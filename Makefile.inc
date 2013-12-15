#########################
## Library Directories ##
########################

# External libraries
BOOST_DIR = /usr/include

# Optional libraries
# Only needed if using CUDA
CUDA_DIR = /usr/local/cuda-5.5

# This FMM library
FMMTL_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
KERNEL_DIR = $(FMMTL_DIR)kernel

####################
## Makefile Setup ##
####################

# Get the shell name to determine the OS
UNAME := $(shell uname)

# Define the C++ compiler to use
CXX := $(shell which g++) -std=gnu++0x -fPIC
# Check the version number of CXX
CXX_VERSION_LT_47 := $(shell expr `$(CXX) -dumpversion | cut -f1,2 -d.` \< 4.7)
ifeq ($(CXX_VERSION_LT_47),1)
  $(error Makefile using g++ v$(shell $(CXX) -dumpversion), please use >= v4.7)
endif
# Check for CUDA compiler
NVCC := $(shell which nvcc)

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files
INCLUDES += -I. -I$(FMMTL_DIR) -I$(KERNEL_DIR)
INCLUDES += -I$(BOOST_DIR)

# Define cxx compile flags
CXXFLAGS  = -fopenmp -funroll-loops -O3 -W -Wall -Wextra -Wfatal-errors

# Define nvcc compile flags   TODO: Detect and generate appropriate sm_XX
NVCCFLAGS := --compiler-options "$(CXXFLAGS)" -w -arch=sm_20 -O3 #-Xptxas="-v"

# Define any directories containing libraries
#   To include directories use -Lpath/to/files
LDFLAGS +=

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
LDLIBS  +=

######################
## Makefile Options ##
######################

ifeq ($(NDEBUG),1)
CXXFLAGS += -DFMMTL_NDEBUG
endif
ifeq ($(DEBUG),1)
CXXFLAGS += -g -fno-inline -DFMMTL_DEBUG
endif
ifeq ($(PROFILE),1)
CXXFLAGS += -g -pg
endif
ifeq ($(LOG),1)
CXXFLAGS += -DFMMTL_LOGGING
endif
ifeq ($(NO_CUDA),1)
NVCC :=
endif
# Set up for CUDA if available
ifeq ($(NVCC),)    # NVCC is not available
else               # NVCC is available
CXXFLAGS += -DFMMTL_WITH_CUDA
INCLUDES += -I$(CUDA_DIR)/include
LDFLAGS  += -L$(CUDA_DIR)/lib64
LDLIBS   += -lcudart
endif

####################
## Makefile Rules ##
####################

# Suffix replacement rules
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule
.SUFFIXES:                                        # Delete the default suffixes
.SUFFIXES: .hpp .cpp .kern .kern.cpp .kern.cu .o  # Define our suffix list

# 'make' - default rule
all: $(EXEC)

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Default rule for creating a .o file from a .cu file
%.o: %.cu
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c -o $@ $<
	@$(NVCC) $(NVCCFLAGS) $(INCLUDES) -M -o $(DEPSFILE) $<
	@sed -i "1s,.*.o,$@," $(DEPSFILE)  # Fix dep file file.o -> /path/file.o

# Default rule for creating a .o file from a .kern.cpp file
%.o: %.kern.cpp
	$(CXX) -DFMMTL_KERNEL $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Default rule for creating a .o file from a .kern.cu file
%.o: %.kern.cu
	$(NVCC) -DFMMTL_KERNEL $(NVCCFLAGS) $(INCLUDES) -c -o $@ $<
	@$(NVCC) -DFMMTL_KERNEL $(NVCCFLAGS) $(INCLUDES) -M -o $(DEPSFILE) $<
	@sed -i "1s,.*.kern.o,$@," $(DEPSFILE)  # Fix dep file file.o -> /path/file.o

# Default rule for creating a .o file from a .kern file
%.o: %.kern
ifeq ($(NVCC),)  # NVCC isn't available
	ln -s $< $(<:%.kern=%.kern.cpp)
	@$(MAKE) --no-print-directory $@
else             # NVCC is availble
	ln -s $< $(<:%.kern=%.kern.cu)
	@$(MAKE) --no-print-directory $@
endif

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o $(KERNEL_DIR)/*.o
	-$(RM) $(EXEC)
	find $(KERNEL_DIR) -maxdepth 1 -type l -delete
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
-include $(wildcard $(DEPSDIR)/*.d)
