#########################
## Library Directories ##
########################

# This FMM library
FMMTL_DIR := $(dir $(lastword $(MAKEFILE_LIST)))
KERNEL_DIR = $(FMMTL_DIR)kernel

# External libraries
THRUST_DIR = /usr/local/cuda/include

# Optional libraries
# Only needed if using CUDA
CUDA_PATH = /usr/local/cuda


####################
## Makefile Setup ##
####################

# Get the host-name if empty
ifeq ($(host-name),)
	host-name := $(shell hostname)
endif
# Get the kernel-name if empty
ifeq ($(kernel-name),)
	kernel-name := $(shell uname -s)
endif
# Get the arch-name if empty
ifeq ($(arch-name),)
	arch-name := $(shell uname -p)
endif

# Define the C++ compiler to use
CXX := $(shell which g++) -std=c++11
# Check the version number of CXX
CXX_VERSION_LT_48 := $(shell expr `$(CXX) -dumpversion | cut -f1,2 -d.` \< 4.8)
ifeq ($(CXX_VERSION_LT_48),1)
  $(error Makefile using g++ v$(shell $(CXX) -dumpversion), please use >= v4.8)
endif
# Check for CUDA compiler
USE_NVCC := $(shell which nvcc)
NVCC     := $(USE_NVCC) -std=c++11 -ccbin=g++-4.8 #$(shell which g++)

# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPSFILE = $(DEPSDIR)/$(notdir $*.d)
DEPSFLAGS = -MD -MF $(DEPSFILE) #-MP

# Define any directories containing header files
#   To include directories use -Ipath/to/files
INCLUDES += -I. -I$(THRUST_DIR) -I$(FMMTL_DIR) -I$(KERNEL_DIR)

# Define cxx compile flags
CXXFLAGS  = -fopenmp -funroll-loops -O3 -Wall -Wextra -Wno-unused-local-typedefs #-Wfatal-errors

# Define nvcc compile flags   TODO: Detect and generate appropriate sm_XX
NVCCFLAGS := -arch=sm_35 -O3 --compiler-options "$(CXXFLAGS)" -Xcompiler -Wno-unused-parameter #-Xptxas="-v"

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
	CXXFLAGS += -DNDEBUG
	NVCCFLAGS += -DNDEBUG
endif
ifeq ($(FMMTL_NDEBUG),1)
	CXXFLAGS += -DFMMTL_NDEBUG
endif
ifeq ($(DEBUG),1)
	CXXFLAGS += -DFMMTL_DEBUG -g -fno-inline
endif
ifeq ($(PROFILE),1)
	CXXFLAGS += -g -pg
endif
ifeq ($(LOG),1)
	CXXFLAGS += -DFMMTL_LOGGING
endif
ifeq ($(NO_CUDA),1)
	USE_NVCC :=
endif
# Set up for CUDA if available
ifeq ($(USE_NVCC),)    # NVCC is not available
else                   # NVCC is available
	CXXFLAGS += -DFMMTL_WITH_CUDA
	# Use cuda lib64 if it exists, else cuda lib
	ifneq ($(wildcard $(CUDA_PATH)/lib64/.*),)
		LDFLAGS += -L$(CUDA_PATH)/lib64
	else
		LDFLAGS += -L$(CUDA_PATH)/lib
	endif
	LDLIBS += -lcudart
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

# Default rule for creating a .o file from a .kern.cpp file
%.o: %.kern.cpp
	$(CXX) -DFMMTL_KERNEL $(CXXFLAGS) $(INCLUDES) $(DEPSFLAGS) -c -o $@ $<

# Default rule for creating a .o file from a .kern.cu file
%.o: %.kern.cu
	$(NVCC) -DFMMTL_KERNEL $(NVCCFLAGS) $(INCLUDES) -c -o $@ $<
	@$(NVCC) -DFMMTL_KERNEL $(NVCCFLAGS) $(INCLUDES) -M -o $(DEPSFILE) $<
	@sed -i -e "1s,.*.kern.o,$@," $(DEPSFILE)  # Fix dep file file.o -> /path/file.o

# Default rule for creating a .o file from a .kern file
%.o: %.kern
ifeq ($(USE_NVCC),)  # NVCC isn't available
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
