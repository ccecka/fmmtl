#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#
include Makefile.inc

# Get the shell name to determine the OS
UNAME := $(shell uname)

# Define the C compiler to use
ifeq ($(UNAME), Linux)
CXX := g++-4.7 -std=gnu++0x
endif
ifeq ($(UNAME), Darwin)
CXX := g++-mp-4.7 -std=gnu++0x
endif
LINK := $(CXX)
# Check for CUDA compiler
NVCC := $(shell which nvcc)

# Define any compile-time flags
CXXFLAGS := $(ADDED_FLAGS) -fopenmp -funroll-loops -O3 -W -Wall -Wextra
ifeq ($(DEBUG),1)
CXXFLAGS += -g -fno-inline -DDEBUG
endif
ifeq ($(PROFILE),1)
CXXFLAGS += -g -pg
endif
# Dependency directory and flags
DEPSDIR := $(shell mkdir -p .deps; echo .deps)
# MD: Dependency as side-effect of compilation
# MF: File for output
# MP: Include phony targets
DEPCFLAGS = -MD -MF $(DEPSDIR)/$(notdir $*.d) -MP

# Define any directories containing header files other than /usr/include
#   To include directories use -Ipath/to/files
INCLUDES := -I. -I$(FMMW_DIR)/include -I$(FMMW_DIR)/kernel
INCLUDES += -I$(BOOST_DIR)

# Define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX
LIBS :=

# Set up for CUDA if available
ifeq ($(NVCC),)
# NVCC isn't available
CXXFLAGS += -DFMMTL_NO_CUDA
else
# NVCC is available
NVCFLAGS := $(ADDED_FLAGS) -arch=sm_30 -w --compiler-options "$(CXXFLAGS)"
INCLUDES += -I$(THRUST_DIR)
LIBS += -L/usr/local/cuda-5.0/lib64 -lcudart
endif

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and
# adding executables (and any additional dependencies) to EXEC below
##################

EXEC += scaling

# 'make' - default rule
all: $(EXEC)

# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?
KERNEL_DIR := $(FMMW_DIR)/kernel
# List the Kernel dependency (.kern) for each executable
gpup2p:    $(KERNEL_DIR)/UnitKernel.o
serialrun: $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o

# Suffix replacement rules
#   $^: the name of the prereqs of the rule
#   $<: the name of the first prereq of the rule
#   $@: the name of the target of the rule
.SUFFIXES:                    # Delete the default suffixes
.SUFFIXES: .hpp .cpp .kern .o # Define our suffix list

# Default rule for creating an exec of $(EXEC) from a .o file
$(EXEC): % : %.o
	$(LINK) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Default rule for creating a .o file from a .cpp file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPCFLAGS) $(INCLUDES) -c -o $@ $<

# Default rule for creating a .o file from a .cu file
%.o: %.cu
	$(NVCC) $(NVCFLAGS) $(INCLUDES) -c -o $@ $<

# Default rule for creating a .o file from a .kern file
%.o: %.kern
ifeq ($(NVCC),)  # NVCC isn't available
	ln -s $< $(<:%.kern=%.cpp)
	$(MAKE) $@ ADDED_FLAGS='-DFMMTL_KERNEL'
	rm $(<:%.kern=%.cpp)
else          # NVCC is availble
	ln -s $< $(<:%.kern=%.cu)
	$(MAKE) $@ ADDED_FLAGS='-DFMMTL_KERNEL'
	rm $(<:%.kern=%.cu)
endif

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o *~ */*~ $(OBJS) $(FMMW_DIR)/kernel/*.o
	-$(RM) $(EXEC)
	$(RM) -r $(DEPSDIR)

# Define rules that do not actually generate the corresponding file
.PHONY: clean all

# Include the dependency files
DEPFILES := $(wildcard $(DEPSDIR)/*.d)
ifneq ($(DEPFILES),)
-include $(DEPFILES)
endif
