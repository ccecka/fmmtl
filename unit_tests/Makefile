EXEC = $(basename $(wildcard *.cpp) $(wildcard *.cu))

include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

compiles:         $(KERNEL_DIR)/KernelSkeleton.o

correctness:      $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o
dual_correctness: $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o

test_s2t:         $(KERNEL_DIR)/Laplace.o #$(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))

test_bbfmm:       $(KERNEL_DIR)/Barycentric.o
