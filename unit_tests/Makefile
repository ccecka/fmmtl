EXEC = compiles

EXEC += correctness
EXEC += dual_correctness

EXEC += multi_level
EXEC += single_level

EXEC += test_kernel
EXEC += test_expansion

EXEC += test_p2p
EXEC += test_tree
EXEC += test_direct
EXEC += test_vec

EXEC += version


include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

compiles:         $(KERNEL_DIR)/KernelSkeleton.o

correctness:      $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o
dual_correctness: $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o

test_p2p:         $(KERNEL_DIR)/Laplace.o #$(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))
