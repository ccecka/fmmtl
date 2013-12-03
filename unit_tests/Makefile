EXEC = compiles

EXEC += correctness
EXEC += dual_correctness

EXEC += multi_level
EXEC += single_level

EXEC += test_p2p
EXEC += test_kern
EXEC += test_tree
EXEC += test_direct
EXEC += test_traits

EXEC += test_vec

EXEC += boost_version
EXEC += thrust_version


include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
compiles:         $(KERNEL_DIR)/KernelSkeleton.o

correctness:      $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o
dual_correctness: $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o

multi_level:      $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o
single_level:     $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o

test_p2p:         $(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))
test_kern:        $(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))
