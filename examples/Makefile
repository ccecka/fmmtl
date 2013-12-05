EXEC += scaling
EXEC += error
EXEC += error_img


include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

scaling:  $(KERNEL_DIR)/Laplace.o
error:    $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o

error_img: LDLIBS += -lpng  # Extra png dependency for making images
error_img: $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o


