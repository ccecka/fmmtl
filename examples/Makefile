EXEC += scaling
EXEC += error
EXEC += error_img


include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
scaling:  $(KERNEL_DIR)/Laplace.o
error:    $(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))

error_img: LDLIBS += -lpng
error_img: $(patsubst %.kern,%.o,$(wildcard $(KERNEL_DIR)/*.kern))


