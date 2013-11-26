EXEC += scaling
EXEC += error
EXEC += error_img

include ../Makefile.inc

# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
scaling:  $(KERNEL_DIR)/Laplace.o
error:    $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o $(KERNEL_DIR)/Laplace.o $(KERNEL_DIR)/Yukawa.o
error_img: $(KERNEL_DIR)/UnitKernel.o $(KERNEL_DIR)/ExpKernel.o $(KERNEL_DIR)/Laplace.o

# Custom rule for error_img
error_img: % : %.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS) -lpng
