EXEC = $(basename $(wildcard *.cpp) $(wildcard *.cu))

include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

svd_test: LDLIBS += -llapack
probe:    LDLIBS += -llapack

plr:      LDLIBS += -llapack
plr_test: LDLIBS += -llapack
plr_implicit: CXXFLAGS += -DWITH_OPENBLAS
plr_implicit: LDLIBS += -L/opt/OpenBLAS/lib -lopenblas -llapack

test_hodlr: CXXFLAGS += -DWITH_OPENBLAS
test_hodlr: LDLIBS += -L/opt/OpenBLAS/lib -lopenblas -llapack
