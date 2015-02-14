EXEC = $(basename $(wildcard *.cpp) $(wildcard *.cu))

include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

plr:      LDLIBS += -llapack
plr_test: LDLIBS += -llapack
plr_implicit: LDLIBS += -L/usr/local/Cellar/openblas/0.2.12/lib -lopenblas -llapack
svd_test: LDLIBS += -llapack
probe:    LDLIBS += -llapack

hodlr_solver: LDLIBS += -llapack
