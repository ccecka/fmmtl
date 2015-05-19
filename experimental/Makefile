EXEC = $(basename $(wildcard *.cpp) $(wildcard *.cu))

include ../Makefile.inc


# Extra dependencies for executables
# TODO: Automatically detect a .kern file in the dependencies?

# List the Kernel dependency (.kern) for each executable
## A kernel dependency should be listed if a full kernel-matrix vector product
## is executed with that kernel

ifeq ($(host-name),crisco-mac)
	OPENBLAS_LIB = -L/usr/local/Cellar/openblas/0.2.14/lib/
endif
ifeq ($(host-name),crisco-nix)
	OPENBLAS_LIB = -L/opt/OpenBLAS/lib/
endif

svd_test: LDLIBS += -llapack
probe:    LDLIBS += -llapack

plr:      LDLIBS += -llapack
plr_test: LDLIBS += -llapack
plr_implicit: CXXFLAGS += -DWITH_OPENBLAS
plr_implicit: LDLIBS += $(OPENBLAS_LIB) -lopenblas

test_hodlr: CXXFLAGS += -DWITH_OPENBLAS
test_hodlr: LDLIBS += $(OPENBLAS_LIB) -lopenblas

test_hodlr_gpu: CXXFLAGS += -DWITH_OPENBLAS -DWITH_CUBLAS -DWITH_CUSOLVER
test_hodlr_gpu: LDLIBS += $(OPENBLAS_LIB) -lopenblas -lcublas -lcusolver

profile_hodlr_gpu: CXXFLAGS += -DWITH_OPENBLAS -DWITH_CUBLAS -DWITH_CUSOLVER
profile_hodlr_gpu: LDLIBS += $(OPENBLAS_LIB) -lopenblas -lcublas -lcusolver

test_plr: CXXFLAGS += -DWITH_OPENBLAS
test_plr: LDLIBS += $(OPENBLAS_LIB) -lopenblas
