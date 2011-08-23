# -*- Makefile -*-
#
# Makefile for flac
#
# Author: Eh Tan <tan2@mail.utexas.edu>
#

## Execute "make" if making production run. Or "make debug=1" for debugging run.
##
## debug = 0: optimized build; 1: debugging build
## omp = 1: enable openmp support
## coverage = 1: enable code coverage analysis (used for debugging)
## gprof = 1: enable profiling (used for tuning optimization)
## cuda = 0: f90 only; 1: cuda runtime; 2: cuda emulation

debug = 0
omp = 1
coverage = 0
gprof = 0
cuda = 0


## any customized flags
CFLAGS =
LDFLAGS =


## Select default fortran 90 compiler

## GNU version >= 4.2, otherwise openmp won't work
F90 = gfortran

## Intel, tested with v11.1
#F90 = ifort

## PGI group
#F90 = pgf90

## Sun
#F90 = f90


########################################################################
## Select compiler and linker flags
## (Usually you won't need to modify anything below)
########################################################################

## for GNU gfortran >= 4.2
ifeq ($(F90), gfortran)
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -O3 -ffast-math
	else
		CFLAGS += -O0 -g -Wall -fbounds-check -ftrapv -fbacktrace -fdump-core \
			-Waliasing -Walign-commons -Wampersand -Wintrinsic-shadow \
			-Wline-truncation -Wsurprising -Wunderflow \
			-finit-integer=-99 -finit-real=nan -frange-check -Werror
	endif

	ifeq ($(coverage), 1)
		CFLAGS += -fprofile-arcs -ftest-coverage
		LDFLAGS += -fprofile-arcs -ftest-coverage
	endif

	ifeq ($(gprof), 1)
		CFLAGS += -pg -g
		LDFLAGS += -pg -g
	endif
endif

## for Intel v11.1
ifeq ($(F90), ifort)
	CFLAGS += -assume byterecl
	ifeq ($(omp), 1)
		CFLAGS += -fopenmp
		LDFLAGS += -fopenmp
	endif
	ifeq ($(debug), 0)
		CFLAGS += -fast
	else
	       CFLAGS += -O0 -g -no-opt-subscript-in-range -ftrapuv -check all -traceback
	       LDFLAGS +=
	endif
endif

## for PGI
ifeq ($(F90), pgf90)
	ifeq ($(omp), 1)
		CFLAGS += -mp
		LDFLAGS += -mp
	endif
        ifeq ($(debug), 0)
                CFLAGS += -O3 -fast -fastsse -Mnoframe -Mipa=fast
        else
               CFLAGS += -O0 -g -Mbounds -Mchkstk -traceback
               LDFLAGS +=
        endif
endif

## for Solaris
ifeq ($(F90), f90)
	ifeq ($(omp), 1)
		CFLAGS += -xopenmp=parallel
		LDFLAGS += -xopenmp
	endif
	# -r8const: Promote single precision const to double precision.
	# -C: Enable runtime subscript range checking
	# -fpover=yes: Provide run-time overflow check during READ
	# -xcheck=%all: Enable all run-time check
	ifeq ($(debug), 0)
		CFLAGS += -O3 -r8const -fast -xloopinfo -s -xipo
	else
		CFLAGS += -g -w3 -r8const -C -fpover=yes -xcheck=%all
	endif
endif


#### CUDA C compiler, compiler flags, linker flags
ifneq ($(cuda), 0)

	ifdef TACC_CUDA_DIR
		CUDAHOME = $(TACC_CUDA_DIR)
		CUDALIB = $(TACC_CUDA_LIB)
	else
		CUDAHOME = $(HOME)/cuda
		CUDALIB = $(CUDAHOME)/lib64
	endif

	CFLAGS += -DUSE_CUDA
	LDFLAGS += -L$(CUDALIB) -lstdc++ -Xlinker -rpath $(CUDALIB)

	NVCC = $(CUDAHOME)/bin/nvcc
	NVCCFLAGS = -use_fast_math --ptxas-options=-v -arch sm_13 # for double precision

	ifeq ($(debug), 0)
		NVCCFLAGS += -O2
	else
		NVCCFLAGS += -g
	endif

	ifeq ($(cuda), 2)
		NVCCFLAGS += -g -deviceemu
		LDFLAGS += -lcudartemu
	else
		LDFLAGS += -lcudart
	endif
endif


## command to invoke mercurial VCS
HG = hg


MODULE_SRCS = \
	marker_data.f90

# f90 source code with preprocessor directives (#if ... #endif)
CPP_SRCS = \
	arrays.F90 \
	flac.F90 \
	par.F90

# f90 source code without preprocessor directives
SRCS =	\
	bar2euler.f90 \
	bc_update.f90 \
	dt_adjust.f90 \
	dt_mass.f90 \
	euler2bar.f90 \
	fl_move.f90 \
	fl_node.f90 \
	fl_rheol.f90 \
	fl_srate.f90 \
	fl_therm.f90 \
	functions.f90 \
	init_areas.f90 \
	init_bc.f90 \
	init_cord.f90 \
	init_marker.f90 \
	init_phase.f90 \
	init_stress.f90 \
	init_temp.f90 \
	init_tracer.f90 \
	init_visc.f90 \
	lpeuler2bar.f90 \
	marker2elem.f90 \
	matprops.f90 \
	newphase2marker.f90 \
	outflac.f90 \
	outmarker.f90 \
	outtracer.f90 \
	read_params.f90 \
	rem_cord.f90 \
	remesh.f90 \
	rem_test.f90 \
	rh_elastic.f90 \
	rh_maxwell.f90 \
	rh_plastic.f90 \
	rmasses.f90 \
	rsflac.f90 \
	saveflac.f90 \
	setflac.f90 \
	user_ab.f90 \
	user_luc.f90

CUDAFLAC_SRC = \
	cu_flac.cu

MODULES = $(MODULE_SRCS:.f90=.mod)
MODULE_OBJS = $(MODULE_SRCS:.f90=.o)
OBJS = $(SRCS:.f90=.o)
OBJS2 = $(CPP_SRCS:.F90=.o)
CUDAFLAC_OBJS = $(CUDAFLAC_SRC:.cu=.o)

EXE = flac
CUEXE = cuflac

ifeq ($(cuda), 0)
all: $(EXE)
else
all: $(CUEXE)
endif

$(EXE): $(MODULE_OBJS) $(OBJS2) $(OBJS)
	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) $(OBJS2) -o $@
	@which hg &> /dev/null && (hg summary && hg diff) > 0.patch

$(CUEXE): $(MODULES_OBJS) $(OBJS2) $(OBJS) $(CUDAFLAC_OBJS)
	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) $(OBJS2) $(CUDAFLAC_OBJS) -o $@
	@which hg &> /dev/null && (hg summary && hg diff) > 0.patch

$(MODULE_OBJS): %.o : %.f90
	$(F90) $(CFLAGS) -c $<

$(OBJS): %.o : %.f90
	$(F90) $(CFLAGS) -c $<

$(OBJS2): %.o : %.F90
	$(F90) $(CFLAGS) -c $<

$(CUDAFLAC_OBJS): %.o : %.cu
	$(NVCC) $(NVCCFLAGS) -L$(CUDALIB) -c $<

touch:
	@rm -f $(OBJS2)

clean:
	@rm -f $(MODULES) $(MODULE_OBJS) $(OBJS) $(OBJS2) $(CUDAFLAC_OBJS) $(EXE)

clean-data:
	@rm -f *.0 *.rs sys.msg marker_init.asc vbc.s area.dat output.asc temp0.dat
