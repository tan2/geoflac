# -*- Makefile -*-
#
# Makefile for flac
#
# Author: Eh Tan <tan2@mail.utexas.edu>
#


#### fortran 90 compiler, compiler flags, linker flags

## for GNU gfortran < 4.2
#F90 = gfortran
#CFLAGS = -O3
#LDFLAGS =

## for GNU gfortran >= 4.2
F90 = gfortran
CFLAGS = -O3 -fopenmp
LDFLAGS = -lgomp

## for Solaris
#F90 = f90
#CFLAGS = -O3 -openmp
#LDFLAGS =

MODULE_SRCS = \
	marker_data.f90

SRCS =	\
	bar2euler.f90 \
	bc_update.f90 \
	dt_adjust.f90 \
	dt_mass.f90 \
	elem2marker.f90 \
	euler2bar.f90 \
	flac.f90 \
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
	par.f90 \
	read_more.f90 \
	read_params.f90 \
	rem_cord.f90 \
	remesh.f90 \
	rem_test.f90 \
	rh_elastic.f90 \
	rh_maxwell.f90 \
	rh_plastic.f90 \
	rh_visc.f90 \
	rmasses.f90 \
	rsflac.f90 \
	saveflac.f90 \
	setflac.f90 \
	user_ab.f90 \
	user_luc.f90

MODULES = $(MODULE_SRCS:.f90=.mod)
MODULE_OBJS = $(MODULE_SRCS:.f90=.o)
OBJS = $(SRCS:.f90=.o)
EXE = flac

all: $(EXE)

$(EXE): $(MODULES) $(OBJS)
	$(F90) $(LDFLAGS) $(MODULE_OBJS) $(OBJS) -o $@

$(MODULES): %.mod : %.f90
	$(F90) $(CFLAGS) -c $<

$(OBJS): %.o : %.f90
	$(F90) $(CFLAGS) -c $<

clean:
	@rm -f $(MODULES) $(MODULE_OBJS) $(OBJS) $(EXE)
