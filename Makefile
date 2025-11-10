# Specify a particular compiler with "make COMPILER=pgi", etc.
# Specify debugging flags with "make MODE=debug"
# Alternatively, these options can be specified as environment variables
# eg. "export COMPILER=gfortran" can be added to $HOME/.bashrc


D = -D
FFLAGS = -O3

  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -Wall -Wextra -pedantic -fbounds-check \
             -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
						 -Wrealloc-lhs-all

    GNUVER := $(shell gfortran -dumpversion | head -1 \
        | sed 's/[^0-9\.]*\([0-9\.]\+\).*/\1/')
    GNUMAJOR := $(shell echo $(GNUVER) | cut -f1 -d\.)
    GNUMINOR := $(shell echo $(GNUVER) | cut -f2 -d\.)

    # gfortran-4.3
    GNUGE43 := $(shell expr $(GNUMAJOR) \>= 4 \& $(GNUMINOR) \>= 3)
    ifeq "$(GNUGE43)" "1"
      FFLAGS += -fbacktrace -fdump-core

      # gfortran-4.6
      GNUGE46 := $(shell expr $(GNUMINOR) \>= 6)
      ifeq "$(GNUGE46)" "1"
        FFLAGS += -Wno-unused-dummy-argument

        # gfortran-4.8
        GNUGE48 := $(shell expr $(GNUMINOR) \>= 8)
        ifeq "$(GNUGE48)" "1"
          FFLAGS += -Wno-target-lifetime
        endif
      endif
    endif
  endif
  MODULEFLAG = -I/usr/include -I$(OBJDIR) -J$(OBJDIR)

ifeq ($(shell hostname),hamilton8.dur.ac.uk)
	MPIF90 ?= mpiifort
	NETCDF = -I /usr/local/Cluster-Apps/netcdf-fortran/ompi/gcc/4.4.4/include
	NETCDFLIB = -L/usr/local/Cluster-Apps/netcdf-fortran/ompi/gcc/4.4.4/lib  -lnetcdff
else
	MPIF90 ?= /usr/lib64/openmpi/bin/mpif90
	NETCDF = -I /usr/lib64/gfortran/modules
	NETCDFLIB = -L/usr/lib64/libnetcdff.so.7 -lnetcdff
endif

FFLAGS += $(MODULEFLAG)
LDFLAGS = $(FFLAGS)
# Set some of the build parameters
TARGET = lare3d

# Set pre-processor defines
DEFINES := $(DEFINE)

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

all: main

SRCDIR = src
OBJDIR = obj
BINDIR = bin
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"'

SRCFILES = boundary.f90 conduct.f90 radiative.f90 control.f90 diagnostics.F90 \
  initial_conditions.f90 lagran.F90 lare3d.f90 mpi_routines.F90 \
  mpiboundary.f90 neutral.f90 normalise.f90 openboundary.f90 remap.f90 \
  setup.F90 shared_data.F90 version_data.F90 welcome.f90 xremap.f90 yremap.f90 \
  zremap.f90 random_generator.f90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/core:$(SDF)/src:$(OBJDIR)

-include $(SRCDIR)/COMMIT

ifeq ($(DONE_COMMIT),)
main: commit
else
main: $(FULLTARGET)
endif

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) $(NETCDF) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) $(NETCDF) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS) $(NETCDFLIB)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)

cleanall: tidy

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log
	$(MAKE) -C $(SDF) cleanall

datatidy:
	@rm -rf Data/*

tarball:
	@sh $(SRCDIR)/make_tarball.sh


$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

commit: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh && $(MAKE) $(MAKECMDGOALS) DONE_COMMIT=1

FORCE:

.PHONY: commit clean cleanall tidy datatidy visit visitclean main FORCE

# All the dependencies

boundary.o: boundary.f90 mpiboundary.o random_generator.o shared_data.o
conduct.o: conduct.f90 boundary.o neutral.o shared_data.o
control.o: control.f90 normalise.o shared_data.o
radiative.o: radiative.f90 boundary.o shared_data.o
diagnostics.o: diagnostics.F90 boundary.o conduct.o shared_data.o \
  version_data.o $(SDFMOD)
initial_conditions.o: initial_conditions.f90 neutral.o shared_data.o \
  boundary.o
lagran.o: lagran.F90 boundary.o openboundary.o conduct.o neutral.o shared_data.o \
  openboundary.o remap.o
lare3d.o: lare3d.f90 boundary.o control.o diagnostics.o initial_conditions.o \
  lagran.o mpi_routines.o neutral.o normalise.o openboundary.o remap.o setup.o \
  shared_data.o welcome.o
mpi_routines.o: mpi_routines.F90 shared_data.o
mpiboundary.o: mpiboundary.f90 shared_data.o
neutral.o: neutral.f90 boundary.o shared_data.o
normalise.o: normalise.f90 shared_data.o
openboundary.o: openboundary.f90 shared_data.o
random_generator.o: random_generator.f90
remap.o: remap.f90 boundary.o shared_data.o xremap.o yremap.o zremap.o
setup.o: setup.F90 diagnostics.o shared_data.o version_data.o welcome.o \
  $(SDFMOD)
shared_data.o: shared_data.F90 $(SDFMOD)
version_data.o: version_data.F90 COMMIT
welcome.o: welcome.f90 shared_data.o version_data.o
xremap.o: xremap.f90 boundary.o
yremap.o: yremap.f90 boundary.o
zremap.o: zremap.f90 boundary.o
