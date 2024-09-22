
# set default fortran compiler:

FC   =  gfortran
#FC  =  ifort


# set mode
MODE = 
# set to "EMPIRE" for EMPIRE mode or leave empty for pure OPTMAN

# set formatted output for *.tlj
OUTMODE =
# set to "formatted" for formatted output of *.tlj or leave blank
# for unformatted

# set parallel
PARALLEL =OPENMP
# set to "OPENMP" for parallelization or leave blank for single thread

# set matrix inversion
MATRIX = 
# set to "LAPACK" for LAPACK matrix inversion or leave blank for
# subroutine



# various flags are set based on compiler FC:
# FFLAGS are the normal complier options for production code
# DFLAGS are the options used when debugging (except for ECIS)
# EFLAGS are the options used for ECIS (no debug allowed)
# OFLAGS are the options used for OPTMAN

LIBS =
DFLAGS = 
FFLAGS = 
EFLAGS = 
OFLAGS = 

ifeq ($(FC),gfortran)

  #---------------------------------
  #----GNU gfortran FORTRAN compiler
  #---------------------------------
  #----flags for production compilation with gfortran
  FFLAGS = -O3 -std=legacy -ftree-vectorize -ffast-math -cpp
  ifeq ($(PARALLEL),OPENMP) 
    FFLAGS += -fopenmp
  endif
  #FFLAGS = -O3 -pg -std=legacy
  #----flags for debuging
  DFLAGS =  -O0 -g --bounds-check -std=legacy -ftrapv 
  # -pg shoudl be added for profiling
  #----flags for OPTMAN
  OFLAGS = -O2 -std=legacy -ffast-math

else ifeq ($(FC),ifort)
 
  #---------------------------
  #----INTEL f95 compiler
  #---------------------------
  #----flags for production compilation using ifort
  # ***  Please note that ifort v11.1 does not work properly with -ipo !!!!!  
  #----flags for debuging
  DFLAGS = -O0 -g -debug all -check all -warn unused -fp-stack-check -ftrapuv -trace -logo
  #------------------------------------------------------------------------------------------------
  # FFLAGS =  -O3 -x=host -logo -parallel -openmp-report1 -par-threshold2 -openmp -vec-report1
  # flags for automatic & openMP directives
  ifeq ($(MATRIX),LAPACK) 
    LIBS =  -mkl#-openmp-lib compat
  endif
  
  #----flags for automatic parallelization
  FFLAGS = -O3 -fpp 
  ifeq ($(PARALLEL),OPENMP) 
    FFLAGS += -qopenmp
  endif
  # Flags for OPTMAN
  OFLAGS = -O2 


endif

ifeq ($(MATRIX),LAPACK)
   FFLAGS = $(FFLAGS) -DLAPACK
endif

ifeq ($(MODE),EMPIRE)
   FFLAGS = $(FFLAGS) -DEMPMODE
endif

ifeq ($(OUTMODE),formatted)
   FFLAGS = $(FFLAGS) -DFORMATTEDOUT
endif

# make sure MAKE knows f90 extension
%.o : src/%.f
	$(FC) $(FFLAGS) -c $<

OBJF = OPTMAND.o SHEMSOFD.o dispers.o KNDITD.o ccrd.o ABCTpar.o DATETpar.o LU_matrix_inv.o opt2Rmatrix.o

all: 
	$(MAKE) optmand

optmand: $(OBJF) 
	$(FC) $(FFLAGS) -o optmand $(OBJF) $(LIBS)

clean:
	rm -f *.o *.mod

cleanall:
	rm -f *.o *.mod empire optmand *.optrpt
