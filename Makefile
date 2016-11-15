#$Rev: 3410 $
#$Author: dbrown $
#$Date: 2013-05-01 03:44:39 +0200 (Wed, 01 May 2013) $

# set default fortran compiler:

#FC   =  gfortran
FC  =  ifort
#FC  =  af95
#FC  =  pfg90
#FC  =  lf95

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
  FFLAGS = -O3 -std=legacy -ftree-vectorize -ffast-math -fopenmp
  #FFLAGS = -O3 -pg -std=legacy
  #----flags for debuging
  DFLAGS =  -O0 -g --bounds-check -std=legacy -ftrapv 
  # -pg shoudl be added for profiling
  #----flags for ECIS
  EFLAGS = -O2 -std=legacy -ffast-math
  #----flags for OPTMAN
  OFLAGS = -O2 -std=legacy -ffast-math
  #----flags for debuging using gfortran compiler
  #----From http://www.macresearch.org/best-fortran-compiler-leopard
  #gfortran -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic

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
  LIBS =  -mkl#-openmp-lib compat
  #----flags for automatic parallelization
  FFLAGS = -O3 -qopenmp#-x=host -parallel -par-report1
  # Flags for ECIS
  EFLAGS = -O2 -x=host
  # Flags for OPTMAN
  OFLAGS = -O2 -x=host -parallel -par-report1
  #----flags for use on NNDC cluster
  ifeq ($(HOSTNAME),nlc.nndc.bnl.gov)
      FFLAGS = -O3
      EFLAGS = -O2
      OFLAGS = -O2
  endif

else ifeq ($(FC),af95)
 
  #---------------------------
  #----Absoft FORTRAN compiler
  #---------------------------
  #----We need to include Unix and Vax libraries to add etime,
  #    getenv, and system functions in Absoft
  LIBS =  -lU77 -lV77 -L/Applications/Absoft11.1/lib64
  #----flags for production
  FFLAGS = -march=host -m64 -O3 -s
  #----flags for debuging
  DFLAGS = -O0 -m64 -g -s
  #----flags for ECIS
  EFLAGS = -O2 -march=host -s -m64
  #----flags for OPTMAN
  OFLAGS = -O2 -march=host -s -m64

else ifeq ($(FC),pfg90)
 
  #---------------------------
  #  Portland Group Fortran Compiler
  #---------------------------
  #----flags for production
  FFLAGS =  -O3 
  #----flags for debuging
  DFLAGS = -O0 -g
  #----flags for ECIS
  EFLAGS = -O2
  #----flags for OPTMAN
  OFLAGS = -O2

else ifeq ($(FC),lf90)
 
  #--------------------------------------
  #----Lahey/Fujitsu f95 compiler
  #--------------------------------------
  #----flag for Lahey/Fujitsu production 
  FFLAGS = -O3
  #----flags for debuging
  DFLAGS = -O0 -g --chk
  #----flags for ECIS
  EFLAGS = -O2
  #----flags for OPTMAN
  OFLAGS = -O2

endif

# make sure MAKE knows f90 extension
%.o : %.f
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
