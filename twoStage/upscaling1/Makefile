###
###  Makefile for some tests of fortran90/95 compilers
###

#
# notes:
# OSTYPE : defined by system
# F90    : defined by user in shell;  nag, lahey
# DEBUG 0/1
#

SRCDIR  = LAPACK77

SOURCES = upscaling.f90 bloco.f90

INCLUDE_FILES = 
OBJECTS       = $(SOURCES:.f90=.o)
MODULES	     = *.mod
OUTPUT        = run
MAKEFILE      = Makefile

RM = rm -f

# setting ostype 
OSTYPE    = $(shell uname)

# define standard Fortran90 compiler
#F90= gfortran
FC = gfortran
#FC = ifort
LD = $(FC)
# define libs
#BLAS      = -L/home/mrborges/LAPACK/blas_LINUX
#LAPACK77  = -L/home/mrborges/LAPACK/lapack_LINUX
BLAS      = -lblas
LAPACK77  = -llapack
LIBS      = $(LAPACK77) $(BLAS) $(SPECIAL_LIBG2C)

ifeq ($(FC),ifort)
	FFLAGS    = -fast
endif

ifeq ($(FC),gfortran)
	FFLAGS    = -O2
endif

ifeq ($(OSTYPE), Linux)
##### linux #####
  ifeq ($(F90), nag)
    ### NAGf95 ###
    ## disable warnings from license manager
    export NAG_LM_OPTS=nowarn
    S_BLAS    = libblas.a libblas.so
    S_LAPACK  = liblapack.a liblapack.so
    FC        = f95-nag
    LD	      = $(FC)
    ifeq ($(DEBUG), 1)
      DEBUGFLAG = -C -g -g90 -gline
    else
      DEBUGFLAG = 
    endif
    FFLAGS    = -v -u -O4 -nan
    LDFLAGS   = 
    BLAS      = -lblas
    LAPACK77  = -llapack
    ifeq ($(DEBUG), 1)
      LIBS      = $(LAPACK77) $(BLAS) $(SPECIAL_LIBG2C) -lefence
    else
      LIBS      = $(LAPACK77) $(BLAS) $(SPECIAL_LIBG2C)
    endif
    CUT_ASM_WARN = 
  endif

  ifeq ($(F90), lahey)
    ### Lahey ###
    S_BLAS    = libblas.a libblas.so
    S_LAPACK  = liblapack.a liblapack.so
    FC        = f95-lah
    LD        = $(FC)
    ifeq ($(DEBUG), 1)
      DEBUGFLAG = --chk -g --trace
    else
      DEBUGFLAG = 
    endif
    FFLAGS    = --wo --warn --f95 -O --tpp --ap
    LDFLAGS   = 
    BLAS      = -lblas
    LAPACK77  = -llapack
    ifeq ($(DEBUG), 1)
      LIBS      = $(LAPACK77) $(BLAS) $(SPECIAL_LIBG2C) -lefence
    else
      LIBS      = $(LAPACK77) $(BLAS) $(SPECIAL_LIBG2C)
    endif
    CUT_ASM_WARN = 2>&1 | grep -v "/tmp/asm"
  endif
endif


$(OUTPUT): $(OBJECTS) 

	$(LD) $(OBJECTS) -o $(OUTPUT) $(LDFLAGS) $(LIBS) $(DEBUGFLAG)

%.o : %.f90 $(MAKEFILE) 
	$(FC) $(@:.o=.f90) -c -o $@ $(DEBUGFLAG) $(FFLAGS) $(CUT_ASM_WARN)


install: $(OUTPUT) 
	cp $(OUTPUT) $(HOME)/bin/$(OUTPUT)


basic-clean:
	$(RM) $(OBJECTS) $(MODULES)

clean: basic-clean
	( $(RM) *.o *.mod *.g90 *~ core \#* $(OUTPUT) )

dist:
	( cd ..; tar -zcvf $(SRCDIR)-`date +"%Y-%m-%d"`.tar.gz \
	  $(SRCDIR)/*.f90 $(SRCDIR)/test_lapack_in.dat \
	  $(SRCDIR)/Doc \
	  $(SRCDIR)/Makefile $(SRCDIR)/README \
	  $(SRCDIR)/no_arch -X $(SRCDIR)/no_arch )

distclean: clean dist
