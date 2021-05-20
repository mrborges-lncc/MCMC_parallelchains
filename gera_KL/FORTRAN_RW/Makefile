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

SOURCES = RANDOM_NUM.f90 KL_condRW.f90

INCLUDE_FILES = 
OBJECTS       = $(SOURCES:.f90=.o)
MODULES	     = *.mod
OUTPUT        = run
MAKEFILE      = Makefile

RM = rm -f

# setting ostype 
OSTYPE    = $(shell uname)

# define standard Fortran90 compiler
FC = ifort
FC = gfortran
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
