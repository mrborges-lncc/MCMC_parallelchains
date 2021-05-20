#diretorios
DIRFONTES = ./fontes
DIRBIN    = ./bin
DIRINC    = ./include
PARDISO_DIR="/home/mrborges/Dropbox/pardiso/"
export PARDISO_LIC_PATH=/home/mrborges/Dropbox/pardiso/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/mrborges/Dropbox/pardiso/

#opcoes de compilacao
#PETSC     = withPetsc
OMP       = withOMP
UMFPACK   = #withumfpack
PARDISO   = #withpardiso
LAPACK    = #withlapack
#HYPRE     = withHYPRE

OPT       =-O0 -g #-fast -g 


TEMPOS    = mostrarTempos
DEBUG     = #debug

#TIPOSOLUCAO = withLagrangeano
TIPOSOLUCAO = withRT

#compiladores
CC        = gcc
FC        = gfortran
CC        = icc
#FC        = ifort
#FC        = mpif90

ifeq ($(FC),ifort)
   FFLAGSI  = -module ${DIRINC} ${OPT} 
   FOMP     = -openmp -fast # -g  #-fast #-g #-warn  #para ifort
   LIBBLAS  = -mkl
endif

ifeq ($(FC),gfortran)
   FFLAGSI  = -J ${DIRINC}  
   FOMP     = -fopenmp -O3 #-g #-finit-local-zero  #-Wall #para gfortran
   LIBBLAS  = #-lblas -llapack
endif

ifeq ($(OMP),withOMP)
   FPPFLAGS  = -D${OMP} 
   FFLAGS    = ${FFLAGS0} ${FOMP}
   LFLAGS    = ${LFLAGS0} ${FOMP} 
else
   FFLAGS    = ${FFLAGS0} 
   LFLAGS    = ${LFLAGS0}
endif

NOMEEXECUTAVEL=escImiscivel.exe

#NOMEEXECUTAVEL=escImiscivelSkyLine.exe

ifeq ($(LAPACK),withlapack)
   FPPFLAGSL = -D${LAPACK}   
   LFLAGS0   = $(LIBBLAS) 
#   NOMEEXECUTAVEL=escImiscivelLAPACK.exe
endif

LIB=$(LIBBLAS)

ifeq ($(UMFPACK),withumfpack)
   UMFPACK_DIR  = /usr/local/lib/
   UMFPACK_INC  = /local/include/
   UMFPACKCMP  = -DLP64 -I ${UMFPACK_DIR} 
   #UMFPACKLNK  = ${DIRBIN}/umf4_f77wrapper.o -L ${UMFPACK_DIR} -lumfpack -lm -lamd -lsuitesparseconfig -lrt
   UMFPACKLNK  = -L ${UMFPACK_DIR} -lumfpack -lm -lamd -lsuitesparseconfig -lrt
   FPPFLAGSU = -D${UMFPACK}   
endif


ifeq ($(HYPRE),withHYPRE)
   HYPRE_DIR  = /usr/local/hypre2.9
   HYPRE_INC  = /local/include/
   HYPRELNK  = "-L${HYPRE_DIR}/lib -L${HYPRE_DIR}/../babel-runtime/sidl/.libs -lbHYPREClient-C -lbHYPREClient-CX -lbHYPREClient-F -lbHYPRE  -lHYPRE "
   FPPFLAGSU = -D${HYPRE}   
endif


ifeq ($(TEMPOS),mostrarTempos)
   FPPFLAGSTime = -D${TEMPOS}
endif

ifeq ($(DEBUG),debug)
   FPPFLAGSDebug = -D${DEBUG}
endif

ifeq ($(UMFPACK),withumfpack)
   NSYMCRS=withcrs
   FPPFLAGSN = -D${NSYMCRS}
endif

ifeq ($(PARDISO),withpardiso)
   FPPFLAGSP = -D${PARDISO}
   NSYMCRS=withcrs
   FPPFLAGSN = -D${NSYMCRS}
endif

ifeq ($(HYPRE),withHYPRE)
   FPPFLAGSH = -D${HYPRE}
   NSYMCRS=withcrs
   FPPFLAGSN = -D${NSYMCRS}
endif


ifeq ($(TIPOSOLUCAO),withLagrangeano)
FPPFLAGST  = -D${TIPOSOLUCAO}
else
ifeq ($(TIPOSOLUCAO),withRT)
FPPFLAGST  = -D${TIPOSOLUCAO}
endif
endif


FFLAGS0   = ${FPPFLAGS}  ${FPPFLAGSL}  ${FPPFLAGSU} ${FPPFLAGSP} ${FPPFLAGSN} ${FPPFLAGSU}  ${FFLAGSI} ${FPPFLAGST} ${FPPFLAGSTime} ${FPPFLAGSDebug}

NOMEMODULO01=mGlobais
NOMEMODULO02=algMatricial
NOMEMODULO03=malha
NOMEMODULO04=propGeoFisica
NOMEMODULO05=leituraEscrita
NOMEMODULO06=funcoesDeForma
NOMEMODULO07=hidrodinamica
NOMEMODULO08=geomecanica
NOMEMODULO09=mcmc
NOMEMODULO10=transporte
NOMEMODULO11=utilitarios
NOMEMODULO12=driverRigido

OBJECTS0 = ${DIRBIN}/${NOMEMODULO01}.o ${DIRBIN}/${NOMEMODULO02}.o ${DIRBIN}/${NOMEMODULO03}.o \
	${DIRBIN}/${NOMEMODULO04}.o ${DIRBIN}/${NOMEMODULO05}.o ${DIRBIN}/${NOMEMODULO06}.o  \
	${DIRBIN}/${NOMEMODULO07}.o ${DIRBIN}/${NOMEMODULO08}.o ${DIRBIN}/${NOMEMODULO09}.o  \
	${DIRBIN}/${NOMEMODULO10}.o ${DIRBIN}/${NOMEMODULO11}.o ${DIRBIN}/${NOMEMODULO12}.o

SOURCES0 = ${DIRFONTES}/${NOMEMODULO01}.F90 ${DIRFONTES}/${NOMEMODULO02}.F90 ${DIRFONTES}/${NOMEMODULO03}.F90 \
	${DIRFONTES}/${NOMEMODULO04}.F90 ${DIRFONTES}/${NOMEMODULO05}.F90 ${DIRFONTES}/${NOMEMODULO06}.F90 \
	${DIRFONTES}/${NOMEMODULO07}.F90 ${DIRFONTES}/${NOMEMODULO08}.F90 ${DIRFONTES}/${NOMEMODULO09}.F90 \
	${DIRFONTES}/${NOMEMODULO10}.F90 ${DIRFONTES}/${NOMEMODULO11}.F90 ${DIRFONTES}/${NOMEMODULO12}.F90

ifeq ($(UMFPACK),withumfpack)
  OBJECTS = ${OBJECTS0} ${DIRBIN}/umf4_f77wrapper.o
  SOURCES = ${SOURCES0} ${DIRFONTES}/umf4_f77wrapper.c
else
  OBJECTS = ${OBJECTS0}
  SOURCES = ${SOURCES0}
endif


all: ${DIRBIN}/${NOMEEXECUTAVEL} 


ifeq ($(UMFPACK),withumfpack)
${DIRBIN}/umf4_f77wrapper.o: ${DIRFONTES}/umf4_f77wrapper.c
	@echo  0, compilando ${DIRFONTES}/umf4_f77wrapper.c 
	-$(CC) -c -I ${UMFPACK_DIR} ${DIRFONTES}/umf4_f77wrapper.c -o ${DIRBIN}/umf4_f77wrapper.o
endif

${DIRBIN}/${NOMEMODULO01}.o: ${DIRFONTES}/${NOMEMODULO01}.F90
	@echo  1, compilando ${NOMEMODULO01}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO01}.o  ${DIRFONTES}/${NOMEMODULO01}.F90

${DIRBIN}/${NOMEMODULO02}.o: ${DIRFONTES}/${NOMEMODULO02}.F90
	@echo  2, compilando ${NOMEMODULO02}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO02}.o  ${DIRFONTES}/${NOMEMODULO02}.F90

${DIRBIN}/${NOMEMODULO03}.o: ${DIRFONTES}/${NOMEMODULO03}.F90
	@echo  3, compilando ${NOMEMODULO03}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO03}.o  ${DIRFONTES}/${NOMEMODULO03}.F90

${DIRBIN}/${NOMEMODULO04}.o: ${DIRFONTES}/${NOMEMODULO04}.F90
	@echo  4, compilando ${NOMEMODULO04}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO04}.o  ${DIRFONTES}/${NOMEMODULO04}.F90

${DIRBIN}/${NOMEMODULO05}.o: ${DIRFONTES}/${NOMEMODULO05}.F90
	@echo  5, compilando ${NOMEMODULO05}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO05}.o  ${DIRFONTES}/${NOMEMODULO05}.F90

${DIRBIN}/${NOMEMODULO06}.o: ${DIRFONTES}/${NOMEMODULO06}.F90
	@echo  6, compilando ${NOMEMODULO06}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO06}.o  ${DIRFONTES}/${NOMEMODULO06}.F90

${DIRBIN}/${NOMEMODULO07}.o: ${DIRFONTES}/${NOMEMODULO07}.F90
	@echo  7, compilando ${NOMEMODULO07}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO07}.o  ${DIRFONTES}/${NOMEMODULO07}.F90

${DIRBIN}/${NOMEMODULO08}.o: ${DIRFONTES}/${NOMEMODULO08}.F90
	@echo  8, compilando ${NOMEMODULO08}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO08}.o  ${DIRFONTES}/${NOMEMODULO08}.F90

${DIRBIN}/${NOMEMODULO09}.o: ${DIRFONTES}/${NOMEMODULO09}.F90
	@echo  9, compilando ${NOMEMODULO09}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO09}.o  ${DIRFONTES}/${NOMEMODULO09}.F90

${DIRBIN}/${NOMEMODULO10}.o: ${DIRFONTES}/${NOMEMODULO10}.F90
	@echo 10, compilando ${NOMEMODULO10}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO10}.o  ${DIRFONTES}/${NOMEMODULO10}.F90

${DIRBIN}/${NOMEMODULO11}.o: ${DIRFONTES}/${NOMEMODULO11}.F90
	@echo 11, compilando ${NOMEMODULO11}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO11}.o  ${DIRFONTES}/${NOMEMODULO11}.F90

${DIRBIN}/${NOMEMODULO12}.o: ${DIRFONTES}/${NOMEMODULO12}.F90
	@echo 11, compilando ${NOMEMODULO12}.F90
	-${FC} -c ${FFLAGS} -o ${DIRBIN}/${NOMEMODULO12}.o  ${DIRFONTES}/${NOMEMODULO12}.F90




${DIRBIN}/${NOMEEXECUTAVEL}: ${OBJECTS}
	@echo  gerando executavel: ${DIRBIN}/${NOMEEXECUTAVEL} 
	-${FC} ${LFLAGS} -o ${DIRBIN}/${NOMEEXECUTAVEL}  ${OBJECTS} ${LIB} ${UMFPACKLNK} 
clean:
	@echo  apagando:  ${DIRBIN}/${NOMEEXECUTAVEL} ${OBJECTS} ${DIRINC}/*.mod
	rm -f ${DIRBIN}/${NOMEEXECUTAVEL} ${OBJECTS} ${DIRINC}/*.mod

