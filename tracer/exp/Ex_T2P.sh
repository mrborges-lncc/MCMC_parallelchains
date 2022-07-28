#!/bin/bash

export _LMFILES_=/scratch/app/modulos/openmpi/gnu/3.1.4
export MPI_ROOT=/scratch/app/openmpi/3.1.4_gnu
export PKG_CONFIG_PATH=/scratch/app/openmpi/3.1.4_gnu/lib/pkgconfig
export PATH=/scratch/app/openmpi/3.1.4_gnu/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/prj/simulreserv/marcio.borges/bin
export LD_LIBRARY_PATH=/scratch/app/openmpi/3.1.4_gnu/lib
export CPPFLAGS=-I/scratch/app/openmpi/3.1.4_gnu/include
export LDFLAGS=-L/scratch/app/openmpi/3.1.4_gnu/lib
export MANPATH=/scratch/app/openmpi/3.1.4_gnu/share/man:/usr/share/man:/usr/local/share/man
export SHELL=/bin/bash
export MPI_MCA_mca_base_component_show_load_errors=0
export OMPI_PATH=/scratch/app/openmpi/3.1.4_gnu/bin
export OMPI_LD_LIBRARY_PATH=/scratch/app/openmpi/3.1.4_gnu/lib
export MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles:/scratch/app/modulos
export LOADEDMODULES=
export MODULESHOME=/usr/share/Modules


INSTDIR=/scratch/simulreserv/rafael.guiraldello/piritech_proj1

#exibe informaçoes sobre o executavel
EXEC=${INSTDIR}/build/linux_gcc/partrans_linux_gcc_opt
GLOBAL_SOL="-ksp_type gmres -ksp_gmres_restart 20 -ksp_rtol 1e-8 -pc_type hypre"
CAP="-capPress_ksp_type cg -capPress_ksp_gmres_restart 20 -capPress_ksp_rtol 1e-8 -capPress_pc_type hypre"


#PATH=/scratch/app/openmpi/3.1.4_gnu/bin/:/bin mpirun -n 4  --oversubscribe $EXEC $GLOBAL_SOL $CAP
mpirun -n 4 --oversubscribe $EXEC $GLOBAL_SOL $CAP
