#set -x
#SCRIPT PARA AUXILIAR NA EXECUCAO DE PROGRAMAS 
# DE SIMULACAO DE ESCOAMENTOS EM RESERVATORIOS 
#
#
#AUTORES: bidu@lncc.br e tuane@lncc.br
# 06.2013
#
# foma de utilizacao . rodarSimulador.sh exp10x10x30
#### DEFINICAO DO NUMERO DE PROCESSOS (MPI) ####
npPadrao=1
numProcs=${1:-${npPadrao}}
export NP=${numProcs} # atribuir valor padrao aa variavel NP


#### EXPERIMENTOS ####

# diretorio do experimento:
# dirPadrao ou lido de primeiro argumento da linha de comando
#dirExp="$(pwd)"/${1:-${dirPadrao}} 
dirExp="$(pwd)"
#rm -rf ${dirExp}/out/*.vtk
arqTela=${dirExp}/output.out
#
#### DEFINICAO DO EXECUTAVEL ####
LOCAL="$(pwd)"
NOMEEXECUTAVEL=runMCMC

DIRBIN="${LOCAL}"

EXECUTAVEL=${DIRBIN}/${NOMEEXECUTAVEL}

#### definicao do comando a ser executado
#comando="(export  OMP_NUM_THREADS=${numThreads} ; cd ${dirExp}; mpirun -np ${numProcs} ${EXECUTAVEL})"
comando=" cd ${dirExp}; mpirun -np  ${numProcs} ${EXECUTAVEL}"

if [ -e ${EXECUTAVEL} ] 
then
  printf "\n diretorio do experimento.: %s\n" ${dirExp}  
  printf "\n nome do executavel.......: %s\n" ${EXECUTAVEL} 
  printf "\n numero de processos......: %d\n" ${NP}
  printf "\n comando .................: %s\n" "${comando}"
  eval ${comando}  |tee  ${arqTela}
else
  printf "\n EXECUTAVEL NAO ENCONTRADO \n"
  printf "\n comando .................: %s\n" "${comando}"
fi
