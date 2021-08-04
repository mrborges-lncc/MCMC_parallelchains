!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 METODO MCMC DOIS ESTAGIOS
!         PARA SELECIONAR CAMPOS DE PERMEABILIDADES,
!               POROSIDADES E MODULO DE YOUNG
!                             MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NAME..................: MARCIO RENTES BORGES
! DATA..................: 09/03/2018
! DATA VERSIONAMENTO....: 04/10/2018
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================================================!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
!
  USE MPI
  USE VARIAVEIS
  USE RANDOM
  USE STATFUNCTIONS
  IMPLICIT NONE
!  INCLUDE 'mpif.h'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LIST OF LOCAL VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER            :: NERROREAD,K,RC,LEN,NK,NCONT
  REAL, DIMENSION(2) :: TARRAY,AUX
  CHARACTER(LEN=128) :: COMMAND
  INTEGER            :: NPROCS, IERROR, NRANK, I, J, TAG
  INTEGER            :: MPI_NSIMUL,MPISTATUS(MPI_STATUS_SIZE)
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) :: RANKNAME
  INTEGER            :: NDIM,NDIMR
  CHARACTER(LEN=1)   :: NSINCRO
  INTEGER, EXTERNAL  :: INTMIN
  INTEGER            :: CONTMIN,ROOT
  REAL               :: PCRLOCAL,LOGERROR_LOCAL
  REAL   ,ALLOCATABLE,DIMENSION(:) :: LOGERROR
  INTEGER,ALLOCATABLE,DIMENSION(:) :: INDX
  LOGICAL            :: STARTSTATUS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ROOT = 0
  NSINCRO = 'Y'
  STARTSTATUS = .FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INICIO DA PARALELIZACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL MPI_INIT(IERROR)
  IF(IERROR.NE.MPI_SUCCESS)THEN
     WRITE(*,1313)
     CALL MPI_ABORT(MPI_COMM_WORLD,RC,IERROR)
  END IF
!
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,NRANK,IERROR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPROCS,IERROR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NRANK.EQ.0)WRITE(*,200)NPROCS
  CALL MPI_GET_PROCESSOR_NAME(RANKNAME,LEN,IERROR)
!
  WRITE(*,201)NRANK,ADJUSTL(TRIM(RANKNAME))
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD = 0
  CALL READ_INPUT(NERROREAD,MPI_NSIMUL,NPROCS,NRANK)
!
  IF(NERROREAD.NE.1)THEN
     WRITE(*,*)'#### LEITURA ERRADA ####'
     WRITE(*,*)'### DO ARQUIVO INPUT ###'
     CALL MPI_FINALIZE(IERROR)
     STOP 1
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(LOGERROR(NPROCS))
  ALLOCATE(INDX(NPROCS))
  LOGERROR  = 0.0
  NERROREAD = 0
  CALL INITIAL(NERROREAD,NRANK,NPROCS,STARTSTATUS)
!
  IF(NERROREAD.NE.1)THEN
     WRITE(*,*)'######### PROBLEMA NA ###########'
     WRITE(*,*)'### INICIALIZACAO DO PROGRAMA ###'
     CALL MPI_FINALIZE(IERROR)
     STOP 2
  END IF
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL DREAMTYPEINITIALIZATION(STARTSTATUS,NRANK,NSINCRO)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COPIA DOS DIRETORIOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(STARTSTATUS.EQV..FALSE.)THEN
     CALL COPY_DIR(NRANK,NSIMUL,NSTAGE)
     CALL COPY_DIRUP(NRANK,NSTAGE)
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!#####################################################!!
!!#####################################################!!
!!#####################################################!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LOOP SOBRE AS REALIZACOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CONTMIN = 0
  VETCONT = 0
  DO WHILE(CONTMIN.LT.NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SALVA O STATUS ATUAL DA SIMULACAO !!!!!!!!!!!!!!!!!!!!!
     CALL SAVE_INITSTATUS(CONTADORC,CONTADORF,&
          CONT,NCHAIN,NRANK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     CONTADORC=CONTADORC+1
!
999  CONTINUE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5&
          .OR.NPROPOSAL(1).EQ.6)THEN
        CALL READMPICONT(CONT,VETCONT,NPROCS,NRANK,NSINCRO)
        CONTMIN = INTMIN(VETCONT,NPROCS)
     ELSE
        CONTMIN = CONT
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL DREAMTYPEMETHODS(CONTMIN,LOGERROR_LOCAL,&
          LOGERROR,NRANK,INDX,NPROCS,PCRLOCAL,NSINCRO)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(CONTADORC.GT.NRTOTAL)THEN
        IF(NRANK.EQ.0)WRITE(*,*)'#################################'
        IF(NRANK.EQ.0)WRITE(*,*)'### SAIDA SEM CONVERGENCIA ######'
        IF(NRANK.EQ.0)WRITE(*,*)'#################################'
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        GO TO 500
     ENDIF
!
     IF(NRANK.EQ.NRANK)WRITE(*,4444)CONTADORC,NRANK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CRIACAO DO INPUT PARA OS CAMPOS !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GERACAO DO CAMPO ALEATORIO !!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO K=1,NPRIORR
        CALL GERAFILEIN(NERROREAD,NRANK,NPROCS,K)
        CALL GERADOR(CONTADORC,NRANK,K)
        IF(NERROREAD.NE.1)THEN
           IF(NRANK.EQ.0)WRITE(*,*)'######### PROBLEMA NA ###########'
           IF(NRANK.EQ.0)WRITE(*,*)'#### CRIACAO DAS ENTRADAS DO ####'
           IF(NRANK.EQ.0)WRITE(*,*)'######  GERADOR DE CAMPOS  ######'
           CALL MPI_FINALIZE(IERROR)
           STOP 600
        END IF
     END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULO DO PRIOR RATIO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL PRIOR_RATIO()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! TWO STAGE ALGORITHM !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(NSTAGE.EQ.2)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  COARSE PROBLEM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPSCALING DO CAMPO GERADO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO K=1,NPRIORR
           CALL UPSCALING(K,NRANK)
        END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESOLUCAO DO PROBLEMA DE TRANSPORTE (COARSE SCALE) !!!
        CALL SIMULADOR_C(NRANK)
        IF(NSIMUL.EQ.5.AND.ERROSIMULADOR.EQ.1)THEN
           WRITE(*,1001)NRANK
           GOTO 999
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DA AMOSTRA !!!!!!!!!!!!!!!!!!!!!!!!!
        CALL READ_AMOSTRA_C(NRANK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULO DO LIKELIHOOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL LIKELIHOOD_C()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROCESSO DE ACEITACAO-REJEICAO !!!!!!!!!!!!!!!!!!!!!!!
        CALL MCMC_C(CONTADORC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ELSE
        ACCEPT = 1
     END IF
!!#####################################################!!
!! END OF FIRST STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FINE SCALE PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(ACCEPT.GT.0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESOLUCAO DO PROBLEMA DE TRANSPORTE (FINE SCALE)!!!!!!
        CALL SIMULADOR_F(NRANK)
        IF(NSIMUL.EQ.5.AND.ERROSIMULADOR.EQ.1)THEN
           WRITE(*,1000)NRANK
           GOTO 999
        END IF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DA AMOSTRA (FINE SCALE)!!!!!!!!!!!!!
        CALL READ_AMOSTRA_F(NRANK)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULO DO LIKELIHOOD (FINE SCALE)!!!!!!!!!!!!!!!!!!!!
        CALL LIKELIHOOD_F()
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROCESSO DE ACEITACAO-REJEICAO !!!!!!!!!!!!!!!!!!!!!!!
        CALL MCMC_F(CONTADORC,NRANK)
!!
        CONTADORF=CONTADORF+1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COPIA DOS CAMPOS ACEITOS E REJEITADOS !!!!!!!!!!!!!!!!
     CALL COPYFILES(CONTADORC,NAMEPRIOR,NRANK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(ACCEPT.GT.0)THEN
        CONT=CONT+1
        IF(NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5.OR.&
             NPROPOSAL(1).EQ.6)THEN
           VETCONT(NRANK+1)=CONT
           CALL SAVEMPICONT(CONT,1,NRANK,NSINCRO)
        END IF
        WRITE(*,555)NRANK,JUMP
     ELSE
        CONTREJ=CONTREJ+1
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END DO ! LOOP PRINCIPAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
500 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,332)NRANK,REAL(CONT)/REAL(CONTADORF)*100.d0
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(ALLOCATED(INDX))DEALLOCATE(INDX)
  IF(ALLOCATED(LOGERROR))DEALLOCATE(LOGERROR)
  IF(ALLOCATED(LOGERRORCF))DEALLOCATE(LOGERRORCF)
  IF(ALLOCATED(NID))DEALLOCATE(NID)
  IF(ALLOCATED(NIDALL))DEALLOCATE(NIDALL)
  IF(ALLOCATED(JVET))DEALLOCATE(JVET)
  IF(ALLOCATED(JVETALL))DEALLOCATE(JVETALL)
  IF(ALLOCATED(PCRALL))DEALLOCATE(PCRALL)
  IF(ALLOCATED(VETCONT))DEALLOCATE(VETCONT)
  IF(ALLOCATED(VET))DEALLOCATE(VET)
  IF(ALLOCATED(DIMX))DEALLOCATE(DIMX)
  IF(ALLOCATED(DIMY))DEALLOCATE(DIMY)
  IF(ALLOCATED(DIMZ))DEALLOCATE(DIMZ)
  IF(ALLOCATED(NX))DEALLOCATE(NX)
  IF(ALLOCATED(NY))DEALLOCATE(NY)
  IF(ALLOCATED(NZ))DEALLOCATE(NZ)
  IF(ALLOCATED(NFLAG))DEALLOCATE(NFLAG)
  IF(ALLOCATED(MKL))DEALLOCATE(MKL)
  IF(ALLOCATED(NTOTAL))DEALLOCATE(NTOTAL)
  IF(ALLOCATED(NCOND))DEALLOCATE(NCOND)
  IF(ALLOCATED(NPROPOSAL))DEALLOCATE(NPROPOSAL)
  IF(ALLOCATED(ALPHA))DEALLOCATE(ALPHA)
  IF(ALLOCATED(VARIAN))DEALLOCATE(VARIAN)
  IF(ALLOCATED(SIG))DEALLOCATE(SIG)
  IF(ALLOCATED(FBETA))DEALLOCATE(FBETA)
  IF(ALLOCATED(NWELLS))DEALLOCATE(NWELLS)
  IF(ALLOCATED(LIKETIPO))DEALLOCATE(LIKETIPO)
  IF(ALLOCATED(NDADOS))DEALLOCATE(NDADOS)
  IF(ALLOCATED(NDADOSI))DEALLOCATE(NDADOSI)
  IF(ALLOCATED(GERATIPO))DEALLOCATE(GERATIPO)
  IF(ALLOCATED(REFDATA))DEALLOCATE(REFDATA)
  IF(ALLOCATED(REFAMOS))DEALLOCATE(REFAMOS)
  IF(ALLOCATED(THETA))DEALLOCATE(THETA)
  IF(ALLOCATED(THETAN))DEALLOCATE(THETAN)
  IF(ALLOCATED(MAXDATA))DEALLOCATE(MAXDATA)
  IF(ALLOCATED(NORMA_L2_REF))DEALLOCATE(NORMA_L2_REF)
  IF(ALLOCATED(NLPRIOR))DEALLOCATE(NLPRIOR)
  IF(ALLOCATED(NUM_AVET))DEALLOCATE(NUM_AVET)
  IF(ALLOCATED(NUM_AVETOLD))DEALLOCATE(NUM_AVETOLD)
  IF(ALLOCATED(XRMVN_MEAN))DEALLOCATE(XRMVN_MEAN)
  IF(ALLOCATED(XRMVN))DEALLOCATE(XRMVN)
  IF(ALLOCATED(XRMVN_VET))DEALLOCATE(XRMVN_VET)
  IF(ALLOCATED(COVMATUP))DEALLOCATE(COVMATUP)
  IF(ALLOCATED(COVMATDOWN))DEALLOCATE(COVMATDOWN)
  IF(ALLOCATED(COVMATID))DEALLOCATE(COVMATID)
  IF(ALLOCATED(PRATIO))DEALLOCATE(PRATIO)
  IF(ALLOCATED(LRATIOF))DEALLOCATE(LRATIOF)
  IF(ALLOCATED(LRATIOC))DEALLOCATE(LRATIOC)
  IF(ALLOCATED(ERRORCC))DEALLOCATE(ERRORCC)
  IF(ALLOCATED(ERRORNC))DEALLOCATE(ERRORNC)
  IF(ALLOCATED(ERRORCF))DEALLOCATE(ERRORCF)
  IF(ALLOCATED(ERRORNF))DEALLOCATE(ERRORNF)
  IF(ALLOCATED(VARPROP))DEALLOCATE(VARPROP)
  IF(ALLOCATED(NFREQ))DEALLOCATE(NFREQ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL MPI_FINALIZE(IERROR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  stop
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
200  FORMAT(/,'#####################################################',/,&
              'NUMERO DE CADEIAS...............................: ',I3,/,&
              '#####################################################',/)
201  FORMAT(/,'#####################################################',/,&
              '### NOME DO PROCESSADOR',I4, ':.....: ',A,/,&
              '#####################################################',/)
1313 FORMAT(/,'#####################################################',/,&
              '#### PROBLEMA NA INICIALIZACAO DO PROTOCOLO MPI #####',/,&
              '#####################################################',/)
300  FORMAT(/,'#####################################################',/,&
              'PROCESSADOR NUMERO..............................: ',I3,/,&
              'VALOR DE NSIMUL.................................: ',I3,/,&
              '#####################################################',/)
4444 FORMAT(/,'#####################################################',/,&
              '## REALIZACAO....: ',I7,'      #####################',/,&
              '## RANK..........: ',I7,'      #####################',/,&
              '#####################################################')

331  FORMAT('TIME TO GENERATE RANDOM FIELD (s) =',F8.2,/,&
         'USER TIME           (s) =',F8.2,/,&
         'SYSTEM TIME         (s) =',F8.2,/)
332  FORMAT('================================',/,&
         'RANK................. =',I7,/,&
         'TAXA DE ACEITACAO (%) =',F8.2,/,&
         '================================') 
333  FORMAT('TIME TO SIMULATION  (s) =',F8.2,/,&
         'USER TIME           (s) =',F8.2,/,&
         'SYSTEM TIME         (s) =',F8.2,/)
334  FORMAT(I7,2X,E15.7)
1000 FORMAT('((----- REFAZENDO SIMULACAO FINE NO RANK: ',I2,'-----))')
1001 FORMAT('((----- REFAZENDO SIMULACAO COARSE NO RANK: ',I2,'-----))')
555  FORMAT('ACCEPTED_JUMP AT RANK',I3,': ',F8.2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!#######################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DE ENTRADA !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_INPUT(NERROREAD,MPI_NSIMU,NPROCS,NK)
!
  USE VARIAVEIS
  USE RANDOM
  IMPLICIT NONE
!
  INTEGER              :: IN_FILE,ISTAT,I,N,J
  INTEGER              :: MPI_NSIMU,NPROCS,NK,DIM
  CHARACTER(LEN=128)   :: FILEINPUT,FILEOUT
  INTEGER, INTENT(OUT) :: NERROREAD
  LOGICAL              :: TFILE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! OPEN INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FILEINPUT = ADJUSTL(TRIM('./in/entra.in'))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IN_FILE = 100+NK
  INQUIRE(FILE=FILEINPUT,EXIST=TFILE)
!
  OPEN(UNIT=IN_FILE,FILE=FILEINPUT,ACTION='READ',&
       STATUS='OLD',FORM='FORMATTED',IOSTAT=ISTAT)
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',ADJUSTL(TRIM(FILEINPUT))
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ INPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  IF(NK.EQ.0)WRITE(*,*)'############ READING INPUT DATA FROM: ##############'
  IF(NK.EQ.0)WRITE(*,*)ADJUSTL(TRIM(FILEINPUT))
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  IF(NK.EQ.0)WRITE(*,*)''
!
  READ(IN_FILE,100)NSTAGE
  IF(NSTAGE.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,108)NSTAGE
  ELSE
     IF(NSTAGE.EQ.2)THEN
        IF(NK.EQ.0)WRITE(*,109)NSTAGE
     ELSE
        IF(NK.EQ.0)WRITE(*,110)
        STOP
     END IF
  END IF
!
  READ(IN_FILE,100)NPRIORR
  IF(NK.EQ.0)WRITE(*,98)NPRIORR
!
  DO I=1,NPRIORR
     READ(IN_FILE,200)NAMEPRIOR(I)
     FILEAM(I)=ADJUSTL(TRIM(NAMEPRIOR(I)))
     IF(NK.EQ.0)WRITE(*,206)I,ADJUSTL(TRIM(FILEAM(I)))
  END DO
  !
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  IF(NSTAGE.GT.0)THEN
     ALLOCATE(NUPSC(NPRIORR))
     READ(IN_FILE,94)(NUPSC(I),I=1,NPRIORR)
     IF(NK.EQ.0)WRITE(*,95)(NUPSC(I),I=1,NPRIORR)
  END IF
!
  READ(IN_FILE,100)NRTOTAL
  IF(NK.EQ.0)WRITE(*,101)NRTOTAL
!
  READ(IN_FILE,100)NR
  IF(NK.EQ.0)WRITE(*,102)NR
!
  READ(IN_FILE,100)NSIMUL
  IF(NK.EQ.0)THEN
     WRITE(*,107)NSIMUL
  END IF
  MPI_NSIMU = NSIMUL
!
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  IF(NK.EQ.0)WRITE(*,*)'DEFINICAO DOS DADOS OBSERVADOS'
  IF(NK.EQ.0)WRITE(*,*)'####################################################'
  READ(IN_FILE,93)NDTYPE,NNORMA
  IF(NK.EQ.0)WRITE(*,103)NDTYPE
  IF(NNORMA.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,123)NNORMA
  ELSE
     IF(NK.EQ.0)WRITE(*,133)NNORMA
  END IF
  !
  ALLOCATE(LOGERRORCF(NRTOTAL))
  ALLOCATE(VETCONT(NPROCS))
  VETCONT=0D0
  ALLOCATE(NORMA_L2_REF(NDTYPE))
  ALLOCATE(MAXDATA(NDTYPE))
  ALLOCATE(NWELLS(NDTYPE))
  ALLOCATE(LIKETIPO(NDTYPE))
  ALLOCATE(NDADOS(NDTYPE))
  ALLOCATE(NDADOSI(NDTYPE))
  ALLOCATE(GERATIPO(NPRIORR))
  ALLOCATE(NPROPOSAL(NPRIORR))
  ALLOCATE(NLPRIOR(NPRIORR))
  ALLOCATE(NUM_AVET(NPRIORR))
  ALLOCATE(NUM_AVETOLD(NPRIORR))
  ALLOCATE(NEWLBD(NPRIORR))
  ALLOCATE(OLDLBD(NPRIORR))
  ALLOCATE(YNEWLBD(NPRIORR))
  ALLOCATE(YOLDLBD(NPRIORR))
  ALLOCATE(NKIND(NPRIORR))
  ALLOCATE(NZMIN(NPRIORR))
  ALLOCATE(ZMAX(NPRIORR))
  ALLOCATE(LBDFIXO(NPRIORR))
  ALLOCATE(XLBDMEAN(NPRIORR))
  ALLOCATE(XLBDSD(NPRIORR))
  ALLOCATE(XLBDBETA(NPRIORR))
  ALLOCATE(NLBD(NPRIORR))
  ALLOCATE(SIG(NPRIORR))
  ALLOCATE(MKL(NPRIORR))
  ALLOCATE(FILE_OUTEXP(NPRIORR,NPROCS))
  ALLOCATE(VARPROP(NPRIORR))
  ALLOCATE(NFREQ(NPRIORR))
!
  NLBD=1
!
  DO I=1,NDTYPE
     READ(IN_FILE,300)NWELLS(I),LIKETIPO(I),NDADOSI(I),NDADOS(I)
     IF(NK.EQ.0)WRITE(*,*)'-----------------------------------------'
     IF(NK.EQ.0)WRITE(*,*)'-----------------------------------------'
     IF(NK.EQ.0)WRITE(*,*)'### CONJUNTO ',I,' DE DADOS'
     IF(NK.EQ.0)WRITE(*,301)NWELLS(I),LIKETIPO(I),NDADOS(I)-NDADOSI(I)+1
     IF(LIKETIPO(I).EQ.1)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< DADOS DE PRODUCAO >>>'
     END IF
     IF(LIKETIPO(I).EQ.2)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< DADOS DE CONCENTRACAO - NORMA L^2 >>>'
     END IF
     IF(LIKETIPO(I).EQ.3)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< DADOS DE CONCENTRACAO - NORMA Inf >>>'
     END IF
  END DO
!
  DO I=1,NPRIORR
     READ(IN_FILE,92)GERATIPO(I),NPROPOSAL(I),NLPRIOR(I),MKL(I),VARPROP(I),NFREQ(I),SIG(I)
     IF(NK.EQ.0)WRITE(*,151)I,GERATIPO(I),NPROPOSAL(I),MKL(I),VARPROP(I),NFREQ(I),SIG(I)
     IF(NPROPOSAL(I).EQ.2)THEN
        IF(NK.EQ.0)WRITE(*,152)
     END IF
     IF(NPROPOSAL(I).EQ.3)THEN
        IF(NK.EQ.0)WRITE(*,153)
     END IF
     IF(NPROPOSAL(I).EQ.4)THEN
        IF(NK.EQ.0)WRITE(*,154)
     END IF
     IF(NPROPOSAL(I).EQ.5)THEN
        IF(NK.EQ.0)WRITE(*,155)
     END IF
     IF(NPROPOSAL(I).EQ.6)THEN
        IF(NK.EQ.0)WRITE(*,155)
     END IF
     IF(NLPRIOR(I).EQ.0)THEN
        IF(NK.EQ.0)WRITE(*,600)
     ELSE
        IF(NK.EQ.0)WRITE(*,601)NLPRIOR(I)
        ALLOCATE(CLBDS(NLPRIOR(I)))
        DO J=1,NLPRIOR(I)
           READ(IN_FILE,200)NAME_AVET(I,J)
           IF(NK.EQ.0)WRITE(*,602)J,ADJUSTL(TRIM(NAME_AVET(I,J)))
        END DO
     END IF
     IF(GERATIPO(I).EQ.33)THEN
        READ(IN_FILE,100)TIPOVAR
        IF(TIPOVAR.EQ.0)THEN
           IF(NK.EQ.0)WRITE(*,700)TIPOVAR
        ELSE
           IF(NK.EQ.0)WRITE(*,701)TIPOVAR
        END IF
     END IF
     IF(VARPROP(I).EQ.0)THEN
        IF(NK.EQ.0)WRITE(*,800)
     ELSE
        IF(VARPROP(I).EQ.1)THEN
           IF(NK.EQ.0)WRITE(*,801)
        ELSE
           WRITE(*,*)'PROBLEMA NA ESCOLHA DO TIPO DE JUMP'
           STOP 'JUMP'
        END IF
     END IF
  END DO
!
  IF(NPROPOSAL(1).EQ.3)THEN
     READ(IN_FILE,302)NSTART,NFREQAM,NPCHAIN
     IF(NK.EQ.0)WRITE(*,303)NSTART,NFREQAM,NPCHAIN
  END IF
!
  IF(NPROPOSAL(1).EQ.5.OR.NPROPOSAL(1).EQ.6)THEN
     READ(IN_FILE,920)NDELTA,NCR,CR_DREAM
     IF(NK.EQ.0)WRITE(*,921)NDELTA,NCR,CR_DREAM
  END IF
!
  ALLOCATE(SIGMA2F(NDTYPE))
  DO I=1,NDTYPE
     READ(IN_FILE,400)SIGMA2F(I)
     IF(NK.EQ.0)WRITE(*,401)I,SIGMA2F(I)
  END DO
!
  ALLOCATE(SIGMA2C(NDTYPE))
  DO I=1,NDTYPE
     READ(IN_FILE,400)SIGMA2C(I)
     IF(NK.EQ.0)WRITE(*,402)I,SIGMA2C(I)
  END DO
!
  DO I=1,NDTYPE
     READ(IN_FILE,200)FILEREF(I)
     FILEREF(I)=ADJUSTL(TRIM(FILEREF(I)))
     IF(NK.EQ.0)WRITE(*,202)ADJUSTL(TRIM(FILEREF(I)))
!
     READ(IN_FILE,200)FILEAM(I)
     FILEAM(I)=ADJUSTL(TRIM(FILEAM(I)))
     IF(NK.EQ.0)WRITE(*,203)ADJUSTL(TRIM(FILEAM(I)))
  END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  READ(IN_FILE,200)NAME_OUT
  NAME_OUT=ADJUSTL(TRIM(NAME_OUT))
  IF(NK.EQ.0)WRITE(*,205)ADJUSTL(TRIM(NAME_OUT))
!
  FILERROR=TRIM('./error/erros_')//TRIM(NAME_OUT)
!  READ(IN_FILE,200)FILERROR
  FILERROR=ADJUSTL(TRIM(FILERROR))
  IF(NK.EQ.0)WRITE(*,204)ADJUSTL(TRIM(FILERROR))
!
  CALL RANDOM_SEED(SIZE=NS)
  ALLOCATE(NSEED(NS))
  READ(IN_FILE,97)(NSEED(I),I=1,NS)
  IF(NK.EQ.0)  CALL RANDOM_SEED(PUT=NSEED)
  IF(NK.EQ.0)  WRITE(*,96)NK,NS,(NSEED(I),I=1,NS)
!
  READ(IN_FILE,100)NPRINT
  IF(NK.EQ.0)WRITE(*,106)NPRINT
  IF(NPRINT.EQ.0)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< NAO SALVA CAMPOS REJEITADOS >>>'
  END IF
  IF(NPRINT.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<<  SALVANDO CAMPOS REJEITADOS >>>'
  END IF
!
  CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD=1D0
!
700 FORMAT('COMPRIMENTO DE CORREL. TEMPORAL ~ NORM.. = ',I10)
701 FORMAT('COMPRIMENTO DE CORREL. TEMPORAL ~ UNIF.. = ',I10)
92  FORMAT(6I10,E12.5)
920 FORMAT(2I10,E12.5)
921 FORMAT(/,&
         '########### PARAMETROS DO METODO DREAM ############',/,&
         'NUMERO DE PARES USADOS NAS PROPOSTAS.... = ',I10,/,&
         'NCR..................................... = ',I10,/,&
         'PROBABILIDADE DE TROCAS DOS ELEMENTOS DA   ',/,&
         'PROPOSTA (PARAMETRO CR DO MET. DREAM)... = ',E12.5,/,&
         '###################################################',/)

93 FORMAT(2I10)
94 FORMAT(8I10)
95 FORMAT('TIPO DE UPSCALING....................... = ',8I10)
96 FORMAT('############################',/,&
          'MAQUINA................: ',I4,/,&
          'NUMERO DE SEMENTES.....: ',I3,/,&
          'SEMENTES...............: ',40I10)
97 FORMAT(40I12)
98 FORMAT('NUMERO DE PRIORS (VARIAVEIS DE ENTRADA). = ',I10)
100 FORMAT(I10)
101 FORMAT('NUMERO MAXIMO DE SIMULACOES............. = ',I10)
102 FORMAT('NUMERO MAXIMO DE CAMPOS SELECIONADOS.... = ',I10)
103 FORMAT('QUANTIDADE DE TIPOS DE DADOS OBSERVADOS. = ',I10)
123 FORMAT('   DADOS SERAO NORMALIZADOS PARA O LIKE. = ',I10)
133 FORMAT('   NAO NORMALIZA OS DADOS NO LIKELIHOOD. = ',I10)
151 FORMAT(/,'-----------------------------------------',/,&
           'DEFINE O GERADOR DOS CAMPOS ALEATORIOS',I2,' = ',I10,/,&
           'OPCAO DE PROPOSAL........................= ',I10,/,&
           'DIMENSAO ESTOCASTICA (KL)................= ',I10,/,&
           'JUMP FIXO (0) ou JUMP VARIAVEL (1).......= ',I10,/,&
           'FREQUENCIA DE AJUSTE DO JUMP = 1.........= ',I10,/,&
           'PARAMETRO DO RANDOM WALK.................= ',E12.5)
152 FORMAT('-----------------------------------------',/,&
         '---------------RANDON WALK---------------',/,&
         '-----------------------------------------')
153 FORMAT('-----------------------------------------',/,&
         '----------------METODO AM----------------',/,&
         '-----------------------------------------')
154 FORMAT('-----------------------------------------',/,&
         '----------------METODO DE----------------',/,&
         '-----------------------------------------')
155 FORMAT('-----------------------------------------',/,&
         '--------------METODO DREAM---------------',/,&
         '-----------------------------------------')
600 FORMAT('SEM MUDANCA NO COMPRIMENTO DE CORRELACAO')
601 FORMAT('LEITURA DAS',I2,' MATRIZES KL')
602 FORMAT('ARQUIVO MATRIZ T.',I2,'...: ',A)
106 FORMAT('SALVAR CAMPOS REJEITADOS? (1)=SIM E (0)=NAO',I10)
107 FORMAT('DEFINE O SIMULADOR QUE SERA USADO....... = ',I10)
108 FORMAT('SINGLE STAGE............................ = ',I10)
109 FORMAT('TWO STAGES.............................. = ',I10)
110 FORMAT('##### ERROR #####',/,&
         'DEFINA O TIPO DE ALGORITMO: (1) UM ESTAGIO (2) DOIS ESTAGIOS')
200 FORMAT(A)
202 FORMAT('NOME DO ARQUIVO DE DADOS REFERENCIA..... = ',A)
203 FORMAT('NOME DO ARQUIVO DE DADOS DA AMOSTRA..... = ',A)
204 FORMAT('NOME DO ARQUIVO DE SAIDA DOS ERROS...... = ',A)
205 FORMAT('NOME BASE DOS ARQUIVOS DE SAIDA......... = ',A)
206 FORMAT('NOME BASE DO ARQUIVO DO PRIOR ',I2,' ....... = ',A)
300 FORMAT(4(I10))
302 FORMAT(3(I10))
301 FORMAT('NUMERO DE POCOS......................... = ',I10,/,&
           'TIPO DE LIKELIHOOD USADO................ = ',I10,/,&
           'NUMERO DE AVALIACOES NO TEMPO........... = ',I10)
303 FORMAT(/,'######## PARAMETROS DO METODO AM ##########',/,&
           'INICIO DA ADAPTACAO..................... = ',I10,/,&
           'FREQUENCIA DA ADAPTACAO................. = ',I10,/,&
           'NUMERO DE ELEMENTOS DA CADEIA USADOS.... = ',I10,/,&
           '###########################################',/)
400 FORMAT(E12.5)
401 FORMAT('S^2 FINE SCALE........................',I2,' = ',E12.5)
402 FORMAT('S^2 COARSE SCALE......................',I2,' = ',E12.5)
500 FORMAT(4I10,2e12.5)
800 FORMAT('JUMP FIXO')
801 FORMAT('JUMP VARIAVEL')
!
END SUBROUTINE READ_INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INICIALIZACAO DO PROGRAMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIAL(NERROREAD,NK,NPROC,TFILE)
  !
  USE MPI
  USE RANDOM
  USE VARIAVEIS
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: FOUT
  CHARACTER(LEN=5)   :: CHARAC
  CHARACTER(LEN=4)   :: NUMB
  INTEGER :: NERROREAD,I,ISTAT,NSINAL,NSINAL2,NSINALAM
  INTEGER :: OUT_FILE,OUT_FILEL,K,IOs,DIM,N,LOC,J,NK,NPROC
  LOGICAL :: TFILE,TFILEB
  INTEGER :: TAG,IERROR,clock
  INTEGER :: STATUS(MPI_STATUS_SIZE)
  INTEGER :: ROOT
  INTEGER, DIMENSION(3) :: NVAL
  INTEGER, DIMENSION(NS):: NKSEED
  REAL,    DIMENSION(NS):: SEED
  INTEGER, PARAMETER    :: IA=843314861,IB=453816693,M=1073741824
  INTEGER :: IMG
  DATA ROOT/0/
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  ALLOCATE(ERRORCC(NDTYPE))
  ALLOCATE(ERRORNC(NDTYPE))
  ALLOCATE(ERRORCF(NDTYPE))
  ALLOCATE(ERRORNF(NDTYPE))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LEITURA DO ARQUIVO DE STATUS INICIAL !!!!!!!!!!!!!!!!!!
  SGN = 1
  FILEINI='./in/init_stat_'
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)     
  FILEINI= ADJUSTL(TRIM(FILEINI)//TRIM(NAME_OUT)//('_RK')//TRIM(CHARAC)//TRIM('.in'))
  OUT_FILE = 200+NK
  OUT_FILEL= 300+NK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INQUIRE(FILE=FILEINI,EXIST=TFILE)
  IF(TFILE)THEN
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     IF(NK.EQ.0)WRITE(*,*)'CONTINUAR EXPERIMENTO ANTERIOR'
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     OPEN(UNIT=OUT_FILE,FILE=FILEINI,STATUS='UNKNOWN',&
          ACTION='READ',IOSTAT=ISTAT)
     READ(UNIT=OUT_FILE,FMT=500)CONTADORC,CONTADORF,CONT,NCHAIN
     READ(UNIT=OUT_FILE,FMT=550)(ERRORNF(I),I=1,NDTYPE)
     READ(UNIT=OUT_FILE,FMT=550)(ERRORNC(I),I=1,NDTYPE)
     READ(UNIT=OUT_FILE,FMT=503)(NSEED(I),I=1,NS)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL RANDOM_SEED(put=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FOUT = TRIM('./out/lambda_')//TRIM(NAME_OUT)//TRIM('.dat')
     FOUT = ADJUSTL(TRIM(FOUT))
     INQUIRE(FILE=FOUT,EXIST=TFILEB)
     IF(TFILEB)THEN
        OPEN(UNIT=OUT_FILEL,FILE=FOUT,STATUS='UNKNOWN',&
             ACTION='READ',IOSTAT=ISTAT)
        READ(OUT_FILEL,505,IOSTAT=IOs)NVAL(1),NVAL(2),NVAL(3),OLDLBD(NVAL(2))
        IF(IOs>0)THEN
           IF(NK.EQ.0)WRITE(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
           IF(NK.EQ.0)WRITE(*,*)"@    ALGO ERRADO NA LEITURA      @"
           IF(NK.EQ.0)WRITE(*,*)"@ PARA A CONTINUACAO DOS LAMBDAS @"
           IF(NK.EQ.0)WRITE(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
           STOP 342
        ELSE IF(IOs<0)THEN
           IF(NK.EQ.0)WRITE(*,*)"ULTIMO LAMBDA:",OLDLBD(NVAL(2))
           STOP 243
        ELSE
           NUM_AVETOLD(NVAL(2))=NVAL(3)
        END IF
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE
!! SEMENTES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%% INICIANDO NOVO EXPERIMENTO %%%%%%%%'
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
     IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!
     NKSEED = NSEED
     TAG = 0
     IF(NK.EQ.0)THEN
        CALL RANDOM_SEED(put=NSEED)
        DO I=1,NPROC-1
           CALL RANDOM_NUMBER(SEED)
           NKSEED = FLOOR(1E+08*SEED)
           NKSEED = IB +IA*NKSEED
           DO K=1,NS
              IF(NKSEED(K).LT.0) NKSEED(K) = (NKSEED(K)+M)+M
           END DO
           CALL MPI_SEND(NKSEED,NS,MPI_INTEGER,I,TAG,MPI_COMM_WORLD,IERROR)
        END DO
     ELSE
        CALL MPI_RECV(NSEED,NS,MPI_INTEGER,ROOT,TAG,MPI_COMM_WORLD,&
             STATUS,IERROR)
        CALL RANDOM_SEED(put=NSEED)
     END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CONT     = 0
     CONTREJ  = 0
     CONTADORC= 0d0
     CONTADORF= 0D0
     TERRORNC = 1.D0
     TERRORNF = 1.D0
     ERRORNC  = 1.D0
     ERRORNF  = 1.D0
     ERRORCC  = 1.D0
     ERRORCF  = 1.D0
     NCHAIN   = 1
     NUM_AVET = 1
  END IF
  IF(NK.EQ.0)WRITE(*,FMT=501)CONTADORC,CONTADORF,CONT,NCHAIN,&
       TERRORNC,TERRORNF
  WRITE(*,FMT=502)NK,NS,(NSEED(I),I=1,NS)
!
  CLOSE(OUT_FILE)        
  CLOSE(OUT_FILEL)        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALOCACAO DE MEMORIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  NCONDMAX  = 20
  ALLOCATE(VET(NPRIORR,0:3,0:NCONDMAX-1))
  ALLOCATE(DIMX(NPRIORR))
  ALLOCATE(DIMY(NPRIORR))
  ALLOCATE(DIMZ(NPRIORR))
  ALLOCATE(NX(NPRIORR))
  ALLOCATE(NY(NPRIORR))
  ALLOCATE(NZ(NPRIORR))
  ALLOCATE(NFLAG(NPRIORR))
  ALLOCATE(NTOTAL(NPRIORR))
  ALLOCATE(NCOND(NPRIORR))
  ALLOCATE(ALPHA(NPRIORR))
  ALLOCATE(VARIAN(NPRIORR))
  ALLOCATE(FBETA(NPRIORR))
  ALLOCATE(PRATIO(NPRIORR))
  ALLOCATE(LRATIOF(NPRIORR))
  ALLOCATE(LRATIOC(NPRIORR))
!
  NTOTAL=0D0
  NCOND=0D0
  NX=0D0
  NY=0D0
  NZ=0D0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEFINICAO DO TIPO DE GERADOR !!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NK.EQ.0)WRITE(*,*)'DEFINICAO DO GERADOR DE CAMPOS ALEATORIOS'
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!
  CALL  COPY_DIRGER(NK,NSIMUL,NSTAGE)
!  
  DO I=1,NPRIORR
     IF(NK.EQ.0)WRITE(*,*)'|||||||||||||||||||||||||||||||||||||||||'
     IF(NK.EQ.0)WRITE(*,*)'|||||||||||||||||||||||||||||||||||||||||'
     IF(NK.EQ.0)WRITE(*,*)'|||||||   VARIAVEL:',I,        '  |||||||'
     IF(NK.EQ.0)WRITE(*,*)'|||||||||||||||||||||||||||||||||||||||||'
     IF(GERATIPO(I).EQ.0)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR LABTRANGEO (RANDOM WALK) 3 >>>'
        FILE_INFIELD(I) = '../gera_LABTRANGEORW3/in'
     END IF
     IF(GERATIPO(I).EQ.1)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN/in'
     END IF
     IF(GERATIPO(I).EQ.2)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR LABTRANGEO >>>'
        FILE_INFIELD(I) = '../gera_LABTRANGEO/in'
     END IF
     IF(GERATIPO(I).EQ.3)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN (RANDOM WALK) >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_RW/in'
     END IF
     IF(GERATIPO(I).EQ.31)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN (RANDOM WALK) >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_RW1/in'
     END IF
     IF(GERATIPO(I).EQ.34)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN (RANDOM WALK) >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_RW2/in'
     END IF
     IF(GERATIPO(I).EQ.32)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN 3D (RANDOM WALK) >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_RW3D/in'
     END IF
     IF(GERATIPO(I).EQ.33)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR KL_FORTRAN 3D (RANDOM WALK) >>>'
        IF(NK.EQ.0)WRITE(*,*)'<<<         DIFERENTES LAMBDAS          >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_LBD/in'
     END IF
     IF(GERATIPO(I).EQ.4)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR LABTRANGEO (RANDOM WALK) >>>'
        FILE_INFIELD(I) = '../gera_LABTRANGEORW/in'
     END IF
     IF(GERATIPO(I).EQ.5)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR LABTRANGEO (RANDOM WALK) 2 >>>'
        FILE_INFIELD(I) = '../gera_LABTRANGEORW2/in'
     END IF
     IF(GERATIPO(I).EQ.6)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN/in'
     END IF
     IF(GERATIPO(I).EQ.7)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN1/in'
     END IF
     IF(GERATIPO(I).EQ.8)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN2/in'
     END IF
     IF(GERATIPO(I).EQ.13)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3D/in'
     END IF
     IF(GERATIPO(I).EQ.14)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3D_2/in'
     END IF
!
     WRITE(NUMB,'(I4.3)')NK
     FILE_INFIELD(I)=ADJUSTL(TRIM(FILE_INFIELD(I)))//TRIM(ADJUSTL(NUMB))//&
          ('/entrada.in')
!
     IF(NK.EQ.0)WRITE(*,504)FILE_INFIELD(I)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT DATA OF RANDOM FIELDS !!!!!!!!!!!!!!!!!!!!!!
     NERROREAD = 0
     CALL READ_FIELD(NERROREAD,I,NK)
!
     IF(NERROREAD.NE.1)THEN
        IF(NK.EQ.0)WRITE(*,*)'####### ERRO NA LEITURA #######'
        IF(NK.EQ.0)WRITE(*,*)'### DO ARQUIVO INPUT FIELDS ###'
!        STOP
     END IF
!
  END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEFINICAO DO SIMULADOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NK.EQ.0)WRITE(*,*)'SIMULADOR UTILIZADO PARA O PROBLEMA DE ESCOAMENTO'
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NSIMUL.EQ.0)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< BLACKBOX >>>'
  END IF
  IF(NSIMUL.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMULADOR >>>'
  END IF
  IF(NSIMUL.EQ.2)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< simuadorRigido >>>'
  END IF
  IF(NSIMUL.EQ.3)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMUL_COMP >>>'
  END IF
  IF(NSIMUL.EQ.4)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< UW SIMULATOR >>>'
  END IF
  IF(NSIMUL.EQ.5)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMULADOR_VISCOELASTICO >>>'
  END IF
  IF(NSIMUL.EQ.6)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< TWOPHASEFLOW MATLAB >>>'
  END IF
  IF(NSIMUL.EQ.7)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< TWOPHASEFLOW OCTAVE >>>'
  END IF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALOCACAO DE MEMORIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%'
  IF(NK.EQ.0)WRITE(*,*)'ALOCACAO DE MEMORIA'
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%'
!
  ALLOCATE(REFDATA(NDTYPE,MAXVAL(NWELLS)+1,MAXVAL(NDADOS)))
  ALLOCATE(REFAMOS(NDTYPE,MAXVAL(NWELLS)+1,MAXVAL(NDADOS)))
  ALLOCATE(ERRORF(NR))
  REFAMOS=0.0D0
  REFDATA=0.0D0
!
  NSINAL=0
  NSINAL2=0
  DO I=1,NPRIORR
     IF(GERATIPO(I).EQ.1.OR.GERATIPO(I).EQ.3.OR.&
          GERATIPO(I).EQ.0.OR.GERATIPO(I).EQ.31.OR.&
          GERATIPO(I).EQ.32.OR.GERATIPO(I).EQ.33.OR.&
          GERATIPO(I).EQ.34.OR.GERATIPO(I).EQ.6.OR.&
          GERATIPO(I).EQ.7.OR.GERATIPO(I).EQ.8.OR.&
          GERATIPO(I).EQ.13.OR.GERATIPO(I).EQ.14)THEN
        IF(NSINAL.EQ.0)THEN
           ALLOCATE(THETA(NPRIORR,MAX(MAXVAL(MKL),MAXVAL(NTOTAL))))
           ALLOCATE(THETAN(NPRIORR,MAX(MAXVAL(MKL),MAXVAL(NTOTAL))))
           NSINAL=1
        END IF
     ELSE
        IF(NSINAL.EQ.0)THEN
           ALLOCATE(LOGPERM(NPRIORR,MAXVAL(NX)*MAXVAL(NY)))
           ALLOCATE(LOGPERMN(NPRIORR,MAXVAL(NX)*MAXVAL(NY)))
           NSINAL2=1
        END IF
     END IF
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ALOCACAO DE MEMORIA PARA O METODO AM !!!!!!!!!!!!!!!!!!
  IF(NPROPOSAL(1).EQ.3)THEN
     DIM=0
     NSINALAM=0
     DO I=1,NPRIORR
        DIM = MAX(DIM,MKL(I))
        NSINALAM = 1
     END DO
     DIMR   = DIM*(DIM+1)/2
     IF(NSINALAM.GT.0)THEN
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%% METODO AM %%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DO VETOR ALEATORIO: ',DIM
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DA MATRIZ DE COVAR: ',DIMR
     !
        ALLOCATE(XRMVN(DIM))
        ALLOCATE(XRMVN_VET(DIM,NPCHAIN))
        ALLOCATE(XRMVN_MEAN(DIM))
        ALLOCATE(COVMATUP(DIMR))
        ALLOCATE(COVMATDOWN(DIMR))
        ALLOCATE(COVMATID(DIMR))
        COVMATUP   = 0.0E0
        COVMATDOWN = 0.0E0
        COVMATID   = 0.0E0
        XRMVN      = 0.0E0
        XRMVN_VET  = 0.0E0
        XRMVN_MEAN = 0.0E0
     !! MATRIZ IDENTIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,DIM
           DO J=I,DIM
              LOC = J*(J-1)/2+I
              IF(I.EQ.J) COVMATID(LOC) = 1.0E+00
           END DO
        END DO
     END IF
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALOCACAO DE MEMORIA PARA O METODO DE OU DREAM !!!!!!!!!
  IF(NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5.OR.&
       NPROPOSAL(1).EQ.6)THEN
     SINALDREAM = .FALSE.
     DIM=0
     NSINALAM=0
     DO I=1,NPRIORR
        DIM = MAX(DIM,MKL(I))
        NSINALAM = 1
     END DO
     DIMR   = DIM*(DIM+1)/2
     IF(NSINALAM.GT.0)THEN
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)THEN
           IF(NPROPOSAL(1).EQ.4)WRITE(*,*)'%%%%%%%%%%%%%%%%%% METODO DE %%%%%%%%%%%%%%%%%%'
           IF(NPROPOSAL(1).EQ.5)WRITE(*,*)'%%%%%%%%%%%%%%%%%% METODO DREAM %%%%%%%%%%%%%%%'
           IF(NPROPOSAL(1).EQ.6)WRITE(*,*)'%%%%%%%%%%%%%%%%%% METODO DREAM mod %%%%%%%%%%%'
        END IF
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DA MATRIZ MATDE: ',DIM*NPROC
     !
        ALLOCATE(XRMVN(DIM))
        ALLOCATE(MATDE(DIM,NPROC))
        ALLOCATE(COVMATDOWN(DIMR))
        ALLOCATE(COVMATID(DIMR))
        ALLOCATE(XRMVN_MEAN(DIM))
        IF(NPROPOSAL(1).EQ.5.OR.NPROPOSAL(1).EQ.6)THEN
           ALLOCATE(JVET(NPRIORR,NCR))
           JVET = 0.0D0
           ALLOCATE(JVETALL(NPRIORR,NCR))
           JVETALL = 0.0D0
           ALLOCATE(NID(NPRIORR,NCR))
           NID = 0
           ALLOCATE(NIDALL(NPRIORR,NCR))
           NIDALL = 0
           ALLOCATE(PCRALL(NPRIORR,NCR))
           DO I=1,NCR
              DO J=1,NPRIORR
                 PCRALL(J,I)= REAL(1)/REAL(NCR)
              END DO
           END DO
        END IF
        COVMATDOWN = 0.0E0
        COVMATID   = 0.0E0
        XRMVN      = 0.0E0
        XRMVN_MEAN = 0.0E0
        MATDE      = 0.0E0
        SIGDE      = SIG(1)
     !! MATRIZ IDENTIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,DIM
           DO J=I,DIM
              LOC = J*(J-1)/2+I
              IF(I.EQ.J) COVMATID(LOC) = 1.0E+00
           END DO
        END DO
     END IF
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DE REFERENCIA !!!!!!!!!!!!!!!!!!!!!!
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NK.EQ.0)WRITE(*,*)'LEITURA DOS DADOS DE REFERENCIA'
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!
  CALL READ_REF(NK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CRIACAO DOS THETAS INICIAIS !!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(TFILE)THEN
     IF(NK.EQ.0)WRITE(*,*)'/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/'
     IF(NK.EQ.0)WRITE(*,*)'CONTINUACAO (NAO GERA NOVO THETA)'
     IF(NK.EQ.0)WRITE(*,*)'---------------------------------'
  ELSE
     DO K=1,NPRIORR
        IF(GERATIPO(K).EQ.0)THEN
           CALL GERATHETA(NTOTAL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.1)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.3)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.31)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.32)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.34)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.33)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
           NUM_AVET=0
           NUM_AVETOLD=0
        ENDIF
!
        IF(GERATIPO(K).EQ.6)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.7)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.8)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.13)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
        IF(GERATIPO(K).EQ.14)THEN
           CALL GERATHETA(MKL(K),GERATIPO(K),K,NK)
        ENDIF
!
     END DO
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD = 1
  !
112 FORMAT(I5)
505 FORMAT(3I7,E10.3)
504 FORMAT('NOME DO ARQUIVO DE INPUT DOS CAMPOS: ',A)
502 FORMAT('MAQUINA.............: ',I4,/,&
           'NUMERO DE SEMENTES..: ',I4,/,'SEMENTES: ',40I12)
503 FORMAT(40I12)
500 FORMAT(4I10)
550 FORMAT(10E12.5)
501 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,&
           'N. DE TESTES MALHA GROSSA.....=',I10,/,&
           'N. DE TESTES MALHA FINA.......=',I10,/,&
           'N. DE CAMPOS SELECIONADOS.....=',I10,/,&
           'N. DE REPETICOES NA CADEIA....=',I10,/,&
           'ERRO ATUAL NA MALHA FINA......=',E12.5,/,&
           'ERRO ATUAL NA MALHA GROSSA....=',E12.5,/,&
           '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
!
END SUBROUTINE INITIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COPY_DIRGER(NPR,NUMSIMUL,NS)
  !
  USE VARIAVEIS, ONLY: GERATIPO,NPRIORR
  IMPLICIT NONE
!  
  INTEGER            :: NPR,NUMSIMUL,NS,J,TIPO
  CHARACTER (LEN=128):: COMMAND
  CHARACTER (LEN=4)  :: NUMB
!
  write(*,100)NPR
  WRITE(NUMB,'(I4.3)')NPR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO J=1,NPRIORR
     TIPO = GERATIPO(J)
!    IF(TIPO.EQ.0) COMMAND='../gera_LABTRANGEORW3/out/xpsi.dat'
     IF(TIPO.EQ.6)THEN
        COMMAND=('cp -r ../gera_KL/MVN/in/ ')//&
             (' ../gera_KL/MVN/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.7)THEN
        COMMAND=('cp -r ../gera_KL/MVN1/in/ ')//&
             (' ../gera_KL/MVN1/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.8)THEN
        COMMAND=('cp -r ../gera_KL/MVN2/in/ ')//&
             (' ../gera_KL/MVN2/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.1)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN/in/ ')//&
             (' ../gera_KL/FORTRAN/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.3)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW/in/ ')//&
             (' ../gera_KL/FORTRAN_RW/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.13)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_KL3D/in/ ')//&
             (' ../gera_KL/FORTRAN_KL3D/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.14)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_KL3D_2/in/ ')//&
             (' ../gera_KL/FORTRAN_KL3D_2/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.31)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW1/in/ ')//&
             (' ../gera_KL/FORTRAN_RW1/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.34)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW2/in/ ')//&
             (' ../gera_KL/FORTRAN_RW2/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.32)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW3D/in/ ')//&
             (' ../gera_KL/FORTRAN_RW3D/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
     IF(TIPO.EQ.33)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_LBD/in/ ')//&
             (' ../gera_KL/FORTRAN_LBD/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
  100 FORMAT('EXPEREMENT COPY: ',I3)
!
END SUBROUTINE COPY_DIRGER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROTINA PARA SALVAR OS ERROS EM CADA ELEMENTO DA CADEIA !!!!!!!!
SUBROUTINE SAVE_INITSTATUS(N0,N1,N2,N3,NK)
!
  USE RANDOM
  USE VARIAVEIS, only: FILEINI,NSEED,NS
  USE VARIAVEIS, only: ERRORNF,ERRORNC,NDTYPE
  IMPLICIT NONE
!
  INTEGER            :: OUT_FILE,ISTAT
  INTEGER            :: N0,N1,N2,N3,I,NK,J
  REAL               :: ERF,ERC
!
  OUT_FILE = 200+NK
!  WRITE(CHARAC,112)NK
!  CHARAC=ADJUSTL(CHARAC)     
!  FOUT = TRIM(FILEINI)//('_RK')//TRIM(CHARAC)//TRIM('.in')
!
  OPEN(UNIT=OUT_FILE,FILE=FILEINI,STATUS='UNKNOWN',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEINI
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(UNIT=OUT_FILE,FMT=500)N0,N1,N2,N3
  WRITE(UNIT=OUT_FILE,FMT=501)(ERRORNF(I),I=1,NDTYPE)
  WRITE(UNIT=OUT_FILE,FMT=501)(ERRORNC(J),J=1,NDTYPE)
!
  CALL RANDOM_SEED(GET=NSEED)
!  WRITE(*,FMT=101)(NSEED(I),I=1,NS)
  WRITE(UNIT=OUT_FILE,FMT=101)(NSEED(I),I=1,NS)
!
  CLOSE(OUT_FILE)
  RETURN
!
500 FORMAT(4I10)
501 FORMAT(200E12.5)
101 FORMAT(40I12)
112 FORMAT(I5)
!
END SUBROUTINE SAVE_INITSTATUS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ THE INPUT DATA OF RANDOM FIELDS !!!!!!!!!!!!!!!!!!
SUBROUTINE READ_FIELD(NERROREAD,K,NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER :: IN_FILE,IN_FILE2,IN_FILE3,ISTAT,I,J,K
  INTEGER :: NLIXO,NFILE,NFILES,NK,LNGT,MMKL
  INTEGER, DIMENSION(NS) :: NSE
  INTEGER :: INIF,NPROP
  INTEGER, INTENT(OUT) :: NERROREAD
  CHARACTER(LEN=128)   :: NAME,FILENT,FOUT
  CHARACTER(LEN=3)     :: FEXP
  REAL                 :: FDX,FDY,SIGG
!
  IN_FILE = 500+NK
!
  OPEN(UNIT=IN_FILE,FILE=FILE_INFIELD(K),STATUS='OLD',&
       ACTION='READ',FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_INFIELD(K)
     STOP
  END IF
!
  IF(GERATIPO(K).EQ.6.OR.GERATIPO(K).EQ.7.OR.&
       GERATIPO(K).EQ.8)THEN
!
     READ(IN_FILE,3000)MMKL
     IF(NK.EQ.0)WRITE(*,2005)MMKL
!
     READ(IN_FILE,4000)(NSE(I),I=1,NS)
     IF(NK.EQ.0)WRITE(*,4001)(NSE(I),I=1,NS)
!
     READ(IN_FILE,7000)FILE_OUT(K)
     IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
!
     READ(IN_FILE,7000)FILE_VET(K)
     IF(NK.EQ.0)WRITE(*,7003)FILE_VET(K)
!
     READ(IN_FILE,7000)FILE_THE(K)
     IF(NK.EQ.0)WRITE(*,7004)FILE_THE(K)
!
     READ(IN_FILE,12)NPROP,SIGG
     IF(NK.EQ.0)WRITE(*,13)NPROP,SIGG
!     SIG(K) = SIGG
!    
  ELSE
!     
     READ(IN_FILE,8000)INIF,NFILES
     IF(NK.EQ.0)WRITE(*,8001)INIF,NFILES
     NFILE=NFILES-INIF+1
     IF(NK.EQ.0)WRITE(*,8002)NFILE
!
     IF(GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.&
          GERATIPO(K).EQ.13.OR.GERATIPO(K).EQ.14)THEN
        READ(IN_FILE,1000)DIMX(K),DIMY(K),DIMZ(K)
        IF(NK.EQ.0)WRITE(*,8801)DIMX(K),DIMY(K),DIMZ(K)
!
        READ(IN_FILE,2000)NX(K),NY(K),NZ(K)
        IF(NK.EQ.0)WRITE(*,9901)NX(K),NY(K),NZ(K)
!
        FDX = DIMX(K)/NX(K)
        FDY = DIMY(K)/NY(K)
        IF(NK.EQ.0)WRITE(*,5001)FDX,FDY
     ELSE
        READ(IN_FILE,1000)DIMX(K),DIMY(K)
        IF(NK.EQ.0)WRITE(*,1001)DIMX(K),DIMY(K)
!
        READ(IN_FILE,2000)NX(K),NY(K)
        IF(NK.EQ.0)WRITE(*,2001)NX(K),NY(K)
!
        FDX = DIMX(K)/NX(K)
        FDY = DIMY(K)/NY(K)
        IF(NK.EQ.0)WRITE(*,5001)FDX,FDY
     END IF
!
     READ(IN_FILE,3000)NFLAG(K)
!
     IF(NFLAG(K)==2)THEN
        IF(NK.EQ.0)WRITE(*,6002)NFLAG(K)
     END IF
     IF(NFLAG(K)==1)THEN
        IF(NK.EQ.0)WRITE(*,6001)NFLAG(K)
     END IF
!
     READ(IN_FILE,9000)FBETA(K),VARIAN(K)
     IF(NFLAG(K)==1)THEN
        IF(NK.EQ.0)WRITE(*,9001)FBETA(K),VARIAN(K)
     END IF
     IF(NFLAG(K)==2)THEN
        IF(NK.EQ.0)WRITE(*,9002)FBETA(K),VARIAN(K)
     END IF
!
     READ(IN_FILE,4000)(NSE(I),I=1,NS)
     IF(NK.EQ.0)WRITE(*,4001)(NSE(I),I=1,NS)
!
     IF(GERATIPO(K).EQ.1)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2005)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
!
        READ(IN_FILE,7000)FILE_VAL(K)
        IF(NK.EQ.0)WRITE(*,7002)FILE_VAL(K)
!
        READ(IN_FILE,7000)FILE_VET(K)
        IF(NK.EQ.0)WRITE(*,7003)FILE_VET(K)
     ENDIF
!
     IF(GERATIPO(K).EQ.2)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2006)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
     ENDIF
!
     IF(GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31&
          .OR.GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.34&
          .OR.GERATIPO(K).EQ.13.OR.GERATIPO(K).EQ.14)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2005)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
!
        READ(IN_FILE,7000)FILE_VAL(K)
        IF(NK.EQ.0)WRITE(*,7002)FILE_VAL(K)
!
        READ(IN_FILE,7000)FILE_VET(K)
        IF(NK.EQ.0)WRITE(*,7003)FILE_VET(K)
!
        READ(IN_FILE,7000)FILE_THE(K)
        IF(NK.EQ.0)WRITE(*,7004)FILE_THE(K)
!
!     READ(IN_FILE,12)NPROPOSAL(K),SIG(K)
     !     IF(NK.EQ.0)WRITE(*,13)NPROPOSAL(K),SIG(K)
        READ(IN_FILE,12)NPROP,SIGG
        IF(NK.EQ.0)WRITE(*,13)NPROP,SIGG
     ENDIF
!
     IF(GERATIPO(K).EQ.33)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2005)MMKL
!
        READ(IN_FILE,1000)LBDFIXO(K)
        IF(NK.EQ.0)WRITE(*,105)LBDFIXO(K)
!
        READ(IN_FILE,1000)NEWLBD(K)
        IF(NK.EQ.0)WRITE(*,100)NEWLBD(K)
!
        READ(IN_FILE,3003)NZMIN(K),ZMAX(K)
        IF(NK.EQ.0)WRITE(*,106)NZMIN(K),ZMAX(K)
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
!
        READ(IN_FILE,7000)FILE_VET(K)
        IF(NK.EQ.0)WRITE(*,7003)FILE_VET(K)
!
        READ(IN_FILE,7000)FILE_THE(K)
        IF(NK.EQ.0)WRITE(*,7004)FILE_THE(K)
!
!     READ(IN_FILE,12)NPROPOSAL(K),SIG(K)
!     IF(NK.EQ.0)WRITE(*,13)NPROPOSAL(K),SIG(K)
        READ(IN_FILE,12)NPROP,SIGG
        IF(NK.EQ.0)WRITE(*,13)NPROP,SIGG
!
     ENDIF
!
     IF(GERATIPO(K).EQ.4)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2006)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
     ENDIF
!
     IF(GERATIPO(K).EQ.5)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2006)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
     ENDIF
!
     IF(GERATIPO(K).EQ.0)THEN
        READ(IN_FILE,3000)MMKL
        IF(NK.EQ.0)WRITE(*,2006)MMKL
!
        READ(IN_FILE,7000)FILE_OUT(K)
        IF(NK.EQ.0)WRITE(*,7001)FILE_OUT(K)
!
        IN_FILE2 = 600+NK
        NAME = '../gera_LABTRANGEORW3/in/entra_niveis.in'
        OPEN(UNIT=IN_FILE2,FILE=NAME,STATUS='UNKNOWN',&
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',NAME
           STOP
        END IF
        READ(IN_FILE2,4000)(NSE(I),I=1,NS)
!
        READ(IN_FILE2,2010)NLIXO,NTOTAL(K)
!     WRITE(*,*)NLIXO,NTOTAL(K)
!
        READ(IN_FILE2,1003)ALPHA(K)
        !
        READ(IN_FILE2,7000)NAME
!
        CLOSE(UNIT=IN_FILE2)
!
     ENDIF
!
     READ(IN_FILE,3000)NCOND(K)
     IF(NK.EQ.0)WRITE(*,3002)NCOND(K)
     IF(NCOND(K).GT.NCONDMAX)THEN
        IF(NK.EQ.0)WRITE(*,505)NCOND(K),NCONDMAX
        STOP
     END IF
!
     IF(GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.GERATIPO(K).EQ.13&
          .OR.GERATIPO(K).EQ.14)THEN
        DO I=0,NCOND(K)-1
           READ(IN_FILE,1011)VET(K,0,I),VET(K,1,I),VET(K,2,I),&
                VET(K,3,I)
           IF(NK.EQ.0)WRITE(*,1012)I+1,VET(K,0,I),VET(K,1,I),VET(K,2,I),&
                VET(K,3,I)
           DO J=0,2
              VET(K,J,I)=VET(K,J,I)+1E-6
           ENDDO
        ENDDO
     ELSE
        DO I=0,NCOND(K)-1
           READ(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
           IF(NK.EQ.0)WRITE(*,1002)I+1,VET(K,0,I),VET(K,1,I),VET(K,2,I)
           DO J=0,1
              VET(K,J,I)=VET(K,J,I)+1E-6
           ENDDO
        ENDDO
     END IF
!
     NERROREAD = 1
!
     IF(GERATIPO(K).EQ.33.OR.GERATIPO(K).EQ.32)THEN
        IF(NK.EQ.0)WRITE(*,*)'((((((((((((((.)))))))))))))))))'
        IF(NK.EQ.0)WRITE(*,*)'DETALHES PARA VARIACAO EM LAMBDA'
        FILENT='in/entra_campo.in'
        IN_FILE3=700+NK
        OPEN(UNIT=IN_FILE3,FILE=FILENT,STATUS='OLD',&
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILENT
           STOP
        END IF
        READ(IN_FILE3,*)XLBDMEAN(K),XLBDSD(K),XLBDBETA(K)
        IF(NK.EQ.0)WRITE(*,104)XLBDMEAN(K),XLBDSD(K),XLBDBETA(K)
        READ(IN_FILE3,*)NKIND(K)
        IF(NKIND(K).EQ.0)THEN
           IF(NK.EQ.0)WRITE(*,*)'((((((((((((((.)))))))))))))))))'
           IF(NK.EQ.0)WRITE(*,*)'RANDOM WALK'
           IF(NK.EQ.0)WRITE(*,*)'((((((((((((((.)))))))))))))))))'
        END IF
        IF(NKIND(K).EQ.1)THEN
           IF(NK.EQ.0)WRITE(*,*)'((((((((((((((.)))))))))))))))))'
           IF(NK.EQ.0)WRITE(*,*)'AUTOREGRESSIVO'
           IF(NK.EQ.0)WRITE(*,*)'((((((((((((((.)))))))))))))))))'
        END IF
        READ(IN_FILE3,*)I
        IF(I.NE.NLPRIOR(K))THEN
           IF(NK.EQ.0)WRITE(*,*)'PROBLEMAS EM NLPRIOR'
           STOP
        END IF
        DO I=1,NLPRIOR(K)
           READ(IN_FILE3,*)CLBDS(I)
           IF(NK.EQ.0)WRITE(*,*)CLBDS(I)
        END DO
        CLOSE(IN_FILE3)
!     
        IF(CONT.EQ.0)THEN
           IF(TIPOVAR.EQ.0)THEN
              YOLDLBD(K)=(NEWLBD(K)-XLBDMEAN(K))/XLBDSD(K)
              OLDLBD(K)=NEWLBD(K)
           ELSE
              YOLDLBD(K)=(NEWLBD(K)-XLBDMEAN(K))/(XLBDSD(K)-XLBDMEAN(K))
              OLDLBD(K)=NEWLBD(K)
           END IF
           IF(NK.EQ.0)WRITE(*,*)"PRIMEIRO LAMBDA (NORMAL):",YOLDLBD(K)
        ELSE
           YOLDLBD(K)=(OLDLBD(K)-XLBDMEAN(K))/XLBDSD(K)
           IF(NK.EQ.0)WRITE(*,*)"ULTIMO  LAMBDA (NORMAL):",YOLDLBD(K)
        END IF
     END IF
  END IF
!
  CLOSE(UNIT=IN_FILE)
!
  IF(NK.EQ.0)WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
  IF(NK.EQ.0)WRITE(*,*)FILE_OUT(K)
  IF(NK.EQ.0)WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MULTIPLOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LNGT = LEN(FILE_OUT(K))
  FOUT = TRIM(FILE_OUT(K))
  DO I=1,LNGT-2
     FEXP = FOUT(I:I+2)
     IF(FEXP.EQ.'exp')THEN
        GOTO 111
     END IF 
  END DO
111 CONTINUE
  DO J=I,LNGT-2
     FEXP = FOUT(J:J)
     IF(FEXP.EQ.'/')THEN
        GOTO 112
     END IF 
  END DO
112 CONTINUE
  FILE_OUT(K) = ADJUSTL(TRIM(FOUT(1:I-1)))
  FILE_EXT(K) = ADJUSTL(TRIM(FOUT(J:LNGT)))
!
!  FOUT = TRIM(FILE_THE(K))
!  DO I=1,LNGT-2
!     FEXP = FOUT(I:I+2)
!     IF(FEXP.EQ.'dat')THEN
!        GOTO 110
!     END IF 
!  END DO
  !110 CONTINUE
  FOUT = ('../out/theta')
  FILE_THE(K) = TRIM(ADJUSTL(FOUT(1:I-2)))
  NERROREAD = 1
!
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
106 FORMAT('MALHA NO TEMPO (PREVISAO)..:',I7,/,&
           'TEMPO PARA PREVISAO........:',F10.3)
105 FORMAT('LAMBDA FIXO................:',F10.3)
104 FORMAT('LAMBDA MEDIO...............:',F10.3,/,&
           'DESVIO DE LAMBDA...........:',F10.3,/,&
           'BETA PARA PASSO NA CADEIA..:',F10.3)
101 FORMAT('LAMBDA.....................:',F10.3)
102 FORMAT('NUMERO DE LAMBDAS..........:',I7)
103 FORMAT('LAMBDAS USADOS.............:',10F10.3)
100 FORMAT('NOVO LAMBDA................:',F10.3)
12 FORMAT(I7,E10.3)
13 FORMAT('TIPO DE PROPOSAL         =',I7,/,&
          'PARAMETRO DO RANDOM WALK =',E10.3)
1000 FORMAT(3F12.5)
1001 FORMAT('DIMENSION X =',F12.5,/,&
            'DIMENSION Y =',F12.5,/)
8801 FORMAT('DIMENSION X =',F12.5,/,&
            'DIMENSION Y =',F12.5,/,&
            'DIMENSION Z =',F12.5,/)
1011 FORMAT(4F12.5)
1012 FORMAT('POINT    =',I6,/,&
          'COORD. X =',F12.5,/,&
          'COORD. Y =',F12.5,/,&
          'COORD. Z =',F12.5,/,&
          'VALUE    =',F12.5,/)
1002 FORMAT('POINT    =',I6,/,&
          'COORD. X =',F12.5,/,&
          'COORD. Y =',F12.5,/,&
          'VALUE    =',F12.5,/)
1003 FORMAT(F12.5)
2000 FORMAT(3I7)
2001 FORMAT('N. OF ELEMENTS IN X =',I7,/,&
            'N. OF ELEMENTS IN Y =',I7,/)
9901 FORMAT('N. OF ELEMENTS IN X =',I7,/,&
            'N. OF ELEMENTS IN Y =',I7,/,&
            'N. OF ELEMENTS IN Z =',I7,/)
2005 FORMAT('M (terms) =',I7,/)
2006 FORMAT('M (NIVEIS)=',I7,/)
2010 FORMAT(2I7)
3000 FORMAT(I7)
3003 FORMAT(I7,F10.3)
3001 FORMAT('FLAG  =',I7)
3002 FORMAT('NCOND =',I7)
4000 FORMAT(40I12)
4001 FORMAT(/,'SEEDS =',40I12,/)
5001 FORMAT('DELTA X =',F12.4,/,&
          'DELTA Y =',F12.4,/)
6001 FORMAT(I7,': EXPONENTIAL')
6002 FORMAT(I7,': FRACTAL')
7000 FORMAT(A)
7001 FORMAT('NAME OF OUTPUT FILE = ',A)
7002 FORMAT('NAME OF INPUT AUTOVALUES FILE  = ',A)
7003 FORMAT('NAME OF INPUT AUTOVECTORS FILE = ',A)
7004 FORMAT('NAME OF INPUT THETA (NORMAL RV)= ',A)
8000 FORMAT(2I7)
8001 FORMAT('N. OF INITIAL FILE =',I7,/,&
            'N. OF FINAL FILE   =',I7)
8002 FORMAT('NUMBER OF FILES    =',I7,/)
9000 FORMAT(2F12.4)
9002 FORMAT('HURST COEF.=',F12.4,/,&
          'VARIANCE   =',F12.4)
9001 FORMAT('CORRELATION LENGTH...=',F12.4,/,&
            'VARIANCE.............=',F12.4)
505 FORMAT('PROBLEMA NA DEFINICAO DE NCONDMAX',/,&
          'NUMERO DE PONTOS CONDICIONADOS SUPERA NCONDMAX',/,&
          'FAVOR ALTERAR NCONDMAX',/,&
          'NCOND=',I5,'NCONDMAX=',I5)
!
END SUBROUTINE READ_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ REFERENCE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_REF(NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER :: I,J,K,IN_FILE,ISTAT,IERR,NK
  REAL   , EXTERNAL :: MAXIMO,NORMA_L2
!
  DO K=1,NDTYPE
!
     IN_FILE = 18+NK
!
     WRITE(*,*)'##########################################'
     WRITE(*,*)'####### LENDO DADOS DE REFERENCIA ########'
     WRITE(*,*)'ARQUIVO REF:', TRIM(ADJUSTL(FILEREF(K)))
!
     OPEN(UNIT=IN_FILE,FILE=FILEREF(K),STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE REF: ',FILEREF(K)
        STOP
     END IF
!
     DO I=1,NDADOS(K)
!        READ(IN_FILE,100)(REFDATA(K,J,I),J=1,NWELLS(K)+1)
        READ(IN_FILE,*)(REFDATA(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
!
     CLOSE(UNIT=IN_FILE)
!     MAXDATA(K)=MAXIMO(K)
     IF(NNORMA.EQ.1)THEN
        REFAMOS = 0.0d0
        NORMA_L2_REF(K)=NORMA_L2(K)
     ELSE
        NORMA_L2_REF(K)=1.D0
     END IF
!
  END DO
!
100  FORMAT(30(e15.8,2x))
!
END SUBROUTINE READ_REF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DEFINE O VALOR MAXIMO DOS DADOS !!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION MAXIMO(K)
  USE VARIAVEIS
  INTEGER, intent(in) :: K
  INTEGER :: I,J
  REAL    :: AUX
!
  AUX=-1E30
  DO I=1,NDADOS(K)
     DO J=2,NWELLS(K)+1
        AUX = MAX(AUX,REFDATA(K,J,I))
     END DO
  END DO
  MAXIMO = AUX
END FUNCTION MAXIMO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DEFINE O VALOR MAXIMO DOS DADOS !!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION LOGERRORMEAN(VET,K,NT,NFREQ)
!
  INTEGER, INTENT(IN)                :: K,NT,NFREQ
  REAL   , INTENT(IN), DIMENSION(NT) :: VET
  INTEGER                            :: I,NI,NF,NN,NL
  REAL                               :: AUX,DH
!
  NL = NFREQ
  NF = K
  NI = K/2
  NN = NF-NI+1
  IF(NN.GT.NL) NI = NF-NL
  NN = NF-NI+1
  DH = REAL(NN)
  AUX = 0.0
!
  DO I=NI,NF
     AUX = VET(I) + AUX
  END DO
!
  LOGERRORMEAN = AUX/DH
!  
END FUNCTION LOGERRORMEAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DEFINE O VALOR MINIMO DOS DADOS !!!!!!!!!!!!!!!!!!!!!
INTEGER FUNCTION INTMIN(VET,K)
  INTEGER, intent(in) :: K
  INTEGER, INTENT(IN),DIMENSION(K) :: VET
  INTEGER :: I,J
!
  INTMIN=VET(1)
  DO I=2,K
     IF(VET(I)<INTMIN) INTMIN=VET(I)
  END DO
END FUNCTION INTMIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DEFINE O VALOR MAXIMO DOS DADOS !!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION NORMA_L2(K)
  USE VARIAVEIS
  INTEGER, intent(in) :: K
  INTEGER :: I,J
  REAL             :: INTEGRAL,AUX,AUXR
!
     INTEGRAL=0.D0
     AUX     =0.D0
     IF(LIKETIPO(K).EQ.1)THEN
        WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        WRITE(*,*)'@@ NAO ESTA BEM DEFINIDO @@'
        WRITE(*,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        STOP
     END IF
!
!     IF(LIKETIPO(K).EQ.1)THEN
!        DO I=1,NWELLS(K)
!           DO J=1,NDADOS(K)-1
!              DTR =(REFDATA(K,1,J+1)-REFDATA(K,1,J))*0.5D0
!              DTA =(REFAMOS(K,1,J+1)-REFAMOS(K,1,J))*0.5D0
!              AUXR=(REFDATA(K,I+1,J+1)-REFAMOS(K,I+1,J+1))
!              AUXA=(REFDATA(K,I+1,J)-REFAMOS(K,I+1,J))
!              AUX = AUX+DTR*(AUXR*AUXR)+DTA*(AUXA*AUXA)
!           ENDDO
!        ENDDO
!     ENDIF
!
     IF(LIKETIPO(K).EQ.2)THEN
        DO I=2,NWELLS(K)+1
           DO J=NDADOSI(K),NDADOS(K)
              AUXR=REFAMOS(K,I,J)-REFDATA(K,I,J)
              AUX =AUX+AUXR*AUXR
           ENDDO
        ENDDO
        NORMA_L2=(AUX)
     ENDIF
!
     IF(LIKETIPO(K).EQ.3)THEN
        AUX = -1.0D+20
        DO I=2,NWELLS(K)+1
           DO J=NDADOSI(K),NDADOS(K)
              AUXR=ABS(REFAMOS(K,I,J)-REFDATA(K,I,J))
              AUX =MAX(AUX,AUXR)
           ENDDO
        ENDDO
        NORMA_L2=AUX
     ENDIF
!
     write(*,100)k,NORMA_L2
!
100  FORMAT('DATA no....:',I3,2X,'NORMA L2........:',E10.3)
!
   END FUNCTION NORMA_L2
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ REFERENCE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_AMOSTRA()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER :: I,J,K,IN_FILE,ISTAT
!
  DO K=1,NDTYPE
!
     IN_FILE = 22
!
     OPEN(UNIT=IN_FILE,FILE=FILEAM(K),STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE SAMPLE: ',FILEAM(K)
        STOP
     END IF
!
     DO I=1,NDADOS(K)
        READ(IN_FILE,100)(REFAMOS(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
!
     CLOSE(UNIT=IN_FILE)
!
  END DO
!
100  FORMAT(30(e15.8,2x))
!
END SUBROUTINE READ_AMOSTRA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WRITE THE NEW INPUT FILE TO FIELDS !!!!!!!!!!!!!!!!!!!!
SUBROUTINE GERAFILEIN(NERROREAD,NRK,NPROC,K)
!
  USE VARIAVEIS
  USE RANDOM
  IMPLICIT NONE
!
  CHARACTER(LEN=128)    :: FILEOUT,FILEVET,FILETHETA
  CHARACTER(LEN=4)      :: NUMB
  INTEGER               :: IN_FILE,ISTAT,I,J,NN,K,INIF,NFILES
  INTEGER, INTENT(IN)   :: NRK,NPROC
  INTEGER               :: NPROP
  INTEGER, INTENT(OUT)  :: NERROREAD
  REAL                  :: TESTE,AUX,BETA,VARUNI
  REAL   ,DIMENSION(NS) :: SEEDS
  INTEGER,DIMENSION(NS) :: N_SEED
  REAL                  :: SIGK
  REAL   ,EXTERNAL      :: RANDOM_WALK_METHOD
!
  INIF   = 0
  NFILES = 0
!
!  DO K=1,NPRIORR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NPROP = NPROPOSAL(K)
  IF(NPROPOSAL(K).EQ.2)THEN
     SIGK = RANDOM_WALK_METHOD(SIG(K),MKL(K),CONTADORC,&
          NFREQ(K),VARPROP(K))
     JUMP = SIGK
     NPROP= 2
  END IF
  IF(NPROPOSAL(K).EQ.3)THEN
     CALL AM_METHOD(K,NRK,SIGK)
     NPROP= 2
     SIGK = SIG(K)
  END IF
  IF(NPROPOSAL(K).EQ.4)THEN
     CALL DE_METHOD(K,NRK,NPROC,MKL(K))
     NPROP = 4
     SIGK  = 0.0
  END IF
  IF(NPROPOSAL(K).EQ.5)THEN
     CALL DREAM_METHOD(K,NRK,NPROC,MKL(K),PCRALL(K,:))
     NPROP = 4
     SIGK  = 0.0
  END IF
  IF(NPROPOSAL(K).EQ.6)THEN
     CALL DREAM_METHOD_MOD(K,NRK,NPROC,MKL(K),PCRALL(K,:))
     NPROP = 4
     SIGK  = 0.0
  END IF
!! SEMENTES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_NUMBER(SEEDS)
  DO I=1,NS
     N_SEED(I) = FLOOR(SEEDS(I)*1.0E+08)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NLPRIOR(K).NE.0)THEN
     IF(GERATIPO(K).EQ.32)THEN
        IF(CONT.EQ.0)THEN
           OLDLBD(K)=XLBDMEAN(K)
           DO I=1,NLPRIOR(K)
              IF(ABS(CLBDS(I)-OLDLBD(K)).LE.1E-06)THEN
                 NUM_AVETOLD(K)=I
              END IF
           END DO
        END IF
        CALL RANDOM_NUMBER(TESTE)
        WRITE(*,*)'************************************'
        WRITE(*,*)'************************************'
        WRITE(*,103)XLBDSD(K)
        WRITE(*,101)OLDLBD(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(NKIND(K).EQ.0)THEN
           IF(TESTE.GT.1.0D0-(1.0D0-XLBDSD(K))*0.5D0)THEN
              NUM_AVET(K)=NUM_AVETOLD(K)+1
              IF(NUM_AVET(K).GT.NLPRIOR(K))NUM_AVET(K)=NLPRIOR(K)
           ELSE IF(TESTE.LT.(1.0D0-XLBDSD(K))*0.5D0)THEN
              NUM_AVET(K)=NUM_AVETOLD(K)-1
              IF(NUM_AVET(K).LT.1)NUM_AVET(K)=1
           ELSE
              NUM_AVET(K)=NUM_AVETOLD(K)
           END IF
        END IF
!
        IF(NKIND(K).EQ.1)THEN
           IF(TESTE.LT.0.05)NUM_AVET(K)=NUM_AVETOLD(K)-2
           IF(TESTE.GT.0.05.AND.TESTE.LT.0.25)NUM_AVET(K)=NUM_AVETOLD(K)-1
           IF(TESTE.GT.0.25.AND.TESTE.LT.0.75)NUM_AVET(K)=NUM_AVETOLD(K)
           IF(TESTE.GT.0.75.AND.TESTE.LT.0.95)NUM_AVET(K)=NUM_AVETOLD(K)+1
           IF(TESTE.GT.0.95)NUM_AVET(K)=NUM_AVETOLD(K)+2
           IF(NUM_AVET(K).LT.1)NUM_AVET(K)=1
           IF(NUM_AVET(K).GT.NLPRIOR(K))NUM_AVET(K)=NLPRIOR(K)
        END IF
!
        NEWLBD(K)=CLBDS(NUM_AVET(K))
        WRITE(*,100)NEWLBD(K)
        WRITE(*,*)'************************************'
        FILE_VET(K)=ADJUSTL(TRIM(NAME_AVET(K,NUM_AVET(K))))
     END IF
  END IF
!
  IF(GERATIPO(K).EQ.33)THEN
     WRITE(*,*)'************************************'
     WRITE(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(*,*)'************************************'
     write(*,*)'BETA (PASSO NA CADEIA)..:',XLBDBETA(K)
     write(*,*)'LAMBDA MEDIO............:',XLBDMEAN(K)
     write(*,*)'DESVIO DE LAMBDA........:',XLBDSD(K)
     WRITE(*,*)'************************************'
     WRITE(*,*)'LAMBDA ATUAL............:',OLDLBD(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(NKIND(K).EQ.0)THEN
        YNEWLBD(K) = YOLDLBD(K)+XLBDBETA(K)*RANDOM_NORMAL()
     END IF
     IF(NKIND(K).EQ.1)THEN
        YNEWLBD(K)=SQRT(1.0D0-XLBDBETA(K)*XLBDBETA(K))*YOLDLBD(K)+ &
             XLBDBETA(K)*RANDOM_NORMAL()
     END IF
     IF(TIPOVAR.EQ.0)THEN
        NEWLBD(K)  = XLBDMEAN(K)+XLBDSD(K)*YNEWLBD(K)
     ELSE
        NEWLBD(K)  = XLBDMEAN(K)+(XLBDSD(K)-XLBDMEAN(K))*YNEWLBD(K)
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     NEWLBD(K)=NINT(NEWLBD(K))
     WRITE(*,*)'LAMBDA TESTADO..........:',NEWLBD(K)
     WRITE(*,*)'************************************'
     WRITE(*,*)'////////////////////////////////////'
     WRITE(*,*)'************************************'
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  WRITE(NUMB,'(I4.3)')NRK
  FILEOUT = TRIM(ADJUSTL(FILE_OUT(K)))//('exp')//TRIM(ADJUSTL(NUMB))//&
       TRIM(ADJUSTL(FILE_EXT(K)))
!
  IN_FILE = 800+NRK
  OPEN(UNIT=IN_FILE,FILE=FILE_INFIELD(K),STATUS='REPLACE',&
       ACTION='WRITE',FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE FIELD: ',FILE_INFIELD(K)
     STOP
  END IF
!
  IF(GERATIPO(K).EQ.6.OR.GERATIPO(K).EQ.7.OR.&
       GERATIPO(K).EQ.8)THEN
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VET(K)))
!     WRITE(NUMB,'(I4.3)')NRK
     FILETHETA = ADJUSTL(TRIM(FILE_THE(K)))//TRIM(ADJUSTL(NUMB))//&
          ('.dat')
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILETHETA))
     WRITE(IN_FILE,1200)NPROP,SIGK
  END IF
!
!
  IF(GERATIPO(K).EQ.1)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VAL(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VET(K)))
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.2)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  FILETHETA = TRIM(ADJUSTL(FILE_THE(K)))//TRIM(ADJUSTL(NUMB))//&
       ('.dat')
!
  IF(GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.GERATIPO(K).EQ.34)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VAL(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VET(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILETHETA))
     WRITE(IN_FILE,1200)NPROP,SIGK
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.13.OR.GERATIPO(K).EQ.14)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K),DIMZ(K)
     WRITE(IN_FILE,2000)NX(K),NY(K),NZ(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VAL(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VET(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILETHETA))
     WRITE(IN_FILE,1200)NPROPOSAL(K),SIGK
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1001)VET(K,0,I),&
             VET(K,1,I),VET(K,2,I),VET(K,3,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.33)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K),DIMZ(K)
     WRITE(IN_FILE,2000)NX(K),NY(K),NZ(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,1000)LBDFIXO(K)
     WRITE(IN_FILE,1000)NEWLBD(K)
     WRITE(IN_FILE,3003)NZMIN(K),ZMAX(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILE_VET(K)))
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILETHETA))
     WRITE(IN_FILE,1200)NPROPOSAL(K),SIGK
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1001)VET(K,0,I),&
             VET(K,1,I),VET(K,2,I),VET(K,3,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.4)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.5)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.0)THEN
     WRITE(IN_FILE,8000)INIF,NFILES
     WRITE(IN_FILE,1000)DIMX(K),DIMY(K)
     WRITE(IN_FILE,2000)NX(K),NY(K)
     WRITE(IN_FILE,3000)NFLAG(K)
     WRITE(IN_FILE,9000)FBETA(K),VARIAN(K)
     WRITE(IN_FILE,4000)(N_SEED(I),I=1,NS)
     WRITE(IN_FILE,3000)MKL(K)
     WRITE(IN_FILE,7000)ADJUSTL(TRIM(FILEOUT))
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  NERROREAD = 1d0
  CLOSE(UNIT=IN_FILE)
!  END DO 
!
  RETURN
!
1002 FORMAT(1000F12.5)
100  FORMAT('LAMBDA TESTADO............:',F10.3)
101  FORMAT('LAMBDA ATUAL..............:',F10.3)
103  FORMAT('PROB. DE PERMANECER EM n..:',F10.3)
1200 FORMAT(I7,E10.3)
1000 FORMAT(3F12.5)
1001 FORMAT(4F12.5)
2000 FORMAT(3I7)
3000 FORMAT(I7)
3003 FORMAT(I7,F10.3)
4000 FORMAT(40I12)
7000 FORMAT(A)
8000 FORMAT(2I7)
9000 FORMAT(2F12.4)
!
END SUBROUTINE GERAFILEIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GERADOR DE CAMPO ALEATORIO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GERADOR(NN,NK,K)
!
  USE VARIAVEIS, ONLY: THETAN,THETA,MKL,GERATIPO,NTOTAL
  USE VARIAVEIS, ONLY: NSIMUL,LOGPERM,NPRIORR,FILE_OUT
!
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: COMMAND
  INTEGER            :: OUT_FILE,ISTAT,I,FLAG,MK
  INTEGER, INTENT(IN):: NN,K,NK
  CHARACTER(LEN=128) :: NOME,FOUT
  CHARACTER(LEN=4)   :: NUMB
!
  WRITE(NUMB,'(I4.3)')NK
!  
  IF(GERATIPO(K).EQ.0)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(NTOTAL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_LABTRANGEORW3/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.1)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(GERATIPO(K).EQ.2)THEN
     COMMAND=('cd ../gera_LABTRANGEO/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.3)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.31)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW1/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.34)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW2/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.32)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW3D/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.33)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_LBD/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.4)THEN
     COMMAND=('cd ../gera_LABTRANGEORW/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.5)THEN
     COMMAND=('cd ../gera_LABTRANGEORW2/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.6)THEN
     COMMAND=('cd ../gera_KL/MVN/; sh rodarMVN.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.7)THEN
     COMMAND=('cd ../gera_KL/MVN1/; sh rodarMVN.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.8)THEN
     COMMAND=('cd ../gera_KL/MVN2/; sh rodarMVN.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.13)THEN
     COMMAND=('cd ../gera_KL/FORTRAN_KL3D/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.14)THEN
     COMMAND=('cd ../gera_KL/FORTRAN_KL3D_2/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output_')//TRIM(ADJUSTL(NUMB))//('.out')
  ENDIF
!
  WRITE(*,*)COMMAND
  CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(GERATIPO(K).EQ.0.OR.GERATIPO(K).EQ.1.OR.&
       GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.&
       GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.&
       GERATIPO(K).EQ.34.OR.GERATIPO(K).EQ.13.OR.&
       GERATIPO(K).EQ.14)THEN
!!!!!!!!!!!!!!!! LEITURA DO THETA ANTERIOR !!!!!!!!!!!!!!!!!!!!!!
     OUT_FILE = 900+NK
     IF(GERATIPO(K).EQ.0)THEN
        NOME=ADJUSTL(TRIM('../gera_LABTRANGEORW3/out/xpsi.dat'))
        MK=NTOTAL(K)
     END IF
     IF(GERATIPO(K).EQ.1)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.3)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.31)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.34)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.32)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.33)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.6)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.7)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN1/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.8)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.13)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.14)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D_2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
!
     NOME = TRIM(ADJUSTL(NOME))
     OPEN(UNIT=OUT_FILE,FILE=NOME,FORM='FORMATTED',STATUS='OLD',&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE THETA: ',NOME
        STOP
     END IF
!
     DO I=1,MK
        READ(OUT_FILE,111)THETAN(K,I)
     ENDDO
           !
     CLOSE(UNIT=OUT_FILE)
!          
!!!!!!!!!!!!!!!! LEITURA DO THETA ATUAL !!!!!!!!!!!!!!!!!!!!!!!!!
     OUT_FILE = 1000+NK
     IF(GERATIPO(K).EQ.0)THEN
        NOME=ADJUSTL(TRIM('../gera_LABTRANGEORW3/out/xpsi.dat'))
        MK=NTOTAL(K)
     END IF
     IF(GERATIPO(K).EQ.1)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.3)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.31)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.34)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.32)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.33)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.6)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.7)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN1/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.8)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN2/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.13)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
     IF(GERATIPO(K).EQ.14)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D_2/out/thetanew'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        MK=MKL(K)
     END IF
!
     NOME = TRIM(ADJUSTL(NOME))
     OPEN(UNIT=OUT_FILE,FILE=NOME,FORM='FORMATTED',STATUS='OLD',&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE (THETANEW): ',NOME
        STOP
     END IF
!
     DO I=1,MK
        READ(OUT_FILE,111)THETA(K,I)
     ENDDO
!
     CLOSE(UNIT=OUT_FILE)
!           
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IF(NSIMUL.EQ.1)THEN
!        NOME='../simulador/fields/amostra_0.dat'
!     END IF
!     IF(NSIMUL.EQ.2)THEN
!        NOME='../simuladorBTMM/exp/fields/amostra_0.dat'
!     END IF
!     IF(NSIMUL.EQ.3)THEN
!        NOME='../simul_comp/exp/fields/amostra_0.dat'
!     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(GERATIPO(K).EQ.2)THEN
     CALL READPERM(NPRIORR,NOME,LOGPERM,K)
  ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(GERATIPO(K).EQ.4)THEN
     CALL READPERM(NPRIORR,NOME,LOGPERM,K)
  ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(GERATIPO(K).EQ.5)THEN
     CALL READPERM(NPRIORR,NOME,LOGPERM,K)
  ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE GERADOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GERATHETA(MK,TIPO,K,NK)
!
  USE RANDOM
  USE VARIAVEIS, ONLY: THETAN
  IMPLICIT NONE
!
  INTEGER            :: I,N,MK,TIPO,K,NP,NK
  INTEGER            :: ISTAT,OUT_FILE
  CHARACTER(LEN=128) :: FILETHE,FOUT
  CHARACTER(LEN=4)   :: NUMB
  REAL               :: TESTE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  OUT_FILE = 1100+NK
  IF(TIPO.EQ.0)  FILETHE='../gera_LABTRANGEORW3/out/xpsi'
  IF(TIPO.EQ.1)  FILETHE='../gera_KL/FORTRAN/out/theta'
  IF(TIPO.EQ.3)  FILETHE='../gera_KL/FORTRAN_RW/out/theta'
  IF(TIPO.EQ.13) FILETHE='../gera_KL/FORTRAN_KL3D/out/theta'
  IF(TIPO.EQ.14) FILETHE='../gera_KL/FORTRAN_KL3D_2/out/theta'
  IF(TIPO.EQ.31) FILETHE='../gera_KL/FORTRAN_RW1/out/theta'
  IF(TIPO.EQ.34) FILETHE='../gera_KL/FORTRAN_RW2/out/theta'
  IF(TIPO.EQ.32) FILETHE='../gera_KL/FORTRAN_RW3D/out/theta'
  IF(TIPO.EQ.33) FILETHE='../gera_KL/FORTRAN_LBD/out/theta'
  IF(TIPO.EQ.6)  FILETHE='../gera_KL/MVN/out/theta'
  IF(TIPO.EQ.7)  FILETHE='../gera_KL/MVN1/out/theta'
  IF(TIPO.EQ.8)  FILETHE='../gera_KL/MVN2/out/theta'
!
  WRITE(NUMB,'(I4.3)')NK
  FOUT = TRIM(ADJUSTL(FILETHE))//TRIM(ADJUSTL(NUMB))//&
       ('.dat')
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE THETA: ',FILETHE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,MK
     THETAN(K,I) = RANDOM_NORMAL()
     WRITE(OUT_FILE,111)THETAN(K,I)
  ENDDO
!
  CLOSE(UNIT=OUT_FILE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  OUT_FILE = 1200+NK
  IF(TIPO.EQ.0) FILETHE='../gera_LABTRANGEORW3/out/xpsi'
  IF(TIPO.EQ.1) FILETHE='../gera_KL/FORTRAN/out/thetanew'
  IF(TIPO.EQ.3) FILETHE='../gera_KL/FORTRAN_RW/out/thetanew'
  IF(TIPO.EQ.13)FILETHE='../gera_KL/FORTRAN_KL3D/out/thetanew'
  IF(TIPO.EQ.14)FILETHE='../gera_KL/FORTRAN_KL3D_2/out/thetanew'
  IF(TIPO.EQ.31)FILETHE='../gera_KL/FORTRAN_RW1/out/thetanew'
  IF(TIPO.EQ.34)FILETHE='../gera_KL/FORTRAN_RW2/out/thetanew'
  IF(TIPO.EQ.32)FILETHE='../gera_KL/FORTRAN_RW3D/out/thetanew'
  IF(TIPO.EQ.33)FILETHE='../gera_KL/FORTRAN_LBD/out/thetanew'
  IF(TIPO.EQ.6) FILETHE='../gera_KL/MVN/out/thetanew'
  IF(TIPO.EQ.7) FILETHE='../gera_KL/MVN1/out/thetanew'
  IF(TIPO.EQ.8) FILETHE='../gera_KL/MVN2/out/thetanew'
!
  WRITE(NUMB,'(I4.3)')NK
  FOUT = TRIM(ADJUSTL(FILETHE))//TRIM(ADJUSTL(NUMB))//&
       ('.dat')
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE THETANEW: ',FILETHE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,MK
     THETAN(K,I) = RANDOM_NORMAL()
     WRITE(OUT_FILE,111)THETAN(K,I)
  ENDDO
!
  CLOSE(UNIT=OUT_FILE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE GERATHETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             
SUBROUTINE READPERM(nnp,fname,PERM,NK)
!
!     Objetivo: le as permeabilidades de um arquivo
!
!--------------------------------------------------------
! 
  implicit none
!      
  integer :: nelemx,nelemy,ntype,inperm,nline,nflag,i,j
  integer :: nr,nelx,nely,NK,nnp
  real   :: xlx,xly,beta
  character(len=128) :: NOME,FILEINN
  character(len=128) :: fname
  character(len=4)   :: EXT,tipo
  character(len=5)   :: C
  REAL   ,DIMENSION(nnp,*) :: PERM
!
!     abre o arquivo de permeabilidades
!   
  inperm=901
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  FILEINN=fname
  NOME=TRIM(FILEINN)
  NOME=ADJUSTL(TRIM(NOME))
  WRITE(*,111)NOME(1:LEN_TRIM(NOME))
  open(inperm, file= NOME)
!
!     dimensoes do dominio
!
  read(inperm,*) xlx
  read(inperm,*) xly
!
!     numero de elementos em cada direcao
!      
  read(inperm,*) nelemx
  read(inperm,*) nelemy
!
!     ntype = 1; campo exponencial
!     ntype = 2; campo fractal
!      
  read(inperm,*) ntype
!
!     beta: coeficiente de Hurst
!      
  read(inperm,*) beta
!
!     leituras vazias
!      
  read(inperm,*) 
  read(inperm,*) 
!
!     inicio da leitura do campo
!      
  do j=1,nelemy
!
     read(inperm,*) nline
     if(nline+1.ne.j) then
        write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
        stop
     end if
!
     read(inperm,*) (PERM(NK,i+(j-1)*nelemx),i=1,nelemx)
!      
     read(inperm,*) nflag
     if(nflag.ne.192837465) then
        write(*,*) 'Erro na leitura do campo de permeabilidade, nflag'
        stop
     end if
!      
  end do ! nelemy     
!     
  close(inperm)
!     calculando a permeabilidade lognormal
!
!     do j=1,nelemy*nelemx
!        PERM(NK,j)=kg*dexp(rho*PERM(NK,j))
!     end do
!
  tipo='K'
  call estatistica(nnp,PERM,nelemx,nelemy,tipo,NK)
!
111 FORMAT('NAME OF INPUT FILE (RANDOM FIELD): ',A)
113 FORMAT(I5)
115 FORMAT(                                      &
         '######################################',/, &
         'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
         'LEITURA DO CAMPO DE PERMEABILIDADES   ',/, &
         'NUMERO DO CAMPO:',I5,/,                    &
         '######################################',/)
END SUBROUTINE READPERM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ESTATISTICA(nnp,prm,nmx,nmy,sc,NK)
!
!     calculando a media e variancia
  implicit none
  REAL   , dimension(nnp,*) :: prm
  integer :: nmx,nmy,i,j,NK,nnp
  REAL             :: xlx,xly,xm,xmm,xv,XMAX,XMIN
  character*4 sc
!
  xm  = 0.d0
  xmm = 0.d0
  xv  = 0.d0
  XMAX = -1.d14
  XMIN = 1.d14
!
  do i=1,nmy*nmx
     xm = xm+prm(NK,i)
     xmm= xmm+prm(NK,i)*prm(NK,i)
     if(prm(NK,i).gt.XMAX) XMAX = prm(NK,i)
     if(prm(NK,i).lt.XMIN) XMIN = prm(NK,i)
  end do
  xm  = xm/(nmx*nmy)
  xmm = xmm/(nmx*nmy)
  xv  = xmm - xm*xm
!
!
  write(*,123)sc,xm,xv,XMAX,XMIN
123 FORMAT(A,                          &
         & ' ######################### ',/,  &
         & 'MEDIA     =',f10.5,/,            &
         & 'VARIANCIA =',f10.5,/,            &
         & 'MAXIMO    =',f10.5,/,            &
         & 'MINIMO    =',f10.5,/,            &
         & '##############################')
!
END SUBROUTINE ESTATISTICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PRIOR RATIO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PRIOR_RATIO()
!
  USE VARIAVEIS, only: THETA,THETAN,MKL,GERATIPO,PRATIO
  USE VARIAVEIS, only: NX,NY,LOGPERM,LOGPERMN,NTOTAL
  USE VARIAVEIS, only: NPRIORR,TPRATIO,NPROPOSAL
  USE VARIAVEIS, ONLY: YNEW,YOLD
  IMPLICIT NONE
!
  INTEGER          :: I,J,K
  REAL             :: INTEG,INTEGN
  REAL             :: AUX,AUXN
  REAL             :: MEAN,SIGMA
!
  IF(NPROPOSAL(1).EQ.2.OR.NPROPOSAL(1).EQ.3.OR.&
     NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5.OR.&
     NPROPOSAL(1).EQ.6)THEN
     DO K=1,NPRIORR
        PRATIO(K)=0.D0
     END DO
  ELSE     
     DO K=1,NPRIORR
!!! INICIALIZANDO AS VARIAVEIS !!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEG =0.D0
        INTEGN=0.D0
        AUX   =0.D0
        AUXN  =0.D0
        MEAN  =0.D0
        SIGMA =2.D0*1.D0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        PRATIO(K)=1.D0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.0)THEN
           DO I=1,NTOTAL(K)
              AUX   = THETA(K,I)-MEAN
              AUXN  = THETAN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.5)THEN
           DO I=1,NX(K)*NY(K)
              AUX   = LOGPERM(K,I)-MEAN
              AUXN  = LOGPERMN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.4)THEN
           DO I=1,NX(K)*NY(K)
              AUX   = LOGPERM(K,I)-MEAN
              AUXN  = LOGPERMN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.&
             GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.34.OR.&
             GERATIPO(K).EQ.13.OR.GERATIPO(K).EQ.14)THEN
           DO I=1,MKL(K)
              AUX   = THETA(K,I)-MEAN
              AUXN  = THETAN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.33)THEN
           DO I=1,MKL(K)
              AUX   = THETA(K,I)-MEAN
              AUXN  = THETAN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
!           AUXN= (YOLD**2)*0.5
!           AUX = (YNEW**2)*0.5
!           write(*,*)'################## prior:', exp(auxn-aux)
!
           PRATIO(K) = (INTEGN-INTEG)!+(AUXN-AUX)
!
        ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.2)THEN
           DO I=1,NX(K)*NY(K)
              AUX   = LOGPERM(K,I)-MEAN
              AUXN  = LOGPERMN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.1)THEN
           DO I=1,MKL(K)
              AUX = THETA(K,I)-MEAN
              AUXN= THETAN(K,I)-MEAN
              INTEG = INTEG + (AUX*AUX)
              INTEGN= INTEGN+ (AUXN*AUXN)
           ENDDO
!
           INTEG = INTEG/SIGMA
           INTEGN= INTEGN/SIGMA
!
           PRATIO(K) = (INTEGN-INTEG)
!
        ENDIF
!
     END DO
  END IF
!
  TPRATIO=0.0
  DO K=1,NPRIORR
     write(*,*)'ATENCAO',k,PRATIO(K)
     TPRATIO = TPRATIO+PRATIO(K)
  END DO
  TPRATIO = EXP(TPRATIO)
  WRITE(*,*)'###############################'
  WRITE(*,100)TPRATIO
  WRITE(*,*)'###############################'
!
  RETURN
!
100 FORMAT('PRIOR ERROR    = ',E10.3)
!
END SUBROUTINE PRIOR_RATIO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPSCALING DOS CAMPOS ALEATORIOS !!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UPSCALING(K,NK)
!
  USE VARIAVEIS, ONLY: NPRIORR,NUPSC
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: COMMAND
  CHARACTER(LEN=4)   :: NUMB
  CHARACTER(LEN=3)   :: FEXP
  INTEGER            :: K,NK,LNGT
  INTEGER            :: I,J,IN_FILE,ISTAT
  INTEGER            :: NA,NB,NC,ND,NE
  CHARACTER(LEN=128) :: FILEAMC,FIN,FOUT,FAUX
  REAL               :: XA,XB
  IN_FILE = 1300+NK
!
  WRITE(NUMB,'(I4.3)')NK
!
  IF(NUPSC(K).EQ.0)THEN
     FILEAMC = ('./upscaling/in')//TRIM(ADJUSTL(NUMB))//('/entrada.in')
  END IF
  IF(NUPSC(K).EQ.1)THEN
     FILEAMC = ('./upscaling1/in')//TRIM(ADJUSTL(NUMB))//('/entrada.in')
  END IF
  IF(NUPSC(K).EQ.2)THEN
     FILEAMC = ('./upscaling2/in')//TRIM(ADJUSTL(NUMB))//('/entrada.in')
  END IF
!
     WRITE(*,*)'#### RANK: ',NK
     WRITE(*,*)'#### LEITURA DO ARQUIVO DE ENTRADA UPSCALING ####'
     WRITE(*,*)FILEAMC
     OPEN(UNIT=IN_FILE,FILE=FILEAMC,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAMC
        STOP
     END IF
!
     READ(IN_FILE,100)NA
     READ(IN_FILE,101)NB,NC
     READ(IN_FILE,102)FIN
     READ(IN_FILE,*)XA,XB
     READ(IN_FILE,102)FOUT
     READ(IN_FILE,101)ND,NE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     LNGT = LEN(FIN)
     FAUX = TRIM(FIN)
     DO I=1,LNGT-2
        FEXP = FAUX(I:I+2)
        IF(FEXP.EQ.'exp')THEN
           GOTO 111
        END IF
     END DO
111  CONTINUE
     DO J=I,LNGT-2
        FEXP = FAUX(J:J)
        IF(FEXP.EQ.'/')THEN
           GOTO 115
        END IF
     END DO
115  CONTINUE
     FIN = TRIM(ADJUSTL(FIN(1:I-1)))//('exp')//TRIM(ADJUSTL(NUMB))//&
          TRIM(ADJUSTL(FIN(J:LNGT)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     LNGT = LEN(FOUT)
     FAUX = TRIM(FOUT)
     DO I=1,LNGT-2
        FEXP = FAUX(I:I+2)
        IF(FEXP.EQ.'exp')THEN
           GOTO 112
        END IF
     END DO
112  CONTINUE
     DO J=I,LNGT-2
        FEXP = FAUX(J:J)
        IF(FEXP.EQ.'/')THEN
           GOTO 116
        END IF
     END DO
116  CONTINUE
     FOUT = TRIM(ADJUSTL(FOUT(1:I-1)))//('exp')//TRIM(ADJUSTL(NUMB))//&
          TRIM(ADJUSTL(FOUT(J:LNGT)))
!
     CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE(*,*)'#### LEITURA DO ARQUIVO DE ENTRADA UPSCALING ####'
     WRITE(*,*)FILEAMC
     OPEN(UNIT=IN_FILE,FILE=FILEAMC,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAMC
        STOP
     END IF
!
     WRITE(IN_FILE,100)NA
     WRITE(IN_FILE,101)NB,NC
     WRITE(IN_FILE,102)FIN
     WRITE(IN_FILE,*)XA,XB
     WRITE(IN_FILE,102)FOUT
     WRITE(IN_FILE,101)ND,NE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     IF(NUPSC(K).EQ.0)THEN
        COMMAND = ('cd ./upscaling/; sh rodarUP.sh in')//TRIM(ADJUSTL(NUMB))
     END IF
     IF(NUPSC(K).EQ.1)THEN
        COMMAND = ('cd ./upscaling1/; sh rodarUP.sh in')//TRIM(ADJUSTL(NUMB))
     END IF
     IF(NUPSC(K).EQ.2)THEN
        COMMAND = ('cd ./upscaling2/; sh rodarUP.sh in')//TRIM(ADJUSTL(NUMB))
     END IF
     WRITE(*,*)'RANK:',NK,TRIM(COMMAND)
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT(I10)
101 FORMAT(2I10)
102 FORMAT(A)
200 FORMAT(30(e15.8,2x))
!
END SUBROUTINE UPSCALING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ UPSCALING INPUT DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_UPSC(K)
!
  USE VARIAVEIS, ONLY: NUPSC
  IMPLICIT NONE
!
  INTEGER            :: I,J,K,IN_FILE,ISTAT
  INTEGER            :: NA,NB,NC,ND,NE
  CHARACTER(LEN=128) :: FILEAMC,FIN,FOUT
  REAL               :: XA,XB
  IN_FILE = 1400
!
  IF(NUPSC(K).EQ.0) FILEAMC = './upscaling/in/entrada.in'
  IF(NUPSC(K).EQ.1) FILEAMC = './upscaling1/in/entrada.in'
  IF(NUPSC(K).EQ.2) FILEAMC = './upscaling2/in/entrada.in'
     
  WRITE(*,*)'#### LEITURA DO ARQUIVO DE ENTRADA UPSCALING ####'
  WRITE(*,*)FILEAMC
  OPEN(UNIT=IN_FILE,FILE=FILEAMC,STATUS='UNKNOWN',&
       FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAMC
     STOP
  END IF
     !
  READ(IN_FILE,100)NA
  READ(IN_FILE,101)NB,NC
  READ(IN_FILE,102)FIN
  READ(IN_FILE,*)XA,XB
  READ(IN_FILE,102)FOUT
  READ(IN_FILE,101)ND,NE
  WRITE(*,100)NA
  WRITE(*,101)NB,NC
  WRITE(*,102)FIN
  WRITE(*,*)XA,XB
  WRITE(*,102)FOUT
  WRITE(*,101)ND,NE
  
!
  CLOSE(UNIT=IN_FILE)
  !
100 FORMAT(I10)
101 FORMAT(2I10)
102 FORMAT(A)
200 FORMAT(30(e15.8,2x))
!
END SUBROUTINE READ_UPSC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SIMULADOR USADO NA MALHA GROSSA !!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SIMULADOR_C(NK)
!
  USE VARIAVEIS, ONLY: NSIMUL,ERROSIMULADOR
  CHARACTER*128 :: COMMAND
  INTEGER            :: NFILEV,ISTATV,SERROR,NK
  CHARACTER(LEN=128) :: FILENV
  CHARACTER(LEN=4)   :: NUMB
!
  NFILEV = 1400+NK
  WRITE(NUMB,'(I4.3)')NK
!
  WRITE(*,*)'#### SIMULACAO ####'
!
  IF(NSIMUL.EQ.0)THEN
     COMMAND='cd ./blackbox/; bash rodarBK.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > outputC')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.1)THEN
     COMMAND='cd ./simulador/; ./run > output_'//TRIM(ADJUSTL(NUMB))//'.out'
  END IF
!
  IF(NSIMUL.EQ.2)THEN
     COMMAND='cd ./simuladorRigido/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.3)THEN
     COMMAND='cd ./simul_comp/; bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.4)THEN
     COMMAND='cd ./co2injection/src/; bash rodarSimuladorPsh'
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     COMMAND='cd ./co2injection/src/mrb_dados/; ./run >output_'//TRIM(ADJUSTL(NUMB))//'.out'
  END IF
!
  IF(NSIMUL.EQ.5)THEN
     COMMAND='cd ./SIMULADOR_VISCOELASTICO/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > outputC')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.6)THEN
     COMMAND='cd ./twophaseflow/; matlab -nodisplay -nosplash -nodesktop -r "Simulator('
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)// ('); exit;"  > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.7)THEN
     COMMAND='cd ./twophaseflow/; octave --no-gui -q -H -W SimulatorOCT.m'
     COMMAND=TRIM(COMMAND)//' '//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)// (' ;')
  END IF
!
  WRITE(*,200)NK,COMMAND
  CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
  IF(NSIMUL.EQ.5)THEN
     FILENV = ('./SIMULADOR_VISCOELASTICO/exp')//&
          TRIM(ADJUSTL(NUMB))//('/error.out')
     OPEN(UNIT=NFILEV,FILE=FILENV,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTATV)
     IF(ISTATV.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILENV
        STOP
     END IF
     READ(NFILEV,100)SERROR
     CLOSE(UNIT=NFILEV)
     ERROSIMULADOR = SERROR
     IF(ERROSIMULADOR.EQ.0)THEN
        WRITE(*,300)NK,ERROSIMULADOR
     ELSE
        WRITE(*,400)NK,ERROSIMULADOR
     END IF
  END IF
!
100 FORMAT(I7)
200 FORMAT('RANK.....: ',I3,/,'COARSE...: ',A)
300 FORMAT('---SIMULACAO COARSE RANK: ',I3,'   -FINALIZADA CORRETAMENTE: ',I3,/)
400 FORMAT('------------------SIMULACAO COARSE RANK: ',I3,'   -ABORTADA: ',I3,/)
!
END SUBROUTINE SIMULADOR_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ REFERENCE SAMPLE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_AMOSTRA_C(NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: I,J,K,IN_FILE,ISTAT,NK,LNGT
  CHARACTER(LEN=128) :: FILEAMC,FAUX
  CHARACTER(LEN=256) :: CHARAUX,CHARAUX2
  CHARACTER(LEN=4)   :: NUMB
  CHARACTER(LEN=3)   :: FEXP
  IN_FILE = 22
  WRITE(NUMB,'(I4.3)')NK
!
  DO K=1,NDTYPE
     FILEAMC = FILEAM(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     LNGT = LEN(FILEAMC)
     FAUX = TRIM(FILEAMC)
     DO I=1,LNGT-2
        FEXP = FAUX(I:I+2)
        IF(FEXP.EQ.'exp')THEN
           GOTO 111
        END IF
     END DO
111  CONTINUE
     DO J=I,LNGT-2
        FEXP = FAUX(J:J)
        IF(FEXP.EQ.'/')THEN
           GOTO 115
        END IF
     END DO
115  CONTINUE
     FILEAMC = TRIM(ADJUSTL(FILEAMC(1:I-1)))//('exp')//TRIM(ADJUSTL(NUMB))//&
          TRIM(ADJUSTL(FILEAMC(J:LNGT)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     FILEAMC = TRIM(ADJUSTL(FILEAMC(2:LEN_TRIM(FILEAMC))))
     WRITE(*,*)'#### LEITURA DOS DADOS DA AMOSTRA (COARSE) ####'
     WRITE(*,*)FILEAMC
     OPEN(UNIT=IN_FILE,FILE=FILEAMC,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAMC
        STOP
     END IF
!
     DO I=1,NDADOS(K)
        READ(IN_FILE,*)(REFAMOS(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
!
     CLOSE(UNIT=IN_FILE)
  END DO
!
100 FORMAT(30(e15.8,2x))
!
END SUBROUTINE READ_AMOSTRA_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! CALCULA O LIKELIHOOD DA MALHA GROSSA !!!!!!!!!!!!!!!!
SUBROUTINE LIKELIHOOD_C()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER        :: K
  REAL, EXTERNAL :: NORMA_L2
!
  DO K=1,NDTYPE
     ERRORCC(K)=NORMA_L2(K)/NORMA_L2_REF(K)
     LRATIOC(K)=((1.0d0/SIGMA2C(K))*&
          (ERRORNC(K)-ERRORCC(K)))
     WRITE(*,100)K,ERRORCC(K)
  END DO
!
  TERRORCC=1.D0
  TLRATIOC=0.D0
  DO K=1,NDTYPE
     TERRORCC = ERRORCC(K)*TERRORCC
     TLRATIOC = LRATIOC(K)+TLRATIOC
  END DO
  TLRATIOC = EXP(TLRATIOC)
!
  RETURN
!
100 FORMAT('DATA no....:',I3,2X,'ERROR COARSE = ',E10.3)
!
END SUBROUTINE LIKELIHOOD_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MCMC_C(NN)
!
  USE VARIAVEIS
  USE RANDOM
  IMPLICIT NONE
!
  INTEGER   :: I,J,K,NN,TEST
  REAL(4)   :: AUX,AUX2
  LOGICAL   :: FIRST
!
  IF(SGN.EQ.1)THEN
     SGN  =0D0
     FIRST=.true.
     TEST =RANDOM_BINOMIAL2(1,AUX,FIRST)
  END IF
  IF(NN.EQ.1)THEN
     AUX     =1.D0
     AUX2    =1.D0
     TERRORNC=TERRORCC
     DO K=1,NDTYPE
        ERRORNC(K)=ERRORCC(K)
     END DO
     ACCEPT  =1
     FIRST   =.true.
     TEST    =RANDOM_BINOMIAL2(1,AUX,FIRST)
  ELSE
     FIRST =.false.
     AUX2  =TLRATIOC*TPRATIO
     AUX   =MIN(1.D0,AUX2)
     TEST  =RANDOM_BINOMIAL2(1,AUX,FIRST)
     ACCEPT=TEST
     IF(ACCEPT.GE.1)THEN
        TERRORNC=TERRORCC
        DO K=1,NDTYPE
           ERRORNC(K)=ERRORCC(K)
        END DO
     END IF
  ENDIF
!
  WRITE(*,*)'LRATIOC=',TLRATIOC,'PRATIO=',TPRATIO
  WRITE(*,*)'TESTE=',TEST,'ACCEPT=',ACCEPT
  IF(TEST.GT.0)THEN
     WRITE(*,*)'#########################################################'
     WRITE(*,*)'##-------------------- CAMPO ACEITO -------------------##'
     WRITE(*,*)'##-------------------- MALHA GROSSA -------------------##'
     WRITE(*,*)'#########################################################'
  ELSE
     NCHAIN = NCHAIN+1
     WRITE(*,*)'#########################################################'
     WRITE(*,*)'#################### CAMPO REJEITADO ####################'
     WRITE(*,*)'####################   MALHA GROSSA  ####################'
     WRITE(*,*)'#########################################################'
  END IF
!
  RETURN
!
END SUBROUTINE MCMC_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SIMULADOR USADO NA MALHA FINA !!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SIMULADOR_F(NK)
!
  USE VARIAVEIS, ONLY: NSIMUL,ERROSIMULADOR
  CHARACTER*128      :: COMMAND
  INTEGER            :: NFILEV,ISTATV,SERROR,NK
  CHARACTER(LEN=128) :: FILENV
  CHARACTER(LEN=4)   :: NUMB
!
  WRITE(NUMB,'(I4.3)')NK
  NFILEV = 1500+NK
!
  WRITE(*,*)'#### SIMULACAO ####'
!
  IF(NSIMUL.EQ.0)THEN
     COMMAND='cd ../blackbox/; bash rodarBK.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.1)THEN
     COMMAND='cd ../simulador/; ./run > output_'//TRIM(ADJUSTL(NUMB))//'.out'
  END IF
!
  IF(NSIMUL.EQ.2)THEN
     COMMAND='cd ../simuladorRigido/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.3)THEN
     COMMAND='cd ../simul_comp/;  bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.4)THEN
     COMMAND='cd ../co2injection/src/; bash rodarSimuladorPsh'
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     COMMAND='cd ../co2injection/src/mrb_dados/; ./run >output_'//TRIM(ADJUSTL(NUMB))//'.out'
  END IF
!
  IF(NSIMUL.EQ.5)THEN
     COMMAND='cd ../SIMULADOR_VISCOELASTICO/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.6)THEN
     COMMAND='cd ../twophaseflow/; matlab -nodisplay -nosplash -nodesktop -r "Simulator('
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)// ('); exit;"  > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  IF(NSIMUL.EQ.7)THEN
     COMMAND='cd ../twophaseflow/; octave --no-gui -q -H -W SimulatorOCT.m '
     COMMAND=TRIM(COMMAND)//' '//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)// (' ;')
  END IF
!
  
  WRITE(*,*)COMMAND
  CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
  IF(NSIMUL.EQ.5)THEN
     FILENV = ('../SIMULADOR_VISCOELASTICO/exp')//&
          TRIM(ADJUSTL(NUMB))//('/error.out')
     WRITE(*,*)NUMB, FILENV
     OPEN(UNIT=NFILEV,FILE=FILENV,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTATV)
     IF(ISTATV.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILENV
        STOP
     END IF
     READ(NFILEV,100)SERROR
     CLOSE(UNIT=NFILEV)
     ERROSIMULADOR = SERROR
     IF(ERROSIMULADOR.EQ.0)THEN
        WRITE(*,300)NK,ERROSIMULADOR
     ELSE
        WRITE(*,400)NK,ERROSIMULADOR
     END IF
  END IF
!
100 FORMAT(I7)
200 FORMAT('RANK.....: ',I3,/,'FINE...: ',A)
300 FORMAT('---SIMULACAO FINE RANK: ',I3,'   -FINALIZADA CORRETAMENTE: ',I3,/)
400 FORMAT('------------------SIMULACAO FINE RANK: ',I3,'   -ABORTADA: ',I3,/)
!
END SUBROUTINE SIMULADOR_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! CALCULA O LIKELIHOOD NA MALHA FINA !!!!!!!!!!!!!!!!!
SUBROUTINE LIKELIHOOD_F()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER          :: K,I,J
  REAL   ,EXTERNAL :: NORMA_L2
!
  WRITE(*,*)'####################################################'
  WRITE(*,*)'##================= LIKELIHOOD ===================##'
  DO K=1,NDTYPE
     ERRORCF(K)=NORMA_L2(K)/NORMA_L2_REF(K)
     LRATIOF(K)=((1.0d0/SIGMA2F(K))*(ERRORNF(K)-ERRORCF(K)))
     write(*,100)K,ERRORCF(K)
  END DO
  WRITE(*,*)'####################################################'
!
  TERRORCF=1.D0
  TLRATIOF=0.D0
  DO K=1,NDTYPE
     TERRORCF = ERRORCF(K)*TERRORCF
     TLRATIOF = LRATIOF(K)+TLRATIOF
  END DO
!
  TLRATIOF = EXP(TLRATIOF)
!
  RETURN
!
100 FORMAT('DATAr no...:',I3,2X,'RELATIVE ERROR..:',E10.3)
!
END SUBROUTINE LIKELIHOOD_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MCMC_F(NN,NK)
!
  USE VARIAVEIS
  USE RANDOM
  IMPLICIT NONE
!
  INTEGER :: I,J,K,NN,TEST,NK,NT
  REAL    :: AUX
  REAL    :: AUX2
  LOGICAL :: FIRST
!
  NT = 1
  IF(SGN.EQ.1)THEN
     SGN     =0D0
     FIRST   =.true.
     TEST    =RANDOM_BINOMIAL2(1,AUX,FIRST)
  END IF
!
  IF(NN.EQ.1)THEN
     AUX     =1.D0
     AUX2    =1.D0
     TERRORNF=TERRORCF
     DO K=1,NDTYPE
        ERRORNF(K)=ERRORCF(K)
     END DO
     ACCEPT  =1
     FIRST   =.true.
     TEST    =RANDOM_BINOMIAL2(NT,AUX,FIRST)
  ELSE
     FIRST   =.true.
     IF(NSTAGE.EQ.1)THEN
        AUX2 =TLRATIOF*TPRATIO
     ELSE
        AUX2 =0.D0
        DO K=1,NDTYPE
           AUX2=AUX2+(0.5D0*SIGMA2C(K))*(ERRORCC(K)-ERRORNC(K))&
               + (0.5D0*SIGMA2F(K))*(ERRORNF(K)-ERRORCF(K))
        END DO
        AUX2=EXP(AUX2)
     END IF
     AUX =MIN(1.D0,AUX2)
     TEST=RANDOM_BINOMIAL2(NT,AUX,FIRST)
     ACCEPT=TEST
     WRITE(*,100)AUX
     IF(ACCEPT.GE.1)THEN
        TERRORNF=TERRORCF
        DO K=1,NDTYPE
           ERRORNF(K)=ERRORCF(K)
        END DO
     END IF
  ENDIF
!
  IF(TEST.GT.0)THEN
     IF(NN.EQ.1)THEN
        NCHAIN = 1
     ELSE
        CALL CHAIN_COUNTER(NCHAIN,NK)
        NCHAIN = 1
     END IF
     WRITE(*,*)'############################'
     WRITE(*,*)'#------ CAMPO ACEITO ------#'
     WRITE(*,*)'#------- MALHA FINA -------#'
     WRITE(*,*)'############################'
     LOGERRORCF(CONTADORC) = LOG(TERRORCF)
  ELSE
     IF(NSTAGE.EQ.1)NCHAIN = NCHAIN+1
     WRITE(*,*)'############################'
     WRITE(*,*)'###### CAMPO REJEITADO #####'
     WRITE(*,*)'######## MALHA FINA ########'
     WRITE(*,*)'############################'
     LOGERRORCF(CONTADORC) = LOGERRORCF(CONTADORC-1)
  END IF
!
  RETURN
!
100 FORMAT('()()()()()()()()()()()()()()()()()()() ',/,&
           '(-) PROBABILITY........=',E10.3,/,&
           '()()()()()()()()()()()()()()()()()()()')
!
END SUBROUTINE MCMC_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COPYFILES(NN,NOMEP,NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: I,J,K,M,NN,IV,NK
  CHARACTER(LEN=128) :: NAME,CHARAUX,CHARAUX2,NAMEV
  CHARACTER(LEN=128) :: NAMEOUT
  CHARACTER(LEN=5)   :: CHA
  CHARACTER(LEN=4)   :: EXT,NUMB
  CHARACTER(LEN=6)   :: EXT2
  CHARACTER(LEN=5)   :: CHARAC
  CHARACTER(LEN=7)   :: CHARRANK
  CHARACTER(LEN=5)   :: C,D,SAMP,CVAR
  CHARACTER(LEN=256) :: COMMAND
  CHARACTER(LEN=128),DIMENSION(1:10):: NOMEP
  D=TRIM('_')
  WRITE(CHARAC,113)NK
  WRITE(NUMB,'(I4.3)')NK
  CHARAC=ADJUSTL(CHARAC)
  CHARRANK = ADJUSTL(('_RK')//TRIM(CHARAC))
!
  IF(ACCEPT.GT.0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO K=1,NPRIORR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(NLPRIOR(K).NE.0.OR.GERATIPO(K).EQ.33)THEN
           CALL SAVE_LAMBDAS(CONT,K,NK)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        EXT2='_0.dat'
        IF(NSIMUL.EQ.0)THEN
           COMMAND='../blackbox/exp/fields/'
        END IF
        IF(NSIMUL.EQ.1)THEN
           COMMAND='../simulador/fields/'
        END IF
        IF(NSIMUL.EQ.2)THEN
           COMMAND='../simuladorRigido/exp/fields/'
        END IF
        IF(NSIMUL.EQ.3)THEN
           COMMAND='../simul_comp/exp/fields/'
        END IF
        IF(NSIMUL.EQ.4)THEN
           COMMAND='../co2injection/src/mrb_dados/fields/'
        END IF
        IF(NSIMUL.EQ.5)THEN
           COMMAND='../SIMULADOR_VISCOELASTICO/exp/fields/'
        END IF
        IF(NSIMUL.EQ.6)THEN
           COMMAND='../twophaseflow/exp/fields/'
        END IF
        IF(NSIMUL.EQ.7)THEN
           COMMAND='../twophaseflow/exp/fields/'
        END IF

        COMMAND=TRIM(COMMAND)//TRIM(NOMEP(K))//TRIM(EXT2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.2)THEN
           COMMAND=ADJUSTL(TRIM(COMMAND))
           CALL READPERM(NPRIORR,COMMAND,LOGPERMN,K)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.4)THEN
           COMMAND=ADJUSTL(TRIM(COMMAND))
           CALL READPERM(NPRIORR,COMMAND,LOGPERMN,K)
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.5)THEN
           COMMAND=ADJUSTL(TRIM(COMMAND))
           CALL READPERM(NPRIORR,COMMAND,LOGPERMN,K)
        ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.0)THEN
           COMMAND=ADJUSTL(TRIM(COMMAND))
        ENDIF             
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        EXT='.dat'
        M=0
        WRITE(C,113)M
        C=ADJUSTL(C)
        NAME=TRIM(FILE_OUT(K))//('exp')//TRIM(ADJUSTL(NUMB))//&
             TRIM(FILE_EXT(K))//TRIM(C)//TRIM(EXT)
        IF(GERATIPO(K).EQ.1)THEN
           NAME=ADJUSTL(TRIM(NAME(4:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.2)THEN
           NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.&
             GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.&
             GERATIPO(K).EQ.34.OR.GERATIPO(K).EQ.13.OR.&
             GERATIPO(K).EQ.14)THEN
           NAME=ADJUSTL(TRIM(NAME(7:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.4)THEN
           NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.5)THEN
           NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.0)THEN
           NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.6)THEN
           NAME=ADJUSTL(TRIM(NAME(7:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.7)THEN
           NAME=ADJUSTL(TRIM(NAME(7:LEN_TRIM(NAME))))
        ENDIF
        IF(GERATIPO(K).EQ.8)THEN
           NAME=ADJUSTL(TRIM(NAME(7:LEN_TRIM(NAME))))
        ENDIF
        WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHA=TRIM(ADJUSTL('mv  '))
        WRITE(CVAR,113)K
        CVAR=ADJUSTL(CVAR)    
! 
        IF(GERATIPO(K).EQ.33)THEN
           NAMEV=NAME(1:LEN_TRIM(NAME)-5)
           DO IV=0,NZMIN(K)-1
              WRITE(SAMP,113)IV
              SAMP=ADJUSTL(SAMP)
              WRITE(C,113)CONT
              C=ADJUSTL(C)
              NAMEOUT=TRIM(' ./select_fields/field_v')//&
                   TRIM(CVAR)//TRIM(D)//TRIM(NAME_OUT)//TRIM(CHARRANK)&
                   //TRIM(D)//TRIM('S')//TRIM(C)//TRIM(D)
              COMMAND=(CHA)//TRIM(NAMEV)//TRIM(SAMP)//TRIM(EXT)
              COMMAND=TRIM(COMMAND)//NAMEOUT
              COMMAND=ADJUSTL(TRIM(COMMAND))
              COMMAND=TRIM(COMMAND)//TRIM(SAMP)//TRIM(EXT)
              COMMAND=ADJUSTL(TRIM(COMMAND))
              WRITE(*,*)'VARIAVEL: ',CVAR,COMMAND
              CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
           END DO
        ELSE
           NAMEOUT=TRIM(' ./select_fields/field_v')//&
                TRIM(CVAR)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
           COMMAND=CHA//NAME
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
           WRITE(*,*)'VARIAVEL:',CVAR,COMMAND
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
        END IF
     END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO K=1,NDTYPE
        CHA=TRIM(ADJUSTL('mv  '))
        WRITE(C,113)K
        C=ADJUSTL(C)     
        NAMEOUT=TRIM(' ./select_prod/prod_D')//&
             TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
             TRIM(CHARRANK)//TRIM(D)
        WRITE(C,113)CONT
        C=ADJUSTL(C)
!
        COMMAND=CHA//FILEAM(K)
        COMMAND=TRIM(COMMAND)//NAMEOUT
        COMMAND=ADJUSTL(TRIM(COMMAND))
        COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
        COMMAND=ADJUSTL(TRIM(COMMAND))
!           
        WRITE(*,*)COMMAND
!
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END DO
!
!!! SALVANDO O ERRO DO FLUXO (PRODUCAO) ACEITOS !!!!!!!!!
     ERRORF(CONT)=TERRORCF
     CALL SAVE_ERROR(CONT,TERRORCF,NK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO K=1,NPRIORR
        IF(GERATIPO(K).EQ.0)THEN
           CHARAUX=TRIM('mv ../gera_LABTRANGEORW3/out/xpsinew.dat ')
           CHARAUX2=TRIM(' ../gera_LABTRANGEORW3/out/xpsi.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.1)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.3)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_RW/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_RW/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_RW/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.13)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_KL3D/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_KL3D/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_KL3D/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.14)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_KL3D_2/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_KL3D_2/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_KL3D_2/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.6)THEN
           CHARAUX=TRIM('mv ../gera_KL/MVN/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/MVN/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/MVN/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.7)THEN
           CHARAUX=TRIM('mv ../gera_KL/MVN1/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/MVN1/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/MVN1/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.8)THEN
           CHARAUX=TRIM('mv ../gera_KL/MVN2/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/MVN2/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)
!
           CHARAUX =TRIM(' ../gera_KL/MVN2/out/theta')//&
           TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.31)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_RW1/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_RW1/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)     
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_RW1/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.34)THEN
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_RW2/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_RW2/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)     
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_RW2/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.32)THEN
           OLDLBD(K) =NEWLBD(K)
           YOLDLBD(K)=YNEWLBD(K)
           NUM_AVETOLD(K)=NUM_AVET(K)
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_RW3D/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_RW3D/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)     
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_RW3D/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.33)THEN
           YOLD=YNEW
           OLDLBD(K)=NEWLBD(K)
           YOLDLBD(K)=YNEWLBD(K)
           CHARAUX=TRIM('mv ../gera_KL/FORTRAN_LBD/out/thetanew')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           CHARAUX2=TRIM(' ../gera_KL/FORTRAN_LBD/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CHA=TRIM(ADJUSTL('cp  '))
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./select_thetas/theta_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONT
           C=ADJUSTL(C)     
!
           CHARAUX =TRIM(' ../gera_KL/FORTRAN_LBD/out/theta')//&
                TRIM(ADJUSTL(NUMB))//('.dat ')
           COMMAND=CHA//CHARAUX
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.4)THEN
           CHARAUX=TRIM('mv ../gera_LABTRANGEORW/out/niveisnew.dat ')
           CHARAUX2=TRIM(' ../gera_LABTRANGEORW/out/niveis.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(GERATIPO(K).EQ.5)THEN
           CHARAUX=TRIM('mv ../gera_LABTRANGEORW2/out/niveisnew.dat ')
           CHARAUX2=TRIM(' ../gera_LABTRANGEORW2/out/niveis.dat')
           COMMAND=TRIM(CHARAUX)//TRIM(CHARAUX2)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!
        ENDIF
     END DO
!              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     IF(NPRINT.EQ.1)THEN
!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        DO K=1,NPRIORR
           EXT='.dat'
           M=0
           WRITE(C,113)M
           C=ADJUSTL(C)
           NAME=TRIM(FILE_OUT(K))//TRIM(C)//TRIM(EXT)
           IF(GERATIPO(K).EQ.1)THEN
              NAME=ADJUSTL(TRIM(NAME(4:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.2)THEN
              NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.&
                GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.&
                GERATIPO(K).EQ.34.OR.GERATIPO(K).EQ.13.OR.&
                GERATIPO(K).EQ.14)THEN
              NAME=ADJUSTL(TRIM(NAME(4:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.4)THEN
              NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.5)THEN
              NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.0)THEN
              NAME=ADJUSTL(TRIM(NAME(1:LEN_TRIM(NAME))))
           ENDIF
           IF(GERATIPO(K).EQ.6.OR.GERATIPO(K).EQ.7.OR.&
                GERATIPO(K).EQ.8)THEN
              NAME=ADJUSTL(TRIM(NAME(4:LEN_TRIM(NAME))))
           ENDIF
!
           WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
           CHA=TRIM(ADJUSTL('mv  '))
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./reject_fields/field_v')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
           WRITE(C,113)CONTREJ
           C=ADJUSTL(C)
!
           COMMAND=CHA//NAME
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
!            
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
        END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO K=1,NDTYPE
           WRITE(C,113)K
           C=ADJUSTL(C)     
           NAMEOUT=TRIM(' ./reject_prod/prod_D')//&
                TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
                TRIM(CHARRANK)//TRIM(D)
!
           COMMAND=CHA//FILEAM(K)
           COMMAND=TRIM(COMMAND)//NAMEOUT
           COMMAND=ADJUSTL(TRIM(COMMAND))
           COMMAND=TRIM(COMMAND)//TRIM(C)//TRIM(EXT)
           COMMAND=ADJUSTL(TRIM(COMMAND))
!
           WRITE(*,*)COMMAND
!
           CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
!  
        END DO
!
     ENDIF
!
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GUARDAR NIVEL ESCOLHIDO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IF(GERATIPO(K).EQ.4)THEN
!     CALL SAVE_LEVEL(NAME_OUT,ACCEPT,GERATIPO(K))
!  ENDIF
!  IF(GERATIPO(K).EQ.5)THEN
!     CALL SAVE_LEVEL(NAME_OUT,ACCEPT,GERATIPO(K))
!  ENDIF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
113 FORMAT(I5)
!
END SUBROUTINE COPYFILES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROTINA PARA SALVAR OS ERROS EM CADA ELEMENTO DA CADEIA !!!!!!!!
SUBROUTINE SAVE_ERROR(NUMB,ER,NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: OUT_FILE,ISTAT,K,NK
  CHARACTER(LEN=128) :: FOUT
  CHARACTER(LEN=5)   :: CHARAC
  INTEGER            :: NUMB
  REAL               :: ER
  LOGICAL            :: TFILE
!
  OUT_FILE = 1600+NK
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)     
  FOUT = TRIM(FILERROR)//('_RK')//TRIM(CHARAC)//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
  INQUIRE(FILE=FOUT,EXIST=TFILE)
!  WRITE(*,*)'ARQUIVO EXISTE:',TFILE
  IF(TFILE)THEN
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='OLD',&
          ACTION='WRITE',IOSTAT=ISTAT,POSITION='APPEND')
  ELSE
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='NEW',&
          ACTION='WRITE',IOSTAT=ISTAT)
  END IF
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FOUT
     STOP
  END IF
!
  WRITE(OUT_FILE,FMT=100)NUMB,ER,(ERRORCF(K),K=1,NDTYPE)
!
  CLOSE(OUT_FILE)
!
100 FORMAT(I7,1X,E15.5,1X,200E15.5)
112 FORMAT(I5)
!
END SUBROUTINE SAVE_ERROR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ REFERENCE SAMPLE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_AMOSTRA_F(NK)
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: I,J,K,IN_FILE,ISTAT,NK,LNGT
  CHARACTER(LEN=128) :: FILEAMC,FAUX
  CHARACTER(LEN=4)   :: NUMB
  CHARACTER(LEN=3)   :: FEXP
!
  WRITE(NUMB,'(I4.3)')NK
!
  IN_FILE = 1700+NK
!
  DO K=1,NDTYPE
     FILEAMC = FILEAM(K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     LNGT = LEN(FILEAMC)
     FAUX = TRIM(FILEAMC)
     DO I=1,LNGT-2
        FEXP = FAUX(I:I+2)
        IF(FEXP.EQ.'exp')THEN
           GOTO 111
        END IF
     END DO
111  CONTINUE
     DO J=I,LNGT-2
        FEXP = FAUX(J:J)
        IF(FEXP.EQ.'/')THEN
           GOTO 115
        END IF
     END DO
115  CONTINUE
     FILEAMC = TRIM(ADJUSTL(FILEAMC(1:I-1)))//('exp')//TRIM(ADJUSTL(NUMB))//&
          TRIM(ADJUSTL(FILEAMC(J:LNGT)))
     FILEAM(K) = TRIM(FILEAMC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE(*,*)'#### LEITURA DOS DADOS DA AMOSTRA (FINE) ####'
     WRITE(*,*)FILEAMC
     OPEN(UNIT=IN_FILE,FILE=FILEAMC,STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAMC
        STOP
     END IF
!
     DO I=1,NDADOS(K)
        READ(IN_FILE,*)(REFAMOS(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
  !
     CLOSE(UNIT=IN_FILE)
  END DO
!
100 FORMAT(30(E15.8,2x))
!
END SUBROUTINE READ_AMOSTRA_F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROTINA CONTAR NUMERO DE VEZES QUE UM ELEMENTO SE !!!!!!
! REPETE NA CADEIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CHAIN_COUNTER(NUMB,NK)
!
  USE VARIAVEIS, ONLY: CONT,NAME_OUT
  IMPLICIT NONE
!
  INTEGER            :: OUT_FILE,ISTAT,NK
  CHARACTER(LEN=128) :: FOUT
  INTEGER            :: NUMB
  CHARACTER(LEN=5)   :: CHARAC
  LOGICAL            :: TFILE
!
  OUT_FILE = 1800+NK
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)
  FOUT = TRIM('./out/nchain_')//TRIM(NAME_OUT)//TRIM('_RK')//&
       TRIM(ADJUSTL(CHARAC))//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
  INQUIRE(FILE=FOUT,EXIST=TFILE)
!
  IF(TFILE)THEN
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='OLD',&
          ACTION='WRITE',IOSTAT=ISTAT,POSITION='APPEND')
  ELSE
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='NEW',&
          ACTION='WRITE',IOSTAT=ISTAT)
  END IF
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FOUT
     STOP "CHAIN_COUNTER SUBROUTINE"
  END IF
!
  WRITE(UNIT=OUT_FILE,FMT=100)CONT,NUMB
!
  CLOSE(OUT_FILE)      
!
100 FORMAT(I7,I7)
112 FORMAT(I5)
!
END SUBROUTINE CHAIN_COUNTER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ROTINA PARA SALVAR OS COMPRIMENTOS DE CORRELACCAO !!!!!
! SELECIONADOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SAVE_LAMBDAS(NUMB1,NUMB2,NK)
!
  USE VARIAVEIS, ONLY: NUM_AVET,NAME_OUT,NEWLBD
  IMPLICIT NONE
!
  INTEGER            :: OUT_FILE,ISTAT,NK
  CHARACTER(LEN=128) :: FOUT
  CHARACTER(LEN=5)   :: CHARAC
  INTEGER            :: NUMB1,NUMB2
  LOGICAL            :: TFILE
!
  OUT_FILE = 1900+NK
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)
  FOUT = TRIM('./out/lambda_')//TRIM(NAME_OUT)//('_RK')//&
       TRIM(CHARAC)//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
  INQUIRE(FILE=FOUT,EXIST=TFILE)
!
  IF(TFILE)THEN
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='OLD',&
          ACTION='WRITE',IOSTAT=ISTAT,POSITION='APPEND')
  ELSE
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='NEW',&
          ACTION='WRITE',IOSTAT=ISTAT)
  END IF
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FOUT
     STOP
  END IF
!
  WRITE(UNIT=OUT_FILE,FMT=100)NUMB1,NUMB2,&
       NUM_AVET(NUMB2),NEWLBD(NUMB2)
!
  CLOSE(OUT_FILE)      
!
100 FORMAT(3I7,E10.3)
112 FORMAT(I5)
!
END SUBROUTINE SAVE_LAMBDAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INTERPOLATION_PREP(LBD,K)
!
  USE VARIAVEIS, ONLY:  CLBDS,XLBDS,NLPRIOR,NLBDPOSI
  USE VARIAVEIS, ONLY:  NLBD,NLBDVET
  REAL               :: LBD
  INTEGER            :: K,I,J,LIMA,LIMB
!
  write(*,*)'(((((((((((((((((())))))))))))))))))'
  IF(LBD.LE.CLBDS(1))THEN
     LBD=CLBDS(1)
     NLBDPOSI=1
  ELSE IF(LBD.GE.CLBDS(NLPRIOR(K)))THEN
     LBD=CLBDS(NLPRIOR(K))
     NLBDPOSI=NLPRIOR(K)
  ELSE
     DO I=1,NLPRIOR(K)-1
        IF(LBD.GT.CLBDS(I).AND.LBD.LE.CLBDS(I+1))THEN
           NLBDPOSI=I
        END IF
     END DO
  END IF
!
  LIMB=NLBD(K)/2
  LIMA=NLPRIOR(K)-LIMB
!
  IF(NLBDPOSI.LE.LIMB)THEN
     DO I=1,NLBD(K)
        NLBDVET(I)=I
     END DO
  ELSE IF(NLBDPOSI.GE.LIMA)THEN
     DO I=1,NLBD(K)
        NLBDVET(I)=NLPRIOR(K)-NLBD(K)+I
     END DO
  ELSE
     DO I=1,NLBD(K)
        NLBDVET(I)=NLBDPOSI-LIMB+I
     END DO
  END IF
!
  write(*,*)'(((((((((((((((((())))))))))))))))))'
  WRITE(*,*)'LAMBDAS USADOS:'
  DO I=1,NLBD(K)
     WRITE(*,*)CLBDS(NLBDVET(I))
  END DO
  write(*,*)'(((((((((((((((((())))))))))))))))))'
END SUBROUTINE INTERPOLATION_PREP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOADTHETAS(K,DM,N,NCONT,NK)
!
  USE VARIAVEIS, ONLY: XRMVN_VET,NAME_OUT
  IMPLICIT NONE
!
  INTEGER            :: IN_FILE,ISTAT,I,J,K,N,DM,NCONT,II
  INTEGER            :: NK
  CHARACTER(LEN=128) :: NOME,FILEINPUT
  CHARACTER(LEN=5)   :: C,D,M,CHARAC
  CHARACTER(LEN=7)   :: CHARRANK
  CHARACTER(LEN=4)   :: EXT
!
  IN_FILE = 2000+NK
  WRITE(*,205)NAME_OUT
  D   = TRIM('_')
  EXT = '.dat'
  WRITE(C,111)K
  C = ADJUSTL(C)
  WRITE(CHARAC,111)NK
  CHARAC   = ADJUSTL(CHARAC)
  CHARRANK = ADJUSTL(('_RK')//TRIM(CHARAC))
!
  II=0
  DO I=NCONT-N+1,NCONT-1
     II=II+1
     WRITE(M,111)I
     M=ADJUSTL(M)
     NOME = ADJUSTL(TRIM(' ./select_thetas/theta_v')//&
          TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
          TRIM(CHARRANK)//TRIM(D)//&
          TRIM(M)//TRIM(EXT))
     FILEINPUT = ADJUSTL(NOME)
     WRITE(*,206)NK,FILEINPUT
!
     OPEN(UNIT=IN_FILE,FILE=FILEINPUT,ACTION='READ',&
          STATUS='OLD',FORM='FORMATTED',IOSTAT=ISTAT)
!
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEINPUT
        STOP
     END IF
!
     DO J=1,DM
        READ(IN_FILE,112)XRMVN_VET(J,II)
!        WRITE(*,*)J,II,XRMVN_VET(J,II)
     END DO
     CLOSE(UNIT=IN_FILE)
  END DO
!
112 FORMAT(E15.8)
111 FORMAT(I5)
205 FORMAT('NOME BASE DOS ARQUIVOS.... : ',A)
206 FORMAT('RANK..:',I5,2x,' - LEITURA DO ARQUIVO.........: ',A)
!
END SUBROUTINE LOADTHETAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOADTHETASDE(NK,K,DM,NVETCONT,NPROC)
!
  USE VARIAVEIS, ONLY: XRMVN,NAME_OUT,MATDE
  IMPLICIT NONE
!
  INTEGER                 :: IN_FILE,ISTAT,NCONT,II
  INTEGER                 :: NK,NPROC,I,J,K,N,DM,NKRANK
  INTEGER,DIMENSION(NPROC):: NVETCONT
  CHARACTER(LEN=128)      :: NOME,FILEINPUT
  CHARACTER(LEN=7)        :: CHARRANK
  CHARACTER(LEN=5)        :: C,D,M,CHARAC
  CHARACTER(LEN=4)        :: EXT
!
  IN_FILE = 2200+NK
!  WRITE(*,205)NAME_OUT
  D  =TRIM('_')
  EXT='.dat'
  WRITE(C,111)K
  C  =ADJUSTL(C)
  DO NKRANK=0,NPROC-1
     FILEINPUT = ''
     WRITE(M,111)NVETCONT(NKRANK+1)-1
     M=ADJUSTL(M)
     WRITE(CHARAC,111)NKRANK
     CHARAC=ADJUSTL(CHARAC)
     CHARRANK = ADJUSTL(('_RK')//TRIM(CHARAC))
!
     NOME = ADJUSTL(TRIM(' ./select_thetas/theta_v')//&
          TRIM(C)//TRIM(D)//TRIM(NAME_OUT)//&
          TRIM(CHARRANK)//TRIM(D)//&
          TRIM(M)//TRIM(EXT))
     FILEINPUT = ADJUSTL(NOME)
     WRITE(*,206)NK,FILEINPUT
!
     OPEN(UNIT=IN_FILE,FILE=FILEINPUT,ACTION='READ',&
          STATUS='OLD',FORM='FORMATTED',IOSTAT=ISTAT)
!
     IF(ISTAT.NE.0)THEN
        WRITE(*,207)FILEINPUT,ISTAT
        STOP 44444
     END IF
!
     DO J=1,DM
        READ(IN_FILE,112)XRMVN(J)
!        WRITE(*,*)J,XRMVN(J)
     END DO
     MATDE(1:DM,NKRANK+1) = XRMVN
     CLOSE(UNIT=IN_FILE)
  END DO
!
  RETURN
!
112 FORMAT(E15.8)
111 FORMAT(I5)
205 FORMAT('NOME BASE DOS ARQUIVOS....: ',A)
206 FORMAT('RANK:',I4,' <-->  LEITURA DO ARQUIVO....: ',A)
207 FORMAT('ERROR ON OPENING INPUT FILE IN LOADTHETASDE: ',A,/,&
           'IOSTAT:',I4)
!
END SUBROUTINE LOADTHETASDE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GERA_AMTHETA(MK,TIPO,K,NK)
!
  USE RANDOM
  USE VARIAVEIS, ONLY: THETAN
  IMPLICIT NONE
!
  INTEGER            :: I,N,MK,TIPO,K,NK
  INTEGER            :: ISTAT,OUT_FILE
  CHARACTER(LEN=128) :: FILETHE
  CHARACTER(LEN=4)   :: NUMB
  REAL               :: TESTE
!
  WRITE(NUMB,'(I4.3)')NK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  OUT_FILE = 500+NK
  IF(TIPO.EQ.0) FILETHE='../gera_LABTRANGEORW3/out/xpsi.dat'
  IF(TIPO.EQ.1) FILETHE=('../gera_KL/FORTRAN/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.3) FILETHE=('../gera_KL/FORTRAN_RW/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.13) FILETHE=('../gera_KL/FORTRAN_KL3D/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.14) FILETHE=('../gera_KL/FORTRAN_KL3D_2/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.31)FILETHE=('../gera_KL/FORTRAN_RW1/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.34)FILETHE=('../gera_KL/FORTRAN_RW2/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.32)FILETHE=('../gera_KL/FORTRAN_RW3D/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.33)FILETHE=('../gera_KL/FORTRAN_LBD/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.6) FILETHE=('../gera_KL/MVN/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.7) FILETHE=('../gera_KL/MVN1/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.8) FILETHE=('../gera_KL/MVN2/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
!
  OPEN(UNIT=OUT_FILE,FILE=FILETHE,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE IN GERA_AMTHETA: ',FILETHE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,MK
     WRITE(OUT_FILE,111)THETAN(K,I)
  ENDDO
!
  CLOSE(UNIT=OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE GERA_AMTHETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DREAM_METHOD(K,NK,NPROC,DIM,PCR)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE MPI
  USE VARIAVEIS, ONLY:  SIG,MKL,XRMVN,MATDE,XRMVN_MEAN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,CONT,VETCONT
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,DIMR,SIGDE
  USE VARIAVEIS, ONLY:  CONTADORC,JVET,NID
  USE VARIAVEIS, ONLY:  NDELTA,NFREQ,CR_DREAM,NCR,JUMP
  USE VARIAVEIS, ONLY:  VARPROP
!
  REAL                  :: GAMMA_DREAM,AUX,CRDREAM
  REAL                  :: GAMMA_DRAUX
  REAL                  :: VUNIF,B,A,BSTAR,COEFF
  INTEGER               :: K,I,J,N,DIM,NK,NPROC,DESTART
  REAL   ,DIMENSION(DIM):: XVET,XSUM
  LOGICAL               :: FIRST
  INTEGER,DIMENSION(2)  :: R
  INTEGER,DIMENSION(DIM):: DVET
  INTEGER, EXTERNAL     :: INTMIN
  INTEGER               :: TEST,DLINHA
  INTEGER               :: IDCR,DELTA,ID
  REAL   ,DIMENSION(NCR):: PCR,CR
  REAL   ,DIMENSION(DIM):: DXVET,STDX
  REAL   ,EXTERNAL      :: GAMMA_DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DXVET   = 0.0
  DVET    = 0
  FIRST   = .TRUE.
  B       = 1.0E-01
  A       = 0.0E+00
  BSTAR   = 1.0E-06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL RANDOM_SEED() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NDELTA.GT.NPROC)NDELTA=NPROC-1
  IF(NDELTA.EQ.0)    NDELTA=1
  DELTA = UNIFDISC(1,NDELTA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NCR
     CR(I) = REAL(I)/REAL(NCR)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IDCR = MULTINOM(NCR,PCR)
  CRDREAM = CR(IDCR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DESTART = INTMIN(VETCONT,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(DESTART.GT.0)THEN
!     SIG(K) = 0.0E+00
     MATDE  = 0.0E+00
!! carrega matriz MATDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL LOADTHETASDE(NK,K,DIM,VETCONT,NPROC)
     STDX = STDV(MATDE,DIM,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     XSUM = 0.0
     DO J=1,DELTA
        R    = CHOOSER(NK+1,NPROC)
        XSUM = XSUM + (MATDE(1:DIM,R(1))-MATDE(1:DIM,R(2)))
        write(*,*)'NK',NK,'D',J,'paresR',R
     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! escolha de quais dimensoes serao atualizadas !!!!!!!!!
!! obs: se CRDREAM = 1.0 todas as dimensoes sao atualiz !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DLINHA = DIM
     ID     = 0
     DO J=1,DIM
        FIRST  = .TRUE.
        TEST   = RANDOM_BINOMIAL2(1,1.0-CRDREAM,FIRST)
        DLINHA = DLINHA-TEST
        IF(TEST.EQ.1)THEN
           THETAN(K,J) = MATDE(J,NK+1)
        ELSE
           ID       = ID+1
           DVET(ID) = J
        END IF
     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(DLINHA.EQ.0)THEN
        DLINHA  = 1
        DVET(1) = UNIFDISC(1,DIM)
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! proposal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     GAMMA_DREAM = 2.38/SQRT(2.0*REAL(DLINHA*DELTA))
     GAMMA_DREAM = SIG(K)
     GAMMA_DREAM = GAMMA_DREAM_METHOD(SIG(K),CONTADORC,&
          NFREQ(K),VARPROP(K))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE(*,100)DELTA,DLINHA,NCR,GAMMA_DREAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO I=1,DLINHA
        J = DVET(I)
        CALL RANDOM_NUMBER(VUNIF)
        VUNIF       = (B-A)*VUNIF+A
        COEFF       = (1.0-VUNIF)*GAMMA_DREAM
        DXVET(J)    = COEFF*XSUM(J) + BSTAR*RANDOM_NORMAL()
        THETAN(K,J) = MATDE(J,NK+1) + DXVET(J)
     END DO
     JUMP = COEFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL GERA_AMTHETA(MKL(K),GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     JVET(K,IDCR) = JVET(K,IDCR) + SUMSQR(DXVET,STDX,DIM,DLINHA)
     NID(K,IDCR)  = NID(K,IDCR) + 1
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT('NUMERO DE PARES USADOS........: ',I4,/,&
           'NUMERO DE DIMENSOES TROCADAS..: ',I4,/,&
           'NUMERO NCR....................: ',I4,/,&
           'LAMBDA METODO DREAM...........: ',F6.3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DREAM_METHOD_MOD(K,NK,NPROC,DIM,PCR)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE MPI
  USE VARIAVEIS, ONLY:  SIG,MKL,XRMVN,MATDE,XRMVN_MEAN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,CONT,VETCONT
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,DIMR,SIGDE
  USE VARIAVEIS, ONLY:  CONTADORC,JVET,NID
  USE VARIAVEIS, ONLY:  NDELTA,NFREQ,CR_DREAM,NCR,JUMP
  USE VARIAVEIS, ONLY:  VARPROP
!
  REAL                  :: GAMMA_DREAM,AUX,CRDREAM
  REAL                  :: GAMMA_DRAUX
  REAL                  :: VUNIF,B,A,BSTAR,COEFF
  INTEGER               :: K,I,J,N,DIM,NK,NPROC,DESTART
  REAL   ,DIMENSION(DIM):: XVET,XSUM
  LOGICAL               :: FIRST
  INTEGER,DIMENSION(2)  :: R
  INTEGER,DIMENSION(DIM):: DVET
  INTEGER, EXTERNAL     :: INTMIN
  INTEGER               :: TEST,DLINHA
  INTEGER               :: IDCR,DELTA,ID
  REAL   ,DIMENSION(NCR):: PCR,CR
  REAL   ,DIMENSION(DIM):: DXVET,STDX
  REAL   ,EXTERNAL      :: GAMMA_DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DXVET   = 0.0
  DVET    = 0
  FIRST   = .TRUE.
  B       = 1.0E-01
  A       = 0.0E+00
  BSTAR   = 1.0E-06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL RANDOM_SEED() !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NDELTA.GT.NPROC)NDELTA=NPROC-1
  IF(NDELTA.EQ.0)    NDELTA=1
  DELTA = UNIFDISC(1,NDELTA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NCR
     CR(I) = REAL(I)/REAL(NCR)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IDCR = MULTINOM(NCR,PCR)
  CRDREAM = CR(IDCR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DESTART = INTMIN(VETCONT,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(DESTART.GT.0)THEN
     SIG(K) = 0.0E+00
     MATDE  = 0.0E+00
!! carrega matriz MATDE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL LOADTHETASDE(NK,K,DIM,VETCONT,NPROC)
     STDX = STDV(MATDE,DIM,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     XSUM = 0.0
     DO J=1,DELTA
        R    = CHOOSER(NK+1,NPROC)
        XSUM = XSUM + (MATDE(1:DIM,R(1))-MATDE(1:DIM,R(2)))
        write(*,*)'NK',NK,'D',J,'paresR',R
     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! escolha de quais dimensoes serao atualizadas !!!!!!!!!
!! obs: se CRDREAM = 1.0 todas as dimensoes sao atualiz !     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DLINHA = DIM
     ID     = 0
     DO J=1,DIM
        FIRST  = .TRUE.
        TEST   = RANDOM_BINOMIAL2(1,1.0-CRDREAM,FIRST)
        DLINHA = DLINHA-TEST
        IF(TEST.EQ.1)THEN
           THETAN(K,J) = MATDE(J,NK+1)
        ELSE
           ID       = ID+1
           DVET(ID) = J
        END IF
     END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(DLINHA.EQ.0)THEN
        DLINHA  = 1
        DVET(1) = UNIFDISC(1,DIM)
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! proposal !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     GAMMA_DREAM = 2.38/SQRT(2.0*REAL(DLINHA*DELTA))
     GAMMA_DREAM = CR_DREAM
     GAMMA_DREAM = GAMMA_DREAM_METHOD(CR_DREAM,CONTADORC,NFREQ(K),VARPROP(K))
     IF(GAMMA_DREAM.GT.SQRT(0.5)) GAMMA_DREAM = SQRT(0.5)-1.0E-03
!! GERACAO DE UMA PROPOSTA INTEIRAMENTE NOVA !!!!!!!!!!!!
!! INTRODUCAO DE UM INDIVIDUO EXTERNO !(theta novo)!!!!!!
!     IF(MOD(CONTADORC,NFREQ).EQ.0)THEN
!        GAMMA_DREAM = SQRT(0.5)-1.0E-04
!        DO I=1,DIM
!           XSUM(I) = (1.0/GAMMA_DREAM)*RANDOM_NORMAL()
!           DVET(I) = I
!        END DO
!        DLINHA = DIM
!     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE(*,100)DELTA,DLINHA,NCR,GAMMA_DREAM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     DO I=1,DLINHA
        J = DVET(I)
        CALL RANDOM_NUMBER(VUNIF)
        VUNIF       = (B-A)*VUNIF+A
        COEFF       = (1.0-VUNIF)*GAMMA_DREAM
        GAMMA_DRAUX = SQRT(1.0-2.0*COEFF*COEFF)
        DXVET(J)    = COEFF*XSUM(J) + BSTAR*RANDOM_NORMAL()
        THETAN(K,J) = GAMMA_DRAUX*MATDE(J,NK+1) + DXVET(J)
     END DO
     JUMP = COEFF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL GERA_AMTHETA(MKL(K),GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     JVET(K,IDCR) = JVET(K,IDCR) + SUMSQR(DXVET,STDX,DIM,DLINHA)
     NID(K,IDCR)  = NID(K,IDCR) + 1
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT('NUMERO DE PARES USADOS........: ',I4,/,&
           'NUMERO DE DIMENSOES TROCADAS..: ',I4,/,&
           'NUMERO NCR....................: ',I4,/,&
           'LAMBDA METODO DREAM...........: ',F6.3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE DREAM_METHOD_MOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DE_METHOD(K,NK,NPROC,DIM)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE MPI
  USE VARIAVEIS, ONLY:  SIG,MKL,XRMVN,MATDE,XRMVN_MEAN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,CONT,VETCONT
  USE VARIAVEIS, ONLY:  CONTADORC,NFREQ,JUMP,VARPROP
!
  REAL                 :: XEPS,XCS,TOL
  INTEGER              :: K,I,N,DIM,NK,DESTART
  INTEGER              :: IERROR,TAG,NPROC
  REAL,DIMENSION(DIM)  :: XVET
  LOGICAL              :: FIRST
  INTEGER,DIMENSION(2) :: R
  INTEGER, EXTERNAL    :: INTMIN
  REAL, EXTERNAL       :: GAMMA_DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TOL    = 1E-08
  FIRST  = .TRUE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(SIG(K).LE.TOL)THEN
     XEPS = 2.38/SQRT(2.0*REAL(DIM))
  ELSE
     XEPS = SIG(K)
  END IF
  XEPS   = GAMMA_DREAM_METHOD(XEPS,&
       CONTADORC,NFREQ(K),VARPROP(K))
  IERROR = 0
  XCS    = SQRT(1.0E-03)
  XVET   = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DESTART = INTMIN(VETCONT,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(DESTART.GT.0)THEN
     MATDE  = 0.0E+00
     CALL LOADTHETASDE(NK,K,DIM,VETCONT,NPROC)
     R = CHOOSER(NK+1,NPROC)
     DO I=1,DIM
        XVET(I) = RANDOM_NORMAL()*XCS
     END DO
     THETAN(K,1:DIM) = (SQRT(1.0 - 2*XEPS*XEPS))*&
          MATDE(1:DIM,NK+1) + &
          XEPS*(MATDE(1:DIM,R(1))-MATDE(1:DIM,R(2)))+XVET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL GERA_AMTHETA(MKL(K),GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT('LAMBDA METODO DE........: ',F6.3)
!
END SUBROUTINE DE_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION RANDOM_WALK_METHOD(SIGK,MKLK,CT,NFREQ,VP)
!
  USE RANDOM
!
  REAL,    INTENT(IN) :: SIGK
  INTEGER, INTENT(IN) :: CT,MKLK
  REAL                :: AUX,TOL
  INTEGER             :: NFREQ,VP
!
  !  CALL RANDOM_SEED()
  TOL = 1.0E-08
  IF(SIGK.LT.TOL)THEN
     AUX = 2.38/SQRT(REAL(MKLK))
  ELSE
     AUX = SIGK
  END IF
  IF(CT.LE.20) AUX = AUX*2.0
  IF(VP.EQ.1)THEN
     TOL = (5.0D-01)*SIGK
     CALL RANDOM_NUMBER(AUX)
     AUX = (SIGK-TOL)*AUX + TOL
  END IF
  IF(MOD(CT,NFREQ).EQ.0)THEN
     IF(AUX.LT.1.0)THEN
        AUX = 1.0
     ELSE
        AUX = SIGK * 2.0
     END IF
  END IF
!
  WRITE(*,100)SIGK,AUX
  RANDOM_WALK_METHOD = AUX
!
100 FORMAT('RANDOM WALK JUMP (',F5.4,')...: ',F5.4)
!
END FUNCTION RANDOM_WALK_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION GAMMA_DREAM_METHOD(SIGK,CT,NFREQ,VP)
!
  USE RANDOM
!
  REAL,    INTENT(IN) :: SIGK
  INTEGER, INTENT(IN) :: CT
  REAL                :: AUX,TOL
  INTEGER             :: NFREQ,VP
!
  !  CALL RANDOM_SEED()
  AUX = SIGK
  IF(CT.LE.20) AUX = AUX*2.0
!
  IF(VP.EQ.1)THEN
     TOL   = (5.0D-01)*SIGK
     CALL RANDOM_NUMBER(AUX)
     AUX = (SIGK-TOL)*AUX + TOL
  END IF
  IF(MOD(CT,NFREQ).EQ.0)THEN
     IF(AUX.LT.1.0)THEN
        AUX = 1.0
     ELSE
        AUX = 2.0 * SIGK
     END IF
  END IF
!
  WRITE(*,100)SIGK,AUX
  GAMMA_DREAM_METHOD = AUX
!
100 FORMAT('DE JUMP (',F5.4,')...: ',F5.4)
!
END FUNCTION GAMMA_DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AM_METHOD(K,NK,SIGK)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE VARIAVEIS, ONLY:  MKL,COVMATUP,XRMVN_VET,DIMR
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,XRMVN_MEAN,XRMVN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,SIG,CONT,NINICIO
  USE VARIAVEIS, ONLY:  NSTART,NFREQAM,NPCHAIN
!
  REAL                          :: XEPS,XCS,XCSININT,XAUX
  INTEGER                       :: K,I,J,N,DIM,NK
  INTEGER                       :: IERROR
  LOGICAL                       :: FIRST
  REAL,DIMENSION(MKL(K))        :: XVET
  REAL                          :: SIGK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NINICIO= NSTART
  DIM    = MKL(K)
  IERROR = 0D0
  FIRST  = .TRUE.
  XCS    = 2.50E-02
  XCSINIT= 1.00E+00
  XEPS   = ((2.38**2)/(REAL(DIM)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(CONT.GE.NINICIO)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     SIG(K) = 0.0E+00
     CALL LOADTHETAS(K,DIM,NPCHAIN,CONT,NK)
     XRMVN_MEAN = XRMVN_VET(1:DIM,NPCHAIN-1)
     COVMATUP   = COVAR(XRMVN_VET,DIM,DIMR,NPCHAIN)
     COVMATUP   = XEPS*COVMATUP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! VERIFICAR SE "COVMATUP" EH SIMETRICA + DEFINIDA !!!!!!
     CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,COVMATUP,COVMATDOWN,FIRST,XRMVN,IERROR)
     IF(IERROR.EQ.1)THEN
        WRITE(*,*)'MATRIZ NAO EH SIMETRICA POSITIVA DEFINIDA'
        CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,XEPS*COVMATID,COVMATDOWN,&
             FIRST,XRMVN,IERROR)
        THETAN(K,1:DIM) = XRMVN
     ELSE
        WRITE(*,*)'MATRIZ SIMETRICA POSITIVA DEFINIDA'
        CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,XEPS*XCS*COVMATID,COVMATDOWN,&
             FIRST,XVET ,IERROR)
        THETAN(K,1:DIM) = XRMVN + XVET
     END IF
     CALL GERA_AMTHETA(DIM,GERATIPO(K),K,NK)
!     write(*,*)'covid',COVMATID
!     write(*,*)'cov',COVMATUP
  END IF
!  ELSE
!     IF(CONT.GT.0)THEN
!        CALL LOADTHETAS(K,DIM,2,CONT,NK)
!        XRMVN_MEAN = XRMVN_VET(1:DIM,1)
!        XAUX       = XCSINIT*SIGK
!     ELSE
!        XAUX       = XCSINIT*SIGK
!     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,XAUX*COVMATID,COVMATDOWN,&
!          FIRST,XVET ,IERROR)
!     THETAN(K,1:DIM) = XVET
!  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL GERA_AMTHETA(DIM,GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!
END SUBROUTINE AM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COPY_DIRUP(NPR,NS)
  !
  USE VARIAVEIS, ONLY: NUPSC,NPRIORR
  IMPLICIT NONE
!  
  INTEGER            :: NPR,NS,K
  CHARACTER (LEN=128):: COMMAND
  CHARACTER (LEN=4)  :: NUMB
!
  write(*,100)NPR
  WRITE(NUMB,'(I4.3)')NPR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  IF(NS.EQ.2)THEN
     DO K=1,NPRIORR
        IF(NUPSC(K).EQ.0)THEN
           COMMAND = ('cp -r ./upscaling/in ')
           COMMAND = TRIM(ADJUSTL(COMMAND))//(' ./upscaling/in')//TRIM(ADJUSTL(NUMB))
        END IF
        IF(NUPSC(K).EQ.1)THEN
           COMMAND = ('cp -r ./upscaling1/in ')
           COMMAND = TRIM(ADJUSTL(COMMAND))//(' ./upscaling1/in')//TRIM(ADJUSTL(NUMB))
        END IF
        IF(NUPSC(K).EQ.2)THEN
           COMMAND = ('cp -r ./upscaling2/in ')
           COMMAND = TRIM(ADJUSTL(COMMAND))//(' ./upscaling2/in')//TRIM(ADJUSTL(NUMB))
        END IF
!
        WRITE(*,*)TRIM(COMMAND)  
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END DO
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
  100 FORMAT('EXPEREMENT COPY: ',I3)
!
END SUBROUTINE COPY_DIRUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE COPY_DIR(NPR,NUMSIMUL,NS)
  !
  USE VARIAVEIS, ONLY: GERATIPO,NPRIORR
  IMPLICIT NONE
!  
  INTEGER            :: NPR,NUMSIMUL,NS,J,TIPO
  CHARACTER (LEN=128):: COMMAND
  CHARACTER (LEN=4)  :: NUMB
!
  write(*,100)NPR
!
  IF(NUMSIMUL.EQ.0)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../blackbox/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../blackbox/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./blackbox/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./blackbox/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.2)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../simuladorRigido/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../simuladorRigido/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./simuladorRigido/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./simuladorRigido/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.3)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../simul_comp/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../simul_comp/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./simul_comp/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./simul_comp/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.5)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../SIMULADOR_VISCOELASTICO/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../SIMULADOR_VISCOELASTICO/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./SIMULADOR_VISCOELASTICO/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./SIMULADOR_VISCOELASTICO/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.6)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../twophaseflow/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../twophaseflow/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./twophaseflow/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./twophaseflow/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.7)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../twophaseflow/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../twophaseflow/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./twophaseflow/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./twophaseflow/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     END IF
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
  100 FORMAT('EXPEREMENT COPY: ',I3)
!
END SUBROUTINE COPY_DIR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CORRIGE_SIMULC(NRANK,NPROCS,K)
  USE VARIAVEIS
  USE MPI
!
  IMPLICIT NONE  
  INTEGER :: NRANK,NCONT,NPROCS,K,NERROREAD,NK,J
  INTEGER :: IERROR
!
  NCONT = 0
  DO WHILE(ERROSIMULADOR.EQ.1)
     NCONT = NCONT+1
     WRITE(*,*)'RANK: ',NRANK, 'REPETINDO "COARSE" PELA:',NCONT,' VEZES'
     DO K=1,NPRIORR
        CALL GERAFILEIN(NERROREAD,NRANK,NPROCS,K)
        CALL GERADOR(CONTADORC,NRANK,K)
        IF(NERROREAD.NE.1)THEN
           IF(NRANK.EQ.0)WRITE(*,*)'######### PROBLEMA NA ###########'
           IF(NRANK.EQ.0)WRITE(*,*)'#### CRIACAO DAS ENTRADAS DO ####'
           IF(NRANK.EQ.0)WRITE(*,*)'######  GERADOR DE CAMPOS  ######'
           STOP
        END IF
     END DO
     CALL PRIOR_RATIO()
     DO K=1,NPRIORR
        CALL UPSCALING(K,NRANK)
     END DO 
     CALL SIMULADOR_C(NRANK)
  END DO
!
END SUBROUTINE CORRIGE_SIMULC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CORRIGE_SIMULF(NRANK,NPROCS,K,NST)
  USE VARIAVEIS
  USE MPI
!
  IMPLICIT NONE  
  INTEGER :: NRANK,NCONT,NPROCS
  INTEGER :: K,NERROREAD,NST,NK,J,IERROR
!
  IF(NST.EQ.2)THEN
100  CONTINUE
     CALL CORRIGE_SIMULC(NRANK,NPROCS,K)
     CALL READ_AMOSTRA_C(NRANK)
     CALL LIKELIHOOD_C()
     CALL MCMC_C(CONTADORC)
     write(*,*)'()()()()()()()()()()()()()()()()()()()()()()()()()()()',accept,errosimulador
     IF(ACCEPT.GT.0)THEN
        write(*,*)'(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)',accept,errosimulador        
        CALL SIMULADOR_F(NRANK)
        write(*,*)'(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)',accept,errosimulador        

        IF(ERROSIMULADOR.EQ.1) GOTO 100
     ELSE
        write(*,*)'(.)-(.)-(.)-(.)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)-(-)',accept,errosimulador        
        GOTO 100
     END IF
  END IF
  IF(NST.EQ.1)THEN
     NCONT = 0
     DO WHILE(ERROSIMULADOR.EQ.1)
        NCONT = NCONT+1
        WRITE(*,*)'RANK: ',NRANK, 'REPETINDO "FINE" PELA:',NCONT,' VEZES'
        DO K=1,NPRIORR
           CALL GERAFILEIN(NERROREAD,NK,NPROCS,K)
           CALL GERADOR(CONTADORC,NK,K)
        END DO
        CALL PRIOR_RATIO()
        CALL SIMULADOR_F(NRANK)
     END DO
  END IF
!
END SUBROUTINE CORRIGE_SIMULF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SAVE MPI FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SAVEMPICONT(NCONT,NPROC,NK,NSINCRO)
  !
  USE MPI
  USE VARIAVEIS, ONLY : VETCONT
  IMPLICIT NONE
  INTEGER            :: IERR,FID,NPROC,NK,I,ISTAT
  INTEGER            :: NCONT
  CHARACTER(LEN=128) :: FILENAME
  CHARACTER(LEN=5)   :: CHARAC
  CHARACTER(LEN=1)   :: NSINCRO
!
  IF(NSINCRO.EQ.'N')THEN
     FID = 3000+NK
     WRITE(CHARAC,'(I4.3)')NK
     CHARAC=TRIM(ADJUSTL(CHARAC))
     FILENAME = TRIM(ADJUSTL('out/vetcont'))//TRIM(ADJUSTL(CHARAC))//&
          TRIM(ADJUSTL('.dat'))
!
     OPEN(UNIT=FID,FILE=FILENAME,STATUS='REPLACE',&
          ACTION='WRITE',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILENAME
        STOP
     END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     WRITE(FID,100) NCONT
!
     CLOSE(UNIT=FID)
  END IF
  IF(NSINCRO.EQ.'Y')THEN
!     CALL MPI_ALLGATHER(NCONT,1,MPI_INTEGER,VETCONT,1,&
!          MPI_INTEGER,MPI_COMM_WORLD,IERR)
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT(I7)
111 FORMAT(I5)
!
END SUBROUTINE SAVEMPICONT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SAVE MPI FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READMPICONT(NCONT,VET,NPROC,NK,NSINCRO)
  !
  use mpi
  USE VARIAVEIS, ONLY      : VETCONT
  IMPLICIT NONE
  INTEGER                  :: IERR,FID,NPROC,NK
  INTEGER                  :: NCONT,I,ISTAT,K
  INTEGER,DIMENSION(NPROC) :: VET
  CHARACTER(LEN=128)       :: FILENAME
  CHARACTER(LEN=5)         :: CHARAC
  INTEGER                  :: I_NUMBER,IERROR
  LOGICAL                  :: I_EXIST,I_OPEN
  CHARACTER(LEN=1)         :: NSINCRO
!
  IF(NSINCRO.EQ.'N')THEN
     FID = 4000+NK
     DO K=1,NPROC
        WRITE(CHARAC,'(I4.3)')K-1
        CHARAC=TRIM(ADJUSTL(CHARAC))
        FILENAME = TRIM(ADJUSTL('out/vetcont'))//TRIM(ADJUSTL(CHARAC))//&
             TRIM(ADJUSTL('.dat'))
900     CONTINUE
        INQUIRE(FILE=FILENAME, OPENED=I_OPEN, EXIST=I_EXIST, NUMBER=I_NUMBER)
        IF(I_OPEN)THEN
        !        GOTO 900
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)I_OPEN, I_EXIST,I_NUMBER
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           WRITE(*,*)'#___________________________________#'
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
           CALL MPI_FINALIZE(IERROR)
           stop 5265
        ELSE
           OPEN(UNIT=FID,FILE=FILENAME,STATUS='UNKNOWN',&
                ACTION='READ',IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILENAME
              STOP
           END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           READ(FID,100)NCONT
           VET(K)=NCONT
!
           CLOSE(UNIT=FID)
        END IF
     END DO
  END IF
  IF(NSINCRO.EQ.'Y')THEN
     CALL MPI_ALLGATHER(NCONT,1,MPI_INTEGER,VETCONT,1,&
          MPI_INTEGER,MPI_COMM_WORLD,IERR)
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT(I7)
111 FORMAT(I5)
!
END SUBROUTINE READMPICONT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION SUMSQR(DX,ST,DM,DL)
!
  IMPLICIT NONE
  INTEGER,INTENT(IN):: DM,DL
  INTEGER           :: I
  REAL,DIMENSION(DM):: DX
  REAL              :: SUM,AUX
  REAL,DIMENSION(DM),INTENT(IN)   :: ST
!
  SUM = 0.0
  DO I=1,DM
     AUX = DX(I)/ST(I)
     SUM = SUM + AUX*AUX
  END DO
  SUMSQR = SUM
END FUNCTION SUMSQR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MVNCASE(NK,NPROC,ND,NDR)
  USE VARIAVEIS
  USE STATFUNCTIONS
  IMPLICIT NONE
!
  REAL, DIMENSION(NDR) :: MCOVAR,MCOVARID,MCOVARD
  REAL, DIMENSION(ND)  :: XMEAN,XNEW
  INTEGER              :: K,NERROREAD,NK,NPROC,ND,NDR
  INTEGER              :: J,I,LOC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!! MATRIZ IDENTIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,ND
     DO J=I,ND
        LOC = J*(J-1)/2+I
        IF(I.EQ.J) MCOVARID(LOC) = 1.0E+00
     END DO
  END DO
!!#####################################################!!
!! MATRIZ DE COVARIANCIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MCOVARID = 0.5
  DO I=1,ND
     DO J=I,ND
        LOC = J*(J-1)/2+I
        IF(I.EQ.J) MCOVARID(LOC) = DFLOAT(J)
     END DO
  END DO
  MCOVARID = 0.0
  XMEAN    = 0.0
  XNEW     = 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  NPRIORR = 1
  K       = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LOOP SOBRE AS REALIZACOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO WHILE(CONT.LT.NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SALVA O STATUS ATUAL DA SIMULACAO !!!!!!!!!!!!!!!!!!!!!
     CALL RANDOM_SEED(put=NSEED)
!
     CONTADORC=CONTADORC+1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(CONTADORC.GT.NRTOTAL)THEN
        IF(NK.EQ.0)WRITE(*,*)'#################################'
        IF(NK.EQ.0)WRITE(*,*)'### SAIDA SEM CONVERGENCIA ######'
        IF(NK.EQ.0)WRITE(*,*)'#################################'
        IF(NK.EQ.0)WRITE(*,332)REAL(CONT)/REAL(CONTADORF)*100.d0
        EXIT
     ENDIF
!
     IF(NK.EQ.0)WRITE(*,*)'#########################################################'
     IF(NK.EQ.0)WRITE(*,*)'######################## REALIZACAO:',CONTADORC,' #######'
     IF(NK.EQ.0)WRITE(*,*)'######################## RANK......:',NK    ,' #######'
     IF(NK.EQ.0)WRITE(*,*)'#########################################################'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! GERACAO DO CAMPO ALEATORIO !!!!!!!!!!!!!!!!!!!!!!!!!!!
     GERATIPO = 99
     MKL      = 2
!     CALL GERAFILEIN(NERROREAD,NK,NPROCS,K)
!        CALL GERADOR(CONTADORC,NK,K)
!        IF(NERROREAD.NE.1)THEN
!           IF(NK.EQ.0)WRITE(*,*)'######### PROBLEMA NA ###########'
!           IF(NK.EQ.0)WRITE(*,*)'#### CRIACAO DAS ENTRADAS DO ####'
!           IF(NK.EQ.0)WRITE(*,*)'######  GERADOR DE CAMPOS  ######'
!           STOP
!        END IF
!     END DO
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! TWO STAGE ALGORITHM !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IF(NSTAGE.EQ.2)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  COARSE PROBLEM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPSCALING DO CAMPO GERADO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        DO K=1,NPRIORR
!           CALL UPSCALING(K,NK)
!        END DO
!       CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESOLUCAO DO PROBLEMA DE TRANSPORTE (COARSE SCALE) !!!
!        CALL SIMULADOR_C(NK)
!        IF(NSIMUL.EQ.5.AND.ERROSIMULADOR.EQ.1)THEN
!           WRITE(*,1001)NK
!           GOTO 999
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DA AMOSTRA !!!!!!!!!!!!!!!!!!!!!!!!!
!        CALL READ_AMOSTRA_C(NK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULO DO LIKELIHOOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        CALL LIKELIHOOD_C()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROCESSO DE ACEITACAO-REJEICAO !!!!!!!!!!!!!!!!!!!!!!!
!        CALL MCMC_C(CONTADORC)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ELSE
!        ACCEPT=1D0
!     END IF
!!#####################################################!!
!! END TWO STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! FINE SCALE PROBLEM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IF(ACCEPT.GT.0)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESOLUCAO DO PROBLEMA DE TRANSPORTE (FINE SCALE)!!!!!!
!        CALL SIMULADOR_F(NK)
!        IF(NSIMUL.EQ.5.AND.ERROSIMULADOR.EQ.1)THEN
!           WRITE(*,1000)NK
!           GOTO 999
!        END IF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LEITURA DOS DADOS DA AMOSTRA (FINE SCALE)!!!!!!!!!!!!!
!        CALL READ_AMOSTRA_F(NK)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CALCULO DO LIKELIHOOD (FINE SCALE)!!!!!!!!!!!!!!!!!!!!
!        CALL LIKELIHOOD_F()
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! PROCESSO DE ACEITACAO-REJEICAO !!!!!!!!!!!!!!!!!!!!!!!
!        CALL MCMC_F(CONTADORC,NK)
!!
        CONTADORF=CONTADORF+1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! COPIA DOS CAMPOS ACEITOS E REJEITADOS !!!!!!!!!!!!!!!!
!    CALL COPYFILES(CONTADORC,NAMEPRIOR,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    CALL RANDOM_SEED(get=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IF(ACCEPT.GT.0)THEN
!        CONT=CONT+1
!        IF(NPROPOSAL(1).EQ.4)THEN
!           VETCONT(NK+1)=CONT
!           CALL SAVEMPICONT(CONT,1,NK)
!        END IF
!     ELSE
!        CONTREJ=CONTREJ+1
!     END IF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
332  FORMAT('================================',/,&
         'RANK................. =',I3,/,&
         'TAXA DE ACEITACAO (%) =',F8.2,/,&
         '================================') 
END SUBROUTINE MVNCASE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ENCONTRAR CADEIAS OUTLIERS !!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTLIER_FINDER(LOGERR,NPR,INDX,INDXT,CONTOUT)
!
  IMPLICIT NONE
  INTERFACE
     SUBROUTINE SASM(x,quart)
       REAL*4 :: x(:),quart
     END SUBROUTINE SASM
  END INTERFACE
!
  REAL   ,INTENT(IN),DIMENSION(NPR)  :: LOGERR
  REAL   ,DIMENSION(NPR)             :: VX
  INTEGER,INTENT(OUT),DIMENSION(NPR) :: INDX,INDXT
  INTEGER,INTENT(OUT)                :: CONTOUT
  INTEGER :: NPR,I,NK
  REAL    :: EPS,QUART1,QUART3,IQR,TOLS,TOLI,MEDIAN
!
  VX     = LOGERR
  EPS    = 1.0E-06
  INDX   = 0
  INDXT  = 0
  CONTOUT= 0
!
  CALL SORT(NPR,VX,INDX,EPS)
!
  QUART1 = 0.25D0
  CALL SASM(VX,QUART1)
  QUART3 = 0.75D0
  CALL SASM(VX,QUART3)
  MEDIAN = 0.50D0
  CALL SASM(VX,MEDIAN)
!
  IQR = QUART3-QUART1
  TOLI= QUART1 - 2.0*IQR
  TOLS= QUART3 + 2.0*IQR
  INDX= 0
!
  DO I=1,NPR
!     IF(LOGERR(I).GT.TOLS.OR.LOGERR(I).LT.TOLI)THEN
     IF(LOGERR(I).GT.TOLS)THEN
        CONTOUT = CONTOUT + 1
        INDX(CONTOUT) = I - 1
     END IF
  END DO
!
  DO I=1,CONTOUT
     NK = INDX(I)
     CALL SWITCH(NK,NPR,INDX,INDXT,CONTOUT,I)
  END DO
!
END SUBROUTINE OUTLIER_FINDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SWITCH(NK,NPR,INDX,INDXT,CONTOUT,LOC)
!
  USE VARIAVEIS    , ONLY: GERATIPO,NPRIORR,LOGERRORCF
  USE VARIAVEIS    , ONLY: CONTADORC
  USE STATFUNCTIONS, ONLY: UNIDRND
!
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: COMMAND
  INTEGER            :: OUT_FILE,ISTAT,J
  INTEGER,INTENT(IN) :: LOC
  INTEGER            :: I,FLAG,MK,K
  LOGICAL            :: SINAL
  INTEGER, INTENT(IN):: NK,NPR,CONTOUT
  CHARACTER(LEN=128) :: NOME,FOUT,NAME
  CHARACTER(LEN=4)   :: NUMB,NUMA
  INTEGER, INTENT(IN) ,DIMENSION(NPR) :: INDX
  INTEGER, INTENT(OUT),DIMENSION(NPR) :: INDXT
!
  SINAL = .TRUE.
  DO WHILE(SINAL)
     MK = UNIDRND(NPR-1)
     SINAL = .FALSE.
     DO J=1,CONTOUT
        IF(MK.EQ.INDX(J)) SINAL=.TRUE.
     END DO
  END DO
  INDXT(LOC) = MK
!
  WRITE(NUMB,'(I4.3)')NK
  WRITE(NUMA,'(I4.3)')MK
!
  DO K=1,NPRIORR
     IF(GERATIPO(K).EQ.0)THEN
        NOME=ADJUSTL(TRIM('../gera_LABTRANGEORW3/out/xpsi.dat'))
     END IF
     IF(GERATIPO(K).EQ.1)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.3)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.13)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.14)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D_2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_KL3D_2/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.31)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.34)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.32)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.33)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.6)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/MVN/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.7)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN1/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/MVN1/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     IF(GERATIPO(K).EQ.8)THEN
        NOME=ADJUSTL(TRIM('../gera_KL/MVN2/out/theta'))//&
             TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
        NAME=ADJUSTL(TRIM('../gera_KL/MVN2/out/theta'))//&
             TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
     END IF
     COMMAND = ('cp -f ')//TRIM(ADJUSTL(NAME))//' '//TRIM(ADJUSTL(NOME))
     CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
     WRITE(*,*)COMMAND
     WRITE(*,100)CONTADORC,NK,MK
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DO K=1,0!NPRIORR
  !    IF(GERATIPO(K).EQ.0)THEN
  !       NOME=ADJUSTL(TRIM('../gera_LABTRANGEORW3/out/xpsi.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.1)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.3)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.31)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW1/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.34)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW2/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.32)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_RW3D/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.33)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/FORTRAN_LBD/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.6)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/MVN/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/MVN/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.7)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/MVN1/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/MVN1/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    IF(GERATIPO(K).EQ.8)THEN
  !       NOME=ADJUSTL(TRIM('../gera_KL/MVN2/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMB))//TRIM(ADJUSTL('.dat'))
  !       NAME=ADJUSTL(TRIM('../gera_KL/MVN2/out/thetanew'))//&
  !            TRIM(ADJUSTL(NUMA))//TRIM(ADJUSTL('.dat'))
  !    END IF
  !    COMMAND = ('cp -f ')//TRIM(ADJUSTL(NAME))//' '//TRIM(ADJUSTL(NOME))
  !    CALL EXECUTE_COMMAND_LINE(COMMAND,WAIT=.TRUE.)
  !    WRITE(*,*)COMMAND
  !    WRITE(*,100)NK,MK
  ! END DO
!
100 FORMAT('CONTADORC NO ROOT: ',I4,'   CHANGE:',I4,'   BY',I4)
END SUBROUTINE SWITCH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SORT(n, ra, ind, eps)
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
! use kinds, ONLY : DP
  implicit none  
  !-input/output variables
  integer, intent(in)   :: n  
  real(4), intent(in)  :: eps
  integer :: ind (n)  
  real(4) :: ra (n)
  !-local variables
  integer :: i, ir, j, l, iind  
  real(4) :: rra  
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  

  END DO sorting    
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
    REAL(4) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt

  !
END SUBROUTINE SORT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SASM(x,quart)
! Calculate Quartiles using SAS Method 5
! This method is the default method of SAS and is based on the empirical distribution function. 
! Based on discussion in this paper http://www.haiweb.org/medicineprices/manual/quartiles_iTSS.pdf
!
implicit none
real*4 :: x(:),quart,a,b,c,tol,diff
integer :: n,ib

tol=1.e-8
n=size(x)

a=n*quart
call getgp(a,b,c)

ib=int(c)
!print *,n,a,b,c,ib

diff=b-0.0d0
if(diff <=tol) then
  quart=(x(ib+1)+x(ib))/2.0d0
else
 quart=x(ib+1)
end if

end subroutine SASM

subroutine getgp(a,b,c)
! Subroutine to that returns the Right hand and Left hand side digits of a decimal number
real*4 :: a,b,c

b=mod(a,1.0d0)
c=a-b

end subroutine getgp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOAD_OUTLIER_INFO(NK)
!
  USE VARIAVEIS, ONLY: CONT,NAME_OUT,CONTADORC,NCHAIN
  USE VARIAVEIS, ONLY: LOGERRORCF,NDTYPE
  IMPLICIT NONE
!
  INTEGER            :: OUT_FILE,ISTAT,NK
  CHARACTER(LEN=128) :: FOUT,FERR
  INTEGER            :: NUMB,I,J,NAUX,SOMA,K
  INTEGER,DIMENSION(CONT)   :: NREP
  REAL   ,DIMENSION(NDTYPE) :: AUX
  CHARACTER(LEN=5)   :: CHARAC
  REAL               :: ER
!
  OUT_FILE = 3800+NK
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)
  FOUT = TRIM('./out/nchain_')//TRIM(NAME_OUT)//TRIM('_RK')//&
       TRIM(ADJUSTL(CHARAC))//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
!
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='OLD',&
       ACTION='READ',IOSTAT=ISTAT)
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FOUT
     STOP
  END IF
!
  SOMA=0.0
  DO I=1,CONT-1
     READ(UNIT=OUT_FILE,FMT=100)NAUX,NREP(I)
     SOMA = SOMA+NREP(I)
  END DO
  NREP(CONT) = NCHAIN
  SOMA = SOMA+NREP(CONT)
!
  CLOSE(OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FERR=TRIM('./error/erros_')//TRIM(NAME_OUT)//TRIM('_RK')//&
       TRIM(ADJUSTL(CHARAC))//TRIM('.dat')
  FERR=ADJUSTL(TRIM(FERR))
!
  OPEN(UNIT=OUT_FILE,FILE=FERR,STATUS='OLD',&
       ACTION='READ',IOSTAT=ISTAT)
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FERR
     STOP
  END IF
!
  K = 0
  DO I=1,CONT
     READ(OUT_FILE,FMT=110)NUMB,ER,(AUX(J),J=1,NDTYPE)
     DO J=1,NREP(I)
        K=K+1
        LOGERRORCF(K) = LOG(ER)
     END DO
  END DO
  CLOSE(OUT_FILE) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
100 FORMAT(I7,I7)
110 FORMAT(I7,1X,E15.5,1X,200E15.5)
112 FORMAT(I5)
END SUBROUTINE LOAD_OUTLIER_INFO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DREAMTYPEINITIALIZATION(STARTS,NK,NSINC)
  USE MPI
  USE VARIAVEIS, ONLY : VETCONT,NPROPOSAL,CONT,SINALDREAM
  IMPLICIT NONE
!
  INTEGER :: NK, IERROR
  LOGICAL :: STARTS
  CHARACTER(1) :: NSINC
!
  IF(NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5.OR.&
       NPROPOSAL(1).EQ.6)THEN
     IF(NPROPOSAL(1).EQ.5)THEN
        SINALDREAM = .TRUE.
        NPROPOSAL(1) = 6
     END IF
     VETCONT(NK+1) = CONT
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
     CALL MPI_ALLGATHER(CONT,1,MPI_INTEGER,VETCONT,1,&
          MPI_INTEGER,MPI_COMM_WORLD,IERROR)
     CALL SAVEMPICONT(CONT,1,NK,NSINC)
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
     IF(STARTS)THEN
        CALL LOAD_OUTLIER_INFO(NK)
     END IF
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
END SUBROUTINE DREAMTYPEINITIALIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DREAMTYPEMETHODS(CONTMIN,LOGERROR_LOCAL,&
     LOGERROR,NK,INDX,NPROCS,PCRLOCAL,NSINCRO)
!
  USE MPI
  USE VARIAVEIS, ONLY: NPROPOSAL,LOGERRORCF,CONTADORC,NRTOTAL
  USE VARIAVEIS, ONLY: JVET,JVETALL,NIDALL,NID,NPRIORR
  USE VARIAVEIS, ONLY: NCR,PCRALL,CONT,NR,VETCONT,SINALDREAM
!
  IMPLICIT NONE
!
  INTEGER :: NK,IERROR,ROOT,NPROCS,J,K,CONTOUT,NLIMIT
  INTEGER :: TAG,STATUS(MPI_STATUS_SIZE),NFREQOUTL
  INTEGER, INTENT(OUT)         :: CONTMIN
  REAL   , EXTERNAL            :: LOGERRORMEAN
  INTEGER, EXTERNAL            :: INTMIN
  INTEGER, DIMENSION(NPROCS)   :: INDX,INDXT
  REAL   , DIMENSION(NPROCS)   :: LOGERROR
  REAL   , INTENT(OUT)         :: LOGERROR_LOCAL,PCRLOCAL
  CHARACTER(LEN=1), INTENT(IN) :: NSINCRO
!
  ROOT      =   0
  TAG       = 100
  NFREQOUTL =   1
  NLIMIT    =   1
!
  IF(NPROPOSAL(1).EQ.5.OR.NPROPOSAL(1).EQ.6)THEN
     !! OUTLIERS IDENTIFICATION !!!!!!!!!!!!!!!!!!!!!!
     IF(CONTMIN.GE.NLIMIT.AND.MOD(CONTADORC,NFREQOUTL).EQ.0)THEN
        IF(SINALDREAM)THEN
           NPROPOSAL(1) = 5
           SINALDREAM   = .FALSE.
        END IF
        LOGERROR_LOCAL = LOGERRORMEAN(LOGERRORCF,CONTADORC,NRTOTAL,NFREQOUTL)
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLGATHER(LOGERROR_LOCAL,1,MPI_REAL,LOGERROR,1,&
             MPI_REAL,MPI_COMM_WORLD,IERROR)
        IF(NK.EQ.ROOT)THEN
           CALL OUTLIER_FINDER(LOGERROR,NPROCS,INDX,INDXT,CONTOUT)
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_BCAST(CONTOUT,1,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERROR)
        CALL MPI_BCAST(INDX,NPROCS,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERROR)
        CALL MPI_BCAST(INDXT,NPROCS,MPI_INTEGER,ROOT,MPI_COMM_WORLD,IERROR)
        DO J=1,CONTOUT
           IF(NK.EQ.INDXT(J))THEN
              CALL MPI_SEND(LOGERRORCF,NRTOTAL,MPI_REAL,INDX(J),&
                   TAG,MPI_COMM_WORLD,IERROR)
           END IF
           IF(NK.EQ.INDX(J))THEN
              CALL MPI_RECV(LOGERRORCF,NRTOTAL,MPI_REAL,INDXT(J),&
                   TAG,MPI_COMM_WORLD,STATUS,IERROR)
           END IF
        END DO
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(CONTADORC.GT.20.AND.CONTADORC.LT.(NR/10))THEN
        CALL MPI_REDUCE(JVET,JVETALL,NPRIORR*NCR,MPI_REAL,&
             MPI_SUM,ROOT,MPI_COMM_WORLD,IERROR)
        CALL MPI_REDUCE(NID,NIDALL,NPRIORR*NCR,MPI_INTEGER,&
             MPI_SUM,ROOT,MPI_COMM_WORLD,IERROR)
        IF(NK.EQ.ROOT)THEN
           DO K=1,NPRIORR
              PCRLOCAL = 0.0
              DO J=1,NCR
                 PCRALL(K,J) = JVETALL(K,J)/DBLE(NIDALL(K,J))
                 PCRLOCAL = PCRLOCAL + PCRALL(K,J)
              END DO
              PCRALL(K,:) = PCRALL(K,:)/PCRLOCAL
           END DO
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_BCAST(PCRALL,NPRIORR*NCR,MPI_REAL,&
             ROOT,MPI_COMM_WORLD,IERROR)
        DO K=1,NPRIORR
           WRITE(*,*)'PCRALL2=',nk,K,PCRALL(K,:),CONTADORC
        END DO
     END IF
     CALL READMPICONT(CONT,VETCONT,NPROCS,NK,NSINCRO)
     CONTMIN = INTMIN(VETCONT,NPROCS)
  ELSE
     CONTMIN = CONT
  END IF
!
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
END SUBROUTINE DREAMTYPEMETHODS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
