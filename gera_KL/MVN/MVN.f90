!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!                 RANDOM VARIABLE GENERATOR
!          GENERATE MULTIVARIATE NORMAL DISTRIBUTION
! #################  CONDITIONAL VERSION  ###################
!
! NAME: MARCIO RENTES BORGES
! DATA: 05/10/2018
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------ 
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULO COM AS VARIAVEIS GLOBAIS DO PROGRAMA !!!!!!!!!!!!!!
MODULE VARIAVEIS
  CHARACTER(LEN=128):: FILE_IN
  CHARACTER(LEN=128):: FILE_VAL,FILE_VET,FILE_THE,FILE_OUT
  REAL              :: SIG
  INTEGER           :: NFILE,NDIM,NERROR
  INTEGER           :: NPROPOSAL
  INTEGER, ALLOCATABLE,DIMENSION(:) :: NSEED
  REAL,    ALLOCATABLE,DIMENSION(:) :: THETAN,THETA,COVAR,COVARD
END MODULE VARIAVEIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
!
  USE VARIAVEIS
  USE RANDOM
!
  IMPLICIT NONE
! LIST OF LOCAL VARIABLES
  INTEGER                    :: I,J,ERROR_READ,NUM_ELEM,NDIMR,IERROR
  INTEGER, ALLOCATABLE       :: SEED(:)
  DOUBLE PRECISION, EXTERNAL :: POWERLAW,EXPONENTIAL
  DOUBLE PRECISION           :: TOL
  REAL                       :: TIME
  REAL, DIMENSION(2)         :: TARRAY
  LOGICAL                    :: FIRST
!
  ERROR_READ = 0
  TOL = 1.0E-12
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL READ_DATA(ERROR_READ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MEMORY ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NDIMR = NDIM*(NDIM+1)/2
  ALLOCATE(COVAR(NDIMR))
  ALLOCATE(COVARD(NDIMR))
  ALLOCATE(THETA(NDIM))
  ALLOCATE(THETAN(NDIM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FIRST  = .TRUE.
  COVAR = 0.0
  COVARD= 0.0
  THETAN= 0.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'LOADING COVARIANCE'
  CALL LOAD_COVAR()
  CALL LOAD_THETA()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETS THE SEEDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_SEED(put=NSEED)      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM FIELDS GENERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL ETIME(TARRAY,TIME)
!
  IF(NPROPOSAL.EQ.2)THEN
     SIG = SIG+TOL
     CALL RANDOM_MVNORM(NDIM,THETA,SIG*COVAR,&
          COVARD,FIRST,THETAN,IERROR)
  END IF
!  
  IF(NPROPOSAL.EQ.4)THEN
     THETAN = THETA
  END IF
!
  CALL WRITE_THETA()
  CALL WRITE_FIELD()
  CALL ETIME(TARRAY,TIME)
  WRITE(*,331)TIME,TARRAY(1),TARRAY(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPRESSAO DA SEMENTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL PRINTSEED()
! MEMORY FREE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DEALLOCATE(COVAR,COVARD,THETA,THETAN)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROR = 1
!
333 FORMAT('TIME TO LOAD EIGENPAIRS (s) =',F8.2,/,&
           'USER TIME               (s) =',F8.2,/,&
           'SYSTEM TIME             (s) =',F8.2,/)
!
331 FORMAT('TIME TO GENERATE        (s) =',F8.2,/,&
           'USER TIME               (s) =',F8.2,/,&
           'SYSTEM TIME             (s) =',F8.2,/)
!
END PROGRAM MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ THE INPUT DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_DATA(NERROREAD)
!
  USE VARIAVEIS
  USE RANDOM
!
  INTEGER :: IN_FILE,ISTAT
  INTEGER, INTENT(OUT) :: NERROREAD
  INTEGER,DIMENSION(2) :: LIXO
!
  IN_FILE = 123
  FILE_IN = './entrada.in'
  OPEN(UNIT=IN_FILE,FILE=FILE_IN,STATUS='UNKNOWN',&
       FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_IN
     STOP
  END IF
!
!        READ(IN_FILE,300)LIXO(1),LIXO(2)
  READ(IN_FILE,100)NDIM
  WRITE(*,200)NDIM
!        READ(IN_FILE,100)LIXO(1)
!
  CALL RANDOM_SEED(size=NS)
  ALLOCATE(NSEED(NS))
  READ(IN_FILE,97)(NSEED(I),I=1,NS)
  WRITE(*,96)NS,(NSEED(I),I=1,NS)
!
  READ(IN_FILE,7000)FILE_OUT
  FILE_OUT = TRIM(ADJUSTL(FILE_OUT))
  WRITE(*,7003)TRIM(ADJUSTL(FILE_OUT))
!
  READ(IN_FILE,7000)FILE_VET
  FILE_VET = TRIM(ADJUSTL(FILE_VET))
  WRITE(*,7001)TRIM(ADJUSTL(FILE_VET))
!
  READ(IN_FILE,7000)FILE_THE
  FILE_THE = TRIM(ADJUSTL(FILE_THE))
  WRITE(*,7002)TRIM(ADJUSTL(FILE_THE))
!
  READ(IN_FILE,12)NPROPOSAL,SIG
  WRITE(*,13)NPROPOSAL,SIG
  IF(NPROPOSAL.EQ.1)WRITE(*,*)'RANDOM WALK'
  IF(NPROPOSAL.EQ.2)WRITE(*,*)'pCN'
!
  NERROREAD = 1d0
  CLOSE(UNIT=IN_FILE)
!
  RETURN
!
300 FORMAT(2I7)
100 FORMAT(I7)
200 FORMAT('DIMENSAO DA MVN: ',I2)
96 FORMAT('NUMERO DE SEMENTES: ',I2,/,'SEMENTES:',40I12)
97 FORMAT(40I12)
12 FORMAT(I7,E10.3)
13 FORMAT('TIPO DE PROPOSAL         =',I7,/,&
        'PARAMETRO DO RANDOM WALK =',E10.3)
4000 FORMAT(3I12)
4001 FORMAT(/,'SEEDS =',I12,I12,I12,I12,I12,I12,I12,I12,/)
7000 FORMAT(A)
7001 FORMAT('NAME OF COVARIANCE MATRIX..= ',A)
7002 FORMAT('NAME OF MEAN VECTOR........= ',A)
7003 FORMAT('NAME OF OUTPUT FILE........= ',A)
!
END SUBROUTINE READ_DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION VMEAN(VEC,NNX,NNY)
!
  DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
  DOUBLE PRECISION             :: SUM
  INTEGER, INTENT(IN)          :: NNX,NNY
  INTEGER                      :: I,J
!
  SUM=0.0
  DO I=0,NNX-1
     DO J=0,NNY-1
        SUM=SUM+VEC(I,J)
     ENDDO
  ENDDO
!
  VMEAN=SUM/DBLE(NNX*NNY)
!
END FUNCTION VMEAN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION VAR(VEC,NNX,NNY)
!
  DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
  DOUBLE PRECISION             :: SUM,AUX
  DOUBLE PRECISION, EXTERNAL   :: VMEAN
  INTEGER, INTENT(IN)          :: NNX,NNY
  INTEGER                      :: I,J
!
  AUX=VMEAN(VEC,NNX,NNY)
  SUM=0.0
  DO I=0,NNX-1
     DO J=0,NNY-1
        SUM=SUM+VEC(I,J)*VEC(I,J)
     ENDDO
  ENDDO
!
  VAR=SUM/DBLE(NNX*NNY)-AUX*AUX
!
END FUNCTION VAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION FUNCTION VARIANCIA(VEC,NNX)
!
  DOUBLE PRECISION, INTENT(IN) :: VEC(NNX)
  DOUBLE PRECISION             :: SUM,ME
  INTEGER, INTENT(IN)          :: NNX
  INTEGER                      :: I
!
  SUM=0.0
  ME =0.0
  DO I=1,NNX
     ME =ME+VEC(I)
  ENDDO
  ME=ME/DBLE(NNX-1)
!
  DO I=1,NNX
     SUM=SUM+(VEC(I)-ME)**2
  END DO
  VARIANCIA=SUM/DBLE(NNX-1)
!
END FUNCTION VARIANCIA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PRINTSEED()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  CHARACTER(LEN=128)             :: NAME
  INTEGER                        :: I,OUT_FILE
  INTEGER, DIMENSION(8)          :: SEED 
!
  OUT_FILE=341
  NAME='./out/seed.dat'
  NAME=ADJUSTL(TRIM(NAME))
  WRITE(*,111)NAME(1:LEN_TRIM(NAME))
!
  OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
  CALL RANDOM_SEED(get=NSEED)
  WRITE(OUT_FILE,103)(SEED(I),I=1,8)
!
  CLOSE(OUT_FILE)
!
103 FORMAT(I10)
111 FORMAT('NAME OF SEED OUTPUT FILE:',A)
!
END SUBROUTINE PRINTSEED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOAD_COVAR()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER          :: ISTAT,IN_FILE
  INTEGER          :: I,J,NDIMR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IN_FILE = 325
  OPEN(UNIT=IN_FILE,FILE=FILE_VET,STATUS='UNKNOWN', &
       FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_THE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NDIMR = NDIM*(NDIM+1)/2
  DO I=1,NDIMR
     READ(IN_FILE,*)COVAR(I)
  ENDDO
  CLOSE(UNIT=IN_FILE)
  return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE LOAD_COVAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOAD_THETA()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER          :: ISTAT,IN_FILE
  INTEGER          :: I,J
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IN_FILE = 326
  OPEN(UNIT=IN_FILE,FILE=FILE_THE,STATUS='UNKNOWN', &
       FORM='FORMATTED',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_THE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NDIM
     READ(IN_FILE,*)THETA(I)
  ENDDO
!
  CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE LOAD_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_THETA()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: ISTAT,OUT_FILE
  INTEGER            :: I,J
  CHARACTER(LEN=128) :: FILETHE
  CHARACTER(LEN=5)   :: FEXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,128-4
     FEXP = FILE_THE(I:I+4)
     IF(FEXP.EQ.'theta')GOTO 123
  END DO
123 CONTINUE
  FILETHE = TRIM(ADJUSTL(FILE_THE(1:I+4)))//('new')//&
       TRIM(ADJUSTL(FILE_THE(I+5:128)))
  OUT_FILE = 326
  FILETHE=TRIM(ADJUSTL(FILETHE))
  OPEN(UNIT=OUT_FILE,FILE=FILETHE,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NDIM
     WRITE(OUT_FILE,111)THETAN(I)
  ENDDO
!
  CLOSE(UNIT=OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE WRITE_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_FIELD()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER            :: ISTAT,OUT_FILE
  INTEGER            :: I,J
  CHARACTER(LEN=128) :: FILETHE
  CHARACTER(LEN=3)   :: FEXP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        DO I=1,128-2
!           FEXP = FILE_OUT(I:I+2)
!           IF(FEXP.EQ.'exp')GOTO 123
!        END DO
!        123 CONTINUE
!        FILETHE = TRIM(ADJUSTL(FILE_OUT(1:I+2)))//&
!             TRIM(ADJUSTL(FILE_THE(I+5:128)))
  FILETHE = TRIM(ADJUSTL(FILE_OUT))//'0.dat'
  WRITE(*,*)'FILE=',FILETHE
  OUT_FILE = 327
  FILETHE=TRIM(ADJUSTL(FILETHE))
  OPEN(UNIT=OUT_FILE,FILE=FILETHE,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NDIM
     WRITE(OUT_FILE,111)THETAN(I)
  ENDDO
!
  CLOSE(UNIT=OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111 FORMAT(e15.8)
!
END SUBROUTINE WRITE_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
