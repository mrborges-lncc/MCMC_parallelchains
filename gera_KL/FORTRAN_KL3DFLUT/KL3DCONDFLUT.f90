!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!                 3D RANDOM FIELDS GENERATOR
!                USING KARHUNEN-LOEVE METHOD
! #################  CONDITIONAL VERSION  ################### 
! CONDITIONING BASED ON OSSIANDER et al. (2014) -
! Conditional Stochastic Simulations of Flow and Transport
! with Karhunen-Loève Expansions, Stochastic Collocation,
! and Sequential Gaussian Simulation, Journal of 
! Applied Mathematics Volume 2014, Article ID 652594, 21 pages
! http://dx.doi.org/10.1155/2014/652594
! AUTHOR: MARCIO RENTES BORGES
! DATE: 29/04/2021
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!     Metodo baseado na expansao de Karhunen-Loeve         
!     
!     ESTA VERSAO PERMITE CONDICIONAR VALORES EM ALGUNS
!                   PONTOS DO DOMINIO
!                                               
!    
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MODULO COM AS VARIAVEIS GLOBAIS DO PROGRAMA !!!!!!!!!!!!!!
MODULE VARIAVEIS
  CHARACTER(LEN=128):: FILE_IN,FILE_OUT,FILE_COND
  CHARACTER(LEN=128):: FILE_VAL,FILE_VET,FILE_THE
  CHARACTER(LEN=128):: FILE_MU,FILE_MM,FILE_FIELD
  REAL(4)  :: FDX,FDY,FDZ,BASE,SIG
  REAL(4)  :: FBETA,FVAR,FMEAN,DIMX,DIMY,DIMZ
  REAL(4)  :: VARIAN,XMEAN,XVAR,FATCOND
  INTEGER  :: NX,NY,NZ,NFLAG,NERROR,INIF,NFILES
  INTEGER  :: NFILE,NDIM,NS
  LOGICAL  :: TFCOND
  INTEGER  :: NCOND,MKL,NPROPOSAL
  INTEGER, ALLOCATABLE,DIMENSION(:)    :: NSEED,PCOND
  REAL(4), ALLOCATABLE,DIMENSION(:,:)  :: AVET,VET,MM
  REAL(4), ALLOCATABLE,DIMENSION(:)    :: X,Y,Z,AVAL,THETA,MU
  REAL(4), ALLOCATABLE,DIMENSION(:)    :: FIELDM
  REAL(4)  :: BETA,RHO
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
  INTEGER              :: I,J,K,ERROR_READ,NUM_ELEM
  INTEGER, ALLOCATABLE :: SEED(:)
  REAL(4), EXTERNAL    :: POWERLAW,EXPONENTIAL
  REAL                 :: TIME,TSTART,TFINAL
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CPU_TIME(TSTART)
  ERROR_READ = 0
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL READ_DATA(ERROR_READ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MEMORY ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NUM_ELEM = NX*NY*NZ
  ALLOCATE(AVAL(MKL))
  ALLOCATE(AVET(NUM_ELEM,MKL))
  ALLOCATE(MU(MKL))
  ALLOCATE(MM(MKL,MKL))
  ALLOCATE(PCOND(MKL))
  ALLOCATE(THETA(MKL))
  ALLOCATE(FIELDM(NUM_ELEM))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLEAR VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CLEAR(THETA,MKL-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETS THE SEEDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_SEED(put=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM FIELDS GENERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  CALL LOAD_FIELDM(NUM_ELEM)
!
  CALL LOAD_THETA()
!
  CALL LOAD_MM()
!
  CALL LOAD_AUTOV()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM CONDITIONING AND GENERATION !!!!!!!!!!!!!!!!!!!!
!
  CALL KLCOND3D()
!
  CALL CPU_TIME(TFINAL)
  TIME = TFINAL - TSTART
  WRITE(*,331)TIME
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPRESSAO DA SEMENTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL PRINTSEED()
! MEMORY FREE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DEALLOCATE(AVAL,AVET,MU,MM,FIELDM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROR = 1
!
331 FORMAT('TOTAL TIME TO GENERATE  (seconds) =',F8.2)
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
        READ(IN_FILE,8000)INIF,NFILES
        WRITE(*,8001)INIF,NFILES
        NFILE=NFILES-INIF+1
        WRITE(*,8002)NFILE
!
        READ(IN_FILE,1000)DIMX,DIMY,DIMZ
        WRITE(*,1001)DIMX,DIMY,DIMZ
!
        READ(IN_FILE,2000)NX,NY,NZ
        WRITE(*,2001)NX,NY,NZ
!
        FDX = DIMX/NX
        FDY = DIMY/NY
        FDZ = DIMZ/NZ
        WRITE(*,5001)FDX,FDY,FDZ
!
        READ(IN_FILE,3000)NFLAG
!
        IF(NFLAG==2)THEN
           WRITE(*,6002)NFLAG
        END IF
        IF(NFLAG==1)THEN
           WRITE(*,6001)NFLAG
        END IF
!
        READ(IN_FILE,100)BETA,RHO
        WRITE(*,101)BETA,RHO
!
        FBETA  = 0.5
        VARIAN = 1.0
        IF(NFLAG==1)THEN
           WRITE(*,9001)FBETA,VARIAN
        END IF
        IF(NFLAG==2)THEN
           WRITE(*,9002)FBETA,VARIAN
        END IF
!
        CALL RANDOM_SEED(size=NS)
        ALLOCATE(NSEED(NS))
        READ(IN_FILE,97)(NSEED(I),I=1,NS)
        WRITE(*,96)NS,(NSEED(I),I=1,NS)
!
        READ(IN_FILE,3000)MKL
        WRITE(*,2005)MKL
!
        READ(IN_FILE,7000)FILE_OUT
        WRITE(*,7001)FILE_OUT
!
        READ(IN_FILE,7000)FILE_FIELD
        WRITE(*,7007)FILE_FIELD
!
        READ(IN_FILE,7000)FILE_VAL
        WRITE(*,7002)FILE_VAL
!
        READ(IN_FILE,7000)FILE_VET
        WRITE(*,7003)FILE_VET
!
        READ(IN_FILE,7000)FILE_MU
        WRITE(*,7005)FILE_MU
!
        READ(IN_FILE,7000)FILE_MM
        WRITE(*,7006)FILE_MM
!
        READ(IN_FILE,7000)FILE_THE
        WRITE(*,7004)FILE_THE
!
        READ(IN_FILE,12)NPROPOSAL,SIG
        WRITE(*,13)NPROPOSAL,SIG
        IF(NPROPOSAL.EQ.0)WRITE(*,*)'RANDOM'
        IF(NPROPOSAL.EQ.1)WRITE(*,*)'RANDOM WALK'
        IF(NPROPOSAL.EQ.2)WRITE(*,*)'pCN'
        IF(NPROPOSAL.EQ.4)WRITE(*,*)'DE or DREAM'
!
        READ(IN_FILE,3000)NCOND
        WRITE(*,3002)NCOND
!
        ALLOCATE(VET(0:3,0:NCOND-1))
        DO I=0,NCOND-1
           READ(IN_FILE,1003)VET(0,I),VET(1,I),VET(2,I),VET(3,I)
           WRITE(*,1002)I+1,VET(0,I),VET(1,I),VET(2,I),VET(3,I)
           DO J=0,2
              VET(J,I)=VET(J,I)+1E-6
           ENDDO
        ENDDO
!
        IF(NCOND.EQ.0)THEN
           NCOND = 1
        END IF
!
        NERROREAD = 1d0
        CLOSE(UNIT=IN_FILE)
!
        RETURN

100     FORMAT(E12.4,E12.4)
101     FORMAT('BETA:',E12.4,/,'RHO:',E12.4)
96      FORMAT('NUMERO DE SEMENTES:',I5,/,'SEMENTES:',40I12)
97      FORMAT(40I12)
12      FORMAT(I7,E10.3)
13      FORMAT('TIPO DE PROPOSAL         =',I7,/,&
               'PARAMETRO DO RANDOM WALK =',E10.3)
1000    FORMAT(3F12.5)
1003    FORMAT(4F12.5)
1001    FORMAT('DIMENSION X =',F12.5,/,&
               'DIMENSION Y =',F12.5,/,&
               'DIMENSION Z =',F12.5,/)
1002    FORMAT('POINT    =',I6,/,&
               'COORD. X =',F12.5,/,&
               'COORD. Y =',F12.5,/,&
               'COORD. Z =',F12.5,/,&
               'VALUE    =',F12.5,/)
2000    FORMAT(3I7)
2001    FORMAT('N. OF ELEMENTS IN X =',I7,/,&
               'N. OF ELEMENTS IN Y =',I7,/,&
               'N. OF ELEMENTS IN Z =',I7,/)
2005    FORMAT('M (terms) =',I7,/)
3000    FORMAT(I7)
3001    FORMAT('FLAG  =',I7)
3002    FORMAT('NCOND =',I7)
4000    FORMAT(3I12)
4001    FORMAT(/,'SEEDS =',I12,I12,I12,I12,I12,I12,I12,I12,/)
5001    FORMAT('DELTA X =',F12.4,/,&
               'DELTA Y =',F12.4,/,&
               'DELTA Z =',F12.4,/)
6001    FORMAT(I7,': EXPONENTIAL')
6002    FORMAT(I7,': FRACTAL')
7000    FORMAT(A)
7001    FORMAT('NAME OF OUTPUT FILE = ',A)
7002    FORMAT('NAME OF INPUT AUTOVALUES FILE  = ',A)
7003    FORMAT('NAME OF INPUT AUTOVECTORS FILE = ',A)
7004    FORMAT('NAME OF INPUT THETA (NORMAL RV)= ',A)
7005    FORMAT('NAME OF INPUT MU TILDE FILE    = ',A)
7006    FORMAT('NAME OF INPUT PROJECT MATRIX   = ',A)
7007    FORMAT('NAME OF MEAN FIELD             = ',A)
8000    FORMAT(2I7)
8001    FORMAT('N. OF INTIAL FILE =',I7,/,&
               'N. OF FINAL FILE  =',I7)
8002    FORMAT('NUMBER OF FILES   =',I7,/)
9000    FORMAT(2F12.4)
9002    FORMAT('HURST COEF.=',F12.4,/,&
               'VARIANCE   =',F12.4)
9001    FORMAT('CORRELATION LENGTH =',F12.4,/,&
               'VARIANCE   =',F12.4) 
!
      END SUBROUTINE READ_DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR2D(VEC,NNX,NNY)
        INTEGER, INTENT(IN):: NNX,NNY
        INTEGER            :: I,J
        REAL(4)            :: VEC(0:NNX-1,0:NNY-1)
!
        !$OMP PARALLEL PRIVATE(I,J)
!$OMP DO
        DO I=0,NNX-1
           DO J=0,NNY-1
             VEC(I,J)=0.d0
          ENDDO
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR3D(VEC,NNX,NNY,NNZ)
        INTEGER, INTENT(IN):: NNX,NNY,NNZ
        INTEGER            :: I,J,K
        REAL(4)            :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
!
!$OMP PARALLEL PRIVATE(I,J)
!$OMP DO
        DO I=0,NNX-1
           DO J=0,NNY-1
              DO K=0,NNZ-1
                 VEC(I,J,K)=0.d0
              END DO
           END DO
        END DO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE CLEAR(VEC,NN)
        INTEGER, INTENT(IN):: NN
        INTEGER            :: I,J
        REAL(4)            :: VEC(0:NN)
!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
        DO I=0,NN
          VEC(I)=0.d0
        ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
      END SUBROUTINE CLEAR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINTASC(VEC,VX,VY,NNX,NNY,NAME)
        INTEGER,   INTENT(IN) :: NNX,NNY
        REAL(4),   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE
!
        BAND=192837465
        OUT_FILE=322
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
          DO J=0,NNY-1
             WRITE(OUT_FILE,100)J+1
             DO I=0,NNX-1
                WRITE(OUT_FILE,101)VX(I),VY(J),VEC(I,NNY-(1+J))
             ENDDO
             WRITE(OUT_FILE,102)BAND
          ENDDO
        CLOSE(OUT_FILE)
!
 100  FORMAT(I7)
 101  FORMAT(F12.4,2X,F12.4,2X,F12.4)
 102  FORMAT(I9)
      END SUBROUTINE PRINTASC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINTBIN(VEC,VX,VY,NNX,NNY,NAME)
        INTEGER, INTENT(IN) :: NNX,NNY
        REAL(4), INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE
!
        BAND=192837465
        OUT_FILE=21
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')!,FORM='UNFORMATTED')
        WRITE(OUT_FILE,*)VEC
        CLOSE(OUT_FILE)
!
 100  FORMAT(I7)
 101  FORMAT(F12.4,2X,F12.4,2X,F12.4)
 102  FORMAT(I9)
      END SUBROUTINE PRINTBIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MENOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION MENOR(A,B)
        INTEGER :: A,B
!
        IF(A<B)THEN
           MENOR = A
        ELSE
           MENOR = B
        END IF
!
      END FUNCTION MENOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MAIOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION MAIOR(A,B)
        INTEGER :: A,B
!
        IF(A>B)THEN
           MAIOR = A
        ELSE
           MAIOR = B
        END IF
        !
      END FUNCTION MAIOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION NARRED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER FUNCTION NARRED(LD,C)
!
        REAL(4), INTENT(IN) :: C,LD
!        REAL(4)             :: AUX,AUX2
!        INTEGER                      :: N
!!
!        N=INT(LD/C)
!        AUX=(N*C)
!        AUX2=AUX+C*0.5
!        write(*,*)LD,REAL(NINT(LD))
!        IF(LD.LT.AUX2)THEN
        !           NARRED=AUX
!        ELSE
!           NARRED=AUX+C
!        END IF
        NARRED = (NINT(LD))
        !
      END FUNCTION NARRED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(4) FUNCTION VMEAN(VEC,NNX,NNY)
        !
        REAL(4), INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
        REAL(4)             :: SUM
        INTEGER, INTENT(IN) :: NNX,NNY
        INTEGER             :: I,J
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(4) FUNCTION VMEAN3D(VEC,NNX,NNY,NNZ)
!
        REAL(4), INTENT(IN) :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
        REAL(4)             :: SUM
        INTEGER, INTENT(IN) :: NNX,NNY,NNZ
        INTEGER             :: I,J,K
!
        SUM=0.0
        DO I=0,NNX-1
           DO J=0,NNY-1
              DO K=0,NNZ-1
                 SUM=SUM+VEC(I,J,K)
              ENDDO
           ENDDO
        END DO
!
        VMEAN3D=SUM/DBLE(NNX*NNY*NNZ)
!
      END FUNCTION VMEAN3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(4) FUNCTION VAR(VEC,NNX,NNY)
!
        REAL(4), INTENT(IN) :: VEC(0:NNX-1,0:NNY-1)
        REAL(4)             :: SUM,AUX
        REAL(4), EXTERNAL   :: VMEAN
        INTEGER, INTENT(IN) :: NNX,NNY
        INTEGER             :: I,J
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
      REAL(4) FUNCTION VAR3D(VEC,NNX,NNY,NNZ)
!
        REAL(4), INTENT(IN) :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
        REAL(4)             :: SUM,AUX
        REAL(4), EXTERNAL   :: VMEAN3D
        INTEGER, INTENT(IN) :: NNX,NNY,NNZ
        INTEGER             :: I,J,K
!
        AUX=VMEAN3D(VEC,NNX,NNY,NNZ)
        SUM=0.0
        DO K=0,NNZ-1
           DO I=0,NNX-1
              DO J=0,NNY-1
                 SUM=SUM+VEC(I,J,K)*VEC(I,J,K)
              ENDDO
           ENDDO
        END DO
!
        VAR3D=SUM/DBLE(NNX*NNY*NNZ)-AUX*AUX
        !
      END FUNCTION VAR3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION VARIANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(4) FUNCTION VARIANCIA(VEC,NNX)
        !
        REAL(4), INTENT(IN) :: VEC(NNX)
        REAL(4)             :: SUM,ME
        INTEGER, INTENT(IN) :: NNX
        INTEGER             :: I
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
      SUBROUTINE LOAD_FIELDM(N)
        USE VARIAVEIS, ONLY: FILE_FIELD,FIELDM,BETA,RHO
        IMPLICIT NONE
        !
        INTEGER :: N,FID,ISTAT,I
        CHARACTER(LEN=128) :: FILEF
        LOGICAL :: EFILE
        !
        FILEF = TRIM(ADJUSTL(FILE_FIELD))
        INQUIRE(FILE=FILEF,EXIST=EFILE)
        FID = 434
        OPEN(UNIT=FID,FILE=FILEF,STATUS='UNKNOWN', &
             ACTION='READ',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEF
           STOP 135
        END IF
        DO I=1,N
           READ(FID,*)FIELDM(I)
        ENDDO
        CLOSE(FID)    
      END SUBROUTINE LOAD_FIELDM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD_AUTOV()
!
        USE VARIAVEIS, ONLY: AVET,NX,NY,NZ,MKL,FILE_VET
        IMPLICIT NONE
!
        INTEGER          :: NELEM
        INTEGER          :: ISTAT,IN_FILEA,IN_FILEV
        INTEGER          :: I,J,K
        INTEGER(8)       :: NREC=0 ! NREC < 2,147,483,647
        REAL             :: TIME,T_START,T_FINAL
        CHARACTER(LEN=3) :: TFILE
        LOGICAL          :: EFILE
        REAL(4),ALLOCATABLE,DIMENSION(:) :: V
!
        NELEM = NX * NY * NZ
        ALLOCATE(V(NELEM * MKL))
        INQUIRE(IOLENGTH=NREC)V
!        NREC = NELEM * MKL * 4
!!!!! FILE VERIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TFILE = 'ERR'
        DO I=1,LEN(FILE_VET)
           IF(FILE_VET(I:I+3).EQ.'.dat')THEN
              J = I
              TFILE = 'ASC'
           END IF
           IF(FILE_VET(I:I+3).EQ.'.bin')THEN
              J = I
              TFILE = 'BIN'
           END IF
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL CPU_TIME(T_START)
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'ERR')THEN
           WRITE(*,*)'ERROR ON READING INPUT FILE: ',FILE_VET
           STOP 100
        END IF
        IF(TFILE.EQ.'ASC')THEN
           IN_FILEV = 132
           FILE_VET = TRIM(ADJUSTL(FILE_VET))
           INQUIRE(FILE=FILE_VET,EXIST=EFILE)
           OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='UNKNOWN', &
                ACTION='READ',IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
              STOP 101
           END IF
           DO I=1,NELEM
              READ(IN_FILEV,*)AVET(I,1:MKL)
           ENDDO
           CLOSE(IN_FILEV)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'BIN')THEN
           IN_FILEV = 326
           OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='OLD', &
                ACCESS='DIRECT',RECL=NREC,IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
              STOP 103
           END IF
           READ(IN_FILEV,REC=1)V
           K = 0
           DO J=1,MKL
              DO I=1,NELEM
                 K = k + 1
                 AVET(I,J) = V(K)
              END DO
           END DO
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL CPU_TIME(T_FINAL)
        TIME = T_FINAL-T_START
        WRITE(*,334)TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
334     FORMAT('TIME TO LOAD AUT (seconds)     =',F8.2)
111     FORMAT(e25.8)
112     FORMAT(40000e15.7)
!
      END SUBROUTINE LOAD_AUTOV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD_MM()
!
        USE VARIAVEIS, ONLY: MM,MU,MKL,FILE_MM,FILE_MU
        IMPLICIT NONE
!
        INTEGER          :: NELEM
        INTEGER          :: ISTAT,IN_FILEA,IN_FILEV
        INTEGER          :: I,J,K
        INTEGER(8)       :: NREC=0
        REAL             :: TIME,T_START,T_FINAL
        CHARACTER(LEN=3) :: TFILE
        LOGICAL          :: EFILE
        REAL(4),ALLOCATABLE,DIMENSION(:) :: V,VU
!
        NELEM = MKL
        ALLOCATE(V(NELEM * MKL))
        INQUIRE(IOLENGTH=NREC)V
!!!!! FILE VERIFICATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        TFILE = 'ERR'
        DO I=1,LEN(FILE_MM)
           IF(FILE_MM(I:I+3).EQ.'.dat')THEN
              J = I
              TFILE = 'ASC'
           END IF
           IF(FILE_MM(I:I+3).EQ.'.bin')THEN
              J = I
              TFILE = 'BIN'
           END IF
        END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL CPU_TIME(T_START)
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'ERR')THEN
           WRITE(*,*)'ERROR ON READING INPUT FILE: ',FILE_MM
           STOP 301
        END IF
        IF(TFILE.EQ.'ASC')THEN
           IN_FILEV = 137
           FILE_MM = TRIM(ADJUSTL(FILE_MM))
           INQUIRE(FILE=FILE_MM,EXIST=EFILE)
           OPEN(UNIT=IN_FILEV,FILE=FILE_MM,STATUS='UNKNOWN', &
                ACTION='READ',IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_MM
              STOP
           END IF
           DO I=1,NELEM
              READ(IN_FILEV,*)MM(I,1:MKL)
           ENDDO
           CLOSE(IN_FILEV)
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'BIN')THEN
           IN_FILEV = 26
           OPEN(UNIT=IN_FILEV,FILE=FILE_MM,STATUS='OLD', &
                ACCESS='DIRECT',RECL=NREC,IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_MM
              STOP 1544
           END IF
           READ(IN_FILEV,REC=1)V
           K = 0
           DO J=1,MKL
              DO I=1,NELEM
                 K = k + 1
                 MM(I,J) = V(K)
              END DO
           END DO
           CLOSE(IN_FILEV)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IN_FILEV = 637
        FILE_MU = TRIM(ADJUSTL(FILE_MU))
        INQUIRE(FILE=FILE_MU,EXIST=EFILE)
        OPEN(UNIT=IN_FILEV,FILE=FILE_MU,STATUS='UNKNOWN', &
             ACTION='READ',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_MU
           STOP 1978
        END IF
        DO I=1,MKL
           READ(IN_FILEV,*)MU(I)
        ENDDO
        CLOSE(IN_FILEV)        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL CPU_TIME(T_FINAL)
        TIME = T_FINAL-T_START
        WRITE(*,334)TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
334     FORMAT('TIME TO LOAD MMatrix (seconds) =',F8.2)
111     FORMAT(e25.8)
112     FORMAT(40000e15.7)
!
      END SUBROUTINE LOAD_MM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KLCOND3D()
!
        USE VARIAVEIS
        USE RANDOM
        IMPLICIT NONE
!
        REAL(4) :: AUX1,VARR,AUXSQ
        REAL(4) :: AUXM,AUXV
        INTEGER :: I,J,M,MX,MY,MZ
        REAL(4) :: POSIX,POSIY,POSIZ,EPSILON
        REAL(4) :: AUX, SIG2,MEANK,VARK
        REAL(4),ALLOCATABLE,DIMENSION(:) :: KPERM,XI
        INTEGER :: MAX,MIN,MIN2
        CHARACTER(LEN=256) :: NAME
        CHARACTER(LEN=4)   :: EXT
        CHARACTER(LEN=5)   :: C
        REAL               :: TIME,T_START,T_FINAL
        REAL(4), EXTERNAL  :: MAIOR,MENOR
        REAL(4), EXTERNAL  :: VMEAN,VAR,VARIANCIA
        REAL(4), EXTERNAL  :: VMEAN3D,VAR3D
        REAL(4),DIMENSION(MKL) :: THETAUX
!        EXTERNAL           :: DGESV
        REAL(4)            :: WT
        INTEGER            :: ISTAT,IN_FILE,NELEM
!
! SET VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        SIG2 = SIG * SIG
        NELEM= NX * NY * NZ
        ALLOCATE(KPERM(NELEM))
        ALLOCATE(XI(NELEM))
        EPSILON = 1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IDENTIFICACAO DOS PONTOS NO VETOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MX=NX-1
        MY=NY-1
        MZ=NZ-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        EXT='.dat'
        FILE_OUT=ADJUSTL(FILE_OUT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUMBER OF FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        DO M=INIF,NFILES+INIF
!
           CALL CPU_TIME(T_START)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GERACAO DA VA GAUSSIANA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.EQ.0)THEN
              DO I=1,MKL
                 THETA(I) = RANDOM_NORMAL()
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.EQ.1)THEN
              DO I=1,MKL
                 THETA(I)=THETA(I)+SIG*RANDOM_NORMAL()
              ENDDO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.EQ.2)THEN
              DO I=1,MKL
                 THETA(I)=SQRT(1.0-SIG2)*THETA(I)+SIG*RANDOM_NORMAL()
!                         +MEANTHETA*(1.D0-SQRT(1.0-SIG2))
              END DO
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.EQ.4)THEN
              WRITE(*,*)'DO NOT CHANGE LOADED THETA'
           END IF
! NORMALIZACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           AUX=SQRT(1.d0/VARIANCIA(THETA,MKL))
!           DO I=1,MKL
!              THETA(I)=AUX*THETA(I)
!           END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONDICIONAMENTO DO CAMPO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NCOND .GT. 0)THEN
              DO I=1,MKL
                 AUX = 0.0D0
                 DO J=1,MKL
                    AUX = AUX + MM(I,J) * THETA(J)
                 END DO
                 THETAUX(I) = MU(I) + AUX
              END DO
              THETA = THETAUX
           END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SAVE THE NEW THETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CALL WRITE_THETA()
!
           DO I=1,NELEM
              AUX=0.D0
              DO J=1,MKL
                 AUX=AUX+AVET(I,J)*THETA(J)
              ENDDO
              XI(I)=AUX
           ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           KPERM = BETA * EXP(RHO * (FIELDM + EPSILON * XI))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           WRITE(C,113)M
           C=ADJUSTL(C)
           NAME=TRIM(FILE_OUT)//TRIM(C)//TRIM(EXT)
           NAME=ADJUSTL(TRIM(NAME))
           WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRINT FIELD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CALL PRINT_OLD3D(XI,NAME)
           CALL PRINT_UT(KPERM,NELEM,NAME)
           CALL CPU_TIME(T_FINAL)
           TIME = T_FINAL-T_START
           WRITE(*,333)TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           MEANK = 0.0
           DO I = 1,NELEM
              MEANK = MEANK + KPERM(I)
           ENDDO
           MEANK = MEANK/REAL(NELEM)
           VARK = 0
           DO I = 1,NELEM
              VARK = VARK + (KPERM(I) - MEANK)**2
           ENDDO
           WRITE(*,444)MEANK,SQRT(VARK/DBLE(NELEM-1))
        ENDDO
        DEALLOCATE(KPERM)
!
444     FORMAT('MEAN OF PERM....: ',E10.4,/,'STD OF PERM.....: ',E10.4)
333     FORMAT('TIME TO GENERATE SAMPLE (seconds) =',F8.2)
111     FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
113     FORMAT(I5)
!
      END SUBROUTINE KLCOND3D
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
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',&
             ACTION='WRITE')
        CALL RANDOM_SEED(get=NSEED)
        WRITE(OUT_FILE,103)(SEED(I),I=1,8)
!
        CLOSE(OUT_FILE)
!
103     FORMAT(I10)
111     FORMAT('NAME OF SEED OUTPUT FILE:',A)
!
      END SUBROUTINE PRINTSEED
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
        IN_FILE = 325
        OPEN(UNIT=IN_FILE,FILE=FILE_THE,STATUS='UNKNOWN', &
             FORM='FORMATTED',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)THEN
           WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_THE
           STOP
        END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO I=1,MKL
           READ(IN_FILE,*)THETA(I)
        ENDDO
!
        CLOSE(UNIT=IN_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111     FORMAT(e15.8)
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
        DO I=1,MKL
           WRITE(OUT_FILE,111)THETA(I)
        ENDDO
!
        CLOSE(UNIT=OUT_FILE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111     FORMAT(e15.8)
!
      END SUBROUTINE WRITE_THETA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINT_OLD3D(XI,NOME)
        USE VARIAVEIS
!
        INTEGER             :: J,I,K,BAND,OUT_FILE,M,L,N
        REAL(4), DIMENSION(NX*NY*NZ) :: XI
        CHARACTER(LEN=256)  :: NAMEF,NOME
        CHARACTER(LEN=4)    :: EXT
        CHARACTER(LEN=5)    :: C
!
        BAND=192837465
        OUT_FILE=323
        NAMEF=TRIM(NOME)
        NAMEF=ADJUSTL(TRIM(NAMEF))
!        WRITE(*,111)N,NAMEF(1:LEN_TRIM(NAMEF))
        OPEN(UNIT=OUT_FILE,FILE=NAMEF,STATUS='REPLACE',ACTION='WRITE')
        WRITE(OUT_FILE,101)DIMX
        WRITE(OUT_FILE,101)DIMY
        WRITE(OUT_FILE,101)DIMZ
        WRITE(OUT_FILE,100)NX
        WRITE(OUT_FILE,100)NY
        WRITE(OUT_FILE,100)NZ
        WRITE(OUT_FILE,100)NFLAG
        WRITE(OUT_FILE,101)FBETA
        WRITE(OUT_FILE,100)2
        WRITE(OUT_FILE,100)2
        !
        N=1
        DO K=1,NZ
           L=K-1
           WRITE(OUT_FILE,100)L
           DO J=1,NY
              M=J-1
              WRITE(OUT_FILE,100)M
              WRITE(OUT_FILE,101)(XI(I),I=N,N+NX-1)
              WRITE(OUT_FILE,102)BAND
              N = N+NX
           END DO
        END DO
!
        CLOSE(OUT_FILE)
!
100     FORMAT(I7)
101     FORMAT(100000F10.5)
102     FORMAT(I9)
113     FORMAT(I5)
111     FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
!
      END SUBROUTINE PRINT_OLD3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE PRINT_OLD(VEC,VX,VY,NNX,NNY,NAME)
      USE VARIAVEIS
      INTEGER,   INTENT(IN) :: NNX,NNY
      REAL(4),   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
      CHARACTER(LEN=128)             :: NAME
      INTEGER                        :: J,I,BAND,OUT_FILE,MMY
!
      BAND=192837465
      OUT_FILE=321
!
      OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
      WRITE(OUT_FILE,103)DIMX
      WRITE(OUT_FILE,103)DIMY
      WRITE(OUT_FILE,100)NX
      WRITE(OUT_FILE,100)NY
      WRITE(OUT_FILE,100)NFLAG
      WRITE(OUT_FILE,103)FBETA
      WRITE(OUT_FILE,100)2
      WRITE(OUT_FILE,100)2
      
      MMY=NNY-1
      DO J=0,MMY
         WRITE(OUT_FILE,100)J
         WRITE(OUT_FILE,101)(VEC(I,J),I=0,NNX-1)
         WRITE(OUT_FILE,102)BAND
      ENDDO
!
      CLOSE(OUT_FILE)
!
100   FORMAT(I7)
101   FORMAT(10000F12.4)
102   FORMAT(I9)
103   FORMAT(1F12.4)
!
      END SUBROUTINE PRINT_OLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE PRINT_UT(XI,NELEM,NAME)
        INTEGER,   INTENT(IN)    :: NELEM
        REAL(4),DIMENSION(NELEM) :: XI
        CHARACTER(LEN=128)       :: NAME
        INTEGER                  :: J,OUT_FILE
!
        OUT_FILE=325
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
        DO J=1,NELEM
           WRITE(OUT_FILE,101)XI(J)
        END DO
!
        CLOSE(OUT_FILE)
!
101     FORMAT(e15.8)
!
      END SUBROUTINE PRINT_UT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
