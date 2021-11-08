!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!                 3D RANDOM FIELDS GENERATOR
!                USING KARHUNEN-LOEVE METHOD
! #################  CONDITIONAL VERSION  ###################
!
! NAME: MARCIO RENTES BORGES
! DATA: 21/04/2021
! 
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
  REAL(4)  :: FDX,FDY,FDZ,BASE,SIG
  REAL(4)  :: FBETA,FVAR,FMEAN,DIMX,DIMY,DIMZ
  REAL(4)  :: VARIAN,XMEAN,XVAR,FATCOND
  INTEGER  :: NX,NY,NZ,NFLAG,NERROR,INIF,NFILES
  INTEGER  :: NFILE,NDIM,NS
  LOGICAL  :: TFCOND
  INTEGER  :: NCOND,MKL,NPROPOSAL
  INTEGER, ALLOCATABLE,DIMENSION(:)    :: NSEED,PCOND
  REAL(4), ALLOCATABLE,DIMENSION(:,:,:):: PSI
  REAL(4), ALLOCATABLE,DIMENSION(:,:)  :: AVET,VET,PS
  REAL(4), ALLOCATABLE,DIMENSION(:)    :: X,Y,Z,AVAL,THETA
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
  REAL, DIMENSION(2)   :: TARRAY
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL ETIME(TARRAY,TSTART)
  ERROR_READ = 0
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL READ_DATA(ERROR_READ)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MEMORY ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NUM_ELEM = NX*NY*NZ
  ALLOCATE(PSI(0:NX-1,0:NY-1,0:NZ-1))
  ALLOCATE(PS(0:NX-1,0:NY-1))
  ALLOCATE(X(0:NX))
  ALLOCATE(Y(0:NY))
  ALLOCATE(Z(0:NZ))
  ALLOCATE(AVAL(MKL))
  ALLOCATE(AVET(NUM_ELEM,MKL))
  ALLOCATE(PCOND(MKL))
  ALLOCATE(THETA(MKL))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLEAR VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CLEAR3D(PSI,NX,NY,NZ)
  CALL CLEAR(X,NX)
  CALL CLEAR(Y,NY)
  CALL CLEAR(Z,NZ)
  CALL CLEAR(THETA,MKL-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONSTRUCT VECTORS OF POSITIONS !!!!!!!!!!!!!!!!!!!!!!!!
!$OMP PARALLEL PRIVATE(I)
!$OMP DO
  DO I=0,NX
     X(I)=(I*FDX)+FDX*.5
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
!
!$OMP PARALLEL PRIVATE(J)
!$OMP DO
  DO J=0,NY
     Y(J)=(J*FDY)+FDY*.5
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
!$OMP PARALLEL PRIVATE(K)
!$OMP DO
  DO K=0,NZ
     Z(K)=(K*FDZ)+FDZ*.5
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETS THE SEEDS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_SEED(put=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM FIELDS GENERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'LOADING EINGENPAIRS'
!
  CALL LOAD_AUTOV()
!
  CALL LOAD_THETA()
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM CONDITIONING AND GENERATION !!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'CONSTRUCTING       '
!
  CALL KLCOND3D()
!
  CALL ETIME(TARRAY,TFINAL)
  TIME = TFINAL - TSTART
  WRITE(*,331)TIME,TARRAY(1),TARRAY(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IMPRESSAO DA SEMENTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALL PRINTSEED()
! MEMORY FREE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DEALLOCATE(AVAL,AVET,PSI)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROR = 1
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
        READ(IN_FILE,9000)FBETA,VARIAN
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
        READ(IN_FILE,7000)FILE_VAL
        WRITE(*,7002)FILE_VAL
!
        READ(IN_FILE,7000)FILE_VET
        WRITE(*,7003)FILE_VET
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
        IF(NCOND > 0)THEN
           WRITE(*,*)'=========================='
           WRITE(*,*)'=========================='
           WRITE(*,*)'CONDITIONING OFF'
           WRITE(*,*)'=========================='
           WRITE(*,*)'=========================='
           STOP 345
        END IF
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
        NERROREAD = 1d0
        CLOSE(UNIT=IN_FILE)
!
        RETURN

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LOAD_AUTOV()
!
        USE VARIAVEIS, ONLY: AVET,NX,NY,NZ,MKL,FILE_VET
        IMPLICIT NONE
!
        INTEGER          :: NELEM
        INTEGER          :: ISTAT,IN_FILEA,IN_FILEV
        INTEGER          :: I,J,K,NREC
        REAL             :: TIME,T_START,T_FINAL
        REAL,DIMENSION(2):: TARRAY
        CHARACTER(LEN=3) :: TFILE
        LOGICAL          :: EFILE
        REAL(4),ALLOCATABLE,DIMENSION(:) :: V
!
        NELEM = NX*NY*NZ
        ALLOCATE(V(NELEM * MKL))
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
        CALL ETIME(TARRAY,T_START)
        T_START = TARRAY(1)
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'ERR')THEN
           WRITE(*,*)'ERROR ON READING INPUT FILE: ',FILE_VET
           STOP 100
        END IF
        IF(TFILE.EQ.'ASC')THEN
           IN_FILEV = 132
           FILE_VET = TRIM(ADJUSTL(FILE_VET))
           INQUIRE(FILE=FILE_VET,EXIST=EFILE)
           OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='UNKNOWN', ACTION='READ',IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
              STOP
           END IF
           DO I=1,NELEM
              READ(IN_FILEV,*)AVET(I,1:MKL)
           ENDDO
           CLOSE(IN_FILEV)
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(TFILE.EQ.'BIN')THEN
           NREC = NELEM*MKL*4
           IN_FILEV = 326
           OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='OLD', &
                ACCESS='DIRECT',RECL=NREC,IOSTAT=ISTAT)
           IF(ISTAT.NE.0)THEN
              WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
              STOP
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
        CALL ETIME(TARRAY,T_FINAL)
        T_FINAL = TARRAY(1)
        TIME = T_FINAL-T_START
        WRITE(*,334)TIME,TARRAY(1),TARRAY(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
334     FORMAT('TIME TO LOAD AUT (s) =',F12.2,/,&
             'USER TIME        (s) =',F12.2,/,&
             'SYSTEM TIME      (s) =',F12.2,/)
111     FORMAT(e25.8)
112     FORMAT(40000e15.7)
!
      END SUBROUTINE LOAD_AUTOV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KLCOND3D()
!
        USE VARIAVEIS
        USE RANDOM
        IMPLICIT NONE
!
        INTEGER :: PNODE(1:NCOND),K,KK
        INTEGER :: NRHS,LDA,LDB,INFO
        INTEGER :: IPIV(NCOND,1)
        REAL(4) :: B(NCOND,1),MAT(NCOND,NCOND)
        REAL(4) :: AUX1,VARR,AUXSQ
        REAL(4) :: AUXM,AUXV
        REAL(4) :: XI(1:NX*NY*NZ)
        INTEGER :: I,J,M,MX,MY,MZ
        REAL(4) :: POSIX,POSIY,POSIZ
        REAL(4) :: AUX, SIG2
        INTEGER :: MAX,MIN,MIN2
        CHARACTER(LEN=256) :: NAME
        CHARACTER(LEN=4)   :: EXT
        CHARACTER(LEN=5)   :: C
        REAL               :: TIME,T_START,T_FINAL
        REAL, DIMENSION(2) :: TARRAY
        REAL(4), EXTERNAL  :: MAIOR,MENOR
        REAL(4), EXTERNAL  :: VMEAN,VAR,VARIANCIA
        REAL(4), EXTERNAL  :: VMEAN3D,VAR3D
!        EXTERNAL           :: DGESV
        REAL(4)            :: WT
        INTEGER            :: ISTAT,IN_FILE
!
!        EXTERNAL :: ETIME
! SET VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        SIG2 = SIG * SIG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IDENTIFICACAO DOS PONTOS NO VETOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO M=0,NCOND-1
           K=0
           DO KK=1,NZ
              POSIZ = REAL(KK)*FDZ
              IF((VET(2,M).LT.POSIZ).AND.(VET(2,M).GT.POSIZ-FDZ))THEN
                 DO J=NY,1,-1
                    POSIY = REAL(NY-J)*FDY
                    IF((VET(1,M).LT.POSIY+FDY.AND.VET(1,M).GT.POSIY))THEN
                       DO I=1,NX
                          POSIX = REAL(I)*FDX
                          IF((VET(0,M).LT.POSIX.AND.VET(0,M).GT.&
                               POSIX-FDX))THEN
                             K=K+1
                             PNODE(M+1)=K
                          ELSE
                             K=K+1
                          END IF
                       END DO
                    ELSE
                       K=K+NX
                    END IF
                 END DO
              ELSE
                 K=K+NX*NY
              END IF
           END DO
        END DO
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
           CALL ETIME(TARRAY,T_START)
           T_START = TARRAY(1)
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
!
           VARR=1.d0
           NRHS=1
           LDA =NCOND
           LDB =MAX(1,NCOND)
           IF(NCOND>0)THEN
              DO I=1,NCOND
                 B(I,1)=0.D0
              ENDDO
              DO I=1,NCOND
                 DO J=NCOND+1,MKL
                    B(I,1)=B(I,1)-&
                         AVET(PNODE(I),PCOND(J))*THETA(PCOND(J))
                 ENDDO
                 B(I,1)=B(I,1)+VET(3,I-1)+AUXSQ-VARR*0.5
                 DO J=1,NCOND
                    MAT(I,J)=AVET(PNODE(I),PCOND(J))
                 ENDDO
              ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! resolucao do sitema linear ! Aa=b !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              CALL DGESV(NCOND,NRHS,MAT,LDA,IPIV,B,LDB,INFO)
!              IF(INFO.EQ.0)WRITE(*,*)'successful exit'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              DO I=1,NCOND
                 THETA(PCOND(I))=B(I,1)
              ENDDO
           ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! SAVE THE NEW THETA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CALL WRITE_THETA()
!
           DO I=1,NX*NY*NZ
              AUX=0.D0
              DO J=1,MKL
                 AUX=AUX+AVET(I,J)*THETA(J)
              ENDDO
              XI(I)=AUX
           ENDDO
!
           K=0
           DO KK=0,MZ
              DO J=0,MY
                 DO I=0,MX
                    K=K+1
                    PSI(I,J,KK)=SQRT(VARIAN)*XI(K)
                 ENDDO
              ENDDO
           END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           WRITE(C,113)M
           C=ADJUSTL(C)
           NAME=TRIM(FILE_OUT)//TRIM(C)//TRIM(EXT)
           NAME=ADJUSTL(TRIM(NAME))
           WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRINT FIELD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CALL PRINTASC(PS,X,Y,NX,NY,NAME)
           CALL PRINT_OLD3D(XI,NAME)
!              CALL PRINT_OLDUW(PS,X,Y,NX,NY,NAME)
!
!          CALL PRINTBIN(PSI,X,Y,NX,NY,NAME)
!              WRITE(*,112)VMEAN(PS,NX,NY),VAR(PS,NX,NY)
           CALL ETIME(TARRAY,T_FINAL)
           T_FINAL = TARRAY(1)
           TIME = T_FINAL-T_START
           WRITE(*,333)TIME,TARRAY(1),TARRAY(2)
        ENDDO
        WRITE(*,112)VMEAN3D(PSI,NX,NY,NZ),VAR3D(PSI,NX,NY,NZ)
!
333     FORMAT('TIME TO GENERATE (s) =',F12.2,/,&
               'USER TIME        (s) =',F12.2,/,&
               'SYSTEM TIME      (s) =',F12.2,/)
111     FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
112     FORMAT('MEAN=',F12.4,'   VARIANCE=',F12.4)
113     FORMAT(I5)
!
      END SUBROUTINE KLCOND3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
      CALL RANDOM_SEED(get=NSEED)
      WRITE(OUT_FILE,103)(SEED(I),I=1,8)
!
      CLOSE(OUT_FILE)
!
103   FORMAT(I10)
111   FORMAT('NAME OF SEED OUTPUT FILE:',A)
!
      END SUBROUTINE PRINTSEED
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
      INTEGER,            INTENT(IN) :: NNX,NNY
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
      SUBROUTINE PRINT_OLDUW(VEC,VX,VY,NNX,NNY,NAME)
        USE VARIAVEIS
        INTEGER,            INTENT(IN) :: NNX,NNY
        REAL(4),   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE,MMY
!
        OUT_FILE=325
!
        OPEN(UNIT=OUT_FILE,FILE=NAME,STATUS='REPLACE',ACTION='WRITE')
        WRITE(OUT_FILE,*)'# FIELD'
        WRITE(OUT_FILE,100)NX,NY
      
        MMY=NNY-1
        DO J=0,MMY
           DO I=0,NNX-1
              WRITE(OUT_FILE,101)VEC(I,J)
           END DO
        END DO
!
        CLOSE(OUT_FILE)
!
100     FORMAT(I4,' x ',I4,' x 1')
101     FORMAT(e15.8)
!
      END SUBROUTINE PRINT_OLDUW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
