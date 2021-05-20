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
! DATA: 04/11/2013
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
  DOUBLE PRECISION  :: FDX,FDY,FDZ,BASE,SIG
  DOUBLE PRECISION  :: XLBD,XNEWLBD,ZMAX
  DOUBLE PRECISION  :: FBETA,FVAR,FMEAN,DIMX,DIMY,DIMZ
  DOUBLE PRECISION  :: VARIAN,XMEAN,XVAR,FATCOND
  INTEGER           :: NX,NY,NZ,NFLAG,NERROR,INIF,NFILES
  INTEGER           :: NFILE,NDIM,NS,NZMIN
  LOGICAL           :: TFCOND
  INTEGER           :: NCOND,MKL,NPROPOSAL
  INTEGER         , ALLOCATABLE,DIMENSION(:)    :: NSEED,PCOND
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:,:):: PSI
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:,:)  :: AVET,VET,PS
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:)    :: X,Y,Z,AVAL,THETA
  DOUBLE PRECISION, ALLOCATABLE,DIMENSION(:)    :: OTIME,XTIME
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
  INTEGER                    :: I,J,K,ERROR_READ,NUM_ELEM
  INTEGER, ALLOCATABLE       :: SEED(:)
  DOUBLE PRECISION, EXTERNAL :: POWERLAW,EXPONENTIAL
  REAL                       :: TIME
  REAL, DIMENSION(2)         :: TARRAY
!
  ERROR_READ = 0
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL READ_DATA(ERROR_READ)
! VERICACAO DA INTERPOLACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL VERIFICA()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MEMORY ALLOCATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NUM_ELEM = NX*NY*NZ
  ALLOCATE(PSI(0:NX-1,0:NY-1,0:NZ-1))
  ALLOCATE(PS(0:NX-1,0:NY-1))
  ALLOCATE(X(0:NX))
  ALLOCATE(Y(0:NY))
  ALLOCATE(Z(0:NZ))
  ALLOCATE(AVAL(1:MKL))
  ALLOCATE(AVET(1:NUM_ELEM,1:MKL))
  ALLOCATE(PCOND(1:MKL))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLEAR VECTORS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CLEAR3D(PSI,NX,NY,NZ)
  CALL CLEAR(X,NX)
  CALL CLEAR(Y,NY)
  CALL CLEAR(Z,NZ)
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
  CALL ETIME(TARRAY,TIME)
!
  CALL LOAD_AUTOV(AVAL,AVET,NUM_ELEM,MKL)
!
  CALL LOAD_THETA()
!
  CALL ETIME(TARRAY,TIME)
  WRITE(*,333)TIME,TARRAY(1),TARRAY(2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RANDOM CONDITIONING AND GENERATION !!!!!!!!!!!!!!!!!!!!
  WRITE(*,*)'CONSTRUCTING       '
  CALL ETIME(TARRAY,TIME)
!
  CALL KLCOND()
!
  CALL ETIME(TARRAY,TIME)
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
!
        IN_FILE = 123
        FILE_IN = 'in/entrada.in'
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
        READ(IN_FILE,*)XLBD
        WRITE(*,778)XLBD
!
        READ(IN_FILE,*)XNEWLBD
        WRITE(*,777)XNEWLBD
!
        READ(IN_FILE,*)NZMIN,ZMAX
        WRITE(*,779)NZMIN,ZMAX
!
        READ(IN_FILE,7000)FILE_OUT
        WRITE(*,7001)FILE_OUT
!
!        READ(IN_FILE,7000)FILE_VAL
!        WRITE(*,7002)FILE_VAL
!
        READ(IN_FILE,7000)FILE_VET
        WRITE(*,7003)FILE_VET
!
        READ(IN_FILE,7000)FILE_THE
        WRITE(*,7004)FILE_THE
!
        READ(IN_FILE,12)NPROPOSAL,SIG
        WRITE(*,13)NPROPOSAL,SIG
        IF(NPROPOSAL.EQ.1)WRITE(*,*)'RANDOM WALK'
        IF(NPROPOSAL.EQ.2)WRITE(*,*)'pCN'
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
        NERROREAD = 1d0
        CLOSE(UNIT=IN_FILE)
!
        RETURN
!
777     FORMAT('COMPRIMENTO DE CORRELACAO NOVO.....:',F10.3)
778     FORMAT('COMPRIMENTO DE CORRELACAO ORIGINAL.:',F10.3)
779     FORMAT('MALHA MINIMA NO TEMPO..............:',I10,/,&
               'DIMENSAO NO TEMPO..................:',F10.3)
96      FORMAT('NUMERO DE SEMENTES:',I2,/,'SEMENTES:',20I12)
97      FORMAT(20I12)
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
        DOUBLE PRECISION   :: VEC(0:NNX-1,0:NNY-1)
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
        DOUBLE PRECISION   :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
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
        DOUBLE PRECISION   :: VEC(0:NN)
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
        INTEGER,            INTENT(IN) :: NNX,NNY
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE
!
        BAND=192837465
        OUT_FILE=321
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
        INTEGER,            INTENT(IN) :: NNX,NNY
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
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
        DOUBLE PRECISION, INTENT(IN) :: C,LD
!        DOUBLE PRECISION             :: AUX,AUX2
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION MEAN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DOUBLE PRECISION FUNCTION VMEAN3D(VEC,NNX,NNY,NNZ)
!
        DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
        DOUBLE PRECISION             :: SUM
        INTEGER, INTENT(IN)          :: NNX,NNY,NNZ
        INTEGER                      :: I,J,K
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
      DOUBLE PRECISION FUNCTION VAR3D(VEC,NNX,NNY,NNZ)
!
        DOUBLE PRECISION, INTENT(IN) :: VEC(0:NNX-1,0:NNY-1,0:NNZ-1)
        DOUBLE PRECISION             :: SUM,AUX
        DOUBLE PRECISION, EXTERNAL   :: VMEAN3D
        INTEGER, INTENT(IN)          :: NNX,NNY,NNZ
        INTEGER                      :: I,J,K
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
      SUBROUTINE LOAD_AUTOV(VAL,VE,NELEM,M)
!
      USE VARIAVEIS
      IMPLICIT NONE
!
      DOUBLE PRECISION :: VAL(1:M)
      DOUBLE PRECISION :: VE(1:NELEM,1:M)
      INTEGER          :: NELEM,M
      INTEGER          :: ISTAT,IN_FILEA,IN_FILEV
      INTEGER          :: I,J
!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IN_FILEA = 321
!      OPEN(UNIT=IN_FILEA,FILE=FILE_VAL,STATUS='UNKNOWN', &
! FORM='FORMATTED',IOSTAT=ISTAT)
!      IF(ISTAT.NE.0)THEN
!         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VAL
!         STOP
!      END IF
!
      IN_FILEV = 323
      OPEN(UNIT=IN_FILEV,FILE=FILE_VET,STATUS='UNKNOWN', &
 FORM='FORMATTED',IOSTAT=ISTAT)
      IF(ISTAT.NE.0)THEN
         WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_VET
         STOP
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      DO I=1,M
!         READ(IN_FILEA,111)VAL(I)
!      ENDDO
!
      DO I=1,NELEM
         !READ(IN_FILEV,112)(VE(I,J),J=1,M)
         READ(IN_FILEV,*)(VE(I,J),J=1,M)
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
111   FORMAT(e25.8)
112   FORMAT(40000e15.7)
!
      END SUBROUTINE LOAD_AUTOV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE KLCOND()
!
        USE VARIAVEIS
        USE RANDOM
        IMPLICIT NONE
!
        INTEGER                       :: PNODE(1:NCOND),K,KK
        INTEGER                       :: NRHS,LDA,LDB,INFO
        INTEGER                       :: IPIV(NCOND,1)
        DOUBLE PRECISION              :: B(NCOND,1),MAT(NCOND,NCOND)
        DOUBLE PRECISION              :: AUX1,VARR,AUXSQ
        DOUBLE PRECISION              :: AUXM,AUXV
        DOUBLE PRECISION              :: XI(1:NX*NY*NZ)
        INTEGER                       :: I,J,M,MX,MY,MZ
        DOUBLE PRECISION              :: POSIX,POSIY,POSIZ
        DOUBLE PRECISION              :: AUX
        INTEGER                       :: MAX,MIN,MIN2
        CHARACTER(LEN=256)            :: NAME
        CHARACTER(LEN=4)              :: EXT
        CHARACTER(LEN=5)              :: C
        REAL                          :: TIME,T_START,T_FINAL
        REAL            , DIMENSION(2):: TARRAY
        DOUBLE PRECISION, EXTERNAL    :: MAIOR,MENOR
        DOUBLE PRECISION, EXTERNAL    :: VMEAN,VAR,VARIANCIA
        DOUBLE PRECISION, EXTERNAL    :: VMEAN3D,VAR3D
        EXTERNAL                      :: DGESV
        DOUBLE PRECISION              :: WT
        INTEGER                       :: ISTAT,IN_FILE
!
!        EXTERNAL :: ETIME
! SET VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IDENTIFICACAO DOS PONTOS NO VETOR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL ETIME(TARRAY,T_START)
        T_START = TARRAY(1)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        AUX=0.D0
        DO I=1,MKL
           AUX1=AVET(1,I)
           AUX=AUX+0.5*AUX1*AUX1
        ENDDO
        AUXSQ=AUX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL ETIME(TARRAY,T_FINAL)
        T_FINAL = TARRAY(1)
        TIME = T_FINAL-T_START
        WRITE(*,334)TIME,TARRAY(1),TARRAY(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        MX=NX-1
        MY=NY-1
        MZ=NZ-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        EXT='.dat'
        FILE_OUT=ADJUSTL(FILE_OUT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CONDITIONING ==> BEST CONDITIONING !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF(NCOND.GT.0)THEN
           FILE_COND = 'in/cond_seq.in'
           FILE_COND = ADJUSTL(TRIM(FILE_COND))
           INQUIRE(FILE=FILE_COND,EXIST=TFCOND)
           WRITE(*,*)'#########################'
           IF(TFCOND)THEN
              WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FILE_COND))
              IN_FILE = 22
              OPEN(UNIT=IN_FILE,FILE=FILE_COND,STATUS='OLD',&
                   ACTION='READ',FORM='FORMATTED',IOSTAT=ISTAT)
              IF(ISTAT.NE.0)THEN
                 WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_COND
                 STOP
              END IF
              DO I=1,MKL
                 READ(IN_FILE,*)PCOND(I)
              END DO
              WRITE(*,*)'SEQUENCY FOR CONDITIONING: ',PCOND(1:NCOND)
              CLOSE(IN_FILE)
           ELSE
              WRITE(*,*)'CONVENTIONAL CONDITIONING'
              DO I=1,MKL
                 PCOND(I)=I
              END DO
           END IF
           WRITE(*,*)'#########################'
        END IF
! MAIN LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NUMBER OF FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
        DO M=INIF,NFILES+INIF
!
           CALL ETIME(TARRAY,T_START)
           T_START = TARRAY(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GERACAO DA VA GAUSSIANA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF(NPROPOSAL.EQ.1)THEN
              DO I=1,MKL
                 THETA(I)=THETA(I)+SIG*RANDOM_NORMAL()
              ENDDO
! NORMALIZACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              AUX=SQRT(1.d0/VARIANCIA(THETA,MKL))
              DO I=1,MKL
                 THETA(I)=AUX*THETA(I)
              END DO
           END IF
           IF(NPROPOSAL.EQ.2)THEN
              DO I=1,MKL
                 THETA(I)=SQRT(1.0-SIG*SIG)*THETA(I)+SIG*RANDOM_NORMAL()
!                         +MEANTHETA*(1.D0-SQRT(1.0-SIG*SIG))
              END DO
           END IF
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
              CALL DGESV(NCOND,NRHS,MAT,LDA,IPIV,B,LDB,INFO)
              IF(INFO.EQ.0)THEN
                 WRITE(*,*)'#########################'
                 WRITE(*,*)'#### successful exit ####'
                 WRITE(*,*)'#########################'
              ELSE
                 CALL INFORMACAO(MAT,INFO)
                 STOP
              END IF
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
              XI(I)=0.5*VARR+AUX-AUXSQ
           ENDDO
!$OMP PARALLEL PRIVATE(I,J,KK)
!$OMP DO
           K=0
           DO KK=0,MZ
              DO J=0,MY
                 DO I=0,MX
                    K=K+1
                    PSI(I,J,KK)=SQRT(VARIAN)*XI(K)
                 ENDDO
              ENDDO
           END DO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           CALL INTERPOLATION()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NORMALIZACAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           AUXM=VMEAN3D(PSI,NX,NY,NZ)
!           WRITE(*,*)AUXM
!           AUXV=1.D0/SQRT(VAR3D(PSI,NX,NY,NZ))
!$OMP PARALLEL PRIVATE(I,J,KK)
!$OMP DO
!           DO KK=0,MZ
!              DO J=0,MY
!                 DO I=0,MX
!                    PSI(I,J,KK)=AUXV*(PSI(I,J,KK)-AUXM)
!                 ENDDO
!              ENDDO
!           END DO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           DO KK=1,NZMIN
              WRITE(C,113)KK-1
              C=ADJUSTL(C)
              NAME=TRIM(FILE_OUT)//TRIM(C)//TRIM(EXT)
              NAME=ADJUSTL(TRIM(NAME))
!              WRITE(*,111)M,NAME(1:LEN_TRIM(NAME))
!$OMP PARALLEL PRIVATE(I,J)
!$OMP DO
              DO J=0,MY
                 DO I=0,MX
                    PS(I,J)=PSI(I,J,KK-1)
                 ENDDO
              ENDDO
!$OMP END DO
!$OMP END PARALLEL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PRINT FIELD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           CALL PRINTASC(PS,X,Y,NX,NY,NAME)
              CALL PRINT_OLD(PS,X,Y,NX,NY,NAME)
!              CALL PRINT_OLDUW(PS,X,Y,NX,NY,NAME)
!
!          CALL PRINTBIN(PSI,X,Y,NX,NY,NAME)
!              WRITE(*,112)VMEAN(PS,NX,NY),VAR(PS,NX,NY)
           END DO
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
334     FORMAT('TIME TO LOAD AUT (s) =',F12.2,/,&
               'USER TIME        (s) =',F12.2,/,&
               'SYSTEM TIME      (s) =',F12.2,/)
111     FORMAT('NAME OF OUTPUT FILE',I5,': ',A)
112     FORMAT('MEAN=',F12.4,'   VARIANCE=',F12.4)
113     FORMAT(I5)
!
      END SUBROUTINE KLCOND
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
!! ALOCACAO DA MEMORIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ALLOCATE(THETA(MKL))
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        OUT_FILE = 326
        FILETHE='./out/thetanew.dat'
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
      SUBROUTINE PRINT_OLD3D(XI,N)
        USE VARIAVEIS
!
        INTEGER, INTENT(IN) :: N
        INTEGER             :: J,I,K,BAND,OUT_FILE,M,L
        DOUBLE PRECISION, DIMENSION(NX*NY*NZ) :: XI
        CHARACTER(LEN=256)  :: NAMEF,NOME
        CHARACTER(LEN=4)    :: EXT
        CHARACTER(LEN=5)    :: C
!
        BAND=192837465
        OUT_FILE=321
!
        NOME=FILE_OUT
        EXT='.dat'
        NAMEF=ADJUSTL(NAMEF)
        WRITE(C,113)N
        C=ADJUSTL(C)
        NAMEF=TRIM(NOME)//TRIM(C)//TRIM(EXT)
        NAMEF=ADJUSTL(TRIM(NAMEF))
        WRITE(*,111)N,NAMEF(1:LEN_TRIM(NAMEF))
        OPEN(UNIT=OUT_FILE,FILE=NAMEF,STATUS='REPLACE',ACTION='WRITE')
        WRITE(OUT_FILE,101)FDX
        WRITE(OUT_FILE,101)FDY
        WRITE(OUT_FILE,101)FDZ
        WRITE(OUT_FILE,100)NX
        WRITE(OUT_FILE,100)NY
        WRITE(OUT_FILE,100)NZ
        WRITE(OUT_FILE,100)NFLAG
        WRITE(OUT_FILE,101)FBETA
        WRITE(OUT_FILE,100)2
        WRITE(OUT_FILE,100)2
!
        DO K=1,NZ
           DO J=1,NY
              M=J-1
              L=K-1
              WRITE(OUT_FILE,100)L,M
              WRITE(OUT_FILE,101)(XI(I),I=L*NY*NX+M*NX+1,J*NX+L*NY*NZ)
              WRITE(OUT_FILE,102)BAND
           END DO
        END DO
!
        CLOSE(OUT_FILE)
!
100     FORMAT(2I7)
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
      DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
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
        DOUBLE PRECISION,   INTENT(IN) :: VEC(0:NNX-1,0:NNY-1),VX(0:NNX),VY(0:NNY)
        CHARACTER(LEN=128)             :: NAME
        INTEGER                        :: J,I,BAND,OUT_FILE,MMY
!
        OUT_FILE=321
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
      SUBROUTINE INTERPOLATION()
!
        USE VARIAVEIS, ONLY:PSI,XNEWLBD,XLBD,FDZ,DIMZ
        USE VARIAVEIS, ONLY:NX,NY,NZ,NZMIN,ZMAX
        DOUBLE PRECISION :: SCALINGF,Y
        INTEGER          :: I,J,K,KMAX
        DOUBLE PRECISION,DIMENSION(0:NZ-1):: YN,B,C,D
        DOUBLE PRECISION,EXTERNAL         :: ISPLINE
        DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:):: OTIME,XTIME 
!
        SCALINGF=XNEWLBD/XLBD
!
        ALLOCATE(OTIME(0:NZ-1))
        ALLOCATE(XTIME(0:NZMIN-1))
        DO K=0,NZ-1
           OTIME(K)=DBLE(K)*SCALINGF*FDZ
        END DO
        DO K=0,NZMIN-1
           XTIME(K)=DBLE(K)*ZMAX/(DBLE(NZMIN))
        END DO
        DO J=0,NY-1
           DO I=0,NX-1
              DO K=0,NZ-1
                 YN(K)=PSI(I,J,K)
              ENDDO
              CALL SPLINE(OTIME,YN,B,C,D,NZ)
              DO K=0,NZMIN-1
                 Y=ISPLINE(XTIME(K),OTIME,YN,B,C,D,NZ)
                 PSI(I,J,K)=Y
              END DO
           ENDDO
        END DO
!
        DEALLOCATE(OTIME,XTIME)
!
      END SUBROUTINE INTERPOLATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
        implicit none
        integer n
        double precision x(n), y(n), b(n), c(n), d(n)
        integer i, j, gap
        double precision h

        gap = n-1
        ! check input
        if ( n < 2 ) return
        if ( n < 3 ) then
           b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
           c(1) = 0.
           d(1) = 0.
           b(2) = b(1)
           c(2) = 0.
           d(2) = 0.
           return
        end if
!
! step 1: preparation
!
        d(1) = x(2) - x(1)
        c(2) = (y(2) - y(1))/d(1)
        do i = 2, gap
           d(i) = x(i+1) - x(i)
           b(i) = 2.0*(d(i-1) + d(i))
           c(i+1) = (y(i+1) - y(i))/d(i)
           c(i) = c(i+1) - c(i)
        end do
!
! step 2: end conditions 
!
        b(1) = -d(1)
        b(n) = -d(n-1)
        c(1) = 0.0
        c(n) = 0.0
        if(n /= 3) then
           c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
           c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
           c(1) = c(1)*d(1)**2/(x(4)-x(1))
           c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
        end if
!
! step 3: forward elimination 
!
        do i = 2, n
           h = d(i-1)/b(i-1)
           b(i) = b(i) - h*d(i-1)
           c(i) = c(i) - h*c(i-1)
        end do
!
! step 4: back substitution
!
        c(n) = c(n)/b(n)
        do j = 1, gap
           i = n-j
           c(i) = (c(i) - d(i)*c(i+1))/b(i)
        end do
!
! step 5: compute spline coefficients
!
        b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
        do i = 1, gap
           b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
           d(i) = (c(i+1) - c(i))/d(i)
           c(i) = 3.*c(i)
        end do
        c(n) = 3.0*c(n)
        d(n) = d(n-1)
      end subroutine spline
!======================================================================
!======================================================================
!======================================================================
      function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
        implicit none
        double precision ispline
        integer n
        double precision  u, x(n), y(n), b(n), c(n), d(n)
        integer i, j, k
        double precision dx

! if u is ouside the x() interval take a boundary value (left or right)
        if(u <= x(1)) then
           ispline = y(1)
           return
        end if
        if(u >= x(n)) then
           ispline = y(n)
           return
        end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
        i = 1
        j = n+1
        do while (j > i+1)
           k = (i+j)/2
           if(u < x(k)) then
              j=k
           else
              i=k
           end if
        end do
!*
!  evaluate spline interpolation
!*
        dx = u - x(i)
        ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      end function ispline
!=============================================================================
!=============================================================================
!=============================================================================
      SUBROUTINE INFORMACAO(A,inform)
        USE VARIAVEIS, ONLY: NCOND
        IMPLICIT NONE
        DOUBLE PRECISION :: A(NCOND,NCOND)
        INTEGER          :: I,J,inform
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'### PROBLEMA NA RESOLUCAO DO SISTEMA LINEAR - INFO= ',inform
        WRITE(*,*)'A positive value of INFO on return from an LAPACK routine'
        WRITE(*,*)'indicates a failure in the course of the algorithm. Common'
        WRITE(*,*)'causes are:'
        WRITE(*,*)'--> a matrix is singular (to working precision);'
        WRITE(*,*)'--> a symmetric matrix is not positive definite;'
        WRITE(*,*)'--> an iterative algorithm for computing eigenvalues or '
        WRITE(*,*)'    eigenvectors fails to converge in the permitted number '
        WRITE(*,*)'    of iterations.'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'### MATRIZ  ###############################################'
        DO I=1,NCOND
           WRITE(*,*)(A(I,J),J=1,NCOND)
        END DO
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'DETALHES: http://www.netlib.org/lapack/lug/node138.html'
        WRITE(*,*)'###########################################################'
        WRITE(*,*)'###########################################################'
      END SUBROUTINE INFORMACAO
!=============================================================================
!=============================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE VERIFICA()
!
        USE VARIAVEIS, ONLY:XNEWLBD,XLBD,FDZ,DIMZ
        USE VARIAVEIS, ONLY:NZ,NZMIN,ZMAX
        DOUBLE PRECISION :: SCALINGF
!
        SCALINGF=XNEWLBD/XLBD
        WRITE(*,100)SCALINGF
!
        IF(SCALINGF*DIMZ.LT.ZMAX)THEN
           WRITE(*,*)'###########################'
           WRITE(*,101)XNEWLBD
           WRITE(*,200)
           WRITE(*,300)(ZMAX/DIMZ)*XLBD
           WRITE(*,400)(ZMAX/(FDZ*0.5*DBLE(NZMIN)))*XLBD
           WRITE(*,*)'###########################'
           STOP
        END IF
!
        IF(FDZ*DBLE(NZMIN)*SCALINGF*0.5.GE.ZMAX)THEN
           WRITE(*,*)'###########################'
           WRITE(*,101)XNEWLBD
           WRITE(*,201)
           WRITE(*,300)(ZMAX/DIMZ)*XLBD
           WRITE(*,400)(ZMAX/(FDZ*0.5*DBLE(NZMIN)))*XLBD
           WRITE(*,*)'###########################'
           STOP
        END IF
!
100     FORMAT('FATOR DE ESCALA....................:',F10.1)
101     FORMAT('COMPRIMENTO DE CORRELACAO ATUAL....:',F10.1)
200     FORMAT('PROBLEMA NA INTERPOLACAO (comp. corr. menor)')
201     FORMAT('PROBLEMA NA INTERPOLACAO (comp. corr. maior)')
300     FORMAT('COMPRIMENTO DE CORRELACAO MINIMO...:',F10.1)
400     FORMAT('COMPRIMENTO DE CORRELACAO MAXIMO...:',F10.1)
!
      END SUBROUTINE VERIFICA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
