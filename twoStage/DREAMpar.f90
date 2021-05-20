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
! NAME: MARCIO RENTES BORGES
! DATA: 09/03/2018
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!------------------------------------------------------------
!------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE VARIAVEIS
  CHARACTER(LEN=128)                :: FILEIN,NAME_OUT,FILEINI
  CHARACTER(LEN=128)                :: FILERROR
  CHARACTER(LEN=128),DIMENSION(1:10):: FILE_OUT,FILE_EXT,FILE_VAL
  CHARACTER(LEN=128),DIMENSION(1:10):: FILE_VET,FILE_THE
  CHARACTER(LEN=128),DIMENSION(1:10):: FILE_INFIELD,FILEREF
  CHARACTER(LEN=128),DIMENSION(1:10):: FILEAM,NAMEPRIOR
  CHARACTER(LEN=128),DIMENSION(1:10,1:40):: NAME_AVET
  CHARACTER(LEN=128),ALLOCATABLE,DIMENSION(:,:):: FILE_OUTEXP
  INTEGER            :: NSTAGE,NPRIORR,NPRINT,NDTYPE,ACCEPT
  INTEGER            :: NR,NRTOTAL,NSIMUL,NS,NCONDMAX,NCHAIN
  INTEGER            :: SGN,CONTADORC,CONTADORF,CONTREJ,CONT
  INTEGER            :: NNORMA,NCUT
  INTEGER,ALLOCATABLE,DIMENSION(:) :: VETCONT
  INTEGER,ALLOCATABLE,DIMENSION(:) :: GERATIPO,LIKETIPO,NLPRIOR
  INTEGER,ALLOCATABLE,DIMENSION(:) :: NWELLS,NDADOS,NDADOSI,NSEED,NUPSC
  INTEGER,ALLOCATABLE,DIMENSION(:) :: NCOND,NX,NY,NZ,NFLAG,MKL,NTOTAL
  INTEGER,ALLOCATABLE,DIMENSION(:) :: NPROPOSAL,NUM_AVET,NUM_AVETOLD
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: ERRORCC
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: ERRORNC,ERRORCF,ERRORNF
  DOUBLE PRECISION   :: TERRORCC,TERRORNC,TERRORCF,TERRORNF
  DOUBLE PRECISION   :: TLRATIOF,TLRATIOC,TPRATIO
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: SIGMA2C,SIGMA2F
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: VET
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: ALPHA,DIMX
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: DIMY,DIMZ,SIG,VARIAN,FBETA
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: REFDATA,REFAMOS
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: ERRORF
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: LOGPERM,LOGPERMN
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: PRATIO,LRATIOF,LRATIOC
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: MAXDATA,NORMA_L2_REF
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: THETA,THETAN
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: NEWLBD,OLDLBD,LBDFIXO
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: YNEWLBD,YOLDLBD
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: XLBDMEAN,XLBDSD,XLBDBETA
  INTEGER         ,ALLOCATABLE,DIMENSION(:)     :: NLBD,NLBDVET,NZMIN,NKIND
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: XLBDS,CLBDS,ZMAX
  INTEGER                                       :: NLBDPOSI
  DOUBLE PRECISION                              :: YNEW,YOLD
  INTEGER                                       :: TIPOVAR, ERROSIMULADOR
  REAL,ALLOCATABLE,DIMENSION(:)                 :: XRMVN_MEAN,XRMVN
  REAL,ALLOCATABLE,DIMENSION(:,:)               :: XRMVN_VET,MATDE
  REAL,ALLOCATABLE,DIMENSION(:)                 :: COVMATUP,COVMATDOWN,COVMATID
  INTEGER                                       :: NDIM,NDIMR,NINICIO,NSTART
  INTEGER                                       :: NFREQ,NPCHAIN,NDELTA
  DOUBLE PRECISION                              :: SIGDE
  REAL                                          :: CR_DREAM
!
END MODULE VARIAVEIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE STATFUNCTIONS
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCAO PARA CALCULO DO VALOR DA PDF DE MVN !!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL function mv_normal_pdf(k,x,x_mean,sig)
!    use mkl95_lapack, only : potrf,trtrs
!    use mkl95_blas, only : dot
    implicit none
    external :: dpotrf,dtrtrs
    integer             :: kr,j,loc
    integer, intent(in) :: k
    double precision,intent(in),dimension(k) :: x, x_mean
    double precision,dimension(k) :: xminusmean
    double precision,intent(in out),dimension(k*(k+1)/2) :: sig
    double precision,dimension(k,k) :: sigma
    double precision :: rtdet_sigma
    doubleprecision :: c, inside_exp
    integer :: i,info,lda,ldb
!
    lda = k
    ldb = 2
    kr  = k*(k+1)/2
!
! Changing matrix format
!
    do i=1,k
       do j=i,k
          loc = j*(j-1)/2 + i
          sigma(i,j) = sig(loc)
          sigma(j,i) = sig(loc)
       end do
    end do
!
! Form Cholesky factors, sigma = GG^T
!
!    call dpotrf( sigma, 'L' ,info )
    call dpotrf('L',k,sigma,k ,info )
    if(info /= 0)then
       write(*,*)' INFO from DPOTRF = ',info
       stop
    endif
!
! det(G) is product of diagonal terms, and = sqrt(det(sigma))
!
    rtdet_sigma=product([(sigma(i,i),i=1,k)])
    xminusmean=x-x_mean
!
! compute y = inverse(G)(x-m)
!
!    call dtrtrs(sigma,xminusmean,'L','N','N',info)
    call dtrtrs('L','N','N',k,k,sigma,lda,xminusmean,ldb,info)
    if(info /= 0)then
       write(*,*)' INFO from TRTRS = ',info
       stop
    endif
!
! argument of exp() is -(1/2) y^T . y
!
    inside_exp=-0.5d0*mrbdot(xminusmean,xminusmean,k)
    c = (atan(1d0)*8d0)**(k/2d0)
    mv_normal_pdf= max(exp(inside_exp)/(c*rtdet_sigma),1d-16)
  end function mv_normal_pdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCAO PARA CALCULO DO PRODUTO INTERNO !!!!!!!!!!!!!!!  
  double precision function mrbdot(X,Y,K)
    IMPLICIT NONE
    INTEGER             :: I
    INTEGER, INTENT(IN) :: K
    DOUBLE PRECISION,DIMENSION(K),INTENT(IN):: X
    DOUBLE PRECISION,DIMENSION(K),INTENT(IN):: Y
    DOUBLE PRECISION    :: PRODI
!
    PRODI = 0.0
    DO I=1,K
       PRODI = PRODI + X(I)*Y(I)
    END DO
    MRBDOT = PRODI
  END function mrbdot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCAO PARA CALCULO DA MEDIA DE VETORES !!!!!!!!!!!!!!  
  FUNCTION MEANV(VET,D,TAM)
    INTEGER, INTENT(IN) :: D,TAM
    REAL, DIMENSION(D,TAM),INTENT(IN)    :: VET
    REAL, DIMENSION(D)  :: MEANV
    INTEGER             :: I,J
    !
    MEANV = 0.0D0
    DO J=1,D
       DO I=1,TAM
          MEANV(J) = MEANV(J)+VET(J,I)
       END DO
       MEANV(J)=MEANV(J)/REAL(TAM)
    END DO
    RETURN
  END FUNCTION MEANV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION CHOOSER(NK,NP)
    INTEGER, INTENT(IN)   :: NK,NP
    INTEGER               :: R1,R2
    INTEGER, DIMENSION(2) :: CHOOSER
    !
    R1 = NK
    R2 = NK
    DO WHILE(R1.EQ.NK)
       R1 = UNIDRND(NP)
    END DO
    DO WHILE(R2.EQ.NK.OR.R1.EQ.R2)
       R2 = UNIDRND(NP)
    END DO
    CHOOSER(1) = R1
    CHOOSER(2) = R2
    RETURN
  END FUNCTION CHOOSER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION UNIDRND(N)
    USE RANDOM
    IMPLICIT NONE
!
    INTEGER :: N
    REAL    :: AUX
!
    CALL RANDOM_NUMBER(AUX)
    AUX = REAL(N)*AUX+1.0
    UNIDRND = INT(AUX)
    RETURN
  END FUNCTION UNIDRND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION UNIFDISC(A,B)
    USE RANDOM
    IMPLICIT NONE
!
    INTEGER :: A,B,NA,NB,N
    INTEGER :: ID,I
    REAL    :: AUX,FA,FB
!
    NB = MAX(A,B)
    NA = MIN(A,B)
    N  = NB-NA+1
    CALL RANDOM_NUMBER(AUX)
    AUX = DFLOAT(NB-NA)*AUX+DFLOAT(NA) 
    UNIFDISC = NINT(AUX)
!
  END FUNCTION UNIFDISC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION MULTINOM(NC,P)
    USE RANDOM
    IMPLICIT NONE
!
    INTEGER :: NC,ID,I
    REAL    :: AUX,FA,FB
    REAL,DIMENSION(NC) :: P
!
    CALL RANDOM_NUMBER(AUX)
    FA = 0.0
    FB = 0.0
    ID = 1
    DO I=1,NC-1
       FA=FA+P(I)
       FB=FA+P(I+1)
       IF(AUX.GT.FA.AND.AUX.LE.FB)ID=I+1
    END DO
    MULTINOM = ID
    RETURN
!
  END FUNCTION MULTINOM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! FUNCAO PARA CALCULO DA COVARIANCIA !!!!!!!!!!!!!!!!!!!
  FUNCTION COVAR(VET,D,DR,TAM)
    INTEGER, INTENT(IN)   :: D,DR,TAM
    REAL, DIMENSION(D,TAM),INTENT(IN) :: VET
    REAL, DIMENSION(DR)   :: COVAR
    INTEGER               :: I,J,K,LOCJ,LOC
    REAL, DIMENSION(D)    :: XMEAN
    REAL                  :: SOMA, AUX
    !
    XMEAN = MEANV(VET,D,TAM)
    DO J=1,D
       LOCJ = J*(J-1)/2
       DO I=J,D
          LOC = LOCJ+I
          SOMA = 0.0E0
          DO K=1,TAM
             SOMA = SOMA + VET(I,K)*VET(J,K)
          END DO
          AUX = SOMA/REAL(TAM) - XMEAN(J)*XMEAN(I)
          COVAR(LOC) = AUX
       END DO
    END DO
    RETURN
  END FUNCTION COVAR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! FUNCAO PARA CALCULO DA COVARIANCIA !!!!!!!!!!!!!!!!!!!
  REAL FUNCTION STDV(VET,D,TAM)
    INTEGER, INTENT(IN)   :: D,TAM
    REAL, DIMENSION(D,TAM),INTENT(IN) :: VET
    INTEGER               :: I,J,K
    REAL, DIMENSION(D)    :: XMEAN
    REAL                  :: SOMA,SOMASQ, AUX
    !
    SOMA  = 0.0
    SOMASQ= 0.0
    DO J=1,D
       DO K=1,TAM
          SOMA   = SOMA + VET(J,K)
       END DO
    END DO
    SOMA  = SOMA/DFLOAT(D*TAM)
    DO J=1,D
       DO K=1,TAM
          SOMASQ = SOMASQ + (VET(J,K)-SOMA)**2
       END DO
    END DO
    SOMASQ= SOMASQ/DFLOAT(D*TAM-1)
    STDV = SQRT(SOMASQ)
!
  END FUNCTION STDV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE STATFUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
!
  USE VARIAVEIS
  USE RANDOM
  USE STATFUNCTIONS
  IMPLICIT NONE
  INCLUDE 'mpif.h'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LIST OF LOCAL VARIABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER            :: NERROREAD,K,RC,LEN,NK,NCONT
  REAL, DIMENSION(2) :: TARRAY,AUX
  CHARACTER(LEN=128) :: COMMAND
  INTEGER            :: NPROCS, IERROR, NRANK, I, J, TAG
  INTEGER            :: MPI_NSIMUL,MPISTATUS(MPI_STATUS_SIZE)
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) :: RANKNAME
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!  IF(NRANK.EQ.0)
  WRITE(*,201)NRANK,ADJUSTL(TRIM(RANKNAME))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD = 0
  CALL READ_INPUT(NERROREAD,MPI_NSIMUL,NPROCS,NRANK)
!
  IF(NERROREAD.NE.1)THEN
     WRITE(*,*)'#### LEITURA ERRADA ####'
     WRITE(*,*)'### DO ARQUIVO INPUT ###'
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NERROREAD = 0
  CALL INITIAL(NERROREAD,NRANK,NPROCS)
!
  IF(NERROREAD.NE.1)THEN
     WRITE(*,*)'######### PROBLEMA NA ###########'
     WRITE(*,*)'### INICIALIZACAO DO PROGRAMA ###'
     STOP
  END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NPROPOSAL(1).EQ.4)THEN
     VETCONT(NRANK+1) = CONT
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
     CALL MPI_ALLGATHER(CONT,1,MPI_INTEGER,VETCONT,1,&
          MPI_INTEGER,MPI_COMM_WORLD,IERROR)
     CALL SAVEMPICONT(CONT,1,NRANK)
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL COPY_DIR(NRANK,NSIMUL,NSTAGE)
  CALL COPY_DIRUP(NRANK,NSTAGE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NPRINT.EQ.6)GOTO 2233
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!#####################################################!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LOOP SOBRE AS REALIZACOES !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO WHILE(CONT.LT.NR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SALVA O STATUS ATUAL DA SIMULACAO !!!!!!!!!!!!!!!!!!!!!
     CALL RANDOM_SEED(put=NSEED)
!
     CALL SAVE_INITSTATUS(CONTADORC,CONTADORF,&
          CONT,NCHAIN,NRANK)
!
     CONTADORC=CONTADORC+1
!
999  CONTINUE
!
     IF(NPROPOSAL(1).EQ.4)CALL READMPICONT(VETCONT,NPROCS,NRANK)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(CONTADORC.GT.NRTOTAL)THEN
        IF(NRANK.EQ.0)WRITE(*,*)'#################################'
        IF(NRANK.EQ.0)WRITE(*,*)'### SAIDA SEM CONVERGENCIA ######'
        IF(NRANK.EQ.0)WRITE(*,*)'#################################'
        IF(NRANK.EQ.0)WRITE(*,332)REAL(CONT)/REAL(CONTADORF)*100.d0
        EXIT
     ENDIF
!
     IF(NRANK.EQ.NRANK)WRITE(*,*)'#########################################################'
     IF(NRANK.EQ.NRANK)WRITE(*,*)'######################## REALIZACAO:',CONTADORC,' #######'
     IF(NRANK.EQ.NRANK)WRITE(*,*)'######################## RANK......:',NRANK    ,' #######'
     IF(NRANK.EQ.NRANK)WRITE(*,*)'#########################################################'
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
           STOP
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
!       CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESOLUCAO DO PROBLEMA DE TRANSPORTE (COARSE SCALE) !!!
        CALL SIMULADOR_C(NRANK)
        IF(NSIMUL.EQ.5.AND.ERROSIMULADOR.EQ.1)THEN
           WRITE(*,1001)NRANK
           GOTO 999
        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
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
        ACCEPT=1D0
     END IF
!!#####################################################!!
!! END TWO STAGE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    CALL RANDOM_SEED(get=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     IF(ACCEPT.GT.0)THEN
        CONT=CONT+1
        IF(NPROPOSAL(1).EQ.4)THEN
           VETCONT(NRANK+1)=CONT
           CALL SAVEMPICONT(CONT,1,NRANK)
        END IF
     ELSE
        CONTREJ=CONTREJ+1
     END IF
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  WRITE(*,332)NRANK,REAL(CONT)/REAL(CONTADORF)*100.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  DEALLOCATE(VETCONT)
  DEALLOCATE(VET)
  DEALLOCATE(DIMX)
  DEALLOCATE(DIMY)
  DEALLOCATE(DIMZ)
  DEALLOCATE(NX)
  DEALLOCATE(NY)
  DEALLOCATE(NZ)
  DEALLOCATE(NFLAG)
  DEALLOCATE(MKL)
  DEALLOCATE(NTOTAL)
  DEALLOCATE(NCOND)
  DEALLOCATE(NPROPOSAL)
  DEALLOCATE(ALPHA)
  DEALLOCATE(VARIAN)
  DEALLOCATE(SIG)
  DEALLOCATE(FBETA)
  DEALLOCATE(NWELLS)
  DEALLOCATE(LIKETIPO)
  DEALLOCATE(NDADOS)
  DEALLOCATE(NDADOSI)
  DEALLOCATE(GERATIPO)
  DEALLOCATE(REFDATA)
  DEALLOCATE(REFAMOS)
  DEALLOCATE(THETA)
  DEALLOCATE(THETAN)
  DEALLOCATE(MAXDATA)
  DEALLOCATE(NORMA_L2_REF)
  DEALLOCATE(NLPRIOR)
  DEALLOCATE(NUM_AVET)
  DEALLOCATE(NUM_AVETOLD)
  DEALLOCATE(XRMVN_MEAN,XRMVN)
  DEALLOCATE(XRMVN_VET,COVMATUP,COVMATDOWN,COVMATID)
!!  DEALLOCATE(PRATIO)
!!  DEALLOCATE(LRATIOF)
!!  DEALLOCATE(LRATIOC)
!!  DEALLOCATE(ERRORCC)
!!  DEALLOCATE(ERRORNC)
!!  DEALLOCATE(ERRORCF)
!!  DEALLOCATE(ERRORNF)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL MPI_FINALIZE(IERROR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
331  FORMAT('TIME TO GENERATE RANDOM FIELD (s) =',F8.2,/,&
         'USER TIME           (s) =',F8.2,/,&
         'SYSTEM TIME         (s) =',F8.2,/)
332  FORMAT('================================',/,&
         'RANK................. =',I3,/,&
         'TAXA DE ACEITACAO (%) =',F8.2,/,&
         '================================') 
333  FORMAT('TIME TO SIMULATION  (s) =',F8.2,/,&
         'USER TIME           (s) =',F8.2,/,&
         'SYSTEM TIME         (s) =',F8.2,/)
334  FORMAT(I7,2X,E15.7)
1000 FORMAT('((----- REFAZENDO SIMULACAO FINE NO RANK: ',I2,'-----))')
1001 FORMAT('((----- REFAZENDO SIMULACAO COARSE NO RANK: ',I2,'-----))')
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
  FILEINPUT = './in/entraDREAM.in'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IN_FILE = 100+NK
  INQUIRE(FILE=FILEINPUT,EXIST=TFILE)
!
  OPEN(UNIT=IN_FILE,FILE=FILEINPUT,ACTION='READ',&
       STATUS='OLD',FORM='FORMATTED',IOSTAT=ISTAT)
!
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEINPUT
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
  READ(IN_FILE,100)NDIM
  NDIMR = NDIM*(NDIM+1)/2
  IF(NK.EQ.0)WRITE(*,108)NDIM,NDIMR
!
  ALLOCATE(XRMVN(NDIM))
  ALLOCATE(XRMVN_VET(NDIM,NPCHAIN))
  ALLOCATE(XRMVN_MEAN(NDIM))
  ALLOCATE(COVMATUP(NDIMR))
  ALLOCATE(COVMATDOWN(NDIMR))
  ALLOCATE(COVMATID(NDIMR))
  COVMATUP   = 0.0E0
  COVMATDOWN = 0.0E0
  COVMATID   = 0.0E0
  XRMVN      = 0.0E0
  XRMVN_VET  = 0.0E0
  XRMVN_MEAN = 0.0E0
  !! MATRIZ IDENTIDADE !!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NDIM
     DO J=I,NDIM
        LOC = J*(J-1)/2+I
        IF(I.EQ.J) COVMATID(LOC) = 1.0E+00
     END DO
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NDIR
     READ(IN_FILE,92)COVMATUP(I)
  END DO



  READ(IN_FILE,100)NPRIORR
  IF(NK.EQ.0)WRITE(*,98)NPRIORR
!
  DO I=1,NPRIORR
     READ(IN_FILE,200)NAMEPRIOR(I)
     FILEAM(I)=ADJUSTL(TRIM(NAMEPRIOR(I)))
     IF(NK.EQ.0)WRITE(*,206)I,FILEAM(I)
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
  IF(NK.EQ.0)WRITE(*,107)NSIMUL
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
        IF(NK.EQ.0)WRITE(*,*)'<<< DADOS DE CONCENTRACAO >>>'
     END IF
  END DO
!
  DO I=1,NPRIORR
     READ(IN_FILE,92)GERATIPO(I),NPROPOSAL(I),NLPRIOR(I),MKL(I),SIG(I)
     IF(NK.EQ.0)WRITE(*,151)I,GERATIPO(I),NPROPOSAL(I),MKL(I),SIG(I)
     IF(NPROPOSAL(I).EQ.2)THEN
        IF(NK.EQ.0)WRITE(*,152)
     END IF
     IF(NPROPOSAL(I).EQ.3)THEN
        IF(NK.EQ.0)WRITE(*,153)
     END IF
     IF(NPROPOSAL(I).EQ.4)THEN
        IF(NK.EQ.0)WRITE(*,154)
     END IF
     IF(NLPRIOR(I).EQ.0)THEN
        IF(NK.EQ.0)WRITE(*,600)
     ELSE
        IF(NK.EQ.0)WRITE(*,601)NLPRIOR(I)
        ALLOCATE(CLBDS(NLPRIOR(I)))
        DO J=1,NLPRIOR(I)
           READ(IN_FILE,200)NAME_AVET(I,J)
           IF(NK.EQ.0)WRITE(*,602)J,NAME_AVET(I,J)
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
  END DO
!
  IF(NPROPOSAL(1).EQ.3)THEN
     READ(IN_FILE,302)NSTART,NFREQ,NPCHAIN
     IF(NK.EQ.0)WRITE(*,303)NSTART,NFREQ,NPCHAIN
  END IF
!
  IF(NPROPOSAL(1).EQ.5)THEN
     READ(IN_FILE,920)NDELTA,NFREQ,CR_DREAM
     IF(NK.EQ.0)WRITE(*,921)NDELTA,NFREQ,CR_DREAM
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
     IF(NK.EQ.0)WRITE(*,202)FILEREF(I)
!
     READ(IN_FILE,200)FILEAM(I)
     FILEAM(I)=ADJUSTL(TRIM(FILEAM(I)))
     IF(NK.EQ.0)WRITE(*,203)FILEAM(I)
  END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  READ(IN_FILE,200)NAME_OUT
  NAME_OUT=ADJUSTL(TRIM(NAME_OUT))
  IF(NK.EQ.0)WRITE(*,205)NAME_OUT
!
  READ(IN_FILE,200)FILERROR
  FILERROR=ADJUSTL(TRIM(FILERROR))
  IF(NK.EQ.0)WRITE(*,204)FILERROR
!
  CALL RANDOM_SEED(size=NS)
  ALLOCATE(NSEED(NS))
  READ(IN_FILE,97)(NSEED(I),I=1,NS)
  CALL RANDOM_SEED(PUT=NSEED)
!  IF(NK.EQ.0)
  WRITE(*,96)NK,NS,(NSEED(I),I=1,NS)
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
92 FORMAT(E12.5)
920 FORMAT(2I10,E12.5)
921 FORMAT(/,'####### PARAMETROS DO METODO DREAM ########',/,&
           'NUMERO DE PARES USADOS NAS PORPOSTAS.... = ',I10,/,&
           'FREQUENCIA DE AJUSTE DO JUMP............ = ',I10,/,&
           'PROBABILIDADE DE TROCAS DOS ELEMENTOS DA   ',/,&
           'PROPOSTA (PARAMETRO CR DO MET. DREAM.... = ',E12.5,/,&
           '###########################################',/)

93 FORMAT(2I10)
94 FORMAT(8I10)
95 FORMAT('TIPO DE UPSCALING....................... = ',8I10)
96 FORMAT('############################',/,&
          'MAQUINA................: ',I3,/,&
          'NUMERO DE SEMENTES.....: ',I3,/,&
          'SEMENTES...............: ',20I10)
97 FORMAT(20I12)
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
           'PARAMETRO DO RANDOM WALK.................= ',E12.5)
152 FORMAT('---------------RANDON WALK---------------')
153 FORMAT('----------------METODO AM----------------')
154 FORMAT('----------------METODO DE----------------')
600 FORMAT('SEM MUDANCA NO COMPRIMENTO DE CORRELACAO')
601 FORMAT('LEITURA DAS',I2,' MATRIZES KL')
602 FORMAT('ARQUIVO n.',I2,' = ',A)
106 FORMAT('SALVAR CAMPOS REJEITADOS? (1)=SIM E (0)=NAO',I10)
107 FORMAT('DEFINE O SIMULADOR QUE SERA USADO....... = ',I10)
108 FORMAT('VARIABLE DIMENSION...................... = ',I10,/,&
           'COVARIANCE DIMENSION (REDUCED).......... = ',I10)
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
!
END SUBROUTINE READ_INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! INICIALIZACAO DO PROGRAMA !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIAL(NERROREAD,NK,NPROC)
  !
  USE MPI
  USE RANDOM
  USE VARIAVEIS
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: FOUT
  CHARACTER(LEN=5)   :: CHARAC
  INTEGER :: NERROREAD,I,ISTAT,NSINAL,NSINAL2,NSINALAM
  INTEGER :: OUT_FILE,OUT_FILEL,K,IOs,DIM,N,LOC,J,NK,NPROC
  INTEGER, DIMENSION(3) :: NVAL
  INTEGER, DIMENSION(NS):: NKSEED
  CHARACTER(LEN=4)      :: NUMB
  REAL,    DIMENSION(NS):: SEED
  LOGICAL :: TFILE,TFILEB
  INTEGER :: TAG,IERROR
  INTEGER :: STATUS(MPI_STATUS_SIZE)
  INTEGER :: ROOT
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
  SGN = 1D0
  FILEINI='./in/init_stat'
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)     
  FILEINI= TRIM(FILEINI)//('_RK')//TRIM(CHARAC)//TRIM('.in')
  OUT_FILE = 200+NK
  OUT_FILEL= 300+NK
!
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
!     CALL RANDOM_SEED(put=NSEED)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     FOUT = TRIM('./out/lambda_')//TRIM(NAME_OUT)//TRIM('.dat')
     FOUT = ADJUSTL(TRIM(FOUT))
     INQUIRE(FILE=FOUT,EXIST=TFILEB)
     IF(TFILEB)THEN
        OPEN(UNIT=OUT_FILEL,FILE=FOUT,STATUS='UNKNOWN',&
             ACTION='READ',IOSTAT=ISTAT)
        DO
           READ(OUT_FILEL,505,IOSTAT=IOs)NVAL(1),NVAL(2),NVAL(3),OLDLBD(NVAL(2))
           IF(IOs>0)THEN
              IF(NK.EQ.0)WRITE(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
              IF(NK.EQ.0)WRITE(*,*)"@    ALGO ERRADO NA LEITURA      @"
              IF(NK.EQ.0)WRITE(*,*)"@ PARA A CONTINUACAO DOS LAMBDAS @"
              IF(NK.EQ.0)WRITE(*,*)"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
              EXIT
           ELSE IF(IOs<0)THEN
              IF(NK.EQ.0)WRITE(*,*)"ULTIMO LAMBDA:",OLDLBD(NVAL(2))
              EXIT
           ELSE
              NUM_AVETOLD(NVAL(2))=NVAL(3)
           END IF
        END DO
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
!        CALL RANDOM_SEED(put=NSEED)
        DO I=1,NPROC-1
           CALL RANDOM_NUMBER(SEED)
           NKSEED = INT(1E+08*SEED)
           CALL MPI_SEND(NKSEED,NS,MPI_INTEGER,I,TAG,MPI_COMM_WORLD,IERROR)
        END DO
     ELSE
        CALL MPI_RECV(NSEED,NS,MPI_INTEGER,ROOT,TAG,MPI_COMM_WORLD,STATUS,IERROR)
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
  WRITE(*,FMT=502)NS,(NSEED(I),I=1,NS)
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
  IF(NSIMUL.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMULADOR >>>'
  END IF
  IF(NSIMUL.EQ.2)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMULADOR_BTMM >>>'
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
          GERATIPO(I).EQ.34)THEN
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
     NDIMR   = DIM*(DIM+1)/2
     IF(NSINALAM.GT.0)THEN
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%% METODO AM %%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DO VETOR ALEATORIO: ',DIM
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DA MATRIZ DE COVAR: ',NDIMR
     !
        ALLOCATE(XRMVN(DIM))
        ALLOCATE(XRMVN_VET(DIM,NPCHAIN))
        ALLOCATE(XRMVN_MEAN(DIM))
        ALLOCATE(COVMATUP(NDIMR))
        ALLOCATE(COVMATDOWN(NDIMR))
        ALLOCATE(COVMATID(NDIMR))
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
  IF(NPROPOSAL(1).EQ.4.OR.NPROPOSAL(1).EQ.5)THEN  
     DIM=0
     NSINALAM=0
     DO I=1,NPRIORR
        DIM = MAX(DIM,MKL(I))
        NSINALAM = 1
     END DO
     NDIMR   = DIM*(DIM+1)/2
     IF(NSINALAM.GT.0)THEN
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%% METODO DE %%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        IF(NK.EQ.0)WRITE(*,*)'DIMENSAO DA MATRIZ MATDE: ',DIM*NPROC
     !
        ALLOCATE(XRMVN(DIM))
        ALLOCATE(MATDE(DIM,NPROC))
        ALLOCATE(COVMATDOWN(NDIMR))
        ALLOCATE(COVMATID(NDIMR))
        ALLOCATE(XRMVN_MEAN(DIM))
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
  CALL READ_REF()
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
502 FORMAT('NUMERO DE SEMENTES: ',I2,/,'SEMENTES: ',20I12)
503 FORMAT(20I12)
500 FORMAT(4I10)
550 FORMAT(10E12.5)
501 FORMAT('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',/,&
           'N. DE TESTES MALHA GROSSA =',I10,/,&
           'N. DE TESTES MALHA FINA   =',I10,/,&
           'N. DE CAMPOS SELECIONADOS =',I10,/,&
           'N. DE REPETICOES NA CADEIA=',I10,/,&
           'ERRO ATUAL NA MALHA FINA  =',E12.5,/,&
           'ERRO ATUAL NA MALHA GROSSA=',E12.5,/,&
           '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
!
END SUBROUTINE INITIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  DOUBLE PRECISION   :: ERF,ERC
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
101 FORMAT(20I12)
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
  DOUBLE PRECISION     :: FDX,FDY,SIGG
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
  READ(IN_FILE,8000)INIF,NFILES
  IF(NK.EQ.0)WRITE(*,8001)INIF,NFILES
  NFILE=NFILES-INIF+1
  IF(NK.EQ.0)WRITE(*,8002)NFILE
!
  IF(GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33)THEN
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
       .OR.GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.34)THEN
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
  IF(GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33)THEN
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
  NERROREAD = 1d0
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
4000 FORMAT(20I12)
4001 FORMAT(/,'SEEDS =',20I12,/)
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
9001 FORMAT('CORRELATION LENGTH =',F12.4,/,&
          'VARIANCE   =',F12.4)
505 FORMAT('PROBLEMA NA DEFINICAO DE NCONDMAX',/,&
          'NUMERO DE PONTOS CONDICIONADOS SUPERA NCONDMAX',/,&
          'FAVOR ALTERAR NCONDMAX',/,&
          'NCOND=',I5,'NCONDMAX=',I5)
!
END SUBROUTINE READ_FIELD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! READ REFERENCE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READ_REF()
!
  USE VARIAVEIS
  IMPLICIT NONE
!
  INTEGER :: I,J,K,IN_FILE,ISTAT
  DOUBLE PRECISION, EXTERNAL :: MAXIMO,NORMAL2
!
  DO K=1,NDTYPE
!
     IN_FILE = 17
!
     WRITE(*,*)'##########################################'
     WRITE(*,*)'####### LENDO DADOS DE REFERENCIA ########'
     WRITE(*,*)TRIM(ADJUSTL(FILEREF(K)))
!
     OPEN(UNIT=IN_FILE,FILE=FILEREF(K),STATUS='UNKNOWN',&
          FORM='FORMATTED',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEREF(K)
        STOP
     END IF
!
     DO I=1,NDADOS(K)
        READ(IN_FILE,100)(REFDATA(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
!
     CLOSE(UNIT=IN_FILE)
!
!     MAXDATA(K)=MAXIMO(K)
     IF(NNORMA.EQ.1)THEN
        NORMA_L2_REF(K)=NORMAL2(K)
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
DOUBLE PRECISION FUNCTION MAXIMO(K)
  USE VARIAVEIS
  INTEGER, intent(in) :: K
  INTEGER :: I,J
!
  MAXIMO=-1E30
  DO I=1,NDADOS(K)
     DO J=2,NWELLS(K)+1
        MAXIMO = MAX(MAXIMO,REFDATA(K,J,I))
     END DO
  END DO
END FUNCTION MAXIMO
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
DOUBLE PRECISION FUNCTION NORMAL2(K)
  USE VARIAVEIS
  INTEGER, intent(in) :: K
  INTEGER :: I,J
  DOUBLE PRECISION :: INTEGRAL,AUX,AUXR
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
     ENDIF
!
!     write(*,100)k,SQRT(AUX)
     NORMAL2=SQRT(AUX)
!
100  FORMAT('DATA no.:',I3,2X,'NORMA L2 = ',E10.3)
!
   END FUNCTION NORMAL2
  
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
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEAM(K)
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
  CHARACTER(LEN=128)   :: FILEOUT,FILEVET,FILETHETA
  CHARACTER(LEN=4)     :: NUMB
  INTEGER              :: IN_FILE,ISTAT,I,J,NN,K,INIF,NFILES
  INTEGER, INTENT(IN)  :: NRK,NPROC
  INTEGER              :: NPROP
  INTEGER, INTENT(OUT) :: NERROREAD
  DOUBLE PRECISION     :: TESTE,AUX,BETA,VARUNI
  DOUBLE PRECISION,DIMENSION(NS) :: SEEDS
  INTEGER         ,DIMENSION(NS) :: N_SEED
  DOUBLE PRECISION     :: SIGK
!
  INIF=0D0
  NFILES=0D0
!
!  DO K=1,NPRIORR
  IF(CONT.EQ.0)THEN
     SIGK = 1.0E+00
  ELSE
     SIGK = SIG(K)
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NPROP = NPROPOSAL(K)
  IF(NPROPOSAL(K).EQ.3)THEN
     CALL AM_METHOD(K,NRK)
     NPROP = 2
  END IF
  IF(NPROPOSAL(K).EQ.4)THEN
     CALL DE_METHOD(K,NRK,NPROC,MKL(K))
     NPROP = 2
  END IF
  IF(NPROPOSAL(K).EQ.5)THEN
     CALL DREAM_METHOD(K,NRK,NPROC,MKL(K),3)
     NPROP = 2
  END IF
!! SEMENTES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RANDOM_NUMBER(SEEDS)
  DO I=1,NS
     N_SEED(I) = INT(SEEDS(I)*1.0E+08)
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
     NEWLBD(K)=DNINT(NEWLBD(K))
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
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILE_INFIELD(K)
     STOP
  END IF
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
     WRITE(IN_FILE,7000)FILETHETA
     WRITE(IN_FILE,1200)NPROP,SIGK
     WRITE(IN_FILE,3000)NCOND(K)
!
     DO I=0,NCOND(K)-1
        WRITE(IN_FILE,1000)VET(K,0,I),VET(K,1,I),VET(K,2,I)
     ENDDO
  ENDIF
!
  IF(GERATIPO(K).EQ.32)THEN
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
4000 FORMAT(20I12)
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
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.1)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  END IF
!
  IF(GERATIPO(K).EQ.2)THEN
     COMMAND=('cd ../gera_LABTRANGEO/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.3)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.31)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW1/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.34)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW2/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.32)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_RW3D/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.33)THEN
!       IF(NN.LE.1)THEN
!          CALL GERATHETA(MKL(K),GERATIPO(K),K)
!       ENDIF
     COMMAND=('cd ../gera_KL/FORTRAN_LBD/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.4)THEN
     COMMAND=('cd ../gera_LABTRANGEORW/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  IF(GERATIPO(K).EQ.5)THEN
     COMMAND=('cd ../gera_LABTRANGEORW2/; sh rodarKL.sh in')//&
          TRIM(ADJUSTL(NUMB))//(' > output.out')
  ENDIF
!
  WRITE(*,*)COMMAND
  CALL SYSTEM(COMMAND)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(GERATIPO(K).EQ.0.OR.GERATIPO(K).EQ.1.OR.&
       GERATIPO(K).EQ.3.OR.GERATIPO(K).EQ.31.OR.&
       GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.33.OR.&
       GERATIPO(K).EQ.34)THEN
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
!
     OPEN(UNIT=OUT_FILE,FILE=NOME,&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',NOME
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
     OPEN(UNIT=OUT_FILE,FILE=NOME,&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',NOME
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
  DOUBLE PRECISION   :: TESTE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ABERTURA DOS ARQUIVOS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  OUT_FILE = 1100+NK
  IF(TIPO.EQ.0) FILETHE='../gera_LABTRANGEORW3/out/xpsi'
  IF(TIPO.EQ.1) FILETHE='../gera_KL/FORTRAN/out/theta'
  IF(TIPO.EQ.3) FILETHE='../gera_KL/FORTRAN_RW/out/theta'
  IF(TIPO.EQ.31)FILETHE='../gera_KL/FORTRAN_RW1/out/theta'
  IF(TIPO.EQ.34)FILETHE='../gera_KL/FORTRAN_RW2/out/theta'
  IF(TIPO.EQ.32)FILETHE='../gera_KL/FORTRAN_RW3D/out/theta'
  IF(TIPO.EQ.33)FILETHE='../gera_KL/FORTRAN_LBD/out/theta'
!
  WRITE(NUMB,'(I4.3)')NK
  FOUT = TRIM(ADJUSTL(FILETHE))//TRIM(ADJUSTL(NUMB))//&
       ('.dat')
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
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
  IF(TIPO.EQ.31)FILETHE='../gera_KL/FORTRAN_RW1/out/thetanew'
  IF(TIPO.EQ.34)FILETHE='../gera_KL/FORTRAN_RW2/out/thetanew'
  IF(TIPO.EQ.32)FILETHE='../gera_KL/FORTRAN_RW3D/out/thetanew'
  IF(TIPO.EQ.33)FILETHE='../gera_KL/FORTRAN_LBD/out/thetanew'
!
  WRITE(NUMB,'(I4.3)')NK
  FOUT = TRIM(ADJUSTL(FILETHE))//TRIM(ADJUSTL(NUMB))//&
       ('.dat')
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
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
  double precision   :: xlx,xly,beta
  character(len=128) :: NOME,FILEINN
  character(len=128) :: fname
  character(len=4)   :: EXT,tipo
  character(len=5)   :: C
  DOUBLE PRECISION,DIMENSION(nnp,*) :: PERM
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
  DOUBLE PRECISION, dimension(nnp,*) :: prm
  integer :: nmx,nmy,i,j,NK,nnp
  DOUBLE PRECISION :: xlx,xly,xm,xmm,xv,XMAX,XMIN
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
  DOUBLE PRECISION :: INTEG,INTEGN
  DOUBLE PRECISION :: AUX,AUXN
  DOUBLE PRECISION :: MEAN,SIGMA
!
  IF(NPROPOSAL(1).EQ.2.OR.NPROPOSAL(1).EQ.3.OR.NPROPOSAL(1).EQ.4)THEN
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
             GERATIPO(K).EQ.32.OR.GERATIPO(K).EQ.34)THEN
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
     CALL SYSTEM(COMMAND)
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
  IF(NSIMUL.EQ.1)THEN
     COMMAND='cd ./simulador/; ./run > output.out'
  END IF
!
  IF(NSIMUL.EQ.2)THEN
     COMMAND='cd ./simuladorBTMM/; bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.3)THEN
     COMMAND='cd ./simul_comp/; bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.4)THEN
     COMMAND='cd ./co2injection/src/; bash rodarSimuladorPsh'
     WRITE(*,*)COMMAND
     CALL SYSTEM(COMMAND)
     COMMAND='cd ./co2injection/src/mrb_dados/; ./run >output.out'
  END IF
!
  IF(NSIMUL.EQ.5)THEN
     COMMAND='cd ./SIMULADOR_VISCOELASTICO/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > outputC')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  WRITE(*,200)NK,COMMAND
  CALL SYSTEM(COMMAND)
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

     FILEAMC = FILEAMC(2:LEN_TRIM(FILEAMC))
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
        READ(IN_FILE,100)(REFAMOS(K,J,I),J=1,NWELLS(K)+1)
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
  INTEGER :: K
  DOUBLE PRECISION,EXTERNAL :: NORMAL2
!
  DO K=1,NDTYPE
     ERRORCC(K)=NORMAL2(K)/NORMA_L2_REF(K)
     LRATIOC(K)=((0.5d0/SIGMA2C(K))*&
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
100 FORMAT('DATA no.:',I3,2X,'ERROR COARSE = ',E10.3)
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
  IF(NSIMUL.EQ.1)THEN
     COMMAND='cd ../simulador/; ./run > output.out'
  END IF
!
  IF(NSIMUL.EQ.2)THEN
     COMMAND='cd ../simuladorBTMM/; bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.3)THEN
     COMMAND='cd ../simul_comp/;  bash rodarSimuladorPsh'
  END IF
!
  IF(NSIMUL.EQ.4)THEN
     COMMAND='cd ../co2injection/src/; bash rodarSimuladorPsh'
     WRITE(*,*)COMMAND
     CALL SYSTEM(COMMAND)
     COMMAND='cd ../co2injection/src/mrb_dados/; ./run >output.out'
  END IF
!
  IF(NSIMUL.EQ.5)THEN
     COMMAND='cd ../SIMULADOR_VISCOELASTICO/; bash rodarSimulador.sh exp'
     COMMAND=TRIM(COMMAND)//TRIM(ADJUSTL(NUMB))
     COMMAND=TRIM(COMMAND)//('/ > output')//TRIM(ADJUSTL(NUMB))//('.out')
  END IF
!
  
  WRITE(*,*)COMMAND
  CALL SYSTEM(COMMAND)
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
  DOUBLE PRECISION,EXTERNAL :: NORMAL2
!
  WRITE(*,*)'####################################################'
  WRITE(*,*)'################### LIKELIHOOD #####################'
  DO K=1,NDTYPE
     ERRORCF(K)=NORMAL2(K)/NORMA_L2_REF(K)
     LRATIOF(K)=((0.5d0/SIGMA2F(K))*(ERRORNF(K)-ERRORCF(K)))
     write(*,100)k,ERRORCF(K)
  END DO
  WRITE(*,*)'####################################################'
!
  TERRORCF=1.D0
  TLRATIOF=0.D0
  DO K=1,NDTYPE
     TERRORCF = ERRORCF(K)*TERRORCF
     TLRATIOF = LRATIOF(K)+TLRATIOF
  END DO
  WRITE(*,101)TERRORCF
!
  TLRATIOF = EXP(TLRATIOF)
!
  RETURN
!
100 FORMAT('DATA no.:',I3,2X,'ERROR....= ',E10.3)
101 FORMAT('ERROR TOTAL............= ',E10.3)
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
  INTEGER          :: I,J,K,NN,TEST,NK
  REAL(4)          :: AUX
  DOUBLE PRECISION :: AUX2
  LOGICAL          :: FIRST
  INTEGER,DIMENSION(NS) :: N_SEED
!
  IF(SGN.EQ.1)THEN
     CALL RANDOM_SEED(get=N_SEED)
     SGN     =0D0
     FIRST   =.true.
     TEST    =RANDOM_BINOMIAL2(1,AUX,FIRST)
     CALL RANDOM_SEED(put=N_SEED)
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
     TEST    =RANDOM_BINOMIAL2(1,AUX,FIRST)
  ELSE
     FIRST   =.false.
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
     WRITE(*,100)AUX
     TEST=RANDOM_BINOMIAL2(1,AUX,FIRST)
     ACCEPT=TEST
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
  ELSE
     IF(NSTAGE.EQ.1)NCHAIN = NCHAIN+1
     WRITE(*,*)'############################'
     WRITE(*,*)'###### CAMPO REJEITADO #####'
     WRITE(*,*)'######## MALHA FINA ########'
     WRITE(*,*)'############################'
  END IF
!
  RETURN
!
100 FORMAT('(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)',/,&
           '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)',/,&
           'PROBABILITY............=',E10.3,/,&
           '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)',/,&
           '(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)(-)')
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
        IF(NSIMUL.EQ.1)THEN
           COMMAND='../simulador/fields/'
        END IF
        IF(NSIMUL.EQ.2)THEN
           COMMAND='../simuladorBTMM/exp/fields/'
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
             GERATIPO(K).EQ.34)THEN
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
              CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
        CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
                GERATIPO(K).EQ.34)THEN
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
           CALL SYSTEM(COMMAND)
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
           CALL SYSTEM(COMMAND)
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
  DOUBLE PRECISION   :: ER
  LOGICAL            :: TFILE
!
  OUT_FILE = 1600+NK
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)     
  FOUT = TRIM(FILERROR)//('_RK')//TRIM(CHARAC)//TRIM('.dat')
  INQUIRE(FILE=FOUT,EXIST=TFILE)
!  WRITE(*,*)'ARQUIVO EXISTE:',TFILE
  IF(TFILE)THEN
     OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='UNKNOWN',&
          ACTION='WRITE',IOSTAT=ISTAT,POSITION='APPEND')
  ELSE
     OPEN(UNIT=OUT_FILE,FILE=FOUT,&
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
        READ(IN_FILE,100)(REFAMOS(K,J,I),J=1,NWELLS(K)+1)
     ENDDO
  !
     CLOSE(UNIT=IN_FILE)
  END DO
!
100 FORMAT(30(e15.8,2x))
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
  CHARACTER(LEN=4)   :: CHARAC
!
!  WRITE(CHARAC,'(I4.3)')NK
  WRITE(CHARAC,'(I4)')NK
!
  OUT_FILE = 1800+NK
  FOUT = TRIM('./out/nchain_')//TRIM(NAME_OUT)//TRIM('_RK')//&
       TRIM(ADJUSTL(CHARAC))//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
!
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='UNKNOWN',&
       ACTION='READWRITE',IOSTAT=ISTAT,ACCESS='APPEND')
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FOUT
     STOP
  END IF
!
  WRITE(UNIT=OUT_FILE,FMT=100)CONT,NUMB
!
  CLOSE(OUT_FILE)      
!
100 FORMAT(I7,I7)
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
!
  WRITE(CHARAC,112)NK
  CHARAC=ADJUSTL(CHARAC)
  OUT_FILE = 1900+NK
  FOUT = TRIM('./out/lambda_')//TRIM(NAME_OUT)//('_RK')//&
       TRIM(CHARAC)//TRIM('.dat')
  FOUT = ADJUSTL(TRIM(FOUT))
!
  OPEN(UNIT=OUT_FILE,FILE=FOUT,STATUS='UNKNOWN',&
       ACTION='READWRITE',IOSTAT=ISTAT,ACCESS='APPEND')
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
  DOUBLE PRECISION   :: LBD
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
SUBROUTINE AM_METHOD(K,NK)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE VARIAVEIS, ONLY:  SIG,MKL,COVMATUP,XRMVN_VET,NDIMR
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,XRMVN_MEAN,XRMVN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,SIG,CONT,NINICIO
  USE VARIAVEIS, ONLY:  NSTART,NFREQ,NPCHAIN
!
  REAL                          :: XEPS,XCS,XCSININT,XAUX
  INTEGER                       :: K,I,J,N,DIM,NK
  INTEGER                       :: IERROR
  LOGICAL                       :: FIRST
  REAL,DIMENSION(MKL(K))        :: XVET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  NINICIO= NSTART
  DIM    = MKL(K)
  IERROR = 0D0
  FIRST  = .TRUE.
  XCS    = 2.50E-02
  XCSINIT= 1.00E+00
  XEPS   = ((2.38**2)/(REAL(DIM)))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SIG(K) = 0.0E+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(CONT.GE.NINICIO)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL LOADTHETAS(K,DIM,NPCHAIN,CONT,NK)
     XRMVN_MEAN = XRMVN_VET(1:DIM,NPCHAIN)
     COVMATUP   = COVAR(XRMVN_VET,DIM,NDIMR,NPCHAIN)
     COVMATUP   = XEPS*COVMATUP
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
  ELSE
     IF(CONT.GT.0)THEN
        CALL LOADTHETAS(K,DIM,CONT,CONT,NK)
        XRMVN_MEAN = XRMVN_VET(1:DIM,CONT-1)
        XAUX       = XCSINIT*XEPS
     ELSE
        XAUX       = 1.0E+00
     END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,XAUX*COVMATID,COVMATDOWN,&
          FIRST,XVET ,IERROR)
     THETAN(K,1:DIM) = XVET
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL GERA_AMTHETA(DIM,GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
!
END SUBROUTINE AM_METHOD
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
  D=TRIM('_')
  EXT='.dat'
  WRITE(C,111)K
  C=ADJUSTL(C)
  WRITE(CHARAC,111)NK
  CHARAC=ADJUSTL(CHARAC)
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
     WRITE(*,206)FILEINPUT
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
206 FORMAT('LEITURA DO ARQUIVO.........: ',A)
!
END SUBROUTINE LOADTHETAS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOADTHETASDE(K,DM,NVETCONT,NPROC)
!
  USE VARIAVEIS, ONLY: XRMVN,NAME_OUT,MATDE
  IMPLICIT NONE
!
  INTEGER                 :: IN_FILE,ISTAT,NCONT,II
  INTEGER                 :: NK,NPROC,I,J,K,N,DM
  INTEGER,DIMENSION(NPROC):: NVETCONT
  CHARACTER(LEN=128)      :: NOME,FILEINPUT
  CHARACTER(LEN=7)        :: CHARRANK
  CHARACTER(LEN=5)        :: C,D,M,CHARAC
  CHARACTER(LEN=4)        :: EXT
!
  IN_FILE = 2200
!  WRITE(*,205)NAME_OUT
  D=TRIM('_')
  EXT='.dat'
  WRITE(C,111)K
  C=ADJUSTL(C)
  DO NK=0,NPROC-1
     WRITE(M,111)NVETCONT(NK+1)-1
     M=ADJUSTL(M)
     WRITE(CHARAC,111)NK
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
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILEINPUT
        STOP
     END IF
!
     DO J=1,DM
        READ(IN_FILE,112)XRMVN(J)
!        WRITE(*,*)J,XRMVN(J)
     END DO
     MATDE(1:DM,NK+1) = XRMVN
     CLOSE(UNIT=IN_FILE)
  END DO
!
  RETURN 
!
112 FORMAT(E15.8)
111 FORMAT(I5)
205 FORMAT('NOME BASE DOS ARQUIVOS.... : ',A)
206 FORMAT('RANK:',I3,' -->  LEITURA DO ARQUIVO.........: ',A)
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
  DOUBLE PRECISION   :: TESTE
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
  IF(TIPO.EQ.31)FILETHE=('../gera_KL/FORTRAN_RW1/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.34)FILETHE=('../gera_KL/FORTRAN_RW2/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.32)FILETHE=('../gera_KL/FORTRAN_RW3D/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  IF(TIPO.EQ.33)FILETHE=('../gera_KL/FORTRAN_LBD/out/theta')//&
       TRIM(ADJUSTL(NUMB))//('.dat')
  OPEN(UNIT=OUT_FILE,FILE=FILETHE,STATUS='REPLACE',&
       ACTION='WRITE',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FILETHE
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
SUBROUTINE DREAM_METHOD(K,NK,NPROC,DIM,NCR)
!
  USE STATFUNCTIONS
  USE RANDOM
  USE MPI
  USE VARIAVEIS, ONLY:  SIG,MKL,XRMVN,MATDE,XRMVN_MEAN
  USE VARIAVEIS, ONLY:  GERATIPO,THETAN,CONT,VETCONT
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,NDIMR,SIGDE
  USE VARIAVEIS, ONLY:  CONTADORC
  USE VARIAVEIS, ONLY:  NDELTA,NFREQ,CR_DREAM
!
  REAL                  :: GAMMA_DREAM,AUX,STDX,AUXX,AUXXX
  REAL                  :: VUNIF,B,A,BSTAR,COEFF
  INTEGER               :: K,I,J,N,DIM,NK,NPROC,DESTART
  REAL,DIMENSION(DIM)   :: XVET,XSUM
  LOGICAL               :: FIRST
  INTEGER,DIMENSION(2)  :: R
  INTEGER,DIMENSION(DIM):: DVET
  INTEGER, EXTERNAL     :: INTMIN
  INTEGER               :: TEST,DLINHA
  INTEGER               :: NCR,IDCR,DELTA,ID
  REAL,DIMENSION(NCR)   :: PCR,CR,JVET
  REAL,DIMENSION(DIM)   :: DXVET
  INTEGER,DIMENSION(NCR):: NID
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DXVET   = 0.0
  JVET    = 0.0
  NID     = 0
  DVET    = 0
  FIRST   = .TRUE.
  B       = 1.0E-01
  A       = -B
  BSTAR   = 1.0E-06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(NDELTA.GT.NPROC)NDELTA=NPROC-1
  IF(NDELTA.EQ.0)    NDELTA=1
  DELTA = UNIFDISC(1,NDELTA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i=1,dim
     do j=1,nproc
        aux = random_normal()
        matde(i,j) = aux
     end do
  end do
  AUX = RANDOM_NORMAL()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO I=1,NCR
     CR(I) = REAL(I)/REAL(NCR)
     PCR(I)= REAL(1)/REAL(NCR)
  END DO
!
  IDCR = MULTINOM(NCR,PCR)
  CR_DREAM = CR(IDCR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DESTART = 10!INTMIN(VETCONT,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(DESTART.GT.0.OR.NPROC.GT.1)THEN
     SIG(K) = 0.0E+00
!     MATDE  = 0.0E+00
!
!     CALL LOADTHETASDE(K,DIM,VETCONT,NPROC)
!
     STDX     = STDV(MATDE,DIM,NPROC)
!
     XSUM = 0.0
     DO J=1,DELTA
        R    = CHOOSER(NK+1,NPROC)
        XSUM = XSUM + (MATDE(1:DIM,R(1))-MATDE(1:DIM,R(2)))
     END DO
!
     DLINHA = DIM
     ID   = 0
     DO J=1,DIM
        TEST   = RANDOM_BINOMIAL2(1,1.0-CR_DREAM,FIRST)
        DLINHA = DLINHA-TEST
        IF(TEST.EQ.1)THEN
           THETAN(K,J) = XRMVN(J)
        ELSE
           ID= ID+1
           DVET(ID) = J
        END IF
     END DO
!
     IF(DLINHA.EQ.0)THEN
        DLINHA = 1
        DVET(1) = UNIFDISC(1,DIM)
     END IF
!
     GAMMA_DREAM = 2.38/SQRT(2.0*REAL(DLINHA*DELTA))
     IF(MOD(CONTADORC,NFREQ).EQ.0)GAMMA_DREAM=1.0E+00
     WRITE(*,100)DELTA,DLINHA,GAMMA_DREAM
!
     DO I=1,DLINHA
        J = DVET(I)
        CALL RANDOM_NUMBER(VUNIF)
        VUNIF       = (B-A)*VUNIF+A
        COEFF       = (1.0-VUNIF)*GAMMA_DREAM
        DXVET(J)    =  COEFF*XSUM(J) + BSTAR*RANDOM_NORMAL()
        THETAN(K,J) = XRMVN(J) + DXVET(J)
     END DO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     CALL GERA_AMTHETA(MKL(K),GERATIPO(K),K,NK)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     JVET(IDCR) = JVET(IDCR) + SUMSQR(DXVET,STDX,DIM,DLINHA)
!     NID(IDCR)  = NID(IDCR)  + 1
!     WRITE(*,*)NK,IDCR,NID(IDCR),JVET(IDCR)
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  stop 'final de metodo dream'
  RETURN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT('NUMERO DE PARES USADOS........: ',I3,/,&
           'NUMERO DE DIMENSOES USADAS....: ',I3,/,&
           'LAMBDA METODO DREAM...........: ',F6.3)
!
END SUBROUTINE DREAM_METHOD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  USE VARIAVEIS, ONLY:  COVMATID,COVMATDOWN,NDIMR,SIGDE
  USE VARIAVEIS, ONLY:  CONTADORC
!
  REAL                 :: XEPS,XCS,XCSININT,XAUX
  INTEGER              :: K,I,J,N,DIM,NK,DESTART
  INTEGER              :: IERROR,TAG,NPROC
  INTEGER              :: ROOT,SOMACONT,NFREQ
  REAL,DIMENSION(DIM)  :: XVET
  LOGICAL              :: FIRST
  INTEGER,DIMENSION(2) :: R
  INTEGER, EXTERNAL    :: INTMIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ROOT   = NK
  FIRST  = .TRUE.
  NFREQ = 10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  XEPS   = 2.38/SQRT(2.0*REAL(DIM))
!  XEPS   = SIGDE
  IF(MOD(CONTADORC,NFREQ).EQ.0)XEPS=1.0E+00
  WRITE(*,100)XEPS
  IERROR = 0D0
  XCS    = 1.0E-03
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DESTART = INTMIN(VETCONT,NPROC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(DESTART.GT.0)THEN
     SIG(K) = 0.0E+00
     MATDE  = 0.0E+00
!
     CALL LOADTHETASDE(K,DIM,VETCONT,NPROC)
!
     R = CHOOSER(NK+1,NPROC)
!
     CALL RANDOM_MVNORM(DIM,XRMVN_MEAN,XCS*COVMATID,&
          COVMATDOWN,FIRST,XVET,IERROR)
     THETAN(K,1:DIM) = XRMVN + &
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     IF(TIPO.EQ.1)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN/in/ ')//&
             (' ../gera_KL/FORTRAN/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
     IF(TIPO.EQ.3)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW/in/ ')//&
             (' ../gera_KL/FORTRAN_RW/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
     IF(TIPO.EQ.31)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW1/in/ ')//&
             (' ../gera_KL/FORTRAN_RW1/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
     IF(TIPO.EQ.34)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW2/in/ ')//&
             (' ../gera_KL/FORTRAN_RW2/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
     IF(TIPO.EQ.32)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_RW3D/in/ ')//&
             (' ../gera_KL/FORTRAN_RW3D/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
     IF(TIPO.EQ.33)THEN
        COMMAND=('cp -r ../gera_KL/FORTRAN_LBD/in/ ')//&
             (' ../gera_KL/FORTRAN_LBD/in')//TRIM(ADJUSTL(NUMB))//('/')
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
  100 FORMAT('EXPEREMENT COPY: ',I3)
!
END SUBROUTINE COPY_DIRGER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
        CALL SYSTEM(COMMAND)
     END DO
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  RETURN
  100 FORMAT('EXPEREMENT COPY: ',I3)
!
END SUBROUTINE COPY_DIRUP
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
  IF(NUMSIMUL.EQ.2)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../simuladorBTMM/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../simuladorBTMM/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL SYSTEM(COMMAND)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ../simuladorBTMM/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ../simuladorBTMM/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.3)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../simul_comp/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../simul_comp/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL SYSTEM(COMMAND)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./simul_comp/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./simul_comp/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
     END IF
  END IF
!
  IF(NUMSIMUL.EQ.5)THEN
     WRITE(NUMB,'(I4.3)')NPR
     COMMAND=TRIM('cp -r ../SIMULADOR_VISCOELASTICO/exp')
     COMMAND=TRIM(COMMAND)//TRIM(' ../SIMULADOR_VISCOELASTICO/exp') &
          //TRIM(ADJUSTL(NUMB))
     WRITE(*,*)COMMAND
     CALL SYSTEM(COMMAND)
     IF(NS.EQ.2)THEN
        COMMAND=TRIM('cp -r ./SIMULADOR_VISCOELASTICO/exp')
        COMMAND=TRIM(COMMAND)//TRIM(' ./SIMULADOR_VISCOELASTICO/exp') &
             //TRIM(ADJUSTL(NUMB))
        WRITE(*,*)COMMAND
        CALL SYSTEM(COMMAND)
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
SUBROUTINE SAVEMPICONT(NCONT,NPROC,NK)
  !
  IMPLICIT NONE
  INTEGER            :: IERR,FID,NPROC,NK,I,ISTAT
  INTEGER            :: NCONT
  CHARACTER(LEN=128) :: FILENAME
  CHARACTER(LEN=5)   :: CHARAC
!
  FID = 3000+NK
  WRITE(CHARAC,'(I4.3)')NK
  CHARAC=TRIM(ADJUSTL(CHARAC))
  FILENAME = TRIM(ADJUSTL('out/vetcont'))//TRIM(ADJUSTL(CHARAC))//&
       TRIM(ADJUSTL('.dat'))
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 FORMAT(I7)
111 FORMAT(I5)
!
END SUBROUTINE SAVEMPICONT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SAVE MPI FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE READMPICONT(VET,NPROC,NK)
  !
  IMPLICIT NONE
  INTEGER                  :: IERR,FID,NPROC,NK
  INTEGER                  :: NCONT,I,ISTAT,K
  INTEGER,DIMENSION(NPROC) :: VET
  CHARACTER(LEN=128)       :: FILENAME
  CHARACTER(LEN=5)         :: CHARAC
!
  FID = 4000+NK
  DO K=1,NPROC
     WRITE(CHARAC,'(I4.3)')K-1
     CHARAC=TRIM(ADJUSTL(CHARAC))
     FILENAME = TRIM(ADJUSTL('out/vetcont'))//TRIM(ADJUSTL(CHARAC))//&
          TRIM(ADJUSTL('.dat'))
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
  END DO
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
  REAL,INTENT(IN)   :: ST
!
  SUM = 0.0
  DO I=1,DM
     AUX = DX(I)/ST
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
