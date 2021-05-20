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
    external                     :: dpotrf,dtrtrs
    integer                      :: kr,j,loc
    integer, intent(in)          :: k
    real,intent(in),dimension(k) :: x, x_mean
    real,dimension(k)            :: xminusmean
    real,intent(in out),dimension(k*(k+1)/2) :: sig
    real,dimension(k,k)                      :: sigma
    real    :: rtdet_sigma
    real    :: c, inside_exp
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
!
  end function mv_normal_pdf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FUNCAO PARA CALCULO DO PRODUTO INTERNO !!!!!!!!!!!!!!!  
  REAL FUNCTION MRBDOT(X,Y,K)
    IMPLICIT NONE
    INTEGER             :: I
    INTEGER,INTENT(IN)  :: K
    REAL   ,DIMENSION(K),INTENT(IN):: X
    REAL   ,DIMENSION(K),INTENT(IN):: Y
    REAL                :: PRODI
!
    PRODI = 0.0
    DO I=1,K
       PRODI = PRODI + X(I)*Y(I)
    END DO
    MRBDOT = PRODI
  END FUNCTION MRBDOT
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
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: NK,NP
    INTEGER               :: R1,R2
    INTEGER, DIMENSION(2) :: CHOOSER
    !
    R1 = NK
    R2 = NK
    DO WHILE(R1.EQ.NK.OR.R2.EQ.NK.OR.R1.EQ.R2)
       R2 = UNIDRND(NP)
       R1 = UNIDRND(NP)
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
    UNIDRND = FLOOR(REAL(N)*AUX)+1
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
    REAL    :: AUX,FA,FB,P
!
    NB = MAX(A,B)
    NA = MIN(A,B)
    N  = NB-NA+1
    P  = 1.0/DBLE(N)
    ID = 1
!    CALL RANDOM_SEED()
    CALL RANDOM_NUMBER(AUX)
    FA = 0.0
    FB = 0.0
    DO I=1,N-1
       FA = FA+P
       FB = FA+P
       IF(AUX.GT.FA.AND.AUX.LE.FB)ID=I+1
    END DO
    UNIFDISC = ID + NA - 1
    RETURN
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
!    CALL RANDOM_SEED()
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!! FUNCAO PARA CALCULO DO DESVIO PADRAO !!!!!!!!!!!!!!!!!
  FUNCTION STDV(VET,D,TAM)
    INTEGER, INTENT(IN)   :: D,TAM
    REAL, DIMENSION(D,TAM),INTENT(IN) :: VET
    INTEGER               :: I,J,K
    REAL, DIMENSION(D)    :: XMEAN,SOMA,SOMASQ,STDV
!
    SOMA  = 0.0D0
    SOMASQ= 0.0D0
    DO J=1,D
       DO K=1,TAM
          SOMA(J) = SOMA(J) + VET(J,K)
       END DO
    END DO
    SOMA  = SOMA/DFLOAT(TAM)
    DO J=1,D
       DO K=1,TAM
          SOMASQ(J) = SOMASQ(J) + (VET(J,K)-SOMA(J))**2
       END DO
    END DO
    SOMASQ= SOMASQ/DFLOAT(TAM-1)
    STDV = SQRT(SOMASQ)
!
  END FUNCTION STDV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE STATFUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
