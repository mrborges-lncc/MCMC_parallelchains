!
!**** new module **********************************************************************
!
      module mGeomecanica

      use mGlobaisEscalares, only: novaMalha

!ESCALARES
      REAL(8) :: KGYNG, RHOYNG ! YOUNG's GEOMETRIC MEAN, STRENGHT
      INTEGER :: NDAUX, NDT
      REAL*8  :: TMANDEL, GEOTIME
      INTEGER :: NROWB, NINTD, NESD, NED, NED2, NSTR
      INTEGER :: NEESQ2, NEESQ, NEE, NEE2,  IOPT, IBBAR, NLOOPS

      INTEGER, ALLOCATABLE :: IDDIS(:,:)
      REAL(8), ALLOCATABLE :: SIGMAT(:),SIGMA0(:),DIVU(:),DIVU0(:)
      REAL(8), ALLOCATABLE :: DIS(:,:), STRSS(:,:,:), VDP(:,:), DTSIGM(:),BALANC(:)
! CREEP
!      INTEGER, ALLOCATABLE, DIMENSION(:)    :: MECLAW
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: DTRL, AVCREP, AVSTRS, STRSS_LIN
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HMTTG, ECREEP, AUXM
      INTEGER :: NROWB2, NWTNITER,NELDOMO
      REAL*8  :: RESMAX, RESIDUAL
!
      real*8,  allocatable :: fDesloc(:,:)
!
      CHARACTER(len=128) :: YNG_IN
      real(8) :: kg,kgphi    ! media geometrica
      real(8) :: rho,rhophi ! coeficiente da variancia (strenght)
!
      contains
!
!**** new **********************************************************************
!
      subroutine montarSistEqAlgGeo(tipo)
!
      use mAlgMatricial,     only: idDesloc, ALHSD, BRHSD, CLHSD, NEQD, LMD, IDIAGD, NALHSD, LOAD, FTOD,FACTOR
      use mMalha,            only: conecNodaisElem, listaDosElemsPorNo, x, numnp, nen, nsd
      use mGlobaisEscalares, only: NDOFD, NLVECTD, creep, optSolver
      use mHidroDinamicaRT,  only: pressaoElem
!
      implicit none
!
      character(len=8) :: tipo
! 
!.... ALLOCATE MEMORY FOR GLOBAL EQUATION SYSTEM 
! 
      print*, "tipo=", tipo
!
      IF(TIPO=='bbarmtrx') THEN
!
      if(.not.allocated(ALHSD)) allocate(ALHSD(NALHSD))
      if(.not.allocated(BRHSD)) allocate(BRHSD(NEQD))
!
!.... CLEAR LEFT AND RIGHT HAND SIDE FOR THE ELASTIC  PROBLEM 
! 
      ALHSD=0.0d00
      BRHSD=0.0d00
!
!.... ACCOUNT THE NODAL FORCES IN THE R.H.S. 
!
         IF (NLVECTD.GT.0) &
            CALL LOAD(idDesloc,fDesloc,brhsd,NDOFD,NUMNP,NLVECTD)
!
         IF (NLVECTD.GT.0) &
            CALL FTOD(idDesloc,DIS,fDesloc,NDOFD,NUMNP,NLVECTD)

!
!.... MOUNT THE L.H.S AT ELEMENT LEVEL 
!....       (S U B R O U T I N E:     W E A K -- --  B I O T 2) 
! 
       if(creep.eqv..false.) then
         CALL BBARMTRX_LIN(x, conecNodaisElem, alhsd, brhsd, idiagD, lmD)
       else
         if(nsd==2)CALL BBARMTRX   (x, pressaoElem, conecNodaisElem, alhsd, brhsd, idiagD, lmD)
         if(nsd==3)CALL BBARMTRX_3D(x, pressaoElem, conecNodaisElem, alhsd, brhsd, idiagD, lmD)
       endif
!
       return
!
      ENDIF
!
      IF(TIPO=='vectrscr'.and.creep.eqv..false.) THEN
!
!.... MOUNT VECTOR SOURCE FOR ELASTIC PROBLEM 
!
         CALL VECTRSCR_LIN(x, conecNodaisElem,pressaoElem, brhsd, lmD)
!
         RETURN
!
      ENDIF
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine BBARMTRX(x, p, conecNodaisElem,alhsd, brhsd, idiagD, lmD)
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE 
!        STOKE'S DISPLACEMENT  ELEMENT AND 
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX 
!        AND RIGHT-HAND SIDE VECTOR 
 
      use mGlobaisEscalares, only: ndofD, ndofP,  nrowsh,nnp
      use mGlobaisArranjos,  only: c
      use mPropGeoFisica,    only: YOUNG, POISBTTM, YNGBOTTM, BULKSOLID, REGION, MECLAW
      use mAlgMatricial,     only: rowdot, coldot, addrhs, addlhs, neqD, nalhsD
      use mFuncoesDeForma,   only: shgq, shlq
      use mMalha,            only: local, numnp, multab
      use mMalha,            only: nsd, numel, numelReserv, nen
      use mSolucoesExternas, only: addlhsCRS, ApGeo, AiGeo
!
      implicit none
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
      real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqD)
      integer, intent(in)    :: idiagD(*)
      integer, intent(in)    :: lmD(ned2,nen,numel)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!  
      LOGICAL DIAG,QUAD,ZERODL,LSYM
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ, BBARI, QIXIBBAR
      REAL(8), DIMENSION(NROWB,NROWB)   :: QIXI
      REAL(8), DIMENSION(NROWB)         :: TENSAO, PRESSURE
!
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL(8) :: C1
      INTEGER :: NEL, I, J, L, K
!
      integer :: tid, omp_get_thread_num,omp_get_num_threads
      integer :: numPrimeiroElemento, numUltimoElemento, numThreads, inicioSol, fimSol
!
      REAL(8) :: BULKROCK, BIOTCOEF 
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!      CONSISTENT MATRIX
! 
      DIAG       = .FALSE. 
      TENSAO     = 0.0D0
      tid        = 1
      numThreads = 1

! !$OMP PARALLEL FIRSTPRIVATE(tid) &
! !$OMP PRIVATE (numPrimeiroElemento, numUltimoElemento, inicioSol,fimSol) &
! !$OMP PRIVATE (I,J,L,K,NEL,BBARI,BBARJ,QIXIBBAR,QIXI,PRESSURE,XL,DISL,TENSAO,C1,IOPT,QUAD,LSYM) &
! !$OMP PRIVATE (R,SHGD,SHGBR,DETD) &
! !$OMP REDUCTION(+:ELEFFMD,ELRESFD) !,BRHSD,ALHSD)
! 
! 
! #ifdef withOMP
!        tid=tid+omp_get_thread_num()
!        numThreads=omp_get_num_threads()
! #endif
!
       if(tid==1) print*, "Em  BBARMTRX, numThreads=",numThreads

   
      numPrimeiroElemento = 1
      numUltimoElemento   = numel
      call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, numThreads, tid-1, inicioSol, fimSol)
!

! !$OMP PARALLEL DO ORDERED SCHEDULE(DYNAMIC)
      DO 500 NEL=inicioSol,fimSol
!      DO 500 NEL=1,numel

      BULKROCK = BULK(YNGBOTTM,POISBTTM)
!
      BIOTCOEF = 1.0D0 - BULKROCK/BULKSOLID
!
      BBARI    = 0.0D0
      BBARJ    = 0.0D0
      QIXIBBAR = 0.0D0
      QIXI     = 0.0D0
!
!.... SETUP ELEMENT PRESSURE VECTOR
!
      IF (REGION(NEL).EQ.1) THEN
         PRESSURE(1)=BIOTCOEF*P(1,NEL)
         PRESSURE(2)=BIOTCOEF*P(1,NEL)
         PRESSURE(3)=0.0D0
         PRESSURE(4)=BIOTCOEF*P(1,NEL) 
       ELSE
         PRESSURE = 0.0D0
      ENDIF
!

!.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY 
! 
      ELEFFMD=0.0D0 
      ELRESFD=0.0D0 
! 
!.... LOCALIZE COORDINATES AND DIRICHLET B.C. 
! 
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD) 
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2) 
!
      QUAD = .TRUE. 
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.  
      CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN) 
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
        DO 100 L=1,NINTD
           R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
           DETD(L) = DETD(L)*R(L) 
 100    CONTINUE
      ENDIF
!
!.... MOUNT STIFFNESS MATRIX 
!
!.... CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
      CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... LOOP OVER INTEGRATIONN POINTS
!
      DO 400 L=1,NINTD 
!
!... ... SETUP TANGENT MATRIX QIXI: ORDER 4X4 FOR MULTIPLICATION
!
         CALL TANG2QX(HMTTG(NEL,L,1:16),QIXI)
!
!...... SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
!
         C1=DETD(L)*WD(L)
!
!.... ...UPLOAD B-BAR MATRIX AT NODE J
!
         DO 200 J=1,NEN
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),  &
     &           SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MULTIPLY QIXI*BBARJ ===>QIXIBBAR
!
            QIXIBBAR=matmul(QIXI,BBARJ)
!             CALL MULTAB(QIXI,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
!
!.... .... UPLOAD B-BAR MATRIX AT NODE I
!
            DO 200 I=1,NEN
               CALL SETBB(BBARI,SHGD(1:NROWSH,I,L),  &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .......  MOUNT ELEMNT STIFFNESS NODAL MATRIX:
!.... .......  K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
!
               ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1)  &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J)      &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4)*C1
!
               ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1)      &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J)          &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4)*C1
 200        CONTINUE
      BBARI=0.d0
!
      DO 350 I=1,NEN
         CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
         DO 340 K=1,NROWB
            TENSAO(K) = STRSS(NEL,L,K)-PRESSURE(K)

 340     CONTINUE
         ELRESFD(NED2*I-1)= ELRESFD(NED2*I-1)  &
     &                      - COLDOT(BBARI(1:4,1),TENSAO,4)*C1
         ELRESFD(NED2*I)  = ELRESFD(NED2*I)    &
     &                      - COLDOT(BBARI(1:4,2),TENSAO,4)*C1
 350  CONTINUE
!
 400  CONTINUE
!  
!.... COMPUTATION OF DIRICHLET B.C. CONTRIBUTION 
! 
     CALL ZTEST(DISL,NEE2,ZERODL) 
!
     IF(.NOT.ZERODL)CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL) 
!
!.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL 
!....     LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR 
! 
      LSYM=.TRUE.

! !$OMP ORDERED

#ifdef withcrs
      call addlhsCRS(ALHSD,ELEFFMD,LMD(1,1,nel),ApGeo,AiGeo,NEE2) 
#else
     CALL ADDLHS(ALHSD,ELEFFMD,IDIAGD,LMD(1,1,NEL),NEE2, &
     &            DIAG,LSYM) 
#endif
!
      CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2) 

! !$OMP END ORDERED
       
 500   CONTINUE 
! !OMP END DO
! 
! !$OMP END PARALLEL
!
      RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
 2222  FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X   &
     &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                    &
     &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                    &
     &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                    &
     &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine BBARMTRX_3D(x, p, conecNodaisElem,alhsd, brhsd, idiagD, lmD)
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE 
!        STOKE'S DISPLACEMENT  ELEMENT AND 
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX 
!        AND RIGHT-HAND SIDE VECTOR 
 
      use mGlobaisEscalares, only: ndofD, ndofP,  nrowsh,nnp
      use mGlobaisArranjos,  only: c
      use mPropGeoFisica,    only: YOUNG, POISBTTM, YNGBOTTM, BULKSOLID, REGION, MECLAW
      use mAlgMatricial,     only: rowdot, coldot, addrhs, addlhs, neqD, nalhsD
      use mFuncoesDeForma,   only: SHG3D, SHLQ3D
      use mMalha,            only: local, numnp, multab
      use mMalha,            only: nsd, numel, numelReserv, nen
      use mSolucoesExternas, only: addlhsCRS, ApGeo, AiGeo
!
      implicit none
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
      real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqD)
      integer, intent(in)    :: idiagD(*)
      integer, intent(in)    :: lmD(ned2,nen,numel)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!  
      LOGICAL DIAG,QUAD,ZERODL,LSYM
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ, BBARI, QIXIBBAR
      REAL(8), DIMENSION(NROWB,NROWB)   :: QIXI
      REAL(8), DIMENSION(NROWB)         :: TENSAO, PRESSURE
!
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL(8) :: C1
      INTEGER :: NEL, I, J, L, K
!
      integer :: tid, omp_get_thread_num,omp_get_num_threads
      integer :: numPrimeiroElemento, numUltimoElemento, numThreads, inicioSol, fimSol
!
      REAL(8) :: BULKROCK, BIOTCOEF 

!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
!       CALL SHLQ(SHLD,WD,NINTD,NEN)
      CALL SHLQ3D  (SHLD,WD,NINTD,NEN)
!
!      CONSISTENT MATRIX
! 
      DIAG       = .FALSE. 
      TENSAO     = 0.0D0
      tid        = 1
      numThreads = 1

      DO 500 NEL=1,numel

      BULKROCK = BULK(YNGBOTTM,POISBTTM)
!
      BIOTCOEF = 1.0D0 - BULKROCK/BULKSOLID
!
      BBARI    = 0.0D0
      BBARJ    = 0.0D0
      QIXIBBAR = 0.0D0
      QIXI     = 0.0D0
!
!.... SETUP ELEMENT PRESSURE VECTOR
!
      IF (REGION(NEL).EQ.1) THEN
         PRESSURE(1)=BIOTCOEF*P(1,NEL)
         PRESSURE(2)=BIOTCOEF*P(1,NEL)
         PRESSURE(3)=0.0D0
         PRESSURE(4)=BIOTCOEF*P(1,NEL) 
       ELSE
         PRESSURE = 0.0D0
      ENDIF
!

!.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY 
! 
      ELEFFMD=0.0D0 
      ELRESFD=0.0D0 
! 
!.... LOCALIZE COORDINATES AND DIRICHLET B.C. 
! 
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD) 
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2) 
!
      QUAD = .TRUE. 
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.  
!       CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
       CALL SHG3D(XL,DETD,SHLD,SHGD,NINTD,NEL,NEN) 
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
        DO 100 L=1,NINTD
           R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
           DETD(L) = DETD(L)*R(L) 
 100    CONTINUE
      ENDIF
!
!.... MOUNT STIFFNESS MATRIX 
!
!.... CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
      CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... LOOP OVER INTEGRATIONN POINTS
!
      DO 400 L=1,NINTD 
!
!... ... SETUP TANGENT MATRIX QIXI: ORDER 4X4 FOR MULTIPLICATION
!
         CALL TANG2QX(HMTTG(NEL,L,1:16),QIXI)
!
!...... SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
!
         C1=DETD(L)*WD(L)
!
!.... ...UPLOAD B-BAR MATRIX AT NODE J
!
         DO 200 J=1,NEN
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),  &
     &           SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MULTIPLY QIXI*BBARJ ===>QIXIBBAR
!
            QIXIBBAR=matmul(QIXI,BBARJ)
!             CALL MULTAB(QIXI,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
!
!.... .... UPLOAD B-BAR MATRIX AT NODE I
!
            DO 200 I=1,NEN
               CALL SETBB(BBARI,SHGD(1:NROWSH,I,L),  &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .......  MOUNT ELEMNT STIFFNESS NODAL MATRIX:
!.... .......  K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
!
               ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1)  &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J)      &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4)*C1
!
               ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1)      &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J)          &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4)*C1
 200        CONTINUE
      BBARI=0.d0
!
      DO 350 I=1,NEN
         CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
         DO 340 K=1,NROWB
            TENSAO(K) = STRSS(NEL,L,K)-PRESSURE(K)

 340     CONTINUE
         ELRESFD(NED2*I-1)= ELRESFD(NED2*I-1)  &
     &                      - COLDOT(BBARI(1:4,1),TENSAO,4)*C1
         ELRESFD(NED2*I)  = ELRESFD(NED2*I)    &
     &                      - COLDOT(BBARI(1:4,2),TENSAO,4)*C1
 350  CONTINUE
!
 400  CONTINUE
!  
!.... COMPUTATION OF DIRICHLET B.C. CONTRIBUTION 
! 
     CALL ZTEST(DISL,NEE2,ZERODL) 
!
     IF(.NOT.ZERODL)CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL) 
!
!.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL 
!....     LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR 
! 
      LSYM=.TRUE.


#ifdef withcrs
      call addlhsCRS(ALHSD,ELEFFMD,LMD(1,1,nel),ApGeo,AiGeo,NEE2) 
#else
     CALL ADDLHS(ALHSD,ELEFFMD,IDIAGD,LMD(1,1,NEL),NEE2, &
     &            DIAG,LSYM) 
#endif
!
      CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2) 
!       
 500   CONTINUE 
!
      RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
 2222  FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X   &
     &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                    &
     &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                    &
     &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                    &
     &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
!
      END SUBROUTINE BBARMTRX_3D
!

!
!**** new **********************************************************************
!
      subroutine BBARMTRX_LIN(x, conecNodaisElem,alhsd, brhsd, idiagD, lmD)
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE 
!        STOKE'S DISPLACEMENT  ELEMENT AND 
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX 
!        AND RIGHT-HAND SIDE VECTOR 

      use mGlobaisEscalares, only: ndofD, nrowsh
      use mGlobaisArranjos,  only: c
      use mPropGeoFisica,    only: YOUNG, POISBTTM, POISMIDL, POISTOP, GEOFORM
      use mAlgMatricial,     only: rowdot, coldot, addrhs, addlhs, neqd,nalhsD
      use mFuncoesDeForma,   only: shgq, shlq
      use mMalha,            only: local, numnp, multab
      use mMalha,            only: nsd, numnp, numel, nen, numelReserv
      use mSolucoesExternas, only: addlhsCRS, ApGeo, AiGeo
!
      implicit none
!
      real*8,  intent(in)    :: x(nsd,numnp)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
      real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqd)
      integer, intent(in)    :: idiagD(*)
      integer, intent(in)    :: lmD(ned2,nen,numel)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!  
      LOGICAL DIAG,QUAD,ZERODL,LSYM
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL*8 :: BBARJ(NROWB,NESD), BBARI(NROWB,NESD),&
     &          QIXI(NROWB,NROWB), QIXIBBAR(NROWB,NESD), C1, POISSON
      INTEGER :: NEL, I, J, L, II, JJ
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: CBBAR(NROWB, NROWB)
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!.... CONSISTENT MATRIX
! 
      DIAG = .FALSE.

!       do nel=1, numel
!       write(852,*) nel,young(nel)
!       enddo
! 
!      stop
!
      DO 500 NEL=1,NUMEL
! 
!.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY 
! 
      ELEFFMD=0.D0
      ELRESFD=0.D0
! 
!.... LOCALIZE COORDINATES AND DIRICHLET B.C. 
! 
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
      QUAD = .TRUE.
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
! 
       CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
         DO 100 L=1,NINTD
            R(L)   = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
            DETD(L) = DETD(L)*R(L)
 100    CONTINUE
      ENDIF
!
      POISSON = POISBTTM

      IF ((NEL.GT.numelReserv).AND.(NEL.LE.NELDOMO)) POISSON = POISMIDL
      IF (NEL.GT.NELDOMO) POISSON = POISTOP
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
      CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
!
!.... FORM STIFFNESS MATRIX 
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... .. FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
          CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
 !
!.... .. LOOP OVER INTEGRATIONN POINTS
!
         DO 400  L=1,NINTD

!
!.... ...SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
!
          C1=DETD(L)*WD(L)
!
!.... ... UPLOAD B-BAR MATRIX AT NODE J
!
          DO 200 J=1,NEN
!
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L), &
     &            SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MULTIPLY CBBAR*BBARJ ===>QIXIBBAR
!
            CALL MULTAB(CBBAR,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
!
!.... .... UPLOAD B-BAR MATRIX AT NODE I
!
            DO 200 I=1,NEN
!
            CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), &
     &            SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MOUNT ELEMNT STIFFNESS NODAL MATRIX:
!.... ...        K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
!
             ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1) &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4)*C1
!
             ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J) &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4)*C1
!
             ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1) &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4)*C1
!
             ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J)  &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4)*C1
!
 200      CONTINUE
 400  CONTINUE

!  
!      COMPUTATION OF DIRICHLET B.C. CONTRIBUTION 
! 
      CALL ZTEST(DISL,NEE2,ZERODL)
!            write(852,*)disl(1:ned2,1:nen), NEE2
!
      IF(.NOT.ZERODL) &
     &CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL)
!            write(852,*)nel, ELRESFD(1:nee2)
!
!.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL 
!        LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR 
! 
      LSYM=.TRUE.
!
#ifdef withcrs
     call addlhsCRS(ALHSD,ELEFFMD,LMD(1,1,nel),ApGeo,AiGeo,NEE2) 
#else
      CALL ADDLHS(ALHSD,ELEFFMD,IDIAGD,LMD(1,1,NEL),NEE2, &
     &            DIAG,LSYM) 
#endif

!
      CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2)

        DO 450 II=1,NEE2
        DO 450 JJ=1,NEE2
            AUXM(NEL,II,JJ)=ELEFFMD(II,JJ)
!             print*, auxm(nel,ii,jj)
 450   CONTINUE
!
 500   CONTINUE
!
!        do nel=1, nalhsD
!       write(541,*), nel, ALHSD(nel) !brhsd(nel)
!        enddo

!        stop

      RETURN
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE SETUPC(CBBAR,YOUNG,POISSON,NROWB,IOPT)
! 
!     PROGRAM TO SETUP ELASTICITY TENSOR 
! 
      IMPLICIT REAL*8 (A-H,O-Z) 
! 
      REAL*8 :: YOUNG,POISSON
      INTEGER :: NROWB,IOPT
      REAL(8) :: CBBAR(NROWB, NROWB)
!
      REAL(8) :: AMU2,ALAM

!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
!
!..... SET MATERIAL PARAMETERS FOR OUT-OF-PLANE COMPONENTS
!
!      C(5,M) == YOUNG   e   C(6,M)==POISSON
!          AMU2 = YOUNG/(ONE+POISSON)
!          ALAM = AMU2*POISSON/(ONE-TWO*POISSON)
!
          CBBAR = 0.0
          AMU2 = YOUNG/(1.0D0+POISSON)
          ALAM = AMU2*POISSON/(1.0D0-2.0D0*POISSON)
!
!..... COLUMN MATRIX D3 
!
          CBBAR(1,4) = ALAM
          CBBAR(2,4) = ALAM
          CBBAR(3,4) = 0.0D0
          CBBAR(4,4) = ALAM + AMU2
!
!..... TRANSPOSE OF D3
!
          CBBAR(4,1) = CBBAR(1,4)
          CBBAR(4,2) = CBBAR(2,4)
          CBBAR(4,3) = CBBAR(3,4)
!
!..... SET MATERIAL PARAMETERS FOR IN-PLANE COMPONENTS
!
!..... .. SET STRESS PLANE CONDITION IOPT.EQ.0
!
          IF (IOPT.EQ.0) ALAM = ALAM*AMU2/(ALAM + AMU2)
!
!..... .. MATRIX D_33: UP TRIANGULAR
!
          CBBAR(1,1) = ALAM + AMU2
          CBBAR(1,2) = ALAM
          CBBAR(2,2) = CBBAR(1,1)
          CBBAR(1,3) = 0.0D0
          CBBAR(2,3) = 0.0D0
          CBBAR(3,3) = 0.5D0*AMU2
!
!..... .. MATRIX D_33: LOW TRIANGULAR
! 
          CBBAR(2,1) = CBBAR(1,2)
          CBBAR(3,1) = CBBAR(1,3)
          CBBAR(3,2) = CBBAR(2,3)
!
 100   CONTINUE
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
       END SUBROUTINE
!
!**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 760 *************** 
!
      SUBROUTINE XMEANSH(SHGBAR,W,DET,R,SHG,NEN,NINT,IOPT,NESD,NROWSH)
!
!.... PROGRAM TO CALCULATE MEAN VALUES OF SHAPE FUNCTION
!        GLOBAL DERIVATIVES FOR B-BAR METHOD
!
!        NOTE: IF IOPT.EQ.2, DET(L)=DET(L)*R(L) UPON ENTRY
!
      use mAlgMatricial, only: coldot
!      use mGlobaisEscalares, only: one

      IMPLICIT REAL*8 (A-H,O-Z) 
!
!.... DEACTIVATE ABOVE CARD(S) FOR SINGLE-PRECISION OPERATION
!
      INTEGER :: NEN, NROWSH, NINT, IOPT, NESD
      REAL(8) :: SHGBAR(3,NEN),W(*),DET(*),R(*),SHG(NROWSH,NEN,*)
!
      REAL(8) :: VOLINV, TEMP1, TEMP2
      INTEGER :: L, I, J
!
      CALL CLEAR(SHGBAR,3*NEN)
!
      VOLINV = 1.0D0/COLDOT(W,DET,NINT)
!
      DO 300 L=1,NINT
         TEMP1 = W(L)*DET(L)*VOLINV
         IF (IOPT.EQ.2) TEMP2 = TEMP1/R(L)
!
         DO 200 J=1,NEN
            DO 100 I=1,NESD
               SHGBAR(I,J) = SHGBAR(I,J) + TEMP1*SHG(I,J,L)
 100        CONTINUE
 200     CONTINUE
         IF (IOPT.EQ.2) SHGBAR(3,J)=SHGBAR(3,J)+TEMP2*SHG(3,J,L)
 300  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 780 **************** 
!
      SUBROUTINE SETBB(BBAR,SHG,SHGBAR,R,NROWSH,NROWB,IOPT,IBBAR) 
!
!..... PROGRAM TO SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
!         TWO-DIMENSIONAL CONTINUUM ELEMENTS
!
!         IBBAR = 0,  STANDARD B-MATRIX
!
!         IBBAR = 1,  MEAN-DILATATIONAL B-MATRIX
!

      IMPLICIT NONE
!
!.... DEACTIVATE ABOVE CARD(S) FOR SINGLE-PRECISION OPERATION
!
      integer :: NROWSH,NROWB,IOPT,IBBAR
      real*8 ::  SHG(NROWSH), SHGBAR(NROWSH), BBAR(NROWB,2), R
!
      real*8 :: temp1, temp2, temp3, constb
!
       BBAR(1,1) = SHG(1)
       BBAR(2,1) = 0.0d00
       BBAR(3,1) = SHG(2)
       BBAR(4,1) = 0.0d00
!
       BBAR(1,2) = 0.0d00 
       BBAR(2,2) = SHG(2)
       BBAR(3,2) = SHG(1)
       BBAR(4,2) = 0.0d00
!
!.... AXISYMMETRIC CASE
!
      IF (IOPT.EQ.2) THEN 
        BBAR(4,1) = SHG(3)/R
        BBAR(4,2) = 0.0d00
      ENDIF
!
      IF (IBBAR.EQ.0) RETURN
!
      CONSTB = 1.0d0/3.0d0
!
!.... ADD CONTRIBUTIONS TO FORM B-BAR
!
!.... DEFINE VALUES FOR ROW=4 CASES: PLAIN STRAIN OR AXISYMMETRIC
!
      IF (IOPT.EQ.2) THEN
         TEMP3 = CONSTB*(SHGBAR(3)-SHG(3)/R)
         BBAR(1,1) = BBAR(1,1)+TEMP3
         BBAR(2,1) = BBAR(2,1)+TEMP3
         BBAR(4,1) = BBAR(4,1)+TEMP3
       ELSE
         BBAR(4,1) = 0.0d00
         BBAR(4,1) = 0.0d00
      ENDIF
!
      TEMP1 = CONSTB*(SHGBAR(1)-SHG(1))
      TEMP2 = CONSTB*(SHGBAR(2)-SHG(2))
!
      BBAR(1,1) = BBAR(1,1) + TEMP1
      BBAR(2,1) = BBAR(2,1) + TEMP1
      BBAR(4,1) = BBAR(4,1) + TEMP1
!
      BBAR(1,2) = BBAR(1,2) + TEMP2
      BBAR(2,2) = BBAR(2,2) + TEMP2
      BBAR(4,2) = BBAR(4,2) + TEMP2
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW **************************************************************** 
!
      SUBROUTINE KDBCGEO(ELEFFM,ELRESF,DL,NEE,LM,NEL)
! 
!.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT 
!     BOUNDARY CONDITION 
! 
      USE mGlobaisEscalares
      IMPLICIT REAL*8 (A-H,O-Z)
! 
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION 
! 
      INTEGER :: NEE, NEL
      REAL*8  :: ELEFFM(NEE,*),ELRESF(*),DL(*)
      INTEGER :: LM(*)
!
      INTEGER :: I,J,L
      REAL(8) :: VAL
! 
!    THIS VERSION OF KDBC IS ONLY VALID FOR THERMOELASTIC CONS. 
! 
      DO 200 J=1,NEE
         L   = LM(J) 
         VAL = DL(J)     
         IF(L.GT.0)      GO TO 200
         IF(VAL.EQ.0.d0) GO TO 200 
         DO 100 I=1,NEE
            ELRESF(I)=ELRESF(I)-ELEFFM(I,J)*VAL
100      CONTINUE
! 
200   CONTINUE
! 
      RETURN
!
      END SUBROUTINE
!      
!**** NEW ****************************************************************
!
      SUBROUTINE SIGMAINIT(REGION,SIGMAT,P,NUMEL,NUMELRESERV)
!
      use mPropGeoFisica,    only: POISBTTM, YNGBOTTM, BULKSOLID
!
!.... PROGRAM TO SETUP INITIAL STRESS 
!
      IMPLICIT NONE
!
      INTEGER :: NEL,NUMEL,NUMELRESERV
      REAL(8) :: BULKROCK, BIOTCOEF
!
      INTEGER, DIMENSION(NUMEL)       :: REGION
      REAL(8), DIMENSION(NUMEL)       :: SIGMAT
      REAL(8), DIMENSION(NUMELRESERV) :: P
!
      BULKROCK = BULK(YNGBOTTM,POISBTTM)
!
      BIOTCOEF = 1.0D0 - BULKROCK/BULKSOLID
!
      DO 500 NEL=1,NUMEL
         IF (REGION(NEL).EQ.1) THEN 
               SIGMAT(NEL) = -3.0D0*BIOTCOEF*P(NEL)
           ELSE
               SIGMAT(NEL) = 0.0D0
         ENDIF
 500  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!     
!**** NEW ****************************************************************
!
      SUBROUTINE CALCVDP(ID,VDP,D,DTRL,NDOF,NUMNP,DELTAT)
!
!.... PROGRAM TO COMPUTE DISPLACEMENT VELOCITY
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      INTEGER :: NDOF, NUMNP
      INTEGER :: ID(NDOF,*)
      REAL*8  :: VDP(NDOF,*),D(NDOF,*),DTRL(NDOF,*)
      REAL*8  :: DELTAT
!
      INTEGER :: I,J
!
!      write(*,*) '  '
!      WRITE(*,*) 'AQUI DELTA T --> ',DELTAT
!      write(*,*) '  '

         DO 200 I=1,NDOF
!
         DO 100 J=1,NUMNP
!         K = ID(I,J)
!           IF (K.GT.0) 
           VDP(I,J)= (D(I,J)- DTRL(I,J))/DELTAT
  100    CONTINUE
!
  200    CONTINUE
!
      RETURN
!
      END SUBROUTINE
!      
!**** NEW ****************************************************************
!
      SUBROUTINE CALCVDP_LIN(ID,VDP,D,BRHS,NDOF,NUMNP,DELTAT)
!
!.... PROGRAM TO COMPUTE DISPLACEMENT VELOCITY
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      DIMENSION ID(NDOF,*),VDP(NDOF,*),D(NDOF,*),BRHS(*)
!
!      write(*,*) '  '
!      WRITE(*,*) 'AQUI DELTA T --> ',DELTAT
!      write(*,*) '  '

         DO 200 I=1,NDOF
!
         DO 100 J=1,NUMNP
            K = ID(I,J)
            IF (K.GT.0) then
               VDP(I,J)= (D(I,J)- BRHS(K))/DELTAT
            ENDIF
  100    CONTINUE
!
  200    CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!*** NEW ***** FUNCTION TO COMPUTE ROCK BULK MODULUS  *******************
!
      FUNCTION BULK(YOUNG,POISSON)
!
!... COMPUTE BULK MODULUS OF ROCK BULK MODULUS
!     
      IMPLICIT NONE
!
      REAL(8) :: AMU2, ALAM, YOUNG, POISSON, BULK
!
      AMU2  = YOUNG/(1.0D0 + POISSON)
!
      ALAM  = AMU2*POISSON/(1.0D0-2.0D0*POISSON)
!
      BULK  = ALAM + AMU2/3.0D0
!
      END FUNCTION
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE POS4ITER(x, conecNodaisElem, p, p0, sat)
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
! 
      use mMalha,            only: nsd, numnp, numel, numelReserv, nen, LOCAL, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, PORE, PORE0, REGION, BULKSOLID
      use mPropGeoFisica,    only: POISBTTM, POISRIFT, POISMIDL, POISTOP
      use mPropGeoFisica,    only: MECLAW, GEOFORM
!
      use mPropGeoFisica,    ONLY: BULKWATER, BULKOIL, BULKSOLID
      use mPropGeoFisica,    ONLY: RHOW, RHOO, MASSCONT, MASSCONT0
!
      IMPLICIT NONE
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv), p0(numelReserv), sat(numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
!
      LOGICAL QUAD
!
      INTEGER :: I,J,K,L,NEL
!
!       INTEGER :: NEN,NSD,NESD,NDOF2,NED2,NEE2,NUMEL,NEESQ2,NINTD
!       INTEGER :: NEG,NROWSH,NUMNP,NROWB,IOPTD,IBBAR
        INTEGER :: IOPTD
! 
! 
       REAL(8), DIMENSION(NESD,NEN)        :: XL
       REAL(8), DIMENSION(NED2,NEN)        :: DISL
       REAL(8), DIMENSION(NINTD)           :: WD, DETD, R
       REAL(8), DIMENSION(NROWSH,NEN,NINTD)  :: SHLD,SHGD
       REAL(8), DIMENSION(NROWSH,NEN)    :: SHGBR
       REAL(8), DIMENSION(NROWB,NROWB)   :: CBBAR
!       REAL(8), DIMENSION(NUMEL)         :: PORE,PORE0,DIVU,DIVU0
!       REAL(8), DIMENSION(NUMEL)         :: P,P0,YOUNG
!       REAL(8), DIMENSION(NUMEL)         :: MASCN,MASCN0,PHIEULER
!       REAL(8), DIMENSION(NROWB,*)       :: STRSS
! !
! !.... ..LOCAL VECTORS AND MATRIZES
! !
       REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ
       REAL(8), DIMENSION(NROWB)         :: TENSAO,STRAIN,DEVSTRS
       REAL(8), DIMENSION(4)             :: UNITVEC
! !
       REAL(8) :: POISSON,AREA,C1,TRCSTRS
!       REAL(8) :: ROWDOT,COLDOT,DIFFDIVU,DIFFPRES
! !
!       REAL(8) :: BULK, BIOTMOD, BULKROCK, BIOTCOEF
!       REAL(8) :: BULKFLUID, DEFNM1, DEFMM1, JACOBIAN
!
      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0  
!
      POISSON = POISBTTM


!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)

!       do i=1, numnp
!          write(847,*), dis(1:NDOFD, i)
!       enddo
!        stop 'xx'

!
      DO 500 NEL=1,NUMEL
!
!.... ..LOCALIZE GEOMECHANICAL REGION TO SETUP ELASTIC PARAMETERS
!
         if(novaMalha.eqv..true.) then
!             POISSON = GEODIC('POISSON',GEOFORM(NEL))
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIFT_CAPROCK') POISSON = POISRIFT
            IF (GEOFORM(NEL).EQ.'LEFT__RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIGHT_RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'SALT_CAPROCK') POISSON = POISMIDL
            IF (GEOFORM(NEL).EQ.'POS_SAL_OVER') POISSON = POISTOP
            
         else
            IF (REGION(NEL).EQ.2) POISSON = POISMIDL
            IF (REGION(NEL).EQ.3) POISSON = POISTOP
         endif
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
        CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPTD)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!         write(521,*), "aaaa, NEL=", nel, DISL(1:NED2, 1:NEN)
!
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)

!         write(521,*), "aaaa, NEL=", nel, DETD(1:NINTD)
!
        DIVU(NEL) = 0.0D0
!
!...  ..CLEAR INITIAL STRAIN
!
        STRAIN=0.d0 
!
!.... ..DEFINE ELEMENT AREA
!
         AREA = 0.0D0
!
!.... ..SETUP FOR AXISYMMETRIC OPTION
!
        IF (IOPTD.EQ.2) THEN
          DO 10 L=1,NINTD
             R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
             DETD(L) = DETD(L)*R(L) 
 10       CONTINUE
        ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPTD,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
        DO 300 L=1,NINTD
!
        C1=WD(L)*DETD(L)
!
        AREA = AREA + C1
!
        DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
        CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &             R(L),NROWSH,NROWB,IOPTD,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
        DO 200 K=1,NROWB
           STRAIN(K)=STRAIN(K)+COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200    CONTINUE
!
 300    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION OVER ELEMENT 
!
        DO 350 K=1,NROWB 
           STRAIN(K)=STRAIN(K)/AREA
 350    CONTINUE
!..
         DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
!
!.... ..COMPUTE MEAN VOLUMETRIC STRESS 
!
        DO 420 K=1,NROWB
           STRSS_LIN(K,NEL)=COLDOT(CBBAR(K,1:4),STRAIN(1:4),4)
!             write(521,*), K, NEL, STRSS_LIN(K,NEL)
 420    CONTINUE

        
!
 500  CONTINUE

!          stop "aquiiiiiiii"
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      subroutine POS4TRNS(x, conecNodaisElem, p, p0, sat)
! 
      use mMalha,            only: nsd, numnp, numel, numelReserv, nen, LOCAL, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, PORE, PORE0, REGION, BULKSOLID
      use mPropGeoFisica,    only: POISBTTM, POISRIFT, POISMIDL, POISTOP
      use mPropGeoFisica,    only: MECLAW, GEOFORM
!
      use mPropGeoFisica,    ONLY: BULKWATER, BULKOIL, BULKSOLID
      use mPropGeoFisica,    ONLY: RHOW, RHOO, MASSCONT, MASSCONT0
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
!
      IMPLICIT NONE
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv), p0(numelReserv), sat(numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL QUAD
      real(8) :: xl(nesd,nen), disl(ned2,nen), AREA, xdivu, c1
      REAL*8  :: POISSON
      integer :: nel, l, j, k
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: UNITVEC(NROWB),BBARJ(NROWB,NESD)
      REAL(8) :: CBBAR(NROWB, NROWB)
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
!
      REAL(8) :: PHIEULER(NUMELRESERV)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL*8 :: STRAIN(NROWB), DEVSTRS(NROWB), TENSAO(NROWB) 
      REAL*8 :: ROOT3D2, TRCSTRS, QVM, JACOBIAN
!
      REAL(8) :: POROEXP, BIOTMOD, BULKROCK, BIOTCOEF, DEFNM1
      REAL(8) :: SATDENSW, SATDENSO, RHOEQUIV, BULKFLUID
      REAL(8) :: ALFADIVU, DIFFPRES
!
      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0  
      ROOT3D2=1.224744871391589D0
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
      DO 500 NEL=1,NUMEL
!
         if(novaMalha.eqv..true.) then
!             POISSON = GEODIC('POISSON',GEOFORM(NEL))
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIFT_CAPROCK') POISSON = POISRIFT
            IF (GEOFORM(NEL).EQ.'LEFT__RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIGHT_RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'SALT_CAPROCK') POISSON = POISMIDL
            IF (GEOFORM(NEL).EQ.'POS_SAL_OVER') POISSON = POISTOP
            
         else
            IF (REGION(NEL).EQ.2) POISSON = POISMIDL
            IF (REGION(NEL).EQ.3) POISSON = POISTOP
         endif
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
         CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
         DIVU(NEL) = 0.0D0
!..
!.... ..CLEAR INITIAL STRAIN
!..
         CALL CLEAR(STRAIN,NROWB) 
         CALL CLEAR(TENSAO,NROWB) 
!..
!.... ..DEFINE ELEMENT AREA
!..
         AREA = 0.0D0
!..
!.... ..SETUP FOR AXISYMMETRIC OPTION
!..
         IF (IOPT.EQ.2) THEN
            DO 10 L=1,NINTD
               R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
               DETD(L) = DETD(L)*R(L) 
 10         CONTINUE
         ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
         DO 300 L=1,NINTD
            C1=WD(L)*DETD(L)
            AREA = AREA + C1
            DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &              R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
            DO 200 K=1,NROWB
                   STRAIN(K)=STRAIN(K)+                          &
     &              COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200        CONTINUE
!..
!.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION 
!..
!.... ..TO COMPUTE MEAN VOLUMETRIC CREEP 
!.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
!
!.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
!
            QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!
            DO 220 K=1,NROWB
                   ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+            &
     &                ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
 220        CONTINUE
            DO 250 K=1,NROWB
               TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
 250        CONTINUE
 300    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT 
!
        DO 350 K=1,NROWB       
           STRAIN(K)      = STRAIN(K)/AREA
           AVSTRS(NEL,K)  = TENSAO(K)/AREA
 350    CONTINUE
!        if(nel.eq.77) write(*,*) strain
!        if(nel.eq.77) write(*,*) dis(1,conecNodaisElem(1,77)), dis(2,conecNodaisElem(1,77)) 
!        if(nel.eq.77) write(*,*) dis(1,conecNodaisElem(2,77)), dis(2,conecNodaisElem(1,77)) 
!        if(nel.eq.77) write(*,*) dis(1,conecNodaisElem(3,77)), dis(3,conecNodaisElem(1,77)) 
!        if(nel.eq.77) write(*,*) dis(1,conecNodaisElem(4,77)), dis(4,conecNodaisElem(1,77)) 
!..
!.... ..COMPUTE MEAN VOLUMETRIC DEFORMATION FOR RESERVOIR
!..
        IF (REGION(NEL).EQ.1) THEN
!..
!.... ..COMPUTE BULK MODULUS AND BIOT COEFICIENT(ALPHA)
!..
           BULKROCK  = BULK(YOUNG(NEL),POISBTTM)
           BIOTCOEF  = 1.0D0 - BULKROCK/BULKSOLID   !also called ALPHA
!
           DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
!
!... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27 
!
           JACOBIAN  = 1.0D0 + DIVU(NEL)
!
!... COMPUTE DENSITY AND BULK OF EQUIVALENT FLUID
! 
           SATDENSW = Sat(NEL)*RHOW
           SATDENSO = (1.0D0-Sat(NEL))*RHOO
           RHOEQUIV = SATDENSW + SATDENSO
!
           BULKFLUID = (SATDENSW/BULKWATER + SATDENSO/BULKOIL)/RHOEQUIV
!... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35 
           DEFNM1    = (BIOTCOEF-PORE0(NEL))/BULKSOLID
!... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
           BIOTMOD = DEFNM1 + PORE0(NEL)*BULKFLUID 
!
           ALFADIVU  = BIOTCOEF*(DIVU(NEL)-DIVU0(NEL))
           DIFFPRES  = P(1,NEL)-P0(NEL)
!
!.... ..COMPUTE POROSITY: linearized form Ref: Coussy 2.Edt. Eqs.4.19b
! 
!            PORE(NEL) = PORE0(NEL)+ALFADIVU+DEFNM1*DIFFPRES
           PORE(NEL) = PORE0(NEL)+BIOTCOEF*(DIVU(NEL)-DIVU0(NEL))+DEFNM1*DIFFPRES
!
!.... ..COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
!
           MASSCONT(NEL)=PORE(NEL) !!MASSCONT0(NEL)+RHOEQUIV*(alfaDIVU+BIOTMOD*DIFFPRES)

           if (masscont(nel).lt.0.0d0) then
               write(*,*) 'nel=',nel
               write(*,*) 'RHOEQUIV=',RHOEQUIV
                write(*,*) 'alfa=',BIOTCOEF
               write(*,*) 'divu=',divu(nel)
               write(*,*) 'divu0=',divu0(nel)
               write(*,*) 'biotmod=',biotmod
               write(*,*) 'diffpress=',diffpres
               write(*,*) '1/m(p-p0)=',BIOTMOD*DIFFPRES
               write(*,*) 'masscont=',masscont(nel), masscont0(nel) 
               write(*,*) 'masscont=',masscont(nel-1), masscont0(nel-1) 
!               stop
           endif
!
!.... RELATION OF EULERIAN (PHIEULER) AND LAGRANGIAN (PORE) POROSITIES
!....    REFERENCE COUSSY EQ. 1.17
!
           PHIEULER(NEL) = PORE(NEL)  !/JACOBIAN
!
        ENDIF
!..
 500  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE POS4TRNS_LIN(DIVU,DIVU0, P, P0,PORE,PORE0,YOUNG,&
     &                    MASCN,MASCN0,REGION,PHIEULER,NUMEL)
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
! 
  
      use mPropGeoFisica,    only: POISBTTM, BULKSOLID, BULKWATER    
!
      IMPLICIT NONE
!
      LOGICAL QUAD
!
      INTEGER :: I,J,K,L,NEL, NUMEL
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      INTEGER, DIMENSION(*)         :: REGION
!
      REAL(8), DIMENSION(*)         :: PORE,PORE0,DIVU,DIVU0
      REAL(8), DIMENSION(*)         :: P,P0, YOUNG
      REAL(8), DIMENSION(*)         :: MASCN,MASCN0,PHIEULER
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: POISSON, DIFFDIVU, DIFFPRES
!
      REAL(8) :: BULKROCK, BIOTCOEF
      REAL(8) :: BULKFLUID, DEFNM1, DEFMM1, JACOBIAN
!
      DO 500 NEL=1,NUMEL

         IF (REGION(NEL).EQ.1) THEN

           BULKROCK  = BULK(YOUNG(NEL),POISBTTM)
           BIOTCOEF  = 1.0D0 - BULKROCK/BULKSOLID   !also called ALPHA
!
!... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35 
!
           DEFNM1    = (BIOTCOEF-PORE0(NEL))/BULKSOLID
!
!... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
!
           DEFMM1 = DEFNM1 + PORE0(NEL)/BULKWATER 
!
!... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27 
!
           JACOBIAN  = 1.0D0 + DIVU(NEL)
!
           DIFFDIVU  = DIVU(NEL)-DIVU0(NEL)
           DIFFPRES  = P(NEL)-P0(NEL)
!
! ..COMPUTE LAGRANGIAN POROSITY
!           PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-diffdivu) 
!
!new line: linearized See: Coussy 2.Edt. Eqs. 4.19
           PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES
!
!...COMPUTE EULERAIN POROSITY
           PHIEULER(NEL) = PORE(NEL) !/JACOBIAN
!
!...COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
!
           MASCN(NEL)= MASCN0(NEL)+(BIOTCOEF*DIFFDIVU+DEFMM1*DIFFPRES)
           if (mascn(nel).lt.0.0d0) then
               write(*,*) 'nel=',nel
!               write(*,*) 'RHOEQUIV=',RHOEQUIV
                write(*,*) 'alfa=',BIOTCOEF
               write(*,*) 'divu=',divu(nel)
               write(*,*) 'divu0=',divu0(nel)
               write(*,*) 'biotmod=',defmm1
               write(*,*) 'diffpress=',diffpres
               write(*,*) '1/m(p-p0)=',defmm1*DIFFPRES
               write(*,*) 'masscont=',mascn(nel), mascn0(nel) 
!               write(*,*) 'masscont=',mascn(nel-1), mascn0(nel-1) 
               stop
           endif
!
!             ELSE

!             PORE(NEL)     = 0.0D0
!             MASCN(NEL)    = 0.0D0
!             PHIEULER(NEL) = 0.0D0
         ENDIF
!
 500  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE POS4CREEP(x, conecNodaisElem)
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
! 
      use mMalha, only:nsd, numnp, numel, nen, LOCAL, numelReserv, multab
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh, nnp, NCREEP
      use mAlgMatricial, only: rowdot, coldot
      use mFuncoesDeForma, only: shgq, shlq
      use mPropGeoFisica, only: YOUNG, POISSCAP, PORE, PORE0, POISBTTM, POISMIDL, POISTOP
      use mPropGeoFisica, only: REGION, MECLAW
!
      IMPLICIT NONE
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!

      real*8,  intent(in)    :: x(nsd,numnp)
      integer, intent(in)    :: conecNodaisElem(nen,numel)

      LOGICAL QUAD
!
      INTEGER :: J,K,L,NEL      
!
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      REAL(8), DIMENSION(NESD,NEN)      :: XL
      REAL(8), DIMENSION(NED2,NEN)      :: DLTRL
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ
      REAL(8), DIMENSION(NROWB,NROWB)   :: QIXI, CBBAR
      REAL(8), DIMENSION(NROWB)         :: DEVSTRS,TENSAO
      REAL(8), DIMENSION(4)             :: STRAIN
!
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)

      REAL(8) :: POISSON
      REAL(8) :: GTIMES2, GTIMES3, BULK, TRCSTRS, QTRIAL, ROOT3D2 
      REAL(8) :: GAMMA,FX,DERFX,DELTAG
!
      DATA ROOT3D2/1.224744871391589D0/

! C
! C      write(*,*) 'nnp = ',nnp, 'ncreep =',ncreep,'mod(nnp,ncreep)',
! C     &            mod(nnp,ncreep)
! c      return
! c
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)

      POISSON = POISBTTM
      TENSAO=0.d0
!
      DO 500 NEL=1,NUMEL
!
         IF (REGION(NEL).EQ.2) POISSON = POISMIDL
         IF (REGION(NEL).EQ.3) POISSON = POISTOP
!
         GTIMES2 = YOUNG(NEL)/(1.0D0+POISSON)
         GTIMES3 = 1.5D0*GTIMES2
!
         BULK = YOUNG(NEL)/(3.0D0*(1.0D0-2.0D0*POISSON))    
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
         CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!..... ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DTRL,DLTRL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.

         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!..
!..... ..SETUP FOR AXISYMMETRIC OPTION
!..
         IF (IOPT.EQ.2) THEN
            DO 100 L=1,NINTD
               R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
               DETD(L) = DETD(L)*R(L) 
 100        CONTINUE
         ENDIF
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!....     FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!..
!.... ...LOOP OVER INTEGRATION POINTS
!..

         DO 400 L=1,NINTD
!..
!..... ..CLEAR INITIAL STRAIN
!..
            STRAIN=0.d0
!
            DO 200 J=1,NEN  
!..
!.... ..... UPLOAD B-BAR MATRIX AT NODE J
!..
               CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &                    R(L),NROWSH,NROWB,IOPT,IBBAR)
!.. 
!.... ..... COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!..
               DO 200 K=1,NROWB
                  STRAIN(K)=STRAIN(K)+ &
     &                    COLDOT(BBARJ(K,1:2),DLTRL(1:2,J),2)
 200        CONTINUE
!..
!.... ..... COMPUTE STRESS 
!..
!            CALL MULTAB(CBBAR,STRAIN,TENSAO,4,4,4,4,4,1,1)
             TENSAO=matmul(CBBAR,STRAIN)
!..
!.... ..... COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,TENSAO,NROWB)
!
!.... ..... COMPUTE EFFECTIVE TRIAL STRESS
!
            QTRIAL=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!..
!... SOLVE NON LINEAR EQUATION FOR GAMMA FRAMEWORK: NEWTON METHOD
!..
            GAMMA=0.0D0
!
            IF (MOD(NNP,NCREEP).EQ.0) THEN
               IF (MECLAW(NEL).EQ.2) THEN 
 220              CONTINUE
                  FX =  FNOTLIN(QTRIAL,GTIMES3,1,GAMMA)
                  DERFX = FNOTLIN(QTRIAL,GTIMES3,2,GAMMA) 
                  DELTAG = -FX/DERFX
                  GAMMA = GAMMA+DELTAG
                  IF (DABS(DELTAG).GT.1.0D-6) GOTO 220
               ENDIF
            ENDIF
!..
!.... ...  UPDATE DEVIATOR STRESS
!..
            DO 230 K=1,NROWB
               DEVSTRS(K)=(1.0D0-GAMMA*GTIMES3/QTRIAL)*DEVSTRS(K)
 230        CONTINUE
!
            STRSS(NEL,L,1) = DEVSTRS(1)+TRCSTRS
            STRSS(NEL,L,2) = DEVSTRS(2)+TRCSTRS
            STRSS(NEL,L,3) = DEVSTRS(3)
            STRSS(NEL,L,4) = DEVSTRS(4)+TRCSTRS

!..
!.... .... UPDATE LOCAL TANGENT MATRIX
!..
! 
            CALL SETTGMX(QIXI,DEVSTRS,GTIMES2,GTIMES3,BULK,QTRIAL, &
     &                   GAMMA,MECLAW(NEL))
!..
!.... .... TRANSFER 4X4-ORDER MATRIX TO GLOBAL TANGENT ARRAY
!..
            CALL QX2TANG(QIXI,HMTTG(NEL,L,1:16))
!
 400     CONTINUE
!
 500  CONTINUE

      RETURN 
!
2222  FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X  &
     &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                  &
     &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                  &
     &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                  &
     &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
 4000 FORMAT(2X,40(1PE15.8,2X)) 
 5000 FORMAT(I4,2X,I1,2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE

!
!**** new **********************************************************************
!
      subroutine UPDATSAT(s, p, phi)
!
      use mMalha,            only: nsd, numelReserv
      use mGlobaisEscalares, only: ndofP
      use mPropGeoFisica,    only: phi0
!
      implicit none
!
      real*8 :: s(numelReserv), p(ndofP,numelReserv), phi(numelReserv)
      integer :: nel

!
!.... SATURATION CORRECTION
!
      do nel=1,numelReserv
         s(nel)=phi0(nel)/phi(nel)*s(nel)!*dexp(-beta*(p(nel)-p0(nel)))
      END DO ! nel

      end subroutine
!
!**** NEW **** FOR STOCHASTIC AND NON-LINEAR FORMULATION ***************** 
!
      SUBROUTINE GEOSETUP(YOUNG,REGION, &
     &                    NUMEL,NROWB,NINTD,NROWB2,IOPT)
!

     use mPropGeoFisica, only:  POISBTTM, POISMIDL, POISTOP

      IMPLICIT NONE
!
!.... PROGRAM TO SETUP INITIAL INELASTIC STOCHASTIC PARAMETERS 
! 
      INTEGER :: NEL,NUMEL,NROWB,NINTD,NROWB2,IOPT
!
!.... NEXT LINES: GLOBAL VECTORS AND MATRIZES
!
      INTEGER, DIMENSION(NUMEL)              :: REGION
!
      REAL(8), DIMENSION(NUMEL)              :: YOUNG
      REAL(8) :: CBBAR(NROWB, NROWB)
! 
      REAL(8) :: POISSON
      integer :: L
!
      DO 500 NEL=1,NUMEL 
!
      IF (REGION(NEL).EQ.1) POISSON = POISBTTM
      IF (REGION(NEL).EQ.2) POISSON = POISMIDL
      IF (REGION(NEL).EQ.3) POISSON = POISTOP
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
      CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!..
!.... SETUP INITIAL TANGENT MATRIX AT ELEMENT GAUSS POINT 
!..
         DO 100 L=1,NINTD
            CALL QX2TANG(CBBAR,HMTTG(NEL,L,1:16))
 100  CONTINUE
!
 500  CONTINUE
!
!.... SETUP INITIAL VISCO-ELASTIC DEFORMATION AT GAUSS POINTS
!
      CALL CLEAR(ECREEP,NUMEL*NINTD*NROWB)
!
      RETURN
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR STOCHASTIC AND NON-LINEAR FORMULATION ***************** 
!
      SUBROUTINE GEOSETUP1(YOUNG,GEOFORM, &
     &                     NUMEL,NROWB,NINTD,NROWB2,IOPT)
!
      use mPropGeoFisica,    only: POISBTTM, POISRIFT, POISMIDL, POISTOP

      IMPLICIT NONE
!
!.... PROGRAM TO SETUP INITIAL INELASTIC STOCHASTIC PARAMETERS 
! 
      INTEGER :: NEL,NUMEL,NROWB,NINTD,NROWB2,IOPT
!
!.... NEXT LINES: GLOBAL VECTORS AND MATRIZES
!
      CHARACTER*12, DIMENSION(NUMEL)         :: GEOFORM
      REAL(8), DIMENSION(NUMEL)              :: YOUNG
      REAL(8) :: CBBAR(NROWB, NROWB)
! 
      REAL(8) :: POISSON
      integer :: L
!
      DO 500 NEL=1,NUMEL 
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIFT_CAPROCK') POISSON = POISRIFT
            IF (GEOFORM(NEL).EQ.'LEFT__RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'RIGHT_RESERV') POISSON = POISBTTM
            IF (GEOFORM(NEL).EQ.'SALT_CAPROCK') POISSON = POISMIDL
            IF (GEOFORM(NEL).EQ.'POS_SAL_OVER') POISSON = POISTOP
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
!
         CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!.... SETUP INITIAL TANGENT MATRIX AT ELEMENT GAUSS POINT 
!
         DO 100 L=1,NINTD
            CALL QX2TANG(CBBAR,HMTTG(NEL,L,1:16))
 100     CONTINUE
 500  CONTINUE
!
!.... SETUP INITIAL VISCO-ELASTIC DEFORMATION AT GAUSS POINTS
!
      CALL CLEAR(ECREEP,NUMEL*NINTD*NROWB)
!
      RETURN
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW ****************************************************************
!
      FUNCTION POWERLAW(X,NDERHO)
!
      USE mPropGeoFisica, only: DTCREEP, CREEPZERO, POWERN, SIGMAREF
!
!.....PROGRAM TO COMPUTE POWER LAW LIKE FUCNTION 
!
      IMPLICIT NONE
!
      REAL*8 :: X
      INTEGER NDERHO
!
      REAL*8 :: POWERLAW, XIMAGE
!
      XIMAGE=(DTCREEP*CREEPZERO)*(X/SIGMAREF)**POWERN
!
      GOTO(100,200) NDERHO
!
!.... FUNCTION ONLY
!
 100  CONTINUE
        POWERLAW = XIMAGE
      RETURN
!..
!.... FIRST DERIVATIVE OF POWER LAW 
!..
 200  CONTINUE
        POWERLAW = POWERN*XIMAGE/X
      RETURN
!
      END FUNCTION
!
!**** NEW ****************************************************************
!
      FUNCTION FNOTLIN(QTRIAL,THREEG,NDERHO,XINPUT)
!
!
!.....PROGRAM TO COMPUTE NON LINEAR FUNCTION FROM POWER LAW
!
      IMPLICIT NONE
!
      INTEGER :: NDERHO
      REAL(8) :: QTRIAL, THREEG, X, XINPUT, FNOTLIN
!
      X = QTRIAL-THREEG*XINPUT
!
      GOTO(100,200) NDERHO
!
!... FUNCTION
!
 100  CONTINUE
        FNOTLIN = XINPUT-POWERLAW(X,NDERHO)
      RETURN
!
!... FIRST DERIVATIVE FUNCTION
!
 200  CONTINUE
        FNOTLIN = 1.0D0+POWERLAW(X,NDERHO)*THREEG
      RETURN
!
      END FUNCTION
!
!**** NEW ****************************************************************
!
      FUNCTION DEVNORM2(X,NROWB)
!
!.....PROGRAM TO COMPUTE SQUARE NORM OF DEVIATOR TENSOR
!
      IMPLICIT NONE
!
      INTEGER :: NROWB
      REAL*8  :: X(NROWB)
!
      REAL*8  :: DEVNORM2
!
      DEVNORM2 = X(1)*X(1)+X(2)*X(2)+X(4)*X(4)+2.0D0*X(3)*X(3)
!
      RETURN
!
      END FUNCTION
!      
!**** NEW ****************************************************************
!
      SUBROUTINE COMPTRACE(REGION,P,AVSTRS,SIGMAT,NROWB,NUMEL,numelReserv)
!
!.... PROGRAM TO COMPUTE TRACE OF TOTAL STRESS
!.... TOTAL STRESS = SOLID_STRSS_TRACE+FLUID_PRESSURE
!
      use mPropGeoFisica,    only: POISBTTM, YNGBOTTM, BULKSOLID
!
      IMPLICIT NONE
!
!
      INTEGER :: NEL,NUMEL,NROWB,numelReserv
!
      INTEGER, DIMENSION(NUMEL)       :: REGION
      REAL(8), DIMENSION(numelReserv) :: P
      REAL(8), DIMENSION(NUMEL)       :: SIGMAT
      REAL(8), DIMENSION(NUMEL,NROWB) :: AVSTRS
!      
      REAL(8) :: BULKROCK, BIOTCOEF
!
      BULKROCK = BULK(YNGBOTTM,POISBTTM)
!
      BIOTCOEF = 1.0D0 - BULKROCK/BULKSOLID
!      

!      SIGMAT(NEL)=STRSS(NEL,1,1)+STRSS(NEL,1,2)+STRSS(NEL,1,4)
!     &            -3.0D0*P(NEL)
!
      DO 500 NEL=1,NUMEL
!
      IF (REGION(NEL).EQ.1) THEN
             SIGMAT(NEL)=AVSTRS(NEL,1)+AVSTRS(NEL,2)+AVSTRS(NEL,4) &
     &                  -3.0D0*BIOTCOEF*P(NEL)
         ELSE
             SIGMAT(NEL) = 0.0D0
      ENDIF
!
500   CONTINUE

      RETURN
!
      END SUBROUTINE
!      
!**** NEW ****************************************************************
!
      SUBROUTINE MASSINIT(MASSCONT,PHI,RHOagua,RHOoleo,SAT,NUMELRESERV)
!
!.... PROGRAM TO SETUP INITIAL MASS CONTENT 
!
      IMPLICIT NONE
!
      INTEGER :: NEL,NUMELRESERV
      REAL(8) :: RHOagua, RHOoleo
!
      REAL(8), DIMENSION(NUMELRESERV)   :: MASSCONT, PHI, SAT
!      
      DO 500 NEL=1,NUMELRESERV
         MassCont(NEL) = Phi(NEL)*(RHOagua*Sat(nel)+RHOoleo*(1.0D0-Sat(nel)))
500   CONTINUE
!
      RETURN
!
      END SUBROUTINE
!      
!**** NEW ****************************************************************
!
      SUBROUTINE COMPTRACE_LIN(P,SIGMAT,NUMEL)
!
!.... PROGRAM TO COMPUTE TRACE OF TOTAL STRESS
!
      USE mMalha, only: numelReserv
!
      IMPLICIT NONE
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      REAL*8 :: P(numelReserv),SIGMAT(*)
      integer :: nel, numel
!      
!... COMPUTE COMPRESSIBILITY 
!
      DO 200 NEL=1,numelReserv
!
!.... COMPUTE TOTAL STRESS = SOLID_STRSS_TRACE+FLUID_PRESSURE
!
      SIGMAT(NEL)=STRSS_LIN(1,NEL)+STRSS_LIN(2,NEL)+STRSS_LIN(4,NEL)-3.0D0*P(NEL)
!
  200 CONTINUE
!
      DO 300 NEL=numelReserv+1,NUMEL
         SIGMAT(NEL)=0.0D0
 300  CONTINUE

      RETURN
!
      END SUBROUTINE
!
!**** NEW ****************************************************************
!
      SUBROUTINE VELSIGMAT(DTSIGM,SIGMAT,SIGMA0,NUMEL,DELTAT)
!
!.... PROGRAM TO COMPUTE STRESS TRACE VELOCITY
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION
!
      REAL*8  :: DTSIGM(*),SIGMAT(*),SIGMA0(*), DELTAT
      INTEGER :: NUMEL
!
      INTEGER :: I
!
!      write(*,*) '  '
!      WRITE(*,*) 'AQUI DELTA T --> ',DELTAT
!      write(*,*) '  '

         DO 200 I=1,NUMEL
            DTSIGM(I) = (SIGMAT(I) - SIGMA0(I))/DELTAT
  200    CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW **** SUBROUTINE MODIFIED FOR GEOMECHANIC  ************* 
!
      SUBROUTINE PRTGLOBAL(NSD,NUMEL,T0,X,U,IUNIT)
!
!.... PROGRAM TO PRINT GLOBAL MASS BALANCE FOR MATLAB
!     
      IMPLICIT NONE 
!
      INTEGER                   :: NUMEL,NSD
      REAL(8), DIMENSION(NSD,*) :: X
      REAL(8), DIMENSION(*)     :: U
      REAL(8)                   :: T0
!     
      INTEGER :: NEL
!     
      INTEGER :: IUNIT
!     
      WRITE(IUNIT,"('#TIMESTEP PRINT OUT = ',F15.8)") T0
!     
      DO 100 NEL=1,NUMEL
         WRITE(IUNIT,"(5(F25.15,2X))") X(1:NSD,NEL),U(NEL)
 100  CONTINUE
!
      END SUBROUTINE
!
!**** NEW **** SUBROUTINE MODIFIED FOR GEOMECHANIC  ************* 
!
      SUBROUTINE PRTSIGMAT(NSD,NUMEL,T0,X,U,IUNIT)
!
!.... PROGRAM TO PRINT GLOBAL MASS BALANCE FOR MATLAB
!     
      IMPLICIT NONE 
!
      INTEGER                   :: NUMEL,NSD
      REAL(8), DIMENSION(NSD,*) :: X
      REAL(8), DIMENSION(*)     :: U
      REAL(8)                   :: T0
!     
      INTEGER :: NEL
!     
      INTEGER :: IUNIT
!     
      WRITE(IUNIT,"('#TIMESTEP PRINT OUT = ',F15.8)") T0
!     
      DO 100 NEL=1,NUMEL
         WRITE(IUNIT,"(5(F25.15,2X))") X(1:NSD,NEL),U(NEL)
 100  CONTINUE
!   
      END SUBROUTINE
!
!**** NEW **** FOR STOCHASTIC AND NON-LINEAR FORMULATION ***************** 
!
      SUBROUTINE GEOREGION(REGION, MECLAW, NUMEL, numelReserv, NELDOMO)
!
!.... PROGRAM TO SETUP GEOMECHANICAL REGIONS: RESERVOIR-PRE-SAL, DOMO, POS-SAL
! 
      IMPLICIT NONE 
!  
      INTEGER :: NEL,NUMEL,numelReserv,NELDOMO
!
      INTEGER, DIMENSION(NUMEL) :: REGION, MECLAW
!
      DO 500 NEL=1,NUMEL
      IF (NEL.LE.numelReserv) THEN
             REGION(NEL) = 1
             MECLAW(NEL) = 1
      ENDIF
      IF ((NEL.GT.numelReserv).AND.(NEL.LE.NELDOMO)) THEN 
             REGION(NEL) = 2
             MECLAW(NEL) = 2
      ENDIF
      IF (NEL.GT.NELDOMO) THEN 
             REGION(NEL) = 3
             MECLAW(NEL) = 1
      ENDIF
!
 500  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW **** MODIFIED 4 HIERARCHICAL MESH  ***************** 
!
      SUBROUTINE GEOREGION1(REGION,MECLAW,GEOFORM,NUMEL)
!
!.... PROGRAM TO SETUP GEOMECHANICAL REGIONS: RESERVOIR-PRE-SAL, DOMO, POS-SAL
! 
      IMPLICIT NONE 
!  
      INTEGER :: NEL, NUMEL, INREGION, INMECLAW, INGEOFOR
!
      CHARACTER*30 :: NINREGION, NINMECLAW, NINGEOFOR
!
      INTEGER, DIMENSION(NUMEL) :: REGION,MECLAW
      INTEGER, DIMENSION(NUMEL) :: IREG, IMEC
      CHARACTER*12, DIMENSION(NUMEL) :: GEOLOC, GEOFORM
!
      INREGION = 752
      INMECLAW = 753
      INGEOFOR = 754
      NINREGION = 'data_elmnt_region.in'
      NINMECLAW = 'data_elmnt_meclaw.in'
      NINGEOFOR = 'elmnt_geoformation.in'
      OPEN(UNIT=INREGION,FILE=NINREGION)
      OPEN(UNIT=INMECLAW,FILE=NINMECLAW)
      OPEN(UNIT=INGEOFOR,FILE=NINGEOFOR)
!
      DO 300 NEL=1,NUMEL
         READ(INREGION,1000) IREG(NEL)
!         WRITE(*,1000) IREG(NEL)
         READ(INMECLAW,1000) IMEC(NEL)
!         WRITE(*,1000) IMEC(NEL)
         READ(INGEOFOR,1010) GEOLOC(NEL)
!         WRITE(555,1010) GEOLOC(NEL)         
 300  CONTINUE
!
      CLOSE(INREGION)
      CLOSE(INMECLAW)
      CLOSE(INGEOFOR)
!
      REGION  = IREG
      MECLAW  = IMEC
      GEOFORM = GEOLOC 
!
      RETURN
!
1000  FORMAT(I10)
1010  FORMAT(A12)
!
      END SUBROUTINE
!
!**** NEW ***** FOR SISMIC REPRESENTATION ******************************
!
!       FUNCTION GEODIC(TASK,GEOFORM)
! !
!       use mPropGeoFisica, only: POISBTTM, POISMIDL, POISTOP, POISRIFT
!       use mPropGeoFisica, only: YNGBOTTM, YNGMIDLE, YNGTOP, YNGRIFT
! !
!       IMPLICIT NONE
! !
!       INTEGER      :: I
!       CHARACTER*7  :: TASK
!       CHARACTER*12 :: GEOFORM
!       CHARACTER*12, DIMENSION(6) :: REFGEO
!       REAL*8, DIMENSION(6)       :: POISSON, MODYOUNG
!       REAL*8                     :: GEOLOC, GEODIC
! !
!       DATA REFGEO(1)     ,     REFGEO(2),     REFGEO(3)/ &
!      &     'RESERVATORIO','RIFT_CAPROCK','RIGHT_RESERV'/ &
!      &     REFGEO(4)     ,     REFGEO(5),     REFGEO(6)/&
!      &     'SALT_CAPROCK','LEFT__RESERV','POS_SAL_OVER'/
! !
!         POISSON(1) = POISBTTM
!         POISSON(2) = POISRIFT
!         POISSON(3) = POISBTTM
!         POISSON(4) = POISMIDL
!         POISSON(5) = POISBTTM
!         POISSON(6) = POISTOP
! !
!         MODYOUNG(1) = YNGBOTTM
!         MODYOUNG(2) = YNGRIFT
!         MODYOUNG(3) = YNGBOTTM
!         MODYOUNG(4) = YNGMIDLE
!         MODYOUNG(5) = YNGBOTTM
!         MODYOUNG(6) = YNGTOP
! !
!       IF (TASK.EQ.'POISSON') THEN 
!          DO 100 I=1,6
!             IF (REFGEO(I).EQ.GEOFORM) GEOLOC = POISSON(I)    
! 100      CONTINUE
!       ENDIF
! !
!       IF (TASK.EQ.'YOUNGMD') THEN
!          DO 200 I=1,6
!             IF (REFGEO(I).EQ.GEOFORM) GEOLOC = MODYOUNG(I)    
! 200      CONTINUE
!       ENDIF
! !
!       GEODIC = GEOLOC
! !
!       END FUNCTION
!
!**** NEW ***** FOR CREEP MODELING *************************************
!
      SUBROUTINE DEVTENSOR(DEVIATOR,TRACED3,TENSORIN,NROWB)
!..
!.... PROGRAM TO COMPUTE DEVIATOR FOR PLANE DEFORMATIONS STATE  
!..
      IMPLICIT NONE
!
      INTEGER :: NROWB
      REAL*8 :: DEVIATOR(NROWB), TENSORIN(NROWB), TRACED3
!
      TRACED3 = (TENSORIN(1)+TENSORIN(2)+TENSORIN(4))/3.0D0
!
      DEVIATOR(1)= TENSORIN(1)-TRACED3
      DEVIATOR(2)= TENSORIN(2)-TRACED3
      DEVIATOR(3)= TENSORIN(3)
      DEVIATOR(4)= TENSORIN(4)-TRACED3
!
      RETURN
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE SETTGMX(QIXI,DEVSTRS,GTIMES2,GTIMES3,BULK,QTRIAL,GAMMA,MECLAW)
! 
!     PROGRAM TO SETUP TANGENT MATRIX 
      USE mGlobaisEscalares, only: nnp, NCREEP
! 
      IMPLICIT NONE
! 
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
! 
      REAL*8  :: GTIMES2, GTIMES3, BULK, QTRIAL, GAMMA
      INTEGER :: MECLAW
!
      REAL*8 :: UNITVECT(4),UNITTENS(4,4),DEVPROJ(4,4),QIXI(4,4)
      REAL*8 :: DEVSTRS(4), VECTNORM, X, BLOCOA, BLOCO2G, BLOCO6G2, BLOCOG2
      INTEGER :: I, J
     
!
      UNITVECT(1)=1.0D0
      UNITVECT(2)=1.0D0
      UNITVECT(3)=0.0D0
      UNITVECT(4)=1.0D0
!
      UNITTENS(1,1)=1.0D0
      UNITTENS(1,2)=0.0D0
      UNITTENS(1,3)=0.0D0
      UNITTENS(1,4)=0.0D0
      UNITTENS(2,1)=0.0D0
      UNITTENS(2,2)=1.0D0
      UNITTENS(2,3)=0.0D0
      UNITTENS(2,4)=0.0D0
      UNITTENS(3,1)=0.0D0
      UNITTENS(3,2)=0.0D0
      UNITTENS(3,3)=0.5D0
      UNITTENS(3,4)=0.0D0 
      UNITTENS(4,1)=0.0D0 
      UNITTENS(4,2)=0.0D0 
      UNITTENS(4,3)=0.0D0 
      UNITTENS(4,4)=1.0D0 
!..
!... MOUNT DEVIATOR PROJETOR
!..
      DO 20 I=1,4
         DO 10 J=1,4
            DEVPROJ(I,J)=UNITTENS(I,J)-UNITVECT(I)*UNITVECT(J)/3.0D0
 10      CONTINUE
 20   CONTINUE
!
!      WRITE(102,2222) 1,1, ((DEVPROJ(I,J),I=1,4),J=1,4)
!..
!.... COMPUTE DEVIATOR NORM
!
       VECTNORM = DEVNORM2(DEVSTRS,4)
!..
!.... COMPUTE FACTORS THAT MULTIPLY FOURTH-ORDER MATRICES
!..
       X = QTRIAL-GTIMES3*GAMMA
       BLOCOA   = POWERLAW(X,2)/(1.0D0+GTIMES3*POWERLAW(X,2))
       BLOCO2G  = GTIMES2*(1.0D0-GAMMA*GTIMES3/QTRIAL)
       BLOCO6G2 = GTIMES2*GTIMES3*(GAMMA/QTRIAL-BLOCOA)/VECTNORM
!
       IF ((MECLAW.EQ.1).OR.(MOD(NNP,NCREEP).GT.0)) THEN
           BLOCO6G2 = 0.0D0
           BLOCOG2  = GTIMES2
       ENDIF
       DO 400 I=1,4
        DO 300 J=1,4
           QIXI(I,J)=BLOCO2G*DEVPROJ(I,J)+BULK*UNITVECT(I)*UNITVECT(J) &
     &               +BLOCO6G2*DEVSTRS(I)*DEVSTRS(J)
 300    CONTINUE
 400   CONTINUE
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
       END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE QX2TANG(QMATR4X4,TANGENT)
! 
!     PROGRAM TO TRANSFER FROM QIXI 4X4 MATRIX TO TANGENT ARRAY
! 
      IMPLICIT NONE
! 
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
! 
      REAL*8 :: TANGENT(16),QMATR4X4(4,4)
!
            TANGENT(1)  = QMATR4X4(1,1)
            TANGENT(2)  = QMATR4X4(1,2)
            TANGENT(3)  = QMATR4X4(1,3)
            TANGENT(4)  = QMATR4X4(1,4)
            TANGENT(5)  = QMATR4X4(2,1)
            TANGENT(6)  = QMATR4X4(2,2)
            TANGENT(7)  = QMATR4X4(2,3)
            TANGENT(8)  = QMATR4X4(2,4)
            TANGENT(9)  = QMATR4X4(3,1)
            TANGENT(10) = QMATR4X4(3,2)
            TANGENT(11) = QMATR4X4(3,3)
            TANGENT(12) = QMATR4X4(3,4)
            TANGENT(13) = QMATR4X4(4,1)
            TANGENT(14) = QMATR4X4(4,2)
            TANGENT(15) = QMATR4X4(4,3)
            TANGENT(16) = QMATR4X4(4,4)
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW *** FOR PLASTICITY ********************************************* 
!
      SUBROUTINE UPDATEINCR(IDDIS,BRHSD,NDOFD,NUMNP)
! 
!.... UPDATE TRIAL DISPLACEMENT ARRAY 
! 
      IMPLICIT NONE
! 
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION 
! 
      INTEGER :: NDOFD,NUMNP
      INTEGER :: IDDIS(NDOFD,*)
      REAL*8  :: BRHSD(*)
!
      INTEGER :: I, J, K
! 
!      write(32,*) 'ENTRANDO update iter =',nwtniter
!      call xxwrite(dtrl,ndof2,numnp)
!      call xxwrite(brhsd,1,10000)

      DO 200 I=1,NDOFD
! 
        DO 100 J=1,NUMNP
           K = IDDIS(I,J)
           IF (K.GT.0) then
           DTRL(I,J) = DTRL(I,J) + BRHSD(K)
           ENDIF
  100   CONTINUE
! 
  200 CONTINUE
!
!      write(32,*) 'SAINDO update iter =',nwtniter
!      call xxwrite(dtrl,ndof2,numnp)
!      call xxwrite(brhsd,1,10000)
      RETURN
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE TANG2QX(TANGENT,QMATR4X4)
! 
!     PROGRAM TO TRANSFER FROM TANGENT ARRAY TO QIXI 4X4 MATRIX 
!
      IMPLICIT NONE
! 
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
! 
      real*8 :: TANGENT(16),QMATR4X4(4,4)
!
         QMATR4X4(1,1) = TANGENT(1)
         QMATR4X4(1,2) = TANGENT(2)
         QMATR4X4(1,3) = TANGENT(3)
         QMATR4X4(1,4) = TANGENT(4)
         QMATR4X4(2,1) = TANGENT(5)
         QMATR4X4(2,2) = TANGENT(6)
         QMATR4X4(2,3) = TANGENT(7)
         QMATR4X4(2,4) = TANGENT(8)
         QMATR4X4(3,1) = TANGENT(9)
         QMATR4X4(3,2) = TANGENT(10)
         QMATR4X4(3,3) = TANGENT(11)
         QMATR4X4(3,4) = TANGENT(12)
         QMATR4X4(4,1) = TANGENT(13)
         QMATR4X4(4,2) = TANGENT(14)
         QMATR4X4(4,3) = TANGENT(15)
         QMATR4X4(4,4) = TANGENT(16)
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE PRINT_DX(DIS ,PORE  ,P     ,S     ,PERM  ,YOUNG ,BALANC, &
     &                    VC  ,AVSTRS,REGION,NDOF2 ,NUMEL ,NROWB ,NUMNP)
!
      USE mGlobaisEscalares, only: NUMDX, NITGEO, NNP
      use mLeituraEscrita,   only: PRINT_DXINFO
!
!.... PROGRAM TO PRINT DATA FROM GEOMECHANIC MODEL
! 
      IMPLICIT NONE
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      CHARACTER*30 NIDISP,NIPRSR,NIPORE,NICREP,NISATR,NIVELT
      CHARACTER*30 NISIGX,NISIGY,NISGTA,NISIGZ
      CHARACTER*30 NIYUNG,NIPERM,NIBLNC
!
      CHARACTER*3 ASTEP
      CHARACTER*2 AITER
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      INTEGER :: I,J,K,NEL,NDOF2,NUMEL,NROWB,NUMNP
      INTEGER :: IDISP, IPRSR, IPORE, ICREP, ISATR, IVELT
      INTEGER :: IPERM, ISIGX, ISIGY, ISGTA, ISIGZ, IYUNG, IBLNC
!
!      DIMENSION DIS(NDOF2,*),PORE(*),P(*),S(*),PERM(*),YOUNG(*),
!     &          VC(NDOF2,*),AVSTRS(NUMEL,*)
!
      INTEGER, DIMENSION(NUMEL)       :: REGION
      REAL(8), DIMENSION(NDOF2,NUMNP) :: DIS
      REAL(8), DIMENSION(NDOF2,*)     :: VC
      REAL(8), DIMENSION(NUMEL)       :: PORE,P,S,PERM,YOUNG,BALANC
      REAL(8), DIMENSION(NUMEL,NROWB) :: AVSTRS
!
      REAL(8), DIMENSION(NROWB)       :: TENSAO, DEVSTRS
!
      REAL(8) :: ROOT3D2,TRCSTRS,QTRIAL,STRIAL,SPORE,SBLNC
!
      DATA ROOT3D2/1.224744871391589D0/
!
!      COMMON /IOUNIT2/ IFEDX,IFNOD,IFMSH,IRBAL
!
!... OUTPUT DATA FILES 
!
      IF (NUMDX.EQ.0) RETURN
!
      IDISP = 631
      IPRSR = 632
      IPORE = 633
      ICREP = 634
      ISATR = 635
      IVELT = 636
      IPERM = 637
      ISIGX = 638
      ISIGY = 639
      ISGTA = 640
      ISIGZ = 641
      IYUNG = 642
      IBLNC = 643
!
!.... SETUP FILES PATH
!
      WRITE(AITER,'(I2.2)') NITGEO 
!
!.... SETUP FILES COUNTER
!
!      WRITE(ASTEP,'(I3.3)') IDINT(TPRT_PHI/DTPRT_PHI)
      WRITE(ASTEP,'(I3.3)') NNP/NUMDX
!
!      WRITE(ASTEP,'(I3.3)') IDINT(T0/TPRT_PHI)
!      WRITE(ASTEP,'(I4.4)') 100*NNP+NLOOPS
!
!.... OUT-PUT FILES FOR OPEN-DX VISUALIZATION
!.... DATA VALUES AT NODAL FEM POINTS
!
      NIDISP = 'dxcreep'//AITER//'/disp'//ASTEP//'.stoc'
      NIPRSR = 'dxcreep'//AITER//'/prsr'//ASTEP//'.stoc'
      NIPORE = 'dxcreep'//AITER//'/pore'//ASTEP//'.stoc'
      NICREP = 'dxcreep'//AITER//'/crep'//ASTEP//'.stoc'
      NISATR = 'dxcreep'//AITER//'/satr'//ASTEP//'.stoc'
      NIVELT = 'dxcreep'//AITER//'/velt'//ASTEP//'.stoc'
      NIPERM = 'dxcreep'//AITER//'/perm'//ASTEP//'.stoc'
      NISIGX = 'dxcreep'//AITER//'/sigx'//ASTEP//'.stoc'
      NISIGY = 'dxcreep'//AITER//'/sigy'//ASTEP//'.stoc'
      NISGTA = 'dxcreep'//AITER//'/sigt'//ASTEP//'.stoc'
      NISIGZ = 'dxcreep'//AITER//'/sigz'//ASTEP//'.stoc'
      NIYUNG = 'dxcreep'//AITER//'/yung'//ASTEP//'.stoc'
      NIBLNC = 'dxcreep'//AITER//'/blnc'//ASTEP//'.stoc'
!
!.....OPEN NODAL DATA FILES
!
      OPEN(UNIT=IDISP, FILE= NIDISP)
      OPEN(UNIT=IPRSR, FILE= NIPRSR)
      OPEN(UNIT=IPORE, FILE= NIPORE)
      OPEN(UNIT=ICREP, FILE= NICREP)
      OPEN(UNIT=ISATR, FILE= NISATR)
      OPEN(UNIT=IVELT, FILE= NIVELT)
      OPEN(UNIT=IPERM, FILE= NIPERM)
      OPEN(UNIT=ISIGX, FILE= NISIGX)
      OPEN(UNIT=ISIGY, FILE= NISIGY)
      OPEN(UNIT=ISGTA, FILE= NISGTA)
      OPEN(UNIT=ISIGZ, FILE= NISIGZ)
      OPEN(UNIT=IYUNG, FILE= NIYUNG)
      OPEN(UNIT=IBLNC, FILE= NIBLNC)
!
!.... PRINT DISPLACEMENTS
!
       DO 30 I=1,NUMNP
          WRITE(IDISP,4000) (DIS(J,I),J=1,NDOF2) 
 30    CONTINUE
!
      DO 500 NEL=1,NUMEL 
!
!.... PRINT PRESSURE
!
      IF (REGION(NEL).EQ.1) SPORE = P(NEL)
      IF (REGION(NEL).GT.1) SPORE = -10.0D0
      WRITE(IPRSR,4000) SPORE
!
!.... PRINT GEOMECHANIC POROSITY: "PORE" 
! 
      IF (REGION(NEL).EQ.1) SPORE = PORE(NEL)
      IF (REGION(NEL).GT.1) SPORE = -10.0D0
      WRITE(IPORE,4000) SPORE
!
!.... PRINT STRESS CREEP VALUE:
! 
      DO 100 K=1,NROWB
!         TENSAO(KK) = STRSS(NEL,1,KK)
         TENSAO(K) = AVSTRS(NEL,K)
 100  CONTINUE
!
      IF (REGION(NEL).EQ.1) QTRIAL = -10.0D0
!      IF (REGION(NEL).GE.2) QTRIAL = 0.0D0
      IF (REGION(NEL).GE.2) THEN 
!.... ... COMPUTE DEVIATOR STRESS TENSOR
         CALL DEVTENSOR(DEVSTRS,TRCSTRS,TENSAO,NROWB)
!.... ... COMPUTE EFFECTIVE TRIAL STRESS
         QTRIAL=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
      ENDIF
!
      WRITE(ICREP,4000) QTRIAL
!
!.... PRINT SATURATION
!
      IF (REGION(NEL).EQ.1) STRIAL =  S(NEL)
      IF (REGION(NEL).GT.1) STRIAL = -10.0D0
      WRITE(ISATR,3000) STRIAL
!
!.... PRINT PERMEABILITY
!
      IF (REGION(NEL).EQ.1) STRIAL = PERM(NEL)
      IF (REGION(NEL).GT.1) STRIAL = -10.0D0
!      WRITE(IPERM,4000) XKKC(PORE(NEL))*PERM(NEL)
      WRITE(IPERM,4000) STRIAL
!
!.... PRINT TOTAL VELOCITY 
!
      IF (REGION(NEL).EQ.1) THEN 
             STRIAL = VC(1,NEL)
             SPORE  = VC(2,NEL)
         ELSE
             STRIAL = 0.0D0
             SPORE  = 0.0D0
      ENDIF
!       WRITE(IVELT,4000) (VC(J,NEL),J=1,NDOF2)
      WRITE(IVELT,4000)  STRIAL, SPORE
!      WRITE(IVELT,4000) (STRSS(NEL,4,K),K=1,4)
!
!.... PRINT STRESS_X COMPONENT
! 
      WRITE(ISIGX,4000) AVSTRS(NEL,1)
!
!.... PRINT STRESS_Y COMPONENT
      IF (REGION(NEL).EQ.1) QTRIAL = -9.0D20
!      IF (REGION(NEL).GE.2) QTRIAL = 0.0D0
      IF (REGION(NEL).GE.2) QTRIAL =  AVSTRS(NEL,2)
!
      WRITE(ISIGY,4000) qtrial   !AVSTRS(NEL,2)
!
!.... PRINT STRESS_XY COMPONENT 
!
      WRITE(ISGTA,4000) AVSTRS(NEL,3)
!
!.... PRINT STRESS_Z COMPONENT 
!
      WRITE(ISIGZ,4000) AVSTRS(NEL,2)
!
!.... PRINT YOUNG MODULUS 
!
!      WRITE(IYUNG,4000) dfloat(REGION(NEL))
      WRITE(IYUNG,4000) YOUNG(NEL)
!
!.... PRINT ELEMENT MASS BALANCE 
!
      IF (REGION(NEL).EQ.1) SBLNC = BALANC(NEL)
      IF (REGION(NEL).GT.1) SBLNC = 0.0D0
      WRITE(IBLNC,4000) SBLNC
!
 500  CONTINUE
! C
! c      IF (NNP.EQ.NVEL) THEN
! c      PI=4.0D0*DATAN(1.0D0)
! c      DO 600 NEL=51,NUMEL,100
! c         radio=DSQRT(XC(1,NEL)**2+XC(2,NEL)**2)
! c         angle=DATAN(XC(2,NEL)/XC(1,NEL))
! 
! C         sgxxcos2 = (dcos(angle))**2*STRSS(NEL,1,1)
! C         sgyysin2 = (dsin(angle))**2*STRSS(NEL,1,2)
! C         sgxysin2 = (dsin(2.0d0*angle))*STRSS(NEL,1,3)
! c         sgxxcos2 = (dcos(angle))**2*AVSTRS(NEL,1)
! c         sgyysin2 = (dsin(angle))**2*AVSTRS(NEL,2)
! c         sgxysin2 = (dsin(2.0d0*angle))*AVSTRS(NEL,3)
! C
! c         sigmar = sgxxcos2+sgyysin2+sgxysin2 
! C
! c         WRITE(101,4000) RADIO,ANGLE*180.0D0/PI,SIGMAR
 600  CONTINUE
!      ENDIF
      CLOSE(IDISP)
      CLOSE(IPRSR)
      CLOSE(IPORE)
      CLOSE(ICREP)
      CLOSE(ISATR)
      CLOSE(IPERM)
      CLOSE(ISIGX)
      CLOSE(ISIGY)
      CLOSE(ISIGZ)
      CLOSE(IYUNG)
      CLOSE(IVELT)
!
!.... PRINT INFORMATION ON NODAL DX FILE "nodestoc.dx"
!
      CALL PRINT_DXINFO('WRITE_FEDX_DATA',NUMNP,NUMEL)
!
      RETURN 
!
 2000 FORMAT(I5,2X,40(1PE15.8,2X)) 
 3000 FORMAT(2X,5(F25.15,2x)) 
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE

!
!**** NEW **** SUBROUTINE BASED ON FOR GEOMECHANIC  ************* 
! 
      SUBROUTINE RSRVR_MASS(nsd,numel,nedfl,numedg,nedg,ieedg, &
     &           fluxo,temp,IEN,VDP,NEN,NDOF2,NED2,NUMNP)
!
!...  COMPUTE GLOBAL MASS BALANCE OVER ALL THE RESERVOIR
!
!...  FRAMEWORK: GAUSS THEOREM APPLIED ON RESERVOIR GLOBAL MASS EQUATION 
!
      USE mLeituraEscrita, only: IRBAL
      use mMalha,          only: local
      use mPropGeoFisica,  only: NELX
!
      IMPLICIT NONE
!
      integer :: numel,nsd,nedfl,numedg,nedg
      integer :: n2,n3,n4,j
      real(8) :: temp
      real(8), dimension(nedfl,*) :: fluxo
      integer, dimension(nedg,*)  :: ieedg
      real(8), dimension(4)       :: flux,xindj
      real(8), dimension(2,4)     :: vnl
!      
      integer :: nel
!
!... BEGIN CATALOG OF GEOMECHANIC 
!
!
      INTEGER, DIMENSION(NEN,*)   :: IEN
      REAL(8), DIMENSION(NDOF2,*) :: VDP
      REAL(8), DIMENSION(NED2,NEN) :: VDPL
      INTEGER :: NEN, NDOF2, NED2, NUMNP
      REAL(8), DIMENSION(4) :: VDPEDG
      REAL(8) :: XXINJC, XXDFRM, XXOUTF
!
!... END CATALOG OF GEOMECHANIC
!
      XXINJC = 0.0D0
      XXDFRM = 0.0D0
      XXOUTF = 0.0D0
!
      vnl=0.d0
      vnl(2,1)=-1.d0
      vnl(1,2)= 1.d0
      vnl(2,3)= 1.d0
      vnl(1,4)=-1.d0
!
      do j=1,nedg
         xindj(j) = vnl(1,j)+vnl(2,j)
      end do
!         
!     EDGE LOCAL ORDER
!
!               3
!            _______
!           |       |
!         4 |       | 2
!           |_______|
!
!               1
!..
!...  INJECTION BOUNDARY
!..
      DO 100 NEL=1,NUMEL,NELX
!
!         n1=ieedg(1,nel)
!         n2=ieedg(2,nel)
!         n3=ieedg(3,nel)
!         n4=ieedg(4,nel)
         N4 = IEEDG(4,NEL)
!     
!....LOCALIZE SOLID'S NODAL VELOCITIES ON ELEMENT
!
       CALL LOCAL(IEN(1,NEL),VDP,VDPL,NEN,NDOF2,NED2) 
!
!.... FROM LOCAL DISPLACEMENT VELOCITIES (VDPL[NDOF2,NEN]) 
!....  TO  LOCAL EDGES DISPLACEMENT VELOCITIES (VDPEDG[EDGE])
!
!       VDPEDG(1)=(VDPL(1,1)+VDPL(1,4))*0.5D0
!       VDPEDG(2)=(VDPL(2,1)+VDPL(2,2))*0.5D0
!       VDPEDG(3)=(VDPL(1,2)+VDPL(1,3))*0.5D0
!       VDPEDG(4)=(VDPL(2,3)+VDPL(2,4))*0.5D0
       VDPEDG(4)=(VDPL(1,1)+VDPL(1,4))*0.5D0
!
!       FLUX(1)=XINDJ(1)*fluxo(1,n1)+VDPEDG(1)
!       FLUX(2)=XINDJ(2)*fluxo(1,n2)+VDPEDG(2)
!       FLUX(3)=XINDJ(3)*fluxo(1,n3)+VDPEDG(3)
!       FLUX(4)=XINDJ(4)*fluxo(1,n4)+VDPEDG(4)
       XXINJC = XXINJC + XINDJ(4)*FLUXO(1,N4)+VDPEDG(4)
!
 100    CONTINUE
!
!     EDGE LOCAL ORDER
!
!               3
!            _______
!           |       |
!         4 |       | 2
!           |_______|
!
!               1
!..
!...  UPPER BOUNDARY
!..
      DO 200 NEL=NUMEL-NELX+1,NUMEL
!
!         n1=ieedg(1,nel)
!         n2=ieedg(2,nel)
!         n3=ieedg(3,nel)
!         n4=ieedg(4,nel)
         N3 = IEEDG(3,NEL)
!     
!....LOCALIZE SOLID'S NODAL VELOCITIES ON ELEMENT
!
       CALL LOCAL(IEN(1,NEL),VDP,VDPL,NEN,NDOF2,NED2) 

!
!.... FROM LOCAL DISPLACEMENT VELOCITIES (VDPL[NDOF2,NEN]) 
!....  TO  LOCAL EDGES DISPLACEMENT VELOCITIES (VDPEDG[EDGE])
!
!       VDPEDG(1)=(VDPL(1,1)+VDPL(1,4))*0.5D0
!       VDPEDG(2)=(VDPL(2,1)+VDPL(2,2))*0.5D0
!       VDPEDG(3)=(VDPL(1,2)+VDPL(1,3))*0.5D0
!       VDPEDG(4)=(VDPL(2,3)+VDPL(2,4))*0.5D0
       VDPEDG(3)=(VDPL(2,3)+VDPL(2,4))*0.5D0


!      WRITE(IMASL,1000) NEL, TEMP,VDPEDG(3)
!
!       FLUX(1)=XINDJ(1)*fluxo(1,n1)+VDPEDG(1)
!       FLUX(2)=XINDJ(2)*fluxo(1,n2)+VDPEDG(2)
!       FLUX(3)=XINDJ(3)*fluxo(1,n3)+VDPEDG(3)
!       FLUX(4)=XINDJ(4)*fluxo(1,n4)+VDPEDG(4)
        XXDFRM = XXDFRM + XINDJ(3)*FLUXO(1,N3)+VDPEDG(3)
!
 200    CONTINUE
!
!     EDGE LOCAL ORDER
!
!               3
!            _______
!           |       |
!         4 |       | 2
!           |_______|
!
!               1
!..
!...  OUT-FLOW BOUNDARY
!..
      DO 300 NEL=NELX,NUMEL,NELX
!
!         n1=ieedg(1,nel)
!         n2=ieedg(2,nel)
!         n3=ieedg(3,nel)
!         n4=ieedg(4,nel)
         N2 = IEEDG(2,NEL)
!     
!....LOCALIZE SOLID'S NODAL VELOCITIES ON ELEMENT
!
       CALL LOCAL(IEN(1,NEL),VDP,VDPL,NEN,NDOF2,NED2) 
!
!.... FROM LOCAL DISPLACEMENT VELOCITIES (VDPL[NDOF2,NEN]) 
!....  TO  LOCAL EDGES DISPLACEMENT VELOCITIES (VDPEDG[EDGE])
!
!       VDPEDG(1)=(VDPL(1,1)+VDPL(1,4))*0.5D0
!       VDPEDG(2)=(VDPL(2,1)+VDPL(2,2))*0.5D0
!       VDPEDG(3)=(VDPL(1,2)+VDPL(1,3))*0.5D0
!       VDPEDG(4)=(VDPL(2,3)+VDPL(2,4))*0.5D0
       VDPEDG(2)=(VDPL(1,2)+VDPL(1,3))*0.5D0
!
!       FLUX(1)=XINDJ(1)*fluxo(1,n1)+VDPEDG(1)
!       FLUX(2)=XINDJ(2)*fluxo(1,n2)+VDPEDG(2)
!       FLUX(3)=XINDJ(3)*fluxo(1,n3)+VDPEDG(3)
!       FLUX(4)=XINDJ(4)*fluxo(1,n4)+VDPEDG(4)
        XXOUTF = XXOUTF + XINDJ(2)*FLUXO(1,N2)+VDPEDG(2)
!
 300    CONTINUE
!
!      WRITE(IMASL,2000) TEMP, XXINJC, XXDFRM, XXOUTF, xxinjc+xxdfrm
      WRITE(IRBAL,2000) TEMP, XXINJC, XXDFRM, XXOUTF, xxinjc+xxdfrm
!
      RETURN
!
 1000 FORMAT(I7,2X,7(1PE15.8,2X))
!
 2000 FORMAT(2X,7(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW **** FORCES APPLIED ON IRREGULAR MESH **************************
!
      SUBROUTINE F4IRRMSH(FDIS,X,NSD,NDOF2,NUMNP,NLVECT,NELX,PRSRINIT)
!
!**** PROGRAM TO SETUP NODAL FORCES OVER IRREGULAR MESH 
!
      IMPLICIT NONE
!
      INTEGER       :: NODE, NREF, NLAST
      INTEGER       :: NSD, NDOF2, NUMNP, NLVECT, NELX
      REAL(8), DIMENSION(NDOF2,NUMNP,NLVECT) :: FDIS
      REAL(8), DIMENSION(NSD,NUMNP)          :: X
!
      REAL(8)       :: LENLEFT, LENRIGHT, PRSRINIT
!
      NREF  = NUMNP-NELX
      NLAST = NUMNP-1
!
      LENRIGHT        =  0.5D0*(X(1,NREF+1) - X(1,NREF))
      FDIS(2,NREF,1)  = -PRSRINIT*LENRIGHT
      LENLEFT         =  0.5D0*(X(1,NUMNP) - X(1,NLAST))
      FDIS(2,NUMNP,1) = -PRSRINIT*LENLEFT
!
      NREF  = NREF+1
!
      DO 500 NODE=NREF,NLAST
         LENLEFT     = 0.5D0*(X(1,NODE)- X(1,NODE-1))
         LENRIGHT    = 0.5D0*(X(1,NODE+1) - X(1,NODE))
         FDIS(2,NODE,1) = -PRSRINIT*(LENLEFT+LENRIGHT)
 500  CONTINUE
!
      RETURN
!  
      END SUBROUTINE
!
!**** NEW **** FORCES APPLIED ON IRREGULAR MESH **************************
!
      SUBROUTINE F4IRRMSH2(FDIS,X,NSD,NDOF2,NUMNP,NLVECT,NELX,NUMNPR,PRSRINIT)
!
!**** PROGRAM TO SETUP NODAL FORCES OVER IRREGULAR MESH 
!
      IMPLICIT NONE
!
      INTEGER       :: NODE, NREF, NLAST, NTOPLEFT, NTOPCENT, NTOPRIGHT
      INTEGER       :: NSD, NDOF2, NUMNP, NLVECT, NELX, NUMNPR
      REAL(8), DIMENSION(NDOF2,NUMNP,NLVECT) :: FDIS
      REAL(8), DIMENSION(NSD,NUMNP)          :: X
!
      REAL(8) :: LENLEFT, LENRIGHT, LENTOP, LENBOTTOM, PRSRINIT
      REAL(8) :: XLEFT,XRIGHT,YTOP,YBOTTOM 
!
      XLEFT   = 0.0D0
      YTOP    = 0.0D0
      XRIGHT  = 0.0D0
      YBOTTOM = 0.0D0
!
      NTOPLEFT  = 1
      NTOPCENT  = 1
      NTOPRIGHT = 1
!
      DO 100 NODE=NUMNPR,NUMNP
         XLEFT   = MIN(X(1,NODE),XLEFT)
         XRIGHT  = MAX(X(1,NODE),XRIGHT)
         YTOP    = MAX(X(2,NODE),YTOP)
         YBOTTOM = MIN(X(2,NODE),YBOTTOM)
100   CONTINUE
!
      WRITE(*,*) 'XLEFT = ',XLEFT,'XRIGHT  = ',XRIGHT
      WRITE(*,*) 'YTOP  = ', YTOP,'YBOTTOM = ',YBOTTOM 
!
!... IF-THEN-ELSE: FOR SMALL CAPS DOMAINS
!
      IF ((XLEFT.EQ.0.0D0).OR.(XLEFT.EQ.-540.0D0)) THEN 
         NREF  = NUMNP-NELX
         NLAST = NUMNP-1
!
         LENRIGHT        =  0.5D0*(X(1,NREF+1) - X(1,NREF))
         FDIS(2,NREF,1)  = -PRSRINIT*LENRIGHT
         LENLEFT         =  0.5D0*(X(1,NUMNP) - X(1,NLAST))
         FDIS(2,NUMNP,1) = -PRSRINIT*LENLEFT
!
         NREF  = NREF+1
!
         DO 120 NODE=NREF,NLAST
            LENLEFT        = 0.5D0*(X(1,NODE)- X(1,NODE-1))
            LENRIGHT       = 0.5D0*(X(1,NODE+1) - X(1,NODE))
            FDIS(2,NODE,1) = -PRSRINIT*(LENLEFT+LENRIGHT)
 120    CONTINUE
!
        RETURN
!
      ENDIF

      DO 200 NODE=NUMNPR,NUMNP
         LENLEFT   = XLEFT   - X(1,NODE)
         LENRIGHT  = XRIGHT  - X(1,NODE)
         LENTOP    = YTOP    - X(2,NODE)
         LENBOTTOM = YBOTTOM - X(2,NODE)
!
         IF ((LENLEFT  .EQ.0.0D0).AND.(LENTOP.EQ.0.0D0)) NTOPLEFT  = NODE
         IF ((LENRIGHT .EQ.0.0D0).AND.(LENTOP.EQ.0.0D0)) NTOPRIGHT = NODE 
         IF ((X(1,NODE).EQ.0.0D0).AND.(LENTOP.EQ.0.0D0)) NTOPCENT  = NODE-2
200   CONTINUE
!
      WRITE(*,*) 'NODE LEFT=',NTOPLEFT,'NODE CENT=',NTOPCENT,'NODE RIGHT=',NTOPRIGHT

      write(*,*) 'point left=',X(1,NTOPLEFT),X(2,NTOPLEFT)
      write(*,*) 'point CENT=',X(1,NTOPCENT),X(2,NTOPCENT)
!      write(*,*) 'point resr=',X(1,NTOPCENT+2),X(2,NTOPCENT+2)
      write(*,*) 'point RIGH=',X(1,NTOPRIGHT),X(2,NTOPRIGHT)
      write(*,*) 'point r+1=',X(1,NTOPright+1),X(2,NTOPright+1)
!
!      IF ((XLEFT.EQ.0.0D0).OR.(XLEFT.EQ.-540.0D0)) THEN 
!         STOP
!      ENDIF
!
!.... FORCE ON TOP LEFT POINT
      LENRIGHT            =  0.5D0*(DABS(X(1,NTOPLEFT-1) - X(1,NTOPLEFT)))
      FDIS(2,NTOPLEFT,1)  = -PRSRINIT*LENRIGHT
!.... FORCE ON TOP RIGHT POINT
      LENLEFT             =  0.5D0*(DABS(X(1,NTOPRIGHT) - X(1,NTOPRIGHT-1)))
      FDIS(2,NTOPRIGHT,1) = -PRSRINIT*LENLEFT
!.... FORCE ON TOP LEFT CENTER POINT
      LENRIGHT              =  0.5D0*(DABS(X(1,NTOPRIGHT+1) - X(1,NTOPCENT)))
      LENLEFT               =  0.5D0*(DABS(X(1,NTOPRIGHT+2) - X(1,NTOPRIGHT+1)))
      FDIS(2,NTOPRIGHT+1,1) = -PRSRINIT*(LENRIGHT+LENLEFT)
!.... FORCE ON TOP RIGHT POINT
      LENRIGHT            =  0.5D0*(DABS(X(1,NTOPCENT+1) - X(1,NTOPCENT)))
      LENLEFT             =  0.5D0*(DABS(X(1,NTOPRIGHT+1) - X(1,NTOPCENT)))
      FDIS(2,NTOPCENT,1)  = -PRSRINIT*(LENRIGHT+LENLEFT)
!.... FORCE ON RIGHT TOP CAP SIDE
      NREF  = NTOPCENT+1
      NLAST = NTOPRIGHT-1
!
      DO 300 NODE=NREF,NLAST
         LENLEFT        = 0.5D0*(DABS(X(1,NODE)- X(1,NODE-1)))
         LENRIGHT       = 0.5D0*(DABS(X(1,NODE+1) - X(1,NODE)))
         FDIS(2,NODE,1) = -PRSRINIT*(LENLEFT+LENRIGHT)
 300  CONTINUE
!.... FORCE ON LEFT TOP CAP SIDE
      NREF  = NTOPRIGHT+2
      NLAST = NTOPLEFT-1
!
      DO 400 NODE=NREF,NLAST
         LENLEFT        = 0.5D0*(DABS(X(1,NODE)- X(1,NODE-1)))
         LENRIGHT       = 0.5D0*(DABS(X(1,NODE+1) - X(1,NODE)))
         FDIS(2,NODE,1) = -PRSRINIT*(LENLEFT+LENRIGHT)
 400  CONTINUE
!
      RETURN
!  
      END SUBROUTINE
! 
!**** NEW **** MODIFIED SUBROUTINE FOR GEOMECHANIC  ********* 
! 
      SUBROUTINE GLOBAL_MASS(nsd,numel,nedfl,numedg,nedg,ieedg,&
     & fluxo,IEN,NEN,NED2,NUMNP)
!
!...  PROGRAM TO COMPUTE GLOBAL MASS BALANCE ON EVERY MESH ELEMENT 
!
!...  FRAMEWORK: GAUSS THEOREM APPLIED ON EACH ELEMENT
!
      use mPropGeoFisica, only: hx, hy
      use mGlobaisEscalares, only: NDOFD
      use mMalha, only: local
!
      implicit none
!
      integer :: numel,nsd,nedfl,numedg,nedg
      integer :: n1,n2,n3,n4,j
      real(8), dimension(1,*)     :: fluxo
      integer, dimension(nedg,*)  :: ieedg
      real(8), dimension(4)       :: flux,xindj
      real(8), dimension(2,4)     :: vnl
!      
      integer :: nel
!
!... BEGIN CATALOG OF GEOMECHANIC 
!
      INTEGER, DIMENSION(NEN,*)   :: IEN
      INTEGER :: NEN, NDOF2, NED2, NUMNP
      REAL(8), DIMENSION(4) :: VDPEDG
      REAL(8) :: TFLUX
      REAL(8), DIMENSION(NED2,NEN) :: VDPL
!
!.... MODIFICATION FOR GEOMECHANICS IN CAPITAL LETTERS
!
!... END CATALOG OF GEOMECHANIC
!
      vnl=0.d0
      vnl(2,1)=-1.d0
      vnl(1,2)= 1.d0
      vnl(2,3)= 1.d0
      vnl(1,4)=-1.d0
!
      do j=1,nedg
         xindj(j) = vnl(1,j)+vnl(2,j)
      end do
!
      do nel=1,numel
!
         n1=ieedg(1,nel)
         n2=ieedg(2,nel)
         n3=ieedg(3,nel)
         n4=ieedg(4,nel)
!         
!... EDGE LOCAL ORDER 
!
!               3
!            _______
!           |       |
!         4 |       | 2
!           |_______|
!
!               1
!
!     
!....LOCALIZE SOLID'S NODAL VELOCITIES ON ELEMENT
!
       CALL LOCAL(IEN(1,NEL),VDP,VDPL,NEN,NDOFD,NED2) 
!
!.... FROM LOCAL DISPLACEMENT VELOCITIES (VDPL[NDOF2,NEN]) 
!....  TO  LOCAL EDGES DISPLACEMENT VELOCITIES (VDPEDG[EDGE])
!
       VDPEDG(1)=(VDPL(2,1)+VDPL(2,2))*0.5D0
       VDPEDG(2)=(VDPL(1,2)+VDPL(1,3))*0.5D0
       VDPEDG(3)=(VDPL(2,3)+VDPL(2,4))*0.5D0
       VDPEDG(4)=(VDPL(1,1)+VDPL(1,4))*0.5D0

!         flux(1)=xindj(1)*f(s(nel),fluxo(1,n1))
!         flux(2)=xindj(2)*g(s(nel),fluxo(1,n2))
!         flux(3)=xindj(3)*f(s(nel),fluxo(1,n3))
!         flux(4)=xindj(4)*g(s(nel),fluxo(1,n4))
!         tflux  =flux(1)+flux(3)+flux(4)+flux(2)
!
       FLUX(1)=XINDJ(1)*FLUXO(1,n1)+VDPEDG(1)
       FLUX(2)=XINDJ(2)*FLUXO(1,n2)+VDPEDG(2)
       FLUX(3)=XINDJ(3)*FLUXO(1,n3)+VDPEDG(3)
       FLUX(4)=XINDJ(4)*FLUXO(1,n4)+VDPEDG(4)
!
       TFLUX  = FLUX(1)+FLUX(3)+FLUX(4)+FLUX(2)
!
       BALANC(NEL) = TFLUX
!
      END DO
!
      END SUBROUTINE


!
!**** NEW **** FOR INCOMPRESSIBILITY ***************************************** 
!
      subroutine VECTRSCR_LIN(x, conecNodaisElem,p, brhsd, lmD)
!
      use mGlobaisEscalares, only: ndofD, nrowsh
      use mAlgMatricial,     only: rowdot, coldot, addrhs, neqD
      use mFuncoesDeForma,   only: shgq,shlq
      use mMalha,            only: local, multab, numel, nen, nsd, numelReserv
      use mPropGeoFisica,    only: POISBTTM, YNGBOTTM, BULKSOLID
!
!.... PROGRAM TO COMPUTE PLASTIC DEFORMATION OVER THE MACRO DOMAIN 
! 
      IMPLICIT NONE
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL DIAG,QUAD,ZERODL
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      INTEGER :: conecNodaisElem(NEN,*),  LMD(NED2,NEN,*)

      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL*8 :: X(NSD,*), ELEFFMD(NEE2,NEE2),ELRESFD(NEE2), &
     &          BRHSD(*), &
     &          SHGBR(NROWSH,NEN), P(NUMELRESERV)

!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL*8  :: BBARI(NROWB,NESD), PRESSURE(NROWB)
      real(8) :: xl(nesd,nen), disl(ned2,nen), c1
      INTEGER :: NEL, I, L, II, JJ
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)

      REAL(8) :: BULKROCK, BIOTCOEF 
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!... CONSISTENT MATRIX
!
      DIAG = .FALSE.

      DO 500 NEL=1,NUMEL

      BULKROCK = BULK(YNGBOTTM,POISBTTM)
!
      BIOTCOEF = 1.0D0 - BULKROCK/BULKSOLID
!       BIOTCOEF = 1.d0
!      
!... CLEAR STIFFNESS MATRIX AND FORCE ARRAY
!
      ELEFFMD=0.D0
      ELRESFD=0.D0
!
!... SETUP ELEMENT PRESSURE VECTOR
!
      IF (NEL.GT.numelReserv) THEN
         PRESSURE=0.D0
      ELSE
         PRESSURE(1)=BIOTCOEF*P(NEL)
         PRESSURE(2)=BIOTCOEF*P(NEL)
         PRESSURE(3)=0.0D0
         PRESSURE(4)=BIOTCOEF*P(NEL) 
      ENDIF

!
!... LOCALIZE COORDINATES AND DIRICHLET B.C.
!
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
      QUAD = .TRUE.
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
      CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
        DO 100 L=1,NINTD
           R(L)   = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
           DETD(L) = DETD(L)*R(L)
 100    CONTINUE
      ENDIF
!
!.... FORM STIFFNESS MATRIX 
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... .. FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... .. LOOP OVER INTEGRATIONN POINTS
!
         DO 400 L=1,NINTD
!
         C1=DETD(L)*WD(L)
!
!**** **** MOUNT FORCE VECTOR ****************** 
!
         CALL CLEAR(BBARI,NROWB*NESD)
!
         DO 300 I=1,NEN
!
         CALL SETBB(BBARI,SHGD(1:NROWSH,I,L),&
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
          ELRESFD(NED2*I-1)= ELRESFD(NED2*I-1) &
     &                       + COLDOT(BBARI(1:4,1),PRESSURE,4)*C1
!
          ELRESFD(NED2*I)  = ELRESFD(NED2*I) &
     &                       + COLDOT(BBARI(1:4,2),PRESSURE,4)*C1
!
  300    CONTINUE
  400    CONTINUE


!
       DO 450 II=1,NEE2
       DO 450 JJ=1,NEE2
!
         ELEFFMD(II,JJ) = AUXM(NEL,II,JJ)
!
 450  CONTINUE
!
!     COMPUTATION OF DIRICHLET B.C. CONTRIBUTION
!   
      CALL ZTEST(DISL,NEE2,ZERODL)
!
       IF(.NOT.ZERODL) &
     & CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL)
!
!.... ASSEMBLE ELEMENT FORCE ARRAY INTO GLOBAL RIGHT-HAND SIDE VECTOR
!
      CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2)
!         write(766,*) nel, ELRESFD(1:nee2)
!
 500  CONTINUE
! 

!         stop "vectrscr"
      RETURN
!
      END SUBROUTINE

!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      subroutine POS4TRNS1(x, conecNodaisElem, p, p0)
! 
      use mMalha, only:nsd, numnp, numel, nen, LOCAL, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh
      use mAlgMatricial, only: rowdot, coldot
      use mFuncoesDeForma, only: shgq, shlq
      use mPropGeoFisica, only: YOUNG, POISSCAP, PORE, PORE0, POISBTTM, POISMIDL, POISTOP, REGION, BULKSOLID
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
!
      IMPLICIT NONE
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numel), p0(numel)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL QUAD
      real(8) :: xl(nesd,nen), disl(ned2,nen), AREA, xdivu, c1
      REAL*8  :: POISSON
      integer :: nel, l, j, k
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: DIVULOC(NUMEL)
      REAL(8) :: UNITVEC(NROWB),BBARJ(NROWB,NESD)
      REAL(8) :: CBBAR(NROWB, NROWB)
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD),WD(NINTD)
!
      REAL*8 :: STRAIN(NROWB),DEVSTRS(NROWB)
      REAL*8, DIMENSION(NROWB) :: TENSAO
      REAL*8 :: ROOT3D2, TRCSTRS, QVM

      REAL(8) :: POROEXP, BIOTMOD, BULKROCK, BIOTCOEF, DEFNM1,DIFFDIVU,DIFFPRES

      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0  
      ROOT3D2=1.224744871391589D0

!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)

      POISSON = POISBTTM
!
      DO 500 NEL=1,NUMEL
!
         IF (REGION(NEL).EQ.2) POISSON = POISMIDL
         IF (REGION(NEL).EQ.3) POISSON = POISTOP
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
        CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
        DIVULOC(NEL) = 0.0D0
!..
!.... ..CLEAR INITIAL STRAIN
!..
        CALL CLEAR(STRAIN,NROWB) 
        CALL CLEAR(TENSAO,NROWB) 
!..
!.... ..DEFINE ELEMENT AREA
!..
         AREA = 0.0D0
!..
!.... ..SETUP FOR AXISYMMETRIC OPTION
!..
        IF (IOPT.EQ.2) THEN
          DO 10 L=1,NINTD
             R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
             DETD(L) = DETD(L)*R(L) 
 10       CONTINUE
        ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
        DO 300 L=1,NINTD
!
        C1=WD(L)*DETD(L)
        AREA = AREA + C1
!
        DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
        CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &             R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
        DO 200 K=1,NROWB
           STRAIN(K)=STRAIN(K)+                                &
     &                COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200    CONTINUE
!..
!.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION 
!..
!.... ..TO COMPUTE MEAN VOLUMETRIC CREEP 
!.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
!
!.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
!
            QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!
        DO 220 K=1,NROWB
           ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+                    &
     &                     ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
 220    CONTINUE
!
        DO 250 K=1,NROWB
           TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
 250    CONTINUE
!
 300    CONTINUE

!
        DO 350 K=1,NROWB       
           STRAIN(K)=STRAIN(K)/AREA
           AVSTRS(NEL,K)=TENSAO(K)/AREA
 350    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT 
        DIVULOC(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)

!        write(888,*) divuloc(nel)
 500  CONTINUE
!..
      DO 700 NEL=1,NUMEL
         DIVU(NEL)=DIVULOC(NEL)
700   CONTINUE

!      DO 800 NEL=1,NUMEL
!.... ..COMPUTE MEAN VOLUMETRIC DEFORMATION FOR RESERVOIR
!..
!        IF (REGION(NEL).EQ.1) THEN

!           DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
!           XDIVU = DIVU(NEL)-DIVU0(NEL)
! !.... ..COMPUTE POROSITY
!           PORE(NEL)=1.0D0-(1.0D0-PORE0(NEL))*DEXP(-XDIVU)

!..
!.... ..COMPUTE BULK MODULUS AND BIOT COEFICIENT(ALPHA)
!..
!           BULKROCK  = BULK(YOUNG(NEL),POISBTTM)
!           BIOTCOEF  = 1.0D0 - BULKROCK/BULKSOLID
!..
!..new line: linearized See: Coussy 2.Edt. Eqs.4.61 
!..
!           DEFNM1    = BIOTMOD(PORE0(NEL),BIOTCOEF)-PORE0(NEL)/BULKOIL
!           DEFNM1    = (BIOTCOEF-PORE0(NEL))/BULKSOLID
!
!           DIFFDIVU  = DIVULOC(NEL)-DIVU0(NEL)
!           DIFFPRES  = P(1,NEL)-P0(NEL)
!
!.... ..COMPUTE POROSITY
!old line  PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-XDIVU)
!
!new line: linearized See: Coussy 2.Edt. Eqs. 4.19
!
!           PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES
!        ENDIF
!..
!800   CONTINUE
!
!      DO 550 NEL=1,NELX1*NELY1
!         DO 550 JJ=1,NROWB
!            STRSS(JJ,NEL)=0.0D0
! 550  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      subroutine POS4TRNS2(x, conecNodaisElem, p, p0)
! 
      use mMalha, only:nsd, numnp, numel, nen, LOCAL, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh
      use mAlgMatricial, only: rowdot, coldot
      use mFuncoesDeForma, only: shgq, shlq
      use mPropGeoFisica, only: YOUNG, POISSCAP, PORE, PORE0, POISBTTM, POISMIDL, POISTOP, REGION, BULKSOLID
      use mPropGeoFisica, only: BULKSOLID, BULKWATER, MASCN, MASCN0, PHIEULER
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
!
      IMPLICIT NONE
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numel), p0(numel)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL QUAD
      real(8) :: xl(nesd,nen), disl(ned2,nen), AREA, xdivu, c1
      REAL*8  :: POISSON
      integer :: nel, l, j, k
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: DIVULOC(NUMEL)
      REAL(8) :: UNITVEC(NROWB),BBARJ(NROWB,NESD)
      REAL(8) :: CBBAR(NROWB, NROWB)
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD),WD(NINTD)
!
      REAL*8 :: STRAIN(NROWB),DEVSTRS(NROWB)
      REAL*8, DIMENSION(NROWB) :: TENSAO
      REAL*8 :: ROOT3D2, TRCSTRS, QVM, JACOBIAN

      REAL(8) :: BIOTMOD, BULKROCK, BIOTCOEF, DEFNM1, DEFMM1, DIFFDIVU, DIFFPRES

      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0  
      ROOT3D2=1.224744871391589D0

!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)

      POISSON = POISBTTM
!
      DO 500 NEL=1,NUMEL
!
         IF (REGION(NEL).EQ.2) POISSON = POISMIDL
         IF (REGION(NEL).EQ.3) POISSON = POISTOP
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
        CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
        DIVU(NEL) = 0.0D0
!..
!.... ..CLEAR INITIAL STRAIN
!..
        CALL CLEAR(STRAIN,NROWB) 
        CALL CLEAR(TENSAO,NROWB) 
!..
!.... ..DEFINE ELEMENT AREA
!..
         AREA = 0.0D0
!..
!.... ..SETUP FOR AXISYMMETRIC OPTION
!..
        IF (IOPT.EQ.2) THEN
          DO 10 L=1,NINTD
             R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
             DETD(L) = DETD(L)*R(L) 
 10       CONTINUE
        ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
        DO 300 L=1,NINTD
!
        C1=WD(L)*DETD(L)
        AREA = AREA + C1
!
        DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
        CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &             R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
        DO 200 K=1,NROWB
           STRAIN(K)=STRAIN(K)+                                &
     &                COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200    CONTINUE
!..
!.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION 
!..
!.... ..TO COMPUTE MEAN VOLUMETRIC CREEP 
!.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
!
!.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
!
            QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!
        DO 220 K=1,NROWB
           ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+                    &
     &                     ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
 220    CONTINUE
!
        DO 250 K=1,NROWB
           TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
 250    CONTINUE
!
 300    CONTINUE

!
        DO 350 K=1,NROWB       
           STRAIN(K)=STRAIN(K)/AREA
           AVSTRS(NEL,K)=TENSAO(K)/AREA
 350    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT 
        DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)

!        write(888,*) divuloc(nel)
 500  CONTINUE
!..
!      DO 700 NEL=1,NUMEL
!         DIVU(NEL)=DIVULOC(NEL)
!700   CONTINUE

      DO 800 NEL=1,NUMELRESERV
!.... ..COMPUTE MEAN VOLUMETRIC DEFORMATION FOR RESERVOIR
!..
           BULKROCK  = BULK(YOUNG(NEL),POISBTTM)

           BIOTCOEF  = 1.0D0 - BULKROCK/BULKSOLID   !also called ALPHA
!
!... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35 
!
           DEFNM1    = (BIOTCOEF-PORE0(NEL))/BULKSOLID
!
!... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
!
           DEFMM1 = DEFNM1 + PORE0(NEL)/BULKWATER 
!
!... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27 
!
           JACOBIAN  = 1.0D0 + DIVU(NEL)
!
           DIFFDIVU  = DIVU(NEL)-DIVU0(NEL)
           DIFFPRES  = P(1,NEL)-P0(NEL)
!
! ..COMPUTE LAGRANGIAN POROSITY
!           PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-diffdivu) 
!
!new line: linearized See: Coussy 2.Edt. Eqs. 4.19
           PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES
!
!...COMPUTE EULERAIN POROSITY
           PHIEULER(NEL) = PORE(NEL)/JACOBIAN

!
!...COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
!
           MASCN(NEL)= MASCN0(NEL)+(BIOTCOEF*DIFFDIVU+DEFMM1*DIFFPRES)

           if (mascn(nel).lt.0.0d0) then
               write(*,*) 'nel=',nel
!               write(*,*) 'RHOEQUIV=',RHOEQUIV
                write(*,*) 'alfa=',BIOTCOEF
               write(*,*) 'divu=',divu(nel)
               write(*,*) 'divu0=',divu0(nel)
               write(*,*) 'biotmod=',defmm1
               write(*,*) 'biotmod*diffDivu=',biotmod*diffdivu
               write(*,*) 'diffpress=',diffpres
               write(*,*) '1/m(p-p0)=',defmm1*DIFFPRES
               write(*,*) 'masscont=',pore(nel), pore0(nel) 
!               write(*,*) 'masscont=',mascn(nel-1), mascn0(nel-1) 
               mascn(nel) = mascn0(nel)
!               stop
           endif
!..
800   CONTINUE
!
!      DO 550 NEL=1,NELX1*NELY1
!         DO 550 JJ=1,NROWB
!            STRSS(JJ,NEL)=0.0D0
! 550  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE



     END MODULE MGEOMECANICA
