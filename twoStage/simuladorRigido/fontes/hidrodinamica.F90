!         programa de elementos finitos em fortran 90 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!
!         Eduardo Garcia e Tuane Lopes
!         bidu@lncc.br, tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013

!
!**** new module **********************************************************************
!
      module mHidroDinamicaRT

      implicit none

      real*8,  allocatable :: pressaoElem(:,:),  pressaoElemAnt(:)
      real*8,  allocatable :: velocLadal(:,:), fVeloc(:,:), vc(:,:), ve(:,:,:)
     
      public  :: hidroGeomecanicaRT, pinit
      private :: montarMatrizVelocidadeRT,  limparSistEqAlgVelocidade, montarSistEqAlgVelocidade
      private :: calcularPressaoRT

      contains
!
!**** new **********************************************************************
!
      subroutine hidroGeomecanicaRT(satElem, satElemAnt)

      use mAlgMatricial,     only : neqV, nAlhsV, alhsV, brhsV, idiagV, idVeloc
      use mAlgMatricial,     only : btod, solverDiretoSkyLine,numCoefPorLinhaVel
      use mGlobaisArranjos,  only : mat, c, beta
      use mGlobaisEscalares, only : ndofV, optSolver, simetriaVel, tempoSolverVel
      use mGlobaisEscalares, only : dtBlocoTransp, tTransporte, nvel, nnp,ttv, ttp, tmVel
      use mMalha,            only : conecNodaisElem, conecLadaisElem, numel, numelReserv, nsd, nen
      use mMalha,            only : x, xc, listaDosElemsPorFace, numLadosReserv, numLadosElem
      use mPropGeoFisica,    only : phi, phi0, nelx, nely, nelz, calcphi, hx, hy, permkx, xkkc, xltGeo
      use mSolucoesExternas, only : solverLapack_RT, solverUMFPack, solverPardiso,solverHYPRE
      use mSolucoesExternas, only : ApVel, AiVel
      use mSolucoesExternas, only : initialGuessVel, myid, num_procs, mpi_comm
      use mLeituraEscrita,   only : iflag_vel, iflag_pres, printd
!
      implicit none
!
      real*8 :: satElem(numelReserv), satElemAnt(numelReserv)

      real*8  :: t1, t2, t3
      integer :: i, nel, arq
      character (len=80) :: nomeArqVel, nomeArqPres
      real*8 :: difP, difV, difS
      real(8) :: xlambda,xmu,bal,maxBal
      real(8), save :: tTeste
      real(8) :: elemAnt, elemDep, dx
      real(8) :: difVel, maxDifVel, minDifVel
!
#ifdef debug
      nomeArqVel ="velocidadeLadal_RT.txt"
      nomeArqPres="pressao_RT.txt"
      open(unit=206, file=nomeArqVel)
      open(unit=207, file=nomeArqPres)
#endif

! 
!******************* VELOCIDADE **********************
! 

      if(.not.allocated(alhsV)) allocate(alhsV(nalhsV))
      if(.not.allocated(brhsV)) allocate(brhsV(neqV))
!
      call limparSistEqAlgVelocidade()
! 
      call timing(t1)
      call montarSistEqAlgVelocidade(nsd, satElem, PHI, PHI0)
      call timing(t2)
!
#ifdef mostrarTempos
      write(*,*) ' tempo montagem do sistema da velocidade ', t2-t1 
#endif
      tmVel=tmVel+(t2-t1)
!
       write(*,'(a)', ADVANCE='NO') 'solucao do sistema de equacoes'
!
      if(optSolver=='lapack') then
!          write(*,'(a)') ', solver direto LAPACK'
         call solverLapack_RT(alhsV, brhsV, neqV, idiagV, numCoefPorLinhaVel, conecLadaisElem, &
                                 listaDosElemsPorFace, numLadosReserv, numLadosElem, nelx, nely,nelz)
      end if   
!
      if(optSolver=='umfpack') then
!          write(*,'(a)') ', solver direto UMFPACK, VELOCITY'
           call solverUMFPack(alhsV, brhsV, ApVel, AiVel,neqV, nalhsV)
      end if

      if(optSolver=='pardiso') then
!          write(*,'(a)') ', solver direto PARDISO, VELOCITY'
           call solverPardiso(alhsV, brhsV, ApVel, AiVel,neqV, nalhsV, simetriaVel, 'vel', 'fact')
           call solverPardiso(alhsV, brhsV, ApVel, AiVel,neqV, nalhsV, simetriaVel, 'vel', 'back')
      end if
!
      if(optSolver=='skyline') then
          write(*,'(a)') ', solver direto SKYLINE, VELOCITY'
         call solverDiretoSkyLine(alhsV, brhsV, idiagV, nalhsV, neqV, 'vel')
      end if
!
      if(optSolver=='HYPRE') then
          write(*,'(a)') ', solver iterativo HYPRE, VELOCITY'
       if(.not.allocated(initialGuessVel)) then
            write(*,'(a)') ', allocate(initialGuessVel(neqV)); initialGuessVel=0.0 '
            allocate(initialGuessVel(neqV)); initialGuessVel=0.0 
       endif
         call solverHYPRE(alhsV, brhsV, initialGuessVel, ApVel, AiVel, neqV, nalhsV, myid, num_procs, mpi_comm)
         initialGuessVel=brhsV
      endif

!
      call btod(idVeloc,velocLadal,brhsV,ndofV,numLadosReserv)

  
      call timing(t3)
#ifdef mostrarTempos
      write(*,*) ' tempo do solver velocidade', t3-t2 
      tempoSolverVel=tempoSolverVel+(t3-t2)
#endif
!
!          call printd(' velocidade  ',  &
!      &                     velocLadal,ndofV,numLados,853) 
!

      ttv=ttv+(t3-t1)
! 
      call timing(t1)

!
      print*, "calcular Pressao"
      call calcularPressaoRT (x, conecNodaisElem, conecLadaisElem,   &
                    pressaoElem, pressaoElemAnt, velocLadal,satElem,dtBlocoTransp)
!
!.... calcula a nova porosidade
!        
!      call calcphi(nsd,numelReserv,mat,c,pressaoElem,pressaoElemAnt,phi,phi0)
! 
      call calcve(numLadosElem,numelReserv,nsd,ndofV,nen, &
     &            velocLadal,ve,conecLadaisElem)
      call calcvc(nsd,nen,numelReserv,vc,ve)

#ifdef debug
         write(206,'(a,i0)') "Realizacao", nnp
         do nel=1, numelReserv
            if(nsd==2) then
               write(206,'(a,i5,a, 6(e15.5))') &
               'nel =', nel,', vel =', velocLadal(1,conecLadaisElem(4,nel)), velocLadal(1,conecLadaisElem(1,nel)), &
               velocLadal(1,conecLadaisElem(2,nel)), velocLadal(1,conecLadaisElem(3,nel))
            else
               write(206,'(a,i5,a, 6(e15.5))') &
               'nel =', nel,', vel =', velocLadal(1,conecLadaisElem(4,nel)), velocLadal(1,conecLadaisElem(1,nel)), &
               velocLadal(1,conecLadaisElem(2,nel)), velocLadal(1,conecLadaisElem(3,nel)),  &
               velocLadal(1,conecLadaisElem(5,nel)), velocLadal(1,conecLadaisElem(6,nel))
            endif
         end do
! 
         write(207,'(a,i0)') "Realizacao", nnp
         do i = 1, numelReserv
            write(207,'(i0, 5x, (e25.10))') i, pressaoElem(1,i)
         end do 
!
#endif
!
!     stop ' depois de hidrodinamica'

      call timing(t2)
      !
      ttp=ttp+(t2-t1)
!
      end subroutine hidroGeomecanicaRT
!
!**** new **********************************************************************
!
      subroutine limparSistEqAlgVelocidade()
         use mAlgMatricial,    only: alhsV, brhsV
        
         velocLadal  = 0.0d00
         if(allocated(brhsV))  brhsV = 0.0d00 
         if(allocated(alhsV)) alhsV = 0.0d00

      end subroutine
!
!**** new **********************************************************************
!
      subroutine montarSistEqAlgVelocidade(nsd, satElem, PHI, PHI0)
!
      use mGlobaisEscalares,    only: nlvectV, ndofV, dtBlocoTransp,nnp, nvel
      use mAlgMatricial,        only: load, ftod, load2
      use mAlgMatricial,        only: alhsV, brhsV, idiagV, idVeloc, lmV, neqV
      use mMalha,               only: x, conecNodaisElem, conecLadaisElem, numLadosElem
      use mMalha,               only: numLadosReserv, listaDosElemsPorFace, numelReserv
      use mSolucoesExternas,    only: solverUMFPack
!4Compressibility
      use mMalha,               only: numel
!
      implicit none
!
      integer, intent(in) :: nsd
      integer :: i
      real*8  :: satElem(numelReserv), PHI(numelReserv), PHI0(numelReserv)
!
      if (nlvectV.gt.0) CALL LOAD(idVeloc,fVeloc,     brhsV,ndofV,numLadosReserv,nlvectV)

      if (nlvectV.gt.0) call ftod(idVeloc,velocLadal,fVeloc,ndofV,numLadosReserv,nlvectV)
! 
      call montarMatrizVelocidadeRT(x, conecNodaisElem, conecLadaisElem, &
          alhsV, brhsV, idiagV, lmV, velocLadal, pressaoElemAnt, satElem,    &
          dtBlocoTransp, PHI, PHI0) 

      end subroutine montarSistEqAlgVelocidade
!
!**** new **********************************************************************
!
      subroutine montarMatrizVelocidadeRT (x, conecNodaisElem, conecLadaisElem, &
                            alhs, brhs,idiagV,lmV,          &
                            velocLadal, pressao, sw,        &
                            dtBlocoTransp, PHI, PHI0 ) 

      use mAlgMatricial,     only: addrhs, addlhs, kdbc, nALHSV, neqV,idVeloc,nedV
      use mGlobaisEscalares, only: one, ndofV, nrowsh, npint, nnp, ligarBlocosHetBeta, iflag_beta
      use mGlobaisArranjos,  only: grav, c, mat, beta
      use mfuncoesDeForma,   only: shlq, shlqrt,shgq, shgqrt
      use mfuncoesDeForma,   only: shlq3d, shg3d, shlqrt3d, shgqrt3d
      use mMalha,            only: local, nsd, numnp, numel, numelReserv, nen
      use mMalha,            only: numLadosElem, numLadosReserv
      use mPropGeoFisica,    only: xkkc,xkkcGeo,xlt,xltGeo,perm, permkx, permky, permkz, YOUNG, xlo, xlw, rhow,rhoo 
      use mPropGeoFisica,    only: gf1, gf2, gf3
      use mSolucoesExternas, only: addlhsCRS, ApVel, AiVel
!..4COMPRESSIBILITY:
      use mPropGeoFisica,    only: YOUNG, POISBTTM
      use mPropGeoFisica,    ONLY: FORMVOL, BULKWATER, BULKOIL, BULKSOLID  
      use mPropGeoFisica,    only: hx,hy,hz
      use mPropGeoFisica,    only: iflag_KC,s1KC,s2KC
!
!.... program to calculate stifness matrix and force array for a
!     singular problem element in one dimension
!     form and assemble into the global left-hand-side matrix
!                  and right-hand side vector
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!       
      real*8,  intent(in)    :: x(nsd,numnp)
      integer, intent(in)    :: conecNodaisElem(nen,numel), conecLadaisElem(numLadosELem,numelReserv)
      real(8), intent(inout) :: alhs(nalhsV), brhs(neqV)
      integer, intent(in)    :: idiagV(neqV)
      integer, intent(in)    :: lmV(ndofV,numLadosElem,numelReserv)
      real*8,  intent(inout) :: velocLadal(ndofV,numLadosReserv)
      real*8,  intent(in)    :: pressao(numelReserv), sw(numelReserv)
      real*8,  intent(in)    :: dtBlocoTransp 
!4COMPRESSIBILITY:
      real*8,  intent(in)    :: PHI(numelReserv), PHI0(numelReserv)
!
      real*8 :: xl(nsd,nen), dl(ndofV,nen), fluxl(ndofV,nen), vnl(nsd,numLadosElem)
      real*8 :: shgrt(nrowsh,numLadosElem,npint), shlrt(nrowsh,numLadosElem,npint)
      real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      real*8 :: det(npint), w(npint)
!
      real*8  :: elfrt(numLadosElem), elert(numLadosElem,numLadosElem)
      integer :: nee, neesq
!
      integer :: nel, m 
      integer :: l, i, j
      real*8  :: pi
      real*8  :: pix, piy, piz, sx, sy, sz, cx, cy, cz
!
      real(8) :: xk,yk,zk,delta,theta
      real(8) :: h,h2,xx,yy,zz,du,dux,duy,duz,xindi,xindj,ff
      real(8) :: c1,phii1,phii2,phii3,phij1,phij2,phij3,a11,a22,a33,divphij,divphii   
      real(8) :: swint,pp
      real(8) :: xlambda,xmu
      real(8) :: dt

      integer :: tid, omp_get_thread_num,omp_get_num_threads
      integer :: numPrimeiroElemento, numUltimoElemento, numThreads, inicioSol, fimSol
! 
      logical :: diag,zerol, quad,lsym
      real(8) :: xxlt,xxlo,xxlw,lambdab
!...
!... BEGIN CATALOG FOR ITERATIVE FORMULATION FOR COMPRESSIBILTY
!
!      REAL(8) :: BIOTMOD, BULK, BETASTAR
      REAL(8) :: BULKROCK, BETACOM, BIOTCOEF, BETAT, BETA0, COEFSIGM
!
      REAL(8) :: ALAM, AMU2, XBULK, POISSON,BETAGEO
!
!... END CATALOG FOR FOR ITERATIVE FORMULATION
!
!      consistent matrix
!
      nee = numLadosElem*ndofV; neesq = nee*nee
      shlrt=0.d0
      tid=1
      numThreads=1
!
      diag = .false.
      quad = .true.
      pi   = 4.0d0*datan(1.0d0)
      dt   = dtBlocoTransp
!
      if(nsd==2) then
         call shlq  (shl,w,npint,nen)
         call shlqrt(numLadosElem,npint,w,shlrt)
      else
         call shlq3d  (shl,w,npint,nen)
         call shlqrt3d(numLadosElem,npint,w,shlrt)
      endif
!      
! !$OMP PARALLEL FIRSTPRIVATE(tid) &
! !$OMP PRIVATE (numPrimeiroElemento, numUltimoElemento, inicioSol,fimSol) &
! !$OMP PRIVATE (xxlw,xxlo,xxlt,lambdab,xk,yk,h,h2,hx,hy,delta,theta, xl,fluxl,shg,shgrt,a11,a22,a33) &
! !$OMP PRIVATE (divphij,divphii ,phii1, phii2, phij1, phij2, AMU2, ALAM, XBULK,det,c1,xindi,xindj,vnl,nel,i,j,l,pp,m) &
! !$OMP REDUCTION(+:elert,elfrt) !,brhs,alhs)
! !
! #ifdef withOMP
!        tid=tid+omp_get_thread_num()
!        numThreads=omp_get_num_threads()
! #endif
!
       if(tid==1) print*, "Em Velocidade, numThreads=",numThreads
!
      numPrimeiroElemento = 1
      numUltimoElemento   = numelReserv
      call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, numThreads, tid-1, inicioSol, fimSol)
!
      do 500 nel=inicioSol,fimSol
!
!      clear stiffness matrix and force array
!
      elfrt=0.0d00; elert=0.0d00
!
!      localize coordinates and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
      call local(conecLadaisElem(1,nel),velocLadal,fluxl,numLadosElem,nedV,nedV)
!
      m = mat(nel)
      quad = .true.
!       if (nen.eq.4.and.conecLadaisElem(3,nel).eq.conecLadaisElem(4,nel)) quad = .false.
!
!     chama a shg. eh necessario pra calcular a fonte e o determinante
!
      if(nsd==2) then
         call shgq (xl,det,shl,shg,npint,nel,quad,nen)
      else
         call shg3d(xl,det,shl,shg,npint,nel,nen)
      endif
!
!... length of the element
      h2=0                                         
      h2=h2+(xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,4))**2 
      if(nsd==3)h2=h2+(xl(3,1)-xl(3,5))**2                      
      h=dsqrt(h2)/2.d00
      h2=h*h   
!     
!rt   calcula as normais externas locais
!      
      vnl=0.d0
      vnl(2,1)=-1.d0
      vnl(1,2)= 1.d0
      vnl(2,3)= 1.d0
      vnl(1,4)=-1.d0
      if(nsd==3) then
         vnl(3,5)=-1.d0
         vnl(3,6)= 1.d0
      endif
!
      if(nsd==2) then
         call shgqrt(numLadosElem,npint,hx,hy,shlrt,shgrt)
      else
         call shgqrt3d(numLadosElem,npint,hx,hy,hz,shlrt,shgrt)  
      endif
!
!..... form stiffness matrix
!
      pp=pressao(nel) 
!
!...  compressibilidade
!
!... coeficientes de Lame
!
      xlambda=c(2,m) 
      xmu=c(3,m)
!       if((iflag_beta==0).and.(ligarBlocosHetBeta.eqv..false.)) then
      beta(nel)=nsd/(nsd*xlambda+2.d0*xmu)
!       endif
!
!     parametros do penalty
      !
      delta=dt/beta(nel)
      theta=1.0d0
      
      do l=1,npint
! 
      c1 = w(l)*det(l)
!
!     calcula x, y, a pressao e suas derivadas no ponto de integracao
!
      xx=0.d0
      yy=0.d0
      if(nsd==3)zz=0.d0
      du=0.d0
      dux=0.d0
      duy=0.d0
      if(nsd==3)duz=0.d0
!
!       swint=0.d0      
!
      do i=1,nen
      xx =xx +shl(nrowsh,i,l)*xl(1,i)
      yy =yy +shl(nrowsh,i,l)*xl(2,i)
      if(nsd==3) zz =zz +shl(nrowsh,i,l)*xl(3,i)
      du =du +shl(nrowsh,i,l)*dl(1,i)
      dux=dux+shg(1,i,l)*dl(1,i)
      duy=duy+shg(2,i,l)*dl(1,i)
      if(nsd==3) duz=duz+shg(3,i,l)*dl(1,i)  
      !swint=swint+shg(3,i,l)*se(i,nel)
      end do
!
!... gravidade
!
      xxlw=xlw(sw(nel))
      xxlo=xlo(sw(nel))
      xxlt=xltGeo(sw(nel))
      lambdab=(xxlw*rhow+xxlo*rhoo)/xxlt
!
!     calculando a permeabilidade
!
!      xk=perm(nel)*xlt(sw(nel))
!      xk=perm(nel)*xlt(swint)
!
!      xk = 1.d0*c(1,m)*xkkc(phi(nel))*permk(nel)*xlt(sw(nel))
      xk = 1.d0*c(1,m)*xkkc(phi(nel),permkx(nel),iflag_KC,s1KC,s2KC)*xxlt
      a11= 1.d0/xk

      yk = 1.d0*c(1,m)*xkkc(phi(nel),permky(nel),iflag_KC,s1KC,s2KC)*xxlt
      a22= 1.d0/yk

      if(nsd==3) then
         zk = 1.d0*c(1,m)*xkkc(phi(nel),permkz(nel),iflag_KC,s1KC,s2KC)*xxlt
         a33= 1.d0/zk
      endif

!
!
!.... vetor de carga - RHS - f  =  gf0  
!
      gf1 = grav(1)
      gf2 = grav(2)
      if(nsd==3)gf3 = grav(3)
!      
      pix=pi*xx
      piy=pi*yy
      if(nsd==3)piz=pi*zz
!      
      sx=dsin(pix)
      sy=dsin(piy)
      if(nsd==3)sz=dsin(piz)
      cx=dcos(pix)
      cy=dcos(piy)
      if(nsd==3)cz=dcos(piz)
! 
!     fonte      
!
      ff=0.d0 !gf1*sx*sy+gf2*cx*cy
!      
      do j=1,numLadosElem 
!
      xindj = vnl(1,j)+vnl(2,j)
      if(nsd==3) xindj = xindj +vnl(3,j)
!
      phij1  =shgrt(1,j,l)*c1*xindj
      phij2  =shgrt(2,j,l)*c1*xindj
      if(nsd==3) phij3  =shgrt(3,j,l)*c1*xindj
      divphij=shgrt(nrowsh,j,l)*c1*xindj

!
!     termo de fonte
!
      elfrt(j)=elfrt(j) &
     &                  +  pp*divphij & ! falta a integral no bordo
     &                  + delta*ff*divphij  &
     &                  + lambdab*(gf1*phij1+gf2*phij2)
!
!     loop nos lados: indice i
!
      do i=1,numLadosElem 
!            
      xindi= vnl(1,i)+vnl(2,i)
      if(nsd==3) xindi=xindi+vnl(3,i)
!      
      phii1  =shgrt(1,i,l)*xindi
      phii2  =shgrt(2,i,l)*xindi
      if(nsd==3)phii3  =shgrt(3,i,l)*xindi 
      divphii=shgrt(nrowsh,i,l)*xindi
!
!     matriz de rigidez
!     
      elert(i,j) = elert(i,j) + a11*phii1*phij1  &
                              + a22*phii2*phij2  &
                              + delta*theta*divphii*divphij
      if(nsd==3) elert(i,j) = elert(i,j)+ a33*phii3*phij3

!          write(985,*) nel, i,j,elert(i,j)  !elert(1:6, 1:6)
!      
      end do ! numLadosElem,i
!
      end do ! numLadosElem,j
!
      end do ! nint,l

!        write(985,*)  a22,phii2,phij2!elert(1:6, 1:6)
!
!
!      computation of Dirichlet b.c. contribution
!
      call ztest(fluxl,nee,zerol)
  
      if(.not.zerol) then
         call kdbc(elert,elfrt,fluxl,nee)
      end if
!
      lsym=.true.
#ifdef withcrs
      call addlhsCRS(alhs,elert,lmV(1,1,nel),ApVel, AiVel,nee) 
#else

      call addlhs(alhs,elert,idiagV,lmV(1,1,nel),nee,diag,lsym) 
#endif
      call addrhs(brhs,elfrt,lmV(1,1,nel),nee) 
! 

  500 continue


!      stop 'depois de velocidade'

! !$OMP END PARALLEL

      return
!
      end subroutine 
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
!*** NEW ***** FUNCTION TO COMPUTE RECIPROCAL OF BIOT MODULUS ***********
!
      FUNCTION BIOTMOD(POROSITY,ALPHA,SAT)
!
!... COMPUTE RECIPROCAL OF BIOT MODULUS,THUS 1/M: REFERENCE CHENG.PDF
!     
      use mPropGeoFisica,    ONLY: BULKWATER, BULKOIL, BULKSOLID
      use mPropGeoFisica,    ONLY: RHOW, RHOO
!
      IMPLICIT NONE
!
      REAL(8) :: ALPHA, SATDENSW, SATDENSO,  DENSITEQ
      REAL(8) :: BIOTMOD, POROSITY, SAT, BULKFLUID
!
      SATDENSW = SAT*RHOW
      SATDENSO = (1.0D0-SAT)*RHOO
      DENSITEQ = SATDENSW + SATDENSO
!
      BULKFLUID = (SATDENSW/BULKWATER + SATDENSO/BULKOIL)/DENSITEQ
!
      BIOTMOD = POROSITY*BULKFLUID + (ALPHA-POROSITY)/BULKSOLID


! 
!       BULKFLUID = SAT/BULKWATER + 1.0D0/BULKOIL
! !
!       BIOTMOD = POROSITY*BULKFLUID + (ALPHA-POROSITY)/BULKSOLID
!     
      END FUNCTION
!
!*** NEW ***** FUNCTION TO COMPUTE BETA_STAR COMPRESSIBILITY ********** 
!
      FUNCTION BETASTAR(POROSITY,ALPHA,BETA,SAT)
!
!... COMPUTE BETA_STAR COMPRESSIBILITY REFERENCE TUTORIAL.PDF
!     
      IMPLICIT NONE
!
!      REAL(8) :: BIOTMOD, ALPHA, BETA, BETASTAR, POROSITY, SAT 
      REAL(8) :: ALPHA, BETA
      REAL(8) :: BETASTAR, POROSITY, SAT 
!	
      BETASTAR = BIOTMOD(POROSITY,ALPHA,SAT) + BETA*ALPHA**2
!     
      END FUNCTION
!
!**** new *************************************************************
!
      subroutine calcularPressaoRT (x, conecNodaisElem, conecLadaisElem,  &
                          pressaoElem, pressaoElemAnt, velocLadal, satElem, dtBlocoTransp)
!
!.... program to calculate stifness matrix and force array for a
!     singular problem element in one dimension
!     form and assemble into the global left-hand-side matrix
!                  and right-hand side vector
!
      use mFuncoesDeForma,   only: shlqrt, shgqrt, shlqrt3d, shgqrt3d
      use mMalha,            only: local
      use mMalha,            only: nsd, numnp, numel,numelReserv,  nen, numLadosElem, numLadosReserv
      use mGlobaisArranjos,  only: c, mat, beta
      use mGlobaisEscalares, only: ndofP, ndofV, nrowsh, npint, ligarBlocosHetBeta, iflag_beta
      use mPropGeoFisica,    only: YOUNG, POISBTTM, BULKSOLID, phi
      use mPropGeoFisica,    only: hx, hy, hz
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!       
      real*8,  intent(in)    :: x(nsd,numnp)
      integer, intent(in)    :: conecNodaisElem(nen,numel), conecLadaisElem(numLadosElem,numelReserv)
      real*8,  intent(inout) :: pressaoElem(ndofP,numelReserv),pressaoElemAnt(numelReserv)
      real*8,  intent(in)    :: velocLadal(ndofV,numLadosReserv), satElem(numelReserv)
      real*8,  intent(in)    :: dtBlocoTransp 
!
      real*8 :: xl(nsd,nen)
      real*8, dimension(nrowsh,numLadosElem,npint) :: shgrt, shlrt
      real*8  :: w(npint)
      real(8) :: dt
!
      integer :: nel,m,i,l
      real(8) :: vnl(nsd,numLadosElem)
      real(8) :: xind,div,fe
      real(8) :: pi,delta,gf1,gf2,gf3
      real(8) :: sx,sy,sz,cx,cy,cz,pix,piy,piz,xg,yg,zg,ff
      real(8) :: xlambda, xmu, dif, teste
!...
!... BEGIN CATALOG FOR ITERATIVE FORMULATION
!
      REAL(8) :: ALAM, AMU2, XBULK, POISSON, XTRACCO
!
!... END CATALOG FOR FOR ITERATIVE FORMULATION
! 
!...
!... BEGIN CATALOG FOR ITERATIVE FORMULATION FOR COMPRESSIBILTY
!
!      REAL(8) :: BIOTMOD, BULK, BETASTAR
      REAL(8) :: BULKROCK, BETACOM, BIOTCOEF, BETAT, BETA0, COEFSIGM, BETAGEO
!
!      consistent matrix
!
      pi=4.d00*datan(1.d00)
!
      dt=dtBlocoTransp
!      
!...  Armazena a solucao anterior
!
!        pressaoElemAnt=pressaoElem(1,:) !? isso é ou não necessário?
!       pressaoElem=0.0d0
!

      if(nsd==2) then
         call shlqrt(numLadosElem,npint,w,shlrt)
      else
         call shlqrt3d(numLadosElem,npint,w,shlrt)
      endif

      do nel=1,numelReserv
!
!     Numeracao local das faces
!
!                 3
!             ________
!            /  6    /|
!           /_______/ |
!           |       |2|
!         4 |   1   | /
!           |_______|/
!
!               5
!
!     vetores normais
!
      vnl=0.d0
      vnl(2,1)=-1.d0
      vnl(1,2)= 1.d0
      vnl(2,3)= 1.d0
      vnl(1,4)=-1.d0
      if(nsd==3) then
         vnl(3,5)=-1.d0
         vnl(3,6)= 1.d0
      endif
! 
!    localize coordinates and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
!      
      if(nsd==2) then
         call shgqrt(nen,npint,hx,hy,shlrt,shgrt) 
      else
         call shgqrt3d(numLadosElem,npint,hx,hy,hz,shlrt,shgrt)
      endif
!
!.... vetor de carga - RHS - f  =  gf0  
!
!       gf1 = grav(1)
!       gf2 = grav(2)
!       gf3 = grav(3)
! 
!...  fonte      
!
      ff=0.0d0 !gf1*sx*sy+gf2*cx*cy
!
!.... material
!
      m = mat(nel) 
!
!.... compressibilidade
!
!...  coeficientes de Lame
!
      xlambda=c(2,m)
      xmu=c(3,m)

      if((iflag_beta==0).and.(ligarBlocosHetBeta.eqv..false.)) then
         beta(nel)=nsd/(nsd*xlambda+2.d0*xmu)
      endif
!
      delta=dt/beta(nel)
!       write(*,'(i0,e20.5)'), nel,beta(nel)

!       
      div=0.d0
      do i=1,numLadosElem
      xind = vnl(1,i)+vnl(2,i)
      if(nsd==3) xind = xind + vnl(3,i)
      fe   = velocLadal(1,conecLadaisElem(i,nel))*xind
      div  = div + shgrt(nrowsh,i,1)*fe 
      end do ! numLadosElem
!
!.... calculo da pressao sobre o elemento
! 
       ff=ff-div
!
       pressaoElem(1,nel)=rk(pressaoElemAnt(nel),delta,ff)
!
!       dif=pressaoElem(1,nel)-pressaoElemAnt(nel)+delta*div
!       if(dif>1.e-05) stop 'balanco de massa não satisfeito'
!       print*, "nel=", nel, "Balanco=", dif

!       write(908,*) nel, pressaoElem(1,nel)
      end do ! nel


!
      return

      end  subroutine 
!
!=======================================================================
!     
      function rk(uc,dt,ff)
!     
      implicit none
!     
      real(8) :: rk,uc,dt,ff
      real(8) :: r1,r2,r3
!     
!     Evolucao no tempo (Runge-Kutta)
!     
!     primeira ordem
!     
!     rk=uc + dt*ff
!     
!     segunda ordem
!     
!       r1 = uc +dt*ff
!       r2 = uc/2.d0 + (r1+dt*ff)/2.d0
!       rk=r2
!     
!     terceira ordem
!     
      r1 = uc+dt*ff
      r2 = uc*3.d0/4.d0 + (r1+dt*ff)/4.d0
      r3 = uc/3.d0      + (r2+dt*ff)*2.d0/3.d0
      rk=r3

      end function
! 
! =======================================================================
!           
       subroutine pinit(numel,p)

       use mPropGeoFisica, only: PRSRINIT
       use mMalha, only: numelReserv 

       implicit none
!      
       integer :: numel
       integer :: i, j
       real(8), dimension(*) :: p
       real(8) :: pi
! 
       pi=4.d0*datan(1.d0)      
! 
!      condicao inicial: pressao nula
! 
!
         do i=1,numel
             p(i)=1.0e+07
          end do
!   
       end subroutine 
!     
!=======================================================================
!
      subroutine calcvc(nsd,nen,numel,v,ve)
!
!     calcula a velocidade no centro dos elementos
!
      implicit none
!      
      integer  :: nen,numel,nel,no,nsd
      real(8)  :: vx,vy
      real(8), dimension (2,*) :: v
      real(8), dimension (nsd,nen,*) :: ve
!     
      do nel=1,numel
!
         vx=0.d0
         vy=0.d0
!
         do no=1,nen
            vx=vx+ve(1,no,nel)
            vy=vy+ve(2,no,nel)
         end do
!
         v(1,nel)=vx/nen
         v(2,nel)=vy/nen
!
      end do
!
      end subroutine

!
!============================================================================
!
      subroutine calcve(nedg,numel,nsd,nedfl,nen,fluxo,vel,ieedg)
!
!     calcula a velocidade nos nos por elemento

      use mGlobaisEscalares, only: nnp
!
      implicit none
!
      integer :: numel,nsd,nen,nedfl,nedg
      real(8), dimension(nedfl,*)   :: fluxo
      real(8), dimension(nsd,nen,*) :: vel
      integer, dimension(nedg,*) :: ieedg
!
      integer :: nel,n1,n2,n3,n4,j
!
!.....  loop on elements
!
      do nel=1,numel
!      
      n1=ieedg(1,nel)
      n2=ieedg(2,nel)
      n3=ieedg(3,nel)
      n4=ieedg(4,nel)
!
      vel(1,1,nel) = fluxo(1,n4)
      vel(2,1,nel) = fluxo(1,n1)
      vel(1,2,nel) = fluxo(1,n2)
      vel(2,2,nel) = fluxo(1,n1)
      vel(1,3,nel) = fluxo(1,n2)
      vel(2,3,nel) = fluxo(1,n3)
      vel(1,4,nel) = fluxo(1,n4)
      vel(2,4,nel) = fluxo(1,n3)
!
      end do

  1001  format('nnp=',i2,1x,'elemnt 602 lado=',i1,1x,'vel_x=',1PE15.8,2X,'vel_y=',1PE15.8)

!
      return
      end subroutine

!**** new **********************************************************************
!
      subroutine velocConst()
!
      use mMalha,            only : conecNodaisElem, conecLadaisElem, numel, numelReserv, nsd, nen
      use mMalha,            only : x, xc, listaDosElemsPorFace, numLadosReserv, numLadosElem
!
      implicit none
!
      integer :: i, j, nel
! 
!******************* VELOCIDADE **********************
! 
      do nel=1, numelReserv
         velocLadal(1,conecLadaisElem(1,nel)) = 0.0
         velocLadal(1,conecLadaisElem(2,nel)) = 0.0
         velocLadal(1,conecLadaisElem(3,nel)) = 0.0
         velocLadal(1,conecLadaisElem(4,nel)) = 0.20
         velocLadal(1,conecLadaisElem(5,nel)) = 0.0
         velocLadal(1,conecLadaisElem(6,nel)) = 0.0
      end do
      !
    end subroutine velocConst

  end module mHidrodinamicaRT
