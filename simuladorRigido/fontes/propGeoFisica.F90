! 
!         programa de elementos finitos em fortran 90 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!
!         Eduardo Garcia e Tuane Lopes
!         bidu@lncc.br, tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!
!=================================================================================
!
      module mPropGeoFisica
!
      integer :: iflag_linear
      real(8) :: xmio,xmiw,sro,srw
      real(8) :: rhow, rhoo
      real(8) :: rmi
      integer :: nsw, nso
      real(8) :: xcbloco,ycbloco,zcbloco,xlbloco,ylbloco,zlbloco
      real(8) :: sbloco,sinicial
      real(8) :: sinj
      real(8) :: gf1,gf2,gf3
      real(8) :: eps_df=1.0d+0
      real(8) :: zer_df=1.0d-2
!
      real*8, allocatable  :: perm(:), permkx(:), permky(:), permkz(:)
      real*8, allocatable  :: phi(:), phi0(:)
      real*8, allocatable  :: YOUNG(:), PORE(:), PORE0(:), PHIEULER(:)
      real*8, allocatable  :: MASCN(:),MASCN0(:)
      INTEGER, ALLOCATABLE :: REGION(:)
      INTEGER, ALLOCATABLE :: MECLAW(:)
!
!... NEW VAR 4 COMPRESSIBILITY: PHYSICAL MASS CONTENT
! 
      REAL(8), ALLOCATABLE :: MASSCONT(:)
      REAL(8), ALLOCATABLE :: MASSCONT0(:)
!
      CHARACTER*12, ALLOCATABLE :: GEOFORM(:)
!
      real(8) :: xcbloco_perm,ycbloco_perm,zcbloco_perm
      real(8) :: xlbloco_perm,ylbloco_perm,zlbloco_perm
      real(8) :: permbloco,perminicial
!
      real(8) :: xcbloco_YNG, ycbloco_YNG, zcbloco_YNG
      real(8) :: xlbloco_YNG, ylbloco_YNG, zlbloco_YNG
      REAL(8) :: YNGINICIAL, YNGBLOCO
!
      real(8) :: xcbloco_phi,ycbloco_phi,zcbloco_phi
      real(8) :: xlbloco_phi,ylbloco_phi,zlbloco_phi
      real(8) :: phibloco,phiinicial
!
!     Descricao: dados dos campos randomicos
! 
      integer :: nrand,nr
!
      integer :: ncont_mix,ncont_prod
      integer :: iflag_mix,iflag_prod
      integer :: iflag_read_phi,iflag_read_perm,IFLAG_READ_YNG
!
      integer :: iflag_KC
      real(8) :: s1KC,s2KC
!
      real(8) :: tprt_mix,tprt_prod,dtprt_mix,dtprt_prod
      real(8) :: tprt_prodF, dtprt_prodF, tprt_conc, dtprt_conc
      integer :: np_rand_prod,np_rand_mix,npmix,npprod,np_rand_conc
      integer :: np_rand_prodF, ninit_prodF, npprodF, npconc 
!
      character(len=128) :: perm_inx,perm_iny,perm_inz,phi_in, beta_in
      character(len=128) :: mixing_out,prod_out
      CHARACTER(len=128) :: YNG_IN
!
      real(8) :: kg,kgphi, kgbeta    ! media geometrica
      real(8) :: rho,rhophi, rhobeta ! coeficiente da variancia (strenght)
      real(8) :: KGYNG, RHOYNG 
!
!.... Geometria
!
      integer :: nelx,nely,nelz
      integer :: nelxReserv, nelyReserv, nelzReserv
      real(8) :: dimx,dimy,dimz
      real(8) :: hx,hy,hz
!
!.... Geomecanica
!

      real*8 :: POISSRSR, POISSCAP
      real*8 :: CREEPZERO, STRSZERO
      real*8 :: PRSRINIT, DTCREEP, POWEREPS, XTEST
!
      REAL(8) :: PERMBOTM, YNGBOTTM
      REAL(8) :: YNGRIFT, YNGMIDLE, YNGTOP
!
      REAL(8) :: POISBTTM, POISMIDL, POISTOP, POISRIFT
      REAL(8) :: SIGMAREF, POWERN, TOLCREEP
!..BEGIN 4 COMPRESSIBILTY
      REAL(8) :: FORMVOL, BULKWATER, BULKOIL, BULKSOLID
!
!..END 4 COMPRESSIBILTY
!
!funcoes e subrotinas
      public :: calcdim, readperm, mapeando, atribuirPermBloco, atribuirPhiBloco, readphi
      public :: calcphi, estatistica,  xkkc, xlt, xlw,  xlo, xkrw, xkro
!
      contains
!     
!======================================================================
!                                              
      subroutine calcphi(nsd,numel,mat,c,p,p0,phi,phi0)
!
      use mGlobaisEscalares, only: geomech
      use mGlobaisArranjos,  only: beta
! 
!     Objetivo: atualiza a porosidade
!
!----------------------------------------------------------------------
! 
      implicit none
      integer :: nsd,numel
      integer, dimension(*)   :: mat
      real(8), dimension(6,*) :: c
      real(8), dimension(*)   :: phi0,phi
      real(8), dimension(*), intent(in) :: p,p0
      integer :: m,nel
      real(8) :: xlambda, xmu
      REAL(8) :: ALAM, AMU2, XBULK, POISSON
!
      do nel=1,numel
!      
!...  Armazena a solucao anterior
!
      phi0(nel)=phi(nel)
!
!.... material
!
      m = mat(nel)     
!
!.... compressibilidade
!
!.... coeficientes de Lame
!
!       xlambda=c(2,m)
!       xmu=c(3,m)
!       beta=nsd/(nsd*xlambda+2.d0*xmu)
        
      if(geomech==1) then
         POISSON = POISBTTM
         AMU2  = YOUNG(NEL)/(1.0D0+POISSON)
         ALAM  = AMU2*POISSON/(1.0D0-2.0D0*POISSON)
         XBULK = ALAM + AMU2/3.0D0
         BETA = 1.0D0/XBULK
      endif
!
!.... calculo da porosidade
! 
!       phi(nel)=1.d0-(1.d0-phi0(nel))*dexp(-beta*(p(nel)-p0(nel)))

!       phi(nel)= beta*(1+phi0(nel))*(p(nel)-p0(nel))
!       phi(nel)= beta*(1-phi0(nel))*(p(nel)-p0(nel))

!       phi(nel) = 1.0 - (1.0-phi0(nel))/(1.0-beta*(p(nel)-p0(nel)))

       phi(nel)=phi0(nel) + beta(nel)*(p(nel)-p0(nel)) !Corrigido 

!        phi(nel)=1.d0-(1.d0-phi0(nel))*dexp(-beta*(p(nel)-p0(nel))) !Antigo
!
      end do ! nel
!
      return
      end subroutine
!
!=======================================================================
!
      function xkkcGeo(phi)
      implicit none
!
!.... Permeabilidade de Kozeny-Carman
!
      real(8) :: xkkcGeo,phi
      real(8) :: c0,cs
!
      c0=1.d0
      cs=48.d0
!
      xkkcGeo=1.d0!*cs*(phi**3)/((1.d0-phi)**1)
!
      end function
!     
!=======================================================================
!  
      function xkkc(phi,perm,iflag,c0,cs)
        implicit none
!     
!.... Permeabilidade de Kozeny-Carman
!
        real(8) :: xkkc,phi,perm
        real(8) :: c0,cs
        integer :: iflag
        !
        if(iflag.EQ.1)then
           xkkc=c0*(phi**3)/(((1.d0-phi)*cs)**2)
        else
           xkkc=perm
        end if
!      
      end function xkkc
!     
!=======================================================================
!     
      function xlt(uu)
!     
!     calcula a mobilidade total
!     
      implicit none
      real(8) :: xlt,uu
!     
!     caso linear
!    
!      xlt=1.d0
!
!       xlt=xlw(uu)+xlo(uu)
!
      xlt=xkrw(uu)/xmiw+xkro(uu)/xmio
!
      end function
!     
!=======================================================================
!     
      function xltGeo(uu)
!     
!     calcula a mobilidade total
!     
      implicit none
      real(8) :: xltGeo,uu
!     
!     caso linear
!    
!      xlt=1.d0
!
       xltGeo=xlw(uu)+xlo(uu)
!
!       xlt=xkrw(uu)/xmiw+xkro(uu)/xmio
!
      end function
!     
!=======================================================================
!     
      function xlw(uu)
!     
!     calcula a mobilidade da agua
!     
      implicit none
      real(8) :: xlw,uu
!     
!.... Agua-Oleo
!
!     mobilidades
!
      xlw=xkrw(uu)/xmiw
!
      end function
!     
!=======================================================================
!     
      function xlo(uu)
!     
!     calcula a mobilidade total
!     
      implicit none
      real(8) :: xlo,uu
!     
!.... Agua-Oleo
!
!     mobilidades
!
      xlo=xkro(uu)/xmio
!
      end function

!     
!=======================================================================
!     
      function xkrw(uu)
!     
!     calcula a permeabilidade relativa da agua
!
      implicit none
!
      real(8) :: xkrw,uu
      real(8) :: uus,zero
!
      zero=0.d0
      xkrw=0.d0
!
!.... Permeabilidades relativas
!
!     Agua-Oleo
!
      uus=uu-srw
      if(uus.lt.zero) uus=zero

!      write(*,*) ' iflag_linear = ',  iflag_linear
      select case(iflag_linear)
      case(1) 
      xkrw=uus/(1.d0-srw)
      case(2)
      xkrw=uus**2/(1.d0-srw)**nsw
      end select
!
      end function
!     
!=======================================================================
!     
      function xkro(uu)
!     
!     calcula a permeavilidade relativa do oleo
!
      implicit none
!
      real(8) :: xkro,uu
      real(8) :: uus,zero
!     
      zero=0.d0
      xkro=0.d0
!
!.... Permeabilidades relativas
!
!     Agua-Oleo
!
      uus=1.d0-sro-uu

!      write(*,*) ' iflag_linear = ',  iflag_linear
      if(uus.lt.zero) uus=zero
      select case(iflag_linear)
      case(1) 
      xkro=uus/(1.d0-sro)
      case(2)
      xkro=uus**2/(1.d0-sro)**nso
      end select
!
      end function
!
!======================================================================
!      
      subroutine calcdim(nsd,numnp,x)
!
      use mGlobaisEscalares, only: geomech
!
      implicit none
!
      integer :: nsd,numnp
      real(8), dimension(nsd,*) :: x   
!
      dimx=dabs(x(1,numnp)-x(1,1))
      dimy=dabs(x(2,numnp)-x(2,1))
      if(nsd==3)dimz=dabs(x(3,numnp)-x(3,1))
   
      hx  =dimx/nelxReserv
      hy  =dimy/nelyReserv
      if(nsd==3)hz  =dimz/nelzReserv

      end subroutine        
!
!**** new **********************************************************************
!
     subroutine lerPropriedadesFisicas()
      use mMalha, only: numel,numnp,nsd,numLadosElem,nen
      use mMalha, only: x, xc, conecNodaisElem
      use mMalha, only: numelReserv       
      use mGlobaisEscalares, only: geomech, ligarBlocosHetBeta, iflag_beta, novaMalha
      use mGlobaisArranjos,  only: beta
!
      implicit none
!
      real*8  :: xlx, xly, xlz
      integer :: nelemx,nelemy,nelemz,nelem,npperm,i

!     cria o vetor de permeabilidades
      npperm = numel*100
      allocate(perm(npperm));           perm=0.d0
      allocate(permkx(numelReserv));    permkx=0.d0
      allocate(permky(numelReserv));    permky=0.d0
      if(nsd==3) then
         allocate(permkz(numelReserv)); permkz=0.d0
      end if
      allocate(phi(numelReserv))
      allocate(phi0(numelReserv))
      if(geomech==1) then
         allocate(PORE(numelReserv));     PORE     = 0.0D0
         allocate(PORE0(numelReserv));    PORE0    = 0.0D0
         allocate(YOUNG(NUMEL));          YOUNG    = 0.0D0
         allocate(PHIEULER(numelReserv)); PHIEULER = 0.0D0
      endif
!
      nelem=numel
      if(geomech==1) nelem=nelxReserv*nelyReserv
    
!.... leitura de dados
!
!...  cria o vetor de porosidade
!
      if(iflag_read_phi==1) then
!
         call readphi(phi_in,phi0,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
         call mapeando(phi,phi0,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
!
      else
!     
         call atribuirPhiBloco(nsd,numelReserv,phi,xc)
!
      end if
!
      if(iflag_read_perm==1.or.iflag_KC==1) then
         if(iflag_KC==1)then
            write(*,*)iflag_KC, s1KC,s2KC
            call KCperm(perm,phi,nelemx,nelemy,nelemz,nelem,nsd,iflag_KC,s1KC,s2KC)
            call mapeando(permkx,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            call mapeando(permky,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            if(nsd==3) then
               call mapeando(permkz,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            end if
         else
            !
            call readperm(perm_inx,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
            call mapeando(permkx,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            call readperm(perm_iny,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
            call mapeando(permky,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            if(nsd==3) then
               call readperm(perm_inz,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
               call mapeando(permkz,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
            endif
         end if
!
      else
!
         IF(GEOMECH==1) THEN        
            call READPERMBLOCO0(numelReserv,permkx)
            call READPERMBLOCO0(numelReserv,permky)
            if(nsd==3) call READPERMBLOCO0(numelReserv,permkz)
         ELSE 
            call atribuirPermBloco(nsd,numelReserv,permkx,permky,permkz,xc)
         ENDIF
!
      end if
!
!.... leitura de dados
!  
!
!       print*, "em leitura", iflag_read_beta
!       stop
      if(iflag_beta==1) then  
!
      call readbeta(beta_in,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)

      call mapeando(beta,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
     
      else

      if(ligarBlocosHetBeta.eqv..true.) call atribuirBetaBloco(nsd,numel,beta,x)
!
      end if

! 
!
!..Bbar.. BEGIN
!...   
      if(geomech==1) then
!...  MOUNT YOUNG MODULUS STOCHASTIC ARRAY 
!
      IF(IFLAG_READ_YNG==1) THEN
!
        CALL READYNG(PERM,XLX,XLY,NELEMX,NELEMY)
!
!...   mapeia o campos de permeabilidades nos elementos
!
        CALL mapeando(young,perm,conecNodaisElem,xc,xlx,xly,xlz,nelemx,nelemy,nelemz, nelem,nen, nsd)
!
      ELSE

        if(novaMalha.eqv..false.) then
           CALL READYNGBLOCO(YOUNG, REGION, NUMEL)
        else
           CALL READYNGBLOCO1(YOUNG, GEOFORM, NUMEL)
        endif

!
      ENDIF 
!
!... MOVE INITIAL STOCHASTIC POROSITY (PHI) TO GEOMECHANICS (PORE) FIELD
!
      PORE=PHI

      endif
!
!..Bbar.. END

     end subroutine lerPropriedadesFisicas
!
!=======================================================================
!                                              
     subroutine readperm(fname,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
!
!     Objetivo: le as permeabilidades de um arquivo
!
       use mGlobaisEscalares, only: geomech
!----------------------------------------------------------------------
! 
       implicit none
!
       real(8), dimension(*) :: perm
       integer :: nsd
!      
       integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j,k
       real(8) :: xlx,xly,xlz,beta,TOL
       character(len=128) :: NAME, FILE_IN
       character(len=128) :: fname
       character(len=4)   :: EXT,tipo
       character(len=5)   :: C
       integer :: cont, contK, fim, numel_
!
!   
       inperm=901
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       WRITE(C,113)(nr-1)
       C=ADJUSTL(C)
!
       file_in=fname
       EXT='.dat'
       NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
       NAME=ADJUSTL(TRIM(NAME))
       WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
       open(inperm, file= NAME)
!
!     dimensoes do dominio
!
       read(inperm,*) xlx
       read(inperm,*) xly
       if(nsd==3) read(inperm,*) xlz

!
!     numero de elementos em cada direcao
!     
       read(inperm,*) nelemx
       read(inperm,*) nelemy
       if(nsd==3)read(inperm,*) nelemz
!
!     verificacao do tamanho do vetor para leitura
!
       if(nsd==2) then
          if(geomech==1) then
             if(nelemx*nelemy.gt.nelxReserv*nelyReserv*100)then
                write(*,115)(nr-1)
                stop
             end if
          else
             if(nelemx*nelemy.gt.nelx*nely*100)then
                write(*,115)(nr-1)
                stop
             end if
          endif
       else
          if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
             write(*,115)(nr-1)
             stop
          end if
       endif

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

!     inicio da leitura do campo
!   
       if(nsd==2) fim=1
       if(nsd==3) fim=nelemz
       contK=0


       do k=1,fim
!
          contK=1+(k-1)*nelemx*nelemy
!
          if(nsd==3) then
             read(inperm,*) nline
             if(nline+1.ne.k) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop
             end if
          endif
!   
          do j=1,nelemy
!
             read(inperm,*) nline
             if(nline+1.ne.j) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop
             end if
!
             cont=contK+(j-1)*nelemx
             read(inperm,*) (perm(i),i=cont,nelemx+cont-1)
             read(inperm,*) nflag
             if(nflag.ne.192837465) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nflag'
                stop
             end if
!      
          end do ! nelemy
       end do ! nelemz    
!     
       close(inperm)
!
!     calculando a permeabilidade lognormal
!
       if(nsd==2) numel_=nelemy*nelemx
       if(nsd==3) numel_=nelemz*nelemy*nelemx
!
       TOL = 1.0e-16
       if(dabs(kg).gt.TOL)then
          do j=1,numel_
             perm(j)=kg*dexp(rho*perm(j))
!             perm(j)=(perm(j))
          end do
       end if
! 
       tipo='K'
       call estatistica(perm,nelemx,nelemy,nelemz,tipo)

!
111    FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113    FORMAT(I5)
115    FORMAT(                                      &
       '######################################',/, &
       'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
       'LEITURA DO CAMPO DE PERMEABILIDADES   ',/, &
       'NUMERO DO CAMPO:',I5,/,                    &
       '######################################',/)
!
     end subroutine readperm

!======================================================================
     subroutine estatistica(prm,nmx,nmy,nmz,sc)
!
       use mMalha, only: nsd
!     calculando a media e variancia
       implicit none
       real(8), dimension(*) :: prm
       integer :: nmx,nmy,nmz,i,numel_
       real(8) :: xm,xmm,xv,XMAX,XMIN
       character*4 sc
!
       xm  = 0.d0
       xmm = 0.d0
       xv  = 0.d0
       XMAX = -1.d18
       XMIN = 1.d18
!
       if(nsd==2) numel_=nmy*nmx
       if(nsd==3) numel_=nmz*nmy*nmx
       do i=1,numel_
          xm = xm+prm(i)
          xmm= xmm+prm(i)*prm(i)
          if(prm(i).gt.XMAX) XMAX = prm(i)
          if(prm(i).lt.XMIN) XMIN = prm(i)
       end do
       xm  = xm/numel_
       xmm = xmm/numel_
       xv  = xmm - xm*xm
!
!
      write(*,123)sc,xm,xv,XMAX,XMIN
 123  FORMAT(A,                          &
     & ' ######################### ',/,  &
     & 'MEDIA     =',1PE15.8,/,          &
     & 'VARIANCIA =',1PE15.8,/,          &
     & 'MAXIMO    =',1PE15.8,/,          &
     & 'MINIMO    =',1PE15.8,/,          &
     & '##############################')
!
      end subroutine estatistica
!
!=================================================================================
!
      subroutine mapeando(permk,perm,conecNodaisElem,x,xlx,xly,xlz,nelemx,nelemy,nelemz,nelem,nen,nsd)
!
!     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO

        use mGlobaisEscalares, only: geomech
!
        implicit none
!
        integer :: nelemx,nelemy,nelemz,nelem,nen,nsd
        integer,dimension(nen,*):: conecNodaisElem
        real(8),dimension(nsd,*)  :: x
        real(8),dimension(*)    :: perm,permk
!      
        integer :: nel,i,no,j,k,ncont
        real(8) :: xi,xf,yi,yf,zi,zf,xx,yy,zz,xlx,xly,xlz
        real(8) :: cdx,cdy,cdz,aux,aux1,aux2,tolx,toly,tolz
!

        if(nsd==2) then
           nelz=1
           nelemz=1
        endif

!     definicacao do tamanlhos dos elementos
!     Campo de permeabilidades
!
        cdx=xlx/nelemx
        cdy=xly/nelemy
        if(nsd==3)cdz=xlz/nelemz
!
!     Impressao dos parametros
!
        if(nsd==2) then
           if(geomech==1) then
              write(*,800)dimx,dimy,nelxReserv,nelyReserv
           else
              write(*,800)dimx,dimy,nelx,nely
           endif
           write(*,890)xlx,xly,nelemx,nelemy
           write(*,891)hx,hy,cdx,cdy
        else 
           write(*,900)dimx,dimy,dimz,nelx,nely,nelz
           write(*,990)xlx,xly,xlz,nelemx,nelemy,nelemz
           write(*,991)hx,hy,hz,cdx,cdy,cdz
        endif
!
!     verificacao da consistencia
!
        aux = (hx-cdx)
        aux1= (hy-cdy)
        if(nsd==3) aux2= (hz-cdz)
        tolx = hx*1.0e-06
        toly = hy*1.0e-06
        if(nsd==3) tolz = hz*1.0e-06
!
        if(nsd==2) then
           if(aux.gt.tolx.or.aux1.gt.toly)then
              write(*,110)
              stop
           end if
        else
           if(aux.gt.tolx.or.aux1.gt.toly.or.aux2.gt.tolz)then
              write(*,110)
              stop
           end if
        end if
!
        aux = (dimx-xlx)
        aux1= (dimy-xly)
        if(nsd==3) aux2= (dimz-xlz)
        tolx = dimx*1.0e-06
        toly = dimy*1.0e-06
        if(nsd==3) tolz = dimz*1.0e-06
!
        if(nsd==2) then
           if(aux.gt.tolx.or.aux1.gt.toly)then
              write(*,111)
              stop
           end if
        else
           if(aux.gt.tolx.or.aux1.gt.toly.or.aux2.gt.tolz)then
              write(*,111)
              stop
           end if
        endif
!
!     loop nos elementos
!             
        if(geomech==1) nelem=nelxReserv*nelyReserv
!
        if((nelx==nelemx).and.(nely==nelemy).and.(nelz==nelemz)) then
           do nel=1,nelem
              xx=0.d0
              yy=0.d0
              if(nsd==3)zz=0.d0
              do i=1,nen
                 no = conecNodaisElem(i,nel)
                 xx=xx+x(1,no)
                 yy=yy+x(2,no)
                 if(nsd==3)zz=zz+x(3,no)
              end do
              xx=xx/nen
              yy=yy/nen
              if(nsd==3)zz=zz/nen
              permk(nel)=perm(nel)
           end do
           return
        end if
!             
        do nel=1,nelem
!       
!     calculo do centroide do elemento
!
           xx=0.d0
           yy=0.d0
           if(nsd==3)zz=0.d0
           do i=1,nen
              no = conecNodaisElem(i,nel)
              xx=xx+x(1,no)
              yy=yy+x(2,no)
              if(nsd==3)zz=zz+x(3,no)
           end do
           xx=xx/nen
           yy=yy/nen
           if(nsd==3)zz=zz/nen
!
!     mapear o centroide no campo de permeabilidade
!
           ncont=0d0
!
           if(nsd==2) then
!
              ncont=0d0
              do j=1,nelemy
                 yi = (j-1)*cdy
                 yf = yi+cdy
!
                 do i=1,nelemx
                    ncont=ncont+1
                    xi = (i-1)*cdx
                    xf = xi+cdx
                    if(xx.ge.xi.and.xx.le.xf)then
                       if(yy.ge.yi.and.yy.le.yf)then
                          permk(nel)=perm(ncont)
                          go to 100
                       end if
                    end if
                 end do !i
              end do !j
!
           else

              do k=1,nelemz
                 zi = (k-1)*cdz
                 zf = zi+cdz
                 do j=1,nelemy
                    yi = (j-1)*cdy
                    yf = yi+cdy
                    do i=1,nelemx
                       ncont=ncont+1
                       xi = (i-1)*cdx
                       xf = xi+cdx
                       if(xx.ge.xi.and.xx.le.xf)then
                          if(yy.ge.yi.and.yy.le.yf)then
                             if(zz.ge.zi.and.zz.le.zf)then
                                permk(nel)=perm(ncont)
                                go to 100
                             end if
                          end if
                       end if
                    end do !i
                 end do !j
              end do !k
           endif
!
100        continue
!
        end do !nel:nelem
!
!      if(geomech==1) then
!      DO I=NELEMY*NELEMX+1,NELX*NELY
!         PERMK(I) = PERMMIDL
!      END DO
!      endif


 800  format(/,                            &
       '##############################',/, &
       'TAMANHO DO RESERVATORIO:',/,'Lx=', &
       f15.7,2x,'Ly=',f15.7,/,             &
       'MALHA:',/,'nx=',i7,' ny=',i7,/,    &
       '##############################')
!
 890  format(/,                              &
       '##############################',/,   &
       'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
       f15.7,'Ly=',f15.7,/,               &
       'MALHA:',/,'nx=',i7,' ny=',i7,/,      &
       '##############################')
!
 891  format(/,                            &
       '##############################',/, &
       'Tamanho do elemento:',/,'hx=',     &
       f10.5,2x,'hy=',f10.5,/, &
       '##############################',/, &
       'Tamanho do bloco geologico:',/,'gx=', &
       f10.5,2x,'gy=',f10.5,/,&
       '##############################',/)

 900  format(/,                            &
       '##############################',/, &
       'TAMANHO DO RESERVATORIO:',/,'Lx=', &
       f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
       'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,    &
       '##############################')
!
 990  format(/,                              &
       '##############################',/,   &
       'TAMANHO DO CAMPO ALEATORIO:',/,'Lx=',&
       f10.5,2x,'Ly=',f10.5,2x,'Lz=',f10.5,/,&
       'MALHA:',/,'nx=',i7,' ny=',i7,' nz=',i7,/,&
       '##############################')
!
 991  format(/,                            &
       '##############################',/, &
       'Tamanho do elemento:',/,'hx=',     &
       f10.5,2x,'hy=',f10.5,2x,'hz=',f10.5,/, &
       '##############################',/, &
       'Tamanho do bloco geologico:',/,'gx=', &
       f10.5,2x,'gy=',f10.5,2x,'gz=',f10.5,/,&
       '##############################',/)

 110  format('####################################',/, &
             'INCONSISTENCIA NO TAMANHO DOS BLOCOS',/, &
             'GEOLOGICOS: MALHA COMPUTACIONAL     ',/, &
             'MAIS GROSSEIRA QUE A MALHA GEOLOGICA',/, &
             '####################################',/)

 111  format('####################################',/, &
             'INCONSISTENCIA NO TAMANHO DOS       ',/, &
             'DOMINIOS: O DOMINIO GEOLOGICO DEVE  ',/, &
             'SER MAIOR OU IGUAL AO DOMINIO       ',/, &
             '####################################',/)
!
      return
!
      end subroutine
!
!=======================================================================
!     
      subroutine atribuirPermBloco(nsd,numel,ux,uy,uz,x)
!
      implicit none
!     
      integer :: nsd,numel
      real*8  :: x(nsd,*) 
      real(8), dimension(*)   :: ux,uy,uz
!
      integer :: i
      real(8) :: xx,yy,zz,xi,xf,yi,yf,zi,zf
!
      print*, "em readpermbloco"
!     
!.... Bloco     
!     
      xi=xcbloco_perm-xlbloco_perm/2.d0
      xf=xcbloco_perm+xlbloco_perm/2.d0
      yi=ycbloco_perm-ylbloco_perm/2.d0
      yf=ycbloco_perm+ylbloco_perm/2.d0
      if(nsd==3) zi=zcbloco_perm-zlbloco_perm/2.d0
      if(nsd==3) zf=zcbloco_perm+zlbloco_perm/2.d0
!     
      do i=1,numel
!
      xx=x(1,i)
      yy=x(2,i)
      if(nsd==3)zz=x(3,i)
!
      ux(i)=perminicial
      uy(i)=perminicial
      if(nsd==3) uz(i)=perminicial

      if(nsd==2) then
         if(xx.gt.xi.and.xx.lt.xf) then
         if(yy.gt.yi.and.yy.lt.yf) then
         ux(i)=permbloco
         uy(i)=permbloco
         end if
         end if
      else
         if(xx.gt.xi.and.xx.lt.xf) then
         if(yy.gt.yi.and.yy.lt.yf) then
         if(zz.gt.zi.and.zz.lt.zf) then
         ux(i)=permbloco
         uy(i)=permbloco
         uz(i)=permbloco
         end if
         end if
         end if
      endif
!     
      end do
!
      end subroutine atribuirPermBloco

!
!=======================================================================
!     
      SUBROUTINE READPERMBLOCO0(numelReserv, U)
!
      implicit none
!     
      INTEGER :: NEL,numelReserv
      REAL(8), DIMENSION(numelReserv) :: U

      DO NEL=1,numelReserv
         U(NEL) = PERMINICIAL 
      END DO
!
!      DO 200 NEL=numelReserv+1,NUMEL
!         U(NEL) = PERMMIDL
! 200  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!======================================================================
!                                              
      subroutine readphi_old(perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
!
!     Objetivo: le as permeabilidades de um arquivo
!
!----------------------------------------------------------------------
!
        implicit none
!
        real(8), dimension(*) :: perm
        integer :: nsd
!      
        integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j,k
        real(8) :: xlx,xly,xlz,beta,TOL
        character(len=128) :: NAME, FILE_IN
        character(len=4)   :: EXT, tipo
        character(len=5)   :: C
!
        if(nsd==2) then
           nelemz=1
           nelz=1
           xlz=1.0
        endif
!      
!     abre o arquivo de permeabilidades
!   
        inperm=903
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        WRITE(C,113)(nr-1)
        C=ADJUSTL(C)
        FILE_IN=phi_in
        EXT='.dat'
        NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
        NAME=ADJUSTL(TRIM(NAME))
        WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
        open(inperm, file= NAME)
!
!     dimensoes do dominio
!
        read(inperm,*) xlx
        read(inperm,*) xly
        if(nsd==3)read(inperm,*) xlz
!
!     numero de elementos em cada direcao
!      
        read(inperm,*) nelemx
        read(inperm,*) nelemy
        if(nsd==3)read(inperm,*) nelemz
!
!     verificacao do tamanho do vetor para leitura
!
        if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
           write(*,115)(nr-1)
           stop
        end if
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
        do k=1,nelemz
!
           if(nsd==3) then
              read(inperm,*) nline
              if(nline+1.ne.k) then
                 write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                 stop
              end if
           endif
!
           do j=1,nelemy
!
              read(inperm,*) nline
              if(nline+1.ne.j) then
                 write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                 stop
              end if
!
              read(inperm,*) (perm(i+(j-1)*nelemx+(k-1)*nelemx*nelemy),i=1,nelemx)
!      
              read(inperm,*) nflag
              if(nflag.ne.192837465) then
                 write(*,*) 'Erro na leitura do campo de permeabilidade, nflag'
                 stop
              end if
!      
           end do ! nelemy     
        end do ! nelemz
!     
        close(inperm)
!
!     calculando a porosidade lognormal
      !
        TOL=1.0e-05
        if(dabs(kgphi).gt.TOL)then
           do j=1,nelemz*nelemy*nelemx
              perm(j)=kgphi*dexp(rhophi*perm(j))
            !print*, j, perm(j)
           end do
        end if
! 
        tipo='PHI'
        call estatistica(perm,nelemx,nelemy,nelemz,tipo)
!
111     FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113     FORMAT(I5)
115     FORMAT(                                     &
             '######################################',/, & 
             'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
             'LEITURA DO CAMPO DE POROSIDADES       ',/, &
             'NUMERO DO CAMPO:',I5,/,                    &
             '######################################',/)
!
      end subroutine readphi_old
!
!=======================================================================
!     
      subroutine atribuirPhiBloco(nsd,numel,u,x)
!
      implicit none
!     
      integer :: nsd,numel
      real*8 :: x(nsd,*) 
      real(8), dimension(*)   :: u
!
      integer :: i
      real(8) :: xx,yy,zz,xi,xf,yi,yf,zi,zf
!
!.... Bloco     
!     
      xi=xcbloco_phi-xlbloco_phi/2.d0
      xf=xcbloco_phi+xlbloco_phi/2.d0
      yi=ycbloco_phi-ylbloco_phi/2.d0
      yf=ycbloco_phi+ylbloco_phi/2.d0
      if(nsd==3) zi=zcbloco_phi-zlbloco_phi/2.d0
      if(nsd==3) zf=zcbloco_phi+zlbloco_phi/2.d0
!     
      do i=1,numel
!     
      xx=x(1,i)
      yy=x(2,i)
      if(nsd==3) zz=x(3,i)
!
      u(i)=phiinicial

      if(nsd==2) then
         if(xx.gt.xi.and.xx.lt.xf) then      
         if(yy.gt.yi.and.yy.lt.yf) then
         u(i)=phibloco
         end if
         end if
      else
         if(xx.gt.xi.and.xx.lt.xf) then      
         if(yy.gt.yi.and.yy.lt.yf) then
         if(zz.gt.zi.and.zz.lt.zf) then
         u(i)=phibloco
         end if
         end if
         end if
      endif

! ! bloquinho1
!       xi=0.0d0; xf=0.5d0
!       yi=0.0d0; yf=0.5d0
!       zi=0.5d0; zf=1.0d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.1
!       end if
!       end if
!       end if
! 
! ! bloquinho2
!       xi=0.5d0; xf=1.0d0
!       yi=0.0d0; yf=0.5d0
!       zi=0.5d0; zf=1.0d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.3
!       end if
!       end if
!       end if
! 
! ! bloquinho3
!       xi=0.0d0; xf=0.5d0
!       yi=0.5d0; yf=1.0d0
!       zi=0.5d0; zf=1.0d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.3
!       end if
!       end if
!       end if
! 
! ! bloquinho4
!       xi=0.5d0; xf=1.0d0
!       yi=0.5d0; yf=1.0d0
!       zi=0.5d0; zf=1.0d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.1
!       end if
!       end if
!       end if
! 
! 
! ! bloquinho1
!       xi=0.0d0; xf=0.5d0
!       yi=0.0d0; yf=0.5d0
!       zi=1.0d0; zf=1.5d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.3
!       end if
!       end if
!       end if
! 
! ! bloquinho2
!       xi=0.5d0; xf=1.0d0
!       yi=0.0d0; yf=0.5d0
!       zi=1.0d0; zf=1.5d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.1
!       end if
!       end if
!       end if
! 
! ! bloquinho3
!       xi=0.0d0; xf=0.5d0
!       yi=0.5d0; yf=1.0d0
!       zi=1.0d0; zf=1.5d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.1
!       end if
!       end if
!       end if
! 
! ! bloquinho4
!       xi=0.5d0; xf=1.0d0
!       yi=0.5d0; yf=1.0d0
!       zi=1.0d0; zf=1.5d0
! 
!       if(xx.gt.xi.and.xx.lt.xf) then      
!       if(yy.gt.yi.and.yy.lt.yf) then
!       if(zz.gt.zi.and.zz.lt.zf) then
!        u(i)=0.3
!       end if
!       end if
!       end if

!     
      end do
!     
      end subroutine atribuirPhiBloco


!
!*** *** NEW: ** READPERM SUBROUTINE MODIFIED FOR YOUNG MODULUS ***********
!                                              
      SUBROUTINE READYNG(perm,xlx,xly,nelemx,nelemy)
!
!....  READ STOCHASTIC FIELD 
!
      implicit none
!
      real(8), dimension(*) :: perm
!      
      integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j
      real(8) :: xlx,xly,beta
      character(len=128) :: NAME, FILE_IN
      character(len=4)   :: EXT , tipo
      character(len=5)   :: C
!

      nelemz=0
!     abre o arquivo de permeabilidades
!   
      inperm=904
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(C,113)(nr-1)
      C=ADJUSTL(C)
!
      file_in=YNG_IN
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
      WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
      open(inperm, file= NAME)
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
!     verificacao do tamanho do vetor para leitura
!
      if(nelemx*nelemy.gt.nelx*nely*100)then
         write(*,115)(nr-1)
         stop
      end if
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
      write(*,*) 'Erro na leitura do campo de YOUNG, nline'
      stop
      end if
!
      read(inperm,*) (perm(i+(j-1)*nelemx),i=1,nelemx)
!     
      read(inperm,*) nflag
      if(nflag.ne.192837465) then
      write(*,*) 'Erro na leitura do campo de YOUNG, nflag'
      stop
      end if
!      
      end do ! nelemy     
!     
      close(inperm)
!     calculando a permeabilidade lognormal
!
      DO J=1,NELEMY*NELEMX
         PERM(j)=KGYNG*DEXP(RHOYNG*PERM(J))
      END DO
!
       tipo='EE'
       call estatistica(perm,nelemx,nelemy,nelemz,tipo)
!
 111  FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
 113  FORMAT(I5)
 115  FORMAT( &
     & '######################################',/, &
     & 'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
     & 'LEITURA DO CAMPO DE MODULO DE YOUNG   ',/, &
     & 'NUMERO DO CAMPO:',I5,/, &
     & '######################################',/)
!
      end SUBROUTINE
!
!=======================================================================
!     
      SUBROUTINE READYNGBLOCO(U,REGION,NUMEL)
!
      IMPLICIT NONE 
!     
      INTEGER :: NEL,NUMEL
      INTEGER, DIMENSION(NUMEL) :: REGION
      REAL(8), DIMENSION(NUMEL) :: U
!     
      DO 500 NEL=1,NUMEL
         IF (REGION(NEL).EQ.1) U(NEL)= YNGBOTTM
         IF (REGION(NEL).EQ.2) U(NEL)= YNGMIDLE
         IF (REGION(NEL).GT.2) U(NEL)= YNGTOP
 500  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!=======================================================================
!     
      SUBROUTINE READYNGBLOCO1(U,GEOFORM,NUMEL)
!
      IMPLICIT NONE 
!     
      INTEGER :: NEL,NUMEL
      REAL(8),      DIMENSION(NUMEL) :: U
      CHARACTER*12, DIMENSION(NUMEL) :: GEOFORM
!     
      DO 500 NEL=1,NUMEL
!          U(NEL) = GEODIC('YOUNGMD',GEOFORM(NEL))

            IF (GEOFORM(NEL).EQ.'RESERVATORIO') U(NEL) = YNGBOTTM
            IF (GEOFORM(NEL).EQ.'RIFT_CAPROCK') U(NEL) = YNGRIFT
            IF (GEOFORM(NEL).EQ.'LEFT__RESERV') U(NEL) = YNGBOTTM
            IF (GEOFORM(NEL).EQ.'RIGHT_RESERV') U(NEL) = YNGBOTTM
            IF (GEOFORM(NEL).EQ.'SALT_CAPROCK') U(NEL) = YNGMIDLE
            IF (GEOFORM(NEL).EQ.'POS_SAL_OVER') U(NEL) = YNGTOP
 500  CONTINUE
!
      RETURN
!
      END SUBROUTINE
! !
! !**** NEW ***** FOR SISMIC REPRESENTATION ******************************
! !
!       FUNCTION GEODIC(TASK,GEOFORM)
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
!=======================================================================
!                                              
      subroutine readbeta(fname,beta,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
!
!     Objetivo: le as permeabilidades de um arquivo
!
      use mGlobaisEscalares, only: geomech
!----------------------------------------------------------------------
! 
      implicit none
!
      real(8), dimension(*) :: beta
      integer :: nsd
!      
      integer :: nelemx,nelemy,nelemz,ntype,inbeta,nline,nflag,i,j,k
      real(8) :: xlx,xly,xlz,betaH
      character(len=128) :: NAME, FILE_IN
      character(len=128) :: fname
      character(len=4)   :: EXT,tipo
      character(len=5)   :: C
      integer :: cont, contK, fim, numel_
!
!   
      inbeta=909
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      WRITE(C,113)(nr-1)
      C=ADJUSTL(C)
!
      file_in=fname
      EXT='.dat'
      NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
      NAME=ADJUSTL(TRIM(NAME))
      WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
      open(inbeta, file= NAME)
!
!     dimensoes do dominio
!
      read(inbeta,*) xlx
      read(inbeta,*) xly
      if(nsd==3) read(inbeta,*) xlz

!
!     numero de elementos em cada direcao
!     
      read(inbeta,*) nelemx
      read(inbeta,*) nelemy
      if(nsd==3)read(inbeta,*) nelemz
!
!     verificacao do tamanho do vetor para leitura
!
      if(nsd==2) then
         if(geomech==1) then
            if(nelemx*nelemy.gt.nelxReserv*nelyReserv*100)then
               write(*,115)(nr-1)
               stop
            end if
         else
            if(nelemx*nelemy.gt.nelx*nely*100)then
               write(*,115)(nr-1)
               stop
            end if
         endif
      else
         if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
            write(*,115)(nr-1)
            stop
         end if
      endif

!
!     ntype = 1; campo exponencial
!     ntype = 2; campo fractal
!      
      read(inbeta,*) ntype
!
!     beta: coeficiente de Hurst
!      
      read(inbeta,*) betaH
!
!     leituras vazias
!      
      read(inbeta,*) 
      read(inbeta,*) 

!     inicio da leitura do campo
!   
      if(nsd==2) fim=1
      if(nsd==3) fim=nelemz
      contK=0


      do k=1,fim

      contK=1+(k-1)*nelemx*nelemy
!
      if(nsd==3) then
      read(inbeta,*) nline
      if(nline+1.ne.k) then
      write(*,*) 'Erro na leitura do campo de beta, nline'
      stop
      end if
      endif
!   
      do j=1,nelemy
!
      read(inbeta,*) nline
      if(nline+1.ne.j) then
      write(*,*) 'Erro na leitura do campo de beta, nline'
      stop
      end if
!
      cont=contK+(j-1)*nelemx
      read(inbeta,*) (beta(i),i=cont,nelemx+cont-1)
      read(inbeta,*) nflag
      if(nflag.ne.192837465) then
      write(*,*) 'Erro na leitura do campo de beta, nflag'
      stop
      end if
!      
      end do ! nelemy
      end do ! nelemz    
!     
      close(inbeta)
!
!     calculando a permeabilidade lognormal
!
      if(nsd==2) numel_=nelemy*nelemx
      if(nsd==3) numel_=nelemz*nelemy*nelemx

      do j=1,numel_
        beta(j)=kgbeta*dexp(rhobeta*beta(j))
      end do
! 
       tipo='K'
       call estatistica(beta,nelemx,nelemy,nelemz,tipo)

!
 111  FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
 113  FORMAT(I5)
 115  FORMAT(                                      &
       '######################################',/, &
       'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
       'LEITURA DO CAMPO DE BETA   ',/, &
       'NUMERO DO CAMPO:',I5,/,                    &
       '######################################',/)
!
      end subroutine readbeta
!
!=======================================================================
!     
      subroutine atribuirBetaBloco(nsd,numel,u,x)
!
      implicit none
!     
      integer :: nsd,numel
      real*8 :: x(nsd,*) 
      real(8), dimension(numel)   :: u
!
      integer :: i
      real(8) :: xx,yy,zz,xi,xf,yi,yf,zi,zf

      do i=1,numel
!
      xx=x(1,i)
      yy=x(2,i)
      if(nsd==3)zz=x(3,i)
!
!.... Bloco     
!     
       u(i)=1.0e-14

! ! bloquinho1
      xi=100.00d0-0.001; xf=250d0+0.001
      yi=5.00d0-0.001;   yf=15.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-12
!         print*, i, u(i)
      end if
      end if

! ! bloquinho2
      xi=100.0d0-0.001; xf=250d0+0.001
      yi=30.0d0-0.001;   yf=45.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-10
!        print*, i, u(i)
      end if
      end if

! ! bloquinho3
      xi=475.0d0-0.001; xf=625d0+0.001
      yi=5.0d0-0.001;   yf=15.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-10
!        print*, i, u(i)
      end if
      end if

! ! bloquinho4
      xi=475.0d0-0.001; xf=625d0+0.001
      yi=30.0d0-0.001;  yf=45.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-11
!        print*, i, u(i)
      end if
      end if! 
! 

! ! bloquinho5
      xi=850.0d0-0.001; xf=1000.0d0+0.001
      yi=5.0d0-0.001;   yf=15.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-18
!        print*, i, u(i)
      end if
      end if

! ! bloquinho6
      xi=850.0d0-0.001; xf=1000.0d0+0.001
      yi=30.0d0-0.001;  yf=45.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-09
!        print*, i, u(i)
      end if
      end if
! 

! ! bloquinho7
      xi=1225.0d0-0.001; xf=1375.0d0+0.001
      yi=5.0d0-0.001;    yf=15.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-09
!        print*, i, u(i)
      end if
      end if

! ! bloquinho8
      xi=1225.0d0-0.001; xf=1375.0d0+0.001
      yi=30.0d0-0.001;   yf=45.0d0+0.001

      if(xx.gt.xi.and.xx.lt.xf) then
      if(yy.gt.yi.and.yy.lt.yf) then
       u(i)=1.0e-10
!        print*, i, u(i)
      end if
      end if

!     
      end do
!     
      end subroutine atribuirBetaBloco


!=================================================================================
!
      subroutine KCperm(permk,perm,nelemx,nelemy,nelemz,nelem,nsd,flag,co,ss)
!
!     MAPEA O CAMPO DE PERMEABILIDADES NO DOMINIO DE SIMULACAO

        use mGlobaisEscalares, only: geomech
!        use mPropGeoFisica,    only: xkkc
        !
        implicit none
!
        character(len=4)   :: tipo
        integer :: nelemx,nelemy,nelemz,nelem,nen,nsd,i,j,k
        real(8),dimension(*)    :: perm,permk
        real(8)                 :: co,ss
!      
        integer :: nel,flag
!
      if(nsd==2) then
         nelemz=1
      endif     
!
      do nel=1,nelemz*nelemy*nelemx
         permk(nel)=xkkc(perm(nel),permk(nel),flag,co,ss)
!         write(*,*)nel,permk(nel)
     end do   
!
     tipo='K'
     call estatistica(permk,nelemx,nelemy,nelemz,tipo)
     return
!
    end subroutine KCperm
!
!=======================================================================
!                                              
     subroutine readphi(fname,perm,xlx,xly,xlz,nelemx,nelemy,nelemz,nsd)
!
!     Objetivo: le as permeabilidades de um arquivo
!
       use mGlobaisEscalares, only: geomech
!----------------------------------------------------------------------
! 
       implicit none
!
       real(8), dimension(*) :: perm
       integer :: nsd
!      
       integer :: nelemx,nelemy,nelemz,ntype,inperm,nline,nflag,i,j,k
       real(8) :: xlx,xly,xlz,beta,TOL
       character(len=128) :: NAME, FILE_IN
       character(len=128) :: fname
       character(len=4)   :: EXT,tipo
       character(len=5)   :: C
       integer :: cont, contK, fim, numel_
!
!   
       inperm=901
! NAME OF OUTPUT FILE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       WRITE(C,113)(nr-1)
       C=ADJUSTL(C)
!
       file_in=fname
       EXT='.dat'
       NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
       NAME=ADJUSTL(TRIM(NAME))
       WRITE(*,111)nr-1,NAME(1:LEN_TRIM(NAME))
       open(inperm, file= NAME)
!
!     dimensoes do dominio
!
       read(inperm,*) xlx
       read(inperm,*) xly
       if(nsd==3) read(inperm,*) xlz

!
!     numero de elementos em cada direcao
!     
       read(inperm,*) nelemx
       read(inperm,*) nelemy
       if(nsd==3)read(inperm,*) nelemz
!
!     verificacao do tamanho do vetor para leitura
!
       if(nsd==2) then
          if(geomech==1) then
             if(nelemx*nelemy.gt.nelxReserv*nelyReserv*100)then
                write(*,115)(nr-1)
                stop
             end if
          else
             if(nelemx*nelemy.gt.nelx*nely*100)then
                write(*,115)(nr-1)
                stop
             end if
          endif
       else
          if(nelemx*nelemy*nelemz.gt.nelx*nely*nelz*100)then
             write(*,115)(nr-1)
             stop
          end if
       endif

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

!     inicio da leitura do campo
!   
       if(nsd==2) fim=1
       if(nsd==3) fim=nelemz
       contK=0


       do k=1,fim
!
          contK=1+(k-1)*nelemx*nelemy
!
          if(nsd==3) then
             read(inperm,*) nline
             if(nline+1.ne.k) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop
             end if
          endif
!   
          do j=1,nelemy
!
             read(inperm,*) nline
             if(nline+1.ne.j) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nline'
                stop
             end if
!
             cont=contK+(j-1)*nelemx
             read(inperm,*) (perm(i),i=cont,nelemx+cont-1)
             read(inperm,*) nflag
             if(nflag.ne.192837465) then
                write(*,*) 'Erro na leitura do campo de permeabilidade, nflag'
                stop
             end if
!      
          end do ! nelemy
       end do ! nelemz    
!     
       close(inperm)
!
!     calculando a permeabilidade lognormal
!
       if(nsd==2) numel_=nelemy*nelemx
       if(nsd==3) numel_=nelemz*nelemy*nelemx
!
       TOL = 1.0e-05
       if(dabs(kgphi).gt.TOL)then
          do j=1,numel_
             perm(j)=kgphi*dexp(rhophi*perm(j))
          end do
       end if
! 
       tipo='PHI'
       call estatistica(perm,nelemx,nelemy,nelemz,tipo)

!
111    FORMAT('NAME OF INPUT FILE (RANDOM FIELD)',I5,': ',A)
113    FORMAT(I5)
115    FORMAT(                                      &
       '######################################',/, &
       'PROBLEMA NO TAMANHO DO VETOR PARA A   ',/, &
       'LEITURA DO CAMPO DE PERMEABILIDADES   ',/, &
       'NUMERO DO CAMPO:',I5,/,                    &
       '######################################',/)
!
     end subroutine readphi


    !=======================================================================

  end module
