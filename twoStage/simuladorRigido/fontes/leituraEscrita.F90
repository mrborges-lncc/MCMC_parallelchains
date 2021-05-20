!=================================================================================
!         programa de elementos finitos  
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
! 
!         + implementacoes de Abimael Loula
!
!         + novo projeto: modular e fortran 90, por
!         Eduardo Garcia,        bidu@lncc.br 
!         Tuane Lopes,           tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!=================================================================================
      module mLeituraEscrita
!
      implicit none
!
      integer :: iin,iecho,icoords,iconects,iconectsL
      integer :: ignuplot,iparaviewS,iparaviewP,iparaviewV
      integer :: isatTransiente
      integer :: nprint
      integer :: qtdImpSat

      integer :: ipres  = 20
      integer :: isat   = 21
      integer :: iphi   = 22
      integer :: ivel   = 23
      integer :: imass  = 24
      integer :: ifiles = 25
      integer :: imasl  = 27
      integer :: IMASG  = 28
      integer :: ISIGT  = 29
      integer :: ifdata = 150
      integer :: ifrand = 151
      integer :: ipress = 155
      integer :: iperm  = 156
      integer :: iporo  = 157
      integer :: iveloc = 158
      integer :: isaturacao = 159
!
!... NEW FOR HIERARCHICAL SETUP FILES COUNTER GT 500
!
      INTEGER :: IFEDX = 501
      INTEGER :: IFNOD = 502
      INTEGER :: IFMSH = 503
      INTEGER :: IRBAL = 504
!
      integer :: iflag_pres,iflag_sat,iflag_phi,iflag_vel
      integer :: iflag_mass, iflag_masl, IFLAG_MASG, IFLAG_SIGT
      integer :: iflag_tipoPrint
!
      integer :: nprtsat = 0
      integer :: nprtphi = 0
      integer :: nprtpres = 0
!
      character(len=128) :: ifpres_out
      character(len=128) :: ifsat_out
      character(len=128) :: ifphi_out
      character(len=128) :: ifperm_out
      character(len=128) :: ifvel_out
      character(len=128) :: ifmass_out
      character(len=128) :: ifmasl_out
      character(len=128) :: IFMASG_OUT
      character(len=128) :: IFSIGT_OUT


      real(8) :: tprt_pres, dtprt_pres
      real(8) :: tprt_sat , dtprt_sat
      real(8) :: tprt_phi , dtprt_phi
      real(8) :: tprt_vel , dtprt_vel
      real(8) :: tprt_mass, dtprt_mass
      real(8) :: tprt_masl, dtprt_masl
      real(8) :: TPRT_MASG, DTPRT_MASG
      real(8) :: TPRT_SIGT, DTPRT_SIGT
      integer :: nppres, npsat, npphi
      integer :: npvel, npmass, npmasl, NPMASG, NPSIGT
!
      logical :: apenasReservatorio=.true.
!
      contains
!
!**** new **********************************************************************
!
      subroutine lerDataIn
!
      use mPropGeoFisica
      use mGlobaisEscalares
!
      implicit none
      character(len=128) :: flag
!
      open(unit=ifdata, file= 'data.in'  )
!
      flag="# linear(1) nao-linear(2)"
      call ireadstat(flag,ifdata,iflag_linear)
!
      flag="# viscosidade da agua"
      call readstat(flag,ifdata,xmiw)
!
      flag="# viscosidade do oleo"
      call readstat(flag,ifdata,xmio)
!
      flag="# massa especifica da agua"
      call readstat(flag,ifdata,rhow)
!
      flag="# massa especifica do oleo"
      call readstat(flag,ifdata,rhoo)

!..BEGIN 4 COMPRESSIBILITY
      IF(GEOMECH==1) THEN
         FLAG="# VOLUME FORMATION"
         CALL READSTAT(FLAG,IFDATA,FORMVOL)
!
         FLAG="# WATER BULK MODULUS"
         CALL READSTAT(FLAG,IFDATA,BULKWATER)
!
         FLAG="# OIL BULK MODULUS"
         CALL READSTAT(FLAG,IFDATA,BULKOIL)
!..
         FLAG="# SOLID BULK MODULUS"
         CALL READSTAT(FLAG,IFDATA,BULKSOLID) 

      ENDIF
!..END 4 COMPRESSIBILITY
!
      flag="# saturacao residual da agua"
      call readstat(flag,ifdata,srw)
!
      flag="# saturacao residual do oleo"
      call readstat(flag,ifdata,sro)
!
      flag="# saturacao inicial da agua"
      call readstat(flag,ifdata,sinicial)

!
      flag="# saturacao na injecao"
      call readstat(flag,ifdata,sinj)
!
      flag="# saturacao do bloco"
      call readstat(flag,ifdata,sbloco)
!
      flag="# x central do bloco de sat"
      call readstat(flag,ifdata,xcbloco)
!
      flag="# y central do bloco de sat"
      call readstat(flag,ifdata,ycbloco)
!
      flag="# z central do bloco de sat"
      call readstat(flag,ifdata,zcbloco)
!
      flag="# lx do bloco de sat"
      call readstat(flag,ifdata,xlbloco)
!
      flag="# ly do bloco de sat"
      call readstat(flag,ifdata,ylbloco)
!
      flag="# lz do bloco de sat"
      call readstat(flag,ifdata,zlbloco)
!
      if(geomech==1) then
         flag="# PERMEABILITY RESERVOIR BOTTOM-REGION"
         call readstat(flag,ifdata,perminicial)
      else
         flag="# permeabilidade inicial"
         call readstat(flag,ifdata,perminicial)
      endif
!
      flag="# permeabilidade do bloco"
      call readstat(flag,ifdata,permbloco)
!
      flag="# x central do bloco de perm"
      call readstat(flag,ifdata,xcbloco_perm)
!
      flag="# y central do bloco de perm"
      call readstat(flag,ifdata,ycbloco_perm)
!
      flag="# z central do bloco de perm"
      call readstat(flag,ifdata,zcbloco_perm)
!
      flag="# lx do bloco de perm"
      call readstat(flag,ifdata,xlbloco_perm)
!
      flag="# ly do bloco de perm"
      call readstat(flag,ifdata,ylbloco_perm)
!
      flag="# lz do bloco de perm"
      call readstat(flag,ifdata,zlbloco_perm)
!
      flag="# porosidade inicial"
      call readstat(flag,ifdata,phiinicial)
!
      flag="# porosidade do bloco"
      call readstat(flag,ifdata,phibloco)
!
      flag="# x central do bloco de phi"
      call readstat(flag,ifdata,xcbloco_phi)
!
      flag="# y central do bloco de phi"
      call readstat(flag,ifdata,ycbloco_phi)
!
      flag="# z central do bloco de phi"
      call readstat(flag,ifdata,zcbloco_phi)
!
      flag="# lx do bloco de phi"
      call readstat(flag,ifdata,xlbloco_phi)
!
      flag="# ly do bloco de phi"
      call readstat(flag,ifdata,ylbloco_phi)
!
      flag="# lz do bloco de phi"
      call readstat(flag,ifdata,zlbloco_phi)
!
      if(geomech==1) then
!
         flag="# YOUNG UNDER-RIFT-REGION"
         call readstat(flag,ifdata,YNGRIFT)
!
         flag="# YOUNG CAP-TOP-REGION"
         call readstat(flag,ifdata,YNGTOP)
!
         flag="# YOUNG CAP-MIDDLE-REGION"
         call readstat(flag,ifdata,YNGMIDLE)
!
         flag="# YOUNG RESERVOIR BOTTOM REGION"
         call readstat(flag,ifdata,YNGBOTTM)
!
         flag="# YOUNG MODULUS BLOCO"
         call readstat(flag,ifdata,YNGBLOCO)
!
         flag="# X CENTRAL DO BLOCO YOUNG MODULUS"
         call readstat(flag,ifdata,xcbloco_YNG)
!
         flag="# Y CENTRAL DO BLOCO YOUNG MODULUS"
         call readstat(flag,ifdata,ycbloco_YNG)
!
         flag="# Z CENTRAL DO BLOCO YOUNG MODULUS"
         call readstat(flag,ifdata,zcbloco_YNG)
!
         flag="# LX DO BLOCO DE YOUNG MODULUS"
         call readstat(flag,ifdata,xlbloco_YNG)
!
         flag="# LY DO BLOCO DE YOUNG MODULUS"
         call readstat(flag,ifdata,ylbloco_YNG)
!
         flag="# LZ DO BLOCO DE YOUNG MODULUS"
         call readstat(flag,ifdata,zlbloco_YNG)
!
         flag="# POISSON UNDER-RIFT-REGION"
         call readstat(flag,ifdata,POISRIFT)
!
         flag="# POISSON CAP-TOP-REGION"
         call readstat(flag,ifdata,POISTOP)
!
         flag="# POISSON CAP-MIDDLE-REGION"
         call readstat(flag,ifdata,POISMIDL)
!
         flag="# POISSON RESERVOIR BOTTOM REGION"
         call readstat(flag,ifdata,POISBTTM)
!
         flag="# INITIAL PRESSURE"
         call readstat(flag,ifdata,PRSRINIT)
!
         flag="# CREEP DEFORMATION REFERENCE" 
         call readstat(flag,ifdata,CREEPZERO)
!
         flag="# TIME INCREMENT FOR CREEP" 
         call readstat(flag,ifdata,DTCREEP)
!
         flag="# STRESS REFERENCE FOR CREEP" 
         call readstat(flag,ifdata,SIGMAREF)
!
         flag="# POWER OF CREEP LAW" 
         call readstat(flag,ifdata,POWERN)
!
         flag="# TOLERANCE FOR ALL ITERATIVE" 
         call readstat(flag,ifdata,TOLCREEP)
!
      endif
!
      flag="# instante inicial"
      call readstat(flag,ifdata,tzero)
      tTransporte = tzero
!
      flag="# tempo total de simulacao"
      call readstat(flag,ifdata,tt)
!
      flag="# numero de calculos da velocidade"
      call ireadstat(flag,ifdata,nvel)
!
      if(geomech==1) then
         flag="# CREEP TIME AS MULTIPLE OF GEOTIME"
         call IREADSTAT(flag,ifdata,NCREEP)
!
         FLAG="# PASSO PARA IMPRESSAO OPENDX-FILES"
         CALL IREADSTAT(flag,ifdata,NUMDX)
!
         FLAG="# ITERATIONS FOR TRANSPORT<-->GEOMEC"
         CALL IREADSTAT(FLAG,IFDATA,NITGEO)
      endif
!  
      flag="# impressao para matlab(0) ou paraview(1) ou paraview(2)"
      call ireadstat(flag,ifdata,iflag_tipoPrint) 
!
      flag="# impressao da saturacao"
      call ireadstat(flag,ifdata,iflag_sat) 
      read(ifdata,"(a)") ifsat_out
      read(ifdata,    *) npsat
! 
!       if(iflag_sat==1)then
!          if(iflag_tipoPrint==0) then
!             ifsat_out=trim(ifsat_out)//'sat_amostra.res'
!          else
!             ifsat_out=trim(ifsat_out)//'resultadoSat.vtk'
!          end if
!          open(unit=isat,file=ifsat_out,status='unknown')
!       end if
!
      flag="# impressao da pressao"
      call ireadstat(flag,ifdata,iflag_pres)
      read(ifdata,"(a)") ifpres_out
      read(ifdata,    *) nppres
! 
!       if(iflag_pres==1)then
!          if(iflag_tipoPrint==0) then
!             ifpres_out=trim(ifpres_out)//'pres_amostra.res'
!          else
!             ifpres_out=trim(ifpres_out)//'resultadoPressao.vtk'
!          end if
!          open(unit=ipres,file=ifpres_out,status='unknown')
!       end if

!
      flag="# impressao da velocidade"
      call ireadstat(flag,ifdata,iflag_vel)
      read(ifdata,"(a)") ifvel_out
      read(ifdata,    *) npvel

!       if(iflag_vel==1)then
!          if(iflag_tipoPrint==0) then
!             ifvel_out=trim(ifvel_out)//'vel_amostra.res'
!          else
!             ifvel_out=trim(ifvel_out)//'resultadoVel.vtk'
!          end if
!          open(unit=ivel,file=ifvel_out,status='unknown')
!       end if
!
      flag="# impressao da porosidade"
      call ireadstat(flag,ifdata,iflag_phi)
      read(ifdata,"(a)") ifphi_out
      read(ifdata,    *) npphi
      ifperm_out = ifphi_out

!       if(iflag_phi==1)then
!          if(iflag_tipoPrint==0) then
!             ifphi_out=trim(ifphi_out)//'phi_amostra.res'
!             ifperm_out=trim(ifphi_out)//'perm_amostra.res'
!          else
!             ifphi_out= trim(ifphi_out)//'resultadoPhi.vtk'
!             ifperm_out=trim(ifphi_out)//'resultadoPerm.vtk'
!          end if
!          open(unit=iphi, file=ifphi_out, status='unknown')
!          open(unit=iperm,file=ifperm_out,status='unknown')
!       end if

!
      flag="# calculo da massa de agua"
      call ireadstat(flag,ifdata,iflag_mass)
      read(ifdata,"(a)") ifmass_out
      read(ifdata,    *) npmass

!       if(iflag_mass==1)then
!          ifmass_out=trim(ifmass_out)//'mass_amostra.res'
!          open(unit=imass,file=ifmass_out,status='unknown')
!       end if
!
      flag="# balanco local da massa de agua"
      call ireadstat(flag,ifdata,iflag_masl)
      read(ifdata,"(a)") ifmasl_out
      read(ifdata,    *) npmasl

!       if(iflag_masl==1)then
!          ifmasl_out=trim(ifmasl_out)//'bal_mass_amostra.res'
!          open(unit=imasl,file=ifmasl_out,status='unknown')
!       end if
!
      if(geomech==1) then
!
!... NEW FOR GLOBAL MASS BALANCE 
!
         FLAG="# balanco global de massa por elemento"
         CALL ireadstat(FLAG,IFDATA,IFLAG_MASG)
         READ(IFDATA,"(a)") IFMASG_OUT
         READ(IFDATA,    *) NPMASG

!          IF(IFLAG_MASG==1) THEN
!             IFMASG_OUT=trim(IFMASG_OUT)//'bal_mass_global.res'
!             OPEN(unit=IMASG,FILE=IFMASG_OUT,STATUS='unknown')
!          ENDIF
! 
!... NEW FOR DERIVATIVE OF TOTAL STRESS TRACE 
!
         FLAG="# derivada do tracco da tensao total"
         CALL ireadstat(FLAG,IFDATA,IFLAG_SIGT)
         READ(IFDATA,"(a)") IFSIGT_OUT
         READ(IFDATA,    *) NPSIGT 
         ENDIF

!          IF(IFLAG_SIGT==1) THEN
!             IFSIGT_OUT=trim(IFSIGT_OUT)//'tracco.res'
!          ENDIF

      !!! provisorio
      rmi=xmiw/xmio
      nsw=1
      nso=1

      if(geomech==1) then
         nsw=2
         nso=2
      endif

      end subroutine lerDataIn
!
!=================================================================================
!
      subroutine lerRandfilesIn
!
      use mPropGeoFisica
      use mMalha, only: nsd
      use mGlobaisEscalares, only: geomech, iflag_beta
!
      implicit none

      character(len=128) :: flag
!
      open(unit=ifrand, file= 'randfiles.in'  )
!
!.... Numero de realizacoes
!
      flag="# realizacoes"
      call ireadstat(flag,ifrand,nrand)
!
!.... Kozeny-Carman 
      flag="# Kozeny-Carman"
      call ireadstat(flag,ifrand,iflag_KC)
      read(ifrand,*)s1KC,s2KC
!
!.... Leitura da permeabilidade
!
      flag="# permeabilidade"
      call ireadstat(flag,ifrand,iflag_read_perm) 
      read(ifrand,"(a)") perm_inx
      read(ifrand,"(a)") perm_iny
      if(nsd==3) read(ifrand,"(a)") perm_inz
      perm_inx=trim(perm_inx)    
      perm_iny=trim(perm_iny) 
      if(nsd==3)  perm_inz=trim(perm_inz)
      read(ifrand,*) kg, rho
!
!.... Leitura da porosidade
!
      flag="# porosidade"
      call ireadstat(flag,ifrand,iflag_read_phi) 
      read(ifrand,"(a)") phi_in
      phi_in=trim(phi_in)
      read(ifrand,*) kgphi, rhophi

!
!.... Leitura da porosidade
!
      flag="# beta"
      call ireadstat(flag,ifrand,iflag_beta) 
      read(ifrand,"(a)") beta_in
      beta_in=trim(beta_in)
      read(ifrand,*) kgbeta, rhobeta
!
!.... Leitura de Modulo de Young
!
      if(geomech==1) then
         flag="# Modulo de Young"
         CALL IREADSTAT(FLAG,IFRAND,IFLAG_READ_YNG) 
         READ(IFRAND,"(A)") YNG_IN
         YNG_IN=trim(YNG_IN)
         READ(IFRAND,*) KGYNG, RHOYNG
      endif
!
!.... Saida do mixing
!
      flag="# mixing"
      call ireadstat(flag,ifrand,iflag_mix)
      read(ifrand,"(a)") mixing_out
      mixing_out=trim(mixing_out)
      read(ifrand,*)npmix
!
!.... Saida da producao
!
      flag="# producao"
      call ireadstat(flag,ifrand,iflag_prod)
      read(ifrand,"(a)") prod_out
      prod_out=trim(prod_out)
      read(ifrand,*)npprod
!
      end subroutine

!**** new **********************************************************************
      subroutine readstat(flag,ifile,x)
!
      implicit none
      integer :: ifile
      real(8) :: x
      character(len=128) :: flag,flag1
!
      read(ifile,"(a)") flag1
      if(trim(flag).eq.trim(flag1)) then
      read(ifile,*) x
      else
      write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
      end subroutine
!
!=================================================================================
!
      subroutine ireadstat(flag,ifile,n)
!
       implicit none

      integer, intent(in)  :: ifile
      integer, intent(out) :: n
      character(len=128), intent(in) :: flag
      character(len=128) :: flag1
!
      read(ifile,"(a)") flag1
      if(trim(flag).eq.trim(flag1)) then
      read(ifile,*) n
      else
      write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
      end subroutine
!
!=================================================================================
!
    subroutine abrirArquivosInformacoesMalha()
!  
!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
      character(len=2) :: nRank
      character(len=20) :: nomeIecho
!
      iin           = 15
      iecho         = 16
      icoords       = 18
      iconects      = 19
      iconectsL     = 25
      isatTransiente= 43
!
      open(unit=iin   , file= 'input.in')
!
      nomeIecho='echo.dat'
      open(unit=iecho , file= nomeIecho)
!
#ifdef debug
      open(unit=icoords    , file= 'coordenadas.dat')
      open(unit=iconects   , file= 'conectsNodais.dat')
      open(unit=iconectsL  , file= 'conectsLadais.dat')
#endif

   end subroutine 

!
!=================================================================================
!
      subroutine abrirArquivosResultados
      use mPropGeoFisica
      use mGlobaisEscalares
!       use mLeituraEscrita
!
      implicit none

      if(iflag_sat==1)then
         if(iflag_tipoPrint==0) then
            ifsat_out=trim(ifsat_out)//'sat_amostra.res'
         end if
         if(iflag_tipoPrint==1) then
            ifsat_out=trim(ifsat_out)//'resultadoSat.vtk'
         end if
         if(iflag_tipoPrint==2) then
            ifsat_out=trim(ifsat_out)//'Sat_'
         end if
         open(unit=isat,file=ifsat_out,status='unknown')
      end if

      if(iflag_pres==1)then
         if(iflag_tipoPrint==0) then
            ifpres_out=trim(ifpres_out)//'pres_amostra.res'
         end if
         if(iflag_tipoPrint==1) then
            ifpres_out=trim(ifpres_out)//'resultadoPressao.vtk'
         end if
         if(iflag_tipoPrint==2) then
            ifpres_out=trim(ifpres_out)//'Pres_'
         end if
         open(unit=ipres,file=ifpres_out,status='unknown')
      end if

      if(iflag_vel==1)then
         if(iflag_tipoPrint==0) then
            ifvel_out=trim(ifvel_out)//'vel_amostra.res'
         else
            ifvel_out=trim(ifvel_out)//'resultadoVel.vtk'
         end if
         open(unit=ivel,file=ifvel_out,status='unknown')
      end if

      if(iflag_phi==1)then
         if(iflag_tipoPrint==0) then
            ifphi_out=trim(ifphi_out)//'phi_amostra.res'
            ifperm_out=trim(ifphi_out)//'perm_amostra.res'
         end if
         if(iflag_tipoPrint==1) then
            ifphi_out= trim(ifphi_out)//'resultadoPhi.vtk'
            ifperm_out=trim(ifphi_out)//'resultadoPerm.vtk'
         end if
         if(iflag_tipoPrint==2) then
            ifphi_out =trim(ifphi_out)//'Phi_'
            ifperm_out=trim(ifperm_out)//'Perm_'
         end if
         open(unit=iphi, file=ifphi_out, status='unknown')
         open(unit=iperm,file=ifperm_out,status='unknown')
      end if

     if(iflag_mass==1)then
         ifmass_out=trim(ifmass_out)//'mass_amostra.res'
         open(unit=imass,file=ifmass_out,status='unknown')
      end if

      if(iflag_masl==1)then
         ifmasl_out=trim(ifmasl_out)//'bal_mass_amostra.res'
         open(unit=imasl,file=ifmasl_out,status='unknown')
      end if

         IF(IFLAG_MASG==1) THEN
            IFMASG_OUT=trim(IFMASG_OUT)//'bal_mass_global.res'
            OPEN(unit=IMASG,FILE=IFMASG_OUT,STATUS='unknown')
         ENDIF

         IF(IFLAG_SIGT==1) THEN
            IFSIGT_OUT=trim(IFSIGT_OUT)//'tracco.res'
         ENDIF

      end subroutine

!**** new **********************************************************************
    subroutine fecharArquivos()
!
      close(iin   )
      close(iecho )
      close(icoords )
      close(ignuplot)
      close(iconectsL)
      close(iconects )
      close(isat)
      close(ipres)
      close(ivel)
      close(iphi)
      close(iperm)
      close(imass)
      close(imasl)
      close(IMASG)
 
    end subroutine fecharArquivos
!      
! =======================================================================
!   
       subroutine inittime
! 
       use mPropGeoFisica
       use mGlobaisEscalares, only: tt, tzero

       implicit none
! 
! .... ajustes para impressao estocastico
! 
!        t0 = tzero
! 
       np_rand_mix  = 1
       np_rand_prod = 1
       np_rand_conc = 1
       np_rand_prodF= 1
       ninit_prodF  = 1
! 
! .... tamanhos dos intervalos de impressoes
! 
       dtprt_pres = (tt-tzero)/nppres
       dtprt_sat  = (tt-tzero)/npsat
       dtprt_phi  = (tt-tzero)/npphi
       dtprt_vel  = (tt-tzero)/npvel
       dtprt_mass = (tt-tzero)/npmass
       dtprt_masl = (tt-tzero)/npmasl
       dtprt_prod = (tt-tzero)/npprod
       dtprt_mix  = (tt-tzero)/npmix
       dtprt_prodF= (tt-tzero)/npprodF
       dtprt_conc = (tt-tzero)/npconc
! 
! .... inicializa os contadores de impressao
! 
       tprt_pres= dtprt_pres
       tprt_sat = dtprt_sat
       tprt_phi = dtprt_phi
       tprt_vel = dtprt_vel
       tprt_mass= dtprt_mass
       tprt_masl= dtprt_masl
       tprt_prod= dtprt_prod
       tprt_mix = dtprt_mix
       tprt_prodF= dtprt_prodF
       tprt_conc= dtprt_conc
! 
       end subroutine
!
!*****************************************************
!
     subroutine imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, perm, satElem)
       use mGlobaisEscalares, only: ndofP, ndofV, tTransporte, numdx
       use mMalha,            only: x, conecNodaisElem, conecLadaisElem, nen, nsd
       use mMalha,            only: numel, numelReserv, numLadosReserv, numLadosElem, numnp


      implicit none

      real*8, intent(in) :: pressaoElem(ndofP, numelReserv), velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: phi(numelReserv), perm(numelReserv), satElem(numelReserv)
      integer :: i, k
!
!
!.... imprime a condicao inicial: pressao
      !
       if(iflag_pres==1) then
          if(iflag_tipoPrint==0) then
             call prt(nsd,numelReserv,tTransporte,pressaoElem,ipres)
          end if
          if(iflag_tipoPrint==1) then
             if(apenasReservatorio.eqv..true.) then
                call escreverArqParaviewReservatorio(ipres, pressaoElem, ndofP, &
                     numelReserv, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
             else
                call escreverArqParaview(ipres, pressaoElem, ndofP, numel, nen, &
                     conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
             endif
          end if
          if(iflag_tipoPrint==2) then
             call PRINT_VTK(ipres,pressaoElem,tTransporte,ifpres_out,nen,ndofV,nprtpres)
          end if
       endif
!
!  
!.... imprime a condicao inicial: velocidade
!
      if(iflag_vel==1) then
         if(iflag_tipoPrint==0) then
            call prtvB(nsd,numelReserv,tTransporte,velocLadal,ndofV, conecLadaisElem,&
                                           numLadosElem,ivel) 
         else
            write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview"
         endif
      end if
!     
!.... imprime a condicao inicial: saturacao
!
      if(iflag_sat==1) then    
         if(iflag_tipoPrint==0) then 
            call prt (nsd,numelReserv,tTransporte,satElem,isat)
         end if
         if(iflag_tipoPrint==1) then 
            if(apenasReservatorio.eqv..true.) then
               call escreverArqParaviewReservatorio(isat, satElem, ndofV,&
                     numelReserv, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
            else
               call escreverArqParaview(isat, satElem, ndofV, numel, nen, &
                                  conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
            endif
         endif
         if(iflag_tipoPrint==2) then
            call PRINT_VTK(isat,satElem,tTransporte,ifsat_out,nen,ndofV,nprtsat)
         end if
      endif
!
!.... imprime a condicao inicial: porosidade
!
      if(iflag_phi==1) then
         if(iflag_tipoPrint==0) then
            call prt(nsd,numelReserv,tTransporte,phi,iphi)
         end if
         if(iflag_tipoPrint==1) then
            if(apenasReservatorio.eqv..true.) then
               call escreverArqParaviewReservatorio(iphi, phi, ndofV, numelReserv, nen, & 
                    conecNodaisElem, 1, 'porosidade', len('porosidade')) 
            else
               call escreverArqParaview(iphi, phi, ndofV, numel, nen, & 
                    conecNodaisElem, 1, 'porosidade', len('porosidade')) 
            endif
         endif
         if(iflag_tipoPrint==2) then
            call PRINT_VTK(iphi,phi,tTransporte,ifphi_out,nen,ndofV,nprtphi)
         end if
       endif
!
!.... imprime a condicao inicial: permeabilidade
!
       if(iflag_phi==1) then
          if(iflag_tipoPrint==0) then
             call prt(nsd,numelReserv,tTransporte,perm,iperm)
          end if
          if(iflag_tipoPrint==1) then
             if(apenasReservatorio.eqv..true.) then
                call escreverArqParaviewReservatorio(iperm, perm, ndofV, numelReserv, nen, &
                     conecNodaisElem, 1, 'permeabilidade', len('permeabilidade'))
             else
                call escreverArqParaview(iperm, perm, ndofV, numel, nen, &
                     conecNodaisElem, 1, 'permeabilidade', len('permeabilidade'))
             endif
          endif
          if(iflag_tipoPrint==2) then
             call PRINT_VTK(iperm,perm,tTransporte,ifperm_out,nen,ndofV,nprtphi)
          end if
       endif
!
      IF (NUMDX.EQ.0) RETURN
! ! 
! ! .... PRINT DISPLACEMENTS DATA FOR OPEN-DX FILE
! ! 
      DO 10 I=1,NUMNP
 	 WRITE(IFNOD,1900) (X(K,I),K=1,NSD)  
 10   CONTINUE 
! ! 
! ! .... PRINT CONECTIVITIES DATA FOR OPEN-DX FILE
! ! 
      DO 20 I=1,NUMEL 
         WRITE(IFMSH,2000) conecNodaisElem(1,I)-1, conecNodaisElem(2,I)-1, &
     &                     conecNodaisElem(4,I)-1, conecNodaisElem(3,I)-1 
 20   CONTINUE 
! !                     
! ! .... OPEN-DX INFORMATION DATA FILES
! ! 
! ! .... NODAL POINTS OF FINITE ELEMENT DISCRETIZATION
       IF (NUMDX.GT.0) CALL PRINT_DXINFO('OPEN__FEDX_FILE',NUMNP,NUMEL)
! !
      RETURN
!
1900  FORMAT(6(1PE11.4,2X)) 
2000  FORMAT(27I6) 
! 
       end subroutine
!
!*****************************************************
!
     subroutine imprimirSolucaoNoTempo(pressaoElem, velocLadal, tempo)
       use mGlobaisEscalares, only : ndofP, ndofV
       use mMalha,            only : nsd, numel, numelReserv, numLadosReserv, numLadosElem, conecLadaisElem

      implicit none

      real*8, intent(in) :: pressaoElem(ndofP,numelReserv), velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: tempo
!
      character(21) :: label
      character(21) :: num
   
! 
!.... imprime a velocidade no centro
!
      if(iflag_vel==1) then
         if(iflag_tipoPrint==0) then
            call prtvB(nsd,numelReserv,tempo,velocLadal,ndofV, conecLadaisElem, numLadosElem,ivel) 
         endif
      end if
!
!.... imprime a pressao no tempo 
!
      if(iflag_pres==1) then
         if(iflag_tipoPrint==0) then
            call prt(nsd,numelReserv,tempo,pressaoElem,ipres)
         else
            write(num,'(f12.5)') tempo
            label="t="//ADJUSTL(num)
            call escreverArqParaviewIntermed(ipres, pressaoElem, ndofP, numelReserv, trim(label), len(trim(label)))
         endif
      endif
!
      end subroutine 
!
!**** new **********************************************************************
!
      subroutine echo

      implicit none
!
!.... program to echo input data
!
      character*4 ia(20)
      integer :: iech, i
!      common /iounit/ iin,iout,iecho,iouter,ignuplot
!
!     cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return
!
      write(iecho,2000) iech
      backspace iin
!
      do 100 i=1,100000
      read(iin,3000,end=200) ia
      if (mod(i,50).eq.1) write(iecho,4000)
      write(iecho,5000) ia

  100 continue
!
  200 continue
      rewind iin
      read(iin,1000) iech
!
      return
!
 500  format('programa de elementos finitos em fortran 90 baseado em:',// &
      'The Finite Element Method, Hughes, T. J. R., (2003)'//)
 1000 format(16i10)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
     ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
     '    eq. 0, no echo of input data                        ',   /5x,&
     '    eq. 1, echo input data                              ',   ///)
 3000 format(20a4)
 4000 format(' ',8('123456789*'),//)
 5000 format(' ',20a4)

      end subroutine echo
!**** new **********************************************************************
      subroutine leituraGeracaoCoordenadas(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: genfl
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8, intent(inout) ::  x(nsd,*)
      integer, intent(in)   :: nsd, numnp, iin, icoords, iprtin
!
      integer :: i, n
!      
      call genfl(x,nsd,iin)
!
      if (iprtin.eq.1) return
!
#ifdef debug
         write(icoords,*) "# Coordenadas ", nsd
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd) 
         end do
#endif
!
      return
!
!  1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,&
!      ' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i12,10x,3(1pe15.8,2x))
      end subroutine
!
!**** NEW ** MODIFIED 4 HIEARARCHICAL ********************************************************
!
      SUBRoutine LEITuraGERAcaoCOORdenadas1(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: genfl
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8, intent(inout) ::  x(nsd,*)
      integer, intent(in)   :: nsd, numnp, iin, icoords, iprtin
!
      REAL(8), DIMENSION(NSD,NUMNP)  :: XLOC
!
      INTEGER  INXCORD
      CHARACTER*30 NAMEIN
!
      INTEGER :: I, N, NODE, HIERARCH
!.... .......   ............ HIERARCH = 0 NOT HIERARCH MESH, HIERARCH = 1 READ HIERARCH MESH 
!      integer :: i, n
!      
      call genfl(x,nsd,iin)
!
      xloc = 0.0d0
!      READ(IIN,1000) HIERARCH
!      write(*,1000) hierarch

!      IF (HIERARCH.EQ.1) THEN
         INXCORD  = 531
         NAMEIN   = 'nodes_rsrv_outsburden.in'
         OPEN(UNIT=INXCORD,FILE=NAMEIN)
         DO 100 NODE=1,NUMNP
            READ(INXCORD,1500) X(1,NODE), X(2,NODE)
!            WRITE(*,1500) X(1,NODE)-XLOC(1,NODE), X(2,NODE)-XLOC(2,NODE)
100      CONTINUE
         CLOSE(INXCORD)
!      ENDIF
!
!...TRANSFER ARRAYS
!

!
      if (iprtin.eq.1) return
!
#ifdef debug
         write(icoords,*) "# Coordenadas ", nsd
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd) 
         end do
#endif
!
      return
!
 1000 FORMAT(I10)
 1500 FORMAT(2X,40(1PE15.8,2X))
 2000 format(6x,i12,10x,3(1pe15.8,2x))
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine leituraCodigosCondContorno(id, ndof, numnp, neq, iin, iecho,iprtin)
!
!.... program to read, generate and write boundary condition data
!        and establish equation numbers
!
      use mMalha, only: igen

      integer, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer:: neq
      integer, intent(inout) :: id(ndof,numnp)
!
      integer :: nn, n, i
      logical pflag
!
      id = 0
      call igen(id,ndof, iin)
!
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
!
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
!
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
!
!.... establish equation numbers
!
      neq = 0
!
      do 400 n=1,numnp
!
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
!
  300 continue
!
  400 continue
!
      return
!
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
      end subroutine
!
!**** new **********************************************************************
!
      subroutine leituraValoresCondContorno(f,ndof,numnp,j,nlvect,iprtin)
!
!.... program to read, generate and write nodal input data
!
!        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
!                             = nodal body forces(j=1)
!
      use mMalha, only : genfl
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: ndof, numnp, j, nlvect, iprtin
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer nlv
      character(len=35) :: rotulo
!
!     call clear(f,nlvect*numnp*ndof)
      f(1:nlvect,1:numnp,1:ndof)=0.0

      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof,iin)
!       
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
! 
      if (iprtin.eq.0) then
!
         if (lzero) then
            if (j.eq.0) write(iecho,1000) nlv
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
!
            if (j.eq.1) then
               rotulo=" n o d a l  b o d y  f o r c e s"
               call printd (rotulo, f,ndof,numnp,iecho)
            end if
!
         endif
      endif
!
  100 continue
!
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
      end subroutine

!**** new **********************************************************************
      subroutine printf(f,ndof,numnp,nlv)
!
!.... program to print prescribed force and boundary condition data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer ndof, numnp, nlv
      real*8 :: f(ndof,numnp,*)
!
      logical lzero
      integer :: nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',&
     ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
     '  b o u n d a r y   c o n d i t i o n s'//5x,&
     ' load vector number = ',i10///5x,&
     ' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end subroutine
      
!**** new **********************************************************************
      subroutine printd(name,dva,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof,numnp, icode
      character (LEN=*) ::  name
      real*8 dva(ndof,*)
!
      logical lzero
      integer nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
!       call ztest(dva(1,n),ndof,lzero)
!       if (.not.lzero) then
!          nn = nn + 1
!          if (mod(nn,50).eq.1) &
!            write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
!       endif
  100 continue
!
      return
!
 1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe30.10,2x))
      end subroutine
!
!**** new **********************************************************************
!
      subroutine printp(a,idiag,neq,nsq,*)
      use mGlobaisEscalares
!
!.... program to print array d after Crout factorization 
!        a = u(transpose) * d * u
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: neq, nsq
      real*8 :: a(*)
      integer :: idiag(*)
!
      integer :: n, i
!
      do 100 n=1,neq
      if (mod(n,50).eq.1) write(iecho,1000) nsq
      write(iecho,1000)
      i = idiag(n)
      write(iecho,2000) n,a(i)
  100 continue
!
      return 1
!
 1000 format('1',' array d of factorization',/&
     ' a = u(transpose) * d * u ',                                //5x,&
     ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end subroutine
!
!**** new **********************************************************************
!
      subroutine prntel(mat,conectElem,nen,numel,tipo)
      implicit none
!
!.... program to print data for element with "nen" nodes
!
!        note: presently the label formats are limited to
!              elements with one to nine nodes
!
      integer :: nen, numel
      integer :: mat(*),conectElem(nen,*)
      integer :: tipo
!
      integer n, i
!
      if(tipo==1) then
      write(iconects,*) "# Conectividades nodais"
      do n=1,numel
        write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if

      if(tipo==2) then
      write(iconectsL,*) "# Conectividades ladais"
      do  n=1,numel
        write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if
!
      return
!
 2000 format(1x,i10,9(2x,i10))
 3000 format(1x,i10,7(2x,i10))
      end subroutine
!
!**** new **********************************************************************
!
      subroutine printResultado(dva, ndof, numnp, inicio, fim, icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof, numnp, inicio, fim, icode
      real*8  :: dva(ndof,numnp)
!
      integer :: n, i
!
      write(icode,*) "# Solucao"
      do 100 n=inicio,fim
         write(icode,2000) n,(dva(i,n),i=1,ndof)
         !write(*,*) n,(dva(i,n),i=1,ndof)
  100 continue
!
      return
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine
!
!**** new **********************************************************************
!
      subroutine prtgnup(name,x,dva,nsd,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: nsd, ndof, numnp, icode
      character*4 name(11)
      real*8 :: x(nsd,*),dva(ndof,*)
!
      integer :: n, j, i
!
      write(icode,*) name
      do 100 n=1,numnp
         write(icode,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
  100 continue
!
      return
!
 2000 format(6(1pe13.6,2x))
      end subroutine
!
!**** new **********************************************************************
!
    subroutine escreverArqParaviewReservatorio(arquivo, campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
    use mMalha, only: x, nsd, numelReserv, numnpReserv

    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)
    integer :: nen
    integer :: conectElem(nen,numelReserv)
    integer :: tipo  !tipo=1 para elemento, e tipo=2 para no
    integer :: tamRot

    character(len=tamRot) :: rotulo

  
    write(arquivo,'(a)')'# vtk DataFile Version 3.0'
    write(arquivo,'(a)')'vtk output'
    write(arquivo,'(a)')'ASCII'
    write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '

    call escreverPontosNodais  (arquivo,x, numnpReserv, nsd)
! 
    write(arquivo,'(a,i10,i10)')'CELLS', numelReserv , (nen+1) * numelReserv
    call escreverConectividades(arquivo,conectElem, numelReserv, nen, nsd) !apenas reservatório
! 
    write(arquivo,'(a,i10)')'CELL_TYPES ', numelReserv
    call escreverTiposElementos(arquivo,numelReserv,nsd)
! 

    if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numelReserv

    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  dim1*dim2

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

    call escreverEscalaresPorElemento(arquivo,campo, dim2,tamRot)

    end subroutine escreverArqParaviewReservatorio
!
!**** new **********************************************************************
!
    subroutine escreverArqParaview(arquivo, campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
    use mMalha, only: x, nsd, numel,numelReserv, numnp, numnpReserv
!
    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)
    integer :: nen
    integer :: conectElem(nen,numel)
    integer :: tipo  !tipo=1 para elemento, e tipo=2 para no
    integer :: tamRot

    character(len=tamRot) :: rotulo

  
    write(arquivo,'(a)')'# vtk DataFile Version 3.0'
    write(arquivo,'(a)')'vtk output'
    write(arquivo,'(a)')'ASCII'
    write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '

    call escreverPontosNodais  (arquivo,x, numnp, nsd)
! 
    write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
    call escreverConectividades(arquivo,conectElem, numel, nen, nsd) !todo o domínio
! 
    write(arquivo,'(a,i10)')'CELL_TYPES ', numel
    call escreverTiposElementos(arquivo,numel,nsd)
! 

    if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numel

    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  dim1*dim2

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

    call escreverEscalaresPorElemento(arquivo,campo, dim2,tamRot)
!
    end subroutine escreverArqParaview

!**** new **********************************************************************
      subroutine escreverPontosNodais  (arquivo,coords, numnp, nsd)
      implicit none
      integer, intent(in) :: arquivo,numnp, nsd
      real*8,  intent(in) :: coords(nsd,numnp)
!
      real*8  :: coordZ = 0.0 
      integer :: d, i
!
      if(nsd==2) then 
        do i=1,numnp
            write(arquivo,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd), coordZ 
        end do
      end if

      if(nsd==3) then
        do i=1,numnp
            write(arquivo,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd)
        end do
      end if
      end subroutine escreverPontosNodais


!**** new **********************************************************************
      subroutine escreverConectividades(arquivo,conectElem, numel, nen, nsd)
      implicit none
      integer, intent(in)  :: arquivo,numel, nen, nsd
      integer, intent(in)  :: conectElem(nen,numel)
!
      integer n, i
!
      if(nsd==2) then
      do  n=1,numel
        write(arquivo,'(i10,9(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

      if(nsd==3) then
      do  n=1,numel
        write(arquivo,'(i10,18(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

 end subroutine escreverConectividades

!**** new **********************************************************************
      subroutine escreverTiposElementos(arquivo,numel, nsd)
      implicit none
      integer, intent(in)   :: arquivo, numel, nsd
!
      integer :: i
!
      if(nsd==2) then
      do  i =1,numel
        write(arquivo,'(a)') '9'!trim(adjustl(tipo))
      end do
      end if 
!
      if(nsd==3) then
      do  i =1,numel
        write(arquivo,'(a)') '12'!trim(adjustl(tipo))
      end do
      end if 

      end subroutine escreverTiposElementos
!
!**** new **********************************************************************
!
      subroutine escreverEscalaresNodais(arquivo,v, tam1, tam2, rotulo, tamRot)
!
      use mMalha, only: numelReserv
!
      implicit none
      integer, intent(in)  :: arquivo,tam1,tam2
      real*8, intent(in)   :: v(tam1,tam2)
      integer :: tamRot
      character(len=tamRot) :: rotulo
!
      character(len=tamRot+5) ::  rotuloN
      integer :: i,j
      character(len=5):: eixo
      real*8 :: limite,zero
!
      limite=1.e-15
      zero=0.0d0
      do i=1,tam1

        if(i>1) then
           write(eixo,'(i0)') i
           rotuloN=trim(rotulo)//trim(eixo)
           write(arquivo,'(3a)')'SCALARS ', trim(rotuloN), ' float '
           write(arquivo,'(a)')'LOOKUP_TABLE default'
        endif
!
      if(apenasReservatorio.eqv..true.) then
        do j=1, tam2
             if(v(i,j).lt.limite) then
                write(arquivo,*) zero
             else 
                write(arquivo,*) v(i,j)
             end if
        end do
!
      else

       do j=1, tam2
          if(tam2<=numelReserv) then
             if(v(i,j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(i,j)
             end if
          else
             write(arquivo,*) zero
          endif
      end do
!
        endif
!
      end do
!
      end subroutine escreverEscalaresNodais
!
!**** new **********************************************************************
!
      subroutine escreverEscalaresPorElemento(arquivo,v, tam, tamRot)
!
      use mMalha, only: numelReserv
!
      implicit none
      integer, intent(in)  :: arquivo,tam
      real*8, intent(in)   :: v(tam)
      integer :: tamRot
!
      integer :: j
      real*8 :: limite,zero

      limite=1.e-15
      zero=0.0d0

      if(apenasReservatorio.eqv..true.) then
      do j=1, tam
              if(v(j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(j)
             end if
      end do
!
      else
!
      do j=1, tam
          if(tam<=numelReserv) then
             if(v(j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(j)
             end if
          else
             write(arquivo,*) zero
          endif
      end do
      endif
!
      end subroutine 
!
!**** new **********************************************************************
!
    subroutine escreverArqParaviewIntermed(arquivo, campo, dim1, dim2, rotulo, tamRot)
    use mMalha, only: x, nsd, numel, numnp

    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)

    integer :: tamRot
    character(len=tamRot) :: rotulo

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

     call escreverEscalaresNodais(arquivo,campo, dim1, dim2, rotulo, tamRot)

     end subroutine escreverArqParaviewIntermed
!
!----------------------------------------------------------------------
!
      subroutine paraview_geraCase(steps)

      implicit none

      integer :: steps
!
      integer :: numInicial, incremento
      real*8  :: incTempo
      integer :: i

      numInicial=0
      incremento=1
      incTempo =0.0

      open(unit=124,file="./out/transiente.case",status="unknown")

      write(124, "('FORMAT',/,'type:',2x,'ensight')")
      write(124, *)
      write(124, "('GEOMETRY',/,'model:',2x,'solucao.geo')")
      write(124, *)
      write(124, "('VARIABLE',/,'scalar per element:', 2x, 'Saturacao', 2x, 'solucao.***' )")
      write(124, *)
      write(124, "('TIME',/,'time set: 1')")
      write(124, "('number of steps:', i10)") steps
      write(124, "('filename start number:', i10)") numInicial
      write(124, "('filename increment:', i10)") incremento
      write(124, "('time values:')")

      do i=1, steps+1
           write(124, *) incTempo
           incTempo=incTempo+1.0
      end do

      end subroutine
!
!----------------------------------------------------------------------
!
    subroutine paraview_geometria(numel,numnp,nsd,x,conecNodaisElem)

    implicit none

      integer :: numel,numnp,nsd
      real(8), dimension(nsd,*) :: x   
      integer :: conecNodaisElem(8,numel)
!
      integer :: i
      open(unit=125,file="./out/solucao.geo",status="unknown")

      write(125,'(a)')'Title1'
      write(125,'(a)')'Title2'
      write(125,'(a)')'node id given'
      write(125,'(a)')'element id given'
      write(125,'(a)')'coordinates'
      write(125,'(i8)')  numnp

      do i = 1, numnp
      WRITE (125,'(I8,3E12.5)') I,x(1,i),x(2,i),x(3,i)
      enddo

      WRITE (125,'(A,/,A,/,A,/,I8)')                     &
                                'part 1'           ,    &
                                'malha'            ,    &
                                'hexa8'            ,    &
                                 numel

      WRITE (125,'(9I8)')  (I,conecNodaisElem(1,i),conecNodaisElem(2,i),conecNodaisElem(3,i), &
                                        conecNodaisElem(4,i),conecNodaisElem(5,i),conecNodaisElem(6,i), & 
                                        conecNodaisElem(7,i),conecNodaisElem(8,i),i=1, numel ) 

     end subroutine paraview_geometria

!
!*****************************************************
!
      subroutine imprimirCaseParaview(x, conecNodaisElem, pressaoElem, satElem, phi, perm)

      use mMalha,               only: numnp,nsd,numel,nen,numelReserv

      implicit none

      real*8 :: x(nsd, numnp)
      integer :: conecNodaisElem(nen,*)
      real*8 ::  pressaoElem(*), satElem(*), phi(*), perm(*)
!
      if(iflag_sat==1) then    
         if(iflag_tipoPrint==1) then 
            open(unit=ipress    , file= './out/solucao.P0001')
            open(unit=iperm     , file= './out/solucao.K0001')
            open(unit=iporo     , file= './out/solucao.PHI0001')
!             open(unit=iveloc      , file= 'solucao.V0001')
            open(unit=isaturacao, file= './out/solucao.S0001')

            call paraview_geraCase(qtdImpSat)
            call paraview_geometria(numel,numnp,nsd,x, conecNodaisElem)
            call paraview_escalarPorElemento(numel, pressaoElem,ipress)
            call paraview_escalarPorElemento(numel, perm,     iperm)
            call paraview_escalarPorElemento(numel, phi,        iporo)
            call paraview_escalarPorElemento(numelReserv, satElem,    isaturacao)
         endif
      endif

      end subroutine

!
!----------------------------------------------------------------------
!
      subroutine paraview_vetorPorElemento(numel,campo,iarq)

! ainda precisa implementar
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElemento(numel,campo,iarq)
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

      print*, "gerando", iarq
!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElementoTransiente(numel,campo,passo,iarq)
      implicit none
!
      integer :: numel,passo,iarq
      real(8), dimension(*) :: campo   
      character(len=128) :: name,sol
      character(len=8)   :: c
      integer :: i
      real(4) :: x     
!      
      x=0.001

      sol="solucao"

      write(c,"(f7.3)") x*passo
      c=adjustl(c)
      name='./out/'//trim(sol)//c 
!      
      open(unit=iarq,file=name,status="unknown")
!
      write(iarq,"('Ensight Scalar passo ',i5)") passo
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!      
!     imprime as coordenadas
!
       write(iarq,"(6(e12.5))") (real(campo(i)),i=1,numel)
!      
      close(iarq)
!
      passo=passo+1
!
      end subroutine          
!
!=======================================================================
!     
      subroutine prt(nsd,numel,t0,u,iunit)
!      
      use mMalha, only: xc
!      
      implicit none
!     
!     imprime campos escalares para o gnuplot ou para o matlab
!
      integer                   :: nsd,numel
      real(8), dimension(*)     :: u
      real(8)                   :: t0
      integer :: iunit
!     
      integer :: nel
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
!     
      do nel=1,numel
         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!     
!=======================================================================
!
      subroutine prtvB(nsd,numel,t0,velocLadal,ndofV,  conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: xc
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV ,numLadosElem
      real(8), dimension(ndofV,*) :: velocLadal
      real(8)                   :: t0
      integer                   :: conecLadaisElem(numLadosElem,numel)
!
      integer :: nel
      real(8) :: vc(nsd)
      real*8 :: mediaCentro
!
      integer :: iunit

      vc=0.0
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        vc(1) = (velocLadal(1,conecLadaisElem(2,nel))+velocLadal(1,conecLadaisElem(4,nel)))/2.0
        vc(2) = (velocLadal(1,conecLadaisElem(1,nel))+velocLadal(1,conecLadaisElem(3,nel)))/2.0
        if(nsd==3) vc(3) = (velocLadal(1,conecLadaisElem(5,nel))+velocLadal(1,conecLadaisElem(6,nel)))/2.0
        mediaCentro=sum(vc)/nsd
        write(iunit,"(6(f25.15,2x))")xc(1:nsd,nel),mediaCentro
!
      end do
!
      write(iunit,*)
!
      end subroutine
!     
!=======================================================================
!     
      subroutine prtv(nsd,numel,ndofV,numLadosReserv,t0,u,iunit)
!      
      use mMalha, only: xc
!
      implicit none
!     
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV, numLadosReserv
      real(8), dimension(ndofV,numLadosReserv) :: u
      real(8)                   :: t0
!
      integer :: nel
!     
      integer :: iunit
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!     
      do nel=1,numel

         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(1,nel),u(2,nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!
!**** NEW **** FOR DATA EXPLORER OUT PUT ************************************* 
!
      SUBROUTINE PRINT_DXINFO(TASK,NUMNP,NUMEL)

      use mGlobaisEscalares, only: NUMDX, NNP, NVEL
!
!..... PROGRAM TO SET-UP AND WRITE DATA ON OPEN-DX FORMAT
!
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      CHARACTER*15 TASK
!
      INTEGER JJ, NUMNP, NUMEL, NINDX, NNSTEP, NTINDX
!
      IF (NUMDX.EQ.0) RETURN

        IF(TASK=='OPEN__FEDX_FILE') THEN
! 
!..... PRINT NODAL AND MESH INFORMATION FOR DATA EXPLORER FILE
! 
	WRITE(IFEDX,1000) '## OpenDX format File' 
	WRITE(IFEDX,1000) '## OutPut Data  at Nodal Points in the' 
	WRITE(IFEDX,1000) '## sense of Finite Element Method'
        WRITE(IFEDX,1000) '##=========================================='
        WRITE(IFEDX,1000) '##=========================================='
! 
!..... Nodes of finite element mesh 
!
	WRITE(IFEDX,1000) '# ' 
	WRITE(IFEDX,1000) '## Nodes locations'
	WRITE(IFEDX,1000) '# '  
	WRITE(IFEDX,1500)  &
     & 'object 1 class array type float rank 1 shape 2 items ', &
     & NUMNP,' data file "fnodes.stoc"' 
!
!..... Conectivity of finite element mesh
!
	WRITE(IFEDX,1000) '# ' 
	WRITE(IFEDX,1000) '## Connectivity' 
	WRITE(IFEDX,1000) '# ' 
	WRITE(IFEDX,1500) &
     & 'object 2 class array type int rank 1 shape 4 items ', &
     & NUMEL,' data file "femesh.stoc"' 
	WRITE(IFEDX,1000) 'attribute "element type" string "quads"' 
	WRITE(IFEDX,1000) 'attribute "ref" string "positions"' 
	WRITE(IFEDX,1000) '#  '
! 
       ENDIF
!
       IF(TASK=='WRITE_FEDX_DATA') THEN
!
!..... PRINT NODAL DATA FOR OPEN-DX FILE
! 
!        NINDX = NUSTEP/10
!
        NINDX = NNP/NUMDX
!
!        NINDX = IDINT(TPRT_PHI/DTPRT_PHI)
!
        NNSTEP=24*NINDX+3
! 
!...... DISPLACEMENTS: VECTOR FIELD
! 
	WRITE(IFEDX,1000)'# Vector field : DISPLACEMENTS' 
	WRITE(IFEDX,2000) NNSTEP, NUMNP, NINDX 
! 
!...... DISPLACEMENT SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   DISPLACEMENT series'
	WRITE(IFEDX,3000) NNSTEP+1, NNSTEP
! 
!.....  PRESSURE: SCALAR FIELD
! 
	WRITE(IFEDX,1000)'# Scalar field : PRESSURE' 
	WRITE(IFEDX,4000) NNSTEP+2, NUMEL, NINDX 
! 
!...... PRESSURE SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar PRESSURE series'
	WRITE(IFEDX,3000) NNSTEP+3, NNSTEP+2
! 
!.....  GEOMECHANIC POROSITY: PORE SCALAR FIELD
! 
	WRITE(IFEDX,1000)'# Scalar field : PORE' 
	WRITE(IFEDX,4010) NNSTEP+4, NUMEL, NINDX 
! 
!...... GEOMECHANIC POROSITY: PORE SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar PORE series'
	WRITE(IFEDX,3000) NNSTEP+5, NNSTEP+4
! 
!.....  CREEP SCALAR FIELD
! 
	WRITE(IFEDX,1000)'# Scalar field : ' 
	WRITE(IFEDX,4020) NNSTEP+6, NUMEL, NINDX 
! 
!...... CREEP SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar CREEP series'
	WRITE(IFEDX,3000) NNSTEP+7, NNSTEP+6
! 
!.....  SATURATION: SCALAR FIELD
! 
	WRITE(IFEDX,1000)'# Scalar field : SATURATION' 
	WRITE(IFEDX,4030) NNSTEP+8, NUMEL, NINDX 
! 
!...... SATURATION SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar SATURATION series'
	WRITE(IFEDX,3000) NNSTEP+9, NNSTEP+8
! 
!.....  VELOCITY VECTOR FIELD
! 
	WRITE(IFEDX,1000)'# Vector field: VELOCITY' 
	WRITE(IFEDX,4040) NNSTEP+10, NUMEL, NINDX 
! 
!...... VELOCITY SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Vector Field VELOCITY series'
	WRITE(IFEDX,3000) NNSTEP+11, NNSTEP+10
! 
!.....  PERMEABILITY: SCALAR FIELD
! 
	WRITE(IFEDX,1000)'# Scalar field : PERMEABILITY' 
	WRITE(IFEDX,4050) NNSTEP+12, NUMEL, NINDX 
! 
!...... PERMEABILITY SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar PERMEABILITY series'
	WRITE(IFEDX,3000) NNSTEP+13, NNSTEP+12
! 
!.....  STRESS ON X DIRECTION
! 
	WRITE(IFEDX,1000)'# Scalar field : Stress_X' 
	WRITE(IFEDX,4060) NNSTEP+14, NUMEL, NINDX 
! 
!...... STRESS_X SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar Stress_X series'
	WRITE(IFEDX,3000) NNSTEP+15, NNSTEP+14
! 
!.....  STRESS ON Y DIRECTION
! 
	WRITE(IFEDX,1000)'# Scalar field : Stress_Y' 
	WRITE(IFEDX,4070) NNSTEP+16, NUMEL, NINDX 
! 
!...... STRESS_Y SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar Stress_Y series'
	WRITE(IFEDX,3000) NNSTEP+17, NNSTEP+16
! 
!.....  SHEAR STRESS XY
! 
	WRITE(IFEDX,1000)'# Scalar field : Stress_XY'
	WRITE(IFEDX,4080) NNSTEP+18, NUMEL, NINDX 
! 
!...... SHEAR STRESS_XY SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar Stress_XY series'
	WRITE(IFEDX,3000) NNSTEP+19, NNSTEP+18
! 
!.....  STRESS Z
! 
	WRITE(IFEDX,1000)'# Scalar field : Stress_Z'
	WRITE(IFEDX,4090) NNSTEP+20, NUMEL, NINDX 
! 
!...... SHEAR STRESS_XY SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar Stress_Z series'
	WRITE(IFEDX,3000) NNSTEP+21, NNSTEP+20
! 
!.....  YOUNG MODULUS 
! 
	WRITE(IFEDX,1000)'# Scalar field : Young Modulus'
	WRITE(IFEDX,4100) NNSTEP+22, NUMEL, NINDX 
! 
!...... YOUNG MODULUS SERIES INFORMATION 
! 
	WRITE(IFEDX,1000)'# Next object is a member of the: '
	WRITE(IFEDX,1000)'#   Scalar Young modulus series'
	WRITE(IFEDX,3000) NNSTEP+23, NNSTEP+22
!
      ENDIF
!
      IF(TASK=='CLOSE_FEDX_FILE') THEN
! C
! C..... PRINT SERIES LINKS INFORMATION FOR DATA EXPLORER FILE
! C
        NTINDX=(NVEL/NUMDX) 

!        NTINDX = IDINT(TPRT_PHI/DTPRT_PHI)+1
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the DISPLACEMENT series object'
      WRITE(IFEDX,1000)'object "displacement" class series'
      DO 301 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-20,JJ-1
 301  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the PRESSURE series object'
      WRITE(IFEDX,1000)'object "pressure" class series'
      DO 302 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-18,JJ-1
 302  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the PORE series object'
      WRITE(IFEDX,1000)'object "pore" class series'
!
      DO 303 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-16,JJ-1
 303  CONTINUE
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the CREEP series object'
      WRITE(IFEDX,1000)'object "creep" class series'
      DO 304 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-14,JJ-1
 304  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the SATURATION series object'
      WRITE(IFEDX,1000)'object "saturation" class series'
      DO 305 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-12,JJ-1
 305  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the VELOCITY series object'
      WRITE(IFEDX,1000)'object "velocity" class series'
      DO 306 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-10,JJ-1
 306  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the PERMEABILITY series object'
      WRITE(IFEDX,1000)'object "permeability" class series'
      DO 307 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-8,JJ-1
 307  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the STRESS X series object'
      WRITE(IFEDX,1000)'object "stress_x" class series'
      DO 308 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-6,JJ-1
 308  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the STRESS Y series object'
      WRITE(IFEDX,1000)'object "stress_y" class series'
      DO 309 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-4,JJ-1
 309  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the SHEAR STRESS series object'
      WRITE(IFEDX,1000)'object "stress_xy" class series'
      DO 310 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ-2,JJ-1
 310  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the STRESS Z series object'
      WRITE(IFEDX,1000)'object "stress_z" class series'
      DO 311 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ,JJ-1
 311  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Here we create the YOUNG MODULUS serie object'
      WRITE(IFEDX,1000)'object "young" class series'
      DO 312 JJ=1,NTINDX
        WRITE(IFEDX,7000) JJ-1,24*JJ+2,JJ-1
 312  CONTINUE
!
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'# Structure of VARIAVEL OF DATA FILE'
      WRITE(IFEDX,1000)'object "campos" class group'
      WRITE(IFEDX,1000)'member "displacement" value "displacement"'
      WRITE(IFEDX,1000)'member "pressure" value "pressure"'
      WRITE(IFEDX,1000)'member "pore" value "pore"'
      WRITE(IFEDX,1000)'member "creep" value "creep"'
      WRITE(IFEDX,1000)'member "saturation" value "saturation"'
      WRITE(IFEDX,1000)'member "velocity" value "velocity"' 
      WRITE(IFEDX,1000)'member "permeability" value "permeability"' 
      WRITE(IFEDX,1000)'member "stress_x" value "stress_x"' 
      WRITE(IFEDX,1000)'member "stress_y" value "stress_y"' 
      WRITE(IFEDX,1000)'member "stress_xy" value "stress_xy"' 
      WRITE(IFEDX,1000)'member "stress_z" value "stress_z"' 
      WRITE(IFEDX,1000)'member "young" value "young"' 
      WRITE(IFEDX,1000)'#  ' 
      WRITE(IFEDX,1000)'end'
      WRITE(IFEDX,1000)'#  ' 
!
      ENDIF
!
      IF(TASK=='OTHERS__ANOTHER') THEN
!
      ENDIF
!
!.... FORMATOS DE SAIDA  OPEN-DX
!
 1000 FORMAT(A) 
 1500 FORMAT(A,I7,A)
 1800 FORMAT(27I6)
 2000 FORMAT('object ',I5,                                 &
     & ' class array type float rank 1 shape 2 items', I8, &
     &' data file "disp',I3.3,'.stoc"'/                    &
     &'attribute "dep" string "positions"'/'#  ') 
!
 3000 FORMAT('object ',I5,' class field'/    &
     &'component "positions" value 1'/       &
     &'component "connections" value 2'/     &
     &'component "data" value ',I5/'#  ')
!
 4000 FORMAT('object ',I5,                     &
     &' class array type float rank 0 items ', &
     & I8,' data file "prsr',I3.3,'.stoc"'/    &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4010 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "pore',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4020 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "crep',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4030 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "satr',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4040 FORMAT('object ',I5,                                 &
     & ' class array type float rank 1 shape 2 items', I8, &
     &' data file "velt',I3.3,'.stoc"'/                    &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4050 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "perm',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4060 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "sigx',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4070 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "sigy',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
4080  FORMAT('object ',I5,                         &
     & ' class array type float rank 0 items', I8, &
     &' data file "sigt',I3.3,'.stoc"'/            &
     &'attribute "dep" string "connections"'/'#  ') 
!
4090  FORMAT('object ',I5,                         &
     & ' class array type float rank 0 items', I8, &
     &' data file "sigz',I3.3,'.stoc"'/            &
     &'attribute "dep" string "connections"'/'#  ') 
!
4100  FORMAT('object ',I5,                          &
     & ' class array type float rank 0 items', I8,  &
     &' data file "yung',I3.3,'.stoc"'/             &
     &'attribute "dep" string "connections"'/'#  ') 
7000  FORMAT('member ',I5,' value ',I5,' position ',I5)
!
      END SUBROUTINE
!
!*** NEW ***************************************************************
!
      SUBROUTINE SETUPDX()
!
      use mGlobaisEscalares, only: NUMDX, NITGEO
!..
!...  PROGRAM TO SETUP MANAGER FILES FOR GRAPHICAL INTERFACE OPEN-DX
!
      CHARACTER*30 NIFEDX,NIFNOD,NIFMSH,NIRBAL
!
      CHARACTER*2 ASTEP
!
      IF (NUMDX.EQ.0) RETURN
!
      WRITE(ASTEP,'(I2.2)') NITGEO
!
      NIFEDX = 'dxcreep'//ASTEP//'/nodestoc.dx'
      NIFNOD = 'dxcreep'//ASTEP//'/fnodes.stoc'
      NIFMSH = 'dxcreep'//ASTEP//'/femesh.stoc'
      NIRBAL = 'dxcreep'//ASTEP//'/global.mass'
!
      OPEN(UNIT=IFEDX, FILE= NIFEDX)
      OPEN(UNIT=IFNOD, FILE= NIFNOD)
      OPEN(UNIT=IFMSH, FILE= NIFMSH)
      OPEN(UNIT=IRBAL, FILE= NIRBAL)
!
      END SUBROUTINE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE PRINT_VTK(arquivo,campo,passo,fname,nen,dim,nprt)
!
    USE mMalha,            only: x,nsd,numel,numnp
    USE mPropGeoFisica,    only: GEOFORM,REGION
    use mMalha,            only: conecNodaisElem, conecLadaisElem
!
    IMPLICIT NONE
!
    INTEGER :: ARQUIVO,ISTAT,d,i,n,nen,dim,j,ptreserv
    REAL(8), DIMENSION(dim,*) :: CAMPO
    CHARACTER(LEN=128)    :: fname,NAME
    CHARACTER(LEN=10)     :: TEMP
    REAL(8)               :: COORDZ = 0.0,PASSO
    integer               :: tipo=1
    CHARACTER(len=4)      :: EXT
    CHARACTER(len=5)      :: C
    INTEGER               :: VARIOSARQ,ZERO,NPRINT
    INTEGER               :: nprt
!
    VARIOSARQ = 1
    ZERO = 0
!
    IF(VARIOSARQ.EQ.1)THEN
       WRITE(C,300)nprt
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT)
    ELSE
       WRITE(C,300)ZERO
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT,ACCESS='APPEND')
    END IF
    WRITE(*,*)'ARQUIVO DE IMPRESSAO:',NAME
!
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',fname
       STOP
    END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(a)')'# vtk DataFile Version 3.0'
       write(arquivo,'(a)')'vtk output'
       write(arquivo,'(a)')'ASCII'
       write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
       write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '
!    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Escreve as coordenadas nodais
!
       if(nsd==2) then 
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd), coordZ 
          end do
       end if
!
       if(nsd==3) then
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd)
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve as conectividades
       write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
       if(nsd==2) then
          do  n=1,numel
             write(arquivo,'(i10,9(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!
       if(nsd==3) then
          do  n=1,numel
             write(arquivo,'(i10,18(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve o tipo de celula
       write(arquivo,'(a,i10)')'CELL_TYPES ', numel
!
       if(nsd==2) then
          do  i =1,numel
             write(arquivo,'(a)') '9'!trim(adjustl(tipo))
          end do
       end if
!
       if(nsd==3) then
          do  i =1,numel
             write(arquivo,'(a)') '12'!trim(adjustl(tipo))
          end do
       end if
!
       if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numel
!       if(tipo==1) write(arquivo,'(a,i10)')'POINT_DATA ', numel

!    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  numel*ndofP
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(TEMP,100)PASSO
    TEMP=ADJUSTL(TRIM(TEMP))
    WRITE(*,*)'TEMPO da IMPRESSAO:',TEMP
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(3a,i10)')'SCALARS ', 'Sw' , ' float ',dim
    ELSE
       write(arquivo,'(3a,i10)')'SCALARS ', 't='//TRIM(TEMP)//'' , ' float ',dim
    END IF
    write(arquivo,'(2a)')'LOOKUP_TABLE ','default'
!
    DO I=1,NUMEL
       WRITE(ARQUIVO,*)(CAMPO(J,I),J=1,DIM)
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CLOSE(ARQUIVO) 
    nprt = nprt + 1
100 FORMAT(F10.5)
200 FORMAT(E15.7)
300 FORMAT(I5)
!
  END SUBROUTINE PRINT_VTK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE
