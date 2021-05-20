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
!     ************************************************************
!     *                                                          *
!     *                                                          *
!     *         A LINEAR STATIC FINITE ELEMENT PROGRAM FOR       *
!     *                                                          *
!     *                 GALERKIN METHOD                     *
!     *                                                          *
!     *                                                          *
!     ************************************************************
!
  program reservoirSimulator
!
      use mGlobaisEscalares, only: exec
      use mLeituraEscrita,   only: abrirArquivosInformacoesMalha, fecharArquivos
!
      implicit none
!
!.... initialization phase
!
      call abrirArquivosInformacoesMalha ()
      print*, ""
      print*, "PREPROCESSAMENTO" 
      call preprocessador()
!
!.... solution phase
!
      if(exec==1) then
           print*, ""
           print*, "Iniciando o PROCESSAMENTO..."
           call processamento()
      endif
!
      call fecharArquivos()
!
end program reservoirSimulator


!**** new **********************************************************************
      subroutine preprocessador()
!
       use mGlobaisArranjos,  only: title
       use mGlobaisEscalares, only: ndofP, ndofV, ndofD, nlvectP, nlvectV, nlvectD, exec, iprtin
       use mGlobaisEscalares, only: simetriaVel, optSolver
       use mAlgMatricial,     only: idVeloc, idDesloc, idiagV, neqV, nalhsV
       use mMalha,            only: nen, nsd, numel, numLadosElem, numLadosReserv, numnp, x
       use mMalha,            only: numnpReserv, numelReserv, numelReserv
       use mMalha,            only: conecLadaisElem,conecNodaisElem,listaDosElemsPorFace,listaDosElemsPorNo
       use mLeituraEscrita,   only: iin, iecho, icoords, echo, leiturageracaocoordenadas, leituracodigoscondcontorno
       use mLeituraEscrita,   only: inittime, leituraValoresCondContorno, lerDataIn, lerRandfilesIn
       use mLeituraEscrita,   only: SETUPDX, abrirArquivosResultados
       use mPropGeoFisica,    only: nelx, nely, nelz,  nr, nrand, lerPropriedadesFisicas, hx, hy, hz
       use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv, PRSRINIT
       use mHidrodinamicaRT,  only: fVeloc, pressaoElem, pinit,pressaoElemAnt
       use mTransporte,       only: satElem, satinit
       use mcmc
!
      implicit none
!
      real*8 :: t1, t2
      character(len=21) :: label
      integer :: n, na
!
!.... input phase
!
      call echo
!
      read(iin,1000) title
      if (title(1).eq.'*end') return
!
      read(iin,'(4i10)') exec,iprtin,nsd
      read(iin,'(7i10)') numnp, numel, nelx, nely, nelz
      read(iin,'(3i10)') nlvectP, nlvectV, nlvectD
!
      write(iecho,1000) title 
      write(iecho,3000) exec, iprtin, nsd
      write(iecho,4000) numnp, numnp, &
                        ndofP, ndofV, &
                        nlvectP, nlvectV
      WRITE(IECHO,5000) nelx,nely,nelz


      nelxReserv=nelx
      nelyReserv=nely
      nelzReserv=nelz
      numelReserv=numel
      if(nsd==2) numLadosReserv=(nelxReserv+1) * (nelyReserv+1) * 2 - ( (nelxReserv+1)+(nelyReserv+1) ) 
      if(nsd==3)then
         numLadosReserv=(nelxReserv*nelyReserv*nelzReserv) &
         +(nelxReserv*nelyReserv) + ((nelxReserv+1)*nelyReserv*nelzReserv) &
         +  (nelyReserv+1)*nelxReserv*nelzReserv
      end if
      print*, "numLadosReserv=", numLadosReserv, nelxReserv, nelyReserv, nelzReserv

!
!.... inicializa os parametros da simulacao
!
      call lerDataIn
      call abrirArquivosResultados
      call lerRandfilesIn
!
!.... initialization phase
!
      ndofP = 1
      ndofV = 1
      ndofD = 2
      if(nsd==2) then 
         numLadosElem = 4
         nen=4
      endif
      if(nsd==3) then
         numLadosElem = 6
         nen=8
      endif
      nr=1
      nrand=1
  
      call alocarMemoria()
!
!.... input coordinate data
!
      call leituraGeracaoCoordenadas(x,nsd,numnp, iin, icoords, iprtin)
!
!.... input boundary condition data and establish equation numbers
!
      call leituraCodigosCondContorno(  idVeloc, ndofV,numLadosReserv,n, iin, iecho, iprtin)
      neqV=n
      allocate(idiagV(neqV));  idiagV=0
!
      print*, "ndofV=", ndofV, "neqV=", neqV
!
!.... input nodal force and prescribed kinematic boundary-value data
!
       if (nlvectV.gt.0) call leituraValoresCondContorno(fVeloc, ndofV, numLadosReserv,1,nlvectV,iprtin)
!
!.... input element data
!
       call topologiaMalhaSistEquacoes(NALHSV, NEQV) 
!
!.... inicializa os tempos de impressao
!
      call inittime
!
!.... estabelece a condicao inicial para a saturacao
! 
      call satinit(nsd,numelReserv,satElem)
!
!.... estabelece a condicao inicial para a pressao
!           
      call pinit(numelReserv,pressaoElem)
      call pinit(numelReserv,pressaoElemAnt)
      call INITMCMC()
!
 1000 format(20a4)
 3000 format(a///&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec ) = ',i10//5x,&
     '    eq. 0, data check                                   ',   /5x,&
     '    eq. 1, execution                                    ',  //5x,&
     ' input data print code . . . . . . . . . . . (iprtin) = ',i10//5x,&
     '    eq. 0, print nodal and element input data           ',   /5x,&
     '    eq. 1, do not print nodal and element input data    ',   /5x, &
     ' number of space dimensions  . . . . . . . . (nsd   ) = ',i10)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . . . (numnpP) = ',i10//5x,&
     ' number of nodal points  . . . . . . . . . . (numnpD) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofP ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofV ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . . (ndofD ) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectP) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectV) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . . (nlvectD) = ',i10//5x)
 5000 FORMAT(5X,  &
     &' MESH DATA FOR RESERVOIR AND OVERBUDEN DOMAINS:        '//5X,     &
     &'    ELEMENTS IN X-DIRECTION GLOBAL DOMAIN. . (NELX ) = ',I10//5X, &
     &'    ELEMENTS IN Y-DIRECTION GLOBAL DOMAIN. . (NELY ) = ',I10//5X, &
     &'    ELEMENTS IN Z-DIRECTION GLOBAL DOMAIN. . (NELZ ) = ',I10//5X)

!
      end subroutine preprocessador
!
!**** new **********************************************************************
!
      subroutine processamento()
        use mAlgMatricial,     only : nAlhsP, neqV, nAlhsV, alhsV, brhsV
        use mGlobaisArranjos,  only : mat, c
        use mGlobaisEscalares
        use mLeituraEscrita,   only : imprimirCondicoesIniciais, imprimirSolucaoNoTempo,PRINT_DXINFO
        use mLeituraEscrita,   only : isat,ipres, escreverArqParaviewIntermed,PRINT_VTK,nprtsat,nprtpres
        use mLeituraEscrita,   only : iflag_sat,iflag_tipoPrint,ifsat_out,iflag_pres,ifpres_out
        use mMalha,            only : nsd, numel, numelReserv, numLadosReserv, nen, numLadosElem, numnp
        use mMalha,            only : x, xc, listaDosElemsPorNo, conecNodaisElem, conecLadaisElem
        use mPropGeoFisica,    only : phi, hx,hy,hz,nelx, nely, nelz,perm, permkx,permky,permkz, phi0
        use mPropGeoFisica,    only : lerPropriedadesFisicas, iflag_prod, TPRT_PROD, DTPRT_PROD, NP_RAND_PROD
        use mTransporte,       only : transport, satElem, satElemAnt
        use mHidrodinamicaRT,  only : hidroGeomecanicaRT,  pressaoElem, velocLadal, pressaoElemAnt, vc, velocConst
        USE MCMC
!
        implicit none
!
!.... solution driver program 
!
!
        integer :: i, nel
        real*8             :: t1, t2, dif, tol
        character(21)      :: label
!
        tol=1e-6
        ttv=0.d00
        ttp=0.d00
        tts=0.d00
        tempoSolverVel=0.d0
        ligarBlocosHetBeta=.false.
!
!.... le dados estocasticos
!    
        call lerPropriedadesFisicas()
!
!.... imprime condicoes iniciais
! 
        call imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, permkx, satElem)
        IF(NIFLAG_CONC.EQ.1) call MONITOR(filecout,pcondc,npcondc,elem_condc,satElem,tTransporte)
        IF(NIFLAG_PRES.EQ.1) call MONITOR(filepout,pcondp,npcondp,elem_condp,pressaoElem,tTransporte)
        IF(NIFLAG_PROD.EQ.1) call MONITORPROD(fileprout,npcondpr,velocLadal,satElem,perm,perm,perm,tTransporte)      
!
!.... tempo de simulacao para cada bloco do transporte
!
        dtBlocoTransp=tt/nvel
!
!-----------------------------------------------------------------------
!
!     SIMULACAO
!
!-----------------------------------------------------------------------
!
        DO NNP=1,NVEL

           write(*,"(//,'tempo inicial: ',f10.5,5x,'Passo: ',i5, 2x,'de ',i5,/)") tTransporte,nnp,nvel
! 
!******************* HIDRODINAMICA********************
!
           print*, "calculando a HIDRODINAMICA"
           call hidroGeomecanicaRT(satElem, satElemAnt)
!
           pressaoElemAnt=pressaoElem(1,:)
!***** simplificacao para velocidade constante igual a 1 ****
!           call velocConst()
! 
!******************* TRANSPORTE **********************
!

           call timing(t1) 
           print*, "calculando o TRANSPORTE:"
           call transport(velocLadal)
           call timing(t2)
#ifdef mostrarTempos
           write(*,'(2(a,f10.5))') 'tempo de parede, transport_ =', t2 - t1
#endif
           tts=tts+(t2-t1)

!
!.... imprime solucao intermediaria no tempo
! 
           call imprimirSolucaoNoTempo(pressaoElem, velocLadal, tTransporte)
!
   end do ! nnp=1,nvel
!
! imprime a solucao final
      if(iflag_sat==1)then
         if(iflag_tipoPrint==1) then
            call gerarLabel(label,tTransporte)
            call escreverArqParaviewIntermed(isat, satElem, ndofV, &
                 numel, trim(label), len(trim(label)))
         endif
         if(iflag_tipoPrint==2) then
            call PRINT_VTK(isat,satElem,tTransporte,ifsat_out,nen,ndofV,nprtsat)
         end if
      endif
      if(iflag_pres==1) then 
         if(iflag_tipoPrint==0) then
!            call prt  (nsd,numelReserv,tTransporte,pressaoElem,ipres)
         end if
         if(iflag_tipoPrint==1) then 
!                  call gerarLabel(label,tTransporte)
!                  call escreverArqParaviewIntermed(ipres, pressaoElem, ndofV, &
!                       numelReserv, trim(label), len(trim(label)))
!!                call paraview_escalarPorElementoTransiente(numel,u,nprint,isatTransiente)
!                  qtdImpSat=qtdImpSat+1
         end if
         if(iflag_tipoPrint==2) then
            call PRINT_VTK(ipres,pressaoElem,tTransporte,ifpres_out,nen,ndofV,nprtpres)
         end if
      end if

!
 1010 format('###########',/,'Fim da realizacao:', &
     & i5,/,'###########',/)

      write(*,*) "Tempo total da velocidade=", ttv
      write(*,*) "      Montagem velocidade=", tmVel
      write(*,*) "      Solver   velocidade=", tempoSolverVel
      write(*,*) "Tempo total da pressao   =", ttp
      
      write(*,*) "    (Tempo total da hidrodinamica = ", ttp+ttv, ")"

      write(*,*) "Tempo total da saturacao =", tts

      print*, " "
      print*, "TEMPO TOTAL DE EXECUCAO=", ttv+ttp+tts

      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine alocarMemoria()
      use mGlobaisArranjos,  only: uTempoN, mat, grav, beta
      use mGlobaisEscalares, only: ndofP, ndofV, nlvectP, nlvectV
      use mMalha,            only: nsd, numel, numelReserv, numnp, numLadosReserv, nen, numLadosElem, x, xc
      use mMalha,            only: listaDosElemsPorNo, listaDosElemsPorFace, conecNodaisElem, conecLadaisElem
      use mAlgMatricial,     only: idVeloc, idDesloc, lmV
      use mTransporte,       only: satElemAnt, satElem
      use mHidrodinamicaRT,  only:pressaoElem, velocLadal, fVeloc, pressaoElemAnt, vc, ve
!
      implicit none
!
!campos
      allocate(pressaoElem(ndofP,numelReserv));      pressaoElem    = 0.0d0
      allocate(pressaoElemAnt(numelReserv));         pressaoElemAnt = 0.0d0 
      allocate(velocLadal (ndofV,numLadosReserv));   velocLadal     = 0.0d0
      allocate(ve (nsd,nen,numelReserv));            ve             = 0.0d0
      allocate(vc (nsd,numelReserv));                vc             = 0.0d0
      allocate(satElem(numelReserv));                satElem        = 0.0d0
      allocate(satElemAnt(numelReserv));             satElemAnt     = 0.0d0
!malha
      allocate(x(nsd,numnp));                        x = 0.0d0
      allocate(xc(nsd,numel));                       xc= 0.0d0
      allocate(conecNodaisElem(nen,numel));          conecNodaisElem=0
      allocate(conecLadaisElem(numLadosElem,numelReserv)); conecLadaisElem=0
      allocate(listaDosElemsPorNo(nen,numnp));       listaDosElemsPorNo=0
      allocate(listaDosElemsPorFace(numLadosElem,numLadosReserv)); listaDosElemsPorFace=0

!contorno
      allocate(idVeloc(ndofV,numLadosReserv));              idVeloc  = 0
      allocate(lmV(ndofV,numLadosElem,numelReserv)); lmV=0
!
      if (nlvectV.ne.0)  then
         allocate(  fVeloc(ndofV,numLadosReserv))
         fVeloc   = 0.0d0
      endif

!material
      allocate(mat(numel)); mat=0.d0

!gravidade
      allocate(grav(3)); grav=0.d0

!tempo
      allocate(uTempoN(numel))


      allocate(beta(numel))


      end subroutine
!
!**** new **********************************************************************
!
      subroutine topologiaMalhaSistEquacoes(NALHSV, NEQV)
       use mGlobaisEscalares
       use mGlobaisArranjos
       use mAlgMatricial,   only: lmV, idVeloc, idiagV, meanbwV, nedV
       use mAlgMatricial,   only: colht,diag,numCoefPorLinhaVel, numCoefPorLinhaGeo
       use mMalha,          only: numel,numnp,nsd,numLadosElem,numLadosReserv,nen,numelReserv
       use mMalha,          only: x, xc, conecNodaisElem, conecLadaisElem, local
       use mMalha,          only: listaDosElemsPorNo, listaDosElemsPorFace, criarListaVizinhos
       use mMalha,          only: genel, genelFaces,  formlm, renumerarMalha
       use mPropGeoFisica,  only: hx,hy,hz,nelx,nely,nelz, REGION, calcdim
       use mPropGeoFisica,  only: nelxReserv, nelyReserv, nelzReserv
       use mLeituraEscrita, only: iin, iecho, nprint, prntel
       use mHidrodinamicaRT, only: fVeloc
       use mSolucoesExternas, only: criarPonteirosMatEsparsa_CRS
!  
      implicit none
!
!.... program to set storage and call tasks for  
!        an one dimensional conveccion-diffusion problem   
!             in non  dimensional form     
!             displacement formulation
!
!
      integer :: NALHSV, NEQV
      integer :: i, m, n, na
      real*8, allocatable :: xl(:,:) 
      real*8 :: t1, t2
!
!
!.... calculate hx, hy, hz
!
       call calcdim(nsd,numnp,x)
       print*, "numnp", numnp
       WRITE(*,'(3(a,f12.5))') " hx=", hx, ", hy=", hy, ", hz=", hz
!
!.... set element parameters
!
      allocate(npar(numParElem))
      read(iin,'(i10,14i10)') (npar(i),i=1,numParElem)
      ntype  = npar( 1)
      numat  = npar( 2)
      nen    = npar( 3)
      nicode = npar( 4)
      if(nsd==2)nrowsh = 3
      if(nsd==3)nrowsh = 4
      nedV   = ndofV
      nprint = 0
      NNP = 0

      tmGeo=0.0
      tmVel=0.0

      if (nicode.eq.0) nicode=nen
      npint   = nicode 
 !
!....... set memory pointers
! 
      allocate(c(6,numat)); c=0.d0
!
!
!.... input element data ('input___')

      write(iecho,1000) ntype,numel,numat,nen,npint
!
!      read material properties
!
      do 400 n=1,numat
      if (mod(n,50).eq.1) write(iecho,4000) numat
      read (iin,  5000) m,(c(i,m),i=1,3)
      write(iecho,6000) m,(c(i,m),i=1,3)
  400 continue
!
!     constant body forces
!
      read  (iin,  7000) (grav(i),i=1,3)
      write (iecho,8000) (grav(i),i=1,3)
!
!    generation of conectivities
!
      call genel(conecNodaisElem,mat,nen,iin)
      if(nsd==2)call genel(conecLadaisElem,mat,numLadosELem,iin)
      if(nsd==3)call genelFaces(conecLadaisElem, numLadosElem, nelxReserv, nelyReserv, nelzReserv, numelReserv,  iin) 
!
      call criarListaVizinhos(nen         ,numnp   ,numel,conecNodaisElem,listaDosElemsPorNo  )
      call criarListaVizinhos(numladosElem,numLadosReserv,numelReserv,conecLadaisElem,listaDosElemsPorFace)

!       do i=1, numLadosReserv
!          write(*,'(i5,a,6i5)'), i, ",",listaDosElemsPorFace(1:6,i)
!       enddo

!       stop

#ifdef debug
      if(iprtin.eq.0) then
         call prntel(mat,conecNodaisElem,nen,         numel,1)
         call prntel(mat,conecLadaisElem,numLadosELem,numel,2)
      end if
#endif
!
!     generation of lm array
!
      call formlm(idVeloc,conecLadaisElem,lmV,nedV,nedV,numLadosElem,numelReserv)
!
!     renumera a malha da velocidade
!
!       if((nelx>nely).and.(nsd==2))call renumerarMalha(ndofV, numLadosElem, numLados, idVeloc, lmV, conecLadaisElem, fVeloc, listaDosElemsPorFace, x, nelx, nely, nelz)
!
!     modification of idiag array
!
      call colht(idiagV,lmV,nedV,numLadosElem,numelReserv,neqV)
!
#ifndef withcrs
      call diag(idiagV,neqV,nalhsV)
#endif

      meanbwV = nalhsV/neqV
      write(iecho,9000) 'Calculo das velocidades nodais', neqV, nalhsV, meanbwV, 8.0*nalhsV/1000/1000

!     calculo de xc (centros dos elementos)
      allocate(xl(nsd,nen)); xl=0.0
      do i=1,numel
         call local(conecNodaisElem(1,i),x,xl,nen,nsd,nsd)
         xc(1,i) = sum(xl(1,1:nen))/nen
         xc(2,i) = sum(xl(2,1:nen))/nen
         if(nsd==3)xc(3,i) = sum(xl(3,1:nen))/nen
      end do

     !Escolha do solver
      optSolver='skyline'
#ifdef withlapack
      optSolver='lapack' 
#endif
#ifdef withumfpack
      optSolver='umfpack'
      simetriaVel=0 ! 0 para não simetrico, unica versão disponivel para o UMFPACK
#endif
#ifdef withpardiso
      optSolver='pardiso'
      simetriaVel=.true. ! 0 para não simetrico ou 1 para simetrico
#endif
      
      print*, "solver escolhido para a solucao do sistema e equacoes: ", optSolver

#ifdef withcrs
      call timing(t1)
      if(nsd==2) numCoefPorLinhaVel=7
      if(nsd==3) numCoefPorLinhaVel=11

      call criarPonteirosMatEsparsa_CRS(nsd, ndofV, neqV, numCoefPorLinhaVel, &
             conecLadaisElem, listaDosElemsPorFace, idVeloc, numLadosReserv, numLadosElem, numLadosElem, na, simetriaVel)

      nalhsV=na
      call timing(t2)
#ifdef mostrarTempos
      print*, "Tempo de pre-processamento da Velocidade para solver externo", t2-t1
#endif
     
#endif
      print*, "NALHSV=", NALHSV
    
      return
!
 1000 format(//,&
     ' two/three-n o d e    e l e m e n t s ',//,5x,&
     ' element type number . . . . . . . . . . . (ntype ) = ',i10,//5x,&
     ' number of elements  . . . . . . . . . . . (numel ) = ',i10,//5x,&
     ' number of element material sets . . . . . (numat ) = ',i10,//5x,&
     ' number of element nodes . . . . . . . . . (nen   ) = ',i10,//5x,&
     ' number of integration points. . . . . . . (npint  ) = ',i10)
 4000  format(///,&
     ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
     ' number of material sets . . . . . . . . . . (numat ) = ',i10///,2x,'set',4x,'Kx ',&
       10x,'Ky',10x,'Kz')
 5000  format(i10,5x,5f10.0)
 6000  format(2x,i3,1x,5(1x,1pe11.4))
 7000 format(8f10.0)
 8000 format(///,&
     ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
     ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
     ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
     ' exemplo 3............................ = ',      1pe15.8,//)
 9000 format('1',a///&
     ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
     ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
     ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
     ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
     ' memoria necessaria para a matriz do sistema  (Mbytes)  = ',e10.2)
!
      end  subroutine 
