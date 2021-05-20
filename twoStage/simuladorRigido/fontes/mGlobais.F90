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
!     *                 GALERKIN METHOD                          *
!     *                                                          *
!     *                                                          *
!     ************************************************************
!
      module mGlobaisArranjos

        integer, allocatable :: npar(:)
        real*8               :: etime(6)
        character*4          :: title(20)

        real*8,  allocatable :: uTempoN(:)

        integer, allocatable :: mat(:)
        real*8, allocatable  :: grav(:), bf(:,:), c(:,:)

        real*8, allocatable  :: beta(:)

      end module mGlobaisArranjos


      module mGlobaisEscalares

        integer :: ntype
        integer :: exec,iprtin
        integer :: numParElem=15
        integer :: ndofP, ndofV, ndofD
        integer :: nlvectP, nlvectV, nlvectD
        real*8, parameter  :: zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0
        real*8, parameter  :: four=4.0d0, five=5.0d0, six=6.0d0
        real*8, parameter  :: pt1667=0.1666666666666667d0, pt25=0.25d0, pt5=0.5d0
        real*8  :: coef
        integer :: nRK, ordemRK
        integer :: optCC
        character(len=10) :: optSolver
        logical :: simetriaVel, simetriaGeo

        integer :: numat
        integer :: nrowsh,nicode,npint

      integer :: nvel,nnp,ns
      integer :: NITGEO, NCREEP, NLOOPS, NUMDX
      real(8) :: tzero,tTransporte,tc,tt,dtBlocoTransp
      real(8) :: tempoNucTrans
      real*8  :: ttv, ttp, tts, tempoSolverVel
      real*8  :: tmVel, tmGeo, tsGeo
      real*8  :: tgeoFase1, tgeoFase2, tgeoFase3, tgeoFase4, ttgeo

      integer :: geomech
      logical :: creep, ligarBlocosHetBeta=.false.
      integer :: iflag_beta
      logical :: novaMalha

      end module mGlobaisEscalares

