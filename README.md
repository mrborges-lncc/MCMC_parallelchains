# MCMC_parallelchains

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
     IF(GERATIPO(I).EQ.6)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN/in'
     END IF
     IF(GERATIPO(I).EQ.7)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN1/in'
     END IF
     IF(GERATIPO(I).EQ.8)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR MVN >>>'
        FILE_INFIELD(I) = '../gera_KL/MVN2/in'
     END IF
     IF(GERATIPO(I).EQ.13)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3D/in'
     END IF
     IF(GERATIPO(I).EQ.14)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3D_2/in'
     END IF
     IF(GERATIPO(I).EQ.15)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3D_3/in'
     END IF
     IF(GERATIPO(I).EQ.10)THEN
        IF(NK.EQ.0)WRITE(*,*)'<<< GERADOR FORTRAN CAMPOS 3D COND >>>'
        FILE_INFIELD(I) = '../gera_KL/FORTRAN_KL3DCOND/in'
     END IF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NK.EQ.0)WRITE(*,*)'SIMULADOR UTILIZADO PARA O PROBLEMA DE ESCOAMENTO'
  IF(NK.EQ.0)WRITE(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  IF(NSIMUL.EQ.0)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< BLACKBOX >>>'
  END IF
  IF(NSIMUL.EQ.1)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< SIMULADOR >>>'
  END IF
  IF(NSIMUL.EQ.2)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< simuadorRigido >>>'
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
  IF(NSIMUL.EQ.6)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< TWOPHASEFLOW MATLAB >>>'
  END IF
  IF(NSIMUL.EQ.7)THEN
     IF(NK.EQ.0)WRITE(*,*)'<<< TWOPHASEFLOW OCTAVE >>>'
  END IF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     
     
     
     
     
     
     
     
     
     
     
