MODULE MCMC
  INTEGER :: NIFLAG_PROD,NIFLAG_CONC,NIFLAG_PRES
  REAL(8) :: TPRT_PRODM,DTPRT_PRODM
  REAL(8) :: TPRT_PRESM,DTPRT_PRESM,TPRT_CONCM,DTPRT_CONCM
  REAL(8), ALLOCATABLE :: PCONDC(:,:),PCONDP(:,:)
  INTEGER, ALLOCATABLE :: ELEM_CONDC(:),FACE_PR(:,:),ELEM_CONDP(:)
  CHARACTER(LEN=128) :: FILEPIN,FILEPRIN,FILECIN
  CHARACTER(LEN=128) :: FILEPOUT,FILEPROUT,FILECOUT
  INTEGER :: NCONDC,NCONDP,NCONDPR
  INTEGER :: NPCONDC,NPCONDP,NPCONDPR,NWELL
  INTEGER, ALLOCATABLE :: NELEMPROD(:,:)
  private :: f, g, h
  contains

SUBROUTINE INITMCMC()

  use mGlobaisEscalares, only: tt,tzero,NVEL
!
  IMPLICIT NONE
!
  CHARACTER(LEN=128) :: FNOME,FILEIN,FILEOUT,NAME,FLAG,FLAG1
  INTEGER :: ISTAT,IDFILE=188
  REAL(8) :: TEMPO, TOL=1E-6
!
  integer :: ierr
!
  FNOME = 'mcmc.in'
  FNOME = ADJUSTL(TRIM(FNOME))
!
  OPEN(UNIT=IDFILE,FILE=FNOME,STATUS='OLD',&
       ACTION='READ',IOSTAT=ISTAT)
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
     STOP
  ELSE
     WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FNOME))
  END IF
!
  WRITE(*,*)'############################################################'
!
!.... Saida da concentracao
!
  FLAG="# concentracao"
  READ(IDFILE,"(a)")FLAG1
  IF(TRIM(FLAG).NE.TRIM(FLAG1)) STOP 1234
  READ(IDFILE,*) NIFLAG_CONC
  READ(IDFILE,"(a)") FILECOUT
  FILECOUT=ADJUSTL(TRIM(FILECOUT))
  WRITE(*,*)'ARQUIVO DADOS PARA A CONCENTRACAO...: ',ADJUSTL(TRIM(FILECOUT))
  READ(IDFILE,*)NCONDC
  READ(IDFILE,"(a)") FILECIN
  FILECIN=TRIM(FILECIN)
  WRITE(*,*)'ARQUIVO DE SAIDA CONCENTRACAO.......: ',ADJUSTL(TRIM(FILECIN))
!
!.... Saida da pressao
!
  FLAG="# pressao"
  READ(IDFILE,"(a)")FLAG1
  IF(TRIM(FLAG).NE.TRIM(FLAG1)) STOP 1235
  READ(IDFILE,*) NIFLAG_PRES
  READ(IDFILE,"(a)") FILEPOUT
  FILEPOUT=ADJUSTL(TRIM(FILEPOUT))
  WRITE(*,*)'ARQUIVO DADOS PARA A PRESSAO........: ',ADJUSTL(TRIM(FILEPOUT))
  READ(IDFILE,*)NCONDP
  READ(IDFILE,"(a)") FILEPIN
  FILEPIN=TRIM(FILEPIN)
  WRITE(*,*)'ARQUIVO DE SAIDA PRESSAO............: ',ADJUSTL(TRIM(FILEPIN))
!
!.... Saida da producao
!
  FLAG="# producao"
  READ(IDFILE,"(a)")FLAG1
  IF(TRIM(FLAG).NE.TRIM(FLAG1)) STOP 1236
  READ(IDFILE,*) NIFLAG_PROD
  READ(IDFILE,"(a)") FILEPROUT
  FILEPROUT=ADJUSTL(TRIM(FILEPROUT))
  WRITE(*,*)'ARQUIVO DADOS PARA A PRODUCAO.......: ',ADJUSTL(TRIM(FILEPROUT))
  READ(IDFILE,*)NCONDPR
  READ(IDFILE,"(a)") FILEPRIN
  FILEPRIN=TRIM(FILEPRIN)
  WRITE(*,*)'ARQUIVO DE SAIDA PRODUCAO...........: ',ADJUSTL(TRIM(FILEPRIN))
!
  CLOSE(IDFILE)
!
  IF(NIFLAG_CONC.NE.1) NCONDC=1
  IF(NIFLAG_PRES.NE.1) NCONDP=1
  IF(NIFLAG_PROD.NE.1) NCONDPR=1
!
  DTPRT_CONCM = tt/REAL(NCONDC)
  DTPRT_PRESM = tt/REAL(NCONDP)
  DTPRT_PRODM = tt/REAL(NCONDPR)
!
  TPRT_PRODM = DTPRT_PRODM
  TPRT_PRESM = DTPRT_PRESM
  TPRT_CONCM = DTPRT_CONCM
!
  WRITE(*,*)'############################################################'
  WRITE(*,*)'###########  NUMERO DE MONITORAMENTOS NO TEMPO  ############'
  WRITE(*,*)'- CONCENTRACAO........:',NCONDC, 'DT.......:',DTPRT_CONCM
  WRITE(*,*)'- PRESSAO.............:',NCONDP, 'DT.......:',DTPRT_PRESM
  WRITE(*,*)'- PRODUCAO............:',NCONDPR,'DT.......:',DTPRT_PRODM
  WRITE(*,*)'############################################################'
  CALL LEITURA()
!
END SUBROUTINE INITMCMC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE IREADSTAT(FLAG,IFILE,N)
!
  IMPLICIT NONE
  INTEGER :: IFILE,N
  CHARACTER(len=128) :: FLAG,FLAG1
!
  READ(IFILE,"(a)") FLAG1
  IF(TRIM(FLAG).EQ.TRIM(FLAG1)) THEN
     READ(IFILE,*) N
  ELSE
     WRITE(*,*) "Erro na leitura de ", FLAG
     STOP
  END IF
!
END SUBROUTINE IREADSTAT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LEITURA()
!
  use mMalha,            only: numelReserv, conecNodaisElem, NSD
  use mMalha,            only: conecLadaisElem, listaDosElemsPorFace
  IMPLICIT NONE
  INTEGER :: IFILE,I,NEL,J,K
  INTEGER, DIMENSION(6) :: N
  CHARACTER(LEN=128) :: FNOME
  INTEGER :: ISTAT,IDFILE=189
!
  WRITE(*,*)'############################################################'
  WRITE(*,*)'###########  NUMERO DE PONTOS DE MONITORAMENTO  ############'
!
! LEITURA DOS DADOS DE CONCENTRACAO
! 
  FNOME = FILECIN
  FNOME = ADJUSTL(TRIM(FNOME))
!
  IF(NIFLAG_CONC.EQ.1)THEN
     OPEN(UNIT=IDFILE,FILE=FNOME,STATUS='UNKNOWN',&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
        STOP
     ELSE
        WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FNOME))
     END IF
     READ(IDFILE,*)NPCONDC
     WRITE(*,*)'- CONCENTRACAO........:',NPCONDC
     ALLOCATE(PCONDC(3,NPCONDC))
     ALLOCATE(ELEM_CONDC(NPCONDC))
     PCONDC = 0.0
     ELEM_CONDC = 0
!
     DO I=1,NPCONDC
        READ(IDFILE,'(3E16.7)')PCONDC(1,I),PCONDC(2,I),PCONDC(3,I)
        WRITE(*,'(3E10.3)')PCONDC(1,I),PCONDC(2,I),PCONDC(3,I)
     END DO
  END IF
!
  CALL MONITORPOSITION(PCONDC,NPCONDC,ELEM_CONDC)
  CLOSE(IDFILE)
!
  WRITE(*,*)'############################################################'
!
! LEITURA DOS DADOS DE PRESSAO
! 
  FNOME = FILEPIN
  FNOME = ADJUSTL(TRIM(FNOME))
!
  IF(NIFLAG_PRES.EQ.1)THEN
     OPEN(UNIT=IDFILE,FILE=FNOME,STATUS='UNKNOWN',&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
        STOP
     ELSE
        WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FNOME))
     END IF
     READ(IDFILE,*)NPCONDP
     WRITE(*,*)'- PRESSAO.............:',NPCONDP
     ALLOCATE(PCONDP(3,NPCONDP))
     ALLOCATE(ELEM_CONDP(NPCONDP))
     PCONDP = 0.0
     ELEM_CONDP = 0
!
     DO I=1,NPCONDP
        READ(IDFILE,'(3E16.7)')PCONDP(1,I),PCONDP(2,I),PCONDP(3,I)
        WRITE(*,'(3E10.3)')PCONDP(1,I),PCONDP(2,I),PCONDP(3,I)
     END DO
  END IF
  CALL MONITORPOSITION(PCONDP,NPCONDP,ELEM_CONDP)
  CLOSE(IDFILE)
!
  WRITE(*,*)'############################################################'
!
! LEITURA DOS DADOS DE PRESSAO
! 
  FNOME = FILEPRIN
  FNOME = ADJUSTL(TRIM(FNOME))
!
  IF(NIFLAG_PROD.EQ.1)THEN
     OPEN(UNIT=IDFILE,FILE=FNOME,STATUS='UNKNOWN',&
          ACTION='READ',IOSTAT=ISTAT)
     IF(ISTAT.NE.0)THEN
        WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
        STOP
     ELSE
        WRITE(*,*)'READING FILE:',ADJUSTL(TRIM(FNOME))
     END IF
     READ(IDFILE,*)NWELL
     WRITE(*,*)'NUMERO de POCOS.......:',NWELL
     READ(IDFILE,*)NPCONDPR
     WRITE(*,*)'- PRODUCAO............:',NPCONDPR
     ALLOCATE(FACE_PR(NWELL,NPCONDPR))
     ALLOCATE(NELEMPROD(NWELL,NPCONDPR))
     FACE_PR   = 0
     NELEMPROD = 0
!
     DO J=1,NWELL
        DO I=1,NPCONDPR
           READ(IDFILE,*)FACE_PR(J,I)
        END DO
     END DO
  END IF
  CLOSE(IDFILE)
!
  WRITE(*,*)'############################################################'
!
  DO K=1,NWELL
     DO I=1,NPCONDPR
        DO J=1,numelReserv
           IF(NSD.EQ.3)THEN
              IF(FACE_PR(K,I).EQ.conecLadaisElem(1,J).OR.FACE_PR(K,I).EQ.conecLadaisElem(2,J).OR. &
                   FACE_PR(K,I).EQ.conecLadaisElem(3,J).OR.FACE_PR(K,I).EQ.conecLadaisElem(4,J).OR. &
                   FACE_PR(K,I).EQ.conecLadaisElem(5,J).OR.FACE_PR(K,I).EQ.conecLadaisElem(6,J))THEN
                 NELEMPROD(K,I) = J
              END IF
           END IF
           IF(NSD.EQ.2)THEN
              IF(FACE_PR(K,I).EQ.conecLadaisElem(1,J).OR.FACE_PR(K,I).EQ.conecLadaisElem(2,J).OR. &
                   FACE_PR(K,I).EQ.conecLadaisElem(3,J).OR.FACE_PR(K,I).EQ.conecLadaisElem(4,J))THEN
                 NELEMPROD(K,I) = J
              END IF
           END IF
        END DO
     END DO
  END DO
!
END SUBROUTINE LEITURA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MONITORPOSITION(PCOND,N,ELEM_COND)
  use mMalha,         only: x,NUMEL,conecNodaisElem,nsd,nen
  IMPLICIT NONE
!
  INTEGER :: N,J,I,NO,NEL
  REAL(8), DIMENSION(3,*) :: PCOND
  REAL(8) :: X1,X2,Y1,Y2,Z1,Z2,H1,H2
  INTEGER, DIMENSION(N),INTENT(OUT) :: ELEM_COND
!
  WRITE(*,*)'########################################'
  WRITE(*,*)'### MONITORAMENTO ######################'
  WRITE(*,*)'### NUMERO DE PONTOS:',N
  WRITE(*,*)'########################################'
  NO = 0
  IF(nsd.eq.2)THEN
     DO 200 I=1,N
        DO NEL=1,NUMEL
           X1=x(1,conecNodaisElem(1,NEL))
           X2=x(1,conecNodaisElem(2,NEL))
           Y1=x(2,conecNodaisElem(1,NEL))
           Y2=x(2,conecNodaisElem(4,NEL))
           IF(PCOND(1,I).LT.X2)THEN
              IF(PCOND(1,I).GE.X1)THEN
                 IF(PCOND(2,I).LT.Y2)THEN
                    IF(PCOND(2,I).GE.Y1)THEN
                       ELEM_COND(I)=NEL
                       NO=NO+1
                       GO TO 200
                    END IF
                 END IF
              END IF
           END IF
        END DO
200     CONTINUE
  END IF
!
  IF(nsd.eq.3)THEN
     DO 300 I=1,N
        DO NEL=1,NUMEL
           X1=x(1,conecNodaisElem(1,NEL))
           X2=x(1,conecNodaisElem(2,NEL))
           Y1=x(2,conecNodaisElem(1,NEL))
           Y2=x(2,conecNodaisElem(4,NEL))
           Z1=x(3,conecNodaisElem(1,NEL))
           Z2=x(3,conecNodaisElem(5,NEL))
           IF(PCOND(1,I).LT.X2)THEN
              IF(PCOND(1,I).GE.X1)THEN
                 IF(PCOND(2,I).LT.Y2)THEN
                    IF(PCOND(2,I).GE.Y1)THEN
                       IF(PCOND(3,I).LT.Z2)THEN
                          IF(PCOND(3,I).GE.Z1)THEN
                             ELEM_COND(I)=NEL
                             NO=NO+1
                             GO TO 300
                          END IF
                       END IF
                    END IF
                 END IF
              END IF
           END IF
        END DO
300     CONTINUE
  END IF
  !
  IF(NO.NE.N)THEN
     WRITE(*,*)'##################################'
     WRITE(*,*)'PROBLEMA NA LEITURA DOS DADOS MCMC'
     WRITE(*,*)'##################################'
     STOP
  END IF
!
END SUBROUTINE MONITORPOSITION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MONITOR(FNOME,PCOND,N,ELEM_COND,U,T)
  use mMalha,         only: nsd
  IMPLICIT NONE
!
  INTEGER :: N,INFILE
  REAL(8) :: T,TOL=1E-09
  REAL(8), DIMENSION(nsd,*) :: PCOND
  REAL(8), DIMENSION(*) :: U
  CHARACTER(LEN=128) :: FNOME
  INTEGER :: ISTAT,I,ZERO=0D0
  INTEGER, DIMENSION(N)   :: ELEM_COND
  CHARACTER(len=128) :: NAME,FILE_IN
  CHARACTER(len=4)   :: EXT
  CHARACTER(len=5)   :: C
!
  INFILE = 163
  WRITE(C,113)ZERO
  C=ADJUSTL(C)
  FILE_IN = FNOME
  EXT='.dat'
  NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
  NAME=ADJUSTL(TRIM(NAME))
!  write(*,*)NAME
!
  OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
       ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
     STOP
  END IF
!
!  write(*,*)'##### MONITORANDO #####'
  DO I=1,N
     IF(U(ELEM_COND(I)).LT.TOL) U(ELEM_COND(I))=0.0
  END DO
  WRITE(INFILE,400)T,(U(ELEM_COND(I)),I=1,N)
!
  CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(100(e12.4))
END SUBROUTINE MONITOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MONITORPROD(FNOME,N,U,SAT,XK,YK,ZK,T)
  use mMalha,         only: nsd,numLadosReserv
  use mPropGeoFisica,  only: xlw,xlo,xlt,gf1,rhoo,rhow
  IMPLICIT NONE
!
  INTEGER :: N,INFILE,NO,FACE
  REAL(8) :: T,ZERO=0D0
  REAL(8), ALLOCATABLE :: PROD(:)
  REAL(8) :: SO
  CHARACTER(LEN=128) :: FNOME
  INTEGER :: ISTAT,I,J,NZERO=0D0
  CHARACTER(len=128) :: NAME,FILE_IN
  CHARACTER(len=4)   :: EXT
  CHARACTER(len=5)   :: C
  real(8), DIMENSION(*) :: U,XK,YK,ZK
  REAL(8), DIMENSION(*) :: SAT
!
  ALLOCATE(PROD(NWELL))
!
  INFILE = 164
  WRITE(C,113)NZERO
  C=ADJUSTL(C)
  FILE_IN = FNOME
  EXT='.dat'
  NAME=TRIM(FILE_IN)//TRIM(C)//TRIM(EXT)
  NAME=ADJUSTL(TRIM(NAME))
!  write(*,*)NAME
!
  OPEN(UNIT=INFILE,FILE=NAME,STATUS='UNKNOWN',&
       ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
  IF(ISTAT.NE.0)THEN
     WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',FNOME
     STOP
  END IF
  !
!  write(*,*)'##### MONITORANDO PRODUCAO #####'
  DO J=1,NWELL
     PROD(J) = 0.0
     DO I=1,N
        FACE    = FACE_PR(J,I)
        NO      = NELEMPROD(J,I)
        SO      = xlo(SAT(NO))/xlt(SAT(NO))
        PROD(J) = PROD(J) + SO*DABS(U(FACE))
     END DO
  END DO
  !
  WRITE(INFILE,400)T,(PROD(J),J=1,NWELL)
!
  CLOSE(INFILE)
!
113 FORMAT(I5)
400 FORMAT(3(e15.8,2x))
END SUBROUTINE MONITORPROD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end MODULE MCMC
