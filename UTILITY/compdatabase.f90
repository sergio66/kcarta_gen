! ifort -extend-source 132 compdatabase.f;  mv a.out compdatabase.x                  ifort compiler

! this program reads in a list of compressed data files and outputs 
! this program used to be called compdatabase97.f, but now since the 
! compressed database has both the main gases and the cross section
! gases, we just call it compdatabase.f

! this program is modified version of compdatabase.f since it goes thru the
! subdirectories and produces a 4 column format, 
!   GASID  START-FREQ   STOP-FREQ  TAG
! where TAG tells us if this is a 0.001, 0.0025 or 0.005 spacing, 
! depending on the prefix of the file.
!    TAG = 1  qXXX-gI.dat     0.001  wavenumber spacing   205-605
!    TAG = 2  rXXX-gI.dat     0.0025 wavenumber spacing   0605-2805 cm-1
!    TAG = 3  sXXX-gI.dat     0.0025 wavenumber spacing   2830-3305 cm-1

! %% orig
!      PARAMETER (kNumkCompT=188)
!      DATA kaTag        / 1,   2,   3,   4,   5,   6,   7,  8,   9   /
!      DATA kaNumkComp   / 08, 89,  09,  20,  24,  16,  20,  1,   1   /
!      DATA kaCTag       /'q', 'r', 'm', 'n', 'v', 'w', 's', 'k', 'p' /
!      DATA kaMinFr      /0500.00000,   0605.00000,    4050.00000,  5000.00000,
!     $                   8000.00000,   14000.00000,   2830.00000, 
!     $                   0200.00000,   0355.00000/
!      DATA kaMaxFr     /0605.00000,   2830.00000,    4950.00000,   8000.00000,
!     $                   14000.00000,  22000.00000,   3305.00000,
!     $                    0355.00000,   0500.00000/
!      DATA kaFrStep  /1.500000e-3,  2.500000e-3,   1.000000e-2,  1.500000e-2,
!     $                 2.500000e-2,  5.000000e-2,   2.500000e-3,
!     $                 0.500000e-3,  1.000000e-3/

! Oct 2019 : changed 'q' to be 500 to 880  at 0.0005 cm-1 res
! Oct 2019 : changed 'r' to be 805 to 2830 at 0.0025 cm-1 res

! check n_gas_wt_spectra.f for allowed gases <<<<<<<<<<<<<<-------------
!  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9 
!         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,X PLUS kSelf,kFor 
! check n_gas_wt_spectra.f for allowed gases <<<<<<<<<<<<<<-------------

! 1)file containing the GasID's and frequencies of each file comp1.param
! 2)file containing the summary of GasID's and frequencies comp.param in the
!   4 column format   GASID  START-FREQ   STOP-FREQ  TAG
! comp.param is used by kcarta.x

! to use this code, go to the subdir containing the database
! eg cd /salsify/production3/motteler/sergio/kcarta/h2o.ieee-be
! do ls -1 *.dat > & /salsify/users/sergio/KCARTA/UTILITY/compdatabase
! then run this program
! you will have to do it separately for the water database and the "rest
! of gases" database; combine the two and then separate out the minor or
! cross section gases. eg put the results into 2 files, 
! '../DATA/General/xsec107.param' and '../DATA/General/comp107.param'

      IMPLICIT NONE

      include '../INCLUDE/TempF90/scatterparam.f90'

! the input file that contains the directory listing has all the files 
! that start with ""q" r" "s"  eg r605_g2.dat
      CHARACTER*80 caStr
      REAL raaArrF(kMaxPts,3),raaBlockF(kMaxPts,3)
      REAL raaArrG(kMaxPts,3),raaBlockG(kMaxPts,3)
      REAL raaArrH(kMaxPts,3),raaBlockH(kMaxPts,3)
      REAL raaArrJ(kMaxPts,3),raaBlockJ(kMaxPts,3)
      REAL raaArrK(kMaxPts,3),raaBlockK(kMaxPts,3)
      REAL raaArrP(kMaxPts,3),raaBlockP(kMaxPts,3)
      REAL raaArrQ(kMaxPts,3),raaBlockQ(kMaxPts,3)
      REAL raaArrR(kMaxPts,3),raaBlockR(kMaxPts,3)
      REAL raaArrS(kMaxPts,3),raaBlockS(kMaxPts,3)
      REAL raaArrM(kMaxPts,3),raaBlockM(kMaxPts,3)
      REAL raaArrN(kMaxPts,3),raaBlockN(kMaxPts,3)
      REAL raaArrO(kMaxPts,3),raaBlockO(kMaxPts,3)
      REAL raaArrV(kMaxPts,3),raaBlockV(kMaxPts,3)
      REAL raaArrU(kMaxPts,3),raaBlockU(kMaxPts,3)

      INTEGER iInt,iTag

      INTEGER iCntJ, iCntK, iCntM, iCntN, iCntP, iCntQ, iCntR, iCntS  
      INTEGER iCnt2J,iCnt2K,iCnt2M,iCnt2N,iCnt2P,iCnt2Q,iCnt2R,iCnt2S 
      INTEGER iCntO,  iCntV,  iCntU,  iCntF, iCntG, iCntH
      INTEGER iCnt2O, iCnt2V, iCnt2U, iCnt2F,iCnt2G,iCnt2H

! set GasID,start/stop freq WAY above what they can be!
! update this subroutine when you add on new bands
      CALL UPDATE_DoInit(iCntF,raaArrF,iCntG,raaArrG,iCntH,raaArrH, &
                         iCntJ,raaArrJ,iCntK,raaArrK,iCntP,raaArrP, &
                         iCntM,raaArrM,iCntN,raaArrN,iCntO,raaArrO, &
                         iCntQ,raaArrQ,iCntR,raaArrR,iCntS,raaArrS, &
                         iCntV,raaArrV,iCntU,raaArrU)

! do ls -1 compdatapath/*.dat> & compdatabase and then run this program
      OPEN(UNIT=10,FILE='compdatabase',FORM='FORMATTED',STATUS='OLD')
 11   READ(10,1000,END=100) caStr
 1000 FORMAT(A80)
!      print *,caStr

      IF (caStr(1:1).EQ. 'f') THEN
        CALL process(caStr,raaArrF,iCntF,0.5) !10000.0*0.00005
        END IF
      IF (caStr(1:1).EQ. 'g') THEN
        CALL process(caStr,raaArrG,iCntG,1.0) !10000.0*0.00010
        END IF
      IF (caStr(1:1).EQ. 'h') THEN
        CALL process(caStr,raaArrH,iCntH,1.5) !10000.0*0.00015
        END IF
      IF (caStr(1:1).EQ. 'j') THEN
        CALL process(caStr,raaArrJ,iCntJ,2.5) !10000.0*0.00025
        END IF
      IF (caStr(1:1).EQ. 'k') THEN
        CALL process(caStr,raaArrK,iCntK,5.0) !10000.0*0.0005
        END IF
      IF (caStr(1:1).EQ. 'p') THEN
        CALL process(caStr,raaArrP,iCntP,10.0) !10000.0*0.0010
        END IF
      IF (caStr(1:1).EQ. 'q') THEN
        !CALL process(caStr,raaArrQ,iCntQ,15.0) !10000.0*0.0015
        CALL process(caStr,raaArrQ,iCntQ,05.0) !10000.0*0.0005
        END IF
      IF (caStr(1:1).EQ. 'r') THEN
        CALL process(caStr,raaArrR,iCntR,25.0) !10000.0*0.0025
        END IF
      IF (caStr(1:1).EQ. 's') THEN
        CALL process(caStr,raaArrS,iCntS,25.0) !10000.0*0.0025
        END IF
      IF (caStr(1:1).EQ. 'm') THEN
        CALL process(caStr,raaArrM,iCntM,100.0) !10000.0*0.01
        END IF
      IF (caStr(1:1).EQ. 'n') THEN
        CALL process(caStr,raaArrN,iCntN,150.0) !10000.0*0.015
        END IF
      IF (caStr(1:1).EQ. 'o') THEN
        CALL process(caStr,raaArrO,iCntO,250.0) !10000.0*0.025
        END IF
      IF (caStr(1:1).EQ. 'v') THEN
        CALL process(caStr,raaArrV,iCntV,500.0) !10000.0*0.050
        END IF
      IF (caStr(1:1).EQ. 'u') THEN
        CALL process(caStr,raaArrU,iCntU,1000.0) !10000.0*0.100
        END IF
      GO TO 11

 100  CONTINUE     
      CLOSE(10)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! first sort array raaArr according to frequency indices
! then sort array raaArr according to iGasID
      IF (iCntF .GT. 0) THEN
        CALL dosort(raaArrF,iCntF,2)
        CALL dosort(raaArrF,iCntF,1)
        iTag = 02        !! spacing = 0.00005 cm-1
        END IF 

      IF (iCntG .GT. 0) THEN
        CALL dosort(raaArrG,iCntG,2)
        CALL dosort(raaArrG,iCntG,1)
        iTag = 04        !! spacing = 0.000010 cm-1
        END IF 

      IF (iCntH .GT. 0) THEN
        CALL dosort(raaArrH,iCntH,2)
        CALL dosort(raaArrH,iCntH,1)
        iTag = 06        !! spacing = 0.000015 cm-1
        END IF 

      IF (iCntJ .GT. 0) THEN
        CALL dosort(raaArrJ,iCntJ,2)
        CALL dosort(raaArrJ,iCntJ,1)
        iTag = 08        !! spacing = 0.000025 cm-1
        END IF 

      IF (iCntK .GT. 0) THEN
        CALL dosort(raaArrK,iCntK,2)
        CALL dosort(raaArrK,iCntK,1)
        iTag = 10        !! spacing = 0.00005 cm-1
        END IF 

      IF (iCntP .GT. 0) THEN
        CALL dosort(raaArrP,iCntP,2)
        CALL dosort(raaArrP,iCntP,1)
        iTag = 12        !! spacing = 0.0010 cm-1
        END IF 

      IF (iCntQ .GT. 0) THEN
        CALL dosort(raaArrQ,iCntQ,2)
        CALL dosort(raaArrQ,iCntQ,1)
        iTag = 15        !! spacing = 0.0015 cm-1, changed to 0.0005 cm-1 in Oct 2019
        END IF 

      IF (iCntR .GT. 0) THEN
        CALL dosort(raaArrR,iCntR,2)
        CALL dosort(raaArrR,iCntR,1)
        iTag = 20        !! spacing = 0.0025 cm-1
        END IF

      IF (iCntS .GT. 0) THEN
        CALL dosort(raaArrS,iCntS,2)
        CALL dosort(raaArrS,iCntS,1)
        iTag = 25        !! spacing = 0.0025 cm-1
        END IF

      IF (iCntM .GT. 0) THEN
        CALL dosort(raaArrM,iCntM,2)
        CALL dosort(raaArrM,iCntM,1)
        iTag = 30        !! spacing = 0.01 cm-1
        END IF

      IF (iCntN .GT. 0) THEN
        CALL dosort(raaArrN,iCntN,2)
        CALL dosort(raaArrN,iCntN,1)
        iTag = 35        !! spacing = 0.015 cm-1
        END IF

      IF (iCntO .GT. 0) THEN
        CALL dosort(raaArrO,iCntO,2)
        CALL dosort(raaArrO,iCntO,1)
        iTag = 40        !! spacing = 0.025 cm-1
        END IF

      IF (iCntV .GT. 0) THEN
        CALL dosort(raaArrV,iCntV,2)
        CALL dosort(raaArrV,iCntV,1)
        iTag = 50        !! spacing = 0.050 cm-1
        END IF

      IF (iCntU .GT. 0) THEN
        CALL dosort(raaArrU,iCntU,2)
        CALL dosort(raaArrU,iCntU,1)
        iTag = 55        !! spacing = 0.10 cm-1
        END IF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! first output all the file info that has been found, into comp1.param
! then summarize the sorted indices so that "blocks" are output, instead
! of individual file information, not comp.param

! the tags are stored in ../INCLUDE/pre_defined.param

      OPEN(UNIT=11,FILE='comp1.param',FORM='FORMATTED',STATUS='UNKNOWN')
      OPEN(UNIT=12,FILE='comp.param',FORM='FORMATTED',STATUS='NEW')

      IF (iCntF .GT. 0) THEN
        DO iInt=1,iCntF
          WRITE (11,2000) iInt,IFIX(raaArrF(iInt,1)),raaArrF(iInt,2),raaArrF(iInt,3),iTag
          END DO
        CALL doBlock(raaArrF,iCntF,raaBlockF,iCnt2F)
        DO iInt=1,iCnt2F
          WRITE (12,2010) IFIX(raaBlockF(iInt,1)),raaBlockF(iInt,2),raaBlockF(iInt,3),iTag
          WRITE (*,*) (raaBlockF(iInt,1)),raaBlockF(iInt,2),raaBlockF(iInt,3),iTag
          END DO
        END IF

      IF (iCntG .GT. 0) THEN
        DO iInt=1,iCntG
          WRITE (11,2000) iInt,IFIX(raaArrG(iInt,1)),raaArrG(iInt,2),raaArrG(iInt,3),iTag
          END DO
        CALL doBlock(raaArrG,iCntG,raaBlockG,iCnt2G)
        DO iInt=1,iCnt2G
          WRITE (12,2010) IFIX(raaBlockG(iInt,1)),raaBlockG(iInt,2),raaBlockG(iInt,3),iTag
          WRITE (*,*) (raaBlockG(iInt,1)),raaBlockG(iInt,2),raaBlockG(iInt,3),iTag
          END DO
        END IF

      IF (iCntH .GT. 0) THEN
        DO iInt=1,iCntH
          WRITE (11,2000) iInt,IFIX(raaArrH(iInt,1)),raaArrH(iInt,2),raaArrH(iInt,3),iTag
          END DO
        CALL doBlock(raaArrH,iCntH,raaBlockH,iCnt2H)
        DO iInt=1,iCnt2H
          WRITE (12,2010) IFIX(raaBlockH(iInt,1)),raaBlockH(iInt,2),raaBlockH(iInt,3),iTag
          WRITE (*,*) (raaBlockH(iInt,1)),raaBlockH(iInt,2),raaBlockH(iInt,3),iTag
          END DO
        END IF

      IF (iCntJ .GT. 0) THEN
        DO iInt=1,iCntJ
          WRITE (11,2000) iInt,IFIX(raaArrJ(iInt,1)),raaArrJ(iInt,2),raaArrJ(iInt,3),iTag
          END DO
        CALL doBlock(raaArrJ,iCntJ,raaBlockJ,iCnt2J)
        DO iInt=1,iCnt2J
          WRITE (12,2010) IFIX(raaBlockJ(iInt,1)),raaBlockJ(iInt,2),raaBlockJ(iInt,3),iTag
          WRITE (*,*) (raaBlockJ(iInt,1)),raaBlockJ(iInt,2),raaBlockJ(iInt,3),iTag
          END DO
        END IF

      IF (iCntK .GT. 0) THEN
        DO iInt=1,iCntK
          WRITE (11,2000) iInt,IFIX(raaArrK(iInt,1)),raaArrK(iInt,2),raaArrK(iInt,3),iTag
          END DO
        CALL doBlock(raaArrK,iCntK,raaBlockK,iCnt2K)
        DO iInt=1,iCnt2K
          WRITE (12,2010) IFIX(raaBlockK(iInt,1)),raaBlockK(iInt,2),raaBlockK(iInt,3),iTag
          END DO
        END IF

      IF (iCntP .GT. 0) THEN
        DO iInt=1,iCntP
          WRITE (11,2000) iInt,IFIX(raaArrP(iInt,1)),raaArrP(iInt,2),raaArrP(iInt,3),iTag
          END DO
        CALL doBlock(raaArrP,iCntP,raaBlockP,iCnt2P)
        DO iInt=1,iCnt2P
          WRITE (12,2010) IFIX(raaBlockP(iInt,1)),raaBlockP(iInt,2),raaBlockP(iInt,3),iTag
          END DO
        END IF

      IF (iCntQ .GT. 0) THEN
        DO iInt=1,iCntQ
          WRITE (11,2000) iInt,IFIX(raaArrQ(iInt,1)),raaArrQ(iInt,2),raaArrQ(iInt,3),iTag
          END DO
        CALL doBlock(raaArrQ,iCntQ,raaBlockQ,iCnt2Q)
        DO iInt=1,iCnt2Q
          WRITE (12,2010) IFIX(raaBlockQ(iInt,1)),raaBlockQ(iInt,2),raaBlockQ(iInt,3),iTag
          END DO
        END IF

      IF (iCntR .GT. 0) THEN
        DO iInt=1,iCntR
          WRITE (11,2000) iInt,IFIX(raaArrR(iInt,1)),raaArrR(iInt,2),raaArrR(iInt,3),iTag
          END DO
        CALL doBlock(raaArrR,iCntR,raaBlockR,iCnt2R)
        DO iInt=1,iCnt2R
          WRITE (12,2010) IFIX(raaBlockR(iInt,1)),raaBlockR(iInt,2),raaBlockR(iInt,3),iTag
          END DO
        END IF

      IF (iCntS .GT. 0) THEN
        DO iInt=1,iCntS
          WRITE (11,2000) iInt,IFIX(raaArrS(iInt,1)),raaArrS(iInt,2),raaArrS(iInt,3),iTag
          END DO
        CALL doBlock(raaArrS,iCntS,raaBlockS,iCnt2S)
        DO iInt=1,iCnt2S
          WRITE (12,2010) IFIX(raaBlockS(iInt,1)),raaBlockS(iInt,2),raaBlockS(iInt,3),iTag
          END DO
        END IF

      IF (iCntM .GT. 0) THEN
        DO iInt=1,iCntM
          WRITE (11,2000) iInt,IFIX(raaArrM(iInt,1)),raaArrM(iInt,2),raaArrM(iInt,3),iTag
          END DO
        CALL doBlock(raaArrM,iCntM,raaBlockM,iCnt2M)
        DO iInt=1,iCnt2M
          WRITE (12,2010) IFIX(raaBlockM(iInt,1)),raaBlockM(iInt,2),raaBlockM(iInt,3),iTag
          END DO
        END IF

      IF (iCntN .GT. 0) THEN
        DO iInt=1,iCntN
          WRITE (11,2000) iInt,IFIX(raaArrN(iInt,1)),raaArrN(iInt,2),raaArrN(iInt,3),iTag
          END DO
        CALL doBlock(raaArrN,iCntN,raaBlockN,iCnt2N)
        DO iInt=1,iCnt2N
          WRITE (12,2010) IFIX(raaBlockN(iInt,1)),raaBlockN(iInt,2),raaBlockN(iInt,3),iTag
          END DO
        END IF

      IF (iCntO .GT. 0) THEN
        DO iInt=1,iCntO
          WRITE (11,2000) iInt,IFIX(raaArrO(iInt,1)),raaArrO(iInt,2),raaArrO(iInt,3),iTag
          END DO
        CALL doBlock(raaArrO,iCntO,raaBlockO,iCnt2O)
        DO iInt=1,iCnt2O
          WRITE (12,2010) IFIX(raaBlockO(iInt,1)),raaBlockO(iInt,2),raaBlockO(iInt,3),iTag
          END DO
        END IF

      IF (iCntV .GT. 0) THEN
        DO iInt=1,iCntV
          WRITE (11,2000) iInt,IFIX(raaArrV(iInt,1)),raaArrV(iInt,2),raaArrV(iInt,3),iTag
          END DO
        CALL doBlock(raaArrV,iCntV,raaBlockV,iCnt2V)
        DO iInt=1,iCnt2V
          WRITE (12,2010) IFIX(raaBlockV(iInt,1)),raaBlockV(iInt,2),raaBlockV(iInt,3),iTag
          END DO
        END IF

      IF (iCntU .GT. 0) THEN
        DO iInt=1,iCntU
          WRITE (11,2000) iInt,IFIX(raaArrU(iInt,1)),raaArrU(iInt,2),raaArrU(iInt,3),iTag
          END DO
        CALL doBlock(raaArrU,iCntU,raaBlockU,iCnt2U)
        DO iInt=1,iCnt2U
          WRITE (12,2010) IFIX(raaBlockU(iInt,1)),raaBlockU(iInt,2),raaBlockU(iInt,3),iTag
          END DO
        END IF

      CLOSE(11)
      CLOSE(12)

 2000 FORMAT(I3,'  ',I3,'  ',f14.6,'   ',f14.6,'     ',I3)
 2010 FORMAT(I3,'  ',f14.6,'   ',f14.6,'     ',I3)


      END


!************************************************************************
!************************************************************************
!************************************************************************
! just add on to this subroutine as needed
      SUBROUTINE UPDATE_DoInit( &
                 iCntF,raaArrF,iCntG,raaArrG,iCntH,raaArrH, &
                 iCntJ,raaArrJ,iCntK,raaArrK,iCntP,raaArrP, &
                 iCntM,raaArrM,iCntN,raaArrN,iCntO,raaArrO, &
                 iCntQ,raaArrQ,iCntR,raaArrR,iCntS,raaArrS, &
                 iCntV,raaArrV,iCntU,raaArrU)

      include '../INCLUDE/TempF90/scatterparam.f90'
      REAL raaArrF(kMaxPts,3),raaBlockF(kMaxPts,3)
      REAL raaArrG(kMaxPts,3),raaBlockG(kMaxPts,3)
      REAL raaArrH(kMaxPts,3),raaBlockH(kMaxPts,3)
      REAL raaArrK(kMaxPts,3),raaBlockK(kMaxPts,3)
      REAL raaArrP(kMaxPts,3),raaBlockP(kMaxPts,3)
      REAL raaArrQ(kMaxPts,3),raaBlockQ(kMaxPts,3)
      REAL raaArrR(kMaxPts,3),raaBlockR(kMaxPts,3)
      REAL raaArrS(kMaxPts,3),raaBlockS(kMaxPts,3)
      REAL raaArrM(kMaxPts,3),raaBlockM(kMaxPts,3)
      REAL raaArrN(kMaxPts,3),raaBlockN(kMaxPts,3)
      REAL raaArrO(kMaxPts,3),raaBlockO(kMaxPts,3)
      REAL raaArrJ(kMaxPts,3),raaBlockJ(kMaxPts,3)
      REAL raaArrV(kMaxPts,3),raaBlockV(kMaxPts,3)
      REAL raaArrU(kMaxPts,3),raaBlockU(kMaxPts,3)

      INTEGER iCntJ, iCntK, iCntM, iCntN, iCntP, iCntQ, iCntR, iCntS, iCntO, iCntV, iCntU
      INTEGER iCntF, iCntG, iCntH

      INTEGER iInt

      iCntF = 0
      iCntG = 0
      iCntH = 0
      iCntJ = 0
      iCntK = 0
      iCntM = 0
      iCntN = 0
      iCntO = 0
      iCntP = 0
      iCntQ = 0
      iCntR = 0
      iCntS = 0
      iCntV = 0
      iCntU = 0

      DO iInt = 1,kMaxPts
        raaArrF(iInt,1) = 10000000.0
        raaArrF(iInt,2) = 10000000.0 
        raaArrF(iInt,3) = 10000000.0      

        raaArrG(iInt,1) = 10000000.0
        raaArrG(iInt,2) = 10000000.0 
        raaArrG(iInt,3) = 10000000.0      

        raaArrH(iInt,1) = 10000000.0
        raaArrH(iInt,2) = 10000000.0 
        raaArrH(iInt,3) = 10000000.0      

        raaArrJ(iInt,1) = 10000000.0
        raaArrJ(iInt,2) = 10000000.0 
        raaArrJ(iInt,3) = 10000000.0      

        raaArrK(iInt,1) = 10000000.0
        raaArrK(iInt,2) = 10000000.0 
        raaArrK(iInt,3) = 10000000.0      

        raaArrP(iInt,1) = 10000000.0
        raaArrP(iInt,2) = 10000000.0 
        raaArrP(iInt,3) = 10000000.0      

        raaArrM(iInt,1) = 10000000.0
        raaArrM(iInt,2) = 10000000.0 
        raaArrM(iInt,3) = 10000000.0      

        raaArrN(iInt,1) = 10000000.0
        raaArrN(iInt,2) = 10000000.0 
        raaArrN(iInt,3) = 10000000.0      

        raaArrO(iInt,1) = 10000000.0
        raaArrO(iInt,2) = 10000000.0 
        raaArrO(iInt,3) = 10000000.0      

        raaArrQ(iInt,1) = 10000000.0
        raaArrQ(iInt,2) = 10000000.0 
        raaArrQ(iInt,3) = 10000000.0      

        raaArrR(iInt,1) = 10000000.0
        raaArrR(iInt,2) = 10000000.0 
        raaArrR(iInt,3) = 10000000.0      

        raaArrS(iInt,1) = 10000000.0
        raaArrS(iInt,2) = 10000000.0 
        raaArrS(iInt,3) = 10000000.0      

        raaArrV(iInt,1) = 10000000.0
        raaArrV(iInt,2) = 10000000.0 
        raaArrV(iInt,3) = 10000000.0      

        raaArrU(iInt,1) = 10000000.0
        raaArrU(iInt,2) = 10000000.0 
        raaArrU(iInt,3) = 10000000.0      
        END DO

      RETURN
      END

!************************************************************************
!************************************************************************
!************************************************************************

! this subroutine goes thru raaArr and puts every gas in terms of freq blocks
      SUBROUTINE doBlock(raaArr,iCnt,raaBlock,iCnt2)

! the information in raaArr is summarized and saved to raaBlock
! iCnt  is the total number of r*.dat files read in file comdatabase
! iCnt2 is the summarized number of blocks that are output

      IMPLICIT NONE
      include '../INCLUDE/TempF90/scatterparam.f90'

      REAL raaArr(kMaxPts,3),raaBlock(kMaxPts,3)
      INTEGER iCnt,iCnt2

! local variables
      INTEGER iT,iBlock,iBoo
      REAL rGasID,rGasIDold,rStart,rStop

      rGasIDold=raaArr(1,1)
      rStart=raaArr(1,2)
      rStop=raaArr(1,3) 

      iT=2 
      iBlock=1
             
 70   CONTINUE
      IF (iT .LE. iCnt) THEN
        rGasId=raaArr(iT,1)
! check to see of the gas is the same
        IF (abs(rGasId-rGasIdOld) .LE. 0.1) THEN
! check if the stop freq set frm previous index, and this start, are same
          IF (abs(rStop-raaArr(iT,2)) .LE. 0.1) THEN
            rStop=raaArr(iT,3)
          ELSE
            raaBlock(iBlock,1)=rGasIDold
            raaBlock(iBlock,2)=rStart 
            raaBlock(iBlock,3)=rStop
            iBlock=iBlock+1
            rGasIDold=raaArr(iT,1)
            rStart=raaArr(iT,2)
            rStop=raaArr(iT,3) 
            END IF
        ELSE     
          raaBlock(iBlock,1)=rGasIDold
          raaBlock(iBlock,2)=rStart 
          raaBlock(iBlock,3)=rStop
          iBlock=iBlock+1
          rGasIDold=raaArr(iT,1)
          rStart=raaArr(iT,2)
          rStop=raaArr(iT,3) 
          END IF
        iT=iT+1
        GO TO 70
        END IF
 
      raaBlock(iBlock,1)=rGasIDold
      raaBlock(iBlock,2)=rStart 
      raaBlock(iBlock,3)=rStop
 
      iCnt2=iBlock

      DO iT=1,iCnt2
        print *,iT,raaBlock(iT,1),raaBlock(iT,2),raaBlock(iT,3) 
        END DO

      RETURN
      END 
 
!************************************************************************
! this subroutine takes the string apart, finding the gas id and freq range
! the string is of the form rFREQ_gID.dat eg r605_g1.dat
      SUBROUTINE process(caStr,raaArr,iCnt,rJump)

      IMPLICIT NONE
      include '../INCLUDE/TempF90/scatterparam.f90'

! take string caStr apart, finding the gas id and freq range, saving for each
! file, the gasID in raaArr(1), start freq in raaArr(2) and stop freq in 
! raaArr(3)
! the string is of the form rFREQ_gID.dat eg r605_g1.dat
      CHARACTER*80 caStr
      REAL raaArr(kMaxPts,3),rJump
      INTEGER iCnt

! local variables
      CHARACTER*10 caBlank
      INTEGER i1,i2,iStart,iGasID
      REAL rStart,rStop,rXStart

      iCnt = iCnt+1
!      print *,'The file prefix and iCnt are : ',caStr(1:1),iCnt
      
! first read the start freq ... stop freq is rJump cm away 
      CALL doblank(caBlank)
      i1 = 2
      i2 = i1
 66   CONTINUE
      IF (caStr(i2:i2) .NE. '_') THEN
        i2 = i2+1
        GO TO 66 
        END IF 

      i2 = i2-1
      caBlank(1:i2-i1+1) = caStr(i1:i2)
!      read(caBlank,*) iStart
!      rStart=iStart*1.0
      read(caBlank,*) rStart
      rStop          = rStart + rJump
      raaArr(iCnt,2) = rStart
      raaArr(iCnt,3) = rStop

! now read the gas id
      CALL doblank(caBlank)
      i1 = i2+3
      i2 = i1
 69   CONTINUE
      IF (caStr(i2:i2) .NE. '.')  THEN
        i2 = i2+1
        GO TO 69
        END IF
      i2 = i2-1
      caBlank(1:12-i1+1) = caStr(i1:i2)
      read(caBlank,*) iGasID
      raaArr(iCnt,1)     = iGasID*1.00

      RETURN
      END 

!************************************************************************
! this subroutine sorts the array in terms of index iIndex. 
! this is a bubble sort : refer Dale/Lilly : Data Structures
      SUBROUTINE dosort(raaArr,iCnt,iIndex)

      IMPLICIT NONE
      include '../INCLUDE/TempF90/scatterparam.f90'

! sort rows 1..iCnt of array raaArr, according to the elements in column iCnt
      REAL raaArr(kMaxPts,3)
      INTEGER iCnt,iIndex

      REAL rT(3)
      INTEGER iI,iJ,iTrue,i1,i2

      i1=1
      iTrue=1

 60   CONTINUE
      IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
        i2=iCnt
        iTrue=-1
 70     CONTINUE
        IF (i2 .GT. i1) THEN
          IF (raaArr(i2,iIndex) .LT. raaArr(i2-1,iIndex)) THEN
            iTrue=1
            rT(1)=raaArr(i2-1,1)
            rT(2)=raaArr(i2-1,2)
            rT(3)=raaArr(i2-1,3)
	
            raaArr(i2-1,1)=raaArr(i2,1)
            raaArr(i2-1,2)=raaArr(i2,2)
            raaArr(i2-1,3)=raaArr(i2,3)

            raaArr(i2,1)=rT(1)
            raaArr(i2,2)=rT(2)
            raaArr(i2,3)=rT(3)

            END IF
          i2=i2-1
          GO TO 70
          END IF
      
        i1=i1+1
        GO TO 60
        END IF


       RETURN
       END
!************************************************************************
! this subroutine BLANKS the string
      SUBROUTINE doblank(caBlank)

      IMPLICIT NONE
      include '../INCLUDE/TempF90/scatterparam.f90'

      CHARACTER*10 caBlank
      INTEGER iInt

      DO iInt=1,10
        caBlank(iInt:iInt)=' '
        END DO

      RETURN
      END 
!***********************************************************************

