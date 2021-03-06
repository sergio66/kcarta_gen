c Copyright 2007
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************************************************************************
c************************************************************************
c this subroutine figures out the Rayleigh optical depths
c it also redoes the atmosphere profiles so we have a 3 layer atmosphere, 
c intead of a 100 layer atmosphere : between 35 - 100 km, 
c                                    between 20 - 35 km (ozone in stratosphere)
c                                    between 0 - 20 km  (troposphere)

      SUBROUTINE AddRayleigh(
     $               raFreq,iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               raVTemp,raPressLevels,raThickness,raSunAngles,iNpmix,
     $               iProfileLayers,ICLDTOPKCARTA,ICLDBOTKCARTA,
     $               iaCldTop,iaCldbot,iNclouds,
     $               raMeanT,iSquish,iaMeanLay,
     $               raaExtTemp,raaScatTemp,raaAsymTemp)

      IMPLICIT NONE
 
      include '../INCLUDE/scatter.param'

c input vars
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds),iNclouds
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iNpmix
      INTEGER iProfileLayers,ICLDTOPKCARTA,ICLDBOTKCARTA
      REAL rFracTop,rFracBot,raSunAngles(kProfLayer)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
c input/output
      REAL raaExtTemp(kMaxPts,kMixFilRows),raaScatTemp(kMaxPts,kMixFilRows)
      REAL raaAsymTemp(kMaxPts,kMixFilRows)
c output
      REAL raMeanT(kProfLayer)
      INTEGER iSquish,iaMeanLay(kProfLayer),iS1,iS2,iX

c local
      REAL raaRayleigh(kMaxPts,kProfLayer),rX,rP1,rP2
      REAL raVT1(kMixFilRows),interptemp
      INTEGER iaRadLayer(kProfLayer),iFr,iLay,iL,iP1,iP2,iM,iaC(3)
      INTEGER iNumCldLay
      INTEGER IACLDTOPX(kMaxClouds), IACLDBOTX(kMaxClouds)
      REAL raPZ(kProfLayer),raTZ(kProfLayer)
      REAL raaE1(kMaxPts,kStdkCartaKK),raaE2(kMaxPts,kStdkCartaKK)
      REAL raPZ1(kStdkCartaKK),raTZ1(kStdkCartaKK),raSum(kStdkCartaKK)
      INTEGER iaLayers(kProfLayer)

      DO iL = 1,iNclouds
        iaCldTopX(iL) = kProfLayer - (iaCldTop(iL)+1) + 1
        iaCldBotX(iL) = kProfLayer - iaCldBot(iL) + 1
      END DO
      IF (iNclouds .LT. 2) THEN
        iL = 2       
        iaCldTopX(iL) = -1
        iaCldBotX(iL) = -1
      END IF

c       print *,iaCldTop(1),iaCldbot(1),iaCldTopX(1),iaCldbotX(1)
c       print *,iaCldTop(2),iaCldbot(2),iaCldTopX(2),iaCldbotX(2)
c       print *,ICLDTOPKCARTA,ICLDBOTKCARTA,iNClouds

      IF (iNclouds .EQ. 0) THEN
         write (KstdErr,*) 'Huh, iNclouds = 0 ???'
         CALL DoStop
      ELSEIF (iNclouds .EQ. 1) THEN
        iNumCldLay = ICLDTOPKCARTA - ICLDBOTKCARTA + 1
      ELSEIF (iNclouds .EQ. 2) THEN
        iNumCldLay = iaCldBotX(1)-iaCldTopX(1) + 1 + 
     $               iaCldBotX(2)-iaCldTopX(2) + 1
      ELSEIF (iNclouds .GT. 2) THEN
        write (KstdErr,*) 'Huh, iNclouds > 2 ???'
        CALL DoStop
      END IF

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
      write(kStdErr,*) 'Rayleigh routine HAS A BUG in trmp of lowest layer'
      CALL DOSTOP

c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      IF ((kSolar .GE. 0) .AND. (raFreq(1) .GE. 5000.0) .AND. 
     $     (raSunAngles(iaRadLayer(1)) .LE. 90) .AND. 
     $     (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
        write(kStdWarn,*) 'adding on Rayleigh ...'
        CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $                    raThickness,raPZ,raTZ,raaRayleigh) 
c        DO iM = 1,kProfLayer
c          print *,'Rayleigh ',iM,raPZ(iM),raTZ(iM),raaRayleigh(1,iM)
c        END DO
      END IF

      !! now squish everything down to 3 layers
      !! 0 km ~ gnd, 20 km ~ 60 mb, 35 km  ~ 6 mb, TOA ~ 0 mb
      IF (iaRadLayer(1) .LT. iaRadLayer(2)) THEN      
        iaLayers(1) = iaRadLayer(1)         !!downlook instr  GND
        iaLayers(4) = iaRadLayer(iNumlayer) !!downlook instr  TOA
      ELSE
        iaLayers(1) = iaRadLayer(iNumlayer) !!uplook instr    GND
        iaLayers(4) = iaRadLayer(1)         !!uplook instr    TOA
      END IF
      iX = iaLayers(1)
 10   CONTINUE
      IF (raPZ(iX) .GE. 60.0) THEN
        iX = iX + 1
        GOTO 10
      END IF
      iaLayers(2) = iX

      iX = iaLayers(1)
 20   CONTINUE
      IF (raPZ(iX) .GE. 6.0) THEN
        iX = iX + 1
        GOTO 20
      END IF
      iaLayers(3) = iX

      iSquish = 3
      DO iL = 1,iSquish+1
        raPZ1(iL) = 0.0
        raTZ1(iL) = 0.0
        raSum(iL) = 0.0
        DO iFr = 1,kMaxpts
          raaE1(iFr,iL) = 0.0
          raaE2(iFr,iL) = 0.0
        END DO
      END DO

      DO iM = 1,iSquish
        iS1 = iaLayers(iM)
        IF (iSquish .LE. 2) THEN
          iS2 = iaLayers(iM+1) - 1
        ELSE
          iS2 = iaLayers(iM+1)
        END IF
        DO iL = iS1,iS2
          raSum(iM) = raSum(iM) + raPZ(iL)
          raTZ1(iM) = raTZ1(iM) + raPZ(iL)*raTZ(iL)
          DO iFr = 1,kMaxpts
            raaE1(iFr,iM) = raaE1(iFr,iM) + raaExtTemp(iFr,iL)
            raaE2(iFr,iM) = raaE2(iFr,iM) + raaRayleigh(iFr,iL)
          END DO
        END DO
      END DO

      DO iM = 1,iSquish
        raTZ1(iM) = raTZ1(iM)/raSum(iM)
c        print *,iM,iaLayers(iM)
      END DO

      !! now stuff all this into raaExtTemp and so on
      !! also reset the parameters for TWOSTREAM

      iNclouds = 1

c      print *,iaRadLayer(1),iaRadLayer(2)
c      print *,ICLDTOPKCARTA,ICLDBOTKCARTA,iaCldTop(1),iaCldBot(1)

      IF (iaRadLayer(1) .LT. iaRadLayer(2)) THEN
        !! downlook instr
        ICLDTOPKCARTA = iaRadLayer(1)+(iSquish-1)
        ICLDBOTKCARTA = iaRadLayer(1)
        iaCldTop(1) = kProfLayer - (iaRadLayer(1)+iSquish)+1
        iaCldBot(1) = kProfLayer - (iaRadLayer(1))+1

        DO iLay = 1,iSquish
          iL = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
             raaExtTemp(iFr,iL)  = raaE1(iFr,iLay) + raaE2(iFr,iLay)
             raaScatTemp(iFr,iL) = raaE2(iFr,iLay)
             raaAsymTemp(iFr,iL) = 0.0
           END DO
           raVTemp(iL) = raTZ1(iLay)
c           print *,'sqh ',10000/raFreq(1),iLay,iL,
c     $              raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $              exp(-raaExtTemp(1,iL)),raVTemp(iL)
         END DO
        DO iLay = iSquish+1,iNumLayer
          iL = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
             raaExtTemp(iFr,iL)  = 1.0e-10
             raaScatTemp(iFr,iL) = 0.0
             raaAsymTemp(iFr,iL) = 0.0
           END DO
         END DO

      ELSE
        !! uplook instr
        ICLDTOPKCARTA = iaRadLayer(iNumLayer)+(iSquish-1)
        ICLDBOTKCARTA = iaRadLayer(iNumLayer)
        iaCldTop(1) = kProfLayer - (iaRadLayer(iNumLayer)+iSquish)+1
        iaCldBot(1) = kProfLayer - (iaRadLayer(iNumLayer))+1

        DO iLay = 1,iSquish
          iL = iaRadLayer(iNumLayer-iLay+1)
          DO iFr = 1,kMaxPts
             raaExtTemp(iFr,iL)  = raaE1(iFr,iLay) + raaE2(iFr,iLay)
             raaScatTemp(iFr,iL) = raaE2(iFr,iLay)
             raaAsymTemp(iFr,iL) = 0.0
           END DO
           raVTemp(iL) = raTZ1(iLay)
c           print *,'sqh ',10000/raFreq(1),iLay,iL,
c     $              raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $              exp(-raaExtTemp(1,iL)),raVTemp(iL)
         END DO
        DO iLay = 1,iNumLayer-iSquish
          iL = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
             raaExtTemp(iFr,iL)  = 1.0e-10
             raaScatTemp(iFr,iL) = 0.0
             raaAsymTemp(iFr,iL) = 0.0
           END DO
c           print *,'sqh ',10000/raFreq(1),iLay,iL,
c     $              raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $              exp(-raaExtTemp(1,iL)),raVTemp(iL)
         END DO

      END IF
c      print *,ICLDTOPKCARTA,ICLDBOTKCARTA,iaCldTop(1),iaCldBot(1)

      RETURN
      END

c************************************************************************
      SUBROUTINE AddRayleighWorks(
     $               raFreq,iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               raVTemp,raPressLevels,raThickness,raSunAngles,iNpmix,
     $               iProfileLayers,ICLDTOPKCARTA,ICLDBOTKCARTA,
     $               iaCldTop,iaCldbot,iNclouds,
     $               raMeanT,iSquish,iaMeanLay,
     $               raaExtTemp,raaScatTemp,raaAsymTemp)

      IMPLICIT NONE
 
      include '../INCLUDE/scatter.param'

c input vars
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds),iNclouds
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iNpmix
      INTEGER iProfileLayers,ICLDTOPKCARTA,ICLDBOTKCARTA
      REAL rFracTop,rFracBot,raSunAngles(kProfLayer)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
c input/output
      REAL raaExtTemp(kMaxPts,kMixFilRows),raaScatTemp(kMaxPts,kMixFilRows)
      REAL raaAsymTemp(kMaxPts,kMixFilRows)
c output
      REAL raMeanT(kProfLayer)
      INTEGER iSquish,iaMeanLay(kProfLayer)

c local
      REAL raaRayleigh(kMaxPts,kProfLayer),rX,rP1,rP2
      REAL raVT1(kMixFilRows),interptemp,raTemp(kMaxPts)
      INTEGER iaRadLayer(kProfLayer),iFr,iLay,iL,iP1,iP2,iM,iaC(3)
      INTEGER iNumCldLay
      INTEGER IACLDTOPX(kMaxClouds), IACLDBOTX(kMaxClouds)
      REAL raPZ(kProfLayer),raTZ(kProfLayer)

       DO iL = 1,iNclouds
         iaCldTopX(iL) = kProfLayer - (iaCldTop(iL)+1) + 1
         iaCldBotX(iL) = kProfLayer - iaCldBot(iL) + 1
       END DO
       IF (iNclouds .LT. 2) THEN
         iL = 2       
         iaCldTopX(iL) = -1
         iaCldBotX(iL) = -1
       END IF

c       print *,iaCldTop(1),iaCldbot(1),iaCldTopX(1),iaCldbotX(1)
c       print *,iaCldTop(2),iaCldbot(2),iaCldTopX(2),iaCldbotX(2)
c       print *,ICLDTOPKCARTA,ICLDBOTKCARTA,iNClouds

      IF (iNclouds .EQ. 0) THEN
         write (KstdErr,*) 'Huh, iNclouds = 0 ???'
         CALL DoStop
      ELSEIF (iNclouds .EQ. 1) THEN
        iNumCldLay = ICLDTOPKCARTA - ICLDBOTKCARTA + 1
      ELSEIF (iNclouds .EQ. 2) THEN
        iNumCldLay = iaCldBotX(1)-iaCldTopX(1) + 1 + 
     $               iaCldBotX(2)-iaCldTopX(2) + 1
      ELSEIF (iNclouds .GT. 2) THEN
        write (KstdErr,*) 'Huh, iNclouds > 2 ???'
        CALL DoStop
      END IF

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      IF ((kSolar .GE. 0) .AND. (raFreq(1) .GE. 5000.0) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .LE. 90) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
        write(kStdWarn,*) 'adding on Rayleigh ...'
        CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $                    raThickness,raPZ,raTZ,raaRayleigh) 
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 0.0
        END DO

        !!! no need to delta scale the optical depths, w, g as
        !!! for Rayleigh, assume w = 1, g = 0 ==> w' = 1, g' = 0
        DO iLay = 1,iNumLayer
          iL   = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + raaRayleigh(iFr,iL) 
          END DO
        END DO
        !! divvy equally amongst the layers from ICLDBOTKCARTA to ICLDTOPKCARTA
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raTemp(iFr)/iNumCldLay
        END DO

        !!!!!!!!! -------------------------------------------- !!!!!!!
        IF (iNclouds .GE. 1) THEN
          IF (iaRadLayer(1) .LT. iaRadLayer(2)) THEN
            write(kStdWarn,*) 'adding on Rayleigh for downlook instr, c1'
            !!! now add on the the optical depths
            DO iLay = 1,iNumLayer
              iL   = iaRadLayer(iLay)
              IF (iL .LE. iaCldTopX(1) .AND. iL .GE. iaCldbotX(1)) THEN
                DO iFr = 1,kMaxPts
                  raaExtTemp(iFr,iL)  = raaExtTemp(iFr,iL) + raTemp(iFr)
                  raaScatTemp(iFr,iL) = raTemp(iFr)
                  raaAsymTemp(iFr,iL) = 0.0
                END DO
c                print *,'ray1 ',10000/raFreq(1),iLay,iL,
c     $                 raTemp(1),raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $                 exp(-raaExtTemp(1,iL))
              END IF
            END DO

          ELSEIF (iaRadLayer(1) .GT. iaRadLayer(2)) THEN
            write(kStdWarn,*) 'adding on Rayleigh for uplook instr, c1'
            !!! now add on the the optical depths
            DO iLay = 1,iNumLayer    !! was iNumLayer,iNumLayer
              iL   = iaRadLayer(iLay)
              IF (iL .LE. iaCldTopX(1) .AND. iL .GE. iaCldbotX(1)) THEN
                DO iFr = 1,kMaxPts
                  raaExtTemp(iFr,iL)  = raaExtTemp(iFr,iL) + raTemp(iFr)
                  raaScatTemp(iFr,iL) = raTemp(iFr)
                  raaAsymTemp(iFr,iL) = 0.0
                END DO
c                print *,'ray ',10000/raFreq(1),iLay,iL,
c     $                  raTemp(1),raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $                  exp(-raaExtTemp(1,iL))
              END IF
            END DO
          END IF
        END IF

        !!!!!!!!! -------------------------------------------- !!!!!!!
        IF (iNclouds .EQ. 2) THEN
          IF (iaRadLayer(1) .LT. iaRadLayer(2)) THEN
            write(kStdWarn,*) 'adding on Rayleigh for downlook instr, c2'
            !!! now add on the the optical depths
            DO iLay = 1,iNumLayer
              iL   = iaRadLayer(iLay)
              IF (iL .LE. iaCldTopX(2) .AND. iL .GE. iaCldbotX(2)) THEN
                DO iFr = 1,kMaxPts
                  raaExtTemp(iFr,iL)  = raaExtTemp(iFr,iL) + raTemp(iFr)
                  raaScatTemp(iFr,iL) = raTemp(iFr)
                  raaAsymTemp(iFr,iL) = 0.0
                END DO
c                print *,'ray2 ',10000/raFreq(1),iLay,iL,
c     $                 raTemp(1),raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $                 exp(-raaExtTemp(1,iL))
              END IF
            END DO

          ELSEIF (iaRadLayer(1) .GT. iaRadLayer(2)) THEN
            write(kStdWarn,*) 'adding on Rayleigh for uplook instr, c2'
            !!! now add on the the optical depths
            DO iLay = 1,iNumLayer    !! was iNumLayer,iNumLayer
              iL   = iaRadLayer(iLay)
              IF (iL .LE. iaCldTopX(2) .AND. iL .GE. iaCldbotX(2)) THEN
                DO iFr = 1,kMaxPts
                  raaExtTemp(iFr,iL)  = raaExtTemp(iFr,iL) + raTemp(iFr)
                  raaScatTemp(iFr,iL) = raTemp(iFr)
                  raaAsymTemp(iFr,iL) = 0.0
                END DO
c                print *,'ray ',10000/raFreq(1),iLay,iL,
c     $                  raTemp(1),raaExtTemp(1,iL),raaScatTemp(1,iL),
c     $                  exp(-raaExtTemp(1,iL))
              END IF
            END DO
          END IF
        END IF
        !!!!!!!!! -------------------------------------------- !!!!!!!
      END IF

      RETURN
      END

c************************************************************************
