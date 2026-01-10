c Copyright 2006 
c University of Maryland Baltimore County 
c All Rights Reserved

c (1)JacobGasAmtFM1,JacobTempFM1 : jacobians from the forward model
c    (includes solar contribution and thermal diffusive contribution)
c (2)Surface Reflectivity = 1/pi for thermal
c    Surface Reflectance for solar is defined by user

c the following variables are not size kMaxPtsJac or kProfLayerJac as they
c are well defined in the other non Jacobian routines
c raFreq(kMaxPts),raUseEmissivity(kMaxPts),raVTemp(kMixFilRows),
c iaaRadLayer(kMaxAtm,kProfLayer),raaExtTemp(kMaxPts,kMixFilRows)
c     $              raSurface,raSun,raThermal,raInten,
c     $              raSunRefl,
c     $              raLayAngles,raSunAngles)

c we also allow the user to compute the temperature jacobians in 
c one of three ways
c recall r(v)= sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

c************************************************************************
c****** THESE ARE THE SCATTERING JACOBIANS FOR THE UP LOOK INSTR ******
c************************************************************************
c this subroutine does the Jacobians for downward looking instrument
      SUBROUTINE UpwardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,rDelta,
     $               iNpMix,iTag,iActualTag,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,
     $            iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rDelta is the wavenumber step size
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which 25 cm-1 block being output
c iNumGases is the number of gases to include
c iaGases is the integer array of GasID's
c iNumLayer is the number of layers in the atmosphere # iAtm
c iaaRadLayer is the list of radiating mixed paths to be included
c raVTemp are the layer temperatures

c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
c raaAllDT has the cumulative d/dT coeffs for current freq block
C     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
C        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
C        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

c raaSumAbCoeff is the cumulative absorption coeffs
c raaAmt  has the gas profiles
c raInten has the radiance vector
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi or kSolarRefl
c raLayAngles = layer dependent satellite view angles
c raSunAngles = layer dependent sun view angles
c raaMix is the mixing table
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rDelta,raPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts),rFracTop,rFracBot
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $     raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      CHARACTER*120 caJacobFile
      INTEGER iJacob,iaJacob(kMaxDQ),iProfileLayers,iTag,iActualTag
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
c this is for NLTE weight fcns
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c this is for the scattering parts
      REAL raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
      REAL raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
      REAL raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
      REAL raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
      REAL raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
      REAL raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      INTEGER iNpMix,iCLoudySky

c local variables
      INTEGER iNumGasesTemp,iaGasesTemp(kMaxGas)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneral(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      REAL radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
      REAL radSolardT(kMaxPtsJac),rWeight
      INTEGER iG,iLay,iIOUN,iLowest
      INTEGER DoGasJacob,iGasJacList
      INTEGER WhichGasPosn,iGasPosn
c for cloud stuff!
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1) 
      INTEGER iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot 
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iCloudLayerTop,iCloudLayerBot
      INTEGER iiDiv,iaRadLayer(kProfLayer),iIWPorDME,iL
      INTEGER iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
      REAL r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck
      REAL muSun,muSat,rSunAngle
      INTEGER iSolarRadOrJac,MP2Lay
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer)
      REAL raaSolarScatterCumul(kMaxPts,kProfLayer),raLMm1_toGnd(kMaxPts)
      INTEGER iWhichLayer

      INTEGER iDefault,iWhichJac,iFr
      INTEGER iDoAdd,iErr

      iDefault  = -1     !! do all jacs (Q,T,W,surface)

      iWhichJac = -1     !! do all jacs (Q,T,W,surface)
      iWhichJac = +20    !! only Q jacs
      iWhichJac = +30    !! only T jacs
      iWhichJac = +40    !! only W jacs
      iWhichJac = +50    !! only S jacs

      !! this only uses T(z) contribution from gases in iaJacob{}
      iWhichJac = -2     !! do all jacs (Q,T,W,surface)
      iWhichJac = +32    !! only T jacs

      iWhichJac = kActualJacs

      IF (iDefault .NE. iWhichJac) THEN 
        print *,'iDefault,iWhichJac = ',iDefault,iWhichJac
        END IF 

c calculate cos(SatAngle)  
      muSat = cos(rSatAngle*kPi/180.0)  

c as we are never directly loooking at the sun, there is a geometry factor  
      rSunAngle = raSunAngles(iaaRadlayer(1,iAtm))
      rOmegaSun = kOmegaSun
      IF (iDoSolar .GE. 0) THEN  
        rSunTemp = kSunTemp  
        write(kStdWarn,*) 'upward looking instrument .. daytime'  
      ELSE IF (iDoSolar .LT. 0) THEN  
        rSunTemp = 0.0  
        write(kStdWarn,*)'upward looking instrument .. nitetime'  
        END IF  

      muSun = 1.0       !!!default   
      IF (iDoSolar .GE. 0) THEN  
        muSun = cos(rSunAngle*kPi/180.0)  
        END IF  

      iDoSolar = kSolar  
      IF (iDoSolar .EQ. 0) THEN  
        !use 5700K 
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K' 
        rSunTemp = kSunTemp  
        DO iFr = 1,kMaxPts 
c compute the Plank radiation from the sun  
          rPlanck       = exp(r2*raFreq(iFr)/rSunTemp)-1.0  
          raSunTOA(iFr) = r1*((raFreq(iFr))**3)/rPlanck  
          raSunTOA(iFr) = raSunTOA(iFr)*rOmegaSun*muSun
          END DO  
      ELSEIF (iDoSolar .EQ. 1) THEN  
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files' 
        !read in data from file 
        CALL ReadSolarData(raFreq,raSunTOA,iTag) 
        DO iFr = 1,kMaxPts 
          raSunTOA(iFr) = raSunTOA(iFr)*rOmegaSun*muSun
          END DO  
        END IF 
           
      DO iLay = 1,kProfLayer
        iaCldLayer(iLay) = 0
        iaCldLayerIWPDME(iLay) = 0
        END DO

      iNumGasesTemp = iNumGases
      DO iLay = 1,iNumGases
        iaGasesTemp(iLay) = iaGases(iLay)
        END DO

      iNumGasesTemp = iNumGasesTemp + 1
      iaGasesTemp(iNumGasesTemp) = 201

      iNumGasesTemp = iNumGasesTemp + 1
      iaGasesTemp(iNumGasesTemp) = 202

c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
        END IF 
      DO iLay=1,iNumLayer 
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay) 
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

cccccccccccccccccccc set these all important variables ****************  
      DO iLay = 1,kProfLayer
        iaCldLayer(iLay) = -1
        iaCldLayerIWPDME(iLay) = -1
        END DO

      iCloudLayerTop = -1  
      iCloudLayerBot = -1  
      IF (iaRadLayer(1) .LT. kProfLayer) THEN  
        iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1  
        iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1  
        iiDiv = 0  
      ELSE  
        !!essentially do mod(iaRadLayer(1),kProfLayer)  
        iiDiv = 1            
 1010      CONTINUE  
        IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN  
          iiDiv = iiDiv + 1  
          GOTO 1010  
          END IF  
        iiDiv = iiDiv - 1  
        iLay = iiDiv  
        iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)  
        iLocalCldTop = iiDiv - iCldTopkCarta + 1  
        iLocalCldBot = iiDiv - iCldBotkCarta + 1  
        iiDiv = iLay  
        END IF  
 
      iOffSet = (kProfLayer-iNumLayer)
c      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1 
      DO iLay = iCldBotkCarta,iCldTopkCarta
        !!! this is for d/dq, d/dT
        iaCldLayer(kProfLayer-iLay+1) = 1 
        !!!this is for d/d(DME), d/d(IWP)
        iaCldLayerIWPDME(kProfLayer-iLay+1 - (kProfLayer-iNumLayer)) = 1
        END DO 
 
      IF (kSolar .GE. 0) THEN
        iSolarRadOrJac = +1
        CALL SolarScatterIntensity_Uplook(
     $      iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer,
     $      iNumLayer,iaRadLayer,
     $      raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot,
     $      iTag,iSolarRadOrJac,raaSolarScatter1Lay)
        END IF

c      DO iLay = 1,iNumLayer
c        iL = iaRadlayer(iLay)
c        print *,'<<<<<>>>>',iLay,iL,iaCldLayer(iLay),iaCldLayerIWPDME(iLay),
c     $            raSunAngles(iL),raLayAngles(iL),
c     $            raaExtTemp(1,iL),raaSSAlbTemp(1,iL),raaAsymTemp(1,iL),
c     $            raaSolarScatter1Lay(1,iL),
c     $            raaPhaseJacobASYM(1,iL),raaExtJacobIWP(1,iL),
c     $            raaExtJacobDME(1,iL)
c        END DO

cccccccccccccccccccc set these all important variables ****************  

      iIOUN = kStdJacob

      iLowest = iaaRadLayer(iAtm,1)            !!! this is for DOWN LOOk instr 
      iLowest = iaaRadLayer(iAtm,iNumLayer)    !!!modify for UPLOOK instr 
      iLowest = MOD(iLowest,kProfLayer) 

      IF (kJacobOutPut .GE. 1) THEN 
        kThermal = -1 
        DO iG=1,kMaxPtsJac 
          radBackgndThermdT(iG) = 0.0 
          radSolardT(iG) = 0.0 
          END DO 
        CALL Find_BT_rad(raInten,radBTdr,raFreq, 
     $                   radBackgndThermdT,radSolardT) 
        END IF 
 
      write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...' 
      CALL DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq, 
     $   iAtm,iNumLayer,iaaRadLayer, rSatAngle,raLayAngles,raSun,raaExtTemp, 
     $   raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd, 
     $   iProfileLayers,raPressLevels) 
      write(kStdWarn,*)'initializing Jacobian loops ...' 
      CALL Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,raLayAngles, 
     $      rTSpace,rTSurface,raUseEmissivity,raSurface,raSun,raThermal, 
     $      raaOneMinusTau,raaTau,raaLay2Gnd,raaRad,raaGeneral) 

      DO iLay = 1,iNumLayer
        iL = iaRadlayer(iLay)        
        IF (iLay .EQ. iNumLayer) THEN 
          DO iFr = 1,1
            raLMm1_toGnd(iFr) = 1.0 
            END DO 
        ELSE 
          DO iFr = 1,1
            raLMm1_toGnd(iFr) = raaLay2Gnd(iFr,iLay+1) 
            END DO 
          END IF 
        END DO

      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. -2) 
     $    .OR. (iWhichJac .EQ. 20)) THEN
        DO iG=1,iNumGases 
c for each of the iNumGases whose ID's <= kMaxDQ 
c have to do all the iNumLayer radiances 
          iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob) 
          IF (iGasJacList .GT. 0) THEN 
            iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases) 
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), 
     $         raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer) 
            DO iLay = iNumLayer,1,-1 
              rWeight = raaMix(iaaRadLayer(iAtm,iLay),iG) 
              IF (iLay .EQ. iNumLayer) THEN 
                rWeight = rWeight*rFracBot 
              ELSEIF (iLay .EQ. 1) THEN 
                rWeight = rWeight*rFracTop 
                END IF 
              write(kStdWarn,*)'gas d/dq gas# layer#',iG,iLay,
     $               iaaRadLayer(iAtm,iLay) 
              CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag,
     $          iDoAdd,iErr)
              IF (iDoAdd .GT. 0) THEN             
                CALL JacobGasAmtFM1UP(raFreq,raSun,raaRad,iGasJacList,iLay, 
     $            iNumGases,
     $            iaaRadLayer,iAtm,iNumLayer, 
     $            raaOneMinusTau,raaTau,raaaAllDQ,raResults, 
     $            raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight) 
                IF (kSolar .GE. 0) THEN
                  CALL SolarScatterGasJacobianUp(
     $                  iNumlayer,
     $                  iTag,iLay,iaRadLayer(iLay),iGasJacList,raFreq,raSunTOA,
     $                  raaLay2Gnd,raLayAngles,raSunAngles,
     $                  raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                  raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                  raaSolarScatter1Lay,raaSolarScatterCumul,
     $                  raResults)
                  END IF
                iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1 
                CALL doJacobOutput(iLowest,raFreq,raResults, 
     $             radBTdr,raaAmt,raInten,iaGases(iG),iWhichLayer,iGasPosn) 
              ELSE
                DO iFr = 1,kMaxPts
                  raResults(iFr) = 0.0
                  END DO
                END IF
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
              END DO 
            END IF 
          END DO 
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iG=1,iNumGases
c for each of the iNumGases whose ID's <= kMaxDQ
c have to do all the iNumLayer radiances
          iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
          IF (iGasJacList .GT. 0) THEN
            iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),
     $           raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = iNumLayer,1,-1
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              END DO
            END IF
          END DO
        END IF

      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 20) .OR.
     $    (iWhichJac .EQ. -2)) THEN
        DO iG = iNumGases+1,iNumGases+2
c for each of the iNumGases whose ID's <= kMaxDQ
c have to do all the iNumLayer radiances
          iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
          IF (iGasJacList .GT. 0) THEN
            iGasPosn = -1   !!!for JacobOutput
            IF (iG .EQ. iNumGases+1) THEN
              iIWPorDME = +1
            ELSEIF (iG .EQ. iNumGases+2) THEN
              iIWPorDME = -1
              END IF
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),
     $            raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
            DO iLay=iNumLayer,1,-1
              rWeight = 1.0
              IF (iG .EQ. iNumGases+1) THEN
                write(kStdWarn,*)'IWP d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
              ELSEIF (iG .EQ. iNumGases+2) THEN
                write(kStdWarn,*)'DME d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
                END IF
              IF (iaCldLayer(iLay) .LE. 0) THEN
                !!no cloud here, so indpt of cloud params!!!
                DO iFr = 1,kMaxPts
                  raResults(iFr) = 0.0
                  END DO
              ELSE
                CALL JacobCloudAmtUp(raFreq,raaRad,raaRadDT,
     $           iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $           raaOneMinusTau,raaTau,
     $           raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $           raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $           raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,
     $           raaLay2Gnd,raResults,raThermal,
     $           rSatAngle,raLayAngles,
     $           raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $           iOffSet)
                IF (kSolar .GE. 0) THEN
                  CALL SolarScatterCldAmtJacobianUp(
     $                   iNumlayer,
     $                   iTag,iLay,iaRadLayer(iLay),iIWPorDME,raFreq,raSunTOA,
     $                   raaLay2Gnd,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,raaSolarScatterCumul,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $                   raResults)
                  END IF
                END IF
              iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1 
              CALL doJacobOutput(iLowest,raFreq,raResults, 
     $           radBTdr,raaAmt,raInten,iaGasesTemp(iG),iWhichLayer,iGasPosn) 
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              END DO
            END IF
          END DO
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iG=iNumGases+1,iNumGases+2
c for each of the iNumGases whose ID's <= kMaxDQ
c have to do all the iNumLayer radiances
          iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
          IF (iGasJacList .GT. 0) THEN
            iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),
     $           raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = iNumLayer,1,-1
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              END DO
            END IF
          END DO
        END IF

c then do the temperatures d/dT
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),
     $                       rDelta,iAtm,0,iNumLayer)
      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 30) .OR.
     $    (iWhichJac .EQ. -2) .OR. (iWhichJac .EQ. 32)) THEN
        DO iLay = iNumLayer,1,-1   
          IF (iNatm .GT. 1) THEN 
            rWeight=0.0 
            DO iG=1,iNumGases 
              rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iLay),iG)           
              END DO 
            rWeight = rWeight/(iNumGases*1.0) 
          ELSE 
            rWeight=1.0 
            END IF 
c for each of the iNumLayer radiances, cumulatively add on all  
c iNumGases contributions (this loop is done in JacobTemp) 
          write(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay) 
          CALL JacobTempFM1UP(raFreq,raSun,raaRad,raaRadDT,iLay, 
     $            iaaRadLayer,iAtm,iNumLayer, 
     $            raaOneMinusTau,raaTau,raaAllDT,raResults, 
     $            raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight) 
          IF (kSolar .GE. 0) THEN
            CALL SolarScatterTempJacobianUp(
     $                   iNumlayer,
     $                   iTag,iLay,iaRadLayer(iLay),raFreq,raSunTOA,
     $                   raaLay2Gnd,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,raaSolarScatterCumul,
     $                   raResults)
            END IF
          iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1 
          CALL doJacobOutput(iLowest,raFreq,raResults, 
     $                     radBTdr,raaAmt,raInten,0,iWhichLayer,-1) 
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
          END DO 
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iLay = iNumLayer,1,-1
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF

c do the weighting functions 
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), 
     $                       rDelta,iAtm,-10,iNumLayer) 
      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 40) .OR.
     $    (iWhichJac .EQ. -2)) THEN
        DO iLay = iNumLayer,1,-1 
          write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay) 
          CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles, 
     $     iaaRadLayer,iAtm,raaLay2Gnd,raaExtTemp,raResults,rFracTop,rFracBot) 
c does not make sense to multiply the weighting fcns with gas amounts etc 
c      iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1 
c      CALL doJacobOutput(raFreq,raResults,radBTdr,raaAmt,raInten,0, 
c     $                               iWhichLayer) 
c so just output the weighting functions 
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
          END DO 
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iLay = iNumLayer,1,-1
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF
 
c computing Jacobians wrt surface parameters is meanigless .. output 0's 
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), 
     $                       rDelta,iAtm,-20,4) 
      DO iG=1,kMaxPts 
        raResults(iG) = 0.0 
        END DO 
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults) 
 
      RETURN 
      END 

c************************************************************************
c this subroutine adds on the solar contribution to the Cld Amt Jacobian
      SUBROUTINE SolarScatterCldAmtJacobianUp(
     $                   iNumlayer,
     $                   iTag,iLay,iLM,iIWPorDME,raFreq,raSunTOA,
     $                   raaLay2Gnd,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,raaSolarScatterCumul,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iTag,iLay,iLM,iIWPorDME,iNumLayer
      INTEGER iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
      REAL raaSolarScatterCumul(kMaxPts,kProfLayer) 
      REAL raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME) 
      REAL raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP) 
      REAL raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP) 
      REAL raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME) 
      REAL raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME) 
      REAL raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME) 
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iFr,iDoSolar
      REAL muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)
      REAL hg2_real,hg2_real_deriv_wrt_g

      CALL InitSomeSunScatUp(
     $        iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,
     $        raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,
     $        raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iIWPorDME .EQ. +1) THEN
        !! this is IWP jacobian
        IF (iaCldLayer(iLay) .EQ. 1) THEN
          DO iFr = 1,kMaxPts
            rW = exp(-raaExtTemp(iFr,iLM)/muSun)

            !!this is the change in optical depths wrt iwp
            rT = exp(-raaExtTemp(iFr,iLM)/muSat)/muSat - 
     $           exp(-raaExtTemp(iFr,iLM)/muSun)/muSun
            rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW

            rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*
     $        (exp(-raaExtTemp(iFr,iLM)/muSat)-exp(-raaExtTemp(iFr,iLM)/muSun))

            !!this is the change in ssalb wrt iwp
            rT = exp(-raaExtTemp(iFr,iLM)/muSat) - 
     $           exp(-raaExtTemp(iFr,iLM)/muSun)
            rY = raaSSAlbJacobIWP(iFr,iLM)*rT*rW/raaExtJacobIWP(iFr,iLM)

            !! now add rZ and rY and divide by alpha_j 
            !! to get (1/alpha_j) d(alpha_j))
            rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)

            raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
            END DO
          END IF

        DO iFr = 1,kMaxPts
          raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
csun          raResults(iFr) = raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
          raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
          END DO

      ELSEIF (iIWPorDME .EQ. -1) THEN
        !! this is DME jacobian
        IF (iaCldLayer(iLay) .EQ. 1) THEN
          DO iFr = 1,kMaxPts
            rW = exp(-raaExtTemp(iFr,iLM)/muSun)

            !!this is the change in optical depths wrt dme
            rT = exp(-raaExtTemp(iFr,iLM)/muSat)/muSat - 
     $           exp(-raaExtTemp(iFr,iLM)/muSun)/muSun
            rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW

            rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*
     $        (exp(-raaExtTemp(iFr,iLM)/muSat)-exp(-raaExtTemp(iFr,iLM)/muSun))

            !!this is the change in ssalb wrt dme
            rT = exp(-raaExtTemp(iFr,iLM)/muSat) - 
     $           exp(-raaExtTemp(iFr,iLM)/muSun)
            rY = raaSSAlbJacobDME(iFr,iLM)*rT*rW/raaExtJacobDME(iFr,iLM)

            !! now add rZ and rY and divide by alpha_j 
            !! to get (1/alpha_j) d(alpha_j))
            rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)

            raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
            END DO

          !!also add on the d(hg)/d(dme) contribution!!!!
          DO iFr = 1,kMaxPts
            rZ = 1/hg2_real(-muSun,-muSat,raaAsymTemp(iFr,iLM))*
     $             raaPhaseJacobASYM(iFr,iLM)
            rZ = rZ*raaAsymJacobDME(iFr,iLM)/raaExtJacobDME(iFr,iLM)
            raTemp1(iFr) = raTemp1(iFr) + 
     $                     rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
            END DO
          END IF

        DO iFr = 1,kMaxPts
          raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
csun          raResults(iFr) = raTemp(iFr)*raaExtJacobDME(iFr,iLM)
          raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobDME(iFr,iLM)
          END DO
      END IF

      RETURN
      END 

c************************************************************************
c this subroutine adds on the solar contribution to the Temp Jacobian
      SUBROUTINE SolarScatterTempJacobianUp(
     $                   iNumlayer,
     $                   iTag,iLay,iLM,raFreq,raSunTOA,
     $                   raaLay2Gnd,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,raaSolarScatterCumul,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iTag,iLay,iLM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
      REAL raaSolarScatterCumul(kMaxPts,kProfLayer) 
      INTEGER iNumLayer
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iFr,iDoSolar
      REAL muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)

      CALL InitSomeSunScatUp(
     $        iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,
     $        raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,
     $        raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iaCldLayer(iLay) .EQ. 1) THEN
        DO iFr = 1,kMaxPts
          rW = exp(-raaExtTemp(iFr,iLM)/muSun) 

          !!this is the change in optical depths wrt T
          rT = exp(-raaExtTemp(iFr,iLM)/muSat)/muSat - 
     $         exp(-raaExtTemp(iFr,iLM)/muSun)/muSun
          rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW

          rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun* 
     $      (exp(-raaExtTemp(iFr,iLM)/muSat)-exp(-raaExtTemp(iFr,iLM)/muSun)) 

          !!this is the change in ssalb wrt gas temperature
          rT = exp(-raaExtTemp(iFr,iLM)/muSat) -
     $         exp(-raaExtTemp(iFr,iLM)/muSun)
          rY = (1 + raaAsymTemp(iFr,iLM))/2.0
          rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
          rY = rT*rW*rY/raaExtTemp(iFr,iLM)

          !! now add rZ and rY and divide by alpha_j 
          !! to get (1/alpha_j) d(alpha_j))
          rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)

          raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
          END DO
        END IF

c add on to the raResults tally
      DO iFr = 1,kMaxPts
        raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
csun        raResults(iFr) = raTemp(iFr)*raaAllDT(iFr,iLM)
        raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaAllDT(iFr,iLM)
        END DO

      RETURN
      END 

c************************************************************************
c this subroutine adds on the solar contribution to the Amt Jacobian
      SUBROUTINE SolarScatterGasJacobianUp(
     $                   iNumlayer,
     $                   iTag,iLay,iLM,iG,raFreq,raSunTOA,
     $                   raaLay2Gnd,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,raaSolarScatterCumul,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iTag,iLay,iLM,iG,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
      REAL raaSolarScatterCumul(kMaxPts,kProfLayer) 
      INTEGER iNumLayer
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iFr,iDoSolar
      REAL muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)

      CALL InitSomeSunScatUp(
     $        iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,
     $        raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,
     $        raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iaCldLayer(iLay) .EQ. 1) THEN
        DO iFr = 1,kMaxPts
          rW = exp(-raaExtTemp(iFr,iLM)/muSun) 

          !!this is the change in optical depths wrt T
          rT = exp(-raaExtTemp(iFr,iLM)/muSat)/muSat - 
     $         exp(-raaExtTemp(iFr,iLM)/muSun)/muSun
          rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW

          rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun* 
     $      (exp(-raaExtTemp(iFr,iLM)/muSat)-exp(-raaExtTemp(iFr,iLM)/muSun)) 

          !!this is the change in ssalb wrt gas temperature
          rT = exp(-raaExtTemp(iFr,iLM)/muSat) -
     $         exp(-raaExtTemp(iFr,iLM)/muSun)
          rY = (1 + raaAsymTemp(iFr,iLM))/2.0
          rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
          rY = rT*rW*rY/raaExtTemp(iFr,iLM)

          !! now add rZ and rY and divide by alpha_j 
          !! to get (1/alpha_j) d(alpha_j))
          rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)

          raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
          END DO
        END IF

c add on to the raResults tally
      DO iFr = 1,kMaxPts
        raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
csun        raResults(iFr) = raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
        raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
        END DO

      RETURN
      END 

c************************************************************************
c this subroutine does some initializations common to the Cloud Jac routines
c in particular, it initializes the cumulative jacobian for (cloud) layers
c beneath this layer (iLay) in question
      SUBROUTINE InitSomeSunScatUp(
     $        iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,
     $        raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,
     $        raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iLM,iLay,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer),iNumLayer
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
      REAL raaSolarScatterCumul(kMaxPts,kProfLayer) 
c output vars
      REAL muSun,muSat
      REAL raTemp(kMaxPts),raTemp1(kMaxPts),raLMm1_toGnd(kMaxPts)

c local vars
      INTEGER iFr,iJ,iLJ
      REAL rX

c iLM = iaRadLayer(iLay)
      DO iFr = 1,kMaxPts
        raTemp1(iFr) = 0.0
        END DO

c stuff on Jan 5, 2006
      DO iFr = 1,kMaxPts
        raTemp(iFr) = 0.0
        END DO
      DO iJ = 1,iLay-1
        iLJ = iaRadLayer(iJ)
        IF (iaCldLayer(iJ) .EQ. 1) THEN
          muSat = cos(raLayAngles(iLJ) * kPi/180)
          muSun = cos(raSunAngles(iLJ) * kPi/180)
          rX = (-1/muSat)*(1+muSat/muSun)
          rX = (-1/muSun)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + 
     $                    raaSolarScatter1Lay(iFr,iLJ)*raaLay2Gnd(iFr,iJ+1)*rX
            END DO
          END IF
        END DO
      DO iJ = iLay+1,iNumLayer
        iLJ = iaRadLayer(iJ)
        IF (iaCldLayer(iJ) .EQ. 1) THEN
          muSat = cos(raLayAngles(iLJ) * kPi/180)
          muSun = cos(raSunAngles(iLJ) * kPi/180)
          rX = (-1/muSat)*(1+muSat/muSun)
          rX = (-1/muSat)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + 
     $                    raaSolarScatter1Lay(iFr,iLJ)*raaLay2Gnd(iFr,iJ+1)*rX
            END DO
          END IF
        END DO

c output   these vars for the rest of the routine
      muSat = cos(raLayAngles(iLM) * kPi/180)
      muSun = cos(raSunAngles(iLM) * kPi/180)

      IF (iLM .EQ. iNumLayer) THEN
        DO iFr = 1,kMaxPts
          raLMm1_toGnd(iFr) = 1.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          raLMm1_toGnd(iFr) = raaLay2Gnd(iFr,iLay+1)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this does the jacobian wrt IWP or DME
      SUBROUTINE JacobCloudAmtUp(raFreq,raaRad,raaRadDT,
     $          iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $          raaOneMinusTau,raaTau,
     $          raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $          raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $          raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,
     $          raaLay2Gnd,raResults,raThermal,
     $          rSatAngle,raLayAngles,
     $          raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $          iOffSet)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatAngle is the satellite viewing angle
c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c daaDT,daaDQ are the d/dq,d/dT matrices
c raaLay2Gnd   is the layer-to-gnd abs coeff matrix
c raaRad has the Planck radiances
c raaRadDT has the d/DT (Planck radiances)
c raaOneMinusTau has 1-tau(satellite angle)
c raaOneMinusTauTh has 1-tau(thermal angle)
c iG has the gas number (1 .. iNumGasesTemp)
c iLay has info on how to find the radiating layer number (1..kProfLayerJac)
c raFreq has the frequencies
c raResults has the results
c raThermal are the downwelling Solar,thermal contributions
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
c raaGeneral,raaGeneralTh have the general results (looping over layers)
      REAL raaGeneral(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle
      REAL raLayAngles(kProfLayer)
      REAL raThermal(kMaxPts)
      INTEGER iNumGasesTemp,iAtm
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac),raFreq(kMaxPts)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)
      INTEGER iG,iLay,iNumLayer,iOffSet
c output vars
      REAL raResults(kMaxPtsJac)

c this is for the scattering parts
      REAL raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
      REAL raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)  !scattering d/d(IWP)
      REAL raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
      REAL raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
      REAL raaSSAlbJacobDME(kMaxPts,kProfLayerJac)  !scattering d/d(DME)
      REAL raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac)  !phase fnc jacobs wrt g
      INTEGER iIWPorDME

c local variables
      INTEGER iFr,iJ1,iM1,MP2Lay
      REAL raTemp(kMaxPtsJac),muSat
      REAL raResultsTh(kMaxPtsJac)

c figure out which of 1..100 this current radiating layer corresponds to
c bleh
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the sat angle weight factor
      muSat = 1.0/cos(rSatAngle*kPi/180.0)
      muSat = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

c read the appropriate layer from general results
      DO iFr=1,kMaxPts
        raResults(iFr) = raaGeneral(iFr,iLay)
        END DO

c note that      iLay + iOffset === iaRadlayer(iLay)
c set the constant factor we have to multiply results with
      IF (iIWPorDME .EQ. 1) THEN
        DO iFr=1,kMaxPts
          raTemp(iFr) = raaExtJacobIWP(iFr,iLay+iOffSet)
          raTemp(iFr) = raaExtJacobIWP(iFr,iM1)
          END DO
      ELSEIF (iIWPorDME .EQ. -1) THEN
        DO iFr=1,kMaxPts
          raTemp(iFr) = raaExtJacobDME(iFr,iLay+iOffSet)
          raTemp(iFr) = raaExtJacobDME(iFr,iM1)
          END DO
        END IF
      CALL MinusOne(raTemp,raResults)

c add on the the derivatives wrt radiating layer
      IF (iLay .LT. iNumLayer) THEN
c this is not the topmost layer
        iJ1 = iLay
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)+
     $      raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)
          END DO 
      ELSE IF (iLay .EQ. iNumLayer) THEN
c do the topmost layer correctly
        iJ1 = iLay
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)+
     $          raTemp(iFr)*raaTau(iFr,iJ1)*raaRad(iFr,iJ1)
          END DO 
        END IF

c now multiply results by the 1/cos(viewing angle) factor
      IF (abs(muSat-1.0000000) .GE. 1.0E-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*muSat
          END DO
        END IF

c see if we have to include thermal backgnd
c ignore for now
c      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
c        CALL JacobTHERMALAmtFM1(raFreq,raaRad,
c     $       iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,
c     $       raUseEmissivity,raTemp,raaLay2Gnd,
c     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
cc now add on the effects to raResults
c        DO iFr=1,kMaxPts
c          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
c          END DO
c        END IF

      RETURN
      END 

c************************************************************************
