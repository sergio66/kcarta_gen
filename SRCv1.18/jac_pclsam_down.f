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
c iaaRadLayer(kMaxAtm,kProfLayer),raaAbs(kMaxPts,kMixFilRows)
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

c put in a lot of rEps in the SUBROUTINE solarscatter/gas/temp/iwpdme as things
c could go nuts if raaIWP = 0 in between cloud layers

c************************************************************************
c****** THESE ARE THE SCATTERING JACOBIANS FOR THE DOWN LOOK INSTR ******
c************************************************************************
c this subroutine does the Jacobians for downward looking instrument
      SUBROUTINE DownwardJacobian_Scat(raFreq,
     $            iProfileLayers,raPressLevels,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX,
     $               iNpMix,iTag,iActualTag,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,
     $            iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

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
      CHARACTER*80 caJacobFile
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
      REAL iwpMAX(MAXNZ)

c local variables
      INTEGER iNumGasesTemp,iaGasesTemp(kMaxGas)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
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
      INTEGER iG,iM,iIOUN,iLowest
      INTEGER DoGasJacob,iGasJacList
      INTEGER WhichGasPosn,iGasPosn
c for cloud stuff!
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1) 
      INTEGER iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot 
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iiDiv,iaRadLayer(kProfLayer),iLay,iIWPorDME,iL
      INTEGER iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
      REAL r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck,muSun,rSunAngle
      INTEGER iSolarRadOrJac,MP2Lay
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer)

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

      rSunAngle = raSunAngles(MP2Lay(iaaRadLayer(1,iAtm)))
c change to radians 
      rSunAngle=(rSunAngle*kPi/180.0) 
      muSun = cos(rSunAngle)    

      rSunTemp  = kSunTemp
      rOmegaSun = kOmegaSun
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

      iIOUN = kStdJacob

      iLowest = iaaRadLayer(iAtm,1)
      iLowest = MOD(iLowest,kProfLayer)

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

cccccccccccccccccccc set these all important variables ****************  
        IF (iaRadLayer(1) .LT. kProfLayer) THEN  
          iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1  
          iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1  
          iiDiv = 0  
        ELSE  
          !!essentially do mod(iaRadLayer(1),kProfLayer)  
          iiDiv = 1            
 1010     CONTINUE  
          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN  
            iiDiv = iiDiv + 1  
            GOTO 1010  
            END IF  
          iiDiv = iiDiv - 1  
          iLay = iiDiv  
          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)  
          iLocalCldTop = iCldTopkCarta - iiDiv + 1  
          iLocalCldBot = iCldBotkCarta - iiDiv + 1  
          iiDiv = iLay  
          END IF  

      iOffSet = (kProfLayer-iNumLayer)
c      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1  
      DO iLay = iCldBotkCarta,iCldTopkCarta
        !!! this is for d/dq, d/dT
        iaCldLayer(iLay) = 1  
        !!!this is for d/d(DME), d/d(IWP)
        iaCldLayerIWPDME(iLay - (kProfLayer-iNumLayer)) = 1
        END DO 

      IF (kSolar .GE. 0) THEN
        iSolarRadOrJac = +1
        CALL SolarScatterIntensity_Downlook(
     $      iDoSolar,raFreq,iaCldLayer,
     $      raSunAngles,raLayAngles,0.0,0.0,
     $      iNumLayer,iaRadLayer,
     $      raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot,
     $      iTag,iSolarRadOrJac,raaSolarScatter1Lay)

c        DO iLay = 1,iNumLayer
c          iL = iaRadlayer(iLay)
c          print *,'<<<<<>>>>',iLay,iL,iaCldLayer(iL),iaCldLayerIWPDME(iL),
c     $            raSunAngles(iL),raLayAngles(iL),
c     $            raaExtTemp(1,iL),raaSSAlbTemp(1,iL),raaAsymTemp(1,iL),
c     $            raaSolarScatter1Lay(1,iL),
c     $            raaPhaseJacobASYM(1,iL),raaExtJacobIWP(1,iL),
c     $            raaExtJacobDME(1,iL)
c          END DO

        END IF

cccccccccccccccccccc set these all important variables ****************  

c initialise the layer-to-space matrix
      CALL AtmosLayer2Space(raaLay2Sp,
     $  rSatAngle,raaExtTemp,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

      write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
      CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq,
     $              iAtm,iNumLayer,iaaRadLayer,
     $              rSatAngle,raLayAngles,raaExtTemp,
     $              raaRad,raaRadDT,raaOneMinusTau,raaTau,
     $              raaLay2Gnd,iProfileLayers,raPressLevels)

      write(kStdWarn,*)'initializing Jacobian loops ...'
      CALL Loop_LookDown(iaaRadLayer,
     $        iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl,
     $        raUseEmissivity,raSurface,raSun,raSunAngles,raThermal,
     $        raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
        write(kStdWarn,*)'initializing Jacobian thermal loops ...'
        CALL Loop_thermaldown(raaRad,rSatAngle,raLayAngles,
     $          iProfileLayers,raPressLevels,
     $          iaaRadLayer,iAtm,iNumLayer,raaExtTemp,
     $          raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)

        END IF

      IF (kJacobOutPut .GE. 1) THEN
c have to set the backgnd thermal, solar radiances so that we can do 
c dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
        IF (kThermal .GE. 0) THEN
          DO iG=1,kMaxPts
            radBackgndThermdT(iG) = raThermal(iG)
            END DO
          END IF

        IF (kSolar .GE. 0) THEN           
c compute the Solar contribution
          iM=1
c remember that raSun is that at the gnd -- we have to propagate this back 
c up to the top of the atmosphere
c note that here we are calculating the SOLAR contribution 
          DO iG=1,kMaxPts
            radSolardT(iG) = raUseEmissivity(iG)*
     $                     raaLay2Sp(iG,iM)*raSun(iG)
            END DO
          END IF

        CALL Find_BT_rad(raInten,radBTdr,raFreq,
     $                   radBackgndThermdT,radSolardT)
        END IF

      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. -2)
     $    .OR. (iWhichJac .EQ. 20)) THEN
        DO iG = 1,iNumGases
c for each of the iNumGases whose ID's <= kMaxDQ
c have to do all the iNumLayer radiances
          iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
          IF (iGasJacList .GT. 0) THEN
            iGasPosn = WhichGasPosn(iaGasesTemp(iG),iaGasesTemp,iNumGasesTemp)
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),
     $          raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
            DO iM=1,iNumLayer
              rWeight = raaMix(iaaRadLayer(iAtm,iM),iG)
              IF (iM .EQ. 1) THEN
                rWeight = rWeight*rFracBot
              ELSEIF (iM .EQ. iNumLayer) THEN
                rWeight = rWeight*rFracTop
                END IF
              write(kStdWarn,*)'gas d/dq gas# lay#',iG,iM,iaaRadLayer(iAtm,iM)
              !! see if this gas does exist for this chunk
              CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag,
     $           iDoAdd,iErr)
              IF (iDoAdd .GT. 0) THEN             
                CALL JacobGasAmtFM1(raFreq,raaRad,raaRadDT,iGasJacList,
     $            iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $            raaOneMinusTau,raaTau,raaaAllDQ,
     $            raaLay2Sp,raResults,raThermal,raaLay2Gnd,
     $            rSatAngle,raLayAngles,
     $            raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $            rWeight)
                IF (kSolar .GE. 0) THEN
                  CALL SolarScatterGasJacobian(
     $                     iTag,iM,iaRadLayer(iM),iGasJacList,raFreq,raSunTOA,
     $                     raaLay2Sp,raLayAngles,raSunAngles,
     $                     raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                     raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                     raaSolarScatter1Lay,
     $                     raResults)
                  END IF
                CALL doJacobOutput(iLowest,raFreq,
     $            raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
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
            DO iLay = 1,iNumLayer
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              END DO
            END IF
          END DO
        END IF

      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. -2) .OR. 
     $    (iWhichJac .EQ. 20)) THEN
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
            DO iM=1,iNumLayer
              rWeight = 1.0
              IF (iG .EQ. iNumGases+1) THEN
                write(kStdWarn,*)'IWP d/dq lay#',iM,iaaRadLayer(iAtm,iM)
              ELSEIF (iG .EQ. iNumGases+2) THEN
                write(kStdWarn,*)'DME d/dq lay#',iM,iaaRadLayer(iAtm,iM)
                END IF
c this is more correct as it only puts out d/d(IWP) and d/d(DME) where there
c actually is cloud
c            IF (iaCldLayerIWPDME(iM) .EQ. 0 .OR. 
c     $          iwpMAX(iaaRadLayer(iAtm,iM)) .LE. 1.0e-10) THEN
c very fun thing happens if you let this line compile ... and look at d/d(IWP)
c it mimics T(z)!!!!!
              IF (iaCldLayerIWPDME(iM) .EQ. 0) THEN
                !!no cloud here, so indpt of cloud params!!!
                DO iFr = 1,kMaxPts
                  raResults(iFr) = 0.0
                 END DO
              ELSE
                CALL JacobCloudAmtFM1(raFreq,raaRad,raaRadDT,
     $            iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $            raaOneMinusTau,raaTau,
     $            raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,
     $            raaLay2Sp,raResults,raThermal,raaLay2Gnd,
     $            rSatAngle,raLayAngles,
     $            raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $            iOffSet)
                IF (kSolar .GE. 0) THEN
                  CALL SolarScatterCldAmtJacobian(
     $                   iTag,iM,iaRadLayer(iM),iIWPorDME,raFreq,raSunTOA,
     $                   raaLay2Sp,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,iwpMAX,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $                   raResults)
                  END IF
                END IF
              CALL doJacobOutput(iLowest,raFreq,
     $          raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
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
            DO iLay = 1,iNumLayer
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
        DO iM=1,iNumLayer
c for each of the iNumLayer radiances, cumulatively add on all 
c iNumGases contributions (this loop is done in JacobTemp)
          write(kStdWarn,*)'temp d/dT layer# = ',iM,iaaRadLayer(iAtm,iM)
c don't have to worry about gases 210,202 as their "weight" for d/dT = 0.0
          IF (iNatm .GT. 1) THEN
            rWeight = 0.0
            DO iG=1,iNumGases
              rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iM),iG)            
              END DO
            rWeight = rWeight/(iNumGases*1.0)
          ELSE
            rWeight = 1.0
            END IF
          CALL JacobTempFM1(raFreq,raaRad,raaRadDT,iM,iNumGases,
     $            iaaRadLayer,iAtm,iNumLayer,
     $            raUseEmissivity,
     $            raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults,
     $            raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,
     $            raaGeneralTh,raaOneMinusTauTh,rWeight)
         IF (kSolar .GE. 0) THEN
            CALL SolarScatterTemperatureJacobian(
     $                   iTag,iM,iaRadLayer(iM),raFreq,raSunTOA,
     $                   raaLay2Sp,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,
     $                   raResults)
            END IF
          CALL doJacobOutput(iLowest,raFreq,raResults,
     $                    radBTdr,raaAmt,raInten,0,iM,-1)
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iLay = 1,iNumLayer
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF

c do the weighting functions
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),
     $                       rDelta,iAtm,-10,iNumLayer)
      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 40) .OR.
     $    (iWhichJac .EQ. -2)) THEN
        DO iM=1,iNumLayer
          write(kStdWarn,*)'wgt fcn # = ',iM,iaaRadLayer(iAtm,iM)
          CALL wgtfcndown(iM,iNumLayer,rSatAngle,raLayAngles,
     $      iaaRadLayer,iAtm,raaLay2Sp,raaExtTemp,raResults,rFracTop,rFracBot,
     $      iNLTEStart,raaPlanckCoeff)
c does not make sense to multiply the weighting fcns with gas amounts etc
c        CALL doJacobOutput(iLowest,raFreq,raResults,
c     $                     radBTdr,raaAmt,raInten,0,iM,-1)
c so just output the weighting functions
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iLay = 1,iNumLayer
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF

c finally do the surface sensitivities : d/d(SurfaceTemp), 
c d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),
     $                       rDelta,iAtm,-20,4)
      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 50) .OR. 
     $    (iWhichJac .EQ. -2)) THEN
        iM=1
c computing Jacobians wrt surface parameters are meaningful
        CALL JacobSurfaceTemp(raFreq,iM,
     $      rTSurface,raUseEmissivity,raaLay2Sp,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBTdr,raaAmt,raInten,-1,iM,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobSurfaceEmis(iM,
     $            raSurface,raThermal,raaLay2Sp,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBTdr,raaAmt,raInten,-2,iM,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobBackgndThermal(iM,
     $            raaLay2Sp,raThermal,raResults) 
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBackgndThermdT,raaAmt,raInten,-3,iM,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobSolar(iM,raaLay2Sp,raSun,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radSolardT,raaAmt,raInten,-4,iM,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
          raResults(iFr) = 0.0
          END DO
        DO iLay = 1,4
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF

      
      RETURN
      END

c************************************************************************
c this subroutine adds on the solar contribution to the Cld Amt Jacobian
      SUBROUTINE SolarScatterCldAmtJacobian(
     $                   iTag,iM,iLM,iIWPorDME,raFreq,raSunTOA,
     $                   raaLay2Sp,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,iwpMAX,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      INTEGER iTag,iM,iLM,iIWPorDME
      INTEGER iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts),iwpMAX(MAXNZ)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
      REAL raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME) 
      REAL raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP) 
      REAL raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP) 
      REAL raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME) 
      REAL raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME) 
      REAL raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME) 
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iLay,iFr,iDoSolar
      REAL muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)
      REAL hg2_real,hg2_real_deriv_wrt_g
      REAL rEps   !! so that if w == 0, rz is finite

      rEps = 1.0e-10

      CALL InitSomeSunScat(
     $                   iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,
     $                   raaSolarScatter1Lay,iaRadLayer,
     $                   raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iIWPorDME .EQ. +1) THEN
        !! this is IWP jacobian
        IF (iaCldLayer(iLM) .EQ. 1) THEN
          muSunSat = (muSun * muSat)/(muSun + muSat)
          DO iFr = 1,kMaxPts
            !!this is the transmission
            rT = exp(-raaExtTemp(iFr,iLM)/muSunSat)

            !!this is the change in ssalb wrt iwp
            rY = raaSSAlbJacobIWP(iFr,iLM)*(1.0-rT)/(raaExtJacobIWP(iFr,iLM)+rEps)

            !!this is the change in optical depths wrt iwp
            rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
            rZ = rZ + rEps

            !! now add rZ and rY and divide by alpha_j 
            !! to get (1/alpha_j) d(alpha_j))
            rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))

            rT = exp(-raaExtTemp(iFr,iLM)/muSun)
            raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)

            END DO
          END IF

        DO iFr = 1,kMaxPts
          raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
csun          raResults(iFr) = raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
          raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
          END DO

      ELSEIF (iIWPorDME .EQ. -1) THEN
        !! this is DME jacobian
        IF (iaCldLayer(iLM) .EQ. 1) THEN
          muSunSat = (muSun * muSat)/(muSun + muSat)
          DO iFr = 1,kMaxPts
            !!this is the transmission
            rT = exp(-raaExtTemp(iFr,iLM)/muSunSat)

            !!this is the change in ssalb wrt dme
            rY = raaSSAlbJacobDME(iFr,iLM)*(1.0-rT)/
     $           (raaExtJacobDME(iFr,iLM)+rEps)

            !!this is the change in optical depths wrt dme
            rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
            rZ = rZ + rEps

            !! now add rZ and rY and divide by alpha_j 
            !! to get (1/alpha_j) d(alpha_j))
            rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))
  
            raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
            END DO

          !!also add on the d(hg)/d(dme) contribution!!!!
          DO iFr = 1,kMaxPts
            rZ = 1/hg2_real(-muSun,muSat,raaAsymTemp(iFr,iLM))*
     $             raaPhaseJacobASYM(iFr,iLM)
            rZ = rZ*raaAsymJacobDME(iFr,iLM)/(raaExtJacobDME(iFr,iLM)+rEps)
            raTemp1(iFr) = raTemp1(iFr) + 
     $                     rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
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
      SUBROUTINE SolarScatterTemperatureJacobian(
     $                   iTag,iM,iLM,raFreq,raSunTOA,
     $                   raaLay2Sp,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iTag,iM,iLM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iLay,iFr,iDoSolar
      REAL muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)
      REAL rEps   !! so that if w == 0, rz is finite

      rEps = 1.0-10
      
      CALL InitSomeSunScat(
     $                   iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,
     $                   raaSolarScatter1Lay,iaRadLayer,
     $                   raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iaCldLayer(iLM) .EQ. 1) THEN
        muSunSat = (muSun * muSat)/(muSun + muSat)
        DO iFr = 1,kMaxPts
          !!this is the transmission
          rT = exp(-raaExtTemp(iFr,iLM)/muSunSat)

          !!this is the change in ssalb wrt gas temperature
          rY = (1 + raaAsymTemp(iFr,iLM))/2.0
          rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
          rY = (1.0-rT)*rY/raaExtTemp(iFr,iLM)

          !!this is the change in optical depths wrt temperature
          rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
          rZ = rZ + rEps

          !! now add rZ and rY and divide by alpha_j 
          !! to get (1/alpha_j) d(alpha_j))
          rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))

          raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)

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
      SUBROUTINE SolarScatterGasJacobian(
     $                   iTag,iM,iLM,iG,raFreq,raSunTOA,
     $                   raaLay2Sp,raLayAngles,raSunAngles,
     $                   raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                   raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,
     $                   raaSolarScatter1Lay,
     $                   raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iTag,iM,iLM,iG,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raFreq(kMaxPts),raSunTOA(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
c output vars
      REAL raResults(kMaxPtsJac)

c local vars 
      INTEGER iL,iLay,iFr,iDoSolar
      REAL muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
      REAL rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)
      REAL rEps   !! so that if w == 0, rz is finite

      rEps = 1.0e-10

      CALL InitSomeSunScat(
     $                   iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,
     $                   raaSolarScatter1Lay,iaRadLayer,
     $                   raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

c adjust for change in scattered radn from this layer, if cloudy
      IF (iaCldLayer(iLM) .EQ. 1) THEN
        muSunSat = (muSun * muSat)/(muSun + muSat)
        DO iFr = 1,kMaxPts
          !!this is the transmission
          rT = exp(-raaExtTemp(iFr,iLM)/muSunSat)

          !!this is the change in ssalb wrt gas amount
          rY = (1 + raaAsymTemp(iFr,iLM))/2.0
          rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
          rY = (1.0-rT)*rY/raaExtTemp(iFr,iLM)

          !!this is the change in optical depths wrt gas amount
          rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
          rZ = rZ + rEps

          !! now add rZ and rY and divide by alpha_j to get 
          !! (1/alpha_j) d(alpha_j))
          rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))

          raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
 
          END DO
        END IF

c      IF (raaSSAlbTemp(iFr,iLM) .LT. rEps) THEN 
c        DO iFr = 1,kMaxPts
c          raTemp1(iFr) = 0.0
c          END DO
c        END IF

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
c beneath this layer (iM) in question
      SUBROUTINE InitSomeSunScat(
     $                   iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,
     $                   raaSolarScatter1Lay,iaRadLayer,
     $                   raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iLM,iM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
c output vars
      REAL muSun,muSat,rFac
      REAL raTemp(kMaxPts),raTemp1(kMaxPts),raLMp1_toSpace(kMaxPts)

c local vars
      INTEGER iFr,iJ,iLJ
      REAL rX

c iLM = iaRadLayer(iM)
      DO iFr = 1,kMaxPts
        raTemp1(iFr) = 0.0
        END DO

c stuff on Jan 5, 2006
      muSat = cos(raLayAngles(iLM) * kPi/180)
      muSun = cos(raSunAngles(iLM) * kPi/180)
      rFac = 1.0 + muSat/muSun

      rX = (-1/muSat)*(1+muSat/muSun) 
      !!! oh boy is this wrong? i think it is fine 1/12/06
      DO iFr = 1,kMaxPts
        raTemp(iFr) = 0.0
        END DO
      DO iJ = 1,iM-1
        iLJ = iaRadLayer(iJ)
        IF (iaCldLayer(iLJ) .EQ. 1) THEN
          !muSat = cos(raLayAngles(iLJ) * kPi/180)
          !muSun = cos(raSunAngles(iLJ) * kPi/180)
          !rX = (-1/muSat)*(1+muSat/muSun)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + 
     $                    raaSolarScatter1Lay(iFr,iLJ)*raaLay2Sp(iFr,iJ+1)*rX
            END DO
          END IF
        END DO

c output these vars for the rest of the routine
      muSat = cos(raLayAngles(iLM) * kPi/180)
      muSun = cos(raSunAngles(iLM) * kPi/180)
      rFac = 1.0 + muSat/muSun

      IF (iLM .EQ. kProfLayer) THEN
        DO iFr = 1,kMaxPts
          raLMp1_toSpace(iFr) = 1.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          raLMp1_toSpace(iFr) = raaLay2Sp(iFr,iM+1)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this does the jacobian wrt IWP or DME
      SUBROUTINE JacobCloudAmtFM1(raFreq,raaRad,raaRadDT,
     $          iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $          raaOneMinusTau,raaTau,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,
     $          raaLay2Sp,raResults,raThermal,raaLay2Gnd,
     $          rSatAngle,raLayAngles,
     $          raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $          iOffSet)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatAngle is the satellite viewing angle
c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c daaDT,daaDQ are the d/dq,d/dT matrices
c raaLay2Sp   is the layer-to-space abs coeff matrix
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
      INTEGER iNumGasesTemp,iAtm,iM
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raFreq(kMaxPts)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)
      INTEGER iG,iLay,iNumLayer,iOffSet

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
      REAL raTemp(kMaxPtsJac),rCos
      REAL raResultsTh(kMaxPtsJac)

c figure out which of 1..100 this current radiating layer corresponds to
c bleh
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the sat angle weight factor
      rCos = 1.0/cos(rSatAngle*kPi/180.0)
      rCos = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

c read the appropriate layer from general results
      DO iFr=1,kMaxPts
        raResults(iFr) = raaGeneral(iFr,iLay)
        END DO

c note that      iLay + iOffset === iaRadlayer(iLay)
c set the constant factor we have to multiply results with
      IF (iIWPorDME .EQ. 1) THEN
        DO iFr=1,kMaxPts
          raTemp(iFr) = raaExtJacobIWP(iFr,iLay+iOffSet)
          END DO
      ELSEIF (iIWPorDME .EQ. -1) THEN
        DO iFr=1,kMaxPts
          raTemp(iFr) = raaExtJacobDME(iFr,iLay+iOffSet)
          END DO
        END IF
      CALL MinusOne(raTemp,raResults)

c add on the the derivatives wrt radiating layer
      IF (iLay .LT. iNumLayer) THEN
c this is not the topmost layer
        iJ1 = iLay
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)+
     $      raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Sp(iFr,iJ1)
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
      IF (abs(rCos-1.0000000) .GE. 1.0E-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*rCos
          END DO
        END IF

c see if we have to include thermal backgnd
c ignore for now
c      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
c        CALL JacobTHERMALAmtFM1(raFreq,raaRad,
c     $       iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,
c     $       raUseEmissivity,raTemp,raaLay2Sp,
c     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
cc now add on the effects to raResults
c        DO iFr=1,kMaxPts
c          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
c          END DO
c        END IF

      RETURN
      END 

c************************************************************************
