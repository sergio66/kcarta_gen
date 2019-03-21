! Copyright 2017
! University of Maryland Baltimore County
! All Rights Reserved

! (1)JacobGasAmtFM1,JacobTempFM1 : jacobians from the forward model
!    (includes solar contribution and thermal diffusive contribution)
! (2)Surface Reflectivity = 1/pi for thermal
!    Surface Reflectance for solar is defined by user

! the following variables are not size kMaxPtsJac or kProfLayerJac as they
! are well defined in the other non Jacobian routines
! raFreq(kMaxPts),raUseEmissivity(kMaxPts),raVTemp(kMixFilRows),
! iaaRadLayer(kMaxAtm,kProfLayer),raaExtTemp(kMaxPts,kMixFilRows)
!     $              raSurface,raSun,raThermal,raInten,
!     $              raSunRefl,
!     $              raLayAngles,raSunAngles)

! we also allow the user to compute the temperature jacobians in
! one of three ways
! recall r(v)= sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

MODULE jac_pclsam_up

USE basic_common
use jac_up
use jac_down
use jac_main
use singlescatter
use clear_scatter_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
!****** THESE ARE THE SCATTERING JACOBIANS FOR THE UP LOOK INSTR ******
!************************************************************************
! this subroutine does the Jacobians for downward looking instrument
    SUBROUTINE UpwardJacobian_ScatPCLSAM(raFreq,iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp, ctype2,rFracx,&
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
    iNpMix,iTag,iActualTag, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA, iPrintAllPCLSAMJacs, &
    iNLTEStart,raaPlanckCoeff, &
    raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rDelta is the wavenumber step size
! iJacob,iaJacob tell which gases to do d/dq for
! caJacobFile is the name of the file to save the Jacobian output to
! iFileID is which 25 cm-1 block being output
! iNumGases is the number of gases to include
! iaGases is the integer array of GasID's
! iNumLayer is the number of layers in the atmosphere # iAtm
! iaaRadLayer is the list of radiating mixed paths to be included
! raVTemp are the layer temperatures

! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
! raaAllDT has the cumulative d/dT coeffs for current freq block
!     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
!        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
!        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

! raaSumAbCoeff is the cumulative absorption coeffs
! raaAmt  has the gas profiles
! raInten has the radiance vector
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi or kSolarRefl
! raLayAngles = layer dependent satellite view angles
! raSunAngles = layer dependent sun view angles
! raaMix is the mixing table
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rDelta,raPressLevels(kProfLayer+1)
    REAL :: raSunRefl(kMaxPts),rFracTop,rFracBot
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
       raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts),rFracx
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
! this is to help the cumulative sums over clouds
    REAL :: raaaAllJacQout(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllJacTout(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllWgtOut(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllSurfOut(kMaxPtsJac,4)
    CHARACTER(80) :: caJacobFile
    INTEGER :: iJacob,iaJacob(kMaxDQ),iProfileLayers,iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,ctype2,iaGases(kMaxGas)
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is for the scattering parts
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    INTEGER :: iNpMix,iCLoudySky,iPrintAllPCLSAMJacs

! local variables
    INTEGER :: iNumGasesTemp,iaGasesTemp(kMaxGas)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
    REAL :: radSolardT(kMaxPtsJac),rWeight
    INTEGER :: iG,iLay,iIOUN,iLowest
    INTEGER :: iGasJacList,iGasPosn
! for cloud stuff!
    REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
    INTEGER :: iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iCloudLayerTop,iCloudLayerBot
    INTEGER :: iiDiv,iaRadLayer(kProfLayer),iIWPorDME,iL
    INTEGER :: iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
    REAL :: r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck
    REAL :: muSun,muSat,rSunAngle
    INTEGER :: iSolarRadOrJac
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer),raLMm1_toGnd(kMaxPts)
    INTEGER :: iWhichLayer

    INTEGER :: iDefault,iWhichJac,iFr,iJ
    INTEGER :: iDoAdd,iErr

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

    IF (iDefault /= iWhichJac) THEN
      print *,'iDefault,iWhichJac = ',iDefault,iWhichJac
    END IF

! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)

! as we are never directly loooking at the sun, there is a geometry factor
    rSunAngle = raSunAngles(iaaRadlayer(1,iAtm))
    rOmegaSun = kOmegaSun
    IF (iDoSolar >= 0) THEN
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument .. daytime'
    ELSE IF (iDoSolar < 0) THEN
      rSunTemp = 0.0
      write(kStdWarn,*)'upward looking instrument .. nitetime'
    END IF

    muSun = 1.0       !!!default
    IF (iDoSolar >= 0) THEN
      muSun = cos(rSunAngle*kPi/180.0)
    END IF

    iDoSolar = kSolar
    IF (iDoSolar == 0) THEN
      !use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      raSunTOA = rattorad(raFreq,rSunTemp)*rOmegaSun*muSun
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      !read in data from file
      CALL ReadSolarData(raFreq,raSunTOA,iTag)
      raSunTOA = raSunTOA*rOmegaSun*muSun
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

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
      IF (iaRadLayer(iLay) > iNpmix) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

! cccccccccccccccccc set these all important variables ****************
    DO iLay = 1,kProfLayer
      iaCldLayer(iLay) = -1
      iaCldLayerIWPDME(iLay) = -1
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (iaRadLayer(1) < kProfLayer) THEN
      iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
      iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
      iiDiv = 0
    ELSE
      ! essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1
 1010 CONTINUE
      IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
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
    DO iLay = iCldBotkCarta,iCldTopkCarta
      !!! this is for d/dq, d/dT
      iaCldLayer(kProfLayer-iLay+1) = 1
      !!!this is for d/d(DME), d/d(IWP)
      iaCldLayerIWPDME(kProfLayer-iLay+1 - (kProfLayer-iNumLayer)) = 1
    END DO
     
    IF (kSolar >= 0) THEN
      iSolarRadOrJac = +1
      CALL SolarScatterIntensity_Uplook( &
        iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer, &
        iNumLayer,iaRadLayer, &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot, &
        iTag,iSolarRadOrJac,raaSolarScatter1Lay)
    END IF

! cccccccccccccccccc set these all important variables ****************

    iIOUN = kStdJacob

    iLowest = iaaRadLayer(iAtm,1)            !!! this is for DOWN LOOk instr
    iLowest = iaaRadLayer(iAtm,iNumLayer)    !!!modify for UPLOOK instr
    iLowest = MOD(iLowest,kProfLayer)

    IF (kJacobOutPut >= 1) THEN
      kThermal = -1
      DO iG=1,kMaxPtsJac
        radBackgndThermdT(iG) = 0.0
        radSolardT(iG) = 0.0
      END DO
      CALL Find_BT_rad(raInten,radBTdr,raFreq,radBackgndThermdT,radSolardT)
    END IF
     
    write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
    CALL DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq, &
      iAtm,iNumLayer,iaaRadLayer, rSatAngle,raLayAngles,raSun,raaExtTemp, &
      raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd, &
      iProfileLayers,raPressLevels)
    write(kStdWarn,*)'initializing Jacobian loops ...'
    CALL Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,raLayAngles, &
      rTSpace,rTSurface,raUseEmissivity,raSurface,raSun,raThermal, &
      raaOneMinusTau,raaTau,raaLay2Gnd,raaRad,raaGeneral)

    DO iLay = 1,iNumLayer
      iL = iaRadlayer(iLay)
      IF (iLay == iNumLayer) THEN
        raLMm1_toGnd = 1.0
      ELSE
        raLMm1_toGnd = raaLay2Gnd(:,iLay+1)
      END IF
    END DO

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2).OR. (iWhichJac == 20)) THEN
      iJ = 0
      DO iG=1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iJ = iJ + 1
          iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          IF (iPrintAllPCLSAMJacs > 0) THEN
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
          END IF
          DO iLay = iNumLayer,1,-1
            rWeight = raaMix(iaaRadLayer(iAtm,iLay),iG)
            IF (iLay == iNumLayer) THEN
                rWeight = rWeight*rFracBot
            ELSEIF (iLay == 1) THEN
                rWeight = rWeight*rFracTop
            END IF
            write(kStdWarn,*)'gas d/dq gas# layer#',iG,iLay, &
            iaaRadLayer(iAtm,iLay)
            CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag, &
              iDoAdd,iErr)
            IF (iDoAdd > 0) THEN
              CALL JacobGasAmtFM1UP(raFreq,raSun,raaRad,iGasJacList,iLay, &
                        iNumGases, &
                        iaaRadLayer,iAtm,iNumLayer, &
                        raaOneMinusTau,raaTau,raaaAllDQ,raResults, &
                        raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
              IF (kSolar >= 0) THEN
                CALL SolarScatterGasJacobianUp( &
                            iNumlayer, &
                            iTag,iLay,iaRadLayer(iLay),iGasJacList,raFreq,raSunTOA, &
                            raaLay2Gnd,raLayAngles,raSunAngles, &
                            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
                            raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
                            raaSolarScatter1Lay,raaSolarScatterCumul, &
                            raResults)
              END IF
              iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
              CALL scale_raResults(raResults,rFracx)				
              CALL doJacobOutput(iLowest,raFreq,raResults, &
                    radBTdr,raaAmt,raInten,iaGases(iG),iWhichLayer,iGasPosn)
                  raaaAllJacQOut(iJ,:,iLay) = raResults
            ELSE
              raResults = 0.0
            END IF
            IF (iPrintAllPCLSAMJacs > 0) THEN
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END IF
          END DO
        END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      iJ = 0
      DO iG=1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iJ = iJ + 1
          iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          IF (iPrintAllPCLSAMJacs > 0) THEN
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
              raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = iNumLayer,1,-1
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              raaaAllJacQOut(iJ,:,iLay) = raResults
            END DO
          END IF
        END IF
      END DO
    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == 20) .OR. (iWhichJac == -2)) THEN
      DO iG = iNumGases+1,iNumGases+2
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iJ = iJ + 1
          iGasPosn = -1   !!!for JacobOutput
          IF (iG == iNumGases+1) THEN
            iIWPorDME = +1
          ELSEIF (iG == iNumGases+2) THEN
            iIWPorDME = -1
          END IF
          IF (iPrintAllPCLSAMJacs > 0) THEN
            print *,'wrtout head',iG,iaGasesTemp(iG),iNumLayer
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
              raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
          END IF
          DO iLay=iNumLayer,1,-1
            rWeight = 1.0
            IF (iG == iNumGases+1) THEN
              write(kStdWarn,*)'IWP d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
            ELSEIF (iG == iNumGases+2) THEN
              write(kStdWarn,*)'DME d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
            END IF
            IF (iaCldLayer(iLay) <= 0) THEN
              ! no cloud here, so indpt of cloud params!!!
              raResults = 0.0
              raaaAllJacQOut(iJ,:,iLay) = raResults
            ELSE
              CALL JacobCloudAmtUp(raFreq,raaRad,raaRadDT, &
                iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
                raaOneMinusTau,raaTau, &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
                raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
                raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME, &
                raaLay2Gnd,raResults,raThermal, &
                rSatAngle,raLayAngles, &
                raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
                iOffSet)
              IF (kSolar >= 0) THEN
                CALL SolarScatterCldAmtJacobianUp( &
                    iNumlayer, &
                    iTag,iLay,iaRadLayer(iLay),iIWPorDME,raFreq,raSunTOA, &
                    raaLay2Gnd,raLayAngles,raSunAngles, &
                    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
                    raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
                    raaSolarScatter1Lay,raaSolarScatterCumul, &
                    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
                    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
                    raResults)
              END IF
              CALL scale_raResults(raResults,rFracx)				
                raaaAllJacQOut(iJ,:,iLay) = raResults
            END IF
            iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
            CALL doJacobOutput(iLowest,raFreq,raResults, &
              radBTdr,raaAmt,raInten,iaGasesTemp(iG),iWhichLayer,iGasPosn)
            IF (iPrintAllPCLSAMJacs > 0) THEN
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END IF
          END DO
        END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iG=iNumGases+1,iNumGases+2
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iJ = iJ + 1
          iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          IF (iPrintAllPCLSAMJacs > 0) THEN
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
              raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = iNumLayer,1,-1
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              raaaAllJacQOut(iJ,:,iLay) = raResults
            END DO
          END IF
        END IF
      END DO
    END IF

! then do the temperatures d/dT
    IF (iPrintAllPCLSAMJacs > 0) THEN
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
        rDelta,iAtm,0,iNumLayer)
    END IF
    IF ((iWhichJac == -1) .OR. (iWhichJac == 30) .OR. (iWhichJac == -2) .OR. (iWhichJac == 32)) THEN
      DO iLay = iNumLayer,1,-1
        IF (iNatm > 1) THEN
          rWeight=0.0
          DO iG=1,iNumGases
            rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iLay),iG)
          END DO
          rWeight = rWeight/(iNumGases*1.0)
        ELSE
          rWeight=1.0
        END IF
        ! for each of the iNumLayer radiances, cumulatively add on all
        ! iNumGases contributions (this loop is done in JacobTemp)
          write(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay)
        CALL JacobTempFM1UP(raFreq,raSun,raaRad,raaRadDT,iLay, &
          iaaRadLayer,iAtm,iNumLayer, &
          raaOneMinusTau,raaTau,raaAllDT,raResults, &
          raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
        IF (kSolar >= 0) THEN
          CALL SolarScatterTempJacobianUp( &
              iNumlayer, &
              iTag,iLay,iaRadLayer(iLay),raFreq,raSunTOA, &
              raaLay2Gnd,raLayAngles,raSunAngles, &
              raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
              raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
              raaSolarScatter1Lay,raaSolarScatterCumul, &
              raResults)
        END IF
        iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
        CALL doJacobOutput(iLowest,raFreq,raResults,radBTdr,raaAmt,raInten,0,iWhichLayer,-1)
        CALL scale_raResults(raResults,rFracx)		    
        raaAllJacTOut(:,iLay) = raResults        
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      IF (iPrintAllPCLSAMJacs > 0) THEN
        DO iLay = iNumLayer,1,-1
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
      END IF
    END IF

! do the weighting functions
    IF (iPrintAllPCLSAMJacs > 0) THEN
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
        rDelta,iAtm,-10,iNumLayer)
    END IF
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. (iWhichJac == -2)) THEN
      DO iLay = iNumLayer,1,-1
        write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
        CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles, &
            iaaRadLayer,iAtm,raaLay2Gnd,raaExtTemp,raResults,rFracTop,rFracBot)
        ! does not make sense to multiply the weighting fcns with gas amounts etc
        !      iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
        !      CALL doJacobOutput(raFreq,raResults,radBTdr,raaAmt,raInten,0,
        !     $                               iWhichLayer)
        ! so just output the weighting functions
        CALL scale_raResults(raResults,rFracx)		
        raaAllWgtOut(:,iLay) = raResults
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
        raaAllWgtOut(:,iLay) = raResults
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iLay = iNumLayer,1,-1
        raaAllWgtOut(:,iLay) = 0.0
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    END IF
     
! computing Jacobians wrt surface parameters is meanigless .. output 0's
   IF (iPrintAllPCLSAMJacs > 0) THEN
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
        rDelta,iAtm,-20,4)
    END IF
    DO iLay = 1,4
      raaAllSurfOut(:,iLay) = 0.0
    END DO
    raResults = 0.0
    IF (iPrintAllPCLSAMJacs > 0) THEN
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    END IF

    RETURN
    end SUBROUTINE UpwardJacobian_ScatPCLSAM

!************************************************************************
! this subroutine adds on the solar contribution to the Cld Amt Jacobian
    SUBROUTINE SolarScatterCldAmtJacobianUp( &
    iNumlayer, &
    iTag,iLay,iLM,iIWPorDME,raFreq,raSunTOA, &
    raaLay2Gnd,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay,raaSolarScatterCumul, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iTag,iLay,iLM,iIWPorDME,iNumLayer
    INTEGER :: iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iFr,iDoSolar
    REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: raW(kMaxPts),raX(kMaxPts),raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMm1_toGnd(kMaxPts)

    CALL InitSomeSunScatUp( &
      iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer, &
      raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer, &
      raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iIWPorDME == +1) THEN
      !! this is IWP jacobian
      IF (iaCldLayer(iLay) == 1) THEN
        raW = exp(-raaExtTemp(:,iLM)/muSun)

        ! this is the change in optical depths wrt iwp
        raT = exp(-raaExtTemp(:,iLM)/muSat)/muSat - exp(-raaExtTemp(:,iLM)/muSun)/muSun
        raZ = -raaSSAlbTemp(:,iLM)*raT*raW

        raZ = raZ - raaSSAlbTemp(:,iLM)*raW/muSun*(exp(-raaExtTemp(:,iLM)/muSat)-exp(-raaExtTemp(:,iLM)/muSun))

        ! this is the change in ssalb wrt iwp
        raT = exp(-raaExtTemp(:,iLM)/muSat) - exp(-raaExtTemp(:,iLM)/muSun)
        raY = raaSSAlbJacobIWP(:,iLM)*raT*raW/raaExtJacobIWP(:,iLM)

        !! now add rZ and rY and divide by alpha_j
        !! to get (1/alpha_j) d(alpha_j))
        raZ = (raY + raZ)/(raaSSAlbTemp(:,iLM)*raT*raW)

        raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMm1_toGnd
      END IF

      raTemp = raTemp + raTemp1
      !sun          raResults = raTemp*raaExtJacobIWP(:,iLM)
      raResults = raResults + raTemp*raaExtJacobIWP(:,iLM)

    ELSEIF (iIWPorDME == -1) THEN
      !! this is DME jacobian
      IF (iaCldLayer(iLay) == 1) THEN
        raW = exp(-raaExtTemp(:,iLM)/muSun)

        ! this is the change in optical depths wrt dme
        raT = exp(-raaExtTemp(:,iLM)/muSat)/muSat - exp(-raaExtTemp(:,iLM)/muSun)/muSun
        raZ = -raaSSAlbTemp(:,iLM)*raT*raW

        raZ = raZ - raaSSAlbTemp(:,iLM)*raW/muSun* &
          (exp(-raaExtTemp(:,iLM)/muSat)-exp(-raaExtTemp(:,iLM)/muSun))

        ! this is the change in ssalb wrt dme
        raT = exp(-raaExtTemp(:,iLM)/muSat) - exp(-raaExtTemp(:,iLM)/muSun)
        raY = raaSSAlbJacobDME(:,iLM)*raT*raW/raaExtJacobDME(:,iLM)

        !! now add rZ and rY and divide by alpha_j
        !! to get (1/alpha_j) d(alpha_j))
        raZ = (raY + raZ)/(raaSSAlbTemp(:,iLM)*raT*raW)

        raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMm1_toGnd

        ! also add on the d(hg)/d(dme) contribution!!!!
        raZ = 1/rahg2_real(-muSun,-muSat,raaAsymTemp(:,iLM))*raaPhaseJacobASYM(:,iLM)
        raZ = raZ*raaAsymJacobDME(:,iLM)/raaExtJacobDME(:,iLM)
        raTemp1 = raTemp1 + raZ*raaSolarScatter1Lay(:,iLM)*raLMm1_toGnd
      END IF

      raTemp = raTemp + raTemp1
      !sun          raResults = raTemp*raaExtJacobDME(:,iLM)
      raResults = raResults + raTemp*raaExtJacobDME(:,iLM)
    END IF

    RETURN
    end SUBROUTINE SolarScatterCldAmtJacobianUp

!************************************************************************
! this subroutine adds on the solar contribution to the Temp Jacobian
    SUBROUTINE SolarScatterTempJacobianUp( &
    iNumlayer, &
    iTag,iLay,iLM,raFreq,raSunTOA, &
    raaLay2Gnd,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay,raaSolarScatterCumul, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iTag,iLay,iLM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
    INTEGER :: iNumLayer
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iFr,iDoSolar
    REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: raW(kMaxPts),raX(kMaxPts),raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMm1_toGnd(kMaxPts)

    CALL InitSomeSunScatUp( &
      iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer, &
      raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer, &
      raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iaCldLayer(iLay) == 1) THEN
      raW = exp(-raaExtTemp(:,iLM)/muSun)

      ! this is the change in optical depths wrt T
      raT = exp(-raaExtTemp(:,iLM)/muSat)/muSat - exp(-raaExtTemp(:,iLM)/muSun)/muSun
      raZ = -raaSSAlbTemp(:,iLM)*raT*raW

      raZ = raZ - raaSSAlbTemp(:,iLM)*raW/muSun*(exp(-raaExtTemp(:,iLM)/muSat)-exp(-raaExtTemp(:,iLM)/muSun))

      ! this is the change in ssalb wrt gas temperature
      raT = exp(-raaExtTemp(:,iLM)/muSat) - exp(-raaExtTemp(:,iLM)/muSun)
      raY = (1 + raaAsymTemp(:,iLM))/2.0
      raY = -raaSSAlbTemp(:,iLM)*(1-raaSSAlbTemp(:,iLM)*raY)
      raY = raT*raW*raY/raaExtTemp(:,iLM)

      !! now add rZ and rY and divide by alpha_j
      !! to get (1/alpha_j) d(alpha_j))
      raZ = (raY + raZ)/(raaSSAlbTemp(:,iLM)*raT*raW)

      raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMm1_toGnd
    END IF

! add on to the raResults tally
    raTemp = raTemp + raTemp1
    !sun        raResults = raTemp*raaAllDT(:,iLM)
    raResults = raResults + raTemp*raaAllDT(:,iLM)

    RETURN
    end SUBROUTINE SolarScatterTempJacobianUp

!************************************************************************
! this subroutine adds on the solar contribution to the Amt Jacobian
    SUBROUTINE SolarScatterGasJacobianUp( &
    iNumlayer, &
    iTag,iLay,iLM,iG,raFreq,raSunTOA, &
    raaLay2Gnd,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay,raaSolarScatterCumul, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iTag,iLay,iLM,iG,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
    INTEGER :: iNumLayer
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iFr,iDoSolar
    REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: raW(kMaxPts),raX(kMaxPts),raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMm1_toGnd(kMaxPts)

    CALL InitSomeSunScatUp( &
      iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer, &
      raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer, &
      raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iaCldLayer(iLay) == 1) THEN
      raW = exp(-raaExtTemp(:,iLM)/muSun)

      ! this is the change in optical depths wrt T
      raT = exp(-raaExtTemp(:,iLM)/muSat)/muSat - exp(-raaExtTemp(:,iLM)/muSun)/muSun
      raZ = -raaSSAlbTemp(:,iLM)*raT*raW

      raZ = raZ - raaSSAlbTemp(:,iLM)*raW/muSun*(exp(-raaExtTemp(:,iLM)/muSat)-exp(-raaExtTemp(:,iLM)/muSun))

      ! this is the change in ssalb wrt gas temperature
      raT = exp(-raaExtTemp(:,iLM)/muSat) - exp(-raaExtTemp(:,iLM)/muSun)
      raY = (1 + raaAsymTemp(:,iLM))/2.0
      raY = -raaSSAlbTemp(:,iLM)*(1-raaSSAlbTemp(:,iLM)*raY)
      raY = raT*raW*raY/raaExtTemp(:,iLM)

      !! now add rZ and rY and divide by alpha_j
      !! to get (1/alpha_j) d(alpha_j))
      raZ = (raY + raZ)/(raaSSAlbTemp(:,iLM)*raT*raW)

      raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMm1_toGnd
    END IF

! add on to the raResults tally
    raTemp = raTemp + raTemp1
    !sun        raResults = raTemp*raaaAllDQ(iG,:,iLM)
    raResults = raResults + raTemp*raaaAllDQ(iG,:,iLM)

    RETURN
    end SUBROUTINE SolarScatterGasJacobianUp

!************************************************************************
! this subroutine does some initializations common to the Cloud Jac routines
! in particular, it initializes the cumulative jacobian for (cloud) layers
! beneath this layer (iLay) in question
    SUBROUTINE InitSomeSunScatUp( &
    iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer, &
    raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer, &
    raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iLM,iLay,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer),iNumLayer
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
! output vars
    REAL :: muSun,muSat
    REAL :: raTemp(kMaxPts),raTemp1(kMaxPts),raLMm1_toGnd(kMaxPts)

! local vars
    INTEGER :: iFr,iJ,iLJ
    REAL :: rX

! iLM = iaRadLayer(iLay)
    raTemp1 = 0.0

! stuff on Jan 5, 2006
    raTemp = 0.0
    DO iJ = 1,iLay-1
      iLJ = iaRadLayer(iJ)
      IF (iaCldLayer(iJ) == 1) THEN
        muSat = cos(raLayAngles(iLJ) * kPi/180)
        muSun = cos(raSunAngles(iLJ) * kPi/180)
        rX = (-1/muSat)*(1+muSat/muSun)
        rX = (-1/muSun)
        raTemp = raTemp + raaSolarScatter1Lay(:,iLJ)*raaLay2Gnd(:,iJ+1)*rX
      END IF
    END DO
    DO iJ = iLay+1,iNumLayer
      iLJ = iaRadLayer(iJ)
      IF (iaCldLayer(iJ) == 1) THEN
        muSat = cos(raLayAngles(iLJ) * kPi/180)
        muSun = cos(raSunAngles(iLJ) * kPi/180)
        rX = (-1/muSat)*(1+muSat/muSun)
        rX = (-1/muSat)
        raTemp = raTemp + raaSolarScatter1Lay(:,iLJ)*raaLay2Gnd(:,iJ+1)*rX
      END IF
    END DO

! output   these vars for the rest of the routine
    muSat = cos(raLayAngles(iLM) * kPi/180)
    muSun = cos(raSunAngles(iLM) * kPi/180)

    IF (iLM == iNumLayer) THEN
     raLMm1_toGnd = 1.0
    ELSE
      raLMm1_toGnd = raaLay2Gnd(:,iLay+1)
    END IF

    RETURN
    end SUBROUTINE InitSomeSunScatUp

!************************************************************************
! this does the jacobian wrt IWP or DME
    SUBROUTINE JacobCloudAmtUp(raFreq,raaRad,raaRadDT, &
    iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
    raaOneMinusTau,raaTau, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME, &
    raaLay2Gnd,raResults,raThermal, &
    rSatAngle,raLayAngles, &
    raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
    iOffSet)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! daaDT,daaDQ are the d/dq,d/dT matrices
! raaLay2Gnd   is the layer-to-gnd abs coeff matrix
! raaRad has the Planck radiances
! raaRadDT has the d/DT (Planck radiances)
! raaOneMinusTau has 1-tau(satellite angle)
! raaOneMinusTauTh has 1-tau(thermal angle)
! iG has the gas number (1 .. iNumGasesTemp)
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raThermal are the downwelling Solar,thermal contributions
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
! raaGeneral,raaGeneralTh have the general results (looping over layers)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle
    REAL :: raLayAngles(kProfLayer)
    REAL :: raThermal(kMaxPts)
    INTEGER :: iNumGasesTemp,iAtm
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac),raFreq(kMaxPts)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iG,iLay,iNumLayer,iOffSet
! output vars
    REAL :: raResults(kMaxPtsJac)

! this is for the scattering parts
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)  !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)  !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac)  !phase fnc jacobs wrt g
    INTEGER :: iIWPorDME

! local variables
    INTEGER :: iFr,iJ1,iM1
    REAL :: raTemp(kMaxPtsJac),muSat
    REAL :: raResultsTh(kMaxPtsJac)

! figure out which of 1..100 this current radiating layer corresponds to
! bleh
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
    muSat = 1.0/cos(rSatAngle*kPi/180.0)
    muSat = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
    raResults = raaGeneral(:,iLay)

! note that      iLay + iOffset === iaRadlayer(iLay)
! set the constant factor we have to multiply results with
    IF (iIWPorDME == 1) THEN
      DO iFr=1,kMaxPts
        raTemp = raaExtJacobIWP(:,iLay+iOffSet)
        raTemp = raaExtJacobIWP(:,iM1)
      END DO
    ELSEIF (iIWPorDME == -1) THEN
      DO iFr=1,kMaxPts
        raTemp = raaExtJacobDME(:,iLay+iOffSet)
        raTemp = raaExtJacobDME(:,iM1)
      END DO
    END IF
    CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
    IF (iLay < iNumLayer) THEN
      ! this is not the topmost layer
      iJ1 = iLay
      DO iFr=1,kMaxPts
        raResults = raResults+ &
        raTemp*raaRad(:,iJ1)*raaLay2Gnd(:,iJ1)
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      ! do the topmost layer correctly
      iJ1 = iLay
      DO iFr=1,kMaxPts
        raResults = raResults+ &
        raTemp*raaTau(:,iJ1)*raaRad(:,iJ1)
      END DO
    END IF

! now multiply results by the 1/cos(viewing angle) factor
    IF (abs(muSat-1.0000000) >= 1.0E-5) THEN
      raResults = raResults*muSat
    END IF

    RETURN
    end SUBROUTINE JacobCloudAmtUp

!************************************************************************
END MODULE jac_pclsam_up
