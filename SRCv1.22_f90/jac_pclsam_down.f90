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
! iaaRadLayer(kMaxAtm,kProfLayer),raaAbs(kMaxPts,kMixFilRows)
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

! put in a lot of rEps in the SUBROUTINE solarscatter/gas/temp/iwpdme as things
! could go nuts if raaIWP = 0 in between cloud layers

MODULE jac_pclsam_down

USE basic_common
USE ttorad_common
use jac_up
use jac_down
use jac_main
use singlescatter
use clear_scatter_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
!****** THESE ARE THE SCATTERING JACOBIANS FOR THE DOWN LOOK INSTR ******
!************************************************************************
! this subroutine does the Jacobians for downward looking instrument
    SUBROUTINE DownwardJacobian_ScatPCLSAM(raFreq, &
    iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp,ctype2,rFracx, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX, &
    iNpMix,iTag,iActualTag, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA, iPrintAllPCLSAMJacs, &
    iNLTEStart,raaPlanckCoeff, &
    raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    CHARACTER(160) :: caJacobFile
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
    INTEGER :: iNpMix,iCLoudySky
    REAL :: iwpMAX(MAXNZ)
    INTEGER :: iPrintAllPCLSAMJacs

! local variables
    INTEGER :: iNumGasesTemp,iaGasesTemp(kMaxGas)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
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
    INTEGER :: iG,iM,iJ,iIOUN,iLowest
    INTEGER :: iGasJacList,iGasPosn
! for cloud stuff!
    REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
    INTEGER :: iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iiDiv,iaRadLayer(kProfLayer),iLay,iIWPorDME,iL
    INTEGER :: iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
    REAL :: r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck,muSun,rSunAngle
    INTEGER :: iSolarRadOrJac
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

    INTEGER :: iDefault,iWhichJac,iFr
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

    rSunAngle = raSunAngles(MP2Lay(iaaRadLayer(1,iAtm)))
    ! change to radians
    rSunAngle=(rSunAngle*kPi/180.0)
    muSun = cos(rSunAngle)

    rSunTemp  = kSunTemp
    rOmegaSun = kOmegaSun
    iDoSolar = kSolar
    IF (iDoSolar == 0) THEN
      ! use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      ! compute the Plank radiation from the sun
      raSunTOA = ttorad(raFreq,rSunTemp)*rOmegaSun*muSun
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      ! read in data from file
      CALL ReadSolarData(raFreq,raSunTOA,iTag)
      raSunTOA = raSunTOA*rOmegaSun*muSun
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

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
     iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
     iL = iaRadLayer(iLay)
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
    IF (iaRadLayer(1) < kProfLayer) THEN
      iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
      iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
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
      iLocalCldTop = iCldTopkCarta - iiDiv + 1
      iLocalCldBot = iCldBotkCarta - iiDiv + 1
      iiDiv = iLay
    END IF

    iOffSet = (kProfLayer-iNumLayer)
    !DO iLay = iCldBotkCarta-1,iCldTopkCarta-1
    DO iLay = iCldBotkCarta,iCldTopkCarta
      !!! this is for d/dq, d/dT
      iaCldLayer(iLay) = 1
      !!!this is for d/d(DME), d/d(IWP)
      iaCldLayerIWPDME(iLay - (kProfLayer-iNumLayer)) = 1
    END DO

    IF (kSolar >= 0) THEN
      iSolarRadOrJac = +1
      CALL SolarScatterIntensity_Downlook( &
        iDoSolar,raFreq,iaCldLayer, &
        raSunAngles,raLayAngles,0.0,0.0, &
        iNumLayer,iaRadLayer, &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot, &
        iTag,iSolarRadOrJac,raaSolarScatter1Lay)

    END IF

! cccccccccccccccccc set these all important variables ****************

! initialise the layer-to-space matrix
    CALL AtmosLayer2Space(raaLay2Sp, &
      rSatAngle,raaExtTemp,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

    write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
    CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq, &
      iAtm,iNumLayer,iaaRadLayer, &
      rSatAngle,raLayAngles,raaExtTemp, &
      raaRad,raaRadDT,raaOneMinusTau,raaTau, &
      raaLay2Gnd,iProfileLayers,raPressLevels)

    write(kStdWarn,*)'initializing Jacobian loops ...'
    CALL Loop_LookDown(iaaRadLayer, &
      iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl, &
      raUseEmissivity,raSurface,raSun,raSunAngles,raThermal, &
      raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      write(kStdWarn,*)'initializing Jacobian thermal loops ...'
      CALL Loop_thermaldown(raaRad,rSatAngle,raLayAngles, &
        iProfileLayers,raPressLevels, &
        iaaRadLayer,iAtm,iNumLayer,raaExtTemp, &
        raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)

    END IF

    IF (kJacobOutPut >= 1) THEN
      ! have to set the backgnd thermal, solar radiances so that we can do
      ! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
      IF (kThermal >= 0) THEN
        DO iG=1,kMaxPts
          radBackgndThermdT(iG) = raThermal(iG)
        END DO
      END IF

      IF (kSolar >= 0) THEN
        ! compute the Solar contribution
        iM = 1
        ! remember that raSun is that at the gnd -- we have to propagate this back
        ! up to the top of the atmosphere
        ! note that here we are calculating the SOLAR contribution
        DO iG = 1,kMaxPts
          radSolardT(iG) = raUseEmissivity(iG)*raaLay2Sp(iG,iM)*raSun(iG)
        END DO
      END IF

      CALL Find_BT_rad(raInten,radBTdr,raFreq,radBackgndThermdT,radSolardT)
    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
      iJ = 0
      DO iG = 1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iJ = iJ + 1
          iGasPosn = WhichGasPosn(iaGasesTemp(iG),iaGasesTemp,iNumGasesTemp)
          IF (iPrintAllPCLSAMJacs > 0) THEN
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                           raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
          END IF
          DO iM=1,iNumLayer
            rWeight = raaMix(iaaRadLayer(iAtm,iM),iG)
            IF (iM == 1) THEN
              rWeight = rWeight*rFracBot
            ELSEIF (iM == iNumLayer) THEN
              rWeight = rWeight*rFracTop
            END IF
            write(kStdWarn,*)'gas d/dq gas# lay#',iG,iM,iaaRadLayer(iAtm,iM)
            !! see if this gas does exist for this chunk
            CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag, &
                             iDoAdd,iErr)
            IF (iDoAdd > 0) THEN
              CALL JacobGasAmtFM1(raFreq,raaRad,raaRadDT,iGasJacList, &
                        iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
                        raaOneMinusTau,raaTau,raaaAllDQ, &
                        raaLay2Sp,raResults,raThermal,raaLay2Gnd, &
                        rSatAngle,raLayAngles, &
                        raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
                        rWeight)
              IF (kSolar >= 0) THEN
                CALL SolarScatterGasJacobian( &
                            iTag,iM,iaRadLayer(iM),iGasJacList,raFreq,raSunTOA, &
                            raaLay2Sp,raLayAngles,raSunAngles, &
                            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
                            raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
                            raaSolarScatter1Lay, &
                            raResults)
              END IF
	      CALL scale_raResults(raResults,rFracx)

!!!!!!sergio debug
!!!!              print *,rFracx,radBTdr(1),raInten(1),raResults(1)

              CALL doJacobOutput(iLowest,raFreq, &
                                 raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
              raaaAllJacQOut(iJ,:,iM) = raResults
            ELSE
              raResults = 0.0
            END IF
            IF (iPrintAllPCLSAMJacs > 0) THEN
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END IF
          END DO       !!! iM=1,iNumLayer
!!!!! sergio 
!!    stop

        END IF       !!! if iGasList > 0
      END DO         !!! iG=1,iNumGases
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      iJ = 0
      DO iG=1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          iJ = iJ + 1
          IF (iPrintAllPCLSAMJacs > 0) THEN
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = 1,iNumLayer
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              raaaAllJacQOut(iJ,:,iM) = raResults
            END DO

          END IF     !! IF (iPrintAllPCLSAMJacs > 0) THEN
        END IF       !! IF (iGasJacList > 0)
      END DO         !! iG = 1,iNumGases
    END IF           !! IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
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
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
          END IF
          DO iM=1,iNumLayer
            rWeight = 1.0
            IF (iG == iNumGases+1) THEN
              write(kStdWarn,*)'IWP d/dq lay#',iM,iaaRadLayer(iAtm,iM)
            ELSEIF (iG == iNumGases+2) THEN
              write(kStdWarn,*)'DME d/dq lay#',iM,iaaRadLayer(iAtm,iM)
            END IF
            ! this is more correct as it only puts out d/d(IWP) and d/d(DME) where there
            ! actually is cloud
            !            IF (iaCldLayerIWPDME(iM) .EQ. 0 .OR.
            !     $          iwpMAX(iaaRadLayer(iAtm,iM)) .LE. 1.0e-10) THEN
            ! very fun thing happens if you let this line compile ... and look at d/d(IWP)
            ! it mimics T(z)!!!!!
            IF (iaCldLayerIWPDME(iM) == 0) THEN
              ! no cloud here, so indpt of cloud params!!!
              raResults = 0.0
              raaaAllJacQOut(iJ,:,iM) = raResults
            ELSE
              CALL JacobCloudAmtFM1(raFreq,raaRad,raaRadDT, &
                        iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
                        raaOneMinusTau,raaTau, &
                        raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
                        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
                        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME, &
                        raaLay2Sp,raResults,raThermal,raaLay2Gnd, &
                        rSatAngle,raLayAngles, &
                        raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
                        iOffSet)
              IF (kSolar >= 0) THEN
                CALL SolarScatterCldAmtJacobian( &
                            iTag,iM,iaRadLayer(iM),iIWPorDME,raFreq,raSunTOA, &
                            raaLay2Sp,raLayAngles,raSunAngles, &
                            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
                            raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
                            raaSolarScatter1Lay,iwpMAX, &
                            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
                            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
                            raResults)
              END IF
            CALL scale_raResults(raResults,rFracx)			
            raaaAllJacQOut(iJ,:,iM) = raResults
            END IF  !! IF (iaCldLayerIWPDME(iM) == 0)
            CALL doJacobOutput(iLowest,raFreq, &
                      raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
            IF (iPrintAllPCLSAMJacs > 0) THEN
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END IF
          END DO
        END IF
      END DO         !!! iG=1,iNumGases
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
!            print *,'wrtout head',iG,iaGasesTemp(iG),iNumLayer
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                           raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iLay = 1,iNumLayer
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
              raaaAllJacQOut(iJ,:,iLay) = raResults
            END DO
          END IF     !!! IF (iPrintAllPCLSAMJacs > 0) THEN
        END IF
      END DO         !!! iG=1,iNumGases
    END IF

! then do the temperatures d/dT
    IF (iPrintAllPCLSAMJacs > 0) THEN
!       print *,'wrtout head T',iNumLayer
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
        rDelta,iAtm,0,iNumLayer)
    END IF
    IF ((iWhichJac == -1) .OR. (iWhichJac == 30) .OR. (iWhichJac == -2) .OR. (iWhichJac == 32)) THEN
      DO iM=1,iNumLayer
        ! for each of the iNumLayer radiances, cumulatively add on all
        ! iNumGases contributions (this loop is done in JacobTemp)
        write(kStdWarn,*)'temp d/dT layer# = ',iM,iaaRadLayer(iAtm,iM)
        ! don't have to worry about gases 210,202 as their "weight" for d/dT = 0.0
        IF (iNatm > 1) THEN
          rWeight = 0.0
          DO iG=1,iNumGases
            rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iM),iG)
          END DO
          rWeight = rWeight/(iNumGases*1.0)
        ELSE
          rWeight = 1.0
        END IF
        CALL JacobTempFM1(raFreq,raaRad,raaRadDT,iM,iNumGases, &
            iaaRadLayer,iAtm,iNumLayer, &
            raUseEmissivity, &
            raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults, &
            raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral, &
            raaGeneralTh,raaOneMinusTauTh,rWeight)
        IF (kSolar >= 0) THEN
          CALL SolarScatterTemperatureJacobian( &
                iTag,iM,iaRadLayer(iM),raFreq,raSunTOA, &
                raaLay2Sp,raLayAngles,raSunAngles, &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
                raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
                raaSolarScatter1Lay, &
                raResults)
        END IF
        CALL scale_raResults(raResults,rFracx)	    
        CALL doJacobOutput(iLowest,raFreq,raResults, &
              radBTdr,raaAmt,raInten,0,iM,-1)
         raaAllJacTOut(:,iM) = raResults        
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      IF (iPrintAllPCLSAMJacs > 0) THEN
        raResults = 0.0
        DO iLay = 1,iNumLayer
           raaAllJacTOut(:,iLay) = raResults        
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
      END IF
    END IF

! do the weighting functions
    IF (iPrintAllPCLSAMJacs > 0) THEN
!       print *,'wrtout head WGT',iNumLayer
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
        rDelta,iAtm,-10,iNumLayer)
    END IF
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. (iWhichJac == -2)) THEN
      DO iM=1,iNumLayer
        write(kStdWarn,*)'wgt fcn # = ',iM,iaaRadLayer(iAtm,iM)
        CALL wgtfcndown(iM,iNumLayer,rSatAngle,raLayAngles, &
            iaaRadLayer,iAtm,raaLay2Sp,raaExtTemp,raResults,rFracTop,rFracBot, &
            iNLTEStart,raaPlanckCoeff)
        CALL scale_raResults(raResults,rFracx)	
        raaAllWgtOut(:,iM) = raResults
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iM = 1,iNumLayer
        raaAllWgtOut(:,iM) = raResults
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
    IF (iPrintAllPCLSAMJacs > 0) THEN
!      print *,'wrtout head Surf',4
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
            rDelta,iAtm,-20,4)
    END IF
    IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR. (iWhichJac == -2)) THEN
      iM=1
      ! computing Jacobians wrt surface parameters are meaningful
      CALL JacobSurfaceTemp(raFreq,iM,rTSurface,raUseEmissivity,raaLay2Sp,raResults)	
      CALL doJacobOutput(iLowest,raFreq,raResults,radBTdr,raaAmt,raInten,-1,iM,-1)
      CALL scale_raResults(raResults,rFracx)	
      IF (iPrintAllPCLSAMJacs > 0) THEN
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END IF
      raaAllSurfOut(:,1) = raResults

      CALL JacobSurfaceEmis(iM,raSurface,raThermal,raaLay2Sp,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radBTdr,raaAmt,raInten,-2,iM,-1)
      CALL scale_raResults(raResults,rFracx)		  
      IF (iPrintAllPCLSAMJacs > 0) THEN
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END IF
      raaAllSurfOut(:,2) = raResults

      CALL JacobBackgndThermal(iM,raaLay2Sp,raThermal,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radBackgndThermdT,raaAmt,raInten,-3,iM,-1)
      CALL scale_raResults(raResults,rFracx)		  
      IF (iPrintAllPCLSAMJacs > 0) THEN
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END IF
      raaAllSurfOut(:,3) = raResults

      CALL JacobSolar(iM,raaLay2Sp,raSun,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radSolardT,raaAmt,raInten,-4,iM,-1)
      CALL scale_raResults(raResults,rFracx)		
      IF (iPrintAllPCLSAMJacs > 0) THEN
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END IF
      raaAllSurfOut(:,4) = raResults
	
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iLay = 1,4
        raaAllSurfOut(:,iLay) = raResults
        IF (iPrintAllPCLSAMJacs > 0) THEN
          CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END IF
      END DO
    END IF
    
    RETURN
    end SUBROUTINE DownwardJacobian_ScatPCLSAM

!************************************************************************
! this subroutine adds on the solar contribution to the Cld Amt Jacobian
    SUBROUTINE SolarScatterCldAmtJacobian( &
    iTag,iM,iLM,iIWPorDME,raFreq,raSunTOA, &
    raaLay2Sp,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay,iwpMAX, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! input vars
    INTEGER :: iTag,iM,iLM,iIWPorDME
    INTEGER :: iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts),iwpMAX(MAXNZ)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iLay,iFr,iDoSolar
    REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: rFac,rX,raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMp1_toSpace(kMaxPts)
    REAL :: rEps   !! so that if w == 0, rz is finite

    rEps = 1.0e-10

    CALL InitSomeSunScat( &
      iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer, &
      raaSolarScatter1Lay,iaRadLayer, &
      raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iIWPorDME == +1) THEN
      !! this is IWP jacobian
      IF (iaCldLayer(iLM) == 1) THEN
        muSunSat = (muSun * muSat)/(muSun + muSat)
        ! this is the transmission
        raT = exp(-raaExtTemp(:,iLM)/muSunSat)

        ! this is the change in ssalb wrt iwp
        raY = raaSSAlbJacobIWP(:,iLM)*(1.0-raT)/(raaExtJacobIWP(:,iLM)+rEps)

        ! this is the change in optical depths wrt iwp
        raZ = raaSSAlbTemp(:,iLM)*((1/muSunSat+1/muSun)*raT - 1/muSun)
        raZ = raZ + rEps

        !! now add rZ and rY and divide by alpha_j
        !! to get (1/alpha_j) d(alpha_j))
        raZ = (raY + raZ)/((raaSSAlbTemp(:,iLM)+rEps)*(1-raT+rEps))

        raT = exp(-raaExtTemp(:,iLM)/muSun)
        raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMp1_toSpace

      END IF

      raTemp = raTemp + raTemp1
      !sun          raResults = raTemp*raaExtJacobIWP(:,iLM)
      raResults = raResults + raTemp*raaExtJacobIWP(:,iLM)

    ELSEIF (iIWPorDME == -1) THEN
      !! this is DME jacobian
      IF (iaCldLayer(iLM) == 1) THEN
        muSunSat = (muSun * muSat)/(muSun + muSat)
        ! this is the transmission
        raT = exp(-raaExtTemp(:,iLM)/muSunSat)

        ! this is the change in ssalb wrt dme
        raY = raaSSAlbJacobDME(:,iLM)*(1.0-raT)/(raaExtJacobDME(:,iLM)+rEps)

        ! this is the change in optical depths wrt dme
        raZ = raaSSAlbTemp(:,iLM)*((1/muSunSat+1/muSun)*raT - 1/muSun)
        raZ = raZ + rEps

        !! now add rZ and rY and divide by alpha_j
        !! to get (1/alpha_j) d(alpha_j))
        raZ = (raY + raZ)/((raaSSAlbTemp(:,iLM)+rEps)*(1-raT+rEps))
              
        raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMp1_toSpace

        ! also add on the d(hg)/d(dme) contribution!!!!
        raZ = 1/rahg2_real(-muSun,muSat,raaAsymTemp(:,iLM))*raaPhaseJacobASYM(:,iLM)
        raZ = raZ*raaAsymJacobDME(:,iLM)/(raaExtJacobDME(:,iLM)+rEps)
        raTemp1 = raTemp1 + raZ*raaSolarScatter1Lay(:,iLM)*raLMp1_toSpace
      END IF

      raTemp = raTemp + raTemp1
      !sun          raResults = raTemp*raaExtJacobDME(:,iLM)
      raResults = raResults + raTemp*raaExtJacobDME(:,iLM)
    END IF

    RETURN
    end SUBROUTINE SolarScatterCldAmtJacobian

!************************************************************************
! this subroutine adds on the solar contribution to the Temp Jacobian
    SUBROUTINE SolarScatterTemperatureJacobian( &
    iTag,iM,iLM,raFreq,raSunTOA, &
    raaLay2Sp,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input vars
    INTEGER :: iTag,iM,iLM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iLay,iFr,iDoSolar
    REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: rFac,raX(kMaxPts),raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMp1_toSpace(kMaxPts)
    REAL :: rEps   !! so that if w == 0, rz is finite

    rEps = 1.0-10
          
    CALL InitSomeSunScat( &
      iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer, &
      raaSolarScatter1Lay,iaRadLayer, &
      raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iaCldLayer(iLM) == 1) THEN
      muSunSat = (muSun * muSat)/(muSun + muSat)
      ! this is the transmission
      raT = exp(-raaExtTemp(:,iLM)/muSunSat)

      ! this is the change in ssalb wrt gas temperature
      raY = (1 + raaAsymTemp(:,iLM))/2.0
      raY = -raaSSAlbTemp(:,iLM)*(1-raaSSAlbTemp(:,iLM)*raY)
      raY = (1.0-raT)*raY/raaExtTemp(:,iLM)

      ! this is the change in optical depths wrt temperature
      raZ = raaSSAlbTemp(:,iLM)*((1/muSunSat+1/muSun)*raT - 1/muSun)
      raZ = raZ + rEps

      !! now add rZ and rY and divide by alpha_j
      !! to get (1/alpha_j) d(alpha_j))
      raZ = (raY + raZ)/((raaSSAlbTemp(:,iLM)+rEps)*(1-raT+rEps))

      raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMp1_toSpace

    END IF

! add on to the raResults tally
    raTemp = raTemp + raTemp1
    !sun        raResults = raTemp*raaAllDT(:,iLM)
    raResults = raResults + raTemp*raaAllDT(:,iLM)

    RETURN
    end SUBROUTINE SolarScatterTemperatureJacobian

!************************************************************************
! this subroutine adds on the solar contribution to the Amt Jacobian
    SUBROUTINE SolarScatterGasJacobian( &
    iTag,iM,iLM,iG,raFreq,raSunTOA, &
    raaLay2Sp,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, &
    raaSolarScatter1Lay, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input vars
    INTEGER :: iTag,iM,iLM,iG,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raFreq(kMaxPts),raSunTOA(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
! output vars
    REAL :: raResults(kMaxPtsJac)

! local vars
    INTEGER :: iL,iLay,iFr,iDoSolar
    REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
    REAL :: rFac,raX(kMaxPts),raY(kMaxPts),raZ(kMaxPts),raT(kMaxPts),raLMp1_toSpace(kMaxPts)
    REAL :: rEps   !! so that if w == 0, rz is finite

    rEps = 1.0e-10

    CALL InitSomeSunScat( &
      iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer, &
      raaSolarScatter1Lay,iaRadLayer, &
      raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
    IF (iaCldLayer(iLM) == 1) THEN
      muSunSat = (muSun * muSat)/(muSun + muSat)
      ! this is the transmission
      raT = exp(-raaExtTemp(:,iLM)/muSunSat)

      ! this is the change in ssalb wrt gas amount
      raY = (1 + raaAsymTemp(:,iLM))/2.0
      raY = -raaSSAlbTemp(:,iLM)*(1-raaSSAlbTemp(:,iLM)*raY)
      raY = (1.0-raT)*raY/raaExtTemp(:,iLM)

      ! this is the change in optical depths wrt gas amount
      raZ = raaSSAlbTemp(:,iLM)*((1/muSunSat+1/muSun)*raT - 1/muSun)
      raZ = raZ + rEps

      !! now add rZ and rY and divide by alpha_j to get
      !! (1/alpha_j) d(alpha_j))
      raZ = (raY + raZ)/((raaSSAlbTemp(:,iLM)+rEps)*(1-raT+rEps))

      raTemp1 = raZ*raaSolarScatter1Lay(:,iLM)*raLMp1_toSpace
             
    END IF

! add on to the raResults tally
    raTemp = raTemp + raTemp1
    !sun        raResults = raTemp*raaaAllDQ(iG,:,iLM)
    raResults = raResults + raTemp*raaaAllDQ(iG,:,iLM)

    RETURN
    end SUBROUTINE SolarScatterGasJacobian

!************************************************************************
! this subroutine does some initializations common to the Cloud Jac routines
! in particular, it initializes the cumulative jacobian for (cloud) layers
! beneath this layer (iM) in question
    SUBROUTINE InitSomeSunScat( &
    iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer, &
    raaSolarScatter1Lay,iaRadLayer, &
    raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input vars
    INTEGER :: iLM,iM,iaCldLayer(kProfLayer),iaRadLayer(kProfLayer)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
! output vars
    REAL :: muSun,muSat,rFac
    REAL :: raTemp(kMaxPts),raTemp1(kMaxPts),raLMp1_toSpace(kMaxPts)

! local vars
    INTEGER :: iFr,iJ,iLJ
    REAL :: rX

! iLM = iaRadLayer(iM)
    raTemp1 = 0.0

! stuff on Jan 5, 2006
    muSat = cos(raLayAngles(iLM) * kPi/180)
    muSun = cos(raSunAngles(iLM) * kPi/180)
    rFac = 1.0 + muSat/muSun

    rX = (-1/muSat)*(1+muSat/muSun)
!!! oh boy is this wrong? i think it is fine 1/12/06
    raTemp = 0.0
    DO iJ = 1,iM-1
      iLJ = iaRadLayer(iJ)
      IF (iaCldLayer(iLJ) == 1) THEN
        ! muSat = cos(raLayAngles(iLJ) * kPi/180)
        ! muSun = cos(raSunAngles(iLJ) * kPi/180)
        ! mX = (-1/muSat)*(1+muSat/muSun)
        raTemp = raTemp + raaSolarScatter1Lay(:,iLJ)*raaLay2Sp(:,iJ+1)*rX
      END IF
    END DO

! output these vars for the rest of the routine
    muSat = cos(raLayAngles(iLM) * kPi/180)
    muSun = cos(raSunAngles(iLM) * kPi/180)
    rFac = 1.0 + muSat/muSun

    IF (iLM == kProfLayer) THEN
      raLMp1_toSpace = 1.0
    ELSE
      raLMp1_toSpace = raaLay2Sp(:,iM+1)
    END IF

    RETURN
    end SUBROUTINE InitSomeSunScat

!************************************************************************
! this does the jacobian wrt IWP or DME
    SUBROUTINE JacobCloudAmtFM1(raFreq,raaRad,raaRadDT, &
    iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
    raaOneMinusTau,raaTau, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME, &
    raaLay2Sp,raResults,raThermal,raaLay2Gnd, &
    rSatAngle,raLayAngles, &
    raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
    iOffSet)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! daaDT,daaDQ are the d/dq,d/dT matrices
! raaLay2Sp   is the layer-to-space abs coeff matrix
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
    INTEGER :: iNumGasesTemp,iAtm,iM
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),raFreq(kMaxPts)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iG,iLay,iNumLayer,iOffSet

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
    REAL :: raTemp(kMaxPtsJac),rCos
    REAL :: raResultsTh(kMaxPtsJac)

! figure out which of 1..100 this current radiating layer corresponds to
! bleh
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
    rCos = 1.0/cos(rSatAngle*kPi/180.0)
    rCos = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
    raResults = raaGeneral(:,iLay)

! note that      iLay + iOffset === iaRadlayer(iLay)
! set the constant factor we have to multiply results with
    IF (iIWPorDME == 1) THEN
      raTemp = raaExtJacobIWP(:,iLay+iOffSet)
    ELSEIF (iIWPorDME == -1) THEN
      raTemp = raaExtJacobDME(:,iLay+iOffSet)
    END IF
    CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
    IF (iLay < iNumLayer) THEN
      ! this is not the topmost layer
      iJ1 = iLay
      raResults = raResults+raTemp*raaRad(:,iJ1)*raaLay2Sp(:,iJ1)
    ELSE IF (iLay == iNumLayer) THEN
      ! do the topmost layer correctly
      iJ1 = iLay
      raResults = raResults+raTemp*raaTau(:,iJ1)*raaRad(:,iJ1)
    END IF

! now multiply results by the 1/cos(viewing angle) factor
    IF (abs(rCos-1.0000000) >= 1.0E-5) THEN
      raResults = raResults*rCos
    END IF

    RETURN
    end SUBROUTINE JacobCloudAmtFM1

!************************************************************************
END MODULE jac_pclsam_down
