c Copyright 1997 
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

c************************************************************************
c************* THESE ARE THE JACOBIANS FOR THE DOWN LOOK INSTR***********
c************************************************************************
c this is for clear sky
c this subroutine does the Jacobians for downward looking instrument
      SUBROUTINE DownwardJacobian(raFreq,iTag,iActualTag,
     $            iProfileLayers,raPressLevels,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,rDelta,
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
      REAL raaAbs(kMaxPts,kMixFilRows)
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

c local variables
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
      INTEGER iG,iLay,iIOUN,iLowest
      INTEGER DoGasJacob,iGasJacList
      INTEGER WhichGasPosn,iGasPosn

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

      iLowest = iaaRadLayer(iAtm,1)
      iLowest = MOD(iLowest,kProfLayer)

      iIOUN = kStdJacob
c initialise the layer-to-space matrix
      CALL AtmosLayer2Space(raaLay2Sp,
     $  rSatAngle,raaAbs,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

        IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
          write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
        END IF

        CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq,
     $              iAtm,iNumLayer,iaaRadLayer,
     $              rSatAngle,raLayAngles,raaAbs,
     $              raaRad,raaRadDT,raaOneMinusTau,raaTau,
     $              raaLay2Gnd,iProfileLayers,raPressLevels)

      IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
        write(kStdWarn,*)'initializing Jacobian loops ...'
      END IF

      CALL Loop_LookDown(iaaRadLayer,
     $        iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl,
     $        raUseEmissivity,raSurface,raSun,raSunAngles,raThermal,
     $        raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
        IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
          write(kStdWarn,*)'initializing Jacobian thermal loops ...'
        END IF
        CALL Loop_thermaldown(raaRad,rSatAngle,raLayAngles,
     $          iProfileLayers,raPressLevels,
     $          iaaRadLayer,iAtm,iNumLayer,raaAbs,
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
          iLay=1
c remember that raSun is that at the gnd -- we have to propagate this back 
c up to the top of the atmosphere
c note that here we are calculating the SOLAR contribution 
          DO iG=1,kMaxPts
            radSolardT(iG) = raUseEmissivity(iG)*
     $                     raaLay2Sp(iG,iLay)*raSun(iG)
            END DO
          END IF

        CALL Find_BT_rad(raInten,radBTdr,raFreq,
     $                   radBackgndThermdT,radSolardT)

        END IF

      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. -2) 
     $    .OR. (iWhichJac .EQ. 20)) THEN
        DO iG=1,iNumGases
c for each of the iNumGases whose ID's <= kMaxDQ
c have to do all the iNumLayer radiances
          iGasJacList = DoGasJacob(iaGases(iG),iaJacob,iJacob)
          IF (iGasJacList .GT. 0) THEN
            iGasPosn = WhichGasPosn(iaGases(iG),iaGases,iNumGases)
            CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),
     $          raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            !! see if this gas does exist for this chunk
            CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag,
     $                         iDoAdd,iErr)
            IF (iDoAdd .GT. 0) THEN             
              DO iLay=1,iNumLayer
                rWeight = raaMix(iaaRadLayer(iAtm,iLay),iG)
                IF (iLay .EQ. 1) THEN
                  rWeight = rWeight*rFracBot
                ELSEIF (iLay .EQ. iNumLayer) THEN
                  rWeight = rWeight*rFracTop
                  END IF
                IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $             (abs(kLongOrShort) .LE. 1)) THEN
                  write(kStdWarn,*)'gas d/dq : gas# iaaRadlayer# :',iG,iaaRadLayer(iAtm,iLay)
                END IF
                CALL JacobGasAmtFM1(raFreq,raaRad,raaRadDT,iGasJacList,
     $            iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $            raaOneMinusTau,raaTau,raaaAllDQ,
     $            raaLay2Sp,raResults,raThermal,raaLay2Gnd,
     $            rSatAngle,raLayAngles,
     $            raaGeneral,raaGeneralTh,raaOneMinusTauTh,
     $            rWeight)
                CALL doJacobOutput(iLowest,raFreq,
     $             raResults,radBTdr,raaAmt,raInten,iaGases(iG),iLay,iGasPosn)
                CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
                END DO
            ELSE
              DO iFr = 1,kMaxPts
                raResults(iFr) = 0.0
                END DO
              DO iLay=1,iNumLayer
                CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
                END DO
              END IF
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

c then do the temperatures d/dT
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),
     $                       rDelta,iAtm,0,iNumLayer)
      IF ((iWhichJac .EQ. -1) .OR. (iWhichJac .EQ. 30) .OR.
     $    (iWhichJac .EQ. -2) .OR. (iWhichJac .EQ. 32)) THEN
        DO iLay=1,iNumLayer
c for each of the iNumLayer radiances, cumulatively add on all 
c iNumGases contributions (this loop is done in JacobTemp)
          IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $         (abs(kLongOrShort) .LE. 1)) THEN
            write(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay)
          END IF
          IF (iNatm .GT. 1) THEN
            rWeight = 0.0
            DO iG=1,iNumGases
              rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iLay),iG)            
              END DO
            rWeight = rWeight/(iNumGases*1.0)
            ELSE
            rWeight = 1.0
            END IF
          CALL JacobTempFM1(raFreq,raaRad,raaRadDT,iLay,iNumGases,
     $            iaaRadLayer,iAtm,iNumLayer,
     $            raUseEmissivity,
     $            raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults,
     $            raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,
     $            raaGeneralTh,raaOneMinusTauTh,rWeight)
          CALL doJacobOutput(iLowest,raFreq,raResults,
     $                    radBTdr,raaAmt,raInten,0,iLay,-1)
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
        DO iLay=1,iNumLayer
          IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $       (abs(kLongOrShort) .LE. 1)) THEN
            write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
          END IF
          CALL wgtfcndown(iLay,iNumLayer,rSatAngle,raLayAngles,
     $      iaaRadLayer,iAtm,raaLay2Sp,raaAbs,raResults,rFracTop,rFracBot,
     $      iNLTEStart,raaPlanckCoeff)

c does not make sense to multiply the weighting fcns with gas amounts etc
c        CALL doJacobOutput(iLowest,raFreq,raResults,
c     $                     radBTdr,raaAmt,raInten,0,iLay,-1)
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
        iLay=1
c computing Jacobians wrt surface parameters are meaningful
        CALL JacobSurfaceTemp(raFreq,iLay,
     $      rTSurface,raUseEmissivity,raaLay2Sp,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBTdr,raaAmt,raInten,-1,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobSurfaceEmis(iLay,
     $            raSurface,raThermal,raaLay2Sp,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBTdr,raaAmt,raInten,-2,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobBackgndThermal(iLay,
     $            raaLay2Sp,raThermal,raResults) 
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radBackgndThermdT,raaAmt,raInten,-3,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

        CALL JacobSolar(iLay,raaLay2Sp,raSun,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults,
     $                 radSolardT,raaAmt,raInten,-4,iLay,-1)
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
c this subroutine computes the Layer_to_Space transmission coeffs for the 
c current atmosphere
      SUBROUTINE AtmosLayer2Space(raaLay2Space,rSatAngle,
     $    raaSumAbs,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatAngle  = current satellite view angle
c raaLay2Sp  = layer to space transmission coefficients
c raaSumAbs  = matrix containing the mixed path abs coeffs
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c raLayAngles    = layer dependent angles
      REAL raaSumAbs(kMaxPts,kMixFilRows),rSatAngle
      REAL raaLay2Space(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer)
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer

c local variables
      INTEGER iFr,iLay,iL1,iLtemp,iaRadLayer(kProfLayerJac),MP2Lay
      REAL rN,rCos

      rCos = cos(rSatAngle*kPi/180.0)

      DO iLay=1,kProfLayerJac
        DO iFr=1,kMaxPts
          raaLay2Space(iFr,iLay) = 0.0
          END DO
        END DO

c set the mixed path numbers for this particular atmosphere
c since this is being called during Jacobian calculations ==> it must have 
c worked during radiance calculations ==> no need to error check
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        END DO

c now that we know the mixed path numbers, see what layer they correspond to
c and build up the layer to space transmissions
c the top most layer is easy
      iL1 = iNumLayer
      iLay = iaRadLayer(iNumLayer)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      
      DO iFr=1,kMaxPts
        IF (raaSumAbs(iFr,iLay) .LT. 0.0) THEN
          raaLay2Space(iFr,iL1) = 1.0
        ELSE 
          raaLay2Space(iFr,iL1) = exp(-raaSumAbs(iFr,iLay)/rCos)
          END IF
        END DO

c set "higher" layer number
      iLtemp = iL1

c now iteratively do the rest of the layers
      DO iL1 = iNumLayer-1,1,-1
        iLay = iaRadLayer(iL1)
        rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
        DO iFr=1,kMaxPts
          IF (raaSumAbs(iFr,iLay) .LT. 0.0) THEN
            rN=0.0
          ELSE 
            rN = raaSumAbs(iFr,iLay)/rCos
          END IF
          raaLay2Space(iFr,iL1) = exp(-rN)*raaLay2Space(iFr,iLtemp)
          END DO
c set "higher" layer number
          iLtemp = iL1
        END DO 

      RETURN
      END

c************************************************************************
c this subroutine calculates the Planck radiances, and the derivatives
c for DOWNWARD looking instrument
      SUBROUTINE DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq,
     $    iAtm,iNumLayer,iaaRadLayer,rSatAngle,raLayAngles,raaAbs,
     $    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd,
     $    iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raLayAngles        == array of layer dependent satellite view angles
c rSatAngle          == Satellite View Angle
c raSurface          == planckian emission from surface 
c raaAbs             == total absorption coeffs
c raFreq            == frequency array
c raRad              == Planck radiations for the layers
c raRadDT            == d/dT(Planck radiations) for the layers
c raVTemp            == mix vertical temperature of the layers
c iAtm               == atmosphere number
c iNumLayer          == number of layers in atmosphere
c iaaRadLayer        == radiating atmophere mixed path info
c raaLay2Gnd       == tau (layer-to-gnd) for the thermal backgnd
c raaOneMinusTau     == B(Ti)(tau(I+1)-tau(i))
c raaTau             == B(Ti)(tau(I+1)-tau(i))
c              surface,solar and backgrn thermal at the surface
c rFracTop/rFracBot  == fractional top/bottom layers
      REAL raaAbs(kMaxPts,kMixFilRows),rSatAngle,raPressLevels(kProfLayer+1)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rFracTop,rFracBot
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac)
      REAL raVTemp(kMixFilRows),raFreq(kMaxPts)
      REAL raLayAngles(kProfLayer)
      INTEGER iAtm,iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iProfileLayers

c local variables
      REAL r1,r2,r3,r4,r5,rAngle,rCos,raVT1(kMixFilRows),InterpTemp,ttorad
      INTEGER iLay,iFr,iL,iM2,iMM2,MP2Lay

      ! need these for derivative of Planck
      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

      rCos = cos(rSatAngle*kPi/180.0)

      DO iL=1,kMixFilRows
        raVT1(iL) = raVTemp(iL)
        END DO
c if the bottom layer is fractional, interpolate!!!!!!
      iL = iaaRadLayer(iAtm,1)
      raVT1(iL) = 
     $         interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
c if the top layer is fractional, interpolate!!!!!!
      iL = iaaRadLayer(iAtm,iNumLayer)
      raVT1(iL) = 
     $        interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)

cdebug
c we also allow the user to compute the temperature jacobians in 
c one of three ways, to test strength of terms in d/dT
c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )
c
c kTempJac =  0 uses both dB/dT * tau, d(tau terms)/dT * B
c kTempJac = -2 uses      dB/dT * tau, d(tau terms)/dT = ==0
c remember gas amount jacobians do not need dB/dT, so they should be OK
      IF ((kTempJac .EQ. 0) .OR. (kTempJac .EQ. -2)) THEN
        DO iL=1,iNumLayer
          iLay = iaaRadLayer(iAtm,iL)
          DO iFr=1,kMaxPts
            r3 = r1*(raFreq(iFr)**3)
c note that the data will be stored in "layer" iL, while the temperature
c is that of raVT1(iLay) ... there should be an equivalence
            r4 = r2*raFreq(iFr)/raVT1(iLay)
            r5=exp(r4)
            raaRad(iFr,iL) = r3/(r5-1.0)
            raaRadDT(iFr,iL) = raaRad(iFr,iL)*r4*r5/(r5-1.0)/raVT1(iLay)
            END DO
         END DO

c this turns off dB(T)/dT while it leaves B(T)
c so don't turn off dK/dT  ==> used for D(tau)/DT 
c thus use this if kTempJac=-1
c kTempJac = -1 uses      dB/dT == 0, d(tau terms)/dT * B
c remember gas amount jacobians do not need dB/dT, so they should be OK
      ELSE IF (kTempJac .EQ. -1) THEN
        DO iL=1,iNumLayer
          iLay = iaaRadLayer(iAtm,iL)
          DO iFr=1,kMaxPts
            r3 = r1*(raFreq(iFr)**3)
c note that the data will be stored in "layer" iL, while the temperature
c is that of raVT1(iLay) ... there should be an equivalence
            r4 = r2*raFreq(iFr)/raVT1(iLay)
            r5=exp(r4)
            raaRad(iFr,iL) = r3/(r5-1.0)
            raaRadDT(iFr,iL) = 0.0 
            END DO 
          END DO           
        END IF

c note that before using, we still have to multiply by raaRad or raaRadDT
      DO iL=1,1
c first find the mixed path number
        iLay = iaaRadLayer(iAtm,iL)
        rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaTau(iFr,iL) = exp(-raaAbs(iFr,iLay)*rFracBot/rCos)
          raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
          END DO
        END DO
      DO iL=2,iNumLayer-1
c first find the mixed path number
        iLay = iaaRadLayer(iAtm,iL)
        rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaTau(iFr,iL) = exp(-raaAbs(iFr,iLay)/rCos)
          raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
          END DO
        END DO
      DO iL = iNumLayer,iNumLayer
c first find the mixed path number
        iLay = iaaRadLayer(iAtm,iL)
        rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaTau(iFr,iL) = exp(-raaAbs(iFr,iLay)*rFracTop/rCos)
          raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
          END DO
        END DO

c check to see if we need the thermal contribution for iDownward = 1
      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
        rAngle = kThermalAngle*kPi/180.0
        rCos = cos(rAngle)
c compute Lay2Gnd (needed for the thermal backgnd inclusion in Jacobians)
c initialize bottommost layer
        iLay = iaaRadLayer(iAtm,1)
        iM2 = iaaRadLayer(iAtm,1)
        iMM2=1
c remember r4 is the 1/cos(theta) weighting factor of the satellite 
c viewing angle while we need 1/cos(theta_thermal_diffuse)
        DO iFr=1,kMaxPts
          raaLay2Gnd(iFr,1) = exp(-raaAbs(iFr,iLay)*rFracBot/rCos)
          END DO
c now go layer by layer from the bottom up to build the transmission matrix
        DO iL=2,iNumLayer-1
          iLay = iaaRadLayer(iAtm,iL)
          DO iFr=1,kMaxPts
            raaLay2Gnd(iFr,iL) = raaLay2Gnd(iFr,iMM2)*
     $                         exp(-raaAbs(iFr,iLay)/rCos)
            END DO
          iMM2 = iL
          END DO
        DO iL = iNumLayer,iNumLayer
          iLay = iaaRadLayer(iAtm,iL)
          DO iFr=1,kMaxPts
            raaLay2Gnd(iFr,iL) = raaLay2Gnd(iFr,iMM2)*
     $           exp(-raaAbs(iFr,iLay)*rFracTop/rCos)
            END DO
          iMM2 = iL
          END DO

        END IF

      RETURN
      END

c************************************************************************
c this subroutine does the general looping for the Jacobians, so all that
c has to be called is raaResults with the appropriate raaDT or raaaDQ
      SUBROUTINE Loop_LookDown(iaaRadLayer,
     $        iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl,
     $        raUseEmissivity,raSurface,raSun,raSunAngles,raThermal,
     $        raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raLayAngles is the array of layer dependent satellite view angles
c raSunAngles is the array of layer dependent satellite view angles
c rSatAngle is the satellite viewing angle
c raSurface is the surface emission
c iNumLayer is the number of layers in the atmosphere
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raaRad has the Planck radiances
c raaOneMinusTau has 1-tau
c iG has the gas number (1 .. iNumGases)
c raSun,raThermal are the downwelling Solar,thermal contributions
c raaGeneral has the results
c raSunRefelct has the solar reflectivity
c iaaRadLayer has the mixed paths  <--> layer info
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raSun(kMaxPts),raThermal(kMaxPts),rSatAngle
      INTEGER iAtm
      REAL raUseEmissivity(kMaxPts),raSurface(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raSunRefl(kMaxPts)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneral(kMaxPtsJac,kProfLayerJac)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      INTEGER iNumLayer

c local variables
      INTEGER iFr,iJ,iJ1,iLyr,iLay,MP2Lay
      REAL raTemp(kMaxPtsJac),raTemp1(kMaxPtsJac),rCos,rWsun

      rCos = cos(rSatAngle*kPi/180.0)
      rWsun=cos(kSolarAngle*kPi/180.0)

c do the bottommost layer first
      iLyr=1
c first do the surface term
      iJ1=1
      CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
      DO iFr=1,kMaxPts
        raaGeneral(iFr,iLyr) = 
     $         raUseEmissivity(iFr)*raSurface(iFr)*raTemp(iFr)
        END DO
c recall raTemp is raL2S from gnd to top
      iLay = iaaRadLayer(iAtm,iLyr)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      rWsun=cos(raSunAngles(MP2Lay(iLay))*kPi/180.0)

      IF (kSolar .GE. 0) THEN
        DO iFr=1,kMaxPts
          raaGeneral(iFr,iLyr) = raaGeneral(iFr,iLyr)+
     $        raSunRefl(iFr)*raSun(iFr)*raTemp(iFr)*(1+rCos/rWsun)
          END DO
        END IF   
c include the EASY part of thermal contribution
      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
        DO iFr=1,kMaxPts
          raaGeneral(iFr,iLyr) = raaGeneral(iFr,iLyr)+
     $         (1.0-raUseEmissivity(iFr))/kPi*raThermal(iFr)*raTemp(iFr)
          END DO
        END IF   

cccdebug === this aids in turn off surface term for both DB/DT amd D(tau)/DT 
ccc      IF (kTempJac .NE. 0) THEN
ccc        DO iFr=1,kMaxPts 
ccc          raaSurfLoop(iFr,iLyr) = raaGeneral(iFr,iLyr) 
ccc          END DO 
ccc        END IF

c set raTemp1 to all zeros (since this is the bottom layer, there is no 
c cumulative contribution
      DO iFr=1,kMaxPts
        raTemp1(iFr) = 0.0
        END DO

c loop over the remaining layers
      DO iLyr=2,iNumLayer
        iLay = iaaRadLayer(iAtm,iLyr)
        rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
c first do the surface term
        iJ1=1
        CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
        DO iFr=1,kMaxPts
          raaGeneral(iFr,iLyr) = 
     $         raUseEmissivity(iFr)*raSurface(iFr)*raTemp(iFr)
          END DO

c recall raTemp is raL2S from gnd to top
        IF (kSolar .GE. 0) THEN
          rWsun=cos(raSunAngles(MP2Lay(iLay))*kPi/180.0)
          DO iFr=1,kMaxPts
            raaGeneral(iFr,iLyr) = raaGeneral(iFr,iLyr)+
     $        raSunRefl(iFr)*raSun(iFr)*raTemp(iFr)*(1+rCos/rWsun)
            END DO
          END IF   
c include the EASY part of thermal contribution
        IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
          DO iFr=1,kMaxPts
            raaGeneral(iFr,iLyr) = raaGeneral(iFr,iLyr)+
     $         (1.0-raUseEmissivity(iFr))/kPi*raThermal(iFr)*raTemp(iFr)
            END DO
          END IF   

cccdebug === this aids in turn off surface term for both DB/DT amd D(tau)/DT 
ccc        IF (kTempJac .NE. 0) THEN
ccc          DO iFr=1,kMaxPts 
ccc            raaSurfLoop(iFr,iLyr) = raaGeneral(iFr,iLyr) 
ccc            END DO 
ccc          END IF

c now loop over the layers that contribute (i.e. < iLyr) ....
        iJ = iLyr-1
        iJ1 = iJ+1
        CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
        DO iFr=1,kMaxPts
          raTemp1(iFr) = raTemp1(iFr)+
     $                 raaOneMinusTau(iFr,iJ)*raaRad(iFr,iJ)*raTemp(iFr) 
          raaGeneral(iFr,iLyr) = raaGeneral(iFr,iLyr)+raTemp1(iFr)
          END DO 

        END DO

      RETURN
      END 
c************************************************************************
c the easy part of backgnd thermal Jacobians is done in Loop_lookdown([])
c this subroutine does the hard part of backgnd thermal Jacobians
      SUBROUTINE Loop_thermaldown(raaRad,rSatAngle,raLayAngles,
     $                iProfileLayers,raPressLevels,
     $                iaaRadLayer,iAtm,iNumLayer,raaSumAbsCoeffs,
     $                raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c raaRad has the Planck radiances
c raaGeneralTh has the results
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
c raLayAngles is the array of layer dependent satellite view angles
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iProfileLayers
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raaRad(kMaxPtsJac,kProfLayerJac),rSatAngle
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer)
      REAL raaSumAbsCoeffs(kMaxPts,kMixFilRows)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      REAL raFreq(kMaxPts),raPressLevels(kProfLayer+1)
      INTEGER iNumLayer

c local variables
      INTEGER iFr,iJ,iLay,iL,iLyr,iaRadLayer(kProfLayer)
      REAL raTemp(kMaxPtsJac),rTh

      INTEGER FindBoundary,iB

      DO iL = 1,iNumLayer
        iaRadLayer(iL) = iaaRadLayer(iAtm,iL)
        END DO

c since we are using acos(3/5) approx here, instead of the more accurate 
c diffusive approximation, might as well also approximate these contributions,
c so as to speed up the code
c*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****

      iB = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
      iB = INT((iB*1.0)/2.0)

      rTh=cos(kThermalAngle*kPi/180.0)

      DO iJ=1,iNumLayer
        iL = iaaRadLayer(iAtm,iJ)
        DO iFr=1,kMaxPts
          raaOneMinusTauTh(iFr,iJ) = 1.0-exp(-raaSumAbsCoeffs(iFr,iL)/rTh)
          END DO
        END DO

c set the contribution of the upper layers to 0
c*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****
      DO iLyr = iB+1,iNumLayer
        iLay = iLyr
        DO iFr=1,kMaxPts
          raaGeneralTh(iFr,iLyr) = 0.0
          END DO
        END DO

c this is "hard" part of the thermal, where we loop over lower layers 
c that contribute
c*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****

c first set the matrix elements for layer iB
      iLyr = iB
      DO iFr=1,kMaxPts
        raaGeneralTh(iFr,iLyr) = 0.0
        END DO
      DO iJ = iNumLayer,iLyr+1,-1
        CALL JacobTermGnd(iJ-1,iLyr,raaLay2Gnd,raTemp)
        DO iFr=1,kMaxPts
          raaGeneralTh(iFr,iLyr) = raaGeneralTh(iFr,iLyr)+
     $        raaOneMinusTauTh(iFr,iJ)*raaRad(iFr,iJ)*raTemp(iFr)
          END DO
        END DO

c now do the rest of the layers
      DO iLyr = iB-1,1,-1
        DO iFr=1,kMaxPts
          raaGeneralTh(iFr,iLyr) = 0.0
          END DO
        iJ = iLyr+1
        CALL JacobTermGnd(iJ-1,iLyr,raaLay2Gnd,raTemp)
        DO iFr=1,kMaxPts
          raaGeneralTh(iFr,iLyr) = raaGeneralTh(iFr,iLyr+1)+
     $        raaOneMinusTauTh(iFr,iJ)*raaRad(iFr,iJ)*raTemp(iFr)
          END DO
        END DO

      RETURN
      END 
c************************************************************************
c this subroutine does the Jacobians wrt gas amount
      SUBROUTINE JacobGasAmtFM1(raFreq,raaRad,raaRadDT,
     $    iG,iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,
     $    raUseEmissivity,
     $    raaOneMinusTau,raaTau,raaaAllDQ,raaLay2Sp,raResults,
     $    raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,
     $    raaGeneralTh,raaOneMinusTauTh,rWeight)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rWeight is the wieght of this gas from *WEIGHT
c rSatAngle is the satellite viewing angle
c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c daaDT,daaDQ are the d/dq,d/dT matrices
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raaaAllDQ has the ALL the d/dq coeffs for current freq block
c raaRad has the Planck radiances
c raaRadDT has the d/DT (Planck radiances)
c raaOneMinusTau has 1-tau(satellite angle)
c raaOneMinusTauTh has 1-tau(thermal angle)
c iG has the gas number (1 .. iNumGases)
c iLay has info on how to find the radiating layer number (1..kProfLayerJac)
c raFreq has the frequencies
c raResults has the results
c raThermal are the downwelling Solar,thermal contributions
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
c raaGeneral,raaGeneralTh have the general results (looping over layers)
      REAL raaGeneral(kMaxPtsJac,kProfLayerJac),rWeight
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle
      REAL raLayAngles(kProfLayer)
      REAL raThermal(kMaxPts)
      INTEGER iNumGases,iAtm
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raFreq(kMaxPts)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)
      INTEGER iG,iLay,iNumLayer

c local variables
      INTEGER iFr,iJ1,iM1,MP2Lay
      REAL raTemp(kMaxPtsJac),rSec
      REAL raResultsTh(kMaxPtsJac)

c      IF (iGasORCld .LE. 0) THEN 
c        write(kStdErr,*) 'only do amount jacobians for gases in this routine'
c        CALL DoStop
c        END IF

c figure out which of 1..100 this current radiating layer corresponds to
c bleh
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the sat angle weight factor
      rSec = 1.0/cos(rSatAngle*kPi/180.0)
      rSec = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

c read the appropriate layer from general results
      DO iFr=1,kMaxPts
        raResults(iFr) = raaGeneral(iFr,iLay)
        END DO

c set the constant factor we have to multiply results with
      !!this is a gas amt jacobian
      DO iFr=1,kMaxPts
        raTemp(iFr) = raaaAllDQ(iG,iFr,iM1)
        END DO
c      print *,iG,iAtm,iLay,iM1,rSec,raResults(1),raTemp(1)

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
      IF (abs(rSec-1.0000000) .GE. 1.0E-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*rSec
          END DO
        END IF

c see if we have to include thermal backgnd
      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
        CALL JacobTHERMALAmtFM1(raFreq,raaRad,
     $       iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,
     $       raUseEmissivity,raTemp,raaLay2Sp,
     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
c now add on the effects to raResults
        DO iFr=1,kMaxPts
          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
          END DO
        END IF

c now multiply results by the weight factor
      IF (abs(rWeight-1.0000000) .GE. 1.0E-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*rWeight
          END DO
        END IF

      RETURN
      END 
c************************************************************************
c this subroutine does the Jacobians wrt temperature
      SUBROUTINE JacobTempFM1(raFreq,raaRad,raaRadDT,iLay,iNumGases,
     $  iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,
     $  raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults,
     $  raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,
     $  raaGeneralTh,raaOneMinusTauTh,rWeight)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rWeight is the cumulative average of all weights
c rSatAngle is the satellite viewing angle
c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c daaDT,daaDQ are the d/dq,d/dT matrices
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raaAllDT has the CUMULATIVE the d/dT coeffs for current freq block
c raaRad has the Planck radiances
c raaRadDT has the d/DT (Planck radiances)
c raaOneMinusTau has 1-tau(satellite angle)
c raaOneMinusTauTh has 1-tau(thermal angle)
c iLay has info on how to find the radiating layer number (1..kProfLayerJac)
c raFreq has the frequencies
c raResults has the results
c raThermal are the downwelling Solar,thermal contributions
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
c raaGeneral,raaGeneralTh have the general resulkts (looping over layers)
      REAL raaGeneral(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle,rWeight
      REAL raThermal(kMaxPts),raLayAngles(kProfLayer)
      INTEGER iNumGases,iAtm,iLay,iNumLayer
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
      REAL raaTau(kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raFreq(kMaxPts)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iJ1,iJp1,iM1,MP2Lay
      REAL raTemp(kMaxPtsJac),rSec,rEmittance
      REAL raResultsTh(kMaxPtsJac)

c figure out which of 1..100 this current radiating layer corresponds to
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the sat angle weight factor
      rSec = 1.0/cos(rSatAngle*kPi/180.0)
      rSec = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

c read the appropriate layer from general results
cdebug
c we also allow the user to compute the temperature jacobians in 
c one of three ways, to test strength of terms in d/dT
c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )
c
c this includes all the surface terms 
c      IF (kTempJac .EQ. 0) THEN
c        DO iFr=1,kMaxPts
c          raResults(iFr) = raaGeneral(iFr,iLay)
c          END DO
c       ELSE
cc this turns off the surface terms if we only want to check DB/DT or D(tau)/DT
cc this way we only have the loop over the layers contributing
c        DO iFr=1,kMaxPts 
c          raResults(iFr) = raaGeneral(iFr,iLay)-raaSurfLoop(iFr,iLay) 
c          END DO 
c        END IF

c this includes all the surface terms 
       DO iFr=1,kMaxPts
         raResults(iFr) = raaGeneral(iFr,iLay)
         END DO

c now set the constant factor we have to multiply results with
cdebug
c we also allow the user to compute the temperature jacobians in 
c one of three ways, to test strength of terms in d/dT
c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
c so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
c so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all
c
c this turns off dK/dT  ==> use for kTempJac=-2 (DB/DT )
      IF (kTempJac .EQ. -2) THEN
        DO iFr=1,kMaxPts 
          raTemp(iFr) = 0.0 
          raResults(iFr) = 0.0
          END DO 
c else if we want dK/dT use this 
c ie use for kTempJac = 0 (DB/DT, d(tau)/DT), -1 (D(tau)/DT )
      ELSE IF (kTempJac .NE. -2) THEN
        DO iFr=1,kMaxPts
          raTemp(iFr) = raaAllDT(iFr,iM1)
          END DO
        CALL MinusOne(raTemp,raResults)
        END IF

c now do the derivatives wrt radiating layer
      IF (iLay .LT. iNumLayer) THEN
c this is not the topmost layer
        iJ1 = iLay
        iJp1 = iLay+1
        DO iFr=1,kMaxPts
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Sp(iFr,iJ1) +
     $               raaOneMinusTau(iFr,iJ1)*raaRadDT(iFr,iJ1)/
     $               rSec/rWeight*raaLay2Sp(iFr,iJp1)
          raResults(iFr) = raResults(iFr)+rEmittance
          END DO 
      ELSE IF (iLay .EQ. iNumLayer) THEN
c do the topmost layer correctly
        iJ1 = iLay
        DO iFr=1,kMaxPts
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Sp(iFr,iJ1)+
     $         raaOneMinusTau(iFr,iJ1)*raaRadDT(iFr,iJ1)/rSec/rWeight
          raResults(iFr) = raResults(iFr)+rEmittance
          END DO 
        END IF

c now multiply results by the 1/cos(viewing angle) factor
      IF (abs(rSec-1.00000) .GE. 1.0e-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*rSec
          END DO
        END IF

c now add on the effects to raResults
      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
c this subroutine does the Jacobians wrt temperature
        CALL JacobTHERMALTempFM1(raFreq,raaRad,raaRadDT,
     $       iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,
     $       raUseEmissivity,raTemp,raaLay2Sp,
     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
        DO iFr=1,kMaxPts
          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
          END DO
        END IF

c now multiply results by the weight factor
      IF (abs(rWeight-1.0000000) .GE. 1.0E-5) THEN
        DO iFr=1,kMaxPts
          raResults(iFr) = raResults(iFr)*rWeight
          END DO
        END IF

      RETURN
      END 

c************************************************************************
c this subroutine does the hard part of backgnd thermal Jacobians wrt amt
        SUBROUTINE JacobTHERMALAmtFM1(raFreq,raaRad,
     $       iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,
     $       raUseEmissivity,raTemp,raaLay2Sp,
     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raaaAllDQ has the ALL the d/dq coeffs for current freq block
c raaRad has the Planck radiances
c iLay has info on how to find the radiating layer number (1..kProfLayerJac)
c raFreq has the frequencies
c raResults has the results
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      INTEGER iNumGases,iAtm,iLay,iNumLayer
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raTemp(kMaxPtsJac),raResultsTh(kMaxPtsJac)
      REAL raFreq(kMaxPts),raaRad(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iJ1,iM1,MP2Lay
      REAL rEmittance,rTempTh,rXYZ

c figure out which of 1..100 this current radiating layer corresponds to
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the thermal angle weight factor
      rTempTh=cos(kThermalAngle*kPi/180.0)

c read the appropriate layer from general results
      DO iFr=1,kMaxPts
        raResultsTh(iFr) = raaGeneralTh(iFr,iLay)
        END DO

c we have already set the constant factor we have to multiply results with
c ie raTemp is already the relevant row of raaaDQ for gas Jacobian
c        or is already the relevant row of raaDt for temperature Jacobian
c        depending on the value of kTempJac
cdebug
c we also allow the user to compute the temperature jacobians in 
c one of three ways, to test strength of terms in d/dT
c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

      CALL MinusOne(raTemp,raResultsTh)

c this is the part where we include the effect of the radiating layer
      IF ((iLay. GT. 1) .AND. (iLay .LE. iNumLayer)) THEN
        iJ1 = iLay
        DO iFr=1,kMaxPts
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)
          raResultsTh(iFr) = raResultsTh(iFr)+rEmittance
          END DO 
      ELSE IF (iLay .EQ. 1) THEN
c do the bottommost layer correctly
        iJ1 = iLay
        DO iFr=1,kMaxPts
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)
          raResultsTh(iFr) = raResultsTh(iFr)+rEmittance
          END DO 
        END IF

c now multiply results by the tau(layer_to_space)
c include a diffusivity factor of 0.5 and a factor of 2pi (azimuth integ)
c thus (1-ems)/pi * (2pi) *(0.5) === (1-ems)
      DO iFr=1,kMaxPts
        rXYZ = (1.0-raUseEmissivity(iFr))/rTempTh
        raResultsTh(iFr) = raResultsTh(iFr)*rXYZ*raaLay2Sp(iFr,1)
        END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the hard part of the THERMAL Jacobians wrt temperature
        SUBROUTINE JacobTHERMALTempFM1(raFreq,raaRad,raaRadDT,
     $       iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,
     $       raUseEmissivity,raTemp,raaLay2Sp,
     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNumLayer is the number of layers in the atmosphere
c iaaRadlayer has the radiating layer information for atmospher # iAtm
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raaaAllDQ has the ALL the d/dq coeffs for current freq block
c raaRad has the Planck radiances
c iG has the gas number (1 .. iNumGases)
c iLay has info on how to find the radiating layer number (1..kProfLayerJac)
c raFreq has the frequencies
c raResults has the results
c raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
      REAL raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
      REAL raaGeneralTh(kMaxPtsJac,kProfLayerJac)
      INTEGER iNumGases,iAtm,iLay,iNumLayer
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raUseEmissivity(kMaxPts)
      REAL raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
      REAL raTemp(kMaxPtsJac),raResultsTh(kMaxPtsJac)
      REAL raFreq(kMaxPts),raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raaRad(kMaxPtsJac,kProfLayerJac)
      REAL raaRadDT(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iJ1,iJm1,iM1,MP2Lay
      REAL rEmittance,rTempTh,rW

c figure out which of 1..100 this current radiating layer corresponds to
      iM1 = iaaRadLayer(iAtm,iLay)
      iM1 = MP2Lay(iM1)

c fix the thermal angle weight factor
      rTempTh=cos(kThermalAngle*kPi/180.0)

c read the appropriate layer from general results
      DO iFr=1,kMaxPts
        raResultsTh(iFr) = raaGeneralTh(iFr,iLay)
        END DO

c we have already set the constant factor we have to multiply results with
      CALL MinusOne(raTemp,raResultsTh)

c this is the part where we include the effect of the radiating layer
      IF ((iLay. GT. 1) .AND. (iLay .LE. iNumLayer)) THEN
        iJ1 = iLay
        iJm1 = iJ1-1
        DO iFr=1,kMaxPts
          rEmittance = raaOneMinusTauTh(iFr,iJ1)
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)+
     $     rEmittance/rTempTh*raaRadDT(iFr,iJ1)*raaLay2Gnd(iFr,iJm1)
          raResultsTh(iFr) = raResultsTh(iFr)+rEmittance
          END DO 
      ELSE IF (iLay .EQ. 1) THEN
c do the bottommost layer correctly
        iJ1 = iLay
        DO iFr=1,kMaxPts
          rEmittance = raaOneMinusTauTh(iFr,iJ1)
          rEmittance = raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)+
     $               rEmittance/rTempTh*raaRadDT(iFr,iJ1)
          raResultsTh(iFr) = raResultsTh(iFr)+rEmittance
          END DO 
        END IF

c now multiply results by the tau(layer_to_space)
c include a diffusivity factor of 0.5 
c thus (1-ems)/pi * (2pi) *(0.5) === (1-ems)
      DO iFr=1,kMaxPts
        rW=(1.0-raUseEmissivity(iFr))/rTempTh
        raResultsTh(iFr) = raResultsTh(iFr)*rW*raaLay2Sp(iFr,1)
        END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the weighting functions for downward looking instr 
      SUBROUTINE wgtfcndown(iLay,iNumLayer,rSatAngle,raLayAngles, 
     $                      iaaRadLayer,iAtm,raaLay2Sp,raaAbs,raResults, 
     $                      rFracTop,rFracBot,
     $                      iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
c raaSumAbCoeff is the cumulative absorption coeffs 
c iNumLayer is the number of layers in the atmosphere 
c iaaRadlayer has the radiating layer information for atmospher # iAtm 
c raaLay2Sp   is the layer-to-space abs coeff matrix 
c iLay has info on how to find the radiating layer number (1..kProfLayerJac) 
c raResults has the results 
c raLayAngles has the layer dependent satellite view angles 
      REAL raLayAngles(kProfLayer) 
      INTEGER iAtm,iLay,iNumLayer 
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer) 
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac) 
      REAL raResults(kMaxPtsJac),rSatAngle,rFracTop,rFracBot 
      REAL raaAbs(kMaxPts,kMixFilRows) 
c this is for NLTE weight fcns
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
 
      REAL rCos 
      INTEGER iFr,iM,iM1,MP2Lay
      INTEGER iModKprofLayer,mod
 
      rCos = cos(raLayAngles(1)*kPi/180.0) 
      rCos = cos(rSatAngle*kPi/180.0) 

      DO iFr=1,kMaxPts 
        raResults(iFr) = 0.0 
        END DO 
 
      IF (iLay .GT. iNumLayer) THEN 
        write(kStdErr,*) 'Cannot compute wt fcn for layer ',iLay 
        write(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers' 
        CALL DoSTOP 
        END IF 
 
      IF (iLay .LE. 0) THEN 
        write(kStdErr,*) 'Cannot compute wt fcn for layer ',iLay 
        write(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers' 
        CALL DoSTOP 
        END IF 
 
      IF (iLay .EQ. iNumLayer) THEN 
c use layer to space transmission iM+1 --> infinity == 1.0 
        iM = iaaRadLayer(iAtm,iLay)  
        rCos = cos(raLayAngles(MP2Lay(iM))*kPi/180.0) 
        iModKprofLayer = mod(iM,kProfLayer) 
        IF (iModKprofLayer .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = 1.0-exp(-raaAbs(iFr,iM)*rFracTop/rCos) 
            END DO 
        ELSEIF (iModKprofLayer .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = (1.0-exp(-raaAbs(iFr,iM)*rFracTop/rCos))*
     $                            raaPlanckCoeff(iFr,iM) 
            END DO 
          END IF
      ELSE IF (iLay .EQ. 1) THEN 
        iM1 = iaaRadLayer(iAtm,iLay+1) 
        iM = iaaRadLayer(iAtm,iLay) 
        iModKprofLayer = mod(iM,kProfLayer) 
        IF (iModKprofLayer .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = (1.0-exp(-raaAbs(iFr,iM)*rFracBot/rCos))* 
     $                         raaLay2Sp(iFr,iLay+1) 
          END DO 
        ELSEIF (iModKprofLayer .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = (1.0-exp(-raaAbs(iFr,iM)*rFracTop/rCos))*
     $                         raaLay2Sp(iFr,iLay+1)*raaPlanckCoeff(iFr,iM) 
            END DO 
          END IF
      ELSE 
        iM1 = iaaRadLayer(iAtm,iLay+1) 
        iM = iaaRadLayer(iAtm,iLay) 
        iModKprofLayer = mod(iM,kProfLayer) 
        IF (iModKprofLayer .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = (1.0-exp(-raaAbs(iFr,iM)/rCos))* 
     $                         raaLay2Sp(iFr,iLay+1) 
          END DO 
        ELSEIF (iModKprofLayer .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts 
            raResults(iFr) = (1.0-exp(-raaAbs(iFr,iM)*rFracTop/rCos))*
     $                         raaLay2Sp(iFr,iLay+1)*raaPlanckCoeff(iFr,iM) 
            END DO 
          END IF
        END IF 
 
      RETURN 
      END 

c************************************************************************ 
