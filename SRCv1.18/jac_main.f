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
c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

c************************************************************************
c**************************** GENERIC ROUTINES **************************
c************************************************************************
c this is the main driver subroutine for clear sky Jacobians
c for the current frequency block, this subroutine calculates ALL the 
c jacobians and then outputs them
      SUBROUTINE find_jacobians(raFreq,iTag,iActualTag,
     $            iFileID,caJacobFile,rTSpace,rTSurface,
     $            raUseEmissivity,rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,iProfileLayers,pProf,
     $            iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rDelta is the kComp file Step
c raLayAngles are the layer dependent satellite view angles
c raSunAngles are the layer dependent sun view angles
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which kcomp cm-1 block being output
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
c raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
c raaMix is the mixing table

c these are to do with the arbitrary pressure layers
      REAL raPresslevels(kProfLayer+1),raThickness(kProfLayer)
      REAL pProf(kProfLayer)
      INTEGER iProfileLayers
c FracTop,rFracBot are the upper layer/lower layer fractions
      REAL raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
      REAL raSurFace(kMaxPts),raSun(kMaxPts)
      REAL raThermal(kMaxPts),rDelta
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $      raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      INTEGER iJacob,iaJacob(kMaxDQ),iTag,iActualTag
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
      CHARACTER*80 caJacobFile
c this is for NLTE weight fcns
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)

c local variables
      INTEGER iDownWard,iI

c set the direction of radiation travel --- the checks of iUpper.iLower have
c already been done in radiance.f
c radiation travelling upwards to instrument ==> sat looking down iDownWard = 1
c radiation travelling down to instrument ==> sat looking up iDownWard =-1
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
        iDownWard = 1
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,iNumLayer))THEN
        iDownWard = -1
      END IF
      IF (abs(iDownWard) .NE. 1) THEN
        write(kStdErr,*) 'hmm : jacobian code cannot decide up/down look!'
        write(kStdErr,*) (iaaRadLayer(iAtm,iI),iI=1,iNumLayer)
        write(kStdErr,*) iAtm,iNumLayer,iaaRadLayer(iAtm,1),
     $                   iaaRadLayer(iAtm,iNumLayer)
        CALL DoStop
      END IF

      IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
        write(kStdWarn,*) 'in Jacobian, have set set iDownWard = ',iDownWard
      END IF

      IF (iDownWard .EQ. 1) THEN
        CALL DownWardJacobian(raFreq,iTag,iActualTag,
     $          iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,raSunRefl,rDelta,
     $          iNLTEStart,raaPlanckCoeff)
      ELSE IF (iDownWard .EQ. -1) THEN
        CALL UpWardJacobian(raFreq,iTag,iActualTag,
     $          iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,rDelta)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine multiplies the array by -1.0*constant where constant 
c depends on whether we are doing d/dT or d/dq
      SUBROUTINE MinusOne(raTorQ,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raResults is the array
c raTorQ === relevant element of raaaDq or raaDt
      REAL raTorQ(kMaxPtsJac)
      REAL raResults(kMaxPtsJac)

      INTEGER iFr
 
      DO iFr = 1,kMaxPts
        raResults(iFr) = -raResults(iFr) * raTorQ(iFr)
        END DO
 
      RETURN
      END

c************************************************************************
c cumulatively, using contributions from each gas, find d/dT jacobian 
c for each layer
      SUBROUTINE cumulativeDT(daaDT,raaAllDT,raaMix,iG,iNatm,
     $                        iaaRadLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaDT has the current gas d/dT coeffs for current freq block
c raaAllDT has the cumulative d/dT coeffs for current freq block
c iNatm is the number of atmospheres to do radiance calcs for
c iG is the current gas
c raaMix is the mixing table
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      INTEGER iG,iNatm,iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore)

      INTEGER iL,iFr
      REAL rW

      IF (iNatm .GT. 1) THEN  
c cannot correctly weight the d/dT, so just use unit weight here and then try
c an average weight when JacobTemp is actually called
        IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $       (abs(kLongOrShort) .LE. 1)) THEN
          write(kStdWarn,*)'Gas iG, weight rW = ',iG,1.0
        END IF

        DO iL = 1,kProfLayerJac
          DO iFr = 1,kMaxPtsJac
            raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + daaDT(iFr,iL)
            END DO
          END DO
      ELSE IF (iNatm .EQ. 1) THEN  
c have only one atmosphere and so correctly weight this gas's contribution to
c d/dT matrix ... then use weight of 1.0 when calling JacobTemp
        iL = iaaRadLayer(1,2)  !for atm#1, find which is the second mixed path
                             !as the first,last could have fractional weights
        rW = raaMix(iL,iG)   !find the gas weight in the second radiating layer
        IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
          write(kStdWarn,*)'jacobian d/dT Gas iG, weight rW = ',iG,rW
        END IF
        DO iL = 1,kProfLayerJac
          DO iFr = 1,kMaxPtsJac
            raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + rW*daaDT(iFr,iL)
            END DO
          END DO
        END IF
   
      RETURN
      END

c************************************************************************
c this subroutine does d/dr(tau_layer2space) for gas iG
c where r == gas amount q or temperature T at layer iM
c and  iL is the relevant layer we want tau_layer2space differentiated
c HENCE IF iL > iM, derivative == 0
c i.e. this does d(tau(l--> inf)/dr_m
      SUBROUTINE JacobTerm(iL,iM,raaLay2Sp,raTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Sp is the transmission frm layer to space
c iM has the layer that we differentiate wrt to
c iL has the radiating layer number (1..kProfLayerJac)
c raTemp has the results, apart from the multiplicative constants
c   which are corrected in MinusOne
      INTEGER iL,iM
      REAL raTemp(kMaxPtsJac),raaLay2Sp(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr

      IF (iL .GT. iM) THEN 
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 0.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaLay2Sp(iFr,iL)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine does d/dr(tau_layer2gnd) for gas iG
c where r == gas amount q or temperature T at layer iM
c and  iL is the relevant layer we want tau_layer2space differentiated
c HENCE IF iL < iM, derivative == 0
c i.e. this does d(tau(l--> 0)/dr_m
      SUBROUTINE JacobTermGnd(iL,iM,raaLay2Gnd,raTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Gnd is the transmission frm layer to ground at diffusion angle
c iM has the layer that we differentiate wrt to
c iL has the radiating layer number (1..kProfLayerJac)
c raTemp has the results, apart from the multiplicative constants
c   which are corrected in MinusOne
      INTEGER iL,iM
      REAL raTemp(kMaxPtsJac),raaLay2Gnd(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr

      IF (iL .LT. iM) THEN 
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 0.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaLay2Gnd(iFr,iL)
          END DO
        END IF

      RETURN
      END


c************************************************************************
c************** THESE HAVE TO DO WITH THE OUTPUT STYLE ******************
c************************************************************************
c this subroutine computes d(Brightness Temp)/d(Rad)
      SUBROUTINE Find_BT_rad(raInten,radBTdr,raFreq,
     $                       radBackgndThermdT,radSolardT)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c raInten is the radiance intensity at the instrument
c raFreq are the frequencies
c radBTdr is the derivative result
      REAL raFreq(kMaxPts),raInten(kMaxPts),radBTdr(kMaxPtsJac)
      REAL radBackgndThermdT(kMaxPtsJac),radSolardT(kMaxPtsJac)
      
      INTEGER iFr
      REAL r1,r2,r3,r4

      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

      DO iFr = 1,kMaxPts
        r3 = r1*r2 * (raFreq(iFr)**4)/(raInten(iFr)**2)
        r4 = 1.0+r1 * (raFreq(iFr)**3)/raInten(iFr)
        radBTdr(iFr) = r3/r4/(alog(r4)**2)
        END DO

      IF (kThermal .LT. 0) THEN
        DO iFr = 1,kMaxPts
          radBackgndThermdT(iFr) = 0.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          r3 = r1*r2 * (raFreq(iFr)**4)/(radBackgndThermdT(iFr)**2)
          r4 = 1.0+r1 * (raFreq(iFr)**3)/radBackGndThermdT(iFr)
          radBackgndThermdT(iFr) = r3/r4/(alog(r4)**2)
          END DO
        END IF

      IF (kSolar .LT. 0) THEN
        DO iFr = 1,kMaxPts
          radSolardT(iFr) = 0.0
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          r3 = r1*r2 * (raFreq(iFr)**4)/(radSolardT(iFr)**2)
          r4 = 1.0+r1 * (raFreq(iFr)**3)/radSolardT(iFr)
          radSolardT(iFr) = r3/r4/(alog(r4)**2)
          END DO
        END IF

      RETURN
      END
c************************************************************************
c this subroutine prepares the output Jacobians according to kJacobOutput
      SUBROUTINE doJacobOutput(iLowest,raFreq,raResults,
     $                   radBTdr,raaAmt,raInten,iGasID,iM,iGasPosn)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLowest is the lowest layer in the atmosphere (modulo kProfLayer)
c raFreq are the frequency wavenumbers
c raResults are the raw d(rad)/d(gas amt)
c raaAmt are the gas profiles
c raInten is the radiant intensity at instrument
c iGasID is well duh... actually if iGasID = -1, then we are doing d/dT
c iM is the layer number (1..100)
c iGasPosn is the position of gasID in the gaslist
c radBTdr is the d(brightness temp)/d(Radiance) array
      INTEGER iGasID,iM,iLowest,iGasPosn
      REAL raFreq(kMaxPts),raResults(kMaxPtsJac)
      REAL raInten(kMaxPts)
      REAL raaAmt(kProfLayerJac,kGasStore),radBTdr(kMaxPtsJac)

      INTEGER iFr,iM1

      IF ((iGasID .EQ. 101) .OR. (iGasID .EQ. 102)) THEN
        iGasPosn  = 1   !!!! corresponds to water
        END IF

      iM1=(iLowest-1) + iM

c      IF (kJacobOutPut .EQ. -1) THEN
c basically do nothing! user wants d(rad)/dq
c        END IF

      IF (kJacobOutput .NE. -1) THEN !oh well, do this
        IF ((iGasID .GT. 0) .AND. (iGasID .LE. 200)) THEN
c we are doing d/dq  for a normal gas
          IF (kJacobOutPut .EQ. 0) THEN
c user wants d(rad)/dq * q for a normal gas
            DO iFr = 1,kMaxPts
              raResults(iFr) = raResults(iFr) * raaAmt(iM1,iGasPosn)
              END DO
          ELSE IF (kJacobOutPut .EQ. 1) THEN
c user wants d(BT)/dq * q for a normal gas; this is the default option
            DO iFr = 1,kMaxPts
              raResults(iFr) = 
     $              raResults(iFr) * raaAmt(iM1,iGasPosn) * radBTdr(iFr)
              END DO
          ELSE IF (kJacobOutPut .EQ. 2) THEN
c user wants d(BT)/dq for a normal gas
            DO iFr = 1,kMaxPts
              raResults(iFr) = raResults(iFr) * radBTdr(iFr)
              END DO
            END IF

        ELSEIF (iGasID .GT. 200) THEN
c we are doing d/dq  for IWP or DME
          IF (kJacobOutPut .EQ. 0) THEN
c user wants d(rad)/dq * q for IWP or DME
            DO iFr = 1,kMaxPts
              raResults(iFr) = raResults(iFr)
              END DO
          ELSE IF (kJacobOutPut .EQ. 1) THEN
c user wants d(BT)/dq * q for IWP or DME for a normal gas
            DO iFr = 1,kMaxPts
              raResults(iFr) = raResults(iFr) * radBTdr(iFr)
              END DO
            END IF

        ELSE IF (iGasID .LE. 0) THEN
c we are doing d/dT or cloud amt, size jacobians
          IF (kJacobOutPut .EQ. 0) THEN
            iFr = 1 
c user wants d(rad)/dT so do nothing 
          ELSE IF (kJacobOutPut .EQ. 1) THEN
c user wants d(BT)/dT
            DO iFr = 1,kMaxPts
              raResults(iFr) = raResults(iFr) * radBTdr(iFr)
              END DO
            END IF
          END IF

        END IF   !IF (kJacobOutput .NE. -1) THEN !oh well, do this

      RETURN
      END

c************************************************************************
c************ THESE HAVE TO DO WITH THE SURFACE PARAMETERS **************
c************************************************************************
c this subroutine does Jacobian wrt Surface Temperature
      SUBROUTINE JacobSurfaceTemp(raFreq,iM,
     $      rTSurface,raUseEmissivity,raaLay2Sp,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Sp   is the layer-to-space abs coeff matrix
c raFreq has the frequencies
c raResults has the results
c iM are the layer <-> mixed path associations
      INTEGER iM
      REAL rTSurface,raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raFreq(kMaxPts)

c local variables
      REAL r1,r2,r3,r4,r5,rad,dradDT
      INTEGER iFr

      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

      DO iFr = 1,kMaxPts
        r3 = r1 * (raFreq(iFr)**3)
        r4 = r2 * raFreq(iFr)/rTSurface
        r5 = exp(r4)
        rad = r3/(r5-1.0)
        dRadDT = rad * r4 * r5/(r5-1.0)/rTSurface
        raResults(iFr) = dRadDT*raUseEmissivity(iFr) * raaLay2Sp(iFr,iM)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian wrt Surface Emissivity
      SUBROUTINE JacobSurfaceEmis(iM,raSurface,raThermal,raaLay2Sp,
     $              raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raSurface is the surface emission
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raResults has the results
c raThermal has the downwelling thermal contribs
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac),raSurface(kMaxPts)
      REAL raResults(kMaxPtsJac),raThermal(kMaxPts)

c local variables
      INTEGER iFr

      DO iFr = 1,kMaxPts
        raResults(iFr) = raaLay2Sp(iFr,iM) * 
     $                 (raSurface(iFr)-raThermal(iFr)/kPi)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian of Backgnd Thermal wrt Surface Emissivity
      SUBROUTINE JacobBackgndThermal(iM,raaLay2Sp,raThermal,raResults) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c raaLay2Sp is the layer-to-space abscoeff matrix 
c raResults has the results 
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raThermal(kMaxPts) 

c local variables
      INTEGER iFr

      DO iFr = 1,kMaxPts
        raResults(iFr) = -raaLay2Sp(iFr,iM)/kPi*raThermal(iFr) 
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian of Solar wrt Sun Surface Emissivit
      SUBROUTINE JacobSolar(iM,raaLay2Sp,raSun,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c raaLay2Sp is the layer-to-space abscoeff matrix 
c raResults has the results 
c raSun,raThermal has the downwelling solar,thermal contribs
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac),
     $     raSun(kMaxPts)

c local variables
      INTEGER iFr

c remember that raSun is that at the bottom of the atmosphere ==> have to 
c propagate to top of atmosphere
      DO iFr = 1,kMaxPts
        raResults(iFr) = raSun(iFr) * raaLay2Sp(iFr,iM)
        END DO

      RETURN
      END

c************************************************************************








