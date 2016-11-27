c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c******** THIS FILE CONTAINS VARIOUS USEFUL SUBROUTINES/FUNCTIONS *******
c** such as sorting, setting vertical temperature profiles, checking ****
c** kcarta.param, checking comp.param and xsec.param, splines etc *******
c************************************************************************
c check for isnan
      LOGICAL FUNCTION isnan_real(x)

      REAL x
      LOGICAL isnan

      print *,x,x+1.0e0
      IF (x+1.0e0 .eq. x) THEN
        isnan = .true.
      ELSE
        isnan = .false.
      ENDIF
      isnan = .true.

      IF (x .ge. 0.0) THEN
        print *,'gt',x
        isnan = .false.
      ELSEIF (x .le. 0.0) THEN
        print *,'lt',x
        isnan = .false.
      ENDIF

      isnan_real = isnan
      RETURN
      END

c************************************************************************
      LOGICAL FUNCTION isnan_double(x)

      logical isnan
      double precision x

      print *,x,x+1.0d0
      IF (x+1.0d0 .eq. x) THEN
        isnan = .true.
      ELSE
        isnan = .false.
      ENDIF

      isnan_double = isnan
      RETURN
      END

c************************************************************************
c this subroutine finds the tropopause by looking for the first cold point
c modelled on Scott Hannon's code tropopause_rtp.m which looks for the
c layer within the 50-400 mb range which has the lowest temp 
      INTEGER FUNCTION find_tropopause(raTemp,raPress,iaRadLayer,iNumLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params      
      REAL raTemp(kMixFilRows)         !! temperature structure
      REAL raPress(kProfLayer+1) !! pressure levels
      INTEGER iaRadLayer(kProfLayer)   !! which are the radiating layers
      INTEGER iNumLayer

c local vars
      REAL raT(kMixFilRows),rX,rJunk
      INTEGER iI,iL,i400mb,i50mb,iJL,iJ

c if        iaRadLayer = 001..100, everything ok
c but if eg iaRadLayer = 201..300, everything not ok with raPress
      iI  = 1
      iL  = iaRadLayer(iI)
      iI  = iNumLayer
      iJL = iaRadLayer(iI)
      iL = max(iL,iJL)
      iJ = 0
      IF (iL .GT. kProfLayer) THEN
        iJ = 1
 4      CONTINUE
        IF ((iJ+1)*kProfLayer .GE. iL) THEN
          GOTO 5
        ELSE
          iJ = iJ + 1
          GOTO 4
        END IF
      END IF

 5    CONTINUE
      iJ = iJ*kProfLayer

      DO iI = 1,kProfLayer
        raT(iI) = 0.0
      END DO

      DO iI = 1,iNumLayer
        iL = iaRadLayer(iI)
        raT(iL) = raTemp(iL)   !!note storing into raT(iL) instead of raT(iI)
c        print *,iI,iL,raPress(iL),raT(iL)
      END DO

      !! find i400mb
      rJunk = 1.0e10
      DO iI = 1,iNumLayer
        iL = iaRadLayer(iI)
        iJL = iL - iJ
        rX = abs(400.0 - raPress(iJL))
        IF (rX .LE. rJunk) THEN
          i400mb = iL
          rJunk = rX
        END IF
      END DO

      !! find i50mb
      rJunk = 1.0e10
      DO iI = 1,iNumLayer
        iL = iaRadLayer(iI)
        iJL = iL - iJ
        rX = abs(50.0 - raPress(iJL))
        IF (rX .LE. rJunk) THEN
          i50mb = iL
          rJunk = rX
        END IF
      END DO
      
      !! now look for bottom layer within [i400mb,i50mb] with lowest cold point
      iL = i50mb+1
      rJunk = raT(iL)
      DO iI = i50mb,i400mb,-1
        IF (raT(iI) .LE. rJunk) THEN
          rJunk = raT(iI)
          iL = iI
        END IF
      END DO

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Look for tropopause within AIRS preslays',i50mb,i400mb
      write(kStdWarn,*) 'Found it at AIRS presslayer ',iL
!      find_tropopause = iL

      !! may need to map this back to iaRadLayer
      DO iI = 1,iNumLayer
        IF (iaRadLayer(iI) .EQ. iL) THEN
          GOTO 10
        END IF
      END DO

 10   CONTINUE
      write(kStdWarn,*) 'this is in atmosphere layer (iaRadLayer) ',iI
      find_tropopause = iI
      
      RETURN
      END

c************************************************************************
c this subroutine sets the kCompressed database uncompression type
      SUBROUTINE SetSplineType(raPresslevels,iProfileLayers,
     $           iNumNLTEGases,iNLTE_SlowORFast,iSplineType)
      
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iplev 
      include '../INCLUDE/KCARTA_database.param'

c input vars      
      REAL raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iNumNLTEGases,iNLTE_SlowORFast
c ouput vars
c iSplinetype = 2 if FAST interp can be done
      INTEGER iSplineType

c local vars
      INTEGER iDefault,iJ,iUpDn,iMax
      REAL rDelta,rMin,rMax,rP,raFracDeltaP(kMaxLayer)
      REAL raFracDeltaPav(kMaxLayer),rAv1,rAv2,r1,r2,rMaxAv

c this assumes arbitrary pressure layering; can be slow, as kCARTA needs to 
c interpolate the (pressure levels, abs coeffs) of the database onto the klayers
c pressure layers, before doing the temperature interpolation
c note abs coeff = stored optical depth/default gas amount
      iDefault    = +1        !!! Spline  .... DEFAULT
      iSplineType = -1        !!! Linear
      iSplineType = -2        !!! Matlab uncompression (linear weights)
      iSplineType = +2        !!! Matlab uncompression (linear weights)
      iSplineType = +1        !!! Spline  .... DEFAULT
      iSplineType = iaaOverrideDefault(1,2)
      IF ((abs(iSPlineType) .NE. 1) .AND. (abs(iSPlineType) .NE. 2)) THEN
        write(kStdErr,*) 'invalid iSplineType = ',iSplineType
	CALL DoStop
      END IF      
      IF ((iSplineType .NE. iDefault) .AND. (kOuterLoop .EQ. 1)) THEN
        write(kStdErr,*),'iSplineType,iDefault = ',iSplineType,iDefault
        write(kStdWarn,*),'iSplineType,iDefault = ',iSplineType,iDefault	
      END IF

c recall kMaxLayer = 100 = the kCARTA adtabase
c        kProfLayer = N  = what klayers was compiled for (hopefully)
c      IF ((iProfileLayers .LE. kMaxLayer) .AND. 
c     $     ((iNumNLTEGases .LE. 0) .OR. (iNLTE_SlowORFast .LT. 0))
c     $     .AND. (kMaxLayer .LE. kProfLayer)) THEN
      IF ((iProfileLayers .LE. kMaxLayer) .AND. 
     $     ((iNumNLTEGases .LE. 0) .OR. (iNLTE_SlowORFast .EQ. -1))
     $     .AND. (kMaxLayer .LE. kProfLayer)) THEN
        !!quite possible that the pressure levels are same as kCARTA database
        !!in which case kcoeffSPL and kcoeffSPLandJAC can be sped up, as we can
        !!straightaway do tempr interpolation of stored abs coeffs without 
        !!having to worry about doing the pressure interpolation as well
        !!note abs coeff = stored optical depth/default gas amount
        iMax = 1
        rP = PLEV_KCARTADATABASE_AIRS(1)
        rDelta = 0.0
        rMin   = +1000.0
        rMax   = -1000.0
        rMaxAv = -1000.0
	write(kStdWarn,*) 'computing difference between default pavg and input layer average'
        write(kStdWarn,*) 'iJ p(database) p(klayers) frac_delta +++ pav(database) pav(klayers) frac_delta'
	write(kStdWarn,*) '------------------------------------------------------------------------------'
        DO iJ = (kMaxLayer-iProfileLayers+1),kMaxLayer
           raFracDeltaP(iJ) = 
     $              abs(raPresslevels(iJ)-PLEV_KCARTADATABASE_AIRS(iJ))/
     $                        PLEV_KCARTADATABASE_AIRS(iJ)
           r1 = PLEV_KCARTADATABASE_AIRS(iJ)-PLEV_KCARTADATABASE_AIRS(iJ+1)
           r2 =log(PLEV_KCARTADATABASE_AIRS(iJ)/PLEV_KCARTADATABASE_AIRS(iJ+1))
           rAv1 = r1/r2
           r1 = raPresslevels(iJ)-raPresslevels(iJ+1)
           r2 = log(raPresslevels(iJ)/raPresslevels(iJ+1))
           rAv2 = r1/r2
           raFracDeltaPav(iJ) = abs(rAv1-rAv2)/rAv1
           write(kStdWarn,100) iJ,PLEV_KCARTADATABASE_AIRS(iJ),raPresslevels(iJ),raFracDeltaP(iJ),rAv1,rAv2,raFracDeltaPav(iJ) 
          rDelta = rDelta + raFracDeltaP(iJ)
          IF (rMin .GT. raFracDeltaP(iJ)) THEN
            rMin = raFracDeltaP(iJ)
          END IF
          IF (rMax .LT. raFracDeltaP(iJ)) THEN
            rMax = raFracDeltaP(iJ)
            iMax = iJ
            rP = PLEV_KCARTADATABASE_AIRS(iJ)
          END IF
          IF (rMaxAv .LT. raFracDeltaPav(iJ)) THEN
            rMaxAv = raFracDeltaPav(iJ)
          END IF
        END DO
        rDelta = rDelta/iProfileLayers
        write(kStdWarn,*) 'difference between kCARTA database and klayers ...'
        IF ((rMin .LE. 1.0e-5) .AND. (rMax .LE. 5.0e-4) .AND. 
     $       (rMaxAv .LE. 5.0e-3)) THEN
          write(kStdWarn,*) '  can use kCARTA database levels for lower atm'
          iSplineType = iSplineType * 2  !!so -1 becomes -2, or +1 becomes +2
        ELSE
          write(kStdWarn,*) '  intrp kCARTA databse lvls for lower part of atm'
        END IF
        write(kStdWarn,*) 'MinDiff, MaxDiff, i(MaxDiff), Press(i(MaxDiff)), Sum(abs(diff))/Numlayers = '
        write(kStdWarn,*) rMin,rMax,iMax,rP,rDelta
        IF (abs(iSplineType) .EQ. 1) THEN
          write(kStdWarn,*) 'iSplineType = ',iSplineType, ' slow : interpolate'
          write(kStdWarn,*) 'database (pavglayers,abscoeeff) onto new plevels'
          write(kStdWarn,*) 'before doing the temp interpolation'
          write(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
        ELSEIF (abs(iSplineType) .EQ. 2) THEN
          write(kStdWarn,*) 'iSplineType = ',iSplineType, ' fast : use '
          write(kStdWarn,*) 'database (pavglayers,abscoeeff) for temp interp'
          write(kStdWarn,*) 'as presslevel scheme of klayers = kCARTA database!!'
          write(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
        ELSE
          write(kStdErr,*) 'need abs(iSplineType) = 1 or 2'
          CALL DoStop
        END IF

      END IF
 100  FORMAT(I3,' ',F10.5,' ',F10.5,' ',E10.5,' +++ ',F10.5,' ',F10.5,' ',E10.5)

      RETURN
      END

c************************************************************************
c this subroutine checks to see if the CO2 ppmv is ok
      SUBROUTINE check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix,
     $                         iaGases,rCO2mult)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
c raaMix     = mixing table info from *MIXFIL 
c iGas       = gas # iGas of iNumGases being added to the cumulative raaSum 
c iIpmix     = which of the mixed paths are being considered 
      REAL raaMix(kMixFilRows,kGasStore) 
      INTEGER iaGases(kMaxGas)
      REAL raaAmt(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
c output vars
      REAL rCO2mult       

c local vars

c these are the individual reference profiles, at kMaxLayer layers
c these are what is in kOrigRefPath
      REAL raR100Amt0(kMaxLayer),raR100Temp0(kMaxLayer)
      REAL raR100PartPress0(kMaxLayer),raR100Press0(kMaxLayer)
c these are the individual reference profiles, at kMaxLayer layers
c these are what is in kCO2ppmvFile
      REAL raR100Amt1(kMaxLayer),raR100Temp1(kMaxLayer)
      REAL raR100PartPress1(kMaxLayer),raR100Press1(kMaxLayer)

      CHARACTER*80 caFName
      REAL rMeanT,rMeanA,rMeanP,rMeanPP,rDirectPPMV0,rDirectPPMV1
      REAL rCO2ppmv
      INTEGER iI,iJ,iGasID,iError,iIPMIX,strfind

c find the weight 
      iIPMIX = 1
      iGasID = 2

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Checking the CO2 ppmv in the database profiles'

      IF (kCO2_UMBCorHARTMAN .EQ. -1) THEN
        write(kStdWarn,*) 'Using HARTMANN/NIRO CO2 linemixing, so NO chi fcns'
      ELSEIF (kCO2_UMBCorHARTMAN .EQ. +1) THEN
        write(kStdWarn,*) 'Using Strow/Tobin/Machado CO2 linemixing, need chi fcns'
      ELSE
        write(kStdErr,*) 'kCO2_UMBCorHARTMAN needs to be +/- 1'
        CALL DoStop
      END IF
      IF ((strfind(kCO2Path,'lblrtm') .EQ. 1) .OR. (strfind(kCO2Path,'LBLRTM') .EQ. 1)) THEN
        IF (kCO2_UMBCorHARTMAN .EQ. +1) THEN
	  write(kStdErr,*) 'kCO2Path has LBLRTM/lblrtm in it, but kCO2_UMBCorHARTMAN = -1'
          CALL DoStop
	END IF
      END IF
      
      IF (iaGases(iGasID) .NE. iGasID) THEN
        write(kStdErr,*) 'assumed iaGases(2) =  2, guess not'
        CALL DoStop
      END IF

      !!!input pressures are in atm; lowest layers (ignore) have P = 1000.0
      !!!                            layer before STtart of Atm have P = 0.0
      !!!                            rest of layers have meaningful P
      DO iI = 1,kProfLayer
        IF ((raaPress(iI,iGasID) .GT. 0.0) .AND. 
     $      (raaPress(iI,iGasID) .LT. 1.5)) THEN
          GOTO 10
        END IF
      END DO
 10   CONTINUE

      rCO2ppmv = 0.0
      DO iJ = iI,kProfLayer
        rMeanA   = raaPartPress(iJ,iGasID)/raaPress(iJ,iGasID) *1e6
        rCO2ppmv = max(rMeanA,rCO2ppmv)
      END DO
      write(kStdWarn,*) 'max rCO2ppmv from input profile = ',rCO2ppmv
      write(kStdWarn,*) 'kCARTA compiled for database CO2ppmv = ',kCO2ppmv

      IF (abs(rCO2ppmv-kCO2ppmv) .GE. 10) THEN
        write(kStdErr,*) 'input profile gasamts have max(CO2 ppmv) = ',rCO2ppmv
        write(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
        write(kStdErr,*) '  '
        write(kStdErr,*) '  If running NLTE SLOW suggest make a new kCARTA '
        write(kStdErr,*) '  compressed database with co2ppmv = ',rCO2ppmv
        write(kStdErr,*) '  '
        write(kStdErr,*) '  If running NLTE FAST code should be able to cope'
c        CALL DoStop
      END IF
      
      rCO2Mult = raaMix(iIpmix,iGasID) 

      !! open the true RefProf at kCO2ppmv
      caFName = kCO2ppmvFile
      write(kStdWarn,*) 'Reading CO2 Reference profile from ',caFName
      CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, 
     $         raR100Temp0,raR100Press0,raR100PartPress0,iError) 

      !! open the supposed RefProf
      CALL FindReferenceName(caFName,iGasID,-1) 
      CALL ReadRefProf(caFName,kMaxLayer,raR100Amt1, 
     $         raR100Temp1,raR100Press1,raR100PartPress1,iError) 

c -------------------------
      rMeanA = 0.0
      DO iI = 1,kMaxLayer
        rMeanA = rMeanA + abs(raR100Amt1(iI)/raR100Amt0(iI))
      END DO

      rMeanPP = 0.0
      DO iI = 1,kMaxLayer
        rMeanPP = rMeanPP + abs(raR100PartPress1(iI)/raR100PartPress0(iI))
      END DO

      rMeanP = 0.0
      DO iI = 1,kMaxLayer
        rMeanP = rMeanP + abs(raR100Press1(iI)/raR100Press0(iI))
      END DO

      rMeanT = 0.0
      DO iI = 1,kMaxLayer
        rMeanT = rMeanT + abs(raR100Temp1(iI)-raR100Temp0(iI))
      END DO
c -------------------------

      !! can directly figure out the ppmv in the "what we hope is true" file
      rDirectPPMV1 = 0.0
      DO iI = 1,kMaxLayer/2
        rDirectPPMV1 = rDirectPPMV1+abs(raR100PartPress1(iI)/raR100Press1(iI))
      END DO
      rDirectPPMV1 = rDirectPPMV1/(kMaxLayer/2)*1000000

      !! can directly figure out the ppmv in the "true" file
      rDirectPPMV0 = 0.0
      DO iI = 1,kMaxLayer/2
        rDirectPPMV0 = rDirectPPMV0+abs(raR100PartPress0(iI)/raR100Press0(iI))
      END DO
      rDirectPPMV0 = rDirectPPMV0/(kMaxLayer/2)*1000000

      write(kStdWarn,*) 'Checking CO2 LTE ppmv ...'
      write(kStdWarn,*) 'kCO2ppmv from kcarta.param       = ',kCO2ppmv
      write (kStdWarn,*) '  mean(rMeanAmt Ratio)  = ',rMeanA/kMaxLayer
      write (kStdWarn,*) '  mean(rMeanP   Ratio) = ',rMeanP/kMaxLayer
      write (kStdWarn,*) '  mean(rMeanPP  Ratio) = ',rMeanPP/kMaxLayer
      write (kStdWarn,*) '  mean(rMean   deltaT) = ',rMeanT/kMaxLayer
      write(kStdWarn,*) 'rCO2ppmv from input TRUErefprof2         = ',rDirectPPMV0
      write(kStdWarn,*) 'rCO2ppmv from input refprof2             = ',rDirectPPMV1
      write(kStdWarn,*) 'rCO2Mult from raaMixTable from .nml file = ',rCO2Mult
      write(kStdWarn,*) '  '

      IF (abs(rCO2Mult-1) .GE. 0.01) THEN
        write(kStdErr,*) 'you have set rCO2Mult in mixtable = ',rCO2Mult
        write(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
        write(kStdErr,*) '  suggest not trying NLTE calcs with this mixratio'
        write(kStdErr,*) '  make a new kCARTA compressed database with'
        write(kStdErr,*) '  co2ppmv = ',rCO2Mult*kCO2ppmv
c        CALL DoStop
      END IF

 100  FORMAT(A25,A80)
 101  FORMAT(A65,F12.8) 
      IF ((rMeanA/kMaxLayer .LE. 0.9995) .OR.
     $    (rMeanA/kMaxLayer .GE. 1.0005)) THEN
	write(kStdErr,100) 'v0 : kCO2ppmvFile =   ', kCO2ppmvFile
        write(kStdErr,101) 'mean rCO2ppmv (raCO2pp/raP) v0 from input TRUErefprof2         = ',rDirectPPMV0
	write(kStdErr,100) 'v1 : ref CO2 profile = ', caFName
        write(kStdErr,101) 'mean rCO2ppmv (raCO2PP/raP) v1 from input refprof2             = ',rDirectPPMV1
        write(kStdErr,*) 'oops check_co2ppmv : rMeanA is bad',rMeanA/kMaxLayer
        CALL DoStop
      END IF
      IF ((rMeanP/kMaxLayer .LE. 0.9995) .OR.
     $    (rMeanP/kMaxLayer .GE. 1.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMeanP is bad',rMeanP/kMaxLayer
        CALL DoStop
      END IF
      IF  ((rMeanPP/kMaxLayer .LE. 0.9995) .OR. 
     $    (rMeanPP/kMaxLayer .GE. 1.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMeanPP is bad',rMeanPP/kMaxLayer
        CALL DoStop
      END IF
      IF ((rMeanT/kMaxLayer .LE. -0.0005) .OR.
     $   (rMeanT/kMaxLayer .GE. +0.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMean deltaT is bad',rMeanT/kMaxLayer
        CALL DoStop
      END IF

      !!! now check the NLTE weak background in Upper Layers
      write (kStdWarn,*) '    checking weak backgnd upper atm ppmvs'
      caFName = kCO2ppmvFileBackUA
      CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, 
     $         raR100Temp0,raR100Press0,raR100PartPress0,iError) 
      !! can directly figure out the ppmv in the "what we hope is true" file
      rDirectPPMV0 = 0.0
      DO iI = 1,kMaxLayer/20
        rDirectPPMV0 = rDirectPPMV0+abs(raR100PartPress0(iI)/raR100Press0(iI))
      END DO
      rDirectPPMV0 = rDirectPPMV0/(kMaxLayer/20)*1000000
      write(kStdWarn,*) 'Weak background Upper Layers PPMV = ',rDirectPPMV0
      IF (abs(rDirectPPMV0 - rDirectPPMV1) .GE. 0.25) THEN
        write(kStdErr,*) 'oops check_co2ppmv : Weak Background UA LTE is bad'
        CALL DoStop
      END IF

      !!! now check the NLTE weak background; the profile should be the SAME as
      !!! for the Standard Profile
      write (kStdWarn,*) '    checking weak backgnd usual atm ppmvs and amts'
      caFName = kCO2ppmvFileBack
      CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, 
     $         raR100Temp0,raR100Press0,raR100PartPress0,iError) 
      !! can directly figure out the ppmv in the "what we hope is true" file
      rDirectPPMV0 = 0.0
      rMeanA = 0.0
      DO iI = 1,kMaxLayer
        rDirectPPMV0 = rDirectPPMV0+abs(raR100PartPress0(iI)/raR100Press0(iI))
        rMeanA = rMeanA + abs(raR100Amt1(iI) - raR100Amt0(iI))
      END DO
      rDirectPPMV0 = rDirectPPMV0/(kMaxLayer)*1000000
      write(kStdWarn,*) 'Weak background Standard Layers PPMV = ',rDirectPPMV0
      IF (abs(rDirectPPMV0 - rDirectPPMV1) .GE. 0.1) THEN
        write(kStdErr,*) 'oops check_co2ppmv : Weak Background LTE is bad'  
        CALL DoStop
      END IF
      IF (rMeanA .GE. 0.1) THEN
        write(kStdErr,*) 'oops check_co2ppmv : Weak backgnd amts different from STD amts'  
        CALL DoStop
      END IF
      write(kStdWarn,*) '  '

      RETURN
      END

c************************************************************************
c this subroutine checks to see if the gasID is 1-36 or 101-102 or 103
      INTEGER FUNCTION MainGas(iGasID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iGasID

      INTEGER iI,i1,i2
       
      iI = -1        !assume this is a main gas (1-36, or 101-102)

      i1 = -1
      i2 = -1

      IF ((iGasID .GE. 1) .AND. (iGasID .LE. kGasComp)) THEN
        i1 = 1          !! main gas 
      ELSEIF (iGasID .EQ. kNewGasHi+1) THEN
        i1 = 1          !! heavy water
      ELSEIF ((iGasID .GE. kNewGasLo) .AND. (iGasID .LE. kNewGasHi)) THEN
        i2 = 1          !! water continuum
      END IF
      
      IF ((i1 .EQ. 1) .OR. (i2 .EQ. 1)) THEN
        iI = 1
      END IF

      IF ((i2 .EQ. 1) .AND. kCKD .LT. 0) THEN
        write(kStdWarn,*) 'Cannot have gases 101,102 with CKD turned off'
        write(kStdErr,*) 'Cannot have gases 101,102 with CKD turned off'
        Call DoSTOP
      END IF

      MainGas = iI

      RETURN
      END

c************************************************************************
c this suroutine sets up the current print options
      SUBROUTINE SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,iType,
     $      iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $      iNumGases,iNpmix,iaNumLayer,iAtm)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iNumGases    !number of gases
      INTEGER iNpMix       !number of mix paths
      INTEGER iAtm         !surrent atmosphere
      INTEGER iaNumLayer(kMaxAtm)   !number of layers in atmospheres

      INTEGER iOutNum      !which printing set this is (1..kMaxPrint)
      INTEGER iPrinter     !what type (1,2,3)
      INTEGER iAtmPr       !if gas spectra, which gas; if MP, which MP,
                           !if radiances, which atm
      INTEGER iNp          !how many to output eg 100 gas spectra, or 50 MPs
      INTEGER iaOp(kPathsOut)  !list of paths/MPs/rads to output
      INTEGER iType        !-10 if dumb, 1 if paths, 2 if MPs, 3 if rads
c this is the printing switch,atmosphere# to print,# of layers to print, 
c   list of layers/paths to print (limited to kProfLayer for now) , and the 
c   pressures at which things are output 
      INTEGER iaPrinter(kMaxPrint),iaNp(kMaxPrint) 
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint) 

      INTEGER iDummy
     
      iPrinter = iaPrinter(iOutNum) 
      iAtmPr = iaGPMPAtm(iOutNum) 
      iNp = iaNp(iOutNum) 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 1)) THEN  
        !output all paths for gas
        iNp = kProfLayer*iNumGases 
      END IF 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 2)) THEN  
        !output all MPs
        iNp = iNpmix 
      END IF 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 3)) THEN  
        !output all radiances for atmosphere
        iNp = iaNumlayer(iAtm)
      END IF 

      DO iDummy = 1,iNp 
        iaOp(iDummy) = iaaOp(iOutNum,iDummy) 
      END DO 

      RETURN
      END

c************************************************************************
c this sets the temperature variation after nm_radnce is read in
      SUBROUTINE SetkTemperVary(iTemperVary)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c input param
      INTEGER iTemperVary      !!! from namelist file

c local var
      INTEGER iConstOrVary

c this is TEMPERATURE variation in layer
c       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileUPWELL_LINEAR_IN_TAU
c       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileDNWELL_LINEAR_IN_TAU
c >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
c >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
c      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
c      kTemperVary = +1     !!!temperature in layer varies
c      kTemperVary = +2     !!!temperature in layer varies linearly, simple
c      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
c      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
c      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM,
c                           !!!  no O(tau) approx, very similar to kTemperVary=4
c      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM,
c                           !!!  debugged for small O(tau), used with EliMlawer 12/2015
c      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
c                           !!!  x/6 as x-->0 compared to kTemperVary = +42 *****
c      IF (kFlux .LE. 0) THEN
c        kTemperVary = -1     !!! temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
c      ELSE
c        kTemperVary = +43    !!! temperature in layer varies linearly, ala RRTM, LBLRTM, and has
c	                     !!! x/6 as x-->0 compared to kTemperVary = +42 ****
c      END IF

      kTemperVary = iTemperVary
      
      iConstOrVary = -1   !! if do flux, do linear vary T with tau
      iConstOrVary = +1   !! if only RaDTrans, then do constant T in layer, default SARTA/kCARTA for RT only
      
      IF (kFlux .LE. 0) THEN
        IF (iConstOrVary .GT. 0) THEN
          kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
  	  write(kStdWarn,*) 'kFlux .LE. 0 so set kTemperVary = -1'
        ELSEIF (iConstOrVary .LT. 0) THEN    
          kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
	                       !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
      	  write(kStdWarn,*) 'kFlux .LT. 0 but set kTemperVary = 43'	
        END IF
      ELSEIF (kFlux .GT. 0) THEN
        kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
	                     !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
	write(kStdWarn,*) 'kFlux .GT. 0 so set kTemperVary = 43'	
      END IF

      !!! new, do what the user wishes!!!
      IF ((kFlux .LE. 0) .AND. (iTemperVary .GT. 0)) THEN
        kTemperVary = +43
      END IF
      
      !!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
      ! kTemperVary = +43           
      !IF (iTemperVary .NE. kTemperVary) THEN
      !  write(kStdWarn,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
      !  write(kStdErr,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary	
      !  kTemperVary = iTemperVary      
      !END IF
      !!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
            
      IF (kTemperVary .EQ. -1) THEN
        write(kStdWarn,*) 'kTemperVary = -1     !layer temp constant (SARTA DEFAULT)'
      ELSEIF (kTemperVary .EQ. +1) THEN
        write(kStdWarn,*) 'kTemperVary = +1     !layer temp varies'
      ELSEIF (kTemperVary .EQ. +2) THEN
        write(kStdWarn,*) 'kTemperVary = +2     !layer temp varies linearly, simple v2'
      ELSEIF (kTemperVary .EQ. +3) THEN
        write(kStdWarn,*) 'kTemperVary = +3     !layer temp varies linearly, ala LBLRTM v3'
      ELSEIF (kTemperVary .EQ. +4) THEN
        write(kStdWarn,*) 'kTemperVary = +4     !layer temp varies linearly, ala LBLRTM v4 O(tau^2)'
      ELSEIF (kTemperVary .EQ. +41) THEN
        write(kStdWarn,*) 'kTemperVary = +41    !layer temp varies linearly, ala LBLRTM v4 (Pade)'
      ELSEIF (kTemperVary .EQ. +42) THEN
        write(kStdWarn,*) 'kTemperVary = +42    !layer temp varies linearly, ala LBLRTM v4 O(tau)'
      ELSEIF (kTemperVary .EQ. +43) THEN
        write(kStdWarn,*) 'kTemperVary = +43    !layer temp varies linearly, ala LBLRTM v4 O(tau) -> tau/6'
      ELSE
        write(kStdErr,*)'kTemperVary = ',kTemperVary,'unknown option'
        CALL DoStop
      END IF 

      iaaOverrideDefault(2,1) = kTemperVary
      
      RETURN
      END

c************************************************************************
c this subroutine does some more initializations
      SUBROUTINE SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iMixFileLines,iVertTempSet,iNpmix
      REAL raaMix(kMixFilRows,kGasStore)       !mixing table
 
      INTEGER iFileIDLO,iFileIDhi

c these are needd to be set if kRTP = -10 (user supplied levels profile) or -5,-6 (TAPE5/TAPE6 from LBLRTM)
c   1,2,3 = SurfPress/SurfTemp/SurfHgt
c   4,5   = ViewHgt and ViewAngle/direction
c   6     = Based on (4=ViewHgt) the code figures out highest pressure at which rad is to be output
      DO iFileIDLO = 1,10
        raRTP_TxtInput(iFileIDLO) = -9999
      END DO

c these are for seeing how cloud params are overwritten
c raaRTPCloudParams0(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   from rtpfile
c raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   after kcarta resets
c this gets set in rtp_interface.f
      DO iFileIDLO = 1,7
        raaRTPCloudParams0(1,iFileIDLO) = -1.0
	raaRTPCloudParamsF(1,iFileIDLO) = -1.0
        raaRTPCloudParams0(2,iFileIDLO) = -1.0
	raaRTPCloudParamsF(2,iFileIDLO) = -1.0
      END DO
      
c this is really for Mie scattering and VIS/UV ocean reflectance
      kSolAzi = 0.0
      kSatAzi = 0.0
      kWindSpeed = 0.0
      kLatitude = 0.0
      kMonth = 1.0
c this is to stop the code flux calcs at LBLRTM toa, default = -1 so keep doing calcs till 0.005 mb
c else if kLBLRTM_toa > 0 then the flux code
c   finds the highest layer whose pressure is greater than this
c   sets all ODS to 0 above this layer so essentially
c     rU = rUp0 exp(-k) + B(T) (1-exp(-k)) --> r0 for upwelling rad
c     rD = 0 since the downwelling rad is initialized with ttorad(f,2.7) ~ 0
      kLBLRTM_toa = -1.0
      
c assume no *mixfil section 
      iMixFileLines = -1 
 
c the vertical temperature profile has not been set yet 
      iVertTempSet = -1 
 
c initialize the mixing table to weights of 0.0 
      iNpMix = 1 
      DO iFileIDLo = 1,kMixFilRows 
        DO iFileIDHi = 1,kGasStore 
          raaMix(iFileIDLo,iFileIDHi) = 0.0 
        END DO 
      END DO 

c initialize kaaNumVectors      
      DO iFileIDLo = 1,kMaxGas 
        DO iFileIDHi = 1,kMaxLayer
          kaaNumVectors(iFileIDLo,iFileIDHi) = 0
        END DO 
      END DO 

      RETURN
      END

c************************************************************************
c this subroutine inits file unit numbers
      SUBROUTINE InitializeFileUnits 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c set up file unit numbers so at run time they default to STDIN,STDOUT,STDOUT  
      kStdDriver = 5 
      kStdkCarta = 6 
      kStdJacob  = 6  
 
c set up common block parameters indicating all units closed 
      kStdErrOpen  = +1    !! logical unit 0 
      kStdWarnOpen = -1 
 
      kStdDriverOpen   = -1 
      kStdkCartaOpen   = -1 
      kStdJacobOpen    = -1 
      kStdFluxOpen     = -1 
      kStdPlanckOpen   = -1 
 
      kCompUnitOpen    = -1 
      kProfileUnitOpen = -1 
 
      kTempUnitOpen = -1 

      kBloatPlanckOpen = -1
      kBloatOutOpen    = -1

      kStdPlanckUAOpen = -1
      kNLTEOutUAOpen   = -1

      kBloatPlanckUAOpen   = -1
      kBloatNLTEOutUAOpen  = -1

      RETURN
      END

c************************************************************************
c this subroutine opens the driver namelist file and figures out the name
c  of the error/warning log file
      SUBROUTINE ErrorLogName(caDriverName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c this is the driver file name from the command line arguments
      CHARACTER*80 caDriverName

      CHARACTER*30 namecomment

c this is for OUTPUT
c caLogFile     = name of success/warning log file 'warning.msg'
c caComment     = comment the user writes
c iOutTypes     = number of printing options specified
c iaPrinter     = for each option, which output type specified
c iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output 
c iaNp          = for each option, how many paths/MPs/layers to be output
c iaaOp         = for each option, list of paths/MP/layers to be output
c raaOp         = for option 3, list fract of layers used for radiance output
c raaUserPress  = for option 3, list of pressures for output radiances
c iNatm2        = number of atmospheres that *OUTPUT thinks there is
      INTEGER iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
      INTEGER iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
      INTEGER iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
      INTEGER iIOUN,iErr
      CHARACTER*80 caComment,caComment1
      CHARACTER*80 caLogFile,caLogFile1
      REAL raaOp(kMaxPrint,kProfLayer),raaOp1(kMaxPrint,kProfLayer)

      NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter,
     $           iaGPMPAtm,iaNp,iaaOp,raaOp

      iIOun=kStdDriver 
      IF (iIOUN .NE. 5) THEN 
        OPEN(UNIT = iIOun,FILE = caDriverName,STATUS='OLD',IOSTAT=iErr) 
        IF (iErr .NE. 0) THEN 
          write (kStdErr,*) 'in subroutine ErrorLogName, error reading'
          write (kStdErr,*) 'namelist file to find name of logfile ... '
          WRITE(kStdErr,1070) iErr, caDriverName 
 1070     FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80) 
          CALL DoSTOP 
        ENDIF 
      END IF 
      kStdDriverOpen = 1 

c      write(kStdWarn,*) 'grepping input nml file for caLogFile name'
c      print *,'translating x...'

      namecomment = '******* OUTPUT section *******'
      caLogFile = 'warning.msg'     !this is the default name
      read (iIOUN,nml = nm_output)
      caLogFile1  =  caLogFile

      kWarnFile  =  caLogFile

      close (iIOUN)
      kStdDriverOpen = -1 

      RETURN
      END

c************************************************************************
c this subroutine summarizs the output options
      SUBROUTINE SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp, 
     $                  raaOp,raaUserPress) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c this is the printing switch,atmosphere# to print,# of layers to print, 
c   list of layers/paths to print (limited to kProfLayer for now) , and the 
c   pressures at which things are output 
      INTEGER iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint) 
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint) 
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer) 

      INTEGER iDummy,iOutnum

      write(kStdWarn,*) '# of printing options selected = ',iOuttypes 
      write(kStdWarn,*) '     index     option type      atm#    numpaths' 
      write(kStdWarn,*) '------------------------------------------------' 
      DO iDummy = 1,iOuttypes 
        write(kStdWarn,30) iDummy,iaPrinter(iDummy),iaGPMPAtm(iDummy), 
     $                    iaNp(iDummy) 
        write(kStdWarn,*) 'paths to be printed : (if numpaths=-1,print all)' 
        write(kStdWarn,*)(iaaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy)) 
        IF (iaPrinter(iDummy) .EQ. 3) THEN 
          write(kStdWarn,*)(raaOp(iDummy,iOutNum), 
     $                       iOutNum=1,iaNp(iDummy)) 
          write(kStdWarn,*)(raaUserPress(iDummy,iOutNum), 
     $                       iOutNum=1,iaNp(iDummy)) 
        END IF 
        write(kStdWarn,*) '    ' 
      END DO
 
 30   FORMAT('     ',4(I3,'          '))

      RETURN
      END

c************************************************************************
c this subroutine checks the MixPath Vertical Temps
      SUBROUTINE CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp, 
     $                        iNpmix,iCO2,iaGases) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
 
      INTEGER iCo2              !which gas used to mimic CO2 temps
      INTEGER iNumGases         !how many gases
      INTEGER iNpMix            !number of mixed paths
      REAL raaTemp(kProfLayer,kGasStore)       !profile temp
      REAL raaMix(kMixFilRows,kGasStore)       !mixing table
      REAL raMixVertTemp(kMixFilRows)          !temperatures of MP layers
      INTEGER iaGases(kMaxGas)               !gasIDs stored in order

      INTEGER iDummy,iFileIDLo

      iCO2 = -1

      IF (iNpmix .GT. 0) THEN 
c search for the CO2 gas === gasID 2 
c since the gases have to be entered in ascending order, either gas1 or gas2 
c is CO2 
        iCO2 = -1 
        IF (iaGases(1) .EQ. 2) THEN 
          iCO2 = 1 
          write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!' 
        ELSE IF (iaGases(2) .EQ. 2) THEN 
          iCO2 = 2 
          write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!' 
        ELSE !!!for some strange reason, no CO2 present 
          iCO2 = 1 
          write(kStdWarn,*) 'Temperature of Gas number 1 will mimic CO2!!' 
      END IF 
 
        CALL GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, 
     $                      iNpmix,iCO2) 
 
        write(kStdWarn,*) 'Checking Mixed Path Temp' 
        iFileIDLo = 0 
        DO iDummy = 1,iNpmix 
          IF (raMixVertTemp(iDummy) .LT. 0.0) THEN 
            write(kStdWarn,*) 'Negative MP Temp in Mixed Path',iDummy 
            iFileIDLo = iFileIDLo+1 
          END IF                  
        END DO 
        IF (iFileIDLo .GT. 0) THEN 
          write(kStdWarn,*) 'Warning! negative MP temperatures found!' 
          write(kStdErr,*) 'Warning! negative MP temperatures found!' 
          CALL DoSTOP 
        END IF 
      END IF 

 1111 FORMAT(A1) 

      RETURN
      END

c************************************************************************
c this subroutine does the command line stuff
      SUBROUTINE DoCommandLine(iMicrosoft,caDriverName,caOutName,
     $                         caJacobFile,iOutFileName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caDriverName is the name of the driver file to be processed 
      CHARACTER*80 caDriverName 
c caOutName is the name of the unformatted output file name 
c integer iOutFileName tells whether or not there is a driver name, or 
c dump output to Unit 6 
      INTEGER iOutFileName 
      CHARACTER*80 caOutName 
c caJacobFile is the name of the unformatted output file name for Jacobians 
      CHARACTER*80 caJacobFile 
c this tells if we have MS Product ie no command line stuff!
      INTEGER iMicroSoft
c this is the number of args
      INTEGER iargc 

      INTEGER iDummy,iError

      IF (iMicroSoft .GT. 0) THEN 
         
        !no command line options .. do this! 
        print *,'Enter (1) if standard kcarta computation ' 
        print *,'Enter (2) if kcarta + jacobian computation ' 
        read *,iDummy 
        IF ((iDummy .GT. 2) .OR. (iDummy .LT. 1)) THEN  
          write(kStdErr,*) 'Microsoft user : please enter 1 or 2' 
          CALL DoSTOP 
        END IF 
       
        print *,'Enter driver namelist filename (enclose in quotes) : ' 
        read *,caDriverName 
        kStdDriver = kStdDriverKK 
 
        print *,'Enter output standard filename  (enclose in quotes) : ' 
        read *,caOutName 
        kStdkCarta = kStdkCartaKK 
        iOutFileName  =  1 
 
        IF (iDummy .EQ. 2) THEN 
          print *,'Enter output jacobian filename  (enclose in quotes) : ' 
          read *,caJacobFile 
          kStdJacob = kStdJacobKK 
        END IF 

        iMicrosoft = iDummy    !tells number of output files (w/o flux, planck)
       
       ELSE 
         !use command line stuff 
         iDummy = iargc() 
 
         IF (iDummy .GT. 3) THEN  
           write(kStdErr,*) 'more than three arguments in command line' 
           write(kStdErr,*) 'is NOT allowed' 
           CALL DoSTOP 
         END IF 

         iOutFileName = -1         !assume no name 
         DO iError = 1,iDummy 
           IF (iError .EQ. 1) THEN 
             CALL getarg(1,caDriverName) 
             IF (caDriverName(1:1) .NE. '-') THEN 
               kStdDriver = kStdDriverKK 
             END IF 
           END IF 
 
           IF (iError .EQ. 2) THEN 
             CALL getarg(2,caOutName) 
             IF (caOutName(1:1) .NE. '-') THEN 
               iOutFileName = 1 
               kStdkCarta = kStdkCartaKK 
             END IF 
           END IF 
 
           IF (iError .EQ. 3) THEN 
             CALL getarg(3,caJacobFile) 
             IF (caJacobFile(1:1) .NE. '-') THEN 
               kStdJacob = kStdJacobKK 
             END IF 
           END IF 
         END DO 
 
        iMicrosoft = iDummy-1  !tells number of output files (w/o flux,planck)
      END IF 

      RETURN
      END 
 
c************************************************************************
c this subroutine stores the reference gas amts/temps etc
      SUBROUTINE StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,
     $   raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c this is the gasnumber (not gasID!!!!!!!!!!!!!!)
      INTEGER iGas
      INTEGER iaGases(kMaxGas)
c these are the individual reference profiles 
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer) 
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c these are the reference profiles stored in matrices 
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore) 
      REAL raaRPress(kProfLayer,kGasStore) 
      REAL raaRPartPress(kProfLayer,kGasStore) 

      INTEGER iI

      DO iI = 1,kProfLayer 
        raaRAmt(iI,iGas) = raRAmt(iI)             !amts 
        raaRTemp(iI,iGas) = raRTemp(iI)           !temps 
        raaRPress(iI,iGas) = raRPress(iI)         !press 
        raaRPartPress(iI,iGas) = raRPartPress(iI) !part press 

        IF (raRAmt(iI) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error in Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iI,' has negative amount',raRAmt(iI) 
          CALL DoStop 
        END IF 
 
        IF (raRTemp(iI) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error in Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iI,' has negative tempr',raRTemp(iI) 
          CALL DoStop 
        END IF 
 
        IF (raRPress(iI) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error in Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iI,' has negative pressure',raRPress(iI) 
          CALL DoStop 
        END IF 
 
        IF (raRPartPress(iI) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error in Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iI,' has negative part press',
     $ raRPartPress(iI) 
          CALL DoStop 
        END IF 
      END DO 

      RETURN
      END

c************************************************************************
c this subroutine sets the reference gas amts/temps etc
      SUBROUTINE SetReference(raRAmt,raRTemp,raRPress,raRPartPress, 
     $   raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c this is the gasnumber (not gasID!!!!!!!!!!!!!!)
      INTEGER iGas
c these are the individual reference profiles 
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer) 
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c these are the reference profiles stored in matrices 
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore) 
      REAL raaRPress(kProfLayer,kGasStore) 
      REAL raaRPartPress(kProfLayer,kGasStore) 

      INTEGER iInt

      DO iInt = 1,kProfLayer 
        raRAmt(iInt) = raaRAmt(iInt,iGas) 
        raRTemp(iInt) = raaRTemp(iInt,iGas) 
        raRPress(iInt) = raaRPress(iInt,iGas) 
        raRPartPress(iInt) = raaRPartPress(iInt,iGas) 
      END DO 

      RETURN
      END 

c************************************************************************
c this subroutine close units
      SUBROUTINE TheEnd(iaGases,iNumGases,iaList,raFiles)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iNumGases,iaGases(kMaxGas),iaList(kNumkCompT)
      REAL raFiles(kNumkCompT)
      INTEGER iI,iJ,iG,iaSum(kMaxGas),iaCount(kMaxGas),iSum,iaJunk(kMaxGas)

      write(kStdWarn,*) '*****************************************'      
      write(kStdWarn,*) 'kComp stats ..............'
      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Freq chunks processed'
      write(kStdWarn,*) (raFiles(iaList(iJ)),iJ=1,kOuterLoop)

      write(kStdWarn,*) 'Stats of compressed vecs per gas per chunk'
      DO iI = 1,kMaxGas
        iaSum(kMaxGas) = 0
      END DO
      DO iI = 1,iNumGases
        iG = iaGases(iI)
        write(kStdWarn,*) 'num compressed vector stats for gas #, gasID ',iI,iG
        write(kStdWarn,*) (kaaNumVectors(iG,iJ),iJ=1,kOuterLoop)
        write(kStdWarn,*) ' '
        DO iJ=1,kOuterLoop
          IF (kaaNumVectors(iG,iJ) .GT. 0) THEN
            iaSum(iG) = iaSum(iG) + kaaNumVectors(iG,iJ)
            iaCount(iG) = iaCount(iG) + 1
          END IF
        END DO
      END DO

      write(kStdWarn,*) 'gas#   GASID   NumVecs   NumChunks  AvgVecs'
      DO iI = 1,iNumGases
        iG = iaGases(iI)
        IF (iaCount(iG) .GT. 0) THEN
          write(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG),
     $                      1.0*iaSum(iG)/iaCount(iG)
        ELSE
          write(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG),0.0
        END IF
      END DO
 15   FORMAT(2(I5,' '),'   ',I5,'      ',I5,'     ',F8.4)

      write(kStdWarn,*) ' '
      write(kStdWarn,*) '  Chunk   NumGases'
      DO iJ = 1,kOuterLoop
        iSum = 0
        DO iI = 1,iNumGases
          iG = iaGases(iI)
          IF (kaaNumVectors(iG,iJ) .GT. 0) iSum = iSum + 1
        END DO
        write(kStdWarn,16) raFiles(iaList(iJ)),iSum
      END DO
  16  FORMAT(F9.2,'    ',I3)

      write(kStdWarn,*) ' '
      DO iJ = 1,kOuterLoop
        iSum = 0
        DO iI = 1,iNumGases
          iG = iaGases(iI)
          IF (kaaNumVectors(iG,iJ) .GT. 0) THEN
            iSum = iSum + 1
            iaJunk(iSum) = iG
          END IF
        END DO
        write(kStdWarn,*) '  Chunk = ',raFiles(iaList(iJ)), ' numgases = ',iSum,' gasIDs are .... '
        write(kStdWarn,*) (iaJunk(iI),iI=1,iSum)
      END DO
  17  FORMAT(I3,' ')

      IF (kStdkCarta .NE. 6) THEN  
        write(kStdWarn,*)'closing binary output file'  
        CLOSE(UNIT = kStdkCarta)      !close file where kCARTA binary goes to  
        kStdkCartaOpen  =  -1 
      END IF  
  
      IF ((kJacobian .GT. 0) .AND. (kStdJacob .NE. 6)) THEN  
        write(kStdWarn,*)'closing jacobian binary file'  
        CLOSE(UNIT = kStdJacob)      !close file where Jacob binary goes to  
        kStdJacobOpen  =  -1 
      END IF  

      IF (kJacobian .GT. 0) THEN
        IF ((kStdJacobOpen .EQ. 1)  .AND. (kStdJacob .NE. 6)) THEN
          write(kStdWarn,*)'closing jacobian binary file'
          CLOSE(UNIT = kStdJacob)       !close file where Jacob binary goes to
        END IF
        IF (kStdJacob2Open .EQ. 1) THEN
          write(kStdWarn,*)'closing jacobian2 (column) binary file'
          CLOSE(UNIT = kStdJacob2)       !close file where Jacob binary goes to
        END IF
      END IF

c      IF (kJacobian .GT. 0) THEN  
c        write(kStdWarn,*)'closing jacobian2 (column) binary file'  
c        CLOSE(UNIT = kStdJacob2)      !close file where Jacob binary goes to  
c        kStdJacob2Open  =  -1 
c      END IF  
 
      IF (kFlux .GT. 0) THEN 
        write(kStdWarn,*)'closing flux binary file'  
        CLOSE(UNIT = kStdFlux)       !close file where flux binary goes to  
        kStdFluxOpen  =  -1 
      END IF  

      IF (kStdPlanckOpen .GT. 0) THEN 
        write(kStdWarn,*)'closing planck binary file'  
        CLOSE(UNIT = kStdPlanck)    !close file where planck binary goes to
        kStdPlanckOpen  =  -1 
      END IF  

      IF (kBloatPlanckOpen .GT. 0) THEN 
        write(kStdWarn,*)'closing bloated planck binary file'  
        CLOSE(UNIT = kBloatNLTEPlanck)    !close file 
        kBloatPlanckOpen  =  -1 
      END IF  

      IF (kBloatOutOpen .GT. 0) THEN 
        write(kStdWarn,*)'closing bloated binary file'  
        CLOSE(UNIT = kBloatNLTEOut)    !close file 
        kBloatOutOpen  =  -1 
      END IF  

      IF (kStdPlanckUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing UA planck binary file'  
        CLOSE(UNIT = kStdPlanckUAOpen)      !close file 
        kStdPlanckUAOpen = -1 
      END IF  

      IF (kNLTEOutUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing UA binary file'  
        CLOSE(UNIT = kNLTEOutUAOpen)      !close file 
        kNLTEOutUAOpen = -1 
      END IF  

      IF (kBloatPlanckUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloat UA planck binary file'  
        CLOSE(UNIT = kBloatPlanckUAOpen)      !close file 
        kBloatPlanckUAOpen = -1 
      END IF  

      IF (kBloatNLTEOutUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloat UA binary file'  
        CLOSE(UNIT = kBloatNLTEOutUAOpen)      !close file 
        kBloatNLTEOutUAOpen = -1 
      END IF  

      write(kStdWarn,*) 'closed files ... !!!!!!!!!!!'
  
      RETURN
      END

c***********************************************************************
c this subroutine closes all files in case of an emergency stop
c assumes the message ends with '$'
      SUBROUTINE DoSTOPMesg(caMessage)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iI,iFound
      CHARACTER     caMessage*(*)
      CHARACTER*120 caMessage2

      DO iI = 1,80
        caMessage2(iI:iI) = ' '
      END DO

      iI = 80
      iI = len(caMessage)
      IF (iI .GT. 120) THEN
        write(kStdErr,*) 'lengthh of error message is over 120 characters!'
	write(kStdErr,*) caMessage
        CALL DoStop
      END IF

 5    CONTINUE
      IF ((caMessage(iI:iI) .NE. '$') .AND. (iI .GT. 1)) THEN
        iI = iI - 1
        GOTO 5
      END IF
 
      IF (iI .LE. 1) THEN
        write(kStdErr,*) 'caMessage needs "$" to end '
        CALL DoStop
      END IF

c      write(kStdErr,*) 'length of caMessage = ',iI
      caMessage2(1:iI-1) = caMessage(1:iI-1)

      write(kStdErr,10) caMessage2
      CALL DoStop

 10   FORMAT(A120)

      RETURN
      END

c***********************************************************************
c this subroutine closes all files in case of an emergency stop
      SUBROUTINE DoSTOP

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      write(kStdWarn,*)'Fatal Error found : closing all units ..'

      IF ((kStdDriverOpen .EQ. 1) .AND. (kStdDriver .NE. 5)) THEN
        write(kStdWarn,*)'closing driver file'
        CLOSE(UNIT = kStdDriver)          !close driver file 
      END IF

      IF ((kStdkCartaOpen .EQ. 1) .AND. (kStdkCarta .NE. 6)) THEN
        write(kStdWarn,*)'closing binary output file'
        CLOSE(UNIT = kStdkCarta)        !close file where kCARTA binary goes to
      END IF

      IF (kCompUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing kcomp/xsec file'
        CLOSE(UNIT = kCompUnit)         !close kCompressed file/xsec data file
      END IF

      IF (kJacobian .GT. 0) THEN
        IF ((kStdJacobOpen .EQ. 1)  .AND. (kStdJacob .NE. 6)) THEN
          write(kStdWarn,*)'closing jacobian binary file'
          CLOSE(UNIT = kStdJacob)       !close file where Jacob binary goes to
        END IF
        IF (kStdJacob2Open .EQ. 1) THEN
          write(kStdWarn,*)'closing jacobian2 (column) binary file'
          CLOSE(UNIT = kStdJacob2)       !close file where Jacob binary goes to
        END IF
      END IF

      IF (kFlux .GT. 0) THEN
        write(kStdWarn,*)'closing flux binary file'
        CLOSE(UNIT = kStdFlux)         !close file where flux binary goes to
      END IF

      IF (kStdPlanckOpen .GT. 0) THEN
        write(kStdWarn,*)'closing planck binary file'
        CLOSE(UNIT = kStdPlanck)        !close file where planck binary goes to
      END IF

      IF (kProfileUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing profile file '
        CLOSE(UNIT = kProfileUnit)       !close profile file 
      END IF

      IF (kTempUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing temporary param file'
        CLOSE(UNIT = kTempUnit)          !close temporary file eg comp.param
      END IF

      IF (kBloatPlanckOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloated planck binary file'  
        CLOSE(UNIT = kBloatNLTEPlanck)      !close file 
        kBloatOutOpen = -1 
      END IF  

      IF (kBloatOutOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloated binary file'  
        CLOSE(UNIT = kBloatNLTEOut)      !close file 
        kBloatOutOpen = -1 
      END IF  

      IF (kStdPlanckUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing UA planck binary file'  
        CLOSE(UNIT = kStdPlanckUA)      !close file 
        kStdPlanckUAOpen = -1 
      END IF  

      IF (kNLTEOutUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing UA binary file'  
        CLOSE(UNIT = kNLTEOutUA)      !close file 
        kNLTEOutUAOpen = -1 
      END IF  
  
      IF (kBloatPlanckUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloat UA planck binary file'  
        CLOSE(UNIT = kBloatPlanckUAOpen)      !close file 
        kBloatPlanckUAOpen = -1 
      END IF  

      IF (kBloatNLTEOutUAOpen .EQ. 1) THEN 
        write(kStdWarn,*)'closing bloat UA binary file'  
        CLOSE(UNIT = kBloatNLTEOutUAOpen)      !close file 
        kBloatNLTEOutUAOpen = -1 
      END IF  

      write(kStdErr,*) 'bad luck ... emergency exit!'
      write(kStdWarn,*) 'bad luck ... emergency exit!'

      CLOSE(UNIT = kStdErr)             !close error log
      CLOSE(UNIT = kStdWarn)            !close warning log
 
      call exit(1)                    !sad exit so return +1

      STOP

      RETURN
      END

c************************************************************************
c this function, depending on iNp, calls a binary or a sequential search
c to find iLay in iaOp
      INTEGER FUNCTION DoOutputLayer(iLay,iNp,iaOp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLay  = layer number to be looked for
c iaOp  = array containing list of layers
c iNp   = search indices 1..iNp of iaOp, to look for iLay
      INTEGER iLay,iNp,iaOp(*)

c integer functions that do the search
      INTEGER BinarySearch,SequentialSearch

      IF (iNp .LT. 16) THEN
        DoOutputLayer = SequentialSearch(iLay,iNp,iaOp)
      ELSE
        DoOutputLayer = BinarySearch(iLay,iNp,iaOp)
      END IF

      RETURN 
      END 

c************************************************************************
c this function checks to see if current GasID should have its d/dq saved
c if it does, the function result is WHICH gas it is in the *JACOBN wishlist
c else the function result = -1
      INTEGER FUNCTION DoGasJacob(iGasID,iaJacob,iJacob)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = current gasID
c iaJacob  = list of GasID's whose d/dq we want to output
c iJacob   = number of GasID's whose d/dq we want to output
      INTEGER iGasID,iJacob,iaJacob(kMaxDQ)

      INTEGER iI,iFound,iAns

      iFound = -1
      iAns = -1
      iI = 1

 15   CONTINUE
      IF ((iFound .LT. 0) .AND. (iI .LE. iJacob)) THEN
c check to see if iGasID is in iaJacob
        IF (iGasID .EQ. iaJacob(iI)) THEN
          iFound = 1
          iAns = iI
        ELSE 
          iI = iI+1
          GO TO 15
          END IF
      END IF
            
      DoGasJacob = iAns

      RETURN
      END

c************************************************************************
c this function checks to which posn GasID is in, in iaGases
c it mimics the "ismember" function in Matlab
      INTEGER FUNCTION WhichGasPosn(iGasID,iaGases,iNumGases)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = current gasID
c iaGases  = list of GasID's that are being used
c iJacob   = number of GasID's that are being used
      INTEGER iGasID,iNumGases,iaGases(kMaxGas)

      INTEGER iI,iFound,iAns

      iFound = -1
      iAns = -1
      iI = 1

 15   CONTINUE
      IF ((iFound .LT. 0) .AND. (iI .LE. iNumGases)) THEN
c check to see if iGasID is in iaGases
        IF (iGasID .EQ. iaGases(iI)) THEN
          iFound = 1
          iAns = iI
        ELSE 
          iI = iI+1
          GO TO 15
          END IF
      END IF
            
      IF ((iGasID .EQ. 101) .OR. (iGasID .EQ. 102)) THEN
        iAns  =  1
      END IF

      WhichGasPosn = iAns

      RETURN
      END

c************************************************************************
c this subroutine initializes all the rows of the 
c (REAL) array of absorption coeffs
      SUBROUTINE InitializeReal(raaAb)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL raaAb(kMaxPts,kProfLayer)

      INTEGER iLay,iFreq

      DO iLay = 1,kProfLayer
        DO iFreq = 1,kMaxPts
          raaAb(iFreq,iLay) = 0.0
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes all the rows of the 
c (REAL) array of mixed paths
      SUBROUTINE InitializeRealMP(raaAb,iNpmix)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL raaAb(kMaxPts,kMixFilRows)
      INTEGER iNpMix

      INTEGER iLay,iFreq

      DO iLay = 1,iNpMix      !note : initialize only wot is necessary
        DO iFreq = 1,kMaxPts
          raaAb(iFreq,iLay) = 0.0
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes the (DOUBLE) array of absorption coeffs
      SUBROUTINE initialize(daaAb)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      DOUBLE PRECISION daaAb(kMaxPts,kProfLayer)

      INTEGER iLay,iFreq

      DO iLay = 1,kProfLayer
        DO iFreq = 1,kMaxPts
          daaAb(iFreq,iLay) = 0.0
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes the (DOUBLE) array of absorption coeffs
c pretty much the same as above routine
      SUBROUTINE initializeJAC(daaAb)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      DOUBLE PRECISION daaAb(kMaxPtsJac,kProfLayerJac)

      INTEGER iLay,iFreq

      DO iLay = 1,kProfLayerJac
        DO iFreq = 1,kMaxPtsJac
          daaAb(iFreq,iLay) = 0.0
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine converts the abs coeff matrix from 
c double to single daa ---> raa
c and saves it in an overall AbsMatrix raaa
      SUBROUTINE DoSet(daaGasAbCoeff,raaaGasAbCoeff,iCount,iDoAdd)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iCount     = which of the iNumGases are being processed
c daaGasAb   = double precision abs coeffs, from the uncompression
c raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
      DOUBLE PRECISION daaGasAbCoeff(kMaxPtsJac,kProfLayerJac)
      REAL raaaGasAbCoeff(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      INTEGER iCount,iDoAdd

c local variables
      INTEGER iLay,iFr

      IF (iDoAdd .GT. 0) THEN
        DO iLay = 1,kProfLayerJac
          DO iFr = 1,kMaxPtsJac
c          raaaGasAbCoeff(iCount,iFr,iLay) = real(daaGasAbCoeff(iFr,iLay))
            raaaGasAbCoeff(iCount,iFr,iLay) = daaGasAbCoeff(iFr,iLay)
          END DO
        END DO
      ELSEIF (iDoAdd .LE. 0) THEN
        DO iLay = 1,kProfLayerJac
          DO iFr = 1,kMaxPtsJac
            raaaGasAbCoeff(iCount,iFr,iLay) = 0.0
          END DO
        END DO
      END IF

      RETURN
      END
c************************************************************************c
c this subroutine converts the abs coeff matrix from 
c double to single daa ---> raa
      SUBROUTINE DoDtoR(daaGasAbCoeff,raaGasAbCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaGasAb   = double precision abs coeffs, from the uncompression
c raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaGasAbCoeff(kMaxPts,kProfLayer)

c local variables
      INTEGER iLay,iFr

      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaGasAbCoeff(iFr,iLay) = real(daaGasAbCoeff(iFr,iLay))
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine converts the layer iL to 0 
      SUBROUTINE ZeroLayer(raaGasAbCoeff,iL)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaGasAbs = 2d matrix that save ALL abs coeffs for current 25 cm-1 chunk
      REAL raaGasAbCoeff(kMaxPts,kProfLayer)
      INTEGER iL

      INTEGER iFr

      DO iFr = 1,kMaxPts
        raaGasAbCoeff(iFr,iL) = 0.0
      END DO

      RETURN
      END

c************************************************************************
c this subroutine checks parameters in kcarta.param ... abort if they
c do not make sense .. the parameters are set in kcarta.param
      SUBROUTINE CheckKCARTAParameters

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      write(kStdWarn,*) 'checking parameters in kcarta.param'
      write(kStdWarn,*) '  '
 
      CALL Check_kaNum_ka100layerCloudType

c do a quick check of the important parameters set by the user
      IF (kMixFilRows .LT. kProfLayer) THEN
        write(kStdErr,*) 'In kcarta.param, need '
        write(kStdErr,*) 'kMixFilRows >= kProfLayer(=',kProfLayer,')'
        write(kStdErr,*) 'please reset and retry'
        CALL DoSTOP
      END IF

      IF (abs(kXsecFormat) .NE. 1) THEN
        write(kStdErr,*) 'kXsecFormat in kcarta.param must be = +/-1'
        write(kStdErr,*) 'please reset and retry'        
        CALL DoSTOP
      END IF

      write(kStdWarn,*) 'Max #of atmospheres from *RADFIL = ',kMaxAtm
      write(kStdWarn,*) 'Max #of gases from *GAS/XSCFIL = ',kGasStore
      write(kStdWarn,*) 'Max #of mixed paths *MIXFIL = ',kMixFilRows
      write(kStdWarn,*) '  '

      RETURN
      END

c************************************************************************
c set the default parameter values, for those that are not set in *PARAM
c read the parameter file to set parameters that have optional values
      SUBROUTINE SetDefaultParams

c NOTE !!!! also double check subroutine EXTRAPAR in strings2.f
c NOTE !!!! also double check subroutine SETDEFAULTPARAMS in misc.f

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iI,iJ

c acos(3/5) * 180/kPi      
      kThermalAngle = 53.1313010235415598
      
c set default values here
c kLayer2Sp
c     -2     Layer transmittance            t(i)=exp(-k(i))
c     -1     Layer Optical depth            k(i)
c      1     Layer-to-Space Optical Depth   k2s(i)=sum(j=i,n)(k(j))
c      2     Layer-to-Space transmittance   t2s(i)=sum(j=i,n)exp(-k(j))
c      3     Layer-to-Ground Optical Depth  k2s(i)=sum(j=1,i)(k(j))
c      4     Layer-to-Ground transmittance  t2s(i)=sum(j=1,i)exp(-k(j))
      kLayer2Sp    = -1

c kCKD      == -1,00,21,23,24 sets which continuum version calculation to do
c              -1 : no continuum
c ----------------------------------------------------------------------------
c AER-CKD      01 : self, foreign   is by CKD and modified by Mlawer/Tobin
c                                 ... this is the MT_CKD version of Dec 2002
c       AIRS   02 : version 02    ... this is the MT_CKD version of Dec 2002,
c                                     but CS modified by Scott Hannon Aug 2003
c       AIRS   03 : version 03    ... this is the MT_CKD version of Dec 2002,
c                                     but CS modified by Scott Hannon Jan 2004
c       AIRS   04 : version 04    ... this is the MT_CKD version of Dec 2002,
c                                     CS,CF modified by Scott Hannon Jan 2004
c       AIRS   05 : version 05    ... this is the MT_CKD version of Dec 2002,
c                                     CS,CF modified by Scott Hannon Jan 2004
c                                     (so far it looks like CKD v4)
c                                     On top of that it puts Dave Tobin's cs,cf
c                                     between 1300-1800 cm-1
c       AIRS   06 : version 06    ... this is the MT_CKD version of Dec 2002,
c                                     CS,CF modified by Scott Hannon Jan 2004
c                                     extremely similar to CKD 4, except have
c                                     fixed the problem at 850 cm-1 which is
c                                     really HNO3 and not a continuum problem
c In asl:/carrot/s1/strow/Tobin/Tobin_radish/NISTdata2/New/, use
c cs0_tsl0_hr.mat and cf0_tsl0_hr.mat.  makecons.m with flag=7
c looks like it should load things in correctly.
c ----------------------------------------------------------------------------
c AER-STD      25 : self, foreign   is by CKD and modified by Mlawer/Tobin
c                                 ... this is the MT_CKD 2.5 version of Sep 2011
c ----------------------------------------------------------------------------
c AER-STD      27 : self, foreign   is by CKD and modified by Mlawer/Tobin
c                                 ... this is the MT_CKD 2.5 version of Feb 2016
c ----------------------------------------------------------------------------
c old versions of AER CKD
c AER-CKD      00 : version 00 
c AER-CKD      21 : version 2.1
c AER-CKD      23 : version 2.3
c AER-CKD      24 : version 2.4
c ----------------------------------------------------------------------------
c ----------------------------------------------------------------------------
c ----------------------------------------------------------------------------
c these were our interim versions (made by us, basically MT_CKD 1)
c no longer in use!
c
c    RAL       12 : self = version 2.4, foreign = blend of 0.0 and 2.4
c                                       from 0-1575, 1625-3000, use v2.4
c                                       from 1575-1625 linearly blend v2.4,0.0
c                                       using a triangle
c    RAL       13 : self = version 2.3, foreign = dave tobin's v2.3
c    RAL       90 : self = version 2.4, foreign = blend of 2.4,Tobin's thesis
c                                       from 0-1300, 1900-3000, use v2.4
c                                       from 1300-1900 use Tobins thesis
c    RAL       50 : self, foreign       blend of 2.4, RAL data
c                                       from 0-1300, 1900-3000, use v2.4
c                                       from 1300-1900 use RAL data
c                                       mst50.m used linear tempr interpolation
c    RAL       51 : self, foreign       blend of 2.4, RAL data
c                                       from 0-1300, 1900-3000, use v2.4
c                                       from 1300-1900 use RAL data
c                                       mst51.m used power tempr interpolation
c                                       Need to be careful ... I put in Dave 
c                                       Tobin's fixes for CS in the 600-1100
c                                       cm-1 range; see mst51_tobin.m in
c                                       the SPECTRA directory (~0.9*ckd24) 
c    RAL       52 : self, foreign       blend of 2.4, RAL data
c                                       from 0-1300, 1900-3000, use v2.4
c                                       from 1300-1900 use RAL data
c                                       mst52.m used power tempr interpolation
c                                       and says CF(296)=CF(243); so use the
c                                       foreign broadened 243 data to get CS243
c    RAL       55 : self, foreign       same as 51 above, but uses
c                                a) CS fudge factor of 0.87 for v <= 1137 cm-1
c                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
c                                c) some fudge factors for CF in 1400-1700 cm-1
c                                d) CKD 2.4 used upto 1400 cm-1
c    RAL       56 : self, foreign       same as 51 above, but uses
c                                a) John Taylor fudge between 600-1200 cm-1
c                                   instead of Dave Tobin fudge
c                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
c                                c) some fudge factors for CF in 1400-1700 cm-1
c                                   these are better than the above fudges
c                                d) CKD 2.4 used upto 1400 cm-1
c    RAL       60 : self, foreign     is a hybrid of 51,55 done by scott
c                                     so that bias errors are reduced, as
c                                     a function of water column amount. 
c                                     eventually will include tuning coeffs
c ----------------------------------------------------------------------------
c now we only allow 
c %% this is new list
c %% 0,21,23,24 are the 1990s CKD
c %% 1,25       are new MT-CKD
c %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
c %%            are derived from MT-CKD25
c
c origCKD = [0 21 23 24];
c MTCKD1  = [ [1] [4 6]];
c MTCKD25 = [ [25]  [] ];
c allowedCKD = [origCKD MTCKD1 MTCKD25];
c ----------------------------------------------------------------------------

      kCKD = 25

c kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
c              -1 if we just do the weighted average to find the MixVertTemps
      kGasTemp = -1

c kLongOrShort == whether the user wants to save ALL header info (+1) or
c                 just a portion (-1) or output just the basic results (0)
      kLongOrShort = 1

c kActualJacs = -1 if we compute and output ALL profile(z) jacs
c               20 if we compute only Q(z) jacs, output 0's everywhere else
c               30 if we compute only T(z) jacs, output 0's everywhere else
c               40 if we compute only W(z) jacs, output 0's everywhere else
c               50 if we compute only S(z) jacs, output 0's everywhere else
c              100 if we compute only stemp and column gas jacs
c for the following the d/dT only uses gases in iaJacob{}
c kActualJacs = -2 if we compute and output ALL profile(z) jacs
c               32 if we compute only T(z) jacs, output 0's everywhere else
c              102 if we compute only stemp and column gas jacs
      kActualJacs = -1    !!! default

c kActualJacsT = -1, kActualJacsB = -1
c if we set kActualJacs = 100, then we can also set the layers of the column
c that we want to perturb; default = -1/-1 means all layers (kActualJacs = 100)
c else kActualJacsABCXYZ sets these as kActualJacsT=ABC, kActualJacsB=XYZ
      kActualJacsT = -1
      kActualJacsB = -1

c kJacobOutput == -1 if we output d(radiance)/dq,d(radiance)/dT
c                  0 if we output d(radiance)/dq * q, d(radiance)/dT
c                  1 if we output d(BT)/dq * q, d(BT)/dT
      kJacobOutput = 1

c kFlux == -1 if we do not want flux/OLR computations or output NLTE 
c                                                        Planck modifiers
c       ==  1 if we want DNWELL flux at every layer      units = mW m-2 
c       ==  2 if we want flux computations : output      units = kelvin day-1
c       ==  3 if we want UPWELL flux at every layer      units = mW m-2 
c       ==  4 if we want outgoing OLR only at TOA,       units = mW m-2
c       ==  5 if we want outgoing OLR at TOA, ILR at GND units = mW m-2
c       ==  6 if we want DNWELL and UPWELL flux at each layer units = mW m-2
c   --> ==  0 if we want to output NLTE Planck modifiers
      kFlux = -1

c kPlanckOut == -1 if we do not want to output Planck modifiers, 1 if we do
      kPlanckOut = -1

c kSurfTemp = -1 == want to use user supplied surface temp in *RADNCE
c              1 == want to use user supplied surface temp in *RADNCE as 
c                     an offset to pressure interpolated temperature
      kSurfTemp = -1     

c kRTP = -5  : read LBLRTM style LAYERS profile (edited TAPE 6); set atm from namelist
c kRTP = -6  : read LBLRTM style LEVELS profile (       TAPE 5); set atm from namelist
c kRTP = -10 : read TEXT style LEVELS   profile; set atm from namelist
c kRTP = -2  : read GENLN4 style LAYERS profile; set atm from namelist
c kRTP = -1  : read old style kLAYERS   profile; set atm from namelist
c kRTP =  0  : read RTP style kLAYERS   profile; set atm from namelist
c kRTP = +1  : read RTP style kLAYERS   profile; set atm from RTP file
c kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
      kRTP = 1

c kTempJac == -2 if we only want d/dT(planck) in temperature jacobian
c             -1 if we only want d/dT((1-tau)(tau->sp)) in temp jacobian
c              0 if we want complete d/dT(planck) + d/dT(tau) in temp jac
      kTempJac = 0

c the following cannot be controlled by the user using *PARAMS
c all the radiance parameters and the Jacobian parameter
c kJacobian ==  1 if analytic Jacobians are to be computed
c              -1 if we do the standard Genln2 computations w/o Jacobians
      kSolar        = 1     !turn on solar
      kSolarAngle   = 0.0   !solar angle
      kSolarRefl    = -1.0  !use (1-ems)/pi
      kThermal      = 0     !use fast diffusive approx
      kThermalAngle = -1.0  !use acos(3/5) in upper layers
      kThermalJacob = 1     !use thermal backgnd in Jacobians

      kJacobian = -1        !do not do Jacobians
      kScatter  = -1        !do not do scattering computations

      k100layerCloud = -1   !assume rtp file does NOT have 100 layer cloud

c 2016
c allow nm_params to define defaults
c   GENERAL iaDefaults(1,:) = iSARTAChi   iSPlineType iCO2Chi  iMatlabORf77   iWhichScatterCode iMethod
c   RT      iaDefaults(2,:) = iGaussPts iGaussQuad  iSnell      iInterpType  iWhichRT  (kTemperVary set in mn_radnce)
c   NLTE    iaDefaults(3,:) = iCurrent    iTalk       iTestGenln2  iNoPressureShiftCO2 iQtips_H98
c                             iLinearOrSpline iDoCO2Continuum iMethod
c   TAPE5/6     iaDefaults(5,:) = iReplaceZeroProf iAIRS101_or_LBL_levels IPLEV iAddLBLRTM 
c      INTEGER iaaOverrideDefault(8,10)
c      COMMON/comBlockDefault/iaaOverrideDefault
      DO iI = 1,4
        DO iJ = 1,10
	  iaaOverrideDefault(iI,iJ) = -9999
	END DO
      END DO
      
c GENERAL
      caaTextOverrideDefault  = 'notset'
      
      iaaOverrideDefault(1,1) = -1    !!! iSARTAChi = -1  for no tuning, see kcartabasic/kcartamain/kcartajpl
                                      !!!                 kcartaparallel and finally used in kcoeffMAIN.f
      iaaOverrideDefault(1,2) = +1    !!! iSplinetype = +1 for SUBR iSetSplineType in kcartamisc.f
      iaaOverrideDefault(1,3) = +2    !!! iCO2Chi = +2     for SUBR multiply_co2_chi_functions in kcoeffMAIN.f
      iaaOverrideDefault(1,4) = +1    !!! iMatlabORf77 = +1  use Maltab style uncompression,  kcoeffMAIN.f
      iaaOverrideDefault(1,5) = +5    !!! iWHichScatterCode = 5 for PCLSAM in rtp_interface.f
      iaaOverrideDefault(1,6) = +1    !!! iReadP = 1 when assuming GENLN2 style profile in n_pth_mix.f
      iaaOverrideDefault(1,7) = -1    !!! iLogOrLinear = -1 when interp scat tables SUBR INTERP_SCAT_TABLE2
                                      !!!   in clear_scatter_misc.f
      iaaOverrideDefault(1,8) = -1    !!! -1 to keep h.vcmin/h.vcmax as read in from RTPfile, +1 to override with rf1,rf2
      
c RadTrans
      iaaOverrideDefault(2,1) = kTemperVary !!! kTemperVary .... can be reset in nm_radnce, and then subr SetkTemperVary
      iaaOverrideDefault(2,2) = +3    !!! THIS IS LBLRTM STYLE iGaussPts = 3 for flux and downwell gauss quad
                                      !!!   see SUBR IntegrateOverAngles_LinearInTau in rad_quad.f
      iaaOverrideDefault(2,3) = 0     !!! SUBR BackGndThermal in rad_diff.f
                                      !!! iDothermal = kThermal; if iDoThermal = -1, no backgnd thermal computed
      				      !!!                                      =  0, backgndthermal with diffusive approx << DEFAULT >>
				      !!!                                            --->>> control further with iaaOverrideDefault(2,4) <<<---
				      !!!                                      = +1, use integration over angles, const-in-tau  layer T
				      !!!                                            --->>> control further with iaaOverrideDefault(2,5) <<<---
				      !!!                                      = +2, use integration over angles, linear-in-tau layer T
				      !!!   this is the main routine, called by all downwelling RT routines in rad_main.
				      !!!   all of them have -1 for iDoAcos35
				      !!!     calls IntegrateOverAngles_LinearInTau (iDoThermal = 2) 
				      !!!     calls IntegrateOverAngles             (iDoThermal = 1) -- also can set iaaOverrideDefault(2,5)
				      !!!     calls DoDiffusivityApprox             (iDoThermal = 0) << DEFAULT >>
      iaaOverrideDefault(2,4) = -1    !!! SUBR radnce4RTP in rtp_interface.f
                                      !!!   raKThermalAngle(iC) = iaaOverrideDefault(2,4) in rtp_interface.f
                                      !!!     = -1, fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
				      !!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 (acos(3/5))
				      !!!
				      !!!   this sets  kSetThermalAngle = -1 for acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
				      !!!                               = +1 for constant angle (typically acos(3/5)) in all layers
				      !!!                               = -2 for same as -1, except linear-in-tau T variation				      
				      !!! SUBR DoDiffusivityApprox in rad_diff.f uses this info
				      !!!   iDiffMethod = kSetThermalAngle
				      !!!     = -1 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
				      !!!          differs from iaaOverrideDefault(2,5) = 0 since here, upper layers use acos(3/5) lower layers are accurate
				      !!!                                                   while there, use layer-varying accurate acos(rDiffusive)
				      !!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 = acos(3/5) << DEFAULT >>
				      !!!           this can be controlled by kThermalAngle, either in nm_params for kRTP = +1
				      !!!                                                   or rakThermalAngle() for kRTP = 0,-1
				      !!!           so in nm_params : set iaaOverride(2,4) = 1, kThermalAngle = 50.0 and thay works!!!
				      !!!     = -2 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers, linear in tau T
				      !!!     = +2 diffusive background using LBLRTM style 3 exponetial gauss quad, not yet implemented
      iaaOverrideDefault(2,5) = 0     !!! SUBR IntegrateOverAngles in rad_quad.f, called by SUBR BackGndThermal
                                      !!!   iGaussQuad =    -1 for integrate using newton quad 0:90/20:90 (VERY SLOW)
                                      !!!                    0 for accurate diffusivity                   (AT ALL LAYERS << DEFAULT >>)
				      !!!                      so this differs from iaaOverrideDefault(2,3) = 0,iaaOverrideDefault(2,4) = -1
				      !!!                      where acos(3/5) is used in upper ayers, and accurate diffusive angle in lower layers
				      !!!                   +1 for gausslegendre w(i) at theta(i)         (QUITE SLOW)
      iaaOverrideDefault(2,6) = +1    !!! iUsualUpwell = +1 for upwell RT with surface term, << DEFAULT >>
                                      !!!                -1 with no surface,
                                      !!!                -2 to only dump downwelling backgnd
				      !!!   see SUBR find_radiances in rad_main.f
      iaaOverrideDefault(2,7) = -1    !!! iUseSnell = -1 for No Snell law raytrace, similar to SARTA
                                      !!!   see SUBR FindLayerAngles in rad_angles.f				     
      iaaOverrideDefault(2,8) = +1    !!! iInterpType = +1 to turn (pav,Tav) into (plevs,Tlevs), only used if kTemperVary = 43
                                      !!!   see SUBR Get_Temp_Plevs in n_pth_mix.f
      iaaOverrideDefault(2,9) = -1    !!! iLBLRTM_highres = -1 do not estimate/fix problems because use 0.0025 cm-1, when kTemperVary = 43 << DEFAULT>>
                                      !!!   see SUBR rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE_EMISS in rad_main.f				     
  
c n_layers_lblrtm.f and n_pth_mix.f  TAPE5/6
      iaaOverrideDefault(3,1) = -1    !!! iAIRS101_or_LBL_levels use LBLRTM, not AIRS 101 levels, for integration
      iaaOverrideDefault(3,2) = +1    !!! iReplaceZeroProf = +1 to add in profiles TAPE5 does not have
      iaaOverrideDefault(3,3) = -1    !!! iAddLBLRTM = -1 when gas profile missing from TAPE5/6, do not add it in
                                      !!!   in n_pth_mix.f

      RETURN
      END 

c************************************************************************
c now check parameters in *PARAM
c this subroutine checks parameters in *PARAMS ... abort if they
c do not make sense .. 
      SUBROUTINE CheckParams

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

      INTEGER i0,iT,iB,iJ,iGah,iConstOrVary
      CHARACTER*9 iIOUN9

      write(kStdWarn,*) 'checking parameters (from *PARAMS) .... '
      write(kStdWarn,*) '  '

      IF ((iabs(kLayer2Sp) .NE. 1) .AND. (iabs(kLayer2Sp) .NE. 2) .AND. (kLayer2Sp .NE. 3) .AND. (kLayer2Sp .NE. 4)) THEN
        write(kStdErr,*) 'In *PARAMS, need kLayer2Sp = +/-1,+/-2,+3,+4'
        write(kStdErr,*) 'kLayer2Sp == do layer-to-space calc or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF
      IF (iabs(kGasTemp) .NE. 1) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kGasTemp = +/- 1'
        write(kStdErr,*) 'kGasTemp = use CO2 temperature profile or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

c CKD    releases are 0,21,23,24
c MT_CKD releases are 1,25
c our modifications are 2,3,4,5 (derived from 1)
c                       51,55,60 (derived from analysing RAL data)
c and various 12,13,50,52,56 which might be gotten rid of eventually
c ----------------------------------------------------------------------------
c now we only allow 
c %% this is new list
c %% 0,21,23,24 are the 1990s CKD
c %% 1,25       are new MT-CKD
c %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
c %%            are derived from MT-CKD25
c
c origCKD = [0 21 23 24];
c MTCKD1  = [ [1] [4 6]];
c MTCKD25 = [ [25] [] ];
c allowedCKD = [origCKD MTCKD1 MTCKD25];
c ----------------------------------------------------------------------------

      IF ((kCKD .NE. -1) 
     !!!! the std CKD pre-2002 versions
     $    .AND. (kCKD .NE. 0) .AND. (kCKD .NE. 21) .AND. (kCKD .NE. 23) .AND. (kCKD .NE. 24) 
     !!!! these are MT_CKD1 and research versions from AIRS data
     $    .AND. (kCKD .NE. 1) .AND. (kCKD .NE. 4) .AND. (kCKD .NE. 6)
     $    .AND. (kCKD .NE. 25) .AND. (kCKD .NE. 27))
     $ THEN 
        write(kStdErr,*) 'In *PARAMS, need kCKD = [-1] for no continuum OR'
        write(kStdErr,*) '                 CKD    versions 0,21,23 or 24'
        write(kStdErr,*) '              MT_CKD    versions 1,  [4,6]'
        write(kStdErr,*) '              MT_CKD    versions 25  [   ]'
        write(kStdErr,*) '              MT_CKD    versions 27  [   ]'	
        write(kStdErr,*) '       (latest AER versions =  1, released Dec 2002)'
        write(kStdErr,*) '       (latest AER versions = 25, released Dec 2010)'
        write(kStdErr,*) '       (latest AER versions = 27, released Feb 2016)'	
        write(kStdErr,*) '           [ are our modifications ] '
        write(kStdErr,*) 'kCKD is water continuum calculation version'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

      IF (iabs(kLongOrShort) .GT. 2) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kLongOrShort = -2,-1,0,+1,+2'
        write(kStdErr,*) 'kLongOrShort = complete header info (+1) or not (-1),  long warning.msg'
        write(kStdErr,*) 'kLongOrShort = complete header info (+2) or not (-2), short warning.msg'
        write(kStdErr,*) '               or file containing results only (0)'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

c kActualJacs = -1 if we compute and output ALL profile(z) jacs
c               20 if we compute only Q(z) jacs, output 0's everywhere else
c               30 if we compute only T(z) jacs, output 0's everywhere else
c               40 if we compute only W(z) jacs, output 0's everywhere else
c               50 if we compute only S(z) jacs, output 0's everywhere else
c              100 if we compute only stemp and column gas jacs
c for the following the d/dT only uses gases in iaJacob{}
c kActualJacs = -2 if we compute and output ALL profile(z) jacs
c               32 if we compute only T(z) jacs, output 0's everywhere else
c              102 if we compute only stemp and column gas jacs

      IF ((kActualJacs .NE. -1)  .AND. (kActualJacs .NE. 20)   .AND.
     $    (kActualJacs .NE. 30)  .AND. (kActualJacs .NE. 40)   .AND.
     $    (kActualJacs .NE. 50)  .AND. (kActualJacs .NE. 100)  .AND.
     $    (kActualJacs .NE. -2)  .AND. (kActualJacs .NE. 32)   .AND.          
     $    (kActualJacs .NE. 102) .AND. (kActualJacs .LT. 102)) THEN
        write(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
        write(kStdErr,*) '  or 100 or 100ABCXYZ'
        write(kStdErr,*) 'OR -2, 32,102 or 102ABCXYZ'
        write(kStdErr,*) 'kActualJacs = 100(2) ==> column gas/temp jacs '
        write(kStdErr,*) '   for all layers, plus stemp jac'
        write(kStdErr,*) 'kActualJacs = 100(2)ABCXYZ ==> column gas/temp jacs '
        write(kStdErr,*) '   for layers ABC to XYZ, plus stemp jac'
        write(kStdErr,*) 'kActualJacs = actually compute all profile jacs (-1)'
        write(kStdErr,*) '              or Q(z),T(z),W(z) jacs (0s elsewhere)'
        write(kStdErr,*) ' '
        write(kStdErr,*) 'You have set this as ',kActualJacs
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF      
      iT = -1
      iB = -1
      iGah = -1

      IF (kActualJacs .GT. 102) THEN
        IF (kActualJacs .LT. 100000001) THEN
          write(kStdErr,*) 'too few characters in kActualJacs : ',kActualJacs
          write(kStdErr,*) 'check it !'
          CALL DoStop
      END IF
        IF (kActualJacs .GT. 102999999) THEN
          write(kStdErr,*) 'too many characters in kActualJacs : ',kActualJacs
          write(kStdErr,*) 'max possible value is 102999999'
          CALL DoStop
      END IF

        i0 = kActualJacs
        !! i0 better be 9 characters long
        write(iIOUN9,99) kActualJacs
        iJ = 9
 10     CONTINUE
        IF ((iIOUN9(iJ:iJ) .NE. ' ') .AND. (iJ .GE. 1)) THEN
c         write(kStdErr,*) 'good',iJ,iIOUN9(iJ:iJ),kActualJacs
          iJ = iJ - 1
          GOTO 10
        ELSEIF (iJ .GT. 0) THEN
          iGah = iJ
          write(kStdErr,*) 'space in  kActualJacs at ',iJ,iIOUN9(iJ:iJ)
        END IF
        IF (iGah .GT. 0) THEN
          write(kStdErr,*) 9-iGah,' chars in kActualJacs = ',kActualJacs
          write(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
          write(kStdErr,*) '  or 100 or 100ABCXYZ'
          write(kStdErr,*) '  or 102 or 102ABCXYZ'
          write(kStdErr,*) 'need 9 characters 10E ABC XYZ .. try again'
          CALL DoStop
        END IF

        write(kStdWarn,*) 'kActualJacs passed test ... '
        IF (kActualJacs .LE. 102000000) THEN
          iT = i0 - 100000000
          iT = int(iT/1000)
          iB = i0 - int(i0/1000)*1000
          kActualJacs = 100
        ELSE
          i0 = i0 - 2000000
          iT = i0 - 100000000
          iT = int(iT/1000)
          iB = i0 - int(i0/1000)*1000
          kActualJacs = 102
        END IF

        IF (iT .LT. iB) THEN
          iJ = iT
          iT = iB
          iB = iJ
        END IF
        IF (iT .GT. kProfLayer) THEN
          write(kStdWarn,*) 'IT = ',iT,' greater than kProfLayer = ',kProfLayer
          write(kStdWarn,*) 'resetting iT = kProfLayer'
          iT = kProfLayer
        END IF
        IF (iB .GT. kProfLayer) THEN
          write(kStdWarn,*) 'IB = ',iB,' greater than kProfLayer = ',kProfLayer
          write(kStdWarn,*) 'resetting iB = 1'
          iB = 1 
        END IF
        kActualJacsT = iT
        kActualJacsB = iB
      END IF        
 99   FORMAT(I9)

      IF ((iabs(kJacobOutput) .NE. 1) .AND. (kJacobOutput .NE. 0) .AND.
     $     (kJacobOutput .NE. 2))  THEN
        write(kStdErr,*) 'kJacobOutput = ',kJacobOutput
        write(kStdErr,*) 'In *PARAMS, need kJacobOutput =-1,0,+1,+2'
        write(kStdErr,*) 'kJacobOutput = format to output Jacobians'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

      IF ((kFlux .LT. -1) .OR. (kFlux .GT.  6)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kFlux =-1 OR 1,2,3,4,5,6'
        write(kStdErr,*) 'where kFlux = do/do not compute fluxes'
        write(kStdErr,*) 'OR         program needs kFlux =-1,+1,2,3,4,5,6 OR 0'
        write(kStdErr,*) 'where kFlux = do not/do  output NLTE Planck'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

ccccc this has been replaced by kSurfTemp < 0 : use raTSurf(iI)
ccccc                           kSurfTemp > 0 : use raTSurf(iI) + InterpedTemp
c      IF ((abs(kSurfTemp)-1.0) .GE. 1e-5) THEN
c        write(kStdErr,*) 'In *PARAMS, program needs kSurfTemp = +/-1.0'
c        write(kStdErr,*) 'where kSurfTemp tells the program how to use'
c        write(kStdErr,*) 'the surftemperatures in *RADNCE'
c        write(kStdErr,*) 'Please reset and retry'
c        CALL DoSTOP 
c      END IF

      !!!kRTP = -6 : read LBLRTM       LAYERS profile; set atm from namelist
      !!!kRTP = -5 : read LBLRTM       LEVELS profile; set atm from namelist
      !!!kRTP = -10 : read LEVELS            profile; set atm from namelist
      !!!kRTP = -2  : read GENLN2 style LAYERS profile; set atm from namelist
      !!!kRTP = -1  : read old    style LAYERS profile; set atm from namelist
      !!!kRTP =  0  : read RTP   style kLAYERS profile; set atm from namelist
      !!!kRTP = +1  : read RTP   style kLAYERS profile; set atm from RTP file
      !!!kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
      IF ((kRTP .NE. -5) .AND. (kRTP .NE. -6) .AND. (kRTP .NE. +2) .AND. 
     $    (kRTP .NE. -10) .AND. ((kRTP .LT. -2) .OR. (kRTP .GT. 1))) THEN
        write(kStdErr,*) 'Need to set RTP = -10,-6,-5,-2,-1,0,+1,+2'
        write(kStdErr,*) 'Please reset kRTP and retry'
        CALL DoSTOP 
      END IF

      IF ((kSurfTemp .GT. 0) .AND. (kRTP .EQ. 1)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from RTP file'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP 
      END IF

      IF ((kSurfTemp .GT. 0) .AND. (kRTP .EQ. 2)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from JPL/NOAA input'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP 
      END IF

      IF ((kSurfTemp .GT. 0) .AND. ((kRTP .EQ. -5) .OR. (kRTP .EQ. -6))) THEN
        write(kStdErr,*) 'Will read surface temperature info from LBLRTM file'
        write(kStdErr,*) 'and ask kCARTA to add on raTSurf offset from nm_radnces!!!'
c        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
c        CALL DoSTOP 
      END IF

      IF ((kSurfTemp .GT. 0) .AND. (kRTP .EQ. -10)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from LEVELS TXT file'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP 
      END IF
 
      IF ((kTempJac .LT. -2) .OR. (kTempJac .GT. 0)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kTempJac=-2,-1,0'
        write(kStdErr,*) 'where kTempJac = use Planck or tau or both '
        write(kStdErr,*) 'when doing d/dT'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
      END IF

      RETURN
      END


c************************************************************************
c this subroutine sets the mixed path effective tempertures, weighted 
c according to the mixing table
      SUBROUTINE GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp,
     $                          iNpmix,iCO2)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNumGases   = number of gases read in from *GASFIL + *XSCFIL
c iNpMix      = number of mixed paths read in from *MIXFIL
c raaTemp     = temperature profile of ONE of the gases (assume all equal)
c raMixVTemp  = computed vertical temp profile for mixing table
c raaMix      = mixing table from *MIXFIL
c iCO2        = which of the gases is CO2 .. if -1, none are CO2
      INTEGER iNumGases,iNpmix,iCO2
      REAL raMixVertTemp(kMixFilRows)
      REAL raaTemp(kProfLayer,kGasStore)
      REAL raaMix(kMixFilRows,kGasStore)

c local variables
      INTEGER iI,iJ,iL,MP2Lay,iBad
      REAL rT,rW

      iBad = 0
      DO iI = 1,kMixFilRows
        raMixVertTemp(iI) = 0.0
      END DO

c kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
c              -1 if we just do the weighted average to find the MixVertTemps
c default = -1

      rW = 0.0
      DO iJ = 1,iNumGases
        rW = rW+raaMix(50,iJ)
      END DO
      IF (rW .LE. 0.00001) THEN
        write(kStdErr,*) 'arbitrarily tested LAYER 50, mixed path weights = 0 .. perhaps you have nm_weight wrong???'
        write(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
        CALL DoStop
      END IF
        
      IF ((kGasTemp .EQ. 1) .AND. (iCO2 .GT. 0)) THEN
c user wants the CO2 profile to be the temperature profile
        DO iI = 1,iNpmix
          iL = MP2Lay(iI)
          raMixVertTemp(iI) = raaTemp(iL,iCO2)
        END DO
c      ELSEIF ((kGasTemp .EQ. -1) .AND. (iCO2 .GT. 0) .AND. raaMix(50,iCO2) .LE. 0.001) THEN
c user wants the CO2 profile to be the temperature profile, but CO2 wgt == 0, so this is an OOPS moment, quick fix
c        write(kStdErr,*) 'oops in GetMixVertTemp,,kGasTemp,iCO2 = ',kGasTemp,iCO2
c        DO iI = 1,iNpmix
c          iL = MP2Lay(iI)
c          raMixVertTemp(iI) = raaTemp(iL,iCO2)
c          write(kStdErr,*) 'oops in GetMixVertTemp, iI,raMixVertTemp(iI) = ',iI,raMixVertTemp(iI)
c        END DO
      ELSE
c calculate the weights      
        DO iI = 1,iNpmix
          rT = 0.0
          rW = 0.0
          iL = MP2Lay(iI)
          DO iJ = 1,iNumGases
            rT = rT+raaTemp(iL,iJ)*raaMix(iI,iJ)
            rW = rW+raaMix(iI,iJ)
          END DO
          rT = rT/rW
          IF (rW .LE. 0.00001) THEN
            iBad = iBad + 1
            write(kStdErr,*) 'hmm, mixed path weight = 0 in GetMixVertTemp for layer ',iI
          END IF
          raMixVertTemp(iI) = rT
c if the weights are set so that mixed path defines unique layers
c these temperatures should now be equal
        END DO
        IF (iBad .GT. 0) THEN
          write(kStdErr,*) 'had total ',iBad,' mixed path weights = 0 .. perhaps you have nm_weight wrong???'
          write(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
          CALL DoStop
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this subroutine computes avg pressure of the layers
      SUBROUTINE FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers

      INTEGER iI,iJ,iDiff

      DO iI = 1,kProfLayer
        pProf(iI) = 0.0
      END DO

      iDiff = kProfLayer - iProfileLayers
      DO iI = 1,iProfileLayers
        iJ = iI + iDiff
        pProf(iJ) = raPressLevels(iJ+1)-raPressLevels(iJ)
        pProf(iJ) = pProf(iJ)/log(raPressLevels(iJ+1)/raPressLevels(iJ))
      END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
      SUBROUTINE databasestuff(iLowerOrUpper,
     $                 raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)


      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output params
      REAL raDatabaseHeight(kMaxLayer)
      REAL raDATABASELEVHEIGHTS(kMaxLayer+1)
      REAL raDATABASELEV(kMaxLayer+1)
c input params
      INTEGER iLowerOrUpper

      IF (iLowerOrUpper .LT. 0) THEN
        CALL databasestuff_lower(iLowerOrUpper,
     $                 raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
      ELSEIF (iLowerOrUpper .GT. 0) THEN
        CALL databasestuff_upper(iLowerOrUpper,
     $                 raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
c for the lower atm
      SUBROUTINE databasestuff_lower(iLowerOrUpper,
     $                 raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)


      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/airsheights.param'
      include '../INCLUDE/airslevels.param'
      include '../INCLUDE/airslevelheights.param'

c output params
      REAL raDatabaseHeight(kMaxLayer)
      REAL raDATABASELEVHEIGHTS(kMaxLayer+1)
      REAL raDATABASELEV(kMaxLayer+1)
c input params
      INTEGER iLowerOrUpper

c local vars 
      INTEGER iI      

      IF (iLowerOrUpper .GT. -1) THEN
        write(kStdErr,*) 'trying to make default lower atm profile'
        CALL DoStop
      END IF

      DO iI = 1,kMaxLayer
        raDatabaseHeight(iI) = DatabaseHeight(iI)
      END DO
      DO iI = 1,kMaxLayer + 1
        raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
      END DO
      DO iI = 1,kMaxLayer
        raDATABASELEV(iI) = DATABASELEV(iI)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
c for the uuper atm
      SUBROUTINE databasestuff_upper(iLowerOrUpper,
     $                 raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/airsheights_upper.param'
      include '../INCLUDE/airslevels_upper.param'
      include '../INCLUDE/airslevelheights_upper.param'

c output params
      REAL raDatabaseHeight(kMaxLayer)
      REAL raDATABASELEVHEIGHTS(kMaxLayer+1)
      REAL raDATABASELEV(kMaxLayer+1)
c input params
      INTEGER iLowerOrUpper

c local vars 
      INTEGER iI      

      IF (iLowerOrUpper .LT. +1) THEN
        write(kStdErr,*) 'trying to make default upper atm profile'
        CALL DoStop
      END IF

      DO iI = 1,kMaxLayer
        raDatabaseHeight(iI) = DatabaseHeight(iI)
      END DO
      DO iI = 1,kMaxLayer + 1
        raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
      END DO
      DO iI = 1,kMaxLayer
        raDATABASELEV(iI) = DATABASELEV(iI)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine will take in 100 AIRS layering stuff and interpolate to
c the new arbitrary layering
c WARNING : this assumes that the user has not mucked up KLAYERS layering
c           such that highest Z pressure (lowest pressure) is NOT TOA
c           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
c do the lower atm (usual -1) or upper atm (NLTE +1)

c kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
c so at arbitrary pressure layering, it deals with abs coeffs
c so we do not need raRamt
c but we do need the interpolated temp and partial pressures

c see subr AddOnAFGLProfile_arblevels in n_pth_mix.f
      SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $           raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $           raaPress,iGas,iGasID,iNumLayers,
     $           raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

      IMPLICIT NONE

      INTEGER iPLEV
      
      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/KCARTA_database.param'
      include '../INCLUDE/airslevelheights.param'
      
c  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
c  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
c The above equation and 3 data points define the 101 AIRS levels, which
c are in airslevels.param

c input
c do the lower atm (usual -1) or upper atm (NLTE +1)
      INTEGER iLowerOrUpper
c these are the individual reference profiles, at kMaxLayer layers
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
c these are the arbitrary profiles stored in matrices
      REAL raaPress(kProfLayer,kGasStore)
      INTEGER iError,iGas,iGasID,iNumLayers,iSplineType
c these are the kLAYERS pressure levels, layer thick for the current profile
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
c  output
c these are the individual reference profiles, at kProfLayer layers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      REAL pMax100,pMin100

c local variables
      INTEGER iI,iJ,iL,iG,iZbndFinal,iNot,iX,iY
      REAL raWorkP(kMaxLayer),raXgivenP(kMaxLayer),
     $     raYgivenP(kMaxLayer),raY2P(kMaxLayer)
      REAL raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2,rPPWgt
      REAL raSortPressLevels(kMaxLayer+1)
      REAL raSortPressHeights(kMaxLayer+1)
      REAL raPPX2(kProfLayer),raQX2(kProfLayer)

      REAL raDataBaseThickness(kMaxLayer)
c      REAL DatabaseHeight(kMaxLayer)
c      REAL DATABASELEVHEIGHTS(kMaxLayer+1)
c      REAL DATABASELEV(kMaxLayer+1)

      INTEGER iaBnd(kProfLayer+1,2)
      REAL    raBndFrac(kProfLayer+1,2)
      REAL rPP,rWgt,rMR,rFrac,rMolecules,rHeight,rQtot

c pressure variables!!!!! ----------------->
c raaPress in atm

      !!!this tells how many layers are NOT dumped out by kLAYERS, iNot is reset below
      iNot = kProfLayer-iNumLayers

c simply put in the pressures
      DO iI = 1,iNot
        !these are "junk"
        raRPress(iI) = raaPress(iNot+1,iGas)
      END DO
      DO iI = iNot+1,kProfLayer
        raRPress(iI) = raaPress(iI,iGas)
      END DO

c now just happily spline everything on!!!!!! for the temps
C     Assign values for interpolation
C     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
      rYP1 = 1.0E+16
      rYPN = 1.0E+16
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
      END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        raRTemp(iI) = +999.999
      END DO
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        IF (iSplineType .EQ. +1) THEN
          CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)	
        ELSE
          CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
        END IF
        raRTemp(iI) = r
      END DO
      
      DO iL = 1,kProfLayer
        raRAmt(iL) = 0.0
	raRPartPress(iL) = 0.0
      END DO
      
      !!!this tells how many layers are NOT dumped out by kLAYERS
      iZbndFinal = kProfLayer-iNumLayers

 345  FORMAT(I3,2(' ',F10.3),2(' ',I3,F10.3,' ',F10.3))
      !! look at the LAYERS and figure out which PLEV_KCARTADATABASE_AIRS bracket them
      iNot = (kProfLayer) - (iNumLayers)+1
      DO iL = iNot,kProfLayer
        !! find plev_airs which is just ABOVE the top of current layer
        iG = kProfLayer+1
 10     CONTINUE
	IF ((PLEV_KCARTADATABASE_AIRS(iG) .LE. raPressLevels(iL+1)) .AND. (iG .GT. 1)) THEN
	  iG = iG - 1
	  GOTO 10
	ELSE
	  iaBnd(iL,2) = min(iG+1,kMaxLayer+1)   !! top bndry of plevs_database is lower pressure than top bndry of raPressLevels layer iL
	END IF
	raBndFrac(iL,2) = (raPressLevels(iL+1)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))/
     $                  (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))

        !! find plev_airs which is just BELOW the bottom of current layer
        iG = 1
 20     CONTINUE	
	IF (PLEV_KCARTADATABASE_AIRS(iG) .GT. raPressLevels(iL)) THEN
	  iG = iG + 1
	  GOTO 20
	ELSE
	  iaBnd(iL,1) = max(iG-1,1) !! bot boundary of plevs_database is bigger pressure than top bndry of raPressLevels layer iL
	END IF
	raBndFrac(iL,1) = (raPressLevels(iL)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))/
     $                  (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))	
	
c      write (*,345) iL,raPressLevels(iL),raPressLevels(iL+1),iaBnd(iL,1),raBndFrac(iL,1),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)),
c     $                                                       iaBnd(iL,2),raBndFrac(iL,2),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))
      END DO
c      stop 'ooooo'
      
c now that we know the weights and boundaries, off we go!!!
c remember pV = nRT ==> p(z) dz/ r T(z) = dn(z)/V = dq(z) ==> Q = sum(p Z / R T)
c so for these fractional combined layers (i), Qnew = sum(p(i) zfrac(i) / R T(i)) = sum(p(i) zfrac(i)/Z(i) Z(i) / RT(i))
c                                                   = sum(p(i)Z(i)/RT(i) zfrac(i)/Z(i))
c or Qnew = sum(frac(i) Q(i))

      DO iX = iNot,kProfLayer
        raRAmt(iX) = 0.0
        rPP = 0.0
	rPPWgt = 0.0
	rMR = 0.0
	rMolecules = 0.0
	rHeight = 0.0
	DO iY = iaBnd(iX,1),iaBnd(iX,2)-1
	  IF (iY .EQ. iaBnd(iX,2)-1) THEN
	    rFrac = raBndFrac(iX,2)
	  ELSEIF (iY .EQ. iaBnd(iX,1)) THEN
	    !! this also takes care of case when iY .EQ. iaBnd(iX,1) .EQ. iaBnd(iX,2)-1
	    rFrac = raBndFrac(iX,1)
	  ELSE
	    rFrac = 1.0
	  END IF
	  rHeight = rHeight + rFrac*DATABASELEVHEIGHTS(iY)
	  rMolecules = rMolecules + raR100Amt(iY)*rFrac*DATABASELEVHEIGHTS(iY)
	  rPP = rPP + raR100PartPress(iY)*rFrac
	  rPPWgt = rPPWgt + rFrac
	  rMR = rMR + raR100PartPress(iY)/raRPress(iY)	      
	  raRAmt(iX) = raRAmt(iX) + raR100Amt(iY)*rFrac
	END DO
	!! method 1
	raRPartPress(iX) = rPP/rPPWgt

c        !! method 2
c	rMR = rMR/((iaBnd(iX,2)-1)-(iaBnd(iX,1))+1)	    
c      	raRPartPress(iX) = rMR * raPressLevels(iX)/1013.25
c	raRAmt(iX) = rMolecules/rHeight

c bumping raRAmt and raRPartPressup n down
c proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amout
c though WV may be a little more complicated as it depends on pp
c        raRAmt(iX) = raRAmt(iX) * 100.0
c        raRPartPress(iX) = raRPartPress(iX) * 20.0
c this proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amout
c though WV may be a little more complicated as it depends on pp	
c	write(*,1234) iGasID,iX,raPressLevels(iX),raaPress(iX,1)*1013.25,raPressLevels(iX+1),iaBnd(iX,1),iaBnd(iX,2),
c     $     	raRTemp(iX),raRPartPress(iX),raRAmt(iX)

      END DO
 1234 FORMAT(2(' ',I3),3(' ',F10.3),2(' ',I3),3(' ',E10.3))

c      IF (iGasID .EQ. 2) THEN
c        DO iL = 1, 100
c	  print *,iL,raR100Amt(iL),raRAmt(iL)
c	END DO
c      END IF
      
      RETURN
      END

c************************************************************************
c this subroutine finds the partial pressure of the layer
c especially useful for GasID = 1
      SUBROUTINE PPThruLayers(
     $     iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press,
     $     raDataBaseThickness,raSortPressLevels,raSortPressHeights,
     $     raPressLevels,raRTemp,r)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
c input variables
      INTEGER iLowerOrUpper   !!upper or lower atm
      INTEGER iI              !!layer of current profile we are looking at
      INTEGER iGasID
c these are the individual reference profiles, at kMaxLayer layers 
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer) 
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer) 
c these are the kCARTA database pressure levels and heights
      REAL raDataBaseThickness(kMaxLayer)  !!AIRS layer thickness (unsorted) m
      REAL raSortPressLevels(kMaxLayer+1)  !!sorted with increasing p        mb
      REAL raSortPressHeights(kMaxLayer+1) !!sorted with increasing p        m
c these are current profile pressure levels
      REAL raPressLevels(kProfLayer+1)
c this is the interpolated reference temperature
      REAL raRTemp(kProfLayer)
c output variable
      REAL r    !!the partial pressure

c local vars
      INTEGER iL,iLp1, iJ,iJp1,iMMM,iCase
      REAL r0,rH1,rH2,rP,rPp1,rHH,rPP,rDeltaH,r1,rNum,rDenom,rNum1,rDenom1,rXX

      REAL DatabaseHeight(kMaxLayer)
      REAL DATABASELEVHEIGHTS(kMaxLayer+1)
      REAL DATABASELEV(kMaxLayer+1)

      CALL databasestuff(iLowerOrUpper,
     $                   DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

      iMMM = kMaxLayer + 1
      r0 = r               !!! partial pressure from naive interpolation

      iL   = iI            !!!lower pressure level of layer iI
      iLp1 = iI + 1        !!!upper pressure level of layer iI
      rP   = raPressLevels(iL)
      rPp1 = raPressLevels(iL+1)

      iJ = 1
      !!!find databaselev(iJ) <= rP by looping up from GND
 10   CONTINUE
      IF (DATABASELEV(iJ) .GT. rP) THEN
        iJ = iJ + 1
        IF (iJ .LT. kMaxLayer+1) THEN
          GOTO 10
        END IF
      END IF
      
      !!check to see if we are OK or if we went one too far; if so, go down one
      IF ((DATABASELEV(iJ) .LT. rP) .AND. (iJ .GT. 1)) iJ = iJ - 1

      iJp1 = kMaxLayer + 1
      !!!find databaselev(iJp1) >= rPp1 by looping down from TOA
 20   CONTINUE
      IF (DATABASELEV(iJp1) .LT. rPp1) THEN
        iJp1 = iJp1 - 1
        IF (iJp1 .GT. 1) THEN
          GOTO 20
        END IF
      END IF        
      !!check to see if we are OK or if we went one too far; if so, go up one
      IF (DATABASELEV(iJp1) .GT. rPp1) iJp1 = iJp1 + 1

c      print *,iI,' : ',rP,iJ,DATABASELEV(iJ),'; ',rPp1,iJp1,DATABASELEV(iJp1)

c      print *,raR100Press(iJ),raR100PartPress(iJ),raR100Temp(iJ),
c     $        raR100Amt(iJ),raDataBaseThickness(iJ)
c      r = raR100Amt(iJ)*kAvog/(raDataBaseThickness(iJ) * 100) !!!molecules/cm3
c      print *, raR100PartPress(iJ)*kAtm2mb*100,
c     $       r*kBoltzmann*raR100Temp(iJ)*1e6,
c     $       raR100PartPress(iJ)*kAtm2mb*100/(r*kBoltzmann*raR100Temp(iJ)*1e6)

      IF ((iJp1 - iJ) .LE. 0) THEN 
        !!something wrong; we need lower DATABASE level < upper DATABASE level
        write(kStdErr,*) 'doing water amount ref partial pressure : '
        write(kStdErr,*) 'Layer iI : found iJ = iJp1 ',iI,iJ
        CALL DoStop
      END IF

      IF (iJ .LE. 0) THEN
        write(kStdErr,*) '>>layer number ',iI,'gas ID ',iGasID
        Call DoStopMesg(' >>oops iJ <= 0 in PPThruLayers, so will have problems in next line$')
      END IF
    
      !!! case 1 : current layer is INSIDE one database layer
      IF ((iJp1 - iJ) .EQ. 1) THEN       
        iCase = 1
        rDeltaH = (rPp1-rP)/(DATABASELEV(iJ+1)-DATABASELEV(iJ))
        rNum   = raR100Amt(iJ)*rDeltaH !! kmol/cm2
        rDenom = raDataBaseThickness(iJ)*rDeltaH   !! cm
        r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
        r = r*1e6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
        r = r/kAtm2mb                                          !!atm
      END IF

      !!! case 2 : current layer straddles one database level 
      !!!          ie partially occupies two database layers
      IF ((iJp1 - iJ) .EQ. 2) THEN
        iCase = 2

        !!! initialise with contribution from lower layer        
        rDeltaH = (DATABASELEV(iJ+1)-rP)/
     $            (DATABASELEV(iJ+1)-DATABASELEV(iJ))
        rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
        rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

        !!! add on contribution from upper layer        
        rDeltaH = (rPp1-DATABASELEV(iJ+1))/
     $            (DATABASELEV(iJ+2)-DATABASELEV(iJ+1))
        rNum1   = raR100Amt(iJ+1)*rDeltaH                       !! kmol/cm2
        rDenom1 = raDataBaseThickness(iJ+1)*rDeltaH             !! frac

        rNum   = rNum + rNum1
        rDenom = rDenom + rDenom1
        r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
        r = r*1e6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
        r = r/kAtm2mb                                          !!atm
      END IF

      !!! case 3 : current layer straddles more than one database level 
      !!!          ie partially occupies 2 database layers, and some full ones
      IF ((iJp1 - iJ) .GT. 2) THEN
        iCase = 3

        !!! initialise with contribution from lowest layer        
        rDeltaH = (DATABASELEV(iJ+1)-rP)/
     $            (DATABASELEV(iJ+1)-DATABASELEV(iJ))
        rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
        rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

        !!! add on contributions from intermediate layers
        DO iL = iJ+1,iJp1-2
          rDeltaH = 1.0
          rNum = rNum + raR100Amt(iL)*rDeltaH                   !!kmol/cm2
          rDenom = rDenom + raDataBaseThickness(iL)
        END DO

        !!! add on contribution from highest layer        
        rDeltaH = (rPp1-DATABASELEV(iJp1-1))/
     $            (DATABASELEV(iJp1)-DATABASELEV(iJp1-1))
        rNum1   = raR100Amt(iJp1-1)*rDeltaH                         !! kmol/cm2
        rDenom1 = raDataBaseThickness(iJp1-1)*rDeltaH               !! frac

        rNum   = rNum + rNum1
        rDenom = rDenom + rDenom1
        r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
        r = r*1e6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
        r = r/kAtm2mb                                          !!atm
      END IF

c      print *,iI,'(',iCase,')',r0,r,r/r0,raRTemp(iI)
c      print *,'------------------------------'
c      print *,' '

      RETURN
      END

c************************************************************************
c this function finds the pressure layer at which rPressStart is within,
c as well as the fraction of the layer that it occupies
      REAL FUNCTION FindBottomTemp(rP,raProfileTemp,
     $                             raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels = actual pressure levels that come out of kLAYERS
c raProfileTemp = actual profile temp
c rP            = pressure at which we want the temperature
      real rP,raProfileTemp(kProfLayer),raPressLevels(kProfLayer+1)
      integer iProfileLayers

      integer iFound,i1,i2,i3,iLowest,iJ
      real rP1,rP2,T1,T2
      real raP(3),raT(3),Y2A(3),rT,raLogP(3)
      real yp1,ypn,work(3)
      INTEGER iLog,iSpline

      iLog = +1       !!!do log(P) for the x-points 
      iLog = -1       !!!do   P    for the x-points 

      iSpline = -1    !!!use linear interpolations
      iSpline = +1    !!!use spline interpolations

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iLog = +1       !!!do log(P) for the x-points 
      iSpline = -1    !!!use linear interpolations
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      rT=0.0

      iLowest = kProfLayer-iProfileLayers+1

      IF (rP .ge. raPressLevels(iLowest)) THEN
        !this is WHOLE of the bottom layer
        i1=iLowest
      ELSE IF (rP .le. raPressLevels(kProfLayer+1)) THEN
        !this is ludicrous
        write(kStdErr,*) rP,raPressLevels(kProfLayer+1)
        write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
        CALL DoStop

      ELSE
        !first find the AIRS layer within which it lies
        iFound=-1 
        i1 = iLowest
        i2 = iLowest+1
 10     CONTINUE 
        IF ((rP .LE.raPressLevels(i1)).AND.(rP .GT.raPressLevels(i2))) THEN 
          iFound=1 
        END IF 
        IF ((iFound .LT. 0) .AND. (i1 .LT. kProfLayer)) THEN 
          i1=i1+1 
          i2=i2+1 
          GO TO 10 
        END IF 
        IF ((iFound .LT. 0)) THEN 
          IF (abs(rP-raPressLevels(kProfLayer+1)) .LE. delta) THEN 
            i1=kProfLayer 
            iFound=1 
          ELSE 
            write(kStdErr,*) 'could not find pressure ',rP 
            write(kStdErr,*) 'within AIRS pressure levels. Please check' 
            write(kStdErr,*) '*RADNCE and *OUTPUT sections' 
            CALL DoSTOP 
          END IF 
        END IF 
      END IF

      IF ((i1 .gt. kProfLayer) .OR. (i1 .lt. iLowest)) THEN
        write(kStdErr,*) 'sorry : cannot find surface temp for '
        write(kStdErr,*) 'layers outside ',iLowest,' and ',kProfLayer
        write(kStdErr,*) 'Allowed Pressure ranges are from : ',
     $ raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
        write(kStdErr,*) 'Surface Pressure is ',rP,' mb'
        call DoStop
      END IF 
          
      !now find the temperature
      IF (i1 .EQ. iLowest) THEN          !do linear interp
        i1 = iLowest
        i2 = iLowest+1
        i3 = iLowest+2
        rP1 = (raPressLevels(i2)-raPressLevels(i1))/
     $           log(raPressLevels(i2)/raPressLevels(i1))
        rP2 = (raPressLevels(i3)-raPressLevels(i2))/
     $           log(raPressLevels(i3)/raPressLevels(i2))
        T1 = raProfileTemp(i1)
        T2 = raProfileTemp(i2)
        IF (iLog .EQ. -1) THEN
          rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)           !!linear in P
        ELSE
          rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1)) !!log(P)
        END IF

      ELSEIF (i1 .GE. (kProfLayer-1)) THEN          !do linear interp
        rP1 = (raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/
     $          log(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
        rP2 = (raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/
     $          log(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
        T1 = raProfileTemp(kProfLayer-1)
        T2 = raProfileTemp(kProfLayer)
        IF (iLog .EQ. -1) THEN
          rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)            !!linear in P
        ELSE 
          rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1))  !!log(P)
        END IF
      ELSE          !do spline ... note that the pressures have to 
                    !be in ascENDing order for good interpolation
        rP1 = (raPressLevels(i1)-raPressLevels(i1-1))/
     $          log(raPressLevels(i1)/raPressLevels(i1-1))
        raP(3) = rP1
        rP1 = (raPressLevels(i1+1)-raPressLevels(i1))/
     $          log(raPressLevels(i1+1)/raPressLevels(i1))
        raP(2) = rP1
        rP1 = (raPressLevels(i1+2)-raPressLevels(i1+1))/
     $          log(raPressLevels(i1+2)/raPressLevels(i1+1))
        raP(1) = rP1
        IF (iLog .EQ. +1) THEN
          DO iJ = 1,3
            raLogP(iJ) = log(raP(iJ))
          END DO
        END IF

        raT(3) = raProfileTemp(i1-1)
        raT(2) = raProfileTemp(i1)
        raT(1) = raProfileTemp(i1+1)

        yp1=1.0e30
        ypn=1.0e30
        IF (iSpline .EQ. +1) THEN
          IF (iLog .EQ. +1) THEN
            CALL rspl(raLogP,raT,3,log(rP),rT,1) 
          ELSE
            CALL rspl(raP,raT,3,rP,rT,1)
          END IF
        ELSEIF (iSpline .EQ. -1) THEN
          IF (iLog .EQ. +1) THEN
            CALL rlinear(raP,raT,3,rP,rT,1)
          ELSE
            CALL rlinear(raLogP,raT,3,log(rP),rT,1)
          END IF
        END IF
      END IF

      FindBottomTemp = rT

      RETURN
      END

c************************************************************************
c this is called by kcartamain and kcartabasic, to set up the profiles
      SUBROUTINE Set_Ref_Current_Profs(
     $      iJax,rDerivTemp,rDerivAmt,
     $      iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,
     $                   raaAmt,raaTemp,raaPress,raaPartPress,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,
     $                   raTAmt,raTTemp,raTPress,raTPartPress,
     $                   raNumberDensity,pProfNLTE,raMixVertTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params
      INTEGER iJax                   !! when testing jacobians
      REAL rDerivTemp,rDerivAmt      !! when testing jacobians
      INTEGER iGas,iaGases(kMaxGas)  !! gasID stored in order they were read in
c these are the reference profiles stored in matrices
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
      REAL raaRPress(kProfLayer,kGasStore)
      REAL raaRPartPress(kProfLayer,kGasStore)
c these are the user specified layer profiles stored in matrices
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
      REAL raMixVertTemp(kMixFilRows)

c output params
c these are the individual reference profiles, at kProfLayer layers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c these are the user specified layer profiles
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPartPress(kProfLayer),raTPress(kProfLayer)
c this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
c when multiplied by the layer height, gives units of /cm^2
c so when multiplied by Rayleigh scattering cross section, this gives units
c of optical depth (no units)
      REAL raNumberDensity(kProfLayer)
      REAL pProfNLTE(kProfLayer),rDummy

c local vars
      INTEGER iInt

c get the reference profile for the current gas if GAS ID <= kGasXsecHi
      IF ((iaGases(iGas) .LE. kGasXsecHi) .OR. 
     $   (iaGases(iGas) .EQ. kNewGasHi+1)) THEN
        CALL SetReference(raRAmt,raRTemp,raRPress,raRPartPress,
     $            raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
      END IF

cjacob
c get actual profiles for the current gas
cc these are set in kcartamain, kcartabasic
cc          rDerivTemp = 0.01
cc          rDerivTemp = 0.1
cc          rDerivTemp = 0.5
cc          rDerivTemp = 1.0
cc          rDerivAmt  = 0.01
cc          rDerivAmt  = 0.1

c     print *,'should I use pProf or profNLTE in NLTE routines???'
c     print *,'looks like i need to use pprofNLTE for Dave Edwards profiles
c     print *,'but can get away with pprof using klayers stuff'
c     print *, ie raPressLevels is consistent with pProf (avg layer press)
c             but it comes from summing partial pressures of MANY gases
c     print *, while for Dave Edwards tests, he only uses one gas (CO2)
c             and so mebbe the pressure levels are not consistent with pProf
c            print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
c     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)

c      Do iInt = 1,80
c        print *,iInt,raaTemp(90,iInt)
c      end do
c      call dostopmesg('stop in Set_Ref_Current_Profs$')
      
      DO iInt=1,kProfLayer
        raTAmt(iInt)          = raaAmt(iInt,iGas)
        raTTemp(iInt)         = raaTemp(iInt,iGas)
        raTPress(iInt)        = raaPress(iInt,iGas)
        raTPartPress(iInt)    = raaPartPress(iInt,iGas)
        !!compute particle number density in number/cm3 (N/V = P/kT)
        raNumberDensity(iInt) = raTPress(iInt)*kAtm2mb*100.0/(kBoltzmann * 
     $                              raTTemp(iInt))*1e-6
        !!NLTE avg layer press in lower atm = kCARTA profile 
        pProfNLTE(iInt) = raTPress(iInt)*kAtm2mb   

c       print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
c     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cc  this stuff is to test the jacobians
c          raTTemp(iInt)=rDerivTemp+raTTemp(iInt)
c          raMixVertTemp(iInt)=rDerivTemp+raTTemp(iInt)
c          raMixVertTemp(iInt)=raTTemp(iInt)

cc jacob test temperature jacobian ... old
c          IF (iInt .EQ. iJax) THEN
c            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
c            raMixVertTemp(iInt)=raTTemp(iInt)
c          END IF
c          IF ((iInt .EQ. iJax) .AND. (iGas .EQ. 1)) THEN
c            rDummy        = raMixVertTemp(iInt)
c            raMixVertTemp(iInt) = rDummy+rDerivTemp
ccc          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
ccc          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
c            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
c          END IF

cc jacob test column gas  jacobian
c          IF (iGas .EQ. 1) THEN
c            rDummy        = raMixVertTemp(iInt)
c            raMixVertTemp(iInt) = rDummy+rDerivTemp
ccc          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
ccc          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
c            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
c          END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cc jacob test jacobian ... new, could also be column
c          iJax = -1
c          iJax = 65
c          rDerivTemp = 1.0
c          IF ((iInt .EQ. iJax) .OR. (iJax .EQ. -1)) THEN
c            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
c            raMixVertTemp(iInt)=raTTemp(iInt)
c            print *,iJax,raaTemp(iInt,iGas),rDerivTemp,raMixVertTemp(iInt)
c          END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
cc jacob test amount jacobian
c         iJax = 65
c         IF ((iInt .EQ. iJax).AND.(iaGases(iGas) .EQ. 2)) THEN
c           raTAmt(iInt)=raaAmt(iInt,iGas)*(1.0+rDerivAmt)
c           print *,'iJax, dq = ',iJax,raaAmt(iInt,iGas)*rDerivAmt
c         END IF
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

c Matlab code to plot the jacs
c suppose you have a 55 layer atm (from 46 to 100) and a 1 k temperature perturbation in layer 65
c   (which is layer 20 of the Radiationg Atmosphere)
c [d0,w] = readkcstd('junk0.dat');
c [jac,w] = readkcjac('junk.jac');
c [dt,w] = readkcstd('junk.dat');
c [m,n] = size(jac); iL = (n-4)/3
c
c plot(w,(dt-d0)/1,'b.-',w,jac(:,20+55),'r')
c
c [fc,qc] = quickconvolve(w,jac,0.25,0.25);
c pcolor(fc,1:iL,qc(:,(1:iL)+iL*0)'); shading flat; colorbar  %% Q jacobian
c pcolor(fc,1:iL,qc(:,(1:iL)+iL*1)'); shading flat; colorbar  %% T jacobian
c pcolor(fc,1:iL,qc(:,(1:iL)+iL*2)'); shading flat; colorbar  %% WGT fcn

      END DO

      RETURN
      END

c************************************************************************
c this figures out CO2 mixing ratio
c since Antartic surface pressures can be as low as 500 mb, CO2.N2O,CO mixing ratios were failing if iBot=20
c this corresponds to raPresslevls(20) = 596 mb
c so subr Get_Temp_Plevs (in n_pth_mix.f) and subr compute_co2_mixratio (in kcartamisc.f)
c both needed to have iBot tweaked to layer 25 or above, which corresponds to raPresslevls(25) = 496 mmb
      SUBROUTINE compute_co2_mixratio(raaPress,raaPartPress,raaAmt,iNumLayer,rFracBot,rCO2MixRatio)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iNumLayer
      REAL raaPress(kProfLayer,kGasStore),rFracBot
      REAL raaPartPress(kProfLayer,kGasStore),raaAmt(kProfLayer,kGasStore)
c output
      REAL rCO2MixRatio

c local
      INTEGER iI,iBot,iTop,iCount,i1,i2,iLowestLay
      REAL rX,rF,cfac
      
      ! Conversion factor (equivalent mm of liquid water per molecules/cm^2)
      ! cfac = 1 (cm^3/gram) * ~18/Navadagro (grams/molecule) * 10 (mm/cm)
      cfac = 10 * 18.015 / (kAvog / 1000)

c >>>>>>>>>>>>>>>>>>>>>>>>>
      iBot = kProfLayer - iNumLayer + 1 + 1  !! assume profile has started from here!!!!
      iTop = iBot + 10
      write(kStdWarn,*) 'WATER : Compute lower trop <ppmv> between lays ',iBot,iTop

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,1)/raaPress(iI,1)
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rWVMix = ',rX,' ppmv'

      iLowestLay = kProfLayer - iNumLayer + 1
      rX = 0.0
      iCount = 0
      DO iI = iLowestLay,kProfLayer
        rF = 1.0
        IF (iI .EQ. iLowestLay) rF = max(min(rFracBot,1.0),0.0)
        iCount = iCount + 1
        rX = rX + raaAmt(iI,1) * rF
      END DO
      rX = rX * cfac * (kAvog)   !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
      write(kStdWarn,*)' column amount water = ',rX,' mm'
     
c >>>>>>>>>>>>>>>>>>>>>>>>>
      iBot = 20                              !! assume profile has started from here !!!!
      iBot = max(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!            
      iTop = kProfLayer - iBot
      write(kStdWarn,*) 'CO2 : Compute lower trop <ppmv> between lays ',iBot,iTop      

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,2)/raaPress(iI,2)
      END DO
      rX = rX * 1.0e6/iCount
      rCO2MixRatio = rX
      write(kStdWarn,*)'avg rCO2Mix = ',rX,' ppmv'

c >>>>>>>>>>>>>>>>>>>>>>>>>
      iBot = 40  !! assume profile has started from here!!!!
      iTop = 70
      write(kStdWarn,*) 'O3 : Compute lower trop <ppmv> between lays ',iBot,iTop      

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,3)/raaPress(iI,3)
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rO3Mix = ',rX,' ppmv'

      ! Conversion factor
      ! Note: a Dobson unit is 10^-5 meters at 1 atm and 273.15 K
      ! cfac = (1/loschmidt) * (1000 du/cm)
      ! Loschmidt (molecules per cm^3 at 1 atm and 273.15 K)
      cfac=1000/2.6867775E+19

      iLowestLay = kProfLayer - iNumLayer + 1
      rX = 0.0
      iCount = 0
      DO iI = iLowestLay,kProfLayer
        rF = 1.0
        IF (iI .EQ. iLowestLay) rF = rFracBot
        iCount = iCount + 1
        rX = rX + raaAmt(iI,3) * rF
      END DO
      rX = rX * cfac * (kAvog) !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
      write(kStdWarn,*)' column amount ozone = ',rX,' du'

c >>>>>>>>>>>>>>>>>>>>>>>>>
      iBot = 20                              !! assume profile has started from here !!!!
      iBot = max(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!            
      iTop = kProfLayer - iBot
      write(kStdWarn,*) 'N2O/CO/CH4 : Compute lower trop <ppmv> between lays ',iBot,iTop            

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,4)/raaPress(iI,4)
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rN2OMix = ',rX,' ppmv'

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,5)/raaPress(iI,5)
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rCOMix = ',rX,' ppmv'

      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,6)/raaPress(iI,6)
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rCH4Mix = ',rX,' ppmv'

      write(kStdWarn,*) 'assumed fracBot = ',rFracBot,max(min(rFracBot,1.0),0.0)
      write(kStdWarn,*) ' '

      RETURN
      END

c************************************************************************
c this duplicates clear sky atmospheres!
        SUBROUTINE duplicate_clearsky_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c these are the "output" variables
      INTEGER iAtmLoop,iNatm
      REAL raAtmLoop(kMaxAtm)
c these are the "input variables"
      REAL raPressLevels(kProfLayer+1)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      INTEGER iaMPSetForRad(kMaxAtm),iProfileLayers
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iAtm                  !this is the atmosphere number
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
c iSetSolarRefl tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raaPrBdry(kMaxAtm,2)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      REAL rakSolarRefl(kMaxPts)
      REAL raSetEmissivity(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raLayerHeight(kProfLayer)
      REAL rJunk1,rJunk2
      
c local var
      INTEGER iX,iY
      REAL rX,SACONV_SUN,ORIG_SACONV_SUN,VACONV

c first find out how many raAtmLoop the user did set
      iX = 1
      DO WHILE ((iX .LE. kMaxAtm) .AND. (raAtmLoop(iX) .GE. 0.0))
        iX = iX + 1
	IF (iX .GT. kMaxAtm) GOTO 10
      END DO
 10   CONTINUE   
      iX = iX - 1 
      
      write(kStdWarn,*) ' >>>> Duplicate Clear Sky Params from Atm # 1 for ',iX,' atmospheres'
      iNatm = iX
      DO iX = 2,iNatm
        iaMPSetForRad(iX)      = iaMPSetForRad(1)

        iaSetEms(iX)           = iaSetEms(1)
        iaSetSolarRefl(iX)     = iaSetSolarRefl(1)
        caEmissivity(iX)       = caEmissivity(1)
        raSetEmissivity(iX)    = raSetEmissivity(1)
        DO iY = 1,kEmsRegions
          raaaSetEmissivity(ix,iY,1) = raaaSetEmissivity(1,iY,1)
          raaaSetEmissivity(ix,iY,2) = raaaSetEmissivity(1,iY,2)
          raaaSetSolarRefl(ix,iY,1)  = raaaSetSolarRefl(1,iY,1)
          raaaSetSolarRefl(ix,iY,2)  = raaaSetSolarRefl(1,iY,2)
        END DO
        iaKSolar(iX)           = iaKSolar(1)
        rakSolarAngle(iX)      = rakSolarAngle(1)
        rakSolarRefl(iX)       = rakSolarRefl(1)
        iaKThermal(iX)         = iaKThermal(1)
        rakThermalAngle(iX)    = rakThermalAngle(1)
        iakThermalJacob(iX)    = iakThermalJacob(1)
        iaSetThermalAngle(iX)  = iaSetThermalAngle(1)
        raSatHeight(iX)        = raSatHeight(1)
        raSatAngle(iX)         = raSatAngle(1)

        DO iY = 1,2
          raaPrBdry(iX,iY)        = raaPrBdry(1,iY)
        END DO
        raPressStart(iX)       = raPressStart(1)
        raPressStop(iX)        = raPressStop(1)

        raTspace(iX)           = raTSpace(1)
        raTSurf(iX)            = raTSurf(1)

        iaNumLayer(iX)         = iaNumLayer(1)
        DO iY = 1,kProfLayer
          iaaRadLayer(iX,iY)     = iaaRadLayer(1,iY)
        END DO

        iaLimb(iX)     = iaLimb(1)
        raFracTop(iX)  = raFracTop(1)
        raFracBot(iX)  = raFracBot(1)
      END DO

c now set the param you need to set
      IF (iAtmLoop .EQ. 1) THEN
        write(kStdWarn,*) '  Resetting raPressStart for looping'
        write(kStdErr,*)  '  Resetting raPressStart for looping'
        IF ((raaPrBdry(1,1) .GT. raaPrBdry(1,2)) .AND. (iaLimb(1) .LE. 0)) THEN
          write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
          write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
        ELSEIF ((raaPrBdry(1,1) .GT. raaPrBdry(1,2)) .AND. (iaLimb(1) .GT. 0)) THEN
          write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
          write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
        END IF
        IF (raaPrBdry(1,1) .LT. raaPrBdry(1,2)) THEN
          write(kStdWarn,*) '  ---> warning : reset TOA press for uplook instr is odd'
          write(kStdErr,*)  '  ---> warning : reset TOA press for uplook instr is odd'
          CALL DoStop
        END IF
        DO iX = 1,iNatm
          raPressStart(iX) = raAtmLoop(iX)
          raaPrBdry(iX,1)  = raAtmLoop(iX)
        END DO

      ELSEIF (iAtmLoop .EQ. 2) THEN
        write(kStdWarn,*) '  Resetting raPressStop for looping'
        write(kStdErr,*)  '  Resetting raPressStop for looping'
        IF (raaPrBdry(1,1) .LT. raaPrBdry(1,2)) THEN
          write(kStdWarn,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
          write(kStdErr,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
        END IF
        DO iX = 1,iNatm
          raPressStop(iX) = raAtmLoop(iX)
          raaPrBdry(iX,2)  = raAtmLoop(iX)
        END DO

      ELSEIF (iAtmLoop .EQ. 3) THEN	
        write(kStdWarn,*) '  Resetting raSatZen for looping (so need to compute scanang)'
        write(kStdErr,*)  '  Resetting raSatZen for looping (so need to compute scanang)'
        IF (iaLimb(1) .GT. 0) THEN
          write(kStdErr,*) 'Atm 1 set up for Limb sounding'
          write(kStdErr,*) '  so cannot willy nilly reset scanang'
          write(kStdErr,*) 'Go and reset raStartPress instead'
          CALL DoStop
        ELSE
          write(kStdWarn,*) '  changing user input SatZen (angle at gnd)  to       Instr ScanAng '
          write(kStdWarn,*) '  raAtmLoop(iX) --> raSatAngle(iX)     GNDsecant  --> SATELLITEsecant'

          IF (rSatHeightCom .LT. 0) THEN
  	    write(kStdWarn,*) '  WARNING : raSatHeight == -1 so kCARTA uses SATELLITEsecant!!!!'
            DO iX = 1,iNatm
              raSatAngle(iX) = raAtmLoop(iX)
  	      rJunk1 = 1.0/cos(raAtmLoop(iX)*kPi/180)
	      rJunk2 = 1.0/cos(raSatAngle(iX)*kPi/180)	    
	      write(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2	    
	    END DO
	  ELSE
            DO iX = 1,iNatm
              raSatAngle(iX) = raAtmLoop(iX)
              !!!! positive number so this is genuine input angle that will vary with layer height
              raSatAngle(iX) = SACONV_SUN(raAtmLoop(iX),0.0,705.0)
	      rJunk1 = 1.0/cos(raAtmLoop(iX)*kPi/180)
	      rJunk2 = 1.0/cos(raSatAngle(iX)*kPi/180)	    
	      write(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2	    
            END DO
	  END IF
        END IF
 111  FORMAT('   ',F10.5,' ---> ',F10.5,'   +++   ',F10.5,' ---> ',F10.5)
      
      ELSEIF (iAtmLoop .EQ. 4) THEN
        write(kStdWarn,*) '  Resetting raSolZen for looping'
        write(kStdErr,*)  '  Resetting raSolZen for looping'
        DO iX = 1,iNatm
          rakSolarAngle(iX) = raAtmLoop(iX)
        END DO

      ELSEIF (iAtmLoop .EQ. 5) THEN
        write(kStdWarn,*) '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
        write(kStdErr,*)  '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
        DO iX = 1,iNatm
          DO iY = 1,kEmsRegions
            raaaSetEmissivity(iX,iY,2) = raaaSetEmissivity(iX,iY,2) + raAtmLoop(iX)
            raaaSetSolarRefl(iX,iY,2)  = (1-raaaSetEmissivity(iX,iY,2))/kPi
          END DO
        END DO

      ELSEIF (iAtmLoop .EQ. 10) THEN
        write(kStdWarn,*) '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'
        write(kStdErr,*)  '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'

      ELSEIF (iAtmLoop .EQ. 100) THEN
        write(kStdWarn,*) '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'
        write(kStdErr,*)  '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'

      ELSE
        write(kStdErr,*) 'Dont know what to do with iAtmLoop = ',iAtmLoop
        Call DoStop
      END IF

      IF (iAtmLoop .LE. 2) THEN
        CALL Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer,
     $                             iProfileLayers,iaMPSetForRad,
     $                             raSatHeight,raSatAngle,raPressStart,raPressStop,
     $                             raFracTop,raFracBot,raPressLevels,raLayerHeight,
     $                             iakSolar,rakSolarAngle)
      END IF

      RETURN
      END

c ************************************************************************
c this subroutine resets iaaRadLayer and/or raSatAngle, if the Start or Stop 
c Pressures have changed
      SUBROUTINE Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer,
     $                             iProfileLayers,iaMPSetForRad,
     $                             raSatHeight,raSatAngle,raPressStart,raPressStop,
     $                             raFracTop,raFracBot,raPressLevels,raLayerHeight,
     $                             iakSolar,rakSolarAngle)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iNatm,iProfileLayers,iaMPSetForRad(kMaxAtm),iakSolar(kMaxAtm)
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm),raLayerHeight(kProfLayer)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm),rakSolarAngle(kMaxAtm)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm)
c raaPrBdry = pressure start/stop
      REAL raaPrBdry(kMaxAtm,2),raPressLevels(kProfLayer+1)

      INTEGER iC,iX,iStart,iStop,iNlay,iDirection,iInt
      REAL LimbViewScanAng

      DO iC = 1,iNAtm
        CALL StartStopMP(iaMPSetForRad(iC),raPressStart(iC),raPressStop(iC),iC,
     $           raPressLevels,iProfileLayers,
     $           raFracTop,raFracBot,raaPrBdry,iStart,iStop)        

        IF (iStop .GE. iStart) THEN
          iNlay = (iStop-iStart+1)
          iDirection = +1                           !down look instr
        ELSE IF (iStop .LE. iStart) THEN
          iNlay = (iStart-iStop+1)
          iDirection = -1                           !up look instr
        END IF
        IF (iNLay .GT. kProfLayer) THEN
          write(kStdErr,*)'Error for atm # ',iC
          write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
          CALL DoSTOP
        END IF
        iaNumlayer(iC) = iNlay

        DO iInt = 1,iNlay
          iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
        END DO

        write(kStdWarn,*) ' Atm#, Press Start/Stop, iStart,iStop, Nlay = ',
     $     iC,raPressStart(iC),raPressStop(iC),iStart,iStop,iNlay

      END DO

      IF (iaLimb(1). GT. 0) THEN
        !! this is limb sounder, do the angle etc
        DO iC = 1,iNatm
          raSatAngle(iC) = LimbViewScanAng(iC,raPressStart,raSatHeight,
     $                         iaaRadLayer,raPressLevels,raLayerHeight)
          IF (iaKsolar(iC) .GE. 0) THEN
            rakSolarAngle(iC) = raSatAngle(iC)  !! this is scanang at TOA instr
            rakSolarAngle(iC) = 89.9            !! this is sol zenith at "surface"
          END IF
        END DO
      END IF
      
      RETURN
      END
c ************************************************************************
c this duplicates cloud sky 2slab atmospheres!
        SUBROUTINE duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,
     $            cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c these are the "output" variables
      INTEGER iAtmLoop,iNatm
      REAL raAtmLoop(kMaxAtm)
c these are the "input variables"
      REAL raPressLevels(kProfLayer+1)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      INTEGER iaMPSetForRad(kMaxAtm),iProfileLayers
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iAtm                  !this is the atmosphere number
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
c iSetSolarRefl tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raaPrBdry(kMaxAtm,2)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      REAL rakSolarRefl(kMaxPts)
      REAL raSetEmissivity(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raLayerHeight(kProfLayer)

c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
      REAL rAngle
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
c this is info about cloud type, cloud frac
      INTEGER ctype1,ctype2
      REAL cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

c local var
      INTEGER iX,iY,iDebug
      REAL rX

      iDebug = +1
      iDebug = -1
      IF (iDebug .GT. 0) THEN

        print *,' ' 
        print *,'INITIAL Clouds Before duplications'
        print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

        print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
        print *,'ctype1,ctype2 = ',ctype1,ctype2
        print *,'iNclouds = ',iNclouds

        print *,'showing iaCloudNumAtm(iX) : '
        print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
        print *,' '

        print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
          print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
        END DO
        print *,' '

        !! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
        print *,'showing iaaScatTable'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
          print *,' '
        END DO
        print *,' '

        print *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
          print *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
          print *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
          print *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
          print *,' '
        END DO

        print *,' '
        IF (iCldProfile .GT. 0) THEN
          print*,'iCldProfile'
          print *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
          print *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
        END IF
      END IF

      !************************************************************************
      !!! now have to update things
      !!! if there are originally 2 clouds then
      !!!   we are going from 2 clouds in atmosphere #1 to adding on 
      !!!                       cloud1 in atmosphere #2 
      !!!                       cloud2 in atmosphere #3 
      !!!                    NO clouds in atmosphere #4
      !!!                    r5 = clr r4 + c1 r2 + c2 r3 + c12 c1 where clr = 1-c1-c2+c12
      !!! if there are originally 1 clouds then
      !!!   we are going from 1 clouds in atmosphere #1 to adding on 
      !!!                     0  cloud1 in atmosphere #2 
      !!!                     0  cloud2 in atmosphere #3 
      !!!                     O  clouds in atmosphere #4
      !!!                    r5 = clr r4 + c1 r1                  where clr = 1-c1
      !!!  IN OTHER WORDS no need to sweat anything if iNclouds ===== 1 YAYAYAYAYAYAYAYAYAYAYAYA

      IF (iCldProfile .GT. 0) THEN
        write(kStdErr,*) 'Ooops can only duplicate clouds slabs, not profiles'
        CALL DoStop
      END IF

      IF (iNclouds .EQ. 1) THEN
        write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
        !!! just duplicate the clear fields
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

      ELSEIF ((iNclouds .GT. 2) .OR. (iNclouds .LE. 0)) THEN
        write(kStdErr,*) 'iNclouds = ',iNclouds ,' huh?? cant duplicate this !!!'
        CALL DoStop
      ELSEIF (iNclouds .EQ. 2) THEN
        DO iX = 1,iNclouds
          iaCloudNumAtm(iX) = 2
        END DO

        ! no need to upgrade iaaCloudWhichLayers
        ! no need to upgrade iaaScatTable

        ! need to upgrade iaaCloudWhichAtm
        iY = 1
        iaaCloudWhichAtm(iY,2) = 2   !! this means cloud #1 is also going to be used in atm #2
        
        iY = 2
        iaaCloudWhichAtm(iY,2) = 3   !! this means cloud #1 is also going to be used in atm #3

        ! no need to upgrade raaPCloudTop,raaPCloudbot
        ! no need to upgrade raaaCloudParams (cloud loading and <dme>)
      
        !!! finally duplicate the clear fields
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
      END IF

      iDebug = -1
      IF (iDebug .GT. 0) THEN
        print *,' ' 
        print *,'FINAL CLOUD after duplications'

        print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

        print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
        print *,'ctype1,ctype2 = ',ctype1,ctype2
        print *,'iNclouds = ',iNclouds

        print *,'showing iaCloudNumAtm(iX) : '
        print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
        print *,' '

        print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
          print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
        END DO
        print *,' '
      END IF

      RETURN
      END

c************************************************************************
c this duplicates cloud sky 100slab atmospheres!
        SUBROUTINE duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,
     $            cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c these are the "output" variables
      INTEGER iAtmLoop,iNatm
      REAL raAtmLoop(kMaxAtm)
c these are the "input variables"
      REAL raPressLevels(kProfLayer+1)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      INTEGER iaMPSetForRad(kMaxAtm),iProfileLayers
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iAtm                  !this is the atmosphere number
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
c iSetSolarRefl tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raaPrBdry(kMaxAtm,2)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      REAL rakSolarRefl(kMaxPts)
      REAL raSetEmissivity(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raLayerHeight(kProfLayer)

c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
      REAL rAngle
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
c this is info about cloud type, cloud frac
      INTEGER ctype1,ctype2
      REAL cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

c local var
      INTEGER iX,iY,iDebug
      REAL rX

      iDebug = +1
      iDebug = -1
      IF (iDebug .GT. 0) THEN

        print *,' ' 
        print *,'INITIAL Clouds Before duplications'
        print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

        print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
        print *,'ctype1,ctype2 = ',ctype1,ctype2
        print *,'iNclouds = ',iNclouds

        print *,'showing iaCloudNumAtm(iX) : '
        print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
        print *,' '

        print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
          print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
        END DO
        print *,' '

        !! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
        print *,'showing iaaScatTable'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
          print *,' '
        END DO
        print *,' '

        print *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
        DO iY = 1,iNclouds
          print *,'Cloud ',iY
          print *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
          print *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
          print *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
          print *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
          print *,' '
        END DO

        print *,' '
        IF (iCldProfile .GT. 0) THEN
          print*,'iCldProfile'
          print *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
          print *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
        END IF
      END IF

      !************************************************************************
      !!! now have to update things
      !!! if there are originally 2 clouds then
      !!!   just do one 100 layer ice/water cloud and one clear calc

      IF (iCldProfile .LT. 0) THEN
        write(kStdErr,*) 'Ooops can only duplicate 100 layer cloud profiles, not slabs'
        CALL DoStop
      END IF

      write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
      !!! just duplicate the clear fields
      CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

      RETURN
      END

c************************************************************************
c funtion to estimate height (in m), given pressure (in mb)
c  based on US STD atm
      REAL FUNCTION p2h(p)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/airsheights.param'
      include '../INCLUDE/airslevels.param'

      INTEGER iI
      REAL pavg,rH,raY2P(kMaxLayer),p,logpavg(kMaxLayer),raHgt(kMaxLayer)
     
      DO iI = 1,100
        pavg = (DATABASELEV(iI+1)-DATABASELEV(iI))/log(DATABASELEV(iI+1)/DATABASELEV(iI))
        logpavg(kMaxLayer-iI+1) = log(pavg)      !! need this to be smallest to largest
        raHgt(kMaxLayer-iI+1)   = DatabaseHEIGHT(iI)
      END DO

c      print *,(logpavg(iI),iI=1,kMaxLayer)
c      print *,(raHgt(iI),iI=1,kMaxLayer)

      IF (p .GE. DATABASELEV(1)) THEN 
        rH = DatabaseHEIGHT(1)
      ELSEIF (p .LE. DATABASELEV(kMaxLayer+1)) THEN 
        rH = DatabaseHEIGHT(kMaxLayer)
      ELSE
        CALL rlinear_one(logpavg,raHgt,kMaxLayer,log(p),rH)
        CALL rspl(logpavg,raHgt,kMaxLayer,log(p),rH,1)
      END IF

      p2h = rH

      RETURN
      END

c************************************************************************
