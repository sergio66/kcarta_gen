! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

!! use must precede IMPLICIT and INCLUDE statements

MODULE kcartamisc

USE basic_common
USE freqfile
USE spline_and_sort_and_common
USE n_main
USE s_misc
USE rad_misc
USE rad_angles
USE n_rad_jac_scat

IMPLICIT NONE

!************************************************************************
!******** THIS FILE CONTAINS VARIOUS USEFUL SUBROUTINES/FUNCTIONS *******
!** such as sorting, setting vertical temperature profiles, checking ****
!** kcartaparam.f90, checking comp.param and xsec.param, splines etc *******
! note xisnan2, xisinf2, xisinfinite2 are in n_pth_mix.f90
!************************************************************************

CONTAINS
     
!************************************************************************

    LOGICAL FUNCTION isfinite(a)
    REAL :: a
    isfinite = (a-a) == 0
    end FUNCTION isfinite

!************************************************************************
! check for isnan
    LOGICAL FUNCTION isnan_real(x)

    REAL :: x
    LOGICAL :: isnan

!    print *,x,x+1.0e0
    IF (x+1.0e0 == x) THEN
      isnan = .TRUE. 
    ELSE
      isnan = .FALSE. 
    ENDIF
    isnan = .TRUE. 

!    IF (x >= 0.0) THEN
!        print *,'gt',x
!        isnan = .FALSE. 
!    ELSEIF (x <= 0.0) THEN
!        print *,'lt',x
!        isnan = .FALSE. 
!    ENDIF

    isnan_real = isnan
    RETURN
    end FUNCTION isnan_real

!************************************************************************
    LOGICAL FUNCTION isnan_double(x)

    logical :: isnan
    double precision :: x

!    print *,x,x+1.0d0
    IF (x+1.0d0 == x) THEN
      isnan = .TRUE. 
    ELSE
      isnan = .FALSE. 
    ENDIF

    isnan_double = isnan
    RETURN
    end FUNCTION isnan_double

!************************************************************************
! this subroutine sets the kCompressed database uncompression type
    SUBROUTINE SetSplineType(raPresslevels,iProfileLayers, &
    iNumNLTEGases,iNLTE_SlowORFast,iSplineType)
          
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'

! input vars
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iNumNLTEGases,iNLTE_SlowORFast
! ouput vars
! iSplinetype = 2 if FAST interp can be done
    INTEGER :: iSplineType

! local vars
    INTEGER :: iDefault,iJ,iUpDn,iMax
    REAL :: rDelta,rMin,rMax,rP,raFracDeltaP(kMaxLayer)
    REAL :: raFracDeltaPav(kMaxLayer),rAv1,rAv2,r1,r2,rMaxAv

! this assumes arbitrary pressure layering; can be slow, as kCARTA needs to
! interpolate the (pressure levels, abs coeffs) of the database onto the klayers
! pressure layers, before doing the temperature interpolation
! note abs coeff = stored optical depth/default gas amount
    iDefault    = +1        !!! Spline  .... DEFAULT
    iSplineType = -1        !!! Linear
    iSplineType = -2        !!! Matlab uncompression (linear weights)
    iSplineType = +2        !!! Matlab uncompression (linear weights)
    iSplineType = +1        !!! Spline  .... DEFAULT
    iSplineType = iaaOverrideDefault(1,2)
    IF ((abs(iSPlineType) /= 1) .AND. (abs(iSPlineType) /= 2)) THEN
      write(kStdErr,*) 'invalid iSplineType = ',iSplineType
      CALL DoStop
    END IF
    IF ((iSplineType /= iDefault) .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*)  'iSplineType,iDefault = ',iSplineType,iDefault
      write(kStdWarn,*) 'iSplineType,iDefault = ',iSplineType,iDefault
    END IF

! recall kMaxLayer = 100 = the kCARTA adtabase
!        kProfLayer = N  = what klayers was compiled for (hopefully)
!      IF ((iProfileLayers .LE. kMaxLayer) .AND.
!     $     ((iNumNLTEGases .LE. 0) .OR. (iNLTE_SlowORFast .LT. 0))
!     $     .AND. (kMaxLayer .LE. kProfLayer)) THEN
    IF ((iProfileLayers <= kMaxLayer) .AND. &
      ((iNumNLTEGases <= 0) .OR. (iNLTE_SlowORFast == -1)) &
     .AND. (kMaxLayer <= kProfLayer)) THEN
      ! quite possible that the pressure levels are same as kCARTA database
      ! in which case kcoeffSPL and kcoeffSPLandJAC can be sped up, as we can
      ! straightaway do tempr interpolation of stored abs coeffs without
      ! having to worry about doing the pressure interpolation as well
      ! note abs coeff = stored optical depth/default gas amount
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
        raFracDeltaP(iJ) = abs(raPresslevels(iJ)-PLEV_KCARTADATABASE_AIRS(iJ))/PLEV_KCARTADATABASE_AIRS(iJ)
        r1 = PLEV_KCARTADATABASE_AIRS(iJ)-PLEV_KCARTADATABASE_AIRS(iJ+1)
        r2 =log(PLEV_KCARTADATABASE_AIRS(iJ)/PLEV_KCARTADATABASE_AIRS(iJ+1))
        rAv1 = r1/r2
        r1 = raPresslevels(iJ)-raPresslevels(iJ+1)
        r2 = log(raPresslevels(iJ)/raPresslevels(iJ+1))
        rAv2 = r1/r2
        raFracDeltaPav(iJ) = abs(rAv1-rAv2)/rAv1
        write(kStdWarn,100) iJ,PLEV_KCARTADATABASE_AIRS(iJ),raPresslevels(iJ),raFracDeltaP(iJ),rAv1,rAv2,raFracDeltaPav(iJ)
        rDelta = rDelta + raFracDeltaP(iJ)
        IF (rMin > raFracDeltaP(iJ)) THEN
          rMin = raFracDeltaP(iJ)
        END IF
        IF (rMax < raFracDeltaP(iJ)) THEN
          rMax = raFracDeltaP(iJ)
          iMax = iJ
          rP = PLEV_KCARTADATABASE_AIRS(iJ)
        END IF
        IF (rMaxAv < raFracDeltaPav(iJ)) THEN
          rMaxAv = raFracDeltaPav(iJ)
        END IF
      END DO
      rDelta = rDelta/iProfileLayers
      write(kStdWarn,*) 'difference between kCARTA database and klayers ...'
      IF ((rMin <= 1.0e-5) .AND. (rMax <= 5.0e-4) .AND. (rMaxAv <= 5.0e-3)) THEN
        write(kStdWarn,*) '  can use kCARTA database levels for lower atm'
        iSplineType = iSplineType * 2  !!so -1 becomes -2, or +1 becomes +2
      ELSE
        write(kStdWarn,*) '  intrp kCARTA databse lvls for lower part of atm'
      END IF
      write(kStdWarn,*) 'MinDiff, MaxDiff, i(MaxDiff), Press(i(MaxDiff)), Sum(abs(diff))/Numlayers = '
      write(kStdWarn,*) rMin,rMax,iMax,rP,rDelta
      IF (abs(iSplineType) == 1) THEN
        write(kStdWarn,*) 'iSplineType = ',iSplineType, ' slow : interpolate'
        write(kStdWarn,*) 'database (pavglayers,abscoeeff) onto new plevels'
        write(kStdWarn,*) 'before doing the temp interpolation'
        write(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
      ELSEIF (abs(iSplineType) == 2) THEN
        write(kStdWarn,*) 'iSplineType = ',iSplineType, ' fast : use '
        write(kStdWarn,*) 'database (pavglayers,abscoeeff) for temp interp'
        write(kStdWarn,*) 'as presslevel scheme of klayers = kCARTA database!!'
        write(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
      ELSE
        write(kStdErr,*) 'need abs(iSplineType) = 1 or 2'
        CALL DoStop
      END IF
    END IF

  100 FORMAT(I3,' ',F10.5,' ',F10.5,' ',E10.5,' +++ ',F10.5,' ',F10.5,' ',E10.5)

    RETURN
    end SUBROUTINE SetSplineType

!************************************************************************
! this subroutine checks to see if the CO2 ppmv is ok
    SUBROUTINE check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix, &
    iaGases,rCO2mult)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input vars
! raaMix     = mixing table info from *MIXFIL
! iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
! iIpmix     = which of the mixed paths are being considered
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iaGases(kMaxGas)
    REAL :: raaAmt(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
! output vars
    REAL :: rCO2mult

! local vars

! these are the individual reference profiles, at kMaxLayer layers
! these are what is in kOrigRefPath
    REAL :: raR100Amt0(kMaxLayer),raR100Temp0(kMaxLayer)
    REAL :: raR100PartPress0(kMaxLayer),raR100Press0(kMaxLayer)
! these are the individual reference profiles, at kMaxLayer layers
! these are what is in kCO2ppmvFile
    REAL :: raR100Amt1(kMaxLayer),raR100Temp1(kMaxLayer)
    REAL :: raR100PartPress1(kMaxLayer),raR100Press1(kMaxLayer)

    CHARACTER(160) :: caFName
    REAL :: rMeanT,rMeanA,rMeanP,rMeanPP,rDirectPPMV0,rDirectPPMV1
    REAL :: rCO2ppmv
    INTEGER :: iI,iJ,iGasID,iError,iIPMIX,iaIndex(kMaxLayer)

! find the weight
    iIPMIX = 1
    iGasID = 2

    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Checking the CO2 ppmv in the database profiles'

    IF (kCO2_UMBCorHARTMAN == -1) THEN
      write(kStdWarn,*) 'Using HARTMANN/NIRO CO2 linemixing, so NO chi fcns'
    ELSEIF (kCO2_UMBCorHARTMAN == +1) THEN
      write(kStdWarn,*) 'Using Strow/Tobin/Machado CO2 linemixing, need chi fcns'
    ELSE
      write(kStdErr,*) 'kCO2_UMBCorHARTMAN needs to be +/- 1'
      CALL DoStop
    END IF
    IF ((strfind(kCO2Path,'lblrtm') == 1) .OR. (strfind(kCO2Path,'LBLRTM') == 1)) THEN
      IF (kCO2_UMBCorHARTMAN == +1) THEN
        write(kStdErr,*) 'kCO2Path has LBLRTM/lblrtm in it, but kCO2_UMBCorHARTMAN = -1'
        CALL DoStop
      END IF
    END IF
          
    IF (iaGases(iGasID) /= iGasID) THEN
      write(kStdErr,*) 'assumed iaGases(2) =  2, guess not'
      CALL DoStop
    END IF

!!!input pressures are in atm; lowest layers (ignore) have P = 1000.0
!!!                            layer before STtart of Atm have P = 0.0
!!!                            rest of layers have meaningful P
    DO iI = 1,kProfLayer
      IF ((raaPress(iI,iGasID) > 0.0) .AND. (raaPress(iI,iGasID) < 1.5)) THEN
        GOTO 10
      END IF
    END DO
 10 CONTINUE

    rCO2ppmv = 0.0
    DO iJ = iI,kProfLayer
      rMeanA   = raaPartPress(iJ,iGasID)/raaPress(iJ,iGasID) *1e6
      rCO2ppmv = max(rMeanA,rCO2ppmv)
    END DO
    write(kStdWarn,*) 'max rCO2ppmv from input profile = ',rCO2ppmv
    write(kStdWarn,*) 'kCARTA compiled for database CO2ppmv = ',kCO2ppmv

    IF (abs(rCO2ppmv-kCO2ppmv) >= 10) THEN
      write(kStdErr,*) 'input profile gasamts have max(CO2 ppmv) = ',rCO2ppmv
      write(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
      write(kStdErr,*) '  '
      write(kStdErr,*) '  If running NLTE SLOW suggest make a new kCARTA '
      write(kStdErr,*) '  compressed database with co2ppmv = ',rCO2ppmv
      write(kStdErr,*) '  '
      write(kStdErr,*) '  If running NLTE FAST code should be able to cope'
    END IF
          
    rCO2Mult = raaMix(iIpmix,iGasID)

!! open the true RefProf at kCO2ppmv
    caFName = kCO2ppmvFile
    write(kStdWarn,*) 'Reading CO2 Reference profile from ',caFName
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0,raR100Temp0,raR100Press0,raR100PartPress0,iError)

!! open the supposed RefProf
    CALL FindReferenceName(caFName,iGasID,-1)
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt1,raR100Temp1,raR100Press1,raR100PartPress1,iError)

! -------------------------
    DO iI = 1,kMaxLayer
      iaIndex(iI) = iI
    END DO

    rMeanA = 0.0
    rMeanA = sum(abs(raR100Amt1(iaIndex)/raR100Amt0(iaIndex)))

    rMeanPP = 0.0
    rMeanPP = sum(abs(raR100PartPress1(iaIndex)/raR100PartPress0(iaIndex)))    

    rMeanP = 0.0
    rMeanP = sum(abs(raR100Press1(iaIndex)/raR100Press0(iaIndex)))    

    rMeanT = 0.0
    rMeanT = sum(abs(raR100Temp1(iaIndex)-raR100Temp0(iaIndex)))
! -------------------------

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
    write(kStdWarn,*) 'kCO2ppmv from kcartaparam.f90       = ',kCO2ppmv
    write (kStdWarn,*) '  mean(rMeanAmt Ratio)  = ',rMeanA/kMaxLayer
    write (kStdWarn,*) '  mean(rMeanP   Ratio) = ',rMeanP/kMaxLayer
    write (kStdWarn,*) '  mean(rMeanPP  Ratio) = ',rMeanPP/kMaxLayer
    write (kStdWarn,*) '  mean(rMean   deltaT) = ',rMeanT/kMaxLayer
    write(kStdWarn,*) 'rCO2ppmv from input TRUErefprof2         = ',rDirectPPMV0
    write(kStdWarn,*) 'rCO2ppmv from input refprof2             = ',rDirectPPMV1
    write(kStdWarn,*) 'rCO2Mult from raaMixTable from .nml file = ',rCO2Mult
    write(kStdWarn,*) '  '

    IF (abs(rCO2Mult-1) >= 0.01) THEN
      write(kStdErr,*) 'you have set rCO2Mult in mixtable = ',rCO2Mult
      write(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
      write(kStdErr,*) '  suggest not trying NLTE calcs with this mixratio'
      write(kStdErr,*) '  make a new kCARTA compressed database with'
      write(kStdErr,*) '  co2ppmv = ',rCO2Mult*kCO2ppmv
    END IF

 100 FORMAT(A25,A80)
 101 FORMAT(A65,F12.8)
    IF ((rMeanA/kMaxLayer <= 0.9995) .OR. (rMeanA/kMaxLayer >= 1.0005)) THEN
      write(kStdErr,100) 'v0 : kCO2ppmvFile =   ', kCO2ppmvFile
      write(kStdErr,101) 'mean rCO2ppmv (raCO2pp/raP) v0 from input TRUErefprof2         = ',rDirectPPMV0
      write(kStdErr,100) 'v1 : ref CO2 profile = ', caFName
      write(kStdErr,101) 'mean rCO2ppmv (raCO2PP/raP) v1 from input refprof2             = ',rDirectPPMV1
      write(kStdErr,*) 'oops check_co2ppmv : rMeanA is bad',rMeanA/kMaxLayer
      CALL DoStop
    END IF
    IF ((rMeanP/kMaxLayer <= 0.9995) .OR. (rMeanP/kMaxLayer >= 1.0005)) THEN
      write(kStdErr,*) 'oops check_co2ppmv : rMeanP is bad',rMeanP/kMaxLayer
      CALL DoStop
    END IF
    IF  ((rMeanPP/kMaxLayer <= 0.9995) .OR. (rMeanPP/kMaxLayer >= 1.0005)) THEN
      write(kStdErr,*) 'oops check_co2ppmv : rMeanPP is bad',rMeanPP/kMaxLayer
      CALL DoStop
    END IF
    IF ((rMeanT/kMaxLayer <= -0.0005) .OR. (rMeanT/kMaxLayer >= +0.0005)) THEN
      write(kStdErr,*) 'oops check_co2ppmv : rMean deltaT is bad',rMeanT/kMaxLayer
      CALL DoStop
    END IF

!!! now check the NLTE weak background in Upper Layers
    write (kStdWarn,*) '    checking weak backgnd upper atm ppmvs'
    caFName = kCO2ppmvFileBackUA
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0,raR100Temp0,raR100Press0,raR100PartPress0,iError)
!! can directly figure out the ppmv in the "what we hope is true" file
    rDirectPPMV0 = 0.0
    DO iI = 1,kMaxLayer/20
        rDirectPPMV0 = rDirectPPMV0 + abs(raR100PartPress0(iI)/raR100Press0(iI))
    END DO
    rDirectPPMV0 = rDirectPPMV0/(kMaxLayer/20)*1000000
    write(kStdWarn,*) 'Weak background Upper Layers PPMV = ',rDirectPPMV0
    IF (abs(rDirectPPMV0 - rDirectPPMV1) >= 0.25) THEN
      write(kStdErr,*) 'oops check_co2ppmv : Weak Background UA LTE is bad'
      CALL DoStop
    END IF

!!! now check the NLTE weak background; the profile should be the SAME as
!!! for the Standard Profile
    write (kStdWarn,*) '    checking weak backgnd usual atm ppmvs and amts'
    caFName = kCO2ppmvFileBack
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)
!! can directly figure out the ppmv in the "what we hope is true" file
    rDirectPPMV0 = 0.0
    rMeanA = 0.0
    DO iI = 1,kMaxLayer
      rDirectPPMV0 = rDirectPPMV0+abs(raR100PartPress0(iI)/raR100Press0(iI))
      rMeanA = rMeanA + abs(raR100Amt1(iI) - raR100Amt0(iI))
    END DO
    rDirectPPMV0 = rDirectPPMV0/(kMaxLayer)*1000000
    write(kStdWarn,*) 'Weak background Standard Layers PPMV = ',rDirectPPMV0
    IF (abs(rDirectPPMV0 - rDirectPPMV1) >= 0.1) THEN
      write(kStdErr,*) 'oops check_co2ppmv : Weak Background LTE is bad'
      CALL DoStop
    END IF
    IF (rMeanA >= 0.1) THEN
      write(kStdErr,*) 'oops check_co2ppmv : Weak backgnd amts different from STD amts'
      CALL DoStop
    END IF
    write(kStdWarn,*) '  '

    RETURN
    end SUBROUTINE check_co2ppmv

!************************************************************************
! this suroutine sets up the current print options
    SUBROUTINE SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,iType, &
    iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
    iNumGases,iNpmix,iaNumLayer,iAtm)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iNumGases    !number of gases
    INTEGER :: iNpMix       !number of mix paths
    INTEGER :: iAtm         !surrent atmosphere
    INTEGER :: iaNumLayer(kMaxAtm)   !number of layers in atmospheres

    INTEGER :: iOutNum      !which printing set this is (1..kMaxPrint)
    INTEGER :: iPrinter     !what type (1,2,3)
    INTEGER :: iAtmPr       !if gas spectra, which gas; if MP, which MP,
! f radiances, which atm
    INTEGER :: iNp          !how many to output eg 100 gas spectra, or 50 MPs
    INTEGER :: iaOp(kPathsOut)  !list of paths/MPs/rads to output
    INTEGER :: iType        !-10 if dumb, 1 if paths, 2 if MPs, 3 if rads
! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output
    INTEGER :: iaPrinter(kMaxPrint),iaNp(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint)

    INTEGER :: iDummy
         
    iPrinter = iaPrinter(iOutNum)
    iAtmPr = iaGPMPAtm(iOutNum)
    iNp = iaNp(iOutNum)

    IF ((iNp < 0) .AND. (iType == 1)) THEN
    ! utput all paths for gas
        iNp = kProfLayer*iNumGases
    END IF

    IF ((iNp < 0) .AND. (iType == 2)) THEN
    ! utput all MPs
        iNp = iNpmix
    END IF

    IF ((iNp < 0) .AND. (iType == 3)) THEN
    ! utput all radiances for atmosphere
        iNp = iaNumlayer(iAtm)
    END IF

    DO iDummy = 1,iNp
        iaOp(iDummy) = iaaOp(iOutNum,iDummy)
    END DO

    RETURN
    end SUBROUTINE SetUpCurrentPrint

!************************************************************************
! this subroutine does some more initializations
    SUBROUTINE SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iMixFileLines,iVertTempSet,iNpmix
    REAL :: raaMix(kMixFilRows,kGasStore)       !mixing table
     
    INTEGER :: iFileIDLO,iFileIDhi

! these are needd to be set if kRTP = -10 (user supplied levels profile) or -5,-6 (TAPE5/TAPE6 from LBLRTM)
!   1,2,3 = SurfPress/SurfTemp/SurfHgt
!   4,5   = ViewHgt and ViewAngle/direction
!   6     = Based on (4=ViewHgt) the code figures out highest pressure at which rad is to be output
    DO iFileIDLO = 1,10
        raRTP_TxtInput(iFileIDLO) = -9999
    END DO

! these are for seeing how cloud params are overwritten
! raaRTPCloudParams0(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12 iCldoTp iCLdBot  from rtpfile
! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12 iCldTop iCldBot  after kcarta resets
! this gets set in rtp_interface.f
    raaRTPCloudParams0(:,1:9) = -9999.0
    raaRTPCloudParamsF(:,1:9) = -9999.0
! this gets set correctly in SetMieTable_RTSPEC
    iaCloudTypeProfile = -1          
 
! this is really for Mie scattering and VIS/UV ocean reflectance
    kSolAzi = 0.0
    kSatAzi = 0.0
    kWindSpeed = 0.0
    kLatitude = 0.0
    kMonth = 1.0
! this is to stop the code flux calcs at LBLRTM toa, default = -1 so keep doing calcs till 0.005 mb
! else if kLBLRTM_toa > 0 then the flux code
!   finds the highest layer whose pressure is greater than this
!   sets all ODS to 0 above this layer so essentially
!     rU = rUp0 exp(-k) + B(T) (1-exp(-k)) --> r0 for upwelling rad
!     rD = 0 since the downwelling rad is initialized with ttorad(f,2.7) ~ 0
    kLBLRTM_toa = -1.0
          
! assume no *mixfil section
    iMixFileLines = -1
     
! the vertical temperature profile has not been set yet
    iVertTempSet = -1
     
! initialize the mixing table to weights of 0.0
    iNpMix = 1
    raaMix(1:kMixFilRows,1:kGasStore) = 0.0

! initialize kaaNumVectors
    kaaNumVectors(1:kMaxGas,1:kMaxLayer) = 0

    RETURN
    end SUBROUTINE SomeMoreInits

!************************************************************************
! this subroutine inits file unit numbers
    SUBROUTINE InitializeFileUnits
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! set up file unit numbers so at run time they default to STDIN,STDOUT,STDOUT
    kStdDriver = 5
    kStdkCarta = 6
    kStdJacob  = 6
     
! set up common block parameters indicating all units closed
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
    end SUBROUTINE InitializeFileUnits

!************************************************************************
! this subroutine opens the driver namelist file and figures out the name
!  of the error/warning log file
    SUBROUTINE ErrorLogName(caDriverName)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! this is the driver file name from the command line arguments
    CHARACTER(160) :: caDriverName

    CHARACTER(30) :: namecomment

! this is for OUTPUT
! caLogFile     = name of success/warning log file 'warning.msg'
! caComment     = comment the user writes
! iOutTypes     = number of printing options specified
! iaPrinter     = for each option, which output type specified
! iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output
! iaNp          = for each option, how many paths/MPs/layers to be output
! iaaOp         = for each option, list of paths/MP/layers to be output
! raaOp         = for option 3, list fract of layers used for radiance output
! raaUserPress  = for option 3, list of pressures for output radiances
! iNatm2        = number of atmospheres that *OUTPUT thinks there is
    INTEGER :: iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
    INTEGER :: iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
    INTEGER :: iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
    INTEGER :: iIOUN,iErr
    CHARACTER(160) :: caComment,caComment1
    CHARACTER(160) :: caLogFile,caLogFile1
    REAL :: raaOp(kMaxPrint,kProfLayer),raaOp1(kMaxPrint,kProfLayer)

    NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter, &
    iaGPMPAtm,iaNp,iaaOp,raaOp

 1070 FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80)
    iIOun=kStdDriver
    IF (iIOUN /= 5) THEN
      OPEN(UNIT = iIOun,FILE = caDriverName,STATUS='OLD',IOSTAT=iErr)
      IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine ErrorLogName, error reading'
        write (kStdErr,*) 'namelist file to find name of logfile ... '
        WRITE(kStdErr,1070) iErr, caDriverName
        CALL DoSTOP
      ENDIF
    END IF
    kStdDriverOpen = 1

!      write(kStdWarn,*) 'grepping input nml file for caLogFile name'
!      print *,'translating x...'

    namecomment = '******* OUTPUT section *******'
    caLogFile = 'warning.msg'     !this is the default name
    read (iIOUN,nml = nm_output)
    caLogFile1  =  caLogFile

    kWarnFile  =  caLogFile

    close (iIOUN)
    kStdDriverOpen = -1

    RETURN
    end SUBROUTINE ErrorLogName

!************************************************************************
! this subroutine summarizs the output options
    SUBROUTINE SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
    raaOp,raaUserPress)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output
    INTEGER :: iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint)
    REAL :: raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer)

    INTEGER :: iDummy,iOutnum

    write(kStdWarn,*) '# of printing options selected = ',iOuttypes
    write(kStdWarn,*) '     index     option type      atm#    numpaths'
    write(kStdWarn,*) '------------------------------------------------'
    DO iDummy = 1,iOuttypes
      write(kStdWarn,30) iDummy,iaPrinter(iDummy),iaGPMPAtm(iDummy),iaNp(iDummy)
      write(kStdWarn,*) 'paths to be printed : (if numpaths=-1,print all)'
      write(kStdWarn,*)(iaaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy))
      IF (iaPrinter(iDummy) == 3) THEN
        write(kStdWarn,*)(raaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy))
        write(kStdWarn,*)(raaUserPress(iDummy,iOutNum),iOutNum=1,iaNp(iDummy))
      END IF
      write(kStdWarn,*) '    '
    END DO
     
    30 FORMAT('     ',4(I3,'          '))

    RETURN
    end SUBROUTINE SummaryOutputs

!************************************************************************
! this subroutine checks the MixPath Vertical Temps
    SUBROUTINE CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp, &
    iNpmix,iCO2,iaGases)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
     
    INTEGER :: iCo2              !which gas used to mimic CO2 temps
    INTEGER :: iNumGases         !how many gases
    INTEGER :: iNpMix            !number of mixed paths
    REAL :: raaTemp(kProfLayer,kGasStore)       !profile temp
    REAL :: raaMix(kMixFilRows,kGasStore)       !mixing table
    REAL :: raMixVertTemp(kMixFilRows)          !temperatures of MP layers
    INTEGER :: iaGases(kMaxGas)               !gasIDs stored in order

    INTEGER :: iDummy,iFileIDLo

    iCO2 = -1

    IF (iNpmix > 0) THEN
      ! search for the CO2 gas === gasID 2
      ! since the gases have to be entered in ascending order, either gas1 or gas2
      ! is CO2
      iCO2 = -1
      IF (iaGases(1) == 2) THEN
        iCO2 = 1
        write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!'
      ELSEIF (iaGases(2) == 2) THEN
        iCO2 = 2
        write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!'
      ELSE !!!for some strange reason, no CO2 present
        iCO2 = 1
        write(kStdWarn,*) 'Temperature of Gas number 1 will mimic CO2!!'
      END IF
         
      CALL GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp,iNpmix,iCO2)
         
      write(kStdWarn,*) 'Checking Mixed Path Temp'
      iFileIDLo = 0
      DO iDummy = 1,iNpmix
        IF (raMixVertTemp(iDummy) < 0.0) THEN
          write(kStdWarn,*) 'Negative MP Temp in Mixed Path',iDummy
          iFileIDLo = iFileIDLo+1
        END IF
      END DO
      IF (iFileIDLo > 0) THEN
        write(kStdWarn,*) 'Warning! negative MP temperatures found!'
        write(kStdErr,*) 'Warning! negative MP temperatures found!'
        CALL DoSTOP
      END IF
    END IF

    1111 FORMAT(A1)

    RETURN
    end SUBROUTINE CheckMixedPathTemps

!************************************************************************
! this subroutine does the command line stuff
    SUBROUTINE DoCommandLine(iMicrosoft,caDriverName,caOutName, &
    caJacobFile,iOutFileName)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! caDriverName is the name of the driver file to be processed
    CHARACTER(160) :: caDriverName
! caOutName is the name of the unformatted output file name
! integer iOutFileName tells whether or not there is a driver name, or
! dump output to Unit 6
    INTEGER :: iOutFileName
    CHARACTER(160) :: caOutName
! caJacobFile is the name of the unformatted output file name for Jacobians
    CHARACTER(160) :: caJacobFile
! this tells if we have MS Product ie no command line stuff!
    INTEGER :: iMicroSoft
! this is the number of args
    INTEGER :: iargc

    INTEGER :: iDummy,iError

    IF (iMicroSoft > 0) THEN                 
      ! do command line options .. do this!
      print *,'Enter (1) if standard kcarta computation '
      print *,'Enter (2) if kcarta + jacobian computation '
      read *,iDummy
      IF ((iDummy > 2) .OR. (iDummy < 1)) THEN
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
         
      IF (iDummy == 2) THEN
        print *,'Enter output jacobian filename  (enclose in quotes) : '
        read *,caJacobFile
        kStdJacob = kStdJacobKK
      END IF

      iMicrosoft = iDummy    !tells number of output files (w/o flux, planck)
               
    ELSE
      ! use command line stuff
      iDummy = iargc()
         
      IF (iDummy > 3) THEN
        write(kStdErr,*) 'more than three arguments in command line'
        write(kStdErr,*) 'is NOT allowed'
        CALL DoSTOP
      END IF

      iOutFileName = -1         !assume no name
      DO iError = 1,iDummy
        IF (iError == 1) THEN
          CALL getarg(1,caDriverName)
          IF (caDriverName(1:1) /= '-') THEN
            kStdDriver = kStdDriverKK
          END IF
        END IF
             
        IF (iError == 2) THEN
          CALL getarg(2,caOutName)
          IF (caOutName(1:1) /= '-') THEN
            iOutFileName = 1
            kStdkCarta = kStdkCartaKK
          END IF
        END IF
             
        IF (iError == 3) THEN
          CALL getarg(3,caJacobFile)
          IF (caJacobFile(1:1) /= '-') THEN
            kStdJacob = kStdJacobKK
          END IF
        END IF
      END DO
         
      iMicrosoft = iDummy-1  !tells number of output files (w/o flux,planck)
    END IF

    RETURN
    end SUBROUTINE DoCommandLine
     
!************************************************************************
! this subroutine stores the reference gas amts/temps etc
    SUBROUTINE StoreReference(raRAmt,raRTemp,raRPress,raRPartPress, &
    raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! this is the gasnumber (not gasID!!!!!!!!!!!!!!)
    INTEGER :: iGas
    INTEGER :: iaGases(kMaxGas)
! these are the individual reference profiles
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the reference profiles stored in matrices
    REAL :: raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
    REAL :: raaRPress(kProfLayer,kGasStore)
    REAL :: raaRPartPress(kProfLayer,kGasStore)

    INTEGER :: iI

    DO iI = 1,kProfLayer
      raaRAmt(iI,iGas) = raRAmt(iI)             !amts
      raaRTemp(iI,iGas) = raRTemp(iI)           !temps
      raaRPress(iI,iGas) = raRPress(iI)         !press
      raaRPartPress(iI,iGas) = raRPartPress(iI) !part press

      IF (raRAmt(iI) < 0.0) THEN
        WRITE(kStdErr,*) 'Error in Ref Profile for Gas',iaGases(iGas)
        WRITE(kStdErr,*) 'Layer ',iI,' has negative amount',raRAmt(iI)
        CALL DoStop
      END IF
         
      IF (raRTemp(iI) < 0.0) THEN
        WRITE(kStdErr,*) 'Error in Ref Profile for Gas',iaGases(iGas)
        WRITE(kStdErr,*) 'Layer ',iI,' has negative tempr',raRTemp(iI)
        CALL DoStop
      END IF
         
      IF (raRPress(iI) < 0.0) THEN
        WRITE(kStdErr,*) 'Error in Ref Profile for Gas',iaGases(iGas)
        WRITE(kStdErr,*) 'Layer ',iI,' has negative pressure',raRPress(iI)
        CALL DoStop
      END IF
         
      IF (raRPartPress(iI) < 0.0) THEN
        WRITE(kStdErr,*) 'Error in Ref Profile for Gas',iaGases(iGas)
        WRITE(kStdErr,*) 'Layer ',iI,' has negative part press',raRPartPress(iI)
        CALL DoStop
      END IF
    END DO

    RETURN
    end SUBROUTINE StoreReference

!************************************************************************
! this subroutine sets the reference gas amts/temps etc
    SUBROUTINE SetReference(raRAmt,raRTemp,raRPress,raRPartPress, &
    raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! this is the gasnumber (not gasID!!!!!!!!!!!!!!)
    INTEGER :: iGas
! these are the individual reference profiles
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the reference profiles stored in matrices
    REAL :: raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
    REAL :: raaRPress(kProfLayer,kGasStore)
    REAL :: raaRPartPress(kProfLayer,kGasStore)

    INTEGER :: iInt

    raRAmt(1:kProfLayer) = raaRAmt(1:kProfLayer,iGas)
    raRTemp(1:kProfLayer) = raaRTemp(1:kProfLayer,iGas)
    raRPress(1:kProfLayer) = raaRPress(1:kProfLayer,iGas)
    raRPartPress(1:kProfLayer) = raaRPartPress(1:kProfLayer,iGas)

    RETURN
    end SUBROUTINE SetReference

!************************************************************************
! this subroutine close units
    SUBROUTINE TheEnd(iaGases,iNumGases,iaList,raFiles)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
          
    INTEGER :: iNumGases,iaGases(kMaxGas),iaList(kNumkCompT)
    REAL :: raFiles(kNumkCompT)
    INTEGER :: iI,iJ,iG,iaSum(kMaxGas),iaCount(kMaxGas),iSum,iaJunk(kMaxGas)
    CHARACTER*50 FMTX
    CHARACTER*3  str3

    write(kStdWarn,*) '*****************************************'
    write(kStdWarn,*) 'kComp stats ..............'
    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Freq chunks processed'
    write(kStdWarn,*) (raFiles(iaList(iJ)),iJ=1,kOuterLoop)

    write(kStdWarn,*) 'Stats of compressed vecs per gas per chunk'
    iaSum = 0
    iaCount = 0
          
    DO iI = 1,iNumGases
      iG = iaGases(iI)
      write(kStdWarn,*) 'num compressed vector stats for gas #, gasID ',iI,iG
      write(kStdWarn,*) (kaaNumVectors(iG,iJ),iJ=1,kOuterLoop)
      write(kStdWarn,*) ' '
      DO iJ=1,kOuterLoop
        IF (kaaNumVectors(iG,iJ) > 0) THEN
          iaSum(iG) = iaSum(iG) + kaaNumVectors(iG,iJ)
          iaCount(iG) = iaCount(iG) + 1
        END IF
      END DO
    END DO

    write(kStdWarn,*) 'gas#   GASID   NumVecs   NumChunks  AvgVecs'
    DO iI = 1,iNumGases
      iG = iaGases(iI)
      IF (iaCount(iG) > 0) THEN
        write(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG),1.0*iaSum(iG)/iaCount(iG)
      ELSE
        write(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG),0.0
      END IF
    END DO
  15 FORMAT(2(I5,' '),'   ',I5,'      ',I5,'     ',F8.4)

    write(kStdWarn,*) ' '
    write(kStdWarn,*) '  Chunk   NumGases'
    DO iJ = 1,kOuterLoop
      iSum = 0
      DO iI = 1,iNumGases
        iG = iaGases(iI)
        IF (kaaNumVectors(iG,iJ) > 0) iSum = iSum + 1
      END DO
      write(kStdWarn,16) raFiles(iaList(iJ)),iSum
    END DO
    16 FORMAT(F9.2,'    ',I3)

    write(kStdWarn,*) ' '
    DO iJ = 1,kOuterLoop
      iSum = 0
      DO iI = 1,iNumGases
        iG = iaGases(iI)
        IF (kaaNumVectors(iG,iJ) > 0) THEN
          iSum = iSum + 1
          iaJunk(iSum) = iG
        END IF
      END DO
      write (str3, '(i3)') iSum
      !FMTX = ''(' // str3 // '(I3,1X)'')'
      write(kStdWarn,'(A,F12.5,A,I3,A)') '  Chunk = ',raFiles(iaList(iJ)), ' numgases = ',iSum,' gasIDs are .... '
      !write(kStdWarn,FMTX) (iaJunk(iI),iI=1,iSum)
      write(kStdWarn,*) (iaJunk(iI),iI=1,iSum)
    END DO
    17 FORMAT(I3,' ')

    IF (kStdkCarta /= 6) THEN
      write(kStdWarn,*)'closing binary output file'
      CLOSE(UNIT = kStdkCarta)      !close file where kCARTA binary goes to
      kStdkCartaOpen  =  -1
    END IF
      
    IF ((kJacobian > 0) .AND. (kStdJacob /= 6)) THEN
      write(kStdWarn,*)'closing jacobian binary file'
      CLOSE(UNIT = kStdJacob)      !close file where Jacob binary goes to
      kStdJacobOpen  =  -1
    END IF

    IF (kJacobian > 0) THEN
      IF ((kStdJacobOpen == 1)  .AND. (kStdJacob /= 6)) THEN
        write(kStdWarn,*)'closing jacobian binary file'
        CLOSE(UNIT = kStdJacob)       !close file where Jacob binary goes to
      END IF
      IF (kStdJacob2Open == 1) THEN
        write(kStdWarn,*)'closing jacobian2 (column) binary file'
        CLOSE(UNIT = kStdJacob2)       !close file where Jacob binary goes to
      END IF
    END IF

!      IF (kJacobian .GT. 0) THEN
!        write(kStdWarn,*)'closing jacobian2 (column) binary file'
!        CLOSE(UNIT = kStdJacob2)      !close file where Jacob binary goes to
!        kStdJacob2Open  =  -1
!      END IF
     
    IF (kFlux > 0) THEN
      write(kStdWarn,*)'closing flux binary file'
      CLOSE(UNIT = kStdFlux)       !close file where flux binary goes to
      kStdFluxOpen  =  -1
    END IF

    IF (kStdPlanckOpen > 0) THEN
      write(kStdWarn,*)'closing planck binary file'
      CLOSE(UNIT = kStdPlanck)    !close file where planck binary goes to
      kStdPlanckOpen  =  -1
    END IF

    IF (kBloatPlanckOpen > 0) THEN
      write(kStdWarn,*)'closing bloated planck binary file'
      CLOSE(UNIT = kBloatNLTEPlanck)    !close file
      kBloatPlanckOpen  =  -1
    END IF

    IF (kBloatOutOpen > 0) THEN
      write(kStdWarn,*)'closing bloated binary file'
      CLOSE(UNIT = kBloatNLTEOut)    !close file
      kBloatOutOpen  =  -1
    END IF

    IF (kStdPlanckUAOpen == 1) THEN
      write(kStdWarn,*)'closing UA planck binary file'
      CLOSE(UNIT = kStdPlanckUAOpen)      !close file
      kStdPlanckUAOpen = -1
    END IF

    IF (kNLTEOutUAOpen == 1) THEN
      write(kStdWarn,*)'closing UA binary file'
      CLOSE(UNIT = kNLTEOutUAOpen)      !close file
      kNLTEOutUAOpen = -1
    END IF

    IF (kBloatPlanckUAOpen == 1) THEN
      write(kStdWarn,*)'closing bloat UA planck binary file'
      CLOSE(UNIT = kBloatPlanckUAOpen)      !close file
      kBloatPlanckUAOpen = -1
    END IF

    IF (kBloatNLTEOutUAOpen == 1) THEN
      write(kStdWarn,*)'closing bloat UA binary file'
      CLOSE(UNIT = kBloatNLTEOutUAOpen)      !close file
      kBloatNLTEOutUAOpen = -1
    END IF

    write(kStdWarn,*) 'closed files ... !!!!!!!!!!!'
      
    RETURN
    end SUBROUTINE TheEnd

!***********************************************************************
! this subroutine initializes all the rows of the
! (REAL) array of absorption coeffs
    SUBROUTINE InitializeReal(raaAb)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: raaAb(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFreq

    raaAb = 0.0

    RETURN
    end SUBROUTINE InitializeReal

!************************************************************************
! this subroutine initializes all the rows of the
! (REAL) array of mixed paths
    SUBROUTINE InitializeRealMP(raaAb,iNpmix)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: raaAb(kMaxPts,kMixFilRows)
    INTEGER :: iNpMix

    INTEGER :: iLay,iFreq

    DO iLay = 1,iNpMix      !note : initialize only wot is necessary
      raaAb(:,iLay) = 0.0
    END DO

    RETURN
    end SUBROUTINE InitializeRealMP

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs
    SUBROUTINE initialize(daaAb)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    DOUBLE PRECISION :: daaAb(kMaxPts,kProfLayer)

    daaAb = 0.0d0

    RETURN
    end SUBROUTINE initialize

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs
! pretty much the same as above routine
    SUBROUTINE initializeJAC(daaAb)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    DOUBLE PRECISION :: daaAb(kMaxPtsJac,kProfLayerJac)

     daaAb = 0.0d0

    RETURN
    end SUBROUTINE initializeJAC

!************************************************************************
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa
! and saves it in an overall AbsMatrix raaa
    SUBROUTINE DoSet(daaGasAbCoeff,raaaGasAbCoeff,iCount,iDoAdd)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iCount     = which of the iNumGases are being processed
! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPtsJac,kProfLayerJac)
    REAL :: raaaGasAbCoeff(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    INTEGER :: iCount,iDoAdd

! local variables
    INTEGER :: iLay,iFr

    IF (iDoAdd > 0) THEN
      raaaGasAbCoeff(iCount,:,:) = daaGasAbCoeff
    ELSEIF (iDoAdd <= 0) THEN
      raaaGasAbCoeff(iCount,:,:) = 0.0
    END IF

    RETURN
    end SUBROUTINE DoSet
!************************************************************************c
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa
    SUBROUTINE DoDtoR(daaGasAbCoeff,raaGasAbCoeff)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)

    raaGasAbCoeff = real(daaGasAbCoeff)

    RETURN
    end SUBROUTINE DoDtoR

!************************************************************************
! this subroutine converts the layer iL to 0
    SUBROUTINE ZeroLayer(raaGasAbCoeff,iL)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaGasAbs = 2d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)
    INTEGER :: iL

    raaGasAbCoeff(:,iL) = 0.0

    RETURN
    end SUBROUTINE ZeroLayer

!************************************************************************
! this subroutine checks parameters in kcartaparam.f90 ... abort if they
! do not make sense .. the parameters are set in kcartaparam.f90
    SUBROUTINE CheckKCARTAParameters

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'
    include '../INCLUDE/TempF90/airslevelsparam.f90'
    include '../INCLUDE/TempF90/airsheightsparam.f90'
    include '../INCLUDE/TempF90/airslevelheightsparam.f90'

    write(kStdWarn,'(A,I3)') 'checking ../INCLUDE/TempF90/airsHeightsLevels*aparam.f90 are for kPlanet = ',kPlanet
    write(kStdWarn,*) '  '

    IF (kPlanet /= iXPlanet0) THEN
      write(kStdErr,'(A,I3,A,I3)') 'kCARTA compiled for kPlanet = ',kPlanet,'but reading in airs*param for iXPlanet0',iXPlanet0
      Call DoStop
    END IF

    IF (kPlanet /= iXPlanet1) THEN
      write(kStdErr,'(A,I3,A,I3)') 'kCARTA compiled for kPlanet = ',kPlanet,'but reading in airs*param for iXPlanet1',iXPlanet1
      Call DoStop
    END IF
    IF (kPlanet /= iXPlanet2) THEN
      write(kStdErr,'(A,I3,A,I3)') 'kCARTA compiled for kPlanet = ',kPlanet,'but reading in airs*param for iXPlanet2',iXPlanet2
      Call DoStop
    END IF
    IF (kPlanet /= iXPlanet3) THEN
      write(kStdErr,'(A,I3,A,I3)') 'kCARTA compiled for kPlanet = ',kPlanet,'but reading in airs*param for iXPlanet3',iXPlanet3
      Call DoStop
    END IF

    write(kStdWarn,*) 'checking parameters in kcartaparam.f90'
    write(kStdWarn,*) '  '
     
    CALL Check_kaNum_ka100layerCloudType

! do a quick check of the important parameters set by the user
    IF (kMixFilRows < kProfLayer) THEN
      write(kStdErr,*) 'In kcartaparam.f90, need '
      write(kStdErr,*) 'kMixFilRows >= kProfLayer(=',kProfLayer,')'
      write(kStdErr,*) 'please reset and retry'
      CALL DoSTOP
    END IF

    IF (abs(kXsecFormat) /= 1) THEN
      write(kStdErr,*) 'kXsecFormat in kcartaparam.f90 must be = +/-1'
      write(kStdErr,*) 'please reset and retry'
      CALL DoSTOP
    END IF

    write(kStdWarn,*) 'Max #of atmospheres from *RADFIL = ',kMaxAtm
    write(kStdWarn,*) 'Max #of gases from *GAS/XSCFIL = ',kGasStore
    write(kStdWarn,*) 'Max #of mixed paths *MIXFIL = ',kMixFilRows
    write(kStdWarn,*) '  '

    RETURN
    end SUBROUTINE CheckKCARTAParameters

!************************************************************************
! this subroutine sets the mixed path effective tempertures, weighted
! according to the mixing table
    SUBROUTINE GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, &
    iNpmix,iCO2)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNumGases   = number of gases read in from *GASFIL + *XSCFIL
! iNpMix      = number of mixed paths read in from *MIXFIL
! raaTemp     = temperature profile of ONE of the gases (assume all equal)
! raMixVTemp  = computed vertical temp profile for mixing table
! raaMix      = mixing table from *MIXFIL
! iCO2        = which of the gases is CO2 .. if -1, none are CO2
    INTEGER :: iNumGases,iNpmix,iCO2
    REAL :: raMixVertTemp(kMixFilRows)
    REAL :: raaTemp(kProfLayer,kGasStore)
    REAL :: raaMix(kMixFilRows,kGasStore)

! local variables
    INTEGER :: iI,iJ,iL,iBad
    REAL :: rT,rW

    iBad = 0
    raMixVertTemp = 0.0

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
! default = -1

    rW = 0.0
    rW = sum(raaMix(50,1:iNumGases))

!    IF (rW <= 1.0e-5) THEN
    IF (rW < 1.0e-8) THEN
      write(kStdErr,*) 'arbitrarily tested LAYER 50, mixed path weights = 0 .. perhaps you have nm_weight wrong???'
      write(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
      CALL DoStop
    END IF
            
    IF ((kGasTemp == 1) .AND. (iCO2 > 0)) THEN
      ! user wants the CO2 profile to be the temperature profile
      DO iI = 1,iNpmix
        iL = MP2Lay(iI)
        raMixVertTemp(iI) = raaTemp(iL,iCO2)
      END DO
    ELSE
      ! calculate the weights
      DO iI = 1,iNpmix
        rT = 0.0
        rW = 0.0
        iL = MP2Lay(iI)
        DO iJ = 1,iNumGases
          rT = rT+raaTemp(iL,iJ)*raaMix(iI,iJ)
          rW = rW+raaMix(iI,iJ)
        END DO
        rT = rT/rW
!       IF (rW <= 1.0e-5) THEN
        IF (rW < 1.0e-8) THEN
          iBad = iBad + 1
          write(kStdErr,*) 'hmm, mixed path weight = 0 in GetMixVertTemp for layer ',iI
        END IF
        raMixVertTemp(iI) = rT
        ! if the weights are set so that mixed path defines unique layers
        ! these temperatures should now be equal
      END DO
      IF (iBad > 0) THEN
        write(kStdErr,*) 'had total ',iBad,' mixed path weights = 0 .. perhaps you have nm_weight wrong???'
        write(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
        CALL DoStop
      END IF
    END IF

    RETURN
    end SUBROUTINE GetMixVertTemp

!************************************************************************
! this subroutine computes avg pressure of the layers
    SUBROUTINE FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: raPressLevels(kProfLayer+1),pProf(kProfLayer)
    INTEGER :: iProfileLayers

    INTEGER :: iI,iJ,iDiff

     pProf = 0.0

    iDiff = kProfLayer - iProfileLayers
    DO iI = 1,iProfileLayers
      iJ = iI + iDiff
      pProf(iJ) = raPressLevels(iJ+1)-raPressLevels(iJ)
      pProf(iJ) = pProf(iJ)/log(raPressLevels(iJ+1)/raPressLevels(iJ))
    END DO

    RETURN
    end SUBROUTINE FindAvgLayerPressure

!************************************************************************
! this subroutine finds the partial pressure of the layer
! especially useful for GasID = 1
    SUBROUTINE PPThruLayers( &
    iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press, &
    raDataBaseThickness,raSortPressLevels,raSortPressHeights, &
    raPressLevels,raRTemp,r)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
     
! input variables
    INTEGER :: iLowerOrUpper   !!upper or lower atm
    INTEGER :: iI              !!layer of current profile we are looking at
    INTEGER :: iGasID
! these are the individual reference profiles, at kMaxLayer layers
    REAL :: raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
    REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the kCARTA database pressure levels and heights
    REAL :: raDataBaseThickness(kMaxLayer)  !!AIRS layer thickness (unsorted) m
    REAL :: raSortPressLevels(kMaxLayer+1)  !!sorted with increasing p        mb
    REAL :: raSortPressHeights(kMaxLayer+1) !!sorted with increasing p        m
! these are current profile pressure levels
    REAL :: raPressLevels(kProfLayer+1)
! this is the interpolated reference temperature
    REAL :: raRTemp(kProfLayer)
! output variable
    REAL :: r    !!the partial pressure

! local vars
    INTEGER :: iL,iLp1, iJ,iJp1,iMMM,iCase
    REAL :: r0,rH1,rH2,rP,rPp1,rHH,rPP,rDeltaH,r1,rNum,rDenom,rNum1,rDenom1,rXX

    REAL :: DatabaseHeight(kMaxLayer)
    REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: DATABASELEV(kMaxLayer+1)

    CALL databasestuff(iLowerOrUpper,DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

    iMMM = kMaxLayer + 1
    r0 = r               !!! partial pressure from naive interpolation

    iL   = iI            !!!lower pressure level of layer iI
    iLp1 = iI + 1        !!!upper pressure level of layer iI
    rP   = raPressLevels(iL)
    rPp1 = raPressLevels(iL+1)

    iJ = 1
!!!find databaselev(iJ) <= rP by looping up from GND
 10 CONTINUE
    IF (DATABASELEV(iJ) > rP) THEN
      iJ = iJ + 1
      IF (iJ < kMaxLayer+1) THEN
        GOTO 10
      END IF
    END IF
          
! check to see if we are OK or if we went one too far; if so, go down one
    IF ((DATABASELEV(iJ) < rP) .AND. (iJ > 1)) iJ = iJ - 1

    iJp1 = kMaxLayer + 1
!!!find databaselev(iJp1) >= rPp1 by looping down from TOA
 20 CONTINUE
    IF (DATABASELEV(iJp1) < rPp1) THEN
      iJp1 = iJp1 - 1
      IF (iJp1 > 1) THEN
        GOTO 20
      END IF
    END IF
! check to see if we are OK or if we went one too far; if so, go up one
    IF (DATABASELEV(iJp1) > rPp1) iJp1 = iJp1 + 1

    IF ((iJp1 - iJ) <= 0) THEN
      ! something wrong; we need lower DATABASE level < upper DATABASE level
      write(kStdErr,*) 'doing water amount ref partial pressure : '
      write(kStdErr,*) 'Layer iI : found iJ = iJp1 ',iI,iJ
      CALL DoStop
    END IF

    IF (iJ <= 0) THEN
      write(kStdErr,*) '>>layer number ',iI,'gas ID ',iGasID
      Call DoStopMesg(' >>oops iJ <= 0 in PPThruLayers, so will have problems in next line')
    END IF
        
!!! case 1 : current layer is INSIDE one database layer
    IF ((iJp1 - iJ) == 1) THEN
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
    IF ((iJp1 - iJ) == 2) THEN
      iCase = 2

      !!! initialise with contribution from lower layer
      rDeltaH = (DATABASELEV(iJ+1)-rP)/(DATABASELEV(iJ+1)-DATABASELEV(iJ))
      rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
      rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

      !!! add on contribution from upper layer
      rDeltaH = (rPp1-DATABASELEV(iJ+1))/(DATABASELEV(iJ+2)-DATABASELEV(iJ+1))
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
    IF ((iJp1 - iJ) > 2) THEN
      iCase = 3

      !!! initialise with contribution from lowest layer
      rDeltaH = (DATABASELEV(iJ+1)-rP)/(DATABASELEV(iJ+1)-DATABASELEV(iJ))
      rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
      rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

      !!! add on contributions from intermediate layers
      DO iL = iJ+1,iJp1-2
        rDeltaH = 1.0
        rNum = rNum + raR100Amt(iL)*rDeltaH                   !!kmol/cm2
        rDenom = rDenom + raDataBaseThickness(iL)
      END DO

      !!! add on contribution from highest layer
      rDeltaH = (rPp1-DATABASELEV(iJp1-1))/(DATABASELEV(iJp1)-DATABASELEV(iJp1-1))
      rNum1   = raR100Amt(iJp1-1)*rDeltaH                         !! kmol/cm2
      rDenom1 = raDataBaseThickness(iJp1-1)*rDeltaH               !! frac

      rNum   = rNum + rNum1
      rDenom = rDenom + rDenom1
      r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
      r = r*1e6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
      r = r/kAtm2mb                                          !!atm
    END IF

    RETURN
    end SUBROUTINE PPThruLayers

!************************************************************************
! this is called by kcartamain and kcartabasic, to set up the profiles
    SUBROUTINE Set_Ref_Current_Profs( &
    iJax,iJax2,iGasJac,iaNumLayer,iaaRadLayer,rDerivTemp,rDerivAmt, &
    iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress, &
    raaAmt,raaTemp,raaPress,raaPartPress, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    raNumberDensity,pProfNLTE,raMixVertTemp)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input params
    INTEGER :: iJax,iJax2             !! when testing jacobians
    INTEGER :: iGasJac                !! which gas to perturb (-1 for temperature)

    REAL :: rDerivTemp,rDerivAmt      !! when testing jacobians
    INTEGER :: iGas,iaGases(kMaxGas)  !! gasID stored in order they were read in
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
! these are the reference profiles stored in matrices
    REAL :: raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
    REAL :: raaRPress(kProfLayer,kGasStore)
    REAL :: raaRPartPress(kProfLayer,kGasStore)
! these are the user specified layer profiles stored in matrices
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    REAL :: raMixVertTemp(kMixFilRows)

! output params
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
! this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
! when multiplied by the layer height, gives units of /cm^2
! so when multiplied by Rayleigh scattering cross section, this gives units
! of optical depth (no units)
    REAL :: raNumberDensity(kProfLayer)
    REAL :: pProfNLTE(kProfLayer),rDummy

! local vars
    INTEGER :: iInt,iAtm,iX,iX2,dq_Const_between_layers,iStartKeepTrack
    REAL :: rDQ_Track_Cumulative,rDQ_Track,rAdjust_rDerivAmt,rJunk,rUseDQ

    iAtm = 1

! get the reference profile for the current gas if GAS ID <= kGasXsecHi
    IF ((iaGases(iGas) <= kGasXsecHi) .OR. (iaGases(iGas) == kNewGasHi+1)) THEN
      CALL SetReference(raRAmt,raRTemp,raRPress,raRPartPress, &
        raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
    END IF

! ************************************************************************

    !  also see ../MATLAB/compare_sumlayer_jac_dRdq_vs_qdRdq.m
    !  also see ../MATLAB/compare_sumlayer_jac_dRdq_vs_qdRdq.m
    !  also see ../MATLAB/compare_sumlayer_jac_dRdq_vs_qdRdq.m
    
    rDQ_Track_Cumulative = 0.0 

    iStartKeepTrack = -1
    rAdjust_rDerivAmt = 1.0
    IF ((kJacobOutput .EQ. 0) .OR. (kJacobOutput .EQ. 1)) THEN
      !! q dR/dq or q dBT/dq
      !! want dq1 = dq2 = dq3 = .... dqm 
      !!   therefore Q1 gammaQ1 = Q2 gammQ2 = Q3 gammaQ3 = ..... Qm gammaQm
      !!   therefore gammQ2/gammaQ1 = Q1/Q2      gammQ3/gammaQ1 = Q1/Q3 .....      gammQm/gammaQ1 = Q1/Qm   
      !! see detailed notes, Book 46 of my spiral bound notebooks
      dq_Const_between_layers = +1  
      rDQ_Track = 1.0
    ELSEIF ((kJacobOutput .EQ. -1) .OR. (kJacobOutput .EQ. 2)) THEN
      !! dR/dq or dBT/dq
      dq_Const_between_layers = -1
      rDQ_Track = 1.0
    END IF

    !!!!!!!!!!!!!!!!!!!!!!!!!!
    iJax = 600        !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 
    !!!!!!!!!!!!!!!!!!!!!!!!!!

    !! overwrite iJax,iJax2,iGasJac
    iJax = 600        !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 

    iJax = 52
    iJax = 32
    iGasJac = 6

    !! overwrite iJax,iJax2,iGasJac
    iJax = 76
    iJax = 85
    iGasJac = 3

    !! overwrite iJax,iJax2,iGasJac
    iJax = 600        !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 

    iJax = 30
    iJax = 20
    iJax = 10
    iJax = 5
    iJax = 4
    iGasJac = 1

    !! overwrite iJax,iJax2,iGasJac
    iJax = 600        !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 

    iJax = 6
    iJax = 20
    iJax = 5
    iGasJac = 2

    !! *************************
    !! *************************
    !! *************************

    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations
    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations
    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations
!    iGasJac = -9999   !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 
!    iJax = 600        !!! NO dQ PERTURBATIONS for Q jacs        NO dQ PERTURBATIONS for Q jacs 
    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations
    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations
    !! DEFAULT : leave this UNCOMMENTED if you DO NOT WANT PERTURBATIONS; comment if you want to test perturbations

    iJax2 = iJax+5    !!! can alter this to do dQ for many adjacent layers
    iJax2 = iJax+0    !!! can alter this to do dQ for many adjacent layers
! ************************************************************************

    DO iInt=1,iaNumLayer(iAtm)
!      print *,iInt,iaaRadLayer(iAtm,iInt)
      IF (iaaRadLayer(iAtm,iInt) .EQ. iJax) THEN
        iX = iInt
!        print *,'A',iInt,iaaRadLayer(iAtm,iInt),iJax
      END IF
      IF (iaaRadLayer(iAtm,iInt) .EQ. iJax2) THEN
        iX2 = iInt
!        print *,'B',iInt,iaaRadLayer(iAtm,iInt),iJax2
      END IF       
    END DO    
!    print *,'iJax,iJax2,iX,iX2 = ',iJax,iJax2,iX,iX2

    DO iInt=1,kProfLayer
      raTAmt(iInt)          = raaAmt(iInt,iGas)
      raTTemp(iInt)         = raaTemp(iInt,iGas)
      raTPress(iInt)        = raaPress(iInt,iGas)
      raTPartPress(iInt)    = raaPartPress(iInt,iGas)
      ! compute particle number density in number/cm3 (N/V = P/kT)
      raNumberDensity(iInt) = raTPress(iInt)*kAtm2mb*100.0/(kBoltzmann * raTTemp(iInt))*1e-6
      ! NLTE avg layer press in lower atm = kCARTA profile
      pProfNLTE(iInt) = raTPress(iInt)*kAtm2mb

      !!       print *,'BB',iInt,pprofNLTE(iInt),raTPress(iInt), &
      !!                   pprofNLTE(iInt)/(raTPress(iInt)*kAtm2mb)

      !!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !!c  this stuff is to test the temperature jacobians
      !!          raTTemp(iInt)=rDerivTemp+raTTemp(iInt)
      !!          raMixVertTemp(iInt)=rDerivTemp+raTTemp(iInt)
      !!          raMixVertTemp(iInt)=raTTemp(iInt)

      !!c jacob test temperature jacobian ... old
      !!          IF (iInt .EQ. iJax) THEN
      !!            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
      !!            raMixVertTemp(iInt)=raTTemp(iInt)
      !!          END IF
      !!          IF ((iInt .EQ. iJax) .AND. (iGas .EQ. 1)) THEN
      !!            rDummy        = raMixVertTemp(iInt)
      !!            raMixVertTemp(iInt) = rDummy+rDerivTemp
      !! c          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
      !! c          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
      !!            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
      !!          END IF

      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !c jacob test column gas  jacobian
      !          IF (iGas .EQ. 1) THEN
      !            rDummy        = raMixVertTemp(iInt)
      !            raMixVertTemp(iInt) = rDummy+rDerivTemp
      ! c          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
      ! c          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
      !            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
      !          END IF

      !c jacob test temperature jacobian ... new, could also be column
      !          rDerivTemp = 1.0
      !          IF (((iInt .GE. iJax) .AND. (iInt .LE. iJax2)) .OR. (iJax .EQ. -1)) THEN
      !            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
      !            raMixVertTemp(iInt)=raTTemp(iInt)
      !            print *,iJax,raaTemp(iInt,iGas),rDerivTemp,raMixVertTemp(iInt)
      !          END IF

      !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !c jacob test amount jacobian, orig
      !         IF ((iInt .EQ. iJax).AND.(iaGases(iGas) .EQ. iGasJac)) THEN
      !           raTAmt(iInt)=raaAmt(iInt,iGas)*(1.0+rDerivAmt)
      !           print *,'iJax, dq = ',iJax,raaAmt(iInt,iGas)*rDerivAmt
      !         END IF
      ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      !c jacob test amount jacobian, new
      !c jacob test amount jacobian, new
      !c jacob test amount jacobian, new
      !    depending on kJacobOutput (-1,2 for dr/dq/dBT/dq) and (0,1 for q dr/dq, q dBT/dq) 
      !    this varies rAdjust_rDerivAmt so that    rUseDQ = constant for each layer, so dQ/Q   varies layer by layer for the layers iJax..iJax2    if kJacobOutput = -1,2  dr/dq,     dBT/dq
      !                                             dQ/Q = constant for each layer,   so rUseDQ varies layer by layer for the layers iJax..iJax2    if kJacobOutput = 0,1   q dr/dq, q dBT/dq
      !
                IF ((iInt .GE. iJax) .AND. (iInt .LE. iJax2) .AND. (iaGases(iGas) .EQ. iGasJac)) THEN
       
                  IF (dq_Const_between_layers .LT. 0) THEN
                    !! dR/dq or dBT/dq
                    IF (iStartKeepTrack .LT. 0) THEN
                      iStartKeepTrack = +1
                      rDQ_Track = raaAmt(iInt,iGas)*rDerivAmt
                      rAdjust_rDerivAmt = 1.0
                    ELSEIF (iStartKeepTrack .GT. 0) THEN
                      rJunk = raaAmt(iInt,iGas)*rDerivAmt
                      rAdjust_rDerivAmt = rDQ_Track/rJunk
                    END IF
                  ELSE
                    !! q dR/dq or q dBT/dq
                    rAdjust_rDerivAmt = 1.0
                  END IF 
       
                  rUseDQ               = rDerivAmt * rAdjust_rDerivAmt
                  raTAmt(iInt)         = raaAmt(iInt,iGas) * (1.0+rUseDQ)
                  raTPartPress(iInt)   = raaPartPress(iInt,iGas) * (1.0+rUseDQ)
                  rDQ_Track_Cumulative = rDQ_Track_Cumulative + raaAmt(iInt,iGas) * rUseDQ
       
                  !!! recall kMaxLayer = 100 so we look 001-100 ..... but the atmosphere is eg from 1013 mb so we have iNumLayer = 97 and the radiating layers are (4,5,6 ... 100)
                  !!!    so when we read in the 97 layer jacobins, eg here index iInt = 6 will correspond to radiating layer 3, so the dq here is for J3 of 97
                  write(kStdErr,'(A,I3.3,1X,A,I3.3,A,1X,I3.3,A,1X,4(ES12.5))') &
                   'iIndex into (001-100) layers, <iXint> (actual radiating iaRadlayer of iNumlayer), RTP layer (out of 101): ', &
                   iInt,'<',iX + (iInt-iJax),'>', (kProfLayer+1)-iInt+1, &
                   '   Q,dq, dq/q, dq_cumulative :',raaAmt(iInt,iGas),raaAmt(iInt,iGas)*rUseDQ,rUseDQ,rDQ_Track_Cumulative
                END IF

      ! ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    END DO

!    print *,'nananana'
!    if (iaGases(iGas) .EQ. 2) call dostop

    RETURN
    end SUBROUTINE Set_Ref_Current_Profs

!************************************************************************
! this figures out CO2 mixing ratio
! since Antartic surface pressures can be as low as 500 mb, CO2.N2O,CO mixing ratios were failing if iBot=20
! this corresponds to raPresslevls(20) = 596 mb
! so subr Get_Temp_Plevs (in n_pth_mix.f) and subr compute_co2_mixratio (in kcartamisc.f)
! both needed to have iBot tweaked to layer 25 or above, which corresponds to raPresslevls(25) = 496 mmb
    SUBROUTINE compute_co2_mixratio(raaPress,raaPartPress,raaAmt,iNumLayer,rFracBot,rCO2MixRatio)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iNumLayer
    REAL :: raaPress(kProfLayer,kGasStore),rFracBot
    REAL :: raaPartPress(kProfLayer,kGasStore),raaAmt(kProfLayer,kGasStore)
! output
    REAL :: rCO2MixRatio

! local
    INTEGER :: iI,iBot,iTop,iCount,i1,i2,iLowestLay
    REAL :: rX,rF,cfac
          
! Conversion factor (equivalent mm of liquid water per molecules/cm^2)
! cfac = 1 (cm^3/gram) * ~18/Navadagro (grams/molecule) * 10 (mm/cm)
    cfac = 10 * 18.015 / (kAvog / 1000)

! >>>>>>>>>>>>>>>>>>>>>>>>>
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
      IF (iI == iLowestLay) rF = max(min(rFracBot,1.0),0.0)
      iCount = iCount + 1
      rX = rX + raaAmt(iI,1) * rF
    END DO
    rX = rX * cfac * (kAvog)   !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
    write(kStdWarn,*)' column amount water = ',rX,' mm'
         
! >>>>>>>>>>>>>>>>>>>>>>>>>
    iBot = 20                              !! assume profile has started from here !!!!
    iBot = max(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!
    iTop = kProfLayer - iBot
    write(kStdWarn,*) 'CO2 : Compute lower trop <ppmv> between lays ',iBot,iTop

    rX = 0.0
    iCount = 0
    DO iI = iBot,iTop
      iCount = iCount + 1
      rX = rX + raaPartPress(iI,2)/(raaPress(iI,2)-raaPartPress(iI,1))   !! wrt dry air
    END DO
    rX = rX * 1.0e6/iCount
    rCO2MixRatio = rX
    write(kStdWarn,*)'avg rCO2Mix = ',rX,' ppmv'

! >>>>>>>>>>>>>>>>>>>>>>>>>
    IF (kPlanet == 03) THEN
      iBot = 40  !! assume profile has started from here!!!!
      iTop = 70
      write(kStdWarn,*) 'O3 : Compute lower trop <ppmv> between lays ',iBot,iTop
  
      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,3)/(raaPress(iI,2)-raaPartPress(iI,1))   !! wrt dry air
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
        IF (iI == iLowestLay) rF = rFracBot
        iCount = iCount + 1
        rX = rX + raaAmt(iI,3) * rF
      END DO
      rX = rX * cfac * (kAvog) !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
      write(kStdWarn,*)' column amount ozone = ',rX,' du'
    END IF

! >>>>>>>>>>>>>>>>>>>>>>>>>
    IF (kPlanet == 03) THEN
      iBot = 20                              !! assume profile has started from here !!!!
      iBot = max(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!
      iTop = kProfLayer - iBot
      write(kStdWarn,*) 'N2O/CO/CH4 : Compute lower trop <ppmv> between lays ',iBot,iTop
  
      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,4)/(raaPress(iI,4)-raaPartPress(iI,1))   !! wrt dry air
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rN2OMix = ',rX,' ppmv'
    END IF
  
    rX = 0.0
    iCount = 0
    DO iI = iBot,iTop
      iCount = iCount + 1
      rX = rX + raaPartPress(iI,5)/(raaPress(iI,5)-raaPartPress(iI,1))   !! wrt dry air
    END DO
    rX = rX * 1.0e6/iCount
    write(kStdWarn,*)'avg rCOMix = ',rX,' ppmv'

    IF (kPlanet == 03) THEN  
      rX = 0.0
      iCount = 0
      DO iI = iBot,iTop
        iCount = iCount + 1
        rX = rX + raaPartPress(iI,6)/(raaPress(iI,6)-raaPartPress(iI,1))   !! wrt dry air
      END DO
      rX = rX * 1.0e6/iCount
      write(kStdWarn,*)'avg rCH4Mix = ',rX,' ppmv'
    END IF
  
    write(kStdWarn,*) 'assumed fracBot = ',rFracBot,max(min(rFracBot,1.0),0.0)
    write(kStdWarn,*) ' '

    RETURN
    end SUBROUTINE compute_co2_mixratio

!************************************************************************
! ----------------------------------------------------------------
!https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap02/datetime.f90
!  This program uses DATE_AND_TIME() to retrieve the system date
!  and the system time.  Then, it converts the date and time
!  information to a readable format.  This program demonstrates
!  the use of concatenation operator // and substring
! ----------------------------------------------------------------

   SUBROUTINE DateTimeX(callname)
   IMPLICIT   NONE
   include '../INCLUDE/TempF90/kcartaparam.f90'
   
   CHARACTER*(*) :: callname
   CHARACTER(LEN = 8)  :: DateINFO                 ! ccyymmdd
   CHARACTER(LEN = 4)  :: Year, Month*2, Day*2
   CHARACTER(LEN = 10) :: TimeINFO, PrettyTime*12  ! hhmmss.sss
   CHARACTER(LEN = 2)  :: Hour, Minute, Second*6

   CALL  DATE_AND_TIME(DateINFO, TimeINFO)

!  decompose DateINFO into year, month and day.
!  DateINFO has a form of ccyymmdd, where cc = century, yy = year
!  mm = month and dd = day

   Year  = DateINFO(1:4)
   Month = DateINFO(5:6)
   Day   = DateINFO(7:8)

!   WRITE(*,*)  'Date information -> ', DateINFO
!   WRITE(*,*)  '            Year -> ', Year
!   WRITE(*,*)  '           Month -> ', Month
!   WRITE(*,*)  '             Day -> ', Day

!  decompose TimeINFO into hour, minute and second.
!  TimeINFO has a form of hhmmss.sss, where h = hour, m = minute
!  and s = second

   Hour   = TimeINFO(1:2)
   Minute = TimeINFO(3:4)
   Second = TimeINFO(5:10)

   PrettyTime = Hour // ':' // Minute // ':' // Second

   WRITE(*,*)
!   WRITE(*,*)  'Time Information -> ', TimeINFO
!   WRITE(*,*)  '            Hour -> ', Hour
!   WRITE(*,*)  '          Minite -> ', Minute
!   WRITE(*,*)  '          Second -> ', Second
!   WRITE(*,*)  '     Pretty Time -> ', PrettyTime

!  the substring operator can be used on the left-hand side.
 
   PrettyTime = ' '
   PrettyTime( :2) = Hour
   PrettyTime(3:3) = ':'
   PrettyTime(4:5) = Minute
   PrettyTime(6:6) = ':'
   PrettyTime(7: ) = Second

!  WRITE(*,*)
!  WRITE(*,*)  '     Pretty Time -> ', PrettyTime 
   WRITE (kStdWarn,*) callname,' run ended ',Year,'/',Month,'/',Day,' at ',PrettyTime
   WRITE (kStdErr,*) callname,' run ended ',Year,'/',Month,'/',Day,' at ',PrettyTime   
   
   RETURN
   END SUBROUTINE DateTimeX
!************************************************************************
   SUBROUTINE DoDump_REF_ACTUAL_profiles(iNumGases,iaGases,iaProfFromRTP,raaPress,raaTemp,raaAmt,raaRAmt,iaNumLayer,iaaRadLayer)

   IMPLICIT   NONE
   include '../INCLUDE/TempF90/kcartaparam.f90'

   ! input
   INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
   INTEGER :: iNumGases,iaGases(kMaxGas),iaProfFromRTP(kMaxGas)
   REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore),raaRAmt(kProfLayer,kGasStore)
   REAL :: raaPress(kProfLayer,kGasStore)

   ! local
   INTEGER :: iGas,iL,iJ,iIOUN_Junk,iFileErr,i1,i2
   CHARACTER*80 :: caIOUN_Junk
   REAL :: ra1,ra2,rJunk1,rJunk2

   caIOUN_Junk = 'warning_prof.tmp'
   iIOUN_Junk = kTempUnit

   OPEN(UNIT=iIOUN_Junk,FILE=caIOUN_Junk,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=iFileErr)
   IF (iFileErr /= 0) THEN
     write(kStdErr,'(A,I4,A)')'error ',iFileErr,' opening ',caIOUN_Junk
     CALL DoSTOP
   END IF

   kTempUnitOPEN = +1
!    do iGas = 1,kGasStore
!      write(kTempUnit,*) iaGases(iGas),iaProfFromRTP(iGas)
!    END DO
!    DO iL = 1,kProfLayer
!      write(kTempUnit,'(I4,F12.5,F12.5)') ,iL,raaPress(iL,1),raaTemp(iL,1)
!    END DO
!    DO iL = 1,kProfLayer
!      write(kTempUnit,*) (raaAmt(iL,iGas),iGas=1,kGasStore)
!    END DO
!    DO iL = 1,kProfLayer
!      write(kTempUnit,*) (raaRAmt(iL,iGas),iGas=1,kGasStore)
!    END DO

   write(kTempUnit,'(A,I3)') '% numgases = ',iNumGases

   ra1 = 0.0
   ra2 = 0.0
   DO iGas = 1,kGasStore
      write(kTempUnit,'(A)') &
      '%'
      write(kTempUnit,'(A)') &
      '% iL     iJ         iGas  iaGases(iGas)  iaProf   raPress    raTemp      raAmt     raRAmt    raAmt/raaRamt'
      write(kTempUnit,'(A)') &
      '%    iaRadLayer(iL)          iGasID      FromRTP'
      write(kTempUnit,'(A)') &
      '%--------------------------------------------------------------------------------------------------------'

     DO iL = 1,iaNumLayer(1)
       iJ = iaaRadLayer(1,iL)
       IF (iaGases(iGas) .LT. 100) THEN
         rJunk1 = raaRAmt(iJ,iGas)         
         rJunk2 = raaAmt(iJ,iGas)/max(rJunk1,1.0e-15)
       ELSE
         rJunk1 = raaRAmt(iJ,1)         
         rJunk2 = raaAmt(iJ,iGas)/max(rJunk1,1.0e-15)
       END IF
       write(kTempUnit,'(5(I4,5X),2(F12.5),3(ES12.5))') & 
         iL,iJ,iGas,iaGases(iGas),iaProfFromRTP(iGas),& 
        raaPress(iJ,1),raaTemp(iJ,1),raaAmt(iJ,iGas),rJunk1,rJunk2
     END DO

     i1 = minval(iaaRadLayer(1,1:iaNumLayer(1)))+1
     i2 = maxval(iaaRadLayer(1,1:iaNumLayer(1)))
     IF (iaGases(iGas) .LT. 100) THEN
       rJunk1 = sum(raaRAmt(i1:I2,iGas))         
       rJunk2 = sum(raaAmt(i1:i2,iGas))/max(rJunk1,1.0e-15)
     ELSE
       rJunk1 = sum(raaRAmt(i1:I2,1))         
       rJunk2 = sum(raaAmt(i1:i2,iGas))/max(rJunk1,1.0e-15)
     END IF
     write(kTempUnit,'(A,I4,I4,I4,ES12.5,ES12.5,F12.5)') &
       '% iG, iGasID, fromRTP file, column Q, column Qref,ratio :  ',   &
        iGas,iaGases(iGas),iaProfFromRTP(iGas),                         & 
        sum(raaAmt(i1:i2,iGas)),rJunk1,rJunk2

   END DO

   CLOSE(UNIT=iIOUN_Junk)
   kTempUnitOPEN = -1

   RETURN
   END SUBROUTINE DoDump_REF_ACTUAL_profiles
!************************************************************************

END MODULE kcartamisc
