! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

!! use must precede IMPLICIT and INCLUDE statements

MODULE kcartamisc

USE basic_common
USE freqfile
USE spline_and_sort
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
!************************************************************************

CONTAINS

!************************************************************************

! check for isnan
    LOGICAL FUNCTION isnan_real(x)

    REAL :: x
    LOGICAL :: isnan

    print *,x,x+1.0e0
    IF (x+1.0e0 == x) THEN
        isnan = .TRUE. 
    ELSE
        isnan = .FALSE. 
    ENDIF
    isnan = .TRUE. 

    IF (x >= 0.0) THEN
        print *,'gt',x
        isnan = .FALSE. 
    ELSEIF (x <= 0.0) THEN
        print *,'lt',x
        isnan = .FALSE. 
    ENDIF

    isnan_real = isnan
    RETURN
    end FUNCTION isnan_real

!************************************************************************
    LOGICAL FUNCTION isnan_double(x)

    logical :: isnan
    double precision :: x

    print *,x,x+1.0d0
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

    include '../INCLUDE/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'

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
            raFracDeltaP(iJ) = &
            abs(raPresslevels(iJ)-PLEV_KCARTADATABASE_AIRS(iJ))/ &
            PLEV_KCARTADATABASE_AIRS(iJ)
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
        IF ((rMin <= 1.0e-5) .AND. (rMax <= 5.0e-4) .AND. &
        (rMaxAv <= 5.0e-3)) THEN
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

    include '../INCLUDE/kcartaparam.f90'

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

    CHARACTER(80) :: caFName
    REAL :: rMeanT,rMeanA,rMeanP,rMeanPP,rDirectPPMV0,rDirectPPMV1
    REAL :: rCO2ppmv
    INTEGER :: iI,iJ,iGasID,iError,iIPMIX

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
        IF ((raaPress(iI,iGasID) > 0.0) .AND. &
        (raaPress(iI,iGasID) < 1.5)) THEN
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
    !        CALL DoStop
    END IF
          
    rCO2Mult = raaMix(iIpmix,iGasID)

!! open the true RefProf at kCO2ppmv
    caFName = kCO2ppmvFile
    write(kStdWarn,*) 'Reading CO2 Reference profile from ',caFName
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)

!! open the supposed RefProf
    CALL FindReferenceName(caFName,iGasID,-1)
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt1, &
    raR100Temp1,raR100Press1,raR100PartPress1,iError)

! -------------------------
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
    !        CALL DoStop
    END IF

    100 FORMAT(A25,A80)
    101 FORMAT(A65,F12.8)
    IF ((rMeanA/kMaxLayer <= 0.9995) .OR. &
    (rMeanA/kMaxLayer >= 1.0005)) THEN
        write(kStdErr,100) 'v0 : kCO2ppmvFile =   ', kCO2ppmvFile
        write(kStdErr,101) 'mean rCO2ppmv (raCO2pp/raP) v0 from input TRUErefprof2         = ',rDirectPPMV0
        write(kStdErr,100) 'v1 : ref CO2 profile = ', caFName
        write(kStdErr,101) 'mean rCO2ppmv (raCO2PP/raP) v1 from input refprof2             = ',rDirectPPMV1
        write(kStdErr,*) 'oops check_co2ppmv : rMeanA is bad',rMeanA/kMaxLayer
        CALL DoStop
    END IF
    IF ((rMeanP/kMaxLayer <= 0.9995) .OR. &
    (rMeanP/kMaxLayer >= 1.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMeanP is bad',rMeanP/kMaxLayer
        CALL DoStop
    END IF
    IF  ((rMeanPP/kMaxLayer <= 0.9995) .OR. &
    (rMeanPP/kMaxLayer >= 1.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMeanPP is bad',rMeanPP/kMaxLayer
        CALL DoStop
    END IF
    IF ((rMeanT/kMaxLayer <= -0.0005) .OR. &
    (rMeanT/kMaxLayer >= +0.0005)) THEN
        write(kStdErr,*) 'oops check_co2ppmv : rMean deltaT is bad',rMeanT/kMaxLayer
        CALL DoStop
    END IF

!!! now check the NLTE weak background in Upper Layers
    write (kStdWarn,*) '    checking weak backgnd upper atm ppmvs'
    caFName = kCO2ppmvFileBackUA
    CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0, &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)
!! can directly figure out the ppmv in the "what we hope is true" file
    rDirectPPMV0 = 0.0
    DO iI = 1,kMaxLayer/20
        rDirectPPMV0 = rDirectPPMV0+abs(raR100PartPress0(iI)/raR100Press0(iI))
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
! this subroutine checks to see if the gasID is 1-36 or 101-102 or 103
    INTEGER FUNCTION MainGas(iGasID)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iGasID

    INTEGER :: iI,i1,i2
           
    iI = -1        !assume this is a main gas (1-36, or 101-102)

    i1 = -1
    i2 = -1

    IF ((iGasID >= 1) .AND. (iGasID <= kGasComp)) THEN
        i1 = 1          !! main gas
    ELSEIF (iGasID == kNewGasHi+1) THEN
        i1 = 1          !! heavy water
    ELSEIF ((iGasID >= kNewGasLo) .AND. (iGasID <= kNewGasHi)) THEN
        i2 = 1          !! water continuum
    END IF
          
    IF ((i1 == 1) .OR. (i2 == 1)) THEN
        iI = 1
    END IF

    IF ((i2 == 1) .AND. kCKD < 0) THEN
        write(kStdWarn,*) 'Cannot have gases 101,102 with CKD turned off'
        write(kStdErr,*) 'Cannot have gases 101,102 with CKD turned off'
        Call DoSTOP
    END IF

    MainGas = iI

    RETURN
    end FUNCTION MainGas

!************************************************************************
! this suroutine sets up the current print options
    SUBROUTINE SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,iType, &
    iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
    iNumGases,iNpmix,iaNumLayer,iAtm)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
! this sets the temperature variation after nm_radnce is read in
    SUBROUTINE SetkTemperVary(iTemperVary)

    IMPLICIT NONE
    include '../INCLUDE/kcartaparam.f90'

! input param
    INTEGER :: iTemperVary      !!! from namelist file

! local var
    INTEGER :: iConstOrVary

! this is TEMPERATURE variation in layer
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileUPWELL_LINEAR_IN_TAU
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileDNWELL_LINEAR_IN_TAU
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
!      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      kTemperVary = +1     !!!temperature in layer varies
!      kTemperVary = +2     !!!temperature in layer varies linearly, simple
!      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
!      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM,
!                           !!!  no O(tau) approx, very similar to kTemperVary=4
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM,
!                           !!!  debugged for small O(tau), used with EliMlawer 12/2015
!      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!!  x/6 as x-->0 compared to kTemperVary = +42 *****
!      IF (kFlux .LE. 0) THEN
!        kTemperVary = -1     !!! temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      ELSE
!        kTemperVary = +43    !!! temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!! x/6 as x-->0 compared to kTemperVary = +42 ****
!      END IF

    kTemperVary = iTemperVary
          
    iConstOrVary = -1   !! if do flux, do linear vary T with tau
    iConstOrVary = +1   !! if only RaDTrans, then do constant T in layer, default SARTA/kCARTA for RT only
          
    IF (kFlux <= 0) THEN
        IF (iConstOrVary > 0) THEN
            kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
            write(kStdWarn,*) 'kFlux <= 0 so set kTemperVary = -1'
        ELSEIF (iConstOrVary < 0) THEN
            kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
        !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
            write(kStdWarn,*) 'kFlux < 0 but set kTemperVary = 43'
        END IF
    ELSEIF (kFlux > 0) THEN
        kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
    !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
        write(kStdWarn,*) 'kFlux > 0 so set kTemperVary = 43'
    END IF

!!! new, do what the user wishes!!!
    IF ((kFlux <= 0) .AND. (iTemperVary > 0)) THEN
        kTemperVary = +43
    END IF
          
!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
! kTemperVary = +43
! F (iTemperVary .NE. kTemperVary) THEN
!  write(kStdWarn,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  write(kStdErr,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  kTemperVary = iTemperVary
! ND IF
!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
                
    IF (kTemperVary == -1) THEN
        write(kStdWarn,*) 'kTemperVary = -1     !layer temp constant (SARTA DEFAULT)'
    ELSEIF (kTemperVary == +1) THEN
        write(kStdWarn,*) 'kTemperVary = +1     !layer temp varies'
    ELSEIF (kTemperVary == +2) THEN
        write(kStdWarn,*) 'kTemperVary = +2     !layer temp varies linearly, simple v2'
    ELSEIF (kTemperVary == +3) THEN
        write(kStdWarn,*) 'kTemperVary = +3     !layer temp varies linearly, ala LBLRTM v3'
    ELSEIF (kTemperVary == +4) THEN
        write(kStdWarn,*) 'kTemperVary = +4     !layer temp varies linearly, ala LBLRTM v4 O(tau^2)'
    ELSEIF (kTemperVary == +41) THEN
        write(kStdWarn,*) 'kTemperVary = +41    !layer temp varies linearly, ala LBLRTM v4 (Pade)'
    ELSEIF (kTemperVary == +42) THEN
        write(kStdWarn,*) 'kTemperVary = +42    !layer temp varies linearly, ala LBLRTM v4 O(tau)'
    ELSEIF (kTemperVary == +43) THEN
        write(kStdWarn,*) 'kTemperVary = +43    !layer temp varies linearly, ala LBLRTM v4 O(tau) -> tau/6'
    ELSE
        write(kStdErr,*)'kTemperVary = ',kTemperVary,'unknown option'
        CALL DoStop
    END IF

    iaaOverrideDefault(2,1) = kTemperVary
          
    RETURN
    end SUBROUTINE SetkTemperVary

!************************************************************************
! this subroutine does some more initializations
    SUBROUTINE SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
! raaRTPCloudParams0(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   from rtpfile
! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   after kcarta resets
! this gets set in rtp_interface.f
    DO iFileIDLO = 1,7
        raaRTPCloudParams0(1,iFileIDLO) = -1.0
        raaRTPCloudParamsF(1,iFileIDLO) = -1.0
        raaRTPCloudParams0(2,iFileIDLO) = -1.0
        raaRTPCloudParamsF(2,iFileIDLO) = -1.0
    END DO
          
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
    DO iFileIDLo = 1,kMixFilRows
        DO iFileIDHi = 1,kGasStore
            raaMix(iFileIDLo,iFileIDHi) = 0.0
        END DO
    END DO

! initialize kaaNumVectors
    DO iFileIDLo = 1,kMaxGas
        DO iFileIDHi = 1,kMaxLayer
            kaaNumVectors(iFileIDLo,iFileIDHi) = 0
        END DO
    END DO

    RETURN
    end SUBROUTINE SomeMoreInits

!************************************************************************
! this subroutine inits file unit numbers
    SUBROUTINE InitializeFileUnits
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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

    include '../INCLUDE/kcartaparam.f90'

! this is the driver file name from the command line arguments
    CHARACTER(80) :: caDriverName

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
    CHARACTER(80) :: caComment,caComment1
    CHARACTER(80) :: caLogFile,caLogFile1
    REAL :: raaOp(kMaxPrint,kProfLayer),raaOp1(kMaxPrint,kProfLayer)

    NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter, &
    iaGPMPAtm,iaNp,iaaOp,raaOp

    iIOun=kStdDriver
    IF (iIOUN /= 5) THEN
        OPEN(UNIT = iIOun,FILE = caDriverName,STATUS='OLD',IOSTAT=iErr)
        IF (iErr /= 0) THEN
            write (kStdErr,*) 'in subroutine ErrorLogName, error reading'
            write (kStdErr,*) 'namelist file to find name of logfile ... '
            WRITE(kStdErr,1070) iErr, caDriverName
            1070 FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80)
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

    include '../INCLUDE/kcartaparam.f90'

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
        write(kStdWarn,30) iDummy,iaPrinter(iDummy),iaGPMPAtm(iDummy), &
        iaNp(iDummy)
        write(kStdWarn,*) 'paths to be printed : (if numpaths=-1,print all)'
        write(kStdWarn,*)(iaaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy))
        IF (iaPrinter(iDummy) == 3) THEN
            write(kStdWarn,*)(raaOp(iDummy,iOutNum), &
            iOutNum=1,iaNp(iDummy))
            write(kStdWarn,*)(raaUserPress(iDummy,iOutNum), &
            iOutNum=1,iaNp(iDummy))
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

    include '../INCLUDE/kcartaparam.f90'
     
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
        ELSE IF (iaGases(2) == 2) THEN
            iCO2 = 2
            write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!'
        ELSE !!!for some strange reason, no CO2 present
            iCO2 = 1
            write(kStdWarn,*) 'Temperature of Gas number 1 will mimic CO2!!'
        END IF
         
        CALL GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, &
        iNpmix,iCO2)
         
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

    include '../INCLUDE/kcartaparam.f90'

! caDriverName is the name of the driver file to be processed
    CHARACTER(80) :: caDriverName
! caOutName is the name of the unformatted output file name
! integer iOutFileName tells whether or not there is a driver name, or
! dump output to Unit 6
    INTEGER :: iOutFileName
    CHARACTER(80) :: caOutName
! caJacobFile is the name of the unformatted output file name for Jacobians
    CHARACTER(80) :: caJacobFile
! this tells if we have MS Product ie no command line stuff!
    INTEGER :: iMicroSoft
! this is the number of args
    INTEGER :: iargc

    INTEGER :: iDummy,iError

    IF (iMicroSoft > 0) THEN
                 
    ! o command line options .. do this!
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
    ! se command line stuff
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

    include '../INCLUDE/kcartaparam.f90'

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
            WRITE(kStdErr,*) 'Error in Ref Profile for Gas', &
            iaGases(iGas)
            WRITE(kStdErr,*) 'Layer ',iI,' has negative amount',raRAmt(iI)
            CALL DoStop
        END IF
         
        IF (raRTemp(iI) < 0.0) THEN
            WRITE(kStdErr,*) 'Error in Ref Profile for Gas', &
            iaGases(iGas)
            WRITE(kStdErr,*) 'Layer ',iI,' has negative tempr',raRTemp(iI)
            CALL DoStop
        END IF
         
        IF (raRPress(iI) < 0.0) THEN
            WRITE(kStdErr,*) 'Error in Ref Profile for Gas', &
            iaGases(iGas)
            WRITE(kStdErr,*) 'Layer ',iI,' has negative pressure',raRPress(iI)
            CALL DoStop
        END IF
         
        IF (raRPartPress(iI) < 0.0) THEN
            WRITE(kStdErr,*) 'Error in Ref Profile for Gas', &
            iaGases(iGas)
            WRITE(kStdErr,*) 'Layer ',iI,' has negative part press', &
            raRPartPress(iI)
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

    include '../INCLUDE/kcartaparam.f90'

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

    DO iInt = 1,kProfLayer
        raRAmt(iInt) = raaRAmt(iInt,iGas)
        raRTemp(iInt) = raaRTemp(iInt,iGas)
        raRPress(iInt) = raaRPress(iInt,iGas)
        raRPartPress(iInt) = raaRPartPress(iInt,iGas)
    END DO

    RETURN
    end SUBROUTINE SetReference

!************************************************************************
! this subroutine close units
    SUBROUTINE TheEnd(iaGases,iNumGases,iaList,raFiles)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
          
    INTEGER :: iNumGases,iaGases(kMaxGas),iaList(kNumkCompT)
    REAL :: raFiles(kNumkCompT)
    INTEGER :: iI,iJ,iG,iaSum(kMaxGas),iaCount(kMaxGas),iSum,iaJunk(kMaxGas)

    write(kStdWarn,*) '*****************************************'
    write(kStdWarn,*) 'kComp stats ..............'
    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Freq chunks processed'
    write(kStdWarn,*) (raFiles(iaList(iJ)),iJ=1,kOuterLoop)

    write(kStdWarn,*) 'Stats of compressed vecs per gas per chunk'
    DO iI = 1,kMaxGas
        iaSum(iI) = 0
        iaCount(iI) = 0
    END DO
          
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
            write(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG), &
            &                       1.0*iaSum(iG)/iaCount(iG)
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
        write(kStdWarn,*) '  Chunk = ',raFiles(iaList(iJ)), ' numgases = ',iSum,' gasIDs are .... '
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

    include '../INCLUDE/kcartaparam.f90'

    REAL :: raaAb(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFreq

    DO iLay = 1,kProfLayer
        DO iFreq = 1,kMaxPts
            raaAb(iFreq,iLay) = 0.0
        END DO
    END DO

    RETURN
    end SUBROUTINE InitializeReal

!************************************************************************
! this subroutine initializes all the rows of the
! (REAL) array of mixed paths
    SUBROUTINE InitializeRealMP(raaAb,iNpmix)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: raaAb(kMaxPts,kMixFilRows)
    INTEGER :: iNpMix

    INTEGER :: iLay,iFreq

    DO iLay = 1,iNpMix      !note : initialize only wot is necessary
        DO iFreq = 1,kMaxPts
            raaAb(iFreq,iLay) = 0.0
        END DO
    END DO

    RETURN
    end SUBROUTINE InitializeRealMP

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs
    SUBROUTINE initialize(daaAb)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    DOUBLE PRECISION :: daaAb(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFreq

    DO iLay = 1,kProfLayer
        DO iFreq = 1,kMaxPts
            daaAb(iFreq,iLay) = 0.0
        END DO
    END DO

    RETURN
    end SUBROUTINE initialize

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs
! pretty much the same as above routine
    SUBROUTINE initializeJAC(daaAb)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    DOUBLE PRECISION :: daaAb(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iLay,iFreq

    DO iLay = 1,kProfLayerJac
        DO iFreq = 1,kMaxPtsJac
            daaAb(iFreq,iLay) = 0.0
        END DO
    END DO

    RETURN
    end SUBROUTINE initializeJAC

!************************************************************************
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa
! and saves it in an overall AbsMatrix raaa
    SUBROUTINE DoSet(daaGasAbCoeff,raaaGasAbCoeff,iCount,iDoAdd)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iCount     = which of the iNumGases are being processed
! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPtsJac,kProfLayerJac)
    REAL :: raaaGasAbCoeff(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    INTEGER :: iCount,iDoAdd

! local variables
    INTEGER :: iLay,iFr

    IF (iDoAdd > 0) THEN
        DO iLay = 1,kProfLayerJac
            DO iFr = 1,kMaxPtsJac
            !          raaaGasAbCoeff(iCount,iFr,iLay) = real(daaGasAbCoeff(iFr,iLay))
                raaaGasAbCoeff(iCount,iFr,iLay) = daaGasAbCoeff(iFr,iLay)
            END DO
        END DO
    ELSEIF (iDoAdd <= 0) THEN
        DO iLay = 1,kProfLayerJac
            DO iFr = 1,kMaxPtsJac
                raaaGasAbCoeff(iCount,iFr,iLay) = 0.0
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE DoSet
!************************************************************************c
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa
    SUBROUTINE DoDtoR(daaGasAbCoeff,raaGasAbCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iLay,iFr

    DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
            raaGasAbCoeff(iFr,iLay) = real(daaGasAbCoeff(iFr,iLay))
        END DO
    END DO

    RETURN
    end SUBROUTINE DoDtoR

!************************************************************************
! this subroutine converts the layer iL to 0
    SUBROUTINE ZeroLayer(raaGasAbCoeff,iL)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raaGasAbs = 2d matrix that save ALL abs coeffs for current 25 cm-1 chunk
    REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)
    INTEGER :: iL

    INTEGER :: iFr

    DO iFr = 1,kMaxPts
        raaGasAbCoeff(iFr,iL) = 0.0
    END DO

    RETURN
    end SUBROUTINE ZeroLayer

!************************************************************************
! this subroutine checks parameters in kcartaparam.f90 ... abort if they
! do not make sense .. the parameters are set in kcartaparam.f90
    SUBROUTINE CheckKCARTAParameters

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
! set the default parameter values, for those that are not set in *PARAM
! read the parameter file to set parameters that have optional values
    SUBROUTINE SetDefaultParams

! NOTE !!!! also double check subroutine EXTRAPAR in strings2.f
! NOTE !!!! also double check subroutine SETDEFAULTPARAMS in misc.f

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iI,iJ

! acos(3/5) * 180/kPi
    kThermalAngle = 53.1313010235415598
          
! set default values here
! kLayer2Sp
!     -2     Layer transmittance            t(i)=exp(-k(i))
!     -1     Layer Optical depth            k(i)
!      1     Layer-to-Space Optical Depth   k2s(i)=sum(j=i,n)(k(j))
!      2     Layer-to-Space transmittance   t2s(i)=sum(j=i,n)exp(-k(j))
!      3     Layer-to-Ground Optical Depth  k2s(i)=sum(j=1,i)(k(j))
!      4     Layer-to-Ground transmittance  t2s(i)=sum(j=1,i)exp(-k(j))
    kLayer2Sp    = -1

! kCKD      == -1,00,21,23,24 sets which continuum version calculation to do
!              -1 : no continuum
! ----------------------------------------------------------------------------
! AER-CKD      01 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD version of Dec 2002
!       AIRS   02 : version 02    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Aug 2003
!       AIRS   03 : version 03    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Jan 2004
!       AIRS   04 : version 04    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!       AIRS   05 : version 05    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     (so far it looks like CKD v4)
!                                     On top of that it puts Dave Tobin's cs,cf
!                                     between 1300-1800 cm-1
!       AIRS   06 : version 06    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     extremely similar to CKD 4, except have
!                                     fixed the problem at 850 cm-1 which is
!                                     really HNO3 and not a continuum problem
! In asl:/carrot/s1/strow/Tobin/Tobin_radish/NISTdata2/New/, use
! cs0_tsl0_hr.mat and cf0_tsl0_hr.mat.  makecons.m with flag=7
! looks like it should load things in correctly.
! ----------------------------------------------------------------------------
! AER-STD      25 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Sep 2011
! ----------------------------------------------------------------------------
! AER-STD      27 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Feb 2016
! ----------------------------------------------------------------------------
! old versions of AER CKD
! AER-CKD      00 : version 00
! AER-CKD      21 : version 2.1
! AER-CKD      23 : version 2.3
! AER-CKD      24 : version 2.4
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! these were our interim versions (made by us, basically MT_CKD 1)
! no longer in use!

!    RAL       12 : self = version 2.4, foreign = blend of 0.0 and 2.4
!                                       from 0-1575, 1625-3000, use v2.4
!                                       from 1575-1625 linearly blend v2.4,0.0
!                                       using a triangle
!    RAL       13 : self = version 2.3, foreign = dave tobin's v2.3
!    RAL       90 : self = version 2.4, foreign = blend of 2.4,Tobin's thesis
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use Tobins thesis
!    RAL       50 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst50.m used linear tempr interpolation
!    RAL       51 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst51.m used power tempr interpolation
!                                       Need to be careful ... I put in Dave
!                                       Tobin's fixes for CS in the 600-1100
!                                       cm-1 range; see mst51_tobin.m in
!                                       the SPECTRA directory (~0.9*ckd24)
!    RAL       52 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst52.m used power tempr interpolation
!                                       and says CF(296)=CF(243); so use the
!                                       foreign broadened 243 data to get CS243
!    RAL       55 : self, foreign       same as 51 above, but uses
!                                a) CS fudge factor of 0.87 for v <= 1137 cm-1
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       56 : self, foreign       same as 51 above, but uses
!                                a) John Taylor fudge between 600-1200 cm-1
!                                   instead of Dave Tobin fudge
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                   these are better than the above fudges
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       60 : self, foreign     is a hybrid of 51,55 done by scott
!                                     so that bias errors are reduced, as
!                                     a function of water column amount.
!                                     eventually will include tuning coeffs
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25]  [] ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

    kCKD = 25

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
    kGasTemp = -1

! kLongOrShort == whether the user wants to save ALL header info (+1) or
!                 just a portion (-1) or output just the basic results (0)
    kLongOrShort = 1

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs
    kActualJacs = -1    !!! default

! kActualJacsT = -1, kActualJacsB = -1
! if we set kActualJacs = 100, then we can also set the layers of the column
! that we want to perturb; default = -1/-1 means all layers (kActualJacs = 100)
! else kActualJacsABCXYZ sets these as kActualJacsT=ABC, kActualJacsB=XYZ
    kActualJacsT = -1
    kActualJacsB = -1

! kJacobOutput == -1 if we output d(radiance)/dq,d(radiance)/dT
!                  0 if we output d(radiance)/dq * q, d(radiance)/dT
!                  1 if we output d(BT)/dq * q, d(BT)/dT
    kJacobOutput = 1

! kFlux == -1 if we do not want flux/OLR computations or output NLTE
!                                                        Planck modifiers
!       ==  1 if we want DNWELL flux at every layer      units = mW m-2
!       ==  2 if we want flux computations : output      units = kelvin day-1
!       ==  3 if we want UPWELL flux at every layer      units = mW m-2
!       ==  4 if we want outgoing OLR only at TOA,       units = mW m-2
!       ==  5 if we want outgoing OLR at TOA, ILR at GND units = mW m-2
!       ==  6 if we want DNWELL and UPWELL flux at each layer units = mW m-2
!   --> ==  0 if we want to output NLTE Planck modifiers
    kFlux = -1

! kPlanckOut == -1 if we do not want to output Planck modifiers, 1 if we do
    kPlanckOut = -1

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
! kSurfTemp = -1 == want to use user supplied surface temp in *RADNCE
!              1 == want to use user supplied surface temp in *RADNCE as
!                     an offset to pressure interpolated temperature
! so in above RTP cases if kSurfTemp < 0 : use raTSurf(iI)
!                          kSurfTemp > 0 : use raTSurf(iI) + InterpedTemp
    kSurfTemp = -1

! kRTP = -5  : read LBLRTM style LAYERS profile (edited TAPE 6); set atm from namelist
! kRTP = -6  : read LBLRTM style LEVELS profile (       TAPE 5); set atm from namelist
! kRTP = -10 : read TEXT style LEVELS   profile; set atm from namelist
! kRTP = -2  : read GENLN4 style LAYERS profile; set atm from namelist
! kRTP = -1  : read old style kLAYERS   profile; set atm from namelist
! kRTP =  0  : read RTP style kLAYERS   profile; set atm from namelist
! kRTP = +1  : read RTP style kLAYERS   profile; set atm from RTP file
! kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
    kRTP = 1

! kTempJac == -2 if we only want d/dT(planck) in temperature jacobian
!             -1 if we only want d/dT((1-tau)(tau->sp)) in temp jacobian
!              0 if we want complete d/dT(planck) + d/dT(tau) in temp jac
    kTempJac = 0

! the following cannot be controlled by the user using *PARAMS
! all the radiance parameters and the Jacobian parameter
! kJacobian ==  1 if analytic Jacobians are to be computed
!              -1 if we do the standard Genln2 computations w/o Jacobians
    kSolar        = 1     !turn on solar
    kSolarAngle   = 0.0   !solar angle
    kSolarRefl    = -1.0  !use (1-ems)/pi
    kThermal      = 0     !use fast diffusive approx
    kThermalAngle = -1.0  !use acos(3/5) in upper layers
    kThermalJacob = 1     !use thermal backgnd in Jacobians

    kJacobian = -1        !do not do Jacobians
    kScatter  = -1        !do not do scattering computations

    k100layerCloud = -1   !assume rtp file does NOT have 100 layer cloud

! 2016
! allow nm_params to define defaults
!   GENERAL iaDefaults(1,:) = iSARTAChi   iSPlineType iCO2Chi  iMatlabORf77   iWhichScatterCode iMethod
!   RT      iaDefaults(2,:) = iGaussPts iGaussQuad  iSnell      iInterpType  iWhichRT  (kTemperVary set in mn_radnce)
!   NLTE    iaDefaults(3,:) = iCurrent    iTalk       iTestGenln2  iNoPressureShiftCO2 iQtips_H98
!                             iLinearOrSpline iDoCO2Continuum iMethod
!   TAPE5/6     iaDefaults(5,:) = iReplaceZeroProf iAIRS101_or_LBL_levels IPLEV iAddLBLRTM
!      INTEGER iaaOverrideDefault(8,10)
!      COMMON/comBlockDefault/iaaOverrideDefault
    DO iI = 1,4
        DO iJ = 1,10
            iaaOverrideDefault(iI,iJ) = -9999
        END DO
    END DO
          
! GENERAL
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
          
! RadTrans
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
    iaaOverrideDefault(2,7) = -1    !!! iUseSnell = -1 for No  Snell law raytrace plus layer curvature effects, similar to SARTA (default)
!!!           = +1 for Yes Snell law raytrace plus layer curvature effects
!!!           = 0  for No  Snell law raytrace NO   layer curvature effects
!!!   see SUBR FindLayerAngles in rad_angles.f
    iaaOverrideDefault(2,8) = +1    !!! iInterpType = +1 to turn (pav,Tav) into (plevs,Tlevs), only used if kTemperVary = 43
!!!   see SUBR Get_Temp_Plevs in n_pth_mix.f
    iaaOverrideDefault(2,9) = -1    !!! iLBLRTM_highres = -1 do not estimate/fix problems because use 0.0025 cm-1, when kTemperVary = 43 << DEFAULT>>
!!!   see SUBR rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMISS in rad_main.f
    iaaOverrideDefault(2,10) = 5    !!! kWhichScatterCode = 5 for PCLSAM (Default)
!!!   0 for ABS clouds, 2 for RTPSEC, 3 for DISORT
      
! n_layers_lblrtm.f and n_pth_mix.f  TAPE5/6
    iaaOverrideDefault(3,1) = -1    !!! iAIRS101_or_LBL_levels use LBLRTM, not AIRS 101 levels, for integration
    iaaOverrideDefault(3,2) = +1    !!! iReplaceZeroProf = +1 to add in profiles TAPE5 does not have
    iaaOverrideDefault(3,3) = -1    !!! iAddLBLRTM = -1 when gas profile missing from TAPE5/6, do not add it in
!!!   in n_pth_mix.f

    RETURN
    end SUBROUTINE SetDefaultParams

!************************************************************************
! now check parameters in *PARAM
! this subroutine checks parameters in *PARAMS ... abort if they
! do not make sense ..
    SUBROUTINE CheckParams

    IMPLICIT NONE
    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: i0,iT,iB,iJ,iGah,iConstOrVary
    CHARACTER(9) :: iIOUN9

    write(kStdWarn,*) 'checking parameters (from *PARAMS) .... '
    write(kStdWarn,*) '  '

    IF ((iabs(kLayer2Sp) /= 1) .AND. (iabs(kLayer2Sp) /= 2) .AND. (kLayer2Sp /= 3) .AND. (kLayer2Sp /= 4)) THEN
        write(kStdErr,*) 'In *PARAMS, need kLayer2Sp = +/-1,+/-2,+3,+4'
        write(kStdErr,*) 'kLayer2Sp == do layer-to-space calc or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF
    IF (iabs(kGasTemp) /= 1) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kGasTemp = +/- 1'
        write(kStdErr,*) 'kGasTemp = use CO2 temperature profile or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

! CKD    releases are 0,21,23,24
! MT_CKD releases are 1,25
! our modifications are 2,3,4,5 (derived from 1)
!                       51,55,60 (derived from analysing RAL data)
! and various 12,13,50,52,56 which might be gotten rid of eventually
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25] [] ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

    IF ((kCKD /= -1) & &
!!! the std CKD pre-2002 versions
     .AND. (kCKD /= 0) .AND. (kCKD /= 21) .AND. (kCKD /= 23) .AND. (kCKD /= 24) & &
!!! these are MT_CKD1 and research versions from AIRS data
     .AND. (kCKD /= 1) .AND. (kCKD /= 4) .AND. (kCKD /= 6) &
     .AND. (kCKD /= 25) .AND. (kCKD /= 27)) &
    THEN
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

    IF (iabs(kLongOrShort) > 2) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kLongOrShort = -2,-1,0,+1,+2'
        write(kStdErr,*) 'kLongOrShort = complete header info (+1) or not (-1),  long warning.msg'
        write(kStdErr,*) 'kLongOrShort = complete header info (+2) or not (-2), short warning.msg'
        write(kStdErr,*) '               or file containing results only (0)'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs

    IF ((kActualJacs /= -1)  .AND. (kActualJacs /= 20)   .AND. &
    (kActualJacs /= 30)  .AND. (kActualJacs /= 40)   .AND. &
    (kActualJacs /= 50)  .AND. (kActualJacs /= 100)  .AND. &
    (kActualJacs /= -2)  .AND. (kActualJacs /= 32)   .AND. &
    (kActualJacs /= 102) .AND. (kActualJacs < 102)) THEN
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

    IF (kActualJacs > 102) THEN
        IF (kActualJacs < 100000001) THEN
            write(kStdErr,*) 'too few characters in kActualJacs : ',kActualJacs
            write(kStdErr,*) 'check it !'
            CALL DoStop
        END IF
        IF (kActualJacs > 102999999) THEN
            write(kStdErr,*) 'too many characters in kActualJacs : ',kActualJacs
            write(kStdErr,*) 'max possible value is 102999999'
            CALL DoStop
        END IF

        i0 = kActualJacs
    !! i0 better be 9 characters long
        write(iIOUN9,99) kActualJacs
        iJ = 9
        10 CONTINUE
        IF ((iIOUN9(iJ:iJ) /= ' ') .AND. (iJ >= 1)) THEN
        !         write(kStdErr,*) 'good',iJ,iIOUN9(iJ:iJ),kActualJacs
            iJ = iJ - 1
            GOTO 10
        ELSEIF (iJ > 0) THEN
            iGah = iJ
            write(kStdErr,*) 'space in  kActualJacs at ',iJ,iIOUN9(iJ:iJ)
        END IF
        IF (iGah > 0) THEN
            write(kStdErr,*) 9-iGah,' chars in kActualJacs = ',kActualJacs
            write(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
            write(kStdErr,*) '  or 100 or 100ABCXYZ'
            write(kStdErr,*) '  or 102 or 102ABCXYZ'
            write(kStdErr,*) 'need 9 characters 10E ABC XYZ .. try again'
            CALL DoStop
        END IF

        write(kStdWarn,*) 'kActualJacs passed test ... '
        IF (kActualJacs <= 102000000) THEN
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

        IF (iT < iB) THEN
            iJ = iT
            iT = iB
            iB = iJ
        END IF
        IF (iT > kProfLayer) THEN
            write(kStdWarn,*) 'IT = ',iT,' greater than kProfLayer = ',kProfLayer
            write(kStdWarn,*) 'resetting iT = kProfLayer'
            iT = kProfLayer
        END IF
        IF (iB > kProfLayer) THEN
            write(kStdWarn,*) 'IB = ',iB,' greater than kProfLayer = ',kProfLayer
            write(kStdWarn,*) 'resetting iB = 1'
            iB = 1
        END IF
        kActualJacsT = iT
        kActualJacsB = iB
    END IF
    99 FORMAT(I9)

    IF ((iabs(kJacobOutput) /= 1) .AND. (kJacobOutput /= 0) .AND. &
    (kJacobOutput /= 2))  THEN
        write(kStdErr,*) 'kJacobOutput = ',kJacobOutput
        write(kStdErr,*) 'In *PARAMS, need kJacobOutput =-1,0,+1,+2'
        write(kStdErr,*) 'kJacobOutput = format to output Jacobians'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

    IF ((kFlux < -1) .OR. (kFlux >  6)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kFlux =-1 OR 1,2,3,4,5,6'
        write(kStdErr,*) 'where kFlux = do/do not compute fluxes'
        write(kStdErr,*) 'OR         program needs kFlux =-1,+1,2,3,4,5,6 OR 0'
        write(kStdErr,*) 'where kFlux = do not/do  output NLTE Planck'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
    IF ((abs(kSurfTemp)-1.0) >= 1e-5) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kSurfTemp = +/-1.0'
        write(kStdErr,*) 'where kSurfTemp tells the program how to use'
        write(kStdErr,*) 'the surftemperatures in *RADNCE'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

!!!kRTP = -6 : read LBLRTM       LAYERS profile; set atm from namelist
!!!kRTP = -5 : read LBLRTM       LEVELS profile; set atm from namelist
!!!kRTP = -10 : read LEVELS            profile; set atm from namelist
!!!kRTP = -2  : read GENLN2 style LAYERS profile; set atm from namelist
!!!kRTP = -1  : read old    style LAYERS profile; set atm from namelist
!!!kRTP =  0  : read RTP   style kLAYERS profile; set atm from namelist
!!!kRTP = +1  : read RTP   style kLAYERS profile; set atm from RTP file
!!!kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
    IF ((kRTP /= -5) .AND. (kRTP /= -6) .AND. (kRTP /= +2) .AND. &
    (kRTP /= -10) .AND. ((kRTP < -2) .OR. (kRTP > 1))) THEN
        write(kStdErr,*) 'Need to set RTP = -10,-6,-5,-2,-1,0,+1,+2'
        write(kStdErr,*) 'Please reset kRTP and retry'
        CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == 1)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from RTP file'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == 2)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from JPL/NOAA input'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. ((kRTP == -5) .OR. (kRTP == -6))) THEN
        write(kStdErr,*) 'Will read surface temperature info from LBLRTM file'
        write(kStdErr,*) 'and ask kCARTA to add on raTSurf offset from nm_radnces!!!'
    !        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
    !        CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == -10)) THEN
        write(kStdErr,*) 'Cannot read surface temperature info from LEVELS TXT file'
        write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
        CALL DoSTOP
    END IF
     
    IF ((kTempJac < -2) .OR. (kTempJac > 0)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kTempJac=-2,-1,0'
        write(kStdErr,*) 'where kTempJac = use Planck or tau or both '
        write(kStdErr,*) 'when doing d/dT'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE CheckParams


!************************************************************************
! this subroutine sets the mixed path effective tempertures, weighted
! according to the mixing table
    SUBROUTINE GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, &
    iNpmix,iCO2)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    DO iI = 1,kMixFilRows
        raMixVertTemp(iI) = 0.0
    END DO

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
! default = -1

    rW = 0.0
    DO iJ = 1,iNumGases
        rW = rW+raaMix(50,iJ)
    END DO
    IF (rW <= 0.00001) THEN
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
    !      ELSEIF ((kGasTemp .EQ. -1) .AND. (iCO2 .GT. 0) .AND. raaMix(50,iCO2) .LE. 0.001) THEN
    ! user wants the CO2 profile to be the temperature profile, but CO2 wgt == 0, so this is an OOPS moment, quick fix
    !        write(kStdErr,*) 'oops in GetMixVertTemp,,kGasTemp,iCO2 = ',kGasTemp,iCO2
    !        DO iI = 1,iNpmix
    !          iL = MP2Lay(iI)
    !          raMixVertTemp(iI) = raaTemp(iL,iCO2)
    !          write(kStdErr,*) 'oops in GetMixVertTemp, iI,raMixVertTemp(iI) = ',iI,raMixVertTemp(iI)
    !        END DO
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
            IF (rW <= 0.00001) THEN
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

    include '../INCLUDE/kcartaparam.f90'

    REAL :: raPressLevels(kProfLayer+1),pProf(kProfLayer)
    INTEGER :: iProfileLayers

    INTEGER :: iI,iJ,iDiff

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
    end SUBROUTINE FindAvgLayerPressure

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
    SUBROUTINE databasestuff(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)


    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

    IF (iLowerOrUpper < 0) THEN
        CALL databasestuff_lower(iLowerOrUpper, &
        raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
    ELSEIF (iLowerOrUpper > 0) THEN
        CALL databasestuff_upper(iLowerOrUpper, &
        raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
    END IF

    RETURN
    end SUBROUTINE databasestuff

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the lower atm
    SUBROUTINE databasestuff_lower(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)


    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airsheightsparam.f90'
    include '../INCLUDE/airslevelsparam.f90'
    include '../INCLUDE/airslevelheightsparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

! local vars
    INTEGER :: iI

    IF (iLowerOrUpper > -1) THEN
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
    end SUBROUTINE databasestuff_lower

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the uuper atm
    SUBROUTINE databasestuff_upper(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airsheights_upperparam.f90'
    include '../INCLUDE/airslevels_upperparam.f90'
    include '../INCLUDE/airslevelheights_upperparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

! local vars
    INTEGER :: iI

    IF (iLowerOrUpper < +1) THEN
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
    end SUBROUTINE databasestuff_upper

!************************************************************************
! this subroutine finds the partial pressure of the layer
! especially useful for GasID = 1
    SUBROUTINE PPThruLayers( &
    iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press, &
    raDataBaseThickness,raSortPressLevels,raSortPressHeights, &
    raPressLevels,raRTemp,r)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
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

    CALL databasestuff(iLowerOrUpper, &
    DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

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

!      print *,iI,' : ',rP,iJ,DATABASELEV(iJ),'; ',rPp1,iJp1,DATABASELEV(iJp1)

!      print *,raR100Press(iJ),raR100PartPress(iJ),raR100Temp(iJ),
!     $        raR100Amt(iJ),raDataBaseThickness(iJ)
!      r = raR100Amt(iJ)*kAvog/(raDataBaseThickness(iJ) * 100) !!!molecules/cm3
!      print *, raR100PartPress(iJ)*kAtm2mb*100,
!     $       r*kBoltzmann*raR100Temp(iJ)*1e6,
!     $       raR100PartPress(iJ)*kAtm2mb*100/(r*kBoltzmann*raR100Temp(iJ)*1e6)

    IF ((iJp1 - iJ) <= 0) THEN
    ! something wrong; we need lower DATABASE level < upper DATABASE level
        write(kStdErr,*) 'doing water amount ref partial pressure : '
        write(kStdErr,*) 'Layer iI : found iJ = iJp1 ',iI,iJ
        CALL DoStop
    END IF

    IF (iJ <= 0) THEN
        write(kStdErr,*) '>>layer number ',iI,'gas ID ',iGasID
        Call DoStopMesg(' >>oops iJ <= 0 in PPThruLayers, so will have problems in next line$')
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
        rDeltaH = (DATABASELEV(iJ+1)-rP)/ &
        (DATABASELEV(iJ+1)-DATABASELEV(iJ))
        rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
        rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

    !!! add on contribution from upper layer
        rDeltaH = (rPp1-DATABASELEV(iJ+1))/ &
        (DATABASELEV(iJ+2)-DATABASELEV(iJ+1))
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
        rDeltaH = (DATABASELEV(iJ+1)-rP)/ &
        (DATABASELEV(iJ+1)-DATABASELEV(iJ))
        rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
        rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac

    !!! add on contributions from intermediate layers
        DO iL = iJ+1,iJp1-2
            rDeltaH = 1.0
            rNum = rNum + raR100Amt(iL)*rDeltaH                   !!kmol/cm2
            rDenom = rDenom + raDataBaseThickness(iL)
        END DO

    !!! add on contribution from highest layer
        rDeltaH = (rPp1-DATABASELEV(iJp1-1))/ &
        (DATABASELEV(iJp1)-DATABASELEV(iJp1-1))
        rNum1   = raR100Amt(iJp1-1)*rDeltaH                         !! kmol/cm2
        rDenom1 = raDataBaseThickness(iJp1-1)*rDeltaH               !! frac

        rNum   = rNum + rNum1
        rDenom = rDenom + rDenom1
        r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
        r = r*1e6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
        r = r/kAtm2mb                                          !!atm
    END IF

!      print *,iI,'(',iCase,')',r0,r,r/r0,raRTemp(iI)
!      print *,'------------------------------'
!      print *,' '

    RETURN
    end SUBROUTINE PPThruLayers

!************************************************************************
! this function finds the pressure layer at which rPressStart is within,
! as well as the fraction of the layer that it occupies
    REAL FUNCTION FindBottomTemp(rP,raProfileTemp, &
    raPressLevels,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raPressLevels = actual pressure levels that come out of kLAYERS
! raProfileTemp = actual profile temp
! rP            = pressure at which we want the temperature
    real :: rP,raProfileTemp(kProfLayer),raPressLevels(kProfLayer+1)
    integer :: iProfileLayers

    integer :: iFound,i1,i2,i3,iLowest,iJ
    real :: rP1,rP2,T1,T2
    real :: raP(3),raT(3),Y2A(3),rT,raLogP(3)
    real :: yp1,ypn,work(3)
    INTEGER :: iLog,iSpline

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

    IF (rP >= raPressLevels(iLowest)) THEN
    ! his is WHOLE of the bottom layer
        i1=iLowest
    ELSE IF (rP <= raPressLevels(kProfLayer+1)) THEN
    ! his is ludicrous
        write(kStdErr,*) rP,raPressLevels(kProfLayer+1)
        write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
        CALL DoStop

    ELSE
    ! irst find the AIRS layer within which it lies
        iFound=-1
        i1 = iLowest
        i2 = iLowest+1
        10 CONTINUE
        IF ((rP <= raPressLevels(i1)) .AND. (rP > raPressLevels(i2))) THEN
            iFound=1
        END IF
        IF ((iFound < 0) .AND. (i1 < kProfLayer)) THEN
            i1=i1+1
            i2=i2+1
            GO TO 10
        END IF
        IF ((iFound < 0)) THEN
            IF (abs(rP-raPressLevels(kProfLayer+1)) <= delta) THEN
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

    IF ((i1 > kProfLayer) .OR. (i1 < iLowest)) THEN
        write(kStdErr,*) 'sorry : cannot find surface temp for '
        write(kStdErr,*) 'layers outside ',iLowest,' and ',kProfLayer
        write(kStdErr,*) 'Allowed Pressure ranges are from : ', &
        raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
        write(kStdErr,*) 'Surface Pressure is ',rP,' mb'
        call DoStop
    END IF
              
! ow find the temperature
    IF (i1 == iLowest) THEN          !do linear interp
        i1 = iLowest
        i2 = iLowest+1
        i3 = iLowest+2
        rP1 = (raPressLevels(i2)-raPressLevels(i1))/ &
        log(raPressLevels(i2)/raPressLevels(i1))
        rP2 = (raPressLevels(i3)-raPressLevels(i2))/ &
        log(raPressLevels(i3)/raPressLevels(i2))
        T1 = raProfileTemp(i1)
        T2 = raProfileTemp(i2)
        IF (iLog == -1) THEN
            rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)           !!linear in P
        ELSE
            rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1)) !!log(P)
        END IF

    ELSEIF (i1 >= (kProfLayer-1)) THEN          !do linear interp
        rP1 = (raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/ &
        log(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
        rP2 = (raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/ &
        log(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
        T1 = raProfileTemp(kProfLayer-1)
        T2 = raProfileTemp(kProfLayer)
        IF (iLog == -1) THEN
            rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)            !!linear in P
        ELSE
            rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1))  !!log(P)
        END IF
    ELSE          !do spline ... note that the pressures have to
    ! e in ascENDing order for good interpolation
        rP1 = (raPressLevels(i1)-raPressLevels(i1-1))/ &
        log(raPressLevels(i1)/raPressLevels(i1-1))
        raP(3) = rP1
        rP1 = (raPressLevels(i1+1)-raPressLevels(i1))/ &
        log(raPressLevels(i1+1)/raPressLevels(i1))
        raP(2) = rP1
        rP1 = (raPressLevels(i1+2)-raPressLevels(i1+1))/ &
        log(raPressLevels(i1+2)/raPressLevels(i1+1))
        raP(1) = rP1
        IF (iLog == +1) THEN
            DO iJ = 1,3
                raLogP(iJ) = log(raP(iJ))
            END DO
        END IF

        raT(3) = raProfileTemp(i1-1)
        raT(2) = raProfileTemp(i1)
        raT(1) = raProfileTemp(i1+1)

        yp1=1.0e30
        ypn=1.0e30
        IF (iSpline == +1) THEN
            IF (iLog == +1) THEN
                CALL rspl1(raLogP,raT,3,log(rP),rT,1)
            ELSE
                CALL rspl1(raP,raT,3,rP,rT,1)
            END IF
        ELSEIF (iSpline == -1) THEN
            IF (iLog == +1) THEN
                CALL rlinear1(raP,raT,3,rP,rT,1)
            ELSE
                CALL rlinear1(raLogP,raT,3,log(rP),rT,1)
            END IF
        END IF
    END IF

    FindBottomTemp = rT

    RETURN
    end FUNCTION FindBottomTemp

!************************************************************************
! this is called by kcartamain and kcartabasic, to set up the profiles
    SUBROUTINE Set_Ref_Current_Profs( &
    iJax,rDerivTemp,rDerivAmt, &
    iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress, &
    raaAmt,raaTemp,raaPress,raaPartPress, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    raNumberDensity,pProfNLTE,raMixVertTemp)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iJax                   !! when testing jacobians
    REAL :: rDerivTemp,rDerivAmt      !! when testing jacobians
    INTEGER :: iGas,iaGases(kMaxGas)  !! gasID stored in order they were read in
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
    INTEGER :: iInt

! get the reference profile for the current gas if GAS ID <= kGasXsecHi
    IF ((iaGases(iGas) <= kGasXsecHi) .OR. &
    (iaGases(iGas) == kNewGasHi+1)) THEN
        CALL SetReference(raRAmt,raRTemp,raRPress,raRPartPress, &
        raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
    END IF

! acob
! get actual profiles for the current gas
!c these are set in kcartamain, kcartabasic
!c          rDerivTemp = 0.01
!c          rDerivTemp = 0.1
!c          rDerivTemp = 0.5
!c          rDerivTemp = 1.0
!c          rDerivAmt  = 0.01
!c          rDerivAmt  = 0.1

!     print *,'should I use pProf or profNLTE in NLTE routines???'
!     print *,'looks like i need to use pprofNLTE for Dave Edwards profiles
!     print *,'but can get away with pprof using klayers stuff'
!     print *, ie raPressLevels is consistent with pProf (avg layer press)
!             but it comes from summing partial pressures of MANY gases
!     print *, while for Dave Edwards tests, he only uses one gas (CO2)
!             and so mebbe the pressure levels are not consistent with pProf
!            print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
!     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)

!      Do iInt = 1,80
!        print *,iInt,raaTemp(90,iInt)
!      end do
!      call dostopmesg('stop in Set_Ref_Current_Profs$')
          
    DO iInt=1,kProfLayer
        raTAmt(iInt)          = raaAmt(iInt,iGas)
        raTTemp(iInt)         = raaTemp(iInt,iGas)
        raTPress(iInt)        = raaPress(iInt,iGas)
        raTPartPress(iInt)    = raaPartPress(iInt,iGas)
    ! compute particle number density in number/cm3 (N/V = P/kT)
        raNumberDensity(iInt) = raTPress(iInt)*kAtm2mb*100.0/(kBoltzmann * &
        raTTemp(iInt))*1e-6
    ! NLTE avg layer press in lower atm = kCARTA profile
        pProfNLTE(iInt) = raTPress(iInt)*kAtm2mb

    !       print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
    !     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !c  this stuff is to test the jacobians
    !          raTTemp(iInt)=rDerivTemp+raTTemp(iInt)
    !          raMixVertTemp(iInt)=rDerivTemp+raTTemp(iInt)
    !          raMixVertTemp(iInt)=raTTemp(iInt)

    !c jacob test temperature jacobian ... old
    !          IF (iInt .EQ. iJax) THEN
    !            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
    !            raMixVertTemp(iInt)=raTTemp(iInt)
    !          END IF
    !          IF ((iInt .EQ. iJax) .AND. (iGas .EQ. 1)) THEN
    !            rDummy        = raMixVertTemp(iInt)
    !            raMixVertTemp(iInt) = rDummy+rDerivTemp
    ! c          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
    ! c          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
    !            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
    !          END IF

    !c jacob test column gas  jacobian
    !          IF (iGas .EQ. 1) THEN
    !            rDummy        = raMixVertTemp(iInt)
    !            raMixVertTemp(iInt) = rDummy+rDerivTemp
    ! c          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
    ! c          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
    !            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
    !          END IF

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !c jacob test jacobian ... new, could also be column
    !          iJax = -1
    !          iJax = 65
    !          rDerivTemp = 1.0
    !          IF ((iInt .EQ. iJax) .OR. (iJax .EQ. -1)) THEN
    !            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
    !            raMixVertTemp(iInt)=raTTemp(iInt)
    !            print *,iJax,raaTemp(iInt,iGas),rDerivTemp,raMixVertTemp(iInt)
    !          END IF

    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !c jacob test amount jacobian
    !         iJax = 65
    !         IF ((iInt .EQ. iJax).AND.(iaGases(iGas) .EQ. 2)) THEN
    !           raTAmt(iInt)=raaAmt(iInt,iGas)*(1.0+rDerivAmt)
    !           print *,'iJax, dq = ',iJax,raaAmt(iInt,iGas)*rDerivAmt
    !         END IF
    ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    ! Matlab code to plot the jacs
    ! suppose you have a 55 layer atm (from 46 to 100) and a 1 k temperature perturbation in layer 65
    !   (which is layer 20 of the Radiationg Atmosphere)
    ! [d0,w] = readkcstd('junk0.dat');
    ! [jac,w] = readkcjac('junk.jac');
    ! [dt,w] = readkcstd('junk.dat');
    ! [m,n] = size(jac); iL = (n-4)/3
    
    ! plot(w,(dt-d0)/1,'b.-',w,jac(:,20+55),'r')
    
    ! [fc,qc] = quickconvolve(w,jac,0.25,0.25);
    ! pcolor(fc,1:iL,qc(:,(1:iL)+iL*0)'); shading flat; colorbar  %% Q jacobian
    ! pcolor(fc,1:iL,qc(:,(1:iL)+iL*1)'); shading flat; colorbar  %% T jacobian
    ! pcolor(fc,1:iL,qc(:,(1:iL)+iL*2)'); shading flat; colorbar  %% WGT fcn

    END DO

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

    include '../INCLUDE/kcartaparam.f90'

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
        rX = rX + raaPartPress(iI,2)/raaPress(iI,2)
    END DO
    rX = rX * 1.0e6/iCount
    rCO2MixRatio = rX
    write(kStdWarn,*)'avg rCO2Mix = ',rX,' ppmv'

! >>>>>>>>>>>>>>>>>>>>>>>>>
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
        IF (iI == iLowestLay) rF = rFracBot
        iCount = iCount + 1
        rX = rX + raaAmt(iI,3) * rF
    END DO
    rX = rX * cfac * (kAvog) !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
    write(kStdWarn,*)' column amount ozone = ',rX,' du'

! >>>>>>>>>>>>>>>>>>>>>>>>>
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
    end SUBROUTINE compute_co2_mixratio

!************************************************************************
! this duplicates clear sky atmospheres!
    SUBROUTINE duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
    iaSetEms,raaaSetEmissivity,raSetEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rJunk1,rJunk2
          
! local var
    INTEGER :: iX,iY
    REAL :: rX

! first find out how many raAtmLoop the user did set
    iX = 1
    DO WHILE ((iX <= kMaxAtm) .AND. (raAtmLoop(iX) >= 0.0))
        iX = iX + 1
        IF (iX > kMaxAtm) GOTO 10
    END DO
    10 CONTINUE
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

! now set the param you need to set
    IF (iAtmLoop == 1) THEN
        write(kStdWarn,*) '  Resetting raPressStart for looping'
        write(kStdErr,*)  '  Resetting raPressStart for looping'
        IF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) <= 0)) THEN
            write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
            write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
        ELSEIF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) > 0)) THEN
            write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
            write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
        END IF
        IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
            write(kStdWarn,*) '  ---> warning : reset TOA press for uplook instr is odd'
            write(kStdErr,*)  '  ---> warning : reset TOA press for uplook instr is odd'
            CALL DoStop
        END IF
        DO iX = 1,iNatm
            raPressStart(iX) = raAtmLoop(iX)
            raaPrBdry(iX,1)  = raAtmLoop(iX)
        END DO

    ELSEIF (iAtmLoop == 2) THEN
        write(kStdWarn,*) '  Resetting raPressStop for looping'
        write(kStdErr,*)  '  Resetting raPressStop for looping'
        IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
            write(kStdWarn,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
            write(kStdErr,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
        END IF
        DO iX = 1,iNatm
            raPressStop(iX) = raAtmLoop(iX)
            raaPrBdry(iX,2)  = raAtmLoop(iX)
        END DO

    ELSEIF (iAtmLoop == 3) THEN
        write(kStdWarn,*) '  Resetting raSatZen for looping (so need to compute scanang)'
        write(kStdErr,*)  '  Resetting raSatZen for looping (so need to compute scanang)'
        IF (iaLimb(1) > 0) THEN
            write(kStdErr,*) 'Atm 1 set up for Limb sounding'
            write(kStdErr,*) '  so cannot willy nilly reset scanang'
            write(kStdErr,*) 'Go and reset raStartPress instead'
            CALL DoStop
        ELSE
            write(kStdWarn,*) '  changing user input SatZen (angle at gnd)  to       Instr ScanAng '
            write(kStdWarn,*) '  raAtmLoop(iX) --> raSatAngle(iX)     GNDsecant  --> SATELLITEsecant'

            IF (rSatHeightCom < 0) THEN
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
        111 FORMAT('   ',F10.5,' ---> ',F10.5,'   +++   ',F10.5,' ---> ',F10.5)
              
    ELSEIF (iAtmLoop == 4) THEN
        write(kStdWarn,*) '  Resetting raSolZen for looping'
        write(kStdErr,*)  '  Resetting raSolZen for looping'
        DO iX = 1,iNatm
            rakSolarAngle(iX) = raAtmLoop(iX)
        END DO

    ELSEIF (iAtmLoop == 5) THEN
        write(kStdWarn,*) '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
        write(kStdErr,*)  '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
        DO iX = 1,iNatm
            DO iY = 1,kEmsRegions
                raaaSetEmissivity(iX,iY,2) = raaaSetEmissivity(iX,iY,2) + raAtmLoop(iX)
                raaaSetSolarRefl(iX,iY,2)  = (1-raaaSetEmissivity(iX,iY,2))/kPi
            END DO
        END DO

    ELSEIF (iAtmLoop == 10) THEN
        write(kStdWarn,*) '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'
        write(kStdErr,*)  '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'

    ELSEIF (iAtmLoop == 100) THEN
        write(kStdWarn,*) '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'
        write(kStdErr,*)  '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'

    ELSE
        write(kStdErr,*) 'Dont know what to do with iAtmLoop = ',iAtmLoop
        Call DoStop
    END IF

    IF (iAtmLoop <= 2) THEN
        CALL Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer, &
        iProfileLayers,iaMPSetForRad, &
        raSatHeight,raSatAngle,raPressStart,raPressStop, &
        raFracTop,raFracBot,raPressLevels,raLayerHeight, &
        iakSolar,rakSolarAngle)
    END IF

    RETURN
    end SUBROUTINE duplicate_clearsky_atm

! ************************************************************************
! this subroutine resets iaaRadLayer and/or raSatAngle, if the Start or Stop
! Pressures have changed
    SUBROUTINE Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer, &
    iProfileLayers,iaMPSetForRad, &
    raSatHeight,raSatAngle,raPressStart,raPressStop, &
    raFracTop,raFracBot,raPressLevels,raLayerHeight, &
    iakSolar,rakSolarAngle)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iNatm,iProfileLayers,iaMPSetForRad(kMaxAtm),iakSolar(kMaxAtm)
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm),raLayerHeight(kProfLayer)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm),rakSolarAngle(kMaxAtm)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2),raPressLevels(kProfLayer+1)

    INTEGER :: iC,iX,iStart,iStop,iNlay,iDirection,iInt

    DO iC = 1,iNAtm
        CALL StartStopMP(iaMPSetForRad(iC),raPressStart(iC),raPressStop(iC),iC, &
        raPressLevels,iProfileLayers, &
        raFracTop,raFracBot,raaPrBdry,iStart,iStop)

        IF (iStop >= iStart) THEN
            iNlay = (iStop-iStart+1)
            iDirection = +1                           !down look instr
        ELSE IF (iStop <= iStart) THEN
            iNlay = (iStart-iStop+1)
            iDirection = -1                           !up look instr
        END IF
        IF (iNLay > kProfLayer) THEN
            write(kStdErr,*)'Error for atm # ',iC
            write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
            CALL DoSTOP
        END IF
        iaNumlayer(iC) = iNlay

        DO iInt = 1,iNlay
            iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
        END DO

        write(kStdWarn,*) ' Atm#, Press Start/Stop, iStart,iStop, Nlay = ', &
        iC,raPressStart(iC),raPressStop(iC),iStart,iStop,iNlay

    END DO

    IF (iaLimb(1) > 0) THEN
    !! this is limb sounder, do the angle etc
        DO iC = 1,iNatm
            raSatAngle(iC) = LimbViewScanAng(iC,raPressStart,raSatHeight, &
            iaaRadLayer,raPressLevels,raLayerHeight)
            IF (iaKsolar(iC) >= 0) THEN
                rakSolarAngle(iC) = raSatAngle(iC)  !! this is scanang at TOA instr
                rakSolarAngle(iC) = 89.9            !! this is sol zenith at "surface"
            END IF
        END DO
    END IF
          
    RETURN
    end SUBROUTINE Reset_IaaRadLayer
! ************************************************************************
! this duplicates cloud sky 2slab atmospheres!
    SUBROUTINE duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop, &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
    iaSetEms,raaaSetEmissivity,raSetEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm, &
    cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac
    INTEGER :: ctype1,ctype2
    REAL :: cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

! local var
    INTEGER :: iX,iY,iDebug
    REAL :: rX

    iDebug = +1
    iDebug = -1
    IF (iDebug > 0) THEN

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
        IF (iCldProfile > 0) THEN
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

    IF (iCldProfile > 0) THEN
        write(kStdErr,*) 'Ooops can only duplicate clouds slabs, not profiles'
        CALL DoStop
    END IF

    IF (iNclouds == 1) THEN
        write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
    !!! just duplicate the clear fields
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    ELSEIF ((iNclouds > 2) .OR. (iNclouds <= 0)) THEN
        write(kStdErr,*) 'iNclouds = ',iNclouds ,' huh?? cant duplicate this !!!'
        CALL DoStop
    ELSEIF (iNclouds == 2) THEN
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
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
    END IF

    iDebug = -1
    IF (iDebug > 0) THEN
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
    end SUBROUTINE duplicate_cloudsky2slabs_atm

!************************************************************************
! this duplicates cloud sky 100slab atmospheres!
    SUBROUTINE duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop, &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
    iaSetEms,raaaSetEmissivity,raSetEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm, &
    cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac
    INTEGER :: ctype1,ctype2
    REAL :: cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

! local var
    INTEGER :: iX,iY,iDebug
    REAL :: rX

    iDebug = +1
    iDebug = -1
    IF (iDebug > 0) THEN

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
        IF (iCldProfile > 0) THEN
            print*,'iCldProfile'
            print *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
            print *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
        END IF
    END IF

!************************************************************************
!!! now have to update things
!!! if there are originally 2 clouds then
!!!   just do one 100 layer ice/water cloud and one clear calc

    IF (iCldProfile < 0) THEN
        write(kStdErr,*) 'Ooops can only duplicate 100 layer cloud profiles, not slabs'
        CALL DoStop
    END IF

    write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
!!! just duplicate the clear fields
    CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
    iaSetEms,raaaSetEmissivity,raSetEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    RETURN
    end SUBROUTINE duplicate_cloudsky100slabs_atm

!************************************************************************
! funtion to estimate height (in m), given pressure (in mb)
!  based on US STD atm
    REAL FUNCTION p2h(p)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airsheightsparam.f90'
    include '../INCLUDE/airslevelsparam.f90'

    INTEGER :: iI
    REAL :: pavg,rH,raY2P(kMaxLayer),p,logpavg(kMaxLayer),raHgt(kMaxLayer)
         
    DO iI = 1,100
        pavg = (DATABASELEV(iI+1)-DATABASELEV(iI))/log(DATABASELEV(iI+1)/DATABASELEV(iI))
        logpavg(kMaxLayer-iI+1) = log(pavg)      !! need this to be smallest to largest
        raHgt(kMaxLayer-iI+1)   = DatabaseHEIGHT(iI)
    END DO

!      print *,(logpavg(iI),iI=1,kMaxLayer)
!      print *,(raHgt(iI),iI=1,kMaxLayer)

    IF (p >= DATABASELEV(1)) THEN
        rH = DatabaseHEIGHT(1)
    ELSEIF (p <= DATABASELEV(kMaxLayer+1)) THEN
        rH = DatabaseHEIGHT(kMaxLayer)
    ELSE
        CALL rlinear1(logpavg,raHgt,kMaxLayer,log(p),rH,1)
        CALL rspl1(logpavg,raHgt,kMaxLayer,log(p),rH,1)
    END IF

    p2h = rH

    RETURN
    end FUNCTION p2h

!************************************************************************
! ----------------------------------------------------------------
!https://pages.mtu.edu/~shene/COURSES/cs201/NOTES/chap02/datetime.f90
!  This program uses DATE_AND_TIME() to retrieve the system date
!  and the system time.  Then, it converts the date and time
!  information to a readable format.  This program demonstrates
!  the use of concatenation operator // and substring
! ----------------------------------------------------------------

   SUBROUTINE DateTime(callname)
   IMPLICIT   NONE
   INCLUDE '../INCLUDE/kcartaparam.f90'
   
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
   END SUBROUTINE DateTime
!************************************************************************

END MODULE kcartamisc
