! Copyright 2004
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kreadVTprofiles

USE basic_common
USE spline_and_sort_and_common
USE kcoeffSPL
USE n_gas_wt_spectra
USE n_nonlte_common

IMPLICIT NONE

CONTAINS

!************************************************************************
!        these subroutines read in the upper atm profiles
!************************************************************************
! this routine reads in the mixing ratios from GENLN2 files
! sample files are stored in /home/sergio/KCARTADATA/General/NLTE
! assumes the files contain info only for CO2 : <hgt press temp ppmv idunno>
! data for iNumVibLevels is read in, from GND to 0.005 mb to 0.00005 mb
! eventually we can use the data from 0.005 mb to 0.00005 mb for the upper atm

!        !! read in alt/press/temp/upper mix ratio from caaUpperMixRatio
!        !! this is DIFFERENT from the UA NLTE compressed database!

    SUBROUTINE MixRatio(iGasID,rCO2mult,iLTEIn,caaUpperMixRatio, &
    raUpper_Pres,raUpper_MixRatio,iNumVibLevels, &
    raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
    raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input vars
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    INTEGER :: iGasID, iLTEIn
    REAL :: rCO2mult
! output vars
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer)
    INTEGER :: iNumVibLevels
! this is if we want to see what a std US profile looks like about 0.005 mb
! it assumes the lower atm has CO2 ~ 385 ppmv
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! local vars
    CHARACTER(80) :: caFname,caStr
    INTEGER :: iIOUN,iI,iErr
    REAL :: rH,rP,rT,rPPMV,rX,rCO2Mix,raTemp(2*kProfLayer),rMax
    INTEGER :: iDefault,iCurrent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rCO2Mix = 370            !!!as of 2005, this is what was used to
!!!  make the kCARTA database
!!!however, see for2mat2for_CO2_370_385.m in
!!! KCARTA/UTILITY, as the database can easily
!!! be scaled to make eg 385 ppmv the default
    rCO2Mix = kCO2ppmv * rCO2mult  !!! this is how much the input profile
!!! wants CO2 to be enhanced by, relative
!!! to the background

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    caFName = caaUpperMixRatio(iLteIn)

    GOTO 777

    100 write(kStdErr,*) 'Error reading UpperMixRatio info : filename = '
    WRITE(kStdErr,1070) iErr, caFName
    CALL DoStop

    777 CONTINUE
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine upper_mix_ratio, '
        write (kStdErr,*) 'error reading file  ... '
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    END IF
    kTempUnitOpen = +1

    rMax = -1.0
    iI = 0

    555 CONTINUE
    read(iIOUN,123) caStr
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSE
        READ (caStr,*) rH,rP,rT,rPPMV,rX
        iI = iI + 1
        raUpper_Pres(iI) = rP
        raUpper_MixRatio(iI) = rPPMV
        IF (rPPMV > rMax) THEN
            rMax = rPPMV
        END IF
    END IF

    20 CONTINUE
    READ(iIOUN,123,END=199) caStr
    READ (caStr,*) rH,rP,rT,rPPMV,rX
    iI = iI + 1
    raUpper_Pres(iI) = rP
    raUpper_MixRatio(iI) = rPPMV
    IF (rPPMV > rMax) THEN
        rMax = rPPMV
    END IF
    GOTO 20

    199 CONTINUE
    close (iIOUN)
    kTempUnitOpen = -1

    iNumVibLevels = iI

    iCurrent = +1   !! if we are using eg 385 ppmv lower atm database, make sure
!! the UA CO2 mixing ratios read in from caaUpperMixRatio
!! in the namelist file, are approx this value. This is the
!! default case eg when looking at current AIRS/IASI obs
    iCurrent = -1   !! whatever lower atm CO2 database we are using, assume the
!! mixing ratios in caaUpperMixRatio are correct, and do not
!! adjust them. This is the case when eg developing NLTE
!! compressed database
    iDefault = +1
    IF  (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        IF (iDefault /= iCurrent) THEN
            write(kStdErr,*) 'In MixRatio, iDefault,iCurrent are unequal',iDefault,iCurrent
            write(kStdErr,*) 'kCARTA will use CO2 MixRatios specified in caaUpperMixRatio'
        END IF
    END IF

!! this is mainly for testing purposes
!! eventually may want to move this to pre_defined.param or kcartaparam.f90
    CALL GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std, &
    raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    IF (iCurrent == +1) THEN
        DO iI = 1,iNumVibLevels
        !!! rescale stuff to be more current, as scaled by
        !!! kCO2ppmv/mixtable rCO2Mult
            raUpper_MixRatio(iI) = raUpper_MixRatio(iI) * rCO2Mix/rMax
            write(kStdWarn,*) 'new MR',raUpper_Pres(iI),raUpper_MixRatio(iI),rCO2Mix/rMax
        END DO
    ELSEIF (iCurrent == -1) THEN
        IF (abs(rCO2Mult - 1.0) > 1.0e-05) THEN
            DO iI = 1,iNumVibLevels
            !!! rescale the stuff to that specified in mixing table
                raUpper_MixRatio(iI) = raUpper_MixRatio(iI) * rCO2mult
                write(kStdWarn,*) 'new MR',raUpper_Pres(iI),raUpper_MixRatio(iI),rCO2Mix/rMax
            END DO
        END IF
    END IF

! now see if we need to swap values (from small to large)
    IF (raUpper_Pres(1) > raUpper_Pres(iNumVibLevels)) THEN
        DO iI = 1,iNumVibLevels
            raTemp(iI) = raUpper_Pres(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raUpper_Pres(iI) = raTemp(iI)
        END DO
        DO iI = 1,iNumVibLevels
            raTemp(iI) = raUpper_MixRatio(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raUpper_MixRatio(iI) = raTemp(iI)
        END DO
    END IF

    123 FORMAT(A80)
    1070 FORMAT('ERROR! number ',I5,' opening MixRatio atm profile:',/,A80)

    RETURN
    end SUBROUTINE MixRatio

!************************************************************************
! this subroutine reads in the NLTE profiles
    SUBROUTINE read_nlte_file_forUAinfo( &
    raLayerHeight,raThickness,caFName,pProf,iNLTEStart, &
    daJL,daJU,iaJ_UorL,iGasID,iISO,iAllOrSome,raPresslevels, &
    raLayBot1,raLayTop1,raThick1, &
    raPavg1,raTavg1,raNLTavg1,raQTipsAvg1, &
    iNumVibLevels,iUpper,iLay,dVibCenter,raUpperPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL :: pProf(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    CHARACTER(80) :: caFname
    DOUBLE PRECISION :: daJU(kHITRAN),daJL(kHITRAN)
    INTEGER :: iaJ_UorL(kHITRAN),iGasID,iISO,iAllOrSome,iNLTEStart
! output params
    REAL :: raLayTop1(kNLTEProfLayer),raLayBot1(kNLTEProfLayer)
    REAL :: raThick1(kNLTEProfLayer)
    REAL :: raPavg1(kNLTEProfLayer),raTavg1(kNLTEProfLayer)
    REAL :: raNLTavg1(kNLTEProfLayer),raQTipsAvg1(kNLTEProfLayer)
    INTEGER :: iLay               !! layers [1 --> lay-1] are in LTE
!! layer  [lay ...    ] are in NLTE
    INTEGER :: iUpper             !! how many layers above kPROFLAYER are added
!! ie (ilay) --> (ilay+iUpper) are in NLTE
    INTEGER :: iNumVibLevels      !! how many total levels were read in
!! so (iLay-1 + iUpper) +1 == iNumVibLevels
    DOUBLE PRECISION :: dVibCenter            !!! vib center from GENLN2 file
    REAL :: raUpperPressLevels(kProfLayer+1)

! local vars
    REAL :: raLayTop(kProfLayer),raLayBot(kProfLayer)  !!! AIRS info
    INTEGER :: iNumMixRatioLevs
    REAL :: raTTemp1(kNLTEProfLayer),raTPress1(kNLTEProfLayer)
    REAL :: raNLTETemp1(kNLTEProfLayer)
    REAL :: p1,p2,rAvgHgt,rP,rPp1,rT,rTp1,rSlope
    REAL :: kCARTA_STartNLTE
    REAL :: raQtips1(kNLTEProfLayer)
    REAL :: raQtipsSwap(kProfLayer),raQtips1Swap(kNLTEProfLayer)
    REAL :: A1,A2,A3,A4,A

    INTEGER :: iI,iErr,iIOUN,iJ,iUpper0,iVibs,iVibs0
    INTEGER :: iLay_Plevs,iLay_Pavg,iFound_Pavg,iFound_Plevs

!!!these are the AIRS pressure levels, in mb
!!!and the AIRS heights in meters
    DO iI = 1, kProfLayer
        raLayBot(iI) = raLayerHeight(iI)
        raLayTop(iI) = raLayerHeight(iI) + raThickness(iI)
    END DO

! for the usual AIRS layers
! using Matlab, eqn is hgt = log(pavg) * p1 + p2
! or pavg = exp((y-p2)/p1)
! p1 = -6.924823699945052e+03, p2 = 4.819268370827240e+04
! where pavg is in mb, and hgt is in meters
    p1 = -6.924823699945052E+03
    p2 = +4.819268370827240E+04

!! for the upper AIRS layers
!! if x = log(Press) then y = A1 x3 + A2 x2 + A3 x + A4
    A1 = -6.6542e+01
    A2 = -1.2139e+03
    A3 = -1.2836e+04
    A4 =  4.0587e+04

    write(kStdWarn,*) ' '
          
    CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL, &
    iGasID,iISO,iAllOrSome, &
    iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTETemp1,dVibCenter)

    write(kStdWarn,*) 'UA (p < 0.005 mb) : iGASID,ISO,iLSGQ,iUSGQ,vibcntr = ', &
    iGasID,iISO,int(daJL(1)),int(daJU(1)),dVibCenter

    DO iI = 1,kNLTEProfLayer
        raThick1(iI) = 0.0
    END DO

    write(kStdWarn,*) 'doing the LTE profile in STRATOSphere (< 0.005 mb)'

    kCARTA_STartNLTE = pProf(iNLTEStart)
    write(kStdWarn,*) 'NLTE inside kCARTA database (1100 - 0.005 mb) : '
    write(kStdWarn,*) '  kCARTA NLTE starts above layer',iNLTEStart ,'(pressure <= ',kCARTA_STartNLTE,' mb)'
    write(kStdWarn,*) '  kCARTA database ends ',raPresslevels(kProfLayer+1)

    write(kStdWarn,*) 'NLTE above standard kCARTA database ( < 0.005 mb) : '
    iLay = 1
    iLay_Plevs = 1
    iLay_Pavg = 1
    iFound_Plevs = -1
    iFound_Pavg = -1
    DO iI = 1,iNumVibLevels - 1
        IF (raTPress1(iI+1) >= 0.005) THEN
        !! use usual AIRS layers
            raLayBot1(iI) = log(raTPress1(iI))*p1 + p2
            raLayTop1(iI) = log(raTPress1(iI+1))*p1 + p2
        ELSEIF (raTPress1(iI+1) < 0.005) THEN
        !! use new upper AIRS layers
            A = log(raTPress1(iI))
            raLayBot1(iI) = A1 *(A**3) + A2 * (A**2) + A3 * A + A4
            A = log(raTPress1(iI+1))
            raLayTop1(iI) = A1 *(A**3) + A2 * (A**2) + A3 * A + A4
        END IF
        raThick1(iI)  = abs(raLayTop1(iI) - raLayBot1(iI))

        rP   = raTPress1(iI)
        rPp1 = raTPress1(iI+1)
        raPavg1(iI) = (rP-rPp1)/alog(rP/rPp1)

        rT   = raTTemp1(iI)
        rTp1 = raTTemp1(iI+1)
        rSlope = (rT-rTp1)/(alog(rP) - alog(rPp1))
        raTavg1(iI) = raTTemp1(iI) - rSlope*(alog(rP)-alog(raPavg1(iI)))

        rT   = raNLTETemp1(iI)
        rTp1 = raNLTETemp1(iI+1)
        rSlope = (rT-rTp1)/(alog(rP) - alog(rPp1))
        raNLTavg1(iI) = raNLTETemp1(iI) - rSlope*(alog(rP)-alog(raPavg1(iI)))

        rT   = raQTips1(iI)
        rTp1 = raQTips1(iI+1)
        rSlope = (rT-rTp1)/(alog(rP) - alog(rPp1))
        raQTipsAvg1(iI) = raQTips1(iI) - rSlope*(alog(rP)-alog(raPavg1(iI)))

        IF ((iFound_Plevs < 0) .AND. &
        (raTPress1(iI+1) < raPresslevels(kProfLayer+1))) THEN
        !!!!this is above the kCARTA database
            iLay_Plevs = iI
            iFound_Plevs = +1
        END IF
    !       print *,iI,raThick1(iI),raPavg1(iI),raTavg1(iI),raNLTavg1(iI),
    !     $          raQTipsAvg1(iI)

        IF ((iFound_Pavg < 0) .AND. &
        (raPavg1(iI) < pProf(kProfLayer))) THEN
        !!!!this is above the kLAYERS layers
            iLay_Pavg = iI
            iFound_Pavg= +1
        END IF
    !       print *,iI,raThick1(iI),raPavg1(iI),raTavg1(iI),raNLTavg1(iI),
    !     $          raQTipsAvg1(iI)
    END DO

    iLay = iLay_Pavg
    iUpper = (iNumVibLevels-iLay + 1 ) - 1

    iLay = iLay_Plevs
    iUpper = (iNumVibLevels-iLay + 1 ) - 1

    IF (iLay_Plevs /= iLay_Pavg) THEN
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(kStdWarn,*) 'NLTE profile file ',caFName
        write(kStdWarn,*) 'has ',iNumVibLevels,' levels of information ...'
        write(kStdWarn,*) 'oops : iLay_Plevs,iLay_Pavg = ',iLay_Plevs,iLay_Pavg
        write(kStdWarn,*)  raTPress1(iLay_Plevs),raTPress1(iLay_Plevs+1)
        write(kStdWarn,*)  raTPress1(iLay_Pavg),raTPress1(iLay_Pavg+1)
        write(kStdWarn,*) 'kLayers may be using layers around 70-80 km'
        write(kStdWarn,*) 'that are too thick!!!!'
        write(kStdWarn,*) ' ---- > using iLay = ',iLay
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    END IF

    iUpper0 = iUpper
    IF (iUpper > kProfLayer) THEN
        iUpper = kProfLayer
        write(kStdWarn,*) 'reset iUpper to ',iUpper
    END IF

    DO iI = 1,iUpper+1
        iJ = iLay+(iI-1)
        raUpperPressLevels(iI) = raTPress1(iJ)
    END DO

    write(kStdWarn,*) 'in read_nlte_file_forUAinfo : '
    write(kStdWarn,*) '  number of levels in NLTE profile = ',iNumVibLevels
    write(kStdWarn,*) '  number of levels above 0.005 mb  = ',iUpper

    RETURN
    end SUBROUTINE read_nlte_file_forUAinfo

!************************************************************************
! this subroutine reads in the NLTE temperatures as well as the info
! for the (t,p,pp,q) parameters for the upper part of the STRATOSphere
    SUBROUTINE read_upperatm_lte_temperature( &
    iGasID,iNLTEStart,iLTEin,iBand,caaNLTETemp, &
    raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
    pProf,raPresslevels,raLayerHeight,raThickness, &
    iUpper,raTTemp,raTAmt,raTPress,raTPartPress, &
    raUpperPressLevels,raUpperThickness)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
! read file from caaNLTETemp(iLTEin)
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! iGASID          tells the gasID (surprise!!!!!!!!!!!!)
! iNLTEStart    tells where the NLTE starts in the kCARTA atmosphere
! pProf,raPresslevels are the AIRS pressure info
! daJL,daJU are the gas lower, upper quantum numbers
    REAL :: pProf(kProfLayer),raPresslevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    INTEGER :: iLTEIn,iBand,iGasID,iNLTEStart
    CHARACTER(80) :: caaNLTETemp(kGasStore)
! these are the input mixing ratios from Dave Edwards
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer)
    INTEGER :: iNumMixRatioLevs
! output params
    INTEGER :: iUpper             !!how many layers above kPROFLAYER are added
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer)
    REAL :: raTPress(kProfLayer),raTPartPress(kProfLayer)
! this tells the pressure levels and layer thicknesses for upper atm
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)

! local variables
    INTEGER :: iISO,iAllOrSome,iNumVibLevels
    REAL :: rMixRatio,MGC,rAvgHgt,raUAMixRatio(kNLTEProfLayer)
    REAL :: raLayTop1(kNLTEProfLayer),raLayBot1(kNLTEProfLayer)
    REAL :: raThick1(kNLTEProfLayer)
    REAL :: raPavg1(kNLTEProfLayer),raTavg1(kNLTEProfLayer), &
    raNLTavg1(kNLTEProfLayer)
    REAL ::    rDummyVibEnergy
    REAL :: raQtips1(kNLTEProfLayer),raQTipsAvg1(kNLTEProfLayer)

    INTEGER :: iI,iErr,iIOUN,iLay,iFound,iJ,iUpper0,iVibs,iVibs0
    DOUBLE PRECISION :: daJU(kHITRAN),daJL(kHITRAN)
    INTEGER :: iaJ_UorL(kHITRAN)
    CHARACTER(80) :: caFname
    DOUBLE PRECISION :: dVibCenter            !!! vib center from GENLN2 file

    CHARACTER(120) :: ca120A,ca120B,ca120C,ca120D
    
    REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer),raUA_refMR(kMaxLayer)
    REAL :: raMixRatioUA_ref(kNLTEProfLayer),raTUA_ref(kNLTEProfLayer)

    iUpper = -1               !!!assume nothing up on high
    MGC = kMGC

!!! read in information supplied by Dave Edwards
!!! he supplies iNumVibLevels pressure,T(lte),T(nonlte) points
    caFName = caaNLTETemp(iLTEin)

!!!! assign iISO,daJU,daJL = dummy values, to use ReadGENLN2_NLTE_Profile
!!!! to read in the kinetic profile
    iISO    = 1                  !!!just a dummy
    daJU(1) = 9 * 1.0d0          !!!just a dummy
    daJL(1) = 1 * 1.0d0          !!!just a dummy
    iAllOrSome = +1              !!!just a dummy

    CALL read_nlte_file_forUAinfo( &
    raLayerHeight,raThickness,caFName,pProf,iNLTEStart, &
    daJL,daJU,iaJ_UorL,iGasID,iISO,iAllOrSome,raPressLevels, &
    raLayBot1,raLayTop1,raThick1, &
    raPavg1,raTavg1,raNLTavg1,raQTipsAvg1, &
    iNumVibLevels,iUpper,iLay,dVibCenter,raUpperPressLevels)

    CALL getUArefprofile(1,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)
    DO iI = 1,kMaxLayer
      raUA_refMR(iI) = raUA_refPP(iI)/raUA_refP(iI)*1.0e6
    END DO

!!! -------------------------------------------------------------------------->>>>>>
    DO iI = 1,kProfLayer
        raTPress(iI)     = 0.0
        raTPartPress(iI) = 0.0
        raTTemp(iI)      = 0.0
        raTAmt(iI)    = 0.0
    END DO

    IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop <= 2)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        write(kStdWarn,*) 'LOWER atm LTE temperatures, gas amounts ....'
        write(kStdWarn,*) 'read in direct from GENLN2 NLTE profile file ',caFName
        write(kStdWarn,*) ' this part is for jollies, not really used ....'
        write(kStdWarn,*) 'iI     pavg   dz   q        T(iI)  '
        write(kStdWarn,*) '       (mb)   (m) (mol/cm2)  (K)'
        write(kStdWarn,*) '---------------------------------------------------'
    END IF

!!!dope, for the regular kCARTA atmosphere, just a check
    DO iI = 1,iLay-1
        iJ = iI
        raUpperThickness(iI) = raThick1(iJ)  !!! in m
        raTPress(iI)     = raPavg1(iJ)   !!! in mb
        raTTemp(iI)      = raTavg1(iJ)
    ! ppmv = (number of gas molecules)/(number of air molecules) * 1e6
    ! pV = NRT ==> N(g) = p(g) V(a)/RT(g),  N(a) = p(a) V(a)/ RT(a)
    ! now V(g) == V(a), T(g) == T(a) and p(total) = p(a) + p(g)
    ! thus ppmv = N(g)/N(a) = p(g)/p(a) = p(g)/(p(total) - p(g)) * 1e6
    ! thus p(g) = (p(total)/(1e6+ppmv)) ppmv
        rAvgHgt = (raLayTop1(iJ)+raLayBot1(iJ))/2.0
        CALL rspl_one(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs,raTPress(iI),rMixRatio)
        raTPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raTPartPress(iI)=raTPress(iI)*rMixRatio/(1.0e6+rMixRatio) !atm
        raTAmt(iI) = raThick1(iJ) * 100.0           !!! in cm
        raTAmt(iI)=raTAmt(iI)*kAtm2mb*100.0*raTPartPress(iI)
        raTAmt(iI) = raTAmt(iI)/1e9/MGC/raTTemp(iI)

        IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1))) THEN
            write(kStdWarn,6543) iI,raTPress(iI)*kAtm2mb,raUpperThickness(iI), &
            raTAmt(iI),raTTemp(iI)
        END IF
    END DO

!!! -------------------------------------------------------------------------->>>>>>
!! now concentrate on upper atm
    DO iI = 1,kProfLayer
        raTPress(iI)     = 0.0
        raTPartPress(iI) = 0.0
        raTTemp(iI)      = 0.0
        raTAmt(iI)    = 0.0
    END DO

    IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'UPPER atm LTE temperatures, gas amounts ....'
        write(kStdWarn,*) 'read in direct from GENLN2 NLTE profile file ',caFName
	write(kStdWarn,*) ' '
	write(kStdWarn,*) ' >>> NOTE if the CO2 mix ratios do not agree, this routine'
	write(kStdWarn,*) ' >>> will modify the input MR so they are same as US STd by end of this subr'
	write(kStdWarn,*) ' '
        ca120B = &
        ' |                    VT profile                                   | US STd UA reference     |           '
        ca120C = &
        ' lay    p (mb)         pp(mb)        T(K)      Q(mol/cm2)  CO2 ppm |  Tstd(K)   CO2Std(ppm)  | CO2Std/CO2'
        ca120D = &
        '---------------------------------------------------------------------------------------------------------'
        write(kStdWarn,120) ca120B
        write(kStdWarn,120) ca120C
        write(kStdWarn,120) ca120D      
    END IF
     
    DO iI = 1,iUpper
        iJ = iLay+(iI-1)
        raUpperThickness(iI) = raThick1(iJ)  !!! in m
        raTPress(iI)     = raPavg1(iJ)   !!! in mb
        raTTemp(iI)      = raTavg1(iJ)
    ! ppmv = (number of gas molecules)/(number of air molecules) * 1e6
    ! pV = NRT ==> N(g) = p(g) V(a)/RT(g),  N(a) = p(a) V(a)/ RT(a)
    ! now V(g) == V(a), T(g) == T(a) and p(total) = p(a) + p(g)
    ! thus ppmv = N(g)/N(a) = p(g)/p(a) = p(g)/(p(total) - p(g)) * 1e6
    ! thus p(g) = (p(total)/(1e6+ppmv)) ppmv
        rAvgHgt = (raLayTop1(iJ)+raLayBot1(iJ))/2.0
        CALL rspl_one(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs,raTPress(iI),rMixRatio)
        raUAMixRatio(iI)     = rMixRatio		      
        raTPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raTPartPress(iI) = raTPress(iI)*rMixRatio/(1.0e6+rMixRatio) !!atm
        raTAmt(iI)    = raThick1(iJ) * 100.0           !!! in cm
        raTAmt(iI)    = raTAmt(iI)*kAtm2mb*100.0*raTPartPress(iI)
        raTAmt(iI)    = raTAmt(iI)/1e9/MGC/raTTemp(iI)

        CALL r_sort_logspl(raUA_refP,raUA_refT, kMaxLayer,raTPress(iI),raTUA_ref(iI))
        CALL r_sort_logspl(raUA_refP,raUA_refMR,kMaxLayer,raTPress(iI),raMixRatioUA_ref(iI))

        IF ((iBand == 1) .AND. &
          (((abs(kLongOrShort) == 2) .AND. (kOuterLoop <= 2)) .OR. &
          (abs(kLongOrShort) <= 1))) THEN
          write(kStdWarn,1080) kProfLayer+iI,raTPress(iI)*kAtm2mb, &
                  raTPartPress(iI)*kAtm2mb,raTTemp(iI),raTamt(iI),rMixRatio, &
                  raTUA_ref(iI),raMixRatioUA_ref(iI),raMixRatioUA_ref(iI)/rMixRatio
        END IF
	raTPartPress(iI) = raTPartPress(iI) * raMixRatioUA_ref(iI)/rMixRatio
	raTAmt(iI)       = raTAmt(iI) * raMixRatioUA_ref(iI)/rMixRatio
	raUAMixRatio(iI) = raUAMixRatio(iI) * raMixRatioUA_ref(iI)/rMixRatio
	
    END DO

    1070 FORMAT('ERROR! number ',I5,' opening upper atm NONLTE profile:',/,A80)
    1080 FORMAT(I4,' ',2(ES12.5,'  '),1(F10.4,'  '),ES12.5,F10.4,'|',2(F10.4,'  '),F10.5)    
    6543 FORMAT(I4,' ',1(E9.4,' '),1(F9.3,' '),1(E9.4,' '),1(F8.3,' '))
    120  FORMAT(A120)    

    RETURN
    end SUBROUTINE read_upperatm_lte_temperature

!************************************************************************
! this subroutine reads in the NONLTE temperatures as well as the info
! for the (t,p,pp,q) parameters for the upper part of the STRATOSphere
! ie for pavg < 0.005 mb (extent of kCARTA database)
    SUBROUTINE read_upperatm_nonlte_temperature( &
    iGasID,iISO,iNLTEStart,iLTEin,iBand,caaNLTETemp, &
    raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
    pProf,raPresslevels,raLayerHeight,raThickness, &
    iUpper,raTTemp,raTAmt,raTPress,raTPartPress, &
    daJL,daJU,iaJ_UorL,raNLTEtemp,raVibQFT, &
    iAllLayersLTE,dVibCenter, &
    raUpperPressLevels,raUpperThickness, &
    raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
    raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
! read file from caaNLTETemp(iLTEin)
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! iGASID        tells the gasID (surprise!!!!!!!!!!!!)
! iNLTEStart    tells where the NLTE starts in the kCARTA atmosphere
! pProf,raPresslevels are the AIRS pressure info
! daJL,daJU are the gas lower, upper quantum numbers
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    REAL :: pProf(kProfLayer),raPresslevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    INTEGER :: iLTEIn,iBand,iGasID,iNLTEStart,iISO,iAllLayersLTE
    CHARACTER(80) :: caaNLTETemp(kGasStore)
! these are the input mixing ratios from Dave Edwards
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer)
    INTEGER :: iNumMixRatioLevs
! these are to do with the US Std UA profiles
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! output params
    INTEGER :: iUpper             !!how many layers above kPROFLAYER are added
    REAL :: raNLTETemp(kProfLayer),raVibQFT(kProfLayer)
    DOUBLE PRECISION :: dVibCenter
    INTEGER :: iaJ_UorL(kHITRAN)
! input/output params
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
    REAL :: raTPartPress(kProfLayer)
! this tells the pressure levels and layer thicknesses for upper atm
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)

! local variables
    REAL :: rMixRatio,MGC,rAvgHgt
    REAL :: raLayTop1(kNLTEProfLayer),raLayBot1(kNLTEProfLayer)
    REAL :: raThick1(kNLTEProfLayer)
    REAL :: raPavg1(kNLTEProfLayer),raTavg1(kNLTEProfLayer), &
    raNLTavg1(kNLTEProfLayer)
    REAL :: rDummyVibEnergy
    REAL :: raQtips1(kNLTEProfLayer),raQTipsAvg1(kNLTEProfLayer)

    INTEGER :: iI,iErr,iIOUN,iLay,iFound,iJ,iUpper0,iVibs,iVibs0,iAllOrSome
    INTEGER :: iNumVibLevels
    REAL :: p,pp,q,t,raUAMixRatio(kNLTEProfLayer)
    CHARACTER(80) :: caFname
    CHARACTER(80) :: caStr
    CHARACTER(3) ::  ca3
    CHARACTER(120) :: ca120A,ca120B,ca120C,ca120D

    REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer),raUA_refMR(kMaxLayer)
    REAL :: raMixRatioUA_ref(kNLTEProfLayer),raTUA_ref(kNLTEProfLayer)
    
    iUpper = -1               !!!assume nothing up on high
    MGC = kMGC
    iAllOrSome = +1

!!! read in information supplied by Dave Edwards
!!! he supplies iNumVibLevels pressure,T(lte),T(nonlte) points
    caFName = caaNLTETemp(iLTEin)

    CALL read_nlte_file_forUAinfo( &
    raLayerHeight,raThickness,caFName,pProf,iNLTEStart, &
    daJL,daJU,iaJ_UorL,iGasID,iISO,iAllOrSome,raPressLevels, &
    raLayBot1,raLayTop1,raThick1, &
    raPavg1,raTavg1,raNLTavg1,raQTipsAvg1, &
    iNumVibLevels,iUpper,iLay,dVibCenter,raUpperPressLevels)

    CALL getUArefprofile(1,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)
    DO iI = 1,kMaxLayer
      raUA_refMR(iI) = raUA_refPP(iI)/raUA_refP(iI)*1.0e6
    END DO
    
    DO iI = 1,kProfLayer
        raTPress(iI)     = 0.0
        raTPartPress(iI) = 0.0
        raTTemp(iI)      = 0.0
        raTamt(iI)       = 0.0
    END DO

    IF ((iBand == 1) .AND. &
      (((abs(kLongOrShort) == 2) .AND. (kOuterLoop <= 2)) .OR. &
      (abs(kLongOrShort) <= 1))) THEN
        ca120A = '      Comparing VT profile info against US STd UA info'
        ca120B = &
        ' |                    VT profile                                   | US STd UA reference     |           '
        ca120C = &
        ' lay    p (mb)         pp(mb)     Tlte(K)      Q(mol/cm2)  CO2 ppm |  Tstd(K)   CO2Std(ppm)  | CO2Std/CO2'
        ca120D = &
        '---------------------------------------------------------------------------------------------------------'
        write(kStdWarn,120) ca120A
        write(kStdWarn,120) ca120B
        write(kStdWarn,120) ca120C
        write(kStdWarn,120) ca120D      
    END IF
    
    DO iI = 1,iUpper
        iJ = iLay+(iI-1)
        raUpperThickness(iI) = raThick1(iJ)  !!! in m
        raTPress(iI)    = raPavg1(iJ)         !!!in mb
        raTTemp(iI)     = raTavg1(iJ)
    ! for sig sig, pipi, deltdelt use NLTE
        IF (iAllLayersLTE == -1) THEN       !!!upper levels in NLTE
            raVibQFT(iI)   = raQtipsAvg1(iJ)
            raNLTETemp(iI) = raNLTAvg1(iJ)
        ELSE                                  !!!all levels in LTE
        ! aVibQFT(iI)   = 1.0               !!!we have read in correct QV
            raNLTETemp(iI) = raTavg1(iJ)
        END IF
    ! ppmv = (number of gas molecules)/(number of air molecules) * 1e6
    ! pV = NRT ==> N(g) = p(g) V(a)/RT(g),  N(a) = p(a) V(a)/ RT(a)
    ! now V(g) == V(a), T(g) == T(a) and p(total) = p(a) + p(g)
    ! thus ppmv = N(g)/N(a) = p(g)/p(a) = p(g)/(p(total) - p(g)) * 1e6
    ! thus p(g) = (p(total)/(1e6+ppmv) ppmv
        rAvgHgt = (raLayTop1(iJ)+raLayBot1(iJ))/2.0
        CALL rspl_one(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs,raTPress(iI),rMixRatio)
        raUAMixRatio(iI) = rMixRatio
        raTPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raTPartPress(iI) = raTPress(iI) * rMixRatio/(1.0e6 + rMixRatio) !!! atm
        raTAmt(iI)       = raThick1(iJ) * 100.0           !!! in cm
        raTAmt(iI)       = raTAmt(iI)*kAtm2mb*100.0*raTPartPress(iI)
        raTAmt(iI)       = raTAmt(iI)/1e9/MGC/raTTemp(iI)
	
        CALL r_sort_logspl(raUA_refP,raUA_refT, kMaxLayer,raTPress(iI),raTUA_ref(iI))
        CALL r_sort_logspl(raUA_refP,raUA_refMR,kMaxLayer,raTPress(iI),raMixRatioUA_ref(iI))

        IF ((iBand == 1) .AND. &
          (((abs(kLongOrShort) == 2) .AND. (kOuterLoop <= 2)) .OR. &
          (abs(kLongOrShort) <= 1))) THEN
          write(kStdWarn,1080) kProfLayer+iI,raTPress(iI)*kAtm2mb, &
                  raTPartPress(iI)*kAtm2mb,raTTemp(iI),raTamt(iI),rMixRatio, &
                  raTUA_ref(iI),raMixRatioUA_ref(iI),raMixRatioUA_ref(iI)/rMixRatio
        END IF
	raTPartPress(iI) = raTPartPress(iI) * raMixRatioUA_ref(iI)/rMixRatio
	raTAmt(iI)       = raTAmt(iI) * raMixRatioUA_ref(iI)/rMixRatio
	raUAMixRatio(iI) = raUAMixRatio(iI) * raMixRatioUA_ref(iI)/rMixRatio
    END DO
    
! ccc this is typical kProfLayer (~ 100) CO2 amounts
! ccc 200 2 0.93681E-05 0.34536E-08   198.16440  0.14706E-09

    IF ((iBand == 1) .AND. &
    (((abs(kLongOrShort) == 2) .AND. (kOuterLoop <= 2)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        ca120A = &
        '       klayers                                     |            vibtemp file'
        ca120B = &
        'iI        pavg       dz          q         T(iI)   |  Tlte       dT      | Tnlte        dT    | MixRatio'
        ca120C = &
        '          mb         m       kmol/cm2  K           |                     |                    |   ppmv'
        ca120D = &
        '--------------------------------------------------------------------------------------------------------'
	write(kStdWarn,*) ' '
        write(kStdWarn,*) 'this is after readjusting VT profile ppmv to be same as that in the reference profile ...'	
        write(kStdWarn,120) ca120A
        write(kStdWarn,120) ca120B
        write(kStdWarn,120) ca120C
        write(kStdWarn,120) ca120D
        DO iI = 1,iUpper
            write(kStdWarn,9876) iI+kProfLayer,raTPress(iI)*kAtm2mb,raUpperThickness(iI),raTamt(iI),raTTemp(iI), &
            raTTemp(iI),00.0, &
            raNLTETemp(iI), raNLTETemp(iI)-raTTemp(iI),raUAMixRatio(iI)
        END DO

    !! compare what we have to US Std UA
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'Comparison between user supplied and US STD kCompressed ',min(iUpper,iUpperStd_Num),' UA layers'
        ca120A = &
        '             CURRENT PROF (N)  LAYERS                        |            US STD UA (interpd) LAYERS'
        ca120B = &
        'iI     pavg          dz          q         T(iI)   MixRatio  |  T(iI)   MixRatio'
        ca120C = &
        '       mb            m        kmol/cm2        K      ppmv    |    K        ppm'
        ca120D = &
        '--------------------------------------------------------------------------------------------------------------------'
        write(kStdWarn,120) ca120A
        write(kStdWarn,120) ca120B
        write(kStdWarn,120) ca120C
        write(kStdWarn,120) ca120D
        DO iI = 1,min(iUpper,iUpperStd_Num)
            write(kStdWarn,4322) iI+kProfLayer,raTPress(iI)*kAtm2mb,raUpperThickness(iI),raTamt(iI), &
            raTTemp(iI),raUAMixRatio(iI), &
            raTUA_ref(iI),raMixRatioUA_ref(iI)
        END DO
    END IF

    1070 FORMAT('ERROR! number ',I5,' opening upper atm NONLTE profile:',/,A80)
    1080 FORMAT(I4,' ',2(ES12.5,'  '),1(F10.4,'  '),ES12.5,F10.4,'|',2(F10.4,'  '),F10.5)
    9876 FORMAT(I4,' ',ES12.5,' ',F9.4,' ',ES12.5,' ',F9.4,' | ',4(F9.4,' '),' | ',F9.4)
    4322 FORMAT(I4,' ',1(ES12.5,' ',F9.4,' ',ES12.5,' ',2(F9.4,' ')),'|',2(F9.4,' '))    
    120  FORMAT(A120)

    RETURN
    end SUBROUTINE read_upperatm_nonlte_temperature

!************************************************************************
!  these subroutines read in NLTE profiles for upper tropo/lower stratos
!************************************************************************
! this subroutine reads in the NONLTE temperatures for the upper part of
! the standard KCARTA atmosphere; files are from Dave Edwards and only have
! pressure level, kinetic temp, NLTE temp  info
! ie for 1100 <= p <= 0.005 mb (within the kCARTA database)
    SUBROUTINE read_nonlte_temperature(iGasID,iISO,iLTEin,iBand, &
    caaNLTETemp, &
    pProf,raPressLevels,raLayerHeight,raThickness,iProfileLayers, &
    raTPress,raTPartPress,raTTemp,raTAmt,daJL,daJU, &
    iaJ_UorL,raLTETemp,raNLTETemp,raVibQFT,iAllLayersLTE,dVibCenter)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
! read file from caaNLTETemp(iLTEin)
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! iGASID          tells the gasID (surprise!!!!!!!!!!!!)
! iNLTEStart    tells where the NLTE starts in the kCARTA atmosphere
! pProf,raPresslevels are the AIRS pressure info
! daJL,daJU are the gas lower, upper quantum numbers
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    REAL :: pProf(kProfLayer),raPresslevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    REAL :: raTTemp(kProfLayer)    !!the kinetic or LTE temp
    REAL :: raTAmt(kProfLayer)     !!actual gas amounts
    REAL :: raTPress(kProfLayer),raTPartPress(kProfLayer) !!actual p,pp
    INTEGER :: iLTEIn,iBand,iGasID,iISO,iAllLayersLTE,iProfileLayers
    CHARACTER(80) :: caaNLTETemp(kGasStore)
! output params are nonLTE temperatures, and Qtips_vib interpolated to profile
    REAL :: raNLTETemp(kProfLayer)   !!!NLTE temperatures
    REAL :: raLTETemp(kProfLayer)    !!!want to send this is from the vib-temp
!!!files as it might be slightly but
!!!importantly different from LAYERS T(z)
    REAL :: raVibQFT(kProfLayer)     !!!modifier to partition fcns
    DOUBLE PRECISION :: dVibCenter   !!!from Dave Edwards files
    INTEGER :: iaJ_UorL(kHITRAN)     !!!are we matching IUSGQ (+1) or ILSGQ (-1)

! local variables
! read the VibTemp profiles(p,LTE,NLTE,qvib) info into these variables
    REAL :: raPressVT(kNLTEProfLayer),  raKineticVT(kNLTEProfLayer)
    REAL :: raNLTETempVT(kNLTEProfLayer),raQtipsVT(kNLTEProfLayer)
    REAL :: raMixRatio(kNLTEProfLayer)

    INTEGER :: iI,iErr,iIOUN,iLay,iFound,iJ,iVibs,iVibs0,iStart
    INTEGER :: i1,i2,iNumVibLevels,iLogOrLinear,iLP,iDefault
    REAL :: rP,rPP,rQ,rT,rFac
    CHARACTER(80) :: caFname
    CHARACTER(100) :: ca1,ca2
    REAL :: raNLTE_STD(kProfLayer),raLTE_STD(kProfLayer)

!!! read in information supplied by user
!!! [supplies iNumVibLevels pressure,T(lte),T(nonlte) points]
    caFName = caaNLTETemp(iLTEin)

    CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL,iGasID,iISO,+1, &
    iNumVibLevels,raPressVT,raKineticVT,raQtipsVT,raNLTETempVT,dVibCenter)

!      do iI = 1,kProfLayer
!        print *,iI,raPressVT(iI),raKineticVT(iI),raNLTETempVT(iI),'+',pProf(iI),raTTemp(iI)
!      end do
!      call dostop

    iLogOrLinear = -1        !! linear interp in log(P)
    iLogOrLinear = +1        !! spline interp in P
    iLogOrLinear = 0         !! linear interp in P, Dave Edwards way

    iDefault = 0
    IF (iDefault /= iLogOrLinear) THEN
        write(kStdErr,*) 'In read_nonlte_temperature iDefault,iLogOrLinear'
        write(kStdErr,*) 'are different : ',iDefault,iLogOrLinear
    END IF

    CALL MapVTprofile2klayers( &
    iLogOrLinear,iGASID,iISO,daJU(1),dVibCenter,pProf, &
    raPressVT,raKineticVT,raNLTETempVT,raQTipsVT,iNumVibLevels, &
    iProfileLayers,raLTETemp,raNLTETemp,raVibQFT)
          
    iStart = -1
    DO iI = 1,kProfLayer
    !! total molecules in layer
        rFac           = pProf(iI) * 100.0 * raThickness(iI) / kMGC / raTTemp(iI) * 1.0e-4
        raMixRatio(iI) = raTAmt(iI)*1000 / rFac * 1.0e6
    END DO

    IF (iAllLayersLTE == +1) THEN
    !!!all layers in LTE; put all temperatures same as KLAYERS temperatures
    !!!before March 3, 2005 had raLTETemp = raNLTETemp == raTTemp, raVibQFT=1

    !!! after March 3, 2005, just make sure NLTETemp = LTETemp, keep
    !!! the VibQFT as was done by computeQVlte.f
        DO iI = 1,kProfLayer
            raNLTETemp(iI) = raLTETemp(iI)
        END DO
    END IF

!! when printing out this stuff with kLongOrShort == -2,+2, assume you are only
!! doing NLTE ie running 2205-2405 cm-1
    IF ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        write(kStdWarn,120) '        klayers                            ||        vibtemp file               |  '
        write(kStdWarn,120) 'iI     pavg      dz        q         T(iI) || Tlte      dT   |     Tnlte     dT |  MixRatio'
        write(kStdWarn,120) '       mb        m       kmol/cm2          ||       K        |         K        |    ppmv'
        write(kStdWarn,120) '-------------------------------------------------------------------------------------------'
        iStart = (kProfLayer-iProfileLayers+1)	
        DO iI = iStart,kProfLayer
            IF (raLTETemp(iI) < 100.0) THEN
                iStart = iI+1
            END IF
            write(kStdWarn,987) iI, &
            pProf(iI),raThickness(iI),raTAmt(iI),raTTemp(iI), &
            raLTETemp(iI),  raLTETemp(iI)-raTTemp(iI), &
            raNLTETemp(iI), raNLTETemp(iI)-raTTemp(iI), &
            raMixRatio(iI)
        END DO

        write(kStdWarn,*) '  '
        write(kStdWarn,*) '  Comparing NLTE profile to US STD Profile, 2350 Band  '
        write(kStdWarn,*) '  '

        ca1 = 'iI   Pavg    Tk(klayers) | Tk(VibT)          dT  |     Tv           dT '
        ca2 = '-------------------------|-----------------------|----------------------'

        IF (iBand == 1) THEN
            CALL GetUSSTD_2350(pProf,iStart,daJL,daJU,iaJ_UorL,raLTE_STD,raNLTE_STD)
            write(kStdWarn,*) ca1
            write(kStdWarn,*) ca2
            DO iI = iStart, kProfLayer
                write(kStdWarn,1234) iI,pProf(iI),raLTETemp(iI),'|', &
                raTTemp(iI), raTTemp(iI)-raLTETemp(iI),'|', &
                raNLTETemp(iI),raNLTETemp(iI)-raLTETemp(iI)
            END DO
        END iF

        IF (iBand == 1) THEN
            write(kStdWarn,*) ' '
            write(kStdWarn,*) ' Comparisons of Tk,tNLTE vs USSTD : '
            ca1 = 'iI   Pavg    Tk(klayers)      TStd         dT  |     Tv      TvSTD      dTv'
            ca2 = '-----------------------------------------------|------------------------------'

            write(kStdWarn,*) ca1
            write(kStdWarn,*) ca2
            DO iI = iStart, kProfLayer
                write(kStdWarn,1250) iI,pProf(iI),raLTETemp(iI),raLTE_STD(iI),raLTETemp(iI)-raLTE_STD(iI),'|', &
                raNLTETemp(iI), raNLTE_STD(iI),raNLTETemp(iI)-raNLTE_STD(iI)
            END DO
        END iF
    END IF

    1234 FORMAT(I3,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1)
    1250 FORMAT(I3,' ',4(F10.5,' '),A1,' ',3(F10.5,' '))
    1070 FORMAT('ERROR! number ',I5,' opening NONLTE profile:',/,A80)
    987  FORMAT(I4,' ',1(E9.4,' '),1(F9.3,' '),1(E9.4,' '),1(F8.3,' '),'||',1(F8.3,' '), &
    &        1(F6.3,' '),'|',2(F8.3,' '),'|',F8.3)
    120  FORMAT(A120)
    
    RETURN
    end SUBROUTINE read_nonlte_temperature

!************************************************************************
! this subroutine maps the VT profile onto the KLAYERS profile
    SUBROUTINE MapVTprofile2klayers( &
    iLogOrLinear,iGASID,iISO,dJU,dVibCenter,pProf, &
    raPressVT,raKineticVT,raNLTETempVT,raQTipsVT,iNumVibLevels, &
    iProfileLayers,raLTETemp,raNLTETemp,raVibQFT)

    IMPLICIT NONE
      
    include '../INCLUDE/kcartaparam.f90'
     
! input vars
    INTEGER :: iLogOrLinear                   !!! what method to use for map
    INTEGER :: iGasID,iISO
    DOUBLE PRECISION :: dVibCenter,dJU
! the VibTemp profiles(p,LTE,NLTE,qvib) info are in these variables
    INTEGER :: iNumVibLevels                  !!!num of vib levels in NLTE prof
    REAL :: raPressVT(kNLTEProfLayer),  raKineticVT(kNLTEProfLayer)
    REAL :: raNLTETempVT(kNLTEProfLayer),raQtipsVT(kNLTEProfLayer)
! the KLAYERS avg pressures are in here
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers
! output vars
    REAL :: raNLTETemp(kProfLayer)   !!!NLTE temperatures
    REAL :: raLTETemp(kProfLayer)    !!!want to send this is from the vib-temp
!!!files as it might be slightly but
!!!importantly different from LAYERS T(z)
    REAL :: raVibQFT(kProfLayer)     !!!modifier to partition fcns
          
! local vars
    INTEGER :: iI,iJ,iLP
! recall things have to be in "ascending" order for the splines
    REAL :: raPress1Swap(kProfLayer),rP,rFac
    REAL :: raLTETemp1Swap(kProfLayer)
    REAL :: raNLTETemp1Swap(kProfLayer)
    REAL :: raQtips1Swap(kProfLayer)

! c      iLogOrLinear = -1        !! linear interp in log(P)
! c      iLogOrLinear = +1        !! spline interp in P
! c      iLogOrLinear = 0         !! linear interp in P, Dave Edwards way

    write(kStdWarn,*) 'map NLTE profile to kLAYERS profile (>= 0.005 mb)'

! -------------------------------------------------->
    IF (iLogOrLinear > 0) THEN
    ! linear interp in log(P)
        DO iI = 1,iNumVibLevels
            iJ = iNumVibLevels - iI + 1
            raPress1Swap(iI)    = log(raPressVT(iJ))
            raLTETemp1Swap(iI)  = raKineticVT(iJ)
            raNLTETemp1Swap(iI) = raNLTETempVT(iJ)
            raQTips1Swap(iI)    = raQTipsVT(iJ)
        END DO

        DO iI = 1,kProflayer
            pProf(iI) = log(pProf(iI))
        END DO

    !!!now interpolate the Dave Edwards NLTE profiles to the KCARTA profile
    ! for all bands specified by user
        CALL linear(raPress1Swap,raLTETemp1Swap,iNumVibLevels,pProf, &
        raLTETemp,kProfLayer)
        CALL linear(raPress1Swap,raNLTETemp1Swap,iNumVibLevels,pProf, &
        raNLTETemp,kProfLayer)
        CALL linear(raPress1Swap,raQTips1Swap,iNumVibLevels,pProf, &
        raVibQFT,kProfLayer)

        DO iI = 1,kProflayer
            pProf(iI) = exp(pProf(iI))
        END DO

    ! -------------------------------------------------->
    ELSEIF (iLogOrLinear < 0) THEN
    ! spline interp in P
        DO iI = 1,iNumVibLevels
            iJ = iNumVibLevels - iI + 1
            raPress1Swap(iI)    = raPressVT(iJ)
            raLTETemp1Swap(iI)  = raKineticVT(iJ)
            raNLTETemp1Swap(iI) = raNLTETempVT(iJ)
            raQTips1Swap(iI)    = raQTipsVT(iJ)
        END DO

    !!!now interpolate the Dave Edwards NLTE profiles to the KCARTA profile
    ! for all bands specified by user
        CALL linear(raPress1Swap,raLTETemp1Swap,iNumVibLevels,pProf,raLTETemp,kProfLayer)
        call spl(raPress1Swap,raNLTETemp1Swap,iNumVibLevels,pProf,raNLTETemp,kProfLayer)
        call spl(raPress1Swap,raQTips1Swap,iNumVibLevels,pProf,raVibQFT,kProfLayer)

    ! -------------------------------------------------->
    ELSEIF (iLogOrLinear == 0) THEN
    !! linear interp in P, Dave Edwards way
        DO iI = kProfLayer-iProfileLayers+1,kProfLayer
            rP = pProf(iI)
            IF (rP > raPressVT(1)) THEN
                iLP = 1
                rFac = &
                LOG(raPressVT(iLP)/rP)/LOG(raPressVT(iLP+1)/raPressVT(iLP))
            ELSEIF (rP <= raPressVT(iNumVibLevels)) THEN
                iLP = iNumVibLevels - 1
                rFac = &
                LOG(raPressVT(iLP)/rP)/LOG(raPressVT(iLP)/raPressVT(iLP+1))
            ELSE
                DO iJ=1,iNumVibLevels-1
                    IF (rP <= raPressVT(iJ) .AND. rP > raPressVT(iJ+1)) THEN
                        iLP = iJ
                        rFac = LOG(raPressVT(iLP)/rP)
                        rFac = rFac/LOG(raPressVT(iLP)/raPressVT(iLP+1))
                    ENDIF
                END DO
            ENDIF

            raLTETemp(iI)  = raKineticVT(iLP) + &
            rFac*(raKineticVT(iLP+1)-raKineticVT(iLP))
            raNLTETemp(iI) = raNLTETempVT(iLP) + &
            rFac*(raNLTETempVT(iLP+1)-raNLTETempVT(iLP))
            raVibQFT(iI)   = raQTipsVT(iLP) + &
            rFac*(raQTipsVT(iLP+1)-raQTipsVT(iLP))
        END DO
    END IF
! -------------------------------------------------->

    RETURN
    end SUBROUTINE MapVTprofile2klayers

!************************************************************************
END MODULE kreadVTprofiles
