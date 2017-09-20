! Copyright 2004
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kreadVTprofiles

IMPLICIT NONE

CONTAINS

!************************************************************************
! this routine reads in Reference Amts for CO2 LA compressed database
    SUBROUTINE LowerAtmNLTERefs(raRPressX,raRPPressX,raRTempx,raRAmtx)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output
    REAL :: raRPressX(kMaxLayer),raRAmtx(kMaxLayer),raRTempx(kMaxLayer)
    REAL :: raRPPressX(kMaxLayer)

! local
    CHARACTER(80) :: caFname,caStr
    INTEGER :: iIOUN,iI,iErr,iX
    REAL :: r1,r2,r3,r4

    caFname = caLA_US_STD_385ppmv

    GOTO 777

    100 write(kStdErr,*) 'Error reading UpperMixRatio info : filename = '
    WRITE(kStdErr,1070) iErr, caFName
    CALL DoStop

    777 CONTINUE
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine LowerAtmNLTERefs, error reading file  ... '
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    END IF
    kTempUnitOpen = +1

    iI = 0

    555 CONTINUE
    read(iIOUN,123) caStr
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSE
        READ (caStr,*) iX,r1,r2,r3,r4
        iI = iI + 1
        raRPressX(iI)  = r1
        raRPPressX(iI) = r2
        raRTempX(iI)   = r3
        raRAmtX(iI)    = r4
    END IF

    20 CONTINUE
    READ(iIOUN,123,END=199) caStr
    READ (caStr,*) iX,r1,r2,r3,r4
    iI = iI + 1
    raRPressX(iI)  = r1
    raRPPressX(iI) = r2
    raRTempX(iI)   = r3
    raRAmtX(iI)    = r4
    GOTO 20

    199 CONTINUE
    close (iIOUN)
    kTempUnitOpen = -1

    123 FORMAT(A80)
    1070 FORMAT('ERROR! number ',I5,' opening loweAtmNLTE profile:',/,A80)

    RETURN
    end SUBROUTINE LowerAtmNLTERefs

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
! this subroutine gets the UA US STD Ref Profile
    SUBROUTINE GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std, &
    raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! this is if we want to see what a std US profile looks like about 0.005 mb
! it assumes the lower atm has CO2 ~ 385 ppmv
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! local vars
    REAL :: r1,r2,r3,r4,r5,r6,r7,r8,r9
    CHARACTER(120) :: caStr
    INTEGER :: iIOUN,iI,iErr,iReason

    GOTO 777

    write(kStdWarn,*) 'SUBR GetUS_Std_UA : caUA_US_STD_385ppmv = ',caUA_US_STD_385ppmv
    100 write(kStdErr,*) 'Error reading GetUS_Std_UA info : filename = '
    WRITE(kStdErr,1070) iErr,caUA_US_STD_385ppmv
    CALL DoStop

    777 CONTINUE
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caUA_US_STD_385ppmv,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine GetUS_Std_UA, error reading file  ... '
        WRITE(kStdErr,1070) iErr, caUA_US_STD_385ppmv
        CALL DoSTOP
    END IF

    kTempUnitOpen = +1

    iI = 0

    555 CONTINUE
!! for IOSTAT look at http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap04/iostatus.html
!! if iOSTAT = 0, everything fine, if -1 this is EOF, if > 0 this is a problem
    read(iIOUN,123,IOSTAT=iReason,ERR=600) caStr
!      print *,caStr,caUA_US_STD_385ppmv
    IF (iReason < 0) GOTO 600       !!!! end of file
    IF (iReason > 0) THEN
        write(kStdErr,*) 'at iT = ',iI,' IOSTAT=iReason for file ',iReason,caUA_US_STD_385ppmv
        CALL DoStop
    END IF
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSE
        READ (caStr,*) iUpperStd_Num,r1,r2,r3,r4,r5,r6,r7,r8,r9
    !        print *,iI+1,iUpperStd_Num,r1,r2,r3,r4,r5
        iI = iI + 1
        raUpperPress_Std(iI)    = r1
        raUpperDZ_Std(iI)       = r2
        raUpperCO2Amt_Std(iI)   = r3
        raUpperTemp_Std(iI)     = r4
        raUpperMixRatio_Std(iI) = r9
        GOTO 555
    END IF

    600 CONTINUE
    close (iIOUN)
    kTempUnitOpen = -1

    iUpperStd_Num = iI

    123 FORMAT(A120)
    1070 FORMAT('ERROR! number ',I5,' opening US Std UA profile:',/,A80)

    RETURN
    end SUBROUTINE GetUS_Std_UA

!************************************************************************

! this routine sees if we need to add in the 0.005 mb point to the NLTE
! profile that was just read in (else we will either double count the
! CO2 optical depth or mistakenly neglect it; not good in either case!)
    SUBROUTINE Add005mblevel( &
    iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTETemp1)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airslevelsparam.f90'

! input/output
    INTEGER :: iNumVibLevels
    REAL :: raTTemp1(kNLTEProfLayer),raTPress1(kNLTEProfLayer)
    REAL :: raQtips1(kNLTEProfLayer),raNLTETemp1(kNLTEProfLayer)

! local vars
    INTEGER :: iAIRS_TOA,iAddNewPt,iI,iSwap
    REAL :: rAIRS_TOA,rY1,rY2,rY3,raTemp(kNLTEProfLayer)
    REAL :: raTTempX(kNLTEProfLayer),raTPressX(kNLTEProfLayer)
    REAL :: raQtipsX(kNLTEProfLayer),raNLTETempX(kNLTEProfLayer)
    REAL :: raLog(kNLTEProfLayer)

! see if we need to swap values (from small to large)
    iSwap = -1        !!! assume no swap
    IF (raTPress1(1) > raTPress1(iNumVibLevels)) THEN
        iSwap = +1
        DO iI = 1,iNumVibLevels
            raTemp(iI) = raTPress1(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raTPressX(iI) = raTemp(iI)
        END DO

        DO iI = 1,iNumVibLevels
            raTemp(iI) = raTTemp1(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raTTempX(iI) = raTemp(iI)
        END DO

        DO iI = 1,iNumVibLevels
            raTemp(iI) = raQtips1(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raQtipsX(iI) = raTemp(iI)
        END DO

        DO iI = 1,iNumVibLevels
            raTemp(iI) = raNLTETemp1(iNumVibLevels-iI+1)
        END DO
        DO iI = 1,iNumVibLevels
            raNLTETempX(iI) = raTemp(iI)
        END DO
    END IF

    DO iI = 1,iNumVibLevels
        raLog(iI) = log(raTPressX(iI))
    END DO

    IF (raTPressX(1) > DATABASELEV(kMaxLayer+1)) THEN
        write(kStdErr,*) 'usual AIRS layers lowest,highest press level = ', &
        DATABASELEV(1),DATABASELEV(kMaxLayer+1)
        write(kStdErr,*) 'from VT files, lowest,highest press level = ', &
        raTPressX(1),raTPressX(iNumVibLevels)
        write(kStdErr,*) 'hmm rather silly, is it not????'
        write(kStdErr,*) 'expect (raTPressX(iNumVibLevels) < rAIRS_TOA)'
        CALL DOStop
    END IF

! now check to see if explicitly have to add in the levelpoint = 0.005 mb into
!  all this, as that is where the usual AIRS levels end
    rAIRS_TOA = DATABASELEV(kMaxLayer+1)
    iAIRS_TOA = 1
    80 CONTINUE
    IF (raTPressX(iAIRS_TOA) <= rAIRS_TOA) THEN
        iAIRS_TOA = iAIRS_TOA + 1
        GOTO 80
    END IF

! check to see if this boundary [iAIRS_TOA - 1, iAIRS_TOA] is "far" from
! 0.005 mb as we want layers of about 0.5 km thickness when press ~ 0.005 mb
! (see ~/KCARTA/INCLUDE/airslevels_upperparam.f90)
    iAddNewPt = -1      !!! assume it is close
    IF ((abs(raTPressX(iAIRS_TOA)-rAIRS_TOA)/rAIRS_TOA    > 1e-3) .AND. &
    (abs(raTPressX(iAIRS_TOA-1)-rAIRS_TOA)/rAIRS_TOA  > 1e-3)) THEN
        iAddNewPt = +1      !!! AIRS TOA is kinda far from both bracket pts
        CALL rspl1(raLog,raTTempX,   iNumVibLevels,log(rAIRS_TOA),rY1,1)
        CALL rspl1(raLog,raQtipsX,   iNumVibLevels,log(rAIRS_TOA),rY2,1)
        CALL rspl1(raLog,raNLTETempX,iNumVibLevels,log(rAIRS_TOA),rY3,1)

    !! move things over by 1
        DO iI = iNumVibLevels,iAIRS_TOA,-1
            raTPressX(iI+1)   = raTPressX(iI)
            raTTempX(iI+1)    = raTTempX(iI)
            raQTipsX(iI+1)    = raQTipsX(iI)
            raNLTETempX(iI+1) = raNLTETempX(iI)
        END DO

    !! replace elements at this index with the new values
        iI = iAIRS_TOA
        raTPressX(iI)   = rAIRS_TOA
        raTTempX(iI)    = rY1
        raQTipsX(iI)    = rY2
        raNLTETempX(iI) = rY3
         
        iNumVibLevels = iNumVibLevels + 1
    END IF

    IF (iAddNewPt > 0) THEN
    ! we added a point and need to replace raYYY1 with raYYYX
        IF (iSwap < 0) THEN
        ! no swaps
            DO iI = 1,iNumVibLevels
                raTPress1(iI)   = raTPressX(iI)
                raTTemp1(iI)    = raTTempX(iI)
                raQTips1(iI)    = raQTipsX(iI)
                raNLTETemp1(iI) = raNLTETempX(iI)
            END DO
        ELSEIF (iSwap > 0) THEN
            DO iI = 1,iNumVibLevels
                raTemp(iI) = raTPressX(iNumVibLevels-iI+1)
            END DO
            DO iI = 1,iNumVibLevels
                raTPress1(iI) = raTemp(iI)
            END DO

            DO iI = 1,iNumVibLevels
                raTemp(iI) = raTTempX(iNumVibLevels-iI+1)
            END DO
            DO iI = 1,iNumVibLevels
                raTTemp1(iI) = raTemp(iI)
            END DO

            DO iI = 1,iNumVibLevels
                raTemp(iI) = raQtipsX(iNumVibLevels-iI+1)
            END DO
            DO iI = 1,iNumVibLevels
                raQtips1(iI) = raTemp(iI)
            END DO

            DO iI = 1,iNumVibLevels
                raTemp(iI) = raNLTETempX(iNumVibLevels-iI+1)
            END DO
            DO iI = 1,iNumVibLevels
                raNLTETemp1(iI) = raTemp(iI)
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE Add005mblevel

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
    iUpper,raUpperTemp,raUpperGasAmt,raUpperPress,raUpperPartPress, &
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
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
! this tells the pressure levels and layer thicknesses for upper atm
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)

! local variables
    INTEGER :: iISO,iAllOrSome,iNumVibLevels
    REAL :: rMixRatio,MGC,rAvgHgt
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

    IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        write(kStdWarn,*) 'LOWER atm LTE temperatures, gas amounts ....'
        write(kStdWarn,*) 'read in direct from GENLN2 NLTE profile file ',caFName
        write(kStdWarn,*) ' this part is for jollies, not really used ....'
        write(kStdWarn,*) 'iI     pavg   dz   q        T(iI)  '
        write(kStdWarn,*) '       (mb)   (m) (mol/cm2)  (K)'
        write(kStdWarn,*) '---------------------------------------------------'
    END IF

    DO iI = 1,kProfLayer
        raUpperPress(iI)     = 0.0
        raUpperPartPress(iI) = 0.0
        raUpperTemp(iI)      = 0.0
        raUpperGasAmt(iI)    = 0.0
    END DO

!!!dope, for the regular kCARTA atmosphere, just a check
    DO iI = 1,iLay-1
        iJ = iI
        raUpperThickness(iI) = raThick1(iJ)  !!! in m
        raUpperPress(iI)     = raPavg1(iJ)   !!! in mb
        raUpperTemp(iI)      = raTavg1(iJ)
    ! ppmv = (number of gas molecules)/(number of air molecules) * 1e6
    ! pV = NRT ==> N(g) = p(g) V(a)/RT(g),  N(a) = p(a) V(a)/ RT(a)
    ! now V(g) == V(a), T(g) == T(a) and p(total) = p(a) + p(g)
    ! thus ppmv = N(g)/N(a) = p(g)/p(a) = p(g)/(p(total) - p(g)) * 1e6
    ! thus p(g) = (p(total)/(1e6+ppmv)) ppmv
        rAvgHgt = (raLayTop1(iJ)+raLayBot1(iJ))/2.0
        CALL rspl1(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        raUpperPress(iI),rMixRatio,1)
        raUpperPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raUpperPartPress(iI)=raUpperPress(iI)*rMixRatio/(1.0e6+rMixRatio) !atm
        raUpperGasAmt(iI) = raThick1(iJ) * 100.0           !!! in cm
        raUpperGasAmt(iI)=raUpperGasAmt(iI)*kAtm2mb*100.0*raUpperPartPress(iI)
        raUpperGasAmt(iI) = raUpperGasAmt(iI)/1e9/MGC/raUpperTemp(iI)

        IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1))) THEN
            write(kStdWarn,6543) iI,raUpperPress(iI)*kAtm2mb,raUpperThickness(iI), &
            raUpperGasAmt(iI),raUpperTemp(iI)
        END IF
    END DO

!! now concentrate on upper atm
    DO iI = 1,kProfLayer
        raUpperPress(iI)     = 0.0
        raUpperPartPress(iI) = 0.0
        raUpperTemp(iI)      = 0.0
        raUpperGasAmt(iI)    = 0.0
    END DO

    IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'UPPER atm LTE temperatures, gas amounts ....'
        write(kStdWarn,*) 'read in direct from GENLN2 NLTE profile file ',caFName
        write(kStdWarn,*) 'iI     pavg   dz   q        T(iI)  '
        write(kStdWarn,*) '       (mb)   (m) (mol/cm2)  (K)'
        write(kStdWarn,*) '---------------------------------------------------'
    END IF
     
    DO iI = 1,iUpper
        iJ = iLay+(iI-1)
        raUpperThickness(iI) = raThick1(iJ)  !!! in m
        raUpperPress(iI)     = raPavg1(iJ)   !!! in mb
        raUpperTemp(iI)      = raTavg1(iJ)
    ! ppmv = (number of gas molecules)/(number of air molecules) * 1e6
    ! pV = NRT ==> N(g) = p(g) V(a)/RT(g),  N(a) = p(a) V(a)/ RT(a)
    ! now V(g) == V(a), T(g) == T(a) and p(total) = p(a) + p(g)
    ! thus ppmv = N(g)/N(a) = p(g)/p(a) = p(g)/(p(total) - p(g)) * 1e6
    ! thus p(g) = (p(total)/(1e6+ppmv)) ppmv
        rAvgHgt = (raLayTop1(iJ)+raLayBot1(iJ))/2.0
        CALL rspl1(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        raUpperPress(iI),rMixRatio,1)
        raUpperPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raUpperPartPress(iI)=raUpperPress(iI)*rMixRatio/(1.0e6+rMixRatio) !!atm
        raUpperGasAmt(iI) = raThick1(iJ) * 100.0           !!! in cm
        raUpperGasAmt(iI) =raUpperGasAmt(iI)*kAtm2mb*100.0*raUpperPartPress(iI)
        raUpperGasAmt(iI) = raUpperGasAmt(iI)/1e9/MGC/raUpperTemp(iI)

        IF  ((iBand == 1) .AND. (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1))) THEN
            write(kStdWarn,6543) iI,raUpperPress(iI)*kAtm2mb,raUpperThickness(iI), &
            raUpperGasAmt(iI),raUpperTemp(iI)
        END IF
    END DO

    1070 FORMAT('ERROR! number ',I5,' opening upper atm NONLTE profile:',/,A80)
    6543 FORMAT(I4,' ',1(E9.4,' '),1(F9.3,' '),1(E9.4,' '),1(F8.3,' '))

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
    CHARACTER(105) :: ca105A,ca105B,ca105C,ca105D

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

    DO iI = 1,kProfLayer
        raTPress(iI)     = 0.0
        raTPartPress(iI) = 0.0
        raTTemp(iI)      = 0.0
        raTamt(iI)       = 0.0
    END DO

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
        CALL rspl1(raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        raTPress(iI),rMixRatio,1)
        raUAMixRatio(iI) = rMixRatio
        raTPress(iI)     = raPavg1(iJ)/kAtm2mb      !!! in atm
        raTPartPress(iI) = raTPress(iI) * rMixRatio/(1.0e6 + rMixRatio) !!! atm
        raTAmt(iI)       = raThick1(iJ) * 100.0           !!! in cm
        raTAmt(iI)       = raTAmt(iI)*kAtm2mb*100.0*raTPartPress(iI)
        raTAmt(iI)       = raTAmt(iI)/1e9/MGC/raTTemp(iI)
    !        write(kStdWarn,1080) kProfLayer+iI,raTPress(iI)*kAtm2mb,
    !     $            raTPartPress(iI)*kAtm2mb,raTTemp(iI),raTamt(iI),rMixRatio
    END DO

! ccc this is typical kProfLayer (~ 100) CO2 amounts
! ccc 200 2 0.93681E-05 0.34536E-08   198.16440  0.14706E-09

    IF ((iBand == 1) .AND. &
    (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        ca105A = &
        '       klayers                              |            vibtemp file'
        ca105B = &
        'iI        pavg    dz      q         T(iI)   |  Tlte       dT  |     Tnlte        dT | MixRatio '
        ca105C = &
        '          mb       m  kmol/cm2  K           |                 |                     |   ppmv'
        ca105D = &
        '----------------------------------------------------------------------------------------------'
        write(kStdWarn,*) ca105A
        write(kStdWarn,*) ca105B
        write(kStdWarn,*) ca105C
        write(kStdWarn,*) ca105D
        DO iI = 1,iUpper
            write(kStdWarn,9876) iI+kProfLayer,raTPress(iI)*kAtm2mb,raUpperThickness(iI),raTamt(iI),raTTemp(iI), &
            raTTemp(iI),00.0, &
            raNLTETemp(iI), raNLTETemp(iI)-raTTemp(iI),raUAMixRatio(iI)
        END DO

    !! compare what we have to US Std UA
        write(kStdWarn,*) 'Comparison between user supplied and US STD kCompressed 35 UA layers'
        ca105A = &
        '             CURRENT PROF (N)  LAYERS                 |            US STD UA (35) LAYERS'
        ca105B = &
        'iI     pavg       dz      q         T(iI)   MixRatio  |  pavg       dz      q         T(iI)   MixRatio'
        ca105C = &
        '       mb         m  kmol/cm2        K        ppmv    |  mb         m  kmol/cm2        K        ppmv  '
        ca105D = &
        '------------------------------------------------------------------------------------------------------'
        write(kStdWarn,*) ca105A
        write(kStdWarn,*) ca105B
        write(kStdWarn,*) ca105C
        write(kStdWarn,*) ca105D
        DO iI = 1,max(iUpper,iUpperStd_Num)
            write(kStdWarn,4321) iI+kProfLayer,raTPress(iI)*kAtm2mb,raUpperThickness(iI),raTamt(iI), &
            raTTemp(iI),raUAMixRatio(iI), &
            raUpperPress_Std(iI),raUpperDZ_Std(iI),raUpperCO2Amt_Std(iI),raUpperTemp_Std(iI), &
            raUpperMixRatio_Std(iI)
        END DO
    END IF

    1070 FORMAT('ERROR! number ',I5,' opening upper atm NONLTE profile:',/,A80)
    1080 FORMAT(I4,' ',2(E9.4,'  '),1(F10.4,'  '),E12.5,F10.4)
    9876 FORMAT(I4,' ',E9.4,' ',F9.4,' ',E9.4,' ',F9.4,' | ',4(F9.4,' '),' | ',F9.4)
    4321 FORMAT(I4,' ',1(E9.4,' ',F9.4,' ',E9.4,' ',2(F9.4,' ')),'|',1(E9.4,' ',F9.4,' ',E9.4,' ',2(F9.4,' ')))

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
        write(kStdWarn,*) '        klayers                            ||        vibtemp file               |  '
        write(kStdWarn,*) 'iI     pavg      dz        q         T(iI) || Tlte      dT   |     Tnlte     dT |  MixRatio'
        write(kStdWarn,*) '       mb        m       kmol/cm2          ||       K        |         K        |    ppmv'
        write(kStdWarn,*) '-------------------------------------------------------------------------------------------'
        DO iI = 1,kProfLayer
            IF (raLTETemp(iI) < 100.0) THEN
                iStart = iI+1
            END IF
            iStart = (kProfLayer-iProfileLayers+1)
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
            ca2 = '-----------------------------------------------|----------------------------'

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
    987 FORMAT(I4,' ',1(E9.4,' '),1(F9.3,' '),1(E9.4,' '),1(F8.3,' '),'||',1(F8.3,' '), &
    &        1(F6.3,' '),'|',2(F8.3,' '),'|',E8.3)

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
        CALL rlinear(raPress1Swap,raLTETemp1Swap,iNumVibLevels,pProf, &
        raLTETemp,kProfLayer)
        CALL rlinear(raPress1Swap,raNLTETemp1Swap,iNumVibLevels,pProf, &
        raNLTETemp,kProfLayer)
        CALL rlinear(raPress1Swap,raQTips1Swap,iNumVibLevels,pProf, &
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
        CALL rlinear(raPress1Swap,raLTETemp1Swap,iNumVibLevels,pProf, &
        raLTETemp,kProfLayer)
        CALL rspl(raPress1Swap,raNLTETemp1Swap,iNumVibLevels,pProf, &
        raNLTETemp,kProfLayer)
        CALL rspl(raPress1Swap,raQTips1Swap,iNumVibLevels,pProf, &
        raVibQFT,kProfLayer)

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
! this generic routine reads in Dave Edwards NLTE profiles
! if necessary, it adds on into at 0.005 mb by calling Add005mblevel

! as of March 2015, if the code CANNOT find the band in question, it gives
! a warning because it replaces the not-found data with raKineticVT

    SUBROUTINE ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL, &
    iGasID,iISO,iAllORSome, &
    iNumVibLevels,raPressVT,raKineticVT,raQtipsVT,raNLTETempVT, &
    dVibCenter)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    CHARACTER(80) :: caFname                !!! file to read
    DOUBLE PRECISION :: daJU(kHITRAN),daJL(kHITRAN) !!! quantum numbers
    INTEGER :: iGasID, iISO                !!! what to read
    INTEGER :: iAllOrSome                  !!! read NLTE (+1) or LTE (-1) stuff?
! output params
    REAL :: raKineticVT(kNLTEProfLayer),raPressVT(kNLTEProfLayer)
    REAL :: raNLTETempVT(kNLTEProfLayer),raQtipsVT(kNLTEProfLayer)
    INTEGER :: iNumVibLevels,iaJ_UorL(kHITRAN)
    DOUBLE PRECISION :: dVibCenter

! local variables
    REAL :: p,pp,q,t,match_band_energy2
    CHARACTER(80) :: caStr
    CHARACTER(3) ::  ca3
    REAL :: raQtipsXVT(kNLTEProfLayer),rDummyVibEnergy
    INTEGER :: iIOUN,iErr,iInputJU,iInputJL,iSetQVT
    INTEGER :: i1,i2,iGasesInFile,iVibs0,iVibs,iI,iJ,iK
    INTEGER :: iaGasIDs(kGasComp),iaVibTemps(kGasComp)
    INTEGER :: iDummy1,iDummyGasID,iDummyISO,iDummyNum,iDummyQuantum
    INTEGER :: iaDummyGasID(kNLTEProfLayer)
    INTEGER :: iMatchTry

    iMatchTry = 0   !!!nothing found

!!!first try to match the UUPER quantum numbers; if no match found, try
!!!to match the LOWER quantum numbers; if nothing is found, give up!

    DO i1 = 1,kHITRAN
        iaJ_UorL(i1) = 0
    END DO

    111 CONTINUE
    iMatchTry = iMatchTry + 1

    dVibCenter = -1.0d0

    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
    ! his does not get executed as error messgs get handled by line 100
        write (kStdErr,*) 'in subroutine ReadGENLN2_NLTE_Profile, '
        write (kStdErr,*) 'error opening file  ... '
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    ELSE
        GOTO 777
    END IF

    100 write(kStdErr,*) 'In ReadGENLN2_NLTE_Profile, Error reading NLTE data'
    WRITE(kStdErr,1070) iErr, caFName
    CALL DoStop

    777 CONTINUE
    kTempUnitOpen=1

    iInputJU = nint(daJU(1))
    iInputJL = nint(daJL(1))
    write(kStdWarn,*) ' try to get NLTE profile for gid,iso,JL,JU = ',iGASID,iISO,iInPutJL,iInputJU
    write(kStdWarn,*) ' from user supplied profile(s) in ',caFName
         
!! Dave Edwards info does start GND = 1, TOA = 120 so things are fine
    123 FORMAT(A80)

    555 CONTINUE
    read(iIOUN,123) caStr
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSEIF (caStr(1:3) == '***') THEN
    !!!aha .. now we are cooking; see ntepro.f in Genln2
    !!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV
    !!!where FLAG  = ***
    !!!      MODEL = atmosphere model number for data file
    !!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL
    !!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB
    !!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL
    !!! so read in flag, model,nset,ngas,nlev where flag = '***'
        read(caStr,*) ca3,i1,i2,iGasesInFile,iNumVibLevels
    END IF

!!!! now read number of (gas,isotope) pairs in the file
!!!! Genln2 directly sez
!!!!   "READ GAS, ISOTOPE PAIRS FOR WHICH THERE IS VIB. TEM. DATA"
    iVibs0 = 0
    iK = 0
    DO iI = 1,kNLTEProfLayer
        iaDummyGasID(iI) = -1
    END DO
    DO iI = 1,iGasesInFile
        read(iIOUN,*) iaGasIDs(iI),iaVibTemps(iI)
        iVibs0 = iVibs0 + iaVibTemps(iI)
        DO iJ = 1,iaVibTemps(iI)
            iK = iK + 1
            iaDummyGasID(iK) = iaGasIDs(iI)
        END DO
    END DO

!!!read the presssure levels
    read(iIOUN,*) (raPressVT(iI),iI=1,iNumVibLevels)   !!! in mb

!!! read the kinetic temperatures
!!! this might be SLIGHTLY different from KLAYERS LTE temps, but it is
!!! important to use THESE temperatures when computing the population
!!! ratios r1,r2, and associated parameters
    read(iIOUN,*) iDummy1
    IF (iDummy1 /= 0) THEN
        write(kStdErr,*) 'Expected a "0" between pressure levels and LTE temps'
        write(kStdErr,*) 'during reading in the NLTE profile : '
        write(kStdErr,*) caFName
        CALL DOStop
    END IF
!!! read the kinetic temp (LTE)
    read(iIOUN,*) (raKineticVT(iI),iI=1,iNumVibLevels)
          
    IF (iAllOrSome == -1) THEN
        GOTO 999    !!!do not need the NLTE stuff
    END IF

!!!read in the QTIPS_VIB vib partition function profiles
!!!dummy string that says "!QV   Vibrational Partition functions"
    read(iIOUN,123) caStr
    iVibs = 0
    iK = 0
    iSetQVT = -1
    800 CONTINUE
    read(iIOUN,*) iDummyGasID,iDummyISO
    read(iIOUN,*) (raQtipsXVT(iI),iI=1,iNumVibLevels)
    iK = iK + 1
    IF (iDummyGasID /= iaDummyGasID(iK)) THEN
        write(kStdErr,*) 'GasID mismatch when reading in vib part fcn profs'
        write(kStdErr,*) 'Expected ',iaDummyGasID(iK),' but got ',iDummyGasID
        CALL DOStop
    END IF
    IF ((iGASID == iDummyGasID) .AND. (iISO == iDummyISO)) THEN
        iSetQVT = +1
        DO iI = 1,iNumVibLevels
            raQtipsVT(iI) = raQtipsXVT(iI)
        END DO
    END IF
    iVibs = iVibs + 1
    IF (iVibs < iVibs0) THEN
        GOTO 800
    ELSE
        GOTO 820
    END IF

    820 CONTINUE

!!!dummy string that says "!TV   Vibrational Temperatures"
    read(iIOUN,123) caStr

!!! read the NLTE temperatures
    888 CONTINUE
    read(iIOUN,*,END=999) iDummyNum,iDummyGASID,iDummyISO,iDummyQuantum, &
    rDummyVibEnergy
    dVibCenter = rDummyVibEnergy * 1.0d0
    read(iIOUN,*) (raNLTETempVT(iI),iI=1,iNumVibLevels)     !!! in K

    IF (iMatchTry == 1) THEN      !!!trying to match upper quantum numbers
        IF ((iGASID == iDummyGasID) .AND. (iDummyQuantum == iInputJU) &
         .AND. (iISO == iDummyISO)) THEN
            write(kStdWarn,*) 'matched upper quantum number'
            DO iI = 1,kHitran
                iaJ_UorL(iI) = +1         !!!matched the upper quantum number
            END DO
            GOTO 999
        ELSE
            GOTO 888
        END IF
    END IF

    IF (iMatchTry == 2) THEN      !!!trying to match lower quantum numbers
        IF ((iGASID == iDummyGasID) .AND. (iDummyQuantum == iInputJL) &
         .AND. (iISO == iDummyISO)) THEN
            write(kStdWarn,*) 'matched lower quantum number'
            DO iI = 1,kHitran
                iaJ_UorL(iI) = -1         !!!matched the lower quantum number
            END DO
            GOTO 999
        ELSE
            GOTO 888
        END IF
    END IF

    999 CONTINUE

    close (iIOUN)
    kTempUnitOpen=-1

    IF ((iMatchTry == 1) .AND. (iaJ_UorL(1) == 0) .AND. &
    (iAllOrSOme > 0)) THEN
        write(kStdWarn,*) 'Did not find a match for IUSGQ; try ILSGQ ...'
        GOTO 111
    END IF

    IF ((iMatchTry == 2) .AND. (iaJ_UorL(1) == 0) .AND. &
    (iAllOrSome > 0)) THEN
        write(kStdErr,*) 'Did not find a match for IUSGQ or ILSGQ ...'
        write(kStdErr,*) '  iGASID,iISO,iInputJL,iInputJU = ',iGASID,iISO,iInputJL,iInputJU
        write(kStdErr,*) '  replacing with Kintetic (LTE) temp'

        write(kStdWarn,*) 'Did not find a match for IUSGQ or ILSGQ ...'
        write(kStdWarn,*) '  iGASID,iISO,iInputJL,iInputJU = ',iGASID,iISO,iInputJL,iInputJU
        write(kStdWarn,*) '  replacing with Kintetic (LTE) temp'

        dVibCenter = match_band_energy2(iGASID,iISO,iInputJU)*1.0d0
        DO iI = 1,kHitran
            iaJ_UorL(iI) = +1         !!!matched the upper quantum number
        END DO

    !!! did not find what we are looking for, set NLTE output temp to this
        DO iI = 1,iNumVibLevels
            raNLTETempVT(iI) = raKineticVT(iI)
        END DO
        IF (iSetQVT < 0) THEN
            write(kStdErr,*) 'Did not find partition fcns for (GasID ISO), setting to 1.0 ',iGASID,iISO
            write(kStdWarn,*) 'Did not find partition fcns for (GasID ISO), setting to 1.0 ',iGASID,iISO
            DO iI = 1,iNumVibLevels
                raQtipsVT(iI) = 1.0
            END DO
        END IF
    !      CALL DoStop
    END IF

    IF ((iAllORSome > 0) .AND. (dVibCenter < 1.0d-2)) THEN
        write (kStdErr,*) 'need to set dVibCenter, but seem not to have'
        write (kStdErr,*) 'been able to do so!'
        Call DoStop
    END IF

    CALL Add005mblevel(iNumVibLevels, &
    raPressVT,raKineticVT,raQtipsVT,raNLTETempVT)

    1070 FORMAT('ERROR! number ',I5,' opening Read GENLN2 nlte profile:',/,A80)

    RETURN
    end SUBROUTINE ReadGENLN2_NLTE_Profile

!************************************************************************
! this function matches up CO2 band center energies with ISO and HITRAIN id
    REAL FUNCTION match_band_energy2(iGID,iISO,iHITRAN)

! see  NONLTE/driver_make_nlte_prof.m
!   IL    MOL   IDMOL  ISO   IDISO  LEVEL   IDAFGL  ENERGY(cm-1)
!    1 CO2       2    626      1     1101      2     667.38000
!    2 CO2       2    626      1    10002      3    1285.40900
!    3 CO2       2    626      1     2201      4    1335.13200
!    4 CO2       2    626      1    10001      5    1388.18500
!    5 CO2       2    626      1    11102      6    1932.47000
!    6 CO2       2    626      1     3301      7    2003.24600
!    7 CO2       2    626      1    11101      8    2076.85600
!    8 CO2       2    626      1       11      9    2349.14300
!    9 CO2       2    626      1    20003     10    2548.36700
!   10 CO2       2    626      1    12202     11    2585.02200
!   11 CO2       2    626      1    20002     12    2671.14300
!   12 CO2       2    626      1     4401     13    2671.71700
!   13 CO2       2    626      1    12201     14    2760.72500
!   14 CO2       2    626      1    20001     15    2797.13600
!   15 CO2       2    626      1     1111     16    3004.01200
!   16 CO2       2    626      1    10012     23    3612.84200
!   17 CO2       2    626      1     2211     24    3659.27300
!   18 CO2       2    626      1    10011     25    3714.78300
!   19 CO2       2    636      2     1101      2     648.47802
!   20 CO2       2    636      2    10002      3    1265.82800
!   21 CO2       2    636      2     2201      4    1297.26400
!   22 CO2       2    636      2    10001      5    1370.06300
!   23 CO2       2    636      2       11      9    2283.48800
!   24 CO2       2    636      2     1111     16    2920.23900
!   25 CO2       2    636      2    10012     23    3527.73800
!   26 CO2       2    636      2     2211     24    3580.75000
!   27 CO2       2    636      2    10011     25    3632.91000
!   28 CO2       2    628      3     1101      2     662.37300
!   29 CO2       2    628      3    10002      3    1259.42600
!   30 CO2       2    628      3     2201      4    1325.14100
!   31 CO2       2    628      3    10001      5    1365.84400
!   32 CO2       2    628      3    11102      6    1901.73700
!   33 CO2       2    628      3     3301      7    1988.32800
!   34 CO2       2    628      3    11101      8    2049.33910
!   35 CO2       2    628      3       11      9    2332.11300
!   36 CO2       2    627      4     1101      2     664.72900
!   37 CO2       2    627      4    10002      3    1272.28700
!   38 CO2       2    627      4     2201      4    1329.84300
!   39 CO2       2    627      4    10001      5    1376.02700
!   40 CO2       2    627      4    11102      6    1916.69500
!   41 CO2       2    627      4     3301      7    1995.35200
!   42 CO2       2    627      4    11101      8    2062.09910
!   43 CO2       2    627      4       11      9    2340.01400
!   44 CO2       2    626      1    11112     36    4247.70500
!   45 CO2       2    626      1     3311     37    4314.91300
!   46 CO2       2    626      1    11111     38    4390.62900
!   47 CO2       2    628      3     1111     16    2982.11200

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iGID,iISO,iHITRAN
! output vars
    REAL :: rJunk
      
! local vars
    INTEGER :: iaISO_table(47)
    INTEGER :: iaID_hitranID(47)
    REAL ::    raEnergy(47)
    INTEGER :: iFound, iI,iMax

    DATA iaISO_table/1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1,     1, &
    &                  1,     1,     1,     1,     1,     1,     2,     2,     2,     2,     2,     2, &
    &                  2,     2,     2,     3,     3,     3,     3,     3,     3,     3,     3,     4, &
    &                  4,     4,     4,     4,     4,     4,     4,     1,     1,     1,     3/
    DATA iaID_hitranID/2,     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,    13, &
    &                    14,    15,    16,    23,    24,    25,    2,     3,    4,     5,     9,     16, &
    &                    23,    24,    25,     2,     3,     4,    5,     6,    7,     8,     9,     2, &
    &                    3,     4,     5,     6,     7,     8,     9,    36,    37,    38,    16/
    DATA raEnergy/0.6674,    1.2854,    1.3351,    1.3882,    1.9325,    2.0032,    2.0769, &
    &               2.3491,    2.5484,    2.5850,    2.6711,    2.6717,    2.7607,    2.7971, &
    &               3.0040,    3.6128,    3.6593,    3.7148,    0.6485,    1.2658,    1.2973, &
    &               1.3701,    2.2835,    2.9202,    3.5277,    3.5808,    3.6329,    0.6624, &
    &               1.2594,    1.3251,    1.3658,    1.9017,    1.9883,    2.0493,    2.3321, &
    &               0.6647,    1.2723,    1.3298,    1.3760,    1.9167,    1.9954,    2.0621, &
    &               2.3400,    4.2477,    4.3149,    4.3906,    2.9821/
    IF (iGID /= 2) THEN
        write(kStdErr,*) 'need gasID = 2 in REAL FUNCTION match_band_energy2'
        CALL DoStop
    END IF
    IF ((iISO < 1) .OR. (iISO > 4)) THEN
        write(kStdErr,*) 'need iISO = 1,2,3,4 in REAL FUNCTION match_band_energy2'
        CALL DoStop
    END IF
         
    rJunk = -1.0
    iFound = -1

    iMax = 47
    iI = 1
    10 CONTINUE
    IF (iI <= iMax) THEN
        IF ((iaISO_table(iI) == iISO) .AND. (iaID_hitranID(iI) == iHITRAN)) THEN
            iFound = 1
            rJunk = raEnergy(iI)
        ELSE
            iI = iI + 1
            GOTO 10
        END IF
    END IF

    IF (iFound < 0) THEN
        write(kStdErr,*) 'did not find band match in in REAL FUNCTION match_band_energy2'
        CALL DoStop
    END IF

    match_band_energy2 = rJunk * 1000.0

    RETURN
    end FUNCTION match_band_energy2
!************************************************************************
END MODULE kreadVTprofiles
