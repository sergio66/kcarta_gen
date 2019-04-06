! Copyright 2004
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kpredictVT

USE basic_common
USE spline_and_sort_and_common
USE s_writefile
USE n_nonlte_common

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine interpolates the SolarAngle for iRegr, and dumps output
! into a vt file for kCARTA to read
    SUBROUTINE OutputVTFile(caVTFile,iRTP,iRegr,rSolarAngle)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL :: rSolarAngle            !!!solar angle for atm #1
    INTEGER :: iRegr               !!!which regr profile this is closest to
    INTEGER :: iRTP                !!!which of the rtp files processed
! output parameters
    CHARACTER(80) :: caVTFile       !!!temp file name where VT profiles are
!!!created and stored if not given by user

! local vars
! these are used to read in info from the 6 angle profs 0,40,60,80,85,90
! MXXTMP = 6
! kMaxUserSet = 12
! kMaxK = 50
    REAL :: raaPressure(MXXTMP,kNLTEProfLayer)
    REAL :: raaKinetic(MXXTMP,kNLTEProfLayer)
    REAL :: raaaVibPartFcn(MXXTMP,kMaxUserSet,kNLTEProfLayer)
    REAL :: raaaVibTemp(MXXTMP,kMaxK,kNLTEProfLayer)
    CHARACTER(80) :: caFName
    CHARACTER(80) :: caaComments(kMaxLayer)
    INTEGER :: iNumComments
! to dump out the profiles
    REAL :: raaVibPartFcn(kMaxUserSet,kNLTEProfLayer)
    REAL :: raaVibTemp(kMaxK,kNLTEProfLayer)

! other local vars
    REAL :: raSolAngles(MXXTMP)
    INTEGER :: iaSolAngles(MXXTMP)
    INTEGER :: iI,iJ,iK,iAngle,iFound
    INTEGER :: iGasesInFile,iNumVibLevels            !gases,levels in file
    INTEGER :: iGasID,iVibPF                         !gas ID, number of isotopes
    INTEGER :: iNUMVibTempMagic                      !number of vib profiles
    INTEGER :: iGasesInFile1,iNumVibLevels1          !gases,levels in file
    INTEGER :: iGasID1,iVibPF1                       !gas ID, number of isotopes
    INTEGER :: iNUMVibTempMagic1                     !number of vib profiles
    INTEGER :: iNumSolAngles
    INTEGER :: iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
    REAL :: raVibCntr(kMaxK)
    INTEGER :: iaCnt1(kMaxK),iaGasID1(kMaxK),iaISO1(kMaxK),iaAFGL1(kMaxK)
    REAL :: raVibCntr1(kMaxK)
    REAL :: ra6(MXXTMP),rT

    DATA iNumSolAngles/6/
    DATA (raSolAngles(iI),iI=1,6) /0.0,40.0,60.0,80.0,85.0,90.0/
    DATA (iaSolAngles(iI),iI=1,6) /0  ,40  ,60  ,80  ,85  ,90  /

    IF (MXXTMP < 6) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume MXXTMP = 6'
        CALL DoStop
    END IF
    IF (kMaxUserSet < 12) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume kMaxUserSet = 12'
        CALL DoStop
    END IF
    IF (kMaxK < 50) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume kMaxK = 50'
        CALL DoStop
    END IF

    write(kStdWarn,*) 'kCARTA : estimate VibTemps based on following info : '
    write(kStdWarn,*) 'Input RTP profile               = ',iRTP
    write(kStdWarn,*) 'Closest Regression profile      = ',iRegr
    write(kStdWarn,*) 'Solar Angle                     = ',rSolarAngle
    write(kStdWarn,*) 'Output vib temperature filename = ',caVTFile

! start reading in the info!
    DO iAngle = 1,iNumSolAngles
        caFName = kLopezPuertas
        CALL addinfo(caFname,iRegr,iaSolAngles(iAngle))
        CALL VibTempProfRead(iAngle,caFName,caaComments,iNumComments, &
        raaPressure,raaKinetic,raaaVibPartFcn,raaaVibTemp, &
        iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr, &
        iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)
        IF (iAngle == 1) THEN
            iNUMVibTempMagic1  = iNUMVibTempMagic
            iGasesInFile1      = iGasesInFile
            iNumVibLevels1     = iNumVibLevels
            iGasID1            = iGasID
            iVibPF1            = iVibPF
            DO iJ = 1,iNUMVibTempMagic
                iaCnt1(iJ)     = iaCnt(iJ)
                iaGasID1(iJ)   = iaGasID(iJ)
                iaISO1(iJ)     = iaISO(iJ)
                iaAFGL1(iJ)    = iaAFGL(iJ)
                raVibCntr1(iJ) = raVibCntr(iJ)
            END DO
        ELSE
        ! heck to ensure parameters are ok
            IF (iGasesInFile1 /= iGasesInFile) THEN
                write(kStdErr,*) 'Error in iGasesInfile!!!'
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            ELSEIF (iNumVibLevels1 /= iNumVibLevels) THEN
                write(kStdErr,*) 'Error in iNumVibLevels!!!'
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            ELSEIF (iGasID1 /= iGasID) THEN
                write(kStdErr,*) 'Error in iGasID!!!'
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            ELSEIF (iVibPF1 /= iVibPF) THEN
                write(kStdErr,*) 'Error in iVibPF!!!'
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            ELSEIF (iNUMVibTempMagic1 /= iNUMVibTempMagic) THEN
                write(kStdErr,*) 'Error in iNUMVibTempMagic!!!'
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            END IF

        !!!if those tests ok, make sure HITRAN IDS ok
            iFound = -1
            DO iJ = 1,iNUMVibTempMagic
                IF (iaCnt1(iJ)     /= iaCnt(iJ))     iFound = +1
                IF (iaGasID1(iJ)   /= iaGasID(iJ))   iFound = +2
                IF (iaISO1(iJ)     /= iaISO(iJ))     iFound = +3
                IF (iaAFGL1(iJ)    /= iaAFGL(iJ))    iFound = +4
                IF (raVibCntr1(iJ) /= raVibCntr(iJ)) iFound = +5
            END DO
            IF (iFound > 0) THEN
                write(kStdErr,*) 'Error checking out HITRAN GasIDs ',iFound
                write(kStdErr,*) 'in subroutine OutputVTFile'
                CALL DoStop
            END IF
                 
        END IF
    END DO

! interpolate, and then dump out the info!
! first check to see if user SolarAngle == one of the 6 angles we have
    iFound = -1
    DO iI = 1,iNumSolAngles
        IF (abs(rSolarAngle - raSolAngles(iI)) <= 0.1) THEN
            iFound = iI
            GOTO 10
        END IF
    END DO

    10 CONTINUE
    IF (iFound > 0) THEN
    ! no need to interpolate; simply set raaVibPartFcn,raaVibTemp
        DO iJ = 1,iNumVibLevels
            DO iI = 1,iVibPF
                raaVibPartFcn(iI,iJ) =  raaaVibPartFcn(iFound,iI,iJ)
            END DO
        END DO

        DO iJ = 1,iNumVibLevels
            DO iI = 1,iNUMVibTempMagic
                raaVibTemp(iI,iJ) =  raaaVibTemp(iFound,iI,iJ)
            END DO
        END DO

    ELSEIF (iFound < 1) THEN
    ! need to interpolate every layer, every variable, in sunangle
        DO iI = 1,iNumVibLevels !!!scan over the 120 pressure levels
            DO iJ = 1,iVibPf      !!!scan over the different isotopes
                DO iK = 1,MXXTMP    !!!scan over the 6 solar angles
                    ra6(iK) = raaaVibPartFcn(iK,iJ,iI)
                END DO
                CALL rspl_one(raSolAngles,ra6,6,rSolarAngle,rT)
                raaVibPartFcn(iJ,iI) =  rT
            END DO
        END DO

        DO iI = 1,iNumVibLevels       !!!scan over the 120 pressure levels
            DO iJ = 1,iNUMVibTempMagic  !!!scan over the different vib levels
                DO iK = 1,MXXTMP          !!!scan over the 6 solar angles
                    ra6(iK) = raaaVibTemp(iK,iJ,iI)
                END DO
                CALL rspl_one(raSolAngles,ra6,6,rSolarAngle,rT)
                raaVibTemp(iJ,iI) =  rT
            END DO
        END DO
    END IF

!!!dump out the info into caVTFile
    CALL VibTempProfWrite(iAngle,caVTFile,caaComments,iNumComments, &
    raaPressure,raaKinetic,raaVibPartFcn,raaVibTemp, &
    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr, &
    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)

    RETURN
    end SUBROUTINE OutputVTFile

!************************************************************************
! this subroutine writes out the vib temps
! in the format that GENLN2 and kCARTA like
    SUBROUTINE VibTempProfWrite(iAngle,caFName,caaComments,iNumComments, &
    raaPressure,raaKinetic,raaVibPartFcn,raaVibTemp, &
    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr, &
    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! input vars
    CHARACTER(80) :: caFname                          !i/o file name
    INTEGER :: iAngle                                !which angle profile
    REAL :: raaPressure(MXXTMP,kNLTEProfLayer)
    REAL :: raaKinetic(MXXTMP,kNLTEProfLayer)
    REAL :: raaVibPartFcn(kMaxUserSet,kNLTEProfLayer)
    REAL :: raaVibTemp(kMaxK,kNLTEProfLayer)
    CHARACTER(80) :: caaComments(kMaxLayer)
    INTEGER :: iNumComments
    INTEGER :: iGasesInFile,iNumVibLevels            !gases,levels in file
    INTEGER :: iGasID,iVibPF                         !gas ID, number of isotopes
    INTEGER :: iNUMVibTempMagic                      !number of vib profiles
    INTEGER :: iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
    REAL :: raVibCntr(kMaxK)

! local vars
    CHARACTER(3) :: ca3
    CHARACTER(1) :: ca1
    CHARACTER(80) :: caStr
    INTEGER :: i1,i2,iIOUN,iErr,iI,iJ,iDummy1,iDummyGasID,iDummyISO
    INTEGER :: iNV,iZero
    REAL :: raP2(kNLTEProfLayer),raT2(kNLTEProfLayer)

    iZero = 0

    iNV = iNumVibLevels
    DO iI = 1,iNV
        raP2(iI) = raaPressure(1,iI)
        raT2(iI) = raaKinetic(1,iI)
    END DO

    iIOUN = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='UNKNOWN', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write(kStdErr,*) 'error writing VibTemp profile ',iErr,caFName
        CALL DoStop
    END IF
     
    kTempUnitOpen = +1

    caaComments(1) = '! INTERPD (sun angle) vib temp profiles reformatted'
    DO i1 = 1,iNumComments
        caStr = caaComments(i1)
        write(iIOUN,123) caStr
    END DO

!!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV
!!!where FLAG  = ***
!!!      MODEL = atmosphere model number for data file
!!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL
!!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB
!!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL
!!! so read in flag, model,nset,ngas,nlev where flag = '***'
    caStr = '***'
    write(iIOUN,99) caStr,1,1,iGasesInFile,iNumVibLevels
    write(iIOUN,*) iGasID,iVibPF ! how many vib part fcns for the NLTE gas
          
    CALL printarray(iIOUN,iNV,raP2)  !!!print out press levels

    write(iIOUN,*) iZero       !dummy
    CALL printarray(iIOUN,iNV,raT2)  !!!print out kinetic temps
     
    caStr = '!QV   Vibrational Partition functions '
    write(iIOUN,13) caStr
    DO iI = 1,iVibPf
        DO iJ = 1,iNV
            raP2(iJ) = raaVibPartFcn(iI,iJ)
        END DO
        write(iIOUN,*) iGasID,iI
        CALL printarray(iIOUN,iNV,raP2)  !!!print out part fcn
    END DO

    caStr = '!TV   Vibrational Temperatures'
    write(iIOUN,13) caStr
    DO iI = 1,iNUMVibTempMagic
        DO iJ = 1,iNV
            raP2(iJ) = raaVibTemp(iI,iJ)
        END DO
        write(iIOUN,20) iI,iaGasID(iI),iaISO(iI),iaAFGL(iI),raVibCntr(iI)
        CALL printarray(iIOUN,iNV,raP2)
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen = -1

    10 FORMAT(A80)
    11 FORMAT('! ',A70)
    12 FORMAT(A70,I5)
    13 FORMAT(A70)
    20 FORMAT(I4,'    ',I5,'  ',I7,' ',I6,'   ',F11.5)
    99 FORMAT(A3,'   ',4(I3,' '))

    123 FORMAT(A80)

    RETURN
    end SUBROUTINE VibTempProfWrite

!************************************************************************
! this subroutine reads in the vib temps
! in the format that GENLN2 and kCARTA like
    SUBROUTINE VibTempProfRead(iAngle,caFName,caaComments,iNumComments, &
    raaPressure,raaKinetic,raaaVibPartFcn,raaaVibTemp, &
    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr, &
    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! input vars
    CHARACTER(80) :: caFname                          !i/o file name
    INTEGER :: iAngle                                !which angle profile
! output vars
    REAL :: raaPressure(MXXTMP,kNLTEProfLayer)
    REAL :: raaKinetic(MXXTMP,kNLTEProfLayer)
    REAL :: raaaVibPartFcn(MXXTMP,kMaxUserSet,kNLTEProfLayer)
    REAL :: raaaVibTemp(MXXTMP,kMaxK,kNLTEProfLayer)
    CHARACTER(80) :: caaComments(kMaxLayer)
    INTEGER :: iNumComments
    INTEGER :: iGasesInFile,iNumVibLevels            !gases,levels in file
    INTEGER :: iGasID,iVibPF                         !gas ID, number of isotopes
    INTEGER :: iNUMVibTempMagic                      !number of vib profiles
!!!next few are the HITRAN id's
    INTEGER :: iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
    REAL :: raVibCntr(kMaxK)
! local vars
    CHARACTER(3) :: ca3
    CHARACTER(1) :: ca1
    CHARACTER(80) :: caStr
    INTEGER :: i1,i2,iIOUN,iErr,iI,iJ,iDummy1,iDummyGasID,iDummyISO

    iIOUN = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write(kStdErr,*) 'error reading VibTemp profile ',iErr,caFName
        CALL DoStop
    END IF
     
    kTempUnitOpen = +1

    iNumComments = 0
    555 CONTINUE
    read(iIOUN,123) caStr
    IF (caStr(1:1) == '!') THEN
        iNumComments = iNumComments + 1
        caaComments(iNumComments) = caStr
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSEIF (caStr(1:3) == '***') THEN
    !!!aha .. now we are cooking
    !!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV
    !!!where FLAG  = ***
    !!!      MODEL = atmosphere model number for data file
    !!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL
    !!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB
    !!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL
    !!! so read in flag, model,nset,ngas,nlev where flag = '***'
        read(caStr,456) ca3,i1,i2,iGasesInFile,iNumVibLevels
    END IF

    caStr = caaComments(iNumComments-1)
    read(caStr,789) ca1,iNUMVibTempMagic

    IF (iGasesInFile /= 1) THEN
        write(kStdErr,*) 'kCARTA can only handle one NLTE gas (CO2) here!!!!'
        CALL DoStop
    END IF
    IF (iNumVibLevels > kNLTEProfLayer) THEN
        write(kStdErr,*) 'kCARTA can handle max',kNLTEProfLayer,' NLTE levels!'
        CALL DoStop
    END IF

    READ(iIOUN,*) iGasID,iVibPF
    IF (iGasID /= 2) THEN
        write(kStdErr,*) 'kCARTA can only handle CO2 for NLTE here!!!!'
        CALL DoStop
    END IF

    IF (iVibPF > kMaxK) THEN
        write(kStdErr,*) 'kCARTA can only handle max ',kMaxK,'Vib part fcns!!'
        CALL DoStop
    END IF

!!!read the presssure levels
    read(iIOUN,*) (raaPressure(iAngle,iI),iI=1,iNumVibLevels)   !!! in mb
     
!!! read the kinetic temperatures
    read(iIOUN,*) iDummy1
    read(iIOUN,*) (raaKinetic(iAngle,iI),iI=1,iNumVibLevels)   !!!LTE temp
     
!!!read in the QTIPS_VIB modifiers
!!!dummy string that says "!QV   Vibrational Partition functions"
    read(iIOUN,123) caStr
    DO iJ = 1,iVibPF
        read(iIOUN,*) iDummyGasID,iDummyISO
        read(iIOUN,*) (raaaVibPartFcn(iAngle,iJ,iI),iI=1,iNumVibLevels)
    END DO

!!! read the NLTE temperatures
!!!dummy string that says "!TV   Vibrational Temperatures"
    read(iIOUN,123) caStr
    DO iJ = 1,iNUMVibTempMagic
        read(iIOUN,*) iaCnt(iJ),iaGasID(iJ),iaISO(iJ),iaAFGL(iJ),raVibCntr(iJ)
        read(iIOUN,*) (raaaVibTemp(iAngle,iJ,iI),iI=1,iNumVibLevels)
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen = -1

    123 FORMAT(A80)
    456 FORMAT(A3,'   ',4(I3,' '))
    789 FORMAT(A1,'   ',I3)

    RETURN
    end SUBROUTINE VibTempProfRead

!************************************************************************
! this function takes string s1 and adds info to it, based on iProf,iSol
    SUBROUTINE addinfo(caS1,iProf,iSol)
     
    IMPLICIT NONE
     
    CHARACTER(80) :: caS1
    INTEGER :: iProf,iSol
     
    INTEGER :: i1
    CHARACTER(2) :: caString
     
    i1 = 80
    10 CONTINUE
    IF ((caS1(i1:i1) == ' ') .AND. (i1 > 1)) THEN
        i1 = i1 - 1
        GOTO 10
    END IF
     
    IF (iSol >= 0) THEN
    ! we are dealing with vib temp files, so we need "vt" names
        caS1(i1+1:i1+2) = 'vt'
        i1 = i1 + 2
    ELSE
    ! we are dealing with profile files, so we need "pt" names
        caS1(i1+1:i1+2) = 'pt'
        i1 = i1 + 2
    END IF
     
! now process iProf so that we end up with a right padded string
! eg 2 ---> '2 ', 12 ---> '12' etc
    WRITE(caString,15) iProf
    IF (iProf >= 10) THEN
        caS1(i1+1:i1+2) = caString(1:2)
        i1 = i1 + 2
    ELSE
        caS1(i1+1:i1+1) = caString(2:2)
        i1 = i1 + 1
    END IF
     
    IF (iSol >= 0) THEN
    ! we are dealing with vib temp files, so we need to put in sol angle
        caS1(i1+1:i1+2) = '_s'
        i1 = i1 + 2
        WRITE(caString,15) iSol
        IF (iSol >= 10) THEN
            caS1(i1+1:i1+2) = caString(1:2)
            i1 = i1 + 2
        ELSE
            caS1(i1+1:i1+1) = caString(2:2)
            i1 = i1 + 1
        END IF
    END IF
     
    caS1(i1+1:i1+4) = '.prf'
     
    15 FORMAT(I2)
     
    RETURN
    end SUBROUTINE addinfo

!************************************************************************
!               LAPACK routines
!************************************************************************

!************************************************************************
!                          POLYNOM ROUTINES
!************************************************************************

!************************************************************************
!                          PREDICTOR ROUTINES
!************************************************************************
! this subroutine finds the closest profile and dumps it into a file
    SUBROUTINE predict_nlte( &
    caOutName,raTemp,raPressLevels,iNumLayers, &
    rSolarAngle,iRTP,iRegr,caVTFile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
    REAL :: raTemp(kProfLayer)           !!!kinetic temp profile
    REAL :: raPressLevels(kProfLayer+1)  !!!pressure levels
    INTEGER :: iNumLayers                !!!number of layers
    REAL :: rSolarAngle                  !!!solar angle for atm #1
    CHARACTER(80) :: caOutName
    INTEGER :: iRTP                      !!!which of the rtp files processed
! output
    INTEGER :: iRegr               !!!which regr profile this is closest to
    CHARACTER(80) :: caVTFile       !!!temp file name where VT profiles are
!!!created and stored if not given by user

! local vars
    REAL :: raLayPress(kProfLayer),raLayPress3(kProfLayerP3)
    REAL :: raTemp3(kProfLayerP3),raAirsLevels(kNLTEProfLayer)
    REAL :: rMin,rA,rB,rC,rPavg
    REAL :: raAirsPressure(kMaxLayer),raAirsTemps(kNLTEProfLayer)
    CHARACTER(4) :: caVT(10)
    CHARACTER(4) :: caQV(4)

    INTEGER :: iI,iJ,iStart,iNumPts3,dI,iAirsMax

    caVT(1) = 'sig1'
    caVT(2) = 'sig2'
    caVT(3) = 'sig3'
    caVT(4) = 'sig4'
    caVT(5) = 'wow1'
    caVT(6) = 'wow2'
    caVT(7) = 'pii1'
    caVT(8) = 'pii2'
    caVT(9) = 'del1'
    caVT(10)= 'del2'

    caQV(1) = 'qqv1'
    caQV(2) = 'qqv2'
    caQV(3) = 'qqv3'
    caQV(4) = 'qqv4'

    iAirsMax = 102    !!!beyond this, airs_levels = beserk

    CALL airs_levels(raPressLevels,raAirsLevels)
!      CALL closest_regr_lowertrop(raTemp,raPressLevels,iNumLayers,iRegr)

!!! this is the pressure layering for the AIRS JPL instrument
    DO iI = 1,kMaxLayer
        raAirsPressure(iI) = raAirsLevels(iI)-raAirsLevels(iI+1)
        raAirsPressure(iI) = raAirsPressure(iI)/ &
        log(raAirsLevels(iI)/raAirsLevels(iI+1))
    END DO

!!!find the pressure layering for the kCARTA profile
    DO iI = kProfLayer-iNumLayers+1,kProfLayer
        raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1)
        raLayPress(iI) = raLayPress(iI)/ &
        log(raPressLevels(iI)/raPressLevels(iI+1))
    END DO

!!!fill lower layers with some "realistic" increasing values
    DO iI = kProfLayer-iNumLayers,1,-1
       raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) + 10*abs(iI-(kProfLayer-iNumLayers+1))
    END DO

!!!interpolate the kintic temps (logarithmically) onto the AIRS levels
    CALL logrspl(raLayPress,raTemp,kProfLayer,raAirsLevels,raAirsTemps,iAirsMax)
!    DO iI = 1,iAirsMax
!        print *,iI,raAirsLevels(iI),raAirsTemps(iI)
!    END DO
            
    CALL DoStop

    CALL VTName_rtp(iRTP,caOutName,caVTFile)
!      CALL OutputVTFile_Poly(caVTFile,raaCompute,raLayPress3,raTemp3,iNumPts3)

    RETURN
    end SUBROUTINE predict_nlte

!************************************************************************
!                          PREDICTOR ROUTINES
!************************************************************************
END MODULE kpredictVT
