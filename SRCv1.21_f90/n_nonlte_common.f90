! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_nonlte_common

USE basic_common
!USE kpredictVT
USE kcoeff_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile
USE s_writefile

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine finds the closest profile, uses results of polynom fit
! (first or second or third order)  and dumps it into a file
    SUBROUTINE polynom_nlte( &
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
    INTEGER :: iRTP                      !!!which iRTP prof read in
! output
    INTEGER :: iRegr               !!!which regr profile this is closest to
    CHARACTER(80) :: caVTFile       !!!temp file name where VT profiles are
!!!created and stored if not given by user

! local vars
    REAL :: raLayPress(kProfLayer),raLayPress3(kProfLayerP3)
    REAL :: raTemp3(kProfLayerP3),raAirsLevels(kNLTEProfLayer)
    REAL :: rMin,rA,rB,rC,rPavg

    INTEGER :: iI,iJ,iStart,iNumPts3,dI

    CALL airs_levels(raPressLevels,raAirsLevels)

!!!first find the pressure layering for the kCARTA profile
    DO iI = kProfLayer-iNumLayers+1,kProfLayer
      raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1)
      raLayPress(iI) = raLayPress(iI)/ &
      log(raPressLevels(iI)/raPressLevels(iI+1))
    END DO

!!!fill lower layers with some "realistic" increasing values
    DO iI = kProfLayer-iNumLayers,1,-1
      raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) + &
        &                    10*abs(iI-(kProfLayer-iNumLayers+1))
    END DO

    iStart = 90
 10 CONTINUE
    IF (raAirsLevels(iStart) < raLayPress(kProfLayer)) THEN
      GOTO 20
    ELSEIF (iStart < 102) THEN
      iStart = iStart + 1
      GOTO 10
    ELSE
      write(kStdErr,*) 'humph : cannot find'
      write(kStdErr,*) 'raAirsLevels(iStart) < raLayPress(kProfLayer)'
      Call DoStop
    END IF
 20 CONTINUE

!!!fill in the arrays for extended layers
    DO iI = 1,kProfLayer
      raLayPress3(iI) = raLayPress(iI)
      raTemp3(iI)     = raTemp(iI)
    END DO
!!!so now can go slightly above the AIRS levels
! SUBROUTINE quadratic_coeffs(raX,raY,iMidPt,iLogOrLinear,A,B,C)
    CALL quadratic_coeffs(raLayPress3,raTemp3,kProfLayer-1,+1,rA,rB,rC)

    DO iI = 1,103-iStart+1
      raLayPress3(kProfLayer+iI) = raAirsLevels(iStart-1+iI)
      ! guess the kinetic temps here!
      rPavg = raLayPress3(kProfLayer+iI)
      raTemp3(kProfLayer+iI) = rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
      iNumPts3 = kProfLayer+iI
    END DO

!!!flip raLayPress3 so it is in increasing order!
!!!ditto raTemp3
    DO iI = 1,int(iNumPts3*1.0/2.0)
      rMin = raLayPress3(iI)
      raLayPress3(iI) = raLayPress3(iNumPts3-iI+1)
      raLayPress3(iNumPts3-iI+1) = rMin

      rMin = raTemp3(iI)
      raTemp3(iI) = raTemp3(iNumPts3-iI+1)
      raTemp3(iNumPts3-iI+1) = rMin
    END DO

!!!also flip raLayPress so it is in increasing order!
!!!ditto raTemp
!!!this is really NOT necessary as we DO NOT use the arrays here
!!!but whatever, we might use them in n_gas_wt_spectra
    DO iI = 1,int(kProfLayer*1.0/2.0)
      rMin = raLayPress(iI)
      raLayPress(iI) = raLayPress(kProfLayer-iI+1)
      raLayPress(kProfLayer-iI+1) = rMin

      rMin = raTemp(iI)
      raTemp(iI) = raTemp(kProfLayer-iI+1)
      raTemp(kProfLayer-iI+1) = rMin
    END DO

    CALL NLTE_PolyTemp(caOutName,raLayPress3,raTemp3,iNumPts3, &
      rSolarAngle,caVTFile,iRTP,iRegr)

    RETURN
    end SUBROUTINE polynom_nlte

!************************************************************************
! this subroutine reads in the line parameters for the band in question
    SUBROUTINE read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
    iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
    daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR, &
    iLineMixBand,iDoVoigtChi)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params, read from caaNLTETemp(iLTEin)
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iLTEIn,iBand,iDOVoigtChi
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
! output params
    INTEGER :: iGasID,iNum,iISO,iLineMixBand
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)

    INTEGER :: iI,iErr,iIOUN,iJU,iJL
    CHARACTER(80) :: caFname
    DOUBLE PRECISION :: dMax,dL,dR,dC

    caFName = caaaNLTEBands(iLTEin,iBand)

    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='unformatted',STATUS='OLD',IOSTAT=iErr)
    IF (iErr /= 0) THEN
      write (kStdErr,*) 'in subroutine read_lineparameters, error reading'
      write (kStdErr,*) 'file that has HITRAN lineparameters'
      WRITE(kStdErr,1070) iErr, caFName
      CALL DoSTOP
    END IF
    kTempUnitOpen = 1

    read(iIOUN) iGasID,iNum,iISO
    write(kstdWarn,*) ' '
    write(kStdWarn,*) '----> Opening HITRAN parameter file : '
    write(kStdWarn,1080) caFName
    write(kStdWarn,*) 'File has ',iNum,' line parameters for gas ',iGasID

    IF (iNum > kHITRAN) THEN
      write(kStdErr,*) 'File has ',iNum,' line parameters for gas ',iGasID
      write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
      write(kStdErr,*) 'Please check kHITRAN in kcartaparam.f90 and fix'
      CALL DoStop
    END IF

    read(iIOUN) (daElower(iI),iI = 1,iNum)     !lower state energy
    read(iIOUN) (daLineCenter(iI),iI = 1,iNum) !line center
    read(iIOUN) (daJL(iI),iI = 1,iNum)         !lower vib quantum number .. basically SAME for all iI=1,iNum
    read(iIOUN) (daJU(iI),iI = 1,iNum)         !upper vib quantum number .. basically SAME for all iI=1,iNum
    read(iIOUN) (daPshift(iI),iI = 1,iNum)     !pressure shift of linecenter
    read(iIOUN) (daStren296(iI),iI = 1,iNum)   !line strenght
    read(iIOUN) (daW_for(iI),iI = 1,iNum)      !foreign broadening/atm
    read(iIOUN) (daW_self(iI),iI = 1,iNum)     !self broadening/atm
    read(iIOUN) (daW_temp(iI),iI = 1,iNum)     !broadening tempr dependance
!! new since July 2015, comes from line.bslq (see lineparameters.m in
!! SRCv1.18/NONLTE/M_Files_for_kcarta_NLTE_LBL_runs/USUALLAYERS/
    read(iIOUN) (daJLowerQuantumRot(iI),iI = 1,iNum)  !J lower quantum rotation state number

    close (iIOUN)
    kTempUnitOpen = -1
          
    iJU = nint(daJU(1))
    iJL = nint(daJL(1))

    1070 FORMAT('ERROR! number ',I5,' opening HITRAN parameter file:',/,A80)
    1080 FORMAT(A80)

! outside of the weaklines, do linemixing for all the bands the user specifies!
! but be careful about Cousin vs linemix, else code becomes VERY slow
    IF (iGasID == 2) THEN
      IF ((iJL == 4) .AND. (iJU == 24) .AND. (iISO == 1)) THEN
        iLineMixBand = +2
        write(kStdWarn,*) 'very strong CO2 band : deldel (2310) iLineMixBand = +2'
      ELSEIF ((iJL == 2) .AND. (iJU == 16) .AND. (iISO == 1)) THEN
        iLineMixBand = +2
        write(kStdWarn,*) 'very strong CO2 band : pi pi  (2320) iLineMixBand = +2'
      ELSEIF ((iJL == 1) .AND. (iJU == 9) .AND. (iISO == 1)) THEN
        iLineMixBand = +2
        write(kStdWarn,*) 'very strong CO2 band : sigsig (2350) iLineMixBand = +2'
      ELSEIF ((iJL == 1) .AND. (iJU == 9) .AND. (iISO == 2)) THEN
        iLineMixBand = +2
        write(kStdWarn,*) 'very strong CO2 band : sigsig (2351) iLineMixBand = +2'
      ELSE
        write(kStdWarn,*) 'strong CO2 bands : iLineMixBand = +1 (cousin)'
        iLineMixBand = +1
      END IF
    ELSE
      iLineMixBand = -1
      write(kStdWarn,*) 'not a CO2  band : no linemixing (voigt) : iLineMixBand = -1'
    END IF

    IF ((iLineMixBand == 2) .AND. (iJU /= 9)) THEN
      !!!only do linemix for 2350, 2351 lines; else do Cousin
      iLineMixBand = +1
      write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
      write(kStdWarn,*) '   reset iLineMixBand = from +2 to +1 (cousin)'
    ELSEIF ((iLineMixBand == 2) .AND. (iJU == 9)) THEN
      IF (iISO <= 2) THEN
        write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
        write(kStdWarn,*) '   iLineMixBand = +2 (linemix)'
      ELSEIF (iISO > 2) THEN
        write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
        write(kStdWarn,*) '   iLineMixBand = +1 (cousin)'
      END IF
    END IF
            
! this always makes cousin the happening one!
!       IF (iLineMixBand .GT. 0) THEN
!         iLineMixBand = +1
!         write(kStdWarn,*) 'doing cousin everywhere; reset iLineMixBand = +1'
!         print *,' ***** reset iLineMixBand = 1   ie linemix --> cousin ****'
!       END IF

!      IF (iLineMixBand .GT. 0) THEN
!        IF (iDoVoigtChi .GT. 0) THEN
!          write(kStdErr,*) 'Cannot have idoVoigtChi > 0 '
!          write(kStdErr,*) 'AND try to use Cousin everywhere!'
!          CALL DoStop
!        END IF
!        print *,' ***** reset iLineMixBand = 1   ie linemix --> cousin ****'
!        iLineMixBand = 1
!        END IF

    RETURN
    end SUBROUTINE read_lineparameters

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
    REAL :: p,pp,q,t
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
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr,ERR = 100)
    IF (iErr /= 0) THEN
     !this does not get executed as error messgs get handled by line 100
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

    IF ((iMatchTry == 1) .AND. (iaJ_UorL(1) == 0) .AND. (iAllOrSOme > 0)) THEN
      write(kStdWarn,*) 'Did not find a match for IUSGQ; try ILSGQ ...'
      GOTO 111
    END IF

    IF ((iMatchTry == 2) .AND. (iaJ_UorL(1) == 0) .AND. (iAllOrSome > 0)) THEN
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
    END IF

    IF ((iAllORSome > 0) .AND. (dVibCenter < 1.0d-2)) THEN
      write (kStdErr,*) 'need to set dVibCenter, but seem not to have'
      write (kStdErr,*) 'been able to do so!'
      Call DoStop
    END IF

    CALL Add005mblevel(iNumVibLevels,raPressVT,raKineticVT,raQtipsVT,raNLTETempVT)

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
!************************************************************************
!                          POLYNOM ROUTINES
!************************************************************************
! this subroutine computes the AIRS pressure levels
! see /asl/packages/klayers/Doc/sci.txt

! Layers and Levels:
! -----------------
! As we use the terms, we mean different things when we speak of "layers" and
! "levels".  A "level" refers to a point value, while a "layer" refers to a
! finite thickness slab value.  We define our layers using levels as layer
! boundaries. For AIRS, the definition is:
!     Plev(i) = exp( (7/2) * ln( A*i^2 + B*i + C ) )

! where
!   i refers to the level number (an integer counter)
!   and A, B, C are constants which obey the relation (exact):
!      Plev(i=1)=1100, Plev(i=38)=300, Plev(i=101)=5.0E-3 mb
!   with approximate values:
!      A = -1.5508E-4
!      B = -5.5937E-2
!      C =  7.4516

! WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
! and so can only be trusted to iI = 103

    SUBROUTINE airs_levels(raPressLevels,raAirsLevels)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: raAirsLevels(kNLTEProfLayer)
    REAL :: raPressLevels(kProfLayerP1)

    INTEGER :: iI
    REAL :: A,B,C
     
    A = -1.5508E-4
    B = -5.5937E-2
    C =  7.4516

! WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
! and so can only be trusted to iI = 103
          
    DO iI = 1,kNLTEProfLayer
      raAirsLevels(iI) = exp( (7.0/2) * log(A*iI*iI + B*iI + C ) )
    END DO

! WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
! and so can only be trusted to iI = 103

    RETURN
    end SUBROUTINE airs_levels

!************************************************************************
! this finds quadratic coeffs a b c   for y = ax2 + bx + c
    SUBROUTINE quadratic_coeffs(raX,raY,iMidPt,iLogOrLinear,rA,rB,rC)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL :: raX(*),raY(*)   !!arrays to fit
    INTEGER :: iMidPt       !!fit about this point
    INTEGER :: iLogOrLinear !! +1 for raX to be log, -1 to be linear
! output params
    REAL :: rA,rB,rC

! local vars
    REAL :: rP0,rPp1,rPm1
    REAL :: rT0,rTp1,rTm1
    REAL :: rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr

    rPp1 = raX(iMidPt+1)
    rP0  = raX(iMidPt)
    rPm1 = raX(iMidPt-1)

    rTp1 = raY(iMidPt+1)
    rT0  = raY(iMidPt)
    rTm1 = raY(iMidPt-1)

! now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n))
    IF (iLogOrLinear == +1) THEN
      rP0  = alog(rP0)
      rPp1 = alog(rPp1)
      rPm1 = alog(rPm1)
    END IF
            
    rDp1 = rTp1-rT0
    rDm1 = rTm1-rT0
     
    rp1    = rPp1-rP0
    rp1sqr = (rPp1-rP0)*(rPp1+rP0)
    rm1    = rPm1-rP0
    rm1sqr = (rPm1-rP0)*(rPm1+rP0)
     
    rA = (rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1)
    rB = rDp1/rp1-rA*(rp1sqr/rp1)
    rC = rT0-rA*rP0*rP0-rB*rP0
     
! finally compute rT
!        rT=rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC

    RETURN
    end SUBROUTINE quadratic_coeffs

!************************************************************************
! this subroutine reads in the generic files and computes a bunch of things!
! based on imput (raLayPress,raKinetic)
    SUBROUTINE NLTE_PolyTemp(caOutName,raLayPress3,raTemp3,iNumPts3, &
    rSolarAngle,caVTFile,iRTP,iRegr)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    CHARACTER(80) :: caOutName
    REAL :: raTemp3(kProfLayerP3)           !!!kinetic temp profile
    REAL :: raLayPress3(kProfLayerP3)       !!!layering
    INTEGER :: iNumPts3                     !!!number of points in the arrays
    REAL :: rSolarAngle
    INTEGER :: iRegr                        !!! closest regr prof number (1:48)
    INTEGER :: iRTP                         !!!which iRTP prof read in
! output param
    CHARACTER(80) ::  caVTfile

! local vars
    CHARACTER(120) :: caDir,caFnamePoly
    CHARACTER(11) ::  caExt(12)
    REAL :: raaaCoef(kNLTEProfLayer,7,4)
    REAL :: raSolar(7),raPressRegr(kNLTEProfLayer)
    REAL :: raaGood(kNLTEProfLayer,7)

    REAL :: raaCompute(kNLTEProfLayer,12),raOut(kNLTEProfLayer)

    INTEGER :: iI,iOrder0,iNumSolar,iNumPressPts,iCompute

    caExt(1) = 'coeffs_ss1_'
    caExt(2) = 'coeffs_ss2_'
    caExt(3) = 'coeffs_ss3_'
    caExt(4) = 'coeffs_ss4_'
    caExt(5) = 'coeffs_pp1_'
    caExt(6) = 'coeffs_pp2_'
    caExt(7) = 'coeffs_dd1_'
    caExt(8) = 'coeffs_dd2_'
    caExt(9) = 'coeffs_qv1_'
    caExt(10)= 'coeffs_qv2_'
    caExt(11)= 'coeffs_qv3_'
    caExt(12)= 'coeffs_qv4_'

    iOrder0 = 2

    caDir = kNLTEPolyFits

! stringchar = ['ss1','ss2','ss3','ss4','pp1','pp2','dd1','dd2',...
!              'qv1','qv2','qv3','qv4'];

    DO iCompute = 1,12
      CALL makeVT_polyname(iOrder0,iCompute,caExt,caDir,caFnamePoly)
      CALL read_VTpoly(iOrder0,caFnamePoly,iNumPressPts,iNumSolar, &
        raSolar,raPressRegr,raaaCoef,raaGood)
      CALL Compute_NLTE_Poly( &
        iNumPressPts,iNumSolar, &
        raSolar,raPressRegr,raaaCoef,raaGood, &
        iCompute,iOrder0,raLayPress3,raTemp3,iNumPts3,rSolarAngle, &
        raOut)
      DO iI = 1,iNumPts3
        raaCompute(iI,iCompute) = raOut(iI)
      END DO
    END DO

    CALL VTName_rtp(iRTP,caOutName,caVTFile)
    CALL OutputVTFile_Poly(caVTFile,raaCompute,raLayPress3,raTemp3,iNumPts3)

    RETURN
    end SUBROUTINE NLTE_PolyTemp

!************************************************************************
! just dumps stuff to a file
    SUBROUTINE OutputVTFile_Poly(caVTFile,raaCompute, &
    raLayPress3,raTemp3,iNumPts3)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    CHARACTER(80) :: caVTfile,caSummaryFile
    REAL :: raaCompute(kNLTEProfLayer,12)
    REAL :: raTemp3(kProfLayerP3)           !!!kinetic temp profile
    REAL :: raLayPress3(kProfLayerP3)       !!!layering
    INTEGER :: iNumPts3
    CHARACTER(3) :: ca3

! local vars
    INTEGER :: iI,iJ,iIOUN,iVibPf,iGasID,iErr,iZero
    INTEGER :: iaGasID(8),iaISO(8),iaAFGL(8)
    REAL :: raVibCntr(8)
    REAL :: raX(kProfLayerP3)
    CHARACTER(70) :: caStr

    iVibPf = 4
    iGasID = 2
    ca3 = '***'

    DATA (iaGasID(iI),iI=1,8) /2,2,2,2,2,2,2,2/
    DATA (iaISO(iI),iI=1,8)   /1,2,3,4,1,2,1,2/
    DATA (iaAFGL(iI),iI=1,8)/ &
    &               9,          9,          9,          9, &
    &               16,         16,         24,         24/
    DATA (raVibCntr(iI),iI=1,8)/ &
    &              2349.143300, 2283.48800, 2332.11300, 2340.01400, &
    &              3004.01200,  2920.23900, 3659.27300, 3580.75000/

    iI = 80
 100 CONTINUE
    IF (caVTfile(iI:iI) == ' ') THEN
      iI = iI - 1
      GOTO 100
    END IF
    caSummaryFile(1:iI) = caVTfile(1:iI)
    caSummaryFile(iI+1:iI+4) = '.sss'
    iIOUN = kTempUnit
! dump out stuff so it looks like the .sss files made from 48 regr profs
    OPEN(UNIT=iIOun,FILE=caSummaryFile,FORM='formatted',STATUS='UNKNOWN',IOSTAT=iErr)
    IF (iErr /= 0) THEN
      write(kStdErr,*) 'trying to dump NLTE summary '
      write(kStdErr,*) 'error opening polynom NLTE file ',iErr,caSummaryFile
      CALL DoStop
    ELSE
      write(kStdWarn,*) 'dumping summary NLTE guesses to ',caVTFile
      kTempUnitOpen = +1
      DO iI = 1,iNumPts3
        write(iIOUN,111) iI,raLayPress3(iI),raTemp3(iI), &
        raaCompute(iI,9),raaCompute(iI,10),raaCompute(iI,11),raaCompute(iI,12), &
        raaCompute(iI,1),raaCompute(iI,2),raaCompute(iI,3),raaCompute(iI,4), &
        raaCompute(iI,5),raaCompute(iI,6),raaCompute(iI,7),raaCompute(iI,8)
      END DO
      CLOSE(iIOUN)
      kTempUnitOpen = +1
    END IF

    iIOUN = kTempUnit
    OPEN(UNIT=iIOun,FILE=caVTFile,FORM='formatted',STATUS='UNKNOWN', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
      write(kStdErr,*) 'trying to dump NLTE guesses '
      write(kStdErr,*) 'error opening NLTE guesses file ',iErr,caVTFile
      CALL DoStop
    ELSE
        write(kStdWarn,*) 'dumping NLTE guesses to ',caVTFile
    END IF
    kTempUnitOpen = +1

!! ca3  ModelNumber   NumOfVibStates in file NumOfGases  NumLevs
    write(iIOUN,99) ca3,1,8,1,iNumPts3
!! GasID NumOfVibPartFcn
    write(iIOUN,*) 2,4
    DO iI = 1,iNumPts3
      raX(iI) = raLayPress3(iNumPts3-iI+1)   !!!flip
    END DO
    CALL printarray(iIOUN,iNumPts3,raX)  !!!print out press levels

    write(iIOUN,*) iZero       !dummy
    DO iI = 1,iNumPts3
      raX(iI) = raTemp3(iNumPts3-iI+1)   !!!flip
    END DO
    CALL printarray(iIOUN,iNumPts3,raX)  !!!print out kinetic temps
     
    caStr = '!QV   Vibrational Partition functions '
    write(iIOUN,13) caStr
    DO iI = 1,iVibPf
      DO iJ = 1,iNumPts3
        raX(iJ) = raaCompute(iNumPts3-iJ+1,9+(iI-1))
      END DO
      write(iIOUN,*) iGasID,iI
      CALL printarray(iIOUN,iNumPts3,raX)  !!!print out part fcn
    END DO

    caStr = '!TV   Vibrational Temperatures'
    write(iIOUN,13) caStr
    DO iI = 1,8
      DO iJ = 1,iNumPts3
        raX(iJ) = raaCompute(iNumPts3-iJ+1,iI)
      END DO
      write(iIOUN,20) iI,iaGasID(iI),iaISO(iI),iaAFGL(iI),raVibCntr(iI)
      CALL printarray(iIOUN,iNumPts3,raX)  !!!print out NLTE VibTemp
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen = -1

 10 FORMAT(A80)
 11 FORMAT('! ',A70)
 12 FORMAT(A70,I5)
 13 FORMAT(A70)
 20 FORMAT(I4,'    ',I5,'  ',I7,' ',I6,'   ',F11.5)
 99 FORMAT(A3,'   ',4(I3,' '))
 111 FORMAT(I4,' ',14(F12.5,' '))
 123 FORMAT(A80)
     
    RETURN
    end SUBROUTINE OutputVTFile_Poly

!************************************************************************
! this makes the name of the polynomial fit to read in
    SUBROUTINE makeVT_polyname(iOrder0,iCompute,caExt,caDir,caFnamePoly)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iOrder0,iCompute
    CHARACTER(120) :: caDir
    CHARACTER(11) :: caExt(12)
! output params
    CHARACTER(120) :: caFnamePoly

    INTEGER :: iI,iLen
    CHARACTER(11) :: cExt

    DO iI = 1,120
      caFnamePoly(iI:iI) = ' '
    END DO

    cExt = caExt(iCompute)
          
    iLen = 120
 10 CONTINUE
    IF ((caDir(iLen:iLen) == ' ') .AND. (iLen > 1)) THEN
      iLen = iLen - 1
      GOTO 10
    END IF

    DO iI = 1,iLen
      caFnamePoly(iI:iI) = caDir(iI:iI)
    END DO
          
    DO iI = 1,11
      caFnamePoly(iLen+iI:iLen+iI) = cExt(iI:iI)
    END DO

    iLen = 120
 20 CONTINUE
    IF ((caFnamePoly(iLen:iLen) == ' ') .AND. (iLen > 1)) THEN
      iLen = iLen - 1
      GOTO 20
    END IF

    IF (iOrder0 == 2) THEN
      caFnamePoly(iLen+1:iLen+5) = '2.txt'
    ELSEIF (iOrder0 == 2) THEN
      caFnamePoly(iLen+1:iLen+5) = '3.txt'
    ELSE
      write(kStdErr,*) 'hmm incorrect iOrder0'
      CALL DoStop
    END IF

    RETURN
    end SUBROUTINE makeVT_polyname

!************************************************************************
! this subroutine computes the thingies
    SUBROUTINE Compute_NLTE_Poly(iNumPressPts,iNumSolar, &
    raSolar,raPressRegr,raaaCoef,raaGood, &
    iCompute,iOrder0,raLayPress3,raTemp3,iNumPts3,rSolarAngle, &
    raOut)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input from the polyfit files
    INTEGER :: iNumPressPts,iNumSolar,iCompute
    REAL :: raaaCoef(kNLTEProfLayer,7,4)
    REAL :: raSolar(7),raPressRegr(kNLTEProfLayer)
    REAL :: raaGood(kNLTEProfLayer,7)
! input from user profile
    INTEGER :: iOrder0
    REAL :: raTemp3(kProfLayerP3)           !!!kinetic temp profile
    REAL :: raLayPress3(kProfLayerP3)       !!!layering
    INTEGER :: iNumPts3
    REAL :: rSolarAngle
! output
    REAL :: raOut(kNLTEProfLayer)

! local vars
    INTEGER :: iI,iA,iB,iP,iX
    REAL :: raY1(kNLTEProfLayer),raY2(kNLTEProfLayer),raX(kNLTEProfLayer)
    REAL :: raaP1(kNLTEProfLayer,4),raaP2(kNLTEProfLayer,4)
    REAL :: raZin(kProfLayerP3),raZout(kProfLayerP3),raOut2(kNLTEProfLayer)
    REAL :: rSlope

    IF ((rSolarAngle < 0.0) .OR. (rSolarAngle > 120.0)) THEN
      write(kStdWarn,*) 'Hmmm solar ang == ',rSolarAngle
      write(kStdWarn,*) 'so setting raOut = 1.0 or to raTemp3'
      IF ((iCompute >= 9) .AND. (iCompute <= 12)) THEN
        DO iI = 1,kNLTEProfLayer
          raOut(iI) = 1.0          !!! vib part fcn
        END DO
      ELSEIF ((iCompute >= 1) .AND. (iCompute <= 8)) THEN
        DO iI = 1,iNumPts3
          raOut(iI) = raTemp3(iI)          !!! kinetic temp
        END DO
        DO iI = iNumPts3+1,kNLTEProfLayer
          raOut(iI) = raOut(iNumPts3)   !!! kinetic temp
        END DO
        GOTO 50
      END IF
    END IF

    IF ((rSolarAngle >= 0.0) .AND. (rSolarAngle <= 0.1)) THEN
      iA = 1
      iB = 1

    ELSEIF ((rSolarAngle >= 119.9) .AND. (rSolarAngle <= 120.0)) THEN
      iA = 7
      iB = 7

    ELSEIF ((rSolarAngle > 0.1) .AND. (rSolarAngle < 119.9)) THEN
      write(kStdWarn,*) 'solar ang == ',rSolarAngle
      write(kStdWarn,*) 'so computing raOut'
      iA = 7
 10   CONTINUE
      IF ((rSolarAngle < raSolar(iA)) .AND. (iA > 1)) THEN
        iA = iA - 1
        GOTO 10
      END IF

      iB = 1
 20   CONTINUE
      IF ((rSolarAngle > raSolar(iB)) .AND. (iB < 7)) THEN
        iB = iB + 1
        GOTO 20
      END IF
    END IF

!!!!!!!!!!!!!!! have figured out iA,iB so proceed

    IF (iA > iB) THEN
      write(kstdErr,*) 'Ooooooops! iA > iB'
      CALL DOStop
    END IF

    write(kStdWarn,*) 'solang bracket ',raSolar(iA),rSolarAngle,raSolar(iB)

    DO iI = 1,iNumPts3
      raZin(iI) = log10(raLayPress3(iI))
    END DO

    IF (iOrder0 == 2) THEN
      iX = 2
      DO iI = 1,kNLTEProfLayer
        raaP1(iI,iX-1) = 0.0
        raaP2(iI,iX-1) = 0.0
      END DO
    ELSE
      iX = 1
    END IF

    IF (iA == iB) THEN
      !! easy; no need to interp!
      DO iP = iX,4
        !!! put in the interp pts
        DO iI = 1,iNumPressPts
          raX(iI)  = log10(raPressRegr(iI))
          raY1(iI) = raaaCoef(iI,iA,iP)
        END DO
        CALL rspl(raX,raY1,iNumPressPts,raZin,raZout,iNumPts3)
        DO iI = 1,iNumPts3
          raaP1(iI,iP) = raZout(iI)
        END DO
      END DO

      DO iI = 1,iNumPts3
        raOut(iI) =  raaP1(iI,1)*(raTemp3(iI)**3) + &
        raaP1(iI,2)*(raTemp3(iI)**2) + &
        raaP1(iI,3)*raTemp3(iI) + raaP1(iI,4)
      END DO
      GOTO 40
    END IF

    IF (iA /= iB) THEN
      !! need to interp in solar ang!
      DO iP = iX,4
       !!! put in the interp pts
       DO iI = 1,iNumPressPts
         raX(iI)  = log10(raPressRegr(iI))
         raY1(iI) = raaaCoef(iI,iA,iP)
         raY2(iI) = raaaCoef(iI,iB,iP)
       END DO
       CALL rspl(raX,raY1,iNumPressPts,raZin,raZout,iNumPts3)
       DO iI = 1,iNumPts3
         raaP1(iI,iP) = raZout(iI)
       END DO
       CALL rspl(raX,raY2,iNumPressPts,raZin,raZout,iNumPts3)
       DO iI = 1,iNumPts3
         raaP2(iI,iP) = raZout(iI)
       END DO
     END DO

      DO iI = 1,iNumPts3
        raOut(iI) =  raaP1(iI,1)*(raTemp3(iI)**3) + &
        raaP1(iI,2)*(raTemp3(iI)**2) + &
        raaP1(iI,3)*raTemp3(iI) + raaP1(iI,4)
      END DO
      DO iI = 1,iNumPts3
        raOut2(iI) =  raaP2(iI,1)*(raTemp3(iI)**3) + &
        raaP2(iI,2)*(raTemp3(iI)**2) + &
        raaP2(iI,3)*raTemp3(iI) + raaP2(iI,4)
      END DO

      !! now do linear interp in solar angle
      DO iI = 1,iNumPts3
        rSlope = (raOut2(iI) - raOut(iI))/(raSolar(iB)-raSolar(iA))
        raOut(iI) = rslope*(rSolarAngle-raSolar(iA)) + raOut(iI)
      END DO
    END IF

 40 CONTINUE

    IF ((iCompute >= 1) .AND. (iCompute <= 8)) THEN
      !!! remember we did a polyfit of (tK vs (dt = tK - tWhatever))
      !!! so tK - dt = tWhatever
      !!! remember we did a polyfit of (tK vs (dt = tWhatever - tK))
      !!! so tK + dt = tWhatever
      DO iI = 1,iNumPts3
        IF (raLayPress3(iI) >= 800) THEN
          ! don't expect NLTE low in the atm
          raOut(iI) = raTemp3(iI)
        ELSE
          raOut(iI) = raTemp3(iI) + raOut(iI)
        END IF
      END DO
    END IF

    IF ((iCompute >= 9) .AND. (iCompute <= 12)) THEN
      !!! remember we did a polyfit of (tK vs QV)
      !!! so pretty much do nothing more
      DO iI = 1,iNumPts3
        IF ((raOut(iI) <= 0.9) .OR. (raOut(iI) >= 1.2) &
           .OR. (raLayPress3(iI) >= 800)) THEN
          ! don't expect NLTE low in the atm
          raOut(iI) = 1.0
        END IF
      END DO
    END IF

 50 CONTINUE

    RETURN
    end SUBROUTINE Compute_NLTE_Poly

!************************************************************************

! this subroutine reads in the coeffs
    SUBROUTINE read_VTpoly(iOrder0,caFnamePoly,iNumPressPts,iNumSolar, &
    raSolar,raPressRegr,raaaCoef,raaGood)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    CHARACTER(120) :: caFnamePoly
    INTEGER :: iOrder0

! output params
    REAL :: raaaCoef(kNLTEProfLayer,7,4)
    REAL :: raSolar(7),raPressRegr(kNLTEProfLayer)
    REAL :: raaGood(kNLTEProfLayer,7)
    INTEGER :: iNumPressPts,iNumSolar

    INTEGER :: iI,iJ,iK,iS,iOrder,iIOUN,iErr
                  
    write(kStdWarn,*) 'Opening file for NLTE polynom eval : '
    write(kStdWarn,*) caFnamePoly

    iIOUN = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFNamePoly,FORM='formatted',STATUS='OLD',IOSTAT=iErr)
    IF (iErr /= 0) THEN
      write(kStdErr,*) 'trying to open polynom NLTE file'
      write(kStdErr,*) 'error opening polynom NLTE file ',iErr,caFNamePoly
      CALL DoStop
    END IF
    kTempUnitOpen = +1

    read(iIOUN,*) iOrder
    IF (iOrder /= iOrder0) THEN
      write(kStdErr,*) 'iOrder,iOrder0 unequal!',iOrder,iOrder0
      CALL DoStop
    END IF

    read(iIOUN,*) iNumPressPts
    IF (iNumPressPts < 1 .OR. iNumPressPts > kNLTEProfLayer) THEN
      write(kStdErr,*) 'Bad iNumPressPts!',iNumPressPts
      CALL DoStop
    END IF

    DO iI = 1,iNumPressPts
      read(iIOUN,*) raPressRegr(iI)
    END DO

    read(iIOUN,*) iNumSolar
    IF (iNumSolar < 1 .OR. iNumSolar > 7) THEN
      write(kStdErr,*) 'Bad iNumSolar!',iNumSolar
      CALL DoStop
    END IF

    DO iI = 1,iNumSolar
      read(iIOUN,*) raSolar(iI)
    END DO

    IF (iOrder0 == 2) THEN
      DO iI = 1,iNumPressPts
        DO iS = 1,iNumSolar
          raaaCoef(iI,iS,1) = 0.0
          read(iIOUN,*) raaGood(iI,iS),raaaCoef(iI,iS,2),raaaCoef(iI,iS,3), &
          raaaCoef(iI,iS,4)
        END DO
      END DO
    ELSEIF (iOrder0 == 3) THEN
      DO iI = 1,iNumPressPts
        DO iS = 1,iNumSolar
          read(iIOUN,*) raaGood(iI,iS),raaaCoef(iI,iS,1),raaaCoef(iI,iS,2), &
          raaaCoef(iI,iS,3),raaaCoef(iI,iS,4)
        END DO
      END DO
    END IF

    close(iIOUN)
    kTempUnitOpen = -1

    RETURN
    end SUBROUTINE read_VTpoly

!************************************************************************
! this subroutine prints out the arrays in a nice format
    SUBROUTINE printarray(iIOUN,iNV,raX)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'
       
    INTEGER :: iIOUN,iNV
!      REAL raX(kNLTEProfLayer)
    REAL :: raX(*)
     
! local vars
    INTEGER :: iM,iD,i1,i2,ia(5),iJ,iFive
    REAL :: ra(5)
     
    iFive = 5
    iM = iimod(iNv,iFive)
    iD = iDiv(iNv,iFive)
     
    ia(1) = -4
    ia(2) = -3
    ia(3) = -2
    ia(4) = -1
    ia(5) =  0
    DO i1 = 1,iD
      DO i2 = 1,5
        ia(i2) = ia(i2) + 5
        ra(i2) = raX(ia(i2))
      END DO
      write(iIOUN,15) (ra(iJ),iJ=1,5)
    END DO
     
    IF ((iM > 0) .AND. (iM < 5)) THEN
      IF (iM == 1) THEN
        ra(1) = raX(iNV)
        write(iIOUN,11) (ra(iJ),iJ=1,iM)
      ELSE IF (iM == 2) THEN
        ra(2) = raX(iNV)
        ra(1) = raX(iNV-1)
        write(iIOUN,12) (ra(iJ),iJ=1,iM)
      ELSE IF (iM == 3) THEN
        ra(3) = raX(iNV)
        ra(2) = raX(iNV-1)
        ra(1) = raX(iNV-2)
        write(iIOUN,13) (ra(iJ),iJ=1,iM)
      ELSE IF (iM == 4) THEN
        ra(4) = raX(iNV)
        ra(3) = raX(iNV-1)
        ra(2) = raX(iNV-2)
        ra(1) = raX(iNV-3)
        write(iIOUN,14) (ra(iJ),iJ=1,iM)
      END IF
    END IF
     
 11 FORMAT(5(E11.5,' '))
 12 FORMAT(5(E11.5,' '))
 13 FORMAT(5(E11.5,' '))
 14 FORMAT(5(E11.5,' '))
 15 FORMAT(5(E11.5,' '))
     
    RETURN
    end SUBROUTINE printarray

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
      write(kStdErr,*) 'usual AIRS layers lowest,highest press level = ',DATABASELEV(1),DATABASELEV(kMaxLayer+1)
      write(kStdErr,*) 'from VT files, lowest,highest press level = ',raTPressX(1),raTPressX(iNumVibLevels)
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
      CALL rspl_one(raLog,raTTempX,   iNumVibLevels,log(rAIRS_TOA),rY1,1)
      CALL rspl_one(raLog,raQtipsX,   iNumVibLevels,log(rAIRS_TOA),rY2,1)
      CALL rspl_one(raLog,raNLTETempX,iNumVibLevels,log(rAIRS_TOA),rY3,1)

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

END MODULE n_nonlte_common
