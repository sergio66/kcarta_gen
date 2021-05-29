! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_mr_common

USE basic_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile

IMPLICIT NONE

CONTAINS

!************************************************************************
! this reads in the reference levels profile
    SUBROUTINE ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iGasID,iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
          
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! output
    INTEGER :: iNumLevsx
    REAL :: raRx110Temp(kMaxLayer+10),raRx110MR(kMaxLayer+10),raRx110Press(kMaxLayer+10)
! local
    REAL :: PLEVx110(kMaxLayer+10)
    REAL :: rLat,raJunk(kMaxLayer+10)
    CHARACTER(80) :: caPFname,caStr,caComment
    INTEGER :: iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec

! first extend PLEV_KCARTADATABASE_AIRS a little above its lwest value, so can do interps well
    DO iI = 1,kMaxLayer+1
      PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(iI)
    END DO
    rLat = PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)/20
    DO iI = kMaxLayer+2,kMaxLayer+10
      PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(kMaxLayer+1) - rLat*(iI-(kMaxLayer+1))
    END DO

    ! caPFname = '/home/sergio/KCARTA/INCLUDE/glatm_16Aug2010.dat'
    caPFname = kcaLevsRefProf

    IF ((kAFGLProf < 1) .OR. (kAFGLProf > 6)) THEN
      write(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
      CALL DoStop
    END IF

    iErr = 0
    iIOUN2 = kProfileUnit
    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO)
    IF (iErrIO /= 0) THEN
      iErr=1
      WRITE(kStdErr,1070) iErrIO, caPfname
      CALL DoSTOP
    ENDIF
 1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
    kProfileUnitOpen=1

    iAFGL = 0
    iFoundGas = -1
          
 10 CONTINUE
    READ(iIOUN2,1080) caStr
    IF ((caStr(1:1) == '!') .OR. (caStr(1:1) == '%')) THEN
      GOTO 10   !! keep reading comments
    ELSE
      READ (caStr,*) iNumLevsx
    END IF
    IF (iNumLevsx > kMaxLayer+1) THEN
      write(kStdErr,*) 'oops iNumLevsx > kMaxLayer+1 ',iNumLevsx,kMaxLayer+1
      CALL DoStop
    END IF

 20 CONTINUE
    iAFGL = iAFGL + 1

 30 CONTINUE
    READ(iIOUN2,1080) caStr
    IF ((caStr(1:1) == '!') .OR. (caStr(1:1) == '%')) THEN
      GOTO 30   !! keep reading comments
    ELSE
      caComment = caStr
    END IF

 40 CONTINUE
    READ(iIOUN2,1080) caStr
    IF ((caStr(1:1) == '!') .OR. (caStr(1:1) == '%')) THEN    
      GOTO 40   !! keep reading comments
    ELSE
      READ (caStr,*) rLat
    END IF
    !write(kStdWarn,*) 'iAFGL Profile ',iAFGL,' rLat ',rLat, ' : ',caComment

    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! Altitude(km)
    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

    IF (iAFGL == kAFGLProf) THEN
      !write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Press(mb)
      READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Temp(K)
      READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)

      !write(kStdWarn,*) '   need to find profile for gasID ',iGasID
    ELSE
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Press(mb)
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Temp(K)
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
    END IF

    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! density (cm-3)
    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

    DO iG = 1,7
      IF (iG == iGasID) THEN
        ! write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
        iFoundGas = +1
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
      ELSE
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
      END IF
    END DO

    IF ((iAFGL == kAFGLProf) .AND. (iFoundGas > 0)) THEN
      GOTO 60   !! found the AFGL prof and found the gas; done!
    ELSEIF ((iAFGL == 6) .AND. (iFoundGas < 0)) THEN
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! modend
      GOTO 50   !! found the AFGL prof but not found the gas
    ELSEIF ((iAFGL < 6) .OR. (iFoundGas < 0)) THEN
      !! either did not find the gas or the AFGL prof
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! modend
      GOTO 20
    END IF

 50 CONTINUE
    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! constituent profs
    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! mingas

    iG = 7
 55 CONTINUE
    iG = iG + 1
    IF (iG == iGasID) THEN
      ! write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
      iFoundGas = +1
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! gas name
      READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
    ELSE
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! gas name
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
    END IF
    IF ((iG < 28) .AND. (iFoundGas < 0)) THEN
      GOTO 55
    END IF
          
    IF (iFoundGas < 0) THEN
      !! bweh try xsec gases
 123  CONTINUE
      READ(iIOUN2,1080) caStr   !! comment
      IF (caStr(1:1) == '!') THEN
        GOTO 123
      ELSEIF (caStr(1:6) == 'DATEND') THEN
        GOTO 60       !!! oh oh end of file, give up!!!!!
      ELSE
        READ (caStr,*) iXsec
        IF (iXsec == iGasID) THEN
          iFoundGas = 1
          READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
        ELSE
          READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
          GOTO 123
        END IF
      END IF
    END IF

  60 CONTINUE
    CLOSE(iIOUN2)
    kProfileUnitOpen=-1
  1080 FORMAT(A80)

    IF (iFoundGas < 0) THEN
      !! finally give up
      write(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
      CALL DOStop
    END IF

!*************************
! finally interp these onto the AIRS pressure levels
    DO iI = 1,iNumLevsx
      raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
      ! print *,iI,raRx110Press(iI),raRx110Temp(iI),raRx110MR(iI)
    END DO

    CALL r_sort_loglinear(raRx110Press,raRx110Temp,iNumLevsx,PLEVx110,raJunk,kMaxLayer+10)
    DO iI = 1,kMaxLayer+10
      raRx110Temp(iI) = raJunk(iI)
    END DO
    CALL r_sort_loglinear(raRx110Press,raRx110MR,iNumLevsx,PLEVx110,raJunk,kMaxLayer+10)
    DO iI = 1,kMaxLayer+10
      raRx110MR(iI) = raJunk(iI)
      raRx110Press(iI) = PLEVx110(iI)
    END DO
    iNumLevsx = kMaxLayer+10

    RETURN
    end SUBROUTINE ReadRefProf_Levels

!************************************************************************
! this basically loops and reads in ref profile for gas iG, if necessary changes pressures to N/m2
! input gas list (from user)                      = iaG
! output gaslist (preset list = 1 2 3 4 5 6 9 12) = union(user list,preset list)

! called by Tack_on_profile
    SUBROUTINE ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,iaPlanetMolecules,iPlanetMolecules, &
    iaG,iaGasUnits,iNumGases, &
    raR100Press,raR100Temp,raaR100MR,iRefLevels,laysORlevs)

    IMPLICIT NONE
          
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input/output
    INTEGER :: iaG(kMaxGas),iaGasUnits(kMaxGas),iNumGases    !! from user supplied list
    INTEGER :: laysORlevs                                    !! +1 if read 100 US Standard Ref layers, -1 if read 50 AFGL levels
    INTEGER :: iaPlanetMolecules(kMaxGas),iPlanetMolecules   !! important molecule list
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! output
! these are the individual reference profiles, at ~ (kMaxLayer + 10) layers
    REAL :: raR100Temp(kMaxLayer+10),raaR100MR(kMaxLayer+10,kMaxGas),raR100Press(kMaxLayer+10)
    INTEGER :: iRefLevels    !!! iRefLevels ~ 100 + 10 extra levels really high up : see SUBR ReadRefProf_Levels

! local
    REAL :: raRx110Temp(kProfLayer+10),raRx110MR(kProfLayer+10),raRx110Press(kProfLayer+10)

    INTEGER :: iL,iJ,iMid,iErr,iG,iFound,iNumGases0
    REAL :: rX,raR100Amt(kMaxLayer+10),raR100PartPress(kMaxLayer+10)
    CHARACTER(80) :: caRefFName
    CHARACTER(3) :: ca3
    INTEGER :: iNumLevsX,iIndex
    REAL :: raPPX(kMaxProfLayer),raQX(kMaxProfLayer) !!US Std layer ppress, amt
    REAL :: raPX(kMaxProfLayer), raTX(kMaxProfLayer) !!US Std layer press, temp
    REAL :: raMMX(kMaxProfLayer)
    CHARACTER(50) :: FMT
          
!! do a union of iaG and iaPlanetMolecules
    iNumGases0 = iNumGases
    DO iG = 1,iPlanetMolecules
      iFound = -1
      DO iL = 1,iNumGases
        IF (iaG(iL) == iaPlanetMolecules(iG)) iFound = +1
      END DO
      IF (iFound < 0) THEN
        write(kStdWarn,*)'gasID ',iaPlanetMolecules(iG),' not found in user list, adding .... '
        iNumGases = iNumGases + 1
        iaG(iNumGases) = iaPlanetMolecules(iG)
        iaGasUnits(iNumGases) = 12   !! VMR
      END IF
    END DO
    IF (iNumGases0 == iNumGases) THEN
      FMT = '(A,I2,A)'
      write(kStdWarn,*) 'you had all ',iPlanetMolecules,' important molecules in the input profile'
    ELSEIF (iNumGases0 < iNumGases) THEN
      FMT = '(A,I2,A,I2,A)'
      write(kStdWarn,FMT) 'you had ',iPlanetMolecules,' input gas profiles; now have ',iNumGases,' gas profiles'
    END IF
          
    DO iG = 1,iNumGases
      IF ((kPlanet /= 3) .OR. (laysORlevs == +1)) THEN
        !! this is old : use only (US) Standard Profile
        ! read kCARTA kProfLayer reference profile
        CALL FindReferenceName(caRefFName,iaG(iG),-1)
        CALL ReadRefProf(caRefFName,kMaxLayer,raR100Amt, &
          raR100Temp,raR100Press,raR100PartPress,iErr)
        iRefLevels = kMaxLayer
        ! change pressures from atm to N/m2
        DO iL = 1,kMaxLayer
          raR100Press(iL)     = raR100Press(iL) * kAtm2mb*100.0
          raR100PartPress(iL) = raR100PartPress(iL) * kAtm2mb*100.0
          raaR100MR(iL,iG)    = raR100PartPress(iL)/raR100Press(iL)
        END DO

      ELSEIF ((kPlanet == 3) .AND. (laysORlevs == -1)) THEN
        !ca3 = kcaLevsRefProf(1:3)
        ca3 = 'DNE'
        iIndex = index(kcaLevsRefProf,'DNE')
        IF (iIndex > 0) THEN
          write(kStdErr,*) 'oops : do not have a set of Levels Reference Profiles'
          CALL DoStop
        END IF
        !! this is new : can use one of the (AFGL) models
        CALL ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iaG(iG),iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
                    
        iRefLevels = iNumLevsx
        DO iL = 1,iNumLevsx
          raR100Press(iL)     = raRx110Press(iL)       !! already in N/m2
          raR100Temp(iL)      = raRx110Temp(iL)
          raaR100MR(iL,iG)    = raRx110MR(iL) /1.0e6   !! change from ppmv to MR
          raR100PartPress(iL) = raR100Press(iL) * raaR100MR(iL,iG)
          ! test IF (iL .LE. 100) print *,iaG(iG),iL,raaR100MR(iL,iG),raMMX(iL),raaR100MR(iL,iG)/raMMX(iL)
        END DO
        DO iL = iNumLevsx+1,kMaxLayer
          raR100Press(iL)     = 0.0
          raaR100MR(iL,iG)    = 0.0
          raR100PartPress(iL) = 0.0
        END DO
      END IF

      ! make sure pressures are decreasing with index ie layers going higher and higher
      IF (raR100Press(1) < raR100Press(2)) THEN
        !!! need to swap!!!
        iMid = ifloor(kMaxLayer/2.0)
        DO iL = 1,iMid
          iJ = kMaxLayer-iL+1
          rX = raR100Amt(iJ)
          raR100Amt(iJ) = raR100Amt(iL)
          raR100Amt(iL) = rX

          rX = raR100Temp(iJ)
          raR100Temp(iJ) = raR100Temp(iL)
          raR100Temp(iL) = rX

          rX = raR100Press(iJ)
          raR100Press(iJ) = raR100Press(iL)
          raR100Press(iL) = rX

          rX = raR100PartPress(iJ)
          raR100PartPress(iJ) = raR100PartPress(iL)
          raR100PartPress(iL) = rX

          rX = raaR100MR(iJ,iG)
          raaR100MR(iJ,iG) = raaR100MR(iL,iG)
          raaR100MR(iL,iG) = rX
        END DO
      END IF
    END DO
          
    RETURN
    end SUBROUTINE ReadRefProf_Units_laysORlevs

!************************************************************************
! this subroutine takes raaG_VMR and adds on new gas profiles using raaR100MR
! only to rPmin < rP < rPmax
    SUBROUTINE InterpNewGasProfiles_to_InputPlevels(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels, &
    iNumLevs,rPmin,rPmax,raP,raaG_VMR)

    IMPLICIT NONE
          
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iNumGases0,iaG0(kMaxGas),iNumGases,iaG(kMaxGas),iNumLevs,iRefLevels
    REAL :: rPmin,rPmax
    REAL :: raP(2*kProfLayer)
    REAL :: raaR100MR(kMaxLayer+10,kMaxGas),raR100Press(kMaxLayer+10)
! output
    REAL :: raaG_VMR(2*kProfLayer,kMaxGas)

! local
    INTEGER :: iG,iL,i1,i2,iJunk
    REAL :: raTempMR(kMaxLayer+10),raTempP(kMaxLayer+10),raTempMRX(2*kMaxLayer),rJunk
    CHARACTER(50) :: FMT

    INTEGER :: NTOTAL
    REAL :: V1S,V2S,DV
          
! first find pressures that span rPmin < rP < rPmax
    IF (raR100Press(iRefLevels) > rPmin*100.0) THEN
      write(kStdErr,*) 'iRefLevels = ',iRefLevels
      write(kStdErr,*) 'oops min pressure in ref database = ',raR100Press(iRefLevels),' N/m2'
      write(kStdErr,*) 'should be SMALLER than min pressure in supplied user profile = ',rPmin*100.0,' N/m2'
      CALL DoStop
    END IF

    IF (raR100Press(1) < rPmax*100.0) THEN
      write(kStdErr,*) 'iRefLevels = ',iRefLevels
      write(kStdErr,*) 'oops max pressure in ref database = ',raR100Press(1),' N/m2'
      write(kStdErr,*) 'should be GREATER than max pressure in supplied user profile = ',rPmax*100,' N/m2'
      CALL DoStop
    END IF

    write(kStdWarn,*) ' '

    i1 = iRefLevels
 10 CONTINUE
    IF ((raR100Press(i1)/100.0 < rPmin) .AND. (i1 > 1)) THEN
      i1 = i1 - 1
      GOTO 10
    END IF
    i1 = min(i1 + 1,iReflevels)

    i2 = 1
 20 CONTINUE
    IF ((raR100Press(i2)/100.0 > rPmax) .AND. (i2 < iRefLevels)) THEN
      i2 = i2 + 1
      GOTO 20
    END IF
    i2 = max(1,i2-1)
                
    FMT = '(A,F15.7,A,I3,A,F15.7)'
    rJunk = raR100Press(1)/100.0
    write(kStdWarn,FMT) 'Ref Database max press ',rJunk,' (at level ',1,'); user max press (in mb) = ',rPmax
    rJunk = raR100Press(iRefLevels)/100.0
    write(kStdWarn,FMT) 'Ref Database min press ',rJunk,' (at level ',iRefLevels,'); user min press (in mb) = ',rPmin
    FMT = '(A,I3,I3,A)'
    write(kStdWarn,FMT) 'Ref profiles between levels ',i1,i2,' spans rPmin,rPmax'
          
    DO iL = i2,i1
      raTempP(iL-i2+1) = raR100Press(iL)/100.0
    END DO

    DO iG = iNumGases0 + 1,iNumGases
      DO iL = i2,i1
        raTempMR(iL-i2+1) = raaR100MR(iL,iG)
      END DO
      Call r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
      DO iL = 1,iNumLevs
        raaG_VMR(iL,iG) = raTempMRX(iL)
      END DO
    END DO

    RETURN
    end SUBROUTINE InterpNewGasProfiles_to_InputPlevels

!************************************************************************
! for testing just do SUBROUTINE InputMR_profile(caPfName)    !! for testing
    SUBROUTINE InputMR_profile(caPfName,iNumLays,iNumGases,iaG,iLowestLev, &
    raTout,raAmountOut,raZout, &
    raPout,raaQout,raaPartPressout, &
    raPbndFinal,raTbndFinal,iZbndFinal)

    IMPLICIT NONE

    INTEGER :: iplev
    include '../INCLUDE/TempF90/kcartaparam.f90'
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'
    include '../INCLUDE/TempF90/airslevelheightsparam.f90'

! input
    CHARACTER(80) :: caPfName
! output
    REAL :: raTout(kProfLayer)                  !! in K
    REAL :: raAmountOut(kProfLayer)             !! in molecules/cm2
    REAL :: raZout(kProfLayer+1)                !! in meters, notice this is at LEVELS
    REAL :: raPout(kProfLayer)                  !! in N/m2 (even though input TAPE5,TAPE6, kRTP=-10, pressures are in mb)
    REAL :: raaQout(kProfLayer,kMaxGas)         !! in molecules/cm2
    REAL :: raaPartPressout(kProfLayer,kMaxGas) !! in N/m2
    INTEGER :: iNumLays,iNumGases,iaG(kMaxGas)
    INTEGER :: iLowestLev
    REAL :: raPbndFinal(kProfLayer+1)           !! final output pressure boundaries, default = 101 AIRS levels
    REAL :: raTbndFinal(kProfLayer+1)           !! final output temperatures at boundaries
    INTEGER :: iZbndFinal                       !! number of output plev bdries,
!! should be iNumLays+1 if using LBLRTM bdies, else 101
          
! local
    INTEGER :: iG,iL,iJ,iK,iaGasUnits(kMaxGas)
    INTEGER :: iFound,iHighestLay,iNumLevs
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
    REAL :: raX(kMaxGas),rX,rP,rT,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
    REAL :: rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
    INTEGER :: iKCARTADirectionUPorDOWN,iMid,iCnt,iOffSet
    REAL :: gamma,z,dz,slope,amount,junk
    REAL :: rFracBot
    REAL :: raPBnd(2*kProfLayer)  !! do we want user defined pressure level boundaries??
    INTEGER :: iZbnd              !! do we want user defined pressure level boundaries??

    REAL :: raPX(kProfLayer+1),raZX(kProfLayer+1),raTX(kProfLayer+1)
    REAL :: raaG_VMRX(kProfLayer+1,kMaxGas),raLayDensityX(kProfLayer+1)
    REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)
          
    CALL Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS, &
      rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta)

    DO iL = 1,kMaxLayer+1
      raPbndFinal(iL) = PLEV_KCARTADATABASE_AIRS(iL)
      raTbndFinal(iL) = 0.0
    END DO
    iZbndFinal = kMaxLayer+1

    iZbnd = +1   !! default to using the input levels as boundaries when doing levels --> layers integration
    iZbnd = -1   !! default to using AIRS 101 levels  as boundaries when doing levels --> layers integration
          
! >>> read in user supplied prof
! this will be at whatever gas units eg RH g/kg  VMR etc, but then CONVERTED to VMR (gasunit 12)
    IF (kRTP == -10) THEN
      CALL ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
        iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases, &
        raP,raT,iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)
      iZbnd = -1
    ELSEIF (kRTP == -5) THEN
      CALL ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
        iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases, &
        raP,raT,raAlt,iZbnd,raPBnd, &
        iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)
    ELSEIF (kRTP == -6) THEN
      !! edited TAPE6, which contains first few parts of TAPE5, and then just profile info from TAPE6 (mol/cm2)
      CALL ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
        iNumLays,rPSurf,rTSurf,rHSurf,iNumGases, &
        raPX,raTX,raLayDensityX,raZX, &
        raPressLevelsX,raTPressLevelsX,raAltitudesX, &
        iaG,iaGasUnits,raaG_VMRX,rPMin,rPMax,rYear,rLat,rLon)
      iZbnd = -1
    ELSE
      write(kStdErr,*) 'huh?? reading text levels file or text LBLRTM TAPE5/TAPE6 MODIFIED',kRTP
      Call DoStop
    END IF

    IF (kRTP /= -6) THEN
      ! >>> interp onto the kCARTA database profiles
      ! >>> also sort the gasIDs into increasing order

      ! find pseudo temperature at this level
      ! recall pav = (p2-p1)/log(p2/p1)
      !        tav = (T2-T1)/log(p2/p1)
              
      !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
      !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
      !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
              
      IF (iZbnd < 0) THEN
        ! >>>>>>>>>> default AIRS 101 levels

        ! >>>> tack on Standard Profile to User Profile if needed;
        ! >>>> also add on a few gases (so we have gasIDs [1 2 3 4 5 6 9 12]) if needed
        ! old and new gases are ALWAYS going be in units = 12 (volume mixing ratio)
        ! have to be careful when adding in upper level info to old gases, as this could be a mix of units!!!!!!!!!!
        ! taper the tacked on profile so that within about 4 sublevels, you are just tacking on US Std (or default ref prof)
        ! can do the CO2 = 370 +yy-2002), CH4 = 1.7 + (yy-2002)*0.01
        CALL Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR, &
            rPMin,rPMax,rYear)
                      
        ! now check everything already is in VMR (gas unit 12); change to that if needed
        ! also adjust MR (wet water) if this is planet Earth
        CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)

        CALL InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta, &
          PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev)

        rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/(PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'iX = iLowestLev ==>'
        write(kStdWarn,*)' PLEV_KCARTADATABASE_AIRS(iX),rPSurf,PLEV_KCARTADATABASE_AIRS(iX+1) = '
        write(kStdWarn,*) PLEV_KCARTADATABASE_AIRS(iLowestLev),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev+1)
        write(kStdWarn,*) 'rfracBot = ',rFracBot
              
        ! >>> do the integrals using integrals wrt p!!! closer to klayers
        CALL DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon, &
            PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot, &
            raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
            iNumLays = kProflayer-iLowestLev+1  !! this is now LAYER number, max == kProfLayer when iLowestLev = 1

        DO iL = 1,101
          raPbndFinal(iL) = raPbndFinal(iL)/100.0    !! convert N/m2 to mb
        END DO

      ELSEIF (iZbnd >= 0) THEN
        ! >>>>>>>>>> user boundaries
        DO iL = 1,iZbnd
          raPbndFinal(iL) = raPbnd(iL)/100.0    !! convert N/m2 to mb
        END DO
        DO iL = iZbnd+1,kProfLayer+1
          raPbndFinal(iL) = 0.005
          raPbndFinal(iL) = raPbnd(iZbnd)/100.0   !! convert N/m2 to mb
        END DO
        iZbndFinal = iZbnd

        ! now check everything already is in VMR (gas unit 12); change to that if needed
        ! also adjust MR (wet water) if this is planet Earth
        CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)

        CALL InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta, &
          PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev,iZbnd,raPBnd)

        rFracBot = (rPSurf-raPBnd(iLowestLev+1))/(raPBnd(iLowestLev)-raPBnd(iLowestLev+1))
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'iX = iLowestLev ==>'
        write(kStdWarn,*)' raPBnd(iX),rPSurf,raPBnd(iX+1) = '
        write(kStdWarn,*) raPBnd(iLowestLev),rPSurf,raPbnd(iLowestLev+1)
        write(kStdWarn,*) 'rfracBot = ',rFracBot

        ! >>> do the integrals using integrals wrt p!!! closer to klayers
        CALL DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon, &
            raPBnd,raAlt,iZbnd,rfracBot, &
            raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
            iNumLevs = iZbnd-iLowestLev+1
            iNumLays = iNumLevs - 1       !! this is now LAYER number, max == iZbnd-1 when iLowestLev = 1

      END IF
                   
    ELSEIF (kRTP == -6) THEN
      CALL DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon, &
        PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot, &
        raPX,raTX,raLayDensityX,raaG_VMRX, &
        raPressLevelsX,raTPressLevelsX,raAltitudesX, &
        raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)
      rfracBot = 1.0
      iNumLevs = iNumLays+1
      iZbnd = iNumLevs
      iOffSet = (kProfLayer+1)-(iNumLays+1)
      DO iL = 1,iZbnd
        raPbndFinal(iL)    = raPressLevelsX(iL)         !! in mb
        raTbndFinal(iL)    = raTPressLevelsX(iL)        !! in K
        raZout(iL+iOffSet) = raAltitudesX(iL)*1000.0    !! reset raZout ... and use  m
      END DO
      iZbndFinal = iNumLevs
    END IF

    RETURN
    end SUBROUTINE InputMR_profile

!************************************************************************
! >>> read in user supplied prof
! this will be at whatever gas units eg RH g/kg  VMR etc, but then CONVERTED to VMR (gasunit 12)
    SUBROUTINE  ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
    iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,raP,raT, &
    iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)

    IMPLICIT NONE
         
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    CHARACTER(80) :: caPfName
    REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
! output
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon,tAugment
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas)

! local var
    INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,iReadInOK,iIgnore,iYes
    REAL :: raX(kMaxGas),rX,rP,rT,raP_Nm2(2*kProfLayer),rPRefMin,rPRefMax,rPRefMinAugmented
    CHARACTER(80) :: caStr
    CHARACTER(40) :: caaUnit(50)

    DO iL = 1,50
      caaUnit(iL) = '** Not Known **'
    END DO
    caaUnit(01) = 'layer amt molecules/cm2'
    caaUnit(10) = 'vmr in ppmv'
    caaUnit(11) = 'vmr in ppbv'
    caaUnit(20) = 'mmr in g/kg'
    caaUnit(21) = 'mmr in g/g'
    caaUnit(30) = 'part press in mb (not accepted)'
    caaUnit(31) = 'part press in atm (not accepted)'
    caaUnit(40) = 'RH % (not accepted)'
    caaUnit(41) = 'RH fraction (not accepted)'

    rPRefMax = 1100.0 !! 1100  mb is AIRS RefPress LowerBdry      see SUBR ReadRefProf_Levels
    rPRefMin = 3.0e-3 !! 0.005 mb extended to 0.003 mb = 0.3 N/m2 see SUBR ReadRefProf_Levels
          
    rPmin = +1.0e6
    rPmax = -1.0e+6

    write(kStdWarn,*) 'Reading in user supplied TXT LVLS file .....'

 1070 FORMAT('ERROR! number ',I5,' opening PRFILE path LEVELS profile file',/,A80)
    iIgnore = -1
    iIOUN2 = kProfileUnit
    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO)
    IF (iErrIO /= 0) THEN
      iErr = 1
      WRITE(kStdErr,1070) iErrIO, caPfname
      CALL DoSTOP
    ENDIF
    kProfileUnitOpen = 1

    READ (iIOUN2,5030,ERR=13,END=13) caStr

    READ (iIOUN2,*) iNumLevs
    IF (iNumLevs > 2*kProfLayer) THEN
      write(kStdErr,*) 'iNumLevs > 2*kProfLayer',iNumLevs,kProfLayer
      CALL DoStop
    END IF

    READ (iIOUN2,*) rPSurf,rTSurf,rHSurf
    IF ((rHSurf > rHmaxKCarta) .OR. (rHSurf < rHminKCarta)) THEN
      write(kStdErr,*) 'need rHmaxKCarta >= rHSurf >= rHminKCarta but have'
      write(kStdErr,*) '(rHmaxKCarta,rHSurf,rHminKCarta) = ',rHmaxKCarta,rHSurf,rHminKCarta
      CALL DoStop
    END IF
    IF ((rPSurf > rPmaxKCarta) .OR. (rPSurf < rPminKCarta)) THEN
      write(kStdErr,*) 'need rPmaxKCarta >= rPSurf >= rPminKCarta but have'
      write(kStdErr,*) '(rPmaxKCarta,rPSurf,rPminKCarta) = ',rPmaxKCarta,rPSurf,rPminKCarta
      CALL DoStop
    END IF
    IF ((rTSurf > kStempMax) .OR. (rTSurf < kStempMin)) THEN
      write(kStdErr,*) 'need kStempMax >= rTSurf >= kStempMin but have'
      write(kStdErr,*) '(kStempMax,rTSurf,kStempMin) = ',kStempMax,rTSurf,kStempMin
      CALL DoStop
    END IF

    raRTP_TxtInput(1) = rPSurf
    raRTP_TxtInput(2) = rTSurf
    raRTP_TxtInput(3) = rHSurf
          
    READ (iIOUN2,*) rYear,rLat,rLon

    READ (iIOUN2,*) iNumGases
    IF (iNumGases > kMaxGas) THEN
      write(kStdErr,*) 'iNumGases > kMaxGas',iNumGases,kMaxGas
      CALL DoStop
    END IF
           
    READ (iIOUN2,*) (iaG(iJ),iJ=1,iNumGases)
    READ (iIOUN2,*) (iaGasUnits(iJ),iJ=1,iNumGases)

    write(kStdWarn,'(A)') '-------------------------------------------------'
    write(kStdWarn,'(A)') 'Index  GasID  UnitsCode  GasUnit'
    write(kStdWarn,'(A)') '-------------------------------------------------'
    DO iG = 1,iNumGases
      write(kStdWarn,'(3(I5),A,A)') iG,iaG(iG),iaGasUnits(iG),'      ',caaUnit(iaGasUnits(iG))
    END DO
    write(kStdWarn,'(A)') '-------------------------------------------------'

!! this works fine if levels entered so thatpres(1) > pres(2) > pres(3) ie from ground up
!! may be buggy for reverse (from TOA down)
    write(kStdWarn,*) 'Input levels profile : ',iNumLevs,' levels for ',iNumGases,' gases'
    write(kStdWarn,*) 'PSurf = ',rPSurf,' mb;  TSurf = ',rTSurf,' K '
    iReadInOK = 0
    iIgnore = -1
    DO iL = 1,iNumLevs
      READ (iIOUN2,*) rP,rT,(raX(iJ),iJ=1,iNumGases)
      IF ((rP <= rPRefMax) .AND. (iIgnore < 0)) THEN
        IF ((rP < rPRefMin) .AND. (rP < rPmin)) iIgnore = +1
        iReadInOK = iReadInOK + 1
        raP_Nm2(iReadInOK) = rP * 100.0  !! change from mb to N/m2
        raP(iReadInOK) = rP
        raT(iReadInOK) = rT
        IF (rPmax <= raP(iReadInOK)) rPmax = raP(iReadInOK)
        IF (rPmin > raP(iReadInOK)) rPmin = raP(iReadInOK)
        DO iJ = 1,iNumGases
          raaG_VMR(iReadInOK,iJ) = raX(iJ)
        END DO
        write(kStdWarn,5040) iReadInOK,raP(iReadInOK),raT(iReadInOK),raaG_VMR(iReadInOK,1)
      ELSEIF (rP <= rPRefMin) THEN
        !! things are OK, we are ignoring pressures for levels above 0.005 mb
        write(kStdWarn,5041) iL,rP,rT
      ELSEIF ((rP > rPRefMin) .AND. (iIgnore > 0)) THEN
        !! things are BAD, this pressure is for levels below 0.005 mb
        !! but we already found a pressre for above 0.005 mb, so pinput is UNSORTED??
        write(kStdWarn,5042) iL,rP,rT
        CALL DoStop
      END IF
    END DO
          
 13 CONTINUE
    CLOSE(iIOUN2)
    kProfileUnitOpen = -11
          
    IF (iNumLevs /= iReadInOK) THEN
      WRITE(kStdWarn,*) 'Input textfile had ',iNumLevs,' levels of which ',iReadInOK,' lay between ',rPRefMax,rPRefMin,' mb'
    ELSE
      WRITE(kStdWarn,*) 'All ',iNumLevs,' input levels in textfile lay between ',rPRefMax,rPRefMin,' mb'
    END IF
    iNumLevs = iReadInOK
          
 5030 FORMAT(A80)
 5040 FORMAT('iReadInOK, rP (mb), rT (K), raG1 = ',I3,' ',F14.6,' ',F8.3,' ',ES12.5)
 5041 FORMAT('Ignore (too low P) rP (mb), rT (K), raG1 = ',I3,' ',ES12.5,' ',F8.3)
 5042 FORMAT('Bad (press below 0.005 mb (UNSORTED PinLEVS??)) rP (mb), rT (K), raG1 = ',I3,' ',ES12.5,' ',F8.3)
     
    write(kStdWarn,*)'   KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
    write(kStdWarn,*)'   kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
    write(kStdWarn,*)' input file : spres (mb)/Height(km)     = ',rPSurf,rHSurf

    rPSurf = rPSurf * 100.0
    rHSurf = rHSurf * 1000.0

! make sure pressures are decreasing with index ie layers going higher and higher
    IF (raP(1) < raP(2)) THEN
      !!! need to swap!!!
      iMid = ifloor(iNumLevs/2.0)
      DO iL = 1,iMid
        iJ = iNumLevs-iL+1
        rX = raP(iJ)
        raP(iJ) = raP(iL)
        raP(iL) = rX

        rX = raP_Nm2(iJ)
        raP_Nm2(iJ) = raP_Nm2(iL)
        raP_Nm2(iL) = rX

        rX = raT(iJ)
        raT(iJ) = raT(iL)
        raT(iL) = rX

        DO iG = 1,iNumGases
          raX(iG) = raaG_VMR(iJ,iG)
          raaG_VMR(iJ,iG) = raaG_VMR(iL,iG)
          raaG_VMR(iL,iG) = raX(iG)
        END DO
      END DO
    END IF
          
    write(kStdWarn,*) 'input file : Highest altitude (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
    write(kStdWarn,*) ' '

    iYes = +1    !!! reset topmost pressure to rPRefMin = 0.005
    iYes = -1    !!! leave topmost pressure level alone, do not reset to rPRefMin = 0.005
!!! actually we augment the refdatabase to 0.2750000 N/m2 = 0.00275
    rPRefMinAugmented = 0.00275
    IF ((rPmin < rPRefMin) .AND. (iYes > 0)) THEN
      !!! interpolate gas concentrations last point to rPRefMin
      write(kStdWarn,*) 'highest press point is ',rPmin,' mb higher than ref database point ',rPRefMin,'mb fixing!'
      write(kStdWarn,*) ' '
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raX(iL) = raaG_VMR(iL,iG)
        END DO
        CALL r_sort_loglinear(raP,raX,iNumLevs,rPRefMin,rX)
        raaG_VMR(iNumLevs,iG) = rX
      END DO

      !!! interpolate T(z) last point to rPRefMin
      CALL r_sort_loglinear(raP,raT,iNumLevs,rPRefMin,rX)
      raT(iNumLevs) = rX
              
      !!! set last point to 0.005 mb
      rPmin = rPRefMin
      raP(iNumLevs) = rPRefMin
      raP_Nm2(iNumLevs) = rPRefMin * 100.0

    ELSEIF (rPmin < rPRefMinAugmented) THEN
      !!! force interpolate gas concentrations last point to rPRefMinAugmented
      write(kStdWarn,*) 'highest press point is ',rPmin,' mb higher than augmented ref dbase pt ',rPRefMinAugmented,'mb fixing!'
      write(kStdWarn,*) ' '
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raX(iL) = raaG_VMR(iL,iG)
        END DO
        CALL r_sort_loglinear(raP,raX,iNumLevs,rPRefMinAugmented,rX)
        raaG_VMR(iNumLevs,iG) = rX
      END DO

      !!! interpolate T(z) last point to rPRefMin
      CALL r_sort_loglinear(raP,raT,iNumLevs,rPRefMinAugmented,rX)
      raT(iNumLevs) = rX
              
      !!! set last point to 0.00275 mb
      rPmin = rPRefMinAugmented
      raP(iNumLevs) = rPRefMinAugmented
      raP_Nm2(iNumLevs) = rPRefMin * 100.0
    END IF
            
          
! now change all units to PPMV (units 10) and thn to VMR (units 12)
    DO iG = 1,iNumGases
    CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP_Nm2,raT,raaG_VMR,-1)
      DO iL = 1,iNumLevs
        raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0e6
        iaGasUnits(iG) = 12   !! VMR
      END DO
    END DO

    RETURN
    END SUBROUTINE 

!************************************************************************
    SUBROUTINE ReadRefProf_Levels2(iGasID,raP,iNumLevsIN,raX)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID,iNumLevsIN
    REAL :: raP(2*kProfLayer)
! output
    REAL :: raX(2*kProfLayer)

! local
    REAL :: raRx110Temp(2*kProfLayer),raRx110MR(2*kProfLayer),raRx110Press(2*kProfLayer)
    REAL :: rLat,raJunk(kProfLayer*2)
    CHARACTER(80) :: caPFname,caStr,caComment
    INTEGER :: iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec,iNumLevsx

    caPFname = kcaLevsRefProf

    IF ((kAFGLProf < 1) .OR. (kAFGLProf > 6)) THEN
      write(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
      CALL DoStop
    END IF

    iErr = 0
    iIOUN2 = kProfileUnit
    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED', IOSTAT=iErrIO)
    IF (iErrIO /= 0) THEN
      iErr=1
      WRITE(kStdErr,1070) iErrIO, caPfname
      CALL DoSTOP
    ENDIF
 1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
    kProfileUnitOpen=1

    iAFGL = 0
    iFoundGas = -1
          
 10 CONTINUE
    READ(iIOUN2,1080) caStr
    IF (caStr(1:1) == '!') THEN
      GOTO 10   !! keep reading comments
    ELSE
      READ (caStr,*) iNumLevsx
    END IF
    IF (iNumLevsx > 2*kProfLayer) THEN
      write(kStdErr,*) 'oops iNumLevsx > 2*kProfLayer ',iNumLevsx,kProfLayer*2
      CALL DoStop
    END IF

 20 CONTINUE
    iAFGL = iAFGL + 1

 30 CONTINUE
    READ(iIOUN2,1080) caStr
    IF (caStr(1:1) == '!') THEN
      GOTO 30   !! keep reading comments
    ELSE
      caComment = caStr
    END IF

 40 CONTINUE
    READ(iIOUN2,1080) caStr
    IF (caStr(1:1) == '!') THEN
      GOTO 40   !! keep reading comments
    ELSE
      READ (caStr,*) rLat
    END IF

    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! Altitude(km)
    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

    IF (iAFGL == kAFGLProf) THEN
      !write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Press(mb)
      READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Temp(K)
      READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)

      !write(kStdWarn,*) '   need to find profile for gasID ',iGasID
    ELSE
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Press(mb)
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Temp(K)
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
    END IF

    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! density (cm-3)
    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

    DO iG = 1,7
      IF (iG == iGasID) THEN
        !write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
        iFoundGas = +1
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
      ELSE
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
      END IF
    END DO

    IF ((iAFGL == kAFGLProf) .AND. (iFoundGas > 0)) THEN
      GOTO 60   !! found the AFGL prof and found the gas; done!
    ELSEIF ((iAFGL == 6) .AND. (iFoundGas < 0)) THEN
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! modend
      GOTO 50   !! found the AFGL prof but not found the gas
    ELSEIF ((iAFGL < 6) .OR. (iFoundGas < 0)) THEN
      !! either did not find the gas or the AFGL prof
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! modend
      GOTO 20
    END IF

 50 CONTINUE
    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! constituent profs
    READ(iIOUN2,1080) caStr   !! comment
    READ(iIOUN2,1080) caStr   !! mingas

    iG = 7
    55 CONTINUE
    iG = iG + 1
    IF (iG == iGasID) THEN
      ! write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
      iFoundGas = +1
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! gas name
      READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
    ELSE
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! gas name
      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
    END IF
    IF ((iG < 28) .AND. (iFoundGas < 0)) THEN
      GOTO 55
    END IF
          
    IF (iFoundGas < 0) THEN
      !! bweh try xsec gases
      123 CONTINUE
      READ(iIOUN2,1080) caStr   !! comment
      IF (caStr(1:1) == '!') THEN
        GOTO 123
      ELSEIF (caStr(1:6) == 'DATEND') THEN
        GOTO 60       !!! oh oh end of file, give up!!!!!
      ELSE
        READ (caStr,*) iXsec
        IF (iXsec == iGasID) THEN
          iFoundGas = 1
          READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
        ELSE
          READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
          GOTO 123
        END IF
      END IF
    END IF

    60 CONTINUE
    CLOSE(iIOUN2)
    kProfileUnitOpen=-1
 1080 FORMAT(A80)

    IF (iFoundGas < 0) THEN
      !! finally give up
      write(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
      CALL DOStop
    END IF

!*************************
! finally interp these onto the raP pressure levels
    DO iI = 1,iNumLevsx
      raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
    END DO

    CALL r_sort_loglinear(raRx110Press,raRx110MR,iNumLevsx,raP,raJunk,iNumLevsIN)
    DO iI = 1,iNumLevsIN
      raX(iI) = raJunk(iI)
    END DO

    RETURN
    end SUBROUTINE ReadRefProf_Levels2
!************************************************************************
    SUBROUTINE Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS, &
    rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! output
    REAL :: rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
    INTEGER :: iKCARTADirectionUPorDOWN
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1),DATABASELEVHEIGHTS(kMaxLayer+1)

! local
    INTEGER :: iL

    rPminKCarta = +1e10
    rPmaxKCarta = -1e10
    rHminKCarta = +1e10
    rHmaxKCarta = -1e10
    rPmin = +1e10
    rPmax = -1e10
    iKCARTADirectionUPorDOWN = +1   !! assume increading index == lower pressure ==> going UP
    DO iL = 1,kMaxLayer+1
      PLEV_KCARTADATABASE_AIRS(iL) = PLEV_KCARTADATABASE_AIRS(iL) * 100.0  !! change mb --> N/m2
      DATABASELEVHEIGHTS(iL)       = DATABASELEVHEIGHTS(iL) * 1000.0       !! change km --> m
      IF (PLEV_KCARTADATABASE_AIRS(iL) > rPmaxKCarta) rPmaxKCarta = PLEV_KCARTADATABASE_AIRS(iL)
      IF (PLEV_KCARTADATABASE_AIRS(iL) < rPminKCarta) rPminKCarta = PLEV_KCARTADATABASE_AIRS(iL)
      IF (DATABASELEVHEIGHTS(iL) > rHmaxKCarta) rHmaxKCarta = DATABASELEVHEIGHTS(iL)
      IF (DATABASELEVHEIGHTS(iL) < rHminKCarta) rHminKCarta = DATABASELEVHEIGHTS(iL)
    END DO
    IF ((PLEV_KCARTADATABASE_AIRS(1) > PLEV_KCARTADATABASE_AIRS(2)) .AND. &
      (DATABASELEVHEIGHTS(1) < DATABASELEVHEIGHTS(2))) THEN
      iKCARTADirectionUPorDOWN = +1   !! increasing index == lower pressure ==> going UP
    ELSEIF ((PLEV_KCARTADATABASE_AIRS(1) < PLEV_KCARTADATABASE_AIRS(2)) .AND. &
      (DATABASELEVHEIGHTS(1) > DATABASELEVHEIGHTS(2))) THEN
      iKCARTADirectionUPorDOWN = -1   !! decreasing index == higher pressure ==> going DOWN
      write (kStdErr,*) 'Bummer : being lazy, wanted this to be going UP'
      CALL DoStop
    ELSE
      write (kStdErr,*) 'Inconsistency : database pressures in one dir, heights in another'
      CALL DoStop
    END IF

    RETURN
    end SUBROUTINE Init_n_layers

!************************************************************************
! http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
! this one takes the ppmv etc and applies to TOP bdry
!  (ie record 2.1.1 if TOA is level 1 ... GND at level L, this CORRECTLY applies ppmv to (L)
! ALTZ(L-1),  PZ(L-1),  TZ(L-1),  ATLZ(L),  PZ(L),  TZ(L)
!    ALTZ(L-1), PZ(L-1) and TZ(L-1) are only required for the first layer.
!    LBLRTM assumes that these quantites are equal to the top of the previous
!    layer for L > 1.
                                                          
    SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
    iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases, &
    raP,raT,raAlt,iZbnd,raPBnd, &
    iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)

    IMPLICIT NONE
         
    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iPLEV
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'
          
! input
    CHARACTER(80) :: caPfName
    REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
! output
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer) !! levels info
    REAL :: raPavg(2*kProfLayer),raTavg(2*KProfLayer)  !!! TAPE5 can have layers average p and T
    REAL :: raPBnd(2*kProfLayer)  !! do we want user defined pressure level boundaries??
    INTEGER :: iZbnd              !! do we want user defined pressure level boundaries??

! local var
    INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
    REAL :: rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,raSumCheck(kMaxGas)
    CHARACTER(80) :: caStr,caStrX,caStrY
    CHARACTER(30) :: caStr30
    CHARACTER(1) ::  c1
    CHARACTER(2) ::  c2
    CHARACTER(3) ::  c3
    CHARACTER(4) ::  c4
    CHARACTER(5) ::  c5
    INTEGER :: iWriteRTP,iNumGasesBAD,iaBadGasProfile(kMaxGas),iReplaceZeroProf
    INTEGER :: IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS
    INTEGER :: XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL
    REAL :: rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,rDVOUT
    REAL :: raEmiss(3),raRefl(3),rSecnto
    INTEGER :: ILNFLG,iForm,iBmax,iOffSet
    REAL :: raX(KMaxGas),rX,rSatZen
    REAL :: raZbnd(2*kProfLayer)
    INTEGER :: iDefault,iAIRS101_or_LBL_levels

    iDefault = +1   !! use AIRS101 levels for the integration
    iAIRS101_or_LBL_levels = +1 !! use AIRS101 levels for the integration
    iAIRS101_or_LBL_levels = -1 !! use LBLRTM  levels for the integration
    iAIRS101_or_LBL_levels = iaaOverrideDefault(3,1)
    IF (abs(iAIRS101_or_LBL_levels) > 1) THEN
     write(kStdErr,*) 'invalid iAIRS101_or_LBL_levels ',iAIRS101_or_LBL_levels
     CALL DoStop
    END IF
    IF (iDefault /= iAIRS101_or_LBL_levels) THEN
        write(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 when doing integration from levels to layers'
        write(kStdErr,*) 'iDefault, iAIRS101_or_LBL_levels = ',iDefault,iAIRS101_or_LBL_levels
    END IF

    DO iL = 1,2*kProfLayer
      raPavg(iL) = -9999.0
      raTavg(iL) = -9999.0
    END DO
          
    rPmin = +1.0e6
    rPmax = -1.0e+6
    iNXsec = -1

    IF (kProfLayer > kMaxLayer) THEN
      write(kStdErr,*) 'oops, want to save LBLRTM temperatures into kLBLRTM_levelT'
      write(kStdErr,*) 'but kProfLayer > kMaLayer+1'
      CALL DoStop
    ELSE
      DO iG = 1,kMaxLayer
        kLBLRTM_levelT(iG) = 0.0
        kLBLRTM_layerTavg(iG) = 0.0
      END DO
      kLBLRTM_levelT(kMaxLayer+1) = 0.0
    END IF
          
    iZbnd = -1   !!! assume we want to use default 101 AIRS levels
    DO iL = 1,kProfLayer+1
      raPbnd(iL) = PLEV_KCARTADATABASE_AIRS(iL)
    END DO
          
    write(kSTdWarn,*) 'Reading in LBLRTM TAPE5 .....'

 1070 FORMAT('ERROR! number ',I5,' opening PRFILE path LBLRTM profile file' ,/,A80)
    iIOUN2 = kProfileUnit
    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO)
    IF (iErrIO /= 0) THEN
      iErr=1
      WRITE(kStdErr,1070) iErrIO, caPfname
      CALL DoSTOP
    ENDIF
    kProfileUnitOpen=1

!! read till you get the $ character, RECORD 1.1
 10 CONTINUE
    READ (iIOUN2,111,ERR=13,END=13) caStr
    IF (caStr(1:1) /= '$') GOTO 10

    READ (iIOUN2,111,ERR=13,END=13) caStr
    CALL read_record_1p2(caStr, &
    IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)
    IF (ICNTNM == 6) THEN
      READ (iIOUN2,111,ERR=13,END=13) caStr
      CALL read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)
    END IF
    IF (IEMIT == 2) THEN
      CALL DOSTOPMesg('oops cannot handle solar file$')
    END IF

    IF ((IHIRAC > 0) .OR. (IAERSL > 0) .OR. (IEMIT == 1) .OR. (IATM == 1) .OR. (ILAS > 0)) THEN
      READ (iIOUN2,111,ERR=13,END=13) caStr
      CALL read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)
    END IF

!! iotflg comes from solar data
!! IF ((IEMIT .EQ. 1) .OR. (IEMIT .EQ. 2 .AND.  IOTFLG .EQ. 1)) THEN
    IF ((IEMIT == 1) .OR. (IEMIT == 2)) THEN
      READ (iIOUN2,111,ERR=13,END=13) caStr
      CALL read_record_1p4(caStr,rTSurf,raEmiss,raRefl)
    END IF

    111 FORMAT(A80)
    write(kStdWarn,*) 'TAPE 5 has iAtm = ',iAtm
    IF (IATM == 0) THEN
      !! need to read in profile
      CALL read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)
      CALL read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raPavg,raTavg,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck, &
        rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta, &
        iaG,iaGasUnits,iForm,iNumGases)
      IF (iAIRS101_or_LBL_levels == +1) THEN
        iZbnd = -1            !! use default AIRS 101 levels, they should span past the rPSurf
      ELSE
        iZbnd = iNumLevs  !! use the LBLRTM pressure levels
      END IF
      IF (iZbnd > 0) THEN
        raPbnd(1) = rPSurf
        iOffSet = kProfLayer-(iZbnd-1)
        DO iG = 1,iNumLevs
          IF (iG < iNumLevs) kLBLRTM_layerTavg(iG+iOffset) = raTavg(iG)
          kLBLRTM_levelT(iG+iOffSet) = raT(iG)
          raPbnd(iG) = raP(iG) * 100.0 !! change from mb to N/m2, as Psurf is finally in N/m2
          ! raPbnd(iG) = raP(iG)         !! keep in mb
        END DO
      END IF
      IF (iXsect == 1) THEN
        CALL read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)
      END IF
    ELSE
      rViewAngle = -9999.0
      rTopHgt    = 705000.0
      CALL read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)
      IF (iBmax > 0) THEN
        iZbnd = +iBmax
        CALL read_record_3p3b(iIOUN2,iBmax,raZbnd)
      END IF
      CALL read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf, &
        rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
        raAlt,raP,raT,raSumCheck, &
        iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon, &
        iZbnd,raZbnd,raPbnd)
    END IF
!!! TAPE5 can also include scanangle info and output file name info for flux, see Record 6,6.1
!!! Unless you want to read in rCos(flux angles) for the three angles, as yu can set them in nm_radnce, using
!!!   iAtmLoop = 3, raAtmLoop(1,2,3)

    write(kStdWarn,*) 'After exiting the equivalent of nolgas and xsecgas in TAPE5 LBLRTM, found',iNumGases,' gases'
    write(kStdWarn,*) '    these are the gasIDs'
    write(kStdWarn,*) (iaG(iJ),iJ=1,iNumGases)
    write(kStdWarn,*) '    these are the gas units'
    write(kStdWarn,*) (iaGasUnits(iJ),iJ=1,iNumGases)
    write(kStdWarn,*) ' '
          
!! now see if data was acutally read in
    iNumGasesBAD = 0
    DO iJ = 1,iNumGases
      IF (raSumCheck(iJ) < 1.0e-20) THEN
        write(kStdWarn,*) ' read in TAPE5, following gas had    zero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
        iNumGasesBAD = iNumGasesBAD + 1
        iaBadGasProfile(iNumGasesBAD) = iJ
      ELSE
        write(kStdWarn,*) ' readin TAPE5, following gas had nonzero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
      END IF
    END DO

!! need to fix the BAD gases
    iDefault = +1
    iReplaceZeroProf = -1    !! assume user knows why there is a ZERO everywhere gas profile
    iReplaceZeroProf = +1    !! assume user wants to replace ZERO everywhere gas profile with climatology
    iReplaceZeroProf = iaaOverrideDefault(3,2)
    IF (abs(iReplaceZeroProf) > 1) THEN
      write(kStdErr,*) 'invalid iReplaceZeroProf ',iReplaceZeroProf
      CALL DoStop
    END IF
    IF ((iDefault /= iReplaceZeroProf) .AND. (iNumGasesBAD > 0)) THEN
      write(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
      write(kStdErr,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf
      write(kStdWarn,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
      write(kStdWarn,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf
    END IF
          
    IF ((iNumGasesBAD > 0) .AND. (iReplaceZeroProf > 0)) THEN
      DO iG = 1,iNumGasesBAD
        CALL substitute_tape5_profile_for_climatology(iaBadGasProfile(iG),iNumLevs,raP,raaG_MR)
        iaGasUnits(iaBadGasProfile(iG)) = 10
        raSumCheck(iaBadGasProfile(iG)) = 0.0
        DO iJ = 1,iNumLevs
          raSumCheck(iaBadGasProfile(iG)) = raSumCheck(iaBadGasProfile(iG)) + &
          raaG_MR(iJ,iaBadGasProfile(iG))
        END DO
        write(kStdWarn,*) 'reset ZERO gas prof for gasID ',iaBadGasProfile(iG),' ppmv sum = ',raSumCheck(iaBadGasProfile(iG))
      END DO
    END IF

    13 CONTINUE
    CLOSE(iIOUN2)
    kProfileUnitOpen = -1

    raRTP_TxtInput(1) = rPSurf
    raRTP_TxtInput(2) = rTSurf
    raRTP_TxtInput(3) = rHSurf    !! km
    raRTP_TxtInput(4) = rTophgt   !! km
    raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, downwelling radiation to instr, else upwelling rad to instr

    rPSurf = rPSurf * 100.0       !! change from mb to N/m2
    rHSurf = rHSurf * 1000.0      !! change from km to m

    write(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
    write(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
    write(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf
         
! make sure pressures are decreasing with index ie layers going higher and higher
    IF (raP(1) < raP(2)) THEN
      !!! need to swap!!!
      iMid = ifloor(iNumLevs/2.0)
      DO iL = 1,iMid
        iJ = iNumLevs-iL+1
        rX = raP(iJ)
        raP(iJ) = raP(iL)
        raP(iL) = rX

        rX = raT(iJ)
        raT(iJ) = raT(iL)
        raT(iL) = rX

        DO iG = 1,iNumGases
          raX(iG) = raaG_MR(iJ,iG)
          raaG_MR(iJ,iG) = raaG_MR(iL,iG)
          raaG_MR(iL,iG) = raX(iG)
        END DO
      END DO
    END IF
    write(kStdWarn,*) 'input file : Top alt (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
    write(kStdWarn,*) ' '

    iWriteRTP = +21   !!!! raP,raT,raaG_MR,raAlt are in levels (g/g or ppmv or whatever) and we do have some layer avg info
    IF (iWriteRTP > 0) THEN
      CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf, &
        raP,raT,raaG_MR,raAlt,raPavg,raTavg,iWriteRTP)
    END IF
          
! now change all units first to PPMV (unit 10) and then to Volume Mix Ratio (unit 12)
! actually this is a waste of time as we go from VMR (units 12) to PPMV (units 10) back to VMR (unit 12)
    DO iG = 1,iNumGases
      CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR,+1)
        iaGasUnits(iG) = 12
      DO iL = 1,iNumLevs
        raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
      END DO
    END DO

    RETURN
    END SUBROUTINE 

!************************************************************************
! this reads in TAPE6 which has the integrated layer amounts
    SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
    iNumLays,rPSurf,rTSurf,rHSurf,iNumGases,raPX,raTX,raLayDensityX,raZX, &
    raPressLevelsX,raTPressLevelsX,raAltitudesX, &
    iaG,iaGasUnits,raaG_MRX,rPMin,rPMax,rYear,rLat,rLon)

    IMPLICIT NONE
         
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    CHARACTER(80) :: caPfName
    REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
! output
    INTEGER :: iNumLays,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
!!! these are LAYERS
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas),raLayDensityX(kProfLayer+1)
    REAL :: raZX(kProfLayer+1)
!!! this is LEVELS
    REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)

! local var
    INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,iNumLevs,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
    REAL :: raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,r1,r2,r3,r4,r5,r6,rPTOA
    CHARACTER(80) :: caStrX,caStrY
    CHARACTER(120) :: caStr120,caStr
    CHARACTER(30) :: caStr30
    CHARACTER(1) ::  c1
    CHARACTER(2) ::  c2a,c2b
    CHARACTER(128) :: cX
    INTEGER :: iWriteRTP,iCountGases,iGCnt,iPass,iOffSet
    INTEGER :: i1,i2,i3,i4,i5,i6
! for printing to lblrtm2rtp
    REAL :: raPY(2*kProfLayer),raTY(2*kProfLayer),raaG_MRY(2*kProfLayer,kMaxGas),raZY(2*kProfLayer)
    REAL :: raPavgJunk(2*kProfLayer),raTavgJunk(2*kProfLayer),raNavgJunk(2*kProfLayer),raSumGasAmtTAPE6(2*kProfLayer)
    REAL :: raPZ(2*kProfLayer),raTZ(2*kProfLayer)
          
    DO iL = 1,2*kProfLayer
      raSumGasAmtTAPE6(iL) = 0.0
    END DO
          
    rPmin = +1.0e6
    rPmax = -1.0e+6
    iNXsec = -1

 1070 FORMAT('ERROR! number ',I5,' opening PRFILE path LBLRTM profile file' ,/,A80)
    write(kSTdWarn,*) 'Reading in modified LBLRTM TAPE6 = 5 line summary of TAPE5 + MID TAPE6.....'
    iIOUN2 = kProfileUnit
    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED', IOSTAT=iErrIO)
    IF (iErrIO /= 0) THEN
      iErr=1
      WRITE(kStdErr,1070) iErrIO, caPfname
      CALL DoSTOP
    ENDIF
    kProfileUnitOpen=1

! this reads a modified version of TAPE6 so that the following info is there
! >>>>>>>>>>>>>>>>>>> simple header made by YOU
! short description of your yada yada yada
! rF1 rF2
! rTSurf,rPSurf,rPTOA
! iNumLevs iNumGases iNXsec
! rViewAngle
! >>>>>>>>>>>>>>>>> cut from TAPE5
! TAPE5 INFO
! for i = 1 : iNumlayes
!   read string1
!   read string2
! end do
! END TAPE5 INFO next line is start of TAPE6 info
! >>>>>>>>>>>>>>>>> cut from TAPE6
! caStr starting with 0LAYER and then all the subsequent info

!---------------------- read in header ------------------------
    write(kStdWarn,*) '   ... reading in 5 line summary of TAPE5/TAPE6 file ....'
    READ (iIOUN2,111,ERR=222,END=222) caStr
    READ (iIOUN2,*) rF1,rF2
    write(kStdWarn,*) 'LBLRTM input indicates start/stop wavenumbers are ',rF1,rF2

    READ (iIOUN2,*) rTSurf,rPSurf,rPTOA
    write(kStdWarn,*) 'LBLRTM input indicates STEMP/SPRES and TOA_Pres are ',rTSurf,rPSurf,rPTOA

    READ (iIOUN2,*) iNumLevs,iNumGases,iNXSec
    write(kStdWarn,*) 'TAPE5 implies kCARTA compute ODs at ',iNumLevs,' press/heights, for iNumGases = ',iNumGases
    IF (iNumLevs > kProfLayer) THEN
      write(kStdErr,*) 'though irrelevant for kCARTA calcs, it will not be able to read in the pressures/heights'
      write(kStdErr,*) iNumLevs,kProfLayer
      CALL DoStop
    END IF
    IF (iNumGases+iNXsec > kMaxGas) THEN
      write(kStdErr,*) 'iNumGases+iNXsec > kMaxGas',iNumGases,iNXsec,kMaxGas
      CALL DoStop
    END IF

    READ (iIOUN2,*) rViewAngle
    iLBROutBdryHorP = +1   !! assume start/stop heigt info given

    iNumLays = iNumLevs - 1
    IF (iNumLays > kProfLayer) THEN
      write(kStdErr,*) 'iNumLays > kProfLayer',iNumLays,kProfLayer
      CALL DoStop
    END IF

!---------------------- read in TAPE5 info (for plevs)  -----------
! see subr read_record_2p1p1
    write(kSTdWarn,*) ' '
    write(kStdWarn,*) '   ... reading in relevant edited TAPE 5 info ....'
    READ (iIOUN2,111,ERR=222,END=222) caStr   !!! start TAPE5 info
    DO iL = 1,iNumLevs-1
      READ (iIOUN2,111,ERR=222,END=222) caStr120
      caStr = caStr120(1:40)
      read(caStr,*) raPAvgJunk(iL),raTavgJunk(iL)
      caStr = caStr120(41:120)
      IF (iL == 1) THEN
        read (caStr,*) raAltitudesX(1),raPressLevelsX(1),raTPressLevelsX(1),raAltitudesX(2),raPressLevelsX(2),raTPressLevelsX(2)
      ELSE
        read (caStr,*) raAltitudesX(iL+1),raPressLevelsX(iL+1),raTPressLevelsX(iL+1)
      END IF
      READ (iIOUN2,111,ERR=222,END=222) caStr  !! mixing ratios, not relevant
      !! pV = nRT ==> pdz = n/Vdz RT ==> p dz/RT = q = moles/m3, which can be compared to what is in TAPE6!!!!
      !! change mb to N/m2 and km to m
      raNAvgJunk(iL) = raPavgJunk(iL)*100.0 * (raAltitudesX(iL+1)-raAltitudesX(iL)) * 1000.0
      raNAvgJunk(iL) = raNAvgJunk(iL) / kMGC /raTavgJunk(iL)    !!! this is in mles/m2
      raNAvgJunk(iL) = raNAvgJunk(iL) * 1e-4 * kAvog/1000.0
    END DO
          
    kLBLRTM_TOA = raPressLevelsX(iNumLevs)
    write(kStdWarn,*) 'kLBLRTM_toa = ',kLBLRTM_toa,' mb (when used with flux calcs)'
    write(kStdWarn,*) ' '
          
    iOffSet = kProfLayer+1 - (iNumLevs)
    DO iL = 1,iNumLevs
      kLBLRTM_levelT(iL+iOffSet) = raTPressLevelsX(iL)
    END DO
          
    READ (iIOUN2,111,ERR=222,END=222) caStr   !!! end TAPE5 info
                
!---------------------- read in TAPE6 info ------------------------

! now start reading in the TAPE6 info
! first set gasunits = 1 = molecules/cm2
    write(kSTdWarn,*) ' '
    write(kStdWarn,*) '   ... reading in relevant edited TAPE 6 info ....'
!! we now read every gas 1 .. iNumGases  in edited TAPE6 == molecules/cm2
    DO iG = 1,iNumGases
      iaG(iG) = iG
      iaGasUnits(iG) = 1   !!! assume hardcoded molecules/cm2
    END DO

    READ(iIOUN2,111) caStr
    IF (caStr(1:6) /= '0LAYER') THEN
      write(kStdErr,*) 'oops, expecting string : 0LAYER                          P(MB)       T(K)    ALPHL    ALPHD    ALPHV ...'
      write(kStdErr,*) 'but instead got ',caSTr
      CALL DOStop
    END IF

    DO iL = 1,iNumLays
      READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT
      raZX(iL)   = r1
      raZX(iL+1) = r2
      IF (iL == 1)        rHSurf = r1
      IF (iL == iNumLays) rTopHgt = r2
      !print *,iL,iNumLays,raZX(iL),rP,rT
    END DO

    rPmax = rPSurf
    rPMin = rPTOA
          
    iPass = 1
    iCountGases = 0
    iGCnt = min(7,iNumGases)  !! in first pass, read at most 7 gases
    iGCnt = 8   !! reads first gases, then "all other gases"
!!! might need to re-code this because of OTHER
!!!  H2O           CO2            O3           N2O            CO           CH4            O2            >>> OTHER <<<
!!! because of OTHER need to read at least iNumGases+1  profiles
!!! might need to re-code this because of OTHER
    iGCnt = min(8,iNumGases+1)   !! reads first gases, then "all other gases"
    write(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at FIRST PASS ',iPass
    write(kStdWarn,*) '  because the first 7 are regular profile and eight is "OTHER" '
    READ(iIOUN2,111) caStr      !! blank except for 1 at beginning
    READ(iIOUN2,111) caStr      !! eg LBLRTM    14/07/04  22:24:32
    READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
    READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O
    DO iL = 1,iNumLays
      READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
      DO iG=1,iGCnt
        raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
      END DO
      raZX(iL)   = r1
      raZX(iL+1) = r2
      !raPX(iL) = rP * 100.0  !! change from mb to N/m2
      raPX(iL) = rP          !! keep in mb
      raTX(iL) = rT
      IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
      IF (rPmin > raPX(iL)) rPmin = raPX(iL)
      DO iG=1,iGCnt-1
        raaG_MRX(iL,iG+iCountGases) = raX(iG)
      END DO
      !print *,iL,'read in molgas',raaG_MRX(iL,1),raaG_MRX(iL,2),raaG_MRX(iL,3)
      raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
    END DO
    READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
    READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM
    iCountGases = iCountGases + (iGCnt-1)   !!! because of >>>> OTHER <<<<

    iOffset = kProfLayer - iNumLays
    DO iL = 1,iNumLays
      kLBLRTM_layerTavg(iL+iOffSet) = raTX(iL)
    END DO

 15 CONTINUE
    IF ((iNumGases == iCountGases) .AND. (iNXsec == 0)) GOTO 222     !!! done, just close file!!!!
    IF ((iNumGases == iCountGases) .AND. (iNXsec > 0)) GOTO 23      !!! done with main gases, need to read xsec gases

! need to do this if need to read in more gases
    iPass = iPass + 1
    iGCnt = min(7,iNumGases+1-iCountGases)
!!! might need to re-code this because of OTHER
!!!  H2O           CO2            O3           N2O            CO           CH4            O2            >>> OTHER <<<
!!! because of OTHER need to read at least iNumGases+1  profiles
!!! might need to re-code this because of OTHER
    write(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at pass ',iPass
    READ(iIOUN2,111) caStr      !! eg 1blank
    READ(iIOUN2,111) caStr      !! LBLRTM    14/07/04  22:24:32
    READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
    READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O
    DO iL = 1,iNumLays
      READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
      DO iG=1,iGCnt
        raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
      END DO
      raZX(iL)   = r1
      raZX(iL+1) = r2
      !raPX(iL) = rP * 100.0  !! change from mb to N/m2
      raPX(iL) = rP          !! keep in mb
      raTX(iL) = rT
      IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
      IF (rPmin > raPX(iL)) rPmin = raPX(iL)
      DO iG=1,iGCnt
        raaG_MRX(iL,iG+iCountGases) = raX(iG)
      END DO
      !!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
    END DO
    iCountGases = iCountGases + iGCnt
    READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
    READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM
          
    GOTO 15

 23 CONTINUE    !! finished reading in LAYER amounts for main gases, now read in AUX CUMULATIVE INFO
    DO i1 = 1,iPass
      READ(iIOUN2,111) caStr         !! blank
      READ(iIOUN2,111) caStr      !! ----
      READ(iIOUN2,111) caStr      !!  MIXING RATIOS BY LAYER
      READ(iIOUN2,111) caStr      !!   P(MB)      T(K)   IPATH         H2O           CO2
      DO i2 = 1,iNumLays
        READ(iIOUN2,111) caStr      !! mix ratios
      END DO
      READ(iIOUN2,111) caStr      !! blank
      READ(iIOUN2,111) caStr      !! LBLRTM    14/07/04  22:24:32
    END DO

!!! now ready to read in xsec gases, assume you only need 1 pass
    READ(iIOUN2,111) caStr      !! blank
    READ(iIOUN2,111) caStr      !!  *****  CROSS SECTIONS  *****
    READ(iIOUN2,111) caStr      !! MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
    READ(iIOUN2,111) caStr      !! P(MB)      T(K)   IPATH     F11           F14           F22
    caStr = caStr(56:120)
    CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
    IF (iNXsec > 7) THEN
      write(kStdWarn,*) 'oops assumed at most 7 xsec gases, so only 1 pass ... need to modify code for > 7 gases'
      CALL DOStop
    END IF
    DO iL = 1,iNumLays
      READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iNXsec)
      DO iG=1,iGCnt
        raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
      END DO
      raZX(iL)   = r1
      raZX(iL+1) = r2
      !raPX(iL) = rP * 100.0  !! change from mb to N/m2
      raPX(iL) = rP          !! keep in mb
      raTX(iL) = rT
      IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
      IF (rPmin > raPX(iL)) rPmin = raPX(iL)
      DO iG=1,iNXsec
        raaG_MRX(iL,iG+iNumGases) = raX(iG)
      END DO
      !!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
    END DO
          
 222 CONTINUE
    CLOSE(iIOUN2)
    kProfileUnitOpen = -11
    iNumGases = iNumGases + iNXsec
          
 111 FORMAT(A120)
 112 FORMAT(I3,2(' ',F10.3), 2(' ',F10.3),5(' ',ES12.5),1(' ',F12.5))
 113 FORMAT(A122)
     
    write(kStdWarn,*) '... finshed reading in TAPE5/TAPE6 combo ....'
    write(kStdWarn,*) ' '
    write(kSTdWarn,*) '     rTSurf,rPSurf,rPTOA = ',rTSurf,rPSurf,rPTOA
    cX = &
    ' Index   Plev      Tlev       Pavg       Tavg     WV(amt)       CO2(amt)   O3(amt)    || q=pdz/RT   totalQ(TAPE6)    Ratio'
    write(kStdWarn,113) cX
    cX = &
    '--------------------------------------------------------------------------------------||-----------------------------------'
         
    write(kStdWarn,113) cX
    DO iL = 1,iNumLays+1
      write(kSTdWarn,112) iL,raPressLevelsX(iL),raTPressLevelsX(iL),raPX(iL),raTX(iL),(raaG_MRX(iL,iJ),iJ=1,3), &
        raNAvgJunk(iL),raSumGasAmtTAPE6(iL),raNAvgJunk(iL)/raSumGasAmtTAPE6(iL)
    END DO
    write(kStdWarn,113) cX
          
    raRTP_TxtInput(1) = rPSurf
    raRTP_TxtInput(2) = rTSurf
    raRTP_TxtInput(3) = rHSurf    !! km
    raRTP_TxtInput(4) = rTophgt   !! km
    raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, then downwell rad, else upwell radn

    rPSurf = rPSurf * 100.0       !! chnage from mb to N/m2
    rHSurf = rHSurf * 1000.0      !! change from km to m

    write(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
    write(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
    write(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf
         
! make sure pressures are decreasing with index ie layers going higher and higher
    IF (raPX(1) < raPX(2)) THEN
      !!! need to swap!!!
      iMid = ifloor(iNumLays/2.0)
      DO iL = 1,iMid
        iJ = iNumLays-iL+1
        rX = raPX(iJ)
        raPX(iJ) = raPX(iL)
        raPX(iL) = rX

        rX = raTX(iJ)
        raTX(iJ) = raTX(iL)
        raTX(iL) = rX

        DO iG = 1,iNumGases
          raX(iG) = raaG_MRX(iJ,iG)
          raaG_MRX(iJ,iG) = raaG_MRX(iL,iG)
          raaG_MRX(iL,iG) = raX(iG)
        END DO
      END DO
    END IF
    write(kStdWarn,*) 'input file : Highest altitude (lowest <press>) = ',raPX(iNumLays),' mb at layer ',iNumLays
    write(kStdWarn,*) 'input file : Highest altitude (lowest plev)    = ',rPTOA,' mb at level ',iNumLays+1
    write(kStdWarn,*) ' '

    iWriteRTP = +1   !!! this is in layers
    IF (iWriteRTP > 0) THEN
      DO iL = 1,iNumLevs
        raPZ(IL) = raPX(iL)             !!! average pressure
        raTZ(iL) = raTPressLevelsX(iL)  !!! levels pressure
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        raTY(IL) = raTX(iL)             !!! layers temp
        raZY(IL) = raZX(iL)             !!! layers alt
        raPY(IL) = raPressLevelsX(iL)        !!! always need levels info for p.plevs
        raZY(IL) = raAltitudesX(iL)        !!! levels alt
        DO iG = 1,iNumGases
          raaG_MRY(iL,iG) = raaG_MRX(iL,iG)
        END DO
      END DO
      CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf, &
        raPY,raTY,raaG_MRY,raZY,raPZ,raTZ,iWriteRTP)
    END IF

! now change all units to MR
!      DO iG = 1,iNumGases
!        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR,+1)
!        DO iL = 1,iNumLevs
!          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
!          iaGasUnits(iG) = 12
!        END DO
!      END DO

    RETURN
    END SUBROUTINE 

!************************************************************************
! this subroutine tacks on Standard Profile above what the user has given, to fill in amounts till 0.005 mb
! also if kCARTA has earlier determined there should be more gases than in user input profile (eg suppose
! user only gave WV/O3 but kCARTA also wants other gases), then this is the place additional profiles added in
! also can do the CO2 = 370 ppm +(yy-2002)*2.2
! also can do the CH4 = 1.7 ppm +(yy-2002)*0.01
    SUBROUTINE Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR, &
    rPMin,rPMax,rYear)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input var
    INTEGER :: iNumLevs
    REAL :: rYear
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! input/output var
    REAL :: rPmin,rPmax
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas)
    INTEGER :: iaG(kMaxGas),iaGasUnits(kMaxGas),iNumGases

! local var
!! individual reference profiles, at kMaxLayer layers, in terms of MR = PPress/Press
    INTEGER :: iaG0(kMaxGas),iNumGases0,i2,iFound,iaGasUnits0(kMaxGas),iIndex
    REAL :: raaR100MR(kMaxLayer+10,kMaxGas),raR100Temp(kMaxLayer+10),raR100Press(kMaxLayer+10)
    INTEGER :: iL,iG,iK,iAbove,iMerge,iRefLevels,laysORlevs
    REAL :: rToffset,rJunk,raOffset(kMaxGas),raJunk(kMaxLayer+10)
    REAL :: rCO2ppmv,rAdjust,rCH4ppbv
    CHARACTER(3) :: ca3
    CHARACTER(50) :: FMT
    INTEGER :: iEarth,iMars,iaEarth(8),iaMars(6),iaPlanetMolecules(kMaxGas),iPlanetMolecules
          
    DATA (iaEarth(iG),iG=1,8) /01,02,03,04,05,06,09,12/      !! 8 important Earth atm molecules = HITRAN
                                                             !! WV,CO2,O3,N2O,CO,CH4,SO2,HNO3
    DATA (iaMars(iG),iG=1,6)  /01,02,03,04,05,06/            !! 6 important Mars  atm molecules = UMBC
                                                             !! WV,CO2,N2,O2,CO,NO

    iEarth = 8    
    iMars = 6

    iaPlanetMolecules = -1
    IF (kPlanet == 03) THEN
      iPlanetMolecules = iEarth
      DO iG = 1,iEarth
        iaPlanetMolecules(iG) = iaEarth(iG)
      END DO
    ELSEIF (kPlanet == 04) THEN
      iPlanetMolecules = iMars
      DO iG = 1,iMars
        iaPlanetMolecules(iG) = iaMars(iG)
      END DO
    ELSE
      write (kStdErr,*) 'oops right now can only handle Earth or Mars. Sorry, Mork'
      Call DoStop
    END IF
          
! first read in reference profile, and if needed info for other gases
!     read in additional gases eg from ERA/ECM we get gas 1,3 but we would like gas 1,2,3,4,5,6,8,12
!     read in ref profiles (eg for Gas2 since Mars, Earth, Venus all have CO2 in the atm; get the Standard VMR)
    iNumGases0 = iNumGases
    DO iL = 1,iNumGases0
      iaG0(iL) = iaG(iL)
      iaGasUnits0(iL) = iaGasUnits(iL)
    END DO
    laysORlevs = +1   !! orig, read in reference P,PP and get mix ratio using PP/P for each gas
    laysORlevs = -1   !! new,  read in one of 6 AFGL P/T/ppmv level profiles
    ! a3 = kcaLevsRefProf(1:3)
    ca3 = 'DNE'
    iIndex = index(kcaLevsRefProf,'DNE')
    IF ((laysORlevs == -1) .AND. (iIndex > 0)) THEN
      write(kStdWarn,*) 'oops : do not have a set of Levels Reference Profiles; code is setting laysORlevs = +1'
      laysORlevs = +1
    END IF
          
    Call ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,iaPlanetMolecules,iPlanetMolecules, &
      iaG,iaGasUnits,iNumGases,raR100Press,raR100Temp,raaR100MR,iRefLevels,laysORlevs)

    IF (iNumGases0 /= iNumGases) THEN
      write(kStdWarn,*)'Orig num of gases, new num of gases = ',iNumGases0,iNumGases
      !! new gases are all at VMR (dimensionless fraction) = gasunit 12
      !! these new profiles are AFGL profiles from glatm, and NOT at sonde or AIRS 101 levels
      !! so they need to be interpolated from (raR100Press,raaR100MR) into raaG_VMR (at sonde P levels between rPmax and rPmin)
      CALL InterpNewGasProfiles_to_InputPlevels(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels, &
        iNumLevs,rPmin,rPmax,raP,raaG_VMR)
    END IF

! then find KCARTA levels above which there is NO user supplied info
    FMT = '(A,F12.5,F12.5,I3)'
    write(kStdWarn,FMT)' user supplied profile : Pmax(mb),Pmin(mb),numlevs = ',rPMax,rPMin,iNumLevs
    write(kStdWarn,FMT)' kCARTA Pav Database   : Pmax(mb),Pmin(mb),numlevs = ', &
    raR100Press(1)/100.0,raR100Press(iRefLevels)/100.0,iRefLevels
          
    iAbove = kMaxLayer
    iAbove = iRefLevels
    10 CONTINUE
    IF ((raR100Press(iAbove)/100.0 < rPMin) .AND. (iAbove > 1)) THEN
      iAbove = iAbove - 1
      GOTO 10
    END IF
    iAbove = iAbove + 1
    IF (iAbove > iRefLevels) THEN
      FMT = '(A,I3,A,I3)'
      write(kStdWarn,FMT) 'tacking on ref  profile info from level ',iAbove,' to ',iRefLevels
      write(kStdWarn,*) 'Hmm ... no need to do this ... looks like user supplied levels profile past 0.005 mb'
      GOTO 123
    END IF

    FMT = '(A,I3,A,I3)'
    write(kStdWarn,FMT) 'tacking on ref profile info from level ',iAbove,' to ',iRefLevels
    write(kStdWarn,FMT) 'which should extend user supplied level info from level ',iNumLevs,' to ', &
    iNumLevs + (iRefLevels-iAbove+1)
          
! now do linear interpolation from iAbove pressure, down to min(raP)
    Call r_sort_loglinear(raR100Press,raR100Temp,iRefLevels,rPMin*100.0,rJunk)
    rToffset = raT(iNumLevs) - rJunk
    FMT = '(A,I3,A,I3)'
    write(kStdWarn,FMT)'tacking on info from kCARTA Pav Dtabase, layers ',iAbove,' to ',iRefLevels
    write(kStdWarn,*)'ToffSet = ',rToffset

! also need to figure out the gas multiplier offset
    DO iG = 1, iNumGases
      DO iL = 1,iRefLevels
        raJunk(iL) = raaR100MR(iL,iG)
      END DO
      Call r_sort_loglinear(raR100Press,raJunk,iRefLevels,rPMin*100.0,rJunk)
      !! need to do a CHANGE OF UNITS to ppmv!!!!!
      raoffset(iG) = raaG_VMR(iNumLevs,iG)/rJunk
      IF (raoffset(iG) < 0) THEN
        raoffset(iG) = 1.0
      END IF
      write(kStdWarn,*) 'Multiplier for gasID = ',iaG(iG),' units = ',iaGasUnits(iG),' = ',raoffset(iG)
    END DO
          
    iMerge = +1    !! in other words, 3 points away from iMerge, we bump up/down the multiplier to 1.000
    DO iL = iAbove,iRefLevels
      iNumLevs = iNumLevs + 1
      raP(iNumLevs) = raR100Press(iL)/100.0

      !aT(iNumLevs) = raR100Temp(iL) + rToffset
      IF ((iL-iAbove) > iMerge) THEN
        rAdjust = 0.0  !! make the temp adjustment effectively 0
      ELSE
        rJunk = (0 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 0, inside 4 levels; this is slope
        rAdjust = rJunk * (iL-iAbove) + rToffset
        rAdjust = rToffset  !!!!!!!! ????
      END IF
      raT(iNumLevs) = raR100Temp(iL) + rAdjust

      DO iG = 1,iNumGases
        !raaG_VMR(iNumLevs,iG) = raaR100MR(iL,iG) * raOffset(iG)
        IF ((iL-iAbove) > iMerge) THEN
          rAdjust = 1  !! make the MR adjustment effectively 1
        ELSE
          rJunk = (1 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 1, inside 4 levels; this is slope
          rAdjust = rJunk * (iL-iAbove) + raOffset(iG)
          rAdjust = raOffset(iG)   !!!! ????
          ! print *,iG,iL,iL-iAbove,raOffset(iG),rAdjust
        END IF
        raaG_VMR(iNumLevs,iG) = raaR100MR(iL,iG) * rAdjust
      END DO
    END DO
    write(kStdWarn,*) 'Have extended user supplied info to level',iNumLevs
          
    123 CONTINUE
          
! >>>>>>>>>>>>>>>>>>>>>>>>> CO2 adjust; growth rate = 2 ppmv/yr
    rCO2ppmv = 370.0
    IF ((rYear > 1970) .AND. (kPlanet == 03)) THEN
      !! should be dt = (yy-2002) + mm/12;  rCO2ppmv = 370 + 2*dt;
      rCO2ppmv = 370.0 + (rYear-2002)*2.0

      !! now see if CO2 was originally there
      iFound = -1
      i2 = 1
 20   CONTINUE
      IF (iaG0(i2) == 2) THEN
        iFound = +1
      ELSEIF (i2 < iNumGases0) THEN
        i2 = i2 + 1
        GOTO 20
      END IF

      IF (iFound > 0) THEN
        !! user had CO2 in supplied prof, no need to adjust MR
        write(kStdWarn,*) 'User has CO2 in supplied profile, so not adjusting bottom levels'
      ELSE
        !! user did not have CO2 in there, yes need to adjust MR
        i2 = 1
        iFound = -1
 25     CONTINUE
        IF (iaG(i2) == 2) THEN
          iFound = +1
        ELSEIF (i2 < iNumGases) THEN
          i2 = i2 + 1
          GOTO 25
        END IF
        IF (iFound < 0) THEN
          write(kStdErr,*) 'hmm, expected to find CO2!!! ooops'
          CALL DoStop
        END IF
        rJunk = 0.0
        DO iL = 3,7
          rJunk = rJunk + raaG_VMR(iL,i2)
        END DO
        rJunk = rJunk/5.0
        write(kStdWarn,*) '  reference CO2 profile that was read in has CO2 MR = ',rJunk*1.0e6,' ppmv'
        write(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCO2ppmv,' ppmv'
        DO iL = 1,iNumLevs
          raaG_VMR(iL,i2) = raaG_VMR(iL,i2) * rCO2ppmv/(rJunk*1.0e6)
        END DO
      END IF
    END IF

! >>>>>>>>>>>>>>>>>>>>>>>>> CH4 adjust; growth rate = 0.01 ppmb/yr
    rCH4ppbv = 1.7
    IF ((rYear > 1970) .AND. (kPlanet == 03)) THEN
      rCH4ppbv = 1.7 + (rYear-2002)*0.01
      rCH4ppbv = 1.8   !! hard coded into klayers right now

      !! now see if CH4 was originally there
      iFound = -1
      i2 = 1
 30   CONTINUE
      IF (iaG0(i2) == 6) THEN
        iFound = +1
      ELSEIF (i2 < iNumGases0) THEN
        i2 = i2 + 1
        GOTO 30
      END IF

      IF (iFound > 0) THEN
        !! user had CH4 in supplied prof, no need to adjust MR
        write(kStdWarn,*) 'User has CH4 in supplied profile, so not adjusting bottom levels'
      ELSE
        !! user did not have CH4 in there, yes need to adjust MR
        i2 = 1
        iFound = -1
 35     CONTINUE
        IF (iaG(i2) == 6) THEN
          iFound = +1
        ELSEIF (i2 < iNumGases) THEN
          i2 = i2 + 1
          GOTO 35
        END IF
        IF (iFound < 0) THEN
          write(kStdErr,*) 'hmm, expected to find CH4!!! ooops'
          CALL DoStop
        END IF
        rJunk = 0.0
        DO iL = 3,7
          rJunk = rJunk + raaG_VMR(iL,i2)
        END DO
        rJunk = rJunk/5.0
        write(kStdWarn,*) '  reference CH4 profile that was read in has CH4 MR = ',rJunk*1.0e6,' ppbv'
        write(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCH4ppbv,' ppbv'
        DO iL = 1,iNumLevs
          raaG_VMR(iL,i2) = raaG_VMR(iL,i2) * rCH4ppbv/(rJunk*1.0e6)
        END DO
      END IF
    END IF

!>>>>>>>>>>>>>>>>>>>>>>>>>

!! check all MR larger than 0
    DO iG = 1, iNumGases
      DO iL = 1,kMaxLayer
        raaG_VMR(iL,iG) = max(raaG_VMR(iL,iG),0.0)
      END DO
    END DO
          
    RETURN
    end SUBROUTINE Tack_on_profile

!************************************************************************
! this subroutine changes the input levels profiles units from VMR to ppmv
!    >>>> and if planet Earth, adjusts mix ratios (wet water) <<<<
! then convers back to VMR
    SUBROUTINE adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer)
! input/output
    REAL :: raaG_VMR(2*kProfLayer,kMaxGas)

! local
    INTEGER :: iL,iG,iDry2Wet,iWat,iFound,iQuiet
    REAL :: rJunk1,rJunk2
    CHARACTER(40) :: FMT
          
    iFound = -1
    iWat = 1
 10 CONTINUE
    IF (iaG(iWat) == 1) THEN
      iFound = +1
    ELSE
      iWat = iWat + 1
      IF (iWat <= iNumGases) THEN
        GOTO 10
      END IF
    END IF

    IF ((kPlanet == 3) .AND. (iFound == 1)) THEN
      write(kStdWarn,*) 'Planet = 3, EARTH .. gas ',iWat,' corresponds to water'
      IF (iWat /= 1) THEN
        write(kStdErr,*) 'Planet = 3, EARTH .. found water but it should have been g1, not g',iWat
        CALL DoStop
      END IF
    ELSEIF ((kPlanet == 3) .AND. (iFound /= 1)) THEN
      write(kStdErr,*) 'Planet = 3, EARTH .. no input gas corresponds to water'
      Call DoStop
    END IF

!! change to ppmv
    FMT = '(A,/A,/A,/A)'
    write(kStdWarn,FMT) 'FINAL PROCESSING : Read in reference gas profiles, text input levels profile etc', &
    '     (ref. profiles to augment input prof to 0.005 mb and/or add standard panet gases)', &
    '  converting from VMR to ppmv so we can do wet water fix', &
    'After this we can do the levels2layers integration'
    iQuiet = +1
    DO iG = 1,iNumGases
      rJunk1 = raaG_VMR(1,iG)
      IF (iaGasUnits(iG) /= 10) THEN
        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_VMR,iQuiet)
        rJunk2 = raaG_VMR(1,iG)
        write(kStdWarn,800) iG,iaG(iG),iaGasUnits(iG),rJunk1,rJunk2
      ELSE
        write(kStdWarn,801) iG,iaG(iG),iaGasUnits(iG),rJunk1
      END IF
    END DO
 800 FORMAT('index=',I3,' gasID=',I3,' InputGasUnits=',I3,' origQ(level1)=',ES12.6,' finalQ_ppmv(level1)=',ES12.6)
 801 FORMAT('index=',I3,' gasID=',I3,' InputGasUnits=',I3,' origQ(level1)=',ES12.6,' already in ppmv')
     
    iDry2Wet = -1 !! do not adjust dry air to wet air MR
    iDry2Wet = +1 !!        adjust dry air to wet air MR  DEFAULT

    ! IF (kRTP == -20) iDry2Wet = -1

    IF (iDry2Wet < 0) THEN
      !! simple change from ppmv (gas units 10) to VMR (which is gas units 12)

      write(kStdWarn,*) 'converting back to ppmv from VMR, not doing wet water fix'
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0e6
        END DO
      END DO

    !! need to make sure gas units code == 10.11.12.20.21 before trying this; see toppmv.f in klayers
    ELSEIF ((iDry2Wet > 0) .AND. (kPlanet == 3) .AND. (iaG(1) == 1)) THEN

      !!! first convert WV mixing ratios
      !!! however, klayars assumes the user has provided WET mixing ratios, so NO NEED to do this
      !        DO iG = 1,1
      !          DO iL = 1,iNumLevs
      !            print *,iL,raaG_VMR(iL,iG),1e6*raaG_VMR(iL,iG) / (1.0e6+raaG_VMR(iL,iG))
      !            raaG_VMR(iL,iG) =  1.0e6 * raaG_VMR(iL,iG) / (1.0e6 + raaG_VMR(iL,iG))
      !          END DO
      !        END DO

      !!! then use WV mixing ratios to fix other gas mixing ratios if they are 10.11.12.20.21
      DO iG = 2,iNumGases
        IF ((iaGasUnits(iG) >= 10) .AND. (iaGasUnits(iG) <= 21)) THEN
          DO iL = 1,iNumLevs
            raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) * (1.0e6 - raaG_VMR(iL,1))/1e6
          END DO
        END IF
      END DO

      !! finally convert ppmv --> VMR
      write(kStdWarn,*) 'converting back to ppmv from VMR, after doing wet water fix'
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0e6
        END DO
      END DO

    ELSE
      write(kStdErr,*) 'Need iDry2Wet == +/- 1'
      Call DoStop
    END IF

!! check all MR larger than 0
    DO iG = 1,iNumGases
      DO iL = 1,kMaxLayer
        raaG_VMR(iL,iG) = max(raaG_VMR(iL,iG),0.0)
      END DO
    END DO
     
    RETURN
    end SUBROUTINE adjustVMR_Earth

!************************************************************************
! this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
    SUBROUTINE InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta, &
    PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input var
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas)
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rPmaxKCarta
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas)
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! output var
    INTEGER :: iLowestLev
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_VMRX(kProfLayer+1,kMaxGas)
    REAL :: rInJunk,rOutJunk

! local var
    INTEGER :: iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
    REAL :: raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_VMRX2(kProfLayer+1,kMaxGas)

    write(kStdWarn,*) ' ----->>>> using default AIRS 101 levels for layers integration'
          
! see which kCARTA pressure level is equal to or greater than surface pressure
    write(kStdWarn,*) '--------------------------------------------------------------------------'
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' >>> prep : klayers within kcarta : interp usergrid onto 101 levels grid >>>'
          
    iLowestLev = kMaxLayer
    10 CONTINUE
    IF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) < rPSurf) .AND. (iLowestLev > 1)) THEN
      iLowestLev = iLowestLev - 1
      GOTO 10
    ELSEIF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) < rPSurf) .AND. (iLowestLev == 1)) THEN
      write (kStdErr,*) 'PSurf = ',rPSurf,' and Max Kcarta Database Plev = ',rPmaxKCarta
      Call DoStop
    ELSEIF (PLEV_KCARTADATABASE_AIRS(iLowestLev) >= rPSurf) THEN
      write(kStdWarn,*)'lowest level = ',iLowestLev
    END IF

! fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
    DO iL = 1,kMaxLayer-iLowestLev+1
      iNumUse = iL
      raPX(iL) = PLEV_KCARTADATABASE_AIRS(iLowestLev+iL)/100.0  ! input press in mb, so change plev_kcarta to mb
    END DO
          
! now interpolate T and MR onto this grid
    Call r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
    DO iG = 1,iNumGases
      DO iL = 1,iNumLevs
        raTemp(iL) = raaG_VMR(iL,iG)
      END DO
      Call r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
      DO iL = 1,iNumUse
        raaG_VMRX(iL,iG) = raTempX(iL)
      END DO
    END DO
          
    1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
     
! need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
    iHigh = iNumUse
 20 CONTINUE
    IF (raP(iNumLevs) > raPX(iHigh)) THEN
      iHigh = iHigh - 1
      GOTO 20
    END IF
    iHigh = iHigh + 1
    IF (iHigh <= iNumUse) THEN
      write(kSTdWarn,*) ' '
      write(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
      write(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
      write(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' to ',iNumUse
      DO iL = iHigh,iNumUse
        rInJunk = raPX(iL)
        Call r_sort_loglinear(raP,raT,iNumLevs,rInJunk,rOutJunk)
        raTX(iL) = rOutJunk

        DO iG = 1,iNumGases
          DO iX = 1,iNumLevs
            raTemp(iX) = raaG_VMR(iX,iG)
          END DO
          Call r_sort_loglinear(raP,raTemp,iNumLevs,rInJunk,rOutJunk)
          raaG_VMRX(iL,iG) = rOutJunk
        END DO
      END DO
    END IF

! now sort the gasIDs into increasing order, and do the same for the mixing ratios
    DO iL = 1,iNumGases
      raArr(iL) = iaG(iL) * 1.0
    END DO
    CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)

! assign temp arrays
    DO iG = 1,iNumGases
      iaG2(iG) = iaG(iG)
    END DO
    DO iL = 1,iNumUse
      DO iG = 1,iNumGases
        raaG_VMRX2(iL,iG) = raaG_VMRX(iL,iG)
      END DO
    END DO

! sort according to iaIndex
    DO iG = 1,iNumGases
      iaG(iG) = iaG2(iaIndex(iG))
    END DO
    DO iL = 1,iNumUse
      DO iG = 1,iNumGases
        raaG_VMRX(iL,iG) = raaG_VMRX2(iL,iaIndex(iG))
      END DO
    END DO

    write(kStdWarn,*) '--------------------------------------------------------------------------'
          
    RETURN
    end SUBROUTINE InterpUser2kCARTA

!************************************************************************
! this subroutine does the klayers integrations, integrals over P
! so this is more similar to kLAYERs, CRTM
! this is newer code; it ignores surface as "lowest" level and only uses the user level info
! this is for AIRS 101 levels

! >>>>>>>> raPout comes in mb
! >>>>>>>> and goes out in N/m2

    SUBROUTINE DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon, &
    PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot, &
    raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: rHSurf,rPSurf,rTSurf,PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    REAL :: DATABASELEVHEIGHTS(kMaxLayer+1),rfracBot,rLat,rLon
    INTEGER :: iLowestLev,iNumGases,iaG(kMaxGas)
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_VMRX(kProfLayer+1,kMaxGas)
! output
    REAL :: raTout(kProfLayer),raAmountOut(kProfLayer),raZout(kProfLayer+1),raPout(kProfLayer)
    REAL :: raaQout(kProfLayer,kMaxGas),raaPartPressOut(kProfLayer,kMaxGas)

! local
    INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
    REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
    REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
    REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
    REAL :: raXYZPress(kMaxlayer+1+1)
    REAL :: raXYZTemp(kMaxlayer+1+1)
    REAL :: raaXYZ_VMR(kMaxlayer+1+1,kMaxGas)
    REAL :: raXYZ_VMRwater(kMaxlayer+1+1)
    REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
    REAL :: zWoo,rConvertQ,rRH
    INTEGER :: iWoo

    rConvertQ = 6.022141e23/1e4     ! multiply by this to change moles/m2 to molecules/cm2
    rConvertQ = kAvog/1000/1.0e4    ! multiply by this to change moles/m2 to molecules/cm2
    rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor

    raaQout = 0.0

    IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
      write (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
      Call DoStop
    END IF

    iNFine = 10
    iNFine = 200
    iNFine = 50  !! default

    write(kStdWarn,*) '  '
    write(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers integration >>>>>>>>>>>>>>>>>>>'

! >>>>>>>>>>
!! set up temporary arrays/matrices for interping onto sublevels
    iUpperLev = kMaxLayer-iLowestLev+1
!! set levels 1 .. iNumLevs
    DO iL = 1,iUpperLev
      raXYZPress(iL) = raPX(iL)
      raXYZTemp(iL) = raTX(iL)
      DO iG = 1,iNumGases
        raaXYZ_VMR(iL,iG) = raaG_VMRX(iL,iG)
        IF (iG ==  1) raXYZ_VMRwater(iL) = raaG_VMRX(iL,iG)
      END DO
    END DO
! >>>>>>>>>>

! now start integrating
! lowest layer is special case

! http://cimss.ssec.wisc.edu/~paulv/
! http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
!     /home/sergio/KCARTA/DOC/sci_klayers.txt

    1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
     
    z = rHSurf !! need to kludge this
    zWoo = rHSurf
    iWoo = -1

    DO iL = iLowestLev,iLowestLev
      raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
      raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
      raPout(iL) = raPout(iL)/100.0    !!! change N/m2 to mb
           
      z0 = z
      iJ = iL - iLowestLev + 1

      amount = 0.0
      rTsum  = 0.0
      rWgt = 0.0

      !!! >>>>>>  this is starting out at         rP = PLEV_KCARTADATABASE_AIRS(iL)        rT = rTSurfx  >>>>>>>>
      rP = PLEV_KCARTADATABASE_AIRS(iL)
      rP = rP/100.0    !!! change N/m2 to mb
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,iUpperLev,rP,rT)

      !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = log(a/c)
      dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))
      dlnp = dlnp / (iNFine)

      !! information for current (sub)level
      rP_n = rP                             !! current sublev press
      rT_n = rT                             !! current sublev temp
      rMR_water_n = raXYZ_VMRwater(1)       !! current water MR
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP,rMR_water_n) !! current water MR

      !! information for next (sub)level
      rP_np1 = log(rP_n) - dlnp
      rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)  !! next sublev MRw

      IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
        iWoo = +1
        zWoo   = z
      END IF

      DO iLoop = 1,iNFine

        rMolarMass_n   = kAtmMolarMass
        rMolarMass_np1 = kAtmMolarMass
        !!! if Earth, adjust molecular mass for displacement by water
        IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
        IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
        rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
        rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

        density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
        density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

        IF (kPlanet /= 3) THEN
          gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
          gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
        ELSE
          gravity = grav_earth(z,0.0,0.0,rLat,rLon)
        END IF

        Tav   = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
        SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
        dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

        rTsum = rTsum + Tav * (density_n + density_np1)/2
        rWgt  = rWgt + (density_n + density_np1)/2
                                 
        Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
        damount = Pav/Tav/kMGC * dz
        ! print *,iLoop,PLEV_KCARTADATABASE_AIRS(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
        amount = amount + damount

        DO iG = 1,iNumGases
          DO iCnt = 1,iUpperLev
            raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
          END DO
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n  )
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1)
          raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
        END DO

        !!! now update for next iteration
        z = z + dz
        rP_n = rP_np1
        rT_n = rT_np1
        rMR_water_n = rMR_water_np1
        rP = rP_n

        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                               !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

        IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
          iWoo = +1
          zWoo   = z
        END IF

      END DO

      raTout(iL)      = rTsum/rWgt                 ! <<< need TAV
      raAmountOut(iL) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
      raZout(iL+1)    = z

    END DO

! displace everything  z
    iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
    raZout(iL)   = -(zWoo - z0)
    raZout(iL+1) = +(z-zWoo)
    write(kStdWarn,*) '  need to displace lowest layer heights'
    write(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
    write(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',PLEV_KCARTADATABASE_AIRS(iL),rPSurf, &
    PLEV_KCARTADATABASE_AIRS(iL+1)
    z = z - zWoo

! c      !no need to adjust this amount for book-keeping, make consistent with eg n_rad_jac_scat.f since we did FULL layer
! c      raAmountOut(iL) = raAmountOut(iL)/rFracBot
! c      DO iG = 1,iNumGases
! c        raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
! c      END DO
          
! go to TOA
    DO iL = iLowestLev + 1,kProfLayer
      z0 = z

      !! compute avg pressure of layer
      raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
      raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
      raPout(iL) = raPout(iL)/100.0    !!! change N/m2 to mb
            
      iJ = iL - iLowestLev + 1
      rP = PLEV_KCARTADATABASE_AIRS(iL)
      rP = rP/100.0    !!! change N/m2 to mb
      rT = raXYZTemp(iJ)
              
      amount = 0.0
      rTsum  = 0.0
      rWgt = 0.0

      !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
      dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))
      dlnp = dlnp / (iNFine)

      !! information for current (sub)level
      rP_n = rP                              !! current sublev press
      rT_n = rT                              !! current sublev temp
      rMR_water_n = raXYZ_VMRwater(iJ)        !! current water MR

      !! information for next (sub)level
      rP_np1 = log(rP_n) - dlnp
      rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

      DO iLoop = 1,iNFine

        rMolarMass_n   = kAtmMolarMass
        rMolarMass_np1 = kAtmMolarMass
        !!! if Earth, adjust molecular mass for displacement by water
        IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
        IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
        rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
        rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

        density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
        density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

        IF (kPlanet /= 3) THEN
          gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
          gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
        ELSE
          gravity = grav_earth(z,0.0,0.0,rLat,rLon)
        END IF

        Tav = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
        SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
        dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

        rTsum = rTsum + Tav * (density_n + density_np1)/2
        rWgt  = rWgt + (density_n + density_np1)/2
                                 
        Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
        damount = Pav/Tav/kMGC * dz
        amount = amount + damount
                  
        DO iG = 1,iNumGases
          DO iCnt = 1,iUpperLev
              raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
          END DO
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n  )
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1)
          raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
          ! q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
          ! q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
          ! raaQout(iL,iG) = raaQout(iL,iG) + (q_n + q_np1)/2 * rConvertQ
        END DO

        !!! now update for next iteration
        z = z + dz
        rP_n = rP_np1
        rT_n = rT_np1
        rMR_water_n = rMR_water_np1
        rP = rP_n

        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

      END DO

      raTout(iL)      = rTsum/rWgt             ! <<< need TAV
      raAmountOut(iL) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
      raZout(iL+1)    = z

      !write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
      !$        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
    END DO
          
    write(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
    !! raZout is in meters
    !! raaQout is in moles/m2
    DO iL = iLowestLev,kProfLayer
      rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
      rT = raTout(iL)
      raPout(iL) = raPout(iL) * 100.0    !! change mb to N/m2
      DO iG = 1,iNumGases
        raaPartPressOut(iL,iG) = raaQout(iL,iG) / raAmountOut(iL) * raPout(iL)
      END DO
    END DO

! testing
!          !!! this automatically puts partial pressure in ATM, assuming
!          !!! gas amount in kilomoles/cm2, length in cm, T in kelvin
!          !!!note "j"!!!
!      DO iL = iLowestLev,kProfLayer
!        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
!        rT = raTout(iL)
!        DO iG = 1,iNumGases
!          rAmt = raaQout(iL,iG)/kAvog
!          rPP  = rAmt*1.0e9*kMGC*rT / (rThick*kAtm2mb*100.0)
!          IF (iG .EQ. 1) THEN
!            print *,iL,iG,raaPartPressOut(iL,iG)/(kAtm2mb*100).0,rPP,rThick,
!     $             raaPartPressOut(iL,iG)/rPP/(kAtm2mb*100.0),raAmountOut(iL),raaQout(iL,iG)
!          END IF
!        print *,iL,raPout(iL)/100,raTout(iL),raAmountOut(iL),
!     $        raaQout(iL,1),raaQout(iL,2),raaQout(iL,3),raaQout(iL,4),raaQout(iL,5),
!     $        raaG_VMRX(iL-iLowestLev+1,1),raaG_VMRX(iL-iLowestLev+1,2),raaG_VMRX(iL-iLowestLev+1,3),
!     $        raaG_VMRX(iL-iLowestLev+1,4),raaG_VMRX(iL-iLowestLev+1,5)
!        END DO
!      END DO
! testing

    IF ((iLowestLev-1) >= 1) THEN
      DO iL = 1,iLowestLev-1
        raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
        raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
        raTout(iL) = kStempMin
        DO iG = 1,iNumGases
          raaPartPressOut(iL,iG) = 0.0
          raaQout(iL,iG) = 0.0
          raaG_VMRX(iL,iG) = 0.0
        END DO
      END DO
    END IF

    IF (kPlanet == 3) THEN
      !! show the RH for gas 1
      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Checking Relative Humidity (integrate over AIRS 101 levels)'
      write(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Q6'
      write(kStdWarn,113) '------------------------------------------------------------------------------------------'
      DO iL = iLowestLev,kProfLayer
        z    = (raZout(iL) + raZout(iL+1))/2/1000
        zWoo = (raZout(iL+1) - raZout(iL))/1000
        rP   = raPout(iL)            !! N/m2
        rT   = raTout(iL)
        rPP  = raaPartPressOut(iL,1) !! N/m2
        rSVP = wexsvp(rT) * 100      !! mb --> N/m2
        rRH  = rPP/rSVP*100.0
        IF (rRH <= 100.0) THEN
          write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL), &
            (raaQout(iL,iG),iG=1,6)
        ELSE
          write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL), &
            (raaQout(iL,iG),iG=1,6),' ***** '
        END IF
      END DO
    END IF
 111 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '))
 112 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '),A7)
 113 FORMAT(A90)

! now find pressure output corresponding to HGT output from LBLRTM
    IF (kRTP == -20) THEN
      IF (raRTP_TxtInput(6) > 0) THEN
        !!! input level boundaries in km, change to mb
        DO iL = iLowestLev,kProfLayer
          raJunk(iL-iLowestLev+1) = raZout(iL)
          raJunk2(iL-iLowestLev+1) = raPout(iL)
        END DO
        CALL r_sort_linear(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo)
        write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
        raRTP_TxtInput(6) = zWoo/100.0  !! mb
      ELSEIF (raRTP_TxtInput(6) < 0) THEN
        !!! input level boundaries in mb
        raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
        write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
      END IF
    END IF

    write(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
          
    RETURN
    end SUBROUTINE DoIntegrateLevels2Layers_wrtP

!************************************************************************
! this subroutine does the klayers integrations, integrals over P
! so this is more similar to kLAYERs, CRTM
! this is newer code; it ignores surface as "lowest" level and only uses the user level info
! this is for USER DEFINED levels

! >>>>>>>> raPout comes in mb
! >>>>>>>> and goes out in N/m2

    SUBROUTINE DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon, &
    raPBnd,raAlt,iZbnd,rfracBot, &
    raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: rHSurf,rPSurf,rTSurf,raPBnd(2*kProfLayer)
    REAL :: raAlt(2*kProfLayer),rfracBot,rLat,rLon
    INTEGER :: iLowestLev,iNumGases,iaG(kMaxGas),iZbnd
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_VMRX(kProfLayer+1,kMaxGas)
! output
    REAL :: raTout(kProfLayer),raAmountOut(kProfLayer),raZout(kProfLayer+1),raPout(kProfLayer)
    REAL :: raaQout(kProfLayer,kMaxGas),raaPartPressOut(kProfLayer,kMaxGas)

! local
    INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
    REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
    REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
    REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
    REAL :: raXYZPress(kMaxlayer+1+1)
    REAL :: raXYZTemp(kMaxlayer+1+1)
    REAL :: raaXYZ_VMR(kMaxlayer+1+1,kMaxGas)
    REAL :: raXYZ_VMRwater(kMaxlayer+1+1)
    REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
    REAL :: zWoo,rConvertQ,rRH
    INTEGER :: iOffset,iWoo

    rConvertQ = 6.022141e23/1e4     ! multiply by this to change moles/m2 to molecules/cm2
    rConvertQ = kAvog/1000/1.0e4    ! multiply by this to change moles/m2 to molecules/cm2
    rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor
          
    raaQout = 0.0

    IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
      write (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
      Call DoStop
    END IF

    iNFine = 10
    iNFine = 200
    iNFine = 50

    write(kStdWarn,*) '  '
    write(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'

! >>>>>>>>>>
!! set up temporary arrays/matrices for interping onto sublevels
    iUpperLev = iZbnd
!! set levels 1 .. iNumLevs
    DO iL = 1,iUpperLev
      raXYZPress(iL) = raPX(iL)
      raXYZTemp(iL) = raTX(iL)
      DO iG = 1,iNumGases
        raaXYZ_VMR(iL,iG) = raaG_VMRX(iL,iG)
        IF (iG ==  1) raXYZ_VMRwater(iL) = raaG_VMRX(iL,iG)
      END DO
      ! print *,iNumGases,iL,raPX(iL),raTX(iL),raaG_VMRX(iL,1),raaG_VMRX(iL,7),raaG_VMRX(iL,22)
    END DO

! >>>>>>>>>>

! now start integrating
! lowest layer is special case

! http://cimss.ssec.wisc.edu/~paulv/
! http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
!     /home/sergio/KCARTA/DOC/sci_klayers.txt

 1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
     
    z = rHSurf !! need to kludge this
    zWoo = rHSurf
    iWoo = -1
          
    iOffSet = kProfLayer-(iZbnd-1)
    DO iL = iLowestLev,iLowestLev
      raPout(iL+iOffset) = log(raPBnd(iL)/raPBnd(iL+1))
      raPout(iL+iOffSet) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL+iOffset)
      raPout(iL+iOffSet) = raPout(iL+iOffSet)/100.0    !!! change N/m2 to mb
            
      z0 = z
      iJ = iL - iLowestLev + 1

      amount = 0.0
      rTsum  = 0.0
      rWgt = 0.0

      !!! >>>>>>  this is starting out at         rP = raPBnd(iL)        rT = rTSurfx  >>>>>>>>
      rP = raPBnd(iL)
      rP = rP/100.0    !!! change N/m2 to mb
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,iUpperLev,rP,rT)

      !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
      dlnp = log(raPBnd(iL)) - log(raPBnd(iL+1))
      dlnp = dlnp / (iNFine)

      !! information for current (sub)level
      rP_n = rP                             !! current sublev press
      rT_n = rT                             !! current sublev temp
      rMR_water_n = raXYZ_VMRwater(1)        !! current water MR
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP,rMR_water_n) !! current water MR

      !! information for next (sub)level
      rP_np1 = log(rP_n) - dlnp
      rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

      IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
        iWoo = +1
        zWoo   = z
      END IF

      DO iLoop = 1,iNFine

        rMolarMass_n   = kAtmMolarMass
        rMolarMass_np1 = kAtmMolarMass
        !!! if Earth, adjust molecular mass for displacement by water
        IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
        IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
        rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
        rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

        density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
        density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

        IF (kPlanet /= 3) THEN
          gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
          gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
        ELSE
          gravity = grav_earth(z,0.0,0.0,rLat,rLon)
        END IF

        Tav   = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
        SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
        dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

        rTsum = rTsum + Tav * (density_n + density_np1)/2
        rWgt  = rWgt + (density_n + density_np1)/2
                                 
        Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
        damount = Pav/Tav/kMGC * dz
        !print *,iLoop,raPBnd(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
        amount = amount + damount

        DO iG = 1,iNumGases
          DO iCnt = 1,iUpperLev
            raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
          END DO
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n  )
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1)
          raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
        END DO

       !!! now update for next iteration
        z = z + dz
        rP_n = rP_np1
        rT_n = rT_np1
        rMR_water_n = rMR_water_np1
        rP = rP_n

        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

        IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
          iWoo = +1
          zWoo   = z
        END IF

      END DO

      raTout(iL+iOffSet)      = rTsum/rWgt                 ! <<< need TAV
      raAmountOut(iL+iOffSet) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
      raZout(iL+1+iOffSet)    = z
                
      ! write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL+iOffSet),
      ! raAmountOut(iL+iOffSet),raZout(iL+1+iOffSet),raPout(iL+iOffSet),z/1000
    END DO

! displace everything  z
    iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
    raZout(iL+iOffSet)   = -(zWoo - z0)
    raZout(iL+1+iOffSet) = +(z-zWoo)
    write(kStdWarn,*) '  need to displace lowest layer heights'
    write(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
    write(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',raPBnd(iL),rPSurf,raPBnd(iL+1)
    z = z - zWoo

! go to TOA as defined by user

    DO iL = iLowestLev + 1,iZbnd-1
      z0 = z
      !! compute avg pressure of layer
      raPout(iL+iOffSet) = log(raPBnd(iL)/raPBnd(iL+1))
      raPout(iL+iOffSet) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL+iOffset)
      raPout(iL+iOffSet) = raPout(iL+iOffSet)/100.0    !!! change N/m2 to mb
            
      iJ = iL - iLowestLev + 1
      rP = raPBnd(iL)
      rP = rP/100.0    !!! change N/m2 to mb
      rT = raXYZTemp(iJ)
              
      amount = 0.0
      rTsum  = 0.0
      rWgt = 0.0

      !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
      dlnp = log(raPBnd(iL)) - log(raPBnd(iL+1))
      dlnp = dlnp / (iNFine)

      !! information for current (sub)level
      rP_n = rP                              !! current sublev press
      rT_n = rT                              !! current sublev temp
      rMR_water_n = raXYZ_VMRwater(iJ)        !! current water MR

      !! information for next (sub)level
      rP_np1 = log(rP_n) - dlnp
      rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
      CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
      CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

      DO iLoop = 1,iNFine

        rMolarMass_n   = kAtmMolarMass
        rMolarMass_np1 = kAtmMolarMass
        !!! if Earth, adjust molecular mass for displacement by water
        IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
        IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
        rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
        rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

        density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
        density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

        IF (kPlanet /= 3) THEN
          gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
          gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
        ELSE
          gravity = grav_earth(z,0.0,0.0,rLat,rLon)
        END IF

        Tav = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
        SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
        dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

        rTsum = rTsum + Tav * (density_n + density_np1)/2
        rWgt  = rWgt + (density_n + density_np1)/2
                                 
        Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
        damount = Pav/Tav/kMGC * dz
        amount = amount + damount
                  
        DO iG = 1,iNumGases
          DO iCnt = 1,iUpperLev
            raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
          END DO
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n  )
          CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1)
          raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
          ! q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
          ! q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
          ! raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + (q_n + q_np1)/2 * rConvertQ
        END DO

        !!! now update for next iteration
        z = z + dz
        rP_n = rP_np1
        rT_n = rT_np1
        rMR_water_n = rMR_water_np1
        rP = rP_n

        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1        )   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1)   !! next sublev MRw

      END DO

      raTout(iL+iOffSet)      = rTsum/rWgt             ! <<< need TAV
      raAmountOut(iL+iOffSet) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
      raZout(iL+1+iOffSet)    = z

    ! write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL),
    ! $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
    END DO

    write(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
    DO iL = 1,iOffSet-1
      raPout(iL) = 0.0
      raTout(iL) = 0.0
      raZout(iL+1) = raZout(iZbnd)
      raAmountOut(iL) = 0.0
    END DO

!! raZout is in meters
!! raaQout is in moles/m2
    DO iL = iLowestLev,iZbnd-1
      rThick = abs(raZout(iL+iOffSet) - raZout(iL+1+iOffSet)) * 100.0
      rT = raTout(iL+iOffSet)
      raPout(iL+iOffSet) = raPout(iL+iOffSet) * 100.0    !! change mb to N/m2
      DO iG = 1,iNumGases
        raaPartPressOut(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) / raAmountOut(iL+iOffSet) * raPout(iL+iOffSet)
      END DO
    END DO
          

    !IF ((iLowestLev-1) .GE. 1) THEN
    DO iL = 1,iOffSet-1
      raPout(iL) = log(raPBnd(iL)/raPBnd(iL+1))
      raPout(iL) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL)
      raTout(iL) = kStempMin
      DO iG = 1,iNumGases
        raaPartPressOut(iL,iG) = 0.0
        raaQout(iL,iG) = 0.0
        raaG_VMRX(iL,iG) = 0.0
      END DO
    END DO

    IF (iZbnd+iOffSet <= kProfLayer) THEN
      DO iL = iZbnd,kProfLayer
        raPout(iL) = 0.0
        raTout(iL) = kStempMin
        DO iG = 1,iNumGases
          raaPartPressOut(iL,iG) = 0.0
          raaQout(iL,iG) = 0.0
          raaG_VMRX(iL,iG) = 0.0
        END DO
      END DO
    END IF
            
    IF (kPlanet == 3) THEN
      !! show the RH for gas 1
      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Checking Relative Humidity (integrate over user levels)'
      write(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
      write(kStdWarn,113) '------------------------------------------------------------------------------------------'
      DO iL = iLowestLev,iZbnd-1
        z    = (raZout(iL+iOffset) + raZout(iL+1+iOffset))/2/1000
        zWoo = (raZout(iL+1+iOffset) - raZout(iL+iOffset))/1000
        rP   = raPout(iL+iOffset)            !! N/m2
        rT   = raTout(iL+iOffset)
        rPP  = raaPartPressOut(iL+iOffset,1) !! N/m2
        rSVP = wexsvp(rT) * 100      !! mb --> N/m2
        rRH  = rPP/rSVP*100.0
        IF (rRH <= 100.0) THEN
          write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet), &
          (raaQout(iL+iOffSet,iG),iG=1,6)
        ELSE
          write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet), &
          (raaQout(iL+iOffSet,iG),iG=1,6),' BAD RH!'
        END IF
      END DO
    END IF
 111 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '))
 112 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '),A7)
 113 FORMAT(A90)

 ! now find pressure output corresponding to HGT output from LBLRTM
    IF (kRTP == -20) THEN
      IF (raRTP_TxtInput(6) > 0) THEN
        !!! input level boundaries in km, change to mb
        DO iL = iLowestLev,kProfLayer
          raJunk(iL-iLowestLev+1) = raZout(iL)
          raJunk2(iL-iLowestLev+1) = raPout(iL)
        END DO
        CALL r_sort_linear(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo)
        write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
        raRTP_TxtInput(6) = zWoo/100.0  !! mb
      ELSEIF (raRTP_TxtInput(6) < 0) THEN
        !!! input level boundaries in mb
        raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
        write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
      END IF
    END IF

    write(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
          
    RETURN
    end SUBROUTINE DoIntegrateLevels2UserLayers_wrtP

!************************************************************************
! this subroutine just takes the LBLRTM mol/cm2 TAPE6 profile and does some simple calcs
    SUBROUTINE DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon, &
    PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot, &
    raPX,raTX,raLayDensityX,raaG_MRX, &
    raPressLevelsX,raTPressLevelsX,raAltitudesX, &
    raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: rHSurf,rPSurf,rTSurf,PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    REAL :: DATABASELEVHEIGHTS(kMaxLayer+1),rfracBot,rLat,rLon
    INTEGER :: iNumLays,iNumGases,iaG(kMaxGas)
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raLayDensityX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas)
! output
    REAL :: raTout(kProfLayer),raAmountOut(kProfLayer),raZout(kProfLayer+1),raPout(kProfLayer)
    REAL :: raaQout(kProfLayer,kMaxGas),raaPartPressOut(kProfLayer,kMaxGas)
    REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)
    INTEGER :: iLowestLev

! local
    INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
    REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
    REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
    REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
    REAL :: raXYZPress(kMaxlayer+1+1)
    REAL :: raXYZTemp(kMaxlayer+1+1)
    REAL :: raaXYZ_MR(kMaxlayer+1+1,kMaxGas)
    REAL :: raXYZ_MRwater(kMaxlayer+1+1)
    REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
    REAL :: zWoo,rConvertQ,rRH
    INTEGER :: iWoo

    raaQout = 0.0

    IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
      write (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
      Call DoStop
    END IF

    iLowestLev = kProfLayer - iNumLays + 1

! d onot need this, as we are using TAPE6 levels
!      rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/
!     $         (PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
! just readjust the layer info read in, and tack on additional gas profiles from refgas
    rFracBot = 1.0
    write(kSTdWarn,*) 'iLowestLev = ',iLowestLev,' rfracBot = ',rFracBot,'( if less than 1, reset to 1)'
    write(kStdWarn,*) 'PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev) = ', &
    PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev)

    DO iL = kProfLayer+1,kProfLayer+1
      raZout(iL) = raAltitudesX(iL-iLowestLev+1)*1000   !! change to meters
    END DO
    DO iL = iLowestLev,kProfLayer
      raZout(iL) = raAltitudesX(iL-iLowestLev+1)*1000   !! change to meters
      raPout(iL) = raPX(iL-iLowestLev+1) * 100.0        !! change mb to N/m2
      raTout(iL) = raTX(iL-iLowestLev+1)
      raAmountOut(iL) = raLayDensityX(iL-iLowestLev+1)
      IF (iL == iLowestLev) raAmountout(iL) = raAmountout(iL)/rFracBot
      DO iG = 1,iNumGases
        raaQout(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)
        IF (iL == iLowestLev) raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
        raaPartPressOut(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)/raLayDensityX(iL-iLowestLev+1) * raPX(iL-iLowestLev+1)*100
      END DO
    END DO

    IF ((iLowestLev-1) >= 1) THEN
      DO iL = 1,iLowestLev-1
        raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
        raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
        raTout(iL) = kStempMin
        DO iG = 1,iNumGases
          raaPartPressOut(iL,iG) = 0.0
          raaQout(iL,iG) = 0.0
          raaG_MRX(iL,iG) = 0.0
        END DO
      END DO
    END IF

          
    IF (kPlanet == 3) THEN
      !! show the RH for gas 1
      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Checking Relative Humidity'
      write(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
      write(kStdWarn,113) '------------------------------------------------------------------------------------------'
      DO iL = iLowestLev,kProfLayer
        z    = (raZout(iL) + raZout(iL+1))/2/1000
        zWoo = (raZout(iL+1) - raZout(iL))/1000
        rP   = raPout(iL)            !! N/m2
        rT   = raTout(iL)
        rPP  = raaPartPressOut(iL,1) !! N/m2
        rSVP = wexsvp(rT) * 100      !! mb --> N/m2
        rRH  = rPP/rSVP*100.0
        IF (rRH <= 100.0) THEN
          write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL), &
          (raaQout(iL,iG),iG=1,6)
        ELSE
          write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL), &
          (raaQout(iL,iG),iG=1,6),' ***** '
        END IF
      END DO
    END IF
    111 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '))
    112 FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '),A7)
    113 FORMAT(A90)

! now find pressure output corresponding to HGT output from LBLRTM
    IF (kRTP == -20) THEN
      IF (raRTP_TxtInput(6) > 0) THEN
        !!! input level boundaries in km, change to mb
        DO iL = iLowestLev,kProfLayer
          raJunk(iL-iLowestLev+1) = raZout(iL)
          raJunk2(iL-iLowestLev+1) = raPout(iL)
        END DO
        CALL r_sort_linear(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo)
        write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
        !raRTP_TxtInput(6) = zWoo/100.0  !! mb
        raRTP_TxtInput(6) = zWoo        !! keep in mb
      ELSEIF (raRTP_TxtInput(6) < 0) THEN
        !!! input level boundaries in mb
        !          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
        !          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
        raRTP_TxtInput(6) = abs(raRTP_TxtInput(6))
        write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' mb'
      END IF
    END IF

    RETURN
    end SUBROUTINE DoLBLRTMLayers2KCARTALayers

!************************************************************************
! copied from klayers/grav.f
! ALL PROTOCOL:
!    GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)

! NPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      LAT     latitude                    degrees
!    REAL      LON     longitude                   degrees
!    REAL      WINDE   wind velecity east          m/s
!    REAL      WINDN   wind velocity north         m/s
!    REAL      Z       altitude                    m

! UTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      GRAV    Earth gravity               m/s^2

! NPUT/OUTPUT PARAMETERS: none

! ETURN VALUES: none

! ARENT(S):
!    INTLEV

! OUTINES CALLED: none

! ILES ACCESSED: none

! OMMON BLOCKS: none

! ESCRIPTION:
!    Function to calculate Earth gravity (gravitation plus
!    centripetal acceleration) for points in the atmosphere.

!    It calculates surface gravity using an equation given in
!    "American Institute of Physics Handbook", 1963, page 2-102.
!    This equation is essentially a variant of the International
!    Gravity Formula with an extra term for longitude (which is
!    very minor).

!    Centripetal acceleration is tangental velocity squared over the
!    radius.

!    Gravitation is the gravitational constant times the Earth's mass
!    divided by the square of the radius.

!    Gravity at any point in the atmosphere is simply surface gravity
!    plus the change in gravitation and centripetal acceleration at
!    that point compared to the surface.

! LGORITHM REFERENCES: see DESCRIPTION

! NOWN BUGS AND LIMITATIONS: none

! OUTINE HISTORY:
!    Date     Programmer        Comments
!------------ ----------------- ----------------------------------------
! Mar  8 1995 Scott Hannon/UMBC created
! Jun 23 1995 Scott Hannon      Correct some comments
! 28 Sep 2007 Scott Hannon      Add more centripetel term comments

! ND====================================================================

!      =================================================================
    REAL FUNCTION GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)
!      =================================================================

!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
    IMPLICIT NONE

!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none

!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none

!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------

    REAL :: Z, WINDE, WINDN, LAT, LON

!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------

    REAL :: G_SUR, R, COSLT, COSLT2, SINLT2, SIN2LT, COSLON, &
    LTRAD, C_SUR, C_Z, RTOT, GRAVZ

!      Constants for (1 - b^2/a^2) with
!      a = 6.378388E+6 m = equatorial radius, and
!      b = 6.356911E+6 m = polar radius.
    REAL :: B2, ABTERM

!      Constants for normal gravity equation
!      (see "American Institute of Physics Handbook", 1963, pg 2-102)
    REAL :: G0
    REAL :: C1, C2, C3

!      Constants for pi/180, 2*pi, and Earth's rotational speed
!      of w=1/86400 rev/s
    REAL :: PI180,PI2, W

    DATA B2, ABTERM /4.041031E+13, 6.724285E-3/
    DATA G0 /9.780455/
    DATA C1,C2,C3 /5.30157E-3, -5.85E-6, 6.40E-6/
    DATA PI180,PI2,W /1.7453293E-2, 6.28318531, 1.1574074E-5/

!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none

!***********************************************************************
!***********************************************************************
!      EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************

!      Calculate longitude term
!      Add offset of 18 degrees, double it, convert to radians, and
!      take the cosine
    COSLON=COS( PI180*2.0*( LON + 18.0 ) )

!      Calculate the latitude terms
!      Convert Latitude into radians
    LTRAD=PI180*LAT
!      Calculate sine and cosine terms
    COSLT = COS(LTRAD)
    COSLT2 = COSLT**2
    SINLT2 = ( SIN(LTRAD ) )**2
    SIN2LT = ( SIN( 2.0*LTRAD ) )**2

!      Calculate the Earth's radius at this latitude
    R = SQRT( B2/( 1.0 - COSLT2*ABTERM ) )

!      Calculate total distance from Earth's center
    RTOT = R + Z

!      Calculate gravity at the Earth's surface
    G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON )

!      Calculate the centripetal term at the Earth's surface
!      Note: the centripetal acceleration due to Earth's rotation
!      is in a direction perpendicular to the Earth's rational axis,
!      but for this gravity function we are only interested in the
!      component parallel to the radial direction.
    C_SUR = COSLT2*R*(PI2*W)**2

!      Calculate the centripetal term at altitude z (with wind)
    C_Z = ( ( PI2*RTOT*COSLT*W + WINDE )**2 + (WINDN)**2 )/RTOT

!      Calculate the change in gravitation with altitude
    GRAVZ=(G_SUR + C_SUR)*(1.0 - R**2/RTOT**2)

    GRAV_EARTH = G_SUR + (C_SUR - C_Z) - GRAVZ

    RETURN
    END FUNCTION GRAV_EARTH

!************************************************************************
! function to compute vapor pressure wexp
    REAL FUNCTION WEXSVP(T)

    IMPLICIT NONE

    REAL :: T
    DOUBLE PRECISION :: TEMP,LOGP

    TEMP = DBLE(T)
    LOGP = -0.58002206D+4 / TEMP &
    + 0.13914993D+1 &
    - 0.48640239D-1 * TEMP &
    + 0.41764768D-4 * TEMP**2 &
    - 0.14452093D-7 * TEMP**3 &
    + 0.65459673D+1 * DLOG(TEMP)
    WEXSVP = SNGL( 1.0D-2 * DEXP( LOGP ) )

!      Hyland & Wexler, 1983, equation for SVP over ice
!           log Pi =  -0.56745359E+4 / T
!                    + 0.63925247E+1
!                    - 0.96778430E-2  T
!                    + 0.62215701E-6  T^2
!                    + 0.20747825E-8  T^3
!                    - 0.94840240E-12 T^4
!                    + 0.41635019E+1 log(T)
!       with T in [K] and Pi in [Pa]

    RETURN
    END FUNCTION WEXSVP

!************************************************************************
! this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
    SUBROUTINE InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta, &
    PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev,iZbnd,raPbnd)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input var
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas)
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rPmaxKCarta
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas)
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    REAL :: raPbnd(2*kProfLayer)  !! user defined boundaries
    INTEGER :: iZbnd              !! user defined boundaries
! output var
    INTEGER :: iLowestLev
    REAL :: raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_VMRX(kProfLayer+1,kMaxGas),rInJunk,rOutJunk

! local var
    INTEGER :: iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
    REAL :: raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_VMRX2(kProfLayer+1,kMaxGas)

    write(kStdWarn,*) ' ----->>>> using user defined ',iZbnd,' pressure levels for layers integration'
          
! see which user defined pressure level is equal to or greater than surface pressure
    write(kStdWarn,*) '--------------------------------------------------------------------------'
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' >>> prep : user defined levels : interp usergrid onto N levels grid >>>'

    rPmaxKCarta = -1.0
    DO iL = 1,iZbnd
      print *,'in InterpUser2UserBnd',iL,raPbnd(iL)
      IF (raPbnd(iL) > rPmaxKCarta) rPmaxKCarta = raPbnd(iL)
    END DO
    write(kStdWarn,*) 'User defined levels for integration : maxP = ',rPmaxKCarta,' Psurf = ',rPSurf,' mb'

!! assume raPbnd(1) > raPbnd(2) > raPbnd(3) > ... > raPbnd(iZbnd) and look
!! for raPbnd greater than raPSurf
    iLowestLev = iZbnd
 10 CONTINUE
    IF ((raPbnd(iLowestLev) < rPSurf) .AND. (iLowestLev > 1)) THEN
      iLowestLev = iLowestLev - 1
      GOTO 10
    ELSEIF ((raPbnd(iLowestLev) < rPSurf) .AND. (iLowestLev == 1)) THEN
      write (kStdErr,*) 'Error!!!  PSurf = ',rPSurf,' and Max UserLev Plev = ',rPmaxKCarta
      Call DoStop
    ELSEIF (raPbnd(iLowestLev) >= rPSurf) THEN
      write(kStdWarn,*)'lowest level = ',iLowestLev
    END IF

! fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
! remember this comes from LBLRTM input profile, so lowest layer will probably CLOSELY or EXACTLY match pSurf
    DO iL = 1,iZbnd
      raPX(iL) = raPbnd(iL)/100.0  ! input raPbnd press in N/m2, so change plev_kcarta to mb
    END DO
    iNumUse = iZbnd
          
! now interpolate T and MR onto this grid
    Call r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
    DO iG = 1,iNumGases
      DO iL = 1,iNumLevs
        raTemp(iL) = raaG_VMR(iL,iG)
      END DO
      Call r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
      DO iL = 1,iNumUse
        raaG_VMRX(iL,iG) = raTempX(iL)
      END DO
    END DO
          
 1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
     
! need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
    iHigh = iNumUse
 20 CONTINUE
    IF (raP(iNumLevs) > raPX(iHigh)) THEN
      print *,'wah iHigh',iNumLevs,iHigh,raP(iNumLevs),raPX(iHigh)
      iHigh = iHigh - 1
      GOTO 20
    END IF

    iHigh = iHigh + 1
    IF (iHigh <= iNumUse) THEN
      write(kSTdWarn,*) ' '
      write(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
      write(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
      write(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' to ',iNumUse
      DO iL = iHigh,iNumUse
        rInJunk = raPX(iL)
        Call r_sort_loglinear(raP,raT,iNumLevs,rInJunk,rOutJunk)
        raTX(iL) = rOutJunk

        DO iG = 1,iNumGases
          DO iX = 1,iNumLevs
            raTemp(iX) = raaG_VMR(iX,iG)
          END DO
          Call r_sort_loglinear(raP,raTemp,iNumLevs,rInJunk,rOutJunk)
          raaG_VMRX(iL,iG) = rOutJunk
        END DO
      END DO
    END IF

! now sort the gasIDs into increasing order, and do the same for the mixing ratios
    DO iL = 1,iNumGases
      raArr(iL) = iaG(iL) * 1.0
    END DO
    CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)

! assign temp arrays
    DO iG = 1,iNumGases
      iaG2(iG) = iaG(iG)
    END DO
    DO iL = 1,iNumUse
      DO iG = 1,iNumGases
        raaG_VMRX2(iL,iG) = raaG_VMRX(iL,iG)
      END DO
    END DO

! sort according to iaIndex
    DO iG = 1,iNumGases
      iaG(iG) = iaG2(iaIndex(iG))
    END DO
    DO iL = 1,iNumUse
      DO iG = 1,iNumGases
        raaG_VMRX(iL,iG) = raaG_VMRX2(iL,iaIndex(iG))
      END DO
    END DO

    write(kStdWarn,*) '--------------------------------------------------------------------------'
          
    RETURN
    end SUBROUTINE InterpUser2UserBnd

!************************************************************************
! this does the actual change for raaG_VMR(iL,iG)
! klayers/Doc/gas_units_code.txt
! klayers/Doc/toppmv.f
!              ------- 10-19 = volume mixing ratio ---------------------
!         10   parts per million volume mixing ratio (ppmv)
!              Number of gas X molecules per 1E+6 "air" molecules

!         11   parts per billion volume mixing ratio (ppbv)
!              Number of gas X molecules per 1E+9 "air" molecules

!         12   volume mixing ratio (unitless fraction)
!              Number of gas X molecules per "air" molecule

!              ------- 20-29 = mass mixing ratio -----------------------
!         20   mass mixing ratio in (g/kg)
!              Grams of gas X per kilogram of "air"

!         21   mass mixing ratio in (g/g)
!              Grams of gas X per gram of "air"

!              ------- 30-39 = partial pressure ------------------------

!         30   partial pressure in millibars (mb) Note: mb=hPa
!              Pressure of gas X as it exists in the atmosphere.
!              To clarify, this means at the corresponding profile
!              temperature and pressure level total pressure.

!         31   partial pressure in atmospheres (atm)
!              Pressure of gas X as it exists in the atmosphere

!              ------- 40-49 = water vapor humidity units --------------

!         40   relative humidity in (percent)
!              100 times actual vapor pressure divided by saturation
!              vapor pressure

!         41   relative humidity (unitless fraction)
!              Actual vapor pressure divided by saturation vapor
!              pressure

!         42   dew point temperature (Kelvin)
!              Temperaure at which vapor will start to condense out

!         43   dew point temperature (Celcius)
!              Temperaure at which vapor will start to condense out

!              Possible additional units might be
!                 inches or centimenters of water vapor
!                 grams of water vapor


! pIN is in mb

    SUBROUTINE changeLVLS_2_ppmv(iGasID,iGasUnits,iNumLevs,iG,PIN,TIN,raaG_VMR,iQuiet)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID,iGasUnits,iNumLevs,iG,iQuiet
    REAL :: PIN(2*kProfLayer),TIN(2*kProfLayer)
! input/output
    REAL :: raaG_VMR(2*kProfLayer,kMaxGas)

! local
    INTEGER :: iL,iCode,NLEV
    REAL :: MRIN(2*kProfLayer,kMaxGas),RJUNK,MDAIR
    INTEGER :: iaLocalGasID(kMaxGas)
    REAL :: MASSF(kMaxGas)
    CHARACTER(50) :: FMT
    CHARACTER(40) :: caaUnit(50)
    CHARACTER(20) :: cID

    include '../INCLUDE/TempF90/gasIDnameparam.f90'
! see cbgids.f in klayers
    DATA (iaLocalGasID(iL),iL=1,12) /01,02,03,04,05,06,07,08,09,10,11,12/
    DATA (MASSF(iL),iL=1,12) /18.015,44.010,47.9982,44.013,28.011,16.043,31.999,30.006,64.063,46.006,17.031,63.013/

    DO iL = 1,50
      caaUnit(iL) = '** Not Known **'
    END DO
    caaUnit(01) = 'layer amt molecules/cm2'
    caaUnit(10) = 'vmr in ppmv'
    caaUnit(11) = 'vmr in ppbv'
    caaUnit(20) = 'mmr in g/kg'
    caaUnit(21) = 'mmr in g/g'
    caaUnit(30) = 'part press in mb (not accepted)'
    caaUnit(31) = 'part press in atm (not accepted)'
    caaUnit(40) = 'RH % (not accepted)'
    caaUnit(41) = 'RH fraction (not accepted)'

    IF (iGasID > 12) THEN
      MASSF(iGasID) = kAtmMolarMass
    ENDIF

    MDAIR = kAtmMolarMass

    DO iL = 1,iNumLevs
      MRIN(iL,iG) = raaG_VMR(iL,iG)
    END DO

! convert to ppmv (gas unit 10)
    iCode = iGasUnits
    NLEV  = iNumLevs

    IF ((ICODE /= 10) .AND. (iQuiet < 0)) THEN
      FMT = '(I2,A,I2,A,A,A,I2,A,A,A)'
      cID = caGID(iGasID)
      write(kStdWarn,FMT) iG,' : gID=',iGasID,'(',cID,') units code ',ICODE,'(',caaUnit(ICODE),') convert to ppmv (unit code 10)'
    END IF
          
    IF ((ICODE == 10) .AND. (iQuiet < 0)) THEN
      FMT = '(I2,A,I2,A,A,A)'
      cID = caGID(iGasID)
      write(kStdWarn,FMT) iG,' : gID=',iGasID,'(',cID,') already in ppmv (unit code 10)'
                   
    ELSEIF (ICODE == 11) THEN
      !parts per billion volume mixing ratio
      !               PPMV = PPBV*1E-3
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*1E-3
     ENDDO
    
    ELSEIF (ICODE == 12) THEN
      !volume mixing ratio
      !PPMV = VMR*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*1E+6
      ENDDO
    
    ELSEIF (ICODE == 20) THEN
      !mass mixing ratio in g/kg
      !PPMV = MRgkg*((MDAIR*1E-3)/MASSF)*1E+6
      RJUNK=1E+3*MDAIR/MASSF(IGasID)
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
      ENDDO
    
    ELSEIF (ICODE == 21) THEN
      !mass mixing ratio in g/g
      !PPMV = MRgg*(MDAIR/MASSF)*1E+6
      RJUNK=1E+6*MDAIR/MASSF(IGasID)
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
      ENDDO
    
    ELSEIF (ICODE == 30) THEN
      ! partial pressure in mb  
      !PPMV = (PPmb/PIN)*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=(MRIN(IL,IG)/PIN(IL))*1E+6
     ENDDO
    
    ELSEIF (ICODE == 31) THEN
      !partial pressure in atm
      !PPMV = (PPatm*1013.25/PIN)*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=(MRIN(IL,IG)*1013.25/PIN(IL))*1E+6
      ENDDO
    
    ELSEIF (ICODE == 40 .AND. iGasID == 1) THEN
      ! relative humidy in percent
      ! note we need to change PN fron N/m2 to mb, so divide by 100
      ! PPMV = (RH%/100)*(SVP/PIN)*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*(WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+4
      ENDDO
    
    ELSEIF (ICODE == 41 .AND. iGasID == 1) THEN
      !relative humidity (fraction)
      !   PPMV = RH*(SVP/PIN)*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*(WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+6
      ENDDO
    
    ELSE
      WRITE(kStdErr,*) 'gasID,icode = ',iGasID,ICODE,' UNKNOWN COMBO to comvert to ppmv!!!'
      CALL DoStop
    ENDIF

    DO iL = 1,iNumLevs
      raaG_VMR(iL,iG) = MRIN(iL,iG)
    END DO

    RETURN
    end SUBROUTINE changeLVLS_2_ppmv

!************************************************************************
! this reads record 1.2 of a LBLRTM TAPE5
    SUBROUTINE read_record_1p2(caStr, &
    IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)

    IMPLICIT NONE
          
! input
    CHARACTER(80) :: caStr
! output
    INTEGER :: IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS

! local
    CHARACTER(1) ::  c1
    CHARACTER(2) ::  c2
    CHARACTER(3) ::  c3
    CHARACTER(4) ::  c4
    CHARACTER(5) ::  c5

! RECORD 1.2
!      IHIRAC, ILBLF4, ICNTNM, IAERSL,  IEMIT,  ISCAN, IFILTR, IPLOT, ITEST,  IATM,  IMRG,  ILAS,   IOD, IXSECT,  MPTS,  NPTS
!           5,     10,     15,     20,     25,     30,     35,    40,    45,    50, 54-55,    60,    65,     70, 72-75, 77-80
!       4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4
! this should read
! /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
!! reads HI=1 F4=1 CN=5 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0
!! reads HI=1 F4=1 CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0

    c1 = caStr(5:5)
    read(c1,'(I1)') IHIRAC
    c1 = caStr(10:10)
    read(c1,'(I1)') ILBLF4
    c1 = caStr(15:15)
    read(c1,'(I1)') ICNTNM
    c1 = caStr(20:20)
    read(c1,'(I1)') IAERSL
    c1 = caStr(25:25)
    read(c1,'(I1)') IEMIT
    c1 = caStr(30:30)
    read(c1,'(I1)') ISCAN
    c1 = caStr(35:35)
    read(c1,'(I1)') IFILTR
    c1 = caStr(40:40)
    read(c1,'(I1)') IPLOT
    c1 = caStr(45:45)
    read(c1,'(I1)') ITEST
    c1 = caStr(50:50)
    read(c1,'(I1)') IATM
    c1 = caStr(55:55)
    read(c1,'(I1)') IMRG
    c1 = caStr(60:60)
    read(c1,'(I1)') ILAS
    c1 = caStr(65:65)
    read(c1,'(I1)') IOD
    c1 = caStr(70:70)
    read(c1,'(I1)') IXSECT
    c4 = caStr(72:75)
    read(c4,'(I4)') MPTS
    c4 = caStr(77:80)
    read(c4,'(I4)') NPTS

    2 FORMAT(I2)
    4 FORMAT(I4)
     
    RETURN
    end SUBROUTINE read_record_1p2

!************************************************************************
! this reads record 1.2a of a LBLRTM TAPE5
    SUBROUTINE read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)

    IMPLICIT NONE

!!   XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL in free format
!!     XSELF  H2O self broadened continuum absorption multiplicative factor
!!     XFRGN  H2O foreign broadened continuum absorption multiplicative factor
!!     XCO2C  CO2 continuum absorption multiplicative factor
!!     XO3CN  O3 continuum absorption multiplicative factor
!!     XO2CN  O2 continuum absorption multiplicative factor
!!     XN2CN  N2 continuum absorption multiplicative factor
!!      XRAYL Rayleigh extinction multiplicative factor

! input
    CHARACTER(80) :: caStr
! output
    INTEGER :: XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

    READ (caStr,*) XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
          
    RETURN
    end SUBROUTINE read_record_1p2a
          
!************************************************************************
! this reads record 1.3 of a LBLRTM TAPE5
    SUBROUTINE read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)

    IMPLICIT NONE
          
! input
    CHARACTER(80) :: caStr
! output
    REAL :: rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,rDVOUT
    INTEGER :: ILNFLG

! local
    CHARACTER(1) ::  c1
    CHARACTER(10) ::  c10

! RECORD 1.3    (required if IHIRAC > 0; IAERSL > 0; IEMIT = 1; IATM = 1; or ILAS > 0; otherwise omit)
!             V1,     V2,   SAMPLE,   DVSET,  ALFAL0,   AVMASS,   DPTMIN,   DPTFAC,   ILNFLG,     DVOUT
!           1-10,  11-20,    21-30,   31-40,   41-50,    51-60,    61-70,    71-80,     85,      90-100
!          E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3

    10 FORMAT(E10.3)

    rSAMPLE = 4.0
    rALFAL0 = 0.04
    rAVMASS = 36
    rDPTMIN = 0.0002
    rDPTFAC = 0.001
    ILNFLG = 0
    rDVOUT = 0.0
          
    c10 = caStr(01:10)
    read(c10,10) rF1
    c10 = caStr(11:20)
    read(c10,10) rF2

    c10 = caStr(21:30)
    read(c10,10) rSAMPLE
    c10 = caStr(31:40)
    read(c10,10) rDVSET
    c10 = caStr(41:50)
    read(c10,10) rALFAL0
    c10 = caStr(51:60)
    read(c10,10) rAVMASS
    c10 = caStr(61:70)
    read(c10,10) rDPTMIN
    c10 = caStr(71:80)
    read(c10,10) rDPTFAC

!    c1 = caStr(85:85)
!    read(c1,'(I1)') ILFLG
!    c1 = caStr(91:100)
!    read(c10,10) rDVOUT

    RETURN
    end SUBROUTINE read_record_1p3

!************************************************************************
! this reads record 1.4 of a LBLRTM TAPE5
    SUBROUTINE read_record_1p4(caStr,rTSurf,raEmiss,raRefl)

    IMPLICIT NONE
          
! RECORD 1.4    (required if IEMIT = 1, or both IEMIT=2 and IOTFLG=2; otherwise omit)
!         TBOUND, SREMIS(1), SREMIS(2), SREMIS(3), SRREFL(1), SRREFL(2), SRREFL(3)
!           1-10,     11-20,     21-30,     31-40,     41-50,     51-60,     61-70
!          E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3

! input
    CHARACTER(80) :: caStr
! output
    REAL :: rTSurf,raEmiss(3),raRefl(3)

! local
    CHARACTER(10) ::  c10
    INTEGER :: iI

    10 FORMAT(E10.3)
     
    c10 = caStr(01:10)
    read(c10,10) rTSurf
    DO iI = 1,3
      raEmiss(iI) = -999.0
      c10 = caStr(iI*10+1:iI*10+10)
      read(c10,10) raEmiss(iI)
    END DO
    IF (raEmiss(1) < 0.0) THEN
      !! expecting emiss file
      raEmiss(2) = -1.0
      raEmiss(3) = -1.0
    END IF
          
    DO iI = 4,6
      raRefl(iI-3) = -999.0
      c10 = caStr(iI*10+1:iI*10+10)
      read(c10,10) raRefl(iI-3)
    END DO
    IF (raRefl(1) < 0.0) THEN
      !! expecting emiss file
      raRefl(2) = -1.0
      raRefl(3) = -1.0
    END IF

    RETURN
    end SUBROUTINE read_record_1p4
          
!************************************************************************
! this reads record 2.1 of a LBLRTM TAPE5
    SUBROUTINE read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)

    IMPLICIT NONE
          
!         IFORM, NLAYRS, NMOL, SECNTO,        ZH1,       ZH2,    ZANGLE
!            2     3-5,   6-10,  11-20,      41-48,     53-60,     66-73
!          1X,I1    I3,    I5,   F10.2,  20X, F8.2,  4X, F8.2,  5X, F8.3

!              IFORM      (0,1) column amount format flag
!                           = 0  read PAVE, WKL(M,L), WBROADL(L) in F10.4, E10.3, E10.3 formats (default)
!                        = 1  read PAVE, WKL(M,L), WBROADL(L) in E15.7 format

!              NLAYRS      number of layers (maximum of 200)

!                NMOL      value of highest molecule number used (default = 7; maximum of 35)
!                                             See Table I for molecule numbers.
!              SECNTO      user entered scale factor for the column amount for the layers defined by NLAYRS
!                                             if positive, looking up
!                                              if negative, looking down
!                                                normal value = 1.0
!                 ZH1      observer altitude
!                 ZH2      end point altitude
!              ZANGLE      zenith angle at H1 (degrees)

! input
    INTEGER :: iIOUN2
! output
    INTEGER :: iForm,iNumLevs,iNumGases
    REAL :: rSecnto,rTopHgt,rHSurf,rViewAngle

! local
    CHARACTER(80) :: caStr
    CHARACTER(10) :: c10
    CHARACTER(1) ::  c1
    CHARACTER(3) ::  c3
    CHARACTER(5) ::  c5
    CHARACTER(8) ::  c8

    REAL :: rJunk
    INTEGER :: iI

    8 FORMAT(E8.3)
    10 FORMAT(E10.3)
    111 FORMAT(A80)
     
    READ (iIOUN2,111,ERR=13,END=13) caStr

    c1 = caStr(2:2)
    read(c1,*) iFORM
    c3 = caStr(3:5)
    read(c3,*) iNumLevs

!! this is really number of layers, so need to increment by 1
!! iNumLevs = iNumLevs + 1
    c5 = caStr(6:10)
    read(c5,*) iNumGases
    c10 = caStr(11:20)
    read(c10,10) rSecnto   !! scale factor for layer ODS : positive if looking up, -ve if looking down
          
    c8 = caStr(41:48)
    read(c8,8) rTopHgt     !! observer altitude
    c8 = caStr(53:60)
    read(c8,8) rHSurf      !! end point altitude
    c8 = caStr(66:73)
    read(c8,8) rViewAngle  !! zenith angle at observer altitude (rTopHgt) so this is scanang

    IF (rHSurf > rTopHgt) THEN
      rJunk = rTopHgt
      rTopHgt = rHSurf
      rHSurf = rJunk
    END IF
          
!      print *,'record 2.1 : ',iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle
    13 CONTINUE
     
    RETURN
    end SUBROUTINE read_record_2p1
          
!************************************************************************
! this reads record 2.2 of a LBLRTM TAPE5
    SUBROUTINE read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raPavg,raTavg,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck, &
    rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta, &
    iaG,iaGasUnits,iForm,iNumGases)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'
       
! input
    INTEGER :: iNumLevs,iIOUN2,iForm,iNumGases
    REAL :: rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta
! output
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),rPSurf,rHSurf,raSumCheck(kMaxGas)
    REAL :: rTSurf,rPmin,rPmax,raAlt(2*kProfLayer),raTavg(2*kProfLayer),raPavg(2*kProfLayer)
    INTEGER :: iaG(kMaxGas),iaGasUnits(kMaxGas)
          
! local
    INTEGER :: iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2
    CHARACTER(120) :: caStrX,caStrY
    CHARACTER(30) ::  caStr30
    CHARACTER(15) ::  c15
    CHARACTER(10) ::  c10
    CHARACTER(8) ::  c8
    CHARACTER(7) ::  c7
    CHARACTER(3) ::  c3
    CHARACTER(2) ::  c2
    CHARACTER(1) ::  c1
    CHARACTER(80) :: ca80
    REAL :: rH,rP,rT,raX(kMaxGas),rJunk,rHSurfJunk

    120 FORMAT(A120)
    15 FORMAT(E15.7)
    10 FORMAT(F10.4)
    8 FORMAT(F8.3)
    7 FORMAT(F7.2)
     
    DO iG = 1,iNumGases
      iaG(iG) = iG
      iaGasUnits(iG) = 12   !!! assume hardcoded VMR
    END DO

    DO iJ = 1,iNumGases
      raSumCheck(iJ) = 0.0
    END DO
            
    DO iL = 1,iNumLevs              !!! >>>>>>>>>>>>>>>>>>>> start reading the MOLGAS profiles
      READ (iIOUN2,120) caStrY

      !! first we want to read Pave and Tave
      IF (iForm == 0) THEN
        c10 = caStrY(01:10)
        read(c10,10) rP
        c10 = caStrY(11:20)
        read(c10,10) rT
        c10 = caStrY(21:30)
        read(c10,10) rJunk
        c3 = caStrY(31:33)
        c2 = caStrY(34:35)
        ca80 = caStrY(37:120)
      ELSEIF (iForm == 1) THEN
        c15 = caStrY(01:15)
        read(c15,15) rP
        c10 = caStrY(16:25)
        read(c10,10) rT
        c10 = caStrY(26:35)
        read(c10,10) rJunk
        c3 = caStrY(36:38)
        c2 = caStrY(39:40)
        ca80 = caStrY(41:120)
      END IF

      raPavg(iL) = rP
      raTavg(iL) = rT
              
      IF (iL == 1) THEN
        !! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
        !!                      which are basically
        !!          SurfAlt      Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
        READ (ca80,*) rHSurfJunk,rPSurf,rTSurf,raAlt(iL+1),raP(iL+1),raT(iL+1)   !! p is in mb
        raAlt(iL) = rHSurfJunk
        raP(iL)   = rPSurf
        raT(iL)   = rTSurf
        IF (rPmax <= raP(iL)) rPmax = raP(iL)
        IF (rPmin > raP(iL)) rPmin = raP(iL)
        IF (rPmax <= rPSurf) rPmax = rPSurf
        IF (rPmin > rPSurf) rPmin = rPSurf
                
        IF (abs(rHSurfJunk-rHSurf) > 0.001) THEN
          write(kStdWarn,*) 'resetting rHSurf from first levl info in TAPE5 ',rHSurfJunk,rHSurf
          rHSurf = rHSurfJunk
        END IF

        IF (rPmax <= raP(iL+1)) rPmax = raP(iL+1)
        IF (rPmin > raP(iL+1)) rPmin = raP(iL+1)

        IF (raP(iL+1) > raP(iL)) THEN
          write(kStdErr,*) 'huh reading LBLRTM TAPE5  iL,raP(iL),raP(iL+1) = ',iL,raP(iL),raP(iL+1)
          CALL DoStop
        END IF
                    
        !print *,iL,raAlt(iL),raP(iL),raT(iL),rHSurfJunk,rPSurf,rTSurf

      ELSE
        !! now read next "BLANK"  and A(z) P(z) and T(z)
        !!                         which are basically
        !!                         A(z) Pz)  and T(z) for next level, ad continuum
        READ (ca80,*) raAlt(iL+1),raP(iL+1),raT(iL+1)    !! p in mb
        IF (rPmax <= raP(iL+1)) rPmax = raP(iL+1)
        IF (rPmin > raP(iL+1)) rPmin = raP(iL+1)
        !print *,iL,raAlt(iL+1),raP(iL+1),raT(iL+1)
      END IF
              
      !! now read MixRatio(z)
      IF (iNumGases <= 7) THEN
        READ (iIOUN2,120) caStrY
        read (caStrY,*) (raX(iJ),iJ=1,iNumGases)
      ELSEIF (iNumGases > 7) THEN
        iGasCntLoop = 0
 100    CONTINUE
        iX1 = iGasCntLoop*7 + 1
        iX2 = iX1 + 7
        iX2 = min(iX2,iNumGases)
        READ (iIOUN2,120) caStrY
        read (caStrY,*) (raX(iJ),iJ=iX1,iX2)
        IF (iX2 < iNumGases) THEN
          iGasCntLoop = iGasCntLoop + 1
          GOTO 100
        END IF
      END IF

      DO iJ = 1,iNumGases
       raaG_MR(iL+1,iJ) = raX(iJ)
       raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
      END DO
      IF (iL == 1) THEN
        DO iJ = 1,iNumGases
          raaG_MR(iL,iJ) = raX(iJ)
          raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
        END DO
      END IF
              
      IF (iL == 1) THEN
        IF ((rHSurf > rHmaxKCarta) .OR. (rHSurf < rHminKCarta)) THEN
          write(kStdErr,*) 'need rHmaxKCarta >= rHSurf >= rHminKCarta but have'
          write(kStdErr,*) '(rHmaxKCarta,rHSurf,rHminKCarta) = ',rHmaxKCarta,rHSurf,rHminKCarta
          CALL DoStop
        END IF
        IF ((rPSurf > rPmaxKCarta) .OR. (rPSurf < rPminKCarta)) THEN
          write(kStdErr,*) 'need rPmaxKCarta >= rPSurf >= rPminKCarta but have'
          write(kStdErr,*) '(rPmaxKCarta,rPSurf,rPminKCarta) = ',rPmaxKCarta,rPSurf,rPminKCarta
          CALL DoStop
        END IF
        IF ((rTSurf > kStempMax) .OR. (rTSurf < kStempMin)) THEN
          write(kStdErr,*) 'need kStempMax >= rTSurf >= kStempMin but have'
          write(kStdErr,*) '(kStempMax,rTSurf,kStempMin) = ',kStempMax,rTSurf,kStempMin
          CALL DoStop
        END IF
      END IF
              
    END DO       !!!! <<<<<<<<<<<<<<<<<<<< done reading the MOLGAS profiles, loop over levs

    iNumLevs = iNumLevs+1   !!! since we added on surface level
          
    kLBLRTM_toa = raP(iNumLevs)
    write(kStdWarn,*) 'kLBLRTM_toa = ',kLBLRTM_toa,' mb (when used with flux calcs)'
    write(kStdWarn,*) ' '
          
    RETURN
    end SUBROUTINE read_record_2p1p1

!************************************************************************
! read the xsect
! this reads record 2.2 of LBLRTM TAPE5
    SUBROUTINE read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'
       
! input
    INTEGER :: iIOUN2,iNumLevs,iNumGases,iNXsec
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer)
! output
    REAL :: raaG_MR(2*kProfLayer,kMaxGas),raSumCheck(kMaxGas)
    INTEGER :: iaG(kMaxGas),iaGasUnits(kMaxGas)
          
! local
    INTEGER :: iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2,iXSBIN,iFormX,iXmol,iNumLevsXsec
    CHARACTER(120) :: caStr,caStrX,caStrY
    CHARACTER(30) :: caStr30
    CHARACTER(15) ::  c15
    CHARACTER(10) ::  c10
    CHARACTER(8) ::  c8
    CHARACTER(7) ::  c7
    CHARACTER(3) ::  c3
    CHARACTER(2) ::  c2
    CHARACTER(1) ::  c1
    CHARACTER(80) :: ca80
    REAL :: raPX(2*kProfLayer),raTX(2*kProfLayer),raAltX(2*kProfLayer),raX(kMaxGas),rP,rT,rJunk
    REAL :: rTSurfX,rHSUrfX,rPSurfX

 111 FORMAT(A80)
 120 FORMAT(A120)
 15  FORMAT(E15.7)
 10  FORMAT(F10.4)
 8   FORMAT(F8.3)
 7   FORMAT(F7.2)

!! now see if there are xsec gases
    READ (iIOUN2,111,ERR=13,END=13) caStr
    READ(caStr,*) iNXsec,iXSBIN
    IF (iNXsec > 0) THEN    !!!! >>>>>>>>>>>>>>> start reading XSCGAS profiles

      DO iG = iNumGases+1,iNumGases+iNXsec
        iaG(iG) = iG
        iaGasUnits(iG) = 12   !!! assume hardcoded VMR
      END DO

      DO iJ = iNumGases+1,iNumGases+iNXsec
        raSumCheck(iJ) = 0.0
      END DO

      READ (iIOUN2,111,ERR=13,END=13) caStr     !!!! xsec names
      CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
            
      READ (iIOUN2,*) iFormX,iNumLevsXsec,iXmol
      IF (iNumLevsXsec > 2*kProfLayer) THEN
        write(kStdErr,*) 'iNumLevsXsec > 2*kProfLayer',iNumLevsXsec,kProfLayer
        CALL DoStop
      END IF
      IF (iNumLevsXsec /= (iNumLevs-1)) THEN
        write(kStdErr,*) 'iNumLevsXsec /= iNumLevs',iNumLevsXsec,iNumLevs-1
        CALL DoStop
      END IF

      !! now ready to read in the profiles, same as in molgas above!
      DO iL = 1,iNumLevsXsec
        READ (iIOUN2,111) caStrY

        !! first we want to read Pave and Tave
        IF (iFormX == 0) THEN
          c10 = caStrY(01:10)
          read(c10,10) rP
          c10 = caStrY(11:20)
          read(c10,10) rT
          c10 = caStrY(21:30)
          read(c10,10) rJunk
          c3 = caStrY(31:33)
          c2 = caStrY(34:35)
          ca80 = caStrY(37:120)
        ELSEIF (iFormX == 1) THEN
          c15 = caStrY(01:15)
          read(c15,15) rP
          c10 = caStrY(16:25)
          read(c10,10) rT
          c10 = caStrY(26:35)
          read(c10,10) rJunk
          c3 = caStrY(36:38)
          c2 = caStrY(39:40)
          ca80 = caStrY(41:120)
        END IF

        IF (iL == 1) THEN
          !! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
          !!                      which are basically
          !!          SurfAlt      Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
          READ (ca80,*) rHSurfX,rPSurfX,rTSurfX,raAltX(iL+1),raPX(iL+1),raTX(iL+1)   !! p is in mb
        ELSE
          READ (ca80,*) raAltX(iL+1),raPX(iL+1),raTX(iL+1)
        END IF

        IF (iNXsec <= 7) THEN
          READ (iIOUN2,120) caStrY
          read (caStrY,*) (raX(iJ),iJ=1,iNXsec)
        ELSEIF (iNXsec > 7) THEN
          iGasCntLoop = 0
 100      CONTINUE
          iX1 = iGasCntLoop*7 + 1
          iX2 = iX1 + 7
          iX2 = min(iX2,iNXsec)
          READ (iIOUN2,120) caStrY
          read (caStrY,*) (raX(iJ),iJ=iX1,iX2)
          IF (iX2 < iNXsec) THEN
            iGasCntLoop = iGasCntLoop + 1
            GOTO 100
          END IF
        END IF

        IF (iL == 1) THEN
          DO iJ = 1,iNXsec
            raaG_MR(iL,iJ+iNumGases) = raX(iJ)
            raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
          END DO
        END IF
                   
        DO iJ = 1,iNXsec
          raaG_MR(iL+1,iJ+iNumGases) = raX(iJ)
          raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
        END DO

      END DO   !!! loop over levels
    END IF     !!! end reading xsc gas profiles

    iNumGases = iNumGases + iNXsec

 13 CONTINUE
     
    RETURN
    end SUBROUTINE read_record_2p2
!************************************************************************
! this reads record 3.1 of LBLRTM TAPE5
    SUBROUTINE read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
           
! input
    INTEGER :: iIOUN2
! output
    INTEGER :: iNumGases,iBmax
    REAL :: rTopHgt,rSatZen,rHSurf

! local
    INTEGER :: iJunk,iI,iY
    CHARACTER(80) :: caStr
    CHARACTER(1) :: c1
    CHARACTER(2) :: c2
    CHARACTER(5) :: c5
    CHARACTER(10) :: c10
    REAL :: rCO2MIX,rJUNK

! 3.1
!      MODEL,  ITYPE, IBMAX,  NOZERO,  NOPRNT,  NMOL, IPUNCH, IFXTYP,   MUNITS,    RE, HSPACE,  VBAR, CO2MX
!          5,     10,    15,      20,      25,    30,     35,  36-37,    39-40, 41-50,  51-60, 61-70, 71-80
!         I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3, F10.3
    READ(iIOUN2,111) caStr
          
    c5 = caStr(1:5)
    read(c5,*) iJunk
    IF (iJunk /= 0) THEN
      write(kStdErr,*) 'in record 3.1 LBLRTM specifies one of its own internal profile OOPS'
      CALL dostop
    END IF

    c5 = caStr(6:10)
    read(c5,*) iJunk
    IF (iJunk /= 3) THEN
      write(kStdWarn,*) 'in record 3.1 LBLRTM specifies iType = ',iJunk,' OOPS only want 3 (gnd to space), resetting'
      iJunk = 3
    END IF

    c5 = caStr(11:15)
    read(c5,*) iBmax
    IF (iBmax == 0) THEN
      write(kStdWarn,*) 'in record 3.1 LBLRTM will use its own boundaries'
    ELSE
      write(kStdWarn,*) 'in record 3.1 LBLRTM will read in user specified boundaries'
    END IF

    c5 = caStr(16:20)
    read(c5,*) iJunk
    c5 = caStr(21:25)
    read(c5,*) iJunk
    c5 = caStr(26:30)
    read(c5,*) iNumGases
    c5 = caStr(31:35)
    read(c5,*) iJunk

    rCO2MIX = 330.0
    c10 = caStr(71:80)
    iY = -1
    iI = 1
    DO iI = 1,10
      IF (c10(iI:iI) /= ' ') iY = +1
    END DO
    IF  (iY > 0) read(c10,*) rCO2mix

! 3.2
!         H1,    H2,   ANGLE,   RANGE,   BETA,   LEN,     HOBS
!       1-10, 11-20,   21-30,   31-40,  41-50, 51-55,    61-70
!      F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3
          
    READ(iIOUN2,111) caStr
    c10 = caStr(1:10)
    read (c10,*) rHSurf
    c10 = caStr(11:20)
    read (c10,*) rTopHgt
    c10 = caStr(21:30)
    read (c10,*) rSatZen

    IF (rHSurf > rTopHgt) THEN
      rJunk = rTopHgt
      rTopHgt = rHSurf
      rHSurf = rJunk
    END IF

    111 FORMAT(A80)
     
    RETURN
    end SUBROUTINE read_record_3p1_and_3p2
!************************************************************************
! this reads in record 3.3b of LBLRTM TAPE5
! altitudes of LBLRTM layer boundaries
    SUBROUTINE read_record_3p3b(iIOUN2,iBmax,raZbnd)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iBmax,iIOUN2
! output
    REAL :: raZbnd(2*kProfLayer)

! local
    INTEGER :: iI,iCnt,i1,i2,iMax

    iCnt = 1
 20 CONTINUE
    i1 = (iCnt-1)*8 + 1
    i2 = iCnt*8
    IF (i2 > iBmax) i2 = iBmax
    READ(iIOUN2,*) (raZbnd(iI),iI=i1,i2)
    IF (i2 < iBmax) THEN
      iCnt = iCnt + 1
      GOTO 20
    END IF
     
    RETURN
    end SUBROUTINE read_record_3p3b

!************************************************************************
! this reads in record 3.4 of LBLRTM TAPE5
! user defined profile
    SUBROUTINE read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf, &
    rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, &
    raAlt,raP,raT,raSumCheck, &
    iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon, &
    iZbnd,raZbnd,raPbnd)
         
    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iNumGases,iIOUN2
! output
    INTEGER :: iNumLevs
    REAL :: raSumCheck(kMaxGas),rPSurf,rHSurf,rTSurf
    REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
    INTEGER :: iaG(kMaxGas),iaGasUnits(kMaxGas)
    REAL :: rPmin,rPmax,rYear,rLat,rLon
    REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
    REAL :: raZbnd(2*kProfLayer),raPbnd(2*kProfLayer)
    INTEGER :: iZbnd   !!! are we using default 101 levels, or LBLRTM defined?
         
! local var
    CHARACTER(5) ::  c5
    INTEGER :: iErr,iErrIO,iL,iJ,iG,iMid,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
    REAL :: raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH
    CHARACTER(80) :: caStr,caStrX,caStrY
    CHARACTER(30) :: caStr30
    CHARACTER(1) ::  c1

    raSumCheck = 0.0

    rPmin = +1.0e6
    rPmax = -1.0e+6
    iNXsec = -1

! record 3.4
    111 FORMAT(A80)
    read(iIOUN2,111) caStr
    c5 = caStr(1:5)
    read(c5,*) iNumLevs

! record 3.5
! see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
!       JCHAR = 1-6           - default to value for specified model atmosphere
!              = " ",A         - volume mixing ratio (ppmv)
!              = B             - number density (cm-3)
!              = C             - mass mixing ratio (gm/kg)
!              = D             - mass density (gm m-3)
!              = E             - partial pressure (mb)
!              = F             - dew point temp (K) *H2O only*
!              = G             - dew point temp (C) *H2O only*
!              = H             - relative humidity (percent) *H2O only*
!              = I             - available for user definition
    DO iG = 1,iNumGases
      iaG(iG) = iG
      iaGasUnits(iG) = 10   !!! assume hardcoded ppmv
    END DO
            
    DO iL = 1,iNumLevs
      ! READ (iIOUN2,*) rH,rP,rT,caStrX
      READ (iIOUN2,111) caStrY
      READ(caStrY,*) rH,rP,rT
      IF (iL == 1) THEN
        rPSurf = rP
        rHSurf = rH
        IF (abs(rTSurf-rT) >= 0.001) THEN
          !! hmm looks like info in first item of Record 3.5 is inconsistent with info in Record 1.4
          write(kStdWarn,*) 'rT first item of Record 3.5 is inconsistent with TBound info in Record 1.4'
          write(kSTdWarn,*) rTSurf,rT
          write(kStdWarn,*) 'artificially set rPSurf to be a little higher than lowest "p" entry'
          rPSurf = rPSurf-0.125
        END IF
      ! TSurf = rT
      END IF
      caStrX = caStrY(31:80)
      caStr30 = caStrX(11:40)

      IF (iL == 1) THEN
        DO iJ=1,iNumGases
          c1 = caStr30(iJ:iJ)
          IF (c1 == 'A') iaGasUnits(iJ) = 10
          IF (c1 == ' ') iaGasUnits(iJ) = 10
          IF (c1 == 'B') iaGasUnits(iJ) = -1
          IF (c1 == 'C') iaGasUnits(iJ) = 20
          IF (c1 == 'D') iaGasUnits(iJ) = -1
          IF (c1 == 'E') iaGasUnits(iJ) = -1
          IF (c1 == 'F') iaGasUnits(iJ) = 42
          IF (c1 == 'G') iaGasUnits(iJ) = 43
          IF (c1 == 'H') iaGasUnits(iJ) = 40
          IF (c1 == 'I') iaGasUnits(iJ) = -1
                      
          IF (iaGasUnits(iJ) < 0) THEN
            write(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
            CALL DOStop
          ELSE
            write(kStdWarn,*) 'LBLRTM gas number ',iJ,' "units" ',c1,' = kCARTA levels units of ',iaGasUnits(iJ)
          END IF
        END DO
      END IF

      !raP(iL) = rP * 100.0  !! change from mb to N/m2
      raP(iL) = rP          !! keep in mb
      raT(iL) = rT
      raAlt(iL) = rH
      IF (rPmax <= raP(iL)) rPmax = raP(iL)
      IF (rPmin > raP(iL)) rPmin = raP(iL)
      READ (iIOUN2,*) (raX(iG),iG=1,iNumGases)
      DO iJ = 1,iNumGases
        raaG_MR(iL,iJ) = raX(iJ)
        raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
      END DO

      IF (iL == 1) THEN
        !! rPSurf already set a few lines above, and need to be careful since we want bdry levels set
        !! so do NOT play with it here
        !          rPSurf = raP(1)/100 !! because we redo this below
        !          rPSurf = raP(1)     !! keep in mb
        IF ((rHSurf > rHmaxKCarta) .OR. (rHSurf < rHminKCarta)) THEN
          write(kStdErr,*) 'need rHmaxKCarta >= rHSurf >= rHminKCarta but have'
          write(kStdErr,*) '(rHmaxKCarta,rHSurf,rHminKCarta) = ',rHmaxKCarta,rHSurf,rHminKCarta
          CALL DoStop
        END IF
        IF ((rPSurf > rPmaxKCarta) .OR. (rPSurf < rPminKCarta)) THEN
          write(kStdErr,*) 'need rPmaxKCarta >= rPSurf >= rPminKCarta but have'
          write(kStdErr,*) '(rPmaxKCarta,rPSurf,rPminKCarta) = ',rPmaxKCarta,rPSurf,rPminKCarta
          CALL DoStop
        END IF
        IF ((rTSurf > kStempMax) .OR. (rTSurf < kStempMin)) THEN
          write(kStdErr,*) 'need kStempMax >= rTSurf >= kStempMin but have'
          write(kStdErr,*) '(kStempMax,rTSurf,kStempMin) = ',kStempMax,rTSurf,kStempMin
          CALL DoStop
        END IF
      END IF
    END DO

!! now see if there are xsec gases
    READ (iIOUN2,111,ERR=13,END=13) caStr
    READ(caStr,*) iNXsec
    IF (iNXsec > 0) THEN
      READ (iIOUN2,111,ERR=13,END=13) caStr     !!!! xsec names
      CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
      READ (iIOUN2,*) iNumLevsXsec
      IF (iNumLevsXsec > 2*kProfLayer) THEN
        write(kStdErr,*) 'iNumLevsXsec > 2*kProfLayer',iNumLevsXsec,kProfLayer
        CALL DoStop
      END IF
      IF (iNumLevsXsec /= iNumLevs) THEN
        write(kStdErr,*) 'iNumLevsXsec /= iNumLevs',iNumLevsXsec,iNumLevs
        CALL DoStop
      END IF
      DO iL = 1,iNumLevsXsec
        READ (iIOUN2,111,ERR=13,END=13) caStrY
        caStr30 = caStrY(16:45)
        IF (iL == 1) THEN
          DO iJ=1,iNXsec
            c1 = caStr30(iJ:iJ)
            IF (c1 == 'A') iaGasUnits(iJ+iNumGases) = 10
            IF (c1 == ' ') iaGasUnits(iJ+iNumGases) = 10
            IF (c1 == 'B') iaGasUnits(iJ+iNumGases) = -1
            IF (c1 == 'C') iaGasUnits(iJ+iNumGases) = 20
            IF (c1 == 'D') iaGasUnits(iJ+iNumGases) = -1
            IF (c1 == 'E') iaGasUnits(iJ+iNumGases) = -1
            IF (c1 == 'F') iaGasUnits(iJ+iNumGases) = 42
            IF (c1 == 'G') iaGasUnits(iJ+iNumGases) = 43
            IF (c1 == 'H') iaGasUnits(iJ+iNumGases) = 40
            IF (c1 == 'I') iaGasUnits(iJ+iNumGases) = -1
                                  
            IF (iaGasUnits(iJ+iNumGases) < 0) THEN
              write(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
              CALL DOStop
            ELSE
              write(kStdWarn,*) 'LBLRTM xsec number ',iJ,' is of "units" ',c1,' = kCARTA input levels units of ', &
              iaGasUnits(iJ+iNumGases)
            END IF
          END DO
        END IF

        READ (iIOUN2,*) (raX(iG),iG=1,iNXsec)
        DO iJ = 1,iNXsec
          raaG_MR(iL,iJ+iNumGases) = raX(iJ)
          raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
        END DO
      END DO
      iNumGases = iNumGases + iNXsec
    END IF

    IF (iZbnd > 0) THEN
      ! need to convert the user heights to user defined pressure levels
      call spl(raAlt,raP,iNumLevs,raZbnd,raPbnd,iZbnd)
      DO iJ = 1,iZbnd
        raPbnd(iJ) = raPbnd(iJ) * 100.0  !! change from mb to N/m2, as Psurf is finally in N/m2
      END DO
    END IF
          
 13 CONTINUE
    RETURN
    end SUBROUTINE read_record_3p4_and_3p5_and_3p7
!************************************************************************
! if TAPE5 has zeros everywhere for a certain gas, this reads in a US Std profile and stuffs it in
    SUBROUTINE substitute_tape5_profile_for_climatology(iG,iNumLevs,raP,raaG_MR)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iG             ! which gasID to put in profile
    INTEGER :: iNumLevs       ! how many levs in profile
    REAL :: raP(2*kProfLayer) ! pressure levels
! input/output
    REAL :: raaG_MR(2*kProfLayer,kMaxGas)

! local vars
    REAL :: raX(2*kProfLayer)
    INTEGER :: iL

    CALL ReadRefProf_Levels2(iG,raP,iNumLevs,raX)

    DO iL = 1,iNumLevs
      raaG_MR(iL,iG) = raX(iL)
    END DO

    RETURN
    end SUBROUTINE substitute_tape5_profile_for_climatology

!************************************************************************
! this writes the LBLRTM input so you can cut and paste into Matlab, and save RTP file
! now write the klayers stuff
    SUBROUTINE lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf, &
    raP,raT,raaG_MR,raAlt,raJunkP,raJunkT,iWriteRTP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
    INTEGER :: iWriteRTP    !! +1 for levels, +2 for layers profiles
    REAL :: rPmin,rPmax,rPSurf,rTSurf,rHSurf,rF1,rF2
    REAL :: raP(2*kProfLayer)    !! this is always LEVELS press
    REAL :: raAlt(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas) !! these could be layers or levels
    REAL :: raJunkP(2*kProfLayer),raJunkT(2*kProfLayer)                         !! these change as iWriteRTP changes
!! from 1 (layers) to 21 (levels)
          
! local var
    REAL :: raX(2*kProfLayer)
    INTEGER :: iL,iG,iLeftjust_lenstr
    CHARACTER(8) :: ca8
    CHARACTER(6) :: ca6
    CHARACTER(2) :: ca2

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' '

    write(kStdWarn,*) 'h.vcmin = ',rF1,';'
    write(kStdWarn,*) 'h.vcmax = ',rF2,';'
    ! write(kStdWarn,*) 'h.pmin = ',rPmin/100.0,';      %% rPmin in N/m2 --> mb'
    ! write(kStdWarn,*) 'h.pmax = ',rPmax/100.0,';      %% rPmax in N/m2 -->  mb'
    write(kStdWarn,*) 'h.pmin = ',rPmin,';      %% in mb'
    write(kStdWarn,*) 'h.pmax = ',rPmax,';      %% im mb'
    write(kStdWarn,*) 'h.ngas = ',iNumGases,';'
    IF (iaGasUnits(1) > 1) THEN
      IF (iWriteRTP <= 1) THEN
        write(kStdErr,*) 'trying to write out LEVELS profile but found inconsistency ...'
        CALL DoStop
      END IF
      write(kStdWarn,*) 'h.ptype = 0;'
      write(kStdWarn,*) 'h.pfields = 1;'
    ELSEIF (iaGasUnits(1) == 1) THEN
      IF (iWriteRTP /= 1) THEN
        write(kStdErr,*) 'trying to write out LAYERS profile but found inconsistency ...'
        CALL DoStop
      END IF
      write(kStdWarn,*) 'h.ptype = 1;'
      write(kStdWarn,*) 'h.pfields = 1;'
    END IF
    write(kStdWarn,*) 'h.glist = [',(iaG(iG),iG=1,iNumGases),']'';'
    write(kStdWarn,*) 'h.gunit = [',(iaGasUnits(iG),iG=1,iNumGases),']'';'

    write(kStdWarn,*) 'p.nlevs = ',iNumLevs,';'
    write(kStdWarn,*) 'p.spres = ',rPsurf/100.0,';  %% in mb'
    write(kStdWarn,*) 'p.stemp = ',rTSurf,'; %%%%% if kSurfTemp < 0, o/w use raStemp from nm_radnce'
    write(kStdWarn,*) 'p.salti = ',rHSurf,'; %%%% WOWOWOWOWOW'

    write(kStdWarn,*) 'p.satzen = 0.0;'
    write(kStdWarn,*) 'p.scanang = 0.0;'
    write(kStdWarn,*) 'p.solzen = 130.0;'
    write(kStdWarn,*) 'p.upwell = 1;'
    write(kStdWarn,*) 'p.zobs = 705000.0;'
    write(kStdWarn,*) 'p.nemis = 2;'
    write(kStdWarn,*) 'p.efreq = [600 3000]'';  %% check this with LBLRTM efreq'
    write(kStdWarn,*) 'p.emis = [1.0 1.0]'';    %% and against nm_radnce EmisFile and raEmiss'
    write(kStdWarn,*) 'p.rho = [0.0 0.0]'';     %% and against nm_radnce raRho'

    write(kStdWarn,*) ' '
    write(kStdWarn,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    IF (iWriteRTP == 1) THEN
      write(kStdWarn,*) '% we have read in TAPE6 (layers averages) so have average lay pressures as well'
      ca8 = 'p.plays'
      CALL write_stringnice(ca8,raJunkP,1.00,8,iNumLevs,-1)       !! write out LEVELS T, as you have written out raT which is layers T
            
      write(kStdWarn,*) '% we have read in TAPE6 (layers averages) but also have TAPE5 (levels T info)'
      ca8 = 'p.tlevs'
      CALL write_stringnice(ca8,raJunkT,1.00,8,iNumLevs,-1)       !! write out LEVELS T, as you have written out raT which is layers T
    ELSEIF (iWriteRTP == 21) THEN
      write(kStdWarn,*) '% we have read in TAPE5 (level info) but have some lay pressures as well'
      ca8 = 'p.plays'
      CALL write_stringnice(ca8,raJunkP,1.00,8,iNumLevs,-1)

      write(kStdWarn,*) '% we have read in TAPE5 (level info) but have some lay temps as well'
      ca8 = 'p.tlays'
      CALL write_stringnice(ca8,raJunkT,1.00,8,iNumLevs,-1)
    END IF
    write(kStdWarn,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(kStdWarn,*) ' '
          
    ca8 = 'p.plevs'
    !CALL write_stringnice(ca8,raP,0.01,8,iNumLevs,-1)       !! raP in N/m2, convert to mb for rtp
    CALL write_stringnice(ca8,raP,1.00,8,iNumLevs,-1)       !! raP in mb, keep in mb for rtp
    ca8 = 'p.ptemp'
    CALL write_stringnice(ca8,raT,1.00,8,iNumLevs,+1)
    ca8 = 'p.palts'
    CALL write_stringnice(ca8,raAlt,1000.00,8,iNumLevs,+1)  !!raAlt in km, convert to m for rtp

    DO iG = 1,iNumGases
      CALL int2str(iaG(iG),ca2)
      ca6 = 'p.gas_'
      ca8 = ca6 // ca2
      DO iL = 1,iNumLevs
        raX(iL) = raaG_MR(iL,iG)
      END DO
      CALL write_stringnice(ca8,raX,1.00,8,iNumLevs,-1)
    END DO

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' '

    RETURN
    end SUBROUTINE lblrtm2rtp

!************************************************************************
! parses LBLRTM xsec names and finds gas ids
    SUBROUTINE XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,iLBLTapeType)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    CHARACTER(80) :: caStr
    INTEGER :: iNXsec,iLBLTapeType
! input/output
    INTEGER :: iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)

! local
    INTEGER :: iaX(kMaxGas),iG,i1,i2
    CHARACTER(10) :: caStr10
    CHARACTER(15) :: caStr15

    IF (iLBLTapeType == -5) THEN
      DO iG = 1,iNXsec
        i1 = 01 + (iG-1)*10
        i2 = 10 + (iG-1)*10
        caStr10 = caStr(i1:i2)
        CALL mapXsecname_to_XsecID(caStr10,iG,iaX)
      END DO

      DO iG = iNumGases+1,iNumGases+iNXsec
        iaG(iG) = iaX(iG-iNumGases)
        iaGasUnits(iG) = 12   !!! assume hardcoded VMR
      END DO

    ELSEIF (iLBLTapeType == -6) THEN
      DO iG = 1,iNXsec
        i1 = 01 + (iG-1)*15
        i2 = 10 + (iG-1)*15
        caStr15 = caStr(i1:i2)
        caStr10 = caSTR15(1:10)
        CALL mapXsecname_to_XsecID(caStr10,iG,iaX)
      END DO

      DO iG = iNumGases+1,iNumGases+iNXsec
        iaG(iG) = iaX(iG-iNumGases)
        iaGasUnits(iG) = 1   !!! assume hardcoded molecules/cm2
      END DO
    END IF

    RETURN
    end SUBROUTINE XsecNamesLBL

!************************************************************************
! this maps xsecnames to xsec IDs
    SUBROUTINE mapXsecname_to_XsecID(caStr10,iG,iaX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    CHARACTER(10) :: caStr10  !! str we need to compare
    INTEGER :: iG            !! index
! input/output
    INTEGER :: iaX(kMaxGas)  !! modify xsec ID according to strcmp

! local
    INTEGER :: iI,iFound,iChecked
    CHARACTER(10) :: caaNamesA(13)
    CHARACTER(10) :: caaNamesB(13)
    CHARACTER(10) :: caaNamesC(13)
    CHARACTER(10) :: caX
          
    IF ((caStr10(1:1) /= ' ') .AND. (caStr10(1:1) /= '%')) THEN
      caX = caStr10
    ELSE
      CALL adjustleftstr(caStr10,caX)
    END IF

    caaNamesA(1) = 'CCL3F'
    caaNamesB(1) = 'F11'
    caaNamesC(1) = 'CFC11'
    caaNamesA(2) = 'CCL2F2'
    caaNamesB(2) = 'F12'
    caaNamesC(2) = 'CFC12'
    caaNamesA(3) = 'CCLF3'
    caaNamesB(3) = 'F13'
    caaNamesC(3) = 'CFC13'
    caaNamesA(4) = 'CF4'
    caaNamesB(4) = 'F14'
    caaNamesC(4) = 'CFC14'
    caaNamesA(5) = 'CHCL2F'
    caaNamesB(5) = 'CFC21'
    caaNamesC(5) = 'F21'
    caaNamesA(6) = 'CHC2F2'
    caaNamesB(6) = 'CFC22'
    caaNamesC(6) = 'F22'
    caaNamesA(7) = 'C2CL3F3'
    caaNamesB(7) = 'CFC113'
    caaNamesC(7) = 'F113'
    caaNamesA(8) = 'C2CL2F4'
    caaNamesB(8) = 'CFC114'
    caaNamesC(8) = 'F114'
    caaNamesA(9) = 'C2CLF5'
    caaNamesB(9) = 'CFC115'
    caaNamesC(9) = 'F115'
    caaNamesA(10) = 'CCL4'
    caaNamesB(10) = 'CCL4'
    caaNamesC(10) = ' '
    caaNamesA(11) = 'CLONO2'
    caaNamesB(11) = 'CLONO2'
    caaNamesC(11) = 'CLNO3'
    caaNamesA(12) = 'N2O5'
    caaNamesB(12) = 'N2O5'
    caaNamesC(12) = ' '
    caaNamesA(13) = 'HNO4'
    caaNamesB(13) = 'HNO4'
    caaNamesC(13) = ' '

    iaX(iG) = -1
    iFound = -1
    iChecked = 1
    10 CONTINUE
    IF ((caaNamesA(iChecked) == caX) .OR. (caaNamesB(iChecked) == caX) .OR. (caaNamesC(iChecked) == caX)) THEN
      iFound = iChecked
    END IF

    IF ((iFound < 0) .AND. (iChecked < 13)) THEN
      iChecked = iChecked + 1
      GOTO 10
    END IF

    IF (iFound < 0) THEN
      write(kStdErr,*) 'Reading LBLRTM file; could not find xsec match for ',caX
      CALL DOStop
    ELSE
      iaX(iG) = 50 + iChecked
      write(kStdWarn,*) '  LBLRTM xsec gas = ',caX,' corresponds to gasID ',iaX(iG)
    END IF

    RETURN
    end SUBROUTINE mapXsecname_to_XsecID

!                           -------------------------------------------------
!                             Alias(1)           Alias(2)           Alias(3)           Alias(4)
!                            ----------         ----------         ----------         ----------
!                            CLONO2              CLNO3
!                            HNO4
!                            CHCL2F                                 CFC21              F21
!                            CCL4
!                            CCL3F               CFCL3              CFC11              F11
!                            CCL2F2              CF2CL2             CFC12              F12
!                            C2CL2F4             C2F4CL2            CFC114             F114
!                            C2CL3F3             C2F3CL3            CFC113             F113
!                            N2O5
!                            HNO3
!                            CF4                                    CFC14              F14
!                            CHCLF2              CHF2CL             CFC22              F22
!                            CCLF3                                  CFC13              F13
!                            C2CLF5                                 CFC115             F115

!************************************************************************
! this simply writes the strings nicely
    SUBROUTINE write_stringnice(ca8,raX,rScale,iNumItemsPerLine,iNumLevs,ForE)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    CHARACTER(8) :: ca8
    INTEGER :: iNumLevs,iNumItemsPerLine,ForE
    REAL :: raX(2*kProfLayer),rScale

! local var
    INTEGER :: iL,iY,iFull,iPart,i1,i2,iLen
    REAL :: raY(2*kProfLayer)

    iFull = ifloor(iNumLevs*1.0/iNumItemsPerLine)
          
    DO iL = 1,iNumLevs
      raY(iL) = raX(iL) * rScale
    END DO

    IF (iNumItemsPerLine /= 8) iNumItemsPerLine = 8
          
    write(kStdWarn,*) ca8,' = [...'
    DO iL = 1,iFull
      i1 = 1 + (iL-1)* iNumItemsPerLine
      i2 = i1 + iNumItemsPerLine - 1
      !print *,iL,i1,i2
      IF (ForE == +1) THEN
        write(kStdWarn,1108) (raY(iY),iY=i1,i2)
      ELSEIF (ForE == -1) THEN
        write(kStdWarn,1308) (raY(iY),iY=i1,i2)
      END IF
    END DO
    IF (i2 < iNumLevs) THEN
      i1 = i2 + 1
      i2 = iNumLevs
      !print *,9999,i1,i2
      iLen = i2-i1+1
      IF (iLen > iNumItemsPerLine) THEN
        write(kSTdErr,*) 'cannot have more than 8 array members per line'
        CALL DoStop
      END IF
      IF (ForE == +1) THEN
        IF (iLen == 1) THEN
          write(kStdWarn,1101) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 2) THEN
          write(kStdWarn,1102) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 2) THEN
          write(kStdWarn,1102) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 3) THEN
          write(kStdWarn,1103) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 4) THEN
          write(kStdWarn,1104) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 5) THEN
          write(kStdWarn,1105) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 6) THEN
          write(kStdWarn,1106) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 7) THEN
          write(kStdWarn,1107) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 8) THEN
          write(kStdWarn,1108) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 9) THEN
          write(kStdWarn,1109) (raY(iY),iY=i1,i2)
        ELSE
          write(kStdErr,*) 'oops so many numbers in the line????'
          Call DoStop
        END IF
      ELSEIF (ForE == -1) THEN
        IF (iLen == 1) THEN
          write(kStdWarn,1301) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 2) THEN
          write(kStdWarn,1302) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 2) THEN
          write(kStdWarn,1302) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 3) THEN
          write(kStdWarn,1303) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 4) THEN
          write(kStdWarn,1304) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 5) THEN
          write(kStdWarn,1305) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 6) THEN
          write(kStdWarn,1306) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 7) THEN
          write(kStdWarn,1307) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 8) THEN
          write(kStdWarn,1308) (raY(iY),iY=i1,i2)
        ELSEIF (iLen == 9) THEN
          write(kStdWarn,1309) (raY(iY),iY=i1,i2)
        ELSE
          write(kStdErr,*) 'oops so many numbers in the line????'
          Call DoStop
        END IF
      END IF
    END IF
    write(kStdWarn,*) ']'';'

 1101 FORMAT(1(' ',F12.5),' ...')
 1102 FORMAT(2(' ',F12.5),' ...')
 1103 FORMAT(3(' ',F12.5),' ...')
 1104 FORMAT(4(' ',F12.5),' ...')
 1105 FORMAT(5(' ',F12.5),' ...')
 1106 FORMAT(6(' ',F12.5),' ...')
 1107 FORMAT(7(' ',F12.5),' ...')
 1108 FORMAT(8(' ',F12.5),' ...')
 1109 FORMAT(9(' ',F12.5),' ...')

 1301 FORMAT(1(' ',ES12.5),' ...')
 1302 FORMAT(2(' ',ES12.5),' ...')
 1303 FORMAT(3(' ',ES12.5),' ...')
 1304 FORMAT(4(' ',ES12.5),' ...')
 1305 FORMAT(5(' ',ES12.5),' ...')
 1306 FORMAT(6(' ',ES12.5),' ...')
 1307 FORMAT(7(' ',ES12.5),' ...')
 1308 FORMAT(8(' ',ES12.5),' ...')
 1309 FORMAT(9(' ',ES12.5),' ...')
           
    RETURN
    end SUBROUTINE write_stringnice

!************************************************************************

END MODULE n_mr_common
