! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE basic_common

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine closes all files in case of an emergency stop
! assumes the message ends with '$'
    SUBROUTINE DoSTOPMesg(caMessage)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iI,iFound
    CHARACTER     caMessage*(*)
    CHARACTER(120) :: caMessage2

    DO iI = 1,80
        caMessage2(iI:iI) = ' '
    END DO

    iI = 80
    iI = len(caMessage)
    IF (iI > 120) THEN
        write(kStdErr,*) 'lengthh of error message is over 120 characters!'
        write(kStdErr,*) caMessage
        CALL DoStop
    END IF

    5 CONTINUE
    IF ((caMessage(iI:iI) /= '$') .AND. (iI > 1)) THEN
        iI = iI - 1
        GOTO 5
    END IF
     
    IF (iI <= 1) THEN
        write(kStdErr,*) 'caMessage needs "$" to end '
        CALL DoStop
    END IF

!      write(kStdErr,*) 'length of caMessage = ',iI
    caMessage2(1:iI-1) = caMessage(1:iI-1)

    write(kStdErr,10) caMessage2
    CALL DoStop

    10 FORMAT(A120)

    RETURN
    end SUBROUTINE DoSTOPMesg

!***********************************************************************
! this subroutine closes all files in case of an emergency stop
    SUBROUTINE DoSTOP

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    write(kStdWarn,*)'Fatal Error found : closing all units ..'

    IF ((kStdDriverOpen == 1) .AND. (kStdDriver /= 5)) THEN
        write(kStdWarn,*)'closing driver file'
        CLOSE(UNIT = kStdDriver)          !close driver file
    END IF

    IF ((kStdkCartaOpen == 1) .AND. (kStdkCarta /= 6)) THEN
        write(kStdWarn,*)'closing binary output file'
        CLOSE(UNIT = kStdkCarta)        !close file where kCARTA binary goes to
    END IF

    IF (kCompUnitOpen == 1) THEN
        write(kStdWarn,*)'closing kcomp/xsec file'
        CLOSE(UNIT = kCompUnit)         !close kCompressed file/xsec data file
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

    IF (kFlux > 0) THEN
        write(kStdWarn,*)'closing flux binary file'
        CLOSE(UNIT = kStdFlux)         !close file where flux binary goes to
    END IF

    IF (kStdPlanckOpen > 0) THEN
        write(kStdWarn,*)'closing planck binary file'
        CLOSE(UNIT = kStdPlanck)        !close file where planck binary goes to
    END IF

    IF (kProfileUnitOpen == 1) THEN
        write(kStdWarn,*)'closing profile file '
        CLOSE(UNIT = kProfileUnit)       !close profile file
    END IF

    IF (kTempUnitOpen == 1) THEN
        write(kStdWarn,*)'closing temporary param file'
        CLOSE(UNIT = kTempUnit)          !close temporary file eg comp.param
    END IF

    IF (kBloatPlanckOpen == 1) THEN
        write(kStdWarn,*)'closing bloated planck binary file'
        CLOSE(UNIT = kBloatNLTEPlanck)      !close file
        kBloatOutOpen = -1
    END IF

    IF (kBloatOutOpen == 1) THEN
        write(kStdWarn,*)'closing bloated binary file'
        CLOSE(UNIT = kBloatNLTEOut)      !close file
        kBloatOutOpen = -1
    END IF

    IF (kStdPlanckUAOpen == 1) THEN
        write(kStdWarn,*)'closing UA planck binary file'
        CLOSE(UNIT = kStdPlanckUA)      !close file
        kStdPlanckUAOpen = -1
    END IF

    IF (kNLTEOutUAOpen == 1) THEN
        write(kStdWarn,*)'closing UA binary file'
        CLOSE(UNIT = kNLTEOutUA)      !close file
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

    write(kStdErr,*) 'bad luck ... emergency exit!'
    write(kStdWarn,*) 'bad luck ... emergency exit!'

    CLOSE(UNIT = kStdErr)             !close error log
    CLOSE(UNIT = kStdWarn)            !close warning log
     
    call exit(1)                    !sad exit so return +1

    STOP

    RETURN
    end SUBROUTINE DoSTOP

!************************************************************************
! this function checks to which posn GasID is in, in iaGases
! it mimics the "ismember" function in Matlab
    INTEGER FUNCTION WhichGasPosn(iGasID,iaGases,iNumGases)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iGasID   = current gasID
! iaGases  = list of GasID's that are being used
! iJacob   = number of GasID's that are being used
    INTEGER :: iGasID,iNumGases,iaGases(kMaxGas)

    INTEGER :: iI,iFound,iAns

    iFound = -1
    iAns = -1
    iI = 1

    15 CONTINUE
    IF ((iFound < 0) .AND. (iI <= iNumGases)) THEN
    ! check to see if iGasID is in iaGases
        IF (iGasID == iaGases(iI)) THEN
            iFound = 1
            iAns = iI
        ELSE
            iI = iI+1
            GO TO 15
        END IF
    END IF
                
    IF ((iGasID == 101) .OR. (iGasID == 102)) THEN
        iAns  =  1
    END IF

    WhichGasPosn = iAns

    RETURN
    end FUNCTION WhichGasPosn

!************************************************************************
! this function checks to see if current GasID should have its d/dq saved
! if it does, the function result is WHICH gas it is in the *JACOBN wishlist
! else the function result = -1
! originally in kcartamisc.f90, first needed in s_writefile

    INTEGER FUNCTION DoGasJacob(iGasID,iaJacob,iJacob)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iGasID   = current gasID
! iaJacob  = list of GasID's whose d/dq we want to output
! iJacob   = number of GasID's whose d/dq we want to output
    INTEGER :: iGasID,iJacob,iaJacob(kMaxDQ)

    INTEGER :: iI,iFound,iAns

    iFound = -1
    iAns = -1
    iI = 1

    15 CONTINUE
    IF ((iFound < 0) .AND. (iI <= iJacob)) THEN
    ! check to see if iGasID is in iaJacob
        IF (iGasID == iaJacob(iI)) THEN
            iFound = 1
            iAns = iI
        ELSE
            iI = iI+1
            GO TO 15
        END IF
    END IF
                
    DoGasJacob = iAns

    RETURN
    end FUNCTION DoGasJacob

!************************************************************************
! this subroutine reads the profile files (for the references)
! for either the lower or the upper atm
! it flags an error if kProfLayers layers are not read in
! ProX === A=amount,T=temperature,P=pressure,PP=partial pressure
! from n_pth_mix, first needed in kcoeffSPL

    SUBROUTINE ReadRefProf(caFname,iNlayIn,raProA,raProT, &
    raProP,raProPP,iErrOut)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! caFName   = name of file that has the profile
! iNlay     = number of layers that are read in
! raProA/T/P/PP = profile amout, temperature,pressure,partial pressure
! iErrOut   = error count (usually associated with file I/O)
    CHARACTER(80) :: caFname
    INTEGER :: iNlayIn,iErrOut
    REAL :: raProA(*),raProT(*)
    REAL :: raProP(*),raProPP(*)

! local variables
    INTEGER :: iErr,iJ,iNlay,iIOUN
    CHARACTER(100) :: caLine

    iIOUN = kProfileUnit

    OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='FORMATTED', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        WRITE(kStdErr,1080) iErr, caFname
        1080 FORMAT('ERROR! number ',I5,' opening profile file:',/,A82)
        CALL DoSTOP
    ENDIF
    kProfileUnitOpen=1
           
!      Read the file (skip comment lines)
    iNlay=0
    20 READ(iIOUN,5020,END=199) caLine
    5020 FORMAT(A100)
    IF (caLine(1:1) /= '!') THEN
        iNlay=iNlay+1
        READ(caLine,*) iJ,raProP(iNlay),raProPP(iNlay), &
        raProT(iNlay),raProA(iNlay)
    ENDIF
    GOTO 20
!      Close the file
    199 CLOSE(iIOUN)
    kProfileUnitOpen=-1

! check to see that ALL layers have been read in (could be greater than 100)
    IF (iNlay /= iNlayIN) THEN
        iErrOUt=1
        WRITE(kStdErr,500) caFName,iNLay,iNLayIN
        500 FORMAT ('Profile File',/,A82,/,' has ',I4,' layers (need ',I4,')')
        CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE ReadRefProf

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    DOUBLE PRECISION FUNCTION dttorad(df,dBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    DOUBLE PRECISION :: dI
    DOUBLE PRECISION :: dF,dBT
          
! local variables
    DOUBLE PRECISION :: d1,d2,dPlanck,dexp
    INTEGER :: iInt
     
    d1 = (kPlanck1)
    d2 = (kPlanck2)

!! 10^10 = e^23.03
    dPlanck = d2*df/dBT
    IF (dPlanck > 23.03) THEN
        dPlanck = 1.0e10
    ELSE
        dPlanck = dexp(dPlanck) - 1
    END IF

    dI = d1*(df**3)/dPlanck

    dttorad = dI

    RETURN
    END FUNCTION dttorad

!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for one point
    REAL function radtot(rf,rI)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    rPlanck = alog(1.0+(r1*(rf**3))/rI)
    rBT     = r2*rf/rPlanck

    radtot = rBT

    RETURN
    end function radtot

!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for an array
    SUBROUTINE radtot_array(raFreq,raInten,raBrightTemp)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq        = array containing wavenumbers
! raInten        = intensity from forward model
! raBrightTemp   = brightness temperatures associated with raInten
    REAL :: raFreq(kMaxPts),raInten(kMaxPts),raBrightTemp(kMaxPts)

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    DO iInt=1,kMaxPts
        rPlanck = alog(1.0+(r1*(raFreq(iInt)**3))/raInten(iInt))
        raBrightTemp(iInt) = r2*raFreq(iInt)/rPlanck
    END DO

    RETURN
    end SUBROUTINE radtot_array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one array point
    SUBROUTINE ttorad_oneBT2array(raF,rBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    DO iFr = 1,kMaxPts
        rPlanck = r2*raF(iFr)/rBT
        IF (rPlanck > 87.49) THEN
            rPlanck = 1.0e38
        ELSE
            rPlanck = exp(rPlanck) - 1.0
        END IF
        raInten(iFr) = r1*(raF(iFr)**3)/rPlanck
    END DO

    RETURN
    end SUBROUTINE ttorad_oneBT2array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities for array
    SUBROUTINE ttorad_array(raF,raBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),raBT(kMaxPts)

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    DO iFr = 1,kMaxPts
        rPlanck = r2*raF(iFr)/raBT(iFr)
        IF (rPlanck > 87.49) THEN
            rPlanck = 1.0e38
        ELSE
            rPlanck = exp(rPlanck) - 1.0
        END IF
        raInten(iFr) = r1*(raF(iFr)**3)/rPlanck
    END DO

    RETURN
    end SUBROUTINE ttorad_array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    REAL function ttorad(rf,rBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
! Constants; values from NIST (CODATA98)
!   c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
!   h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
!   k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1
!   c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
!   c2 = (h*c/k) * 100;

! at small T, exp(c2 fr/T) >>> 1
!   rad --> c1 fr^3  exp(-c2 fr/T)
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
          
    rPlanck = r2*rf/rBT
    IF (rPlanck > 87.49) THEN
        rPlanck = 1.0e38
    ELSE
        rPlanck = exp(rPlanck) - 1
    END IF

    rI = r1*(rf**3)/rPlanck

    ttorad = rI

    RETURN
    end function ttorad

!************************************************************************
! this function converts the Mixed Path number to a layer number
    INTEGER FUNCTION MP2Lay(iNum)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNum === mixed path that we want to convert to a layer
! eg 110 --> 10       200 --> 100
    INTEGER :: iNum
     
    INTEGER :: iT

    iT = MOD(iNum,kProfLayer)
    IF (iT == 0) THEN
        iT = kProfLayer
    END IF
     
    MP2Lay = iT
     
    RETURN
    end FUNCTION MP2Lay
!************************************************************************
! this subroutine finds the tropopause by looking for the first cold point
! modelled on Scott Hannon's code tropopause_rtp.m which looks for the
! layer within the 50-400 mb range which has the lowest temp
    INTEGER FUNCTION find_tropopause(raTemp,raPress,iaRadLayer,iNumLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL :: raTemp(kMixFilRows)         !! temperature structure
    REAL :: raPress(kProfLayer+1) !! pressure levels
    INTEGER :: iaRadLayer(kProfLayer)   !! which are the radiating layers
    INTEGER :: iNumLayer

! local vars
    REAL :: raT(kMixFilRows),rX,rJunk
    INTEGER :: iI,iL,i400mb,i50mb,iJL,iJ

! if        iaRadLayer = 001..100, everything ok
! but if eg iaRadLayer = 201..300, everything not ok with raPress
    iI  = 1
    iL  = iaRadLayer(iI)
    iI  = iNumLayer
    iJL = iaRadLayer(iI)
    iL = max(iL,iJL)
    iJ = 0
    IF (iL > kProfLayer) THEN
        iJ = 1
        4 CONTINUE
        IF ((iJ+1)*kProfLayer >= iL) THEN
            GOTO 5
        ELSE
            iJ = iJ + 1
            GOTO 4
        END IF
    END IF

    5 CONTINUE
    iJ = iJ*kProfLayer

    DO iI = 1,kProfLayer
        raT(iI) = 0.0
    END DO

    DO iI = 1,iNumLayer
        iL = iaRadLayer(iI)
        raT(iL) = raTemp(iL)   !!note storing into raT(iL) instead of raT(iI)
    !        print *,iI,iL,raPress(iL),raT(iL)
    END DO

!! find i400mb
    rJunk = 1.0e10
    DO iI = 1,iNumLayer
        iL = iaRadLayer(iI)
        iJL = iL - iJ
        rX = abs(400.0 - raPress(iJL))
        IF (rX <= rJunk) THEN
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
        IF (rX <= rJunk) THEN
            i50mb = iL
            rJunk = rX
        END IF
    END DO
          
!! now look for bottom layer within [i400mb,i50mb] with lowest cold point
    iL = i50mb+1
    rJunk = raT(iL)
    DO iI = i50mb,i400mb,-1
        IF (raT(iI) <= rJunk) THEN
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
        IF (iaRadLayer(iI) == iL) THEN
            GOTO 10
        END IF
    END DO

    10 CONTINUE
    write(kStdWarn,*) 'this is in atmosphere layer (iaRadLayer) ',iI
    find_tropopause = iI
          
    RETURN
    end FUNCTION find_tropopause

!************************************************************************

END MODULE basic_common
