! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved
!
!************************************************************************
! how to declare a function that outputs an array (or a matrix) see for example
! http://www.pcc.qub.ac.uk/tec/courses/f90/stu-notes/F90_notesMIF_5.html     or
! https://stackoverflow.com/questions/26347090/how-to-declare-the-type-of-a-function-that-returns-an-array-in-fortran
! see example of making functions to output arrays
!
!    module so_func
!    INTEGER, PARAMETER :: MAX_SIZE = 5
!    TYPE MY_DATA
!        INTEGER :: SIZE
!        REAL, DIMENSION(MAX_SIZE) :: DATA
!    ENDTYPE
! contains
!
!    FUNCTION f1(A,N) RESULT(X)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    REAL :: X(N)
!    ! ....
!    X = 1.0+A
!    END FUNCTION f1
!
!    TYPE(MY_DATA) FUNCTION f2(A,N)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    ! ....
!    f2%SIZE = N
!    f2%DATA(1:N) = 1.0+A
!    END FUNCTION f2
!
!    FUNCTION f3(A,N)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    REAL :: f3(N)
!    ! ....
!    f3 = 1.0+A
!    END FUNCTION f3
!
!end module
!
!program SO_RESULT
!    use so_func
!    implicit none
!    integer, parameter :: n=5
!    REAL :: A(n), y1(n), y3(n)    
!    TYPE(MY_DATA) :: y2
!    INTEGER :: i
!
!    ! Variables
!    A =(/ (i, i=1,n) /)
!    y1 = f1(A,n)
!    y2 = f2(A,n)
!    y3 = f3(A,n)
!end program SO_RESULT

!************************************************************************
MODULE basic_common

IMPLICIT NONE

CONTAINS

!************************************************************************
! see https://stackoverflow.com/questions/11691925/fortran-equivalent-to-matlab-find-application-to-slicing-matrix-without-memory
! example usage
!     REAL*8 A(8)
!     INTEGER n, npos, pos(8)
!     n=8
!     A = (/ 19, 20, 21, 22, 23, 24, 25, 26 /)
!     ! Find the positions of vector A that is equal to 22 
!     CALL FindInVector(n,A==22,npos,pos)
!     WRITE(*,*) pos(1:npos)
!
!     ! Find the positions of vector A that contains even numbers 
!     CALL FindInVector(n,ABS(A/2.d0-INT(A/2.d0))<1.d-2,npos,pos)
!     WRITE(*,*) pos(1:npos)

     SUBROUTINE FindInVector(n,TF,npos,pos)
    ! Inlet variables
    INTEGER,INTENT(IN):: n      ! Dimension of logical vector
    LOGICAL,INTENT(IN):: TF(n)  ! Logical vector (True or False)

    ! Outlet variables
    INTEGER npos                ! number of "true" conditions
    INTEGER pos(n)              ! position of "true" conditions

    ! Internal variables
    INTEGER i                   ! counter
    INTEGER v(n)                ! vector of all positions

    pos = 0                     ! Initialize pos
    FORALL(i=1:n)   v(i) = i    ! Enumerate all positions
    npos  = COUNT(TF)           ! Count the elements of TF that are .True.
    pos(1:npos)= pack(v, TF)    ! With Pack function, verify position of true conditions

    END SUBROUTINE FindInVector

!************************************************************************
! this integer function is "floor" -- assume rX > 0
    INTEGER FUNCTION iFloor(rX)
     
    IMPLICIT NONE

    REAL :: rX
     
    iFloor = int(floor(rX))

    RETURN
    end FUNCTION iFloor
!************************************************************************
! this integer function is "ceil" -- assume rX > 0
    INTEGER FUNCTION iCeil(rX)
     
    IMPLICIT NONE

    REAL :: rX
     
    iCeil = int(ceiling(rX))
     
    RETURN
    end FUNCTION iCeil

!************************************************************************
! this function does a temperature interpolation on a fractional layer
! this uses modified Scott Hannon's method of doing a quad fit to the layer,
! layer above, layer below  of the form
!     T = a (ln P(avg))^2 + b (ln P(avg)) + c
    REAL FUNCTION InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac, &
    iTopORBot,iL)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raVTemp  = array containing the original 1.0 fraction temps
! rFrac    = frac of layer that we need
! iTopORBot= do we need top or bottom of layer (+1/-1)
! iL       = which of the mixed paths

! for a down looking instrument, we need bottom frac
! for a   up looking instrument, we need top frac
! for bottommost layer, we need top frac
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raVTemp(kMixFilRows),rFrac
    INTEGER :: iTopORBot,iL,iProfileLayers

    REAL :: rT,rP         !user spedfd pressure, temp calculated at this press
    REAL :: rPavg         !given rP,rP1, need to compute rPavg
    REAL :: rT0,rTm1,rTp1 !avg temps of 3 adjacent layers
    REAL :: rP0,rPm1,rPp1 !avg pressures of 3 adjacent layers
    REAL :: rA,rB,rC      !need to find eqn of quadratic
    REAL :: rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr  !temporary variables
    INTEGER :: i0,im1,ip1,iW
    INTEGER :: iLowest

    iLowest = kProfLayer - iProfileLayers + 1

    IF (abs(rFrac-1.00) <= delta) THEN
        rT = raVTemp(iL)       !use the original temp .. no need to intrp
    ! thse next three lines are to debug the function, for iTopBot = +1
        rP = raPressLevels(MP2Lay(iL))
        rPp1 = raPressLevels(MP2Lay(iL)+1)
        rPavg=(rP-rPp1)/alog(rP/rPp1)

    ELSE   !oh boy .. have to interp!!!!!!!!

        iW = iCeil(iL*1.0/(kProfLayer*1.0))    !which set of mxd paths this is
        i0 = MP2Lay(iL) !lower pressure level .. rP is within this press layer
        ip1 = i0+1      !upper pressure leve1 .. this is one press layer above
        im1 = i0-1      !                     .. this is one press layer below

    ! have to recompute what the user specified pressure was!!
        IF (iTopORBot == 1) THEN          !top frac of layer
        ! ressure specified by user
            rP = raPressLevels(ip1)+rFrac*(raPressLevels(i0)-raPressLevels(ip1))
        ELSE                                !bot frac of layer
        ! ressure specified by user
            rP = -rFrac*(raPressLevels(i0)-raPressLevels(ip1))+raPressLevels(i0)
        END IF

    ! compute the average pressure of the fractional layer
        IF (iTopOrBot == 1) THEN
            IF (abs(rP-raPressLevels(ip1)) >= delta) THEN
                rPavg = (rP-raPressLevels(ip1))/alog(rP/raPressLevels(ip1))
            ELSE
                rPavg = rP
            END IF
        ELSE
            IF (abs(rP-raPressLevels(i0)) >= delta) THEN
                rPavg = (raPressLevels(i0)-rP)/alog(raPressLevels(i0)/rP)
            ELSE
                rPavg = rP
            END IF
        END IF

        IF ((i0 <= (kProfLayer-1)) .AND. (i0 >= (iLowest+1)))  THEN
        ! can safely look at layer i0, and layer above/below it
        ! avg press of layer i0+1
            rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/ &
            alog(raPressLevels(ip1)/raPressLevels(ip1+1))
        ! avg press of layer i0
            rP0 = (raPressLevels(i0)-raPressLevels(ip1))/ &
            alog(raPressLevels(i0)/raPressLevels(ip1))
        ! avg press of layer i0-1
            rPm1 = (raPressLevels(im1)-raPressLevels(i0))/ &
            alog(raPressLevels(im1)/raPressLevels(i0))
        ! temperatures of these levels from raVTemp
            rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
            rT0 = raVTemp(i0+(iW-1)*kProfLayer)
            rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
        ELSE IF (i0 == kProfLayer) THEN
        ! first redefine i0,ip1,im1
            i0 = kProfLayer-1
            ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
            im1 = i0-1    !                     .. this is one press layer below
        ! can now safely look at layer i0, and layer above/below it
        ! avg press of layer i0+1
            rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/ &
            alog(raPressLevels(ip1)/raPressLevels(ip1+1))
        ! avg press of layer i0
            rP0 = (raPressLevels(i0)-raPressLevels(ip1))/ &
            alog(raPressLevels(i0)/raPressLevels(ip1))
        ! avg press of layer i0-1
            rPm1 = (raPressLevels(im1)-raPressLevels(i0))/ &
            alog(raPressLevels(im1)/raPressLevels(i0))
        ! temperatures of these levels from raVTemp
            rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
            rT0 = raVTemp(i0+(iW-1)*kProfLayer)
            rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
                    
        ELSE IF (i0 == iLowest) THEN
        ! first redefine i0,ip1,im1
            i0 = iLowest+1
            ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
            im1 = i0-1    !                     .. this is one press layer below
        ! can now safely look at layer i0, and layer above/below it
        ! avg press of layer i0+1
            rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/ &
            alog(raPressLevels(ip1)/raPressLevels(ip1+1))
        ! avg press of layer i0
            rP0 = (raPressLevels(i0)-raPressLevels(ip1))/ &
            alog(raPressLevels(i0)/raPressLevels(ip1))
        ! avg press of layer i0-1
            rPm1 = (raPressLevels(im1)-raPressLevels(i0))/ &
            alog(raPressLevels(im1)/raPressLevels(i0))
        ! temperatures of these levels from raVTemp
            rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
            rT0 = raVTemp(i0+(iW-1)*kProfLayer)
            rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
        END IF
              
    ! now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n))
        rP0  = alog(rP0)
        rPp1 = alog(rPp1)
        rPm1 = alog(rPm1)
    !        print *,rpp1,rp0,rPm1
    !      print *,rTp1,rT0,rTm1
              
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
        rT = rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
    END IF

    InterpTemp = rT

    RETURN
    end FUNCTION InterpTemp

!************************************************************************
! this subroutine closes all files in case of an emergency stop
! assumes the message ends with '$'
    SUBROUTINE DoSTOPMesg(caMessage)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iI,iFound
    CHARACTER     caMessage*(*)
    CHARACTER(160) :: caMessage2

    DO iI = 1,160
        caMessage2(iI:iI) = ' '
    END DO

    iI = 160
    iI = len(caMessage)
    IF (iI > 160) THEN
        write(kStdErr,*) 'length of error message is over 160 characters!'
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

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    include '../INCLUDE/TempF90/kcartaparam.f90'

! caFName   = name of file that has the profile
! iNlay     = number of layers that are read in
! raProA/T/P/PP = profile amout, temperature,pressure,partial pressure
! iErrOut   = error count (usually associated with file I/O)
    CHARACTER(160) :: caFname
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
    IF ((caLine(1:1) /= '!') .AND. (caLine(1:1) /= '%')) THEN    
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
! this function converts the Mixed Path number to a layer number
    INTEGER FUNCTION MP2Lay(iNum)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input params
    REAL :: raTemp(kMixFilRows)         !! temperature structure
    REAL :: raPress(kProfLayer+1) !! pressure levels
    INTEGER :: iaRadLayer(kProfLayer)   !! which are the radiating layers
    INTEGER :: iNumLayer

! local vars
    REAL :: raT(kMixFilRows),rX,rJunk,radTdz(kProfLayer)
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
        ! print *,iI,iL,raPress(iL),raT(iL),raThickness(iL),radTdz(iL)
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
    write(kStdWarn,*) 'old tropopause (T = Tmin) is in atmosphere layer (iaRadLayer) ',iL
    write(kStdWarn,*) '    which is pressure = ',raPress(iaRadLayer(iL))

    find_tropopause = iI
          
    RETURN
    end FUNCTION find_tropopause

!************************************************************************

! Tropopause altitude determination from temperature profile measurements of reduced vertical resolution
! Nils König, Peter Braesicke, and Thomas von Clarmann, AMT 
! Accepted: 26 Jun 2019 – Published: 29 Jul 2019
!
! The tropopause constitutes a vertical separation in the atmosphere
! that segregates the lower weather active region, viz., the
! troposphere, from an upper, steadier region, the
! stratosphere. High-altitude temperature soundings that became possible
! at the end of the 19th century showed an – at that time – unexpected
! temperature behavior, where temperatures would stagnate or even
! increase with height (see Hoinka, 1997, for a historical
! overview). Once it was established that this observation was no
! measurement error, and that above the troposphere another region of
! the atmosphere exists, namely the stratosphere, an unambiguous
! definition for the height of the boundary, the tropopause, had to be
! agreed on. The earliest comprehensive definition provided by the
! British Meteorological Office was based on either the existence of a
! temperature inversion or an abrupt transition to a temperature
! gradient below 2 K km−1. If the first two criteria were not met, a
! more general vertical temperature gradient criterion was applied: “at
! the point where the mean fall of temperature for the kilometer next
! above is 2 K or less provided that it does not exceed 2 K for any
! subsequent kilometer” (Dines, 1919, cited after Hoinka, 1997). A
! similar definition, focusing solely on the lapse rate of 2 K km−1 was
! adapted by the World Meteorological Organization (WMO) in later years
! (World Meteorological Organization, 1957). Since then additional
! definitions of the tropopause have emerged, focusing on the behavior
! of dynamical quantities (e.g., Hoerling et al., 1991) or of trace gas
! changes (e.g., Pan et al., 2004). However, the most commonly used
! method to define the position of the tropopause is still the WMO
! criterion.

! This subroutine finds the tropopause by looking for the first cold point
! modelled on Scott Hannon's code tropopause_rtp.m which looks for the
! layer within the 50-400 mb range which has the lowest temp
    INTEGER FUNCTION find_tropopauseNew(raTemp,raPress,raThickness,raLayerHeight,iaRadLayer,iNumLayer)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input params
    REAL :: raThickNess(kProfLayer),raLayerHeight(kProfLayer)
    REAL :: raTemp(kMixFilRows)         !! temperature structure
    REAL :: raPress(kProfLayer+1) !! pressure levels
    INTEGER :: iaRadLayer(kProfLayer)   !! which are the radiating layers
    INTEGER :: iNumLayer

! local vars
    REAL :: raT(kMixFilRows),rX,rJunk,pProf(kProfLayer),pN(kProfLayer),pD(kProfLayer)
    REAL :: raZalts(kProfLayer+1),raZCenterAlts(kProfLayer+1),rH
    INTEGER :: iI,iL,i400mb,i50mb,iJL,iJ

    REAL :: radTdz(kProfLayer),r2Kperkm,r2Kperkm_press,rStretch,rTotal
    INTEGER i2km,iK,iStillOK,iFound,iTempInversion

! if        iaRadLayer = 001..100, everything ok
! but if eg iaRadLayer = 201..300, everything not ok with raPress

! check that thicknesses are making sense
    IF (iaRadLayer(1) < kProfLayer) THEN
      !!! this is downlook instr
      iI = iaRadLayer(1) + 2
    ELSEIF (iaRadLayer(1) .EQ. kProfLayer) THEN
      !!! this is uplook instr
      iI = iaRadLayer(iNumLayer) - 2
    END IF

    if (raThickness(iI) < 0) then
      write(kSTdErr,*) 'in find_tropopauseNew looks like layer thickness are wrong'
      write(kStdErr,*) raThickness
      call doStop
    end if

    pN = raPress(1:kProfLayer)-raPress(2:kProfLayer+1)
    pD = log(raPress(1:kProfLayer)/raPress(2:kProfLayer+1))
    pProf = pN/pD
    
    raZalts = 0.0
    raZAlts(1:kProfLayer) = raLayerHeight(1:kProfLayer)
    raZCenterAlts(1:kProfLayer) = raZalts(1:kProfLayer) + 0.5*(raZalts(2:kProfLayer+1)-raZalts(1:kProfLayer))

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
      IF (iL .GT. 1) THEN 
        radTdz(iL) = -(raTemp(iL)-raTemp(iL-1))/((raThickness(iL))/1000.0)
        radTdz(iL) = -(raTemp(iL)-raTemp(iL-1))/((raZCenterAlts(iL)-raZCenterAlts(iL-1))/1000.0)
      ELSE
        radTdz(iL) = -(raTemp(iL+1)-raTemp(iL))/((raThickness(iL))/1000.0)
        radTdz(iL) = -(raTemp(iL+1)-raTemp(iL))/((raZCenterAlts(iL+1)-raZCenterAlts(iL))/1000.0)
      END IF
      !!!! OLD write(*,'(2(I4,1x),4(F12.5,1x))') ,iI,iL,raPress(iL),raT(iL),raThickness(iL),radTdz(iL)
!      write(*,'(2(I4,1x),4(F12.5,1x))') ,iI,iL,pProf(iL),raZCenterAlts(iL),raT(iL),radTdz(iL)
    END DO

    !! check to see if there is a surface inversion
    iFound = -1
    iI = 1
    iTempInversion = 1
    do while ((iFound < 0) .and. (iI < iNumLayer))
      iL = iaRadLayer(iI)
      if (radTdz(iL) > 0.2) then
        i2Km = 0
        rStretch = raThickness(iL)/1000.0
        !! see how many layers we have to go to stretch to 2 km above this one
 18     continue
        rStretch = rStretch  + raThickness(iL+i2Km)/1000.0
        if ( (rStretch < 2.0) .and. (iI + i2Km < iNumLayer)) then
!          print *,i2km,iL,iL+i2Km,raThickness(iL),rStretch
          i2Km = i2Km + 1
          goto 18
        end if
      
        iStillOK = 0        
        do iK = 0,i2Km
          if (radTdz(iL + iK) > 0.2) then
             iStillOK = iStillOK + 1
           end if
        end do
        IF (iStillOK .EQ. i2km+1) THEN
          iFound = +1 
          iTempInversion = iI                 
!          print *,'yay invers',iI
        else
          iI = iI + 1
        end if
      else
        iI = iI + 1
      end if
    end do

    if (iTempInversion > 1) write(kStdWarn,*) 'tempinversion till iL = ',iTempInversion,' hgt = ',raZCenterAlts(iL),' km'

!    TROPOPAUSE Determines the tropopause height The WMO definition of the
!    tropopause: "The lowest level at which the lapse rate decreases to
!    2 C/km or less, provided that the average lapse rate between this
!    level and all higher levels within 2 km does not exceed 2 C/km."
    iFound = -1
    iI = iTempInversion+1  !! make sure you avoid first few km in case there is a temp inversion above surface
    iI = iTempInversion+2  !! make sure you avoid first few km in case there is a temp inversion above surface
    do while ((iFound < 0) .and. (iI < iNumLayer))
      iL = iaRadLayer(iI)
      if (radTdz(iL) < 2.0) then
        i2Km = 0
        rStretch = raThickness(iL)/1000.0
        !! see how many layers we have to go to strecth to 2 km above this one
 19     continue
        rStretch = rStretch  + raThickness(iL+i2Km)/1000.0
        if ( (rStretch < 2.0) .and. (iI + i2Km < iNumLayer)) then
!          print *,i2km,iL,iL+i2Km,raThickness(iL),rStretch
          i2Km = i2Km + 1
          goto 19
        end if
      
        iStillOK = 0        
        do iK = 0,i2Km
          if (radTdz(iL + iK) < 2.0) then
             iStillOK = iStillOK + 1
           end if
        end do
        IF (iStillOK .EQ. i2km+1) THEN
          iFound = +1                  
          r2Kperkm = radTdz(iL)
          r2Kperkm_press = pProf(iL)
!          print *,'yay ',iI,iL,r2Kperkm_press
        else
          iI = iI + 1
        end if
      else
        iI = iI + 1
      end if
    end do

    write(kStdWarn,'(A,I3,I3)') 'new tropopause (dT/dz < 2 K/km) is in atmosphere layer iI,iaRadLayer(iI) ',iI,iL
    write(kStdWarn,'(A,F12.5,A,F12.5,A)') '    which is layer avg pressure = ',pProf(iaRadLayer(iI)), &
                                        ' mb; altitude ',raZCenteralts(iaRadLayer(iI))/1000,' km'
    find_tropopauseNew = iL
          
    RETURN
    end FUNCTION find_tropopauseNew

!************************************************************************
! this subroutine computes d(Brightness Temp)/d(Rad) for jacobian output
    SUBROUTINE Find_BT_rad(raInten,radBTdr,raFreq, &
    radBackgndThermdT,radSolardT)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
! raInten is the radiance intensity at the instrument
! raFreq are the frequencies
! radBTdr is the derivative result
    REAL :: raFreq(kMaxPts),raInten(kMaxPts),radBTdr(kMaxPtsJac)
    REAL :: radBackgndThermdT(kMaxPtsJac),radSolardT(kMaxPtsJac)
          
    INTEGER :: iFr
    REAL :: r1,r2,ra3(kMaxPts),ra4(kMaxPts)

!! need these for derivatives of Planck
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    ra3 = r1*r2 * (raFreq**4)/(raInten**2)
    ra4 = 1.0 + r1 * (raFreq**3)/raInten
    radBTdr = ra3/ra4/(alog(ra4)**2)

    IF (kThermal < 0) THEN
      radBackgndThermdT = 0.0
    ELSE
      ra3 = r1*r2 * (raFreq**4)/(radBackgndThermdT**2)
      ra4 = 1.0 + r1 * (raFreq**3)/radBackGndThermdT
      radBackgndThermdT = ra3/ra4/(alog(ra4)**2)
    END IF

    IF (kSolar < 0) THEN
      radSolardT = 0.0
    ELSE
      ra3 = r1*r2 * (raFreq**4)/(radSolardT**2)
      ra4 = 1.0 + r1 * (raFreq**3)/radSolardT
      radSolardT = ra3/ra4/(alog(ra4)**2)
    END IF

    RETURN
    end SUBROUTINE Find_BT_rad
!************************************************************************
! **********************  generic jacobian stuff ************************
!************************************************************************
! this subroutine does d/dr(tau_layer2space) for gas iG
! where r == gas amount q or temperature T at layer iM
! and  iL is the relevant layer we want tau_layer2space differentiated
! HENCE IF iL > iM, derivative == 0
! i.e. this does d(tau(l--> inf)/dr_m
    SUBROUTINE JacobTerm(iL,iM,raaLay2Sp,raTemp)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaLay2Sp is the transmission frm layer to space
! iM has the layer that we differentiate wrt to
! iL has the radiating layer number (1..kProfLayerJac)
! raTemp has the results, apart from the multiplicative constants
!   which are corrected in MinusOne
    INTEGER :: iL,iM
    REAL :: raTemp(kMaxPtsJac),raaLay2Sp(kMaxPtsJac,kProfLayerJac)

    IF (iL > iM) THEN
      raTemp = 0.0
    ELSE
      raTemp = raaLay2Sp(:,iL)
    END IF

    RETURN
    end SUBROUTINE JacobTerm

!************************************************************************
! this subroutine multiplies the array by -1.0*constant where constant
! depends on whether we are doing d/dT or d/dq
    SUBROUTINE MinusOne(raTorQ,raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raResults is the array
! raTorQ === relevant element of raaaDq or raaDt
    REAL :: raTorQ(kMaxPtsJac)
    REAL :: raResults(kMaxPtsJac)

    raResults = -raResults * raTorQ
     
    RETURN
    end SUBROUTINE MinusOne

!************************************************************************

END MODULE basic_common
