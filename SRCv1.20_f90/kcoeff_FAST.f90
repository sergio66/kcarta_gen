! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeff_FAST

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine determines the weights for temp and pressure interpolation
! duplicating the Matlab 2011 version
!************************************************************************
    SUBROUTINE xWeights(raPPart,raPTemp,pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)


    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    INTEGER :: iPlev
    include '../INCLUDE/KCARTA_databaseparam.f90'

! input
! pProf           = actual layers from kLAYERS avg pressure
! iProfileLayers  = number of layers
! iSPlineType     = if we have previously determined fast or slow interp
! raPTemp         = actual temperature profile
! raPPart         = actual water vapor partial pressure profile

! output
!     iaP1,iaP2      are the indices          for pressure    offsets
!     raP1,raP2      are the weights          for pressure    offsets
!     iaT11,iaT12    are the indices          for temperature offsets
!     raT11,raT12    are the weights          for temperature offsets
!     raJT11,raJT12  are the jacobian weights for temperature offsets
!     iaQ11,iaQ12    are the indices          for WV    offsets
!     raQ11,raQ12    are the weights          for WV    offsets
                         
    REAL ::    raPTemp(kProfLayer),pProf(kProfLayer),raPPart(kProfLayer)
    INTEGER :: iProfileLayers,iSplineType
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

    REAL :: raToffSet(kMaxTemp),raTSpan(kMaxTemp)
    REAL :: raPPoffSet(kMaxWater),raPPSpan(kMaxWater)

! local vars
    INTEGER :: iL,iLr,iGasID,iLowest,iSwap
    REAL :: rP,p1,p2,rSumWgt1,rSumWgt2,rSwap,p1LEV,p2LEV
    REAL :: rT,t11,t12,t21,t22,wsum
    REAL :: rQ,q11,q12,q21,q22,qsum
    INTEGER :: i1,i2,iFindMaxMin
    INTEGER :: ix1,ix2,ix3
    REAL ::    rx1,rx2,rx3

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName,caStr
    INTEGER :: iE,iX
    REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer)
    REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
    REAL :: xPLEV_KCARTADATABASE_AIRS(kMaxLayer+1)

    iGasID = 1           !!! this will give all the info we need
    CALL FindReferenceName(caFName,iGasID,-1)
    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T, &
    raOrig100P,raOrig100PP,iE)

    iLowest = kProfLayer - iProfileLayers + 1

    rsumWgt1 = 0.0
    rsumWgt2 = 0.0

    DO iL = 1,kMaxLayer+1
        xPLEV_KCARTADATABASE_AIRS(iL) = PLEV_KCARTADATABASE_AIRS(iL)/1013.25
    END DO

    i1 = 0
    DO iL = -5,5
        i1 = i1 + 1
        raToffSet(i1) = iL*10.0
    END DO

    raPPoffSet(1) = 0.1
    raPPoffSet(2) = 1.0
    raPPoffSet(3) = 3.3
    raPPoffSet(4) = 6.7
    raPPoffSet(5) = 10.0
    456 FORMAT(I3,6(' ',ES12.5))
     
    write(kStdWarn,*) 'Figuring out the kComp Corner Weights ....'
    caStr = &
    'Layer  P(mb)  span   wgtP1    T(K)  indx    wgtT1    water span   wgtQ1'
    write(kStdWarn,*) caStr
    caStr = &
    '              indx1                                  amt   indx1'
    write(kStdWarn,*) caStr
    caStr = &
    '-----------------------------|------------------------|------------------------'
    write(kStdWarn,*) caStr
    caStr = &
    'iL     rP       iaP1    raP1 |  rT     iaT11   raT11  |   rQ       iaQ11  raQ11'
    write(kStdWarn,*) caStr
    caStr = &
    '-----------------------------|------------------------|------------------------'
    write(kStdWarn,*) caStr

!      DO iL = 1,kProfLayer
!        print *,iL,xPLEV_KCARTADATABASE_AIRS(iL),raOrig100P(iL),xPLEV_KCARTADATABASE_AIRS(iL+1)
!      END DO
!      print *,p1LEV,rP,p2LEV,' ',p1,p2
!      call dostop
          
    DO iL = iLowest,kProfLayer
    !!!!!!!!!!!!!!!!!!!!!! pressure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rP = pProf(iL)/1013.25
    !! get pressure bounding interval [p1 p2]
    !! replaced .... raOrig100P,kMaxLayer) with
    !!          .... xPLEV_KCARTADATABASE_AIRS,kMaxLayer+1)
    !! so

    !! this was code till Dec 2015
        iaP1(iL) = iFindMaxMin(+1,rP,xPLEV_KCARTADATABASE_AIRS,kMaxLayer+1)
        iaP2(iL) = iFindMaxMin(-1,rP,xPLEV_KCARTADATABASE_AIRS,kMaxLayer+1)

    !! this is newer code, Jan 2016
        iaP1(iL) = iFindMaxMin(+1,rP,raOrig100P,kMaxLayer)
        iaP2(iL) = iFindMaxMin(-1,rP,raOrig100P,kMaxLayer)

    !! >>>>>>>>>>> now switch from LEVELS to average LAY pressure <<<<<<<<<
        IF (iaP1(iL) > kMaxLayer) iaP1(iL) = kMaxLayer
        IF (iaP2(iL) > kMaxLayer) iaP2(iL) = kMaxLayer
        p1 = raOrig100P(iaP1(iL))  !! lower avg press bound (so upper in hgt)
    !! if iSplineType ~ 2 weight should be ~ 0
    !! as should be iL=iaP2(iL) should be 1
    !! ie they are OFFSET and NOT the same layer
        p2 = raOrig100P(iaP2(iL))  !! upper avg press bound (so lower in hgt)
    !! if iSplineType ~ 2 weight should be ~ 1
    !! as should be one-> one correspondance
    !! between iL and iaP2(iL)
        p1LEV = xPLEV_KCARTADATABASE_AIRS(iaP1(iL))
        p2LEV = xPLEV_KCARTADATABASE_AIRS(iaP2(iL))
              
    !! pressure interp weight
    !! oops this originally was 1e-3. but the <p> and rp are in atm, so can range from 1013/1013 to 0.005/1013
    !! or from 1 to 1e-6
        IF (abs(p1-p2) >= 1.0e-6) THEN
            raP2(iL) = (rP - p1) / (p2 - p1)
            raP1(iL) = 1.0 - raP2(iL)
        ELSE
            raP2(iL) = 1.0
            raP1(iL) = 0.0
        END IF
              
    !      write(kStdErr,456) ,iL,p1,rp,p2,p1-p2,raP1(iL),raP2(iL)

    !! want raP2(iL) > raP1(iL)
        IF (raP1(iL) < raP2(iL)) THEN
            rSwap = raP1(iL)
            raP1(iL) = raP2(iL)
            raP2(iL) = rSwap
            iSwap = iaP1(iL)
            iaP1(iL) = iaP2(iL)
            iaP2(iL) = iSwap
        END IF
        IF (abs(iSplinetype) == +2) THEN
            CALL Check_xWeights_iSplinetype(iL,iaP1,iaP2,raP1,raP2)
        END IF

        rSumWgt1 = rSumWgt1 + raP1(iL)
        rSumWgt2 = rSumWgt2 + raP2(iL)

    !!!!!!!!!!!!!!!!!!!!!! temperature !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        rT = raPTemp(iL)
    !! get temperature tabulation bounding interval [t11, t12] at p1
        DO i2 = 1,kMaxTemp
            raTSpan(i2) = raTOffSet(i2) + raOrig100T(iaP1(iL))
        END DO
        iaT11(iL) = iFindMaxMin(+1,rT,raTSpan,kMaxTemp)
        iaT12(iL) = iFindMaxMin(-1,rT,raTSpan,kMaxTemp)
    !! need iaT11 and iaT12 to be separate, for the jacobian
        IF (iaT11(iL) == iaT12(iL)) THEN
            IF (iaT11(iL) == kMaxTemp) THEN
                iaT11(iL) = iaT11(iL)-1   !!move iaT11(iL) down by one
            ELSE
                iaT12(iL) = iaT12(iL)+1   !!move iaT12(iL) up by one
            END IF
        END IF
        t11 = raTSpan(iaT11(iL))  ! lower temperature bound at p1
        t12 = raTSpan(iaT12(iL))  ! upper temperature bound at p1
        IF (abs(t11 - t12) > 1.0e-3) THEN
        ! temperature interpolation weight
            raT12(iL)  = (rT - t11) / (t12 - t11)
            raT11(iL)  = 1.0 - raT12(iL)
            raJT12(iL) = +1/(t12-t11)
            raJT11(iL) = -1/(t12-t11)
        ELSE
            raT12(iL)  = 1.0
            raT11(iL)  = 0.0
            raJT12(iL) = +1/(t12-t11)
            raJT11(iL) = -1/(t12-t11)
        END IF

    ! get temperature tabulation bounding interval [t21, t22] at p2
        DO i2 = 1,kMaxTemp
            raTSpan(i2) = raTOffSet(i2) + raOrig100T(iaP2(iL))
        END DO
        iaT21(iL) = iFindMaxMin(+1,rT,raTSpan,kMaxTemp)
        iaT22(iL) = iFindMaxMin(-1,rT,raTSpan,kMaxTemp)
    !! need iaT21 and iaT22 to be separate, for the jacobian
        IF (iaT21(iL) == iaT22(iL)) THEN
            IF (iaT21(iL) == kMaxTemp) THEN
                iaT21(iL) = iaT21(iL)-1   !!move iaT21(iL) down by one
            ELSE
                iaT22(iL) = iaT22(iL)+1   !!move iaT22(iL) up by one
            END IF
        END IF
        t21 = raTSpan(iaT21(iL))  ! lower temperature bound at p2
        t22 = raTSpan(iaT22(iL))  ! upper temperature bound at p2
        IF (abs(t21 - t22) > 1.0e-3) THEN
        ! temperature interpolation weight
            raT22(iL)  = (rT - t21) / (t22 - t21)
            raT21(iL)  = 1.0 - raT22(iL)
            raJT22(iL) = +1/(t22-t21)
            raJT21(iL) = -1/(t22-t21)
        ELSE
            raT22(iL)  = 1.0
            raT21(iL)  = 0.0
            raJT22(iL) = +1/(t22-t21)
            raJT21(iL) = -1/(t22-t21)
        END IF

        wsum = raP1(iL)*raT11(iL) + raP1(iL)*raT12(iL) + &
        raP2(iL)*raT21(iL) + raP2(iL)*raT22(iL)
        IF (abs(1.0-wsum) >= 1.0e-6) THEN
            print *,iL,raP1(iL),raP2(iL)
            write(kStdErr,*) 'error: PT wsum should be 1, not ',wsum
            CALL DoStop
        END IF

    !!!!!!!!!!!!!!!!!!!!!! water partial pressure !!!!!!!!!!!!!!!!!!!!!!!!
        rQ = raPPart(iL)

    ! get partpressure tabulation bounding interval [t11, t12] at p1
        DO i2 = 1,kMaxWater
            raPPSpan(i2) = raPPoffSet(i2) * raOrig100PP(iaP1(iL))
        END DO
        iaQ11(iL) = iFindMaxMin(+1,rQ,raPPSpan,kMaxWater)
        iaQ12(iL) = iFindMaxMin(-1,rQ,raPPSpan,kMaxWater)
        q11 = raPPSpan(iaQ11(iL))  ! lower partpressure bound at p1
        q12 = raPPSpan(iaQ12(iL))  ! upper partpressure bound at p1
    !! this delta(q) was originally 1e-3, but the q of water can get pretty small!!!
        IF (abs(q11 - q12) > 1.0e-11) THEN
        ! partpressure interpolation weight
            raQ12(iL) = (rQ - q11) / (q12 - q11)
            raQ11(iL) = 1.0 - raQ12(iL)
        ELSE
            raQ12(iL) = 1.0
            raQ11(iL) = 0.0
        END IF

    ! get partpressure tabulation bounding interval [t21, t22] at p2
        DO i2 = 1,kMaxWater
            raPPSpan(i2) = raPPoffSet(i2) * raOrig100PP(iaP2(iL))
        END DO
        iaQ21(iL) = iFindMaxMin(+1,rQ,raPPSpan,kMaxWater)
        iaQ22(iL) = iFindMaxMin(-1,rQ,raPPSpan,kMaxWater)
        q21 = raPPSpan(iaQ21(iL))  ! lower partpressure bound at p2
        q22 = raPPSpan(iaQ22(iL))  ! upper partpressure bound at p2
    !! this delta(q) was originally 1e-3, but the q of water can get pretty small!!!
        IF (abs(q21 - q22) > 1.0e-11) THEN
        ! partpressure interpolation weight
            raQ22(iL) = (rQ - q21) / (q22 - q21)
            raQ21(iL) = 1.0 - raQ22(iL)
        ELSE
            raQ22(iL) = 1.0
            raQ21(iL) = 0.0
        END IF

        qsum = raP1(iL)*raT11(iL)*raQ11(iL) + raP1(iL)*raT12(iL)*raQ11(iL) + &
        raP2(iL)*raT21(iL)*raQ21(iL) + raP2(iL)*raT22(iL)*raQ21(iL) + &
        raP1(iL)*raT11(iL)*raQ12(iL) + raP1(iL)*raT12(iL)*raQ12(iL) + &
        raP2(iL)*raT21(iL)*raQ22(iL) + raP2(iL)*raT22(iL)*raQ22(iL)
        IF (abs(1-qsum) >= 1.0e-6) THEN
            write(kStdErr,*) 'error: PTQ wsum should be 1, not ',qsum
            CALL DoStop
        END IF

        IF (raP1(iL) >= raP2(iL)) THEN
        ! this is the dominant weight for pressure
            iX1 = iaP1(iL)
            iX2 = iaT11(iL)
            iX3 = iaQ11(iL)
            rX1 = raP1(iL)
            rX2 = raT11(iL)
            rX3 = raQ11(iL)
        ELSE
            iX1 = iaP2(iL)
            iX2 = iaT12(iL)
            iX3 = iaQ12(iL)
            rX1 = raP2(iL)
            rX2 = raT12(iL)
            rX3 = raQ12(iL)
        END IF
    !        write(kStdWarn,100) iL,rP,iaP1(iL),raP1(iL),rT,iaT11(iL),raT11(iL),rQ,iaQ11(iL),raQ11(iL)
        write(kStdWarn,100) iL,rP,iX1,rX1,rT,iX2,rX2,rQ,iX3,rX3
    END DO

    iL = (kProfLayer-iLowest+1)
    write(kStdWarn,*) ' '
    write(kStdWarn,*) 'Mean pressure 1 weight = ',rSumWgt1/iL
    write(kStdWarn,*) 'Mean pressure 2 weight = ',rSumWgt2/iL
    IF ((rSumWgt1/iL <= 1.0e-3) .AND. (rSumWgt2/iL >= 1.0-1.0e-3)) THEN
        iE = 2
        write(kStdWarn,*) 'xWeights says do fast interp (iSplinetype = 2)'
    ELSE
        iE = 1
        write(kStdWarn,*) 'xWeights says do slow interp (iSplinetype = 1)'
    END IF
    write(kStdWarn,*) 'iSplineType was already determined to be ',iSPlineType
    IF (abs(iE) /= abs(iSplinetype)) THEN
        write(kStdErr,*) 'OOPS : Determining whether SLOW or FAST interp!'
        write(kStdErr,*) '     : SetSplineType determined ',iSplinetype
        write(kStdErr,*) '     : xWeights      determined ',iE
        write(kStdErr,*) ' ..... run proceeding with SetSplineType = ',iE
        iSplinetype = iE
    END IF
    write(kStdWarn,*) ' '

    write(kStdWarn,*) '**********************************************************************'
    write(kStdWarn,*) ' '
    write (kStdWarn,*) 'caVersion (tells you template sedded to make kcartaparam.f90) = '
    write(kStdWarn,222) caVersion(1:80)
    write(kStdWarn,*) ' '
    write(kStdWarn,*) 'The following paths hold for kCompressed Database : '
    write(kStdWarn,*) '  (unless explicitly overridden by nm_spectra)'
    write(kStdWarn,120) kWaterPath(1:80)
    write(kStdWarn,121) kWaterIsotopePath(1:80)
    write(kStdWarn,122) kCO2Path(1:80)
    IF (kCO2_UMBCorHARTMAN == +1) THEN
        write(kStdWarn,*) '  kCO2_UMBCorHARTMAN = ',kCO2_UMBCorHARTMAN,' ==> UMBC CO2 linemixing ;;; ChiFile = '
        write(kStdWarn,222) kChiFile(1:80)
    ELSEIF (kCO2_UMBCorHARTMAN == -1) THEN
        write(kStdWarn,*) '  kCO2_UMBCorHARTMAN = ',kCO2_UMBCorHARTMAN,' ==> LBLRTM CO2 linemixing'
    END IF
    write(kStdWarn,123) kCompPath(1:80)
    write(kStdWarn,124) kCKDPath(1:80)
    write(kStdWarn,125) kCKD_Compr_Path(1:80)
    write(kStdWarn,*) ' '
    write(kStdWarn,126) kOrigRefPath(1:80)
          
    120 FORMAT('Water MolGas  ',A80)
    121 FORMAT('Water HDO     ',A80)
    122 FORMAT('CO2           ',A80)
    222 FORMAT('     ',A80)
    123 FORMAT('Other gases   ',A80)
    124 FORMAT('kCKDPath      ',A80)
    125 FORMAT('kCKDCompPath  ',A80)
    126 FORMAT('RefGas Profs  ',A80)
    write(kStdWarn,*) ' '
    write(kStdWarn,*) '**********************************************************************'
          
    write(kStdWarn,*) ' '
          
    90 FORMAT(F12.5)
    111 FORMAT(I3,' ',3(F9.5,' ',I3,' ',F9.5,' '))
    100 FORMAT(I3,' ',1(ES12.4,' ',I3,' ',F9.5,' '),1(F9.5,' ',I3,' ',F9.5,' '), &
    &  1(ES12.4,' ',I3,' ',F9.5,' '))

    RETURN
    end SUBROUTINE xWeights

!************************************************************************
! this subroutine checks to see if we can keep the iSPlinetype = +2 going
    SUBROUTINE Check_xWeights_iSplinetype(iL,iaP1,iaP2,raP1,raP2)

    IMPLICIT NONE
    include '../INCLUDE/kcartaparam.f90'

! inputs
    INTEGER :: iL
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)

    INTEGER :: iI
    REAL ::    rP,rQ
! if iSplinetype == +/-2 then
!   iL = iaP1(iL) = iaP2(iL), raP1(iL) = 0, raP2(iL) = 1

!! if iSplinetype == 2, asuming we can do FAST interp
!! so either iaP1 or iaP2 must equal iL
    IF ((iaP1(iL) /= iL) .AND. (iaP2(iL) /= iL)) THEN
        write(kStdErr,*) 'iaP1 =  '
        write(kStdErr,*) (iaP1(iI),iI=1,kProfLayer)
        write(kStdErr,*) 'iaP2 =  '
        write(kStdErr,*) (iaP2(iI),iI=1,kProfLayer)
        write(kStdErr,*) 'raP1 =  '
        write(kStdErr,*) (raP1(iI),iI=1,kProfLayer)
        write(kStdErr,*) 'raP2 =  '
        write(kStdErr,*) (raP2(iI),iI=1,kProfLayer)
        write(kStdErr,*) 'need iaP1(iL) or iaP2(iL) == iL',iL,iaP1(iL),iaP2(iL)
        Call DoStop
    END IF

    IF (raP1(iL) > raP2(iL)) THEN
        rP = raP2(iL)
        rQ = raP1(iL)
    ELSE
        rQ = raP2(iL)
        rP = raP1(iL)
    END IF

!! if kRTP == 6, the pavg is coming straight from LBLRTM layers code,
!! which may not agree with KLAYERS code
    IF (abs(rP-0.0) >= 1.0e-3) THEN
        write(kStdErr,*) 'WARNING need the smaller weight ~ 0, not ',rP
        IF (kRTP /= -6) THEN
            CALL DoStop
        END IF
    END IF

    IF (abs(rQ-1.0) >= 1.0e-3) THEN
        write(kStdErr,*) 'WARNING need the larger weight ~ 1, not ',rQ
        IF (kRTP /= -6) THEN
            CALL DoStop
        END IF
    END IF

!! if all this is satisfied, go ahead
    iaP1(iL) = iL
    iaP2(iL) = iL
    raP1(iL) = rP
    raP2(iL) = rQ

    RETURN
    end SUBROUTINE Check_xWeights_iSplinetype

!************************************************************************
! this function finds bounding intervals
    INTEGER FUNCTION iFindMaxMin(iWhich,rV,raArray,iLen)

    IMPLICIT NONE
    include '../INCLUDE/kcartaparam.f90'

! input
!   iWhich   = +1 for finding iX so that raArray(iX) < rV
!              -1 for finding iX so that raArray(iX) > rV
!   rV       = this is input value; find values in raArray which bracket this
!   raArray  = array in which to search for values
!   iLen     = array length

    INTEGER :: iWhich,iLen
    REAL :: raArray(*),rV

    INTEGER :: iOut,iAscOrDsc,iJ

    IF ((raArray(1) <= raArray(2)) .AND. (raArray(2) <= raArray(3)) &
     .AND. (raArray(iLen-2) <= raArray(iLen-1)) &
     .AND. (raArray(iLen-1) <= raArray(iLen))) THEN
    !!! array values incrs with index eg temperatures increase with offset
    !!  raArray(iX,iWhich = -1) < rV < raArray(iY,iWhich = +1); iX < iY
        iAscorDsc = +1
    ELSEIF ((raArray(1) >= raArray(2)) .AND. (raArray(2) >= raArray(3)) &
         .AND. (raArray(iLen-2) >= raArray(iLen-1)) &
         .AND. (raArray(iLen-1) >= raArray(iLen))) THEN
    !!! array values decrs with index eg pressure decreasing with height
    !!  raArray(iX,iWhich = -1) < rV < raArray(iY,iWhich = +1); iX > iY
        iAscorDsc = -1
    ELSE
        write(kStdErr,*) 'iFindMaxMin : array values must increase OR decrease'
        CALL DoStop
    END IF

    iOut = -1
    IF (iAscOrDsc == -1) THEN
        IF (iWhich == +1) THEN
            iOut = 1
            10 CONTINUE
            IF ((raArray(iOut) > rV) .AND. (iOut < iLen)) THEN
                iOut = iOut + 1
                GOTO 10
            END IF
        ELSEIF (iWhich == -1) THEN
            iOut = iLen
            20 CONTINUE
            IF ((raArray(iOut) < rV) .AND. (iOut > 1)) THEN
                iOut = iOut - 1
                GOTO 20
            END IF
        END IF
    ELSEIF (iAscOrDsc == +1) THEN
        IF (iWhich == -1) THEN
            iOut = 1
            30 CONTINUE
            IF ((raArray(iOut) < rV) .AND. (iOut < iLen)) THEN
                iOut = iOut + 1
                GOTO 30
            END IF
        ELSEIF (iWhich == +1) THEN
            iOut = iLen
            40 CONTINUE
            IF ((raArray(iOut) > rV) .AND. (iOut > 1)) THEN
                iOut = iOut - 1
                GOTO 40
            END IF
        END IF
    END IF

    IF (iOut == -1) THEN
        write(kStdErr,*) 'iFindMaxMin not set?? iWhich requires +/-1 ',iWhich
        Call DoStop
    END IF

!      DO iJ = 1,iLen
!        if (rV .LT. raArray(iJ)) print *,iJ,iWhich,rV,raArray(iJ),'---',iOut
!        if (rV .GE. raArray(iJ)) print *,iJ,iWhich,rV,raArray(iJ),'+++',iOut
!      END DO
          
    iFindMaxMin = iOut
     
    RETURN
    end FUNCTION iFindMaxMin

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! iGasID = 1 (WATER), rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
! the spline is done wrt partial pressure, while the scaling is done
! wrt partial pressure
    SUBROUTINE xwater(iGasID,rFileStartFr,iTag,iActualTag,iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPPart,raRPart,raPTemp,raRTemp, &
    iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE
    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! iGasID     = GASID ==1 for water
! rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
! iTag       = current k-comp block of 25 cm-1 that is being processed
! iProfLayer = number of layers in profile === kProfLayer
! iL,iU      = min/max layer number (=1,kMaxlayer)
! daaAbs     = final uncompressed abs coefficient for gas iGasID
! iErr       = errors (mainly associated with file I/O, could be associated
!              with incorrect number of layers in compresse database etc)
! raP/RAmt   = arrays containing actual/reference gas amounts
! raP/RPart  = arrays containing actual/reference gas partial pressures
! raP/RTemp  = arrays containing actual/reference gas temperatures
! daaDQ      = analytic Jacobian wrt amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ      = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    INTEGER :: iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
    INTEGER :: iProfileLayers,iSPlineType
    REAL :: raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer)
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! the weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

! local variables associated with uncompressing the water database files
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp),iLowerOrUpper
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX1(kMaxK,kMaxTemp,kMaxLayer), &
    daaaKX2(kMaxK,kMaxTemp,kMaxLayer), &
    daaaKX3(kMaxK,kMaxTemp,kMaxLayer), &
    daaaKX4(kMaxK,kMaxTemp,kMaxLayer), &
    daaaKX5(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    INTEGER :: iDefault,iMultiplyHeavyWater

    IF ((iGasID /= 1) .AND. (iGasID /= kNewGasHi+1)) THEN
        write(kStdErr,*) 'Expecting to read in water profile/database'
        iErr=1
        CALL DoSTOP
    END IF

    iIOUN = kCompUnit
    CALL CompFileName(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    CALL rdcompwater(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts, &
    iNLay,iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, &
    daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUX)

! check that the file has the data for the correct gas
    IF (iFileGasID /= iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
    END IF

! check that the data file has the right number of layers ===== AIRS layers
    IF (iNLay /= kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
        1010 FORMAT('Error! file : ',/,A120,/, &
        'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
    END IF

! kGenln2Water   = self broadening correction for water, using interpolation
!                  in water partial pressure (+1)
!                = just do what Genln2 does (which would be the same as the
!                  uncompresssion for CO2 (-1)
!      INTEGER kGenln2Water
!      PARAMETER (kGenln2Water=+1)

! interpolate compressed data in temperature, and then in partial pressure,
! to get abs coeff matrix
    IF (kGenln2Water > 0) THEN
    ! e worry about the self broadening corrections
    ! his is pretty good
        IF (kJacobian > 0) THEN
            IF (abs(iSplineType) == 2) THEN
            !! very fast
                CALL x2GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset, &
                daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx, &
                raPTemp,raRTemp,raPPart,raRPart,iaTsort, &
                iNk,iKm,iKn,iUm,iUn,daaDQ,daaDT,iDoDQ,iGasID,pProf, &
                iProfileLayers,iSPlineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            ELSE
            !! fast
                CALL xGetAbsCoeffWaterJAC(daaAbsCoeff,daToffset, &
                daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx, &
                raPTemp,raRTemp,raPPart,raRPart,iaTsort, &
                iNk,iKm,iKn,iUm,iUn,daaDQ,daaDT,iDoDQ,iGasID,pProf, &
                iProfileLayers,iSPlineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            END IF
        ELSE
            iLowerOrUpper = -1
            IF (abs(iSplineType) == 2) THEN
            !! very fast
                CALL x2GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset, &
                daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,raPTemp, &
                raRTemp,raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
                pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            ELSE
            !! fast
                CALL xGetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset, &
                daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,raPTemp, &
                raRTemp,raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
                pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            END IF
        END IF

    ! because iKtype=2 do any necessary jacobians calcs HERE!
        IF (kJacobian > 0) THEN
            IF (iDoDQ > 0)  THEN
                IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                    CALL FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDQ,raPAmt)
                END IF
            END IF
            IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
            (kActualJacs == 32) .OR. &
            (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
                CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
            END IF
        END IF
            
    ELSE
    ! Genln2Water .LT. 0  ==> do same uncompression as for CO2
    ! e do not worry about the self broadening corrections
    ! his is not very good at all!
    ! interpolate compressed data in temperature, to get abs coeff matrix
        IF (kJacobian >= 0) THEN
            IF (abs(iSplineType) == 2) THEN
            !! very fast
                CALL x2GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
                raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
                daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)
            ELSE
            !! fast
                CALL xGetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
                raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
                daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            END IF
        ELSE
            iLowerOrUpper = -1
            IF (abs(iSplineType) == 2) THEN
            !! very fast
                CALL x2GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
                raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
                pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)
            ELSE
            !! fast
                CALL xGetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
                raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
                pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)
            END IF
        END IF

    ! because of iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
        IF (kJacobian >= 0) THEN
            IF (iDoDQ > 0) THEN
                IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                    CALL FinalAmtDeriv(daaDQ,iKtype)
                END IF
            END IF
            IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
            (kActualJacs == 32) .OR. &
            (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
                CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
            END IF
        END IF
    END IF

! convert absorption coefficient correctly if necessary
    IF (iKtype == 2) THEN
        CALL RaisePower(daaAbsCoeff)
    END IF

! now compute optical depth = gas amount * abs coeff
    CALL AmtScale(daaAbsCoeff,raPAmt)

    RETURN
    end SUBROUTINE xwater

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! for the COUSIN CO2 files
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE xCousinContribution(iGasID, &
    rFileStartFr,iTag,iActualTag,iProfileLayers,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,pProf, &
    rLTEstrength,iNLTEStart,iLowerOrUpper,daaAbsCoeff,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! iGasID     = GASID
! rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
! iTag       = current k-comp block of 25 cm-1 that is being processed
! iProfLayer = number of layers in profile === kProfLayer
! iL,iU      = min/max layer number (=1,kMaxlayer)
! daaAbs     = final uncompressed abs coefficient for gas iGasID
! iErr       = errors (mainly associated with fOAile I/O, could be associated
!              with incorrect number of layers in compresse database etc)
! raP/RAmt   = arrays containing actual/reference gas amounts
! raP/RPart  = arrays containing actual/reference gas partial pressures
! raP/RTemp  = arrays containing actual/reference gas temperatures
! new stuff
! rLTEStrength = tells the weight ~ 1.1212
! iNLTEStart = tells where to replace LineMix spectra with Cousin spectra
    INTEGER :: iGasID,iErr,iL,iU,iTag,iActualTag
    INTEGER :: iProfileLayers,iLowerOrUpper,iSPlineType
    REAL :: raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    REAL :: rLTEStrength
    INTEGER :: iNLTEStart
! the weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)
        
! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp),iFr
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    DOUBLE PRECISION :: daaCousin(kMaxPts,kProfLayer)

    IF (iGasID /= 2) THEN
        write(kStdErr,*) 'This is only for CO2!!!'
        CALL DoStop
    END IF

    iIOUN = kCompUnit
    CALL CompFileName(-1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay, &
    iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, &
    daaaKX,daaUX)

! check that the file has the data for the correct gas
    IF (iFileGasID /= iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
    END IF

! check that the data file has the right number of layers
    IF (iNLay /= kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
        1010 FORMAT('Error! file : ',/,A120,/, &
        'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
    END IF

! interpolate compressed data in temperature, to get abs coeff matrix
    IF (kJacobian >= 0) THEN
        write(kStdErr,*) 'Cannot do Jacobians and willy nilly replace '
        write(kStdErr,*) 'linemix spectroscopy with cousin spectroscopy'
        CALL DoStop
    ELSE
        iLowerOrUpper = -1    !!!!only use 100 AIRS layers; nuthin above
        IF (abs(iSplineType) == 2) THEN
        !! very fast
            CALL x2GetAbsCoeffNOJAC(daaCousin,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)

        ELSE
        !! fast
            CALL xGetAbsCoeffNOJAC(daaCousin,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)

        END IF
    END IF

! convert absorption coefficient correctly if necessary
    IF (iKtype == 2) THEN
        CALL RaisePower(daaCousin)
    END IF

! now compute optical depth = gas amount * abs coeff
    CALL AmtScale(daaCousin,raPAmt)

! now multiply all spectra by scaling factor, and replace daaGasAb as required
    write (kStdWarn,*) 'Replacing LINEMIX spectra with COUSIN spectra'
    write (kStdWarn,*) 'start layer, strength = ',iNLTEStart, &
    abs(rLTEStrength)
    DO iL = iNLTEStart,kProfLayer
        DO iFr = 1,kMaxPts
            daaAbsCoeff(iFr,iL) = daaCousin(iFr,iL) * abs(rLTEStrength)
        END DO
    END DO

    RETURN
    end SUBROUTINE xCousinContribution

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE xothergases(iGasID,rFileStartFr,iTag,iActualTag, &
    iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! iGasID     = GASID
! rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
! iTag       = current k-comp block of 25 cm-1 that is being processed
! iProfLayer = number of layers in profile === kProfLayer
! iL,iU      = min/max layer number (=1,kMaxlayer)
! daaAbs     = final uncompressed abs coefficient for gas iGasID
! iErr       = errors (mainly associated with file I/O, could be associated
!              with incorrect number of layers in compresse database etc)
! raP/RAmt   = arrays containing actual/reference gas amounts
! raP/RPart  = arrays containing actual/reference gas partial pressures
! raP/RTemp  = arrays containing actual/reference gas temperatures
! daaDT      = analytic Jacobian wrt temperature
! daaDQ      = analytic Jacobian wrt amount
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    INTEGER :: iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
    INTEGER :: iProfileLayers,iSplineType
    REAL :: raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)
        
! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp)
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    INTEGER :: iLowerOrUpper,iJ,iJunk

    iIOUN = kCompUnit
    IF ((iGasID /= kNewGasLo) .AND. (iGasID /= kNewGasHi)) THEN
        CALL CompFileName(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    ELSE
        CALL CompFileNameCKD(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    END IF

    CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay, &
    iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, &
    daaaKX,daaUX)

! check that the file has the data for the correct gas
    IF (iFileGasID /= iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
    END IF

! check that the data file has the right number of layers
    IF (iNLay /= kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
        1010 FORMAT('Error! file : ',/,A120,/, &
        'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
    END IF

! interpolate compressed data in temperature, to get abs coeff matrix
    IF (kJacobian >= 0) THEN
        IF (abs(iSplineType) == 2) THEN
        !! very fast
            CALL x2GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
            daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
        !! fast
            CALL xGetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
            daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        END IF
    ELSE
        iLowerOrUpper = -1
        IF (abs(iSplineType) == 2) THEN
        !! very fast
            CALL x2GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
        !! fast
            CALL xGetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        END IF
    END IF

! because of the iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
    IF (kJacobian >= 0) THEN
        IF (iDoDQ > 0) THEN
            IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                CALL FinalAmtDeriv(daaDQ,iKtype)
            END IF
        END IF
        IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
        (kActualJacs == 32) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
            CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
        END IF
    END IF

! convert absorption coefficient correctly if necessary
    IF (iKtype == 2) THEN
        CALL RaisePower(daaAbsCoeff)
    END IF

! now compute optical depth = gas amount * abs coeff
    CALL AmtScale(daaAbsCoeff,raPAmt)
          
! print *,iGasID,int(rFileStartFr),kAltDirs

    IF (iGasID == 2) THEN
        CALL multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)
    END IF
           
    RETURN
    end SUBROUTINE xothergases

!************************************************************************
! this subroutine calls the routines to read CKD continmuum k-compressed data
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE xCKDgases(iGasID,rFileStartFr,iTag,iActualTag, &
    iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! iGasID     = GASID
! rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
! iTag       = current k-comp block of 25 cm-1 that is being processed
! iProfLayer = number of layers in profile === kProfLayer
! iL,iU      = min/max layer number (=1,kMaxlayer)
! daaAbs     = final uncompressed abs coefficient for gas iGasID
! iErr       = errors (mainly associated with file I/O, could be associated
!              with incorrect number of layers in compresse database etc)
! raP/RAmt   = arrays containing actual/reference gas amounts
! raP/RPart  = arrays containing actual/reference gas partial pressures
! raP/RTemp  = arrays containing actual/reference gas temperatures
! daaDT      = analytic Jacobian wrt temperature
! daaDQ      = analytic Jacobian wrt amount
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    INTEGER :: iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
    INTEGER :: iProfileLayers,iSplineType
    REAL :: raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)
        
! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp)
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge,WhichGasPosn
    INTEGER :: iDefault,iLowerOrUpper,iJ,iJunk
    REAL :: raRXAmt(kProfLayer),raPXAmt(kProfLayer)
            
    DO iJ = 1,kProfLayer
        raRXAmt(iJ) = 1.0
        raPXAmt(iJ) = 1.0
    END DO

    iIOUN = kCompUnit
    IF ((iGasID /= kNewGasLo) .AND. (iGasID /= kNewGasHi)) THEN
        write (kStdErr,*) 'xCKDgases is only for gases 101,102'
        CALL DoStop
    ELSE
        CALL CompFileNameCKD(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    END IF
    CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay, &
    iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, &
    daaaKX,daaUX)

! check that the file has the data for the correct gas
    IF (iFileGasID /= iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
    END IF

! check that the data file has the right number of layers
    IF (iNLay /= kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
        1010 FORMAT('Error! file : ',/,A120,/, &
        'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
    END IF

! interpolate compressed data in temperature, to get abs coeff matrix
    IF (kJacobian >= 0) THEN
        IF (abs(iSplineType) == 2) THEN
        !! very fast
            CALL x2GetAbsCoeffJAC_CKD(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
            daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
            write(kStdErr,*) 'klayers different from DataBase'
            write(kStdErr,*) 'Prefer using older CKD code'
            CALL DoStop
        !! fast
            CALL xGetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
            daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        END IF
    ELSE
        iLowerOrUpper = -1
        IF (abs(iSplineType) == 2) THEN
        !! very fast
            CALL x2GetAbsCoeffNOJAC_CKD(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
            write(kStdErr,*) 'klayers different from DataBase'
            write(kStdErr,*) 'Prefer using older CKD code'
            CALL DoStop
        !! fast
            CALL xGetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        END IF
    END IF

! because of the iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
    IF (kJacobian >= 0) THEN
        IF (iDoDQ > 0) THEN
            IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                CALL FinalAmtDeriv(daaDQ,iKtype)
            END IF
        END IF
        IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
        (kActualJacs == 32) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
            CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPXAmt)
        END IF
    END IF

! convert absorption coefficient correctly if necessary
    IF (iKtype == 2) THEN
        CALL RaisePower(daaAbsCoeff)
    END IF

! now compute optical depth = gas amount * abs coeff
!      CALL AmtScale(daaAbsCoeff,raPXAmt)
! no need as raPXAmt == 1
           
    RETURN
    end SUBROUTINE xCKDgases

!************************************************************************
END MODULE kcoeff_FAST
