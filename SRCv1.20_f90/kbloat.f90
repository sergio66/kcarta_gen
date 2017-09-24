! Copyright 2002
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kbloat

USE basic_common
USE spline_and_sort_and_common
USE rad_misc
USE rad_main
USE s_writefile

IMPLICIT NONE

CONTAINS

! bloated output routines for high res NLTE
!************************************************************************
! this subroutine adds on current LTE gas optical depths to cumulative mixed
! path optical depth
! almost the same as SUBROUTINE Accumulate(raaSum,raaGas,raaMix,iGas,iIpmix)
    SUBROUTINE AccumulateForBloat(raaSum,raaGas,raaMix,iGas,iIpmix)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raaSum     = cumulative spectra associated with the mixed paths
! raaGas     = current gas absorption spectra
! raaMix     = mixing table info from *MIXFIL
! iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
! iIpmix     = which of the mixed paths are being considered
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iIpmix,iGas
    REAL :: raaSum(kMaxPts,kProfLayer)
    REAL :: raaGas(kMaxPts,kProfLayer)
     
    INTEGER :: iFreq,iL
    REAL :: rL

    IF (iIPMIX > kProfLayer) THEN
        write(kStdErr,*) 'AccumulateForBloat assumes 1..kProfLayer mix paths'
        CALL DoStop
    END IF

! find out which of the 100 layers is associated with this mixed path
    iL = MP2Lay(iIpmix)
     
! find the weight
    rL = raaMix(iIpmix,iGas)
     
! add on contribution of the iGas th gas to the iIpmix th row of raaSum
    DO iFreq=1,kMaxPts
        raaSum(iFreq,iIpmix)=raaSum(iFreq,iIpmix)+rL*raaGas(iFreq,iL)
    END DO
     
    RETURN
    end SUBROUTINE AccumulateForBloat

!************************************************************************
! this adds together stuff so as to output Mixed paths at high resolution
    SUBROUTINE SumBloatMP( &
    daFreqBloat,raFreq,raaCo2_LTE,raaRestOfLTEGases,iNLTEStart, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)
    REAL :: raaCo2_LTE(kMaxPts,kProfLayer)
    REAL :: raFreq(kMaxPts)
    INTEGER :: iNLTEStart

! input/output params
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
    DOUBLE PRECISION :: daFreqT(kMaxPts),daTemp(kMaxPts)
    DOUBLE PRECISION :: dY,daY(kBloatPts)
    INTEGER :: iFr,iL

    DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
    END DO

! in SetPlanckCoeffBloat we have already done, from iNLTEStart to kProflayer
!       daaNLTEGasAbCoeffBloat(iFr,iL)+blah_from_raaRestOfLTEGases(iFr,iL)
! so no need to redo it here
! however we ned to do it from 1 t- NLTEStart-1

    DO iL = 1,iNLTEStart-1
        DO iFr = 1,kMaxPts
        !! we have already stored daaNLTEGasAbCoeffBloat= NLTE_CO2 + LTE_CO2
        !! so we need to add on LTE_others to get total mixed path
            daTemp(iFr)  = raaRestOfLTEGases(iFr,iL)*1.0d0
        END DO
        CALL dlinear_smart(daFreqT,daTemp, kMaxPts,daFreqBloat,daY, kBloatPts)
        DO iFr = 1,kBloatPts
            daaSumNLTEGasAbCoeffBloat(iFr,iL) = &
            max(daaNLTEGasAbCoeffBloat(iFr,iL)+daY(iFr),0.0d0)
        END DO
    END DO

    DO iL = iNLTEStart,kProfLayer
        DO iFr = 1,kBloatPts
            daaSumNLTEGasAbCoeffBloat(iFr,iL) = &
            max(daaNLTEGasAbCoeffBloat(iFr,iL),0.0d0)
        END DO
    END DO

    RETURN
    end SUBROUTINE SumBloatMP

!************************************************************************
! this subroutine changes double --> real for the KCARTA atmosphere
! this computes the cumulative multiplication modification to Planck
! this is for the usual 100 kCARTA database layers
! however, this is at the bloated HIGH resolution!!!!!!!
    SUBROUTINE SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer, &
    raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    raaSumAbCoeff,raaPlanckCoeff, &
    raaRestOfLTEGases,raaCO2_LTE, &
    daFreqBloat,daaSumNLTEGasAbCoeffBloat, &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL :: raFreq(kMaxPts)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
! cumulative nonlte gas optical depths
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer) !!planck mod, so far
    REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)             !!mixed path abscoeff
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)          !!rest of LTE gases
    REAL :: raaCO2_LTE(kMaxPts,kProfLayer)                 !!LTE part of CO2
    INTEGER :: iAtm                                        !!which atmosphere
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)             !!set of mixed paths
    INTEGER :: iNLTEStart                                !!where nonLTE starts
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
! for NLTE gas, this contains cumulative optical depths = LTE + NLTE
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
! input/output vars
! at input, only contains NLTE opt depths for CO2
! at output, contains (LTE+NLTE) opt depths for CO2, plus other LTE gases
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)

! local variables
    DOUBLE PRECISION :: daFreqT(kMaxPts)
    DOUBLE PRECISION :: dY, daY(kBloatPts), daTemp(kMaxPts)
    DOUBLE PRECISION :: dY1,daY1(kBloatPts),daTemp1(kMaxPts)
    INTEGER :: iFr,iL,iA,iAtmTemp,iM,iStart

    DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
    END DO

    iAtmTemp = 1
    IF (iAtmTemp /= iatm) THEN
        write(kStdErr,*) 'oh oh iAtmTemp /= iatm ',iAtmTemp,iatm
        write(kStdErr,*) 'kCARTA assumes the NLTE is for and only for atm #1'
        CALL DoStop
    END IF
    iA = 1    !!!assume current atmospher uses mixed path layers 1-100
    iM = iaaRadLayer(iAtmTemp,1)
           
! ind which set of Mixed Paths current atm uses eg 1-100, 101-200 etc
    DO iL = 1,kProfLayer
        IF (iM <= kProfLayer*iL) THEN
            iA = iL
            GOTO 10
        END IF
    END DO
    10 CONTINUE
           
! initialize to 1.0 below where NLTE starts!!!!
    DO iL = 1,iNLTEStart-1
        DO iFr = 1,kBloatPts
            daaPlanckCoeffBloat(iFr,iL) = 1.0
        END DO
    END DO

    DO iL = iNLTEStart,kProfLayer
        iM = iL + (iA-1)*kProfLayer
        DO iFr = 1,kMaxPts
        !!!the LTE part is easily gotten from raaRestOfLTEGases + raaCo2_LTE
            daTemp(iFr)  = raaRestOfLTEGases(iFr,iM)*1.0d0 + &
            raaCo2_LTE(iFr,iM)*1.0d0
            daTemp1(iFr) = raaRestOfLTEGases(iFr,iM)*1.0d0
        END DO
    !!!since abs(LTE) << abs(CO2_NLTE), can get away with linear interp
        CALL dlinear_smart(daFreqT,daTemp, kMaxPts,daFreqBloat,daY, kBloatPts)
        CALL dlinear_smart(daFreqT,daTemp1,kMaxPts,daFreqBloat,daY1,kBloatPts)
        DO iFr = 1,kBloatPts
        ! dY is thus the LTE component of numerator
            dY = daY(iFr)
        ! update this var so it gives the TOTAL optical depth at the layer
        !! or Sum(NLTE_CO2 + LTE_CO2 + LTE_others)
        !! we have already stored daaNLTEGasAbCoeffBloat = NLTE_CO2 + LTE_CO2
        !! so we need to add on LTE_others
            dY = max(daaNLTEGasAbCoeffBloat(iFr,iL)+daY1(iFr),0.0d0)
            daaNLTEGasAbCoeffBloat(iFr,iL) = dY + dDeltaNLTE
        ! for Planck, add on LTE component
            dY = max(daY(iFr)  + daaPlanckCoeffBloat(iFr,iL),0.0d0)
            dY = max(daY1(iFr) + daaPlanckCoeffBloat(iFr,iL),0.0d0)
            dY = dY + dDeltaNLTE
            daaPlanckCoeffBloat(iFr,iL) = dY/daaNLTEGasAbCoeffBloat(iFr,iL)
        END DO
    END DO

    1010 FORMAT(I6,' ',3(F20.7,' '))
    1111 FORMAT(I6,' ',3(F20.7,' '))

    RETURN
    end SUBROUTINE SetPlanckCoeffBloat

!************************************************************************
! this subroutine simply bloats stuff, on a 1-1 correspondance using i1,i2
    SUBROUTINE AccumulateBloat(daHighX,daHighY,iHigh,daLowX,daLowY,iLow, &
    iWideMeshLoop,i1,i2,iType,daBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
! iType = +1 ==> use daHighY at kMaxPtsBox resolution
!                  use the iHigh points from this vector
! iType = -1 ==> use daLowY  at kMaXpts    resolution
!                  use the iCount points from this vector
    DOUBLE PRECISION :: daHighY(kMaxPtsBox),daHighX(kMaxPtsBox)
    DOUBLE PRECISION :: daLowY(kMaxPtsBox),daLowX(kMaxPtsBox)
    INTEGER :: i1,i2,iType,iHigh,iLow,iWideMeshLoop
! output params
    DOUBLE PRECISION :: daBloat(kBloatPts)

! local vars
    INTEGER :: iI,iJ,iK,iS,iE,iOffset
    DOUBLE PRECISION :: daTemp(kBloatPts)

    IF (iType == 1) THEN
        iOffSet = (iWideMeshLoop-1)*2000
        DO iJ = 1,iHigh-1
            daBloat(iJ+iOffSet) = daBloat(iJ+iOffSet) + daHighY(iJ)
        END DO
    ELSEIF (iType == -1) THEN
        CALL dlinear_smart(daLowX,daLowY,iLow,daHighX,daTemp,5*(i2-i1+1))
        iS = -4
        iE = 0
        DO iI = 1,i2-i1+1
            iS = iS + 5
            iE = iE + 5
            iOffSet = (i1-1)*5
            DO iJ = iS,iE
                daBloat(iJ+iOffSet) = daBloat(iJ+iOffSet) + daTemp(iJ)
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE AccumulateBloat

!************************************************************************
! this is the driver
! to bloat up the weak LTE lines (and the stronger LTE lines low down
! in the atm) to higher resolution

! this subroutine bloats up the optical depths, at double precision
! can bloat up all layers, or just select layers
! usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
! at the higher altitudes
    SUBROUTINE BloatCoeffsDriver(iTag,iActualTag,rFileStartFr,raFreq,daaGasAbCoeff, &
    raaRestOfLTEGases,raaCO2_LTE, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat,iNLTEStart,rNLTEStrength, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
    iGas,iGasID,iL,iU,iUseWeakBackGnd, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    pProf,iProfileLayers, &
    raPAmt,raPTemp,raPPress,raPPartPress,iSplineType)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL ::             raFreq(kMaxPts)                        !!kCARTA res
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
! at kCARTA res
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)            !!high res pts
    REAL :: rNLTEStrength                                 !!mixed path strength
    INTEGER :: iNLTEStart                                 !!where NLTE starts
    INTEGER :: iTag,iActualTag                            !!tells info about chunk
    REAL :: rFileStartFr                                  !!start point of chunk

! pProf       = actual layers (from kLAYERS) avg pressure, in iProfileLayers
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (= kProfLayer)
! iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
! iaCont    = whether or not to do continuum calculation .. iaCont(iCount)
! caXecF    = file name of cross section data
! daaTemp   = matrix containing the uncompressed k-spectra
! raVtemp   = vertical temperature profile for the Mixed paths
! iVTSet    = has the vertical temp been set, to check current temp profile
! rFileStartFr   = which k-comp file chunk to uncompress
! iTag      = which k-comp file chunk to uncompress
! raFreq    = wavenumber array
! iErr      = errors (mainly associated with file I/O)
! daaDQ     = analystic Jacobian wrt water amount
! daaDT     = analystic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPress(kProfLayer),raPPartPress(kProfLayer)
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raVTemp(kProfLayer)
    INTEGER :: iCount,iGas,iGasID,iL,iU,iErr,iRefLayer,iSplineType
    INTEGER :: iUseWeakBackGnd
          
! output params
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)
    REAL ::  raaCO2_LTE(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iLL

    write(kStdWarn,*) 'Using std 0.0025 -> 0.0005 code for "highres" bloat'
    CALL BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff, &
    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd, &
    raPAmt,raPTemp,raPPress,raPPartPress, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat,iNLTEStart,rNLTEStrength, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

    RETURN
    end SUBROUTINE BloatCoeffsDriver

!************************************************************************
! this just bloats up the weak LTE lines (and the stronger LTE lines low down
! in the atm) to higher resolution
! this subroutine bloats up the optical depths, at double precision
! can bloat up all layers, or just select layers
! usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
! at the higher altitudes

    SUBROUTINE BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff, &
    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd, &
    raPAmt,raPTemp,raPPress,raPPartPress, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat,iNLTEStart,rNLTEStrength, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    REAL ::             raFreq(kMaxPts)                        !!kCARTA res
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
! at kCARTA res
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)            !!high res pts
    REAL :: rNLTEStrength                                 !!mixed path strength
    INTEGER :: iNLTEStart                               !!where NLTE starts
    INTEGER :: iTag                                     !!tells info about chunk
    REAL ::    rFileStartFr                             !!start point of chunk
    INTEGER :: iUseWeakBackGnd                          !!go thru with this?
    REAL :: raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPress(kProfLayer),raPPartPress(kProfLayer)

! output params
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)
    REAL ::  raaCO2_LTE(kMaxPts,kProfLayer)

! local variables
    DOUBLE PRECISION :: daTemp(kMaxPts),raaTemp(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daFreqT(kMaxPts)
    DOUBLE PRECISION :: daY(kBloatPts)
    DOUBLE PRECISION :: df,dfFine,dK
    INTEGER :: iFr,iL,iM,iA,iIpMix,iS,iE
    REAL :: raFreqBloat(kBloatPts)
    INTEGER :: iDo,iLinearOrSpline,iDefault,iStart

    iDefault = -1
    iLinearOrSpline = +1    !!do linear .. not good at all
    iLinearOrSpline = -1    !!do spline .. much better, but overall still sux

    IF (iDefault /= iLinearOrSpline) THEN
        write(kStdErr,*) 'iDefault,iLinearOrSpline = ',iDefault,iLinearOrSpline
    END IF

    DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
    END DO

!!! at this point, the code has uncompressed the database, so the lower
!!!   layers contain the lower trop LTE optical depths
!!! in addition the code has replaced the (weak+strong) upper atm LTE
!!!    optical depths with the (weak) upper trop LTE depths
!!! So save the CO2 LTE stuff
    DO iL = 1,kProfLayer
        DO iFr = 1,kMaxPts
            raaCO2_LTE(iFr,iL) = sngl(daaNLTEGasAbCoeff(iFr,iL))
        END DO
    END DO

! <---!!!! --------------------------->
!!! update stuff with the weak LTE lines that have been computed in
!!! by calling LTE_spectra. This is for upper atmosphere
! set coeffs of upper layers
    iS = iNLTEStart
    iE = kProfLayer
    write (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iS,' up'
    DO iL = iS,iE
        DO iFr = 1,kMaxPts
            daTemp(iFr)  = daaNLTEGasAbCoeff(iFr,iL)
        END DO
        IF (iLinearOrSpline == 1) THEN
            CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        ELSE
            CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        END IF
        DO iFr = 1,kMaxPtsBox
            dK = max(daY(iFr),0.0d0)
            daaNLTEGasAbCoeffBloat(iFr,iL)    = dK
            daaSumNLTEGasAbCoeffBloat(iFr,iL) = dK
        END DO
    END DO

    DO iL = iS,iE
        DO iFr = 1,kMaxPts
            daTemp(iFr)  = daaPlanckCoeff(iFr,iL)
        END DO
        IF (iLinearOrSpline == 1) THEN
            CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        ELSE
            CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        END IF
        DO iFr = 1,kMaxPtsBox
            daaPlanckCoeffBloat(iFr,iL) = max(daY(iFr),0.0d0)
        END DO
    END DO

! <---!!!! --------------------------->
!!! update stuff with all LTE lines that have been computed from
!!! uncompressing the kCompressed Database. This is for lower  atmosphere
! set coeffs of lower layers
    iS = 1
    iE = iNLTEStart-1
    write (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iE,' down'
    DO iL = iS,iE
        DO iFr = 1,kMaxPts
            daTemp(iFr)  = daaGasAbCoeff(iFr,iL)
        END DO
        IF (iLinearOrSpline == 1) THEN
            CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        ELSE
            CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        END IF
        DO iFr = 1,kMaxPtsBox
            dK = max(daY(iFr),0.0d0)
            daaNLTEGasAbCoeffBloat(iFr,iL)    = dK
            daaSumNLTEGasAbCoeffBloat(iFr,iL) = dK
        !          daaPlanckCoeffBloat(iFr,iL)       = 1.0d0
        END DO
    END DO

    RETURN
    end SUBROUTINE BloatCoeffsStandard

!************************************************************************
! this subroutine bloats up the LTE coeffs from the upper atmosphere database
    SUBROUTINE bloatUAstuff( &
    daaUpperNLTEGasAbCoeff,daaUpperPlanckCoeff,iUpper, &
    raFreq,daFreqBloat, &
    daaUpperPlanckCoeffBloat,daaUpperNLTEGasAbCoeffBloat, &
    daaUpperSumNLTEGasAbCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iUpper
    REAL :: raFreq(kMaxPts)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
! output params
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local vars
    INTEGER :: iFr,iL,iS,iE,iDefault,iLinearOrSpline
    DOUBLE PRECISION :: dK,daTemp(kMaxPts),daFreqT(kMaxPts),daY(kBloatPts)

    iDefault = +1
    iLinearOrSpline = -1    !!do spline .. worse!
    iLinearOrSpline = +1    !!do linear .. better here!

    IF (iDefault /= iLinearOrSpline) THEN
        write(kStdErr,*) 'iDefault,iLinearOrSpline = ',iDefault,iLinearOrSpline
    END IF

    iS = 1
    iE = iUpper

    DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
    END DO

    DO iL = iS,iE
        DO iFr = 1,kMaxPts
            daTemp(iFr)  = daaUpperNLTEGasAbCoeff(iFr,iL)
        END DO
        IF (iLinearOrSpline == 1) THEN
            CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        ELSE
            CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        END IF
        DO iFr = 1,kMaxPtsBox
            dK = max(daY(iFr),0.0d0)
            daaUpperNLTEGasAbCoeffBloat(iFr,iL)    = dK
            daaUpperSumNLTEGasAbCoeffBloat(iFr,iL) = dK
        END DO
    END DO

    DO iL = iS,iE
        DO iFr = 1,kMaxPts
            daTemp(iFr)  = daaUpperPlanckCoeff(iFr,iL)
        END DO
        IF (iLinearOrSpline == 1) THEN
            CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        ELSE
            CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts)
        END IF
        DO iFr = 1,kMaxPtsBox
            dK = max(daY(iFr),0.0d0)
            daaUpperPlanckCoeffBloat(iFr,iL) = dK
        END DO
    END DO

    RETURN
    end SUBROUTINE bloatUAstuff

!************************************************************************
! this subroutine reads in the precomputed LBL weakbackgnd data
! as you can see this is hardcoded for
!   /home/sergio/SPECTRA/IPFILES/layop_md_va0_usethis_co2'
    SUBROUTINE ReadFixedDataWeakBack_Bloat(iGasID,iStart,iTag,rFileStartFr, &
    raPAmt,raPTemp,raPPress,raPPartPress, &
    rNLTEStrength,daaWeakOptDepthBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPress(kProfLayer),raPPartPress(kProfLayer)
    INTEGER :: iStart,iTag,iGasID
    REAL :: rFileSTartFr
    DOUBLE PRECISION :: daaWeakOptDepthBloat(kBloatPts,kProfLayer)
    REAL :: rNLTEStrength

! local variables
    INTEGER :: iIoun,I,J,iErr, iL, iFr
    INTEGER :: IDGAS, NPTS, NLAY
    DOUBLE PRECISION :: SFREQ, FSTEP
    CHARACTER(80) :: FNAM
    CHARACTER(4) :: ca4
    DOUBLEPRECISION df,sf,dK
    INTEGER :: iDefault
    REAL :: raFreqBloat(kBloatPts)

    df = kaFineFrStep(iTag)
    sf = 1.0d0*nint(rFileStartFr)
    DO iFr = 1,kBloatPts
        raFreqBloat(iFr) = sf + (iFr-1)*df
    END DO

    FNAM = '/carrot/s1/sergio/AIRSCO2/BACKGND_vt_mdH96/hehe_highres_weaklte_'
    I = 80
    10 CONTINUE
    IF (FNAM(I:I) == ' ') THEN
        I = I - 1
        GOTO 10
    END IF
    IF ((nint(rFileStartFr) < 1000) .OR. &
    (nint(rFileStartFr) > 9999)) THEN
        write(kStdErr,*) 'Ooooooerr cannot accept freqs outside 1000 ... 9999'
        CALL DoStop
    END IF
     
    write(ca4,40) nint(rFileStartFr)
    40 FORMAT(I4)
    FNAM(I+1:I+4) = ca4(1:4)
    FNAM(I+5:I+8) = '.dat'

    write(kStdWarn,*) 'In ReadFixedDataWeakBack_Bloat, open file for GasID'
    write(kStdWarn,*) iGasID
    write(kStdWarn,*) FNAM

    iIOUN = kCompUnit
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine ReadFixedDataWeakBack_Bloat'
        WRITE(kStdErr,1010) IERR, FNAM
        1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
        CALL DoSTOP
    ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, NPTS, NLAY
    READ(iIOUN) SFREQ, FSTEP

    IF (IDGAS /= iGasID) THEN
        write(kStdErr,*) 'Wanted gas ',iGasID,'but new data is for gas ',IDGAS
        CALL DoStop
    END IF

    IF (NLAY /= kProfLayer) THEN
        write(kStdErr,*) ' new data has ',NLAY,'layers instead of kProfLayer'
        CALL DoStop
    END IF

    IF (NPTS /= kBloatPts) THEN
        write(kStdErr,*) ' new data has ',NPTS,' freq pts instead of kBloatPts'
        CALL DoStop
    END IF

    IF (abs(df-FSTEP) > 1e-4) THEN
        write(kStdErr,*) ' new data has ',FSTEP,'freq spacing instead of ',df
        CALL DoStop
    END IF

    IF (abs(sf-SFREQ) > 1e-4) THEN
        write(kStdErr,*) ' new data starts at ',SFREQ,' instead of ',sf
        CALL DoStop
    END IF
          
! passed all tests, so read data
!     Read in the optical depths
    DO I=1,kProfLayer
        READ(iIOUN) (daaWeakOptDepthBloat(J,I),J=1,kBloatPts)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

    write (kStdWarn,*) 'Read in HIGH RES WEAK BKGND spectra for gasID ',IDGAS

    RETURN
    end SUBROUTINE ReadFixedDataWeakBack_Bloat

!************************************************************************
! this opens the bloated file four output
    SUBROUTINE OpenOutputBloatFile( &
    iType,iNumLayers,caBloatFile,rFrLow,rFrHigh, &
    iFileIDLo,iFileIDHi,iTag,iTotalStuff, &
    iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    CHARACTER(80) :: caBloatFile
    REAL :: rFrLow,rFrHigh
    INTEGER :: iFileIDLo,iFileIDHi,iTag,iTotalStuff

    INTEGER :: iNumLayers !!number of layers in atm #1
    INTEGER :: iType      !!+1 if regular output; -1 if planck output
    INTEGER :: iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

! local vars
    INTEGER :: iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks
    REAL ::    r1,r2,r3

    iNLTEChunks = -1
    i1 = +123456
    i2 = -123456
    DO iI = 1,iNumNLTEGases
        IF (iaNLTEChunks(iI) >= iNLTEChunks) THEN
            iNLTEChunks = iaNLTEChunks(iI)
        END IF
        DO iJ = 1,iaNLTEChunks(iI)
            IF (iaaNLTEChunks(iI,iJ) <= i1) THEN
                i1 = iaaNLTEChunks(iI,iJ)
            END IF
            IF (iaaNLTEChunks(iI,iJ) >= i2) THEN
                i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
            END IF
        END DO
    END DO
    r1 = i1*1.0
    r2 = i2*1.0
    iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
     
!!!now check to see we actually start from r1, while running kCARTA!
    IF (r1 < rFrLow) THEN
        r1 = rFrLow
        i1 = iFloor(r1)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
    END IF

!!!now check to see we actually do go to r2, while running kCARTA!
!      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
!        r2 = rFrHigh-kaBlSize(iTag)
    IF (r2 > rFrHigh) THEN
        r2 = rFrHigh
        i2 = iFloor(r2)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
    END IF

    IF (iType > 0) THEN
        iIoun1 = kBloatNLTEOut
    ELSEIF (iType < 0) THEN
        iIoun1 = kBloatNLTEPlanck
    ELSE
        write(kStdErr,*) 'Unknown print option for bloated files'
        CALL DoStop
    END IF

! open unformatted file as a fresh file to be written to
    OPEN(UNIT=iIoun1,FILE=caBloatFile,STATUS='NEW', &
    FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
    IF (iFileErr /= 0) THEN
        WRITE(kStdErr,304) iFileErr,iIOUN1,caBloatFile
        write(kStdErr,*)'make sure the file does not exist!'
        CALL DoSTOP
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,*) 'Printing BLOAT FILE summary stats at open : '
    write(kStdWarn,*) caBloatFile
    write(kStdWarn,*) 'Regular freqs etc .... '
    write(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
    write(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
    write(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
    write(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff
    write(kStdWarn,*) 'Bloated  freqs etc .... '
    write(kStdWarn,*) '  iType = = ',iType
    write(kStdWarn,*) '  iNumlayers = ',iNumLayers
    write(kStdWarn,*) '  kBoxCarUse = ',kBoxCarUse
    write(kStdWarn,*) '  df_fine = ',kaFineFrStep(iTag)
    write(kStdWarn,*) 'start,stop high res integers : ',i1,i2
    write(kStdWarn,*) 'number of highres chunks = ',iNLTEChunks
    write(kStdWarn,*) 'start,stop high res doubles : ',r1,r2
    write(kStdWarn,*) ' '

    304 FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file : &
    ',/,A80)

    IF (iType > 0) THEN
        kBloatOutOpen = 1
        write(kStdWarn,*) 'Opened following file for bloated general output :'
        write(kStdWarn,*) caBloatFile
    ELSEIF (iType < 0) THEN
        kBloatPlanckOpen = 1
        write(kStdWarn,*) 'Opened following file for bloated planck output :'
        write(kStdWarn,*) caBloatFile
    END IF

    WRITE(iIOUN1) kProfLayer
    WRITE(iIOUN1) rFrLow,rFrHigh
    WRITE(iIOUN1) iFileIDLo,iFileIDHi     !!user should expect x5 of these!
    WRITE(iIOUN1) kaFrStep(iTag)          !!regular freq step size
    WRITE(iIOUN1) kaBlSize(iTag)          !!regular 10k point freq block size
    WRITE(iIOUN1) iTotalStuff             !!regular number of outputs
!!!!!now do the fine stuff
    WRITE(iIOUN1) iType                   !!tells reader if reg +1 or plck -1
    WRITE(iIOUN1) iNumLayers              !!tells #of layers in atm 1
    WRITE(iIOUN1) kBoxCarUse              !!number of boxcar pts
    WRITE(iIOUN1) sngl(kaFineFrStep(iTag))!!fine freq step size
    WRITE(iIOUN1) i1,i2,iNLTEChunks       !!which chunks we are dealing with
! and how many to expect
    WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks
    WRITE(iIOUN1) kaBlSize(iTag)/kBoxCarUse
! fine 10000 point freq block size
    WRITE(iIOUN1) iTotalStuff*kBoxCarUse  !!fine number of outputs
          
    RETURN
    end SUBROUTINE OpenOutputBloatFile

!************************************************************************
! this subroutine dumps out gas optical depths at high resolution
    SUBROUTINE out_bloat(raFreq,rFreqStart,rFreqEnd, &
    iType,daFreqBloat, &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter, &
    caPlanckBloatFile,caOutBloatFile, &
    iFileID, &
    iaPath,iNp,iaOp)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! iWHich         = +1 for optical depths, -1 for planck modifiers
! raaGasAbs  = single gas abs coeffs
! iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to
! iaPath     = list of the paths corresponding to the current gas
! iNp        = total number of paths to be output
! iaOp       = list of the paths to be output
    REAL :: raFreq(kMaxPts),rFreqStart,rFreqEnd
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    INTEGER :: iPrinter,iFileID,iType
    INTEGER :: iNp,iaOp(kPathsOut),iaPath(kProfLayer)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile

! local variables
    INTEGER :: iIntOffset,iStart
    INTEGER :: iInt,iDp,iPath,iLay,iIOUN,iLoop,iL
    REAL :: raL2S(kMaxPts),raF(kMaxPts)
    CHARACTER(80) :: caOut

    IF (iType > 0) THEN
        iIOUN = kBloatNLTEOut
        caOut = caOutBloatFile
    ELSE
        write(kStdErr,*) 'Should be calling out_bloat_planck'
        CALL DoStop
    END IF

! thius certainly works for gas paths; might be wierd for mixed paths
    iStart = iDiv(iaPath(1),kProfLayer)

! write spectra to unformatted file
! if iPrinter=1 then have to check for valid paths
    DO iLay = 1,kProfLayer
    ! check to see if this path should be output
        iPath = iStart*kProfLayer + iLay
        iDp = DoOutputLayer(iPath,iNp,iaOp)

        IF ((iDp > 0) .AND. (iType > 0) .AND. (kBloatOutOpen > 0)) THEN
            IF (iPrinter == 1) THEN
                write(kStdWarn,*)'output high res GAS opt depths at layer = ',iPath
            ELSEIF (iPrinter == 2) THEN
                write(kStdWarn,*)'output high res MIXPATH depths at layer = ',iPath
            ELSEIF (iPrinter == 3) THEN
                write(kStdWarn,*)'output high res RADIANCES at layer = ',iPath
            END IF

            DO iLoop = 1,kBoxCarUse

                IF (kLayer2Sp == -1) THEN
                    DO iInt=1,kMaxPts
                        iIntOffSet = iInt + (iLoop-1)*kMaxPts
                        raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                        raL2S(iInt) = sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
                    END DO
                    CALL  wrtout(iIOUN,caOut,raF,raL2S)
                ELSE IF (kLayer2Sp == -2) THEN
                    DO iInt=1,kMaxPts
                        iIntOffSet = iInt + (iLoop-1)*kMaxPts
                        raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                        raL2S(iInt)=sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
                        raL2S(iInt) = exp(-raL2S(iInt))
                    END DO
                    CALL wrtout(iIOUN,caOut,raF,raL2S)
                ELSE IF (kLayer2Sp == 1) THEN
                    DO iInt=1,kMaxPts
                        raL2S(iInt) = 0.0
                    END DO
                    DO iL = iLay,kProfLayer
                        DO iInt=1,kMaxPts
                            iIntOffSet = iInt + (iLoop-1)*kMaxPts
                            raL2S(iInt) = raL2s(iInt) + &
                            sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
                        END DO
                    END DO
                    DO iInt=1,kMaxPts
                        iIntOffSet = iInt + (iLoop-1)*kMaxPts
                        raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                    END DO
                    CALL wrtout(iIOUN,caOut,raF,raL2S)
                ELSE IF (kLayer2Sp == 2) THEN
                    DO iInt=1,kMaxPts
                        raL2S(iInt) = 0.0
                    END DO
                    DO iL = iLay,kProfLayer
                        DO iInt=1,kMaxPts
                            iIntOffSet = iInt + (iLoop-1)*kMaxPts
                            raL2S(iInt) = raL2s(iInt) + &
                            sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
                        END DO
                    END DO
                    DO iInt=1,kMaxPts
                        iIntOffSet = iInt + (iLoop-1)*kMaxPts
                        raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                        raL2S(iInt) = exp(-raL2S(iInt))
                    END DO
                    CALL wrtout(iIOUN,caOut,raF,raL2S)

                END IF       !IF (iKlayerSp = -2,-1,1,2)
            END DO         !DO iLoop = 1,kBoxCarUse
        END IF           !IF ((iDp > 0) .AND. (iType .EQ -1,1)) THEN
    END DO             !DO iLay = 1,kProfLayer

    RETURN
    end SUBROUTINE out_bloat

!************************************************************************
! this subroutine dumps out gas optical depths at high resolution
    SUBROUTINE out_bloat_planck(raFreq,rFreqStart,rFreqEnd, &
    iType,daFreqBloat, &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter, &
    caPlanckBloatFile,caOutBloatFile, &
    iFileID, &
    iAtm,iNumLayers,iaaRadLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! iWHich         = +1 for optical depths, -1 for planck modifiers
! raaGasAbs  = single gas abs coeffs
! iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to
! iaPath     = list of the paths corresponding to the current gas
! iNp        = total number of paths to be output
! iaOp       = list of the paths to be output
    REAL :: raFreq(kMaxPts),rFreqStart,rFreqEnd
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    INTEGER :: iPrinter,iFileID,iType
    INTEGER :: iNp,iaOp(kPathsOut),iaPath(kProfLayer)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile
    INTEGER :: iAtm,iNumLayers,iaaRadLayer(kMaxAtm,kProfLayer)

! local variables
    INTEGER :: iIntOffset,iaRadLayer(kProfLayer)
    INTEGER :: iInt,iDp,iPath,iLay,iIOUN,iLoop,iL
    REAL :: raL2S(kMaxPts),raF(kMaxPts)
    CHARACTER(80) :: caOut

    IF (iType > 0) THEN
        write(kStdErr,*) 'Should be calling out_bloat'
        CALL DoStop
    ELSE
        iIOUN = kBloatNLTEPlanck
        caOut = caPlanckBloatFile
    END IF
     
    DO iInt = 1,iNumLayers
        iaRadLayer(iInt) = iaaRadLayer(iAtm,iInt)
    END DO

! write spectra to unformatted file
! if iPrinter=1 then have to check for valid paths
    DO iLay = 1,iNumLayers
    ! check to see if this path should be output
        iPath = iaRadLayer(iLay)
        IF ((iType < 0) .AND. (kBloatPlanckOpen > 0)) THEN
            write(kStdWarn,*)'output high res planck coeffs at layer = ',iPath
            DO iLoop = 1,kBoxCarUse
                DO iInt=1,kMaxPts
                    iIntOffSet = iInt + (iLoop-1)*kMaxPts
                    raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                    raL2S(iInt) = sngl(daaPlanckCoeffBloat(iIntOffset,iPath))
                END DO
                CALL wrtout(iIOUN,caOut,raF,raL2S)
            END DO         !DO iLoop = 1,kBoxCarUse

        END IF           !IF ((iDp > 0) .AND. (iType .EQ -1,1)) THEN
    END DO             !DO iLay = 1,kProfLayer

    RETURN
    end SUBROUTINE out_bloat_planck

!************************************************************************
!************************************************************************
!************************************************************************
END MODULE kbloat
