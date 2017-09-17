! Copyright 2002
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:39
 
! University of Maryland Baltimore County
! All Rights Reserved

! bloated output routines for high res NLTE
!************************************************************************
! this subroutine adds on current LTE gas optical depths to cumulative mixed
! path optical depth
! almost the same as SUBROUTINE Accumulate(raaSum,raaGas,raaMix,iGas,iIpmix)

SUBROUTINE AccumulateForBloat(raaSum,raaGas,raaMix,iGas,iIpmix)


REAL, INTENT(OUT)                        :: raaSum(kMaxPts,kProfLayer)
REAL, INTENT(IN)                         :: raaGas(kMaxPts,kProfLayer)
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN OUT)                  :: iIpmix
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaSum     = cumulative spectra associated with the mixed paths
! raaGas     = current gas absorption spectra
! raaMix     = mixing table info from *MIXFIL
! iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
! iIpmix     = which of the mixed paths are being considered





INTEGER :: iFreq,iL,MP2Lay
REAL :: rL

IF (iIPMIX > kProfLayer) THEN
  WRITE(kStdErr,*) 'AccumulateForBloat assumes 1..kProfLayer mix paths'
  CALL DoStop
END IF

! find out which of the 100 layers is associated with this mixed path
iL=MP2Lay(iIpmix)

! find the weight
rL=raaMix(iIpmix,iGas)

! add on contribution of the iGas th gas to the iIpmix th row of raaSum
DO iFreq=1,kMaxPts
  raaSum(iFreq,iIpmix)=raaSum(iFreq,iIpmix)+rL*raaGas(iFreq,iL)
END DO

RETURN
END SUBROUTINE AccumulateForBloat

!************************************************************************
! this adds together stuff so as to output Mixed paths at high resolution

SUBROUTINE SumBloatMP(  &
    daFreqBloat,raFreq,raaCo2_LTE,raaRestOfLTEGases,iNLTEStart,  &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat)


NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaCo2_LTE(kMaxPts,kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaRestOfL
INTEGER, INTENT(IN)                      :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daFreqBloat(kBloatPts)
REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)




! input/output params
DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
DOUBLE PRECISION :: daFreqT(kMaxPts),daTemp(kMaxPts)
DOUBLE PRECISION :: dY,daY(kBloatPts)
INTEGER :: iFr,iL

DO iFr = 1,kMaxPts
  daFreqT(iFr) = raFreq(iFr)*1.0D0
END DO

! in SetPlanckCoeffBloat we have already done, from iNLTEStart to kProflayer
!       daaNLTEGasAbCoeffBloat(iFr,iL)+blah_from_raaRestOfLTEGases(iFr,iL)
! so no need to redo it here
! however we ned to do it from 1 t- NLTEStart-1

DO iL = 1,iNLTEStart-1
  DO iFr = 1,kMaxPts
!! we have already stored daaNLTEGasAbCoeffBloat= NLTE_CO2 + LTE_CO2
!! so we need to add on LTE_others to get total mixed path
    daTemp(iFr)  = raaRestOfLTEGases(iFr,iL)*1.0D0
  END DO
  CALL dlinear_smart(daFreqT,daTemp, kMaxPts,daFreqBloat,daY, kBloatPts)
  DO iFr = 1,kBloatPts
    daaSumNLTEGasAbCoeffBloat(iFr,iL) =  &
        MAX(daaNLTEGasAbCoeffBloat(iFr,iL)+daY(iFr),0.0D0)
  END DO
END DO

DO iL = iNLTEStart,kProfLayer
  DO iFr = 1,kBloatPts
    daaSumNLTEGasAbCoeffBloat(iFr,iL) =  &
        MAX(daaNLTEGasAbCoeffBloat(iFr,iL),0.0D0)
  END DO
END DO

RETURN
END SUBROUTINE SumBloatMP

!************************************************************************
! this subroutine changes double --> real for the KCARTA atmosphere
! this computes the cumulative multiplication modification to Planck
! this is for the usual 100 kCARTA database layers
! however, this is at the bloated HIGH resolution!!!!!!!

SUBROUTINE SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer,  &
    raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff, raaSumAbCoeff,raaPlanckCoeff,  &
    raaRestOfLTEGases,raaCO2_LTE, daFreqBloat,daaSumNLTEGasAbCoeffBloat,  &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)


INTEGER, INTENT(IN)                      :: iNLTEStart
INTEGER, INTENT(IN OUT)                  :: iAtm
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: raaRestOfL
REAL, INTENT(IN OUT)                     :: raaCO2_LTE(kMaxPts,kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params

DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
!!cumulative nonlte gas optical depths
DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer) !!planck mod, so far
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)             !!mixed path abscoeff
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)          !!rest of LTE gases
REAL :: !!LTE part of CO2

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)             !!set of mixed paths

DOUBLE PRECISION :: daFreqBloat(kBloatPts)
!!for NLTE gas, this contains cumulative optical depths = LTE + NLTE
DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
! input/output vars
!!at input, only contains NLTE opt depths for CO2
!!at output, contains (LTE+NLTE) opt depths for CO2, plus other LTE gases
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)

! local variables
DOUBLE PRECISION :: daFreqT(kMaxPts)
DOUBLE PRECISION :: dY, daY(kBloatPts), daTemp(kMaxPts)
DOUBLE PRECISION :: dY1,daY1(kBloatPts),daTemp1(kMaxPts)
INTEGER :: iFr,iL,iA,iAtmTemp,iM,iStart

DO iFr = 1,kMaxPts
  daFreqT(iFr) = raFreq(iFr)*1.0D0
END DO

iAtmTemp = 1
IF (iAtmTemp /= iatm) THEN
  WRITE(kStdErr,*) 'oh oh iAtmTemp .NE. iatm ',iAtmTemp,iatm
  WRITE(kStdErr,*) 'kCARTA assumes the NLTE is for and only for atm #1'
  CALL DoStop
END IF
iA = 1    !!!assume current atmospher uses mixed path layers 1-100
iM = iaaRadLayer(iAtmTemp,1)

!find which set of Mixed Paths current atm uses eg 1-100, 101-200 etc
DO iL = 1,kProfLayer
  IF (iM <= kProfLayer*iL) THEN
    iA = iL
    GO TO 10
  END IF
END DO
10   CONTINUE

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
    daTemp(iFr)  = raaRestOfLTEGases(iFr,iM)*1.0D0 + raaCo2_LTE(iFr,iM)*1.0D0
    daTemp1(iFr) = raaRestOfLTEGases(iFr,iM)*1.0D0
  END DO
!!!since abs(LTE) << abs(CO2_NLTE), can get away with linear interp
  CALL dlinear_smart(daFreqT,daTemp, kMaxPts,daFreqBloat,daY, kBloatPts)
  CALL dlinear_smart(daFreqT,daTemp1,kMaxPts,daFreqBloat,daY1,kBloatPts)
  DO iFr = 1,kBloatPts
!!dY is thus the LTE component of numerator
    dY = daY(iFr)
!!update this var so it gives the TOTAL optical depth at the layer
!! or Sum(NLTE_CO2 + LTE_CO2 + LTE_others)
!! we have already stored daaNLTEGasAbCoeffBloat = NLTE_CO2 + LTE_CO2
!! so we need to add on LTE_others
    dY = MAX(daaNLTEGasAbCoeffBloat(iFr,iL)+daY1(iFr),0.0D0)
    daaNLTEGasAbCoeffBloat(iFr,iL) = dY + dDeltaNLTE
!!for Planck, add on LTE component
    dY = MAX(daY(iFr)  + daaPlanckCoeffBloat(iFr,iL),0.0D0)
    dY = MAX(daY1(iFr) + daaPlanckCoeffBloat(iFr,iL),0.0D0)
    dY = dY + dDeltaNLTE
    daaPlanckCoeffBloat(iFr,iL) = dY/daaNLTEGasAbCoeffBloat(iFr,iL)
  END DO
END DO

1010 FORMAT(I6,' ',3(F20.7,' '))
1111 FORMAT(I6,' ',3(F20.7,' '))

RETURN
END SUBROUTINE SetPlanckCoeffBloat

!************************************************************************
! this subroutine simply bloats stuff, on a 1-1 correspondance using i1,i2

SUBROUTINE AccumulateBloat(daHighX,daHighY,iHigh,daLowX,daLowY,iLow,  &
    iWideMeshLoop,i1,i2,iType,daBloat)


DOUBLE PRECISION, INTENT(IN OUT)         :: daHighX(kMaxPtsBox)
DOUBLE PRECISION, INTENT(IN)             :: daHighY(kMaxPtsBox)
INTEGER, INTENT(IN)                      :: iHigh
DOUBLE PRECISION, INTENT(IN OUT)         :: daLowX(kMaxPtsBox)
DOUBLE PRECISION, INTENT(IN OUT)         :: daLowY(kMaxPtsBox)
INTEGER, INTENT(IN OUT)                  :: iLow
NO TYPE, INTENT(IN OUT)                  :: iWideMeshL
INTEGER, INTENT(IN)                      :: i1
INTEGER, INTENT(IN)                      :: i2
INTEGER, INTENT(IN OUT)                  :: iType
DOUBLE PRECISION, INTENT(OUT)            :: daBloat(kBloatPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
! iType = +1 ==> use daHighY at kMaxPtsBox resolution
!                  use the iHigh points from this vector
! iType = -1 ==> use daLowY  at kMaXpts    resolution
!                  use the iCount points from this vector


INTEGER :: iWideMeshLoop
! output params


! local vars
INTEGER :: iI,iJ,iK,iS,iE,iOffset
DOUBLE PRECISION :: daTemp(kBloatPts)

IF (iType == 1) THEN
  iOffSet = (iWideMeshLoop-1)*2000
  DO iJ = 1,iHigh-1
    daBloat(iJ+iOffSet) = daBloat(iJ+iOffSet) + daHighY(iJ)
  END DO
ELSE IF (iType == -1) THEN
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
END SUBROUTINE AccumulateBloat

!************************************************************************
! this is the driver
! to bloat up the weak LTE lines (and the stronger LTE lines low down
! in the atm) to higher resolution

! this subroutine bloats up the optical depths, at double precision
! can bloat up all layers, or just select layers
! usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
! at the higher altitudes

SUBROUTINE BloatCoeffsDriver(iTag,iActualTag,rFileStartFr,raFreq,daaGasAbCoeff,  &
    raaRestOfLTEGases,raaCO2_LTE,  &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,  &
    daFreqBloat,iNLTEStart,rNLTEStrength,  &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,  &
    iGas,iGasID,iL,iU,iUseWeakBackGnd, raRAmt,raRTemp,raRPress,raRPartPress,  &
    pProf,iProfileLayers, raPAmt,raPTemp,raPPress,raPPartPress,iSplineType)


INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
NO TYPE, INTENT(IN OUT)                  :: rFileStart
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: daaGasAbCo
NO TYPE, INTENT(IN OUT)                  :: raaRestOfL
REAL, INTENT(IN OUT)                     :: raaCO2_LTE(kMaxPts,kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rNLTEStren
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN OUT)                  :: iU
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
REAL, INTENT(IN OUT)                     :: raRAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: raPAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPPartPre
NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
REAL :: !!kCARTA res
DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
!!at kCARTA res
DOUBLE PRECISION :: daFreqBloat(kBloatPts)            !!high res pts
REAL :: rNLTEStrength                                 !!mixed path strength


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

REAL :: raPPartPress(kProfLayer)

REAL :: raRPartPress(kProfLayer)
INTEGER :: iProfileLayers
REAL :: raVTemp(kProfLayer)
INTEGER :: iCount, iErr,iRefLayer,iSplineType
INTEGER :: iUseWeakBackGnd

! output params
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)


! local variables
INTEGER :: iLL

WRITE(kStdWarn,*) 'Using std 0.0025 -> 0.0005 code for "highres" bloat'

CALL BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff,  &
    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd,  &
    raPAmt,raPTemp,raPPress,raPPartPress,  &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,  &
    daFreqBloat,iNLTEStart,rNLTEStrength,  &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

RETURN
END SUBROUTINE BloatCoeffsDriver

!************************************************************************
! this just bloats up the weak LTE lines (and the stronger LTE lines low down
! in the atm) to higher resolution
! this subroutine bloats up the optical depths, at double precision
! can bloat up all layers, or just select layers
! usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
! at the higher altitudes

SUBROUTINE BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff,  &
    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd,  &
    raPAmt,raPTemp,raPPress,raPPartPress,  &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,  &
    daFreqBloat,iNLTEStart,rNLTEStrength,  &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)


INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: rFileStart
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: daaGasAbCo
NO TYPE, INTENT(IN OUT)                  :: raaRestOfL
REAL, INTENT(OUT)                        :: raaCO2_LTE(kMaxPts,kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
REAL, INTENT(IN OUT)                     :: raPAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPPartPre
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
INTEGER, INTENT(IN)                      :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rNLTEStren
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
REAL :: !!kCARTA res
DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
!!at kCARTA res
DOUBLE PRECISION :: daFreqBloat(kBloatPts)            !!high res pts
REAL :: rNLTEStrength                                 !!mixed path strength


REAL :: rFileStartFr                             !!start point of chunk
INTEGER :: iUseWeakBackGnd                          !!go thru with this?
REAL :: pProf(kProfLayer)
REAL :: raPPartPress(kProfLayer)

! output params
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)


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
  WRITE(kStdErr,*) 'iDefault,iLinearOrSpline = ',iDefault,iLinearOrSpline
END IF

DO iFr = 1,kMaxPts
  daFreqT(iFr) = raFreq(iFr)*1.0D0
END DO

!!! at this point, the code has uncompressed the database, so the lower
!!!   layers contain the lower trop LTE optical depths
!!! in addition the code has replaced the (weak+strong) upper atm LTE
!!!    optical depths with the (weak) upper trop LTE depths
!!! So save the CO2 LTE stuff
DO iL = 1,kProfLayer
  DO iFr = 1,kMaxPts
    raaCO2_LTE(iFr,iL) = SNGL(daaNLTEGasAbCoeff(iFr,iL))
  END DO
END DO

! <---!!!! --------------------------->
!!! update stuff with the weak LTE lines that have been computed in
!!! by calling LTE_spectra. This is for upper atmosphere
!!set coeffs of upper layers
iS = iNLTEStart
iE = kProfLayer
WRITE (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iS,' up'
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
    dK = MAX(daY(iFr),0.0D0)
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
    daaPlanckCoeffBloat(iFr,iL) = MAX(daY(iFr),0.0D0)
  END DO
END DO

! <---!!!! --------------------------->
!!! update stuff with all LTE lines that have been computed from
!!! uncompressing the kCompressed Database. This is for lower  atmosphere
!!set coeffs of lower layers
iS = 1
iE = iNLTEStart-1
WRITE (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iE,' down'
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
    dK = MAX(daY(iFr),0.0D0)
    daaNLTEGasAbCoeffBloat(iFr,iL)    = dK
    daaSumNLTEGasAbCoeffBloat(iFr,iL) = dK
!          daaPlanckCoeffBloat(iFr,iL)       = 1.0d0
  END DO
END DO

RETURN
END SUBROUTINE BloatCoeffsStandard

!************************************************************************
! this subroutine bloats up the LTE coeffs from the upper atmosphere database

SUBROUTINE bloatUAstuff( daaUpperNLTEGasAbCoeff,daaUpperPlanckCoeff,iUpper,  &
    raFreq,daFreqBloat, daaUpperPlanckCoeffBloat,daaUpperNLTEGasAbCoeffBloat,  &
    daaUpperSumNLTEGasAbCoeffBloat)


NO TYPE, INTENT(IN OUT)                  :: daaUpperNL
NO TYPE, INTENT(IN OUT)                  :: daaUpperPl
INTEGER, INTENT(IN)                      :: iUpper
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaUpperPl
NO TYPE, INTENT(IN OUT)                  :: daaUpperNL
NO TYPE, INTENT(IN OUT)                  :: daaUpperSu
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params


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
  WRITE(kStdErr,*) 'iDefault,iLinearOrSpline = ',iDefault,iLinearOrSpline
END IF

iS = 1
iE = iUpper

DO iFr = 1,kMaxPts
  daFreqT(iFr) = raFreq(iFr)*1.0D0
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
    dK = MAX(daY(iFr),0.0D0)
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
    dK = MAX(daY(iFr),0.0D0)
    daaUpperPlanckCoeffBloat(iFr,iL) = dK
  END DO
END DO

RETURN
END SUBROUTINE bloatUAstuff

!************************************************************************
! this subroutine reads in the precomputed LBL weakbackgnd data
! as you can see this is hardcoded for
!   /home/sergio/SPECTRA/IPFILES/layop_md_va0_usethis_co2'

SUBROUTINE ReadFixedDataWeakBack_Bloat(iGasID,iStart,iTag,rFileStartFr,  &
    raPAmt,raPTemp,raPPress,raPPartPress, rNLTEStrength,daaWeakOptDepthBloat)


INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN OUT)                  :: iStart
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: rFileStart
REAL, INTENT(IN OUT)                     :: raPAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPPartPre
NO TYPE, INTENT(IN OUT)                  :: rNLTEStren
NO TYPE, INTENT(IN OUT)                  :: daaWeakOpt
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL :: pProf(kProfLayer)
REAL :: raPPartPress(kProfLayer)

REAL :: rFileSTartFr
DOUBLE PRECISION :: daaWeakOptDepthBloat(kBloatPts,kProfLayer)
REAL :: rNLTEStrength

! local variables
INTEGER :: iIoun,I,J,iErr, iL, iFr
INTEGER :: IDGAS, NPTS, NLAY
DOUBLE PRECISION :: SFREQ, FSTEP
CHARACTER (LEN=80) :: FNAM
CHARACTER (LEN=4) :: ca4
DOUBLEPRECISION df,sf,dK
INTEGER :: iDefault
REAL :: raFreqBloat(kBloatPts)

df = kaFineFrStep(iTag)
sf = 1.0D0*nint(rFileStartFr)
DO iFr = 1,kBloatPts
  raFreqBloat(iFr) = sf + (iFr-1)*df
END DO

FNAM = '/carrot/s1/sergio/AIRSCO2/BACKGND_vt_mdH96/hehe_highres_weaklte_'
I = 80
10   CONTINUE
IF (FNAM(I:I) == ' ') THEN
  I = I - 1
  GO TO 10
END IF
IF ((nint(rFileStartFr) < 1000) .OR. (nint(rFileStartFr) > 9999)) THEN
  WRITE(kStdErr,*) 'Ooooooerr cannot accept freqs outside 1000 ... 9999'
  CALL DoStop
END IF

WRITE(ca4,40) nint(rFileStartFr)
40   FORMAT(I4)
FNAM(I+1:I+4) = ca4(1:4)
FNAM(I+5:I+8) = '.dat'

WRITE(kStdWarn,*) 'In ReadFixedDataWeakBack_Bloat, open file for GasID'
WRITE(kStdWarn,*) iGasID
WRITE(kStdWarn,*) FNAM

iIOUN = kCompUnit
OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', IOSTAT=IERR)
IF (IERR /= 0) THEN
  WRITE(kStdErr,*) 'In subroutine ReadFixedDataWeakBack_Bloat'
  WRITE(kStdErr,1010) IERR, FNAM
  1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
  CALL DoSTOP
END IF
kCompUnitOpen=1

!     Read in the header
READ(iIOUN) IDGAS, NPTS, NLAY
READ(iIOUN) SFREQ, FSTEP

IF (IDGAS /= iGasID) THEN
  WRITE(kStdErr,*) 'Wanted gas ',iGasID,'but new data is for gas ',IDGAS
  CALL DoStop
END IF

IF (NLAY /= kProfLayer) THEN
  WRITE(kStdErr,*) ' new data has ',NLAY,'layers instead of kProfLayer'
  CALL DoStop
END IF

IF (NPTS /= kBloatPts) THEN
  WRITE(kStdErr,*) ' new data has ',NPTS,' freq pts instead of kBloatPts'
  CALL DoStop
END IF

IF (ABS(df-FSTEP) > 1E-4) THEN
  WRITE(kStdErr,*) ' new data has ',FSTEP,'freq spacing instead of ',df
  CALL DoStop
END IF

IF (ABS(sf-SFREQ) > 1E-4) THEN
  WRITE(kStdErr,*) ' new data starts at ',SFREQ,' instead of ',sf
  CALL DoStop
END IF

! passed all tests, so read data
!     Read in the optical depths
DO I=1,kProfLayer
  READ(iIOUN) (daaWeakOptDepthBloat(J,I),J=1,kBloatPts)
  ENDDO
    
    CLOSE(iIOUN)
    kCompUnitOpen=-1
    
    WRITE (kStdWarn,*) 'Read in HIGH RES WEAK BKGND spectra for gasID ',IDGAS
    
    RETURN
  END SUBROUTINE ReadFixedDataWeakBack_Bloat
  
!************************************************************************
! this opens the bloated file four output
  
  SUBROUTINE OpenOutputBloatFile(  &
      iType,iNumLayers,caBloatFile,rFrLow,rFrHigh,  &
      iFileIDLo,iFileIDHi,iTag,iTotalStuff,  &
      iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
  
  
  INTEGER, INTENT(OUT)                     :: iType
  INTEGER, INTENT(IN)                      :: iNumLayers
  NO TYPE, INTENT(IN OUT)                  :: caBloatFil
  REAL, INTENT(IN OUT)                     :: rFrLow
  REAL, INTENT(IN OUT)                     :: rFrHigh
  INTEGER, INTENT(OUT)                     :: iFileIDLo
  INTEGER, INTENT(OUT)                     :: iFileIDHi
  INTEGER, INTENT(OUT)                     :: iTag
  NO TYPE, INTENT(IN OUT)                  :: iTotalStuf
  NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
  NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
  NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  
  CHARACTER (LEN=80) :: caBloatFile
  
  INTEGER :: iTotalStuff
  
  
  
  INTEGER :: iaNLTEChunks(kGasStore)
  INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
  
! local vars
  INTEGER :: iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor
  REAL :: r1,r2,r3
  
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
  ELSE IF (iType < 0) THEN
    iIoun1 = kBloatNLTEPlanck
  ELSE
    WRITE(kStdErr,*) 'Unknown print option for bloated files'
    CALL DoStop
  END IF
  
! open unformatted file as a fresh file to be written to
  OPEN(UNIT=iIoun1,FILE=caBloatFile,STATUS='NEW',  &
      FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,304) iFileErr,iIOUN1,caBloatFile
    WRITE(kStdErr,*)'make sure the file does not exist!'
    CALL DoSTOP
  END IF
  
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*) 'Printing BLOAT FILE summary stats at open : '
  WRITE(kStdWarn,*) caBloatFile
  WRITE(kStdWarn,*) 'Regular freqs etc .... '
  WRITE(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
  WRITE(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
  WRITE(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
  WRITE(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff
  WRITE(kStdWarn,*) 'Bloated  freqs etc .... '
  WRITE(kStdWarn,*) '  iType = = ',iType
  WRITE(kStdWarn,*) '  iNumlayers = ',iNumLayers
  WRITE(kStdWarn,*) '  kBoxCarUse = ',kBoxCarUse
  WRITE(kStdWarn,*) '  df_fine = ',kaFineFrStep(iTag)
  WRITE(kStdWarn,*) 'start,stop high res integers : ',i1,i2
  WRITE(kStdWarn,*) 'number of highres chunks = ',iNLTEChunks
  WRITE(kStdWarn,*) 'start,stop high res doubles : ',r1,r2
  WRITE(kStdWarn,*) ' '
  
  304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file :  &
      ',/,A80)
  
  IF (iType > 0) THEN
    kBloatOutOpen = 1
    WRITE(kStdWarn,*) 'Opened following file for bloated general output :'
    WRITE(kStdWarn,*) caBloatFile
  ELSE IF (iType < 0) THEN
    kBloatPlanckOpen = 1
    WRITE(kStdWarn,*) 'Opened following file for bloated planck output :'
    WRITE(kStdWarn,*) caBloatFile
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
  WRITE(iIOUN1) SNGL(kaFineFrStep(iTag))!!fine freq step size
  WRITE(iIOUN1) i1,i2,iNLTEChunks       !!which chunks we are dealing with
!!and how many to expect
  WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks
  WRITE(iIOUN1) kaBlSize(iTag)/kBoxCarUse
!!fine 10000 point freq block size
  WRITE(iIOUN1) iTotalStuff*kBoxCarUse  !!fine number of outputs
  
  RETURN
END SUBROUTINE OpenOutputBloatFile

!************************************************************************
! this subroutine writes out the data header at each standard 10000 pt chunk

SUBROUTINE HeaderBloatFile(caOutFile,rFrLow,rFrHigh,daFreqBloat, iTag,iType)


CHARACTER (LEN=80), INTENT(IN)           :: caOutFile
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iType
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

DOUBLE PRECISION :: daFreqBloat(kBloatPts)



! local vars
INTEGER :: iIOUN1,iI,iJump,ifloor,i10000
CHARACTER (LEN=80) :: caOut
DOUBLE PRECISION :: daStart(kBoxCarUse),daEnd(kBoxCarUse)
DOUBLE PRECISION :: daStartBox(kBoxCarUse),daEndBox(kBoxCarUse)

iJump = ifloor(kBoxCarUse*1.0/2.0)

caOut = caOutFile

IF (iType == +2) THEN
  iIoun1 = kBloatNLTEOutUA
ELSE IF (iType == +1) THEN
  iIoun1 = kBloatNLTEOut
ELSE IF (iType == -1) THEN
  iIoun1 = kBloatNLTEPlanck
!      ELSEIF (iType .EQ. -2) THEN
!        iIoun1 = kBloatNLTEUAPlanck
ELSE
  WRITE(kStdErr,*) 'Unknown print option for bloated files'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'These are the chunks for high resolution file ',iIOUN1
DO iI = 1,kBoxCarUse
  i10000 = kMaxPts*(iI-1)
!!!these are the high resolution start/stop points for the chunk
  daStart(iI) = daFreqBloat(i10000 + 1)
  daEnd(iI)   = daFreqBloat(i10000 + kMaxPts)
!!!after boxcar avg, these are the kCARTA resolution start/stop points
  daStartBox(iI) = daFreqBloat(i10000 + 1 + iJump)
  daEndBox(iI)   = daFreqBloat(i10000 + kMaxPts - iJump)
  WRITE(kStdWarn,*) iI,daStart(iI),daEnd(iI),daStartBox(iI),daEndBox(iI)
END DO

WRITE(iIOUN1) kMaxPts,rFrLow,rFrHigh,kaFrStep(iTag)   !usual kCARTA stuff
WRITE(iIOUN1) kBloatPts,kaFineFrStep(iTag),kBoxCarUse !high res stuff
WRITE(iIOUN1) (daStart(iI),iI=1,kBoxCarUse)
WRITE(iIOUN1) (daEnd(iI),iI=1,kBoxCarUse)
WRITE(iIOUN1) (daStartBox(iI),iI=1,kBoxCarUse)
WRITE(iIOUN1) (daEndBox(iI),iI=1,kBoxCarUse)

RETURN
END SUBROUTINE HeaderBloatFile

!************************************************************************
! this subroutine just dumps out kBoxCarUse worth of chunks of data

SUBROUTINE wrtout_bloated_rad(iIOUN,caOutBloatFile, raFreqBloat,raIntenBloat)


INTEGER, INTENT(IN OUT)                  :: iIOUN
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
NO TYPE, INTENT(IN OUT)                  :: raFreqBloa
NO TYPE, INTENT(IN OUT)                  :: raIntenBlo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars
REAL :: raFreqBloat(kBloatPts),raIntenBloat(kBloatPts)

CHARACTER (LEN=80) :: caOutBloatFile

! local vars
INTEGER :: iLoop,iInt,iIntOffSet
REAL :: raF(kMaxPts),raInten(kMaxPts)

DO iLoop = 1,kBoxCarUse
  DO iInt=1,kMaxPts
    iIntOffSet    = iInt + (iLoop-1)*kMaxPts
    raF(iInt)     = raFreqBloat(iIntOffSet)
    raInten(iInt) = raIntenBloat(iIntOffset)
  END DO
  CALL wrtout(iIOUN,caOutBloatFile,raF,raInten)
END DO

RETURN
END SUBROUTINE wrtout_bloated_rad

!************************************************************************
! this subroutine dumps out gas optical depths at high resolution

SUBROUTINE out_bloat(raFreq,rFreqStart,rFreqEnd, iType,daFreqBloat,  &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,  &
    caPlanckBloatFile,caOutBloatFile, iFileID,  &
    iaPath,iNp,iaOp)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFreqStart
REAL, INTENT(IN OUT)                     :: rFreqEnd
INTEGER, INTENT(IN OUT)                  :: iType
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
INTEGER, INTENT(IN OUT)                  :: iPrinter
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iaPath(kProfLayer)
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

DOUBLE PRECISION :: daFreqBloat(kBloatPts)
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)


CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile

! local variables
INTEGER :: iIntOffset,iStart
INTEGER :: iInt,iDiv,iDp,iPath,iLay,DoOutputLayer,iIOUN,iLoop,iL
REAL :: raL2S(kMaxPts),raF(kMaxPts)
CHARACTER (LEN=80) :: caOut

IF (iType > 0) THEN
  iIOUN = kBloatNLTEOut
  caOut = caOutBloatFile
ELSE
  WRITE(kStdErr,*) 'Should be calling out_bloat_planck'
  CALL DoStop
END IF

! thius certainly works for gas paths; might be wierd for mixed paths
iStart=iDiv(iaPath(1),kProfLayer)

! write spectra to unformatted file
! if iPrinter=1 then have to check for valid paths
DO iLay = 1,kProfLayer
! check to see if this path should be output
  iPath = iStart*kProfLayer + iLay
  iDp = DoOutputLayer(iPath,iNp,iaOp)
  
  IF ((iDp > 0).AND.(iType > 0) .AND. (kBloatOutOpen > 0)) THEN
    IF (iPrinter == 1) THEN
      WRITE(kStdWarn,*)'output high res GAS opt depths at layer = ',iPath
    ELSE IF (iPrinter == 2) THEN
      WRITE(kStdWarn,*)'output high res MIXPATH depths at layer = ',iPath
    ELSE IF (iPrinter == 3) THEN
      WRITE(kStdWarn,*)'output high res RADIANCES at layer = ',iPath
    END IF
    
    DO iLoop = 1,kBoxCarUse
      
      IF (kLayer2Sp == -1) THEN
        DO iInt=1,kMaxPts
          iIntOffSet = iInt + (iLoop-1)*kMaxPts
          raF(iInt)   = SNGL(daFreqBloat(iIntOffSet))
          raL2S(iInt) = SNGL(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
        END DO
        CALL  wrtout(iIOUN,caOut,raF,raL2S)
      ELSE IF (kLayer2Sp == -2) THEN
        DO iInt=1,kMaxPts
          iIntOffSet = iInt + (iLoop-1)*kMaxPts
          raF(iInt)   = SNGL(daFreqBloat(iIntOffSet))
          raL2S(iInt)=SNGL(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
          raL2S(iInt) = EXP(-raL2S(iInt))
        END DO
        CALL wrtout(iIOUN,caOut,raF,raL2S)
      ELSE IF (kLayer2Sp == 1) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt) = 0.0
        END DO
        DO iL = iLay,kProfLayer
          DO iInt=1,kMaxPts
            iIntOffSet = iInt + (iLoop-1)*kMaxPts
            raL2S(iInt) = raL2s(iInt) +  &
                SNGL(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
          END DO
        END DO
        DO iInt=1,kMaxPts
          iIntOffSet = iInt + (iLoop-1)*kMaxPts
          raF(iInt)   = SNGL(daFreqBloat(iIntOffSet))
        END DO
        CALL wrtout(iIOUN,caOut,raF,raL2S)
      ELSE IF (kLayer2Sp == 2) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt) = 0.0
        END DO
        DO iL = iLay,kProfLayer
          DO iInt=1,kMaxPts
            iIntOffSet = iInt + (iLoop-1)*kMaxPts
            raL2S(iInt) = raL2s(iInt) +  &
                SNGL(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
          END DO
        END DO
        DO iInt=1,kMaxPts
          iIntOffSet = iInt + (iLoop-1)*kMaxPts
          raF(iInt)   = SNGL(daFreqBloat(iIntOffSet))
          raL2S(iInt) = EXP(-raL2S(iInt))
        END DO
        CALL wrtout(iIOUN,caOut,raF,raL2S)
        
      END IF       !IF (iKlayerSp = -2,-1,1,2)
    END DO         !DO iLoop = 1,kBoxCarUse
  END IF           !IF ((iDp .GT. 0) .AND. (iType .EQ -1,1)) THEN
END DO             !DO iLay = 1,kProfLayer

RETURN
END SUBROUTINE out_bloat

!************************************************************************
! this subroutine dumps out gas optical depths at high resolution

SUBROUTINE out_bloat_planck(raFreq,rFreqStart,rFreqEnd, iType,daFreqBloat,  &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,  &
    caPlanckBloatFile,caOutBloatFile, iFileID,  &
    iAtm,iNumLayers,iaaRadLayer)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFreqStart
REAL, INTENT(IN OUT)                     :: rFreqEnd
INTEGER, INTENT(IN OUT)                  :: iType
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaNLTEGas
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
INTEGER, INTENT(IN OUT)                  :: iPrinter
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayers
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

DOUBLE PRECISION :: daFreqBloat(kBloatPts)
DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)

INTEGER :: iNp,iaOp(kPathsOut),iaPath(kProfLayer)
CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! local variables
INTEGER :: iIntOffset,iaRadLayer(kProfLayer)
INTEGER :: iInt,iDiv,iDp,iPath,iLay,DoOutputLayer,iIOUN,iLoop,iL
REAL :: raL2S(kMaxPts),raF(kMaxPts)
CHARACTER (LEN=80) :: caOut

IF (iType > 0) THEN
  WRITE(kStdErr,*) 'Should be calling out_bloat'
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
    WRITE(kStdWarn,*)'output high res planck coeffs at layer = ',iPath
    DO iLoop = 1,kBoxCarUse
      DO iInt=1,kMaxPts
        iIntOffSet = iInt + (iLoop-1)*kMaxPts
        raF(iInt)   = SNGL(daFreqBloat(iIntOffSet))
        raL2S(iInt) = SNGL(daaPlanckCoeffBloat(iIntOffset,iPath))
      END DO
      CALL wrtout(iIOUN,caOut,raF,raL2S)
    END DO         !DO iLoop = 1,kBoxCarUse
    
  END IF           !IF ((iDp .GT. 0) .AND. (iType .EQ -1,1)) THEN
END DO             !DO iLay = 1,kProfLayer

RETURN
END SUBROUTINE out_bloat_planck

!************************************************************************
! the main dog! radiative transfer! DOUBLE precision
! BIG Assumption, is that we only dump out radiances at TOP of layers
! compared to rad_main, where we could dump radiances in midlayer

SUBROUTINE radiances_bloat( raFreq,raVTemp,caOutBloatFile,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    daFreqBloat,daaSumNLTEOptDepthBloat,daaPlanckCoeffBloat,  &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,  &
    daaUpperNLTEGasAbCoeffBloat)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperSu
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
NO TYPE, INTENT(IN OUT)                  :: daaUpperPl
NO TYPE, INTENT(IN OUT)                  :: daaUpperSu
NO TYPE, INTENT(IN OUT)                  :: daaUpperNL
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iNLTEStart  = which layer NONLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutBloatFile  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)



INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iDumpAllUARads
CHARACTER (LEN=80) :: caOutBloatFile
! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! these are the bloated lower atmosphere matrices
DOUBLE PRECISION :: daaSumNLTEOptDepthBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daFreqBloat(kBloatPts)
DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
DOUBLE PRECISION :: daaLayTrans(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaEmission(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaLay2Sp(kBloatPts,kProfLayer)
REAL :: raSumLayEmission(kBloatPts),raSurfaceEmissionToSpace(kBloatPts)
DOUBLE PRECISION :: daSumLayEmission(kBloatPts),  &
    daSurfaceEmissionToSpace(kBloatPts)
REAL :: rDum1,rDum2,ttorad,rOmegaSun,rCO2
! to do the thermal,solar contribution
REAL :: raIntenBloat(kBloatPts)
DOUBLE PRECISION :: daIntenBloat(kBloatPts)
REAL :: raSurfaceBloat(kBloatPts),raEmissivityBloat(kBloatPts)
DOUBLE PRECISION :: daSurfaceBloat(kBloatPts),daEmissivityBloat(kBloatPts)
REAL :: raSunBloat(kBloatPts)
DOUBLE PRECISION :: daSunBloat(kBloatPts)
REAL :: raThermalBloat(kBloatPts),raSunReflBloat(kBloatPts)
DOUBLE PRECISION :: daThermalBloat(kBloatPts),daSunReflBloat(kBloatPts)
REAL :: rThermalRefl,rCos,rPlanck,rMPTemp
DOUBLE PRECISION :: dThermalRefl,DCOS,dPlanck,dMPTemp,dttorad
REAL :: raVT1(kMixFilRows),rXYZ
DOUBLE PRECISION :: daVT1(kMixFilRows),dXYZ
INTEGER :: iDoThermal,iDoSolar,MP2Lay

DOUBLE PRECISION :: dEmission, dTrans
REAL :: raOutFrac(kProfLayer),raFreqBloat(kBloatPts),InterpTemp
INTEGER :: iIOUN
REAL :: bt2rad,t2s
INTEGER :: iFr1

INTEGER :: i1,i2,iFloor,iDownWard,iSTopNormalRadTransfer

! set the direction of radiation travel
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = 1
  i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
  i2=iaaRadLayer(iAtm,iNumLayer)-1
  i2=iFloor(i2*1.0/kProfLayer)
  IF (rTSpace > 5.0) THEN
    WRITE(kStdErr,*) 'you want satellite to be downward looking'
    WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
    WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
    WRITE(kStdErr,*) 'Please retry'
    CALL DoSTOP
  END IF
ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
! radiation travelling downwards to instrument ==> sat looking up
! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = -1
  i1=iaaRadLayer(iAtm,1)-1
  i1=iFloor(i1*1.0/(1.0*kProfLayer))
  i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
  WRITE(kStdErr,*) 'Huh ? NLTE for UPLOOK instrument?????'
  CALL DoStop
END IF
WRITE(kStdWarn,*) ' ------------------> <----------------------'
WRITE(kStdWarn,*) 'Doing BLOATED nlte rad transfer ...'
WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUp=90, iLow=1 is acceptable
! eg iUp=140,iLow=90 is NOT acceptable
IF (i1 /= i2) THEN
  WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
  WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
  WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
  CALL DoSTOP
END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
IF (i1 > kProfLayer) THEN
  WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
  CALL DoSTOP
END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

! here we go .....

iIOUN = kBloatNLTEOut

rThermalRefl=1.0/kPi
dThermalRefl=1.0/kPi

DO iFr = 1,kBloatPts
  raFreqBloat(iFr) = SNGL(daFreqBloat(iFr))
END DO

! calculate cos(SatAngle)
rCos = COS(rSatAngle*kPi/180.0)
DCOS = DBLE(rCos)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop


! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) > kProfLayer) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Err, assume only 1 .. kProfLayer mixedpaths'
    WRITE(kStdErr,*) '  in the high res matrices'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
  daVT1(iFr) = DBLE(raVTemp(iFr))
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
daVT1(iL)=DBLE(raVT1(iL))
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
daVT1(iL)=DBLE(raVT1(iL))
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh=iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS = DBLE(rCos)
  DO iFr=1,kBloatPts
    dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)*DBLE(rFracBot)/DCOS
    daaLayTrans(iFr,iLay) = DEXP(-dXYZ)
    daaEmission(iFr,iLay) = 0.0D0
  END DO
END DO

DO iLay=2,iNumLayer-1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS = DBLE(rCos)
  DO iFr=1,kBloatPts
    dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)/DCOS
    daaLayTrans(iFr,iLay) = DEXP(-dXYZ)
    daaEmission(iFr,iLay) = 0.0D0
  END DO
END DO
DO iLay=iNumLayer,iNumLayer
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS = DBLE(rCos)
  DO iFr=1,kBloatPts
    dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)*DBLE(rFracTop)/DCOS
    daaLayTrans(iFr,iLay) = DEXP(-dXYZ)
    daaEmission(iFr,iLay) = 0.0D0
  END DO
END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NONLTE computations done!
! so usually only the usual LTE computations are done!!
IF (iNLTEStart > kProfLayer) THEN
  iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
  WRITE (kStdWarn,*) 'Normal rad transfer .... no NLTE'
  WRITE (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
ELSE
  iLay = 1
  987    CONTINUE
  iL=iaRadLayer(iLay)
  iLModKprofLayer = MOD(iL,kProfLayer)
  IF (iLModKprofLayer == 0) THEN
    iLModKprofLayer = kProfLayer
  END IF
  IF ((iLModKprofLayer < iNLTEStart).AND.(iLay < iNumLayer)) THEN
    iLay = iLay + 1
    GO TO 987
  END IF
  iSTopNormalRadTransfer = iLay
  WRITE (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
  WRITE (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
END IF

DO iLay = 1,iNumLayer
  iL      = iaRadLayer(iLay)
! first get the Mixed Path temperature for this radiating layer
  dMPTemp = raVT1(iL)
  iLModKprofLayer = MOD(iL,kProfLayer)
  IF (iLModKprofLayer == 0) THEN
    iLModKprofLayer = kProfLayer
    iL = kProfLayer
  END IF
  IF (iLModKprofLayer < iNLTEStart) THEN
!normal, no LTE emission stuff
    DO iFr=1,kBloatPts
      dPlanck = dttorad(daFreqBloat(iFr),dMPTemp)
      daaEmission(iFr,iLay) = (1.0-daaLayTrans(iFr,iLay))*dPlanck
    END DO
  ELSE IF (iLModKprofLayer >= iNLTEStart) THEN
!new; LTE emission stuff
    DO iFr=1,kBloatPts
      dPlanck = dttorad(daFreqBloat(iFr),dMPTemp)
      dPlanck = dPlanck*daaPlanckCoeffBloat(iFr,iL)
      daaEmission(iFr,iLay) = (1.0-daaLayTrans(iFr,iLay))*dPlanck
    END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL rspl(raFreq,raThermal,kMaxPts,raFreqBloat, raThermalBloat,kBloatPts)
ELSE
  DO iFr = 1,kBloatPts
    raThermalBloat(iFr) = 0.0
  END DO
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL rspl(raFreq,raSun,kMaxPts,raFreqBloat,raSunBloat,kBloatPts)
  CALL rspl(raFreq,raSunRefl,kMaxPts,raFreqBloat, raSunReflBloat,kBloatPts)
ELSE
  DO iFr = 1,kBloatPts
    raSunBloat(iFr) = 0.0
  END DO
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

CALL rlinear(raFreq,raSurface,kMaxPts, raFreqBloat,raSurfaceBloat,kBloatPts)
CALL rlinear(raFreq,raUseEmissivity,kMaxPts,  &
    raFreqBloat,raEmissivityBloat,kBloatPts)

DO iFr=1,kBloatPts
  raIntenBloat(iFr)=raSurfaceBloat(iFr)*raEmissivityBloat(iFr)+  &
      raThermalBloat(iFr)*(1.0-raEmissivityBloat(iFr))*rThermalRefl+  &
      raSunBloat(iFr)*raSunReflBloat(iFr)
  daIntenBloat(iFr) = DBLE(raIntenBloat(iFr))
END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL      = iaRadLayer(iLay)
  rCos    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS    = DBLE(rCos)
  rMPTemp = raVT1(iL)
  dMPTemp = DBLE(raVT1(iL))
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
! do the radiative transfer thru this bottom layer
  DO iFr=1,kBloatPts
    daIntenBloat(iFr) = daaEmission(iFr,iLay) +  &
        daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
    raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
  END DO
  
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
    CALL wrtout_bloated_rad(iIOUN, caOutBloatFile,raFreqBloat,raIntenBloat)
  ELSE IF (iDp >= 2) THEN
    WRITE(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
    WRITE(kStdErr,*) 'more than one radiance per layer'
    CALL DoStop
  END IF
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL      = iaRadLayer(iLay)
  rCos    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS    = DBLE(rCos)
  rMPTemp = raVT1(iL)
  dMPTemp = DBLE(raVT1(iL))
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
! now do the radiative transfer thru this complete layer
  DO iFr=1,kBloatPts
    daIntenBloat(iFr) = daaEmission(iFr,iLay) +  &
        daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
    raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
  END DO
  
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
    CALL wrtout_bloated_rad(iIOUN, caOutBloatFile,raFreqBloat,raIntenBloat)
  ELSE IF (iDp >= 2) THEN
    WRITE(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
    WRITE(kStdErr,*) 'more than one radiance per layer'
    CALL DoStop
  END IF
  
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
777  CONTINUE
DO iLay=iHigh,iHigh
  iL      = iaRadLayer(iLay)
  rCos    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DCOS    = DBLE(rCos)
  rMPTemp = raVT1(iL)
  dMPTemp = DBLE(raVT1(iL))
  
  IF (iUpper >= 1) THEN
!!! need to compute stuff at extra layers (100-200 km)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp >= 1) THEN
      
      WRITE(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer'
      WRITE(kStdWarn,*) 'This is the top of the usual AIRS atmosphere'
      WRITE(kStdWarn,*) '   you have iDoUpperATM > 0'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output stuff into the blah_UA file'
      WRITE(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of'
      WRITE(kStdWarn,*) 'stratosphere into both this and the UA file'
      
!do radiative transfer thru this layer
      DO iFr=1,kBloatPts
        daIntenBloat(iFr) =  &
            daaEmission(iFr,iLay)+daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
        raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
      END DO
      
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTransBloatDouble(raFreqBloat,  &
          daIntenBloat,daFreqBloat,raLayAngles(MP2Lay(iL)),  &
          iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,  &
          raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
!!! forget about interpolation thru the layers, just dump out the
!!! radiance at the top of stratosphere (120-200 km)
      
      WRITE(kStdWarn,*) 'outputting bloated stuff at TOTAL Complete TOA into'
      WRITE(kStdWarn,*) 'usual BLOATED binary file (iLay = ',iLay,')'
      
      DO iFr=1,kBloatPts
        raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
      END DO
      DO iFr=1,iDp
        CALL wrtout_bloated_rad(iIOUN,  &
            caOutBloatFile,raFreqBloat,raIntenBloat)
      END DO
    END IF
  END IF
  
  IF (iUpper < 1) THEN
!!! no need to compute stuff at extra layers (100-200 km)
!!! so just do usual stuff
!!! see if this mixed path layer is in the list iaOp to be output
!!! since we might have to do fractions!
    DO iFr = 1,kBloatPts
      daIntenBloat(iFr) = daaEmission(iFr,iLay) +  &
          daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
      raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
    END DO
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp == 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' rad layer'
      CALL wrtout_bloated_rad(iIOUN, caOutBloatFile,raFreqBloat,raIntenBloat)
    ELSE IF (iDp >= 2) THEN
      WRITE(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
      WRITE(kStdErr,*) 'more than one radiance per layer'
      CALL DoStop
    END IF
  END IF
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

RETURN
END SUBROUTINE radiances_bloat

!************************************************************************
! this subroutine very quickly does the radiative transfer (DOUBLE PREC)
! since the optical depths are soooooooooo small, use double precision

SUBROUTINE UpperAtmRadTransBloatDouble( raFreqBloat,  &
    daIntenBloat,daFreqBloat,rSatAngle,  &
    iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)


NO TYPE, INTENT(IN OUT)                  :: raFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daIntenBlo
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
REAL, INTENT(IN OUT)                     :: rSatAngle
INTEGER, INTENT(IN)                      :: iUpper
NO TYPE, INTENT(IN OUT)                  :: daaUpperPl
NO TYPE, INTENT(IN OUT)                  :: daaUpperSu
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
!   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!   upper atm abs coeff and planck coeff
DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
!   input wavevector and integer stating which layer to stop rad transfer at
DOUBLE PRECISION :: daFreqBloat(kBloatPts)
REAL :: raFreqBloat(kBloatPts)

! do we want to do upper atm NLTE computations?
INTEGER :: iDoUpperAtmNLTE
! do we dump all or some rads?
INTEGER :: iDumpAllUARads
! input/output pararameter
!   this contains the radiance incident at kCARTA TOA (0.005 mb)
!   it will finally contain the radiance exiting from TOP of UPPER ATM
DOUBLE PRECISION :: daIntenBloat(kBloatPts)

! local variables
INTEGER :: iFr,iL,iIOUN
DOUBLE PRECISION :: dEmission,dTrans,rMu,daInten0(kBloatPts)
REAL :: raIntenBloat(kBloatPts),ttorad
CHARACTER (LEN=80) :: caOutName

caOutName = 'DumDum'
iIOUN = kBloatNLTEOutUA

IF (iDoUpperAtmNLTE <= 0) THEN
  WRITE (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
  CALL DoStop
ELSE
  WRITE(kStdWarn,*) 'Doing UA (NLTE) radtransfer at highres 0.0005 cm-1 '
END IF

!!compute radiance intensity thru NEW uppermost layers of atm
DO iFr = 1,kBloatPts
  daInten0(iFr) = daIntenBloat(iFr)
  raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
END DO

iL = 0
IF (kNLTEOutUAOpen > 0) THEN
  WRITE(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
!!always dump out the 0.005 mb TOA radiance if the UA file is open
  CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
END IF

rMu = COS(rSatAngle*kPi/180.0)

DO iL = 1,iUpper-1
  
  DO iFr = 1,kBloatPts
    dTrans = daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)/(rMu*1.0D0)
    dTrans = EXP(-dTrans)
    dEmission = daaUpperPlanckCoeffBloat(iFr,iL) *  &
        DBLE(ttorad(raFreqBloat(iFr),raUpperTemp(iL))*1.0D0)* (1.0D0 - dTrans)
    daIntenBloat(iFr) = dEmission + daIntenBloat(iFr)*dTrans
    
    raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
  END DO
  
  IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
    WRITE(kStdWarn,*) 'dumping out UA highres rads at iL = ',iL
!!dump out the 0.000025 mb TOA radiance
    CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
  END IF
END DO

DO iL = iUpper,iUpper
  DO iFr = 1,kBloatPts
    dTrans = daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)/(rMu*1.0D0)
    dTrans = EXP(-dTrans)
    dEmission = daaUpperPlanckCoeffBloat(iFr,iL) *  &
        DBLE(ttorad(raFreqBloat(iFr),raUpperTemp(iL))*1.0D0)* (1.0D0 - dTrans)
    daIntenBloat(iFr) = dEmission + daIntenBloat(iFr)*dTrans
    raIntenBloat(iFr) = SNGL(daIntenBloat(iFr))
  END DO
  
  IF (kNLTEOutUAOpen > 0) THEN
!!always dump out 0.000025 mb TOA radiance if the UA file is open
    WRITE(kStdWarn,*) 'dumping out 0.000025 mb highres UA rads iL = ',iL
!!dump out the 0.000025 mb TOA radiance
    CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
  END IF
  
END DO

11   FORMAT(I3,' ',9(E15.8,' '))
3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

RETURN
END SUBROUTINE UpperAtmRadTransBloatDouble

!************************************************************************
!************************************************************************
!************************************************************************
