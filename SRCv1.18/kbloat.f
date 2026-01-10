c Copyright 2002  
c University of Maryland Baltimore County  
c All Rights Reserved 

c bloated output routines for high res NLTE
c************************************************************************
c this subroutine adds on current LTE gas optical depths to cumulative mixed
c path optical depth
c almost the same as SUBROUTINE Accumulate(raaSum,raaGas,raaMix,iGas,iIpmix) 
      SUBROUTINE AccumulateForBloat(raaSum,raaGas,raaMix,iGas,iIpmix) 

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c raaSum     = cumulative spectra associated with the mixed paths 
c raaGas     = current gas absorption spectra 
c raaMix     = mixing table info from *MIXFIL 
c iGas       = gas # iGas of iNumGases being added to the cumulative raaSum 
c iIpmix     = which of the mixed paths are being considered 
      REAL raaMix(kMixFilRows,kGasStore) 
      INTEGER iIpmix,iGas 
      REAL raaSum(kMaxPts,kProfLayer) 
      REAL raaGas(kMaxPts,kProfLayer) 
 
      INTEGER iFreq,iL,MP2Lay 
      REAL rL 

      IF (iIPMIX .GT. kProfLayer) THEN
        write(kStdErr,*) 'AccumulateForBloat assumes 1..kProfLayer mix paths'
        CALL DoStop
        END IF 

c find out which of the 100 layers is associated with this mixed path 
      iL=MP2Lay(iIpmix) 
 
c find the weight 
      rL=raaMix(iIpmix,iGas) 
 
c add on contribution of the iGas th gas to the iIpmix th row of raaSum 
      DO iFreq=1,kMaxPts 
         raaSum(iFreq,iIpmix)=raaSum(iFreq,iIpmix)+rL*raaGas(iFreq,iL) 
         END DO 
 
      RETURN
      END

c************************************************************************
c this adds together stuff so as to output Mixed paths at high resolution
      SUBROUTINE SumBloatMP(
     $         daFreqBloat,raFreq,raaCo2_LTE,raaRestOfLTEGases,iNLTEStart,
     $                   daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      REAL raaRestOfLTEGases(kMaxPts,kProfLayer)
      REAL raaCo2_LTE(kMaxPts,kProfLayer)
      REAL raFreq(kMaxPts)
      INTEGER iNLTEStart

c input/output params
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

c local variables
      DOUBLE PRECISION daFreqT(kMaxPts),daTemp(kMaxPts)
      DOUBLE PRECISION dY,daY(kBloatPts)
      INTEGER iFr,iL

      DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
        END DO      

c in SetPlanckCoeffBloat we have already done, from iNLTEStart to kProflayer
c       daaNLTEGasAbCoeffBloat(iFr,iL)+blah_from_raaRestOfLTEGases(iFr,iL)   
c so no need to redo it here
c however we ned to do it from 1 t- NLTEStart-1

      DO iL = 1,iNLTEStart-1
        DO iFr = 1,kMaxPts 
          !! we have already stored daaNLTEGasAbCoeffBloat= NLTE_CO2 + LTE_CO2
          !! so we need to add on LTE_others to get total mixed path
          daTemp(iFr)  = raaRestOfLTEGases(iFr,iL)*1.0d0
          END DO
        CALL dlinear_smart(daFreqT,daTemp, kMaxPts,daFreqBloat,daY, kBloatPts) 
        DO iFr = 1,kBloatPts 
          daaSumNLTEGasAbCoeffBloat(iFr,iL) = 
     $                  max(daaNLTEGasAbCoeffBloat(iFr,iL)+daY(iFr),0.0d0)
          END DO 
        END DO 

      DO iL = iNLTEStart,kProfLayer 
        DO iFr = 1,kBloatPts 
          daaSumNLTEGasAbCoeffBloat(iFr,iL) = 
     $                  max(daaNLTEGasAbCoeffBloat(iFr,iL),0.0d0)
          END DO 
        END DO 

      RETURN
      END

c************************************************************************
c this subroutine changes double --> real for the KCARTA atmosphere
c this computes the cumulative multiplication modification to Planck
c this is for the usual 100 kCARTA database layers
c however, this is at the bloated HIGH resolution!!!!!!!
      SUBROUTINE SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer,
     $              raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $              raaSumAbCoeff,raaPlanckCoeff,
     $              raaRestOfLTEGases,raaCO2_LTE,
     $              daFreqBloat,daaSumNLTEGasAbCoeffBloat,
     $              daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      REAL raFreq(kMaxPts)
      DOUBLE PRECISION daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) 
                           !!cumulative nonlte gas optical depths
      DOUBLE PRECISION daaPlanckCoeff(kMaxPts,kProfLayer) !!planck mod, so far
      REAL raaSumAbCoeff(kMaxPts,kMixFilRows)             !!mixed path abscoeff
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raaRestOfLTEGases(kMaxPts,kProfLayer)          !!rest of LTE gases
      REAL raaCO2_LTE(kMaxPts,kProfLayer)                 !!LTE part of CO2
      INTEGER iAtm                                        !!which atmosphere
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)             !!set of mixed paths
      INTEGER iNLTEStart                                !!where nonLTE starts
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      !!for NLTE gas, this contains cumulative optical depths = LTE + NLTE
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) 
c input/output vars
      !!at input, only contains NLTE opt depths for CO2
      !!at output, contains (LTE+NLTE) opt depths for CO2, plus other LTE gases
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)

c local variables
      DOUBLE PRECISION daFreqT(kMaxPts)
      DOUBLE PRECISION dY, daY(kBloatPts), daTemp(kMaxPts)
      DOUBLE PRECISION dY1,daY1(kBloatPts),daTemp1(kMaxPts)
      INTEGER iFr,iL,iA,iAtmTemp,iM,iStart

      DO iFr = 1,kMaxPts
        daFreqT(iFr) = raFreq(iFr)*1.0d0
        END DO

      iAtmTemp = 1
      IF (iAtmTemp .NE. iatm) THEN
        write(kStdErr,*) 'oh oh iAtmTemp .NE. iatm ',iAtmTemp,iatm
        write(kStdErr,*) 'kCARTA assumes the NLTE is for and only for atm #1'
        CALL DoStop
        END IF
      iA = 1    !!!assume current atmospher uses mixed path layers 1-100 
      iM = iaaRadLayer(iAtmTemp,1) 
       
      !find which set of Mixed Paths current atm uses eg 1-100, 101-200 etc 
      DO iL = 1,kProfLayer 
        IF (iM .LE. kProfLayer*iL) THEN 
          iA = iL       
          GOTO 10 
          END IF 
        END DO 
 10   CONTINUE 
       
c initialize to 1.0 below where NLTE starts!!!!
      DO iL = 1,iNLTEStart-1
        DO iFr = 1,kBloatPts
          daaPlanckCoeffBloat(iFr,iL) = 1.0
          END DO
        END DO

      DO iL = iNLTEStart,kProfLayer 
        iM = iL + (iA-1)*kProfLayer 
        DO iFr = 1,kMaxPts 
          !!!the LTE part is easily gotten from raaRestOfLTEGases + raaCo2_LTE
          daTemp(iFr)  = raaRestOfLTEGases(iFr,iM)*1.0d0 + 
     $                   raaCo2_LTE(iFr,iM)*1.0d0
          daTemp1(iFr) = raaRestOfLTEGases(iFr,iM)*1.0d0
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
          dY = max(daaNLTEGasAbCoeffBloat(iFr,iL)+daY1(iFr),0.0d0)
          daaNLTEGasAbCoeffBloat(iFr,iL) = dY + dDeltaNLTE
          !!for Planck, add on LTE component 
          dY = max(daY(iFr)  + daaPlanckCoeffBloat(iFr,iL),0.0d0) 
          dY = max(daY1(iFr) + daaPlanckCoeffBloat(iFr,iL),0.0d0) 
          dY = dY + dDeltaNLTE
          daaPlanckCoeffBloat(iFr,iL) = dY/daaNLTEGasAbCoeffBloat(iFr,iL)
          END DO 
        END DO 

 1010 FORMAT(I6,' ',3(F20.7,' '))
 1111 FORMAT(I6,' ',3(F20.7,' '))

      RETURN
      END

c************************************************************************
c this subroutine simply bloats stuff, on a 1-1 correspondance using i1,i2
      SUBROUTINE AccumulateBloat(daHighX,daHighY,iHigh,daLowX,daLowY,iLow,
     $                           iWideMeshLoop,i1,i2,iType,daBloat)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
c iType = +1 ==> use daHighY at kMaxPtsBox resolution
c                  use the iHigh points from this vector
c iType = -1 ==> use daLowY  at kMaXpts    resolution
c                  use the iCount points from this vector
      DOUBLE PRECISION daHighY(kMaxPtsBox),daHighX(kMaxPtsBox)
      DOUBLE PRECISION daLowY(kMaxPtsBox),daLowX(kMaxPtsBox)
      INTEGER i1,i2,iType,iHigh,iLow,iWideMeshLoop
c output params
      DOUBLE PRECISION daBloat(kBloatPts)

c local vars
      INTEGER iI,iJ,iK,iS,iE,iOffset
      DOUBLE PRECISION daTemp(kBloatPts)

      IF (iType .EQ. 1) THEN
        iOffSet = (iWideMeshLoop-1)*2000
        DO iJ = 1,iHigh-1
          daBloat(iJ+iOffSet) = daBloat(iJ+iOffSet) + daHighY(iJ)
          END DO
      ELSEIF (iType .EQ. -1) THEN
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
      END

c************************************************************************
c this is the driver 
c to bloat up the weak LTE lines (and the stronger LTE lines low down
c in the atm) to higher resolution

c this subroutine bloats up the optical depths, at double precision
c can bloat up all layers, or just select layers
c usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
c at the higher altitudes
      SUBROUTINE BloatCoeffsDriver(iTag,iActualTag,rFileStartFr,raFreq,daaGasAbCoeff,
     $    raaRestOfLTEGases,raaCO2_LTE,
     $    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $    daFreqBloat,iNLTEStart,rNLTEStrength,
     $    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $    iGas,iGasID,iL,iU,iUseWeakBackGnd,
     $    raRAmt,raRTemp,raRPress,raRPartPress,
     $    pProf,iProfileLayers,
     $    raPAmt,raPTemp,raPPress,raPPartPress,iSplineType)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      REAL             raFreq(kMaxPts)                        !!kCARTA res
      DOUBLE PRECISION daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
      DOUBLE PRECISION daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
      DOUBLE PRECISION daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
                                                         !!at kCARTA res
      DOUBLE PRECISION daFreqBloat(kBloatPts)            !!high res pts
      REAL rNLTEStrength                                 !!mixed path strength
      INTEGER iNLTEStart                                 !!where NLTE starts
      INTEGER iTag,iActualTag                            !!tells info about chunk
      REAL rFileStartFr                                  !!start point of chunk

c pProf       = actual layers (from kLAYERS) avg pressure, in iProfileLayers
c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (= kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c iaCont    = whether or not to do continuum calculation .. iaCont(iCount)
c caXecF    = file name of cross section data
c daaTemp   = matrix containing the uncompressed k-spectra
c raVtemp   = vertical temperature profile for the Mixed paths
c iVTSet    = has the vertical temp been set, to check current temp profile
c rFileStartFr   = which k-comp file chunk to uncompress
c iTag      = which k-comp file chunk to uncompress
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ     = analystic Jacobian wrt water amount
c daaDT     = analystic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
      REAL raPPress(kProfLayer),raPPartPress(kProfLayer)
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      INTEGER iProfileLayers
      REAL raVTemp(kProfLayer)
      INTEGER iCount,iGas,iGasID,iL,iU,iErr,iRefLayer,iSplineType
      INTEGER iUseWeakBackGnd
      
c output params
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      REAL raaRestOfLTEGases(kMaxPts,kProfLayer)
      REAL  raaCO2_LTE(kMaxPts,kProfLayer)

c local variables
      INTEGER iLL

      write(kStdWarn,*) 'Using std 0.0025 -> 0.0005 code for "highres" bloat'
      CALL BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff,
     $    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd,
     $    raPAmt,raPTemp,raPPress,raPPartPress,
     $    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $    daFreqBloat,iNLTEStart,rNLTEStrength,
     $    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

      RETURN
      END

c************************************************************************
c this just bloats up the weak LTE lines (and the stronger LTE lines low down
c in the atm) to higher resolution
c this subroutine bloats up the optical depths, at double precision
c can bloat up all layers, or just select layers
c usually used for the LTE stuff (lower altitudes), plus the weak LTE lines
c at the higher altitudes

      SUBROUTINE BloatCoeffsStandard(iTag,rFileStartFr,raFreq,daaGasAbCoeff,
     $    raaRestOfLTEGases,raaCO2_LTE,iUseWeakBackGnd,
     $    raPAmt,raPTemp,raPPress,raPPartPress,
     $    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $    daFreqBloat,iNLTEStart,rNLTEStrength,
     $    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      REAL             raFreq(kMaxPts)                        !!kCARTA res
      DOUBLE PRECISION daaNLTEGasAbCoeff(kMaxPts,kProfLayer)    !!kCARTA res
      DOUBLE PRECISION daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer) !!kCARTA res
      DOUBLE PRECISION daaPlanckCoeff(kMaxPts,kProfLayer)       !!kCARTA res
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer) !!uncompressed LTE,
                                                         !!at kCARTA res
      DOUBLE PRECISION daFreqBloat(kBloatPts)            !!high res pts
      REAL rNLTEStrength                                 !!mixed path strength
      INTEGER iNLTEStart                               !!where NLTE starts
      INTEGER iTag                                     !!tells info about chunk
      REAL    rFileStartFr                             !!start point of chunk
      INTEGER iUseWeakBackGnd                          !!go thru with this?
      REAL raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
      REAL raPPress(kProfLayer),raPPartPress(kProfLayer)

c output params
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer) !!high res
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      REAL raaRestOfLTEGases(kMaxPts,kProfLayer)
      REAL  raaCO2_LTE(kMaxPts,kProfLayer)

c local variables
      DOUBLE PRECISION daTemp(kMaxPts),raaTemp(kMaxPts,kProfLayer)
      DOUBLE PRECISION daFreqT(kMaxPts)
      DOUBLE PRECISION daY(kBloatPts)
      DOUBLE PRECISION df,dfFine,dK
      INTEGER iFr,iL,iM,iA,iIpMix,iS,iE
      REAL raFreqBloat(kBloatPts)
      INTEGER iDo,iLinearOrSpline,iDefault,iStart

      iDefault = -1
      iLinearOrSpline = +1    !!do linear .. not good at all
      iLinearOrSpline = -1    !!do spline .. much better, but overall still sux

      IF (iDefault .NE. iLinearOrSpline) THEN
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

c <---!!!! --------------------------->
      !!! update stuff with the weak LTE lines that have been computed in 
      !!! by calling LTE_spectra. This is for upper atmosphere
      !!set coeffs of upper layers
      iS = iNLTEStart
      iE = kProfLayer
      write (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iS,' up'
      DO iL = iS,iE
        DO iFr = 1,kMaxPts
          daTemp(iFr)  = daaNLTEGasAbCoeff(iFr,iL)
          END DO
        IF (iLinearOrSpline .EQ. 1) THEN
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
        IF (iLinearOrSpline .EQ. 1) THEN
          CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts) 
        ELSE
          CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts) 
          END IF
        DO iFr = 1,kMaxPtsBox
          daaPlanckCoeffBloat(iFr,iL) = max(daY(iFr),0.0d0)
          END DO
        END DO

c <---!!!! --------------------------->
      !!! update stuff with all LTE lines that have been computed from
      !!! uncompressing the kCompressed Database. This is for lower  atmosphere
      !!set coeffs of lower layers
      iS = 1
      iE = iNLTEStart-1
      write (kStdWarn,*) 'Bloat LTE optD, planck mod from layer ',iE,' down'
      DO iL = iS,iE
        DO iFr = 1,kMaxPts
          daTemp(iFr)  = daaGasAbCoeff(iFr,iL)
          END DO
        IF (iLinearOrSpline .EQ. 1) THEN
          CALL dlinear(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts) 
        ELSE
          CALL dspl(daFreqT,daTemp,kMaxPts,daFreqBloat,daY,kBloatPts) 
          END IF
        DO iFr = 1,kMaxPtsBox
          dK = max(daY(iFr),0.0d0)
          daaNLTEGasAbCoeffBloat(iFr,iL)    = dK
          daaSumNLTEGasAbCoeffBloat(iFr,iL) = dK
c          daaPlanckCoeffBloat(iFr,iL)       = 1.0d0
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this subroutine bloats up the LTE coeffs from the upper atmosphere database
      SUBROUTINE bloatUAstuff(
     $           daaUpperNLTEGasAbCoeff,daaUpperPlanckCoeff,iUpper,
     $           raFreq,daFreqBloat,
     $           daaUpperPlanckCoeffBloat,daaUpperNLTEGasAbCoeffBloat,
     $           daaUpperSumNLTEGasAbCoeffBloat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c input params
      INTEGER iUpper
      REAL raFreq(kMaxPts)
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaUpperPlanckCoeff(kMaxPts,kProfLayer)
c output params
      DOUBLE PRECISION daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

c local vars
      INTEGER iFr,iL,iS,iE,iDefault,iLinearOrSpline
      DOUBLE PRECISION dK,daTemp(kMaxPts),daFreqT(kMaxPts),daY(kBloatPts)

      iDefault = +1
      iLinearOrSpline = -1    !!do spline .. worse!
      iLinearOrSpline = +1    !!do linear .. better here!

      IF (iDefault .NE. iLinearOrSpline) THEN
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
        IF (iLinearOrSpline .EQ. 1) THEN
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
        IF (iLinearOrSpline .EQ. 1) THEN
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
      END

c************************************************************************
c this subroutine reads in the precomputed LBL weakbackgnd data 
c as you can see this is hardcoded for 
c   /home/sergio/SPECTRA/IPFILES/layop_md_va0_usethis_co2'
      SUBROUTINE ReadFixedDataWeakBack_Bloat(iGasID,iStart,iTag,rFileStartFr,
     $           raPAmt,raPTemp,raPPress,raPPartPress,
     $           rNLTEStrength,daaWeakOptDepthBloat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

      REAL raPAmt(kProfLayer),raPTemp(kProfLayer),pProf(kProfLayer)
      REAL raPPress(kProfLayer),raPPartPress(kProfLayer)
      INTEGER iStart,iTag,iGasID
      REAL rFileSTartFr
      DOUBLE PRECISION daaWeakOptDepthBloat(kBloatPts,kProfLayer) 
      REAL rNLTEStrength

c local variables
      INTEGER iIoun,I,J,iErr, iL, iFr
      INTEGER IDGAS, NPTS, NLAY
      DOUBLE PRECISION SFREQ, FSTEP
      CHARACTER*120 FNAM
      CHARACTER*4 ca4
      DOUBLEPRECISION df,sf,dK
      INTEGER iDefault
      REAL raFreqBloat(kBloatPts)

      df = kaFineFrStep(iTag)
      sf = 1.0d0*nint(rFileStartFr)
      DO iFr = 1,kBloatPts
        raFreqBloat(iFr) = sf + (iFr-1)*df
        END DO

      FNAM = '/carrot/s1/sergio/AIRSCO2/BACKGND_vt_mdH96/hehe_highres_weaklte_'
      I = 80
 10   CONTINUE
      IF (FNAM(I:I) .EQ. ' ') THEN
        I = I - 1
        GOTO 10 
        END IF
      IF ((nint(rFileStartFr) .LT. 1000) .OR. 
     $    (nint(rFileStartFr) .GT. 9999)) THEN
        write(kStdErr,*) 'Ooooooerr cannot accept freqs outside 1000 ... 9999'
        CALL DoStop
        END IF
 
      write(ca4,40) nint(rFileStartFr)
 40   FORMAT(I4)
      FNAM(I+1:I+4) = ca4(1:4)
      FNAM(I+5:I+8) = '.dat'

      write(kStdWarn,*) 'In ReadFixedDataWeakBack_Bloat, open file for GasID'
      write(kStdWarn,*) iGasID
      write(kStdWarn,*) FNAM

      iIOUN = kCompUnit
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine ReadFixedDataWeakBack_Bloat'
          WRITE(kStdErr,1010) IERR, FNAM
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
          CALL DoSTOP
          ENDIF
      kCompUnitOpen=1

C     Read in the header
      READ(iIOUN) IDGAS, NPTS, NLAY
      READ(iIOUN) SFREQ, FSTEP

      IF (IDGAS .NE. iGasID) THEN
        write(kStdErr,*) 'Wanted gas ',iGasID,'but new data is for gas ',IDGAS
        CALL DoStop
        END IF

      IF (NLAY .NE. kProfLayer) THEN
        write(kStdErr,*) ' new data has ',NLAY,'layers instead of kProfLayer'
        CALL DoStop
        END IF

      IF (NPTS .NE. kBloatPts) THEN
        write(kStdErr,*) ' new data has ',NPTS,' freq pts instead of kBloatPts'
        CALL DoStop
        END IF

      IF (abs(df-FSTEP) .gt. 1e-4) THEN
        write(kStdErr,*) ' new data has ',FSTEP,'freq spacing instead of ',df
        CALL DoStop
        END IF

      IF (abs(sf-SFREQ) .gt. 1e-4) THEN
        write(kStdErr,*) ' new data starts at ',SFREQ,' instead of ',sf
        CALL DoStop
        END IF
      
c passed all tests, so read data
C     Read in the optical depths
      DO I=1,kProfLayer
        READ(iIOUN) (daaWeakOptDepthBloat(J,I),J=1,kBloatPts)
        ENDDO

      CLOSE(iIOUN)
      kCompUnitOpen=-1

      write (kStdWarn,*) 'Read in HIGH RES WEAK BKGND spectra for gasID ',IDGAS

      RETURN
      END

c************************************************************************
c this opens the bloated file four output
      SUBROUTINE OpenOutputBloatFile(
     $                         iType,iNumLayers,caBloatFile,rFrLow,rFrHigh,
     $                         iFileIDLo,iFileIDHi,iTag,iTotalStuff,
     $                         iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

      CHARACTER*120 caBloatFile
      REAL rFrLow,rFrHigh
      INTEGER iFileIDLo,iFileIDHi,iTag,iTotalStuff

      INTEGER iNumLayers !!number of layers in atm #1
      INTEGER iType      !!+1 if regular output; -1 if planck output
      INTEGER iaNLTEChunks(kGasStore)
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

c local vars
      INTEGER iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor
      REAL    r1,r2,r3

      iNLTEChunks = -1
      i1 = +123456
      i2 = -123456
      DO iI = 1,iNumNLTEGases
        IF (iaNLTEChunks(iI) .GE. iNLTEChunks) THEN
          iNLTEChunks = iaNLTEChunks(iI)
          END IF
        DO iJ = 1,iaNLTEChunks(iI)
          IF (iaaNLTEChunks(iI,iJ) .LE. i1) THEN
            i1 = iaaNLTEChunks(iI,iJ)
            END IF
          IF (iaaNLTEChunks(iI,iJ) .GE. i2) THEN
            i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
            END IF
          END DO
        END DO
      r1 = i1*1.0
      r2 = i2*1.0 
      iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
 
      !!!now check to see we actually start from r1, while running kCARTA!
      IF (r1 .LT. rFrLow) THEN
        r1 = rFrLow
        i1 = iFloor(r1)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
        END IF

      !!!now check to see we actually do go to r2, while running kCARTA!
c      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
c        r2 = rFrHigh-kaBlSize(iTag)
      IF (r2 .GT. rFrHigh) THEN
        r2 = rFrHigh
        i2 = iFloor(r2)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
        END IF

      IF (iType .GT. 0) THEN
        iIoun1 = kBloatNLTEOut
      ELSEIF (iType .LT. 0) THEN
        iIoun1 = kBloatNLTEPlanck
      ELSE
        write(kStdErr,*) 'Unknown print option for bloated files'
        CALL DoStop
        END IF

c open unformatted file as a fresh file to be written to
      OPEN(UNIT=iIoun1,FILE=caBloatFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
      IF (iFileErr .NE. 0) THEN
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

 304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file : 
     $      ',/,A120)

      IF (iType .GT. 0) THEN
        kBloatOutOpen = 1
        write(kStdWarn,*) 'Opened following file for bloated general output :'
        write(kStdWarn,*) caBloatFile
      ELSEIF (iType .LT. 0) THEN
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
                                            !!and how many to expect
      WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks 
      WRITE(iIOUN1) kaBlSize(iTag)/kBoxCarUse
                                            !!fine 10000 point freq block size
      WRITE(iIOUN1) iTotalStuff*kBoxCarUse  !!fine number of outputs
      
      RETURN
      END

c************************************************************************
c this subroutine writes out the data header at each standard 10000 pt chunk
      SUBROUTINE HeaderBloatFile(caOutFile,rFrLow,rFrHigh,daFreqBloat,
     $                           iTag,iType)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c iWHich         = +1 for optical depths, -1 for planck modifiers
c raaGasAbs  = single gas abs coeffs
c iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
c iaPath     = list of the paths corresponding to the current gas 
c iNp        = total number of paths to be output
c iaOp       = list of the paths to be output
      REAL             rFrLow,rFrHigh
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      CHARACTER*120 caOutFile
      INTEGER iTag,iType

c local vars
      INTEGER iIOUN1,iI,iJump,ifloor,i10000
      CHARACTER*120 caOut
      DOUBLE PRECISION daStart(kBoxCarUse),daEnd(kBoxCarUse)
      DOUBLE PRECISION daStartBox(kBoxCarUse),daEndBox(kBoxCarUse)

      iJump = ifloor(kBoxCarUse*1.0/2.0)

      caOut = caOutFile

      IF (iType .EQ. +2) THEN
        iIoun1 = kBloatNLTEOutUA
      ELSEIF (iType .EQ. +1) THEN
        iIoun1 = kBloatNLTEOut
      ELSEIF (iType .EQ. -1) THEN
        iIoun1 = kBloatNLTEPlanck
c      ELSEIF (iType .EQ. -2) THEN
c        iIoun1 = kBloatNLTEUAPlanck
      ELSE
        write(kStdErr,*) 'Unknown print option for bloated files'
        CALL DoStop
        END IF

      write(kStdWarn,*) 'These are the chunks for high resolution file ',iIOUN1
      DO iI = 1,kBoxCarUse
        i10000 = kMaxPts*(iI-1)
        !!!these are the high resolution start/stop points for the chunk
        daStart(iI) = daFreqBloat(i10000 + 1)
        daEnd(iI)   = daFreqBloat(i10000 + kMaxPts)
        !!!after boxcar avg, these are the kCARTA resolution start/stop points
        daStartBox(iI) = daFreqBloat(i10000 + 1 + iJump)
        daEndBox(iI)   = daFreqBloat(i10000 + kMaxPts - iJump)
        write(kStdWarn,*) iI,daStart(iI),daEnd(iI),daStartBox(iI),daEndBox(iI)
        END DO

      WRITE(iIOUN1) kMaxPts,rFrLow,rFrHigh,kaFrStep(iTag)   !usual kCARTA stuff
      WRITE(iIOUN1) kBloatPts,kaFineFrStep(iTag),kBoxCarUse !high res stuff
      WRITE(iIOUN1) (daStart(iI),iI=1,kBoxCarUse)
      WRITE(iIOUN1) (daEnd(iI),iI=1,kBoxCarUse)
      WRITE(iIOUN1) (daStartBox(iI),iI=1,kBoxCarUse)
      WRITE(iIOUN1) (daEndBox(iI),iI=1,kBoxCarUse)

      RETURN
      END

c************************************************************************
c this subroutine just dumps out kBoxCarUse worth of chunks of data
      SUBROUTINE wrtout_bloated_rad(iIOUN,caOutBloatFile,
     $                              raFreqBloat,raIntenBloat)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input vars
      REAL         raFreqBloat(kBloatPts),raIntenBloat(kBloatPts)
      INTEGER      iIOUN
      CHARACTER*120 caOutBloatFile

c local vars
      INTEGER iLoop,iInt,iIntOffSet
      REAL raF(kMaxPts),raInten(kMaxPts)

      DO iLoop = 1,kBoxCarUse
        DO iInt=1,kMaxPts
          iIntOffSet    = iInt + (iLoop-1)*kMaxPts
          raF(iInt)     = raFreqBloat(iIntOffSet)
          raInten(iInt) = raIntenBloat(iIntOffset)
          END DO
        CALL wrtout(iIOUN,caOutBloatFile,raF,raInten)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine dumps out gas optical depths at high resolution
      SUBROUTINE out_bloat(raFreq,rFreqStart,rFreqEnd,
     $                  iType,daFreqBloat,
     $                  daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,
     $                  caPlanckBloatFile,caOutBloatFile,
     $                  iFileID,
     $                  iaPath,iNp,iaOp)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c iWHich         = +1 for optical depths, -1 for planck modifiers
c raaGasAbs  = single gas abs coeffs
c iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
c iaPath     = list of the paths corresponding to the current gas 
c iNp        = total number of paths to be output
c iaOp       = list of the paths to be output
      REAL raFreq(kMaxPts),rFreqStart,rFreqEnd
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      INTEGER iPrinter,iFileID,iType
      INTEGER iNp,iaOp(kPathsOut),iaPath(kProfLayer)
      CHARACTER*120 caPlanckBloatFile,caOutBloatFile

c local variables
      INTEGER iIntOffset,iStart
      INTEGER iInt,iDiv,iDp,iPath,iLay,DoOutputLayer,iIOUN,iLoop,iL
      REAL raL2S(kMaxPts),raF(kMaxPts)
      CHARACTER*120 caOut

      IF (iType .GT. 0) THEN
        iIOUN = kBloatNLTEOut
        caOut = caOutBloatFile
      ELSE
        write(kStdErr,*) 'Should be calling out_bloat_planck'
        CALL DoStop
        END IF

c thius certainly works for gas paths; might be wierd for mixed paths
      iStart=iDiv(iaPath(1),kProfLayer)

c write spectra to unformatted file
c if iPrinter=1 then have to check for valid paths
      DO iLay = 1,kProfLayer
c check to see if this path should be output
        iPath = iStart*kProfLayer + iLay
        iDp = DoOutputLayer(iPath,iNp,iaOp)

        IF ((iDp .GT. 0).AND.(iType .GT. 0) .AND. (kBloatOutOpen .GT. 0)) THEN
          IF (iPrinter .EQ. 1) THEN
            write(kStdWarn,*)'output high res GAS opt depths at layer = ',iPath
          ELSEIF (iPrinter .EQ. 2) THEN
            write(kStdWarn,*)'output high res MIXPATH depths at layer = ',iPath
          ELSEIF (iPrinter .EQ. 3) THEN
            write(kStdWarn,*)'output high res RADIANCES at layer = ',iPath
            END IF

          DO iLoop = 1,kBoxCarUse

            IF (kLayer2Sp .EQ. -1) THEN
              DO iInt=1,kMaxPts
                iIntOffSet = iInt + (iLoop-1)*kMaxPts
                raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                raL2S(iInt) = sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
                END DO
              CALL  wrtout(iIOUN,caOut,raF,raL2S)
            ELSE IF (kLayer2Sp .EQ. -2) THEN
              DO iInt=1,kMaxPts
                iIntOffSet = iInt + (iLoop-1)*kMaxPts
                raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                raL2S(iInt)=sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iLay))
                raL2S(iInt) = exp(-raL2S(iInt))
                END DO
              CALL wrtout(iIOUN,caOut,raF,raL2S)
            ELSE IF (kLayer2Sp .EQ. 1) THEN
              DO iInt=1,kMaxPts
                raL2S(iInt) = 0.0
                END DO
              DO iL = iLay,kProfLayer
                DO iInt=1,kMaxPts
                  iIntOffSet = iInt + (iLoop-1)*kMaxPts
                  raL2S(iInt) = raL2s(iInt) + 
     $                          sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
                  END DO
                END DO
              DO iInt=1,kMaxPts
                iIntOffSet = iInt + (iLoop-1)*kMaxPts
                raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
                END DO
              CALL wrtout(iIOUN,caOut,raF,raL2S)
            ELSE IF (kLayer2Sp .EQ. 2) THEN
              DO iInt=1,kMaxPts
                raL2S(iInt) = 0.0
                END DO
              DO iL = iLay,kProfLayer
                DO iInt=1,kMaxPts
                  iIntOffSet = iInt + (iLoop-1)*kMaxPts
                  raL2S(iInt) = raL2s(iInt) + 
     $                          sngl(daaNLTEGasAbCoeffBloat(iIntOffset,iL))
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
          END IF           !IF ((iDp .GT. 0) .AND. (iType .EQ -1,1)) THEN
        END DO             !DO iLay = 1,kProfLayer

      RETURN
      END

c************************************************************************
c this subroutine dumps out gas optical depths at high resolution
      SUBROUTINE out_bloat_planck(raFreq,rFreqStart,rFreqEnd,
     $                  iType,daFreqBloat,
     $                  daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,
     $                  caPlanckBloatFile,caOutBloatFile,
     $                  iFileID,
     $                  iAtm,iNumLayers,iaaRadLayer)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c iWHich         = +1 for optical depths, -1 for planck modifiers
c raaGasAbs  = single gas abs coeffs
c iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
c iaPath     = list of the paths corresponding to the current gas 
c iNp        = total number of paths to be output
c iaOp       = list of the paths to be output
      REAL raFreq(kMaxPts),rFreqStart,rFreqEnd
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      INTEGER iPrinter,iFileID,iType
      INTEGER iNp,iaOp(kPathsOut),iaPath(kProfLayer)
      CHARACTER*120 caPlanckBloatFile,caOutBloatFile
      INTEGER iAtm,iNumLayers,iaaRadLayer(kMaxAtm,kProfLayer)

c local variables
      INTEGER iIntOffset,iaRadLayer(kProfLayer)
      INTEGER iInt,iDiv,iDp,iPath,iLay,DoOutputLayer,iIOUN,iLoop,iL
      REAL raL2S(kMaxPts),raF(kMaxPts)
      CHARACTER*120 caOut

      IF (iType .GT. 0) THEN
        write(kStdErr,*) 'Should be calling out_bloat'
        CALL DoStop
      ELSE
        iIOUN = kBloatNLTEPlanck
        caOut = caPlanckBloatFile
        END IF
 
      DO iInt = 1,iNumLayers
        iaRadLayer(iInt) = iaaRadLayer(iAtm,iInt)
        END DO

c write spectra to unformatted file
c if iPrinter=1 then have to check for valid paths
      DO iLay = 1,iNumLayers
c check to see if this path should be output
        iPath = iaRadLayer(iLay)
        IF ((iType .LT. 0) .AND. (kBloatPlanckOpen .GT. 0)) THEN
          write(kStdWarn,*)'output high res planck coeffs at layer = ',iPath
          DO iLoop = 1,kBoxCarUse
            DO iInt=1,kMaxPts
              iIntOffSet = iInt + (iLoop-1)*kMaxPts
              raF(iInt)   = sngl(daFreqBloat(iIntOffSet))
              raL2S(iInt) = sngl(daaPlanckCoeffBloat(iIntOffset,iPath))
              END DO
            CALL wrtout(iIOUN,caOut,raF,raL2S)
            END DO         !DO iLoop = 1,kBoxCarUse

          END IF           !IF ((iDp .GT. 0) .AND. (iType .EQ -1,1)) THEN
        END DO             !DO iLay = 1,kProfLayer

      RETURN
      END

c************************************************************************
c the main dog! radiative transfer! DOUBLE precision
c BIG Assumption, is that we only dump out radiances at TOP of layers
c compared to rad_main, where we could dump radiances in midlayer
      SUBROUTINE radiances_bloat(
     $              raFreq,raVTemp,caOutBloatFile,
     $              iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $              rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $              rSatAngle,rFracTop,rFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              daFreqBloat,daaSumNLTEOptDepthBloat,daaPlanckCoeffBloat,
     $    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $    daaUpperNLTEGasAbCoeffBloat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNLTEStart  = which layer NONLTE calcs start
c raaPlanckCoeff = how to affect the Planck computation
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutBloatFile  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for computing radiances
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag,iDumpAllUARads
      CHARACTER*120 caOutBloatFile
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c these are the bloated lower atmosphere matrices
      DOUBLE PRECISION daaSumNLTEOptDepthBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      DOUBLE PRECISION daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
      DOUBLE PRECISION daaLayTrans(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaEmission(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaLay2Sp(kBloatPts,kProfLayer)
      REAL raSumLayEmission(kBloatPts),raSurfaceEmissionToSpace(kBloatPts)
      DOUBLE PRECISION daSumLayEmission(kBloatPts),
     $                 daSurfaceEmissionToSpace(kBloatPts)
      REAL rDum1,rDum2,ttorad,rOmegaSun,rCO2
c to do the thermal,solar contribution
      REAL             raIntenBloat(kBloatPts)
      DOUBLE PRECISION daIntenBloat(kBloatPts)
      REAL             raSurfaceBloat(kBloatPts),raEmissivityBloat(kBloatPts)
      DOUBLE PRECISION daSurfaceBloat(kBloatPts),daEmissivityBloat(kBloatPts)
      REAL             raSunBloat(kBloatPts)
      DOUBLE PRECISION daSunBloat(kBloatPts)
      REAL             raThermalBloat(kBloatPts),raSunReflBloat(kBloatPts)
      DOUBLE PRECISION daThermalBloat(kBloatPts),daSunReflBloat(kBloatPts)
      REAL             rThermalRefl,rCos,rPlanck,rMPTemp
      DOUBLE PRECISION dThermalRefl,dCos,dPlanck,dMPTemp,dttorad
      REAL             raVT1(kMixFilRows),rXYZ
      DOUBLE PRECISION daVT1(kMixFilRows),dXYZ
      INTEGER iDoThermal,iDoSolar,MP2Lay

      DOUBLE PRECISION dEmission, dTrans               
      REAL raOutFrac(kProfLayer),raFreqBloat(kBloatPts),InterpTemp
      INTEGER iIOUN
      REAL bt2rad,t2s
      INTEGER iFr1

      INTEGER i1,i2,iFloor,iDownWard,iSTopNormalRadTransfer
     
c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2=iaaRadLayer(iAtm,iNumLayer)-1
        i2=iFloor(i2*1.0/kProfLayer)
        IF (rTSpace .GT. 5.0) THEN
          write(kStdErr,*) 'you want satellite to be downward looking'
          write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
          write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
          write(kStdErr,*) 'Please retry'
          CALL DoSTOP
          END IF
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,iNumLayer))THEN
c radiation travelling downwards to instrument ==> sat looking up
c i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = -1
        i1=iaaRadLayer(iAtm,1)-1
        i1=iFloor(i1*1.0/(1.0*kProfLayer))
        i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
        write(kStdErr,*) 'Huh ? NLTE for UPLOOK instrument?????'
        CALL DoStop
        END IF
      write(kStdWarn,*) ' ------------------> <----------------------'
      write(kStdWarn,*) 'Doing BLOATED nlte rad transfer ...'
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUp=90, iLow=1 is acceptable
c eg iUp=140,iLow=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
        END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
        END IF

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

c here we go .....

      iIOUN = kBloatNLTEOut

      rThermalRefl=1.0/kPi
      dThermalRefl=1.0/kPi
      
      DO iFr = 1,kBloatPts
        raFreqBloat(iFr) = sngl(daFreqBloat(iFr))
        END DO

c calculate cos(SatAngle)
      rCos = cos(rSatAngle*kPi/180.0)
      dCos = dble(rCos)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop


c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
        END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
          END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
          END IF
        IF (iaRadLayer(iLay) .GT. kProfLayer) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Err, assume only 1 .. kProfLayer mixedpaths'
          write(kStdErr,*) '  in the high res matrices'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
          END IF
        END DO

c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
        daVT1(iFr) = dble(raVTemp(iFr))
        END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      daVT1(iL)=dble(raVT1(iL))
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      daVT1(iL)=dble(raVT1(iL))
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh=iaOp(iLay)
          END IF
        END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay=1,1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         dCos = DBLE(rCos)
         DO iFr=1,kBloatPts
           dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)*dble(rFracBot)/dCos
           daaLayTrans(iFr,iLay) = dexp(-dXYZ)
           daaEmission(iFr,iLay) = 0.0d0
           END DO
         END DO

       DO iLay=2,iNumLayer-1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         dCos = dble(rCos)
         DO iFr=1,kBloatPts
           dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)/dCos
           daaLayTrans(iFr,iLay) = dexp(-dXYZ)
           daaEmission(iFr,iLay) = 0.0d0
           END DO
         END DO
       DO iLay=iNumLayer,iNumLayer
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         dCos = dble(rCos)
         DO iFr=1,kBloatPts
           dXYZ = daaSumNLTEOptDepthBloat(iFr,iL)*dble(rFracTop)/dCos
           daaLayTrans(iFr,iLay) = dexp(-dXYZ)
           daaEmission(iFr,iLay) = 0.0d0
           END DO
         END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NONLTE computations done!
c so usually only the usual LTE computations are done!!
      IF (iNLTEStart .GT. kProfLayer) THEN
        iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
        write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
      ELSE
        iLay = 1
 987    CONTINUE
        iL=iaRadLayer(iLay)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer .EQ. 0) THEN
          iLModKprofLayer = kProfLayer
          END IF
        IF ((iLModKprofLayer .LT. iNLTEStart).AND.(iLay .LT. iNumLayer)) THEN
          iLay = iLay + 1
          GOTO 987
          END IF
        iSTopNormalRadTransfer = iLay
        write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
        END IF

      DO iLay = 1,iNumLayer
        iL      = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        dMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer .EQ. 0) THEN
          iLModKprofLayer = kProfLayer
          iL = kProfLayer
          END IF
        IF (iLModKprofLayer .LT. iNLTEStart) THEN   
          !normal, no LTE emission stuff
          DO iFr=1,kBloatPts
            dPlanck = dttorad(daFreqBloat(iFr),dMPTemp)
            daaEmission(iFr,iLay) = (1.0-daaLayTrans(iFr,iLay))*dPlanck
            END DO
        ELSEIF (iLModKprofLayer .GE. iNLTEStart) THEN
          !new; LTE emission stuff
          DO iFr=1,kBloatPts
            dPlanck = dttorad(daFreqBloat(iFr),dMPTemp)	  
            dPlanck = dPlanck*daaPlanckCoeffBloat(iFr,iL)
            daaEmission(iFr,iLay) = (1.0-daaLayTrans(iFr,iLay))*dPlanck
            END DO
          END IF
        END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL rspl(raFreq,raThermal,kMaxPts,raFreqBloat,
     $             raThermalBloat,kBloatPts)
      ELSE
        DO iFr = 1,kBloatPts
          raThermalBloat(iFr) = 0.0
          END DO
        write(kStdWarn,*) 'no thermal backgnd to calculate'
        END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL rspl(raFreq,raSun,kMaxPts,raFreqBloat,raSunBloat,kBloatPts)
        CALL rspl(raFreq,raSunRefl,kMaxPts,raFreqBloat,
     $             raSunReflBloat,kBloatPts)
      ELSE
        DO iFr = 1,kBloatPts
          raSunBloat(iFr) = 0.0
          END DO
        write(kStdWarn,*) 'no solar backgnd to calculate'
        END IF

      CALL rlinear(raFreq,raSurface,kMaxPts,
     $             raFreqBloat,raSurfaceBloat,kBloatPts)
      CALL rlinear(raFreq,raUseEmissivity,kMaxPts,
     $             raFreqBloat,raEmissivityBloat,kBloatPts)

      DO iFr=1,kBloatPts
        raIntenBloat(iFr)=raSurfaceBloat(iFr)*raEmissivityBloat(iFr)+
     $          raThermalBloat(iFr)*(1.0-raEmissivityBloat(iFr))*rThermalRefl+
     $          raSunBloat(iFr)*raSunReflBloat(iFr)
        daIntenBloat(iFr) = dble(raIntenBloat(iFr))
        END DO

c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
        iL      = iaRadLayer(iLay)
        rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        dCos    = dble(rCos)
        rMPTemp = raVT1(iL)
        dMPTemp = dble(raVT1(iL))

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
c do the radiative transfer thru this bottom layer
        DO iFr=1,kBloatPts
          daIntenBloat(iFr) = daaEmission(iFr,iLay) +
     $                        daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
          raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
          END DO

        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .EQ. 1) THEN
          write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
          CALL wrtout_bloated_rad(iIOUN,
     $                            caOutBloatFile,raFreqBloat,raIntenBloat)
        ELSEIF (iDp .GE. 2) THEN
          write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out' 
          write(kStdErr,*) 'more than one radiance per layer'
          CALL DoStop
          END IF
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
        iL      = iaRadLayer(iLay)
        rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        dCos    = dble(rCos)
        rMPTemp = raVT1(iL)
        dMPTemp = dble(raVT1(iL))

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
c now do the radiative transfer thru this complete layer
        DO iFr=1,kBloatPts
          daIntenBloat(iFr) = daaEmission(iFr,iLay) +
     $                        daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
          raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
          END DO

        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .EQ. 1) THEN
          write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
          CALL wrtout_bloated_rad(iIOUN,
     $                            caOutBloatFile,raFreqBloat,raIntenBloat)
        ELSEIF (iDp .GE. 2) THEN
          write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out' 
          write(kStdErr,*) 'more than one radiance per layer'
          CALL DoStop
          END IF

        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
 777  CONTINUE
      DO iLay=iHigh,iHigh
        iL      = iaRadLayer(iLay)
        rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        dCos    = dble(rCos)
        rMPTemp = raVT1(iL)
        dMPTemp = dble(raVT1(iL))

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN

            write(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer' 
            write(kStdWarn,*) 'This is the top of the usual AIRS atmosphere' 
            write(kStdWarn,*) '   you have iDoUpperATM > 0' 
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere' 
            write(kStdWarn,*) 'and output stuff into the blah_UA file' 
            write(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of' 
            write(kStdWarn,*) 'stratosphere into both this and the UA file' 
 
            !do radiative transfer thru this layer 
            DO iFr=1,kBloatPts 
              daIntenBloat(iFr) =  
     $          daaEmission(iFr,iLay)+daIntenBloat(iFr)*daaLayTrans(iFr,iLay) 
              raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
              END DO 

            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTransBloatDouble(raFreqBloat,
     $        daIntenBloat,daFreqBloat,raLayAngles(MP2Lay(iL)),
     $        iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of stratosphere (120-200 km)

            write(kStdWarn,*) 'outputting bloated stuff at TOTAL Complete TOA into' 
            write(kStdWarn,*) 'usual BLOATED binary file (iLay = ',iLay,')' 

            DO iFr=1,kBloatPts 
              raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
              END DO
            DO iFr=1,iDp
                CALL wrtout_bloated_rad(iIOUN,
     $                            caOutBloatFile,raFreqBloat,raIntenBloat)
              END DO
            END IF
          END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           DO iFr = 1,kBloatPts
             daIntenBloat(iFr) = daaEmission(iFr,iLay) +
     $                           daIntenBloat(iFr)*daaLayTrans(iFr,iLay)
             raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
             END DO
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .EQ. 1) THEN
            write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' rad layer'
            CALL wrtout_bloated_rad(iIOUN,
     $                              caOutBloatFile,raFreqBloat,raIntenBloat)
          ELSEIF (iDp .GE. 2) THEN
            write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out' 
            write(kStdErr,*) 'more than one radiance per layer'
            CALL DoStop
            END IF
          END IF

        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

 3579 FORMAT(I4,' ',F10.5,' ',5(E12.5,' '))

      RETURN
      END

c************************************************************************
c this subroutine very quickly does the radiative transfer (DOUBLE PREC)
c since the optical depths are soooooooooo small, use double precision
      SUBROUTINE UpperAtmRadTransBloatDouble(
     $        raFreqBloat,
     $        daIntenBloat,daFreqBloat,rSatAngle,
     $        iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters 
c   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c   upper atm abs coeff and planck coeff
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
c   input wavevector and integer stating which layer to stop rad transfer at
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      REAL raFreqBloat(kBloatPts),rSatAngle
      INTEGER iUpper
c do we want to do upper atm NLTE computations?
      INTEGER iDoUpperAtmNLTE
c do we dump all or some rads?
      INTEGER iDumpAllUARads
c input/output pararameter
c   this contains the radiance incident at kCARTA TOA (0.005 mb) 
c   it will finally contain the radiance exiting from TOP of UPPER ATM
      DOUBLE PRECISION daIntenBloat(kBloatPts)

c local variables
      INTEGER iFr,iL,iIOUN
      DOUBLE PRECISION dEmission,dTrans,rMu,daInten0(kBloatPts)
      REAL raIntenBloat(kBloatPts),ttorad
      CHARACTER*120 caOutName

      caOutName = 'DumDum'
      iIOUN = kBloatNLTEOutUA
  
      IF (iDoUpperAtmNLTE .LE. 0) THEN
        write (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
        CALL DoStop
      ELSE
        write(kStdWarn,*) 'Doing UA (NLTE) radtransfer at highres 0.0005 cm-1 '
      END IF

      !!compute radiance intensity thru NEW uppermost layers of atm 
      DO iFr = 1,kBloatPts
        daInten0(iFr) = daIntenBloat(iFr)
        raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
        END DO

      iL = 0
      IF (kNLTEOutUAOpen .GT. 0) THEN
        write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
        !!always dump out the 0.005 mb TOA radiance if the UA file is open
        CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
        END IF

      rMu = cos(rSatAngle*kPi/180.0)

      DO iL = 1,iUpper-1

        DO iFr = 1,kBloatPts
          dTrans = daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)/(rMu*1.0d0)
          dTrans = exp(-dTrans)
          dEmission = daaUpperPlanckCoeffBloat(iFr,iL) *
     $                  dble(ttorad(raFreqBloat(iFr),raUpperTemp(iL))*1.0d0)*
     $                  (1.0d0 - dTrans)
          daIntenBloat(iFr) = dEmission + daIntenBloat(iFr)*dTrans

          raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
          END DO

        IF ((iDumpAllUARads .GT. 0) .AND. (kNLTEOutUAOpen .GT. 0)) THEN
          write(kStdWarn,*) 'dumping out UA highres rads at iL = ',iL
          !!dump out the 0.000025 mb TOA radiance
          CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
          END IF
        END DO

      DO iL = iUpper,iUpper
        DO iFr = 1,kBloatPts
          dTrans = daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)/(rMu*1.0d0)
          dTrans = exp(-dTrans)
          dEmission = daaUpperPlanckCoeffBloat(iFr,iL) *
     $                  dble(ttorad(raFreqBloat(iFr),raUpperTemp(iL))*1.0d0)*
     $                  (1.0d0 - dTrans)
          daIntenBloat(iFr) = dEmission + daIntenBloat(iFr)*dTrans
          raIntenBloat(iFr) = sngl(daIntenBloat(iFr))
          END DO

        IF (kNLTEOutUAOpen .GT. 0) THEN
          !!always dump out 0.000025 mb TOA radiance if the UA file is open
          write(kStdWarn,*) 'dumping out 0.000025 mb highres UA rads iL = ',iL
          !!dump out the 0.000025 mb TOA radiance
          CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
          END IF

        END DO

 11   FORMAT(I3,' ',9(E15.8,' '))
 3579 FORMAT(I4,' ',F10.5,' ',5(E12.5,' '))

      RETURN
      END

c************************************************************************
c************************************************************************
c************************************************************************
