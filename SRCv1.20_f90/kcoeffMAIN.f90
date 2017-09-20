! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeffMAIN

use kcoeffSPL
use kcoeffSPLJAC
use kcoeff_FAST
use kcoeff_FAST_details
use kcoeff_FAST_details2
use kcont_xsec

IMPLICIT NONE

CONTAINS


!************************************************************************
!********* this file has the main k-compressed routines *****************
!** which include reading in the data, doing the spline interpolations **
!********** and finding the absorption matrix for the relevant gas ******
!************************************************************************

!************************************************************************
!*********** NOTE THESE VARIABLES ARE DOUBLE PRECISION ******************
!************************************************************************

!  if kJacobian == 1 then call the spline-type routines
!  if kJacobian == -1 then call OLD routines
! where JAC( indicates routines for calculating jacobians, and NOJAC(
! indicates original routines w/o jacobians

! Also, if iGasID > kMaxDQ then we do not need to calculate d/dq, only d/dT
! contribution. Since kMaxDQ >= 1, water d/dq,d/dT always calculated
!************************************************************************
! this is the MAIN routine
    SUBROUTINE UsualLTEUncompress(iGas,iaGases, &
    raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
    raTAmt,raTTemp,raTPress,raTPartPress,iaCont, &
    pProf,iProfileLayers, &
    raVertTemp,iVertTempSet, &
    rFileStartFr,iTag,iActualTag,raFreq,iError,iDoDQ,iSplineType, &
    iNumNewGases,iaNewGasID,caaaNewChunks,iaNewData,iaaNewChunks, &
    iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr, &
    daaDQ,daaDT,daaGasAbCoeff, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iaNewGasID(kGasStore),iaNewData(kGasStore)
    INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)
    INTEGER :: iaAltComprDirs(kGasStore),iNumAltComprDirs
    CHARACTER(80) :: caaAltComprDirs(kGasStore)
    INTEGER :: iGas,iaGases(kMaxGas)
    INTEGER :: iProfileLayers,iVertTempSet,iError,iDoDQ
    INTEGER :: iSplineType,iL_low,iL_high,iTag,iActualTag
    INTEGER :: iaCONT(kMaxGas)
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
    REAL :: pProf(kProfLayer),raVertTemp(kProfLayer),raFreq(kMaxPts)
    REAL :: rFileStartFr
    REAL :: rAltMinFr,rAltMaxFr
! the Matlab weights
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

! output vars
! daaDT,daaDQ are the d/dq,d/dT matrices
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! daaGasAbCoeff has the uncompressed gas absorption coeff
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)

! local vars
    INTEGER :: iNewIn,OutSideSpectra,NewDataChunk,iWhichChunk
    INTEGER :: ix

    kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath

    iNewIn = -1
    iNewIn = OutsideSpectra(iaGases(iGas),iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs, &
    rFileStartFr,rAltMinFr,rAltMaxFr,iTag)
    IF (iNewIn < 0) THEN
    ! se kCompressed Database w/o worrying
        kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
        CALL GasContribution(iGas,iaGases(iGas),kProfLayer, &
        raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
        pProf,iProfileLayers, &
        raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile, &
        raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,raFreq, &
        iError,iDoDQ, &
        daaDQ,daaDT,daaGasAbCoeff,iSplineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    END IF

    IF ((iNewIn >= 1) .AND. (iNewIn < 1000)) THEN
    !!!! new monochromatic spectra from file
        iWhichChunk=NewDataChunk(iNewIn,iaNewData,iaaNewChunks, &
        rFileStartFr)
        IF (iWhichChunk > 0) THEN
        ! ead in new spectra
            CALL ReadNewData(iGas,iaGases(iGas),kProfLayer, &
            iL_low,iL_high,raTAmt,raTTemp,raTPress,raTPartPress, &
            iaCont,iTag,iActualTag,raFreq,daaDQ,daaDT,iDoDQ, &
            daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks, &
            kaFrStep(iTag),rFileStartFr*1.00000)
        ELSE
        ! se kCompressed Database w/o worrying
            kAltComprDirs = -1  !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
            CALL GasContribution(iGas,iaGases(iGas),kProfLayer, &
            raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
            pProf,iProfileLayers, &
            raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile, &
            raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag, &
            raFreq,iError,iDoDQ, &
            daaDQ,daaDT,daaGasAbCoeff,iSplineType, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        END IF
    END IF

    IF ((iNewIn >= 1001) .AND. (iNewIn < 10000)) THEN
    !!!! use alternate compressed database
        kAltComprDirs = +1 !!overwrite one of kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
        CALL GasContributionAlternateDataBase(iGas,iaGases(iGas),kProfLayer, &
        raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
        pProf,iProfileLayers, &
        raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile, &
        raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag, &
        raFreq,iError,iDoDQ, &
        daaDQ,daaDT,daaGasAbCoeff,iSplineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22, &
        iNewIN-1000,iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr)
    END IF

!      if (iGas .EQ. 1) then
!        print *, '%kcoeffmain : check me out '
!c        do ix = 1,100
!c          print *, 'kcoeffmain : check me out ',ix,daaDT(4602,ix)
!c          end do
!        do ix = 1,10000
!          print *, ix,daaDT(ix,4)
!         end do
!       end if

    RETURN
    end SUBROUTINE UsualLTEUncompress

!************************************************************************
! this subroutine gets the contribution of the i-th gas to the
! absorbtion coefficients. Using the reference profile, it scales the
! contribution accordingly
! also, the gases are weighted by raaMix (from MIXFIL)
! iGasID is the gas ID, while rFileStartFr identifies the wavenumber block
    SUBROUTINE GasContribution(iCount,iGasID,iRefLayer, &
    raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,pProf,iProfileLayers, &
    raAmt,raTemp,raPress,raPartPress,iaCont,caXsecF, &
    raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,raFreq,iErr,iDoDQ, &
    daaDQ,daaDT,daaTemp, &
    iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf     = actual layers (from kLAYERS) avg pressure, in iProfileLayers
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raAmt(kProfLayer),raTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPress(kProfLayer),raPartPress(kProfLayer)
    INTEGER :: iaCont(kMaxGas),iVTSet,iTag,iProfileLayers
    REAL :: raVTemp(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iCount,iGasID,iL,iU,iErr,iRefLayer
    INTEGER :: iActualTag,iDoDQ,iSplineType
    CHARACTER(80) :: caXsecF
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer),kFrStep,rFileStartFr
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the Matlab weights
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

! local variables
    INTEGER :: iFr,iLay

    iErr = -1

    kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath

    IF (((1 <= iGasID) .AND. (iGasID <= kGasComp)) .OR. &
    (iGasID == kNewGasHi+1)) THEN
        kFrStep = kaFrStep(iTag)
        CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress, &
        raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag, &
        raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont, &
        pProf,iProfileLayers,iDoDQ, &
        daaDQ,daaDT,daaTemp,iSPlineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    END IF

    IF ((kGasXsecLo <= iGasID) .AND. (iGasID <= kGasXsecHi)) THEN
        kFrStep = kaFrStep(iTag)
        IF (KXsecFormat > 0) THEN
            write (kStdWarn,*) ' xsec gas : using kCompressed Database format'
            CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress, &
            raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag, &
            raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont, &
            pProf,iProfileLayers,iDoDQ, &
            daaDQ,daaDT,daaTemp,iSPlineType, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
            write (kStdWarn,*) ' xsec gas : using old style binary file format'
            CALL CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep,daaTemp, &
            raVTemp,iVTSet,raFreq,iErr,caXsecF, &
            raAmt,raTemp,raPress,raPartPress,iaCont, &
            daaDQ,daaDT,iDoDQ)
        END IF
    END IF

    IF ((kNewGasLo <= iGasID) .AND. (iGasID <= kNewGasHi)) THEN
        IF (kCKD >= 0) THEN           !add on continuum
            kFrStep = kaFrStep(iTag)
            CALL driver_continuum(iCount,iGasID,iRefLayer,iProfileLayers, &
            raRAmt,raRTemp,raRPress,raRPartPress, &
            iL,iU,daaTemp,raVTemp,iVTSet, &
            rFileStartFr,iTag,iActualTag, &
            raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont, &
            daaDQ,daaDT,iDoDQ,iSplineType,pProf, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
            write(kStdWarn,*)'even if 101,102 are valid gases, '
            write(kStdWarn,*)'not adding on continuum as kCKD < 0'
            DO iLay=1,kProfLayer
                DO iFr=1,kMaxPts
                    daaTemp(iFr,iLay)=0.0d0
                END DO
            END DO
            IF (iDoDQ > -2) THEN
                DO iLay=1,kProfLayer
                    DO iFr=1,kMaxPts
                        daaDQ(iFr,iLay)=0.0d0
                    END DO
                END DO
            END IF
        END IF
    END IF

    RETURN
    end SUBROUTINE GasContribution

!************************************************************************
! this subroutine gets the contribution of the i-th gas to the
! absorbtion coefficients. Using the reference profile, it scales the
! contribution accordingly
! also, the gases are weighted by raaMix (from MIXFIL)
! iGasID is the gas ID, while rFileStartFr identifies the wavenumber block
! same as GasContribution except it substitudes COMPRESSED DATABASE
    SUBROUTINE GasContributionAlternateDatabase(iCount,iGasID,iRefLayer, &
    raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,pProf,iProfileLayers, &
    raAmt,raTemp,raPress,raPartPress,iaCont,caXsecF, &
    raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,raFreq,iErr,iDoDQ, &
    daaDQ,daaDT,daaTemp, &
    iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22, &
    iNewIN,iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf     = actual layers (from kLAYERS) avg pressure, in iProfileLayers
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raAmt(kProfLayer),raTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPress(kProfLayer),raPartPress(kProfLayer)
    INTEGER :: iaCont(kMaxGas),iVTSet,iTag,iProfileLayers
    REAL :: raVTemp(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iCount,iGasID,iL,iU,iErr,iRefLayer
    INTEGER :: iActualTag,iDoDQ,iSplineType
    CHARACTER(80) :: caXsecF
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer),kFrStep,rFileStartFr
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the Matlab weights
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
! the alt database
    INTEGER :: iNewIN,iNumAltComprDirs,iaAltComprDirs(kGasStore)
    CHARACTER(80) :: caaAltComprDirs(kGasStore)
    REAL :: rAltMinFr,rAltMaxFr

! local variables
    INTEGER :: iFr,iLay,strfind

    iErr = -1

    DO iFr = 1,80
        kcaAltComprDirs(iFr:iFr) = ' '
    END DO
    kcaAltComprDirs = caaAltComprDirs(iNewIN)
    write(kStdWarn,*) '>>> substituting caCompressedDataPath for gasID ',iGasID
    write(kStdWarn,80) kcaAltComprDirs

    IF (iGasID == 2) THEN
        IF ((strfind(kcaAltComprDirs,'lblrtm') == 1) .OR. (strfind(kcaAltComprDirs,'LBLRTM') == 1)) THEN
            IF (kCO2_UMBCorHARTMAN > 0) THEN
                write(kStdWarn,*) ' kCO2_UMBCorHARTMAN is +1 for UMBC CO2 linemix',kCO2_UMBCorHARTMAN
                write(kStdWarn,*) ' ignore so NO chi fcns, since you are using LBLRTM CO2 database'
                write(kStdErr,*)  ' ignore kCO2_UMBCorHARTMAN (+1) so NO chi fcns, for LBLRTM CO2 database'
            ! CO2_UMBCorHARTMAN = -1
            END IF
        END IF
    END IF
            
    IF ( ((1 <= iGasID) .AND. (iGasID <= kGasComp)) .OR. (iGasID == kNewGasHi+1)) THEN
        kFrStep = kaFrStep(iTag)
        CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress, &
        raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag, &
        raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont, &
        pProf,iProfileLayers,iDoDQ, &
        daaDQ,daaDT,daaTemp,iSPlineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    END IF

    IF ((kGasXsecLo <= iGasID) .AND. (iGasID <= kGasXsecHi)) THEN
        kFrStep = kaFrStep(iTag)
        IF (KXsecFormat > 0) THEN
            write (kStdWarn,*) ' xsec gas : using kCompressed Database format'
            CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress, &
            raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag, &
            raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont, &
            pProf,iProfileLayers,iDoDQ, &
            daaDQ,daaDT,daaTemp,iSPlineType, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
        ELSE
            write(kStdWarn,*) ' xsec gas : using old style binary file format'
            write(kStdErr,*)  'not supported here in GasContributionAlternateDatabase'
            CALL DoStop
        END IF
    END IF

    IF ((kNewGasLo <= iGasID) .AND. (iGasID <= kNewGasHi)) THEN
        IF (kCKD >= 0) THEN           !add on continuum
            kFrStep = kaFrStep(iTag)
            write(kStdErr,*)  'not supported here in GasContributionAlternateDatabase'
            CALL DoStop
        END IF
    END IF

    80 FORMAT(A80)

    RETURN
    end SUBROUTINE GasContributionAlternateDatabase

!************************************************************************
! this subroutine basically computes water continuum
    SUBROUTINE driver_continuum(iCount,iGasID,iRefLayer,iProfileLayers, &
    raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,daaTemp,raVTemp,iVTSet, &
    rFileStartFr,iTag,iActualTag, &
    raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont, &
    daaDQ,daaDT,iDoDQ,iSplineType,pProf, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers (from kLAYERS) avg pressure, in iProfileLayers
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
! daaDQ      = analytic Jacobian wrt amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raTPart(kProfLayer),raTPress(kProfLayer),rFileStartFr
    INTEGER :: iaCont(kMaxGas),iVTSet,iDoDQ
    INTEGER :: iTag,iSPlineType,iActualTag,iProfileLayers
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    REAL :: raVTemp(kProfLayer),raFreq(KmaxPts)
    INTEGER :: iGasID,iL,iU,iErr,iCount,iRefLayer
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the Matlab weights
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

    INTEGER :: iAns
    REAL :: kFrStep

    iErr=0

    IF (iGasID == kNewGasLo) THEN
    ! check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep = kaFrStep(iTag)
        IF ((iAns > 0) .AND. (kCKD >= 0)) THEN
            CALL AddContinuum(iGasID,iTag,iActualTag,rFileStartFr, &
            iRefLayer,iProfileLayers,raFreq,raTAmt,raTTemp, &
            kFrStep,raTPress,raTPart,iL,iU,daaTemp, &
            daaDQ,daaDT,iDoDQ,iSPlineType, &
            raRAmt,raRTemp,raRPress,raRPartPress,pProf, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
            write(kStdWarn,*) 'added self water continuum'
        END IF
    END IF

    IF (iGasID == kNewGasHi) THEN
    ! check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep = kaFrStep(iTag)
        IF ((iAns > 0) .AND. (kCKD >= 0)) THEN
            CALL AddContinuum(iGasID,iTag,iActualTag,rFileStartFr, &
            iRefLayer,iProfileLayers,raFreq,raTAmt,raTTemp, &
            kFrStep,raTPress,raTPart,iL,iU,daaTemp, &
            daaDQ,daaDT,iDoDQ,iSPlineType, &
            raRAmt,raRTemp,raRPress,raRPartPress,pProf, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)
            write(kStdWarn,*) 'added foreign water continuum'
        END IF
    END IF

    RETURN
    end SUBROUTINE driver_continuum

!************************************************************************
! this subroutine computes the contribution of Gases 1-kGasComp (if present)
    SUBROUTINE compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress, &
    raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag, &
    raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont, &
    pProf,iProfileLayers,iDoDQ, &
    daaDQ,daaDT,daaTemp, &
    iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! pProf       = actual layers from kLAYERS avg pressure
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
! daaDQ      = analytic Jacobian wrt amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raTPart(kProfLayer),raTPress(kProfLayer),rCheckTemp
    INTEGER :: iaCont(kMaxGas),iVTSet,iDoDQ,iTag,iSPlineType,iActualTag
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    REAL :: raVTemp(kProfLayer),raFreq(KmaxPts),rFileStartFr
    INTEGER :: iGasID,iL,iU,iErr,iCount,iRefLayer,iProfileLayers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! the Matlab weights
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

    INTEGER :: iErr1,iCnt,iMatlabORf77,iDefault

    iErr1=0    !this sees if data successfully uncompressed

    IF (iErr1 <= 0) THEN
    ! set the vertical temperature profile if iGasID < 51 (first profile read)
        IF ((iGasID <= kGasComp) .AND. (iVTSet < 0)) THEN
            write(kStdWarn,* )'Setting the vertical tempr profile ...'
            DO iCnt=1,kProfLayer
                raVTemp(iCnt)=raTTemp(iCnt)
            END DO
            iVTset=1
        END IF
    ! if previous profiles have been read in, check to make sure the
    ! temperature profiles are the same!!!!
        IF ((iGasID <= kGasComp) .AND. (iVTSet > 0)) THEN
            IF ((abs(kLongOrShort) > 1) .AND. (kOuterLoop == 1)) THEN
                write(kStdWarn,*) 'Regular gas : Checking the vertical tempr profile ...'
            ELSE
                write(kStdWarn,*) 'Regular gas : Checking the vertical tempr profile ...'
            ENDIF
            DO iCnt=1,kProfLayer
                rCheckTemp=raTTemp(iCnt)-raVTemp(iCnt)
                IF (abs(rCheckTemp) >= 1.0e-3) THEN
                    write(kStdErr,*)'Warning!!Tempr profiles do not match!!'
                    write(kStdErr,*)'Gas#,layer, gastemp, vertical temp = '
                    write(kStdErr,*)iCount,iCnt,raTTemp(iCnt),raVTemp(iCnt)

                    write(kStdWarn,*)'Warning!!Tempr profiles do not match!!'
                    write(kStdWarn,*)'Gas#,layer, gastemp, vertical temp = '
                    write(kStdWarn,*)iCount,iCnt,raTTemp(iCnt),raVTemp(iCnt)

                END IF
            END DO
        END IF
    END IF

    iDefault = +1
    iMatlabORf77 = -1   !!! use original (pre 2011) f77 uncompression routines
    iMatlabORf77 = +1   !!! use Matlab based (2011)     uncompression routines
    iMatlabORf77 = iaaOverrideDefault(1,4)
    IF (abs(iMatlabORf77) /= 1) THEN
        write(kStdErr,*) 'invalid iMatlabORf77 = ',iMatlabORf77
        CALL DoStop
    END IF
    IF ((iMatlabORf77 /= iDefault) .AND. (kOuterLoop == 1)) THEN
        write(kStdErr,*) 'using iMatlab/f77 = ',iMatlabORf77,' not ',iDefault
        write(kStdWarn,*) 'using iMatlab/f77 = ',iMatlabORf77,' not ',iDefault
    END IF

! then read in the compressed data ... if water do the uncompression
! differently than if the gas is one of the others
    iErr1=0
    IF (iErr <= 0) THEN
        IF ((iGasID == 1) .OR. (iGasID == (kNewGasHi + 1))) THEN
            IF (iMatlabORf77 == -1) THEN
                CALL water(iGasID,rFileStartFr,iTag,iActualTag, &
                iRefLayer,iL,iU, &
                raTAmt,raRAmt,raTPart,raRPartPress,raTTemp,raRTemp, &
                iErr1,iDoDQ,pProf,iProfileLayers, &
                daaDQ,daaDT,daaTemp,iSplineType)
            ELSEIF (iMatlabORf77 == +1) THEN
                CALL xwater(iGasID,rFileStartFr,iTag,iActualTag, &
                iRefLayer,iL,iU, &
                raTAmt,raRAmt,raTPart,raRPartPress,raTTemp,raRTemp, &
                iErr1,iDoDQ,pProf,iProfileLayers, &
                daaDQ,daaDT,daaTemp,iSplineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)
            END IF
        ELSE
            IF (iMatlabORf77 == -1) THEN
                CALL othergases(iGasID,rFileStartFr,iTag,iActualTag, &
                iRefLayer,iL,iU, &
                raTAmt,raRAmt,raTTemp,raRTemp, &
                iErr1,iDoDQ,pProf,iProfileLayers, &
                daaDQ,daaDT,daaTemp,iSplineType)
            ELSEIF (iMatlabORf77 == +1) THEN
                CALL xothergases(iGasID,rFileStartFr,iTag,iActualTag, &
                iRefLayer,iL,iU, &
                raTAmt,raRAmt,raTTemp,raRTemp, &
                iErr1,iDoDQ,pProf,iProfileLayers, &
                daaDQ,daaDT,daaTemp,iSplineType, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)
            END IF
        END IF
        IF (iErr1 <= 0) THEN
            WRITE(kStdWarn,1005) iGasID
            Write(kStdWarn,*) '     '
        ELSE
            iErr=1
            WRITE(kStdWarn,1006) iGasID
            Write(kStdWarn,*) '     '
        END IF
    END IF

    1005 FORMAT('Successfully read in k-comp data for GAS ID ',I3)
    1006 FORMAT('Error ... Could not read k-comp data for GAS ID ',I3)


    RETURN
    end SUBROUTINE compressed

!************************************************************************
! this function checks to see if there is NEW data for the  current gas
!   if so, it returns a positive integer that tells the code which spectra
!   dataset to use (from *SPECTRA) ** would be between 1 to 110 **
! this function then checks to see if there is ALT COMPRESSED DATABASE data for the  current gas
!   if so, it returns a positive integer that tells the code which spectra
!   dataset to use (from *SPECTRA) ** would be between 1001 to 1110 **
! else it returns -1
    INTEGER FUNCTION OutsideSpectra(iGasID,iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs, &
    rFileStartFr,rAltMinFr,rAltMaxFr,iTag)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iGasID       tells current gasID
! iNumNewGases tells how many new gases to use
! iaNewGasID   tells the gasIDs of the new gases
! iTag         tells the kcompressed band code is looking at
    INTEGER :: iGasID, iNumNewGases, iaNewGasID(kGasStore),iTag
! iNumAltComprDirs  tells if we need to use alternate compressed gas dirs
! iaAltComprDirs    tells which gases gave compressed data stored in alternate dirs
    INTEGER :: iNumAltComprDirs,iaAltComprDirs(kGasStore)
! these tell the start/stop wavenumbers for  alternate database (and we are looking at chunk rFileStartFr)
    REAL :: rAltMinFr,rAltMaxFr,rFileStartFr

! local vars
    INTEGER :: iI,iJ

    iI = -1

!      print *,iNumNewGases
!      print *,(iaNewGasID(iI),iI=1,iNumNewGases)
!      call dostopmesg('subr OutsideSpectra$')
          
    IF (iNumNewGases > 0) THEN
        iJ = 1
    ! search to see if there is new data!
        10 CONTINUE
        IF (iaNewGasID(iJ) == iGasID) THEN
            iI = iJ
        ELSEIF (iJ  < iNumNewGases) THEN
            iJ = iJ + 1
            GOTO 10
        END IF
        IF (iI > 0) THEN
            write(kStdWarn,*) '>>> found alternate monochromatic SPECTRA for gasID ',iGasID
            IF (iGASID == 2) write(kStdWarn,*) '  >>> gasID = 2, so could simply be NLTE check ...'
        END IF

    !      ELSEIF ((iNumAltComprDirs .GT. 0) .AND. (rFileStartFr+0.05 .GE. rAltMinFr-0.05)
    !     $                                  .AND. (rFileStartFr-0.05 .LE. rAltMaxFr+0.05)) THEN
    ELSEIF ((iNumAltComprDirs > 0) .AND. (rFileStartFr+0.05 >= rAltMinFr-0.05) &
         .AND. (rFileStartFr+kaBlSize(iTag)-0.05 <= rAltMaxFr+0.05)) THEN
        iJ = 1
    ! search to see if there is new data!
        20 CONTINUE
        IF (iaAltComprDirs(iJ) == iGasID) THEN
            iI = iJ
        ELSEIF (iJ  < iNumAltComprDirs) THEN
            iJ = iJ + 1
            GOTO 20
        END IF
        IF (iI > 0) THEN
            write(kStdWarn,*) '>>> found alternate COMPRESSED DATABASE for gasID ',iGasID
        END IF
        IF (iI > 0) iI = iI + 1000
    END IF

!      print *,iNumNewGases,iNumAltComprDirs,rFileStartFr,rAltMinFr,rAltMaxFr,iI

    OutsideSpectra = iI

    RETURN
    end FUNCTION OutsideSpectra

!************************************************************************
! this function checks to see if there is NEW data for the  current chunk
    INTEGER FUNCTION NewDataChunk(iNewIn,iaNewData,iaaNewChunks, &
    rFileStartFr)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNewIn        tells which NewSpectra set this gas corresponds to
! iaNewData     tells how many new data sets were read in for this set
! iaaNewChunks  tells which data chunks were read in for this set
! rFileStartFr is the integer ID of the relevant k-comp file
    INTEGER :: iaNewData(kGasStore),iaaNewChunks(kGasStore,kNumKCompT)
    INTEGER :: iNewIn
    REAL :: rFileStartFr

    INTEGER :: iI,iJ,iK

    iI = -1

    iJ = 1
! search to see if there is new data!
    10 CONTINUE
    IF (iaaNewChunks(iNewIn,iJ) == nint(rFileStartFr)) THEN
        iI = iJ
    ELSEIF (iJ  < iaNewData(iNewIn)) THEN
        iJ = iJ + 1
        GOTO 10
    END IF

    NewDataChunk = iI

    RETURN
    end FUNCTION NewDataChunk

!************************************************************************
! this subroutine reads in the new data
    SUBROUTINE ReadNewData(iCount,iGasID, &
! these parameters are mainly junk, but required for water continuum
    iRefLayer,iL,iU,raTAmt,raTTemp,raTPress,raTPart,iaCont, &
    iTag,iActualTag,raFreq, &
    daaDQ,daaDT,iDoDQ, &
! these parameters are to read in the externally computed spectra
    daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks,df,sf)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iWhichChunk   tells which file number to read in
! iNewIn        tells which NewSpectra set this gas corresponds to
! caaaNewChunks tells the name of the files associated with the chunks
! daaGasAbCoeff has the uncompressed gas absorption coeff
! df is the frequency step, sf is the start freq
! iGasID is the GasID
    INTEGER :: iWhichChunk,iGasID,iNewIn
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)
    REAL :: df,sf
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
! iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
! raFreq    = wavenumber array
! iErr      = errors (mainly associated with file I/O)
! daaDQ      = analytic Jacobian wrt amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! iTag      = which k-comp file chunk to uncompress
! iCount    = which of the iNumGases is being processed
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPart(kProfLayer),raTPress(kProfLayer)
    INTEGER :: iaCont(kMaxGas),iDoDQ,iCount
    REAL :: raFreq(KmaxPts)
    INTEGER :: iL,iU,iErr,iRefLayer,iTag,iActualTag
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iIoun,I,J
    INTEGER :: IDGAS, NPTS, NLAY
    DOUBLE PRECISION :: SFREQ, FSTEP
    CHARACTER(80) :: FNAM

    FNAM=caaaNewChunks(iNewIn,iWhichChunk)

    write(kStdWarn,*) 'In ReadNewData, opening data file for GasID'
    write(kStdWarn,*) iGasID,FNAM

    iIOUN = kCompUnit
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine ReadNewData'
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

    IF (NPTS /= kMaxPts) THEN
        write(kStdErr,*) ' new data has ',NPTS,' freq pts instead of kMaxPts'
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
        READ(iIOUN) (daaGasAbCoeff(J,I),J=1,kMaxPts)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

    write (kStdWarn,*) 'Read in NEW spectra for gasID ',IDGAS

    RETURN
    end SUBROUTINE ReadNewData

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! iGasID = 1 (WATER), rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
! the spline is done wrt partial pressure, while the scaling is done
! wrt partial pressure
    SUBROUTINE water(iGasID,rFileStartFr,iTag,iActualTag,iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPPart,raRPart,raPTemp,raRTemp, &
    iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSPlineType)

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
            CALL GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset, &
            daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx, &
            raPTemp,raRTemp,raPPart,raRPart,iaTsort, &
            iNk,iKm,iKn,iUm,iUn,daaDQ,daaDT,iDoDQ,iGasID,pProf, &
            iProfileLayers,iSPlineType)
        ELSE
            iLowerOrUpper = -1
            CALL GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset, &
            daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,raPTemp, &
            raRTemp,raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
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
            CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
            daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType)
        ELSE
            iLowerOrUpper = -1
            CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx, &
            raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
            pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
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
    end SUBROUTINE water

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! for the COUSIN CO2 files
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE CousinContribution(iGasID, &
    rFileStartFr,iTag,iActualTag,iProfileLayers,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,pProf, &
    rLTEstrength,iNLTEStart,iLowerOrUpper,daaAbsCoeff,iSPlineType)

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
        
! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp),iFr
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    DOUBLE PRECISION :: daaCousin(kMaxPts,kProfLayer)

    IF (iGasID /= 2) THEN
        write(kStdErr,*) 'This subroutine is onlty for CO2!!!'
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
        CALL GetAbsCoeffNOJAC(daaCousin,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
        pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
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
    end SUBROUTINE CousinContribution

!************************************************************************
! this subroutine calls the routines to read in the k-compressed data
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE othergases(iGasID,rFileStartFr,iTag,iActualTag, &
    iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSplineType)

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
        
! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp)
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    INTEGER :: iLowerOrUpper,iJ

    iIOUN = kCompUnit
    CALL CompFileName(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
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
        CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
        daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype)
    ELSE
        iLowerOrUpper = -1
        CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
        pProf,iProfileLayers,iSplineType,iLowerOrUpper)
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
          
!      print *,iGasID,int(rFileStartFr),kAltComprDirs

    IF (iGasID == 2) THEN
        CALL multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)
    END IF

    RETURN
    end SUBROUTINE othergases

!************************************************************************
! this subroutine mutiplies the daaGasAbsCoeff by CO2 chi functions
    SUBROUTINE multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
    REAL :: rFileStartFr
! input/output
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

! local vars
    INTEGER :: iCO2Chi,iDefault
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge,WhichGasPosn

    iCO2Chi = 0  !!no chi fixes applied .. with database being
! PARAMETER (kCO2Path = '/asl/data/kcarta/v20.ieee-le/etc.ieee-le/')
! this would be the original linemix from RAL, before AIRS was
! launced in April 2002
    iCO2Chi = 3 !!new after March 2004; fixes 4 + 15 um;
    iCO2Chi = 2 !!default prior to Mar 2004; only fixes 4um; leaves wiggles

    iDefault = 2
    iCO2Chi = 0
    iCO2Chi = 2    !!! DEFAULT
    iCO2Chi = iaaOverrideDefault(1,3)
    IF ((iCO2Chi /= 0) .AND. (iCO2Chi /= 2)) THEN
        write(kStdErr,*) 'invalid iCO2Chi = ',iCO2Chi
        CALL DoStop
    END IF

!      IF (kCO2_UMBCorHARTMAN .EQ. +1) THEN
!        iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
!      ELSEIF (kCO2_UMBCorHARTMAN .EQ. -1) THEN
!        iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing
!      END IF
    IF ((kCO2_UMBCorHARTMAN == +1) .AND. (kAltComprDirs == -1)) THEN
        iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
    ELSEIF ((kCO2_UMBCorHARTMAN == -1) .OR. (kAltComprDirs == +1)) THEN
        iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing, or other databases
    END IF

!!!! iCO2Chi = 2   TESTING

    10 FORMAT(I2,I2,I2)
    IF ((iCO2Chi /= iDefault) .AND. (kAltComprDirs == -1) .AND. (kOuterLoop == 1)) THEN
        write(kStdErr,*) ' CO2 chi fudge iDefault,iCO2Chi = ',iDefault,iCO2Chi,'  kCO2_UMBCorHARTMAN = ',kCO2_UMBCorHARTMAN
    !      ELSEIF ((iCO2Chi .NE. iDefault) .AND. (kAltComprDirs .EQ. +1) .AND. (kOuterLoop .EQ. 1)) THEN
    !        write(kStdErr,*) ' CO2 chi fudge iDefault,iCO2Chi = ',iDefault,iCO2Chi,' but using other user suppl CO2 dir'
    END IF

    IF (iCO2Chi == 2) THEN
    ! this is old; prior to March 2004
    ! notice how we only need to fudge fix 2255,2280 and 2305,2405 chunks here!
        iChiChunks = 4
        iaChiChunks(1) = 2255
        iaChiChunks(2) = 2280
        iaChiChunks(3) = 2380
        iaChiChunks(4) = 2405
    ELSE IF (iCO2Chi == 3) THEN
    ! this is new; after March 2004
    ! notice we fix 15 um and 4 um here
        iChiChunks = 18
        iaChiChunks(1)  =  630
        iaChiChunks(2)  =  655
        iaChiChunks(3)  =  680
        iaChiChunks(4)  =  705
        iaChiChunks(5)  =  730
        iaChiChunks(6)  =  755
        iaChiChunks(7)  = 2180
        iaChiChunks(8)  = 2205
        iaChiChunks(9)  = 2230
        iaChiChunks(10) = 2255
        iaChiChunks(11) = 2280
        iaChiChunks(12) = 2355
        iaChiChunks(13) = 2380
        iaChiChunks(14) = 2405
        iaChiChunks(15) = 2430
        iaChiChunks(16) = 2530
        iaChiChunks(17) = 2555
        iaChiChunks(18) = 2580
    END IF

    IF (iCO2Chi > 0) THEN
        iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
    !        print *,int(rFileStartFr),iCO2Chi,iDoFudge
        IF (iDoFudge > 0) THEN
            CALL co2_4um_fudge(daaAbsCoeff,rFileStartFr, &
            iCO2Chi,iaChiChunks,iChiChunks)
        END IF
    END IF
           
    RETURN
    end SUBROUTINE multiply_co2_chi_functions

!************************************************************************
! baed on LLStrow's analysis of the AERI data, this subroutine further scales
! the CO2 absorption coefficients in the 2380-2405 regions
    SUBROUTINE AERI_LLS_fudge(daaAbsCoeff,rFileStartFr,iWhichGas)

    include '../INCLUDE/kcartaparam.f90'

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! rFileStartFr = which chunk
    INTEGER :: iWhichGas
    REAL :: rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay,iChi
    REAL :: raFreq(kMaxPts),raMystery(kMaxPts)

    iChi = -1
    IF (iWHichGas == 1) THEN
        write(kStdWarn,*) 'need H20 AERI chifunction for ',rFileStartFr
        iChi = +1
    ELSEIF (iWHichGas == 2) THEN
        write(kStdWarn,*) 'need CO2 AERI chifunction for ',rFileStartFr
        iChi = +1
    END IF

    IF (iChi > 0) THEN
        DO iFr=1,kMaxPts
            raFreq(iFr) = rFileStartFr*1.0 + (iFr-1)*0.0025
        END DO
        Call FindMysteryFactor(raFreq,raMystery,iWhichGas)
        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*raMystery(iFr)
            END DO
        END DO
    ELSE
        write(kStdWarn,*) 'do not need CO2 chifunction for chunk ',rFileStartFr
    END IF

    RETURN
    end SUBROUTINE AERI_LLS_fudge
         
!************************************************************************
! this subroutine scales the CO2 absorption coefficients in the 2355,2280
! and 2380,2405 chunks
! look at chi.f for the fudge for 2380,2405 cm-1
    SUBROUTINE Co2_4um_fudge(daaAbsCoeff,rFileStartFr, &
    iCO2Chi,iaChiChunks,iChiChunks)

    include '../INCLUDE/kcartaparam.f90'

! iCO2Chi is for the following :
! the blah.txt  files are the orig,      done in late june/early july 2002
! the blah_a.txt files are refinement 1, done in early aug 2002
! the blah_b.txt files are refinement 2, done in mid nov 2002       iCO2Chi = 2
!   Scott did further refinemnents in 2003 and Jan 2004; see file   iCO2Chi = 3
!   /home/sergio/SPECTRA/CKDLINUX/co2_tune.m

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! rFileStartFr = which chunk
    REAL :: rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iCO2Chi

    INTEGER :: iDoFudge,WhichGasPosn,iI,iJ,iK
    INTEGER :: iFr,iLay,iIOUN,iERR,iChi
    REAL :: raF(kMaxPts),raChi(kMaxPts),m,c,rF,rX,rChi,r1,r2,r3,rAmp,rAmp0
    CHARACTER(120) :: FNAME
    CHARACTER(3) :: ca3
    CHARACTER(4) :: ca4

    IF ((iCO2Chi /= 2) .AND. (iCO2Chi /= 3)) THEN
        write(kStdErr,*) 'Illegal type for co2 chi function',iCO2Chi
        CALL DoStop
    END IF

    iChi = -1
    iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
    IF (iDoFudge > 0) THEN
        iChi = +1
        DO iI = 1,120
            fname(iI:iI) = ' '
        END DO
        FNAME = 'co2_4um_fudge_'
        iJ = 1
        11 CONTINUE
        IF ((fname(iJ:iJ) /= ' ') .AND. (iJ < 120)) THEN
            iJ = iJ + 1
            GOTO 11
        END IF
        IF (rFileStartFr < 1000) THEN
            write(ca3,30) nint(rFileStartFr)
            fname(iJ:iJ+2) = ca3(1:3)
            iJ = iJ+3
        ELSEIF (rFileStartFr >= 1000) THEN
            write(ca4,40) nint(rFileStartFr)
            fname(iJ:iJ+3) = ca4(1:4)
            iJ = iJ+4
        END IF
    END IF

    IF (iCO2Chi == 2) THEN
        fname(iJ:iJ+5) = '_b.txt'
    ELSEIF (iCO2Chi == 3) THEN
        fname(iJ:iJ+5) = '_c.txt'
    END IF

    30 FORMAT(I3)
    40 FORMAT(I4)

!      IF (iFileStartFR .EQ. 2255) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2255 chunk ....'
!        FNAME = 'co2_4um_fudge_2255.txt'
!        FNAME = 'co2_4um_fudge_2255_a.txt'
!        FNAME = 'co2_4um_fudge_2255_b.txt'  !!!'a' and 'b are copies
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2280) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2280 chunk ....'
!        FNAME = 'co2_4um_fudge_2280.txt'
!        FNAME = 'co2_4um_fudge_2280_a.txt'
!        FNAME = 'co2_4um_fudge_2280_b.txt'  !!!'a' and 'b are copies
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2380) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2380 chunk ....'
!        FNAME = 'co2_4um_fudge_2380.txt'
!        FNAME = 'co2_4um_fudge_2380_b.txt'
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2405) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2405 chunk ....'
!        FNAME = 'co2_4um_fudge_2405.txt'
!        FNAME = 'co2_4um_fudge_2405_b.txt'
!        iChi = +1
!      END IF
     
    CALL FindChiFileName(fname)

    IF (iChi > 0) THEN
        write(kStdWarn,*) '   Reading in CO2 chifile ',fname
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
        IF (IERR /= 0) THEN
            WRITE(kStdErr,*) 'In subroutine co2_4um_fudge'
            WRITE(kStdErr,1010) IERR, FNAME
            1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
            CALL DoSTOP
        ENDIF
        kTempUnitOpen=1
        READ(iIOUN,*) (raF(iFr),raChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen=-1

    ! to print out the chifcn
    !          DO iFr=1,kMaxPts
    !           print *,raF(iFr),raChi(iFr)
    !           end do

        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*raChi(iFr)
            END DO
        END DO
    ELSE
        write(kStdWarn,*) 'do not need CO2 chifunction for chunk ',rFileStartFr
    END IF

    RETURN
    end SUBROUTINE Co2_4um_fudge

!************************************************************************
! this subroutine reads in a generic Scott Hannon tuning file, and applies
! the tuning to INDIVIDUAL gas optical depths (ie daaAbsCoeff)
! basically copied from SUBROUTINE water_sarta_fudge
    SUBROUTINE generic_sarta_tunmult(iGasID,raFreq,raaAbsCoeff,iSARTAChi)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
          
! input
!   iSARTAChi = -1 for no chi, +1 for AIRS, +2 for IASI, +3 for CrIS
!   iGasID    = gasID
!   raFreq    = wavenumber array
!  raaAbsCoeff = gas OD

    INTEGER :: iGasID,iSARTAChi
    REAL :: rFileStartFr,raFreq(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay,iIOUN,iERR,iMin,iMax,iNpts
    REAL :: raF(kMaxPts),raChi(kMaxPts),raChi2(kMaxPts)
    REAL :: rF,fixed,water,watercon,o3,co,ch4,nte
    CHARACTER(80) :: Fname,caStr
    REAL :: raFX(kMaxPts),raChiX(kMaxPts)

    REAL :: raBeginBand(100),raEndBand(100)
    INTEGER :: iBand,iX,iaOK(kMaxPts)

! see /asl/packages/sartaV108_PGEv6/Src_AIRS_PGEv6/tunmlt_df.f
!    The tuning multiplier file must consist of MXCHAN lines of data
!    sorted (in ascending order) by channel ID and containing the
!    following 9 columns:
!       1    2   3    4     5    6   7   8    9
!       ID RJUNK XF XH2OL XH2OC XO3 XCO XCH4 XNTE

    Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.txt'
    IF (iSARTAChi == +1) THEN
    !! /asl/packages/sartaV108_PGEv6/Src_AIRS_PGEv6/tunmlt_df.f
        Fname = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_PGEv6.txt'
    ELSEIF (iSARTAChi == +2) THEN
    !! /asl/packages/sartaV108/Src_rtpV201_pclsam_slabcloud_hg3/incFTC_iasi_may09_wcon_nte_swch4.f
        Fname = '/asl/data/sarta_database/Data_IASI_may09/Coef/tunmlt_wcon_nte.txt'
    ELSE
        write(kStdErr,*) 'only doing SARTA AIRS tunings'
        CALL DoStop
    END IF

    iIOUN = kTempUnit

    rFileStartFr = raFreq(1)

    1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)

    IF ((rFileStartFr >= 605.0) .AND. (rFileStartFr <= 2805.0)) THEN
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
        IF (IERR /= 0) THEN
            WRITE(kStdErr,*) 'In subroutine generic_sarta_tunmult'
            WRITE(kStdErr,1010) IERR, FNAME
            CALL DoSTOP
        ENDIF

        iNpts = 0
        kTempUnitOpen=1
        10 CONTINUE
        READ(iIOUN,80,ERR=20) caStr
        IF (caStr(1:1) == '%') THEN
            GOTO 10
        ELSE
            iNpts = iNpts + 1
            read(caStr,*) iFr,rF,fixed,water,watercon,o3,co,ch4,nte
            raF(iNpts)   = rF
            raChi(iNpts) = 1.0
            IF ((iGasID == 1) .OR. (iGasID == 103)) THEN
                raChi(iNpts) = water
            !          ELSEIF (iGasID .EQ. 2) THEN
            !            raChi(iNpts) = co2
            ELSEIF (iGasID == 3) THEN
                raChi(iNpts) = o3
            ELSEIF (iGasID == 5) THEN
                raChi(iNpts) = co
            ELSEIF (iGasID == 6) THEN
                raChi(iNpts) = ch4
            ELSEIF (iGasID <= 100) THEN
                raChi(iNpts) = fixed
            END IF
        END IF
        GOTO 10

        20 CONTINUE
        CLOSE(iIOUN)
        kTempUnitOpen = -1

        IF (iSARTAChi == +1) THEN
        ! add on two boundary points to cover the huge 1600-2050 cm-1 gap
        ! 2696  1614.607  1.00000  1.13000  1.00000  1.00000  1.00000  1.00000  1.00000
        ! 2697  2165.877  0.99000  1.00000  1.53116  1.00000  1.00000  1.00000  1.00000
            iNpts = iNpts + 1
            raF(iNpts) = 1614.607 + 0.1
            raChi(iNpts) = 1.0

            iNpts = iNpts + 1
            raF(iNpts) = 2165.877 - 0.1
            raChi(iNpts) = 1.0
        END IF

    !x        call r_sort_spl(raF,raChi,2378,raFreq,raChi2,kMaxPts)
        call r_sort_spl(raF,raChi,iNpts,raFreq,raChi2,kMaxPts)

    !        !! now check simple boundaries, between [raF(1) and raF(end)]
    !        IF (raFreq(1) .LT. raF(1)) THEN
    !          iMin = 1
    !          IF (raFreq(kMaxPts) .LT. raF(1)) THEN
    !            iMax = kMaxPts
    !          ELSE
    !            iMax = int((raF(1)-raFreq(1))/0.0025)
    !          END IF
    !          DO iFr=iMin,iMax
    !            raChi2(iFr) = 1.0
    !          END DO
    !        END IF
    !        IF (raFreq(kMaxPts) .GT. raF(2378)) THEN
    !          iMax = kMaxPts
    !          IF (raFreq(1) .GT. raF(2378)) THEN
    !            iMin = 1
    !          ELSE
    !            iMin = int((raF(2378)-raFreq(1))/0.0025)
    !          END IF
    !          DO iFr=iMin,iMax
    !            raChi2(iFr) = 1.0
    !          END DO
    !        END IF

    !! now check more complicated boundaries
        CALL DoSort2Real(raF,raChi,iNpts,1)
        iBand = 0
        iX = 1
        30 CONTINUE
        IF (abs(raChi(iX)-1.0) >= 1.0e-5) THEN
        !! found  the start of a "band" at which to multiply ODs
            iBand = iBand + 1
            raBeginBand(iBand) = raF(iX)
            iX = iX + 1
            40 CONTINUE
        !! now look for the band end
            IF ((abs(raChi(iX)-1.0) >= 1.0e-5) .AND. (iX < iNpts)) THEN
                iX = iX + 1
                GOTO 40
            ELSE
                raEndBand(iBand) = raF(iX)
            END IF
        END IF
        IF (iX < iNpts) THEN
        !! did not find start of new band, so go to next point
            iX = iX + 1
            GOTO 30
        END IF

    !! now adjust raChi2 so that anything outside the bands gets set to 1
        DO iFr = 1,kMaxPts
            iaOK(iFr) = -1
        END DO
        DO iX = 1,iBand
            DO iFr = 1,kMaxPts
                IF ((raFreq(iFr) >= raBeginBand(iX)) .AND. (raFreq(iFr) <= raEndBand(iX))) THEN
                    iaOK(iFr) = +1
                END IF
            END DO
        END DO
        DO iFr = 1,kMaxPts
            IF (iaOK(iFr) == -1) raChi2(iFr) = 1.0
        END DO

    !        IF (iGasID .LE. 6) THEN
    !          print *,iGasID,iBand
    !          print *,(raBeginBand(iX),iX=1,iBand)
    !          print *,(raEndBand(iX),iX=1,iBand)
    !        END IF

        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                raaAbsCoeff(iFr,iLay) = raaAbsCoeff(iFr,iLay)*max(0.0,raChi2(iFr))
            END DO
        END DO

    END IF

    80 FORMAT(A80)

    RETURN
    end SUBROUTINE generic_sarta_tunmult
         
!************************************************************************
! this subroutine scales water absorption coefficients in the >= 2380 chunks
    SUBROUTINE water_4um_fudge(daaAbsCoeff,rFileStartFr)

    include '../INCLUDE/kcartaparam.f90'

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! rFileStartFr = which chunk
    REAL :: rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay,iIOUN,iERR,iChi,i1,i2,iFileEnd
    REAL :: raF(kMaxPts),raChi(kMaxPts)
    DOUBLE PRECISION :: dX,daX(kMaxPts),dSlope,xA,x

    iChi = -1
    IF (abs(rFileStartFR-2830.0) <= 0.0001) THEN
        write(kStdWarn,*) 'need water chifunction for ',nint(rFileStartFR),' &
        chunk ....'
        iChi = +2
    ELSEIF (rFileStartFR >= 2405.0) THEN
        write(kStdWarn,*) 'need water chifunction for ',nint(rFileStartFR),' &
        chunk ....'
        iChi = +1
    END IF

    IF (iChi == 1) THEN
    !!! for all freqs, for all layers, multiply lines by 0.93
        dX = 0.93d0
        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*dX
            END DO
        END DO
    ELSEIF (iChi == 2) THEN
    !!! for all freqs < 2392, multiply lines by 1.00
        dX = 1.0d0
        iFileEnd = 2405
        i1 = nint((2392.0-2380)*10000.0/25.0)
        DO iFr = 1,i1
            daX(iFr) = dX
        END DO

    !!! for all freqs > 2393, multiply lines by 0.93
        dX = 0.93d0
        iFileEnd = 2405
        i2 = nint((2393.0-2380)*10000.0/25.0)
        DO iFr = i2,kMaxPts
            daX(iFr) = dX
        END DO

    !!!smoothly join the mulitplication factors
        dSlope = (dX-1.0d0)/((2393-2392)*1.0d0)
        xA = 2392.0d0
        DO iFr = i1+1,i2-1
            x = 2380.0d0 + (iFr-1)*0.0025*1.0d0
            daX(iFr) = dSlope*(x-xA) + 1.0d0
        END DO

        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*daX(iFr)
            END DO
        END DO

    ELSE
        write(kStdWarn,*) 'do not need H2O chifunction for chunk ',rFileStartFr
    END IF

    RETURN
    end SUBROUTINE water_4um_fudge
         
!************************************************************************
! this subroutine scales water absorption coefficients in the >= 2380 chunks
    SUBROUTINE water_sarta_fudge(daaAbsCoeff,rFileStartFr)

    include '../INCLUDE/kcartaparam.f90'

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! rFileStartFr = which chunk
    REAL :: rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay,iIOUN,iERR,iChi
    REAL :: raFreq(kMaxPts)
    REAL :: raF(kMaxPts),raChi(kMaxPts),raChi2(kMaxPts)
    REAL :: rF,fixed,water,watercon,o3,co,ch4,co2
    CHARACTER(80) :: Fname,caStr

    iChi = -1
    IF (rFileStartFR >= 1380.0 .AND. rFileStartFR <= 1680.0) THEN
        write(kStdWarn,*) 'need water_sarta_fudge for ',nint(rFileStartFR),' &
        chunk ....'
        write(kStdWarn,*)' h20_sarta_fudge for ',nint(rFileStartFR),' chunk ..'
        iChi = +1
    END IF

    IF (iChi == 1) THEN
        Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.txt'
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
        IF (IERR /= 0) THEN
            WRITE(kStdErr,*) 'In subroutine water_sarta_fudge'
            WRITE(kStdErr,1010) IERR, FNAME
            1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
            CALL DoSTOP
        ENDIF
        kTempUnitOpen=1
        10 CONTINUE
        READ(iIOUN,80) caStr
        IF (caStr(1:1) == '%') THEN
            GOTO 10
        ELSE
            read(caStr,*) iFr,rF,fixed,water,watercon,o3,co,ch4,co2
            raF(iFr)   = rF
            raChi(iFr) = water
            IF (iFr < 2378) THEN
                GOTO 10
            ELSE
                GOTO 20
            END IF
        END IF
        20 CONTINUE
        CLOSE(iIOUN)
        kTempUnitOpen = -1

        DO iFr = 1,kMaxPts
            raFreq(iFr) = rFileStartFr*1.0 + (iFr-1)*0.0025
        END DO

        call rspl(raF,raChi,2378,raFreq,raChi2,kMaxPts)

        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*(raChi2(iFr)*1.0d0)
            END DO
        END DO

    ELSE
        write(kStdWarn,*) 'do not need H2O chifunction for chunk ',rFileStartFr
    END IF

    80 FORMAT(A80)

    RETURN
    end SUBROUTINE water_sarta_fudge
         
!************************************************************************
! this subroutine scales the absorption coefficients by looking at the
! amounts in the gas profiles, thereby computing the optical depth
! compute optical depth = gas amount * abs coeff
    SUBROUTINE AmtScale(daaAbsCoeff,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! raP/RAmt    = actual/reference gas amount profiles
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay
    REAL :: rScale
    DOUBLE PRECISION :: dZero

    dZero=DBLE(0.0)

    DO iLay=1,kProfLayer
        rScale=raPAmt(iLay)
        DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=max(daaAbsCoeff(iFr,iLay),dZero)*rScale
        END DO
    END DO

    RETURN
    end SUBROUTINE AmtScale
         
!************************************************************************
! this subroutine finishes the calculation of d/dq(abscoeff)
    SUBROUTINE FinalAmtDeriv(daaDQ,iType)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaDT   = d/dT
! iType   = compression type
    INTEGER :: iType
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iLay,iFr

    IF (iType == 2) THEN
        DO iLay = 1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
                daaDQ(iFr,iLay)=(daaDQ(iFr,iLay)**4)
            END DO
        END DO
    ELSE
        DO iLay = 1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
                daaDQ(iFr,iLay)=(daaDQ(iFr,iLay))
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE FinalAmtDeriv

!************************************************************************
! this subroutine finishes the computation of d/dT (absCoeff)
    SUBROUTINE FinalTempDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iKtype      = compression type (power = 1 or 4)
! daaAbsCoeff = uncompressed gas abs coeffs
! daaDT       = d/dT
! ra(P/R)Amt  = reference/actual gas amounts
    INTEGER :: iKtype
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDA(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iLay
    REAL :: rScale

    IF (iKtype == 1) THEN
        DO iLay=1,kProfLayer
            rScale = raPAmt(iLay)
            DO iFr=1,kMaxPts
                daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale
            END DO
        END DO
    ELSE
    ! remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
    ! what comes from the spline interpolation routines is d/dT (K(v,T))^(1/4)
    ! or in other words, we get = (1/4) K(v,T)^(-3/4) dK/dT = daaDA
    ! we need d/dT(optical depth) = d/dT q(actual) K(v,T) = q(actual) d/dT K(v,T)
    ! so we get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3
        DO iLay=1,kProfLayer
            rScale = raPAmt(iLay)
            DO iFr=1,kMaxPts
                daaDA(iFr,iLay)=daaDA(iFr,iLay)*rscale* &
                &                       4.0*(daaAbsCoeff(iFr,iLay)**3)
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE FinalTempDeriv

!************************************************************************
! this subroutine finishes the computation of d/dq (absCoeff) for water
    SUBROUTINE FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iKtype      = compression type (power = 1 or 4)
! daaAbsCoeff = uncompressed gas abs coeffs
! daaDT       = d/dT
! ra(P/R)Amt  = reference/actual gas amounts
    INTEGER :: iKtype
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDA(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iLay
    REAL :: rScale

! similar to FinalTempDeriv above
! remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
! what comes from the spline interpolation routines is d/dq (K(v,T))^(1/4)
! or Zget = (1/4) K(v,T)^(-3/4) dK/dq = daaDA
! we need d/dq(optical depth) = d/dq q(actual) K(v,T) = q(actual) d/dq K(v,T)
!                                                        + K(v,T)
! get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3 +  daaAbsCoeff^4
    DO iLay=1,kProfLayer
        rScale=raPAmt(iLay)
        DO iFr=1,kMaxPts
            daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale* &
            &                       4.0*(daaAbsCoeff(iFr,iLay)**3) &
            + (daaAbsCoeff(iFr,iLay)**4)
        END DO
    END DO

    RETURN
    end SUBROUTINE FinalWaterAmtDeriv

!************************************************************************
! this subroutine raises the compressed matrix elements to the 4th power
    SUBROUTINE RaisePower(daaAbsCoeff)

    IMPLICIT NONE
          
    include '../INCLUDE/kcartaparam.f90'

! daaGasAbsCoeff = uncompressed, scaled gas abs coeffs
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFr

    DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=(daaAbsCoeff(iFr,iLay))**4
        END DO
    END DO

    RETURN
    end SUBROUTINE RaisePower

!************************************************************************
! CCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
!      Reads a compressed K and U data file for water
    SUBROUTINE RDCOMPWATER(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, &
    NLAY, KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT, &
    KX1, KX2, KX3, KX4, KX5, UX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

!      Calling parameters
! fnam     = name of relevant file that contains compressed database
! iIOUN    = file IOUNIT number
! IDGAS    = GAS ID
! sfreq    = start freq
! fstep    = frequency point spacing
! npts     = number of freq points === kMaxPts
! nlay     = number of layers == kProfLayer
! ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
! nk       = number of singular vectors (<= kMaxK)
! kt       = number of temperature profiles (=11)
! kn       = number of layers == kMaxLayers
! um       = number of freq points == kMaxPts (in UX)
! un       = number of singular vectors (in UX)
! KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
! UX       = uncompression matrix
! ITO,ITSORT = base points to do temperature interpolation
    CHARACTER(120) :: FNAM
    INTEGER :: iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0, &
    ITSORT(kMaxTemp)
    DOUBLE PRECISION :: SFREQ,FSTEP,TOFF(kMaxTemp), &
    KX1(kMaxK,kMaxTemp,kMaxLayer), &
    KX2(kMaxK,kMaxTemp,kMaxLayer), &
    KX3(kMaxK,kMaxTemp,kMaxLayer), &
    KX4(kMaxK,kMaxTemp,kMaxLayer), &
    KX5(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

!     Local variables
    INTEGER :: IERR,I,J,K,RESORT,IHOLD
    REAL :: rTemp

    1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMPWATER'
        WRITE(kStdErr,1010) IERR, FNAM
        CALL DoSTOP
    ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, &
    UM, UN

    1110 FORMAT('Error! Compressed data array dimension exceeds ', &
    'max size')
    1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
!     Make sure the array sizes are <= to the declared sizes
    IF (NK > kMaxK) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1120) NK, kMaxK
        CALL DoSTOP
    ENDIF

    1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
    IF (KT > kMaxTemp) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1130) KT, kMaxTemp
        CALL DoSTOP
    ENDIF

    1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
    IF (KN > kMaxLayer) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1140) KN, kMaxLayer
        CALL DoSTOP
    ENDIF

    1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
    IF (UM > kMaxPts) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1150) UM, kMaxPts
        CALL DoSTOP
    ENDIF

    1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
    IF (UN > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
        CALL DoSTOP
    ENDIF

!     Read in the temperature offsets
    READ(iIOUN) (TOFF(I),I=1,KT)

! check to make sure the offsets differ by kTempOffSet_database
    DO I=1,KT-1
        rTemp = abs(TOFF(I)-TOFF(I+1))
        rTemp = abs(rTemp - kTempOffSet_database)
        IF (rTemp > 0.001) THEN
            write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
            write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
            write(kStdErr,*) (TOFF(J),J=1,KT)
            write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
            write(kStdErr,*) I,kTempOffSet_database
            CALL DoStop
        END IF
    ENDDO

!      Find which of the offsets is 0 (error if none).
    IT0=0
    DO I=1,KT
        IF (TOFF(I) == 0.0) IT0=I
        ITSORT(I)=I
    ENDDO
    IF (IT0 == 0) THEN
        WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
        1180 FORMAT('ERROR! One of the temperature offsets must be 0',/, &
        'offsets =',20(' ',F5.1))
        CALL DoSTOP
    ENDIF

!     Sort the indices of the temperature offsets in ascending order
    RESORT=1
    10 IF (RESORT == 1) THEN
        RESORT=0
        DO I=1,KT-1
            IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
                IHOLD=ITSORT(I)
                ITSORT(I)=ITSORT(I+1)
                ITSORT(I+1)=IHOLD
                RESORT=1
            ENDIF
        ENDDO
        GOTO 10
    ENDIF

!      Read in the five K matrices
! for the old mat2for files READ(iIOUN) ((KX1(I,J,K),J=1,KT),K=1,KN)
! for the WATER mat2for files, READ(iIOUN) ((KX1(I,J,K),K=1,KN),J=1,KT)
    DO I=1,NK
        READ(iIOUN) ((KX1(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
        READ(iIOUN) ((KX2(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
        READ(iIOUN) ((KX3(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
        READ(iIOUN) ((KX4(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
        READ(iIOUN) ((KX5(I,J,K),K=1,KN),J=1,KT)
    ENDDO

!     Read in the U matrix
    DO I=1,NK
        READ(iIOUN) (UX(J,I),J=1,NPTS)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

    RETURN
    END SUBROUTINE RDCOMPWATER

!************************************************************************
! CCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
!      Reads a compressed K and U data file for gases other than water
    SUBROUTINE RDCOMP(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, NLAY, &
    KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT, &
    KX, UX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! fnam     = name of relevant file that contains compressed database
! iiIOUN    = file IOUNIT number
! IDGAS    = GAS ID
! sfreq    = start freq
! fstep    = frequency point spacing
! npts     = number of freq points === kMaxPts
! nlay     = number of layers == kMaxLayer
! ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
! nk       = number of singular vectors (<= kMaxK)
! kt       = number of temperature profiles (=11)
! kn       = number of layers == kMaxLayers
! um       = number of freq points == kMaxPts (in UX)
! un       = number of singular vectors (in UX)
! c un should be equal to nk
! KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
! UX       = uncompression matrix
! ITO,ITSORT = base points to do temperature interpolation
!      Calling parameters
    CHARACTER(120) :: FNAM
    INTEGER :: iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0, &
    ITSORT(kMaxTemp)
    DOUBLE PRECISION :: SFREQ,FSTEP,TOFF(kMaxTemp), &
    KX(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

!      Local variables
    INTEGER :: IERR,I,J,K,RESORT,IHOLD,iDebugMatlab
    REAL :: rTemp
!-----------------------------------------------------------------------

    1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMP'
        WRITE(KSTDERR,1010) IERR, FNAM
        CALL DoSTOP
    ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, &
    UM, UN

    1110 FORMAT('Error! Compressed data array dimension exceeds ', &
    'max size')
    1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
!     Make sure the array sizes are <= to the declared sizes
    IF (NK > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1120) NK, kMaxK
        CALL DoSTOP
    ENDIF

    1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
    IF (KT > kMaxTemp) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1130) KT, kMaxTemp
        CALL DoSTOP
    ENDIF

    1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
    IF (KN > kMaxLayer) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1140) KN, kMaxLayer
        CALL DoSTOP
    ENDIF

    1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
    IF (UM > kMaxPts) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1150) UM, kMaxPts
        CALL DoSTOP
    ENDIF

    1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
    IF (UN > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
        CALL DoSTOP
    ENDIF

!      Read in the temperature offsets
    READ(iIOUN) (TOFF(I),I=1,KT)

! check to make sure the offsets differ by kTempOffSet_database
    DO I=1,KT-1
        rTemp = abs(TOFF(I)-TOFF(I+1))
        rTemp = abs(rTemp - kTempOffSet_database)
        IF (rTemp > 0.001) THEN
            write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
            write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
            write(kStdErr,*) (TOFF(J),J=1,KT)
            write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
            write(kStdErr,*) I,kTempOffSet_database
            CALL DoStop
        END IF
    ENDDO

!      Find which of the offsets is 0 (error if none).
    1180 FORMAT('ERROR! One of the temperature offsets must be 0',/, &
    'offsets =',20(' ',F5.1))
    IT0=0
    DO I=1,KT
        IF (TOFF(I) == 0.0) IT0=I
        ITSORT(I)=I
    ENDDO
    IF (IT0 == 0) THEN
        WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
        CALL DoSTOP
    ENDIF

!     Sort the indices of the temperature offsets in ascending order
    RESORT=1
    10 IF (RESORT == 1) THEN
        RESORT=0
        DO I=1,KT-1
            IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
                IHOLD=ITSORT(I)
                ITSORT(I)=ITSORT(I+1)
                ITSORT(I+1)=IHOLD
                RESORT=1
            ENDIF
        ENDDO
        GOTO 10
    ENDIF

    iDebugMatlab = +1   !!! do     debug using print statements
    iDebugMatlab = -1   !!! do not debug using print statements
    IF (iDebugMatlab > 0) THEN
        print *,idgas,nk,kn,kt,ktype
    !! kx is a matrix(nk,kt,kn) = matrix(nk,11,100)
    !!   where nk = number of basis vectors
    !! this is matrix "kcomp" in the cg5v4250.mat files kcomp(nk,100,11)
    !! where the code is cg"IDGAS"v"FREQ".mat -- see abscmp/ReadmeVIS
    END IF
    DO I=1,NK
        READ(iIOUN) ((KX(I,J,K),K=1,KN),J=1,KT)
    ENDDO

!     Read in the U matrix
!! this is matrix "B" in the cg5v4250.mat files  B(10000,nk)
    DO I=1,NK
        READ(iIOUN) (UX(J,I),J=1,NPTS)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

!      iDebugMatlab = 1
    IF (iDebugMatlab > 0) THEN
        IF (idgas == 2) THEN
            print *,(UX(J,1),J=1,5)     !! should equal B(1:5,1)
            print *,(KX(1,J,6),J=1,6)   !! should equal kcomp(1,6,1:6)
        END IF
    END IF

    RETURN
    END SUBROUTINE RDCOMP

!************************************************************************
END MODULE kcoeffMAIN
