! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeffMAIN

use basic_common
use s_misc
use kcoeff_common
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
    iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr, &
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
    REAL :: rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)
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
    INTEGER :: iNewIn,iWhichChunk
    INTEGER :: ix

    kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath

    iNewIn = -1
    iNewIn = OutsideSpectra(iaGases(iGas),iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs, &
    rFileStartFr,rAltMinFr,rAltMaxFr,iTag)
    IF (iNewIn < 0) THEN
    ! use kCompressed Database w/o worrying
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
        iWhichChunk = NewDataChunk(iNewIn,iaNewData,iaaNewChunks, &
        rFileStartFr)
        IF (iWhichChunk > 0) THEN
            ! read in new spectra
	    write(kStdWarn,*) ' Reading in external spectra for iGas =',iGas,' iNewIn = ',iNewIn,' iWhichChunk = ',iWhichChunk
            CALL ReadNewData(iGas,iaGases(iGas),kProfLayer, &
            iL_low,iL_high,raTAmt,raTTemp,raTPress,raTPartPress, &
            iaCont,iTag,iActualTag,raFreq,daaDQ,daaDT,iDoDQ, &
            daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks, &
            kaFrStep(iTag),rFileStartFr*1.00000)
        ELSE
        ! use kCompressed Database w/o worrying
            kAltComprDirs = -1  !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
	    write(kStdWarn,*) ' still using usual spectra for iGas =',iGas,' iNewIn = ',iNewIn,' iWhichChunk = ',iWhichChunk	    
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
        write(kStdWarn,*) ' Reading in alternate compressed database for iGas =',iGas,' iNewIn = ',iNewIn
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
        iNewIN-1000,iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr)
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
    iNewIN,iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr)

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
    REAL :: rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)

! local variables
    INTEGER :: iFr,iLay

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

      IF (abs(raAltComprDirsScale(iNewIN)-1.0) .GE. 1.0e-8) THEN
        write(kStdWarn,*) '  scale factor for ',iGasID,' is ',raAltComprDirsScale(iNewIN)
        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            daaTemp(iFr,iLay) = daaTemp(iFr,iLay) * raAltComprDirsScale(iNewIN)
          END DO
        END DO
      END IF												    
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

            IF (abs(raAltComprDirsScale(iNewIN)-1.0) .GE. 1.0e-8) THEN
              write(kStdWarn,*) '  scale factor for ',iGasID,' is ',raAltComprDirsScale(iNewIN)
              DO iLay = 1,kProfLayer
                DO iFr = 1,kMaxPts
                  daaTemp(iFr,iLay) = daaTemp(iFr,iLay) * raAltComprDirsScale(iNewIN)
                END DO
              END DO
            END IF												    

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
 120 FORMAT(A80)
    IF (iFileGasID /= iGasID) THEN
      IF ((iFileGasID == 110) .AND. (iGasID == 1)) THEN
        write(kStdWarn,*) 'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
        write(kStdErr,*)  'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
      ELSE
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF
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
      IF ((iFileGasID == 110) .AND. (iGasID == 1)) THEN
        write(kStdWarn,*) 'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
        write(kStdErr,*)  'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
      ELSE
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF
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
      IF ((iFileGasID == 110) .AND. (iGasID == 1)) THEN
        write(kStdWarn,*) 'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
        write(kStdErr,*)  'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
      ELSE
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF
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
    Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_ones.txt'
    IF (iSARTAChi == +1) THEN
      !! /asl/packages/sartaV108_PGEv6/Src_AIRS_PGEv6/tunmlt_df.f    
      Fname = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_PGEv6.txt'
      Fname = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt' 
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
        READ(iIOUN,80,ERR=20,IOSTAT=IERR) caStr
        IF ((caStr(1:1) == '%') .AND. (IERR == 0)) THEN
          GOTO 10
        ELSEIF (IERR == 0) THEN
          iNpts = iNpts + 1
          read(caStr,*) iFr,rF,fixed,water,watercon,o3,co,ch4,nte
          raF(iNpts)   = rF
          raChi(iNpts) = 1.0
          IF ((iGasID == 1) .OR. (iGasID == 103)) THEN
            raChi(iNpts) = water
          !ELSEIF (iGasID .EQ. 2) THEN
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
	ELSEIF (iERR /= 0) THEN
	  GOTO 20
        END IF
        GOTO 10

20      CONTINUE
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
END MODULE kcoeffMAIN
