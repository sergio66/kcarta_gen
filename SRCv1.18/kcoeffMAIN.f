c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c********* this file has the main k-compressed routines *****************
c** which include reading in the data, doing the spline interpolations **
c********** and finding the absorption matrix for the relevant gas ******
c************************************************************************

c************************************************************************
c*********** NOTE THESE VARIABLES ARE DOUBLE PRECISION ******************
c************************************************************************
c
c  if kJacobian == 1 then call the spline-type routines
c  if kJacobian == -1 then call OLD routines
c where JAC( indicates routines for calculating jacobians, and NOJAC(
c indicates original routines w/o jacobians
c
c Also, if iGasID > kMaxDQ then we do not need to calculate d/dq, only d/dT
c contribution. Since kMaxDQ >= 1, water d/dq,d/dT always calculated
c************************************************************************
c this is the MAIN routine
      SUBROUTINE UsualLTEUncompress(iGas,iaGases,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,
     $          pProf,iProfileLayers,
     $          raVertTemp,iVertTempSet,
     $          rFileStartFr,iTag,iActualTag,raFreq,iError,iDoDQ,iSplineType,
     $          iNumNewGases,iaNewGasID,caaaNewChunks,iaNewData,iaaNewChunks,
     $          iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr,
     $          daaDQ,daaDT,daaGasAbCoeff,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input vars
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore)
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
      CHARACTER*120 caaaNewChunks(kGasStore,kNumkCompT) 
      INTEGER iaAltComprDirs(kGasStore),iNumAltComprDirs
      CHARACTER*120 caaAltComprDirs(kGasStore) 
      INTEGER iGas,iaGases(kMaxGas)
      INTEGER iProfileLayers,iVertTempSet,iError,iDoDQ
      INTEGER iSplineType,iL_low,iL_high,iTag,iActualTag
      INTEGER iaCONT(kMaxGas)
c these are the individual reference profiles, at kProfLayer layers 
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer) 
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer) 
c these are the user specified layer profiles 
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer) 
      REAL raTPartPress(kProfLayer),raTPress(kProfLayer) 
      REAL pProf(kProfLayer),raVertTemp(kProfLayer),raFreq(kMaxPts)
      REAL rFileStartFr
      REAL rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)

c output vars
c daaDT,daaDQ are the d/dq,d/dT matrices 
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac) 
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac) 
c daaGasAbCoeff has the uncompressed gas absorption coeff 
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer) 

c local vars
      INTEGER iNewIn,OutSideSpectra,NewDataChunk,iWhichChunk
      INTEGER ix

      kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath

      iNewIn = -1
      iNewIn = OutsideSpectra(iaGases(iGas),iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs,
     $                        rFileStartFr,rAltMinFr,rAltMaxFr,iTag)
      IF (iNewIn .LT. 0) THEN
        !use kCompressed Database w/o worrying
        kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
        CALL GasContribution(iGas,iaGases(iGas),kProfLayer,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          pProf,iProfileLayers,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile,
     $          raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,raFreq,
     $          iError,iDoDQ,
     $          daaDQ,daaDT,daaGasAbCoeff,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          END IF

      IF ((iNewIn .GE. 1) .AND. (iNewIn .LT. 1000)) THEN
        !!!! new monochromatic spectra from file
        iWhichChunk=NewDataChunk(iNewIn,iaNewData,iaaNewChunks,
     $                               rFileStartFr)
        IF (iWhichChunk .GT. 0) THEN
          !read in new spectra
          CALL ReadNewData(iGas,iaGases(iGas),kProfLayer,
     $            iL_low,iL_high,raTAmt,raTTemp,raTPress,raTPartPress,
     $            iaCont,iTag,iActualTag,raFreq,daaDQ,daaDT,iDoDQ,
     $            daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks,
     $            kaFrStep(iTag),rFileStartFr*1.00000)
        ELSE
          !use kCompressed Database w/o worrying
          kAltComprDirs = -1  !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath
          CALL GasContribution(iGas,iaGases(iGas),kProfLayer,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          pProf,iProfileLayers,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile,
     $          raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,
     $          raFreq,iError,iDoDQ,
     $          daaDQ,daaDT,daaGasAbCoeff,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
        END IF
      END IF

      IF ((iNewIn .GE. 1001) .AND. (iNewIn .LT. 10000)) THEN
        !!!! use alternate compressed database
        kAltComprDirs = +1 !!overwrite one of kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath        
        CALL GasContributionAlternateDataBase(iGas,iaGases(iGas),kProfLayer,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          pProf,iProfileLayers,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile,
     $          raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,
     $          raFreq,iError,iDoDQ,
     $          daaDQ,daaDT,daaGasAbCoeff,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22,
     $          iNewIN-1000,iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr)
      END IF

c      if (iGas .EQ. 1) then
c        print *, '%kcoeffmain : check me out '
cc        do ix = 1,100
cc          print *, 'kcoeffmain : check me out ',ix,daaDT(4602,ix)
cc          end do
c        do ix = 1,10000
c          print *, ix,daaDT(ix,4)
c         end do
c       end if

      RETURN
      END

c************************************************************************
c this subroutine gets the contribution of the i-th gas to the
c absorbtion coefficients. Using the reference profile, it scales the
c contribution accordingly
c also, the gases are weighted by raaMix (from MIXFIL)
c iGasID is the gas ID, while rFileStartFr identifies the wavenumber block
      SUBROUTINE GasContribution(iCount,iGasID,iRefLayer,
     $  raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,pProf,iProfileLayers,
     $  raAmt,raTemp,raPress,raPartPress,iaCont,caXsecF,
     $  raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,raFreq,iErr,iDoDQ,
     $  daaDQ,daaDT,daaTemp,
     $  iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf     = actual layers (from kLAYERS) avg pressure, in iProfileLayers
c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raAmt(kProfLayer),raTemp(kProfLayer),pProf(kProfLayer)
      REAL raPress(kProfLayer),raPartPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iVTSet,iTag,iProfileLayers
      REAL raVTemp(kProfLayer),raFreq(kMaxPts)
      INTEGER iCount,iGasID,iL,iU,iErr,iRefLayer
      INTEGER iActualTag,iDoDQ,iSplineType
      CHARACTER*120 caXsecF
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer),kFrStep,rFileStartFr
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)

c local variables
      INTEGER iFr,iLay

      iErr = -1

      kAltComprDirs = -1   !! stick to kWaterPath,kCKD_Compr_Path,kWaterIsotopePath,kCO2Path,kCompPath

      IF (((1 .LE. iGasID) .AND. (iGasID .LE. kGasComp)) .OR. 
     $    (iGasID .EQ. kNewGasHi+1)) THEN
        kFrStep = kaFrStep(iTag)
        CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,
     $    raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp,iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
      END IF

      IF ((kGasXsecLo .LE. iGasID) .AND. (iGasID .LE. kGasXsecHi)) THEN
        kFrStep = kaFrStep(iTag)
        IF (KXsecFormat .GT. 0) THEN
          write (kStdWarn,*) ' xsec gas : using kCompressed Database format'
          CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $      raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,
     $      raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $      pProf,iProfileLayers,iDoDQ,
     $      daaDQ,daaDT,daaTemp,iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
        ELSE
          write (kStdWarn,*) ' xsec gas : using old style binary file format'
          CALL CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep,daaTemp,
     $         raVTemp,iVTSet,raFreq,iErr,caXsecF,
     $         raAmt,raTemp,raPress,raPartPress,iaCont,
     $         daaDQ,daaDT,iDoDQ)
        END IF
      END IF

      IF ((kNewGasLo .LE. iGasID) .AND. (iGasID .LE. kNewGasHi)) THEN
        IF (kCKD .GE. 0) THEN           !add on continuum
          kFrStep = kaFrStep(iTag)
          CALL driver_continuum(iCount,iGasID,iRefLayer,iProfileLayers,
     $      raRAmt,raRTemp,raRPress,raRPartPress,
     $      iL,iU,daaTemp,raVTemp,iVTSet,
     $      rFileStartFr,iTag,iActualTag,
     $      raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $                   daaDQ,daaDT,iDoDQ,iSplineType,pProf,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
        ELSE
          write(kStdWarn,*)'even if 101,102 are valid gases, '
          write(kStdWarn,*)'not adding on continuum as kCKD < 0'         
          DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
              daaTemp(iFr,iLay)=0.0d0
            END DO
          END DO
          IF (iDoDQ .GT. -2) THEN
            DO iLay=1,kProfLayer
              DO iFr=1,kMaxPts
                daaDQ(iFr,iLay)=0.0d0
              END DO
            END DO
          END IF
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this subroutine gets the contribution of the i-th gas to the
c absorbtion coefficients. Using the reference profile, it scales the
c contribution accordingly
c also, the gases are weighted by raaMix (from MIXFIL)
c iGasID is the gas ID, while rFileStartFr identifies the wavenumber block
c same as GasContribution except it substitudes COMPRESSED DATABASE
      SUBROUTINE GasContributionAlternateDatabase(iCount,iGasID,iRefLayer,
     $  raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,pProf,iProfileLayers,
     $  raAmt,raTemp,raPress,raPartPress,iaCont,caXsecF,
     $  raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,raFreq,iErr,iDoDQ,
     $  daaDQ,daaDT,daaTemp,
     $  iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22,
     $          iNewIN,iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,
     $          caaAltComprDirs,rAltMinFr,rAltMaxFr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf     = actual layers (from kLAYERS) avg pressure, in iProfileLayers
c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raAmt(kProfLayer),raTemp(kProfLayer),pProf(kProfLayer)
      REAL raPress(kProfLayer),raPartPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iVTSet,iTag,iProfileLayers
      REAL raVTemp(kProfLayer),raFreq(kMaxPts)
      INTEGER iCount,iGasID,iL,iU,iErr,iRefLayer
      INTEGER iActualTag,iDoDQ,iSplineType
      CHARACTER*120 caXsecF
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer),kFrStep,rFileStartFr
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)
c the alt database
      INTEGER iNewIN,iNumAltComprDirs,iaAltComprDirs(kGasStore)
      CHARACTER*120 caaAltComprDirs(kGasStore)
      REAL rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)

c local variables
      INTEGER iFr,iLay,strfind

      iErr = -1

      DO iFr = 1,80
        kcaAltComprDirs(iFr:iFr) = ' '
      END DO
      kcaAltComprDirs = caaAltComprDirs(iNewIN)
      write(kStdWarn,*) '>>> substituting caCompressedDataPath for gasID ',iGasID
      write(kStdWarn,80) kcaAltComprDirs

      IF (iGasID .EQ. 2) THEN
        IF ((strfind(kcaAltComprDirs,'lblrtm') .EQ. 1) .OR. (strfind(kcaAltComprDirs,'LBLRTM') .EQ. 1)) THEN
	  IF (kCO2_UMBCorHARTMAN .GT. 0) THEN
  	    write(kStdWarn,*) ' kCO2_UMBCorHARTMAN is +1 for UMBC CO2 linemix',kCO2_UMBCorHARTMAN
	    write(kStdWarn,*) ' ignore so NO chi fcns, since you are using LBLRTM CO2 database'
	    write(kStdErr,*)  ' ignore kCO2_UMBCorHARTMAN (+1) so NO chi fcns, for LBLRTM CO2 database'
	    !kCO2_UMBCorHARTMAN = -1
	  END IF
	END IF
      END IF
        
      IF ( ((1 .LE. iGasID) .AND. (iGasID .LE. kGasComp)) .OR. (iGasID .EQ. kNewGasHi+1)) THEN
        kFrStep = kaFrStep(iTag)
        CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,
     $    raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp,iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
        IF (abs(raAltComprDirsScale(iNewIN)-1.0) .GE. 1.0e-8) THEN
	  write(kStdWarn,*) '  scale factor for ',iGasID,' is ',raAltComprDirsScale(iNewIN)
	  DO iLay = 1,kProfLayer
	    DO iFr = 1,kMaxPts
              daaTemp(iFr,iLay) = daaTemp(iFr,iLay) * raAltComprDirsScale(iNewIN)
	    END DO
	  END DO
	END IF
      END IF

      IF ((kGasXsecLo .LE. iGasID) .AND. (iGasID .LE. kGasXsecHi)) THEN
        kFrStep = kaFrStep(iTag)
        IF (KXsecFormat .GT. 0) THEN
          write (kStdWarn,*) ' xsec gas : using kCompressed Database format'
          CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $      raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,
     $      raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $      pProf,iProfileLayers,iDoDQ,
     $      daaDQ,daaDT,daaTemp,iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
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

      IF ((kNewGasLo .LE. iGasID) .AND. (iGasID .LE. kNewGasHi)) THEN
        IF (kCKD .GE. 0) THEN           !add on continuum
          kFrStep = kaFrStep(iTag)
          write(kStdErr,*)  'not supported here in GasContributionAlternateDatabase'
          CALL DoStop
        END IF
      END IF

 80   FORMAT(A120)

      RETURN
      END

c************************************************************************
c this subroutine basically computes water continuum
      SUBROUTINE driver_continuum(iCount,iGasID,iRefLayer,iProfileLayers,
     $   raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,daaTemp,raVTemp,iVTSet,
     $   rFileStartFr,iTag,iActualTag,
     $   raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont,
     $                   daaDQ,daaDT,iDoDQ,iSplineType,pProf,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers (from kLAYERS) avg pressure, in iProfileLayers
c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
c daaDQ      = analytic Jacobian wrt amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer),pProf(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer),rFileStartFr
      INTEGER iaCont(kMaxGas),iVTSet,iDoDQ
      INTEGER iTag,iSPlineType,iActualTag,iProfileLayers
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      REAL raVTemp(kProfLayer),raFreq(KmaxPts)
      INTEGER iGasID,iL,iU,iErr,iCount,iRefLayer
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)

      INTEGER iAns
      REAL kFrStep

      iErr=0

      IF (iGasID .EQ. kNewGasLo) THEN
c check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep = kaFrStep(iTag)
        IF ((iAns .GT. 0) .AND. (kCKD .GE. 0)) THEN
          CALL AddContinuum(iGasID,iTag,iActualTag,rFileStartFr,
     $                  iRefLayer,iProfileLayers,raFreq,raTAmt,raTTemp,
     $                        kFrStep,raTPress,raTPart,iL,iU,daaTemp,
     $                        daaDQ,daaDT,iDoDQ,iSPlineType,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,pProf,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          write(kStdWarn,*) 'added self water continuum'
        END IF 
      END IF

      IF (iGasID .EQ. kNewGasHi) THEN
c check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep = kaFrStep(iTag)
        IF ((iAns .GT. 0) .AND. (kCKD .GE. 0)) THEN
          CALL AddContinuum(iGasID,iTag,iActualTag,rFileStartFr,
     $                   iRefLayer,iProfileLayers,raFreq,raTAmt,raTTemp,
     $                        kFrStep,raTPress,raTPart,iL,iU,daaTemp,
     $                        daaDQ,daaDT,iDoDQ,iSPlineType,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,pProf,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          write(kStdWarn,*) 'added foreign water continuum'
        END IF 
      END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the contribution of Gases 1-kGasComp (if present)
      SUBROUTINE compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,rFileStartFr,iTag,iActualTag,
     $    raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp,
     $    iSPlineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c pProf       = actual layers from kLAYERS avg pressure
c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
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
c daaDQ      = analytic Jacobian wrt amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer),pProf(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer),rCheckTemp
      INTEGER iaCont(kMaxGas),iVTSet,iDoDQ,iTag,iSPlineType,iActualTag
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      REAL raVTemp(kProfLayer),raFreq(KmaxPts),rFileStartFr
      INTEGER iGasID,iL,iU,iErr,iCount,iRefLayer,iProfileLayers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)

      INTEGER iErr1,iCnt,iMatlabORf77,iDefault

      iErr1=0    !this sees if data successfully uncompressed

      IF (iErr1 .LE. 0) THEN
c set the vertical temperature profile if iGasID < 51 (first profile read)
        IF ((iGasID .LE. kGasComp) .AND. (iVTSet .LT. 0)) THEN 
          write(kStdWarn,* )'Setting the vertical tempr profile ...'
          DO iCnt=1,kProfLayer
            raVTemp(iCnt)=raTTemp(iCnt)
          END DO
          iVTset=1
        END IF
c if previous profiles have been read in, check to make sure the 
c temperature profiles are the same!!!!
        IF ((iGasID .LE. kGasComp) .AND. (iVTSet .GT. 0)) THEN 
          IF ((abs(kLongOrShort) .GT. 1) .AND. (kOuterLoop .EQ. 1)) THEN
            write(kStdWarn,*) 'Regular gas : Checking the vertical tempr profile ...'
          ELSE
            write(kStdWarn,*) 'Regular gas : Checking the vertical tempr profile ...'
          ENDIF
          DO iCnt=1,kProfLayer
            rCheckTemp=raTTemp(iCnt)-raVTemp(iCnt)
            IF (abs(rCheckTemp) .GE. 1.0e-3) THEN
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
      IF (abs(iMatlabORf77) .NE. 1) THEN
        write(kStdErr,*) 'invalid iMatlabORf77 = ',iMatlabORf77
	CALL DoStop
      END IF
      IF ((iMatlabORf77 .NE. iDefault) .AND. (kOuterLoop .EQ. 1)) THEN
        write(kStdErr,*) 'using iMatlab/f77 = ',iMatlabORf77,' not ',iDefault
        write(kStdWarn,*) 'using iMatlab/f77 = ',iMatlabORf77,' not ',iDefault	
      END IF

c then read in the compressed data ... if water do the uncompression 
c differently than if the gas is one of the others
      iErr1=0
      IF (iErr .LE. 0) THEN
        IF ((iGasID .EQ. 1) .OR. (iGasID .EQ. (kNewGasHi + 1))) THEN
          IF (iMatlabORf77 .EQ. -1) THEN
            CALL water(iGasID,rFileStartFr,iTag,iActualTag,
     $        iRefLayer,iL,iU,
     $        raTAmt,raRAmt,raTPart,raRPartPress,raTTemp,raRTemp,
     $        iErr1,iDoDQ,pProf,iProfileLayers,
     $        daaDQ,daaDT,daaTemp,iSplineType)
          ELSEIF (iMatlabORf77 .EQ. +1) THEN
            CALL xwater(iGasID,rFileStartFr,iTag,iActualTag,
     $        iRefLayer,iL,iU,
     $        raTAmt,raRAmt,raTPart,raRPartPress,raTTemp,raRTemp,
     $        iErr1,iDoDQ,pProf,iProfileLayers,
     $        daaDQ,daaDT,daaTemp,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          END IF
        ELSE 
          IF (iMatlabORf77 .EQ. -1) THEN
            CALL othergases(iGasID,rFileStartFr,iTag,iActualTag,
     $         iRefLayer,iL,iU,
     $         raTAmt,raRAmt,raTTemp,raRTemp,
     $         iErr1,iDoDQ,pProf,iProfileLayers,
     $         daaDQ,daaDT,daaTemp,iSplineType)
          ELSEIF (iMatlabORf77 .EQ. +1) THEN
            CALL xothergases(iGasID,rFileStartFr,iTag,iActualTag,
     $         iRefLayer,iL,iU,
     $         raTAmt,raRAmt,raTTemp,raRTemp,
     $         iErr1,iDoDQ,pProf,iProfileLayers,
     $         daaDQ,daaDT,daaTemp,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          END IF
        END IF
        IF (iErr1 .LE. 0) THEN
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
      END

c************************************************************************
c this function checks to see if there is NEW data for the  current gas
c   if so, it returns a positive integer that tells the code which spectra
c   dataset to use (from *SPECTRA) ** would be between 1 to 110 ** 
c this function then checks to see if there is ALT COMPRESSED DATABASE data for the  current gas
c   if so, it returns a positive integer that tells the code which spectra
c   dataset to use (from *SPECTRA) ** would be between 1001 to 1110 **  
c else it returns -1
      INTEGER FUNCTION OutsideSpectra(iGasID,iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs,
     $                                rFileStartFr,rAltMinFr,rAltMaxFr,iTag) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iGasID       tells current gasID
c iNumNewGases tells how many new gases to use
c iaNewGasID   tells the gasIDs of the new gases
c iTag         tells the kcompressed band code is looking at
      INTEGER iGasID, iNumNewGases, iaNewGasID(kGasStore),iTag
c iNumAltComprDirs  tells if we need to use alternate compressed gas dirs
c iaAltComprDirs    tells which gases gave compressed data stored in alternate dirs
      INTEGER iNumAltComprDirs,iaAltComprDirs(kGasStore)
c these tell the start/stop wavenumbers for  alternate database (and we are looking at chunk rFileStartFr)
      REAL rAltMinFr,rAltMaxFr,rFileStartFr

c local vars
      INTEGER iI,iJ

      iI = -1

c      print *,iNumNewGases
c      print *,(iaNewGasID(iI),iI=1,iNumNewGases)
c      call dostopmesg('subr OutsideSpectra$')
      
      IF (iNumNewGases .GT. 0) THEN
        iJ = 1
c search to see if there is new data!     
 10     CONTINUE
        IF (iaNewGasID(iJ) .EQ. iGasID) THEN
          iI = iJ
        ELSEIF (iJ  .LT. iNumNewGases) THEN
          iJ = iJ + 1
          GOTO 10
        END IF       
        IF (iI .GT. 0) THEN
          write(kStdWarn,*) '>>> found alternate monochromatic SPECTRA for gasID ',iGasID
          IF (iGASID .EQ. 2) write(kStdWarn,*) '  >>> gasID = 2, so could simply be NLTE check ...'
        END IF

c      ELSEIF ((iNumAltComprDirs .GT. 0) .AND. (rFileStartFr+0.05 .GE. rAltMinFr-0.05) 
c     $                                  .AND. (rFileStartFr-0.05 .LE. rAltMaxFr+0.05)) THEN
      ELSEIF ((iNumAltComprDirs .GT. 0) .AND. (rFileStartFr+0.05 .GE. rAltMinFr-0.05) 
     $                                  .AND. (rFileStartFr+kaBlSize(iTag)-0.05 .LE. rAltMaxFr+0.05)) THEN
        iJ = 1
c search to see if there is new data!     
 20     CONTINUE
        IF (iaAltComprDirs(iJ) .EQ. iGasID) THEN
          iI = iJ
        ELSEIF (iJ  .LT. iNumAltComprDirs) THEN
          iJ = iJ + 1
          GOTO 20
        END IF        
        IF (iI .GT. 0) THEN
          write(kStdWarn,*) '>>> found alternate COMPRESSED DATABASE for gasID ',iGasID
        END IF
        IF (iI .GT. 0) iI = iI + 1000
      END IF

c      print *,iNumNewGases,iNumAltComprDirs,rFileStartFr,rAltMinFr,rAltMaxFr,iI

      OutsideSpectra = iI

      RETURN
      END

c************************************************************************
c this function checks to see if there is NEW data for the  current chunk
      INTEGER FUNCTION NewDataChunk(iNewIn,iaNewData,iaaNewChunks,
     $                             rFileStartFr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iNewIn        tells which NewSpectra set this gas corresponds to
c iaNewData     tells how many new data sets were read in for this set
c iaaNewChunks  tells which data chunks were read in for this set
c rFileStartFr is the integer ID of the relevant k-comp file 
      INTEGER iaNewData(kGasStore),iaaNewChunks(kGasStore,kNumKCompT)
      INTEGER iNewIn
      REAL rFileStartFr

      INTEGER iI,iJ,iK

      iI = -1

      iJ = 1
c search to see if there is new data!     
 10   CONTINUE
      IF (iaaNewChunks(iNewIn,iJ) .EQ. nint(rFileStartFr)) THEN
        iI = iJ
      ELSEIF (iJ  .LT. iaNewData(iNewIn)) THEN
        iJ = iJ + 1
        GOTO 10
      END IF        

      NewDataChunk = iI

      RETURN
      END

c************************************************************************
c this subroutine reads in the new data 
      SUBROUTINE ReadNewData(iCount,iGasID,
c these parameters are mainly junk, but required for water continuum
     $ iRefLayer,iL,iU,raTAmt,raTTemp,raTPress,raTPart,iaCont,
     $ iTag,iActualTag,raFreq,
     $ daaDQ,daaDT,iDoDQ,
c these parameters are to read in the externally computed spectra
     $     daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks,df,sf)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iWhichChunk   tells which file number to read in
c iNewIn        tells which NewSpectra set this gas corresponds to
c caaaNewChunks tells the name of the files associated with the chunks 
c daaGasAbCoeff has the uncompressed gas absorption coeff 
c df is the frequency step, sf is the start freq
c iGasID is the GasID
      INTEGER iWhichChunk,iGasID,iNewIn
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer) 
      CHARACTER*120 caaaNewChunks(kGasStore,kNumkCompT) 
      REAL df,sf
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ      = analytic Jacobian wrt amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c iTag      = which k-comp file chunk to uncompress
c iCount    = which of the iNumGases is being processed
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iDoDQ,iCount
      REAL raFreq(KmaxPts)
      INTEGER iL,iU,iErr,iRefLayer,iTag,iActualTag
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iIoun,I,J
      INTEGER IDGAS, NPTS, NLAY
      DOUBLE PRECISION SFREQ, FSTEP
      CHARACTER*120 FNAM

      FNAM=caaaNewChunks(iNewIn,iWhichChunk)

      write(kStdWarn,*) 'In ReadNewData, opening data file for GasID'
      write(kStdWarn,*) iGasID,FNAM

      iIOUN = kCompUnit
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine ReadNewData'
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

      IF (NPTS .NE. kMaxPts) THEN
        write(kStdErr,*) ' new data has ',NPTS,' freq pts instead of kMaxPts'
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
        READ(iIOUN) (daaGasAbCoeff(J,I),J=1,kMaxPts)
      ENDDO

      CLOSE(iIOUN)
      kCompUnitOpen=-1

      write (kStdWarn,*) 'Read in NEW spectra for gasID ',IDGAS

      RETURN
      END

c************************************************************************
c this subroutine calls the routines to read in the k-compressed data
c iGasID = 1 (WATER), rFileStartFr identifies the frequency range
c have to send in ref temp, amount and profile temp,amount
c the spline is done wrt partial pressure, while the scaling is done
c wrt partial pressure
      SUBROUTINE water(iGasID,rFileStartFr,iTag,iActualTag,iProfLayer,iL,iU,
     $      raPAmt,raRAmt,raPPart,raRPart,raPTemp,raRTemp,
     $      iErr,iDoDQ,pProf,iProfileLayers,
     $      daaDQ,daaDT,daaAbsCoeff,iSPlineType)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
c iGasID     = GASID ==1 for water
c rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
c iTag       = current k-comp block of 25 cm-1 that is being processed
c iProfLayer = number of layers in profile === kProfLayer
c iL,iU      = min/max layer number (=1,kMaxlayer)
c daaAbs     = final uncompressed abs coefficient for gas iGasID
c iErr       = errors (mainly associated with file I/O, could be associated
c              with incorrect number of layers in compresse database etc)
c raP/RAmt   = arrays containing actual/reference gas amounts
c raP/RPart  = arrays containing actual/reference gas partial pressures
c raP/RTemp  = arrays containing actual/reference gas temperatures
c daaDQ      = analytic Jacobian wrt amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ      = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      INTEGER iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
      INTEGER iProfileLayers,iSPlineType
      REAL raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)

c local variables associated with uncompressing the water database files
      CHARACTER*120 caFName
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
      INTEGER iT0,iaTsort(kMaxTemp),iLowerOrUpper
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp)
      DOUBLE PRECISION daaaKX1(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX2(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX3(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX4(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX5(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK)
      INTEGER iDefault,iMultiplyHeavyWater

      IF ((iGasID .NE. 1) .AND. (iGasID .NE. kNewGasHi+1)) THEN
        write(kStdErr,*) 'Expecting to read in water profile/database'
        iErr=1
        CALL DoSTOP
      END IF

      iIOUN = kCompUnit
      CALL CompFileName(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
      CALL rdcompwater(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,
     $        iNLay,iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort,
     $        daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUX)

c check that the file has the data for the correct gas
      IF (iFileGasID .NE. iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
 1000   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF

c check that the data file has the right number of layers ===== AIRS layers
      IF (iNLay .NE. kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
 1010   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
      END IF

c kGenln2Water   = self broadening correction for water, using interpolation  
c                  in water partial pressure (+1)  
c                = just do what Genln2 does (which would be the same as the 
c                  uncompresssion for CO2 (-1) 
c      INTEGER kGenln2Water 
c      PARAMETER (kGenln2Water=+1) 

c interpolate compressed data in temperature, and then in partial pressure, 
c to get abs coeff matrix        
      IF (kGenln2Water .GT. 0) THEN
        !ie worry about the self broadening corrections
        !this is pretty good 
        IF (kJacobian .GT. 0) THEN
          CALL GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset,
     $       daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,
     $       raPTemp,raRTemp,raPPart,raRPart,iaTsort,
     $       iNk,iKm,iKn,iUm,iUn,daaDQ,daaDT,iDoDQ,iGasID,pProf,
     $       iProfileLayers,iSPlineType)
        ELSE
          iLowerOrUpper = -1  
          CALL GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset,
     $      daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,raPTemp,
     $      raRTemp,raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
        END IF

c because iKtype=2 do any necessary jacobians calcs HERE!
        IF (kJacobian .GT. 0) THEN
          IF (iDoDQ .GT. 0)  THEN
            IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 20)) THEN
              CALL FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDQ,raPAmt)
            END IF
          END IF
          IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 30) .OR.
     $        (kActualJacs .EQ. 32) .OR.
     $        (kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102)) THEN
            CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
          END IF
        END IF
    
      ELSE      
        !kGenln2Water .LT. 0  ==> do same uncompression as for CO2
        !ie do not worry about the self broadening corrections
        !this is not very good at all!
c interpolate compressed data in temperature, to get abs coeff matrix
        IF (kJacobian .GE. 0) THEN
          CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $      daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType)
        ELSE
          iLowerOrUpper = -1
          CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
        END IF

c because of iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
        IF (kJacobian .GE. 0) THEN
          IF (iDoDQ .GT. 0) THEN
            IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 20)) THEN
              CALL FinalAmtDeriv(daaDQ,iKtype)
            END IF
          END IF
          IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 30) .OR. 
     $        (kActualJacs .EQ. 32) .OR. 
     $        (kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102)) THEN
            CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
          END IF
        END IF
      END IF

c convert absorption coefficient correctly if necessary
      IF (iKtype .eq. 2) THEN
        CALL RaisePower(daaAbsCoeff)
      END IF

c now compute optical depth = gas amount * abs coeff
      CALL AmtScale(daaAbsCoeff,raPAmt)

      RETURN
      END

c************************************************************************
c this subroutine calls the routines to read in the k-compressed data
c for the COUSIN CO2 files
c iGasID tells which gas type, rFileStartFr identifies the frequency range
c have to send in ref temp, amount and profile temp,amount
      SUBROUTINE CousinContribution(iGasID,
     $       rFileStartFr,iTag,iActualTag,iProfileLayers,iL,iU,
     $      raPAmt,raRAmt,raPTemp,raRTemp,pProf,
     $      rLTEstrength,iNLTEStart,iLowerOrUpper,daaAbsCoeff,iSPlineType)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
c iGasID     = GASID
c rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
c iTag       = current k-comp block of 25 cm-1 that is being processed
c iProfLayer = number of layers in profile === kProfLayer
c iL,iU      = min/max layer number (=1,kMaxlayer)
c daaAbs     = final uncompressed abs coefficient for gas iGasID
c iErr       = errors (mainly associated with fOAile I/O, could be associated
c              with incorrect number of layers in compresse database etc)
c raP/RAmt   = arrays containing actual/reference gas amounts
c raP/RPart  = arrays containing actual/reference gas partial pressures
c raP/RTemp  = arrays containing actual/reference gas temperatures
c new stuff
c rLTEStrength = tells the weight ~ 1.1212
c iNLTEStart = tells where to replace LineMix spectra with Cousin spectra
      INTEGER iGasID,iErr,iL,iU,iTag,iActualTag
      INTEGER iProfileLayers,iLowerOrUpper,iSPlineType
      REAL raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      REAL rLTEStrength
      INTEGER iNLTEStart
    
c local variables associated with uncompressing the database
      CHARACTER*120 caFName
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
      INTEGER iT0,iaTsort(kMaxTemp),iFr
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp)
      DOUBLE PRECISION daaaKX(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK)
      DOUBLE PRECISION daaCousin(kMaxPts,kProfLayer)

      IF (iGasID .NE. 2) THEN
        write(kStdErr,*) 'This subroutine is onlty for CO2!!!'
        CALL DoStop
      END IF

      iIOUN = kCompUnit
      CALL CompFileName(-1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
      CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay,
     $              iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort,
     $              daaaKX,daaUX)

c check that the file has the data for the correct gas
      IF (iFileGasID .NE. iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
 1000   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF

c check that the data file has the right number of layers
      IF (iNLay .NE. kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
 1010   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
      END IF

c interpolate compressed data in temperature, to get abs coeff matrix
      IF (kJacobian .GE. 0) THEN
        write(kStdErr,*) 'Cannot do Jacobians and willy nilly replace '
        write(kStdErr,*) 'linemix spectroscopy with cousin spectroscopy'
        CALL DoStop
      ELSE
        iLowerOrUpper = -1    !!!!only use 100 AIRS layers; nuthin above
        CALL GetAbsCoeffNOJAC(daaCousin,daToffset,daaaKx,daaUx,
     $         raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $         pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
      END IF

c convert absorption coefficient correctly if necessary
      IF (iKtype .eq. 2) THEN
        CALL RaisePower(daaCousin)
      END IF

c now compute optical depth = gas amount * abs coeff
      CALL AmtScale(daaCousin,raPAmt)

c now multiply all spectra by scaling factor, and replace daaGasAb as required
      write (kStdWarn,*) 'Replacing LINEMIX spectra with COUSIN spectra'
      write (kStdWarn,*) 'start layer, strength = ',iNLTEStart,
     $     abs(rLTEStrength)
      DO iL = iNLTEStart,kProfLayer
        DO iFr = 1,kMaxPts
          daaAbsCoeff(iFr,iL) = daaCousin(iFr,iL) * abs(rLTEStrength)
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine calls the routines to read in the k-compressed data
c iGasID tells which gas type, rFileStartFr identifies the frequency range
c have to send in ref temp, amount and profile temp,amount
      SUBROUTINE othergases(iGasID,rFileStartFr,iTag,iActualTag,
     $      iProfLayer,iL,iU,
     $      raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers,
     $      daaDQ,daaDT,daaAbsCoeff,iSplineType)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
c iGasID     = GASID
c rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
c iTag       = current k-comp block of 25 cm-1 that is being processed
c iProfLayer = number of layers in profile === kProfLayer
c iL,iU      = min/max layer number (=1,kMaxlayer)
c daaAbs     = final uncompressed abs coefficient for gas iGasID
c iErr       = errors (mainly associated with file I/O, could be associated
c              with incorrect number of layers in compresse database etc)
c raP/RAmt   = arrays containing actual/reference gas amounts
c raP/RPart  = arrays containing actual/reference gas partial pressures
c raP/RTemp  = arrays containing actual/reference gas temperatures
c daaDT      = analytic Jacobian wrt temperature
c daaDQ      = analytic Jacobian wrt amount
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      INTEGER iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
      INTEGER iProfileLayers,iSplineType
      REAL raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
    
c local variables associated with uncompressing the database
      CHARACTER*120 caFName
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
      INTEGER iT0,iaTsort(kMaxTemp)
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp)
      DOUBLE PRECISION daaaKX(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK)
      INTEGER iLowerOrUpper,iJ

      iIOUN = kCompUnit
      CALL CompFileName(+1,iGasID,rFileStartFr,iTag,iActualTag,caFName)
      CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay,
     $              iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort,
     $              daaaKX,daaUX)

c check that the file has the data for the correct gas
      IF (iFileGasID .NE. iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
 1000   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
      END IF

c check that the data file has the right number of layers
      IF (iNLay .NE. kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
 1010   FORMAT('Error! file : ',/,A120,/,
     $         'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
      END IF

c interpolate compressed data in temperature, to get abs coeff matrix
      IF (kJacobian .GE. 0) THEN
        CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $  raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype)
      ELSE
        iLowerOrUpper = -1
        CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $         raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $         pProf,iProfileLayers,iSplineType,iLowerOrUpper)
      END IF

c because of the iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
      IF (kJacobian .GE. 0) THEN
        IF (iDoDQ .GT. 0) THEN
          IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 20)) THEN
            CALL FinalAmtDeriv(daaDQ,iKtype)
          END IF
        END IF
        IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 30) .OR. 
     $      (kActualJacs .EQ. 32) .OR. 
     $      (kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102)) THEN
          CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
        END IF
      END IF

c convert absorption coefficient correctly if necessary
      IF (iKtype .eq. 2) THEN
        CALL RaisePower(daaAbsCoeff)
      END IF

c now compute optical depth = gas amount * abs coeff
      CALL AmtScale(daaAbsCoeff,raPAmt)
      
c      print *,iGasID,int(rFileStartFr),kAltComprDirs

      IF (iGasID .EQ. 2) THEN
        CALL multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine mutiplies the daaGasAbsCoeff by CO2 chi functions
      SUBROUTINE multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL rFileStartFr
c input/output
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

c local vars
      INTEGER iCO2Chi,iDefault
      INTEGER iaChiChunks(kMaxGas),iChiChunks,iDoFudge,WhichGasPosn

      iCO2Chi = 0  !!no chi fixes applied .. with database being
        !!PARAMETER (kCO2Path = '/asl/data/kcarta/v20.ieee-le/etc.ieee-le/') 
        !!this would be the original linemix from RAL, before AIRS was
        !!launched in April 2002
      iCO2Chi = 3 !!new after March 2004; fixes 4 + 15 um;
      iCO2Chi = 2 !!default prior to Mar 2004; only fixes 4um; leaves wiggles

      iDefault = 2
      iCO2Chi = 0
      iCO2Chi = 2    !!! DEFAULT
      iCO2Chi = iaaOverrideDefault(1,3)      
      IF ((iCO2Chi .NE. 0) .AND. (iCO2Chi .NE. 2)) THEN
        write(kStdErr,*) 'invalid iCO2Chi = ',iCO2Chi
	CALL DoStop
      END IF

c      IF (kCO2_UMBCorHARTMAN .EQ. +1) THEN
c        iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
c      ELSEIF (kCO2_UMBCorHARTMAN .EQ. -1) THEN
c        iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing
c      END IF
      IF ((kCO2_UMBCorHARTMAN .EQ. +1) .AND. (kAltComprDirs .EQ. -1)) THEN
        iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
      ELSEIF ((kCO2_UMBCorHARTMAN .EQ. -1) .OR. (kAltComprDirs .EQ. +1)) THEN
        iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing, or other databases
      END IF

      !!!! iCO2Chi = 2   TESTING

 10   FORMAT(I2,I2,I2)
 999  FORMAT('CO2 chi fudge iDefault = ',I2,' iCO2Chi = ',I2,' kCO2_UMBCorHARTMAN = ',I2)
      IF ((iCO2Chi .NE. iDefault) .AND. (kAltComprDirs .EQ. -1) .AND. (kOuterLoop .EQ. 1)) THEN      
        write(kStdWarn,999) iDefault,iCO2Chi,kCO2_UMBCorHARTMAN
        write(kStdErr,999)  iDefault,iCO2Chi,kCO2_UMBCorHARTMAN	
c      ELSEIF ((iCO2Chi .NE. iDefault) .AND. (kAltComprDirs .EQ. +1) .AND. (kOuterLoop .EQ. 1)) THEN
c        write(kStdErr,*) ' CO2 chi fudge iDefault,iCO2Chi = ',iDefault,iCO2Chi,' but using other user suppl CO2 dir'
      END IF

      IF (iCO2Chi .EQ. 2) THEN
c this is old; prior to March 2004
c notice how we only need to fudge fix 2255,2280 and 2305,2405 chunks here!
        iChiChunks = 4
        iaChiChunks(1) = 2255
        iaChiChunks(2) = 2280
        iaChiChunks(3) = 2380
        iaChiChunks(4) = 2405
      ELSE IF (iCO2Chi .EQ. 3) THEN
c this is new; after March 2004
c notice we fix 15 um and 4 um here
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

      IF (iCO2Chi .GT. 0) THEN
        iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
c        print *,int(rFileStartFr),iCO2Chi,iDoFudge
        IF (iDoFudge .GT. 0) THEN
          write(kStdWarn,*) ' doing CO2 chi fudge for ',int(rFileStartFr),iCO2Chi,iDoFudge
          write(kStdErr,*) ' doing CO2 chi fudge for ',int(rFileStartFr),iCO2Chi,iDoFudge
          CALL co2_4um_fudge(daaAbsCoeff,rFileStartFr,
     $                         iCO2Chi,iaChiChunks,iChiChunks)
        END IF 
      END IF
       
      RETURN
      END

c************************************************************************
! this subroutine mutiplies the daaGasAbsCoeff by CO2/WV chi functions
! reference : Measurements and modeling of absorption by CO2+ H2O mixtures in
!    the spectral region beyond the CO2 ν3-band head : H. Tran, M. Turbet,
!    P. Chelin, X. Landsheere, Icarus 306 (2018) 116–121
!
! od ~ rho^2 xco2 xwv L CA   wheo rho = denisity in amagat, xco2/xwv are the mix ratios. L = path length of layer
! Ha Tran told me no experimentally measured T dependance, no d(CA)/dT = 0; however the density changes with
! T since rho ~ P/RT so d(rho^2)/dT ~ -2P/RT^3

      SUBROUTINE add_co2_wv_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness,
     &                                    daaDQ,daaDT,iDoDQ)
     
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

! input
      INTEGER :: iGasID,iDoDQ
      REAL :: raFreq(kMaxPts)
      REAL :: raPress(kProfLayer),raaPartPress(kProfLayer,kGasStore),raTemp(kProfLayer)
      REAL :: raThickness(kProfLayer)
! input/output
      DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

! local vars
      LOGICAL :: isfinite2
      INTEGER :: iCO2Chi,iDefault,iI,iFr,iNumPts
      INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge

      INTEGER :: iIOUN,iERR
      CHARACTER(120) :: FNAME
      REAL :: raFChi(kMaxPts),raChi(kMaxPts),raX(kMaxPts),rScale
      REAL :: raXCO2(kProfLayer),raXWV(kprofLayer),raRho(kProfLayer)

      FNAME = '//home/sergio/SPECTRA/CKDLINUX/ca_wv_co2_forkcarta_2000_2900.dat'
    
      iCO2Chi = 0  !!no continuum chi fixes applied .. with database being
      iCO2Chi = 2  !!default July 2018; only fixes 4um

      iDefault = 2
      iCO2Chi = 2    !!! DEFAULT
      iCO2Chi = iaaOverrideDefault(1,9)
      IF ((iCO2Chi /= 0) .AND. (iCO2Chi /= 2)) THEN
        write(kStdErr,*) 'invalid iCO2/WV Chi = ',iCO2Chi
        CALL DoStop
      END IF

      iCO2Chi = -1    !! testing

 10   FORMAT(I2,I2,I2)
 999  FORMAT('CO2 chi fudge iDefault = ',I2,' iCO2Chi = ',I2,' kCO2_UMBCorHARTMAN = ',I2)
 
      IF (iGasID .EQ. 2 .and. iCO2Chi > 0 .and. (raFreq(1) >= 2355.0 .and. raFreq(1) <= 2805.0))  THEN

        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
        IF (IERR /= 0) THEN
          WRITE(kStdErr,*) 'In subroutine multiply_co2_wv_continuum'
          WRITE(kStdErr,1010) IERR, FNAME
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
          CALL DoSTOP
        ENDIF

        kTempUnitOpen=1
        READ(iIOUN) iNumPts
        READ(iIOUN) (raFChi(iFr),iFr=1,iNumPts)
        READ(iIOUN) (raChi(iFr),iFr=1,iNumPts)      
        CLOSE(iIOUN)
        kTempUnitOpen=-1

        !! get the abs coeff in cm-1, normalized by (rho)^2 Xco2 Xwv
        CALL rspl(raFChi,raChi,iNumPts,raFreq,raX,kMaxPts)
      
!      print *,-9999,raFChi(1),raChi(1),raFreq(1),raX(1)
!      DO iI = 1,iNumPts
!        print *,iI,iNumPts,raFChi(iI),raChi(iI)
!      END DO
!      DO iI = 1,kMaxPts,100
!        print *,iI,raFreq(iI),raX(iI)
!      END DO
      
        DO iI = 1,kProfLayer
          raXWV(iI)  = raaPartPress(iI,1)/raPress(iI)      ! fraction of WV
          raXCO2(iI) = raaPartPress(iI,2)/raPress(iI)      ! fraction of CO2
	  !! now compute density in amagats
	  raRho(iI) = (raPress(iI) * kAtm2mb * 100)/(kMGC * raTemp(iI))*kAvog/1000.0   !! molecules/m3
	  raRho(iI) = raRho(iI)/1.0e6/2.6867805e19                                     !! amagat
	
          rScale = raRho(iI)*raRho(iI)*raXWV(iI)*raXCO2(iI)*raThickness(iI)*100.0      !! thickness m --> cm	
!          print *,iI,raXWV(iI),raXCO2(iI),raRho(iI),raThickness(iI),rScale,raX(1),daaAbsCoeff(1,iI)
        END DO

        DO iI = 1,kProfLayer
          rScale = raRho(iI)*raRho(iI)*raXWV(iI)*raXCO2(iI)*raThickness(iI)*100.0      !! thickness m --> cm
	  if (isfinite2(rScale) .EQV. .false. ) rScale = 0.0
          DO iFr = 1,kMaxPts
	    daaAbsCoeff(iFr,iI) = daaAbsCoeff(iFr,iI) + rScale*raX(iFr)
	  END DO
        END DO
      END IF
           
      RETURN
      end SUBROUTINE add_co2_wv_continuum

!************************************************************************

      LOGICAL FUNCTION isfinite2(a)
      REAL :: a
      isfinite2 = (a-a) == 0
      end FUNCTION isfinite2

!************************************************************************
c baed on LLStrow's analysis of the AERI data, this subroutine further scales 
c the CO2 absorption coefficients in the 2380-2405 regions
      SUBROUTINE AERI_LLS_fudge(daaAbsCoeff,rFileStartFr,iWhichGas)

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c rFileStartFr = which chunk
      INTEGER iWhichGas
      REAL rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iChi
      REAL raFreq(kMaxPts),raMystery(kMaxPts)

      iChi = -1
      IF (iWHichGas .EQ. 1) THEN
        write(kStdWarn,*) 'need H20 AERI chifunction for ',rFileStartFr
        iChi = +1
      ELSEIF (iWHichGas .EQ. 2) THEN
        write(kStdWarn,*) 'need CO2 AERI chifunction for ',rFileStartFr
        iChi = +1
      END IF

      IF (iChi .GT. 0) THEN
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
      END
     
c************************************************************************
c this subroutine scales the CO2 absorption coefficients in the 2355,2280 
c and 2380,2405 chunks
c look at chi.f for the fudge for 2380,2405 cm-1
      SUBROUTINE Co2_4um_fudge(daaAbsCoeff,rFileStartFr,
     $                       iCO2Chi,iaChiChunks,iChiChunks)

      include '../INCLUDE/kcarta.param'

c iCO2Chi is for the following :
c the blah.txt  files are the orig,      done in late june/early july 2002
c the blah_a.txt files are refinement 1, done in early aug 2002
c the blah_b.txt files are refinement 2, done in mid nov 2002       iCO2Chi = 2
c   Scott did further refinemnents in 2003 and Jan 2004; see file   iCO2Chi = 3
c   /home/sergio/SPECTRA/CKDLINUX/co2_tune.m

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c rFileStartFr = which chunk
      REAL rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      INTEGER iaChiChunks(kMaxGas),iChiChunks,iCO2Chi

      INTEGER iDoFudge,WhichGasPosn,iI,iJ,iK
      INTEGER iFr,iLay,iIOUN,iERR,iChi
      REAL raF(kMaxPts),raChi(kMaxPts),m,c,rF,rX,rChi,r1,r2,r3,rAmp,rAmp0
      CHARACTER*120 FNAME
      CHARACTER*3 ca3
      CHARACTER*4 ca4

      IF ((iCO2Chi .NE. 2) .AND. (iCO2Chi .NE. 3)) THEN
        write(kStdErr,*) 'Illegal type for co2 chi function',iCO2Chi
        CALL DoStop
      END IF

      iChi = -1
      iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
      IF (iDoFudge .GT. 0) THEN
        iChi = +1
        DO iI = 1,120
          fname(iI:iI) = ' '
        END DO
        FNAME = 'co2_4um_fudge_'
        iJ = 1
 11     CONTINUE
        IF ((fname(iJ:iJ) .NE. ' ') .AND. (iJ .LT. 120)) THEN
          iJ = iJ + 1
          GOTO 11
        END IF
        IF (rFileStartFr .LT. 1000) THEN
          write(ca3,30) nint(rFileStartFr)
          fname(iJ:iJ+2) = ca3(1:3)
          iJ = iJ+3
        ELSEIF (rFileStartFr .GE. 1000) THEN
          write(ca4,40) nint(rFileStartFr)
          fname(iJ:iJ+3) = ca4(1:4)
          iJ = iJ+4
        END IF
      END IF

      IF (iCO2Chi .EQ. 2) THEN
        fname(iJ:iJ+5) = '_b.txt'
      ELSEIF (iCO2Chi .EQ. 3) THEN
        fname(iJ:iJ+5) = '_c.txt'
      END IF

 30    FORMAT(I3)
 40    FORMAT(I4)

c      IF (iFileStartFR .EQ. 2255) THEN
c        write(kStdWarn,*) 'need CO2 chifunction for 2255 chunk ....'
c        FNAME = 'co2_4um_fudge_2255.txt'
c        FNAME = 'co2_4um_fudge_2255_a.txt'
c        FNAME = 'co2_4um_fudge_2255_b.txt'  !!!'a' and 'b are copies 
c        iChi = +1
c      ELSEIF (iFileStartFR .EQ. 2280) THEN
c        write(kStdWarn,*) 'need CO2 chifunction for 2280 chunk ....'
c        FNAME = 'co2_4um_fudge_2280.txt'
c        FNAME = 'co2_4um_fudge_2280_a.txt'
c        FNAME = 'co2_4um_fudge_2280_b.txt'  !!!'a' and 'b are copies 
c        iChi = +1
c      ELSEIF (iFileStartFR .EQ. 2380) THEN
c        write(kStdWarn,*) 'need CO2 chifunction for 2380 chunk ....'
c        FNAME = 'co2_4um_fudge_2380.txt'
c        FNAME = 'co2_4um_fudge_2380_b.txt'
c        iChi = +1
c      ELSEIF (iFileStartFR .EQ. 2405) THEN
c        write(kStdWarn,*) 'need CO2 chifunction for 2405 chunk ....'
c        FNAME = 'co2_4um_fudge_2405.txt'
c        FNAME = 'co2_4um_fudge_2405_b.txt'
c        iChi = +1
c      END IF
 
      CALL FindChiFileName(fname)

      IF (iChi .GT. 0) THEN
        write(kStdWarn,*) '   Reading in CO2 chifile ',fname
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine co2_4um_fudge'
          WRITE(kStdErr,1010) IERR, FNAME
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
          CALL DoSTOP
        ENDIF
        kTempUnitOpen=1
        READ(iIOUN,*) (raF(iFr),raChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen=-1

c to print out the chifcn
c          DO iFr=1,kMaxPts
c           print *,raF(iFr),raChi(iFr)
c           end do

        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*raChi(iFr)
          END DO
        END DO
      ELSE
        write(kStdWarn,*) 'do not need CO2 chifunction for chunk ',rFileStartFr
      END IF        

      RETURN
      END

c************************************************************************
c this subroutine reads in a generic Scott Hannon tuning file, and applies
c the tuning to INDIVIDUAL gas optical depths (ie daaAbsCoeff)
c basically copied from SUBROUTINE water_sarta_fudge
      SUBROUTINE generic_sarta_tunmult(iGasID,raFreq,raaAbsCoeff,iSARTAChi)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      
c input
c   iSARTAChi = -1 for no chi, +1 for AIRS, +2 for IASI, +3 for CrIS
c   iGasID    = gasID
c   raFreq    = wavenumber array
c  raaAbsCoeff = gas OD

      INTEGER iGasID,iSARTAChi
      REAL rFileStartFr,raFreq(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iIOUN,iERR,iMin,iMax,iNpts
      REAL raF(kMaxPts),raChi(kMaxPts),raChi2(kMaxPts)
      REAL rF,fixed,water,watercon,o3,co,ch4,nte
      CHARACTER*120 Fname,caStr
      REAL raFX(kMaxPts),raChiX(kMaxPts)

      REAL raBeginBand(100),raEndBand(100)
      INTEGER iBand,iX,iaOK(kMaxPts)

c see /asl/packages/sartaV108_PGEv6/Src_AIRS_PGEv6/tunmlt_df.f
C    The tuning multiplier file must consist of MXCHAN lines of data
C    sorted (in ascending order) by channel ID and containing the
C    following 9 columns:
C       1    2   3    4     5    6   7   8    9
C       ID RJUNK XF XH2OL XH2OC XO3 XCO XCH4 XNTE

c      write(kStdWarn,*) 'inside generic_sarta_tunmult for GasID = ',iGasID,' raFreq(1) = ',raFreq(1),' iSARTAChi = ',iSARTAChi
c      write(kStdErr,*)  'inside generic_sarta_tunmult for GasID = ',iGasID,' raFreq(1) = ',raFreq(1),' iSARTAChi = ',iSARTAChi      
      
      Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.txt'
      Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_ones.txt'
      IF (iSARTAChi .EQ. +1) THEN
        !! /asl/packages/sartaV108_PGEv6/Src_AIRS_PGEv6/tunmlt_df.f
        Fname = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_PGEv6.txt'
        Fname = '/asl/data/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt'
      ELSEIF (iSARTAChi .EQ. +2) THEN
        !! /asl/packages/sartaV108/Src_rtpV201_pclsam_slabcloud_hg3/incFTC_iasi_may09_wcon_nte_swch4.f
        Fname = '/asl/data/sarta_database/Data_IASI_may09/Coef/tunmlt_wcon_nte.txt'
      ELSE
        write(kStdErr,*) 'only doing SARTA AIRS/IASI tunings'
        CALL DoStop
      END IF

      iIOUN = kTempUnit

      rFileStartFr = raFreq(1)

 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)

      IF ((rFileStartFr .GE. 605.0) .AND. (rFileStartFr .LE. 2805.0)) THEN
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine generic_sarta_tunmult'
          WRITE(kStdErr,1010) IERR, FNAME
          CALL DoSTOP
        ENDIF

        iNpts = 0
        kTempUnitOpen=1
 10     CONTINUE
        READ(iIOUN,80,ERR=20,IOSTAT=IERR) caStr
        IF ((caStr(1:1) .EQ. '%') .AND. (IERR .EQ. 0)) THEN
          GOTO 10
        ELSEIF (IERR .EQ. 0) THEN
          iNpts = iNpts + 1
          read(caStr,*) iFr,rF,fixed,water,watercon,o3,co,ch4,nte
          raF(iNpts)   = rF
          raChi(iNpts) = 1.0
          IF ((iGasID .EQ. 1) .OR. (iGasID .EQ. 103)) THEN
            raChi(iNpts) = water
c          ELSEIF (iGasID .EQ. 2) THEN
c            raChi(iNpts) = co2
          ELSEIF (iGasID .EQ. 3) THEN
            raChi(iNpts) = o3
          ELSEIF (iGasID .EQ. 5) THEN
            raChi(iNpts) = co
          ELSEIF (iGasID .EQ. 6) THEN
            raChi(iNpts) = ch4
          ELSEIF (iGasID .LE. 100) THEN
            raChi(iNpts) = fixed
          END IF
	ELSEIF (iERR .NE. 0) THEN
	  GOTO 20	  
        END IF
        GOTO 10

 20     CONTINUE
        CLOSE(iIOUN)
        kTempUnitOpen = -1

        IF (iSARTAChi .EQ. +1) THEN
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

cx        call r_sort_spl(raF,raChi,2378,raFreq,raChi2,kMaxPts)
        call r_sort_spl(raF,raChi,iNpts,raFreq,raChi2,kMaxPts)

c        !! now check simple boundaries, between [raF(1) and raF(end)]
c        IF (raFreq(1) .LT. raF(1)) THEN
c          iMin = 1
c          IF (raFreq(kMaxPts) .LT. raF(1)) THEN
c            iMax = kMaxPts
c          ELSE
c            iMax = int((raF(1)-raFreq(1))/0.0025)
c          END IF
c          DO iFr=iMin,iMax
c            raChi2(iFr) = 1.0
c          END DO
c        END IF
c        IF (raFreq(kMaxPts) .GT. raF(2378)) THEN
c          iMax = kMaxPts
c          IF (raFreq(1) .GT. raF(2378)) THEN
c            iMin = 1
c          ELSE
c            iMin = int((raF(2378)-raFreq(1))/0.0025)
c          END IF
c          DO iFr=iMin,iMax
c            raChi2(iFr) = 1.0
c          END DO
c        END IF

        !! now check more complicated boundaries
        CALL DoSort2Real(raF,raChi,iNpts,1)
        iBand = 0
        iX = 1
 30     CONTINUE
        IF (abs(raChi(iX)-1.0) .GE. 1.0e-5) THEN
          !! found  the start of a "band" at which to multiply ODs
          iBand = iBand + 1
          raBeginBand(iBand) = raF(iX)
          iX = iX + 1
  40      CONTINUE
          !! now look for the band end
          IF ((abs(raChi(iX)-1.0) .GE. 1.0e-5) .AND. (iX .LT. iNpts)) THEN
            iX = iX + 1
            GOTO 40
          ELSE
            raEndBand(iBand) = raF(iX)
          END IF
        END IF
        IF (iX .LT. iNpts) THEN
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
            IF ((raFreq(iFr) .GE. raBeginBand(iX)) .AND. (raFreq(iFr) .LE. raEndBand(iX))) THEN
              iaOK(iFr) = +1
            END IF
          END DO
        END DO
        DO iFr = 1,kMaxPts
          IF (iaOK(iFr) .EQ. -1) raChi2(iFr) = 1.0
        END DO

c        IF (iGasID .LE. 6) THEN
c          print *,iGasID,iBand
c          print *,(raBeginBand(iX),iX=1,iBand)
c          print *,(raEndBand(iX),iX=1,iBand)
c        END IF

        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            raaAbsCoeff(iFr,iLay) = raaAbsCoeff(iFr,iLay)*max(0.0,raChi2(iFr))
          END DO
        END DO

      END IF

 80   FORMAT(A120)

      RETURN
      END
     
c************************************************************************
c this subroutine scales water absorption coefficients in the >= 2380 chunks
      SUBROUTINE water_4um_fudge(daaAbsCoeff,rFileStartFr)

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c rFileStartFr = which chunk
      REAL rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iIOUN,iERR,iChi,i1,i2,iFileEnd
      REAL raF(kMaxPts),raChi(kMaxPts)
      DOUBLE PRECISION dX,daX(kMaxPts),dSlope,xA,x

      iChi = -1
      IF (abs(rFileStartFR-2830.0) .LE. 0.0001) THEN
        write(kStdWarn,*) 'need water chifunction for ',nint(rFileStartFR),' 
     $ chunk ....'
        iChi = +2
      ELSEIF (rFileStartFR .GE. 2405.0) THEN
        write(kStdWarn,*) 'need water chifunction for ',nint(rFileStartFR),' 
     $ chunk ....'
        iChi = +1
      END IF

      IF (iChi .EQ. 1) THEN
        !!! for all freqs, for all layers, multiply lines by 0.93
        dX = 0.93d0
        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*dX
          END DO
        END DO
      ELSEIF (iChi .EQ. 2) THEN
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
      END
     
c************************************************************************
c this subroutine scales water absorption coefficients in the >= 2380 chunks
      SUBROUTINE water_sarta_fudge(daaAbsCoeff,rFileStartFr)

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c rFileStartFr = which chunk
      REAL rFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iIOUN,iERR,iChi
      REAL raFreq(kMaxPts)
      REAL raF(kMaxPts),raChi(kMaxPts),raChi2(kMaxPts)
      REAL rF,fixed,water,watercon,o3,co,ch4,co2
      CHARACTER*120 Fname,caStr

      iChi = -1
      IF (rFileStartFR .GE. 1380.0 .AND. rFileStartFR .LE. 1680.0) THEN
        write(kStdWarn,*) 'need water_sarta_fudge for ',nint(rFileStartFR),' 
     $ chunk ....'
        write(kStdWarn,*)' h20_sarta_fudge for ',nint(rFileStartFR),' chunk ..'
        iChi = +1
      END IF

      IF (iChi .EQ. 1) THEN
        Fname = '/home/sergio/SPECTRA/CKDLINUX/tunmlt_jan04deliv.txt'
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine water_sarta_fudge'
          WRITE(kStdErr,1010) IERR, FNAME
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
          CALL DoSTOP
        ENDIF
        kTempUnitOpen=1
 10     CONTINUE
        READ(iIOUN,80) caStr
        IF (caStr(1:1) .EQ. '%') THEN
          GOTO 10
        ELSE 
          read(caStr,*) iFr,rF,fixed,water,watercon,o3,co,ch4,co2
          raF(iFr)   = rF
          raChi(iFr) = water
          IF (iFr .LT. 2378) THEN
            GOTO 10
          ELSE
            GOTO 20
          END IF
        END IF
 20     CONTINUE
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

 80   FORMAT(A120)

      RETURN
      END
     
c************************************************************************
c this subroutine scales the absorption coefficients by looking at the 
c amounts in the gas profiles, thereby computing the optical depth
c compute optical depth = gas amount * abs coeff
      SUBROUTINE AmtScale(daaAbsCoeff,raPAmt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c raP/RAmt    = actual/reference gas amount profiles
      REAL raPAmt(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay
      REAL rScale
      DOUBLE PRECISION dZero

      dZero=DBLE(0.0)

      DO iLay=1,kProfLayer
        rScale=raPAmt(iLay)
        DO iFr=1,kMaxPts
          daaAbsCoeff(iFr,iLay)=max(daaAbsCoeff(iFr,iLay),dZero)*rScale
        END DO
      END DO

      RETURN
      END
     
c************************************************************************
c this subroutine finishes the calculation of d/dq(abscoeff)
      SUBROUTINE FinalAmtDeriv(daaDQ,iType)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaDT   = d/dT
c iType   = compression type
      INTEGER iType
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

      INTEGER iLay,iFr

      IF (iType .EQ. 2) THEN
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
      END 

c************************************************************************
c this subroutine finishes the computation of d/dT (absCoeff) 
      SUBROUTINE FinalTempDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iKtype      = compression type (power = 1 or 4)
c daaAbsCoeff = uncompressed gas abs coeffs
c daaDT       = d/dT
c ra(P/R)Amt  = reference/actual gas amounts
      INTEGER iKtype
      REAL raPAmt(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDA(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iLay
      REAL rScale

      IF (iKtype .EQ. 1) THEN
        DO iLay=1,kProfLayer
          rScale = raPAmt(iLay)
          DO iFr=1,kMaxPts
            daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale
          END DO
        END DO
      ELSE 
c remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
c what comes from the spline interpolation routines is d/dT (K(v,T))^(1/4)
c or in other words, we get = (1/4) K(v,T)^(-3/4) dK/dT = daaDA
c we need d/dT(optical depth) = d/dT q(actual) K(v,T) = q(actual) d/dT K(v,T)
c so we get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3 
        DO iLay=1,kProfLayer
          rScale = raPAmt(iLay)
          DO iFr=1,kMaxPts
            daaDA(iFr,iLay)=daaDA(iFr,iLay)*rscale*
     $                      4.0*(daaAbsCoeff(iFr,iLay)**3)
          END DO
        END DO
      END IF        

      RETURN
      END

c************************************************************************
c this subroutine finishes the computation of d/dq (absCoeff) for water
      SUBROUTINE FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iKtype      = compression type (power = 1 or 4)
c daaAbsCoeff = uncompressed gas abs coeffs
c daaDT       = d/dT
c ra(P/R)Amt  = reference/actual gas amounts
      INTEGER iKtype
      REAL raPAmt(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDA(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iLay
      REAL rScale

c similar to FinalTempDeriv above
c remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
c what comes from the spline interpolation routines is d/dq (K(v,T))^(1/4)
c or Zget = (1/4) K(v,T)^(-3/4) dK/dq = daaDA
c we need d/dq(optical depth) = d/dq q(actual) K(v,T) = q(actual) d/dq K(v,T)
c                                                        + K(v,T)
c get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3 +  daaAbsCoeff^4
      DO iLay=1,kProfLayer
        rScale=raPAmt(iLay)
        DO iFr=1,kMaxPts
          daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale*
     $                      4.0*(daaAbsCoeff(iFr,iLay)**3)
     $        + (daaAbsCoeff(iFr,iLay)**4)
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine raises the compressed matrix elements to the 4th power
      SUBROUTINE RaisePower(daaAbsCoeff)

      IMPLICIT NONE
      
      include '../INCLUDE/kcarta.param'

c daaGasAbsCoeff = uncompressed, scaled gas abs coeffs
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iLay,iFr

      DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
          daaAbsCoeff(iFr,iLay)=(daaAbsCoeff(iFr,iLay))**4
        END DO
      END DO

      RETURN
      END

c************************************************************************
CCCCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
C      Reads a compressed K and U data file for water
       SUBROUTINE RDCOMPWATER(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, 
     $               NLAY, KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT,
     $               KX1, KX2, KX3, KX4, KX5, UX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

C      Calling parameters
c fnam     = name of relevant file that contains compressed database
c iIOUN    = file IOUNIT number
c IDGAS    = GAS ID
c sfreq    = start freq
c fstep    = frequency point spacing
c npts     = number of freq points === kMaxPts
c nlay     = number of layers == kProfLayer
c ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
c nk       = number of singular vectors (<= kMaxK)
c kt       = number of temperature profiles (=11)
c kn       = number of layers == kMaxLayers
c um       = number of freq points == kMaxPts (in UX)
c un       = number of singular vectors (in UX)
c KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
c UX       = uncompression matrix
c ITO,ITSORT = base points to do temperature interpolation
       CHARACTER*120 FNAM
       INTEGER iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0,
     $    ITSORT(kMaxTemp)
       DOUBLE PRECISION SFREQ,FSTEP,TOFF(kMaxTemp),
     $    KX1(kMaxK,kMaxTemp,kMaxLayer),
     $    KX2(kMaxK,kMaxTemp,kMaxLayer),
     $    KX3(kMaxK,kMaxTemp,kMaxLayer),
     $    KX4(kMaxK,kMaxTemp,kMaxLayer),
     $    KX5(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

C     Local variables
      INTEGER IERR,I,J,K,RESORT,IHOLD
      REAL rTemp

 1010   FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMPWATER'
        WRITE(kStdErr,1010) IERR, FNAM
        CALL DoSTOP
      ENDIF
      kCompUnitOpen=1

C     Read in the header
      READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN,
     $    UM, UN

 1110 FORMAT('Error! Compressed data array dimension exceeds ',
     $       'max size')
 1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
C     Make sure the array sizes are <= to the declared sizes
      IF (NK .GT. kMaxK) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1120) NK, kMaxK
        CALL DoSTOP
      ENDIF

 1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
      IF (KT .GT. kMaxTemp) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1130) KT, kMaxTemp
        CALL DoSTOP
      ENDIF

 1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
      IF (KN .GT. kMaxLayer) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1140) KN, kMaxLayer
        CALL DoSTOP
      ENDIF

 1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
      IF (UM .GT. kMaxPts) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1150) UM, kMaxPts
        CALL DoSTOP
      ENDIF

 1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
      IF (UN .GT. kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
        CALL DoSTOP
      ENDIF

C     Read in the temperature offsets
      READ(iIOUN) (TOFF(I),I=1,KT)

c check to make sure the offsets differ by kTempOffSet_database
      DO I=1,KT-1
        rTemp = abs(TOFF(I)-TOFF(I+1))
        rTemp = abs(rTemp - kTempOffSet_database)
        IF (rTemp .GT. 0.001) THEN
          write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
          write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
          write(kStdErr,*) (TOFF(J),J=1,KT)
          write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
          write(kStdErr,*) I,kTempOffSet_database         
          CALL DoStop
        END IF
      ENDDO

C      Find which of the offsets is 0 (error if none).
      IT0=0
      DO I=1,KT
        IF (TOFF(I) .EQ. 0.0) IT0=I
        ITSORT(I)=I
       ENDDO
       IF (IT0 .EQ. 0) THEN
         WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
 1180    FORMAT('ERROR! One of the temperature offsets must be 0',/,
     $       'offsets =',20(' ',F5.1))
         CALL DoSTOP
       ENDIF
C
C     Sort the indices of the temperature offsets in ascending order
      RESORT=1
 10   IF (RESORT .EQ. 1) THEN
        RESORT=0
        DO I=1,KT-1
          IF (TOFF( ITSORT(I) ) .GT. TOFF( ITSORT(I+1) )) THEN
            IHOLD=ITSORT(I)
            ITSORT(I)=ITSORT(I+1)
            ITSORT(I+1)=IHOLD
            RESORT=1
          ENDIF
        ENDDO
      GOTO 10
      ENDIF

C      Read in the five K matrices
c for the old mat2for files READ(iIOUN) ((KX1(I,J,K),J=1,KT),K=1,KN)
c for the WATER mat2for files, READ(iIOUN) ((KX1(I,J,K),K=1,KN),J=1,KT)
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

C     Read in the U matrix
      DO I=1,NK
         READ(iIOUN) (UX(J,I),J=1,NPTS)
      ENDDO

      CLOSE(iIOUN)
      kCompUnitOpen=-1

      RETURN
      END

c************************************************************************
CCCCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
C      Reads a compressed K and U data file for gases other than water
       SUBROUTINE RDCOMP(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, NLAY,
     $                   KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT,
     $                   KX, UX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c fnam     = name of relevant file that contains compressed database
c iiIOUN    = file IOUNIT number
c IDGAS    = GAS ID
c sfreq    = start freq
c fstep    = frequency point spacing
c npts     = number of freq points === kMaxPts
c nlay     = number of layers == kMaxLayer
c ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
c nk       = number of singular vectors (<= kMaxK)
c kt       = number of temperature profiles (=11)
c kn       = number of layers == kMaxLayers
c um       = number of freq points == kMaxPts (in UX)
c un       = number of singular vectors (in UX)
ccc un should be equal to nk
c KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
c UX       = uncompression matrix
c ITO,ITSORT = base points to do temperature interpolation
C      Calling parameters
       CHARACTER*120 FNAM
       INTEGER iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0,
     $    ITSORT(kMaxTemp)
       DOUBLE PRECISION SFREQ,FSTEP,TOFF(kMaxTemp),
     $    KX(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)
C
C      Local variables
       INTEGER IERR,I,J,K,RESORT,IHOLD,iDebugMatlab
       REAL rTemp
C-----------------------------------------------------------------------
C
 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMP'
        WRITE(KSTDERR,1010) IERR, FNAM
        CALL DoSTOP
      ENDIF
      kCompUnitOpen=1
C
C     Read in the header
      READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN,
     $    UM, UN
C
 1110 FORMAT('Error! Compressed data array dimension exceeds ',
     $       'max size')
 1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
C     Make sure the array sizes are <= to the declared sizes
      IF (NK .GT. kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1120) NK, kMaxK
        CALL DoSTOP
      ENDIF
C
 1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
      IF (KT .GT. kMaxTemp) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1130) KT, kMaxTemp
        CALL DoSTOP
      ENDIF
C
 1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
      IF (KN .GT. kMaxLayer) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1140) KN, kMaxLayer
        CALL DoSTOP
      ENDIF
C
 1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
      IF (UM .GT. kMaxPts) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1150) UM, kMaxPts
        CALL DoSTOP
      ENDIF
C
 1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
      IF (UN .GT. kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
        CALL DoSTOP
      ENDIF
C
C      Read in the temperature offsets
      READ(iIOUN) (TOFF(I),I=1,KT)

c check to make sure the offsets differ by kTempOffSet_database
      DO I=1,KT-1
        rTemp = abs(TOFF(I)-TOFF(I+1))
        rTemp = abs(rTemp - kTempOffSet_database)
        IF (rTemp .GT. 0.001) THEN
          write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
          write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
          write(kStdErr,*) (TOFF(J),J=1,KT)
          write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
          write(kStdErr,*) I,kTempOffSet_database         
          CALL DoStop
        END IF
      ENDDO
C
C      Find which of the offsets is 0 (error if none).
 1180  FORMAT('ERROR! One of the temperature offsets must be 0',/,
     $       'offsets =',20(' ',F5.1))
      IT0=0
      DO I=1,KT
        IF (TOFF(I) .EQ. 0.0) IT0=I
        ITSORT(I)=I
      ENDDO
      IF (IT0 .EQ. 0) THEN
        WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
          CALL DoSTOP
       ENDIF
C
C     Sort the indices of the temperature offsets in ascending order
      RESORT=1
 10   IF (RESORT .EQ. 1) THEN
        RESORT=0
        DO I=1,KT-1
          IF (TOFF( ITSORT(I) ) .GT. TOFF( ITSORT(I+1) )) THEN
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
      IF (iDebugMatlab .gt. 0) THEN
        print *,idgas,nk,kn,kt,ktype
        !! kx is a matrix(nk,kt,kn) = matrix(nk,11,100) 
        !!   where nk = number of basis vectors
        !! this is matrix "kcomp" in the cg5v4250.mat files kcomp(nk,100,11)
        !! where the code is cg"IDGAS"v"FREQ".mat -- see abscmp/ReadmeVIS
       END IF
       DO I=1,NK
         READ(iIOUN) ((KX(I,J,K),K=1,KN),J=1,KT)
       ENDDO

C     Read in the U matrix
      !! this is matrix "B" in the cg5v4250.mat files  B(10000,nk)
      DO I=1,NK
        READ(iIOUN) (UX(J,I),J=1,NPTS)
      ENDDO
C
      CLOSE(iIOUN)
      kCompUnitOpen=-1

c      iDebugMatlab = 1
      IF (iDebugMatlab .gt. 0) THEN
        IF (idgas .EQ. 2) THEN
          print *,(UX(J,1),J=1,5)     !! should equal B(1:5,1)
          print *,(KX(1,J,6),J=1,6)   !! should equal kcomp(1,6,1:6)
        END IF
      END IF

      RETURN
      END

c************************************************************************
