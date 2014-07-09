c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c********* this file has the k-compressed routines **********************
c** which include reading in the data, doing the spline interpolations **
c********** and finding the absorption matrix for the relevant gas ******
c************************************************************************

c************************************************************************
c*********** NOTE THESE VARIABLES ARE DOUBLE PRECISION ******************
c************************************************************************
c
c  if kJacobian == 1 then call the spline-type routines
c  if kJacobian == -1 then call OLD routines
c where JAC( indicates routines for calculating jacobians, and OLD( indicates 
c original routines w/o jacobians
c
c Also, if iGasID > kMaxDQ then we do not need to calculate d/dq, only d/dT
c contribution. Since kMaxDQ >= 1, water d/dq,d/dT always calculated
c************************************************************************
c this is the MAIN routine
c this subroutine gets the contribution of the i-th gas to the
c absorbtion coefficients. Using the reference profile, it scales the
c contribution accordingly
c also, the gases are weighted by raaMix (from MIXFIL)
c iGasID is the gas ID, while iFileStartFr identifies the wavenumber block
      SUBROUTINE GasContribution(iCount,iGasID,iRefLayer,
     $  raRAmt,raRTemp,raRPress,raRPartPress,iL,iU,pProf,iProfileLayers,
     $  raAmt,raTemp,raPress,raPartPress,iaCont,caXsecF,
     $  raVTemp,iVTSet,iFileStartFr,iTag,raFreq,iErr,iDoDQ,
     $  daaDQ,daaDT,daaTemp)

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
c iFileStartFr   = which k-comp file chunk to uncompress
c iTag      = which k-comp file chunk to uncompress
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ     = analystic Jacobian wrt water amount
c daaDT     = analystic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raAmt(kProfLayer),raTemp(kProfLayer),pProf(kProfLayer)
      REAL raPress(kProfLayer),raPartPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iVTSet,iTag,iProfileLayers
      REAL raVTemp(kProfLayer),raFreq(kMaxPts)
      INTEGER iCount,iGasID,iL,iU,iErr,iRefLayer,iFileStartFr,iDoDQ
      CHARACTER*80 caXsecF
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer),kFrStep
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr,iLay

      iErr = -1
  
      IF ((1 .LE. iGasID) .AND. (iGasID .LE. kGasComp)) THEN
        kFrStep=kaFrStep(iTag)
        CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,iFileStartFr,iTag,
     $    raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp)
        END IF

      IF ((kGasXsecLo .LE. iGasID) .AND. (iGasID .LE. kGasXsecHi)) THEN
        kFrStep=kaFrStep(iTag)
        IF (KXsecFormat .GT. 0) THEN
          write (kStdWarn,*) 'using kCompressed Database format'
          CALL compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,iFileStartFr,iTag,
     $    raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp)
        ELSE
          write (kStdWarn,*) 'using binary file format'
          CALL CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep,daaTemp,
     $         raVTemp,iVTSet,raFreq,iErr,caXsecF,
     $         raAmt,raTemp,raPress,raPartPress,iaCont,
     $         daaDQ,daaDT,iDoDQ)
          END IF
        END IF

      IF ((kNewGasLo .LE. iGasID) .AND. (iGasID .LE. kNewGasHi)) THEN
        IF (kCKD .GE. 0) THEN           !add on continuum
         write(kStdWarn,*)'adding on continuum'
          kFrStep=kaFrStep(iTag)
          CALL newgases(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,daaTemp,raVTemp,iVTSet,iFileStartFr,iTag,
     $    raFreq,iErr,raAmt,raTemp,raPress,raPartPress,iaCont,
     $    daaDQ,daaDT,iDoDQ)
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

c this subroutine basically computes water continuum
      SUBROUTINE newgases(iCount,iGasID,iRefLayer,raRAmt,raRTemp,
     $ raRPress,raRPartPress,iL,iU,daaTemp,raVTemp,iVTSet,iFileStartFr,
     $ iTag,raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont,
     $ daaDQ,daaDT,iDoDQ)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iCount    = which of the iNumGases is being processed
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c iaCont    = whether or not to do continuum calculation .. iaCont(iCount)
c caXecF    = file name of cross section data
c daaTemp   = matrix containing the uncompressed k-spectra
c raVtemp   = vertical temperature profile for the Mixed paths
c iVTSet    = has the vertical temp been set, to check current temp profile
c iFileStartFr   = which k-comp file chunk to uncompress
c iTag      = which k-comp file chunk to uncompress
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ      = analytic Jacobian wrt water amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iVTSet,iDoDQ,iTag
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      REAL raVTemp(kProfLayer),raFreq(KmaxPts)
      INTEGER iGasID,iL,iU,iErr,iCount,iRefLayer,iFileStartFr
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

      INTEGER iAns
      REAL kFrStep

      iErr=0

      IF (iGasID .EQ. kNewGasLo) THEN
c check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep=kaFrStep(iTag)
        IF ((iAns .GT. 0) .AND. (kCKD .GE. 0)) THEN
          CALL AddContinuum(iGasID,iTag,iRefLayer,raFreq,raTAmt,raTTemp,
     $                        kFrStep,raTPress,raTPart,iL,iU,daaTemp,
     $                        daaDQ,daaDT,iDoDQ)
          write(kStdWarn,*) 'added self water continuum'
          END IF 
        END IF

      IF (iGasID .EQ. kNewGasHi) THEN
c check if user wants to include continuum .. if so, compute and include
        iAns=iaCont(1)        !check to see water continuum ON or OFF
        kFrStep=kaFrStep(iTag)
        IF ((iAns .GT. 0) .AND. (kCKD .GE. 0)) THEN
          CALL AddContinuum(iGasID,iTag,iRefLayer,raFreq,raTAmt,raTTemp,
     $                        kFrStep,raTPress,raTPart,iL,iU,daaTemp,
     $                        daaDQ,daaDT,iDoDQ)
          write(kStdWarn,*) 'added foreign water continuum'
          END IF 
        END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the contribution of Gases 1-kGasComp (if present)
      SUBROUTINE compressed(iCount,iGasID,iRefLayer,raRAmt,raRTemp,raRPress,
     $    raRPartPress,iL,iU,raVTemp,iVTSet,iFileStartFr,iTag,
     $    raFreq,iErr,raTAmt,raTTemp,raTPress,raTPart,iaCont,
     $    pProf,iProfileLayers,iDoDQ,
     $    daaDQ,daaDT,daaTemp)

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
c iFileStartFr   = which k-comp file chunk to uncompress
c iTag      = which k-comp file chunk to uncompress
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ      = analytic Jacobian wrt water amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer),pProf(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer),rCheckTemp
      INTEGER iaCont(kMaxGas),iVTSet,iDoDQ,iTag
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      REAL raVTemp(kProfLayer),raFreq(KmaxPts)
      INTEGER iGasID,iL,iU,iErr,iCount,iRefLayer,iFileStartFr,iProfileLayers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

      INTEGER iErr1,iCnt

      iErr1=0    !this sees if data successfully uncompressed

      IF (iErr1 .LE. 0) THEN
c set the vertical temperature profile if iGasID < 51 (first profile read)
        IF ((iGasID .LE. kGasComp).AND.(iVTSet .LT. 0)) THEN 
          write(kStdWarn,* )'Setting the vertical tempr profile ...'
          DO iCnt=1,kProfLayer
            raVTemp(iCnt)=raTTemp(iCnt)
            END DO
          iVTset=1
          END IF
c if previous profiles have been read in, check to make sure the 
c temperature profiles are the same!!!!
        IF ((iGasID .LE. kGasComp).AND.(iVTSet .GT. 0)) THEN 
          write(kStdWarn,*) 'Checking the vertical tempr profile ...'
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

c then read in the compressed data ... if water do the uncompression 
c differently than if the gas is one of the others
      iErr1=0
      IF (iErr .LE. 0) THEN
        IF (iGasID .EQ. 1) THEN
          CALL water(iGasID,iFileStartFr,iTag,iRefLayer,iL,iU,
     $      raTAmt,raRAmt,raTPart,raRPartPress,raTTemp,raRTemp,
     $      iErr1,iDoDQ,pProf,iProfileLayers,daaDQ,daaDT,daaTemp)
        ELSE 
          CALL othergases(iGasID,iFileStartFr,iTag,iRefLayer,iL,iU,
     $       raTAmt,raRAmt,raTTemp,raRTemp,iErr1,iDoDQ,pProf,iProfileLayers,
     $       daaDQ,daaDT,daaTemp)
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


 1000 FORMAT('Successfully read in profile for GAS ID',I2)
 1005 FORMAT('Successfully read in k-comp data for GAS ID ',I2)

 1001 FORMAT('Error ... Could not read in profile for GAS ID',I2)
 1006 FORMAT('Error ... Could not read k-comp data for GAS ID ',I2)


      RETURN
      END

c************************************************************************
c this function checks to see if there is NEW data for the  current gas
c if so, it returns a positive integer that tells the code which spectra
c dataset to use (from *SPECTRA) , else it returns -1
      INTEGER FUNCTION OutsideSpectra(iGasID,iNumNewGases,iaNewGasID) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iGasID       tells current gasID
c iNumNewGases tells how many new gases to use
c iaNewGasID   tells the gasIDs of the new gases
      INTEGER iGasID, iNumNewGases, iaNewGasID(kGasStore)

      INTEGER iI,iJ

      iI = -1

      IF (iNumNewGases .GT. 0) THEN
        iJ=1
c search to see if there is new data!     
 10     CONTINUE
        IF (iaNewGasID(iJ) .EQ. iGasID) THEN
          iI=iJ
        ELSEIF (iJ  .LT. iNumNewGases) THEN
          iJ = iJ + 1
          GOTO 10
          END IF        
        END IF

      OutsideSpectra = iI

      RETURN
      END

c************************************************************************

c this function checks to see if there is NEW data for the  current chunk
      INTEGER FUNCTION NewDataChunk(iNewIn,iaNewData,iaaNewChunks,
     $                             iFileStartFr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iNewIn        tells which NewSpectra set this gas corresponds to
c iaNewData     tells how many new data sets were read in for this set
c iaaNewChunks  tells which data chunks were read in for this set
c iFileStartFr is the integer ID of the relevant k-comp file 
      INTEGER iaNewData(kGasStore),iaaNewChunks(kGasStore,kNumKCompT)
      INTEGER iFileStartFr,iNewIn

      INTEGER iI,iJ

      iI = -1

      iJ=1
c search to see if there is new data!     
 10   CONTINUE
      IF (iaaNewChunks(iNewIn,iJ) .EQ. iFileStartFr) THEN
        iI=iJ
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
     $ iRefLayer,iL,iU,raTAmt,raTTemp,raTPress,raTPart,iaCont,iTag,raFreq,
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
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT) 
      REAL df,sf
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c daaDQ      = analytic Jacobian wrt water amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c iTag      = which k-comp file chunk to uncompress
c iCount    = which of the iNumGases is being processed
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPart(kProfLayer),raTPress(kProfLayer)
      INTEGER iaCont(kMaxGas),iDoDQ,iCount
      REAL raFreq(KmaxPts)
      INTEGER iL,iU,iErr,iRefLayer,iTag
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iIoun,I,J
      INTEGER IDGAS, NPTS, NLAY
      DOUBLE PRECISION SFREQ, FSTEP
      CHARACTER*80 FNAM

      FNAM=caaaNewChunks(iNewIn,iWhichChunk)

      write(kStdWarn,*) 'In ReadNewData, opening data file for GasID'
      write(kStdWarn,*) iGasID,FNAM

      iIOUN=kCompUnit
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine ReadNewData'
          WRITE(kStdErr,1010) IERR, FNAM
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
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
C     Read in the U matrix
      DO I=1,kProfLayer
        READ(iIOUN) (daaGasAbCoeff(J,I),J=1,kMaxPts)
        ENDDO

      CLOSE(iIOUN)
      kCompUnitOpen=-1

      write (kStdWarn,*) 'Read in NEW spectra for gasID ',IDGAS

      RETURN
      END

c************************************************************************
c this adds on the continuum
      SUBROUTINE AddContinuum(iGasID,iTag,iRefLayer,raFreq,raAmt,raTemp,
     $ kFrStep,raPress,raPartPress,iL,iU,daaTemp,daaDQ,daaDT,iDoDQ)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c daaTemp   = matrix containing the uncompressed k-spectra
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      INTEGER iGasID,iL,iU,iErr,iRefLayer,iDoDQ,iTag
      REAL raFreq(kMaxPts),kFrStep
      REAL raAmt(kProfLayer),raTemp(kProfLayer)
      REAL raPartPress(kProfLayer),raPress(kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iLay,iFr,iIOUN,iLoc
      REAL rFStep,rMin
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer),der1,der2
      REAL rEC,rECCount
      REAL raMult(kProfLayer)
      CHARACTER*80 caFName

c these are to read in the binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)
      INTEGER iCKD,iM,iN

      INTEGER idQT
      REAL a1,a2

c first check to see if exact calculations of d/dq,d/dT should be enabled
      IF (kJacobian .GE. 0) THEN 
        idQT=1
      ELSE
        idQT=-1
        END IF
 
c initialize raaCon to all zeros
      DO iFr=1,kMaxPts
        DO iLay=1,kProfLayer
          daaCon(iFr,iLay)=0.0
          END DO
        END DO

      iLay=kProfLayer
      iFr=kMaxPts
      rFStep=kFrStep

c calculate the continuum contribution
c figure out the filename
      CALL CKDFileName(caFName,iGasID,iTag)

      iIOUN = kTempUnit 
c open file and load in data
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=iErr)
      IF (iErr .NE. 0) THEN
        WRITE(kStdErr,1080) iErr, caFname
        WRITE(kStdWarn,1080) iErr, caFname
 1080   FORMAT('ERROR! number ',I5,' opening CKD binary file:',/,A82)
        CALL DoSTOP
        ENDIF

      kTempUnitOpen=1
      READ(iIOUN) d1,d2,df       !read start/stop freq, df
      IF (d1 .GT. raFreq(1)) THEN
        write(kStdErr,*) 'CKD file has freqs that start too high!'
        CALL DoStop
        END IF
      IF (d2 .LT. raFreq(kMaxPts)) THEN
        write(kStdErr,*) 'CKD file has freqs end start too low!'
        CALL DoStop
        END IF

c&&&&&&&&
      READ(iIOUN) iLoc,iCKD,iM,iN 
      !read line type, CKD vers, # of temp,freq pts
      IF (iLoc .LT. 0) THEN
        write(kStdErr,*) 'need correct continuum for local lineshape !'
        CALL DoStop
        END IF
      IF (iCKD .NE. kCKD) THEN
        write(kStdErr,*) 'CKD versions do not match!'
        CALL DoStop
        END IF
      IF (iM .NE. kTempCKD) THEN
        write(kStdErr,*) 'Need more Temp offsets in CKD binary file'
        CALL DoStop
        END IF
      IF (iN .NE. kFreqCKD) THEN
        write(kStdErr,*) 'Need more Freqs in CKD binary file'
        CALL DoStop
        END IF

      READ(iIOUN) (daTemprt(iLay),iLay=1,iM)   !read the temps

      DO iLay=1,iM
        !for current temp, read the continuum data for each freq
        READ(iIOUN) (daaCKD(iLay,iFr),iFr=1,iN) 
        END DO

      CLOSE(iIOUN)
      kTempUnitOpen=-1
c&&&&&&&&

c interpolate
c this is what was used in SRCv1.07, which is what Scott Hannon used
c for the Fast Models
c      CALL ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,
c     $                 kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID)

c this is the new and improved version, which should be quicker than Quad
      CALL ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,
     $                 kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID)

      rMin=1.0e10 
      rECCount=0.0 
      DO iLay=1,kProfLayer 
        DO iFr=1,kMaxPts 
          IF (daaCon(iFr,iLay) .lt. 0.0) THEN
            write(kStdWarn,*) 'continuum < 0 for iGasID,ilay,iFr',
     $                         iGasID,iLay,iFr,daaCon(iFr,iLay) 
            rECCount = rECCount + 1.0 
            IF (daaCon(iFr,iLay) .lt. rMin) THEN
              rMin = daaCon(iFr,iLay) 
              END IF 
            END IF
          END DO
        END DO

c add the continuum abs coeffs to the k-compressed coeffs raaCon = daaTemp
c keeping track of whether the added values are all greater than 0.0

      IF (iGasID .EQ. kNewGasLo) THEN      !CSelf
        DO iLay =1,kProfLayer
          raMult(iLay)=raAmt(iLay)*kAvog*raPartPress(iLay)
          END DO
      ELSEIF (iGasID .EQ. kNewGasHi) THEN      !CFor
        DO iLay =1,kProfLayer
          raMult(iLay)=raAmt(iLay)*kAvog*(raPress(iLay)-raPartPress(iLay))
          END DO
        END IF

      DO iLay =1,kProfLayer
        a1=raMult(iLay)*296.0/raTemp(iLay)
        a2=kPlanck2/2/raTemp(iLay)
        DO iFr=1,kMaxPts
          daaTemp(iFr,iLay)=daaCon(iFr,iLay)*a1*raFreq(iFr)*
     $                       tanh(a2*raFreq(iFr))
          END DO
        END DO

      IF (rECCount .gt. 0.5) THEN
        rEC=rEC/rECCount
        write(kStdWarn,*)'Error in CKD data!!! Some values negative!!!'
        write(kStdWarn,*)'and reset to 0.0'
        write(kStdWarn,*) rECCount,' values have avg value ',rEC
        write(kStdWarn,*) 'min negative value in 10000*100 = ',rMin
        IF (abs(rMin) .GT. 1.0e-7) THEN
          iErr=1
          write(kStdErr,*)'Error in CKD data!!! Some values negative!!!'
          CALL DoSTOP
          END IF
        END IF          

      IF (iDoDQ .GE. -1) THEN
        !!the gas amount jacobians
        IF (iGasID .EQ. 101) THEN   !!!there is a factor of 2
          DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
              daaDQ(iFr,iLay) = 2 * daaTemp(iFr,iLay)/raAmt(iLay) 
              END DO
            END DO
        ELSEIF (iGasID .EQ. 102) THEN
          DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
              daaDQ(iFr,iLay) = daaTemp(iFr,iLay)/raAmt(iLay) 
              END DO
            END DO
          END IF

        !!this is the temperature jacobian
        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            der1=-(296/raTemp(iLay)**2)* 
     $           tanh(kPlanck2*raFreq(iFr)/2/raTemp(iLay))
            der1=der1-(296*kPlanck2*raFreq(iFr)/2/(raTemp(iLay)**3))/
     $           (cosh(kPlanck2*raFreq(iFr)/2/raTemp(iLay))**2)
            der1=der1*daaCon(iFr,iLay)

            der2=(296/raTemp(iLay))* 
     $           tanh(kPlanck2*raFreq(iFr)/2/raTemp(iLay))
            der2=der2*daaDT(iFr,iLay)

            daaDT(iFr,iLay)=raMult(iLay)*raFreq(iFr)*(der1+der2)
            END DO
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, linear in temperature
c does temperature interpolation for SELF and for foreign
c as this is the new CKD UMBC thingies
      SUBROUTINE ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
        END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=1,kProfLayer
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
        END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)

              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
 
              END DO
            END DO
          END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)

              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!
  
              END DO
            END DO
          END IF
        END IF

      RETURN
      END
c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, linear in temperature
c only does temperature interpolation for SELF, and not for foreign
c as this was the CKD 0,2.1,2.3,2.4 thingies
      SUBROUTINE ComputeCKD_Linear_March2002(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
        END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=1,kProfLayer
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
        END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
  
              END DO
            END DO
          END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
              daaDT(iFr,iLay)=0.0d0  !this is temp jacobian!
  
              END DO
            END DO
          END IF
        END IF

      RETURN
      END
c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, quadratic in temperature
      SUBROUTINE ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
        END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=1,kProfLayer
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
        END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y1=daaCKD(iL-1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL-1,iF+1)-daaCKD(iL-1,iF))/df
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
              x1 = daTemprt(iL-1)
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)
              z1=y1-y3
              z2=y2-y3
              t1=x1*x1-x3*x3
              t2=x1-x3
              t3=x2*x2-x3*x3
              t4=x2-x3
              b=(z1*t3-z2*t1)/(t3*t2-t1*t4)
              a=(z2-b*t4)/t3
              c=y3-a*x3*x3-b*x3

              x=raTemp(iLay)
              daaCon(iFr,iLay)=a*x*x + b*x + c     !this is temp dependance!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
  
              END DO
            END DO
          END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y1=daaCKD(iL-1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL-1,iF+1)-daaCKD(iL-1,iF))/df
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
              x1 = daTemprt(iL-1)
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)
              z1=y1-y3
              z2=y2-y3
              t1=x1*x1-x3*x3
              t2=x1-x3
              t3=x2*x2-x3*x3
              t4=x2-x3
              b=(z1*t3-z2*t1)/(t3*t2-t1*t4)
              a=(z2-b*t4)/t3
              c=y3-a*x3*x3-b*x3

              x=raTemp(iLay)
              daaCon(iFr,iLay)=a*x*x + b*x + c     !this is temp dependance!
              daaDT(iFr,iLay) =2*a*x + b           !this is temp jacobian!

              END DO
            END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=1,kProfLayer
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
              daaDT(iFr,iLay)=0.0d0  !this is temp jacobian!
  
              END DO
            END DO
          END IF
        END IF

      RETURN
      END
c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses spline interpolations ........ blooody slow
      SUBROUTINE ComputeCKD_SlowSpline(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp,daC(kTempCKD)
      DOUBLE PRECISION dyp1,dypn,daY2(kTempCKD),daWork(kTempCKD)

C     Assign values for interpolation
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
        END DO   

      DO iFr=1,kMaxPts
        iF=iaFrIndex(iFr)   !index of closest freq lower than raFreq(iFr)
        DO iL=1,kTempCKD
          daC(iL)=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          END DO
        CALL DSPLY2(daTemprt,daC,kTempCKD,dYP1,dYPN,daY2,daWork) 
        DO iL=1,kProfLayer
          CALL DSPLIN(daTemprt,daC,daY2,kTempCKD,raTemp(iL),daaCon(iFr,iL)) 
          END DO
        END DO

      RETURN
      END
c************************************************************************

c this subroutine calls the routines to read in the k-compressed data
c iGasID = 1 (WATER), iFileStartFr identifies the frequency range
c have to send in ref temp, amount and profile temp,amount
c the spline is done wrt partial pressure, while the scaling is done
c wrt partial pressure
      SUBROUTINE water(iGasID,iFileStartFr,iTag,iProfLayer,iL,iU,
     $      raPAmt,raRAmt,raPPart,raRPart,raPTemp,raRTemp,
     $      iErr,iDoDQ,pProf,iProfileLayers,
     $      daaDQ,daaDT,daaAbsCoeff)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
c iGasID     = GASID ==1 for water
c iFileStartFr    = current k-comp block of 25 cm-1 that is being processed
c iTag       = current k-comp block of 25 cm-1 that is being processed
c iProfLayer = number of layers in profile === kProfLayer
c iL,iU      = min/max layer number (=1,kMaxlayer)
c daaAbs     = final uncompressed abs coefficient for gas iGasID
c iErr       = errors (mainly associated with file I/O, could be associated
c              with incorrect number of layers in compresse database etc)
c raP/RAmt   = arrays containing actual/reference gas amounts
c raP/RPart  = arrays containing actual/reference gas partial pressures
c raP/RTemp  = arrays containing actual/reference gas temperatures
c daaDQ      = analytic Jacobian wrt water amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ      = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      INTEGER iGasID,iErr,iFileStartFr,iProfLayer,iL,iU,iDoDQ,iTag
      INTEGER iProfileLayers
      REAL raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)

c local variables associated with uncompressing the water database files
      CHARACTER*80 caFName
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
      INTEGER iT0,iaTsort(kMaxTemp)
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp)
      DOUBLE PRECISION daaaKX1(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX2(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX3(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX4(kMaxK,kMaxTemp,kMaxLayer),
     $                 daaaKX5(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK)

      IF (iGasID .NE. 1) THEN
        write(kStdErr,*) 'Expecting to read in water profile'
        iErr=1
        CALL DoSTOP
        END IF

      iIOUN=kCompUnit
      CALL CompFileName(+1,iGasID,iFileStartFr,iTag,caFName)
      CALL rdcompwater(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,
     $        iNLay,iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort,
     $        daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUX)

c check that the file has the data for the correct gas
      IF (iFileGasID .NE. iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
 1000   FORMAT('Error! file : ',/,A80,/,
     $         'contains data for Gas ',i2,' not desired gas',I2)
        CALL DoSTOP
        END IF

c check that the data file has the right number of layers ===== AIRS layers
      IF (iNLay .NE. kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
 1010   FORMAT('Error! file : ',/,A80,/,
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
     $       iNk,iKm,iKn,iUm,iUn,daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers)
      ELSE
          CALL GetAbsCoeffWaterOLD(daaAbsCoeff,daToffset,
     $      daaaKX1,daaaKX2,daaaKX3,daaaKX4,daaaKX5,daaUx,raPTemp,
     $      raRTemp,raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers)
          END IF

c because iKtype=2 do any necessary jacobians calcs HERE!
        IF (kJacobian .GT. 0) THEN
          IF (iDoDQ .GT. 0)  THEN
            CALL FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDQ,raPAmt)
            END IF
          CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
          END IF
          
      ELSE      
        !kGenln2Water .LT. 0  ==> do same uncompression as for CO2
        !ie do not worry about the slef broadening corrections
        !this is not very good at all!
c interpolate compressed data in temperature, to get abs coeff matrix
        IF (kJacobian .GE. 0) THEN
          CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx2,daaUx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $      daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers)
        ELSE
          CALL GetAbsCoeffOLD(daaAbsCoeff,daToffset,daaaKx2,daaUx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers)
          END IF

c because of iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
        IF (kJacobian .GE. 0) THEN
          CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
          IF (iDoDQ .GT. 0) THEN
            CALL FinalAmtDeriv(daaDQ,iKtype)
            END IF
          END IF
        END IF

c convert absorption coefficient correctly if necessary
      IF (iKtype .eq. 2) THEN
        CALL RaisePower(daaAbsCoeff)
        END IF

c now compute optical depth = gas amount * abs coeff
      CALL AmtScale(daaAbsCoeff,raPAmt)

cc Scott will tune the reflected thermal in the Fast Model himself, so only  
cc put fudges in water continuum!  
cc Eventually we will hve to modify line strengths in the HITRAN database 
cc notice how we only need to fudge fix 2380 and above chunks here!
cc      IF (iGasID .EQ. 1) THEN
cc        IF (iFileStartFr .GE. 2380) THEN
cc          CALL water_4um_fudge(daaAbsCoeff,iFileStartFr)
cc          END IF
cc        END IF

      RETURN
      END

c************************************************************************
c this subroutine calls the routines to read in the k-compressed data
c iGasID tells which gas type, iFileStartFr identifies the frequency range
c have to send in ref temp, amount and profile temp,amount
      SUBROUTINE othergases(iGasID,iFileStartFr,iTag,iProfLayer,iL,iU,
     $      raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers,
     $      daaDQ,daaDT,daaAbsCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
c iGasID     = GASID
c iFileStartFr    = current k-comp block of 25 cm-1 that is being processed
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
      INTEGER iGasID,iErr,iFileStartFr,iProfLayer,iL,iU,iDoDQ,iTag
      INTEGER iProfileLayers
      REAL raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
    
c local variables associated with uncompressing the database
      CHARACTER*80 caFName
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
      INTEGER iT0,iaTsort(kMaxTemp)
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp)
      DOUBLE PRECISION daaaKX(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK)

      iIOUN=kCompUnit
      CALL CompFileName(+1,iGasID,iFileStartFr,iTag,caFName)
      CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay,
     $              iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort,
     $              daaaKX,daaUX)

c check that the file has the data for the correct gas
      IF (iFileGasID .NE. iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
 1000   FORMAT('Error! file : ',/,A80,/,
     $         'contains data for Gas ',i2,' not desired gas',I2)
        CALL DoSTOP
        END IF

c check that the data file has the right number of layers
      IF (iNLay .NE. kMaxLayer) THEN
        iErr=1
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer
 1010   FORMAT('Error! file : ',/,A80,/,
     $         'contains data for ',i3,' layers but kMaxLayer = ',I3)
        CALL DoSTOP
        END IF

c interpolate compressed data in temperature, to get abs coeff matrix
      IF (kJacobian .GE. 0) THEN
        CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $  raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers)
      ELSE
        CALL GetAbsCoeffOLD(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $         raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $         pProf,iProfileLayers)
        END IF

c because of the iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
      IF (kJacobian .GE. 0) THEN
        CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
        IF (iDoDQ .GT. 0) THEN
          CALL FinalAmtDeriv(daaDQ,iKtype)
          END IF
        END IF

c convert absorption coefficient correctly if necessary
      IF (iKtype .eq. 2) THEN
        CALL RaisePower(daaAbsCoeff)
        END IF

c now compute optical depth = gas amount * abs coeff
      CALL AmtScale(daaAbsCoeff,raPAmt)

c notice how we only need to fudge fix 2255,2280 and 2305,2405 chunks here!
      IF (iGasID .EQ. 2) THEN 
        IF ((iFileStartFr .EQ. 2255) .OR. (iFileStartFr .EQ. 2280) .OR.
     $      (iFileStartFr .EQ. 2380) .OR. (iFileStartFr .EQ. 2405)) THEN 
          CALL co2_4um_fudge(daaAbsCoeff,iFileStartFr) 
          END IF 
        END IF 
 
      RETURN
      END
c************************************************************************

c this subroutine scales the CO2 absorption coefficients in the 2355,2280 
c and 2380,2405 chunks
c look at chi.f for the fudge for 2380,2405 cm-1
      SUBROUTINE Co2_4um_fudge(daaAbsCoeff,iFileStartFr)

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c iFileStartFr = which chunk
      INTEGER iFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iIOUN,iERR,iChi
      REAL raF(kMaxPts),raChi(kMaxPts),m,c,rF,rX,rChi,r1,r2,r3,rAmp,rAmp0
      CHARACTER*120 FNAME

c the blah.txt  files are the orig,      done in late june/early july 2002
c the blah_a.txt files are refinement 1, done in early aug 2002
c the blah_b.txt files are refinement 2, done in mid nov 2002

      iChi = -1
      IF (iFileStartFR .EQ. 2255) THEN
        write(kStdWarn,*) 'need CO2 chifunction for 2255 chunk ....'
        FNAME = 'co2_4um_fudge_2255.txt'
        FNAME = 'co2_4um_fudge_2255_a.txt'
        iChi = +1
      ELSEIF (iFileStartFR .EQ. 2280) THEN
        write(kStdWarn,*) 'need CO2 chifunction for 2280 chunk ....'
        FNAME = 'co2_4um_fudge_2280.txt'
        FNAME = 'co2_4um_fudge_2280_a.txt'
        iChi = +1
c      ELSEIF (iFileStartFR .EQ. 2355) THEN
c        write(kStdWarn,*) 'need CO2 chifunction for 2355 chunk ....'
c        FNAME = 'co2_4um_fudge_2355.txt'
c        iChi = +1
      ELSEIF (iFileStartFR .EQ. 2380) THEN
        write(kStdWarn,*) 'need CO2 chifunction for 2380 chunk ....'
        FNAME = 'co2_4um_fudge_2380.txt'
        FNAME = 'co2_4um_fudge_2380_b.txt'
        iChi = +1
      ELSEIF (iFileStartFR .EQ. 2405) THEN
        write(kStdWarn,*) 'need CO2 chifunction for 2405 chunk ....'
        FNAME = 'co2_4um_fudge_2405.txt'
        FNAME = 'co2_4um_fudge_2405_b.txt'
        iChi = +1
        END IF

      CALL FindChiFileName(fname)

      IF (iChi .GT. 0) THEN
        iIOUN=kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine co2_4um_fudge'
          WRITE(kStdErr,1010) IERR, FNAME
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
          ENDIF
        kTempUnitOpen=1
        READ(iIOUN,*) (raF(iFr),raChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen=-1

        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*raChi(iFr)
            END DO
          END DO
      ELSE
        write(kStdWarn,*) 'do not need CO2 chifunction for chunk ',iFileStartFr
        END IF        

      RETURN
      END
     
c************************************************************************
c this subroutine scales water absorption coefficients in the >= 2380 chunks
      SUBROUTINE water_4um_fudge(daaAbsCoeff,iFileStartFr)

      include '../INCLUDE/kcarta.param'

c daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
c iFileStartFr = which chunk
      INTEGER iFileStartFr
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iLay,iIOUN,iERR,iChi,i1,i2,iFileEnd
      REAL raF(kMaxPts),raChi(kMaxPts)
      DOUBLE PRECISION dX,daX(kMaxPts),dSlope,xA,x

      iChi = -1
      IF (iFileStartFR .EQ. 2380) THEN
        write(kStdWarn,*) 'need water chifunction for ',iFileStartFR,' 
     $ chunk ....'
        iChi = +2
      ELSEIF (iFileStartFR .GE. 2405) THEN
        write(kStdWarn,*) 'need water chifunction for ',iFileStartFR,' 
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
        write(kStdWarn,*) 'do not need H2O chifunction for chunk ',iFileStartFr
        END IF        

      RETURN
      END
     
c************************************************************************

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
      DOUBLE PRECISION daaDA(kMaxPts,kProfLayer)

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
c or Zget = (1/4) K(v,T)^(-3/4) dK/dT = daaDA
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
      DOUBLE PRECISION daaDA(kMaxPts,kProfLayer)

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
       CHARACTER*80 FNAM
       INTEGER iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0,
     $    ITSORT(kMaxTemp)
       DOUBLE PRECISION SFREQ,FSTEP,TOFF(kMaxTemp),
     $    KX1(kMaxK,kMaxTemp,kMaxLayer),
     $    KX2(kMaxK,kMaxTemp,kMaxLayer),
     $    KX3(kMaxK,kMaxTemp,kMaxLayer),
     $    KX4(kMaxK,kMaxTemp,kMaxLayer),
     $    KX5(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

C      Local variables
       INTEGER IERR,I,J,K,RESORT,IHOLD
C-----------------------------------------------------------------------
c Larry M. suggested ACTION='READ' option, but SGI does not seem to support this
       OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine RDCOMPWATER'
          WRITE(kStdErr,1010) IERR, FNAM
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
       ENDIF
       kCompUnitOpen=1

C      Read in the header
       READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN,
     $    UM, UN

C      Make sure the array sizes are <= to the declared sizes
       IF (NK .GT. kMaxK) THEN
          WRITE(kStdErr,1110)
 1110     FORMAT('Error! Compressed data array dimension exceeds ',
     $       'max size')
          WRITE(kStdErr,1120) NK, kMaxK
 1120     FORMAT('NK = ',I3,', kMaxK = ',I3)
          CALL DoSTOP
       ENDIF

       IF (KT .GT. kMaxTemp) THEN
          WRITE(kStdErr,1110)
          WRITE(kStdErr,1130) KT, kMaxTemp
 1130     FORMAT('KT = ',I2,', kMaxTemp = ',I2)
          CALL DoSTOP
       ENDIF

       IF (KN .GT. kMaxLayer) THEN
          WRITE(kStdErr,1110)
          WRITE(kStdErr,1140) KN, kMaxLayer
 1140     FORMAT('KN = ',I3,', kMaxLayer = ',I3)
          CALL DoSTOP
       ENDIF

       IF (UM .GT. kMaxPts) THEN
          WRITE(kStdErr,1110)
          WRITE(kStdErr,1150) UM, kMaxPts
 1150     FORMAT('UM = ',I5,', kMaxPts = ',I5)
          CALL DoSTOP
       ENDIF

       IF (UN .GT. kMaxK) THEN
          WRITE(KSTDERR,1110)
          WRITE(KSTDERR,1160) UN, kMaxK
 1160     FORMAT('UN = ',I3,', kMaxK = ',I3)
          CALL DoSTOP
       ENDIF

C      Read in the temperature offsets
       READ(iIOUN) (TOFF(I),I=1,KT)

C      Find which of the offsets is 0 (error if none).
       IT0=0
       DO I=1,KT
          IF (TOFF(I) .EQ. 0.0) IT0=I
          ITSORT(I)=I
       ENDDO
       IF (IT0 .EQ. 0) THEN
          WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
 1180     FORMAT('ERROR! One of the temperature offsets must be 0',/,
     $       'offsets =',20(' ',F5.1))
          CALL DoSTOP
       ENDIF
C
C      Sort the indices of the temperature offsets in ascending order
       RESORT=1
 10    IF (RESORT .EQ. 1) THEN
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

C      Read in the U matrix
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
c KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
c UX       = uncompression matrix
c ITO,ITSORT = base points to do temperature interpolation
C      Calling parameters
       CHARACTER*80 FNAM
       INTEGER iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0,
     $    ITSORT(kMaxTemp)
       DOUBLE PRECISION SFREQ,FSTEP,TOFF(kMaxTemp),
     $    KX(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)
C
C      Local variables
       INTEGER IERR,I,J,K,RESORT,IHOLD
C-----------------------------------------------------------------------
C
       OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=IERR)
       IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine RDCOMP'
          WRITE(KSTDERR,1010) IERR, FNAM
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
       ENDIF
       kCompUnitOpen=1
C
C      Read in the header
       READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN,
     $    UM, UN
C
C      Make sure the array sizes are <= to the declared sizes
       IF (NK .GT. kMaxK) THEN
          WRITE(KSTDERR,1110)
 1110     FORMAT('Error! Compressed data array dimension exceeds ',
     $       'max size')
          WRITE(KSTDERR,1120) NK, kMaxK
 1120     FORMAT('NK = ',I3,', kMaxK = ',I3)
          CALL DoSTOP
       ENDIF
C
       IF (KT .GT. kMaxTemp) THEN
          WRITE(KSTDERR,1110)
          WRITE(KSTDERR,1130) KT, kMaxTemp
 1130     FORMAT('KT = ',I2,', kMaxTemp = ',I2)
          CALL DoSTOP
       ENDIF
C
       IF (KN .GT. kMaxLayer) THEN
          WRITE(KSTDERR,1110)
          WRITE(KSTDERR,1140) KN, kMaxLayer
 1140     FORMAT('KN = ',I3,', kMaxLayer = ',I3)
          CALL DoSTOP
       ENDIF
C
       IF (UM .GT. kMaxPts) THEN
          WRITE(KSTDERR,1110)
          WRITE(KSTDERR,1150) UM, kMaxPts
 1150     FORMAT('UM = ',I5,', kMaxPts = ',I5)
          CALL DoSTOP
       ENDIF
C
       IF (UN .GT. kMaxK) THEN
          WRITE(KSTDERR,1110)
          WRITE(KSTDERR,1160) UN, kMaxK
 1160     FORMAT('UN = ',I3,', kMaxK = ',I3)
          CALL DoSTOP
       ENDIF
C
C      Read in the temperature offsets
       READ(iIOUN) (TOFF(I),I=1,KT)
C
C      Find which of the offsets is 0 (error if none).
       IT0=0
       DO I=1,KT
          IF (TOFF(I) .EQ. 0.0) IT0=I
          ITSORT(I)=I
       ENDDO
       IF (IT0 .EQ. 0) THEN
          WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
 1180     FORMAT('ERROR! One of the temperature offsets must be 0',/,
     $       'offsets =',20(' ',F5.1))
          CALL DoSTOP
       ENDIF
C
C      Sort the indices of the temperature offsets in ascending order
       RESORT=1
 10    IF (RESORT .EQ. 1) THEN
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

cWARNING!!!!!!!!!!!!! 
C      Read in the K matrices
c       DO I=1,NK
c          READ(iIOUN) ((KX(I,J,K),J=1,KT),K=1,KN)
c       ENDDO

       DO I=1,NK
         READ(iIOUN) ((KX(I,J,K),K=1,KN),J=1,KT)
         ENDDO

C
C      Read in the U matrix
       DO I=1,NK
          READ(iIOUN) (UX(J,I),J=1,NPTS)
       ENDDO
C
       CLOSE(iIOUN)
       kCompUnitOpen=-1

       RETURN
       END

c************************************************************************
c************************************************************************
c************************************************************************
c**  this is the old code (cross sections and continuum gateway calls) **
c************************************************************************
c************************************************************************
c************************************************************************
c this subroutine is the gateway call to XSEC/calxsc 
c compute the contribution of Gases 29-63 (if present) 
      SUBROUTINE CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep, 
     $      daaTemp,raVTemp,iVTSet,raFreq,iErr,caXsecName, 
     $      raTAmt,raTTemp,raTPress,raTPart,iaCont, 
     $      daaDQ,daaDT,iDoDQ) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
c iCount    = which of the iNumGases is being processed 
c iGasID    = iaGasID(iCount) = gas ID of current gas 
c iRefLayer = number of layers in the reference profiles (=kProfLayer) 
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer) 
c iaCont    = whether or not to do continuum calculation .. iaCont(iCount) 
c caXecF    = file name of cross section data 
c daaTemp   = matrix containing the uncompressed k-spectra 
c raVtemp   = vertical temperature profile for the Mixed paths 
c iVTSet    = has the vertical temp been set, to check current temp profile 
c iFileStartFr   = which k-comp file chunk to uncompress 
c raFreq    = wavenumber array 
c daaDQ     = analytic jacobian wrt gas amount 
c daaDT     = analytic jacobian wrt temperature 
c iErr      = errors (mainly associated with file I/O) 
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs 
c kFrStep  = kFreqStep 
      INTEGER iaCont(kMaxGas),iDoDQ 
      CHARACTER*80 caXsecName 
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer) 
      INTEGER iGasID,iL,iU,iErr,iRefLayer,iCount,iVTSet 
      REAL raTPress(kProfLayer),raTPart(kProfLayer),kFrStep 
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer) 
      REAL raFreq(kMaxPts),raVTemp(kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac) 
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac) 
 
c local variables 
      REAL raXamnt(kProfLayer),raAdjustProf(kProfLayer),rMin 
      REAL raaXsec(kMaxPts,kProfLayer),rFStep,rCheckTemp 
      INTEGER iNmol, iaMolid(MXXMOL),iLay,iFr,iWhichUsed
 
c used in d/dq, d/dT if kSpline = 1 
      REAL raaSplineDQ(kMaxPtsJac,kProfLayerJac) 
      REAL raaSplineDT(kMaxPtsJac,kProfLayerJac) 
      INTEGER idQT 
 
      REAL AVOG 
      DATA AVOG/6.022045E+26/ 
 
c to see if all the abs coefficients are less than zero ...tau=exp(-abs) <=1 
      REAL rEC,rECCount 
 
c first check to see if exact calculations of d/dq,d/dT should be enabled 
      IF (kJacobian .GT. 0) THEN  
        idQT=1 
      ELSE 
        idQT=-1 
        END IF 

      iLay=kProfLayer 
      iFr=kMaxPts 
      rFStep=kFrStep 
 
c only calculate the current gas cross-section 
      DO iLay=1,MXXMOL 
        iaMolid(iLay)=0 
        END DO 
      iNmol=1 
      iaMolid(1)=iGasID 
 
      iLay=kProfLayer 
 
      IF (iErr .LE. 0) THEN 
c set the vertical temperature profile if iCount=1 (first profile read) 
        IF ((iVTSet .LT. 0) .AND. (iGasID .LE. kGasComp)) THEN  
          write(kStdWarn,*) 'Setting vertical temp profile ...' 
          DO iLay=1,kProfLayer 
            raVTemp(iLay)=raTTemp(iLay) 
            END DO 
          END IF 
c if previous profiles have been read in, check to make sure the  
c temperature profiles are the same!!!! 
        IF ((iVTSet .GT. 0) .AND. (iGasID .LE. kGasComp)) THEN  
          write(kStdWarn,*) 'Checking the vertical temp profile ...' 
          DO iLay=1,kProfLayer 
            rCheckTemp=raTTemp(iLay)-raVTemp(iLay) 
            IF (abs(rCheckTemp) .GE. 1.0e-3) THEN 
              write(kStdWarn,*) 'Warning!!Temp profiles do not match!!!' 
              write(kStdWarn,*) 'Gas#,layer, gastemp, vertical temp = ' 
              write(kStdWarn,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay) 
              write(kStdErr,*) 'Warning!!Temp profiles do not match!!!' 
              write(kStdErr,*) 'Gas#,layer, gastemp, vertical temp = ' 
              write(kStdErr,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay) 
              END IF 
            END DO 
          END IF 
        END IF 
 
      IF (iErr .LT. 0) THEN 
c convert amount from k.moles/cm^2 to molecules/cm^2 
        DO iLay=1,kProfLayer 
          raXamnt(iLay)=raTAmt(iLay)*AVOG 
          raAdjustProf(iLay)=raTTemp(iLay) 
          ENDDO 
 
c make sure that the relevant area in the cross section matrix is zeroed 
        DO iLay=1,kProfLayer 
          DO iFr=1,kMaxPts 
            raaXsec(iFr,iLay)=0.0 
            END DO 
          END DO 
 
c if we need to do jacobian calculations using splines, 
c initialize the matrices here  
        IF (kJacobian .GT. 0) THEN 
          IF (iDoDQ .GT. 0) THEN 
            DO iLay=1,kProfLayerJac 
              DO iFr=1,kMaxPtsJac 
                raaSplineDQ(iFr,iLay)=0.0 
                raaSplineDT(iFr,iLay)=0.0 
                END DO 
              END DO 
          ELSE IF (iDoDQ .LT. 0) THEN 
            DO iLay=1,kProfLayerJac 
              DO iFr=1,kMaxPtsJac 
                raaSplineDT(iFr,iLay)=0.0 
                END DO 
              END DO 
            END IF   
          END IF   
 
 
c this call calculates ONLY the cross sections  
c            OR 
c as idQT=1 this call calculates the cross sections, AND d/dq,d/dT!!!, saving 
c them in raaSplineDQ ,raaSplineDT respectively 

        iLay=kProfLayer 
        iFr=kMaxPts 
        rFStep=kFrStep 

        CALL CALXSC(caXsecName,iFr,raFreq,rFStep, 
     $      iLay, raAdjustProf,raXamnt,iNmol,iaMolid, 
     $      raaXsec,iWhichUsed, 
     $      idQT, raaSplineDQ ,raaSplineDT, iGasID, iDoDQ) 
 
c since CALXSC computed the cross-section and saved it in a 3d matrix, 
c pull out the relevant cross section data into a 2d matrix and add it to  
c the cumulative abscoeff matrix daaTemp 
c first check whether the CALXSC routine actually found a matching  
c gas+wavenumber region 
        rMin=1.0e10 
        rEC=0.0 
        rECCount=0.0 
        IF ((iWhichUsed .GT. 0) .AND. (iWhichUsed .LE. iNmol)) THEN 
          WRITE(kStdWarn,1000) iCount,iGasID 
 1000          FORMAT('adding on cross section for GAS(',I2,') = ',I2) 
          DO iLay=1,kProfLayer 
            DO iFr=1,kMaxPts 
              IF (raaXsec(iFr,iLay) .LT. 0.0) THEN 
c flag this error, as the cross section values should all be > 0 
                iErr=1 
                rEC=rEC+raaXsec(iFr,iLay) 
                rECCount=rECCount+1 
                IF (raaXsec(iFr,iLay) .LT. rMin) THEN 
                  rMin = raaXsec(iFr,iLay) 
                  END IF 
                raaXsec(iFr,iLay)=0.0 
                END IF 
              daaTemp(iFr,iLay)=raaXsec(iFr,iLay) 
              END DO 
            END DO 
          ELSE 
            WRITE(kStdWarn,1010) iCount,iGasID 
 1010              FORMAT ('no XSEC contribution due to GAS(',I2,') = ',I2) 
          END IF 
        IF (rECCount .gt. 0.5) THEN 
          rEC=rEC/rECCount 
          write(kStdWarn,*) 'Error in XSEC data!!! Some values negative!' 
          write(kStdWarn,*) 'and reset to 0.0' 
          write(kStdWarn,*) 'rECCount values have avg value ',rEC 
          write(kStdWarn,*) 'min negative value in 10000*100 = ',rMin 
          IF (abs(rMin) .GT. 1.0e-7) THEN 
            iErr=1 
            CALL DoSTOP 
            END IF 
          END IF           
c end main if statement 
        END IF 
 
c do the inclusion of the exact derivatives here 
      IF (kJacobian .GT. 0) THEN 
c exact calculations have already been performed in calXSC,calq ...  
c so apart from the avog factor in raaDQ, no need to do much here!!! 
        IF (iDoDQ .GT. 0) THEN           
          DO iLay=1,kProfLayerJac 
            DO iFr=1,kMaxPtsJac 
              daaDQ(iFr,iLay)=raaSplineDQ(iFr,iLay)*avog 
              daaDT(iFr,iLay)=raaSplineDT(iFr,iLay) 
              END DO 
            END DO 
        ELSE IF (iDoDQ .LT. 0) THEN 
          DO iLay=1,kProfLayerJac 
            DO iFr=1,kMaxPtsJac 
              daaDT(iFr,iLay)=raaSplineDT(iFr,iLay) 
              END DO 
            END DO 
          END IF 
        END IF 
         
      RETURN 
      END 
 
c************************************************************************ 
