c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c this file reads *MOLGAS,*XSCGAS,*FREQ (very easy),*MIXFIL

c************************************************************************
c this subroutine deals with the 'MOLGAS' keyword
c and basically deals with GasID 1--28
c
c main output parameter is iaMOLgases
c synopsis IN : iNGas    = number of gasIDs to be read in (if -1, all)
c               iaGeseNL = list of gasIDs
c         OUT : iaMOLgases(iI) = set to 1 if gasID was found in this section
      SUBROUTINE molgas4(iNGas,iaGasesNL,iaMOLgases)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNGas     = number of gases stored in this namelist
c iaGasesNL = gasIDs stored in this namelist
c             the elements contain the gasIDs, or -1 ow
c                so eg if gases1,2,5,101 needed then
c            iaGasesNL(1)=1  iaGasesNL(2)=2  iaGasesNL(3)=5 iaGasesNL(4)=101 
c iaMOLgases   = overall gasIDs to be used by kCARTA
c                the elements are set to 0 if gas not required, 1 if required
c                so eg if gases1,2,5,101 needed then
c            iaMOLgases(1)=iaMOLgases(2)=iaMOLgases(5)=iaMOLgases(101)=1
c            all else 0
c MOLGAS or XSCGAS ... array iaInputOrder is built up from this array, each
c time it is returned from a call to either subroutine MOLGAS or xscfil
      INTEGER iNGas,iaGasesNL(kGasComp),iaMOLgases(kMaxGas)

c local variables
      INTEGER iTag,iErr,iNotMainGas,MainGas
      INTEGER iInt,iaTemp(kMaxGas),iC,iCC,iaInDataBase(kGasComp)
      INTEGER iNgasesCheck
      INTEGER iCheckCompDataBase
      CHARACTER*7 caWord

      caWord='*MOLGAS'

      DO iInt=1,kMaxGas
        iaMOLgases(iInt)=0
        iaTemp(iInt)=0
        END DO

c read the no of gases stored in the file
      iErr=-1      
      IF (iNgas .GT. 0) THEN
c use the molecular ID's from the namelist file 
        DO iInt=1,iNgas
          iNotMainGas=MainGas(iaGasesNL(iInt))
          IF (iNotMainGas .LT. 0) THEN
            write(kStdErr,*) 'Invalid MOLGAS GasID',iaGasesNL(iInt),' entered'
            write(kStdErr,*) 'Please check *MOLGAS and retry'
            write(kStdErr,*) 'Note : If the GasID printed was -100, you '
            write(kStdErr,*) 'probably entered less gasIDs than you promised '
            CALL DoSTOP 
          ELSE
            iaTemp(iInt)=iaGasesNL(iInt)
            END IF
          END DO      
      ELSE IF (iNgas .LT. 0) THEN
c use all gases in the compressed database
c  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9
c         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,X PLUS kSelf,kFor
        iNgas=0
        DO iC=1,kGasComp
          iaInDataBase(iC)=-1
          END DO
        DO iC=1,kGasComp
          iTag=-1
          iCC=iCheckCompDataBase(iC,-100.0,-100.0,iTag,iErr)
          IF (iCC .GT. 0) THEN
            iNgas=iNgas+1
            iaInDataBase(iC)=1
            END IF
          END DO
c now based on which gases were found, reset array iaTemp
        iCC=1
        DO iC=1,kGasComp
          IF (iaInDataBase(iC) .GT. 0) THEN
            write(kStdWarn,*)'Including gasID ',iC,' in comp database'
            iaTemp(iCC)=iC
            iCC=iCC+1
            END IF
          END DO

        IF (kCKD .GE. 0) THEN
          DO iC=kNewGasLo,kNewGasHi
            write(kStdWarn,*)'Including gasID ',iC,' in comp database'
            iaTemp(iCC)=iC
            iCC=iCC+1
            iNgas=iNgas+1
            END DO
          END IF

        END IF

c check the molecular ID's are between 1 and kGasComp , and
c kNewGasLo and kNewGasHi
c (should be in the compressed data base)
c if iNgas = -1, everything should have been set correctly above; if user
c enetered in the GasIDs him/her self, there could be mistakes
      DO iC=1,iNgas
        iNotMainGas=MainGas(iaTemp(iC))
          IF (iNotMainGas .LT. 0) THEN
c gas does not exist in the compressed base ... stop
          WRITE(kStdErr,777) iaTemp(iC),1,kGasComp
  777       FORMAT('Error in MOLGAS!! found Gas ID  = ',I2,'(
     $ MOLGAS ID''s should be between ',I2,' and ',I2 ,')')
          CALL DoSTOP
          END IF
c check to see if kcomp files do exist
        iTag=-1
        IF (iaTemp(IC) .LE. kGasComp) THEN
          iCC=iCheckCompDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
          IF (iCC .LT. 0) THEN
            WRITE(kStdWarn,780) iaTemp(iC)
            END IF
        ELSE
          WRITE(kStdWarn,785) iaTemp(iC)
          END IF
        END DO
 780  FORMAT('Warning! GasID ',I3,' not in compressed data base')
 785  FORMAT('Warning! GasID ',I3,' is a new continuum gas')
      
      iNgasesCheck=iNgas

c set the identities of the gases whose abs spectra are in files
c iaMOLgases keeps track of how many times the GasID has been found, so that
c no double counting of the gases (between MOLGAS and XSCGAS) is done
      DO iInt=1,iNgas
        iaMOLgases(iaTemp(iInt))=iaMOLgases(iaTemp(iInt))+1
        END DO

c check to see that the same gas ID has not been entered more than once
      DO iInt=1,iNgas
        IF ((iaMOLgases(iaTemp(iInt)).GT.1) .OR. 
     $      (iaMOLgases(iaTemp(iInt)).LT.0)) THEN
          write(kStdErr,*) 'Gas ID',iaTemp(iInt),' entered more than'
          write(kStdErr,*)'once. Please check *MOLGAS and retry'
          CALL DoSTOP 
          END IF
        END DO

      iNgas=iNgasesCheck
      write(kStdWarn,*) 'MOLGAS ... gases stored are ',iNgas
      DO iInt=1,kGasComp
        IF (iaMOLgases(iInt) .GT. 0) THEN
          write(kStdWarn,*) '     going to use compressed database gas ',iInt
          END IF
        END DO

      DO iInt=kNewGasLo,kNewGasHi
        IF (iaMOLgases(iInt) .GT. 0) THEN
          write(kStdWarn,*) '     going to use new continuum gas ',iInt
          END IF
        END DO

      IF (kCKD .GE. 0) THEN
        DO iInt=kNewGasLo,kNewGasHi
          IF (iaMOLgases(iInt) .LE. 0) THEN
            write(kStdWarn,*) 'Cannot have CKD on and gasIDs 101/102 unused'
            write(kStdErr,*) 'Cannot have CKD on and gasIDs 101/102 unused'
            Call DoSTOP
            END IF
          END DO
        END IF
        
      RETURN
      END

c************************************************************************

c this subroutine deals with the 'XSCGAS' keyword
c and basically deals with GasID 51--63
c synopsis IN : iNXsec    = number of gasIDs to be read in (if -1, all)
c               iaLXsecNL = list of gasIDs
c         OUT : iaXSCgases(iI) = set to 1 if gasID was found in this section
      SUBROUTINE xscgas4(iNXsec,iaLXsecNL,iaXSCgases)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNXSec     = number of gases stored in this namelist
c iaLXGasesNL = gasIDs stored in this namelist
c             the elements contain the gasIDs, or -1 ow
c                so eg if gases 51,52 needed then
c            iaLXSecNL(1)=51  iaLXSecNL(2)=52 
c iaXSCgases   = overall gasIDs to be used by kCARTA
c                the elements are set to 0 if gas not required, 1 if required
c                so eg if gases 51,52 needed then
c            iaXSCgases(51)=iaXSCgases(52)=1
c            all else 0
c MOLGAS or XSCGAS ... array iaInputOrder is built up from this array, each
c time it is returned from a call to either subroutine MOLGAS or xscfil
      INTEGER iNXsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1),iaXSCgases(kMaxGas)

c local variables
      CHARACTER*7 caWord
      INTEGER iInt,iErr,iC,iCC,iCheckXsecDataBase,iNXsecCheck
      INTEGER iaTemp(kMaxGas),iTag
      INTEGER iaInDataBase(kProfLayer)      
      caWord='*XSCGAS'

      DO iInt=1,kMaxGas
        iaXSCgases(iInt)=0
        iaTemp(iInt)=0
        END DO

      iErr=-1
      IF (iNxsec .GT. 0) THEN
c use the xsec gas ID's from the namelist file 
        DO iInt=1,iNXsec
          IF ((iaLXsecNL(iInt) .LT. kGasXsecLo) .OR. 
     $        (iaLXsecNL(iInt) .GT. kGasXsecHi)) THEN
            write(kStdErr,*) 'Invalid XSCGAS GasID',iaLXsecNL(iInt),' entered'
            write(kStdErr,*) 'Please check *XSCGAS and retry'
            write(kStdErr,*) 'Note : If the GasID printed was -100, you '
            write(kStdErr,*) 'probably entered less gasIDs than you promised '
            CALL DoSTOP 
          ELSE
            iaTemp(iInt)=iaLXsecNL(iInt)
            END IF
          END DO      
      ELSE IF (iNxsec .LT. 0) THEN
c use all gases in the xsec database
c  check to see the following reference profiles exist : 51..63
        iNxsec=0
        DO iC=kGasXsecLo,kGasXsecHi
          iaInDataBase(iC)=-1
          END DO
        DO iC=kGasXsecLo,kGasXsecHi
          iTag=-1
          iCC=iCheckXsecDataBase(iC,-100.0,-100.0,iTag,iErr)
          IF (iCC .GT. 0) THEN
            iNxsec=iNxsec+1
            iaInDataBase(iC)=1
            END IF
          END DO
c now based on which gases were found, reset array iaTemp
        iCC=1
        DO iC=kGasXsecLo,kGasXsecHi
          IF (iaInDataBase(iC) .GT. 0) THEN
            write(kStdWarn,*) 'Including gasID ',iC,' in xsec database'
            iaTemp(iCC)=iC
            iCC=iCC+1
            END IF
          END DO
        END IF
 
c check the molecular ID's are between kGasXsecLo and kXsecGasHi
c (should be in the cross sec data base)
c if iNXsec = -1, everything should have been set correctly above; if user
c enetered in the GasIDs him/her self, there could be mistakes
      DO iC=1,iNXsec
        IF ((iaTemp(iC).LT.kGasXsecLo).OR.(iaTemp(iC).GT.kGasXsecHi)) THEN
c gas does not exist in the xsec data base ... stop
          WRITE(kStdErr,777) iaTemp(iC),kGasXsecLo,kGasXsecHi
 777      FORMAT('Error in XSCGAS!! found Gas ID  = ',I2,'(
     $ XSCGAS ID''s should be between ',I2,' and ',I2 ,')')
           CALL DoSTOP
           END IF
c check to see if cross section data files for this gas does exist
         iTag=-1
         iCC=iCheckXsecDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
         IF (iCC .LT. 0) THEN
           WRITE(kStdWarn,780) iaTemp(iC)
 780       FORMAT('Warning! GasID ',I2,' not in xsec data base')
           END IF
         END DO        

      iNXsecCheck=iNXsec
c set the identities of the gases whose molecular cross sections are in file
c iaXSCgases keeps track of how many times the GasID has been found, so that
c no double counting of the gases (between MOLGAS and XSCGAS) is done
      DO iInt=1,iNXsec
        iaXSCgases(iaTemp(iInt))=iaXSCgases(iaTemp(iInt))+1
        END DO

c check to see that the same gas ID has not been entered more than once
      DO iInt=1,iNXsec
        IF ((iaXSCgases(iaTemp(iInt)).GT.1) .OR. 
     $      (iaXSCgases(iaTemp(iInt)).LT.0)) THEN
          write(kStdErr,*) 'Gas ID ',iaTemp(iInt),' entered more than'
          write(kStdErr,*) 'once Please check *XSCGAS and retry'
          CALL DoSTOP 
          END IF
        END DO

      iNXsec=iNXsecCheck
      write(kStdWarn,*)'cross section gases stored are ',iNXsec
      DO iInt=kGasXsecLo,kGasXsecHi
        IF (iaXSCgases(iInt) .GT. 0) THEN
          write(kStdWarn,*) 'CrossSect database gas ',iInt
          END IF
        END DO

      RETURN
      END

c************************************************************************
c this subroutine is for the new spectra
c this reads in the info for the new spectra for ONE gas
      SUBROUTINE spectra4(iNumNewGases,
     $             iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNumNewGases   tells number of new gases 
c iaNewGasID     tells which gases we want to update spectroscopy 
c iaNewData      tells how many new data sets to read in for each gas 
c iaaNewChunks   tells which data chunks to read in 
c caaaNewChunks  tells the name of the files associated with the chunks 
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore)
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT) 
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT) 

c local variables
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iCount,iErr

      caWord='*SPECTR'
      iErr=-1

 5030 FORMAT(A130) 

      iNumLinesRead=0 
 13   IF (iNumLinesRead .GT. 0) THEN 
        iErr=1 
        WRITE(kStdErr,5010) caWord 
        CALL DoSTOP 
        END IF 
 5010 FORMAT('Error reading section ',A7) 

      iNumLinesRead=1
c read how many gases will have new spectroscopy
      iErr=-1  
      IF ((iNumNewGases .LT. 1) .OR. (iNumNewGases .GT. kGasStore)) THEN  
        write(kStdErr,*)'need a valid number of gases in *SPECTRA!!'  
        write(kStdErr,*)'must be > 0, < kGasStore!!'
        write(kStdErr,*)'please check and retry!' 
        CALL DoSTOP  
        END IF

      DO iCount=1,iNumNewGases
        IF ((iaNewGasID(iCount) .LT. 1) .OR. 
     $      (iaNewGasID(iCount) .GT. kMaxGas)) THEN 
          write(kStdErr,*)'need a valid gas ID in *SPECTRA!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        IF ((iaNewData(iCount) .LT. 1) .OR. 
     $      (iaNewData(iCount) .GT. kNumkCompT)) THEN 
          write(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        END DO
        
c really no more error checking possible, as we have to actually go ahead and
c open the files as necessary

      RETURN
      END

c************************************************************************
c this subroutine is for nonlte
c this reads in the info for the new spectra for ONE gas
      SUBROUTINE nonlte(rLTEstrength,iNumLTEGases,
     $             iaLTEGasID,iaLTEData,iaaLTEChunks,caaaLTEChunks)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rLTEstrength   tells how strongly to add on all the files (default 1.0)
c iNumLTEGases   tells number of new gases 
c iaLTEGasID     tells which gases we want to update spectroscopy 
c iaLTEData      tells how many new data sets to read in for each gas 
c iaaLTEChunks   tells which data chunks to read in 
c caaaLTEChunks  tells the name of the files associated with the chunks 
      REAL rLTEstrength
      INTEGER iaLTEGasID(kGasStore),iaLTEData(kGasStore)
      INTEGER iNumLTEGases,iaaLTEChunks(kGasStore,kNumkCompT) 
      CHARACTER*80 caaaLTEChunks(kGasStore,kNumkCompT) 

c local variables
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iCount,iErr

      caWord='*NONLTE'
      iErr=-1

 5030 FORMAT(A130) 

      iNumLinesRead=0 
 13   IF (iNumLinesRead .GT. 0) THEN 
        iErr=1 
        WRITE(kStdErr,5010) caWord 
        CALL DoSTOP 
        END IF 
 5010 FORMAT('Error reading section ',A7) 

      IF (abs(rLTEstrength) .GT. 10.0) THEN
        write(kStdWarn,*)'wow! you multiply your nonLTE data by',rLTEstrength
        END IF

      iNumLinesRead=1
c read how many gases will have new nonLTE spectroscopy
      iErr=-1  
      IF ((iNumLTEGases .LT. 1) .OR. (iNumLTEGases .GT. kGasStore)) THEN  
        write(kStdErr,*)'need a valid number of gases in *NONLTE!!'  
        write(kStdErr,*)'must be > 0, < kGasStore!!'
        write(kStdErr,*)'please check and retry!' 
        CALL DoSTOP  
        END IF

      DO iCount=1,iNumLTEGases
        IF ((iaLTEGasID(iCount) .LT. 1) .OR. 
     $      (iaLTEGasID(iCount) .GT. kMaxGas)) THEN 
          write(kStdErr,*)'need a valid gas ID in *SPECTRA!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        IF ((iaLTEData(iCount) .LT. 1) .OR. 
     $      (iaLTEData(iCount) .GT. kNumkCompT)) THEN 
          write(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        END DO
        
c really no more error checking possible, as we have to actually go ahead and
c open the files as necessary

      RETURN
      END

c************************************************************************
