c Copyright 2000
c University of Maryland Baltimore County 
c All Rights Reserved

c this file deals with reading in the user file
c main routines for PTHFIL are here 
c also, the REFERENCE profile reading subroutine is here

c************************************************************************
c this subroutine deals with the 'PTHFIL' keyword 
      SUBROUTINE pthfil4(raaAmt,raaTemp,raaPress,raaPartPress,caPFName,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,iNumLayers)

      include 'kcarta.param'

c caPFName = character*80 profile name
c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gas ID's should be read in
c iaWhichGasRead = array that tracks whch gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iNumLayers= number of layers per gas profile (=kProfLayer)
c raLayerHeight = heights of layers in km
      REAL raLayerHeight(kProfLayer)
      INTEGER iNumLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      INTEGER iNpath
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

c local variables
      CHARACTER*7 caWord
      INTEGER iNumLinesRead

      caWord='*PTHFIL'

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('xError reading section ',A7,' of main user/PTHFIL file')

      write(kStdWarn,*)'Profile file to be read is  : '
      write(kStdWarn,*)caPfname
      CALL readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,iNumLayers,caPfName)

 5030 FORMAT(A130)

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword
c this differs from GENLN2 format in that
c (1) instead of veloctiy, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

      SUBROUTINE readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,iNumLayers,caPfName)

      include 'kcarta.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iNumLayers= number of layers per gas profile (=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
      INTEGER iNumLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

c local variables
      REAL raaHeight(kProfLayer,kGasStore)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iIOUN2,iNpath,iaNpathcounter(kProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn=0

c this variable keeps track of how many gases should be read in
      iNeed2Read=iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi
        IF (iaGases(iIDGas) .EQ. 1) THEN
          iNeed2Read=iNeed2Read-1
          END IF
        END DO

c this keeps track of the GasID used for the temperature .. hopefully water
      iTempFound=-1
c this keeps track of if we need to read in reference profiles
      iNeedMoreProfiles=-1

      caWord='*PTHFIL'
      iErr=-1

      iNumberOfGasesRead=0
c set all individual gas paths to zero        
      DO iNpath=1,kProfLayer
        iaNpathcounter(iNpath)=0
        END DO

c set this temp varaiable
      DO iNpath=1,kMaxGas
        iaAlreadyIn(iNpath)=-1
        END DO

c set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
c so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
c iaGases(1)=iaGases(2)=iaGases(22)=iaGases(51)=1
c so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
      DO iNpath=1,kMaxGas
        iaInputOrder(iNpath)=-1
        END DO
      iErr=1
      DO iNpath=1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr)=iNpath
          iErr=iErr+1
          END IF
        END DO

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        print *,iNumLinesRead
        print *,'profile name = ',caPfname        
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('yError reading section ',A7,' of main user/PTHFIL file')

      iIOUN2=kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
 1070     FORMAT('ERROR! number ',I5,' opening PRFILE path profile file'
     $            ,/,A80)
          CALL DoSTOP
        ENDIF
      kProfileUnitOpen=1

 30   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr
      CALL rightpad130(caStr)
      IF (caStr(1:1) .EQ. '!') THEN
        GO TO 30
      ELSE
        iNumLinesRead=iNumLinesRead+1
        READ (caStr,*,ERR=13,END=13) iNpath
        WRITE(kStdWarn,3000) iNpath,iNeed2Read*kProfLayer
 3000   FORMAT('input file has ',I5,' paths; kCARTA needs ',I5,' paths')
        iGasesInProfile=INT(iNpath*1.0/(kProfLayer*1.0))
        END IF

      IF (iNeed2Read*kProfLayer .GT. iNpath) THEN
        write(kStdWarn,*)'kCARTA/kLAYERS mismatch!!!!'
        write(kStdWarn,*)'Not enough gas paths in your supplied profile'
        iNeedMoreProfiles = 1
        CALL DoStop
        END IF

c now loop iNpath/iNumGases  times for each gas in the user supplied profile
 35   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr
      CALL rightpad130(caStr)
      IF (caStr(1:1) .EQ. '!') THEN
        GO TO 35
      ELSE
         iNumLinesRead=iNumLinesRead+1
         READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
         iaNpathCounter(iIDgas)=iaNpathCounter(iIDgas)+1

         IF (rAmt .lt. 0.0) THEN
           WRITE(kStdErr,1080)
           WRITE(kStdErr,1111) iIDgas,iaNpathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
           END IF

         IF (rT .lt. 0.0) THEN
           WRITE(kStdErr,1081)
           WRITE(kStdErr,1111) iIDgas,iaNpathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
           END IF

         IF (rP .lt. 0.0) THEN
           WRITE(kStdErr,1082)
           WRITE(kStdErr,1111) iIDgas,iaNpathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
           END IF

         IF (rPP .lt. 0.0) THEN
           WRITE(kStdErr,1083)
           WRITE(kStdErr,1111) iIDgas,iaNpathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
           END IF

         END IF

 1111      FORMAT('gasID, layer = ',I5,I5)
 1080      FORMAT('negative gas amount in PRFILE profile file')
 1081      FORMAT('negative gas temp in PRFILE profile file')
 1082      FORMAT('negative layer pressure in PRFILE profile file')
 1083      FORMAT('negative gas partial press in PRFILE profile file')


 40   CONTINUE
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
      IF (iaGases(iIDgas) .GT. 0) THEN
c always reset these variables at the beginning of this IF loop
         iFound=-1
         iGasIndex=1
 999     CONTINUE
         IF (iaInputOrder(iGasIndex) .EQ. iIDgas) THEN
           iFound=1
           END IF
         IF ((iFound .LT. 0) .AND. (iGasIndex .LT. iNumGases)) THEN
           iGasIndex=iGasIndex+1
           GO TO 999
           END IF
     		         
         IF (iFound .GT. 0) THEN 
           raaAmt(iaNpathCounter(iIDgas),iGasIndex)=rAmt
           raaTemp(iaNpathCounter(iIDgas),iGasIndex)=rT
           raaPress(iaNpathCounter(iIDgas),iGasIndex)=rP
           raaPartPress(iaNpathCounter(iIDgas),iGasIndex)=rPP
           raaHeight(iaNpathCounter(iIDgas),iGasIndex)=rH

           iaWhichGasRead(iIDgas)=1
           iaCont(iIDgas)=1       !continuum "always" included
           !water is a special case
           IF ((iIDGas .EQ. 1) .AND. (kCKD .GE. 0)) THEN
             iaCont(iIDgas)=1
           ELSE IF ((iIDGas .EQ. 1) .AND. (kCKD .LT. 0)) THEN
             iaCont(iIDgas)=-1
             END IF
           END IF
        END IF

c check to see if for the current gas (iGasID) we have read iNpath layers
      IF (iaNpathCounter(iIDgas) .LT. kProfLayer) THEN
        GOTO 35
        END IF

      WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
 4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3) 
      iFileGasesReadIn=iFileGasesReadIn+1

c this checks to see if we have read the profiles for all iNumGases required
c note that the gases read in MUST have been entered in GASFIL or XSCFIL 
c to count toward the tally ...
      IF (iaGases(iIDgas) .GT. 0) THEN
        iNumberOfGasesRead=iNumberOfGasesRead+1
        !so that this gas is not "reread" in from ReadOtherGases
        iaAlreadyIn(iNumberOfGasesRead) = iIDGas      
      ELSE
        write(kStdWarn,6000) iIDgas
 6000   FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
        END IF

      IF (iFileGasesReadIn .LT. iGasesInProfile) THEN
        GOTO 35
        END IF

      CLOSE(iIOUN2)
      kProfileUnitOpen=-1

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      DO iIDGas = kNewGasLo,kNewGasHi
        IF ((iaGases(iIDGas) .EQ. 1) .AND. (iaGases(1) .EQ. 1)) THEN
          write(kStdWarn,*)'Using water profile for gasID ',iIDGas
          iNumberOfGasesRead=iNumberOfGasesRead+1
          iaWhichGasRead(iIDgas)=1
          iFound=-1
          iGasIndex=1
 777      CONTINUE
          IF (iaInputOrder(iGasIndex) .EQ. iIDgas) THEN
            iFound=1
            END IF
          IF ((iFound .LT. 0) .AND. (iGasIndex .LT. iNumGases)) THEN
            iGasIndex=iGasIndex+1
            GO TO 777
            END IF
          !gasID=1 (water) has to be the first gas stuck in there!!!
          DO iP=1,kProfLayer
            raaAmt(iP,iGasIndex)       = raaAmt(iP,1)
            raaTemp(iP,iGasIndex)      = raaTemp(iP,1)
            raaPress(iP,iGasIndex)     = raaPress(iP,1)
            raaPartPress(iP,iGasIndex) = raaPartPress(iP,1)
            raaHeight(iP,iGasIndex)    = raaHeight(iP,1)
            END DO
        ELSEIF ((iaGases(iIDGas) .EQ. 1) .AND. (iaGases(1) .LT. 1)) THEN
          write(kStdErr,*) 'Cannot have continuum gas (101,102) w/o water'
          write(kStdErr,*) 'If you need to turn off water, but have continuum'
          write(kStdErr,*) 'you need to use the mixing table, not MOLGAS'
          CALL DoStop
          END IF
        END DO

c first check to see if all required gases found in the user supplied profile
      IF (iNumberOfGasesRead .LT. iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdErr,*) 'your profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        CALL DoStop
        END IF

 5030 FORMAT(A130)

c now set raLayerHeight
      DO iFound=1,kProfLayer
        raLayerHeight(iFound)=raaHeight(iFound,1)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine reads the profile files (for the references)
c it flags an error if kProfLayers layers are not read in
c ProX === A=amount,T=temperature,P=pressure,PP=partial pressure
       SUBROUTINE ReadRefProf(caFname,iIOUN,iNlay,raProA,raProT,
     $      raProP,raProPP,raQ,iL,iU,iErrOut)

      include 'kcarta.param'

c caFName   = name of file that has the profile
c iIOUN     = file IO UNIT number
c iNlay     = number of layers that are read in
c raProA/T/P/PP = profile amout, temperature,pressure,partial pressure
c raQ           = ratio of KLAYERS/AIRS reference gas amounts
c iL,iU     = lower/upper layer numbers (=1,kProfLayer)
c iErrOut   = error count (usually associated with file I/O)
       CHARACTER*80 caFname
       INTEGER iIOUN,iNlay,iL,iU,iErrOut
       REAL raProA(kProfLayer),raProT(kProfLayer),raQ(kProfLayer)
       REAL raProP(kProfLayer),raProPP(kProfLayer)

c local variables
       INTEGER iErr,iJ(kProfLayer),iX
       CHARACTER*100 caLine

       OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErr)
       IF (iErr .NE. 0) THEN
          WRITE(kStdErr,1080) iErr, caFname
 1080     FORMAT('ERROR! number ',I5,' opening profile file:',/,A82)
          CALL DoSTOP
       ENDIF
       kProfileUnitOpen=1
       
C      Read the file (skip comment lines)
       iNlay=0
 20    READ(iIOUN,5020,END=199) caLine
 5020  FORMAT(A100)
       IF (caLine(1:1) .NE. '!') THEN
          iNlay=iNlay+1
c old, before June 2011
c          READ(caLine,*) iJ(iNlay),raProP(iNlay),raProPP(iNlay),
c     $                   raProT(iNlay),raProA(iNlay),raQ(iNlay)

c new
          READ(caLine,*) iX,raProP(iNlay),raProPP(iNlay),
     $                   raProT(iNlay),raProA(iNlay)
          raQ(iNlay) = 1.0
          iJ(iNlay) = iNlay
       ENDIF
       GOTO 20
C      Close the file
 199   CLOSE(iIOUN)
       kProfileUnitOpen=-1

c check to see if the profile file includes layers iL to iU inclusive
       IF ((iJ(1) .gt. iL) .or. (iJ(iNlay) .lt. iU)) THEN
         iErrOut=1
         write(kStdErr,*)'Error while reading reference profile',caFname
         write(kStdErr,*)'not enough layers in the profile !! '
         write(kStdErr,*)'Low = ',iL,' ; bottom layer read in is ',iJ(1)
         write(kStdErr,*)'Top = ',iU,'; top layer read in is ',iJ(iNlay)
         CALL DoSTOP
         END IF

c check to see that ALL layers have been read in (could be greater than 100)
       IF (iNlay .NE. kProfLayer) THEN
         iErrOUt=1
         WRITE(kStdErr,500) caFName,iNlay
 500     FORMAT ('Profile File',/,A82,/,' has ',I2,' layers (need 100)')
         CALL DoSTOP
         END IF

       RETURN
       END

c************************************************************************
      SUBROUTINE rightpad130(caName)

      include 'kcarta.param' 

      CHARACTER*130 caName

      CHARACTER*130 caTempName
      INTEGER iR,iL,iInt

c find the "right" length of the input root name
      iR=len(caName)
 11      continue
      IF (caName(iR:iR) .eq. ' ') THEN
        iR=iR-1
        GO TO 11      
        END IF

c find the "left" length of the input root name
      iL=1
 12      continue
      IF (caName(iL:iL) .eq. ' ') THEN
        iL=iL+1
        GO TO 12      
        END IF
c thus the entered word exists between iL:iR

      IF (iL .GT. iR) THEN
        write(kStdErr,*) 'Fatal Error! Invalid (blank) string in file !'
        CALL DoSTOP 
        END IF

c now rearrange the string, so that it is right padded with blanks
c this is basically equivalent to  ADJUSTR(caName)
      DO iInt=1,130
        caTempName(iInt:iInt)=' '
        END DO
      caTempName(1:iR-iL+1)=caName(iL:iR)
      caName(1:130)=caTempName(1:130)

      RETURN
      END

c************************************************************************

c this subroutine deals with the 'WEIGHT' keyword == mixed path settings
c recall iNumGases    = total # of gases read in from MOLGAS and XSCGAS
c        iNpath       = iNumGases*iNumLayers
c        iNumLayers   = (iHigh-iLow+1) == kMaxLayer (max)
c        iMixFileLines= total number of uncommented lines in mixtable section

      SUBROUTINE mixfil4(raaMix,iNpmix,iNumGases,iNumLayers,
     $ iaGases,iNpath,caaMixFileLines,iMixFileLines)

      include 'kcarta.param'

c raaMix     = final mixing table
c iNpmix     = number of mixed paths that are read in = kProfLayer*Sets of MPs
c iNumGases  = number of gases from *MOLGAS + *XSCGAS
c iNumLayers = number of layers per gas profile (=kProfLayer)
c iaGases = allowed gas ID's
c iNpath     = total number of paths = iNumGases*kProfLayer
c caaMixFileLines = lines contining character description of mixtable
c iMixFileLines = number of lines in mixfile
      INTEGER iNpath,iNumLayers,iNpmix,iMixFileLines
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iNumGases,iaGases(kMaxGas)
      CHARACTER*130 caaMixFileLines(kProfLayer)

      CHARACTER*7 caWord

c local variables
      INTEGER iIpmix,iNpMixTemp,iGas

      caWord='*WEIGHT'

      iMixFileLines=-1

c initialize the mixing table to all 0.0's
      DO iGas=1,kGasStore
        DO iIpmix=1,kMixFilRows
          raaMix(iIpmix,iGas)=0.0
          END DO
        END DO

      IF (iNpmix .LT. 1) THEN
        WRITE(kStdErr,7000)
 7000   FORMAT('Error ! Need iNpmix > 0 (sets of mixed paths)')
        CALL DoSTOP
        END IF

c recall the data is for eack kProfLayer chunk ==> total number of mixed paths
c read in each time is really ......
      iNpMixTemp=iNpmix
      iNpMix=iNpMix*kProfLayer

      IF (iNpmix .GT. kMixFilRows) THEN
        WRITE(kStdErr,7010) iNpmix,kMixFilRows
 7010   FORMAT('iNpmix*kProfLayer = ',I4,' is greater than max allowed
     $ (from kcarta.param,',/,' kMixFilRows = ',I4,')')
        write(kStdErr,*)'Reset iNpmix (# of 100chunks in mixfile)'
        write(kStdErr,*)'or kMixFilRows, in file kcarta.param '
        CALL DoSTOP
        END IF

c now have to keep on parsing caaMixFileLines till we have the necessary info
c for the iNpMixTemp sets of mixed paths
      CALL ReadMixfil4(raaMix,iNpmixTemp,caaMixFileLines,
     $  iNumGases,iNumLayers,iMixFileLines,iaGases)

      RETURN
      END

c************************************************************************
c this subroutine reads in the specified number of free format reals
c from the input string ... if there are not enough, it will read in 
c the next line from the input file.

c NOTE THAT EACH MIXING TABLE IS READ IN AS A 
c     1 X NPATH TABLE,
c WHERE NPATH = iNumGases*iNumLayers. THE TABLE IS STORED AS AN
c      iNpmix X iNumGases TABLE,
c iNpmixTemp*kProfLayer = iNpmix == num of mixing paths that can be read in 
      SUBROUTINE ReadMixfil4(raaMix,iNpmixTemp,caaMixFileLines,
     $   iNumGases,iNumLayers,iMixFileLines,iaGases)

      include 'kcarta.param'

c caStr130     = current line to be processed. If not enough info found in it,
c                additional lines will be read in as necessary
c iNpmixTemp   = total number of mixed paths to be read in / kProfLayer
c iNumGases    = number of gases read in from *MOLGAS + *XSCGAS
c iNumLayers   = number of layers === kProfLayers
c iMixFileLines= number of relevant lines in mixfile that contains the mix info
c caaMixFileLines = lines contining character description of mixtable
c iaGases = allowed gas ID's
      CHARACTER*130 caaMixFileLines(kProfLayer)
      INTEGER iNumLayers,iNpmixTemp,iMixFileLines
      INTEGER iNumGases,iaGases(kMaxGas)
      REAL raaMix(kMixFilRows,kGasStore)

      CHARACTER*130 caaLinesTemp(kProfLayer)
      REAL raStore(kMaxGas),raStore2(kMaxGas)
      INTEGER iInt,iErr, iLay,iGas,iNpath
      CHARACTER*130 caStr130
c iaInputOrder = gas ID's in order they were read in
      INTEGER iaInputOrder(kMaxGas),iaStore(kMaxGas)
      INTEGER iMPReadIn, iMpSet, iSpecial, iW
      REAL rC

c do initializations ...

      DO iNpath=1,kMaxGas
        iaInputOrder(iNpath)=-1
        END DO
      iErr=1
      DO iNpath=1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr)=iNpath
          iErr=iErr+1
          END IF
        END DO

      iNpath=iNumGases*iNumLayers
      iErr=-1

 99   IF (iErr .GT. 0) THEN
        write(kStdErr,*)'Error while reading in *WEIGHT ... '
        CALL DoSTOP
        END IF

 30   CONTINUE
      DO iInt=1,kMaxGas
        raStore(iInt)=0.0
        raStore2(iInt)=0.0
        END DO

c remember there are three possibilities when defining a mixed path set
c  (a) iN  list of weights       ... individually define weights
c  (b) iN  -1 rW -1              ... default weights for all gases to rW
c  (c) iN  -1 rW iG              ... default weights for all gases to rW
c        i1 r1 i2 r2 ..... iiG riG      except for a few gases

      iMPReadIn=0
      iMixFileLines=1
 50   CONTINUE
      caStr130=caaMixFileLines(iMixFileLines)
c     if user wants to separately list weight, then the input format is 
c          iNatm    list of weights (> 0)
c     if user wants to default all/most values, then the input format is 
c          iNatm    -1   rWeight -1
      read (caStr130,*) iMPSet,rC
      IF (rc .gt. 0.0) THEN         !we know user is listing gas weights
        !iaStore will be irrelevant
        !info will be in raStore2
        CALL readlist(caaMixFileLines,-1,iNumgases,iMixFileLines,
     $                raStore2,iaStore)
        DO iInt=1,iNumGases
          raStore(iaInputOrder(iInt))=raStore2(iInt)
          END DO
      ELSE
        read (caStr130,*) iMPSet,iW,rC,iSpecial
        IF (iSpecial .lt. 0) THEN            !all gases to have same weight
          iMixFileLines=iMixFileLines+1      !consider one more line read 
          DO iInt=1,kMaxGas
            raStore(iInt)=rC
            END DO
        ELSE
          !first default all gas weights to rC
          DO iInt=1,kMaxGas
            raStore(iInt)=rC
            END DO
          CALL readlist(caaMixFileLines,1,iSpecial,iMixFileLines,
     $                  raStore2,iaStore)
          !then add on the special weights
          DO iLay=1,iSpecial 
            IF (iaGases(iaStore(iLay)) .GT. 0) THEN 
              raStore(iaStore(iLay))=raStore2(iLay) 
            ELSE 
              write(kStdErr,*)'Error while reading in *WEIGHT ... ' 
              write(kStdErr,*)'Gas ID',(iaStore(iLay)),' not in  
     $ *MOLGAS or *XSCGAS' 
              CALL DoSTOP 
              END IF   
            END DO
          END IF
        END IF

c having successfully read in the current 100 mp block of the mixing table,  
c set raaStore from the elements stored in raStore 
      write(kStdWarn,*) 'read in info for weighted set number ',iMPReadIn+1
      write(kStdWarn,*) 'the weights are : ' 
      DO iGas=1,iNumGases 
        write(kStdWarn,*)iGas,iaInputOrder(iGas),raStore(iaInputOrder(iGas)) 
        END DO 
      write(kStdWarn,*) ' ' 
      DO iGas=1,iNumGases 
        DO iLay=1,iNumLayers 
          raaMix(iMPReadIn*iNumLayers+iLay,iGas)=
     $             raStore(iaInputOrder(iGas))
          END DO 
        END DO 

      iMPReadIn=iMPReadIN+1

      IF (iMPReadIn .lt. iNpMixTemp) THEN       !read in next set of MPs
        GOTO 50
        END IF

      iMixFileLines=iMixFileLines-1   !this is correct number of caaMix read
                               !as code is ready to read in next unread line


c but we have to add on the iNpmix as the first line
      DO iGas=1,iMixFileLines
        caaLinesTemp(iGas+1)=caaMixFileLines(iGas)
        END DO

c with iNpmix added on, this is the number that will get dumped to output file
      iMixFileLines=iMixFileLines+1   

      DO iGas=2,iMixFileLines
        caaMixFileLines(iGas)=caaLinesTemp(iGas)
        END DO
      write(caStr130,37) iNpmixTemp
      caaMixFileLines(1)=caStr130
 37   FORMAT(I5)

      RETURN
      END

c************************************************************************
c this subroutine reads stuff from caaM, and stores the resuls either in
c raStore alone or iaStore,raStore depending on value of iType
      SUBROUTINE readlist(caaM,iType,iNum,iLines,raStore,iaStore)

      include 'kcarta.param'

c INPUT PARAMS
c caaM   = input array of character strings to be parsed
c iType  = -1 if want to read in one integer and bunch of reals
c        = +1 if want to read in integer,real combos
c iNum   = how many reals, or how many (int, real) combos to read in
c iLines = where we are in array index before we start looking for info
c OUTPUT/MODIFIED PARAMS
c iLines = where we finally are in array index before all info found
c raStore = output list of reals that are read in
c iaStore = output list of integers that are read in
      CHARACTER*130 caaM(kProfLayer)
      INTEGER iaStore(kMaxGas),iType,iNum,iLines
      REAL raStore(kMaxGas)

      INTEGER iI,iJ,iAdd
      REAL raStore2(kMaxGas)
      
      DO iI=1,kMaxGas
        iaStore(iI)=-1
        raStore(iI)=-1.0
        raStore2(iI)=-1.0
        END DO

      IF (iType .LT. 0) THEN
        !read in list of integer (discard) followed by iNum reals, using 
        !current line as the line to start reading stuff from
        iLines=iLines
      ELSEIF (iType .GT. 0) THEN
        !read in list of (integer,real) pairs,  using next line as the 
        !line to start reading stuff from
        iLines=iLines+1
        END IF

      iAdd=0
 5    IF (iAdd .GT. 0) THEN 
        write (kStdErr,*) 'whoops! wierd error parsing in caaMixFileInfo'
        CALL DoSTOP
        END IF
 10   IF (iAdd .GT. 0) THEN
         CLOSE(unit=kTempUnit)
         END IF
      iAdd=iAdd+1
      OPEN(unit=kTempUnit,file=kDumbFile)
      DO iJ=0,iAdd-1
        write (kTempUnit,*) caaM(iLines+iJ)
        END DO
      CLOSE(unit=kTempUnit)

      IF (iType. LT. 0) THEN
        !read in list of integer followed by iNum reals
        OPEN(unit=kTempUnit,file=kDumbFile)
        read (kTempUnit,*,end=10,err=5) iJ,(raStore(iI),iI=1,iNum)
        CLOSE(unit=kTempUnit)
        iLines=iLines+iAdd
        END IF

      IF (iType. GT. 0) THEN
        !read in list of iNum (integer,real) pairs
        OPEN(unit=kTempUnit,file=kDumbFile)
        read (kTempUnit,*,end=10,err=5) (raStore2(iI),iI=1,2*iNum)
        CLOSE(unit=kTempUnit)
        DO iI=1,iNum
          iaStore(iI)=int(raStore2(iI*2-1))
          raStore(iI)=    raStore2(iI*2)
          END DO
        iLines=iLines+iAdd
        END IF

      RETURN
      END

c************************************************************************
