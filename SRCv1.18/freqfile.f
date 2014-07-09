c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c******** THIS FILE CONTAINS ROUTINES THAT FIND FREQ CHUNKS ***** *******
c************************************************************************
c this function checks to see if the cross section data base exists for
c gasID, between freqs rL,rH
c if rL,rH < 0 then checking basic presence of gas in data base
      INTEGER FUNCTION iCheckXsecDataBase(iGasID,rL,rH,iTagIn,iErr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c iGasID = GAS ID to be searched for in the database
c rL,rH  = freq start/stop points
c iErr   = error status (mainly associated with not finding the relevant file)
c iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
c          (not used if kXsecFormat < 0)
      INTEGER iGasID,iErr,iTagIn
      REAL rL,rH

c local variables
      INTEGER iIOUN,iFileErr,iID,iTag
      INTEGER iLine,iNpts,iTemps
      CHARACTER*80 caLine,caFName
      REAL rLower,rHigher

c assume GASID , freqs are wrong
      iCheckXsecDataBase=-1

      caFName = kXsecParamFile
      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=caFName,STATUS='old',
     $    FORM='FORMATTED',IOSTAT=iFileErr)

      IF (kXsecFormat .LT. 0) THEN      
ccccccccccc this is the original format : read old style xsec.param file
        IF (iFileErr .NE. 0) THEN
          iErr=0
          WRITE(kStdErr,103) iFileErr,caFName
 103      FORMAT('ERROR! number ',I5,' opening xsec database file : 
     $    ',/,A84)
          CALL DoSTOP
        END IF
        kTempUnitOpen=1

c read file util GASID, freq bounds match found or EOF
 20     READ(iIOUN,5020,END=777) caLine
        READ(caLine,*) iLine,iID,iNpts,iTemps,rLower,rHigher

        IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is the -1 option in XSCGAS
          iCheckXsecDataBase=1
        END IF

        IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $    (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
          iCheckXsecDataBase=1
        END IF

        IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
     $     (rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
          iCheckXsecDataBase=1
        END IF

        IF (iCheckXsecDataBase .LT. 0) THEN
          GOTO 20
        END IF
      
 777    CONTINUE
      ELSE
ccccccccccc this is the new format : read comp.param style xsec.param file
c read file util GASID, freq bounds match found or EOF
 30     READ(iIOUN,5020,END=888) caLine
        READ(caLine,*) iID,rLower,rHigher,iTag

        IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is basically the -1 option in XSCGAS
          iCheckXSecDataBase=1
        END IF

        IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $            (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
c this is when we WANT to UNCOMPRESS files!
          IF (iTag .EQ. iTagIn) THEN   !actually uncompressing stuff
            iCheckXsecDataBase=1
          END IF
        END IF

        IF (iCheckXSecDataBase .LT. 0) THEN
          GOTO 30
        END IF
 888    CONTINUE

      END IF

      CLOSE(iIOUN)
      kTempUnitOpen=-1

 5020 FORMAT(A80)
      RETURN
      END
c************************************************************************
c this function checks to see if the comp data file exists for
c gasID, between freqs rL,rH
c if rL,rH < 0 then checking basic presence of gas in data base
      INTEGER FUNCTION iCheckCompDataBase(iGasID,rL,rH,iTagIn,iErr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID = GAS ID to be searched for in the database
c rL,rH  = freq start/stop points
c iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
c iErr   = error status (mainly associated with not finding the relevant file)
      INTEGER iGasID,iErr,iTagIn
      REAL rL,rH

c local variables
      INTEGER iIOUN,iFileErr,iID,iTag
      REAL rLower,rHigher
      CHARACTER*80 caLine,caFname

c assume GASID , freqs are wrong
      iCheckCompDataBase=-1

      caFName = kCompParamFile
      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='old',
     $    FORM='FORMATTED',IOSTAT=iFileErr)

      IF (iFileErr .NE. 0) THEN
        iErr=0
        WRITE(kStdErr,103) iFileErr,caFname
 103    FORMAT('ERROR! number ',I5,' opening comp database file : 
     $  ',/,A84)
        CALL DoSTOP
      END IF
      kTempUnitOpen=1

c read file util GASID, freq bounds match found or EOF
 20   READ(iIOUN,5020,END=777) caLine
 5020 FORMAT(A80)
      READ(caLine,*) iID,rLower,rHigher,iTag

      IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is basically the -1 option in MOLGAS
        iCheckCompDataBase=1
      END IF

      IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $            (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
c this is when we WANT to UNCOMPRESS files!
        IF (iTag .EQ. iTagIn) THEN   !actually uncompressing stuff
          iCheckCompDataBase=1
        END IF
      END IF

      IF (iCheckCompDataBase .LT. 0) THEN
        GOTO 20
      END IF
      
 777  CONTINUE
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      RETURN
      END
c************************************************************************
c this function checks to see if the GAS ID/frequency combination are in
c the compressed data base or in the xsec database
      SUBROUTINE DataBaseCheck(iGasID,raFreq,iTag,iActualTag,iDoAdd,iErr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = GAS ID to be searched for in CompDataBase or XsecDataBase
c raFreq   = wavenumber array that contains present chunk of 25 cm-1
c iErr     = errors (associated with file I/O)
c iTag     = 1,2,3 telling us wavenumber spacing 0.001 0.0025 0.005 cm-1
c iDoAdd   = -1,+1 tells us if we add on current gas to current 10000 pts
      INTEGER iGasID,iErr,iTag,iActualTag,iDoAdd
      REAL raFreq(kMaxPts)

      INTEGER iCheckCompDataBase,iCheckXsecDataBase

      iDoAdd=-1

 333  FORMAT('---> Warning! No comprsd data gasID ',I3,
     $ ' in ',f10.4,' cm-1 chunk; tagIN = ', I3)
 334  FORMAT('---> Warning! No comprsd data heavywater gasID ',I3,
     $ ' in ',f10.4,' cm-1 chunk')
 414  FORMAT('Comp data exists for HeavyWater (gas#) ',I3,f10.4,' cm-1 chunk')

c check to see if the k-comp file exists for the (GAS ID/freq) combination
      IF ((1 .LE. iGasID) .AND. (iGasID .LE. kGasComp)) THEN
        iDoAdd=
     $   iCheckCompDataBase(iGasID,raFreq(1),raFreq(kMaxPts),iActualTag,iErr)
        IF (iDoAdd .LT. 0) THEN
          WRITE(kStdWarn,333) iGasID,raFreq(1),iActualTag
          WRITE(kStdWarn,*) ' '
        END IF
        GOTO 2000
      END IF

      IF (iGasID .EQ. kMaxGas) THEN 
        !! this is heavy water
        IF (((raFreq(1) + 1.0) .GE. kWaterIsobandStart1*1.0)  .AND. 
     $      ((raFreq(kMaxPts) - 1.0) .LE. kWaterIsobandStop1*1.0)) THEN
          iDoAdd = +1
          WRITE(kStdWarn,414) iGasID,raFreq(1)
        ELSEIF (((raFreq(1) + 1.0) .GE. kWaterIsobandStart2*1.0)  .AND. 
     $      ((raFreq(kMaxPts) - 1.0) .LE. kWaterIsobandStop2*1.0)) THEN
          iDoAdd = +1
          WRITE(kStdWarn,414) iGasID,raFreq(1)
        ELSE
          iDoAdd = -1
          WRITE(kStdWarn,334) iGasID,raFreq(1)
        END IF
      END IF

c check to see if the xsec file exists for the (GAS ID/freq) combination
 444  FORMAT('---> Warning! No xsec data for gasID ',I3,
     $ ' in ',f10.4,' cm-1 chunk')

      IF ((kGasXsecLo .LE. iGasID) .AND. (iGasID .LE. kGasXsecHi)) THEN
        iDoAdd=iCheckXsecDataBase(iGasID,raFreq(1),raFreq(kMaxPts),
     $                            iActualTag,iErr)
        IF (iDoAdd .LT. 0) THEN
          WRITE(kStdWarn,444) iGasID,raFreq(1)
          WRITE(kStdWarn,*) ' '
        END IF
        GOTO 2000
      END IF

ccc   IF ((29 .LE. iGasID) .AND. (iGasID .LE. 50)) THEN
      IF ((kGasComp+1 .LE. iGasID) .AND. (iGasID .LE. kGasXsecLo-1)) THEN
        iDoAdd=-1
        WRITE(kStdWarn,1000) iGasID
 1000   FORMAT('No contribution to absorption cross sections from GAS ID
     $ = ',I2,/,'... skipping to next gas ...')
        GOTO 2000
      END IF

c check to see if we need to add on the water continuum
      IF (kCKD .ge. 0) THEN
        IF ((iGasID .GE. kNewGasLo) .AND. (iGasID .LE. kNewGasHi)) THEN
          iDoAdd = 1
          GOTO 2000
        END IF
      END IF

      IF ((iGasID .LT. 1) .OR. (iGasID .GT. kMaxGas)) THEN
        iErr=1
        WRITE(kStdWarn,1010) iGasID
 1010   FORMAT('Cannot add contribution of GAS ID = ',I2)
        CALL DoSTOP
      END IF

 2000 CONTINUE
      RETURN
      END 

c************************************************************************
c this subroutine checks that kNumKcompT, as set in kcarta.param
c agrees with what one would expect from kaMinFr,kaMaxFr

      SUBROUTINE Check_kaNum_ka100layerCloudType

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

      INTEGER iI,iJ
      REAL rF,rG

      !! in pre_defined we use ka100layerCloudType(kMaxUserSet)
      !! but then in kcarta.param we define kMaxClouds
      !! beter have kMaxUserSet > kMaxClouds else array bound problems!!!
      IF (kMaxUserSet .LT. kMaxClouds) THEN
        write(kStdErr,*) 'because of comblock6, need kMaxUserSet > kMaxClouds'
        write(kStdErr,*) 'kMaxUserSet,kMaxClouds = ',kMaxUserSet,kMaxClouds
        CALL DoStop
      END IF

      DO iJ=1,kW
        rF = kMaxPts*kaFrStep(iJ)
        rG = kaBlSize(iJ)
        IF (abs(rF-kaBlSize(iJ)) .GT. kaFrStep(iJ)/2.0) THEN
          write(kStdErr,*) 'iJ= ',iJ
          write(kStdErr,*) 'kcarta.param claims kaBlSize(iJ)=',rG
          write(kStdErr,*) 'while it should be ',rF
          CALL DoSTOP
        END IF
        iI = NINT((kaMaxFr(iJ)-kaMinFr(iJ))/rF)
        IF (iI .NE. kaNumKComp(iJ)) THEN
          write(kStdErr,*) 'iJ = ',iJ
          write(kStdErr,*) 'kcarta.param says that the number of '
          write(kStdErr,*) 'kCompressed files = kaNumKComp(iJ) = ',kaNumKComp(iJ)
          write(kStdErr,*) 'based on following parameters, there should be iI = ',iI
          write(kStdErr,*) kaMinFr(iJ),kaMaxFr(iJ),kaFrStep(iJ),kMaxPts
          write(kStdErr,*) 'iI = iNT((kMaxFreq-kMinFreq)/
     $(kMaxPts*kFreqStep))'
          CALL DoSTOP
        END IF
      END DO

      iI=0
      DO iJ=1,kW
        iI = iI+(kaNumkComp(iJ))
      END DO
      IF (iI .NE. kNumkCompT) THEN
        write(kStdErr,*) 'kcarta.param says kNumkCompT = ',kNumkCompT
        write(kStdErr,*) 'while it should be ',iI
        CALL DoSTOP
      END IF

      RETURN
      END
       
c************************************************************************
c this subroutine uses the input file to set the file  number for the
c compressed data
c with the following three constraints : 
c   1) the wavenumbers have to be between rFreqMin1 cm-1 and rFreqMax3 cm-1.
c   2) since each k-compressed file is XYZ cm-1 long, (iU-iL) <= XYZ
c   3) iL and iU have to fall in the same XYZ cm-1 block
c If all checks out, the relevant file block is set up, as is the
c wavenumber array 
      SUBROUTINE GetFreq(rFrLow,rFrHigh,
     $                   iFileIDLo,iFileIDHi,
     $                   raBlock,raFiles,
     $                   iaTagIndex,iaActualTag,
     $                   raFileStep,iaList,iTotal) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFrLow    = lower frequency bound      
c rFrHigh   = upper frequency bound      
c iFileIDLo    = lower file index ID (between 1--kNumKcompT)
c iFileIDLo    = upper file index ID (between 1--kNumKcompT)
c iaTagIndex tells which TagIndex (1 .. 10) is associated with which file 
c     (1,2,3,4 for k,p,q,r etc
c iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
c   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
c   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
c raFiles  tells which is the current kComp "file number" eg 630, 805 etc 
c          very useful so program easily knows r630_g1.dat etc 
c raBlock  tells the current kComp wavenumber block 
c raFileStep tells you the current wavenumber step size*10000 
c iaList   has the final list of iTotal files that should be uncompressed
      INTEGER iaList(kNumkCompT),iTotal
      REAL raFileStep(kNumkCompT),raFiles(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT) 
      REAL raBlock(kNumkCompT)
      REAL rFrLow,rFrHigh
      INTEGER iFileIDLo,iFileIDHi

c local variables
c rFileStartFrLo = lower file ID (605 -- 2805)
c rFileStartFrHi = upper file ID 
      REAL rFileStartFrLo,rFileStartFrHi
      REAL rLow,rHigh,rTemp
      INTEGER iDummy

      write(kStdWarn,*) '*********************************************'
      write(kStdWarn,*) '**Setting freqs, fileID lower/upper bounds***'

c      print *,'input (rFrLow,rFrHigh) : '
c      read *,rFrLow,rFrHigh

      IF (rFrLow .ge. rFrHigh) THEN
        write(kStdWarn,*) 'Swapping frequencies found in *FRQNCY'
        rTemp = rFrLow
        rFrLow = rFrHigh
        rFrHigh = rTemp
      END IF        
      write(kStdWarn,*) 'input freqs are  : ',rFrLow,rFrHigh

c remember : int(x)  === trunc(x)
c remember : nint(x) === round(x)
c if the bands are FIR3,FIR2,FIR1,IR,NIR1,NIR2 onwards then
      IF (rFrlow .GE. 140.0) THEN
        !! start/stop chunks are integer numbers, so use this fact
        rLow = rFrLow
        rHigh = rFrHigh 
        !first round down rFrLow, and round up rFrHigh
        rFrLow=1.0*INT(rFrLow)
        IF (abs(rFrHigh-1.0*INT(rFrHigh)) .GT. 0.000000) THEN
          rFrHigh = rFrHigh+1.0
        END IF
        rFrHigh=1.0*INT(rFrHigh)
        write(kStdWarn,*) 'rounded down/up rFrLow/rFrhigh to ',rFrLow,rFrHigh
      ELSE
        write(kStdWarn,*) ' --> for bands with wavenumbers < 140 cm-1, kCARTA'
        write(kStdWarn,*) ' --> compressed files are not on integer bdries'
        write(kStdWarn,*) ' --> so we hope you got things correct!'
      END IF

      rLow = rFrLow
      rHigh = rFrHigh

      IF (rLow .ge. rHigh) THEN
        rTemp = rLow
        rLow = rHigh
        rHigh = rTemp
      END IF        

      CALL filebounds(rLow,rHigh,
     $             rFileStartFrLo,rFileStartFrHi,iFileIDLo,iFileIDHi,
     $             raBlock,raFiles,
     $             iaTagIndex,iaActualTag,
     $             raFileStep,iaList,iTotal)     

      rFrLow = rLow
      rFrHigh = rHigh
      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'rFrLow,rFrHigh after = ',rFrLow,rFrHigh
      write(kStdWarn,*) 'num chunks = ',iTotal
      write(kStdWarn,*) 'Start/Stop FileTags = ',
     $   iaActualTag(iaList(1)),iaActualTag(iaList(iTotal))
      write(kStdWarn,*) 'compfile 1: freq,FileID ',rFileStartFrLo,iFileIDLo
      write(kStdWarn,*) 'compfile N: freq,FileID ',rFileStartFrHi,iFileIDHi

      write (kStdWarn,*) 'File info : '
      write (kStdWarn,*) '#,rBlock,iFile,iTagIndex,iActualTag,iFileStep,iList'
      write (kStdWarn,*) '----------------------------------------------'
      do iDummy = iFileIDLo,iFileIDHi
        write(kStdWarn,*) iDummy,raBlock(iDummy),raFiles(iDummy),
     $    iaTagIndex(iDummy),iaActualTag(iDummy),
     $    raFileStep(iDummy),iaList(iDummy-iFileIDLo+1)
      END DO

      IF (iaActualTag(iaList(1)) .NE. iaActualTag(iaList(iTotal))) THEN
        write(kStdWarn,*) 'Start file tag = ',iaActualTag(iaList(1))
        write(kStdWarn,*) 'Stop  file tag = ',iaActualTag(iaList(iTotal))

        write(kStdErr,*) 'program requires you choose start/stop freqs'
        write(kStdErr,*) 'that only span one wavenumber spacing, as '
        write(kStdErr,*) 'defined at the bottom of kcarta.param'
        DO iDummy = 1,kW
          write(kStdErr,333) kaMinFr(iDummy),kaMaxFr(iDummy),
     $                 kaFrStep(iDummy),kaTag(iDummy)
        END DO
        write(kStdErr,*) 'please reset *FRQNCY and retry',rFrLow,rFrHigh
        CALL DoSTOP
      END IF

 333  FORMAT('rF1,rF2,d_f,Tag = ',2(f8.2,'  '),f12.7,' ',i4)

      write(kStdWarn,*)'****Finished freqs, fileID lower/upper bounds**'

      RETURN
      END

c************************************************************************
c this subroutine calculates which of the relevant k-comp
c blocks are required for the problem ... this is for the user input files
c also makes sure the fres lie between kMinFreq1,kMaxFreq3 ... else it sets 
c them to default values and asks if the user wishes to continue
c artificially set up the "imaginary" kNumkComp th point to account for 2805

c note we have to be smart about possibility of overlapping q,r,s eg
      SUBROUTINE filebounds(rL,rH,
     $             rFileStartFrLo,rFileStartFrHi,iFileIDLo,iFileIDHi,
     $             raBlock,raFiles,
     $             iaTagIndex,iaActualTag,raFileStep,iaList,iTotal) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rL,rH     = from namelist file, these are the lower and upper wavenumber
c             bounds : could be reset here
c rFileStartFrLo = from rL, lower file ID bound (in wavenumbers e.g. 605)
c rFileStartFrHi = from rH, upper file ID bound (in wavenumbers e.g. 2780)
c iFileIDLo    = from rL, lower fileID (e.g. 605 ==> 1+kNumKfile(1))
c iFileIDHi    = from rH, upper fileID(e.g. 2780 ==>kNumKcompT(1)+
c                                                            kNumKcompT(2))
c iaTagIndex tells which TagIndex (1 .. 10) is associated with which file 
c     (1,2,3,4 for k,p,q,r etc
c iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
c   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
c   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
c raFiles  tells which is the current kComp "file number" eg 630, 805 etc 
c          very useful so program easily knows r630_g1.dat etc 
c raFileStep tells you the current wavenumber step size*10000 
c raBlock  tells the current kComp wavenumber block (ie start freq)

c iaList   has the final list of iTotal files that should be uncompressed
      INTEGER iaList(kNumkCompT),iTotal
      REAL raFileStep(kNumkCompT),raFiles(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT) 
      REAL raBlock(kNumkCompT),rL,rH
      REAL rFileStartFrLo,rFileStartFrHi
      INTEGER iFileIDLo,iFileIDHi

c local variables
      REAL rFileStartFrLo1,rFileStartFrHi1
      REAL rFileStartFrLo2,rFileStartFrHi2
      INTEGER iFileIDLo1,iFileIDHi1,iFileIDLo2,iFileIDHi2
      REAL raBlockEnd(kNumkCompT),rL1,rH1,rL2,rH2,rEnd,raDelta(kNumkCompT)
      INTEGER iInt,iFound,iTruthLo,iTruthHi
      INTEGER iErr,iDummy,iJunk

c first check lower, upper limits
      CALL CheckLimits(rL,rH)

      iErr=0
      DO iDummy=1,kW            !compute total number of kCOMP chunks present
        iErr = iErr+kaNumkComp(iDummy)
      END DO
      IF (iErr .NE. kNumkCompT) THEN
        write(kStdErr,*) 'Error in kNumkCompT',iErr,kNumkCompT
        CALL DoStop
      END IF

c note there could be "overlaps" between the q r s files eg
c     < 1=500-550> <2=550-600>  <3=600-650> <4=650-700>
c                                  <5=605-630> <6=630-655> <7=655-680> .....
      iDummy=0 
      DO iFound=1,kW 
        DO iInt=1,kaNumkComp(iFound) 
          iDummy = iDummy+1 
          !this is needed by kcartamain.f
          raBlock(iDummy)     = kaMinFr(iFound)+(iInt-1)*kaBlSize(iFound) 
          raFiles(iDummy)     = raBlock(iDummy) 
          iaActualTag(iDummy) = kaTag(iFound) 
          iaTagIndex(iDummy)  = iFound
          raFileStep(iDummy)  = kaBlSize(iFound) 

          !this is needed within misc.f
          raBlockEnd(iDummy) = raBlock(iDummy)+kaBlSize(iFound)     
          raBlockEnd(iDummy) = raBlockEnd(iDummy)-kaFrStep(iFound)  
          raDelta(iDummy)    = kaFrStep(iFound)
        END DO 
      END DO 

c go thru raBlock and see where we think the start file ID, stop file ID
c should be set at
      CALL LowerLimits(rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $     rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,raFiles,raDelta)
      CALL UpperLimits(rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $     rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,raFiles,raDelta)

c now set the fileID's to be used
      CALL LowerFileID(rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $  rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,
     $  iaActualTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,rFileStartFrLo)
      CALL UpperFileID(rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $  rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,
     $  iaActualTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,rFileStartFrHi)

      IF ((iFileIDLo .GE. 1) .AND. (iFileIDLo .LE. kNumkCompT) .AND.
     $    (iFileIDHi .GE. 1) .AND. (iFileIDHi .LE. kNumkCompT) .AND.
     $    (iFileIDLo .LE. iFileIDHi)  .AND. 
     $    (iTruthLo .GT. 0) .AND. (iTruthHi .GT. 0)) THEN
        iErr=-1
      ELSE
        write(kStdErr,*) 'Error in setting file '
        write(kStdErr,*) 'lower/upper bounds from *FRQNCY'
        CALL DoSTOP
      END IF

      iJunk = kLongOrShort
      iJunk = 2
      IF  (abs(iJunk) .LE. 1) THEN  !! verbose printing
        iFound = iaActualTag(1)
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'iCnt  raBlock  raBlockEnd  iaTag'
        write(kStdWarn,*) '--------------------------------'
        DO iDummy = 1,kNumkCompT-1
          IF ((iDummy. NE. iFileIDLo) .AND. (iDummy .NE. iFileIDHi)) THEN
            write(kStdWarn,10) iDummy,raBlock(iDummy),raBlockEnd(iDummy),
     $       iaActualTag(iDummy)
          ELSE
            write(kStdWarn,11) iDummy,raBlock(iDummy),raBlockEnd(iDummy),
     $       iaActualTag(iDummy)
          END IF
          IF (iFound .ne.  iaActualTag(iDummy+1)) THEN
            write(kStdWarn,*) '--------------------------------'        
            iFound =  iaActualTag(iDummy+1)
          END IF
        END DO
        iDummy = kNumkCompT
        IF ((iDummy. NE. iFileIDLo) .AND. (iDummy .NE. iFileIDHi)) THEN
          write(kStdWarn,10) iDummy,raBlock(iDummy),raBlockEnd(iDummy),
     $       iaActualTag(iDummy)
        ELSE
          write(kStdWarn,11) iDummy,raBlock(iDummy),raBlockEnd(iDummy),
     $       iaActualTag(iDummy)
        END IF
        write(kStdWarn,*) '--------------------------------'
        write(kStdWarn,*) ' '
      ELSE  !! quiet printing
        iFound = iaActualTag(1)
        write(kStdWarn,*) ' '
        DO iDummy = 1,kNumkCompT-1
          IF (iFound .ne.  iaActualTag(iDummy+1)) THEN
            iFound =  iaActualTag(iDummy+1)
          END IF
        END DO
        iDummy = kNumkCompT
      END IF
        
 10   FORMAT(I4,' ',F12.5,' ',F12.5,' ',I3)
 11   FORMAT(I4,' ',F12.5,' ',F12.5,' ',I3,' <<<<<<<<<<<<<<<')

c now create iaList as necessary
c remember that iFileID is basically an index setting equal to iDummy above
c thus if we know iFileID we know everything!!!!!
      IF (iaActualTag(iFileIDLo) .EQ. iaActualTag(iFileIDHi)) THEN  
        !very easy everything is in either q or r or s database
        iTotal=0
        DO iInt=1,(iFileIDHi-iFileIDLo+1)
          iTotal = iTotal+1
          iaList(iTotal) = iFileIDLo+(iInt-1)     !save file ID
          IF  (abs(kLongOrShort) .LE. 1) THEN  !! verbose printing
            write(kStdWarn,*) iTotal,iaActualTag(iaList(iTotal)),
     $          raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
          END IF
        END DO
      ELSE
        !very hard : mixing of q,r,s databases
        iTotal=0
        iTotal = iTotal+1
        iaList(iTotal) = iFileIDLo                !save file ID
        rEnd = raBlockEnd(iaList(iTotal))
        write(kStdWarn,*) iTotal,iaActualTag(iaList(iTotal)),
     $              raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
        DO iInt=2,(iFileIDHi-iFileIDLo+1)
          IF (rEnd .LT. raBlockEnd(iFileIDLo+iInt-1)) THEN
            iTotal = iTotal+1
            iaList(iTotal) = iFileIDLo+(iInt-1)   !save file ID
            write(kStdWarn,*) iTotal,iaActualTag(iaList(iTotal)),
     $             raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
            rEnd = raBlockEnd(iaList(iTotal))
          END IF
        END DO
      END IF

c now that we have the ActualTag eg 20 for r=0605-2830 cm-1, or 
c                                   25 for s=2805-3305 cm-1, or 
c                                   30 for m=4000-5000 cm-1,
c we have also mapped this back into an TagIndex

c      call dostop

      RETURN
      END

c************************************************************************
c this subroutine checks to make sure rL,rH lie within the lower/upper bounds
c of wavenumbers in the kComp files
      SUBROUTINE CheckLimits(rL,rH)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rl,rH are the lower and upper limits the user sends in; this subroutine 
c                       can change them

      REAL rL,rH

      IF (rL .LT. kaMinFr(1)) THEN
        write(kStdWarn,*) 'Error!!Setting min wavenumber to ',kaMinFr(1)
        rL = kaMinFr(1)
      END IF
      IF (rL .GT. (kaMaxFr(kW)-kaBlSize(kW))) THEN
        !! input start wavenumber waaaaaay tooooooooo high
        write(kStdWarn,*) 'Error!!Setting minimum wavenumber to ',
     $        kaMaxFr(kW)-kaBlSize(kW)
        rL = kaMaxFr(kW)-kaBlSize(kW)
      END IF
      IF (rH .GT. kaMaxFr(kW)) THEN
        write(kStdWarn,*)'Error!!Setting max wavenumber to ',kaMaxFr(kW)
        rH = kaMaxFr(kW)
      END IF
      IF (rH .LT. (kaMinFr(1)+kaBlSize(1))) THEN
        !! input stop wavenumber waaaaaay tooooooooo low
        write(kStdWarn,*) 'Error!!Setting maximum wavenumber to ',
     $         kaMinFr(1)+kaBlSize(1)
        rH = kaMinFr(1)+kaBlSize(1)
      END IF

 1111 FORMAT(A1)

      RETURN
      END

c************************************************************************
c this subroutine sets the lower file ID limits, by going thru the list top
c to bottom and bottom to top
      SUBROUTINE LowerLimits(rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $ rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,raFiles,raDelta)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c INPUT
c rL         = user set lower bound
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c raDelta    = freq spacing
c raFiles    = int(file start freq)
c OUTPUT
c rL1,rL2    = freq in which rL is set down to, wne checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c rFileStartFrLo1 = file lower freq where  rL is found, when checking 
c                   list from t->b
c rFileStartFrLo2 = file lower freq where  rL is found, when checking 
c                    list from b->t
      REAL raFiles(kNumkCompT),raDelta(kNumkCompT)
      INTEGER iFileIDLo1,iFileIDLo2
      REAL rFileStartFrLo1,rFileStartFrLo2
      REAL rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

      INTEGER iInt

c now check which block rL falls in, starting from the extreme highest file ID
c original code
      iInt = kNumkCompT
      rL1 = rL
 11   CONTINUE
      IF ((iInt .GT. 1).AND.(rL1 .LT. raBlock(iInt))) THEN
        iInt = iInt-1
        GO TO 11
      END IF
      iFileIDLo1 = iInt
      rFileStartFrLo1 = raFiles(iInt)
      rL1 = raBlock(iInt)
      write(kStdWarn,*) 'set lower freq to start of kcomp block',rL1

c now check which block rL falls in, starting from the extreme lowest file ID
c new code
      rL2 =  rL
      iInt = 1
 12   CONTINUE
      IF ((iInt .LT. kNumkCompT).AND.(rL2 .GE. raBlockEnd(iInt))) THEN
        iInt = iInt+1
        GO TO 12
      END IF
      iFileIDLo2 = iInt
      rFileStartFrLo2 = raFiles(iInt)
      rL2 = raBlock(iInt)
      write(kStdWarn,*) 'reset lower freq to start of kcomp block',rL2

      IF (abs(rL1 - rL2) .GT. 0.01) THEN
        write(kStdWarn,*) 'hmmm, lower freq needs to be fixed ...'
      END IF

      RETURN
      END
c************************************************************************
c this subroutine sets the upper file ID limits, by going thru the list top
c to bottom and bottom to top
      SUBROUTINE UpperLimits(rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $     rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,raFiles,raDelta)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c INPUT
c rH         = user set lower bound
c raDelta    = file freq spacing
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c raFiles    = int(file start freq)
c OUTPUT
c rH1,rH2    = freq in which rL is set down to, when checking list t->b, b->t
c iFileIDHi1    = file ID in which rL is found, when checking list from t->b
c iFileIDHi2    = file ID in which rL is found, when checking list from b->t
c rFileStartFrHi1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c rFileStartFrHi2 = file lower freq where  rL is found, when checking list 
c                   from b->t
      REAL raFiles(kNumkCompT),raDelta(kNumkCompT)
      INTEGER iFileIDHi1,iFileIDHi2
      REAL rFileStartFrHi1,rFileStartFrHi2
      REAL rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

      INTEGER iInt

c now check which block rH falls in, starting from the extreme lowest file ID
c original code
      iInt=1
      rH1 = rH
 22   CONTINUE
      IF (iInt .EQ. kNumkCompT) THEN
         GOTO 24
ccc was raBlock(iInt+1)
      ELSEIF ((iInt .LT. kNumkCompT) .AND. 
     $ (rH1-2*raDelta(iInt) .GT. raBlockEnd(iInt))) THEN
c        print *,iInt,rH1-2*raDelta(iInt),raBlock(iInt),raBlockEnd(iInt)
        iInt = iInt+1
        GO TO 22
      END IF
 24   CONTINUE
      iFileIDHi1 = iInt
      rFileStartFrHi1 = raFiles(iInt)
      rH1 = raBlockEnd(iInt)
      write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH1

c now check which block rH falls in, starting from the extreme highest file ID
c new code
      iInt = kNumkCompT
      rH2 = rH
 23   CONTINUE
      IF ((iInt .GT. 1).AND.(rH2+2*raDelta(iInt) .LE. raBlock(iInt))) THEN
        iInt = iInt-1
        GO TO 23
      END IF
      iFileIDHi2 = iInt
      rFileStartFrHi2 = raFiles(iInt)
      rH2 = raBlockEnd(iInt)
      write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH2

      IF (abs(rH1 - rH2) .GT. 0.01) THEN
        write(kStdWarn,*) 'hmmm, upper freqs need to be fixed ...'
      END IF

      RETURN
      END

c************************************************************************
c this subroutine sets which file ID should be used as the lower limit
c remember tags depend on q r s ... = 1 2 3 .....
      SUBROUTINE LowerFileID(rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $  rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,
     $  iaActualTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,rFileStartFrLo)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c OUTPUT
c iFileIDLo     = final iFileIDLo
c iFinalIDLo = final rFileStartFrLo
c iTruthLo   = have we certainly set the file ID???
c INPUT
c rL         = user set lower bound
c rL1,rL2    = freq in which rL is set down to, when checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c iFileIDHi1    = file ID in which rH is found, when checking list from b->t
c iFileIDHi2    = file ID in which rH is found, when checking list from t->b
c rFileStartFrLo1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c rFileStartFrLo2 = file lower freq where  rL is found, when checking list 
c                   from b->t
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
c   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
c   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
      INTEGER iFileIDLo1,iFileIDLo2
      REAL rFileStartFrLo1,rFileStartFrLo2,rFileStartFrLo
      REAL rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iTruthLo,iFileIDHi1,iFileIDHi2
      INTEGER iFileIDLo

      INTEGER iDiff11,iDiff12,iDiff21,iDiff22

      iTruthLo=-1
      IF (iFileIDLo1 .EQ. iFileIDLo2) THEN !lower limit simple:equal file ID's
         iFileIDLo = iFileIDLo1
         rFileStartFrLo = rFileStartFrLo1
         iTruthLo=1
         rL = rL1
         write(kStdWarn,*) 'Equal lower limit : iFileIDLo = ',iFileIDLo
       END IF

      IF (iFileIDLo1 .NE. iFileIDLo2) THEN    !lower limit hard
        write(kStdWarn,*) 'Unequal lower limit : iFileIDLo1,2 = ',
     $                                  iFileIDLo1,iFileIDLo2
         write(kStdWarn,*)'thinking ... '
        write(kStdWarn,*) 'loTag',iFileIDLo1,iaActualTag(iFileIDLo1), 
     $                            iFileIDLo2,iaActualTag(iFileIDLo2)
        write(kStdWarn,*) 'hiTag',iFileIDHi1,iaActualTag(iFileIDHi1),
     $                            iFileIDHi2,iaActualTag(iFileIDHi2)
        IF (iaActualTag(iFileIDHi1) .EQ. iaActualTag(iFileIDHi2)) THEN
          !for UPPER limit, equal tags; so now see if either of the 
          !lower iFileIDLo1 or iFileIDLo2 fall in the same tag block
          IF (iaActualTag(iFileIDLo1) .EQ. iaActualTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE IF (iaActualTag(iFileIDLo2) .EQ. iaActualTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          ELSE IF ((iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)) .LT. 
     $               (iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2))) THEN
            !look for min diff between tags
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE
            !look for min diff between tags
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          END IF
        END IF

        IF (iaActualTag(iFileIDHi1) .NE. iaActualTag(iFileIDHi2)) THEN
          !for UPPER limit, unequal file tags
          !so look for minimum difference between TAGS
          iDiff11 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo1)
          iDiff12 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo2)
          iDiff21 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)
          iDiff22 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2)

          IF (iaActualTag(iFileIDLo1) .EQ. iaActualTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE IF (iaActualTag(iFileIDLo2) .EQ. iaActualTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          ELSE IF (iaActualTag(iFileIDLo1) .EQ. iaActualTag(iFileIDHi1)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE IF (iaActualTag(iFileIDLo2) .EQ. iaActualTag(iFileIDHi1)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          !tags are rather different : come to the first acceptable combination
          ELSE IF (iDiff11.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE IF (iDiff21.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo1
            rFileStartFrLo = rFileStartFrLo1
            iTruthLo=1
            rL = rL1
          ELSE IF (iDiff12.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          ELSE IF (iDiff22.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo2
            rFileStartFrLo = rFileStartFrLo2
            iTruthLo=1
            rL = rL2
          END IF
        END IF
      END IF

      iFileIDlo1 = iFileIDlo
      iFileIDlo2 = iFileIDlo
      rFileStartFrLo1 = rFileStartFrLo
      rFileStartFrLo2 = rFileStartFrLo

      write(kStdWarn,*) 'Low bound:itruthlo,iTag,iFileIDlo,rFileStartFrLo= '
      write(kStdWarn,*) itruthlo,kaTag(iaActualTag(itruthlo)),iFileIDlo,rFileStartFrLo

      RETURN
      END

c************************************************************************
c this subroutine sets which file ID should be used as the upper limit
      SUBROUTINE UpperFileID(rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $   rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,
     $   iaActualTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,rFileStartFrHi)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c OUTPUT
c iFileIDHi     = final iFileIDHi
c iFinalIDHi = final rFileStartFrHi
c iTruthHi   = have we certainly set the file ID???
c INPUT
c rL         = user set lower bound
c rL1,rL2    = freq in which rL is set doiwn to, wne checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c iFileIDHi1    = file ID in which rH is found, when checking list from b->t
c iFileIDHi2    = file ID in which rH is found, when checking list from t->b
c rFileStartFrLo1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c rFileStartFrLo2 = file lower freq where  rL is found, when checking list 
c                   from b->t
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
c   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
c   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
      INTEGER iFileIDHi1,iFileIDHi2
      REAL rFileStartFrHi1,rFileStartFrHi2,rFileStartFrHi
      REAL rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iTruthHi,iFileIDLo1,iFileIDLo2
      INTEGER iFileIDHi

      INTEGER iDiff11,iDiff12,iDiff21,iDiff22

      iTruthHi=-1
      IF (iFileIDHi1 .EQ. iFileIDHi2) THEN    !upper limit simple
        iFileIDHi = iFileIDHi1
        rFileStartFrHi = rFileStartFrHi1
        iTruthHi=1
        rH = rH1
        write(kStdWarn,*)'Equal upper limit : iFileIDHi = ',iFileIDHi
      END IF

      IF (iFileIDHi1 .NE. iFileIDHi2) THEN    !upper limit hard
        write(kStdWarn,*) 'Unequal upper limit : iFileIDHi1,2 = ',
     $            iFileIDHi1,iFileIDHi2
        write(kStdWarn,*)'thinking ... '
        write(kStdWarn,*) 'loTag',iFileIDLo1,iaActualTag(iFileIDLo1), 
     $                            iFileIDLo2,iaActualTag(iFileIDLo2)
        write(kStdWarn,*) 'hiTag',iFileIDHi1,iaActualTag(iFileIDHi1),
     $                            iFileIDHi2,iaActualTag(iFileIDHi2)

        IF (iaActualTag(iFileIDLo1) .EQ. iaActualTag(iFileIDLo2)) THEN
          !for LOWER limit, equal tags; so now see if either of the 
          !upper iFileIDHi1 or iFileIDHi2 fall in the same tag block
          IF (iaActualTag(iFileIDHi1) .EQ. iaActualTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE IF (iaActualTag(iFileIDHi2) .EQ. iaActualTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
          ELSE IF ((iaActualTag(iFileIDLo2)-iaActualTag(iFileIDHi1)) .LT. 
     $              (iaActualTag(iFileIDLo2)-iaActualTag(iFileIDHi2))) THEN
            !look for min diff between tags
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE
            !look for min diff between tags
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
          END IF
        END IF

        IF (iaActualTag(iFileIDLo1) .NE. iaActualTag(iFileIDLo2)) THEN
          !for LOWER limit, unequal file tags
          !so look for minimum difference between TAGS
          iDiff11 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo1)
          iDiff12 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo2)
          iDiff21 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)
          iDiff22 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2)

          IF (iaActualTag(iFileIDHi1) .EQ. iaActualTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE IF (iaActualTag(iFileIDHi2) .EQ. iaActualTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
          ELSE IF (iaActualTag(iFileIDHi1) .EQ. iaActualTag(iFileIDLo1)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE IF (iaActualTag(iFileIDHi2) .EQ. iaActualTag(iFileIDLo1)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
         !tags are rather different : come to the first acceptable combination
          ELSE IF (iDiff11.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE IF (iDiff21.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi1
            rFileStartFrHi = rFileStartFrHi1
            iTruthHi=1
            rH = rH1
          ELSE IF (iDiff12.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
          ELSE IF (iDiff22.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi2
            rFileStartFrHi = rFileStartFrHi2
            iTruthHi=1
            rH = rH2
          END IF
        END IF
      END IF

c now subtract delta(wavenumber) so that we don't have to do an additional
c set of calculations because of the additional wavenumber point
      write(kStdWarn,*)'High bound:itruthhi,iTag,iFileIDhi,rFileStartFrHi = '
      write(kStdWarn,*) itruthhi,kaTag(iaActualTag(itruthhi)),iFileIDhi,rFileStartFrHi
      write(kStdWarn,*) ' '

      RETURN
      END
c************************************************************************
