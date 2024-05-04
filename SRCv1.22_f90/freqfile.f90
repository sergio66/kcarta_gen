! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!******** THIS FILE CONTAINS ROUTINES THAT FIND FREQ CHUNKS ***** *******
!************************************************************************

MODULE freqfile

USE basic_common

IMPLICIT NONE

CONTAINS

! this function checks to see if the cross section data base exists for
! gasID, between freqs rL,rH
! if rL,rH < 0 then checking basic presence of gas in data base
    INTEGER FUNCTION iCheckXsecDataBase(iGasID,rL,rH,iTagIn,iErr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iGasID = GAS ID to be searched for in the database
! rL,rH  = freq start/stop points
! iErr   = error status (mainly associated with not finding the relevant file)
! iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
!          (not used if kXsecFormat < 0)
    INTEGER :: iGasID,iErr,iTagIn
    REAL :: rL,rH

! local variables
    INTEGER :: iIOUN,iFileErr,iID,iTag
    INTEGER :: iLine,iNpts,iTemps
    CHARACTER(160) :: caLine,caFName
    REAL :: rLower,rHigher

! assume GASID , freqs are wrong
    iCheckXsecDataBase = -1

 103 FORMAT('ERROR! number ',I5,' opening xsec database file : ',/,A84)
    caFName = kXsecParamFile
    iIOUN = kTempUnit
    OPEN(UNIT=iIOUN,FILE=caFName,STATUS='old',FORM='FORMATTED',IOSTAT=iFileErr)
    IF (kXsecFormat < 0) THEN
      ! ccccccccc this is the original format : read old style xsec.param file
      IF (iFileErr /= 0) THEN
        iErr = 0
        WRITE(kStdErr,103) iFileErr,caFName
        CALL DoSTOP
      END IF
      kTempUnitOpen = 1
              
      ! read file util GASID, freq bounds match found or EOF
 20   READ(iIOUN,5020,END=777) caLine
      IF (caLine(1:1) .EQ. '!') GOTO 20
 
      READ(caLine,*) iLine,iID,iNpts,iTemps,rLower,rHigher

      IF ((iID == iGasID) .AND. (rL < 0) .AND. (rH < 0)) THEN
        ! basic presence of gas tested for, and found
        ! this is the -1 option in XSCGAS
        iCheckXsecDataBase = 1
      END IF

      IF ((iID == iGasID) .AND. (rL >= rLower) .AND. (rH <= rHigher)) THEN
        ! presence of gas tested for, with freq bounds, and found well within
         iCheckXsecDataBase = 1
      END IF

      IF ((iID == iGasID) .AND. (rL <= rHigher) .AND. (rH >= rLower)) THEN
        ! presence of gas tested for, with freq bounds, and found
        iCheckXsecDataBase = 1
      END IF

      IF (iCheckXsecDataBase < 0) THEN
        GOTO 20
      END IF
              
  777 CONTINUE
    ELSE
      ! ccccccccc this is the new format : read comp.param style xsec.param file
      ! read file util GASID, freq bounds match found or EOF
 30   READ(iIOUN,5020,END=888) caLine
      READ(caLine,*) iID,rLower,rHigher,iTag

      IF ((iID == iGasID) .AND. (rL < 0) .AND. (rH < 0)) THEN
        ! basic presence of gas tested for, and found
        ! this is basically the -1 option in XSCGAS
        iCheckXSecDataBase = 1
      END IF

      IF ((iID == iGasID) .AND. (rL >= rLower) .AND. (rH <= rHigher)) THEN
        ! presence of gas tested for, with freq bounds, and found well within
        ! this is when we WANT to UNCOMPRESS files!
        IF (iTag == iTagIn) THEN   !actually uncompressing stuff
          iCheckXsecDataBase = 1
        END IF
      END IF

      IF (iCheckXSecDataBase < 0) THEN
        GOTO 30
      END IF
 888 CONTINUE

    END IF

    CLOSE(iIOUN)
    kTempUnitOpen = -1

    5020 FORMAT(A80)
    RETURN
    end FUNCTION iCheckXsecDataBase

!************************************************************************
! this function checks to see if the comp data file exists for
! gasID, between freqs rL,rH
! if rL,rH < 0 then checking basic presence of gas in data base
    INTEGER FUNCTION iCheckCompDataBase(iGasID,rL,rH,iTagIn,iErr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iGasID = GAS ID to be searched for in the database
! rL,rH  = freq start/stop points
! iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
! iErr   = error status (mainly associated with not finding the relevant file)
    INTEGER :: iGasID,iErr,iTagIn
    REAL :: rL,rH

! local variables
    INTEGER :: iIOUN,iFileErr,iID,iTag
    REAL :: rLower,rHigher
    CHARACTER(160) :: caLine,caFname

! assume GASID , freqs are wrong
    iCheckCompDataBase = -1

 103 FORMAT('ERROR! number ',I5,' opening comp database file : ',/,A84)
    caFName = kCompParamFile
    iIOUN = kTempUnit
    OPEN(UNIT=iIOUN,FILE=caFname,STATUS='old',FORM='FORMATTED',IOSTAT=iFileErr)
    IF (iFileErr /= 0) THEN
      iErr = 0
      WRITE(kStdErr,103) iFileErr,caFname
      CALL DoSTOP
    END IF
    kTempUnitOpen = 1

! read file until GASID, freq bounds match found or EOF
 20   READ(iIOUN,5020,END=777) caLine
      IF (caLine(1:1) .EQ. '!') GOTO 20

 5020 FORMAT(A80)
    READ(caLine,*) iID,rLower,rHigher,iTag

    IF ((iID == iGasID) .AND. (rL < 0) .AND. (rH < 0)) THEN
      ! basic presence of gas tested for, and found
      ! this is basically the -1 option in MOLGAS
      iCheckCompDataBase = 1
    END IF

    IF ((iID == iGasID) .AND. (rL >= rLower) .AND. (rH <= rHigher)) THEN
      ! presence of gas tested for, with freq bounds, and found well within
      ! this is when we WANT to UNCOMPRESS files!
      IF (iTag == iTagIn) THEN   !actually uncompressing stuff
        iCheckCompDataBase = 1
      END IF
    END IF

    IF (iCheckCompDataBase < 0) THEN
      GOTO 20
    END IF
          
  777 CONTINUE
    CLOSE(iIOUN)
    kTempUnitOpen = -1

    RETURN
    end FUNCTION iCheckCompDataBase

!************************************************************************
! this function checks to see if the GAS ID/frequency combination are in
! the compressed data base or in the xsec database
    SUBROUTINE DataBaseCheck(iGasID,raFreq,iTag,iActualTag,iDoAdd,iErr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iGasID   = GAS ID to be searched for in CompDataBase or XsecDataBase
! raFreq   = wavenumber array that contains present chunk of 25 cm-1
! iErr     = errors (associated with file I/O)
! iTag     = 1,2,3 telling us wavenumber spacing 0.001 0.0025 0.005 cm-1
! iDoAdd   = -1,+1 tells us if we add on current gas to current 10000 pts
    INTEGER :: iGasID,iErr,iTag,iActualTag,iDoAdd
    REAL :: raFreq(kMaxPts)

    iDoAdd = -1

  333 FORMAT('---> Warning! No comprsd data gasID ',I3,' in ',f10.4,' cm-1 chunk; tagIN = ', I3)
  334 FORMAT('---> Warning! No comprsd data heavywater gasID ',I3,' in ',f10.4,' cm-1 chunk')
  414 FORMAT('Comp data exists for HeavyWater (gas#) ',I3,f10.4,' cm-1 chunk')

    ! check to see if the k-comp file exists for the (GAS ID/freq) combination
    IF ((1 <= iGasID) .AND. (iGasID <= kGasComp)) THEN
      iDoAdd = iCheckCompDataBase(iGasID,raFreq(1),raFreq(kMaxPts),iActualTag,iErr)
      IF (iDoAdd < 0) THEN
        WRITE(kStdWarn,333) iGasID,raFreq(1),iActualTag
        WRITE(kStdWarn,*) ' '
      END IF
      GOTO 2000
    END IF

    IF (iGasID == kMaxGas) THEN
      !! this is heavy water
      IF (((raFreq(1) + 1.0) >= kWaterIsobandStart1*1.0)  .AND. &
         ((raFreq(kMaxPts) - 1.0) <= kWaterIsobandStop1*1.0)) THEN
        iDoAdd = +1
        WRITE(kStdWarn,414) iGasID,raFreq(1)
      ELSEIF (((raFreq(1) + 1.0) >= kWaterIsobandStart2*1.0)  .AND. &
          ((raFreq(kMaxPts) - 1.0) <= kWaterIsobandStop2*1.0)) THEN
        iDoAdd = +1
        WRITE(kStdWarn,414) iGasID,raFreq(1)
      ELSE
        iDoAdd = -1
        WRITE(kStdWarn,334) iGasID,raFreq(1)
      END IF
    END IF

!   check to see if the xsec file exists for the (GAS ID/freq) combination
    444 FORMAT('---> Warning! No xsec data for gasID ',I3, ' in ',f10.4,' cm-1 chunk')

    IF ((kGasXsecLo <= iGasID) .AND. (iGasID <= kGasXsecHi)) THEN
      iDoAdd = iCheckXsecDataBase(iGasID,raFreq(1),raFreq(kMaxPts),iActualTag,iErr)
      IF (iDoAdd < 0) THEN
         WRITE(kStdWarn,444) iGasID,raFreq(1)
         WRITE(kStdWarn,*) ' '
      END IF
      GOTO 2000
    END IF

 1000 FORMAT('No contribution to absorption cross sections from GAS ID = ',I2,/,'... skipping to next gas ...')
    IF ((kGasComp+1 <= iGasID) .AND. (iGasID <= kGasXsecLo-1)) THEN
      iDoAdd = -1
      WRITE(kStdWarn,1000) iGasID
      GOTO 2000
    END IF

! check to see if we need to add on the water continuum
    IF (kCKD >= 0) THEN
      IF ((iGasID >= kNewGasLo) .AND. (iGasID <= kNewGasHi)) THEN
        iDoAdd = 1
        GOTO 2000
      END IF
    END IF

 1010 FORMAT('Cannot add contribution of GAS ID = ',I2)
    IF ((iGasID < 1) .OR. (iGasID > kMaxGas)) THEN
      iErr = 1
      WRITE(kStdWarn,1010) iGasID
      CALL DoSTOP
    END IF

 2000 CONTINUE
    RETURN
    end SUBROUTINE DataBaseCheck

!************************************************************************
! this subroutine checks that kNumKcompT, as set in kcartaparam.f90
! agrees with what one would expect from kaMinFr,kaMaxFr

    SUBROUTINE Check_kaNum_ka100layerCloudType

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iI,iJ
    REAL :: rF,rG

!! in pre_defined we use ka100layerCloudType(kMaxUserSet)
!! but then in kcartaparam.f90 we define kMaxClouds
!! beter have kMaxUserSet > kMaxClouds else array bound problems!!!
    IF (kMaxUserSet < kMaxClouds) THEN
      write(kStdErr,*) 'because of comblock6, need kMaxUserSet > kMaxClouds'
      write(kStdErr,*) 'kMaxUserSet,kMaxClouds = ',kMaxUserSet,kMaxClouds
      CALL DoStop
    END IF

! kW is set in pre_definedVERYHIGHRES_IR.param (or pre_definedVERYHIGHRES_IR.param.f90)
    DO iJ = 1,kW
      rF = kMaxPts*kaFrStep(iJ)
      rG = kaBlSize(iJ)
      IF (abs(rF-kaBlSize(iJ)) > kaFrStep(iJ)/2.0) THEN
        write(kStdErr,*) 'iJ =  ',iJ
        write(kStdErr,*) 'kcartaparam.f90 claims kaBlSize(iJ) = ',rG
        write(kStdErr,*) 'while it should be ',rF
        CALL DoSTOP
      END IF
      iI = NINT((kaMaxFr(iJ)-kaMinFr(iJ))/rF)
      IF (iI /= kaNumKComp(iJ)) THEN
        write(kStdErr,*) 'iJ = ',iJ
        write(kStdErr,*) 'kcartaparam.f90 says that the number of '
        write(kStdErr,*) 'kCompressed files = kaNumKComp(iJ) = ',kaNumKComp(iJ)
        write(kStdErr,*) 'but based on following parameters, there should be iI = ',iI
        write(kStdErr,*) kaMinFr(iJ),kaMaxFr(iJ),kaFrStep(iJ),kMaxPts
        write(kStdErr,*) 'iI = iNT((kMaxFreq-kMinFreq)/(kMaxPts*kFreqStep))'
        CALL DoSTOP
      END IF
    END DO

    iI = 0
    DO iJ = 1,kW
      iI = iI+(kaNumkComp(iJ))
    END DO
    IF (iI /= kNumkCompT) THEN
      write(kStdErr,*) 'kcartaparam.f90 says kNumkCompT = ',kNumkCompT
      write(kStdErr,*) 'while it should be ',iI
      CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE Check_kaNum_ka100layerCloudType
           
!************************************************************************
! this subroutine uses the input file to set the file  number for the
! compressed data
! with the following three constraints :
!   1) the wavenumbers have to be between rFreqMin1 cm-1 and rFreqMax3 cm-1.
!   2) since each k-compressed file is XYZ cm-1 long, (iU-iL) <= XYZ
!   3) iL and iU have to fall in the same XYZ cm-1 block
! If all checks out, the relevant file block is set up, as is the
! wavenumber array
    SUBROUTINE GetFreq(rFrLow,rFrHigh, &
    iFileIDLo,iFileIDHi, &
    raBlock,raFiles, &
    iaTagIndex,iaActualTag, &
    raFileStep,iaList,iChunkTotal)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rFrLow    = lower frequency bound
! rFrHigh   = upper frequency bound
! iFileIDLo    = lower file index ID (between 1--kNumKcompT)
! iFileIDLo    = upper file index ID (between 1--kNumKcompT)
! iaTagIndex tells which TagIndex (1 .. 10) is associated with which file
!     (1,2,3,4 for k,p,q,r etc
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
! raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!          very useful so program easily knows r630_g1.dat etc
! raBlock  tells the current kComp wavenumber block
! raFileStep tells you the current wavenumber step size*10000
! iaList   has the final list of iChunkTotal files that should be uncompressed
    INTEGER :: iaList(kNumkCompT),iChunkTotal
    REAL :: raFileStep(kNumkCompT),raFiles(kNumkCompT)
    INTEGER :: iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT)
    REAL :: raBlock(kNumkCompT)
    REAL :: rFrLow,rFrHigh
    INTEGER :: iFileIDLo,iFileIDHi

! local variables
! rFileStartFrLo = lower file ID (605 -- 2805)
! rFileStartFrHi = upper file ID
    REAL :: rFileStartFrLo,rFileStartFrHi
    REAL :: rLow,rHigh,rTemp
    INTEGER :: iDummy

    write(kStdErr,*) ' '
    write(kStdWarn,*) ' '    
    
    write(kStdWarn,*) '*********************************************'
    write(kStdWarn,*) '**Setting freqs, fileID lower/upper bounds***'

    IF (rFrLow >= rFrHigh) THEN
      write(kStdWarn,*) 'Swapping frequencies found in *FRQNCY'
      rTemp = rFrLow
      rFrLow = rFrHigh
      rFrHigh = rTemp
    END IF
    write(kStdWarn,*) 'input freqs are  : ',rFrLow,rFrHigh

! remember : int(x)  === trunc(x)
! remember : nint(x) === round(x)
! if the bands are FIR3,FIR2,FIR1,IR,NIR1,NIR2 onwards then
    IF (rFrlow >= 140.0) THEN
      !! start/stop chunks are integer numbers, so use this fact
      rLow = rFrLow
      rHigh = rFrHigh
      !first round down rFrLow, and round up rFrHigh
      rFrLow = 1.0*INT(rFrLow)
      IF (abs(rFrHigh-1.0*INT(rFrHigh)) > 0.000000) THEN
        rFrHigh = rFrHigh+1.0
      END IF
      rFrHigh = 1.0*INT(rFrHigh)
      write(kStdWarn,*) 'rounded down/up rFrLow/rFrhigh to ',rFrLow,rFrHigh
    ELSE
      write(kStdWarn,*) ' --> for bands with wavenumbers < 140 cm-1, kCARTA'
      write(kStdWarn,*) ' --> compressed files are not on integer bdries'
      write(kStdWarn,*) ' --> so we hope you got things correct!'
    END IF

    rLow = rFrLow
    rHigh = rFrHigh

    IF (rLow >= rHigh) THEN
      rTemp = rLow
      rLow = rHigh
      rHigh = rTemp
    END IF

    CALL filebounds(rLow,rHigh, &
      rFileStartFrLo,rFileStartFrHi,iFileIDLo,iFileIDHi, &
      raBlock,raFiles, &
      iaTagIndex,iaActualTag, &
      raFileStep,iaList,iChunkTotal)

    rFrLow = rLow
    rFrHigh = rHigh
    write(kStdWarn,*) ' '
    write(kStdWarn,*) 'rFrLow,rFrHigh after = ',rFrLow,rFrHigh
    write(kStdWarn,*) 'num chunks = ',iChunkTotal
    write(kStdWarn,*) 'Start/Stop FileTags = ', &
    iaActualTag(iaList(1)),iaActualTag(iaList(iChunkTotal))
    write(kStdWarn,*) 'compfile 1: freq,FileID ',rFileStartFrLo,iFileIDLo
    write(kStdWarn,*) 'compfile N: freq,FileID ',rFileStartFrHi,iFileIDHi

    write (kStdWarn,*) 'File info : '
    write (kStdWarn,*) '#,rBlock,iFile,iTagIndex,iActualTag,iFileStep,iList'
    write (kStdWarn,*) '----------------------------------------------'
    DO iDummy = iFileIDLo,iFileIDHi
      write(kStdWarn,111) iDummy,raBlock(iDummy),raFiles(iDummy), &
        iaTagIndex(iDummy),iaActualTag(iDummy), &
        raFileStep(iDummy),iaList(iDummy-iFileIDLo+1)
    END DO
    111 FORMAT(I4,2(' ',F10.4),2(' ',I3),' ',F10.4,' ',I5)
     
    IF (iaActualTag(iaList(1)) /= iaActualTag(iaList(iChunkTotal))) THEN
      write(kStdWarn,*) 'Start file tag = ',iaActualTag(iaList(1))
      write(kStdWarn,*) 'Stop  file tag = ',iaActualTag(iaList(iChunkTotal))

      write(kStdErr,*) 'program requires you choose start/stop freqs'
      write(kStdErr,*) 'that only span one wavenumber spacing, as '
      write(kStdErr,*) 'defined at the bottom of kcartaparam.f90'
      DO iDummy = 1,kW
        write(kStdErr,333) kaMinFr(iDummy),kaMaxFr(iDummy),kaFrStep(iDummy),kaTag(iDummy)
      END DO
      write(kStdErr,*) 'please reset *FRQNCY and retry',rFrLow,rFrHigh
      CALL DoSTOP
    END IF

    333 FORMAT('rF1,rF2,d_f,Tag = ',2(f8.2,'  '),f12.7,' ',i4)

    write(kStdWarn,*)'****Finished freqs, fileID lower/upper bounds**'

    RETURN
    end SUBROUTINE GetFreq

!************************************************************************
! this subroutine calculates which of the relevant k-comp
! blocks are required for the problem ... this is for the user input files
! also makes sure the fres lie between kMinFreq1,kMaxFreq3 ... else it sets
! them to default values and asks if the user wishes to continue
! artificially set up the "imaginary" kNumkComp th point to account for 2805

! note we have to be smart about possibility of overlapping q,r,s eg
    SUBROUTINE filebounds(rL,rH, &
    rFileStartFrLo,rFileStartFrHi,iFileIDLo,iFileIDHi, &
    raBlock,raFiles, &
    iaTagIndex,iaActualTag,raFileStep,iaList,iChunkTotal)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rL,rH     = from namelist file, these are the lower and upper wavenumber
!             bounds : could be reset here
! rFileStartFrLo = from rL, lower file ID bound (in wavenumbers e.g. 605)
! rFileStartFrHi = from rH, upper file ID bound (in wavenumbers e.g. 2780)
! iFileIDLo    = from rL, lower fileID (e.g. 605 ==> 1+kNumKfile(1))
! iFileIDHi    = from rH, upper fileID(e.g. 2780 ==>kNumKcompT(1)+
!                                                            kNumKcompT(2))
! iaTagIndex tells which TagIndex (1 .. 10) is associated with which file
!     (1,2,3,4 for k,p,q,r etc
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
! raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!          very useful so program easily knows r630_g1.dat etc
! raFileStep tells you the current wavenumber step size*10000
! raBlock  tells the current kComp wavenumber block (ie start freq)

! iaList   has the final list of iChunkTotal files that should be uncompressed
    INTEGER :: iaList(kNumkCompT),iChunkTotal
    REAL :: raFileStep(kNumkCompT),raFiles(kNumkCompT)
    INTEGER :: iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT)
    REAL :: raBlock(kNumkCompT),rL,rH
    REAL :: rFileStartFrLo,rFileStartFrHi
    INTEGER :: iFileIDLo,iFileIDHi

! local variables
    REAL :: rFileStartFrLo1,rFileStartFrHi1
    REAL :: rFileStartFrLo2,rFileStartFrHi2
    INTEGER :: iFileIDLo1,iFileIDHi1,iFileIDLo2,iFileIDHi2
    REAL :: raBlockEnd(kNumkCompT),rL1,rH1,rL2,rH2,rEnd,raDelta(kNumkCompT)
    INTEGER :: iInt,iFound,iTruthLo,iTruthHi
    INTEGER :: iErr,iDummy,iJunk

! first check lower, upper limits
    CALL CheckLimits(rL,rH)

    iErr = 0
    iErr = sum(kaNumkComp(1:kW))
    IF (iErr /= kNumkCompT) THEN
      write(kStdErr,*) 'Error in kNumkCompT',iErr,kNumkCompT
      CALL DoStop
    END IF

! note there could be "overlaps" between the q r s files eg
!     < 1=500-550> <2=550-600>  <3=600-650> <4=650-700>
!                                  <5=605-630> <6=630-655> <7=655-680> .....
    iDummy = 0
    DO iFound = 1,kW
      DO iInt = 1,kaNumkComp(iFound)
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

! go thru raBlock and see where we think the start file ID, stop file ID
! should be set at
    CALL LowerLimits(rL,rL1,rL2,iFileIDLo1,iFileIDLo2, &
      rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,raFiles,raDelta)
    CALL UpperLimits(rH,rH1,rH2,iFileIDHi1,iFileIDHi2, &
      rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,raFiles,raDelta)

! now set the fileID's to be used
    CALL LowerFileID(rL,rL1,rL2,iFileIDLo1,iFileIDLo2, &
      rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd, &
      iaActualTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,rFileStartFrLo)
    CALL UpperFileID(rH,rH1,rH2,iFileIDHi1,iFileIDHi2, &
      rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd, &
      iaActualTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,rFileStartFrHi)

    IF ((iFileIDLo >= 1) .AND. (iFileIDLo <= kNumkCompT) .AND. &
      (iFileIDHi >= 1) .AND. (iFileIDHi <= kNumkCompT) .AND. &
      (iFileIDLo <= iFileIDHi)  .AND. &
      (iTruthLo > 0) .AND. (iTruthHi > 0)) THEN
      iErr = -1
    ELSE
      write(kStdErr,*) 'Error in setting file '
      write(kStdErr,*) 'lower/upper bounds from *FRQNCY'
      CALL DoSTOP
    END IF

    iJunk  =  kLongOrShort
    iJunk  =  2
    IF  (abs(iJunk) <= 1) THEN  !! verbose printing
      iFound = iaActualTag(1)
      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'iCnt  raBlock  raBlockEnd  iaTag'
      write(kStdWarn,*) '--------------------------------'
      DO iDummy = 1,kNumkCompT-1
        IF ((iDummy /= iFileIDLo) .AND. (iDummy /= iFileIDHi)) THEN
          write(kStdWarn,10) iDummy,raBlock(iDummy),raBlockEnd(iDummy),iaActualTag(iDummy)
        ELSE
          write(kStdWarn,11) iDummy,raBlock(iDummy),raBlockEnd(iDummy),iaActualTag(iDummy)
        END IF
        IF (iFound /=  iaActualTag(iDummy+1)) THEN
          write(kStdWarn,*) '--------------------------------'
          iFound =  iaActualTag(iDummy+1)
        END IF
      END DO
      iDummy = kNumkCompT
      IF ((iDummy /= iFileIDLo) .AND. (iDummy /= iFileIDHi)) THEN
        write(kStdWarn,10) iDummy,raBlock(iDummy),raBlockEnd(iDummy),iaActualTag(iDummy)
      ELSE
        write(kStdWarn,11) iDummy,raBlock(iDummy),raBlockEnd(iDummy),iaActualTag(iDummy)
      END IF
      write(kStdWarn,*) '--------------------------------'
      write(kStdWarn,*) ' '
    ELSE  !! quiet printing
      iFound = iaActualTag(1)
      write(kStdWarn,*) ' '
      DO iDummy = 1,kNumkCompT-1
        IF (iFound /=  iaActualTag(iDummy+1)) THEN
          iFound =  iaActualTag(iDummy+1)
        END IF
      END DO
      iDummy = kNumkCompT
    END IF
            
    10 FORMAT(I4,' ',F12.5,' ',F12.5,' ',I3)
    11 FORMAT(I4,' ',F12.5,' ',F12.5,' ',I3,' <<<<<<<<<<<<<<<')

! now create iaList as necessary
! remember that iFileID is basically an index setting equal to iDummy above
! thus if we know iFileID we know everything!!!!!
    IF (iaActualTag(iFileIDLo) == iaActualTag(iFileIDHi)) THEN
      !very easy everything is in either q or r or s database
      iChunkTotal = 0
      DO iInt = 1,(iFileIDHi-iFileIDLo+1)
        iChunkTotal = iChunkTotal+1
        iaList(iChunkTotal) = iFileIDLo+(iInt-1)     !save file ID
        IF  (abs(kLongOrShort) <= 1) THEN  !! verbose printing
          write(kStdWarn,*) iChunkTotal,iaActualTag(iaList(iChunkTotal)), &
              raBlock(iaList(iChunkTotal)),raBlockEnd(iaList(iChunkTotal))
        END IF
      END DO
    ELSE
      !very hard : mixing of q,r,s databases
      iChunkTotal = 0
      iChunkTotal = iChunkTotal+1
      iaList(iChunkTotal) = iFileIDLo                !save file ID
      rEnd = raBlockEnd(iaList(iChunkTotal))
      write(kStdWarn,*) iChunkTotal,iaActualTag(iaList(iChunkTotal)), &
      raBlock(iaList(iChunkTotal)),raBlockEnd(iaList(iChunkTotal))
      DO iInt = 2,(iFileIDHi-iFileIDLo+1)
        IF (rEnd < raBlockEnd(iFileIDLo+iInt-1)) THEN
          iChunkTotal = iChunkTotal+1
          iaList(iChunkTotal) = iFileIDLo+(iInt-1)   !save file ID
          write(kStdWarn,*) iChunkTotal,iaActualTag(iaList(iChunkTotal)), &
            raBlock(iaList(iChunkTotal)),raBlockEnd(iaList(iChunkTotal))
          rEnd = raBlockEnd(iaList(iChunkTotal))
        END IF
      END DO
    END IF

! now that we have the ActualTag eg 20 for r=0605-2830 cm-1, or
!                                   25 for s=2805-3305 cm-1, or
!                                   30 for m=4000-5000 cm-1,
! we have also mapped this back into an TagIndex

    RETURN
    end SUBROUTINE filebounds

!************************************************************************
! this subroutine checks to make sure rL,rH lie within the lower/upper bounds
! of wavenumbers in the kComp files
    SUBROUTINE CheckLimits(rL,rH)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rl,rH are the lower and upper limits the user sends in; this subroutine
!                       can change them

    REAL :: rL,rH

    IF (rL < kaMinFr(1)) THEN
      write(kStdWarn,*) 'Error!!Setting min wavenumber to ',kaMinFr(1)
      rL = kaMinFr(1)
    END IF
    IF (rL > (kaMaxFr(kW)-kaBlSize(kW))) THEN
      !! input start wavenumber waaaaaay tooooooooo high
      write(kStdWarn,*) 'Error!!Setting minimum wavenumber to ',kaMaxFr(kW)-kaBlSize(kW)
      rL = kaMaxFr(kW)-kaBlSize(kW)
    END IF
    IF (rH > kaMaxFr(kW)) THEN
      write(kStdWarn,*)'Error!!Setting max wavenumber to ',kaMaxFr(kW)
      rH = kaMaxFr(kW)
    END IF
    IF (rH < (kaMinFr(1)+kaBlSize(1))) THEN
      !! input stop wavenumber waaaaaay tooooooooo low
      write(kStdWarn,*) 'Error!!Setting maximum wavenumber to ',kaMinFr(1)+kaBlSize(1)
      rH = kaMinFr(1)+kaBlSize(1)
    END IF

    1111 FORMAT(A1)

    RETURN
    end SUBROUTINE CheckLimits

!************************************************************************
! this subroutine sets the lower file ID limits, by going thru the list top
! to bottom and bottom to top
    SUBROUTINE LowerLimits(rL,rL1,rL2,iFileIDLo1,iFileIDLo2, &
    rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd,raFiles,raDelta)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! INPUT
! rL         = user set lower bound
! raBlock    = file start freqs
! raBlockEnd = file stop freqs
! raDelta    = freq spacing
! raFiles    = int(file start freq)
! OUTPUT
! rL1,rL2    = freq in which rL is set down to, wne checking list t->b, b->t
! iFileIDLo1    = file ID in which rL is found, when checking list from t->b
! iFileIDLo2    = file ID in which rL is found, when checking list from b->t
! rFileStartFrLo1 = file lower freq where  rL is found, when checking
!                   list from t->b
! rFileStartFrLo2 = file lower freq where  rL is found, when checking
!                    list from b->t
    REAL :: raFiles(kNumkCompT),raDelta(kNumkCompT)
    INTEGER :: iFileIDLo1,iFileIDLo2
    REAL :: rFileStartFrLo1,rFileStartFrLo2
    REAL :: rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

    INTEGER :: iInt

! now check which block rL falls in, starting from the extreme highest file ID
! original code
    iInt = kNumkCompT
    rL1 = rL
  11 CONTINUE
    IF ((iInt > 1) .AND. (rL1 < raBlock(iInt))) THEN
      iInt = iInt-1
      GO TO 11
    END IF
    iFileIDLo1 = iInt
    rFileStartFrLo1 = raFiles(iInt)
    rL1 = raBlock(iInt)
    write(kStdWarn,*) 'set lower freq to start of kcomp block',rL1

! now check which block rL falls in, starting from the extreme lowest file ID
! new code
    rL2 =  rL
    iInt = 1
    12 CONTINUE
    IF ((iInt < kNumkCompT) .AND. (rL2 >= raBlockEnd(iInt))) THEN
      iInt = iInt+1
      GO TO 12
    END IF
    iFileIDLo2 = iInt
    rFileStartFrLo2 = raFiles(iInt)
    rL2 = raBlock(iInt)
    write(kStdWarn,*) 'reset lower freq to start of kcomp block',rL2

    IF (abs(rL1 - rL2) > 0.01) THEN
      write(kStdWarn,*) 'hmmm, lower freq needs to be fixed ...'
    END IF

    RETURN
    end SUBROUTINE LowerLimits
!************************************************************************
! this subroutine sets the upper file ID limits, by going thru the list top
! to bottom and bottom to top
    SUBROUTINE UpperLimits(rH,rH1,rH2,iFileIDHi1,iFileIDHi2, &
    rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd,raFiles,raDelta)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! INPUT
! rH         = user set lower bound
! raDelta    = file freq spacing
! raBlock    = file start freqs
! raBlockEnd = file stop freqs
! raFiles    = int(file start freq)
! OUTPUT
! rH1,rH2    = freq in which rL is set down to, when checking list t->b, b->t
! iFileIDHi1    = file ID in which rL is found, when checking list from t->b
! iFileIDHi2    = file ID in which rL is found, when checking list from b->t
! rFileStartFrHi1 = file lower freq where  rL is found, when checking list
!                   from t->b
! rFileStartFrHi2 = file lower freq where  rL is found, when checking list
!                   from b->t
    REAL :: raFiles(kNumkCompT),raDelta(kNumkCompT)
    INTEGER :: iFileIDHi1,iFileIDHi2
    REAL :: rFileStartFrHi1,rFileStartFrHi2
    REAL :: rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

    INTEGER :: iInt

! now check which block rH falls in, starting from the extreme lowest file ID
! original code
    iInt=1
    rH1 = rH
    22 CONTINUE
    IF (iInt == kNumkCompT) THEN
      GOTO 24
      ! c was raBlock(iInt+1)
    ELSEIF ((iInt < kNumkCompT) .AND. (rH1-2*raDelta(iInt) > raBlockEnd(iInt))) THEN
      !  print *,iInt,rH1-2*raDelta(iInt),raBlock(iInt),raBlockEnd(iInt)
      iInt = iInt+1
      GO TO 22
    END IF
  24 CONTINUE
    iFileIDHi1 = iInt
    rFileStartFrHi1 = raFiles(iInt)
    rH1 = raBlockEnd(iInt)
    write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH1

! now check which block rH falls in, starting from the extreme highest file ID
! new code
    iInt = kNumkCompT
    rH2 = rH
    23 CONTINUE
    IF ((iInt > 1) .AND. (rH2+2*raDelta(iInt) <= raBlock(iInt))) THEN
      iInt = iInt-1
      GO TO 23
    END IF
    iFileIDHi2 = iInt
    rFileStartFrHi2 = raFiles(iInt)
    rH2 = raBlockEnd(iInt)
    write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH2

    IF (abs(rH1 - rH2) > 0.01) THEN
      write(kStdWarn,*) 'hmmm, upper freqs need to be fixed ...'
    END IF

    RETURN
    end SUBROUTINE UpperLimits

!************************************************************************
! this subroutine sets which file ID should be used as the lower limit
! remember tags depend on q r s ... = 1 2 3 .....
    SUBROUTINE LowerFileID(rL,rL1,rL2,iFileIDLo1,iFileIDLo2, &
    rFileStartFrLo1,rFileStartFrLo2,raBlock,raBlockEnd, &
    iaActualTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,rFileStartFrLo)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! OUTPUT
! iFileIDLo     = final iFileIDLo
! iFinalIDLo = final rFileStartFrLo
! iTruthLo   = have we certainly set the file ID???
! INPUT
! rL         = user set lower bound
! rL1,rL2    = freq in which rL is set down to, when checking list t->b, b->t
! iFileIDLo1    = file ID in which rL is found, when checking list from t->b
! iFileIDLo2    = file ID in which rL is found, when checking list from b->t
! iFileIDHi1    = file ID in which rH is found, when checking list from b->t
! iFileIDHi2    = file ID in which rH is found, when checking list from t->b
! rFileStartFrLo1 = file lower freq where  rL is found, when checking list
!                   from t->b
! rFileStartFrLo2 = file lower freq where  rL is found, when checking list
!                   from b->t
! raBlock    = file start freqs
! raBlockEnd = file stop freqs
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
    INTEGER :: iFileIDLo1,iFileIDLo2
    REAL :: rFileStartFrLo1,rFileStartFrLo2,rFileStartFrLo
    REAL :: rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
    INTEGER :: iaActualTag(kNumkCompT),iTruthLo,iFileIDHi1,iFileIDHi2
    INTEGER :: iFileIDLo

    INTEGER :: iDiff11,iDiff12,iDiff21,iDiff22

    iTruthLo = -1
    IF (iFileIDLo1 == iFileIDLo2) THEN !lower limit simple:equal file ID's
      iFileIDLo = iFileIDLo1
      rFileStartFrLo = rFileStartFrLo1
      iTruthLo = 1
      rL = rL1
      write(kStdWarn,*) 'Equal lower limit : iFileIDLo = ',iFileIDLo
    END IF

    IF (iFileIDLo1 /= iFileIDLo2) THEN    !lower limit hard
      write(kStdWarn,*) 'Unequal lower limit : iFileIDLo1,2 = ',iFileIDLo1,iFileIDLo2
      write(kStdWarn,*)'thinking ... '
      write(kStdWarn,*) 'loTag',iFileIDLo1,iaActualTag(iFileIDLo1),iFileIDLo2,iaActualTag(iFileIDLo2)
      write(kStdWarn,*) 'hiTag',iFileIDHi1,iaActualTag(iFileIDHi1),iFileIDHi2,iaActualTag(iFileIDHi2)
      IF (iaActualTag(iFileIDHi1) == iaActualTag(iFileIDHi2)) THEN
        ! for UPPER limit, equal tags; so now see if either of the
        ! lower iFileIDLo1 or iFileIDLo2 fall in the same tag block
        IF (iaActualTag(iFileIDLo1) == iaActualTag(iFileIDHi2)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSEIF (iaActualTag(iFileIDLo2) == iaActualTag(iFileIDHi2)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
          rL = rL2
        ELSEIF ((iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)) < &
                (iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2))) THEN
          !look for min diff between tags
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSE
          !look for min diff between tags
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
          rL = rL2
        END IF
      END IF

      IF (iaActualTag(iFileIDHi1) /= iaActualTag(iFileIDHi2)) THEN
        !for UPPER limit, unequal file tags
        !so look for minimum difference between TAGS
        iDiff11 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo1)
        iDiff12 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo2)
        iDiff21 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)
        iDiff22 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2)

        IF (iaActualTag(iFileIDLo1) == iaActualTag(iFileIDHi2)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSEIF (iaActualTag(iFileIDLo2) == iaActualTag(iFileIDHi2)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
          rL = rL2
        ELSEIF (iaActualTag(iFileIDLo1) == iaActualTag(iFileIDHi1)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSEIF (iaActualTag(iFileIDLo2) == iaActualTag(iFileIDHi1)) THEN
          !these two tags are the same
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
          rL = rL2
        ELSEIF (iDiff11 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          !tags are rather different : come to the first acceptable combination
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSEIF (iDiff21 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDLo = iFileIDLo1
          rFileStartFrLo = rFileStartFrLo1
          iTruthLo = 1
          rL = rL1
        ELSEIF (iDiff12 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
          rL = rL2
        ELSEIF (iDiff22 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDLo = iFileIDLo2
          rFileStartFrLo = rFileStartFrLo2
          iTruthLo = 1
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
    end SUBROUTINE LowerFileID

!************************************************************************
! this subroutine sets which file ID should be used as the upper limit
    SUBROUTINE UpperFileID(rH,rH1,rH2,iFileIDHi1,iFileIDHi2, &
    rFileStartFrHi1,rFileStartFrHi2,raBlock,raBlockEnd, &
    iaActualTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,rFileStartFrHi)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! OUTPUT
! iFileIDHi     = final iFileIDHi
! iFinalIDHi = final rFileStartFrHi
! iTruthHi   = have we certainly set the file ID???
! INPUT
! rL         = user set lower bound
! rL1,rL2    = freq in which rL is set doiwn to, wne checking list t->b, b->t
! iFileIDLo1    = file ID in which rL is found, when checking list from t->b
! iFileIDLo2    = file ID in which rL is found, when checking list from b->t
! iFileIDHi1    = file ID in which rH is found, when checking list from b->t
! iFileIDHi2    = file ID in which rH is found, when checking list from t->b
! rFileStartFrLo1 = file lower freq where  rL is found, when checking list
!                   from t->b
! rFileStartFrLo2 = file lower freq where  rL is found, when checking list
!                   from b->t
! raBlock    = file start freqs
! raBlockEnd = file stop freqs
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
    INTEGER :: iFileIDHi1,iFileIDHi2
    REAL :: rFileStartFrHi1,rFileStartFrHi2,rFileStartFrHi
    REAL :: rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
    INTEGER :: iaActualTag(kNumkCompT),iTruthHi,iFileIDLo1,iFileIDLo2
    INTEGER :: iFileIDHi

    INTEGER :: iDiff11,iDiff12,iDiff21,iDiff22

    iTruthHi=-1
    IF (iFileIDHi1 == iFileIDHi2) THEN    !upper limit simple
      iFileIDHi = iFileIDHi1
      rFileStartFrHi = rFileStartFrHi1
      iTruthHi = 1
      rH = rH1
      write(kStdWarn,*)'Equal upper limit : iFileIDHi = ',iFileIDHi
    END IF

    IF (iFileIDHi1 /= iFileIDHi2) THEN    !upper limit hard
      write(kStdWarn,*) 'Unequal upper limit : iFileIDHi1,2 = ',iFileIDHi1,iFileIDHi2
      write(kStdWarn,*)'thinking ... '
      write(kStdWarn,*) 'loTag',iFileIDLo1,iaActualTag(iFileIDLo1),iFileIDLo2,iaActualTag(iFileIDLo2)
      write(kStdWarn,*) 'hiTag',iFileIDHi1,iaActualTag(iFileIDHi1),iFileIDHi2,iaActualTag(iFileIDHi2)

      IF (iaActualTag(iFileIDLo1) == iaActualTag(iFileIDLo2)) THEN
        ! or LOWER limit, equal tags; so now see if either of the
        ! pper iFileIDHi1 or iFileIDHi2 fall in the same tag block
        IF (iaActualTag(iFileIDHi1) == iaActualTag(iFileIDLo2)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSEIF (iaActualTag(iFileIDHi2) == iaActualTag(iFileIDLo2)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        ELSEIF ((iaActualTag(iFileIDLo2)-iaActualTag(iFileIDHi1)) < &
              (iaActualTag(iFileIDLo2)-iaActualTag(iFileIDHi2))) THEN
          !look for min diff between tags
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSE
          !look for min diff between tags
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        END IF
      END IF

      IF (iaActualTag(iFileIDLo1) /= iaActualTag(iFileIDLo2)) THEN
        ! or LOWER limit, unequal file tags
        ! so look for minimum difference between TAGS
        iDiff11 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo1)
        iDiff12 = iaActualTag(iFileIDHi1)-iaActualTag(iFileIDLo2)
        iDiff21 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo1)
        iDiff22 = iaActualTag(iFileIDHi2)-iaActualTag(iFileIDLo2)

        IF (iaActualTag(iFileIDHi1) == iaActualTag(iFileIDLo2)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSEIF (iaActualTag(iFileIDHi2) == iaActualTag(iFileIDLo2)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        ELSEIF (iaActualTag(iFileIDHi1) == iaActualTag(iFileIDLo1)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSEIF (iaActualTag(iFileIDHi2) == iaActualTag(iFileIDLo1)) THEN
          !these two tags are the same
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        ELSEIF (iDiff11 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          !tags are rather different : come to the first acceptable combination
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSEIF (iDiff21 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDHi = iFileIDHi1
          rFileStartFrHi = rFileStartFrHi1
          iTruthHi = 1
          rH = rH1
        ELSEIF (iDiff12 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        ELSEIF (iDiff22 == min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
          iFileIDHi = iFileIDHi2
          rFileStartFrHi = rFileStartFrHi2
          iTruthHi = 1
          rH = rH2
        END IF
      END IF
    END IF

! now subtract delta(wavenumber) so that we don't have to do an additional
! set of calculations because of the additional wavenumber point
    write(kStdWarn,*)'High bound:itruthhi,iTag,iFileIDhi,rFileStartFrHi = '
    write(kStdWarn,*) itruthhi,kaTag(iaActualTag(itruthhi)),iFileIDhi,rFileStartFrHi
    write(kStdWarn,*) ' '

    RETURN
    end SUBROUTINE UpperFileID
!************************************************************************

END MODULE FreqFile
