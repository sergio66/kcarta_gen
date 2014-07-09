c Copyright 2000
c University of Maryland Baltimore County 
c All Rights Reserved

c this file deals with reading in the user file
c main routines for PTHFIL are here 
c also, the REFERENCE profile reading subroutine is here

c************************************************************************
c this subroutine deals with the 'PTHFIL' keyword 
      SUBROUTINE pthfil4(raaAmt,raaTemp,raaPress,raaPartPress,caPFName,iRTP,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath,
     $      iProfileLayers,raPressLevels,raThickness)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caPFName = character*80 profile name
c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gas ID's should be read in
c iaWhichGasRead = array that tracks whch gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c raLayerHeight = heights of layers in KM!!!!!!!
c iRTP  = if RTP KLAYERS profile, which one of the profiles to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
c iProfileLayers = tells how many layers read in from RTP or KLAYERS file
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayerHeight(kProfLayer)
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      INTEGER iNpath,iRTP
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
 5010 FORMAT('Error reading section ',A7,' of main user/PTHFIL file')

      !!!kRTP = -1 : read old style kLAYERS profile; set atm from namelist
      !!!kRTP =  0 : read RTP style kLAYERS profile; set atm from namelist
      !!!kRTP = +1 : read RTP style kLAYERS profile; set atm from RTP file
      IF (kRTP .LT. 0) THEN
        write(kStdWarn,*) 'old style Klayers Profile file to be read is  : '
        write(kStdWarn,*) caPfname
        CALL readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,raPressLevels,raThickness)
        iProfileLayers = kProfLayer !!!!!old style always has kProfLayer layers
      ELSEIF (kRTP .GE. 0) THEN
        write(kStdWarn,*) 'new style RTP profile to be read is  : '
        write(kStdWarn,*) caPfname
        write(kStdWarn,*) 'within this file, we will read profile # ',iRTP
        CALL readRTP(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPressLevels,raThickness)
        END IF

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
     $      iNpath,caPfName,raPressLevels,raThickness)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in KM!!!!!!
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
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

      INTEGER iL1,iL2,length130,iI
      CHARACTER*130 ca1,ca2,caTemp

      DO iI = 1,kProfLayer
        raPressLevels(iI) = 0.0
        raThickness(iI) = 0.0
        END DO
      raPressLevels(kProfLayer+1) = 0.0

      ca1 = '! TOTAL NO. OF LAYERS IN ATM:'
      ca2 = '! TOTAL NO. OF GASES:'
      iL1 = length130(ca1) 
      iL2 = length130(ca2) 

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
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error reading section ',A7,' of main user/PTHFIL file')

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

      IF (caStr(1:iL1) .EQ. ca1) THEN
        DO iP=1,130
          caTemp(iP:iP) = ' '
          END DO
        caTemp(1:10) = caStr(iL1+1:iL1+10)
        read (caTemp,*) iL1         !!!now we know how many KLAYERS in profile
        END IF

      IF (caStr(1:iL2) .EQ. ca2) THEN
        DO iP=1,130
          caTemp(iP:iP) = ' '
          END DO
        caTemp(1:10) = caStr(iL2+1:iL2+10)
        read (caTemp,*) iL2         !!!now we know how many gases in profile
        END IF

      IF (caStr(1:1) .EQ. '!') THEN
        GO TO 30
      ELSE
        iNumLinesRead=iNumLinesRead+1
        READ (caStr,*,ERR=13,END=13) iNpath
        WRITE(kStdWarn,3000) iNpath,iNeed2Read*kProfLayer
 3000   FORMAT('input file has ',I5,' paths; kCARTA needs ',I5,' paths')
        iGasesInProfile=INT(iNpath*1.0/(kProfLayer*1.0))
        END IF

      !!!now check if this agrees with iL1,iL2 above
      IF (kProfLayer .NE. iL1) THEN
        write (kStdErr,*) 'KLAYERS profile has ',iL2,' gases in atm'
        write (kStdErr,*) 'KLAYERS profile has ',iL1,' layers in atm'
        write (kStdErr,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdErr,*) 'Need kProfLayer == Layers in Profile'
        write (kStdErr,*) 'Either redo klayers.x or kcarta.x'
        CALL DoStop
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
           raaHeight(iaNpathCounter(iIDgas),iGasIndex)=rH*1000  !!change to m!

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

c now reread the profile file, so that we can get info about presslevels
c and layer thicknesses
c just turn this off to read the real old klayers files (eg that Scott uses for
c his kcartav1.07 runs ....)
      CALL GetMoreInfo(raPressLevels,raThickness,caPfName)

      RETURN
      END

c************************************************************************
c this subroutine reads the info from the KLAYERS profile, storing info
c about layer thicknesses and presslevels
      SUBROUTINE GetMoreInfoOld(raPressLevels,raThickness,caPfName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c parameters
      CHARACTER*80 caPFname
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)

c local variables      
      INTEGER iIOUN2,iErrIO,iErr
      CHARACTER*130 caTemp,caStr
      CHARACTER c1a,c1b
      CHARACTER*3 c3
      INTEGER iI,iJ,iNumP,iNumH,iTrue
      REAL p1,p2,raP1(kProfLayer+1),raP2(kProfLayer+1)
      REAL h1,raH1(kProfLayer)

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

      iNumP = 0
      iNumH = 0
 40   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr

      IF ((iNumP .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        IF (caStr(3:6) .EQ. 'PATH') THEN
          iNumP = iNumP + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            END DO
          iI = 3
 10       CONTINUE
          IF (caStr(iI:iI+2) .NE. 'km,') THEN
            iI = iI+1
            GOTO 10
            END IF
          iI = iI + 3
          iJ = 130-iI+1
          caTemp(1:iJ) = caStr(iI:130)
          read(caTemp,*) p1,c1b,p2
          raP1(iNumP) = p1
          raP2(iNumP) = p2
          END IF
        END IF

      IF ((iNumH .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        IF (caStr(3:5) .EQ. 'RAY') THEN
          iNumH = iNumH + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            END DO
          iI = 3
 20       CONTINUE
          IF (caStr(iI:iI) .NE. '=') THEN
            iI = iI+1
            GOTO 20
            END IF
          iI = iI + 1
          iJ = 130-iI+1
          caTemp(1:iJ) = caStr(iI:130)
          read(caTemp,*) h1
          raH1(iNumH) = h1
          END IF
        END IF

      IF ((iNumP .LT. kProfLayer) .OR. (iNumH .LT. kProfLayer)) THEN
        GOTO 40
        END IF

      CLOSE(iIOUN2)

      kProfileUnitOpen=-1
      raP1(kProfLayer+1) = raP2(kProfLayer)  !!set highest pressure level

c check to see that raP1(1) > raP1(2) >> raP1(kProfLayer)
      IF (raP1(1) .LT. raP1(2)) THEN         
        !swap pressure levels
        DO iI = 1,kProfLayer+1
          raP2(iI) = raP1(kProfLayer+1 - iI + 1) 
          END DO
        DO iI = 1,kProfLayer+1
          raP1(iI) = raP2(iI)
          END DO
        !swap layer thickness
        DO iI = 1,kProfLayer
          raP2(iI) = raH1(kProfLayer - iI + 1) 
          END DO
        DO iI = 1,kProfLayer
          raH1(iI) = raP2(iI)
          END DO
        END IF

c double check to see that raP1(1) > raP1(2) >> raP1(kProfLayer)
      iTrue = 1       !!!assume all is ok
      DO iI = 1,kProfLayer
        IF (raP1(iI) .LT. raP1(iI+1)) THEN
          iTrue = -1
          END IF
        END DO
      IF (iTrue .LT. 0) THEN
        write(kStdErr,*) 'Pressure levels from TEXT klayers file are not in'
        write(kStdErr,*) 'monotonically decreasing or increasing order!!!!!!'
        CALL DoStop
        END IF

      DO iI = 1,kProfLayer+1
        raPresslevels(iI) = raP1(iI)
        END DO

      DO iI = 1,kProfLayer
        raThickness(iI) = raH1(iI)*1000           !change from km to m
        END DO

      GOTO 14

 5030 FORMAT(A130)
 13   write(kStdErr,*) 'error reading text kLAYERS profile for layer thickness'
      write(kStdErr,*) '   or reading text kLAYERS profile for pressure levels'
      CALL DOStop

 14   CONTINUE
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (KLAYERS file)'
      write (kStdWarn,*) '---------------------------------------------------'
      write (kStdWarn,*) 'Lowest  z : press levels (mb) = ',raP1(1),raP1(2)
      write (kStdWarn,*) 'Highest z : press levels (mb) = ',raP1(kProfLayer),
     $                                                  raP1(kProfLayer+1)
      write (kStdWarn,*) 'Lowest  z : layer thick (km) = ',raH1(1),raH1(2)
      write (kStdWarn,*) 'Highest z : layer thick (km) = ',raH1(kProfLayer-1),
     $                                                    raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(1)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this subroutine reads the info from the KLAYERS profile, storing info
c about layer thicknesses and presslevels
      SUBROUTINE GetMoreInfo(raPressLevels,raThickness,caPfName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c parameters
      CHARACTER*80 caPFname
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)

c local variables      
      INTEGER iIOUN2,iErrIO,iErr,SubString,iS
      CHARACTER*130 caTemp,caTemp2,caStr
      CHARACTER c1a,c1b
      CHARACTER*10 ca10a,ca10b,ca10c,ca10d,ca10e
      CHARACTER*3 c3
      INTEGER iI,iJ,iNumP,iNumH,iTrue,iK,iL
      REAL p1,p2,raP1(kProfLayer+1),raP2(kProfLayer+1)
      REAL h1,raH1(kProfLayer)

      DO iI=1,5
        ca10a(iI:iI) = ' '
        ca10b(iI:iI) = ' '
        ca10c(iI:iI) = ' '
        ca10d(iI:iI) = ' '
        END DO
      ca10a(1:4) = 'PATH'
      ca10b(1:3) = 'km,'
      ca10c(1:3) = ' - '
      ca10d(1:3) = 'RAY'
      ca10e(1:1) = '='

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

      iNumP = 0
      iNumH = 0
 40   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr

      IF ((iNumP .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        iS = SubString(caStr,ca10a,3,6)
        IF (iS .EQ. 1) THEN
          iNumP = iNumP + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            caTemp2(iI:iI) = ' '
            END DO
          iI = 3
 10       CONTINUE
          iS = SubString(caStr,ca10b,iI,iI+2)
          IF (iS .EQ. -1) THEN
            iI = iI+1
            GOTO 10
            END IF
          iI = iI + 3
          iJ = 130-iI
          iL = iI
          DO iK = 1,iJ
            caTemp(iK:iK) = caStr(iL:iL)
            iL=iL+1
            END DO
          read(caTemp,*) p1
          iI = 1
 15       CONTINUE
          iS = SubString(caTemp,ca10c,iI,iI+2)
          IF (iS .EQ. -1) THEN
            iI = iI+1
            GOTO 15
            END IF
          iI = iI + 3
          iJ = 130-iI
          iL = iI
          !caTemp(1:iJ) = caStr(iI:130)
          DO iK = 1,iJ
            caTemp2(iK:iK) = caTemp(iL:iL)
            iL=iL+1
            END DO
          read(caTemp2,*) p2
          raP1(iNumP) = p1
          raP2(iNumP) = p2
          END IF
        END IF

      IF ((iNumH .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        iS = SubString(caStr,ca10d,3,5)
        IF (iS .EQ. 1) THEN
          iNumH = iNumH + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            END DO
          iI = 3
 20       CONTINUE
          iS = SubString(caStr,ca10e,iI,iI)
          IF (iS .EQ. -1) THEN
            iI = iI+1
            GOTO 20
            END IF
          iI = iI + 1
          iJ = 130-iI+1
          caTemp(1:iJ) = caStr(iI:130)
          read(caTemp,*) h1
          raH1(iNumH) = h1
          END IF
        END IF

      IF ((iNumP .LT. kProfLayer) .OR. (iNumH .LT. kProfLayer)) THEN
        GOTO 40
        END IF

      CLOSE(iIOUN2)

      kProfileUnitOpen=-1
      raP1(kProfLayer+1) = raP2(kProfLayer)  !!set highest pressure level

c check to see that raP1(1) > raP1(2) >> raP1(kProfLayer)
      IF (raP1(1) .LT. raP1(2)) THEN         
        !swap pressure levels
        DO iI = 1,kProfLayer+1
          raP2(iI) = raP1(kProfLayer+1 - iI + 1) 
          END DO
        DO iI = 1,kProfLayer+1
          raP1(iI) = raP2(iI)
          END DO
        !swap layer thickness
        DO iI = 1,kProfLayer
          raP2(iI) = raH1(kProfLayer - iI + 1) 
          END DO
        DO iI = 1,kProfLayer
          raH1(iI) = raP2(iI)
          END DO
        END IF

c double check to see that raP1(1) > raP1(2) >> raP1(kProfLayer)
      iTrue = 1       !!!assume all is ok
      DO iI = 1,kProfLayer
        IF (raP1(iI) .LT. raP1(iI+1)) THEN
          iTrue = -1
          END IF
        END DO
      IF (iTrue .LT. 0) THEN
        write(kStdErr,*) 'Pressure levels from TEXT klayers file are not in'
        write(kStdErr,*) 'monotonically decreasing or increasing order!!!!!!'
        CALL DoStop
        END IF

      DO iI = 1,kProfLayer+1
        raPresslevels(iI) = raP1(iI)
        END DO

      DO iI = 1,kProfLayer
        raThickness(iI) = raH1(iI)*1000           !change from km to m
        END DO

      GOTO 14

 5030 FORMAT(A130)
 13   write(kStdErr,*) 'error reading text kLAYERS profile for layer thickness'
      write(kStdErr,*) '   or reading text kLAYERS profile for pressure levels'
      CALL DOStop

 14   CONTINUE
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (KLAYERS file)'
      write (kStdWarn,*) '---------------------------------------------------'
      write (kStdWarn,*) 'Lowest  z : press levels (mb) = ',raP1(1),raP1(2)
      write (kStdWarn,*) 'Highest z : press levels (mb) = ',raP1(kProfLayer),
     $                                                  raP1(kProfLayer+1)
      write (kStdWarn,*) 'Lowest  z : layer thick (km) = ',raH1(1),raH1(2)
      write (kStdWarn,*) 'Highest z : layer thick (km) = ',raH1(kProfLayer-1),
     $                                                    raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(1)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this function checks to see if substring is in string
      INTEGER FUNCTION SubString(caX,ca10,i1,i2)

      CHARACTER*130 caX
      CHARACTER*10  ca10
      INTEGER i1,i2

      INTEGER K,iI,iJ,iLen

      K = -1           !!!asuume ca10 is in caX not at specified locations

      iLen = i2-i1 + 1
      iI = 1
      IJ = i1
 10   CONTINUE
      IF (caX(iJ:iJ) .EQ. ca10(iI:iI)) THEN
        K = 1
        iI = iI + 1
        iJ = iJ + 1
        IF (iJ .LT. iLen) THEN
          GOTO 10
          END IF
      ELSE
        K = -1
        END IF

      SubString = K

      RETURN
      END
c************************************************************************

c this subroutine reads the profile files (for the references)
c it flags an error if kProfLayers layers are not read in
c ProX === A=amount,T=temperature,P=pressure,PP=partial pressure
       SUBROUTINE ReadRefProf(caFname,iIOUN,iNlayIn,raProA,raProT,
     $      raProP,raProPP,iErrOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caFName   = name of file that has the profile
c iIOUN     = file IO UNIT number
c iNlay     = number of layers that are read in
c raProA/T/P/PP = profile amout, temperature,pressure,partial pressure
c iErrOut   = error count (usually associated with file I/O)
       CHARACTER*80 caFname
       INTEGER iIOUN,iNlayIn
       INTEGER iErrOut 
       REAL raProA(*),raProT(*)
       REAL raProP(*),raProPP(*)

c local variables
       INTEGER iErr,iJmin,iJmax,iJ,iNlay
       CHARACTER*100 caLine

       OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErr)
       IF (iErr .NE. 0) THEN
          WRITE(kStdErr,1080) iErr, caFname
 1080     FORMAT('ERROR! number ',I5,' opening profile file:',/,A82)
          CALL DoSTOP
       ENDIF
       kProfileUnitOpen=1
       
       iJMin = +1000000
       iJMax = -1000000

C      Read the file (skip comment lines)
       iNlay=0
 20    READ(iIOUN,5020,END=199) caLine
 5020  FORMAT(A100)
       IF (caLine(1:1) .NE. '!') THEN
         iNlay=iNlay+1
         READ(caLine,*) iJ,raProP(iNlay),raProPP(iNlay),
     $                   raProT(iNlay),raProA(iNlay)         
         ENDIF
       GOTO 20
C      Close the file
 199   CLOSE(iIOUN)
       kProfileUnitOpen=-1

c check to see that ALL layers have been read in (could be greater than 100)
       IF (iNlay .NE. iNlayIN) THEN
         iErrOUt=1
         WRITE(kStdErr,500) caFName,iNLay,iNLayIN
 500     FORMAT ('Profile File',/,A82,/,' has ',I4,' layers (need ',I4,')')
         CALL DoSTOP
         END IF

       RETURN
       END

c************************************************************************
      SUBROUTINE rightpad130(caName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

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
c this finds the length of a string
      INTEGER FUNCTION length130(caName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

      CHARACTER*130 caName

      CHARACTER*130 caTempName
      INTEGER iR,iL,iInt

      CALL rightpad130(caName)

c find the "right" length of the input root name
      iR=len(caName)
 11   continue
      IF (caName(iR:iR) .eq. ' ') THEN
        iR=iR-1
        GO TO 11      
        END IF

      length130 = iR

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

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caStr130     = current line to be processed. If not enough info found in it,
c                additional lines will be read in as necessary
c iNpmixTemp   = total number of mixed paths to be read in / kProfLayer
c              = number of sets of mixed paths to be read in
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

c      DO iInt = 1,10
c        write(*,6000) caaMixFileLines(iInt)
c        END DO
c 6000 FORMAT(A80)

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
      IF (rC .gt. 0.0) THEN         !we know user is listing gas weights
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
c this subroutine does NOT use kDumbFile
      SUBROUTINE readlist(caaM,iType,iNum,iLines,raStore,iaStore)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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

      INTEGER iI,iJ,iAdd,iLinesStop
      REAL raStore2(kMaxGas)
      CHARACTER*130 caaTemp(kProfLayer)
      
c      DO iI=1,kProfLayer
c        caaTemp(iI) = caaM(iI)
c        END DO
c 30   CONTINUE
c      print *,'enter how many to read? '
c      read *,iAdd
c      read(caaTemp,*)(raStore(iI),iI=1,iAdd)
c      print *,(raStore(iI),iI=1,iAdd)     
c      print *,'again (-1 to end)  '
c      read *,iAdd
c      IF (iADD .GT. 0) GOTO 30
c      stop

      DO iI=1,kMaxGas
        iaStore(iI)=-1
        raStore(iI)=-1.0
        raStore2(iI)=-1.0
        END DO

      IF (iType .LT. 0) THEN
        !read in list of integer (discard) followed by iNum reals, using 
        !current line as the line to start reading stuff from
        iLines     = iLines
        iLinesStop = iLines-1
      ELSEIF (iType .GT. 0) THEN
        !read in list of (integer,real) pairs,  using next line as the 
        !line to start reading stuff from
        iLines     = iLines+1
        iLinesStop = iLines-1
        END IF

      iAdd = -1
c initialize the temporary ccaArray, to read just the required number of things
c increment iLinesStop by 1 each time this loop is called
 5    iLinesStop=iLinesStop+1
      DO iI=1,(iLinesStop-iLines+1)
        caaTemp(iI) = caaM(iLines+(iI-1))
        END DO

 10   IF (iAdd .GT. 0) THEN 
        write (kStdErr,*) 'whoops! wierd error parsing in caaMixFileInfo'
        CALL DoSTOP
        END IF

      IF (iType. LT. 0) THEN
        !read in list of integer followed by iNum reals
        read (caaTemp,*,err=5,end=10) iJ,(raStore(iI),iI=1,iNum)
        END IF

      IF (iType. GT. 0) THEN
        !read in list of iNum (integer,real) pairs
        read (caaTemp,*,err=5,end=10) (raStore2(iI),iI=1,2*iNum)
        DO iI=1,iNum
          iaStore(iI)=int(raStore2(iI*2-1))
          raStore(iI)=    raStore2(iI*2)
          END DO
        END IF

      iLines = iLinesStop+1

 6000 FORMAT(I3,' ',A80)
      RETURN
      END

c************************************************************************
c this subroutine reads stuff from caaM, and stores the resuls either in
c raStore alone or iaStore,raStore depending on value of iType
c this subroutine works very well. unfortunately, it uses kDumbFile which might
c be unwise as many possible kcarta.x might be trying to use the same file
c simultanelosuly, thus mucking up the all important MIXING TABLE
      SUBROUTINE readlistWORKS(caaM,iType,iNum,iLines,raStore,iaStore)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
      CHARACTER*80 kDumbFile

      kDumbFile='dumdum.dum'

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
