c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Resered

c this file deals with reading in the RTP file

c************************************************************************
c this subroutine deals with figuring out the wavenumbers
c very bloody simple, as Howard Motteler very nicely specifies these numbers
c for me in the RTP file!
      SUBROUTINE  IdentifyChannels(rf1,rf2,iRTP,caPFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c output
      REAL rf1,rf2           !the low and high end wavenumbers
c input
      CHARACTER*130 caPFname !the RTP file to peruse
      INTEGER iRTP           !which profile to read
      INTEGER inpath
c local variables for RTP file
      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

c other local variables
      integer i

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
c      write(kStdWarn,*) 'read open status = ', status

c      DO i = 1, iRTP
c        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
c        status = rtpread(rchan, prof)
c      END DO
      status = rtpclose(rchan)
c      write(kStdWarn,*) 'read close status = ', status
      kProfileUnitOpen = -1

      rf1 = head%vcmin
      rf2 = head%vcmax

      IF ((rf1 .LT. 0) .AND. (rf2 .LT. 0)) THEN
        write(kStdWarn,*) 'resetting head%vcmin from ',rf1,' to 605.0'
        rf1 = 605.0
        write(kStdWarn,*) 'resetting head%vcmax from ',rf2,' to 2830.0'
        rf2 = 2830.0
      ENDIF 

      IF (rf1 .LT. 0) THEN
        write(kStdErr,*) 'head%vcmin = ',rf1
        CALL DoStop
      END IF

      IF (rf2 .LT. 0) THEN
        write(kStdErr,*) 'head%vcmax = ',rf2
        CALL DoStop
      END IF

c      rf1 = 2105
c      rf2 = 2405
c      print *,'*********************************************************'
c      print *, 'Set rf1,rf2 =  ',rf1,rf2,' cm-1 for testing'
c      print *,'*********************************************************'

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format

c the kLAYERS format already differs from GENLN2 format in that
c (1) instead of velocity, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

c now we have an additional format to deal with, which should be MUCH simpler
      SUBROUTINE READRTP(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iRTP,iProfileLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      INTEGER inpath  ! added ESM inpath
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname
c for 100 layer clouds
      INTEGER iaCld100Read(3)
      REAL raaCld100Amt(kProfLayer,3)

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer iPtype

      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      fname(1:80) = caPFName(1:80)

      DO status = 1,3
        iaCld100Read(status) = -1
      END DO

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      iPtype = head%ptype
      write(kStdWarn,*) 'head%ptype = ',iPtype
      write(kStdWarn,*) 'head%ngas  = ',head%ngas
      status = rtpclose(rchan)

      IF (iPtype .EQ. 0) THEN
        write(kStdErr,*) 'KCARTA expects layers or pseudolevels profile'
        write(kStdErr,*) 'h.ptype == 1 or 2'
        CALL DOStop
      ELSEIF (iPtype .EQ. 1) THEN
        !! layers profile
        write(kStdWarn,*) 'Expecting ',iNumGases,' gases in rtp profile'
        IF (head%ngas .GE. iNumGases) THEN
          !read in rtp profile; hope all gas profiles are there
          CALL READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iaCld100Read,raaCld100Amt,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)
        ELSEIF (head%ngas .LT. iNumGases) THEN
          !read in rtp profile; augment profiles using US Std
          CALL READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iaCld100Read,raaCld100Amt,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)
        END IF
      ELSEIF (iPtype .EQ. 2) THEN
        !! pseudolevels profile
        CALL READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 1
c ie these are the AIRS layers
c EXPECTS to find ALL necessary gases here

c the kLAYERS format already differs from GENLN2 format in that
c (1) instead of velocity, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

c now we have an additional format to deal with, which should be MUCH simpler
      SUBROUTINE READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iaCld100Read,raaCld100Amt,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iRTP,iProfileLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

      REAL raaHeight(kProfLayer,kGasStore),MGC,delta1
      REAL raH1(kProfLayer),raP1(kProfLayer+1)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iNpath,iaNpathcounter(kMaxProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
       
      INTEGER iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
      CHARACTER*130 ca1,ca2,caTemp

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,iPtype
      REAL raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname
      logical isfinite

      MGC = 8.314674269981136  

      DO i = 1,kProfLayer
        pProf(i) = 0.0
      END DO

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      iPtype = head%ptype
      write(kStdWarn,*) 'head%ptype = ',iPtype

      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
      write(kStdWarn,*) 'read open status = ', status

      DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status .EQ. -1) THEN
          write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdWarn,*) fname

          write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdErr,*) fname
          CALL DoStop
        END IF
      END DO

      write (kStdWarn,*) 'success : read in RTP profile ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status

      kProfileUnitOpen = -1

      IF (prof%plevs(1) .lt. prof%plevs(prof%nlevs)) THEN
        !layers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof%plevs(prof%nlevs)
        kRTP_pTop = prof%plevs(1)
        kRTPBot   = kProfLayer - (prof%nlevs-1) + 1
        kRTPTop   = kProfLayer
      ELSE
        !layers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof%plevs(prof%nlevs)
        kRTP_pBot  = prof%plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof%nlevs-1
      END IF

      iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
      iProfileLayers = iL1
      iGasInRTPFile = head%ngas              !!! number of gases

      IF (prof%nlevs .GT. kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
      END IF
 
      write(kStdWarn,*) 'Reading profile from RTP file... '
      write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
      write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
      write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

      !!!now check if this agrees with iL1,iGasInRTPFile above
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
      END IF
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
      END IF

      DO i = 1,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof%palts(i)                     !!!!in meters
        raPressLevels(j) = prof%plevs(i)                !!!!in mb
      END DO

      DO i = 1,prof%nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
      END DO

      !check that spres lies withn plevs(nlevs) and plevs(nlevs-1)
      IF ((prof%plevs(prof%nlevs) .GT. prof%spres) .AND. 
     $    (prof%plevs(prof%nlevs-1) .GT. prof%spres)) THEN
        write(kStdErr,*) 'p.nlevs | p.plevs(p.nlevs) p.spres p.plevs(p.nlevs-1)'
        write(kStdErr,*) prof%nlevs,prof%plevs(prof%nlevs),prof%spres,prof%plevs(prof%nlevs-1)
        write(kStdErr,*) 'spres not between p.plevs(nlevs) and p.plevs(nlevs-1)'
        print *,'i       raP(i)          raPavg(i)        raP(i+1)    spres'
        print *,'----------------------------------------------------------'
        DO i = 1,prof%nlevs-1
          print *,i,raPressLevels(i),pProf(i),raPressLevels(i+1),prof%spres
        END DO

        CALL DoStop
      END IF

      IF (iDownWard .EQ. -1) THEN
        !!!add on dummy stuff
        !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof%nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
      ELSE
        !!!add on dummy stuff
        !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
      END IF

      DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
        IF (raThickness(i) .LE. 100.00) THEN
          write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          CALL DoStop
        END IF
      END DO

c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn = 0

c this variable keeps track of how many gases should be read in
      iNeed2Read = iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) .EQ. 1) THEN
          iNeed2Read = iNeed2Read-1
        END IF
      END DO

c this keeps track of the GasID used for the temperature .. hopefully water
c this keeps track of if we need to read in more gas profiles
      iTempFound        = -1
      iNeedMoreProfiles = -1

      caWord = '*PTHFIL'
      iErr   = -1

      iNumberOfGasesRead = 0
c set all individual gas paths to zero        
      DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
      END DO

c set this temp varaiable
      DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
      END DO

c set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
c so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
c iaGases(1) = iaGases(2) = iaGases(22)=iaGases(51)=1
c so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
      DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
      END DO
      iErr = 1
      DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr) = iNpath
          iErr = iErr+1
        END IF
      END DO

c now loop iNpath/iNumGases  times for each gas in the user supplied profile
c make sure you only do things for gases 1- 63
      DO iG = 1, iGasInRTPFile
        iIDGas = head%glist(iG)
        IF ((iIDGas .GT. kGasXsecHi) .AND. (iIDGAS. LT. kNewCloudLo)) THEN
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
          write(kStdWarn,*) 'this is something we may ignore for "gas" profs'
          write(kStdWarn,*) 'as it looks like continuum (101,102)'
        ELSEIF (iIDGas .LE. kGasXsecHi) THEN
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
          !!! first fill things out with stuff from the RTP file
          DO i = 1, prof%nlevs - 1
            j = iFindJ(kProfLayer,I,iDownWard)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

            rAmt = prof%gamnt(i,iG) / kAvog
            IF (isfinite(rAmt) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
              CALL dostop
            END IF              
            rT   = prof%ptemp(i)
            IF (isfinite(rT) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
              CALL dostop
            END IF              

            plays(i) = (prof%plevs(i)-prof%plevs(i+1))/
     $                 log(prof%plevs(i)/prof%plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            IF (iDownWard .EQ. -1) THEN
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "j"!!! 
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0) 
            ELSE 
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "i"!!!
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0) 
            END IF
            rH   = prof%palts(i)
            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                              iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                !write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH       !lalready in meters
                iaWhichGasRead(iIDgas) = 1
              END IF
            END IF
          END DO              !DO i = 1, prof%nlevs - 1 for klayers info

          !!! then fill bottom of atm with zeros for gas amt, partial pressure
          DO i = prof%nlevs, kProfLayer
            j = iFindJ(kProfLayer,I,iDownWard)
            iIDGas = head%glist(iG)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
            IF (iDownWard .EQ. -1) THEN
              delta1 = (300-prof%ptemp(prof%nlevs-1))/(1-(kProfLayer-prof%nlevs))
              rT   = 300.0  + delta1*j
              rT = 300.0
            ELSE
              delta1 = (200-prof%ptemp(prof%nlevs-1))/(kProfLayer-prof%nlevs)
              rT   = prof%ptemp(prof%nlevs-1) + delta1*j
              rT   = 300.0
            END IF
            rAmt = 0.0
            rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
            rPP  = 0.0
            rH   = raHeight(j)
            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                            iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                write(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,
     $ 'gindx,layer ',iGasIndex,i
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH
                iaWhichGasRead(iIDgas) = 1
              END IF
            END IF
          END DO    !DO i = prof%nlevs, kProfLayer for zeros
          CALL ContinuumFlag(iIDGas,iaCont)
    
          iFileGasesReadIn = iFileGasesReadIn+1
          WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
c this checks to see if we have read the profiles for all iNumGases required
c note that the gases read in MUST have been entered in GASFIL or XSCFIL 
c to count toward the tally ...
          IF (iaGases(iIDgas) .GT. 0) THEN
            iNumberOfGasesRead = iNumberOfGasesRead+1
            iaAlreadyIn(iNumberOfGasesRead) = iIDGas      
          ELSE
            write(kStdWarn,6000) iIDgas
          END IF

        ELSEIF ((iIDGAS. GE. kNewCloudLo) .AND. (iIDGAS. LE. kNewCloudHi)) THEN
          k100layerCloud = +1
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
          !!! first fill things out with stuff from the RTP file
          iNpathCounterJunk = 0
          DO i = 1, prof.nlevs - 1
            j = iFindJ(kProfLayer,I,iDownWard)
            iNpathCounterJunk = iNpathCounterJunk + 1

            rAmt = prof.gamnt(i,iG)
            IF (isfinite(rAmt) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
              CALL dostop
            END IF              
            rT   = prof.ptemp(i)
            IF (isfinite(rT) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
              CALL dostop
            END IF              

            plays(i) = (prof.plevs(i)-prof.plevs(i+1))/
     $                 log(prof.plevs(i)/prof.plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            IF (iDownWard .EQ. -1) THEN
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "j"!!! 
              rPP  = 0
            ELSE 
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "i"!!!
              rPP  = 0
            END IF
            rH   = prof.palts(i)

            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iNpathCounterJunk)
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            iGasIndex = iIDgas-kNewCloudLo+1
            raaCld100Amt(j,iGasIndex) = rAmt
            iaCld100Read(iGasIndex)   = 1
          END DO              !DO i = 1, prof.nlevs - 1 for klayers info

          !!! then fill bottom of atm with zeros for gas amt, partial pressure
          DO i = prof.nlevs, kProfLayer
            j = iFindJ(kProfLayer,I,iDownWard)
            iIDGas = head.glist(iG)
            iNpathCounterJunk = iNpathCounterJunk + 1
            IF (iDownWard .EQ. -1) THEN
              delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
              rT   = 300.0  + delta1*j
              rT = 300.0
            ELSE
              delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
              rT   = prof.ptemp(prof.nlevs-1) + delta1*j
              rT   = 300.0
            END IF
            rAmt = 0.0
            rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
            rPP  = 0.0
            rH   = raHeight(j)
            raaCld100Amt(j,iGasIndex)       = rAmt
            iaCld100Read(iGasIndex)      = 1
          END DO    !DO i = prof.nlevs, kProfLayer for zeros

        END IF      !if iGasID <= 63
      END DO

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,
     $          iaInputOrder,iNumGases,
     $          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

c first check to see if all required gases found in the user supplied profile
      IF (iNumberOfGasesRead .LT. iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'your profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'your profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        CALL DoStop
      END IF

 4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3) 
 6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
 5030 FORMAT(A130)
 4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

c now set raLayerHeight
      DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
      END DO

!change layer thickness to meters, because this is what rad_* routines need
      DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
      END DO
      DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
      END DO

      i = prof%nlevs - 1     !!!!!!number of layers in RTP file
      i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',
     $ raP1(i),raP1(i+1)
      write (kStdWarn,*) 'Highest layer : press levels (mb) = ',
     $ raP1(kProfLayer),raP1(kProfLayer+1)
      write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
      write (kStdWarn,*) '2 Highest layers thickness (km) = ',
     $  raH1(kProfLayer-1),raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
c see ../INCLUDE/KCARTA_database.param for the kCARTA database definitions
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 1
c ie these are the AIRS layers
c EXPECTS to find MOST gases here, but then goes off to augment the 
c remaining gases from US Std

c the kLAYERS format already differs from GENLN2 format in that
c (1) instead of velocity, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

c now we have an additional format to deal with, which should be MUCH simpler
      SUBROUTINE READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iaCld100Read,raaCld100Amt,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iRTP,iProfileLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

      REAL raaHeight(kProfLayer,kGasStore),MGC,delta1
      REAL raH1(kProfLayer),raP1(kProfLayer+1)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iNpath,iaNpathcounter(kMaxProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
       
      INTEGER iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
      CHARACTER*130 ca1,ca2,caTemp

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,iPtype
      REAL raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

      INTEGER iNpathCounterJunk,iaCld100Read(3)
      REAL raaCld100Amt(kProfLayer,3)

      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname
      logical isfinite

      MGC = 8.314674269981136  

      DO i = 1,kProfLayer
        pProf(i) = 0.0
      END DO

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      iPtype = head%ptype
      write(kStdWarn,*) 'head%ptype = ',iPtype

      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
      write(kStdWarn,*)  'read open status = ', status

      DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status .EQ. -1) THEN
          write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdWarn,*) fname

          write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdErr,*) fname
          CALL DoStop
        END IF
      END DO

      write (kStdWarn,*) 'success : read in RTP profile ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status

      kProfileUnitOpen = -1

      IF (prof%plevs(1) .lt. prof%plevs(prof%nlevs)) THEN
        !layers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof%plevs(prof%nlevs)
        kRTP_pTop = prof%plevs(1)
        kRTPBot   = kProfLayer - (prof%nlevs-1) + 1
        kRTPTop   = kProfLayer
      ELSE
        !layers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof%plevs(prof%nlevs)
        kRTP_pBot  = prof%plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof%nlevs-1
      END IF

      iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
      iProfileLayers = iL1
      iGasInRTPFile = head%ngas              !!! number of gases

      IF (prof%nlevs .GT. kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
      END IF
 
      write(kStdWarn,*) 'Reading profile from RTP file... '
      write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
      write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
      write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

      !!!now check if this agrees with iL1,iGasInRTPFile above
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
      END IF
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
      END IF

      DO i = 1,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof%palts(i)                     !!!!in meters
        raPressLevels(j) = prof%plevs(i)                !!!!in mb
      END DO

      DO i = 1,prof%nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
      END DO

      IF (iDownWard .EQ. -1) THEN
        !!!add on dummy stuff
        !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof%nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
      ELSE
        !!!add on dummy stuff
        !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
      END IF

      DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
        IF (raThickness(i) .LE. 100.00) THEN
          write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          CALL DoStop
        END IF
      END DO

c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn=0

c this variable keeps track of how many gases should be read in
      iNeed2Read = iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) .EQ. 1) THEN
          iNeed2Read = iNeed2Read-1
        END IF
      END DO

c this keeps track of the GasID used for the temperature .. hopefully water
c this keeps track of if we need to read in more gas profiles
      iTempFound        = -1
      iNeedMoreProfiles = -1

      caWord = '*PTHFIL'
      iErr   = -1

      iNumberOfGasesRead = 0
c set all individual gas paths to zero        
      DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
      END DO

c set this temp varaiable
      DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
      END DO

c set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
c so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
c iaGases(1) = iaGases(2) = iaGases(22) = iaGases(51) = 1
c so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
      DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
      END DO
      iErr = 1
      DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr) = iNpath
          iErr = iErr+1
        END IF
      END DO

c now loop iNpath/iNumGases  times for each gas in the user supplied profile
c make sure you only do things for gases 1- 63
      DO iG = 1, iGasInRTPFile
        iIDGas = head%glist(iG)
        IF ((iIDGas .GT. kGasXsecHi) .AND. (iIDGAS. LT. kNewCloudLo)) THEN
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
          write(kStdWarn,*) 'this is something we may ignore for "gas" profs'
          write(kStdWarn,*) 'as it looks like continuum (101,102)'
        ELSEIF (iIDGas .LE. kGasXsecHi) THEN
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
          !!! first fill things out with stuff from the RTP file
          DO i = 1, prof%nlevs - 1
            j = iFindJ(kProfLayer,I,iDownWard)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

            rAmt = prof%gamnt(i,iG) / kAvog
            IF (isfinite(rAmt) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
              CALL dostop
            END IF              
            rT   = prof%ptemp(i)
            IF (isfinite(rT) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
              CALL dostop
            END IF              

            plays(i) = (prof%plevs(i)-prof%plevs(i+1))/
     $                 log(prof%plevs(i)/prof%plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            IF (iDownWard .EQ. -1) THEN
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "j"!!! 
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0) 
            ELSE 
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "i"!!!
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0) 
            END IF
            rH   = prof%palts(i)

            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                              iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                !write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH       !lalready in meters
                iaWhichGasRead(iIDgas) = 1
              END IF
            END IF
          END DO              !DO i = 1, prof%nlevs - 1 for klayers info

          !!! then fill bottom of atm with zeros for gas amt, partial pressure
          DO i = prof%nlevs, kProfLayer
            j = iFindJ(kProfLayer,I,iDownWard)
            iIDGas = head%glist(iG)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
            IF (iDownWard .EQ. -1) THEN
              delta1 = (300-prof%ptemp(prof%nlevs-1))/(1-(kProfLayer-prof%nlevs))
              rT   = 300.0  + delta1*j
              rT = 300.0
            ELSE
              delta1 = (200-prof%ptemp(prof%nlevs-1))/(kProfLayer-prof%nlevs)
              rT   = prof%ptemp(prof%nlevs-1) + delta1*j
              rT   = 300.0
            END IF
            rAmt = 0.0
            rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
            rPP  = 0.0
            rH   = raHeight(j)
            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                            iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                write(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,
     $ 'gindx,layer ',iGasIndex,i
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH
                iaWhichGasRead(iIDgas) = 1
              END IF
            END IF
          END DO    !DO i = prof%nlevs, kProfLayer for zeros
          CALL ContinuumFlag(iIDGas,iaCont)
    
          iFileGasesReadIn = iFileGasesReadIn+1
          WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
c this checks to see if we have read the profiles for all iNumGases required
c note that the gases read in MUST have been entered in GASFIL or XSCFIL 
c to count toward the tally ...
          IF (iaGases(iIDgas) .GT. 0) THEN
            iNumberOfGasesRead = iNumberOfGasesRead+1
            iaAlreadyIn(iNumberOfGasesRead) = iIDGas      
          ELSE
            write(kStdWarn,6000) iIDgas
          END IF

        ELSEIF ((iIDGAS. GE. kNewCloudLo) .AND. (iIDGAS. LE. kNewCloudHi)) THEN
          k100layerCloud = +1
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
          !!! first fill things out with stuff from the RTP file
          iNpathCounterJunk = 0
          DO i = 1, prof.nlevs - 1
            j = iFindJ(kProfLayer,I,iDownWard)
            iNpathCounterJunk = iNpathCounterJunk + 1

            rAmt = prof.gamnt(i,iG)
            IF (isfinite(rAmt) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
              CALL dostop
            END IF              
            rT   = prof.ptemp(i)
            IF (isfinite(rT) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
              CALL dostop
            END IF              

            plays(i) = (prof.plevs(i)-prof.plevs(i+1))/
     $                 log(prof.plevs(i)/prof.plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            IF (iDownWard .EQ. -1) THEN
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "j"!!! 
              rPP  = 0
            ELSE 
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "i"!!!
              rPP  = 0
            END IF
            rH   = prof.palts(i)

            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iNpathCounterJunk)
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            iGasIndex = iIDgas-kNewCloudLo+1
            raaCld100Amt(j,iGasIndex) = rAmt
            iaCld100Read(iGasIndex)   = 1
          END DO              !DO i = 1, prof.nlevs - 1 for klayers info

          !!! then fill bottom of atm with zeros for gas amt, partial pressure
          DO i = prof.nlevs, kProfLayer
            j = iFindJ(kProfLayer,I,iDownWard)
            iIDGas = head.glist(iG)
            iNpathCounterJunk = iNpathCounterJunk + 1
            IF (iDownWard .EQ. -1) THEN
              delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
              rT   = 300.0  + delta1*j
              rT = 300.0
            ELSE
              delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
              rT   = prof.ptemp(prof.nlevs-1) + delta1*j
              rT   = 300.0
            END IF
            rAmt = 0.0
            rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
            rPP  = 0.0
            rH   = raHeight(j)
            raaCld100Amt(j,iGasIndex)       = rAmt
            iaCld100Read(iGasIndex)      = 1
          END DO    !DO i = prof.nlevs, kProfLayer for zeros

        END IF      !if iGasID <= 63
      END DO

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,
     $          iaInputOrder,iNumGases,
     $          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

c first check to see if all required gases found in the user supplied profile
      IF (iNumberOfGasesRead .LT. iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'head%ptype = 1 profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'head%ptype = 1 profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        CALL AddOnAFGLProfile(kAFGLProf,
     $         iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $         raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)
      END IF

 4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3) 
 6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
 5030 FORMAT(A130)
 4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

c now set raLayerHeight
      DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
      END DO

!change layer thickness to meters, because this is what rad_* routines need
      DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
      END DO
      DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
      END DO

      i = prof%nlevs - 1     !!!!!!number of layers in RTP file
      i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',
     $ raP1(i),raP1(i+1)
      write (kStdWarn,*) 'Highest layer : press levels (mb) = ',
     $ raP1(kProfLayer),raP1(kProfLayer+1)
      write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
      write (kStdWarn,*) '2 Highest layers thickness (km) = ',
     $  raH1(kProfLayer-1),raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
c see ../INCLUDE/KCARTA_database.param for the kCARTA database definitions
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 2
c ie these are the AIRS pseudolayers

c the kLAYERS format already differs from GENLN2 format in that
c (1) instead of velocity, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

c now we have an additional format to deal with, which should be MUCH simpler
      SUBROUTINE READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'
      include '../INCLUDE/airslevelheights.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iRTP,iProfileLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

      REAL raaHeight(kProfLayer,kGasStore),MGC,delta1
      REAL raH1(kProfLayer),raP1(kProfLayer+1)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iNpath,iaNpathcounter(kProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
       
      INTEGER iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
      CHARACTER*130 ca1,ca2,caTemp

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,iPtype,LBOT
      REAL raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)
      REAL raActualLayTemps(kProfLayer),TSURFA

      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname
      logical isfinite

      MGC = 8.314674269981136  

      DO i = 1,kProfLayer
        pProf(i) = 0.0
      END DO

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      iPtype = head%ptype
      write(kStdWarn,*) 'head%ptype = ',iPtype

      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
      write(kStdWarn,*)  'read open status = ', status

      DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status .EQ. -1) THEN
          write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdWarn,*) fname

          write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdErr,*) fname
          CALL DoStop
        END IF
      END DO

      write (kStdWarn,*) 'success : read in RTP profile ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status  =  ', status

      kProfileUnitOpen = -1

      prof%nlevs = prof%nlevs + 1   !!this really was number of LAYERS

      IF (prof%plevs(1) .lt. prof%plevs(prof%nlevs)) THEN
        !!reset prof%plevs so it has ALL the AIRS levels(1:101), rather than
        !!AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
        DO i = 1,kProfLayer+1
          prof%plevs(i) = PLEV_KCARTADATABASE_AIRS(kProfLayer+1-i+1)
          prof%palts(i) = DATABASELEVHEIGHTS(kProfLayer+1-i+1)*1000
        END DO
        DO i = 1,kProfLayer
          plays(i) = PAVG_KCARTADATABASE_AIRS(kProfLayer-i+1)
        END DO
        !layers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof%plevs(prof%nlevs)
        kRTP_pTop = prof%plevs(1)
        kRTPBot   = kProfLayer - (prof%nlevs-1) + 1
        kRTPTop   = kProfLayer
      ELSE
        !!reset prof%plevs so it has ALL the AIRS levels(1:101), rather than
        !!AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
        DO i = 1,kProfLayer+1
          prof%plevs(i) = PLEV_KCARTADATABASE_AIRS(i)
          prof%palts(i) = DATABASELEVHEIGHTS(i)*1000
        END DO
        DO i = 1,kProfLayer
          plays(i) = PAVG_KCARTADATABASE_AIRS(i)
        END DO
        !layers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof%plevs(prof%nlevs)
        kRTP_pBot  = prof%plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof%nlevs-1
      END IF

      iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
      iProfileLayers = iL1
      iGasInRTPFile = head%ngas              !!! number of gases

      IF (prof%nlevs .GT. kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
      END IF
 
      write(kStdWarn,*) 'Reading profile from RTP file... '
      write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
      write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
      write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

      !!!now check if this agrees with iL1,iGasInRTPFile above
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
      END IF
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
      END IF

      DO i = 1,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof%palts(i)                     !!!!in meters
        raPressLevels(j) = prof%plevs(i)                !!!!in mb
      END DO

      DO i = 1,prof%nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
      END DO

      IF (iDownWard .EQ. -1) THEN
        !!!add on dummy stuff
        !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof%nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
      ELSE
        !!!add on dummy stuff
        !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof%nlevs)
        DO i = prof%nlevs+1, kProfLayer + 1
          j = iFindJ(kProfLayer+1,I,iDownWard)
          raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
      END IF

      DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
        IF (raThickness(i) .LE. 100.00) THEN
          write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
          CALL DoStop
        END IF
      END DO
     
c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn = 0

c this variable keeps track of how many gases should be read in
      iNeed2Read = iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) .EQ. 1) THEN
          iNeed2Read = iNeed2Read-1
        END IF
      END DO

c this keeps track of the GasID used for the temperature .. hopefully water
c this keeps track of if we need to read in more gas profiles
      iTempFound        = -1
      iNeedMoreProfiles = -1

      caWord = '*PTHFIL'
      iErr   = -1

      iNumberOfGasesRead = 0
c set all individual gas paths to zero        
      DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
      END DO

c set this temp varaiable
      DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
      END DO

c set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
c so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then 
c         iaGases(1) = iaGases(2) = iaGases(22) = iaGases(51) = 1
c so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
      DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
      END DO
      iErr = 1
      DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr) = iNpath
          iErr = iErr+1
        END IF
      END DO

c**********
      !now map the pseudolevel temps to the layer temps
      !see /asl/packages/sartaV105/Src/mean_t.f
      write(kStdWarn,*) 'replacing pseudolevel temps with layer temps'
      LBOT = prof%nlevs-1
      !Do top layer (special case) 
      raActualLayTemps(1) = prof%ptemp(1) 
      !Loop down over the layers 
      DO i = 2,LBOT-1 
        raActualLayTemps(i) = 0.5*( prof%ptemp(i-1) + prof%ptemp(i) ) 
      ENDDO 
      ! Interpolate to get air temperature at the surface 
      TSURFA = prof%ptemp(LBOT-1) + ( prof%ptemp(LBOT) - prof%ptemp(LBOT-1) )* 
     $    (prof%spres - plays(LBOT))/(plays(LBOT+1)-plays(LBOT))
      !Do bottom layer (special case) 
      raActualLayTemps(LBOT) = 0.5*( prof%ptemp(LBOT-1) + TSURFA) 
c**********

c      DO i = 1, prof%nlevs
c        print *,i,prof%plevs(i),prof%ptemp(i),raActualLayTemps(i),
c     $          prof%gamnt(i,1),prof%gamnt(i,3)
c      END DO

c now loop iNpath/iNumGases  times for each gas in the user supplied profile
c make sure you only do things for gases 1- 63
      DO iG = 1, iGasInRTPFile
        iIDGas = head%glist(iG)
        IF (iIDGas .GT. kGasXsecHi) THEN
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
          write(kStdWarn,*) 'this is something we will ignore for "gas" profs'
          write(kStdWarn,*) 'either cloud (201,202,203) or cont (101,102)'
        ELSE
          write(kStdWarn,*) ' ---------------------------------------------'
          write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile
          !!! first fill things out with stuff from the RTP file
          DO i = 1, prof%nlevs - 1
            j = iFindJ(kProfLayer,I,iDownWard)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

            rAmt = prof%gamnt(i,iG) / kAvog
            IF (isfinite(rAmt) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
              CALL dostop
            END IF              

            rT   = prof%ptemp(i)
            rT   = raActualLayTemps(i)         !use this instead of prof%ptemp
            IF (isfinite(rT) .EQ. .false.) THEN
              write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
              CALL dostop
            END IF              

            plays(i) = (prof%plevs(i)-prof%plevs(i+1))/
     $                 log(prof%plevs(i)/prof%plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            IF (iDownWard .EQ. -1) THEN
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "j"!!! 
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0) 
            ELSE 
              !!! this automatically puts partial pressure in ATM, assuming 
              !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
              !!!note "i"!!!
              rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0) 
            END IF
            rH   = prof%palts(i)
            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                              iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                !write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH       !lalready in meters
                iaWhichGasRead(iIDgas)    = 1
              END IF
            END IF
          END DO              !DO i = 1, prof%nlevs - 1 for klayers info

          !!! then fill bottom of atm with zeros for gas amt, partial pressure
          DO i = prof%nlevs, kProfLayer
            j = iFindJ(kProfLayer,I,iDownWard)
            iIDGas = head%glist(iG)
            iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
            IF (iDownWard .EQ. -1) THEN
              delta1=(300-prof%ptemp(prof%nlevs-1))/(1-(kProfLayer-prof%nlevs))
              rT   = 300.0  + delta1*j
              rT = 300.0
            ELSE
              delta1 = (200-prof%ptemp(prof%nlevs-1))/(kProfLayer-prof%nlevs)
              rT   = prof%ptemp(prof%nlevs-1) + delta1*j
              rT   = 300.0
            END IF
            rAmt = 0.0
            rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
            rPP  = 0.0
            rH   = raHeight(j)
            !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
            CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
            IF (iaGases(iIDgas) .GT. 0) THEN
              Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                            iFound,iGasIndex)
              IF (iFound .GT. 0) THEN 
                write(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,
     $ 'gindx,layer ',iGasIndex,i
                raaAmt(j,iGasIndex)       = rAmt
                raaTemp(j,iGasIndex)      = rT
                raaPress(j,iGasIndex)     = rP
                raaPartPress(j,iGasIndex) = rPP
                raaHeight(j,iGasIndex)    = rH
                iaWhichGasRead(iIDgas)    = 1
              END IF
            END IF
          END DO    !DO i = prof%nlevs, kProfLayer for zeros
          CALL ContinuumFlag(iIDGas,iaCont)
    
          iFileGasesReadIn = iFileGasesReadIn+1
          WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
c this checks to see if we have read the profiles for all iNumGases required
c note that the gases read in MUST have been entered in GASFIL or XSCFIL 
c to count toward the tally ...
          IF (iaGases(iIDgas) .GT. 0) THEN
            iNumberOfGasesRead = iNumberOfGasesRead+1
            iaAlreadyIn(iNumberOfGasesRead) = iIDGas      
          ELSE
            write(kStdWarn,6000) iIDgas
          END IF
        END IF      !if iGasID <= 63
      END DO

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,
     $          iaInputOrder,iNumGases,
     $          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

c first check to see if all required gases found in the user supplied profile
      IF (iNumberOfGasesRead .LT. iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'head%ptype = 2 profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'head%ptype = 2 profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        CALL AddOnAFGLProfile(kAFGLProf,
     $         iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $         raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)
      END IF

 4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3) 
 6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
 5030 FORMAT(A130)
 4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

c now set raLayerHeight
      DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
      END DO

!change layer thickness to meters, because this is what rad_* routines need
      DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
      END DO
      DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
      END DO

      i = prof%nlevs - 1     !!!!!!number of layers in RTP file
      i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',
     $ raP1(i),raP1(i+1)
      write (kStdWarn,*) 'Highest layer : press levels (mb) = ',
     $ raP1(kProfLayer),raP1(kProfLayer+1)
      write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
      write (kStdWarn,*) '2 Highest layers thickness (km) = ',
     $  raH1(kProfLayer-1),raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
c see ../INCLUDE/KCARTA_database.param for the kCARTA database definitions
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this subroutine adds on US STandard Profile gas amounts for those gases NOT
c in the RTP profile (using h.ptype = 2 or even 1)
      SUBROUTINE AddOnStandardProfile(
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

      implicit none

      include '../INCLUDE/kcarta.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      INTEGER iNumberofGasesRead
      INTEGER igasindex  ! added ESM
      REAL    raaHeight(kProfLayer,kGasStore)
      INTEGER iaWhichGasRead(kMaxGas),iNumGases,iaInputOrder(kMaxGas)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)

c local variables
      INTEGER iI,iJ,iaNeed(kMaxGas),iNewRead,iFound,iFoundX,iIDgas,iLay,iG0,iK
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer)  !!US Std layer ppress, amt
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer)  !!US Std layer press, temp
      CHARACTER*1 cY,cN

      cY = 'Y'
      cN = 'N'

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'read STD profiles for ',iNumberOfGasesRead, ' gases ...'

c this is the list of gases for which we need profiles
      iK = 0
      write(kStdWarn,*) '            count gasID  found?'
      write(kStdWarn,*) '-------------------------------'
      DO iI = 1,iNumGases
        IF (iaInputOrder(iI) .LT. 100) THEN
          IF (iaWhichGasRead(iI) .EQ. +1) THEN
             iK = iK + 1
c            write(kStdWarn,200) iI,iaInputOrder(iI),cY
            write(kStdWarn,200) iK,iI,cY
          ELSE
c            write(kStdWarn,200) iI,iaInputOrder(iI),cN
            write(kStdWarn,200) -1,iI,cN
          END IF
        END IF
      END DO


      !!! now do gasids 101,102,103
      iI = 1
      IF (iaInputOrder(iI) .EQ. 1) THEN
        IF (iaWhichGasRead(iI) .EQ. +1) THEN
          write(kStdWarn,200) iI,101,cY
          write(kStdWarn,200) iI,102,cY
          write(kStdWarn,200) iI,103,cY
        ELSE  
          write(kStdWarn,200) iI,101,cN
          write(kStdWarn,200) iI,102,cN
          write(kStdWarn,200) iI,103,cN
        END IF
      END IF
      write(kStdWarn,*) ' '

 200  FORMAT('yi, idgas = ',I5,I5,'  ',A1)

c this is the list of gases for which we have read in the profiles
c      DO iI = 1,kMaxGas
c        IF (iaWhichGasRead(iI) .EQ. +1) THEN
c          write(kStdWarn,*) 'RTP file had profile for gasID ',iI
c        END IF
c      END DO
c      write(kStdWarn,*) ' '

      iFound = -1
      DO iI = 1,kMaxGas
        IF (iaWhichGasRead(iI) .EQ. +1) THEN
          GOTO 10
        END IF
      END DO
 10   CONTINUE
      iG0 = iI     !!this gas filled out from RTP file, so use T(z),h(z)
      
c thus the difference between the lists is what we need US Std profiles for
      DO iI = 1,kMaxGas
        iaNeed(iI) = -1
      END DO

      iNewRead = 0
      DO iI = 1,iNumGases
        iFound = -1
        iIDgas = iaInputOrder(iI)
        IF (iaWhichGasRead(iIDgas) .EQ. +1) THEN
          iFound = +1
        END IF
        IF (iFound .LT. 0) THEN
          !! gas not found in rtp file, so need the US Std profile
          iNewRead = iNewRead + 1
          iaNeed(iNewRead) = iIDgas
          CALL get_us_std(iIDgas,raPX,raPPX,raTX,raQX)
          Call FindIndexPosition(iIDgas,iNumGases,iaInputOrder,
     $                           iFoundX,iGasIndex)
          IF (iFoundX .GT. 0) THEN 
            !write(kStdWarn,4321) iIDGas,iLay,rAmt,rT,rP,rPP
            DO iLay = 1,kProfLayer
              raaAmt(iLay,iGasIndex)       = raQX(iLay)
              raaTemp(iLay,iGasIndex)      = raaTemp(iLay,iG0)
              raaPress(iLay,iGasIndex)     = raaPress(iLay,iG0)
              raaPartPress(iLay,iGasIndex) = raPPX(iLay)
              raaHeight(iLay,iGasIndex)    = raaHeight(iLay,iG0)
            END DO
            iaWhichGasRead(iIDgas)    = 1
          ELSE 
            write (kStdErr,*) 'huh? FindIndexPosition failed for ',iIDgas
            CALL DoStop
          END IF
        END IF
      END DO

      write(kStdWarn,*) 'Before entering "AddOnStandardProfile" '
      write(kStdWarn,*) '  had read in profiles for ',iNumberofGasesRead
      write(kStdWarn,*) '  out of ',iNumGases
      write(kStdWarn,*) 'Read in ',iNewRead,' more in "AddOnStandardProfile"'
 
      IF (iNumberofGasesRead + iNewRead .NE. iNumGases) THEN
        write(kStdErr,*) 'need iNumberofGasesRead + iNewRead = iNumGases'
        CALL DoStop
      END IF

      RETURN
      END

c************************************************************************
c this subroutine reads in the US Std Profile name and returns the gas amount
      SUBROUTINE get_us_std(iIDgas,raPX,raPPX,raTX,raQX)

      implicit none
      include '../INCLUDE/kcarta.param'

c input var
      INTEGER iIDgas
c output var
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer) !!US Std layer ppress, amt
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer) !!US Std layer press, temp

c local vars
      INTEGER i1,i2,iLenX,iLen,iLay,iX, iioun, ierr  ! added ESM iioun, ierr
      CHARACTER*80 caFname0,caFname,cnameX
      CHARACTER c1,c2
      CHARACTER*2 c12
      CHARACTER*100 caLine

c      caFname0 = kUSStd
      caFname0 = kOrigRefPath

      iLen = 80
 100  CONTINUE
      IF (caFname0(iLen:iLen) .EQ. ' ') THEN
        iLen = iLen - 1
        GOTO 100
      END IF

      cnameX = 'us_std_gas_'
      cnameX = 'refgas'
      iLenX = 80
 200  CONTINUE
      IF (cnameX(iLenX:iLenX) .EQ. ' ') THEN
        iLenX = iLenX - 1
        GOTO 200
      END IF

      DO i1 = 1,iLenX
        caFname0(iLen+i1:iLen+i1) = cnameX(i1:i1)
      END DO
      iLen = iLen + iLenX

c      caFname0 = '/home/sergio/KCARTADATA/USSTD/us_std_gas_' 

      IF (iIDgas .LT. 10) THEN
        c1 = ' '
        i2 = iIDgas
        c2 = CHAR(i2+48)
        c12 = c2//c1
      ELSE
        i2 = iIDgas/10
        i1 = (iIDgas-10*i2)
        c1 = CHAR(i1+48)
        c2 = CHAR(i2+48)
        c12 = c2//c1
      END IF
      DO i1 = 1,80
        caFname(i1:i1) = ' '
      END DO
      DO i1 = 1,iLen
        caFname(i1:i1) = caFname0(i1:i1)
      END DO
      DO i1 = 1,2
        caFname(iLen+i1:iLen+i1) = c12(i1:i1)
      END DO
      write(kStdWarn,*) 'need to add on US Std for ',iIDgas

 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)  
      iIOUN = kTempUnit  
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='FORMATTED',  
     $    IOSTAT=IERR)  
      IF (IERR .NE. 0) THEN  
        WRITE(kStdErr,*) 'In subroutine get_us_std have file I/O error'  
        write(kStdErr,*) 'reference path = ',caFname0
        write(kStdErr,*) 'gasID = ',iIDGas
        WRITE(kStdErr,1010) IERR, caFname
        CALL DoSTOP  
      ENDIF  
 
ccc this was orig code
c      kTempUnitOpen = 1  
c      iLay = 0 
c 1020 CONTINUE 
c      iLay = iLay + 1 
c      READ(iIOUN,*,END=1030) iX,raPx(iLay),raPPx(iLay),raTx(iLay),raQx(iLay) 
c      GOTO 1020 
c 
c 1030  CONTINUE 
c      CLOSE(iIOUN)  
c      kTempUnitOpen = -1  

c this is new code, same as in subr ReadRefProf (in n_pth_mix.f)
      kTempUnitOpen = 1  
      iLay = 0
 20   READ(iIOUN,5020,END=199) caLine
 5020 FORMAT(A100)
      IF (caLine(1:1) .NE. '!') THEN
        iLay=iLay+1
        READ(caLine,*) iX,raPx(iLay),raPPx(iLay),raTx(iLay),raQx(iLay)
      ENDIF
      GOTO 20
 199  CLOSE(iIOUN)
      kTempUnitOpen = -1

      RETURN
      END

c************************************************************************
c this subroutine reads in the US Std Profile name and returns the gas amount
      SUBROUTINE getAFGL(iProfileNum,iIDgas,raPX,raPPX,raTX,raQX)

      implicit none
      include '../INCLUDE/kcarta.param'

c input var
      INTEGER iIDgas,iProfileNum
c output var
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer) !!US Std layer ppress, amt
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer) !!US Std layer press, temp

c local vars
      INTEGER i1,i2,iLenX,iLen,iLay,iX, iioun, ierr  ! added ESM iioun, ierr
      CHARACTER*80 caFname0,caFname,cnameX
      CHARACTER c1,c2
      CHARACTER*2 c12
      CHARACTER*100 caLine

      caFname0 = kOrigRefPath

      iLen = 80
 100  CONTINUE
      IF (caFname0(iLen:iLen) .EQ. ' ') THEN
        iLen = iLen - 1
        GOTO 100
      END IF
      caFName0(iLen:iLen) = '/'
      iLen = iLen+1

      IF (kAFGLProf .EQ. 1) caFName0(iLen:iLen) = '1'
      IF (kAFGLProf .EQ. 2) caFName0(iLen:iLen) = '2'
      IF (kAFGLProf .EQ. 3) caFName0(iLen:iLen) = '3'
      IF (kAFGLProf .EQ. 4) caFName0(iLen:iLen) = '4'
      IF (kAFGLProf .EQ. 5) caFName0(iLen:iLen) = '5'      
      IF (kAFGLProf .EQ. 6) caFName0(iLen:iLen) = '6'
      iLen = iLen+1

      caFName0(iLen:iLen) = '/'

      cnameX = 'us_std_gas_'
      cnameX = 'refgas'
      cnameX = 'afglgas'
      iLenX = 80
 200  CONTINUE
      IF (cnameX(iLenX:iLenX) .EQ. ' ') THEN
        iLenX = iLenX - 1
        GOTO 200
      END IF

      DO i1 = 1,iLenX
        caFname0(iLen+i1:iLen+i1) = cnameX(i1:i1)
      END DO
      iLen = iLen + iLenX

      IF (iIDgas .LT. 10) THEN
        c1 = ' '
        i2 = iIDgas
        c2 = CHAR(i2+48)
        c12 = c2//c1
      ELSE
        i2 = iIDgas/10
        i1 = (iIDgas-10*i2)
        c1 = CHAR(i1+48)
        c2 = CHAR(i2+48)
        c12 = c2//c1
      END IF
      DO i1 = 1,80
        caFname(i1:i1) = ' '
      END DO
      DO i1 = 1,iLen
        caFname(i1:i1) = caFname0(i1:i1)
      END DO
      DO i1 = 1,2
        caFname(iLen+i1:iLen+i1) = c12(i1:i1)
      END DO
      write(kStdWarn,*) 'need to add on AFGL for ',iIDgas

 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)  
      iIOUN = kTempUnit  
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='FORMATTED',  
     $    IOSTAT=IERR)  
      IF (IERR .NE. 0) THEN  
        WRITE(kStdErr,*) 'In subroutine getAFGL have file I/O error'  
        write(kStdErr,*) 'looking to fill in profile with = ',caFname0
        write(kStdErr,*) 'gasID = ',iIDGas
        WRITE(kStdErr,1010) IERR, caFname
        CALL DoSTOP  
      ENDIF  

c this is new code, same as in subr ReadRefProf (in n_pth_mix.f)
      kTempUnitOpen = 1  
      iLay = 0
 20   READ(iIOUN,5020,END=199) caLine
 5020 FORMAT(A100)
      IF (caLine(1:1) .NE. '!') THEN
        iLay=iLay+1
        READ(caLine,*) iX,raPx(iLay),raPPx(iLay),raTx(iLay),raQx(iLay)
      ENDIF
      GOTO 20
 199  CLOSE(iIOUN)
      kTempUnitOpen = -1

      RETURN
      END

c************************************************************************
c this subroutine deals with the 'RADNCE' keyword, but for new .rtp files
      SUBROUTINE radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,
     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $   raPressLevels,iProfileLayers,
     $   raFracTop,raFracBot,raaPrBdry,
     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $   raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,
     $   iakSolar,rakSolarAngle,rakSolarRefl,
     $   raSatAzimuth,raSolAzimuth,
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iaNumLayer,iaaRadLayer,raProfileTemp,
     $   cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP,
     $   raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c caSetEmissivity= array that gives name of emissivity files (if any) 
c raSetEmissivity= array that gives constant emissivity value (if set)
c iNpmix     = number of mixed paths read in from mixfile
c iaMPSetForRad = array telling which MP set to associate with which atm
c iNatm       = number of atmospheres
c raPressStart = start pressure for radiating atmos
c raPressStop  = stop pressure for radiating atmos
c raTSpace    = array containing background temperature for each atmosphere
c raTSurf    = array contianing surface temperature for each atmosphere
c raSatAngle = array containing satellite view angle for each atmosphere
c raSatHeight= array containing satellite height for each atmosphere
c iaNumLayer = array containing number of layers in each atmosphere
c iaaRadLayer= matrix containing list of layers in each atmosphere
c iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
c raaaSetEmissivity = array containing the wavenumber dependent emissivities
c raFracTop  = top fraction
c raFracBot  = bottom fraction
c raaPrBdry  = matrix that keeps start/stop pressures
c the next few only work for DOWNWARD LOOK instr
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c rakSolarRefl   =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c raProfileTemp = array containing CO2 gas profile temperature
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
c this is for the cloud, if any, that is associated with the atmosphere
c raPressLevels are the actual pressure levels from the KLAYERS file
      INTEGER iProfileLayers,ctype1,ctype2
      REAL raPressLevels(kProfLayer+1)

      REAL cfrac12,cFrac1,cFrac2,cngwat1,cngwat2,raCemis(kMaxClouds)
      REAL raCprtop(kMaxClouds), raCprbot(kMaxClouds)
      REAL raCngwat(kMaxClouds), raCpsize(kMaxClouds)
      INTEGER iaCtype(kMaxClouds),iMPSetForRadRTP
      INTEGER iNclouds_RTP     !!!tells how many clouds
      INTEGER iaNML_Ctype(kMaxClouds)

      CHARACTER*80 caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm) 
      INTEGER iaMPSetForRad(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm),iNpmix
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      INTEGER iRTP     !!!tells which profile info, radiance info, to read
      CHARACTER*130  caPFName !!!tells which profile 

c local variables
      CHARACTER*7 caWord
      INTEGER iNlay,iStart,iStop,iErr
      REAL rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight,rT
      INTEGER iDirection,iW,iInt
      INTEGER iC,iNumLinesRead
      REAL FindSurfaceTemp,rSize1,rSize2
      INTEGER iInstrType

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,upwell,iOKscanang,iOKzobs,iOKsatzen
      REAL raHeight(kProfLayer+1),raThickness(kProfLayer),pobs,pobs1,pTemp,rSURFaltitude
      REAL r1,rEms,rAngleX,rAngleY,saconv_sun,orig_saconv_sun
      INTEGER*4 i4CTYPE1,i4CTYPE2
 
      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname
      real rf1,rf2

      print *,'booboo A ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      fname(1:80) = caPFName(1:80)

      write(kStdWarn,*) 'Using RTP file to set atm info ....'
      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF

      print *,'booboo B ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      kProfileUnitOpen = +1
      !since we succesffuly read in profile, no need to do as many checks here
      DO i = 1, iRTP
        status = rtpread(rchan, prof)
      END DO
      status = rtpclose(rchan)
      kProfileUnitOpen = -1

      print *,'booboo c ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      rf1 = head%vcmin
      rf2 = head%vcmax

      caWord = '*RADNCE'
      iErr = -1

      iNatm = 1        !can only read ONE atmosphere per RTP profile

c now get the relevant info from rchan,prof
      iC = 1
      iW = iMPSetForRadRTP
      iaMPSetForRad(1) = iMPSetForRadRTP

      !!!assume the instrument is downlooking, from TOA
      pobs = 0.0
     
      rSURFaltitude = prof%salti
      
      rPressStart = prof%spres
      rPressStop  = 0.000            ! ----------> assume TOA

      !!!then go ahead and look at variables prof%pobs
      !!!note that variable pobs is reset only if prof%pobs > 0, else it
      !!!stays at a value of 0.0000 (TOA)
      pobs1 = prof%pobs
      IF (pobs1 .GT. 0) THEN
        IF (pobs1  .lt. raPressLevels(iProfileLayers+1)) THEN
          write(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
          write(kStdWarn,*) 'frm ',pobs1,' to ',raPressLevels(iProfileLayers+1)
          pobs1 = raPressLevels(iProfileLayers+1)
          pobs = pobs1
        END IF
        upwell = ((kProfLayer + 1) - prof%nlevs) + 1
        IF (pobs1  .gt. raPressLevels(upwell)) THEN
          write(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
          write(kStdWarn,*) 'from ',pobs1,' to ',raPressLevels(upwell)
          pobs1 = raPressLevels(upwell)
          pobs = pobs1
        END IF
      END IF
      pobs = pobs1

c testing
c      prof%satzen = -abs(prof%satzen)
c      prof%zobs   = -1000
c       prof%scanang = -abs(prof%scanang) * 1000

      !!!assume the instrument is downlooking, from TOA
      upwell = 1
      !!!then go ahead and look at variables prof%upwell
      IF (prof%upwell .EQ. 1) THEN
        !!radiation is travelling upwards so this is a downlook instrument
        upwell = prof%upwell
      ELSEIF (prof%upwell .EQ. 2) THEN
        !!radiation is travelling downwards so this is a uplook instrument
        upwell = prof%upwell
      ELSE
        write(kStdErr,*) 'need prof%upwell = 1 (downlook) or 2 (uplook)'
        write(kStdErr,*) 'prof%upwell = ',prof%upwell
        CALL DoStop
      END IF

      !now that we have rPressStart,rPressStop (defining pressure boundaries 
      !for the atm) and upwell (direction of radiation travel), check to see 
      !if things need to be reset
      IF (upwell .EQ. 1) THEN
        !need rPressStart > rPressStop
        rPressStart = prof%spres
        rPressStop  = pobs
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof%pobs,prof%upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a downlook instr, with Surf,OBS press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart .LE. rPressStop) THEN
          write(kStdErr,*) 'For downlook instr, need rPressStart > rPressStop'
          CALL DoStop
        END IF
      ELSEIF (upwell .EQ. 2) THEN
        !need rPressStart < rPressStop
        rPressStart = 0.0
        rPressStop  = prof%spres
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof%pobs,prof%upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a uplook instr, with TOA,Surf press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart .GE. rPressStop) THEN
          write(kStdErr,*) 'For uplook instr, need rPressStart < rPressStop'
          CALL DoStop
        END IF
      ELSE
        write(kStdErr,*) 'Need to set prof%upwell = 1 or 2'
        Call DOStop
      END IF

      rTbdy       = kTSpace          ! -------> assume deep space
      rTSurf      = prof%stemp

c      DO i = 1, prof%nlevs - 1
c        rT   = prof%ptemp(i)
c      END DO 
c      print *,rTSurf,rT,rT+1
c      rTSurf      = rT+1   

c!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c scott has this in sarta.f (for AIRS!!)
C      Convert SATZEN or SATANG to viewing angle
c       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
C         Convert zenith angle at surface to view angle at satellite
c          SVA=SACONV( SATZEN, SALT*1000 )/CONV
c       ELSE
C         Check if scan angle is valid
c          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
C            View angle should be within a few degrees of scan angle
c             SVA=ABS( SATANG )
c          ELSE
c             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
c 1030        FORMAT('Error! Profile',I5,
c     $          ': invalid angles for SATZEN ',1PE11.4,
c     $          ' and SATANG ',E11.4)
c             STOP
c        ENDIF
c     ENDIF
c       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
c       IF (SVA .GT. ANGMAX) THEN
C         Truncate angle if too big
c          WRITE(IOINFO,1040) IPROF, SVA
c 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
c     $       1PE11.4,' to 53 degrees')
c          SVA=ANGMAX
c     ENDIF
C      Convert SATZEN or SATANG to viewing angle
c       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
C         Convert zenith angle at surface to view angle at satellite
c          SVA=SACONV( SATZEN, SALT*1000 )/CONV
c       ELSE
cC         Check if scan angle is valid
c          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
C            View angle should be within a few degrees of scan angle
c             SVA=ABS( SATANG )
c          ELSE
c             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
c 1030        FORMAT('Error! Profile',I5,
c     $          ': invalid angles for SATZEN ',1PE11.4,
c     $          ' and SATANG ',E11.4)
c             STOP
c        ENDIF
c     ENDIF
c
c       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
c       IF (SVA .GT. ANGMAX) THEN
C         Truncate angle if too big
c          WRITE(IOINFO,1040) IPROF, SVA
c 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
c     $       1PE11.4,' to 53 degrees')
c          SVA=ANGMAX
c     ENDIF

c assume scanang, satzen, zobs make sense
      write(kStdWarn,*) ' '

      iOKscanang = +1
      IF (abs(prof%scanang) .gt. 89.99) THEN
        write(kStdWarn,*) 'whoops : prof%scanang = ',prof%scanang
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKscanang = -1
      END IF

      iOKsatzen = +1
      IF (abs(prof%satzen) .gt. 89.99) THEN
        write(kStdWarn,*) 'whoops : prof%satzen = ',prof%satzen
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKsatzen = -1
      END IF

      iOKzobs = +1
      IF ((prof%zobs .lt. 2.00*1000) .AND. (prof%upwell .EQ. 1)) THEN
        write(kStdWarn,*) 'whoops : prof%zobs = ',prof%zobs
        write(kStdWarn,*) '         zobs < 2 km in height!'
        write(kStdWarn,*) '         does not make sense for downlook instr'
        iOKzobs = -1
      END IF

      rAngle      = prof%scanang     ! -------> ignore satzen,satazi
      rHeight     = prof%zobs        ! -------> use raSatHeight(iC) in m (changed in July 2013)

      IF ((iOKscanang .EQ. -1) .AND. (iOKsatzen .EQ. -1) .AND. 
     $      (iOKzobs .EQ. -1)) THEN
        write(kStdWarn,*) 'scanang,satzen,zobs do not make sense!' 
        write(kStdErr,*) 'scanang,satzen,zobs do not make sense!' 
        CALL DOStop
      END IF

      !! Larrabee prefers to use satzen, so if it exists, use that!!!!!
c > Secondly, I was looking into Larrabee's request to get kCARTA to handle
c > p.satzen rather than p.scanang. It'll take some doing, but in Matlab,
c > suppose I have p.satzen. Do I use
c >     [zang]=saconv( p.satzen,705000);
c > to convert that into scanang, for kCARTA to use?
c Yes, except use prof%zobs rather than a hardcoded 705000.
      rAngleX     = prof%satzen      ! -------> rtp_interface originally
                                     !          ignored satzen,satazi
c so the conversion is  p.scanang = orig_saconv_sun( p.satzen,prof%zobs);  %% by Scott
c so the conversion is  p.scanang = saconv( p.satzen,prof%zobs);           %% by Sergio
      write(kStdWarn,*) ' '
      IF (prof%upwell .EQ. 1) THEN
        IF (rHeight .GT. 2.0) THEN
          rAngleY = ORIG_SACONV_SUN(rAngleX, rHeight)
          rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, rHeight/1000)
        ELSE
          rAngleY = ORIG_SACONV_SUN(rAngleX, 705.00)                !! default AIRS hgt, dangerous
          rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, 705.00) !! default AIRS hgt, dangerous
          rAngleY = -9999
        ENDIF
        write(kStdWarn,*) 'downlook instr : satellite hgt, view angle info : '
        write(kStdWarn,*) 'scanang, zobs(km), satzen, saconv(satzen,zobs) = ',
     $    rAngle,rHeight/1000,rAngleX,rAngleY
      END IF

      IF (prof%upwell .EQ. 2) THEN
        !! uplook instrument
        IF (iOKscanang .EQ. 1) THEN 
          write(kStdWarn,*) 'Uplook instr : use prof%scanang'
        ELSEIF (iOKscanang .EQ. -1) THEN 
          write(kStdWarn,*) 'Uplook instr : incorrect prof%scanang'
          write(kStdErr,*) 'Uplook instr : incorrect prof%scanang',rAngle
          CALL DoStop
        END IF
      ELSEIF (prof%upwell .EQ. 1) THEN
        !! downlook instr
        IF ((iOKscanang .EQ. 1) .AND. (iOKsatzen .EQ. 1) .AND. 
     $      (iOKzobs .EQ. 1)) THEN
          !! all angles seem reasonable; now check consistency between them
          IF (abs(abs(rAngle)-abs(rAngleY)) .LE. 1.0e-2) THEN
            write(kStdWarn,*) 'scanang,satzen,zobs present in rtp file'
            write(kStdWarn,*) 'scanang and saconv(satzen,zobs) agree'
            rAngle = rAngleY   !!! no need to do anything
          ELSEIF (abs(abs(rAngle)-abs(rAngleY)) .GT. 1.0e-2) THEN
            write(kStdWarn,*) 'scanang,satzen,zobs present in rtp file'
            write(kStdWarn,*) 'scanang and saconv(satzen,zobs) disagree'
            write(kSTdWarn,*) 'using satzen (AIRS preference!!!)'
            IF (prof%zobs .LT. 2000.0) THEN
              write(kStdErr,*) 'used 705 km as satellite height'
            ELSE
              write(kStdErr,*) 'used ',prof%zobs/1000 ,' km as satellite hght'
            END IF
            rAngle = rAngleY   !!! replace input scanang with that derived
                               !!! from satzen,zobs
          END IF
        ELSEIF ((iOKscanang .EQ. 1) .AND. 
     $          ((iOKsatzen .EQ. -1) .AND. (iOKzobs .EQ. +1))) THEN
          write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
          rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang .EQ. 1) .AND. 
     $          ((iOKsatzen .EQ. +1) .AND. (iOKzobs .EQ. -1))) THEN
          write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
          rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang .EQ. 1) .AND. 
     $          ((iOKsatzen .EQ. -1) .AND. (iOKzobs .EQ. -1))) THEN
          write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
          rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang .EQ. -1) .AND. 
     $          ((iOKsatzen .EQ. +1) .AND. (iOKzobs .EQ. +1))) THEN
          !! satzen or zobs or both make sense, scanang is incorrect
          rAngle = rAngleY   !!! no need to do anything
          write(kStdWarn,*) 'scanang does not make sense, but satzen,zobs do'
          write(kSTdWarn,*) 'using satzen to derive scanang (AIRS preference!!)'
        END IF
      END IF

      write(kStdWarn,*) 'using kCARTA scanang = ',rAngle
      write(kStdWarn,*) ' '
c!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (rTbdy .GT. 3.0) THEN
        write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
        CALL DoStop
      END IF

      CALL StartStopMP(iW,rPressStart,rPressStop,iC,
     $                 raPressLevels,iProfileLayers,
     $                 raFracTop,raFracBot,raaPrBdry,iStart,iStop)
      raPressStart(1) = raaPrBdry(1,1)
      raPressStop(1)  = raaPrBdry(1,2)

c figure out if the start/stop MixedPath numbers are legitimate
      IF ((iStart .GT. iNpmix).OR.(iStart .LT. 1) .OR.
     $    (iStop .GT. iNpmix).OR.(iStop.LT. 1)) THEN
          write(kStdErr,*)'Error while setting Start/Stop Mixed Path '
          write(kStdErr,*)'numbers for atmosphere # ',iC
          write(kStdErr,*)'Must be between 1 and ',iNpmix
          write(kStdErr,*)'Instead they are ',iStart,iStop
          CALL DoSTOP
        END IF

c figure out how many radiating layers (or MPs) in this atmosphere, and check 
c that that it is less than or equal to kProfLayer
      IF (iStop .GE. iStart) THEN
        iNlay = (iStop-iStart+1)
        iDirection = +1                           !down look instr
      ELSE IF (iStop .LE. iStart) THEN
        iNlay = (iStart-iStop+1)
        iDirection = -1                           !up look instr
      END IF
      IF (iNLay .GT. kProfLayer) THEN
        write(kStdErr,*)'Error for atm # ',iC
        write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
        CALL DoSTOP
      END IF

c set the B.C.'s
      raTSpace(iC) = rTbdy
      raTSurf(iC)  = FindSurfaceTemp(rPressStart,rPressStop,rTSurf,
     $                     raProfileTemp,raPressLevels,iProfileLayers)

      raSatAngle(iC)=rAngle
      IF (abs(rAngle) .LE. 1.0e-4) THEN !nadir view
        rHeight = -1.0
        raSatHeight(iC) = -1.0
      ELSE
        IF (rHeight .gt. 0.0) THEN
          raSatHeight(iC) = rHeight   !height in m
        ELSE
          rHeight = -1.0
          raSatHeight(iC) = rHeight   !height in m
        END IF
      END IF
      raSatAzimuth(iC) = prof%satazi
      raSolAzimuth(iC) = prof%solazi
      iaNumLayer(iC) = iNlay

      write(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
      write(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
      write(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted= ',
     $                         rTsurf,raTSurf(iC)

c set the mixed path numbers for the current atmosphere, in direction of
c radiation travel
      DO iInt = 1,iNlay
        iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
      END DO

c use the solar on/off, thermal on/off etc. 
c sun is only on if 0 < prof%solzen < 90
      !!rakSolarAngle(iC) = abs(prof%sunang)   !!!RTP v 097-
      rakSolarAngle(iC) = prof%solzen          !!!RTP v 098+
      IF ((prof%solzen .GE. 0.0) .AND. (prof%solzen .LE. 90.0)) THEN
        IF ((rf1 .GE. 605.0) .AND. (rf2 .LE. 2830.0)) THEN
          iakSolar(iC) = +1 
        ELSEIF ((rf1 .LT. 605.0) .OR. (rf2 .GT. 2830.0)) THEN
          iakSolar(iC) = +0   !!! do not have a solar database yet
        END IF 
      ELSE 
        iakSolar(iC) = -1
      END IF
      raKSolarRefl(iC) = -1.0
      iaKThermal(iC)   = 0
      raKThermalAngle(iC) = -1.0

c      print *,'----> warning : set raKthermalangle = 53.3 (acos(3/5))'
c      raKThermalAngle(iC) = +53.13
c      print *,'----> so this will be used at all layers '
c      print *,'----> instead of varying the diffusivity angle'

      iakThermalJacob(iC) = 1
c use the solar on/off, thermal on/off etc. 
      kSolar        = iaKSolar(iC)
      IF (abs(raKSolarAngle(iC) - 90.0) .le. 1.0e-5) then
        write(kStdWarn,*) 'resetting solar angle = 90 to 89.9, iAtm = ',iC
        raKSolarAngle(iC) = 89.9
      END IF
      kSolarAngle   = raKSolarAngle(iC)
      kSolarRefl    = raKSolarRefl(iC)
      kThermal      = iaKThermal(iC)
      kThermalAngle = raKThermalAngle(iC)
      kThermalJacob = iakThermalJacob(iC)

      IF (kThermalAngle  .LT. 0) THEN
        kSetThermalAngle = -1   !use accurate angles lower down in atm
      ELSE
        kSetThermalAngle = +1   !use user specified angle everywhere
      END IF

      IF (iDirection .GT. 0) THEN
        !check things make sense for downlook in
        IF ((kSolarAngle .LT. 0.0) .OR. (kSolarAngle .GT. 90.0)) THEN
          write(kStdWarn,*) 'Warning! Resetting Solar Angle from ',kSolarAngle,' to 150.0'
          write(kStdWarn,*) 'and setting kSolar from ',kSolar, ' to -1 (solar = off)'
          kSolarAngle = 150.0
          kSolar      = -1
        END IF
        IF ((abs(kSolar) .NE. 1) .AND. (kSolar .NE. 0)) THEN
          write(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
          CALL DoSTOP 
        END IF
        IF (abs(kThermal) .GT. 1) THEN
          write(kStdErr,*)'need Thermal on/off parameter = -1/0/1'
          CALL DoSTOP 
        END IF
        IF (abs(kThermalJacob) .NE. 1) THEN
          write(kStdErr,*)'need ThermalJacob on/off parameter = -1/1'
          CALL DoSTOP 
        END IF
        !set the diffusivity angle in degrees
        !IF ((kThermalAngle .LT. 0.0).OR.(kThermalAngle .GT. 90.0)) THEN
        IF (kThermalAngle .GT. 90.0) THEN
          write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
          kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
        END IF
      END IF

      IF (iDirection .LT. 0) THEN
        IF ((kWhichScatterCode .EQ. 2) .OR. (kWhichScatterCode .EQ. 4)) THEN
          !set to nonsense values for uplooking instrument RTSPEC SCAT
          !as these CANNOT handle solar
          kSolar = -1          !!!RTPSEC, FIRST ORDER PERTURB cannot handle sun
          IF (kSolar .NE. iaKSolar(iC)) THEN
            write(kStdErr,*) 'in radnce4RTP, kSolar = -1 but you have a solar'
            write(kStdErr,*) 'angle for the profile!!!',kWhichScatterCode
            CALL DoStop
          END IF
          kSolarAngle   = -90.0
          kSolarRefl    = 0.0
          kThermal      = 0
          kThermalAngle = -45.0
          kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
          kThermalJacob = -1
        ELSE
          !set to nonsense values for uplooking instr kCARTA, DISORT,TWOSTR
          !!!kSolar = -1        !!!kCARTA nonscatter can handle this
          !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
          !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
          kSolarRefl    = 0.01     
          kThermal      = 0
          kThermalAngle = -45.0
          kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
          kThermalJacob = -1
        END IF
      END IF

c So if {\sf iakThermal(iI) = 0}, then {\sf rakThermalAngle(iI)} should be used
c with care.  If it is set at a negative value $x$, then for the upper
c layers the diffusive angle $acos(3/5)$ is used for the reflected
c thermal, while for the lower layers, a parameterized optimum
c diffusivity angle is used.  If it is set at a positive value, then for
c all layers, the diffusive angle $acos(x)$ is used for the
c reflected thermal. acos(3/5) = 53.1301 degrees

c      print *,'**************************************************'
c      kThermal = -1
c      kSolar = -1
c      print *,'*******  in rtp_interface,  kThermal = -1 ********'
c      print *,'*******  in rtp_interface,  kSolar   = -1 ********'
c      print *,'**************************************************'
      iakSolar(iC)          = kSolar
      rakSolarAngle(iC)     = kSolarAngle
      rakSolarRefl(iC)      = kSolarRefl
      iakThermal(iC)        = kThermal
      rakThermalAngle(iC)   = kThermalAngle
      iakThermalJacob(iC)   = kThermalJacob
      iaSetThermalAngle(iC) = kSetThermalAngle

      write(kStdWarn,*)'Solar on/off, Solar angle, Solar emiss = ',
     $             kSolar,kSolarAngle,kSolarRefl
      write(kStdWarn,*)'Thermal on/off,Thermal angle,Thermal Jacob =',
     $              kThermal,kThermalAngle,kThermalJacob

c now read in the emissivity values 
      iaSetEms(iC) = prof%nemis
      IF (iaSetEms(iC) .GT. kEmsRegions) THEN 
        write(kStdErr,*)'Cannot set so many emiss regions. Change' 
        write(kStdErr,*)'kEmsRegions in kcarta.param and recompile' 
        CALL DoSTOP 
      END IF 

      IF (iaSetEms(iC) .GE. 2) THEN
        DO i=1,iaSetEms(iC) 
          r1   = prof%efreq(i)
          rEms = prof%emis(i)
          write(kStdWarn,*) r1,rEms 
          raaaSetEmissivity(iC,i,1) = r1 
          raaaSetEmissivity(iC,i,2) = rEms 
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN 
            write(kStdErr,*)'Need emissivity between 0 and 1' 
            write(kStdErr,*)'check your emissivity values in file' 
            CALL DoSTOP 
          END IF 
        END DO 
      END IF

      IF ((iaSetEms(iC) .EQ. 1) .AND. (prof%emis(1) .GT. 1.0)) THEN 
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'rtpfile : has ONE emiss point, and emiss > 1'
        write(kStdWarn,*) '  assuming LIMB VIEW : RESET emissivity = 0'
        iaLimb(iC) = +1
        i = 1
        r1 = 6.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1 
        raaaSetEmissivity(iC,i,2) = rEms
        i = 2
        r1 = 36000000.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1 
        raaaSetEmissivity(iC,i,2) = rEms
      ELSEIF ((iaSetEms(iC) .EQ. 1) .AND. (prof%emis(1) .LE. 1.0)) THEN 
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'Fooling emissivity file so that it uses two points'
        write(kStdWarn,*) 'with constant emissivity',prof%emis(1)
        iaSetEms(iC) = 2
        i = 1
        r1 = 6.0
        rEms = prof%emis(1)
        IF (rEms .LT. 0.0) THEN 
          write(kStdErr,*)'Need emissivity between 0 and 1' 
          write(kStdErr,*)'   RTP file has ',prof%nemis,' emissivity points'
          write(kStdErr,*)'   first point (rf,rEms) = ',prof%efreq(1),prof%emis(1)
          CALL DoSTOP 
        END IF 
        raaaSetEmissivity(iC,i,1) = r1 
        raaaSetEmissivity(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms  
        i = 2
        r1 = 36000000.0
        rEms = prof%emis(1)
        raaaSetEmissivity(iC,i,1) = r1 
        raaaSetEmissivity(iC,i,2) = rEms 
        write(kStdWarn,*) r1,rEms  
      END IF
 
c now read in the solar refl values 
      iaSetSolarRefl(iC) = prof%nemis !! new, before it was nrho
      IF (iaSetSolarRefl(iC) .GT. kEmsRegions) THEN 
        write(kStdErr,*)'Cannot set so many solar refl regions. Change' 
        write(kStdErr,*)'kEmsRegions in kcarta.param and recompile' 
        CALL DoSTOP 
      END IF 
      IF (iaSetSolarRefl(iC) .LT. 1) THEN 
        write(kStdWarn,*)'No points in the solar refl file' 
        write(kStdWarn,*)'Will assume that we are using refl = (1-ems)/pi)' 
        iaSetSolarRefl(iC) = iaSetEms(iC) 
        DO i = 1,iaSetEms(iC) 
          !first is wavenumber, second is emissivity --> reflectance
          raaaSetSolarRefl(iC,i,1)  = raaaSetEmissivity(iC,i,1)
          raaaSetSolarRefl(iC,i,2)  = (1-raaaSetEmissivity(iC,i,2))/kPi
        END DO 
      ELSE IF (iaSetSolarRefl(iC) .LT. 2) THEN 
        write(kStdWarn,*)'For solar refl, Need > 1 point for interpolation'
        write(kStdWarn,*)'Fooling reflectivity file so that it uses two '
        write(kStdWarn,*)'points, with constant reflectivity',prof%rho(1)
        iaSetSolarRefl(iC) = 2
        i = 1
        r1 = 3.6
        rEms = prof%rho(1)
        raaaSetSolarRefl(iC,i,1) = r1 
        raaaSetSolarRefl(iC,i,2) = rEms 
        write(kStdWarn,*) r1,rEms  
        i = 2
        r1 = 3600.0
        rEms = prof%rho(1)
        raaaSetSolarRefl(iC,i,1) = r1 
        raaaSetSolarRefl(iC,i,2) = rEms 
        write(kStdWarn,*) r1,rEms  
      ELSE
        DO i=1,iaSetSolarRefl(iC) 
          r1   = prof%efreq(i)   !!new rfreq = efreq
          rEms = prof%rho(i)
          write(kStdWarn,*) r1,rEms 
          raaaSetSolarRefl(iC,i,1) = r1 
          raaaSetSolarRefl(iC,i,2) = rEms 
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN 
            write(kStdErr,*)'Need reflectance between 0 and 1' 
            write(kStdErr,*)'check your reflectance values in file' 
            CALL DoSTOP 
          END IF 
        END DO 
      END IF       

      print *,'booboo d ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      !now see if there is a cloud to be used with this atmosphere

      cfrac12 = prof%cfrac12 

      print *,'booboo d1 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      ctype1  = int(prof%ctype)
      cfrac1  = prof%cfrac
      cngwat1 = prof%cngwat
      rSize1  = prof%cpsize

      print *,'booboo d2 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      ctype2  = int(prof%ctype2)
      cfrac2  = prof%cfrac2
      cngwat2 = prof%cngwat2
      rSize2  = prof%cpsize2

      print *,'booboo d3 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      i4ctype1 = prof%ctype
      i4ctype2 = prof%ctype2

      print *,'booboo e ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

      print *,'xmama0',i4ctype1,i4ctype2
      print *,'xmama1',ctype1,cfrac1,cngwat1,rSize1
      print *,'xmama2',ctype2,cfrac2,cngwat2,rSize2
      print *,'xmama12',cfrac12
      print *,prof%clrflag,prof%ctype,prof%ctype2,prof%upwell,prof%findex,prof%atrack,prof%xtrack,prof%ifov,prof%robsqual,prof%itype

      IF ((cfrac1 .LE. 0) .AND. (cfrac2 .GT. 0)) THEN
        write(kStdErr,*) 'kCARTA assumes if cfrac1 > 0 then cfrac2 >= 0'
        write(kStdErr,*) 'kCARTA assumes if cfrac1 = 0 then cfrac2  = 0'
        write(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        CALL DOSTOP
      END IF

      IF (cfrac12 .GT. max(cfrac1,cfrac2)) THEN
        write(kStdErr,*) 'kCARTA assumes cfac12 <= max(cfrac1,cfrac2)'
        write(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        CALL DOSTOP
      END IF

      IF (prof%cfrac .gt. 0.0) THEN
        !!!first cloud is easy to do
        i = 1
        iaCtype(i) =  prof%ctype       !cloud type 1=cirrus 2=water etc
        raCemis(i) = 1.0               !assume cloud totally emissive
        raCngwat(i) = prof%cngwat*1.00 !IWP
        raCpsize(i) = prof%cpsize*1.00 !in microns
        raCprtop(i) = prof%cprtop
        raCprbot(i) = prof%cprbot

        i = 2
        iaCtype(i) =  prof%ctype2        !cloud type 1=cirrus 2=water etc
        raCemis(i) = 1.0                 !assume cloud totally emissive
        raCngwat(i) = prof%cngwat2*1.00  !IWP
        raCpsize(i) = prof%cpsize2*1.00  !in microns
        raCprtop(i) = prof%cprtop2
        raCprbot(i) = prof%cprbot2

       ELSE
        cfrac1  = 0.0            !assume clear sky, use dummy values
        cfrac2  = 0.0            !assume clear sky, use dummy values
        cngwat1 = 0.0
        cngwat2 = 0.0
        DO i = 1,kMaxClouds
          iaCtype(i) =  0                !cloud type
          raCemis(i)  = 0.0              !assume cloud totally emissive
          raCngwat(i) = 0.0              !IWP
          raCpsize(i) = 1.0              !in microns
          raCprtop(i) = 500.0
          raCprbot(i) = 600.0
        END DO
      END IF

      RETURN
      END

c************************************************************************
c this subroutine quickly sets up stuff for ONE atmosphere
c there could be more than one cloud
      SUBROUTINE SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop,
     $    cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,iNclouds_RTP,iaKsolar,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
     $    iBinORasc,caaCloudFile,iaNML_Ctype,
     $    iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
     $    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
     $    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input params
      INTEGER iakSolar(kMaxAtm),ctype1,ctype2
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      !these are the cloud parameters read in from the RTP file
      REAL Cfrac,cfrac1,cfrac2,cfrac12,raCemis(kMaxClouds),cngwat1,cngwat2
      REAL raCprtop(kMaxClouds),  raCprbot(kMaxClouds)
      REAL raCngwat(kMaxClouds),  raCpsize(kMaxClouds)
      INTEGER iaCtype(kMaxClouds),iBinORasc,iNclouds_RTP
      CHARACTER*80 caaCloudFile(kMaxClouds)
      INTEGER iaNML_Ctype(kMaxClouds)
c output params, 
c     above set into the cloud parameters .....
c iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaCloudName(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c raaPrBdry is the pressure boundaries for the atms
      REAL raaPrBdry(kMaxAtm,2)      
c raPresslevls are the KLAYERS pressure levels
c iProfileLayers = tells how many layers read in from RTP or KLAYERS file
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iProfileLayers
c this tells if the cloud, when "expanded", has same IWP or exponentially
c decreasing IWP
      REAL raExp(kMaxClouds)
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
      INTEGER iNatm
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iJ1,iI,iIn,iJ,iScat,iaTemp(kMixFilRows),iTop,iBot,iNum,iErr
      REAL rPT,rPB,rP1,rP2,rSwap,r1,r2
      CHARACTER*80 caName
      INTEGER FindCloudLayer
c these are to check that the scattering table names are unique
      INTEGER iaTable(kCloudLayers*kMaxClouds),iWhichScatterCode,iDefault
      CHARACTER*80 caaTable(kCloudLayers*kMaxClouds)
c these are to match iaCtype vs iaNML_Ctype
      INTEGER iFound1,iFound2,iNclouds_RTPX

      IF (kAllowScatter .EQ. -1) THEN
        write(kStdErr,*) 'bkcarta.x (basic) version does not allow scattering'
        write(kStdErr,*) 'Please either use Makefile to compile/run kcarta.x'
        write(kStdErr,*) 'which allows scattering, or modify our .nml and/or'
        write(kStdErr,*) '.rtp files, so that only basic kCARTA is used'
        CALL DoStop
      END IF    

      iWhichScatterCode = 6         !!RAYLEIGH in CLEAR SKY, nir/vis/uv
      iWhichScatterCode = 5         !!PCLSAM
      iWhichScatterCode = 4         !!r = r0 + r1 = perturb (not yet done)
      iWhichScatterCode = 3         !!DISORT
      iWhichScatterCode = 2         !!RTSPEC
      iWhichScatterCode = 1         !!TWOSTREAM  DEFAULT
      iWhichScatterCode = 0         !!simple absorb; directly goes to rad_main

      iDefault = 5

      iWhichScatterCode = 6         !!RAYLEIGH in clear sky
      iWhichScatterCode = 5         !!PCLSAM

      IF (iDefault .NE. iWhichScatterCode) THEN
        print *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
      END IF

      IF (iWhichScatterCode .EQ. 6) THEN
        kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 5) THEN
        kWhichScatterCode = 5        !use PCLSAM
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 4) THEN
        kWhichScatterCode = 4        !use r = r0 + r1 = perturb
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 3) THEN
        kWhichScatterCode = 3        !use Disort
        kScatter          = 1        !use this setting
        kDis_Pts          = 400      !do 1 every 400 pts
      ELSEIF (iWhichScatterCode .EQ. 2) THEN
        kWhichScatterCode = 2        !use RTSPEC
        kScatter          = 1        !use this setting  SingleScatter
        kScatter          = 3        !use this setting  Hybrid
        IF (kScatter .NE. 3) THEN
          write (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
        END IF
      ELSEIF (iWhichScatterCode .EQ. 1) THEN
        kWhichScatterCode = 1        !use TwoStream
        kScatter          = 1        !use one run of TwoStream
      ELSEIF (iWhichScatterCode .EQ. 0) THEN
        kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
        kScatter          = 1        !   
      END IF

      IF ((iakSolar(1) .GE. 0)  .AND. (kWhichScatterCode .EQ. 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
        CALL DoStop
      END IF

      IF ((iakSolar(1) .GE. 0)  .AND. (kWhichScatterCode .EQ. 4)) THEN
        write(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
        CALL DoStop
      END IF
 
      IF (ctype1 .NE. iaCtype(1)) THEN
        write(kStdErr,*) 'hmm ctype1,iaCtype(1) = ',ctype1,iaCtype(1)
        CALL DoStop
      END IF
      IF (ctype2 .NE. iaCtype(2)) THEN
        write(kStdErr,*) 'hmm ctype2,iaCtype(2) = ',ctype2,iaCtype(2)
        CALL DoStop
      END IF

      iNclouds_RTPX = iNclouds_RTP

      IF ((ctype1 .LE. 0) .AND. (ctype2 .LE. 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 0'
        iNclouds_RTP = 0
        DO iI = 1,kMaxClouds
          raaaCloudParams(iI,1,1) = 0.0
        END DO
        ctype1 = -9999
        ctype2 = -9999
        iaCtype(1) = -9999 
        iaCtype(2) = -9999 
        RETURN
      END IF

      IF ((ctype1 .LE. 0) .AND. (ctype2 .GT. 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' swapping info for clouds1,2'
 
        ctype1       = ctype2
        cfrac1       = cfrac2
        cngwat1      = cngwat2
        iaCtype(1)   = ctype1

        raCprTop(1)  = raCprTop(2)
        raCprBot(1)  = raCprBot(2)
        raCpSize(1)  = raCpSize(2)
        raCngwat(1)  = raCngwat(2)
        iNclouds_RTP = 1

        ctype2       = -9999
        iaCtype(2)   = -9999 
        raCngwat(2)  = 0.0
        cngwat2      = 0.0
        cfrac2       = 0.0
        cfrac12      = 0.0
      END IF

      IF ((ctype1 .GT. 0) .AND. (ctype2 .LE. 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
        iNclouds_RTP = 1
      END IF
 
      iFound1 = -1
      DO iI = 1,iNclouds_RTPX
        IF (iaNML_Ctype(iI) .EQ. ctype1) THEN
          iFound1 = iI
          GOTO 1234
        END IF
      END DO
 1234 CONTINUE
      IF ((iFound1 .LT. 0) .AND. (ctype1 .GT. 0)) THEN
        write(kStdErr,*) 'ctype1 cfrac1 cngwat1 = ',ctype1,cfrac1,cngwat1
        write(kStdErr,*) 'Could not find a match between ctype1 = ',ctype1,' and iaNML_Ctype'
        CALL DoStop
      END IF

      iFound2 = -1
      DO iI = 1,iNclouds_RTPX
        IF (iaNML_Ctype(iI) .EQ. ctype2) THEN
          iFound2 = iI
          GOTO 2345
        END IF
      END DO
 2345 CONTINUE
      IF ((iFound2 .LE. 0) .AND. (ctype2 .GT. 0)) THEN
        write(kStdErr,*) 'ctype2 cfrac2 cngwat2 cfrac12 = ',ctype2,cfrac2,cngwat2, cfrac12
        write(kStdErr,*) 'Could not find a match between ctype2 = ',ctype2,' and iaNML_Ctype'
        CALL DoStop
      END IF

      IF ((iNclouds_RTP .EQ. 3) .AND. (iFound1 .GT. 0) .AND. (iFound2 .GT. 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found two clouds, resetting'
        iNclouds_RTP = 2
      ELSEIF ((iNclouds_RTP .EQ. 3) .AND. (iFound1 .GT. 0) .AND. (iFound2 .LE. 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
        iNclouds_RTP = 1
      ELSEIF ((iNclouds_RTP .EQ. 3) .AND. (iFound1 .LE. 0) .AND. (iFound2 .GE. 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
        iNclouds_RTP = 1
      ELSEIF ((iNclouds_RTP .EQ. 3) .AND. (iFound1 .LE. 0) .AND. (iFound2 .LE. 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found zero clouds, resetting'
        iNclouds_RTP = 0
        ctype1 = -9999
        iaCtype(1) = -9999
        raCngwat(1) = 0.0
        cngwat1     = 0.0
        cfrac1       = 0.0

        ctype2       = -9999
        iaCtype(2)   = -9999 
        raCngwat(2)  = 0.0
        cngwat2      = 0.0
        cfrac2       = 0.0
        cfrac12      = 0.0
      END IF

      IF (iWhichScatterCode .EQ. 0) THEN
        write(kStdWarn,*) 'ONE Purely absorptive cloud is set from rtp!!!!!'
        iScatBinaryFile         = iBinORasc
        iNClouds                = 1
        caaScatter(1)           = caaCloudFile(iFound1)
        raaScatterPressure(1,1) = raCprTop(1)
        raaScatterPressure(1,2) = raCprBot(1)
        raScatterDME(1)         = raCpSize(1)
        raScatterIWP(1)         = raCngWat(1)
        write(kStdWarn,*) 'Cloud datafile  is : '
        write(kStdWarn,222) caaScatter(1)
        write(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1),
     $      raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
        GOTO 333
      END IF

c from /asl/packages/sartaV105/yukyung_readme.txt
c      prof%udef(11,:) = cngwat
c      prof%udef(12,:) = cpsize
c      prof%udef(13,:) = cprtop
c      prof%udef(14,:) = cprbot
c      prof%udef(15,:) = cfrac
c      prof%udef(16,:) = "cfrac12", fraction of FOV containing both clouds
c      prof%udef(17,:) = ctype {currently not used}
c      prof%udef(18,:) = cemis for cloud 2

c      raExp(1)               = 0.0      !same amount in all layers
c      iScatBinaryFile        = iBinORasc
c      iNClouds               = 1
c      iaCloudNumLayers(1)    = 1
c      iaaScatTable(1,1)      = 1
c      caaCloudName(1)        = 'RTP cloud'
c      caaaScatTable(1,1)     = cfile xxxx caaCloudFile(iI)
c      raaaCloudParams(1,1,1) = cngwat
c      raaaCloudParams(1,1,2) = cpsize
c      raaPCloudTop(1,1)      = cprtop
c      raaPCloudBot(1,1)      = cprbot
c      iaCloudNumAtm(1)       = 1
c      iaaCloudWhichAtm(1,1)  = 1
c      iaPhase(1)             = -1       !default to HG phase function

      write(kStdWarn,*) ' '
      iScatBinaryFile        = iBinORasc
      iNClouds               = iNclouds_RTP

      DO iI = 1,iNclouds_RTP
        raExp(iI)               = 0.0      !same amount in all layers
        iaCloudNumLayers(iI)    = 1
        iaaScatTable(iI,1)      = iI
        caaCloudName(iI)        = 'RTP cloud'
        IF (iI .EQ. 1) THEN
          caaaScatTable(iI,1)     = caaCloudFile(iFound1)
        ELSEIF (iI .EQ. 2) THEN
          caaaScatTable(iI,1)     = caaCloudFile(iFound2)
        END IF
        raaaCloudParams(iI,1,1) = raCngwat(iI)
        raaaCloudParams(iI,1,2) = raCpsize(iI)
        raaPCloudTop(iI,1)      = raCprtop(iI)
        raaPCloudBot(iI,1)      = raCprbot(iI)
        iaCloudNumAtm(iI)       = 1
        iaaCloudWhichAtm(iI,1)  = 1
        iaPhase(iI)             = -1       !default to HG phase function

        write(kStdWarn,*)    'cloud info for RTP cloud # ',iI
        write (KStdWarn,222) 'cloud file    = ',caaCloudFile(iI)
        write (kStdWarn,*)   'cloud top     = ',raCprtop(iI),' mb'
        write (kStdWarn,*)   'cloud bot     = ',raCprbot(iI),' mb'
        write (kStdWarn,*)   'cloud IWP     = ',raCngwat(iI),' gm m-2'
        write (kStdWarn,*)   'particle size = ',raCpsize(iI),' um'
        IF (iI .EQ. 1) THEN
          write (kStdWarn,*)   'cloud frac    = ',cfrac1
        ELSEIF (iI .EQ. 2) THEN
          write (kStdWarn,*)   'cloud frac    = ',cfrac2
        END IF
      END DO

      !!now have to stretch out the cloud if necessary
      CALL ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,
     $     caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable,
     $     iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp,
     $     raPressLevels,iProfileLayers,
     $     raFracTop,raFracBot,raPressStart,raPressStop,iNatm)

c now start checking the info
      DO iIn=1,iNclouds
        caName = caaCloudName(iIn)
        iJ     = iaCloudNumLayers(iIn)
        iaCloudNumLayers(iIn) = iJ
        write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '

c set individual cloud layer parameters but STRETCH the cloud out as necessary
c from pressure level rPT to pressure level rPB
c note it will occupy the entire layer
        DO iJ1=1,iJ
          !top and bottom pressures CloudName/Type  IWP/LWP DME          
          rPT = raaPCloudTop(iIn,iJ1)
          rPB = raaPCloudBot(iIn,iJ1)

          IF (rPT .GT. rPB) THEN
            rSwap = rPT
            rPT   = rPB
            rPB   = rSwap
            write (kStdWarn,*) 'Swapped cloud top & bottom pressures'
          END IF

          iTop = FindCloudLayer(rPT,raPressLevels,iProfileLayers)
          iBot = FindCloudLayer(rPB,raPressLevels,iProfileLayers)
          iNum = iTop
          IF ((iTop - iBot) .LT. 0) THEN
            write (kStdErr,*) 'the top of your cloud is below the bottom'
            CALL DoStop
          END IF

          iaaCloudWhichLayers(iIn,iJ1) = iNum  !layer number wrt 1 ..kProfLayer
          rP1 = raaaCloudParams(iIn,iJ1,1)     !IWP
          rP2 = raaaCloudParams(iIn,iJ1,2)     !mean size

          iScat=iaaScatTable(iIn,iJ1)
          caName=caaaScatTable(iIn,iJ1)
          write(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
          write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
          write(kStdWarn,*) '   mean particle size (um)  = ',rP2
          write(kStdWarn,*) '   scatter table number = ',iScat
          write(kStdWarn,222) caName
        END DO 

c set how many, and which atmospheres to use with this cloud
        iNum=iaCloudNumAtm(iIn)
        IF (iNum .GT. iNatm) THEN
          write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!' 
          write(kStdErr,*)'*SCATTR wants to use',iNum,' atmospheres!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DOStop
        END IF
        !check which atmospheres to use this cloud         
        DO iJ = 1,iNum
          iaTemp(iJ)=iaaCloudWhichAtm(iIn,iJ)
        END DO
        DO iJ = 1,iNum
          IF (iaTemp(iJ) .GT. iNatm) THEN
            write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!' 
            write(kStdErr,*)'*SCATTR wants to use atmosphere #',iaTemp(iJ) 
            write(kStdErr,*)'please check and retry!'
            CALL DOStop
          END IF
        END DO

        write(kStdWarn,*) 'number of atms for cloud is ',iNum
        write(kStdWarn,*) '  atmospheres to be used with this cloud  : '
        write(kStdWarn,*)(iaTemp(iJ),iJ=1,iNum)        
        write(kStdWarn,*) '  '

      END DO
cccccccccccccccccccccc now check the info
 
      write(kStdWarn,*) 'finished preprocessing *SCATTR .. checking info ...'

      write(kStdWarn,*) 'checking cloud boundaries lies in start/stop press...'
      !check that cloud boundaries lie within those defined for atmosphere
      DO iIn = 1,iNClouds
        !these would be cloud top and bottom pressures
        r1 = raPressLevels(iaaCloudWhichLayers(iIn,1)+1)
        r2 = raPressLevels(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
        !check top pressure
        DO iJ = 1,iaCloudNumAtm(iIn)
          iI  = iaaCloudWhichAtm(iIn,iJ)
          rPT = raaPrBdry(iI,1)         !start pressure
          rPB = raaPrBdry(iI,2)         !stop pressure
          IF (rPT .GT. rPB) THEN      !atm is for down look instr
            rP1 = rPT
            rPT = rPB
            rPB = rP1
          END IF
          !check top pressure
          IF (r1 .LT. rPT) THEN
            write(kStdErr,*)'*RADNCE defines top pressure for atmosphere'
            write(kStdErr,*)'number ',iI,' as ',rPT
            write(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
            write(kStdErr,*)'that atmosphere; cloud top at ',r1
            iErr = 1
            CALL DOStop
          END IF
          !check bot pressure
          IF (r2 .GT. rPB) THEN
            write(kStdWarn,*)'*RADNCE defines bottom pressure for atmosphere'
            write(kStdWarn,*)'number ',iI,' as ',rPB
            write(kStdWarn,*)'*SCATTR says to use cloud number ',iIn,' in'
            write(kStdWarn,*)'that atmosphere; cloud bottom at',r2
            write(kStdWarn,*)'Resetting r2 ...'
            r2 = rPB
          END IF
        END DO
      END DO

      write(kStdWarn,*) 'checking cloud layers sequential ...'
      !check that the layers for a cloud are sequential eg 16,15,14
      DO iIn=1,iNclouds
        !if there is only one layer in the cloud, things OK, else
        IF (iaCloudNumLayers(iIn)  .GT. 1) THEN
          iJ = 1
          iJ1 = iaaCloudWhichLayers(iIn,iJ)
          DO iJ = 2,iaCloudNumLayers(iIn)
            iScat = iaaCloudWhichLayers(iIn,iJ)
            IF (iScat .GE. iJ1) THEN
              write(kStdErr,*) 'checking cloud # ',iIn
              write(kStdErr,*) 'layer ',iJ,' is not below preceding layer'
              write(kStdErr,*) 'please check and retry'
              CALL DoStop
            END IF
            IF ((iJ1-iScat) .GT. 1) THEN
              write(kStdErr,*) 'checking cloud # ',iIn
              write(kStdErr,*) 'layers not sequential!!',iJ1,iScat
              write(kStdErr,*) 'please check and retry'
              CALL DoStop
            END IF
            iJ1 = iScat
          END DO
        END IF
      END DO

c check that the scattering tables are unique within a cloud
      write(kStdWarn,*) 'checking scattering tables unique within a cloud ...'
      DO iIn = 1,iNclouds
        DO iJ = 1,iaCloudNumLayers(iIn)
          iI = iaaScatTable(iIn,iJ)
          caName = caaaScatTable(iIn,iJ)
          DO iJ1 = iJ+1,iaCloudNumLayers(iIn)
            IF (iI .EQ. iaaScatTable(iIn,iJ1)) THEN
              write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
              write(kStdWarn,*) 'found nonunique scattering table numbers'
              write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE
     $ profile layer temperature'
            END IF
            IF (caName .EQ. caaaScatTable(iIn,iJ1)) THEN
              write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
              write(kStdWarn,*) 'found nonunique scattering table file names'
              write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE
     $ profile layer temperature'
            END IF
          END DO
        END DO
      END DO 

c if this test is successfully passed, then do the next check!!!
c check across all clouds that the scattering tables are unique
c map this code to rtspec.f
c these are to check that the scattering table names are unique
      DO iIn = 1,kMaxClouds*kCloudLayers
        iaTable(iIn) = -1
        caaTable(iIn) = '                                                     '
      END DO
      write(kStdWarn,*) 'checking scattering tables unique thru all clouds ...'
      DO iIn = 1,iNclouds
        DO iJ = 1,iaCloudNumLayers(iIn)
          iI = iaaScatTable(iIn,iJ)
          caName = caaaScatTable(iIn,iJ)
          IF (iaTable(iI) .LT. 0) THEN  !nothing associated with this yet
            iaTable(iI) = 1
            caaTable(iI) = caName
          ELSE                          !check to see file names are the same
            IF (caaTable(iI) .NE. caName) THEN
              write(kStdErr,*)'Scattering table #',iI,' <-> ',caaTable(iI)
              write(kStdErr,*)'for same scattering table, new cloud in
     $ *SCATTR is associating file ',caName
              write(kStdErr,*)'please check and retry'
              CALL DoStop
            END IF
          END IF
        END DO
      END DO 

 222  FORMAT(A80) 
 333  CONTINUE

      RETURN
      END

c************************************************************************
c this subroutine reads 48 regression profiles to see which is closest to the
c current profile, so that kCARTA can figure out which NLTE temps to use!!!
      SUBROUTINE NLTE_RegrTemp(raa48Temp,raa48Press,ia48layers)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c the 48 regression profile kinetic temperatures
      REAL raa48Temp(kMaxLayer,kRegrProf)
      REAL raa48Press(kMaxLayer,kRegrProf)
      INTEGER ia48layers(kRegrProf)

      REAL rAmt,rT,rP,rPP,rH,rdP,rdT,raPressLevels(kMaxLayer+1)
      INTEGER iDownWard

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,iK,iJ,iI,iRTP,iL1,iProfileLayers,iGasInRTPFile,iFindJ
      REAL pProf(kMaxLayer),pTemp(kMaxLayer),MGC,plays(kMaxLayer)
      
      integer rtpopen, rtpread, rtpwrite, rtpclose
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      MGC = 8.314674269981136  

      DO iI = 1,kRegrProf
        ia48layers(iI) = 0
        DO iJ = 1,kMaxLayer
          raa48Temp(iJ,iI)  = 0.0
          raa48Press(iJ,iI) = 0.0
        END DO
      END DO

      fname = kRegrFile

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
      write(kStdWarn,*)  'read open status = ', status

      iRTP = 48
      DO iJ = 1, iRTP      
c        write(kStdWarn,*) 'Reading temperature regression profile ',iJ

        status = rtpread(rchan, prof)
        IF (prof%plevs(1) .lt. prof%plevs(prof%nlevs)) THEN
          !layers are from TOA to the bottom
          iDownWard = -1
        ELSE
          !layers are from GND to the top
          iDownWard = +1
        END IF

        iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
        iProfileLayers = iL1
        ia48layers(iJ) = iL1
        iGasInRTPFile = head%ngas              !!! number of gases

        IF (prof%nlevs .GT. kMaxLayer+1) THEN
          write(kStdErr,*) 'this routine compiled for ',kMaxLayer,' layers'
          write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
          write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
          CALL DoStop
        END IF

        DO i = 1,prof%nlevs
          j = iFindJ(kMaxLayer+1,I,iDownWard)        !!!!notice the kProf+1
          raPressLevels(j) = prof%plevs(i)            !!!!in mb
        END DO

        DO i = 1,prof%nlevs-1
          pProf(i) = raPressLevels(i) - raPressLevels(i+1)
          pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
        END DO

c now loop (for water only) in the supplied profile
        DO i = 1, prof%nlevs - 1
          j = iFindJ(kMaxLayer,I,iDownWard)
          rT   = prof%ptemp(i)
          plays(i) = (prof%plevs(i)-prof%plevs(i+1))/
     $               log(prof%plevs(i)/prof%plevs(i+1))
          rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
          rP   = plays(i)               !need pressure in mb
          raa48Temp(j,iJ)  = rT
          raa48Press(j,iJ) = rP
        END DO

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
        DO i = prof%nlevs, kMaxLayer
          j = iFindJ(kMaxLayer,I,iDownWard)
          raa48Temp(j,iJ)  = 300.0
          raa48Press(j,iJ) = 1200.0 + i
        END DO

      END DO

      write (kStdWarn,*) 'success : read in 48 regression profiles ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status
      kProfileUnitOpen = -1

c      DO iI = 1,kRegrProf
cc        DO iJ = kMaxLayer-ia48layers(iI)+1,kMaxLayer
c        DO iJ = 1,kMaxLayer
c          print *,iI,ia48layers(iI),iJ,raa48Press(iJ,iI),raa48Temp(iJ,iI)
c        END DO
c        print *,' '
c      END DO

      RETURN
      END

c************************************************************************
c this subroutine prints out error messages
      SUBROUTINE FindError(rAmt,rT,rP,rPP,iIDgas,iCnt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables, to process
      REAL rAmt,rT,rP,rPP
      INTEGER iIDGas,iCnt

c local vars 
      INTEGER iError
 
      iError = -1

      IF ((rAmt .lt. 0.0) .OR. (rAmt. gt. 1.0e3)) THEN
        WRITE(kStdWarn,1080)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rAmt
        iError = 1
        rAmt = 0.0
        !CALL DoStop
      END IF

      IF ((rT .lt. 0.0) .OR. (rT. gt. 1.0e3)) THEN
        WRITE(kStdWarn,1081)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rT
        iError = 1
        rT = 0.0
        !CALL DoStop
      END IF

      IF ((rP .lt. 0.0) .OR. (rP .GT. 1.0e5)) THEN
        WRITE(kStdWarn,1082)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rP
        iError = 1
        rP = 0.0
        !CALL DoStop
      END IF

      IF ((rPP .lt. 0.0) .OR. (rPP .GT. 1.0e5)) THEN
        WRITE(kStdWarn,1083)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rPP
        iError = 1
        rPP = 0.0
        !CALL DoStop
      END IF

      IF (iError .EQ. 1) THEN
        rP = 1.0e3
        rPP = 1.0e-3
        rT = 300.0
        rAmt = 0.000000
        write(kStdWarn,4321) iIDGas,iCnt,rAmt,rT,rP,rPP
      END IF

        
 1111 FORMAT('gasID, layer = ',I5,I5,F12.5)
 1080 FORMAT('negative or bad gas amount in PRFILE profile file')
 1081 FORMAT('negative or bad gas temp in PRFILE profile file')
 1082 FORMAT('negative or bad layer pressure in PRFILE profile file')
 1083 FORMAT('negative or bad gas partial press in PRFILE profile file')
 4321 FORMAT('Reset RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

      RETURN
      END

c************************************************************************
c this subroutine finds the position where we wanna store stuff
      SUBROUTINE FindIndexPosition(iID,iNumGases,iaInputOrder,iFnd,iGasIndex)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iID,iNumGases        !current gasID, number gases from MOL/XSC
      INTEGER iaInputOrder(kMaxGas)!list of gases from MOL/XSCGAS
c output variables
      INTEGER iFnd,iGasIndex       !is RTP gas in the MOL/XSC gas list?
                                   !if so, where is it in the list?

      iFnd = -1
      iGasIndex = 1
 999  CONTINUE
      IF (iaInputOrder(iGasIndex) .EQ. iID) THEN
        iFnd = 1
      END IF
      IF ((iFnd .LT. 0) .AND. (iGasIndex .LT. iNumGases)) THEN
        iGasIndex = iGasIndex+1
        GO TO 999
      END IF
     	
      RETURN
      END

c************************************************************************
c this subroutine adds on the continuum flag
      SUBROUTINE ContinuumFlag(iIDGas,iaCont)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iaCont(kMaxGas),iIDGas
	         
      iaCont(iIDgas)  =  1       !continuum "always" included
      !water is a special case
      IF ((iIDGas .EQ. 1) .AND. (kCKD .GE. 0)) THEN
        iaCont(iIDgas) = 1
      ELSE IF ((iIDGas .EQ. 1) .AND. (kCKD .LT. 0)) THEN
        iaCont(iIDgas) = -1
      END IF

      RETURN
      END

c************************************************************************
c this adds on the water profile info (gasID = 1) to gasID 101,102
c and to gasID 103 = heavy water
c if the user wants the effects of water continuum added on	         
      SUBROUTINE AddWaterContinuumProfile(iaGases,iNumberofGasesRead,
     $        iaWhichGasRead,iaInputOrder,iNumGases,
     $        raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raaHeight(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
      INTEGER iNumberOFGasesRead,iaInputOrder(kMaxGas)          
     
      INTEGER iIDGas,iFound,iP,iGasIndex

      DO iIDGas = kNewGasLo,kNewGasHi+1
        IF ((iaGases(iIDGas) .EQ. 1) .AND. (iaGases(1) .EQ. 1)) THEN
          write(kStdWarn,*)'Using water profile for gasID ',iIDGas
          iNumberOfGasesRead     = iNumberOfGasesRead + 1
          iaWhichGasRead(iIDgas) = 1
          iFound    = -1
          iGasIndex = 1
 777      CONTINUE
          IF (iaInputOrder(iGasIndex) .EQ. iIDgas) THEN
            iFound = 1
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

      RETURN
      END

c************************************************************************
c this function finds what index to put
      INTEGER FUNCTION iFindJ(iL,I,iDownWard)

      IMPLICIT NONE

      INTEGER iL,i,iDownWard

      INTEGER j

      IF (iDownWard .EQ. -1) THEN
        j = iL - i + 1 
      ELSE
        j = i
      END IF

      iFindJ = j

      RETURN
      END

c************************************************************************
c           these functions deal with reading CLOUD PROFILES             
c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format
c same as READRTP except now it looks for 
c   h.ngas > 4 : for 4 gases + (1,2 or 3) cloud types
c   h.glist = cobinations of (201,202,203)
c     water drop           : ctype=101     in prof%gas_201
c     Baran ice aggragates : ctype=201     in prof%gas_202
c     andesite dust        : ctype=301     in prof%gas_203
c   h.gunit = 1            : molecules/cm2 for gases  01:63
c           = 5            : g/cm2         for clouds 201+
c
c ---------------> no units conversion required <---------------
      SUBROUTINE READRTP_CLD100LAYER(iRTP,iProfileLayers,
     $    caPFName,caCloudPfName,iNclouds_RTP,
     $    caaCloudFile,iaCloudScatType,
     $    raPresslevels,iBinOrAsc,
     $    iaaRadLayer,iNumLayer,iaKsolar,
c these are the outputs, just relevant for the 100 profile layers clouds
     $    iNclouds2,   ! added ESM
     $    iaCldTypes,raaKlayersCldAmt,
     $    ctype1,ctype2,cfrac1,cfrac2,
c these are the outputs, as also set from SetRTPCloud
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
     $    iScatBinaryFile,iNclouds3, ! added ESM
     $    iaCloudNumLayers,caaCloudName,
     $    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
     $    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
     $    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels2, ! added ESM
     $    iProfileLayers2 ) ! added ESM

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'
      include 'rtpdefs.f'
      INTEGER iplev
      INTEGER inatm   ! Added ESM 
      INTEGER iProfilelayers2   ! Added ESM 
      
      REAL raaPrBdry(kMaxAtm,2)    ! Added ESM 
      include '../INCLUDE/KCARTA_database.param'

c input params ---------------------------------------------------->
c   raaTemp/Press  = current gas profile parameters
c   iNpathClds     = total number of paths to be read in 
c                    (iNclouds*kProfLayers)
c   iProfileLayers = actual number of layers per gas profile (<=kProfLayer)
c   caCloudfPfName = name of file containing user supplied profiles
c   raLayerHeight  = heights of layers in km
c   iRTP           = which profile to read in
c   iNclouds_RTP   = how many clouds are claimed to be in the new .rtp file
      REAL raPressLevels(kProfLayer+1)
      REAL raPressLevels2(kProfLayer+1)
      INTEGER ctype1,ctype2
      INTEGER iakSolar(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer
      INTEGER iRTP,iProfileLayers,iNclouds_RTP,iBinOrAsc
      INTEGER iNclouds2, iNclouds3  ! added ESM
      INTEGER iaCloudScatType(kMaxCLouds)
      CHARACTER*80 caPFname,caCloudPfname,caaCloudFile(kMaxClouds)
c output params  -------------------------------------------------->
c   kMaxClouds == 5
c   raaKlayersCldAmt = cloud profiles(s) in g/m2
c   iNumCLds       = number of clouds
c   iaCldTypes     = type(s) of cloud
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
      INTEGER iNclouds,iaCldTypes(kMaxClouds)
c output params, above set into the cloud parameters ---------------->
c iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaCloudName(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c this tells if the cloud, when "expanded", has same IWP or exponentially
c decreasing IWP
      REAL raExp(kMaxClouds)
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds),cfrac1,cfrac2
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      INTEGER iaCtype(kMaxClouds)
      REAL    raCemis(kMaxClouds),raCngwat(kMaxCLouds),raCpsize(kMaxClouds)
      REAL    raCprtop(kMaxClouds),raCprbot(kMaxClouds)
 
c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer   i,j,k,iG,iH,iGX,iDownward,iGasIndex,iCldIndex,iMapIndex,iDiv
      REAL      raHeight(kProfLayer+1),pProf(kProfLayer)
      REAL      kcraPressLevels(kProfLayer+1)
      REAL      kcRTP_pTop,kcRTP_pBot, deltaP
      INTEGER   kcRTPTop,kcRTPBot,iL1,iGasInRTPFile,iaCldGasID(3)
      INTEGER   iDefault,iWhichScatterCode,iaLoopGasIndex(3),iFound,iMax
      CHARACTER caCldGasID(3)
      REAL      raSumCldAmt(3)
      REAL      cfrac   ! added ESM
      REAL      raCloudDME(kMaxCLouds)
      INTEGER   iactype_rtp(kMaxClouds)

      integer   ii      ! added ESM
      integer rtpopen, rtpread, rtpwrite, rtpclose, ifindj ! added ESM
      type(RTPHEAD) head
      type(RTPPROF) prof
      type(RTPATTR) hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      iWhichScatterCode = 6         !!RAYLEIGH in CLEAR SKY, nir/vis/uv
      iWhichScatterCode = 5         !!PCLSAM
      iWhichScatterCode = 4         !!r = r0 + r1 = perturb (not yet done)
      iWhichScatterCode = 3         !!DISORT
      iWhichScatterCode = 2         !!RTSPEC
      iWhichScatterCode = 1         !!TWOSTREAM
      iWhichScatterCode = 0         !!simple absorb; directly goes to rad_main

      iWhichScatterCode = 5         !!PCLSAM
      iDefault = 5

      IF (iDefault .NE. iWhichScatterCode) THEN
        print *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
      END IF

      IF (iWhichScatterCode .EQ. 6) THEN
        kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 5) THEN
        kWhichScatterCode = 5        !use PCLSAM
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 4) THEN
        kWhichScatterCode = 4        !use r = r0 + r1 = perturb
        kScatter          = 1        !
      ELSEIF (iWhichScatterCode .EQ. 3) THEN
        kWhichScatterCode = 3        !use Disort
        kScatter          = 1        !use this setting
        kDis_Pts          = 400      !do 1 every 400 pts
      ELSEIF (iWhichScatterCode .EQ. 2) THEN
        kWhichScatterCode = 2        !use RTSPEC
        kScatter          = 1        !use this setting  SingleScatter
        kScatter          = 3        !use this setting  Hybrid
        IF (kScatter .NE. 3) THEN
          write (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
        END IF
      ELSEIF (iWhichScatterCode .EQ. 1) THEN
        kWhichScatterCode = 1        !use TwoStream
        kScatter          = 1        !use one run of TwoStream
      ELSEIF (iWhichScatterCode .EQ. 0) THEN
        kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
        kScatter          = 1        !   
      END IF

      IF ((iakSolar(1) .GE. 0)  .AND. (kWhichScatterCode .EQ. 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
        CALL DoStop
      END IF

      IF ((iakSolar(1) .GE. 0)  .AND. (kWhichScatterCode .EQ. 4)) THEN
        write(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
        CALL DoStop
      END IF

       IF (iWhichScatterCode .EQ. 0) THEN
         write(kStdWarn,*) 'Purely absorptive cloud is set from rtp!!!!!'
         iScatBinaryFile         = iBinORasc
         iNClouds                = 1
         caaScatter(1)           = caaCloudFile(1)
         raaScatterPressure(1,1) = raCprTop(1)
         raaScatterPressure(1,2) = raCprBot(1)
         raScatterDME(1)         = raCpSize(1)
         raScatterIWP(1)         = raCngWat(1)
         write(kStdWarn,*) 'Cloud datafile  is : '
         write(kStdWarn,222) caaScatter(1)
         write(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1),
     $      raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
         GOTO 333
       END IF

      iScatBinaryFile         = iBinORasc

c <----------------------------------------------------------------------->
      ! now read the cloud parameters  --------------------------->

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen=+1
      write(kStdWarn,*)  'read open status = ', status

      DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status .EQ. -1) THEN
          write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdErr,*) fname
          CALL DoStop
        END IF
      END DO

      write (kStdWarn,*) 'success : read in RTP with gas profile ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status

      kProfileUnitOpen = -1

      !now see if there is a cloud to be used with this atmosphere
      IF (prof%cfrac .gt. 0.0) THEN
        cfrac   =  1.0            !assume total cloud cover
        cfrac   =  prof%cfrac     !be more honest
        cfrac1  =  prof%cfrac     !be more honest
        !!!first cloud is easy to do
        DO i = 1,1
          ctype1    = prof%ctype
          iaCtype(i) =  prof%ctype       !cloud type 1=cirrus 2=water etc
          raCemis(i) = 1.0               !assume cloud totally emissive
          raCngwat(i) = -9999            !IWP set by cloud profile
          raCpsize(i) = prof%cpsize*1.00 !in microns
          raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
          raCprbot(i) = prof%spres
        END DO

        cfrac2  =  prof%cfrac2     !be more honest
        DO i = 2,iNclouds_RTP
          ctype2    = prof%ctype2
          iaCtype(i) =  prof%ctype2       !cloud type 1=cirrus 2=water etc
          raCemis(i) = 1.0               !assume cloud totally emissive
          raCngwat(i) = -9999            !IWP set by cloud profile
c          raCpsize(i) = prof%udef(12)    !in microns
          raCpsize(i) = prof%cpsize2     !in microns
          raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
          raCprbot(i) = prof%spres
        END DO
       ELSE
        cfrac  =  0.0            !assume clear sky, use dummy values
        cfrac1 =  0.0
        cfrac2 =  0.0
        ctype1 = -101
        ctype2 = -101
        DO i = 1,kMaxClouds
          iaCtype(i) =  -101             !cloud type
          raCemis(i)  = 0.0              !assume cloud totally emissive
          raCngwat(i) = 0.0              !IWP
          raCpsize(i) = 1.0              !in microns
          raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
          raCprbot(i) = prof%spres
        END DO
      END IF

c <----------------------------------------------------------------------->
      !now read the actual cloud profile ------------------->

      fname(1:80) = caCloudPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      IF (status .eq. -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
      END IF
      kProfileUnitOpen = +1
      write(kStdWarn,*)  'read open status  =  ', status

      DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status .EQ. -1) THEN
          write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
          write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
          write(kStdErr,*) fname
          CALL DoStop
        END IF
      END DO
      write (kStdWarn,*) 'success : read in RTP with cloud profile ',iRTP
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status

      kProfileUnitOpen = -1
c <----------------------------------------------------------------------->

      IF (prof%plevs(1) .lt. prof%plevs(prof%nlevs)) THEN
        !layers are from TOA to the bottom
        iDownWard = -1
        kcRTP_pBot = prof%plevs(prof%nlevs)
        kcRTP_pTop = prof%plevs(1)
        kcRTPBot   = kProfLayer - (prof%nlevs-1) + 1
        kcRTPTop   = kProfLayer
      ELSE
        !layers are from GND to the top
        iDownWard = +1
        kcRTP_pTop = prof%plevs(prof%nlevs)
        kcRTP_pBot  = prof%plevs(1)
        kcRTPTop   = 1
        kcRTPBot   = prof%nlevs-1
      END IF

      iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
      IF (iProfileLayers .NE. iL1) THEN
        write(kStdErr,*) 'Oops : gas profiles have ',iProfileLayers,' layers'
        write(kStdErr,*) '     : cld profiles have ',iL1            ,' layers'
        CALL DoStop
      END IF

      iGasInRTPFile = head%ngas              !!! number of gases >=5, <= 7
  
      IF (prof%nlevs .GT. kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
      END IF
 
      write(kStdWarn,*) 'Reading profile from RTP file... '
      write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
      write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
      write(kStdWarn,*) '  top,bot = ',kcRTPBot,kcRTPTop,kcRTP_pBot,kcRTP_pTop

      !!!now check if this agrees with iL1,iGasInRTPFile above
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
      END IF
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
      END IF

      DO i = 1,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        kcraPressLevels(j) = prof%plevs(i)              !!!!in mb
      END DO

      DO i = 1,1
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        deltaP = kcraPressLevels(j) - raPressLevels(j)
        IF (abs(deltaP)/kcraPressLevels(j) .gt. 0.001) THEN
          write(kStdWarn,*) 'comparing pressure levels of gas,cld profiles'
          write(kStdWarn,*) 'TOA oops at i,j = ',i,j
          write(kStdWarn,*) 'kcraPressLevels,raPressLevels = ',
     $         kcraPressLevels(j),raPressLevels(j) 
        END IF
      END DO
      DO i = 2,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        deltaP = kcraPressLevels(j) - raPressLevels(j)
        IF (abs(deltaP)/kcraPressLevels(j) .gt. 0.001) THEN
          write(kStdErr,*) 'comparing pressure levels of gas,cld profiles'
          write(kStdErr,*) 'oops at i,j = ',i,j
          write(kStdErr,*) 'kcraPressLevels,raPressLevels = ',
     $         kcraPressLevels(j),raPressLevels(j) 
          CALL DoStop
        END IF
      END DO
        
      write(kStdWarn,*) 'checked input gas   profile pressure levels vs'
      write(kStdWarn,*) '        input cloud profile pressure levels ...'

c this variable keeps track of how many gases should be read in
c check that the gases are indeed 201,202,203

      iaCldGasID(1) = 201
      iaCldGasID(2) = 202
      iaCldGasID(3) = 203

      !! check for clouds
      iNclouds = 0
      DO j = 1,3
        iaLoopGasIndex(j) = -1
        DO i = 1,head%ngas
          IF (head%glist(i) .eq. iaCldGasID(j)) THEN
            write(kStdWarn,*) 'in  the rtp file, found cloudID = ',head%glist(i)
            iNclouds                 = iNclouds + 1
            iaCldTypes(iNclouds)     = iaCldGasID(j)
            iaLoopGasIndex(iNclouds) = i
            GOTO 10
          END IF
        END DO
 10     CONTINUE
      END DO

      IF (iNclouds .NE. iNclouds_RTP) THEN
        write(kStdErr,*) 'In SUBROUTINE READRTP_CLD100LAYER we have'
        write(kStdErr,*) 'iNclouds,iNclouds_RTP = ',iNclouds,iNclouds_RTP
        CALL DoStop
      END IF

c      prof%udef(11,:) = cngwat
c      prof%udef(12,:) = cpsize
c      prof%udef(13,:) = cprtop
c      prof%udef(14,:) = cprbot
c      prof%udef(15,:) = cfrac
c      prof%udef(16,:) = "cfrac12", fraction of FOV containing both clouds
c      prof%udef(17,:) = ctype {currently not used}
c SARTA has
c   000 - 099 = black   clouds
c   100 - 199 = water   clouds
c   200 - 299 = ice     clouds
c   300 - 399 = mineral clouds
c   400 - 499 = seasalt clouds
c   500 - 599 = soot/smoke         clouds
c   600 - 699 = sulphate/pollution clouds
c right now all kCARTA worries about is water, ice or mineral = (300-699) clds

      !! check that ctype and cpsizes match between nml and rtp info
      !! can only check upto 2 clouds from RTP file
      iactype_rtp(1) = prof%ctype
      IF (iNclouds .GE. 2) THEN
c        iactype_rtp(2) = prof%udef(17)
        iactype_rtp(2) = prof%ctype2
      END IF

      raCloudDME(1) = prof%cpsize
      raCloudDME(2) = prof%cpsize2

      iMax = min(iNclouds,2)
      DO j = 1,iMax
 
        print *,'yaya',j,iaCloudScatType(j),iactype_rtp(j)
        iG = int(iactype_rtp(j)*1.0/100.0)
        iH = iDiv(iaCloudScatType(j),100)

        IF ((iG .LT. 1) .OR. (iH .LT. 1)) THEN
          write(kStdErr,*) 'kCARTA cannot do "black clouds" here '
          write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',
     $    j,iaCloudScatType(j),iactype_rtp(j)
          CALL DoStop
        END IF

        IF (iG .NE. iH) THEN
          write(kStdErr,*) 'RTP and kCARTA have different cloud types here '
          write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',
     $    j,iaCloudScatType(j),iactype_rtp(j)
          CALL DoStop
        ELSE
          write(kStdWarn,*) 'j,iaCloudScatType(j),iactype_rtp(j),p.gas_ID(j) = ', j,iaCloudScatType(j),iactype_rtp(j),iaCldTypes(j)
        END IF
        
c        IF ((j .EQ. 1) .AND. (iG .NE. 1)) THEN
c          write(kStdErr,*) 'oh oh kCARTA assumes gas201 == water'
c          CALL DOStop
c        ELSEIF ((j .EQ. 2) .AND. (iG .NE. 2)) THEN
c          write(kStdErr,*) 'oh oh kCARTA assumes gas202 == ice'
c          CALL DOStop
c        ELSEIF ((j .EQ. 3) .AND. (iG .NE. 3)) THEN
c          write(kStdErr,*) 'oh oh kCARTA assumes gas203 == dust'
c          CALL DOStop
c        END IF

      END DO

      IF (iMax .LT. iNclouds) THEN
        write(kStdWarn,*) 'cannot check nml/rtp for cpsize,ctype beyond cld2'
      END IF

      iFound = iNclouds

      DO iG = 1,iNclouds
        iGasIndex = iG
        iCldIndex = iaCldTypes(iGasIndex) - iaCldGasID(1) + 1
        !!! first fill things out with stuff from the RTP file
        !!! loop till you find the gas
        iGX = iaLoopGasIndex(iG)
        DO i = 1, prof%nlevs - 1
          j  = iFindJ(kProfLayer,I,iDownWard)
          raaKlayersCldAmt(j,iGasIndex) = max(prof%gamnt(i,iGX),0.0)
        END DO      
        !!! then fill bottom of atm with zeros for cld amt
        DO i = prof%nlevs, kProfLayer
          j  = iFindJ(kProfLayer,I,iDownWard)
          raaKlayersCldAmt(j,iGasIndex) = 0.0
        END DO
      END DO

      DO i=1,3
        raSumCldAmt(i) = 0.0
      END DO
      write(kStdWarn,*) '  '
      write(kStdWarn,*) ' Lay    GasID=       GasID=      GasID='
      write(kStdWarn,*) '        ',(iaCldTypes(iG),iG=1,iNclouds)
      write(kStdWarn,*) '----------------------------------------'
      DO iI = 1,iFound
        DO i = kProfLayer-(prof%nlevs-1)+1,kProfLayer
          !write(kStdWarn,*) i,(raaKlayersCldAmt(i,iG),iG=1,iNclouds)
          raSumCldAmt(iI) =  raSumCldAmt(iI) + raaKlayersCldAmt(iI,1)
        END DO
      END DO
      write(kStdWarn,*) '----------------------------------------'
      write(kStdWarn,*) 'total cldamt in g/m2 = ',(raSumCldAmt(i),i=1,3)

      write(kStdWarn,*) ' '
      DO iI = 1,iNclouds_RTP
        raExp(iI)               = 0.0      !same amount in all layers
        iaCloudNumLayers(iI)    = iNumLayer
        !iaaScatTable(iI,1)      = iI
        !caaCloudName(iI)        = 'RTP cloud'
        !caaaScatTable(iI,1)     = caaCloudFile(iI)
        !!!these next params are actually set from the nm_profile combos
        !!!raaaCloudParams(iI,1,1) = raCngwat(iI)
        !!!raaaCloudParams(iI,1,2) = raCpsize(iI)
        IF (raCloudDME(iI) .LT. 0) THEN
          write(kStdErr,*) 'rtp file gives eff diam < 0 ',iI,raCloudDME(iI)
          CALL DoStop
        ELSE
          DO j = 1,kCloudLayers
            raaaCloudParams(iI,j,1) = raSumCldAmt(iI)
            raaaCloudParams(iI,j,2) = raCloudDME(iI)
          END DO
        END IF
         
        raaPCloudTop(iI,1)      = raCprtop(iI)
        raaPCloudBot(iI,1)      = raCprbot(iI)
        iaCloudNumAtm(iI)       = 1
        iaaCloudWhichAtm(iI,1)  = 1
        iaPhase(iI)             = -1       !default to HG phase function

        write (KStdWarn,222) 'cloud file    = ',caaCloudFile(iI)
        write (kStdWarn,*)   'cloud top     = ',raCprtop(iI),' mb'
        write (kStdWarn,*)   'cloud bot     = ',raCprbot(iI),' mb'
        write (kStdWarn,*)   'cloud IWP     = ',raSumCldAmt(iI),' gm m-2'
        write (kStdWarn,*)   'particle size = ',raCloudDME(iI),' um'
      END DO

      DO iG = 1,iNclouds
        IF ((iaCloudScatType(iG) .GE. 100) .AND. (iaCloudScatType(iG) .LE. 199)) caCldGasID(iG) = 'W'
        IF ((iaCloudScatType(iG) .GE. 200) .AND. (iaCloudScatType(iG) .LE. 299)) caCldGasID(iG) = 'I'
        IF ((iaCloudScatType(iG) .GE. 300) .AND. (iaCloudScatType(iG) .LE. 699)) caCldGasID(iG) = 'A'
      END DO

      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' CLD | TYPE | KLAYERSID | NMLID |           datafile'
      write(kStdWarn,*) ' ----------------------------------------------------'

c now set the auxiliary info as in SetRTPCloud      
      DO iG = 1,iNclouds
        iGasIndex = iG
        iCldIndex = iaCldTypes(iGasIndex) - iaCldGasID(1) + 1
        iMapIndex = iG
        iaCloudNumLayers(iGasIndex)   =  iProfileLayers      
        iaCloudNumAtm(iGasIndex)      =  1
        iaaCloudWhichAtm(iGasIndex,1) =  1
        iaPhase(iGasIndex)            = -1
        caaCloudName(iGasIndex)       = 'RTP cloud'
        caaaScatTable(iGasIndex,1)    = caaCloudFile(iMapIndex)
        raaPCloudTop(iGasIndex,1)     = raCprTop(iGasINdex)    !!!TOA
        raaPCloudBot(iGasIndex,1)     = raCprBot(iGasIndex)    !!!GND
        raExp(iGasIndex)              = 0.0           !!!default "same" 
        DO i = 1,iProfileLayers
          j = iaaRadLayer(1,i)
          iaaCloudWhichLayers(iGasIndex,i) = 
     $         iaaRadLayer(1,iProfileLayers)-iaaRadLayer(1,i) + 
     $         (kProfLayer - iProfileLayers + 1)
          iaaCloudWhichLayers(iGasIndex,i) = j
          iaaScatTable(iGasIndex,i) = iMapIndex
          raaaCloudParams(iG,j,1)   = raaKlayersCldAmt(j,iGasIndex)
          raaaCloudParams(iG,j,2)   = raCloudDME(iGasIndex)
        END DO
        write(kStdWarn,110) iG,caCldGasID(iCldIndex),iaCldTypes(iG),
     $          iaaScatTable(iGasIndex,1),caaaScatTable(iGasIndex,1)
      END DO
      write(kStdWarn,*) ' ---------------------------------------------'

 222  FORMAT(A80) 
 100  FORMAT('  ',I3,'    ',A1,'      ',I3,'      ',A80)
 110  FORMAT('  ',I3,'    ',A1,'      ',I3,'       ',I3,'   ',A80)
 333  CONTINUE

c! added ESM
      iNclouds2 = iNclouds   
      iNclouds3 = iNclouds   
      iProfileLayers = iProfileLayers2
      DO i = 1,prof%nlevs
        raPresslevels2(i) = raPresslevels(i)
      ENDDO
      
      RETURN
      END

c************************************************************************
c this subroutine adds on RefProfile 1,2,3,4,5,6 gas amounts for those gases NOT
c in the RTP profile (using h.ptype = 2 or even 1)
c the afgl profiles generated by /home/sergio/KCARTA/UTILITY/AFGLprofs.m
      SUBROUTINE AddOnAFGLProfile(iProfileNum,
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

      implicit none

      include '../INCLUDE/kcarta.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      INTEGER iNumberofGasesRead,iProfileNum
      INTEGER igasindex  ! added ESM
      REAL    raaHeight(kProfLayer,kGasStore)
      INTEGER iaWhichGasRead(kMaxGas),iNumGases,iaInputOrder(kMaxGas)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)

c local variables
      INTEGER iI,iJ,iaNeed(kMaxGas),iNewRead,iFound,iFoundX,iIDgas,iLay,iG0,iK
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer)  !!US Std layer ppress, amt
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer)  !!US Std layer press, temp
      CHARACTER*1 cY,cN

      cY = 'Y'
      cN = 'N'

      IF ((iProfileNum .LT. 1) .OR. (iProfileNum .GT. 6)) THEN
        write(kStdErr,*) 'Can only substitute profs from AFGL '
        write(kStdErr,*) '1=STD, 2=TRP, 3=MLS, 4=MLW, 5=SAS, 6=SAW'
        CALL DoStop
      END IF

      IF (iProfileNum .EQ. 1) THEN
        CALL AddOnStandardProfile(
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)
        RETURN
      END IF

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'read profiles for ',iNumberOfGasesRead, ' gases ...'
      write(kStdWarn,*) 'for AFLG profile ',iProfileNum

c this is the list of gases for which we need profiles
      iK = 0
      write(kStdWarn,*) '            count gasID  found?'
      write(kStdWarn,*) '-------------------------------'
      DO iI = 1,iNumGases
        IF (iaInputOrder(iI) .LT. 100) THEN
          IF (iaWhichGasRead(iI) .EQ. +1) THEN
             iK = iK + 1
c            write(kStdWarn,200) iI,iaInputOrder(iI),cY
            write(kStdWarn,200) iK,iI,cY
          ELSE
c            write(kStdWarn,200) iI,iaInputOrder(iI),cN
            write(kStdWarn,200) -1,iI,cN
          END IF
        END IF
      END DO

      !!! now do gasids 101,102,103
      iI = 1
      IF (iaInputOrder(iI) .EQ. 1) THEN
        IF (iaWhichGasRead(iI) .EQ. +1) THEN
          write(kStdWarn,200) iI,101,cY
          write(kStdWarn,200) iI,102,cY
          write(kStdWarn,200) iI,103,cY
        ELSE  
          write(kStdWarn,200) iI,101,cN
          write(kStdWarn,200) iI,102,cN
          write(kStdWarn,200) iI,103,cN
        END IF
      END IF
      write(kStdWarn,*) ' '

 200  FORMAT('xi, idgas = ',I5,I5,'  ',A1)

c this is the list of gases for which we have read in the profiles
c      DO iI = 1,kMaxGas
c        IF (iaWhichGasRead(iI) .EQ. +1) THEN
c          write(kStdWarn,*) 'RTP file had profile for gasID ',iI
c        END IF
c      END DO
c      write(kStdWarn,*) ' '

      iFound = -1
      DO iI = 1,kMaxGas
        IF (iaWhichGasRead(iI) .EQ. +1) THEN
          GOTO 10
        END IF
      END DO
 10   CONTINUE
      iG0 = iI     !!this gas filled out from RTP file, so use T(z),h(z)
      
c thus the difference between the lists is what we need AFGL  profiles for
      DO iI = 1,kMaxGas
        iaNeed(iI) = -1
      END DO

      iNewRead = 0
      DO iI = 1,iNumGases
        iFound = -1
        iIDgas = iaInputOrder(iI)
        IF (iaWhichGasRead(iIDgas) .EQ. +1) THEN
          iFound = +1
        END IF
        IF (iFound .LT. 0) THEN
          !! gas not found in rtp file, so need the US Std profile
          iNewRead = iNewRead + 1
          iaNeed(iNewRead) = iIDgas
          !!CALL get_us_std(iIDgas,raPX,raPPX,raTX,raQX)
          CALL getAFGL(iProfileNum,iIDgas,raPX,raPPX,raTX,raQX)
          Call FindIndexPosition(iIDgas,iNumGases,iaInputOrder,
     $                           iFoundX,iGasIndex)
          IF (iFoundX .GT. 0) THEN 
            !write(kStdWarn,4321) iIDGas,iLay,rAmt,rT,rP,rPP
            DO iLay = 1,kProfLayer
              raaAmt(iLay,iGasIndex)       = raQX(iLay)
              raaTemp(iLay,iGasIndex)      = raaTemp(iLay,iG0)
              raaPress(iLay,iGasIndex)     = raaPress(iLay,iG0)
              raaPartPress(iLay,iGasIndex) = raPPX(iLay)
              raaHeight(iLay,iGasIndex)    = raaHeight(iLay,iG0)
            END DO
            iaWhichGasRead(iIDgas)    = 1
          ELSE 
            write (kStdErr,*) 'huh? FindIndexPosition failed for ',iIDgas
            CALL DoStop
          END IF
        END IF
      END DO

      write(kStdWarn,*) 'Before entering "AddOnAFGLProfile" '
      write(kStdWarn,*) '  had read in profiles for ',iNumberofGasesRead
      write(kStdWarn,*) '  out of ',iNumGases
      write(kStdWarn,*) 'Read in ',iNewRead,' more in "AddOnAFGLProfile"'
 
      IF (iNumberofGasesRead + iNewRead .NE. iNumGases) THEN
        write(kStdErr,*) 'need iNumberofGasesRead + iNewRead = iNumGases'
        CALL DoStop
      END IF

      RETURN
      END

c************************************************************************
