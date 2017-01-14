c Copyright 2014
c University of Maryland Baltimore County 
c All Rights Reserved

c this file deals with reading in the user file
c main routines for PTHFIL are here 
c also, the REFERENCE profile reading subroutine is here

c************************************************************************
c this subroutine computes the temperatures at the pressure levels
c since Antartic surface pressures can be as low as 500 mb, CO2.N2O,CO mixing ratios were failing if iBot=20
c this corresponds to raPresslevls(20) = 596 mb
c so subr Get_Temp_Plevs (in n_pth_mix.f) and subr compute_co2_mixratio (in kcartamisc.f)
c both needed to have iBot tweaked to layer 25 or above, which corresponds to raPresslevls(25) = 496 mmb
      SUBROUTINE Get_Temp_Plevs(iProfileLayers,iaGases,raaTemp,raaPress,raThickness,
     $  raPressLevels,raTPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL raaPress(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iaGases(kMaxGas),iProfileLayers
c output
      REAL raTPressLevels(kProfLayer+1)           

c local
      INTEGER iI,iJ,iOffSet,iCO2_ind,iDefault,iInterpType
      REAL raT(kProfLayer),raP(kProfLayer),rX,rY,logP(kProfLayer),rPmin,rPmax,dx
      REAL grav0,grav,Re,Rd,raGeopotentialThick1(kProfLayer),raGeopotentialThick2(kProfLayer)
      REAL rInt,raTX(10),raPX(10),p2h

      iCO2_ind = 2                 !! assume we have CO2 in profile

      Re    = kPlanetRadius        !! Planet radius
      grav0 = kGravity             !! m/s2
      Rd    = kSpecificGasDryAir   !!  specific gas constant for dry air J K−1kg−1  

      IF (iaGases(2) .NE. 1) THEN
        write(kStdErr,*) 'No CO2 in your profile??????????? have to "check" T(z) some other way'
        iCO2_ind = 1     !!! hopefully there is water!!!!
      END IF

      iOffSet = kProfLayer-iProfileLayers

      DO iI = 1,kProfLayer
        if (raThickness(iI) .LT. 0) THEN
	  write(kSTdErr,*) iI,raThickness(iI)
	  call dostopmesg('rathickness < 0 in Get_Temp_Plevs$')
	end if
      END DO
      
      rPmin = +1.0e10
      rPmax = -1.0e10
      iI = 75   !!! make this really high up
      IF ((raaPress(iI,iCO2_ind) .GT. raaPress(iI+1,iCO2_ind)) .AND.
     $    (raaPress(iI+1,iCO2_ind) .GT. raaPress(iI+2,iCO2_ind))) THEN
        !! pressures decreasing with index; use CO2 temps and pressures
        DO iI = kProfLayer-iProfileLayers+1,kProfLayer
          iJ = iI-iOffSet
          raP(iJ) = raaPress(kProfLayer-iJ+1,iCO2_ind)*kAtm2mb
          logP(iJ) = log(raP(iJ))
          raT(iJ) = raaTemp(kProfLayer-iJ+1,iCO2_ind)
          IF (logP(iJ) .GT. rPmax) rPmax = logP(iJ)
          IF (logP(iJ) .LT. rPmin) rPmin = logP(iJ)
c           print *,'a',iI,iJ,kProfLayer-iJ+1,raP(iJ),raT(iJ)
        END DO
      ELSEIF 
     $((raaPress(iI,iCO2_ind) .LT. raaPress(iI+1,iCO2_ind)) .AND.
     $ (raaPress(iI+1,iCO2_ind) .LT. raaPress(iI+2,iCO2_ind))) THEN
        !! pressures increasing with index; use CO2 temps and pressures
        DO iI = kProfLayer-iProfileLayers+1,kProfLayer
          iJ = iI-iOffSet
          raP(iJ) = raaPress(iJ,iCO2_ind)*kAtm2mb
          logP(iJ) = log(raP(iJ))
          raT(iJ) = raaTemp(iJ,iCO2_ind)
          IF (logP(iJ) .GT. rPmax) rPmax = logP(iJ)
          IF (logP(iJ) .LT. rPmin) rPmin = logP(iJ)
c          print *,'b',iJ,raP(iJ),raT(iJ)
        END DO
      ELSE
        DO iJ = iI-2,iI+2
	  write(kStdErr,*) 'Get_Temp_Plevs ',iJ,iCO2_ind,raaPress(iJ,iCO2_ind)
	END DO
        write(kStdErr,*) 'In Get_Temp_Plevs, Pressures neither increasing nor decreasing'
        CALL DoStop
      END IF

      iDefault = +1       ! spline is default
      iInterpType = -1    ! linear
      iInterpType = +1    ! spline
      iInterpType = iaaOverrideDefault(2,8)
      IF (abs(iInterpType) .NE. 1) THEN
        write(kStdErr,*) 'need iInterpType = +1 or -1 not ',iInterpType
        CALL DoStop
      END IF		       
      
      IF ((iDefault .NE. iInterpType) .AND. (kOuterLoop .EQ. 1)) THEN
        write(kStdErr,*) 'in Get_Temp_Plevs, default is take layer temps and spline interp them to plevs'
        write(kStdErr,*) '                   we are   taking layer temps and linear interp them to plevs'
        write(kStdWarn,*) 'in Get_Temp_Plevs, default is take layer temps and spline interp them to plevs'
        write(kStdWarn,*) '                   we are   taking layer temps and linear interp them to plevs'
      END IF
      
      DO iI = 1,kProfLayer+1
        IF (raPressLevels(iI) .GT. 0) THEN
          IF ((log(raPressLevels(iI)) .GE. rPmin) .AND. (log(raPressLevels(iI)) .LE. rPmax)) THEN
	    !! most of the points
	    IF (iInterpType .EQ. +1) THEN
              CALL rspl1(logP,raT,iProfileLayers,log(raPressLevels(iI)),rY,1)
c              CALL rspl_diffyp1n(logP,raT,iProfileLayers,log(raPressLevels(iI)),rY,1)	      
	    ELSE
              CALL rlinear1(logP,raT,iProfileLayers,log(raPressLevels(iI)),rY,1)
	    END IF
          ELSE
	    !! couple or so points at the top or bottom boundaries
            CALL rlinear1(logP,raT,iProfileLayers,log(raPressLevels(iI)),rY,1)
          END IF
          raTPressLevels(iI) = rY
        ELSE
          raTPressLevels(iI) = -9999
        END IF
      END DO

      write(kStdWarn,*) 'iI    pLower    pLayAvg     pUpper  |  tLower   tLayAvg     tUpper'
      write(kStdWarn,*) '-----------------------------------------------------------------------'
      DO iI = kProfLayer-iProfileLayers+1,kProfLayer
        iJ = iI-iOffSet
	iJ = iI
        write(kStdWarn,2345) iI,raPressLevels(iI),raP(kProfLayer-iJ+1),raPressLevels(iI+1),
     $                          raTPressLevels(iI),raT(kProfLayer-iJ+1),raTPressLevels(iI+1)
      END DO
      write(kStdWarn,*) '-----------------------------------------------------------------------'      
      
 1234 FORMAT(I3,5(' ',F10.4))
 2345 FORMAT(I3,6(' ',F10.4)) 

c >>>>>>>>>>> this is to FORCE tape5 temperatures into the LevelTemperatures and LayerAvgTemps
      IF ((kRTP .EQ. -5) .OR. (kRTP .EQ. -6))THEN
        write(kStdWarn,*) 'kRTP = -5,-6 so resetting raTPressLevels (levelsT)    with TAPE5 levels info'
        write(kStdWarn,*) 'kRTP = -5,-6 so resetting raVTemp        (layersTavg) with TAPE5 levels info'
	write(kStdWarn,*) '     iI          raTPressLevels(iI)            kLBLRTM_levelT(iI)'
	write(kStdWarn,*) '-----------------------------------------------------------------'	
	DO iI = 1,kProflayer+1
	  write(kStdWarn,*) iI,raTPressLevels(iI),kLBLRTM_levelT(iI),raTPressLevels(iI)-kLBLRTM_levelT(iI)
          raTPressLevels(iI) = kLBLRTM_levelT(iI)
        END DO
	write(kStdWarn,*) '-----------------------------------------------------------------'
	DO iI = 1,kProflayer
	  write(kStdWarn,*) iI,raaTemp(iI,1),kLBLRTM_layerTavg(iI),raaTemp(iI,1)-kLBLRTM_layerTavg(iI)
	  IF (kRTP .EQ. -5) THEN
  	    write(kStdWarn,*) 'kRTP = -5, use these T in uncompress',iI,raaTemp(iI,1),kLBLRTM_layerTavg(iI),
     $	    raaTemp(iI,1)-kLBLRTM_layerTavg(iI)
	  END IF
c          raVTemp(iI) = kLBLRTM_layerTavg(iI)
        END DO	
      END IF
      
c      DO iI = kProfLayer-iProfileLayers+1,kProfLayer
c        write(*,1234) iI,raPressLevels(iI),raP(kProfLayer-iI+1),raPressLevels(iI+1),
c     $             raT(kProfLayer-iI+1),raTPressLevels(iI)
c      END DO
c      iI = kProfLayer + 1
c      write(*,1234) iI,raPressLevels(iI),raPressLevels(iI),raPressLevels(iI),
c     $             raT(1),raTPressLevels(iI)
cc      DO iI = kProfLayer-iProfileLayers+1,kProfLayer+1
cc        print *,iI,iProfileLayers,raPressLevels(iI),raTPressLevels(iI)
cc      END DO
c      call dostop
      
c hyposometric equation
c now find geopotential thickness = Rd/g integral_p1^p2 T(p) d(ln(p)) = 
c   Rd/g 0.5(T1+T2) (ln(p2)-ln(p1)) = Rd/g 0.5(T1+T2) ln(p2/p1)
c
c we know mg = GMm/(Re+h)^2 = GM/Re^2 (1+h/Re)^-2  = grav(0) (1-2h/Re)
c method 1 : just use slab endpoints
      DO iI = kProfLayer-iProfileLayers+1,kProfLayer
        rInt = 0.5*(raTPressLevels(iI)+raTPressLevels(iI+1))
        rInt = rInt*log(raPressLevels(iI+1)/raPressLevels(iI))
        grav = grav0 * (1 - 2*p2h(raPressLevels(iI))/1000/Re)
        raGeopotentialThick1(iI) = Rd/grav * abs(rInt)
c        print *,iI,raThickness(iI),raGeopotentialThick1(iI),
c     $             (1-raGeopotentialThick1(iI)/raThickness(iI))*100
      END DO
c method 2 : use 10 points including slab ends

      write(kStdWarn,*) 'geopotential height vs layer thickness'
      write(kStdWarn,*) ' iI      p2h      laythick    GeoHgt1       GeoHgt2      Error1   Error2'
      write(kStdWarn,*) '------------------------------------------------------------------------'

      DO iI = kProfLayer-iProfileLayers+1,kProfLayer
        raTX(01) = raTPressLevels(iI)
        raTX(10) = raTPressLevels(iI+1)
        raPX(01) = raPressLevels(iI)
        raPX(10) = raPressLevels(iI+1)
        dx = log(raPressLevels(iI+1)/raPressLevels(iI)) * 1/9  !! pressure spacing in (ln(p)) space
        DO iJ = 2,9
          raPX(iJ) = exp(log(raPX(01)) + (iJ-1)*dx)
          CALL rspl1(logP,raT,iProfileLayers,log(raPX(iJ)),rY,1)
          raTX(iJ) = rY
        END DO
c        DO iJ = 1,10
c          print *,iI,iJ,raPX(01),raTX(01),raPX(iJ),raTX(iJ)
c        END DO
        raGeopotentialThick2(iI) = 0.0
        grav = grav0 * (1 - 2*p2h(raPressLevels(iI))/1000/Re)
        DO iJ = 1,9        
          rInt = 0.5*(raTX(iJ)+raTX(iJ+1))
          rInt = rInt*log(raPX(iJ+1)/raPX(iJ))
          raGeopotentialThick2(iI) = raGeopotentialThick2(iI) + Rd/grav * abs(rInt)
        END DO
        write(kStdWarn,111) iI,p2h(raPressLevels(iI)),abs(raThickness(iI)),
     $             raGeopotentialThick1(iI),raGeopotentialThick2(iI),
     $             (1-raGeopotentialThick1(iI)/abs(raThickness(iI)))*100,
     $             (1-raGeopotentialThick2(iI)/abs(raThickness(iI)))*100
      END DO
 111  FORMAT(I3,6('  ',F10.4))
 
c GeoHgt1/2 and Error1/2 shows Scotts klayers and my hypersometric equation, are consistent

      write(kStdWarn,*) '------------------------------------------------------------'
      write(kStdWarn,*) '  '

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword
c this differs from GENLN2 format in that
c (1) instead of veloctiy, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102,103 so 
c things have to be done slightly differently
c also added on gas 103 (heavy water isotope 4  HOD = 162)
      SUBROUTINE readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,raPressLevels,raThickness)

      IMPLICIT NONE

      INTEGER iPLEV
      
      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/KCARTA_database.param'

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

c local variables for our klayers files
c READ (caStr,*) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      REAL raaHeight(kProfLayer,kGasStore)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iIOUN2,iNpath,iaGasPathCounter(kProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

      INTEGER iL1,iL2,length130,iI,iDefault,iReadP
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
c note we use WATER amts for self and forn continuum, as well as heavy water
c so be careful
      DO iIDGas = kNewGasLo,kNewGasHi+1
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
        iaGasPathCounter(iNpath)=0
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

      iIOUN2 = kProfileUnit
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
        write (kStdWarn,*) ' '
        write (kStdWarn,*) 'KLAYERS profile has ',iL2,' gases in atm'
        write (kStdWarn,*) 'KLAYERS profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'so input text file has ',iNpath, ' paths'
        write (kStdWarn,*) ' '
        write (kStdWarn,*) '  we need paths for ',iNeed2Read,' gases'
        write (kStdWarn,*) '  contained inside ',kProfLayer,' layers'
        write (kStdWarn,*) ' '
      END IF

      IF (kProfLayer .NE. iL1) THEN
        write (kStdErr,*) 'KLAYERS profile has ',iL2,' gases in atm'
        write (kStdErr,*) 'KLAYERS profile has ',iL1,' layers in atm'
        write (kStdErr,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdErr,*) 'for reading in old GENLN2 profs, Need kProfLayer == Layers in Profile'
        write (kStdErr,*) 'Either redo klayers.x or kcarta.x'
        CALL DoStop
      END IF

      IF (iNeed2Read*kProfLayer .GT. iNpath) THEN
        write(kStdWarn,*) 'kCARTA/kLAYERS mismatch!!!!'
        write(kStdWarn,*) 'Not enough gas paths in your supplied profile'
        write(kStdWarn,*) '  may need to add on gases!'
        iNeedMoreProfiles = 1
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
         iaGasPathCounter(iIDgas)=iaGasPathCounter(iIDgas)+1

         IF (rAmt .lt. 0.0) THEN
           WRITE(kStdErr,1080)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rT .lt. 0.0) THEN
           WRITE(kStdErr,1081)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rP .lt. 0.0) THEN
           WRITE(kStdErr,1082)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rPP .lt. 0.0) THEN
           WRITE(kStdErr,1083)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
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
           raaAmt(iaGasPathCounter(iIDgas),iGasIndex)=rAmt
           raaTemp(iaGasPathCounter(iIDgas),iGasIndex)=rT
           raaPress(iaGasPathCounter(iIDgas),iGasIndex)=rP
           raaPartPress(iaGasPathCounter(iIDgas),iGasIndex)=rPP
           raaHeight(iaGasPathCounter(iIDgas),iGasIndex)=rH*1000 !!change to m!

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
      IF (iaGasPathCounter(iIDgas) .LT. kProfLayer) THEN
        GOTO 35
      END IF

      WRITE(kStdWarn,4000) iaGasPathCounter(iIDgas),iIDgas
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
      DO iIDGas = kNewGasLo,kNewGasHi+1
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
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'your profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        IF (iNeedMoreProfiles .EQ. 1) THEN
          write(kStdWarn,*) 'adding on AFGL profile ',kAFGLProf, ' for remaining gases'
          CALL AddOnAFGLProfile(kAFGLProf,
     $         iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $         raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
        ELSE
          !!this is just debugging, and/or to stop, just like KCARTAv1.12-
          write(kStdErr,*) 'your profile did not have all the gases'
          write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
          write(kStdErr,*) 'did not "add" profiles'
          CALL DoStop
        END IF
      END IF

 5030 FORMAT(A130)

c now set raLayerHeight
      DO iFound=1,kProfLayer
        raLayerHeight(iFound)=raaHeight(iFound,1)
      END DO

c now reread the profile file, so that we can get info about presslevels
c and layer thicknesses
c just turn this off to read the real old klayers files 
c -->>> (eg that Scott uses for his kcartav1.07 runs ....) <<<---
      iDefault = +1
      iReadP = +1    !! assume GENN2 style profile has p info
      iReadP = -1    !! assume GENN2 style profile does not have p info
      iReadP = iaaOverrideDefault(1,6)
      IF (abs(iReadP) .NE. 1) THEN
        write(kStdErr,*) 'invalid iReadP = ',iReadP
        CALL DoStop
      END IF		       
      
      IF (iDefault .EQ. iReadP) THEN
        CALL GetMoreInfo(raPressLevels,raThickness,caPfName)
      ELSE
        write(kStdErr,*) 'have turned off GetMoreInfo to read Scott Hannon profiles'
	DO iFound = 1,kProfLayer+1
	  raPressLevels(iFound) = PLEV_KCARTADATABASE_AIRS(iFound)
	END DO
      ENDIF

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
      CHARACTER c1a,c1b,cDumb
      CHARACTER*10 ca10a,ca10b,ca10c,ca10d,ca10e,ca10f
      CHARACTER*3 c3
      CHARACTER*2 c2
      INTEGER iI,iJ,iNumP,iNumH,iTrue,iK,iL
      REAL p1,p2,raP1(kProfLayer+1),raP2(kProfLayer+1)
      REAL raX1(kProfLayer+1),raX2(kProfLayer+1)
      REAL h1,raH1(kProfLayer),r1,r2,r3,r4,r5
      INTEGER iProfileLayers,iPrintfile

      DO iI=1,5
        ca10a(iI:iI) = ' '
        ca10b(iI:iI) = ' '
        ca10c(iI:iI) = ' '
        ca10d(iI:iI) = ' '
        ca10f(iI:iI) = ' '
      END DO
      ca10a(1:4) = 'PATH'
      ca10b(1:3) = 'km,'
      ca10c(1:3) = ' - '
      ca10d(1:3) = 'RAY'
      ca10e(1:1) = '='
      ca10f(1:8) = 'HEIGHTS:'

      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
          CALL DoSTOP
      ENDIF
      kProfileUnitOpen=1
 1070 FORMAT('ERROR! number ',I5,' opening PRFILE path profile file'
     $            ,/,A80)

 25   FORMAT(F8.3,A2,F8.3) 
      iNumP = 0
      iNumH = 0
 40   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr

      IF ((iNumP .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        !!! find PATH
        iS = SubString(caStr,ca10a,3,6)
        IF (iS .EQ. 1) THEN
          iNumP = iNumP + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            caTemp2(iI:iI) = ' '
          END DO

          !!! find HEIGHTS
          iI = 8
 5        CONTINUE
          iS = SubString(caStr,ca10f,iI,iI+7)
          IF (iS .EQ. -1) THEN
            iI = iI+1
            GOTO 5
          END IF
c differs from GetMoreInfo_GENLN4 becuz of this 
c          iI = iI + 8
          iI = iI + 7 + 8
          iJ = 130-iI
          iL = iI
          DO iK = 1,iJ
            caTemp(iK:iK) = caStr(iL:iL)
            iL=iL+1
          END DO
          read(caTemp,25) p1,c2,p2
          raX1(iNumP) = p1
          raX2(iNumP) = p2

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

      iProfileLayers = kProfLayer

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

      raX1(KProfLayer+1) = raX2(kProfLayer)
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
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',
     $ raP1(1),raP1(2)
      write (kStdWarn,*) 'Highest layer : press levels (mb) = ',
     $ raP1(kProfLayer),raP1(kProfLayer+1)
      write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(1),raH1(2)
      write (kStdWarn,*) '2 Highest layers thickness (km) = ',
     $  raH1(kProfLayer-1),raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
c see ../INCLUDE/KCARTA_database.param for the kCARTA database definitions
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(1)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)


c this is to dump 
c 101 pressure levels heights for the DEFINITION of the kCARTA database
c this comes from KCARTADATA/RefProf.For.v107up/mypref.op.new
c run bkcarta.x with the profile to be read in being mypref.op.new
       iPrintfile = -1
       IF (iPrintfile .GT. 0) THEN
         print *,'c from KCARTADATA/RefProf.For.v107up/mypref.op.new'
         print *,'c the heights of the pressure levels for the kCARTA DATABASE'
         print *,'      REAL DATABASELEVHEIGHTS(101)'
         print *,' '
         print *,'C-----------------------------------------------------------'
         print *,'PLEV with the AIRS layer boundary pressure heights (in km)'
         print *,'C-----------------------------------------------------------'
         print *,'      DATA  (DATABASELEVHEIGHTS(IPLEV), IPLEV = 101, 1, -1 )'
         DO iI = 1,20
           iK = kProfLayer - (iI-1)*5 + 1
           iL = iK - 5 + 1
           r1 = raX1(iK)
           r2 = raX1(iK-1)
           r3 = raX1(iK-2)
           r4 = raX1(iK-3)
           r5 = raX1(iK-4)
           IF (iI .EQ. 1) THEN
             write (*,5001) r1,r2,r3,r4,r5
           ELSE
             write (*,5002) r1,r2,r3,r4,r5
           END IF
         END DO
        iJ = 1
        write (*,5003) raX1(iJ)
      END IF

 5001 FORMAT ('     $   / ',5(F11.5 ,', '))
 5002 FORMAT ('     $     ',5(F11.5 ,', '))
 5003 FORMAT ('     $     ',F11.5, '/')

      RETURN
      END

c************************************************************************
c this function checks to see if substring ca10 is in string caX(i1:i2)
c if found, it returns +1, else returns -1
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
c for either the lower or the upper atm
c it flags an error if kProfLayers layers are not read in
c ProX === A=amount,T=temperature,P=pressure,PP=partial pressure
       SUBROUTINE ReadRefProf(caFname,iNlayIn,raProA,raProT,
     $      raProP,raProPP,iErrOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caFName   = name of file that has the profile
c iNlay     = number of layers that are read in
c raProA/T/P/PP = profile amout, temperature,pressure,partial pressure
c iErrOut   = error count (usually associated with file I/O)
       CHARACTER*80 caFname
       INTEGER iNlayIn,iErrOut 
       REAL raProA(*),raProT(*)
       REAL raProP(*),raProPP(*)

c local variables
       INTEGER iErr,iJ,iNlay,iIOUN
       CHARACTER*100 caLine

       iIOUN = kProfileUnit

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
c iNumLayers = number of layers per gas profile ( = kProfLayer)
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
        WRITE(kStdErr,7010) iNpmix
        WRITE(kStdErr,7011) kMixFilRows
        write(kStdErr,*)'Reset iNpmix (# of 100 chunks in mixfile)'
        write(kStdErr,*)'or kMixFilRows, in file kcarta.param '
        CALL DoSTOP
      END IF

 7010 FORMAT('iNpmix*kProfLayer = ',I4,' is greater than max allowed')
 7011 FORMAT('(from kcarta.param, kMixFilRows = ',I4,')')

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
c      END DO
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
          raaMix(iMPReadIn*iNumLayers+iLay,iGas) = raStore(iaInputOrder(iGas))
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
c this subroutine reads stuff from caaM, and stores the results either in
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
c this subroutine deals with 'PTHFIL' keyword
c this differs from GENLN2 format in that
c (1) instead of veloctiy, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently
      SUBROUTINE readGENLN4LAYERS(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,raPressLevels,raThickness,raTPressLevels,
     $      iNumLayers)

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
c iNumLayers = how many layers this profile set puts in (<= kProfLayers)
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases,iNumLayers
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

c local variables for our klayers files
c READ (caStr,*) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      REAL raaHeight(kProfLayer,kGasStore)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr,caStrShape
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iIOUN2,iNpath,iaGasPathCounter(kProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)
c extra variables for GENLN4 "layers" files
c READ (caStr,*) iMfil,iIDgas,iISO,rAmt,rT,rTU,rTL,rP,rPT,rPB,rPP,rV
      REAL rTU,rTL,rPT,rPB,rV  !!upper,lower level temp and pressure, velocity
      INTEGER iMFil,iISO,iInFile,iPathOffset
      CHARACTER*80 caLine
      CHARACTER*6 ex,path
      CHARACTER*9 heights,gas,gas2,gas3
      CHARACTER*4 star,star2,km,mb
      CHARACTER*2 d1,d2
      CHARACTER*1 c1
      REAL h1,h2,p1,p2
      CHARACTER*10 caLineShape,caCon

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

      INTEGER iL1,iL2,iL3,iL4,iL5,iL6,iL7,length130,iIinFile,iI
      CHARACTER*130 ca1,ca2,ca3,ca4,ca5,ca6,ca7,caTemp
      INTEGER iNumPathsPerSet,iNumSet

      DO iI = 1,kProfLayer
        raPressLevels(iI) = 0.0
        raThickness(iI) = 0.0
      END DO
      raPressLevels(kProfLayer+1) = 0.0

      ca1 = '! TOTAL NO. OF LAYERS IN ATM:'
      ca2 = '! TOTAL NO. OF GASES:'
      ca3 = '! PATH'
      ca4 = '! RAY LENGTH'
      ca5 = '! SET'
      ca6 = '!***********'
      ca7 = '! TOTAL NO. OF PATHS FOR EACH SET, NO. OF SETS:'
      iL1 = length130(ca1) 
      iL2 = length130(ca2) 
      iL3 = length130(ca3) 
      iL4 = length130(ca4) 
      iL5 = length130(ca5) 
      iL6 = length130(ca6) 
      iL7 = length130(ca7) 

c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn=0

c this variable keeps track of how many gases should be read in
      iNeed2Read=iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi+1
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
c set all individual gas paths to the offset 
      DO iNpath=1,kProfLayer
        iaGasPathCounter(iNpath)=0
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

      iIOUN2 = kProfileUnit
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

      !!!dumb settings
      iL1 = iNumPathsPerSet
      iL2 = iNumSet

 30   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr
      CALL rightpad130(caStr)

      IF (caStr(1:iL7) .EQ. ca7) THEN
        READ (iIOUN2,5030,ERR=13,END=13) caStr
        CALL rightpad130(caStr)
        READ (caStr,*) c1,iNumPathsPerSet,iNumSet
        iL1 = iNumPathsPerSet
        iL2 = iNumSet
        GOTO 30
      END IF

      IF (caStr(1:iL1) .EQ. ca1) THEN
        DO iP=1,130
          caTemp(iP:iP) = ' '
        END DO
        caTemp(1:10) = caStr(iL1+1:iL1+10)
        read (caTemp,*) iL1         !!!now we know how many LAYERS in profile
      END IF

      IF (caStr(1:iL2) .EQ. ca2) THEN
        DO iP=1,130
          caTemp(iP:iP) = ' '
        END DO
        caTemp(1:10) = caStr(iL2+1:iL2+10)
        read (caTemp,*) iL2         !!!now we know how many gases in profile
      END IF

      IF (iNumPathsPerSet .NE. iL1) THEN
        write(kStdErr,*) 'doofus : need iNumPathsPerSet = iL1'
        write(kStdErr,*) iNumPathsPerSet,iL1
        CALL DoStop
      END IF

      IF (iNumSet .NE. iL2) THEN
        write(kStdErr,*) 'doofus : need iNumSet = iL2'
        write(kStdErr,*) iNumSet,iL2
        CALL DoStop
      END IF

      IF (caStr(1:1) .EQ. '!') THEN
        GO TO 30
      ELSE
        iNumLinesRead=iNumLinesRead+1
        READ (caStr,*,ERR=13,END=13) iNpath
        WRITE(kStdWarn,3000) iNpath,iNeed2Read*kProfLayer
 3000   FORMAT('input file has ',I5,' paths; kCARTA needs ',I5,' paths')
        IF (kProfLayer .EQ. iNumPathsPerSet) THEN
          iGasesInProfile = INT(iNpath*1.0/(kProfLayer*1.0))
        ELSEIF (kProfLayer .GT. iNumPathsPerSet) THEN
          iGasesInProfile = INT(iNpath*1.0/(iNumPathsPerSet*1.0))
          write(kStdWarn,*) 'warning iNumPathsPerSet < kProfLayer'
        ELSEIF (kProfLayer .LT. iNumPathsPerSet) THEN
          write(kStdErr,*) 'OOPS iNumPathsPerSet > kProfLayer'
          CALL DoStop
        END IF
      END IF

      iPathOffset = +1    !!assume kProfLayer == iL1 == no. of layers in prof 
      !!!now check if this agrees with iL1,iL2 above
      IF (kProfLayer .LT. iL1) THEN
        iPathOffset = -1
        write (kStdErr,*) 'GENLN4 LAYERS profile has ',iL2,' gases in atm'
        write (kStdErr,*) 'GENLN4 LAYERS profile has ',iL1,' layers in atm'
        write (kStdErr,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdErr,*) 'Need kProfLayer >= Layers in Profile'
        write (kStdErr,*) 'Either redo g4layers.x or kcarta.x'
        CALL DoStop
      END IF

      write(kStdWarn,*) ' '
      IF (iNeed2Read*kProfLayer .GT. iNpath) THEN
        write(kStdWarn,*)'iGases2BeRead, kProfLayer = ',iNeed2Read,kProfLayer
        write(kStdWarn,*)'  where iGases2BeRead=num gases, modulo 101/102/103'
        write(kStdWarn,*)'Num of paths in your profile iNpath = ',iNpath
        write(kStdWarn,*)'kCARTA/GENLN4 LAYERS mismatch!!!!'
        write(kStdWarn,*)'Mebbe too few gas paths in your supplied profile'
      ELSE
        write(kStdWarn,*)'iGases2BeRead, kProfLayer = ',iNeed2Read,kProfLayer
        write(kStdWarn,*)'where iGases2BeRead=num of gases, modulo 101/102/103'
        write(kStdWarn,*)'Num of paths in your profile iNpath = ',iNpath
      END IF

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Going to start reading genln4 layers .... '
      write(kStdWarn,*) ' '
      write(kStdWarn,*) '  number of layers, gases in profile = ',iL1,iL2
      write(kStdWarn,*) '  number of paths in profile = ',iNpath

      READ (iIOUN2,5030,ERR=13,END=13) caStr   !!!read the lonely "!"

      iPathOffset = kProfLayer-iL1 !!!offset with which to store path info 
      iNumLayers =  kProfLayer - iPathOffSet  !!! === iL1
      IF (iPathOffSet .GT. 0) THEN
        write(kStdWarn,*) 'Need to offset GENLN4 profile by ',iPathOffSet
        write(kStdWarn,*) 'so as to stuff it into kProfLayer'
        write(kStdWarn,*) 'Number of profile layers = ',iNumLayers
      ELSEIF (iPathOffSet .LT. 0) THEN
        write(kStdErr,*) 'Need to offset GENLN4 profile by ',iPathOffSet
        write(kStdWarn,*) 'so as to stuff it into kProfLayer! IMPOSSIBLE!'
        CALL DoStop
      ENDIF

c set all individual gas paths to the offset 
      DO iI = 1,kProfLayer
        iaGasPathCounter(iI) = iPathOffSet
      END DO


c now loop iNpath/iNumGases  times for each gas in the user supplied profile
 35   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr
      CALL rightpad130(caStr)
      IF (caStr(1:iL3) .EQ. ca3) THEN
        !!!this has the hU,hB, pU,pB info, so read this
        READ(caStr,*) ex,path,iIinFile,gas,gas2,gas3,heights,
     $                h1,d1,h2,km,p1,d2,p2,mb
c        print *,ex,path,iIinFile,gas,gas2,gas3,heights,d1,km,d2,mb
c        print *,h1,h2,p1,p2
        iNumLinesRead=iNumLinesRead+1
        GO TO 35
      ELSEIF (caStr(1:iL5) .EQ. ca5) THEN
        !!!this has "set number" info, discard and go on
        GO TO 35
      ELSEIF (caStr(1:iL6) .EQ. ca6) THEN
        !!!this has "***********" info, discard and go on
        GO TO 35
      ELSEIF (caStr(1:iL4) .EQ. ca4) THEN
        !!!this has the ray length, theta info, so discard and read stuff
         iNumLinesRead=iNumLinesRead+1
         !!! since this would be too large a string to read,
         !!!directly read in the info
         READ (iIOUN2,*,ERR=13,END=13) iMfil,iIDgas,iISO,rAmt,rT,rTU,rTL,
     $                                rP,rPT,rPB,rPP
         READ (iIOUN2,5030,ERR=13,END=13) caStrShape         
         READ (caStrShape,*,ERR=13,END=13) caLineShape,caCon
         iaGasPathCounter(iIDgas) = iaGasPathCounter(iIDgas)+1
c         print *,iaGasPathCounter(iIDgas),iIDgas,iISO,rAmt,rT,rP,rPP

         IF (rAmt .lt. 0.0) THEN
           WRITE(kStdErr,1080)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rT .lt. 0.0) THEN
           WRITE(kStdErr,1081)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rP .lt. 0.0) THEN
           WRITE(kStdErr,1082)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

         IF (rPP .lt. 0.0) THEN
           WRITE(kStdErr,1083)
           WRITE(kStdErr,1111) iIDgas,iaGasPathCounter(iIDgas)
           CALL DoSTOP
           iErr=1
         END IF

       END IF

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
           raTPressLevels(iaGasPathCounter(iIDgas))   = rTU
           raTPressLevels(iaGasPathCounter(iIDgas)+1) = rTL
           raaAmt(iaGasPathCounter(iIDgas),iGasIndex)=rAmt
           raaTemp(iaGasPathCounter(iIDgas),iGasIndex)=rT
           raaPress(iaGasPathCounter(iIDgas),iGasIndex)=rP
           raaPartPress(iaGasPathCounter(iIDgas),iGasIndex)=rPP
           raaHeight(iaGasPathCounter(iIDgas),iGasIndex)=(h1+h2)/2*1000  !!in m

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
      IF (iaGasPathCounter(iIDgas) .LT. kProfLayer) THEN
        GOTO 35
      END IF

      WRITE(kStdWarn,4000) iaGasPathCounter(iIDgas)-iPathOffSet,iIDgas
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
      DO iIDGas = kNewGasLo,kNewGasHi+1
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
      CALL GetMoreInfo_Genln4(raPressLevels,raThickness,caPfName,iPathOffSet)

 1111 FORMAT('gasID, layer = ',I5,I5)
 1080 FORMAT('negative gas amount in PRFILE profile file')
 1081 FORMAT('negative gas temp in PRFILE profile file')
 1082 FORMAT('negative layer pressure in PRFILE profile file')
 1083 FORMAT('negative gas partial press in PRFILE profile file')

      RETURN
      END

c************************************************************************
c this subroutine reads the info from the GENLN4 LAYERS profile, storing info
c about layer thicknesses and presslevels
      SUBROUTINE GetMoreInfo_GENLN4(raPressLevels,raThickness,
     $                              caPfName,iPathOffSet)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c parameters
      CHARACTER*80 caPFname
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iPathOffSet

c local variables      
      INTEGER iIOUN2,iErrIO,iErr,SubString,iS
      CHARACTER*130 caTemp,caTemp2,caStr
      CHARACTER c1a,c1b,cDumb
      CHARACTER*10 ca10a,ca10b,ca10c,ca10d,ca10e,ca10f
      CHARACTER*3 c3
      CHARACTER*2 c2
      INTEGER iI,iJ,iNumP,iNumH,iTrue,iK,iL
      REAL p1,p2,raP1(kProfLayer+1),raP2(kProfLayer+1)
      REAL raX1(kProfLayer+1),raX2(kProfLayer+1)
      REAL h1,raH1(kProfLayer),r1,r2,r3,r4,r5
      INTEGER iProfileLayers,iPrintfile

      DO iI=1,5
        ca10a(iI:iI) = ' '
        ca10b(iI:iI) = ' '
        ca10c(iI:iI) = ' '
        ca10d(iI:iI) = ' '
        ca10f(iI:iI) = ' '
      END DO
      ca10a(1:4) = 'PATH'
      ca10b(1:3) = 'km,'
      ca10c(1:3) = ' - '
      ca10d(1:3) = 'RAY'
      ca10e(1:1) = '='
      ca10f(1:8) = 'HEIGHTS:'

      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
          CALL DoSTOP
      ENDIF
      kProfileUnitOpen=1
 1070 FORMAT('ERROR! number ',I5,' opening PRFILE path profile file'
     $            ,/,A80)

 25   FORMAT(F8.3,A2,F8.3) 
      iNumP = iPathOffSet
      iNumH = iPathOffSet
 40   CONTINUE
      READ (iIOUN2,5030,ERR=13,END=13) caStr

      IF ((iNumP .LT. kProfLayer) .AND. (caStr(1:1) .EQ. '!')) THEN
        !!! find PATH
        iS = SubString(caStr,ca10a,3,6)
        IF (iS .EQ. 1) THEN
          iNumP = iNumP + 1
          DO iI=1,130
            caTemp(iI:iI) = ' '
            caTemp2(iI:iI) = ' '
          END DO

          !!! find HEIGHTS
          iI = 8
 5        CONTINUE
          iS = SubString(caStr,ca10f,iI,iI+7)
          IF (iS .EQ. -1) THEN
            iI = iI+1
            GOTO 5
          END IF
ccc differs from GetMoreInfo becuz of this!      iI = iI + 7 + 8
          iI = iI + 8
          iJ = 130-iI
          iL = iI
          DO iK = 1,iJ
            caTemp(iK:iK) = caStr(iL:iL)
            iL=iL+1
          END DO
          read(caTemp,25) p1,c2,p2
          raX1(iNumP) = p1
          raX2(iNumP) = p2

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
c          print *,iNumP,raX1(iNumP),raX2(iNumP),raP1(iNumP),raP2(iNumP)
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

      iProfileLayers = kProfLayer - iPathOffSet

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
      DO iI = 1 + iPathOffSet,kProfLayer
        IF (raP1(iI) .LT. raP1(iI+1)) THEN
          iTrue = -1
        END IF
      END DO
      IF (iTrue .LT. 0) THEN
        write(kStdErr,*) 'Pressure levels from TEXT klayers file are not in'
        write(kStdErr,*) 'monotonically decreasing or increasing order!!!!!!'
        CALL DoStop
      END IF

      raX1(KProfLayer+1) = raX2(kProfLayer)   !!set highest height
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
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',
     $   raP1(1+iPathOffSet),raP1(2+iPathOffSet)
      write (kStdWarn,*) 'Highest layer : press levels (mb) = ',
     $   raP1(kProfLayer),raP1(kProfLayer+1)
      write (kStdWarn,*) '2 Lowest layers thickness (km) = ',
     $   raH1(1+iPathOffSet),raH1(2+iPathOffSet)
      write (kStdWarn,*) '2 Highest layers thickness (km) = ',
     $  raH1(kProfLayer-1),raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
c see ../INCLUDE/KCARTA_database.param for the kCARTA database definitions
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',
     $     raP1(1+iPathOffSet)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $     raP1(kProfLayer+1)

c this is to dump 
c 101 pressure levels heights for the DEFINITION of the kCARTA database
c this comes from KCARTADATA/RefProf.For.v107up/mypref.op.new
c run bkcarta.x with the profile to be read in being mypref.op.new
       iPrintfile = -1
       IF (iPrintfile .GT. 0) THEN
         print *,'c from KCARTADATA/RefProf.For.v107up/mypref.op.new'
         print *,'c the heights of the pressure levels for the kCARTA DATABASE'
         print *,'      REAL DATABASELEVHEIGHTS(101)'
         print *,' '
         print *,'C-----------------------------------------------------------'
         print *,'PLEV with the AIRS layer boundary pressure heights (in km)'
         print *,'C-----------------------------------------------------------'
         print *,'      DATA  (DATABASELEVHEIGHTS(IPLEV), IPLEV = 101, 1, -1 )'
         DO iI = 1,20
           iK = kProfLayer - (iI-1)*5 + 1
           iL = iK - 5 + 1
           r1 = raX1(iK)
           r2 = raX1(iK-1)
           r3 = raX1(iK-2)
           r4 = raX1(iK-3)
           r5 = raX1(iK-4)
           IF (iI .EQ. 1) THEN
             write (*,5001) r1,r2,r3,r4,r5
           ELSE
             write (*,5002) r1,r2,r3,r4,r5
           END IF
         END DO
        iJ = 1
        write (*,5003) raX1(iJ)
      END IF

 5001 FORMAT ('     $   / ',5(F11.5 ,', '))
 5002 FORMAT ('     $     ',5(F11.5 ,', '))
 5003 FORMAT ('     $     ',F11.5, '/')

      RETURN
      END

c************************************************************************
c this subroutine reads in USER SUPPLIED LEVELS profile (in a text file) and
c outputs the results in terms of what kCARTA wants
c
c raaAmt in moles/cm2, raaPress,raaPartPress in atm, raaTemp in K, raPlevs in mb
c raaAmt in moles/cm2, raaPress,raaPartPress in atm, raaTemp in K, raPlevs in mb
c raaAmt in moles/cm2, raaPress,raaPartPress in atm, raaTemp in K, raPlevs in mb
c
      SUBROUTINE UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress,
     $                 raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $                 iNPath,caPfName,iRTP,
     $                 iProfileLayers,raPressLevels,raTPressLevels,raThickness)

      IMPLICIT NONE

      INTEGER iplev
      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/KCARTA_database.param'

c input
      CHARACTER*80 caPfName
      INTEGER iNpath,iRTP
      INTEGER iaGases(kMaxGas)
c output
      REAL raaTemp(kProfLayer,kGasStore)        !! in K
      REAL raaPress(kProfLayer,kGasStore)       !! in atm
      REAL raaAmt(kProfLayer,kGasStore)         !! in moles/m2 --> need to go to molecules/cm2
      REAL raaPartPress(kProfLayer,kGasStore)   !! in atm
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)  !! plevs in mb, temps in K
      REAL raThickness(kProfLayer) 
      INTEGER iProfileLayers,iKnowTP,iAFGLProf
      REAL raLayerHeight(kProfLayer)
      INTEGER iaWhichGasRead(kMaxGas),iNumGases,iZbndFinal

c local var
      REAL raPbndFinal(kProfLayer+1),raTbndFinal(kProfLayer+1)
      INTEGER iLowestLevLVL2LAY,iNumGasesLVL2LAY,iaGasesLVL2LAY(kMaxGas),iaMap(kMaxGas),iaGasIDMap(kMaxGas)
      INTEGER iL,iG,iNumberofGasesRead,iaInputOrder(kMaxGas),iCnt,iOffSet
      REAL raPoutLVL2LAY(kProfLayer)                  !! in N/m2 --> so need to go to mb
      REAL raToutLVL2LAY(kProfLayer)  
      REAL raZoutLVL2LAY(kProfLayer+1)                !! in meters, notice this is at LEVELS
      REAL raAmountOutLVL2LAY(kProfLayer)             !! in molecules/cm2
      REAL raaQoutLVL2LAY(kProfLayer,kMaxGas)         !! in moles/m2 --> need to go to molecules/cm2
      REAL raaPartPressoutLVL2LAY(kProfLayer,kMaxGas) !! in N/m2 --> so need to go to mb
      REAL raaHeight(kProfLayer,kMaxGas)
      REAL rT,rAmt
      INTEGER iaBnd(kProfLayer+1,2)
      REAL    raBndFrac(kProfLayer+1,2)      
      
      iCnt = 0
      DO iL = 1,kMaxGas
        iaMap(iL)      = -1
        iaGasIDMap(iL) = -1
        IF (iaGases(iL) .EQ. +1) THEN
          iCnt = iCnt + 1
          iaMap(iCnt)    = iaGases(iL)
          iaGasIDMap(iL) = iCnt
        END IF
      END DO

      IF ((kRTP .EQ. -10) .OR. (kRTP .EQ. -5) .OR. (kRTP .EQ. -6)) THEN
        !! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
	!! even though pressures in TAPE5, TAPE6, text file are in mb
	!!    <<<< pressures raPoutLVL2LAY,raaPartPressoutLVL2LAY are in N/m2 >>>> <<< raPBndFInal is in mb >>>>
	!!    <<<< pressures raPoutLVL2LAY,raaPartPressoutLVL2LAY are in N/m2 >>>> <<< raPBndFInal is in mb >>>>	
	!!    <<<< pressures raPoutLVL2LAY,raaPartPressoutLVL2LAY are in N/m2 >>>> <<< raPBndFInal is in mb >>>>
	write(kStdErr,*)  '>>> RUNNING INTERNAL KLAYERS to integrate text levels/layers --> layers'
	write(kStdErr,*)  '>>>    this is quick and convenient; highly recommend you use klayers instead!'
	write(kStdWarn,*) '>>> RUNNING INTERNAL KLAYERS to integrate text levels/layers --> layers'
	write(kStdWarn,*) '>>>    this is quick and convenient; highly recommend you use klayers instead!'
	
        CALL InputMR_profile(caPfName,iProfileLayers,iNumGasesLVL2LAY,iaGasesLVL2LAY,iLowestLevLVL2LAY,
     $                     raToutLVL2LAY,raAmountOutLVL2LAY,raZoutLVL2LAY,
     $                     raPoutLVL2LAY,raaQoutLVL2LAY,raaPartPressoutLVL2LAY,
     $                     raPbndFinal,raTBndFinal,iZbndFinal)
        write(kStdWarn,*) 'Read ',iNumGasesLVL2LAY,' gases at ',iProfileLayers,' output layers from the text file'
	IF (iZbndFinal .GT. 0) THEN
	  write(kStdWarn,*) 'hhmmm looks like we have new pressure levels from LBLRTM TAPE5 and/or TAPE6'
	END IF

        write(kStdErr,*)  '>>> FINISHED INTERNAL KLAYERS to integrate text levels/layers --> layers'
        write(kStdWarn,*) '>>> FINISHED INTERNAL KLAYERS to integrate text levels/layers --> layers'	

      ELSE
        write(kStdErr,*) 'huh?? only use this routine if reading text levels file or text LBLRTM TAPE5,6???',kRTP
        Call DoStop
      END IF

      iNumberofGasesRead = 0
      DO iG = 1,kMaxGas
        iaWhichGasRead(iG) = -1
      END DO
      DO iG = 1,iNumGasesLVL2LAY
        iNumberofGasesRead = iNumberofGasesRead + 1
        iaWhichGasRead(iaGasesLVL2LAY(iG)) = +1
      END DO

      !!    <<<< pressures raPoutLVL2LAY,raaPartPressoutLVL2LAY are in N/m2 >>>> <<< raPBndFInal is in mb >>>>
      iOffSet = kProfLayer - (iZbndFinal-1)
      DO iL = 1,iZbndFinal
        !! routine only called in kRPT = -10,-6,-5 so use New Plevs      
        !! raPressLevels(iL) = PLEV_KCARTADATABASE_AIRS(iL)
        raPressLevels(iL+iOffSet) = raPbndFinal(iL)
	!!also set the temps
        raTPressLevels(iL+iOffSet) = raTBndFinal(iL)
c	print *,iL,iOffset,raPressLevels(iL+iOffSet),raTPressLevels(iL+iOffSet)
      END DO

 345  FORMAT(I3,2(' ',F10.3),2(' ',I3,F10.3,' ',F10.3))
      DO iL = 1,iZbndFinal-1
        !! find plev_airs which is just ABOVE the top of current layer
        iG = kProfLayer+1
 10     CONTINUE
	IF ((PLEV_KCARTADATABASE_AIRS(iG) .LE. raPbndFinal(iL+1)) .AND. (iG .GT. 1)) THEN
	  iG = iG - 1
	  GOTO 10
	ELSE
	  iaBnd(iL,2) = min(iG+1,kMaxLayer+1)   !! top bndry of plevs_database is lower pressure than top bndry of raPbndFinal layer iL
	END IF
	raBndFrac(iL,2) = (raPbndFinal(iL+1)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))/
     $                  (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))

        !! find plev_airs which is just BELOW the bottom of current layer
        iG = 1
 20     CONTINUE	
	IF (PLEV_KCARTADATABASE_AIRS(iG) .GT. raPbndFinal(iL)) THEN
	  iG = iG + 1
	  GOTO 20
	ELSE
	  iaBnd(iL,1) = max(iG-1,1) !! bot boundary of plevs_database is bigger pressure than top bndry of raPbndFinal layer iL
	END IF
	raBndFrac(iL,1) = (raPbndFinal(iL)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))/
     $                  (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))	
	
c      write (*,345) iL,raPbndFinal(iL),raPbndFinal(iL+1),iaBnd(iL,1),raBndFrac(iL,1),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)),
c     $                                                   iaBnd(iL,2),raBndFrac(iL,2),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))
      END DO

 999  FORMAT(I3,' ',A20,4(' ',ES10.3))
      DO iL = 1,iOffSet
        raLayerHeight(iL) = 0.0
        raThickness(iL)   = 0.0
        DO iG = 1,iNumGasesLVL2LAY
          raaPress(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))     = 0.0
          raaTemp(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))      = 0.0
          raaAmt(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))       = 0.0
          raaPartPress(iL,iaGasIDMap(iaGasesLVL2LAY(iG))) = 0.0
          raaHeight(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))    = 0.0
        END DO
      END DO
      DO iL = iOffSet+1,kProflayer
        raLayerHeight(iL) = 0.5*(raZoutLVL2LAY(iL) + raZoutLVL2LAY(iL+1))
        raThickness(iL)   = raZoutLVL2LAY(iL+1) - raZoutLVL2LAY(iL)
d	write(*,999) iL,' n pth mix',raLayerHeight(iL),raThickness(iL),raPoutLVL2LAY(iL),raPoutLVL2LAY(iL)/1013.255
        DO iG = 1,iNumGasesLVL2LAY
          raaPress(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))     = raPoutLVL2LAY(iL)/100/1013.255  !! N/m2 --> mb --> atm
          raaTemp(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))      = raToutLVL2LAY(iL)
          raaAmt(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))       = raaQoutLVL2LAY(iL,iG)/kAvog !! change molecules/cm2 to kilomoles/cm2
          raaPartPress(iL,iaGasIDMap(iaGasesLVL2LAY(iG))) = raaPartPressoutLVL2LAY(iL,iG)/100/1013.255    !! N/m2 --> mb --> atm
          raaHeight(iL,iaGasIDMap(iaGasesLVL2LAY(iG)))    = raLayerHeight(iL)
        END DO
      END DO
      
      DO iG = 1,kMaxGas
        iaInputOrder(iG) = -1
      END DO
      iL = 1
      DO iG = 1,kMaxGas
        IF (iaGases(iG) .GT. 0) THEN
          iaInputOrder(iL) = iG
          iL = iL+1
        END IF
      END DO

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      !!! no need to worry about finky plevs, as all the hard work is already done
      !!! in reading in the WaterLayers  Profile from RTP or LBLRTM
      CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,
     $          iaInputOrder,iNumGases,
     $          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

c now got to add in the missing gases, see READRTP_1B
c this was the original code
c      CALL AddOnAFGLProfile(kAFGLProf,
c     $         iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
c     $         raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
      CALL AddOnAFGLProfile_arblevels(kAFGLProf,
     $         iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $         raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness,
     $         iaBnd,raBndFrac,raPbndFinal,iZbndFinal)

      RETURN
      END

c************************************************************************
c this subroutine adds on RefProfile 1,2,3,4,5,6 gas amounts for those gases NOT
c in the RTP profile (using h.ptype = 2 or even 1)
c the afgl profiles generated by /home/sergio/KCARTA/UTILITY/AFGLprofs.m
c see subr MakeRefProf
      SUBROUTINE AddOnAFGLProfile_arblevels(iProfileNum,
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness,
     $      iaBnd,raBndFrac,raPbndFinal,iZbndFinal)

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
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
c these are new, and tell you how to integrate default AITRS 100 layers onto the arb LBLRTM layers
      INTEGER iaBnd(kProfLayer+1,2),iZbndFinal
      REAL    raBndFrac(kProfLayer+1,2)      
      REAL    raPbndFinal(kProfLayer+1)
      
c local variables
      INTEGER iI,iJ,iaNeed(kMaxGas),iNewRead,iFound,iFoundX,iIDgas,iLay,iG0,iK,iX,iY,iOffset,iMaxL
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer)   !!US Std layer ppress, amt at 100 layers
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer)   !!US Std layer press, temp at 100 layers
      REAL raQX2(kMaxProfLayer),raPPX2(kMaxProfLayer) !!need to change to new layers
      REAL rPP,rPPWgt,rFrac,rMR
      CHARACTER*1 cY,cN

      iOffSet = kProfLayer - (iZbndFinal-1)
      
      cY = 'Y'
      cN = 'N'

      IF ((iProfileNum .LT. 1) .OR. (iProfileNum .GT. 6)) THEN
        write(kStdErr,*) 'Can only substitute profs from AFGL '
        write(kStdErr,*) '1=STD, 2=TRP, 3=MLS, 4=MLW, 5=SAS, 6=SAW'
        CALL DoStop
      END IF

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'read profiles for ',iNumberOfGasesRead, ' gases ...'
      write(kStdWarn,*) 'for AFGL profile ',iProfileNum

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
 
      iG0 = 1
      
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
	  
          CALL getAFGL(iProfileNum,iIDgas,raPX,raPPX,raTX,raQX)

c now that we know the weights and boundaries, off we go!!!
c remember pV = nRT ==> p(z) dz/ r T(z) = dn(z)/V = dq(z) ==> Q = sum(p Z / R T)
c so for these fractional combined layers (i), Qnew = sum(p(i) zfrac(i) / R T(i)) = sum(p(i) zfrac(i)/Z(i) Z(i) / RT(i))
c                                                   = sum(p(i)Z(i)/RT(i) zfrac(i)/Z(i))
c or Qnew = sum(frac(i) Q(i))

          !! do the integral from AIRS layers to ARB layers
	  iMaxL = min(kProfLayer,iZbndFinal) !! so if iZbndFinal = 101, we only do 100 layers ....
          DO iX = 1,iMaxL
	    raQX2(iX) = 0.0
	    rPP = 0.0
	    rPPWgt = 0.0
	    rMR = 0.0
	    DO iY = iaBnd(iX,1),iaBnd(iX,2)-1
	      IF (iY .EQ. iaBnd(iX,2)-1) THEN
	        rFrac = raBndFrac(iX,2)
	      ELSEIF (iY .EQ. iaBnd(iX,1)) THEN
	        !! this also takes care of case when iY .EQ. iaBnd(iX,1) .EQ. iaBnd(iX,2)-1
	        rFrac = raBndFrac(iX,1)
	      ELSE
	        rFrac = 1.0
	      END IF
	      rPP = rPP + raPPX(iY)*rFrac
	      rPPWgt = rPPWgt + rFrac
	      rMR = rMR + raPPX(iY)/raPX(iY)	      
	      raQX2(iX) = raQX2(iX) + raQX(iY)*rFrac
	    END DO
	    raPPX2(iX) = rPP/rPPWgt                      !! one way
	    rMR = rMR/((iaBnd(iX,2)-1)-(iaBnd(iX,1))+1)	    
	    raPPX2(iX) = rMR * raPbndFinal(iX)/1013.25   !! another way
c	    write(*,1234) iIDGas,iX,raPbndFinal(iX)/1013.25,raaPress(iX,1),iaBnd(iX,1),iaBnd(iX,2),raQX2(iX),raPPX2(iX),rPP/rPPWgt
	  END DO
 1234     FORMAT(2(' ',I3),2(' ',F10.3),2(' ',I3),3(' ',E10.3))
 
          Call FindIndexPosition(iIDgas,iNumGases,iaInputOrder,
     $                           iFoundX,iGasIndex)
          IF (iFoundX .GT. 0) THEN
            DO iX = 1,kProfLayer
	      raaTemp(iX,iGasIndex) = 00.0	      
	      raaTemp(iX,iGasIndex) = 200.0
	      raaTemp(iX,iGasIndex) = raaTemp(iX,iG0)	      
	    END DO
	  
            DO iLay = 1,iZbndFinal-1
              raaAmt(iLay+iOffSet,iGasIndex)       = raQX2(iLay)
              raaTemp(iLay+iOffSet,iGasIndex)      = raaTemp(iLay+iOffSet,iG0)
              raaPress(iLay+iOffSet,iGasIndex)     = raaPress(iLay+iOffSet,iG0)
              raaPartPress(iLay+iOffSet,iGasIndex) = raPPX2(iLay)
              raaHeight(iLay+iOffSet,iGasIndex)    = raaHeight(iLay+iOffSet,iG0)
c	      IF (iLay  .EQ. 25) print *,iIDGas,iLay,iLay+iOffSet,raaTemp(iLay+iOffSet,iG0)
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
c this subroutine adds on RefProfile 1,2,3,4,5,6 gas amounts for those gases NOT
c in the RTP profile (using h.ptype = 2 or even 1)
c the afgl profiles generated by /home/sergio/KCARTA/UTILITY/AFGLprofs.m
      SUBROUTINE AddOnAFGLProfile(iProfileNum,
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)

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
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      
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
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
        RETURN
      END IF

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'read profiles for ',iNumberOfGasesRead, ' gases ...'
      write(kStdWarn,*) 'for AFGL profile ',iProfileNum

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
      CHARACTER*120 caFname0,caFname,cnameX
      CHARACTER c1,c2
      CHARACTER*2 c12
      CHARACTER*100 caLine

      DO iLen = 1,120
        caFname0(iLen:iLen) = ' '
      END DO
      
      caFname0(1:80) = kOrigRefPath

      iLen = 120
 100  CONTINUE
      IF (caFname0(iLen:iLen) .EQ. ' ') THEN
        iLen = iLen - 1
        GOTO 100
      END IF

      IF (kAFGLProf .NE. 1) THEN
        caFName0(iLen:iLen) = '/'
        iLen = iLen+1

        IF (kAFGLProf .EQ. 2) caFName0(iLen:iLen) = '2'
        IF (kAFGLProf .EQ. 3) caFName0(iLen:iLen) = '3'
        IF (kAFGLProf .EQ. 4) caFName0(iLen:iLen) = '4'
        IF (kAFGLProf .EQ. 5) caFName0(iLen:iLen) = '5'      
        IF (kAFGLProf .EQ. 6) caFName0(iLen:iLen) = '6'
        iLen = iLen+1

        cnameX = 'afglgas'
      ELSE
        ! IF (kAFGLProf .EQ. 1) caFName0(iLen:iLen) = '1'   !! no need to go down a dir for US STD    
        cnameX = 'us_std_gas_'
        cnameX = 'refgas'
      END IF
      
      caFName0(iLen:iLen) = '/'

      iLenX = 120
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
      DO i1 = 1,120
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
c this subroutine adds on US STandard Profile gas amounts for those gases NOT
c in the RTP profile (using h.ptype = 2 or even 1)

c note : have to be careful about 101 AIRS levls versus actual raPressLevels

      SUBROUTINE AddOnStandardProfile(
     $      iNumberofGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,
     $      raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)

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
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      
c local variables
      INTEGER iI,iJ,iaNeed(kMaxGas),iNewRead,iFound,iFoundX,iIDgas,iLay,iG0,iK
      REAL raPPX(kMaxProfLayer),raQX(kMaxProfLayer)  !!US Std layer ppress, amt
      REAL raPX(kMaxProfLayer), raTX(kMaxProfLayer)  !!US Std layer press, temp
      CHARACTER*1 cY,cN

      cY = 'Y'
      cN = 'N'

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'read in USER SUPPLIED profiles for ',iNumberOfGasesRead, ' gases ...'

c this is the list of gases for which we need profiles
      iK = 0
      write(kStdWarn,*) '            count gasID  found?'
      write(kStdWarn,*) '-------------------------------'
      DO iI = 1,kMaxGas
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

c ---->>> this already taken care of in previous loop  July 2014
c      !!! now do gasids 101,102,103
c      iI = 1
c      IF (iaInputOrder(iI) .EQ. 1) THEN
c        IF (iaWhichGasRead(iI) .EQ. +1) THEN
c          write(kStdWarn,200) iI,101,cY
c          write(kStdWarn,200) iI,102,cY
c          write(kStdWarn,200) iI,103,cY
c        ELSE  
c          write(kStdWarn,200) iI,101,cN
c          write(kStdWarn,200) iI,102,cN
c          write(kStdWarn,200) iI,103,cN
c        END IF
c      END IF
c ---->>> this already taken care of in previous loop  July 2014

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

c      DO iI = 1,kProfLayer+1
c        print *,iI,raPressLevels(iI)
c      END DO
c      call dostopmesg('allison$')
      
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
      INTEGER iDefault,iAddLBLRTM

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

      iDefault = -1
      iAddLBLRTM = +1   !! if profile missing and kRTP = -5 or -6, do     add
      iAddLBLRTM = -1   !! if profile missing and kRTP = -5 or -6, do not add
      iAddLBLRTM = iaaOverrideDefault(3,3)
      IF (abs(iAddLBLRTM) .NE. 1) THEN
        write(kStdErr,*) 'invalid iAddLBLRTM = ',iAddLBLRTM
        CALL DoStop
      END IF		       
      
      IF ((kRTP .NE. -5) .AND. (kRTP .NE. -6)) THEN
c       write(kStdWarn,*) 'need to add on US Std for ',iIDgas, ' from ',caFname
        write(kStdWarn,*) 'need to add on US Std for ',iIDgas
      ELSEIF (((kRTP .EQ. -5) .OR. (kRTP .EQ. -6)) .AND. (iAddLBLRTM .EQ. -1)) THEN
        write(kStdWarn,*) 'need to add on US Std for ',iIDgas,' LBLRTM set raQ=0'
      ELSEIF (((kRTP .EQ. -5) .OR. (kRTP .EQ. -6)) .AND. (iAddLBLRTM .EQ. +1)) THEN
        write(kStdWarn,*) 'need to add on US Std for ',iIDgas,' LBLRTM set raQ=0, reset and add profile'
        write(kStdErr,*) 'need to add on US Std for ',iIDgas,' LBLRTM set raQ=0, reset and add profile'
      END IF

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
c      iLay = Ilay + 1 
C      READ(iIOUN,*,END=1030) iX,raPx(iLay),raPPx(iLay),raTx(iLay),raQx(iLay) 
C      Goto 1020 
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
        IF (((kRTP .EQ. -5) .OR. (kRTP .EQ. -6)) .AND. (iAddLBLRTM .EQ. -1)) THEN
          raPPx(iLay) = 0.0
          raQx(iLay) = 0.0
        END IF
      ENDIF
      GOTO 20
 199  CLOSE(iIOUN)
      kTempUnitOpen = -1

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
