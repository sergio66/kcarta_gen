c Copyright 2004  
c University of Maryland Baltimore County  
c All Rights Reserved 

c************************************************************************
c this subroutine interpolates the SolarAngle for iRegr, and dumps output
c into a vt file for kCARTA to read
      SUBROUTINE OutputVTFile(caVTFile,iRTP,iRegr,rSolarAngle)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      REAL rSolarAngle            !!!solar angle for atm #1 
      INTEGER iRegr               !!!which regr profile this is closest to
      INTEGER iRTP                !!!which of the rtp files processed
c output parameters
      CHARACTER*80 caVTFile       !!!temp file name where VT profiles are
                                  !!!created and stored if not given by user

c local vars
c these are used to read in info from the 6 angle profs 0,40,60,80,85,90
c MXXTMP = 6
c kMaxUserSet = 12
c kMaxK = 50
      REAL raaPressure(MXXTMP,kNLTEProfLayer)
      REAL raaKinetic(MXXTMP,kNLTEProfLayer) 
      REAL raaaVibPartFcn(MXXTMP,kMaxUserSet,kNLTEProfLayer)
      REAL raaaVibTemp(MXXTMP,kMaxK,kNLTEProfLayer)
      CHARACTER*80 caFName
      CHARACTER*80 caaComments(kMaxLayer)
      INTEGER iNumComments
c to dump out the profiles
      REAL raaVibPartFcn(kMaxUserSet,kNLTEProfLayer)
      REAL raaVibTemp(kMaxK,kNLTEProfLayer)

c other local vars
      REAL raSolAngles(MXXTMP)
      INTEGER iaSolAngles(MXXTMP)
      INTEGER iI,iJ,iK,iAngle,iFound
      INTEGER iGasesInFile,iNumVibLevels            !gases,levels in file
      INTEGER iGasID,iVibPF                         !gas ID, number of isotopes
      INTEGER iNUMVibTempMagic                      !number of vib profiles
      INTEGER iGasesInFile1,iNumVibLevels1          !gases,levels in file
      INTEGER iGasID1,iVibPF1                       !gas ID, number of isotopes
      INTEGER iNUMVibTempMagic1                     !number of vib profiles
      INTEGER iNumSolAngles
      INTEGER iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
      REAL raVibCntr(kMaxK)
      INTEGER iaCnt1(kMaxK),iaGasID1(kMaxK),iaISO1(kMaxK),iaAFGL1(kMaxK)
      REAL raVibCntr1(kMaxK)
      REAL ra6(MXXTMP),rT

      DATA iNumSolAngles/6/
      DATA (raSolAngles(iI),iI=1,6) /0.0,40.0,60.0,80.0,85.0,90.0/
      DATA (iaSolAngles(iI),iI=1,6) /0  ,40  ,60  ,80  ,85  ,90  /

      IF (MXXTMP .LT. 6) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume MXXTMP = 6'
        CALL DoStop
      END IF
      IF (kMaxUserSet .LT. 12) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume kMaxUserSet = 12'
        CALL DoStop
      END IF
      IF (kMaxK .LT. 50) THEN
        write(kStdErr,*) 'Sorry .. NLTE interpolations assume kMaxK = 50'
        CALL DoStop
      END IF

      write(kStdWarn,*) 'kCARTA : estimate VibTemps based on following info : '
      write(kStdWarn,*) 'Input RTP profile               = ',iRTP
      write(kStdWarn,*) 'Closest Regression profile      = ',iRegr
      write(kStdWarn,*) 'Solar Angle                     = ',rSolarAngle
      write(kStdWarn,*) 'Output vib temperature filename = ',caVTFile

c start reading in the info!
      DO iAngle = 1,iNumSolAngles
        caFName = kLopezPuertas
        CALL addinfo(caFname,iRegr,iaSolAngles(iAngle)) 
        CALL VibTempProfRead(iAngle,caFName,caaComments,iNumComments,
     $    raaPressure,raaKinetic,raaaVibPartFcn,raaaVibTemp,
     $    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr,
     $    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)
        IF (iAngle .EQ. 1) THEN
          iNUMVibTempMagic1  = iNUMVibTempMagic
          iGasesInFile1      = iGasesInFile
          iNumVibLevels1     = iNumVibLevels
          iGasID1            = iGasID
          iVibPF1            = iVibPF
          DO iJ = 1,iNUMVibTempMagic
            iaCnt1(iJ)     = iaCnt(iJ)
            iaGasID1(iJ)   = iaGasID(iJ)
            iaISO1(iJ)     = iaISO(iJ)
            iaAFGL1(iJ)    = iaAFGL(iJ)
            raVibCntr1(iJ) = raVibCntr(iJ)
          END DO
        ELSE
          !check to ensure parameters are ok
          IF (iGasesInFile1 .NE. iGasesInFile) THEN
            write(kStdErr,*) 'Error in iGasesInfile!!!'
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          ELSEIF (iNumVibLevels1 .NE. iNumVibLevels) THEN
            write(kStdErr,*) 'Error in iNumVibLevels!!!'
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          ELSEIF (iGasID1 .NE. iGasID) THEN
            write(kStdErr,*) 'Error in iGasID!!!'
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          ELSEIF (iVibPF1 .NE. iVibPF) THEN
            write(kStdErr,*) 'Error in iVibPF!!!'
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          ELSEIF (iNUMVibTempMagic1 .NE. iNUMVibTempMagic) THEN
            write(kStdErr,*) 'Error in iNUMVibTempMagic!!!'
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          END IF

          !!!if those tests ok, make sure HITRAN IDS ok
          iFound = -1
          DO iJ = 1,iNUMVibTempMagic
            IF (iaCnt1(iJ)     .NE. iaCnt(iJ))     iFound = +1
            IF (iaGasID1(iJ)   .NE. iaGasID(iJ))   iFound = +2
            IF (iaISO1(iJ)     .NE. iaISO(iJ))     iFound = +3
            IF (iaAFGL1(iJ)    .NE. iaAFGL(iJ))    iFound = +4
            IF (raVibCntr1(iJ) .NE. raVibCntr(iJ)) iFound = +5
          END DO
          IF (iFound .GT. 0) THEN
            write(kStdErr,*) 'Error checking out HITRAN GasIDs ',iFound
            write(kStdErr,*) 'in subroutine OutputVTFile'
            CALL DoStop
          END IF
     
        END IF
      END DO

c interpolate, and then dump out the info!
c first check to see if user SolarAngle == one of the 6 angles we have
      iFound = -1
      DO iI = 1,iNumSolAngles
        IF (abs(rSolarAngle - raSolAngles(iI)) .LE. 0.1) THEN
          iFound = iI
          GOTO 10
        END IF
      END DO

 10   CONTINUE
      IF (iFound .GT. 0) THEN 
        !!no need to interpolate; simply set raaVibPartFcn,raaVibTemp
        DO iJ = 1,iNumVibLevels
          DO iI = 1,iVibPF
            raaVibPartFcn(iI,iJ) =  raaaVibPartFcn(iFound,iI,iJ) 
          END DO
        END DO

        DO iJ = 1,iNumVibLevels
          DO iI = 1,iNUMVibTempMagic
            raaVibTemp(iI,iJ) =  raaaVibTemp(iFound,iI,iJ) 
          END DO
        END DO

      ELSEIF (iFound .LT. 1) THEN 
        !!need to interpolate every layer, every variable, in sunangle
        DO iI = 1,iNumVibLevels !!!scan over the 120 pressure levels
          DO iJ = 1,iVibPf      !!!scan over the different isotopes
            DO iK = 1,MXXTMP    !!!scan over the 6 solar angles
              ra6(iK) = raaaVibPartFcn(iK,iJ,iI)
            END DO
            CALL rspl(raSolAngles,ra6,6,rSolarAngle,rT,1)
            raaVibPartFcn(iJ,iI) =  rT
          END DO
        END DO

        DO iI = 1,iNumVibLevels       !!!scan over the 120 pressure levels
          DO iJ = 1,iNUMVibTempMagic  !!!scan over the different vib levels
            DO iK = 1,MXXTMP          !!!scan over the 6 solar angles
              ra6(iK) = raaaVibTemp(iK,iJ,iI)
            END DO
            CALL rspl(raSolAngles,ra6,6,rSolarAngle,rT,1)
            raaVibTemp(iJ,iI) =  rT
          END DO
        END DO
      END IF

      !!!dump out the info into caVTFile
      CALL VibTempProfWrite(iAngle,caVTFile,caaComments,iNumComments,
     $    raaPressure,raaKinetic,raaVibPartFcn,raaVibTemp,
     $    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr,
     $    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)

      RETURN
      END

c************************************************************************
c this subroutine writes out the vib temps
c in the format that GENLN2 and kCARTA like 
      SUBROUTINE VibTempProfWrite(iAngle,caFName,caaComments,iNumComments,
     $    raaPressure,raaKinetic,raaVibPartFcn,raaVibTemp,
     $    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr,
     $    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)

      IMPLICIT NONE 

      include '../INCLUDE/kcarta.param'
 
c input vars 
      CHARACTER*80 caFname                          !i/o file name 
      INTEGER iAngle                                !which angle profile
      REAL raaPressure(MXXTMP,kNLTEProfLayer)
      REAL raaKinetic(MXXTMP,kNLTEProfLayer) 
      REAL raaVibPartFcn(kMaxUserSet,kNLTEProfLayer)
      REAL raaVibTemp(kMaxK,kNLTEProfLayer)
      CHARACTER*80 caaComments(kMaxLayer)
      INTEGER iNumComments
      INTEGER iGasesInFile,iNumVibLevels            !gases,levels in file
      INTEGER iGasID,iVibPF                         !gas ID, number of isotopes
      INTEGER iNUMVibTempMagic                      !number of vib profiles
      INTEGER iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
      REAL raVibCntr(kMaxK)

c local vars 
      CHARACTER*3 ca3
      CHARACTER*1 ca1
      CHARACTER*80 caStr
      INTEGER i1,i2,iIOUN,iErr,iI,iJ,iDummy1,iDummyGasID,iDummyISO
      INTEGER iNV,iZero
      REAL raP2(kNLTEProfLayer),raT2(kNLTEProfLayer)

      iZero = 0

      iNV = iNumVibLevels
      DO iI = 1,iNV
        raP2(iI) = raaPressure(1,iI)
        raT2(iI) = raaKinetic(1,iI)
      END DO

      iIOUN = kTempUnit
      OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='UNKNOWN',  
     $     IOSTAT=iErr)   
      IF (iErr .NE. 0) THEN   
        write(kStdErr,*) 'error writing VibTemp profile ',iErr,caFName
        CALL DoStop
      END IF   
 
      kTempUnitOpen = +1

      caaComments(1) = '! INTERPD (sun angle) vib temp profiles reformatted'
      DO i1 = 1,iNumComments
        caStr = caaComments(i1)
        write(iIOUN,123) caStr 
      END DO

        !!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV 
        !!!where FLAG  = *** 
        !!!      MODEL = atmosphere model number for data file 
        !!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL 
        !!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB      
        !!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL  
        !!! so read in flag, model,nset,ngas,nlev where flag = '***' 
      caStr = '***' 
      write(iIOUN,99) caStr,1,1,iGasesInFile,iNumVibLevels 
      write(iIOUN,*) iGasID,iVibPF ! how many vib part fcns for the NLTE gas 
      
      CALL printarray(iIOUN,iNV,raP2)  !!!print out press levels 

      write(iIOUN,*) iZero       !dummy 
      CALL printarray(iIOUN,iNV,raT2)  !!!print out kinetic temps 
 
      caStr = '!QV   Vibrational Partition functions ' 
      write(iIOUN,13) caStr 
      DO iI = 1,iVibPf 
        DO iJ = 1,iNV
          raP2(iJ) = raaVibPartFcn(iI,iJ)
        END DO
        write(iIOUN,*) iGasID,iI 
        CALL printarray(iIOUN,iNV,raP2)  !!!print out part fcn 
      END DO 

      caStr = '!TV   Vibrational Temperatures' 
      write(iIOUN,13) caStr 
      DO iI = 1,iNUMVibTempMagic
        DO iJ = 1,iNV
          raP2(iJ) = raaVibTemp(iI,iJ)
        END DO
        write(iIOUN,20) iI,iaGasID(iI),iaISO(iI),iaAFGL(iI),raVibCntr(iI)
        CALL printarray(iIOUN,iNV,raP2) 
      END DO 

      CLOSE(iIOUN)
      kTempUnitOpen = -1

 10   FORMAT(A80) 
 11   FORMAT('! ',A70) 
 12   FORMAT(A70,I5) 
 13   FORMAT(A70) 
 20   FORMAT(I4,'    ',I5,'  ',I7,' ',I6,'   ',F11.5) 
 99   FORMAT(A3,'   ',4(I3,' ')) 

 123  FORMAT(A80)

      RETURN
      END

c************************************************************************
c this subroutine reads in the vib temps
c in the format that GENLN2 and kCARTA like 
      SUBROUTINE VibTempProfRead(iAngle,caFName,caaComments,iNumComments,
     $    raaPressure,raaKinetic,raaaVibPartFcn,raaaVibTemp,
     $    iaCnt,iaGasID,iaISO,iaAFGL,raVibCntr,
     $    iGasesInFile,iNumVibLevels,iGasID,iVibPF,iNUMVibTempMagic)
 
      IMPLICIT NONE 

      include '../INCLUDE/kcarta.param'
 
c input vars 
      CHARACTER*80 caFname                          !i/o file name 
      INTEGER iAngle                                !which angle profile
c output vars
      REAL raaPressure(MXXTMP,kNLTEProfLayer)
      REAL raaKinetic(MXXTMP,kNLTEProfLayer) 
      REAL raaaVibPartFcn(MXXTMP,kMaxUserSet,kNLTEProfLayer)
      REAL raaaVibTemp(MXXTMP,kMaxK,kNLTEProfLayer)
      CHARACTER*80 caaComments(kMaxLayer)
      INTEGER iNumComments
      INTEGER iGasesInFile,iNumVibLevels            !gases,levels in file
      INTEGER iGasID,iVibPF                         !gas ID, number of isotopes
      INTEGER iNUMVibTempMagic                      !number of vib profiles
      !!!next few are the HITRAN id's
      INTEGER iaCnt(kMaxK),iaGasID(kMaxK),iaISO(kMaxK),iaAFGL(kMaxK)
      REAL raVibCntr(kMaxK)
c local vars 
      CHARACTER*3 ca3
      CHARACTER*1 ca1
      CHARACTER*80 caStr
      INTEGER i1,i2,iIOUN,iErr,iI,iJ,iDummy1,iDummyGasID,iDummyISO

      iIOUN = kTempUnit
      OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',  
     $     IOSTAT=iErr)   
      IF (iErr .NE. 0) THEN   
        write(kStdErr,*) 'error reading VibTemp profile ',iErr,caFName
        CALL DoStop
      END IF   
 
      kTempUnitOpen = +1

      iNumComments = 0
 555  CONTINUE 
      read(iIOUN,123) caStr 
      IF (caStr(1:1) .EQ. '!') THEN 
        iNumComments = iNumComments + 1
        caaComments(iNumComments) = caStr
        !!!these are comments at beginning of file, so skip 
        GOTO 555 
      ELSEIF (caStr(1:3) .EQ. '***') THEN 
        !!!aha .. now we are cooking 
        !!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV 
        !!!where FLAG  = *** 
        !!!      MODEL = atmosphere model number for data file 
        !!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL 
        !!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB      
        !!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL  
        !!! so read in flag, model,nset,ngas,nlev where flag = '***' 
        read(caStr,456) ca3,i1,i2,iGasesInFile,iNumVibLevels 
      END IF 

      caStr = caaComments(iNumComments-1)
      read(caStr,789) ca1,iNUMVibTempMagic

      IF (iGasesInFile .NE. 1) THEN
        write(kStdErr,*) 'kCARTA can only handle one NLTE gas (CO2) here!!!!'
        CALL DoStop
      END IF
      IF (iNumVibLevels .GT. kNLTEProfLayer) THEN
        write(kStdErr,*) 'kCARTA can handle max',kNLTEProfLayer,' NLTE levels!'
        CALL DoStop
      END IF

      READ(iIOUN,*) iGasID,iVibPF
      IF (iGasID .NE. 2) THEN
        write(kStdErr,*) 'kCARTA can only handle CO2 for NLTE here!!!!'
        CALL DoStop
      END IF

      IF (iVibPF .GT. kMaxK) THEN
        write(kStdErr,*) 'kCARTA can only handle max ',kMaxK,'Vib part fcns!!'
        CALL DoStop
      END IF

      !!!read the presssure levels 
      read(iIOUN,*) (raaPressure(iAngle,iI),iI=1,iNumVibLevels)   !!! in mb 
 
      !!! read the kinetic temperatures 
      read(iIOUN,*) iDummy1 
      read(iIOUN,*) (raaKinetic(iAngle,iI),iI=1,iNumVibLevels)   !!!LTE temp
 
      !!!read in the QTIPS_VIB modifiers 
      !!!dummy string that says "!QV   Vibrational Partition functions" 
      read(iIOUN,123) caStr       
      DO iJ = 1,iVibPF
        read(iIOUN,*) iDummyGasID,iDummyISO 
        read(iIOUN,*) (raaaVibPartFcn(iAngle,iJ,iI),iI=1,iNumVibLevels)
      END DO

      !!! read the NLTE temperatures 
      !!!dummy string that says "!TV   Vibrational Temperatures" 
      read(iIOUN,123) caStr       
      DO iJ = 1,iNUMVibTempMagic
        read(iIOUN,*) iaCnt(iJ),iaGasID(iJ),iaISO(iJ),iaAFGL(iJ),raVibCntr(iJ)
        read(iIOUN,*) (raaaVibTemp(iAngle,iJ,iI),iI=1,iNumVibLevels)
      END DO

      CLOSE(iIOUN)
      kTempUnitOpen = -1

 123  FORMAT(A80)
 456  FORMAT(A3,'   ',4(I3,' ')) 
 789  FORMAT(A1,'   ',I3)

      RETURN 
      END 

c************************************************************************ 
c this subroutine prints out the arrays in a nice format 
      SUBROUTINE printarray(iIOUN,iNV,raX) 

      IMPLICIT NONE  
 
      include '../INCLUDE/kcarta.param' 
   
      INTEGER iIOUN,iNV,iimod,iDiv 
      REAL raX(kNLTEProfLayer) 
 
c local vars  
      INTEGER iM,iD,i1,i2,ia(5),iJ,iFive 
      REAL ra(5) 
 
      iFive = 5 
      iM = iimod(iNv,iFive) 
      iD = iDiv(iNv,iFive) 
 
      ia(1) = -4 
      ia(2) = -3 
      ia(3) = -2 
      ia(4) = -1 
      ia(5) =  0 
      DO i1 = 1,iD 
        DO i2 = 1,5 
          ia(i2) = ia(i2) + 5 
          ra(i2) = raX(ia(i2)) 
        END DO 
        write(iIOUN,15) (ra(iJ),iJ=1,5) 
      END DO 
 
      IF ((iM .GT. 0) .AND. (iM .LT. 5)) THEN 
        IF (iM .EQ. 1) THEN 
          ra(1) = raX(iNV) 
          write(iIOUN,11) (ra(iJ),iJ=1,iM) 
        ELSE IF (iM .EQ. 2) THEN 
          ra(2) = raX(iNV) 
          ra(1) = raX(iNV-1) 
          write(iIOUN,12) (ra(iJ),iJ=1,iM) 
        ELSE IF (iM .EQ. 3) THEN 
          ra(3) = raX(iNV) 
          ra(2) = raX(iNV-1) 
          ra(1) = raX(iNV-2) 
          write(iIOUN,13) (ra(iJ),iJ=1,iM) 
        ELSE IF (iM .EQ. 4) THEN 
          ra(4) = raX(iNV) 
          ra(3) = raX(iNV-1) 
          ra(2) = raX(iNV-2) 
          ra(1) = raX(iNV-3) 
          write(iIOUN,14) (ra(iJ),iJ=1,iM) 
        END IF 
      END IF 
 
 11   FORMAT(5(E11.5,' ')) 
 12   FORMAT(5(E11.5,' ')) 
 13   FORMAT(5(E11.5,' ')) 
 14   FORMAT(5(E11.5,' ')) 
 15   FORMAT(5(E11.5,' ')) 
 
      RETURN 
      END 

c************************************************************************ 
c this function takes string s1 and adds info to it, based on iProf,iSol 
      SUBROUTINE addinfo(caS1,iProf,iSol) 
 
      IMPLICIT NONE 
 
      CHARACTER*80 caS1 
      INTEGER iProf,iSol 
 
      INTEGER i1 
      CHARACTER*2 caString 
 
      i1 = 80 
 10      CONTINUE 
      IF ((caS1(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) THEN 
        i1 = i1 - 1 
        GOTO 10 
      END IF 
 
      IF (iSol .GE. 0) THEN 
        !!we are dealing with vib temp files, so we need "vt" names 
        caS1(i1+1:i1+2) = 'vt' 
        i1 = i1 + 2 
      ELSE 
        !!we are dealing with profile files, so we need "pt" names 
        caS1(i1+1:i1+2) = 'pt' 
        i1 = i1 + 2 
      END IF 
 
c now process iProf so that we end up with a right padded string  
c eg 2 ---> '2 ', 12 ---> '12' etc  
      WRITE(caString,15) iProf 
      IF (iProf .GE. 10) THEN 
        caS1(i1+1:i1+2) = caString(1:2) 
        i1 = i1 + 2 
      ELSE 
        caS1(i1+1:i1+1) = caString(2:2) 
        i1 = i1 + 1 
      END IF 
 
      IF (iSol .GE. 0) THEN 
        !!we are dealing with vib temp files, so we need to put in sol angle 
        caS1(i1+1:i1+2) = '_s' 
        i1 = i1 + 2 
        WRITE(caString,15) iSol 
        IF (iSol .GE. 10) THEN 
          caS1(i1+1:i1+2) = caString(1:2) 
          i1 = i1 + 2 
        ELSE 
          caS1(i1+1:i1+1) = caString(2:2) 
          i1 = i1 + 1 
        END IF 
      END IF 
 
      caS1(i1+1:i1+4) = '.prf' 
 
 15   FORMAT(I2)  
 
      RETURN 
      END 

c************************************************************************ 
c               LAPACK routines
c************************************************************************

c************************************************************************
c                          POLYNOM ROUTINES 
c************************************************************************
c this subroutine computes the AIRS pressure levels
c see /asl/packages/klayers/Doc/sci.txt 

c Layers and Levels:
c -----------------
c As we use the terms, we mean different things when we speak of "layers" and
c "levels".  A "level" refers to a point value, while a "layer" refers to a
c finite thickness slab value.  We define our layers using levels as layer
c boundaries. For AIRS, the definition is:
c     Plev(i) = exp( (7/2) * ln( A*i^2 + B*i + C ) )
c
c where
c   i refers to the level number (an integer counter)
c   and A, B, C are constants which obey the relation (exact):
c      Plev(i=1)=1100, Plev(i=38)=300, Plev(i=101)=5.0E-3 mb
c   with approximate values:
c      A = -1.5508E-4
c      B = -5.5937E-2
c      C =  7.4516

c WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
c and so can only be trusted to iI = 103

      SUBROUTINE airs_levels(raPressLevels,raAirsLevels)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

      REAL raAirsLevels(kNLTEProfLayer)
      REAL raPressLevels(kProfLayerP1)

      INTEGER iI
      REAL A,B,C
 
      A = -1.5508E-4
      B = -5.5937E-2
      C =  7.4516

c WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
c and so can only be trusted to iI = 103
      
      DO iI = 1,kNLTEProfLayer
        raAirsLevels(iI) = exp( (7.0/2) * log(A*iI*iI + B*iI + C ) )       
      END DO

c WARNING it bogs down for iI >= 104 ie it BARELY goes above the AIRS levels
c and so can only be trusted to iI = 103

      RETURN
      END

c************************************************************************
c this finds quadratic coeffs a b c   for y = ax2 + bx + c
      SUBROUTINE quadratic_coeffs(raX,raY,iMidPt,iLogOrLinear,rA,rB,rC)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      REAL raX(*),raY(*)   !!arrays to fit
      INTEGER iMidPt       !!fit about this point
      INTEGER iLogOrLinear !! +1 for raX to be log, -1 to be linear
c output params
      REAL rA,rB,rC

c local vars
      REAL rP0,rPp1,rPm1
      REAL rT0,rTp1,rTm1
      REAL rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr

      rPp1 = raX(iMidPt+1)           
      rP0  = raX(iMidPt)           
      rPm1 = raX(iMidPt-1)           

      rTp1 = raY(iMidPt+1)           
      rT0  = raY(iMidPt)           
      rTm1 = raY(iMidPt-1)           

c now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n)) 
      IF (iLogOrLinear .EQ. +1) THEN
        rP0  = alog(rP0) 
        rPp1 = alog(rPp1) 
        rPm1 = alog(rPm1) 
      END IF
        
      rDp1 = rTp1-rT0 
      rDm1 = rTm1-rT0 
 
      rp1    = rPp1-rP0 
      rp1sqr = (rPp1-rP0)*(rPp1+rP0) 
      rm1    = rPm1-rP0 
      rm1sqr = (rPm1-rP0)*(rPm1+rP0) 
 
      rA = (rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1) 
      rB = rDp1/rp1-rA*(rp1sqr/rp1) 
      rC = rT0-rA*rP0*rP0-rB*rP0 
 
c finally compute rT 
c        rT=rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC 

      RETURN
      END

c************************************************************************
c this subroutine finds the closest profile, uses results of polynom fit
c (first or second or third order)  and dumps it into a file
      SUBROUTINE polynom_nlte(
     $             caOutName,raTemp,raPressLevels,iNumLayers,
     $             rSolarAngle,iRTP,iRegr,caVTFile)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input
      REAL raTemp(kProfLayer)           !!!kinetic temp profile
      REAL raPressLevels(kProfLayer+1)  !!!pressure levels
      INTEGER iNumLayers                !!!number of layers
      REAL rSolarAngle                  !!!solar angle for atm #1 
      CHARACTER*80 caOutName
      INTEGER iRTP                      !!!which iRTP prof read in 
c output
      INTEGER iRegr               !!!which regr profile this is closest to
      CHARACTER*80 caVTFile       !!!temp file name where VT profiles are
                                  !!!created and stored if not given by user

c local vars
      REAL raLayPress(kProfLayer),raLayPress3(kProfLayerP3)
      REAL raTemp3(kProfLayerP3),raAirsLevels(kNLTEProfLayer)
      REAL rMin,rA,rB,rC,rPavg

      INTEGER iI,iJ,iStart,iNumPts3,dI

      CALL airs_levels(raPressLevels,raAirsLevels)

      !!!first find the pressure layering for the kCARTA profile
      DO iI = kProfLayer-iNumLayers+1,kProfLayer
        raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1) 
        raLayPress(iI) = raLayPress(iI)/
     $                        log(raPressLevels(iI)/raPressLevels(iI+1)) 
      END DO

      !!!fill lower layers with some "realistic" increasing values
      DO iI = kProfLayer-iNumLayers,1,-1
        raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) + 
     $                   10*abs(iI-(kProfLayer-iNumLayers+1))
      END DO

      iStart = 90
 10   CONTINUE
      IF (raAirsLevels(iStart) .LT. raLayPress(kProfLayer)) THEN
        GOTO 20
      ELSEIF (iStart .LT. 102) THEN
        iStart = iStart + 1
        GOTO 10
      ELSE
        write(kStdErr,*) 'humph : cannot find'
        write(kStdErr,*) 'raAirsLevels(iStart) .LT. raLayPress(kProfLayer)'
        Call DoStop
      END IF
 20   CONTINUE

      !!!fill in the arrays for extended layers
      DO iI = 1,kProfLayer
        raLayPress3(iI) = raLayPress(iI)
        raTemp3(iI)     = raTemp(iI)
      END DO
      !!!so now can go slightly above the AIRS levels
      !!SUBROUTINE quadratic_coeffs(raX,raY,iMidPt,iLogOrLinear,A,B,C)
      CALL quadratic_coeffs(raLayPress3,raTemp3,kProfLayer-1,+1,rA,rB,rC)

      DO iI = 1,103-iStart+1
        raLayPress3(kProfLayer+iI) = raAirsLevels(iStart-1+iI)
        !!guess the kinetic temps here!
        rPavg = raLayPress3(kProfLayer+iI)
        raTemp3(kProfLayer+iI) = rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
        iNumPts3 = kProfLayer+iI
      END DO

      !!!flip raLayPress3 so it is in increasing order!
      !!!ditto raTemp3
      DO iI = 1,int(iNumPts3*1.0/2.0)
        rMin = raLayPress3(iI)
        raLayPress3(iI) = raLayPress3(iNumPts3-iI+1)
        raLayPress3(iNumPts3-iI+1) = rMin

        rMin = raTemp3(iI)
        raTemp3(iI) = raTemp3(iNumPts3-iI+1)
        raTemp3(iNumPts3-iI+1) = rMin
      END DO

      !!!also flip raLayPress so it is in increasing order!
      !!!ditto raTemp
      !!!this is really NOT necessary as we DO NOT use the arrays here
      !!!but whatever, we might use them in n_gas_wt_spectra
      DO iI = 1,int(kProfLayer*1.0/2.0)
        rMin = raLayPress(iI)
        raLayPress(iI) = raLayPress(kProfLayer-iI+1)
        raLayPress(kProfLayer-iI+1) = rMin

        rMin = raTemp(iI)
        raTemp(iI) = raTemp(kProfLayer-iI+1)
        raTemp(kProfLayer-iI+1) = rMin
      END DO

      CALL NLTE_PolyTemp(caOutName,raLayPress3,raTemp3,iNumPts3,
     $                   rSolarAngle,caVTFile,iRTP,iRegr)

      RETURN
      END

c************************************************************************
c this subroutine reads in the generic files and computes a bunch of things!
c based on imput (raLayPress,raKinetic)
      SUBROUTINE NLTE_PolyTemp(caOutName,raLayPress3,raTemp3,iNumPts3,
     $                         rSolarAngle,caVTFile,iRTP,iRegr)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      CHARACTER*80 caOutName
      REAL raTemp3(kProfLayerP3)           !!!kinetic temp profile
      REAL raLayPress3(kProfLayerP3)       !!!layering
      INTEGER iNumPts3                     !!!number of points in the arrays
      REAL rSolarAngle
      INTEGER iRegr                        !!! closest regr prof number (1:48)
      INTEGER iRTP                         !!!which iRTP prof read in 
c output param
      CHARACTER*80  caVTfile

c local vars
      CHARACTER*120 caDir,caFnamePoly
      CHARACTER*11  caExt(12)
      REAL raaaCoef(kNLTEProfLayer,7,4)
      REAL raSolar(7),raPressRegr(kNLTEProfLayer)
      REAL raaGood(kNLTEProfLayer,7)

      REAL raaCompute(kNLTEProfLayer,12),raOut(kNLTEProfLayer)

      INTEGER iI,iOrder0,iNumSolar,iNumPressPts,iCompute

      caExt(1) = 'coeffs_ss1_'
      caExt(2) = 'coeffs_ss2_'
      caExt(3) = 'coeffs_ss3_'
      caExt(4) = 'coeffs_ss4_'
      caExt(5) = 'coeffs_pp1_'
      caExt(6) = 'coeffs_pp2_'
      caExt(7) = 'coeffs_dd1_'
      caExt(8) = 'coeffs_dd2_'
      caExt(9) = 'coeffs_qv1_'
      caExt(10)= 'coeffs_qv2_'
      caExt(11)= 'coeffs_qv3_'
      caExt(12)= 'coeffs_qv4_'

      iOrder0 = 2

      caDir = kNLTEPolyFits

c stringchar = ['ss1','ss2','ss3','ss4','pp1','pp2','dd1','dd2',... 
c              'qv1','qv2','qv3','qv4']; 

      DO iCompute = 1,12
        CALL makeVT_polyname(iOrder0,iCompute,caExt,caDir,caFnamePoly)
        CALL read_VTpoly(iOrder0,caFnamePoly,iNumPressPts,iNumSolar,
     $                   raSolar,raPressRegr,raaaCoef,raaGood)
        CALL Compute_NLTE_Poly(
     $             iNumPressPts,iNumSolar,
     $             raSolar,raPressRegr,raaaCoef,raaGood,
     $             iCompute,iOrder0,raLayPress3,raTemp3,iNumPts3,rSolarAngle,
     $             raOut)
        DO iI = 1,iNumPts3
          raaCompute(iI,iCompute) = raOut(iI)
        END DO
      END DO

      CALL VTName_rtp(iRTP,caOutName,caVTFile)
      CALL OutputVTFile_Poly(caVTFile,raaCompute,raLayPress3,raTemp3,iNumPts3)

      RETURN
      END

c************************************************************************
c just dumps stuff to a file
      SUBROUTINE OutputVTFile_Poly(caVTFile,raaCompute,
     $                             raLayPress3,raTemp3,iNumPts3)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      CHARACTER*80 caVTfile,caSummaryFile
      REAL raaCompute(kNLTEProfLayer,12)
      REAL raTemp3(kProfLayerP3)           !!!kinetic temp profile
      REAL raLayPress3(kProfLayerP3)       !!!layering
      INTEGER iNumPts3
      CHARACTER*3 ca3

c local vars
      INTEGER iI,iJ,iIOUN,iVibPf,iGasID,iErr,iZero
      INTEGER iaGasID(8),iaISO(8),iaAFGL(8)
      REAL raVibCntr(8)
      REAL raX(kProfLayerP3)
      CHARACTER*70 caStr

      iVibPf = 4
      iGasID = 2
      ca3 = '***'

      DATA (iaGasID(iI),iI=1,8) /2,2,2,2,2,2,2,2/
      DATA (iaISO(iI),iI=1,8)   /1,2,3,4,1,2,1,2/
      DATA (iaAFGL(iI),iI=1,8)/     
     $              9,          9,          9,          9,
     $              16,         16,         24,         24/
      DATA (raVibCntr(iI),iI=1,8)/  
     $             2349.143300, 2283.48800, 2332.11300, 2340.01400,
     $             3004.01200,  2920.23900, 3659.27300, 3580.75000/

      iI = 80
 100  CONTINUE
      IF (caVTfile(iI:iI) .EQ. ' ') THEN
        iI = iI - 1
        GOTO 100
      END IF
      caSummaryFile(1:iI) = caVTfile(1:iI)
      caSummaryFile(iI+1:iI+4) = '.sss'
      iIOUN = kTempUnit
      !!dump out stuff so it looks like the .sss files made from 48 regr profs
      OPEN(UNIT=iIOun,FILE=caSummaryFile,FORM='formatted',STATUS='UNKNOWN',  
     $     IOSTAT=iErr)   
      IF (iErr .NE. 0) THEN   
        write(kStdErr,*) 'trying to dump NLTE summary '
        write(kStdErr,*) 'error opening polynom NLTE file ',iErr,caSummaryFile
        CALL DoStop
      ELSE
        write(kStdWarn,*) 'dumping summary NLTE guesses to ',caVTFile
        kTempUnitOpen = +1
        DO iI = 1,iNumPts3
          write(iIOUN,111) iI,raLayPress3(iI),raTemp3(iI),
     $  raaCompute(iI,9),raaCompute(iI,10),raaCompute(iI,11),raaCompute(iI,12),
     $  raaCompute(iI,1),raaCompute(iI,2),raaCompute(iI,3),raaCompute(iI,4),
     $  raaCompute(iI,5),raaCompute(iI,6),raaCompute(iI,7),raaCompute(iI,8)
        END DO
        CLOSE(iIOUN)
        kTempUnitOpen = +1
      END IF   

      iIOUN = kTempUnit
      OPEN(UNIT=iIOun,FILE=caVTFile,FORM='formatted',STATUS='UNKNOWN',  
     $     IOSTAT=iErr)   
      IF (iErr .NE. 0) THEN   
        write(kStdErr,*) 'trying to dump NLTE guesses '
        write(kStdErr,*) 'error opening NLTE guesses file ',iErr,caVTFile
        CALL DoStop
      ELSE
        write(kStdWarn,*) 'dumping NLTE guesses to ',caVTFile
      END IF   
      kTempUnitOpen = +1

      !! ca3  ModelNumber   NumOfVibStates in file NumOfGases  NumLevs
      write(iIOUN,99) ca3,1,8,1,iNumPts3
      !! GasID NumOfVibPartFcn
      write(iIOUN,*) 2,4
      DO iI = 1,iNumPts3
        raX(iI) = raLayPress3(iNumPts3-iI+1)   !!!flip
      END DO
      CALL printarray(iIOUN,iNumPts3,raX)  !!!print out press levels 

      write(iIOUN,*) iZero       !dummy 
      DO iI = 1,iNumPts3
        raX(iI) = raTemp3(iNumPts3-iI+1)   !!!flip
      END DO
      CALL printarray(iIOUN,iNumPts3,raX)  !!!print out kinetic temps 
 
      caStr = '!QV   Vibrational Partition functions ' 
      write(iIOUN,13) caStr 
      DO iI = 1,iVibPf 
        DO iJ = 1,iNumPts3
          raX(iJ) = raaCompute(iNumPts3-iJ+1,9+(iI-1))
        END DO
        write(iIOUN,*) iGasID,iI 
        CALL printarray(iIOUN,iNumPts3,raX)  !!!print out part fcn 
      END DO 

      caStr = '!TV   Vibrational Temperatures' 
      write(iIOUN,13) caStr 
      DO iI = 1,8
        DO iJ = 1,iNumPts3
          raX(iJ) = raaCompute(iNumPts3-iJ+1,iI)
        END DO
        write(iIOUN,20) iI,iaGasID(iI),iaISO(iI),iaAFGL(iI),raVibCntr(iI)
        CALL printarray(iIOUN,iNumPts3,raX)  !!!print out NLTE VibTemp
      END DO 

      CLOSE(iIOUN)
      kTempUnitOpen = -1

 10   FORMAT(A80) 
 11   FORMAT('! ',A70) 
 12   FORMAT(A70,I5) 
 13   FORMAT(A70) 
 20   FORMAT(I4,'    ',I5,'  ',I7,' ',I6,'   ',F11.5) 
 99   FORMAT(A3,'   ',4(I3,' ')) 
 111  FORMAT(I4,' ',14(F12.5,' '))
 123  FORMAT(A80)
 
      RETURN
      END

c************************************************************************
c this subroutine computes the thingies
      SUBROUTINE Compute_NLTE_Poly(iNumPressPts,iNumSolar,
     $              raSolar,raPressRegr,raaaCoef,raaGood,
     $              iCompute,iOrder0,raLayPress3,raTemp3,iNumPts3,rSolarAngle,
     $              raOut)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input from the polyfit files
      INTEGER iNumPressPts,iNumSolar,iCompute
      REAL raaaCoef(kNLTEProfLayer,7,4)
      REAL raSolar(7),raPressRegr(kNLTEProfLayer)
      REAL raaGood(kNLTEProfLayer,7)
c input from user profile
      INTEGER iOrder0
      REAL raTemp3(kProfLayerP3)           !!!kinetic temp profile
      REAL raLayPress3(kProfLayerP3)       !!!layering
      INTEGER iNumPts3
      REAL rSolarAngle
c output 
      REAL raOut(kNLTEProfLayer)

c local vars
      INTEGER iI,iA,iB,iP,iX
      REAL raY1(kNLTEProfLayer),raY2(kNLTEProfLayer),raX(kNLTEProfLayer)
      REAL raaP1(kNLTEProfLayer,4),raaP2(kNLTEProfLayer,4)
      REAL raZin(kProfLayerP3),raZout(kProfLayerP3),raOut2(kNLTEProfLayer)
      REAL rSlope

      IF ((rSolarAngle .LT. 0.0) .OR. (rSolarAngle .GT. 120.0)) THEN
        write(kStdWarn,*) 'Hmmm solar ang == ',rSolarAngle
        write(kStdWarn,*) 'so setting raOut = 1.0 or to raTemp3'
        IF ((iCompute .GE. 9) .AND. (iCompute .LE. 12)) THEN
          DO iI = 1,kNLTEProfLayer
            raOut(iI) = 1.0          !!! vib part fcn
          END DO
        ELSEIF ((iCompute .GE. 1) .AND. (iCompute .LE. 8)) THEN
          DO iI = 1,iNumPts3
            raOut(iI) = raTemp3(iI)          !!! kinetic temp
          END DO
          DO iI = iNumPts3+1,kNLTEProfLayer
            raOut(iI) = raOut(iNumPts3)   !!! kinetic temp
          END DO
          GOTO 50
        END IF
      END IF

      IF ((rSolarAngle .GE. 0.0) .AND. (rSolarAngle .LE. 0.1)) THEN
        iA = 1
        iB = 1

      ELSEIF ((rSolarAngle .GE. 119.9) .AND. (rSolarAngle .LE. 120.0)) THEN
        iA = 7
        iB = 7

      ELSEIF ((rSolarAngle .GT. 0.1) .AND. (rSolarAngle .LT. 119.9)) THEN
        write(kStdWarn,*) 'solar ang == ',rSolarAngle
        write(kStdWarn,*) 'so computing raOut'
        iA = 7
 10     CONTINUE
        IF ((rSolarAngle .LT. raSolar(iA)) .AND. (iA .GT. 1)) THEN
          iA = iA - 1
          GOTO 10
        END IF

        iB = 1
 20     CONTINUE
        IF ((rSolarAngle .GT. raSolar(iB)) .AND. (iB .LT. 7)) THEN
          iB = iB + 1
          GOTO 20
        END IF
      END IF

      !!!!!!!!!!!!!!! have figured out iA,iB so proceed

      IF (iA. GT. iB) THEN
        write(kstdErr,*) 'Ooooooops! iA > iB'
        CALL DOStop
      END IF

      write(kStdWarn,*) 'solang bracket ',raSolar(iA),rSolarAngle,raSolar(iB) 

      DO iI = 1,iNumPts3
        raZin(iI) = log10(raLayPress3(iI))
      END DO

      IF (iOrder0 .EQ. 2) THEN      
        iX = 2
        DO iI = 1,kNLTEProfLayer
          raaP1(iI,iX-1) = 0.0
          raaP2(iI,iX-1) = 0.0
        END DO
      ELSE
        iX = 1
      END IF

      IF (iA. EQ. iB) THEN
        !! easy; no need to interp! 
        DO iP = iX,4
          !!! put in the interp pts
          DO iI = 1,iNumPressPts
            raX(iI)  = log10(raPressRegr(iI))
            raY1(iI) = raaaCoef(iI,iA,iP) 
          END DO
          CALL rspl(raX,raY1,iNumPressPts,raZin,raZout,iNumPts3)
          DO iI = 1,iNumPts3
            raaP1(iI,iP) = raZout(iI)
          END DO
        END DO

        DO iI = 1,iNumPts3
          raOut(iI) =  raaP1(iI,1)*(raTemp3(iI)**3) + 
     $                 raaP1(iI,2)*(raTemp3(iI)**2) + 
     $                 raaP1(iI,3)*raTemp3(iI) + raaP1(iI,4)
        END DO
        GOTO 40
      END IF

      IF (iA .NE. iB) THEN
        !! need to interp in solar ang! 
        DO iP = iX,4
          !!! put in the interp pts
          DO iI = 1,iNumPressPts
            raX(iI)  = log10(raPressRegr(iI))
            raY1(iI) = raaaCoef(iI,iA,iP) 
            raY2(iI) = raaaCoef(iI,iB,iP) 
          END DO
          CALL rspl(raX,raY1,iNumPressPts,raZin,raZout,iNumPts3)
          DO iI = 1,iNumPts3
            raaP1(iI,iP) = raZout(iI)
          END DO
          CALL rspl(raX,raY2,iNumPressPts,raZin,raZout,iNumPts3)
          DO iI = 1,iNumPts3
            raaP2(iI,iP) = raZout(iI)
          END DO
        END DO

        DO iI = 1,iNumPts3
          raOut(iI) =  raaP1(iI,1)*(raTemp3(iI)**3) + 
     $                 raaP1(iI,2)*(raTemp3(iI)**2) + 
     $                 raaP1(iI,3)*raTemp3(iI) + raaP1(iI,4)
        END DO
        DO iI = 1,iNumPts3
          raOut2(iI) =  raaP2(iI,1)*(raTemp3(iI)**3) + 
     $                  raaP2(iI,2)*(raTemp3(iI)**2) + 
     $                  raaP2(iI,3)*raTemp3(iI) + raaP2(iI,4)
        END DO

        !! now do linear interp in solar angle
        DO iI = 1,iNumPts3
          rSlope = (raOut2(iI) - raOut(iI))/(raSolar(iB)-raSolar(iA))
          raOut(iI) = rslope*(rSolarAngle-raSolar(iA)) + raOut(iI)
        END DO
      END IF

 40   CONTINUE

c      print *,' '
c      print *,' raaaCoef at iA',iA
c      DO iI = iNumPts3,1,-1
c        print *,iI,raPressRegr(iI),
c     $       raaaCoef(iI,iA,2),raaaCoef(iI,iA,3),raaaCoef(iI,iA,4)
c      END DO
c      print *,' '
c      print *,' raaP1',iA
c      DO iI = iNumPts3,1,-1
c        print *,iI,raLayPress3(iI),
c     $       raaP1(iI,1),raaP1(iI,2),raaP1(iI,3),raaP1(iI,4)
c      END DO
c      print *,' '
c      print *,' raaP2',iB
c      DO iI = iNumPts3,1,-1
c        print *,iI,raLayPress3(iI),
c     $        raaP2(iI,1),raaP2(iI,2),raaP2(iI,3),raaP2(iI,4)
c      END DO

c      print *,' '
c      print *,' raTemp,rOut'
c      DO iI = 1,iNumPts3
c        print *,iI,raTemp3(iI),raOut(iI)
c      END DO

      IF ((iCompute .GE. 1) .AND. (iCompute .LE. 8)) THEN
        !!! remember we did a polyfit of (tK vs (dt = tK - tWhatever))
        !!! so tK - dt = tWhatever
        !!! remember we did a polyfit of (tK vs (dt = tWhatever - tK))
        !!! so tK + dt = tWhatever
        DO iI = 1,iNumPts3
          IF (raLayPress3(iI) .GE. 800) THEN
            !!don't expect NLTE low in the atm
            raOut(iI) = raTemp3(iI)
          ELSE
            raOut(iI) = raTemp3(iI) + raOut(iI)
          END IF
        END DO
      END IF

      IF ((iCompute .GE. 9) .AND. (iCompute .LE. 12)) THEN
        !!! remember we did a polyfit of (tK vs QV)
        !!! so pretty much do nothing more
        DO iI = 1,iNumPts3
          IF ((raOut(iI) .LE. 0.9) .OR. (raOut(iI) .GE. 1.2)
     $         .OR. (raLayPress3(iI) .GE. 800)) THEN
            !!don't expect NLTE low in the atm
            raOut(iI) = 1.0
          END IF
        END DO
      END IF

 50   CONTINUE

      RETURN
      END

c************************************************************************
c this subroutine reads in the coeffs
      SUBROUTINE read_VTpoly(iOrder0,caFnamePoly,iNumPressPts,iNumSolar,
     $                   raSolar,raPressRegr,raaaCoef,raaGood)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      CHARACTER*120 caFnamePoly
      INTEGER iOrder0

c output params
      REAL raaaCoef(kNLTEProfLayer,7,4)
      REAL raSolar(7),raPressRegr(kNLTEProfLayer)
      REAL raaGood(kNLTEProfLayer,7)
      INTEGER iNumPressPts,iNumSolar

      INTEGER iI,iJ,iK,iS,iOrder,iIOUN,iErr
              
      write(kStdWarn,*) 'Opening file for NLTE polynom eval : '
      write(kStdWarn,*) caFnamePoly

      iIOUN = kTempUnit
      OPEN(UNIT=iIOun,FILE=caFNamePoly,FORM='formatted',STATUS='OLD',  
     $     IOSTAT=iErr)   
      IF (iErr .NE. 0) THEN   
        write(kStdErr,*) 'trying to open polynom NLTE file'
        write(kStdErr,*) 'error opening polynom NLTE file ',iErr,caFNamePoly
        CALL DoStop
      END IF   
      kTempUnitOpen = +1

      read(iIOUN,*) iOrder
      IF (iOrder .NE. iOrder0) THEN
        write(kStdErr,*) 'iOrder,iOrder0 unequal!',iOrder,iOrder0
        CALL DoStop
      END IF

      read(iIOUN,*) iNumPressPts
      IF (iNumPressPts .LT. 1 .OR. iNumPressPts .GT. kNLTEProfLayer) THEN
        write(kStdErr,*) 'Bad iNumPressPts!',iNumPressPts
        CALL DoStop
      END IF

      DO iI = 1,iNumPressPts
        read(iIOUN,*) raPressRegr(iI)
      END DO

      read(iIOUN,*) iNumSolar
      IF (iNumSolar .LT. 1 .OR. iNumSolar .GT. 7) THEN
        write(kStdErr,*) 'Bad iNumSolar!',iNumSolar
        CALL DoStop
      END IF

      DO iI = 1,iNumSolar
        read(iIOUN,*) raSolar(iI)
      END DO

      IF (iOrder0 .EQ. 2) THEN
        DO iI = 1,iNumPressPts
          DO iS = 1,iNumSolar
            raaaCoef(iI,iS,1) = 0.0
            read(iIOUN,*) raaGood(iI,iS),raaaCoef(iI,iS,2),raaaCoef(iI,iS,3),
     $                    raaaCoef(iI,iS,4)
c            IF (iS .EQ. 2) THEN 
c              print *, iI,iS,raPressRegr(iI),
c     $                  raaGood(iI,iS),raaaCoef(iI,iS,2),
c     $                  raaaCoef(iI,iS,3),raaaCoef(iI,iS,4)
c            END IF

         END DO
       END DO
      ELSEIF (iOrder0 .EQ. 3) THEN
        DO iI = 1,iNumPressPts
          DO iS = 1,iNumSolar
            read(iIOUN,*) raaGood(iI,iS),raaaCoef(iI,iS,1),raaaCoef(iI,iS,2),
     $                    raaaCoef(iI,iS,3),raaaCoef(iI,iS,4)
c            print *, iI,iS,raPressRegr(iI),raaGood(iI,iS),
c     $ raaaCoef(iI,iS,1),raaaCoef(iI,iS,2),raaaCoef(iI,iS,3),raaaCoef(iI,iS,4)
         END DO
       END DO
      END IF

      close(iIOUN)
      kTempUnitOpen = -1

      RETURN
      END

c************************************************************************
c this makes the name of the polynomial fit to read in
      SUBROUTINE makeVT_polyname(iOrder0,iCompute,caExt,caDir,caFnamePoly)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input params
      INTEGER iOrder0,iCompute
      CHARACTER*120 caDir
      CHARACTER*11 caExt(12)
c output params
      CHARACTER*120 caFnamePoly

      INTEGER iI,iLen
      CHARACTER*11 cExt

      DO iI = 1,120
        caFnamePoly(iI:iI) = ' '
      END DO

      cExt = caExt(iCompute)
      
      iLen = 120
 10   CONTINUE
      IF ((caDir(iLen:iLen) .EQ. ' ') .AND. (iLen .GT. 1)) THEN
        iLen = iLen - 1
        GOTO 10
      END IF

      DO iI = 1,iLen
        caFnamePoly(iI:iI) = caDir(iI:iI)
      END DO
      
      DO iI = 1,11
        caFnamePoly(iLen+iI:iLen+iI) = cExt(iI:iI)
      END DO

      iLen = 120
 20   CONTINUE
      IF ((caFnamePoly(iLen:iLen) .EQ. ' ') .AND. (iLen .GT. 1)) THEN
        iLen = iLen - 1
        GOTO 20
      END IF

      IF (iOrder0 .EQ. 2) THEN
        caFnamePoly(iLen+1:iLen+5) = '2.txt'
      ELSEIF (iOrder0 .EQ. 2) THEN
        caFnamePoly(iLen+1:iLen+5) = '3.txt'
      ELSE
        write(kStdErr,*) 'hmm incorrect iOrder0'
        CALL DoStop
      END IF

      RETURN
      END

c************************************************************************
c************************************************************************
c                          POLYNOM ROUTINES 
c************************************************************************

c************************************************************************
c                          PREDICTOR ROUTINES 
c************************************************************************
c this subroutine finds the closest profile and dumps it into a file
      SUBROUTINE predict_nlte(
     $             caOutName,raTemp,raPressLevels,iNumLayers,
     $             rSolarAngle,iRTP,iRegr,caVTFile)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

c input
      REAL raTemp(kProfLayer)           !!!kinetic temp profile
      REAL raPressLevels(kProfLayer+1)  !!!pressure levels
      INTEGER iNumLayers                !!!number of layers
      REAL rSolarAngle                  !!!solar angle for atm #1 
      CHARACTER*80 caOutName
      INTEGER iRTP                      !!!which of the rtp files processed
c output
      INTEGER iRegr               !!!which regr profile this is closest to
      CHARACTER*80 caVTFile       !!!temp file name where VT profiles are
                                  !!!created and stored if not given by user

c local vars
      REAL raLayPress(kProfLayer),raLayPress3(kProfLayerP3)
      REAL raTemp3(kProfLayerP3),raAirsLevels(kNLTEProfLayer)
      REAL rMin,rA,rB,rC,rPavg
      REAL raAirsPressure(kMaxLayer),raAirsTemps(kNLTEProfLayer)
      CHARACTER*4 caVT(10)
      CHARACTER*4 caQV(4)

      INTEGER iI,iJ,iStart,iNumPts3,dI,iAirsMax

      caVT(1) = 'sig1'
      caVT(2) = 'sig2'
      caVT(3) = 'sig3'
      caVT(4) = 'sig4'
      caVT(5) = 'wow1'
      caVT(6) = 'wow2'
      caVT(7) = 'pii1'
      caVT(8) = 'pii2'
      caVT(9) = 'del1'
      caVT(10)= 'del2'

      caQV(1) = 'qqv1'
      caQV(2) = 'qqv2'
      caQV(3) = 'qqv3'
      caQV(4) = 'qqv4'

      iAirsMax = 102    !!!beyond this, airs_levels = beserk

      CALL airs_levels(raPressLevels,raAirsLevels)
c      CALL closest_regr_lowertrop(raTemp,raPressLevels,iNumLayers,iRegr)

      !!! this is the pressure layering for the AIRS JPL instrument
      DO iI = 1,kMaxLayer
        raAirsPressure(iI) = raAirsLevels(iI)-raAirsLevels(iI+1)
        raAirsPressure(iI) = raAirsPressure(iI)/
     $                         log(raAirsLevels(iI)/raAirsLevels(iI+1))
      END DO

      !!!find the pressure layering for the kCARTA profile
      DO iI = kProfLayer-iNumLayers+1,kProfLayer
        raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1) 
        raLayPress(iI) = raLayPress(iI)/
     $                        log(raPressLevels(iI)/raPressLevels(iI+1)) 
      END DO

      !!!fill lower layers with some "realistic" increasing values
      DO iI = kProfLayer-iNumLayers,1,-1
        raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) + 
     $                   10*abs(iI-(kProfLayer-iNumLayers+1))
      END DO

      !!!interpolate the kintic temps (logarithmically) onto the AIRS levels
      CALL logrspl(raLayPress,raTemp,kProfLayer,
     $             raAirsLevels,raAirsTemps,iAirsMax)
      DO iI = 1,iAirsMax
        print *,iI,raAirsLevels(iI),raAirsTemps(iI)
      END DO
        
      CALL DoStop

      CALL VTName_rtp(iRTP,caOutName,caVTFile)
c      CALL OutputVTFile_Poly(caVTFile,raaCompute,raLayPress3,raTemp3,iNumPts3)

      RETURN
      END

c************************************************************************
c                          PREDICTOR ROUTINES 
c************************************************************************
