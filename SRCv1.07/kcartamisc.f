c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c******** THIS FILE CONTAINS VARIOUS USEFUL SUBROUTINES/FUNCTIONS *******
c** such as sorting, setting vertical temperature profiles, checking ****
c** kcarta.param, checking comp.param and xsec.param, splines etc *******
c************************************************************************
c***************** very few string routines in this file ****************
c************************************************************************
c this subroutine checks to see if the gasID is 1-36 or 101-102
      INTEGER FUNCTION MainGas(iGasID)

      include 'kcarta.param'

      INTEGER iGasID
      INTEGER iI,i1,i2
       
      iI=-1        !assume this is a main gas (1-36, or 101-102)

      i1=-1
      i2=-1

      IF ((iGasID .GE. 1) .AND. (iGasID .LE. kGasComp)) THEN
        i1=1  
        END IF

      IF ((iGasID .GE. kNewGasLo) .AND. (iGasID .LE. kNewGasHi)) THEN
        i2=1  
        END IF
      
      IF ((i1 .EQ. 1) .OR. (i2 .EQ. 1)) THEN
        iI=1
        END IF

      IF ((i2 .EQ. 1) .AND. kCKD .LT. 0) THEN
        write(kStdWarn,*) 'Cannot have gases 101,102 with CKD turned off'
        write(kStdErr,*) 'Cannot have gases 101,102 with CKD turned off'
        Call DoSTOP
        END IF

      MainGas=iI

      RETURN
      END

c************************************************************************
c this subroutine close units
      SUBROUTINE TheEnd

      include 'kcarta.param'

      IF (kStdkCarta .NE. 6) THEN  
        write(kStdWarn,*)'closing binary output file'  
        CLOSE(UNIT=kStdkCarta)      !close file where kCARTA binary goes to  
        kStdkCartaOpen=-1 
        END IF  
  
      IF ((kJacobian .GT. 0) .AND. (kStdJacob .NE. 6)) THEN  
        write(kStdWarn,*)'closing jacobian binary file'  
        CLOSE(UNIT=kStdJacob)         !close file where Jacob binary goes to  
        kStdJacobOpen=-1 
        END IF  
 
      IF (kFlux .GT. 0) THEN 
        write(kStdWarn,*)'closing flux binary file'  
        CLOSE(UNIT=kStdFlux)         !close file where flux binary goes to  
        kStdFluxOpen=-1 
        END IF  
  
      RETURN
      END
c************************************************************************

c this suroutine sets up the current print options
      SUBROUTINE SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,iType,
     $      iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $      iNumGases,iNpmix,iaNumLayer,iAtm)

      include 'kcarta.param'

      INTEGER iNumGases    !number of gases
      INTEGER iNpMix       !number of mix paths
      INTEGER iAtm         !surrent atmosphere
      INTEGER iaNumLayer(kMaxAtm)   !number of layers in atmospheres

      INTEGER iOutNum      !which printing set this is (1..kMaxPrint)
      INTEGER iPrinter     !what type (1,2,3)
      INTEGER iAtmPr       !if gas spectra, which gas; if MP, which MP,
                           !if radiances, which atm
      INTEGER iNp          !how many to output eg 100 gas spectra, or 50 MPs
      INTEGER iaOp(kPathsOut)  !list of paths/MPs/rads to output
      INTEGER iType        !-10 if dumb, 1 if paths, 2 if MPs, 3 if rads
c this is the printing switch,atmosphere# to print,# of layers to print, 
c   list of layers/paths to print (limited to kProfLayer for now) , and the 
c   pressures at which things are output 
      INTEGER iaPrinter(kMaxPrint),iaNp(kMaxPrint) 
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint) 

      INTEGER iDummy
     
      iPrinter=iaPrinter(iOutNum) 
      iAtmPr=iaGPMPAtm(iOutNum) 
      iNp=iaNp(iOutNum) 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 1)) THEN  
        !output all paths for gas
        iNp=kProfLayer*iNumGases 
        END IF 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 2)) THEN  
        !output all MPs
        iNp=iNpmix 
        END IF 

      IF ((iNp .LT. 0) .AND. (iType .EQ. 3)) THEN  
        !output all radiances for atmosphere
        iNp=iaNumlayer(iAtm)
        END IF 

      DO iDummy=1,iNp 
        iaOp(iDummy)=iaaOp(iOutNum,iDummy) 
        END DO 

      RETURN
      END

c************************************************************************
c this subroutine does some more initialixations
      SUBROUTINE SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix) 
 
      include 'kcarta.param'

      INTEGER iMixFileLines,iVertTempSet,iNpmix
      REAL raaMix(kMixFilRows,kGasStore)       !mixing table
 
      INTEGER iFileIDLO,iFileIDhi

c assume no *mixfil section 
      iMixFileLines=-1 
 
c the vertical temperature profile has not been set yet 
      iVertTempSet=-1 
 
c initialize the mixing table to weights of 0.0 
      iNpMix=1 
      DO iFileIDLo=1,kMixFilRows 
        DO iFileIDHi=1,kGasStore 
          raaMix(iFileIDLo,iFileIDHi)=0.0 
          END DO 
        END DO 

      RETURN
      END

c************************************************************************
c this subroutine inits file unit numbers
      SUBROUTINE InitializeFileUnits 
 
      include 'kcarta.param'

c set up file unit numbers so at run time they default to STDIN,STDOUT,STDOUT  
      kStdDriver=5 
      kStdkCarta=6 
      kStdJacob=6  
 
c set up common block parameters indicating all units closed 
      kStdErrOpen=-1 
      kStdWarnOpen=-1 
 
      kStdDriverOpen=-1 
      kStdkCartaOpen=-1 
      kStdJacobOpen =-1 
      kStdFluxOpen =-1 
 
      kCompUnitOpen=-1 
      kProfileUnitOpen=-1 
 
      kTempUnitOpen=-1 

      RETURN
      END

c************************************************************************

c this subroutine summarizs the output options
      SUBROUTINE SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp, 
     $                  raaOp,raaUserPress) 
 
      include 'kcarta.param'

c this is the printing switch,atmosphere# to print,# of layers to print, 
c   list of layers/paths to print (limited to kProfLayer for now) , and the 
c   pressures at which things are output 
      INTEGER iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint) 
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint) 
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer) 

      INTEGER iDummy,iOutnum

      write(kStdWarn,*) '# of printing options selected = ',iOuttypes 
      write(kStdWarn,*) '     index     option type      atm#  numpaths' 
      write(kStdWarn,*) '----------------------------------------------' 
      DO iDummy=1,iOuttypes 
        write(kStdWarn,*) iDummy,iaPrinter(iDummy),iaGPMPAtm(iDummy), 
     $                    iaNp(iDummy) 
        write(kStdWarn,*) 'paths to be printed : (if numpaths=-1,  
     $ print all)' 
        write(kStdWarn,*)(iaaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy)) 
        IF (iaPrinter(iDummy) .EQ. 3) THEN 
          write(kStdWarn,*)(raaOp(iDummy,iOutNum), 
     $                       iOutNum=1,iaNp(iDummy)) 
          write(kStdWarn,*)(raaUserPress(iDummy,iOutNum), 
     $                       iOutNum=1,iaNp(iDummy)) 
          END IF 
        write(kStdWarn,*) '    ' 
        END dO 
 
      RETURN
      END

c************************************************************************
c this subroutine checks the MixPath Vertical Temps
      SUBROUTINE CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp, 
     $                        iNpmix,iCO2,iaGases) 
 
      include 'kcarta.param'
 
      INTEGER iCo2              !which gas used to mimic CO2 temps
      INTEGER iNumGases         !how many gases
      INTEGER iNpMix            !number of mixed paths
      REAL raaTemp(kProfLayer,kGasStore)       !profile temp
      REAL raaMix(kMixFilRows,kGasStore)       !mixing table
      REAL raMixVertTemp(kMixFilRows)          !temperatures of MP layers
      INTEGER iaGases(kMaxGas)               !gasIDs stored in order

      INTEGER iDummy,iFileIDLo

      iCO2=-1

      IF (iNpmix .GT. 0) THEN 
c search for the CO2 gas === gasID 2 
c since the gases have to be entered in ascending order, either gas1 or gas2 
c is CO2 
        iCO2=-1 
        IF (iaGases(1) .EQ. 2) THEN 
          iCO2=1 
          write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!' 
        ELSE IF (iaGases(2) .EQ. 2) THEN 
          iCO2=2 
          write(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!' 
        ELSE !!!for some strange reason, no CO2 present 
          iCO2=1 
          write(kStdWarn,*) 'Temperature of Gas number 1 will mimic CO2!!' 
        END IF 
 
        CALL GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, 
     $                      iNpmix,iCO2) 
 
        write(kStdWarn,*) 'Checking Mixed Path Temp' 
        iFileIDLo=0 
        DO iDummy=1,iNpmix 
          IF (raMixVertTemp(iDummy) .LT. 0.0) THEN 
            write(kStdWarn,*) 'Negative MP Temp in Mixed Path',iDummy 
            iFileIDLo=iFileIDLo+1 
            END IF                  
          END DO 
        IF (iFileIDLo .GT. 0) THEN 
          write(kStdWarn,*) 'Warning! negative MP temperatures found!' 
          write(kStdErr,*) 'Warning! negative MP temperatures found!' 
          CALL DoSTOP 
          END IF 
        END IF 

 1111 FORMAT(A1) 

      RETURN
      END

c************************************************************************
c this subroutine does the command line stuff
      SUBROUTINE DoCommandLine(iMicrosoft,caDriverName,caOutName,
     $                         caJacobFile,iOutFileName)

      include 'kcarta.param'

c caDriverName is the name of the driver file to be processed 
      CHARACTER*80 caDriverName 
c caOutName is the name of the unformatted output file name 
c integer iOutFileName tells whether or not there is a driver name, or 
c dump output to Unit 6 
      INTEGER iOutFileName 
      CHARACTER*80 caOutName 
c caJacobFile is the name of the unformatted output file name for Jacobians 
      CHARACTER*80 caJacobFile 
c this tells if we have MS Product ie no command line stuff!
      INTEGER iMicroSoft
c this is the number of args
      INTEGER iargc 

      INTEGER iDummy,iError

      IF (iMicroSoft .GT. 0) THEN 
         
        !no command line options .. do this! 
        print *,'Enter (1) if standard kcarta computation ' 
        print *,'Enter (2) if kcarta + jacobian computation ' 
        read *,iDummy 
        IF ((iDummy .GT. 2) .OR. (iDummy .LT. 1)) THEN  
          write(kStdErr,*) 'more than two arguments, or none ' 
          write(kStdErr,*) 'is NOT allowed' 
          CALL DoSTOP 
          END IF 
       
        print *,'Enter driver namelist filename (enclose in quotes) : ' 
        read *,caDriverName 
        kStdDriver=kStdDriverKK 
 
        print *,'Enter output standard filename  (enclose in quotes) : ' 
        read *,caOutName 
        kStdkCarta=kStdkCartaKK 
        iOutFileName = 1 
 
        IF (iDummy .EQ. 2) THEN 
          print *,'Enter output jacobian filename  (enclose in quotes) : ' 
          read *,caJacobFile 
          kStdJacob=kStdJacobKK 
          END IF 
 
       ELSE 
         !use command line stuff 
         iDummy = iargc() 
 
         IF (iDummy .GT. 3) THEN  
           write(kStdErr,*) 'more than three arguments in command line' 
           write(kStdErr,*) 'is NOT allowed' 
           CALL DoSTOP 
           END IF 

         iOutFileName = -1         !assume no name 
         DO iError=1,iDummy 
           IF (iError .EQ. 1) THEN 
             CALL getarg(1,caDriverName) 
             IF (caDriverName(1:1) .NE. '-') THEN 
               write(kStdWarn,*) 'driver file name is',caDriverName 
               kStdDriver=kStdDriverKK 
               END IF 
             END IF 
 
           IF (iError .EQ. 2) THEN 
             CALL getarg(2,caOutName) 
             IF (caOutName(1:1) .NE. '-') THEN 
               iOutFileName = 1 
               write(kStdWarn,*) 'output file name is',caOutName 
               kStdkCarta=kStdkCartaKK 
               END IF 
             END IF 
 
           IF (iError .EQ. 3) THEN 
             CALL getarg(3,caJacobFile) 
             IF (caJacobFile(1:1) .NE. '-') THEN 
               write(kStdWarn,*) 'jacob file name is',caJacobFile 
               kStdJacob=kStdJacobKK 
               END IF 
             END IF 
           END DO 
 
        END IF 

      RETURN
      END 
 
c************************************************************************

c this subroutine stores the reference gas amts/temps etc
      SUBROUTINE StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,raQ, 
     $   raaRAmt,raaRTemp,raaRPress,raaRPartPress,raaQAirs,iGas,iaGases) 
 
      include 'kcarta.param'

c this is the gasnumber (not gasID!!!!!!!!!!!!!!)
      INTEGER iGas
      INTEGER iaGases(kMaxGas)
c these are the individual reference profiles 
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer) 
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer),raQ(kProfLayer) 
c these are the reference profiles stored in matrices 
      REAL raaQAirs(kProfLayer,kGasStore)
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore) 
      REAL raaRPress(kProfLayer,kGasStore) 
      REAL raaRPartPress(kProfLayer,kGasStore) 

      INTEGER iInt

      DO iInt=1,kProfLayer 
        raaQAirs(iInt,iGas)=raQ(iInt)       !ratio of KPROF/KAIRS amts 
        raaRAmt(iInt,iGas)=raRAmt(iInt)             !amts 
        raaRTemp(iInt,iGas)=raRTemp(iInt)           !temps 
        raaRPress(iInt,iGas)=raRPress(iInt)         !press 
        raaRPartPress(iInt,iGas)=raRPartPress(iInt) !part press 
 
        IF (raRAmt(iInt) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error reading Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iInt,' has negative amount' 
          CALL DoStop 
          END IF 
 
        IF (raRTemp(iInt) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error reading Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iInt,' has negative tempr' 
          CALL DoStop 
          END IF 
 
        IF (raRPress(iInt) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error reading Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iInt,' has negative pressure' 
          CALL DoStop 
          END IF 
 
        IF (raRPartPress(iInt) .LT. 0.0) THEN 
          WRITE(kStdErr,*) 'Error reading Ref Profile for Gas', 
     $                          iaGases(iGas) 
          WRITE(kStdErr,*) 'Layer ',iInt,' has negative part press' 
          CALL DoStop 
          END IF 
        END DO 

      RETURN
      END
c************************************************************************

c this subroutine sets the reference gas amts/temps etc
      SUBROUTINE SetReference(raRAmt,raRTemp,raRPress,raRPartPress,raQ, 
     $   raaRAmt,raaRTemp,raaRPress,raaRPartPress,raaQAirs,iGas) 
 
      include 'kcarta.param'

c this is the gasnumber (not gasID!!!!!!!!!!!!!!)
      INTEGER iGas
c these are the individual reference profiles 
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer) 
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer),raQ(kProfLayer) 
c these are the reference profiles stored in matrices 
      REAL raaQAirs(kProfLayer,kGasStore)
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore) 
      REAL raaRPress(kProfLayer,kGasStore) 
      REAL raaRPartPress(kProfLayer,kGasStore) 

      INTEGER iInt

      DO iInt=1,kProfLayer 
        raQ(iInt)=raaQAirs(iInt,iGas) 
        raRAmt(iInt)=raaRAmt(iInt,iGas) 
        raRTemp(iInt)=raaRTemp(iInt,iGas) 
        raRPress(iInt)=raaRPress(iInt,iGas) 
        raRPartPress(iInt)=raaRPartPress(iInt,iGas) 
        END DO 


      RETURN
      END 
c************************************************************************
 

c this integer function is "floor" -- assume rX > 0 
      INTEGER FUNCTION iFloor(rX) 
 
      REAL rX 
 
      INTEGER iI,iIm1,iIp1,iIm2,iIp2,iF 
       
      iF=nint(rX)

      iI=nint(rX) 
      iIm1=iI-1 
      iIm2=iI-2 
      iIp1=iI+1 
      iIp2=iI+2 
 
      IF (rX .GE. iIm2*1.0) THEN 
        iF=iIm2 
        END IF 
      IF (rX .GE. iIm1*1.0) THEN 
        iF=iIm1 
        END IF 
      IF (rX .GE. iI*1.0) THEN 
        iF=iI 
        END IF 
      IF (rX .GE. iIp1*1.0) THEN 
        iF=iIp1 
        END IF 
      IF (rX .GE. iIp2*1.0) THEN 
        iF=iIp2 
        END IF 
 
      iFloor = iF

      RETURN 
      END 
c************************************************************************ 
c this integer function is "ceil" -- assume rX > 0 
      INTEGER FUNCTION iCeil(rX) 
 
      REAL rX 
 
      INTEGER iI,iIm1,iIp1,iIm2,iIp2,iC
       
      iC=nint(rX)

      iI=nint(rX) 
      iIm1=iI-1 
      iIm2=iI-2 
      iIp1=iI+1 
      iIp2=iI+2 
 
      IF (rX .LE. iIp2*1.0) THEN 
        iC=iIp2 
        END IF 
      IF (rX .LE. iIp1*1.0) THEN 
        iC=iIp1 
        END IF 
      IF (rX .LE. iI*1.0) THEN 
        iC=iI 
        END IF 
      IF (rX .LE. iIm1*1.0) THEN 
        iC=iIm1 
        END IF 
      IF (rX .LE. iIm2*1.0) THEN 
        iC=iIm2 
        END IF 

      iCeil = iC
 
      RETURN 
      END 
c************************************************************************ 
 
c this subroutine closes the files in case of an emergency stop
      SUBROUTINE DoSTOP

      include 'kcarta.param'

      write(kStdWarn,*)'Fatal Error found : closing all units ..'

      IF ((kStdDriverOpen .EQ. 1) .AND. (kStdDriver .NE. 5)) THEN
        write(kStdWarn,*)'closing driver file'
        CLOSE(UNIT=kStdDriver)          !close driver file 
        END IF

      IF ((kStdkCartaOpen .EQ. 1) .AND. (kStdkCarta .NE. 6)) THEN
        write(kStdWarn,*)'closing binary output file'
        CLOSE(UNIT=kStdkCarta)          !close file where kCARTA binary goes to
        END IF

      IF (kCompUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing kcomp/xsec file'
        CLOSE(UNIT=kCompUnit)           !close kCompressed file/xsec data file
        END IF

      IF (kJacobian .GT. 0) THEN
        IF ((kStdJacobOpen .EQ. 1)  .AND. (kStdJacob .NE. 6)) THEN
          write(kStdWarn,*)'closing jacobian binary file'
          CLOSE(UNIT=kStdJacob)         !close file where Jacob binary goes to
          END IF
        END IF

      IF (kFlux .GT. 0) THEN
        write(kStdWarn,*)'closing flux binary file'
        CLOSE(UNIT=kStdFlux)         !close file where flux binary goes to
        END IF

      IF (kProfileUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing profile file '
        CLOSE(UNIT=kProfileUnit)       !close profile file 
        END IF

      IF (kTempUnitOpen .EQ. 1) THEN
        write(kStdWarn,*)'closing temporary param file'
        CLOSE(UNIT=kTempUnit)          !close temporary file eg comp.param
        END IF

      CLOSE(UNIT=kStdErr)             !close error log
      CLOSE(UNIT=kStdWarn)            !close warning log
 
      STOP

      RETURN
      END
c************************************************************************
c this subroutine sorts an integer array
c this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

      SUBROUTINE DoSort(iaArr,iCnt)

      include 'kcarta.param'

c iaArr = integer array to be sorted
c iCnt  = sort indices 1..iCnt of iaArr
      INTEGER iaArr(*),iCnt

c local variables
      INTEGER iTrue,iT,i1,i2

      i1=1
      iTrue=1

 60   CONTINUE
      IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
        i2=iCnt
        iTrue=-1
 70     CONTINUE
        IF (i2 .GT. i1) THEN
          IF (iaArr(i2) .LT. iaArr(i2-1)) THEN
            iTrue=1
            iT=iaArr(i2-1)
            iaArr(i2-1)=iaArr(i2)
            iaArr(i2)=iT
            END IF
          i2=i2-1
          GO TO 70
          END IF
      
        i1=i1+1
        GO TO 60
        END IF

       RETURN
       END
c************************************************************************
c this subroutine sorts a real array raP1 into ascending order or
c descending order depending on value of iUD ... ia1,ra1 are accordingly set
c NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!
c this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

      SUBROUTINE DoSortPress(ia1,ra1,raP1,iCnt,iUD)

      include 'kcarta.param'

c raArr = real array to be sorted
c iCnt  = sort indices 1..iCnt of raArr
c iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order
      REAL ra1(*),raP1(*)
      INTEGER iCnt,iUD,ia1(*)

c local variables
      INTEGER iTrue,iT,i1,i2
      REAL rT

      IF (iUD .GT. 0) THEN
        i1=1
        iTrue=1

 60     CONTINUE
        IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
          i2=iCnt
          iTrue=-1
 70       CONTINUE
          IF (i2 .GT. i1) THEN
            IF (raP1(i2) .LT. raP1(i2-1)) THEN
              iTrue=1
              rT=raP1(i2-1)
              raP1(i2-1)=raP1(i2)
              raP1(i2)=rT
              rT=ra1(i2-1)
              ra1(i2-1)=ra1(i2)
              ra1(i2)=rT
              iT=ia1(i2-1)
              ia1(i2-1)=ia1(i2)
              ia1(i2)=iT
              END IF
            i2=i2-1
            GO TO 70
            END IF
      
           i1=i1+1
           GO TO 60
           END IF

      ELSE IF (iUD .LT. 0) THEN
        DO i1=1,iCnt
          raP1(i1)=-raP1(i1)
          END DO
        i1=1
        iTrue=1

 80     CONTINUE
        IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
          i2=iCnt
          iTrue=-1
 90       CONTINUE
          IF (i2 .GT. i1) THEN
            IF (raP1(i2) .LT. raP1(i2-1)) THEN
              iTrue=1
              rT=raP1(i2-1)
              raP1(i2-1)=raP1(i2)
              raP1(i2)=rT
              rT=ra1(i2-1)
              ra1(i2-1)=ra1(i2)
              ra1(i2)=rT
              iT=ia1(i2-1)
              ia1(i2-1)=ia1(i2)
              ia1(i2)=iT
              END IF
            i2=i2-1
            GO TO 90
            END IF
      
           i1=i1+1
           GO TO 80
           END IF

        DO i1=1,iCnt
          raP1(i1)=-raP1(i1)
          END DO
        i1=1
        END IF
        
 
      RETURN
      END
c************************************************************************
c this subroutine sorts a real array into ascending order or descending order
c depending on value of iUD
c NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!
c this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

      SUBROUTINE DoSortReal(raArr,iCnt,iUD)

      include 'kcarta.param'

c raArr = real array to be sorted
c iCnt  = sort indices 1..iCnt of raArr
c iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order
      REAL raArr(*)
      INTEGER iCnt,iUD

c local variables
      INTEGER iTrue,i1,i2
      REAL rT

      IF (iUD .GT. 0) THEN
        i1=1
        iTrue=1

 60     CONTINUE
        IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
          i2=iCnt
          iTrue=-1
 70       CONTINUE
          IF (i2 .GT. i1) THEN
            IF (raArr(i2) .LT. raArr(i2-1)) THEN
              iTrue=1
              rT=raArr(i2-1)
              raArr(i2-1)=raArr(i2)
              raArr(i2)=rT
              END IF
            i2=i2-1
            GO TO 70
            END IF
      
           i1=i1+1
           GO TO 60
           END IF

      ELSE IF (iUD .LT. 0) THEN
        DO i1=1,iCnt
          raArr(i1)=-raArr(i1)
          END DO
        i1=1
        iTrue=1

 80     CONTINUE
        IF ((i1 .LT. iCnt) .AND. (iTrue .GT. 0)) THEN
          i2=iCnt
          iTrue=-1
 90       CONTINUE
          IF (i2 .GT. i1) THEN
            IF (raArr(i2) .LT. raArr(i2-1)) THEN
              iTrue=1
              rT=raArr(i2-1)
              raArr(i2-1)=raArr(i2)
              raArr(i2)=rT
              END IF
            i2=i2-1
            GO TO 90
            END IF
      
           i1=i1+1
           GO TO 80
           END IF

        DO i1=1,iCnt
          raArr(i1)=-raArr(i1)
          END DO
        i1=1
        END IF
        
 
      RETURN
      END
c************************************************************************
c this function, depending on iNp, calls a binary or a sequential search
c to find iLay in iaOp
      INTEGER FUNCTION DoOutputLayer(iLay,iNp,iaOp)

      include 'kcarta.param'

c iLay  = layer number to be looked for
c iaOp  = array containing list of layers
c iNp   = search indices 1..iNp of iaOp, to look for iLay
c      INTEGER iLay,iNp,iaOp(kPathsOut)
      INTEGER iLay,iNp,iaOp(*)

c integer functions that do the search
      INTEGER BinarySearch,SequentialSearch

      IF (iNp .LT. 16) THEN
        DoOutputLayer=SequentialSearch(iLay,iNp,iaOp)
      ELSE
        DoOutputLayer=BinarySearch(iLay,iNp,iaOp)
        END IF

      RETURN 
      END 
c************************************************************************
c this function searches for the value uses a binary search
c it returns +1 if iLay is found in iaOp(1:iNp), -1 otherwise
c refer Dale/Lilly : Pascal plus Data structures
      INTEGER FUNCTION BinarySearch(iLay,iNp,iaOp)

      include 'kcarta.param'

c iLay  = layer number to be looked for
c iaOp  = array containing list of layers
c iNp   = search indices 1..iNp of iaOp, to look for iLay
c      INTEGER iLay,iNp,iaOp(kPathsOut)
      INTEGER iLay,iNp,iaOp(*)

c local variables
      INTEGER iAns,iFound,iMp,iF,iL,iLocation,IDIV

      iAns=-1
      iFound=-1
      iF=1
      iL=iNp

 15   CONTINUE
      IF ((iF .LE. iL) .AND. (iFound .LT. 0)) THEN
        iMp=IDIV(iF+iL,2)
        IF (iaOp(iMp) .EQ. iLay) THEN
           iFound=1
        ELSE
          IF (iaOp(iMp) .GT. iLay) THEN
            iL=iMp-1
          ELSE
            iF=iMp+1
            END IF
          END IF
        GO TO 15
        END IF
      
      IF (iFound .GT. 0) THEN
        iLocation=iMp
        iAns=1
        END IF

      BinarySearch=iAns
      RETURN
      END
c************************************************************************
c this function checks to see if layer# iLay should be output
c by doing a sequential search : returns +1 if value found, -1 otherwise
      INTEGER FUNCTION SequentialSearch(iLay,iNp,iaOp)

      include 'kcarta.param'

c iLay  = layer number to be looked for
c iaOp  = array containing list of layers
c iNp   = search indices 1..iNp of iaOp, to look for iLay
c      INTEGER iLay,iNp,iaOp(kPathsOut)
      INTEGER iLay,iNp,iaOp(*)

c local variables
      INTEGER iDp,iDpC
      
      iDp=-1
      IF (iNp .LT. 0) THEN
c easy ! print all layers
        iDp=1
        END IF
      IF (iNp .GT. 0) THEN
c actually have to go thru list to see if this layer is to be output
        iDpC=1
 101    CONTINUE
        IF ((iDpc .LE. iNp)  .AND. (iaOp(iDpC) .EQ. iLay)) THEN            
          iDp=1
          END IF 
        IF ((iDpc .LE. iNp)  .AND. (iDp .LT. 0)) THEN
          iDpc=iDpc+1
          GO TO 101
          END IF
        END IF 

      SequentialSearch=iDp

      RETURN
      END
c************************************************************************

c this function checks to see if current GasID should have its d/dq saved
c if it does, the function result is WHICH gas it is in the *JACOBN wishlist
c else the function result = -1
      INTEGER FUNCTION DoGasJacob(iGasID,iaJacob,iJacob)

      include 'kcarta.param'

c iGasID   = current gasID
c iaJacob  = list of GasID's whose d/dq we want to output
c iJacob   = number of GasID's whose d/dq we want to output

      INTEGER iGasID,iJacob,iaJacob(kMaxDQ)

      INTEGER iI,iFound,iAns

      iFound=-1
      iAns=-1
      iI=1

 15   CONTINUE
      IF ((iFound .LT. 0) .AND. (iI .LE. iJacob)) THEN
c check to see if iGasID is in iaJacob
        IF (iGasID .EQ. iaJacob(iI)) THEN
          iFound=1
          iAns=iI
        ELSE 
          iI=iI+1
          GO TO 15
          END IF
        END IF
            
      DoGasJacob=iAns

      RETURN
      END

c************************************************************************
c this function checks to see if the cross section data base exists for
c gasID, between freqs rL,rH
c if rL,rH < 0 then checking basic presence of gas in data base
      INTEGER FUNCTION iCheckXsecDataBase(iGasID,rL,rH,iTagIn,iErr)

      include 'kcarta.param' 

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

      caFName=kXsecParamFile
      iIOUN=kTempUnit
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
          ELSE                         !actually uncompressing stuff
c           print *,'for GasID ',iGasID,'program says iTag  = ',iTagIN
c           print *,'while comp.param database file says iTag  = ',iTag
c           print *,'going on ... '
            END IF
          END IF

c this next option cannot exist, as during the uncompression we take entire 
c 25 cm-1 chunks at a time
c      IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
c     $(rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
c        iCheckCompDataBase=1
c        END IF

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

c modified 3/30 to include iTag

      INTEGER FUNCTION iCheckCompDataBase(iGasID,rL,rH,iTagIn,iErr)

      include 'kcarta.param'

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

      caFName=kCompParamFile
      iIOUN=kTempUnit
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
c        ELSE                         !actually uncompressing stuff
c          print *,'for GasID ',iGasID,'program says iTag  = ',iTagIN
c          print *,'while comp.param database file says iTag  = ',iTag
c          print *,'going on ... '
          END IF
        END IF

c this next option cannot exist, as during the uncompression, we take entire 
c 25 cm-1 chunks at a time
c      IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
c     $(rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
c        iCheckCompDataBase=1
c        END IF

      IF (iCheckCompDataBase .LT. 0) THEN
        GOTO 20
        END IF
      
 777  CONTINUE
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      RETURN
      END
c************************************************************************
c this subroutine initializes all the rows of the 
c (REAL) array of absorption coeffs
      SUBROUTINE InitializeReal(raaAb)

      include 'kcarta.param'

      REAL raaAb(kMaxPts,kProfLayer)

      INTEGER iLay,iFreq

      DO iLay=1,kProfLayer
        DO iFreq=1,kMaxPts
          raaAb(iFreq,iLay)=0.0
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes all the rows of the 
c (REAL) array of mixed paths
      SUBROUTINE InitializeRealMP(raaAb,iNpmix)

      include 'kcarta.param'

      REAL raaAb(kMaxPts,kMixFilRows)
      INTEGER iNpMix

      INTEGER iLay,iFreq

      DO iLay=1,iNpMix      !note : initialize only wot is necessary
        DO iFreq=1,kMaxPts
          raaAb(iFreq,iLay)=0.0
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes the (DOUBLE) array of absorption coeffs
      SUBROUTINE initialize(daaAb)

      include 'kcarta.param'

      DOUBLE PRECISION daaAb(kMaxPts,kProfLayer)

      INTEGER iLay,iFreq

      DO iLay=1,kProfLayer
        DO iFreq=1,kMaxPts
          daaAb(iFreq,iLay)=0.0
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this subroutine initializes the (DOUBLE) array of absorption coeffs
c pretty much the same as above routine
      SUBROUTINE initializeJAC(daaAb)

      include 'kcarta.param'

      DOUBLE PRECISION daaAb(kMaxPtsJac,kProfLayerJac)

      INTEGER iLay,iFreq

      DO iLay=1,kProfLayerJac
        DO iFreq=1,kMaxPtsJac
          daaAb(iFreq,iLay)=0.0
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this function does the integer DIVISION eg 100 div 3 = 33, 10 div 2 = 5
c assuming i1,i2 > 0
      INTEGER FUNCTION idiv(i1,i2)

      INTEGER i1,i2,iInt

c recall INT truncates the real number
      iInt=INT((i1*1.0)/(i2*1.0))
      idiv=iInt 
      RETURN
      END

c************************************************************************
c this function checks to see if the GAS ID/frequency combination are in
c the compressed data base or in the xsec database
      SUBROUTINE DataBaseCheck(iGasID,raFreq,iTag,iDoAdd,iErr)

      include 'kcarta.param'

c iGasID   = GAS ID to be searched for in CompDataBase or XsecDataBase
c raFreq   = wavenumber array that contains present chunk of 25 cm-1
c iErr     = errors (associated with file I/O)
c iTag     = 1,2,3 telling us wavenumber spacing 0.001 0.0025 0.005 cm-1
c iDoAdd   = -1,+1 tells us if we add on current gas to current 10000 pts
      INTEGER iGasID,iErr,iTag,iDoAdd
      REAL raFreq(kMaxPts)

      INTEGER iCheckCompDataBase,iCheckXsecDataBase

      iDoAdd=-1

c check to see if the k-comp file exists for the (GAS ID/freq) combination
      IF ((1 .LE. iGasID) .AND. (iGasID .LE. kGasComp)) THEN
        iDoAdd=
     $   iCheckCompDataBase(iGasID,raFreq(1),raFreq(kMaxPts),iTag,iErr)
        IF (iDoAdd .LT. 0) THEN
          WRITE(kStdWarn,333) iGasID,raFreq(1),raFreq(kMaxPts)
 333      FORMAT('Warning! No compressed data file for gasID = ',I3,
     $ 'between  wavenumbers ',f14.8,' and ',f14.8)
          END IF
        GOTO 2000
        END IF
 
c check to see if the xsec file exists for the (GAS ID/freq) combination
      IF ((kGasXsecLo .LE. iGasID) .AND. (iGasID .LE. kGasXsecHi)) THEN
        iDoAdd=iCheckXsecDataBase(iGasID,raFreq(1),raFreq(kMaxPts),iTag,iErr)
        IF (iDoAdd .LT. 0) THEN
          WRITE(kStdWarn,444) iGasID,raFreq(1),raFreq(kMaxPts)
 444      FORMAT('Warning! No xsec data for gasID = ',I3,
     $ ' between  wavenumbers ',f14.8,' and ',f14.8)
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
c this subroutine converts the abs coeff matrix from 
c double to single daa ---> raa
c and saves it in an overall AbsMatrix raaa
      SUBROUTINE DoSet(daaGasAbCoeff,raaaGasAbCoeff,iCount)

      include 'kcarta.param'

c iCount     = which of the iNumGases are being processed
c daaGasAb   = double precision abs coeffs, from the uncompression
c raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
      DOUBLE PRECISION daaGasAbCoeff(kMaxPtsJac,kProfLayerJac)
      REAL raaaGasAbCoeff(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      INTEGER iCount

c local variables
      INTEGER iLay,iFr

      DO iLay=1,kProfLayerJac
        DO iFr=1,kMaxPtsJac
c          raaaGasAbCoeff(iCount,iFr,iLay)=real(daaGasAbCoeff(iFr,iLay))
          raaaGasAbCoeff(iCount,iFr,iLay)=daaGasAbCoeff(iFr,iLay)
          END DO
        END DO

      RETURN
      END
c************************************************************************
c this subroutine converts the abs coeff matrix from 
c double to single daa ---> raa
      SUBROUTINE DoDtoR(daaGasAbCoeff,raaGasAbCoeff)

      include 'kcarta.param'

c daaGasAb   = double precision abs coeffs, from the uncompression
c raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaGasAbCoeff(kMaxPts,kProfLayer)

c local variables
      INTEGER iLay,iFr

      DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
          raaGasAbCoeff(iFr,iLay)=real(daaGasAbCoeff(iFr,iLay))
          END DO
        END DO

      RETURN
      END
c************************************************************************
C
       SUBROUTINE rSPLIN(XA,YA,Y2A,N,X,Y)

       include 'kcarta.param'
C
C REAL version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      Y2A : I  : DOUB arr : 2nd derivative of points
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),Y2A(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of pints X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! rSPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF
C
       A=(XA(KHI) - X)/H
       B=(X - XA(KLO))/H
C
       Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +
     $    Y2A(KHI)*(B**3 - B) )*(H**2)/6.0
C
       RETURN
       END

c************************************************************************
C
       SUBROUTINE rSPLY2(XA,YA,N,YP1,YPN,Y2A,WORK)
C
C REAL version
C      -----------------------------------------------------------------
C      Calc 2nd derivative as preperation for SPLINT spline routine
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      YP1 : I  : DOUB     : derivative of 1st point
C      YPN : I  : DOUB     : derivative of last point
C      Y2A : O  : DOUB arr : 2nd derivative array(N)
C      WORK: O  : DOUB arr : workspace array(N)
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),Y2A(*),YP1,YPN,WORK(*)
       INTEGER N
C
C      Local Variables
       INTEGER I,K
       REAL P,QN,SIG,UN
C
C      -----------------------------------------------------------------
C
C      Lower boundary
       IF (YP1 .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          Y2A(1)=0.0
          WORK(1)=0.0
       ELSE
C         Set to a specific first derivative
          Y2A(1)=-0.5
          WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/
     $       (XA(2) - XA(1)) - YP1)
       ENDIF
C
C      Decomposition loop of the tridiagonal algorithm
       DO I=2,N-1
c this is from the progas code
c          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
          P=SIG*Y2A(I-1) + 2.0
          Y2A(I)=(SIG - 1.0)/P
          WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -
     $       (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
          WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) -
     $       SIG*WORK(I-1) )/P
       ENDDO
C
C      Upper boundary
       IF (YPN .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          QN=0.0
          UN=0.0
       ELSE
C         Set to a specific first derivative
          QN=0.5
          UN=( 3.0/(XA(N) - XA(N-1)) )*( YPN -
     $       (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
       ENDIF
       Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
C
C      Assign the other 2nd derivatives using the back-substitution
C      loop of the tridiagonal algorithm
       DO K=N-1,1,-1
          Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
       ENDDO
C
       RETURN
       END

c************************************************************************
C  double precision
       SUBROUTINE LINEAR(XA,YA,N,X,Y)

       include 'kcarta.param'

C linear interpolation
C double precision version
C      -----------------------------------------------------------------
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       DOUBLE PRECISION XA(*),YA(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       DOUBLE PRECISION A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of points X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF

       A=(XA(KHI) - X)/H
       B=YA(KHI)-YA(KLO)

       Y=YA(KHI)-A*B

       RETURN
       END

c************************************************************************
c real
       SUBROUTINE LINEAR_REAL(XA,YA,N,X,Y)

       include 'kcarta.param'

C linear interpolation
C double precision version
C      -----------------------------------------------------------------
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of points X falls (bisect loop)
       KLO=1
       KHI=N

 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF

       A=(XA(KHI) - X)/H
       B=YA(KHI)-YA(KLO)

       Y=YA(KHI)-A*B

       RETURN
       END

c************************************************************************
       SUBROUTINE xlinear(XA,YA,N,XOUT,YOUT,NOUT) 
 
c this subroutine directly calls linear_real

       include 'kcarta.param' 
 
C real version 
C      ----------------------------------------------------------------- 
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y 
C      XA  : I  : DOUB arr : x array(N) in increasing order  IN 
C      YA  : I  : DOUB arr : y array(N)                      IN 
C      N   : I  : INT      : number of points in arrays 
C      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline 
C      YOUT   : O  : DOUB ARR     : y points from spline interpolation 
C      NOUT   : I  : INT          : number of points at which to spline 
C      ----------------------------------------------------------------- 
C 
C      Parameters 
       REAL XA(*),YA(*),XOUT(*),YOUT(*) 
       INTEGER N,NOUT 
 
       INTEGER I 
 
       DO I=1,NOUT 
         CALL LINEAR_REAL(XA,YA,N,XOUT(I),YOUT(I))  
         END DO 

       RETURN 
       END 
 
c************************************************************************ 
       SUBROUTINE xspl(XA,YA,N,XOUT,YOUT,NOUT) 
 
c this subroutine directly calls dsply2 and then dspline

       include 'kcarta.param' 
 
C real version 
C      ----------------------------------------------------------------- 
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y 
C      XA  : I  : DOUB arr : x array(N) in increasing order  IN 
C      YA  : I  : DOUB arr : y array(N)                      IN 
C      N   : I  : INT      : number of points in arrays 
C      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline 
C      YOUT   : O  : DOUB ARR     : y points from spline interpolation 
C      NOUT   : I  : INT          : number of points at which to spline 
C      ----------------------------------------------------------------- 
C 
C      Parameters 
       DOUBLE PRECISION XA(*),YA(*),XOUT(*),YOUT(*) 
       INTEGER N,NOUT 
 
       DOUBLE PRECISION Y2A(kMaxPts),work(kMaxPts),Yp1,yPn
       INTEGER I 

       yp1=1.0e16
       ypn=1.0e16

       CALL dSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work) 
       DO I=1,NOUT 
         CALL dSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I))  
         END DO 
 
       RETURN 
       END 
 
c************************************************************************ 
       SUBROUTINE dSPLIN(XA,YA,Y2A,N,X,Y)

       include 'kcarta.param'
C
C double precision version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      Y2A : I  : DOUB arr : 2nd derivative of points
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       DOUBLE PRECISION XA(*),YA(*),Y2A(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       DOUBLE PRECISION A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of pints X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! dSPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF
C
       A=(XA(KHI) - X)/H
       B=(X - XA(KLO))/H
C
       Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +
     $    Y2A(KHI)*(B**3 - B) )*(H**2)/6.0
C
       RETURN
       END

c************************************************************************
C
       SUBROUTINE dSPLY2(XA,YA,N,YP1,YPN,Y2A,WORK)
C
C double precision version
C      -----------------------------------------------------------------
C      Calc 2nd derivative as preperation for SPLINT spline routine
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      YP1 : I  : DOUB     : derivative of 1st point
C      YPN : I  : DOUB     : derivative of last point
C      Y2A : O  : DOUB arr : 2nd derivative array(N)
C      WORK: O  : DOUB arr : workspace array(N)
C      -----------------------------------------------------------------
C
C      Parameters
       DOUBLE PRECISION XA(*),YA(*),Y2A(*),YP1,YPN,WORK(*)
       INTEGER N
C
C      Local Variables
       INTEGER I,K
       DOUBLE PRECISION P,QN,SIG,UN
C
C      -----------------------------------------------------------------
C
C      Lower boundary
       IF (YP1 .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          Y2A(1)=0.0
          WORK(1)=0.0
       ELSE
C         Set to a specific first derivative
          Y2A(1)=-0.5
          WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/
     $       (XA(2) - XA(1)) - YP1)
       ENDIF
C
C      Decomposition loop of the tridiagonal algorithm
       DO I=2,N-1
c this is from the progas code
c          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
          P=SIG*Y2A(I-1) + 2.0
          Y2A(I)=(SIG - 1.0)/P
          WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -
     $       (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
          WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) -
     $       SIG*WORK(I-1) )/P
       ENDDO
C
C      Upper boundary
       IF (YPN .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          QN=0.0
          UN=0.0
       ELSE
C         Set to a specific first derivative
          QN=0.5
          UN=( 3.0/(XA(N) - XA(N-1)) )*( YPN -
     $       (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
       ENDIF
       Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
C
C      Assign the other 2nd derivatives using the back-substitution
C      loop of the tridiagonal algorithm
       DO K=N-1,1,-1
          Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
       ENDDO
C
       RETURN
       END

c************************************************************************
c simple function to see if integer iI is a member of a list iaSet (that
c has iElements in it)  
c output : 1 if in set, -1 if not
      INTEGER FUNCTION  InSet(iI,iaSet,iNumElements)

      include 'kcarta.param'
      INTEGER iI,iaSet(kMaxUserSet),iNumElements

      INTEGER iJ,iAns

c assume not in set
      iAns=-1

      IF ((iNumElements .LT. 1).OR.(iNumElements .GT. kMaxUserSet))THEN
        write(kStdErr,*) 'need a valid number of elements in list'
        CALL DoSTOP
        END IF

      iJ=0
 10   CONTINUE
      iJ=iJ+1
      IF (iaSet(iJ) .EQ. iI) THEN
        iAns = 1
      ELSE IF (iJ .LT. iNumElements) THEN
        GO TO 10
        END IF
  
      InSet=iAns
      RETURN
      END


c************************************************************************
c this subroutine checks that kNumkFile, as set in kcarta.param
c agrees with what one would ex[ect from kaMinFr,kaMaxFr

      SUBROUTINE Check_kaNum

      include 'kcarta.param'

      INTEGER iI,iJ
      REAL rF,rG

      DO iJ=1,kW
        rF=kMaxPts*kaFrStep(iJ)
        rG=kaBlSize(iJ)
        IF (abs(rF-kaBlSize(iJ)) .GT. kaFrStep(iJ)/2.0) THEN
          write(kStdErr,*) 'iJ= ',iJ
          write(kStdErr,*) 'kcarta.param claims kaBlSize(iJ)=',rG
          write(kStdErr,*) 'while it should be ',rF
          CALL DoSTOP
          END IF
        iI=INT((kaMaxFr(iJ)-kaMinFr(iJ))/rF)
        IF (iI .NE. kaNumKComp(iJ)) THEN
          write(kStdErr,*) 'iJ = ',iJ
          write(kStdErr,*) 'kcarta.param says that the number of '
          write(kStdErr,*) 'kCompressed files = ',kaNumKComp(iJ)
          write(kStdErr,*) 'based on folllowing parameters, there should 
     $ be ',iI
          write(kStdErr,*) kaMaxFr(iJ),kaMinFr(iJ),kMaxPts,kaFrStep(iJ)
          write(kStdErr,*) 'iI=INT((kMaxFreq-kMinFreq)/
     $(kMaxPts*kFreqStep))'
          CALL DoSTOP
          END IF
        END DO

      iI=0
      DO iJ=1,kW
        iI=iI+(kaNumkComp(iJ))
        END DO
      IF (iI .NE. kNumkCompT) THEN
        write(kStdErr,*) 'kcarta.param says kNumkCompT = ',kNumkCompT
        write(kStdErr,*) 'while it should be ',iI
        CALL DoSTOP
        END IF

      RETURN
      END
       
c************************************************************************
c this subroutine checks parameters in kcarta.param ... abort if they
c do not make sense .. the parameters are set in kcarta.param
      SUBROUTINE CheckKCARTAParameters

      include 'kcarta.param'
      include 'NewRefProfiles/outincLAY.param'

      write(kStdWarn,*) 'checking parameters in kcarta.param'
      write(kStdWarn,*) '  '
 
      CALL Check_kaNum

c do a quick check of the important parameters set by the user
      IF (kMixFilRows .LT. kProfLayer) THEN
        write(kStdErr,*) 'In kcarta.param, need '
        write(kStdErr,*) 'kMixFilRows >= kProfLayer(=',kProfLayer,')'
        write(kStdErr,*) 'please reset and retry'
        CALL DoSTOP
        END IF

      IF (MYNLAY .NE. kProfLayer) THEN
        write(kStdErr,*) 'kCARTA(kProfLayer) must be = KLAYERS(MYNLAY)'
        write(kStdErr,*) 'please reset kProfLayer and or MYNLAY'
        write(kStdErr,*) 'also check your kLAYERS run!!!!'
        CALL DoSTOP
        END IF

      IF (abs(kXsecFormat) .NE. 1) THEN
        write(kStdErr,*) 'kXsecFormat in kcarta.param must be = +/-1'
        write(kStdErr,*) 'please reset and retry'        
        CALL DoSTOP
        END IF

      write(kStdWarn,*) 'Max #of atmospheres from *RADFIL = ',kMaxAtm
      write(kStdWarn,*) 'Max #of gases from *GAS/XSCFIL = ',kGasStore
      write(kStdWarn,*) 'Max #of mixed paths *MIXFIL = ',kMixFilRows
      write(kStdWarn,*) '  '

      RETURN
      END

c************************************************************************
c set the default parameter values, for those that are not set in *PARAM
c read the parameter file to set some parameters that have optional values
      SUBROUTINE SetDefaultParams

c NOTE !!!! also double check subroutine EXTRAPAR in strings2.f
c NOTE !!!! also double check subroutine SETDEFAULTPARAMS in misc.f

      include 'kcarta.param'

c set default values here
c kLayer2Sp
c     -2     Layer transmittance            t(i)=exp(-k(i))
c     -1     Layer Optical depth            k(i)
c      1     Layer-to-Space Optical Depth   k2s(i)=sum(j=i,n)(k(j))
c      2     Layer-to-Space transmittance   t2s(i)=sum(j=i,n)exp(-k(j))
c kCKD      == -1,00,21,23,24 sets which continuum version calculation to do
c              -1 = no continuum

c kCKD      == -1,00,21,23,24 sets which continuum version calculation to do 
c              -1 : no continuum 
c    STD       01 : self, foreign   is by CKD and modified by Tobin 
c                                 ... this is the MT_CKD version of Dec 2002 
 
c    AIRS      02 : version 02    ... this is the MT_CKD version of Dec 2002, 
c                                     but CS modified by Scott Hannon Aug 2003 
c    AIRS      03 : version 03    ... this is the MT_CKD version of Dec 2002, 
c                                     but CS modified by Scott Hannon Jan 2004 
c    AIRS      04 : version 04    ... this is the MT_CKD version of Dec 2002, 
c                                     CS,CF modified by Scott Hannon Jan 2004 
c    AIRS      05 : version 05    ... this is the MT_CKD version of Dec 2002, 
c                                     CS,CF modified by Scott Hannon Jan 2004 
c                                     (so far it looks like CKD v4) 
c                                    On top of that it puts Dave Tobin's cs,cf 
c                                    between 1300-1800 cm-1 
c    AIRS      06 : version 06    ... this is the MT_CKD version of Dec 2002, 
c                                     CS,CF modified by Scott Hannon Dec 2005 
c                                     except slightly extended to intersect
c                                     multiplier=1 at AIRS module boundaries
c
c    STD       00 : version 00  
c    STD       21 : version 2.1 
c    STD       23 : version 2.3 
c    STD       24 : version 2.4 
c ---------------------------------------------------------------------------- 
c    RAL       12 : self = version 2.4, foreign = blend of 0.0 and 2.4 
c                                       from 0-1575, 1625-3000, use v2.4 
c                                       from 1575-1625 linearly blend v2.4,0.0 
c                                       using a triangle 
c    RAL       13 : self = version 2.3, foreign = dave tobin's v2.3 
c    RAL       90 : self = version 2.4, foreign = blend of 2.4,Tobin's thesis 
c                                       from 0-1300, 1900-3000, use v2.4 
c                                       from 1300-1900 use Tobins thesis 
c    RAL       50 : self, foreign       blend of 2.4, RAL data 
c                                       from 0-1300, 1900-3000, use v2.4 
c                                       from 1300-1900 use RAL data 
c                                       mst50.m used linear tempr interpolation 
c    RAL       51 : self, foreign       blend of 2.4, RAL data 
c                                       from 0-1300, 1900-3000, use v2.4 
c                                       from 1300-1900 use RAL data 
c                                       mst51.m used power tempr interpolation 
c                                       Need to be careful ... I put in Dave  
c                                       Tobin's fixes for CS in the 600-1100 
c                                       cm-1 range; see mst51_tobin.m in 
c                                       the SPECTRA directory (~0.9*ckd24)  
c    RAL       52 : self, foreign       blend of 2.4, RAL data 
c                                       from 0-1300, 1900-3000, use v2.4 
c                                       from 1300-1900 use RAL data 
c                                       mst52.m used power tempr interpolation 
c                                       and says CF(296)=CF(243); so use the 
c                                      foreign broadened 243 data to get CS243 
c    RAL       55 : self, foreign       same as 51 above, but uses 
c                                a) CS fudge factor of 0.87 for v <= 1137 cm-1 
c                                b) CS fudge factor of 3.20 for v >= 2394 cm-1 
c                               c) some fudge factors for CF in 1400-1700 cm-1 
c                                d) CKD 2.4 used upto 1400 cm-1 
c    RAL       56 : self, foreign       same as 51 above, but uses 
c                                a) John Taylor fudge between 600-1200 cm-1 
c                                   instead of Dave Tobin fudge 
c                                b) CS fudge factor of 3.20 for v >= 2394 cm-1 
c                               c) some fudge factors for CF in 1400-1700 cm-1 
c                                   these are better than the above fudges 
c                                d) CKD 2.4 used upto 1400 cm-1 
c    RAL       60 : self, foreign     is a hybrid of 51,55 done by scott 
c                                     so that bias errors are reduced, as 
c                                     a function of water column amount.  
c                                     eventually will include tuning coeffs 

c ---------------------------------------------------------------------------- 
      kCKD         = 1 

c kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
c              -1 if we just do the weighted average to find the MixVertTemps
c kLongOrShort == whether the user wants to save ALL header info (+1) or
c                 just a portion (-1)
      kLayer2Sp=-1
c      kCKD=21
      kGasTemp=-1
      kLongOrShort=1

c kJacobOutput == -1 if we output d(radiance)/dq,d(radiance)/dT
c                  0 if we output d(radiance)/dq * q, d(radiance)/dT
c                  1 if we output d(BT)/dq * q, d(BT)/dT
      kJacobOutput=1
c kFlux == -1 if we do not want flux computations
c       ==  1 if we want flux computations : output units = radiance s-1
c       ==  2 if we want flux computations : output units = kelvin s-1
      kFlux=-1

c kSurfTemp = -1.0 == want to use user supplied surface temp in *RADNCE
c              1.0 == want to use user supplied surface temp in *RADNCE as 
c                     an offset to pressure interpolated temperature
      kSurfTemp=-1.0     

c kTempJac == -2 if we only want d/dT(planck) in temperature jacobian
c             -1 if we only want d/dT((1-tau)(tau->sp)) in temp jacobian
c              0 if we want complete d/dT(planck) + d/dT(tau) in temp jac
      kTempJac=0

c the following cannot be controlled by the user using *PARAMS
c all the radiance parameters and the Jacobian parameter
c kJacobian ==  1 if analytic Jacobians are to be computed
c              -1 if we do the standard Genln2 computations w/o Jacobians
      kSolar = 1          !turn on solar
      kSolarAngle=0.0     !solar angle
      kSolarRefl=-1.0     !use (1-ems)/pi
      kThermal=0          !use fast diffusive approx
      kThermalAngle=-1.0  !use acos(3/5) in upper layers
      kThermalJacob=1     !use thermal backgnd in Jacobians

      kJacobian=-1        !do not do Jacobians
      kScatter=-1         !do not do scattering computations

      RETURN
      END 

c************************************************************************
c now check parameters in *PARAM
c this subroutine checks parameters in *PARAMS ... abort if they
c do not make sense .. 
      SUBROUTINE CheckParams

      include 'kcarta.param'

      write(kStdWarn,*) 'checking parameters (from *PARAMS) .... '
      write(kStdWarn,*) '  '
      IF ((iabs(kLayer2Sp) .NE. 1) .AND. (iabs(kLayer2Sp) .NE. 2)) THEN
        write(kStdErr,*) 'In *PARAMS, need kLayer2Sp = +/-1,+/-2'
        write(kStdErr,*) 'kLayer2Sp == do layer-to-space calc or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF
      IF (iabs(kGasTemp) .NE. 1) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kGasTemp = +/- 1'
        write(kStdErr,*) 'kGasTemp = use CO2 temperature profile or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF

c CKD    releases are 0,21,23,24 
c MT_CKD releases are 1 
c our modifications are 2,3,4,5 (derived from 1) 
c                       51,55,60 (derived from analysing RAL data) 
c and various 12,13,50,52,56 which might be gotten rid of eventually 
      IF ((kCKD .NE. -1)  
     !!!! the std CKD versions 
     $    .AND. (kCKD .NE. 0) .AND. (kCKD .NE. 1)   
     $    .AND. (kCKD .NE. 21) .AND. (kCKD .NE. 23) .AND. (kCKD .NE. 24)  
     !!!! these are the research versions from AIRS data 
     $    .AND. (kCKD .NE. 2)  .AND. (kCKD .NE. 3)  
     $    .AND. (kCKD .NE. 4)  .AND. (kCKD .NE. 5) .AND. (kCKD .NE. 6)
c     $    .AND. (kCKD .NE. 4)  .AND. (kCKD .NE. 5)   
     !!!! these are the research versions from RAL and AIRS data 
     $    .AND. (kCKD .NE. 12) .AND. (kCKD .NE. 13) 
     $    .AND. (kCKD .NE. 50) .AND. (kCKD .NE. 51) .AND. (kCKD .NE. 52) 
     $    .AND. (kCKD .NE. 55) .AND. (kCKD .NE. 56) .AND. (kCKD .NE. 60)) THEN 
        write(kStdErr,*) 'In *PARAMS, need kCKD = [-1] for no continuum OR' 
        write(kStdErr,*) '                 CKD    versions 0,21,23 or 24' 
        write(kStdErr,*) '              MT_CKD    versions 1,2,3,4,(5),6' 
        write(kStdErr,*) '       (latest AER versions = 1, released Dec 2002)' 
        write(kStdErr,*) '       Also have research version 2 and' 
        write(kStdErr,*) '                            12,13,50,51,52,55,56,60' 
        write(kStdErr,*) '       of which v2,51,55,60 are best' 
        write(kStdErr,*) '       (from a combination of RAL and AIRS data)' 
        write(kStdErr,*) 'kCKD is water continuum calculation version' 
        write(kStdErr,*) 'Please reset and retry' 
        CALL DoSTOP  
        END IF 
 
      IF (iabs(kLongOrShort) .NE. 1) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kLongOrShort = +/-1'
        write(kStdErr,*) 'kLongOrShort = print complete header or not'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF
      IF (iabs(kJacobOutput) .GT. 1) THEN
        write(kStdErr,*) 'In *PARAMS, need kJacobOutput =-1,0,+1'
        write(kStdErr,*) 'kJacobOutput = format to output Jacobians'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF
      IF ((iabs(kFlux) .NE. 1) .AND. (kFlux .NE.  2)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kFlux =-1,+1,+2'
        write(kStdErr,*) 'where kFlux = do/do not compute fluxes'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF
      IF ((abs(kSurfTemp)-1.0) .GE. 1e-5) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kSurfTemp = +/-1.0'
        write(kStdErr,*) 'where kSurfTemp tells the program how to use'
        write(kStdErr,*) 'the surftemperatures in *RADNCE'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF
      IF ((kTempJac .LT. -2) .OR. (kTempJac .GT. 0)) THEN
        write(kStdErr,*) 'In *PARAMS, program needs kTempJac=-2,-1,0'
        write(kStdErr,*) 'where kTempJac = use Planck or tau or both '
        write(kStdErr,*) 'when doing d/dT'
        write(kStdErr,*) 'Please reset and retry'
        CALL DoSTOP 
        END IF

      RETURN
      END


c************************************************************************
c this subroutine sets the mixed path effective tempertures, weighted 
c according to the mixing table
      SUBROUTINE GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp,
     $                          iNpmix,iCO2)
 
      include 'kcarta.param'

c iNumGases   = number of gases read in from *GASFIL + *XSCFIL
c iNpMix      = number of mixed paths read in from *MIXFIL
c raaTemp     = temperature profile of ONE of the gases (assume all equal)
c raMixVTemp  = computed vertical temp profile for mixing table
c raaMix      = mixing table from *MIXFIL
c iCO2        = which of the gases is CO2 .. if -1, none are CO2
      INTEGER iNumGases,iNpmix,iCO2
      REAL raMixVertTemp(kMixFilRows)
      REAL raaTemp(kProfLayer,kGasStore)
      REAL raaMix(kMixFilRows,kGasStore)

c local variables
      INTEGER iI,iJ,iL,MP2Lay
      REAL rT,rW

      DO iI=1,kMixFilRows
        raMixVertTemp(iI)=0.0
        END DO

      IF ((kGasTemp .EQ. 1) .AND. (iCO2 .GT. 0)) THEN
c user wants the CO2 profile to be the temperature profile
        DO iI=1,iNpmix
          iL=MP2Lay(iI)
          raMixVertTemp(iI)=raaTemp(iL,iCO2)
c          print *,iI,raMixVertTemp(iI),raaTemp(iL,iCO2)
          END DO
      ELSE
c calculate the weights      
        DO iI=1,iNpmix
          rT=0.0
          rW=0.0
          iL=MP2Lay(iI)
          DO iJ=1,iNumGases
            rT=rT+raaTemp(iL,iJ)*raaMix(iI,iJ)
            rW=rW+raaMix(iI,iJ)
            END DO
          rT=rT/rW
          raMixVertTemp(iI)=rT
c if the weights are set so that mixed path defines unique layers
c these temperatures should now be equal
c          print *,iI,raMixVertTemp(iI),raaTemp(iL,1)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine uses the nput file to set the file  number for the
c compressed data
c with the following three constraints : 
c   1) the wavenumbers have to be bewteen rFreqMin1 cm-1 and rFreqMax3 cm-1.
c   2) since each k-compressed file is 25 cm-1 long, (iU-iL)<=25
c   3) iL and iU have to fall in the same 25cm-1 block
c If all checks out, the relevant file block is set up, as is the
c wavenumber array 
      SUBROUTINE GetFreq(raWaves,rFrLow,rFrHigh,
     $  iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi,
     $  raBlock,iaFiles,iaTag,iaFileStep,iaList,iTotal) 

      include 'kcarta.param'

c raWaves   = frequency wavenumber array for first kCompressed file
c rFrLow    = lower frequency bound      
c rFrHigh   = upper frequency bound      
c iFileStartFrLo = lower file ID (605 -- 2805)
c iFileStartFrHi = upper file ID 
c iFileIDLo    = lower file index ID (1--kNumkFile)
c iFileIDLo    = upper file index ID 
c iaTag    tells which Tag is associated with which file  (1,2,3 for q,r,s..)
c iaFiles  tells which is the current kComp "file number" eg 630, 805 etc 
c          very useful so program easily knows r630_g1.dat etc 
c raBlock  tells the current kComp wavenumber block 
c iaFileStep tells you the current wavenumber step size*10000 
c iaList   has the final list of iTotal files that should be uncompressed
      INTEGER iaFileStep(kNumkCompT),iaList(kNumkCompT),iTotal
      INTEGER iaFiles(kNumkCompT),iaTag(kNumkCompT) 
      REAL raBlock(kNumkCompT)
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      INTEGER iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi

c local variables
      REAL rLow,rHigh,rTemp
      INTEGER iDummy

      write(kStdWarn,*) '**Setting freqs, fileID lower/upper bounds***'

      IF (rFrLow .ge. rFrHigh) THEN
c swap values
        write(kStdWarn,*) 'Swapping frequencies found in *FRQNCY'
        rTemp=rFrLow
        rFrLow=rFrHigh
        rFrHigh=rTemp
        END IF        

c remember : int(x)  === trunc(x)
c remember : nint(x) === round(x)
      rLow=rFrLow
      rHigh=rFrHigh
c first round down rFrLow, and round up rFrHigh
      rFrLow=1.0*INT(rFrLow)
      IF (abs(rFrHigh-1.0*INT(rFrHigh)) .GT. 0.000000) THEN
        rFrHigh=rFrHigh+1.0
        END IF
      rFrHigh=1.0*INT(rFrHigh)
      write(kStdWarn,*) 'rounded down/up rFrLow/rFrhigh as necessary ..'
      write(kStdWarn,*) 'the values are now : ',rFrLow,rFrHigh

      rLow=rFrLow
      rHigh=rFrHigh

      IF (rLow .ge. rHigh) THEN
c swap values
        rTemp=rLow
        rLow=rHigh
        rHigh=rTemp
        END IF        

      CALL filebounds(raWaves,rLow,rHigh,
     $             iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi,
     $             raBlock,iaFiles,iaTag,iaFileStep,iaList,iTotal)     
      rFrLow=rLow
      rFrHigh=rHigh
      write(kStdWarn,*) 'rFrLow,rFrHigh after = ',rFrLow,rFrHigh

      IF (iaTag(iaList(1)) .NE. iaTag(iaList(iTotal))) THEN
        write(kStdWarn,*) 'Start file tag = ',iaTag(iaList(1))
        write(kStdWarn,*) 'Stop  file tag = ',iaTag(iaList(iTotal))

c        print *,'iTotal = ',iTotal
c        print *,(iaList(iDummy),iDummy=1,iTotal)
c        print *,' '
c        print *,(iaFiles(iaList(iDummy)),iDummy=1,iTotal)
c        print *,' '
c        print *,(iaTag(iaList(iDummy)),iDummy=1,iTotal)

        write(kStdErr,*) 'program requires you choose start/stop freqs'
        write(kStdErr,*) 'that only span one wavenumber spacing, as '
        write(kStdErr,*) 'defined at the bottom of kcarta.param'
        DO iDummy = 1,kW
          write(kStdErr,333) kaMinFr(iDummy),kaMaxFr(iDummy),
     $                 kaFrStep(iDummy),kaTag(iDummy)
          END DO
        write(kStdErr,*) 'please reset *FRQNCY and retry'
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
      SUBROUTINE filebounds(raWaves,rL,rH,
     $             iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi,
     $             raBlock,iaFiles,iaTag,iaFileStep,iaList,iTotal) 


      include 'kcarta.param'

c rL,rH     = lower/upper frequency wavenumber bounds : could be reset here
c raWaves   = frequency wavenumbers for the first kcomp chunk to be processed

c iFileStartFrLo = from rL, lower file ID bound (in wavenumbers e.g. 605)
c iFileStartFrHi = from rH, upper file ID bound (in wavenumbers e.g. 2780)
c iFileIDLo    = from rL, lower fileID (e.g. 605 ==> 1+kNumKfile(1))
c iFileIDHi    = from rH, upper fileID(e.g. 2780 ==>kNumkFile(1)+
c                                                            kNumKFile(2))

c iaTag    tells which Tag is associated with which file  (1,2,3 for q,r,s..)
c iaFiles  tells which is the current kComp "file number" eg 630, 805 etc 
c          very useful so program easily knows r630_g1.dat etc 
c iaFileStep tells you the current wavenumber step size*10000 
c raBlock  tells the current kComp wavenumber block (ie start freq)

c iaList   has the final list of iTotal files that should be uncompressed
      INTEGER iaFileStep(kNumkCompT),iaList(kNumkCompT),iTotal
      INTEGER iaFiles(kNumkCompT),iaTag(kNumkCompT) 
      REAL raBlock(kNumkCompT),rL,rH,raWaves(kMaxPts)
      INTEGER iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi

c local variables
      INTEGER iFileStartFrLo1,iFileStartFrHi1,iFileIDLo1,iFileIDHi1
      INTEGER iFileStartFrLo2,iFileStartFrHi2,iFileIDLo2,iFileIDHi2
      REAL raBlockEnd(kNumkCompT),rL1,rH1,rL2,rH2,rEnd
      INTEGER iInt,iFound,iMax,iTruthLo,iTruthHi
      INTEGER iErr,iDummy

c first check lower, upper limits
      CALL CheckLimits(rL,rH)

      iMax=0
      DO iDummy=1,kW            !compute total number of kCOMP chunks present
        iMax=iMax+kaNumkComp(iDummy)
        END DO

c note there could be "overlaps" between the q r s files eg
c     < 1=500-550> <2=550-600>  <3=600-650> <4=650-700>
c                                  <5=605-630> <6=630-655> <7=655-680> .....
      iDummy=0 
      DO iFound=1,kW 
        DO iInt=1,kaNumkComp(iFound) 
          iDummy=iDummy+1 
          !this is needed by kcartamain.f
          raBlock(iDummy)=kaMinFr(iFound)+(iInt-1)*kaBlSize(iFound) 
          iaFiles(iDummy)=NINT(raBlock(iDummy)) 
          iaTag(iDummy)=kaTag(iFound) 
          iaFileStep(iDummy)=kaBlSize(iFound) 

          !this is needed within misc.f
          raBlockEnd(iDummy)=raBlock(iDummy)+kaBlSize(iFound)     
          raBlockEnd(iDummy)=raBlockEnd(iDummy)-kaFrStep(iFound)  
          END DO 
        END DO 

c go thru raBlock and see where we think the start file ID, stop file ID
c should be set at
      CALL LowerLimits(iMax,rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $     iFileStartFrLo1,iFileStartFrLo2,raBlock,raBlockEnd,iaFiles)
      CALL UpperLimits(iMax,rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $     iFileStartFrHi1,iFileStartFrHi2,raBlock,raBlockEnd,iaFiles)

c now set the fileID's to be used
      CALL LowerFileID(iMax,rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $  iFileStartFrLo1,iFileStartFrLo2,raBlock,raBlockEnd,
     $  iaTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,iFileStartFrLo)
      CALL UpperFileID(iMax,rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $  iFileStartFrHi1,iFileStartFrHi2,raBlock,raBlockEnd,
     $  iaTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,iFileStartFrHi)

      IF ((iFileIDLo .GE. 1) .AND. (iFileIDLo .LE. iMax) .AND.
     $    (iFileIDHi .GE. 1) .AND. (iFileIDHi .LE. iMax) .AND.
     $    (iFileIDLo .LE. iFileIDHi)  .AND. 
     $    (iTruthLo .GT. 0) .AND. (iTruthHi .GT. 0)) THEN
        iErr=-1
      ELSE
        write(kStdErr,*) 'Error in setting file '
        write(kStdErr,*) 'lower/upper bounds from *FRQNCY'
        CALL DoSTOP
        END IF

c now create iaList as necessary
c remember that iFileID is basically an index setting equal to iDummy above
c thus if we know iFileID we know everything!!!!!
      IF (iaTag(iFileIDLo) .EQ. iaTag(iFileIDHi)) THEN  
        !very easy everything is in either q or r or s database
        iTotal=0
        DO iInt=1,(iFileIDHi-iFileIDLo+1)
          iTotal=iTotal+1
          iaList(iTotal)=iFileIDLo+(iInt-1)     !save file ID
          write(kStdWarn,*) iTotal,iaTag(iaList(iTotal)),
     $          raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
          END DO
      ELSE
        !very hard : mixing of q,r,s databases
        iTotal=0
        iTotal=iTotal+1
        iaList(iTotal)=iFileIDLo                !save file ID
        rEnd=raBlockEnd(iaList(iTotal))
        write(kStdWarn,*) iTotal,iaTag(iaList(iTotal)),
     $              raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
        DO iInt=2,(iFileIDHi-iFileIDLo+1)
          IF (rEnd .LT. raBlockEnd(iFileIDLo+iInt-1)) THEN
            iTotal=iTotal+1
            iaList(iTotal)=iFileIDLo+(iInt-1)   !save file ID
            write(kStdWarn,*) iTotal,iaTag(iaList(iTotal)),
     $             raBlock(iaList(iTotal)),raBlockEnd(iaList(iTotal))
            rEnd=raBlockEnd(iaList(iTotal))
            END IF
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine checks to make sure rL,rH lie within the lower/upper bounds
c of wavenumbers in the kComp files
      SUBROUTINE CheckLimits(rL,rH)

      include 'kcarta.param'

c rl,rH are the lower and upper limits the user sends in; this subroutine 
c                       can change them

      REAL rL,rH
      CHARACTER cAns

      cAns='y'

      IF (rL .LT. kaMinFr(1)) THEN
        write(kStdWarn,*) 'Error!!Setting min wavenumber to ',kaMinFr(1)
ccc        write(kStdErr,*) 'Error!!Setting min wavenumber to ',kaMinFr(1)
        rL=kaMinFr(1)
c        PRINT *,'Do you wish to continue? (yY/nN) '
c        READ(5,1111)cAns
c        IF ((cAns .EQ. 'n') .OR. (cAns .EQ. 'N')) THEN
c          CALL DoSTOP
c          END IF
        END IF
      IF (rL .GT. (kaMaxFr(kW)-kaBlSize(kW))) THEN
        write(kStdWarn,*) 'Error!!Setting minimum wavenumber to ',
     $        kaMaxFr(kW)-kaBlSize(kW)
ccc        write(kStdErr,*) 'Error!!Setting minimum wavenumber to ',
ccc     $        kaMaxFr(kW)-kaBlSize(kW)
        rL=kaMaxFr(kW)-kaBlSize(kW)
c        PRINT *,'Do you wish to continue? (yY/nN) '
c        READ(5,1111)cAns
c        IF ((cAns .EQ. 'n') .OR. (cAns .EQ. 'N')) THEN
c          CALL DoSTOP
c          END IF
        END IF
      IF (rH .GT. kaMaxFr(kW)) THEN
        write(kStdWarn,*)'Error!!Setting max wavenumber to ',kaMaxFr(kW)
ccc        write(kStdErr,*)'Error!!Setting max wavenumber to ',kaMaxFr(kW)
        rH=kaMaxFr(kW)
c        PRINT *,'Do you wish to continue? (yY/nN) '
c        READ(5,1111)cAns
c        IF ((cAns .EQ. 'n') .OR. (cAns .EQ. 'N')) THEN
c          CALL DoSTOP
c          END IF
        END IF
      IF (rH .LT. (kaMinFr(1)+kaBlSize(1))) THEN
        write(kStdWarn,*) 'Error!!Setting maximum wavenumber to ',
     $         kaMinFr(1)+kaBlSize(1)
ccc        write(kStdErr,*) 'Error!!Setting maximum wavenumber to ',
ccc     $         kaMinFr(1)+kaBlSize(1)
        rH=kaMinFr(1)+kaBlSize(1)
c        PRINT *,'Do you wish to continue? (yY/nN) '
c        READ(5,1111)cAns
c        IF ((cAns .EQ. 'n') .OR. (cAns .EQ. 'N')) THEN
c          CALL DoSTOP
c          END IF
        END IF

 1111 FORMAT(A1)

      RETURN
      END

c************************************************************************
c this subroutine sets the lower file ID limits, by going thru the list top
c to bottom and bottom to top
      SUBROUTINE LowerLimits(iMax,rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $ iFileStartFrLo1,iFileStartFrLo2,raBlock,raBlockEnd,iaFiles)

      include 'kcarta.param'

c INPUT
c iMax       = max number of kcomp files (sum of files in q r s)
c rL         = user set lower bound
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaFiles    = int(file start freq)
c OUTPUT
c rL1,rL2    = freq in which rL is set down to, wne checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c iFileStartFrLo1 = file lower freq where  rL is found, when checking 
c                   list from t->b
c iFileStartFrLo2 = file lower freq where  rL is found, when checking 
c                    list from b->t
      INTEGER iaFiles(kNumkCompT)
      INTEGER iMax,iFileIDLo1,iFileIDLo2,iFileStartFrLo1,iFileStartFrLo2
      REAL rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

      INTEGER iInt

c now check which block rL falls in, starting from the extreme highest file ID
c original code
      iInt=iMax
      rL1=rL
 11   CONTINUE
      IF ((iInt .GT. 1).AND.(rL1 .LT. raBlock(iInt))) THEN
        iInt=iInt-1
        GO TO 11
        END IF
      iFileIDLo1=iInt
      iFileStartFrLo1=iaFiles(iInt)
      rL1=raBlock(iInt)
      write(kStdWarn,*) 'set lower freq to start of kcomp block',rL1
ccc      write(kStdErr,*) 'set lower freq to start of kcomp block',rL1

c now check which block rL falls in, starting from the extreme lowest file ID
c new code
      rL2=rL
      iInt=1
 12   CONTINUE
      IF ((iInt .LT. iMax).AND.(rL2 .GT. raBlockEnd(iInt))) THEN
        iInt=iInt+1
        GO TO 12
        END IF
      iFileIDLo2=iInt
      iFileStartFrLo2=iaFiles(iInt)
      rL2=raBlock(iInt)
      write(kStdWarn,*) 'reset lower freq to start of kcomp block',rL2
ccc      write(kStdErr,*) 'reset lower freq to start of kcomp block',rL2

      RETURN
      END
c************************************************************************
c this subroutine sets the upper file ID limits, by going thru the list top
c to bottom and bottom to top
      SUBROUTINE UpperLimits(iMax,rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $     iFileStartFrHi1,iFileStartFrHi2,raBlock,raBlockEnd,iaFiles)

      include 'kcarta.param'

c INPUT
c iMax       = max number of kcomp files (sum of files in q r s)
c rH         = user set lower bound
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaFiles    = int(file start freq)
c OUTPUT
c rH1,rH2    = freq in which rL is set doiwn to, wne checking list t->b, b->t
c iFileIDHi1    = file ID in which rL is found, when checking list from t->b
c iFileIDHi2    = file ID in which rL is found, when checking list from b->t
c iFileStartFrHi1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c iFileStartFrHi2 = file lower freq where  rL is found, when checking list 
c                   from b->t
      INTEGER iaFiles(kNumkCompT)
      INTEGER iMax,iFileIDHi1,iFileIDHi2,iFileStartFrHi1,iFileStartFrHi2
      REAL rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)

      INTEGER iInt

c now check which block rH falls in, starting from the extreme lowest file ID
c original code
      iInt=1
      rH1=rH
 22   CONTINUE
      IF ((iInt .LT. iMax).AND.(rH1 .GT. raBlock(iInt+1))) THEN
        iInt=iInt+1
        GO TO 22
        END IF
      iFileIDHi1=iInt
      iFileStartFrHi1=iaFiles(iInt)
      rH1=raBlockEnd(iInt)
      write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH1
ccc      write(kStdErr,*) 'reset upper freq to end of kcomp block',rH1

c now check which block rH falls in, starting from the extreme highest file ID
c new code
      iInt=iMax
      rH2=rH
 23   CONTINUE
      IF ((iInt .GT. 1).AND.(rH2 .LE. raBlock(iInt))) THEN
        iInt=iInt-1
        GO TO 23
        END IF
      iFileIDHi2=iInt
      iFileStartFrHi2=iaFiles(iInt)
      rH2=raBlockEnd(iInt)
      write(kStdWarn,*) 'reset upper freq to end of kcomp block',rH2
ccc      write(kStdErr,*) 'reset upper freq to end of kcomp block',rH2

      RETURN
      END

c************************************************************************
c this subroutine sets which file ID should be used as the lower limit
c remember tags depend on q r s ... = 1 2 3 .....
      SUBROUTINE LowerFileID(iMax,rL,rL1,rL2,iFileIDLo1,iFileIDLo2,
     $  iFileStartFrLo1,iFileStartFrLo2,raBlock,raBlockEnd,
     $  iaTag,iFileIDHi1,iFileIDHi2,iTruthLo,iFileIDLo,iFileStartFrLo)

      include 'kcarta.param'

c OUTPUT
c iFileIDLo     = final iFileIDLo
c iFinalIDLo = final iFileStartFrLo
c iTruthLo   = have we certainly set the file ID???
c INPUT
c iMax       = max number of kcomp files (sum of files in q r s)
c rL         = user set lower bound
c rL1,rL2    = freq in which rL is set doiwn to, wne checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c iFileIDHi1    = file ID in which rH is found, when checking list from b->t
c iFileIDHi2    = file ID in which rH is found, when checking list from t->b
c iFileStartFrLo1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c iFileStartFrLo2 = file lower freq where  rL is found, when checking list 
c                   from b->t
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaTag      = 1,2,3 depending on q r s
      INTEGER iMax,iFileIDLo1,iFileIDLo2,iFileStartFrLo1,iFileStartFrLo2
      REAL rL,rL1,rL2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
      INTEGER iaTag(kNumkCompT),iTruthLo,iFileIDHi1,iFileIDHi2
      INTEGER iFileIDLo,iFileStartFrLo

      INTEGER iDiff11,iDiff12,iDiff21,iDiff22

      iTruthLo=-1
      IF (iFileIDLo1 .EQ. iFileIDLo2) THEN !lower limit simple:equal file ID's
         iFileIDLo=iFileIDLo1
         iFileStartFrLo=iFileStartFrLo1
         iTruthLo=1
         rL=rL1
         write(kStdWarn,*) 'Equal lower limit : iFileIDLo = ',iFileIDLo
         END IF

      IF (iFileIDLo1 .NE. iFileIDLo2) THEN    !lower limit hard
        write(kStdWarn,*) 'Unequal lower limit : iFileIDLo1,2 = ',
     $                                  iFileIDLo1,iFileIDLo2
         write(kStdWarn,*)'thinking ... '

        IF (iaTag(iFileIDHi1) .EQ. iaTag(iFileIDHi2)) THEN
          !for UPPER limit, equal tags; so now see if either of the 
          !lower iFileIDLo1 or iFileIDLo2 fall in the same tag block
          IF (iaTag(iFileIDLo1) .EQ. iaTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE IF (iaTag(iFileIDLo2) .EQ. iaTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
          ELSE IF ((iaTag(iFileIDHi2)-iaTag(iFileIDLo1)) .LT. 
     $               (iaTag(iFileIDHi2)-iaTag(iFileIDLo2))) THEN
            !look for min diff between tags
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE
            !look for min diff between tags
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
            END IF
          END IF

        IF (iaTag(iFileIDHi1) .NE. iaTag(iFileIDHi2)) THEN
          !for UPPER limit, unequal file tags
          !so look for minimum difference between TAGS
          iDiff11=iaTag(iFileIDHi1)-iaTag(iFileIDLo1)
          iDiff12=iaTag(iFileIDHi1)-iaTag(iFileIDLo2)
          iDiff21=iaTag(iFileIDHi2)-iaTag(iFileIDLo1)
          iDiff22=iaTag(iFileIDHi2)-iaTag(iFileIDLo2)

          IF (iaTag(iFileIDLo1) .EQ. iaTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE IF (iaTag(iFileIDLo2) .EQ. iaTag(iFileIDHi2)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
          ELSE IF (iaTag(iFileIDLo1) .EQ. iaTag(iFileIDHi1)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE IF (iaTag(iFileIDLo2) .EQ. iaTag(iFileIDHi1)) THEN
            !these two tags are the same
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
          !tags are rather different : come to the first acceptable combination
          ELSE IF (iDiff11.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE IF (iDiff21.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo1
            iFileStartFrLo=iFileStartFrLo1
            iTruthLo=1
            rL=rL1
          ELSE IF (iDiff12.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
          ELSE IF (iDiff22.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDLo = iFileIDLo2
            iFileStartFrLo=iFileStartFrLo2
            iTruthLo=1
            rL=rL2
            END IF
          END IF
        END IF

      iFileIDlo1=iFileIDlo
      iFileIDlo2=iFileIDlo
      iFileStartFrLo1=iFileStartFrLo
      iFileStartFrLo2=iFileStartFrLo

      write(kStdWarn,*) 'Low bound:itruthlo,iFileIDlo,iFileStartFrLo= '
      write(kStdWarn,*) itruthlo,iFileIDlo,iFileStartFrLo

      RETURN
      END

c************************************************************************
c this subroutine sets which file ID should be used as the upper limit
      SUBROUTINE UpperFileID(iMax,rH,rH1,rH2,iFileIDHi1,iFileIDHi2,
     $   iFileStartFrHi1,iFileStartFrHi2,raBlock,raBlockEnd,
     $   iaTag,iFileIDLo1,iFileIDLo2,iTruthHi,iFileIDHi,iFileStartFrHi)

      include 'kcarta.param'

c OUTPUT
c iFileIDHi     = final iFileIDHi
c iFinalIDHi = final iFileStartFrHi
c iTruthHi   = have we certainly set the file ID???
c INPUT
c iMax       = max number of kcomp files (sum of files in q r s)
c rL         = user set lower bound
c rL1,rL2    = freq in which rL is set doiwn to, wne checking list t->b, b->t
c iFileIDLo1    = file ID in which rL is found, when checking list from t->b
c iFileIDLo2    = file ID in which rL is found, when checking list from b->t
c iFileIDHi1    = file ID in which rH is found, when checking list from b->t
c iFileIDHi2    = file ID in which rH is found, when checking list from t->b
c iFileStartFrLo1 = file lower freq where  rL is found, when checking list 
c                   from t->b
c iFileStartFrLo2 = file lower freq where  rL is found, when checking list 
c                   from b->t
c raBlock    = file start freqs
c raBlockEnd = file stop freqs
c iaTag      = 1,2,3 depending on q r s
      INTEGER iMax,iFileIDHi1,iFileIDHi2,iFileStartFrHi1,iFileStartFrHi2
      REAL rH,rH1,rH2,raBlock(kNumkCompT),raBlockEnd(kNumkCompT)
      INTEGER iaTag(kNumkCompT),iTruthHi,iFileIDLo1,iFileIDLo2
      INTEGER iFileIDHi,iFileStartFrHi

      INTEGER iDiff11,iDiff12,iDiff21,iDiff22

      iTruthHi=-1
      IF (iFileIDHi1 .EQ. iFileIDHi2) THEN    !upper limit simple
        iFileIDHi=iFileIDHi1
        iFileStartFrHi=iFileStartFrHi1
        iTruthHi=1
        rH=rH1
        write(kStdWarn,*)'Equal upper limit : iFileIDHi = ',iFileIDHi
        END IF

      IF (iFileIDHi1 .NE. iFileIDHi2) THEN    !upper limit hard
        write(kStdWarn,*) 'Unequal upper limit : iFileIDHi1,2 = ',
     $            iFileIDHi1,iFileIDHi2
        write(kStdWarn,*)'thinking ... '

        IF (iaTag(iFileIDLo1) .EQ. iaTag(iFileIDLo2)) THEN
          !for LOWER limit, equal tags; so now see if either of the 
          !upper iFileIDHi1 or iFileIDHi2 fall in the same tag block
          IF (iaTag(iFileIDHi1) .EQ. iaTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE IF (iaTag(iFileIDHi2) .EQ. iaTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
          ELSE IF ((iaTag(iFileIDLo2)-iaTag(iFileIDHi1)) .LT. 
     $              (iaTag(iFileIDLo2)-iaTag(iFileIDHi2))) THEN
            !look for min diff between tags
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE
            !look for min diff between tags
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
            END IF
          END IF

        IF (iaTag(iFileIDLo1) .NE. iaTag(iFileIDLo2)) THEN
          !for LOWER limit, unequal file tags
          !so look for minimum difference between TAGS
          iDiff11=iaTag(iFileIDHi1)-iaTag(iFileIDLo1)
          iDiff12=iaTag(iFileIDHi1)-iaTag(iFileIDLo2)
          iDiff21=iaTag(iFileIDHi2)-iaTag(iFileIDLo1)
          iDiff22=iaTag(iFileIDHi2)-iaTag(iFileIDLo2)

          IF (iaTag(iFileIDHi1) .EQ. iaTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE IF (iaTag(iFileIDHi2) .EQ. iaTag(iFileIDLo2)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
          ELSE IF (iaTag(iFileIDHi1) .EQ. iaTag(iFileIDLo1)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE IF (iaTag(iFileIDHi2) .EQ. iaTag(iFileIDLo1)) THEN
            !these two tags are the same
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
         !tags are rather different : come to the first acceptable combination
          ELSE IF (iDiff11.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE IF (iDiff21.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi1
            iFileStartFrHi=iFileStartFrHi1
            iTruthHi=1
            rH=rH1
          ELSE IF (iDiff12.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
          ELSE IF (iDiff22.EQ.min(iDiff11,iDiff12,iDiff21,iDiff22)) THEN
            iFileIDHi = iFileIDHi2
            iFileStartFrHi=iFileStartFrHi2
            iTruthHi=1
            rH=rH2
            END IF
          END IF
        END IF

c now subtract delta(wavenumber) so that we don't have to do an additional
c set of calculations because of the additional wavenumber point


      write(kStdWarn,*)'High bound:itruthhi,iFileIDhi,iFileStartFrHi = '
      write(kStdWarn,*) itruthhi,iFileIDhi,iFileStartFrHi

      RETURN
      END
c************************************************************************
