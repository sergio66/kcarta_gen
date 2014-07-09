c program to read in UNFORMATTED jacobian output from kcarta.x and save to  
c easier-to-read UNFORMATTED file, or text file 
 
      IMPLICIT NONE

      include 'convolve.param'

c the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName
      CHARACTER*130 caComment
      INTEGER iIOUN,iNumLayers,iSetLow,iSetHigh,iNumPaths
      INTEGER iaGasID(kGasStore)
      INTEGER iNatm,iBinOrText
      INTEGER iStore,iLoop,ii,ii1,ii2,iJ
      INTEGER iOk,iBlSz,iNBl,iMaxSet,iFreqPts
      INTEGER posn,offset,SizePerDataBlock,status
      REAL raEntire(kMaxEntire),rFr1,rFr2,raFreq(kMaxEntire)
      INTEGER iNumPathsOut,iNumMixPathsOut,iaNumLayersOut(kMaxAtm)
      INTEGER iaArtificial(kProfLayer),iArtificial,iWhichStore
      INTEGER iTotal,iOutNumber,iImportant,iaNumOut(kMaxPrint)

      INTEGER iExtra,iNumGases,iaNumLayers(kMaxAtm)

      INTEGER iInstrType,iLpLp,iOutTypeLoop
    
 90   CONTINUE
      iInstrType = 0

      print *,'Do you want a text or binary output? (-1/1) '  
      read *,iBinOrText  

      print *,'Enter INPUT binary file name '
      read 1000,caInName
 1000 FORMAT(A80)

c first read header to figure out iStore ******************
      iIOUN=10
      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

      CALL readmainjacobheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh) 
      CALL readjacobinfo(iIOUN,iNumGases,iNatm,iaNumLayers) 
      CALL readgasjacob(iIOUN,iNumGases,100) 
      read(iIOUN) iTotal,iOutNumber,iImportant  
      read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)  
  
      print *,'SUCCESSFULLY READ IN THE HEADER ... ON TO THE DATA' 

c this is the number of paths/mixed paths/layers that will be output each chunk
c iStore=(iNumGases+1) as 6 gases .. multiply by iNumLayers, and add 4 
c because we also do d/d(surface Temp), d/d(surfaceemis),
c d(thermal)/d(surace emis)
      iExtra=4
      iStore=0
      DO ii=1,iNatm
        iStore=iStore+iaNumLayers(ii)
        END DO
      iStore=(iNumGases+2)*iStore + (iExtra*iNatm) 

      IF (iBinOrText .LT. 0) THEN  
        print *,'For each k-comp file used in the kcarta.x run ... ' 
        print *,'Total number of gases = ',iNumGases 
        print *,'Total number of layers to output = ',iStore 
        print *,'Total number of surface temp/emiss to output = ',4 
        print *,'GRAND total to be output = ',iStore 
        print *,'Enter which of the above you want to store : '  
        read *,iWhichStore  
        END IF  

      IF (iOutNumber .NE. iNatm) THEN 
        print *,'iOutNumber,iNatm = ',iOutNumber,iNatm 
        print *,'iOutNumber .NE. iNatm' 
        STOP 
        END IF 

      iJ=0  
      DO iI=1,iOutNumber  
        IF (iaNumOut(iI) .NE. iaNumLayers(iI)) THEN  
          print *, 'iaNumOut(iI) .NE. iaNumLayers(iI), iI = ',iI 
          STOP 
          END IF 
        END DO  
   
      iBlSz=kMaxPts
      iNBl=1

      iMaxSet=iSetHigh

      CLOSE(iIOUN)

      IF (iBinOrText .LT. 0) THEN 
        iStore=iWhichStore                 !only need to read upto here!! 
        END IF 
 
c set up the artificial structure to mimic that of readkcarta.f 
      iArtificial=0 
      DO iJ=1,iNatm 
        DO iI=1,iImportant    !these are number of gas jacobs for this atm 
          iArtificial=iArtificial+1 
          iaArtificial(iArtificial)=iaNumOut(iJ) 
          END DO 
        iArtificial=iArtificial+1               !tempr jacob 
        iaArtificial(iArtificial)=iaNumOut(iJ) 
        iArtificial=iArtificial+1               !wgt fcn 
        iaArtificial(iArtificial)=iaNumOut(iJ) 
        iArtificial=iArtificial+1               !surface jacob 
        iaArtificial(iArtificial)=4 
        END DO 

c now iArtificial  <---> iOutNumber  (cf readkcarta.f) 
c now iaArtificial  <---> iaNumOut 

      iLoop=0         
      DO iOutTypeLoop = 1,iArtificial 
        DO iLpLp=1,iaArtificial(iOutTypeLoop) 
          iLoop=iLoop+1 
 
          IF (iLoop .GT. iStore) THEN 
            print *,'no need to read in more data!!' 
            STOP 
            END IF 
 
          print *,'Data # of (total) iStore ',iLoop,iStore 
 
c reread the header over and over again!!
          OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

          CALL readmainjacobheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh) 
          CALL readjacobinfo(iIOUN,iNumGases,iNatm,iaNumLayers) 
          CALL readgasjacob(iIOUN,iNumGases,100) 
          read(iIOUN) iTotal,iOutnumber,iImportant  
          read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)  
 
          CALL readdata(iSetLow,iMaxSet,iBlSz,
     $           iIOUN,iStore,iNatm,iLoop,raEntire,raFreq,iFreqPts,
     $           iTotal,iArtificial,iaArtificial,iOutTypeLoop,iLpLp) 

c this is to test the reading functions!
c if iLoop==1 then print the header info at top of output daya file
          CALL output2(raFreq,iFreqPts,raEntire,iStore,
     $                  caInName,iLoop,iBinOrText,iWhichStore)

          CLOSE(iIOUN)
          END DO
        END DO
         
      END 

c************************************************************************
      SUBROUTINE output2(raFreq,iFreqPts,raEntire,  
     $       iStore,caInName,iLoop,iBinOrText,iWhichStore)  
  
      include 'convolve.param'  
  
      REAL raEntire(kMaxEntire),raFreq(kMaxEntire)  
      INTEGER iFreqPts,iStore,iLoop,iBinOrText  
      INTEGER iWhichStore  
      CHARACTER*80 caInName  
  
c local variables  
      CHARACTER*80 caOutName  
      INTEGER iI,iJ,iIOUN  
  
      iIOUN=11  
        
c get the name of output file  
      CALL GetOutName2(caInName,caOutName,iBinOrText)  
  
      IF (iBinOrText .EQ. 1) THEN            !binary format for rdairs.m  
c put in the header info ... iI tells how many spectra to expect  
c only need to save this header when we are outputting the VERY first spectra  
        IF (iLoop .EQ. 1) THEN  
          ii=iStore  
          OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='NEW',  
     $         FORM='UNFORMATTED')  
          write(iIOUN) 1,iFreqPts,iI  
          write(iIOUN) (raFreq(iI),iI=1,iFreqPts)  
          CLOSE(iIOUN)  
          END IF  
c now output!!!  
        OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='OLD',  
     $     FORM='UNFORMATTED',ACCESS='APPEND')  
          write(iIOUN) (raEntire(iJ),iJ=1,iFreqPts)  
        CLOSE(iIOUN)  
        END IF  
  
      IF (iBinOrText .EQ. -1) THEN     !column text format for rdairs.m  
c put in the header info ... iI tells how many spectra to expect  
c only need to save this header when we are outputting the VERY first spectra  
        IF (iLoop .EQ. iWhichStore) THEN  
          OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='NEW',  
     $         FORM='FORMATTED')  
          DO iJ=1,iFreqPts  
            write(iIOUN,5050) raFreq(iJ),raEntire(iJ)  
            END DO  
          CLOSE(iIOUN)  
          END IF  
        END IF  
 5050       FORMAT(f10.4,'  ',1pe12.5)  
  
      RETURN  
      END     
  
c************************************************************************  
c this takes the input file name, finds the root (before .dat or .datJAC)   
c and appends .datCON  
      SUBROUTINE GetOutName2(caInName,caOutName,iBinOrText)  
  
      CHARACTER*80 caInName,caOutName  
      INTEGER iBinOrText  
        
      INTEGER iI  
  
      DO iI=1,80  
        caOutName(iI:iI)=' '  
        END DO  
  
c search for the .dat in caInName  
      iI=80  
 15         CONTINUE  
      IF ((caInName(iI:iI) .NE. '.') .AND. (iI .GT. 1)) THEN  
        iI=iI-1  
        GO TO 15  
        END IF  
      IF (iI .EQ. 1) THEN  
        print *,'wierd input name!!'  
        STOP  
        END IF  
  
      caOutName(1:iI-1)=caInName(1:iI-1)  
      IF (iBinOrText .GT. 0) THEN  
        caOutName(iI:iI+9)='JAC.datCON'  
      ELSE  
        caOutName(iI:iI+9)='JAC.datTXT'  
        END IF  
  
      RETURN  
      END  
c************************************************************************  
