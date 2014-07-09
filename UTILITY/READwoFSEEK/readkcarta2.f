c program to read in UNFORMATTED output from kcarta.x and save to 
c easier-to-read UNFORMATTED file, or text file

      IMPLICIT NONE

      include 'convolve.param'

c the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName
      INTEGER iIOUN,iNumLayers,iSetLow,iSetHigh,iNumPaths
      INTEGER iNumPathsOut,iaOutPaths(kPathsOut)
      INTEGER iNpmix,iMixFileLines,iNumMixPathsOut,iNatm
      INTEGER iaOutMixPaths(kPathsOut),iaAtmNum(kMaxAtm)
      INTEGER iaNumLayersInAtm(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iaaOutRadPaths(kMaxAtm,kPathsOut)
      INTEGER iaNumLayersOut(kMaxAtm),ii1,ii2
      INTEGER iStore,ii,iaEms(kMaxAtm),iLpLp
      INTEGER iOk,iBlSz,iNBl,iMaxSet,iFreqPts
      INTEGER posn,offset,SizePerDataBlock,status,iLoop,iOutTypeLoop
      REAL raEntire(kMaxEntire)
      REAL rFr1,rFr2,raFreq(kMaxEntire)
      INTEGER iLorS,iNumRadiances

      INTEGER iInstrType,iInstr,rInPointSp,rInStart,rFWHM
      INTEGER iBinOrText,iWhichStore

      INTEGER iTotal,iOutNumber,iaNumOut(kMaxPrint)

 90   CONTINUE
      iInstrtype=0

      print *,'Do you want a text or binary output? (-1/1) '
      read *,iBinOrText

      print *,'Enter INPUT binary file name '
      read 1000,caInName
 1000 FORMAT(A80)

c first read header to figure out iStore ******************
      iIOUN=10
      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

      CALL readmainheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh,iLorS) 
      CALL readgaspaths(iIOUN,iNumPathsOut,iLorS) 
      CALL readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS) 
      CALL readatmospheres(iIOUN,iNpmix,iNatm,iaNumLayersInAtm, 
     $                     iaNumLayersOut) 
c this is the summary info that has to be read in each time!
      read(iIOUN) iTotal,iOutNumber
      read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)
 
      IF (iNpmix .EQ. 0) THEN
        iNumMixPathsOut=0
        END IF

      IF (iNatm .EQ. 0) THEN
        DO ii=1,iNatm
          iaNumLayersOut(ii)=0
          END DO
        END IF

      IF (iNatm .GT.  0) THEN
        DO ii=1,iNatm
          IF (iaNumLayersOut(ii) .EQ. -1) THEN
            iaNumLayersOut(ii)=iaNumLayersInAtm(ii)
            END IF
          END DO
        END IF

c      print *,'SUCCESSFULLY READ IN THE HEADER ... ON TO THE DATA' 

c this is the number of paths/mixed paths/layers that will be output each chunk
      iStore=0
      DO ii=1,iNatm
        iStore=iStore+iaNumLayersOut(ii)
        END DO                                 !num of radiances to output

      iWhichStore=iNumPathsOut+iNumMixPathsOut !num of paths/MP to output
      iStore=iStore+iWhichStore                !total num to output

      iNumRadiances=iStore-iWhichStore

      IF (iBinOrText .LT. 0) THEN
        print *,'For each k-comp file used in the atmos.x run '
        print *,'Total # of paths to output = ',iNumPathsOut
        print *,'Total # of mixed paths to output = ',iNumMixPathsOut
        print *,'Total # of radiances to output = ',iStore-iWhichStore
        print *,'GRAND total to be output = ',iStore
        print *,'Enter which of the above you want to store : '
        read *,iWhichStore
        END IF

c the following info is for each k-comp file
      iBlSz=kMaxPts
      iNBl=1

      iMaxSet=iSetHigh

      CLOSE(iIOUN)

      IF (iBinOrText .LT. 0) THEN
        iStore=iWhichStore                 !only need to read upto here!!
        END IF

      iLoop=0        
      DO iOutTypeLoop = 1,iOutNumber
        DO iLpLp=1,iaNumOut(iOutTypeLoop)
          iLoop=iLoop+1

          IF (iLoop .GT. iStore) THEN
            print *,'no need to read in more data!!'
            STOP
            END IF

          print *,'Data # of (total) iStore ',iLoop,iStore

c reread the header over and over again!!
          OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

          CALL readmainheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh,iLorS) 
          CALL readgaspaths(iIOUN,iNumPathsOut,iLorS) 
          CALL readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS) 
          CALL readatmospheres(iIOUN,iNpmix,iNatm,iaNumLayersInAtm, 
     $                     iaNumLayersOut) 
c this is the summary info that has to be read in each time!
          read(iIOUN) iTotal,iOutNumber
          read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)

          CALL readdata(iSetLow,iMaxSet,iBlSz,
     $           iIOUN,iStore,iNumPathsOut,iNumMixPathsOut,
     $           iaNumLayersOut,iNatm,iLoop,raEntire,raFreq,iFreqPts,
     $           iTotal,iOutNumber,iaNumOut,iOutTypeLoop,iLpLp)

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

      iIOUN=10
      
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
 5050 FORMAT(f10.4,'  ',1pe12.5)

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
 15   CONTINUE
      IF ((caInName(iI:iI) .NE. '.') .AND. (iI .GT. 1)) THEN
        iI=iI-1
        GO TO 15
        END IF
      IF (iI .EQ. 1) THEN
        print *,'wierd input name!!'
        STOP
        END IF

      caOutName(1:iI)=caInName(1:iI)
      IF (iBinOrText .GT. 0) THEN
        caOutName(iI+1:iI+6)='datCON'
      ELSE
        caOutName(iI+1:iI+6)='datTXT'
        END IF

      RETURN
      END
c************************************************************************
