! this program takes in UNFORMATTED jacob output from atmos.x and allows the
! user to read the output

! the main rudiments of the program are the SAME AS readkcarta.f

! NOTE : Absoft behaves very wierdly with fseek,ftell
! so for Linux machines, use g77 or PDF Fortran
! wherever you see g77 in this file and readbinary.f, comment out 
!      iDummy = fseek(iIOUN,idatasize,1)    !!!!if using PDF(Linux) or SGI
! and use the 
!      CALL fseek(iIOUN,idatasize,1)         !!!!if using g77

      IMPLICIT NONE

      include 'convolve.param'

! the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName
      INTEGER iIOUN
      INTEGER iNatm,iSave,iOk,iChoice
      REAL rFr1,rFr2

      INTEGER iExtra,iNumGases,iaNumLayers(kMaxAtm),iNumLayers

      INTEGER iInstr,rInPointSp,rInStart,rFWHM,iBinOrText
      INTEGER iLorS,iNOAA,iInstrType

      INTEGER iaArtificial(kMaxLayer),iArtificial
      INTEGER iOutNumber,iaNumOut(kMaxPrint),iDummy,iDatapoints
      INTEGER FTELL
      INTEGER iTotal,iStart,iStore,iFreqPts,iImportant
      REAL rTemp,raFreq(kMaxEntire) 
      REAL raEntire(kMaxEntire),raConvolve(kMaxEntire) 
      INTEGER iI,iJ,iK,iLoop,iStartOffset,iJump,iDataSize,iWhichStore 
    
 90   CONTINUE
      iInstrType = 0

      print *,'Do you want a text or binary output? (-1/1) ' 
      read *,iBinOrText 

      print *,'Enter INPUT binary file name '
      read 1000,caInName
 1000 FORMAT(A80)

! first read header to figure out iStore ******************
      iIOUN=10
      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

      CALL readmainjacobheader(iIOUN,rFr1,rFr2,iNumLayers)

      CALL readjacobinfo(iIOUN,iNumGases,iNatm,iaNumLayers)

      CALL readgasjacob(iIOUN,iNumGases,iNumLayers)

      print *,'SUCCESSFULLY READ IN THE HEADER ... ON TO THE DATA' 

! this is the number of paths/mixed paths/layers that will be output each chunk

! this is the number of jacobian stuff that will be output each chunk
! iStore=(iNumGases+2) as d/dq for each gas, plus temp jac, plus wgt fcn ... 
! multiply by iNumLayers, and add 3  because we also do d/d(surface Temp), 
! d/d(surfaceemis), d(thermal)/d(surace emis)
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
      
! the following info is for each k-comp file 

      read(iIOUN) iTotal,iOutNumber,iImportant 
      read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber) 
      iStartOffset=ftell(iIOUN)  !figure out EXACTLY where data starts 

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
 
! figure out the frequencies 
      CALL findfreqsjacob(iIOUN,iFreqPts,iDataSize,iImportant,iaNumOut,iOutNumber,iTotal,raFreq) 

! set up the artificial structure to mimic that of readkcarta.f
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
 
! now  iArtificial  <---> iOutNumber
! now iaArtificial  <---> iaNumOut

      iLoop=0 
      DO iJ=1,iArtificial
 
! jump to point in file where main data blocks start 
        CALL fseek(iIOUN,iStartOffSet,0)         !!!!if using g77
 
        iJump=0 
        DO iI=1,iJ-1 
! jump to point in file where data for output option iJ starts  
!        READ(iIOUN) iMainType,iSubMainType,iNumberOut   
!        READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta   
          iJump=iJump+4+(3*4)+4           
          iJump=iJump+(4+(4*4)+4)        !this reads the above 2 header lines 
          !this reads all the data in iaArtificial(iI) 
          iJump=iJump+(4+(kMaxPts*4)+4)*iaArtificial(iI) 
          END DO 
        CALL fseek(iIOUN,iJump,1)         !!!!if using g77

! now skip over frequency header 
!      READ(iIOUN) iMainType,iSubMainType,iNumberOut   
!      READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta   
        iJump=4+(3*4)+4           
        iJump=iJump+(4+(4*4)+4)        !this reads the above 2 header lines 
        CALL fseek(iIOUN,iJump,1)         !!!!if using g77

        iStart=ftell(iIOUN) 
 
        print *,'outoption  ',iJ, ' of ',iArtificial
 
        DO iK = 1,iaArtificial(iJ) 
          CALL fseek(iIOUN,iStart,0)     !!!!if using g77
          iJump=(4+(kMaxPts*4)+4)*(iK-1) 
          CALL fseek(iIOUN,iJump,1)     !!!!if using g77
 
          CALL ReadData(iIOUN,iDataSize,raEntire,iDataPoints,iTotal) 
 
          IF (iDataPoints .NE. iFreqPts) THEN 
            print *,'iDataPoints .NE. iFreqPts',iDataPoints,iFreqPts 
            STOP 
            END IF 
              
! this is to test the reading functions! 
! if iLoop==1 then print the header info at top of output data file 
          iLoop=iLoop+1 
          CALL output2(raFreq,iFreqPts,raEntire,iStore,caInName,iLoop,iBinOrText,iWhichStore) 
 
          END DO 
        END DO 
      CLOSE(iIOUN) 

      END 

!************************************************************************
! this subroutine reads the freqs 
      SUBROUTINE findfreqsjacob(iIOUN,iFreqPts,iDataSize,iImportant,iaNumOut,iOutNumber,iTotal,raFreq) 
 
      include 'convolve.param' 
 
      REAL raFreq(kMaxEntire) 
      INTEGER iTotal,iIOUN,iOutNumber,iDataSize,iaNumOut(kMaxPrint) 
      INTEGER iFreqPts,iImportant
 
      INTEGER iI,iJ,iL,iMainType,iSubMainType,ikMaxPts,iNumberOut 
      INTEGER offset,iK,iF
      INTEGER FTELL
      INTEGER iDummy,iWhere 
      REAL rFrLow,rFrHigh,rDelta  
 
      iF=0 
      DO iI=1,iTotal 
        iWhere=FTELL(iIOUN) 
        DO iJ=1,iOutnumber          !loop over atmospheres 

          DO iL=1,iImportant+2      !iImportant Gases,one temp, one weight fcn
            READ(iIOUN) iMainType,iSubMainType,iNumberOut 
            READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta 
            IF (iNumberOut .NE. iaNumOut(iJ)) THEN
              print *,'iNumberOut .NE. iaNumOut(iJ), iJ = ',iJ
              STOP
              END IF
            offset=(ikMaxPts*4 + 2*4)*iaNumOut(iJ) 
            CALL fseek(iIOUN,offset,1)         !!!!if using g77
            END DO
          !surface jacobians
          READ(iIOUN) iMainType,iSubMainType,iNumberOut 
          READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta 
          IF (iNumberOut .NE. 4) THEN
            print *,'iNumberOut .NE. 4, iJ = ',iJ
            STOP
            END IF
          offset=(ikMaxPts*4 + 2*4)*4 
          CALL fseek(iIOUN,offset,1)         !!!!if using g77

          END DO 

        iDataSize=FTELL(iIOUN) 
        !set up raFreq 
        print *,rFrLow,rFrHigh,rDelta 
        if (abs(rDelta-0.0025) .le. 0.00001) then
          rDelta=0.0025
          end if       
        DO iL=1,ikMaxPts 
          iF=iF+1 
          raFreq(iF)=rFrLow+(iL-1)*rDelta 
          END DO 
        END DO 
 
      iFreqPts=iF 
 
      iDataSize=iDataSize-iWhere 
      iDataSize=iDataSize-(ikMaxPts*4 + 2*4)    !this is one data block) 

      RETURN 
      END 
 
!************************************************************************ 
      SUBROUTINE output2(raFreq,iFreqPts,raEntire,iStore,caInName,iLoop,iBinOrText,iWhichStore) 
 
      include 'convolve.param' 
 
      REAL raEntire(kMaxEntire),raFreq(kMaxEntire) 
      INTEGER iFreqPts,iStore,iLoop,iBinOrText 
      INTEGER iWhichStore 
      CHARACTER*80 caInName 
 
! local variables 
      CHARACTER*80 caOutName 
      INTEGER iI,iJ,iIOUN 
 
      iIOUN=11 
       
! get the name of output file 
      CALL GetOutName2(caInName,caOutName,iBinOrText) 
 
      IF (iBinOrText .EQ. 1) THEN            !binary format for rdairs.m 
! put in the header info ... iI tells how many spectra to expect 
! only need to save this header when we are outputting the VERY first spectra 
        IF (iLoop .EQ. 1) THEN 
          ii=iStore 
          OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='NEW',FORM='UNFORMATTED') 
          write(iIOUN) 1,iFreqPts,iI 
          write(iIOUN) (raFreq(iI),iI=1,iFreqPts) 
          CLOSE(iIOUN) 
          END IF 
! now output!!! 
        OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='OLD',FORM='UNFORMATTED',ACCESS='APPEND') 
          write(iIOUN) (raEntire(iJ),iJ=1,iFreqPts) 
        CLOSE(iIOUN) 
        END IF 
 
      IF (iBinOrText .EQ. -1) THEN     !column text format for rdairs.m 
! put in the header info ... iI tells how many spectra to expect 
! only need to save this header when we are outputting the VERY first spectra 
        IF (iLoop .EQ. iWhichStore) THEN 
          OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='NEW',FORM='FORMATTED') 
          DO iJ=1,iFreqPts 
            write(iIOUN,5050) raFreq(iJ),raEntire(iJ) 
            END DO 
          CLOSE(iIOUN) 
          END IF 
        END IF 
 5050    FORMAT(f10.4,'  ',1pe12.5) 
 
      RETURN 
      END    
 
!************************************************************************ 
! this takes the input file name, finds the root (before .dat or .datJAC) 
! and appends .datCON
      SUBROUTINE GetOutName(caInName,caOutName)

      CHARACTER*80 caInName,caOutName
      
      INTEGER iI

      DO iI=1,80
        caOutName(iI:iI)=' '
        END DO

! search for the .dat in caInName
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
      caOutName(iI+1:iI+6)='datCON'

      RETURN
      END
!************************************************************************
! this takes the input file name, finds the root (before .dat or .datJAC)  
! and appends .datCON 
      SUBROUTINE GetOutName2(caInName,caOutName,iBinOrText) 
 
      CHARACTER*80 caInName,caOutName 
      INTEGER iBinOrText 
       
      INTEGER iI 
 
      DO iI=1,80 
        caOutName(iI:iI)=' ' 
        END DO 
 
! search for the .dat in caInName 
      iI=80 
 15      CONTINUE 
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
!************************************************************************ 







