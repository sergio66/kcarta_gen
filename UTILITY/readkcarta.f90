! this program takes in the UNFORMATTED output from atmos.x and allows the
! user to read the final radiance and save to a FORMATTED or UNFORMATTED file

! either do Gaussian SRF,AIRS SRF or NOAA SRF -- 
! see paramSRF.f for AIRS details, paramNOAA for NOAA details

! NOTE : Absoft behaves very wierdly with fseek,ftell
! so for Linux machines, use g77 or PDF Fortran
! this file is for g77 (which uses subroutines fseek instead of functions)
!          iDummy=fseek(iIOUN,offset,1)       !! if using absoft
!          CALL fseek(iIOUN,OffSet,1)         !! if using g77
!
! example : ../BIN/readkcarta.x fin=../WORK/junk.dat
!
! https://community.intel.com/t5/Intel-Fortran-Compiler/Errors-using-GETPID/m-p/1094827?profile.language=es
! ifort -fpp -D IFORT ex_getpid.f90 &&  ./a.out
! my pid=       29421@
!
! gfortran -cpp ex_getpid.f90 && ./a.out
! my pid=       29430

      PROGRAM readkcarta

#ifdef IFORT
    USE IFPORT
#endif

      IMPLICIT NONE

      include 'convolve.param'

! the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName
      CHARACTER*130 caaMixedPathInfo(kMixFilRows)
      INTEGER iIOUN,iNumPathsOut
      INTEGER iNpmix,iNumMixPathsOut,iNatm
      INTEGER iaNumLayersOut(kMaxAtm),iaNumLayersInAtm(kMaxAtm)

      REAL rFr1,rFr2

      INTEGER iLorS,iNOAA,iInstrType
      INTEGER iInstr,rInPointSp,rInStart,rFWHM
      INTEGER iBinOrText

      INTEGER iOutNumber,iaNumOut(kMaxPrint),iDummy,iDatapoints
      INTEGER iTotal,iStart,iStore,iFreqPts
      REAL rTemp,raFreq(kMaxEntire)
      REAL raEntire(kMaxEntire),raConvolve(kMaxEntire)
      INTEGER iI,iJ,iK,iLoop,iStartOffset,iJump,iDataSize,iWhichStore,iFileWhere

!      INTEGER ftell

! for vararg
      CHARACTER*80 FIN, FOUT, mode
      INTEGER IARGC,IOERR
      INTEGER J
      INTEGER NARGS
 
      CHARACTER*80 BUF
      CHARACTER*80 VAL
      CHARACTER*80 VAR

!************************************************************************
! see /asl/packages/rtpV201/utils/rdinfo_subset.f
!      ------------
!      Set defaults
!      ------------
       IOERR = 6
       FIN =  'kcarta.dat'           ! input filename
       FOUT = 'kcarta.txt'           ! output filename

!      -----------------------------------------------------------------
!      Loop on program parameters
!      --------------------------
!      Determine the number of command-line arguments
       NARGS=IARGC()

       IF (NARGS .EQ. 0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('readkcarta.f90 must be run with at least 1 argument')
          WRITE(IOERR,1011)
 1011     FORMAT('   fin  = <filename> txt input file  {mandatory}','//'    &
                 '   fout = <filename> rtp output file {optional, set to kcarta.txt}')
          STOP
       ENDIF

!      Loop over the command-line arguments
       DO iI = 1, NARGS

!         Pull out the ith argument
          CALL GETARG(iI, BUF)

!         Find the "=" character in the command-line argument string
          J=INDEX(BUF, '=')

          IF (J .NE. 0) THEN

!            Name of variable
             VAR = BUF(1:J-1)
             CALL UPCASE(VAR)

!            Specified value
             VAL = BUF(J+1:LEN(BUF))

!            Big "IF" to set parameters
!            ----------------------------
             IF (VAR(1:3) .EQ. 'FIN') THEN
                FIN=VAL
             ELSEIF (VAR(1:4) .EQ. 'FOUT') THEN
                FOUT=VAL
             ELSE
                WRITE(IOERR,1020) VAR
 1020           FORMAT('Unknown command-line argument: ',A6)
                STOP

             ENDIF

          ENDIF
       ENDDO  ! end of loop over command-line arguments
       write(*,'(A,A)') 'input  file name = ',fin
       write(*,'(A,A)') 'output file name = ',fout

!************************************************************************

 90   CONTINUE
      iInstrtype=0

      print *,'Do you want a text or binary output? (-1/1) '
      read *,iBinOrText

!      print *,'Enter INPUT binary file name '
!      read 1000,caInName
! 1000 FORMAT(A80)
       caInName = fin

! first read header to figure out iStore ******************
      iIOUN=10
!      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED',
!     $     RECL=8)
      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

      iFileWhere = ftell(iIOUN)
      print *,'location 0 : just opened file iFileWhere = ',iFileWhere

      CALL readmainheaderLONG(iIOUN,rFr1,rFr2,iLorS)

      iFileWhere = ftell(iIOUN)
      print *,'location 1 : readmainheader -  iFileWhere = ',iFileWhere

      CALL readgaspaths(iIOUN,iNumPathsOut,iLorS)

      iFileWhere = ftell(iIOUN)
      print *,'location 2 : readgaspaths -  iFileWhere = ',iFileWhere

      CALL readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS)

      iFileWhere = ftell(iIOUN)
      print *,'location 3 : readmixedpaths -  iFileWhere = ',iFileWhere

      IF (iNpmix .EQ. 0) THEN
        iNumMixPathsOut=0
      END IF

      CALL readatmospheres(iIOUN,iNpmix,iNatm,iaNumLayersInAtm,iaNumLayersOut)

      iFileWhere = ftell(iIOUN)
      print *,'location 4 : readatmospheres -  iFileWhere = ',iFileWhere

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

! this is the number of paths/mixed paths/layers that will be output each chunk
      iStore=0
      DO ii=1,iNatm
        iStore=iStore+iaNumLayersOut(ii)
      END DO                                 !num of radiances to output

      iWhichStore=iNumPathsOut+iNumMixPathsOut !num of paths/MP to output
      iStore=iStore+iWhichStore                !total num to output

      IF (iBinOrText .LT. 0) THEN
        print *,'For each k-comp file used in the atmos.x run '
        print *,'Total # of paths to output       = ',iNumPathsOut
        print *,'Total # of mixed paths to output = ',iNumMixPathsOut
        print *,'Total # of radiances to output   = ',iStore-iWhichStore
        print *,'GRAND total to be output         = ',iStore
        print *,'Enter which of the above you want to store : '
        read *,iWhichStore
        print *,'  ... wiil read output NUMBER ',iWhichStore
      END IF

! the following info is for each k-comp file

! cccccccc this is the summary info that HAS to be read in!!!!!!!!
      read(iIOUN) iTotal,iOutNumber
      read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)
      iStartOffset = ftell(iIOUN)  !figure out EXACTLY where data starts
      iJ=0
      DO iI=1,iOutNumber
        iJ=iJ+iaNumOut(iI)
        END DO
      IF (iJ .NE. iStore) THEN
        print *,'hmm : number of things to read  aintcha jivin'
        STOP
      END IF

!      print *,'jhs',iStartOffset

      iFileWhere = ftell(iIOUN)
      print *,'location 5 : before call findfreqs iFileWhere = ',iFileWhere

! figure out the frequencies
      CALL findfreqs(iIOUN,iFreqPts,iDataSize,iaNumOut,iOutNumber,iTotal,raFreq)
      print *,'klm'

!      IF (abs(rFr1-raFreq(1)) .GT. 0.02) THEN
!        print *,'lower freq bounds make no sense'
!        STOP
!        END IF       
!      IF (abs(rFr2-raFreq(iFreqPts)) .GT. 0.02) THEN
!        print *,'upper freq bounds make no sense'
!        STOP
!        END IF       

      print *,'start/stop wavenumbers are : ',rFr1,rFr2

      iLoop=0
      DO iJ=1,iOutnumber

! jump to point in file where main data blocks start
        CALL fseek_local(iIOUN,iSTartOffSet,0)         !!!!if using g77

        iJump=0
        DO iI=1,iJ-1
! jump to point in file where data for output option iJ starts 
!c        READ(iIOUN) iMainType,iSubMainType,iNumberOut  
!cc        READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta  
          !this reads above 2 header lines
          iJump=iJump+kByteSize+(3*kRealSize)+kByteSize          
          iJump=iJump+(kByteSize+(4*kRealSize)+kByteSize) 
          !this reads all the data in iaNumOut(iI)
          iJump=iJump+(kByteSize+(kMaxPts*kRealSize)+kByteSize)*iaNumOut(iI)
          END DO

        CALL fseek_local(iIOUN,iJump,1)         !!!!if using g77

! now skip over frequency header
!c      READ(iIOUN) iMainType,iSubMainType,iNumberOut  
!c      READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta  
        !this reads above 2 header lines
        iJump=kByteSize+(3*kRealSize)+kByteSize          
        iJump=iJump+(kByteSize+(4*kRealSize)+kByteSize)        
        CALL fseek_local(iIOUN,iJump,1)         !!!!if using g77

        iStart = ftell(iIOUN)

        print *,'outoption  ',iJ, ' of ',iOutnumber

        DO iK = 1,iaNumOut(iJ)
          CALL fseek_local(iIOUN,iStart,0)         !!!!if using g77
          iJump=(kByteSize+(kMaxPts*kRealSize)+kByteSize)*(iK-1)
          CALL fseek_local(iIOUN,iJump,1)         !!!!if using g77

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

      END PROGRAM

!************************************************************************
! this subroutine reads the freqs
      SUBROUTINE findfreqs(iIOUN,iFreqPts,iDataSize,iaNumOut,iOutNumber,iTotal,raFreq)

#ifdef IFORT
    USE IFPORT
#endif

      implicit none

      include 'convolve.param'

      REAL raFreq(kMaxEntire)
      INTEGER iTotal,iIOUN,iOutNumber,iDataSize,iaNumOut(kMaxPrint)
      INTEGER iFreqPts

      INTEGER iI,iJ,iL,iMainType,iSubMainType,ikMaxPts,iNumberOut
      INTEGER offset,iK,iF,iDummy,iFileWhere
      REAL rFrLow,rFrHigh,rDelta 
!      INTEGER FTELL

!     print *,'aaaa ',0,0,iTotal

      iFileWhere = ftell(iIOUN)
      print *,'just entered  findfreqs iFileWhere = ',iFileWhere

      iF=0
      DO iI=1,iTotal
        iFileWhere = FTELL(iIOUN)
        write(*,'(A,3(I10))') '    findfreqs starting loop --- ii,iTotal,iFileWhere = ',iI,iTotal,iFileWhere
        DO iJ=1,iOutnumber
!          print *,'aaaa ',iI,iJ,iTotal,iFileWhere
          READ(iIOUN) iMainType,iSubMainType,iNumberOut
          READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta
          write(*,'(6X,3(I10))') iMainType,iSubMainType,iNumberOut
          write(*,'(6X,I10,3(F20.12))') ikMaxPts,rFrLow,rFrHigh,rDelta

          offset = (ikMaxPts*kByteSize + 2*kRealSize)*iaNumOut(iJ)
          iFileWhere = FTELL(iIOUN)
          print *,'    findfreqs starting inner loop --- iJ,iOutnumber,offset,iFileWhere = ',iJ,iOutNumber,offset,iFileWhere

!          CALL fseek_local(iIOUN,OffSet,1)             !!!!if using g77
          CALL fseek_local(iIOUN,iOutNumber,1)         !!!!if using g77
        END DO

        iDataSize = FTELL(iIOUN)
        !set up raFreq
        if (abs(rDelta-0.0025) .le. 0.00001) then
          rDelta=0.0025
        end if       
        DO iL=1,ikMaxPts
          iF=iF+1
          raFreq(iF)=rFrLow+(iL-1)*rDelta
        END DO
      END DO

      iFreqPts=iF

      iDataSize = iDataSize-iFileWhere
      iDataSize = iDataSize-(ikMaxPts*kRealSize + 2*kByteSize)  !one data block

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
 5050 FORMAT(f10.4,'  ',1pe12.5)

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

!************************************************************************
       SUBROUTINE UPCASE(BUF)

!      Convert a string to upper case

       CHARACTER*(*) BUF

       INTEGER I
       INTEGER IC
       INTEGER J

       DO I=1,LEN(BUF)
          IC=ICHAR( BUF(I:I) )
          IF (IC .GE. ICHAR('a') .AND. IC .LE. ICHAR('z')) THEN
             J=IC + (ICHAR('A') - ICHAR('a'))
             BUF(I:I)=CHAR(J)
          ENDIF
       ENDDO

       RETURN
       END

!************************************************************************
