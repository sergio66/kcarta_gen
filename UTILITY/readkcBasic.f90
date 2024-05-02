! this program takes in the UNFORMATTED output from kcarta.x and allows the
! user to convolve the final radiance and save to a FORMATTED file

! usage #1 : readkcBasic.x infile.dat outfile.bin
!       #2 : readkcBasic.x infile.dat outfile.txt iWhich
! usage    : readkcBasic.x arg1 arg2 [optional integer]
! where arg1,arg2 are names of input and output files; if third argument is
! specified, the output file will be a text file that stores the iWhich th row
! of data; else the output file is a binary file

! NOTE : Absoft behaves very wierdly with fseek,ftell
! so for Linux machines, use g77 or PDF Fortran
! while for SGI, you can use the supplied software

! wherever you see g77 in this file and readbinary.f, comment out 
!      iDummy = fseek(iIOUN,idatasize,1)    !!!!if using PDF (Linux) or SGI
! and use the 
!      CALL fseek(iIOUN,idatasize,1)         !!!!if using g77

      IMPLICIT NONE

      include 'convolve.param'

! the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName,caOutName
      INTEGER iIOUN,iE,iB,iTotalStuff
      REAL rFr1,rFr2,rDv

      INTEGER iI,iJ,iK,iIOUN2,iTextOrBinary,iWhich,iChunks
      REAL raFreq(kMaxEntire),rF
      REAL raData(kMaxPts),raEntire(kMaxEntire)

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
!      CALL DoCommandLine(caInName,caOutName,iTextOrBinary,iWhich)

      print *,'Do you want a text or binary output? (-1/1) '
      read *,iTextOrBinary

!      print *,'Enter INPUT binary file name '
!      read 1000,caInName
! 1000 FORMAT(A80)
      caInName  = fin
      caOutName = fout
       
      iIOUN = 10
      iIOUN2 = 11

      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')
      CALL readmainheaderShort(iIOUN,rFr1,rFr2,rDv,iB,iE,iTotalStuff)

      IF (iTextOrBinary .EQ. -1) THEN
        IF (iWhich .GT. iTotalStuff) THEN
          write(*,'(A,I3,A,I3,A)') 'You want output # ',iWhich,' from file which has ',iTotalStuff,' outputs!!! Impossible!!!'
          STOP
          END IF
        OPEN(UNIT=iIOUN2,FILE=caOutName,STATUS='NEW',FORM='FORMATTED')
      ELSE
        OPEN(UNIT=iIOUN2,FILE=caOutName,STATUS='NEW',FORM='UNFORMATTED')
        WRITE (iIOUN2) 1,(iE-iB+1)*kMaxPts,iTotalStuff
        END IF

      print *,'reading in data .. might take a while .... '

      iChunks = iE - iB + 1

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      IF (iTextOrBinary .EQ. -1) THEN       !!dump selected text output
        DO iI = 1,iChunks
          DO iJ = 1,iTotalStuff
            read(iIOUN) (raData(iK),iK=1,kMaxPts)
            IF (iJ .EQ. iWhich) THEN 
              DO iK = 1,kMaxpts   
                rF = rFr1 + ((iI-1)*kMaxPts + iK-1)*rDv
                write(iIOUN2,*) rF,raData(iK)
                END DO
              END IF
            END DO
          END DO
        END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (iTextOrBinary .NE. -1) THEN       !!dump all binary output
        DO iK = 1,iChunks*kMaxPts
          raFreq(iK) = rFr1 + (ik-1)*rDv
          END DO
        write(iIOUN2) (raFreq(iK),iK=1,iChunks*kMaxPts)

        iChunks = iE - iB + 1
        DO iJ = 1,iTotalStuff
          CALL ReWinderSEEK(iIOUN,iJ)
          CALL ReaderSEEK(iIOUN,raEntire,iJ,iTotalStuff,iChunks)
          write(iIOUN2) (raEntire(iK),iK=1,iChunks*kMaxPts)
          END DO
        END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CLOSE(iIOUN)
      CLOSE(iIOUN2)

 1000 FORMAT(A80)


      END
!************************************************************************
! this subroutine simply rewinds to beginning of file
! and then repositions to the data we want to start reading
      SubRoutine ReWinderSEEK(iIOUN,iJ)

      IMPLICIT NONE
 
      include 'convolve.param'

      INTEGER iIOUN,iJ

! local variables
      INTEGER iE,iB,iTotalStuff,iI,iK,iDummy,iDataSize
      REAL rFr1,rFr2,rDv,raData(kMaxPts)

!     iDummy=ftell(iIOUN)

      REWIND(iIOUN)

      CALL readmainheaderShort(iIOUN,rFr1,rFr2,rDv,iB,iE,iTotalStuff)

      iDataSize = (kMaxPts*4 + 2*4)*(iJ-1)
      CALL fseek_local(iIOUN,idatasize,1)         !!!!if using g77

      RETURN
      END
!************************************************************************
! this subroutine reads in the data we need!!!!!!! 
! more easily done using fseeks
       SUBROUTINE ReaderSEEK(iIOUN,raEntire,iJ,iTotalStuff,iChunks)

      IMPLICIT NONE
      include 'convolve.param'

      INTEGER iIOUN,iJ,iChunks,iTotalStuff
      REAL raEntire(kMaxEntire)

! local variables
      REAL raData(kMaxPts)
      INTEGER iI,iK,iL,iDummy,iDataSize

      iDataSize = (kMaxPts*4 + 2*4)*(iTotalStuff-1)
      DO iI = 1,iChunks
        !!data reader should already be pointing at relevant data set of the
        !!first KCARTA chunk
        read(iIOUN) (raData(iK),iK=1,kMaxPts)        
        DO iK = 1,kMaxPts
          raEntire((iI-1)*kMaxPts + iK) = raData(iK)
          END DO
        !!now jump to next data chunk, if it exists
        IF ((iI .LT. iChunks) .AND. (iDataSize .GT. 0)) THEN
          CALL fseek_local(iIOUN,idatasize,1)          !!!!if using g77
          END IF
        END DO

      RETURN
      END
!***********************************************************************
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
