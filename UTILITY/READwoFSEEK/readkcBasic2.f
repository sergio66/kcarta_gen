c this program takes in the UNFORMATTED output from kcarta.x and allows the
c user to convolve the final radiance and save to a FORMATTED file
c usage #1 : readkcBasic.x infile.dat outfile.bin
c       #2 : readkcBasic.x infile.dat outfile.txt iWhich
c usage    : readkcBasic.x arg1 arg2 [optional integer]
c where aerg1,arg2 are names of input and output files; if third argument is
c specified, the output file will be a text file that stores the iWhich th row
c of data; else the output file is a binary file

c NON fSEEK version

      IMPLICIT NONE

      include 'convolve.param'

c the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName,caOutName
      INTEGER iIOUN,iE,iB,iTotalStuff
      REAL rFr1,rFr2,rDv

      INTEGER iI,iJ,iK,iIOUN2,iTextOrBinary,iWhich,iChunks
      REAL raFreq(kMaxEntire),rF
      REAL raData(kMaxPts),raEntire(kMaxEntire)

      print *,'Enter INPUT binary file name '
      read 1000,caInName
 

      iIOUN = 10
      iIOUN2 = 11

      CALL DoCommandLine(caInName,caOutName,iTextOrBinary,iWhich)

      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')
      CALL readmainheaderShort(iIOUN,rFr1,rFr2,rDv,iB,iE,iTotalStuff)

      IF (iTextOrBinary .EQ. -1) THEN
        IF (iWhich .GT. iTotalStuff) THEN
          print *,'You want output # ',iWhich,' from file which has ',
     $      iTotalStuff,' outputs!!! Impossible!!!'
          STOP
          END IF
        OPEN(UNIT=iIOUN2,FILE=caOutName,STATUS='NEW',FORM='FORMATTED')
      ELSE
        OPEN(UNIT=iIOUN2,FILE=caOutName,STATUS='NEW',FORM='UNFORMATTED')
        WRITE (iIOUN2) 1,(iE-iB+1)*kMaxPts,iTotalStuff
        END IF

      print *,'reading in data .. might take a while .... '

      iChunks = iE - iB + 1
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (iTextOrBinary .NE. -1) THEN       !!dump all binary output
        DO iK = 1,iChunks*kMaxPts
          raFreq(iK) = rFr1 + (ik-1)*rDv
          END DO
        write(iIOUN2) (raFreq(iK),iK=1,iChunks*kMaxPts)

        DO iJ = 1,iTotalStuff
          CALL ReWinder(iIOUN,iJ)
          print *,iJ,iTotalStuff
          CALL Reader(iIOUN,raEntire,iJ,iTotalStuff,iChunks)
          write(iIOUN2) (raEntire(iK),iK=1,iChunks*kMaxPts)
          END DO
        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      CLOSE(iIOUN)
      CLOSE(iIOUN2)

 1000 FORMAT(A80)

      END
c************************************************************************
c this subroutine simply rewinds to beginning of file
c and then repositions to the data we want to start reading
      SubRoutine ReWinder(iIOUN,iJ)
 
      include 'convolve.param'

      INTEGER iIOUN,iJ

c local variables
      INTEGER iE,iB,iTotalStuff,iI,ik
      REAL rFr1,rFr2,rDv,raData(kMaxPts)

      REWIND(UNIT = iIOUN)
      CALL readmainheaderShort(iIOUN,rFr1,rFr2,rDv,iB,iE,iTotalStuff)

      DO iI = 1,iJ-1
         read(iIOUN) (raData(iK),iK=1,kMaxPts)
         END DO

      RETURN
      END
c************************************************************************
c this subroutine reads in the data we need!!!!!!! 
c more easily done using fseeks
       SUBROUTINE Reader(iIOUN,raEntire,iJ,iTotalStuff,iChunks)

      include 'convolve.param'

      INTEGER iIOUN,iJ,iChunks,iTotalStuff
      REAL raEntire(kMaxEntire)

c local variables
      REAL raData(kMaxPts)
      INTEGER iI,iK,iL

      DO iI = 1,iChunks
        !!data reader should already be pointing at relevant data set of the
        !!first KCARTA chunk
        read(iIOUN) (raData(iK),iK=1,kMaxPts)        
        DO iK = 1,kMaxPts
          raEntire((iI-1)*kMaxPts + iK) = raData(iK)
          END DO
        !!now jump to next data chunk, if it exists
        IF (iI .LT. iChunks) THEN
          DO iL = 1,(iTotalStuff-1)
            read(iIOUN) (raData(iK),iK=1,kMaxPts)        
            END DO
          END IF
        END DO

      RETURN
      END
c************************************************************************
      
