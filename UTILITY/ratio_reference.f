      !this file loads in the reference profiles from different directories,
      !and spits out the ratio (to check klAYERS/kCARTA Src v1.09 vs 1.10)

      IMPLICIT NONE

      CHARACTER*80 caF1,caF2, caF
      REAL raP1(100),raP2(100),raPP1(100),raPP2(100)
      REAL raT1(100),raT2(100),raA1(100),raA2(100)
      INTEGER iProf

      caF1 = '../SRCv1.09/NewRefProfiles/refgas'  !infile1
      caF2 = '../SRCv1.10/NewRefProfiles/refgas'  !infile2
      caF  = '../WORK/ratio_refgas'               !results

      caF1 = '../DATA/RefProf.For.v107up/refgas'  !infile1
      caF2 = '../SRCv1.10/NewRefProfiles/refgas'  !infile2
      caF  = '../WORK/ratio_refgas'               !results

      print *,'Enter which profile to process ...'
      read *,iProf

      CALL SetName(caF1,iProf)
      CALL SetName(caF2,iProf)
      CALL SetName(caF ,iProf)

      CALL ReadFile(caF1,raP1,raPP1,raT1,raA1)      
      CALL ReadFile(caF2,raP2,raPP2,raT2,raA2)      
      CALL OutFile(caF,raP1,raPP1,raT1,raA1,raP2,raPP2,raT2,raA2)

      END

c************************************************************************
      SUBROUTINE SetName(caF,iI)
 
      IMPLICIT NONE
  
      CHARACTER*80 caF
      INTEGER iI
 
      CHARACTER*2 caString2, caTemp2
      INTEGER iInt,iLen

      iInt = 80 
 11   CONTINUE 
      IF ((caF(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
        END IF 
      iLen = iInt 
 
      DO iInt=1,2 
        caString2(iInt:iInt)=' ' 
        caTemp2(iInt:iInt)=' ' 
        END DO 

      WRITE(caString2,12) iI
 12   FORMAT(I2) 
c this is right justified ... change to left justified 
      iInt=1 
 14   CONTINUE 
      IF (caString2(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 14 
        END IF 
      caTemp2(1:2-iInt+1)=caString2(iInt:2) 

      caF(iLen+1:iLen+2) = caTemp2(1:2)
      print *,caF

      RETURN
      END
    
c************************************************************************
      SUBROUTINE ReadFile(caF,raP,raPP,raT,raA)

      IMPLICIT NONE

      CHARACTER*80 caF
      REAL raP(100),raPP(100),raT(100),raA(100)

      CHARACTER*80 caLine
      INTEGER iI,iJ,iIO

      iIO = 10
      OPEN(unit = iIO, file = caF, form = 'formatted')

 10   CONTINUE
      read(iIO,5030) caLine
      IF (caLine(1:1) .EQ. '!') THEN
        GOTO 10
        END IF

      iI = 1
      read (caLine,*) iJ,raP(iI),raPP(iI),raT(iI),raA(iI)
      DO iI = 2,100
        read (iIO,*) iJ,raP(iI),raPP(iI),raT(iI),raA(iI)
        END DO
      CLOSE(iIO)

 5030 FORMAT(A80)

      RETURN
      END

c************************************************************************
      SUBROUTINE OutFile(caF,raP1,raPP1,raT1,raA1,raP2,raPP2,raT2,raA2)

      IMPLICIT NONE

      CHARACTER*80 caF
      REAL raP1(100),raPP1(100),raT1(100),raA1(100)
      REAL raP2(100),raPP2(100),raT2(100),raA2(100)

      INTEGER iI,iIO

      iIO = 10
      OPEN(unit = iIO, file = caF, form = 'formatted')

      DO iI = 1,100
        write (iIO,30) iI,raP1(iI),raPP1(iI),raT1(iI),raA1(iI),
     $                    raP1(iI)/raP2(iI),raPP1(iI)/raPP2(iI),
     $                    raT1(iI)/raT2(iI),raA1(iI)/raA2(iI)
        END DO
      CLOSE(iIO)

 30   FORMAT(I3,' ',8(E12.6,'  '))

      RETURN
      END

c************************************************************************
      
