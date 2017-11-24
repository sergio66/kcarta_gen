      !ifort -extend-source 132 f77_nchar_to_120char.f -o f77_nchar_to_120char


PROGRAM f77_nchar_to_120char  

IMPLICIT none  
INTEGER :: iIOUN1, iIOUN2, IERR, iCnt, iLong, iLen, xlen, iLenCut  
CHARACTER (LEN=80) :: caF1  
CHARACTER (LEN=81) :: caF2  

CHARACTER (LEN=132) :: caIn, caX  

INTEGER :: iI  
iLenCut = 100  

iLenCut = 80  
print * , ' input caF1 : '  
read *, caF1  

caF2 = 'x'//caF1  
iIOUN1 = 11  
iIOUN2 = 12  
OPEN (UNIT = iIOUN1, FILE = caF1, STATUS = 'OLD', FORM = &
 'FORMATTED', IOSTAT = IERR)

OPEN (UNIT = iIOUN2, FILE = caF2, STATUS = 'UNKNOWN', FORM = &
 'FORMATTED', IOSTAT = IERR)
iCnt = 0  
iLong = 0  
  101 CONTINUE  
read (iIOUN1, '(A)', END = 199) caIn  
iCnt = iCnt + 1  
caIn = trim (caIn)  
ilen = xlen (caIn)  
IF (iLen.GT.iLenCut) THEN  
!        write(*,'(I6,I6,A)') iCnt,iLen,caIn
  iLong = iLong + 1  
  CALL CutStr (iLenCut, iCnt, caIn, caX)  
  write (iIOUN2, '(A)') caIn (1:iLenCut)  
  write (iIOUN2, '(A)') caX (1:iLenCut)  
ELSE  
        !!! no problemo
  write (iIOUN2, '(A)') caIn (1:iLenCut)  

ENDIF  

GOTO 101  
  199 CONTINUE  
CLOSE (iIOUN1)  

CLOSE (iIOUN2)  
write ( * , '(A,A,A,I4,A,I4,A,I3,A)') 'file ', trim (caF1) , ' had &
& ', iCnt, linesofwhich',iLong,'werelongerthan',iLenCut,'char'


end PROGRAM f77_nchar_to_120char
!************************************************************************

INTEGER FUNCTION xlen (caStrx)  

IMPLICIT NONE  
CHARACTER (LEN=132) :: caStrx  

INTEGER :: iI, iFound  
iFound = - 1  
iI = len (caStrx)  
DO WHILE ( (iI.GT.1) .AND. (iFound.LT.0) )  
IF (caStrx (iI:iI) .NE.' ') THEN  
  iFound = iI  
ELSE  
  iI = iI - 1  
ENDIF  

ENDDO  

xlen = iI  
RETURN  


END function xlen
!************************************************************************

SUBROUTINE CutStr (iLenCut, iCnt, ca1, ca2)  

IMPLICIT NONE  
INTEGER :: iLenCut, iCnt  

CHARACTER (LEN=132) :: ca1, ca2  

CHARACTER (LEN=132) :: ca0  

INTEGER :: iI, iLen, xlen, iFound  
ca0 = ca1  

iLen = xlen (ca0)  
DO iI = 1, 132  
ca2 (iI:iI) = ' '  

ENDDO  
iFound = - 1  
iI = iLenCut - 1  
DO WHILE ( (iI.GT.1) .AND. (iFound.LT.0) )  
IF ( (ca1 (iI:iI) .EQ.' ') .OR. (ca1 (iI:iI) .EQ.',') ) THEN  
  iFound = iI  
ELSE  
  iI = iI - 1  
ENDIF  

ENDDO  
ca2 (1:7)  = '     c '  

ca2 (8:8 + iLen - iFound) = ca1 (iFound:iLen)  
DO iI = iFound+1, 132  
ca1 (iI:iI) = ' '  

ENDDO  
print * , ' '  
write ( * , '(A,I4,A,A,A)') 'Line ', iCnt, ' Split the input strin &
&g <', , ca0 (1:iLen) , '> into the &
&           following'
write ( * , '(A,A,A)') ca1, ' and ', ca2  

print * , ' '  
RETURN  

END SUBROUTINE CutStr
