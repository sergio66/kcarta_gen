! Copyright 2000
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:46
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!******** THIS FILE CONTAINS VARIOUS USEFUL SUBROUTINES/FUNCTIONS *******
!******** such as linear interp and sorting routines                 ****
!************************************************************************
! simple function to see if integer iI is a member of a list iaSet (that
! has iElements in it)
! output : index into set if in set, -1 if not

INTEGER FUNCTION  InSet(iI,iaSet,iNumElements)


INTEGER, INTENT(IN OUT)                  :: iI
INTEGER, INTENT(IN OUT)                  :: iaSet(*)
NO TYPE, INTENT(IN OUT)                  :: iNumElemen
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INTEGER :: iNumElements

INTEGER :: iJ,iAns

! assume not in set
iAns = -1

!      IF ((iNumElements .LT. 1).OR.(iNumElements .GT. kMaxUserSet)) THEN
IF (iNumElements < 1) THEN
  WRITE(kStdErr,*) 'fcn INSET : need valid number of elements',iNumElements
  CALL DoSTOP
END IF

iJ = 0
10   CONTINUE
iJ = iJ+1
IF (iaSet(iJ) == iI) THEN
  iAns = iJ
ELSE IF (iJ < iNumElements) THEN
  GO TO 10
END IF

InSet = iAns
RETURN
END FUNCTION  InSet

!************************************************************************
! this function does the integer DIVISION eg 100 div 3 = 33, 10 div 2 = 5
! assuming i1,i2 > 0

INTEGER FUNCTION idiv(i1,i2)


INTEGER, INTENT(IN OUT)                  :: i1
INTEGER, INTENT(IN OUT)                  :: i2
IMPLICIT NONE

INTEGER :: iInt

! recall INT truncates the real number
iInt=INT((i1*1.0)/(i2*1.0))
idiv=iInt
RETURN
END FUNCTION idiv

!************************************************************************
! this function does the integer MOD eg 100 div 3 = 1, 10 div 2 = 0
! assuming i1,i2 > 0  it does mod(i1,i2) = i1 - int(i1/i2)*i2
! so eg mod(100,3) = 100 - int(100/3)*3 = 100 - 33*3 = 1
! so eg mod(100,5) = 100 - int(100/5)*5 = 100 - 20*5 = 0
! so eg mod(5,100) = 5 - int(5/100)*100 = 5 - 0*100 = 5
! so eg mod(10,10) = 10 - int(10/10)*10 = 10 - 1*10 = 0

INTEGER FUNCTION iimod(i1,i2)


INTEGER, INTENT(IN)                      :: i1
INTEGER, INTENT(IN)                      :: i2
IMPLICIT NONE

INTEGER :: iInt

! recall INT truncates the real number
iInt = INT((i1*1.0)/(i2*1.0))
iimod = i1 - iInt*i2

RETURN
END FUNCTION iimod

!************************************************************************
! this integer function is "floor" -- assume rX > 0

INTEGER FUNCTION iFloor(rX)


REAL, INTENT(IN OUT)                     :: rX
IMPLICIT NONE



INTEGER :: iI,iIm1,iIp1,iIm2,iIp2,IF

IF=nint(rX)

iI=nint(rX)
iIm1=iI-1
iIm2=iI-2
iIp1=iI+1
iIp2=iI+2

IF (rX >= iIm2*1.0) THEN
  IF=iIm2
END IF
IF (rX >= iIm1*1.0) THEN
  IF=iIm1
END IF
IF (rX >= iI*1.0) THEN
  IF=iI
END IF
IF (rX >= iIp1*1.0) THEN
  IF=iIp1
END IF
IF (rX >= iIp2*1.0) THEN
  IF=iIp2
END IF

iFloor = IF

RETURN
END FUNCTION iFloor
!************************************************************************
! this integer function is "ceil" -- assume rX > 0

INTEGER FUNCTION iCeil(rX)


REAL, INTENT(IN OUT)                     :: rX
IMPLICIT NONE



INTEGER :: iI,iIm1,iIp1,iIm2,iIp2,iC

iC=nint(rX)

iI=nint(rX)
iIm1=iI-1
iIm2=iI-2
iIp1=iI+1
iIp2=iI+2

IF (rX <= iIp2*1.0) THEN
  iC=iIp2
END IF
IF (rX <= iIp1*1.0) THEN
  iC=iIp1
END IF
IF (rX <= iI*1.0) THEN
  iC=iI
END IF
IF (rX <= iIm1*1.0) THEN
  iC=iIm1
END IF
IF (rX <= iIm2*1.0) THEN
  iC=iIm2
END IF

iCeil = iC

RETURN
END FUNCTION iCeil

!************************************************************************
!                  SEARCHING and SORTING routines
!************************************************************************
! this subroutine sorts the arr values and returns indx values
! there are No <= kMaxPts values to be sorted

SUBROUTINE NumericalRecipesIndexer(indx,arr,No)


INTEGER, INTENT(OUT)                     :: indx(*)
REAL, INTENT(IN)                         :: arr(*)
INTEGER, INTENT(IN)                      :: No
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'


INTEGER :: !integer array that gives the indices
REAL :: !array of abs coeffs of layer closest to gnd


INTEGER, PARAMETER :: M=7
INTEGER, PARAMETER :: NSTACK=50
INTEGER :: n,i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
REAL :: a

! sort the k values, saving the index
! refer Numerical Recipes in F77, Second Edition, pg 330 (subroutine indexx)

IF (No > kMaxPts) THEN
  WRITE (kStdErr,*) 'can only sort and index <= kMaxPts points'
  CALL DoStop
END IF

n=No
DO j=1,n
  indx(j) = j
END DO

jstack=0
l=1
ir=n
1    IF ((ir-l) < M) THEN
  DO j=l+1,ir
    indxt=indx(j)
    a=arr(indxt)
    DO i=j-1,1,-1
      IF (arr(indx(i)) <= a) GO TO 2
      indx(i+1)=indx(i)
    END DO
    i=0
    2        indx(i+1)=indxt
  END DO
  IF (jstack == 0) GO TO 1000
  ir = istack(jstack)
  l=istack(jstack-1)
  jstack=jstack-2
ELSE
  k=(l+ir)/2
  itemp=indx(k)
  indx(k)=indx(l+1)
  indx(l+1)=itemp
  IF (arr(indx(l+1)) > arr(indx(ir))) THEN
    itemp=indx(l+1)
    indx(l+1)=indx(ir)
    indx(ir)=itemp
  END IF
  IF (arr(indx(l)) > arr(indx(ir))) THEN
    itemp=indx(l)
    indx(l)=indx(ir)
    indx(ir)=itemp
  END IF
  IF (arr(indx(l+1)) > arr(indx(l))) THEN
    itemp=indx(l+1)
    indx(l+1)=indx(l)
    indx(l)=itemp
  END IF
  i=l+1
  j=ir
  indxt=indx(l)
  a=arr(indxt)
  3      CONTINUE
  i=i+1
  IF (arr(indx(i)) < a) GO TO 3
  4      CONTINUE
  j=j-1
  IF (arr(indx(j)) > a) GO TO 4
  IF (j < i) GO TO 5
  itemp = indx(i)
  indx(i) = indx(j)
  indx(j) = itemp
  GO TO 3
  5      indx(l)=indx(j)
  indx(j)=indxt
  jstack=jstack+2
  IF (jstack > nstack) THEN
    WRITE(kStdErr,*) 'nstack too small in indexx (sorting k)'
    CALL DoStop
  END IF
  IF (ir-i+1 >= j-l) THEN
    istack(jstack) = ir
    istack(jstack-1) = i
    ir=j-1
  ELSE
    istack(jstack) = j-1
    istack(jstack-1) = l
    l=i
  END IF
END IF
GO TO 1
1000 CONTINUE

RETURN
END SUBROUTINE NumericalRecipesIndexer

!************************************************************************
! this subroutine sorts an integer array
! this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

SUBROUTINE DoSort(iaArr,iCnt)


INTEGER, INTENT(IN OUT)                  :: iaArr(*)
INTEGER, INTENT(IN)                      :: iCnt
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iaArr = integer array to be sorted
! iCnt  = sort indices 1..iCnt of iaArr


! local variables
INTEGER :: iTrue,iT,i1,i2

i1=1
iTrue=1

60   CONTINUE
IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
  i2=iCnt
  iTrue=-1
  70     CONTINUE
  IF (i2 > i1) THEN
    IF (iaArr(i2) < iaArr(i2-1)) THEN
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
END SUBROUTINE DoSort

!************************************************************************
! this subroutine sorts a real array raP1 into ascending order or
! descending order depending on value of iUD ... ia1,ra1 are accordingly set
! NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!
! this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

SUBROUTINE DoSortPress(ia1,ra1,raP1,iCnt,iUD)


INTEGER, INTENT(IN OUT)                  :: ia1(*)
REAL, INTENT(IN OUT)                     :: ra1(*)
REAL, INTENT(IN OUT)                     :: raP1(*)
INTEGER, INTENT(IN)                      :: iCnt
INTEGER, INTENT(IN OUT)                  :: iUD
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raArr = real array to be sorted
! iCnt  = sort indices 1..iCnt of raArr
! iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order



! local variables
INTEGER :: iTrue,iT,i1,i2
REAL :: rT

IF (iUD > 0) THEN
  i1=1
  iTrue=1
  
  60     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2=iCnt
    iTrue=-1
    70       CONTINUE
    IF (i2 > i1) THEN
      IF (raP1(i2) < raP1(i2-1)) THEN
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
  
ELSE IF (iUD < 0) THEN
  DO i1=1,iCnt
    raP1(i1)=-raP1(i1)
  END DO
  i1=1
  iTrue=1
  
  80     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2=iCnt
    iTrue=-1
    90       CONTINUE
    IF (i2 > i1) THEN
      IF (raP1(i2) < raP1(i2-1)) THEN
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
END SUBROUTINE DoSortPress

!************************************************************************
! this subroutine looks for uniqueness (by first sorting)
! depending on value of iUD
! NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!

SUBROUTINE DoUniqueReal(raArr,iCnt,iUD,rEps)


REAL, INTENT(IN OUT)                     :: raArr(*)
INTEGER, INTENT(IN OUT)                  :: iCnt
INTEGER, INTENT(IN OUT)                  :: iUD
REAL, INTENT(IN OUT)                     :: rEps
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raArr = real array to be sorted
! iCnt  = sort indices 1..iCnt of raArr
! iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order
! rEps = tolerance



! local
INTEGER :: iTrack,iI
REAL :: raArrNew(2*kProfLayer)

IF (iCnt > 2*kProfLayer) THEN
  WRITE(kStdErr,*)  'Oops, DoUniqueReal only allows iCnt <= 2*kProfLayer ',iCnt,2*kProfLayer
  WRITE(kStdWarn,*) 'Oops, DoUniqueReal only allows iCnt <= 2*kProfLayer ',iCnt,2*kProfLayer
  CALL DoStop
END IF

CALL DoSortReal(raArr,iCnt,iUD)

iTrack = 0
DO iI = 2,iCnt
  IF ((raArr(iI)-raArr(iI-1)) > rEps) THEN
    iTrack = iTrack + 1
    raArrNew(iTrack) = raArr(iI-1)
  END IF
END DO

iI = iCnt
!! could be the last few elements were the same, in which case we NEED to add last element
!      IF ((raArr(iI)-raArr(iI-1)) .GT. rEps) THEN
!        iTrack = iTrack + 1
!        raArrNew(iTrack) = raArr(iI)
!      END IF
IF ((raArr(iCnt)-raArrNew(iTrack)) > rEps) THEN
  iTrack = iTrack + 1
  raArrNew(iTrack) = raArr(iI)
END IF

IF (iTrack < iCnt) THEN
  WRITE(kStdWarn,*) 'sent in array with ',iCnt,' entries of which ',iTrack,' were unique'
  DO iI = 1,iTrack
    raArr(iI) = raArrNew(iI)
  END DO
  iCnt = iTrack
END IF

RETURN
END SUBROUTINE DoUniqueReal

!************************************************************************
! this subroutine sorts a real array into ascending order or descending order
! depending on value of iUD
! NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!
! this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

SUBROUTINE DoSortReal(raArr,iCnt,iUD)


REAL, INTENT(IN OUT)                     :: raArr(*)
INTEGER, INTENT(IN)                      :: iCnt
INTEGER, INTENT(IN OUT)                  :: iUD
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raArr = real array to be sorted
! iCnt  = sort indices 1..iCnt of raArr
! iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order



! local variables
INTEGER :: iTrue,i1,i2
REAL :: rT

IF (iUD > 0) THEN
  i1=1
  iTrue=1
  
  60     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2=iCnt
    iTrue=-1
    70       CONTINUE
    IF (i2 > i1) THEN
      IF (raArr(i2) < raArr(i2-1)) THEN
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
  
ELSE IF (iUD < 0) THEN
  DO i1=1,iCnt
    raArr(i1)=-raArr(i1)
  END DO
  i1=1
  iTrue=1
  
  80     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2=iCnt
    iTrue=-1
    90       CONTINUE
    IF (i2 > i1) THEN
      IF (raArr(i2) < raArr(i2-1)) THEN
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
END SUBROUTINE DoSortReal

!************************************************************************
! this subroutine sorts a real array into ascending order or descending order
! depending on value of iUD
! NOTE!!!!!! iUD=-1 ONLY WORKS FOR POSITIVE VALUES in raARR!!!!!!!!!!
! this is a bubble sort : refer Dale/Lilly : Pascal plus Data Structures

! so for example can send in (freq,tuning) from Scott AIRS tuning oceffs, which
! are sorted (1:2378) into AIRS channels and then (2379:2840) into arb channels
! so sort everything correctly!

SUBROUTINE DoSort2Real(raArr1,raArr2,iCnt,iUD)


REAL, INTENT(IN OUT)                     :: raArr1(*)
REAL, INTENT(IN OUT)                     :: raArr2(*)
INTEGER, INTENT(IN)                      :: iCnt
INTEGER, INTENT(IN OUT)                  :: iUD
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raArr = real array to be sorted
! iCnt  = sort indices 1..iCnt of raArr
! iUD = 1 .. sort into ascending order, iUD = -1 .. sort into decending order



! local variables
INTEGER :: iTrue,i1,i2
REAL :: rT1,rT2

IF (iUD > 0) THEN
  i1 = 1
  iTrue = 1
  
  60     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2 = iCnt
    iTrue = -1
    70       CONTINUE
    IF (i2 > i1) THEN
      IF (raArr1(i2) < raArr1(i2-1)) THEN
        iTrue = 1
        rT1 = raArr1(i2-1)
        raArr1(i2-1) = raArr1(i2)
        raArr1(i2) = rT1
        
        rT2 = raArr2(i2-1)
        raArr2(i2-1) = raArr2(i2)
        raArr2(i2) = rT2
        
      END IF
      i2 = i2-1
      GO TO 70
    END IF
    
    i1 = i1+1
    GO TO 60
  END IF
  
ELSE IF (iUD < 0) THEN
  DO i1 = 1,iCnt
    raArr1(i1) = -raArr1(i1)
    raArr2(i1) = -raArr2(i1)
  END DO
  i1 = 1
  iTrue = 1
  
  80     CONTINUE
  IF ((i1 < iCnt) .AND. (iTrue > 0)) THEN
    i2 = iCnt
    iTrue = -1
    90       CONTINUE
    IF (i2 > i1) THEN
      IF (raArr1(i2) < raArr1(i2-1)) THEN
        iTrue = 1
        
        rT1 = raArr1(i2-1)
        raArr1(i2-1) = raArr1(i2)
        raArr1(i2) = rT1
        
        rT2 = raArr2(i2-1)
        raArr2(i2-1) = raArr2(i2)
        raArr2(i2) = rT2
        
      END IF
      i2 = i2-1
      GO TO 90
    END IF
    
    i1 = i1+1
    GO TO 80
  END IF
  
  DO i1 = 1,iCnt
    raArr1(i1) = -raArr1(i1)
    raArr2(i1) = -raArr2(i1)
  END DO
  i1 = 1
END IF

RETURN
END SUBROUTINE DoSort2Real

!************************************************************************
! this function searches for the value uses a binary search
! it returns +1 if iLay is found in iaOp(1:iNp), -1 otherwise
! refer Dale/Lilly : Pascal plus Data structures

INTEGER FUNCTION BinarySearch(iLay,iNp,iaOp)


INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(*)
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay  = layer number to be looked for
! iaOp  = array containing list of layers
! iNp   = search indices 1..iNp of iaOp, to look for iLay
!      INTEGER iLay,iNp,iaOp(kPathsOut)


! local variables
INTEGER :: iAns,iFound,iMp,IF,iL,iLocation,IDIV

iAns=-1
iFound=-1
IF=1
iL=iNp

15   CONTINUE
IF ((IF <= iL) .AND. (iFound < 0)) THEN
  iMp=IDIV(IF+iL,2)
  IF (iaOp(iMp) == iLay) THEN
    iFound=1
  ELSE
    IF (iaOp(iMp) > iLay) THEN
      iL=iMp-1
    ELSE
      IF=iMp+1
    END IF
  END IF
  GO TO 15
END IF

IF (iFound > 0) THEN
  iLocation=iMp
  iAns=1
END IF

BinarySearch=iAns
RETURN
END FUNCTION BinarySearch

!************************************************************************
! this function checks to see if layer# iLay should be output
! by doing a sequential search : returns +1 if value found, -1 otherwise

INTEGER FUNCTION SequentialSearch(iLay,iNp,iaOp)


INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(*)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay  = layer number to be looked for
! iaOp  = array containing list of layers
! iNp   = search indices 1..iNp of iaOp, to look for iLay
!      INTEGER iLay,iNp,iaOp(kPathsOut)


! local variables
INTEGER :: iDp,iDpC

iDp=-1
IF (iNp < 0) THEN
! easy ! print all layers
  iDp=1
END IF
IF (iNp > 0) THEN
! actually have to go thru list to see if this layer is to be output
  iDpC=1
  101    CONTINUE
  IF ((iDpc <= iNp)  .AND. (iaOp(iDpC) == iLay)) THEN
    iDp=1
  END IF
  IF ((iDpc <= iNp)  .AND. (iDp < 0)) THEN
    iDpc=iDpc+1
    GO TO 101
  END IF
END IF

SequentialSearch=iDp

RETURN
END FUNCTION SequentialSearch

!************************************************************************
!                     SPLINES
!       SUBROUTINE rspl(rXA,rYA,N,rXOUT,rYOUT,NOUT)
!       SUBROUTINE logrspl(rXA,rYA,N,rXOUT,rYOUT,NOUT)
!       SUBROUTINE dspl(dXA,dYA,N,dXOUT,dYOUT,NOUT)
!************************************************************************

! this subroutine first sorts, then splines

SUBROUTINE r_sort_spl(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN)                         :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN OUT)                  :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox)
INTEGER :: iI

DO iI = 1,N
  XAA(iI) = XA(iI)
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

CALL rspl(XAA,YAA,N,XOUT,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_spl

!************************************************************************
! this subroutine first sorts, then splines
! changes input X, output XOUT, to log(X) and log(XOUT)

SUBROUTINE r_sort_logspl(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox),XBB(kMaxPtsBox)
INTEGER :: iI

DO iI = 1,N
  XAA(iI) = LOG(XA(iI))
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

DO iI = 1,NOUT
  XBB(iI) = LOG(XOUT(iI))
END DO

CALL rspl(XAA,YAA,N,XBB,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_logspl

!************************************************************************
! this subroutine directly calls rsply2 and then rspline

SUBROUTINE rspl1(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN OUT)                  :: N
REAL, INTENT(IN OUT)                     :: XOUT
REAL, INTENT(IN OUT)                     :: YOUT
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! real version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------

!      Parameters



REAL :: Y2A(kMaxPtsBox),work(kMaxPtsBox),Yp1,yPn
INTEGER :: I

IF (NOUT /= 1) THEN
  WRITE(kStdErr,*) 'rspl1 assumes you only want 1 calc'
  CALL DoStop
END IF

yp1=1.0E16
ypn=1.0E16

CALL rSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work)
DO I=1,NOUT
  CALL rSPLIN(XA,YA,Y2A,N,XOUT,YOUT)
END DO

RETURN
END SUBROUTINE rspl1

!************************************************************************
! this subroutine directly calls rsply2 and then rspline

SUBROUTINE rspl(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN OUT)                  :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! real version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------

!      Parameters



REAL :: Y2A(kMaxPtsBox),work(kMaxPtsBox),Yp1,yPn
INTEGER :: I

yp1=1.0E16
ypn=1.0E16

CALL rSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work)
DO I=1,NOUT
  CALL rSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I))
END DO

RETURN
END SUBROUTINE rspl

!************************************************************************
! this subroutine directly calls rsply2 and then rspline
! same as rspl except it changes XA, XOUT to log(XA),log(XOUT)

SUBROUTINE logrspl(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN)                         :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! real version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------

!      Parameters



REAL :: Y2A(kMaxPtsBox),work(kMaxPtsBox),Yp1,yPn
REAL :: LOGXA(kMaxPtsBox),LOGXOUT(kMaxPtsBox),rTemp
REAL :: raX(kMaxPtsBox),raY(kMaxPtsBox)
INTEGER :: I,iFlip,iStart

yp1=1.0E16
ypn=1.0E16

!      write(kStdErr,*) 'logrspl does not seem to work ???'
!      CALL DoStop

DO I = 1,N
  IF (XA(I) < 0.0) THEN
    WRITE(kStdErr,*) 'Error .. cannot do log(x) if x < 0'
    CALL DoStop
  ELSE
    LOGXA(I) = LOG(XA(I))
  END IF
END DO

DO I = 1,NOUT
  IF (xout(I) < 0.0) THEN
    WRITE(kStdErr,*) 'Error .. cannot do log(xout) if xout < 0'
    CALL DoStop
  ELSE
    LOGXOUT(I) = LOG(xout(I))
  END IF
END DO

iFlip = -1    !!assume everything ordered correctly
IF ((XA(N-1) > XA(N)) .OR. (XA(N-2) > XA(N-1))) THEN
  iFlip = +1
  DO I = 1,N
    raX(N-I+1) = LOGXA(I)
  END DO
  DO I = 1,N
    logxa(I) = raX(I)
  END DO
  DO I = 1,N
    raY(N-I+1) = ya(I)
  END DO
ELSE
  iFlip = -1
  DO I = 1,N
    raY(I) = ya(I)
  END DO
END IF

!!! find out where we can start from
iStart = 1
10   CONTINUE
IF (XOUT(iStart) > XA(1)) THEN
  iSTart = iStart + 1
  GO TO 10
END IF

CALL rSPLY2(logxa,raY,N,Yp1,Ypn,Y2A,work)
DO I=iStart,NOUT
  CALL rSPLIN(logxa,raY,Y2A,N,logXOUT(I),YOUT(I))
!        print *,I,logXOUT(I),YOUT(I)
END DO

RETURN
END SUBROUTINE logrspl

!************************************************************************
! this subroutine directly calls rsply2 and then rspline
! same as SUBROUTINE rspl except yp1,ypn are set differently

SUBROUTINE rspl_diffyp1n(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN OUT)                  :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! real version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------

!      Parameters



REAL :: Y2A(kMaxPtsBox),work(kMaxPtsBox),Yp1,yPn
INTEGER :: I

yp1=1.0
ypn=1.0

CALL rSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work)
DO I=1,NOUT
  CALL rSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I))
END DO

RETURN
END SUBROUTINE rspl_diffyp1n

!************************************************************************

SUBROUTINE rSPLIN(XA,YA,Y2A,N,X,Y)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
REAL, INTENT(IN)                         :: Y2A(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: X
REAL, INTENT(OUT)                        :: Y
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! REAL version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      Y2A : I  : DOUB arr : 2nd derivative of points
!      N   : I  : INT      : number of points in arrays
!      X   : I  : DOUB     : x point at which to evaluate spline
!      Y   : O  : DOUB     : y point from spline interpolation
!      -----------------------------------------------------------------

!      Parameters



!      Local Variables
INTEGER :: K,KLO,KHI
REAL :: A,B,H

!      -----------------------------------------------------------------

!      Determine between which pair of pints X falls (bisect loop)
KLO=1
KHI=N
20    IF ( (KHI - KLO) > 1) THEN
  K=(KHI + KLO)/2
  IF (XA(K) > X) THEN
    KHI = k
  ELSE
    KLO = k
  END IF
  GO TO 20
END IF

H=XA(KHI) - XA(KLO)
IF (H <= 0.0) THEN
  WRITE(kStdErr,1010) KLO,KHI,XA(KLO),XA(KHI)
  1010     FORMAT('ERROR! rSPLINT: bad XA array.',/,  &
      'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,  &
      ', XA(KHI)=',E12.5,'. Quitting.')
  CALL DoStop
END IF

A=(XA(KHI) - X)/H
B=(X - XA(KLO))/H

Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +  &
    Y2A(KHI)*(B**3 - B) )*(H**2)/6.0

RETURN
END SUBROUTINE rSPLIN

!************************************************************************

SUBROUTINE rSPLY2(XA,YA,N,YP1,YPN,Y2A,WORK)


REAL, INTENT(IN)                         :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: YP1
REAL, INTENT(IN)                         :: YPN
REAL, INTENT(OUT)                        :: Y2A(*)
REAL, INTENT(OUT)                        :: WORK(*)
IMPLICIT NONE


! REAL version
!      -----------------------------------------------------------------
!      Calc 2nd derivative as preperation for SPLINT spline routine
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      N   : I  : INT      : number of points in arrays
!      YP1 : I  : DOUB     : derivative of 1st point
!      YPN : I  : DOUB     : derivative of last point
!      Y2A : O  : DOUB arr : 2nd derivative array(N)
!      WORK: O  : DOUB arr : workspace array(N)
!      -----------------------------------------------------------------

!      Parameters



!      Local Variables
INTEGER :: I,K
REAL :: P,QN,SIG,UN

!      -----------------------------------------------------------------

!      Lower boundary
IF (YP1 > 1.0E+15) THEN
!         "Natural" boundary condition
  Y2A(1)=0.0
  WORK(1)=0.0
ELSE
!         Set to a specific first derivative
  Y2A(1)=-0.5
  WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/ (XA(2) - XA(1)) - YP1)
END IF

!      Decomposition loop of the tridiagonal algorithm
DO I=2,N-1
! this is from the progas code
!          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
  SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
  P=SIG*Y2A(I-1) + 2.0
  Y2A(I)=(SIG - 1.0)/P
  WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -  &
      (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
  WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) - SIG*WORK(I-1) )/P
  ENDDO
    
!      Upper boundary
    IF (YPN > 1.0E+15) THEN
!         "Natural" boundary condition
      QN=0.0
      UN=0.0
    ELSE
!         Set to a specific first derivative
      QN=0.5
      UN=( 3.0/(XA(N) - XA(N-1)) )*( YPN -  &
          (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
    END IF
    Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
    
!      Assign the other 2nd derivatives using the back-substitution
!      loop of the tridiagonal algorithm
    DO K=N-1,1,-1
      Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
      ENDDO
        
        RETURN
      END SUBROUTINE rSPLY2
      
!************************************************************************
      
      SUBROUTINE dspl(XA,YA,N,XOUT,YOUT,NOUT)
      
      
      DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
      DOUBLE PRECISION, INTENT(IN OUT)         :: YA(*)
      INTEGER, INTENT(IN OUT)                  :: N
      DOUBLE PRECISION, INTENT(IN OUT)         :: XOUT(*)
      DOUBLE PRECISION, INTENT(IN OUT)         :: YOUT(*)
      INTEGER, INTENT(IN)                      :: NOUT
      IMPLICIT NONE
      
! this subroutine directly calls dsply2 and then dspline
      
      INCLUDE '../INCLUDE/kcartaparam.f90'
      
! double precision version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------
      
!      Parameters
      
      
      
      DOUBLE PRECISION :: Y2A(kMaxPtsBox),work(kMaxPtsBox),Yp1,yPn
      INTEGER :: I
      
      yp1=1.0E16
      ypn=1.0E16
      
      CALL dSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work)
      DO I=1,NOUT
        CALL dSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I))
      END DO
      
      RETURN
    END SUBROUTINE dspl
    
!************************************************************************
    
    SUBROUTINE dSPLIN(XA,YA,Y2A,N,X,Y)
    
    
    DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
    DOUBLE PRECISION, INTENT(IN)             :: YA(*)
    DOUBLE PRECISION, INTENT(IN)             :: Y2A(*)
    INTEGER, INTENT(IN)                      :: N
    DOUBLE PRECISION, INTENT(IN OUT)         :: X
    DOUBLE PRECISION, INTENT(OUT)            :: Y
    IMPLICIT NONE
    
    INCLUDE '../INCLUDE/kcartaparam.f90'
    
! double precision version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      Y2A : I  : DOUB arr : 2nd derivative of points
!      N   : I  : INT      : number of points in arrays
!      X   : I  : DOUB     : x point at which to evaluate spline
!      Y   : O  : DOUB     : y point from spline interpolation
!      -----------------------------------------------------------------
    
!      Parameters
    
    
    
!      Local Variables
    INTEGER :: K,KLO,KHI
    DOUBLE PRECISION :: A,B,H
    
!      -----------------------------------------------------------------
    
!      Determine between which pair of pints X falls (bisect loop)
    KLO=1
    KHI=N
    20    IF ( (KHI - KLO) > 1) THEN
      K=(KHI + KLO)/2
      IF (XA(K) > X) THEN
        KHI = k
      ELSE
        KLO = k
      END IF
      GO TO 20
    END IF
    
    H=XA(KHI) - XA(KLO)
    IF (H <= 0.0) THEN
      WRITE(kStdErr,1010) KLO,KHI,XA(KLO),XA(KHI)
      1010     FORMAT('ERROR! dSPLINT: bad XA array.',/,  &
          'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,  &
          ', XA(KHI)=',E12.5,'. Quitting.')
      CALL DoStop
    END IF
    
    A=(XA(KHI) - X)/H
    B=(X - XA(KLO))/H
    
    Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +  &
        Y2A(KHI)*(B**3 - B) )*(H**2)/6.0
    
    RETURN
  END SUBROUTINE dSPLIN
  
!************************************************************************
  
  SUBROUTINE dSPLY2(XA,YA,N,YP1,YPN,Y2A,WORK)
  
  
  DOUBLE PRECISION, INTENT(IN)             :: XA(*)
  DOUBLE PRECISION, INTENT(IN)             :: YA(*)
  INTEGER, INTENT(IN)                      :: N
  DOUBLE PRECISION, INTENT(IN OUT)         :: YP1
  DOUBLE PRECISION, INTENT(IN)             :: YPN
  DOUBLE PRECISION, INTENT(OUT)            :: Y2A(*)
  DOUBLE PRECISION, INTENT(OUT)            :: WORK(*)
  IMPLICIT NONE
  
  
! double precision version
!      -----------------------------------------------------------------
!      Calc 2nd derivative as preperation for SPLINT spline routine
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      N   : I  : INT      : number of points in arrays
!      YP1 : I  : DOUB     : derivative of 1st point
!      YPN : I  : DOUB     : derivative of last point
!      Y2A : O  : DOUB arr : 2nd derivative array(N)
!      WORK: O  : DOUB arr : workspace array(N)
!      -----------------------------------------------------------------
  
!      Parameters
  
  
  
!      Local Variables
  INTEGER :: I,K
  DOUBLE PRECISION :: P,QN,SIG,UN
  
!      -----------------------------------------------------------------
  
!      Lower boundary
  IF (YP1 > 1.0E+15) THEN
!         "Natural" boundary condition
    Y2A(1)=0.0
    WORK(1)=0.0
  ELSE
!         Set to a specific first derivative
    Y2A(1)=-0.5
    WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/ (XA(2) - XA(1)) - YP1)
  END IF
  
!      Decomposition loop of the tridiagonal algorithm
  DO I=2,N-1
! this is from the progas code
!          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
    SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
    P=SIG*Y2A(I-1) + 2.0
    Y2A(I)=(SIG - 1.0)/P
    WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -  &
        (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
    WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) - SIG*WORK(I-1) )/P
    ENDDO
      
!      Upper boundary
      IF (YPN > 1.0E+15) THEN
!         "Natural" boundary condition
        QN=0.0
        UN=0.0
      ELSE
!         Set to a specific first derivative
        QN=0.5
        UN=( 3.0/(XA(N) - XA(N-1)) )*( YPN -  &
            (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
      END IF
      Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
      
!      Assign the other 2nd derivatives using the back-substitution
!      loop of the tridiagonal algorithm
      DO K=N-1,1,-1
        Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
        ENDDO
          
          RETURN
        END SUBROUTINE dSPLY2
        
!************************************************************************
!                       linear interps
!************************************************************************
! real linear interpolator for ONE "x" point
        
        SUBROUTINE RLINEAR1(XA,YA,N,X,Y,NOUT)
        
        
        REAL, INTENT(IN OUT)                     :: XA(*)
        REAL, INTENT(IN)                         :: YA(*)
        INTEGER, INTENT(IN)                      :: N
        REAL, INTENT(IN OUT)                     :: X
        REAL, INTENT(OUT)                        :: Y
        INTEGER, INTENT(IN OUT)                  :: NOUT
        IMPLICIT NONE
        
        INCLUDE '../INCLUDE/kcartaparam.f90'
        
! linear interpolation
! double precision version
!      -----------------------------------------------------------------
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      N   : I  : INT      : number of points in arrays
!      X   : I  : DOUB     : x point at which to evaluate spline
!      Y   : O  : DOUB     : y point from spline interpolation
!      -----------------------------------------------------------------
        
!      Parameters
        
        
        
!      Local Variables
        INTEGER :: K,KLO,KHI
        REAL :: A,B,H
        
        IF (NOUT /= 1) THEN
          WRITE(kStdErr,*) 'RLINEAR1 only gives one output'
          CALL DoStop
        END IF
        
!      -----------------------------------------------------------------
        
!      Determine between which pair of points X falls (bisect loop)
        KLO=1
        KHI=N
        
        20    IF ( (KHI - KLO) > 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) > X) THEN
            KHI = k
          ELSE
            KLO = k
          END IF
          GO TO 20
        END IF
        
        H=XA(KHI) - XA(KLO)
        IF (H <= 0.0) THEN
          WRITE(kStdErr,1010) KLO,KHI,XA(KLO),XA(KHI)
          1010     FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
              'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,  &
              ', XA(KHI)=',E12.5,'. Quitting.')
          CALL DoStop
        END IF
        
        A=(XA(KHI) - X)/H
        B=YA(KHI)-YA(KLO)
        
        Y=YA(KHI)-A*B
        
        RETURN
      END SUBROUTINE RLINEAR1
      
!************************************************************************
! this subroutine directly calls rlinear, for MANY "x" points
      
      SUBROUTINE rlinear(XA,YA,N,XOUT,YOUT,NOUT)
      
      
      REAL, INTENT(IN OUT)                     :: XA(*)
      REAL, INTENT(IN OUT)                     :: YA(*)
      INTEGER, INTENT(IN OUT)                  :: N
      REAL, INTENT(IN OUT)                     :: XOUT(*)
      REAL, INTENT(IN OUT)                     :: YOUT(*)
      INTEGER, INTENT(IN)                      :: NOUT
      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/kcartaparam.f90'
      
! real version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------
      
!      Parameters
      
      
      
      INTEGER :: I
      
      DO I=1,NOUT
        CALL RLINEAR1(XA,YA,N,XOUT(I),YOUT(I),1)
      END DO
      
      RETURN
    END SUBROUTINE rlinear
    
!************************************************************************
!  double precision - smart version; assumes XA in increasing order
!  and also assumes "x" is sent in as increasing values each time
! "x" comes in singly
    
    SUBROUTINE dLINEAR_SMART_ONE(NSMART_LO,NSMART_HI,XA,YA,N,X,Y)
    
    
    INTEGER, INTENT(OUT)                     :: NSMART_LO
    INTEGER, INTENT(OUT)                     :: NSMART_HI
    DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
    DOUBLE PRECISION, INTENT(IN)             :: YA(*)
    INTEGER, INTENT(IN)                      :: N
    DOUBLE PRECISION, INTENT(IN OUT)         :: X
    DOUBLE PRECISION, INTENT(OUT)            :: Y
    IMPLICIT NONE
    
    INCLUDE '../INCLUDE/kcartaparam.f90'
    
! linear interpolation
! double precision version
!      -----------------------------------------------------------------
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      N   : I  : INT      : number of points in arrays
!      X   : I  : DOUB     : x point at which to evaluate spline
!      Y   : O  : DOUB     : y point from spline interpolation
!  NSMART  : I/O : INT     : guess at where the bisect occurs
!                              if set at -1, then do usual trial
!                              if set elsewehre, see if this works!
!      -----------------------------------------------------------------
    
!      Parameters
    
    
    
!      Local Variables
    INTEGER :: K,KLO,KHI
    DOUBLE PRECISION :: A,B,H
    
!      -----------------------------------------------------------------
    
!      Determine between which pair of points X falls (bisect loop)
    IF ((NSMART_LO < 0) .AND. (NSMART_HI < 0)) THEN
!!!start looking from beginning
      KLO=1
      KHI=N
      20      IF ( (KHI - KLO) > 1) THEN
        K=(KHI + KLO)/2
        IF (XA(K) > X) THEN
          KHI = k
        ELSE
          KLO = k
        END IF
        GO TO 20
      END IF
      NSMART_LO = KLO
      NSMART_HI = KHI
    ELSE
!!!start looking from previous try!
      KLO=MAX(NSMART_LO - 1,1)
      KHI=MIN(NSMART_HI + 1,N)
      30      IF ( (KHI - KLO) > 1) THEN
        K=(KHI + KLO)/2
        IF (XA(K) > X) THEN
          KHI = k
        ELSE
          KLO = k
        END IF
        GO TO 30
      END IF
      NSMART_LO = KLO
      NSMART_HI = KHI
    END IF
    
    H=XA(KHI) - XA(KLO)
    IF (H <= 0.0) THEN
      WRITE(kStdErr,1010) KLO,KHI,XA(KLO),XA(KHI)
      1010     FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
          'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,  &
          ', XA(KHI)=',E12.5,'. Quitting.')
      CALL DoStop
    END IF
    
    A=(XA(KHI) - X)/H
    B=YA(KHI)-YA(KLO)
    
    Y=YA(KHI)-A*B
    
    RETURN
  END SUBROUTINE dLINEAR_SMART_ONE
  
!************************************************************************
! this subroutine directly calls linear, in a smart way
  
  SUBROUTINE dlinear_smart(XA,YA,N,XOUT,YOUT,NOUT)
  
  
  DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
  DOUBLE PRECISION, INTENT(IN OUT)         :: YA(*)
  INTEGER, INTENT(IN OUT)                  :: N
  DOUBLE PRECISION, INTENT(IN OUT)         :: XOUT(*)
  DOUBLE PRECISION, INTENT(IN OUT)         :: YOUT(*)
  INTEGER, INTENT(IN)                      :: NOUT
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  
! double precision version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------
  
!      Parameters
  
  INTEGER :: NSMART_LO,NSMART_HI
  
  INTEGER :: I
  
  I = 1
  NSMART_LO = -1
  NSMART_HI = -1
  CALL dLINEAR_smart_one(NSMART_LO,NSMART_HI,XA,YA,N,XOUT(I),YOUT(I))
  
  DO I=2,NOUT
    CALL dLINEAR_smart_one(NSMART_LO,NSMART_HI,XA,YA,N,XOUT(I),YOUT(I))
  END DO
  
  RETURN
END SUBROUTINE dlinear_smart

!************************************************************************
!  double precision

SUBROUTINE dLINEAR_ONE(XA,YA,N,X,Y)


DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
DOUBLE PRECISION, INTENT(IN)             :: YA(*)
INTEGER, INTENT(IN)                      :: N
DOUBLE PRECISION, INTENT(IN OUT)         :: X
DOUBLE PRECISION, INTENT(OUT)            :: Y
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
! double precision version
!      -----------------------------------------------------------------
!      XA  : I  : DOUB arr : x array(N) in increasing order
!      YA  : I  : DOUB arr : y array(N)
!      N   : I  : INT      : number of points in arrays
!      X   : I  : DOUB     : x point at which to evaluate spline
!      Y   : O  : DOUB     : y point from spline interpolation
!      -----------------------------------------------------------------

!      Parameters



!      Local Variables
INTEGER :: K,KLO,KHI
DOUBLE PRECISION :: A,B,H

!      -----------------------------------------------------------------

!      Determine between which pair of points X falls (bisect loop)
KLO=1
KHI=N
20    IF ( (KHI - KLO) > 1) THEN
  K=(KHI + KLO)/2
  IF (XA(K) > X) THEN
    KHI = k
  ELSE
    KLO = k
  END IF
  GO TO 20
END IF

H=XA(KHI) - XA(KLO)
IF (H <= 0.0) THEN
  WRITE(kStdErr,1010) KLO,KHI,XA(KLO),XA(KHI)
  1010     FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
      'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,  &
      ', XA(KHI)=',E12.5,'. Quitting.')
  CALL DoStop
END IF

A=(XA(KHI) - X)/H
B=YA(KHI)-YA(KLO)

Y=YA(KHI)-A*B

RETURN
END SUBROUTINE dLINEAR_ONE

!************************************************************************
! this subroutine directly calls linear

SUBROUTINE dlinear(XA,YA,N,XOUT,YOUT,NOUT)


DOUBLE PRECISION, INTENT(IN OUT)         :: XA(*)
DOUBLE PRECISION, INTENT(IN OUT)         :: YA(*)
INTEGER, INTENT(IN OUT)                  :: N
DOUBLE PRECISION, INTENT(IN OUT)         :: XOUT(*)
DOUBLE PRECISION, INTENT(IN OUT)         :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! double precision version
!      -----------------------------------------------------------------
!      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
!      XA  : I  : DOUB arr : x array(N) in increasing order  IN
!      YA  : I  : DOUB arr : y array(N)                      IN
!      N   : I  : INT      : number of points in arrays
!      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline
!      YOUT   : O  : DOUB ARR     : y points from spline interpolation
!      NOUT   : I  : INT          : number of points at which to spline
!      -----------------------------------------------------------------

!      Parameters



INTEGER :: I

DO I=1,NOUT
  CALL dLINEAR_ONE(XA,YA,N,XOUT(I),YOUT(I))
END DO

RETURN
END SUBROUTINE dlinear

!************************************************************************
! this subroutine first sorts, then interps

SUBROUTINE r_sort_linear(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN)                         :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN OUT)                  :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox)
INTEGER :: iI

DO iI = 1,N
  XAA(iI) = XA(iI)
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

CALL rlinear(XAA,YAA,N,XOUT,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_linear

!************************************************************************
! this subroutine first sorts, then linearines
! changes input X, output XOUT, to log(X) and log(XOUT)

SUBROUTINE r_sort_loglinear(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT(*)
REAL, INTENT(IN OUT)                     :: YOUT(*)
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox),XBB(kMaxPtsBox)
INTEGER :: iI

DO iI = 1,N
  XAA(iI) = LOG(XA(iI))
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

DO iI = 1,NOUT
  XBB(iI) = LOG(XOUT(iI))
END DO

CALL rlinear(XAA,YAA,N,XBB,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_loglinear

!************************************************************************
! this subroutine first sorts, then interps

SUBROUTINE r_sort_linear1(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN)                         :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT
REAL, INTENT(IN OUT)                     :: YOUT
INTEGER, INTENT(IN OUT)                  :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox)
INTEGER :: iI

IF (NOUT /= 1) THEN
  WRITE(kStdErr,*) 'r_sort_linear1 is only for 1 point'
  CALL DoStop
END IF

DO iI = 1,N
  XAA(iI) = XA(iI)
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

CALL rlinear1(XAA,YAA,N,XOUT,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_linear1

!************************************************************************
! this subroutine first sorts, then linearines
! changes input X, output XOUT, to log(X) and log(XOUT)

SUBROUTINE r_sort_loglinear1(XA,YA,N,XOUT,YOUT,NOUT)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN OUT)                     :: XOUT
REAL, INTENT(IN OUT)                     :: YOUT
INTEGER, INTENT(IN)                      :: NOUT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!      Parameters



REAL :: XAA(kMaxPtsBox),YAA(kMaxPtsBox),XBB
INTEGER :: iI

IF (NOUT /= 1) THEN
  WRITE(kStdErr,*) 'r_sort_loglinear1 is only for 1 point'
  CALL DoStop
END IF

DO iI = 1,N
  XAA(iI) = LOG(XA(iI))
  YAA(iI) = YA(iI)
END DO
CALL DoSort2Real(XAA,YAA,N,1)

DO iI = 1,NOUT
  XBB = LOG(XOUT)
END DO

CALL rlinear1(XAA,YAA,N,XBB,YOUT,NOUT)

RETURN
END SUBROUTINE r_sort_loglinear1

!************************************************************************
