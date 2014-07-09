c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************ 
      SUBROUTINE printstar

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      write(kStdWarn,*)
     $'**********************************************************************'

      RETURN
      END

c************************************************************************ 
      SUBROUTINE printpound

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      write(kStdWarn,*)
     $'######################################################################'

      RETURN
      END

c************************************************************************ 
c this subroutine conactenates two ca80 strings together
c ca1 + ca2 = ca2       (storing the result in ca2)
      SUBROUTINE ConcatCA80(ca1,ca2)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 ca1,ca2
      INTEGER iGasID

c local variables
      INTEGER iI,i1,i2,iJ
      CHARACTER*80 caSum

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)

c do some initializations
      iI=80
 11   CONTINUE
      IF ((ca1(iI:iI) .EQ. ' ') .AND. (iI .GE. 1)) THEN
        iI=iI-1
        GO TO 11
        END IF
      i1=iI

      iI=80
 12   CONTINUE
      IF ((ca2(iI:iI) .EQ. ' ') .AND. (iI .GE. 1)) THEN
        iI=iI-1
        GO TO 12
        END IF
      i2=iI

      DO iI=1,80
        caSum(iI:iI) = ' '
        END DO

      DO iI=1,i1
        caSum(iI:iI) = ca1(iI:iI)
        END DO

      iJ = i1+1
      DO iI=1,i2
        caSum(iJ:iJ) = ca2(iI:iI)
        iJ = iJ + 1
        END DO

      DO iI=1,80
        ca2(iI:iI) = caSum(iI:iI)
        END DO

      RETURN
      END
c************************************************************************ 
c this subroutine finds the name of the CKD binary file
      SUBROUTINE CKDFileName(caFName,iGasID,iTag)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID      = current ID of gas to be processed (101=self,102=foreign)
c caFName     = name of file that contains CKD data
      CHARACTER*80 caFName
      INTEGER iGasID,iTag

c local variables
      CHARACTER*2 caString,caTemp
      CHARACTER*4 caTemp4
      CHARACTER*80 caDir
      INTEGER iInt,iLenDir,iL

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)
      IF (iTag .EQ. 2) THEN 
        !! 605-2830 cm-1   'r' prefix 
        caDir = kCKDPath          !open AIRS 100 layer reference profile 
      ELSEIF (iTag .EQ. 3) THEN 
        !! 4050-4950 cm-1   'm' prefix 
        caDir = kCKDPathm          !open AIRS 100 layer reference profile 
      ELSE 
        write(kStdErr,*) 'So far have only coded up iTag = 2,3' 
        write(kStdErr,*) 'Cannot find compressed datafilename for other bands' 
        CALL DoStop          
        END IF 

c do some initializations
      iInt=80
 11   CONTINUE
      IF ((caDir(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN
        iInt=iInt-1
        GO TO 11
        END IF
      iLenDir=iInt

      DO iInt=1,80
        caFName=' '
        END DO

      caString='  '
      caTemp='00'

c now process iGasID so that we end up with a right padded string
c eg iGasID=101 ==> 'Self', 102 ==> For
      caTemp4='    '
      IF (iGasID .EQ. kNewGasLo) THEN
        caTemp4='Self'
        iL=4
      ELSE
        caTemp4='For'
        iL=3
        END IF

c CKD         = CKD version 0,21,23,24
c now process CKD ver so that we end up with a right padded string
c eg 2 ---> '2 ', 12 ---> '12' etc
      WRITE(caString,15) kCKD
 15      FORMAT(I2)
c this is right justified ... change to left justified
      iInt=1
 16      continue
      IF (caString(iInt:iInt) .eq. ' ') THEN
        iInt=iInt+1
        GO TO 16
        END IF
      caTemp(1:2-iInt+1)=caString(iInt:2)

c now concatenate all this together in caFname
      caFName(1:iLendir)=caDir(1:iLenDir)
      caFName(iLenDir+1:iLenDir+3)='CKD'
      iLenDir=iLenDir+3                           !we have added 3 characters
      caFName(iLenDir+1:iLenDir+iL)=caTemp4(1:iL) !add on 'Self' or 'For'
      iLenDir=iLenDir+iL                          !we have added iL characters
      IF (kCKD .GE. 10) THEN
        caFName(iLenDir+1:iLenDir+2)=caTemp(1:2)
        iLenDir=iLenDir+2                           !we have added 2 characters
      ELSEIF (kCKD .LT. 10) THEN
        caFName(iLenDir+1:iLenDir+1)=caTemp(1:1)
        iLenDir=iLenDir+1                           !we have added 1 character
        END IF
      caFName(iLenDir+1:iLenDir+4)='.bin'

      WRITE (kStdWarn,2000) caFName
 2000  FORMAT('CKD datafile name is : ',/,A80)

      RETURN
      END
c************************************************************************ 

c this subroutine finds the name of the reference profile for the current gas
      SUBROUTINE FindReferenceName(caFName,iGasID,iNewOrAirs)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNewOrAirs  = open new klayers (kProfLayer) or AIRS (100 layer) file (+/-1)
c iGasID      = current ID of gas to be processed
c caFName     = name of file that contains reference profile
      CHARACTER*80 caFName
      INTEGER iGasID,iNewOrAirs

c local variables
      CHARACTER*2 caString,caTemp
      CHARACTER*80 caDir
      INTEGER iInt,iLenDir

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)
      IF (iNewOrAirs .GT. 0) THEN
        !!!caDir=kNewRefPath       !open new kLAYERS created reference profile
        write(kStdErr,*) 'Only using original 100 AIRS profiles'
        CALL DoStop
      ELSE
        caDir=kOrigRefPath      !open original AIRS100 layer reference profile
        END IF

c do some initializations
      iInt=80
 11   CONTINUE
      IF ((caDir(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN
        iInt=iInt-1
        GO TO 11
        END IF
      iLenDir=iInt

      DO iInt=1,80
        caFName=' '
        END DO

      caString='  '
      caTemp='  '

c now process iGasID so that we end up with a right padded string
c eg 2 ---> '2 ', 12 ---> '12' etc
      WRITE(caString,15) iGasID
 15      FORMAT(I2)
c this is right justified ... change to left justified
      iInt=1
 16      continue
      IF (caString(iInt:iInt) .eq. ' ') THEN
        iInt=iInt+1
        GO TO 16
        END IF
      caTemp(1:2-iInt+1)=caString(iInt:2)

c now concatenate all this together in caFname
      caFName(1:iLendir)=caDir(1:iLenDir)
      caFName(iLenDir+1:iLenDir+6)='refgas'
      caFName(iLenDir+7:iLenDir+8)=caTemp(1:2)

      IF (iNewORAirs .LT. 0) THEN
        WRITE (kStdWarn,2000) caFName
      ELSE
        WRITE (kStdWarn,2010) caFName
        ENDIF
 2000 FORMAT('100 AIRS layer Reference datafile name is : ',/,A80)
 2010 FORMAT('kProfLayer     Reference datafile name is : ',/,A80)

      RETURN
      END

c************************************************************************
c for 1 <= iGasID  <= kGasComp
c thus if we send in (2,430,1) then caFname='q430_g2.dat   '
      SUBROUTINE CompFileNameOld(iGasID,iFileStartFr,iTag,caFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = ID of gas to be processed
c iFileStartFr  = start freq of file (integer)
c caFName  = final name of file we wish to read in 
c iTag     = depending on (1,2,3), this is q r or s
      CHARACTER*80 caFName
      INTEGER iGasID,iFileStartFr,iTag
      
c local variables
      CHARACTER*13 caTemp
      CHARACTER*67 caDir
      CHARACTER*2 caString2,caTemp2
      CHARACTER*4 caString4,caTemp4
      CHARACTER*6 caTemp6
      CHARACTER*8 caTemp8
      INTEGER iInt,iLen2,iLen4,iLenDir

ccccc user supplied info
c here give the directory path to where the compressed files are ...
c keeping in mind that the water files may be in a directory by themselves
      IF (iGasID .EQ. 1) THEN
         caDir=kWaterPath
      ELSEIF (iGasID .EQ. 2) THEN
         caDir=kCO2Path
      ELSE
         caDir=kCompPath
        END IF

      iInt=67
 11      CONTINUE
      IF ((caDir(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN
        iInt=iInt-1
        GO TO 11
        END IF
      iLenDir=iInt

      DO iInt=1,80
        caFName=' '
        END DO
      DO iInt=1,15
        caTemp=' '
        END DO
      DO iInt=1,2
        caString2(iInt:iInt)=' '
        caTemp2(iInt:iInt)=' '
        END DO
      DO iInt=1,4
        caString4(iInt:iInt)=' '
        caTemp4(iInt:iInt)=' '
        END DO
      caTemp6='     '
      caTemp8='       '
 
c now process iFileStartFr so that we end up with a right padded string 
c eg 605 ---> '605 ', 1200 ---> '1200' etc
      WRITE(caString4,12) iFileStartFr
 12      FORMAT(I4)
c this is right justified ... change to left justified
      iInt=1
 14   continue
      IF (caString4(iInt:iInt) .eq. ' ') THEN
        iInt=iInt+1
        GO TO 14
        END IF
      caTemp4(1:4-iInt+1)=caString4(iInt:4)
      iLen4=4-iInt+1

c now process iGasID so that we end up with a right padded string to which
c we append '.dat'
c eg 2 ---> '2.dat ', 12 ---> '12.dat' etc
      WRITE(caString2,15) iGasID
 15      FORMAT(I2)
c this is right justified ... change to left justified
      iInt=1
 16      continue
      IF (caString2(iInt:iInt) .eq. ' ') THEN
        iInt=iInt+1
        GO TO 16
        END IF
      caTemp2(1:2-iInt+1)=caString2(iInt:2)
      iLen2=2-iInt+1
      IF (iLen2 .EQ. 1) THEN
        caTemp6(1:1)=caTemp2(1:1)
        caTemp6(2:5)='.dat'
      ELSE
        caTemp6(1:2)=caTemp2(1:2)
        caTemp6(3:6)='.dat'
        END IF
c now prepend '_g'
      caTemp8(1:2)='_g'
      caTemp8(3:8)=caTemp6

c using caTemp, put together the compressed data file name, starting with 'r'
      IF (iTag .EQ. 1) THEN 
        caTemp(1:1)='q'
      ELSE IF (iTag .EQ. 2) THEN 
        caTemp(1:1)='r'
      ELSE IF (iTag .EQ. 3) THEN 
        caTemp(1:1)='s'
        END IF

c followed by the iFileStartFr
      caTemp(2:5)=caTemp4(1:4)
c followed by _g + IGASID + .dat
      IF (iFileStartFr .GT. 1000) THEN
        caTemp(6:13)=caTemp8(1:8)
      ELSE
        caTemp(5:12)=caTemp8(1:8)
      END IF

c now concatenate directory path and filename : caDir+caTemp
      IF (iLenDir .GT. 0) THEN
        caFName(1:iLenDir)=caDir(1:iLenDir)
        END IF
      caFName(iLenDir+1:iLenDir+1+13)=caTemp(1:13)      

c      WRITE (6,2000) caFName
c 2000 FORMAT('Compressed datafile name is : ',/,A80)
      RETURN
      END

c************************************************************************
c for 1 <= iGasID  <= kGasComp 
c thus if we send in (+1,2,430,1) then caFname='q430_g2.dat   ' 
c also has to add on the DIRECTORY 
c only worry about linemixing/COUSIN/NLTE for CO2 in the 605-2830 cm-1 tag  
c (iTag=2,prefix='r') 
      SUBROUTINE CompFileName(iLineMix,iGasID,iFileStartFr,iTag,caFName) 
 
      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 
 
c iGasID   = ID of gas to be processed 
c iFileStartFr  = start freq of file (integer) 
c caFName  = final name of file we wish to read in  
c iTag     = depending on (1,2,3), this is q r or s 
c iLineMix = for CO2, can do for 1100 - 0.005 mb :  
c                            iLineMix = +1 for the linemixing files, or 
c                            iLineMix = -1 for the Cousin files, or 
c                            iLineMix = -2 for the weak backgnd 4 um CO2 files 
c                     can do for 0.005 - 0.00025 mb 
c                            iLineMix = -3 for the weak backgnd 4 um CO2 files 
c            for all other gases, iLineMix = irrelevant (set to +1 usually) 
 
      CHARACTER*80 caFName 
      INTEGER iGasID,iFileStartFr,iTag,iLineMix 
       
c local variables 
      CHARACTER*13 caTemp 
      CHARACTER*67 caDir,caDirUse 
      CHARACTER*2 caString2,caTemp2 
      CHARACTER*4 caString4,caTemp4 
      CHARACTER*6 caTemp6 
      CHARACTER*8 caTemp8 
      INTEGER iInt,iLen2,iLen4,iLenDir 

ccccc user supplied info 
c here give the directory path to where the compressed files are ... 
c keeping in mind that the water files may be in a directory by themselves 
      IF (iTag .EQ. 3) THEN 
        !! -------> 4050-4950 cm-1  <---------------------   'm' prefix 
        IF (iGasID .EQ. 1) THEN 
          caDir = kWaterPathm 
        ELSE 
          caDir = kCompPathm 
          END IF 
      ELSEIF (iTag .EQ. 2) THEN 
        !! -------> 0605-2830 cm-1  <---------------------   'r' prefix 
        IF (iGasID .EQ. 1) THEN 
          caDir = kWaterPath 
        ELSEIF (iGasID .EQ. 2) THEN 
          IF (iLineMix .EQ. +1) THEN 
            !!use linemixing 
            caDir = kCO2Path 
          ELSEIF (iLineMix .EQ. -1) THEN 
            !!use cousin lineshapes from GENLN2 
            caDir = kCousin_CO2Path 
          ELSEIF (iLineMix .EQ. -2) THEN 
            !!use weak background, plus 2353,2354 and 2321 lines 
            caDir = caWeakCO2Path 
          ELSEIF (iLineMix .EQ. -3) THEN 
            !!use weak background, plus 2353,2354 and 2321 lines for upperatm 
            caDir = caWeakUpperAtmCO2Path 
          ELSE 
            write(kStdErr,*) 'Wrong option for CO2 in CompFileName',iLineMix 
            CALL DoStop 
            END IF 
        ELSE 
          caDir = kCompPath 
          END IF             
      ELSE 
        write(kStdErr,*) 'So far have only coded up iTag = 2,3' 
        write(kStdErr,*) 'Cannot find compressed datafilenames for other bands' 
        CALL DoStop          
        END IF 
 
      iInt=67 
 11      CONTINUE 
      IF ((caDir(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
        END IF 
      iLenDir=iInt 
 
      DO iInt=1,80 
        caFName=' ' 
        END DO 
      DO iInt=1,13 
        caTemp=' ' 
        END DO 
      DO iInt=1,2 
        caString2(iInt:iInt)=' ' 
        caTemp2(iInt:iInt)=' ' 
        END DO 
      DO iInt=1,4 
        caString4(iInt:iInt)=' ' 
        caTemp4(iInt:iInt)=' ' 
        END DO 
      caTemp6='     ' 
      caTemp8='       ' 

c now process iFileStartFr so that we end up with a right padded string  
c eg 605 ---> '605 ', 1200 ---> '1200' etc 
      WRITE(caString4,12) iFileStartFr 
 12         FORMAT(I4) 
c this is right justified ... change to left justified 
      iInt=1 
 14      continue 
      IF (caString4(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 14 
        END IF 
      caTemp4(1:4-iInt+1)=caString4(iInt:4) 
      iLen4=4-iInt+1 
 
c now process iGasID so that we end up with a right padded string to which 
c we append '.dat' 
c eg 2 ---> '2.dat ', 12 ---> '12.dat' etc 
      WRITE(caString2,15) iGasID 
 15         FORMAT(I2) 
c this is right justified ... change to left justified 
      iInt=1 
 16         continue 
      IF (caString2(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 16 
        END IF 
      caTemp2(1:2-iInt+1)=caString2(iInt:2) 
      iLen2=2-iInt+1 
      IF (iLen2 .EQ. 1) THEN 
        caTemp6(1:1)=caTemp2(1:1) 
        caTemp6(2:5)='.dat' 
      ELSE 
        caTemp6(1:2)=caTemp2(1:2) 
        caTemp6(3:6)='.dat' 
        END IF 
c now prepend '_g' 
      caTemp8(1:2)='_g' 
      caTemp8(3:8)=caTemp6 
 
c using caTemp, put together the compressed data file name, starting with 'r' 
      IF (iTag .EQ. 1) THEN  
        caTemp(1:1)='q' 
      ELSE IF (iTag .EQ. 2) THEN  
        caTemp(1:1)='r' 
      ELSE IF (iTag .EQ. 3) THEN  
        caTemp(1:1)='m' 
      ELSE IF (iTag .EQ. 4) THEN  
        caTemp(1:1)='n' 
      ELSE IF (iTag .EQ. 5) THEN  
        caTemp(1:1)='v' 
      ELSE IF (iTag .EQ. 6) THEN  
        caTemp(1:1)='w' 
        END IF 
 
      caTemp(1:1) = kaCtag(iTag) 
 
c followed by the iFileStartFr 
      caTemp(2:5)=caTemp4(1:4) 
c followed by _g + IGASID + .dat 
      IF (iFileStartFr .GT. 1000) THEN 
        caTemp(6:13)=caTemp8(1:8) 
      ELSE 
        caTemp(5:12)=caTemp8(1:8) 
      END IF 

c now concatenate directory path and filename : caDir+caTemp 
      IF (iLenDir .GT. 0) THEN 
        caFName(1:iLenDir)=caDir(1:iLenDir) 
        END IF 
      caFName(iLenDir+1:iLenDir+1+13)=caTemp(1:13)       
 
c      WRITE (6,2000) caFName 
c 2000 FORMAT('Compressed datafile name is : ',/,A80) 
 
      RETURN 
      END 

c************************************************************************
c this finds the overall CHI function file name
      SUBROUTINE FIndChiFileName(caFile)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*120 caFile    !!!input and output parameter

      CHARACTER*120 caTemp,caDir
      INTEGER iI,iJ

      DO iI = 1,120
        caTemp(iI:iI) = ' '
        END DO

      caDir = kChiFIle
      iI = 120
 10   CONTINUE
      IF ((caDir(iI:iI) .EQ. ' ') .AND. (iI .GT. 1)) THEN
        iI = iI - 1
        GOTO 10
        ENDIF
      caTemp(1:iI) = caDir(1:iI)

      iJ = 120
 20   CONTINUE
      IF ((caFile(iJ:iJ) .EQ. ' ') .AND. (iJ .GT. 1)) THEN
        iJ = iJ - 1
        GOTO 20
        ENDIF
      caTemp(iI+1:iI+1+iJ) = caFile(1:iJ)

      DO iI = 1,120
        caFile(iI:iI) = caTemp(iI:iI)
        END DO

      write (kStdWarn,30) caFile
 30   FORMAT('chifile = ',A120)

      RETURN
      END

c************************************************************************
c this function finds the pressure layer at which rPressStart is within,
c as well as the fraction of the layer that it occupies
      REAL FUNCTION FindBottomTemp(rP,raProfileTemp,
     $                             raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels = actual pressure levels that come out of kLAYERS
c raProfileTemp = actual profile temp
c rP            = pressure at which we want the temperature
      real rP,raProfileTemp(kProfLayer),raPressLevels(kProfLayer+1)
      integer iProfileLayers

      integer iFound,i1,i2,i3,iLowest
      real rP1,rP2,T1,T2
      real raP(3),raT(3),Y2A(3),rT
      real yp1,ypn,work(3)

      rT=0.0

      iLowest = kProfLayer-iProfileLayers+1

      if (rP .ge. raPressLevels(iLowest)) then
        !this is WHOLE of the bottom layer
        i1=iLowest
      else if (rP .le. raPressLevels(kProfLayer+1)) then
        !this is ludicrous
        write(kStdErr,*) rP,raPressLevels(kProfLayer+1)
        write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
        CALL DoStop

      else
        !first find the AIRS layer within which it lies
        iFound=-1 
        i1 = iLowest
        i2 = iLowest+1
 10     CONTINUE 
        IF ((rP .LE.raPressLevels(i1)).AND.(rP .GT.raPressLevels(i2))) THEN 
          iFound=1 
          END IF 
        IF ((iFound .LT. 0) .AND. (i1 .LT. kProfLayer)) THEN 
          i1=i1+1 
          i2=i2+1 
          GO TO 10 
          END IF 
        IF ((iFound .LT. 0)) THEN 
          IF (abs(rP-raPressLevels(kProfLayer+1)) .LE. delta) THEN 
            i1=kProfLayer 
            iFound=1 
          ELSE 
            write(kStdErr,*) 'could not find pressure ',rP 
            write(kStdErr,*) 'within AIRS pressure levels. Please check' 
            write(kStdErr,*) '*RADNCE and *OUTPUT sections' 
            CALL DoSTOP 
            END IF 
          END IF 
        END IF

      IF ((i1 .gt. kProfLayer) .OR. (i1 .lt. iLowest)) THEN
        write(kStdErr,*) 'sorry : cannot find surface temp for '
        write(kStdErr,*) 'layers outside ',iLowest,' and ',kProfLayer
        write(kStdErr,*) 'Allowed Pressure ranges are from : ',
     $ raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
        write(kStdErr,*) 'Surface Pressure is ',rP,' mb'
        call DoStop
        END IF 
          
      !now find the temperature
      if (i1 .EQ. iLowest) then          !do linear interp
        i1 = iLowest
        i2 = iLowest+1
        i3 = iLowest+2
        rP1=(raPressLevels(i2)-raPressLevels(i1))/
     $          log(raPressLevels(i2)/raPressLevels(i1))
        rP2=(raPressLevels(i3)-raPressLevels(i2))/
     $          log(raPressLevels(i3)/raPressLevels(i2))
        T1=raProfileTemp(i1)
        T2=raProfileTemp(i2)
        rT=T2-(rP2-rP)*(T2-T1)/(rP2-rP1)
      elseif (i1 .GE. (kProfLayer-1)) then          !do linear interp
        rP1=(raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/
     $          log(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
        rP2=(raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/
     $          log(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
        T1=raProfileTemp(kProfLayer-1)
        T2=raProfileTemp(kProfLayer)
        rT=T2-(rP2-rP)*(T2-T1)/(rP2-rP1)
      else          !do spline ... note that the pressures have to 
                    !be in ascending order for good interpolation
        rP1=(raPressLevels(i1)-raPressLevels(i1-1))/
     $          log(raPressLevels(i1)/raPressLevels(i1-1))
        raP(3)=rP1
        rP1=(raPressLevels(i1+1)-raPressLevels(i1))/
     $          log(raPressLevels(i1+1)/raPressLevels(i1))
        raP(2)=rP1
        rP1=(raPressLevels(i1+2)-raPressLevels(i1+1))/
     $          log(raPressLevels(i1+2)/raPressLevels(i1+1))
        raP(1)=rP1
        raT(3)=raProfileTemp(i1-1)
        raT(2)=raProfileTemp(i1)
        raT(1)=raProfileTemp(i1+1)

        yp1=1.0e30
        ypn=1.0e30
        CALL rSPLY2(raP,raT,3,YP1,YPN,Y2A,WORK)
        call rSplin(raP,raT,Y2A,3,rP,rT)

        end if

      FindBottomTemp=rT

      RETURN
      END
c************************************************************************
c this subroutine will take in 100 AIRS layering stuff and interpolate to
c the new arbitrary layering
c WARNING : this assumes that the user has not mucked up KLAYERS layering
c           such that higehst Z pressure (lowest pressure) is NOT TOA
c           ie still need lowest pressure (highest z) = 0.005 mb!!!!!

      SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $           raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $           raaPress,iGas,iGasID,iNumLayers,iError)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c these are the individual reference profiles, at kProfLayer layers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c these are the individual reference profiles, at kMaxLayer layers
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
c these are the arbitrary profiles stored in matrices
      REAL raaPress(kProfLayer,kGasStore)
      INTEGER iError,iGas,iGasID,iNumLayers

c local variables
      INTEGER iI,iNot
      REAL raWorkP(kMaxLayer),raXgivenP(kMaxLayer),
     $     raYgivenP(kMaxLayer),raY2P(kMaxLayer)
      REAL raWork(kMaxTemp),rYP1,rYPN,rXPT,r

C     Assign values for interpolation
C     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
      rYP1=1.0E+16
      rYPN=1.0E+16

      !!!this tells how many layers are NOT dumped out by kLAYERS
      iNot = kProfLayer-iNumLayers

c now just happily spline everything on!!!!!! for the amts
c recall you need raXgivenP to be increasing
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100Amt(kMaxLayer-iI+1)
        END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        r = 0.0
        raRAmt(iI) = r
        END DO
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        raRAmt(iI) = r
        IF ((r .LT. 0.0)) THEN
          write (kStdErr,*) 'In creating reference profile, negative'
          write (kStdErr,*) 'gas amount found!!!! GasID = ',iGasID
          CALL DoStop
          END IF
        END DO

c now just happily spline everything on!!!!!! for the temps
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
        END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        r = 273.15
        raRTemp(iI) = r
        END DO
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        raRTemp(iI) = r
        END DO

c now just happily spline everything on!!!!!! for the partial pressures
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100PartPress(kMaxLayer-iI+1)
        END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        r = 0.0
        raRPartPress(iI) = r
        END DO
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        IF ((r .LT. 0.0) .AND. (iGasID .EQ. 1)) THEN
          write (kStdErr,*) 'In creating water reference profile, negative'
          write (kStdErr,*) 'gas partial pressure found!!!!'
          CALL DoStop
        ELSEIF ((r .LT. 0.0) .AND. (iGasID .NE. 1)) THEN
          write (kStdWarn,*) 'Warning!!In creating reference profile,negative'
          write (kStdWarn,*) 'gas partial pressure found!!!! Reset to 0'
          write (kStdWarn,*) 'Gas ID, layer, PP = ',iGasID,iI,r
          r = 0.0
          END IF
        raRPartPress(iI) = r
        END DO

c simply put in the pressures
      DO iI = 1,iNot
        raRPress(iI) = raaPress(iNot+1,iGas)
        END DO
      DO iI = iNot+1,kProfLayer
        raRPress(iI) = raaPress(iI,iGas)
        END DO

c      do iI = iNot+1,kProfLayer
c        print *,iGas,iI,raRAmt(iI)/raR100Amt(iI),raRTemp(iI)/raR100Temp(iI),
c     $                    raRPress(iI)/raR100Press(iI),
c     $                    raRPartPress(iI)/raR100PartPress(iI)
c        end do
c      print *,' '

      RETURN
      END
c************************************************************************
