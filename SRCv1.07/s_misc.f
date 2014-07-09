c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************ 
      SUBROUTINE printstar

      include 'kcarta.param'

      write(kStdWarn,*)
     $'**********************************************************************'

      RETURN
      END

c************************************************************************ 
      SUBROUTINE printpound

      include 'kcarta.param'

      write(kStdWarn,*)
     $'######################################################################'

      RETURN
      END

c************************************************************************ 
c this subroutine finds the name of the CKD binary file
      SUBROUTINE CKDFileName(caFName,iGasID)

      include 'kcarta.param'

c iGasID      = current ID of gas to be processed (101=self,102=foreign)
c caFName     = name of file that contains CKD data
      CHARACTER*80 caFName
      INTEGER iGasID

c local variables
      CHARACTER*2 caString,caTemp
      CHARACTER*4 caTemp4
      CHARACTER*80 caDir
      INTEGER iInt,iLenDir,iL

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)
      caDir=kCKDPath          !open AIRS 100 layer reference profile

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
c eg iGasID=101 ==> 'Self', 102 ==> For
      caTemp='    '
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
c      caFName(iLenDir+1:iLenDir+2)=caTemp(1:2)
c      iLenDir=iLenDir+2                           !we have added 2 characters
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

      include 'kcarta.param'

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
        caDir=kNewRefPath       !open kLAYERS created reference profile
      ELSE
        caDir=kRefPath          !open AIRS 100 layer reference profile
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

      WRITE (kStdWarn,2000) caFName
 2000  FORMAT('Reference datafile name is : ',/,A80)

      RETURN
      END

c************************************************************************
c for 1 <= iGasID  <= kGasComp
c thus if we send in (2,430,1) then caFname='q430_g2.dat   '
      SUBROUTINE CompFileName(iGasID,iFileStartFr,iTag,caFName)

      include 'kcarta.param'

c iGasID   = ID of gas to be processed
c iFileStartFr  = start freq of file (integer)
c caFName  = final name of file we wish to read in 
c iTag     = depending on (1,2,3), this is q r or s
      CHARACTER*120 caFName
      INTEGER iGasID,iFileStartFr,iTag
      
c local variables
      CHARACTER*13 caTemp
      CHARACTER*107 caDir
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

      iInt=107
 11      CONTINUE
      IF ((caDir(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN
        iInt=iInt-1
        GO TO 11
        END IF
      iLenDir=iInt

      DO iInt=1,120
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
      ELSE IF (iTag .EQ. 20) THEN 
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
c 2000 FORMAT('Compressed datafile name is : ',/,A120)
      RETURN
      END

c************************************************************************
c this function finds the pressure layer at which rPressStart is within,
c as well as the fraction of the layer that it occupies
      REAL FUNCTION FindBottomTemp(rP,raProfileTemp)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      real rP,raProfileTemp(kProfLayer)
      integer iFound,i1,i2
      real rP1,rP2,T1,T2
      real raP(3),raT(3),Y2A(3),rT

      rT=0.0

      if (rP .ge. plev(1)) then
        !this is WHOLE of the bottom layer
        i1=1
      else if (rP .le. plev(kProfLayer+1)) then
        !this is ludicrous
        write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
        CALL DoStop

      else
        !first find the AIRS layer within which it lies
        iFound=-1 
        i1=1 
        i2=2 
 10          CONTINUE 
        IF ((rP .LE.plev(i1)).AND.(rP .GT.plev(i2))) THEN 
          iFound=1 
          END IF 
        IF ((iFound .LT. 0) .AND. (i1 .LT. kProfLayer)) THEN 
          i1=i1+1 
          i2=i2+1 
          GO TO 10 
          END IF 
        IF ((iFound .LT. 0)) THEN 
          IF (abs(rP-plev(kProfLayer+1)) .LE. delta) THEN 
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

      IF ((i1 .gt. kProfLayer) .or. (i1 .lt. 0)) then
        write(kStdErr,*) 'sorry : cannot find surface temp for '
        write(kStdErr,*) 'layers outside 1 and ',kProfLayer
        call DoStop
        end if
          

      !now find the temperature
      if (i1 .EQ. 1) then          !do linear interp
        rP1=(plev(2)-plev(1))/log(plev(2)/plev(1))
        rP2=(plev(3)-plev(2))/log(plev(3)/plev(2))
        T1=raProfileTemp(1)
        T2=raProfileTemp(2)
        rT=T2-(rP2-rP)*(T2-T1)/(rP2-rP1)
      elseif (i1 .GE. (kProfLayer-1)) then          !do linear interp
        rP1=(plev(kProfLayer)-plev(kProfLayer-1))/
     $        log(plev(kProfLayer)/plev(kProfLayer-1))
        rP2=(plev(kProfLayer+1)-plev(kProfLayer))/
     $        log(plev(kProfLayer+1)/plev(kProfLayer))
        T1=raProfileTemp(kProfLayer-1)
        T2=raProfileTemp(kProfLayer)
        rT=T2-(rP2-rP)*(T2-T1)/(rP2-rP1)
      else          !do spline ... note that the pressures have to 
                    !be in ascending order for good interpolation
        rP1=(plev(i1)-plev(i1-1))/log(plev(i1)/plev(i1-1))
        raP(3)=rP1
        rP1=(plev(i1+1)-plev(i1))/log(plev(i1+1)/plev(i1))
        raP(2)=rP1
        rP1=(plev(i1+2)-plev(i1+1))/log(plev(i1+2)/plev(i1+1))
        raP(1)=rP1
        raT(3)=raProfileTemp(i1-1)
        raT(2)=raProfileTemp(i1)
        raT(1)=raProfileTemp(i1+1)
        call rSplin(raP,raT,Y2A,3,rP,rT)
        end if

      FindBottomTemp=rT

      RETURN
      END
c************************************************************************
