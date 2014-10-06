c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c note : adjustl is a f77 standard function .. that does not seem to be
c        implemented on all compilers! (eg it does not work on g77)
c so can instead do call adjustleftstr(caIn,caOut)

c************************************************************************ 
      LOGICAL FUNCTION isfinite(a)
      REAL a
      isfinite = (a-a).EQ.0
      END

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
c this subroutine takes an integer (0 <= iInt <= 99) and converts to a string
      SUBROUTINE int2str(iInt,caStrOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*(*) caStrOut
      CHARACTER*2 caStr
      INTEGER iInt
 
      IF (iInt .LT. 0) THEN
        write(kStdErr,*) 'int2str expects iInt >= 0'
        Call Dostop
      END IF
    
      IF (iInt .LT. 10) THEN
        write(caStr,1) iInt
      ELSEIF (iInt .LT. 100) THEN
        write(caStr,2) iInt
c      ELSEIF (iInt .LT. 1000) THEN
c        write(caStr,3) iInt
c      ELSEIF (iInt .LT. 10000) THEN
c        write(caStr,4) iInt
c      ELSEIF (iInt .LT. 100000) THEN
c        write(caStr,5) iInt
      ELSE
        write(kStdErr,*) 'hmm iInt >= 100 int2str unhappy'
        Call DoStop
      END IF

      CALL adjustleftstr(caStr,caStrOut)

 1    FORMAT(I1)
 2    FORMAT(I2)
 3    FORMAT(I3)
 4    FORMAT(I4)
 5    FORMAT(I5)

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
      INTEGER iI,i1,i2,iJ,iLeftjust_lenstr
      CHARACTER*80 caSum

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)

c do some initializations
      i1 = iLeftjust_lenstr(ca1,len(ca1))
      i2 = iLeftjust_lenstr(ca2,len(ca1))

      CALL blankstr(caSum,len(caSum))

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
      SUBROUTINE CKDFileName(caFName,iGasID,iTag,iActualTag)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID      = current ID of gas to be processed (101=self,102=foreign)
c caFName     = name of file that contains CKD data
      CHARACTER*120 caFName
      INTEGER iGasID,iTag,iActualTag

c local variables
      CHARACTER*2 caString,caTemp
      CHARACTER*4 caTemp4
      CHARACTER*120 caDir
      INTEGER iInt,iLenDir,iL,iLeftjust_lenstr
c      CHARACTER*(*) adjustl

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)
      IF (kaTag(iTag) .EQ. 20) THEN
        !! 605-2830 cm-1   'r' prefix
        caDir = kCKDPath          
      ELSEIF (kaTag(iTag) .EQ. 55) THEN
        !! -------> 25000-45000 cm-1  <---------------------   'u' prefix
        caDir = kCKDPathu          
      ELSEIF (kaTag(iTag) .EQ. 50) THEN
        !! -------> 12000-25000 cm-1  <---------------------   'v' prefix
        caDir = kCKDPathv          
      ELSEIF (kaTag(iTag) .EQ. 40) THEN
        !! -------> 8250-12500 cm-1  <---------------------   'o' prefix
        caDir = kCKDPatho          
      ELSEIF (kaTag(iTag) .EQ. 35) THEN
        !! -------> 5500-8200 cm-1  <---------------------   'n' prefix
        caDir = kCKDPathn          
      ELSEIF (kaTag(iTag) .EQ. 30) THEN
        !! -------> 3350-5550 cm-1  <---------------------   'm' prefix
        caDir = kCKDPathm          
      ELSEIF (kaTag(iTag) .EQ. 25) THEN
        !! -------> 2830-3350 cm-1  <---------------------   's' prefix
        caDir = kCKDPaths
      ELSEIF (kaTag(iTag) .EQ. 15) THEN
        !! -------> 500-605 cm-1  <---------------------   'q' prefix
        caDir = kCKDPathq   
      ELSEIF (kaTag(iTag) .EQ. 12) THEN
        !! -------> 300-510 cm-1  <---------------------   'p' prefix
        caDir = kCKDPathp   
      ELSEIF (kaTag(iTag) .EQ. 10) THEN
        !! -------> 140-310 cm-1  <---------------------   'k' prefix
        caDir = kCKDPathk   
      ELSEIF (kaTag(iTag) .EQ. 08) THEN
        !! -------> 080-150 cm-1  <---------------------   'j' prefix
        caDir = kCKDPathj   
      ELSEIF (kaTag(iTag) .EQ. 06) THEN
        !! -------> 050-080 cm-1  <---------------------   'h' prefix
        caDir = kCKDPathh   
      ELSEIF (kaTag(iTag) .EQ. 04) THEN
        !! -------> 030-050 cm-1  <---------------------   'g' prefix
        caDir = kCKDPathg   
      ELSEIF (kaTag(iTag) .EQ. 02) THEN
        !! -------> 020-030 cm-1  <---------------------   'f' prefix
        caDir = kCKDPathf   
      ELSE
        write(kStdErr,*) 'Current iTag = ',iTag
        write(kStdErr,*) 'No compressed datafilename(s) for this band'
        CALL DoStop         
      END IF
        
c do some initializations
      iLenDir = iLeftjust_lenstr(caDir,len(caDir))

      CALL blankstr(caFname,len(caFName))

      caString = '  '
      caTemp   = '00'

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

c CKD         = CKD version 0,1,21,23,24,50+
c now process CKD ver so that we end up with a right padded string
c eg 2 ---> '2 ', 12 ---> '12' etc
      WRITE(caString,15) kCKD
 15   FORMAT(I2)
c this is right justified ... change to left justified
c      caTemp = adjustl(caString)
      CALL adjustleftstr(caString,caTemp)

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
 2000 FORMAT('CKD datafile name is : ',/,A120)

      RETURN
      END

c************************************************************************ 
c this subroutine finds the name of the reference profile for the current gas
      SUBROUTINE FindReferenceName(caFName,iGasID,iLowerOrUpper)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLowerOrUpper  = open lower atm AIRS (100 layer) file (-1)
c                   or  upper atm AIRS (100 layer) file (+1)
c iGasID         = current ID of gas to be processed
c caFName        = name of file that contains reference profile
      CHARACTER*80 caFName
      INTEGER iGasID,iLowerOrUpper

c local variables
      CHARACTER*2 caString,caTemp
      CHARACTER*80 caDir
      INTEGER iInt,iLenDir,iLeftjust_lenstr
c      CHARACTER*(*) adjustl

ccccc user supplied info
c here give the directory path to where the reference profiles are
c note we need all layers read in (kProfLayer >= kMaxLayer)
      IF (iLowerOrUpper .EQ. +1) THEN
        !the CO2 profile in NLTE/UA = caUpperAtmRefPath
        !should be DIFFERENT from that in NLTE/USUALLAYERS
        !so open upper atm AIRS100 layer reference profile 
        caDir = caUpperAtmRefPath  
      ELSEIF (iLowerOrUpper .EQ. -1) THEN
        !open original AIRS100 layer reference profile
        !INTERESTING NOTE : as expected, the CO2 profile in KOrigRefPath
        !should be the same as that in NLTE/USUALLAYERS
        caDir = kOrigRefPath      
      ELSE
        write(kStdErr,*) 'need iLowerOrUpper = -/+ 1 or 2'
        CALL DoStop
      END IF

c do some initializations
      iLenDir = iLeftjust_lenstr(caDir,len(caDir))

      CALL blankstr(caFName,len(caFname))
      CALL blankstr(caString,len(caString))
      CALL blankstr(caTemp,len(caTemp))

c now process iGasID so that we end up with a right padded string
c eg 2 ---> '2 ', 12 ---> '12' etc
 15   FORMAT(I2)
      WRITE(caString,15) iGasID
c this is right justified ... change to left justified
c      caTemp = adjustl(caString) 
      CALL adjustleftstr(caString,caTemp)

c now concatenate all this together in caFname
      caFName(1:iLendir)           = caDir(1:iLenDir)
      caFName(iLenDir+1:iLenDir+6) = 'refgas'
      caFName(iLenDir+7:iLenDir+8) = caTemp(1:2)

      iLenDir = iLeftjust_lenstr(caFName,len(caFname))

      IF (iLowerOrUpper .EQ. +1) THEN
        IF (kCO2ppmv .EQ. 370) THEN
          !!!no need to do nothing, this is the name of the ref profile
          iLendir = iLendir
        ELSEIF (kCO2ppmv .EQ. 378) THEN
          caFName(iLenDir+1:iLenDir+8) = '_378ppmv'
        ELSEIF (kCO2ppmv .EQ. 385) THEN
          caFName(iLenDir+1:iLenDir+8) = '_385ppmv'
        ELSEIF (kCO2ppmv .EQ. 400) THEN
          caFName(iLenDir+1:iLenDir+8) = '_400ppmv'
        ELSE
          write(kStdErr,*) 'in  subroutine FindReferenceName '
          write(kStdErr,*) '?? CO2 ppmv = 370/378/385/400, not ',kCO2ppmv
          CALL DoStop
        END IF
      END IF

      IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $     (abs(kLongOrShort) .LE. 1)) THEN
        IF (iLowerOrUpper .LT. 0) THEN
          WRITE (kStdWarn,2000) caFName
        ELSE
          WRITE (kStdWarn,2010) caFName
        ENDIF
      END IF

 2000 FORMAT('lower atm 100 AIRS layer Reference datafile name is : ',/,A80)
 2010 FORMAT('upper atm 100 AIRS layer Reference datafile name is : ',/,A80)

      RETURN
      END

c************************************************************************
c for 1 <= iGasID  <= kGasComp
c so if we send in (+1,2,430,1) then caFname='q430_g2.dat   '
c also has to add on the DIRECTORY
c only worry about linemixing/COUSIN/NLTE for CO2 in the 605-2830 cm-1 tag 
c (iTag=2,prefix='r')
      SUBROUTINE CompFileName(iLineMix,iGasID,rFileStartFr,
     $                        iTag,iActualTag,caFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = ID of gas to be processed
c rFileStartFr  = start freq of file 
c caFName  = final name of file we wish to read in 
c iTag     = depending on (1,2,3), this is q r or s
c
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c iLineMix = for CO2, can do for 1100 - 0.005 mb : 
c                            iLineMix = +1   linemixing files, or
c                            iLineMix = -1   Cousin files, or
c                            iLineMix = -2   weak backgnd 4 um CO2 files
c                     can do for 0.005 - 0.00025 mb
c                            iLineMix = -3   weak backgnd 4 um CO2 files
c NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
c   for CO2 4 um NLTE (1100-0.005 mb) 
c                            iLineMix = +100+x  ODs          as fcn of solzen 
c                            iLineMix = +200+x  PlanckCoeffs as fcn of solzen
c   for CO2 4 um NLTE (0.005-0.0005 mb) 
c                            iLineMix = +300+x  ODs          as fcn of solzen 
c                            iLineMix = +400+x  PlanckCoeffs as fcn of solzen
c where the +x = +0,40,60,80,85,90 and depends on solzen
c NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW NEW 
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c            for all other gases, iLineMix = irrelevant (set to +1 usually)
c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      CHARACTER*120 caFName
      INTEGER iGasID,iTag,iActualTag,iLineMix
      REAL rFileStartFr
      
c local variables
      CHARACTER*16 caTemp
      CHARACTER*120 caDir,caDirUse
      CHARACTER*2 caString2,caTemp2
      CHARACTER*5 caString3,caTemp3
      CHARACTER*5 caString5,caTemp5
      CHARACTER*9 caTemp9
      CHARACTER*16 caTemp16
      CHARACTER*2  ca2
      INTEGER iInt,iLen2,iLen3,iLen5,iLen16,iLenDir,iLenX,iLeftjust_lenstr
c      CHARACTER*(*) adjustl  !not implemented in g77!!!
      
ccccc user supplied info
c here give the directory path to where the compressed files are ...
c keeping in mind that the water files may be in a directory by themselves
      IF (kaTag(iTag) .EQ. 55) THEN
        !! -------> 25000-45000 cm-1  <---------------------   'u' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathu
        ELSE
          caDir = kCompPathu
        END IF
      ELSEIF (kaTag(iTag) .EQ. 50) THEN
        !! -------> 12000-25000 cm-1  <---------------------   'v' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathv
        ELSE
          caDir = kCompPathv
        END IF
      ELSEIF (kaTag(iTag) .EQ. 40) THEN
        !! -------> 8250-12250 cm-1  <---------------------   'o' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPatho
        ELSE
          caDir = kCompPatho
        END IF
      ELSEIF (kaTag(iTag) .EQ. 35) THEN
        !! -------> 5550-8200 cm-1  <---------------------   'n' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathn
        ELSE
          caDir = kCompPathn
        END IF
      ELSEIF (kaTag(iTag) .EQ. 30) THEN
        !! -------> 3550-5550 cm-1  <---------------------   'm' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathm
        ELSE
          caDir = kCompPathm
        END IF
      ELSEIF (kaTag(iTag) .EQ. 25) THEN
        !! -------> 2830-3580 cm-1  <---------------------   's' prefix
        IF ((iGasID .EQ. 1) .AND. (rFileSTartFr+1 .GT. kWaterIsobandStop2))THEN
          caDir = kWaterPaths
        ELSEIF ((iGasID.EQ.1).AND.(rFileSTartFr+1.LE.kWaterIsobandStop2)) THEN
          caDir = kWaterIsotopePath
        ELSEIF (iGasID .EQ. 103) THEN
          !!heavy water : isotope 4; iDoAdd = -1 unless
          !!  kWaterIsoBandStart2 <= f < kWaterIsoBandStop2
          caDir = kWaterIsotopePath
        ELSE
          caDir = kCompPaths
        END IF
      ELSEIF (kaTag(iTag) .EQ. 15) THEN
        !! -------> 500-605 cm-1  <---------------------   'q' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathq
        ELSE
          caDir = kCompPathq
        END IF
      ELSEIF (kaTag(iTag) .EQ. 12) THEN
        !! -------> 300-510 cm-1  <---------------------   'p' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathp
        ELSE
          caDir = kCompPathp
        END IF
      ELSEIF (kaTag(iTag) .EQ. 10) THEN
        !! -------> 140-310 cm-1  <---------------------   'k' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathk
        ELSE
          caDir = kCompPathk
        END IF
      ELSEIF (kaTag(iTag) .EQ. 08) THEN
        !! -------> 080-140 cm-1  <---------------------   'j' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathj
        ELSE
          caDir = kCompPathj
        END IF
      ELSEIF (kaTag(iTag) .EQ. 06) THEN
        !! -------> 050-080 cm-1  <---------------------   'h' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathh
        ELSE
          caDir = kCompPathh
        END IF
      ELSEIF (kaTag(iTag) .EQ. 04) THEN
        !! -------> 030-050 cm-1  <---------------------   'g' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathg
        ELSE
          caDir = kCompPathg
        END IF
      ELSEIF (kaTag(iTag) .EQ. 02) THEN
        !! -------> 020-030 cm-1  <---------------------   'f' prefix
        IF (iGasID .EQ. 1) THEN
          caDir = kWaterPathf
        ELSE
          caDir = kCompPathf
        END IF
      !!!!--->>> this is default kCARTA database from 605-2830 cm-1 <<<--!!!!
      ELSEIF (kaTag(iTag) .EQ. 20) THEN
        !! -------> 0605-2830 cm-1  <---------------------   'r' prefix
        IF (iGasID .EQ. 1) THEN
          IF (rFileStartFr+1.0 .LT. kWaterIsobandStart1) THEN
            !! usual water : all isotopes 1 2 3 4 5 6   f< 1105 cm-1
            caDir = kWaterPath
          ELSEIF ((rFileStartFr+1. GE. kWaterIsobandStart1) .AND. 
     $            (rFileStartFr+1. LE. kWaterIsobandStop1)) THEN
            !!special water : isotopes 1 2 3   5 6   1105 < f < 1705 cm-1
            caDir = kWaterIsotopePath
          ELSEIF ((rFileStartFr+1. GE. kWaterIsobandStop1) .AND. 
     $            (rFileStartFr+1. LE. kWaterIsobandStart2)) THEN
            !! usual water : all isotopes 1 2 3 4 5 6   1730 < f <= 2405 cm-1
            caDir = kWaterPath
          ELSEIF ((rFileStartFr+1. GE. kWaterIsobandStart2) .AND. 
     $            (rFileStartFr+1. LE. kWaterIsobandStop2)) THEN
            !!special water : isotopes 1 2 3   5 6   2405 < f < 3305 cm-1
            caDir = kWaterIsotopePath
          END IF
        ELSEIF (iGasID.EQ.103) THEN
          IF ((rFileStartFr+1. GE. kWaterIsobandStart1) .AND. 
     $        (rFileStartFr+1. LE. kWaterIsobandStop1)) THEN
            !!heavy water : isotope 4
            caDir = kWaterIsotopePath
          ELSEIF ((rFileStartFr+1. GE. kWaterIsobandStart2) .AND. 
     $          (rFileStartFr+1. LE. kWaterIsobandStop2)) THEN
            !!heavy water : isotope 4
            caDir = kWaterIsotopePath
        END IF
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
          ELSEIF ((iLineMix .GE. 100) .AND. (iLineMix .LE. 290)) THEN
            !!use compressed tables for NLTE in LA
            caDir = kCompressedNLTE_LA
          ELSEIF ((iLineMix .GE. 300) .AND. (iLineMix .LE. 490)) THEN
            !!use compressed tables for NLTE in UA
            caDir = kCompressedNLTE_UA
          ELSE
            write(kStdErr,*) 'Wrong option for CO2 in CompFileName',iLineMix
            CALL DoStop
          END IF
        ELSE
          caDir = kCompPath
        END IF            
      ELSE
        write(kStdErr,*) 'iTag,iGasID,rFileStartFr = ',iTag,iGasID,rFileStartFr
        write(kStdErr,*) 'have only coded up iTag for finite number of bands'
        write(kStdErr,*) 'No compressed datafilename(s) for this band'
        CALL DoStop         
      END IF

      IF (kAltDir .EQ. +1) caDir = kcaAltDir   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      iLenDir = iLeftjust_lenstr(caDir,len(caDir))

      CALL blankstr(caFName,len(caFName))
      CALL blankstr(caTemp,len(caTemp))

      CALL blankstr(caString2,len(caString2))
      CALL blankstr(caTemp2,len(caTemp2))

      CALL blankstr(caString3,len(caString3))
      CALL blankstr(caTemp3,len(caTemp3))

      CALL blankstr(caString5,len(caString5))
      CALL blankstr(caTemp5,len(caTemp5))

      CALL blankstr(caTemp9,len(caTemp9))
      CALL blankstr(caTemp16,len(caTemp16))
 
c now process rFileStartFr so that we end up with a right padded string 
c eg 605 ---> '605 ', 1200 ---> '1200' etc
 12   FORMAT(I5)
 13   FORMAT(F5.1)
      IF (abs(rFileStartFr - anint(rFileStartFr)) .LE. 0.001) THEN
        WRITE(caString5,12) nint(rFileStartFr)
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      ELSE
        WRITE(caString5,13) rFileStartFr
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      END IF

c now process iGasID so that we end up with a right padded string to which
c we append '.dat'
c eg 2 ---> '2.dat ', 12 ---> '12.dat' etc
 15   FORMAT(I3)
      WRITE(caString3,15) iGasID
c this is right justified ... change to left justified
c      caTemp3 = adjustl(caString3)
      CALL adjustleftstr(caString3,caTemp3)
      iLen3 = iLeftjust_lenstr(caTemp3,3)
      !gasIDs range from 1-->103 = 3 char

c this is new : if gases are [1],[3..103], then we don;t have to 
c worry about NLTE
      IF (iLineMix .LT. 100) THEN
        ca2 = 'XY'    !! who cares
c else we do need to worry about NLTE, which depends on solzen
c as given by 00,40,60,80,85,90 
      ELSEIF ((iLineMix .EQ. 100) .OR. (iLineMix .EQ. 200) .OR. 
     $        (iLineMix .EQ. 300) .OR. (iLineMix .EQ. 400)) THEN
        ca2 = '00'    !! solzen = 00
      ELSEIF ((iLineMix .EQ. 140) .OR. (iLineMix .EQ. 240) .OR. 
     $        (iLineMix .EQ. 340) .OR. (iLineMix .EQ. 440)) THEN
        ca2 = '40'    !! solzen = 40
      ELSEIF ((iLineMix .EQ. 160) .OR. (iLineMix .EQ. 260) .OR. 
     $        (iLineMix .EQ. 360) .OR. (iLineMix .EQ. 460)) THEN
        ca2 = '60'    !! solzen = 60
      ELSEIF ((iLineMix .EQ. 180) .OR. (iLineMix .EQ. 280) .OR. 
     $        (iLineMix .EQ. 380) .OR. (iLineMix .EQ. 480)) THEN
        ca2 = '80'    !! solzen = 80
      ELSEIF ((iLineMix .EQ. 185) .OR. (iLineMix .EQ. 285) .OR. 
     $        (iLineMix .EQ. 385) .OR. (iLineMix .EQ. 485)) THEN
        ca2 = '85'    !! solzen = 85
      ELSEIF ((iLineMix .EQ. 190) .OR. (iLineMix .EQ. 290) .OR. 
     $        (iLineMix .EQ. 390) .OR. (iLineMix .EQ. 490)) THEN
        ca2 = '90'    !! solzen = 90
      END IF

      IF (iLineMix .LT. 100) THEN
c no need to worry : only LTE gases in any band
c filenames eg r2205_g3.dat means "r" prefix, 2205 cm-1 chunk, g3
        caTemp9 = '_g' // caTemp3(1:iLen3) // '.dat'  

c else we need to worry about NLTE ODs and Planck coeffs in 
c lower or upper atm
c filenames eg r2205n60g2.dat means 
c    "r" prefix, 2205 cm-1 chunk, solzen 60, n = NLTE LA ods
      ELSEIF ((iLineMix .GE. 100) .AND. (iLineMix .LE. 190)) THEN
        !NLTE LA ods
        caTemp9 = 'n' // ca2 // 'g' // caTemp3(1:iLen3) // '.dat'  
      ELSEIF ((iLineMix .GE. 200) .AND. (iLineMix .LE. 290)) THEN
        !NLTE LA planck coeffs
        caTemp9 = 'p' // ca2 // 'g' // caTemp3(1:iLen3) // '.dat'  
      ELSEIF ((iLineMix .GE. 300) .AND. (iLineMix .LE. 390)) THEN
        !NLTE LA ods
        caTemp9 = 'N' // ca2 // 'g' // caTemp3(1:iLen3) // '.dat'  
      ELSEIF ((iLineMix .GE. 400) .AND. (iLineMix .LE. 490)) THEN
        !NLTE LA planck coeffs
        caTemp9 = 'P' // ca2 // 'g' // caTemp3(1:iLen3) // '.dat'  
      END IF

c using caTemp, put together the compressed data file name, starting with 'r'
      iLenX = 1
 100  CONTINUE
      IF (kaTag(iLenX) .EQ. iActualTag) THEN
        GOTO 101
      ELSEIF (iLenX .EQ. kW) THEN
        write (kStdErr,*) 'Ooops cannot match iTag with kaTag'
        CALL DoStop
      ELSEIF (iLenX .LT. kW) THEN
        iLenX = iLenX + 1
        GOTO 100         
      END IF
 101  CONTINUE
      caTemp(1:1) = kaCtag(iLenX)
      
c followed by rFileStartFr
      caTemp(2:6)=caTemp5(1:5)

      iLenX = iLeftjust_lenstr(caTemp,len(caTemp)) + 1
      caTemp(iLenX:iLenX+9-1)=caTemp9(1:9)
      !!now concatenate directory path and filename : caDir+caTemp
      IF (iLenDir .GT. 0) THEN
        caFName(1:iLenDir)=caDir(1:iLenDir)
      END IF
      caFName(iLenDir+1:iLenDir+1+15)=caTemp(1:14)      

c      WRITE (kStdWarn,2000) caFName
c 2000 FORMAT('Compressed datafile name is : ',/,A120)      

      RETURN
      END

c************************************************************************
c for 1 <= iGasID  <= kGasComp
c thus if we send in (+1,2,430,1) then caFname='q430_g2.dat   '
c also has to add on the DIRECTORY
c only worry about linemixing/COUSIN/NLTE for CO2 in the 605-2830 cm-1 tag 
c (iTag=2,prefix='r')
      SUBROUTINE CompFileNameCKD(iLineMix,iGasID,rFileStartFr,
     $                        iTag,iActualTag,caFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iGasID   = ID of gas to be processed
c rFileStartFr  = start freq of file 
c caFName  = final name of file we wish to read in 
c iTag     = depending on (1,2,3), this is q r or s
c iLineMix = for CO2, can do for 1100 - 0.005 mb : 
c                            iLineMix = +1 for the linemixing files, or
c                            iLineMix = -1 for the Cousin files, or
c                            iLineMix = -2 for the weak backgnd 4 um CO2 files
c                     can do for 0.005 - 0.00025 mb
c                            iLineMix = -3 for the weak backgnd 4 um CO2 files
c            for all other gases, iLineMix = irrelevant (set to +1 usually)

      CHARACTER*120 caFName
      INTEGER iGasID,iTag,iActualTag,iLineMix
      REAL rFileStartFr
      
c local variables
      CHARACTER*26 caTemp
      CHARACTER*120 caDir,caDirUse
      CHARACTER*2 caString2,caTemp2
      CHARACTER*5 caString3,caTemp3
      CHARACTER*5 caString5,caTemp5
      CHARACTER*9 caTemp9
      CHARACTER*16 caTemp16
      INTEGER iInt,iLen2,iLen3,iLen5,iLen16,iLenDir,iLenX,iLeftjust_lenstr
c      CHARACTER*(*) adjustl  !not implemented in g77!!!
      
ccccc user supplied info
c here give the directory path to where the compressed files are ...
c keeping in mind that the water files may be in a directory by themselves
      IF ((iGasID .EQ. kNewGasLo) .OR. (iGasID .EQ. kNewGasHi)) THEN
        caDir = kCKD_Compr_Path
      ELSE 
        write(kStdErr,*) 'this routine is only for gases 101,102'
        write(kStdErr,*) 'iTag,iGasID,rFileStartFr = ',iTag,iGasID,rFileStartFr
        CALL DoStop         
      END IF

      iLenDir = iLeftjust_lenstr(caDir,len(caDir))

      CALL blankstr(caFName,len(caFName))
      CALL blankstr(caTemp,len(caTemp))

      CALL blankstr(caString2,len(caString2))
      CALL blankstr(caTemp2,len(caTemp2))

      CALL blankstr(caString3,len(caString3))
      CALL blankstr(caTemp3,len(caTemp3))

      CALL blankstr(caString5,len(caString5))
      CALL blankstr(caTemp5,len(caTemp5))

      CALL blankstr(caTemp9,len(caTemp9))
      CALL blankstr(caTemp16,len(caTemp16))
 
c now process rFileStartFr so that we end up with a right padded string 
c eg 605 ---> '605 ', 1200 ---> '1200' etc
 12   FORMAT(I5)
 13   FORMAT(F5.1)
      IF (abs(rFileStartFr - anint(rFileStartFr)) .LE. 0.001) THEN
        WRITE(caString5,12) nint(rFileStartFr)
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      ELSE
        WRITE(caString5,13) rFileStartFr
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      END IF

c now process CKD so that we end up with a right padded string 
c eg 1 ---> '1 ', 24 ---> '24' etc
 112   FORMAT(I2)
      WRITE(caString2,112) kCKD
c this is right justified ... change to left justified
c      caTemp2 = adjustl(caString2)
        CALL adjustleftstr(caString2,caTemp2)
      iLen2   = iLeftjust_lenstr(caTemp2,2)

c now process iGasID so that we end up with a right padded string to which
c we append '.dat'
c eg 2 ---> '2.dat ', 12 ---> '12.dat' etc
 15   FORMAT(I3)
      WRITE(caString3,15) iGasID
c this is right justified ... change to left justified
c      caTemp3 = adjustl(caString3)
      CALL adjustleftstr(caString3,caTemp3)
      iLen3 = iLeftjust_lenstr(caTemp3,3)
      !gasIDs from 1-->103 = 3 char. CKD = 1:99 so at most 2 char
      caTemp16 = '_g' // caTemp3(1:iLen3) // '_CKD_' // caTemp2(1:iLen2) // '.dat'  
      iLen16   = iLeftjust_lenstr(caTemp16,16)

c using caTemp, put together the compressed data file name, starting with 'r'
      iLenX = 1
 100  CONTINUE
      IF (kaTag(iLenX) .EQ. iActualTag) THEN
        GOTO 101
      ELSEIF (iLenX .EQ. kW) THEN
        write (kStdErr,*) 'Ooops cannot match iTag with kaTag'
        CALL DoStop
      ELSEIF (iLenX .LT. kW) THEN
        iLenX = iLenX + 1
        GOTO 100         
      END IF
 101  CONTINUE
      caTemp(1:1) = kaCtag(iLenX)
      
c followed by rFileStartFr
      caTemp(2:6)=caTemp5(1:5)

      iLenX = iLeftjust_lenstr(caTemp,len(caTemp)) + 1
      caTemp(iLenX:iLenX+iLen16-1)=caTemp16(1:iLen16)

      !!now concat the CKD number in front of this
      caTemp = '/' // caTemp2(1:iLen2) // '/' // caTemp

      !!now concatenate directory path and filename : caDir+caTemp
      IF (iLenDir .GT. 0) THEN
        caFName(1:iLenDir) = caDir(1:iLenDir)
      END IF
      iLenX  = iLeftjust_lenstr(caTemp,len(caTemp))
      caFName(iLenDir+1:iLenDir+1+iLenX)=caTemp(1:iLenX)      
c      print *,'A : ',caFName

c      WRITE (kStdWarn,2000) caFName
c 2000 FORMAT('Compressed datafile name is : ',/,A120)      

      RETURN
      END

c************************************************************************
c this finds the overall CHI function file name
      SUBROUTINE FindChiFileName(caFile)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*120 caFile    !!!input and output parameter

      CHARACTER*120 caTemp,caDir
      INTEGER iI,iJ,iLeftjust_lenstr

      CALL BlankStr(caTemp,len(caTemp))

      caDir = kChiFile

      iI = iLeftjust_lenstr(caDir,len(caDir))
      caTemp(1:iI) = caDir(1:iI)

      iJ = iLeftjust_lenstr(caFile,len(caFile))
      caTemp(iI+1:iI+1+iJ) = caFile(1:iJ)

      DO iI = 1,120
        caFile(iI:iI) = caTemp(iI:iI)
      END DO

c      write (kStdWarn,30) caFile
c 30   FORMAT('chifile = ',A120)

      RETURN
      END

c************************************************************************
c this subroutine gets the solar file name
      SUBROUTINE GetSolarFileName(fname,rFileStartFr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 fname
      REAL rFileStartFr

      CHARACTER*5 caString5,caTemp5
      INTEGER iInt,iLen5,iLenDir,iLeftjust_lenstr
c      CHARACTER*(*) adjustl

      CALL blankstr(fname,len(fname))

      fname = kSolarPath

      iLenDir = iLeftjust_lenstr(fname,len(fname))
      fname(iLenDir+1:iLenDir+9)='rad_solar'
      iLenDir = iLeftjust_lenstr(fname,len(fname))

c now process rFileStartFr so that we end up with a right padded string  
c eg 605 ---> '605 ', 1200 ---> '1200', 10500 --> '10500' etc 
      CALL blankstr(caString5,len(caString5))
      CALL blankstr(caTemp5,len(caTemp5))

 12   FORMAT(I5)
 13   FORMAT(F5.1)
      IF (abs(rFileStartFr - anint(rFileStartFr)) .LE. 0.001) THEN
        WRITE(caString5,12) nint(rFileStartFr)
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      ELSE
        WRITE(caString5,13) rFileStartFr
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      END IF

c now concat the 3 parts together fname+caTemp5+'.dat'
      fname(iLenDir+1:iLenDir+1+iLen5)=caTemp5(1:iLen5)
      iLenDir=iLenDir+iLen5
      fname(iLenDir+1:iLenDir+1+4) = '.dat'

      RETURN
      END

c************************************************************************
c this appends the start wavenumber to caWeak
c basically the same as GetSolarFileName
      SUBROUTINE Get_WeakOptDepthFileName(iGasID,caWeak,raWaves,caWeakFr)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

      CHARACTER*80 caWeak,caWeakFr
      REAL raWaves(kMaxPts)
      INTEGER iGasID

c local variables
      INTEGER iFloor,iLeftjust_lenstr
      REAL rFileStartFr
      CHARACTER*2 caString2,caTemp2
      CHARACTER*5 caString5,caTemp5
      INTEGER iInt,iLen2,iLen5,iLenDir
c      CHARACTER*(*) adjustl

      rFileStartFr = raWaves(1)

      CALL BlankStr(caWeakFr,len(caWeakFr))
      caWeakFr = caWeak

      iLenDir = iLeftjust_lenstr(caWeakFr,len(caWeakFr))
      caWeakFr(iLenDir+1:iLenDir+3) = 'gas'
      iLenDir = iLenDir + 3
 
c now process iGasID so that we end up with a right padded string  
c eg 605 ---> '605 ', 1200 ---> '1200' etc 
      CALL BlankStr(caString2,len(caString2))
      CALL BlankStr(caTemp2,len(caTemp2))

      WRITE(caString2,11) iGasID
 11   FORMAT(I2) 
c this is right justified ... change to left justified 
c      caTemp2 = adjustl(caString2)
      CALL adjustleftstr(caString2,caTemp2)
      iLen2   = iLeftjust_lenstr(caTemp2,len(caTemp2))
      caWeakFr(iLenDir+1:iLenDir+1+iLen2)=caTemp2(1:iLen2)

      iLenDir = iLeftjust_lenstr(caWeakFr,len(caWeakFr))
      caWeakFr(iLenDir+1:iLenDir+1) ='_'
      iLenDir = iLenDir + 1

c now process rFileStartFr so that we end up with a right padded string  
c eg 605 ---> '605 ', 1200 ---> '1200' etc 
 12   FORMAT(I5)
 13   FORMAT(F5.1)
      IF (abs(rFileStartFr - anint(rFileStartFr)) .LE. 0.001) THEN
        WRITE(caString5,12) nint(rFileStartFr)
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      ELSE
        WRITE(caString5,13) rFileStartFr
c this is right justified ... change to left justified
c        caTemp5 = adjustl(caString5)
        CALL adjustleftstr(caString5,caTemp5)
        iLen5   = iLeftjust_lenstr(caTemp5,5)
      END IF

c now concat the 3 parts together caWeakFr+caTemp4+caString4
      caWeakFr(iLenDir+1:iLenDir+1+iLen5) = caTemp5(1:iLen5)
      iLenDir = iLenDir+iLen5
      caWeakFr(iLenDir+1:iLenDir+1+4) = '.dat'

      RETURN
      END

c************************************************************************
c this subtroutine finds the mystery absorption factor
c used by Co2_4um_AERI_LLS_fudge
c WARNING : for iWhichGas = 0, raMystery is used additively (so no effect = 0)
c         : for iWhichGas=1,2, raMystery is used factorly (so no effect = 1)

      SUBROUTINE FindMysteryFactor(raFreq,raMystery,iWhichGas)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input 
      REAL raFreq(kMaxPts)
      INTEGER iWhichGas   ! 0 ==> mystery, +1 = water, +2 = CO2 etc
c output
      REAL raMystery(kMaxPts)

c local vars
      INTEGER iNumMystery,iIOUN,iErr,iFr
      REAL raXadj(kMaxPts),raYadj(kMaxPts)
      CHARACTER*80 caMystery

      IF (iWhichGas .EQ. 0) THEN
        caMystery = '/asl/data/kcarta/KCARTADATA/General/mystery.dat'
      ELSEIF (iWhichGas .EQ. 1) THEN
        caMystery = '/asl/data/kcarta/KCARTADATA/General/aerifudge_gas1.dat'
      ELSEIF (iWhichGas .EQ. 2) THEN
        caMystery = '/asl/data/kcarta/KCARTADATA/General/aerifudge_gas2.dat'
      ELSE
        write(kStdErr,*) 'Invalid code for iWhichGas in FindMysteryFactor'
        CALL DoStop
      END IF

 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80) 
      iIOUN = kTempUnit 
      OPEN(UNIT=iIOUN,FILE=caMystery,STATUS='OLD',FORM='FORMATTED', 
     $    IOSTAT=IERR) 
        IF (IERR .NE. 0) THEN 
          WRITE(kStdErr,*) 'In subroutine FindMysteryFactor' 
          WRITE(kStdErr,1010) IERR, caMystery
          CALL DoSTOP 
          ENDIF 

      kTempUnitOpen = 1 
      iFr = 0
 1020 CONTINUE
      iFr = iFr + 1
      READ(iIOUN,*,END=1030) raXadj(iFr),raYadj(iFr)
      GOTO 1020

 1030 CONTINUE
      CLOSE(iIOUN) 
      kTempUnitOpen = -1 
 
      Call rlinear(raXadj,raYadj,iFr,raFreq,raMystery,kMaxPts)

c      DO iFr = 1,kMaxPts
c        print *,raFreq(iFr),raMystery(iFr)
c      END DO

      RETURN
      END

c************************************************************************
c this funtion returns the "left justified length" of a string by 
c looking for last blank char
c assumes you've done all the work beforehand
c         and max number of characters = iSizeStr
c this is pretty much the same idea as f95 lnblank (standard fcn)
      INTEGER FUNCTION iLeftjust_lenstr(caStr,iSizeStr)

      INTEGER iL,iSizeStr
      CHARACTER*(*) caStr

      iL = iSizeStr
 5    CONTINUE
      IF ((caStr(iL:iL) .EQ. ' ') .AND. (iL .GT. 1)) THEN
        iL = iL - 1
        GOTO 5
      END IF

 20   iLeftjust_lenstr = iL

      RETURN
      END
c************************************************************************
c this subroutine blanks a string
      SUBROUTINE blankstr(caStr,iSizeStr)

      INTEGER iSizeStr
      CHARACTER*(*) caStr

      INTEGER iL

      DO iL = 1,iSizeStr
        caStr(iL:iL) = ' '
      END DO

      RETURN
      END
c************************************************************************
c this subroutine Adjusts string to the left, removing leading blanks and inserting trailing blanks.
      SUBROUTINE adjustleftstr(caStr,caStrOut)
      !! call adjustleftstr(caString5,caTemp5) should be equivalent to ABSOFT caTemp5 = adjustl(caString5) 

      include '../INCLUDE/kcarta.param'       

      CHARACTER*(*) caStr
      CHARACTER*(*) caStrOut

      INTEGER iL,i1,i2,iSizeStr

      iSizeStr = len(caStr)

      IF (len(caStr) .NE. len(caStrOut)) THEN
        write(kStdErr,*) 'hmm, len(caStr) .NE. len(caStrOut)'
        CALL DoStop
      END IF

      DO iL = 1,iSizeStr
        caStrOut(iL:iL) = ' '
      END DO

      i1 = 1
      i2 = 0
 10   CONTINUE
      IF ((caStr(i1:i1) .EQ. ' ') .AND. (i1 .LT. iSizeStr)) THEN
        i2 = i1
        i1 = i1 + 1
        GOTO 10
      END IF

      IF (i2 .EQ. iSizeStr) THEN
        !! all are blank, yaya!!!!! no need to do anything
        i1 = i2
      ELSEIF (i2 .EQ. 0) THEN
        !!nothing is blank, just copy entire string
        DO iL = 1,iSizeStr
          caStrOut(iL:iL) = caStr(iL:iL)
        END DO
      ELSE
        i2 = i2 + 1
        !!first few are blak, copy starting from i2+1
        DO iL = i2,iSizeStr
          caStrOut(iL-i2+1:iL-i2+1) = caStr(iL:iL)
        END DO
      END IF
         
      
      RETURN
      END


c************************************************************************
