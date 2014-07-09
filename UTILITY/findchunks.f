c compile with f77 -u findchunks.f
c      this file finds which gases have kCOEFF data for the chunks
c      copied from kcartamisc.f

      CHARACTER*80 kCompParamFile,kXsecParamFile
      INTEGER iCount,iGasID,iChunk,iI,iJ,iDummy,iErr
      REAL rL,rH
      INTEGER iaLists(63),iCheckXsecDataBase,iCheckCompDataBase
      INTEGER kXsecFormat

c      PARAMETER (kXsecParamFile = '../DATA/General/xsec107.param')
c      PARAMETER (kCompParamFile = '../DATA/General/comp107.param')
c      PARAMETER (kXsecFormat = 1)      

      PARAMETER (kXsecParamFile = '../DATA/General/xsec.param')
      PARAMETER (kCompParamFile = '../DATA/General/comp97.param')
      PARAMETER (kXsecFormat = -1)      

      iChunk = 1
      rL     = 605.0
      rH     = 630.0

      do iChunk=1,89

        !initialize list to 0
        do iGasID=1,63
          iaLists(iGasID)=0
          end do

        iCount=0

        !check comp files
        do iGasID=1,49
          iDummy=iCheckCompDataBase(kCompParamFile,
     $                                    iGasID,rL,rH,2,iErr)
          if (iDummy .GT. 0) THEN
            iCount=iCount+1
            iaLists(iCount)=iGasID
            end if
          end do

        !check xsec files
        do iGasID=50,63
          iDummy=iCheckXsecDataBase(kXsecParamFile,kXsecFormat,
     $                                    iGasID,rL,rH,2,iErr)
          if (iDummy .GT. 0) THEN
            iCount=iCount+1
            iaLists(iCount)=iGasID
            end if
          end do

          !print out the results
          WRITE(*,103) rL,rH,iCount
 103      FORMAT('start,stop freqs = ',F10.5,'  ',F10.5,' count = ',3I)
          print *, (iaLIsts(iJ),iJ=1,iCount)
          print *,'      '
           
          rL=rL+25.0
          rH=rH+25.0
          end do

        end

                
c************************************************************************
c this function checks to see if the cross section data base exists for
c gasID, between freqs rL,rH
c if rL,rH < 0 then checking basic presence of gas in data base
      INTEGER FUNCTION iCheckXsecDataBase(kXsecParamFile,kXsecFormat,
     $                                    iGasID,rL,rH,iTagIn,iErr)

c iGasID = GAS ID to be searched for in the database
c rL,rH  = freq start/stop points
c iErr   = error status (mainly associated with not finding the relevant file)
c iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
c          (not used if kXsecFormat < 0)
      INTEGER iGasID,iErr,iTagIn
      REAL rL,rH

c local variables
      INTEGER iIOUN,iFileErr,iID,iTag
      INTEGER iLine,iNpts,iTemps
      CHARACTER*80 caLine,caFName,kXsecParamFile
      REAL rLower,rHigher

      INTEGER kTempUnit,kTempUnitOpen,kXsecFormat
      kTempUnit=10

c assume GASID , freqs are wrong
      iCheckXsecDataBase=-1

      caFName=kXsecParamFile
      iIOUN=kTempUnit
      OPEN(UNIT=iIOUN,FILE=caFName,STATUS='old',
     $    FORM='FORMATTED',IOSTAT=iFileErr)

      IF (kXsecFormat .LT. 0) THEN      
ccccccccccc this is the original format : read old style xsec.param file
        IF (iFileErr .NE. 0) THEN
          iErr=0
          WRITE(*,103) iFileErr,caFName
 103      FORMAT('ERROR! number ',I5,' opening xsec database file : 
     $    ',/,A84)
          END IF
        kTempUnitOpen=1

c read file util GASID, freq bounds match found or EOF
 20     READ(iIOUN,5020,END=777) caLine
        READ(caLine,*) iLine,iID,iNpts,iTemps,rLower,rHigher

        IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is the -1 option in XSCGAS
          iCheckXsecDataBase=1
          END IF

        IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $    (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
          iCheckXsecDataBase=1
          END IF

        IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
     $     (rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
          iCheckXsecDataBase=1
          END IF

        IF (iCheckXsecDataBase .LT. 0) THEN
          GOTO 20
          END IF
      
 777    CONTINUE
      ELSE
ccccccccccc this is the new format : read comp.param style xsec.param file
c read file util GASID, freq bounds match found or EOF
 30     READ(iIOUN,5020,END=888) caLine
        READ(caLine,*) iID,rLower,rHigher,iTag

        IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is basically the -1 option in XSCGAS
          iCheckXSecDataBase=1
          END IF

        IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $            (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
c this is when we WANT to UNCOMPRESS files!
          IF (iTag .EQ. iTagIn) THEN   !actually uncompressing stuff
            iCheckXsecDataBase=1
          ELSE                         !actually uncompressing stuff
c           print *,'for GasID ',iGasID,'program says iTag  = ',iTagIN
c           print *,'while comp.param database file says iTag  = ',iTag
c           print *,'going on ... '
            END IF
          END IF

c this next option cannot exist, as during the uncompression we take entire 
c 25 cm-1 chunks at a time
c      IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
c     $(rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
c        iCheckCompDataBase=1
c        END IF

        IF (iCheckXSecDataBase .LT. 0) THEN
          GOTO 30
          END IF
 888    CONTINUE

        END IF

      CLOSE(iIOUN)
      kTempUnitOpen=-1

 5020 FORMAT(A80)
      RETURN
      END
c************************************************************************
c this function checks to see if the comp data file exists for
c gasID, between freqs rL,rH
c if rL,rH < 0 then checking basic presence of gas in data base

c modified 3/30 to include iTag

      INTEGER FUNCTION iCheckCompDataBase(kCompParamFile,
     $                 iGasID,rL,rH,iTagIn,iErr)

c iGasID = GAS ID to be searched for in the database
c rL,rH  = freq start/stop points
c iTagIn = 1,2,3 tells expected wavenumber spacing 0.001 0.0025 0.005 spacing
c iErr   = error status (mainly associated with not finding the relevant file)
      INTEGER iGasID,iErr,iTagIn
      REAL rL,rH

c local variables
      INTEGER iIOUN,iFileErr,iID,iTag
      REAL rLower,rHigher
      CHARACTER*80 caLine,caFname,kCompParamFile

      INTEGER kTempUnit,kTempUnitOpen
      kTempUnit=10

c assume GASID , freqs are wrong
      iCheckCompDataBase=-1

      caFName=kCompParamFile
      iIOUN=kTempUnit
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='old',
     $    FORM='FORMATTED',IOSTAT=iFileErr)

      IF (iFileErr .NE. 0) THEN
        iErr=0
        WRITE(*,103) iFileErr,caFname
 103    FORMAT('ERROR! number ',I5,' opening comp database file : 
     $  ',/,A84)
        END IF
      kTempUnitOpen=1

c read file util GASID, freq bounds match found or EOF
 20   READ(iIOUN,5020,END=777) caLine
 5020 FORMAT(A80)
      READ(caLine,*) iID,rLower,rHigher,iTag

      IF ((iID .EQ. iGasID) .AND. (rL .LT. 0) .AND. (rH .LT. 0)) THEN
c basic presence of gas tested for, and found
c this is basically the -1 option in MOLGAS
        iCheckCompDataBase=1
        END IF

      IF ((iID .EQ. iGasID) .AND. (rL .GE. rLower) .AND. 
     $            (rH .LE. rHigher)) THEN
c presence of gas tested for, with freq bounds, and found well within
c this is when we WANT to UNCOMPRESS files!
        IF (iTag .EQ. iTagIn) THEN   !actually uncompressing stuff
          iCheckCompDataBase=1
        ELSE                         !actually uncompressing stuff
c          print *,'for GasID ',iGasID,'program says iTag  = ',iTagIN
c          print *,'while comp.param database file says iTag  = ',iTag
c          print *,'going on ... '
          END IF
        END IF

c this next option cannot exist, as during the uncompression, we take entire 
c 25 cm-1 chunks at a time
c      IF ((iID .EQ. iGasID) .AND. (rL .LE. rHigher) .AND. 
c     $(rH .GE. rLower)) THEN
c presence of gas tested for, with freq bounds, and found
c        iCheckCompDataBase=1
c        END IF

      IF (iCheckCompDataBase .LT. 0) THEN
        GOTO 20
        END IF
      
 777  CONTINUE
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      RETURN
      END
c************************************************************************
