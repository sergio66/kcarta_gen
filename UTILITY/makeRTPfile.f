c      this takes in a TestPointProf text file and changes it to a RTP file
c      Sergio De Souza-Machado
c      very closely linked to the ftest1.f file that Howard Motteler wrote

        program makeRTPfile
        implicit none

c RTP declarations
        include 'rtpdefs.f'
        integer rtpopen, rtpread, rtpwrite, rtpclose
        record /RTPHEAD/ head
        record /RTPPROF/ prof
        record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)

c other variables
        integer status, rchan, iPreSet
        character*80 fname, mode

        iPreset = -1
        if (iPreset .gt. 0) then
c set only one header attribute
          hatt(1).fname = 'glist'//char(0)
          hatt(1).aname = 'note'//char(0)
          hatt(1).atext = 'translated from text sonde profiles'//char(0)
c set some sample profile attributes
          patt(1).fname = 'ptemp'//char(0)
          patt(1).aname = 'units'//char(0)
          patt(1).atext = 'degrees K'//char(0)
        else
          CALL RTPINIT(HEAD,PROF) 
          endif

c read in the point profile
        CALL ReadProfile(prof,head)
        
c create the file, write the attributes, write the (mostly unfilled) 
c profile structure, and close the file
c
        mode = 'c'
        print *,'Writing output in running dir ..'
        print *,'output name = makeRTPfile.ip.rtp'
        fname = 'makeRTPfile.ip.rtp'

        status = rtpopen(fname, mode, head, hatt, patt, rchan)
        print *, 'open status,rchan = ', status,rchan
        status = rtpwrite(rchan, prof)
        !print *, 'write status = ', status
        status = rtpclose(rchan)
        !print *, 'close status = ', status

        stop
        end

c************************************************************************
c this subroutine prompts the user for the name of the point profile,
c opens it, reads in the info and saves it in the "prof" record

c the info in the levels profile is as follows 
c !!! first few lines are comments
c cType,mlevs,ngas    where cType = 'p' or 'P'
c gasID list
c lat lon
c rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)

      SUBROUTINE ReadProfile(prof,head)

      implicit none

      include 'rtpdefs.f'

      record /RTPHEAD/ head
      record /RTPPROF/ prof

      INTEGER ngas,mlevs,glist(MaxGas)
      REAL pmin,pmax
      CHARACTER*80 caFName
      CHARACTER*280 caStr
      INTEGER iIOUN2,iErrIO
      INTEGER iI,iZ
      REAL rLat,rLon
      REAL rPress,rTemp,rDensity,raGasAmts(MAXGAS)
      CHARACTER cType

      print *,'Enter file name to process : '
      read *,caFName

      pmin = 1.0e30
      pmax = -1.0e30

      iIOUN2 = 15

c now loop iNpath/iNumGases  times for each gas in the user supplied profile 
      print *,caFName

      OPEN(UNIT=iIOun2,FILE=caFName,STATUS='OLD',FORM='FORMATTED', 
     $    IOSTAT=iErrIO) 
      IF (iErrIO .NE. 0) THEN 
          WRITE(*,1070) iErrIO
        ENDIF 
 1070 FORMAT('Error number ',I5,' opening your profile')

      !!!!!!!! go thru the comments at the top of the file
 30   CONTINUE 
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      IF (caStr(1:1) .EQ. '!') THEN 
        GO TO 30 
        END IF

      !!read type of sonde profile (pressure or height), num levels, num gases
      READ (caStr,*,ERR=13,END=13) cType,mlevs,ngas 
      IF ((cType .NE. 'p') .AND. (cType .NE. 'P'))  THEN
        print *,'Error : makeRTPfile assumes pressure info in the file'
        CLOSE(iIOUN2)
        STOP
        ENDIF

      IF (ngas .GT. MAXGAS) THEN
        print *,'Error : makeRTPfile can only handle ',MAXGAS,' gases'
        CLOSE(iIOUN2)
        STOP
        ENDIF

      IF (mlevs .GT. MAXLEV) THEN
        print *,'Error : makeRTPfile can only handle ',MAXLEV,' levels'
        CLOSE(iIOUN2)
        STOP
        ENDIF

      !! read the gasID list
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) (glist(iI),iI = 1,ngas)
 
      !! read the latitude and longitude
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) rLat,rLon

      prof.plat = rLat
      prof.nlevs = mlevs
      prof.spres = 1100.0      !!start from BOTTOM

      print *,'the gasIDs are : '
      print *, (glist(iI),iI = 1,ngas)
      print *, 'rLat,rLon = ', rLat,rLon
      print *, 'number of levels, number of gases =', mlevs,ngas
      print *, 'WARNING : surf pres temporarily set at ', prof.spres

      !!read the profile info
      DO iZ = 1,mlevs
        READ (iIOUN2,5030,ERR=13,END=13) caStr 
        READ (caStr,*,ERR=13,END=13) rPress,rTemp,rDensity,
     $                               (raGasAmts(iI),iI = 1,ngas)
        !print *,rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)        
        prof.plevs(iZ) = rPress
        prof.pTemp(iZ) = rTemp
        DO iI = 1,ngas
          prof.gamnt(iZ,iI) = raGasAmts(iI)
          END DO
        IF (rPress .LT. pmin) pmin = rPress
        IF (rPress .GT. pmax) pmax = rPress

        ENDDO

 13   CONTINUE
      CLOSE(iIOUN2)
 5030 FORMAT(A280)

c set HEAD values
      head.ptype = LEVPRO
      head.pfields = PROFBIT
c      head.mrho = 0
      head.memis = 0
      head.ngas = ngas
      do iI = 1,ngas
        head.glist(iI) = glist(iI)
        head.gunit(iI) = 10
        end do
      head.mlevs = mlevs
      head.nchan = 0
      head.pmin = pmin
      head.pmax = pmax
c      head.udef1 = -9999
c      head.udef2 = -9999

      RETURN
      END
c************************************************************************
