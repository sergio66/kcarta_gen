!      this takes in a TestPointProf text file and changes it to a RTP file
!      Sergio De Souza-Machado
!      very closely linked to the ftest1.f file that Howard Motteler wrote
!
! so eg      makeRTPfile.x fin=day_py4cats.atmX                      will produce  makeRTPfile.ip.rtp
! so eg      makeRTPfile.x fin=day_py4cats.atmX fout=junk.ip.rtp     will produce  junk.ip.rtp

      program makeRTPfile
      implicit none

! RTP declarations
      include 'rtpdefs.f90'
      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof,profR,profR1
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)

! other variables
      integer status, IOPCI, IOPCO, iPreSet, iI, iP

! for vararg
      CHARACTER*80 FIN, FOUT, mode
      INTEGER IARGC, IOERR
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
       FIN =  'profile.txt'           ! input filename
       FOUT = 'makeRTPfile.ip.rtp'    ! output filename
       FOUT = 'profile_ip.rtp'        ! output filename

!      -----------------------------------------------------------------
!      Loop on program parameters
!      --------------------------
!      Determine the number of command-line arguments
       NARGS=IARGC()

       IF (NARGS .EQ. 0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('makeRTPfile.f90 must be run with at least 1 argument')
          WRITE(IOERR,1011)
 1011     FORMAT('   fin  = <filename> txt input file  {mandatory}','//'    &
                 '   fout = <filename> rtp output file {optinal, set to profile_ip.rtp}')
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
      iPreset = -1
      IF (iPreset .gt. 0) THEN
! set only one header attribute
        hatt(1).fname = 'glist'//char(0)
        hatt(1).aname = 'note'//char(0)
        hatt(1).atext = 'translated from text sonde profiles'//char(0)
! set some sample profile attributes
        patt(1).fname = 'ptemp'//char(0)
        patt(1).aname = 'units'//char(0)
        patt(1).atext = 'degrees K'//char(0)
      ELSE
        CALL RTPINIT(head,prof)
      END IF

!      print *,FIN
!      print *,FOUT

! read in the point profile
      CALL ReadProfile(prof,head,FIN)
      print *,'head ngas  = ',head.ngas
      print *,'     glist = ',head.glist(1:head.ngas)
      print *,'     gunit = ',head.gunit(1:head.ngas)
      write(*,'(A,I3,2(F12.4),ES12.4)') 'prof  nlevs,stemp,spres,sum(WV) = ', &
        prof.nlevs,prof.stemp,prof.spres,sum(prof.gamnt(1:prof.nlevs,1))
      print *,'     rlat,rlon                 = ',prof.rlat,prof.rlon
      print *,'        iI      P(z)           T(z)           WV(z)'
      print *,'-----------------------------------------------------------'
      do iI = 1,prof.nlevs
        print *,'write',iI,prof.plevs(iI),prof.ptemp(iI),prof.gamnt(iI,1)
      end do
      print *,'-----------------------------------------------------------'

! create the file, write the attributes, write the (mostly unfilled) 
! profile structure, and close the file
! 'c'=create, 'r'=read : see eg /home/sergio/OTHERSTUFF/rtpV201/src/rtpopen.c
      mode = 'c'
      write(*,'(A,A)') 'Writing output in running dir .. output name = ',FOUT

!      IOPCO = 1
      status = rtpopen(FOUT, mode, head, hatt, patt, IOPCO)
      print *, 'write open status, IOPCO  = ', status,IOPCO
      status = rtpwrite(IOPCO, prof)
      print *, 'write status = ', status
      status = rtpclose(IOPCO)
      print *, 'write close status = ', status

      print *,' '
      print *,' '
      print *,' '

!************************************************************************

      IP = +1
      if (iP > 0) THEN
        mode = 'r'
!       FOUT = 'pin_feb2002_sea_airs45deg_op.so2.latlon.const_emiss_0.9.rtp'
        write(*,'(A,A)') 'now reading in file ',FOUT
        status = rtpopen(FOUT, mode, head, hatt, patt, IOPCI)
        print *, 'read open status, IOPCI  = ', status,IOPCI

        iP = 0
        do while (.true.)
          status = rtpread(IOPCI, profR1)
          if (status .eq. -1) goto 22

          profR = profR1
          iP = iP + 1
          write(*,'(A,I3,2(F12.4),ES12.4)') 'xprof nlevs,xstemp,xspres,sum(WV) = ', &
           profR.nlevs,profR.stemp,profR.spres,sum(profR.gamnt(1:profR.nlevs,1))
          print *,'     rlat,rlon                 = ',profR.rlat,profR.rlon

!           do iI = 1,profR.nlevs
!             print *,'read',iI,profR.plevs(iI),profR.ptemp(iI),profR.gamnt(iI,1)
!           end do

        end do
  22    continue
        write(*,'(A,I3,A)') 'successfully read ',iP,' profiles .... here is P(z),T(z),WV(z) for last one'
        write(*,'(A,I3,2(F12.4),ES12.4)') 'xprof nlevs,xstemp,xspres,sum(WV) = ', & 
          profR.nlevs,profR.stemp,profR.spres,sum(profR.gamnt(1:profR.nlevs,1))
        print *,'     rlat,rlon                 = ',profR.rlat,profR.rlon
        print *,'        iI      P(z)           T(z)           WV(z)'
        print *,'-----------------------------------------------------------'
        do iI = 1,profR.nlevs
          print *,'read',iP,iI,profR.plevs(iI),profR.ptemp(iI),profR.gamnt(iI,1)
        end do
        print *,'-----------------------------------------------------------'

        status = rtpclose(IOPCI)
        print *, 'read close status = ', status
      end if

      stop
      end

!************************************************************************
! this subroutine prompts the user for the name of the point profile,
! opens it, reads in the info and saves it in the "prof" record

! the info in the levels profile is as follows 
! !!! first few lines are comments, starting with ! or % or #
! cType,nlevs,ngas    where cType = 'p' or 'P'
! gasID list
! lat lon
! spres stemp         -9999 if just set it to max pres, temp at max pres
! loop   ii = 1 : nlevs
!   rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)

      SUBROUTINE ReadProfile(prof,head,caFname)

      implicit none

      include 'rtpdefs.f90'

      record /RTPHEAD/ head
      record /RTPPROF/ prof

      INTEGER ngas,nlevs,glist(MaxGas),thetype,ispres,istemp
      REAL pmin,pmax
      CHARACTER*80 caFName
      CHARACTER*280 caStr
      INTEGER iIOUN2,iErrIO
      INTEGER iI,iZ
      REAL rLat,rLon,spres,stemp
      REAL rPress,rTemp,rDensity,raGasAmts(MAXGAS),rZ
      CHARACTER*1 cType

      pmin = 1.0e30
      pmax = -1.0e30

      iIOUN2 = 15

! now loop iNpath/iNumGases  times for each gas in the user supplied profile 
!      print *,caFName

      OPEN(UNIT=iIOun2,FILE=caFName,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO) 
      IF (iErrIO .NE. 0) THEN 
        WRITE(*,1070) iErrIO
      END IF 
 1070 FORMAT('Error number ',I5,' opening your profile')

      !!!!!!!! go thru the comments at the top of the file
 30   CONTINUE 
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      IF ((caStr(1:1) .EQ. '!') .OR. (caStr(1:1) .EQ. '%') .OR. (caStr(1:1) .EQ. '#')) THEN 
        GO TO 30 
      END IF

      !!read type of sonde profile (pressure or height), num levels, num gases
      READ (caStr,*,ERR=13,END=13) cType,nlevs,ngas 
      IF ((cType .NE. 'p') .AND. (cType .NE. 'P') .AND. (cType .NE. 'h') .AND. (cType .NE. 'H'))  THEN
        print *,'Error : makeRTPfile assumes pressure info in the file'
        CLOSE(iIOUN2)
        STOP
      END IF
      IF ((cType .EQ. 'p') .OR. (cType .EQ. 'P')) THEN
        thetype = 1
      ELSEIF ((cType .EQ. 'h') .OR. (cType .EQ. 'H')) THEN
        thetype = 2
      END IF
      
      IF (ngas .GT. MAXGAS) THEN
        print *,'Error : makeRTPfile can only handle ',MAXGAS,' gases'
        CLOSE(iIOUN2)
        STOP
      END IF

      IF (nlevs .GT. MAXLEV) THEN
        print *,'Error : makeRTPfile can only handle ',MAXLEV,' levels'
        CLOSE(iIOUN2)
        STOP
      END IF

      !! read the gasID list
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) (glist(iI),iI = 1,ngas)
 
      !! read the latitude and longitude
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) rLat,rLon

      ispres = -1
      istemp = -1
      !! read the spres and stemp
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) spres,stemp
      if (spres > 0) ispres = +1
      if (stemp > 0) istemp = +1
        
      prof.rlat = rLat
      prof.rlon = rLon
      prof.nlevs = nlevs

      !!read the profile info
      DO iZ = 1,nlevs
        READ (iIOUN2,5030,ERR=13,END=13) caStr
        IF (thetype .EQ. 1) THEN
          READ (caStr,*,ERR=13,END=13) rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)
!          print *,rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)        
        ELSEIF (thetype .EQ. 2) THEN
          READ (caStr,*,ERR=13,END=13) rZ,rPress,rTemp,(raGasAmts(iI),iI = 1,ngas)
!          print *,rZ,rPress,rTemp,(raGasAmts(iI),iI = 1,ngas)        
        END IF
        
        prof.plevs(iZ) = rPress
        prof.ptemp(iZ) = rTemp
        IF (thetype .EQ. 2) THEN
          prof.palts(iZ) = rZ
        END IF

        DO iI = 1,ngas
          prof.gamnt(iZ,iI) = raGasAmts(iI)
        END DO
        IF (rPress .LT. pmin) THEN
          pmin = rPress
        END IF
        IF (rPress .GT. pmax) THEN
          pmax = rPress
          IF (ispres < 0) spres = rPress
          IF (istemp < 0) stemp = rTemp
        END IF

      ENDDO

 13   CONTINUE
      CLOSE(iIOUN2)
 5030 FORMAT(A280)

      prof.spres = spres
      prof.stemp = stemp

! set HEAD values
      head.ptype = LEVPRO
      head.pfields = PROFBIT
!      head.mrho = 0
      head.memis = 0
      head.ngas = ngas
! see /asl/packages/klayersV205/Doc/gas_units_code.txt     10 = ppmv     
      DO iI = 1,ngas
        head.glist(iI) = glist(iI)
        head.gunit(iI) = 10
      END DO

      prof.nlevs = nlevs
      head.nchan = 0
      head.pmin = pmin
      head.pmax = pmax
!      head.udef = -9999
!      head.udef = -9999

      print*,'<----------- SUMMARY of INPUT TXT file -----------> '
      print *,' the gasIDs are : '
      print *, (glist(iI),iI = 1,ngas)
      print *,'  rLat,rLon = ', prof.rLat,prof.rLon
      print *,'  number of levels, number of gases =', prof.nlevs,head.ngas
      print *,'  surf pres, surf temp = ', prof.spres,prof.stemp
      print*,'<----------- SUMMARY of INPUT TXT file -----------> '
      print *,'  '

      RETURN
      END
!************************************************************************
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

!----------------------------------------------------------------------
