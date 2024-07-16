!      this takes in a TestPointProf text file and changes it to a RTP file
!      Sergio De Souza-Machado
!      very closely linked to the ftest1.f file that Howard Motteler wrote
!
! make clean; make -f Makefile_f77 makeRTPfile.x
!
        program makeRTPfile
        implicit none

! RTP declarations
        include 'rtpdefs.f'
        integer rtpopen, rtpread, rtpwrite, rtpclose
        record /RTPHEAD/ head
        record /RTPPROF/ prof
        record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)

! other variables
        integer status, rchan, iPreSet
        character*80 mode
        CHARACTER*80 FIN, FOUT

        INTEGER NARGS, iI, J
        INTEGER IARGC, IOERR
        CHARACTER*80 BUF
        CHARACTER*80 VAL
        CHARACTER*80 VAR

        iPreset = -1
        if (iPreset .gt. 0) then
! set only one header attribute
          hatt(1).fname = 'glist'//char(0)
          hatt(1).aname = 'note'//char(0)
          hatt(1).atext = 'translated from text sonde profiles'//char(0)
! set some sample profile attributes
          patt(1).fname = 'ptemp'//char(0)
          patt(1).aname = 'units'//char(0)
          patt(1).atext = 'degrees K'//char(0)
        else
          CALL rtpinit(head,prof) 
        endif

!************************************************************************
! command line argument processing, from Scott's earlier codes
!      ------------
!      Set defaults
!      ------------
       IOERR = 6
       FIN =  'profile.txt'           ! input filename
       FOUT = 'profile_ip90.rtp'      ! output filename

!      -----------------------------------------------------------------
!      Loop on program parameters
!      --------------------------
!      Determine the number of command-line arguments
       NARGS=IARGC()

       IF (NARGS .EQ. 0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('makeRTPfile.f90 must be run with at least 1 argument')
          WRITE(IOERR,1011)
 1011     FORMAT('   fin  = <filename> txt input file  {mandatory}','//'  &
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

! read in the point profile
        CALL ReadProfile(prof,head,fin)
        print *,'head ngas  = ',head.ngas
        print *,'     glist = ',head.glist(1:head.ngas)
        print *,'     gunit = ',head.gunit(1:head.ngas)
        write(*,'(A,I3,2(F12.4),ES12.4)') 'prof  nlevs,stemp,spres,sum(WV) = ', &
              prof.nlevs,prof.stemp,prof.spres,sum(prof.gamnt(1:prof.nlevs,1))
        print *,'     rlat,rlon                 = ',prof.rlat,prof.rlon
        print *,'        iI      P(z)           T(z)           WV(z)'
        print *,'-----------------------------------------------------------'
        do iI = 1,prof.nlevs
          print *,'before write',iI,prof.plevs(iI),prof.ptemp(iI),prof.gamnt(iI,1)
        end do
        print *,'before write : nlevs = ',prof.nlevs
!        do iI = 1+prof.nlevs,MAXLEV
!          print *,'before write',iI,prof.plevs(iI),prof.ptemp(iI),prof.gamnt(iI,1)
!        end do
        print *,'-----------------------------------------------------------'
        
! create the file, write the attributes, write the (mostly unfilled) 
! profile structure, and close the file
!
        mode = 'c'
        print *,'Writing output in running dir ..'
        print *,'output name = ',fout

        status = rtpopen(fout, mode, head, hatt, patt, rchan)
        print *, 'open status f77 code,rchan = ', status,rchan
        status = rtpwrite(rchan, prof)
        !print *, 'write status F77 code = ', status
        status = rtpclose(rchan)
        !print *, 'close status F77 code = ', status

        stop
        end

!************************************************************************
! this subroutine prompts the user for the name of the point profile,
! opens it, reads in the info and saves it in the "prof" record

! the info in the levels profile is as follows 
! !!! first few lines are comments
! cType,mlevs,ngas    where cType = 'p' or 'P'
! gasID list
! lat lon
! rPress,rTemp,rDensity,(raGasAmts(iI),iI = 1,ngas)

      SUBROUTINE ReadProfile(prof,head,caFname)

      implicit none

      include 'rtpdefs.f'

      record /RTPHEAD/ head
      record /RTPPROF/ prof

      INTEGER ngas,nlevs,glist(MaxGas),ispres,istemp,thetype
      REAL pmin,pmax,spres,stemp
      CHARACTER*80 caFname
      CHARACTER*280 caStr
      INTEGER iIOUN2,iErrIO
      INTEGER iI,iZ,iL
      REAL rLat,rLon,rZ
      REAL rPress,rTemp,rDensity,raGasAmts(MAXGAS)
      CHARACTER cType

      pmin = 1.0e30
      pmax = -1.0e30

      iIOUN2 = 15

! now loop iNpath/iNumGases  times for each gas in the user supplied profile 
!      print *,caFname

      OPEN(UNIT=iIOun2,FILE=caFname,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO) 
      IF (iErrIO .NE. 0) THEN 
          WRITE(*,1070) iErrIO
        ENDIF 
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
      ENDIF

      IF ((cType .EQ. 'p') .OR. (cType .EQ. 'P')) THEN
        thetype = 1
      ELSEIF ((cType .EQ. 'h') .OR. (cType .EQ. 'H')) THEN
        thetype = 2
      END IF
      print *,'reading in textfile : thetype = ',thetype
 
      IF (ngas .GT. MAXGAS) THEN
        print *,'Error : makeRTPfile can only handle ',MAXGAS,' gases'
        CLOSE(iIOUN2)
        STOP
      ENDIF

      IF (nlevs .GT. MAXLEV) THEN
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

      ispres = -1
      istemp = -1
      !! read the spres and stemp
      READ (iIOUN2,5030,ERR=13,END=13) caStr 
      READ (caStr,*,ERR=13,END=13) spres,stemp
      if (spres > 0) ispres = +1
      if (stemp > 0) istemp = +1

      prof.nlevs = nlevs
      head.mlevs = nlevs
        
      prof.rlat = rlat
      prof.rlon = rlon
      prof.rtime = 2.0e9
      prof.plat = rlat
      prof.plon = rlon
      prof.ptime = 2.0e9

      !!read the profile info
      DO iL = 1,nlevs
        iZ = nlevs + 1 - iL  ! reverse order of loop index

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
          prof.palts(iZ) = rZ * 1000.0
          print *,'palts',iZ,prof.palts(iZ),rPress,rTemp
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

      prof.scanang = -10.0
      prof.satzen  = 15.0
      prof.solzen  = 130.0

      prof.zobs = 705000.0
      prof.scanang = 0.0
      prof.satzen  = 0.0
      prof.solzen  = 130.0

!!!!! blah blah
!      prof.nemis = 2
!      prof.efreq(1) = 600.0
!      prof.efreq(2) = 3000.0
!      prof.emis(1) = 0.9
!      prof.emis(2) = 0.9

!!! realistic ocean emiss : [h,ha,p,pa] = rtpread('/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.rtp');
      prof.nemis = 19
      prof.efreq(1:19) = (/0.7692,0.8000,0.8333,0.8696,0.9091,0.9524,1.0000,1.0417,1.0870,1.1364,1.1905,1.2500,2.4390,2.5000,2.5641,2.6316,2.7027,2.7778,2.8571/) * 1000;
      prof.emis(1:19)  = (/0.9767,0.9829,0.9889,0.9919,0.9925,0.9913,0.9891,0.9875,0.9863,0.9851,0.9840,0.9831,0.9776,0.9770,0.9764,0.9756,0.9745,0.9733,0.9715/)

      prof.rho = (1-prof.emis)/3.1415
      prof.rho = (1-prof.emis)/3.1415

! set HEAD values
      head.ptype = LEVPRO
      head.pfields = PROFBIT
!      head.mrho = 0
      head.memis = 19
      
      if (head.memis .NE. prof.nemis) then
        write(*,'(A,I3,I3)') 'head.memis .NE. prof.nemis : STOPPING ',head.memis,prof.nemis
        STOP
      END IF

      head.ngas = ngas
! see /asl/packages/klayersV205/Doc/gas_units_code.txt     10 = ppmv     
      do iI = 1,ngas
        head.glist(iI) = glist(iI)
        head.gunit(iI) = 10
      end do

      head.nchan = 2
      head.ichan(1) = 1520
      head.vchan(1) = 1231.22
      head.ichan(2) = 1820
      head.vchan(2) = 1419.9
      head.pmin = pmin
      head.pmax = pmax
!      head.iudef = -9999
!      head.itype = -9999

      prof.salti = 0.0
      write(*,'(A,I3,4(F12.5))') 'inside ReadProf : nlevs, rlat, rlon, stemp, spres = ', &
             prof.nlevs,prof.rlon,prof.rlat,prof.stemp,prof.spres  

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
