C=======================================================================
C=======================================================================
C
C        GODDARD SPACE FLIGHT CENTER
C
C        AIRS
C
C        SUSSKIND RETRIEVALS
C
!F77====================================================================

!ROUTINE NAME:  SET_SOLAR
C
!ABSTRACT:
C
C    Read in the AIRS solar radiance data
C
!CALL PROTOCOL:
C
C    SUBROUTINE set_solar(iounit, filename, Nchan, freq, Tsun, Rsun )

!INPUT PARAMETERS:
C
C   Name        Type     DESCRIPTION
C   ----        ----     -----------
C   iounit      integer  unit number of open and read file, must be >=10
C   filename    char*80  file name w/ path
C   Nchan       integer  number of IR channels in FREQ()
C   freq        real     frequency list  (1:Nchan)   cm-1
C   Tsun        real     optional parameter for using Planck Rsun(i) = function

!OUTPUT PARAMETERS:
C
C   Name        Type     DESCRIPTION
C   ----        ----     -----------
C   Rsun        real     solar radiance  mW/(M^2.str.cm-1)
C
!INPUT/OUTPUT PARAMETERS:  None
C
!RETURN VALUES:  
C
!PARENT(S):
C
C   trfinit
C
!ROUTINES CALLED:  None
C
!FILES ACCESSED:  None
C
!COMMON BLOCKS:  None
C
!DESCRIPTION:
C
C    Sets the solar radiance data for each AIRS channel.  This is the
C    effective solar radiation at the "surface" of the Sun.
C    It is  convolved with the instrument channel spectral response function.
C
C    If the unit number, iounit, is less than 10 then an error
c      message is generated
c
c    if the filename is equal to 'not.set' then no file is loaded
c      and the Planck function values are used.
C
C    The frequencies will be tested against those that are loaded for
C    consistency.  If the microwave channels are below the IR channels,
C    then the call should be made as (Nall_IR = # of IR channels, offset
C    equals index of first IR channel))
C
C      call set_solar(iounit,filename, Nall_IR, FREQ(offset), Rsun(offset))
C    
C
C    If the file exists:
C      Channels above 1900 cm-1 MUST be set via the file, otherwise an
C      error occurs.
C      Channels between 350 < FREQ(i) <= 1900 will be  set equal to the
C      Planck value if it is not set within the file.
C      Channels where FREQ(i) <= 350 will be set to zero if they are not
C      set within the file.
C
C    How the file was generated:
C    ---------------------------
C
C    The effective solar surface radiance is converted from irradiance
C    measurements via Earth based instruments (e.g., see Thekaehara, 1973).
C    Therefore, it is converted from irradiance at the time of observation
C    to an effective solar radiance at the surface of the Sun for the mean
C    Earth-Sun distance of 1 AU.  For example, Thekaekara has an irradiance
C    value  of 0.530 mW/cm^2/um at 4.6 um == 2173.9125 cm-1.
C    This is converted to W/M^2/cm^-1/steradian as follows:
C       Irradiance = Om*Blam
C         Om=steradians = pi*RADIUS_SUN^2/DISTANCE_SUN^2
C         Blam=Planck function w.r.t. wavelength  mW/cm^2/um/steradian
C    We want Bnu, which is Planck function w.r.t. wavenumber.  The two
C    Planck functions are equal in the sense of Flux (milli-Watts/cm^2)
C       Om*Bnu*dnu = Om*Blam*dlam
C         dnu = wavenumber interval, cm^-1
C         dlam = wavelength interval, microns
C       Bnu = (dlam/dnu)*(Om*Blam)/Om  ! milli-Watts/cm^2/steradian/cm^-1
C         lam = 10^4/nu    (actually wavenumber in vacuum, wavelength in air)
C         dlam/dnu = 10^4/nu^2 = lam^2/10^4
C         Om = pi*(6.951e+10/1.4957e+13)^2 = 6.78509e-5
C       Bnu(2173.9) = (4.6^2/10^4)*(0.53)/6.78509e-5 = 16.5286 mW/cm^2/st/cm^-1
C       Bnu(2173.9) = (10^4 cm^2/M^2) * 16.5286 = 1.65286E+5 mW/M^2/st/cm^-1
C       Bnu(2173.9) = (10^-3 W/mW) * 1.65286E+5 = 165.286 W/M^2/st/cm-1
C
!ALGORITHM REFERENCES:
C
C       Thekaekara, M.P. 1973.  Solar energy outside the Earth's atmosphere. 
C         Solar Energy v.14 p.109-127.
C
!KNOWN BUGS AND LIMITATIONS:  None
C
C
C    The file must have the same channel numbering system of the IR RTA,
C    which is also the index within FREQ(i).
C
C    Microwave channels are NOT set with this function.
C
C    The wrong constants for RADIUS_SUN and DISTANCE_SUN are used
C    to be compatible with previous versions of the code and with
C    the routine set_orbit.F.  The  value of Omega using these
C    values is 6.78509e-05.  These should be changed to the values below.
C      DISTANCE_SUN == 1 AU == 1.49600E+11 meters  (1994 AENA pg. K-5)
C      RADIUS_SUN (IAU 1994 definition, 1994 AENA pg. L-5) = 0.696000E9 meters.
C    which will yield a value of Omega = 6.79994e-05.  This error translates
C    to an error in AU of about 0.001, and is much larger than the other errors
C    in AU (see sun_distance.F) if an inconsistentcy between set_orbit.F
C    and the file read by set_solar.F were to exist.
C
!ROUTINE HISTORY:
C   Date      Programmer      Description
C   8/98      Scott Hannon    original program
C  10/98      C. Barnet       converted program for GSFC "wrapper"
C   1/99      C.Barnet        added error if iounit < 10
C
!END====================================================================
C
C
C=======================================================================

       subroutine set_solar(iounit, filename, Nchan, freq, Tsun, Rsun )
       implicit none

       real*4 Planck

C Input Vars
       integer*4    iounit     ! unit to use for I/O
       character*80 filename   ! name of UMBC solar file
       integer*4    Nchan      ! number of channels to fill
       real*4       freq(*)    ! frequency of channels
       real*4       Tsun       ! Planck temperature of Sun, Kelvin [option]
       real*4       Rsun(*)    ! solar radiance mW/M^2/steradian/cm^-1

C Local Vars
       CHARACTER*80 CLINE
       INTEGER*4    I
       integer*4    icnt
       INTEGER*4    IERR
       INTEGER*4    ICHAN
       integer*4    NIR
       REAL*4       FRQCHN
       real*4       FRQ_IR           ! lowest IR frequency, highest uW freq
       parameter    (FRQ_IR=350.0)
       real*4       FRQLIM
       parameter    (FRQLIM=0.01)    ! freq's must be accurate to this value
       REAL*4       SUNCHN

       character*80 filex

C***********************************************************************
C                    EXECUTABLE CODE
C***********************************************************************
C
C      -------------------------------------
C      Set default radiance for all channels
C      -------------------------------------

      do i = 1, Nchan
        Rsun(i) = 0.0
      enddo

      if(iounit.lt.10) then
        print 1011, iounit
        i = 8
        call softexit('SET_SOLA',i) 
      endif

      call strippath(filename,filex)
      if(filex.eq.'not.set') then
        do i = 1, Nchan
          if(freq(i).gt.350.0) then
            Rsun(i) = Planck(freq(i),Tsun)
          endif
        enddo
        goto 9999
      endif

C     --------------------------------------------------------------
C     Open the solar radiance file and overwrite default values with
C                 values from file
C     --------------------------------------------------------------

      OPEN(UNIT=iounit,FILE=filename,FORM='FORMATTED',STATUS='OLD',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        print 1020, IERR, filename
        i = 4
        call softexit('SET_SOLA',i)
      ENDIF
C
C     Initialize the channel counter

C     -----------------------
C     Read the solar rad file
C     -----------------------
C     Read a line of text from the file
      ICNT = 0
      ICHAN = 0
   10 ICHAN = ICHAN + 1
      IF(ICHAN.gt.NCHAN) goto 910       ! got them all
   
      if(FREQ(ICHAN).lt.FRQ_IR) then
        ICHAN = ICHAN + 1
        goto 10           ! look for IR channel
      endif

c     OK, FREQ(ICHAN)   is a valid IR channel, now find it in file

   20 READ(iounit,9000,END=910) CLINE
C     Determine if the text line is data or a comment
      IF (CLINE(1:1) .EQ. '!') goto 20

      ICNT = ICNT + 1

C     Read the data from the text line
      READ(CLINE,*)  I, FRQCHN, SUNCHN

      if( ABS(FRQCHN-freq(ICHAN)) .GT. FRQLIM ) goto 20  ! not in user list

      Rsun(ICHAN) = SUNCHN*1000.0   !convert Watts to mW

      goto 10

  910 CLOSE(iounit)   ! either EOF or all FREQ()'s are matched

C      -------------------------------------------------------------------
C      test to make sure than the channels that matter have a value loaded
C      -------------------------------------------------------------------

   30  ICNT = 0
       NIR  = 0 
       do i = 1, Nchan
         FRQCHN = FREQ(i)
         if(FRQCHN.gt.1900.0) NIR = NIR + 1
 
         if(Rsun(i).le.0.0) then
C
C          first, ensure that the short-wave channels were set
           if(FRQCHN.gt.1900.0) then
             ICNT = ICNT + 1
             print 4510, i, FRQCHN
C
C          for long-wave, the Planck approximation will be substituted
           elseif(FRQCHN.gt.FRQ_IR)  then
             Rsun(i) = Planck(FRQCHN,TSUN)
C
C          anything less than FRQ_IR (cm-1 IR or GHz microwave) will have zero
           endif
         endif
       enddo

       if( ICNT.gt.0 ) then
         print 4500, ICNT, NIR
         i = 6
         call softexit('SET_SOLA',i)
       endif

 9999  RETURN
 1011  format('set_solar: Error, iounit=',i3,' .. must be >= 10')
 1020  FORMAT('set_solar: Error ',I5,' opening file:',/,A80)
 1040  FORMAT('set_solar: Error! Channel number is out of range.',/,
     $       'Range is 1 to ',I4,', but channel list has ',I7,'.')
 9000  FORMAT(A80)
 3000  format('set_solar: file and unit=',i3,' MUST exist, file=',/,
     1     a80)
 4400 format('set_solar: incorrect freq for channel!',/,
     1       '  FRQCHN: ',f10.4,'  FREQ: ',f10.4)  
 4500 format('set_solar: # of short wave number channels NOT set =',
     $       i6,'/',i6)
 4510 format('set_solar: FREQ(',i5,')=',f10.4,' is NOT set')
       END
