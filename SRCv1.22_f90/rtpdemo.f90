! pgf90  -o rtpdemo.x  -Mextend -Mnoupcase -O -Mbounds rtpdemo.f
! asl/opt/absoft/absoft10.0/bin/af77 -o rtpdemo.x -I/asl/packages/rtpV201_g80/include -w -W -C -O -A -N3 -s -f -N15 rtpdemo.f -L/asl/packages/rtpV201_g80/lib -lrtp -L/asl/opt/lib -lmfhdf -ldf -ljpeg -lz
!************************************************************************

    include '../INCLUDE/kcartaparam.f90'
    include '/asl/packages/rtpV201_g80/include/rtpdefs.f'

! caSetEmissivity= array that gives name of emissivity files (if any)
! raSetEmissivity= array that gives constant emissivity value (if set)
! iNpmix     = number of mixed paths read in from mixfile
! iaMPSetForRad = array telling which MP set to associate with which atm
! iNatm       = number of atmospheres
! raPressStart = start pressure for radiating atmos
! raPressStop  = stop pressure for radiating atmos
! raTSpace    = array containing background temperature for each atmosphere
! raTSurf    = array contianing surface temperature for each atmosphere
! raSatAngle = array containing satellite view angle for each atmosphere
! raSatHeight= array containing satellite height for each atmosphere
! iaNumLayer = array containing number of layers in each atmosphere
! iaaRadLayer= matrix containing list of layers in each atmosphere
! iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
! raaaSetEmissivity = array containing the wavenumber dependent emissivities
! raFracTop  = top fraction
! raFracBot  = bottom fraction
! raaPrBdry  = matrix that keeps start/stop pressures
! the next few only work for DOWNWARD LOOK instr
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! rakSolarRefl   =solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! raProfileTemp = array containing CO2 gas profile temperature
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
! this is for the cloud, if any, that is associated with the atmosphere
! raPressLevels are the actual pressure levels from the KLAYERS file
    INTEGER :: iProfileLayers,ctype1,ctype2
    REAL :: raPressLevels(kProfLayer+1)

    REAL :: cfrac12,cFrac1,cFrac2,cngwat1,cngwat2,raCemis(kMaxClouds)
    REAL :: raCprtop(kMaxClouds), raCprbot(kMaxClouds)
    REAL :: raCngwat(kMaxClouds), raCpsize(kMaxClouds)
    INTEGER :: iaCtype(kMaxClouds),iMPSetForRadRTP
    INTEGER :: iNclouds_RTP     !!!tells how many clouds
    INTEGER :: iaNML_Ctype(kMaxClouds)

    CHARACTER(160) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm),iNpmix
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
    INTEGER :: iRTP     !!!tells which profile info, radiance info, to read
    CHARACTER(130) ::  caPFName !!!tells which profile

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNlay,iStart,iStop,iErr
    REAL :: rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight,rT
    INTEGER :: iDirection,iW,iInt
    INTEGER :: iC,iNumLinesRead
    REAL :: FindSurfaceTemp,rSize1,rSize2
    INTEGER :: iInstrType

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: i,j,k,iG,upwell,iOKscanang,iOKzobs,iOKsatzen
    REAL :: raHeight(kProfLayer+1),raThickness(kProfLayer),pobs,pobs1,pTemp,rSURFaltitude
    REAL :: r1,rEms,rAngleX,rAngleY
!      REAL saconv_sun,orig_saconv_sun
    INTEGER*4 :: i4CTYPE1,i4CTYPE2
     
    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(32) :: mode
    character(160) :: fname
    real :: rf1,rf2

    iRTP = 1
    caPFName = '/home/sergio/MATLABCODE/RTPMAKE/100LAYER/test_1chan.op.rtp'

    print *,'booboo A ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    fname(1:80) = caPFName(1:80)

    write(kStdWarn,*) 'Using RTP file to set atm info ....'
    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Stop    !Call Dostop
    END IF

    print *,'X1A prof.ctype,prof.ctype2 = ',prof.ctype,prof.ctype2
    print *,'booboo B ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    kProfileUnitOpen = +1
! ince we succesffuly read in profile, no need to do as many checks here
    DO i = 1, iRTP
        status = rtpread(rchan, prof)
    END DO
    status = rtpclose(rchan)
    kProfileUnitOpen = -1

    print *,'X1B prof.ctype,prof.ctype2 = ',prof.ctype,prof.ctype2
    print *,'booboo c ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    rf1 = head.vcmin
    rf2 = head.vcmax

    caWord = '*RADNCE'
    iErr = -1

    iNatm = 1        !can only read ONE atmosphere per RTP profile

! now get the relevant info from rchan,prof
    iC = 1
    iW = iMPSetForRadRTP
    iaMPSetForRad(1) = iMPSetForRadRTP

!!!assume the instrument is downlooking, from TOA
    pobs = 0.0
         
    rSURFaltitude = prof.salti
          
    rPressStart = prof.spres
    rPressStop  = 0.000            ! ----------> assume TOA

!!!then go ahead and look at variables prof.pobs
!!!note that variable pobs is reset only if prof.pobs > 0, else it
!!!stays at a value of 0.0000 (TOA)
    pobs1 = prof.pobs
    IF (pobs1 > 0) THEN
        IF (pobs1  < raPressLevels(iProfileLayers+1)) THEN
            write(kStdWarn,*) 'From reading info in RTP file, reset prof.pobs'
            write(kStdWarn,*) 'frm ',pobs1,' to ',raPressLevels(iProfileLayers+1)
            pobs1 = raPressLevels(iProfileLayers+1)
            pobs = pobs1
        END IF
        upwell = ((kProfLayer + 1) - prof.nlevs) + 1
        IF (pobs1  > raPressLevels(upwell)) THEN
            write(kStdWarn,*) 'From reading info in RTP file, reset prof.pobs'
            write(kStdWarn,*) 'from ',pobs1,' to ',raPressLevels(upwell)
            pobs1 = raPressLevels(upwell)
            pobs = pobs1
        END IF
    END IF
    pobs = pobs1

! testing
!      prof.satzen = -abs(prof.satzen)
!      prof.zobs   = -1000
!       prof.scanang = -abs(prof.scanang) * 1000

!!!assume the instrument is downlooking, from TOA
    upwell = 1
!!!then go ahead and look at variables prof.upwell
    IF (prof.upwell == 1) THEN
    ! radiation is travelling upwards so this is a downlook instrument
        upwell = prof.upwell
    ELSEIF (prof.upwell == 2) THEN
    ! radiation is travelling downwards so this is a uplook instrument
        upwell = prof.upwell
    ELSE
        write(kStdErr,*) 'need prof.upwell = 1 (downlook) or 2 (uplook)'
        write(kStdErr,*) 'prof.upwell = ',prof.upwell
        Stop    !Call Dostop
    END IF

! ow that we have rPressStart,rPressStop (defining pressure boundaries
! or the atm) and upwell (direction of radiation travel), check to see
! f things need to be reset
    IF (upwell == 1) THEN
    ! eed rPressStart > rPressStop
        rPressStart = prof.spres
        rPressStop  = pobs
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof.pobs,prof.upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a downlook instr, with Surf,OBS press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart <= rPressStop) THEN
            write(kStdErr,*) 'For downlook instr, need rPressStart > rPressStop'
            Stop    !Call Dostop
        END IF
    ELSEIF (upwell == 2) THEN
    ! eed rPressStart < rPressStop
        rPressStart = 0.0
        rPressStop  = prof.spres
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof.pobs,prof.upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a uplook instr, with TOA,Surf press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart >= rPressStop) THEN
            write(kStdErr,*) 'For uplook instr, need rPressStart < rPressStop'
            Stop    !Call Dostop
        END IF
    ELSE
        write(kStdErr,*) 'Need to set prof.upwell = 1 or 2'
        Stop    !Call Dostop
    END IF

    rTbdy       = kTSpace          ! -------> assume deep space
    rTSurf      = prof.stemp

!      DO i = 1, prof.nlevs - 1
!        rT   = prof.ptemp(i)
!      END DO
!      print *,rTSurf,rT,rT+1
!      rTSurf      = rT+1

!!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scott has this in sarta.f (for AIRS!!)
!      Convert SATZEN or SATANG to viewing angle
!       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
!         Convert zenith angle at surface to view angle at satellite
!          SVA=SACONV( SATZEN, SALT*1000 )/CONV
!       ELSE
!         Check if scan angle is valid
!          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
!            View angle should be within a few degrees of scan angle
!             SVA=ABS( SATANG )
!          ELSE
!             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
! 1030        FORMAT('Error! Profile',I5,
!     $          ': invalid angles for SATZEN ',1PE11.4,
!     $          ' and SATANG ',E11.4)
!             STOP
!        ENDIF
!     ENDIF
!       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
!       IF (SVA .GT. ANGMAX) THEN
!         Truncate angle if too big
!          WRITE(IOINFO,1040) IPROF, SVA
! 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
!     $       1PE11.4,' to 53 degrees')
!          SVA=ANGMAX
!     ENDIF
!      Convert SATZEN or SATANG to viewing angle
!       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
!         Convert zenith angle at surface to view angle at satellite
!          SVA=SACONV( SATZEN, SALT*1000 )/CONV
!       ELSE
!C         Check if scan angle is valid
!          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
!            View angle should be within a few degrees of scan angle
!             SVA=ABS( SATANG )
!          ELSE
!             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
! 1030        FORMAT('Error! Profile',I5,
!     $          ': invalid angles for SATZEN ',1PE11.4,
!     $          ' and SATANG ',E11.4)
!             STOP
!        ENDIF
!     ENDIF

!       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
!       IF (SVA .GT. ANGMAX) THEN
!         Truncate angle if too big
!          WRITE(IOINFO,1040) IPROF, SVA
! 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
!     $       1PE11.4,' to 53 degrees')
!          SVA=ANGMAX
!     ENDIF

! assume scanang, satzen, zobs make sense
    write(kStdWarn,*) ' '

    iOKscanang = +1
    IF (abs(prof.scanang) > 89.99) THEN
        write(kStdWarn,*) 'whoops : prof.scanang = ',prof.scanang
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKscanang = -1
    END IF

    iOKsatzen = +1
    IF (abs(prof.satzen) > 89.99) THEN
        write(kStdWarn,*) 'whoops : prof.satzen = ',prof.satzen
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKsatzen = -1
    END IF

    iOKzobs = +1
    IF ((prof.zobs < 2.00*1000) .AND. (prof.upwell == 1)) THEN
        write(kStdWarn,*) 'whoops : prof.zobs = ',prof.zobs
        write(kStdWarn,*) '         zobs < 2 km in height!'
        write(kStdWarn,*) '         does not make sense for downlook instr'
        iOKzobs = -1
    END IF

    rAngle      = prof.scanang     ! -------> ignore satzen,satazi
    rHeight     = prof.zobs        ! -------> use raSatHeight(iC) in m (changed in July 2013)

    IF ((iOKscanang == -1) .AND. (iOKsatzen == -1) .AND. &
    (iOKzobs == -1)) THEN
        write(kStdWarn,*) 'scanang,satzen,zobs do not make sense!'
        write(kStdErr,*) 'scanang,satzen,zobs do not make sense!'
        Stop    !Call Dostop
    END IF

!! Larrabee prefers to use satzen, so if it exists, use that!!!!!
! > Secondly, I was looking into Larrabee's request to get kCARTA to handle
! > p.satzen rather than p.scanang. It'll take some doing, but in Matlab,
! > suppose I have p.satzen. Do I use
! >     [zang]=saconv( p.satzen,705000);
! > to convert that into scanang, for kCARTA to use?
! Yes, except use prof.zobs rather than a hardcoded 705000.
    rAngleX     = prof.satzen      ! -------> rtp_interface originally
!          ignored satzen,satazi
! so the conversion is  p.scanang = orig_saconv_sun( p.satzen,prof.zobs);  %% by Scott
! so the conversion is  p.scanang = saconv( p.satzen,prof.zobs);           %% by Sergio
    write(kStdWarn,*) ' '
    IF (prof.upwell == 1) THEN
        IF (rHeight > 2.0) THEN
        !          rAngleY = ORIG_SACONV_SUN(rAngleX, rHeight)
        !          rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, rHeight/1000)
            rAngleY = rAngleX
        ELSE
        !          rAngleY = ORIG_SACONV_SUN(rAngleX, 705.00)                !! default AIRS hgt, dangerous
        !          rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, 705.00) !! default AIRS hgt, dangerous
            rAngleY = -9999
        ENDIF
        write(kStdWarn,*) 'downlook instr : satellite hgt, view angle info : '
        write(kStdWarn,*) 'scanang, zobs(km), satzen, saconv(satzen,zobs) = ', &
        rAngle,rHeight/1000,rAngleX,rAngleY
    END IF

    IF (prof.upwell == 2) THEN
    !! uplook instrument
        IF (iOKscanang == 1) THEN
            write(kStdWarn,*) 'Uplook instr : use prof.scanang'
        ELSEIF (iOKscanang == -1) THEN
            write(kStdWarn,*) 'Uplook instr : incorrect prof.scanang'
            write(kStdErr,*) 'Uplook instr : incorrect prof.scanang',rAngle
            Stop    !Call Dostop
        END IF
    ELSEIF (prof.upwell == 1) THEN
    !! downlook instr
        IF ((iOKscanang == 1) .AND. (iOKsatzen == 1) .AND. &
        (iOKzobs == 1)) THEN
        !! all angles seem reasonable; now check consistency between them
            IF (abs(abs(rAngle)-abs(rAngleY)) <= 1.0e-2) THEN
                write(kStdWarn,*) 'scanang,satzen,zobs present in rtp file'
                write(kStdWarn,*) 'scanang and saconv(satzen,zobs) agree'
                rAngle = rAngleY   !!! no need to do anything
            ELSEIF (abs(abs(rAngle)-abs(rAngleY)) > 1.0e-2) THEN
                write(kStdWarn,*) 'scanang,satzen,zobs present in rtp file'
                write(kStdWarn,*) 'scanang and saconv(satzen,zobs) disagree'
                write(kSTdWarn,*) 'using satzen (AIRS preference!!!)'
                IF (prof.zobs < 2000.0) THEN
                    write(kStdErr,*) 'used 705 km as satellite height'
                ELSE
                    write(kStdErr,*) 'used ',prof.zobs/1000 ,' km as satellite hght'
                END IF
                rAngle = rAngleY   !!! replace input scanang with that derived
            !!! from satzen,zobs
            END IF
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == -1) .AND. (iOKzobs == +1))) THEN
            write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
            rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == +1) .AND. (iOKzobs == -1))) THEN
            write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
            rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == -1) .AND. (iOKzobs == -1))) THEN
            write(kStdWarn,*) 'satzen or zobs or both do not make sense',rAngleX,rHeight*1000
            rAngle = rAngle   !!! no need to do anything
        ELSEIF ((iOKscanang == -1) .AND. &
            ((iOKsatzen == +1) .AND. (iOKzobs == +1))) THEN
        !! satzen or zobs or both make sense, scanang is incorrect
            rAngle = rAngleY   !!! no need to do anything
            write(kStdWarn,*) 'scanang does not make sense, but satzen,zobs do'
            write(kSTdWarn,*) 'using satzen to derive scanang (AIRS preference!!)'
        END IF
    END IF

    write(kStdWarn,*) 'using kCARTA scanang = ',rAngle
    write(kStdWarn,*) ' '
!!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (rTbdy > 3.0) THEN
        write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
        Stop    !Call Dostop
    END IF

!      CALL StartStopMP(iW,rPressStart,rPressStop,iC,
!     $                 raPressLevels,iProfileLayers,
!     $                 raFracTop,raFracBot,raaPrBdry,iStart,iStop)
    raPressStart(1) = raaPrBdry(1,1)
    raPressStop(1)  = raaPrBdry(1,2)

! figure out if the start/stop MixedPath numbers are legitimate
    IF ((iStart > iNpmix) .OR. (iStart < 1) .OR. &
    (iStop > iNpmix) .OR. (iStop < 1)) THEN
        write(kStdErr,*)'Error while setting Start/Stop Mixed Path '
        write(kStdErr,*)'numbers for atmosphere # ',iC
        write(kStdErr,*)'Must be between 1 and ',iNpmix
        write(kStdErr,*)'Instead they are ',iStart,iStop
    !          Stop    !Call Dostop
    END IF

! figure out how many radiating layers (or MPs) in this atmosphere, and check
! that that it is less than or equal to kProfLayer
    IF (iStop >= iStart) THEN
        iNlay = (iStop-iStart+1)
        iDirection = +1                           !down look instr
    ELSE IF (iStop <= iStart) THEN
        iNlay = (iStart-iStop+1)
        iDirection = -1                           !up look instr
    END IF
    IF (iNLay > kProfLayer) THEN
        write(kStdErr,*)'Error for atm # ',iC
        write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
        Stop    !Call Dostop
    END IF

! set the B.C.'s
    raTSpace(iC) = rTbdy
!      raTSurf(iC)  = FindSurfaceTemp(rPressStart,rPressStop,rTSurf,
!     $                     raProfileTemp,raPressLevels,iProfileLayers)
    raTSurf(iC) = rTsurf

    raSatAngle(iC)=rAngle
    IF (abs(rAngle) <= 1.0e-4) THEN !nadir view
        rHeight = -1.0
        raSatHeight(iC) = -1.0
    ELSE
        IF (rHeight > 0.0) THEN
            raSatHeight(iC) = rHeight   !height in m
        ELSE
            rHeight = -1.0
            raSatHeight(iC) = rHeight   !height in m
        END IF
    END IF
    raSatAzimuth(iC) = prof.satazi
    raSolAzimuth(iC) = prof.solazi
    iaNumLayer(iC) = iNlay

    write(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
    write(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
    write(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted= ', &
    rTsurf,raTSurf(iC)

! set the mixed path numbers for the current atmosphere, in direction of
! radiation travel
    DO iInt = 1,iNlay
        iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
    END DO

! use the solar on/off, thermal on/off etc.
! sun is only on if 0 < prof.solzen < 90
! rakSolarAngle(iC) = abs(prof.sunang)   !!!RTP v 097-
    rakSolarAngle(iC) = prof.solzen          !!!RTP v 098+
    IF ((prof.solzen >= 0.0) .AND. (prof.solzen <= 90.0)) THEN
        IF ((rf1 >= 605.0) .AND. (rf2 <= 2830.0)) THEN
            iakSolar(iC) = +1
        ELSEIF ((rf1 < 605.0) .OR. (rf2 > 2830.0)) THEN
            iakSolar(iC) = +0   !!! do not have a solar database yet
        END IF
    ELSE
        iakSolar(iC) = -1
    END IF
    raKSolarRefl(iC) = -1.0
    iaKThermal(iC)   = 0
    raKThermalAngle(iC) = -1.0

!      print *,'----> warning : set raKthermalangle = 53.3 (acos(3/5))'
!      raKThermalAngle(iC) = +53.13
!      print *,'----> so this will be used at all layers '
!      print *,'----> instead of varying the diffusivity angle'

    iakThermalJacob(iC) = 1
! use the solar on/off, thermal on/off etc.
    kSolar        = iaKSolar(iC)
    IF (abs(raKSolarAngle(iC) - 90.0) <= 1.0e-5) then
        write(kStdWarn,*) 'resetting solar angle = 90 to 89.9, iAtm = ',iC
        raKSolarAngle(iC) = 89.9
    END IF
    kSolarAngle   = raKSolarAngle(iC)
    kSolarRefl    = raKSolarRefl(iC)
    kThermal      = iaKThermal(iC)
    kThermalAngle = raKThermalAngle(iC)
    kThermalJacob = iakThermalJacob(iC)

    IF (kThermalAngle  < 0) THEN
        kSetThermalAngle = -1   !use accurate angles lower down in atm, const  in tau temp variation
        IF ((kFlux > 0) .OR. (kTemperVary >= 4)) THEN
        ! kSetThermalAngle = -2   !use accurate angles lower down in atm, linear in tau temp variation
            kThermal = +2           !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
            kSetThermalAngle = +2   !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
        END IF
    ELSE
        kSetThermalAngle = +1   !use user specified angle everywhere
    END IF
    write(kStdWarn,*) 'in rtpdemo.f --> kFlux,kTemperVary,kSetThermalAngle = ',kFlux,kTemperVary,kSetThermalAngle
          
    IF (iDirection > 0) THEN
    ! heck things make sense for downlook in
        IF ((kSolarAngle < 0.0) .OR. (kSolarAngle > 90.0)) THEN
            write(kStdWarn,*) 'Warning! Resetting Solar Angle from ',kSolarAngle,' to 150.0'
            write(kStdWarn,*) 'and setting kSolar from ',kSolar, ' to -1 (solar = off)'
            kSolarAngle = 150.0
            kSolar      = -1
        END IF
        IF ((abs(kSolar) /= 1) .AND. (kSolar /= 0)) THEN
            write(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
            Stop    !Call Dostop
        END IF
        IF ((abs(kThermal) > 1) .AND. (kThermal /= 2)) THEN
            write(kStdErr,*)'need Thermal on/off parameter = -1/0/1/2'
            Stop    !Call Dostop
        END IF
        IF (abs(kThermalJacob) /= 1) THEN
            write(kStdErr,*)'need ThermalJacob on/off parameter = -1/1'
            Stop    !Call Dostop
        END IF
    ! et the diffusivity angle in degrees
    ! F ((kThermalAngle .LT. 0.0).OR.(kThermalAngle .GT. 90.0)) THEN
        IF (kThermalAngle > 90.0) THEN
            write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
        END IF
    END IF

    IF (iDirection < 0) THEN
        IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 4)) THEN
        ! et to nonsense values for uplooking instrument RTSPEC SCAT
        ! s these CANNOT handle solar
            kSolar = -1          !!!RTPSEC, FIRST ORDER PERTURB cannot handle sun
            IF (kSolar /= iaKSolar(iC)) THEN
                write(kStdErr,*) 'in radnce4RTP, kSolar = -1 but you have a solar'
                write(kStdErr,*) 'angle for the profile!!!',kWhichScatterCode
                Stop    !Call Dostop
            END IF
            kSolarAngle   = -90.0
            kSolarRefl    = 0.0
            kThermal      = 0
            kThermalAngle = -45.0
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
            kThermalJacob = -1
        ELSE
        ! et to nonsense values for uplooking instr kCARTA, DISORT,TWOSTR
        !!!kSolar = -1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
            kSolarRefl    = 0.01
            kThermal      = 0
            kThermalAngle = -45.0
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
            kThermalJacob = -1
        END IF
    END IF

! So if {\sf iakThermal(iI) = 0}, then {\sf rakThermalAngle(iI)} should be used
! with care.  If it is set at a negative value $x$, then for the upper
! layers the diffusive angle $acos(3/5)$ is used for the reflected
! thermal, while for the lower layers, a parameterized optimum
! diffusivity angle is used.  If it is set at a positive value, then for
! all layers, the diffusive angle $acos(x)$ is used for the
! reflected thermal. acos(3/5) = 53.1301 degrees

!      print *,'**************************************************'
!      kThermal = -1
!      kSolar = -1
!      print *,'*******  in rtp_interface,  kThermal = -1 ********'
!      print *,'*******  in rtp_interface,  kSolar   = -1 ********'
!      print *,'**************************************************'
    iakSolar(iC)          = kSolar
    rakSolarAngle(iC)     = kSolarAngle
    rakSolarRefl(iC)      = kSolarRefl
    iakThermal(iC)        = kThermal
    rakThermalAngle(iC)   = kThermalAngle
    iakThermalJacob(iC)   = kThermalJacob
    iaSetThermalAngle(iC) = kSetThermalAngle

    write(kStdWarn,*)'Solar on/off, Solar angle, Solar emiss = ', &
    kSolar,kSolarAngle,kSolarRefl
    write(kStdWarn,*)'Thermal on/off,Thermal angle,Thermal Jacob =', &
    kThermal,kThermalAngle,kThermalJacob

! now read in the emissivity values
    iaSetEms(iC) = prof.nemis
    IF (iaSetEms(iC) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many emiss regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        Stop    !Call Dostop
    END IF

    IF (iaSetEms(iC) >= 2) THEN
        DO i=1,iaSetEms(iC)
            r1   = prof.efreq(i)
            rEms = prof.emis(i)
            write(kStdWarn,*) r1,rEms
            raaaSetEmissivity(iC,i,1) = r1
            raaaSetEmissivity(iC,i,2) = rEms
            IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
                write(kStdErr,*)'Need emissivity between 0 and 1'
                write(kStdErr,*)'check your emissivity values in file'
                Stop    !Call Dostop
            END IF
        END DO
    END IF

    IF ((iaSetEms(iC) == 1) .AND. (prof.emis(1) > 1.0)) THEN
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'rtpfile : has ONE emiss point, and emiss > 1'
        write(kStdWarn,*) '  assuming LIMB VIEW : RESET emissivity = 0'
        iaLimb(iC) = +1
        i = 1
        r1 = 6.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        i = 2
        r1 = 36000000.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
    ELSEIF ((iaSetEms(iC) == 1) .AND. (prof.emis(1) <= 1.0)) THEN
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'Fooling emissivity file so that it uses two points'
        write(kStdWarn,*) 'with constant emissivity',prof.emis(1)
        iaSetEms(iC) = 2
        i = 1
        r1 = 6.0
        rEms = prof.emis(1)
        IF (rEms < 0.0) THEN
            write(kStdErr,*)'Need emissivity between 0 and 1'
            write(kStdErr,*)'   RTP file has ',prof.nemis,' emissivity points'
            write(kStdErr,*)'   first point (rf,rEms) = ',prof.efreq(1),prof.emis(1)
            Stop    !Call Dostop
        END IF
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
        i = 2
        r1 = 36000000.0
        rEms = prof.emis(1)
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
    END IF
     
! now read in the solar refl values
    iaSetSolarRefl(iC) = prof.nemis !! new, before it was nrho
    IF (iaSetSolarRefl(iC) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many solar refl regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        Stop    !Call Dostop
    END IF
    IF (iaSetSolarRefl(iC) < 1) THEN
        write(kStdWarn,*)'No points in the solar refl file'
        write(kStdWarn,*)'Will assume that we are using refl = (1-ems)/pi)'
        iaSetSolarRefl(iC) = iaSetEms(iC)
        DO i = 1,iaSetEms(iC)
        ! irst is wavenumber, second is emissivity --> reflectance
            raaaSetSolarRefl(iC,i,1)  = raaaSetEmissivity(iC,i,1)
            raaaSetSolarRefl(iC,i,2)  = (1-raaaSetEmissivity(iC,i,2))/kPi
        END DO
    ELSE IF (iaSetSolarRefl(iC) < 2) THEN
        write(kStdWarn,*)'For solar refl, Need > 1 point for interpolation'
        write(kStdWarn,*)'Fooling reflectivity file so that it uses two '
        write(kStdWarn,*)'points, with constant reflectivity',prof.rho(1)
        iaSetSolarRefl(iC) = 2
        i = 1
        r1 = 3.6
        rEms = prof.rho(1)
        raaaSetSolarRefl(iC,i,1) = r1
        raaaSetSolarRefl(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
        i = 2
        r1 = 3600.0
        rEms = prof.rho(1)
        raaaSetSolarRefl(iC,i,1) = r1
        raaaSetSolarRefl(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
    ELSE
        DO i=1,iaSetSolarRefl(iC)
            r1   = prof.efreq(i)   !!new rfreq = efreq
            rEms = prof.rho(i)
            write(kStdWarn,*) r1,rEms
            raaaSetSolarRefl(iC,i,1) = r1
            raaaSetSolarRefl(iC,i,2) = rEms
            IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
                write(kStdErr,*)'Need reflectance between 0 and 1'
                write(kStdErr,*)'check your reflectance values in file'
                Stop    !Call Dostop
            END IF
        END DO
    END IF

    print *,'booboo d ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

! ow see if there is a cloud to be used with this atmosphere

    cfrac12 = prof.cfrac12

    print *,'booboo d1 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    ctype1  = int(prof.ctype)
    cfrac1  = prof.cfrac
    cngwat1 = prof.cngwat
    rSize1  = prof.cpsize

    print *,'booboo d2 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    ctype2  = int(prof.ctype2)
    cfrac2  = prof.cfrac2
    cngwat2 = prof.cngwat2
    rSize2  = prof.cpsize2

    print *,'booboo d3 ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    i4ctype1 = prof.ctype
    i4ctype2 = prof.ctype2

    print *,'booboo e ',cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2,iNclouds_RTP

    print *,'xmama0',i4ctype1,i4ctype2
    print *,'xmama1',ctype1,cfrac1,cngwat1,rSize1
    print *,'xmama2',ctype2,cfrac2,cngwat2,rSize2
    print *,'xmama12',cfrac12
    print *,prof.clrflag,prof.ctype,prof.ctype2,prof.upwell,prof.findex,prof.atrack,prof.xtrack,prof.ifov,prof.robsqual,prof.itype

    IF ((cfrac1 <= 0) .AND. (cfrac2 > 0)) THEN
        write(kStdErr,*) 'kCARTA assumes if cfrac1 > 0 then cfrac2 >= 0'
        write(kStdErr,*) 'kCARTA assumes if cfrac1 = 0 then cfrac2  = 0'
        write(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        STOP    !CALL DOSTOP
    END IF

    IF (cfrac12 > max(cfrac1,cfrac2)) THEN
        write(kStdErr,*) 'kCARTA assumes cfac12 <= max(cfrac1,cfrac2)'
        write(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        STOP    !CALL DOSTOP
    END IF

    IF (prof.cfrac > 0.0) THEN
    !!!first cloud is easy to do
        i = 1
        iaCtype(i) =  prof.ctype       !cloud type 1=cirrus 2=water etc
        raCemis(i) = 1.0               !assume cloud totally emissive
        raCngwat(i) = prof.cngwat*1.00 !IWP
        raCpsize(i) = prof.cpsize*1.00 !in microns
        raCprtop(i) = prof.cprtop
        raCprbot(i) = prof.cprbot

        i = 2
        iaCtype(i) =  prof.ctype2        !cloud type 1=cirrus 2=water etc
        raCemis(i) = 1.0                 !assume cloud totally emissive
        raCngwat(i) = prof.cngwat2*1.00  !IWP
        raCpsize(i) = prof.cpsize2*1.00  !in microns
        raCprtop(i) = prof.cprtop2
        raCprbot(i) = prof.cprbot2

    ELSE
        cfrac1  = 0.0            !assume clear sky, use dummy values
        cfrac2  = 0.0            !assume clear sky, use dummy values
        cngwat1 = 0.0
        cngwat2 = 0.0
        DO i = 1,kMaxClouds
            iaCtype(i) =  0                !cloud type
            raCemis(i)  = 0.0              !assume cloud totally emissive
            raCngwat(i) = 0.0              !IWP
            raCpsize(i) = 1.0              !in microns
            raCprtop(i) = 500.0
            raCprbot(i) = 600.0
        END DO
    END IF

    END PROGRAM

!************************************************************************
