! Copyright 2015
! University of Maryland Baltimore County
! All Rights Resered

! this file deals with reading in the RTP file
! scientific format ESW.n   http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html

! see Makefile for rtpdefs.f90 ... Makefile_intel_hdf_rtp

MODULE rtp_interface_f90

USE basic_common
USE spline_and_sort_and_common
USE n_pth_mix
USE s_misc
USE freqfile
USE rad_angles
USE n_rad_jac_scat
USE n_pth_mix

IMPLICIT NONE

CONTAINS

!************************************************************************
!           these functions deal with reading CLOUD PROFILES
!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format
! same as READRTP except now it looks for
!   h.ngas > 4 : for 4 gases + (1,2 or 3) cloud types
!   h.glist = cobinations of (201,202,203)
!     water drop           : ctype=101     in prof%gas_201
!     Baran ice aggragates : ctype=201     in prof%gas_202
!     andesite dust        : ctype=301     in prof%gas_203
!   h.gunit = 1            : molecules/cm2 for gases  01:63
!           = 5            : g/cm2         for clouds 201+

! OLD NLTE routines, need to clean these up
! this subroutine reads 48 regression profiles to see which is closest to the
! current profile, so that kCARTA can figure out which NLTE temps to use!!!
    SUBROUTINE NLTE_RegrTemp(raa48Temp,raa48Press,ia48layers)

    implicit none

    include '../INCLUDE/TempF90/kcartaparam.f90'
    include 'rtpdefs.f90'

! the 48 regression profile kinetic temperatures
    REAL :: raa48Temp(kMaxLayer,kRegrProf)
    REAL :: raa48Press(kMaxLayer,kRegrProf)
    INTEGER :: ia48layers(kRegrProf)

    REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT,raPressLevels(kMaxLayer+1)
    INTEGER :: iDownWard

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: i,j,k,iG,iK,iJ,iI,iRTP,iL1,iProfileLayers,iGasInRTPFile
    REAL :: pProf(kMaxLayer),pTemp(kMaxLayer),MGC,plays(kMaxLayer)
          
    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(160) :: fname

    MGC = kMGC

    ia48layers = 0
    raa48Temp = 0.0
    raa48Press = 0.0

    fname = kRegrFile

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status

    iRTP = 48
    DO iJ = 1, iRTP
    !        write(kStdWarn,*) 'Reading temperature regression profile ',iJ

        status = rtpread(rchan, prof)
        IF (prof%plevs(1) < prof%plevs(prof%nlevs)) THEN
        ! ayers are from TOA to the bottom
            iDownWard = -1
        ELSE
        ! ayers are from GND to the top
            iDownWard = +1
        END IF

        iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
        iProfileLayers = iL1
        ia48layers(iJ) = iL1
        iGasInRTPFile = head.ngas              !!! number of gases

        IF (prof%nlevs > kMaxLayer+1) THEN
            write(kStdErr,*) 'this routine compiled for ',kMaxLayer,' layers'
            write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
            write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
            CALL DoStop
        END IF

        DO i = 1,prof%nlevs
            j = iFindJ(kMaxLayer+1,I,iDownWard)        !!!!notice the kProf+1
            raPressLevels(j) = prof%plevs(i)            !!!!in mb
        END DO

        DO i = 1,prof%nlevs-1
!            j = iFindJ(kMaxLayer+1,I,iDownWard)        !!!!notice the kProf+1
!            pProf(j) = raPressLevels(i) - raPressLevels(i+1)
!            pProf(j) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
            pProf(i) = raPressLevels(i) - raPressLevels(i+1)
            pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
        END DO

    ! now loop (for water only) in the supplied profile
        DO i = 1, prof%nlevs - 1
            j = iFindJ(kMaxLayer,I,iDownWard)
            rT   = prof%ptemp(i)
            plays(i) = (prof%plevs(i)-prof%plevs(i+1))/ &
            log(prof%plevs(i)/prof%plevs(i+1))
            rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
            rP   = plays(i)               !need pressure in mb
            raa48Temp(j,iJ)  = rT
            raa48Press(j,iJ) = rP
        END DO

    !!! then fill bottom of atm with zeros for gas amt, partial pressure
        DO i = prof%nlevs, kMaxLayer
            j = iFindJ(kMaxLayer,I,iDownWard)
            raa48Temp(j,iJ)  = 300.0
            raa48Press(j,iJ) = 1200.0 + i
        END DO

    END DO

    write (kStdWarn,*) 'success : read in 48 regression profiles ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status
    kProfileUnitOpen = -1

    RETURN
    end SUBROUTINE NLTE_RegrTemp

!************************************************************************

! this subroutine finds the closest regression profile
    SUBROUTINE closest_regr_lowertrop(raTempIn, &
    raPressLevels,iNumLayers,iRegr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: raTempIn(kProfLayer)         !!!input kinetic temp profile
    REAL :: raPressLevels(kProfLayer+1)  !!!pressure levels
    INTEGER :: iNumLayers                !!!number of layers
! output
    INTEGER :: iRegr               !!!which regr profile this is closest to

! local vars
    REAL :: raa48Temp(kMaxLayer,kRegrProf)
    REAL :: raa48Press(kMaxLayer,kRegrProf)
    INTEGER :: ia48layers(kRegrProf),iN,i800,i25
    REAL :: raTemp(kProfLayer)

    REAL :: raLayPress(kProfLayer),raX(kMaxLayer),raT(kMaxLayer),rMin,rChi
    INTEGER :: iI,iJ

    raTemp = raTempIn

!!!first find the pressure layering for the kCARTA profile
    DO iI = kProfLayer-iNumLayers+1,kProfLayer
        raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1)
        raLayPress(iI) = raLayPress(iI)/log(raPressLevels(iI)/raPressLevels(iI+1))
    END DO
!!!fill lower layers with some "realistic" increasing values
    DO iI = kProfLayer-iNumLayers,1,-1
        raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) + &
        &                    10*abs(iI-(kProfLayer-iNumLayers+1))
    END DO

!!!now need to flip raLayPress so it is in increaing order!
!!!and do the same flip to raTemp
    DO iI = 1,int(kProfLayer*1.0/2.0)
        rMin = raLayPress(iI)
        raLayPress(iI) = raLayPress(kProfLayer-iI+1)
        raLayPress(kProfLayer-iI+1) = rMin

        rMin = raTemp(iI)
        raTemp(iI) = raTemp(kProfLayer-iI+1)
        raTemp(kProfLayer-iI+1) = rMin
    END DO

!!!! read in the temperature profiles for the 48 regression profiles
    call NLTE_RegrTemp(raa48Temp,raa48Press,ia48layers)
!!!! all 48 profiles should have the same pressure layering
    i800 = 1           !!!this is roughly the 800 mb pressure = 2 km
    i25  = kProfLayer  !!!this is roughly the  25 mb pressure = 25 km
    DO iI = 1,kMaxLayer
        raX(iI) = raa48Press(iI,1)
        IF (raX(iI) >= 800.0) THEN
            i800 = iI
        END IF
        IF (raX(iI) >= 25.0) THEN
            i25 = iI
        END IF
    END DO

!!!! all 48 profiles have the same pressure layering
!!!! so interpolate the kCARTA (pressure,temp) profile on the regression
!!!! profile grid so we can do the comparisons
    call spl(raLayPress,raTemp,kProfLayer,raX,raT,kMaxLayer)

!!!!do the comparisons
    rMin = 1.0e10
    DO iI = 1, 48
        rChi = 0.0
        DO iJ = i800,i25
            rChi = rChi + (raT(iJ)-raa48Temp(iJ,iI))**2
        END DO
    !        print *,iI,rChi
        IF (rChi < rMin) THEN
            iRegr = iI
            rMin = rChi
        END IF
    END DO

    write(kStdWarn,*) 'Closest lower-atm regr profile = ',iRegr
!      print *, 'Closest regr = ',iRegr
!      CALL DoStop

    RETURN
    end SUBROUTINE closest_regr_lowertrop

!************************************************************************
END MODULE rtp_interface_f90
