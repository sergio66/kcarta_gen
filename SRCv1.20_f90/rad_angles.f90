! Copyright 2014
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:42
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the misc angle  routines    ***************
!************************************************************************
!************************************************************************
! viewing geometry angles summary
! define conv = pi/180, Re = 6400 km

!   [zang]=vaconv( sva, salt, alt );
!     ra = Re + alt     = radius at atmosphere height of interest ~ 20 km
!     rs = Re + salt    = radius at satellite orbit               ~ 705 km
!     zang = arcsin ((rs/ra) * sin(conv*sva)) / conv              convert SCANAG (at satellite) to SATZEN (at surface)

!   [zang]=saconv( surfzang, alt );
!     ra = Re + alt     = radius at atmosphere height of interest ~ 20 km
!     rs = Re + salt    = radius at satellite orbit               ~ 705 km
!     zang = arcsin ((re/ra) * sin(conv*surfzang)) / conv          convert SATZEN (at surface) to SCANANG (at satellite)


!        !! test going from surface satzen to sat scanang (at 705 km) back to surface satzen
!      !! should have satzen > scanang
!        print *,'orig satzen(h=0) = 45',saconv_sun(45.0,0.0,705.0),vaconv(saconv_sun(45.0,0.0,705.0),0.0,705.0)
! >>  orig satzen(h=0) = 45   39.54217       45.00000

!        !! test going from sat scanang (at 705 km) to surface satzen (0 km) back to sat scanang (at 705 km)
!      !! should have satzen > scanang
!        print *,'orig scanang(h=705) = 45',vaconv(45.0,0.0,705.0),saconv_sun(vaconv(45.0,0.0,705.0),0.0,705.0)
! >>  orig scanang(h=705) = 45   51.75453       45.00000

!************************************************************************

! this function does the solar viewing angle conversion, by Sergio
! given solar zenith angle == angle wrt nadir at Earth surface,
! can compute solar angle at height h : from satzen(h=0) to scanang(h=H)

! ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
! SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

! directly compare to saconv.m
! directly compare to saconv.m
! directly compare to saconv.m

REAL FUNCTION SACONV_SUN(LSZA, SURFALT, ALT )


NO TYPE, INTENT(IN OUT)                  :: LSZA
NO TYPE, INTENT(IN OUT)                  :: SURFALT
NO TYPE, INTENT(IN OUT)                  :: ALT
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input param
REAL :: LSZA         !! solar/satellite zenith angle at local surface (which is not necessarily 0)
REAL :: SURFALT      !! surface altitude ABOVE EARTH surface (km) (eg [6400 Km + ] aircraft at 05 km)
REAL :: ALT          !! layer altitude   ABOVE EARTH surface (km) (eg [6400 Km + ] layer    at 15 km)

REAL :: rX,rPi,rEarth

!       rPi = 3.1415927
!       rEarth = 6370.0

rPi = kPi
rEarth = kPlanetRadius

!! p = Snell's law = n(i) r(i) sin(theta(i)) = constant =>
!!     (R+h(i)) sin (theta(i)) = constant (assume n(i) = 1)
!! so Re sin(theta_earth) = (Re+h(i))sin (theta(i))

rX = (SURFALT + rEarth)/(ALT + rEarth) * SIN(LSZA * rPi/180.0)
!! SURFALT < ALT ==> ratio in front of sin(LSZA) < 1
!!   ==> local angles get "smaller" the higher you go

IF (rX > 0.99999) rX = 0.99999
IF (rX < 0.00000) rX = 0.00000

SACONV_SUN = ASIN(rX) * 180/rPi

RETURN
END FUNCTION SACONV_SUN

!************************************************************************

! ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
! SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

! directly compare to saconv.m
! directly compare to saconv.m
! directly compare to saconv.m

! this function does the surface --> arb height viewing angle conversion, by Scott
! copied from saconv.f in SARTA package, modified to
!   have input lay height in km
!   have output angle in degrees
! Function to convert the surface solar/satellite zenith angle SZA into the
! local solar/satellite angle at altitude ALT (sunang/scanang)
!    from surface satzen(h=0) to scanang(h=H) at arb altitude H

!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SZA     Surface Zenith Angle        degrees


!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL fun  SACONV  local (height H) angle    degrees

!    Function to convert the Zenith Angle SZA at the Earth's
!    surface into a the local angle at altitude ALT.
!    The local angle generally varies slightly with altitude
!    due to the curvature of the Earth and its atmosphere.
!    The effect is largest at the maximum zenith angle, and
!    disappears as the zenith angle approaches 0 degrees.
!    Currently this function only considers the geometry of the
!    situation, and no refractive effects are included.

!    The layers of the atmosphere may be considered as concentric
!    rings with some average altitude. A ray traced thru these rings
!    at any viewing angle other than nadir will have a slightly
!    different angle (relative to the outward radial at the point
!    of intersection) in each ring.

!    The local angle may be calculated (using The Law of
!    Sines) if we know:
!       The solar/satellire zenith angle, SZA, at the Earth's surface (ALT=0)
!       The layer altitude, ALT.
!       The radius of the Earth, RE.

!    The solution uses the law of sines and sin(180 - x) = sin(x)

REAL FUNCTION ORIG_SACONV_SUN( SZA, ALT )


NO TYPE, INTENT(IN OUT)                  :: SZA
NO TYPE, INTENT(IN)                      :: ALT
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

REAL :: sza,alt

! local variables
REAL :: saconv,conv, re, ra

!      ------------------
!      Assign some values
!      ------------------
!      CONV = pi/180 = degrees to radians conversion factor
CONV=1.7453292E-02

!      RE = radius of the Earth (in km)
!       RE=6.37E+03
RE = kPlanetRadius

!      RA = radius of the point to calc the angle at (in km)
!      Note: layer altitude already in kilometers
RA = rE + ALT

!      -----------------
!      Do the conversion
!      -----------------

SACONV=ASIN( (RE/RA) * SIN(CONV*SZA) )

! change back to degrees
ORIG_SACONV_SUN = SACONV/conv

RETURN
END FUNCTION ORIG_SACONV_SUN

!************************************************************************
!    Function to convert the AIRS satellite viewing angle into the
!    local path angle (      Viewing Angle CONVersion   )
! for downward looking instrument : from scanang(h=H) to satzen(h=0)

! directly compare to vaconv.m
! directly compare to vaconv.m
! directly compare to vaconv.m
! can Matlab debug version, but WARNING its input arguments are switched thusly
!   [zang]=vaconv( sva, salt, alt );
! cd /asl/matlab2012/science
! >> alt = [0:10:80]
! >> [zang]=vaconv(50, 705*1000, alt'*1000 )

! uses law of sines

!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SVA     Satellite viewing angle     degrees
!    REAL      SALT    Satellite altitude          km

REAL FUNCTION VACONV( SVA, ALT, SALT )


REAL, INTENT(IN OUT)                     :: SVA
REAL, INTENT(IN)                         :: ALT
REAL, INTENT(IN)                         :: SALT
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'



!      LOCAL VARIABLES
REAL :: CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink

!      CONV = pi/180 = degrees to radians conversion factor
CONV = 1.7453292E-02
theta = sva*conv

!      RE = radius of the Earth (in km)
!       RE=6.37E+03
RE = kPlanetRadius

!      RA = radius of the point to calc the angle at (in km)
RA = rE + ALT

!      RS = radius of the satellite orbit (in km)
RS = rE + SALT


!      -----------------
!      Do the conversion
!      -----------------

rBoink =  (RS/RA) * SIN(CONV*SVA)
IF (ABS(rBoink) > 1.0) rBoink = 0.9999

RTEMP=ASIN(rBoink)

! change back to degrees
VACONV = RTEMP/conv

RETURN
END FUNCTION VACONV

!************************************************************************
!    Function to convert the AIRS satellite viewing angle into the
!    local path angle (      Viewing Angle CONVersion   ) WITH RefIndex
! for downward looking instrument

! directly compare to vaconv.m, but this uses Snell
! directly compare to vaconv.m, but this uses Snell
! directly compare to vaconv.m, but this uses Snell
! can Matlab debug version, but WARNING its input arguments are switched thusly
!   [zang]=vaconv( sva, salt, alt );
! cd /asl/matlab2012/science
! >> alt = [0:10:80]
! >> [zang]=vaconv(50, 705*1000, alt'*1000 )

! uses law of sines


!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SVA     Satellite viewing angle     degrees
!    REAL      SALT    Satellite altitude          km

REAL FUNCTION VACONV_Snell( SVA, ALT, SALT, rNumberDensity )


REAL, INTENT(IN OUT)                     :: SVA
REAL, INTENT(IN)                         :: ALT
REAL, INTENT(IN)                         :: SALT
NO TYPE, INTENT(IN OUT)                  :: rNumberDen
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL :: rNumberDensity

!      LOCAL VARIABLES
REAL :: CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink
REAL :: lambda,n0,gamma,ref_ind,mu
REAL :: rLosch             !Loschmidt number

rLosch = (kAtm2mb*100/kBoltzmann/273 * 1E-6)    ! (p0/k T0)

!      CONV = pi/180 = degrees to radians conversion factor
CONV = 1.7453292E-02
theta = sva*conv

!      RE = radius of the Earth (in km)
RE = 6367.512

! http://mintaka.sdsu.edu/GF/explain/atmos_refr/horizon.html
! http://www.ess.uci.edu/~cmclinden/link/xx/node45.html for TONS of detail
!      RA = radius of the point to calc the angle at (in km)
RA = rE + ALT

!      RS = radius of the satellite orbit (in km)
RS = rE + SALT

! ref ind of air at 10 um
lambda = 10
n0 = 64.328 + 29498.1/(146-1/lambda/lambda) + 255.4/(41 - 1/lambda/lambda)
n0 = 1 + n0/1E6
gamma = (n0*n0-1)/(n0*n0+2)

IF (ALT >= 200) THEN
  ref_ind = 1.0
ELSE
  mu = rNumberDensity/rLosch * gamma
  ref_ind = SQRT((2*mu+1)/(1-mu))
END IF


!      -----------------
!      Do the conversion
!      -----------------

rBoink =  ABS((RS/RA) * SIN(CONV*SVA)/ref_ind)
IF (ABS(rBoink) > 1.0) rBoink = 0.9999

RTEMP = ASIN(rBoink)

! change back to degrees
VACONV_SNELL = RTEMP/conv

RETURN
END FUNCTION VACONV_Snell

!************************************************************************
! this subroutine finds the layer dependent satellite viewing angle

SUBROUTINE FindSunLayerAngles(rSatHeight,rSurfHeight,iAtm,iaNumLayer,iaaRadLayer,  &
    raLayHgt,rPrBdry1,rPrBdry2,rSunAngle,raSunAngles)


REAL, INTENT(IN OUT)                     :: rSatHeight
NO TYPE, INTENT(IN OUT)                  :: rSurfHeigh
INTEGER, INTENT(IN OUT)                  :: iAtm
NO TYPE, INTENT(IN OUT)                  :: iaNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN)                         :: raLayHgt(kProfLayer)
REAL, INTENT(IN OUT)                     :: rPrBdry1
REAL, INTENT(IN OUT)                     :: rPrBdry2
REAL, INTENT(IN)                         :: rSunAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rSatHeight = tells us if we need to do ray tracing (> 0) or not (< 0)
! raLayHgt   = height of individual layers, read in from profile, in meters
! rPrBdry1   = start pressure boundary
! rPrBdry2   = stop pressure boundary
! rSunAngle  = sun angle at TOA
! raSunAngles = layer dependent sun viewing angle
REAL :: rSurfHeight
REAL :: raSunAngles(kProfLayer)
INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: salt,saconv_sun,rSunHeight
INTEGER :: iI,iX,iMin,iMax,iaRadLayer(kProfLayer)

rSunHeight = 150.0E9  !! in km

! as default, set all angles to be the sun view angle
DO iI=1,kProfLayer
  raSunAngles(iI) = rSunAngle
END DO

iMin = +1000000
iMax = -1000000
DO iI=1,iaNumlayer(iAtm)
  iaRadLayer(iI) = iaaRadLayer(iAtm,iI)
  IF (iaRadLayer(iI) > iMax) iMax = iaRadLayer(iI)
  IF (iaRadLayer(iI) < iMin) iMin = iaRadLayer(iI)
END DO

! *********************** rSatHeight, raLayHgt is in METERS ***************

IF ((rSatHeight > 0.0) .AND. (ABS(rSunAngle) > 1.0E-4)) THEN
!have to find layer dependent angles
  IF (rPrBdry1 > rPrBdry2) THEN !downward looking instr
    DO iI=1,kProfLayer
      IF (rSatHeight > raLayHgt(iI)) THEN
        raSunAngles(iI) = saconv_sun(ABS(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
        IF (rSunAngle < 0.0) raSunAngles(iI) = -raSunAngles(iI)
        IF (kOuterLoop == 1) THEN
          IF (iI >= iMin .AND. iX <= iMax) THEN
            iX = (iI - iMin + 1)
          ELSE
            iX = -1
          END IF
          IF (iX == 1) WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
          WRITE(kStdWarn,*)'sun: lay#/rad# lay/surfhgt, localzen/traced angle ',  &
              iI,iX,raLayHgt(iI)/1000,rSurfHeight/1000,rSunAngle,raSunAngles(iI)
        END IF
      END IF
    END DO
  END IF
  IF (rPrBdry2 > rPrBdry1) THEN !upward looking instr
    salt=705.00
    DO iI=1,kProfLayer
!            write(kStdWarn,*)'sun up lay hgt ',iI,raLayHgt(iI)/1000
      IF (rSatHeight > raLayHgt(iI)) THEN
        raSunAngles(iI) = saconv_sun(ABS(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
        IF (rSunAngle < 0.0) raSunAngles(iI) = -raSunAngles(iI)
        IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
              (ABS(kLongOrShort) <= 1)) THEN
          WRITE(kStdWarn,*)'sun up lay hgt, orig/traced angle ',  &
              iI,raLayHgt(iI)/1000,rSunAngle,raSunAngles(iI)
        END IF
      END IF
    END DO
  END IF
END IF

!      call dostop

RETURN
END SUBROUTINE FindSunLayerAngles

!************************************************************************
! this subroutine finds the layer dependent satellite viewing angle

SUBROUTINE FindLayerAngles(rSatHeight,raLayHgt,  &
    rPrBdry1,rPrBdry2,rSatAngle,raLayAngles,  &
    iAtm,iaaRadLayer,iaNumlayer,raNumberDensity)


REAL, INTENT(IN)                         :: rSatHeight
REAL, INTENT(IN)                         :: raLayHgt(kProfLayer)
REAL, INTENT(IN OUT)                     :: rPrBdry1
REAL, INTENT(IN OUT)                     :: rPrBdry2
REAL, INTENT(IN)                         :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
INTEGER, INTENT(IN OUT)                  :: iAtm
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
NO TYPE, INTENT(IN)                      :: iaNumlayer
NO TYPE, INTENT(IN OUT)                  :: raNumberDe
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rSatHeight = satellite height in meters
!   (from spec : http://asl.umbc.edu/pub/motteler/hdf/rtpV105/doc/rtpspec.html)
! [sergio@tara-fe1 WORK]$ grep -in zobs   /asl/packages/sartaV108/Src/*.f
! /asl/packages/sartaV108/Src/rdrtp_df.f:329:       ZSAT=PROF.zobs/1000
!  /asl/packages/sartaV108/Src/rdrtp.f:328:       ZSAT=PROF.zobs/1000

! raLayHgt = height of individual layers, read in from profile, in meters

! rPrBdry1 = start pressure boundary
! rPrBdry2= stop pressure boundary
! rSatAngle = satellite viewing angle (scanang)
! raLayAngles = layer dependent satellite viewing angle (satzen)

REAL :: raLayAngles(kProfLayer)
INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raNumberDensity(kProfLayer)

REAL :: vaconv,vaconv_Snell,saconv_sun
REAL :: raLayAnglesSnell(kProfLayer),raLayAnglesNoSnell(kProfLayer)
INTEGER :: iI,iaRadLayer(kProfLayer),iMin,iMax,iX,iDefault,iUseSnell

! as default, set all angles to be the satellite view angle
DO iI=1,kProfLayer
  raLayAngles(iI)      = rSatAngle
  raLayAnglesSnell(iI) = rSatAngle
END DO

iDefault = -1   !! no  Snell, yes layer curvature, similar to SARTA

iUseSnell = +1  !! yes Snell, yes layer curvature
iUseSnell = -1  !! no  Snell, yes layer curvature, similar to SARTA
iUseSnell = 0   !! no  Snell, no  layer curvature

iUseSnell = iaaOverrideDefault(2,7)
IF (ABS(iUseSnell) > 1) THEN
  WRITE(kStdErr,*) 'invalid iUseSnell ',iUseSnell
  CALL DoStop
END IF
IF ((iDefault /= iUseSnell) .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdWarn,*) 'not/using Snell law in FindLayerAngles (raytrace thru layers)',iUseSnell
  WRITE(kStdErr,*) 'not/using Snell law in FindLayerAngles (raytrace thru layers)',iUseSnell
END IF

iMin = +1000000
iMax = -1000000
DO iI = 1,iaNumlayer(iAtm)
  iaRadLayer(iI) = iaaRadLayer(iAtm,iI)
  IF (iaRadLayer(iI) > iMax) iMax = iaRadLayer(iI)
  IF (iaRadLayer(iI) < iMin) iMin = iaRadLayer(iI)
END DO

! *********** vaconv assumes layheight,satheight is in KM ***************
! *********** SatHeight is in METERS (from RTP), and raLayHgt is in METERS ****

! orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)
! orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)
! orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)

IF ((rSatHeight > 0.0) .AND. (ABS(rSatAngle) > 1.0E-4) .AND. (ABS(iUseSnell) == 0)) THEN
  DO iI=1,kProfLayer
    raLayAnglesNoSnell(iI) = ABS(rSatAngle)
    raLayAnglesSnell(iI) = ABS(rSatAngle)
    IF (kOuterLoop == 1) THEN
      IF (iI >= iMin .AND. iX <= iMax) THEN
        iX = (iI - iMin + 1)
      ELSE
        iX = -1
      END IF
      IF (iX == 1) THEN
        WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
        WRITE(kStdWarn,*) 'Downlook Instr'
        WRITE(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
        WRITE(kStdWarn,*) 'TOA scanang, GND satzen = ',ABS(rSatAngle),vaconv(ABS(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
        WRITE(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
      END IF
      WRITE(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI),  &
          1.0/COS(raLayAngles(iI)*kPi/180.0)
    END IF      !! if kOuterLoop .EQ. 1
    
  END DO
ELSE IF ((rSatHeight > 0.0) .AND. (ABS(rSatAngle) > 1.0E-4) .AND. (ABS(iUseSnell) == 1)) THEN
!have to find layer dependent angles
  IF (rPrBdry1 > rPrBdry2) THEN !downward looking instr
    DO iI=1,kProfLayer
!            print *,iI,rSatHeight,raLayHgt(iI)
      IF (rSatHeight > raLayHgt(iI)) THEN
        raLayAnglesNoSnell(iI) = vaconv(ABS(rSatAngle),raLayHgt(iI)/1000,  &
            rSatHeight/1000)
        raLayAnglesSnell(iI) = vaconv_Snell(ABS(rSatAngle),raLayHgt(iI)/1000,  &
            rSatHeight/1000,raNumberDensity(iI))
! diff between Snell and noSnell is less than 2e-2
!              print *,iI,raLayAnglesSnell(iI),raLayAnglesNoSnell(iI),raLayAnglesNoSnell(iI)-raLayAnglesSnell(iI)
! but Scott/SARTA uses NoSnell
        IF (iUseSnell == 1) THEN
          raLayAngles(iI) = raLayAnglesSnell(iI)    !!! sergio
        ELSE IF (iUseSnell == -1) THEN
          raLayAngles(iI) = raLayAnglesNoSnell(iI)      !!! scott
        END IF
        IF (rSatAngle < 0.0) raLayAngles(iI) = -raLayAngles(iI)
        IF (kOuterLoop == 1) THEN
          IF (iI >= iMin .AND. iX <= iMax) THEN
            iX = (iI - iMin + 1)
          ELSE
            iX = -1
          END IF
          IF (iX == 1) THEN
            WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
            WRITE(kStdWarn,*) 'Downlook Instr'
            WRITE(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
            WRITE(kStdWarn,*) 'TOA scanang, GND satzen = ',ABS(rSatAngle),vaconv(ABS(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
            WRITE(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
          END IF
          WRITE(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI),  &
              1.0/COS(raLayAngles(iI)*kPi/180.0)
        END IF      !! if kOuterLoop .EQ. 1
      END IF
    END DO
  ELSE
!no need to do anything much, the angles are so close to nadir
    iX = -1
    WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
    WRITE(kStdWarn,*) 'Downlook Instr'
    WRITE(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
    WRITE(kStdWarn,*) 'TOA scanang, GND satzen = ',ABS(rSatAngle),vaconv(ABS(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
    WRITE(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
    DO iI=1,kProfLayer
      WRITE(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI),1.0
    END DO
  END IF
  
  IF ((rPrBdry2 > rPrBdry1) .AND. (ABS(iUseSnell) == 0)) THEN !upward looking instr
    DO iI=1,kProfLayer
      raLayAnglesNoSnell(iI) = ABS(rSatAngle)
      raLayAngles(iI) = raLayAnglesSnell(iI)
      IF (kOuterLoop == 1) THEN
        IF (iI >= iMin .AND. iX <= iMax) THEN
          iX = (iI - iMin + 1)
        ELSE
          iX = -1
        END IF
        IF (iX == 1) THEN
          WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
          WRITE(kStdWarn,*)'up : lay#/rad# lay/sat hgt, local satzen/layer iI angle'
        END IF
        WRITE(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
      END IF
    END DO
    
  ELSE IF ((rPrBdry2 > rPrBdry1) .AND. (ABS(iUseSnell) == 1)) THEN !upward looking instr
!          print *,'FindLayerAngles',iMin,iMax,iAtm,rSatHeight,iaNumLayer(iAtm),
!     $               raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000
    DO iI=1,kProfLayer
      IF (rSatHeight > raLayHgt(iI)) THEN
        raLayAnglesNoSnell(iI) = saconv_sun(ABS(rSatAngle),  &
            raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000, raLayHgt(iI)/1000)
        raLayAngles(iI) = raLayAnglesSnell(iI)
        raLayAngles(iI) = raLayAnglesNoSnell(iI)
        IF (rSatAngle < 0.0) raLayAngles(iI) = -raLayAngles(iI)
        IF (kOuterLoop == 1) THEN
          IF (iI >= iMin .AND. iX <= iMax) THEN
            iX = (iI - iMin + 1)
          ELSE
            iX = -1
          END IF
          IF (iX == 1) THEN
            WRITE(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
            WRITE(kStdWarn,*)'up : lay#/rad# lay/sat hgt, local satzen/layer iI angle'
          END IF
          WRITE(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
        END IF
      END IF
    END DO
  END IF
END IF


999  FORMAT(I3,' ',I3,'  ',5(F10.4,' '))

RETURN
END SUBROUTINE FindLayerAngles

!************************************************************************
! this subroutine computes the gauss-legendre abscissa weights and points
! from Numerical Recipes
! nn is the number of points <= kProfLayer, x,w are output abscissa and wts

SUBROUTINE FindGauss(nn,daX,daW)


INTEGER, INTENT(IN)                      :: nn
DOUBLE PRECISION, INTENT(OUT)            :: daX(kGauss)
DOUBLE PRECISION, INTENT(OUT)            :: daW(kGauss)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'




DOUBLE PRECISION :: daX1(2*kGauss),daW1(2*kGauss)
DOUBLE PRECISION :: x1,x2
INTEGER :: m,j,i,n
DOUBLE PRECISION :: z1,z,xm,xl,pp,p1,p2,p3,epss

epss = 3.0E-11

x1 = -1.0D0
x2 = +1.0D0

IF ((nn > kGauss) .OR. (nn < 0)) THEN
  WRITE (kStdErr,*) 'need 0 < nn <= kGauss'
  CALL DoStop
END IF

n=nn*2

IF (MOD(n,2) == 1) THEN
  WRITE (kStdErr,*) 'need n to be even'
  CALL DoSTOP
END IF

IF (x2 < x1) THEN
  xm = x1
  x1 = x2
  x2 = xm
END IF

m  = (n+1)/2

m  = n/2
xm = 0.5*(x2+x1)
xl = 0.5*(x2-x1)

DO i = 1,m                    !loop over desired roots
  z = COS(kPi*(i-0.25)/(n+0.5))
  20     CONTINUE
  p1 = 1.0
  p2 = 0.0
  DO j = 1,n
    p3 = p2
    p2 = p1
    p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
  END DO
  pp = n*(z*p1-p2)/(z*z-1)
  z1 = z
  z  = z1-p1/pp
  IF (ABS(z-z1) > epss) THEN
    GO TO 20
  END IF
  
  daX(i)     = xm+xl*z
  daW(i)     = 2*xl/((1-z*z)*pp*pp)
  
  daX1(i)     = xm-xl*z
  daX1(n+1-i) = xm+xl*z
  daW1(i)     = 2*xl/((1-z*z)*pp*pp)
  daW1(n+1-i) = daW(i)
END DO

RETURN
END SUBROUTINE FindGauss

!************************************************************************
! this subroutine reads in the first moment (l=1) quadrature points and weights
! using results in Table 1, second column (m=1) of
! "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
! Journal of Atmospheric Sciences, v57, pg 753 (2000)

SUBROUTINE FindGauss2old(nn,daX,daW)


INTEGER, INTENT(OUT)                     :: nn
DOUBLE PRECISION, INTENT(OUT)            :: daX(kGauss)
DOUBLE PRECISION, INTENT(OUT)            :: daW(kGauss)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'




IF (nn == 1) THEN
  daX(1) = 2.0/3.0
  daW(1) = 1.0/2.0
ELSE IF (nn == 2) THEN
  daX(1) = 0.3550510
  daX(2) = 0.8449489
  daW(1) = 0.1819856
  daW(2) = 0.3180414
ELSE IF (nn == 3) THEN
  daX(1) = 0.2123405
  daX(2) = 0.5905331
  daX(3) = 0.9114120
  daW(1) = 0.0698270
  daW(2) = 0.2292411
  daW(3) = 0.2009319
ELSE IF (nn == 4) THEN
  daX(1) = 0.1397599
  daX(2) = 0.4164096
  daX(3) = 0.7231570
  daX(4) = 0.9428958
  daW(1) = 0.0311810
  daW(2) = 0.1298475
  daW(3) = 0.2034646
  daW(4) = 0.1355069
ELSE
  WRITE(kStdErr,*) 'FindGauss2 : need nn = 1,2,3 or 4 not ',nn
  CALL DoStop
END IF

RETURN
END SUBROUTINE FindGauss2old

!************************************************************************
! this subroutine reads in the first moment (l=1) quadrature points and weights
! using results in Table 1, second column (m=1) of
! "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
! Journal of Atmospheric Sciences, v57, pg 753 (2000)
! also see RRTM code v3.3 rtreg.ffor more decimal points

SUBROUTINE FindGauss2(nn,daX,daW)


INTEGER, INTENT(OUT)                     :: nn
DOUBLE PRECISION, INTENT(OUT)            :: daX(kGauss)
DOUBLE PRECISION, INTENT(OUT)            :: daW(kGauss)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: ii


IF (nn == 1) THEN
  daX(1) = 3.0/2.0
  daW(1) = 1.0/2.0
ELSE IF (nn == 2) THEN
  daX(1) = 2.81649655
  daX(2) = 1.18350343
  daW(1) = 0.1819586183
  daW(2) = 0.3180413817
ELSE IF (nn == 3) THEN
  daX(1) = 4.70941630
  daX(2) = 1.69338507
  daX(3) = 1.09719858
  daW(1) = 0.0698269799
  daW(2) = 0.2292411064
  daW(3) = 0.2009319137
ELSE IF (nn == 4) THEN
  daX(1) = 7.15513024
  daX(2) = 2.40148179
  daX(3) = 1.38282560
  daX(4) = 1.06056257
  daW(1) = 0.0311809710
  daW(2) = 0.1298475476
  daW(3) = 0.2034645680
  daW(4) = 0.1355069134
ELSE
  WRITE(kStdErr,*) 'FindGauss2 : need nn = 1,2,3 or 4 not ',nn
  CALL DoStop
END IF

DO ii = 1,nn
!change from secant to cosine
  daX(ii) = 1.0/daX(ii)
END DO

RETURN
END SUBROUTINE FindGauss2

!************************************************************************
! this function does the expintegral raY = E3(raX), using recursion

SUBROUTINE expint3(raX,raY)


REAL, INTENT(OUT)                        :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI

CALL expint(raX,raY)               !do y = expint(x)
DO iI = 1,kMaxPts
  IF (raX(iI) < 1.0E-28) THEN
    raY(iI) = 0.5
  ELSE IF (raX(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint3 iI,raX(iI) = ',iI,raX(iI)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE
    raY(iI) = (EXP(-raX(iI))*(1-raX(iI)) + (raX(iI)**2)*raY(iI))/2.0
  END IF
END DO

RETURN
END SUBROUTINE expint3

!************************************************************************
! this function does the expintegral raY = E2(raX), using recursion

SUBROUTINE expint2(raX,raY)


REAL, INTENT(OUT)                        :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI

CALL expint(raX,raY)               !do y = expint(x)
DO iI = 1,kMaxPts
  IF (raX(iI) < 1.0E-28) THEN
    raY(iI) = 1.0
  ELSE IF (raX(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint2 iI,raX(iI) = ',iI,raX(iI)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE
    raY(iI) = EXP(-raX(iI)) - raX(iI)*raY(iI)
  END IF
END DO

RETURN
END SUBROUTINE expint2

!************************************************************************
! this function does the expintegral raY = E1(raX), ala Matlab
! really, this can be found in any Math handbook, or on the web
! assumes input raX has no elements == 0, and all are real

SUBROUTINE expint(raX,raY)


REAL, INTENT(IN)                         :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI,iJ,iK,iFr,iFound,iN
INTEGER :: iaPositive(kMaxPts),iPos,iaNegative(kMaxPts),iNeg
DOUBLE PRECISION :: p(9),egamma,daX(kMaxPts),daY(kMaxPts),DIMAG
COMPLEX*16 daZ(kMaxPts)
DOUBLE PRECISION :: daPterm(kMaxPts),daTerm(kMaxPts)
DOUBLE PRECISION :: am2(kMaxPts),am1(kMaxPts),bm2(kMaxPts),bm1(kMaxPts)
DOUBLE PRECISION :: f(kMaxPts),oldf(kMaxPts),alpha,beta,a,b

DATA (P(iJ),iJ=1,9)/  &
    -3.602693626336023D-09, -4.819538452140960D-07, -2.569498322115933D-05,  &
    -6.973790859534190D-04, -1.019573529845792D-02, -7.811863559248197D-02,  &
    -3.012432892762715D-01, -7.773807325735529D-01,  8.267661952366478D+00/

DO iFr = 1,kMaxPts
  daY(iFr) = 0.0D0
  DO iJ = 1,9
    iI = 9-iJ
    daY(iFr) = daY(iFr) + p(iJ)*((raX(iFr)*1.0D0)**iI)
  END DO
END DO

! polyv = polyval(p,real(x))
! k = find( abs(imag(x)) <= polyv )   this is iPos, except imag(x) == 0 always
iPos = 0
iNeg = 0
DIMAG = 0.0D0
DO iFr = 1,kMaxPts
  IF (DIMAG <= daY(iFr)) THEN
    iPos = iPos + 1
    iaPositive(iPos) = iFr
  ELSE
    iNeg = iNeg + 1
    iaNegative(iNeg) = iFr
  END IF
END DO

! ---------- these are for x >= 0 ------------------------------------------
IF (iPos >= 1) THEN
  egamma = 5.7721566490153286061D-1
  DO iK = 1,iPos
    iFr = iaPositive(iK)
    daX(iFr) = raX(iFr)*1.0D0
    daZ(iFr) = -egamma - CDLOG(DCMPLX(daX(iFr),0.0D0))
    daPterm(iFr) = daX(iFr)
    daterm(iFr)  = daX(iFr)
  END DO
  
  iJ = 1
  
  10     CONTINUE
  iJ = iJ + 1
  
  DO iK = 1,iPos
    iFr = iaPositive(iK)
    daZ(iFr) = daZ(iFr) + daTerm(iFr)
    daPterm(iFr) = -daX(iFr) * daPterm(iFr)/iJ
    daTerm(iFr)  = daPterm(iFr)/iJ
  END DO
  
! check for convergence
  iFound = -1
  DO iK = 1,iPos
    iFr = iaPositive(iK)
    IF (ABS(daTerm(iFr)) > 1.0D-16*CDABS(daZ(iFr))) THEN
      iFound = +1
      GO TO 11
    END IF
  END DO
  11      CONTINUE
  
  IF (iFound > 0) THEN
    GO TO 10
  END IF
  
  DO iK = 1,iPos
    iFr = iaPositive(iK)
    raY(iFr) = SNGL(dreal(daZ(iFr)))
  END DO
END IF

! ---------- these are for x < 0  ------------------------------------------
IF (iNeg >= 1) THEN
  iN = 1  !calculating E1(x)
  
  DO iK = 1,iNeg
    iFr = iaNegative(iK)
    daX(iFr)  = raX(iFr)*1.0D0
    am2(iFr)  = 0.0D0
    bm2(iFr)  = 1.0D0
    am1(iFr)  = 1.0D0
    bm1(iFr)  = daX(iFr)
    f(iFr)    = am1(iFr)/bm1(iFr)
    oldf(iFr) = 1.0D100
  END DO
  
  iJ = 2
  20     CONTINUE
  
! calculate the coefficients of the recursion formulas for j even
  alpha = iN - 1+(iJ/2) ! note: beta= 1
  
  DO iK = 1,iNeg
    iFr = iaNegative(iK)
    
!calculate A(j), B(j), and f(j)
    a = am1(iFr) + alpha * am2(iFr)
    b = bm1(iFr) + alpha * bm2(iFr)
    
! save new normalized variables for next pass through the loop
!  note: normalization to avoid overflow or underflow
    am2(iFr) = am1(iFr) / b
    bm2(iFr) = bm1(iFr) / b
    am1(iFr) = a / b
    bm1(iFr) = 1.0D0
    
    oldf(iFr) =f(iFr)
    f(iFr) = am1(iFr)
  END DO
  
  iJ = iJ+1
  
! calculate the coefficients for j odd
  alpha = (iJ-1)/2
  DO iK = 1,iNeg
    iFr = iaNegative(iK)
    beta = daX(iFr)
    a = beta * am1(iFr) + alpha * am2(iFr)
    b = beta * bm1(iFr) + alpha * bm2(iFr)
    am2(iFr) = am1(iFr) / b
    bm2(iFr) = bm1(iFr) / b
    am1(iFr) = a / b
    bm1(iFr) = 1.0D0
    oldf(iFr) = f(iFr)
    f(iFr) = am1(iFr)
  END DO
  
  iJ = iJ+1
  
! check for convergence
  iFound = -1
  DO iK = 1,iNeg
    iFr = iaNegative(iK)
    IF (ABS(f(iFr)-oldf(iFr)) > 1.0D-14*ABS(f(iFr))) THEN
      iFound = +1
      GO TO 21
    END IF
  END DO
  21     CONTINUE
  
  IF (iFound > 0) THEN
    GO TO 20
  END IF
  
  DO iK = 1,iNeg
    iFr = iaNegative(iK)
    daY(iFr) =  EXP(-daX(iFr)) * f(iFr)
!!! daY is really complex, but ignore the imag part :
!!!     - i*pi*((real(xk)<0)&(imag(xk)==0))
    raY(iFr) = SNGL(daY(iFr))
  END DO
  
END IF

RETURN
END SUBROUTINE expint

!************************************************************************
! calls function to do expintegral raY = E(n,raX), ala Numerical Recipes
! this should be faster than the Matlab version

SUBROUTINE expintfast3(raX,raY)


REAL, INTENT(OUT)                        :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI
REAL :: fastExpint

DO iI = 1,kMaxPts
  IF (raX(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint iI,raX(iI) = ',iI,raX(iI)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE
    raY(iI) = fastExpint(3,raX(iI))
  END IF
END DO

RETURN
END SUBROUTINE expintfast3

!************************************************************************
! this is fast exponential integrals, from Numerical recipes
! called by SUBROUTINE expintfast3(n,raX)

REAL FUNCTION fastExpint(n,x)


INTEGER, INTENT(IN)                      :: n
REAL, INTENT(IN)                         :: x
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars



! local vars
INTEGER :: MaxIT
REAL :: Euler,FPmin,Eps

INTEGER :: i,ii,nm1
REAL :: r,fact,del,psi,a,b,c,d,h

r      = 0.0
MaxIT = 100
Euler = 0.5772156649
FPMin = 1.0E-30
Eps   = 1.0E-7

nm1 = n-1
IF (n == 0) THEN
  r = EXP(-x)/x
ELSE IF (x <= Eps) THEN
  r = 1.0/nm1
ELSE IF (x > 1.0) THEN
  b = x + n
  c = 1/FPmin
  d = 1.0/b
  h = d
  DO i = 1,MAXIT
    a = -i*(nm1+i)
    b = b + 2.0
    d = 1/(a*d+b)
    c = b + a/c
    del = c*d
    h = h*del
    IF (ABS(del-1) <= eps) THEN
      r = h*EXP(-x)
      GO TO 10
    END IF
  END DO
  
ELSE IF ((x >= Eps) .AND. (x <= 1.0)) THEN
!set first term ... wow
  IF (nm1 /= 0) THEN
    r = 1.0/nm1
  ELSE
    r = -LOG(x) - Euler
  END IF
  fact = 1.0
  DO i = 1,MaxIT
    fact = fact * (-x/i)
    IF (i /= nm1) THEN
      del = -fact/(i-nm1)
    ELSE
      psi = -euler
      DO ii = 1,nm1
        psi = psi + 1.0/ii
      END DO
      del = fact*(-LOG(x) + psi)
    END IF
    r = r + del
    IF (ABS(del) < ABS(r)*Eps) GO TO 10
  END DO
END IF

10   CONTINUE
fastExpint = r

RETURN
END FUNCTION fastExpint

!************************************************************************
! does expintegral raY = E3(raX), using fifth order polynom fits to expint3

SUBROUTINE expintsuperfast3_6(raX,raY)


REAL, INTENT(OUT)                        :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI,i,j,iWhich
REAL :: fastExpint,raaPout(5,6)

!!this is for x > 0.00 and x <= 0.10
DATA ((raaPout(i,j),j=1,6),i=1,1)/  &
    4.999997490759554E-01, -9.952726790169388E-01,  2.469512034408480E+00,  &
    -1.666132334280038E+01,  1.120106392322948E+02, -3.420299925382391E+02/

!!this is for x > 0.10 and x <= 1.00
DATA ((raaPout(i,j),j=1,6),i=2,2)/  &
    4.964930888860016E-01,  -9.002236958326287E-01,  1.058882468914909E+00,  &
    -9.590227134465121E-01,   5.595049802866501E-01, -1.460130149043994E-01/

!!this is for x > 1.00 and x <= 6.00
DATA ((raaPout(i,j),j=1,6),i=3,3)/  &
    9.348983139232728E-03,   8.153540216831002E-01, -1.318262705075364E+00,  &
    1.507590471592263E+00,  -9.895617788618024E-01,  2.739113519528383E-01/

!!this is for x > 6.00 and x <= 10.00
DATA ((raaPout(i,j),j=1,6),i=4,4)/  &
    3.700147480617722E-04,   9.814904610819430E-01, -2.593611383300170E+00,  &
    6.727199620937857E+00,  -1.257888244164997E+01,  1.144472478235488E+01/

!!this is for x > 10.00 and x <= 20.00
DATA ((raaPout(i,j),j=1,6),i=5,5)/  &
    3.360502391024836E-05,  9.969795018873142E-01, -2.884279567510683E+00,  &
    9.510785069482322E+00, -2.619205947978609E+01,  3.863540676528253E+01/

DO iI = 1,kMaxPts
  raY(iI) = 0.0
  IF (raX(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint iI,raX(iI) = ',iI,raX(iI)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE IF (raX(iI) <= 0.1) THEN
    iWhich = 1
    DO j = 1,6
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raX(iI)**(j-1))
    END DO
  ELSE IF (raX(iI) <= 1.0) THEN
    iWhich = 2
    DO j = 1,6
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raX(iI)**(j-1))
    END DO
  ELSE IF (raX(iI) <= 6.0) THEN
    iWhich = 3
    DO j = 1,6
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE IF (raX(iI) <= 10.0) THEN
    iWhich = 4
    DO j = 1,6
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE IF (raX(iI) <= 20.0) THEN
    iWhich = 5
    DO j = 1,6
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE !!   expint3(x) <= 1e-10 for x >= 20
    iWhich = 6
    raY(iI) = 0.0
  END IF
END DO

RETURN
END SUBROUTINE expintsuperfast3_6

!************************************************************************
! does expintegral raY = E3(raX), using 2 - 5 order polynom fits to expint3

SUBROUTINE expintsuperfast3(raX,raY)


REAL, INTENT(OUT)                        :: raX(kMaxPts)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

! output vars


! local vars
INTEGER :: iI,i,j,iWhich,iaMax(5),iCnt
REAL :: fastExpint,raaPout(5,6)

DATA (iaMax(i),i=1,5) /3, 6, 6, 3, 3/

!!this is for x > 0.00 and x <= 0.10
DATA ((raaPout(i,j),j=1,3),i=1,1)/  &
    4.999958864847588E-01, -9.704714999886551E-01,  1.357944294418206E+00/

!!this is for x > 0.10 and x <= 1.00
DATA ((raaPout(i,j),j=1,6),i=2,2)/  &
    4.964930888860016E-01,  -9.002236958326287E-01,  1.058882468914909E+00,  &
    -9.590227134465121E-01,   5.595049802866501E-01, -1.460130149043994E-01/

!!this is for x > 1.00 and x <= 6.00
DATA ((raaPout(i,j),j=1,6),i=3,3)/  &
    9.348983139232728E-03,   8.153540216831002E-01, -1.318262705075364E+00,  &
    1.507590471592263E+00,  -9.895617788618024E-01,  2.739113519528383E-01/

!!this is for x > 6.00 and x <= 10.00
DATA ((raaPout(i,j),j=1,3),i=4,4)/  &
    6.850600578733759E-03,   8.121694943009118E-01, -9.875485935560442E-01/

!!this is for x > 10.00 and x <= 20.00
DATA ((raaPout(i,j),j=1,3),i=5,5)/  &
    1.873706336320716E-03,   9.115323790575932E-01, -1.489123902774130E+00/

DO iI = 1,kMaxPts
  raY(iI) = 0.0
  IF (raX(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint iI,raX(iI) = ',iI,raX(iI)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE IF (raX(iI) <= 0.1) THEN
    iWhich = 1
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raX(iI)**(j-1))
    END DO
  ELSE IF (raX(iI) <= 1.0) THEN
    iWhich = 2
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raX(iI)**(j-1))
    END DO
  ELSE IF (raX(iI) <= 6.0) THEN
    iWhich = 3
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE IF (raX(iI) <= 10.0) THEN
    iWhich = 4
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE IF (raX(iI) <= 20.0) THEN
    iWhich = 5
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raX(iI))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raX(iI))
  ELSE !!   expint3(x) <= 1e-10 for x >= 20
    iWhich = 6
    raY(iI) = 0.0
  END IF
END DO

RETURN
END SUBROUTINE expintsuperfast3

!************************************************************************
! does expintegral raY = E3(raX), using 2 - 5 order polynom fits to expint3

SUBROUTINE expintsuperfast3matrix(iL,raaX,raY)


INTEGER, INTENT(IN OUT)                  :: iL
REAL, INTENT(OUT)                        :: raaX(kMaxPts,*)
REAL, INTENT(OUT)                        :: raY(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars


! output vars


! local vars
INTEGER :: iI,i,j,iWhich,iaMax(5),iCnt
REAL :: fastExpint,raaPout(5,6)

DATA (iaMax(i),i=1,5) /3, 6, 6, 3, 3/

!!this is for x > 0.00 and x <= 0.10
DATA ((raaPout(i,j),j=1,3),i=1,1)/  &
    4.999958864847588E-01, -9.704714999886551E-01,  1.357944294418206E+00/

!!this is for x > 0.10 and x <= 1.00
DATA ((raaPout(i,j),j=1,6),i=2,2)/  &
    4.964930888860016E-01,  -9.002236958326287E-01,  1.058882468914909E+00,  &
    -9.590227134465121E-01,   5.595049802866501E-01, -1.460130149043994E-01/

!!this is for x > 1.00 and x <= 6.00
DATA ((raaPout(i,j),j=1,6),i=3,3)/  &
    9.348983139232728E-03,   8.153540216831002E-01, -1.318262705075364E+00,  &
    1.507590471592263E+00,  -9.895617788618024E-01,  2.739113519528383E-01/

!!this is for x > 6.00 and x <= 10.00
DATA ((raaPout(i,j),j=1,3),i=4,4)/  &
    6.850600578733759E-03,   8.121694943009118E-01, -9.875485935560442E-01/

!!this is for x > 10.00 and x <= 20.00
DATA ((raaPout(i,j),j=1,3),i=5,5)/  &
    1.873706336320716E-03,   9.115323790575932E-01, -1.489123902774130E+00/

DO iI = 1,kMaxPts
  raY(iI) = 0.0
  IF (raaX(iI,iL) < 0.0) THEN
    WRITE(kStdErr,*) 'in expint iI,raaX(iI,iL) = ',iI,raaX(iI,iL)
    WRITE(kStdErr,*) 'cannot have negative values!'
    CALL DoStop
  ELSE IF (raaX(iI,iL) <= 0.1) THEN
    iWhich = 1
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raaX(iI,iL)**(j-1))
    END DO
  ELSE IF (raaX(iI,iL) <= 1.0) THEN
    iWhich = 2
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*(raaX(iI,iL)**(j-1))
    END DO
  ELSE IF (raaX(iI,iL) <= 6.0) THEN
    iWhich = 3
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raaX(iI,iL))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raaX(iI,iL))
  ELSE IF (raaX(iI,iL) <= 10.0) THEN
    iWhich = 4
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raaX(iI,iL))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raaX(iI,iL))
  ELSE IF (raaX(iI,iL) <= 20.0) THEN
    iWhich = 5
    iCnt = iaMax(iWhich)
    DO j = 1,iCnt
      raY(iI) = raY(iI) + raaPout(iWhich,j)*((1.0/raaX(iI,iL))**(j-1))
    END DO
    raY(iI) = raY(iI)*EXP(-raaX(iI,iL))
  ELSE !!   expint3(x) <= 1e-10 for x >= 20
    iWhich = 6
    raY(iI) = 0.0
  END IF
END DO

RETURN
END SUBROUTINE expintsuperfast3matrix

!************************************************************************