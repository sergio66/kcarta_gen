c Copyright 2014
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the misc angle  routines    ***************
c************************************************************************
c************************************************************************
c viewing geometry angles summary 
c define conv = pi/180, Re = 6400 km
c
c   [zang]=vaconv( sva, salt, alt );
c     ra = Re + alt     = radius at atmosphere height of interest ~ 20 km
c     rs = Re + salt    = radius at satellite orbit               ~ 705 km 
c     zang = arcsin ((rs/ra) * sin(conv*sva)) / conv              convert SCANAG (at satellite) to SATZEN (at surface)
c
c   [zang]=saconv( surfzang, alt );
c     ra = Re + alt     = radius at atmosphere height of interest ~ 20 km
c     rs = Re + salt    = radius at satellite orbit               ~ 705 km 
c     zang = arcsin ((re/ra) * sin(conv*surfzang)) / conv          convert SATZEN (at surface) to SCANANG (at satellite)
c
c
c        !! test going from surface satzen to sat scanang (at 705 km) back to surface satzen
c	!! should have satzen > scanang
c        print *,'orig satzen(h=0) = 45',saconv_sun(45.0,0.0,705.0),vaconv(saconv_sun(45.0,0.0,705.0),0.0,705.0)
c >>  orig satzen(h=0) = 45   39.54217       45.00000
c
c        !! test going from sat scanang (at 705 km) to surface satzen (0 km) back to sat scanang (at 705 km)
c	!! should have satzen > scanang
c        print *,'orig scanang(h=705) = 45',vaconv(45.0,0.0,705.0),saconv_sun(vaconv(45.0,0.0,705.0),0.0,705.0)	
c >>  orig scanang(h=705) = 45   51.75453       45.00000
  
c************************************************************************

c this function does the solar viewing angle conversion, by Sergio
c given solar zenith angle == angle wrt nadir at Earth surface,
c can compute solar angle at height h : from satzen(h=0) to scanang(h=H)

c ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
c SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

c directly compare to saconv.m
c directly compare to saconv.m
c directly compare to saconv.m

       REAL FUNCTION SACONV_SUN(LSZA, SURFALT, ALT )
                        
       IMPLICIT NONE

       ! input param
       REAL LSZA         !! solar/satellite zenith angle at local surface (which is not necessarily 0)
       REAL SURFALT      !! surface altitude ABOVE EARTH surface (km) (eg [6400 Km + ] aircraft at 05 km)
       REAL ALT          !! layer altitude   ABOVE EARTH surface (km) (eg [6400 Km + ] layer    at 15 km)

       REAL rX,rPi,rEarth

       rPi = 3.1415927
       rEarth = 6370.0
     
       !! p = Snell's law = n(i) r(i) sin(theta(i)) = constant => 
       !!     (R+h(i)) sin (theta(i)) = constant (assume n(i) = 1)
       !! so Re sin(theta_earth) = (Re+h(i))sin (theta(i))

       rX = (SURFALT + rEarth)/(ALT + rEarth) * sin(LSZA * rPi/180.0)
       !! SURFALT < ALT ==> ratio in front of sin(LSZA) < 1 
       !!   ==> local angles get "smaller" the higher you go 

       if (rX .GT. 0.99999) rX = 0.99999
       if (rX .LT. 0.00000) rX = 0.00000

       SACONV_SUN = ASIN(rX) * 180/rPi

       RETURN
       END

c************************************************************************

c ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
c SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

c directly compare to saconv.m
c directly compare to saconv.m
c directly compare to saconv.m

c this function does the surface --> arb height viewing angle conversion, by Scott
c copied from saconv.f in SARTA package, modified to 
c   have input lay height in km
c   have output angle in degrees
C Function to convert the surface solar/satellite zenith angle SZA into the
C local solar/satellite angle at altitude ALT (sunang/scanang)
c    from surface satzen(h=0) to scanang(h=H) at arb altitude H

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SZA     Surface Zenith Angle        degrees


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL fun  SACONV  local (height H) angle    degrees

C    Function to convert the Zenith Angle SZA at the Earth's
C    surface into a the local angle at altitude ALT.
C    The local angle generally varies slightly with altitude
C    due to the curvature of the Earth and its atmosphere.
C    The effect is largest at the maximum zenith angle, and
C    disappears as the zenith angle approaches 0 degrees.
C    Currently this function only considers the geometry of the
C    situation, and no refractive effects are included.
C
C    The layers of the atmosphere may be considered as concentric
C    rings with some average altitude. A ray traced thru these rings
C    at any viewing angle other than nadir will have a slightly
C    different angle (relative to the outward radial at the point
C    of intersection) in each ring. 
C
C    The local angle may be calculated (using The Law of
C    Sines) if we know:
C       The solar/satellire zenith angle, SZA, at the Earth's surface (ALT=0)
C       The layer altitude, ALT.
C       The radius of the Earth, RE.
C
C    The solution uses the law of sines and sin(180 - x) = sin(x)
       REAL FUNCTION ORIG_SACONV_SUN( SZA, ALT )
 
       implicit none 
       real sza,alt

c local variables
       real saconv,conv, re, ra

C      ------------------
C      Assign some values
C      ------------------
C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
C
C      RE = radius of the Earth (in km)
       RE=6.37E+03
C
C      RA = radius of the point to calc the angle at (in km)
C      Note: layer altitude already in kilometers
       RA = rE + ALT
C
C      -----------------
C      Do the conversion
C      -----------------
C
       SACONV=ASIN( (RE/RA) * SIN(CONV*SZA) )
       
c change back to degrees
       ORIG_SACONV_SUN = SACONV/conv

       RETURN
       END
       
c************************************************************************
C    Function to convert the AIRS satellite viewing angle into the
C    local path angle (      Viewing Angle CONVersion   )
c for downward looking instrument : from scanang(h=H) to satzen(h=0)

c directly compare to vaconv.m
c directly compare to vaconv.m
c directly compare to vaconv.m
c can Matlab debug version, but WARNING its input arguments are switched thusly
c   [zang]=vaconv( sva, salt, alt );
c cd /asl/matlab2012/science
c >> alt = [0:10:80]
c >> [zang]=vaconv(50, 705*1000, alt'*1000 )

c uses law of sines

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SVA     Satellite viewing angle     degrees
c    REAL      SALT    Satellite altitude          km

       REAL FUNCTION VACONV( SVA, ALT, SALT )

      IMPLICIT NONE

       REAL SVA, ALT, SALT

C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink

C      CONV = pi/180 = degrees to radians conversion factor
       CONV = 1.7453292E-02
       theta = sva*conv

C      RE = radius of the Earth (in km)
       RE=6.37E+03

C      RA = radius of the point to calc the angle at (in km)
       RA = rE + ALT

C      RS = radius of the satellite orbit (in km)
       RS = rE + SALT

C
C      -----------------
C      Do the conversion
C      -----------------
C
       rBoink =  (RS/RA) * SIN(CONV*SVA)
       IF (abs(rBoink) .GT. 1.0) rBoink = 0.9999
       
       RTEMP=ASIN(rBoink)

c change back to degrees
       VACONV = RTEMP/conv
        
       RETURN
       END

c************************************************************************
C    Function to convert the AIRS satellite viewing angle into the
C    local path angle (      Viewing Angle CONVersion   ) WITH RefIndex
c for downward looking instrument

c directly compare to vaconv.m, but this uses Snell
c directly compare to vaconv.m, but this uses Snell
c directly compare to vaconv.m, but this uses Snell
c can Matlab debug version, but WARNING its input arguments are switched thusly
c   [zang]=vaconv( sva, salt, alt );
c cd /asl/matlab2012/science
c >> alt = [0:10:80]
c >> [zang]=vaconv(50, 705*1000, alt'*1000 )

c uses law of sines


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SVA     Satellite viewing angle     degrees
c    REAL      SALT    Satellite altitude          km

       REAL FUNCTION VACONV_Snell( SVA, ALT, SALT, rNumberDensity )

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

       REAL SVA, ALT, SALT, rNumberDensity

C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink
       REAL lambda,n0,gamma,ref_ind,mu
       REAL rLosch             !Loschmidt number

       rLosch = (kAtm2mb*100/kBoltzmann/273 * 1e-6)    ! (p0/k T0)

C      CONV = pi/180 = degrees to radians conversion factor
       CONV = 1.7453292E-02
       theta = sva*conv

C      RE = radius of the Earth (in km)
       RE = 6367.512

c http://mintaka.sdsu.edu/GF/explain/atmos_refr/horizon.html
c http://www.ess.uci.edu/~cmclinden/link/xx/node45.html for TONS of detail
C      RA = radius of the point to calc the angle at (in km)
       RA = rE + ALT

C      RS = radius of the satellite orbit (in km)
       RS = rE + SALT

c ref ind of air at 10 um
       lambda = 10
       n0 = 64.328 + 29498.1/(146-1/lambda/lambda) + 255.4/(41 - 1/lambda/lambda)
       n0 = 1 + n0/1e6
       gamma = (n0*n0-1)/(n0*n0+2)

       IF (ALT .GE. 200) THEN
         ref_ind = 1.0
       ELSE
         mu = rNumberDensity/rLosch * gamma
         ref_ind = sqrt((2*mu+1)/(1-mu))
       END IF

C
C      -----------------
C      Do the conversion
C      -----------------
C
       rBoink =  abs((RS/RA) * SIN(CONV*SVA)/ref_ind)
       if (abs(rBoink) .gt. 1.0) rBoink = 0.9999

       RTEMP = ASIN(rBoink)
       
c change back to degrees
       VACONV_SNELL = RTEMP/conv
        
       RETURN
       END

c************************************************************************
c this subroutine finds the layer dependent satellite viewing angle
      SUBROUTINE FindSunLayerAngles(rSatHeight,rSurfHeight,iAtm,iaNumLayer,iaaRadLayer,
     $           raLayHgt,rPrBdry1,rPrBdry2,rSunAngle,raSunAngles) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatHeight = tells us if we need to do ray tracing (> 0) or not (< 0)
c raLayHgt   = height of individual layers, read in from profile, in meters
c rPrBdry1   = start pressure boundary
c rPrBdry2   = stop pressure boundary
c rSunAngle  = sun angle at TOA
c raSunAngles = layer dependent sun viewing angle
      REAL rPrBdry1,rPrBdry2,rSunAngle,rSatHeight,rSurfHeight
      REAL raLayHgt(kProfLayer),raSunAngles(kProfLayer)
      INTEGER iAtm,iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

      REAL salt,saconv_sun,rSunHeight
      INTEGER iI,iX,iMin,iMax,iaRadLayer(kProfLayer)

      rSunHeight = 150.0e9  !! in km

c as default, set all angles to be the sun view angle
      DO iI=1,kProfLayer
        raSunAngles(iI) = rSunAngle
      END DO

      iMin = +1000000
      iMax = -1000000
      DO iI=1,iaNumlayer(iAtm)
        iaRadLayer(iI) = iaaRadLayer(iAtm,iI)
        IF (iaRadLayer(iI) .GT. iMax) iMax = iaRadLayer(iI)
        IF (iaRadLayer(iI) .LT. iMin) iMin = iaRadLayer(iI)
      END DO

c *********************** rSatHeight, raLayHgt is in METERS ***************

      IF ((rSatHeight .GT. 0.0) .AND. (abs(rSunAngle) .GT. 1.0e-4)) THEN   
        !have to find layer dependent angles
        IF (rPrBdry1 .GT. rPrBdry2) THEN !downward looking instr
          DO iI=1,kProfLayer
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raSunAngles(iI) = saconv_sun(abs(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
              IF (rSunAngle .lt. 0.0) raSunAngles(iI) = -raSunAngles(iI)
              IF (kOuterLoop .EQ. 1) THEN
                IF (iI .GE. iMin .AND. iX .LE. iMax) THEN
                  iX = (iI - iMin + 1)
                ELSE
                  iX = -1
                END IF
                IF (iX .EQ. 1) write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
                write(kStdWarn,*)'sun: lay#/rad# lay/surfhgt, localzen/traced angle ',
     $            iI,iX,raLayHgt(iI)/1000,rSurfHeight/1000,rSunAngle,raSunAngles(iI)
              END IF
            END IF
          END DO
        END IF
        IF (rPrBdry2 .GT. rPrBdry1) THEN !upward looking instr
          salt=705.00
          DO iI=1,kProfLayer
c            write(kStdWarn,*)'sun up lay hgt ',iI,raLayHgt(iI)/1000
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raSunAngles(iI) = saconv_sun(abs(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
              IF (rSunAngle .lt. 0.0) raSunAngles(iI) = -raSunAngles(iI)
              IF (((abs(kLongOrShort) .EQ. 2) .AND. (kOuterLoop .EQ. 1)) .OR. 
     $            (abs(kLongOrShort) .LE. 1)) THEN
                write(kStdWarn,*)'sun up lay hgt, orig/traced angle ',
     $            iI,raLayHgt(iI)/1000,rSunAngle,raSunAngles(iI)
              END IF
            END IF
          END DO
        END IF
      END IF

c      call dostop

      RETURN
      END

c************************************************************************ 
c this subroutine finds the layer dependent satellite viewing angle
      SUBROUTINE FindLayerAngles(rSatHeight,raLayHgt, 
     $              rPrBdry1,rPrBdry2,rSatAngle,raLayAngles,
     $              iAtm,iaaRadLayer,iaNumlayer,raNumberDensity)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatHeight = satellite height in meters 
c   (from spec : http://asl.umbc.edu/pub/motteler/hdf/rtpV105/doc/rtpspec.html)
c [sergio@tara-fe1 WORK]$ grep -in zobs   /asl/packages/sartaV108/Src/*.f
c /asl/packages/sartaV108/Src/rdrtp_df.f:329:       ZSAT=PROF.zobs/1000
c  /asl/packages/sartaV108/Src/rdrtp.f:328:       ZSAT=PROF.zobs/1000
c
c raLayHgt = height of individual layers, read in from profile, in meters
c
c rPrBdry1 = start pressure boundary
c rPrBdry2= stop pressure boundary
c rSatAngle = satellite viewing angle (scanang)
c raLayAngles = layer dependent satellite viewing angle (satzen)
      REAL rSatHeight,rPrBdry1,rPrBdry2,rSatAngle
      REAL raLayHgt(kProfLayer),raLayAngles(kProfLayer)
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iAtm
      REAL raNumberDensity(kProfLayer)

      REAL vaconv,vaconv_Snell,saconv_sun
      REAL raLayAnglesSnell(kProfLayer),raLayAnglesNoSnell(kProfLayer)
      INTEGER iI,iaRadLayer(kProfLayer),iMin,iMax,iX

c as default, set all angles to be the satellite view angle
      DO iI=1,kProfLayer
        raLayAngles(iI)      = rSatAngle
        raLayAnglesSnell(iI) = rSatAngle
      END DO

      iMin = +1000000
      iMax = -1000000
      DO iI = 1,iaNumlayer(iAtm)
        iaRadLayer(iI) = iaaRadLayer(iAtm,iI)
        IF (iaRadLayer(iI) .GT. iMax) iMax = iaRadLayer(iI)
        IF (iaRadLayer(iI) .LT. iMin) iMin = iaRadLayer(iI)
      END DO

c *********** vaconv assumes layheight,satheight is in KM ***************
c *********** SatHeight is in METERS (from RTP), and raLayHgt is in METERS ****

c orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)
c orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)
c orig/traced angle  <<< --------- >>> scanang(at satellite)/satzen(local angle)

      IF ((rSatHeight .GT. 0.0) .AND. (abs(rSatAngle) .GT. 1.0e-4)) THEN   
        !have to find layer dependent angles
        IF (rPrBdry1 .GT. rPrBdry2) THEN !downward looking instr
          DO iI=1,kProfLayer
c            print *,iI,rSatHeight,raLayHgt(iI)
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raLayAnglesNoSnell(iI) = vaconv(abs(rSatAngle),raLayHgt(iI)/1000,
     $                                rSatHeight/1000)
              raLayAnglesSnell(iI) = vaconv_Snell(abs(rSatAngle),raLayHgt(iI)/1000,
     $                                rSatHeight/1000,raNumberDensity(iI))
              raLayAngles(iI) = raLayAnglesSnell(iI)
              IF (rSatAngle .lt. 0.0) raLayAngles(iI) = -raLayAngles(iI)
              IF (kOuterLoop .EQ. 1) THEN
                IF (iI .GE. iMin .AND. iX .LE. iMax) THEN
                  iX = (iI - iMin + 1)
                ELSE
                  iX = -1
                END IF
              IF (iX .EQ. 1) write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
              write(kStdWarn,*)'dn : lay#/rad# lay/sat hgt, satellite scanang/local satzen angle ',
     $              iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
              END IF
            END IF
          END DO
      ELSE  
        !no need to do anything much, the angles are so close to nadir
        DO iI=1,kProfLayer
          write(kStdWarn,*)'dn : lay#/rad# lay/sat hgt, satellite scanang/local satzen angle ',
     $              iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
        END DO

      END IF

        IF (rPrBdry2 .GT. rPrBdry1) THEN !upward looking instr
c          print *,'FindLayerAngles',iMin,iMax,iAtm,rSatHeight,iaNumLayer(iAtm),
c     $               raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000
          DO iI=1,kProfLayer
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raLayAnglesNoSnell(iI) = saconv_sun(abs(rSatAngle),
     $                                    raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000,
     $                                    raLayHgt(iI)/1000)
              raLayAngles(iI) = raLayAnglesSnell(iI)
              raLayAngles(iI) = raLayAnglesNoSnell(iI)
              IF (rSatAngle .lt. 0.0) raLayAngles(iI) = -raLayAngles(iI)
              IF (kOuterLoop .EQ. 1) THEN
                IF (iI .GE. iMin .AND. iX .LE. iMax) THEN
                  iX = (iI - iMin + 1)
                ELSE
                   iX = -1
                 END IF
                IF (iX .EQ. 1) write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
                write(kStdWarn,*)'up : lay#/rad# lay/sat hgt, local satzen/layer iI angle ',
     $            iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
              END IF
            END IF
          END DO
        END IF
      END IF

      RETURN
      END

c************************************************************************ 
c this subroutine computes the gauss-legendre abscissa weights and points
c from Numerical Recipes
c nn is the number of points <= kProfLayer, x,w are output abscissa and wts
      SUBROUTINE FindGauss(nn,daX,daW)
       
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER nn
      DOUBLE PRECISION daX(kGauss),daW(kGauss)

      DOUBLE PRECISION daX1(2*kGauss),daW1(2*kGauss)
      DOUBLE PRECISION x1,x2
      INTEGER m,j,i,n
      DOUBLE PRECISION z1,z,xm,xl,pp,p1,p2,p3,epss

      epss = 3.0e-11
 
      x1 = -1.0D0
      x2 = +1.0D0
      
      IF ((nn .gt. kGauss) .or. (nn .lt. 0)) THEN
        write (kStdErr,*) 'need 0 < nn <= kGauss' 
        CALL DoStop
      END IF

      n=nn*2

      IF (MOD(n,2) .EQ. 1) THEN
        write (kStdErr,*) 'need n to be even' 
        CALL DoSTOP
      END IF

      IF (x2 .LT. x1) THEN
        xm = x1
        x1 = x2
        x2 = xm
      END IF

      m  = (n+1)/2

      m  = n/2
      xm = 0.5*(x2+x1)
      xl = 0.5*(x2-x1)

      DO i = 1,m                    !loop over desired roots
        z = cos(kPi*(i-0.25)/(n+0.5))
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
        IF (abs(z-z1) .gt. epss) THEN
          GOTO 20
        END IF
           
        daX(i)     = xm+xl*z
        daW(i)     = 2*xl/((1-z*z)*pp*pp)

        daX1(i)     = xm-xl*z
        daX1(n+1-i) = xm+xl*z
        daW1(i)     = 2*xl/((1-z*z)*pp*pp)
        daW1(n+1-i) = daW(i)
      END DO
 
      RETURN
      END

c************************************************************************ 
c this subroutine reads in the first moment (l=1) quadrature points and weights
c using results in Table 1, second column (m=1) of
c "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
c Journal of Atmospheric Sciences, v57, pg 753 (2000)
      SUBROUTINE FindGauss2old(nn,daX,daW)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER nn
      DOUBLE PRECISION daX(kGauss),daW(kGauss)

      IF (nn .EQ. 1) THEN
        daX(1) = 2.0/3.0
        daW(1) = 1.0/2.0
      ELSEIF (nn .EQ. 2) THEN
        daX(1) = 0.3550510
        daX(2) = 0.8449489
        daW(1) = 0.1819856
        daW(2) = 0.3180414
      ELSEIF (nn .EQ. 3) THEN
        daX(1) = 0.2123405
        daX(2) = 0.5905331
        daX(3) = 0.9114120
        daW(1) = 0.0698270
        daW(2) = 0.2292411
        daW(3) = 0.2009319
      ELSEIF (nn .EQ. 4) THEN
        daX(1) = 0.1397599
        daX(2) = 0.4164096
        daX(3) = 0.7231570
        daX(4) = 0.9428958
        daW(1) = 0.0311810
        daW(2) = 0.1298475
        daW(3) = 0.2034646
        daW(4) = 0.1355069
      ELSE
        write(kStdErr,*) 'FindGauss2 : need nn = 1,2,3 or 4 not ',nn
        CALL DoStop
      END IF

      RETURN
      END

c************************************************************************ 
c this subroutine reads in the first moment (l=1) quadrature points and weights
c using results in Table 1, second column (m=1) of
c "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
c Journal of Atmospheric Sciences, v57, pg 753 (2000)
c also see RRTM code v3.3 rtreg.ffor more decimal points
      SUBROUTINE FindGauss2(nn,daX,daW)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER nn,ii
      DOUBLE PRECISION daX(kGauss),daW(kGauss)

      IF (nn .EQ. 1) THEN
        daX(1) = 3.0/2.0
        daW(1) = 1.0/2.0
      ELSEIF (nn .EQ. 2) THEN
        daX(1) = 2.81649655
        daX(2) = 1.18350343
        daW(1) = 0.1819586183
        daW(2) = 0.3180413817
      ELSEIF (nn .EQ. 3) THEN
        daX(1) = 4.70941630
        daX(2) = 1.69338507
        daX(3) = 1.09719858
        daW(1) = 0.0698269799
        daW(2) = 0.2292411064
        daW(3) = 0.2009319137
      ELSEIF (nn .EQ. 4) THEN
        daX(1) = 7.15513024
        daX(2) = 2.40148179
        daX(3) = 1.38282560
        daX(4) = 1.06056257
        daW(1) = 0.0311809710
        daW(2) = 0.1298475476
        daW(3) = 0.2034645680
        daW(4) = 0.1355069134
      ELSE
        write(kStdErr,*) 'FindGauss2 : need nn = 1,2,3 or 4 not ',nn
        CALL DoStop
      END IF

      DO ii = 1,nn
        !change from secant to cosine
        daX(ii) = 1.0/daX(ii)
      END DO
      
      RETURN
      END

c************************************************************************ 
