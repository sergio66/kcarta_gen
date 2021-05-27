! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_angles

USE basic_common
USE n_misc

IMPLICIT NONE

CONTAINS

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

! this subroutine finds the layer dependent satellite viewing angle
    SUBROUTINE FindSunLayerAngles(rSatHeight,rSurfHeight,iAtm,iaNumLayer,iaaRadLayer, &
    raLayHgt,rPrBdry1,rPrBdry2,rSunAngle,raSunAngles)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rSatHeight = tells us if we need to do ray tracing (> 0) or not (< 0)
! raLayHgt   = height of individual layers, read in from profile, in meters
! rPrBdry1   = start pressure boundary
! rPrBdry2   = stop pressure boundary
! rSunAngle  = sun angle at TOA
! raSunAngles = layer dependent sun viewing angle
    REAL :: rPrBdry1,rPrBdry2,rSunAngle,rSatHeight,rSurfHeight
    REAL :: raLayHgt(kProfLayer),raSunAngles(kProfLayer)
    INTEGER :: iAtm,iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

    REAL :: salt,rSunHeight
    INTEGER :: iI,iX,iMin,iMax,iaRadLayer(kProfLayer)

    rSunHeight = 150.0e9  !! in km

! as default, set all angles to be the sun view angle
    raSunAngles = rSunAngle

    iMin = +1000000
    iMax = -1000000
    DO iI=1,iaNumlayer(iAtm)
      iaRadLayer(iI) = iaaRadLayer(iAtm,iI)
      IF (iaRadLayer(iI) > iMax) iMax = iaRadLayer(iI)
      IF (iaRadLayer(iI) < iMin) iMin = iaRadLayer(iI)
    END DO

! *********************** rSatHeight, raLayHgt is in METERS ***************

    IF ((rSatHeight > 0.0) .AND. (abs(rSunAngle) > 1.0e-4)) THEN
      !have to find layer dependent angles
      IF (rPrBdry1 > rPrBdry2) THEN !downward looking instr
        DO iI=1,kProfLayer
          IF (rSatHeight > raLayHgt(iI)) THEN
            raSunAngles(iI) = saconv_sun(abs(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
            IF (rSunAngle < 0.0) raSunAngles(iI) = -raSunAngles(iI)
            IF (kOuterLoop == 1) THEN
              IF (iI >= iMin .AND. iX <= iMax) THEN
                iX = (iI - iMin + 1)
              ELSE
                iX = -1
              END IF
              IF (iX == 1) write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
              write(kStdWarn,*)'sun: lay#/rad# lay/surfhgt, localzen/traced angle ', &
                        iI,iX,raLayHgt(iI)/1000,rSurfHeight/1000,rSunAngle,raSunAngles(iI)
            END IF
          END IF
        END DO
      END IF
      IF (rPrBdry2 > rPrBdry1) THEN !upward looking instr
        salt=705.00
        DO iI=1,kProfLayer
          ! write(kStdWarn,*)'sun up lay hgt ',iI,raLayHgt(iI)/1000
          IF (rSatHeight > raLayHgt(iI)) THEN
            raSunAngles(iI) = saconv_sun(abs(rSunAngle),rSurfHeight/1000,raLayHgt(iI)/1000)
            IF (rSunAngle < 0.0) raSunAngles(iI) = -raSunAngles(iI)
            IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. (abs(kLongOrShort) <= 1)) THEN
              write(kStdWarn,*)'sun up lay hgt, orig/traced angle ',iI,raLayHgt(iI)/1000,rSunAngle,raSunAngles(iI)
            END IF
          END IF
        END DO
      END IF
    END IF

    RETURN
    end SUBROUTINE FindSunLayerAngles

!************************************************************************
! this subroutine finds the layer dependent satellite viewing angle
    SUBROUTINE FindLayerAngles(rSatHeight,raLayHgt, &
    rPrBdry1,rPrBdry2,rSatAngle,raLayAngles, &
    iAtm,iaaRadLayer,iaNumlayer,raNumberDensity)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: rSatHeight,rPrBdry1,rPrBdry2,rSatAngle
    REAL :: raLayHgt(kProfLayer),raLayAngles(kProfLayer)
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iAtm
    REAL :: raNumberDensity(kProfLayer)

    REAL :: raLayAnglesSnell(kProfLayer),raLayAnglesNoSnell(kProfLayer)
    INTEGER :: iI,iaRadLayer(kProfLayer),iMin,iMax,iX,iDefault,iUseSnell,iStartLayer

! as default, set all angles to be the satellite view angle

    write(kStdWarn,*) 'FindLayerAngles : iAtm,rSatAngle,rSatHeight = ',iAtm,rSatAngle,rSatHeight

    raLayAngles      = rSatAngle
    raLayAnglesSnell = rSatAngle

    iDefault = -1   !! no  Snell, yes layer curvature, similar to SARTA
          
    iUseSnell = +1  !! yes Snell, yes layer curvature
    iUseSnell = -1  !! no  Snell, yes layer curvature, similar to SARTA
    iUseSnell = 0   !! no  Snell, no  layer curvature
          
    iUseSnell = iaaOverrideDefault(2,7)
    IF (abs(iUseSnell) > 1) THEN
      write(kStdErr,*) 'invalid iUseSnell ',iUseSnell
      CALL DoStop
    END IF
    if ((iDefault /= iUseSnell) .AND. (kOuterLoop == 1)) THEN
      write(kStdWarn,*) 'not(-1)/using(+1) Snell law in FindLayerAngles (raytrace thru layers)',iUseSnell
      write(kStdErr,*)  'not(-1)/using(+1) Snell law in FindLayerAngles (raytrace thru layers)',iUseSnell
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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DOWNLOOK INSTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF ((rSatHeight > 0.0) .AND. (abs(rSatAngle) > 1.0e-4) .AND. (abs(iUseSnell) == 0)) THEN
      iStartLayer = iaRadLayer(1)
      !! no  Snell, no  layer curvature
      DO iI=1,kProfLayer
        raLayAnglesNoSnell(iI) = abs(rSatAngle)
        raLayAnglesSnell(iI) = abs(rSatAngle)
        IF (kOuterLoop == 1) THEN
          IF (iI >= iMin .AND. iX <= iMax) THEN
            iX = (iI - iMin + 1)
          ELSE
            iX = -1
          END IF
          IF (iX == 1) THEN
            write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
            write(kStdWarn,*) 'Downlook Instr'
            write(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
            write(kStdWarn,*) 'TOA scanang, GND satzen = ',abs(rSatAngle),vaconv(abs(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
            write(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
          END IF
          write(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI), &
                &  1.0/cos(raLayAngles(iI)*kPi/180.0)
        END IF      !! if kOuterLoop == 1
      END DO

    ELSEIF ((rSatHeight > 0.0) .AND. (abs(rSatAngle) > 1.0e-4) .AND. (abs(iUseSnell) == 1)) THEN

      !have to find layer dependent angles as there is layer curvature
      IF (rPrBdry1 > rPrBdry2) THEN !downward looking instr
        iStartLayer = iaRadLayer(1)
        DO iI=1,kProfLayer
          ! print *,iI,rSatHeight,raLayHgt(iI)
          IF (rSatHeight > raLayHgt(iI)) THEN
            raLayAnglesNoSnell(iI) = vaconv(abs(rSatAngle),raLayHgt(iI)/1000-raLayHgt(iStartLayer)/1000,rSatHeight/1000)
            raLayAnglesSnell(iI) = vaconv_Snell(abs(rSatAngle),raLayHgt(iI)/1000-raLayHgt(iStartLayer)/1000, &
                                                rSatHeight/1000,raNumberDensity(iI))

            
            IF (iI == 1) THEN
              raLayAnglesNoSnell(iI) = vaconv(abs(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
              raLayAnglesSnell(iI) = vaconv_Snell(abs(rSatAngle),kSurfAlt/1000, &
                                                  rSatHeight/1000,raNumberDensity(iI))
            ELSE
              raLayAnglesNoSnell(iI) = vaconv(abs(rSatAngle),raLayHgt(iI-1)/1000,rSatHeight/1000)
              raLayAnglesSnell(iI) = vaconv_Snell(abs(rSatAngle),raLayHgt(iI-1)/1000, &
                                                  rSatHeight/1000,raNumberDensity(iI))
            END IF
            ! diff between Snell and noSnell is less than 2e-2
            !              print *,iI,raLayAnglesSnell(iI),raLayAnglesNoSnell(iI),raLayAnglesNoSnell(iI)-raLayAnglesSnell(iI)
            ! but Scott/SARTA uses NoSnell
            IF (iUseSnell == 1) THEN
              raLayAngles(iI) = raLayAnglesSnell(iI)    !!! sergio, use layer curvature+Snell
            ELSEIF (iUseSnell == -1) THEN
              raLayAngles(iI) = raLayAnglesNoSnell(iI)  !!! scott/SARTA, only use layer curvature
            END IF
            IF (rSatAngle < 0.0) raLayAngles(iI) = -raLayAngles(iI)
            IF (kOuterLoop == 1) THEN
              IF (iI >= iMin .AND. iX <= iMax) THEN
                iX = (iI - iMin + 1)
              ELSE
                iX = -1
              END IF
              IF (iX == 1) THEN
                write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
                write(kStdWarn,*) 'Downlook Instr'
                write(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
                write(kStdWarn,*) 'TOA scanang, GND satzen = ',abs(rSatAngle),vaconv(abs(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
                write(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
              END IF
              write(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI), &
                        &  1.0/cos(raLayAngles(iI)*kPi/180.0)
            END IF      !! if kOuterLoop == 1
          END IF
        END DO
      ELSE
        ! no need to do anything much, the angles are so close to nadir
        iX = -1
        write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
        write(kStdWarn,*) 'Downlook Instr'
        write(kStdWarn,*) 'Surf Pressure (mb), Surf Altitude (m) = ',kSurfPress,kSurfAlt
        write(kStdWarn,*) 'TOA scanang, GND satzen = ',abs(rSatAngle),vaconv(abs(rSatAngle),kSurfAlt/1000,rSatHeight/1000)
        write(kStdWarn,*)'lay#/rad#   lay hgt   sat hgt  sat.scanang loc.satzen sec(satzen)'
        DO iI=1,kProfLayer
          write(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI),1.0
        END DO
      END IF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UPLOOK INSTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      IF ((rPrBdry2 > rPrBdry1) .AND. (abs(iUseSnell) == 0)) THEN !upward looking instr
        iStartLayer = iaRadLayer(kProfLayer)
        DO iI=1,kProfLayer
          raLayAnglesNoSnell(iI) = abs(rSatAngle)
          raLayAngles(iI) = raLayAnglesSnell(iI)
          IF (kOuterLoop == 1) THEN
            IF (iI >= iMin .AND. iX <= iMax) THEN
              iX = (iI - iMin + 1)
            ELSE
              iX = -1
            END IF
            IF (iX == 1) THEN
              write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
              write(kStdWarn,*)'up : lay#/rad# lay/sat hgt, local satzen/layer iI angle'
            END IF
            write(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
          END IF
        END DO
                    
      ELSEIF ((rPrBdry2 > rPrBdry1) .AND. (abs(iUseSnell) == 1)) THEN !upward looking instr
        !          print *,'FindLayerAngles',iMin,iMax,iAtm,rSatHeight,iaNumLayer(iAtm),
        !     $               raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000
        iStartLayer = iaRadLayer(kProfLayer)
        DO iI=1,kProfLayer
          IF (rSatHeight > raLayHgt(iI)) THEN
            raLayAnglesNoSnell(iI) = saconv_sun(abs(rSatAngle), &
                    raLayHgt(iaaRadLayer(iAtm,iaNumLayer(iAtm)))/1000, &
                    raLayHgt(iI)/1000)
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
                write(kStdWarn,*) '------------>>> these are used by Atmosphere ',iAtm
                write(kStdWarn,*)'up : lay#/rad# lay/sat hgt, local satzen/layer iI angle'
              END IF
              write(kStdWarn,999) iI,iX,raLayHgt(iI)/1000,rSatHeight/1000,rSatAngle,raLayAngles(iI)
            END IF
          END IF
        END DO
      END IF
    END IF
    
    999 FORMAT(I3,' ',I3,'  ',5(F10.4,' '))
     
    RETURN
    end SUBROUTINE FindLayerAngles

!************************************************************************
! this subroutine computes the gauss-legendre abscissa weights and points
! from Numerical Recipes
! nn is the number of points <= kProfLayer, x,w are output abscissa and wts
    SUBROUTINE FindGauss(nn,daX,daW)
           
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: nn
    DOUBLE PRECISION :: daX(kGauss),daW(kGauss)

    DOUBLE PRECISION :: daX1(2*kGauss),daW1(2*kGauss)
    DOUBLE PRECISION :: x1,x2
    INTEGER :: m,j,i,n
    DOUBLE PRECISION :: z1,z,xm,xl,pp,p1,p2,p3,epss

    epss = 3.0e-11
     
    x1 = -1.0D0
    x2 = +1.0D0
          
    IF ((nn > kGauss) .OR. (nn < 0)) THEN
      write (kStdErr,*) 'need 0 < nn <= kGauss'
      CALL DoStop
    END IF

    n=nn*2

    IF (MOD(n,2) == 1) THEN
      write (kStdErr,*) 'need n to be even'
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
      z = cos(kPi*(i-0.25)/(n+0.5))
 20   CONTINUE
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
      IF (abs(z-z1) > epss) THEN
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
    end SUBROUTINE FindGauss

!************************************************************************
! this subroutine reads in the first moment (l=1) quadrature points and weights
! using results in Table 1, second column (m=1) of
! "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
! Journal of Atmospheric Sciences, v57, pg 753 (2000)
    SUBROUTINE FindGauss2old(nn,daX,daW)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: nn
    DOUBLE PRECISION :: daX(kGauss),daW(kGauss)

    IF (nn == 1) THEN
      daX(1) = 2.0/3.0
      daW(1) = 1.0/2.0
    ELSEIF (nn == 2) THEN
      daX(1) = 0.3550510
      daX(2) = 0.8449489
      daW(1) = 0.1819856
      daW(2) = 0.3180414
    ELSEIF (nn == 3) THEN
      daX(1) = 0.2123405
      daX(2) = 0.5905331
      daX(3) = 0.9114120
      daW(1) = 0.0698270
      daW(2) = 0.2292411
      daW(3) = 0.2009319
    ELSEIF (nn == 4) THEN
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
    end SUBROUTINE FindGauss2old

!************************************************************************
! this subroutine reads in the first moment (l=1) quadrature points and weights
! using results in Table 1, second column (m=1) of
! "Gaussian Quadrature and Application to Infrared Radiation" by J.Li
! Journal of Atmospheric Sciences, v57, pg 753 (2000)
! also see RRTM code v3.3 rtreg.ffor more decimal points
    SUBROUTINE FindGauss2(nn,daX,daW)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: nn,ii
    DOUBLE PRECISION :: daX(kGauss),daW(kGauss)

    IF (nn == 1) THEN
      daX(1) = 3.0/2.0
      daW(1) = 1.0/2.0
    ELSEIF (nn == 2) THEN
      daX(1) = 2.81649655
      daX(2) = 1.18350343
      daW(1) = 0.1819586183
      daW(2) = 0.3180413817
    ELSEIF (nn == 3) THEN
      daX(1) = 4.70941630
      daX(2) = 1.69338507
      daX(3) = 1.09719858
      daW(1) = 0.0698269799
      daW(2) = 0.2292411064
      daW(3) = 0.2009319137
    ELSEIF (nn == 4) THEN
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
      !c hange from secant to cosine
      daX(ii) = 1.0/daX(ii)
    END DO
          
    RETURN
    end SUBROUTINE FindGauss2

!************************************************************************
! this function does the expintegral raY = E3(raX), using recursion
    SUBROUTINE expint3(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI,iaPos(kMaxPts),npos
    REAL :: eps

    eps = 1.0e-28

    CALL expint(raX,raY)               !do y = expint(x)

    CALL FindInVector(kMaxPts,raX < eps,npos,iaPos)
    IF (npos > 0) THEN
      raY(iaPos(1:npos)) = 0.5
    END IF

    CALL FindInVector(kMaxPts,raX >= eps,npos,iaPos)
    IF (npos > 0) THEN
      raY(iaPos(1:npos)) = (exp(-raX(iaPos(1:npos)))*(1-raX(iaPos(1:npos))) + (raX(iaPos(1:npos))**2)*raY(iaPos(1:npos)))/2.0
    END IF


    RETURN
    end SUBROUTINE expint3

!************************************************************************
! this function does the expintegral raY = E2(raX), using recursion
    SUBROUTINE expint2(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI,iaPos(kMaxPts),npos
    REAL :: eps

    eps = 1.0e-28

    CALL expint(raX,raY)               !do y = expint(x)

    CALL FindInVector(kMaxPts,raX < eps,npos,iaPos)
    IF (npos > 0) THEN
      raY(iaPos(1:npos)) = 1.0
    END IF

    CALL FindInVector(kMaxPts,raX >= eps,npos,iaPos)
    IF (npos > 0) THEN
      raY(iaPos(1:npos)) = exp(-raX(iaPos(1:npos))) - raX(iaPos(1:npos))*raY(iaPos(1:npos))
    END IF

    RETURN
    end SUBROUTINE expint2

!************************************************************************
! this function does the expintegral raY = E1(raX), ala Matlab
! really, this can be found in any Math handbook, or on the web
! assumes input raX has no elements == 0, and all are real
    SUBROUTINE expint(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI,iJ,iK,iFr,iFound,iN
    INTEGER :: iaPositive(kMaxPts),iPos,iaNegative(kMaxPts),iNeg
    DOUBLE PRECISION :: p(9),egamma,daX(kMaxPts),daY(kMaxPts),dImag
    COMPLEX*16 daZ(kMaxPts)
    DOUBLE PRECISION :: daPterm(kMaxPts),daTerm(kMaxPts)
    DOUBLE PRECISION :: am2(kMaxPts),am1(kMaxPts),bm2(kMaxPts),bm1(kMaxPts)
    DOUBLE PRECISION :: f(kMaxPts),oldf(kMaxPts),alpha

    REAL, DIMENSION(:), ALLOCATABLE :: raBeta,raA,raB
    INTEGER :: AllocateStatus,DeAllocateStatus
    INTEGER :: iaPos(kMaxPts),npos

    DATA (P(iJ),iJ=1,9)/ &
    -3.602693626336023d-09, -4.819538452140960d-07, -2.569498322115933d-05, &
    -6.973790859534190d-04, -1.019573529845792d-02, -7.811863559248197d-02, &
    -3.012432892762715d-01, -7.773807325735529d-01,  8.267661952366478d+00/

    daY = 0.0d0
    DO iJ = 1,9
      iI = 9-iJ
      daY = daY + p(iJ)*((raX*1.0D0)**iI)
    END DO

! polyv = polyval(p,real(x))
! k = find( abs(imag(x)) <= polyv )   this is iPos, except imag(x) == 0 always
    iPos = 0
    iNeg = 0
    dImag = 0.0D0
    DO iFr = 1,kMaxPts
      IF (dImag <= daY(iFr)) THEN
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
      daX(iaPositive(1:iPos)) = raX(iaPositive(1:iPos))*1.0d0
      daZ(iaPositive(1:iPos)) = -egamma - cdlog(dcmplx(daX(iaPositive(1:iPos)),0.0D0))
      daPterm(iaPositive(1:iPos)) = daX(iaPositive(1:iPos))
      daterm(iaPositive(1:iPos))  = daX(iaPositive(1:iPos))

      iJ = 1

 10   CONTINUE
      iJ = iJ + 1

      daZ(iaPositive(1:iPos)) = daZ(iaPositive(1:iPos)) + daTerm(iaPositive(1:iPos))
      daPterm(iaPositive(1:iPos)) = -daX(iaPositive(1:iPos)) * daPterm(iaPositive(1:iPos))/iJ
      daTerm(iaPositive(1:iPos))  = daPterm(iaPositive(1:iPos))/iJ

      ! check for convergence : any of the terms fail then go back to 10
      iFound = -1
      CALL FindInVector(iPos,abs(daTerm(iaPositive(1:iPos))) > 1.0d-16*cdabs(daZ(iaPositive(1:iPos))),npos,iaPos)
      IF (npos > 0) THEN
        iFound = +1
        GOTO 11
      END IF
!      DO iK = 1,iPos
!        iFr = iaPositive(iK)
!        IF (abs(daTerm(iFr)) > 1.0d-16*cdabs(daZ(iFr))) THEN
!          iFound = +1
!          GOTO 11
!        END IF
!      END DO
 11   CONTINUE

      IF (iFound > 0) THEN
        GOTO 10
      END IF

      raY(iaPositive(1:iPos)) = sngl(dreal(daZ(iaPositive(1:iPos))))
    END IF

! ---------- these are for x < 0  ------------------------------------------
    IF (iNeg >= 1) THEN

      ALLOCATE ( raBeta(iNeg), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory for raBeta ***"      
      ALLOCATE ( raA(iNeg), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory for raA ***"      
      ALLOCATE ( raB(iNeg), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory for raB ***"      

      iN = 1  !calculating E1(x)

      daX(iaNegative(1:iNeg))  = raX(iaNegative(1:iNeg))*1.0d0
      am2(iaNegative(1:iNeg))  = 0.0D0
      bm2(iaNegative(1:iNeg))  = 1.0D0
      am1(iaNegative(1:iNeg))  = 1.0D0
      bm1(iaNegative(1:iNeg))  = daX(iaNegative(1:iNeg))
      f(iaNegative(1:iNeg))    = am1(iaNegative(1:iNeg))/bm1(iaNegative(1:iNeg))
      oldf(iaNegative(1:iNeg)) = 1.0D100

      iJ = 2
   20 CONTINUE

      ! calculate the coefficients of the recursion formulas for j even
      alpha = iN - 1+(iJ/2) ! note: beta = 1

      !calculate A(j), B(j), and f(j)
      raA = am1(iaNegative(1:iNeg)) + alpha * am2(iaNegative(1:iNeg))
      raB = bm1(iaNegative(1:iNeg)) + alpha * bm2(iaNegative(1:iNeg))
         
      ! save new normalized variables for next pass through the loop
      !  note: normalization to avoid overflow or underflow
      am2(iaNegative(1:iNeg)) = am1(iaNegative(1:iNeg)) / raB
      bm2(iaNegative(1:iNeg)) = bm1(iaNegative(1:iNeg)) / raB
      am1(iaNegative(1:iNeg)) = raA / raB
      bm1(iaNegative(1:iNeg)) = 1.0D0
          
      oldf(iaNegative(1:iNeg)) =f(iaNegative(1:iNeg))
      f(iaNegative(1:iNeg)) = am1(iaNegative(1:iNeg))

      iJ = iJ+1

      ! calculate the coefficients for j odd
      alpha = (iJ-1)/2
      raBeta = daX(iaNegative(1:iNeg))
      raA = raBeta * am1(iaNegative(1:iNeg)) + alpha * am2(iaNegative(1:iNeg))
      raB = raBeta * bm1(iaNegative(1:iNeg)) + alpha * bm2(iaNegative(1:iNeg))
      am2(iaNegative(1:iNeg)) = am1(iaNegative(1:iNeg)) / raB
      bm2(iaNegative(1:iNeg)) = bm1(iaNegative(1:iNeg)) / raB
      am1(iaNegative(1:iNeg)) = raA / raB
      bm1(iaNegative(1:iNeg)) = 1.0D0
      oldf(iaNegative(1:iNeg)) = f(iaNegative(1:iNeg))
      f(iaNegative(1:iNeg)) = am1(iaNegative(1:iNeg))

      iJ = iJ+1
           
      ! check for convergence, any fails go back to 20
      iFound = -1
      CALL FindInVector(iNeg,abs(f(iaNegative(1:iNeg))-oldf(iaNegative(1:iNeg))) > 1.0d-14*abs(f(iaNegative(1:iNeg))),npos,iaPos)
      IF (npos > 0) THEN
        iFound = +1
        GOTO 21
      END IF
!      DO iK = 1,iNeg
!        iFr = iaNegative(iK)
!        IF (abs(f(iFr)-oldf(iFr)) > 1.0d-14*abs(f(iFr))) THEN
!          iFound = +1
!          GOTO 21
!        END IF
!      END DO
  21  CONTINUE

      IF (iFound > 0) THEN
        GOTO 20
      END IF

      daY(iaNegative(1:iNeg)) =  exp(-daX(iaNegative(1:iNeg))) * f(iaNegative(1:iNeg))
      !!! daY is really complex, but ignore the imag part :
      !!!     - i*pi*((real(xk)<0)&(imag(xk)==0))
      raY(iaNegative(1:iNeg)) = sngl(daY(iaNegative(1:iNeg)))

      DEALLOCATE (raBeta, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Error while deallocating raBeta ***"
      DEALLOCATE (raA, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Error while deallocating raA ***"
      DEALLOCATE (raB, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Error while deallocating raB ***"

    END IF

    RETURN
    end SUBROUTINE expint

!************************************************************************
! calls function to do expintegral raY = E(n,raX), ala Numerical Recipes
! this should be faster than the Matlab version
    SUBROUTINE expintfast3(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI
    INTEGER :: iaPos(kMaxPts),npos

    CALL FindInVector(kMaxPts,raX < 0,npos,iaPos)
    IF (npos > 0) THEN
      write(kStdWarn,'(I5)') 'found ',npos,' elements less than zero in expintfast3, reset to 0 '
      raX(iaPos(1:npos)) = 0.0
    END IF

    raY = rafastExpint(3,raX)

    RETURN
    end SUBROUTINE expintfast3

!************************************************************************
! this is fast exponential integrals, from Numerical recipes
! called by eg SUBROUTINE expintfast3(n,raX)

    FUNCTION rafastExpint(n,raX)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    INTEGER :: n
    REAL :: raX(kMaxPts)
! output vars
    REAL :: rafastExpint(kMaxPts)

! local vars
    INTEGER :: MaxIT
    REAL :: Euler,FPmin

    INTEGER :: i,ii,nm1,iStop
    REAL :: r(kMaxPts),fact(kMaxPts),del(kMaxPts),psi(kMaxPts),a(kMaxPts),b(kMaxPts),c(kMaxPts),d(kMaxPts),h(kMaxPts)
    INTEGER :: iaPos(kMaxPts),npos,iaPosx(kMaxPts),nposx
    REAL :: eps

    eps = 1.0e-28
       
    r      = 0.0

    MaxIT = 100
    Euler = 0.5772156649
    FPMin = 1.0e-30
!    Eps   = 1.0e-7

    nm1 = n-1
    IF (n == 0) THEN
      r = 1.0e28
      CALL FindInVector(kMaxPts,raX > eps,npos,iaPos)
      IF (npos > 0) THEN
        r(iaPos(1:npos)) = exp(-raX(iaPos(1:npos)))/raX(iaPos(1:npos))
      END IF
    ELSEIF (n > 0) THEN
      CALL FindInVector(kMaxPts,raX < eps,npos,iaPos)
      IF (npos > 0) THEN
        r(iaPos(1:npos)) = 1.0/nm1
      END IF     !! raX < eps

      CALL FindInVector(kMaxPts,raX > 1.0,npos,iaPos)
      IF (npos > 0) THEN
        b(iaPos(1:nPos)) = raX(iaPos(1:nPos)) + n
        c(iaPos(1:nPos)) = 1/FPmin
        d(iaPos(1:nPos)) = 1.0/b(iaPos(1:nPos))
        h(iaPos(1:nPos)) = d(iaPos(1:nPos))
        iStop = -1
        i = 0
        DO WHILE ((i < MAXIT) .AND. (iStop < 0))
          i = i + 1
          a(iaPos(1:nPos)) = -i*(nm1+i)
          b(iaPos(1:nPos)) = b(iaPos(1:nPos)) + 2.0
          d(iaPos(1:nPos)) = 1/(a(iaPos(1:nPos))*d(iaPos(1:nPos))+b(iaPos(1:nPos)))
          c(iaPos(1:nPos)) = b(iaPos(1:nPos)) + a(iaPos(1:nPos))/c(iaPos(1:nPos))
          del(iaPos(1:nPos)) = c(iaPos(1:nPos))*d(iaPos(1:nPos))
          h(iaPos(1:nPos)) = h(iaPos(1:nPos))*del(iaPos(1:nPos))

          CALL FindInVector(npos,abs(del(iaPos(1:nPos))-1) < eps,nposx,iaPosx)
          IF (nposx == npos) THEN
            !! all done
            r(iaPos(1:nPos)) = h(iaPos(1:nPos))*exp(-raX(iaPos(1:nPos)))
            iStop = +1
          END IF
        END DO
      END IF      !! raX > 1.0

      CALL FindInVector(kMaxPts,(raX > eps .AND. raX <= 1.0),npos,iaPos)
      IF (npos > 0) THEN
        !get first term ... wow
        IF (nm1 /= 0) THEN
          r(iaPos(1:nPos)) = 1.0/nm1
        ELSE
          r(iaPos(1:nPos)) = -log(raX(iaPos(1:nPos))) - Euler
        END IF

        fact(iaPos(1:nPos)) = 1.0
        iStop = -1
        i = 0
        DO WHILE ((i < MAXIT) .AND. (iStop < 0))
          i = i + 1
          fact(iaPos(1:nPos)) = fact(iaPos(1:nPos)) * (-raX(iaPos(1:nPos))/i)
          IF (i /= nm1) THEN
            del(iaPos(1:nPos)) = -fact(iaPos(1:nPos))/(i-nm1)
          ELSE
            psi(iaPos(1:nPos)) = -euler
            DO ii = 1,nm1
              psi(iaPos(1:nPos)) = psi(iaPos(1:nPos)) + 1.0/ii
            END DO
            del(iaPos(1:nPos)) = fact*(-log(raX(iaPos(1:nPos))) + psi(iaPos(1:nPos)))
          END IF
          r(iaPos(1:nPos)) = r(iaPos(1:nPos)) + del(iaPos(1:nPos))

          CALL FindInVector(npos,abs(del(iaPos(1:nPos))) < eps*abs(r(iaPos(1:nPos))),nposx,iaPosx)
          IF (nposx == npos) THEN
            !! all done
            iStop = +1
          END IF

        END DO
      END IF    !! raX > eps .AND. raX <= 1.0
    END IF      !! if n == 0
      
    rafastExpint = r

    RETURN
    end FUNCTION rafastExpint

!************************************************************************
! does expintegral raY = E3(raX), using fifth order polynom fits to expint3
    SUBROUTINE expintsuperfast3_6(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI,i,j,iWhich
    REAL :: raaPout(5,6)
    INTEGER :: iaPos(kMaxPts),npos,iaPosx(kMaxPts),nposx

! this is for x > 0.00 and x <= 0.10
    DATA ((raaPout(i,j),j=1,6),i=1,1)/ &
    &   4.999997490759554e-01, -9.952726790169388e-01,  2.469512034408480e+00, &
    -1.666132334280038e+01,  1.120106392322948e+02, -3.420299925382391e+02/

! this is for x > 0.10 and x <= 1.00
    DATA ((raaPout(i,j),j=1,6),i=2,2)/ &
    &   4.964930888860016e-01,  -9.002236958326287e-01,  1.058882468914909e+00, &
    -9.590227134465121e-01,   5.595049802866501e-01, -1.460130149043994e-01/

! this is for x > 1.00 and x <= 6.00
    DATA ((raaPout(i,j),j=1,6),i=3,3)/ &
    &   9.348983139232728e-03,   8.153540216831002e-01, -1.318262705075364e+00, &
    &   1.507590471592263e+00,  -9.895617788618024e-01,  2.739113519528383e-01/

! this is for x > 6.00 and x <= 10.00
    DATA ((raaPout(i,j),j=1,6),i=4,4)/ &
    &   3.700147480617722e-04,   9.814904610819430e-01, -2.593611383300170e+00, &
    &   6.727199620937857e+00,  -1.257888244164997e+01,  1.144472478235488e+01/

! this is for x > 10.00 and x <= 20.00
    DATA ((raaPout(i,j),j=1,6),i=5,5)/ &
    &    3.360502391024836e-05,  9.969795018873142e-01, -2.884279567510683e+00, &
    &    9.510785069482322e+00, -2.619205947978609e+01,  3.863540676528253e+01/

    raY = 0.0

    IF (minval(raX) < 0.0) THEN
      write(kStdErr,*) 'in expintsuperfast3_6 found minval(raX) < 0'
      write(kStdErr,*) 'cannot have negative values, resetting to 0!'
      WHERE (raX < 0)
        raX = 0.0
      END WHERE
    END IF

    CALL FindInVector(kMaxPts,raX <= 0.1,npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 1
      DO j = 1,6
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)*(raX(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 0.1 .AND. raX <= 1.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 2
      DO j = 1,6
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)*(raX(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 1.0 .AND. raX <= 6.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 3
      DO j = 1,6
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 6.0 .AND. raX <= 10.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 4
      DO j = 1,6
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 10.0 .AND. raX <= 20.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 5
      DO j = 1,6
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 20.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 6
      raY(iaPos(1:npos)) = 0.0
    END IF

    RETURN
    end SUBROUTINE expintsuperfast3_6

!************************************************************************
! does expintegral raY = E3(raX), using 2 - 5 order polynom fits to expint3
    SUBROUTINE expintsuperfast3(raX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raX(kMaxPts)
! output vars
    REAL :: raY(kMaxPts)

! local vars
    INTEGER :: iI,i,j,iWhich,iaMax(5),iCnt
    REAL :: raaPout(5,6)
    INTEGER :: iaPos(kMaxPts),npos,iaPosx(kMaxPts),nposx

    DATA (iaMax(i),i=1,5) /3, 6, 6, 3, 3/

! this is for x > 0.00 and x <= 0.10
    DATA ((raaPout(i,j),j=1,3),i=1,1)/ &
    &   4.999958864847588e-01, -9.704714999886551e-01,  1.357944294418206e+00/

! this is for x > 0.10 and x <= 1.00
    DATA ((raaPout(i,j),j=1,6),i=2,2)/ &
    &   4.964930888860016e-01,  -9.002236958326287e-01,  1.058882468914909e+00, &
    -9.590227134465121e-01,   5.595049802866501e-01, -1.460130149043994e-01/

! this is for x > 1.00 and x <= 6.00
    DATA ((raaPout(i,j),j=1,6),i=3,3)/ &
    &   9.348983139232728e-03,   8.153540216831002e-01, -1.318262705075364e+00, &
    &   1.507590471592263e+00,  -9.895617788618024e-01,  2.739113519528383e-01/

! this is for x > 6.00 and x <= 10.00
    DATA ((raaPout(i,j),j=1,3),i=4,4)/ &
    &   6.850600578733759e-03,   8.121694943009118e-01, -9.875485935560442e-01/

! this is for x > 10.00 and x <= 20.00
    DATA ((raaPout(i,j),j=1,3),i=5,5)/ &
    &   1.873706336320716e-03,   9.115323790575932e-01, -1.489123902774130e+00/

    raY = 0.0

    IF (minval(raX) < 0.0) THEN
      write(kStdErr,*) 'in expintsuperfast3 found minval(raX) < 0'
      write(kStdErr,*) 'cannot have negative values, resetting to 0!'
      WHERE (raX < 0)
        raX = 0.0
      END WHERE
    END IF

    CALL FindInVector(kMaxPts,raX <= 0.1,npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 1
      iCnt = iaMax(iWhich)
      do j = 1,iCnt
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)*(raX(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 0.1 .AND. raX <= 1.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 2
      iCnt = iaMax(iWhich)
      do j = 1,iCnt
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)*(raX(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 1.0 .AND. raX <= 6.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 3
      iCnt = iaMax(iWhich)
      do j = 1,iCnt
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 6.0 .AND. raX <= 10.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 4
      iCnt = iaMax(iWhich)
      do j = 1,iCnt
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 10.0 .AND. raX <= 20.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 5
      iCnt = iaMax(iWhich)
      do j = 1,iCnt
          raY(iaPos(1:npos)) = raY(iaPos(1:npos)) + raaPout(iWhich,j)/(rax(iaPos(1:npos))**(j-1))
      END DO
    END IF

    CALL FindInVector(kMaxPts,(raX > 20.0),npos,iaPos)
    IF (npos >= 1) THEN
      iWhich = 6
      raY(iaPos(1:npos)) = 0.0
    END IF

    RETURN
    end SUBROUTINE expintsuperfast3

!************************************************************************
! does expintegral raY = E3(raX), using 2 - 5 order polynom fits to expint3
    SUBROUTINE expintsuperfast3matrix(iL,raaX,raY)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars
    REAL :: raaX(kMaxPts,*)
    INTEGER :: iL
! output vars
    REAL :: raY(kMaxPts)

! local vars
    REAL :: raX(kMaxPts)

    raX = raaX(:,iL)
    raY = 0.0
    CALL expintsuperfast3(raX,raY)

    RETURN
    end SUBROUTINE expintsuperfast3matrix

!************************************************************************
END MODULE rad_angles
