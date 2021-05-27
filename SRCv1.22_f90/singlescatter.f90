MODULE singlescatter

USE clear_scatter_misc
USE basic_common
USE ttorad_common

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine calculates the scattered solar intensity at the bottom of
! each layer for a downlook instr
! pretty much the same as  SUBROUTINE Solar in rad_misc.f
! except it does the [w P exp(1-exp(-x))] stuff at the end!

! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))

! if iSolarRadOrJac == +1 ==> compute stuff for rads for rad_DOWN_pclsam_solar
! if iSolarRadOrJac == -1 ==> compute stuff for jacs for
!                                                    AddSolarScatterGasJacobian
    SUBROUTINE SolarScatterIntensity_Downlook( &
    iDoSolar,raFreq,iaCldLayer, &
    raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth, &
    iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
    iTag,iSolarRadOrJac,raaSolarScatterX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! input vars
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raLayAngles(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag,iDoSolar
    INTEGER :: iaCldLayer(kProfLayer),iSolarRadOrJac
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows), &
    raaAsym(kMaxPts,kMixFilRows)
    REAL :: rFracTop,rFracBot,rSatAzimuth,rSolAzimuth
! output variable
    REAL :: raaSolarScatterX(kMaxPts,kProfLayer)

! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),raU(kMaxPts),muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: muSat,raKabs(kMaxPts),rSilly
    INTEGER :: iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: raNoScale(kMaxPts)

    raNoScale = 1.0

    IF (iDoSolar == 0) THEN
      !use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      !read in data from file
      CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
     raSun  = raSun*muSun*rOmegaSun
    raKAbs = 0.0

! note raExtraSun is initialized to all zeros
    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
      iaRadLayerTemp,iT,iExtraSun,raExtraSun)

! note how raaSolarScatter is calculated at the TOP of the layer!!!
! ie by computing raaSolarScatter before updating raKAbs, we do solar
!    radiance at the top-of-the-layer
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
      ! the current defined atmosphere used all Gnd-kProfLayer layers
      DO iLay=iNumLayer,2,-1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs + raaExt(:,iL)/muSun*raNoScale
      END DO
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
        !!!solar intensity at GND
        raSun = raSun*exp(-raKAbs)
      END DO
      raExtraSun = 0.0

    ELSE IF (iExtraSun > 0) THEN
      ! all upper layers not used eg instrument could be on a low flying aircraft
      IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
        write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
        write(kStdWarn,*)'but posn of instrument is at middle of '
        write(kStdWarn,*)'layer ==> need to add extra term'
        !first do the highest layer .. make it "full"
        iI = iNumLayer
        write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
        DO iLay=iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
          raExtraSun = raSun*exp(-rakAbs)
        END DO
        !now do remaining layers, all the way to the ground-1
        DO iLay=iNumLayer-1,2,-1
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKAbs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        DO iLay=1,1
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKAbs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
          !!!solar intensity at GND
          raSun = raSun*exp(-raKAbs)
        END DO
      END IF
         
      IF (iT > iNumLayer) THEN
        write(kStdWarn,*)'need to do the upper layers as well!!'
        !now do top layers, all the way to the instrument
        DO iLay=iT,iNumLayer+1,-1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        !now do the layer instrument is in
        DO iLay=iNumLayer,iNumLayer
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
          raExtraSun = raSun*(exp(-raKabs))
        END DO
        !now do all the way to the ground-1
        DO iLay=iNumLayer-1,2,-1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        !now do ground
        DO iLay=1,1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
          !!!solar intensity at GND
          raSun = raSun*exp(-raKAbs)
        END DO
      END IF
         
    END IF

!!! ---------------->                     <------------------------
!!! ---------------->                     <------------------------
!!! now do the actual scattering intensity computation
    IF (iSolarRadOrJac == +1) THEN  !!! do stuff for rads
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        muSun = cos(raSunAngles(iL)*kPi/180.0)
        IF (iaCldLayer(iL) == 1) THEN
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raU   = raaExt(:,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
          raU   = 1.0-exp(-raU*raNoScale)
          raU   = raU * muSun * raaSSAlb(:,iL)/kForP/(muSat+muSun)
          raU   = raU * rahg2_azimuth_real(-muSun,muSat,rSolAzimuth,rSatAzimuth,raaAsym(:,iL))
          raaSolarScatterX(:,iL) = raaSolarScatterX(:,iL)*raU
        ELSE
          ! set contribution to solar scattered radiation to zero
          raaSolarScatterX(:,iL) = 0.0
        END IF
      END DO

    ELSEIF (iSolarRadOrJac == -1) THEN  !!!
      !accumulate scattered solar rad. Not used anywhere ...
      !initialize raSun == cumulative solar scattered contribution
      raSun = 0.0

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        muSun = cos(raSunAngles(iL)*kPi/180.0)
        IF (iaCldLayer(iL) == 1) THEN
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raU   = raaExt(:,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
          raU   = 1.0-exp(-raU*raNoScale)
          raU   = raU * muSun*raaSSAlb(:,iL)/kForP/(muSat+muSun)
          raU   = raU * rahg2_azimuth_real(-muSun,muSat,rSolAzimuth,rSatAzimuth,raaAsym(:,iL))
          !! this is contribution from current layer
          raaSolarScatterX(:,iL) = raaSolarScatterX(:,iL)*raU
          !! add on effects from previous layer(s)
          raU   = exp(-raaExt(:,iL)/abs(muSat))
          raaSolarScatterX(:,iL) = raSun*raU+raaSolarScatterX(:,iL)
        ELSEIF (iaCldLayer(iL) == -1) THEN
          ! just use attenuated cumulative solar scatter from previous calc
          raU   = exp(-raaExt(:,iL)/abs(muSat))
          raaSolarScatterX(:,iL) = raSun*raU
        END IF

        ! update raSun
        raSun = raaSolarScatterX(:,iL)
      END DO

    END IF

    RETURN
    end SUBROUTINE SolarScatterIntensity_Downlook
      
!************************************************************************
! this function computes the Henyey Greenstein function, assuming the
! cos(phi1-phi2) factor = 1
    REAL Function hg2_azimuth_real(mu1,mu2,phi1,phi2,g)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry
    REAL :: phi1,phi2

    REAL :: normB,mu0,yexact,delphi

    delphi = (phi1-phi2)*kPi/180.0

! normB is normalisation of mu from -1 to 1 and works out to be 2.0
!! we also know that (1/2) integral P(-1,1) = 1
! ormB = 1/sqrt(1+g*g - 2.0*g) - 1/sqrt(1+g*g + 2.0*g)
! ormB = (1-g*g)/g * normB
! ormB = 2.0

!!!compute mu0 = cos ofangle between the two
    mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)*cos(delphi)
     
    yexact = (1 + g*g - 2.0*g*mu0) * sqrt(1 + g*g - 2.0*g*mu0)
    yexact = (1-g*g)/yexact
! exact = yexact/normB * 2.0
    yexact = yexact

    hg2_azimuth_real = yexact

    RETURN
    end Function hg2_azimuth_real

!************************************************************************
! this function computes the Henyey Greenstein function, assuming the
! cos(phi1-phi2) factor = 1
    function rahg2_azimuth_real(mu1,mu2,phi1,phi2,raG)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: mu1,mu2,raG(kMaxPts)       !mu1,mu2 are the two angles, g is the asymmetry
    REAL :: phi1,phi2
    REAL :: rahg2_azimuth_real(kMaxPts)

    REAL :: normB,mu0,yexact(kMaxPts),delphi

    delphi = (phi1-phi2)*kPi/180.0

! normB is normalisation of mu from -1 to 1 and works out to be 2.0
!! we also know that (1/2) integral P(-1,1) = 1
! ormB = 1/sqrt(1+g*g - 2.0*g) - 1/sqrt(1+g*g + 2.0*g)
! ormB = (1-g*g)/g * normB
! ormB = 2.0

!!!compute mu0 = cos ofangle between the two
    mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)*cos(delphi)
     
    yexact = (1 + raG*raG - 2.0*raG*mu0) * sqrt(1 + raG*raG - 2.0*raG*mu0)
    yexact = (1-raG*raG)/yexact
! exact = yexact/normB * 2.0
    yexact = yexact

    rahg2_azimuth_real = yexact

    RETURN
    end Function rahg2_azimuth_real

!************************************************************************
! each layer for an uplook instr, mainly for use in Jacobian code
! almost same as SUBROUTINE SolarScatterIntensity_Downlook

! if iSolarRadOrJac == +1 ==> compute stuff for rads for rad_DOWN_pclsam_solar
! if iSolarRadOrJac == -1 ==> compute stuff for jacs for
!                                                    AddSolarScatterGasJacobian
    SUBROUTINE SolarScatterIntensity_Uplook( &
    iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer, &
    iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
    iTag,iSolarRadOrJac,raaSolarScatterX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! input vars
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raLayAngles(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag,iDoSolar
    INTEGER :: iaCldLayer(kProfLayer),iSolarRadOrJac
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows), &
    raaAsym(kMaxPts,kMixFilRows)
    REAL :: rFracTop,rFracBot
! output variable
    REAL :: raaSolarScatterX(kMaxPts,kProfLayer)

! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raTau(kMaxPts),raSun(kMaxPts),rU,muSun,raTauSum(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,rSunAngle,rFrac
    REAL :: muSat,raKabs(kMaxPts),rSilly
    INTEGER :: iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay,iLow
    REAL :: raSolarScatter(kMaxPts),raNoScale(kMaxPts)

    IF (iDoSolar == 0) THEN
      !use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      !read in data from file
      CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = 1
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)

! now adjust raSun by cos(rSunAngle) * rSolidAngle
    raSun = raSun*muSun*rOmegaSun

! need cumulative optical depth to TOP, not BOTTOM, of the layer, so initialize
! as follows
    iLay = 1
    iL  = iaRadLayer(iLay)
    raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
    raTauSum = -raaExt(:,iL)*rFracTop*raNoScale

    iLow = iNumLayer
    DO iLay = 1,iLow
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      muSun    = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      IF (iLay == 1) THEN
        rFrac = rFracTop
      ELSE IF (iLay == iNumLayer) THEN
        rFrac = rFracBot
      ELSE
        rFrac = 1.0
      END IF
      !cumulative optical depth to TOP, not BOTTOM, of the layer
      raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
      raTauSum  = raTauSum + raaExt(:,iL)*rFrac*raNoScale
      !!! now see if we need the solar scatter term
      IF (iaCldLayer(iLay) == 1) THEN
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raTau = raaExt(:,iL)*rFrac*raNoScale
        raSolarScatter = rahg2_real(-muSun,-muSat,raaAsym(:,iL)) * &
            (exp(-raTau/muSat)-exp(-raTau/muSun))
        raSolarScatter = raSolarScatter*raaSSAlb(:,iL)
        raSolarScatter = raSolarScatter*raSun*exp(-raTauSum/muSun)
        raSolarScatter = raSolarScatter*muSun/(muSat-muSun)/kForP
        raaSolarScatterX(:,iL) = raSolarScatter
      END IF
    END DO

    RETURN
    end SUBROUTINE SolarScatterIntensity_Uplook

!************************************************************************
! this subroutine calculates the scattered solar intensity at the bottom of
! each layer for a downlook instr

! used for debugging mysterious N2O discrepancies with SARTA daytime!

! pretty much the same as  SUBROUTINE Solar in rad_misc.f
! except it does the [w P exp(1-exp(-x))] stuff at the end!

! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))

! if iSolarRadOrJac == +1 ==> compute stuff for rads for rad_DOWN_pclsam_solar
! if iSolarRadOrJac == -1 ==> compute stuff for jacs for
!                                                    AddSolarScatterGasJacobian
    SUBROUTINE xSolarScatterIntensity_Downlook_debug( &
    iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer, &
    iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
    iTag,iSolarRadOrJac,raaSolarScatterX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! input vars
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raLayAngles(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag,iDoSolar
    INTEGER :: iaCldLayer(kProfLayer),iSolarRadOrJac
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows), &
    raaAsym(kMaxPts,kMixFilRows)
    REAL :: rFracTop,rFracBot
! output variable
    REAL :: raaSolarScatterX(kMaxPts,kProfLayer)

! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),raU(kMaxPts),muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: muSat,raKabs(kMaxPts),rSilly
    INTEGER :: iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: raNoScale(kMaxPts)

    raNoScale = 1.0

    IF (iDoSolar == 0) THEN
      !use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      !read in data from file
      CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    raSun  = raSun*muSun*rOmegaSun
    raKAbs = 0.0

! note raExtraSun is initialized to all zeros
    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
      iaRadLayerTemp,iT,iExtraSun,raExtraSun)

! note how raaSolarScatter is calculated at the TOP of the layer!!!
! ie by computing raaSolarScatter before updating raKAbs, we do solar
!    radiance at the top-of-the-layer
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
      ! the current defined atmosphere used all Gnd-kProfLayer layers
      DO iLay=iNumLayer,2,-1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs + raaExt(:,iL)/muSun*raNoScale
      END DO
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
        !!!solar intensity at GND
        raSun = raSun*exp(-raKAbs)
      END DO
      raExtraSun = 0.0
    END IF

!!! ---------------->                     <------------------------
!!! ---------------->                     <------------------------
!!! now do the actual scattering intensity computation
    IF (iSolarRadOrJac == +1) THEN  !!! do stuff for rads
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        muSun = cos(raSunAngles(iL)*kPi/180.0)
        IF (iaCldLayer(iL) == 1) THEN
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raU   = raaExt(:,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
          raU   = 1.0 - exp(-raU*raNoScale)
          raU   = muSun/(muSat+muSun)/kForP * raaSSAlb(:,iL) * raU
          raU   = raU * rahg2_real(-muSun,muSat,raaAsym(:,iL))
          raaSolarScatterX(:,iL) = raaSolarScatterX(:,iL)*raU
        ELSE
          ! set contribution to solar scattered radiation to zero
          raaSolarScatterX(:,iL) = 0.0
        END IF
      END DO
    END IF

    RETURN
    end SUBROUTINE xSolarScatterIntensity_Downlook_debug
      
!************************************************************************
! this subroutine calculates the scattered solar intensity at the bottom of
! each layer for a downlook instr
! pretty much the same as  SUBROUTINE Solar in rad_misc.f
! except it does the [w P exp(1-exp(-x))] stuff at the end!

! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))

! if iSolarRadOrJac == +1 ==> compute stuff for rads for rad_DOWN_pclsam_solar
! if iSolarRadOrJac == -1 ==> compute stuff for jacs for
!                                                    AddSolarScatterGasJacobian
    SUBROUTINE xSolarScatterIntensity_DownlookXWORKS( &
    iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer, &
    iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
    iTag,iSolarRadOrJac,raaSolarScatterX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! input vars
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raLayAngles(kProfLayer),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag,iDoSolar
    INTEGER :: iaCldLayer(kProfLayer),iSolarRadOrJac
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows), &
    raaAsym(kMaxPts,kMixFilRows)
    REAL :: rFracTop,rFracBot
! output variable
    REAL :: raaSolarScatterX(kMaxPts,kProfLayer)

! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),raU(kMaxPts),muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: muSat,raKabs(kMaxPts),rSilly
    INTEGER :: iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: raNoScale(kMaxPts)

    raNoScale = 1.0

    IF (iDoSolar == 0) THEN
      !use 5700K
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      rSunTemp = kSunTemp
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      !read in data from file
      CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    raSun  = raSun*muSun*rOmegaSun
    raKAbs = 0.0

! note raExtraSun is initialized to all zeros
    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
      iaRadLayerTemp,iT,iExtraSun,raExtraSun)

! note how raaSolarScatter is calculated at the TOP of the layer!!!
! ie by computing raaSolarScatter before updating raKAbs, we do solar
!    radiance at the top-of-the-layer
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
      ! the current defined atmosphere used all Gnd-kProfLayer layers
      DO iLay = iNumLayer,2,-1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs + raaExt(:,iL)/muSun*raNoScale
      END DO
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
        raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
        raKAbs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
        !!!solar intensity at GND
        raSun = raSun*exp(-raKAbs)
      END DO
      raExtraSun = 0.0

    ELSE IF (iExtraSun > 0) THEN
      ! all upper layers not used eg instrument could be on a low flying aircraft
      IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
        write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
        write(kStdWarn,*)'but posn of instrument is at middle of '
        write(kStdWarn,*)'layer ==> need to add extra term'
        !first do the highest layer .. make it "full"
        iI=iNumLayer
        write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
        DO iLay=iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
          raExtraSun = raSun*exp(-rakAbs)
        END DO
        !now do remaining layers, all the way to the ground-1
        DO iLay=iNumLayer-1,2,-1
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKAbs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        DO iLay=1,1
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKAbs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
          !!!solar intensity at GND
          raSun = raSun*exp(-raKAbs)
        END DO
      END IF
         
      IF (iT > iNumLayer) THEN
        write(kStdWarn,*)'need to do the upper layers as well!!'
        !now do top layers, all the way to the instrument
        DO iLay=iT,iNumLayer+1,-1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        !now do the layer instrument is in
        DO iLay=iNumLayer,iNumLayer
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
          raExtraSun = raSun*(exp(-raKabs))
        END DO
        !now do all the way to the ground-1
        DO iLay=iNumLayer-1,2,-1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)/muSun*raNoScale
        END DO
        !now do ground
        DO iLay=1,1
          iL = iaRadLayerTemp(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raaSolarScatterX(:,iL) = raSun*exp(-raKAbs)
          raKabs = raKAbs+raaExt(:,iL)*rFracBot/muSun*raNoScale
          !!!solar intensity at GND
          raSun = raSun*exp(-raKAbs)
        END DO
      END IF         
    END IF

!!! ---------------->                     <------------------------
!!! ---------------->                     <------------------------
!!! now do the actual scattering intensity computation
    IF (iSolarRadOrJac == +1) THEN  !!! do stuff for rads
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        muSun = cos(raSunAngles(iL)*kPi/180.0)
        IF (iaCldLayer(iL) == 1) THEN
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raU   = raaExt(:,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
          raU   = 1.0-exp(-raU*raNoScale)
          raU   = raU * muSun * raaSSAlb(:,iL)/kForP/(muSat+muSun)
          raU   = raU * rahg2_real(-muSun,muSat,raaAsym(:,iL))
          raaSolarScatterX(:,iL) = raaSolarScatterX(:,iL)*raU
        ELSE
          ! set contribution to solar scattered radiation to zero
          raaSolarScatterX(:,iL) = 0.0
        END IF
      END DO

    ELSEIF (iSolarRadOrJac == -1) THEN  !!!
      ! accumulate scattered solar rad. Not used anywhere ...
      ! initialize raSun == cumulative solar scattered contribution
      raSun = 0.0

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        muSun = cos(raSunAngles(iL)*kPi/180.0)
        IF (iaCldLayer(iL) == 1) THEN
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raU   = raaExt(:,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
          raU   = 1.0 - exp(-raU*raNoScale)
          raU   = raU * muSun*raaSSAlb(:,iL)/kForP/(muSat+muSun)
          raU   = raU * rahg2_real(-muSun,muSat,raaAsym(:,iL))
          !! this is contribution from current layer
          raaSolarScatterX(:,iL) = raaSolarScatterX(:,iL)*raU
          !! add on effects from previous layer(s)
          raU   = exp(-raaExt(:,iL)/abs(muSat))
          raaSolarScatterX(:,iL) = raSun*raU+raaSolarScatterX(:,iL)
        ELSEIF (iaCldLayer(iL) == -1) THEN
          ! just use attenuated cumulative solar scatter from previous calc
          raU   = exp(-raaExt(:,iL)/abs(muSat))
          raaSolarScatterX(:,iL) = raSun*raU
        END IF

        ! update raSun
        raSun = raaSolarScatterX(:,iL)
      END DO

    END IF

    RETURN
    end SUBROUTINE xSolarScatterIntensity_DownlookXWORKS
      
!************************************************************************
END MODULE singlescatter
