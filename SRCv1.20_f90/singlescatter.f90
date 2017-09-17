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
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),rU,muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: ttorad,muSat,raKabs(kMaxPts),hg2_azimuth_real,rSilly
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: rNoScale

    rNoScale = 1.0

    IF (iDoSolar == 0) THEN
    ! se 5700K
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
    ELSEIF (iDoSolar == 1) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
    ! ead in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    DO iFr=1,kMaxPts
        raSun(iFr)  = raSun(iFr)*muSun*rOmegaSun
        raKAbs(iFr) = 0.0
    END DO

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
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr) + raaExt(iFr,iL)/muSun*rNoScale
            END DO
        END DO
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
            !!!solar intensity at GND
                raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
            END DO
        END DO
        DO iFr=1,kMaxPts
            raExtraSun(iFr) = 0.0
        END DO

    ELSE IF (iExtraSun > 0) THEN
    ! all upper layers not used eg instrument could be on a low flying aircraft
        IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
            write(kStdWarn,*)'but posn of instrument is at middle of '
            write(kStdWarn,*)'layer ==> need to add extra term'
        ! irst do the highest layer .. make it "full"
            iI=iNumLayer
            write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
            DO iLay=iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                    raExtraSun(iFr) = raSun(iFr)*exp(-rakAbs(iFr))
                END DO
            END DO
        ! ow do remaining layers, all the way to the ground-1
            DO iLay=iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
            DO iLay=1,1
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
                !!!solar intensity at GND
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
        IF (iT > iNumLayer) THEN
            write(kStdWarn,*)'need to do the upper layers as well!!'
        ! ow do top layers, all the way to the instrument
            DO iLay=iT,iNumLayer+1,-1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
        ! ow do the layer instrument is in
            DO iLay=iNumLayer,iNumLayer
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                    raExtraSun(iFr) = raSun(iFr)*(exp(-raKabs(iFr)))
                END DO
            END DO
        ! ow do all the way to the ground-1
            DO iLay=iNumLayer-1,2,-1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
        ! ow do ground
            DO iLay=1,1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
                !!!solar intensity at GND
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
    END IF

!      DO iLay = 1,iNumLayer
!        iL = iaRadLayer(iLay)
!        print *,iLay,iL,iaCldLayer(iL),raaSolarScatterX(1,iL)
!        END DO
!      CALL DOStop

!!! ---------------->                     <------------------------
!!! ---------------->                     <------------------------
!!! now do the actual scattering intensity computation
    IF (iSolarRadOrJac == +1) THEN  !!! do stuff for rads
        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            muSun = cos(raSunAngles(iL)*kPi/180.0)
            IF (iaCldLayer(iL) == 1) THEN
                DO iFr = 1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    rU   = raaExt(iFr,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
                    rU   = 1.0-exp(-rU*rNoScale)
                    rU   = rU * muSun * raaSSAlb(iFr,iL)/kForP/(muSat+muSun)
                    rU   = rU * hg2_azimuth_real(-muSun,muSat, &
                    rSolAzimuth,rSatAzimuth,raaAsym(iFr,iL))
                    raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)*rU
                END DO
            ELSE
            ! set contribution to solar scattered radiation to zero
                DO iFr = 1,kMaxPts
                    raaSolarScatterX(iFr,iL) = 0.0
                END DO
            END IF
        END DO

    ELSEIF (iSolarRadOrJac == -1) THEN  !!!
    ! accumulate scattered solar rad. Not used anywhere ...
    ! nitialize raSun == cumulative solar scattered contribution
        DO iFr = 1,kMaxPts
            raSun(iFr) = 0.0
        END DO

        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            muSun = cos(raSunAngles(iL)*kPi/180.0)
            IF (iaCldLayer(iL) == 1) THEN
                DO iFr = 1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    rU   = raaExt(iFr,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
                    rU   = 1.0-exp(-rU*rNoScale)
                    rU   = rU * muSun*raaSSAlb(iFr,iL)/kForP/(muSat+muSun)
                    rU   = rU * hg2_azimuth_real(-muSun,muSat, &
                    rSolAzimuth,rSatAzimuth,raaAsym(iFr,iL))
                !! this is contribution from current layer
                    raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)*rU
                !! add on effects from previous layer(s)
                    rU   = exp(-raaExt(iFr,iL)/abs(muSat))
                    raaSolarScatterX(iFr,iL) = &
                    raSun(iFr)*rU+raaSolarScatterX(iFr,iL)
                END DO
            ELSEIF (iaCldLayer(iL) == -1) THEN
            ! just use attenuated cumulative solar scatter from previous calc
                DO iFr = 1,kMaxPts
                    rU   = exp(-raaExt(iFr,iL)/abs(muSat))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*rU
                END DO
            END IF

        ! update raSun
            DO iFr = 1,kMaxPts
                raSun(iFr) = raaSolarScatterX(iFr,iL)
            END DO
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
    REAL :: ttorad,muSat,raKabs(kMaxPts),hg2_real,rSilly
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay,iLow
    REAL :: rSolarScatter,rNoScale

    IF (iDoSolar == 0) THEN
    ! se 5700K
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
    ELSEIF (iDoSolar == 1) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
    ! ead in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = 1
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)

! now adjust raSun by cos(rSunAngle) * rSolidAngle
    DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*muSun*rOmegaSun
    END DO

! need cumulative optical depth to TOP, not BOTTOM, of the layer, so initialize
! as follows
    iLay = 1
    iL  = iaRadLayer(iLay)
    DO iFr=1,kMaxPts
        rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
        raTauSum(iFr) = -raaExt(iFr,iL)*rFracTop*rNoScale
    END DO

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
    ! umulative optical depth to TOP, not BOTTOM, of the layer
        DO iFr=1,kMaxPts
            rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
            raTauSum(iFr)  = raTauSum(iFr) + raaExt(iFr,iL)*rFrac*rNoScale
        END DO
    !!! now see if we need the solar scatter term
        IF (iaCldLayer(iLay) == 1) THEN
            DO iFr = 1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raTau(iFr) = raaExt(iFr,iL)*rFrac*rNoScale
                rSolarScatter = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) * &
                (exp(-raTau(iFr)/muSat)-exp(-raTau(iFr)/muSun))
                rSolarScatter = rSolarScatter*raaSSAlb(iFr,iL)
                rSolarScatter = rSolarScatter*raSun(iFr)*exp(-raTauSum(iFr)/muSun)
                rSolarScatter = rSolarScatter*muSun/(muSat-muSun)/kForP
                raaSolarScatterX(iFr,iL) = rSolarScatter
            END DO
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
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),rU,muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: ttorad,muSat,raKabs(kMaxPts),hg2_real,rSilly
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: rNoScale

    rNoScale = 1.0

    IF (iDoSolar == 0) THEN
    ! se 5700K
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
    ELSEIF (iDoSolar == 1) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
    ! ead in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    DO iFr=1,kMaxPts
    !        raSun(iFr)  = 1000.0
        raSun(iFr)  = raSun(iFr)*muSun*rOmegaSun
        raKAbs(iFr) = 0.0
    END DO

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
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr) + raaExt(iFr,iL)/muSun*rNoScale
            END DO
        END DO
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
            !!!solar intensity at GND
                raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
            END DO
        END DO
        DO iFr=1,kMaxPts
            raExtraSun(iFr) = 0.0
        END DO
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
                DO iFr = 1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    rU   = raaExt(iFr,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
                    rU   = 1.0 - exp(-rU*rNoScale)
                    rU   = muSun/(muSat+muSun)/kForP * raaSSAlb(iFr,iL) * rU
                    rU   = rU * hg2_real(-muSun,muSat,raaAsym(iFr,iL))
                    raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)*rU

                !              rU   = muSun/(muSat+muSun)/kForP *
                !     $               hg2_real(-muSun,muSat,raaAsym(iFr,iL))
                !c               rU   = raaSSAlb(iFr,iL)
                !c              raaSolarScatterX(iFr,iL) = rU*1000
                !              raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)
                !              if (ifr .eq. 1) then
                !                   print *,iNumLayer-iL+1,rOmegaSun,raaSolarScatterX(iFr,iL)
                !                   end if
                !              raaSolarScatterX(iFr,iL) = 1000.0
                END DO
            ELSE
            ! set contribution to solar scattered radiation to zero
                DO iFr = 1,kMaxPts
                    raaSolarScatterX(iFr,iL) = 0.0
                END DO
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
    REAL :: raExtraSun(kMaxPts),raSun(kMaxPts),rU,muSun
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: ttorad,muSat,raKabs(kMaxPts),hg2_real,rSilly
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: rNoScale

    rNoScale = 1.0

    IF (iDoSolar == 0) THEN
    ! se 5700K
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
    ELSEIF (iDoSolar == 1) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
    ! ead in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    iLay      = iNumLayer
    iL        = iaRadLayer(iLay)
    rSunAngle = raSunAngles(MP2Lay(iL))*kPi/180
    muSun     = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    DO iFr=1,kMaxPts
        raSun(iFr)  = raSun(iFr)*muSun*rOmegaSun
        raKAbs(iFr) = 0.0
    END DO

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
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr) + raaExt(iFr,iL)/muSun*rNoScale
            END DO
        END DO
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
                rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
            !!!solar intensity at GND
                raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
            END DO
        END DO
        DO iFr=1,kMaxPts
            raExtraSun(iFr) = 0.0
        END DO

    ELSE IF (iExtraSun > 0) THEN
    ! all upper layers not used eg instrument could be on a low flying aircraft
        IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
            write(kStdWarn,*)'but posn of instrument is at middle of '
            write(kStdWarn,*)'layer ==> need to add extra term'
        ! irst do the highest layer .. make it "full"
            iI=iNumLayer
            write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
            DO iLay=iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                    raExtraSun(iFr) = raSun(iFr)*exp(-rakAbs(iFr))
                END DO
            END DO
        ! ow do remaining layers, all the way to the ground-1
            DO iLay=iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
            DO iLay=1,1
                iL = iaRadLayer(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKAbs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
                !!!solar intensity at GND
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
        IF (iT > iNumLayer) THEN
            write(kStdWarn,*)'need to do the upper layers as well!!'
        ! ow do top layers, all the way to the instrument
            DO iLay=iT,iNumLayer+1,-1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
        ! ow do the layer instrument is in
            DO iLay=iNumLayer,iNumLayer
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                    raExtraSun(iFr) = raSun(iFr)*(exp(-raKabs(iFr)))
                END DO
            END DO
        ! ow do all the way to the ground-1
            DO iLay=iNumLayer-1,2,-1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)/muSun*rNoScale
                END DO
            END DO
        ! ow do ground
            DO iLay=1,1
                iL = iaRadLayerTemp(iLay)
                muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*exp(-raKAbs(iFr))
                    raKabs(iFr) = raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/muSun*rNoScale
                !!!solar intensity at GND
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
    END IF

!      DO iLay = 1,iNumLayer
!        iL = iaRadLayer(iLay)
!        print *,iLay,iL,iaCldLayer(iL),raaSolarScatterX(1,iL)
!        END DO
!      CALL DOStop

!!! ---------------->                     <------------------------
!!! ---------------->                     <------------------------
!!! now do the actual scattering intensity computation
    IF (iSolarRadOrJac == +1) THEN  !!! do stuff for rads
        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            muSun = cos(raSunAngles(iL)*kPi/180.0)
            IF (iaCldLayer(iL) == 1) THEN
                DO iFr = 1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    rU   = raaExt(iFr,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
                    rU   = 1.0-exp(-rU*rNoScale)
                    rU   = rU * muSun * raaSSAlb(iFr,iL)/kForP/(muSat+muSun)
                    rU   = rU * hg2_real(-muSun,muSat,raaAsym(iFr,iL))
                    raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)*rU
                END DO
            ELSE
            ! set contribution to solar scattered radiation to zero
                DO iFr = 1,kMaxPts
                    raaSolarScatterX(iFr,iL) = 0.0
                END DO
            END IF
        END DO

    ELSEIF (iSolarRadOrJac == -1) THEN  !!!
    ! accumulate scattered solar rad. Not used anywhere ...
    ! nitialize raSun == cumulative solar scattered contribution
        DO iFr = 1,kMaxPts
            raSun(iFr) = 0.0
        END DO

        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            muSun = cos(raSunAngles(iL)*kPi/180.0)
            IF (iaCldLayer(iL) == 1) THEN
                DO iFr = 1,kMaxPts
                    rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
                    rU   = raaExt(iFr,iL)*(1.0/abs(muSat) + 1.0/abs(muSun))
                    rU   = 1.0-exp(-rU*rNoScale)
                    rU   = rU * muSun*raaSSAlb(iFr,iL)/kForP/(muSat+muSun)
                    rU   = rU * hg2_real(-muSun,muSat,raaAsym(iFr,iL))
                !! this is contribution from current layer
                    raaSolarScatterX(iFr,iL) = raaSolarScatterX(iFr,iL)*rU
                !! add on effects from previous layer(s)
                    rU   = exp(-raaExt(iFr,iL)/abs(muSat))
                    raaSolarScatterX(iFr,iL) = &
                    raSun(iFr)*rU+raaSolarScatterX(iFr,iL)
                END DO
            ELSEIF (iaCldLayer(iL) == -1) THEN
            ! just use attenuated cumulative solar scatter from previous calc
                DO iFr = 1,kMaxPts
                    rU   = exp(-raaExt(iFr,iL)/abs(muSat))
                    raaSolarScatterX(iFr,iL) = raSun(iFr)*rU
                END DO
            END IF

        ! update raSun
            DO iFr = 1,kMaxPts
                raSun(iFr) = raaSolarScatterX(iFr,iL)
            END DO
        END DO

    END IF

    RETURN
    end SUBROUTINE xSolarScatterIntensity_DownlookXWORKS
      
!************************************************************************
