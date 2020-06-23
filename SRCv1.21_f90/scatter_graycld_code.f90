! Copyright 2007
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_graycld_code

USE basic_common
USE ttorad_common
USE jac_main
USE jac_pclsam_up
USE jac_pclsam_down
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_diff_and_quad
USE spline_and_sort_and_common

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
!************************************************************************
! stripped down version from rad_main.f
!************************************************************************
!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth

! NO NLTE allowed here!

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raaInten(kMaxPts,kProfLayer),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP,iProfileLayers
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2)
    REAL :: raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

! this is for Rayleigh
    REAL :: raaRayleigh(kMaxPts,kProfLayer)
    REAL :: raPZ(kProfLayer),raTZ(kProfLayer)

! local variables
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar
! need this
    REAL :: raInten(kMaxPts)

! for the NLTE which is not used in this routine
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),r0
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN
    REAL :: bt2rad,t2s
    INTEGER :: iFr1,troplayer
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular

! for printing
    INTEGER :: iPrintX,iMax
          
    IF ((raFreq(1) >= 10000) .AND. (raSunAngles(50) <= 90)) THEN
        write(kStdWarn,*) 'daytime downlook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV'
        write(kStdWarn,*) 'oops not yet coded or black clouds'
        CALL DoStop
    !        CALL rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV(raFreq,raaInten,raVTemp,
    !     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
    !     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
    !     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
    !     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
    !     $    raThickness,raPressLevels,iProfileLayers,pProf,
    !     $    raTPressLevels,iKnowTP,
    !     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
    !        RETURN
    END IF

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,4) == 3) rThermalRefl = 1.0   !! nick nalli
          
! calculate cos(SatAngle)
    rCos=cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),fFracBot,rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracBot,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iLay,iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm), &
        raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer iL = ',iL,' temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer iL = ',iL,' temp : orig, interp ',raVTemp(iL),raVT1(iL)

!      DO iFr = 1,100
!        print *,'cloud ',iFr,raVTemp(iFr),raVT1(iFr)
!      END DO
!      print *,'here debug cloud'
!      Call DoStop

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
    troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iMax = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iMax) THEN
            iMax = iaOp(iLay)
        END IF
    END DO

    iHigh = -1
    DO iLay=1,iNumLayer
        IF (iaRadLayer(iLay) > iHigh) THEN
            iHigh = iLay
        END IF
    END DO

    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

          
! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay)  = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO
          
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr) = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
    iNLTEStart = kProfLayer + 1
    iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
    iUpper = -1
    write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
    write (kStdWarn,*) 'stop normal radtransfer after ',iSTopNormalRadTransfer,' iterations'

    DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
    ! ormal, no LTE emission stuff
        DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    iSpecular = +1    !some specular refl, plus diffuse
    iSpecular = -1    !no   specular refl, only diffuse

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    IF (iSpecular > 0) THEN
        write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
        CALL loadspecular(raFreq,raSpecularRefl)
        DO iFr=1,kMaxPts
        ! aSpecularRefl(iFr) = 0.0272   !!! smooth water
            raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+ &
            raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
            raSun(iFr)*(raSpecularRefl(iFr) + raSunRefl(iFr))
        END DO
    ELSE
        DO iFr=1,kMaxPts
            raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+ &
            raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
            raSun(iFr)*raSunRefl(iFr)
        END DO
    END IF

!      print *,iNp
!      print *,(iaOp(iL),iL = 1,iNp)

    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        iPrintX = Inset(iL,iaOp,iNp)
        IF (iPrintX > 0) THEN
            write(kStdWarn,*) 'GND output rads at',iLay,' th rad layer = ',iL,' kcarta layer to ind ',iPrintX
            DO iFr=1,kMaxPts
                raaInten(iFr,iPrintX) = raInten(iFr)
            END DO
        END IF
    ! xxxxx            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr) = raaEmission(iFr,iLay) + &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay = 2,iHigh-1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    !        print *,iLay,iL,raLayAngles(MP2Lay(iL)),rMPTemp
        iPrintX = Inset(iL,iaOp,iNp)
        IF (iPrintX > 0) THEN
            write(kStdWarn,*) 'MID ATM output rads at',iLay,' th rad layer = ',iL,' kcarta layer to ind ',iPrintX
            DO iFr=1,kMaxPts
                raaInten(iFr,iPrintX) = raInten(iFr)
            END DO
        END IF
    ! xxxxxxx            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr) = raaEmission(iFr,iLay) + &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
    !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
            iL = iaRadLayer(iLay)
            rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            rMPTemp = raVT1(iL)
        !          print *,iLay,iL,raLayAngles(MP2Lay(iL)),rMPTemp
        ! xxxxxxxx  CALL wrtout(iIOUN,caOutName,raFreq,raInten2)

        !c need to do radiative transfer thru this layer
            DO iFr=1,kMaxPts
                raInten(iFr) = raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO

            iPrintX = Inset(iL,iaOp,iNp)
            IF (iPrintX > 0) THEN
                write(kStdWarn,*) 'TOA output rads at',iLay,' th rad layer  = ',iL,' kcarta layer to ind ',iPrintX
                DO iFr=1,kMaxPts
                    raaInten(iFr,iPrintX) = raInten(iFr)
                END DO
            END IF

        END DO
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_GRAY

!************************************************************************
END MODULE scatter_graycld_code
