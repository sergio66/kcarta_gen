! Copyright 2006
! University of Maryland Baltimore County
! All Rights Reserved

MODULE clear_scatter_basic

USE basic_common
USE spline_and_sort_and_common
USE s_writefile
USE s_misc
!USE rad_misc
!USE kbloat
!USE knonlte

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine computes the UPWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile, and NO scattering
    SUBROUTINE RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)             !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: temp(maxnz)                  !temperature profile (1+kProfLayer)
    REAL :: rCos                         !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| > 1.0)
    INTEGER :: iVary                     !should we model temp dependance???
! 1 yes EXP, +2 yes LINEAR, -1 no    !!! ORIGINALLY 0 = linear so sheck this
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at bottom of layer
! utput : intensity at top of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaP1,iVaryVary
    REAL :: rBeta,rTT,rZeta,rBooga,radtot,rad1
    INTEGER :: iRTSPEC
    REAL :: planck1,planck0,del,gwak,tau0,trans

    iVaryVary = iVary

    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileUPWELL, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer+1
    END IF

    IF ((iBeta >= kProfLayer-15) .AND. (iVaryVary >= 2)) THEN
    !!!! if we use RTSPEC, we get junky results close to TOA because of
    !!!! real vs double precision
        iVaryVary = -1
    !!!! but i've changed this so it mimics GASRT2 tau->0 approx
        iVaryVary = iVary
    END IF

!!! model the variation as B(x) = Bo exp(rBooga x)
!!! recall B(bottom) = Bb = Bo exp(rBooga tau) = Bo (since tau = 0)  ==> Bo = Bb
!!! recall B(top)    = Bt = Bo exp(rBooga tau)
!!! this rBooga = 1/tau ln(Bt/Bb) which varies with wavenumber as tau varies with wavenumber
!!!
!!!!this is how temperature in layer varies with tau
    IF (iVaryVary == +1) THEN        !!!!exponential in tau dependance of T
        rBooga = log(TEMP(iBeta+1)/TEMP(iBeta))
    ELSEIF (iVaryVary >= 2) THEN       !!!!linear in tau dependance of T
        rBooga = 0.0
    ELSEIF (iVaryVary == -1) THEN       !!!!no tau dependance of T
        rBooga = 0.0
    END IF

    IF (iVaryVary >= 2) THEN
        iRTSPEC = 1
    ELSE
        iRTSPEC = -1    !!!RTSPEC does a simple "exponential in rad" way
    END IF

    IF (iVary == -2) THEN
    !!!NO LAYER EMISSION csun debug!!!
        DO iFr=1,kMaxPts
            raInten(iFr) = raInten(iFr)*exp(-raaAbs(iFr,iL)/rCos)
        END DO

    ELSEIF ((iVary > -1) .AND. (iRTSPEC < 0)) THEN
    !!!either exp temperature dependance or none; rBooga carries this info
    !!! >>>>>>>>>>>>>>> this is basically exp in tau <<<<<<<<<<<<<<<<<<<<<<
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rbeta = 1/raaAbs(iFr,iL) * rBooga
                rTT   = ttorad(raFreq(iFr),TEMP(iBeta))/(1 + rbeta*rCos)
                rZeta = (raInten(iFr) - rTT) * exp(-raaAbs(iFr,iL)/rCos)
                raInten(iFr) = rZeta + rTT * exp(raaAbs(iFr,iL) * rbeta)
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rbeta = 1/(raaAbs(iFr,iL)*rFrac) * rBooga
                rTT   = ttorad(raFreq(iFr),TEMP(iBeta))/(1 + rbeta*rCos)
                rZeta = (raInten(iFr) - rTT) * exp(-raaAbs(iFr,iL)*rFrac/rCos)
                raInten(iFr) = rZeta + rTT * exp(raaAbs(iFr,iL)*rFrac * rbeta)
            END DO
        END IF

    ELSEIF ((iVary > -1) .AND. (iRTSPEC >= 0)) THEN
    !!!!do the RTSPEC way  .... see GASRT2 in RTSPEC
        IF (rFrac >= 0.9999) THEN !!!full layer
            gwak = 1.0
            iBetaP1 = iBeta + 1
        ELSE IF (rFrac < 0.9999) THEN !!!partial layer
            gwak = rFrac
            IF ((TEMP(iBeta+1) < 150) .OR. (TEMP(iBeta+1) > 350)) THEN
                iBetaP1 = ibeta
            ELSE
                iBetaP1 = ibeta + 1
            END IF
        END IF
        IF (rFrac >= 1.0000) gwak = 1.0
        IF (rFrac < 0.9999) gwak = rFrac
        rad1=raInten(1)
        DO iFr=1,kMaxPts
            planck1 = ttorad(raFreq(iFr),TEMP(iBeta))
            planck0 = ttorad(raFreq(iFr),TEMP(iBetaP1))
            tau0 = (raaAbs(iFr,iL)*gwak)/rCos
            IF (tau0 < 0.001) THEN
                raInten(iFr) = raInten(iFr)*(1-tau0) + tau0*0.5*(PLANCK0+PLANCK1)
            ELSE
                del = (planck1-planck0)/tau0
                trans = exp(-tau0)
                raInten(iFr) = raInten(iFr)*trans + (planck0+del &
                - trans*(planck0+del*(1.0+tau0)))
            END IF
        END DO

    END IF

    RETURN
    end SUBROUTINE RT_ProfileUPWELL

!************************************************************************
! this subroutine computes the DNWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile and NO scattering
    SUBROUTINE RT_ProfileDNWELL(raFreq,raaAbs,iL,TEMP,rCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)             !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: temp(maxnz)                  !temperature profile (1+kProfLayer)
    REAL :: rCos                         !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| = 1.0)
    INTEGER :: iVary                     !should we model temp dependance???
! 1 yes EXP, 2 yes LINEAR, -1 no     !!! originally 0 = linear so sheck this
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at top of layer
! utput : intensity at bottom of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaM1
    REAL :: rBeta,rTT,rZeta,rBooga,mu
    INTEGER :: iRTSPEC
    REAL :: planck1,planck0,del,gwak,tau0,trans

    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileDNWELL, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer+1
    END IF

!!!!this is how temperature in layer varies with tau
    IF (iVary == +1) THEN          !!!!exponential in tau dependance of T
        rBooga = log(TEMP(iBeta+1)/TEMP(iBeta))
    ELSEIF (iVary >= 2) THEN       !!!!linear in tau dependance of T
        rBooga = 0.0
    ELSEIF (iVary == -1) THEN       !!!!no tau dependance of T
        rBooga = 0.0
    END IF

    IF (iVary >= 2) THEN
        iRTSPEC = 1             !!!RTSPEC does a simple "linear in tau" way
    ELSE
        iRTSPEC = -1
    END IF

    mu = abs(rCos)

    IF (iVary == -2) THEN
    !!!NO LAYER EMISSION csun debug!!!
        DO iFr=1,kMaxPts
            raInten(iFr) = raInten(iFr)*exp(-raaAbs(iFr,iL)/rCos)
        END DO
    ELSEIF ((iVary >= -1) .AND. (iRTSPEC < 0)) THEN
    !!!either exp temperature dependace or none; rBooga carries this info
        IF (rFrac >= 0.9999) THEN
            DO iFr=1,kMaxPts
                rbeta = 1/raaAbs(iFr,iL) * rBooga
                rTT   = ttorad(raFreq(iFr),TEMP(iBeta))/(rbeta*mu - 1)
                rZeta = exp(-raaAbs(iFr,iL)/mu) * exp(rBeta*raaAbs(iFr,iL)) - 1.0
                raInten(iFr) = raInten(iFr)* exp(-raaAbs(iFr,iL)/mu) + rTT*rZeta
            END DO
        ELSE
            DO iFr=1,kMaxPts
                rbeta = 1/(raaAbs(iFr,iL)*rFrac) * rBooga
                rTT   = ttorad(raFreq(iFr),TEMP(iBeta))/(rbeta*mu - 1)
                rZeta = &
                exp(-raaAbs(iFr,iL)*rFrac/mu)*exp(rBeta*raaAbs(iFr,iL)*rFrac)-1.0
                raInten(iFr) = raInten(iFr)*exp(-raaAbs(iFr,iL)*rFrac/mu)+rTT*rZeta
            END DO
        END IF

    ELSEIF ((iVary >= -1) .AND. (iRTSPEC >= 0)) THEN
    !!!!do the RTSPEC way  .... see GASRT2 in RTSPEC
        IF (rFrac >= 0.9999) THEN !!!full layer
            gwak = 1.0
            iBetaM1 = iBeta - 1
        ELSE IF (rFrac > 0.0) THEN !!!partial layer
            gwak = rFrac
            IF ((TEMP(iBeta-1) < 150) .OR. (TEMP(iBeta-1) > 350)) THEN
                iBetaM1 = ibeta
            ELSE
                iBetaM1 = ibeta - 1
            END IF
        END IF
        DO iFr=1,kMaxPts
            planck0 = ttorad(raFreq(iFr),TEMP(iBeta))
            planck1 = ttorad(raFreq(iFr),TEMP(iBetaM1))
            tau0 = (raaAbs(iFr,iL)*gwak)/rCos
            IF (tau0 < 0.001) THEN
                raInten(iFr) = raInten(iFr)*(1-tau0) + tau0*0.5*(PLANCK0+PLANCK1)
            ELSE
                del = (planck1-planck0)/tau0
                trans = exp(-tau0)
                raInten(iFr) = raInten(iFr)*trans + (PLANCK1-DEL &
                - TRANS*(PLANCK1-DEL*(1.0+tau0)))
            END IF
        END DO
    END IF


    RETURN
    end SUBROUTINE RT_ProfileDNWELL

!************************************************************************
! this subroutine computes the UPWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile, and NO scattering
!!! ref : IEEE TRANSACTIONS ON GEO AND REMOTE SENSING,
!!!   VOL. 44, NO. 5, MAY 2006, Forward Model and Jacobians for Tropospheric
!!!   Emission Spectrometer Retrievals
!!!   Shepard A. Clough, Mark W. Shephard, John Worden, Patrick D. Brown,
!!!   Helen M. Worden, Mingzhao Luo, Clive D. Rodgers, Curtis P. Rinsland,
!!!   Aaron Goldman, Linda Brown, Susan S. Kulawik, Annmarie Eldering, Michael
!!!   Lampel, Greg Osterman, Reinhard Beer, Kevin Bowman, Karen E. Cady-Pereira,
!!!   and Eli J. Mlawer

!!! ref : MODTRAN Cloud and Multiple Scattering Upgrades with Application to AVIRIS
!!!   A. Berk,* L. S. Bernstein,* G. P. Anderson,† P. K. Acharya,*
!!!   D. C. Robertson,* J. H. Chetwynd,† and S. M. Adler-Golden*
!!!   REMOTE SENS. ENVIRON. 65:367–375 (1998)
!!!   Elsevier Science Inc., 1998 0034-4257/98/$19.00
!!!   655 Avenue of the Americas, New York, NY 10010

!!! or do simplie linear in tau YAY

    SUBROUTINE RT_ProfileUPWELL_LINEAR_IN_TAU( &
    raFreq,raaAbs,iL,TEMPLEV,TEMPLAY,rCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)              !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: tempLEV(maxnz)               !level temperature profile (1+kProfLayer)
    REAL :: tempLAY(kMixFilRows)         !layer temperature profile (0+kProfLayer)
    REAL :: rCos                         !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| > 1.0)
    INTEGER :: iVary                     !should we model temp dependance??? +2,+3,+4
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at bottom of layer
! utput : intensity at top of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaP1
    REAL :: rBeff,rFcn
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: rZeta,rZeta2,rAbs,rTrans

    IF (iVary < 2) THEN
        write(kStdErr,*) 'this is upwell for linear in tau .. need iVary = 2 or 3 or 4'
        CALL DoStop
    END IF

    IF (iVary == 41) iVary = 43     !!! have debugged 04, 42, 43 for small tau O(tau^2)
          
    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileUPWELL_LINTAU, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer
    END IF

    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level
    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level  XXXXX this is the one we want XXXXX
    CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
!! (which is NOT necessarily average of above 2)
    IF (kOuterLoop == 1) THEN
        write(kStdWarn,2345) iL,TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    END IF
          
    1234 FORMAT(I3,3(' ',F10.3))
    2345 FORMAT('up [iLDN=iL iLay iLUP=iLp1]',I3,3(' ',F10.3))
     
!      IF (iVary .EQ. 4) THEN
!        ! new option
!        DO iFr = 1,kMaxPts
!          raIntenAvg(iFr) = 0.5 * (raIntenP(iFr) + raIntenP1(iFr))
!        END DO
!      END IF

    IF (iVary == 2) THEN
    !!! lim tau --> 0 , rFcn --> 0
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP0)
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        END IF

    ELSEIF (iVary == +3) THEN
    !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        END IF

    ELSEIF (iVary == +40) THEN
    !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
    !        print *,iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP1(iFr)*(1-rTrans) + (raIntenP(iFr) - raIntenP1(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + rZeta
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP1(iFr)*(1-rTrans) + (raIntenP(iFr) - raIntenP1(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr)*exp(-rAbs/rCos) + rZeta
            END DO
        END IF

    ELSEIF (iVary == +41) THEN
    !        print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! PADE APPROX, two term (combo of GENLN2 and LBLRTM)
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rTrans = exp(-rAbs)
            rZeta = 0.2*rAbs    !! pade one
            rFcn = (raIntenAvg(iFr) + rZeta*raIntenP1(iFr))/(1+rZeta)
            rZeta = 0.193*rAbs    !! pade two
            rZeta2 = 0.013*rAbs*rAbs    !! pade two
            rFcn = (raIntenAvg(iFr) + (rZeta + rZeta2)*raIntenP1(iFr))/(1+rZeta+rZeta2)
            rFcn = (1-rTrans)*rFcn
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +42) THEN
    !        print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, GENLN2 style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP1(iFr))
            IF (rAbs >= 0.05) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP1(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs
                rFcn = rAbs*raIntenP1(iFr) + rZeta*(1-rAbs/2) - rTrans * rZeta
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'up',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +43) THEN
    !         http://www.wolframalpha.com/input/?i=1-2*%281%2Fx-exp%28-x%29%2F%281-exp%28-x%29%29%29
    !        print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = raIntenP1(iFr) - raIntenAvg(iFr)
            IF (rAbs >= 0.06) THEN
                rTrans = exp(-rAbs)
                rZeta2 = 1.0 - 2.0*(1/rAbs - rTrans/(1-rTrans))
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            ELSE
                rTrans = 1 - rAbs + 0.5*(rAbs * rAbs)
                rZeta2 = rAbs/6.0 - (rAbs**3)/360.0 + (rAbs**5)/15120.0   !! mathematica
                rZeta2 = rAbs/6.0
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'up',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    ELSEIF (iVary == +4) THEN
    !        print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, MY style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP1(iFr))
            IF (rAbs > 0.1) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP1(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs + 0.5*rAbs**2
                rZeta2 = rZeta*(rAbs/2-(rAbs**2)/3+(rAbs**3)/6)
                rFcn   = (1-rTrans)*raIntenP1(iFr) + rZeta2
            END IF
        !          IF (iFr .EQ. 1) THEN
        !            print *,'<<up>>',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    END IF

    RETURN
    end SUBROUTINE RT_ProfileUPWELL_LINEAR_IN_TAU

!************************************************************************
! this subroutine computes the DNWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile, and NO scattering
!!! ref : IEEE TRANSACTIONS ON GEO AND REMOTE SENSING,
!!!   VOL. 44, NO. 5, MAY 2006, Forward Model and Jacobians for Tropospheric
!!!   Emission Spectrometer Retrievals
!!!   Shepard A. Clough, Mark W. Shephard, John Worden, Patrick D. Brown,
!!!   Helen M. Worden, Mingzhao Luo, Clive D. Rodgers, Curtis P. Rinsland,
!!!   Aaron Goldman, Linda Brown, Susan S. Kulawik, Annmarie Eldering, Michael
!!!   Lampel, Greg Osterman, Reinhard Beer, Kevin Bowman, Karen E. Cady-Pereira,
!!!   and Eli J. Mlawer

!!! ref : MODTRAN Cloud and Multiple Scattering Upgrades with Application to AVIRIS
!!!   A. Berk,* L. S. Bernstein,* G. P. Anderson,† P. K. Acharya,*
!!!   D. C. Robertson,* J. H. Chetwynd,† and S. M. Adler-Golden*
!!!   REMOTE SENS. ENVIRON. 65:367–375 (1998)
!!!   Elsevier Science Inc., 1998 0034-4257/98/$19.00
!!!   655 Avenue of the Americas, New York, NY 10010

!!! or do simplie linear in tau YAY

    SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU( &
    raFreq,raaAbs,iL,TEMPLEV,TEMPLAY,rCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)             !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: tempLEV(maxnz)               !level temperature profile (1+kProfLayer)
    REAL :: tempLAY(kMixFilRows)         !layer temperature profile (0+kProfLayer)
    REAL :: rCos                         !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| > 1.0)
    INTEGER :: iVary                     !should we model temp dependance??? +2,+3,+4
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at top of layer
! utput : intensity at bottom of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaP1
    REAL :: rBeff,rFcn
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: rZeta,rZeta2,rAbs,rTrans

    IF (iVary < 2) THEN
        write(kStdErr,*) 'this is downwell for linear in tau .. need iVary = 2 or 3 or 4'
        CALL DoStop
    END IF

    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileDNWELL_LINTAU, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    IF (iVary == 41) iVary = 43     !!! have debugged 04, 42, 43 for small tau O(tau^2)

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer
    END IF

    IF (iVary < 4) THEN
        IF (iBeta > 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta-1),raIntenP1)
        ELSEIF (iBeta == 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP1)
        END IF
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)
    END IF

    IF (iVary >= 4) THEN
    !! new option
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX this is the one we want XXXXX
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
    !! (which is NOT necessarily average of above 2)

        IF (kOuterLoop == 1) THEN
            write(kStdWarn,2345) iL,TEMPLEV(iBeta+1),TEMPLAY(iBeta),TEMPLEV(iBeta)
        END IF
    END IF
          
    1234 FORMAT(I3,3(' ',F10.3))
    2345 FORMAT('dn [iLUP=iLp1 iLay=iL iLDN=iL]',I3,3(' ',F10.3))
     
    IF (iVary == 2) THEN
    !!! lim tau --> 0 , rFcn --> 0
        write(kStdErr,*) 'huh iVary = 2 is a little buggy'
        CALL DoStop
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP0)
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        END IF

    ELSEIF (iVary == +3) THEN
    !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        END IF

    ELSEIF (iVary == +40) THEN
    !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + rZeta
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + rZeta
            END DO
        END IF

    ELSEIF (iVary == +41) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rTrans = exp(-rAbs)
            rZeta = 0.2*rAbs    !! pade one
            rFcn = (raIntenAvg(iFr) + rZeta*raIntenP(iFr))/(1+rZeta)
            rZeta = 0.193*rAbs    !! pade two
            rZeta2 = 0.013*rAbs*rAbs    !! pade two
            rFcn = (raIntenAvg(iFr) + (rZeta + rZeta2)*raIntenP(iFr))/(1+rZeta+rZeta2)
            rFcn = (1-rTrans)*rFcn
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +42) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, GENLN2 style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs >= 0.05) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs
                rFcn = rAbs*raIntenP(iFr) + rZeta*(1-rAbs/2) - rTrans * rZeta
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'down',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +43) THEN
    !        print *,'dn flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = raIntenP(iFr) - raIntenAvg(iFr)
            IF (rAbs >= 0.06) THEN
                rTrans = exp(-rAbs)
                rZeta2 = 1.0 - 2.0*(1/rAbs - rTrans/(1-rTrans))
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            ELSE
                rTrans = 1 - rAbs + 0.5*(rAbs * rAbs)
                rZeta2 = rAbs/6.0 - (rAbs**3)/360.0 + (rAbs**5)/15120.0  !! mathematica
                rZeta2 = rAbs/6.0
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'up',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO
              
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    ELSEIF (iVary == +4) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, MY style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs > 0.1) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs + 0.5*rAbs**2
                rZeta2 = rZeta*(rAbs/2-(rAbs**2)/3+(rAbs**3)/6)
                rFcn   = (1-rTrans)*raIntenP(iFr) + rZeta2
            END IF
        !          IF (iFr .EQ. 1) THEN
        !            print *,'>>down<<',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    END IF
          
    RETURN
    end SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU

!************************************************************************
! this subroutine computes the DNWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile, and NO scattering
! assumes ONE angle for all freq points
!!! ref : IEEE TRANSACTIONS ON GEO AND REMOTE SENSING,
!!!   VOL. 44, NO. 5, MAY 2006, Forward Model and Jacobians for Tropospheric
!!!   Emission Spectrometer Retrievals
!!!   Shepard A. Clough, Mark W. Shephard, John Worden, Patrick D. Brown,
!!!   Helen M. Worden, Mingzhao Luo, Clive D. Rodgers, Curtis P. Rinsland,
!!!   Aaron Goldman, Linda Brown, Susan S. Kulawik, Annmarie Eldering, Michael
!!!   Lampel, Greg Osterman, Reinhard Beer, Kevin Bowman, Karen E. Cady-Pereira,
!!!   and Eli J. Mlawer

!!! ref : MODTRAN Cloud and Multiple Scattering Upgrades with Application to AVIRIS
!!!   A. Berk,* L. S. Bernstein,* G. P. Anderson,† P. K. Acharya,*
!!!   D. C. Robertson,* J. H. Chetwynd,† and S. M. Adler-Golden*
!!!   REMOTE SENS. ENVIRON. 65:367–375 (1998)
!!!   Elsevier Science Inc., 1998 0034-4257/98/$19.00
!!!   655 Avenue of the Americas, New York, NY 10010

!!! or do simplie linear in tau YAY

! this is SAME as RT_ProfileDNWELL_CONST_IN_TAU_FORFLUX EXCEPT very importantly, since the
! atmosphere was defined for downlook instrument, that means we have to be VERY CAREFUL with directions
! so as to ensure
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX this is the one we want XXXXXXXX
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level
! CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
!! (which is NOT necessarily average of above 2)
! is changed to
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX this is the one we want XXXXXXXX
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta-1),raIntenP1)   !! ttorad of upper level
! CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
!! (which is NOT necessarily average of above 2)
                                                                       
    SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX( &
    raFreq,raaAbs,iL,TEMPLEV,TEMPLAY,rCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)             !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: tempLEV(maxnz)               !level temperature profile (1+kProfLayer)
    REAL :: tempLAY(kMixFilRows)         !layer temperature profile (0+kProfLayer)
    REAL :: rCos                         !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| > 1.0)
    INTEGER :: iVary                     !should we model temp dependance??? +2,+3,+4
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at top of layer
! utput : intensity at bottom of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaP1
    REAL :: rBeff,rFcn
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: rZeta,rZeta2,rAbs,rTrans

    IF (iVary < 2) THEN
        write(kStdErr,*) 'this is downwell for linear in tau .. need iVary = 2 or 3 or 4'
        CALL DoStop
    END IF

    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileDNWELL_LINTAU, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    IF (iVary == 41) iVary = 43     !!! have debugged 04, 42, 43 for small tau O(tau^2)

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer
    END IF

    IF (iVary < 4) THEN
        IF (iBeta > 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)
        ELSEIF (iBeta == 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP1)
        END IF
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)
    END IF

! RT_ProfileUPWELL_LINEAR_IN_TAU
!     iBeta = MOD(iL,kProfLayer)
!     IF (iBeta .EQ. 0) THEN
!       iBeta = kProfLayer
!     END IF
!     IF (iL .EQ. kProfLayer+1) THEN
!      iBeta = kProfLayer
!    END IF
!    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level
!    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level  XXXXX this is the one we want XXXXX
!    CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
                                                                 
    IF (iVary >= 4) THEN
    !! new option
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)    !! ttorad of lower level XXXX this is the one we want XXXXXXXX
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)  !! ttorad of upper level
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
    !! (which is NOT necessarily average of above 2)
        IF (kOuterLoop == 1) THEN
            write(kStdWarn,2345) iL,TEMPLEV(iBeta+1),TEMPLAY(iBeta),TEMPLEV(iBeta)
        END IF
    END IF
          
    1234 FORMAT(I3,3(' ',F10.3))
    2345 FORMAT('dn [iLUP=iLp1 iLay=iL iLDN=iL]',I3,3(' ',F10.3))
          
    IF (iVary == 2) THEN
    !!! lim tau --> 0 , rFcn --> 0
        write(kStdErr,*) 'huh iVary = 2 is a little buggy'
        CALL DoStop
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP0)
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/rCos))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*rCos*(rAbs/rCos-1) + &
                rFcn*rCos*exp(-rAbs/rCos)
            END DO
        END IF

    ELSEIF (iVary == +3) THEN
    !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + &
                rFcn * (1 - exp(-rAbs/rCos))
            END DO
        END IF

    ELSEIF (iVary == +40) THEN
    !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + rZeta
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/rCos)
                    rFcn = rCos/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/rCos) + rZeta
            END DO
        END IF

    ELSEIF (iVary == +41) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rTrans = exp(-rAbs)
            rZeta = 0.2*rAbs    !! pade one
            rFcn = (raIntenAvg(iFr) + rZeta*raIntenP(iFr))/(1+rZeta)
            rZeta = 0.193*rAbs    !! pade two
            rZeta2 = 0.013*rAbs*rAbs    !! pade two
            rFcn = (raIntenAvg(iFr) + (rZeta + rZeta2)*raIntenP(iFr))/(1+rZeta+rZeta2)
            rFcn = (1-rTrans)*rFcn
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +42) THEN
    !        print *,'fluxybuyxy down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta-1)
    !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, GENLN2 style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs >= 0.05) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs
                rFcn = rAbs*raIntenP(iFr) + rZeta*(1-rAbs/2) - rTrans * rZeta
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'down',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +43) THEN
    !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = raIntenP(iFr) - raIntenAvg(iFr)
            IF (rAbs >= 0.06) THEN
                rTrans = exp(-rAbs)
                rZeta2 = 1.0 - 2.0*(1/rAbs - rTrans/(1-rTrans))
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            ELSE
                rTrans = 1 - rAbs + 0.5*(rAbs * rAbs)
                rZeta2 = rAbs/6.0 - (rAbs**3)/360.0 + (rAbs**5)/15120.0  !! mathematica
                rZeta2 = rAbs/6.0
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            !          print *,rAbs,rTrans,(1-rTrans),raIntenAvg(iFr),rZeta,rZeta2,rFcn,rCos,rFrac
            !          call dostop
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'up',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO
    !        print *,'dn flux ',iL,iBeta,rFrac,raaAbs(1,iL),rAbs,rTrans,TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1),rFcn,raInten(1)

    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    ELSEIF (iVary == +4) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, MY style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/rCos*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs > 0.1) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs + 0.5*rAbs**2
                rZeta2 = rZeta*(rAbs/2-(rAbs**2)/3+(rAbs**3)/6)
                rFcn   = (1-rTrans)*raIntenP(iFr) + rZeta2
            END IF
        !          IF (iFr .EQ. 1) THEN
        !            print *,'>>down<<',iL,iBeta,rCos,rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    END IF
          
    RETURN
    end SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX

!************************************************************************
! this subroutine computes the DNWARD rad transfer thru an atmospheric layer,
! assuming there is a temperature profile, and NO scattering
! assumes DIFFERENT angle for all freq points
!!! ref : IEEE TRANSACTIONS ON GEO AND REMOTE SENSING,
!!!   VOL. 44, NO. 5, MAY 2006, Forward Model and Jacobians for Tropospheric
!!!   Emission Spectrometer Retrievals
!!!   Shepard A. Clough, Mark W. Shephard, John Worden, Patrick D. Brown,
!!!   Helen M. Worden, Mingzhao Luo, Clive D. Rodgers, Curtis P. Rinsland,
!!!   Aaron Goldman, Linda Brown, Susan S. Kulawik, Annmarie Eldering, Michael
!!!   Lampel, Greg Osterman, Reinhard Beer, Kevin Bowman, Karen E. Cady-Pereira,
!!!   and Eli J. Mlawer

!!! ref : MODTRAN Cloud and Multiple Scattering Upgrades with Application to AVIRIS
!!!   A. Berk,* L. S. Bernstein,* G. P. Anderson,† P. K. Acharya,*
!!!   D. C. Robertson,* J. H. Chetwynd,† and S. M. Adler-Golden*
!!!   REMOTE SENS. ENVIRON. 65:367–375 (1998)
!!!   Elsevier Science Inc., 1998 0034-4257/98/$19.00
!!!   655 Avenue of the Americas, New York, NY 10010

!!! or do simplie linear in tau YAY

! this is SAME as RT_ProfileDNWELL_CONST_IN_TAU_FORFLUX EXCEPT very importantly, since the
! atmosphere was defined for downlook instrument, that means we have to be VERY CAREFUL with directions
! so as to ensure
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX this is the one we want XXXXXXXX
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level
! CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
!! (which is NOT necessarily average of above 2)
! is changed to
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX this is the one we want XXXXXXXX
! CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta-1),raIntenP1)   !! ttorad of upper level
! CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
!! (which is NOT necessarily average of above 2)
                                                                       
    SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang( &
    raFreq,raaAbs,iL,TEMPLEV,TEMPLAY,raCos,rFrac,iVary,raInten)
          
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: raFreq(kMaxPts)             !wavenumbers
    REAL :: raaAbs(kMaxPts,kMixFilRows)  !mixing table
    INTEGER :: iL                        !which row of mix table
    REAL :: tempLEV(maxnz)               !level temperature profile (1+kProfLayer)
    REAL :: tempLAY(kMixFilRows)         !layer temperature profile (0+kProfLayer)
    REAL :: raCos(kMaxPts)               !satellite view angle
    REAL :: rFrac                        !fractional (0<f<1) or full (|f| > 1.0)
    INTEGER :: iVary                     !should we model temp dependance??? +2,+3,+4
! output parameters
    REAL :: raInten(kMaxPts)             !input  : intensity at top of layer
! utput : intensity at bottom of layer

! local variables
    INTEGER :: iFr,iBeta,iBetaP1
    REAL :: rBeff,rFcn
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: rZeta,rZeta2,rAbs,rTrans

    IF (iVary < 2) THEN
        write(kStdErr,*) 'this is downwell for linear in tau .. need iVary = 2 or 3 or 4'
        CALL DoStop
    END IF

    IF (rFrac < 0) THEN
        write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileDNWELL_LINTAU, reset to > 0'
        rFrac = abs(rFrac)
    END IF

    IF (iVary == 41) iVary = 43     !!! have debugged 04, 42, 43 for small tau O(tau^2)

    iBeta = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF

    IF (iL == kProfLayer+1) THEN
        iBeta = kProfLayer
    END IF

    IF (iVary < 4) THEN
        IF (iBeta > 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)
        ELSEIF (iBeta == 1) THEN
            CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP1)
        END IF
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)
    END IF

! RT_ProfileUPWELL_LINEAR_IN_TAU
!     iBeta = MOD(iL,kProfLayer)
!     IF (iBeta .EQ. 0) THEN
!       iBeta = kProfLayer
!     END IF
!     IF (iL .EQ. kProfLayer+1) THEN
!      iBeta = kProfLayer
!    END IF
!    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level
!    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level  XXXXX this is the one we want XXXXX
!    CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
                                                                 
    IF (iVary >= 4) THEN
    !! new option
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)    !! ttorad of lower level XXXX this is the one we want XXXXXXXX
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)  !! ttorad of upper level
        CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
    !! (which is NOT necessarily average of above 2)
        IF (kOuterLoop == 1) THEN
            write(kStdWarn,2345) iL,TEMPLEV(iBeta+1),TEMPLAY(iBeta),TEMPLEV(iBeta)
        END IF
    END IF
          
    1234 FORMAT(I3,3(' ',F10.3))
    2345 FORMAT('dn [iLUP=iLp1 iLay=iL iLDN=iL]',I3,3(' ',F10.3))
     
    IF (iVary == 2) THEN
    !!! lim tau --> 0 , rFcn --> 0
        write(kStdErr,*) 'huh iVary = 2 is a little buggy'
        CALL DoStop
        CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP0)
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/raCos(iFr)))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*raCos(iFr)*(rAbs/raCos(iFr)-1) + &
                rFcn*raCos(iFr)*exp(-rAbs/raCos(iFr))
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = (raIntenP1(iFr) - raIntenP0(iFr) + 1.0e-10)/(rAbs + 1.0e-10)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + &
                raIntenP0(iFr) * (1 - exp(-rAbs/raCos(iFr)))
                IF (rAbs >= 0.001) &
                raInten(iFr) = raInten(iFr) + rFcn*raCos(iFr)*(rAbs/raCos(iFr)-1) + &
                rFcn*raCos(iFr)*exp(-rAbs/raCos(iFr))
            END DO
        END IF

    ELSEIF (iVary == +3) THEN
    !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/raCos(iFr))
                    rFcn = raCos(iFr)/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + &
                rFcn * (1 - exp(-rAbs/raCos(iFr)))
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                rFcn = 1.0
                IF (rAbs >= 0.001) THEN
                    rFcn = exp(-rAbs/raCos(iFr))
                    rFcn = raCos(iFr)/rAbs - rFcn/(1-rFcn)
                END IF
                rFcn = raIntenP1(iFr) + 2*(raIntenAvg(iFr)-raIntenP1(iFr))*rFcn
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + &
                rFcn * (1 - exp(-rAbs/raCos(iFr)))
            END DO
        END IF

    ELSEIF (iVary == +40) THEN
    !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
    !!! lim tau --> 0 , rFcn --> 1
        IF (rFrac >= 0.9999) THEN
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/raCos(iFr))
                    rFcn = raCos(iFr)/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + rZeta
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                rAbs = raaAbs(iFr,iL)*rFrac
                IF (rAbs >= 0.0001) THEN
                    rTrans = exp(-rAbs/raCos(iFr))
                    rFcn = raCos(iFr)/rAbs * (1 - rTrans)
                ELSE
                    rFcn = 1.0
                    rTrans = 1.0
                END IF
                rZeta = raIntenP(iFr)*(1-rTrans) + (raIntenP1(iFr) - raIntenP(iFr))*(rFcn - rTrans)
                raInten(iFr) = raInten(iFr) * exp(-rAbs/raCos(iFr)) + rZeta
            END DO
        END IF

    ELSEIF (iVary == +41) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/raCos(iFr)*rFrac
            rTrans = exp(-rAbs)
            rZeta = 0.2*rAbs    !! pade one
            rFcn = (raIntenAvg(iFr) + rZeta*raIntenP(iFr))/(1+rZeta)
            rZeta = 0.193*rAbs    !! pade two
            rZeta2 = 0.013*rAbs*rAbs    !! pade two
            rFcn = (raIntenAvg(iFr) + (rZeta + rZeta2)*raIntenP(iFr))/(1+rZeta+rZeta2)
            rFcn = (1-rTrans)*rFcn
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +42) THEN
    !        print *,'fluxybuyxy down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta-1)
    !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, GENLN2 style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/raCos(iFr)*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs >= 0.05) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs
                rFcn = rAbs*raIntenP(iFr) + rZeta*(1-rAbs/2) - rTrans * rZeta
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'down',iL,iBeta,raCos(iFr),rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    ELSEIF (iVary == +43) THEN
    !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/raCos(iFr)*rFrac
            rZeta = raIntenP(iFr) - raIntenAvg(iFr)
            IF (rAbs >= 0.06) THEN
                rTrans = exp(-rAbs)
                rZeta2 = 1.0 - 2.0*(1/rAbs - rTrans/(1-rTrans))
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            ELSE
                rTrans = 1 - rAbs + 0.5*(rAbs * rAbs)
                rZeta2 = rAbs/6.0 - (rAbs**3)/360.0 + (rAbs**5)/15120.0  !! mathematica
                rZeta2 = rAbs/6.0
                rFcn = (1-rTrans)*(raIntenAvg(iFr) + rZeta * rZeta2)
            !          print *,rAbs,rTrans,(1-rTrans),raIntenAvg(iFr),rZeta,rZeta2,rFcn,raCos(iFr),rFrac
            !          call dostop
            END IF
        !          if (iFr .EQ. 1) THEN
        !            print *,'up',iL,iBeta,raCos(iFr),rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO
    !        print *,'dn flux ',iL,iBeta,rFrac,raaAbs(1,iL),rAbs,rTrans,TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1),rFcn,raInten(1)

    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    ELSEIF (iVary == +4) THEN
    !        print *,'down flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
    !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
    !!! LINEAR IN TAU, MY style
        DO iFr = 1,kMaxPts
            rAbs = raaAbs(iFr,iL)/raCos(iFr)*rFrac
            rZeta = 2*(raIntenAvg(iFr)-raIntenP(iFr))
            IF (rAbs > 0.1) THEN
                rTrans = exp(-rAbs)
                rFcn = (1-rTrans)*(raIntenP(iFr) + rZeta/rAbs) - rTrans * rZeta
            ELSE
                rTrans = 1 - rAbs + 0.5*rAbs**2
                rZeta2 = rZeta*(rAbs/2-(rAbs**2)/3+(rAbs**3)/6)
                rFcn   = (1-rTrans)*raIntenP(iFr) + rZeta2
            END IF
        !          IF (iFr .EQ. 1) THEN
        !            print *,'>>down<<',iL,iBeta,raCos(iFr),rAbs,rTrans,rZeta,rFcn,raInten(iFr)
        !          end if
            raInten(iFr) = raInten(iFr)*rTrans + rFcn
        END DO

    END IF
          
    RETURN
    end SUBROUTINE RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang

!************************************************************************
! this function does a temperature interpolation on a fractional layer
! this uses modified Scott Hannon's method of doing a quad fit to the layer,
! layer above, layer below  of the form
!     T = a (ln P(avg))^2 + b (ln P(avg)) + c

! this is almost the same as REAL FUNCTION InterpTemp, but it tries to
! account for large temperature difference between surface and air temp
! of the lowest layer, by using the surface parameters
    REAL FUNCTION InterpTempSurf(iProfileLayers,raPressLevels,raVTemp,rFrac, &
    iTopORBot,iL,rSurfTemp,rSurfPress)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raVTemp  = array containing the original 1.0 fraction temps
! rFrac    = frac of layer that we need
! iTopORBot= do we need top or bottom of layer (+1/-1)
! iL       = which of the mixed paths

! for a down looking instrument, we need bottom frac
! for a   up looking instrument, we need top frac
! for bottommost layer, we need top frac
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raVTemp(kMixFilRows),rFrac
    INTEGER :: iTopORBot,iL,iProfileLayers
    REAL :: rSurfPress,rSurfTemp

    REAL :: rT,rP         !user specfd pressure, temp calculated at this press
    REAL :: rPavg         !given rP,rP1, need to compute rPavg
    REAL :: rT0,rTm1,rTp1 !avg temps of 3 adjacent layers
    REAL :: rP0,rPm1,rPp1 !avg pressures of 3 adjacent layers
    REAL :: rA,rB,rC      !need to find eqn of quadratic
    REAL :: rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr  !temporary variables
    REAL :: xa(3),ya(3)
    INTEGER :: i0,im1,ip1,iW, i00,im11,ip11
    INTEGER :: iLowest

    iLowest = kProfLayer - iProfileLayers + 1

    iW = iCeil(iL*1.0/(kProfLayer*1.0))  !from which set of mxd paths this is
    i0 = MP2Lay(iL) !lower pressure level .. rP is within this press layer
    ip1 = i0+1      !upper pressure leve1 .. this is one press layer above
    im1 = i0-1      !                     .. this is one press layer below

! have to recompute what the user specified pressure was!!
    IF (iTopORBot == 1) THEN          !top frac of layer
    ! ressure specified by user
        rP=raPressLevels(ip1)+rFrac*(raPressLevels(i0)-raPressLevels(ip1))
    ELSE                                !bot frac of layer
    ! ressure specified by user
        rP=-rFrac*(raPressLevels(i0)-raPressLevels(ip1))+raPressLevels(i0)
    END IF

! compute the average pressure of the fractional layer
    IF (iTopOrBot == 1) THEN
        IF (abs(rP-raPressLevels(ip1)) >= delta) THEN
            rPavg=(rP-raPressLevels(ip1))/alog(rP/raPressLevels(ip1))
        ELSE
            rPavg=rP
        END IF
    ELSE
        IF (abs(rP-raPressLevels(i0)) >= delta) THEN
            rPavg=(raPressLevels(i0)-rP)/alog(raPressLevels(i0)/rP)
        ELSE
            rPavg=rP
        END IF
    END IF

! avg press,temperature of layer i0
    rP0=(raPressLevels(i0)-raPressLevels(ip1))/ &
    alog(raPressLevels(i0)/raPressLevels(ip1))
    rT0=raVTemp(i0+(iW-1)*kProfLayer)
! avg press, temperature of layer i0+1
    rPp1=(raPressLevels(ip1)-raPressLevels(ip1+1))/ &
    alog(raPressLevels(ip1)/raPressLevels(ip1+1))
    rTp1=raVTemp(ip1+(iW-1)*kProfLayer)
! surface parameters
    rPm1 = rSurfPress
    rTm1 = rSurfTemp

! now compute the fit for rT(n) = ax(n)^2 + bx(n) + c where x(n) = alog(P(n))
    rPavg = alog(rPavg)

    rP0=alog(rP0)
    rPp1=alog(rPp1)
    rPm1=alog(rPm1)
           
!      rDp1=rTp1-rT0
!      rDm1=rTm1-rT0

!      rp1=rPp1-rP0
!      rp1sqr=(rPp1-rP0)*(rPp1+rP0)
!      rm1=rPm1-rP0
!      rm1sqr=(rPm1-rP0)*(rPm1+rP0)

!      rA=(rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1)
!      rB=rDp1/rp1-rA*(rp1sqr/rp1)
!      rC=rT0-rA*rP0*rP0-rB*rP0

! finally compute rT
!      rT=rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
!      print *,'rPavg,rT = ',rPavg,rT

! use rSpl
    xa(1) = rPp1
    xa(2) = rP0
    xa(3) = rPm1
    ya(1) = rTp1
    ya(2) = rT0
    ya(3) = rTm1
     
    CALL rspl1(xa,ya,3,rPavg,rT,1)
!      print *,'rPavg,rT = ',exp(rPavg),rT

    InterpTempSurf=rT
    RETURN
    end FUNCTION InterpTempSurf

!************************************************************************
! this is quick clear sky downlook radT, based on rad_main.f : rad_trans_SAT_LOOK_DOWN
    SUBROUTINE quick_clear_radtrans_downlook( &
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP,rCO2MixRatio, &
    raaRadsX,iNumOutX,iWriteToOutputFile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
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
    INTEGER :: iWriteToOutputFile
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP,iProfileLayers
    REAL :: raThickness(kProfLayer),pProf(kProfLayer),rCO2MixRatio, &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

! raaRadsX,iNumOutX are to keep up with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: rCos,raInten2(kMaxPts),rMPTemp
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: rDum1,rDum2
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

! for the NLTE which is not used in this routine
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),rT
    REAL :: raVT1(kMixFilRows)
    REAL :: bt2rad,t2s,rPlanck
    INTEGER :: iFr1,troplayer
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular,iFrX

! for NLTE
    REAL :: suncos,scos1,vsec1
          
    iNumOutX = 0

    rThermalRefl=1.0/kPi
          
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
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

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
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh=-1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh = iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh
          
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
    write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

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

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    !         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF ((iDp > 0) .AND. (iWriteToOutputFile > 0)) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaRadsX)   !!don't worry about raaRadsX here
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                iNumOutX = iNumOutX + 1
                DO iFrX = 1,kMaxPts
                    raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                END DO
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            rT = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raInten(iFr) = rPlanck*(1-rT) + raInten(iFr)*rT
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF ((iDp > 0) .AND. (iWriteToOutputFile > 0)) THEN
            write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaRadsX)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                iNumOutX = iNumOutX + 1
                DO iFrX = 1,kMaxPts
                    raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                END DO
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            rT = exp(-raaAbs(iFr,iL)/rCos)
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raInten(iFr) = rPlanck*(1-rT) + raInten(iFr)*rT
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
    !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
            iL = iaRadLayer(iLay)
            rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            rMPTemp = raVT1(iL)

            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)

            IF (iDoSolar < 0) THEN
                IF ((iDp > 0) .AND. (iWriteToOutputFile > 0)) THEN
                    write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                    DO iFr=1,iDp
                        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                        raSun,-1,iNumLayer,rFracTop,rFracBot, &
                        iProfileLayers,raPressLevels, &
                        iNLTEStart,raaRadsX)
                        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                        iNumOutX = iNumOutX + 1
                        DO iFrX = 1,kMaxPts
                            raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                        END DO
                    END DO
                END IF
            ELSE
                IF (iDp == 1) THEN
                    write(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer'

                    suncos = raSunAngles(iaRadLayer(1))           !! at surface
                    scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
                    vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

                    suncos = cos(suncos*kPi/180.0)
                    scos1  = cos(scos1*kPi/180.0)
                    vsec1  = 1/cos(vsec1*kPi/180.0)

                    DO iFr=1,kMaxPts
                        rT = exp(-raaAbs(iFr,iL)/rCos)
                        rPlanck = ttorad(raFreq(iFr),rMPTemp)
                        raInten2(iFr) = rPlanck*(1-rT) + raInten(iFr)*rT
                    END DO
                                
                    CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
                    iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
                    iNumOutX = iNumOutX + 1
                    DO iFrX = 1,kMaxPts
                        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                    END DO
                !              print *,'abcde',raSunAngles(iaRadLayer(1)),suncos,raFreq(1),raInten(1),raInten2(1)
                    IF ((iDp == 1) .AND. (iWriteToOutputFile > 0)) THEN
                        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                    END IF
                                
                ELSEIF (iDp > 1) THEN
                    write(kStdErr,*) 'oops in scatter_pclsam_cpde, at NLTE, dump more than 1 rad at TOA???'
                    CALL DoStop
                END IF
            END IF            !! if iDoSolar
        END DO              !! do iLay = iHigh,iHigh
    END IF                !! if iHigh > 0

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE quick_clear_radtrans_downlook

!************************************************************************
! this is quick clear sky uplook radT, based on rad_main.f : rad_trans_SAT_LOOK_DOWN
    SUBROUTINE quick_clear_radtrans_uplook( &
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    raaRadsX,iNumOutX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP,iProfileLayers
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

! raaRadsX,iNumOutX are to keep up with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: rCos,raInten2(kMaxPts),rMPTemp
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: rDum1,rDum2
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

! for the NLTE which is not used in this routine
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),rT
    REAL :: raVT1(kMixFilRows)
    REAL :: bt2rad,t2s,rPlanck
    INTEGER :: iFr1,troplayer
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular,iFrX

    iNumOutX = 0

    rThermalRefl=1.0/kPi
          
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
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

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
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the lowest layer that we need to output radiances for
    iHigh = +100000000
    DO iLay=1,iNp
        IF (iaOp(iLay) < iHigh) THEN
            iHigh = iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh
          
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

    iLay = 1
    iL = iaRadLayer(iLay)
    rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = sngl(kTSpace)
    DO iFr=1,kMaxPts
        rT = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
        raInten(iFr) = ttorad(raFreq(iFr),rMPTemp)
    END DO

! now we can compute the downwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the top layer (could be fractional)
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    !         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaRadsX)   !!don't worry about raaRadsX here
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                iNumOutX = iNumOutX + 1
                DO iFrX = 1,kMaxPts
                    raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                END DO
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            rT = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raInten(iFr) = rPlanck*(1-rT) + raInten(iFr)*rT
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then go down thru the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaRadsX)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                iNumOutX = iNumOutX + 1
                DO iFrX = 1,kMaxPts
                    raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                END DO
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            rT = exp(-raaAbs(iFr,iL)/rCos)
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raInten(iFr) = rPlanck*(1-rT) + raInten(iFr)*rT
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the bottommost layer (could be fractional)
    777 CONTINUE
    IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
    !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
            iL = iaRadLayer(iLay)
            rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            rMPTemp = raVT1(iL)

            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaRadsX)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                    iNumOutX = iNumOutX + 1
                    DO iFrX = 1,kMaxPts
                        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
                    END DO
                END DO
            END IF
        END DO
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE quick_clear_radtrans_uplook

!************************************************************************
! this function loads in the specular reflection file
    SUBROUTINE loadspecular(raFreq,raSpecularRefl)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input var
    REAL :: raFreq(kMaxPts)
! output var
    REAL :: raSpecularRefl(kMaxPts)

! local vars
    INTEGER :: iIOUN,iL,iFlip,iI
    CHARACTER(80) :: fname,caLine
    REAL :: raW(kMaxPts),raR(kMaxPts),r1,r2,r3,r4,r5,raTemp(kMaxPts)

! name = sun view d(phi) windspeed
    fname = '/home/sergio/SBDART/V2.4/rho_22.12_23.86_212.97_7.3'
    iIOUN = kTempUnit
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='FORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
        WRITE(kStdErr,1010) IL, FNAME
        1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)
        CALL DoSTOP
    ENDIF

    iL = 0
    kTempUnitOpen = +1
    20 READ(iIOUN,5020,END=777) caLine
    iL = iL + 1
    READ(caLine,*) r1,r2,r3,r4,r5   !!wavenumber sun satellite d(phi) rho
    raW(iL) = r1
    raR(iL) = r5
    GOTO 20

    777 CONTINUE
    CLOSE(iIOUN)
    kTempUnitOpen = -1

    iFlip = -1    !!assume everything ordered correctly
    IF ((raW(iL-1) > raW(iL)) .OR. (raW(iL-2) > raW(iL-1))) THEN
        iFlip = +1
        DO iI = 1,iL
            raTemp(iL-iI+1) = raW(iI)
        END DO
        DO iI = 1,iL
            raW(iI) = raTemp(iI)
        END DO
        DO iI = 1,iL
            raTemp(iL-iI+1) = raR(iI)
        END DO
        DO iI = 1,iL
            raR(iI) = raTemp(iI)
        END DO
    END IF

    CALL rspl(raW,raR,iL, raFreq,raSpecularRefl,kMaxPts)
    write(kStdWarn,*) 'specular refl (1) = ',raFreq(1),raSpecularRefl(1)
    write(kStdWarn,*) 'for mu_angles ',r2*180/kPi,r3*180/kPi,r4*180/kPi

    5020 FORMAT(A80)

    RETURN
    end SUBROUTINE loadspecular

!************************************************************************
! this subroutine calculates the solar contribution AT GND
! ie take solar radiance incident from space at TOA, and then propagate down to surface
! then adjst by cos(SunAngle) * Earth-Sun solid angle
    SUBROUTINE Solar(iDoSolar,raSun,raFreq,raSunAngles, &
    iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raSun(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))
     
! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: r1,r2,rPlanck,rCos,raKabs(kMaxPts)
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!!! raSun will be in units of mW/m2/sr/cm-1 with NO sun solidangle correction
    IF (iDoSolar == 0) THEN
        write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
        ! ompute the Plank radiation from the sun
            rPlanck=exp(r2*raFreq(iFr)/rSunTemp)-1.0
            raSun(iFr) = r1*((raFreq(iFr))**3)/rPlanck
        END DO
    ELSEIF (iDoSolar == 1) THEN
        IF (raFreq(1) >= 605) THEN
            write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        ! ead in data from file
            CALL ReadSolarData(raFreq,raSun,iTag)
        ELSEIF (raFreq(1) < 605) THEN
        !! who cares, solar contribution is so small
            write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
            rSunTemp = kSunTemp
            DO iFr=1,kMaxPts
            ! compute the Plank radiation from the sun
                rPlanck=exp(r2*raFreq(iFr)/rSunTemp)-1.0
                raSun(iFr) = r1*((raFreq(iFr))**3)/rPlanck
            END DO
        END IF
    END IF

!! now do the solid angle correction
! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    rSunAngle = raSunAngles(MP2Lay(iaRadLayer(1)))
! change to radians
    rSunAngle = rSunAngle*kPi/180.0
    rCos      = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*rCos*rOmegaSun      !!!!this is correct
        raKAbs(iFr) = 0.0
    END DO

    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtraSun,raExtraSun)
      
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
    ! the current defined atmosphere used all Gnd-100 layers
        DO iLay = iNumLayer,2,-1
            iL = iaRadLayer(iLay)
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
            !!!!this is wrong!! raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)
                raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
            END DO
        END DO
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
                raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
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
            iI = iNumLayer
            write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                    raExtraSun(iFr) = raSun(iFr)*exp(-rakAbs(iFr))
                END DO
            END DO
        ! ow do remaining layers, all the way to the ground-1
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
            DO iLay=1,1
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO

        END IF
         
        IF (iT > iNumLayer) THEN
            write(kStdWarn,*)'need to do the upper layers as well!!'
        ! ow do top layers, all the way to the instrument
            DO iLay = iT,iNumLayer+1,-1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
        ! ow do the layer instrument is in
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                    raExtraSun(iFr) = raSun(iFr)*(exp(-raKabs(iFr)))
                END DO
            END DO
        ! ow do all the way to the ground-1
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
        ! ow do ground
            DO iLay=1,1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
    END IF

    RETURN
    end SUBROUTINE Solar
      
!************************************************************************
! this subroutine reads in the solar data files
!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;
!!!! then kCARTA  internally multiplies by 1000
!!!!  eventually also multiples by solidangle

    SUBROUTINE ReadSolarData(raFreq,raSun,iTag)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! raSun    = final solar contr in mW/m2/sr/cm-1
! raFreq  = frequency array
! iTag     = 1,2,3 and tells what the wavenumber spacing is
    INTEGER :: iTag
    REAL :: raSun(kMaxPts),raFreq(kMaxPts)
     
    CHARACTER(80) :: fname
    INTEGER :: iIOUN,iL,iU,iFr
    DOUBLE PRECISION :: fs,fe,df,daSun(kMaxPts)

    iIOUN = kTempUnit
    CALL GetSolarFileName(fname,raFreq(1))
    write(kStdWarn,*) 'solar data file = ',fname
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
        WRITE(kStdErr,1010) IL, FNAME
        1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)
        CALL DoSTOP
    ENDIF

    READ(iIOUN) fs,fe,df

    IF (abs(fs - raFreq(1)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1011) fs,raFreq(1)
        1011 FORMAT('ERROR! solar data file has start freq ',f10.5,' while the &
        start wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

    IF (abs(fe - raFreq(kMaxPts)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1012) fe,raFreq(kMaxPts)
        1012 FORMAT('ERROR! solar data file has stop freq ',f10.5,' while the &
        stop wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

    IF (abs(df - kaFrStep(iTag)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1013) df,kaFrStep(iTag)
        1013 FORMAT('ERROR! solar data file has delta freq ',f10.5,' while the &
        wavenumber spacing of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

! now map the data
    iL = iNT((fs-raFreq(1))/kaFrStep(iTag))
    iU = iL+kMaxPts

!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;
    READ(iIOUN) (daSun(iFr),iFr=1,kMaxPts)
    CLOSE(iIOUN)
!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;

!!! data files units of W cm-2 sr-1
!!! so need to multiply by 1000 to change to mW/m2/sr/cm-1
    DO iFr=1,kMaxPts
        raSun(iFr) = daSun(iFr)*1000.0
    !        write (6,2000) iFr,raFreq(iFr),raSun(iFr)
    END DO

    2000 FORMAT(I6,' ',f10.5,' ',e10.5)
    RETURN
    end SUBROUTINE ReadSolarData

!************************************************************************
! this subroutine checks to see how many radiances are to be output at this
! pressure layer
    SUBROUTINE DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp, &
    iOutNum)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iLay       = which of the radiating layers in atmosphere we are processing
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = list of fractions used for output for current atmosphere
! iOutNum    = of all options found in *OUTPUT, this pertains to current atmos
! raOutFrac  = list of fractions used if this layer has radiances to be output
! iDp        = number of fractional radiances to output
    REAL :: raOutFrac(kProfLayer),raaOp(kMaxPrint,kProfLayer)
    INTEGER :: iNp,iaOp(kPathsOut),iDp,iOutNum,iLay

! local variables
    INTEGER :: iDpC

    iDp=-1                   !assume nothing to be output

    IF (iNp < 0) THEN
    ! easy ! print the radiance at the end of this layer
        iDp=1
        raOutFrac(iDp)=1.0
    END IF

    IF (iNp > 0) THEN
        iDp=0
    ! actually have to go thru list to see if this layer is to be output
        iDpC=1
        101 CONTINUE
        IF (iaOp(iDpC) == iLay) THEN
            iDp = iDp+1
            raOutFrac(iDp) = raaOp(iOutNum,iDpc)
        END IF
        IF (iDpc < iNp) THEN
            iDpc = iDpc+1
            GO TO 101
        END IF
        IF (iDp == 0) THEN   !to make things oki doki, set no output to -1
            iDp = -1
        END IF
    END IF

    RETURN
    end SUBROUTINE DoOutPutRadiance

!************************************************************************
! this subroutine does the radiantive transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp
    SUBROUTINE RadianceInterPolate(iDir,rFrac,raFreq,raVTemp,rCos, &
    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun, &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iSun     = for uplooking instr, should we include sun in FOV?
! raSun    = for uplooking instr, should we include sun in FOV?
! raFreq  = wavenumbers
! iLay     = which layer in the list
! iaRadLayer = list of layers in atmosphere
! iDir     = direction of radiation travel (+1=downward look,-1=upward look)
! rFrac    = fraction of layer to use
! raVTemp  = mixed vertical temps
! rCos     = cos(satellite angle)
! raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
! raInten2 = interpolated radiation intensity at pressure level specified
! iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here
    INTEGER :: iDir,iLay,iaRadLayer(KProfLayer),iSun
    INTEGER :: iNumLayer,iProfileLayers
    REAL :: rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
    REAL :: raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

    IF (iDir < 0) THEN            !radiance going down to instr on gnd
        CALL UpLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos, &
        iLay,iaRadLayer,raaAbs,raInten,raSun,iSun, &
        iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)
    ELSE                             !radiance going up to instr in space
        CALL DownLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos, &
        iLay,iaRadLayer,raaAbs,raInten,raSun,iSun, &
        iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels, &
        iNLTEStart,raaPlanckCoeff)
    END IF

    RETURN
    end SUBROUTINE RadianceInterPolate

!************************************************************************
! this subroutine does the radiantive transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp
    SUBROUTINE UpLookInstrInterp(raInten2, &
    iDir,rFrac,raFreq,raVTemp,rCos, &
    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun, &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input parameters
!   iSun     = for uplooking instr, should we include sun in FOV?
!   raSun    = for uplooking instr, should we include sun in FOV?
!   raFreq  = wavenumbers
!   iLay     = which layer in the list
!   iaRadLayer = list of layers in atmosphere
!   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
!   rFrac    = fraction of layer to use
!   raVTemp  = mixed vertical temps
!   rCos     = cos(satellite angle)
!   raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
!   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here

! output parameters
!   raInten2 = interpolated radiation intensity at pressure level specified

    INTEGER :: iDir,iLay,iaRadLayer(KProfLayer),iSun
    INTEGER :: iNumLayer,iProfileLayers
    REAL :: rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
    REAL :: raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: raPressLevels(kProfLayer+1)
          
    INTEGER :: iFr,iL
    REAL :: rPlanck,rTrans,rEmis,rT,rFrac_k,rFrac_T
     
! iDir < 0  !radiance going down to instr on earth surface

! n all layers except bottommost layer, rFrac_k == rFrac_T
    rFrac_k=0.0            !use this much in k interpolation
    rFrac_T=0.0            !use this much in T interpolation

    iL = iaRadLayer(iLay)
    IF ((iLay > 1) .AND. (iLay < iNumLayer)) THEN
        rFrac_k = rFrac
        rFrac_T = rFrac
    ELSE IF (iLay == 1) THEN !!!topmost layer
    
    !====================== presslev(i1+1)

    ! --------------------- TopOfAtm
    ! XXXXXXXXXXXXXXXXXXXXX                         fraction rFracTop of full layer
    ! ---------------------  up look instr posn
    ! /////////////////////
    ! /////////////////////                        fraction rFrac of full layer
    !====================== presslev(i1)
        write(kStdWarn,*) 'recomputing fraction for top layer ...'
        rFrac_k = rFracTop-rFrac
    !!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
    !!!!!!!at surface, then we don't have to interpolate anything
        rFrac_T=(rFrac+rFracTop)/2.0  !!sort of do an average
        IF (rFrac/rFracTop > (1.0+1000*delta)) THEN
            write(kStdErr,*) rFrac,rFracTop
            write(kStdErr,*)'Cannot output radiance at such low'
            write(kStdErr,*)'pressure (topmost layer)'
            CALL DoStop
        END IF
    ELSE IF (iLay == iNumLayer) THEN !problem!!bottommost layer
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracBot > (1.0+1000*delta)) THEN
            write(kStdErr,*) rFrac,rFracBot
            write(kStdErr,*)'Cannot output radiance at such high'
            write(kStdErr,*)'pressure (bottommost layer)'
            CALL DoStop
        END IF
    ELSE
        write(kStdErr,*)'Cannot output radiance at this layer; not'
        write(kStdErr,*)'within atmosphere defined by user!!'
        CALL DoSTOP
    END IF

    IF (rFrac_k < 100*delta) THEN
        rFrac_k=0.00
    END IF
            
    write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

    IF (rFrac_k <= delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
            raInten2(iFr) = raInten(iFr)
        END DO
    ELSE                           !interpolate
        IF (iLay /= 1) THEN
        ! op part of most layers
            rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
        ELSE
        ! ottom  part of top layer
            rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
        END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
        DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rT)
            rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
            rEmis=(1.0-rTrans)*rPlanck
            raInten2(iFr) = rEmis+raInten(iFr)*rTrans
        END DO
        IF (iSun >= 0) THEN
            DO iFr=1,kMaxPts
                rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
                raInten2(iFr) = raInten2(iFr)+raSun(iFr)*rTrans
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE UpLookInstrInterp

!************************************************************************
! this subroutine does the radiative transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp
    SUBROUTINE DownLookInstrInterp(raInten2, &
    iDir,rFrac,raFreq,raVTemp,rCos, &
    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun, &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input parameters
!   raFreq  = wavenumbers
!   iLay     = which layer in the list
!   iaRadLayer = list of layers in atmosphere
!   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
!   rFrac    = fraction of layer to use
!   raVTemp  = mixed vertical temps
!   rCos     = cos(satellite angle)
!   raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
!   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here
! output parameters
!   raInten2 = interpolated radiation intensity at pressure level specified

    INTEGER :: iDir,iLay,iaRadLayer(KProfLayer),iSun
    INTEGER :: iNumLayer,iProfileLayers
    REAL :: rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
    REAL :: raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iL
    REAL :: rPlanck,rTrans,rEmis,rT,rFrac_k,rFrac_T

! iDir > 0  !radiance going up to instr in space
     
! n all layers except bottommost layer, rFrac_k == rFrac_T
    rFrac_k = 0.0            !use this much in k interpolation
    rFrac_T = 0.0            !use this much in T interpolation

    iL = iaRadLayer(iLay)

    IF ((iLay > 1) .AND. (iLay < iNumLayer)) THEN
    ! o problem; full layer in mixtable
        rFrac_k = rFrac
        rFrac_T = rFrac
    ELSE IF (iLay == 1) THEN !!!bottommost layer
    
    !====================== presslev(i1+1)
    ! XXXXXXXXXXXXXXXXXXXXX                         fraction rFrac of full layer
    ! ---------------------  down look instr posn
    ! /////////////////////
    ! /////////////////////                        fraction rFracBot of full layer
    ! --------------------- surface
    
    
    !====================== presslev(i1)
        write(kStdWarn,*)'recomputing fraction for bottom layer ...'
    ! ->> this is old
    ! ->>        rFrac_k = rFracBot-rFrac
    !!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
    !!!!!!!at surface, then we don't have to interpolate anything
    ! ->>        rFrac_T = (rFrac+rFracBot)/2.0  !!sort of do an average
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracBot > (1.0+1000*delta)) THEN
            write(kStdErr,*) rFrac,rFracBot
            write(kStdErr,*)'Cannot output radiance at such high'
            write(kStdErr,*)'pressure (bottommost layer)'
            CALL DoStop
        END IF
    ELSE IF (iLay == iNumLayer) THEN !problem!!top most layer
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracTop > (1.0+1000*delta)) THEN
            write(kStdErr,*) rFrac,rFracTop
            write(kStdErr,*)'Cannot output radiance at such low'
            write(kStdErr,*)'pressure (topmost layer)'
            CALL DoStop
        END IF
    ELSE
        write(kStdErr,*)'Cannot output radiance at this layer; not'
        write(kStdErr,*)'within atmosphere defined by user!!'
        CALL DoSTOP
    END IF

    IF (rFrac_k < 100*delta) THEN
        rFrac_k = 0.00
    END IF
            
    write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

    IF (rFrac_k <= delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
            raInten2(iFr) = raInten(iFr)
        END DO
    ELSE                           !interpolate
        IF (iLay /= 1) THEN
        ! ottom part of most layers
            rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
        ELSE
        ! op part of bottom layer
            rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
        END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
    ! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
    ! so usually only the usual LTE computations are done!!
        IF (iNLTEStart > kProfLayer) THEN    !!!normal, no emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rT)
                rTrans = exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
                rEmis  = (1.0-rTrans)*rPlanck
                raInten2(iFr) = rEmis+raInten(iFr)*rTrans
            END DO
        ELSE IF (iNLTEStart <= kProfLayer) THEN
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rT)
                rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
                rEmis=(1.0-rTrans)*rPlanck*raaPlanckCoeff(iFr,iL)
                raInten2(iFr) = rEmis+raInten(iFr)*rTrans
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE DownLookInstrInterp

!************************************************************************
! this subroutine mimics the SARTA NLTE model
! see /asl/packages/sartaV106/Src/calnte.f
    SUBROUTINE Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
    iaRadLayer,iNumLayer,raInten,rCO2MixRatio)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: SUNCOS ! solar zenith angle cosine at surface
    REAL ::  SCOS1 ! solar zenith angle cosine at layer1
    REAL ::  VSEC1 ! satellite view zenith angle secant at layer1
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    REAL :: rCO2MixRatio !!! CO2 mixing ratio

! input/output
    REAL :: raInten(kMaxPts)

! local vars
    INTEGER :: iI,iJ,iC
    INTEGER :: NCHNTE,iERR,iIOUN,ICHAN

!      !!/asl/packages/sartaV106/Src/incFTC_apr05_nte.f
!      INTEGER MXCHNN,NNCOEF
!      !MXCHNN = 201
!      !NNCOEF = 6     !/asl/packages/sartaV106/Src/incFTC_apr05_nte.f
!      PARAMETER(MXCHNN = 201,NNCOEF = 6)
!      REAL COEFN(NNCOEF,MXCHNN)

!     !!/asl/packages/sartaV108/Src/incFTC_airs_apr08_m130_m150_template_cal.f
    INTEGER :: MXCNTE ! max # of channels for non-LTE (203)
    INTEGER :: NNCOEF ! # of coefs for non-LTE
    INTEGER :: NTEBOT ! bottom layer for CO2TOP calc
    REAL :: CO2NTE ! ref CO2 mixing ratio for non-LTE coefs (ppmv)
    PARAMETER(MXCNTE = 203)
    PARAMETER(NNCOEF = 7)
    PARAMETER(NTEBOT = 10)
    PARAMETER(CO2NTE = 370.0)
    REAL :: COEFN(NNCOEF,MXCNTE)
         
    CHARACTER(120) :: FNCOFN

    REAL :: tHigh,pred1,pred2,pred3,pred4,pred5,pred6,FRQCHN
    REAL :: raDrad(kMaxPts),raFrad(kMaxPts),raDkcarta(kMaxPts)
    REAL :: CO2TOP

    tHigh = 0.0
    DO iI = iNumLayer-4,iNumLayer
        tHigh = tHigh + raVTemp(iaRadLayer(iI))
    END DO
    tHigh = tHigh/5

    PRED1 = 1.0
    PRED2 = SCOS1
    PRED3 = SCOS1*SCOS1
    PRED4 = SCOS1*VSEC1
    PRED5 = SCOS1*THIGH
    PRED6 = SUNCOS

!      !/asl/packages/sartaV106/Src/incFTC_apr05_nte.f
!      FNCOFN= '/asl/data/sarta_database/Data_jan04untun/Coef/setnte_oct05.dat'
!      need to process this : from be to le
!      !matlab cd : /home/sergio/KCARTADATA/NLTE/SARTA_COEFS/
!      !fname='/asl/data/sarta_database/Data_jan04untun/Coef/setnte_oct05.dat';
!      ![idchan, freq, coef] = rd_nte_be(fname);
!      !/wrt_nte_le.m
!      FNCOFN= '/home/sergio/KCARTADATA/NLTE/SARTA_COEFS/setnte_oct05.le.dat'
!      FNCOFN= '/asl/data/kcarta/KCARTADATA/NLTE/SARTA_COEFS/nonLTE7_m150.le.dat'

! /asl/packages/sartaV108/Src_rtpV201/incFTC_iasi_sep08_wcon_nte.f :
!       PARAMETER(FNCOFN=
!            $ '/asl/data/sarta_database/Data_IASI_sep08/Coef/nte_7term.dat')

    FNCOFN = kSartaNLTE
    iIOUN = kTempUnit

    OPEN(UNIT=iIOUN,FILE=FNCOFN,FORM='UNFORMATTED',STATUS='OLD', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'error reading SARTA NLTE coeffs file',IERR, FNCOFN
        CALL DoStop
    ENDIF

    kTempUnitOpen = +1
    iJ=1
    DO iI=1,MXCNTE
    !       Read data for this frequency/channel
        READ(iIOUN) ICHAN, FRQCHN, (COEFN(IC,iJ),IC=1,NNCOEF)
        raFrad(iI) = FRQCHN
        iJ = iJ + 1
    ENDDO
    NCHNTE = iJ - 1
    CLOSE(iIOUN)
    kTempUnitOpen = -1

!      print *,kSartaNLTE
!      print *,MXCNTE,NCHNTE,NNCOEF
!      call dostop

    CO2TOP = kCO2ppmv      !! this is the expected CO2 at TOA
    CO2TOP = rCO2MixRatio  !! better to use this, as it comes from profile

! asl/packages/sartaV106or8/Src/calnte.f
! ote we assume the freqs raFrad are SORTED so we don't have problems
! oing the spline interpolation onto kCARTA freqs
    DO iI = 1,NCHNTE
        raDrad(iI)=( COEFN(1,iI)*PRED1 ) + &
        ( COEFN(2,iI)*PRED2 ) + &
        ( COEFN(3,iI)*PRED3 ) + &
        ( COEFN(4,iI)*PRED4 ) + &
        ( COEFN(5,iI)*PRED5 ) + &
        ( COEFN(6,iI)*PRED6 )

    ! this is new; see calnte.f in eg sartaV108/Src_rtpV201
    ! this accounts for changing concentration of CO2
    !       Adjust DRAD for CO2 mixing ratio
    !       DRAD=DRAD*(COEFN(7,I)*(CO2TOP - CO2NTE) + 1.0)
                
        raDrad(iI) = raDrad(iI)*(COEFN(7,iI)*(CO2TOP - CO2NTE) + 1.0)
    !       print *,iI,raFrad(iI),raDrad(iI)
    END DO

!      print *,'sarta nlte oops'
!      print *,'sartanlte pred : ',PRED1,PRED2,PRED3,PRED4,PRED5,PRED6,COEFN(7,iI-1),rCO2MixRatio,CO2TOP,CO2NTE
!      CALL DoStop

    CALL rspl(raFrad,raDrad,NCHNTE,raFreq,raDkcarta,kMaxPts)    !! too dangerous,   small 4 um lte rads, wiggly NLTE correction
    CALL rlinear(raFrad,raDrad,NCHNTE,raFreq,raDkcarta,kMaxPts) !! hopefully safer, small 4 um lte rads, straightline NLTE correction
    DO iI = 1,kMaxPts
        IF ((raFreq(iI) < raFrad(1)) .OR. (raFreq(iI) > raFrad(NCHNTE))) THEN
            raDkcarta(iI) = 0.0
        END IF
    !       print *,iI,raFreq(iI),raDkcarta(iI),raInten(iI)
    !! LTE rad is raInten(iI); NLTE correction raDkcarta(iI) should be positive; if negative then just stick to the LTE rad
        raInten(iI) = max(raInten(iI) + raDkcarta(iI),raInten(iI))
    END DO

    RETURN
    end SUBROUTINE Sarta_NLTE
          
!************************************************************************
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
    SUBROUTINE AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtra,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! input params
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
!                  temporarily define atm from GND to TOP of atmosphere
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    REAL :: rFracTop
! output params
! raExtra = array initialized to all zeros
! iExtra = -1 if no layeres added on, +1 if layers added on
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
! iT             = number of layers in this temporary atmosphere
    INTEGER :: iaRadLayerTemp(kMixFilRows),iExtra,iT
    REAL :: raExtra(kMaxPts)
     
    INTEGER :: iI,iFr
     
    iExtra = -1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4)) THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
        write(kStdWarn,*) 'top most layer is fractional layer. Some'
        write(kStdWarn,*) 'portion needed above instrument to calculate'
        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
            write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF
     
    RETURN
    end SUBROUTINE AddUppermostLayers

!************************************************************************
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
! same as above routine, except that it is quiet!!!!!!!! (for scatttering code)
    SUBROUTINE  AddUppermostLayersQ(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtra,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
!                  temporarily define atm from GND to TOP of atmosphere
! iT             = number of layers in this temporary atmosphere
! iExtra = -1 if no layeres added on, +1 if layers added on
! raExtra = array initialized to all zeros
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows),iExtra
    REAL :: raExtra(kMaxPts),rFracTop
     
    INTEGER :: iI,iFr
     
    iExtra=-1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4))THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    !        write(kStdWarn,*) 'top most layer is fractional layer. Some'
    !        write(kStdWarn,*) 'portion needed above instrument to calculate'
    !        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
        !          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
    !        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
    !        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF
     
    RETURN
    END SUBROUTINE 

!************************************************************************
! allows a NLTE atmosphere, does SLOW ACCURATE GENLN2 calc

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles, &
    iTag,raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer),rCO2MixRatio
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2,rOmegaSun
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    INTEGER :: iCloudLayerTop,iCloudLayerBot
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,troplayer
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
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
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar

    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO

! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'slow nlte radmodel ...'
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh=iaOp(iLay)
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
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
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
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
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
    END DO
          
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    IF (iNLTEStart > kProfLayer) THEN
        iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
        write (kStdErr,*) 'Normal rad transfer .... no NLTE'
        write (kStdErr,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
        write (kStdErr,*) 'should be calling rad_trans_SAT_LOOK_DOWN'
        CALL DoStop
    ELSE
        iLay = 1
        987 CONTINUE
        iL=iaRadLayer(iLay)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF ((iLModKprofLayer < iNLTEStart) .AND. (iLay < iNumLayer)) THEN
            iLay = iLay + 1
            GOTO 987
        END IF
        iSTopNormalRadTransfer = iLay
        write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
    END IF

    DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF (iLModKprofLayer < iNLTEStart) THEN
        ! ormal, no LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iLModKprofLayer >= iNLTEStart) THEN
        ! ew; LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
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

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    DO iFr=1,kMaxPts
        raInten(iFr)=raSurface(iFr)*raUseEmissivity(iFr)+ &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
        raSun(iFr)*raSunRefl(iFr)
    END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    DO iLay=iHigh,iHigh
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

        IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp >= 1) THEN

                write(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer'
                write(kStdWarn,*) 'This is the top of the usual AIRS atmosphere'
                write(kStdWarn,*) '   you have iDoUpperATM > 0'
                write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
                write(kStdWarn,*) 'and output stuff into the blah_UA file'
                write(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of'
                write(kStdWarn,*) 'stratosphere into both this and the UA file'

            ! o radiative transfer thru this layer
                DO iFr=1,kMaxPts
                    raInten(iFr) = &
                    raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
                END DO

            ! ow do complete rad transfer thru upper part of atmosphere
                CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)), &
                iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of stratosphere (120-200 km)

                write(kStdWarn,*) 'finally outputting radiances at TOTAL Complete TOA into'
                write(kStdWarn,*) 'usual binary file (iLay = ',iLay,')'

                DO iFr=1,iDp
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten)
                END DO
            END IF
        END IF

        IF (iUpper < 1) THEN
        !!! no need to compute stuff at extra layers (100-200 km)
        !!! so just do usual stuff
        !!! see if this mixed path layer is in the list iaOp to be output
        !!! since we might have to do fractions!
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        END IF

    !c no need to do radiative transfer thru this layer
    !c        DO iFr=1,kMaxPts
    !c          raInten(iFr)=raaEmission(iFr,iLay)+
    !c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
    !c          END DO

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_SLOW

!************************************************************************
! allows a NLTE atmosphere, does FAST SARTA NLTE calc

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FAST(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles, &
    iTag,raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2,rOmegaSun
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    INTEGER :: iCloudLayerTop,iCloudLayerBot
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,troplayer
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    REAL :: SUNCOS ! solar zenith angle cosine at surface
    REAL ::  SCOS1 ! solar zenith angle cosine at layer1
    REAL ::  VSEC1 ! satellite view zenith angle secant at layer1

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
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
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'fast nlte radmodel ...'
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh=iaOp(iLay)
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
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
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
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
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
    END DO
          
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    IF (iNLTEStart <= kProfLayer) THEN
        write (kStdWarn,*) 'this routine expects SARTA NLTE',iNLTEStart
        CALL DoStop
    ELSE
    !        write (kStdWarn,*) 'Normal rad transfer .... '
    !        write (kStdWarn,*) 'adding on SARTA NLTE if DAYTIME'
        iLay = 1
        iSTopNormalRadTransfer = kProfLayer + 1
    END IF

    DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
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
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
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

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    iSpecular = +1    !some specular refl, plus diffuse
    iSpecular = -1    !no   specular refl, only diffuse
    IF (iSpecular > 0) THEN
        write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN_NLTE_FAST'
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

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    DO iLay=iHigh,iHigh
    ! o rad transfer to TOA (80 km)
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+ &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

        suncos = raSunAngles(iaRadLayer(1))           !! at surface
        scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
        vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

        suncos = cos(suncos*kPi/180.0)
        scos1  = cos(scos1*kPi/180.0)
        vsec1  = 1/cos(vsec1*kPi/180.0)

        IF (iDoSolar >= 0) THEN
        !         do iFr = 1,iNumlayer
        !           write(kStdWarn,*) iFr,iaRadLayer(iFr),raSunAngles(iaRadLayer(iFr))
        !           end do
            write(kStdWarn,*)'day .... add SARTA_NLTE at solangle for chunk ',raSunAngles(iaRadLayer(1)),raFreq(1)
            CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
            iaRadLayer,iNumlayer,raInten,rCO2MixRatio)
        ELSE
            write(kStdWarn,*)'nighttime ... do not need SARTA_NLTE'
        END IF

        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FAST

!************************************************************************
! this file reads a binary made from the ASCII sscatmie.x file
! and returns the extinction, absm asymmetry coeffs
! this is a combination of subroutines
!      INTERP_SCAT_TABLE2 and READ_SSCATTAB_BINARY
    SUBROUTINE FIND_ABS_ASY_EXT(SCATFILE,DME,IWP,pT,pB,raPLevels,RAFREQ, &
    iaRadLayer,iNumLayer,EXTINCT,ABSC,ASYM,ILT,ILB)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

!       Input parameters:
!     SCATFILE   file name of scattering file
!     DME        particle size to interpolate for
!     IWP        iwp normalization
!     WAVES      wavenumbers
!     pT,pB      pressure (top + bottom) where the cloud layer is
!     raPLevels        AIRS pressure levels
!     iaRadLayer current atmosphere layers
!       Output parameters:
!     EXTINCT, ABS, ASYM  : the particle scattering coefficients for each layer

    INTEGER :: iaRadlayer(kProfLayer),iNumLayer
    CHARACTER*(*) SCATFILE
    REAL :: raFreq(kMaxPts), DME, IWP, pT, pB, raPLevels(kProfLayer+1)
    REAL :: extinct(kMaxPts),absc(kMaxPts),asym(kMaxPts)

! output layer
!     IL                  : which AIRS layer the cloud is in
    INTEGER :: iLT,iLB

! local variables
    CHARACTER(1) :: caScale(MAXSCAT)
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
          
    INTEGER :: I,IF,iMod,iS,iL
    REAL :: ee,aa,gg, waveno

    I = 1

    CALL READ_SSCATTAB_BINARY(SCATFILE,   & !!!!!!MAXTAB, MAXGRID,
    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
    NWAVETAB(I), WAVETAB(1,I), &
    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
    TABPHI1UP(1,I), TABPHI1DN(1,I), &
    TABPHI2UP(1,I), TABPHI2DN(1,I))

!       !!!get rid of delta scaling
!      CALL UnScaleMie(
!     $        caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
!     $        ndme(i)*nwavetab(i))
                 
    DO iF = 1,kMaxPts
        waveno = raFreq(iF)
    !  here we only need the simpler first choice as we are not messing
    !  around with the phase functions
        CALL INTERP_SCAT_TABLE2 (WAVENO, DME, ee, aa, gg, &
        NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
        TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
        EXTINCT(iF) = ee * iwp/1000.0
        ABSC(iF)    = ee * iwp/1000.0 * (1.0 - aa)
        ASYM(iF)    = gg
    END DO

!     figure out what AIRS layers the cloud is in between

! do the top layer --------------------------------->
    iL = 1
    10 CONTINUE
    IF (raPLevels(iL) <= 1.0e-3) THEN
        iL = iL + 1
        GOTO 10
    END IF

    IF (pT > raPLevels(iL)) THEN
        write(kStdErr,*) 'cloud top pressure (',pT,' mb) is too large!!!'
        CALL DoStop
    END IF
    IF (pT < raPLevels(kProfLayer+1)) THEN
        write(kStdErr,*) 'cloud top pressure (',pT,' mb) is too small!!!'
        CALL DoStop
    END IF

    iL = 1
    20 CONTINUE
    IF ((pT <= raPLevels(iL)) .AND. (pT >= raPLevels(iL+1))) THEN
        GOTO 30
    ELSE
        iL = iL + 1
        GOTO 20
    END IF
          
    30 CONTINUE

    IF ((iL < 1) .OR. (iL > kProfLayer)) THEN
        write(kStdErr,*) 'iL = ',iL,' ... out of range!!!'
        CALL DoStop
    END IF

!!!now see how this can be put into iaRadLayer
! figure out maximum Mixed Path Layer in the atmosphere
    IF (iaRadlayer(1) > iaRadLAyer(iNumLayer)) THEN
        iS = iaRadlayer(1)
    ELSE
        iS = iaRadLAyer(iNumLayer)
    END IF
    iMod = 1
    40 CONTINUE
    IF ((iMod * kProfLayer) < iS) THEN
        iMod = iMod + 1
        GOTO 40
    END IF
!!!so, this is the Mixed Path Layer with Cloud in it
    iL = (iMod-1)*kProfLayer + iL
!!!now see which iaRadLayer this corresponds to
    iS = 1
    50 CONTINUE
    IF ((iaRadLayer(iS) /= iL) .AND. (iS <= iNumLayer)) THEN
        iS = iS + 1
        GOTO 50
    END IF

    iL = iS
    write(kStdWarn,*) '  Putting top of abs cloud into iaRadLayer(',iL,')'

    iL = iaRadLayer(1) + iL - 1
    write(kStdWarn,*) '    which is MP radiating layer ',iL

    write(kStdWarn,*) '  This is for cloud pressure = ',pT
    write(kStdWarn,*) '  Corresponding AIRS levels are : ',raPLevels(iL), &
    raPLevels(iL+1)

    iLT = iL

! do the bottom layer --------------------------------->
    iL = 1
    15 CONTINUE
    IF (raPLevels(iL) <= 1.0e-3) THEN
        iL = iL + 1
        GOTO 15
    END IF

    IF (pB > raPLevels(iL)) THEN
        write(kStdErr,*) 'cloud bot pressure (',pB,' mb) is too large!!!'
        CALL DoStop
    END IF
    IF (pB < raPLevels(kProfLayer+1)) THEN
        write(kStdErr,*) 'cloud bot pressure (',pB,' mb) is too small!!!'
        CALL DoStop
    END IF

    iL = 1
    25 CONTINUE
    IF ((pB <= raPLevels(iL)) .AND. (pB >= raPLevels(iL+1))) THEN
        GOTO 35
    ELSE
        iL = iL + 1
        GOTO 25
    END IF
          
    35 CONTINUE

    IF ((iL < 1) .OR. (iL > kProfLayer)) THEN
        write(kStdErr,*) 'iL = ',iL,' ... out of range!!!'
        CALL DoStop
    END IF

!!!now see how this can be put into iaRadLayer
! figure out maximum Mixed Path Layer in the atmosphere
    IF (iaRadlayer(1) > iaRadLAyer(iNumLayer)) THEN
        iS = iaRadlayer(1)
    ELSE
        iS = iaRadLAyer(iNumLayer)
    END IF
    iMod = 1
    45 CONTINUE
    IF ((iMod * kProfLayer) < iS) THEN
        iMod = iMod + 1
        GOTO 45
    END IF
!!!so, this is the Mixed Path Layer with Cloud in it
    iL = (iMod-1)*kProfLayer + iL
!!!now see which iaRadLayer this corresponds to
    iS = 1
    55 CONTINUE
    IF ((iaRadLayer(iS) /= iL) .AND. (iS <= iNumLayer)) THEN
        iS = iS + 1
        GOTO 55
    END IF

    iL = iS
    write(kStdWarn,*) '  Putting bot of abs cloud into iaRadLayer(',iL,')'

    iL = iaRadLayer(1) + iL - 1
    write(kStdWarn,*) '    which is MP radiating layer ',iL

    write(kStdWarn,*) '  This is for cloud pressure = ',pB
    write(kStdWarn,*) '  Corresponding AIRS levels are : ',raPLevels(iL), &
    raPLevels(iL+1)

    iLB = iL

! see if the layers make sense
    IF (iLB > iLT) THEN
        write(kStdErr,*) 'oops in FIND_ABS_ASY_EXT iLB > iLT',iLB,iLT
        CALL DOStop
    END IF

! see if we need to adjust the individual cloud opt depths
    IF (iLB /= iLT) THEN
        write(kStdWarn,*) 'adjusting the cld abs depths for each layer'
        DO iF = 1,kMaxPts
            EXTINCT(iF) = EXTINCT(iF)/(iLT-iLB+1)
            ABSC(iF)    = ABSC(iF)/(iLT-iLB+1)
            ASYM(iF)    = ASYM(iF)
        END DO
    END IF

    RETURN
    END SUBROUTINE FIND_ABS_ASY_EXT

!************************************************************************
    SUBROUTINE INTERP_SCAT_TABLE2_modified (WAVENO, DME, &
    EXTINCT, SSALB, ASYM, &
    NDME, DMETAB, NWAVE, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM)
!       Interpolates the scattering properties from the table for
!     a particular wavenumber and particle size.  Does a bilinear
!     interpolation, but optimized for fixed particle size and slowly
!     varying wavenumber.  If the DME is the same as last time then we
!     can just linearly interpolate in wavenumber between stored
!     scattering values.  If the DME has changed then we linearly
!     interpolate between the DMETAB grid lines for each of the two
!     wavenumber grid lines.

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL ::     WAVENO, DME
    REAL ::     EXTINCT, SSALB, ASYM
    INTEGER ::  NDME, NWAVE
    REAL ::     DMETAB(NDME), WAVETAB(NWAVE)
    REAL ::     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME)
    REAL ::     TABASYM(NWAVE,NDME)
    INTEGER ::  IW0, IW1, ID, IL, IU, IM
    LOGICAL ::  NEWIW
    REAL ::     FWAV, FDME, FLDME, F
    REAL ::     OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
! sergio do not save iw0,iw1, olddme
!      SAVE     IW0, IW1, ID, OLDDME, FDME, FLDME
!      SAVE     ID, FDME, FLDME
!      SAVE     EXT0,EXT1, ALB0,ALB1, ASYM0,ASYM1
    DATA     IW0/1/, IW1/2/

    iw0 = 1
    iw1 = 2
    olddme = 0.0

    iw0 = 1
    iw1 = nwave
    olddme = -10.0
          
!         Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
        write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
        IF (WAVENO < WAVETAB(1)) THEN
            WAVENO = WAVETAB(1)
        ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
            WAVENO = WAVETAB(NWAVE)
        END IF
    ! ALL DoStop
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
        write(kStdErr,*) DME,' outside ',DMETAB(1),':',DMETAB(NDME)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: particle Dme out of range ... RESET'
        IF (DME < DMETAB(1)) THEN
            DME = DMETAB(1)
        ELSEIF (DME > DMETAB(NDME)) THEN
            DME = DMETAB(NDME)
        END IF
    ! ALL DoStop
    END IF

!         See if wavenumber is within last wavenumber grid, otherwise
!           find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
!      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (WAVENO >= WAVETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE. 
    ENDIF

    IF (DME /= OLDDME) THEN
    !         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (DME >= DMETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
    ENDIF

    IF (DME /= OLDDME .OR. NEWIW) THEN
    !         If not the same Dme or a new wavenumber grid, then
    !           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)
    ENDIF

!         Linearly interpolate the scattering properties in wavenumber
    FWAV    = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
    F       = 1-FWAV
    EXTINCT = F*EXT0 + FWAV*EXT1
    SSALB   = F*ALB0 + FWAV*ALB1
    ASYM    = F*ASYM0 + FWAV*ASYM1

    OLDDME = DME

    RETURN
    end SUBROUTINE INTERP_SCAT_TABLE2_modified

!************************************************************************
! set the vertical temperatures of the atmosphere for TWOSTREAM
! this sets the temperatures at the pressure level boundaries, using the
! temperatures of the pressure layers that have been supplied by kLayers
    SUBROUTINE SetTWOSTRTemp(TEMP,iaRadLayer,raVTemp,iNumLayer, &
    iDownWard,iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are variables that come in from kcartamain.f
    REAL :: raVTemp(kMixFilRows),raPressLevels(kProfLayer+1)
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,iDownWard,iProfileLayers
! these are variables that we have to set
    REAL ::    TEMP(*)

! local variables
    INTEGER :: iL,iLay,iM,iaRadLayerTemp(kMixFilRows),iOffSet,iJump,iLowest
    REAL :: Temp1(maxnz)
    REAL :: pavg(kProfLayer),rP,raProfileTemp(kProfLayer)

    iLowest = kProfLayer - iProfileLayers + 1

    DO iLay=1,MAXNZ
        Temp1(iLay) = -10.0
        Temp(iLay)  = -10.0
    END DO

    DO iLay = iLowest,kProfLayer
        pavg(iLay) = raPressLevels(iLay+1)-raPressLevels(iLay)
        pavg(iLay) = pavg(iLay)/log(raPressLevels(iLay+1)/raPressLevels(iLay))
    END DO

! now set iaRadLayerTemp the same as  iaRadLayer if downlook instr
!     set iaRadLayerTemp flipped from iaRadLayer if uplook   instr
    IF (iDownWard == 1) THEN      !!!!keep everything the same
        DO iLay = 1,iNumLayer
            iaRadLayerTemp(iLay) = iaRadLayer(iLay)
        END DO
    ELSE            !!!gotta do a bit of reverse logic for uplook instr
        DO iLay = 1,iNumLayer
            iaRadLayerTemp(iLay) = iaRadLayer(iNumLayer-iLay+1)
        END DO
    END IF

! see which set of Mixed Paths the current atmosphere occupies eg
! set 1 = 1..100, set2= 101..200 etc
! eg if current atmosphere is from MixfilPath 110 to 190, and kProfLayer = 100,
! then we set iMod as 2      idiv(150,100) = 1  === 2nd set of mixed paths
! assume each atmosphere has at least 25 layers in it!!!
    iM = idiv(iaRadLayerTemp(25),kProfLayer)+1
    DO iLay=1,kProfLayer
        raProfileTemp(iLay) = raVTemp(iLay+(iM-1)*kProfLayer)
    END DO

    DO iLay=1,iNumLayer
        iL = iaRadLayerTemp(iLay)
    ! ap this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
        iL = iL-idiv(iL,kProfLayer)*kProfLayer
        IF (iL == 0) THEN
            iL = kProfLayer
        END IF
        rP=raPressLevels(iL+1)-10000*delta
        if (rp < raPressLevels(kProfLayer+1)) then
            rp = raPressLevels(kProfLayer+1)+10000*delta
        end if
        TEMP1(iNumLayer-iLay+1) = FindBottomTemp(rP,raProfileTemp, &
        raPressLevels,iProfileLayers)
    END DO

    rP = raPressLevels(iLowest)
    rP = DISORTsurfPress          !!!from scatterparam.f90
    TEMP1(iNumLayer+1) = FindBottomTemp(rP,raProfileTemp, &
    raPressLevels,iProfileLayers)

    IF (iDownWard == 1) THEN
        DO iLay=1,iNumLayer+1
            temp(iLay) = temp1(iLay)
        END DO
    ELSE
        DO iLay=1,iNumLayer+1
            temp(iLay) = temp1((iNumLayer+1)-iLay+1)
        END DO
    END IF

    IF (iDownWard == -1) THEN
    !!!suppose atm is in kCARTA layers 5 -- 100 (96 layers ==> 97 levels)
    !!!this is same as RTSPEC levels 1 -- 97 ...
    !!!   so temp(iI) is filled from levels 1 ..97
    !!!so push it up so that it occupies KLAYERS levels 5 ... 101

    !!!set up the temp1 array
        DO iLay=1,kProfLayer+1
            temp1(iLay) = temp(iLay)
            temp(iLay)  = -10.0
        END DO

    !!!push up the stuff so it occupies kLAYERS levels 5 .. 80
        iOffSet = (iaRadLayer(iNumLayer) - 1) - (iM-1)*kProfLayer
        DO iLay = 1,iNumLayer+1
            temp(iLay+iOffSet) = temp1(iLay)
        END DO
    END IF

    IF (iDownWard == 1) THEN
    !!!suppose atm is in kCARTA layers 5 -- 79 (75 layers ==> 76 levels)
    !!!this is same as RTSPEC levels 1 -- 76 ...
    !!!   so temp(iI) is filled from levels 1 ..76
    !!!so now push this down so it fills RTSPEC levels 26 .. 101
    !!!then flip it so it occupies KLAYERS levels 1 ... 76
    !!!and then push it up so it occupies KLAYERS levels 5 -- 80
    !!!   (which is same as KLAYERS layers  5-79!!!!)

    !!!set up the temp1 array
        DO iLay=1,kProfLayer+1
            temp1(iLay) = temp(iLay)
            temp(iLay)  = -10.0
        END DO

    !!!push it down so it occupies RTSPEC levels 26 .. 101
        iOffSet = kProfLayer-iNumLayer
        DO iLay=1,iNumLayer+1
            temp(iLay+iOffSet) = temp1(iLay)
        END DO

    !!!now flip it so it occupies KLAYERS levels 1 ..76
        DO iLay = 1,kProfLayer + 1
            TEMP1(iLay) = TEMP(iLay)
        END DO
        DO iLay = 1,kProfLayer + 1
            TEMP(iLay) = TEMP1((kProfLayer+1)-iLay+1)
        END DO
        DO iLay=1,kProfLayer+1
            temp1(iLay) = temp(iLay)
            temp(iLay)  = -10.0
        END DO

    !!!push up the stuff sp it occupies kLAYERS levels 5 .. 80
        iOffSet = iaRadLayer(1) - 1 - (iM-1)*kProfLayer
        DO iLay = 1,iNumLayer+1
            temp(iLay+iOffSet) = temp1(iLay)
        END DO
    END IF

    RETURN
    end SUBROUTINE SetTWOSTRTemp

!************************************************************************
! this subroutine resets the TEMPerature array that comes out of
! GetAbsProfileRTSPEC : so that it is same as raVertTemp
! this is because
! RTSPEC will want stuff from RTSPEC layerM --> layerN and ignore N+1 to 100
! so this code is a little bit smart and reset temps so they are ok
    SUBROUTINE ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
    iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! output variable
    REAL :: TEMP(MAXNZ)     !temperature of layers, in kCARTA layering style
!1 = GND, 100 = TOA
! input variables
    REAL :: raVTemp(kMixFilRows),rSurfaceTemp,raPressLevels(kProfLayer+1)
    INTEGER :: iDownWard,iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iProfileLayers

    INTEGER :: iii,iaRadLayer(kProfLayer)
    REAL :: TEMP1(MAXNZ)

    DO iii = 1,iNumLayer
        iaRadLayer(iii) = iaaRadLayer(iAtm,iii)
    END DO

    CALL SetTWOSTRTemp(TEMP,iaRadLayer,raVTemp,iNumLayer, &
    iDownWard,iProfileLayers,raPressLevels)

    RETURN
    end SUBROUTINE ResetTemp_Twostream
!************************************************************************
! this function finds the pressure layer at which rPressStart is within,
! as well as the fraction of the layer that it occupies
    REAL FUNCTION FindBottomTemp(rP,raProfileTemp, &
    raPressLevels,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raPressLevels = actual pressure levels that come out of kLAYERS
! raProfileTemp = actual profile temp
! rP            = pressure at which we want the temperature
    real :: rP,raProfileTemp(kProfLayer),raPressLevels(kProfLayer+1)
    integer :: iProfileLayers

    integer :: iFound,i1,i2,i3,iLowest,iJ
    real :: rP1,rP2,T1,T2
    real :: raP(3),raT(3),Y2A(3),rT,raLogP(3)
    real :: yp1,ypn,work(3)
    INTEGER :: iLog,iSpline

    iLog = +1       !!!do log(P) for the x-points
    iLog = -1       !!!do   P    for the x-points

    iSpline = -1    !!!use linear interpolations
    iSpline = +1    !!!use spline interpolations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iLog = +1       !!!do log(P) for the x-points
    iSpline = -1    !!!use linear interpolations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rT=0.0

    iLowest = kProfLayer-iProfileLayers+1

    IF (rP >= raPressLevels(iLowest)) THEN
    ! his is WHOLE of the bottom layer
        i1=iLowest
    ELSE IF (rP <= raPressLevels(kProfLayer+1)) THEN
    ! his is ludicrous
        write(kStdErr,*) rP,raPressLevels(kProfLayer+1)
        write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
        CALL DoStop

    ELSE
    ! first find the AIRS layer within which it lies
        iFound=-1
        i1 = iLowest
        i2 = iLowest+1
        10 CONTINUE
        IF ((rP <= raPressLevels(i1)) .AND. (rP > raPressLevels(i2))) THEN
            iFound=1
        END IF
        IF ((iFound < 0) .AND. (i1 < kProfLayer)) THEN
            i1=i1+1
            i2=i2+1
	    GOTO 10
        END IF
        IF ((iFound < 0)) THEN
            IF (abs(rP-raPressLevels(kProfLayer+1)) <= delta) THEN
                i1=kProfLayer
                iFound=1
            ELSE
                write(kStdErr,*) 'could not find pressure ',rP
                write(kStdErr,*) 'within AIRS pressure levels. Please check'
                write(kStdErr,*) '*RADNCE and *OUTPUT sections'
                CALL DoSTOP
            END IF
        END IF
    END IF

    IF ((i1 > kProfLayer) .OR. (i1 < iLowest)) THEN
        write(kStdErr,*) 'sorry : cannot find surface temp for '
        write(kStdErr,*) 'layers outside ',iLowest,' and ',kProfLayer
        write(kStdErr,*) 'Allowed Pressure ranges are from : ', &
        raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
        write(kStdErr,*) 'Surface Pressure is ',rP,' mb'
        call DoStop
    END IF
              
! now find the temperature
    IF (i1 == iLowest) THEN          !do linear interp
        i1 = iLowest
        i2 = iLowest+1
        i3 = iLowest+2
        rP1 = (raPressLevels(i2)-raPressLevels(i1))/ &
        log(raPressLevels(i2)/raPressLevels(i1))
        rP2 = (raPressLevels(i3)-raPressLevels(i2))/ &
        log(raPressLevels(i3)/raPressLevels(i2))
        T1 = raProfileTemp(i1)
        T2 = raProfileTemp(i2)
        IF (iLog == -1) THEN
            rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)           !!linear in P
        ELSE
            rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1)) !!log(P)
        END IF

    ELSEIF (i1 >= (kProfLayer-1)) THEN          !do linear interp
        rP1 = (raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/ &
        log(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
        rP2 = (raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/ &
        log(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
        T1 = raProfileTemp(kProfLayer-1)
        T2 = raProfileTemp(kProfLayer)
        IF (iLog == -1) THEN
            rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)            !!linear in P
        ELSE
            rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1))  !!log(P)
        END IF
	
    ELSE
         !do spline ... note that the pressures have to
	 !be in ascENDing order for good interpolation
			     
        rP1 = (raPressLevels(i1)-raPressLevels(i1-1))/ &
        log(raPressLevels(i1)/raPressLevels(i1-1))
        raP(3) = rP1
        rP1 = (raPressLevels(i1+1)-raPressLevels(i1))/ &
        log(raPressLevels(i1+1)/raPressLevels(i1))
        raP(2) = rP1
        rP1 = (raPressLevels(i1+2)-raPressLevels(i1+1))/ &
        log(raPressLevels(i1+2)/raPressLevels(i1+1))
        raP(1) = rP1
        IF (iLog == +1) THEN
            DO iJ = 1,3
                raLogP(iJ) = log(raP(iJ))
            END DO
        END IF

        raT(3) = raProfileTemp(i1-1)
        raT(2) = raProfileTemp(i1)
        raT(1) = raProfileTemp(i1+1)

        yp1=1.0e30
        ypn=1.0e30
        IF (iSpline == +1) THEN
            IF (iLog == +1) THEN
                CALL rspl1(raLogP,raT,3,log(rP),rT,1)
            ELSE
                CALL rspl1(raP,raT,3,rP,rT,1)
            END IF
        ELSEIF (iSpline == -1) THEN
            IF (iLog == +1) THEN
                CALL rlinear1(raP,raT,3,rP,rT,1)
            ELSE
                CALL rlinear1(raLogP,raT,3,log(rP),rT,1)
            END IF
        END IF
    END IF

    FindBottomTemp = rT

    RETURN
    end FUNCTION FindBottomTemp

!************************************************************************
! if param kSurfTemp = +1, this computes surface temp by interpolating across
! pressure layers, and adds on offet given by rTSurf (which is the usual
! parameter normally used for surface temp in *RADNCE)
! else if kSurfTemp = -1, it just returns the user specified temperature
    REAL FUNCTION FindSurfaceTemp(rPressStart,rPressStop, &
    rTSurf,raProfileTemp, &
    raPresslevels,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    REAL :: rPressStart,rPressStop,rTSurf,raProfileTemp(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers

! local variables
    REAL :: rT
    INTEGER :: iI

    rT = rTSurf   !! this is the temp set in nm_radnce; logic below determines if it is offset or actual stemp
          
! this ORIGINAL code, allowed user to add in an offset
    IF ((kSurfTemp > 0) .AND. ((kRTP == -1) .OR. (kRTP == 0))) THEN
    ! ave to adjust temperature .. do this for down AND up look instr
        IF (rPressStart > rPressStop) THEN  ! for down looking instr
            rT = FindBottomTemp(rPressStart,raProfileTemp, &
            raPressLevels,iProfileLayers)
            rT = rT+rTSurf
        ELSEIF (rPressStart < rPressStop) THEN  ! for up looking instr
            rT = FindBottomTemp(rPressStop,raProfileTemp, &
            raPressLevels,iProfileLayers)
            rT = rT+rTSurf
        END IF
    END IF
!      ELSEIF ((kSurfTemp .gt. 0) .AND. ((kRTP .EQ. -5) .OR. (kRTP .EQ. -6))) THEN
!        !just state this has already been taken care of in subr radnce4
!       write(kStdWarn,*) 'kSurfTemp > 0 and kRTP = -5 or -6'
!       write(kStdWarn,*) 'so we already added in raTSurf offset to stemp from TAPE5/6'
!        rT = rT
!        END IF
!      END IF

! this was the code in Dec 2016, get rid of it
! this was the code in Dec 2016, get rid of it
! this was the code in Dec 2016, get rid of it
! why make life complicated, just directly give USER defined STEMP
!      IF ((kSurfTemp .gt. 0) .AND. ((kRTP .GE. -6) .AND. (kRTP .LE. 0))) THEN
!        rT = rTSurf
!      END IF
! replace with this why make life complicated, just directly give USER defined STEMP
    IF ((kSurfTemp > 0) .AND. ((kRTP >= -6) .AND. (kRTP <= -5))) THEN
        rT = rTSurf
    END IF

    FindSurfaceTemp = rT
            
    IF (rT < 190.0) THEN
        write(kStdErr,*)'Surface Temperature = ',rT-273,' deg C (',rT,' K)'
        write(kStdErr,*)'brrrrrrrrrrrrrrrrrrrrrrrrrr!!!!!!!'
        write(kStdErr,*)'kCARTA allows surface temps between 190 and 350K'
        CALL DoSTOP
    END IF

    IF (rT > 350.0) THEN
        write(kStdErr,*)'Surface Temperature = ',rT-273,' deg C (',rT,' K)'
        write(kStdErr,*)'whew!!!!! bloody hot!!!!!!!'
        write(kStdErr,*)'kCARTA allows temps between 210 and 350K'
        CALL DoSTOP
    END IF
            
    RETURN
    end FUNCTION FindSurfaceTemp

!************************************************************************

    SUBROUTINE INTERP_SCAT_TABLE2 (WAVENO, DME, &
    EXTINCT, SSALB, ASYM, &
    NDME, DMETAB, NWAVE, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM)
!       Interpolates the scattering properties from the table for
!     a particular wavenumber and particle size.  Does a bilinear
!     interpolation, but optimized for fixed particle size and slowly
!     varying wavenumber.  If the DME is the same as last time then we
!     can just linearly interpolate in wavenumber between stored
!     scattering values.  If the DME has changed then we linearly
!     interpolate between the DMETAB grid lines for each of the two
!     wavenumber grid lines.

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL ::     WAVENO, DME
    REAL ::     EXTINCT, SSALB, ASYM
    INTEGER ::  NDME, NWAVE
    REAL ::     DMETAB(NDME), WAVETAB(NWAVE)
    REAL ::     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME)
    REAL ::     TABASYM(NWAVE,NDME)
    INTEGER ::  IW0, IW1, ID, IL, IU, IM
    LOGICAL ::  NEWIW
    REAL ::     FWAV, FDME, FLDME, F
    REAL ::     OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
! sergio do not save iw0,iw1, olddme
!      SAVE     IW0, IW1, ID, OLDDME, FDME, FLDME
!      SAVE     ID, FDME, FLDME
!      SAVE     EXT0,EXT1, ALB0,ALB1, ASYM0,ASYM1
    DATA     IW0/1/, IW1/2/

    INTEGER :: iLogORLinear,iDefault

    iw0 = 1
    iw1 = 2
    olddme = 0.0

    iw0 = 1
    iw1 = nwave
    olddme = -10.0

    iDefault = -1       !! do linear for w,g and log for e
    iLogOrLinear = +1    !! do linear for w,g and log for e; default RTSPEC
    iLogOrLinear = -1    !! do log for w,g    and log for e; default SARTA
    iLogOrLinear = iaaOverrideDefault(1,7)
    IF (abs(iLogOrLinear) /= 1) THEN
        write(kStdErr,*) 'invalid iLogOrLinear = ',iLogOrLinear
        CALL DoStop
    END IF
    IF ((iDefault /= iLogOrLinear) .AND. (kOuterLoop == 1)) THEN
        write (kStdErr,*) 'in INTERP_SCAT_TABLE2'
        write (kStdErr,*)  'iDefault,iLogOrLinear = ',iDefault,iLogOrLinear
    END IF

!         Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
        write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
        IF (WAVENO < WAVETAB(1)) THEN
            WAVENO = WAVETAB(1)
        ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
            WAVENO = WAVETAB(NWAVE)
        END IF
    ! ALL DoStop
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
    !        write(kStdErr,*) DME,' outside ',DMETAB(1),':',DMETAB(NDME)
    !        write(kStdErr,*) 'INTERP_SCAT_TABLE: particle Dme out of range ... RESET'
        IF (DME < DMETAB(1)) THEN
            DME = DMETAB(1)
        ELSEIF (DME > DMETAB(NDME)) THEN
            DME = DMETAB(NDME)
        END IF
    ! ALL DoStop
    END IF

!         See if wavenumber is within last wavenumber grid, otherwise
!           find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
!      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (WAVENO >= WAVETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE. 
    ENDIF

    IF (DME /= OLDDME) THEN
    !         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (DME >= DMETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
    ENDIF

    IF ((DME /= OLDDME .OR. NEWIW) .AND. (iLogOrLinear == +1)) THEN
    !         If not the same Dme or a new wavenumber grid, then
    !           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)

    ! looking at sarta code, Scott Hannon ALWAYS does a log interp
    ELSEIF ((DME /= OLDDME .OR. NEWIW) .AND. (iLogOrLinear == -1)) THEN
    !         If not the same Dme or a new wavenumber grid, then
    !           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = EXP( (1-FLDME)*LOG(TABSSALB(IW0,ID)) &
        + FLDME*LOG(TABSSALB(IW0,ID+1)) )
        ALB1 = EXP( (1-FLDME)*LOG(TABSSALB(IW1,ID)) &
        + FLDME*LOG(TABSSALB(IW1,ID+1)) )
        ASYM0 = EXP( (1-FLDME)*LOG(TABASYM(IW0,ID)) &
        + FLDME*LOG(TABASYM(IW0,ID+1)) )
        ASYM1 = EXP( (1-FLDME)*LOG(TABASYM(IW1,ID)) &
        + FLDME*LOG(TABASYM(IW1,ID+1)) )
    ENDIF

!         Linearly interpolate the scattering properties in wavenumber
    FWAV    = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
    F       = 1-FWAV
    EXTINCT = F*EXT0 + FWAV*EXT1
    SSALB   = F*ALB0 + FWAV*ALB1
    ASYM    = F*ASYM0 + FWAV*ASYM1

    OLDDME = DME

    RETURN
    END SUBROUTINE INTERP_SCAT_TABLE2

!************************************************************************

! this file reads a binary made from the ASCII sscatmie.x file
    SUBROUTINE READ_SSCATTAB_BINARY(SCATFILE,    & !!!  MAXTAB, MAXGRID,
    cScale, NMUOBS, MUTAB, NDME, DMETAB, NWAVE, WAVETAB, &
    MUINC, TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

!       Reads in the single scattering table for a number of wavenumbers,
!     particle sizes, and viewing angles.  The scattering properties are
!     computed for a IWC/LWC of 1 g/m^3.
!       Input parameters:
!     SCATFILE   file name of scattering file
!     MAXTAB     maximum array size for whole table
!     MAXGRID    maximum array size for table grids
!       Output parameters:
!     NMUOBS     number of viewing angle mu grid values
!     MUTAB      viewing angle grid values
!     NDME       number of Dme grid values
!     DMETAB     Dme grid values
!     NWAVE      number of wavenumber grid values
!     WAVETAB    wavenumber grid values
!     MUINC(2)   cosine zenith of two incident angles
!     TABEXTINCT tabulated extinction (km^-1)
!     TABSSALB   tabulated single scattering albedo
!     TABASYM    tabulated asymmetry parameter
!     TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN
!                tabulated phase function terms for incident radiance angles

!!!!      INTEGER  MAXTAB, MAXGRID
    INTEGER ::  NMUOBS, NDME, NWAVE
    REAL ::     MUTAB(*), DMETAB(*), WAVETAB(*)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(*), TABSSALB(*), TABASYM(*)
    REAL ::     TABPHI1UP(*), TABPHI1DN(*)
    REAL ::     TABPHI2UP(*), TABPHI2DN(*)
    CHARACTER*(*) SCATFILE
    INTEGER ::  IMU, ID, IW, K2, K3
    INTEGER :: IERR
    CHARACTER cScale

    OPEN (UNIT = kTempUnit, STATUS='OLD', FORM='UNFORMATTED', &
    FILE=SCATFILE, IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,1010) IERR, SCATFILE
        CALL DoSTOP
    ENDIF
    1010 FORMAT('ERROR! number ',I5,' opening scatter data file:',/,A120)

    kTempUnitOpen=1
    READ(kTempUnit) NMUOBS
    READ(kTempUnit) NDME
    READ(kTempUnit) NWAVE

    IF (MAX(NMUOBS,NDME,NWAVE) > MAXGRID) THEN
        write(kStdErr,*) '(MAX(NMUOBS,NDME,NWAVE) > MAXGRID) '
        write(kStdErr,*)  NMUOBS,NDME,NWAVE,MAXGRID
        write(kStdErr,*) 'READ_SSCATTAB_BINARY: MAXGRID exceeded'
        CALL DoStop
    END IF
    IF (NMUOBS*NDME*NWAVE > MAXTAB) THEN
        write(kStdErr,*) '(NMUOBS*NDME*NWAVE > MAXTAB)'
        write(kStdErr,*)  NMUOBS,NDME,NWAVE,MAXTAB
        write(kStdErr,*) 'READ_SSCATTAB_BINARY: MAXTAB exceeded'
        CALL DoStop
    END IF

    READ(kTempUnit) MUINC(1), MUINC(2)
    READ(kTempUnit) cScale
    READ(kTempUnit) (MUTAB(IMU), IMU = 1, NMUOBS)

!      print *,NMUOBS,NDME,NWAVE
!      print *,MUINC(1),MUINC(2)
!      print *,cScale
!      print *,(MUTAB(IMU), IMU = 1, NMUOBS)

    DO IW = 1, NWAVE
        DO ID = 1, NDME
            K2 = IW-1 + NWAVE*(ID-1)
            K3 = NMUOBS*K2
        !          print *,IW,ID,K2,K3
            READ(kTempUnit) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1), &
            TABSSALB(K2+1), TABASYM(K2+1)
        !          print *,K2,DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1),
        !     $       TABSSALB(K2+1), TABASYM(K2+1)
            READ(kTempUnit) (TABPHI1UP(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI2UP(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI1DN(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI2DN(IMU+K3), IMU = 1, NMUOBS)
        ENDDO
    ENDDO
          
    CLOSE (kTempUnit)
    kTempUnitOpen=-1

    write(kStdWarn,*)'success : read in binary scattr data from file = '
    write(kStdWarn,1020) scatfile

!      call dostop
          
    1020 FORMAT(A70)

    RETURN
    END SUBROUTINE READ_SSCATTAB_BINARY

! ************************************************************************
    SUBROUTINE READ_SSCATTAB(SCATFILE,                & !!!  MAXTAB, MAXGRID,
    cScale, NMUOBS, MUTAB, NDME, DMETAB, NWAVE, WAVETAB, &
    MUINC, TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

!       Reads in the single scattering table for a number of wavenumbers,
!     particle sizes, and viewing angles.  The scattering properties are
!     computed for a IWC/LWC of 1 g/m^3.
!       Input parameters:
!     SCATFILE   file name of scattering file
!     MAXTAB     maximum array size for whole table
!     MAXGRID    maximum array size for table grids
!       Output parameters:
!     cScale     Scaling (n,y,h,g) ... needed for DISORT
!     NMUOBS     number of viewing angle mu grid values
!     MUTAB      viewing angle grid values
!     NDME       number of Dme grid values
!     DMETAB     Dme grid values
!     NWAVE      number of wavenumber grid values
!     WAVETAB    wavenumber grid values
!     MUINC(2)   cosine zenith of two incident angles
!     TABEXTINCT tabulated extinction (km^-1)
!     TABSSALB   tabulated single scattering albedo
!     TABASYM    tabulated asymmetry parameter
!     TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN
!                tabulated phase function terms for incident radiance angles

!!!!      INTEGER  MAXTAB, MAXGRID
    INTEGER ::  NMUOBS, NDME, NWAVE
    REAL ::     MUTAB(*), DMETAB(*), WAVETAB(*)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(*), TABSSALB(*), TABASYM(*)
    REAL ::     TABPHI1UP(*), TABPHI1DN(*)
    REAL ::     TABPHI2UP(*), TABPHI2DN(*)
    CHARACTER*(*) SCATFILE
    INTEGER ::  IMU, ID, IW, K2, K3

    CHARACTER(80) :: caLine
    CHARACTER cScale

    OPEN (UNIT=2, STATUS='OLD', FILE=SCATFILE)
    READ (2,*)
    READ (2,*) NMUOBS
    READ (2,*) NDME
    READ (2,*) NWAVE
    IF (MAX(NMUOBS,NDME,NWAVE) > MAXGRID) &
    STOP 'READ_SSCATTAB: MAXGRID exceeded'
    IF (NMUOBS*NDME*NWAVE > MAXTAB) &
    STOP 'READ_SSCATTAB: MAXTAB exceeded'
    READ (2,*) MUINC(1), MUINC(2)
    READ (2,*)
    READ (2,*)
    READ (2,30) caLine
    READ (2,*)
    READ (2,*)
    READ (2,*) (MUTAB(IMU), IMU = 1, NMUOBS)
    DO IW = 1, NWAVE
        DO ID = 1, NDME
            K2 = IW-1 + NWAVE*(ID-1)
            K3 = NMUOBS*K2
            READ(2,*)
            READ(2,*)
            READ(2,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1), &
            TABSSALB(K2+1), TABASYM(K2+1)
            READ(2,*) (TABPHI1UP(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI2UP(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI1DN(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI2DN(IMU+K3), IMU = 1, NMUOBS)
        ENDDO
    ENDDO

    30 FORMAT(A80)

    CLOSE (UNIT=2)
     
    CALL FindScalingParameter(caLine,cScale)

    RETURN
    END SUBROUTINE READ_SSCATTAB

!************************************************************************
! this is for reading in stuff from eg Baran's files
    SUBROUTINE READ_SSCATTAB_SPECIAL(SCATFILE, &
    cScale, NMUOBS, MUTAB, NDME, DMETAB, NWAVE, WAVETAB, &
    MUINC, TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

!       Reads in the single scattering table for a number of wavenumbers,
!     particle sizes, and viewing angles.  The scattering properties are
!     computed for a IWC/LWC of 1 g/m^3.
!       Input parameters:
!     SCATFILE   file name of scattering file
!     MAXTAB     maximum array size for whole table
!     MAXGRID    maximum array size for table grids
!       Output parameters:
!     cScale     Scaling (n,y,h,g) ... needed for DISORT
!     NMUOBS     number of viewing angle mu grid values
!     MUTAB      viewing angle grid values
!     NDME       number of Dme grid values
!     DMETAB     Dme grid values
!     NWAVE      number of wavenumber grid values
!     WAVETAB    wavenumber grid values
!     MUINC(2)   cosine zenith of two incident angles
!     TABEXTINCT tabulated extinction (km^-1)
!     TABSSALB   tabulated single scattering albedo
!     TABASYM    tabulated asymmetry parameter
!********* these are dummy!
!     TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN
!                tabulated phase function terms for incident radiance angles
!********* these are dummy!

    INTEGER ::  NMUOBS, NDME, NWAVE
    REAL ::     MUTAB(*), DMETAB(*), WAVETAB(*)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(*), TABSSALB(*), TABASYM(*)
    REAL ::     TABPHI1UP(*), TABPHI1DN(*)
    REAL ::     TABPHI2UP(*), TABPHI2DN(*)
    CHARACTER*(*) SCATFILE
    INTEGER ::  IMU, ID, IW, K2, K3

    CHARACTER(80) :: caLine
    CHARACTER cScale

    NMUOBS = -9999

    OPEN (UNIT=2, STATUS='OLD', FILE=SCATFILE)
    READ (2,*)
    READ (2,*) NDME
    READ (2,*) NWAVE
    IF (MAX(NMUOBS,NDME,NWAVE) > MAXGRID) &
    STOP 'READ_SSCATTAB_SPECIAL: MAXGRID exceeded'
    IF (NMUOBS*NDME*NWAVE > MAXTAB) &
    STOP 'READ_SSCATTAB_SPECIAL: MAXTAB exceeded'
    READ (2,*)
    READ (2,*)
    READ (2,*)
    READ (2,*)

    DO IW = 1, NWAVE
        DO ID = 1, NDME
            K2 = IW-1 + NWAVE*(ID-1)
            READ(2,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1), &
            TABSSALB(K2+1), TABASYM(K2+1)
            TABEXTINCT(K2+1) = TABEXTINCT(K2+1) * 1000.0   !!!to be like sscatmie
        !          write(kStdWarn,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1),
        !     $       TABSSALB(K2+1), TABASYM(K2+1)
        ENDDO
    ENDDO

    30 FORMAT(A80)

    CLOSE (UNIT=2)

    cScale =  'H'

    RETURN
    END SUBROUTINE READ_SSCATTAB_SPECIAL

!************************************************************************
! this subroutine parses the line and finds first nonzero character
    SUBROUTINE FindScalingParameter(caLine,cScale)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
      
    CHARACTER(80) :: caLine
    CHARACTER cScale
      
    INTEGER :: iI,iFound
      
    iFound = -1
    iI = 1
    10 CONTINUE
    IF (caLine(iI:iI) == ' ') THEN
        iI = iI + 1
    ELSE
        iFound = 1
        cScale = caLine(iI:iI)
    END IF
    IF ((iI <= 80) .AND. (iFound < 0)) THEN
        GOTO 10
    END IF
            
    IF ((caLine(iI:iI) /= 'N') .AND. (caLine(iI:iI) /= 'Y') .AND. &
    (caLine(iI:iI) /= 'G') .AND. (caLine(iI:iI) /= 'H')) THEN
        iFound = -1
        iI = iI + 1
        GOTO 10
    END IF
            
    IF ((iI == 80) .AND. (iFound < 0)) THEN
        write (kStdErr,*) 'never found scaling parameter (n,y,g,h)!!!!'
        CALL DoSTOP
    ELSE
        write (kStdWarn,*) 'scaling parameter in Mie Tables = ',cScale
    END IF
     
    RETURN
    end SUBROUTINE FindScalingParameter
    
!************************************************************************
! this subroutine very quickly does the radiative transfer
! since the optical depths are soooooooooo small, use double precision
    SUBROUTINE UpperAtmRadTrans(raInten,raFreq,rSatAngle, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
!   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!   upper atm abs coeff and planck coeff
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!   input wavevector and integer stating which layer to stop rad transfer at
    REAL :: raFreq(kMaxPts),rSatAngle
    INTEGER :: iUpper
! do we want to do upper atm NLTE computations?
    INTEGER :: iDoUpperAtmNLTE
! do we dump all or some rads?
    INTEGER :: iDumpAllUARads
! input/output pararameter
!   this contains the radiance incident at kCARTA TOA (0.005 mb)
!   it will finally contain the radiance exiting from TOP of UPPER ATM
    REAL :: raInten(kMaxPts)

! local variables
    INTEGER :: iFr,iL,iIOUN
    REAL :: rEmission,rTrans,rMu,raInten0(kMaxPts)
    DOUBLE PRECISION :: daInten(kMaxPts),dTrans,dEmission
    CHARACTER(80) :: caOutName

    caOutName = 'DumDum'
    iIOUN = kNLTEOutUA
      
    IF (iDoUpperAtmNLTE <= 0) THEN
        write (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
        CALL DoStop
    ELSE
        write(kStdWarn,*) 'Doing UA (NLTE) radtransfer at 0.0025 cm-1 '
    END IF

! compute radiance intensity thru NEW uppermost layers of atm
    DO iFr = 1,kMaxPts
        raInten0(iFr) = raInten(iFr)
        daInten(iFr)  = dble(raInten(iFr))
    END DO

    iL = 0
    IF (kNLTEOutUAOpen > 0) THEN
        write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
    ! always dump out the 0.005 mb TOA radiance if the UA file is open
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF

    rMu = cos(rSatAngle*kPi/180.0)

    DO iL = 1,iUpper - 1

        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = (raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = (raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            (ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans

            raInten(iFr) = sngl(daInten(iFr))
        END DO

        IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
            write(kStdWarn,*) 'dumping out UA rads at iL = ',iL
        ! dump out the radiance at this HIGH pressure level
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    DO iL = iUpper,iUpper
        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = dble(raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = dble(raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            dble(ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans
            raInten(iFr) = sngl(daInten(iFr))

        END DO

        IF (kNLTEOutUAOpen > 0) THEN
        ! always dump out the 0.000025 mb TOA radiance if the UA file is open
            write(kStdWarn,*) 'dumping out 0.000025 mb UA rads iL = ',iL
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE UpperAtmRadTrans

!************************************************************************

END MODULE clear_scatter_basic
