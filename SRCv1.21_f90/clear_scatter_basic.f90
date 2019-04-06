! Copyright 2006
! University of Maryland Baltimore County
! All Rights Reserved

MODULE clear_scatter_basic

USE basic_common
USE ttorad_common
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
    INTEGER :: iRTSPEC
    REAL :: raBeta(kMaxPts),raTT(kMaxPts),raZeta(kMaxPts),rBooga
    REAL :: raPlanck1(kMaxPts),raPlanck0(kMaxPts),raDel(kMaxPts),raTau0(kMaxPts),raTrans(kMaxPts),gwak

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
      raInten = raInten*exp(-raaAbs(:,iL)/rCos)

    ELSEIF ((iVary > -1) .AND. (iRTSPEC < 0)) THEN
      !!!either exp temperature dependance or none; rBooga carries this info
      !!! >>>>>>>>>>>>>>> this is basically exp in tau <<<<<<<<<<<<<<<<<<<<<<
      IF (rFrac >= 0.9999) THEN
        raBeta = 1/raaAbs(:,iL) * rBooga
        raTT   = ttorad(raFreq,TEMP(iBeta))/(1 + raBeta*rCos)
        raZeta = (raInten - raTT) * exp(-raaAbs(:,iL)/rCos)
        raInten = raZeta + raTT * exp(raaAbs(:,iL) * raBeta)
      ELSE
        raBeta = 1/(raaAbs(:,iL)*rFrac) * rBooga
        raTT   = ttorad(raFreq,TEMP(iBeta))/(1 + raBeta*rCos)
        raZeta = (raInten - raTT) * exp(-raaAbs(:,iL)*rFrac/rCos)
        raInten = raZeta + raTT * exp(raaAbs(:,iL)*rFrac * raBeta)
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

      raPlanck1 = ttorad(raFreq,TEMP(iBeta))
      raPlanck0 = ttorad(raFreq,TEMP(iBetaP1))
      raTau0 = (raaAbs(:,iL)*gwak)/rCos
      WHERE (ratau0 < 0.001)
        raInten = raInten*(1-raTau0) + ratau0*0.5*(raPLANCK0+raPLANCK1)
      ELSEWHERE (raTau0 >= 0.001)
        raDel = (raPlanck1-raPlanck0)/raTau0
        raTrans = exp(-raTau0)
        raInten = raInten*raTrans + (raPlanck0 + raDel - raTrans*(raPlanck0+raDel*(1.0+raTau0)))
      ENDWHERE

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
    REAL :: mu
    INTEGER :: iRTSPEC
    REAL :: raBeta(kMaxPts),raTT(kMaxPts),raZeta(kMaxPts),rBooga
    REAL :: raPlanck1(kMaxPts),raPlanck0(kMaxPts),raDel(kMaxPts),raTau0(kMaxPts),raTrans(kMaxPts),gwak

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
      raInten = raInten*exp(-raaAbs(:,iL)/rCos)
    ELSEIF ((iVary >= -1) .AND. (iRTSPEC < 0)) THEN
      !!!either exp temperature dependace or none; rBooga carries this info
      IF (rFrac >= 0.9999) THEN
        raBeta = 1/raaAbs(:,iL) * rBooga
        raTT   = ttorad(raFreq,TEMP(iBeta))/(raBeta*mu - 1)
        raZeta = exp(-raaAbs(:,iL)/mu) * exp(raBeta*raaAbs(:,iL)) - 1.0
        raInten = raInten*exp(-raaAbs(:,iL)/mu) + raTT*raZeta
      ELSE
        raBeta = 1/(raaAbs(:,iL)*rFrac) * rBooga
        raTT   = ttorad(raFreq,TEMP(iBeta))/(raBeta*mu - 1)
        raZeta = exp(-raaAbs(:,iL)*rFrac/mu) * exp(raBeta*raaAbs(:,iL)*rFrac) - 1.0
        raInten = raInten*exp(-raaAbs(:,iL)*rFrac/mu) + raTT*raZeta
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

      IF (rFrac >= 1.0000) gwak = 1.0
      IF (rFrac < 0.9999) gwak = rFrac

      raPlanck1 = ttorad(raFreq,TEMP(iBeta))
      raPlanck0 = ttorad(raFreq,TEMP(iBetaM1))
      raTau0 = (raaAbs(:,iL)*gwak)/rCos
      WHERE (ratau0 < 0.001)
        raInten = raInten*(1-raTau0) + ratau0*0.5*(raPLANCK0+raPLANCK1)
      ELSEWHERE (raTau0 >= 0.001)
        raDel = (raPlanck1-raPlanck0)/raTau0
        raTrans = exp(-raTau0)
        raInten = raInten*raTrans + (raPlanck1 - raDel - raTrans*(raPlanck1-raDel*(1.0+raTau0)))
      ENDWHERE
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
    REAL :: rBeff,raFcn(kMaxPts)
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: raZeta(kMaxPts),raZeta2(kMaxPts),raAbs(kMaxPts),raTrans(kMaxPts),raCos(kMaxPts)

    raCos = rCos

    IF (iVary < 2) THEN
      write(kStdErr,*) 'this is upwell for linear in tau .. need iVary = 2 or 3 or 4'
      CALL DoStop
    END IF

    IF (rFrac < 0) THEN
      write(kStdErr,*) 'Warning rFrac < 0 in RT_ProfileUPWELL_LINTAU, reset to > 0'
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

    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level
    CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1)   !! ttorad of upper level  XXXXX this we want XXXXX
    CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)    !! ttorad of Tlayer
    !! (which is NOT necessarily average of above 2)
    IF (kOuterLoop == 1) THEN
      write(kStdWarn,2345) iL,TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
    END IF
          
 1234 FORMAT(I3,3(' ',F10.3))
 2345 FORMAT('up [iLDN=iL iLay iLUP=iLp1]',I3,3(' ',F10.3))
     
    IF (iVary == 2) THEN
      !!! lim tau --> 0 , rFcn --> 0
      CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP0)
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = (raIntenP1 - raIntenP0 + 1.0e-10)/(raAbs + 1.0e-10)
      raInten = raInten * exp(-raAbs/raCos) + raIntenP0 * (1 - exp(-raAbs/raCos))
      WHERE (raAbs >= 0.001) 
        raInten = raInten + raFcn*raCos*(raAbs/raCos-1) + raFcn*raCos*exp(-raAbs/raCos)
      END WHERE

    ELSEIF (iVary == +3) THEN
      !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = 1.0
      WHERE (raAbs >= 0.001)
        raFcn = exp(-raAbs/rCos)
        raFcn = raCos/raAbs - raFcn/(1-raFcn)
      END WHERE
      raFcn = raIntenP1 + 2*(raIntenAvg-raIntenP1)*raFcn
      raInten = raInten * exp(-raAbs/raCos) + raFcn * (1 - exp(-raAbs/raCos))

    ELSEIF (iVary == +40) THEN
      !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = 1.0
      raTrans = 1.0
      WHERE (raAbs >= 0.0001)
        raTrans = exp(-raAbs/raCos)
        raFcn = raCos/raAbs * (1 - raTrans)
      END WHERE
      raZeta = raIntenP1*(1-raTrans) + (raIntenP - raIntenP1)*(raFcn - raTrans)
      raInten = raInten * exp(-raAbs/raCos) + raZeta

    ELSEIF (iVary == +41) THEN
      !!! print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! PADE APPROX, two term (combo of GENLN2 and LBLRTM)
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raTrans = exp(-raAbs)
      raZeta = 0.2*raAbs             !! pade one
      raFcn = (raIntenAvg + raZeta*raIntenP1)/(1+raZeta)
      raZeta = 0.193*raAbs           !! pade two
      raZeta2 = 0.013*raAbs*raAbs    !! pade two
      raFcn = (raIntenAvg + (raZeta + raZeta2)*raIntenP1)/(1+raZeta+raZeta2)
      raFcn = (1-raTrans)*raFcn
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +42) THEN
      !!!      print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
      !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, GENLN2 style
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP1)
      WHERE (raAbs >= 0.05) 
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP1 + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs
        raFcn = raAbs*raIntenP1 + raZeta*(1-raAbs/2) - raTrans * raZeta
      END WHERE
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +43) THEN
      !!!         http://www.wolframalpha.com/input/?i=1-2*%281%2Fx-exp%28-x%29%2F%281-exp%28-x%29%29%29
      !!!        print *,'up flux ',iL,iBeta,rFrac,raaAbs(1,iL),TEMPLEV(iBeta),TEMPLAY(iBeta),TEMPLEV(iBeta+1)
      !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = raIntenP1 - raIntenAvg
      WHERE (raAbs >= 0.06) 
        raTrans = exp(-raAbs)
        raZeta2 = 1.0 - 2.0*(1/raAbs - raTrans/(1-raTrans))
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*(raAbs * raAbs)
        raZeta2 = raAbs/6.0 - (raAbs**3)/360.0 + (raAbs**5)/15120.0   !! mathematica
        raZeta2 = raAbs/6.0
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      END WHERE
      raInten = raInten*raTrans + raFcn

    !y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')

    ELSEIF (iVary == +4) THEN
      !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, MY style
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP1)
      WHERE (raAbs > 0.1) 
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP1 + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*raAbs**2
        raZeta2 = raZeta*(raAbs/2-(raAbs**2)/3+(raAbs**3)/6)
        raFcn   = (1-raTrans)*raIntenP1 + raZeta2
      END WHERE
      raInten = raInten*raTrans + raFcn

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
    REAL :: rBeff,raFcn(kMaxPts)
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: raZeta(kMaxPts),raZeta2(kMaxPts),raAbs(kMaxPts),raTrans(kMaxPts)

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
      CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)      !! ttorad of lower level  XXXX what we want XXXXX
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
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = (raIntenP1 - raIntenP0 + 1.0e-10)/(raAbs + 1.0e-10)
      raInten = raInten * exp(-raAbs/rCos) + raIntenP0 * (1 - exp(-raAbs/rCos))
      WHERE (raAbs >= 0.001)
          raInten = raInten + raFcn*rCos*(raAbs/rCos-1) + raFcn*rCos*exp(-raAbs/rCos)
      END WHERE

    ELSEIF (iVary == +3) THEN
      !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL) * rFrac
      END IF
      raFcn = 1.0
      WHERE (raAbs >= 0.001) 
        raFcn = exp(-raAbs/rCos)
        raFcn = rCos/raAbs - raFcn/(1-raFcn)
      END WHERE
      raFcn = raIntenP1 + 2*(raIntenAvg-raIntenP1)*raFcn
      raInten = raInten * exp(-raAbs/rCos) + raFcn * (1 - exp(-raAbs/rCos))

    ELSEIF (iVary == +40) THEN
      !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = 1.0
      raTrans = 1.0
      WHERE (raAbs >= 0.0001) 
        raTrans = exp(-raAbs/rCos)
        raFcn = rCos/raAbs * (1 - raTrans)
      END WHERE
      raZeta = raIntenP*(1-raTrans) + (raIntenP1 - raIntenP)*(raFcn - raTrans)
      raInten = raInten * exp(-raAbs/rCos) + raZeta

    ELSEIF (iVary == +41) THEN
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raTrans = exp(-raAbs)
      raZeta = 0.2*raAbs             !! pade one
      raFcn = (raIntenAvg + raZeta*raIntenP)/(1+raZeta)
      raZeta = 0.193*raAbs           !! pade two
      raZeta2 = 0.013*raAbs*raAbs    !! pade two
      raFcn = (raIntenAvg + (raZeta + raZeta2)*raIntenP)/(1+raZeta+raZeta2)
      raFcn = (1-raTrans)*raFcn
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +42) THEN
      !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, GENLN2 style
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP)
      WHERE (raAbs >= 0.05) 
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs
        raFcn = raAbs*raIntenP + raZeta*(1-raAbs/2) - raTrans * raZeta
      END WHERE
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +43) THEN
      !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = raIntenP - raIntenAvg
      WHERE (raAbs >= 0.06)
        raTrans = exp(-raAbs)
        raZeta2 = 1.0 - 2.0*(1/raAbs - raTrans/(1-raTrans))
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*(raAbs * raAbs)
        raZeta2 = raAbs/6.0 - (raAbs**3)/360.0 + (raAbs**5)/15120.0  !! mathematica
        raZeta2 = raAbs/6.0
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      END WHERE
      raInten = raInten*raTrans + raFcn
              
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')

    ELSEIF (iVary == +4) THEN
      !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, MY style
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = 2*(raIntenAvg - raIntenP)
      WHERE (raAbs > 0.1)
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*raAbs**2
        raZeta2 = raZeta*(raAbs/2-(raAbs**2)/3+(raAbs**3)/6)
        raFcn   = (1-raTrans)*raIntenP + raZeta2
      END WHERE
      raInten = raInten*raTrans + raFcn

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
    REAL :: rBeff,raFcn(kMaxPts)
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: raZeta(kMaxPts),raZeta2(kMaxPts),raAbs(kMaxPts),raTrans(kMaxPts)

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
                                                                 
    IF (iVary >= 4) THEN
      !! new option
      CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta),raIntenP)    !! ttorad of lower level XXXX what we want XXXXXXXX
      CALL ttorad_oneBT2array(raFreq,TEMPLEV(iBeta+1),raIntenP1) !! ttorad of upper level
      CALL ttorad_oneBT2array(raFreq,TEMPLAY(iBeta),raIntenAvg)  !! ttorad of Tlayer
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
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = (raIntenP1 - raIntenP0 + 1.0e-10)/(raAbs + 1.0e-10)
      raInten = raInten * exp(-raAbs/rCos) + raIntenP0 * (1 - exp(-raAbs/rCos))
      WHERE (raAbs >= 0.001)
          raInten = raInten + raFcn*rCos*(raAbs/rCos-1) + raFcn*rCos*exp(-raAbs/rCos)
      END WHERE

    ELSEIF (iVary == +3) THEN
      !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL) * rFrac
      END IF
      raFcn = 1.0
      WHERE (raAbs >= 0.001) 
        raFcn = exp(-raAbs/rCos)
        raFcn = rCos/raAbs - raFcn/(1-raFcn)
      END WHERE
      raFcn = raIntenP1 + 2*(raIntenAvg-raIntenP1)*raFcn
      raInten = raInten * exp(-raAbs/rCos) + raFcn * (1 - exp(-raAbs/rCos))

    ELSEIF (iVary == +40) THEN
      !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = 1.0
      raTrans = 1.0
      WHERE (raAbs >= 0.0001) 
        raTrans = exp(-raAbs/rCos)
        raFcn = rCos/raAbs * (1 - raTrans)
      END WHERE
      raZeta = raIntenP*(1-raTrans) + (raIntenP1 - raIntenP)*(raFcn - raTrans)
      raInten = raInten * exp(-raAbs/rCos) + raZeta

    ELSEIF (iVary == +41) THEN
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raTrans = exp(-raAbs)
      raZeta = 0.2*raAbs             !! pade one
      raFcn = (raIntenAvg + raZeta*raIntenP)/(1+raZeta)
      raZeta = 0.193*raAbs           !! pade two
      raZeta2 = 0.013*raAbs*raAbs    !! pade two
      raFcn = (raIntenAvg + (raZeta + raZeta2)*raIntenP)/(1+raZeta+raZeta2)
      raFcn = (1-raTrans)*raFcn
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +42) THEN
      !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, GENLN2 style
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP)
      WHERE (raAbs >= 0.05) 
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs
        raFcn = raAbs*raIntenP + raZeta*(1-raAbs/2) - raTrans * raZeta
      END WHERE
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +43) THEN
      !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = raIntenP - raIntenAvg
      WHERE (raAbs >= 0.06)
        raTrans = exp(-raAbs)
        raZeta2 = 1.0 - 2.0*(1/raAbs - raTrans/(1-raTrans))
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*(raAbs * raAbs)
        raZeta2 = raAbs/6.0 - (raAbs**3)/360.0 + (raAbs**5)/15120.0  !! mathematica
        raZeta2 = raAbs/6.0
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      END WHERE
      raInten = raInten*raTrans + raFcn
              
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')

    ELSEIF (iVary == +4) THEN
      !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, MY style
      raAbs = raaAbs(:,iL)/rCos*rFrac
      raZeta = 2*(raIntenAvg - raIntenP)
      WHERE (raAbs > 0.1)
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*raAbs**2
        raZeta2 = raZeta*(raAbs/2-(raAbs**2)/3+(raAbs**3)/6)
        raFcn   = (1-raTrans)*raIntenP + raZeta2
      END WHERE
      raInten = raInten*raTrans + raFcn

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
    REAL :: rBeff,raFcn(kMaxPts)
    REAL :: raIntenP(kMaxPts),raIntenP1(kMaxPts),raIntenP0(kMaxPts)
    REAL :: raIntenAvg(kMaxPts)
    REAL :: raZeta(kMaxPts),raZeta2(kMaxPts),raAbs(kMaxPts),raTrans(kMaxPts)

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
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = (raIntenP1 - raIntenP0 + 1.0e-10)/(raAbs + 1.0e-10)
      raInten = raInten * exp(-raAbs/raCos) + raIntenP0 * (1 - exp(-raAbs/raCos))
      WHERE (raAbs >= 0.001)
        raInten = raInten + raFcn*raCos*(raAbs/raCos-1) + raFcn*raCos*exp(-raAbs/raCos)
      END WHERE

    ELSEIF (iVary == +3) THEN
      !!! this was done on June 24, 2013 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 13
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL)*rFrac
      END IF
      raFcn = 1.0
      WHERE (raAbs >= 0.001)
        raFcn = exp(-raAbs/raCos)
        raFcn = raCos/raAbs - raFcn/(1-raFcn)
      END WHERE
      raFcn = raIntenP1 + 2*(raIntenAvg - raIntenP1)*raFcn
      raInten = raInten * exp(-raAbs/raCos) + raFcn * (1 - exp(-raAbs/raCos))

    ELSEIF (iVary == +40) THEN
      !!! orig code uptil Oct 2015, buggy as it used raIntenP instead of raIntenAvg
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 9
      !!! lim tau --> 0 , rFcn --> 1
      IF (rFrac >= 0.9999) THEN
        raAbs = raaAbs(:,iL)
      ELSE
        raAbs = raaAbs(:,iL) * rFrac
      END IF
      raFcn = 1.0
      raTrans = 1.0
      WHERE (raAbs >= 0.0001)
        raTrans = exp(-raAbs/raCos)
        raFcn = raCos/rAabs * (1 - raTrans)
     END WHERE
     raZeta = raIntenP*(1-raTrans) + (raIntenP1 - raIntenP)*(raFcn - raTrans)
     raInten = raInten * exp(-raAbs/raCos) + raZeta

    ELSEIF (iVary == +41) THEN
      !!! this was done on Nov 04, 2014 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! PADE APPROX two term (combo of GENLN2 and LBLRTM)
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raTrans = exp(-raAbs)
      raZeta = 0.2*raAbs             !! pade one
      raFcn = (raIntenAvg + raZeta*raIntenP)/(1+raZeta)
      raZeta = 0.193*raAbs           !! pade two
      raZeta2 = 0.013*raAbs*raAbs    !! pade two
      raFcn = (raIntenAvg + (raZeta + raZeta2)*raIntenP)/(1+raZeta+raZeta2)
      raFcn = (1-raTrans)*raFcn
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +42) THEN
      !!! this was done on Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, GENLN2 style
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP)
      WHERE (raAbs >= 0.05)
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs
        raFcn = raAbs*raIntenP + raZeta*(1-raAbs/2) - raTrans * raZeta
      END WHERE
      raInten = raInten*raTrans + raFcn

    ELSEIF (iVary == +43) THEN
      !!! this was done on jan 2016 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, LBLRTM style, where for small OD (x)  means the function --> x/6
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = raIntenP - raIntenAvg
      WHERE (raAbs >= 0.06) 
        raTrans = exp(-raAbs)
        raZeta2 = 1.0 - 2.0*(1/raAbs - raTrans/(1-raTrans))
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*(raAbs * raAbs)
        raZeta2 = raAbs/6.0 - (raAbs**3)/360.0 + (raAbs**5)/15120.0  !! mathematica
        raZeta2 = raAbs/6.0
        raFcn = (1-raTrans)*(raIntenAvg + raZeta * raZeta2)
      END WHERE
      raInten = raInten*raTrans + raFcn

    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    !  y=1e-3; x = y : y : 250*y; T = exp(-x); plot(x,(1-T).*(1./x-T./(1-T)),'b.-',x,x.*(1/2-2*x/6+x.*x/6),'r',x,x/2,'k')
    ELSEIF (iVary == +4) THEN
      !!! this was done Oct 2015 .. looking at Clough et al, JGR 1992 v97
      !!! pg 15761, LBL calcs of atmospheric fluxed and cooling rates, Eqn 12
      !!! LINEAR IN TAU, MY style
      raAbs = raaAbs(:,iL)/raCos*rFrac
      raZeta = 2*(raIntenAvg-raIntenP)
      WHERE (raAbs > 0.1) 
        raTrans = exp(-raAbs)
        raFcn = (1-raTrans)*(raIntenP + raZeta/raAbs) - raTrans * raZeta
      ELSEWHERE
        raTrans = 1 - raAbs + 0.5*raAbs**2
        raZeta2 = raZeta*(raAbs/2-(raAbs**2)/3+(raAbs**3)/6)
        raFcn   = (1-raTrans)*raIntenP + raZeta2
      END WHERE
      raInten = raInten*raTrans + raFcn

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
      !pressure specified by user
      rP = raPressLevels(ip1)+rFrac*(raPressLevels(i0)-raPressLevels(ip1))
    ELSE                                !bot frac of layer
      !pressure specified by user
      rP = -rFrac*(raPressLevels(i0)-raPressLevels(ip1))+raPressLevels(i0)
    END IF

! compute the average pressure of the fractional layer
    IF (iTopOrBot == 1) THEN
      IF (abs(rP-raPressLevels(ip1)) >= delta) THEN
        rPavg = (rP-raPressLevels(ip1))/alog(rP/raPressLevels(ip1))
      ELSE
        rPavg = rP
      END IF
    ELSE
      IF (abs(rP-raPressLevels(i0)) >= delta) THEN
        rPavg = (raPressLevels(i0)-rP)/alog(raPressLevels(i0)/rP)
      ELSE
        rPavg = rP
      END IF
    END IF

! avg press,temperature of layer i0
    rP0 = (raPressLevels(i0)-raPressLevels(ip1))/alog(raPressLevels(i0)/raPressLevels(ip1))
    rT0 = raVTemp(i0+(iW-1)*kProfLayer)
! avg press, temperature of layer i0+1
    rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/alog(raPressLevels(ip1)/raPressLevels(ip1+1))
    rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
! surface parameters
    rPm1 =  rSurfPress
    rTm1 = rSurfTemp

! now compute the fit for rT(n) = ax(n)^2 + bx(n) + c where x(n) = alog(P(n))
    rPavg = alog(rPavg)

    rP0  = alog(rP0)
    rPp1 = alog(rPp1)
    rPm1 = alog(rPm1)
           
    xa(1) = rPp1
    xa(2) = rP0
    xa(3) = rPm1
    ya(1) = rTp1
    ya(2) = rT0
    ya(3) = rTm1
     
    CALL rspl_one(xa,ya,3,rPavg,rT,1)

    InterpTempSurf=rT
    RETURN
    end FUNCTION InterpTempSurf

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

    INTEGER, DIMENSION(:), ALLOCATABLE :: iaIndexAlloc
    INTEGER :: AllocateStatus,DeAllocateStatus

! name = sun view d(phi) windspeed
    fname = '/home/sergio/SBDART/V2.4/rho_22.12_23.86_212.97_7.3'
    iIOUN = kTempUnit
 1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='FORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
      WRITE(kStdErr,1010) IL, FNAME
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

      ALLOCATE ( iaIndexAlloc(iL), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory for iaIndexAlloc ***"      
      iaIndexAlloc = (/ (iI, iI = 1, iL) /)

      raTemp(iL-iaIndexAlloc+1) = raW(1:iL)
      raW(1:iL) = raTemp(1:iL)
      raTemp(iL-iaIndexAlloc+1) = raR(1:iL)

      DEALLOCATE (iaIndexAlloc, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Error while deallocating iaIndexAlloc ***"

      raR(1:iL) = raTemp(1:iL)
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
      !compute the Plank radiation from the sun
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      IF (raFreq(1) >= 605) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        !read in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
      ELSEIF (raFreq(1) < 605) THEN
        !! solar contribution is so small at these wavenumbers
        write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
        rSunTemp = kSunTemp
        !compute the Plank radiation from the sun
        raSun = ttorad(raFreq,rSunTemp)
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
    raSun = raSun*rCos*rOmegaSun      !!!!this is correct
    raKAbs = 0.0

    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
      iaRadLayerTemp,iT,iExtraSun,raExtraSun)
      
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
      ! the current defined atmosphere used all Gnd-100 layers
      DO iLay = iNumLayer,2,-1
        iL = iaRadLayer(iLay)
        rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raKAbs = raKAbs+raaAbs(:,iL)/rCos
      END DO
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raKAbs = raKAbs+raaAbs(:,iL)*rFracBot/rCos
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
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKabs = raKAbs+raaAbs(:,iL)/rCos
          raExtraSun = raSun*exp(-rakAbs)
        END DO
        !now do remaining layers, all the way to the ground-1
        DO iLay = iNumLayer-1,2,-1
          iL = iaRadLayer(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKAbs = raKAbs+raaAbs(:,iL)/rCos
        END DO
        DO iLay=1,1
          iL = iaRadLayer(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKAbs = raKAbs+raaAbs(:,iL)*rFracBot/rCos
          raSun = raSun*exp(-raKAbs)
        END DO
      END IF
         
      IF (iT > iNumLayer) THEN
        write(kStdWarn,*)'need to do the upper layers as well!!'
        !now do top layers, all the way to the instrument
        DO iLay = iT,iNumLayer+1,-1
          iL = iaRadLayerTemp(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKabs = raKAbs+raaAbs(:,iL)/rCos
        END DO
        !now do the layer instrument is in
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayerTemp(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKabs = raKAbs+raaAbs(:,iL)/rCos
          raExtraSun = raSun*(exp(-raKabs))
        END DO
        !now do all the way to the ground-1
        DO iLay = iNumLayer-1,2,-1
          iL = iaRadLayerTemp(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKabs = raKAbs+raaAbs(:,iL)/rCos
        END DO
        !now do ground
        DO iLay=1,1
          iL = iaRadLayerTemp(iLay)
          rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raKabs = raKAbs+raaAbs(:,iL)*rFracBot/rCos
          raSun = raSun*exp(-raKAbs)
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

 1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)

    iIOUN = kTempUnit
    CALL GetSolarFileName(fname,raFreq(1))
    write(kStdWarn,*) 'solar data file = ',fname
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
      WRITE(kStdErr,1010) IL, FNAME
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
    raSun = daSun*1000.0

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

    iDp = -1                   !assume nothing to be output

    IF (iNp < 0) THEN
      ! easy ! print the radiance at the end of this layer
      iDp = 1
      raOutFrac(iDp) = 1.0
    END IF

    IF (iNp > 0) THEN
      iDp = 0
      ! actually have to go thru list to see if this layer is to be output
      iDpC = 1
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
    REAL :: raPlanck(kMaxPts),raTrans(kMaxPts),raEmis(kMaxPts),rT,rFrac_k,rFrac_T
     
! iDir < 0  !radiance going down to instr on earth surface

! n all layers except bottommost layer, rFrac_k == rFrac_T
    rFrac_k = 0.0            !use this much in k interpolation
    rFrac_T = 0.0            !use this much in T interpolation

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
      raInten2 = raInten
    ELSE                           !interpolate
      IF (iLay /= 1) THEN
        !top part of most layers
        rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
      ELSE
        !bottom  part of top layer
        rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
      END IF
      write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
      raPlanck = ttorad(raFreq,rT)
      raTrans  = exp(-raaAbs(:,iL)*rFrac_k/rCos)
      raEmis   = (1.0-raTrans)*raPlanck
      raInten2 = raEmis+raInten*raTrans
      IF (iSun >= 0) THEN
        raTrans = exp(-raaAbs(:,iL)*rFrac_k/rCos)
        raInten2 = raInten2+raSun*raTrans
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
    REAL :: raPlanck(kMaxPts),raTrans(kMaxPts),raEmis(kMaxPts),rT,rFrac_k,rFrac_T

! iDir > 0  !radiance going up to instr in space
     
! n all layers except bottommost layer, rFrac_k == rFrac_T
    rFrac_k = 0.0            !use this much in k interpolation
    rFrac_T = 0.0            !use this much in T interpolation

    iL = iaRadLayer(iLay)

    IF ((iLay > 1) .AND. (iLay < iNumLayer)) THEN
      !no problem; full layer in mixtable
      rFrac_k = rFrac
      rFrac_T = rFrac
    ELSEIF (iLay == 1) THEN !!!bottommost layer
    
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
      raInten2 = raInten
    ELSE                           !interpolate
      IF (iLay /= 1) THEN
        !bottom part of most layers
        rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
      ELSE
        !top part of bottom layer
        rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
      END IF
      write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
      ! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
      ! so usually only the usual LTE computations are done!!
      IF (iNLTEStart > kProfLayer) THEN    !!!normal, no emission stuff
        raPlanck = ttorad(raFreq,rT)
        raTrans = exp(-raaAbs(:,iL)*rFrac_k/rCos)
        raEmis  = (1.0-raTrans)*raPlanck
        raInten2 = raEmis + raInten*raTrans
      ELSE IF (iNLTEStart <= kProfLayer) THEN
        raPlanck = ttorad(raFreq,rT)
        raTrans = exp(-raaAbs(:,iL)*rFrac_k/rCos)
        raEmis = (1.0-raTrans)*raPlanck*raaPlanckCoeff(:,iL)
        raInten2 = raEmis + raInten*raTrans
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

    OPEN(UNIT=iIOUN,FILE=FNCOFN,FORM='UNFORMATTED',STATUS='OLD',IOSTAT=IERR)
    IF (IERR /= 0) THEN
      WRITE(kStdErr,*) 'error reading SARTA NLTE coeffs file',IERR, FNCOFN
      CALL DoStop
    ENDIF

    kTempUnitOpen = +1
    iJ=1
    DO iI=1,MXCNTE
      ! Read data for this frequency/channel
      READ(iIOUN) ICHAN, FRQCHN, (COEFN(IC,iJ),IC=1,NNCOEF)
      raFrad(iI) = FRQCHN
      iJ = iJ + 1
    ENDDO
    NCHNTE = iJ - 1
    CLOSE(iIOUN)
    kTempUnitOpen = -1

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

    CALL rspl(raFrad,raDrad,NCHNTE,raFreq,raDkcarta,kMaxPts)    !! too dangerous,   small 4 um lte rads, wiggly NLTE correction
    CALL rlinear(raFrad,raDrad,NCHNTE,raFreq,raDkcarta,kMaxPts) !! hopefully safer, small 4 um lte rads, straightline NLTE correction
    DO iI = 1,kMaxPts
      IF ((raFreq(iI) < raFrad(1)) .OR. (raFreq(iI) > raFrad(NCHNTE))) THEN
        raDkcarta(iI) = 0.0
      END IF
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

    raExtra = 0.0
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
      ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
      iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4)) THEN
      ! even though the current defined atmosphere has all g-100 layers,
      ! 100th layer had frac 0 < f < 1
      iExtra=1
      ! extend the defined atmosphere so it includes all upper layers
      ! copy the currently defined atmosphere
      iT = 0
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
      iT = 0
      DO iI=1,iNumLayer
        iT = iT+1
        iaRadLayerTemp(iI) = iaRadLayer(iI)
      END DO
      ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
 15   CONTINUE
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
     
    iExtra = -1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    raExtra = 0.0
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
      ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
      iExtra = -1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4))THEN
      ! even though the current defined atmosphere has all g-100 layers,
      ! 100th layer had frac 0 < f < 1
      iExtra = 1
      ! extend the defined atmosphere so it includes all upper layers
      ! copy the currently defined atmosphere
      iT = 0
      DO iI=1,iNumLayer
        iT = iT+1
        iaRadLayerTemp(iI) = iaRadLayer(iI)
      END DO
      !        write(kStdWarn,*) 'top most layer is fractional layer. Some'
      !        write(kStdWarn,*) 'portion needed above instrument to calculate'
      !        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
      ! current defined atmosphere does not have all g-100 layers
      iExtra = 1
      ! extend the defined atmosphere so it includes all upper layers
      ! copy the currently defined atmosphere
      iT = 0
      DO iI=1,iNumLayer
        iT = iT+1
        iaRadLayerTemp(iI) = iaRadLayer(iI)
      END DO
      ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
 15   CONTINUE
      IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
        iT = iT+1
        iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
        !write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
        GO TO 15
      END IF
      ! write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
      ! write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF
     
    RETURN
    END SUBROUTINE 

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
          
    ! Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
      write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
      write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
      IF (WAVENO < WAVETAB(1)) THEN
        WAVENO = WAVETAB(1)
      ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
        WAVENO = WAVETAB(NWAVE)
      END IF
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
      write(kStdErr,*) DME,' outside ',DMETAB(1),':',DMETAB(NDME)
      write(kStdErr,*) 'INTERP_SCAT_TABLE: particle Dme out of range ... RESET'
      IF (DME < DMETAB(1)) THEN
        DME = DMETAB(1)
      ELSEIF (DME > DMETAB(NDME)) THEN
        DME = DMETAB(NDME)
      END IF
    END IF

    ! See if wavenumber is within last wavenumber grid, otherwise
    ! find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
    ! IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
      IL = 1
      IU = NWAVE
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
      !Find the grid location and interpolation factor for DME
      IL = 1
      IU = NDME
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
      ! If not the same Dme or a new wavenumber grid, then
      ! linearly interpolate omega and g and log interpolate extinction
      EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
      EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
      ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
      ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
      ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
      ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)
    ENDIF

    ! Linearly interpolate the scattering properties in wavenumber
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

    Temp1 = -10.0
    Temp(1:kProfLayer) = -10.0

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
      !map this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
      iL = iL-idiv(iL,kProfLayer)*kProfLayer
      IF (iL == 0) THEN
        iL = kProfLayer
      END IF
      rP=raPressLevels(iL+1)-10000*delta
      if (rp < raPressLevels(kProfLayer+1)) then
        rp = raPressLevels(kProfLayer+1)+10000*delta
      end if
      TEMP1(iNumLayer-iLay+1) = FindBottomTemp(rP,raProfileTemp,raPressLevels,iProfileLayers)
    END DO

    rP = raPressLevels(iLowest)
    rP = DISORTsurfPress          !!!from scatterparam.f90
    TEMP1(iNumLayer+1) = FindBottomTemp(rP,raProfileTemp,raPressLevels,iProfileLayers)

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
      !this is WHOLE of the bottom layer
      i1 = iLowest
    ELSE IF (rP <= raPressLevels(kProfLayer+1)) THEN
      !this is ludicrous
      write(kStdErr,*) rP,raPressLevels(kProfLayer+1)
      write(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
      CALL DoStop

    ELSE
      ! first find the AIRS layer within which it lies
      iFound = -1
      i1 = iLowest
      i2 = iLowest+1
      10 CONTINUE
      IF ((rP <= raPressLevels(i1)) .AND. (rP > raPressLevels(i2))) THEN
        iFound=1
      END IF
      IF ((iFound < 0) .AND. (i1 < kProfLayer)) THEN
        i1 = i1+1
        i2 = i2+1
        GOTO 10
      END IF
      IF ((iFound < 0)) THEN
        IF (abs(rP-raPressLevels(kProfLayer+1)) <= delta) THEN
          i1 = kProfLayer
          iFound = 1
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
      write(kStdErr,*) 'Allowed Pressure ranges are from : ',&
        raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
      write(kStdErr,*) 'Surface Pressure is ',rP,' mb'
      call DoStop
    END IF
              
! now find the temperature
    IF (i1 == iLowest) THEN          !do linear interp
      i1 = iLowest
      i2 = iLowest+1
      i3 = iLowest+2
      rP1 = (raPressLevels(i2)-raPressLevels(i1))/log(raPressLevels(i2)/raPressLevels(i1))
      rP2 = (raPressLevels(i3)-raPressLevels(i2))/log(raPressLevels(i3)/raPressLevels(i2))
      T1 = raProfileTemp(i1)
      T2 = raProfileTemp(i2)
      IF (iLog == -1) THEN
        rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)           !!linear in P
      ELSE
        rT = T2-(log(rP2/rP))*(T2-T1)/(log(rP2/rP1)) !!log(P)
      END IF

    ELSEIF (i1 >= (kProfLayer-1)) THEN          !do linear interp
      rP1 = (raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/log(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
      rP2 = (raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/log(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
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
      		     
      rP1 = (raPressLevels(i1)-raPressLevels(i1-1))/log(raPressLevels(i1)/raPressLevels(i1-1))
      raP(3) = rP1
      rP1 = (raPressLevels(i1+1)-raPressLevels(i1))/log(raPressLevels(i1+1)/raPressLevels(i1))
      raP(2) = rP1
      rP1 = (raPressLevels(i1+2)-raPressLevels(i1+1))/log(raPressLevels(i1+2)/raPressLevels(i1+1))
      raP(1) = rP1
      IF (iLog == +1) THEN
        DO iJ = 1,3
          raLogP(iJ) = log(raP(iJ))
        END DO
      END IF

      raT(3) = raProfileTemp(i1-1)
      raT(2) = raProfileTemp(i1)
      raT(1) = raProfileTemp(i1+1)

      yp1 = 1.0e30
      ypn = 1.0e30
      IF (iSpline == +1) THEN
        IF (iLog == +1) THEN
          CALL rspl_one(raLogP,raT,3,log(rP),rT,1)
        ELSE
          CALL rspl_one(raP,raT,3,rP,rT,1)
        END IF
      ELSEIF (iSpline == -1) THEN
        IF (iLog == +1) THEN
          CALL rlinear_one(raP,raT,3,rP,rT,1)
        ELSE
          CALL rlinear_one(raLogP,raT,3,log(rP),rT,1)
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
      !have to adjust temperature .. do this for down AND up look instr
      IF (rPressStart > rPressStop) THEN  ! for down looking instr
        rT = FindBottomTemp(rPressStart,raProfileTemp,raPressLevels,iProfileLayers)
        rT = rT+rTSurf
      ELSEIF (rPressStart < rPressStop) THEN  ! for up looking instr
        rT = FindBottomTemp(rPressStop,raProfileTemp,raPressLevels,iProfileLayers)
        rT = rT+rTSurf
      END IF
    END IF
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

    ! Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
      write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
      write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
      IF (WAVENO < WAVETAB(1)) THEN
        WAVENO = WAVETAB(1)
      ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
        WAVENO = WAVETAB(NWAVE)
      END IF
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
      IF (DME < DMETAB(1)) THEN
        DME = DMETAB(1)
      ELSEIF (DME > DMETAB(NDME)) THEN
        DME = DMETAB(NDME)
      END IF
    END IF

    ! See if wavenumber is within last wavenumber grid, otherwise
    !   find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
    !IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
      IL = 1
      IU = NWAVE
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
      !Find the grid location and interpolation factor for DME
      IL = 1
      IU = NDME
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
      !If not the same Dme or a new wavenumber grid, then
      !  linearly interpolate omega and g and log interpolate extinction
      EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
      EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
      ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
      ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
      ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
      ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)

      ! looking at sarta code, Scott Hannon ALWAYS does a log interp
    ELSEIF ((DME /= OLDDME .OR. NEWIW) .AND. (iLogOrLinear == -1)) THEN
      !If not the same Dme or a new wavenumber grid, then
      !  linearly interpolate omega and g and log interpolate extinction
      EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
      EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
      ALB0 = EXP( (1-FLDME)*LOG(TABSSALB(IW0,ID))   + FLDME*LOG(TABSSALB(IW0,ID+1)) )
      ALB1 = EXP( (1-FLDME)*LOG(TABSSALB(IW1,ID))   + FLDME*LOG(TABSSALB(IW1,ID+1)) )
      ASYM0 = EXP( (1-FLDME)*LOG(TABASYM(IW0,ID))   + FLDME*LOG(TABASYM(IW0,ID+1)) )
      ASYM1 = EXP( (1-FLDME)*LOG(TABASYM(IW1,ID))   + FLDME*LOG(TABASYM(IW1,ID+1)) )
    ENDIF

    ! Linearly interpolate the scattering properties in wavenumber
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

    OPEN (UNIT = kTempUnit, STATUS='OLD', FORM='UNFORMATTED', FILE=SCATFILE, IOSTAT=IERR)
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
 
    !print *,NMUOBS,NDME,NWAVE
    !print *,MUINC(1),MUINC(2)
    !print *,cScale
    !print *,(MUTAB(IMU), IMU = 1, NMUOBS)

    DO IW = 1, NWAVE
      DO ID = 1, NDME
        K2 = IW-1 + NWAVE*(ID-1)
        K3 = NMUOBS*K2
        READ(kTempUnit) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1), &
        TABSSALB(K2+1), TABASYM(K2+1)
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
        READ(2,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1),TABSSALB(K2+1), TABASYM(K2+1)
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
    IF (MAX(NMUOBS,NDME,NWAVE) > MAXGRID) STOP 'READ_SSCATTAB_SPECIAL: MAXGRID exceeded'
    IF (NMUOBS*NDME*NWAVE > MAXTAB) STOP 'READ_SSCATTAB_SPECIAL: MAXTAB exceeded'
    READ (2,*)
    READ (2,*)
    READ (2,*)
    READ (2,*)

    DO IW = 1, NWAVE
      DO ID = 1, NDME
        K2 = IW-1 + NWAVE*(ID-1)
        READ(2,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1), TABSSALB(K2+1), TABASYM(K2+1)
        TABEXTINCT(K2+1) = TABEXTINCT(K2+1) * 1000.0   !!!to be like sscatmie
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
    REAL :: raEmission(kMaxPts),raTrans(kMaxPts),rMu,raInten0(kMaxPts)
    DOUBLE PRECISION :: daInten(kMaxPts),daTrans(kMaxPts),daEmission(kMaxPts)
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
    raInten0 = raInten
    daInten  = dble(raInten)

    iL = 0
    IF (kNLTEOutUAOpen > 0) THEN
      write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
      ! always dump out the 0.005 mb TOA radiance if the UA file is open
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF

    rMu = cos(rSatAngle*kPi/180.0)

    DO iL = 1,iUpper - 1

      raTrans = raaUpperSumNLTEGasAbCoeff(:,iL)/rMu
      raTrans = exp(-raTrans)
      raEmission = (1.0 - raTrans) * raaUpperPlanckCoeff(:,iL) * ttorad(raFreq,raUpperTemp(iL))
      raInten = raEmission + raInten*raTrans

      daTrans = (raaUpperSumNLTEGasAbCoeff(:,iL)*1.0d0/(rMu*1.0d0))
      daTrans = exp(-daTrans)
      daEmission = (raaUpperPlanckCoeff(:,iL)*1.0d0) * dble(ttorad(raFreq,raUpperTemp(iL))*1.0d0)*(1.0d0 - daTrans)
      daInten = daEmission + daInten*daTrans

      raInten = sngl(daInten)

      IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
        write(kStdWarn,*) 'dumping out UA rads at iL = ',iL
        ! dump out the radiance at this HIGH pressure level
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF

    END DO

    DO iL = iUpper,iUpper
      raTrans = raaUpperSumNLTEGasAbCoeff(:,iL)/rMu
      raTrans = exp(-raTrans)
      raEmission = (1.0 - raTrans) * raaUpperPlanckCoeff(:,iL) * ttorad(raFreq,raUpperTemp(iL))
      raInten = raEmission + raInten*raTrans

      daTrans = dble(raaUpperSumNLTEGasAbCoeff(:,iL)*1.0d0/(rMu*1.0d0))
      daTrans = exp(-daTrans)
      daEmission = dble(raaUpperPlanckCoeff(:,iL)*1.0d0) * dble(ttorad(raFreq,raUpperTemp(iL))*1.0d0)*(1.0d0 - daTrans)
      daInten = daEmission + daInten*daTrans
      raInten = sngl(daInten)

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
