! Copyright 2001
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:46
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
!************************************************************************

! this does the scattering radiative transfer for a downlook instrument
! the only difference between Cloud_DownLook and Cloud_UpLook are that the
! ordering of the layers in iaRadLayer is inverted, and so
! iLocalCldTop and iLocalCldBot are inverted

! same as /home/sergio/MATLAB/RADTrans/DOWN/cloudyD25_best_prof_solar.m

SUBROUTINE Cloud_UpOrDownLook(iNumLayer,iDir,iLocalCldTop,iLocalCldBot,  &
    rFracTop,rFracBot, iaRadLayer,raLayAngles,TEMP,raFreq,  &
    raaExt,raaScat,raaAsym,radSolarCld,muSun,muSat,  &
    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,  &
    raTrDown12,raReDown12,raEmissDown12,raSunDown12, raW0,raAsym,  &
    iPhase,raPhasePoints,raComputedPhase)


INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN)                      :: iDir
NO TYPE, INTENT(IN OUT)                  :: iLocalCldT
NO TYPE, INTENT(IN OUT)                  :: iLocalCldB
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaScat(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: radSolarCl
REAL, INTENT(IN OUT)                     :: muSun
REAL, INTENT(IN OUT)                     :: muSat
REAL, INTENT(IN OUT)                     :: raTau12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTrUp12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raReUp12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raEmissUp1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunUp12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTrDown12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raReDown12(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
REAL, INTENT(IN OUT)                     :: raSunDown1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raW0(kMaxPts)
REAL, INTENT(IN OUT)                     :: raAsym(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input parameters


INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
INTEGER :: !atmosphere layering
REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
REAL :: !temperature profile (levels)
REAL :: !wavenumbers
REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
!these next three are self explanatory



INTEGER :: iPhase                  !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)


! output parameters

!these next few are self explanatory : optical depth, and the
!cumulative up/down transmission, reflection, emission, solar trans

REAL :: raEmissUp12(kMaxPts)
REAL :: raEmissDown12(kMaxPts)
REAL :: raSunDown12(kMaxPts)
! these are the lowest cloud layer asymmetry and single scattering


! local variables
INTEGER :: N,iFr,iLay,iL,iDp

! general coeffs for the layers
REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raTb(kMaxPts),raTt(kMaxPts)
REAL :: raBeta(kMaxPts),raSun0(kMaxPts)
REAL :: raA(kMaxPts),raB(kMaxPts),raE(kMaxPts),raF(kMaxPts)
REAL :: raG(kMaxPts),raH(kMaxPts)
REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts)

! arbitrary angle stuff
REAL :: raTau1(kMaxPts),raTau2(kMaxPts)
REAL :: raTrUp1(kMaxPts),raReUp1(kMaxPts)
REAL :: raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
REAL :: raTrDown1(kMaxPts),raReDown1(kMaxPts),raEmissDown1(kMaxPts)
REAL :: raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
REAL :: raSunUp1(kMaxPts),raSunUp2(kMaxPts)
REAL :: raSunDown2(kMaxPts)

! stream angle stuff
REAL :: raT1(kMaxPts),raT1star(kMaxPts),raT2(kMaxPts),raT2star(kMaxPts)
REAL :: raR1(kMaxPts),raR1star(kMaxPts),raR2(kMaxPts),raR2star(kMaxPts)
REAL :: raE1p(kMaxPts),raE1m(kMaxPts),raE2p(kMaxPts),raE2m(kMaxPts)
REAL :: raF1p(kMaxPts),raF1m(kMaxPts),raF2p(kMaxPts),raF2m(kMaxPts)
REAL :: raT12(kMaxPts),raT12star(kMaxPts)
REAL :: raR12(kMaxPts),raR12star(kMaxPts)
REAL :: raE12p(kMaxPts),raE12m(kMaxPts)
REAL :: raS1p(kMaxPts),raS2p(kMaxPts),raS12p(kMaxPts)
REAL :: raS1m(kMaxPts),raS2m(kMaxPts),raS12m(kMaxPts)
REAL :: raCumSum1(kMaxPts),raCumSum2(kMaxPts),raDet(kMaxPts)

IF (iDir > 0) THEN
  N = iLocalCldTop - iLocalCldBot + 1
ELSE IF (iDir < 0) THEN
  N = iLocalCldBot - iLocalCldTop + 1
ELSE IF (iDir == 0) THEN
  WRITE (kStdErr,*) 'iDir = ',iDir,' invalid number!!!'
  CALL DoStop
END IF

DO iFr = 1,kMaxPts
  raCumSum1(iFr) = 0.0
  raCumSum2(iFr) = 0.0
END DO

IF (N < 1) THEN
  WRITE(kStdErr,*) 'Huh? negative number of cld lays in Cld_UpDnLook'
  WRITE(kStdErr,*) 'Local CldTop,CldBot = ',iLocalCldTop,iLocalCldBot
  CALL DoStop
END IF

! -------------- cloud has only one layer ---------------------------------
IF (N == 1) THEN             !bloody simple
  iLay    = iLocalCldTop
  
  CALL OneLayerScatProp(iLay,iDir,N,2,1,
! input params  &
  iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP,  &
      raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld,  &
      iPhase,raPhasePoints,raComputedPhase,
! output params  &
  raCumSum1,raCumSum2,                         !AccumulateSolarDept  &
      raTau1,raTau2,muSat,                         !optical depths,angle  &
      raW0,raAsym,raTb,raBeta,                     !LayerScatteringProp  &
      raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, !asymcoeffs_solar  &
      raTrUp2,raReUp2,raEmissUp2,raSunUp2,         !t_r_e_arb_up_Solar_prof  &
      raTrDown2,raReDown2,raEmissDown2,raSunDown2, !t_r_e_arb_dn_Solar_prof  &
      raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
!if more than one layer, need t_r_e_streams_Solar_prof
  
!need to assign raXXX12 --> raXXX2 for output purposes
  CALL update_streams_arb_solar(  &
      raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12,  &
      raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,  &
      raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2,  &
      raEmissUp2, raEmissDown2, raSunup2, raSunDown2)
  
END IF

! -------------- cloud has many layers ---------------------------------
IF (N > 1) THEN             !have to add the layers together
!because the sun takes away up down symmetry, we have to do this
!from the top down
  
!do the stuff for the arbitrary angles, for CLOUD LAYER 2
  iLay    = iLocalCldTop
  CALL OneLayerScatProp(iLay,iDir,N,2,1,
! input params  &
  iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP,  &
      raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld,  &
      iPhase,raPhasePoints,raComputedPhase,
! output params  &
  raCumSum1,raCumSum2,                         !AccumulateSolarDepth  &
      raTau1,raTau2,muSat,                         !optical depths,angle  &
      raW0,raAsym,raTb,raBeta,                     !LayerScatteringProp  &
      raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, !asymcoeffs_solar  &
      raTrUp2,raReUp2,raEmissUp2,raSunUp2,         !t_r_e_arb_up_Solar_prof  &
      raTrDown2,raReDown2,raEmissDown2,raSunDown2, !t_r_e_arb_dn_Solar_prof  &
      raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
!if more than one layer, need t_r_e_streams_Solar_prof
  
!do the stuff for the arbitrary angles, for CLOUD LAYER 1
  iLay = iLay - iDir
  CALL OneLayerScatProp(iLay,iDir,N,1,2,
! input params  &
  iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP,  &
      raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld,  &
      iPhase,raPhasePoints,raComputedPhase,
! output params  &
  raCumSum1,raCumSum2,                         !AccumulateSolarDepth  &
      raTau1,raTau2,muSat,                         !optical depths,angle  &
      raW0,raAsym,raTb,raBeta,                     !LayerScatteringProp  &
      raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, !asymcoeffs_solar  &
      raTrUp1,raReUp1,raEmissUp1,raSunUp1,         !t_r_e_arb_up_Solar_prof  &
      raTrDown1,raReDown1,raEmissDown1,raSunDown1, !t_r_e_arb_dn_Solar_prof  &
      raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
!if more than one layer, need t_r_e_streams_Solar_prof
  
!loop over the layers
  DO  iDp = 1,N-1
!!! add layers together
    
    CALL addstar_arb_solar(  &
        raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1,  &
        raSunUp1,raSunDown1,raTau1,  &
        raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2,  &
        raSunUp2,raSunDown2,raTau2,  &
        raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m,  &
        raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,  &
        raReUp12,raReDown12,raTrUp12,raTrDown12,  &
        raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12,  &
        raCumSum1,raCumSum2,muSat,muSun)
    
    CALL addstar_solar(  &
        raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m,  &
        raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m,  &
        raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,  &
        raCumSum1,raCumSum2,muSun)
    
    iLay = iLay - iDir
    IF (iDp < (N-1)) THEN          !!!!!add on new layer stuff
      CALL OneLayerScatProp(iLay,iDir,N,1,3,
! input params  &
      iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP,  &
          raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld,  &
          iPhase,raPhasePoints,raComputedPhase,
! output params  &
      raCumSum1,raCumSum2,                        !AccumulateSolarDepth  &
          raTau1,raTau2,muSat,                        !optical depths  &
          raW0,raAsym,raTb,raBeta,                    !LayerScatteringProp  &
          raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, !asym_solar  &
          raTrUp1,raReUp1,raEmissUp1,raSunUp1,         !arb_up_Solar_prof  &
          raTrDown1,raReDown1,raEmissDown1,raSunDown1, !arb_dn_Solar_prof  &
          raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
!if more than one layer, need t_r_e_streams_Solar_prof
      
      CALL update_streams_solar(  &
          raR2, raR2star, raT2, raT2star, raE2p, raE2m, raS2p, raS2m,  &
          raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
      CALL update_streams_arb_solar(  &
          raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2,  &
          raEmissUp2, raEmissDown2, raSunUp2, raSunDown2,  &
          raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12,  &
          raEmissUp12,raEmissDown12,raSunup12,raSunDown12)
      
    END IF
  END DO
END IF

RETURN
END SUBROUTINE Cloud_UpOrDownLook

!************************************************************************
! this finds out if we use complete HG function, or first order approx

SUBROUTINE FindIHG(avg,iHG,iArbHG,raW0)


REAL, INTENT(OUT)                        :: avg
INTEGER, INTENT(OUT)                     :: iHG
INTEGER, INTENT(OUT)                     :: iArbHG
REAL, INTENT(IN)                         :: raW0(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'





avg = (raW0(1)+raW0(kMaxPts))/2.0

!c      IF ((kSolar .GE. 0) .AND. (avg .LE. 0.75)) THEN
!c        !!use full HG phase function
!c        iHG    = 1
!c        iArbHG = 1
!c        !!use full HG phase function for twostream only
!c        iHG    = 1
!c        iArbHG = 1
!c      ELSE
!c        !!use first order HG phase function
!c        iHG    = -1
!c        iArbHG = -1
!c        END IF

!c ******* this was sept 2001 *********
IF (kSolar >= 0) THEN
!!use full HG phase function
  iHG    = 1
  iArbHG = 1
ELSE
!!use first order HG phase function
  iHG    = -1
  iArbHG = -1
END IF
!c ******* this was sept 2001 *********

!c ******* try this dec 2001 *********
!!use full HG phase function
iHG    = 1
iArbHG = 1
!c ******* try this dec 2001 *********

RETURN
END SUBROUTINE FindIHG

!************************************************************************
! this checks to see if the emissions are positive or not
! flag the negatives as "bad" by keeping a running count of wot's bad

SUBROUTINE check_raE_raF(raE,raF,IF,iEminus,iFminus,iaEMinus,iaFMinus,  &
    rEmin_positive_for_chunk,rFmin_positive_for_chunk)


REAL, INTENT(IN)                         :: raE(kMaxPts)
REAL, INTENT(IN)                         :: raF(kMaxPts)
INTEGER, INTENT(IN)                      :: IF
INTEGER, INTENT(OUT)                     :: iEminus
INTEGER, INTENT(OUT)                     :: iFminus
INTEGER, INTENT(OUT)                     :: iaEMinus(kMaxPts)
INTEGER, INTENT(OUT)                     :: iaFMinus(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: rEmin_posi
NO TYPE, INTENT(IN OUT)                  :: rFmin_posi
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input params
REAL :: !!the up/down emission

! output params
REAL :: rEmin_positive_for_chunk,rFmin_positive_for_chunk



!      raE(iF) = max(raE(iF),0.0)
!      raF(iF) = max(raF(iF),0.0)
IF ((raE(IF) > 0).AND.(raE(IF) < rEmin_positive_for_chunk))THEN
  rEmin_positive_for_chunk = raE(IF)
END IF
IF ((raF(IF) > 0).AND.(raF(IF) < rFmin_positive_for_chunk))THEN
  rFmin_positive_for_chunk = raF(IF)
END IF

IF (raE(IF) < 0) THEN
  iEminus = iEMinus + 1
  iaEMinus(iEMinus) = IF
END IF
IF (raF(IF) < 0) THEN
  iFminus = iFMinus + 1
  iaFMinus(iFMinus) = IF
END IF

RETURN
END SUBROUTINE check_raE_raF

!************************************************************************
! this subroutine computes the coefficients
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/asymcoeffs_solar.m
! here, view angle === stream angle

SUBROUTINE asymcoeffs_solar(raAsym,raW0,raBo,raRadBb,raRadBt,  &
    raTau,raBeta,radSolarCld,raCumSum,muSun,  &
    iPhase,raPhasePoints,raComputedPhase,  &
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)


REAL, INTENT(IN)                         :: raAsym(kMaxPts)
REAL, INTENT(IN)                         :: raW0(kMaxPts)
REAL, INTENT(IN)                         :: raBo(kMaxPts)
REAL, INTENT(IN)                         :: raRadBb(kMaxPts)
REAL, INTENT(IN)                         :: raRadBt(kmaxPts)
REAL, INTENT(IN)                         :: raTau(kMaxPts)
REAL, INTENT(IN)                         :: raBeta(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: radSolarCl
REAL, INTENT(IN)                         :: raCumSum(kMaxPts)
NO TYPE, INTENT(OUT)                     :: muSun
NO TYPE, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
REAL, INTENT(OUT)                        :: raA(kMaxPts)
REAL, INTENT(OUT)                        :: raB(kMaxPts)
REAL, INTENT(OUT)                        :: raDet(kMaxPts)
REAL, INTENT(OUT)                        :: raE(kMaxPts)
REAL, INTENT(OUT)                        :: raF(kMaxPts)
REAL, INTENT(OUT)                        :: raG(kMaxPts)
REAL, INTENT(OUT)                        :: raH(kMaxPts)
REAL, INTENT(OUT)                        :: raKp(kMaxPts)
REAL, INTENT(OUT)                        :: raKm(kMaxPts)
REAL, INTENT(OUT)                        :: raAp(kMaxPts)
REAL, INTENT(OUT)                        :: raAm(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! really we do not need raA or raB, so we do not really need RadBb or RadBt
! so we do not need to compute rho or theta

! input parameters
REAL :: !asym coeff, single scatter alb
REAL :: !1/k * log(T(n+1)/T(n))
REAL :: !bottom level temperature
REAL :: !planck, top,bottom rad
REAL :: !total extinction of layer
REAL :: radSolarCld(kMaxPts),muSun      !sun intensity at cloud top, angle
REAL :: !adjusts sun intensity as necessary
INTEGER :: iPhase                  !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! output parameters
REAL :: !!det is more important!
REAL :: !boundary condtions for soln
REAL :: !boundary condtions for soln
REAL :: !eigen stuff

! local variables
REAL :: mu_pm                             !two stream angle
INTEGER :: IF,iHG,iArbHG,iReset_raEF
REAL :: b,det,alpha,theta,rho,avg,betaa,angle,p_plus,p_minus,hg2,psi,rSun
REAL :: avgmin,rEmin_positive_for_chunk,rFmin_positive_for_chunk
INTEGER :: iaEMinus(kMaxPts),iaFMinus(kMaxPts),iEMinus,iFMinus

INTEGER :: i1,i2,iDebug
REAL :: z1,z2,z3,z4

iDebug = -1
i1 = 9222
i2 = 9224

rEmin_positive_for_chunk = +1.0*(100.0**8)
rFmin_positive_for_chunk = +1.0*(100.0**8)

iEMinus = 0
iFMinus = 0

avgmin = kAvgMin

mu_pm = 1/SQRT(3.0)
muSun = ABS(muSun)

CALL FindIHG(avg,iHG,iArbHG,raW0)

IF ((avg >= avgmin) .AND. (kSolar >= 0)) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
    b = (1-raAsym(IF))/2.0
    raKp(IF) = 1/mu_pm * SQRT((1-raW0(IF))*(1-raW0(IF)*raAsym(IF)))
    raKm(IF) =-1/mu_pm * SQRT((1-raW0(IF))*(1-raW0(IF)*raAsym(IF)))
    
    alpha = (raW0(IF)*(1-b)-1)
    raAp(IF) = -(mu_pm/(raW0(IF)*b))*(raKp(IF) + alpha/mu_pm)
    raAm(IF) = -(mu_pm/(raW0(IF)*b))*(raKm(IF) + alpha/mu_pm)
    
!          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
    det = (raW0(IF)*b)*(raW0(IF)*b)-alpha*alpha+(betaa*mu_pm)*(betaa*mu_pm)
    
    raE(IF) = raBo(IF)*(1-raW0(IF))/det
    raE(IF) = raE(IF)*(betaa*mu_pm + alpha - raW0(IF)*b)
    
    raF(IF) = raBo(IF)*(1-raW0(IF))/det
    raF(IF) = raF(IF)*(-betaa*mu_pm + alpha - raW0(IF)*b)
    
    CALL check_raE_raF(raE,raF,IF,iEminus,iFminus,iaEMinus,iaFMinus,  &
        rEmin_positive_for_chunk,rFmin_positive_for_chunk)
    
    rSun = radSolarCld(IF) * EXP(-raCumSum(IF)/muSun)
    
    IF (iHG <= 0) THEN
!! before this loop was if (iHG .GT. 0) THEN
!! but gave lousy answers when we fix the minus sign (see below)
      angle = mu_pm/muSun
      det = angle**2 - alpha**2 + (raW0(IF)*b)**2
      p_plus  = hg2(-ABS(muSun),+ABS(mu_pm),raAsym(IF),  &
          iPhase,raPhasePoints,raComputedPhase)
      p_minus = hg2(-ABS(muSun),-ABS(mu_pm),raAsym(IF),  &
          iPhase,raPhasePoints,raComputedPhase)
!!! minus sign is incorrect, but gives best answers if we use this!
!!! raG(iF) = (alpha+angle)*p_plus - raW0(iF)*b*p_minus
      raG(IF) = (alpha+angle)*p_plus + raW0(IF)*b*p_minus
      raG(IF) = raW0(IF)*rSun/kForP/det * raG(IF)
      raH(IF) = raW0(IF)*b*(-p_plus) + (alpha-angle)*p_minus
      raH(IF) = raW0(IF)*rSun/kForP/det * raH(IF)
!z1 = raG(iF)
!z2 = raH(iF)
      
    ELSE
      angle = mu_pm/muSun
      det = angle**2 - alpha**2 + (raW0(IF)*b)**2
      psi = 3.0*raAsym(IF)*muSun*mu_pm
      z3 = -alpha + raW0(IF)*b
      z4 =  alpha + raW0(IF)*b
      raG(IF) =  z3 + psi*z4 + angle*(psi-1)
      raG(IF) = -raW0(IF)*rSun/kForP/det * raG(IF)
      raH(IF) =  z3 - psi*z4 + angle*(psi+1)
      raH(IF) = -raW0(IF)*rSun/kForP/det * raH(IF)
      
!            print *,iF,z1,z2,raG(iF),raH(iF),p_plus,p_minus,1-psi,1+psi
    END IF
    
    theta = raRadBb(IF) - raE(IF) - raG(IF)*EXP(-raTau(IF)/muSun)
    rho   = raRadBt(IF) - raF(IF)*EXP(betaa*raTau(IF)) - raH(IF)
    
    raDet(IF) = raAp(IF)*EXP(raKm(IF)*raTau(IF)) -  &
        raAm(IF)*EXP(raKp(IF)*raTau(IF))
    
! these last two are really not needed, but what the heck!!!!!
    raA(IF) = ( theta*EXP(raKm(IF)*raTau(IF)) - rho*raAm(IF))
    raA(IF) = raA(IF)/raDet(IF)
    raB(IF) = (-theta*EXP(raKp(IF)*raTau(IF)) + rho*raAp(IF))
    raB(IF) = raB(IF)/raDet(IF)
  END DO
  
ELSE IF ((avg >= avgmin) .AND. (kSolar < 0)) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
    b = (1.0-raAsym(IF))/2.0
    raKp(IF) = 1/mu_pm * SQRT((1-raW0(IF))*(1-raW0(IF)*raAsym(IF)))
    raKm(IF) =-1/mu_pm * SQRT((1-raW0(IF))*(1-raW0(IF)*raAsym(IF)))
    
    alpha = raW0(IF)*(1.0-b)-1.0
    raAp(IF) = -(mu_pm/(raW0(IF)*b))*(raKp(IF) + alpha/mu_pm)
    raAm(IF) = -(mu_pm/(raW0(IF)*b))*(raKm(IF) + alpha/mu_pm)
    
!          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
    det = (raW0(IF)*b)*(raW0(IF)*b)-alpha*alpha+(betaa*mu_pm)*(betaa*mu_pm)
!          print *,iF,det
    
    raE(IF) = raBo(IF)*(1-raW0(IF))/det
    raE(IF) = raE(IF)*(betaa*mu_pm + alpha - raW0(IF)*b)
    
    raF(IF) = raBo(IF)*(1-raW0(IF))/det
    raF(IF) = raF(IF)*(-betaa*mu_pm + alpha - raW0(IF)*b)
    
    CALL check_raE_raF(raE,raF,IF,iEminus,iFminus,iaEMinus,iaFMinus,  &
        rEmin_positive_for_chunk,rFmin_positive_for_chunk)
    
    IF (iDebug > 0) THEN
      IF ((IF >= i1) .AND. (IF <= i2)) THEN
        z1 = raBo(IF)*(1-raW0(IF))/det
        z2 = (+betaa*mu_pm + alpha - raW0(IF)*b)
        z3 = (-betaa*mu_pm + alpha - raW0(IF)*b)
!              print *,iF,det,z1,z2,z3,raE(iF),raF(iF)
        PRINT *,IF,raW0(IF),raAsym(IF),b,alpha,betaa,det,raE(IF)
      END IF
    END IF
    
    raG(IF) = 0.0
    raH(IF) = 0.0
    
    theta = raRadBb(IF) - raE(IF) - raG(IF)*EXP(-raTau(IF)/muSun)
    rho   = raRadBt(IF) - raF(IF)*EXP(betaa*raTau(IF)) - raH(IF)
    
    raDet(IF) = raAp(IF)*EXP(raKm(IF)*raTau(IF)) -  &
        raAm(IF)*EXP(raKp(IF)*raTau(IF))
    
!!!! these last two are really not needed, but what the heck!!!!!
    raA(IF) = ( theta*EXP(raKm(IF)*raTau(IF)) - rho*raAm(IF))
    raA(IF) = raA(IF)/raDet(IF)
    raB(IF) = (-theta*EXP(raKp(IF)*raTau(IF)) + rho*raAp(IF))
    raB(IF) = raB(IF)/raDet(IF)
  END DO
  
ELSE
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
    b = (1-raAsym(IF))/2
    raKp(IF) = 1/mu_pm
    raKm(IF) =-1/mu_pm
    
    raAp(IF) = 0.0
    raAm(IF) = 1.0
    
    raE(IF) = raBo(IF)/(1+betaa*mu_pm)
    raF(IF) = raBo(IF)/(1-betaa*mu_pm)
    
    raG(IF) = 0.0
    raH(IF) = 0.0
    
    theta = raRadBb(IF) - raE(IF)
    rho   = raRadBt(IF) - raF(IF)*EXP(raTau(IF)*betaa)
    raDet(IF) = raAp(IF)*EXP(raKm(IF)*raTau(IF)) -  &
        raAm(IF)*EXP(raKp(IF)*raTau(IF))
!!! these last two are really not needed, but what the heck!!!!!
    raA(IF) = rho*EXP(raKm(IF)*raTau(IF))
    raB(IF) = theta
    
  END DO
END IF

!      DO iF = 1,kMaxPts
!        print *,iF,raE(iF),raG(iF)
!        END DO
!      CALL DoStop

! if there are "bad" emissions, set them to min(good emissions)
iReset_raEF = +1     !!go ahead and reset things
iReset_raEF = -1     !!leave things as they are
IF ((avg >= avgmin) .AND. (iReset_raEF > 0)) THEN
  IF (iEminus > 0) THEN
    WRITE(kStdWarn,*) 'oops, layer has negative E emiss!! "fixing" this'
    WRITE(kStdWarn,*) 'by setting raE(bad) = ',rEmin_positive_for_chunk
    WRITE(kStdWarn,*) 'for ',iEMinus, ' points'
    DO IF = 1,iEminus
      raE(iaEMinus(IF)) = rEmin_positive_for_chunk
      raE(iaEMinus(IF)) = 0.0                  !!before Apr 07
      raE(iaEMinus(IF)) = raBo(iaEMinus(IF))   !!Apr 07
    END DO
  END IF
  IF (iFminus > 0) THEN
    WRITE(kStdWarn,*) 'oops, layer has negative F emiss!! "fixing" this'
    WRITE(kStdWarn,*) 'by setting raF(bad) = ',rFmin_positive_for_chunk
    WRITE(kStdWarn,*) 'for ',iFMinus, ' points'
    DO IF = 1,iFminus
      raF(iaFMinus(IF)) = rFmin_positive_for_chunk
      raF(iaFMinus(IF)) = 0.0                  !!before Apr 07
      raF(iaEMinus(IF)) = raBo(iaEMinus(IF))   !!Apr 07
    END DO
  END IF
END IF

RETURN
END SUBROUTINE asymcoeffs_solar

!************************************************************************
! this subroutine computes upwards R,T,Emiss for arbitrary view angle
! same as /home/sergio/MATLAB/RADTrans/UP/t_r_e_arb_up_Solar_prof.m

SUBROUTINE t_r_e_arb_up_Solar_prof(
! input params  &
raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
    raTau,muSat,muSun,raAsym,raW0,raCumSum,radSolarCld,raBo,  &
    iPhase,raPhasePoints,raComputedPhase,
! output params  &
tr_up,re_up,emiss_up,sun_up)

IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input parameters
REAL :: muSat,muSun,raTau(kMaxPts)   !view,solar angle and extinction
REAL :: raCumSum(kMaxPts),radSolarCld(kMaxPts)!direct solar extinction
REAL :: raAsym(kMaxPts),raW0(kMaxPts)   !asymmetry and albedo
REAL :: raDet(kMaxPts)                  !determinant
REAL :: raE(kMaxPts),raF(kMaxPts)       !used for emission
REAL :: raG(kMaxPts),raH(kMaxPts)       !used for sun
REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff
REAL :: raBeta(kMaxPts),raBo(kMaxPts)   !temp of lower level, variation
INTEGER :: iPhase                  !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! output parameters
REAL :: tr_up(kMaxPts),re_up(kMaxPts),emiss_up(kMaxPts),sun_up(kMaxPts)

! local variables
INTEGER :: IF
REAL :: mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
REAL :: rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
REAL :: eup,fup,gup,hup,sup,hg2,betaa,mu,avgmin

INTEGER :: iArbHG,iHG

avgmin = kAvgMin

mu_pm  = 1/SQRT(3.0)
mu     = ABS(muSat)
muSun  = ABS(muSun)

CALL FindIHG(avg,iHG,iArbHG,raW0)

IF (avg >= avgmin) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    rSun = radSolarCld(IF)*EXP(-raCumSum(IF)/muSun)
    det = raDet(IF)
    
    b = (1-raAsym(IF))/2
    alpha = (raW0(IF)*(1-b)-1)
    
    alpha_p = raW0(IF)/2/mu * (1 + 3.0*raAsym(IF)*mu_pm*mu)
    alpha_m = raW0(IF)/2/mu * (1 - 3.0*raAsym(IF)*mu_pm*mu)
    
    vp = 1/mu + raKp(IF)
    vm = 1/mu + raKm(IF)
    vb = betaa + 1/mu
    vs = 1/muSun + 1/mu
    
!evaluate at raTau(iF)
    evm = EXP(vm*raTau(IF)) - 1
    evp = EXP(vp*raTau(IF)) - 1
    evb = EXP(vb*raTau(IF)) - 1
    evs = (EXP(vs*raTau(IF))-1)/vs
    
    omega = evb/vb
    gamma_a = (alpha_p * raAp(IF) + alpha_m) * evp/vp
    gamma_b = (alpha_p * raAm(IF) + alpha_m) * evm/vm
    
    kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(IF)*b)**2
    kappa = (betaa*mu_pm + alpha - raW0(IF)*b)/kappa
    
    tr_up(IF) = (gamma_a*EXP(raKm(IF)*raTau(IF)) -  &
        gamma_b*EXP(raKp(IF)*raTau(IF)))/det
    re_up(IF) = (raAp(IF) *gamma_b - raAm(IF)*gamma_a)/det
    
    eup = -tr_up(IF) + omega*(alpha_p + 1/(mu*kappa))
    fup = -re_up(IF) + omega*alpha_m
    emiss_up(IF) = raE(IF)*eup + raF(IF)*fup
    
    IF (kSolar >= 0) THEN
      gup = EXP(-raTau(IF)/muSun)*(evs*alpha_p - tr_up(IF))
      hup = EXP(-raTau(IF)/muSun)*evs*alpha_m - re_up(IF)
      IF (iarbhg > 0) THEN
        sup = EXP(-raTau(IF)/muSun)*raW0(IF)*evs/kForP/mu *  &
            hg2(-ABS(muSun),ABS(mu),raAsym(IF),  &
            iPhase,raPhasePoints,raComputedPhase)
      ELSE
        sup = EXP(-raTau(IF)/muSun)*raW0(IF)*evs/kForP/mu *  &
            (1-3.0*raAsym(IF)*mu*muSun)
      END IF
      sun_up(IF) = (raG(IF)*gup + raH(IF)*hup + rSun*sup)/rSun
!sun            sun_up(iF) = sup
    ELSE
      sun_up(IF) = 0.0
    END IF
    
  END DO
  
ELSE IF ((avg < avgmin) .AND. (kTemperVary < 0)) THEN
  DO IF = 1,kMaxPts
!!!!!!!!!no scattering, clear sky solns
    betaa = raBeta(IF)
    
! evaluate at tau0
    tr_up(IF) = 0.0        !remember, no 2stream radiation knocked here
    re_up(IF) = 0.0        !remember, no 2stream radiation knocked here
    
    eup = raBo(IF)*(1 - EXP(-raTau(IF)/mu))
    eup = eup*EXP(raTau(IF)/mu)
    fup = 0.0
    gup = 0.0          !remember, no scattering of solar beam to here
    hup = 0.0          !remember, no scattering of solar beam to here
    
    emiss_up(IF) = eup
    sun_up(IF)   = 0.0
  END DO
  
ELSE IF ((avg < avgmin) .AND. (kTemperVary > 0)) THEN
  DO IF = 1,kMaxPts
!!!!!!!!!no scattering, clear sky solns
    betaa = raBeta(IF)
    
! evaluate at tau0
    tr_up(IF) = 0.0        !remember, no 2stream radiation knocked here
    re_up(IF) = 0.0        !remember, no 2stream radiation knocked here
    
    eup = raBo(IF)/(1+betaa*mu)* (EXP(betaa*raTau(IF))-EXP(-raTau(IF)/mu))
    eup = eup*EXP(raTau(IF)/mu)
    
    fup = 0.0
    gup = 0.0          !remember, no scattering of solar beam to here
    hup = 0.0          !remember, no scattering of solar beam to here
    
    emiss_up(IF) = eup
    sun_up(IF)   = 0.0
  END DO
END IF

RETURN
END SUBROUTINE t_r_e_arb_up_Solar_prof

!************************************************************************
! this subroutine computes downwards R,T,Emiss for arbitrary view angle
! same as /home/sergio/MATLAB/RADTrans/DOWN/t_r_e_arb_down_Solar_prof.m

SUBROUTINE t_r_e_arb_down_Solar_prof(  &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, radSolarCld,raBo,  &
    iPhase,raPhasePoints,raComputedPhase, tr_dn,re_dn,emiss_dn,sun_dn)


REAL, INTENT(IN)                         :: raDet(kMaxPts)
REAL, INTENT(IN)                         :: raE(kMaxPts)
REAL, INTENT(IN)                         :: raF(kMaxPts)
REAL, INTENT(IN OUT)                     :: raG(kMaxPts)
REAL, INTENT(IN)                         :: raH(kMaxPts)
REAL, INTENT(IN)                         :: raKp(kMaxPts)
REAL, INTENT(IN)                         :: raKm(kMaxPts)
REAL, INTENT(IN)                         :: raAp(kMaxPts)
REAL, INTENT(IN)                         :: raAm(kMaxPts)
REAL, INTENT(IN)                         :: raBeta(kMaxPts)
REAL, INTENT(IN)                         :: raTau(kMaxPts)
REAL, INTENT(IN OUT)                     :: muSat
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(IN)                         :: raAsym(kMaxPts)
REAL, INTENT(IN)                         :: raW0(kMaxPts)
REAL, INTENT(IN)                         :: raCumSum(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: radSolarCl
REAL, INTENT(IN)                         :: raBo(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
REAL, INTENT(OUT)                        :: tr_dn(kMaxPts)
REAL, INTENT(OUT)                        :: re_dn(kMaxPts)
REAL, INTENT(OUT)                        :: emiss_dn(kMaxPts)
REAL, INTENT(OUT)                        :: sun_dn(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input parameters
REAL :: !view,solar angle and extinction
REAL :: radSolarCld(kMaxPts)!direct solar extinction
REAL :: !asymmetry and albedo
REAL :: !determinant
REAL :: !used for emission
REAL :: !used for sun
REAL :: !eigen stuff
REAL :: !temp of lower level, variation
INTEGER :: iPhase                  !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! output parameters


! local variables
INTEGER :: IF
REAL :: mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
REAL :: rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
REAL :: edn,fdn,gdn,hdn,sdn,hg2,betaa,mu,avgmin
INTEGER :: iArbHG,iHG

avgmin = kAvgMin

mu_pm  = 1/SQRT(3.0)
mu     = ABS(muSat)
muSun  = ABS(muSun)

CALL FindIHG(avg,iHG,iArbHG,raW0)

IF (avg >= avgmin) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    rSun = radSolarCld(IF)*EXP(-raCumSum(IF)/muSun)
    det = raDet(IF)
    
    b = (1-raAsym(IF))/2
    alpha = (raW0(IF)*(1-b)-1)
    
    alpha_p = raW0(IF)/2/mu * (1 - 3.0*raAsym(IF)*mu_pm*mu)
    alpha_m = raW0(IF)/2/mu * (1 + 3.0*raAsym(IF)*mu_pm*mu)
    
    vp = -1/mu + raKp(IF)
    vm = -1/mu + raKm(IF)
    vb = betaa - 1/mu
    vs = 1/muSun - 1/mu
    
!evaluate at 0.0
    evm = EXP(vm*raTau(IF)) - 1
    evp = EXP(vp*raTau(IF)) - 1
    evb = EXP(vb*raTau(IF)) - 1
    evs = (EXP(vs*raTau(IF))-1)/vs
    
    omega = evb/vb
    gamma_a = (alpha_p * raAp(IF) + alpha_m) * evp/vp
    gamma_b = (alpha_p * raAm(IF) + alpha_m) * evm/vm
    
    kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(IF)*b)**2
    kappa = (betaa*mu_pm + alpha - raW0(IF)*b)/kappa
    
    re_dn(IF) = (gamma_a*EXP(raKm(IF)*raTau(IF)) -  &
        gamma_b*EXP(raKp(IF)*raTau(IF)))/det
    tr_dn(IF) = (raAp(IF) *gamma_b - raAm(IF)*gamma_a)/det
    
    edn = -re_dn(IF) + omega*(alpha_p + 1/(mu*kappa))
    fdn = -tr_dn(IF)*EXP(betaa*raTau(IF)) + omega*alpha_m
    emiss_dn(IF) = raE(IF)*edn + raF(IF)*fdn
    
    IF (kSolar >= 0) THEN
      gdn = EXP(-raTau(IF)/muSun)*(evs*alpha_p - re_dn(IF))
      hdn = EXP(-raTau(IF)/muSun)*evs*alpha_m - tr_dn(IF)
      
      IF (iarbhg > 0) THEN
        sdn = EXP(-raTau(IF)/muSun)*raW0(IF)*evs/kForP/mu *  &
            hg2(-ABS(muSun),-ABS(mu),raAsym(IF),  &
            iPhase,raPhasePoints,raComputedPhase)
      ELSE
        sdn = EXP(-raTau(IF)/muSun)*raW0(IF)*evs/kForP/mu *  &
            (1+3.0*raAsym(IF)*mu*muSun)
      END IF
      sun_dn(IF) = (raG(IF)*gdn + raH(IF)*hdn + rSun*sdn)/rSun
      
    ELSE
      sun_dn(IF) = 0.0
    END IF
    
  END DO
  
ELSE IF ((avg < avgmin) .AND. (kTemperVary < 0)) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
!!!!!!!!!no scattering, clear sky solns
    
! evaluate at tau0
    tr_dn(IF) = 0.0        !remember, no 2stream radiation knocked here
    re_dn(IF) = 0.0        !remember, no 2stream radiation knocked here
    
    edn = 0.0
    fdn = raBo(IF)*(1-EXP(-raTau(IF)/mu))
    
    gdn = 0.0          !remember, no scattering of solar beam to here
    hdn = 0.0          !remember, no scattering of solar beam to here
    
    emiss_dn(IF) = fdn
    sun_dn(IF)   = 0.0
  END DO
  
ELSE IF ((avg < avgmin) .AND. (kTemperVary > 0)) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
!!!!!!!!!no scattering, clear sky solns
    
! evaluate at tau0
    tr_dn(IF) = 0.0        !remember, no 2stream radiation knocked here
    re_dn(IF) = 0.0        !remember, no 2stream radiation knocked here
    
    edn = 0.0
    fdn = (1-EXP(-raTau(IF)/mu))*EXP(betaa*raTau(IF))
    fdn = raBo(IF)/(1 - betaa*mu) * fdn
    
    gdn = 0.0          !remember, no scattering of solar beam to here
    hdn = 0.0          !remember, no scattering of solar beam to here
    
    emiss_dn(IF) = fdn
    sun_dn(IF)   = 0.0
  END DO
END IF

DO IF = 1,kMaxPts
  tr_dn(IF)    = tr_dn(IF)    * EXP(raTau(IF)/mu)
  re_dn(IF)    = re_dn(IF)    * EXP(raTau(IF)/mu)
  emiss_dn(IF) = emiss_dn(IF) * EXP(raTau(IF)/mu)
  sun_dn(IF)   = sun_dn(IF)   * EXP(raTau(IF)/mu)
END DO

RETURN
END SUBROUTINE t_r_e_arb_down_Solar_prof

!************************************************************************
! this is for the twostream radiation
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/t_r_e_streams_Solar_prof.m

SUBROUTINE t_r_e_streams_Solar_prof(  &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
    raTau1,muSat,muSun,raAsym,raW0,raCumSum1,radSolarCld,raTb,  &
    raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)


REAL, INTENT(IN)                         :: raDet(kMaxPts)
REAL, INTENT(IN)                         :: raE(kMaxPts)
REAL, INTENT(IN)                         :: raF(kMaxPts)
REAL, INTENT(IN OUT)                     :: raG(kMaxPts)
REAL, INTENT(IN)                         :: raH(kMaxPts)
REAL, INTENT(IN)                         :: raKp(kMaxPts)
REAL, INTENT(IN OUT)                     :: raKm(kMaxPts)
REAL, INTENT(IN)                         :: raAp(kMaxPts)
REAL, INTENT(IN)                         :: raAm(kMaxPts)
REAL, INTENT(IN)                         :: raBeta(kMaxPts)
REAL, INTENT(IN)                         :: raTau1(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: muSat
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(IN OUT)                     :: raAsym(kMaxPts)
REAL, INTENT(IN OUT)                     :: raW0(kMaxPts)
REAL, INTENT(IN)                         :: raCumSum1(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: radSolarCl
REAL, INTENT(IN OUT)                     :: raTb(kMaxPts)
REAL, INTENT(OUT)                        :: raT1(kMaxPts)
REAL, INTENT(OUT)                        :: raT1star(kMaxPts)
REAL, INTENT(OUT)                        :: raR1(kMaxPts)
REAL, INTENT(OUT)                        :: raR1star(kMaxPts)
REAL, INTENT(OUT)                        :: raE1p(kMaxPts)
REAL, INTENT(OUT)                        :: raE1m(kMaxPts)
REAL, INTENT(OUT)                        :: raS1p(kMaxPts)
REAL, INTENT(OUT)                        :: raS1m(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! output params
REAL :: !transmissions
REAL :: !reflections
REAL :: !emissions
REAL :: !sun stuff

! input params
REAL :: muSat                     !solar,view angle (not needed)
REAL :: !asymmetry
REAL :: radSolarCld(kMaxPts), !solar stuff
REAL :: !radiance at bottom
REAL :: !albedo and temperature variation
REAL :: !eigenvalues
REAL :: !eigenvectors
REAL :: !optical depth, determinant
REAL :: !up and down emission
REAL :: !up and down sun stuff


! local variables
INTEGER :: IF,iHG,iArbHG
REAL :: avg,det,bw1,bw2,hg2,rSun,tau,mu_pm,tau0
REAL :: tr_up,re_up,eup,fup,gup,hup,sun_up,emiss_up
REAL :: tr_down,re_down,edn,fdn,gdn,hdn,sun_down,emiss_down,betaa
REAL :: avgmin

avgmin = kAvgMin

CALL FindIHG(avg,iHG,iArbHG,raW0)

mu_pm   = 1/SQRT(3.0)
muSun  = ABS(muSun)

IF (avg > avgmin) THEN
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    rSun = radSolarCld(IF)*EXP(-raCumSum1(IF)/muSun)
    
    tau  = raTau1(IF)
    tau0 = raTau1(IF)
    tr_up = raAp(IF)*EXP(raKm(IF)*tau0)*EXP(raKp(IF)*tau) -  &
        raAm(IF)*EXP(raKp(IF)*tau0)*EXP(raKm(IF)*tau)
    tr_up = tr_up/raDet(IF)
    re_up = raAp(IF)*raAm(IF)*(EXP(raKm(IF)*tau) - EXP(raKp(IF)*tau))
    re_up = re_up/raDet(IF)
    eup = EXP(betaa*tau) - tr_up
    fup = -re_up*EXP(betaa*tau0)
    gup = EXP(-tau0/muSun)*(EXP(tau/muSun) - tr_up)
    hup = -re_up
    emiss_up = raE(IF)*eup + raF(IF)*fup
    IF (kSolar >= 0) THEN
      sun_up = (raG(IF)*gup+ raH(IF)*hup)/rSun
    ELSE
      sun_up = 0.0
    END IF
    raT1(IF) = tr_up
    raR1(IF) = re_up
    raE1p(IF) = emiss_up
    raS1p(IF) = sun_up
    
    tau = 0.0
    tau0 = raTau1(IF)
    tr_down = raAp(IF)*EXP(raKm(IF)*tau)-raAm(IF)*EXP(raKp(IF)*tau)
    tr_down = tr_down/raDet(IF)
    re_down = EXP(raKm(IF)*tau0)*EXP(raKp(IF)*tau) -  &
        EXP(raKp(IF)*tau0)*EXP(raKm(IF)*tau)
    re_down = re_down/raDet(IF)
    edn = -re_down
    fdn = -tr_down*EXP(betaa*tau0) + EXP(betaa*tau)
    emiss_down = raE(IF)*edn + raF(IF)*fdn
    IF (kSolar >= 0) THEN
      gdn = EXP(-tau0/muSun)*(-re_down)
      hdn = EXP(-(tau0-tau)/muSun) - tr_down
      sun_down = (raG(IF)*gdn + raH(IF)*hdn)/rSun
    ELSE
      sun_down = 0.0
    END IF
    
    raT1star(IF) = tr_down
    raR1star(IF) = re_down
    raE1m(IF) = emiss_down
    raS1m(IF) = sun_down
    
  END DO
  
ELSE
  DO IF = 1,kMaxPts
    betaa = raBeta(IF)
    
    tau  = raTau1(IF)
    tau0 = raTau1(IF)
    
    tr_up = EXP(-raKp(IF)*tau)
    re_up = 0.0
    eup = EXP(betaa*tau) - tr_up
    fup = 0.0
    gup = 0.0
    hup = 0.0
    emiss_up = raE(IF)*eup + raF(IF)*fup
    sun_up = 0.0
    raT1(IF) = tr_up
    raR1(IF) = re_up
    raE1p(IF) = emiss_up
    raS1p(IF) = sun_up
    
    tau = 0.0
    tr_down = EXP(-raKp(IF)*(tau0-tau))
    re_down = 0.0
    edn = 0.0
    fdn = 1 - EXP(-tau0/mu_pm)*EXP(betaa*tau0)
    gdn = 0.0
    hdn = 0.0
    emiss_down = raE(IF)*edn + raF(IF)*fdn
    sun_down   = 0.0
    raT1star(IF) = tr_down
    raR1star(IF) = re_down
    raE1m(IF) = emiss_down
    raS1m(IF) = sun_down
    
  END DO
END IF

RETURN
END SUBROUTINE t_r_e_streams_Solar_prof

!************************************************************************
! this subroutine adds together stuff for the arbitrary angles
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/addstar_arb_solar.m

SUBROUTINE  addstar_arb_solar(  &
    raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1,  &
    raSunUp1,raSunDown1,raTau1,  &
    raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2,  &
    raSunUp2,raSunDown2,raTau2,  &
    raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m,  &
    raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,  &
    raReUp12,raReDown12,raTrUp12,raTrDown12,  &
    raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12,  &
    raCumSum1,raCumSum2,muSat,muSun)


REAL, INTENT(IN)                         :: raReUp1(kMaxPts)
REAL, INTENT(IN)                         :: raReDown1(kMaxPts)
REAL, INTENT(IN)                         :: raTrUp1(kMaxPts)
REAL, INTENT(IN)                         :: raTrDown1(kMaxPts)
REAL, INTENT(IN)                         :: raEmissUp1(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
REAL, INTENT(IN)                         :: raSunUp1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunDown1(kMaxPts)
NO TYPE, INTENT(IN)                      :: raTau1
REAL, INTENT(IN)                         :: raReUp2(kMaxPts)
REAL, INTENT(IN)                         :: raReDown2(kMaxPts)
REAL, INTENT(IN)                         :: raTrUp2(kMaxPts)
REAL, INTENT(IN)                         :: raTrDown2(kMaxPts)
REAL, INTENT(IN)                         :: raEmissUp2(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
REAL, INTENT(IN OUT)                     :: raSunUp2(kMaxPts)
REAL, INTENT(IN)                         :: raSunDown2(kMaxPts)
NO TYPE, INTENT(IN)                      :: raTau2
REAL, INTENT(IN)                         :: raR1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raT1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raR1star(kMaxPts)
REAL, INTENT(IN)                         :: raT1star(kMaxPts)
REAL, INTENT(IN OUT)                     :: raE1p(kMaxPts)
REAL, INTENT(IN OUT)                     :: raE1m(kMaxPts)
REAL, INTENT(IN OUT)                     :: raS1p(kMaxPts)
REAL, INTENT(IN OUT)                     :: raS1m(kMaxPts)
REAL, INTENT(IN OUT)                     :: raR2(kMaxPts)
REAL, INTENT(IN)                         :: raT2(kMaxPts)
REAL, INTENT(IN)                         :: raR2star(kMaxPts)
REAL, INTENT(IN OUT)                     :: raT2star(kMaxPts)
REAL, INTENT(IN OUT)                     :: raE2p(kMaxPts)
REAL, INTENT(IN OUT)                     :: raE2m(kMaxPts)
REAL, INTENT(IN OUT)                     :: raS2p(kMaxPts)
REAL, INTENT(IN OUT)                     :: raS2m(kMaxPts)
REAL, INTENT(OUT)                        :: raReUp12(kMaxPts)
REAL, INTENT(OUT)                        :: raReDown12(kMaxPts)
REAL, INTENT(OUT)                        :: raTrUp12(kMaxPts)
REAL, INTENT(OUT)                        :: raTrDown12(kMaxPts)
NO TYPE, INTENT(IN)                      :: raEmissUp1
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
REAL, INTENT(OUT)                        :: raSunUp12(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSunDown1
REAL, INTENT(OUT)                        :: raTau12(kMaxPts)
REAL, INTENT(IN OUT)                     :: raCumSum1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raCumSum2(kMaxPts)
REAL, INTENT(IN OUT)                     :: muSat
REAL, INTENT(OUT)                        :: muSun
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
!! layer stuff for the arb angles, layer n


REAL :: raEmissDown1(kMaxPts),ratau1(kMaxPts)

!! layer stuff for the arb angles, layer n+1


REAL :: raEmissDown2(kMaxPts),ratau2(kMaxPts)

!! layer stuff for the stream angles, layer n




!! layer stuff for the stream angles, layer n+1




!!arbitrary stuff


! output params
!! layer stuff for the arb angles, layer n AND n+1


REAL :: raEmissUp12(kMaxPts),raEmissDown12(kMaxPts)
REAL :: raSunDown12(kMaxPts)

! local variables
INTEGER :: IF
REAL :: mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat11,mat12,mat21,mat22
REAL :: mu,det,a1,a2,tau12,e1,e2

mu = ABS(muSat)
muSun = ABS(muSun)

DO IF = 1,kMaxPts
  
  det = 1-raR1(IF)*raR2star(IF)
  
!!!!!!!reflection and transmission
  mat1 = 0.0
  mat2 = 0.0
  mat3 = raReUp2(IF)*EXP(raTau1(IF)/mu)
  mat4 = 0.0
  mat5 = raT2(IF)*raReUp1(IF)/det
  mat6 = 0.0
  mat7 = raR1(IF)*raT2(IF)*raTrUp2(IF)/det*EXP(raTau1(IF)/mu)
  mat11 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
  
  mat1 = raTrUp1(IF)
  mat2 = 0.0
  mat3 = 0.0
  mat4 = 0.0
  mat5 = raR2star(IF)*raT1star(IF)*raReUp1(IF)/det
  mat6 = 0.0
  mat7 = raT1star(IF)*raTrUp2(IF)/det*EXP(raTau1(IF)/mu)
  mat12 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
  
  mat1 = raTrDown2(IF)
  mat2 = 0.0
  mat3 = 0.0
  mat4 = raR1(IF)*raT2(IF)*raReDown2(IF)/det
  mat5 = 0.0
  mat6 = raT2(IF)*raTrDown1(IF)/det*EXP(raTau2(IF)/mu)
  mat7 = 0.0
  mat21 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
  
  mat1 = 0.0
  mat2 = raReDown1(IF)*EXP(raTau2(IF)/mu)
  mat3 = 0.0
  mat4 = raT1star(IF)*raReDown2(IF)/det
  mat5 = 0.0
  mat6 = raR2star(IF)*raT1star(IF)*raTrDown1(IF)/det*EXP(raTau2(IF)/mu)
  mat7 = 0.0
  mat22 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
  
!!!!emission
  a1 = raEmissUp1(IF) + 0.0 + EXP(raTau1(IF)/mu)*raEmissUp2(IF) + 0.0 +  &
      raTrUp1(IF)/det*(raR2star(IF)*raE1p(IF)+raE2m(IF)) + 0.0 +  &
      raTrUp2(IF)/det*(raR1(IF)*raE2m(IF)+raE1p(IF))*EXP(raTau1(IF)/mu)
  
  a2 = raEmissDown2(IF) + EXP(raTau2(IF)/mu)*raEmissDown1(IF) + 0.0 +  &
      raReDown2(IF)/det*(raR1(IF)*raE2m(IF) + raE1p(IF)) + 0.0 +  &
      raTrDown1(IF)/det*(raR2star(IF)*raE1p(IF)+raE2m(IF))*EXP(raTau2(IF)/mu)  &
      + 0.0
  
  raReUp12(IF)      = mat11
  raTrUp12(IF)      = mat12
  raTrDown12(IF)    = mat21
  raReDown12(IF)    = mat22
  raEmissUp12(IF)   = a1
  raEmissDown12(IF) = a2
  
  raTau12(IF)  = raTau1(IF) + raTau2(IF)     !total optical depth so far
  
!        print *,raEmissUp1(iF),raEmissUp2(iF),raTrUp1(iF),
!     $         raTrUp2(iF),raE1p(iF),raE2m(iF),raR1(iF),raR2star(iF),det
  
END DO

IF (kSolar >= 0) THEN           !!!do the solar part
  DO IF = 1,kMaxPts
    det = 1 - raR1(IF)*raR2star(IF)
    e1 = EXP(-raCumsum1(IF)/muSun)
    e2 = EXP(-raCumsum2(IF)/muSun)
    
    raSunUp12(IF) = raSunUp1(IF)*e1 + raSunUp2(IF)*e2*EXP(raTau1(IF)/mu) +  &
        1/det*(raReUp1(IF)*(e1*raR2star(IF)*raS1p(IF) + e2*raS2m(IF))) +  &
        EXP(raTau1(IF)/mu)*raTrUp2(IF)/det *  &
        (e2*raR1(IF)*raS2m(IF) + e1*raS1p(IF))
    raSunDown12(IF) = raSunDown2(IF)*e2 +  &
        raSunDown1(IF)*e1*EXP(raTau2(IF)/mu) +  &
        1/det*(raReDown2(IF)*(e2*raR1(IF)*raS2m(IF) + e1*raS1p(IF))) +  &
        EXP(raTau2(IF)/mu)*raTrDown1(IF)/det *  &
        (e1*raR2star(IF)*raS1p(IF) + e2*raS2m(IF))
  END DO
END IF

RETURN
END SUBROUTINE  addstar_arb_solar

!************************************************************************
! this subroutine adds together the layers for the stream angles
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/addstar_solar.m

SUBROUTINE addstar_solar(  &
    raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m,  &
    raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m,  &
    raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,  &
    raCumSum1,raCumSum2,muSun)


REAL, INTENT(OUT)                        :: raR12(kMaxPts)
REAL, INTENT(OUT)                        :: raT12(kMaxPts)
REAL, INTENT(OUT)                        :: raR12star(kMaxPts)
REAL, INTENT(OUT)                        :: raT12star(kMaxPts)
REAL, INTENT(OUT)                        :: raE12p(kMaxPts)
REAL, INTENT(OUT)                        :: raE12m(kMaxPts)
REAL, INTENT(OUT)                        :: raS12p(kMaxPts)
REAL, INTENT(OUT)                        :: raS12m(kMaxPts)
REAL, INTENT(IN)                         :: raR1(kMaxPts)
REAL, INTENT(IN)                         :: raT1(kMaxPts)
REAL, INTENT(IN)                         :: raR1star(kMaxPts)
REAL, INTENT(IN)                         :: raT1star(kMaxPts)
REAL, INTENT(IN)                         :: raE1p(kMaxPts)
REAL, INTENT(IN)                         :: raE1m(kMaxPts)
REAL, INTENT(IN)                         :: raS1p(kMaxPts)
REAL, INTENT(IN)                         :: raS1m(kMaxPts)
REAL, INTENT(IN)                         :: raR2(kMaxPts)
REAL, INTENT(IN)                         :: raT2(kMaxPts)
REAL, INTENT(IN)                         :: raR2star(kMaxPts)
REAL, INTENT(IN)                         :: raT2star(kMaxPts)
REAL, INTENT(IN)                         :: raE2p(kMaxPts)
REAL, INTENT(IN)                         :: raE2m(kMaxPts)
REAL, INTENT(IN)                         :: raS2p(kMaxPts)
REAL, INTENT(IN)                         :: raS2m(kMaxPts)
REAL, INTENT(IN)                         :: raCumSum1(kMaxPts)
REAL, INTENT(IN)                         :: raCumSum2(kMaxPts)
REAL, INTENT(OUT)                        :: muSun
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params
!! layer stuff for the stream angles, layer n



!! layer stuff for the stream angles, layer n+1



! for the sun

! output params
!! layer stuff for the stream angles, layer n AND n+1



! local variables
INTEGER :: IF
REAL :: det,e1,e2

muSun = ABS(muSun)

DO IF = 1,kMaxPts
  det = 1 - raR1(IF)*raR2star(IF)
  
  raR12(IF) = raT2star(IF)*raR1(IF)/det
  raT12star(IF) = raT2star(IF)/det
  raT12(IF) = raT1(IF)/det
  raR12star(IF) = raT1(IF)*raR2star(IF)/det
  
  raE12p(IF) = raE2p(IF) + (raR12(IF)*raE2m(IF)+raT12star(IF)*raE1p(IF))
  raE12m(IF) = raE1m(IF) + (raT12(IF)*raE2m(IF)+raR12star(IF)*raE1p(IF))
  
  raR12(IF)     = raR12(IF)*raT2(IF) + raR2(IF)
  raT12star(IF) = raT1star(IF)*raT12star(IF)
  raT12(IF)     = raT2(IF)*raT12(IF)
  raR12star(IF) = raR1star(IF) + raT1star(IF)*raR12star(IF)
END DO

IF (kSolar >= 0) THEN  !!!!do the solar part
  DO IF = 1,kMaxPts
    det = 1 - raR1(IF)*raR2star(IF)
    e1 = EXP(-raCumSum1(IF)/muSun)
    e2 = EXP(-raCumSum2(IF)/muSun)
    
    raS12p(IF) = (raT2star(IF)/det*raR1(IF)*raS2m(IF) + raS2p(IF))*e2 +  &
        raT2star(IF)/det*raS1p(IF)*e1
    raS12m(IF) = (raT1(IF)/det*raR2star(IF)*raS1p(IF) + raS1m(IF))*e1 +  &
        raT1(IF)/det*raS2m(IF)*e2
  END DO
END IF

RETURN
END SUBROUTINE addstar_solar

!************************************************************************
! this subroutine updates the stream angle added layers
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/update_streams_solar.m

SUBROUTINE update_streams_solar(  &
    raR2,raR2star,raT2,raT2star,raE2p,raE2m,raS2p,raS2m,  &
    raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)


REAL, INTENT(OUT)                        :: raR2(kMaxPts)
REAL, INTENT(OUT)                        :: raR2star(kMaxPts)
REAL, INTENT(OUT)                        :: raT2(kMaxPts)
REAL, INTENT(OUT)                        :: raT2star(kMaxPts)
REAL, INTENT(OUT)                        :: raE2p(kMaxPts)
REAL, INTENT(OUT)                        :: raE2m(kMaxPts)
REAL, INTENT(OUT)                        :: raS2p(kMaxPts)
REAL, INTENT(OUT)                        :: raS2m(kMaxPts)
REAL, INTENT(IN)                         :: raR12(kMaxPts)
REAL, INTENT(IN)                         :: raR12star(kMaxPts)
REAL, INTENT(IN)                         :: raT12(kMaxPts)
REAL, INTENT(IN)                         :: raT12star(kMaxPts)
REAL, INTENT(IN)                         :: raE12p(kMaxPts)
REAL, INTENT(IN)                         :: raE12m(kMaxPts)
REAL, INTENT(IN)                         :: raS12p(kMaxPts)
REAL, INTENT(IN)                         :: raS12m(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! stream angle stuff









! local variables
INTEGER :: IF

!!!only need to update the upper layer
DO IF = 1,kMaxPts
  raR2(IF)     = raR12(IF)
  raR2star(IF) = raR12star(IF)
  raT2(IF)     = raT12(IF)
  raT2star(IF) = raT12star(IF)
  raE2p(IF)    = raE12p(IF)
  raE2m(IF)    = raE12m(IF)
  raS2p(IF)    = raS12p(IF)
  raS2m(IF)    = raS12m(IF)
END DO

RETURN
END SUBROUTINE update_streams_solar

!************************************************************************
! this subroutine updates the arb angle added layers
! this just sets varNew == varOld
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/update_streams_arb_solar.m

SUBROUTINE update_streams_arb_solar(  &
    raTauNew,raReUpNew,raReDownNew,raTrUpNew,raTrDownNew,  &
    raEmissUpNew,raEmissDownNew,raSunUpNew,raSunDownNew,  &
    raTauOld,raReUpOld,raReDownOld,raTrUpOld,raTrDownOld,  &
    raEmissUpOld,raEmissDownOld,raSunupOld,raSunDownOld)


REAL, INTENT(OUT)                        :: raTauNew(kMaxPts)
REAL, INTENT(OUT)                        :: raReUpNew(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raReDownNe
REAL, INTENT(OUT)                        :: raTrUpNew(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raTrDownNe
NO TYPE, INTENT(IN OUT)                  :: raEmissUpN
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
REAL, INTENT(OUT)                        :: raSunUpNew(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSunDownN
REAL, INTENT(IN)                         :: raTauOld(kMaxPts)
REAL, INTENT(IN)                         :: raReUpOld(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raReDownOl
REAL, INTENT(IN)                         :: raTrUpOld(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raTrDownOl
NO TYPE, INTENT(IN OUT)                  :: raEmissUpO
NO TYPE, INTENT(IN OUT)                  :: raEmissDow
NO TYPE, INTENT(IN OUT)                  :: raSunupOld
NO TYPE, INTENT(IN OUT)                  :: raSunDownO
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! arbitrary angle stuff

REAL :: raEmissUpNew(kMaxPts)
REAL :: raEmissUpOld(kMaxPts)
REAL :: raTrDownNew(kMaxPts),raReDownNew(kMaxPts),raEmissDownNew(kMaxPts)
REAL :: raTrDownOld(kMaxPts),raReDownOld(kMaxPts),raEmissDownOld(kMaxPts)
REAL :: raSunDownNew(kMaxPts)
REAL :: raSunUpOld(kMaxPts),raSunDownOld(kMaxPts)

! local variables
INTEGER :: IF

!!!only need to update the upper layer
DO IF = 1,kMaxPts
  raTauNew(IF)       = raTauOld(IF)
  raReUpNew(IF)      = raReUpOld(IF)
  raReDownNew(IF)    = raReDownOld(IF)
  raTrUpNew(IF)      = raTrUpOld(IF)
  raTrDownNew(IF)    = raTrDownOld(IF)
  raEmissUpNew(IF)   = raEmissUpOld(IF)
  raEmissDownNew(IF) = raEmissDownOld(IF)
  raSunUpNew(IF)     = raSunUpOld(IF)
  raSunDownNew(IF)   = raSunDownOld(IF)
END DO

RETURN
END SUBROUTINE update_streams_arb_solar

!************************************************************************
! this subroutine sets the scattering properties of the current layer

SUBROUTINE LayerScatteringProp(rFrac,  &
    raFreq,rMPTemp,iL,raaExt,raaScat,raaAsym, rBeta,raBeta_RadOrTemp,
! output params  &
raW0,raAsym,raTau,raBo,raBeta)


REAL, INTENT(IN)                         :: rFrac
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rMPTemp
INTEGER, INTENT(IN OUT)                  :: iL
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaScat(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAsym(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: rBeta
NO TYPE, INTENT(IN OUT)                  :: raBeta_Rad
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters






REAL :: rBeta                 !log(T(n+1)/T(n)) -->
!log(B(iF,T(n+1))/B(iF,T(n)))
REAL :: raBeta_RadOrTemp(kMaxPts)  !log(T(n+1)/T(n)) -->
!log(B(iF,T(n+1))/B(iF,T(n)))
! output parameters
REAL :: raTau(kMaxPts)        !total layer optical depth (includes scatter)
REAL :: raW0(kMaxPts),raAsym(kMaxPts)   !layer albedo, asymmetry
REAL :: raBo(kMaxPts)         !lower level radiances

REAL :: raBeta(kMaxPts)       !beta = 1/tau log(T(n+1)/T(n))

! local variables
INTEGER :: iFr
REAL :: ttorad   !!!same as ttorad(v,T) but could be double precision (oh well)

IF ((rFrac >= 0.0) .AND. (rFrac <= 1.0)) THEN
  DO iFr = 1,kMaxPts
    raBo(iFr)       = ttorad(raFreq(iFr),rMPTemp)
    raTau(iFr)      = raaExt(iFr,iL) * rFrac
    raW0(iFr)       = raaScat(iFr,iL)/raaExt(iFr,iL)
    raAsym(iFr)     = raaAsym(iFr,iL)
    
!ccc          !!!used only by kTwoStream, apparently!
!ccc          !!!log(T(i)/T(i+1))
!ccc          raBeta(iFr)     = 1/raTau(iFr) * rBeta
    
!!!used by CHARTS, DISORT and RTSPEC
!!!log(B(T(i),vj)/B(T(i+1),vj))
    raBeta(iFr)     = 1/raTau(iFr) * raBeta_RadOrTemp(iFr)
    
  END DO
ELSE
!!!! impossible
  WRITE(kStdErr,*) 'In subroutine LayerScatteringProp :'
  WRITE(kStdErr,*) 'Code claims fraction for this layer is : ',rFrac
  CALL DoStop
END IF

RETURN
END SUBROUTINE LayerScatteringProp

!************************************************************************
! this accumulates the optical depth from top of cloud to top of layer

SUBROUTINE AccumulateSolarDepth(raCumSum,raX,iDo)


REAL, INTENT(OUT)                        :: raCumSum(kMaxPts)
REAL, INTENT(IN)                         :: raX(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iDo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
INTEGER :: iDo                 !1 to add on raX to raCumSum, -1 if not
REAL :: !total optical depth from cloud top

! local variables
INTEGER :: iFr

IF (iDo > 0) THEN
  DO iFr = 1,kMaxPts
    raCumSum(iFr) = raCumSum(iFr) + raX(iFr)
  END DO
ELSE
  DO iFr = 1,kMaxPts
    raCumSum(iFr) = 0.0
  END DO
END IF


RETURN
END SUBROUTINE AccumulateSolarDepth

!************************************************************************
! this function computes the phase function, assuming cos(phi1-phi2) factor = 1

REAL FUNCTION hg2(mu1,mu2,g,iPhase,raPhasePoints,raComputedPhase)


REAL, INTENT(IN)                         :: mu1
REAL, INTENT(IN)                         :: mu2
NO TYPE, INTENT(IN OUT)                  :: g
NO TYPE, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

REAL :: g       !mu1,mu2 are the two angles, g is the asymmetry
INTEGER :: iPhase       !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

REAL :: normB,mu0,yexact,yexactNew,hg2_real
DOUBLE PRECISION :: hg2_double
INTEGER :: iJ

IF (iPhase < 0) THEN
  hg2 = hg2_double(DBLE(mu1),DBLE(mu2),DBLE(g))
END IF

IF (iPhase > 0) THEN
! ********************* warning .. use rspl if the variables are real
!!!compute mu0 = cos ofangle between the two
  mu0 = mu1*mu2 + SQRT(1-mu1*mu1)*SQRT(1-mu2*mu2)
  CALL dSpl(raPhasePoints,raComputedPhase,MaxPhase,mu0,yexactNew,1)
! ********************* warning .. use rspl if the variables are reals
  hg2 = yexactNew
END IF

RETURN
END FUNCTION hg2

!************************************************************************
! this subroutine finds out the partial fraction for layer : default = +1.0
! for use in  LayerScatteringProp

SUBROUTINE FindTheFrac(iUpDown,iL,iaRadLayer,iNumLayer,  &
    rFracTop,rFracBot,rFrac)

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars

INTEGER, INTENT(IN OUT)                  :: iUpDown
INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(OUT)                        :: rFrac


INTEGER :: !!layering of current atm


! output vars


rFrac = 1.0  !!default multiplier
IF (iUpDown > 0) THEN    !!!!down look instr : radiation going up
  IF (iL == iaRadLayer(1)) THEN
    rFrac = rFracBot
  ELSE IF (iL == iaRadLayer(iNumLayer)) THEN
    rFrac = rFracTop
  END IF
ELSE IF (iUpDown < 0) THEN    !!!!uplook instr : radiation going down
  IF (iL == iaRadLayer(1)) THEN
    rFrac = rFracTop
  ELSE IF (iL == iaRadLayer(iNumLayer)) THEN
    rFrac = rFracBot
  END IF
END IF

RETURN
END SUBROUTINE FindTheFrac

!************************************************************************
! just a printer

SUBROUTINE PrintStuffTwoStream(iLay,iL,  &
    raAsym,raW0,raTb,raRadBb,raRadBt,raTau2,  &
    raBeta,raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  &
    raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm,  &
    raTrUp2,raReUp2,raEmissUp2,raSunUp2)


INTEGER, INTENT(OUT)                     :: iLay
INTEGER, INTENT(OUT)                     :: iL
REAL, INTENT(IN OUT)                     :: raAsym(kMaxPts)
REAL, INTENT(IN OUT)                     :: raW0(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTb(kMaxPts)
REAL, INTENT(IN OUT)                     :: raRadBb(kMaxPts)
REAL, INTENT(IN OUT)                     :: raRadBt(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTau2(kMaxPts)
REAL, INTENT(IN OUT)                     :: raBeta(kMaxPts)
REAL, INTENT(OUT)                        :: raA(kMaxPts)
REAL, INTENT(OUT)                        :: raB(kMaxPts)
REAL, INTENT(IN OUT)                     :: raDet(kMaxPts)
REAL, INTENT(OUT)                        :: raE(kMaxPts)
REAL, INTENT(OUT)                        :: raF(kMaxPts)
REAL, INTENT(OUT)                        :: raG(kMaxPts)
REAL, INTENT(OUT)                        :: raH(kMaxPts)
REAL, INTENT(OUT)                        :: raKp(kMaxPts)
REAL, INTENT(OUT)                        :: raKm(kMaxPts)
REAL, INTENT(OUT)                        :: raAp(kMaxPts)
REAL, INTENT(OUT)                        :: raAm(kMaxPts)
REAL, INTENT(OUT)                        :: raT(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTstar(kMaxPts)
REAL, INTENT(OUT)                        :: raR(kMaxPts)
REAL, INTENT(IN OUT)                     :: raRstar(kMaxPts)
REAL, INTENT(OUT)                        :: raEp(kMaxPts)
REAL, INTENT(OUT)                        :: raEm(kMaxPts)
REAL, INTENT(OUT)                        :: raSp(kMaxPts)
REAL, INTENT(OUT)                        :: raSm(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTrUp2(kMaxPts)
REAL, INTENT(IN OUT)                     :: raReUp2(kMaxPts)
REAL, INTENT(IN OUT)                     :: raEmissUp2(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunUp2(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
REAL :: raFreq(kMaxPts)                 !wavenumbers
REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud


REAL :: raTt(kMaxPts)
REAL :: raSun0(kMaxPts)







REAL :: raFp(kMaxPts),raFm(kMaxPts)


! arbitrary angle stuff
REAL :: raSunDown2(kMaxPts)

REAL :: raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)




INTEGER :: iFr,i1,i2

i1 = 1
i2 = 2

i1 = 9223
i2 = 9223

i1 = 9222
i2 = 9224

PRINT *,'          '
PRINT *,' ---->   iLay,iL = ',iLay,iL
PRINT *,'from Asymcoeffs_solar ....'
PRINT *,'  asym solar Asym   ',(raAsym(iFr),iFr=i1,i2)
PRINT *,'  asym solar W0     ',(raW0(iFr),iFr=i1,i2)
PRINT *,'  asym solar Tb     ',(raTb(iFr),iFr=i1,i2)
PRINT *,'  asym solar Radbot ',(raRadBb(iFr),iFr=i1,i2)
PRINT *,'  asym solar Radtop ',(raRadBt(iFr),iFr=i1,i2)
PRINT *,'  asym solar Tau2   ',(raTau2(iFr),iFr=i1,i2)
PRINT *,'  asym solar Beta   ',(raBeta(iFr),iFr=i1,i2)
PRINT *,'  asym solar raDET  ',(raDet(iFr),iFr=i1,i2)
PRINT *,'  asym solar raA    ',(raA(iFr),iFr=i1,i2)
PRINT *,'  asym solar raB    ',(raB(iFr),iFr=i1,i2)
PRINT *,'  asym solar raE    ',(raE(iFr),iFr=i1,i2)
PRINT *,'  asym solar raF    ',(raF(iFr),iFr=i1,i2)
PRINT *,'  asym solar raG    ',(raG(iFr),iFr=i1,i2)
PRINT *,'  asym solar raH    ',(raH(iFr),iFr=i1,i2)
PRINT *,'  asym solar raKp   ',(raKp(iFr),iFr=i1,i2)
PRINT *,'  asym solar raKm   ',(raKm(iFr),iFr=i1,i2)
PRINT *,'  asym solar raAp   ',(raAp(iFr),iFr=i1,i2)
PRINT *,'  asym solar raAm   ',(raAm(iFr),iFr=i1,i2)
PRINT *,'from t_r_e_streams_Solar_prof'
PRINT *,'  asym solar raT    ',(raT(iFr),iFr=i1,i2)
PRINT *,'  asym solar raTstar',(raTstar(iFr),iFr=i1,i2)
PRINT *,'  asym solar raR    ',(raR(iFr),iFr=i1,i2)
PRINT *,'  asym solar raRstar',(raRstar(iFr),iFr=i1,i2)
PRINT *,'  asym solar raEp   ',(raEp(iFr),iFr=i1,i2)
PRINT *,'  asym solar raEm   ',(raEm(iFr),iFr=i1,i2)
PRINT *,'  asym solar raSp   ',(raSp(iFr),iFr=i1,i2)
PRINT *,'  asym solar raSm   ',(raSm(iFr),iFr=i1,i2)
PRINT *,'from  t_r_e_arb_up_Solar_prof'
PRINT *,'  tre T             ',(raTrUp2(iFr),iFr=i1,i2)
PRINT *,'  tre R             ',(raReUp2(iFr),iFr=i1,i2)
PRINT *,'  tre E             ',(raEmissUp2(iFr),iFr=i1,i2)
PRINT *,'  tre S             ',(raSunUp2(iFr),iFr=i1,i2)

RETURN
END SUBROUTINE PrintStuffTwoStream

!************************************************************************
! this subroutine does stuff for one layer

SUBROUTINE OneLayerScatProp(iLay,iDir,iNumLaysInCld,i1or2,iAccumulate,
! input params  &
iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP,  &
    raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld,  &
    iPhase,raPhasePoints,raComputedPhase,
! input/output params  &
raCumSum1,raCumSum2,                         !AccumulateSolarDepth  &
    raTau1,raTau2,muSat,                       !optical depths,angle
! output params  &
raW0,raAsym,raTb,raBeta,                     !LayerScatteringProp  &
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, !asymcoeffs_solar  &
    raTrUp,raReUp,raEmissUp,raSunUp,         !t_r_e_arb_up_Solar_prof  &
    raTrDown,raReDown,raEmissDown,raSunDown, !t_r_e_arb_dn_Solar_prof  &
    raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm)
!t_r_e_streams_Solar_prof


INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iDir
NO TYPE, INTENT(IN OUT)                  :: iNumLaysIn
INTEGER, INTENT(IN)                      :: i1or2
NO TYPE, INTENT(IN OUT)                  :: iAccumulat
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars which determine how to use this subroutine


INTEGER :: iNumLaysInCld                !one or more cloud layers
INTEGER :: iAccumulate                  !do we cumulatively sum opt dep

! input params
INTEGER :: iNumLayer                    !number of layers in atm
INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
REAL :: raFreq(kMaxPts)                 !wavenumbers
REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
!these next three are self explanatory
REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
REAL :: raaAsym(kMaxPts,kMixFilRows)
REAL :: muSun                          !solar angle
INTEGER :: iPhase                       !use supplied phase fcn, or HG
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
REAL :: rFracTop,rFracBot
! input/output params
REAL :: raCumSum1(kMaxPts),raCumSum2(kMaxPts)
REAL :: raTau1(kMaxPts),raTau2(kMaxPts),muSat
! output parameters
REAL :: raW0(kMaxPts),raAsym(kMaxPts)
REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raTb(kMaxPts),raTt(kMaxPts)
REAL :: raBeta(kMaxPts),raSun0(kMaxPts)
REAL :: raA(kMaxPts),raB(kMaxPts),raE(kMaxPts),raF(kMaxPts)
REAL :: raG(kMaxPts),raH(kMaxPts)
REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts)
!!!! arbitrary angle stuff
REAL :: raTrUp(kMaxPts),raReUp(kMaxPts),raEmissUp(kMaxPts)
REAL :: raTrDown(kMaxPts),raReDown(kMaxPts),raEmissDown(kMaxPts)
REAL :: raSunUp(kMaxPts),raSunDown(kMaxPts)
!!! stream angle stuff
REAL :: raT(kMaxPts),raTstar(kMaxPts)
REAL :: raR(kMaxPts),raRstar(kMaxPts)
REAL :: raEp(kMaxPts),raEm(kMaxPts)
REAL :: raFp(kMaxPts),raFm(kMaxPts)
REAL :: raSp(kMaxPts),raSm(kMaxPts)

! local vars
INTEGER :: iL,iBeta,MP2Lay,iFr,iDebug
REAL :: rMPTemp,rFrac,ttorad
REAL :: rBeta,raBeta_RadOrTemp(kMaxPts),r1,r2
REAL :: raCumSum(kMaxPts),raDet(kMaxPts),raTau(kMaxPts)

iL      = iaRadLayer(iLay)
muSat = ABS(COS(raLayAngles(MP2Lay(iL))*kPi/180.0))
iBeta   = MOD(iL,kProfLayer)
IF (iBeta == 0) THEN
  iBeta = kProfLayer
END IF
rMPTemp = TEMP(iBeta)

rBeta   = LOG(TEMP(iBeta+1)/TEMP(iBeta))
DO iFr = 1,kMaxPts
  r1 = ttorad(raFreq(iFr),TEMP(iBeta+1))
  r2 = ttorad(raFreq(iFr),TEMP(iBeta))
  raBeta_RadOrTemp(iFr) = LOG(r1/r2)
END DO

IF (iAccumulate == 1) THEN
!total depth so far = 0.0
  CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1)
ELSE IF (iAccumulate == 2) THEN
!total depth = added on
  CALL AccumulateSolarDepth(raCumSum1,raTau2,+1)
ELSE IF (iAccumulate == 3) THEN
!total depth = added on
  CALL AccumulateSolarDepth(raCumSum1,raTau1,+1)
ELSE
  WRITE(kStdErr,*) 'iAccumulate = ',iAccumulate,' in OneLayerScatProp!!!'
  CALL DoStop
END IF

CALL FindTheFrac(iDir,iL,iaRadLayer,iNumLayer,rFracTop,rFracBot,rFrac)
CALL LayerScatteringProp( rFrac,raFreq,rMPTemp,iL,raaExt,raaScat,raaAsym,  &
    rBeta,raBeta_RadOrTemp, raW0,raAsym,raTau,raTb,raBeta)

IF (i1or2 == 1) THEN
  DO iFr = 1,kMaxPts
    raTau1(iFr)   = raTau(iFr)
    raCumSum(iFr) = raCumSum1(iFr)
  END DO
ELSE IF (i1or2 == 2) THEN
  DO iFr = 1,kMaxPts
    raTau2(iFr) = raTau(iFr)
    raCumSum(iFr) = raCumSum2(iFr)
  END DO
ELSE
  WRITE(kStdErr,*) 'i1or2 = ',i1or2,' in OneLayerScatProp!!!'
  CALL DoStop
END IF

!do the stuff for the stream angles, for up going radiation
CALL asymcoeffs_solar( raAsym,raW0,raTb,raRadBb,raRadBt,raTau,raBeta,  &
    radSolarCld,raCumSum,muSun, iPhase,raPhasePoints,raComputedPhase,  &
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)

IF (iNumLaysInCld > 1) THEN
  CALL t_r_e_streams_Solar_prof(  &
      raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
      raTau,muSat,muSun,raAsym,raW0,raCumSum, radSolarCld,raTb,  &
      raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm)
END IF

!do the stuff for arbirtary angles, for up going radiation
CALL t_r_e_arb_up_Solar_prof(  &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, radSolarCld,raTb,  &
    iPhase,raPhasePoints,raComputedPhase, raTrUp,raReUp,raEmissUp,raSunUp)

!do the stuff for arbirtary angles, for dn going radiation
CALL t_r_e_arb_down_Solar_prof(  &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,  &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, radSolarCld,raTb,  &
    iPhase,raPhasePoints,raComputedPhase, raTrDown,raReDown,raEmissDown,raSunDown)

iDebug = -1
IF (iDebug > 0) THEN
  CALL PrintStuffTwoStream(iLay,iL, raAsym,raW0,raTb,raRadBb,raRadBt,raTau,  &
      raBeta,raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  &
      raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm,  &
      raTrUp,raReUp,raEmissUp,raSunUp)
END IF

RETURN
END SUBROUTINE OneLayerScatProp

!************************************************************************
