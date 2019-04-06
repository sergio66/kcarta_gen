! Copyright 2001
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
    SUBROUTINE Cloud_UpOrDownLook(iNumLayer,iDir,iLocalCldTop,iLocalCldBot, &
    rFracTop,rFracBot, &
    iaRadLayer,raLayAngles,TEMP,raFreq, &
    raaExt,raaScat,raaAsym,radSolarCld,muSun,muSat, &
    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12, &
    raTrDown12,raReDown12,raEmissDown12,raSunDown12, &
    raW0,raAsym, &
    iPhase,raPhasePoints,raComputedPhase)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    INTEGER :: iDir                         !up or down looking instr
    INTEGER :: iNumLayer                    !number of layers in atm
    INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
    REAL :: raFreq(kMaxPts)                 !wavenumbers
    REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
! hese next three are self explanatory
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: muSun                          !solar angle
    INTEGER :: iPhase                  !use supplied phase fcn, or HG
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
    REAL :: rFracTop,rFracBot

! output parameters
    REAL :: muSat                         !view angle at lowest cloud layer
! hese next few are self explanatory : optical depth, and the
! umulative up/down transmission, reflection, emission, solar trans
    REAL :: raTau12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
    REAL :: raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)
! these are the lowest cloud layer asymmetry and single scattering
    REAL :: raW0(kMaxPts),raAsym(kMaxPts)
          
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
    REAL :: raTrUp1(kMaxPts),raReUp1(kMaxPts),raEmissUp1(kMaxPts)
    REAL :: raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
    REAL :: raTrDown1(kMaxPts),raReDown1(kMaxPts),raEmissDown1(kMaxPts)
    REAL :: raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
    REAL :: raSunUp1(kMaxPts),raSunUp2(kMaxPts)
    REAL :: raSunDown1(kMaxPts),raSunDown2(kMaxPts)

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
    ELSEIF (iDir < 0) THEN
        N = iLocalCldBot - iLocalCldTop + 1
    ELSEIF (iDir == 0) THEN
        write (kStdErr,*) 'iDir = ',iDir,' invalid number!!!'
        Call DoStop
    END IF

    DO iFr = 1,kMaxPts
        raCumSum1(iFr) = 0.0
        raCumSum2(iFr) = 0.0
    END DO

    IF (N < 1) THEN
        write(kStdErr,*) 'Huh? negative number of cld lays in Cld_UpDnLook'
        write(kStdErr,*) 'Local CldTop,CldBot = ',iLocalCldTop,iLocalCldBot
        CALL DoStop
    END IF

! -------------- cloud has only one layer ---------------------------------
    IF (N == 1) THEN             !bloody simple
        iLay    = iLocalCldTop

        CALL OneLayerScatProp(iLay,iDir,N,2,1, &
    ! input params
        iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP, &
        raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld, &
        iPhase,raPhasePoints,raComputedPhase, &
    ! output params
        raCumSum1,raCumSum2,                          & !AccumulateSolarDept
        raTau1,raTau2,muSat,                          & !optical depths,angle
        raW0,raAsym,raTb,raBeta,                      & !LayerScatteringProp
        raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  & !asymcoeffs_solar
        raTrUp2,raReUp2,raEmissUp2,raSunUp2,          & !t_r_e_arb_up_Solar_prof
        raTrDown2,raReDown2,raEmissDown2,raSunDown2,  & !t_r_e_arb_dn_Solar_prof
        raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
    ! f more than one layer, need t_r_e_streams_Solar_prof

    ! eed to assign raXXX12 --> raXXX2 for output purposes
        CALL update_streams_arb_solar( &
        raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12, &
        raEmissUp12,raEmissDown12,raSunUp12,raSunDown12, &
        raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2, &
        raEmissUp2, raEmissDown2, raSunup2, raSunDown2)

    ENDIF

! -------------- cloud has many layers ---------------------------------
    IF (N > 1) THEN             !have to add the layers together
    ! ecause the sun takes away up down symmetry, we have to do this
    ! rom the top down

    ! o the stuff for the arbitrary angles, for CLOUD LAYER 2
        iLay    = iLocalCldTop
        CALL OneLayerScatProp(iLay,iDir,N,2,1, &
    ! input params
        iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP, &
        raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld, &
        iPhase,raPhasePoints,raComputedPhase, &
    ! output params
        raCumSum1,raCumSum2,                          & !AccumulateSolarDepth
        raTau1,raTau2,muSat,                          & !optical depths,angle
        raW0,raAsym,raTb,raBeta,                      & !LayerScatteringProp
        raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  & !asymcoeffs_solar
        raTrUp2,raReUp2,raEmissUp2,raSunUp2,          & !t_r_e_arb_up_Solar_prof
        raTrDown2,raReDown2,raEmissDown2,raSunDown2,  & !t_r_e_arb_dn_Solar_prof
        raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
    ! f more than one layer, need t_r_e_streams_Solar_prof

    ! o the stuff for the arbitrary angles, for CLOUD LAYER 1
        iLay = iLay - iDir
        CALL OneLayerScatProp(iLay,iDir,N,1,2, &
    ! input params
        iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP, &
        raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld, &
        iPhase,raPhasePoints,raComputedPhase, &
    ! output params
        raCumSum1,raCumSum2,                          & !AccumulateSolarDepth
        raTau1,raTau2,muSat,                          & !optical depths,angle
        raW0,raAsym,raTb,raBeta,                      & !LayerScatteringProp
        raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  & !asymcoeffs_solar
        raTrUp1,raReUp1,raEmissUp1,raSunUp1,          & !t_r_e_arb_up_Solar_prof
        raTrDown1,raReDown1,raEmissDown1,raSunDown1,  & !t_r_e_arb_dn_Solar_prof
        raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
    ! f more than one layer, need t_r_e_streams_Solar_prof

    ! oop over the layers
        DO  iDp = 1,N-1
        !!! add layers together

            CALL addstar_arb_solar( &
            raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1, &
            raSunUp1,raSunDown1,raTau1, &
            raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2, &
            raSunUp2,raSunDown2,raTau2, &
            raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, &
            raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m, &
            raReUp12,raReDown12,raTrUp12,raTrDown12, &
            raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12, &
            raCumSum1,raCumSum2,muSat,muSun)

            CALL addstar_solar( &
            raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m, &
            raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, &
            raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m, &
            raCumSum1,raCumSum2,muSun)

            iLay = iLay - iDir
            IF (iDp < (N-1)) THEN          !!!!!add on new layer stuff
                CALL OneLayerScatProp(iLay,iDir,N,1,3, &
            ! input params
                iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP, &
                raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld, &
                iPhase,raPhasePoints,raComputedPhase, &
            ! output params
                raCumSum1,raCumSum2,                         & !AccumulateSolarDepth
                raTau1,raTau2,muSat,                         & !optical depths
                raW0,raAsym,raTb,raBeta,                     & !LayerScatteringProp
                raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  & !asym_solar
                raTrUp1,raReUp1,raEmissUp1,raSunUp1,          & !arb_up_Solar_prof
                raTrDown1,raReDown1,raEmissDown1,raSunDown1,  & !arb_dn_Solar_prof
                raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
            ! f more than one layer, need t_r_e_streams_Solar_prof

                CALL update_streams_solar( &
                raR2, raR2star, raT2, raT2star, raE2p, raE2m, raS2p, raS2m, &
                raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
                CALL update_streams_arb_solar( &
                raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2, &
                raEmissUp2, raEmissDown2, raSunUp2, raSunDown2, &
                raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12, &
                raEmissUp12,raEmissDown12,raSunup12,raSunDown12)

            END IF
        END DO
    ENDIF

    RETURN
    end SUBROUTINE Cloud_UpOrDownLook

!************************************************************************
! this finds out if we use complete HG function, or first order approx
    SUBROUTINE FindIHG(avg,iHG,iArbHG,raW0)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iHG,iArbHG
    REAL :: raW0(kMaxPts)
    REAL :: avg

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
    ! use full HG phase function
        iHG    = 1
        iArbHG = 1
    ELSE
    ! use first order HG phase function
        iHG    = -1
        iArbHG = -1
    END IF
!c ******* this was sept 2001 *********

!c ******* try this dec 2001 *********
! use full HG phase function
    iHG    = 1
    iArbHG = 1
!c ******* try this dec 2001 *********

    RETURN
    end SUBROUTINE FindIHG

!************************************************************************
! this checks to see if the emissions are positive or not
! flag the negatives as "bad" by keeping a running count of wot's bad
    SUBROUTINE check_raE_raF(raE,raF,iF,iEminus,iFminus,iaEMinus,iaFMinus, &
    rEmin_positive_for_chunk,rFmin_positive_for_chunk)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input params
    REAL :: raE(kMaxPts),raF(kMaxPts)            !!the up/down emission
    INTEGER :: iF
! output params
    REAL :: rEmin_positive_for_chunk,rFmin_positive_for_chunk
    INTEGER :: iEminus,iFminus,iaEMinus(kMaxPts),iaFMinus(kMaxPts)


!      raE(iF) = max(raE(iF),0.0)
!      raF(iF) = max(raF(iF),0.0)
    IF ((raE(iF) > 0) .AND. (raE(iF) < rEmin_positive_for_chunk))THEN
        rEmin_positive_for_chunk = raE(iF)
    END IF
    IF ((raF(iF) > 0) .AND. (raF(iF) < rFmin_positive_for_chunk))THEN
        rFmin_positive_for_chunk = raF(iF)
    END IF

    IF (raE(iF) < 0) THEN
        iEminus = iEMinus + 1
        iaEMinus(iEMinus) = iF
    END IF
    IF (raF(iF) < 0) THEN
        iFminus = iFMinus + 1
        iaFMinus(iFMinus) = iF
    END IF

    RETURN
    end SUBROUTINE check_raE_raF

!************************************************************************
! this subroutine computes the coefficients
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/asymcoeffs_solar.m
! here, view angle === stream angle
    SUBROUTINE asymcoeffs_solar(raAsym,raW0,raBo,raRadBb,raRadBt, &
    raTau,raBeta,radSolarCld,raCumSum,muSun, &
    iPhase,raPhasePoints,raComputedPhase, &
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! really we do not need raA or raB, so we do not really need RadBb or RadBt
! so we do not need to compute rho or theta

! input parameters
    REAL :: raAsym(kMaxPts),raW0(kMaxPts)    !asym coeff, single scatter alb
    REAL :: raBeta(kMaxPts)                  !1/k * log(T(n+1)/T(n))
    REAL :: raBo(kMaxPts)                    !bottom level temperature
    REAL :: raRadBb(kMaxPts),raRadBt(kmaxPts)!planck, top,bottom rad
    REAL :: raTau(kMaxPts)             !total extinction of layer
    REAL :: radSolarCld(kMaxPts),muSun      !sun intensity at cloud top, angle
    REAL :: raCumSum(kMaxPts)          !adjusts sun intensity as necessary
    INTEGER :: iPhase                  !use supplied phase fcn, or HG
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! output parameters
    REAL :: raDet(kMaxPts),raA(kMaxPts),raB(kMaxPts)  !!det is more important!
    REAL :: raE(kMaxPts),raF(kMaxPts)  !boundary condtions for soln
    REAL :: raG(kMaxPts),raH(kMaxPts)  !boundary condtions for soln
    REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff

! local variables
    REAL :: mu_pm                             !two stream angle
    INTEGER :: iF,iHG,iArbHG,iReset_raEF
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

    mu_pm = 1/sqrt(3.0)
    muSun = abs(muSun)

    CALL FindIHG(avg,iHG,iArbHG,raW0)

    IF ((avg >= avgmin) .AND. (kSolar >= 0)) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

            b = (1-raAsym(iF))/2.0
            raKp(iF) = 1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
            raKm(iF) =-1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
             
            alpha = (raW0(iF)*(1-b)-1)
            raAp(iF) = -(mu_pm/(raW0(iF)*b))*(raKp(iF) + alpha/mu_pm)
            raAm(iF) = -(mu_pm/(raW0(iF)*b))*(raKm(iF) + alpha/mu_pm)

        !          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
            det = &
            (raW0(iF)*b)*(raW0(iF)*b)-alpha*alpha+(betaa*mu_pm)*(betaa*mu_pm)

            raE(iF) = raBo(iF)*(1-raW0(iF))/det
            raE(iF) = raE(iF)*(betaa*mu_pm + alpha - raW0(iF)*b)

            raF(iF) = raBo(iF)*(1-raW0(iF))/det
            raF(iF) = raF(iF)*(-betaa*mu_pm + alpha - raW0(iF)*b)

            CALL check_raE_raF(raE,raF,iF,iEminus,iFminus,iaEMinus,iaFMinus, &
            rEmin_positive_for_chunk,rFmin_positive_for_chunk)

            rSun = radSolarCld(iF) * exp(-raCumSum(iF)/muSun)

            IF (iHG <= 0) THEN
            !! before this loop was if (iHG .GT. 0) THEN
            !! but gave lousy answers when we fix the minus sign (see below)
                angle = mu_pm/muSun
                det = angle**2 - alpha**2 + (raW0(iF)*b)**2
                p_plus  = hg2(-abs(muSun),+abs(mu_pm),raAsym(iF), &
                iPhase,raPhasePoints,raComputedPhase)
                p_minus = hg2(-abs(muSun),-abs(mu_pm),raAsym(iF), &
                iPhase,raPhasePoints,raComputedPhase)
            !!! minus sign is incorrect, but gives best answers if we use this!
            !!! raG(iF) = (alpha+angle)*p_plus - raW0(iF)*b*p_minus
                raG(iF) = (alpha+angle)*p_plus + raW0(iF)*b*p_minus
                raG(iF) = raW0(iF)*rSun/kForP/det * raG(iF)
                raH(iF) = raW0(iF)*b*(-p_plus) + (alpha-angle)*p_minus
                raH(iF) = raW0(iF)*rSun/kForP/det * raH(iF)
            ! 1 = raG(iF)
            ! 2 = raH(iF)

            ELSE
                angle = mu_pm/muSun
                det = angle**2 - alpha**2 + (raW0(iF)*b)**2
                psi = 3.0*raAsym(iF)*muSun*mu_pm
                z3 = -alpha + raW0(iF)*b
                z4 =  alpha + raW0(iF)*b
                raG(iF) =  z3 + psi*z4 + angle*(psi-1)
                raG(iF) = -raW0(iF)*rSun/kForP/det * raG(iF)
                raH(iF) =  z3 - psi*z4 + angle*(psi+1)
                raH(iF) = -raW0(iF)*rSun/kForP/det * raH(iF)

            !            print *,iF,z1,z2,raG(iF),raH(iF),p_plus,p_minus,1-psi,1+psi
            END IF

            theta = raRadBb(iF) - raE(iF) - raG(iF)*exp(-raTau(iF)/muSun)
            rho   = raRadBt(iF) - raF(iF)*exp(betaa*raTau(iF)) - raH(iF)
             
            raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) - &
            raAm(iF)*exp(raKp(iF)*raTau(iF))

        ! these last two are really not needed, but what the heck!!!!!
            raA(iF) = ( theta*exp(raKm(iF)*raTau(iF)) - rho*raAm(iF))
            raA(iF) = raA(iF)/raDet(iF)
            raB(iF) = (-theta*exp(raKp(iF)*raTau(iF)) + rho*raAp(iF))
            raB(iF) = raB(iF)/raDet(iF)
        END DO

    ELSEIF ((avg >= avgmin) .AND. (kSolar < 0)) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

            b = (1.0-raAsym(iF))/2.0
            raKp(iF) = 1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
            raKm(iF) =-1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
             
            alpha = raW0(iF)*(1.0-b)-1.0
            raAp(iF) = -(mu_pm/(raW0(iF)*b))*(raKp(iF) + alpha/mu_pm)
            raAm(iF) = -(mu_pm/(raW0(iF)*b))*(raKm(iF) + alpha/mu_pm)

        !          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
            det = &
            (raW0(iF)*b)*(raW0(iF)*b)-alpha*alpha+(betaa*mu_pm)*(betaa*mu_pm)
        !          print *,iF,det

            raE(iF) = raBo(iF)*(1-raW0(iF))/det
            raE(iF) = raE(iF)*(betaa*mu_pm + alpha - raW0(iF)*b)

            raF(iF) = raBo(iF)*(1-raW0(iF))/det
            raF(iF) = raF(iF)*(-betaa*mu_pm + alpha - raW0(iF)*b)

            CALL check_raE_raF(raE,raF,iF,iEminus,iFminus,iaEMinus,iaFMinus, &
            rEmin_positive_for_chunk,rFmin_positive_for_chunk)
                       
            IF (iDebug > 0) THEN
                IF ((iF >= i1) .AND. (iF <= i2)) THEN
                    z1 = raBo(iF)*(1-raW0(iF))/det
                    z2 = (+betaa*mu_pm + alpha - raW0(iF)*b)
                    z3 = (-betaa*mu_pm + alpha - raW0(iF)*b)
                !              print *,iF,det,z1,z2,z3,raE(iF),raF(iF)
                    print *,iF,raW0(iF),raAsym(iF),b,alpha,betaa,det,raE(iF)
                END IF
            END IF

            raG(iF) = 0.0
            raH(iF) = 0.0

            theta = raRadBb(iF) - raE(iF) - raG(iF)*exp(-raTau(iF)/muSun)
            rho   = raRadBt(iF) - raF(iF)*exp(betaa*raTau(iF)) - raH(iF)
             
            raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) - &
            raAm(iF)*exp(raKp(iF)*raTau(iF))

        !!!! these last two are really not needed, but what the heck!!!!!
            raA(iF) = ( theta*exp(raKm(iF)*raTau(iF)) - rho*raAm(iF))
            raA(iF) = raA(iF)/raDet(iF)
            raB(iF) = (-theta*exp(raKp(iF)*raTau(iF)) + rho*raAp(iF))
            raB(iF) = raB(iF)/raDet(iF)
        END DO

    ELSE
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

            b = (1-raAsym(iF))/2
            raKp(iF) = 1/mu_pm
            raKm(iF) =-1/mu_pm
             
            raAp(iF) = 0.0
            raAm(iF) = 1.0
             
            raE(iF) = raBo(iF)/(1+betaa*mu_pm)
            raF(iF) = raBo(iF)/(1-betaa*mu_pm)
             
            raG(iF) = 0.0
            raH(iF) = 0.0
             
            theta = raRadBb(iF) - raE(iF)
            rho   = raRadBt(iF) - raF(iF)*exp(raTau(iF)*betaa)
            raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) - &
            raAm(iF)*exp(raKp(iF)*raTau(iF))
        !!! these last two are really not needed, but what the heck!!!!!
            raA(iF) = rho*exp(raKm(iF)*raTau(iF))
            raB(iF) = theta

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
            write(kStdWarn,*) 'oops, layer has negative E emiss!! "fixing" this'
            write(kStdWarn,*) 'by setting raE(bad) = ',rEmin_positive_for_chunk
            write(kStdWarn,*) 'for ',iEMinus, ' points'
            DO iF = 1,iEminus
                raE(iaEMinus(iF)) = rEmin_positive_for_chunk
                raE(iaEMinus(iF)) = 0.0                  !!before Apr 07
                raE(iaEMinus(iF)) = raBo(iaEMinus(iF))   !!Apr 07
            END DO
        END IF
        IF (iFminus > 0) THEN
            write(kStdWarn,*) 'oops, layer has negative F emiss!! "fixing" this'
            write(kStdWarn,*) 'by setting raF(bad) = ',rFmin_positive_for_chunk
            write(kStdWarn,*) 'for ',iFMinus, ' points'
            DO iF = 1,iFminus
                raF(iaFMinus(iF)) = rFmin_positive_for_chunk
                raF(iaFMinus(iF)) = 0.0                  !!before Apr 07
                raF(iaEMinus(iF)) = raBo(iaEMinus(iF))   !!Apr 07
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE asymcoeffs_solar

!************************************************************************
! this subroutine computes upwards R,T,Emiss for arbitrary view angle
! same as /home/sergio/MATLAB/RADTrans/UP/t_r_e_arb_up_Solar_prof.m
    SUBROUTINE t_r_e_arb_up_Solar_prof( &
! input params
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
    raTau,muSat,muSun,raAsym,raW0,raCumSum,radSolarCld,raBo, &
    iPhase,raPhasePoints,raComputedPhase, &
! output params
    tr_up,re_up,emiss_up,sun_up)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

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
    INTEGER :: iF
    REAL :: mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
    REAL :: rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
    REAL :: eup,fup,gup,hup,sup,hg2,betaa,mu,avgmin

    INTEGER :: iArbHG,iHG

    avgmin = kAvgMin

    mu_pm  = 1/sqrt(3.0)
    mu     = abs(muSat)
    muSun  = abs(muSun)

    CALL FindIHG(avg,iHG,iArbHG,raW0)

    IF (avg >= avgmin) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)
            rSun = radSolarCld(iF)*exp(-raCumSum(iF)/muSun)
            det = raDet(iF)
                      
            b = (1-raAsym(iF))/2
            alpha = (raW0(iF)*(1-b)-1)
             
            alpha_p = raW0(iF)/2/mu * (1 + 3.0*raAsym(iF)*mu_pm*mu)
            alpha_m = raW0(iF)/2/mu * (1 - 3.0*raAsym(iF)*mu_pm*mu)

            vp = 1/mu + raKp(iF)
            vm = 1/mu + raKm(iF)
            vb = betaa + 1/mu
            vs = 1/muSun + 1/mu
             
        ! valuate at raTau(iF)
            evm = exp(vm*raTau(iF)) - 1
            evp = exp(vp*raTau(iF)) - 1
            evb = exp(vb*raTau(iF)) - 1
            evs = (exp(vs*raTau(iF))-1)/vs
             
            omega = evb/vb
            gamma_a = (alpha_p * raAp(iF) + alpha_m) * evp/vp
            gamma_b = (alpha_p * raAm(iF) + alpha_m) * evm/vm
             
            kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(iF)*b)**2
            kappa = (betaa*mu_pm + alpha - raW0(iF)*b)/kappa

            tr_up(iF) = (gamma_a*exp(raKm(iF)*raTau(iF)) - &
            gamma_b*exp(raKp(iF)*raTau(iF)))/det
            re_up(iF) = (raAp(iF) *gamma_b - raAm(iF)*gamma_a)/det

            eup = -tr_up(iF) + omega*(alpha_p + 1/(mu*kappa))
            fup = -re_up(iF) + omega*alpha_m
            emiss_up(iF) = raE(iF)*eup + raF(iF)*fup

            IF (kSolar >= 0) THEN
                gup = exp(-raTau(iF)/muSun)*(evs*alpha_p - tr_up(iF))
                hup = exp(-raTau(iF)/muSun)*evs*alpha_m - re_up(iF)
                IF (iarbhg > 0) THEN
                    sup = exp(-raTau(iF)/muSun)*raW0(iF)*evs/kForP/mu * &
                    hg2(-abs(muSun),abs(mu),raAsym(iF), &
                    iPhase,raPhasePoints,raComputedPhase)
                ELSE
                    sup = exp(-raTau(iF)/muSun)*raW0(iF)*evs/kForP/mu * &
                    (1-3.0*raAsym(iF)*mu*muSun)
                END IF
                sun_up(iF) = (raG(iF)*gup + raH(iF)*hup + rSun*sup)/rSun
            ! un            sun_up(iF) = sup
            ELSE
                sun_up(iF) = 0.0
            END IF

        END DO

    ELSEIF ((avg < avgmin) .AND. (kTemperVary < 0)) THEN
        DO iF = 1,kMaxPts
        !!!!!!!!!no scattering, clear sky solns
            betaa = raBeta(iF)

        ! evaluate at tau0
            tr_up(iF) = 0.0        !remember, no 2stream radiation knocked here
            re_up(iF) = 0.0        !remember, no 2stream radiation knocked here

            eup = raBo(iF)*(1 - exp(-raTau(iF)/mu))
            eup = eup*exp(raTau(iF)/mu)
            fup = 0.0
            gup = 0.0          !remember, no scattering of solar beam to here
            hup = 0.0          !remember, no scattering of solar beam to here

            emiss_up(iF) = eup
            sun_up(iF)   = 0.0
        END DO

    ELSEIF ((avg < avgmin) .AND. (kTemperVary > 0)) THEN
        DO iF = 1,kMaxPts
        !!!!!!!!!no scattering, clear sky solns
            betaa = raBeta(iF)

        ! evaluate at tau0
            tr_up(iF) = 0.0        !remember, no 2stream radiation knocked here
            re_up(iF) = 0.0        !remember, no 2stream radiation knocked here

            eup = raBo(iF)/(1+betaa*mu)* &
            (exp(betaa*raTau(iF))-exp(-raTau(iF)/mu))
            eup = eup*exp(raTau(iF)/mu)

            fup = 0.0
            gup = 0.0          !remember, no scattering of solar beam to here
            hup = 0.0          !remember, no scattering of solar beam to here

            emiss_up(iF) = eup
            sun_up(iF)   = 0.0
        END DO
    END IF

    RETURN
    end SUBROUTINE t_r_e_arb_up_Solar_prof

!************************************************************************
! this subroutine computes downwards R,T,Emiss for arbitrary view angle
! same as /home/sergio/MATLAB/RADTrans/DOWN/t_r_e_arb_down_Solar_prof.m
    SUBROUTINE t_r_e_arb_down_Solar_prof( &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, &
    radSolarCld,raBo, &
    iPhase,raPhasePoints,raComputedPhase, &
    tr_dn,re_dn,emiss_dn,sun_dn)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: muSat,muSun,raTau(kMaxPts)  !view,solar angle and extinction
    REAL :: raCumSum(kMaxPts),radSolarCld(kMaxPts)!direct solar extinction
    REAL :: raAsym(kMaxPts),raW0(kMaxPts)   !asymmetry and albedo
    REAL :: raDet(kMaxPts)                  !determinant
    REAL :: raE(kMaxPts),raF(kMaxPts)       !used for emission
    REAL :: raG(kMaxPts),raH(kMaxPts)       !used for sun
    REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff
    REAL :: raBeta(kMaxPts),raBo(kMaxPts)    !temp of lower level, variation
    INTEGER :: iPhase                  !use supplied phase fcn, or HG
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
          
! output parameters
    REAL :: tr_dn(kMaxPts),re_dn(kMaxPts),emiss_dn(kMaxPts),sun_dn(kMaxPts)

! local variables
    INTEGER :: iF
    REAL :: mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
    REAL :: rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
    REAL :: edn,fdn,gdn,hdn,sdn,hg2,betaa,mu,avgmin
    INTEGER :: iArbHG,iHG

    avgmin = kAvgMin

    mu_pm  = 1/sqrt(3.0)
    mu     = abs(muSat)
    muSun  = abs(muSun)

    CALL FindIHG(avg,iHG,iArbHG,raW0)

    IF (avg >= avgmin) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)
            rSun = radSolarCld(iF)*exp(-raCumSum(iF)/muSun)
            det = raDet(iF)

            b = (1-raAsym(iF))/2
            alpha = (raW0(iF)*(1-b)-1)
             
            alpha_p = raW0(iF)/2/mu * (1 - 3.0*raAsym(iF)*mu_pm*mu)
            alpha_m = raW0(iF)/2/mu * (1 + 3.0*raAsym(iF)*mu_pm*mu)

            vp = -1/mu + raKp(iF)
            vm = -1/mu + raKm(iF)
            vb = betaa - 1/mu
            vs = 1/muSun - 1/mu
             
        ! valuate at 0.0
            evm = exp(vm*raTau(iF)) - 1
            evp = exp(vp*raTau(iF)) - 1
            evb = exp(vb*raTau(iF)) - 1
            evs = (exp(vs*raTau(iF))-1)/vs
             
            omega = evb/vb
            gamma_a = (alpha_p * raAp(iF) + alpha_m) * evp/vp
            gamma_b = (alpha_p * raAm(iF) + alpha_m) * evm/vm
             
            kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(iF)*b)**2
            kappa = (betaa*mu_pm + alpha - raW0(iF)*b)/kappa

            re_dn(iF) = (gamma_a*exp(raKm(iF)*raTau(iF)) - &
            gamma_b*exp(raKp(iF)*raTau(iF)))/det
            tr_dn(iF) = (raAp(iF) *gamma_b - raAm(iF)*gamma_a)/det

            edn = -re_dn(iF) + omega*(alpha_p + 1/(mu*kappa))
            fdn = -tr_dn(iF)*exp(betaa*raTau(iF)) + omega*alpha_m
            emiss_dn(iF) = raE(iF)*edn + raF(iF)*fdn

            IF (kSolar >= 0) THEN
                gdn = exp(-raTau(iF)/muSun)*(evs*alpha_p - re_dn(iF))
                hdn = exp(-raTau(iF)/muSun)*evs*alpha_m - tr_dn(iF)
                   
                IF (iarbhg > 0) THEN
                    sdn = exp(-raTau(iF)/muSun)*raW0(iF)*evs/kForP/mu * &
                    hg2(-abs(muSun),-abs(mu),raAsym(iF), &
                    iPhase,raPhasePoints,raComputedPhase)
                ELSE
                    sdn = exp(-raTau(iF)/muSun)*raW0(iF)*evs/kForP/mu * &
                    (1+3.0*raAsym(iF)*mu*muSun)
                END IF
                sun_dn(iF) = (raG(iF)*gdn + raH(iF)*hdn + rSun*sdn)/rSun

            ELSE
                sun_dn(iF) = 0.0
            END IF

        END DO

    ELSEIF ((avg < avgmin) .AND. (kTemperVary < 0)) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

        !!!!!!!!!no scattering, clear sky solns
                  
        ! evaluate at tau0
            tr_dn(iF) = 0.0        !remember, no 2stream radiation knocked here
            re_dn(iF) = 0.0        !remember, no 2stream radiation knocked here

            edn = 0.0
            fdn = raBo(iF)*(1-exp(-raTau(iF)/mu))

            gdn = 0.0          !remember, no scattering of solar beam to here
            hdn = 0.0          !remember, no scattering of solar beam to here

            emiss_dn(iF) = fdn
            sun_dn(iF)   = 0.0
        END DO

    ELSEIF ((avg < avgmin) .AND. (kTemperVary > 0)) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

        !!!!!!!!!no scattering, clear sky solns
                  
        ! evaluate at tau0
            tr_dn(iF) = 0.0        !remember, no 2stream radiation knocked here
            re_dn(iF) = 0.0        !remember, no 2stream radiation knocked here

            edn = 0.0
            fdn = (1-exp(-raTau(iF)/mu))*exp(betaa*raTau(iF))
            fdn = raBo(iF)/(1 - betaa*mu) * fdn

            gdn = 0.0          !remember, no scattering of solar beam to here
            hdn = 0.0          !remember, no scattering of solar beam to here

            emiss_dn(iF) = fdn
            sun_dn(iF)   = 0.0
        END DO
    END IF

    DO iF = 1,kMaxPts
        tr_dn(iF)    = tr_dn(iF)    * exp(raTau(iF)/mu)
        re_dn(iF)    = re_dn(iF)    * exp(raTau(iF)/mu)
        emiss_dn(iF) = emiss_dn(iF) * exp(raTau(iF)/mu)
        sun_dn(iF)   = sun_dn(iF)   * exp(raTau(iF)/mu)
    END DO

    RETURN
    end SUBROUTINE t_r_e_arb_down_Solar_prof

!************************************************************************
! this is for the twostream radiation
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/t_r_e_streams_Solar_prof.m
    SUBROUTINE t_r_e_streams_Solar_prof( &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
    raTau1,muSat,muSun,raAsym,raW0,raCumSum1,radSolarCld,raTb, &
    raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output params
    REAL :: raT1(kMaxPts),raT1star(kMaxPts)    !transmissions
    REAL :: raR1(kMaxPts),raR1star(kMaxPts)    !reflections
    REAL :: raE1p(kMaxPts),raE1m(kMaxPts)      !emissions
    REAL :: raS1p(kMaxPts),raS1m(kMaxPts)      !sun stuff

! input params
    REAL :: muSun,muSat                     !solar,view angle (not needed)
    REAL :: raAsym(kMaxPts)                    !asymmetry
    REAL :: radSolarCld(kMaxPts),raCumSum1(kMaxPts)  !solar stuff
    REAL :: raTb(kMaxPts)                      !radiance at bottom
    REAL :: raW0(kMaxPts),raBeta(kMaxPts)      !albedo and temperature variation
    REAL :: raKp(kMaxPts),raKm(kMaxPts)        !eigenvalues
    REAL :: raAp(kMaxPts),raAm(kMaxPts)        !eigenvectors
    REAL :: raTau1(kMaxPts),raDet(kMaxPts)     !optical depth, determinant
    REAL :: raE(kMaxPts),raF(kMaxPts)          !up and down emission
    REAL :: raG(kMaxPts),raH(kMaxPts)          !up and down sun stuff
          

! local variables
    INTEGER :: iF,iHG,iArbHG
    REAL :: avg,det,bw1,bw2,hg2,rSun,tau,mu_pm,tau0
    REAL :: tr_up,re_up,eup,fup,gup,hup,sun_up,emiss_up
    REAL :: tr_down,re_down,edn,fdn,gdn,hdn,sun_down,emiss_down,betaa
    REAL :: avgmin

    avgmin = kAvgMin

    CALL FindIHG(avg,iHG,iArbHG,raW0)

    mu_pm   = 1/sqrt(3.0)
    muSun  = abs(muSun)

    IF (avg > avgmin) THEN
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)
            rSun = radSolarCld(iF)*exp(-raCumSum1(iF)/muSun)

            tau  = raTau1(iF)
            tau0 = raTau1(iF)
            tr_up = raAp(iF)*exp(raKm(iF)*tau0)*exp(raKp(iF)*tau) - &
            raAm(iF)*exp(raKp(iF)*tau0)*exp(raKm(iF)*tau)
            tr_up = tr_up/raDet(iF)
            re_up = raAp(iF)*raAm(iF)*(exp(raKm(iF)*tau) - exp(raKp(iF)*tau))
            re_up = re_up/raDet(iF)
            eup = exp(betaa*tau) - tr_up
            fup = -re_up*exp(betaa*tau0)
            gup = exp(-tau0/muSun)*(exp(tau/muSun) - tr_up)
            hup = -re_up
            emiss_up = raE(iF)*eup + raF(iF)*fup
            IF (kSolar >= 0) THEN
                sun_up = (raG(iF)*gup+ raH(iF)*hup)/rSun
            ELSE
                sun_up = 0.0
            END IF
            raT1(iF) = tr_up
            raR1(iF) = re_up
            raE1p(iF) = emiss_up
            raS1p(iF) = sun_up

            tau = 0.0
            tau0 = raTau1(iF)
            tr_down = raAp(iF)*exp(raKm(iF)*tau)-raAm(iF)*exp(raKp(iF)*tau)
            tr_down = tr_down/raDet(iF)
            re_down = exp(raKm(iF)*tau0)*exp(raKp(iF)*tau) - &
            exp(raKp(iF)*tau0)*exp(raKm(iF)*tau)
            re_down = re_down/raDet(iF)
            edn = -re_down
            fdn = -tr_down*exp(betaa*tau0) + exp(betaa*tau)
            emiss_down = raE(iF)*edn + raF(iF)*fdn
            IF (kSolar >= 0) THEN
                gdn = exp(-tau0/muSun)*(-re_down)
                hdn = exp(-(tau0-tau)/muSun) - tr_down
                sun_down = (raG(iF)*gdn + raH(iF)*hdn)/rSun
            ELSE
                sun_down = 0.0
            END IF

            raT1star(iF) = tr_down
            raR1star(iF) = re_down
            raE1m(iF) = emiss_down
            raS1m(iF) = sun_down

        END DO

    ELSE
        DO iF = 1,kMaxPts
            betaa = raBeta(iF)

            tau  = raTau1(iF)
            tau0 = raTau1(iF)

            tr_up = exp(-raKp(iF)*tau)
            re_up = 0.0
            eup = exp(betaa*tau) - tr_up
            fup = 0.0
            gup = 0.0
            hup = 0.0
            emiss_up = raE(iF)*eup + raF(iF)*fup
            sun_up = 0.0
            raT1(iF) = tr_up
            raR1(iF) = re_up
            raE1p(iF) = emiss_up
            raS1p(iF) = sun_up
             
            tau = 0.0
            tr_down = exp(-raKp(iF)*(tau0-tau))
            re_down = 0.0
            edn = 0.0
            fdn = 1 - exp(-tau0/mu_pm)*exp(betaa*tau0)
            gdn = 0.0
            hdn = 0.0
            emiss_down = raE(iF)*edn + raF(iF)*fdn
            sun_down   = 0.0
            raT1star(iF) = tr_down
            raR1star(iF) = re_down
            raE1m(iF) = emiss_down
            raS1m(iF) = sun_down

        END DO
    END IF

    RETURN
    end SUBROUTINE t_r_e_streams_Solar_prof

!************************************************************************
! this subroutine adds together stuff for the arbitrary angles
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/addstar_arb_solar.m
    SUBROUTINE  addstar_arb_solar( &
    raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1, &
    raSunUp1,raSunDown1,raTau1, &
    raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2, &
    raSunUp2,raSunDown2,raTau2, &
    raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, &
    raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m, &
    raReUp12,raReDown12,raTrUp12,raTrDown12, &
    raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12, &
    raCumSum1,raCumSum2,muSat,muSun)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
!! layer stuff for the arb angles, layer n
    REAL :: raReUp1(kMaxPts),raReDown1(kMaxPts)
    REAL :: raTrUp1(kMaxPts),raTrDown1(kMaxPts)
    REAL :: raEmissUp1(kMaxPts),raEmissDown1(kMaxPts),ratau1(kMaxPts)
    REAL :: raSunUp1(kMaxPts),raSunDown1(kMaxPts)
!! layer stuff for the arb angles, layer n+1
    REAL :: raReUp2(kMaxPts),raReDown2(kMaxPts)
    REAL :: raTrUp2(kMaxPts),raTrDown2(kMaxPts)
    REAL :: raEmissUp2(kMaxPts),raEmissDown2(kMaxPts),ratau2(kMaxPts)
    REAL :: raSunUp2(kMaxPts),raSunDown2(kMaxPts)
!! layer stuff for the stream angles, layer n
    REAL :: raR1(kMaxPts),raT1(kMaxPts),raR1star(kMaxPts)
    REAL :: raT1star(kMaxPts)
    REAL :: raE1p(kMaxPts),raE1m(kMaxPts)
    REAL :: raS1p(kMaxPts),raS1m(kMaxPts)
!! layer stuff for the stream angles, layer n+1
    REAL :: raR2(kMaxPts),raT2(kMaxPts),raR2star(kMaxPts)
    REAL :: raT2star(kMaxPts)
    REAL :: raE2p(kMaxPts),raE2m(kMaxPts)
    REAL :: raS2p(kMaxPts),raS2m(kMaxPts)
! arbitrary stuff
    REAL :: muSun,muSat
    REAL :: raCumSum1(kMaxPts),raCumSum2(kMaxPts)
! output params
!! layer stuff for the arb angles, layer n AND n+1
    REAL :: raReUp12(kMaxPts),raReDown12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raTrDown12(kMaxPts)
    REAL :: raEmissUp12(kMaxPts),raEmissDown12(kMaxPts),raTau12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)

! local variables
    INTEGER :: iF
    REAL :: mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat11,mat12,mat21,mat22
    REAL :: mu,det,a1,a2,tau12,e1,e2

    mu = abs(muSat)
    muSun = abs(muSun)
          
    DO iF = 1,kMaxPts

        det = 1-raR1(iF)*raR2star(iF)

    !!!!!!!reflection and transmission
        mat1 = 0.0
        mat2 = 0.0
        mat3 = raReUp2(iF)*exp(raTau1(iF)/mu)
        mat4 = 0.0
        mat5 = raT2(iF)*raReUp1(iF)/det
        mat6 = 0.0
        mat7 = raR1(iF)*raT2(iF)*raTrUp2(iF)/det*exp(raTau1(iF)/mu)
        mat11 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
         
        mat1 = raTrUp1(iF)
        mat2 = 0.0
        mat3 = 0.0
        mat4 = 0.0
        mat5 = raR2star(iF)*raT1star(iF)*raReUp1(iF)/det
        mat6 = 0.0
        mat7 = raT1star(iF)*raTrUp2(iF)/det*exp(raTau1(iF)/mu)
        mat12 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
         
        mat1 = raTrDown2(iF)
        mat2 = 0.0
        mat3 = 0.0
        mat4 = raR1(iF)*raT2(iF)*raReDown2(iF)/det
        mat5 = 0.0
        mat6 = raT2(iF)*raTrDown1(iF)/det*exp(raTau2(iF)/mu)
        mat7 = 0.0
        mat21 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
         
        mat1 = 0.0
        mat2 = raReDown1(iF)*exp(raTau2(iF)/mu)
        mat3 = 0.0
        mat4 = raT1star(iF)*raReDown2(iF)/det
        mat5 = 0.0
        mat6 = raR2star(iF)*raT1star(iF)*raTrDown1(iF)/det*exp(raTau2(iF)/mu)
        mat7 = 0.0
        mat22 = mat1 +  mat2 +  mat3 +  mat4 +  mat5 +  mat6 + mat7
         
    !!!!emission
        a1 = raEmissUp1(iF) + 0.0 + exp(raTau1(iF)/mu)*raEmissUp2(iF) + &
        &         0.0 + &
        raTrUp1(iF)/det*(raR2star(iF)*raE1p(iF)+raE2m(iF)) + &
        &         0.0 + &
        raTrUp2(iF)/det*(raR1(iF)*raE2m(iF)+raE1p(iF))*exp(raTau1(iF)/mu)
         
        a2 = raEmissDown2(iF) + exp(raTau2(iF)/mu)*raEmissDown1(iF) + 0.0 + &
        raReDown2(iF)/det*(raR1(iF)*raE2m(iF) + raE1p(iF)) + &
        &          0.0 + &
        raTrDown1(iF)/det*(raR2star(iF)*raE1p(iF)+raE2m(iF))*exp(raTau2(iF)/mu) &
        + 0.0
              
        raReUp12(iF)      = mat11
        raTrUp12(iF)      = mat12
        raTrDown12(iF)    = mat21
        raReDown12(iF)    = mat22
        raEmissUp12(iF)   = a1
        raEmissDown12(iF) = a2
                 
        raTau12(iF)  = raTau1(iF) + raTau2(iF)     !total optical depth so far

    !        print *,raEmissUp1(iF),raEmissUp2(iF),raTrUp1(iF),
    !     $         raTrUp2(iF),raE1p(iF),raE2m(iF),raR1(iF),raR2star(iF),det

    END DO

    IF (kSolar >= 0) THEN           !!!do the solar part
        DO iF = 1,kMaxPts
            det = 1 - raR1(iF)*raR2star(iF)
            e1 = exp(-raCumsum1(iF)/muSun)
            e2 = exp(-raCumsum2(iF)/muSun)

            raSunUp12(iF) = raSunUp1(iF)*e1 + &
            raSunUp2(iF)*e2*exp(raTau1(iF)/mu) + &
            &                     1/det*(raReUp1(iF)*(e1*raR2star(iF)*raS1p(iF) + &
            e2*raS2m(iF))) + &
            exp(raTau1(iF)/mu)*raTrUp2(iF)/det * &
            (e2*raR1(iF)*raS2m(iF) + e1*raS1p(iF))
            raSunDown12(iF) = raSunDown2(iF)*e2 + &
            raSunDown1(iF)*e1*exp(raTau2(iF)/mu) + &
            &                     1/det*(raReDown2(iF)*(e2*raR1(iF)*raS2m(iF) + &
            e1*raS1p(iF))) + &
            exp(raTau2(iF)/mu)*raTrDown1(iF)/det * &
            (e1*raR2star(iF)*raS1p(iF) + e2*raS2m(iF))
        END DO
    END IF

    RETURN
    END SUBROUTINE 

!************************************************************************
! this subroutine adds together the layers for the stream angles
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/addstar_solar.m
    SUBROUTINE addstar_solar( &
    raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m, &
    raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, &
    raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m, &
    raCumSum1,raCumSum2,muSun)
         
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
!! layer stuff for the stream angles, layer n
    REAL :: raR1(kMaxPts),raT1(kMaxPts),raR1star(kMaxPts)
    REAL :: raT1star(kMaxPts),raE1p(kMaxPts),raE1m(kMaxPts)
    REAL :: raS1p(kMaxPts),raS1m(kMaxPts)
!! layer stuff for the stream angles, layer n+1
    REAL :: raR2(kMaxPts),raT2(kMaxPts),raR2star(kMaxPts)
    REAL :: raT2star(kMaxPts),raE2p(kMaxPts),raE2m(kMaxPts)
    REAL :: raS2p(kMaxPts),raS2m(kMaxPts)
! for the sun
    REAL :: raCumSum1(kMaxPts),raCumSum2(kMaxPts),muSun
! output params
!! layer stuff for the stream angles, layer n AND n+1
    REAL :: raR12(kMaxPts),raT12(kMaxPts),raR12star(kMaxPts)
    REAL :: raT12star(kMaxPts),raE12p(kMaxPts),raE12m(kMaxPts)
    REAL :: raS12p(kMaxPts),raS12m(kMaxPts)
! local variables
    INTEGER :: iF
    REAL :: det,e1,e2

    muSun = abs(muSun)

    DO iF = 1,kMaxPts
        det = 1 - raR1(iF)*raR2star(iF)

        raR12(iF) = raT2star(iF)*raR1(iF)/det
        raT12star(iF) = raT2star(iF)/det
        raT12(iF) = raT1(iF)/det
        raR12star(iF) = raT1(iF)*raR2star(iF)/det

        raE12p(iF) = raE2p(iF) + (raR12(iF)*raE2m(iF)+raT12star(iF)*raE1p(iF))
        raE12m(iF) = raE1m(iF) + (raT12(iF)*raE2m(iF)+raR12star(iF)*raE1p(iF))
         
        raR12(iF)     = raR12(iF)*raT2(iF) + raR2(iF)
        raT12star(iF) = raT1star(iF)*raT12star(iF)
        raT12(iF)     = raT2(iF)*raT12(iF)
        raR12star(iF) = raR1star(iF) + raT1star(iF)*raR12star(iF)
    END DO

    IF (kSolar >= 0) THEN  !!!!do the solar part
        DO iF = 1,kMaxPts
            det = 1 - raR1(iF)*raR2star(iF)
            e1 = exp(-raCumSum1(iF)/muSun)
            e2 = exp(-raCumSum2(iF)/muSun)

            raS12p(iF) = (raT2star(iF)/det*raR1(iF)*raS2m(iF) + raS2p(iF))*e2 + &
            raT2star(iF)/det*raS1p(iF)*e1
            raS12m(iF) = (raT1(iF)/det*raR2star(iF)*raS1p(iF) + raS1m(iF))*e1 + &
            raT1(iF)/det*raS2m(iF)*e2
        END DO
    END IF

    RETURN
    end SUBROUTINE addstar_solar

!************************************************************************
! this subroutine updates the stream angle added layers
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/update_streams_solar.m
    SUBROUTINE update_streams_solar( &
    raR2,raR2star,raT2,raT2star,raE2p,raE2m,raS2p,raS2m, &
    raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! stream angle stuff
    REAL :: raT2(kMaxPts),raT2star(kMaxPts)
    REAL :: raR2(kMaxPts),raR2star(kMaxPts)
    REAL :: raE2p(kMaxPts),raE2m(kMaxPts)
    REAL :: raS2p(kMaxPts),raS2m(kMaxPts)
    REAL :: raT12(kMaxPts),raT12star(kMaxPts)
    REAL :: raR12(kMaxPts),raR12star(kMaxPts)
    REAL :: raE12p(kMaxPts),raE12m(kMaxPts)
    REAL :: raS12p(kMaxPts),raS12m(kMaxPts)

! local variables
    INTEGER :: iF
     
!!!only need to update the upper layer
    DO iF = 1,kMaxPts
        raR2(iF)     = raR12(iF)
        raR2star(iF) = raR12star(iF)
        raT2(iF)     = raT12(iF)
        raT2star(iF) = raT12star(iF)
        raE2p(iF)    = raE12p(iF)
        raE2m(iF)    = raE12m(iF)
        raS2p(iF)    = raS12p(iF)
        raS2m(iF)    = raS12m(iF)
    END DO

    RETURN
    end SUBROUTINE update_streams_solar

!************************************************************************
! this subroutine updates the arb angle added layers
! this just sets varNew == varOld
! same as /home/sergio/MATLAB/RADTrans/GENERAL_CLOUD/update_streams_arb_solar.m
    SUBROUTINE update_streams_arb_solar( &
    raTauNew,raReUpNew,raReDownNew,raTrUpNew,raTrDownNew, &
    raEmissUpNew,raEmissDownNew,raSunUpNew,raSunDownNew, &
    raTauOld,raReUpOld,raReDownOld,raTrUpOld,raTrDownOld, &
    raEmissUpOld,raEmissDownOld,raSunupOld,raSunDownOld)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! arbitrary angle stuff
    REAL :: raTauNew(kMaxPts),raTauOld(kMaxPts)
    REAL :: raTrUpNew(kMaxPts),raReUpNew(kMaxPts),raEmissUpNew(kMaxPts)
    REAL :: raTrUpOld(kMaxPts),raReUpOld(kMaxPts),raEmissUpOld(kMaxPts)
    REAL :: raTrDownNew(kMaxPts),raReDownNew(kMaxPts),raEmissDownNew(kMaxPts)
    REAL :: raTrDownOld(kMaxPts),raReDownOld(kMaxPts),raEmissDownOld(kMaxPts)
    REAL :: raSunUpNew(kMaxPts),raSunDownNew(kMaxPts)
    REAL :: raSunUpOld(kMaxPts),raSunDownOld(kMaxPts)

! local variables
    INTEGER :: iF

!!!only need to update the upper layer
    DO iF = 1,kMaxPts
        raTauNew(iF)       = raTauOld(iF)
        raReUpNew(iF)      = raReUpOld(iF)
        raReDownNew(iF)    = raReDownOld(iF)
        raTrUpNew(iF)      = raTrUpOld(iF)
        raTrDownNew(iF)    = raTrDownOld(iF)
        raEmissUpNew(iF)   = raEmissUpOld(iF)
        raEmissDownNew(iF) = raEmissDownOld(iF)
        raSunUpNew(iF)     = raSunUpOld(iF)
        raSunDownNew(iF)   = raSunDownOld(iF)
    END DO

    RETURN
    end SUBROUTINE update_streams_arb_solar

!************************************************************************
! this subroutine sets the scattering properties of the current layer
    SUBROUTINE LayerScatteringProp(rFrac, &
    raFreq,rMPTemp,iL,raaExt,raaScat,raaAsym, &
    rBeta,raBeta_RadOrTemp, &
! output params
    raW0,raAsym,raTau,raBo,raBeta)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input parameters
    REAL :: rFrac                 !full or fractional layer
    INTEGER :: iL                 !current MP number
    REAL :: rMPTemp               !current layer temperature
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts)
    REAL :: rBeta                 !log(T(n+1)/T(n)) -->
! og(B(iF,T(n+1))/B(iF,T(n)))
    REAL :: raBeta_RadOrTemp(kMaxPts)  !log(T(n+1)/T(n)) -->
! og(B(iF,T(n+1))/B(iF,T(n)))
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

        ! cc          !!!used only by kTwoStream, apparently!
        ! cc          !!!log(T(i)/T(i+1))
        ! cc          raBeta(iFr)     = 1/raTau(iFr) * rBeta

        !!!used by CHARTS, DISORT and RTSPEC
        !!!log(B(T(i),vj)/B(T(i+1),vj))
            raBeta(iFr)     = 1/raTau(iFr) * raBeta_RadOrTemp(iFr)

        END DO
    ELSE
    !!!! impossible
        write(kStdErr,*) 'In subroutine LayerScatteringProp :'
        write(kStdErr,*) 'Code claims fraction for this layer is : ',rFrac
        CALL DoStop
    END IF

    RETURN
    end SUBROUTINE LayerScatteringProp

!************************************************************************
! this accumulates the optical depth from top of cloud to top of layer
    SUBROUTINE AccumulateSolarDepth(raCumSum,raX,iDo)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iDo                 !1 to add on raX to raCumSum, -1 if not
    REAL :: raCumSum(kMaxPts),raX(kMaxPts)  !total optical depth from cloud top

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
    end SUBROUTINE AccumulateSolarDepth

!************************************************************************
! this function computes the phase function, assuming cos(phi1-phi2) factor = 1
    REAL Function hg2(mu1,mu2,g,iPhase,raPhasePoints,raComputedPhase)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL :: mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry
    INTEGER :: iPhase       !use supplied phase fcn, or HG
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    REAL :: normB,mu0,yexact,yexactNew,hg2_real
    DOUBLE PRECISION :: hg2_double
    INTEGER :: iJ

    IF (iPhase < 0) THEN
        hg2 = hg2_double(dble(mu1),dble(mu2),dble(g))
    END IF

    IF (iPhase > 0) THEN
    ! ********************* warning .. use rspl if the variables are real
    !!!compute mu0 = cos ofangle between the two
        mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)
        CALL spl(raPhasePoints,raComputedPhase,MaxPhase,mu0,yexactNew,1)
    ! ********************* warning .. use rspl if the variables are reals
        hg2 = yexactNew
    END IF

    RETURN
    end Function hg2

!************************************************************************
! this subroutine finds out the partial fraction for layer : default = +1.0
! for use in  LayerScatteringProp
    SUBROUTINE FindTheFrac(iUpDown,iL,iaRadLayer,iNumLayer, &
    rFracTop,rFracBot,rFrac)

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iUpDown                !!up or down look instr
    INTEGER :: iL                     !!layer we are looking at
    INTEGER :: iaRadLayer(kProfLayer) !!layering of current atm
    INTEGER :: iNumLayer              !!number of layers in atm
    REAL :: rFracTop,rFracBot
! output vars
    REAL :: rFrac

    rFrac = 1.0  !!default multiplier
    IF (iUpDown > 0) THEN    !!!!down look instr : radiation going up
        IF (iL == iaRadLayer(1)) THEN
            rFrac = rFracBot
        ELSEIF (iL == iaRadLayer(iNumLayer)) THEN
            rFrac = rFracTop
        END IF
    ELSEIF (iUpDown < 0) THEN    !!!!uplook instr : radiation going down
        IF (iL == iaRadLayer(1)) THEN
            rFrac = rFracTop
        ELSEIF (iL == iaRadLayer(iNumLayer)) THEN
            rFrac = rFracBot
        END IF
    END IF
     
    RETURN
    end SUBROUTINE FindTheFrac

!************************************************************************
! just a printer
    SUBROUTINE PrintStuffTwoStream(iLay,iL, &
    raAsym,raW0,raTb,raRadBb,raRadBt,raTau2, &
    raBeta,raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, &
    raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm, &
    raTrUp2,raReUp2,raEmissUp2,raSunUp2)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
    REAL :: raFreq(kMaxPts)                 !wavenumbers
    REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud

    REAL :: raW0(kMaxPts),raAsym(kMaxPts)
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raTb(kMaxPts),raTt(kMaxPts)
    REAL :: raBeta(kMaxPts),raSun0(kMaxPts)
    REAL :: raA(kMaxPts),raB(kMaxPts),raE(kMaxPts),raF(kMaxPts)
    REAL :: raG(kMaxPts),raH(kMaxPts)
    REAL :: raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts)

    REAL :: raT(kMaxPts),raTstar(kMaxPts)
    REAL :: raR(kMaxPts),raRstar(kMaxPts)
    REAL :: raEp(kMaxPts),raEm(kMaxPts)
    REAL :: raFp(kMaxPts),raFm(kMaxPts)
    REAL :: raSp(kMaxPts),raSm(kMaxPts)
          
! arbitrary angle stuff
    REAL :: raTau2(kMaxPts),raSunUp2(kMaxPts),raSunDown2(kMaxPts)
    REAL :: raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
    REAL :: raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
    REAL :: raDet(kMaxPts)

    INTEGER :: iLay,iL

    INTEGER :: iFr,i1,i2

    i1 = 1
    i2 = 2

    i1 = 9223
    i2 = 9223

    i1 = 9222
    i2 = 9224

    print *,'          '
    print *,' ---->   iLay,iL = ',iLay,iL
    print *,'from Asymcoeffs_solar ....'
    print *,'  asym solar Asym   ',(raAsym(iFr),iFr=i1,i2)
    print *,'  asym solar W0     ',(raW0(iFr),iFr=i1,i2)
    print *,'  asym solar Tb     ',(raTb(iFr),iFr=i1,i2)
    print *,'  asym solar Radbot ',(raRadBb(iFr),iFr=i1,i2)
    print *,'  asym solar Radtop ',(raRadBt(iFr),iFr=i1,i2)
    print *,'  asym solar Tau2   ',(raTau2(iFr),iFr=i1,i2)
    print *,'  asym solar Beta   ',(raBeta(iFr),iFr=i1,i2)
    print *,'  asym solar raDET  ',(raDet(iFr),iFr=i1,i2)
    print *,'  asym solar raA    ',(raA(iFr),iFr=i1,i2)
    print *,'  asym solar raB    ',(raB(iFr),iFr=i1,i2)
    print *,'  asym solar raE    ',(raE(iFr),iFr=i1,i2)
    print *,'  asym solar raF    ',(raF(iFr),iFr=i1,i2)
    print *,'  asym solar raG    ',(raG(iFr),iFr=i1,i2)
    print *,'  asym solar raH    ',(raH(iFr),iFr=i1,i2)
    print *,'  asym solar raKp   ',(raKp(iFr),iFr=i1,i2)
    print *,'  asym solar raKm   ',(raKm(iFr),iFr=i1,i2)
    print *,'  asym solar raAp   ',(raAp(iFr),iFr=i1,i2)
    print *,'  asym solar raAm   ',(raAm(iFr),iFr=i1,i2)
    print *,'from t_r_e_streams_Solar_prof'
    print *,'  asym solar raT    ',(raT(iFr),iFr=i1,i2)
    print *,'  asym solar raTstar',(raTstar(iFr),iFr=i1,i2)
    print *,'  asym solar raR    ',(raR(iFr),iFr=i1,i2)
    print *,'  asym solar raRstar',(raRstar(iFr),iFr=i1,i2)
    print *,'  asym solar raEp   ',(raEp(iFr),iFr=i1,i2)
    print *,'  asym solar raEm   ',(raEm(iFr),iFr=i1,i2)
    print *,'  asym solar raSp   ',(raSp(iFr),iFr=i1,i2)
    print *,'  asym solar raSm   ',(raSm(iFr),iFr=i1,i2)
    print *,'from  t_r_e_arb_up_Solar_prof'
    print *,'  tre T             ',(raTrUp2(iFr),iFr=i1,i2)
    print *,'  tre R             ',(raReUp2(iFr),iFr=i1,i2)
    print *,'  tre E             ',(raEmissUp2(iFr),iFr=i1,i2)
    print *,'  tre S             ',(raSunUp2(iFr),iFr=i1,i2)

    RETURN
    end SUBROUTINE PrintStuffTwoStream

!************************************************************************
! this subroutine does stuff for one layer
    SUBROUTINE OneLayerScatProp(iLay,iDir,iNumLaysInCld,i1or2,iAccumulate, &
! input params
    iNumLayer,iaRadlayer,rFracTop,rFracBot,raLayAngles,TEMP, &
    raFreq,raaExt,raaScat,raaAsym,muSun,radSolarCld, &
    iPhase,raPhasePoints,raComputedPhase, &
! input/output params
    raCumSum1,raCumSum2,                          & !AccumulateSolarDepth
    raTau1,raTau2,muSat,                        & !optical depths,angle
! output params
    raW0,raAsym,raTb,raBeta,                      & !LayerScatteringProp
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,  & !asymcoeffs_solar
    raTrUp,raReUp,raEmissUp,raSunUp,          & !t_r_e_arb_up_Solar_prof
    raTrDown,raReDown,raEmissDown,raSunDown,  & !t_r_e_arb_dn_Solar_prof
    raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm)
! _r_e_streams_Solar_prof
     
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input vars which determine how to use this subroutine
    INTEGER :: iLay                         !which mixed path layer
    INTEGER :: iDir                         !looking at top frac or bot frac
    INTEGER :: iNumLaysInCld                !one or more cloud layers
    INTEGER :: iAccumulate                  !do we cumulatively sum opt dep
    INTEGER :: i1or2                        !looking at lay 1 or lay 2
! input params
    INTEGER :: iNumLayer                    !number of layers in atm
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
    REAL :: raFreq(kMaxPts)                 !wavenumbers
    REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
! hese next three are self explanatory
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
    REAL ::    rMPTemp,rFrac,ttorad
    REAL ::    rBeta,raBeta_RadOrTemp(kMaxPts),r1,r2
    REAL ::    raCumSum(kMaxPts),raDet(kMaxPts),raTau(kMaxPts)

    iL      = iaRadLayer(iLay)
    muSat = abs(cos(raLayAngles(MP2Lay(iL))*kPi/180.0))
    iBeta   = MOD(iL,kProfLayer)
    IF (iBeta == 0) THEN
        iBeta = kProfLayer
    END IF
    rMPTemp = TEMP(iBeta)

    rBeta   = log(TEMP(iBeta+1)/TEMP(iBeta))
    DO iFr = 1,kMaxPts
        r1 = ttorad(raFreq(iFr),TEMP(iBeta+1))
        r2 = ttorad(raFreq(iFr),TEMP(iBeta))
        raBeta_RadOrTemp(iFr) = log(r1/r2)
    END DO

    IF (iAccumulate == 1) THEN
    ! otal depth so far = 0.0
        CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1)
    ELSEIF (iAccumulate == 2) THEN
    ! otal depth = added on
        CALL AccumulateSolarDepth(raCumSum1,raTau2,+1)
    ELSEIF (iAccumulate == 3) THEN
    ! otal depth = added on
        CALL AccumulateSolarDepth(raCumSum1,raTau1,+1)
    ELSE
        write(kStdErr,*) 'iAccumulate = ',iAccumulate,' in OneLayerScatProp!!!'
        CALL DoStop
    END IF

    CALL FindTheFrac(iDir,iL,iaRadLayer,iNumLayer,rFracTop,rFracBot,rFrac)
    CALL LayerScatteringProp( &
    rFrac,raFreq,rMPTemp,iL,raaExt,raaScat,raaAsym, &
    rBeta,raBeta_RadOrTemp, &
    raW0,raAsym,raTau,raTb,raBeta)

    IF (i1or2 == 1) THEN
        DO iFr = 1,kMaxPts
            raTau1(iFr)   = raTau(iFr)
            raCumSum(iFr) = raCumSum1(iFr)
        END DO
    ELSEIF (i1or2 == 2) THEN
        DO iFr = 1,kMaxPts
            raTau2(iFr) = raTau(iFr)
            raCumSum(iFr) = raCumSum2(iFr)
        END DO
    ELSE
        write(kStdErr,*) 'i1or2 = ',i1or2,' in OneLayerScatProp!!!'
        CALL DoStop
    END IF

! o the stuff for the stream angles, for up going radiation
    CALL asymcoeffs_solar( &
    raAsym,raW0,raTb,raRadBb,raRadBt,raTau,raBeta, &
    radSolarCld,raCumSum,muSun, &
    iPhase,raPhasePoints,raComputedPhase, &
    raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)

    IF (iNumLaysInCld > 1) THEN
        CALL t_r_e_streams_Solar_prof( &
        raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
        raTau,muSat,muSun,raAsym,raW0,raCumSum, &
        radSolarCld,raTb, &
        raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm)
    END IF

! o the stuff for arbirtary angles, for up going radiation
    CALL t_r_e_arb_up_Solar_prof( &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, &
    radSolarCld,raTb, &
    iPhase,raPhasePoints,raComputedPhase, &
    raTrUp,raReUp,raEmissUp,raSunUp)

! o the stuff for arbirtary angles, for dn going radiation
    CALL t_r_e_arb_down_Solar_prof( &
    raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta, &
    raTau,muSat,muSun,raAsym,raW0,raCumSum, &
    radSolarCld,raTb, &
    iPhase,raPhasePoints,raComputedPhase, &
    raTrDown,raReDown,raEmissDown,raSunDown)

    iDebug = -1
    IF (iDebug > 0) THEN
        CALL PrintStuffTwoStream(iLay,iL, &
        raAsym,raW0,raTb,raRadBb,raRadBt,raTau, &
        raBeta,raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm, &
        raT,raTstar,raR,raRstar,raEp,raEm,raSp,raSm, &
        raTrUp,raReUp,raEmissUp,raSunUp)
    END IF

    RETURN
    end SUBROUTINE OneLayerScatProp

!************************************************************************
