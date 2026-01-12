c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************************************************************************
c************************************************************************

c this does the scattering radiative transfer for a downlook instrument
c the only difference between Cloud_DownLook and Cloud_UpLook are that the
c ordering of the layers in iaRadLayer is inverted, and so 
c iLocalCldTop and iLocalCldBot are inverted

      SUBROUTINE Cloud_DownLook(iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,raWaves,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c input parameters
      INTEGER iNumLayer                    !number of layers in atm
      INTEGER iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
      INTEGER iaRadLayer(kProfLayer)       !atmosphere layering
      REAL raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
      REAL TEMP(MAXNZ)                     !temperature profile (levels)
      REAL raWaves(kMaxPts)                !wavenumbers
      REAL radSolarCld(kMaxPts)            !solar intensity at top of cloud
      !these next three are self explanatory
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL mu_sun                          !solar angle
c output parameters
      REAL mu_view                         !view angle at lowest cloud layer
      !these next few are self explanatory : optical depth, and the
      !cumulative up/down transmission, reflection, emission, solar trans
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)
      ! these are the lowest cloud layer asymmetry and single scattering
      REAL raW0(kMaxPts),raAsym(kMaxPts)
      
c local variables 
      INTEGER N,iFr,iLay,iL,iBeta,iDp,MP2LAY
      REAL rCos,rMPTemp

c general coeffs for the layers
      REAL rTT,rSunAngle,rBeta,rZeta,raBoo(kMaxPts)
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raTb(kMaxPts),raTt(kMaxPts)
      REAL raBeta(kMaxPts),raSun0(kMaxPts)
      REAL raA(kMaxPts),raB(kMaxPts),raE(kMaxPts),raF(kMaxPts)
      REAL raG(kMaxPts),raH(kMaxPts)
      REAL raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts)
      
c arbitrary angle stuff
      REAL raTau1(kMaxPts),raTau2(kMaxPts)
      REAL raTrUp1(kMaxPts),raReUp1(kMaxPts),raEmissUp1(kMaxPts)
      REAL raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
      REAL raTrDown1(kMaxPts),raReDown1(kMaxPts),raEmissDown1(kMaxPts)
      REAL raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
      REAL raSunUp1(kMaxPts),raSunUp2(kMaxPts)
      REAL raSunDown1(kMaxPts),raSunDown2(kMaxPts)

c stream angle stuff
      REAL raT1(kMaxPts),raT1star(kMaxPts),raT2(kMaxPts),raT2star(kMaxPts)
      REAL raR1(kMaxPts),raR1star(kMaxPts),raR2(kMaxPts),raR2star(kMaxPts)
      REAL raE1p(kMaxPts),raE1m(kMaxPts),raE2p(kMaxPts),raE2m(kMaxPts)
      REAL raF1p(kMaxPts),raF1m(kMaxPts),raF2p(kMaxPts),raF2m(kMaxPts)
      REAL raT12(kMaxPts),raT12star(kMaxPts)
      REAL raR12(kMaxPts),raR12star(kMaxPts)
      REAL raE12p(kMaxPts),raE12m(kMaxPts)
      REAL raS1p(kMaxPts),raS2p(kMaxPts),raS12p(kMaxPts)
      REAL raS1m(kMaxPts),raS2m(kMaxPts),raS12m(kMaxPts)
      REAL raCumSum1(kMaxPts),raCumSum2(kMaxPts),raDet(kMaxPts)

      N = iLocalCldTop - iLocalCldBot + 1

      DO iFr = 1,kMaxPts
        raCumSum1(iFr) = 0.0
        raCumSum2(iFr) = 0.0
        END DO

      IF (N .LT. 1) THEN
        write(kStdErr,*) 'Huh ? negative number of cloud layers'
        CALL DoStop
        END IF 

c -------------- cloud has only one layer ---------------------------------
      IF (N .eq. 1) THEN             !bloody simple  
        iLay    = iLocalCldTop
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))

        CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1)  !!total depth = 0.0
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau2,raTb,raBeta)
        !do the stuff for the arbitrary angles, for up going radiation 
        CALL asymcoeffs_solar(
     $              raAsym,raW0,raTb,raRadBb,raRadBt,raTau2,raBeta,
     $              radSolarCld,raCumSum2,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)

c        print *,'          '
c        print *,'asym solar Asym   ',(raAsym(iF),iF=1,2)
c        print *,'asym solar W0     ',(raW0(iF),iF=1,2)
c        print *,'asym solar Tb     ',(raTb(iF),iF=1,2)
c        print *,'asym solar Radbot ',(raRadBb(iF),iF=1,2)
c        print *,'asym solar Radtop ',(raRadBt(iF),iF=1,2)
c        print *,'asym solar Tau2   ',(raTau2(iF),iF=1,2)
c        print *,'asym solar Beta   ',(raBeta(iF),iF=1,2)
c        print *,'asym solar raA    ',(raA(iF),iF=1,2)
c        print *,'asym solar raB    ',(raB(iF),iF=1,2)
c        print *,'asym solar raDET  ',(raDet(iF),iF=1,2)
c        print *,'asym solar raE    ',(raE(iF),iF=1,2)
c        print *,'asym solar raF    ',(raF(iF),iF=1,2)
c        print *,'asym solar raG    ',(raG(iF),iF=1,2)
c        print *,'asym solar raH    ',(raH(iF),iF=1,2)
c        print *,'asym solar raKp   ',(raKp(iF),iF=1,2)
c        print *,'asym solar raKm   ',(raKm(iF),iF=1,2)
c        print *,'asym solar raAp   ',(raAp(iF),iF=1,2)
c        print *,'asym solar raAm   ',(raAm(iF),iF=1,2)

        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrUp2,raReUp2,raEmissUp2,raSunUp2)

c        print *,'          '
c        print *,'tre T   ',(raTrUp2(iF),iF=1,2)
c        print *,'tre R   ',(raReUp2(iF),iF=1,2)
c        print *,'tre E   ',(raEmissUp2(iF),iF=1,2)

        !do this for completeness
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrDown2,raReDown2,raEmissDown2,raSunDown2)

        DO iFr = 1,kMaxPts
          raTau12(iFr)       = raTau2(iFr)
          raSunDown12(iFr)   = raSunDown2(iFr)
          raTrDown12(iFr)    = raTrDown2(iFr)
          raReDown12(iFr)    = raReDown2(iFr)
          raEmissDown12(iFr) = raEmissDown2(iFr)
          raSunUp12(iFr)   = raSunUp2(iFr)
          raTrUp12(iFr)    = raTrUp2(iFr)
          raReUp12(iFr)    = raReUp2(iFr)
          raEmissUp12(iFr) = raEmissUp2(iFr)
          END DO
        ENDIF

c -------------- cloud has many layers ---------------------------------
      IF (N .GT. 1) THEN             !have to add the layers together
        !because the sun takes away up down symmetry, we have to do this
        !from the top down
        !do the stuff for the arbitrary angles, for CLOUD LAYER 2
        iLay    = iLocalCldTop 
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
        CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1) !!total depth = 0.0
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau2,raTb,raBeta)
        CALL asymcoeffs_solar(raAsym,raW0,raTb,raRadBb,raRadBt,raTau2,raBeta,
     $              radSolarCld,raCumSum2,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
        CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrUp2,raReUp2,raEmissUp2,raSunUp2)
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrDown2,raReDown2,raEmissDown2,raSunDown2)

        !do the stuff for the arbitrary angles, for CLOUD LAYER 1
        iLay = iLay - 1
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
        CALL AccumulateSolarDepth(raCumSum1,raTau2,+1) !!add upper layer depth
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau1,raTb,raBeta)
        CALL asymcoeffs_solar(raAsym,raW0,raTb,raRadBb,raRadBt,raTau1,raBeta,
     $              radSolarCld,raCumSum1,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
        CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrUp1,raReUp1,raEmissUp1,raSunUp1)
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrDown1,raReDown1,raEmissDown1,raSunDown1)

        !loop over the layers
        DO  iDp = 1,N-1
          !!! add layers together

          CALL addstar_arb_solar(
     $      raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1,
     $      raSunUp1,raSunDown1,raTau1,
     $      raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2,
     $      raSunUp2,raSunDown2,raTau2,
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raReUp12,raReDown12,raTrUp12,raTrDown12,
     $      raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12,
     $      raCumSum1,raCumSum2,mu_view,mu_sun)

          CALL addstar_solar(
     $      raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m, 
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raCumSum1,raCumSum2,mu_sun)

          iLay = iLay - 1
          IF (iDp .LT. (N-1)) THEN          !!!!!add on new layer stuff
            iL      = iaRadLayer(iLay)
            rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            mu_view = abs(rCos) 
            iBeta = MOD(iL,kProfLayer)
            IF (iBeta .EQ. 0) THEN
              iBeta = kProfLayer
              END IF
            rMPTemp = TEMP(iBeta)
            rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
            CALL AccumulateSolarDepth(raCumSum1,raTau1,+1) !!add on depths
            CALL LayerScatteringProp(
     $                    raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $                    raW0,raAsym,raTau1,raTb,raBeta)
            CALL asymcoeffs_solar(
     $              raAsym,raW0,raTb,raRadBb,raRadBt,raTau1,raBeta,
     $              radSolarCld,raCumSum1,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
            CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
            CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrUp1,raReUp1,raEmissUp1,raSunUp1)
            CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrDown1,raReDown1,raEmissDown1,raSunDown1)

            CALL update_streams_solar(
     $            raR2, raR2star, raT2, raT2star, raE2p, raE2m, raS2p, raS2m,
     $            raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
            CALL update_streams_arb_solar(            
     $              raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2,
     $              raEmissUp2, raEmissDown2, raSunUp2, raSunDown2,
     $              raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12,
     $              raEmissUp12,raEmissDown12,raSunup12,raSunDown12)

            END IF
          END DO
        ENDIF

      RETURN
      END

c************************************************************************
c this does the scattering radiative transfer for a uplook instrument
      SUBROUTINE Cloud_UpLook(iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,raWaves,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c input parameters
      INTEGER iNumLayer                    !number of layers in atm
      INTEGER iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
      INTEGER iaRadLayer(kProfLayer)       !atmosphere layering
      REAL raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
      REAL TEMP(MAXNZ)                     !temperature profile (levels)
      REAL raWaves(kMaxPts)                !wavenumbers
      REAL radSolarCld(kMaxPts)            !solar intensity at top of cloud
      !these next three are self explanatory
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL mu_sun                          !solar angle
c output parameters
      REAL mu_view                         !view angle at lowest cloud layer
      !these next few are self explanatory : optical depth, and the
      !cumulative up/down transmission, reflection, emission, solar trans
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)
      ! these are the lowest cloud layer asymmetry and single scattering
      REAL raW0(kMaxPts),raAsym(kMaxPts)

c local variables 
      INTEGER N,iFr,iLay,iL,iBeta,iDp,MP2LAY
      REAL rCos,rMPTemp

c general coeffs for the layers
      REAL rTT,rSunAngle,rBeta,rZeta,raBoo(kMaxPts)
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raTb(kMaxPts),raTt(kMaxPts)
      REAL raBeta(kMaxPts),raSun0(kMaxPts)
      REAL raA(kMaxPts),raB(kMaxPts),raE(kMaxPts),raF(kMaxPts)
      REAL raG(kMaxPts),raH(kMaxPts)
      REAL raKp(kMaxPts),raKm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts)
      
c arbitrary angle stuff
      REAL raTau1(kMaxPts),raTau2(kMaxPts)
      REAL raTrUp1(kMaxPts),raReUp1(kMaxPts),raEmissUp1(kMaxPts)
      REAL raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
      REAL raTrDown1(kMaxPts),raReDown1(kMaxPts),raEmissDown1(kMaxPts)
      REAL raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
      REAL raSunUp1(kMaxPts),raSunUp2(kMaxPts)
      REAL raSunDown1(kMaxPts),raSunDown2(kMaxPts)

c stream angle stuff
      REAL raT1(kMaxPts),raT1star(kMaxPts),raT2(kMaxPts),raT2star(kMaxPts)
      REAL raR1(kMaxPts),raR1star(kMaxPts),raR2(kMaxPts),raR2star(kMaxPts)
      REAL raE1p(kMaxPts),raE1m(kMaxPts),raE2p(kMaxPts),raE2m(kMaxPts)
      REAL raF1p(kMaxPts),raF1m(kMaxPts),raF2p(kMaxPts),raF2m(kMaxPts)
      REAL raT12(kMaxPts),raT12star(kMaxPts)
      REAL raR12(kMaxPts),raR12star(kMaxPts)
      REAL raE12p(kMaxPts),raE12m(kMaxPts)
      REAL raS1p(kMaxPts),raS2p(kMaxPts),raS12p(kMaxPts)
      REAL raS1m(kMaxPts),raS2m(kMaxPts),raS12m(kMaxPts)
      REAL raCumSum1(kMaxPts),raCumSum2(kMaxPts),raDet(kMaxPts)


      N = iLocalCldBot - iLocalCldTop + 1

      DO iFr = 1,kMaxPts
        raCumSum1(iFr) = 0.0
        raCumSum2(iFr) = 0.0
        END DO

      IF (N .LT. 1) THEN
        write(kStdErr,*) 'Huh ? negative number of cloud layers'
        CALL DoStop
        END IF 

c -------------- cloud has only one layer ---------------------------------
      IF (N .eq. 1) THEN             !bloody simple  
        iLay    = iLocalCldTop
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))

        CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1)  !!total depth = 0.0
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau2,raTb,raBeta)
        !do the stuff for the arbitrary angles, for up going radiation 
        CALL asymcoeffs_solar(
     $              raAsym,raW0,raTb,raRadBb,raRadBt,raTau2,raBeta,
     $              radSolarCld,raCumSum2,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrDown2,raReDown2,raEmissDown2,raSunDown2)
        !do this for completeness
        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrUp2,raReUp2,raEmissUp2,raSunUp2)
        DO iFr = 1,kMaxPts
          raTau12(iFr)       = raTau2(iFr)
          raSunDown12(iFr)   = raSunDown2(iFr)
          raTrDown12(iFr)    = raTrDown2(iFr)
          raReDown12(iFr)    = raReDown2(iFr)
          raEmissDown12(iFr) = raEmissDown2(iFr)
          raSunUp12(iFr)   = raSunUp2(iFr)
          raTrUp12(iFr)    = raTrUp2(iFr)
          raReUp12(iFr)    = raReUp2(iFr)
          raEmissUp12(iFr) = raEmissUp2(iFr)
          END DO
        ENDIF

c        !!g,w0,RadBb,RadBt,Tb,Tt,tau0,mu_view,mu_sun
c        print *,'---> ',raAsym(1),raW0(1),radtot(raWaves(1),raRadBb(1)),'<---'
c        print *,'---> ',radtot(raWaves(1),raRadBt(1)) ,' <---'
c        print *,'---> ',rMPTemp,rMPTemp,raTau2(1),' <----'
c        print *,'===> ',radtot(raWaves(1),raDiffuseInten(1)) ,' <==='

c -------------- cloud has many layers ---------------------------------
      IF (N .GT. 1) THEN             !have to add the layers together
        !because the sun takes away up down symmetry, we have to do this
        !from the top down
        !do the stuff for the arbitrary angles, for CLOUD LAYER 2
        iLay    = iLocalCldTop 
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
        CALL AccumulateSolarDepth(raCumSum2,raCumSum1,-1) !!total depth = 0.0
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau2,raTb,raBeta)
        CALL asymcoeffs_solar(raAsym,raW0,raTb,raRadBb,raRadBt,raTau2,raBeta,
     $              radSolarCld,raCumSum2,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
        CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raT2,raT2star,raR2,raR2star,raE2p,raE2m,raS2p,raS2m)
        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrUp2,raReUp2,raEmissUp2,raSunUp2)
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau2,mu_view,mu_sun,raAsym,raW0,raCumSum2,
     $              radSolarCld,raTb,
     $              raTrDown2,raReDown2,raEmissDown2,raSunDown2)

        !do the stuff for the arbitrary angles, for CLOUD LAYER 1
        iLay = iLay + 1        !!!!!!+1 becuz of wierd ordering
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        iBeta = MOD(iL,kProfLayer)
        IF (iBeta .EQ. 0) THEN
          iBeta = kProfLayer
          END IF
        rMPTemp = TEMP(iBeta)
        rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
        CALL AccumulateSolarDepth(raCumSum1,raTau2,+1) !!add upper layer depth
        CALL LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau1,raTb,raBeta)
        CALL asymcoeffs_solar(raAsym,raW0,raTb,raRadBb,raRadBt,raTau1,raBeta,
     $              radSolarCld,raCumSum1,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
        CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
        CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrUp1,raReUp1,raEmissUp1,raSunUp1)
        CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrDown1,raReDown1,raEmissDown1,raSunDown1)

        !loop over the layers
        DO  iDp = 1,N-1
          !!! add layers together

          CALL addstar_arb_solar(
     $      raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1,
     $      raSunUp1,raSunDown1,raTau1,
     $      raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2,
     $      raSunUp2,raSunDown2,raTau2,
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raReUp12,raReDown12,raTrUp12,raTrDown12,
     $      raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12,
     $      raCumSum1,raCumSum2,mu_view,mu_sun)

          CALL addstar_solar(
     $      raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m, 
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raCumSum1,raCumSum2,mu_sun)

          iLay = iLay + 1
          IF (iDp .LT. (N-1)) THEN          !!!!!add on new layer stuff
            iL      = iaRadLayer(iLay)
            rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            mu_view = abs(rCos) 
            iBeta = MOD(iL,kProfLayer)
            IF (iBeta .EQ. 0) THEN
              iBeta = kProfLayer
              END IF
            rMPTemp = TEMP(iBeta)
            rBeta = log(TEMP(iBeta+1)/TEMP(iBeta))
            CALL AccumulateSolarDepth(raCumSum1,raTau1,+1) !!add on depths
            CALL LayerScatteringProp(
     $                    raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $                    raW0,raAsym,raTau1,raTb,raBeta)
            CALL asymcoeffs_solar(
     $              raAsym,raW0,raTb,raRadBb,raRadBt,raTau1,raBeta,
     $              radSolarCld,raCumSum1,mu_sun,
     $              raA,raB,raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm)
            CALL t_r_e_streams_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)
            CALL t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrUp1,raReUp1,raEmissUp1,raSunUp1)
            CALL t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,
     $              radSolarCld,raTb,
     $              raTrDown1,raReDown1,raEmissDown1,raSunDown1)

            CALL update_streams_solar(
     $            raR2, raR2star, raT2, raT2star, raE2p, raE2m, raS2p, raS2m,
     $            raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
            CALL update_streams_arb_solar(            
     $              raTau2, raReUp2, raReDown2, raTrUp2, raTrDown2,
     $              raEmissUp2, raEmissDown2, raSunUp2, raSunDown2,
     $              raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12,
     $              raEmissUp12,raEmissDown12,raSunup12,raSunDown12)

            END IF
          END DO
        ENDIF
 
      RETURN
      END

c************************************************************************
c this finds out if we use complete HG function, or first order approx
      SUBROUTINE FindIHG(avg,iHG,iArbHG,raW0)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iHG,iArbHG
      REAL raW0(kMaxPts)
      REAL avg

      avg = (raW0(1)+raW0(kMaxPts))/2.0

cc      IF ((kSolar .GE. 0) .AND. (avg .LE. 0.75)) THEN
cc        !!use full HG phase function
cc        iHG    = 1
cc        iArbHG = 1
cc        !!use full HG phase function for twostream only
cc        iHG    = 1
cc        iArbHG = 1
cc      ELSE
cc        !!use first order HG phase function
cc        iHG    = -1
cc        iArbHG = -1
cc        END IF

cc ******* this was sept 2001 *********
      IF (kSolar .GE. 0) THEN
        !!use full HG phase function
        iHG    = 1
        iArbHG = 1
      ELSE
        !!use first order HG phase function
        iHG    = -1
        iArbHG = -1
        END IF
cc ******* this was sept 2001 *********

cc ******* try this dec 2001 *********
      !!use full HG phase function
      iHG    = 1
      iArbHG = 1
cc ******* try this dec 2001 *********

      RETURN
      END

c************************************************************************
c this subroutine computes the coefficients
      SUBROUTINE asymcoeffs_solar(raAsym,raW0,raBo,raRadBb,raRadBt,
     $            raTau,raBeta,radSolarCld,raCumSum,mu_sun,
     $            raA,raB,raDet,raE,raF,raG,raH,rakp,rakm,raAp,raAm)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c really we do not need raA or raB, so we do not really need RadBb or RadBt
c so we do not need to compute rho or theta

c input parameters
      REAL raAsym(kMaxPts),raW0(kMaxPts)    !asym coeff, single scatter alb
      REAL raBeta(kMaxPts)                  !1/k * log(T(n+1)/T(n))
      REAL raBo(kMaxPts)                    !bottom level temperature
      REAL raRadBb(kMaxPts),raRadBt(kmaxPts)!planck, top,bottom rad
      REAL raTau(kMaxPts)             !total extinction of layer
      REAL radSolarCld(kMaxPts),mu_sun      !sun intensity at cloud top, angle
      REAL raCumSum(kMaxPts)          !adjusts sun intensity as necessary
c output parameters
      REAL mu_pm                             !two stream angle
      REAL raDet(kMaxPts),raA(kMaxPts),raB(kMaxPts)  !!det is more important!
      REAL raE(kMaxPts),raF(kMaxPts)  !boundary condtions for soln
      REAL raG(kMaxPts),raH(kMaxPts)  !boundary condtions for soln
      REAL rakp(kMaxPts),rakm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff

c local variables
      INTEGER iF,iHG,iArbHG
      REAL b,det,alpha,theta,rho,avg,betaa,angle,p_plus,p_minus,hg2,psi,rSun
      REAL avgmin

      avgmin = 1.0e-6
      avgmin = 1.0e-3

      mu_pm = 1/sqrt(3.0)
      mu_sun = abs(mu_sun)

      CALL FindIHG(avg,iHG,iArbHG,raW0)

      IF ((avg .ge. avgmin) .AND. (kSolar .GE. 0)) THEN
        DO iF = 1,kMaxPts
          betaa = raBeta(iF)

          b = (1-raAsym(iF))/2.0
          raKp(iF) = 1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF))) 
          raKm(iF) =-1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
 
          alpha = (raW0(iF)*(1-b)-1) 
          raAp(iF) = -(mu_pm/(raW0(iF)*b))*(raKp(iF) + alpha/mu_pm)
          raAm(iF) = -(mu_pm/(raW0(iF)*b))*(raKm(iF) + alpha/mu_pm)

          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
          raE(iF) = raBo(iF)*(1-raW0(iF))/det
          raE(iF) = raE(iF)*(betaa*mu_pm + alpha - raW0(iF)*b)
          raF(iF) = raBo(iF)*(1-raW0(iF))/det
          raF(iF) = raF(iF)*(-betaa*mu_pm + alpha - raW0(iF)*b)

          rSun = radSolarCld(iF) * exp(-raCumSum(iF)/mu_sun)
          IF (iHG .GT. 0) THEN
            angle = mu_pm/mu_sun
            det = angle**2 - alpha**2 + (raW0(iF)*b)**2
            p_plus  = hg2(-abs(mu_sun),+abs(mu_pm),raAsym(iF)) 
            p_minus = hg2(-abs(mu_sun),-abs(mu_pm),raAsym(iF)) 
            raG(iF) = (alpha+angle)*p_plus - raW0(iF)*b*p_minus 
            raG(iF) = raW0(iF)*rSun/4/det * raG(iF)
            raH(iF) = raW0(iF)*b*(-p_plus) + (alpha-angle)*p_minus 
            raH(iF) = raW0(iF)*rSun/4/det * raH(iF) 

          ELSE
            angle = mu_pm/mu_sun 
            psi = 3*raAsym(iF)*mu_sun*mu_pm 
            det = (mu_pm/mu_sun)**2 - alpha**2 + (raW0(iF)*b)**2 
            raG(iF) = raW0(iF)*rSun/4/det * ((-alpha+raW0(iF)*b) +
     $                           psi*(alpha+raW0(iF)*b)+angle*(psi-1)) 
            raH(iF) = raW0(iF)*rSun/4/det * ((-alpha+raW0(iF)*b) - 
     $                           psi*(alpha+raW0(iF)*b)+angle*(psi+1)) 
            raG(iF) = -raG(iF) 
            raH(iF) = -raH(iF)
            END IF

          theta = raRadBb(iF) - raE(iF) - raG(iF)*exp(-raTau(iF)/mu_sun)
          rho   = raRadBt(iF) - raF(iF)*exp(betaa*raTau(iF)) - raH(iF) 
 
          raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) -
     $                raAm(iF)*exp(raKp(iF)*raTau(iF)) 

c these last two are really not needed, but what the heck!!!!!
          raA(iF) = ( theta*exp(raKm(iF)*raTau(iF)) - rho*raAm(iF)) 
          raA(iF) = raA(iF)/raDet(iF)
          raB(iF) = (-theta*exp(raKp(iF)*raTau(iF)) + rho*raAp(iF)) 
          raB(iF) = raB(iF)/raDet(iF) 
          END DO

      ELSEIF ((avg .ge. avgmin) .AND. (kSolar .LT. 0)) THEN
        DO iF = 1,kMaxPts
          betaa = raBeta(iF)

          b = (1-raAsym(iF))/2.0
          raKp(iF) = 1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF))) 
          raKm(iF) =-1/mu_pm * sqrt((1-raW0(iF))*(1-raW0(iF)*raAsym(iF)))
 
          alpha = (raW0(iF)*(1-b)-1) 
          raAp(iF) = -(mu_pm/(raW0(iF)*b))*(raKp(iF) + alpha/mu_pm)
          raAm(iF) = -(mu_pm/(raW0(iF)*b))*(raKm(iF) + alpha/mu_pm)

          det = (raW0(iF)*b)**2 - alpha*alpha + (betaa*mu_pm)**2
          raE(iF) = raBo(iF)*(1-raW0(iF))/det
          raE(iF) = raE(iF)*(betaa*mu_pm + alpha - raW0(iF)*b)
          raF(iF) = raBo(iF)*(1-raW0(iF))/det
          raF(iF) = raF(iF)*(-betaa*mu_pm + alpha - raW0(iF)*b)

          raG(iF) = 0.0
          raH(iF) = 0.0 

          theta = raRadBb(iF) - raE(iF) - raG(iF)*exp(-raTau(iF)/mu_sun)
          rho   = raRadBt(iF) - raF(iF)*exp(betaa*raTau(iF)) - raH(iF) 
 
          raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) -
     $                raAm(iF)*exp(raKp(iF)*raTau(iF)) 

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
          rakp(iF) = 1/mu_pm 
          rakm(iF) =-1/mu_pm 
 
          raAp(iF) = 0.0
          raAm(iF) = 1.0 
 
          raE(iF) = raBo(iF)/(1+betaa*mu_pm) 
          raF(iF) = raBo(iF)/(1-betaa*mu_pm) 
 
          raG(iF) = 0.0 
          raH(iF) = 0.0 
 
          theta = raRadBb(iF) - raE(iF)
          rho   = raRadBt(iF) - raF(iF)*exp(raTau(iF)*betaa) 
          raDet(iF) = raAp(iF)*exp(raKm(iF)*raTau(iF)) -
     $                raAm(iF)*exp(raKp(iF)*raTau(iF)) 
          !!! these last two are really not needed, but what the heck!!!!!
          raA(iF) = rho*exp(raKm(iF)*raTau(iF)) 
          raB(iF) = theta 

          END DO
        END IF  

      RETURN
      END

c************************************************************************
c this subroutine computes upwards R,T,Emiss for arbitrary view angle
      SUBROUTINE t_r_e_arb_up_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau,mu_view,mu_sun,raAsym,raW0,raCumSum,radSolarCld,raBo,
     $              tr_up,re_up,emiss_up,sun_up)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
      REAL mu_view,mu_sun,raTau(kMaxPts)   !view,solar angle and extinction
      REAL raCumSum(kMaxPts),radSolarCld(kMaxPts)!direct solar extinction 
      REAL raAsym(kMaxPts),raW0(kMaxPts)   !asymmetry and albedo
      REAL raDet(kMaxPts)                  !determinant 
      REAL raE(kMaxPts),raF(kMaxPts)       !used for emission
      REAL raG(kMaxPts),raH(kMaxPts)       !used for sun
      REAL rakp(kMaxPts),rakm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff
      REAL raBeta(kMaxPts),raBo(kMaxPts)   !temp of lower level, variation
      
c output parameters
      REAL tr_up(kMaxPts),re_up(kMaxPts),emiss_up(kMaxPts),sun_up(kMaxPts)

c local variables
      INTEGER iF
      REAL mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
      REAL rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
      REAL eup,fup,gup,hup,sup,hg2,betaa,mu,avgmin

      INTEGER iArbHG,iHG

      avgmin = 1.0e-6
      avgmin = 1.0e-3

      mu_pm   = 1/sqrt(3.0)
      mu      = abs(mu_view)
      mu_sun  = abs(mu_sun)

      CALL FindIHG(avg,iHG,iArbHG,raW0)

      IF (avg .ge. avgmin) THEN
        DO iF = 1,kMaxPts
          betaa = raBeta(iF)
          rSun = radSolarCld(iF)*exp(-raCumSum(iF)/mu_sun)
          det = raDet(iF)
          
          b = (1-raAsym(iF))/2 
          alpha = (raW0(iF)*(1-b)-1) 
 
          alpha_p = raW0(iF)/2/mu * (1 + 3*raAsym(iF)*mu_pm*mu)     
          alpha_m = raW0(iF)/2/mu * (1 - 3*raAsym(iF)*mu_pm*mu)     

          vp = 1/mu + raKp(iF)         
          vm = 1/mu + raKm(iF) 
          vb = betaa + 1/mu              
          vs = 1/mu_sun + 1/mu   
 
          !evaluate at raTau(iF) 
          evm = exp(vm*raTau(iF)) - 1 
          evp = exp(vp*raTau(iF)) - 1 
          evb = exp(vb*raTau(iF)) - 1     
          evs = (exp(vs*raTau(iF))-1)/vs 
 
          omega = evb/vb  
          gamma_a = (alpha_p * raAp(iF) + alpha_m) * evp/vp 
          gamma_b = (alpha_p * raAm(iF) + alpha_m) * evm/vm 
 
          kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(iF)*b)**2
          kappa = (betaa*mu_pm + alpha - raW0(iF)*b)/kappa 

          tr_up(iF) = (gamma_a*exp(raKm(iF)*raTau(iF)) - 
     $                 gamma_b*exp(raKp(iF)*raTau(iF)))/det
          re_up(iF) = (raAp(iF) *gamma_b - raAm(iF)*gamma_a)/det

          eup = -tr_up(iF) + omega*(alpha_p + 1/(mu*kappa)) 
          fup = -re_up(iF) + omega*alpha_m 
          emiss_up(iF) = raE(iF)*eup + raF(iF)*fup 

          IF (kSolar .GE. 0) THEN
            gup = exp(-raTau(iF)/mu_sun)*(evs*alpha_p - tr_up(iF)) 
            hup = exp(-raTau(iF)/mu_sun)*evs*alpha_m - re_up(iF)    
            IF (iarbhg .GT. 0) THEN
              sup = exp(-raTau(iF)/mu_sun)*raW0(iF)*evs/4.0/mu * 
     $              hg2(-abs(mu_sun),abs(mu),raAsym(iF)) 
            ELSE 
              sup = exp(-raTau(iF)/mu_sun)*raW0(iF)*evs/4.0/mu *
     $               (1-3*raAsym(iF)*mu*mu_sun) 
              END IF
            sun_up(iF) = (raG(iF)*gup + raH(iF)*hup + rSun*sup)/rSun

          ELSE 
            sun_up(iF) = 0.0
            END IF

          END DO

      ELSEIF ((avg .le. avgmin) .and. (kTemperVary .LT. 0)) THEN
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

      ELSEIF ((avg .lt. avgmin) .and. (kTemperVary .GT. 0)) THEN
        DO iF = 1,kMaxPts
          !!!!!!!!!no scattering, clear sky solns
          betaa = raBeta(iF)

          ! evaluate at tau0 
          tr_up(iF) = 0.0        !remember, no 2stream radiation knocked here 
          re_up(iF) = 0.0        !remember, no 2stream radiation knocked here 

          eup = raBo(iF)/(1+betaa*mu)*
     $           (exp(betaa*raTau(iF))-exp(-raTau(iF)/mu)) 
          eup = eup*exp(raTau(iF)/mu) 

          fup = 0.0 
          gup = 0.0          !remember, no scattering of solar beam to here
          hup = 0.0          !remember, no scattering of solar beam to here 

          emiss_up(iF) = eup 
          sun_up(iF)   = 0.0 
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine computes downwards R,T,Emiss for arbitrary view angle
      SUBROUTINE t_r_e_arb_down_solar(
     $              raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $              raTau,mu_view,mu_sun,raAsym,raW0,raCumSum,
     $              radSolarCld,raBo,
     $              tr_dn,re_dn,emiss_dn,sun_dn)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
      REAL mu_view,mu_sun,raTau(kMaxPts)  !view,solar angle and extinction
      REAL raCumSum(kMaxPts),radSolarCld(kMaxPts)!direct solar extinction 
      REAL raAsym(kMaxPts),raW0(kMaxPts)   !asymmetry and albedo
      REAL raDet(kMaxPts)                  !determinant 
      REAL raE(kMaxPts),raF(kMaxPts)       !used for emission
      REAL raG(kMaxPts),raH(kMaxPts)       !used for sun
      REAL rakp(kMaxPts),rakm(kMaxPts),raAp(kMaxPts),raAm(kMaxPts) !eigen stuff
      REAL raBeta(kMaxPts),raBo(kMaxPts)    !temp of lower level, variation
      
c output parameters
      REAL tr_dn(kMaxPts),re_dn(kMaxPts),emiss_dn(kMaxPts),sun_dn(kMaxPts)

c local variables
      INTEGER iF
      REAL mu_pm,det,b,alpha,alpha_p,alpha_m,vp,vm
      REAL rSun,avg,evm,evp,evb,evs,omega,gamma_a,gamma_b,vb,vs,kappa
      REAL edn,fdn,gdn,hdn,sdn,hg2,betaa,mu,avgmin
      INTEGER iArbHG,iHG

      avgmin = 1.0e-6
      avgmin = 1.0e-3

      mu_pm   = 1/sqrt(3.0)
      mu      = abs(mu_view)
      mu_sun  = abs(mu_sun)

      CALL FindIHG(avg,iHG,iArbHG,raW0)

      IF (avg .ge. avgmin) THEN
        DO iF = 1,kMaxPts
          betaa = raBeta(iF)
          rSun = radSolarCld(iF)*exp(-raCumSum(iF)/mu_sun)
          det = raDet(iF)

          b = (1-raAsym(iF))/2 
          alpha = (raW0(iF)*(1-b)-1) 
 
          alpha_p = raW0(iF)/2/mu * (1 - 3*raAsym(iF)*mu_pm*mu)     
          alpha_m = raW0(iF)/2/mu * (1 + 3*raAsym(iF)*mu_pm*mu)     

          vp = -1/mu + raKp(iF)         
          vm = -1/mu + raKm(iF) 
          vb = betaa - 1/mu              
          vs = 1/mu_sun - 1/mu   
 
          !evaluate at 0.0
          evm = exp(vm*raTau(iF)) - 1 
          evp = exp(vp*raTau(iF)) - 1 
          evb = exp(vb*raTau(iF)) - 1     
          evs = (exp(vs*raTau(iF))-1)/vs 
 
          omega = evb/vb  
          gamma_a = (alpha_p * raAp(iF) + alpha_m) * evp/vp 
          gamma_b = (alpha_p * raAm(iF) + alpha_m) * evm/vm 
 
          kappa = (betaa*mu_pm)**2 - alpha*alpha + (raW0(iF)*b)**2
          kappa = (betaa*mu_pm + alpha - raW0(iF)*b)/kappa 

          re_dn(iF) = (gamma_a*exp(raKm(iF)*raTau(iF)) - 
     $                 gamma_b*exp(raKp(iF)*raTau(iF)))/det
          tr_dn(iF) = (raAp(iF) *gamma_b - raAm(iF)*gamma_a)/det

          edn = -re_dn(iF) + omega*(alpha_p + 1/(mu*kappa)) 
          fdn = -tr_dn(iF)*exp(betaa*raTau(iF)) + omega*alpha_m 
          emiss_dn(iF) = raE(iF)*edn + raF(iF)*fdn 

          IF (kSolar .GE. 0) THEN
            gdn = exp(-raTau(iF)/mu_sun)*(evs*alpha_p - re_dn(iF)) 
            hdn = exp(-raTau(iF)/mu_sun)*evs*alpha_m - tr_dn(iF) 
   
            IF (iarbhg .GT. 0) THEN
              sdn = exp(-raTau(iF)/mu_sun)*raW0(iF)*evs/4.0/mu * 
     $              hg2(-abs(mu_sun),-abs(mu),raAsym(iF)) 
            ELSE 
              sdn = exp(-raTau(iF)/mu_sun)*raW0(iF)*evs/4.0/mu *
     $               (1+3*raAsym(iF)*mu*mu_sun) 
              END IF
            sun_dn(iF) = (raG(iF)*gdn + raH(iF)*hdn + rSun*sdn)/rSun

          ELSE 
            sun_dn(iF) = 0.0
            END IF

          END DO

      ELSEIF ((avg .lt. avgmin) .AND. (kTemperVary .LT. 0)) THEN 
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

      ELSEIF ((avg .lt. avgmin) .AND. (kTemperVary .GT. 0)) THEN 
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
        tr_dn(iF)    = tr_dn(iF) * exp(raTau(iF)/mu)
        re_dn(iF)    = re_dn(iF) * exp(raTau(iF)/mu)
        emiss_dn(iF) = emiss_dn(iF) * exp(raTau(iF)/mu)
        sun_dn(iF)   = sun_dn(iF) * exp(raTau(iF)/mu)
        END DO

      RETURN
      END

c************************************************************************
c this is for the twostream radiation
      SUBROUTINE t_r_e_streams_solar(
     $         raDet,raE,raF,raG,raH,raKp,raKm,raAp,raAm,raBeta,
     $         raTau1,mu_view,mu_sun,raAsym,raW0,raCumSum1,radSolarCld,raTb,
     $         raT1,raT1star,raR1,raR1star,raE1p,raE1m,raS1p,raS1m)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output params
      REAL raT1(kMaxPts),raT1star(kMaxPts)    !transmissions
      REAL raR1(kMaxPts),raR1star(kMaxPts)   !reflections
      REAL raE1p(kMaxPts),raE1m(kMaxPts)      !emissions
      REAL raS1p(kMaxPts),raS1m(kMaxPts)      !sun stuff

c input params
      REAL mu_sun,mu_view                     !solar,view angle (not needed)
      REAL raAsym(kMaxPts)                    !asymmetry
      REAL radSolarCld(kMaxPts),raCumSum1(kMaxPts)  !solar stuff
      REAL raTb(kMaxPts)                      !radiance at bottom
      REAL raW0(kMaxPts),raBeta(kMaxPts)      !albedo and temperature variation
      REAL raKp(kMaxPts),raKm(kMaxPts)        !eigenvalues
      REAL raAp(kMaxPts),raAm(kMaxPts)        !eigenvectors
      REAL raTau1(kMaxPts),raDet(kMaxPts)     !optical depth, determinant 
      REAL raE(kMaxPts),raF(kMaxPts)          !up and down emission
      REAL raG(kMaxPts),raH(kMaxPts)          !up and down sun stuff
      

c local variables 
      INTEGER iF,iHG,iArbHG
      REAL avg,det,bw1,bw2,hg2,rSun,tau,mu_pm,tau0
      REAL tr_up,re_up,eup,fup,gup,hup,sun_up,emiss_up
      REAL tr_down,re_down,edn,fdn,gdn,hdn,sun_down,emiss_down,betaa,avgmin

      avgmin = 1.0e-6
      avgmin = 1.0e-3

      CALL FindIHG(avg,iHG,iArbHG,raW0)

      mu_pm   = 1/sqrt(3.0)
      mu_sun  = abs(mu_sun)

      IF (avg .gt. avgmin) THEN
        DO iF = 1,kMaxPts
          betaa = raBeta(iF)
          rSun = radSolarCld(iF)*exp(-raCumSum1(iF)/mu_sun)

          tau  = raTau1(iF) 
          tau0 = raTau1(iF)
          tr_up = raAp(iF)*exp(raKm(iF)*tau0)*exp(raKp(iF)*tau) - 
     $            raAm(iF)*exp(raKp(iF)*tau0)*exp(raKm(iF)*tau) 
          tr_up = tr_up/raDet(iF) 
          re_up = raAp(iF)*raAm(iF)*(exp(raKm(iF)*tau) - exp(raKp(iF)*tau)) 
          re_up = re_up/raDet(iF) 
          eup = exp(betaa*tau) - tr_up 
          fup = -re_up*exp(betaa*tau0) 
          gup = exp(-tau0/mu_sun)*(exp(tau/mu_sun) - tr_up) 
          hup = -re_up 
          emiss_up = raE(iF)*eup + raF(iF)*fup 
          IF (kSolar .GE. 0) THEN
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
          re_down = exp(raKm(iF)*tau0)*exp(raKp(iF)*tau) - 
     $              exp(raKp(iF)*tau0)*exp(raKm(iF)*tau)
          re_down = re_down/raDet(iF) 
          edn = -re_down 
          fdn = -tr_down*exp(betaa*tau0) + exp(betaa*tau) 
          emiss_down = raE(iF)*edn + raF(iF)*fdn 
          IF (kSolar .GE. 0) THEN           
            gdn = exp(-tau0/mu_sun)*(-re_down) 
            hdn = exp(-(tau0-tau)/mu_sun) - tr_down 
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
      END

c************************************************************************
c this subroutine adds together stuff for the arbitrary angles
      SUBROUTINE  addstar_arb_solar(
     $      raReUp1,raReDown1,raTrUp1,raTrDown1,raEmissUp1,raEmissDown1,
     $      raSunUp1,raSunDown1,raTau1,
     $      raReUp2,raReDown2,raTrUp2,raTrDown2,raEmissUp2,raEmissDown2,
     $      raSunUp2,raSunDown2,raTau2,
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raReUp12,raReDown12,raTrUp12,raTrDown12,
     $      raEmissUp12,raEmissDown12,raSunUp12,raSunDown12,raTau12,
     $      raCumSum1,raCumSum2,mu_view,mu_sun)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params
      !! layer stuff for the arb angles, layer n
      REAL raReUp1(kMaxPts),raReDown1(kMaxPts)
      REAL raTrUp1(kMaxPts),raTrDown1(kMaxPts)
      REAL raEmissUp1(kMaxPts),raEmissDown1(kMaxPts),ratau1(kMaxPts)
      REAL raSunUp1(kMaxPts),raSunDown1(kMaxPts)
      !! layer stuff for the arb angles, layer n+1
      REAL raReUp2(kMaxPts),raReDown2(kMaxPts)
      REAL raTrUp2(kMaxPts),raTrDown2(kMaxPts)
      REAL raEmissUp2(kMaxPts),raEmissDown2(kMaxPts),ratau2(kMaxPts)
      REAL raSunUp2(kMaxPts),raSunDown2(kMaxPts)
      !! layer stuff for the stream angles, layer n
      REAL raR1(kMaxPts),raT1(kMaxPts),raR1star(kMaxPts)
      REAL raT1star(kMaxPts)
      REAL raE1p(kMaxPts),raE1m(kMaxPts)
      REAL raS1p(kMaxPts),raS1m(kMaxPts)
      !! layer stuff for the stream angles, layer n+1
      REAL raR2(kMaxPts),raT2(kMaxPts),raR2star(kMaxPts)
      REAL raT2star(kMaxPts)
      REAL raE2p(kMaxPts),raE2m(kMaxPts)
      REAL raS2p(kMaxPts),raS2m(kMaxPts)
      !!arbitrary stuff
      REAL mu_sun,mu_view
      REAL raCumSum1(kMaxPts),raCumSum2(kMaxPts)
c output params
      !! layer stuff for the arb angles, layer n AND n+1
      REAL raReUp12(kMaxPts),raReDown12(kMaxPts)
      REAL raTrUp12(kMaxPts),raTrDown12(kMaxPts)
      REAL raEmissUp12(kMaxPts),raEmissDown12(kMaxPts),raTau12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)

c local variables
      INTEGER iF
      REAL mat1,mat2,mat3,mat4,mat5,mat6,mat7,mat11,mat12,mat21,mat22
      REAL mu,det,a1,a2,tau12,e1,e2

      mu = abs(mu_view)
      mu_sun = abs(mu_sun)
      
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
        a1 = raEmissUp1(iF) + 0.0 + exp(raTau1(iF)/mu)*raEmissUp2(iF) + 
     $        0.0 +  
     $        raTrUp1(iF)/det*(raR2star(iF)*raE1p(iF)+raE2m(iF)) +  
     $        0.0 +  
     $   raTrUp2(iF)/det*(raR1(iF)*raE2m(iF)+raE1p(iF))*exp(raTau1(iF)/mu) 
 
        a2 = raEmissDown2(iF) + exp(raTau2(iF)/mu)*raEmissDown1(iF) + 0.0 + 
     $         raReDown2(iF)/det*(raR1(iF)*raE2m(iF) + raE1p(iF)) +  
     $         0.0 +  
     $  raTrDown1(iF)/det*(raR2star(iF)*raE1p(iF)+raE2m(iF))*exp(raTau2(iF)/mu)
     $         + 0.0 
      
        raReUp12(iF)      = mat11       
        raTrUp12(iF)      = mat12 
        raTrDown12(iF)    = mat21       
        raReDown12(iF)    = mat22 
        raEmissUp12(iF)   = a1        
        raEmissDown12(iF) = a2 
 
        raTau12(iF)  = raTau1(iF) + raTau2(iF)     !total optical depth so far 
        END DO

      IF (kSolar .GE. 0) THEN           !!!do the solar part
        DO iF = 1,kMaxPts
          det = 1 - raR1(iF)*raR2star(iF)
          e1 = exp(-raCumsum1(iF)/mu_sun) 
          e2 = exp(-raCumsum2(iF)/mu_sun) 

          raSunUp12(iF) = raSunUp1(iF)*e1 +
     $                    raSunUp2(iF)*e2*exp(raTau1(iF)/mu) + 
     $                    1/det*(raReUp1(iF)*(e1*raR2star(iF)*raS1p(iF) + 
     $                    e2*raS2m(iF))) + 
     $                    exp(raTau1(iF)/mu)*raTrUp2(iF)/det *
     $                    (e2*raR1(iF)*raS2m(iF) + e1*raS1p(iF)) 
          raSunDown12(iF) = raSunDown2(iF)*e2 + 
     $                    raSunDown1(iF)*e1*exp(raTau2(iF)/mu) + 
     $                    1/det*(raReDown2(iF)*(e2*raR1(iF)*raS2m(iF) + 
     $                    e1*raS1p(iF))) +
     $                    exp(raTau2(iF)/mu)*raTrDown1(iF)/det *
     $                    (e1*raR2star(iF)*raS1p(iF) + e2*raS2m(iF)) 
         END DO
       END IF

      RETURN
      END

c************************************************************************
c this subroutine adds together the layers for the stream angles
      SUBROUTINE addstar_solar(
     $      raR12,raT12,raR12star,raT12star,raE12p,raE12m,raS12p,raS12m, 
     $      raR1,raT1,raR1star,raT1star,raE1p,raE1m,raS1p,raS1m, 
     $      raR2,raT2,raR2star,raT2star,raE2p,raE2m,raS2p,raS2m,
     $      raCumSum1,raCumSum2,mu_sun)
     
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params
      !! layer stuff for the stream angles, layer n
      REAL raR1(kMaxPts),raT1(kMaxPts),raR1star(kMaxPts)
      REAL raT1star(kMaxPts),raE1p(kMaxPts),raE1m(kMaxPts)
      REAL raS1p(kMaxPts),raS1m(kMaxPts)
      !! layer stuff for the stream angles, layer n+1
      REAL raR2(kMaxPts),raT2(kMaxPts),raR2star(kMaxPts)
      REAL raT2star(kMaxPts),raE2p(kMaxPts),raE2m(kMaxPts)
      REAL raS2p(kMaxPts),raS2m(kMaxPts)
c for the sun
      REAL raCumSum1(kMaxPts),raCumSum2(kMaxPts),mu_sun
c output params
      !! layer stuff for the stream angles, layer n AND n+1
      REAL raR12(kMaxPts),raT12(kMaxPts),raR12star(kMaxPts)
      REAL raT12star(kMaxPts),raE12p(kMaxPts),raE12m(kMaxPts)
      REAL raS12p(kMaxPts),raS12m(kMaxPts)
c local variables
      INTEGER iF
      REAL det,e1,e2

      mu_sun = abs(mu_sun)

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

      IF (kSolar .GE. 0) THEN  !!!!do the solar part
        DO iF = 1,kMaxPts
          det = 1 - raR1(iF)*raR2star(iF) 
          e1 = exp(-raCumSum1(iF)/mu_sun) 
          e2 = exp(-raCumSum2(iF)/mu_sun) 

          raS12p(iF) = (raT2star(iF)/det*raR1(iF)*raS2m(iF) + raS2p(iF))*e2 +
     $                  raT2star(iF)/det*raS1p(iF)*e1
          raS12m(iF) = (raT1(iF)/det*raR2star(iF)*raS1p(iF) + raS1m(iF))*e1 +
     $                  raT1(iF)/det*raS2m(iF)*e2
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine updates the stream angle added layers
      SUBROUTINE update_streams_solar(
     $          raR2,raR2star,raT2,raT2star,raE2p,raE2m,raS2p,raS2m,
     $          raR12,raR12star,raT12,raT12star,raE12p,raE12m,raS12p,raS12m)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c stream angle stuff
      REAL raT2(kMaxPts),raT2star(kMaxPts)
      REAL raR2(kMaxPts),raR2star(kMaxPts)
      REAL raE2p(kMaxPts),raE2m(kMaxPts)
      REAL raS2p(kMaxPts),raS2m(kMaxPts)
      REAL raT12(kMaxPts),raT12star(kMaxPts)
      REAL raR12(kMaxPts),raR12star(kMaxPts)
      REAL raE12p(kMaxPts),raE12m(kMaxPts)
      REAL raS12p(kMaxPts),raS12m(kMaxPts)

c local variables
      INTEGER iF
 
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
      END

c************************************************************************
c this subroutine updates the arb angle added layers
      SUBROUTINE update_streams_arb_solar(            
     $              raTau2,raReUp2,raReDown2,raTrUp2,raTrDown2,
     $              raEmissUp2,raEmissDown2,raSunUp2,raSunDown2,
     $              raTau12,raReUp12,raReDown12,raTrUp12,raTrDown12,
     $              raEmissUp12,raEmissDown12,raSunup12,raSunDown12)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c arbitrary angle stuff
      REAL raTau2(kMaxPts),raTau12(kMaxPts)
      REAL raTrUp2(kMaxPts),raReUp2(kMaxPts),raEmissUp2(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown2(kMaxPts),raReDown2(kMaxPts),raEmissDown2(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp2(kMaxPts),raSunDown2(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)

c local variables
      INTEGER iF

      !!!only need to update the upper layer 
      DO iF = 1,kMaxPts
        raTau2(iF)       = raTau12(iF)
        raReUp2(iF)      = raReUp12(iF) 
        raReDown2(iF)    = raReDown12(iF)
        raTrUp2(iF)      = raTrUp12(iF) 
        raTrDown2(iF)    = raTrDown12(iF) 
        raEmissUp2(iF)   = raEmissUp12(iF)
        raEmissDown2(iF) = raEmissDown12(iF) 
        raSunUp2(iF)     = raSunUp12(iF)
        raSunDown2(iF)   = raSunDown12(iF)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine sets the scattering properties of the current layer      
      SUBROUTINE LayerScatteringProp(
     $              raWaves,rMPTemp,iL,raaExt,raaScat,raaAsym,rBeta,
     $              raW0,raAsym,raTau,raBo,raBeta)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
      INTEGER iL                 !current MP number
      REAL rMPTemp               !current layer temperature
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts)
      REAL rBeta                 !log(T(n+1)/T(n))
c output parameters
      REAL raTau(kMaxPts)        !total layer optical depth (includes scatter)
      REAL raW0(kMaxPts),raAsym(kMaxPts)   !layer albedo, asymmetry
      REAL raBo(kMaxPts)         !lower level radiances 
      REAL raBeta(kMaxPts)       !beta = 1/tau log(T(n+1)/T(n))

c local variables
      INTEGER iFr
      REAL ttorad_local

      DO iFr = 1,kMaxPts
        raBo(iFr)       = ttorad_local(raWaves(iFr),rMPTemp)  
        raTau(iFr)      = raaExt(iFr,iL)
        raW0(iFr)       = raaScat(iFr,iL)/raaExt(iFr,iL)
        raAsym(iFr)     = raaAsym(iFr,iL)
        raBeta(iFr)     = 1/raTau(iFr) * rBeta
        END DO

      RETURN
      END

c************************************************************************
c this accumulates the optical depth from top of cloud to top of layer
      SUBROUTINE AccumulateSolarDepth(raCumSum,raX,iDo)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
      INTEGER iDo                 !1 to add on raX to raCumSum, -1 if not
      REAL raCumSum(kMaxPts),raX(kMaxPts)  !total optical depth from cloud top 

c local variables
      INTEGER iFr

      IF (iDo .GT. 0) THEN
        DO iFr = 1,kMaxPts
          raCumSum(iFr) = raCumSum(iFr) + raX(iFr)
          END DO
      ELSE
        DO iFr = 1,kMaxPts
          raCumSum(iFr) = 0.0
          END DO
        END IF
       

      RETURN
      END

c************************************************************************
c this function computes the Henyey Greenstein function, assuming the
c cos(phi1-phi2) factor = 1
      REAL Function hg2(mu1,mu2,g)

      IMPLICIT NONE

      REAL mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry

      REAL normB,mu0,yexact

      !! normB is normalisation of mu from -1 to 1 
      !! we also know that (1/2) integral P(-1,1) = 1 
      normB = 1/sqrt(1+g*g - 2*g) - 1/sqrt(1+g*g + 2*g)
      normB = (1-g*g)/g * normB

      !!!compute mu0 = cos ofangle between the two
      mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2) 
 
      yexact = (1 + g*g - 2*g*mu0) * sqrt(1 + g*g - 2*g*mu0)
      yexact = (1-g*g)/yexact
      yexact = yexact/normB * 2 
 
      hg2 = yexact

      RETURN
      END

c************************************************************************
c same as ttorad, except it might need to be promoted to double precision!

      REAL function ttorad_local(rf,rBT)
c rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rf = wavenumber, rI = intensity, rBT = brightness temp
      REAL rf,rI,rBT

c local variables
      REAL r1,r2,rPlanck
      INTEGER iInt
 
      r1=kPlanck1
      r2=kPlanck2

      rPlanck = exp(r2*rf/rBT) - 1 
      rI      = r1*(rf**3)/rPlanck

      ttorad_local = rI

      RETURN
      END

c************************************************************************
