c Copyright 2002  
c University of Maryland Baltimore County  
c All Rights Reserved 

c this has the lineshape I/O stuff, that is used by knonlte.f
c these include : qfcns, line broadening, lte strenghts, population ratios
c it also reads in HITRAN parameters and US standard Optical Depths

c************************************************************************
c                     these are the lineshape subroutines
c************************************************************************
c this calls the qtips fcns
      SUBROUTINE qfcn(dT,iGasID,iISO,dPartitionFcn,dMass)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iGasID,iISO
      DOUBLE PRECISION   dT             !!temperature
c output params
      DOUBLE PRECISION dPartitionFcn,dMass

c local vars
      DOUBLE PRECISION dH98,dH96
      INTEGER iQtips_H98,iDefault

      iDefault = +1   !default; use the most recent qtips values
      iQtips_H98     = -1   !use some other value
      iQtips_H98     = +1   !default; use the most recent qtips values

c      CALL qfcnH98(dT,iGasID,iISO,dH98,dMass)
c      CALL qfcnH96(dT,iGasID,iISO,dH96,dMass)
c      print *,dT,dH98/dH96

      IF (iQtips_H98 .NE. iDefault) THEN
        print *,'iQtips_H98,iDefault = ',iQtips_H98,iDefault
      END IF

      IF (iQtips_H98 .GT. 0) THEN  !!call default (most recent)
        CALL qfcnH98(dT,iGasID,iISO,dPartitionFcn,dMass)
      ELSE
        CALL qfcnH96(dT,iGasID,iISO,dPartitionFcn,dMass)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine finds the broadening coefficients
c sorta mimics run7co2
      SUBROUTINE Broad(iGasID,iNum,daLineShift,
     $           daW_For,daW_Self,daW_temp,dP,dPP,dLTE,daBroad)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iGasID,iNum
      DOUBLE PRECISION daW_For(kHITRAN),daW_Self(kHITRAN),daW_temp(kHITRAN)
      DOUBLE PRECISION dP,dPP,dLTE,daLineShift(kHITRAN)
c output params
      DOUBLE PRECISION daBroad(kHITRAN)

c local variables
      INTEGER iLines
      DOUBLE PRECISION d1,d2,dm

      DO iLines = 1,iNum
        d1 = daW_For(iLines)    !!foreign broadening coeff
        d2 = daW_Self(iLines)   !!self    broadening coeff
        dm = daW_temp(iLines)   !!temperature effects

        d1 = (dP - dPP) * d1 * ((296.0/dLTE) ** dm)

        IF ((d2 .LE. 1.0d-16) .AND. (iGasID .GT. 1)) THEN
          d2 = d1 !!make self coeff == air coeff
        ELSE IF ((d2 .LE. 1.0d-16) .AND. (iGasID .EQ. 1)) THEN
          d2 = d1 * 5.0 !!self coeff
        END IF

        IF (iGasID .EQ. 2) THEN
          d2 = (dPP) * d2 * ((296.0/dLTE) ** (0.685*1.0d0))
        ELSEIF (iGasID .NE. 2) THEN
          d2 = (dPP) * d2 * ((296.0/dLTE) ** dm)
        END IF

        daBroad(iLines) = (d1 + d2)
      END DO

      RETURN
      END

c************************************************************************
c find the linestrengths, at LTE!!!!!!!!!!!!!!!!!!!!!!!!!
c also find the QTIPS modifier, which should be 1 since all at LTE
      SUBROUTINE Strengths_lte_only(
     $                     iNum,iGasID,dGasAmt,daPartitionFcn,dLTE,
     $                     daLineShift,daELower,daStren296,
     $                     daStrenLTE,daQtipsFactor)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iNum,iGasID
      DOUBLE PRECISION dGasAmt,daPartitionFcn(kHITRAN),dLTE
      DOUBLE PRECISION daLineShift(kHITRAN),daElower(kHITRAN)
      DOUBLE PRECISION daStren296(kHITRAN)     !!straight from HITRAN at 296 K
c output params
      DOUBLE PRECISION daStrenLTE(kHITRAN),daQtipsFactor(kHITRAN)

c local variables
      INTEGER iLines
      DOUBLE PRECISION v0,s1,sb,se,s2,s3,TS,sq,jfac,dPartitionFcn

      jfac = 1.0d0
      TS = 296.0 * 1.0d0

      DO iLines = 1,iNum
        v0 = daLineShift(iLines)

        !!!!boltzmann factor
        s1 = kPlanck2 *(dLTE-TS)/(dLTE*TS)
        sb = dexp(s1 * daELower(iLines))

        !!!stimulated emission
        s2 = -kPlanck2 * v0/dLTE
        s3 = -kPlanck2 * v0/TS

        IF (DEXP(s3) .EQ. 1.0) THEN
          se = 1.0d0
        ELSE
          se = (1.0d0 - dexp(s2))/(1.0d0 - dexp(s3))
        END IF

        dPartitionFcn = daPartitionFcn(iLines)
        sq = dPartitionFcn
         
        daStrenLTE(iLines) = dGasAmt*daStren296(iLines)*sb*se*sq * kAvog

c no adjustment to SQ due to NLTE
        daQtipsFactor(iLines) = 1.0d0   

      END DO

      RETURN
      END

c************************************************************************
c find the linestrengths, at NLTE!!!!!!!!!!!!!!!!!!!!!!!!!
c   dPartitionFcn is from LTE qtips.m
c   vibqft        is from NLTE file, and is vib.qtips interpolated to profile
c In other words
      !!find the line strengths at *** LTE *****
      !!also find QTIPS modifier due to vibration partition function at NLTE
      !!so ultimately you have the *** NLTE *** strengths
      !!   *** so you have to be careful with AlphaBetaFactors ******
      SUBROUTINE Strengths_nlte_vibpartfcn(
     $                              iL,iNum,iGasID,iISO,dGasAmt,dPartitionFcn,
     $                              dLTE,dNLTE,dvibqft,
     $                              daLineShift,daELower,daStren296,
     $                              daStrenNLTE,daQtipsFactor)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iNum,iGasID,iISO,iL
      DOUBLE PRECISION dGasAmt,dPartitionFcn,dLTE,dNLTE,dvibqft
      DOUBLE PRECISION daLineShift(kHITRAN),daElower(kHITRAN)
      DOUBLE PRECISION daStren296(kHITRAN)     !!straight from HITRAN at 296 K
c output params
      DOUBLE PRECISION daStrenNLTE(kHITRAN),daQtipsFactor(kHITRAN)

c local variables
      INTEGER iLines
      DOUBLE PRECISION v0,s1,sb,se,s2,s3,TS,sq
      DOUBLE PRECISION qr,qrstd,qv,qvstd,qttem,qtstd,sq_lte,sq_nlte

      TS = 296.0 * 1.0d0

      sq_lte = dPartitionFcn    !!this is partition fcn at LTE         

cc Now figure out how SQ changes due to the fact the molecule is in NLTE
cc this was the OLD version of GENLN2; 
cc see /asl/packages/Genln2/Genln2/zpthadj.f
cc       IF ((NTEFLG) .AND.
cc     + ((BDAT) .OR. (TP .GT. 300.))) THEN
cc         IF (IGAS .LE. 32) THEN
cc           CALL QTFCT(IGAS,ISO,VTBOT,SQ)
cc         ELSEIF (IGAS .GE. 33) THEN
cc           CALL QTINT(IGAS,ISO,VTBOT,SQ)
cc       ENDIF
cc         IF (IGAS .EQ. 2 .OR. IGAS .EQ. 5) JFAC = 1.0
cc         IF (IGAS .EQ. 1 .OR. IGAS .EQ. 3) JFAC = 1.5
cc         SQ = SQ*(TP/VTBOT)**JFAC
cc so no longer need these lines
cc        !!!nonLTE partition fcn modifier
cc        IF ((iGasID .EQ. 2) .OR. (iGasID .EQ. 5)) jfac = 1.0d0
cc        IF ((iGasID .EQ. 1) .OR. (iGasID .EQ. 3)) jfac = 1.5d0
cc        daQtipsFactor(iLines) = ((dLTE/dNLTE) ** jfac)   

c this is from the new GENLN2 that Dave Edwards sent me
      CALL QNLTE(iGasID,iISO,dLTE,QR,QRSTD,QV,QVSTD)
      QTTEM = QR*DVIBQFT       !!qtotal at temp = T  = QR(T)  QV(T)
      QTSTD = QRSTD*QVSTD      !!qtotal at temp = T0 = QT(T0) QV(T0)
      SQ_NLTE = QTSTD/QTTEM

      DO iLines = 1,iNum
        v0 = daLineShift(iLines)

        !!!!boltzmann factor
        s1 = kPlanck2 *(dLTE-TS)/(dLTE*TS)
        sb = dexp(s1 * daELower(iLines))

        !!!stimulated emission
        s2 = -kPlanck2 * v0/dLTE
        s3 = -kPlanck2 * v0/TS

        IF (DEXP(s3) .EQ. 1.0) THEN
          se = 1.0d0
        ELSE
          se = (1.0d0 - dexp(s2))/(1.0d0 - dexp(s3))
        END IF

        sq = sq_nlte          !!this is partition fcn at NLTE         
        daStrenNLTE(iLines) = dGasAmt*daStren296(iLines)*sb*se*sq*kAvog

        !!!nonLTE partition fcn modifier
        IF (iGasID .EQ. 2) THEN
          daQtipsFactor(iLines) = SQ_NLTE/SQ_LTE * 1.0d0
        ELSE
          daQtipsFactor(iLines) = 1.0d0
        END IF
      END DO

c     write(*,1234) dLTE,dNLTE,dPartitionFcn,dvibqft
 1234 FORMAT(4(F15.5,' '))

      RETURN
      END

c************************************************************************
c this is to find how the upper populatios are enhanced/depleted
c             and how the lower populatios are enhanced/depleted
c this mimics pthadj.f in GENLN2
      SUBROUTINE NLTEPopulationRatios(
     $       iL,iNum,dP,dLTE,dLTEr1r2,dNLTE,dVibCenter,
     $       daElower,daLineShift,iaJ_UorL,
     $       daRTop,daRBot,daCFactor,daKFactor)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iNum,iL
      DOUBLE PRECISION dLTE,dLTEr1r2,dNLTE,dP
      DOUBLE PRECISION daLineShift(kHITRAN),daElower(kHITRAN),dVibCenter
      INTEGER iaJ_UorL(kHITRAN)

c output params
c daRTop tells how the upper populations are affected by NONLTE
c daRBot tells how the lower populations are affected by NONLTE
      DOUBLE PRECISION daRTop(kHITRAN),daRBot(kHITRAN)
c daKFactor tells you how the line strength is affected by the nonLTE
c daCFactor tells you how the planck fcn    is affected by the nonLTE
      DOUBLE PRECISION daCFactor(kHITRAN),daKFactor(kHITRAN)

c local variables
      INTEGER iLines,iTalk,iTdat,iBdat,iDefault
      DOUBLE PRECISION dTop,dBot,dEnergyTop,dEnergyBot,vtdef,dDiv,dDeltaT,dTK
c don't need this to be output at all, so leave it local
      DOUBLE PRECISION daDelta(kHITRAN)
      
      iDefault = -1       !!!!just summarize things        DEFAULT
      iTalk = +1          !!!!oh no!!! give out all the info!!!
      iTalk = -1          !!!!just summarize things        DEFAULT
      IF (iTalk .NE. iDefault) THEN
        print *,'iTalk,iDefault = ',iTalk,iDefault
      END IF

      !dLTE     is from the KLAYERS profiles (and/or averaged over mixed paths)
      !dLTEr1r2 is from the VibTemp profiles 

      vtdef = 2.75d+2
      dTK = dLTEr1r2

      DO iLines = 1, iNum
        !!set up the defaults to LTE
        dTop = dTK
        dBot = dTK

c        !genln4 stuff, see pthadj.f
c        IF (dP/kAtm2mb .LT. 1.0d-5) THEN      !!!this has to be in atmospheres
c          print *,'oopsy ... need default low P,Tvib : ',iL,dP,dP/kAtm2mb
c          dTop = min(dTK,vtdef)
c          dBot = min(dTK,vtdef)
c          CALL DoStop
c        END IF

        dEnergytop = daElower(iLines) + daLineShift(iLines)
        dEnergyBot = daElower(iLines)
        iTdat = -1
        iBdat = -1

        !!now put in the NLTE values
        !!from GENLN4, pthadj.f,ntepro.f   ... we are sending in v_upper IUSGQ
        !!KOPRA doc says dBot = dNLTE, dEnergyBot = dVibCenter IGNORE!!
        IF (iaJ_UorL(iLines) .EQ. +1) THEN
          !! we matched the iUSGQ
          dTop = dNLTE
          dBot = dBot
          dEnergytop = dVibCenter
          iTdat = +1
        ELSEIF (iaJ_UorL(iLines) .EQ. -1) THEN        
          !! we matched the iLSGQ
          dBot = dNLTE
          dTop = dTop
          dEnergyBot = dVibCenter
          iBdat = +1
        ELSE
          write(kStdErr,*) 'In NLTEPopulationRatios, had wierd iaJ_UorL'
          write(kStdErr,*) 'iLines,iaJ_UorL(iLines) = ',iLines,iaJ_UorL(iLines)
          CALL DoStop
        END IF

c !! genln4 stuff        
        dDeltaT = (dTK - dTop)/(dTK*dTop)
        daRTop(iLines) = dexp(-kPlanck2 * dEnergytop * dDeltaT)
        dDeltaT = (dTK - dBot)/(dTK*dBot)
        daRBot(iLines) = dexp(-kPlanck2 * dEnergyBot * dDeltaT)

        IF ((iBdat .EQ. +1) .AND. (iTdat .EQ. -1)) THEN
          daRTop(iLines) = daRBot(iLines)
        END IF

        daDelta(iLines) = dexp(-kPlanck2 * daLineShift(iLines)/dLTE)

        dDiv = 1.0d0/(1.0d0 - daDelta(iLines))
        IF (dabs(dTK-dNLTE) .GT. 1.0d-3) THEN
          daCFactor(iLines) = (daRBot(iLines)-daRTop(iLines))*dDiv
          daKFactor(iLines) =
     $                    (daRBot(iLines)-daRTop(iLines)*daDelta(iLines))*dDiv
        ELSE
          daCFactor(iLines) = 0.0d0
          daKFactor(iLines) = 1.0d0
        END IF
c        daANLTE(iLines)   = daKfactor(iLines)
c        daCNLTE(iLines)   = daRTop(iLines)

c        if (iLines .EQ. 41) THEN
c          write(*,1400) iL,iLines,daLineshift(iLines),dVibCenter,
c     $                  dLTE,dLTEr1r2,dNLTE,
c     $                  daRTop(iLines),daKFactor(iLines)
c        END IF

        IF (iTalk .GT. 0) THEN
          write(*,2000) iLines,daRTop(iLines),daRBot(iLines),daDelta(iLines),
     $                1.0d0-daCFactor(iLines)/daKFactor(iLines),
     $                daKFactor(iLines)

          write(*,2000) iLines,daRTop(iLines),daRBot(iLines),daDelta(iLines),
     $                daCFactor(iLines),daKFactor(iLines)
        END IF
      END DO

 1400 FORMAT(2(I3,' '),2(F11.5,' '),3(F7.3,' '),2(F10.5,' '))
 2000 FORMAT(I3,' ',D12.4,' ',D12.4,' ',D12.4,' ',D12.4,' ',D12.4,' ',D12.4)
 
      RETURN
      END

c************************************************************************
c simple subroutine to assign dAlpha,dBeta
c already taken care of daQtipsFactor in subroutine Strengths_nlte_vibpartfcn
      !ALPHA tells how NLTE optical depths modified from LTE optical depths
      !BETA  tells how planck function     modified from LTE Planck fcn
      !            though this is literally one half of the story

      !dAlpha  = daKFactor(iLines)
      !dPlanck = 1-daCFactor(iLines)/daKfactor(iLines)
      !        = daRTop(iLines)/daKfactor(iLines)
      SUBROUTINE AlphaBeta_Factors(
     $                    iLines,iUpdateSUMNLTE,rNLTEstrength,daRTop,daRBot,
     $                    daCFactor,daKFactor,dAlpha,dBeta)      

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iLines          !!!which line we are looking at
      INTEGER iUpdateSUMNLTE  !!!is abs(dNLTE-dLTE) > 0.1
      DOUBLE PRECISION daRTop(kHITRAN),daRBot(kHITRAN)
      DOUBLE PRECISION daCFactor(kHITRAN),daKFactor(kHITRAN)
      REAL rNLTEstrength      !!how much to weight the NLTE by
c output params 
      DOUBLE PRECISION dAlpha,dBeta

c local vars
      DOUBLE PRECISION dNLTEStrength
      INTEGER iAlwaysLTE

      dNLTEStrength = rNLTEStrength * 1.0d0

      IF (iUpdateSumNLTE .EQ. +1) THEN
        dAlpha = dNLTEstrength*daKFactor(iLines)  !!!directly from GENLN2
        dBeta  = dNLTEstrength*daRTop(iLines)     !!!directly from GENLN2
      ELSE
        dBeta  = 1.0d0 !!!if in LTE, rBot = rTop = 1 ==> planck unmodified
        dAlpha = 1.0d0
      END IF

c      !!!! to test this mode, by setting it to LTE ratios
c      iAlwaysLTE = -1
c      IF (iAlwaysLTE .EQ. +1) THEN
c        !this is the equivalent of iAllLayersLTE = +1 in the namelist file
c        dBeta  = 1.00d0
c        dAlpha = 1.00d0
c      END IF

      RETURN
      END

c************************************************************************
c this subroutine reads in the weak LTE optical depths that have been computed
c at the US STandard Atmosphere plus 10 offsets, and dumps it out to the 
c pressure layering specified, at the gas amounts specified
c this is basically a cut and dried version of subroutine "othergases"
c in file kcoeffMAIN.f
      SUBROUTINE read_std_optdepths_new(iGasID,raFreq,iStart,
     $            raTAmt,raTTemp,raTPress,raTPartPress,
     $            raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, 
     $            iProfileLayers,iSplineType,pProf,
     $            raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,
     $            iUnCompressType,
     $            iLowerOrUpper,daaWeakOptDepth)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters self explanatory .. the arrays contain the current profile
      INTEGER iUnCompressType !! -2,-3 for weak backgnd CO2; see CompFileName
      INTEGER iGasID
      INTEGER iStart        !!! layer above which there is NLTE
      REAL raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
      REAL raTPartPress(kProfLayer),raFreq(kMaxPts)
      INTEGER iLowerOrUpper !!-1 if usual KLAYERS studff, +1 if UA
c these are the individual reference profiles
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c this sets the vertical temp profile, to be checked for all gases
      REAL raVertTemp(kProfLayer),rFileStartFr
      INTEGER iVertTempSet,iTag,iActualTag
      INTEGER iL_low,iL_high,iSPlineType
      INTEGER iProfileLayers   !number of layers in RTP file
      REAL pProf(kProfLayer)
c output parameter
c daaWeakOptDepth tells LTE optical depth
      DOUBLE PRECISION daaWeakOptDepth(kMaxPts,kProfLayer)

c local variables
      DOUBLE PRECISION daaDQ(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPts,kProfLayer)
      CHARACTER*120 caFName 
      INTEGER iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn 
      INTEGER iT0,iaTsort(kMaxTemp),iErr
      DOUBLE PRECISION dSfreq,dFStep,daToffset(kMaxTemp) 
      DOUBLE PRECISION daaaKX(kMaxK,kMaxTemp,kMaxLayer) 
      DOUBLE PRECISION daaUX(kMaxPts,kMaxK) 
      INTEGER iaChiChunks(kMaxGas),iChiChunks,iDoFudge,WhichGasPosn 
      INTEGER iLLS,iType,iDefault,iL
        
      IF (iUnCompressType .NE. -2 .AND. iUnCompressType .NE. -3) THEN
        write(kStdErr,*) 'Code is trying to uncompress weakbackgnd CO2 depths'
        write(kStdErr,*) 'either in the usual AIRS layers (1100-0.005 mb) or'
        write(kStdErr,*) 'in the upper AIRS layers (0.005 - 0.000025 mb)'
        write(kStdErr,*) 'need iUnCompressType -2,-3, not ',iUnCompressType
        CALL DoStop
      END IF

      iIOUN = kCompUnit 
      CALL CompFileName(iUnCompressType,iGasID,rFileStartFr,iTag,iActualTag,
     $                  caFName) 
      CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay, 
     $              iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, 
     $              daaaKX,daaUX) 

c check that the file has the data for the correct gas 
      IF (iFileGasID .NE. iGasID) THEN 
        iErr=1 
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID 
 1000   FORMAT('Error! file : ',/,A80,/, 
     $         'contains data for GasID ',I3,' not desired GasID : ',I3) 
        CALL DoSTOP 
      END IF 
 
c check that the data file has the right number of layers 
      IF (iNLay .NE. kMaxLayer) THEN 
        iErr=1 
        WRITE(kStdErr,1010) caFName,iNLay,kMaxLayer 
 1010   FORMAT('Error! file : ',/,A80,/, 
     $         'contains data for ',i3,' layers but kMaxLayer = ',I3) 
        CALL DoSTOP 
      END IF 
 
      CALL GetAbsCoeffNOJAC(daaWeakOptDepth,daToffset,daaaKx,daaUx, 
     $         raTTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, 
     $         pProf,iProfileLayers,iSplineType,iLowerOrUpper) 

c convert absorption coefficient correctly if necessary 
      IF (iKtype .eq. 2) THEN 
        CALL RaisePower(daaWeakOptDepth) 
      END IF 
 
c now compute optical depth = gas amount * abs coeff 
      CALL AmtScale(daaWeakOptDepth,raTAmt) 

c      do iL = 1,kProfLayer
c        print *,'weak',iL,daaWeakOptDepth(1,iL),raTAmt(iL)
c        end do

      RETURN
      END

c************************************************************************
c this subroutine almost same as read_std_optdepths_new above. BUT : it
c uses the new database generated for 70-120 km 
      SUBROUTINE read_std_optdepths_upper_UA(
     $                iGas,iSplineType,raUpperPressLevels,
     $                pProfNLTE_upatm,iNumProfNLTE_upatm,
     $                rFileStartFr,iTag,iActualTag,
     $                iGasID,raFreq,
     $                raTAmt,raTTemp,raTPress,raTPartPress,
     $                iNumUpperLayers,daaWeakOptDepth)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 
      include '../INCLUDE/airslevels_upper.param'
      include '../INCLUDE/airslevelheights_upper.param'

c input parameters self explanatory .. the arrays contain the current profile
      INTEGER iTag,iActualTag,iGasID,iGas,iSplineType
      REAL raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
      REAL raTPartPress(kProfLayer),raFreq(kMaxPts),rFileStartFr
      INTEGER iNumUpperLayers    !!!! how many raTXYZlayers are there?
c note these two correspond to each other
      REAL raUpperPressLevels(kProfLayer+1),pProfNLTE_upatm(kProfLayer)
      INTEGER iNumProfNLTE_upatm
c output parameter
c daaWeakOptDepth   tells LTE optical depth
      DOUBLE PRECISION daaWeakOptDepth(kMaxPts,kProfLayer)

c local vars
      REAL pProf(kProfLayer)       !!! this is the layers profile
      REAL raVertTemp(kProfLayer)  !!! this is the temperature profile
      INTEGER iStart,iI,iUnCompressType,iFr
      REAL raXTemp(kProfLayer),raXAmt(kProfLayer),raXPress(kProfLayer)
      REAL raXPartPress(kProfLayer)
c these are the individual reference profiles
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      INTEGER iL_low,iL_high,iVertTempSet,iLowerOrUpper

      REAL raP2(kMaxLayer),raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)

      IF (iGasID .NE. 2) THEN
        write(kStdErr,*) 'read_std_optdepths_upper_UA only for CO2'
        write(kStdErr,*) 'you have called it with gasID = ',iGasID
        CALL DoStop
      END IF

c units : require raTPress,raTPartPress in atm !!!!!!
c                 raRPress,raXPress     in atm
c                 pProf in mb

c      print *,'>>>>>>>>>>>>> FANCULA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
c      CALL ua_avg_databasepressures(raP2,raPP2,raT2,raA2)
c      iNumUpperLayers = 100
c      DO iI = 1,kProfLayer
c        raTPress(iI) = raP2(iI)/kAtm2mb
c        raTPartPress(iI) = raPP2(iI)/kAtm2mb
c        raTTemp(iI) = raT2(iI)
c        raTAmt(iI) = raA2(iI)
c        pProfNLTE_upatm(iI) = raTPress(iI)*kAtm2mb
c        raUpperPressLevels(iI) = DATABASELEV(iI)
c      END DO
c      iI = 101
c      raUpperPressLevels(iI) = DATABASELEV(iI)
c      DO iI = 1,iNumUpperLayers
c        print *,iI,raTPress(iI),raTPartPress(iI),raTTemp(iI),raTAmt(iI),
c     $          raUpperPressLevels(iI)
c      END DO
c      print *,'>>>>>>>>>>>>> FANCULA >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

c      !need to put stuff into (iLowest --> 100) instead of (1 .. iX)
c      print *,iNumUpperLayers
c      print *,(pProfNLTE_upatm(iI),iI=1,100)
c      print *,(raUpperPressLevels(iI),iI=1,100)
c      print *,(raTPress(iI),iI=1,100)
c      call dostop

      iStart = (kProfLayer - iNumUpperLayers + 1)
      IF ((iStart .LT. 1) .OR. (iStart .GT. kProfLayer)) THEN
        write(kStdErr,*) 'in read_std_optdepths_upper_UA '
        write(kStdErr,*) '  called with iNumUpperLayers = ',iNumUpperLayers
        write(kStdErr,*) '  which gives (bad) iStart = ',iStart
        CALL DoStop
      END IF
      DO iI = 1,iStart-1
        raXTemp(iI)      = 300.0
        raXAmt(iI)       = 0.0
        raXPress(iI)     = 1013.0/kAtm2mb
        raXPartPress(iI) = 0.0
        !! need to do these two
        !! pProfNLTE_upatm may have DIFFERENT numlayers than iNumUpperLayers???
        !! fancula fancula !!
        pProf(iI)        = 1013.0
        raVertTemp(iI)   = 300.0
      END DO
      DO iI = iStart,kProfLayer
        raXTemp(iI)      = raTTemp(iI-iStart+1)
        raXAmt(iI)       = raTAmt(iI-iStart+1)
        raXPress(iI)     = raTPress(iI-iStart+1)
        raXPartPress(iI) = raTPartPress(iI-iStart+1)
        !! need to do these two
        !! pProfNLTE_upatm may have DIFFERENT numlayers than iNumUpperLayers???
        !! fancula fancula !!
        pProf(iI)        = pProfNLTE_upatm(iI-iStart+1)
        raVertTemp(iI)   = raTTemp(iI-iStart+1)
      END DO

      CALL upper_co2_default_profile(iGas,iGasID,pProf,iNumUpperLayers,
     $             raUpperPressLevels,iSplineType,
     $             raRAmt,raRTemp,raRPress,raRPartPress)

      iL_low  = 1
      iL_high = kProfLayer
      iVertTempSet = +1
      iUnCompressType = -3  !!!upper atm weak CO2 backgnd
      iLowerOrUpper = +1
      CALL read_std_optdepths_new(iGasID,raFreq,iStart,
     $              raXAmt,raXTemp,raXPress,raXPartPress,
     $              raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, 
     $              iNumUpperLayers,iSplineType,pProf,
     $              raVertTemp,iVertTempSet,rFileStartFr,
     $              iTag,iActualTag,iUnCompressType,
     $              iLowerOrUpper,daaWeakOptDepth)

      !!now shift the opt depths so that layers 1 .. iUpper have kopticaldepth
      DO iI = iStart,kProfLayer
        DO iFr = 1,kMaxPts
          daaWeakOptDepth(iFr,iI-iStart+1) =  daaWeakOptDepth(iFr,iI)
        END DO
      END DO
       
      !! and the rest are zeroed
      iL_low = kProfLayer + 1
      iL_low = iL_low - iStart + 1
      DO iI = iL_Low,kProfLayer
        DO iFr = 1,kMaxPts
          daaWeakOptDepth(iFr,iI) =  0.0d0
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the upper atm default profile
c and puts it into kProfLayers
      SUBROUTINE upper_co2_default_profile(iGas,iGasID,pProf,iNumLayers,
     $                                  raUpperPressLevels,iSplineType,
     $                                  raRAmt,raRTemp,raRPress,raRPartPress)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input
      REAL pProf(kProfLayer)       !!! this is the layers profile
      REAL raUpperPressLevels(kProfLayer+1)
      INTEGER iGasID,iGas,iNumLayers,iSplineType
c output
c these are the individual reference profiles
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)

c local
      CHARACTER*80 caFName
      INTEGER iI,iLowerOrUpper,iError
      REAL raDummyThickness(kProfLayer)

c these are the individual reference profiles, at kMaxLayer layers 
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer) 
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer) 
c these are the arbitrary profiles (avg press) stored in matrices 
      REAL raaPress(kProfLayer,kGasStore) 

      DO iI = 1,kProfLayer
        raaPress(iI,iGas) = pProf(iI)/kAtm2mb
        raRamt(iI)        = 0.0
        raRTemp(iI)       = 0.0
        raRPress(iI)      = 0.0
        raRPartPress(iI)  = 0.0
      END DO

      iLowerOrUpper = +1
      CALL FindReferenceName(caFName,iGasID,+1)
      CALL ReadRefProf(caFName,kMaxLayer,
     $         raR100Amt,raR100Temp,raR100Press,raR100PartPress,iError) 
      CALL MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $        raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $        raaPress,iGas,iGasID,iNumLayers,
     $        raUpperPressLevels,raDummyThickness,iSplineType,+1,iError)

      RETURN
      END

c************************************************************************
c             these subroutines read in HITRAN parameters
c************************************************************************
c this subroutine reads in the line params for the strong bands that are in LTE
c it reads in all lines, then compares against the list of NLTE lines and 
c removes them, leaving only the stronger LTE lines.
      SUBROUTINE read_stronglineLTE_lineparameters(
     $   iGasID,caaaNLTEBands,iaNLTEBands,iLTEin,dLineStrenMin,caStrong,
     $   iNum,daIso,daElower,daLineCenter,daJL,daJU,daPshift,
     $   daStren296,daW_For,daW_self,daW_temp)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params
c caStrong tells us the name of the file that has the weak bkgnd line params
c dLineStrneMin = -1 if use all lines, some other number if use lines  >= X
c iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
c caaaNLTEBands tells the name of the files containing the line parameters
c iLTEin    tells which set of data to use
      CHARACTER*80 caStrong
      DOUBLE PRECISION dLineStrenMin
      INTEGER iaNLTEBands(kGasStore),iLTEin
      CHARACTER*80 caaaNLTEBands(kGasStore,kNumkCompT) 
c output params 
      INTEGER iGasID,iNum
      DOUBLE PRECISION daElower(kHITRAN),daLineCenter(kHITRAN),daISO(kHITRAN)
      DOUBLE PRECISION daJL(kHITRAN),daJU(kHITRAN)
      DOUBLE PRECISION daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
      DOUBLE PRECISION daW_self(kHITRAN),daW_temp(kHITRAN)

c local variables
      DOUBLE PRECISION daElower0(kHITRAN),daLineCenter0(kHITRAN)
      DOUBLE PRECISION daJL0(kHITRAN),daJU0(kHITRAN),daISO0(kHITRAN)
      DOUBLE PRECISION daPshift0(kHITRAN),daStren2960(kHITRAN)
      DOUBLE PRECISION daW_self0(kHITRAN),daW_for0(kHITRAN),daW_temp0(kHITRAN)

      DOUBLE PRECISION daElowerSum(kHITRAN),daLineCenterSum(kHITRAN)
      DOUBLE PRECISION daJLSum(kHITRAN),daJUSum(kHITRAN),daISOSum(kHITRAN)
      DOUBLE PRECISION daPshiftSum(kHITRAN),daStren296Sum(kHITRAN)
      DOUBLE PRECISION daW_selfSum(kHITRAN),daW_forSum(kHITRAN)
      DOUBLE PRECISION daW_tempSum(kHITRAN)

      DOUBLE PRECISION dMin,dMax

      INTEGER iI,iJ,iK,iErr,iIOUN,iNum0,iBands,iLoop,iCumLineSum,iISO
      INTEGER iFound,iRead
      CHARACTER*80 caFname,caTemp
      CHARACTER*80 caaCheckNamesIn(kGasStore)
      CHARACTER*80 caaCheckNamesMatch(kGasStore)
      CHARACTER*80 caaCheckNamesUse(kGasStore)
      DOUBLE PRECISION dLmin,dLmax
      INTEGER iaBin(10)   !!!to bin the linestrengths
      INTEGER iaOk(kHITRAN),iStrongLTEFound,iStrongNLTEFound

      iStrongLTEFound  = 0   !!!how many strong bands are in LTE
      iStrongNLTEFound = 0   !!!how many strong bands are in NLTE

c first read the names of all strong bands for the molecule
      caFName = caStrong
      iIOun = kTempUnit
      OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',
     $     IOSTAT=iErr) 
      IF (iErr .NE. 0) THEN 
        write (kStdErr,*) 'read_stronglineLTE_lineparameters, error reading'
        write (kStdErr,*) 'datafile that contains names of strong bands'
        WRITE(kStdErr,1070) iErr, caFName 
        CALL DoSTOP 
      END IF 
      kTempUnitOpen = 1 

      iBands = 0
      iRead = 0
 10   CONTINUE
      read(iIOUN,15,ERR = 20, END = 30) caTemp
      IF (caTemp(1:1) .EQ. '!') THEN
        !! this is a comment line, ignore
        GOTO 10
      END IF
      iRead = iRead + 1
      caaCheckNamesIn(iRead) = caTemp
      DO iI = 1,iaNLTEBands(iLTEin)
        IF (caaaNLTEBands(iLTEin,iI) .EQ. caTemp) THEN
          iBands = iBands + 1
          caaCheckNamesMatch(iBands) = caTemp
        END IF
      END DO
      GOTO 10
      
      kTempUnitOpen = -1 
 15   FORMAT(A80)

 20   write (kStdErr,15) 'Error reading info in ',caFName

 30   CONTINUE

      IF (iBands .EQ. 0) THEN
        write(kStdErr,*) ' '
        write(kStdErr,*) 'In subroutine read_stronglineLTE_lineparameters :' 
        write(kStdErr,*) '  your namelist file says "caaStrongLines" contains'
        write(kStdErr,*) '  the list of strong bands that are in LTE or NLTE'
        write(kStdErr,*) 'A typical datafile name in "caaStrongLines" is (from caStrong)'
        write(kStdErr,*) caaCheckNamesIn(1)
        write(kStdErr,*) '  your kcarta.param says "caStrongLineParams" has'
        write(kStdErr,*) '  typical name of parameter file (has eg (line center, width,stren etc)'
        write(kStdErr,*) 'Bandfilename to be found using NLTEBandMapper,caaaNLTEBands(iI,iJ) is : '
        write(kStdErr,*) caaaNLTEBands(iLTEin,1)
        write(kStdErr,*) 'Cannot find a single match!!! '
        write(kStdErr,*)    'in namelist file check caaStrongLines = '
        write(kStdErr,*) caStrong
        write(kStdErr,*)    'in kcarta.param check caStrongLineParams = '
        write(kStdErr,*) caStrongLineParams
        CALL DoStop
      END IF

c find (band) names in caStrong that are NOT in NLTE list caaaNLTEBands!!!!
c thus find names in caaCheckNamesIn that are NOT in caaCheckNamesMatch
      iK = 0
      DO iI = 1,iRead
        iFound = -1
        DO iJ = 1,iBands
          IF (caaCheckNamesIn(iI) .EQ. caaCheckNamesMatch(iJ)) THEN
            iFound = +1        !!!name from caStrong is in the NLTE list
            iStrongNLTEFound  = iStrongNLTEFound + 1
          ENDIF
        END DO
        IF (iFound .EQ. -1) THEN
          !!!name from caStrong is NOT in the NLTE list
          iK = iK + 1
          caaCheckNamesUse(iK) = caaCheckNamesIn(iI)
          iStrongLTEFound  = iStrongLTEFound + 1
        END IF
      END DO

      write(kStdWarn,*) 'in subroutine read_stronglineLTE_lineparameters :'
      write(kStdWarn,*) 'found ',iStrongLTEFound ,' strong bands in  LTE'
      write(kStdWarn,*) 'found ',iStrongNLTEFound,' strong bands in NLTE'

      IF (iStrongNLTEFound .EQ. 0) THEN
        write(kStdErr,*) 'in subroutine read_stronglineLTE_lineparameters :'
        write(kStdErr,*) 'did not find ANY Strong NLTE bands!!!'
        write(kStdErr,*) 'code thinks all your strong bands are in LTE'
        write(kStdErr,*) 'Check info in caaStrongLines (from nm_nonlte) '
        write(kStdErr,*) '      versus that in NLTEBandMapper subroutine'
        CALL DoStop
      END IF

      IF (iStrongNLTEFound .NE. iaNLTEBands(iLTEin)) THEN
        write(kStdErr,*) 'in subroutine read_stronglineLTE_lineparameters :'
        write(kStdErr,*) 'number of Strong NLTE bands <> iaNLTEBands(iGasID)!'
        write(kStdErr,*) 'Check info in caaStrongLines (from nm_nonlte) '
        write(kStdErr,*) '      versus that in NLTEBandMapper subroutine'
        write(kStdErr,*) '      versus that in caaaNLTEBands (from nm_nonlte)'
        write(kStdErr,*) 'iStrongNLTEFound = ',iStrongNLTEFound
        write(kStdErr,*) 'iaNLTEBands(iLTEin) = ',iaNLTEBands(iLTEin),iLTEin
        CALL DoStop
      END IF

      write(kStdWarn,*) 'putting together the Strong LTE bands .....'

      iBands = iK

c now read the line parameters of the files
      iCumLineSum = 0

      DO iLoop = 1,iBands
        caFName = caaCheckNamesUse(iLoop)
        iIOun = kTempUnit

        OPEN(UNIT=iIOun,FILE=caFName,FORM='unformatted',STATUS='OLD',
     $     IOSTAT=iErr) 
        IF (iErr .NE. 0) THEN 
          write (kStdErr,*) 'read_stronglineLTE_lineparameters, error reading'
          write (kStdErr,*) 'datafile that contains parameters of strong band'
          WRITE(kStdErr,1070) iErr, caFName 
          CALL DoSTOP 
        END IF 
        kTempUnitOpen = 1 

        read(iIOUN) iGasID,iNum0,iISO
        write(kstdWarn,*) ' '
        write(kStdWarn,*) '----> Opening HITRAN parameter file : '
        write(kStdWarn,1080) caFName
        write(kStdWarn,*) '  with ',iNum0,' line params gID,iso = ',iGasID,iISO

        IF ((iNum0 + iCumLineSum) .GT. kHITRAN) THEN
          write(kStdErr,*) 'Total lines so far = ',iNum0,' for gas ',iGasID
          write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
          write(kStdErr,*) 'Please check kHITRAN in kcarta.param and fix'
          CALL DoStop
        END IF

        read(iIOUN) (daElower0(iI),iI = 1,iNum0)     !lower state energy
        read(iIOUN) (daLineCenter0(iI),iI = 1,iNum0) !line center    
        read(iIOUN) (daJL0(iI),iI = 1,iNum0)         !lower quantum number
        read(iIOUN) (daJU0(iI),iI = 1,iNum0)         !upper quantum number
        read(iIOUN) (daPshift0(iI),iI = 1,iNum0)     !pressure shift linecenter
        read(iIOUN) (daStren2960(iI),iI = 1,iNum0)   !line strenght
        read(iIOUN) (daW_for0(iI),iI = 1,iNum0)      !foreign broadening/atm
        read(iIOUN) (daW_self0(iI),iI = 1,iNum0)     !self broadening/atm
        read(iIOUN) (daW_temp0(iI),iI = 1,iNum0)     !broaden tempr dependance

        DO iI = 1,iNum0
          daIso0(iI)                      = iISO*1.0d0
          daIsoSum(iCumLineSum+iI)        = daIso0(iI)
          daELowerSum(iCumLineSum+iI)     = daELower0(iI)
          daLineCenterSum(iCumLineSum+iI) = daLineCenter0(iI)
          daJLSum(iCumLineSum+iI)         = daJL0(iI)
          daJUSum(iCumLineSum+iI)         = daJU0(iI)
          daPShiftSum(iCumLineSum+iI)     = daPShift0(iI)
          daStren296Sum(iCumLineSum+iI)   = daStren2960(iI)
          daW_forSum(iCumLineSum+iI)      = daW_for0(iI)
          daW_selfSum(iCumLineSum+iI)     = daW_self0(iI)
          daW_tempSum(iCumLineSum+iI)     = daW_temp0(iI)
        END DO

        iCumLineSUm = iCumLineSum + iNum0

        close (iIOUN)
        kTempUnitOpen=-1 
      END DO        

 1070 FORMAT('ERROR! number ',I5,' opening HITRAN parameter file:',/,A80) 
 1080 FORMAT(A80) 

      !bin the linestrengths
      DO iI = 1,10
        iaBin(iI) = 0
      END DO
      dLMin = +1.0d30
      dLMax = -1.0d30
            
      iNum = 0
      DO iI = 1,iCumLineSum
        IF (daSTren296Sum(iI) .GT. dLmax) dLMax = daSTren296Sum(iI)
        IF (daSTren296Sum(iI) .LT. dLmin) dLMin = daSTren296Sum(iI)
        IF (log10(daSTren296Sum(iI)) .GT. -20.0) THEN
            iaBin(1) = iaBin(1) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -21.0) THEN
          iaBin(2) = iaBin(2) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -22.0) THEN 
          iaBin(3) = iaBin(3) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -23.0) THEN
          iaBin(4) = iaBin(4) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -24.0) THEN
          iaBin(5) = iaBin(5) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -25.0) THEN 
          iaBin(6) = iaBin(6) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -26.0) THEN 
          iaBin(7) = iaBin(7) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -27.0) THEN 
          iaBin(8) = iaBin(8) + 1
        ELSEIF (log10(daSTren296Sum(iI)) .GT. -28.0) THEN 
           iaBin(9) = iaBin(9) + 1
        ELSE 
          iaBin(10) = iaBin(10) + 1
        END IF

        !!!if no need to prune, set all lines as output
        IF (dLineStrenMin .LT. 0.0d0) THEN    !!!use all lines
          iNum = iCumLineSum
          daIso(iI)        = daIsoSum(iI)        !isotope
          daElower(iI)     = daElowerSum(iI)     !lower state energy
          daLineCenter(iI) = daLineCenterSum(iI) !line center    
          daJL(iI)         = daJLSum(iI)         !lower quantum number
          daJU(iI)         = daJUSum(iI)         !upper quantum number
          daPshift(iI)     = daPshiftSum(iI)     !pressure shift linecenter
          daStren296(iI)   = daStren296Sum(iI)   !line strenght
          daW_for(iI)      = daW_forSum(iI)      !foreign broadening/atm
          daW_self(iI)     = daW_selfSum(iI)     !self broadening/atm
          daW_temp(iI)     = daW_tempSum(iI)     !broadening tempr dependance
        ELSEIF (daStren296Sum(iI) .LT. dLineStrenMin) THEN
          iaOK(iI) = -1          !!!do not allow this line to be used
        ELSEIF (daStren296Sum(iI) .GE. dLineStrenMin) THEN
          iaOK(iI) = +1          !!!allow this line to be used
          iNum = iNum + 1
          daIso(iNum)        = daIsoSum(iI)        !isotope
          daElower(iNum)     = daElowerSum(iI)     !lower state energy
          daLineCenter(iNum) = daLineCenterSum(iI) !line center    
          daJL(iNum)         = daJLSum(iI)         !lower quantum number
          daJU(iNum)         = daJUSum(iI)         !upper quantum number
          daPshift(iNum)     = daPshiftSum(iI)     !pressure shift linecenter
          daStren296(iNum)   = daStren296Sum(iI)   !line strength
          daW_for(iNum)      = daW_forSum(iI)      !foreign broadening/atm
          daW_self(iNum)     = daW_selfSum(iI)     !self broadening/atm
          daW_temp(iNum)     = daW_tempSum(iI)     !broadening tempr dependance
        END IF
      END DO

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Backgnd lines : '
      write(kStdWarn,*) 'min, max linestrength = ',dLmin,dLmax
      write(kStdWarn,*) 'number of lines to be used = ',iNum
      IF (iNum .EQ. 0) THEN
        write(kStdWarn,*) 'ALL Strong Band lines are in NLTE !!!!'
        write(kStdWarn,*) 'No  Strong Band lines are in  LTE !!!!'
      ELSE
        write(kStdWarn,*) 'Some Strong Band lines are in LTE !!!!'
        write(kStdWarn,*) 'log10(strength)     num of lines'
        write(kStdWarn,*) '--------------------------------'
        DO iI = 1,10
          write(kStdWarn,300) -19-iI,iaBin(iI)
        END DO
      END IF
      write(kStdWarn,*) ' '

 300  FORMAT('        ',I4,'             ',I5)

      RETURN
      END

c************************************************************************
c this subroutine reads in the line parameters for the band in question 
      SUBROUTINE read_lineparameters(iLTEin,iBand,caaaNLTEBands,
     $     iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift,
     $     daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,
     $     iLineMixBand,iDoVoigtChi)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params, read from caaNLTETemp(iLTEin)
c caaaNONLTETemp  tells the name of the files containing the nonLTE temps      
      INTEGER iLTEIn,iBand,iDOVoigtChi
      CHARACTER*80 caaaNLTEBands(kGasStore,kNumkCompT) 
c output params 
      INTEGER iGasID,iNum,iISO,iLineMixBand
      DOUBLE PRECISION daElower(kHITRAN),daLineCenter(kHITRAN)
      DOUBLE PRECISION daJL(kHITRAN),daJU(kHITRAN)
      DOUBLE PRECISION daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
      DOUBLE PRECISION daW_self(kHITRAN),daW_temp(kHITRAN)
      DOUBLE PRECISION daJLowerQuantumRot(kHITRAN)
      CHARACTER*1      caJPQR(kHITRAN)

      INTEGER iI,iErr,iIOUN,iJU,iJL
      CHARACTER*80 caFname
      DOUBLE PRECISION dMax,dL,dR,dC

      caFName = caaaNLTEBands(iLTEin,iBand)

      iIOun = kTempUnit
      OPEN(UNIT=iIOun,FILE=caFName,FORM='unformatted',STATUS='OLD',
     $     IOSTAT=iErr) 
      IF (iErr .NE. 0) THEN 
        write (kStdErr,*) 'in subroutine read_lineparameters, error reading'
        write (kStdErr,*) 'file that has HITRAN lineparameters'
        WRITE(kStdErr,1070) iErr, caFName 
        CALL DoSTOP 
      END IF 
      kTempUnitOpen=1 

      read(iIOUN) iGasID,iNum,iISO
      write(kstdWarn,*) ' '
      write(kStdWarn,*) '----> Opening HITRAN parameter file : '
      write(kStdWarn,1080) caFName
      write(kStdWarn,*) 'File has ',iNum,' line parameters for gas ',iGasID

      IF (iNum .GT. kHITRAN) THEN
        write(kStdErr,*) 'File has ',iNum,' line parameters for gas ',iGasID
        write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
        write(kStdErr,*) 'Please check kHITRAN in kcarta.param and fix'
        CALL DoStop
      END IF

      read(iIOUN) (daElower(iI),iI = 1,iNum)     !lower state energy
      read(iIOUN) (daLineCenter(iI),iI = 1,iNum) !line center    
      read(iIOUN) (daJL(iI),iI = 1,iNum)         !lower vib quantum number .. basically SAME for all iI=1,iNum
      read(iIOUN) (daJU(iI),iI = 1,iNum)         !upper vib quantum number .. basically SAME for all iI=1,iNum
      read(iIOUN) (daPshift(iI),iI = 1,iNum)     !pressure shift of linecenter
      read(iIOUN) (daStren296(iI),iI = 1,iNum)   !line strenght
      read(iIOUN) (daW_for(iI),iI = 1,iNum)      !foreign broadening/atm
      read(iIOUN) (daW_self(iI),iI = 1,iNum)     !self broadening/atm
      read(iIOUN) (daW_temp(iI),iI = 1,iNum)     !broadening tempr dependance
      !! new since July 2015, comes from line.bslq (see lineparameters.m in
      !! SRCv1.18/NONLTE/M_Files_for_kcarta_NLTE_LBL_runs/USUALLAYERS/
      read(iIOUN) (daJLowerQuantumRot(iI),iI = 1,iNum)  !J lower quantum rotation state
c      read(iIOUN) (caJPQR(iI),iI = 1,iNum)       !P Q or R
      
      close (iIOUN)
      kTempUnitOpen=-1 

      iJU = nint(daJU(1))
      iJL = nint(daJL(1))

 1070 FORMAT('ERROR! number ',I5,' opening HITRAN parameter file:',/,A80) 
 1080 FORMAT(A80) 

c outside of the weaklines, do linemixing for all the bands the user specifies!
c but be careful about Cousin vs linemix, else code becomes VERY slow
      IF (iGasID .EQ. 2) THEN
        IF ((iJL .EQ. 4) .AND. (iJU .EQ. 24) .AND. (iISO .EQ. 1)) THEN
          iLineMixBand = +2
          write(kStdWarn,*) 'strong CO2 band : deldel (2310) iLineMixBand = +2'
        ELSEIF ((iJL .EQ. 2).AND.(iJU .EQ.16) .AND. (iISO .EQ. 1)) THEN
          iLineMixBand = +2
          write(kStdWarn,*) 'strong CO2 band : pi pi  (2320) iLineMixBand = +2'
        ELSEIF ((iJL .EQ. 1).AND. (iJU .EQ.9) .AND. (iISO .EQ. 1)) THEN
          iLineMixBand = +2
          write(kStdWarn,*) 'strong CO2 band : sigsig (2350) iLineMixBand = +2'
        ELSEIF ((iJL .EQ. 1).AND. (iJU .EQ.9) .AND. (iISO .EQ. 2)) THEN
          iLineMixBand = +2
          write(kStdWarn,*) 'strong CO2 band : sigsig (2351) iLineMixBand = +2'
        ELSE
          write(kStdWarn,*) 'strongest CO2 bands : iLineMixBand = +1 (cousin)'
          iLineMixBand = +1
        END IF
      ELSE
        iLineMixBand = -1
        write(kStdWarn,*) 'medium strength CO2 bands : no linemixing (voigt)'
      END IF

      IF ((iLineMixBand .EQ. 2) .AND. (iJU .NE. 9)) THEN
        !!!only do linemix for 2350, 2351 lines; else do Cousin
        iLineMixBand = +1
        write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
        write(kStdWarn,*) '   reset iLineMixBand = from +2 to +1 (cousin)'
      ELSEIF ((iLineMixBand .EQ. 2) .AND. (iJU .EQ. 9)) THEN
        IF (iISO .LE. 2) THEN
          write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
          write(kStdWarn,*) '   iLineMixBand = +2 (linemix)'
        ELSEIF (iISO .GT. 2) THEN
          write(kStdWarn,*) 'iJL, iJU,iSO = ',iJL,iJU,iISO
          write(kStdWarn,*) '   iLineMixBand = +1 (cousin)'
        END IF
      END IF
        
c this always makes cousin the happening one!
c       IF (iLineMixBand .GT. 0) THEN
c         iLineMixBand = +1
c         write(kStdWarn,*) 'doing cousin everywhere; reset iLineMixBand = +1'
c         print *,' ***** reset iLineMixBand = 1   ie linemix --> cousin ****' 
c       END IF

c      IF (iLineMixBand .GT. 0) THEN
c        IF (iDoVoigtChi .GT. 0) THEN
c          write(kStdErr,*) 'Cannot have idoVoigtChi > 0 ' 
c          write(kStdErr,*) 'AND try to use Cousin everywhere!'
c          CALL DoStop
c        END IF
c        print *,' ***** reset iLineMixBand = 1   ie linemix --> cousin ****' 
c        iLineMixBand = 1
c        END IF

      RETURN
      END

c************************************************************************
c               these are the PARTITION FUNCTIONS
c************************************************************************
c this finds the qfcn at LTE, and is mainly for rotational energy
c look at SPECTRA/qtips.m; this is for H96
c the latest GENLN2 version still seems to use these params
      SUBROUTINE qfcnH96(dT,iGasID,iISO,dPartitionFcn,dMass)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iGasID,iISO
      DOUBLE PRECISION   dT             !!temperature
c output params
      DOUBLE PRECISION dPartitionFcn,dMass

c local variables
      DOUBLE PRECISION qt,q296,a,b,c,d,T,rq,qaCO2_296(8),qaO3_296(5)
      INTEGER iFound

      DATA (qaCO2_296(iFound),iFound=1,8)/
c          CO2 626,         636,         628,         627,         
     +      .286219E+03, .576928E+03, .607978E+03, .354389E+04, 
c             638,         637,         828,         728;   
     +      .123528E+04, .714432E+04, .323407E+03, .376700E+04/

      DATA (qaO3_296(iFound),iFound=1,5)/
c        O3 666,          668,         686,         667,         676;
     +  .348186E+04, .746207E+04, .364563E+04, .430647E+05, .212791E+05/

      iFound = -1
      q296 = 1.0d0
      qt = 1.0d0

      IF (iGASID .EQ. 2) THEN
        IF (iISO .EQ. 1) THEN
          iFound = +1
          a = -.13617D+01
          b = .94899D+00
          c = -.69259D-03
          d =  .25974D-05
          dMass = 44.0 * 1.0d0
        ELSEIF (iISO .EQ. 2) THEN
          iFound = +1
          a = -.20631D+01
          b = .18873D+01
          c = -.13669D-02
          d = .54032D-05
          dMass = 45.0 * 1.0d0
        ELSE IF (iISO .EQ. 3) THEN
          iFound = +1
          a = -.29175D+01
          b = .20114D+01
          c = -.14786D-02
          d = .55941D-05
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 4) THEN
          iFound = +1
          a = -.16558D+02
          b = .11733D+02
          c = -.85844D-02
          d = .32379D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 5) THEN
          iFound = +1
          a = -.44685D+01
          b = .40330D+01
          c = -.29590D-02
          d = .11770D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 6) THEN
          iFound = +1
          a = -.26263D+02
          b = .23350D+02
          c = -.17032D-01
          d = .67532D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 7) THEN
          iFound = +1
          a = -.14811D+01
          b = .10667D+01
          c = -.78758D-03
          d = .30133D-05
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 8) THEN
          iFound = +1
          a = -.17600D+02
          b = .12445D+02
          c = -.91837D-02
          d = .34915D-04
          dMass = 46.0 * 1.0d0
        ELSE
          write (kStdErr,*) 'sorry, no qtips data for CO2 isotope ', iISO
        END IF
      END IF

      IF (iGASID .EQ. 3) THEN
        IF (iISO .EQ. 1) THEN
          iFound = +1
          a = -.16443D+03
          b =  .69047D+01
          c =  .10396D-01
          d =  .26669D-04 
          dMass = 48.0 * 1.0d0
        ELSE
          write (kStdErr,*) 'sorry, no qtips data for O3 isotope ', iISO
        END IF
      END IF

      IF (iFound .LT. 0) THEN
        write (kStdErr,*) 'Sorry qtips did not find gasID/isotope combination'
        CALL DoStop
      END IF

      t = dT
      qt = a + b*T + c*T*T + d*T*T*T
      t = 296.0 * 1.0d0
      q296 = a + b*T + c*T*T + d*T*T*T

c      print *,dT,q296,qaco2_296(iISO),q296/qaco2_296(iISO)

      dPartitionFcn = q296/qt

      RETURN
      END

c************************************************************************
c this finds the qfcn at LTE, and is mainly for rotational energy
c look at SPECTRA/qtips.m; this is for H98 and H2k versions
      SUBROUTINE qfcnH98(dT,iGasID,iISO,dPartitionFcn,dMass)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input params 
      INTEGER iGasID,iISO
      DOUBLE PRECISION   dT             !!temperature
c output params
      DOUBLE PRECISION dPartitionFcn,dMass

c local variables
      DOUBLE PRECISION qt,q296,a,b,c,d,T,rq
      INTEGER iFound

      iFound = -1
      q296 = 1.0d0
      qt = 1.0d0

      IF (iGASID .EQ. 2) THEN
        IF (iISO .EQ. 1) THEN
          iFound = +1
          a = -.21995D+01
          b =  .96751D+00
          c = -.80827D-03
          d =  .28040D-05
          dMass = 44.0 * 1.0d0
        ELSEIF (iISO .EQ. 2) THEN
          iFound = +1
          a =  -.38840D+01
          b =  .19263D+01
          c = -.16058D-02
          d =  .58202D-05
          dMass = 45.0 * 1.0d0
        ELSE IF (iISO .EQ. 3) THEN
          iFound = +1
          a = -.47289D+01
          b =  .20527D+01
          c = -.17421D-02
          d = .60748D-05
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 4) THEN
          iFound = +1
          a = -.27475D+02
          b = .11973D+02
          c = -.10110D-01
          d = .35187D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 5) THEN
          iFound = +1
          a = -.84191D+01
          b = .41186D+01
          c = -.34961D-02
          d = .12750D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 6) THEN
          iFound = +1
          a = -.48468D+02
          b = .23838D+02
          c = -.20089D-01
          d = .73067D-04
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 7) THEN
          iFound = +1
          a = -.22278D+01
          b = .10840D+01
          c = -.89718D-03
          d =  .32143D-05
          dMass = 46.0 * 1.0d0
        ELSE IF (iISO .EQ. 8) THEN
          iFound = +1
          a = -.29547D+02
          b = .12714D+02
          c =  -.10913D-01
          d = .38169D-04
          dMass = 46.0 * 1.0d0
        ELSE
          write (kStdErr,*) 'sorry, no qtips data for CO2 isotope ', iISO
        END IF
      END IF

      IF (iGASID .EQ. 3) THEN
        IF (iISO .EQ. 1) THEN
          iFound = +1
          a = -.13459D+03
          b =  .62255D+01
          c =  .14811D-01
          d =  .18608D-04 
          dMass = 48.0 * 1.0d0
        ELSE
          write (kStdErr,*) 'sorry, no qtips data for O3 isotope ', iISO
        END IF
      END IF

      IF (iFound .LT. 0) THEN
        write (kStdErr,*) 'Sorry qtips did not find gasID/isotope combination'
        CALL DoStop
      END IF

      t = dT
      qt = a + b*T + c*T*T + d*T*T*T
      t = 296.0 * 1.0d0
      q296 = a + b*T + c*T*T + d*T*T*T

      dPartitionFcn = q296/qt

      RETURN
      END

c************************************************************************
c this is copied from GENLN3, except that variables are explicitly typed
C  PROGRAM        QNLTE   SUBROUTINE
C-----------------------------------------------------------------------
C  PURPOSE        INTERNAL PARTITION CALCULATION:
C                 SEPARATE FOR VIBRATIONAL AND ROTATIONAL SUMS.
C                 USED IN THE NON-LTE CALCULATION
C-----------------------------------------------------------------------
C  VERSION        4.0   D.P. EDWARDS   96/07/1
C                 CODE FROM BOB GAMACHE
C                 PRIVATE COMMUNICATION JUNE 1996
c                 last changed 98-3-3
C-----------------------------------------------------------------------
C  ARGUMENT       IGAS   I*4  I/P  HITRAN GAS ID
C                 ISO    I*4  I/P  HITRAN ISO INDEX
C                 TP     I*4  I/P  PATH TEMPERATURE [K]
C                 SQR    I*4  O/P  LTE ROTATIONAL PARTITION SUM @ TP
C                 SQRSTD I*4  O/P  LTE ROTATIONAL PARTITION SUM @ 296K
C                 SQR    I*4  O/P  LTE VIBRATIONAL PARTITION SUM @ TP
C                 SQR    I*4  O/P  LTE VIBRATIONAL PARTITION SUM @ 296K
C***********************************************************************      
      SUBROUTINE QNLTE(IGAS,ISO,TP,SQR,SQRSTD,SQV,SQVSTD)

      IMPLICIT NONE
c input and output parameters
      integer  igas, iso
      double precision tp, sqr, sqrstd, sqv, sqvstd

c local variables
      DOUBLE PRECISION tstd,t,gj,qr,gjstd,qrstd,qv,qvstd
      INTEGER mol
C-----------------------------------------------------------------------
      TSTD = 296.D0
c...set temperature, molecule number, and isotope (HITRAN codes)
c...return with vibrational and rotational partition sums and state 
c...independent nuclear factors.
      T = TP
C  3   write(*,*) '   Enter temperature'
C      read(*,*) T
      if(T.lt.100. .or. T.gt.450.) then
        t=450.
c        write(*,*) '  Out of temperature range'
c        STOP 'QNLTE'
C        write(*,*) '  Enter a temperature between 100 and 450K'
C	go to 3
      endif

      MOL = IGAS
C
C 4    write(*,*) '   Enter Molecule number (HITRAN code)'
C      write(*,*) '   1 = H2O' 
C      write(*,*) '   2 = CO2' 
C      write(*,*) '   3 = O3' 
C      write(*,*) '   5 = CO' 
C      write(*,*) '   10 = NO2' 
C      read(*,*) Mol
      if(Mol.ne.1 .and. Mol.ne.2 .and. Mol.ne.3 .and. Mol.ne.4 
     +    .and. Mol.ne.5 .and. mol.ne.10) then
        WRITE(*,*) '  No nonLTE Q-data for molecule'
        STOP 'QNLTE'
C        write(*,*) '  Please enter molecule in list'
C        go to 4
      endif

  5   if(Mol.eq.1) then
C        write(*,*) '   Enter isotope code'
C        write(*,*) '   Enter 1 for 161'
C        write(*,*) '   Enter 2 for 181'
C        write(*,*) '   Enter 3 for 171'
C        write(*,*) '   Enter 4 for 162'
C        read(*,*) iso
        if(Iso.ne.1 .and. Iso.ne.2 .and. Iso.ne.3 .and. Iso.ne.4 
     +  ) then
          write(*,*) '  Please enter isotope in list'
          STOP 'QNLTE'
C          go to 5
        endif
      endif

  6   if(Mol.eq.2) then
C        write(*,*) '   Enter isotope code'
C        write(*,*) '   Enter 1 for 626'
C        write(*,*) '   Enter 2 for 636'
C        write(*,*) '   Enter 3 for 628'
C        write(*,*) '   Enter 4 for 627'
C        write(*,*) '   Enter 5 for 638'
C        write(*,*) '   Enter 6 for 637'
C        write(*,*) '   Enter 7 for 828'
C        write(*,*) '   Enter 8 for 827'
C        read(*,*) iso
        if(Iso.ne.1 .and. Iso.ne.2 .and. Iso.ne.3 .and. Iso.ne.4 
     +    .and. Iso.ne.5 .and. Iso.ne.6 .and. Iso.ne.7 .and. 
     +    Iso.ne.8) then
          write(*,*) '  Please enter isotope in list'
          STOP 'QNLTE'
C          go to 6
        endif
      endif

  7    if(Mol.eq.3) then
C        write(*,*) '   Enter isotope code'
C        write(*,*) '   Enter 1 for 666'
C        write(*,*) '   Enter 2 for 668'
C        write(*,*) '   Enter 3 for 686'
C        write(*,*) '   Enter 4 for 667'
C        write(*,*) '   Enter 5 for 676'
C        read(*,*) iso
        if(Iso.ne.1 .and. Iso.ne.2 .and. Iso.ne.3 .and. Iso.ne.4 
     +    .and. Iso.ne.5) then
          write(*,*) '  Please enter isotope in list'
          STOP 'QNLTE'
C          go to 7
        endif
      endif

  8   if(Mol.eq.5) then
C        write(*,*) '   Enter isotope code'
C        write(*,*) '   Enter 1 for 26'
C        write(*,*) '   Enter 2 for 36'
C        write(*,*) '   Enter 3 for 28'
C        write(*,*) '   Enter 4 for 27'
C        write(*,*) '   Enter 5 for 38'
C        write(*,*) '   Enter 6 for 37'
C        read(*,*) iso
        if(Iso.ne.1 .and. Iso.ne.2 .and. Iso.ne.3 .and. Iso.ne.4 
     +    .and. Iso.ne.5 .and. Iso.ne.6) then
          write(*,*) '  Please enter isotope in list'
          STOP 'QNLTE'
C          go to 8
        endif
      endif

      if(Mol.eq.10) then
        if(Iso.ne.1)then
          write(*,*) '   only the 646 species'
          STOP 'QNLTE'
C          go to 8
        endif
      endif

  900 format(2I5, 2F6.0, 2E12.6)

C        write(*,*) ' Mol, Iso,   T,   gj,     Qr,     Qv'

      if(mol.eq.1) then
        call QrH2O(iso, T, gj, Qr)
        CALL QRH2O(ISO, TSTD, GJSTD, QRSTD)
        call QvH2O(iso, T, Qv)
        CALL QVH2O(ISO, TSTD, QVSTD) 
C        write(*,900) Mol, Iso, T, gj, Qr, Qv
      endif

      if(mol.eq.2) then
        call QrCO2(iso, T, gj, Qr)
        CALL QRCO2(ISO, TSTD, GJSTD, QRSTD)
        call QvCO2(iso, T, Qv)
        CALL QVCO2(ISO, TSTD, QVSTD) 
C        write(*,900) Mol, Iso, T, gj, Qr, Qv
      endif

      if(mol.eq.3) then
        call QrO3(iso, T, gj, Qr)
        CALL QRO3(ISO, TSTD, GJSTD, QRSTD)
        call QvO3(iso, T, Qv)
        CALL QVO3(ISO, TSTD, QVSTD) 
C        write(*,900) Mol, Iso, T, gj, Qr, Qv
      endif

      if(mol.eq.5) then
        call QrCO(iso, T, gj, Qr)
        CALL QRCO(ISO, TSTD, GJSTD, QRSTD)
        call QvCO(iso, T, Qv)
        CALL QVCO(ISO, TSTD, QVSTD) 
C        write(*,900) Mol, Iso, T, gj, Qr, Qv
      endif

      if(mol.eq.10) then
        call QrNO2(iso, T, gj, Qr)
        CALL QRNO2(ISO, TSTD, GJSTD, QRSTD)
        call QvNO2(iso, T, Qv)
        CALL QVNO2(ISO, TSTD, QVSTD) 
C        write(*,900) Mol, Iso, T, gj, Qr, Qv
      endif
      
      SQR = QR
      SQRSTD = QRSTD
      SQV = QV
      SQVSTD = QVSTD

      end


c     *****************
      Subroutine QrH2O(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 161
c		! 2 = 181
c		! 3 = 171
c		! 4 = 162
     I T, 	! temperature in K 
     O gj, 	! state independent nuclear degeneracy factor
     O Qr)	! rotational partition function at T
 
      integer iso
      double precision t,gj,qr

      double precision xgj(4), Qcoef(4,4), Q296(4),tref,eps,qt
      integer j

      data xgj/ 1.,1.,6.,6./
c...        --   161                              
      data (Qcoef( 1,j),j=1,4)/-.53023D+01, .28404D+00,
     +                .12564D-02,-.56017D-06/
c...        --   181                              
      data (Qcoef( 2,j),j=1,4)/-.53509D+01, .28640D+00,
     +                .12671D-02,-.56496D-06/
c...        --   171                              
      data (Qcoef( 3,j),j=1,4)/-.31968D+02, .17117D+01,
     +                .75725D-02,-.33761D-05/
c...        --   162                              
      data (Qcoef( 4,j),j=1,4)/-.27586D+02, .13989D+01,
     +                .62715D-02,-.27996D-05/

c             H2O 161,         181,         171,         162                
      data q296/ .175795D+03, .175795D+03, .105061D+04, .863355D+03/

      DATA Tref,eps / 296., 1.D-05/
c
      gj = xgj(iso)
c...total internal sum at 296K
      IF (DABS(T-Tref) .lt. eps) THEN
         Qr = Q296(iso)
         go to 99
      endif
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qr = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end
c

c     *****************
      Subroutine QrCO2(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 626
c		! 2 = 636
c		! 3 = 628
c		! 4 = 627
c		! 5 = 638
c		! 6 = 637
c		! 7 = 828
c		! 8 = 728
     I T, 	! temperature in K 
     O gj, 	! state independent nuclear degeneracy factor
     O Qr)	! rotational partition function at T

      integer iso
      double precision t,gj,qr

      double precision xgj(8), Qcoef(8,4), Q296(8),tref,eps,qt
      integer j
 
      data xgj/ 1.,2.,1.,6.,2.,12.,1.,6./

c...        --   626                              
      data (Qcoef( 1,j),j=1,4)/ .16706D+00, .89058D+00,
     +                .99971D-05, .23271D-11/
c...        --   636                              
      data (Qcoef( 2,j),j=1,4)/ .33413D+00, .17811D+01,
     +                .19992D-04, .68786D-11/
c...        --   628                              
      data (Qcoef( 3,j),j=1,4)/ .33417D+00, .18877D+01,
     +                .21196D-04, .11156D-12/
c...        --   627                              
      data (Qcoef( 4,j),j=1,4)/ .20037D+01, .11014D+02,
     +                .12437D-03, .98216D-10/
c...        --   638                              
      data (Qcoef( 5,j),j=1,4)/ .66833D+00, .37755D+01,
     +                .42368D-04, .24009D-10/
c...        --   637                              
      data (Qcoef(6,j),j=1,4)/ .40097D+01, .22029D+02,
     +                .24703D-03, .58218D-10/
c...        --   828                              
      data (Qcoef(7,j),j=1,4)/ .16713D+00, .10020D+01,
     +                .11252D-04,-.72034D-12/
c...        --   728                              
      data (Qcoef(8,j),j=1,4)/ .20048D+01, .11683D+02,
     +                .13095D-03, .58038D-10/

c             CO2 626,         636,         628,         627,        
      data q296/ .264654D+03, .529281D+03, .560962D+03, .327317D+04,
c        638,         637,         828,         728
     +  .112194D+04, .654616D+04, .297752D+03, .347174D+04/

      DATA Tref,eps / 296., 1.D-05/
c
      gj = xgj(iso)
c...total internal sum at 296K
      IF (DABS(T-Tref) .lt. eps) THEN
         Qr = Q296(iso)
         go to 99
      endif
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        Qr = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qr = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end

c     *****************
      Subroutine QrO3(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 666
c		! 2 = 668
c		! 3 = 686
c		! 4 = 667
c		! 5 = 676
     I T, 	! temperature in K 
     O gj, 	! state independent nuclear degeneracy factor
     O Qr)	! rotational partition function at T

      integer iso
      double precision t,gj,qr

      double precision xgj(5), Qcoef(5,4), Q296(5),tref,eps,qt
      integer j

      data xgj/ 1.,1.,1.,6.,6./
c...        --   666                              
      data (Qcoef(1,j),j=1,4)/-.11668D+03, .53395D+01,
     +                .24438D-01,-.11317D-04/
c...        --   668                              
      data (Qcoef(2,j),j=1,4)/-.24923D+03, .11401D+02,
     +                .52186D-01,-.24168D-04/
c...        --   686                              
      data (Qcoef(3,j),j=1,4)/-.12178D+03, .55723D+01,
     +                .25504D-01,-.11811D-04/
c...        --   667                              
      data (Qcoef(4,j),j=1,4)/-.14489D+04, .66289D+02,
     +                .30342D+00,-.14052D-03/
c...        --   676                              
      data (Qcoef(5,j),j=1,4)/-.71583D+03, .32756D+02,
     +                .14992D+00,-.69429D-04/

c            O3   666,         668,         686,         667,
      data q296/ .331145D+04, .707103D+04, .345587D+04, .411127D+05, 
c                 676
     +           .203146D+05/

      DATA Tref,eps / 296., 1.D-05/
c
      gj = xgj(iso)
c...total internal sum at 296K
      IF (DABS(T-Tref) .lt. eps) THEN
         Qr = Q296(iso)
         go to 99
      endif
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qr = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end
c
c     *****************
      Subroutine QrCO(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 26
c		! 2 = 36
c		! 3 = 28
c		! 4 = 27
c		! 5 = 38
c		! 6 = 37
     I T, 	! temperature in K 
     O gj, 	! state independent nuclear degeneracy factor
     O Qr)	! rotational partition function at T

      integer iso
      double precision t,gj,qr

      double precision xgj(6), Qcoef(6,4), Q296(6),tref,eps,qt
      integer j

      data xgj/ 1.,2.,1.,6.,2.,12./
c...        --    26                              
      data (Qcoef(1,j),j=1,4)/ .33715D+00, .35986D+00,
     +                .28521D-05,-.80804D-10/
c...        --    36                              
      data (Qcoef(2,j),j=1,4)/ .67401D+00, .75291D+00,
     +                .59556D-05,-.15617D-09/
c...        --    28                              
      data (Qcoef(3,j),j=1,4)/ .33701D+00, .37790D+00,
     +                .29903D-05,-.80827D-10/
c...        --    27                              
      data (Qcoef(4,j),j=1,4)/ .20224D+01, .22153D+01,
     +                .17551D-04,-.50468D-09/
c...        --    38                              
      data (Qcoef(5,j),j=1,4)/ .67363D+00, .79248D+00,
     +                .62561D-05,-.15419D-09/
c...        --    37                              
      data (Qcoef(6,j),j=1,4)/ .40433D+01, .46404D+01,
     +                .36125D-04,-.10135D-08/

c              CO 26,          36,          28,          27,                
      data q296/ .107104D+03, .224055D+03, .112455D+03, .659264D+03,
c                 38,          37
     +           .235793D+03, .138074D+04/

      DATA Tref,eps / 296., 1.D-05/
c
      gj = xgj(iso)
c...total internal sum at 296K
      IF (DABS(T-Tref) .lt. eps) THEN
         Qr = Q296(iso)
         go to 99
      endif
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qr = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end

c     *****************
      Subroutine QrNO2(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 646
     I T, 	! temperature in K 
     O gj, 	! state independent nuclear degeneracy factor
     O Qr)	! rotational partition function at T
 

      integer iso
      double precision t,gj,qr

      double precision Qcoef(4),tref,eps,qt,q296
      integer j

c...        --   646                              
      data (Qcoef(j),j=1,4)/-.46394D+03, .21206D+02,
     +                .97410D-01,-.44482D-04/
c     NO2
      gj = 3.
      q296 =  .131941D+05

      DATA Tref,eps / 296., 1.D-05/
c
c...total internal sum at 296K
      IF (DABS(T-Tref) .lt. eps) THEN
         Qr = Q296
         go to 99
      endif
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qr = Qcoef(1) 
     +       + Qcoef(2)*T
     +       + Qcoef(3)*T*T 
     +       + Qcoef(4)*T*T*T 
c      
   99 return
      end

c     *****************
      Subroutine QvH2O(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 161
c		! 2 = 181
c		! 3 = 171
c		! 4 = 162
     I T, 	! temperature in K 
     O Qv)	! rotational partition function at T

      integer iso
      double precision T,qv

      double precision qcoef(4,4),qt
      integer j

c...        --   161                              
      data (Qcoef( 1,j),j=1,4)/ .99824D+00, .31601D-04,
     +               -.17994D-06, .33068D-09/
c...        --   181                              
      data (Qcoef( 2,j),j=1,4)/ .99823D+00, .31881D-04,
     +               -.18197D-06, .33539D-09/
c...        --   171                              
      data (Qcoef( 3,j),j=1,4)/ .99823D+00, .31751D-04,
     +               -.18101D-06, .33316D-09/
c...        --   162                              
c...        --   162                              
      data (Qcoef( 4,j),j=1,4)/ .99803D+00, .38536D-04,
     +               -.24214D-06, .49662D-09/

c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qv = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end
c

c     *****************
      Subroutine QvCO2(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 626
c		! 2 = 636
c		! 3 = 628
c		! 4 = 627
c		! 5 = 638
c		! 6 = 637
c		! 7 = 828
c		! 8 = 728
     I T, 	! temperature in K 
     O Qv)	! rotational partition function at T

      integer iso
      double precision t,qv

      double precision qcoef(8,4),tref,eps,qt 
      integer j

c...        --   626                              
      data (Qcoef( 1,j),j=1,4)/ .10428D+01,-.69155D-03,
     +                .27905D-05, .80679D-10/
c...        --   636                              
      data (Qcoef( 2,j),j=1,4)/ .10464D+01,-.75718D-03,
     +                .31239D-05,-.12129D-09/
c...        --   628                              
      data (Qcoef( 3,j),j=1,4)/ .10437D+01,-.70716D-03,
     +                .28643D-05, .48899D-10/
c...        --   627                              
      data (Qcoef( 4,j),j=1,4)/ .10432D+01,-.69871D-03,
     +                .28282D-05, .67812D-10/
c...        --   638                              
      data (Qcoef( 5,j),j=1,4)/ .10471D+01,-.77104D-03,
     +                .31995D-05,-.14088D-09/
c...        --   637                              
      data (Qcoef( 6,j),j=1,4)/ .10467D+01,-.76445D-03,
     +                .31635D-05,-.13131D-09/
c...        --   828                              
      data (Qcoef( 7,j),j=1,4)/ .10445D+01,-.72141D-03,
     +                .29381D-05, .24622D-10/
c...        --   728                              
      data (Qcoef( 8,j),j=1,4)/ .10445D+01,-.71842D-03,
     +                .29025D-05, .58059D-10/

c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qv = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end

c     *****************
      Subroutine QvO3(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 666
c		! 2 = 668
c		! 3 = 686
c		! 4 = 667
c		! 5 = 676
     I T, 	! temperature in K 
     O Qv)	! rotational partition function at T
 
      integer iso
      double precision t,qv

      double precision qcoef(5,4),tref,eps,qt
      integer j

c...        --   666                              
      data (Qcoef( 1,j),j=1,4)/ .10220D+01,-.30911D-03,
     +                .87191D-06, .15073D-08/
c...        --   668                              
      data (Qcoef( 2,j),j=1,4)/ .10242D+01,-.34657D-03,
     +                .10436D-05, .14088D-08/
c...        --   686                              
      data (Qcoef( 3,j),j=1,4)/ .10243D+01,-.34434D-03,
     +                .10054D-05, .14929D-08/
c...        --   667                              
      data (Qcoef( 4,j),j=1,4)/ .10220D+01,-.30911D-03,
     +                .87191D-06, .15073D-08/
c...        --   676                              
      data (Qcoef( 5,j),j=1,4)/ .10220D+01,-.30911D-03,
     +                .87191D-06, .15073D-08/
c
c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qv = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end
c
c     *****************
      Subroutine QvCO(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 26
c		! 2 = 36
c		! 3 = 28
c		! 4 = 27
c		! 5 = 38
c		! 6 = 37
     I T, 	! temperature in K 
     O Qv)	! rotational partition function at T

      integer iso
      double precision t,qv

      double precision qcoef(6,4),tref,eps,qt
      integer j

c...        --    26                              
      data (Qcoef( 1,j),j=1,4)/ .99939D+00, .10070D-04,
     +               -.50995D-07, .81541D-10/
c...        --    36                              
      data (Qcoef( 2,j),j=1,4)/ .99931D+00, .11310D-04,
     +               -.57586D-07, .92720D-10/
c...        --    28                              
      data (Qcoef( 3,j),j=1,4)/ .99931D+00, .11433D-04,
     +               -.58223D-07, .93783D-10/
c...        --    27                              
      data (Qcoef( 4,j),j=1,4)/ .99935D+00, .10761D-04,
     +               -.54661D-07, .87746D-10/
c...        --    38                              
      data (Qcoef( 5,j),j=1,4)/ .99923D+00, .12847D-04,
     +               -.65811D-07, .10682D-9/
c...        --    37                              
      data (Qcoef( 6,j),j=1,4)/ .99934D+00, .10926D-04,
     +               -.55541D-07, .89245D-10/

c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qv = Qcoef(iso,1) 
     +       + Qcoef(iso,2)*T
     +       + Qcoef(iso,3)*T*T 
     +       + Qcoef(iso,4)*T*T*T 
c      
   99 return
      end

c     *****************
      Subroutine QvNO2(
     I iso, 	! isotope code (HITRAN notation)
c		! 1 = 646
     I T, 	! temperature in K 
     O Qv)	! rotational partition function at T

      integer iso
      double precision t,qv

      double precision qcoef(4),tref,eps,qt
      integer j
 
c...        --   646                              
      data (Qcoef(j),j=1,4)/ .10145D+01,-.20470D-03,
     +                .59546D-06, .88374D-09/

c...value depends on temperature range
      if(T.lt.100. .OR. T.gt.450.) then
        QT = -1.
        write(*,*) '  OUT OF TEMPERATURE RANGE'
        go to 99
      endif

        Qv = Qcoef(1) 
     +       + Qcoef(2)*T
     +       + Qcoef(3)*T*T 
     +       + Qcoef(4)*T*T*T 
c      
   99 return
      end


c************************************************************************
