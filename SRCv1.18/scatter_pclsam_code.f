c Copyright 2006
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************************************************************************
c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output

      SUBROUTINE find_radiances_pclsam(iRadOrColJac,raFreq,raaAbs,iMRO,tcc,raCC,
     $                          iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,      	  
     $         raaExt,raaSSAlb,raaAsym,iKnowTP,
     $         iPhase,raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $            raaSumAbCoeff,caFluxFile,
     $            caJacobFile,caJacobFile2,
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob,
     $         raaRadsX,iNumOutX)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaExt     = matrix containing the mixed path abs coeffs + cloud ext
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for computing radiances
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
c TEMP        = tempertaure profile in terms of pressure levels
      REAL rSatAzimuth,rSolAzimuth
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),rSurfPress,raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iRadOrColJac
      INTEGER iNpmix,iFileID,iTag,iKnowTP
      REAL raaAbs(kMaxPts,kMixFilRows)
      CHARACTER*80 caOutName
      REAL Temp(MAXNZ)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
c this is to do with flux
      CHARACTER*80 caFluxFile
c this is to do with jacobians
      CHARACTER*80 caJacobFile,caJacobFile2
      INTEGER iNumGases,iaGases(kMaxGas),iNatm
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAmt(kProfLayer,kGasStore)
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ),iIOUN_USE,iIOUN_IN
c this is to do with cloud fracs
      INTEGER iNumOutX,iMRO
      REAL raaRadsX(kMaxPts,kProfLayer),tcc,raCC(kProfLayer)
      INTEGER iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
      REAL    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
      REAL    rCLrFrac               !! clear fraction
      INTEGER iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

      INTEGER i1,i2,iFloor,iDownWard

      !! --------- kAvgMin is a global variable in kcarta.param -------- !!
      !!kAvgMin is a global variable in kcarta.param .. set as required
      !!it is the average of single scattering albedo (w0); if less than some
      !!value, then basically there is no scattering and so can do some 
      !!approximations!!!!!
      kAvgMin = 1.0d-3     !!!before Feb 14, 2003
      kAvgMin = 1.0d-6
      !! --------- kAvgMin is a global variable in kcarta.param -------- !!

      !! --------- kTemperVary is a global variable in kcarta.param -------- !!
      !!  used in scatter_twostream but not in scatter_pclsam                !!

      DO i1=1,kMaxPts
        raInten(i1) = 0.0
        ENDDO

      IF (iRadOrColJac .EQ. 1) THEN 
        iIOUN_IN = kStdkCarta 
      ELSE 
        iIOUN_IN = kStdJacob2 
      END IF 

c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2=iaaRadLayer(iAtm,iNumLayer)-1
        i2=iFloor(i2*1.0/kProfLayer)
        IF (rTSpace .GT. 5.0) THEN
          write(kStdErr,*) 'you want satellite to be downward looking'
          write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
          write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
          write(kStdErr,*) 'Please retry'
          CALL DoSTOP
        END IF
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,iNumLayer))THEN
c radiation travelling downwards to instrument ==> sat looking up
c i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = -1
        i1=iaaRadLayer(iAtm,1)-1
        i1=iFloor(i1*1.0/(1.0*kProfLayer))
        i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
      END IF
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUpper=90,iLower=1 is acceptable
c eg iUpper=140,iLower=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
        CALL DoSTOP
      END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
      END IF

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      ! iRadOrColJac = -1 ==> do coljacs
      ! iRadOrColJac = +1 ==> do rads only

      IF ((iRadOrColJac .EQ. -1) .AND. (k100layerCloud .GT. 0)) THEN
        write(kStdWarn,*) ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'
        write(kStdErr,*)  ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'	
        Call DoStop
      ELSEIF ((iRadOrColJac .EQ. -1) .AND. (k100layerCloud .EQ. -1)) THEN
        write(kStdWarn,*) ' ---> Doing Column Gas AMT, STEMP Jacs ...'
        CALL rad_pclsam_coljac(raFreq,iDownward,
     $        iJacob,iaJacob,raaaColDQ,raaSumAbCoeff,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $            caJacobFile,caJacobFile2,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)

      ELSEIF (iRadOrColJac .EQ. +1) THEN
        IF ((iDownward .EQ. 1) .AND. (k100layerCloud .LE. 1)) THEN
          CALL rad_DOWN_pclsam_solar(raFreq,+1,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
        ELSEIF ((iDownward .EQ. 1) .AND. (k100layerCloud .EQ. 100) .AND. (abs(iMRO) .EQ. 1)) THEN
          CALL rad_DOWN_pclsam_solar100_simplemodel(raFreq,+1,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,
     $                          iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,      	  
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
        ELSEIF ((iDownward .EQ. 1) .AND. (k100layerCloud .EQ. 100) .AND. (abs(iMRO) .EQ. 2)) THEN
          CALL rad_DOWN_pclsam_solar100_MRO_driver(raFreq,+1,iKnowTP,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,
     $                          iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,      	  
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
        ELSE 
          CALL rad_UP_pclsam_solar(raFreq,+1,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,
     $         raaRadsX,iNumOutX)

        END IF
      END IF

      RETURN
      END

c************************************************************************
c this subroutine loops over finding 
c                         column jacobians (rQjac(1),rQjac(2),...)
c                         column temperature jacobian
c                         stemp jacobian (rST)
c for a 10 % gas amount perturbation
      SUBROUTINE rad_pclsam_coljac(raFreq,iDownward,
     $        iJacob,iaJacob,raaaColDQ,raaSumAbCoeff,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $            caJacobFile,caJacobFile2,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL rSatAzimuth,rSolAzimuth,raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iIOUN_USE
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
      CHARACTER*80 caJacobFile,caJacobFile2
c this is to do with cloud fracs
      INTEGER iNumOutX
      REAL raaRadsX(kMaxPts,kProfLayer)

c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ),iDownWard

      INTEGER iI,iL,iFr,iJ
      REAL raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)

      INTEGER iDefault,iColJac,iJacT,iJacB
      REAL rDefaultColMult,raMixVertTemp2(kMixFilRows)

      iDefault  = +1   !! do the (Stemp,col) Jacs

      iColJac = -1     !! skip  the (Stemp,col) Jacs (dump out zeros)
      iColJac = +1     !! do  the (Stemp,col) Jacs
      
      rDefaultColMult = kDefaultColMult

      !! remember we define iJacT and iJacB as radiating layer number with respect to SURFACE
      IF  (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN 
        !!downlook instrument : radiation going UP to instrument, very easy
        iJacT = kActualJacsT
        iJacB = kActualJacsB
      ELSE
        !!uplook instrument : radiation going DOWN to instrument, got to swap things
        iJacT = iNumlayer-kActualJacsB+1
        iJacB = iNumlayer-kActualJacsT+1
      END IF

      IF (iDefault .NE. iColJac) THEN 
        print *,'rad_main : col jacs : calculating numbers (slow) instead of '
        print *,' dumping out zeros (fast)'
        print *,'rad_main : col jacs : iDefault,iColJac = ',
     $    iDefault,iColJac
      END IF 
       
      IF (iColJac .EQ. +1) THEN
          !! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
          DO iJ = 1,iJacob
            write(kStdWarn,*) ' '
            write(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
            DO iL = 1,iNumLayer
              iI = iaaRadLayer(iAtm,iL)
c              IF ((iI .GE. kActualJacsB) .AND. (iI .LE. kActualJacsT)) THEN
              IF ((iL .GE. iJacB) .AND. (iL .LE. iJacT)) THEN
                write(kStdWarn,*) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
                DO iFr = 1,kMaxPts
                  raaTemp(iFr,iI) = raaExt(iFr,iI) + 
     $                              rDefaultColMult*raaaColDQ(iJ,iFr,iI)
                END DO
              ELSE
c                write(kStdWarn,*) 'not perturbing gas layer ',iI
                DO iFr = 1,kMaxPts
                  raaTemp(iFr,iI) = raaExt(iFr,iI)
                END DO
              END IF
            END DO

            IF (iDownward .EQ. 1) THEN
              CALL rad_DOWN_pclsam_solar(raFreq,-1,
     $         raInten,raVTemp,raaTemp,raaSSAlb,raaAsym,
     $         iPhase,raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,TEMP,
     $         iNp,iaOp,raaOp,iNpmix,iFileID,
     $         caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
            ELSE 
              CALL rad_UP_pclsam_solar(raFreq,-1,
     $          raInten,raVTemp,raaTemp,raaSSAlb,raaAsym,
     $          iPhase,raPhasePoints,raComputedPhase,
     $          ICLDTOPKCARTA, ICLDBOTKCARTA,
     $          rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $          rSatAngle,rFracTop,rFracBot,TEMP,
     $          iNp,iaOp,raaOp,iNpmix,iFileID,
     $          caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $          raSurface,raSun,raThermal,raSunRefl,
     $          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $          raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,
     $         raaRadsX,iNumOutX)
            END IF
          END DO

        !! do the column temperature jacobian radiance 
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
c        DO iL = 1,kMixFilRows
c          raMixVertTemp2(iL) = raVTemp(iL) + 1.0
c        END DO
        DO iL = 1,kMixFilRows
          raMixVertTemp2(iL) = 0.0
        END DO
        DO iL = 1,iNumLayer
          iI = iaaRadLayer(iAtm,iL)
c          IF ((iI .GE. kActualJacsB) .AND. (iI .LE. kActualJacsT)) THEN
          IF ((iL .GE. iJacB) .AND. (iL .LE. iJacT)) THEN
            write(kStdWarn,*) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
            raMixVertTemp2(iI) = raVTemp(iI) + 1.0
          ELSE
c            write(kStdWarn,*) 'not perturbing tempr layer ',iI
            raMixVertTemp2(iI) = raVTemp(iI)
          END IF
        END DO
        IF (iDownward .EQ. 1) THEN
           CALL rad_DOWN_pclsam_solar(raFreq,-1,
     $        raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
         ELSE 
           CALL rad_UP_pclsam_solar(raFreq,-1,
     $        raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,
     $         raaRadsX,iNumOutX)
         END IF

        !! do the stemp jacobian radiance 
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
        IF (iDownward .EQ. 1) THEN
           CALL rad_DOWN_pclsam_solar(raFreq,-1,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp+1,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)
         ELSE 
           CALL rad_UP_pclsam_solar(raFreq,-1,
     $        raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp+1,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,
     $         raaRadsX,iNumOutX)
         END IF

        ELSE
c          iIOUN_USE = kStdJacob2
c          write(kStdWarn,*) 'dump out zeros instead of col jacs/stemp jacs'
c          DO iFr = 1,kMaxPts
c            raJunk(iFr) = 0.0
c          END DO
c          DO iI = 1,iJacob+1
c            CALL wrtout(iIOUN_USE,caJacobFile2,raFreq,raJunk)
c          END DO
          !! just do nothing and dump out nothing, to save space
        END IF
        RETURN
      END 

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c allows for temperature variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

c WARNING raaExt,raaSSAlb,raaAsym are really adjusted optical depth,
c                                           single scattering albedo
c                                           asymmetry factor

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c this is for an UPLOOK instrument

c allows for temperature variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis
      SUBROUTINE rad_UP_pclsam_solar(raFreq,iRadOrColJac,raInten,
     $    raVTemp,raaExt,raaSSAlb,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,
     $    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,     
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,
     $         raaRadsX,iNumOutX)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL rSatAzimuth,rSolAzimuth
      REAL raSurFace(kMaxPts),raThermal(kMaxPts)
      REAL raSunScatter(kMaxPts),raSun(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers
c this is to do with NLTE
      INTEGER iNLTEStart,iRadorColJac
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
c this is to do with cloud fracs
      INTEGER iNumOutX
      REAL raaRadsX(kMaxPts,kProfLayer)

c local vars
      INTEGER iaRadLayer(kProfLayer),iFr,iFrX,iLay,iDp,iDoSolar,MP2Lay
      INTEGER iaCldLayer(kProfLayer),iIOUN,iSimple,iiDiv,iL,iLow
      INTEGER iCloudLayerTop,iCloudLayerBot,iLocalCldTop,iLocalCldBot
      REAL ttorad,raVT1(kMixFilRows),raOutfrac(kProfLayer),muSat,muSun
      REAL raInten2(kMaxPts),rAngleTrans,rAngleEmission,rPlanck,rMPTemp
      REAL rSolarScatter,raTau(kMaxPts),rSunAngle,rSunTemp,rOmegaSun
      REAL hg2_real,InterpTemp,rNoScale,rThermalRefl
      INTEGER iDefault,iDebugScatterJacobian

c these are from twostream code, so we know the reflection from surface!!!!!!
c just do the estimate using the twostream angles, rather than actual angle
c this will not be included in the Jacobian (at least, not now)
      REAL raRad1(kMaxPts),raRefl(kMaxPts),raTOA2GND(kMaxPts)
      REAL raSunDirect(kMaxPts)

      iDefault = +1         !!debugging,    use constant temp
      iDefault = -1         !!no debugging, use linear variation in temp

      iDebugScatterJacobian = +1
      iDebugScatterJacobian = -1

      IF (iDebugScatterJacobian .NE. iDefault) THEN
        write(kStdErr,*) 'In DoEmissionLinearInTau_Uplook have '
        write(kStdErr,*) 'iDebugScatterJacobian,iDefault = ',
     $                    iDebugScatterJacobian,iDefault
      END IF

      iNumOutX = 0

      rThermalRefl = 1.0/kPi
 
c calculate cos(SatAngle) 
      muSat = cos(rSatAngle*kPi/180.0) 
 
c if iDoSolar = 1, then include solar contribution from file 
c if iDoSolar = 0 then include solar contribution from T=5700K 
c if iDoSolar = -1, then solar contribution = 0 
      iSimple = nint(kProfLayer/2.0) 
      iDoSolar = kSolar 
      IF (kSolar .GE. 0) THEN 
        rSunAngle = raSunAngles(iSimple) 
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .LE. 1.0e-2) THEN 
          !!!do not want divergences in the code 
          rSunAngle = rSunAngle + 0.1 
          write(kStdWarn,*) 'Uplook instr : For PCLSAM code, reset sun angle'  
          write(kStdWarn,*) 'slightly different from satellite angle'  
          write(kStdErr,*) 'Uplook instr : For PCLSAM code, reset sun angle'  
          write(kStdErr,*) 'slightly different from satellite angle'  
          Call DoStop
        END IF 
      END IF 
 
      iDoSolar = kSolar 
 
c as we are never directly loooking at the sun, there is a geometry factor 
      rOmegaSun = kOmegaSun
      IF (iDoSolar .GE. 0) THEN 
        rSunTemp = kSunTemp 
        write(kStdWarn,*) 'upward looking instrument .. daytime' 
      ELSE IF (iDoSolar .LT. 0) THEN 
        rSunTemp = 0.0 
        write(kStdWarn,*)'upward looking instrument .. nitetime' 
      END IF 

      muSun = 1.0       !!!default  
      iLay  = iNumLayer
      iL    = iaRadLayer(iLay)
      IF (iDoSolar .GE. 0) THEN 
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
      END IF 
 
      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
      write(kStdWarn,*) 'iNumLayer,rTSpace,rTSurf = '
      write(kStdWarn,*)  iNumLayer,rTSpace,rTSurf
 
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      DO iLay=1,iNumLayer 
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay) 
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN 
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay) 
          CALL DoSTOP  
        END IF 
        IF (iaRadLayer(iLay) .LT. 1) THEN 
          write(kStdErr,*)'Error in forward model for atmosphere ',iAtm 
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay) 
          CALL DoSTOP  
        END IF 
      END DO 
 
cccccccccccccccccccc set these all important variables **************** 
      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      DO iLay = 1,kProfLayer
        iaCldLayer(iLay) = -1
      END DO
      IF ((ICLDTOPKCARTA .GT. 0) .AND. (ICLDBOTKCARTA .GT. 0)) THEN
        IF (iaRadLayer(1) .LT. kProfLayer) THEN 
          iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1 
          iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1 
          iiDiv = 0 
        ELSE 
          !!essentially do mod(iaRadLayer(1),kProfLayer) 
          iiDiv = 1           
 1010     CONTINUE 
          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN 
            iiDiv = iiDiv + 1 
            GOTO 1010 
          END IF 
          iiDiv = iiDiv - 1 
          iLay = iiDiv 
          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv) 
          iLocalCldTop = iiDiv - iCldTopkCarta + 1 
          iLocalCldBot = iiDiv - iCldBotkCarta + 1 
          iiDiv = iLay 
        END IF 
        DO iLay = iCldBotkCarta,iCldTopkCarta
          iaCldLayer(kProfLayer-iLay+1) = 1
        END DO
      ELSE
        write(kStdWarn,*)'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
      END IF

c      DO iLay = 1,iNumLayer
c        iL = iaRadLayer(iLay)
c        print *,iLay,iL,raaExt(1,iL),iaCldLayer(iLay)
c      END DO

cccccccccccccccccccc set these all important variables **************** 

c find the lowest layer that we need to output radiances for 
c note that since mixed paths are ordered 100,99,98 .. 1 here, we really 
c need to find the highest integer i.e. if we have to output radiances 
c at the 10,20 and 99 th layers in the atmosphere, we better loop down to 
c the 99th mixed path (which happens to be the layer just above ground) 
      iLow    = -1 
      DO iLay = 1,iNp 
        IF (iaOp(iLay) .GT. iLow) THEN 
          iLow = iaOp(iLay) 
        END IF 
      END DO 
      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers' 
      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer) 
      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow 
 
c set the temperature of the bottommost layer correctly 
      DO iFr=1,kMixFilRows 
        raVT1(iFr) = raVTemp(iFr) 
      END DO 

c if the bottom layer is fractional, interpolate!!!!!! 
      iL = iaRadLayer(iNumLayer) 
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL) 
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)  
c if the top layer is fractional, interpolate!!!!!! 
      iL = iaRadLayer(1) 
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)  
 
      IF (iDoSolar .EQ. 0) THEN 
        DO iFr=1,kMaxPts 
c NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually 
c compute the Plank radiation from the sun 
          raSunScatter(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO 
      ELSEIF (iDoSolar .EQ. 1) THEN 
        CALL ReadSolarData(raFreq,raSunScatter,iTag) 
      ELSE 
        DO iFr=1,kMaxPts 
          raSunScatter(iFr) = 0.0 
        END DO 
      END IF 
      !! keep this at 0, else it will complicate RadianceInterPolate bahaha
      DO iFr=1,kMaxPts 
        raSun(iFr) = 0.0 
      END DO 

      DO iFr=1,kMaxPts 
c initialize the diffuse downward contribution to 0 
c INTIALIZE the emission seen at satellite to 0.0 
        raInten(iFr)        = 0.0 
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual 
        raSurface(iFr) = ttorad(raFreq(iFr),rTSurf)
        raSunScatter(iFr) = raSunScatter(iFr) * rOmegaSun * muSun
        raSunDirect(iFr) = raSunScatter(iFr)
      END DO 

      DO iFr=1,kMaxPts 
c compute emission from the top of atm == eqn 4.26 of Genln2 manual
c initialize the cumulative thermal radiation 
        raThermal(iFr) = ttorad(raFreq(iFr),sngl(kTSpace))
        raRad1(iFr)    = raThermal(iFr)
      END DO 

c using Dave Turner's thesis ideas    
      CALL twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym,
     $                   iNumLayer,raRefl)

c now go from top of atmosphere down to the surface to compute the total 
c radiation from top of layer down to the surface 
c as we go from the top of the atmosphere downto the bottom, we keep the  
c cumulative effects (from layer iNumLayer to iLay) in each of  
c raThermal and raSolar  
 
c note that as direction of radiation travel is defined as 100,99,98,..,1 
c which is what is stored in iaRadLayer, we have to  
c      DO iLay=1,iNumLayer instead of DO iLay=iNumLayer,1,-1 
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer  

c --------------------------------------------------------------------->
c --------------------------------------------------------------------->

c xyz go from TOA to cldtop ------------------------------------------->
      DO iLay = 1,iLocalCldTop-1
        iL      = iaRadLayer(iLay) 
        muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 

c        print *,'above',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
c     $          rMPTemp,raThermal(1)
 
c see if this mixed path layer is in the list iaOp to be output    
c as we might have to do fractional layers!! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
c note this really messes up if it is a scattering layer, as it does not use
c cloudy sky rad transfer. Else it is fine for clear sky
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 

c now do the complete radiative transfer thru this layer 
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,
     $                           iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,
     $                           iaCldLayer,raaExt,raThermal) 

c see if we have to add on the solar contribution to do transmission thru atm 
        IF (iDoSolar .GT. 0) THEN 
c note that the angle is the solar angle = satellite angle 
          IF (iLay .EQ. 1) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracTop/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracTop
            END DO 
          ELSE IF (iLay .EQ. iNumLayer) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracBot/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracBot
            END DO 
          ELSE 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)
            END DO 
          END IF 

          !!! now see if we need the solar scatter term
          IF (iaCldLayer(iLay) .EQ. 1) THEN
            rNoScale = 1.0  !!! before Feb 2, 2006
            DO iFr = 1,kMaxPts
              rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
              rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) * 
     $                       (exp(-rNoScale*raTau(iFr)/muSat) - 
     $                        exp(-rNoScale*raTau(iFr)/muSun))
              rSolarScatter  = 
     $                      rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*
     $                         muSun/(muSat-muSun)/kForP
              raThermal(iFr) = raThermal(iFr) + rSolarScatter
            END DO
          END IF
        END IF  
      END DO 

c ------------------------------------------------------------------------>
      DO iFr = 1,kMaxPts
        raRad1(iFr) = raThermal(iFr)
      END DO

c go from cldtop to gnd, back to cldbot, refl back to gnd
      DO iLay = iLocalCldTop,iLow
        iL      = iaRadLayer(iLay) 
        muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 

c        print *,'haha down1',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,
c     $          rMPTemp,raRad1(1)

        DO iFr = 1,kMaxPts
          !!! simple radtrans
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raRad1(iFr) = raRad1(iFr) * exp(-raaExt(iFr,iL)/muSat) + 
     $             (1-exp(-raaExt(iFr,iL)/muSat))*rPlanck
        END DO
      END DO
        
c see if we have to add on the solar contribution from TOA to GND
      IF (iDoSolar .GT. 0) THEN 
        DO iLay = 1,iLow
          iL = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          DO iFr = 1,kMaxPts
            raTOA2GND(iFr) =  raTOA2GND(iFr) + raaExt(iFr,iL)/muSun
          END DO
        END DO
        DO iFr = 1,kMaxPts
          raSunDirect(iFr) =  raSunDirect(iFr)*exp(-raTOA2GND(iFr))
        END DO
      END IF

c pretend raRad1 is the background reflected thermal, and bounce it back up
      DO iFr=1,kMaxPts 
        raRad1(iFr) = 
     $         raSurface(iFr)*raUseEmissivity(iFr)+ 
     $         raRad1(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ 
     $         raSun(iFr)*raSunRefl(iFr) 
      END DO 

c go from gnd to cldbot
      DO iLay = iLow,iLocalCldBot+1,-1
        iL      = iaRadLayer(iLay) 
        muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 

c        print *,'haha up',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,
c     $          rMPTemp,raRad1(1)

        DO iFr = 1,kMaxPts
          !!! simple radtrans
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raRad1(iFr) = raRad1(iFr) * exp(-raaExt(iFr,iL)/muSat) + 
     $             (1-exp(-raaExt(iFr,iL)/muSat))*rPlanck
        END DO
      END DO

c refl from cldbottom
      DO iFr = 1,kMaxPts
        raRad1(iFr) = raRad1(iFr) * raRefl(iFr)
c        print *,iFr,raRefl(iFr)
      END DO
c ------------------------------------------------------------------------>

c xyz go from cldtop to cldbot ------------------------------------------->
      DO iLay = iLocalCldTop,iLocalCldBot
        iL      = iaRadLayer(iLay) 
        muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 

c        print *,'thru',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
c     $          rMPTemp,raThermal(1)
 
c see if this mixed path layer is in the list iaOp to be output    
c as we might have to do fractional layers!! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
c note this really messes up if it is a scattering layer, as it does not use
c cloudy sky rad transfer. Else it is fine for clear sky
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 

c now do the complete radiative transfer thru this layer 
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,
     $                           iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,
     $                           iaCldLayer,raaExt,raThermal) 

c see if we have to add on the solar contribution to do transmission thru atm 
        IF (iDoSolar .GT. 0) THEN 
c note that the angle is the solar angle = satellite angle 
          IF (iLay .EQ. 1) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracTop/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracTop
            END DO 
          ELSE IF (iLay .EQ. iNumLayer) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracBot/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracBot
            END DO 
          ELSE 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)
            END DO 
          END IF 

          !!! now see if we need the solar scatter term
          IF (iaCldLayer(iLay) .EQ. 1) THEN
            rNoScale = 1.0  !!! before Feb 2, 2006
            DO iFr = 1,kMaxPts
              rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
              rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) * 
     $                       (exp(-rNoScale*raTau(iFr)/muSat) - 
     $                        exp(-rNoScale*raTau(iFr)/muSun))
              rSolarScatter  = 
     $                      rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*
     $                         muSun/(muSat-muSun)/kForP
              raThermal(iFr) = raThermal(iFr) + rSolarScatter
            END DO
          END IF
        END IF  
      END DO 

c add on reflectance from cloud
      DO iFr = 1,kMaxPts
        raThermal(iFr) = raThermal(iFr) + raRad1(iFr)
      END DO

c xyz go from cldbot to gnd  ------------------------------------------->
      DO iLay = iLocalCldBot+1,iLow
        iL      = iaRadLayer(iLay) 
        muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 

c        print *,'below',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
c     $          rMPTemp,raThermal(1)
 
c see if this mixed path layer is in the list iaOp to be output    
c as we might have to do fractional layers!! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
c note this really messes up if it is a scattering layer, as it does not use
c cloudy sky rad transfer. Else it is fine for clear sky
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 

c now do the complete radiative transfer thru this layer 
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,
     $                           iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,
     $                           iaCldLayer,raaExt,raThermal) 

c see if we have to add on the solar contribution to do transmission thru atm 
        IF (iDoSolar .GT. 0) THEN 
c note that the angle is the solar angle = satellite angle 
          IF (iLay .EQ. 1) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracTop/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracTop
            END DO 
          ELSE IF (iLay .EQ. iNumLayer) THEN 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracBot/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)*rFracBot
            END DO 
          ELSE 
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0) 
            DO iFr=1,kMaxPts 
              rAngleTrans = exp(-raaExt(iFr,iL)/muSun) 
              raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans 
              raTau(iFr) = raaExt(iFr,iL)
            END DO 
          END IF 

          !!! now see if we need the solar scatter term
          IF (iaCldLayer(iLay) .EQ. 1) THEN
            rNoScale = 1.0  !!! before Feb 2, 2006
            DO iFr = 1,kMaxPts
              rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
              rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) * 
     $                       (exp(-rNoScale*raTau(iFr)/muSat) - 
     $                        exp(-rNoScale*raTau(iFr)/muSun))
              rSolarScatter  = 
     $                      rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*
     $                         muSun/(muSat-muSun)/kForP
              raThermal(iFr) = raThermal(iFr) + rSolarScatter
            END DO
          END IF
        END IF  
      END DO 

      !!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
      IF (kJacobian .GT. 0) THEN  
        !set raInten to rad at ground (instr) level
        DO iFr=1,kMaxPts 
          raInten(iFr) = raInten2(iFr) 
        END DO 
      END IF 
 
c set raSun = 0.0, and anything else is 
c scattered into the FOV as I do not have the "looking directly in" capability
c codes up yet (even though it is a simple eqn)
      DO iFr=1,kMaxPts 
        raSun(iFr) = 0.0
      END DO

csun      iNumOutX = iNumOutX + 1        
csun      CALL wrtout(iIOUN,caOutName,raFreq,raThermal) 

      RETURN
      END

c************************************************************************
c this subroutince computes "linear in tau" radiance for downlook instr
c uses single precision and a few tricks when the layer optical depth is small
      SUBROUTINE DoEmissionLinearInTau_Downlook(
     $                           iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,raaLayTrans,
     $                           iaCldLayer,raaExt,raaEmission)
      
      IMPLICIT NONE

      include '../INCLUDE/scatter.param' 

c input vars
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iaCldLayer(kProfLayer)
      REAL rFracTop,rFracBot,raLayAngles(kProfLayer)
      REAL raVT1(kMixFilRows),temp(maxnz),raFreq(kMaxPts)
      REAL raaLayTrans(kMaxPts,kProfLayer),raaExt(kMaxPts,kMixFilRows)
c output vars
      REAL raaEmission(kMaxPts,kProfLayer)

c local vars
      INTEGER iLay,iL,iFr,iLModKprofLayer,MP2Lay
      REAL ttorad,rPlanck,rSunTemp,rMPTemp,muSat
      REAL bup,rdn,bdn,rup,rL,rU,rFrac,r0,db
      INTEGER iDefault,iDebugScatterJacobian

      iDefault = +1         !!debugging,    constant temp, as in SARTA
      iDefault = -1         !!no debugging, use linear variation in temp

      iDebugScatterJacobian = +1
      iDebugScatterJacobian = -1

      IF (iDebugScatterJacobian .NE. iDefault) THEN
        write(kStdErr,*) 'In DoEmissionLinearInTau_Downlook have '
        write(kStdErr,*) 'iDebugScatterJacobian,iDefault = ',
     $   iDebugScatterJacobian,iDefault
      END IF

      DO iLay = 1,iNumLayer 
        iL = iaRadLayer(iLay) 
        IF (iLay .EQ. 1) THEN
          rFrac = rFracBot
        ELSEIF (iLay .EQ. iNumLayer) THEN
          rFrac = rFracTop
        ELSE
          rFrac = 1.0
        END IF
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
c first get the Mixed Path temperature for this radiating layer 
        rMPTemp = raVT1(iL) 
        rdn = temp(iL)     !!lower level
        rup = temp(iL+1)   !!upper level
        iLModKprofLayer = mod(iL,kProfLayer) 
        !normal, no LTE emission stuff 

        IF (iaCldLayer(iLay) .EQ. -1) THEN
          DO iFr=1,kMaxPts 
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck 
          END DO
        ELSEIF (iaCldLayer(iLay) .EQ. 1) THEN
          IF (iDebugScatterJacobian .EQ. -1) THEN
            !optical depth varies linearly with radiance
            DO iFr=1,kMaxPts 

              !lower level radiance
              bup = ttorad(raFreq(iFr),rdn)
              !upper level radiance
              bdn = ttorad(raFreq(iFr),rup)

              db = bup - bdn
              r0 = raaExt(iFr,iL)*rFrac/muSat

              IF (r0 .gt. 1.0e-4) THEN
                rL = (1.0 - raaLayTrans(iFr,iLay))
                rL = bup*rL - db*rL/r0 
                rU = db * raaLayTrans(iFr,iLay)
                raaEmission(iFr,iLay) = max(rL + rU,0.0)
              ELSE
                rL = r0
                rL = bup*rL - db
                rU = db * (1-r0)
                raaEmission(iFr,iLay) = max(rL + rU,0.0)
              END IF

            END DO
          ELSEIF (iDebugScatterJacobian .EQ. +1) THEN
            !no variation 
            DO iFr=1,kMaxPts 
              rPlanck = ttorad(raFreq(iFr),rMPTemp)
              raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck 
            END DO
          END IF

        END IF 
      END DO 
 
      RETURN
      END

c************************************************************************
c this subroutine computes "linear in tau" radiance for downlook instr
c uses single precision and a few tricks when the layer optical depth is small
      SUBROUTINE DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,
     $                        iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                        raLayAngles,raVT1,temp,raFreq,
     $                        iaCldLayer,raaExt,raThermal) 

      IMPLICIT NONE

      include '../INCLUDE/scatter.param' 

c input vars
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iaCldLayer(kProfLayer),iLay
      INTEGER iDebugScatterJacobian
      REAL rFracTop,rFracBot,raLayAngles(kProfLayer)
      REAL raVT1(kMixFilRows),temp(maxnz),raFreq(kMaxPts)
      REAL raaExt(kMaxPts,kMixFilRows)
c output vars
      REAL raThermal(kMaxPts)

c local vars
      INTEGER iL,iFr,iLModKprofLayer,MP2Lay
      REAL rAngleEmission,rAngleTrans
      REAL ttorad,rPlanck,rSunTemp,rMPTemp,muSat
      REAL bup,rdn,bdn,rup,rL,rU,rFrac,r0,db

      iL = iaRadLayer(iLay) 
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
 
      rMPTemp = raVT1(iL) 

      iLModKprofLayer = mod(iL,kProfLayer) 
      !normal, no LTE emission stuff 

      IF (iaCldLayer(iLay) .EQ. -1) THEN
        !!usual clear sky calcs
        IF (iLay .EQ. 1) THEN 
          DO iFr=1,kMaxPts 
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            rAngleTrans    = exp(-raaExt(iFr,iL)*rFracTop/muSat) 
            rAngleTrans    = exp(-raaExt(iFr,iL)/muSat) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
csun            rAngleEmission = 0.0
            raThermal(iFr) = raThermal(iFr)*rAngleTrans + rAngleEmission 
          END DO 
        ELSEIF (iLay .EQ. iNumLayer) THEN 
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	  
            rAngleTrans    = exp(-raaExt(iFr,iL)*rFracBot/muSat) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
csun            rAngleEmission = 0.0
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission 
          END DO 
        ELSE 
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	  	  
            rAngleTrans    = exp(-raaExt(iFr,iL)/muSat) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
csun            rAngleEmission = 0.0
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission 
          END DO 
        END IF 

      ELSEIF (iaCldLayer(iLay) .EQ. 1) THEN
        !!cloudy sky calcs
        rFrac = 1.0
        rdn = temp(iL)     !!lower level
        rup = temp(iL+1)   !!upper level

        IF (iLay .EQ. 1) THEN 
          rFrac = rFracTop
        ELSEIF (iLay .EQ. iNumLayer) THEN 
          rFrac = rFracBot
        ELSE
          rFrac = 1.0
        END IF

        IF (iDebugScatterJacobian .EQ. -1) THEN
          !! do the PCLSAM variation in radiance with optical depth
          DO iFr=1,kMaxPts 
            !lower level radiance
            bup = ttorad(raFreq(iFr),rdn)

            !upper level radiance
            bdn = ttorad(raFreq(iFr),rup)	    

            db = bup - bdn
            r0 = raaExt(iFr,iL)*rFrac/muSat

            rAngleTrans = exp(-r0)
  
            IF (r0 .gt. 1.0e-4) THEN
              rL = (1.0 - rAngleTrans)
              rL = bup*rL + db*rL/r0 
              rU = -db 
            ELSE
              rL = r0
              rL = bup*rL - db
              rU = -db
            END IF
            rAngleEmission = rL + rU
csun            rAngleEmission = 0.0
            raThermal(iFr) = raThermal(iFr)*exp(-r0) + rAngleEmission
          END DO

        ELSEIF (iDebugScatterJacobian .EQ. +1) THEN
          !no variation 
          DO iFr=1,kMaxPts 
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            rAngleTrans    = exp(-raaExt(iFr,iL)/muSat) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
csun            rAngleEmission = 0.0
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission 
          END DO
        END IF
      END IF
      RETURN
      END

c************************************************************************
c uses Dave Turner's thesis ideas
      SUBROUTINE twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym,
     $                         iNumLayer,raRefl)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iaCldLayer(kProfLayer)
c output vars
      REAL raRefl(kMaxPts)

c local vars
      INTEGER iLay,iL,iFr,iNumLayer
c these are from twostream code, so we know the reflection from surface!!!!!!
c just do the estimate using the twostream angles, rather than actual angle
c this will not be included in the Jacobian (at least, not now)
      REAL mu2str,rB,raAlpha(kMaxPts),raG(kMaxPts)
      REAL raKp(kMaxPts),raKm(kMaxPts),raW(kMaxPts)
      REAL raCp(kMaxPts),raCm(kMaxPts)
      REAL raG2C(kMaxPts),raDelta(kMaxPts)

      DO iFr=1,kMaxPts 
        raG2C(iFr)  = 0.0
        raRefl(iFr) = 0.0
        raW(iFr)    = 0.0
      END DO 

      DO iLay = 1,iNumLayer 
        iL = iaRadLayer(iLay) 
        IF (iaCldLayer(iLay) .EQ. 1) THEN
          DO iFr = 1,kMaxPts
            raG(iFr) = raaAsym(iFr,iL)
          END DO 
          GOTO 1234
        END IF
      END DO

 1234 CONTINUE

      DO iLay = 1,iNumLayer 
        iL = iaRadLayer(iLay) 
        IF (iaCldLayer(iLay) .EQ. 1) THEN
          print *,iLay,iL,iaCldLayer(iLay),raG(1),raaExt(1,iL)
          DO iFr = 1,kMaxPts
            raG2C(iFr) = raG2C(iFr) + raaExt(iFr,iL)
            raW(iFr)   = max(raW(iFr),raaSSAlb(iFr,iL))
          END DO 
        END IF
      END DO

c compute the reflection coefficient
      mu2str = 1.0/sqrt(3.0)
      DO iFr = 1,kMaxPts
        rB = (1-raG(iFr))/2.0
        raKp(iFr) = +1/mu2str*sqrt((1-raW(iFr))*(1-raW(iFr)*raG(iFr))) 
        raKm(iFr) = -1/mu2str*sqrt((1-raW(iFr))*(1-raW(iFr)*raG(iFr))) 
        raAlpha(iFr) = raW(iFr)*(1-rB) - 1 
        raCp(iFr) = -(raKp(iFr) + raAlpha(iFr)/mu2str)*mu2str/(raW(iFr)*rB) 
        raCm(iFr) = -(raKm(iFr) + raAlpha(iFr)/mu2str)*mu2str/(raW(iFr)*rB) 
        raDelta(iFr) = raCp(iFr)*exp(raKm(iFr)*raG2C(iFr)) - 
     $                   raCm(iFr)*exp(raKp(iFr)*raG2C(iFr)) 
        raRefl(iFr) = 
     $       (exp(raKm(iFr)*raG2C(iFr))-exp(raKp(iFr)*raG2C(iFr)))/raDelta(iFr)
        raRefl(iFr) = max(raRefl(iFr),0.0) 
        raRefl(iFr) = min(raRefl(iFr),1.0) 
      END DO

      RETURN
      END

c************************************************************************
c this is for k100layerCloud <= 1
      SUBROUTINE rad_DOWN_pclsam_solar(raFreq,iRadOrColJac,
     $    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,
     $    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL rSatAzimuth,rSolAzimuth
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iRadorColJac
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
c this is to do with cloud fracs
      INTEGER iNumOutX
      REAL raaRadsX(kMaxPts,kProfLayer)

c local variables
      INTEGER iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer)

c to do the thermal,solar contribution
      REAL rThermalRefl,radtot,rLayT,rEmission,rSunAngle
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

c to do fast NLTE
      REAL suncos,scos1,vsec1

c general
      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot
      INTEGER i1,i2,iLoop,iDebug
      INTEGER iSTopNormalRadTransfer
      REAL rFrac,rL,rU,r0
      REAL raCC(kProfLayer),rC
       
      iNumOutX = 0
      
      rThermalRefl = 1.0/kPi

c calculate cos(SatAngle)
      muSat = cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

cccccccccccccccccccc set these all important variables ****************
      iLocalCldTop = -1
      iLocalCldBot = -1
      DO iLay = 1,kProfLayer
        iaCldLayer(iLay) = -1   !!assume no cld
      END DO
      IF ((ICLDTOPKCARTA .GT. 0) .AND. (ICLDBOTKCARTA .GT. 0)) THEN
        IF (iaRadLayer(1) .LT. kProfLayer) THEN
          iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
          iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
          iiDiv = 0
        ELSE
          !!essentially do mod(iaRadLayer(1),kProfLayer)
          iiDiv = 1          
 1010     CONTINUE
          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
            iiDiv = iiDiv + 1
            GOTO 1010
          END IF
          iiDiv = iiDiv - 1
          iLay = iiDiv
          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
          iLocalCldTop = iCldTopkCarta - iiDiv + 1
          iLocalCldBot = iCldBotkCarta - iiDiv + 1
          iiDiv = iLay
        END IF
        DO iLay = iCldBotkCarta,iCldTopkCarta
          iaCldLayer(iLay) = 1
        END DO
      ELSE
        write(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
      END IF

cccccccccccccccccccc set these all important variables ****************
             
c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
      
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh=iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)*rFracBot/muSat)
          raaEmission(iFr,iLay) = 0.0
        END DO
      END DO
      DO iLay=2,iNumLayer-1
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)/muSat)
          raaEmission(iFr,iLay) = 0.0
        END DO
      END DO
      DO iLay=iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)*rFracTop/muSat)
          raaEmission(iFr,iLay) = 0.0
        END DO
      END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer 
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!  
c note iNLTEStart = kProfLayer + 1, so only LTE is done 
      iNLTEStart = kProfLayer + 1 
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere 
      iUpper = -1 
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE' 
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer 

      CALL DoEmissionLinearInTau_Downlook(
     $                           iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,raaLayTrans,
     $                           iaCldLayer,raaExt,raaEmission)

c now go from top of atmosphere down to the surface to compute the total 
c radiation from top of layer down to the surface 
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance 
c from the top of atmosphere is not reflected 
      IF (iDoThermal .GE. 0) THEN 
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, 
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1) 
      ELSE 
        write(kStdWarn,*) 'no thermal backgnd to calculate' 
      END IF 
 
c see if we have to add on the solar contribution 
      IF (iDoSolar .GE. 0) THEN
        !this figures out the solar intensity at the ground, for reflection up
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag) 
        !this figures backscattered solar intensity
        iSolarRadOrJac = +1   !!! compute rad
        CALL SolarScatterIntensity_Downlook(
     $      iDoSolar,raFreq,iaCldLayer,
     $      raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth,
     $      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot,
     $      iTag,iSolarRadorJac,raaSolarScatter1Lay)
      ELSE 
        write(kStdWarn,*) 'no solar backgnd to calculate' 
        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            raaSolarScatter1Lay(iFr,iLay) = 0.0
          END DO
        END DO
      END IF 

      DO iFr=1,kMaxPts 
        raInten(iFr) = 
     $         raSurface(iFr)*raUseEmissivity(iFr)+ 
     $         raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ 
     $         raSun(iFr)*raSunRefl(iFr) 
      END DO 

 4321 FORMAT(I5,' ',7(F10.4,' ')) 
 
      r0 = raInten(1) 
c now we can compute the upwelling radiation!!!!! 
c compute the total emission using the fast forward model, only looping  
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh 
 
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional) 
      DO iLay=1,1 
         iL = iaRadLayer(iLay) 
         muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
         rMPTemp = raVT1(iL) 
c see if this mixed path layer is in the list iaOp to be output 
c since we might have to do fractions! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 
 
c now do the radiative transfer thru this bottom layer 
        r0 = raInten(9523) 
        DO iFr=1,kMaxPts 
         raInten(iFr) = raaEmission(iFr,iLay) + 
     $                  raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                  raaSolarScatter1Lay(iFr,iL)
csun          raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
csun     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO 
      END DO 


c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full) 
      DO iLay=2,iHigh-1 
         iL = iaRadLayer(iLay) 
         muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
         rMPTemp = raVT1(iL) 
c         print *,iLay,rMPTemp,raaExt(8000,iL),raaLayTrans(8000,iLay) 
c see if this mixed path layer is in the list iaOp to be output 
c since we might have to do fractions! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 
 
c now do the radiative transfer thru this complete layer 

        r0 = raInten(9523) 
        DO iFr=1,kMaxPts 
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                   raaSolarScatter1Lay(iFr,iL)
csun          raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) +
csun     $                   raaSolarScatter1Lay(iFr,iL)
        END DO 
      END DO 
 
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^ 
c then do the topmost layer (could be fractional) 
 777  CONTINUE 
      DO iLay=iHigh,iHigh 
        iL = iaRadLayer(iLay) 
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 
        r0 = raInten(9523) 

        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 

        IF (iDoSolar .LT. 0) THEN
          IF (iDp .GT. 0) THEN 
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
            DO iFr=1,iDp 
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, 
     $            raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $            iProfileLayers,raPressLevels, 
     $            iNLTEStart,raaPlanckCoeff) 
              iNumOutX = iNumOutX + 1           
              DO iFrX = 1,kMaxPts
                raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
              END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
            END DO 
          END IF 
        ELSE
          IF (iDp .EQ. 1) THEN 	  
            write(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer' 

            suncos = raSunAngles(iaRadLayer(1))           !! at surface
            scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
            vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

            suncos = cos(suncos*kPi/180.0)
            scos1  = cos(scos1*kPi/180.0)
            vsec1  = 1/cos(vsec1*kPi/180.0)

            DO iFr=1,kMaxPts
              raInten2(iFr) = raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
	  
            CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1,
     $                  iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
            iNumOutX = iNumOutX + 1           
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
	  
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
	  ELSEIF (iDp .GT. 1) THEN
	    write(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
	    CALL DoStop
	  END IF
        END IF
       
csun          DO iFr=1,kMaxPts 
csun            raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
csun     $                     raInten(iFr)*raaLayTrans(iFr,iLay)
csun          END DO
csun          CALL wrtout(iIOUN,caOutName,raFreq,raInten) 
 
cc no need to do radiative transfer thru this layer 
cc        DO iFr=1,kMaxPts 
cc          raInten(iFr) = raaEmission(iFr,iLay)+ 
cc     $                   raaSolarScatter1Lay(iFr,iL) +
cc     $                   raInten(iFr)*raaLayTrans(iFr,iLay) 
cc        END DO 
      END DO      !!       DO iLay=iHigh,iHigh 
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^ 

      RETURN
      END

c************************************************************************
c this is for k100layerCloud == 100
c so we have to do two simulataneous runs,
c        one for clouds + gas raaExt
c        one for gas only     raaAbs
      SUBROUTINE rad_DOWN_pclsam_solar100_simplemodel(raFreq,iRadOrColJac,
     $    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,
     $                          iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,      	  
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,
     $    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL rSatAzimuth,rSolAzimuth
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      REAL raaAbs(kMaxPts,kMixFilRows)      
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iRadorColJac
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
c this is to do with cloud fracs
      INTEGER iNumOutX,iMRO
      REAL raaRadsX(kMaxPts,kProfLayer),tcc,raCC(KProfLayer)
      INTEGER iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
      REAL    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
      REAL    rCLrFrac               !! clear fraction
      INTEGER iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

c local variables
      INTEGER iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
      REAL raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
      REAL ttorad,rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer)

c to do the thermal,solar contribution
      REAL rThermalRefl,radtot,rLayT,rEmission,rSunAngle
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

c to do fast NLTE
      REAL suncos,scos1,vsec1

c general
      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot
      INTEGER i1,i2,iLoop,iDebug
      INTEGER iSTopNormalRadTransfer
      REAL rFrac,rL,rU,r0
c this is for the cloudy/clear streams
      REAL raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
      REAL raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
      REAL raaExtWeighted(kMaxPts,kMixFilRows)
      REAL raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)

      iNumOutX = 0
      
      rThermalRefl = 1.0/kPi

      IF (abs(iMRO) .EQ. 1) THEN
        write(kStdWarn,*) 'Simple 100 layer cloud model uses raCC and tcc'
      ELSE
        write(kStdErr,*) 'this routine is for simple 100 layer cloud model (iMRO = +/-1) ',iMRO
        CALL DoStop
      END IF
      
c calculate cos(SatAngle)
      muSat = cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

cccccccccccccccccccc set these all important variables ****************
      iLocalCldTop = -1
      iLocalCldBot = -1
      DO iLay = 1,kProfLayer
        iaCldLayer(iLay) = -1   !!assume no cld
      END DO
      IF ((ICLDTOPKCARTA .GT. 0) .AND. (ICLDBOTKCARTA .GT. 0)) THEN
        IF (iaRadLayer(1) .LT. kProfLayer) THEN
          iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
          iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
          iiDiv = 0
        ELSE
          !!essentially do mod(iaRadLayer(1),kProfLayer)
          iiDiv = 1          
 1010     CONTINUE
          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
            iiDiv = iiDiv + 1
            GOTO 1010
          END IF
          iiDiv = iiDiv - 1
          iLay = iiDiv
          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
          iLocalCldTop = iCldTopkCarta - iiDiv + 1
          iLocalCldBot = iCldBotkCarta - iiDiv + 1
          iiDiv = iLay
        END IF
        DO iLay = iCldBotkCarta,iCldTopkCarta
          iaCldLayer(iLay) = 1
        END DO
      ELSE
        write(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
      END IF

cccccccccccccccccccc set these all important variables ****************
             
c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
      
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh=iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
      DO iLay=1,1
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)*rFracBot/muSat)
          raaEmission(iFr,iLay) = 0.0
          raaLayTransGasOnly(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/muSat)
          raaEmissionGasOnly(iFr,iLay) = 0.0
        END DO
      END DO
      DO iLay=2,iNumLayer-1
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)/muSat)
          raaEmission(iFr,iLay) = 0.0
          raaLayTransGasOnly(iFr,iLay) = exp(-raaAbs(iFr,iL)/muSat)
          raaEmissionGasOnly(iFr,iLay) = 0.0
        END DO
      END DO
      DO iLay=iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
          raaLayTrans(iFr,iLay) = exp(-raaExt(iFr,iL)*rFracTop/muSat)
          raaEmission(iFr,iLay) = 0.0
          raaLayTransGasOnly(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/muSat)
          raaEmissionGasOnly(iFr,iLay) = 0.0
        END DO
      END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raSunGasOnly(iFr)     = 0.0
        raThermalGasOnly(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer 
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!  
c note iNLTEStart = kProfLayer + 1, so only LTE is done 
      iNLTEStart = kProfLayer + 1 
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere 
      iUpper = -1 
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE' 
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer 

      CALL DoEmissionLinearInTau_Downlook(
     $                           iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,raaLayTrans,
     $                           iaCldLayer,raaExt,raaEmission)
      CALL DoEmissionLinearInTau_Downlook(
     $                           iNumLayer,iaRadLayer,rFracTop,rFracBot,
     $                           raLayAngles,raVT1,temp,raFreq,raaLayTransGasOnly,
     $                           iaCldLayer,raaAbs,raaEmissionGasOnly)

c now go from top of atmosphere down to the surface to compute the total 
c radiation from top of layer down to the surface 
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance 
c from the top of atmosphere is not reflected 
      IF (iDoThermal .GE. 0) THEN 
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, 
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1) 
        CALL BackGndThermal(raThermalGasOnly,raVT1,rTSpace,raFreq, 
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,-1) 
      ELSE 
        write(kStdWarn,*) 'no thermal backgnd to calculate' 
      END IF 
 
c see if we have to add on the solar contribution 
      IF (iDoSolar .GE. 0) THEN
        !this figures out the solar intensity at the ground, for reflection up
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag) 
        CALL Solar(iDoSolar,raSunGasOnly,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag) 
        !this figures backscattered solar intensity
        iSolarRadOrJac = +1   !!! compute rad
        CALL SolarScatterIntensity_Downlook(
     $      iDoSolar,raFreq,iaCldLayer,
     $      raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth,
     $      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot,
     $      iTag,iSolarRadorJac,raaSolarScatter1Lay)
      ELSE 
        write(kStdWarn,*) 'no solar backgnd to calculate' 
        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            raaSolarScatter1Lay(iFr,iLay) = 0.0
          END DO
        END DO
      END IF 

      DO iFr = 1,kMaxPts
        raSun(iFr)     = tcc * raSun(iFr)     + (1.0 - tcc) * raSunGasOnly(iFr)
        raThermal(iFr) = tcc * raThermal(iFr) + (1.0 - tcc) * raThermalGasOnly(iFr)	
      END DO
      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaExtWeighted(iFr,iLay) = raCC(iLay) * raaExt(iFr,iLay) + (1.0 - raCC(iLay)) * raaAbs(iFr,iLay)
        END DO
      END DO

c at ground raInten (from cloud) = raIntenGasOnly = raaIntenWeighted
      DO iFr=1,kMaxPts 
        raInten(iFr) = 
     $         raSurface(iFr)*raUseEmissivity(iFr)+ 
     $         raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ 
     $         raSun(iFr)*raSunRefl(iFr)
       raIntenGasOnly(iFr) = raInten(iFr)
       raIntenWeighted(iFr) = raInten(iFr)
      END DO 

 4321 FORMAT(I5,' ',7(F10.4,' ')) 
 
      r0 = raInten(1) 
c now we can compute the upwelling radiation!!!!! 
c compute the total emission using the fast forward model, only looping  
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh 
 
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional) 
      DO iLay=1,1 
         iL = iaRadLayer(iLay) 
         muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
         rMPTemp = raVT1(iL) 
c see if this mixed path layer is in the list iaOp to be output 
c since we might have to do fractions! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 
 
c now do the radiative transfer thru this bottom layer 
        DO iFr=1,kMaxPts
c	 IF ((iFr .Eq. 1) .AND. (iLay .EQ. 1)) THEN	 
c           print *,1,iL,'<-',raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
c	 END IF	
         raInten(iFr) = raaEmission(iFr,iLay) + 
     $                  raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                  raaSolarScatter1Lay(iFr,iL)
         raIntenGasOnly(iFr) = raaEmissionGasOnly(iFr,iLay) + 
     $                  raIntenGasOnly(iFr)*raaLayTransGasOnly(iFr,iLay)
         raIntenWeighted(iFr) = raCC(iL) * raInten(iFr) + (1-raCC(iL)) * raIntenGasOnly(iFr)
c	 IF ((iFr .Eq. 1) .AND. (iLay .EQ. 1)) THEN	 
c           print *,1,iL,'->',raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
c	 END IF
	 raInten(iFr) = raIntenWeighted(iFr)
	 raIntenGasOnly(iFr) = raIntenWeighted(iFr) 	 
csun          raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
csun     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
      END DO 

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full) 
      DO iLay=2,iHigh-1 
         iL = iaRadLayer(iLay) 
         muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
         rMPTemp = raVT1(iL) 
c         print *,iLay,rMPTemp,raaExt(8000,iL),raaLayTrans(8000,iLay) 
c see if this mixed path layer is in the list iaOp to be output 
c since we might have to do fractions! 
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 
        IF (iDp .GT. 0) THEN 
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
          DO iFr=1,iDp 
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $        raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, 
     $        raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $        iProfileLayers,raPressLevels, 
     $        iNLTEStart,raaPlanckCoeff) 
            iNumOutX = iNumOutX + 1
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
          END DO 
        END IF 
 
c now do the radiative transfer thru this complete layer 

        r0 = raInten(9523) 
        DO iFr=1,kMaxPts 
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                   raaSolarScatter1Lay(iFr,iL)
         raIntenGasOnly(iFr) = raaEmissionGasOnly(iFr,iLay) + 
     $                  raIntenGasOnly(iFr)*raaLayTransGasOnly(iFr,iLay)
         raIntenWeighted(iFr) = raCC(iL) * raInten(iFr) + (1-raCC(iL)) * raIntenGasOnly(iFr)
c	 IF ((iFr .Eq. 1) .AND. (iLay .EQ. iHigh-1)) THEN
c           print *,iL,raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
c	 END IF	 
	 raInten(iFr) = raIntenWeighted(iFr)
	 raIntenGasOnly(iFr) = raIntenWeighted(iFr) 	      
csun          raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) +
csun     $                   raaSolarScatter1Lay(iFr,iL)
        END DO 
      END DO 

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^ 
c then do the topmost layer (could be fractional) 
 777  CONTINUE 
      DO iLay=iHigh,iHigh 
        iL = iaRadLayer(iLay) 
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
        rMPTemp = raVT1(iL) 
        r0 = raInten(9523) 

        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum) 

        IF (iDoSolar .LT. 0) THEN
          IF (iDp .GT. 0) THEN 
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer' 
            DO iFr=1,iDp 
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, 
     $            raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, 
     $            raSun,-1,iNumLayer,rFracTop,rFracBot, 
     $            iProfileLayers,raPressLevels, 
     $            iNLTEStart,raaPlanckCoeff) 
              iNumOutX = iNumOutX + 1           
              DO iFrX = 1,kMaxPts
                raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
              END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2) 
            END DO 
          END IF 
        ELSE
          IF (iDp .EQ. 1) THEN 	  
            write(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer' 

            suncos = raSunAngles(iaRadLayer(1))           !! at surface
            scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
            vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

            suncos = cos(suncos*kPi/180.0)
            scos1  = cos(scos1*kPi/180.0)
            vsec1  = 1/cos(vsec1*kPi/180.0)

            !!! assume no clouds at TOA, so no need to brseak it into clear/cloudy streams
            DO iFr=1,kMaxPts
              raInten2(iFr) = raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
	  
            CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1,
     $                  iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
            iNumOutX = iNumOutX + 1           
            DO iFrX = 1,kMaxPts
              raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
            END DO
	  
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
	  ELSEIF (iDp .GT. 1) THEN
	    write(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
	    CALL DoStop
	  END IF
        END IF
       
csun          DO iFr=1,kMaxPts 
csun            raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
csun     $                     raInten(iFr)*raaLayTrans(iFr,iLay)
csun          END DO
csun          CALL wrtout(iIOUN,caOutName,raFreq,raInten) 
 
cc no need to do radiative transfer thru this layer 
cc        DO iFr=1,kMaxPts 
cc          raInten(iFr) = raaEmission(iFr,iLay)+ 
cc     $                   raaSolarScatter1Lay(iFr,iL) +
cc     $                   raInten(iFr)*raaLayTrans(iFr,iLay) 
cc        END DO 
      END DO      !!       DO iLay=iHigh,iHigh 
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^ 

      RETURN
      END

c************************************************************************
c this is for k100layerCloud == 100
c so we have to do MANY subpixel simulataneous runs,
c                       for clouds + gas raaExt
c        and finally for gas only        raaAbs
      SUBROUTINE rad_DOWN_pclsam_solar100_MRO_driver(raFreq,iRadOrColJac,iKnowTP,
     $    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,
     $                          iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,      	  
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,
     $    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         raaRadsX,iNumOutX)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL rSatAzimuth,rSolAzimuth
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      REAL raaAbs(kMaxPts,kMixFilRows)      
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iKnowTP
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iRadorColJac
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
c this is to do with cloud fracs
      INTEGER iNumOutX,iMRO
      REAL raaRadsX(kMaxPts,kProfLayer),tcc,raCC(KProfLayer)
      INTEGER iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
      REAL    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
      REAL    rCLrFrac               !! clear fraction
      INTEGER iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

c local variables
      INTEGER iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
      REAL raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
      REAL ttorad,rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer)

c to do the thermal,solar contribution
      REAL rThermalRefl,radtot,rLayT,rEmission,rSunAngle
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

c to do fast NLTE
      REAL suncos,scos1,vsec1

c general
      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot
      INTEGER i1,i2,iLoop,iDebug
      INTEGER iSTopNormalRadTransfer
      REAL rFrac,rL,rU,r0
c this is for the cloudy/clear streams
      REAL raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
      REAL raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
      REAL raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)

c BIG ASSUMPTION : we are only interested in TOA radiances
      REAL raWeightedRadiance(kMaxPts)
      REAL rEps
      REAL raaTempAbs(kMaxPts,kMixFilRows)
      INTEGER iCldSubPixel,iaSwap(kProfLayer),iNumSwap

      IF (iOutNum .GT. 1) THEN
        write(kStdErr,*) 'rad_DOWN_pclsam_solar100_MRO_driver assumes only radiance at TOA will be dumped out'
        CALL DoStop
      END IF
      
      DO iFr = 1,kMaxPts
        raWeightedRadiance(iFr) = 0.0
      END DO

c clear sky
      IF (rClrfrac .GE. rEps) THEN
        IF (kOuterLoop .EQ. 1) write(kStdWarn,*) 'MRO : clrfrac = ',rClrfrac
	print *,'MRO ClearFrac ',rClrfrac,' for ',raFreq(1),' cm-1'
        CALL quick_clear_radtrans_downlook(
     $      raFreq,raInten,raVTemp,
     $      raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $      caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $      raThickness,raPressLevels,iProfileLayers,pProf,
     $      raTPressLevels,iKnowTP,rCO2MixRatio,
     $         raaRadsX,iNumOutX,-1)
	DO iFr = 1,kMaxPts
	  raWeightedRadiance(iFr) = raWeightedRadiance(iFr) + rClrfrac*raInten(iFr)
        END DO
      END IF

c loop over cloud subpixels
      DO iCldSubPixel = 1,iNumSubPixels
        IF (kOuterLoop .EQ. 1) write(kStdWarn,*) 'MRO : index/cldfrac = ',iCldSubPixel,raCfrac(iCldSubPixel)

	IF (iCldSubPixel .EQ. 1) THEN
	  ! set the ODs to gas ODS	
	  DO iLay = 1,kProfLayer
            DO iFr = 1,kMaxPts
              raaTempAbs(iFr,iLay) = raaAbs(iFr,iLay)
            END DO
	  END DO
        END IF
	
        ! find which lays to swap
	iNumSwap = 0
	DO iLay = 1,kProfLayer
	  IF (iaaCldLaySubPixel(iLay,iCldSubPixel) .EQ. 1) THEN
	    iNumSwap = iNumSwap + 1
	    iaSwap(iNumSwap) = iLay
	  END IF
	END DO

        !swap in necessary cldlays
	DO iLay = 1,iNumSwap
	  iL = iaSwap(iLay)
          DO iFr = 1,kMaxPts
            raaTempAbs(iFr,iL) = raaExt(iFr,iL)
          END DO
	END DO	
        CALL quick_clear_radtrans_downlook(
     $      raFreq,raInten,raVTemp,
     $      raaTempAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $      caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $      raThickness,raPressLevels,iProfileLayers,pProf,
     $      raTPressLevels,iKnowTP,rCO2MixRatio,
     $         raaRadsX,iNumOutX,-1)
	DO iFr = 1,kMaxPts
	  raWeightedRadiance(iFr) = raWeightedRadiance(iFr) + raCfrac(iCldSubPixel)*raInten(iFr)
        END DO
        write(kStdErr,111) raFreq(1),iCldSubPixel,iNumSubPixels,iNumSwap,raWeightedRadiance(1)

	IF (iCldSubPixel .LT. iNumSubPixels) THEN
          !swap back in original gasODs
  	  DO iLay = 1,iNumSwap
 	    iL = iaSwap(iLay)
            DO iFr = 1,kMaxPts
              raaTempAbs(iFr,iL) = raaAbs(iFr,iL)
            END DO
	  END DO	
        END IF
      END DO
 111  FORMAT('MRO CloudFrac for ',F10.2,' cm-1; loop N/Tot',I3,I3,' swap numlays ',I3,' rad = ',F10.6)

      CALL wrtout(iIOUN,caOutName,raFreq,raWeightedRadiance)
      
      RETURN
      END
      
c************************************************************************
