c Copyright 2001
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
      SUBROUTINE find_radiances_perturb(
     $         raFreq,raaExt,raaScat,raaAsym,
     $         iPhase,raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
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
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),rSurfPress
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*120 caOutName
      REAL Temp(MAXNZ)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

      INTEGER i1,i2,iFloor,iDownWard

      !! --------- kAvgMin is a global variable in kcarta.param -------- !!
      !!kAvgMin is a global variable in kcarta.param .. set as required
      !!it is the average of single scattering albedo (w0); if less than some
      !!value, then basically there is no scattering and so can do some 
      !!approximations!!!!!
      kAvgMin = 1.0d-3     !!!before Feb 14, 2003
      kAvgMin = 1.0d-6
      !! --------- kAvgMin is a global variable in kcarta.param -------- !!

      DO i1=1,kMaxPts
        raInten(i1)=0.0
        ENDDO
     
c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2 = iaaRadLayer(iAtm,iNumLayer)-1
        i2 = iFloor(i2*1.0/kProfLayer)
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
        i1 = iaaRadLayer(iAtm,1)-1
        i1 = iFloor(i1*1.0/(1.0*kProfLayer))
        i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
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

      IF (iDownward .EQ. 1) THEN
        CALL rad_DOWN_perturb(raFreq,
     $        raInten,raVTemp,raaExt,raaScat,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

      ELSE 
        CALL rad_UP_perturb(raFreq,
     $        raInten,raVTemp,raaExt,raaScat,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp)
        END IF
 
      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

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

      SUBROUTINE rad_DOWN_perturb(raFreq,raInten,
     $    raVTemp,raaExt,raaScat,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

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
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*120 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv
      REAL raaLayTrans(kMaxPts,kProfLayer),rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),muSat
      REAL raInten1(kMaxPts),raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput

      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot

c do we need to output stuff from within the cloud?
      INTEGER iSimple,iLModKprofLayer,iSTopNormalRadTransfer,iF
      REAL raOutFrac(kProfLayer),r0 

      iIOUN = kStdkCarta

      rThermalRefl = 1.0/kPi

c calculate cos(SatAngle)
      muSat = cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
      IF (iDoSolar .GT. -1) THEN
        write(kStdErr,*)'No solar scatteriung in perturb code yet!'
        Call DoStop
        END IF

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
          iHigh = iaOp(iLay)
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
       DO iLay = iNumLayer,iNumLayer
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
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raInten1(iFr)  = 0.0
        raSurface(iFr) = raInten(iFr)
        END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL .LT. iNLTEStart) THEN   !normal, no LTE emission stuff
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iL .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
          END IF
        END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr) = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
        END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
      IF (iNLTEStart .GT. kProfLayer) THEN
        iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
        write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
      ELSE
        iLay = 1
 987    CONTINUE
        iL = iaRadLayer(iLay)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF ((iLModKprofLayer .LT. iNLTEStart).AND.(iLay .LT. iNumLayer)) THEN
          iLay = iLay + 1
          GOTO 987
          END IF
        iSTopNormalRadTransfer = iLay
        write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
        END IF

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer .LT. iNLTEStart) THEN   
          !normal, no LTE emission stuff
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iLModKprofLayer .GE. iNLTEStart) THEN
          !new; LTE emission stuff
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL) 
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
          END IF
        END DO

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
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
        END IF


c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c do the radiation at the surface
      DO iFr=1,kMaxPts
        raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $                 raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $                 raSun(iFr)*raSunRefl(iFr)
        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c now compute the upwelling radiation!!!!! at view angle upto cloud bottom
c DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

      IF (iLocalCldBot .GT. 1) THEN
c first do the bottommost layer (could be fractional) at viewing angle
        DO iLay=1,1
          iL = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            !only output at final TWOSTREAM pass from GND to BOT of CLD
            write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq,
     $          raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $          raSun,-1,iNumLayer,rFracTop,rFracBot,
     $          iProfileLayers,raPressLevels,
     $          iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
              END DO
            END IF
            
c now do the radiative transfer thru this bottom layer
          iSimple = +1
          iL = iaRadLayer(iLay) - iiDiv*kProfLayer
          CALL RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $           muSat,rFracBot,iL,iSimple,raInten1,raInten)
          END DO

c then do the layers till the cloudbot (all will be full) at the viewing angle
        DO iLay = 2,iLocalCldBot-1
          iL      = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            !only output at final TWOSTREAM pass from GND to BOT of CLD
            write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq,
     $          raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $          raSun,-1,iNumLayer,rFracTop,rFracBot,
     $          iProfileLayers,raPressLevels,
     $          iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
              END DO
            END IF
        
c now do the radiative transfer thru each of these complete layers
          iSimple = +1
          iL = iaRadLayer(iLay) - iiDiv*kProfLayer
          CALL RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $                       muSat,1.0,iL,iSimple,raInten1,raInten)
          END DO
        END IF

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        DO iLay = iLocalCldBot,iLocalCldTop
          iL      = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            !only output at final TWOSTREAM pass from GND to BOT of CLD
            write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq,
     $          raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $          raSun,-1,iNumLayer,rFracTop,rFracBot,
     $          iProfileLayers,raPressLevels,
     $          iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
              END DO
            END IF
        
c now do the radiative transfer thru each of these complete layers
        iSimple = -1
        iL = iaRadLayer(iLay) - iiDiv*kProfLayer
        IF (iLay .EQ. 1) THEN
          CALL RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $                       muSat,rFracBot,iL,iSimple,raInten1,raInten)
        ELSE
          CALL RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $                       muSat,1.0,iL,iSimple,raInten1,raInten)
          END IF
        END DO

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c then do the rest of the layers till the last but one (all will be full)
      DO iLay = iLocalCldTop+1,iHigh-1
         iL = iaRadLayer(iLay)
         muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'at iLay = ',iLay,' going up from past cloud'
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq,
     $        raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF

c now do the radiative transfer thru this complete layer
        iSimple = +1
        iL = iaRadLayer(iLay) - iiDiv*kProfLayer
        CALL RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,muSat,1.0,
     $                       iL,iSimple,raInten1,raInten)
        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            write(kStdWarn,*) 'assume you need to output rad at TOA'
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
            write(kStdWarn,*) 'and output everything at the top of this'
            write(kStdWarn,*) 'stratosphere'
            !do radiative transfer thru this layer
            DO iFr=1,kMaxPts
              raInten(iFr) = 
     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
              END DO
            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,-1)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
            DO iFr=1,iDp
              CALL wrtout(iIOUN,caOutName,raFreq,raInten)
              END DO
            END IF
          END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
           IF (iDp .GT. 0) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
             DO iFr=1,iDp
               CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
               CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
              END DO
            END IF
          END IF

cc no need to do radiative transfer thru this layer
cc        DO iFr=1,kMaxPts
cc          raInten(iFr) = raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc          END DO

        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this is for an UPLOOK instrument

c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis

      SUBROUTINE rad_UP_perturb(raFreq,raInten,
     $    raVTemp,raaExt,raaScat,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp)
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
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*120 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
      REAL raaLayTrans(kMaxPts,kProfLayer),rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot,iRepeat

c general coeffs for the layers
      REAL mu_view
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts)
      REAL radSolarCld(kMaxPts),raW0(kMaxPts),raAsym0(kMaxPts)

c arbitrary angle stuff
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)

c do we need to output stuff from within the cloud?
      REAL raTop(kMaxPts),raBot(kMaxPts)
      INTEGER iInsideCloud,iSimple

c other stuff for uplook inst
      REAL raDiffuseInten(kMaxPts),raDownViewAngle(kMaxPts)
      REAL rOmegaSun,rLocalAbs,rFrac,rBotOfCld
      INTEGER iLow,iiDiv

      write(kStdErr,*) 'No uplook code yet for scatter_perturb_uplook'
      CALL DoStop

      RETURN
      END 

c************************************************************************
c this subroutine does the zeroth plus first order radiation computation
c this uses P(u,u') = 1 + 3 g u u' when finding O(1) solution
      SUBROUTINE RT_perturb_UP(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $                         muSat,rFrac,iL,iSimple,raInten1,raInten)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raFreq(kMaxPts)
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL rFrac,muSat,raVT1(kMixFilRows)
      INTEGER iSimple                  !do only zeroth order (clear sky) +1
                                       !or zero+first order (scatter)    -1
      INTEGER iL                       !which kMixFilRow to use from raaBoink
c input/output vars
      REAL raInten(kMaxPts)            !Zeroth + First order intensity
      REAL raInten1(kMaxPts)           !First order intensity

c local vars
      INTEGER iFr
      REAL rBT,rX,rW,r0,r1,rA,rB,raTemp2Rad(kMaxPts)

      IF (iSimple .GT. 0) THEN
        !!!straightforward (non scattering!!) computation
        rBT = raVT1(iL)
        CALL ttorad_oneBT2array(raFreq,rBT,raTemp2Rad)
        DO iFr = 1,kMaxPts
          rX = exp(-raaExt(iFr,iL)/muSat)
          raInten(iFr) = raInten(iFr)*rX + raTemp2Rad(iFr)*(1.0-rX)
          END DO
        END IF

      IF (iSimple .LT. 0) THEN
        !!!do first  order first, then zeroth order
        rBT = raVT1(iL)
        CALL ttorad_oneBT2array(raFreq,rBT,raTemp2Rad)
        CALL Perturb_First(raFreq,rBT,iL,muSat,raaExt,raaScat,raaAsym,
     $                     raTemp2Rad,raInten,raInten1)
        DO iFr = 1,kMaxPts
          rX = exp(-raaExt(iFr,iL)/muSat)
          r0 = raInten(iFr)*rX + raTemp2Rad(iFr)*(1.0-rX) !!zeroth order
          raInten(iFr) = r0 + raInten1(iFr)
          END DO
        END IF

      RETURN
      END 

c************************************************************************
c this is the guts of the perturbation computation
      SUBROUTINE Perturb_First(raFreq,rBT,iL,muSat,raaExt,raaScat,raaAsym,
     $                         raTemp2Rad,raInten,raInten1)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raFreq(kMaxPts),rBT,muSat
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raInten(kMaxPts),raTemp2Rad(kMaxPts)
      INTEGER iL
c output vars
      REAL raInten1(kMaxPts)

c local vars
      INTEGER iFr
      REAL raW(kMaxPts),raX(kMaxPts),raG(kMaxPts),raEmiss(kMaxPts),r1
      REAL raBC1(kMaxPts),raBC2(kMaxPts),ra1(kMaxPts),ra2(kMaxPts)
      REAL raBC(kMaxPts),raBCE(kMaxPts),rMu

      rMu = muSat
      DO iFr = 1,kMaxPts
        raX(iFr) = raaExt(iFr,iL)
        raG(iFr) = raaAsym(iFr,iL)
        raW(iFr) = raaScat(iFr,iL)/raaExt(iFr,iL)
        END DO

      CALL FromIncidentIO(raX,raW,raG,rMu,ra1,raBC1)
      CALL FromRadiatingLayer(raX,raW,raG,rMu,ra1,ra2,raBC2)
      CALL FromEmission(rMu,raX,raW,raTemp2Rad,raEmiss,raBCE)
      CALL BoundaryCondition(rMu,raInten,raTemp2Rad,raX,raW,raG,raBC)

      DO iFr = 1,kMaxPts
c        r1 = raBC(iFr) + 
c     $       raW(iFr)/rMu*(ra1(iFr)*raInten(iFr) + ra2(iFr)*raTemp2Rad(iFr))
c asuume raBC = 0 always
        r1 = raW(iFr)/rMu*((ra1(iFr)-raBC1(iFr))*raInten(iFr) + 
     $                     (ra2(iFr)-raBC2(iFr))*raTemp2Rad(iFr)) +
     $                     (raEmiss(iFr)-raBCE(iFr))
        raInten1(iFr) =  r1*exp(-raX(iFr)/rMu)
        print *,iFr,raX(iFr),raW(iFr),raEmiss(iFr)-raBCE(iFr),
     $                                (ra1(iFr)-raBC(iFr))*raInten(iFr),
     $                                (ra2(iFr)-raBC2(iFr))*raTemp2Rad(iFr)
        END DO

      RETURN
      END 

c************************************************************************
c this function simply does the emission integration
c raBC is the boundary condition at tau = 0
      SUBROUTINE FromEmission(rMu,raX,raW,raTemp2Rad,raEmiss,raBCE)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raW(kMaxPts),raTemp2Rad(kMaxPts),raX(kMaxPts),rMu
c output vars
      REAL raEmiss(kMaxPts),raBCE(kMaxPts)

      INTEGER iFr

      DO iFr = 1,kMaxPts
        raEmiss(iFr) = -raTemp2Rad(iFr)*raW(iFr)*exp(-raX(iFr)/rMu)
        raBCE(iFr) =   +raTemp2Rad(iFr)*raW(iFr)
        END DO

      RETURN
      END 

c************************************************************************
c this function computes scattering of incident radiation into beam
c this is integral [exp(tau/mu))*(E2(tau) + 3 g mu E3(tau))]
c raBC is the boundary condition at tau = 0
      SUBROUTINE FromIncidentIO(raX,raW,raG,rMu,ra1,raBC1)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raW(kMaxPts),raX(kMaxPts),raG(kMaxPts),rMu
c output vars
      REAL ra1(kMaxPts),raBC1(kMaxPts)

c local vars
      REAL r1,r2
      INTEGER iFr

      CALL E2E3(raX,raW,raG,rMu,ra1)

      IF (abs(rMu-1) .GT. 1e-16) THEN
        r1 = 1/(2*(1-rMu)*(1-rMu))*rMu*(1 + (2*rMu-3)*rMu*rMu)
      ELSE
        r1 = 1.5
        END IF

      DO iFr = 1,kMaxPts 
        r2 = rMu + 3*raG(iFr)*rMu*r1
        raBC1(iFr) = r2
        END DO
      
      RETURN
      END 

c************************************************************************
c this function computes scattering of layer planck radiation into beam
c this is integral [exp(tau/mu))*(1 - E2(tau) - 3 g mu E3(tau))]
c so it basically uses results ra1 from subroutine FromIncidentIO
c raBC is the boundary condition at tau = 0
      SUBROUTINE FromRadiatingLayer(raX,raW,raG,rMu,ra1,ra2,raBC2)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raW(kMaxPts),raX(kMaxPts),raG(kMaxPts),rMu
      REAL ra1(kMaxPts)                   !this has E2E3 results
c output vars
      REAL ra2(kMaxPts),raBC2(kMaxPts)

c local vars
      REAL r1,r2
      INTEGER iFr

      DO iFr = 1,kMaxPts
        !!!need to do the extra integral exp(t/mu) dt; else same as ra1
        ra2(iFr) = rMu*exp(raX(iFr)/rMu) - ra1(iFr)
        END DO

      IF (abs(rMu-1) .GT. 1e-16) THEN
        r1 = 1/(2*(1-rMu)*(1-rMu))*rMu*(1 + (2*rMu-3)*rMu*rMu)
      ELSE
        r1 = 1.5
        END IF

      DO iFr = 1,kMaxPts 
        r2 = rMu + 3*raG(iFr)*rMu*r1
        raBC2(iFr) = 1-r2
        END DO

      RETURN
      END 

c************************************************************************
c this is integral [exp(tau/mu))*(E2(tau) + 3 g mu E3(tau))]
c and is used by subroutines FromRadiatingLayer,FromIncidentIO
c also does the evaluation of the BC at tau = 0
      SUBROUTINE E2E3(raX,raW,raG,rMu,ra1)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raW(kMaxPts),raX(kMaxPts),raG(kMaxPts),rMu
c output vars
      REAL ra1(kMaxPts)

c local vars
      INTEGER iFr
      REAL r2,r3,rMuM1,rRecip,raPower1(kMaxPts),raPower2(kMaxPts)

c integral [exp(tau/mu) E2(tau)] = int [exp(tau/mu)(exp(-tau) - tau E1(tau))
      IF (abs(rMu-1) .GE. 1.0e-5) THEN
        !view angle far from nadir
        rMuM1  = 1/rMu - 1
        rRecip = rMu/(1-rMu)
        CALL IntegPower1(1.0,1/rMu,raX,raPower1)
        CALL IntegPower2(1.0,1/rMu,raX,raPower2)
        DO iFr = 1,kMaxPts
          !!this is integral E2(t) exp(t/mu)
          r2 = rRecip*exp(raX(iFr)*rMuM1) - raPower1(iFr)
          !!this is integral E3(t) exp(t/mu)
          r3 = 0.5*(rRecip*exp(raX(iFr)*rMuM1)*(1 - (raX(iFr)-rRecip)) + 
     $         raPower2(iFr))
          ra1(iFr)  = r2 + r3*3*rMu*raG(iFr)
          END DO
      ELSEIF (abs(rMu-1) .LT. 1.0e-5) THEN
        !view angle far almost nadir
        CALL IntegPower1(1.0,1/rMu,raX,raPower1)
        CALL IntegPower2(1.0,1/rMu,raX,raPower2)
        DO iFr = 1,kMaxPts
          !!this is integral E2(t) exp(t/mu)
          r2 = raX(iFr) - raPower1(iFr)
          !!this is integral E3(t) exp(t/mu)
          r3 = 0.5*(raX(iFr) - raX(iFr)*raX(iFr)/2 + raPower2(iFr))
          ra1(iFr)  = r2 + r3*3*rMu*raG(iFr)
          END DO
        END IF

      RETURN
      END 

c************************************************************************
c this evalutes the integrals at z = tau = 0
      SUBROUTINE BoundaryCondition(rMu,raInten,raTemp2Rad,raX,raW,raG,raBC)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL rMu,raInten(kMaxPts),raTemp2Rad(kMaxPts)
      REAL raW(kMaxPts),raX(kMaxPts),raG(kMaxPts)
c output vars
      REAL raBC(kMaxPts)

c local vars
      INTEGER iFr
      REAL r1,r2,r3

      IF (abs(rMu-1) .GT. 1e-16) THEN
        r1 = 1/(2*(1-rMu)*(1-rMu))*rMu*(1 + (2*rMu-3)*rMu*rMu)
      ELSE
        r1 = 1.5
        END IF

      DO iFr = 1,kMaxPts 
        r2 = rMu + 3*raG(iFr)*rMu*r1
        r3 = raInten(iFr)*r2 + raTemp2Rad(iFr)*(1-r2)
        raBC(iFr) = raBC(iFr) + raTemp2Rad(iFr)*raW(iFr) - raW(iFr)/rMu*r3
        END DO

      RETURN
      END 

c************************************************************************
c this function does integral x exp(beta x) E1(gamma x)
c for our applications, gamma = 1, beta = 1/cos(theta) = 1 to infinity
      SUBROUTINE IntegPower1(rGamma0,rBeta0,raX,raPower1)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raX(kMaxPts),rGamma0,rBeta0
c output vars
      REAL raPower1(kMaxPts)

c local vars
      INTEGER iFr
      REAL rX,rY,rZ,rGB,rBG,raTau(kMaxPts),raY1(kMaxPts),raY2(kMaxPts)
      REAL rBeta,rGamma

      IF (abs(rBeta0-rGamma0) .LE. 1.0e-15) THEN
        rBeta = rBeta0 + 1.0e-5
        rGamma = rGamma0
      ELSE
        rBeta = rBeta0
        rGamma = rGamma0
        END IF

      rBG = rBeta - rGamma !!typically 0 or positive
      rGB = rGamma - rBeta !!typically 0 or negative

      !!!rBeta-rGamma > 0; since rGamma = 1, rBeta = 1/cos(theta) not nadir

      DO iFr = 1,kMaxPts
        raTau(iFr) = raX(iFr)*rGB
        END DO
      CALL expint(raTau,raY1)

      DO iFr = 1,kMaxPts
        raTau(iFr) = raX(iFr)*rGamma
        END DO
      CALL expint(raTau,raY2)

      DO iFr = 1,kMaxPts
        rX = rBeta*raX(iFr)
        rY = rGB*raX(iFr)
        rZ = -rBeta*exp(-rY) + rGB*raY1(iFr) + 
     $       rGB*exp(rX)*(rX-1)*raY2(iFr)
        raPower1(iFr) = rZ/rBeta/rBeta/rGB
        END DO

      RETURN
      END 

c************************************************************************
c this function does integral x2 exp(beta x) E1(gamma x)
c for our applications, gamma = 1, beta = 1/cos(theta) = 1 to infinity
      SUBROUTINE IntegPower2(rGamma0,rBeta0,raX,raPower2)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raX(kMaxPts),rGamma0,rBeta0
c output vars
      REAL raPower2(kMaxPts)

c local vars
      INTEGER iFr
      REAL rX,rBG,rGB,raTau(kMaxPts),raY1(kMaxPts),raY2(kMaxPts),r1,r2
      REAL rBeta,rGamma

      IF (abs(rBeta0-rGamma0) .LE. 1.0e-15) THEN
        rBeta = rBeta0 + 1.0e-5
        rGamma = rGamma0
      ELSE
        rBeta = rBeta0
        rGamma = rGamma0
        END IF

      rBG = rBeta - rGamma !!typically 0 or positive
      rGB = rGamma - rBeta !!typically 0 or negative

      !!!this is the actual integral
      !!!rBeta-rGamma > 0; since rGamma = 1, rBeta = 1/cos(theta) not nadir

      DO iFr = 1,kMaxPts
        raTau(iFr) = raX(iFr)*rGB
        END DO
      CALL expint(raTau,raY1)

      DO iFr = 1,kMaxPts
        raTau(iFr) = raX(iFr)*rGamma
        END DO
      CALL expint(raTau,raY2)

      DO iFr = 1,kMaxPts
        rX = rBeta*raX(iFr)
        r1 = rBeta*(rX-3) - rGamma*(rX-2)
        r2 = rBG*raX(iFr)
        raPower2(iFr) = -exp(rX)*(((1-rX)*(1-rX)+1)*raY2(iFr))  +
     $                  (-1)/(rBG*rBG)*(rBeta*exp(r2)*r1-2*raY1(iFr)*rBG*rBG)
        raPower2(iFr) = raPower2(iFr)/rBeta/rBeta/rBeta
        END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the zeroth plus first order radiation computation
c this is kinda a silly call as it assumes that P(u,u') = delta(u,u')
c ie there only is forward scattering ==> really, this is scatter free!
c and so w should be zero
      SUBROUTINE RT_perturb_UP_delta(raFreq,raaExt,raaScat,raaAsym,raVT1,
     $                         muSat,rFrac,iL,iSimple,raInten1,raInten)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c input vars
      REAL raFreq(kMaxPts)
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL rFrac,muSat,raVT1(kMixFilRows)
      INTEGER iSimple                  !do only zeroth order (clear sky) +1
                                       !or zero+first order (scatter)    -1
      INTEGER iL                       !which kMixFilRow to use from raaBoink
c input/output vars
      REAL raInten(kMaxPts)            !Zeroth + First order intensity
      REAL raInten1(kMaxPts)           !First order intensity

c local vars
      INTEGER iFr
      REAL rBT,rX,rW,rFirst,rSecond,raTemp2Rad(kMaxPts)

      IF (iSimple .GT. 0) THEN
        rW = 0.0
        rFirst = 0.0
        !!!straightforward computation
        rBT = raVT1(iL)
        CALL ttorad_oneBT2array(raFreq,rBT,raTemp2Rad)
        DO iFr = 1,kMaxPts
          rX = exp(-raaExt(iFr,iL)/muSat)
          raInten(iFr) = raInten(iFr)*rX + raTemp2Rad(iFr)*(1.0-rX)
          END DO
        END IF

      IF (iSimple .LT. 0) THEN
        !!!do first  order first, then zeroth order
        rBT = raVT1(iL)
        CALL ttorad_oneBT2array(raFreq,rBT,raTemp2Rad)
        DO iFr = 1,kMaxPts
          rX = exp(-raaExt(iFr,iL)/muSat)
          rW = (raaScat(iFr,iL)/raaExt(iFr,iL))/muSat
          rFirst       = rW*(raInten(iFr)-raTemp2Rad(iFr))*raaExt(iFr,iL)*rX
          rSecond      = rFirst * rW/muSat * raaExt(iFr,iL)/2.0
          raInten(iFr) = raInten(iFr)*rX + raTemp2Rad(iFr)*(1.0-rX) + 
     $                   rFirst + rSecond
          END DO
        END IF

      RETURN
      END 

c************************************************************************
