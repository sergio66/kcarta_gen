c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c (1)JacobGasAmtFM1,JacobTempFM1 : jacobians from the forward model
c    (includes solar contribution and thermal diffusive contribution)
c (2)Surface Reflectivity = 1/pi for thermal
c    Surface Reflectance for solar is defined by user

c the following variables are not size kMaxPtsJac or kProfLayerJac as they
c are well defined in the other non Jacobian routines
c raWaves(kMaxPts),raUseEmissivity(kMaxPts),raVTemp(kMixFilRows),
c iaaRadLayer(kMaxAtm,kProfLayer),raaAbs(kMaxPts,kMixFilRows)
c     $              raSurface,raSun,raThermal,raInten,
c     $              raSunRefl,
c     $              raLayAngles,raSunAngles)

c we also allow the user to compute the temperature jacobians in 
c one of three ways
c recall r(v)= sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT(planck)
c so          -1      ==> only use d/dT(1-tau)(tauL2S)
c so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

c************************************************************************
c**************************** GENERIC ROUTINES **************************
c************************************************************************

c this is the main driver subroutine for Jacobians, for SIMPLE SCATTERING
c it first redoes raaAbsTemp so that the clouds are included, and then 
c calls the main jacobian routines
c for the current frequency block, this subroutine calculates ALL the 
c jacobians and then outputs them
      SUBROUTINE find_jacobians_scat(raWaves,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,iNpmix)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c rDelta is the kComp file Step
c raLayAngles are the layer dependent satellite view angles
c raSunAngles are the layer dependent sun view angles
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which kcomp cm-1 block being output
c iNumGases is the number of gases to include
c iaGases is the integer array of GasID's
c iNumLayer is the number of layers in the atmosphere # iAtm
c iaaRadLayer is the list of radiating mixed paths to be included
c raVTemp are the layer temperatures

c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
c raaAllDT has the cumulative d/dT coeffs for current freq block
C     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
C        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
C        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

c raaSumAbCoeff is the cumulative absorption coeffs
c raaAmt  has the gas profiles
c raInten has the radiance vector
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
c raaMix is the mixing table

c these are to do with the arbitrary pressure layering
      REAL raThickness(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c FracTop,rFracBot are the upper layer/lower layer fractions
      REAL raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
      REAL raSurFace(kMaxPts),raSun(kMaxPts)
      REAL raThermal(kMaxPts),rDelta
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $      raVTemp(kMixFilRows),rSatAngle,raWaves(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      INTEGER iJacob,iaJacob(kMaxDQ)
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
      CHARACTER*80 caJacobFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL rAngle
      INTEGER iTag,iBinaryFile,iNpmix

c local variables
      REAL raaAbsTemp(kMaxPts,kMixFilRows)
      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

      INTEGER iDownWard,i1,i2,iFloor
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iaRadLayer(kProfLayer)

      INTEGER iCloudySky,iLayers,iII
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
      REAL raSolarBeam(kMaxPts),rSolarAngle,ttorad
       
      INTEGER NSCATTAB, NCLDLAY, NABSNU, NLEV
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

      DO iF=1,kMaxPts
        raSolarBeam(iF) = 0.0
        END DO

c these are the first few parts of simple_scat

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
          write(kStdErr,*) 'blackbody temp of space >> 2.96K'
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
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
        END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
        END IF

      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDownward .EQ. 1) THEN
        rAngle=rSatAngle
      ELSE
        rAngle=-rSatAngle
        END IF

c now these are the first few lines of interface_simple
      WRITE (kStdWarn,*) 'SIMPLE SCATTER radiative transfer code'

      CALL SetMieTables_RTSPEC(            
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,  
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        +1,             !!!!iSergio = +1 as this is MY code 
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 

      CALL CopyRaaAbs(raaAbs,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer)

      CALL AddCloud(raWaves,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

      CALL find_surface_backgnd_radiances(raWaves,raaAbsTemp,raVTemp, 
     $         iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,
     $         rTSpace,rTSurface,raUseEmissivity, 
     $         iProfileLayers,raPressLevels, 
     $         raSurface,raThermal) 

      CALLfind_jacobians(raWaves,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbsTemp,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,iProfileLayers,pProf)
      RETURN
      END

c************************************************************************

c this is the main driver subroutine for Jacobians
c for the current frequency block, this subroutine calculates ALL the 
c jacobians and then outputs them
      SUBROUTINE find_jacobians(raWaves,
     $            iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $            rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,iProfileLayers,pProf)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rDelta is the kComp file Step
c raLayAngles are the layer dependent satellite view angles
c raSunAngles are the layer dependent sun view angles
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which kcomp cm-1 block being output
c iNumGases is the number of gases to include
c iaGases is the integer array of GasID's
c iNumLayer is the number of layers in the atmosphere # iAtm
c iaaRadLayer is the list of radiating mixed paths to be included
c raVTemp are the layer temperatures

c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
c raaAllDT has the cumulative d/dT coeffs for current freq block
C     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
C        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
C        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

c raaSumAbCoeff is the cumulative absorption coeffs
c raaAmt  has the gas profiles
c raInten has the radiance vector
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
c raaMix is the mixing table

c these are to do with the arbitrary pressure layers
      REAL raPresslevels(kProfLayer+1),raThickness(kProfLayer),pProf(kProfLayer)
      INTEGER iProfileLayers
c FracTop,rFracBot are the upper layer/lower layer fractions
      REAL raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
      REAL raSurFace(kMaxPts),raSun(kMaxPts)
      REAL raThermal(kMaxPts),rDelta
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $      raVTemp(kMixFilRows),rSatAngle,raWaves(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      INTEGER iJacob,iaJacob(kMaxDQ)
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
      CHARACTER*80 caJacobFile

c local variables
      INTEGER iDownWard

c set the direction of radiation travel --- the checks of iUpper.iLower have
c already been done in radiance.f
c radiation travelling upwards to instrument ==> sat looking down iDownWard = 1
c radiation travelling down to instrument ==> sat looking up iDownWard =-1
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
        iDownWard = 1
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,iNumLayer))THEN
        iDownWard = -1
        END IF
      write(kStdWarn,*) 'in Jacobian, have set iDownWard = ',iDownWard

      IF (iDownWard .EQ. 1) THEN
        CALL DownWardJacobian(raWaves,iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,raSunRefl,rDelta)
      ELSE IF (iDownWard .EQ. -1) THEN
        CALL UpWardJacobian(raWaves,iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,rDelta)
        END IF

      RETURN
      END

c************************************************************************
c this subroutine multiplies the array by -1.0*constant where constant 
c depends on whether we are doing d/dT or d/dq
      SUBROUTINE MinusOne(raTorQ,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raResults is the array
c raTorQ === relevant element of raaaDq or raaDt
      REAL raTorQ(kMaxPtsJac)
      REAL raResults(kMaxPtsJac)

      INTEGER iFr
 
      DO iFr=1,kMaxPts
        raResults(iFr)=-raResults(iFr)*raTorQ(iFr)
        END DO
 
      RETURN
      END

c************************************************************************
c cumulatively, using contributions from each gas, find d/dT jacobian 
c for each layer
      SUBROUTINE cumulativeDT(daaDT,raaAllDT,raaMix,iG,iNatm,
     $                        iaaRadLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaDT has the current gas d/dT coeffs for current freq block
c raaAllDT has the cumulative d/dT coeffs for current freq block
c iNatm is the number of atmospheres to do radiance calcs for
c iG is the current gas
c raaMix is the mixing table
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      INTEGER iG,iNatm,iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore)

      INTEGER iL,iFr
      REAL rW

      IF (iNatm .GT. 1) THEN  
c cannot correctly weight the d/dT, so just use unit weight here and then try
c an average weight when JacobTemp is actually called
        write(kStdWarn,*)'Gas iG, weight rW = ',iG,1.0
        DO iL=1,kProfLayerJac
          DO iFr=1,kMaxPtsJac
            raaAllDT(iFr,iL)=raaAllDT(iFr,iL)+daaDT(iFr,iL)
            END DO
          END DO
      ELSE IF (iNatm .EQ. 1) THEN  
c have only one atmosphere and so correctly weight this gas's contribution to
c d/dT matrix ... then use weight of 1.0 when calling JacobTemp
        iL=iaaRadLayer(1,2)  !for atm#1, find which is the second mixed path
                             !as the first,last could have fractional weights
        rW=raaMix(iL,iG)     !find the gas weight in the second radiating layer
        write(kStdWarn,*)'Gas iG, weight rW = ',iG,rW
        DO iL=1,kProfLayerJac
          DO iFr=1,kMaxPtsJac
            raaAllDT(iFr,iL)=raaAllDT(iFr,iL)+rW*daaDT(iFr,iL)
            END DO
          END DO
        END IF
   
      RETURN
      END

c************************************************************************
c this subroutine does d/dr(tau_layer2space) for gas iG
c where r == gas amount q or temperature T at layer iM
c and  iL is the relevant layer we want tau_layer2space differentiated
c HENCE IF iL > iM, derivative == 0
c i.e. this does d(tau(l--> inf)/dr_m
      SUBROUTINE JacobTerm(iL,iM,raaLay2Sp,raTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Sp is the transmission frm layer to space
c iM has the layer that we differentiate wrt to
c iL has the radiating layer number (1..kProfLayerJac)
c raTemp has the results, apart from the multiplicative constants
c   which are corrected in MinusOne
      INTEGER iL,iM
      REAL raTemp(kMaxPtsJac),raaLay2Sp(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr

      IF (iL .GT. iM) THEN 
        DO iFr=1,kMaxPts
          raTemp(iFr)=0.0
          END DO
      ELSE
        DO iFr=1,kMaxPts
          raTemp(iFr)=raaLay2Sp(iFr,iL)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine does d/dr(tau_layer2gnd) for gas iG
c where r == gas amount q or temperature T at layer iM
c and  iL is the relevant layer we want tau_layer2space differentiated
c HENCE IF iL < iM, derivative == 0
c i.e. this does d(tau(l--> 0)/dr_m
      SUBROUTINE JacobTermGnd(iL,iM,raaLay2Gnd,raTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Gnd is the transmission frm layer to ground at diffusion angle
c iM has the layer that we differentiate wrt to
c iL has the radiating layer number (1..kProfLayerJac)
c raTemp has the results, apart from the multiplicative constants
c   which are corrected in MinusOne
      INTEGER iL,iM
      REAL raTemp(kMaxPtsJac),raaLay2Gnd(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iFr

      IF (iL .LT. iM) THEN 
        DO iFr=1,kMaxPts
          raTemp(iFr)=0.0
          END DO
      ELSE
        DO iFr=1,kMaxPts
          raTemp(iFr)=raaLay2Gnd(iFr,iL)
          END DO
        END IF

      RETURN
      END


c************************************************************************
c************** THESE HAVE TO DO WITH THE OUTPUT STYLE ******************
c************************************************************************
c this subroutine computes d(Brightness Temp)/d(Rad)
      SUBROUTINE Find_BT_rad(raInten,radBTdr,raWaves,
     $                       radBackgndThermdT,radSolardT)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c raInten is the radiance intensity at the instrument
c raWaves are the frequencies
c radBTdr is the derivative result
      REAL raWaves(kMaxPts),raInten(kMaxPts),radBTdr(kMaxPtsJac)
      REAL radBackgndThermdT(kMaxPtsJac),radSolardT(kMaxPtsJac)
      
      INTEGER iFr
      REAL r1,r2,r3,r4

      r1=kPlanck1
      r2=kPlanck2

      DO iFr=1,kMaxPts
        r3=r1*r2*(raWaves(iFr)**4)/(raInten(iFr)**2)
        r4=1.0+r1*(raWaves(iFr)**3)/raInten(iFr)
        radBTdr(iFr)=r3/r4/(alog(r4)**2)
        END DO

      IF (kThermal .LT. 0) THEN
        DO iFr=1,kMaxPts
          radBackgndThermdT(iFr)=0.0
          END DO
      ELSE
        DO iFr=1,kMaxPts
          r3=r1*r2*(raWaves(iFr)**4)/(radBackgndThermdT(iFr)**2)
          r4=1.0+r1*(raWaves(iFr)**3)/radBackGndThermdT(iFr)
          radBackgndThermdT(iFr)=r3/r4/(alog(r4)**2)
          END DO
        END IF

      IF (kSolar .LT. 0) THEN
        DO iFr=1,kMaxPts
          radSolardT(iFr)=0.0
          END DO
      ELSE
        DO iFr=1,kMaxPts
          r3=r1*r2*(raWaves(iFr)**4)/(radSolardT(iFr)**2)
          r4=1.0+r1*(raWaves(iFr)**3)/radSolardT(iFr)
          radSolardT(iFr)=r3/r4/(alog(r4)**2)
          END DO
        END IF

      RETURN
      END
c************************************************************************
c this subroutine prepares the output Jacobians according to kJacobOutput
      SUBROUTINE doJacobOutput(iLowest,raWaves,raResults,
     $                   radBTdr,raaAmt,raInten,iGasID,iM,iGasPosn)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLowest is the lowest layer in the atmosphere (modulo kProfLayer)
c raWaves are the frequency wavenumbers
c raResults are the raw d(rad)/d(gas amt)
c raaAmt are the gas profiles
c raInten is the radiant intensity at instrument
c iGasID is well duh... actually if iGasID = -1, then we are doing d/dT
c iM is the layer number (1..100)
c iGasPosn is the position of gasID in the gaslist
c radBTdr is the d(brightness temp)/d(Radiance) array
      INTEGER iGasID,iM,iLowest,iGasPosn
      REAL raWaves(kMaxPts),raResults(kMaxPtsJac)
      REAL raInten(kMaxPts)
      REAL raaAmt(kProfLayerJac,kGasStore),radBTdr(kMaxPtsJac)

      INTEGER iFr,iM1,iGasID_Temp

      iGasID_Temp = iGasID
      IF ((iGasID_temp .EQ. 101) .OR. (iGasID_temp .EQ. 102)) THEN
        iGasPosn  = 1
        END IF

      iM1=(iLowest-1)+iM

c      IF (kJacobOutPut .EQ. -1) THEN
c basically do nothing! user wants d(rad)/dq
c        END IF

      IF (kJacobOutput .NE. -1) THEN !oh well, do this

        IF (iGasID .GT. 0) THEN
c we are doing d/dq   
          IF (kJacobOutPut .EQ. 0) THEN
c user wants d(rad)/dq * q
            DO iFr=1,kMaxPts
              raResults(iFr)=raResults(iFr)*raaAmt(iM1,iGasPosn)
              END DO
          ELSE IF (kJacobOutPut .EQ. 1) THEN
c user wants d(BT)/dq * q
            DO iFr=1,kMaxPts
              raResults(iFr)=raResults(iFr)*raaAmt(iM1,iGasPosn)*radBTdr(iFr)
              END DO
            END IF

        ELSE IF (iGasID .LE. 0) THEN
c we are doing d/dT   
          IF (kJacobOutPut .EQ. 0) THEN
            iFr = 1 
c user wants d(rad)/dT so do nothing 
          ELSE IF (kJacobOutPut .EQ. 1) THEN
c user wants d(BT)/dT
            DO iFr=1,kMaxPts
              raResults(iFr)=raResults(iFr)*radBTdr(iFr)
              END DO
            END IF
          END IF

        END IF   !IF (kJacobOutput .NE. -1) THEN !oh well, do this

      RETURN
      END

c************************************************************************
c************ THESE HAVE TO DO WITH THE SURFACE PARAMETERS***************
c************************************************************************
c this subroutine does Jacobian wrt Surface Temperature
      SUBROUTINE JacobSurfaceTemp(raWaves,iM,
     $      rTSurface,raUseEmissivity,raaLay2Sp,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raaLay2Sp   is the layer-to-space abs coeff matrix
c raWaves has the frequencies
c raResults has the results
c iM are the layer <-> mixed path associations
      INTEGER iM
      REAL rTSurface,raUseEmissivity(kMaxPts)
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raWaves(kMaxPts)

c local variables
      REAL r1,r2,r3,r4,r5,rad,dradDT
      INTEGER iFr

      r1=kPlanck1
      r2=kPlanck2

      DO iFr=1,kMaxPts
        r3=r1*(raWaves(iFr)**3)
        r4=r2*raWaves(iFr)/rTSurface
        r5=exp(r4)
        rad=r3/(r5-1.0)
        dRadDT=rad*r4*r5/(r5-1.0)/rTSurface
        raResults(iFr)=dRadDT*raUseEmissivity(iFr)*raaLay2Sp(iFr,iM)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian wrt Surface Emissivity
      SUBROUTINE JacobSurfaceEmis(iM,raSurface,raThermal,raaLay2Sp,
     $              raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raSurface is the surface emission
c raaLay2Sp   is the layer-to-space abs coeff matrix
c raResults has the results
c raThermal has the downwelling thermal contribs
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac),raSurface(kMaxPts)
      REAL raResults(kMaxPtsJac),raThermal(kMaxPts)

c local variables
      INTEGER iFr

      DO iFr=1,kMaxPts
        raResults(iFr)=raaLay2Sp(iFr,iM)*
     $                 (raSurface(iFr)-raThermal(iFr)/kPi)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian of Backgnd Thermal wrt Surface Emissivity
      SUBROUTINE JacobBackgndThermal(iM,raaLay2Sp,raThermal,raResults) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c raaLay2Sp is the layer-to-space abscoeff matrix 
c raResults has the results 
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac)
      REAL raResults(kMaxPtsJac),raThermal(kMaxPts) 

c local variables
      INTEGER iFr

      DO iFr=1,kMaxPts
        raResults(iFr)=-raaLay2Sp(iFr,iM)/kPi*raThermal(iFr) 
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does Jacobian of Solar wrt Sun Surface Emissivit
      SUBROUTINE JacobSolar(iM,raaLay2Sp,raSun,raResults)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c raaLay2Sp is the layer-to-space abscoeff matrix 
c raResults has the results 
c raSun,raThermal has the downwelling solar,thermal contribs
c iM,rSatAngle are the layer <-> mixed path associations and satellite angle
      INTEGER iM
      REAL raaLay2Sp(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac),
     $     raSun(kMaxPts)

c local variables
      INTEGER iFr

c remember that raSun is that at the bottom of the atmosphere ==> have to 
c propagate to top of atmosphere
      DO iFr=1,kMaxPts
        raResults(iFr)=raSun(iFr)*raaLay2Sp(iFr,iM)
        END DO

      RETURN
      END

c************************************************************************








