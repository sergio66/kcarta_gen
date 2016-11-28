c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with K.F. Evan's rtspec code *************
c************** Any additional routines are also included here **********
c************************************************************************
c note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
c              rtspec, layer 1 == TOA, layer kProfLayer = ground
c                      there are nlev = 1 + iNumlayer  levels
c                      need to set temperature at levels from 1 .. 1+iNumLayer
c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output
      SUBROUTINE doscatter_rtspec(raFreq,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers,iaCldTypes(kmaxclouds)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSUrfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caOutName
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
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL rAngle
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)

      INTEGER i1,i2,iFloor,iDownWard

      DO i1=1,kMaxPts
        raInten(i1)=0.0
        ENDDO
     
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

      CALL interface_rtspec(raFreq,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 
 
      RETURN
      END

c************************************************************************
c the interface call to RTSPEC
c the main difference here is if the cloud layers for 2 different clouds are 
c noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
c artificial cloud of IWP=0.0g/m2 is set for layers 3-40
c kinda based on the interface to DISORT, except that it sets up this 
c intermediate "empty" cloud

      SUBROUTINE interface_rtspec(        
        !first the usual kCARTA variables
     $        raFreq,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
         !then the necessary scattering variables
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
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
c iDownward = +1 ==> downward looking instrument
c             -1 ==> upward looking instrument
    
      INTEGER iaCldTypes(kMaxClouds)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iDownWard,iBinaryFile
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers

c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)

c local variables
ccc      IMPLICIT NONE
C     The scattering tables are read in with READ_SSCATTAB.  The scattering
C     table is 3D: wavenumber, particle size, and viewing angle.
C         Scattering table variables:
C       MUTAB is view angle values (cosine zenith),
C       DMETAB is particle size values (median mass diameter, micron),
C       WAVETAB is wavenumber values (cm^-1).
C       MUINC(2) are the mu values of the two incident angles
C       TABEXTINCT is extinction, TABSSALB is single scattering albedo,
C       TABASYM is the asymmetry parameter
C       TABPHI??? are phase function info for incident directions
     

ccc      INTEGER  MAXTAB, MAXGRID, MAXSCAT
ccc      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

C         Radiative transfer variables: 
      INTEGER NSCATTAB, NCLDLAY, NABSNU, NLEV
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      INTEGER I, JNU1, JNU2
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
      REAL    SFCTEMP, SFCEMIS
      REAL    RADOBS
      REAL    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
      REAL  ABSNU1, ABSNU2, ABSDELNU
      REAL  WAVENO
      CHARACTER*80 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)
c      CHARACTER*24  OUTUNITS, OUTAVERAGING

c new local variables
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER iCldTopkCarta,iCldBotKcarta

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iBdry,FindBoundary,iaRadLayer(kProfLayer)

      INTEGER iCloudySky,iLayers,iII,iiDiv
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle,ttorad

C from original v2 rtspec code, to show differences :   
c      INTEGER  MAXABSNU, MAXNZ, MAXSPEC 
c      PARAMETER (MAXABSNU=100000, MAXNZ=100, MAXSPEC=10)  
c      INTEGER NUMSPEC, NSCATTAB, NCLDLAY, NABSNU, NLEV 
c      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ) 
c      INTEGER I, J, L, M, JNU1, JNU2 
cINTEGER J0, NBAND, IBAND, NWAVENO(MAXSPEC), NOUTNU 
cLOGICAL TWOSIDEBAND ********* not using outtb(maxspec),outtriang(maxspec) 
c      LOGICAL BINARYABSFILE 
c      REAL    ZOBS, MUOBS(MAXSPEC), ZTOP, IWP(MAXNZ), DME(MAXNZ) 
c      REAL    SFCTEMP, SFCEMIS 
cREAL SFCEMIS1(MAXSPEC), SFCEMIS2(MAXSPEC) 
cREAL RADOBS(2,MAXSPEC,MAXABSNU), OUTRAD(MAXABSNU) (real radobs,outrad,onu) 
c      REAL    HEIGHT(MAXNZ), TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU) 
c      REAL*8  ABSNU1, ABSNU2, ABSDELNU 
c      REAL*8  OUTNU1(MAXSPEC), OUTNU2(MAXSPEC), OUTDELNU(MAXSPEC) 
cREAL*8  BANDCENTER(MAXSPEC), BANDOFFSET(MAXSPEC) 
c****** used to have REAL*8  WAVENO, ONU1, ONU2, WT1, WT2, INVOUTDELNU 
c****** used to have REAL*8 sumwt1, sumwt2, sumrad1, sumrad2 
c      REAL*8  OUTNU(MAXABSNU), WAVENO(2,MAXSPEC,MAXABSNU), XO 
c      CHARACTER*100 SCATFILE(MAXSCAT), ABSFILE(MAXSPEC), OUTPUTFILE 
c      CHARACTER*1   RTMODEL, ABSTYPE 
cCHARACTER*24  OUTUNITS(MAXSPEC), OUTAVERAGE(MAXSPEC) (was outunits,outavging) 
      REAL    RAD0UPOBS(MAXNZ), RAD0DNOBS(MAXNZ)

      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
       
      iIOUN = kStdkCarta

c      WRITE (*,*) 'Singlescat, Eddington or hybrid RT model (S, E or H)'
      IF (kScatter .EQ. 1) THEN
        RTMODEL='S'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, SingleScat'
      ELSE IF (kScatter .EQ. 2) THEN
        RTMODEL='E'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Eddington'
      ELSE IF (kScatter .EQ. 3) THEN
        RTMODEL='H'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Hybrid'
        END IF

c      WRITE(*,*) 'Number of scattering tables'
c      NSCATTAB=iNclouds
c      IF (NSCATTAB .GT. MAXSCAT) THEN
c         write(kStdErr,*) 'RTSPEC: MAXSCAT exceeded'
c         CALL DoSTOP 
c         END IF
c      WRITE(*,*) 'Input scattering table file name for each table'
c      DO I = 1, NSCATTAB
c        READ (*,'(A100)') SCATFILE(I)
c      ENDDO

c      WRITE (*,*) 'Absorption file type (A-ascii text, B-binary)'
c      READ (*,'(A1)') ABSTYPE 
c      BINARYABSFILE = ABSTYPE .EQ. 'B' .OR. ABSTYPE .EQ. 'b'
c      WRITE(*,*) 'Output file name'
c      READ (*,'(A)') OUTPUTFILE
C         Read in the first absorption profile file
C           Absorption file may be binary or ascii
c      CALL READ_ABS_PROFILE (ABSFILE(1), BINARYABSFILE,
c     $               ABSNU1, ABSNU2, ABSDELNU, NABSNU, 
c     $               MAXNZ, MAXABSNU, NLEV, HEIGHT, TEMP, ABSPROF)
c      IF (OUTDELNU(1) .LT. ABSDELNU) THEN
c        PRINT *, OUTDELNU(1),ABSDELNU
c        WRITE (kStdErr,*) 'RTSPEC: output wavenumber resolution less than',
c     $                ' in absorption file'
c        CAll DoSTOP
c      ENDIF

      CALL SetMieTables_RTSPEC(raFreq, 
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, 
     $        -1,             !!!!iSergio = -1 as this is ORIG RTSPEC code 
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT,  iCldTopKcarta,iCldBotkCarta)

      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

        IF (iDownWard .GT. 0) THEN    
          !down look instr, in kCARTA layering style
          IOBS   = iNumLayer     
        ELSE IF (iDownWard .LT. 0) THEN    !up look instr
          !up look instr, in KCARTA layering style
          IOBS   = 1             
          END IF
      !change to RTSPEC layering
      iobs=(iNumLayer+1)-iobs+1

      !!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
      IF (iCloudySky .LT. 0) THEN
        rSolarAngle = kSolarAngle
        !!!!note that we do not care about background thermal accurately here
        write (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
        DO iii = 1,iNp
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
              END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $         rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,
     $         raUseEmissivity(iF),raFreq(iF),raInten(iF),iaOp(iii),
     $         iProfileLayers,raPressLevels)
            END DO
           CALL wrtout(iIOUN,caOutName,raFreq,raInten) 
           END DO
        GOTO 9876
        END IF

C     Get spectral range in absorption file for this spectral segment
      JNU1=1
      JNU2 = kMaxPts

      MUOBS=cos(rSatAngle*kPi/180.0)  !viewing angle cosine
      IF (iDownWard .EQ. -1) THEN
        MUOBS=-MUOBS
        END IF
      SFCTEMP = rSurfaceTemp
      WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',muobs,sfctemp

c see if there are any layers between TOA and instr; if there are, set
c TOA_to_instr to the cumulative radiative transfer of the necessary k, else 
c set it to -10.0
      DO iL=1,kProfLayer
        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
        END DO
      iiDiv = 0
 555  CONTINUE
      IF (iiDiv*kProfLayer .LT. iaRadLayer(1)) THEN
        iiDiv = iiDiv + 1
        GOTO 555
        END IF
      iiDiv = iiDiv - 1

      iBdry=-1
      !find the layer where you quit doing background thermal using only 
      !acos(3/5) and start using more accurate approx
      IF (kThermal .GE. 0) THEN
        iBdry=FindBoundary(raFreq,iProfileLayers,raPresslevels,iaRadLayer) 
        END IF

      CALL Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raFreq,raaAbs,TOA_to_instr) 

      DO iii = 1,iNp                     !!!!!!loop over outputs iii=1,iNp
        IF (iDownWard .EQ. +1)  THEN
          iobs = iaOp(iii)+1
          iobs = (iNumLayer+1) - iobs + 1
        ELSE 
          iobs = iaOp(iii)
          END IF
        iObs = iObs + (iiDiv*kProfLayer)

        IF (ICLDTOP.EQ.0 .OR. ICLDBOT.EQ.0 .OR. IOBS.EQ.0) THEN
          WRITE (kStdErr,*) 'RTSPEC: Observer or cloud top height does',
     .                ' not match absorption file levels'
          STOP
          END IF

c        print *,rtmodel,muobs,iobs,nlev,sfctemp,ncldlay,icldtop,icldbot
c        DO iF = JNU1, JNU2            !!!!loop over frequencies
c          raInten(iF) = 0.0
c          END DO
c        jnu1 = 3080
c        jnu2 = 3080

        DO iF = JNU1, JNU2            !!!!loop over frequencies
          IF ((iObs .EQ. 1) .AND. (iDownWard .LT. 0)) THEN
            raInten(iF) = ttorad(raFreq(iF),sngl(kTSpace))
          ELSE
            DO iL = 1,NLEV-1
              taugas(iL) = absprof(iL,iF)
              END DO

            SFCEMIS=raUseEmissivity(iF)
            WAVENO = raFreq(iF)
            CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL,
     $               MUOBS, IOBS, WAVENO,
     $               NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, 
     !!!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID, 
     $               NSCATTAB, MUINC,
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,
     $               TABEXTINCT, TABSSALB, TABASYM,
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN,
     $               RADOBS, 
     $               IBDRY, TOA_to_instr(iF),iiDiv) !these are 3 new parameters
            raInten(iF)=radobs     

            END IF                  !!!!end actually computing rad if things ok
          ENDDO                     !!!!loop over frequencies
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END DO

 9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm

      RETURN
      END

c************************************************************************
c this function quickly computes the background thermal contribution downto
c the layer where the instrument is
      REAL FUNCTION backgnd(raaAbs,waveno,raVTemp,iaaRadLayer,iNumLayer,
     $                      iNpmix,iAtm,rFracTop,j)

      include '../INCLUDE/scatter.param'

      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,raVTemp(kMixFilRows),waveno
      INTEGER iAtm,j,iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iNpmix
c     raFreq(j)=waveno

      INTEGER iLay,iaRadLayer(kProfLayer),iaRadLayerTemp(kProfLayer)
      INTEGER iT,iExtraThermal
      REAL rad,ttorad,rMPTemp,raExtraThermal(kMaxPts),rThetaEff
      REAL rTrans, rEmiss, rAbs

      DO iLay=1,iNumLayer 
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) 
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

      CALL AddUppermostLayersQ(iaRadLayer,iNumLayer,rFracTop, 
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal) 

      rThetaEff=3/5.0

c see  SUBROUTINE Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace, 

c do radiation at space temp
      rMPTemp = sngl(kTSpace)
      rad = ttorad(waveno,rMPTemp)

c go from top of atmosphere to instrument 
      IF (iExtraThermal .LT. 0) THEN
c just have to worry about top layer being fractional 
        rAbs=raaAbs(j,iaRadLayer(iNumLayer))*(1-rFracTop)
        rTrans=exp(-rAbs/rThetaEff)

        rMPTemp=raVTemp(iaRadLayer(iNumLayer))
        rEmiss=ttorad(waveno,rMPTemp)
        rEmiss=(1.0-rTrans)*rEmiss

        rad=rad*rTrans + rEmiss
      ELSE    
c go from top of atmosphere to layer above instrument 
        DO iLay=iT,iNumLayer+1
          rAbs=raaAbs(j,iaRadLayerTemp(iLay))
          rTrans=exp(-rAbs/rThetaEff)
          rMPTemp=raVTemp(iaRadLayerTemp(iNumLayer))
          rEmiss=ttorad(waveno,rMPTemp)	  
          rEmiss=(1.0-rTrans)*rEmiss
          rad=rad*rTrans + rEmiss
          END DO
c do layer where instrument is .. don't worry about interpolating temp
        DO iLay=iNumLayer,iNumLayer
          rAbs=raaAbs(j,iaRadLayerTemp(iLay))*(1-rFracTop)
          rTrans=exp(-rAbs/rThetaEff)
          rMPTemp=raVTemp(iaRadLayerTemp(iNumLayer))
          rEmiss=ttorad(waveno,rMPTemp)	  	  
          rEmiss=(1.0-rTrans)*rEmiss
          rad=rad*rTrans + rEmiss
          END DO
        END IF
 
      backgnd=rad*2.0*kPi

      RETURN
      END
c************************************************************************
c include Evans code here
cccccc      include '../INCLUDE/scatter_rtspec.f'   this is compiled separately

C      PROGRAM RTSPEC   NEW!!!!!!!!!!!!!!! for up and down look!!!!!!
C       Program for computing thermal atmospheric radiative transfer with 
C     a multiple layer cloud in an atmosphere.  Outputs a radiance spectrum,
C     suitably averaged.  Three approximate solutions for the scattering 
C     radiative transfer in the cloud layers are available: Single scattering 
C     approximation, Eddington second approximation and a hybrid single 
C     scattering/Eddington approximation.  Comparison with an "exact" 
C     doubling-adding model for ice clouds has shown that the hybrid model 
C     is superior in the mid IR (300-3000 cm^-1) but not always in the 
C     submm (0-50 cm^-1) (for particle sizes Dme=30 to 300 um).
C     Reference:  Deeter, M. N. and K. F. Evans, 1998: 
C         A hybrid Eddington-single scattering radiative transfer model 
C         for computing radiances from thermally emitting atmospheres. 
C         J. Quant. Spectosc. Radiat. Transfer, 60, 635-648.
C       The scattering approximations do radiative transfer for a single
C     homogeneous layer.  The input to them is the incident radiation on
C     the layer.  The incident radiance is computed by a nonscattering
C     radiative transfer integration of cloud and gaseous absorption profiles
C     that are read in from files.  The scattering routines output
C     upwelling radiance at the top of the layer, and this is propagated
C     to the observation level for output.  The multiple layers are dealt
C     with by computing the appropriate incident radiances/fluxes for each
C     layer, and then running the single layer routines for each layer from
C     the bottom up. For upward looking geometry, the downwelling radiance 
C     is computed by just switching the incident radiation and the top/bottom
C     temperature of the layer, thus the "upwelling" radiance that subroutine
C     outputs is actually the downwelling. Note: It is not necesary to switch
C     the "phi" functions since "+/-" doesn't mean up/down just forward/back.
C     For multiple layers, the single layer routines must be run from top down.
C       The cloud boundaries must be at levels in the input absorption profile.
C     The surface is assumed to reflect specularly with the specified 
C     emissivity. Multiple reflection between the cloud (scattering layer)
C     and surface is ignored.  The input includes the cloud properties,
C     which are ice water path (IWP) and particle size (Dme) for each layer.
C     There may be multiple input scattering table files for different 
C     cloud types.
C       The output is upwelling or downwelling brightness temperature or 
C     radiance as a function of wavenumber at the desired observation level.
C     The radiance is calculated at the absorption file wavenumbers.
C     For output the radiance may be averaged with either a triangular
C     or rectangular bandpass or convolved with sinc or sinc squared.
C     The width of the rectangular bandpass or the full width half max
C     (FHHM) of the triangular bandpass is equal to the output wavenumber 
C     spacing.  The output wavenumbers, which are the center of the 
C     bandpasses, are on the regular grid specified by the user.  The 
C     convolution averaging (for a Fourier Transform Spectrometer) is
C     done by FFTing and filtering with either a rectangular (sinc) or
C     triangular (sinc squared) windowing function. 
C       Multiple spectral segments may be output, with potentially different 
C     absorption profile input files, output units, and output bandpass 
C     averaging for each one.  Multiple "channels" may be output by making 
C     the output spacing for each segment equal to zero, setting the 
C     channel width, and using rectangular averaging. The channel width
C     extends equally from each side of the output wavenumber. If the channel 
C     width is one absorption file wavenumber spacing or less then 
C     monochromatic radiative transfer is done.  A double sideband receiver
C     may be emulated, in which case two rectangular bandpasses equally
C     spaced from the receiver center frequency are simulated.
C
C       The input ascii absorption file has the following format:
C     1) header line
C     2) starting, ending, delta, and number of wavenumbers in file
C     3) number of levels (1 + number of atmosphere layers)
C     4) heights of levels (from top to bottom) (in km)
C     5) temperatures at levels (in Kelvin)
C     wavenumber (cm^-1) and layer optical depths for each Nlev-1 layers
C        . . . 
C     A binary format absorption file may be input for faster loading
C     (the BINARYABSFILE logical is set from user input; see 
C     READ_ABS_PROFILE subroutine).
C
C       The input scattering tables are in the format produced by 
C     sscatmie.f.  This format also accomodates single scattering
C     information for randomly nonspherical oriented particles.
C
C     Author: Frank Evans  University of Colorado   May 1999
C     Co-authors:
C       Aaron Evans    (multilayer code, spectral segments, double sidebands,
C                       upward looking geometry, output convolution)
C       Merritt Deeter (single scattering layer routines)

c************************************************************************
