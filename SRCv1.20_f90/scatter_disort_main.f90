! Copyright 1997
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:44
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with STamnes et al  disort code  *********
!************** Any additional routines are also included here **********
!************************************************************************
!************** All the routines in this file are necessary *************
!************************************************************************

!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              disort, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!************************************************************************


! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

SUBROUTINE doscatter_disort(raFreq,raaAbs,raVTemp,  &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,iProfileLayers,pProf,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raNumberDe
INTEGER, INTENT(IN OUT)                  :: iDoFlux
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iDoFlux     = do radiance or flux computation
! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
! raNumberDensity = P/RT == number of particles in each layer of atm
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)

INTEGER :: iProfileLayers
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raNumberDensity(kProfLayer)
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp


INTEGER :: iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG


INTEGER :: i1,i2,iFloor,iDownWard

DO i1=1,kMaxPts
  raInten(i1) = 0.0
  ENDDO
    
! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer) = 100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = 1
      i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2=iaaRadLayer(iAtm,iNumLayer)-1
      i2=iFloor(i2*1.0/kProfLayer)
      IF (rTSpace > 5.0) THEN
        WRITE(kStdErr,*) 'you want satellite to be downward looking'
        WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        WRITE(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
! radiation travelling downwards to instrument ==> sat looking up
! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer) = 100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1=iaaRadLayer(iAtm,1)-1
      i1=iFloor(i1*1.0/(1.0*kProfLayer))
      i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard
    
! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
      WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
      CALL DoSTOP
    END IF
    
! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF
    
    WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end) = ',  &
        iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)
    
    IF (iDownward == 1) THEN
      rAngle=rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF
    
    CALL interface_disort(raFreq,raInten,raVTemp,  &
        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
        rAngle,rFracTop,rFracBot,
!     $        iNp,iaOp,raaOp,iNpmix,iFileID,  &
    iNp,iaOp,iNpmix,iFileID, caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
!     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
!     $        raLayAngles,raSunAngles,  &
    raLayAngles,iDoFlux, raThickness,raPressLevels,iProfileLayers,pProf,  &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)
    
    RETURN
  END SUBROUTINE doscatter_disort
  
!************************************************************************
! this subroutine sets up the scattering table info from SSCATMIE.F
  
  SUBROUTINE SetMieTables_DISORT(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
      raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT)
  
  
  REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
  INTEGER, INTENT(IN OUT)                  :: iAtm
  NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
  NO TYPE, INTENT(IN)                      :: iNclouds
  NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
  NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
  NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
  NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
  NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
  INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
  NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
  NO TYPE, INTENT(IN OUT)                  :: raComputed
  NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
  NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
  INTEGER, INTENT(IN)                      :: iNumLayer
  NO TYPE, INTENT(IN OUT)                  :: iDownWard
  NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
  INTEGER, INTENT(IN OUT)                  :: NMUOBS(MAXSCAT)
  INTEGER, INTENT(IN OUT)                  :: NDME(MAXSCAT)
  INTEGER, INTENT(IN OUT)                  :: NWAVETAB(MAXSCAT)
  REAL, INTENT(IN OUT)                     :: MUTAB(MAXGRID,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: DMETAB(MAXGRID,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: WAVETAB(MAXGRID,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: MUINC(2)
  REAL, INTENT(IN OUT)                     :: TABEXTINCT(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABSSALB(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABASYM(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABPHI1UP(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABPHI1DN(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABPHI2UP(MAXTAB,MAXSCAT)
  REAL, INTENT(IN OUT)                     :: TABPHI2DN(MAXTAB,MAXSCAT)
  INTEGER, INTENT(OUT)                     :: NSCATTAB
  INTEGER, INTENT(OUT)                     :: NCLDLAY
  INTEGER, INTENT(OUT)                     :: ICLDTOP
  INTEGER, INTENT(OUT)                     :: ICLDBOT
  INTEGER, INTENT(OUT)                     :: IOBS
  INTEGER, INTENT(OUT)                     :: ISCATTAB(MAXNZ)
  REAL, INTENT(OUT)                        :: IWP(MAXNZ)
  REAL, INTENT(OUT)                        :: DME(MAXNZ)
  NO TYPE, INTENT(IN OUT)                  :: iaCloudWit
  NO TYPE, INTENT(IN OUT)                  :: iaScatTabl
  INTEGER, INTENT(OUT)                     :: iCloudySky
  INTEGER, INTENT(OUT)                     :: IACLDTOP(kMaxClouds)
  INTEGER, INTENT(OUT)                     :: IACLDBOT(kMaxClouds)
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
! ---------------- inputs needed to read scattering tables -------------------
  
! this is which atm number is being used, and whether these are binary files
  INTEGER :: iBinaryFile, iDownward
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
  INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
  INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
  INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
  INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
  CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
  REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this is just to set everything about clouds relative to TOA layer
  INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
! this tells if there is phase info associated with the cloud; else use HG
  
  REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
  
! ---------------- outputs from the scattering tables -------------------
! --------------------- produced by Evans Mie code ----------------------
!     The scattering tables are read in with READ_SSCATTAB.  The scattering
!     table is 3D: wavenumber, particle size, and viewing angle.
!         Scattering table variables:
!       MUTAB is view angle values (cosine zenith),
!       DMETAB is particle size values (median mass diameter, micron),
!       WAVETAB is wavenumber values (cm^-1).
!       MUINC(2) are the mu values of the two incident angles
!       TABEXTINCT is extinction, TABSSALB is single scattering albedo,
!       TABASYM is the asymmetry parameter
!       TABPHI??? are phase function info for incident directions
  
!cc      INTEGER  MAXTAB, MAXGRID, MAXSCAT
!cc      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
  CHARACTER (LEN=80) :: SCATFILE(MAXSCAT)
  
  
  
  
  
  
  
  
  
  
  INTEGER :: NLEV, NABSNU
  
  REAL :: !ztop, zobs not needed
  INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
  
  
  
  
! ---------------------------- local variables ----------------------------
  INTEGER :: iaTable(kMaxClouds*kCloudLayers),iIn,iJ,iReadTable,I
  INTEGER :: iCloud,iStep
  REAL :: extinct
  INTEGER :: LL,II,N,M,iB_atm,iT_Atm,iLayers
  REAL :: dmetab_phase(kProfLayer)
  INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
  
  CHARACTER (LEN=80) :: caName
  CHARACTER (LEN=1) :: caScale(MAXSCAT)
  
!initialise all scattering info to null
  
!!!! need to reference the cloud tops and bottoms wrt TOP layer of
!!!! defined atm
!!!! eg if atm from 971 to 150 mb (plane at 150 mb) ==>
!!!!       this occupies kCARTA layers 19-78
!!!!   if cloud from 248 to 214 mb  ==>
!!!!       this occupies kCARTA layers 69 to 72
!!!! Thus the cloud occupies RTSPEC atmosphere "tau" from 6 to 9
  IF (iDownWard == 1) THEN
    iB_Atm=iaaRadLayer(iAtm,1)
    iT_Atm=iaaRadLayer(iAtm,iNumLayer)
  ELSE IF (iDownWard == -1) THEN
    iT_Atm=iaaRadLayer(iAtm,1)
    iB_Atm=iaaRadLayer(iAtm,iNumLayer)
  END IF
  
!!!!!!hwowever we also do fluxes, so even if the atm is defined so it
!!!!!!is for an uplook instrument, RTSPEC will be called in a downlook
!!!!!!fashion, and vice versa
  
  IF (iB_Atm > iT_Atm) THEN
    iCloudySky = iT_Atm
    iT_Atm = iB_Atm
    iB_atm = iCloudySky
  END IF
  
  iCloudySky = -1  !!!!!!! assume NO clouds associated with this atm
  
  iCldTopKcarta = -1
  iCldBotKcarta = kProfLayer+1
  NSCATTAB=-1000   !!!!!!! total of how many scattering tables (files)
  DO iIn=1,kMaxClouds*kCloudLayers
    iaTable(iIn) = -1
  END DO
  DO iIn=1,MAXSCAT
    ScatFile(iIn) = ' '
    iaScatTable_With_Atm(iIn) = -1
  END DO
  DO iIn=1,kMaxClouds
    iaCloudWithThisAtm(iIn) = -1
    iaCldTop(iIn) = -1
    iaCldBot(iIn) = -1
  END DO
  
!!!go thru info and check against whether should be used with this atm
  DO iIn=1,iNclouds
    
!!associate scattering tables with the clouds
!!initialise info for this iIn th cloud
    DO iJ=1,iaCloudNumLayers(iIn)
      iI=iaaScatTable(iIn,iJ)
      IF (iI > MAXSCAT) THEN
        WRITE(kStdErr,*)'in interface_disort, you have set it up so'
        WRITE(kStdErr,*)'MAXSCAT < kMaxClouds*kCloudLayers'
        WRITE(kStdErr,*)'please reset and retry'
        CALL DoSTOP
      END IF
      caName=caaaScatTable(iIn,iJ)
      IF (iaTable(iI) < 0) THEN  !nothing associated with this yet
        IF (iI > NSCATTAB) THEN
          NSCATTAB=iI
        END IF
        iaTable(iI) = 1
        ScatFile(iI) = caName
      END IF
    END DO
    
!!check to see if this cloud is to be used with this atm
    DO iJ=1,iaCloudNumAtm(iIn)
      IF (iaaCloudWhichAtm(iIn,iJ)  == iAtm) THEN
        iCloudySky = iIn         !!!! set this up
        iaCloudWithThisAtm(iIn) = 1
        
!!!!!these are the kCARTA layers 1 to 100 = GND to TOA
        IACLDTOP(iIn) = iaaCloudWhichLayers(iIn,1)+1
        IACLDBOT(iIn) = iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn))
        
        IF (iCldTopkCarta < iaCldTop(iIn)-1) THEN
          iCldTopkCarta = iaCldTop(iIn)-1
        END IF
        IF (iCldBotkCarta > iaCldBot(iIn)) THEN
          iCldBotkCarta = iaCldBot(iIn)
        END IF
        
        WRITE(kStdWarn,*)'cloud # ',iIn,' associated with atm # ',iAtm
        WRITE(kStdWarn,*)'cloud is in KCARTA layers ',  &
            iaCldTop(iIn)-1,' to ',iaCldBot(iIn)
        
!!!!!these are the DISORT layers 100 to 1 = GND to TOA
        iaCldbot(iIn) = iT_Atm - iaCldbot(iIn) + 1
        iaCldtop(iIn) = iT_Atm - iaCldtop(iIn) + 1
        
!            iaCldBot(iIn) = iaCldBot(iIn) + 1
!            iaCldTop(iIn) = iaCldTop(iIn) + 1
        
        WRITE(kStdWarn,*)'cloud is in DISORT layers ',  &
            iaCldTop(iIn)+1,' to ',iaCldBot(iIn)
        
      END IF
    END DO
    
!!check to see which scattering tables to be used with this atm
    DO iJ=1,iaCloudNumLayers(iIn)
      iI=iaaScatTable(iIn,iJ)
      IF (iaCloudWithThisAtm(iIn)  == 1) THEN
        iaScatTable_With_Atm(iI) = 1
        WRITE(kStdWarn,*)'scatter table ',iI,' used with atm # ',iAtm
      END IF
    END DO
    
  END DO      !!!!!!!!main       DO iIn=1,iNclouds
  
!     Only read in scattering tables that are needed for this atm
  iReadTable = 1
  IF (iReadTable > 0) THEN
    IF (iBinaryFile == 1) THEN
      DO I = 1, NSCATTAB
        IF (iaScatTable_With_Atm(I) > 0) THEN
          WRITE(kStdWarn,*) 'Reading binary scatter data for table #',I
          CALL READ_SSCATTAB_BINARY(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,  &
              caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I),  &
              NWAVETAB(I), WAVETAB(1,I),  &
              MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),  &
              TABPHI1UP(1,I), TABPHI1DN(1,I), TABPHI2UP(1,I), TABPHI2DN(1,I))
          IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR.  &
                (ABS(MUINC(2)-0.7887) > 0.001)) THEN
            WRITE(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
            CALL DoStop
          END IF
        END IF
        ENDDO
        ELSE IF (iBinaryFile == -1) THEN
          DO I = 1, NSCATTAB
            IF (iaScatTable_With_Atm(I) > 0) THEN
              WRITE(kStdWarn,*) 'Reading ascii scatter data for table #',I
              CALL READ_SSCATTAB(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,  &
                  caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I),  &
                  NWAVETAB(I), WAVETAB(1,I),  &
                  MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),  &
                  TABPHI1UP(1,I), TABPHI1DN(1,I), TABPHI2UP(1,I), TABPHI2DN(1,I))
              IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR.  &
                    (ABS(MUINC(2)-0.7887) > 0.001)) THEN
                WRITE(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                CALL DoStop
              END IF
            END IF
            ENDDO
            ELSE IF (iBinaryFile == 0) THEN
              DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                  WRITE(kStdWarn,*) 'Reading ascii scatter data for table #',I
                  CALL READ_SSCATTAB_SPECIAL(SCATFILE(I),  &
                      caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I),  &
                      NWAVETAB(I), WAVETAB(1,I),  &
                      MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),  &
                      TABPHI1UP(1,I), TABPHI1DN(1,I),  &
                      TABPHI2UP(1,I), TABPHI2DN(1,I))
                  
                END IF
                ENDDO
                END IF    !iBinaryFile .GT. 0
              END IF      !iReadTable  .GT. 0
              
!check to see that all MIE tables read in had the same nmuobs used
!assuming that this atm does have a cloud associated with it
              IF (iCloudySky > 0) THEN
                iLayers = 0
                LL = 0
                DO i=1,nscattab
                  IF (iaScatTable_With_Atm(I) > 0) THEN
                    iLayers = iLayers + nmuobs(i)
                    LL = LL + 1  !!keep track of how many scattering tables used
                    II = I       !!keep track of which scattering table used with atm
                  END IF
                END DO
                IF (INT(iLayers*1.0/LL) /= nmuobs(II)) THEN
                  WRITE (kStdErr,*) iLayers,LL,INT(iLayers*1.0/LL),nmuobs(II)
                  WRITE (kStdErr,*) 'Some of the Mie Scattering tables had different'
                  WRITE (kStdErr,*) 'number of angles used in computing Mie coeffs'
                  WRITE (kStdErr,*) 'Please recheck  sscatmie.x and rerun'
                  CALL DoStop
                END IF
              ELSE
                WRITE (kStdWarn,*) 'no clouds with atmosphere number ',iAtm,' !!!'
              END IF
              
! Frank Evans code scales the Mie scattering parameters, so we have to
! unscale them!!!!!!!!
!          DO I = 1, NSCATTAB
!            IF (iaScatTable_With_Atm(I).GT. 0) THEN
!              CALL UnScaleMie(
!     $          caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
!     $          ndme(i)*nwavetab(i))
!            END IF
!          END DO
              
! code from n_rad_jac_scat.f
              iCloud = -1
              IF (iCloudySky < 0) THEN
                WRITE(kStdWarn,*)'Could not find a cloud for atmosphere #',iAtm
                WRITE(kStdWarn,*)'setting IWP = -100.0'
                iCloud=1    !say cloud number one associated with this atmosphere
                ncldlay=1   !say fictitious cloud occupies one layer
                IWP(1)      = -100.0   !but ensure cloud has NO particles in it!
                DME(1)      = -10.0    !but ensure cloud has NO particles in it!
                ISCATTAB(1) = -1
                
              ELSE
!!!!!find total number of clouds, and hence total number of layers
                NCLDLAY=0
                iLayers = 0
                DO i=1,kMaxClouds
                  IF (iaCloudWithThisAtm(i) == 1) THEN
                    ncldlay = ncldlay + iaCloudNumLayers(i)
                    WRITE(kStdWarn,*) 'Cloud #, num layers = ',i,iaCloudNumLayers(i)
                    
                    WRITE(kStdWarn,*) 'L   iwp  dme  iscattab   kCARTA Layer  : '
                    
                    DO iStep = 1, iaCloudNumLayers(i)
                      iLayers = iLayers + 1
                      IWP(iLayers) = raaaCloudParams(i,iStep,1)
                      DME(iLayers) = raaaCloudParams(i,iStep,2)
                      ISCATTAB(iLayers) = iaaScatTable(i,iStep)
                      WRITE(kStdWarn,*) iLayers,iwp(iLayers),dme(iLayers),  &
                          iscattab(iLayers),iaaCloudWhichLayers(i,iStep)
                    END DO
                  END IF
                END DO
              END IF
              
!     Find the levels for top of cloud and observation level
!     remember that these numbers are with respect to the KLAYERS pressure
!     levels and layers
!     these will be reset when passed in and out of GetAbsProfileDISORT
!     NOTE : here we are still in kCARTA frame ie
              
!   TOA    --------------
!          layer iNumlayer
!          --------------
!              .....
!          --------------             KCARTA
!             layer 2
!          --------------
!             layer 1
!   GND --------------------------
              
! when we call GetAbsProfile, variables icldtop,icldbot,iobs will be reset
! to reflect the rtspec layering
!   TOA    --------------
!             layer 1
!          --------------
!              .....                 DISORT
!          --------------
!        layer iNumLayer-1
!          --------------
!         layer iNumLayer
!   GND --------------------------
              
              IF (IWP(1) <= 0.0) THEN  !we have no cloud; set up fictitious clouds
                IF (iDownWard > 0) THEN
!down look instr : set cloud BELOW observer, in kCARTA layer #1
                  ICLDTOP = 2
                  ICLDBOT = 1
!down look instr : set cloud BELOW observer, in DISORT layer #2
                  ICLDTOP = 2
                  ICLDBOT = 3
                  IOBS    = iNumLayer
                ELSE IF (iDownWard < 0) THEN    !up look instr
!up look instr : set cloud ABOVE observer, in kCARTA layer #iNumLayer
                  ICLDTOP = iNumLayer+1
                  ICLDBOT = iNumLayer
!up look instr : set cloud ABOVE observer, in DISORT layer #1
                  ICLDTOP = 2
                  ICLDBOT = 1
                  IOBS    = 1
                END IF
              END IF
              
              IF (IWP(1) > 0.0) THEN
!do not really need icldtop/bot, but just set it up
                ICLDTOP=iaaCloudWhichLayers(iCloudySky,1)+1
                ICLDBOT=iaaCloudWhichLayers(iCloudySky,iaCloudNumLayers(iCloudySky))
                
                icldbot = iT_Atm - icldbot + 1
                icldtop = iT_Atm - icldtop + 1
                
                icldbot = icldbot + 1
                icldtop = icldtop + 1
                
                IF (iDownWard > 0) THEN
                  IOBS   = iNumLayer
                ELSE IF (iDownWard < 0) THEN
                  IOBS   = 1
                END IF
              END IF
              
              iobs=(iNumLayer+1)-iobs+1
              
              30   FORMAT(I3,' ',A80)
              
              RETURN
            END SUBROUTINE SetMieTables_DISORT
            
!************************************************************************
! this subtroutine does some initializations for a uplook instr
            
            SUBROUTINE Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,  &
                rFracTop,raFreq,raaAbs,rSatAngle,iTag,  &
                raTopIntensity,raSolarBeam,TOA_to_instr)
            
            
            INTEGER, INTENT(IN OUT)                  :: iAtm
            NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
            INTEGER, INTENT(IN OUT)                  :: iNumLayer
            REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
            REAL, INTENT(IN OUT)                     :: rFracTop
            REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
            REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
            REAL, INTENT(IN OUT)                     :: rSatAngle
            INTEGER, INTENT(IN OUT)                  :: iTag
            NO TYPE, INTENT(IN OUT)                  :: raTopInten
            NO TYPE, INTENT(IN OUT)                  :: raSolarBea
            NO TYPE, INTENT(IN OUT)                  :: TOA_to_ins
            IMPLICIT NONE
            
            INCLUDE '../INCLUDE/kcartaparam.f90'
            
! input variables
            INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
            REAL :: !temperature profile
            
            
            REAL :: !matrix of abs coeffs
            
! output variables
            REAL :: raTopIntensity(kMaxPts)     !intensity at TOA
            REAL :: raSolarBeam(kMaxPts)        !intensity diue to solar rad
            REAL :: TOA_to_instr(kMaxPts)       !abs coeffs
            
! local variables
            INTEGER :: iaRadLayer(kProfLayer),iL,IF
            REAL :: ttorad,rDummy
            
            DO iL=1,kMaxPts
              TOA_to_instr(iL) = 0.0
            END DO
            
! bring incident space radiation down from TOA to instrument
            DO IF=1,kMaxPts
! compute the Plank radiation from space
              raTopIntensity(IF) = ttorad(raFreq(IF),SNGL(kTSpace))
            END DO
            
! set the solar beam intensity ...
! if the sun is ON and is in the FOV of the instrument, then the BC of 5700K
!   to the TOA is set, else the BC of 2.6K to TOA is set; also fbeam set to 0
! if the sun is ON and is NOT in the FOV of the instrument, then the fbeam is
!   set to solar beam, while BC of TOA is 2.6K is set
            
            DO IF=1,kMaxPts     !!!!assume sun is NOT ON
              raSolarBeam(IF) = 0.0
            END DO
            
            IF (kSolar >= 0) THEN
              IF (ABS(ABS(rSatAngle)-ABS(kSolarAngle)) >= 1.0E-3) THEN
!sun is on, but not in FOV of instr
!so set fbeam correctly to that of sun
                CALL SolarBeamDisort(kSolar,raSolarBeam,raFreq,iTag)
                rDummy = ABS(COS(kSolarAngle*kPi/180))
! bring solar radiation down from TOA to instrument
                CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,  &
                    rFracTop,raFreq,raaAbs,TOA_to_instr)
                DO IF=1,kMaxPts
                  raSolarBeam(IF) = raSolarBeam(IF)*EXP(-TOA_to_instr(IF)/rDummy)
                END DO
              END IF
              
              IF (ABS(ABS(rSatAngle)-ABS(kSolarAngle)) <= 1.0E-3) THEN
!sun is on, and in FOV of instr
!so set fbeam correctly to 0.0, and the BC of 5700K
              END IF
              
            END IF     !!IF (kSolar .GE. 0)
            
            RETURN
          END SUBROUTINE Init_UpLook
!************************************************************************
! this subtroutine does some initializations for a downlook instr
          
          SUBROUTINE Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,  &
              rFracTop,raFreq,raaAbs,rSatAngle,iTag,  &
              raTopIntensity,raSolarBeam,TOA_to_instr)
          
          
          INTEGER, INTENT(IN OUT)                  :: iAtm
          NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
          INTEGER, INTENT(IN OUT)                  :: iNumLayer
          REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
          REAL, INTENT(IN OUT)                     :: rFracTop
          REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
          REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
          REAL, INTENT(IN OUT)                     :: rSatAngle
          INTEGER, INTENT(IN OUT)                  :: iTag
          NO TYPE, INTENT(IN OUT)                  :: raTopInten
          NO TYPE, INTENT(IN OUT)                  :: raSolarBea
          NO TYPE, INTENT(IN OUT)                  :: TOA_to_ins
          IMPLICIT NONE
          
          INCLUDE '../INCLUDE/kcartaparam.f90'
          
! input variables
          INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
          REAL :: !temperature profile
          
          
          REAL :: !matrix of abs coeffs
          
! output variables
          REAL :: raTopIntensity(kMaxPts)     !intensity at TOA
          REAL :: raSolarBeam(kMaxPts)        !intensity diue to solar rad
          REAL :: TOA_to_instr(kMaxPts)       !abs coeffs
          
! local variables
          INTEGER :: iaRadLayer(kProfLayer),iL,IF
          REAL :: ttorad,rDummy
          
! see if there are any layers between TOA and instr; if there are, set
! TOA_to_instr to the cumulative k(TOA to aircraft), else set it to 0
          DO iL=1,kProfLayer
            iaRadLayer(iL) = iaaRadLayer(iAtm,iL)
          END DO
          CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,  &
              rFracTop,raFreq,raaAbs,TOA_to_instr)
          
! bring incident space radiation down from TOA to instrument
          DO IF=1,kMaxPts
            raTopIntensity(IF) = ttorad(raFreq(IF),SNGL(kTSpace))
          END DO
! this is technically incorrect ... we really should do rad transfer here
          DO IF=1,kMaxpts
            raTopIntensity(IF) = raTopIntensity(IF)*EXP(-TOA_to_instr(IF))
          END DO
          
          IF (kSolar >= 0) THEN
            CALL SolarBeamDisort(kSolar,raSolarBeam,raFreq,iTag)
            rDummy = ABS(COS(kSolarAngle*kPi/180))
! bring solar radiation down from TOA to instrument
            DO IF=1,kMaxPts
              raSolarBeam(IF) = raSolarBeam(IF)*EXP(-TOA_to_instr(IF)/rDummy)
            END DO
          ELSE
            DO IF=1,kMaxPts
              raSolarBeam(IF) = 0.0
            END DO
          END IF
          
          RETURN
        END SUBROUTINE Init_DownLook
!************************************************************************
! set raCorrelatedK to abs coeffs of layer nearest the ground
        
        SUBROUTINE SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,  &
            iaaRadLayer,iAtm,iNumLayer, ABSPROF)
        
        
        INTEGER, INTENT(IN OUT)                  :: iDownWard
        NO TYPE, INTENT(IN OUT)                  :: raCorrelat
        REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
        NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
        INTEGER, INTENT(IN OUT)                  :: iAtm
        INTEGER, INTENT(IN)                      :: iNumLayer
        REAL, INTENT(IN OUT)                     :: ABSPROF(MAXNZ,MAXABSNU)
        IMPLICIT NONE
        
        INCLUDE '../INCLUDE/scatterparam.f90'
        
        
        REAL :: raCorrelatedK(kMaxPts)  !to see how the "k distributions" are
        REAL :: !matrix of abs coeffs
        INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
        
        
        INTEGER :: iI,iJ,iaIndx(kMaxPts),iMethod
        REAL :: raL2S(kMaxPts),raaWgt(kMaxPts,kProfLayer),raPeak(kMaxPts)
        
        iMethod = -1
        
        IF (iMethod == -1) THEN
!simply store optical depth of layer closest to gnd
          IF (iDownWard == 1) THEN
!! as k(gnd) increases, down look instr only penetrates less into
!! atmosphere, so radiance decreases
            DO iI=1,kMaxPts
              raCorrelatedK(iI) = raaAbs(iI,iaaRadLayer(iAtm,1))
            END DO
          ELSE
!! as k(gnd) increases, uplook instr only penetrates less into
!! atmosphere, so radiance increases as you see hot stuff in
!!your face
            DO iI=1,kMaxPts
              raCorrelatedK(iI) = raaAbs(iI,iaaRadLayer(iAtm,iNumLayer))
            END DO
          END IF
        END IF
        
        IF (iMethod == +1) THEN
!store optical depth of layer that peaks the weight fcn
!!!!!!remember ABSPROF(1,:) = TOA, ABSPROF(NLEV-1,:) = GND
          IF (iDownward == 1) THEN
!do down look weight fcn
            DO iI=1,kMaxPts
              raL2S(iI)  = 0.0        !!!optical depth to space = 0.0
              iaIndx(iI) = 1
              raPeak(iI) = -1.0E10
            END DO
            DO iI = 1,kMaxPts
              DO iJ = iNumLayer,1,-1
                raaWgt(iI,iJ) = absprof(iNumLayer-iJ+1,iI)
                raaWgt(iI,iJ) = (1-EXP(-raaWgt(iI,iJ)))*EXP(-raL2S(iI))
                raL2S(iI)     = raL2S(iI) + absprof(iNumLayer-iJ+1,iI)
              END DO
            END DO
          ELSE IF (iDownward == -1) THEN
!do up look weight fcn
            DO iI=1,kMaxPts
              raL2S(iI) = 0.0        !!!optical depth to gnd = 0.0
              iaIndx(iI) = 1
              raPeak(iI) = -1.0E10
            END DO
            DO iI = 1,kMaxPts
              DO iJ = 1,iNumLayer
                raaWgt(iI,iJ) = absprof(iNumLayer-iJ+1,iI)
                raaWgt(iI,iJ) = (1-EXP(-raaWgt(iI,iJ)))*EXP(-raL2S(iI))
                raL2S(iI)     = raL2S(iI) + absprof(iNumLayer-iJ+1,iI)
              END DO
            END DO
          END IF
!!!having stored the wgt fcns, find where it peaks
          DO iI = 1,kMaxPts
            DO iJ = 1,iNumLayer
              IF (raaWgt(iI,iJ) > raPeak(iI)) THEN
                raPeak(iI) = raaWgt(iI,iJ)
                iaIndx(iI) = iJ
              END IF
            END DO
          END DO
!!!now set this info to raCorrelatedK
          DO iI = 1,kMaxPts
            raCorrelatedK(iI) = absprof(iNumLayer-iaIndx(iI) + 1,iI)
          END DO
        END IF
        
        RETURN
      END SUBROUTINE SetCorrelatedK
      
!************************************************************************
! set up the single scatter albedos etc for the clouds
      
      SUBROUTINE SetUpClouds(nstr, nmuobs, iaCloudWithThisAtm,  &
          iaCldTop,iaCldBot,iaCloudNumLayers,rF, iAtm,iaaRadLayer,iNumLayer,  &
          IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB,  &
          TABEXTINCT, TABSSALB, TABASYM, ISCATTAB, extinct,dtauc,ssalb,asym,pmom)
      
      
      INTEGER, INTENT(IN OUT)                  :: nstr
      INTEGER, INTENT(IN)                      :: nmuobs
      NO TYPE, INTENT(IN OUT)                  :: iaCloudWit
      NO TYPE, INTENT(IN OUT)                  :: iaCldTop
      NO TYPE, INTENT(IN OUT)                  :: iaCldBot
      NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
      REAL, INTENT(IN OUT)                     :: rF
      INTEGER, INTENT(IN OUT)                  :: iAtm
      NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
      INTEGER, INTENT(IN OUT)                  :: iNumLayer
      REAL, INTENT(IN OUT)                     :: IWP(MAXNZ)
      REAL, INTENT(IN OUT)                     :: DME(MAXNZ)
      INTEGER, INTENT(IN OUT)                  :: NDME(MAXSCAT)
      REAL, INTENT(IN OUT)                     :: DMETAB(MAXGRID,MAXSCAT)
      INTEGER, INTENT(IN OUT)                  :: NWAVETAB(MAXSCAT)
      REAL, INTENT(IN OUT)                     :: WAVETAB(MAXGRID,MAXSCAT)
      REAL, INTENT(IN OUT)                     :: TABEXTINCT(MAXTAB,MAXSCAT)
      REAL, INTENT(IN OUT)                     :: TABSSALB(MAXTAB,MAXSCAT)
      REAL, INTENT(IN OUT)                     :: TABASYM(MAXTAB,MAXSCAT)
      INTEGER, INTENT(IN)                      :: ISCATTAB(MAXNZ)
      REAL, INTENT(IN OUT)                     :: extinct
      DOUBLE PRECISION, INTENT(IN OUT)         :: dtauc(maxcly)
      DOUBLE PRECISION, INTENT(IN OUT)         :: ssalb(maxcly)
      NO TYPE, INTENT(OUT)                     :: asym
      DOUBLE PRECISION, INTENT(IN OUT)         :: pmom(0:maxmom,maxcly)
      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/scatterparam.f90'
      
! inputs
      INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
      
      INTEGER :: iaCloudWithThisAtm(kMaxClouds)  !is this cloud in this atm?
      INTEGER :: IACLDTOP(kMaxClouds)            !cloud top layer
      INTEGER :: IACLDBOT(kMaxClouds)            !cloud bot layer
      INTEGER :: iaCloudNumLayers(kMaxClouds)    !number of layers cloud occupies
      REAL :: !iwp, particle size
      
      
      
      
      
      
! outputs
      
      DOUBLE PRECISION :: !optical depths; also used as INPUT
      DOUBLE PRECISION :: ASYM(maxnz)
      DOUBLE PRECISION :: !scattering phase fcn
      DOUBLE PRECISION :: !single scatter albedo of layers
      
! local variables
      INTEGER :: LL,M,ICLDTOP,ICLDBOT,N,L,I,nmom, iiDiv,Nprime
      REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
      DOUBLE PRECISION :: tauC(kProfLayer),tauCG(kProfLayer)
      
      iiDiv = 0
      555  CONTINUE
      IF (iiDiv*kProfLayer < iaaRadLayer(iAtm,3)) THEN
        iiDiv = iiDiv + 1
        GO TO 555
      END IF
      iiDiv = iiDiv - 1
      
!!!!!!!!!!!!!! ********* CLOUD SCATTERING ************* !!!!!!!!!!!!
      LL=0
      DO M = 1,kMaxClouds
        IF (iaCloudWithThisAtm(M) > 0) THEN
          ICLDTOP = IACLDTOP(M)
          ICLDBOT = IACLDBOT(M)
          DO N = ICLDTOP, ICLDBOT - 1
            Nprime = N-iiDiv*kProfLayer
!L is the cloud layer = 1(top)..nlay(bot) in cloud
            L = LL + N-ICLDTOP+1
            
!            write (kstdwarn,*) 'm,tp,bt,l,iwp(l),dme(l) = ',
!     $ m,icldtop,icldbot,l,iwp(l),dme(l)
            
            I = ISCATTAB(L)     !!!!!I is the scattering table info number
            
!      Interpolate to get values of extinction, s.s. albedo, and
!      phi function values at given obs. mu, waveno, and particle size.
!      Note: we don't need phi functions for Eddington-only solution
!      This means that while the rtspec code had a choice of
!            CALL INTERP_SCAT_TABLE2 (rF, DME(L),    versus
!            CALL INTERP_SCAT_TABLE3 (rF, DME(L),
!      here we only need the simpler first choice as we are not messing
!      around with the phase functions
            CALL INTERP_SCAT_TABLE2 (rF, DME(L),  &
                EXTINCT, SSALB_RTSPEC(L), ASYM_RTSPEC(L),  &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I),  &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
            
! but the indices have to be modified so we put info into the correct layers
! also have to set asymmetry parameter
            TAUC(L) = DBLE(IWP(L)*EXTINCT/1000.0)
!!!!!!!!!!!! orig code TAUCG(L) = TAUGAS(N) + TAUC(L)
            TAUCG(L) = dtauc(Nprime) + TAUC(L)
            IF (TAUCG(L) > 1.0D-5) THEN
              SSALB_RTSPEC(L) = SSALB_RTSPEC(L)*TAUC(L)/TAUCG(L)
            ELSE
              SSALB_RTSPEC(L) = 0.0
            END IF
! so now set the indices correctly and do relevant changes of real --> double
!!!orig code dtauc(N)      = blah, ssalb(N)      = blah
!!!onew code dtauc(Nprime) = blah, ssalb(Nprime) = blah
            dtauc(Nprime) = DBLE(TAUCG(L))
            SSALB(Nprime) = DBLE(SSALB_RTSPEC(L))
            IF (ssalb_rtspec(L) > 1.0E-6) THEN
              asym(Nprime)  = DBLE(asym_rtspec(L))
            ELSE
              asym(Nprime)  = DBLE(0.0)
            END IF
            
!!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
            nmom = MAX(2*nmuobs + 1,nstr)
            nmom = MIN(maxmom,nmom)
!! get the phase moments for the HG phase function
            CALL getmom(3,asym(Nprime),nmom,pmom(0,Nprime))
            ENDDO                 !DO N=ICLDTOP,ICLDBOT-1
              LL = LL + iaCloudNumLayers(M)
            END IF                  !IF (iaCloudWithThisAtm(M) .GT. 0) THEN
          END DO                    !DO M = 1,kMaxClouds
          
          RETURN
        END SUBROUTINE SetUpClouds
        
!************************************************************************
! this subroutine sets up the optical depths, single scatter albedo etc for
! a atmosphere where there is only gas + Rayleigh scattering
        
        SUBROUTINE SetUpRayleigh(nlev,nstr,nmuobs, rF,raDensity,raThickness,  &
            dtauc,ssalb,asym,pmom)
        
        
        INTEGER, INTENT(IN OUT)                  :: nlev
        NO TYPE, INTENT(IN OUT)                  :: nstr
        INTEGER, INTENT(IN)                      :: nmuobs
        REAL, INTENT(IN OUT)                     :: rF
        REAL, INTENT(IN)                         :: raDensity(kProfLayer)
        NO TYPE, INTENT(IN OUT)                  :: raThicknes
        DOUBLE PRECISION, INTENT(IN OUT)         :: dtauc(maxcly)
        DOUBLE PRECISION, INTENT(IN OUT)         :: ssalb(maxcly)
        NO TYPE, INTENT(OUT)                     :: asym
        DOUBLE PRECISION, INTENT(IN OUT)         :: pmom(0:maxmom,maxcly)
        IMPLICIT NONE
        
        INCLUDE '../INCLUDE/scatterparam.f90'
! inputs
        INTEGER :: nstr          !number of levels, number of streams
        
        
        REAL :: raThickness(kProfLayer)
! outputs
        DOUBLE PRECISION :: !optical depths; also used as INPUT
        DOUBLE PRECISION :: ASYM(maxnz)
        DOUBLE PRECISION :: !scattering phase fcn
        DOUBLE PRECISION :: !single scatter albedo of layers
        
! local variables
        REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
        DOUBLE PRECISION :: tauC(kProfLayer),tauCG(kProfLayer)
        REAL :: Rayleigh
        INTEGER :: N,nmom
        
        DO N = 1,NLEV-1   !!!!!!!to include scattering
!indices have to be modified so we put info into the correct layers
          TAUC(N) = DBLE(rayleigh(rF,raDensity(N),raThickness(N)))
          TAUCG(N) = dtauc(N) + TAUC(N)
          SSALB_RTSPEC(N) = SNGL(TAUC(N)/TAUCG(N))
! so now set the indices correctly and do relevant changes of real --> double
          dtauc(N) = TAUCG(N)
          SSALB(N) = DBLE(SSALB_RTSPEC(N))
          asym(N)  = DBLE(0.0)
          
!!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
          nmom = MAX(2*nmuobs + 1,nstr)
          nmom = MIN(maxmom,nmom)
!! get the phase moments for Rayleigh scattering
          CALL getmom(2,asym(N),nmom,pmom(0,N))
          ENDDO
            
            RETURN
          END SUBROUTINE SetUpRayleigh
          
!************************************************************************
! this does the final initializations before calling DISORT
          
          SUBROUTINE FinalInitialization( !!!!inputs  &
              iDownWard,rSatAngle,rTopIntensity,rSolarBeam,emiss,  &
              rSurfaceTemp,dtauc,dTotalOpticalDepth,iDoFlux,nlev, iNp,iaOp,  &
!!!!outputs  &
          usrtau,ntau,utau,usrang,numu,umu, nphi,phi,fisot,fbeam,umu0,phi0,  &
              ibcnd,lamber,albedo,  &
              btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)
          
          
          INTEGER, INTENT(IN OUT)                  :: iDownWard
          REAL, INTENT(IN OUT)                     :: rSatAngle
          NO TYPE, INTENT(IN OUT)                  :: rTopIntens
          REAL, INTENT(IN OUT)                     :: rSolarBeam
          REAL, INTENT(IN OUT)                     :: emiss
          NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
          DOUBLE PRECISION, INTENT(IN)             :: dtauc(maxcly)
          NO TYPE, INTENT(IN OUT)                  :: dTotalOpti
          INTEGER, INTENT(IN OUT)                  :: iDoFlux
          INTEGER, INTENT(IN)                      :: nlev
          INTEGER, INTENT(IN)                      :: iNp
          INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
          LOGICAL, INTENT(OUT)                     :: usrtau
          INTEGER, INTENT(OUT)                     :: ntau
          DOUBLE PRECISION, INTENT(OUT)            :: utau(maxulv)
          LOGICAL, INTENT(OUT)                     :: usrang
          INTEGER, INTENT(OUT)                     :: numu
          DOUBLE PRECISION, INTENT(OUT)            :: umu(maxumu)
          INTEGER, INTENT(OUT)                     :: nphi
          DOUBLE PRECISION, INTENT(OUT)            :: phi(maxphi)
          DOUBLE PRECISION, INTENT(OUT)            :: fisot
          DOUBLE PRECISION, INTENT(OUT)            :: fbeam
          DOUBLE PRECISION, INTENT(OUT)            :: umu0
          DOUBLE PRECISION, INTENT(OUT)            :: phi0
          INTEGER, INTENT(OUT)                     :: ibcnd
          LOGICAL, INTENT(OUT)                     :: lamber
          DOUBLE PRECISION, INTENT(OUT)            :: albedo
          DOUBLE PRECISION, INTENT(OUT)            :: btemp
          DOUBLE PRECISION, INTENT(OUT)            :: ttemp
          DOUBLE PRECISION, INTENT(OUT)            :: temis
          LOGICAL, INTENT(OUT)                     :: plank
          LOGICAL, INTENT(OUT)                     :: onlyfl
          DOUBLE PRECISION, INTENT(OUT)            :: accur
          NO TYPE, INTENT(OUT)                     :: prnt
          NO TYPE, INTENT(IN OUT)                  :: header
          IMPLICIT NONE
          
          INCLUDE '../INCLUDE/scatterparam.f90'
          
! input variables
          
          DOUBLE PRECISION :: !optical depths; also used as INPUT
          DOUBLE PRECISION :: dTotalOpticalDepth
          REAL :: rTopIntensity, rSurfaceTemp
          
! output variables
          
          
          DOUBLE PRECISION :: !tau's at which to output results
          DOUBLE PRECISION :: !ang's at which to output results
          DOUBLE PRECISION :: !azimuthal phi's to output radiance
          
          
          
          
          
          
          CHARACTER (LEN=127             !dumb comment) :: header
          
          LOGICAL :: prnt(5)                   !prnt(1) = true, print input variables
          
! local variables
          INTEGER :: iI,iJ
          DOUBLE PRECISION :: d1,d2,dCumulative(maxcly)
          
          d1 = dtauc(1)
          d2 = dTotalOpticalDepth
          
          IF (iDownward == 1) THEN
!down look instrument
            usrtau = .true.   !return intensity at ONE optical depth
!instead of at all levels
            
            ntau    = 1       !return intensity at one user level
            utau(1) = DBLE(0.000000000)                  !!! for TOA
            
!!!!allow more than one level
            ntau = iNp
            dCumulative(1) = dtauc(1)
            DO iI = 2,nlev-1
              dCumulative(iI) = 0.0
              dCumulative(iI) = dtauc(iI) + dCumulative(iI-1)
            END DO
            
            DO iI = 1,ntau
!!we will get out radiance at TOP of layer
              utau(iI) = dCumulative((nlev-1)-iaOp(iI)+1) - d1
            END DO
            
            usrang = .true.   !return intensity at one user polar angle
            numu = 1
            umu(1) = DBLE(ABS(COS(rSatAngle*kPi/180)))
            nphi = 1          !return intensity at one user zenith angle
            phi(1) = DBLE(0.0)
            
            
          ELSE
!up look instrument
            usrtau = .true.   !return intensity at ONE optical depth
!instead of at all levels
            
            ntau    = 1       !return intensity at one user level
            utau(1) = dTotalOpticalDepth                 !!! for bottom
            
!!!!allow more than one level
            ntau = iNp
            dCumulative(1) = dtauc(1)
            DO iI = 2,nlev-1
              dCumulative(iI) = 0.0
              dCumulative(iI) = dtauc(iI) + dCumulative(iI-1)
            END DO
            
            DO iI = 1,ntau
!!we will get out radiance at BOTTOM of layer; flip iaOp stuff
              utau(iI) = dCumulative(iaOp(iI)) - d1
            END DO
            
            usrang = .true.   !return intensity at one user polar angle
            numu = 1
            umu(1) = DBLE(-ABS(COS(rSatAngle*kPi/180)))
            nphi = 1          !return intensity at one user zenith angle
            phi(1) = DBLE(0.0)
          END IF
          
          fisot  = DBLE(rTopIntensity)
          
          IF (ABS(ABS(rSatAngle) - ABS(kSolarAngle)) <= 1E-7) THEN
            kSolarAngle = kSolarAngle+1E-4
          END IF
          fbeam  = DBLE(rSolarBeam)
          umu0   = DBLE(COS(kSolarAngle*kPi/180))
          phi0   = DBLE(0.0)
          
          ibcnd = 0           !general case
          lamber = .false.    !specify bidir reflectance
          lamber = .true.     !isotropic reflect lower bdry ==> specify albedo
          albedo  = DBLE(1.0-emiss)   !from defns in books,papers
          
          btemp  = DBLE(rSurfaceTemp)      !!ground temperature
          IF (iDownward == 1) THEN
!down look instrument
            ttemp  = DBLE(kTSpace)         !!2.7 Kelvin
          ELSE
!up look instrument
            IF (kSolar < 0) THEN
              ttemp  = DBLE(kTSpace)       !!2.7 Kelvin or 5700K!!!!!!
            ELSE IF (kSolar >= 0) THEN
              IF (ABS(ABS(kSolarAngle)-ABS(rSatAngle)) <= 1.0E-3) THEN
                ttemp  = DBLE(kSunTemp)    !!sun in FOV .. 5700K!!!!!!
              ELSE
                ttemp  = DBLE(kTSpace)     !!sun not in FOV .. 2.7!!!!!!
              END IF
            END IF
          END IF
          
          temis  = DBLE(1.0)
          plank  = .true.      !emission from layers
          
          IF (iDoFlux == -1) THEN
            onlyfl = .false.     !only need intensity as well
            header  = 'scattering computations'
          ELSE IF (iDoFlux == +1) THEN
            onlyfl = .true.     !compute flux only
            usrtau = .false.    !compute at every boundary
            ntau = nlev
            numu = 0
            nphi = 0
            header  = 'flux computations'
          END IF
          
          accur   = DBLE(0.000)
          prnt(1) = .false.
          prnt(2) = .false.
          prnt(3) = .false.
          prnt(4) = .false.
          prnt(5) = .false.
          
          RETURN
        END SUBROUTINE FinalInitialization
        
!************************************************************************
! this is the main interface to DISORT
        
        SUBROUTINE interface_disort(
!first the usual kCARTA variables  &
        raFreq,raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,
!     $        iNp,iaOp,raaOp,iNpmix,iFileID,  &
        iNp,iaOp,iNpmix,iFileID, caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
!     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
!     $        raLayAngles,raSunAngles,  &
        raLayAngles,iDoFlux,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
            iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)
        
        IMPLICIT NONE
        
        INCLUDE '../INCLUDE/scatterparam.f90'
        
! iDoFlux     = do flux (+1) or radiance (-1)
! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
        REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
        REAL :: pProf(kProfLayer)
        INTEGER :: iProfileLayers
        REAL :: raNumberDensity(kProfLayer)
        REAL :: raLayAngles(kProfLayer)
        REAL :: raaAbs(kMaxPts,kMixFilRows)
        REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
        REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
        REAL :: rFracTop,rFracBot,rSurfPress
        REAL :: raInten(kMaxPts)
        INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag,iDoFlux
        INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
        INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
        CHARACTER (LEN=80) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
        INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
        INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
        INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
        INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
        CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
        REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
        INTEGER :: iaPhase(kMaxClouds)
        
! local variables
        CHARACTER (LEN=80) :: SCATFILE(MAXSCAT)
        
        INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
        REAL :: MUTAB(MAXGRID,MAXSCAT)
        REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
        REAL :: MUINC(2)
        REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
        REAL :: TABASYM(MAXTAB,MAXSCAT)
        REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
        REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
        
        INTEGER :: NSCATTAB, NCLDLAY, NLEV, NABSNU
        INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
        REAL :: IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
        INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
        INTEGER :: iCloudySky
        INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
        
! ---------------------------------------------------------------------
        
!   RTSPEC Radiative transfer variables:
        REAL :: TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
        
! ---------------------------------------------------------------------
        
! these are direct from DISOTEST.F
        
        CHARACTER (LEN=127          !dumb comment) :: header
        
        LOGICAL :: lamber  !true  ==> isotropic reflecting lower bdry
!          so need to specify albedo
!false ==> bidirectionally reflecting bottom bdry
!      ==> need to specify fcn BDREF()
        LOGICAL :: plank   !use the plank function for local emission
        LOGICAL :: onlyfl  !true ==> return flux, flux divergence, mean intensitys
!falsetrue ==> return flux, flux divergence, mean
!              intensities and intensities
        LOGICAL :: prnt(5) !prnt(1) = true, print input variables
!prnt(2) = true, print fluxes
!prnt(3) = true, intensities at user angles and levels
!prnt(4) = true, print planar transmitivity and albedo
!prnt(5) = true, print phase function moments PMOM
        LOGICAL :: usrtau  !false ==> rad quantities return at every bdry
!true  ==> rad quantities return at NTAU optical depths
!          which will be specified in utau(1:ntau)
        LOGICAL :: usrang  !false ==> rad quantities returned at stream angles
!true  ==> rad quantities returned at user angles
!          which will be specified in umu(1:numu)
        
        INTEGER :: ibcnd            !0 ==> general case beam (fbeam), isotropic
!      top illumination (fisot), thermal top
!      emission (temis,ttemp),internal thermal
!      sources (temper), reflection at bottom
!      (lamber, albedo, bdref), thermal
!      emission from bottom (btemp)
!1 ==> return only albedo/trans of entire
!      medium vs incident beam angle
        INTEGER :: nmom             !number of legendre phase polynoms ie phase fcn
        INTEGER :: nlyr             !number of computational layers in DISORT
!so nlvl = nlyr + 1
        INTEGER :: nstr             !number of radiation streams
        INTEGER :: ntau             !associated with LOGICAL usrtau, print results
!at this many optical depths
        INTEGER :: numu             !associated with LOGICAL usrang, specifies how
!many polar angles results to be printed (umu)
        INTEGER :: nphi             !specifies how many azimuth angles results to
!be printed (phi) .. can only be 0 if
!onlyfl = .true.
        
        DOUBLE PRECISION :: accur   !accuracy convergence criterion for azimuth
!(fourier cosine) series .. set between 0-0.01
        DOUBLE PRECISION :: albedo  !bottom bdry albedo
        DOUBLE PRECISION :: btemp   !bottom surface temp
        DOUBLE PRECISION :: dtauc(maxcly)
!optical depths of computational layers
        DOUBLE PRECISION :: fisot   !intensity of top bdry isotropic illumination
        DOUBLE PRECISION :: fbeam   !intensity of incident // beam at TOA
        DOUBLE PRECISION :: phi(maxphi)
!the azimuthal phi's to output radiance
        DOUBLE PRECISION :: pmom(0:maxmom,maxcly)
!scattering phase fcn in legendre polynoms
        DOUBLE PRECISION :: phi0    !solar azimuth
        DOUBLE PRECISION :: ssalb(maxcly)
!single scatter albedo of computational layers
        DOUBLE PRECISION :: temper(0:maxcly) !temperature of the levels (0=TOA)
        DOUBLE PRECISION :: temis            !emissivity of top bdry
        DOUBLE PRECISION :: ttemp            !temperature of top bdry
        DOUBLE PRECISION :: wvnmhi, wvnmlo   !bounds within which to do computations
        DOUBLE PRECISION :: umu(maxumu)      !ang's at which to output results
        DOUBLE PRECISION :: umu0         !polar angle cosine of incident solar beam
        DOUBLE PRECISION :: utau(maxulv)     !tau's at which to output results
        
!!!!!output variables
        DOUBLE PRECISION :: dfdt(maxulv)  !flux diverge d(net flux)/d(optical depth)
!where 'net flux' includes direct beam
        DOUBLE PRECISION :: flup(maxulv)  !diffuse up flux
        DOUBLE PRECISION :: rfldir(maxulv)!direct beam flux (without delta scaling)
        DOUBLE PRECISION :: rfldn(maxulv) !diffuse down flux = total-direct beam
        DOUBLE PRECISION :: uavg(maxulv)  !mean intensity (including direct beam)
        DOUBLE PRECISION :: UU( MAXUMU, MAXULV, MAXPHI )
!intensity if ONLYFL = .false., 0 o/w
        DOUBLE PRECISION :: albmed(maxumu)!albedo of medium as fcn of cos(umu(*))
!only set if ibcn == 1
        DOUBLE PRECISION :: trnmed(maxumu)!transmission as fcn of cos(umu(*))
!only set if ibcn == 1
! ---------------------------------------------------------------------
        DOUBLE PRECISION :: dTotalOpticalDepth
        DOUBLE PRECISION :: ASYM(maxnz)
! ---------------------------------------------------------------------
        
! these variables are to get the parameters from Frank Evans Mie Code
        REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
        
! this is to do "correlated k" computations
        REAL :: raCorrelatedK(kMaxPts)
        
! more local variables
        INTEGER :: iaStep(kMaxPts),iDiv,iScatter
        CHARACTER (LEN=80) :: caName
        INTEGER :: iIn,iJ,iI,iScat,iIOUN,IF,iFF,iFFMax,iL
        REAL :: TOA_to_instr(kMaxPts)
        INTEGER :: iaRadLayer(kProfLayer)
        REAL :: raSolarBeam(kMaxPts),raTopIntensity(kMaxPts)
        REAL :: ttorad,rDummy
        INTEGER :: iReadTable,iII
        INTEGER :: iDumper,iRayleigh
        REAL :: raDensity(kProfLayer)
        REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle
        
! these parameters are to step thru some of the 10000 pts
        REAL :: raaIntenSTEP(kProfLayer,kMaxPts),raIntenSTEP(kMaxPts)
        REAL :: raaNoscatterSTEP(kProfLayer,kMaxPts),raNoscatterSTEP(kMaxPts)
        REAL :: rakSTEP(kMaxPts),raFreqSTEP(kMaxPts),radtot
        INTEGER :: iSTEP,iStepPts
        
! these variables are to get the parameters from Frank Evans Mie Code
        REAL :: extinct
        INTEGER :: LL,L,I,N,M
        
! this is for RAYLEIGH
        REAL :: raThicknessRayleigh(kProfLayer)
        
! this is to convert from W to mW
        REAL :: raInten2(kMaxPts)
        
! this is if you want a funky dunky phase function
        REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
        
        rSolarAngle = kSolarAngle
        
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  kScatter = 1 works good for both up and down look
!  kScatter = 3 works good for up look
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! Rayleigh scatter has hardly any effect at infrared wavenumbers
        iRayleigh = -1       !!! -1 : do not want Rayleigh Scattering
!!! +1 : do Rayleigh Scattering
        
        iIOUN = kStdkCarta
        
        DO iIn=1,kMaxPts
          raIntenSTEP(iIn) = 0.0
          raFreqSTEP(iIn)  = 0.0
          rakSTEP(iIn)     = 0.0
        END DO
        
        WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',COS(rSatAngle*kPi/180.0),  &
            rSurfaceTemp
        
        CALL SetMieTables_DISORT(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
            raPhasePoints,raComputedPhase,  &
            iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
            TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
            TABPHI2UP, TABPHI2DN,  &
            NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
            IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
            iCloudySky, IACLDTOP, IACLDBOT)
        
        CALL GetAbsProfileDISORT(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
            iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
            NABSNU, NLEV, TEMP, ABSPROF,  &
            ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity,  &
            raDensity,raLayerTemp,  &
            iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)
        
!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
        IF (iCloudySky < 0) THEN
!!!!note that we do not care about background thermal accurately here
          WRITE (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
          DO iii = 1,iNp
            DO IF = 1,kMaxPts
              DO iL = 1,NLEV-1
                raTau(iL)  = absprof(iL,IF)
              END DO
              CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,  &
                  rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,  &
                  raUseEmissivity(IF),raFreq(IF),raInten(IF),iaOp(iii),  &
                  iProfileLayers,raPressLevels)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
          END DO
          GO TO 9876
        END IF
        
        WRITE(kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
        WRITE(kStdWarn,*) 'for the clear <-> cloudy comparisons'
!this fills up raaNoScatterStep
        DO iii = 1,iNumlayer
          DO IF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,IF)
            END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,  &
                rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,  &
                raUseEmissivity(IF),raFreq(IF),raaNoScatterStep(iii,IF),  &
                iaOp(iii),iProfileLayers,raPressLevels)
          END DO
        END DO
        
!!!!!!!!!! if iCloudSky .GT. 0 go thru and do the DISORT stuff !!!!!!!!
        CALL SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,  &
            iaaRadLayer,iAtm,iNumLayer,ABSPROF)
        
! set up some things for the instrument
        IF (iDownward == 1) THEN             !!down look instr
          CALL Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,  &
              rFracTop,raFreq,raaAbs,rSatAngle,iTag,  &
              raTopIntensity,raSolarBeam,TOA_to_instr)
        ELSE
          CALL Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,  &
              rFracTop,raFreq,raaAbs,rSatAngle,iTag,  &
              raTopIntensity,raSolarBeam,TOA_to_instr)
        END IF
        
! set the temperature levels, changing to double precision
        DO iL=1,NLEV
          temper(iL-1) = DBLE(temp(iL))
        END DO
        
! **************************** SPEED UP CODE ***************************
        nstr  = kDis_nstr   ! number of streams used by DISORT (2,4,8,16 ...)
        iStep = kDis_pts    ! number of wavenumber pts to use (1,2,...,10000)
! out of 10000
! **************************** SPEED UP CODE ***************************
        IF (iStep > kMaxPts) THEN
          WRITE(kStdWarn,*) 'Resetting kDis_Pts to kMaxPts'
          iStep = kMaxPts
        END IF
        
        IF (iStep < 20) THEN
          WRITE(kStdWarn,*) 'Resetting kDis_Pts to 20'
          iStep = 20
        END IF
        
!if you want to do 10 pts out of 10000 pts, then you have to do rad
!transfer on points 1,1001,2001,3001,...10000
!ie step over 10000/iStep points
        IF (kScatter /= 2) THEN
          iStep = iDiv(kMaxPts,iStep)
        END IF
        
! set up array of wavenumber indices that we step thru, to do the rad transfer
! (as DISORT is quite slow, we will not use each and every point)
        iScatter = kScatter
        CALL FindWavenumberPoints(iStep,iScatter,raCorrelatedK,iaStep,iFFMax)
        
        iStepPts = 0
        DO iFF = 1,iFFMax
          IF = iaStep(iFF)
          iStepPts = iStepPts + 1
          DO iL=1,NLEV-1
            dtauc(iL) = DBLE(absprof(iL,IF))
            raTau(iL) = absprof(iL,IF)
          END DO
          
          wvnmlo = DBLE(raFreq(IF))
          wvnmhi = DBLE(raFreq(IF)+kaFrStep(iTag))
! ++++++++++++++++ this is to include MIE data from RTSPEC +++++++++++++++++
          
          nlyr = nlev-1
          DO iL=1,kProfLayer
            ssalb(iL) = DBLE(0.0)
            asym(iL)  = DBLE(0.0)
            asym_rtspec(iL) = 0.0
          END DO
          
          DO iL=0,maxmom
            DO I = 1,maxcly
              pmom(iL,I) = DBLE(0.0)
            END DO
          END DO
          
! to test no scattering, just replace following doloops with DO N = 1,-1
          
          IF (iRayleigh == -1) THEN     !want cloud, not Rayleigh scattering
            CALL SetUpClouds(nstr,nmuobs(1),iaCloudWithThisAtm,  &
                iaCldTop,iaCldBot,iaCloudNumLayers,raFreq(IF),  &
                iAtm,iaaRadLayer,iNumLayer,  &
                IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB,  &
                TABEXTINCT, TABSSALB, TABASYM, ISCATTAB,  &
                extinct,dtauc,ssalb,asym,pmom)
          ELSE IF (iRayleigh == +1) THEN   !want Rayleigh, not cloud scattering
            CALL SetUpRayleigh(nlev,nstr,nmuobs(1),raFreq(IF),raDensity,  &
                raThicknessRayleigh,dtauc,ssalb,asym,pmom)
          END IF
          
          dTotalOpticalDepth=DBLE(0.0)
          DO iL=1,NLEV-1
            dTotalOpticalDepth = dTotalOpticalDepth+dtauc(iL)
          END DO
          
! ++++++++++++++++ final initializations ++++++++++++++++++++++++++++++
!     note we do not need flux here!!!!
          CALL FinalInitialization(  &
              iDownWard,rSatAngle,raTopIntensity(IF),raSolarBeam(IF),  &
              raUseEmissivity(IF),rSurfaceTemp,  &
              dtauc,dTotalOpticalDepth,iDoFlux,nlyr+1,iNp,iaOp,  &
              usrtau,ntau,utau,usrang,numu,umu, nphi,phi,fisot,fbeam,umu0,phi0,  &
              ibcnd,lamber,albedo,  &
              btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)
          
!!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
          nmom = MAX(2*nmuobs(1) + 1,nstr)
          nmom = MIN(maxmom,nmom)
          
          CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,  &
              WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,  &
              UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,  &
          FISOT, LAMBER, ALBEDO  BTEMP, TTEMP, TEMIS,  &
                PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN,  &
                FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )
            
            DO iii = 1,iNp
              raKStep(iStepPts)     = raCorrelatedK(IF)
              raFreqSTEP(iStepPts)  = raFreq(IF)
              raaIntenSTEP(iii,iStepPts) = SNGL(uu(1,iii,1))
            END DO
            
          END DO   !!loop over freqs
! ------------------------ END OF DISORT ---------------------------------
          
!interpolate coarse step grid onto fine kCARTA grid
          DO iii = 1,iNp
            DO IF = 1,iStepPts
              raIntenStep(IF)     = raaIntenStep(iii,IF)
              raNoScatterStep(IF) = raaNoScatterStep(iii,IF)
            END DO
            CALL Interpolator(raFreqStep,rakStep,raIntenStep,  &
                raNoScatterSTEP,raaNoScatterStep,iii, iStepPts,iDownWard,nlev,  &
                iProfileLayers,raPressLevels,  &
                raCorrelatedK,raLayerTemp,absprof,raFreq,  &
                rSatAngle,rSolarAngle, rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                raInten)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
!        iF = 1
!        print *,iii,iF,raIntenStep(iF),raNoScatterStep(iF),raInten(iF)
          END DO
          
          9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm
          
          900  FORMAT(8(E12.4,'  '))
          999  FORMAT(I2,'  ',I6,'   ',F9.4,'  ',F8.6,'  ',F8.4,'  ',3(E12.6,'  '))
          
          RETURN
        END SUBROUTINE interface_disort
        
!************************************************************************
! this subroutine figures out how to do the interpolations
        
        SUBROUTINE Interpolator(raFreqStep,rakStep,raIntenStep,  &
            raNoScatterSTEP,raaNoScatterStep,iii, iStepPts,iDownWard,nlev,  &
            iProfileLayers,raPressLevels,  &
            raCorrelatedK,raLayerTemp,absprof,raFreq, rSatAngle,rSolarAngle,  &
            rSurfaceTemp,rSurfPress,raUseEmissivity, raInten)
        
        
        REAL, INTENT(IN OUT)                     :: raFreqStep(kMaxPts)
        NO TYPE, INTENT(IN OUT)                  :: rakStep
        NO TYPE, INTENT(IN OUT)                  :: raIntenSte
        NO TYPE, INTENT(IN OUT)                  :: raNoScatte
        NO TYPE, INTENT(IN OUT)                  :: raaNoScatt
        INTEGER, INTENT(IN OUT)                  :: iii
        INTEGER, INTENT(IN OUT)                  :: iStepPts
        INTEGER, INTENT(IN OUT)                  :: iDownWard
        INTEGER, INTENT(IN OUT)                  :: nlev
        NO TYPE, INTENT(IN OUT)                  :: iProfileLa
        NO TYPE, INTENT(IN OUT)                  :: raPressLev
        NO TYPE, INTENT(IN OUT)                  :: raCorrelat
        NO TYPE, INTENT(IN OUT)                  :: raLayerTem
        NO TYPE, INTENT(IN)                      :: absprof
        REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
        REAL, INTENT(IN OUT)                     :: rSatAngle
        NO TYPE, INTENT(IN OUT)                  :: rSolarAngl
        NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
        REAL, INTENT(IN OUT)                     :: rSurfPress
        NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
        REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
        IMPLICIT NONE
        
        INCLUDE '../INCLUDE/scatterparam.f90'
        
! inputs
        REAL :: raPressLevels(kProfLayer+1)
        REAL :: raKstep(kMaxPts),raIntenStep(kMaxPts)
        REAL :: raNoScatterSTEP(kMaxPts),raaNoScatterStep(kProfLayer,kMaxPts)
        REAL :: raCorrelatedK(kMaxPts), raLayerTemp(kProfLayer)
        REAL :: ABSPROF(MAXNZ,MAXABSNU), rSolarAngle,rSurfaceTemp
        REAL :: raUseEmissivity(kMaxPts)
        INTEGER :: iProfileLayers
! output
        
        
! local variables
        REAL :: DWA(kMaxPts),raTau(kProfLayer),raNoScatterAll(kMaxPts)
        INTEGER :: Indx(kMaxPts),IF,iL
        REAL :: radtot,ttorad
        
        IF (kScatter == 1) THEN
!!!interpolate wavenumber points
          IF (iDownWard == -1) THEN      !!!works better for uplook instr
!          DO iF = 1,kMaxPts
!            DO iL = 1,NLEV-1
!              raTau(iL)  = absprof(iL,iF)
!            END DO
!            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
!     $         rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,
!     $         raUseEmissivity(iF),raFreq(iF),raNoScatterAll(iF),-1,
!     $         iProfileLayers,raPressLevels)
!          END DO
            DO IF = 1,kMaxPts
              raNoScatterAll(IF) = raaNoScatterStep(iii,IF)
            END DO
            CALL INTERP_PLANCK_0(  &
                raFreqStep,raIntenStep,raNoScatterStep,iStepPts,  &
                raFreq,raNoScatterAll,raInten)
          ELSE IF (iDownWard == 1) THEN   !!!works better for downlook instr
            DO IF = 1,kMaxPts
              DO iL = 1,NLEV-1
                raTau(iL)  = absprof(iL,IF)
              END DO
              CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,  &
                  rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,  &
                  raUseEmissivity(IF),raFreq(IF),raNoScatterAll(IF),-1,  &
                  iProfileLayers,raPressLevels)
            END DO
            CALL INTERP_PLANCK_0(  &
                raFreqStep,raIntenStep,raNoScatterStep,iStepPts,  &
                raFreq,raNoScatterAll,raInten)
          END IF
        END IF
        
        IF ((kScatter >= 2) .AND. (iDownWard == -1)) THEN
!!!first interpolate in k
          DO IF = 1,kMaxPts
            CALL rLINEAR_ONE(raKSTEP,raIntenStep,iStepPts,  &
                raCorrelatedK(IF),raInten(IF))
          END DO
          
        ELSE IF ((kScatter >= 2) .AND. (iDownWard == 1)) THEN
!!!interpolate in k
!!!do cumulative distr function to see if things might be too spiky
          DO IF=1,kMaxPts
            DO iL=1,NLEV-1
              raTau(iL)  = absprof(iL,IF)
            END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,  &
                rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,  &
                raUseEmissivity(IF),raFreq(IF),raNoScatterAll(IF),-1,  &
                iProfileLayers,raPressLevels)
          END DO
          CALL INTERP_PLANCK_3(  &
              raFreqStep,raKStep,raIntenStep,raNoScatterStep,iStepPts,  &
              raFreq,raNoScatterAll,raCorrelatedK,raInten)
          
        END IF
        
        RETURN
      END SUBROUTINE Interpolator
!************************************************************************
! this subroutine takes in the input abs coeffs (raaAbsCoeff) where raa(1,:)
! is the lowest layer and raa(kProfLayer,:) is the highest layer .. it then
! outputs these  abs coeffs into absprof, where absprof(1,:) is the top, and
! absprof(iNumLayer,:) = ground
      
! remember kCARTA has -----------------
!                      layer 100 = TOA
!                     -----------------
      
!                           ...
      
!                     -----------------
!                      layer 1 = GND
!                     -----------------
! for a DOWNLOOK instr, the layering is set as 1,2,3,4,5,......,100
!   rFracTop set for the top layer (100) and rFracBot set for the bot layer (1)
! for a   UPLOOK instr, the layering is set as 100,99,98,........,1
!   rFracTop set for the top layer (100) and rFracBot set for the bot layer (1)
! and this is WHAT IS ALWAYS set prior to the rad transfer computations
      
! remember DISORT has -----------------
!      and RTSPEC      layer 1 = TOA
!                     -----------------
      
!                           ...
      
!                     -----------------
!                      layer 100 = GND
!                     -----------------
! so for an downlook instr, everything is easily coded below (ie simply
! flip  1 -->  100, 100 --> 1,
! flip  2 -->   99,  99 --> 2,
! flip  3 -->   98,  98 --> 3 etc)
      
      
!     sets optical depths for NLEV-1 layers and NABSNU wavenumbers.
!     The temperature (K) of the profile is also returned. (no height needed)
      
! same as GetAbsProfile in scatter_rtspec_main,f except that
! 1) calls getverticaltempDISORT
! 2) sets the particle number denisties of the layers
      
      SUBROUTINE GetAbsProfileDISORT(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
          iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
          NABSNU, NLEV, TEMP, ABSPROF,  &
          ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity,  &
          raDensity,raLayerTemp,  &
          iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)
      
      
      REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
      REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
      INTEGER, INTENT(IN)                      :: iNumLayer
      NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
      INTEGER, INTENT(IN OUT)                  :: iAtm
      INTEGER, INTENT(OUT)                     :: iNpmix
      REAL, INTENT(IN)                         :: rFracTop
      REAL, INTENT(IN)                         :: rFracBot
      REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
      NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
      REAL, INTENT(IN OUT)                     :: rSurfPress
      INTEGER, INTENT(IN OUT)                  :: NABSNU
      NO TYPE, INTENT(IN OUT)                  :: NLEV
      REAL, INTENT(IN OUT)                     :: TEMP(*)
      REAL, INTENT(IN OUT)                     :: ABSPROF(MAXNZ,*)
      INTEGER, INTENT(IN OUT)                  :: ICLDTOP
      INTEGER, INTENT(IN OUT)                  :: iCLDBOT
      INTEGER, INTENT(IN OUT)                  :: IOBS
      NO TYPE, INTENT(IN OUT)                  :: iDownward
      REAL, INTENT(IN OUT)                     :: iwp
      NO TYPE, INTENT(IN OUT)                  :: raNumberDe
      REAL, INTENT(OUT)                        :: raDensity(*)
      NO TYPE, INTENT(IN OUT)                  :: raLayerTem
      NO TYPE, INTENT(IN OUT)                  :: iProfileLa
      NO TYPE, INTENT(IN OUT)                  :: raPressLev
      NO TYPE, INTENT(IN OUT)                  :: raThicknes
      NO TYPE, INTENT(IN OUT)                  :: raThicknes
      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/scatterparam.f90'
      
! these are variables that come in from kcartamain.f
      
      REAL :: rSurfaceTemp
      REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER :: iDownWard,iProfileLayers
      REAL :: raNumberDensity(kProfLayer),raLayerTemp(kProfLayer)
! these are variables that we have to set
      INTEGER :: NLEV         !!!!!!!!!!!!!MAXNZ, MAXABSNU
      
      
      
      REAL :: raThicknessRayleigh(kProfLayer)
      
! local variables
      INTEGER :: iaRadLayer(kProfLayer), iaTemp(kProfLayer),iFr, iL, iLay, iiDiv
      REAL :: NU, raVT1(kMixFilRows), InterpTemp, InterpTempSurf
      
! these are to flip the temperature, abs profiles if instr looks up
      REAL :: raTemp(kProfLayer+1),raaTemp(kProfLayer,kMaxPts),rT2
      
      nabsnu = kMaxPts
      nlev=iNumLayer+1           !this is the number of pressure levels
      
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
          WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
          WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
      END DO
      
      iiDiv = 0
      555  CONTINUE
      IF (iiDiv*kProfLayer < iaRadLayer(3)) THEN
        iiDiv = iiDiv + 1
        GO TO 555
      END IF
      iiDiv = iiDiv - 1
      
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
      
! set the temperature of the bottommost layer correctly
      IF (iDownWard == 1) THEN         !downlook instr
! if the bottommost layer is fractional, interpolate!!!!!!
        iL=iaRadLayer(1)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,  &
            1,iL)
        rT2=InterpTempSurf(iProfileLayers,raPressLevels,raVTemp,rFracBot,  &
            1,iL,rSurfaceTemp,rSurfPress)
        WRITE(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
        iL=iaRadLayer(iNumLayer)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,  &
            -1,iL)
        WRITE(kStdWarn,*) 'top temp interped to ',raVT1(iL)
      ELSE IF (iDownWard == -1) THEN       !uplook instr
! if the bottom layer is fractional, interpolate!!!!!!
        iL=iaRadLayer(iNumLayer)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,  &
            1,iL)
        WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
        iL=iaRadLayer(1)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,  &
            -1,iL)
        WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)
      END IF
      
! now set the kCARTA LAYERS temperatures, used with NoScatterRadTransfer
! recall for DISORT , toa = layer 1   while for kCARTA, toa = 100
! recall for DISORT , gnd = layer 100 while for kCARTA, gnd = 1
! but also that for UPLOOK instr, clear sky kCARTA flips the layering!!!!
      IF (iDownward == 1) THEN
        DO iLay=1,iNumLayer
          raLayerTemp(iLay) = raVT1(iaRadLayer(iNumLayer-iLay+1))
        END DO
      ELSE
        DO iLay=1,iNumLayer
          raLayerTemp(iLay) = raVT1(iaRadLayer(iNumLayer-iLay+1))
        END DO
        DO iLay=1,iNumLayer
          raVT1(iLay) = raLayerTemp(iNumLayer-iLay+1)
        END DO
        DO iLay=1,iNumLayer
          raLayerTemp(iLay) = raVT1(iLay)
        END DO
      END IF
      
!because kCARTA flips everything for uplook instrument early on in
!the code, we have to reflip everything for DISORT
      IF (iDownWard == -1) THEN
        DO iLay = 1, iNumLayer
          iaTemp(iLay) = iaRadLayer(iLay)
        END DO
        DO iLay = 1, iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaTemp(iLay)
        END DO
      END IF
      
! set the vertical temperatures of the atmosphere
! instead of temperatures at n layers, it computes temps at (n+1) levels
      CALL SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,  &
          rSurfaceTemp,iProfileLayers,raPressLevels)
      
! now set up the abs coeffs
! initialize array to all zeroes
      DO iFr=1,kMaxPts
        DO iLay=iNumLayer+1,kProfLayer
          absprof(iLay,iFr) = 0.0
          raDensity(iLay) = 0.0
          raThicknessRayleigh(iLay) = 0.0
        END DO
      END DO
      
      DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
        nu=1.0
        IF (iLay == 1) THEN
          nu=rFracBot
        ELSE IF (iLay == iNumLayer) THEN
          nu=rFracTop
        END IF
        raDensity(iNumLayer-iLay+1) = raNumberDensity(iL-iiDiv*kProfLayer)
!!!!!Laythick is from outlayers.param    before Nov 17,2001
!!!!!Laythick is gotten from the KLAYERS profile after Nov 17, 2001
!!!!!only need for Rayleigh scattering, which is small in infrared
        raThicknessRayleigh(iNumLayer-iLay+1) =  &
            raThickness(iL-iiDiv*kProfLayer)*nu
        DO iFr=1,kMaxPts
!absprof wants level 1 == TOA, level iNumLayer= gnd
          absprof(iNumLayer-iLay+1,iFr) = raaAbs(iFr,iL)*nu
        END DO
      END DO
      
! comment this out on Oct 20, 2000 as is taken care of in interface_disort
! now set iobs
!      iobs=(iNumLayer+1)-iobs+1
      
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
! remember there is ONE more level than there are layers
!      icldtop=(iNumLayer+1)-icldtop+1
!      icldbot=(iNumLayer+1)-icldbot+1
      
      RETURN
    END SUBROUTINE GetAbsProfileDISORT
!************************************************************************
    
! set the vertical temperatures of the atmosphere
! this sets the temperatures at the pressure level boundaries, using the
! temperatures of the pressure layers that have been supplied by kLayers
! note that temeprature of bottom level is NOT set to rSurfaceTemp
    
    SUBROUTINE SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,  &
        rSurfaceTemp,iProfileLayers,raPressLevels)
    
    
    REAL, INTENT(IN OUT)                     :: TEMP(*)
    INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
    REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
    INTEGER, INTENT(IN)                      :: iNumLayer
    NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
    NO TYPE, INTENT(IN OUT)                  :: iProfileLa
    NO TYPE, INTENT(IN OUT)                  :: raPressLev
    IMPLICIT NONE
    
    INCLUDE '../INCLUDE/scatterparam.f90'
    
! these are variables that come in from kcartamain.f
    REAL :: rSurfaceTemp,raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! these are variables that we have to set
    
    
! local variables
    INTEGER :: iL, iLay ,iM, idiv
    REAL :: FindBottomTemp,Temp1(maxnz)
    REAL :: pavg(kProfLayer),rP,raProfileTemp(kProfLayer)
    
    DO iLay=1,MAXNZ
      Temp1(iLay) = 0.0
    END DO
    
! see which set of Mixed Paths the current atmosphere occupies eg
! set 1 = 1..100, set2= 101..200 etc
! eg if current atmosphere is from MixfilPath 110-190, and kProfLayer = 100,
! then we set iMod as 2      idiv(150,100) = 1  === 2nd set of mixed paths
! assume each atmosphere has at least 25 layers in it!!!
    iM=idiv(iaRadLayer(25),kProfLayer)+1
    DO iLay=1,kProfLayer
      raProfileTemp(iLay) = raVTemp(iLay+(iM-1)*kProfLayer)
    END DO
    
    DO iLay=1,iNumLayer
      iL=iaRadLayer(iLay)
!map this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
      iL=iL-idiv(iL,kProfLayer)*kProfLayer
      IF (iL == 0) THEN
        iL = kProfLayer
      END IF
      rP=raPressLevels(iL+1)-10000*delta
      IF (rp < raPressLevels(kProfLayer+1)) THEN
        rp = raPressLevels(kProfLayer+1)+10000*delta
      END IF
      TEMP1(iNumLayer-iLay+1) = FindBottomTemp(rP,raProfileTemp,  &
          raPressLevels,iProfileLayers)
    END DO
    
!ccc this is where we differ from RTSPEC
!ccc      TEMP1(iNumLayer+1) = rSurfaceTemp
    rP = DISORTsurfPress
    TEMP1(iNumLayer+1) = FindBottomTemp(rP,raProfileTemp,  &
        raPressLevels,iProfileLayers)
    
    DO iLay=1,iNumLayer+1
      temp(iLay) = temp1(iLay)
    END DO
    
    RETURN
  END SUBROUTINE SetDISORTTemp
  
!************************************************************************
! this subroutine is from LBLRTMv5.10
  
  REAL FUNCTION Rayleigh(fr,wtot,layerthick)
  
  
  REAL, INTENT(IN)                         :: fr
  REAL, INTENT(IN OUT)                     :: wtot
  REAL, INTENT(IN)                         :: layerthick
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
!     ********** Rayleigh Scattering calculation **********
  
!     The formulation, adopted from MODTRAN_3.5 (using approximation
!     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization
!     = 0.0279, output in km-1 for T=273K & P=1 ATM) is as follows:
  
!     The rayleigh "absorption" coefficient (scattering in the direct
!     beam) ray_abs can be defined as
  
!         ray_abs = (fr**4/(9.38076E18-1.08426E09*fr**2))
!     *        *wmol_tot*conv_cm2mol
  
!     where fr is the wavenumber value, wmol_tot is the total
!     molecular amount in the layer, and conv_cm2mol is the conversion
!     factor derived by multiplying air density (2.68675E19 mol/cm3)
!     at 273 K with the number of km per cm (1.e-5 km/cm).
  
!     For numerical purposes, the layer amount of all molecules is
!     calculated above as WTOT, which has been scaled by 1.e-20. We
!     have taken 1.e20 out of the air density in the denominator
!     as well. In addition, a factor of 10,000 (cm-1) has been
!     divided out of fr. Finally, the radiation field is
!     excluded, so xfr**4 is replaced by xfr**3. When
!     JRAD=1, the radiation field is put in by multiplying the
!     absorption by xfr.
  
!     Rayleigh scattering in the direct beam is only calculated for
!     model runs > 3100 cm-1.
  
!        Thus the current formulation is
  
  REAL :: Losch             !Loschmidt number
  
  REAL :: xrayl
!wtot is total number of molecules
!in layer and xrayl is multiplication factor
  
  REAL :: conv_cm2mol,xfr
  REAL :: l,l1,a0,a1,a2,a3,y0,y1,y2,y3
  
  xrayl = 1.0
  conv_cm2mol = 1./(2.68675E-1*1.e5)
  
  l=10000.0/fr        !change from cm-1 to um
  
!      xfr = fr/1.e4
!      y0 = (xfr**4/(9.38076E2-10.8426*xfr**2))*wtot*conv_cm2mol*1e-20
!      y0 = y0*100              !change abs coeff from cm-1 to m-1
!      y0 = y0*layerthick
!c LBLRTM had xfr**3 and not xfr**4
!c      y0 = (xfr**3/(9.38076E2-10.8426*xfr**2))*wtot*conv_cm2mol*XRAYL
  
  
! works quite well in the 300-3000 cm-1 range! use this!
! this is from SBDART
! calculate molecular rayleigh scattering coefficient
! using approximation of shettle et al., 1980 (appl. opt., 2873-4)
! with the depolarization = 0.0279 instead of 0.035
! input: fr = frequency in wavenumbers (cm-1)
! output: y0 = molecular scattering coefficient (km-1)
!               for temperature = 273 k & pressure = 1 atm.
  Losch = (kAtm2mb*100/kBoltzmann/273 * 1E-6)    ! (p0/k T0)
  y1 = fr**4/(9.38076E+18 - 1.08426E+09 * fr ** 2) * (wtot/Losch)
  y1 = y1*layerthick/1000  !the abs coeff above is in km-1
  
! this is from Radiative Transfer in the Ocean, Stamnes&Thomas pg 73
!      a0=3.9729066
!      a1=4.6547659e-2
!      a2=4.5055995e-4
!      a3=2.3229848e-5
!      y0=0.0
!      y0=a0 + a1/(l*l) + a2/(l*l*l*l) + a3/(l*l*l*l*l*l)
!      y0=1e-28/(l*l*l*l)*y0              !!!cross section in cm2
!      y0 = y0*(500*100)                  !!!ty0pical lay0er ~ 500 m
!      y0=y0*wtot                         !!! cm3 * cm-3 = optical depth
  
! try this one, from same book
!      a0=1.000292        !!ref index of air
!      !change wavelength to m (1um = 1e-6 m)
!      !change wtot to m-3 (by multiplying by 1e6)
!      y1=32*(kPi**3)*((a0-1)**2)/(3*l*l*l*l*1e-24*wtot*1e6)
!      y1=y1*250             !!typicl layer ~ 250m
  
! try this one, from same book
!      y2=8*kpi/3*((2*kpi/l*1e-6)**4)*(a0*a0-1)/(a0*a0+2)*1e-9*(wtot*1e6)
  
! try this from Applied Optics
!      y3 = 8*(kPi**3)/3*(wtot*1e6)*((a0*a0-1)/(a0*a0+2)**2)*(l*1e-6)
  
! try this from Liou
!      !!!need wavelength in um to compute refractive index (eqn 3.70)
!      a0 = 1/(l*l)
!      a0 = 6452.8 + 2949810/(146-a0) + 25540/(41-a0)
!      a0 = a0*1.0e-8 + 1       !!this is the ref index
  
!      a1 = 0.035               !!this is the anisotropy factor
  
!      !!!now need wavelength in cm
!      l1 = l*1.0e-4
!      y3 = 8*(kPi**3)/3/(wtot)*((a0*a0-1)**2) *(6+3*a1)/(6-7*a1)/(l1**4)*25000
  
!      !!!need wavelength in cm to compute scattering cross section (eqn 3.72)
!      y3 = 8*(kPi**3)/3/(l1**4)
!      y3 = y3 * ((a0*a0-1)**2)
!      y3 = y3 * (6+3*a1)/(6-7*a1)
!      Losch = (kAtm2mb*100/kBoltzmann/273 * 1e-6)    ! (p0/k T0)
!      y3 = y3/Losch/Losch       !!!!!!!ONLY WAY TO MAKE THINGS WORK
!      y2 = y3   !!!!!!this is eqn 3.72 which is scattering PER molecule
  
!      y3 = y3*wtot*(layerthick*100)  !!each layer ~ 250 m thick ~ 25000 cm
!      print *,fr,wtot,y1,y3
  
  rayleigh = y1
  
  RETURN
END FUNCTION Rayleigh


!************************************************************************
! this subroutine finds out the index of wavenumber points at which to do
! the radiative transfer

SUBROUTINE FindWavenumberPoints(iStep,iScatter,raCorrelatedK, iaStep,iFFMax)


NO TYPE, INTENT(IN)                      :: iStep
NO TYPE, INTENT(IN)                      :: iScatter
NO TYPE, INTENT(IN OUT)                  :: raCorrelat
INTEGER, INTENT(OUT)                     :: iaStep(kMaxPts)
INTEGER, INTENT(OUT)                     :: iFFMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: iStep               !tells how many points to step thru, or find
INTEGER :: iScatter            !1,2 or 3 tells which type to use
INTEGER :: !gives the index of wavenumber points to use

REAL :: raCorrelatedK(kMaxPts) !array of abs coeffs of layer closest to gnd

INTEGER :: IF,iL

DO IF=1,kMaxPts
  iaStep(IF) = 0
END DO

IF ((iScatter < 1) .OR. (iScatter > 3)) THEN
  WRITE (kStdErr,*) 'For DISORT, need to choose kScatter = 1,2 or 3'
  WRITE (kStdErr,*) 'iScatter = ',iScatter
  WRITE (kStdErr,*) 'Please retry'
  CALL DoStop
END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
IF (iScatter == 1) THEN
!!!!! use equal spacing of points; nothing special about chosen points
!!!!! do radiative transfer on freq index pts 1,1+iStep,1+2*iStep,...
!!!!!       1+N*iStep,10000
  iL=1
  IF=1
  10     CONTINUE
  iaStep(iL) = IF
  iL=iL+1
  IF=IF+iStep
  IF (IF <= kMaxPts) THEN
    GO TO 10
  ELSE
    GO TO 20
  END IF
  20     CONTINUE
  IF (iaStep(iL-1) /= kMaxPts) THEN
!make sure last point radiances are computed for, is the 10000th pt
    iaStep(iL) = kMaxPts
    iFFMax = iL
  ELSE
    iL=iL-1
    iFFMax = iL
  END IF
END IF            !IF (kScatter .EQ. 1) THEN
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

IF (iScatter == 2) THEN
!!!!! do radiative transfer on freq index pts a,b,c,...z
!!!!! where the points are chosen such that they have the smallest
!!!!! absorption coeffs
  CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
END IF            !IF (kScatter .EQ. 2) THEN

IF (iScatter == 3) THEN
!!!!! do radiative transfer on freq index pts a,b,c,...z
!!!!! where the points are chosen such that they range from the
!!!!! smallest to largest absorption coeffs
  CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
END IF            !IF (kScatter .EQ. 3) THEN

RETURN
END SUBROUTINE FindWavenumberPoints
!************************************************************************
! this subroutine sorts the k values and returns index values
! if iScatter = 2, sort the k values from smallest to largest
!                  return indices of smallest iStep values
! if iScatter = 3, sort the k values from smallest to largest
!                  return indices of smallest to largest k, in steps of iStep
! the outputs are iaStep and iFFMax

SUBROUTINE FindIndices(iScatter,iStep,arr,iaStep,iFFMax)


NO TYPE, INTENT(IN OUT)                  :: iScatter
NO TYPE, INTENT(IN)                      :: iStep
REAL, INTENT(IN)                         :: arr(kMaxPts)
INTEGER, INTENT(OUT)                     :: iaStep(kMaxPts)
INTEGER, INTENT(OUT)                     :: iFFMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: iStep               !tells how many points to step thru, or find
INTEGER :: iScatter            !1,2 or 3 tells which type to use
INTEGER :: !gives the index of wavenumber points to use

REAL :: !array of abs coeffs of layer closest to gnd

INTEGER :: indx(kMaxPts),i,j,IF,IL
REAL :: TempArr(kMaxPts)

CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

! now we have sorted the k so that indx(i) contains the sorted keys info
IF (iScatter == 2) THEN
!simply give smallest k values for all but last two points.
!for last point, give largest k value
!for second last point, give medium k value
  DO i=1,iStep-2
    iaStep(i) = indx(i)
  END DO
  j = kMaxPts-i
  IF (MOD(j,2) == 0) THEN
    j=j/2 + i
  ELSE
    j=(j+1)/2 + i
  END IF
!        iaStep(iStep-1) = indx(kMaxPts-100)
  iaStep(iStep-1) = indx(j)
  iaStep(iStep)   = indx(kMaxPts)
  iFFMax = iStep
END IF

IF (iScatter == 3) THEN  !give smallest to largest k values
  iL=1
  IF=1
  10     CONTINUE
  iaStep(iL) = indx(IF)
  iL=iL+1
  IF=IF+iStep
  IF (IF <= kMaxPts) THEN
    GO TO 10
  ELSE
    GO TO 20
  END IF
  20     CONTINUE
  IF (iaStep(iL-1) /= indx(kMaxPts)) THEN
!make sure radiances are computed for largest k
    iaStep(iL) = indx(kMaxPts)
    iFFMax = iL
  ELSE
    iL=iL-1
    iFFMax = iL
  END IF
END IF

! very important .. want to keep integrity of the Planck fcn!!!!!
! so make sure the original 1st and 10000th points are always used
IF (iScatter >= 2) THEN
  
!!!!!!make sure (iaStep(i) contains "1")
  j = -1
  iL = 1
  30      CONTINUE
  IF (iaStep(iL) == 1)  THEN
    j = +1      !!!!!yes, it does contain "1"
  ELSE
    iL = iL + 1 !!!!keep on hoping
  END IF
  IF ((iL <= iFFMax) .AND. (j < 0)) THEN
    GO TO 30
  END IF
  IF (j < 0) THEN !!!!!!poooey, not found
    i = -1
    iL = 1
    40       CONTINUE
    IF (indx(iL) == 1) THEN
      i = +1       !!!!!yes, it is found
      iFFMax = iFFMax + 1
      iaStep(iFFMax) = indx(iL)
    ELSE
      iL = iL + 1
      GO TO 40
    END IF
  END IF
  
!!!!!!make sure (iaStep(i) contains "kMaxPts")
  j = -1
  iL = 1
  50      CONTINUE
  IF (iaStep(iL) == kMaxPts)  THEN
    j = +1      !!!!!yes, it does contain "kMaxPts"
  ELSE
    iL = iL + 1 !!!!keep on hoping
  END IF
  IF ((iL <= iFFMax) .AND. (j < 0)) THEN
    GO TO 50
  END IF
  IF (j < 0) THEN !!!!!!poooey, not found
    i = -1
    iL = 1
    60       CONTINUE
    IF (indx(iL) == kMaxPts) THEN
      i = +1       !!!!!yes, it is found
      iFFMax = iFFMax + 1
      iaStep(iFFMax) = indx(iL)
    ELSE
      iL = iL + 1
      GO TO 60
    END IF
  END IF
  
! reset iaStep(iFFMax)
  DO i=1,kMaxPts
    TempArr(i) = 1.0E10
  END DO
  DO i=1,iFFMax
    TempArr(iaStep(i)) = arr(iaStep(i))
  END DO
  CALL NumericalRecipesIndexer(indx,TempArr,kMaxPts)
  DO i=1,iFFMax
    iaStep(i) = indx(i)
  END DO
END IF

RETURN
END SUBROUTINE FindIndices

!************************************************************************
! this subroutine sorts the k values and returns index values
! if iScatter = 2, sort the k values from smallest to largest
!                  return indices of smallest iStep values
! if iScatter = 3, sort the k values from smallest to largest
!                  return indices of smallest to largest k, in steps of iStep
! the outputs are iaStep and iFFMax

SUBROUTINE FindIndices_NotUsed(iScatter,iStep,arr,iaStep,iFFMax)


NO TYPE, INTENT(IN OUT)                  :: iScatter
NO TYPE, INTENT(IN)                      :: iStep
REAL, INTENT(IN)                         :: arr(kMaxPts)
INTEGER, INTENT(OUT)                     :: iaStep(kMaxPts)
INTEGER, INTENT(OUT)                     :: iFFMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: iStep               !tells how many points to step thru, or find
INTEGER :: iScatter            !1,2 or 3 tells which type to use
INTEGER :: !gives the index of wavenumber points to use

REAL :: !array of abs coeffs of layer closest to gnd

INTEGER :: indx(kMaxPts)
REAL :: TempArr(kMaxPts),kmin,kmax,k0,dk

INTEGER :: i,j,IF,iL

CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

! now we have sorted the k so that indx(i) contains the sorted keys info
!      DO i=1,kMaxPts
!        print *,i,indx(i),arr(indx(i))
!      END DO

IF (iScatter == 2) THEN
!simply give smallest k values for all but last two points.
!for last point, give largest k value
!for second last point, give medium k value
  DO i=1,iStep-2
    iaStep(i) = indx(i)
  END DO
  j = kMaxPts-i
  IF (MOD(j,2) == 0) THEN
    j=j/2 + i
  ELSE
    j=(j+1)/2 + i
  END IF
!        iaStep(iStep-1) = indx(kMaxPts-100)
  iaStep(iStep-1) = indx(j)
  iaStep(iStep)   = indx(kMaxPts)
  iFFMax = iStep
END IF

IF (iScatter == 3) THEN
!smallest to largest k values, equally stepped (kmax-kmin)/(iStep-1)
  kmin = arr(indx(1))
  kmax = arr(indx(kMaxPts))
  dk = (kmax-kmin)/(iStep-1)
  k0 = kmin
  
  k0 = k0 + dk
  iaStep(1) = indx(1)
  IF=2
  iL=2
  
  10     CONTINUE
  IF ((arr(indx(iL)) < k0) .AND. (iL < kMaxPts)) THEN
    iL = iL + 1
    GO TO 10
  ELSE
    iaStep(IF) = indx(iL)
    IF = IF + 1
    iL = iL + 1
    k0 = k0 + dk
    IF (k0 < kmax) THEN
      GO TO 10
    END IF
  END IF
  IF (iaStep(IF-1) /= indx(kMaxPts)) THEN
!make sure radiances are computed for largest k
    iaStep(IF) = indx(kMaxPts)
    iFFMax = IF
  ELSE
    IF=IF-1
    iFFMax = IF
  END IF
END IF

! very important .. want to keep integrity of the Planck fcn!!!!!
! so make sure the original 1st and 10000th points are always used
IF (iScatter >= 2) THEN
  
!!!!!!make sure (iaStep(i) contains "1")
  j = -1
  iL = 1
  30      CONTINUE
  IF (iaStep(iL) == 1)  THEN
    j = +1      !!!!!yes, it does contain "1"
  ELSE
    iL = iL + 1 !!!!keep on hoping
  END IF
  IF (iL <= iFFMax) THEN
    GO TO 30
  END IF
  IF (j < 0) THEN !!!!!!poooey, not found
    i = -1
    iL = 1
    40       CONTINUE
    IF (indx(iL) == 1) THEN
      i = +1       !!!!!yes, it is found
      iFFMax = iFFMax + 1
      iaStep(iFFMax) = indx(iL)
    ELSE
      iL = iL + 1
      GO TO 40
    END IF
  END IF
  
!!!!!!make sure (iaStep(i) contains "kMaxPts")
  j = -1
  iL = 1
  50      CONTINUE
  IF (iaStep(iL) == kMaxPts)  THEN
    j = +1      !!!!!yes, it does contain "kMaxPts"
  ELSE
    iL = iL + 1 !!!!keep on hoping
  END IF
  IF (iL <= iFFMax) THEN
    GO TO 50
  END IF
  IF (j < 0) THEN !!!!!!poooey, not found
    i = -1
    iL = 1
    60       CONTINUE
    IF (indx(iL) == kMaxPts) THEN
      i = +1       !!!!!yes, it is found
      iFFMax = iFFMax + 1
      iaStep(iFFMax) = indx(iL)
    ELSE
      iL = iL + 1
      GO TO 60
    END IF
  END IF
  
! reset iaStep(iFFMax)
  DO i=1,kMaxPts
    TempArr(i) = 1.0E10
  END DO
  DO i=1,iFFMax
    TempArr(iaStep(i)) = arr(iaStep(i))
  END DO
  CALL NumericalRecipesIndexer(indx,TempArr,kMaxPts)
  DO i=1,iFFMax
    iaStep(i) = indx(i)
  END DO
END IF

RETURN
END SUBROUTINE FindIndices_NotUsed
!************************************************************************

SUBROUTINE INTERP_PLANCK_0(XA,YA,CA,N,raFO,raNS,raOut)

! this does a very smart interpolation ahem ahem
!        CALL INTERP_PLANCK_0(
!     $          raFreqStep,raIntenStep,raNoScatterStep,iStepPts,
!     $          raFreq,raNoScatterAll,raInten)


REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
REAL, INTENT(IN)                         :: CA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: raFO(*)
REAL, INTENT(IN)                         :: raNS(*)
REAL, INTENT(OUT)                        :: raOut(*)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      XA  : I  : REAL arr : freq array   array(N) in increasing order
!      YA  : I  : REAL arr : intensity y  array(N) from scattering
!      CA  : I  : REAL arr : intensity y  array(N) from clear sky
!      N   : I  : INT      : number of points fom DISORT, in array

!      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
!      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
!      raOut : O  : REAL arr : output radiance             array (kMaxPts)





!      Local Variables
INTEGER :: K,KLO,KHI,iI,iMethod,iDefault
REAL :: A,B,H,radtot,ttorad
REAL :: yn,y0,X,y,w,dx

iMethod = 2        !! uses rads; seems to be a good method, as it
!! will not barf if the radiances are 0.0
!! but it does lead to some noticeable problems
iMethod = 3        !! usesBT(rads); seems to be a good method. needs some
!! tweaks to handle TOA rads for uplook inst (==0)

DO iI=1,kMaxPts
  X = raFO(iI)
  IF (X > xa(n)) X = xa(n) - 1.0E-7
  IF (X < xa(1)) X = xa(1) + 1.0E-7
  
!Determine between which pair of points X falls (bisect loop)
  KLO = 1
  KHI = N
  
  20     IF ( (KHI - KLO) > 1) THEN
    K = (KHI + KLO)/2
    IF (XA(K) > X) THEN
      KHI = k
    ELSE
      KLO = k
    END IF
    GO TO 20
  END IF
  
  dx = XA(KHI) - x
  H  = XA(KHI) - XA(KLO)
  
  IF (H <= 0.0) THEN
    WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI),X,N
    1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
        'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5,  &
        ', XA(KHI) = ',E12.5,',X = ',E12.5,', numpts = ',I4,'. Quitting.')
  END IF
  
  IF (iMethod == 1) THEN
!fix this, as we need to define dh,dl
!B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
!Y=(YA(KHI) + dh) - h*b
    WRITE (kStdErr,*) 'iMethod = 1 invalid in INTERP_PLANCK_0'
    CALL DoStop
    
  ELSE IF (iMethod == 2) THEN
!interpolate radiances
!YA  : I  : REAL arr : intensity y array(N) from scattering
!CA  : I  : REAL arr : intensity y array(N) from clear sky
!!!interpolate in radiance space
    yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
    y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
    B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
    Y = yn - dx*b
    Y = Y + raNS(iI)
    
  ELSE IF (iMethod == 3) THEN
    IF (ya(KHI) > 0 .AND. CA(KHI) > 0 .AND.  ya(KLO) > 0  &
          .AND. CA(KLO) > 0) THEN
!METHOD 3
!interpolate in BT space if rads ~= 0
!!!interpolate in temperature space
!!! get out stuff from clear sky
      yn = radtot(xa(khi),ya(KHI))-radtot(xa(khi),CA(KHI))
      y0 = radtot(xa(klo),ya(Klo))-radtot(xa(klo),CA(Klo))
      B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
      Y = yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
      Y = Y + radtot(x,raNS(iI)) !!add on clear sky stuff
      Y = ttorad(x,y)
    ELSE
!METHOD 3
!interpolate in rad space if rads = 0
      yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
      y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
      B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
      Y = yn - dx*b
      Y = Y + raNS(iI)
    END IF
    
  ELSE IF (iMethod == 4) THEN
!METHOD 4 : directly interpolate in BT space
!interpolate BTs
!fix this
    WRITE (kStdErr,*) 'iMethod = 4 invalid in INTERP_PLANCK_0'
!B = (radtot(WA(KHI),YA(KHI))-radtot(WA(KLO),YA(KLO)))/
!     $        (XA(KHI) - XA(KLO))
!Y=radtot(WA(KHI),YA(KHI))-dx*b !!r(vo+dv,ko+dk)=r(vo,ko)+dr/dko dk
!Y=ttorad(W,Y)+0*dh             !!change bck to rad, add on dr/dvodv
    
  ELSE IF (iMethod == 5) THEN
    WRITE (kStdErr,*) 'iMethod = 5 invalid in INTERP_PLANCK_0'
!Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv
  END IF
  
  raOut(iI) = Y
END DO

RETURN
END SUBROUTINE INTERP_PLANCK_0

!************************************************************************
! this does a linear interpolation, but tries to weight the y axis so that
! information relating to the Planck depenadance on wavenumber, is included

SUBROUTINE INTERP_PLANCK_1(WA,DWA,XA,YA,N,W,X,Y)


REAL, INTENT(IN)                         :: WA(*)
REAL, INTENT(IN)                         :: DWA(*)
REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: W
REAL, INTENT(IN OUT)                     :: X
REAL, INTENT(OUT)                        :: Y
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      WA  : I  : REAL arr : wavenumber w array(N)
!     DWA  : I  : REAL arr : dy/dw so we can weight the interpolations
!      XA  : I  : REAL arr : abs coeff x array(N) in increasing order
!      YA  : I  : REAL arr : intensity y array(N)
!      N   : I  : INT      : number of points in arrays
!      W   : I  : REAL     : wavnumber at which spline is evaluated
!      X   : I  : REAL     : x point at which to evaluate spline
!      Y   : O  : REAL     : y point from spline interpolation
!      -----------------------------------------------------------------

!      Parameters



!      Local Variables
INTEGER :: K,KLO,KHI
REAL :: A,B,H,dh,dl,radtot,ttorad

! can try METHOD 1
! using linear interp (yh-yh)/(xh-xh) = (yh-y)/(xh-x) ==>
! y = yh - (xh-x)(yh-yl)/(xh-xl)
! modify this to
! y = YH - (xh-x)(YH-YL)/(xh-xl)
! where dh,dl take into account the variation of planck wrt wavenumber
! rad(vo+dv) = rad(vo) + d(rad)/dv * dv
! hence YH = yh + dr/dw dh and YL = yl + dr/dw dl

! can try METHOD 2
! r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk + dr/dvo dv

!      -----------------------------------------------------------------

!     Determine between which pair of points X falls (bisect loop)
KLO=1
KHI=N

20   IF ( (KHI - KLO) > 1) THEN
  K=(KHI + KLO)/2
  IF (XA(K) > X) THEN
    KHI = k
  ELSE
    KLO = k
  END IF
  GO TO 20
END IF

H=XA(KHI) - XA(KLO)
dh = W - WA(khi)       !!!wavenumber diff
dl = WA(klo) - W       !!!wavenumber diff
dh = dh * DWA(khi)     !!!first derivative
dl = dl * DWA(klo)     !!!first derivative

IF (H <= 0.0) THEN
  WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
  1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
      'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5,  &
      ', XA(KHI) = ',E12.5,'. Quitting.')
END IF

! METHOD 1
!      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
!      Y=(YA(KHI) + dh) - h*b

! METHOD 2
! interpolate radiances
!      B=(YA(KHI) - YA(KLO))/(XA(KHI) - XA(KLO)) !!slope
!      Y=YA(KHI) - h*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
!      Y=Y+dh            !!add on                            + dr/dvo dv
! interpolate BTs
B=(radtot(WA(KHI),YA(KHI)) - radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
Y=radtot(WA(KHI),YA(KHI)) - h*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
Y=ttorad(W,Y)+0*dh                !!change back to rad, add on  dr/dvo dv

! METHOD 3
!      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv

RETURN
END SUBROUTINE INTERP_PLANCK_1

!************************************************************************

SUBROUTINE INTERP_PLANCK_2(WA,DWA,XA,YA,CA,N, raFO,raNS,raCK,raOut)


REAL, INTENT(IN)                         :: WA(*)
REAL, INTENT(IN)                         :: DWA(*)
REAL, INTENT(IN OUT)                     :: XA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
REAL, INTENT(IN)                         :: CA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: raFO(*)
REAL, INTENT(IN)                         :: raNS(*)
REAL, INTENT(IN)                         :: raCK(*)
REAL, INTENT(OUT)                        :: raOut(*)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      WA  : I  : REAL arr : wavenumber w array(N)
!     DWA  : I  : REAL arr : dy/dw        array(N) to weight the interpolations
!      XA  : I  : REAL arr : corr coeff x array(N) in increasing order
!      YA  : I  : REAL arr : intensity y  array(N) from scattering
!      CA  : I  : REAL arr : intensity y  array(N) from clear sky
!      N   : I  : INT      : number of points fom DISORT, in array

!      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
!      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
!      raCK  : I  : REAL arr : entire correlated k         array (kMaxPts)
!      raOut : O  : REAL arr : output radiance             array (kMaxPts)





!      Local Variables
INTEGER :: K,KLO,KHI,iI
REAL :: A,B,H,dh,dl,radtot,ttorad
REAL :: yn,y0,X,y,w,dx

! can try METHOD 1
! using linear interp (yh-yh)/(xh-xh) = (yh-y)/(xh-x) ==>
! y = yh - (xh-x)(yh-yl)/(xh-xl)
! modify this to
! y = YH - (xh-x)(YH-YL)/(xh-xl)
! where dh,dl take into account the variation of planck wrt wavenumber
! rad(vo+dv) = rad(vo) + d(rad)/dv * dv
! hence YH = yh + dr/dw dh and YL = yl + dr/dw dl

! can try METHOD 2
! r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk + dr/dvo dv


DO iI=1,kMaxPts
  X=raCK(iI)
  W=raFO(iI)
  
!Determine between which pair of points X falls (bisect loop)
  KLO=1
  KHI=N
  
  20     IF ( (KHI - KLO) > 1) THEN
    K=(KHI + KLO)/2
    IF (XA(K) > X) THEN
      KHI = k
    ELSE
      KLO = k
    END IF
    GO TO 20
  END IF
  
  dx = XA(KHI) - x
  H  = XA(KHI) - XA(KLO)
  dh = W - WA(khi)       !!!wavenumber diff
  dl = WA(klo) - W       !!!wavenumber diff
  dh = dh * DWA(khi)     !!!first derivative
  dl = dl * DWA(klo)     !!!first derivative
  
  IF (H <= 0.0) THEN
    WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
    1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
        'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5,  &
        ', XA(KHI) = ',E12.5,'. Quitting.')
  END IF
  
! METHOD 1
!      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
!      Y=(YA(KHI) + dh) - h*b
  
! METHOD 2
! interpolate radiances
!      YA  : I  : REAL arr : intensity y array(N) from scattering
!      CA  : I  : REAL arr : intensity y array(N) from clear sky
!        !!!interpolate in radiance space
!        yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
!        y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
!        B=(yn-y0)/(XA(KHI) - XA(KLO)) !!slope
!        Y=yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
!        Y=Y+raNS(iI) !!add on clear sky stuff
!        !!add on dr/dvo dv
!        !!CANNOT DO THIS as we are modelling differences, not actual numbers
!        !!Y=Y+dh
  
!!!interpolate in temperature space
!!! get out stuff from clear sky
  yn = radtot(wa(khi),ya(KHI))-radtot(wa(khi),CA(KHI))
  y0 = radtot(wa(klo),ya(Klo))-radtot(wa(klo),CA(Klo))
  B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
  Y = yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
  Y = Y + radtot(w,raNS(iI)) !!add on clear sky stuff
  Y = ttorad(w,y)
  
! interpolate BTs
!      B=(radtot(WA(KHI),YA(KHI))-radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
!      Y=radtot(WA(KHI),YA(KHI))-dx*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
!      Y=ttorad(W,Y)+0*dh              !!change back to rad, add on  dr/dvo dv
  
! METHOD 3
!      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv
  
  raOut(iI) = Y
END DO

RETURN
END SUBROUTINE INTERP_PLANCK_2

!************************************************************************
! same as PLANCK2 except we do interpolations in K space

SUBROUTINE INTERP_PLANCK_3(WA,XA,YA,CA,N,raFO,raNS,raCK,raOut)


REAL, INTENT(IN)                         :: WA(*)
REAL, INTENT(OUT)                        :: XA(*)
REAL, INTENT(IN)                         :: YA(*)
REAL, INTENT(IN)                         :: CA(*)
INTEGER, INTENT(IN)                      :: N
REAL, INTENT(IN)                         :: raFO(*)
REAL, INTENT(IN)                         :: raNS(*)
REAL, INTENT(IN)                         :: raCK(*)
REAL, INTENT(OUT)                        :: raOut(*)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      WA  : I  : REAL arr : wavenumber w array(N)
!      XA  : I  : REAL arr : corr coeff x array(N) in increasing order
!      YA  : I  : REAL arr : intensity y  array(N) from scattering
!      CA  : I  : REAL arr : intensity y  array(N) from clear sky
!      N   : I  : INT      : number of points fom DISORT, in array

!      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
!      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
!      raCK  : I  : REAL arr : entire correlated k         array (kMaxPts)
!      raOut : O  : REAL arr : output radiance             array (kMaxPts)





!     Local Variables
INTEGER :: K,KLO,KHI,iI,indx(kMaxPts),iDok,iDoF
REAL :: A,B,H,radtot,ttorad
REAL :: yn,y0,X,Y1,Y2,YTotal,dx
! these have WA,XA,YA,CA sorted according to wavenumber
REAL :: WA_sort(kMaxPts),XA_sort(kMaxPts)
REAL :: YA_sort(kMaxPts),CA_sort(kMaxPts)

1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,  &
    'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5,  &
    ', XA(KHI) = ',E12.5,'. Quitting.')

CALL NumericalRecipesIndexer(indx,WA,N)
DO iI=1,N
  WA_sort(iI) = WA(indx(iI))
  XA_sort(iI) = XA(indx(iI))
  YA_sort(iI) = YA(indx(iI))
  CA_sort(iI) = CA(indx(iI))
END DO

iDoK = -1       !!!!!do not interpolate in K
iDoF = +1       !!!!!do interpolate in F

DO iI=1,kMaxPts
  Y1 = 0.0
  IF (iDoK > 0) THEN
!!!!! -------------  first do dT/dk dk !!!!!!!!!!!!!
    X=raCK(iI)
    
!Determine between which pair of points X falls (bisect loop)
    KLO=1
    KHI=N
    
    20       IF ( (KHI - KLO) > 1) THEN
      K=(KHI + KLO)/2
      IF (XA(K) > X) THEN
        KHI = k
      ELSE
        KLO = k
      END IF
      GO TO 20
    END IF
    
    dx = XA(KHI) - x
    H  = XA(KHI) - XA(KLO)
    
    IF (H <= 0.0) THEN
      WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
    END IF
    
!!!interpolate in temperature space
!!! get out stuff from clear sky
    yn = radtot(wa(khi),ya(KHI))-radtot(wa(khi),CA(KHI))
    y0 = radtot(wa(klo),ya(Klo))-radtot(wa(klo),CA(Klo))
    B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
    Y1 = yn - dx*b   !!this is T(vo+dv,ko+dk) = T(vo,ko) + dT/dko dk
  END IF
  
  Y2=0.0
  IF (iDoF > 0) THEN
!!!!! -------------  then do dT/dv dv !!!!!!!!!!!!!
    X=raFO(iI)
    
!Determine between which pair of points X falls (bisect loop)
    KLO=1
    KHI=N
    
    30       IF ( (KHI - KLO) > 1) THEN
      K=(KHI + KLO)/2
      IF (WA_Sort(K) > X) THEN
        KHI = k
      ELSE
        KLO = k
      END IF
      GO TO 30
    END IF
    
    dx = WA_Sort(KHI) - x
    H  = WA_Sort(KHI) - WA_Sort(KLO)
    
    IF (H <= 0.0) THEN
      WRITE(kStdWarn,1010) KLO,KHI,WA_Sort(KLO),WA_Sort(KHI)
    END IF
    
!!!interpolate in temperature space
!!! get out stuff from clear sky
    yn = radtot(wa_sort(khi),ya_sort(KHI)) - radtot(wa_sort(khi),ca_sort(KHI))
    y0 = radtot(wa_sort(klo),ya_sort(Klo)) - radtot(wa_sort(klo),ca_sort(Klo))
    B = (yn-y0)/(wa_sort(KHI) - wa_sort(KLO)) !!slope
    Y2 = yn - dx*b   !!this is T(vo+dv,ko+dk) = T(vo,ko) + dT/dvo dv
  END IF
  
! remember      iDoK = -1       !!!!!do not interpolate in K, so Y1 == 0.0
! remember      iDoF = +1       !!!!!do interpolate in F    , so Y2 <> 0.0
  
  YTotal = Y1 + Y2        !!this is the total dT/dk dk + dT/dv dv
  YTotal = YTotal + radtot(raFO(iI),raNS(iI)) !!add on clear sky stuff
!!!!!now change from BT to radiance
  YTotal = ttorad(raFO(iI),YTotal)
  raOut(iI) = YTotal
  
END DO

RETURN
END SUBROUTINE INTERP_PLANCK_3

!************************************************************************
! this subrtouine computes drad/dwavenumber for a Planck blackbody radiance

SUBROUTINE drad_dv(WA,DWA,YA,N)


REAL, INTENT(IN OUT)                     :: WA(*)
REAL, INTENT(IN OUT)                     :: DWA(*)
REAL, INTENT(IN OUT)                     :: YA(*)
INTEGER, INTENT(IN)                      :: N
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      WA  : I  : REAL arr : wavenumber w array(N)
!     DWA  : O  : REAL arr : dy/dw so we can weight the interpolations
!      YA  : I  : REAL arr : intensity y array(N)
!      N   : I  : INT      : number of points in arrays
!      -----------------------------------------------------------------

!     Parameters


INTEGER :: i
REAL :: raT(kMaxPts),r1,r2,v,u,dv,du

! we know r(fo+df) = r(fo) + df dr/df
! Let us assume r(fo) = r(fo,T) ===> first we have to find equivalent temp T
! having found this, we then simply diff dr/df where r = radiance at temp T
! so given radiance r(f), do rad2bt(f,r) --> T
!                         then differentiate d(rPlanck(f,T))/df

r1 = SNGL(kPlanck1)  !! need this for rad2bt and dBT/dT
r2 = SNGL(kPlanck2)  !! need this for rad2bt and dBT/dT

! first convert the intensities to equivalent temps T
! bt = c2 * fr / log(1 + c1 * fr^3 / rd)
DO i = 1,N
  raT(i) = r2*wa(i)/(LOG(1+r1*(wa(i)**3)/ya(i)))
END DO

! now do dr/df where r=ttorad = planck radiance at temp T, frequency fr
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
DO i = 1,N
  u  = r1*(wa(i)**3)
  du = r1*3*(wa(i)**2)
  v  = EXP(r2*wa(i)/raT(i))-1
  dv = r2/raT(i) * EXP(r2*wa(i)/raT(i))
!!!! now do the derivative
  dwa(i) = (v*du - u*dv)/(v*v)
END DO

RETURN
END SUBROUTINE drad_dv

!************************************************************************
! this doe a cumulative K distribution fcn, and then for the first 50%
! of the lowest k values, just does an average of intensity

SUBROUTINE CumulativeK(raKall,raF,raK,raI,iN)


NO TYPE, INTENT(IN)                      :: raKall
REAL, INTENT(IN OUT)                     :: raF(kMaxPts)
REAL, INTENT(IN OUT)                     :: raK(kMaxPts)
REAL, INTENT(IN OUT)                     :: raI(kMaxPts)
INTEGER, INTENT(IN)                      :: iN
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL :: raKAll(kMaxPts)        !this is all 10000 pts k(lowest layer)
REAL :: !these are freqs of pts where DISORT used
REAL :: !these are ks of pts where DISORT used
REAL :: !intensities of pts where DISORT used


! note that we will do a sort on k, and then accordingly, average some of the
! points in raK and raI so as to make things smoother

INTEGER :: iCnt50,iCnt75,iInCnt50,iInCnt75
REAL :: rCnt50,rCnt75,rF,rI,rK,radtot,ttorad
INTEGER :: indx(kMaxPts),iI,iJ,iCnt
REAL :: NormalizedK(kMaxPts),kmin,kmax,dk
REAL :: TempF(kMaxPts),TempI(kMaxPts),TempK(kMaxPts)
REAL :: raCumulativeK(kMaxPts),raCumulativeDistr(kMaxPts)

CALL NumericalRecipesIndexer(indx,raKAll,kMaxPts)

!!set up the normalised vector knorm = ksorted/max(k)
DO iI=1,kMaxPts
  NormalizedK(iI) = raKall(indx(iI))/raKall(indx(kMaxPts))
END DO

kmin=0.0
kmax=NormalizedK(kMaxPts)   !!!should be 1.0

iJ = 500
dk = 1.0/iJ

! compute the cumulative distr function
iCnt=1
DO iI=1,iJ
  raCumulativeK(iI) = kmin + iI*dk
  10     CONTINUE
  IF (NormalizedK(iCnt) <= raCumulativeK(iI)) THEN
    iCnt = iCnt + 1
    IF (iCnt <= kMaxPts) THEN
      GO TO 10
    END IF
  END IF
  iCnt = MIN(iCnt,kMaxPts)
  raCumulativeDistr(iI) = iCnt*1.0/kMaxPts
END DO

! now see which index required before 50% of the cumulative function is used
iCnt50 = 1
iCnt75 = 1
30   CONTINUE
IF (raCumulativeDistr(iCnt50) < 0.5) THEN
  iCnt50 = iCnt50 + 1
  GO TO 30
END IF
40   CONTINUE
IF (raCumulativeDistr(iCnt75) < 0.75) THEN
  iCnt75 = iCnt75 + 1
  GO TO 40
END IF

!find values of raK corresponding to the normalised k=0.5,0.75 values
!remember raK is where DISORT is called, and so these should be
!already sorted from smallest to largest
rCnt50 = raCumulativeK(iCnt50)*raKall(indx(kMaxPts))
rCnt75 = raCumulativeK(iCnt75)*raKall(indx(kMaxPts))
iInCnt50 = 1
iInCnt75 = 1
50   CONTINUE
IF (raK(iInCnt50) < rCnt50) THEN
  iInCnt50 = iInCnt50 + 1
  GO TO 50
END IF
60   CONTINUE
IF (raK(iInCnt75) < rCnt75) THEN
  iInCnt75 = iInCnt75 + 1
  GO TO 60
END IF
!      print *,'normalised cml upto 50% = ',iCnt50,raCumulativeK(iCnt50),rCnt50
!      print *,'normalised cml upto 75% = ',iCnt75,raCumulativeK(iCnt75),rCnt75
!      print *,'input cml upto 50% = ',iInCnt50,raK(iInCnt50)
!      print *,'input cml upto 75% = ',iInCnt75,raK(iInCnt75)

!since there is so much variation in the lowest radainces, and they
!account for so much of raK, let us avg the first few values and
!shove the rest around
IF (raCumulativeK(iCnt50) <= 0.1) THEN
!first initialise the temp arrays
  DO iI=1,iN
    TempF(iI) = raF(iI)
    TempK(iI) = raK(iI)
    TempI(iI) = raI(iI)
  END DO
!50% of the (k/kmax) are smaller in magnitude that 0.1
!this will drive the interpolation wrt raKStep nuts!!!!
!so do an average;
!make sure at least FIVE points remain for the interpolations
  IF ((iN - iInCnt50) <= 5) THEN
    iInCnt50 = iN-5
  END IF
  rF=0.0
  rI=0.0
  rK=0.0
  DO iI=2,iInCnt50
    rF=rF+raF(iI)
    rK=rK+raK(iI)
    rI=rI+radtot(raF(iI),raI(iI))      !!!average the temps!!!
  END DO
  rF=rF/(iInCnt50-1)
  rK=rK/(iInCnt50-1)
  rI=rI/(iInCnt50-1)
  rI=ttorad(rF,rI)                      !!!change to radiance!!!
  
!now shove this info into the arrays
!keep the smallest "k" info!!!!!
  raF(1) = raF(1)
  raK(1) = raK(1)
  raI(1) = raI(1)
!you have averaged the next few to get something
  raF(2) = rF
  raK(2) = rK
  raI(2) = rI
!now fill in the larger "k" values
  DO iI=3,iN - iInCnt50 + 2
    raF(iI) = TempF(iInCnt50 + iI - 2)
    rak(iI) = Tempk(iInCnt50 + iI - 2)
    raI(iI) = TempI(iInCnt50 + iI - 2)
  END DO
  IN = iN - iInCnt50 + 2
END IF

RETURN
END SUBROUTINE CumulativeK
!************************************************************************
