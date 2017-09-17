! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with K.F. Evan's rtspec code *************
!************** Any additional routines are also included here **********
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              rtspec, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer
!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output
    SUBROUTINE doscatter_rtspec(raFreq,raaAbs,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers,iaCldTypes(kmaxclouds)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSUrfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(80) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)

    INTEGER :: i1,i2,iFloor,iDownWard

    DO i1=1,kMaxPts
        raInten(i1)=0.0
    ENDDO
         
! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
    ! radiation travelling upwards to instrument ==> sat looking down
    ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
    ! set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2=iaaRadLayer(iAtm,iNumLayer)-1
        i2=iFloor(i2*1.0/kProfLayer)
        IF (rTSpace > 5.0) THEN
            write(kStdErr,*) 'you want satellite to be downward looking'
            write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
            write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
            write(kStdErr,*) 'Please retry'
            CALL DoSTOP
        END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
    ! radiation travelling downwards to instrument ==> sat looking up
    ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
    ! set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = -1
        i1=iaaRadLayer(iAtm,1)-1
        i1=iFloor(i1*1.0/(1.0*kProfLayer))
        i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), &
        i1,i2
        CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
    END IF

    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
        rAngle=rSatAngle
    ELSE
        rAngle=-rSatAngle
    END IF

    CALL interface_rtspec(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)
     
    RETURN
    end SUBROUTINE doscatter_rtspec

!************************************************************************
! the interface call to RTSPEC
! the main difference here is if the cloud layers for 2 different clouds are
! noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
! artificial cloud of IWP=0.0g/m2 is set for layers 3-40
! kinda based on the interface to DISORT, except that it sets up this
! intermediate "empty" cloud

    SUBROUTINE interface_rtspec( &
! irst the usual kCARTA variables
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
! hen the necessary scattering variables
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

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
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
        
    INTEGER :: iaCldTypes(kMaxClouds)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
    CHARACTER(80) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)

! local variables
! c      IMPLICIT NONE
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
         

! c      INTEGER  MAXTAB, MAXGRID, MAXSCAT
! c      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    INTEGER :: I, JNU1, JNU2
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    REAL ::    SFCTEMP, SFCEMIS
    REAL ::    RADOBS
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::  WAVENO
    CHARACTER(80) :: SCATFILE(MAXSCAT)
    CHARACTER(1) ::   RTMODEL
    CHARACTER(1) :: caScale(MAXSCAT)
!      CHARACTER*24  OUTUNITS, OUTAVERAGING

! new local variables
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: iCldTopkCarta,iCldBotKcarta

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName
    INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iBdry,FindBoundary,iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers,iII,iiDiv
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle,ttorad

! from original v2 rtspec code, to show differences :
!      INTEGER  MAXABSNU, MAXNZ, MAXSPEC
!      PARAMETER (MAXABSNU=100000, MAXNZ=100, MAXSPEC=10)
!      INTEGER NUMSPEC, NSCATTAB, NCLDLAY, NABSNU, NLEV
!      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
!      INTEGER I, J, L, M, JNU1, JNU2
! NTEGER J0, NBAND, IBAND, NWAVENO(MAXSPEC), NOUTNU
! OGICAL TWOSIDEBAND ********* not using outtb(maxspec),outtriang(maxspec)
!      LOGICAL BINARYABSFILE
!      REAL    ZOBS, MUOBS(MAXSPEC), ZTOP, IWP(MAXNZ), DME(MAXNZ)
!      REAL    SFCTEMP, SFCEMIS
! EAL SFCEMIS1(MAXSPEC), SFCEMIS2(MAXSPEC)
! EAL RADOBS(2,MAXSPEC,MAXABSNU), OUTRAD(MAXABSNU) (real radobs,outrad,onu)
!      REAL    HEIGHT(MAXNZ), TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)
!      REAL*8  ABSNU1, ABSNU2, ABSDELNU
!      REAL*8  OUTNU1(MAXSPEC), OUTNU2(MAXSPEC), OUTDELNU(MAXSPEC)
! EAL*8  BANDCENTER(MAXSPEC), BANDOFFSET(MAXSPEC)
!****** used to have REAL*8  WAVENO, ONU1, ONU2, WT1, WT2, INVOUTDELNU
!****** used to have REAL*8 sumwt1, sumwt2, sumrad1, sumrad2
!      REAL*8  OUTNU(MAXABSNU), WAVENO(2,MAXSPEC,MAXABSNU), XO
!      CHARACTER*100 SCATFILE(MAXSCAT), ABSFILE(MAXSPEC), OUTPUTFILE
!      CHARACTER*1   RTMODEL, ABSTYPE
! HARACTER*24  OUTUNITS(MAXSPEC), OUTAVERAGE(MAXSPEC) (was outunits,outavging)
    REAL ::    RAD0UPOBS(MAXNZ), RAD0DNOBS(MAXNZ)

    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
           
    iIOUN = kStdkCarta

!      WRITE (*,*) 'Singlescat, Eddington or hybrid RT model (S, E or H)'
    IF (kScatter == 1) THEN
        RTMODEL='S'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, SingleScat'
    ELSE IF (kScatter == 2) THEN
        RTMODEL='E'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Eddington'
    ELSE IF (kScatter == 3) THEN
        RTMODEL='H'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Hybrid'
    END IF

!      WRITE(*,*) 'Number of scattering tables'
!      NSCATTAB=iNclouds
!      IF (NSCATTAB .GT. MAXSCAT) THEN
!         write(kStdErr,*) 'RTSPEC: MAXSCAT exceeded'
!         CALL DoSTOP
!         END IF
!      WRITE(*,*) 'Input scattering table file name for each table'
!      DO I = 1, NSCATTAB
!        READ (*,'(A100)') SCATFILE(I)
!      ENDDO

!      WRITE (*,*) 'Absorption file type (A-ascii text, B-binary)'
!      READ (*,'(A1)') ABSTYPE
!      BINARYABSFILE = ABSTYPE .EQ. 'B' .OR. ABSTYPE .EQ. 'b'
!      WRITE(*,*) 'Output file name'
!      READ (*,'(A)') OUTPUTFILE
!         Read in the first absorption profile file
!           Absorption file may be binary or ascii
!      CALL READ_ABS_PROFILE (ABSFILE(1), BINARYABSFILE,
!     $               ABSNU1, ABSNU2, ABSDELNU, NABSNU,
!     $               MAXNZ, MAXABSNU, NLEV, HEIGHT, TEMP, ABSPROF)
!      IF (OUTDELNU(1) .LT. ABSDELNU) THEN
!        PRINT *, OUTDELNU(1),ABSDELNU
!        WRITE (kStdErr,*) 'RTSPEC: output wavenumber resolution less than',
!     $                ' in absorption file'
!        CAll DoSTOP
!      ENDIF

    CALL SetMieTables_RTSPEC(raFreq, & &
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
    -1,              &  & !!!!iSergio = -1 as this is ORIG RTSPEC code
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT,  iCldTopKcarta,iCldBotkCarta)

    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
    iProfileLayers,raPressLevels)

    IF (iDownWard > 0) THEN
    ! own look instr, in kCARTA layering style
        IOBS   = iNumLayer
    ELSE IF (iDownWard < 0) THEN    !up look instr
    ! p look instr, in KCARTA layering style
        IOBS   = 1
    END IF
! hange to RTSPEC layering
    iobs=(iNumLayer+1)-iobs+1

!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
    IF (iCloudySky < 0) THEN
        rSolarAngle = kSolarAngle
    !!!!note that we do not care about background thermal accurately here
        write (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
        DO iii = 1,iNp
            DO iF = 1,kMaxPts
                DO iL = 1,NLEV-1
                    raTau(iL)  = absprof(iL,iF)
                END DO
                CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
                rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
                raUseEmissivity(iF),raFreq(iF),raInten(iF),iaOp(iii), &
                iProfileLayers,raPressLevels)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END DO
        GOTO 9876
    END IF

!     Get spectral range in absorption file for this spectral segment
    JNU1=1
    JNU2 = kMaxPts

    MUOBS=cos(rSatAngle*kPi/180.0)  !viewing angle cosine
    IF (iDownWard == -1) THEN
        MUOBS=-MUOBS
    END IF
    SFCTEMP = rSurfaceTemp
    WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',muobs,sfctemp

! see if there are any layers between TOA and instr; if there are, set
! TOA_to_instr to the cumulative radiative transfer of the necessary k, else
! set it to -10.0
    DO iL=1,kProfLayer
        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
    END DO
    iiDiv = 0
    555 CONTINUE
    IF (iiDiv*kProfLayer < iaRadLayer(1)) THEN
        iiDiv = iiDiv + 1
        GOTO 555
    END IF
    iiDiv = iiDiv - 1

    iBdry=-1
! ind the layer where you quit doing background thermal using only
! cos(3/5) and start using more accurate approx
    IF (kThermal >= 0) THEN
        iBdry=FindBoundary(raFreq,iProfileLayers,raPresslevels,iaRadLayer)
    END IF

    CALL Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,TOA_to_instr)

    DO iii = 1,iNp                     !!!!!!loop over outputs iii=1,iNp
        IF (iDownWard == +1)  THEN
            iobs = iaOp(iii)+1
            iobs = (iNumLayer+1) - iobs + 1
        ELSE
            iobs = iaOp(iii)
        END IF
        iObs = iObs + (iiDiv*kProfLayer)

        IF (ICLDTOP == 0 .OR. ICLDBOT == 0 .OR. IOBS == 0) THEN
            WRITE (kStdErr,*) 'RTSPEC: Observer or cloud top height does', &
            ' not match absorption file levels'
            STOP
        END IF

    !        print *,rtmodel,muobs,iobs,nlev,sfctemp,ncldlay,icldtop,icldbot
    !        DO iF = JNU1, JNU2            !!!!loop over frequencies
    !          raInten(iF) = 0.0
    !          END DO
    !        jnu1 = 3080
    !        jnu2 = 3080

        DO iF = JNU1, JNU2            !!!!loop over frequencies
            IF ((iObs == 1) .AND. (iDownWard < 0)) THEN
                raInten(iF) = ttorad(raFreq(iF),sngl(kTSpace))
            ELSE
                DO iL = 1,NLEV-1
                    taugas(iL) = absprof(iL,iF)
                END DO

                SFCEMIS=raUseEmissivity(iF)
                WAVENO = raFreq(iF)
                CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
                MUOBS, IOBS, WAVENO, &
                NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
                NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, & &
            !!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID,
                NSCATTAB, MUINC, &
                NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
                TABEXTINCT, TABSSALB, TABASYM, &
                TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
                RADOBS, &
                IBDRY, TOA_to_instr(iF),iiDiv) !these are 3 new parameters
                raInten(iF)=radobs

            END IF                  !!!!end actually computing rad if things ok
        ENDDO                     !!!!loop over frequencies
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END DO

    9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm

    RETURN
    end SUBROUTINE interface_rtspec

!************************************************************************
! this function quickly computes the background thermal contribution downto
! the layer where the instrument is
    REAL FUNCTION backgnd(raaAbs,waveno,raVTemp,iaaRadLayer,iNumLayer, &
    iNpmix,iAtm,rFracTop,j)

    include '../INCLUDE/scatterparam.f90'

    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,raVTemp(kMixFilRows),waveno
    INTEGER :: iAtm,j,iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iNpmix
!     raFreq(j)=waveno

    INTEGER :: iLay,iaRadLayer(kProfLayer),iaRadLayerTemp(kMixFilRows)
    INTEGER :: iT,iExtraThermal
    REAL :: rad,ttorad,rMPTemp,raExtraThermal(kMaxPts),rThetaEff
    REAL :: rTrans, rEmiss, rAbs

    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

    CALL AddUppermostLayersQ(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

    rThetaEff=3/5.0

! see  SUBROUTINE Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace,

! do radiation at space temp
    rMPTemp = sngl(kTSpace)
    rad = ttorad(waveno,rMPTemp)

! go from top of atmosphere to instrument
    IF (iExtraThermal < 0) THEN
    ! just have to worry about top layer being fractional
        rAbs=raaAbs(j,iaRadLayer(iNumLayer))*(1-rFracTop)
        rTrans=exp(-rAbs/rThetaEff)

        rMPTemp=raVTemp(iaRadLayer(iNumLayer))
        rEmiss=ttorad(waveno,rMPTemp)
        rEmiss=(1.0-rTrans)*rEmiss

        rad=rad*rTrans + rEmiss
    ELSE
    ! go from top of atmosphere to layer above instrument
        DO iLay=iT,iNumLayer+1
            rAbs=raaAbs(j,iaRadLayerTemp(iLay))
            rTrans=exp(-rAbs/rThetaEff)
            rMPTemp=raVTemp(iaRadLayerTemp(iNumLayer))
            rEmiss=ttorad(waveno,rMPTemp)
            rEmiss=(1.0-rTrans)*rEmiss
            rad=rad*rTrans + rEmiss
        END DO
    ! do layer where instrument is .. don't worry about interpolating temp
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
    end FUNCTION backgnd
!************************************************************************
! include Evans code here
! cccc      include '../INCLUDE/scatter_rtspec.f'   this is compiled separately

!      PROGRAM RTSPEC   NEW!!!!!!!!!!!!!!! for up and down look!!!!!!
!       Program for computing thermal atmospheric radiative transfer with
!     a multiple layer cloud in an atmosphere.  Outputs a radiance spectrum,
!     suitably averaged.  Three approximate solutions for the scattering
!     radiative transfer in the cloud layers are available: Single scattering
!     approximation, Eddington second approximation and a hybrid single
!     scattering/Eddington approximation.  Comparison with an "exact"
!     doubling-adding model for ice clouds has shown that the hybrid model
!     is superior in the mid IR (300-3000 cm^-1) but not always in the
!     submm (0-50 cm^-1) (for particle sizes Dme=30 to 300 um).
!     Reference:  Deeter, M. N. and K. F. Evans, 1998:
!         A hybrid Eddington-single scattering radiative transfer model
!         for computing radiances from thermally emitting atmospheres.
!         J. Quant. Spectosc. Radiat. Transfer, 60, 635-648.
!       The scattering approximations do radiative transfer for a single
!     homogeneous layer.  The input to them is the incident radiation on
!     the layer.  The incident radiance is computed by a nonscattering
!     radiative transfer integration of cloud and gaseous absorption profiles
!     that are read in from files.  The scattering routines output
!     upwelling radiance at the top of the layer, and this is propagated
!     to the observation level for output.  The multiple layers are dealt
!     with by computing the appropriate incident radiances/fluxes for each
!     layer, and then running the single layer routines for each layer from
!     the bottom up. For upward looking geometry, the downwelling radiance
!     is computed by just switching the incident radiation and the top/bottom
!     temperature of the layer, thus the "upwelling" radiance that subroutine
!     outputs is actually the downwelling. Note: It is not necesary to switch
!     the "phi" functions since "+/-" doesn't mean up/down just forward/back.
!     For multiple layers, the single layer routines must be run from top down.
!       The cloud boundaries must be at levels in the input absorption profile.
!     The surface is assumed to reflect specularly with the specified
!     emissivity. Multiple reflection between the cloud (scattering layer)
!     and surface is ignored.  The input includes the cloud properties,
!     which are ice water path (IWP) and particle size (Dme) for each layer.
!     There may be multiple input scattering table files for different
!     cloud types.
!       The output is upwelling or downwelling brightness temperature or
!     radiance as a function of wavenumber at the desired observation level.
!     The radiance is calculated at the absorption file wavenumbers.
!     For output the radiance may be averaged with either a triangular
!     or rectangular bandpass or convolved with sinc or sinc squared.
!     The width of the rectangular bandpass or the full width half max
!     (FHHM) of the triangular bandpass is equal to the output wavenumber
!     spacing.  The output wavenumbers, which are the center of the
!     bandpasses, are on the regular grid specified by the user.  The
!     convolution averaging (for a Fourier Transform Spectrometer) is
!     done by FFTing and filtering with either a rectangular (sinc) or
!     triangular (sinc squared) windowing function.
!       Multiple spectral segments may be output, with potentially different
!     absorption profile input files, output units, and output bandpass
!     averaging for each one.  Multiple "channels" may be output by making
!     the output spacing for each segment equal to zero, setting the
!     channel width, and using rectangular averaging. The channel width
!     extends equally from each side of the output wavenumber. If the channel
!     width is one absorption file wavenumber spacing or less then
!     monochromatic radiative transfer is done.  A double sideband receiver
!     may be emulated, in which case two rectangular bandpasses equally
!     spaced from the receiver center frequency are simulated.

!       The input ascii absorption file has the following format:
!     1) header line
!     2) starting, ending, delta, and number of wavenumbers in file
!     3) number of levels (1 + number of atmosphere layers)
!     4) heights of levels (from top to bottom) (in km)
!     5) temperatures at levels (in Kelvin)
!     wavenumber (cm^-1) and layer optical depths for each Nlev-1 layers
!        . . .
!     A binary format absorption file may be input for faster loading
!     (the BINARYABSFILE logical is set from user input; see
!     READ_ABS_PROFILE subroutine).

!       The input scattering tables are in the format produced by
!     sscatmie.f.  This format also accomodates single scattering
!     information for randomly nonspherical oriented particles.

!     Author: Frank Evans  University of Colorado   May 1999
!     Co-authors:
!       Aaron Evans    (multilayer code, spectral segments, double sidebands,
!                       upward looking geometry, output convolution)
!       Merritt Deeter (single scattering layer routines)

!************************************************************************
