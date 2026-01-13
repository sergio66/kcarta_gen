c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with K.F. Evan's rtspec code *************
c************** Any additional routines are also included here **********
c************************************************************************
c************** All the routines in this file are necessary *************
c************************************************************************

c************************************************************************
c note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
c              rtspec, layer 1 == TOA, layer kProfLayer = ground
c                      there are nlev = 1 + iNumlayer  levels
c************************************************************************


c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output
      SUBROUTINE doscatter(raWaves,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm)

      include 'kcarta.param'
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
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
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID
      CHARACTER*80 caOutName
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
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

      CALL interface_rtspec(raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward) 
 
      RETURN
      END

c************************************************************************
      SUBROUTINE interface_rtspec(        
        !first the usual kCARTA variables
     $        raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
         !then the necessary scattering variables
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable, 
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward) 

      include 'scatter.param'
      include 'NewRefProfiles/outpresslevels.param'

c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
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

      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iDownWard,iBinaryFile
      CHARACTER*80 caOutName
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 

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
      INTEGER I, L, JNU1, JNU2
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
      REAL    SFCTEMP, SFCEMIS
      REAL    RADOBS
      REAL    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
      REAL  ABSNU1, ABSNU2, ABSDELNU
      REAL  WAVENO
      CHARACTER*80 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
c      CHARACTER*24  OUTUNITS, OUTAVERAGING

c new local variables
      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iBdry,FindBoundary,iaRadLayer(kProfLayer)

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

      iIOUN=kStdkCarta

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
c copied from s_scatter_spectra.f .. all table names etc are unique, so no 
c need to make more checks
      NSCATTAB=-1000
      DO iIn=1,kMaxClouds*kCloudLayers 
        iaTable(iIn)=-1
        END DO 
      DO iIn=1,MAXSCAT
        ScatFile(iIn)='                                                    ' 
        END DO 
      DO iIn=1,iNclouds 
        DO iJ=1,iaCloudNumLayers(iIn) 
          iI=iaaScatTable(iIn,iJ) 
          IF (iI .GT. MAXSCAT) THEN
            write(kStdErr,*)'unfortunately, in interface_rtspec, '
            write(kStdErr,*)'MAXSCAT < kMaxClouds*kCloudLayers'
            write(kStdErr,*)'please reset and retry'             
            CALL DoSTOP
            END IF
          caName=caaaScatTable(iIn,iJ) 
          IF (iaTable(iI) .LT. 0) THEN  !nothing associated with this yet 
            IF (iI .GT. NSCATTAB) THEN
              NSCATTAB=iI
              END IF
            iaTable(iI)=1 
            ScatFile(iI)=caName
            END IF
          END DO 
        END DO  

C         Read in scattering tables
      IF (iBinaryFile .EQ. 1) THEN
        DO I = 1, NSCATTAB 
          CALL READ_SSCATTAB_BINARY(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,
     $          NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), 
     $          NWAVETAB(I), WAVETAB(1,I),
     $          MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
     $          TABPHI1UP(1,I), TABPHI1DN(1,I), 
     $          TABPHI2UP(1,I), TABPHI2DN(1,I))
          IF ((ABS(MUINC(1)-0.2113) .GT. 0.001) .OR.
     $      (ABS(MUINC(2)-0.7887) .GT. 0.001)) THEN
            write (kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
            CALL DoStop
            END IF
          ENDDO
      ELSE IF (iBinaryFile .EQ. -1) THEN
        DO I = 1, NSCATTAB 
          CALL READ_SSCATTAB(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,
     $          NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), 
     $          NWAVETAB(I), WAVETAB(1,I),
     $          MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
     $          TABPHI1UP(1,I), TABPHI1DN(1,I), 
     $          TABPHI2UP(1,I), TABPHI2DN(1,I))
          IF ((ABS(MUINC(1)-0.2113) .GT. 0.001) .OR.
     $      (ABS(MUINC(2)-0.7887) .GT. 0.001)) THEN
            write (kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
            CALL DoStop
            END IF
          ENDDO
        END IF

c      WRITE(*,*) 'Number of cloud layers'
c      READ (*,*) NCLDLAY
c      IF (NCLDLAY .GT. MAXNZ) STOP 'RTSPEC: MAXNZ exceeded by NCLDLAY'
c      WRITE (*,*) 'IWP/LWP (g/m^2), Dme (micron), scat table # for each layer'
c      DO L = 1, NCLDLAY
c        READ (*,*) IWP(L), DME(L), ISCATTAB(L)
c      ENDDO
c code from s_scatter_spectra.f
      iCloud=-1
      DO iJ=1,iNclouds 
        DO iScat=1,iaCloudNumAtm(iJ) 
          IF (iaaCloudWhichAtm(iJ,iScat) .EQ. iAtm) THEN 
            iCloud=iJ             !this one cloud is associated with this atm
            END IF 
          END DO 
        END DO
      IF (iCloud .LT. 0) THEN
        write(kStdWarn,*)'Could not find a cloud for atmosphere #',iAtm
        write(kStdWarn,*)'setting IWP = -100.0'
        iCloud=1    !say cloud number one is associated with this atmosphere
        ncldlay=1   !say fictitious cloud occupies one layer
        IWP(1)      = -100.0   !but make sure cloud has NO particles in it!
        DME(1)      = -10.0    !but make sure cloud has NO particles in it!
        ISCATTAB(1) = -1
      ELSE
        NCLDLAY=iaCloudNumLayers(iCloud)
        write(KStdWarn,*) 'Cloud number, num layers = ',iCloud,NCLDLAY
        DO L = 1, NCLDLAY
          IWP(L) = raaaCloudParams(iCloud,L,1) 
          DME(L) = raaaCloudParams(iCloud,L,2) 
          ISCATTAB(L) = iaaScatTable(iCloud,L)
          ENDDO
        END IF

c not needed
c      WRITE(*,*) 'Observation level (km)'
c      READ (*,*) ZOBS
c      WRITE (*,*) 'Input cloud top height (km)'
c      READ (*,*) ZTOP

C     Find the levels for top of cloud and observation level
c     remember that these numbers are with respect to the KLAYERS pressure
c     levels and layers
c     these will be reset when they are passed in and out of GetAbsProfile
c     NOTE : here we are still in kCARTA frame ie 
c
c   TOA    --------------
c          layer iNumlayer
c          --------------
c              .....
c          --------------             KCARTA
c             layer 2
c          --------------
c             layer 1
c   GND --------------------------
c
c when we call GetAbsProfile, the variables icldtop,icldbot,iobs will be reset
c to reflect the rtspec layering
c   TOA    --------------
c             layer 1
c          --------------
c              .....                 RTSPEC
c          --------------
c        layer iNumLayer-1
c          --------------
c         layer iNumLayer
c   GND --------------------------

      IF (IWP(1) .GT. 0.0) THEN
        ICLDTOP=iaaCloudWhichLayers(iCloud,1)+1
        ICLDBOT=iaaCloudWhichLayers(iCloud,iaCloudNumLayers(iCloud))
        IF (iDownWard .GT. 0) THEN
          IOBS   =iNumLayer       
        ELSE IF (iDownWard .LT. 0) THEN
          IOBS   =1
          END IF
        END IF
       
      IF (IWP(1) .LE. 0.0) THEN  !we have no cloud; set up fictitious clouds
        IF (iDownWard .GT. 0) THEN    
          !down look instr : set cloud BELOW observer, in layer #1
          ICLDTOP=2
          ICLDBOT=1
          IOBS   =iNumLayer     
        ELSE IF (iDownWard .LT. 0) THEN    !up look instr
          !up look instr : set cloud ABOVE observer, in layer #iNumLayer
          ICLDTOP=iNumLayer+1
          ICLDBOT=iNumLayer
          IOBS   =1             
          END IF
        END IF

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
      CALL GetAbsProfile(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1))

C     Get spectral range in absorption file for this spectral segment
      JNU1=1
      JNU2=KMaxPts

      MUOBS=cos(rSatAngle*kPi/180.0)  !viewing angle cosine
      IF (iDownWard .EQ. -1) THEN
        MUOBS=-MUOBS
        END IF
      SFCTEMP = rSurfaceTemp
      WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',muobs,sfctemp

      iBdry=-1
      !find the layer where you quit doing background thermal using only 
      !acos(3/5) and start using more accurate approx
      IF (kThermal .GE. 0) THEN
        iBdry=FindBoundary(raWaves) 
        END IF

c see if there are any layers between TOA and instr; if there are, set
c TOA_to_instr to the cumulative radiative transfer of the necessary k, else 
c set it to -10.0
      DO iL=1,kProfLayer
        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
        END DO
      CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,TOA_to_instr) 

C     Calculate the upwelling observed radiance for this wavenumber
c note that i have rewritten gasrt1 so that it automatically computes the
c backgnd thermal contribution ie do not need to initialise
cc        IF (kThermal .LT. 0) THEN
cc          radobs=0.0
cc        ELSE 
cc          radobs=backgnd(raaAbs,waveno,raVTemp,iaaRadLayer,iNumLayer,
cc     $                      iNpmix,iAtm,rFracTop,iF)
cc          END IF

      DO iF = JNU1, JNU2
        DO iL=1,NLEV-1
          taugas(iL)=absprof(iL,iF)
          END DO
        SFCEMIS=raUseEmissivity(iF)
        WAVENO = raWaves(iF)
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
     $               IBDRY, TOA_to_instr(iF))    !these are 2 new parameters
        raInten(iF)=radobs
        ENDDO

      CALL wrtout(iIOUN,caOutName,raWaves,raInten)

      RETURN
      END

c************************************************************************
c this subroutine takes in the input abs coeffs (raaAbsCoeff) where raa(1,:)
c is the lowest layer and raa(kProfLayer,:) is the highest layer .. it then 
c outputs these  abs coeffs intp absprof, where absprof(1,:) is the top, and
c absprof(iNumLayer,:) = ground

C     sets optical depths for NLEV-1 layers and NABSNU wavenumbers. 
C     The temperature (K) of the profile is also returned. (no height needed)

      SUBROUTINE GetAbsProfile(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS, iDownward, iwp)

      include 'scatter.param' 

c these are variables that come in from kcartamain.f 
      REAL raaAbs(kMaxPts,kMixFilRows),raWaves(kMaxPts),rFracTop,rFracBot
      REAL raVTemp(kMixFilRows),rSurfaceTemp
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNpmix
      INTEGER iDownWard
      REAL iwp
c these are variables that we have to set
      INTEGER  NABSNU, NLEV         !!!!!!!!!!!!!MAXNZ, MAXABSNU 
      REAL   ABSNU1, ABSNU2, ABSDELNU 
      REAL   TEMP(*), ABSPROF(MAXNZ,*) 
      INTEGER  ICLDTOP,iCLDBOT,IOBS

c local variables
      INTEGER iaRadLayer(kProfLayer), iFr, iL, iLay 
      REAL NU

c these are to flip the temperature, abs profiles if instr looks up
      REAL raTemp(kProfLayer+1),raaTemp(kProfLayer,kMaxPts)

      absnu1=raWaves(1)
      absnu2=raWaves(kMaxPts)
      absdelnu=(absnu2-absnu1)/(kMaxPts-1)
      nabsnu=kMaxPts
      nlev=iNumLayer+1           !this is the number of pressure levels

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

c set the vertical temperatures of the atmosphere 
      CALL SetRTSpecTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,rSurfaceTemp)

c now set up the abs coeffs
c initialize array to all zeroes
c      print *,'rFracTop,rFracBot = ',rFracTop,rFracBot
      DO iFr=1,kMaxPts
        DO iLay=iNumLayer+1,kProfLayer
          absprof(iLay,iFr)=0.0
          END DO
        END DO

       DO iLay=1,iNumLayer
         iL=iaRadLayer(iLay) 
         nu=1.0
         IF (iLay .EQ. 1) THEN
           nu=rFracBot           
         ELSE IF (iLay .EQ. iNumLayer) THEN
           nu=rFracTop
           END IF
         DO iFr=1,kMaxPts 
           !absprof wants level 1 == TOA, level iNumLayer= gnd
           absprof(iNumLayer-iLay+1,iFr)=raaAbs(iFr,iL)*nu
           END DO 
         END DO 

c now set icldtop icldbot, iobs
c iDownward = +1 ==> downward looking instrument
c             -1 ==> upward looking instrument
c remember there is ONE more level than there are layers
      icldtop=(iNumLayer+1)-icldtop+1
      icldbot=(iNumLayer+1)-icldbot+1
      iobs=(iNumLayer+1)-iobs+1

cx1      if (iDownWard .gt. 0) then
c        iobs=1
c        icldtop=iNumLayer-icldtop+1
c        icldbot=iNumLayer-icldbot+1
c      else if (iDownWard.lt. 0) then
c        iobs=iNumLayer
c        icldtop=iNumLayer-icldtop+1
c        icldbot=iNumLayer-icldbot+1

      IF  (iDownWard.lt. 0) THEN
        !now flip TEMPerature array
        !remember there is one more level than there are layers
        DO  iLay=1,iNumLayer+1
          raTemp(iLay)=TEMP((iNumLayer+1)-iLay+1)
          END DO
        DO  iLay=1,iNumLayer+1
          TEMP(iLay)=raTemp(iLay)
          END DO
 
        !now flip absprof array
        DO iFr=1,kMaxPts 
          DO  iLay=1,iNumLayer
            raaTemp(iLay,iFr)=absprof(iNumLayer-iLay+1,iFr)
            END DO
          END DO
        DO iFr=1,kMaxPts 
          DO  iLay=1,iNumLayer
            absprof(iLay,iFr)=raaTemp(iLay,iFr)
            END DO
          END DO
        END IF

      RETURN 
      END
c************************************************************************
c set the vertical temperatures of the atmosphere 
c this sets the temperatures at the pressure level boundaries, using the
c temperatures of the pressure layers that have been supplied by kLayers
      SUBROUTINE SetRTSpecTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,rSurfaceTemp)

      include 'scatter.param'
      include 'NewRefProfiles/outpresslevels.param'

c these are variables that come in from kcartamain.f 
      REAL raVTemp(kMixFilRows),rSurfaceTemp
      INTEGER iNumLayer
c these are variables that we have to set
      REAL    TEMP(*)

c local variables
      INTEGER iaRadLayer(kProfLayer), iL, iLay ,iM, idiv
      REAL FindBottomTemp,Temp1(maxnz)
      REAL pavg(kProfLayer),rP,raProfileTemp(kProfLayer)

      DO iLay=1,MAXNZ
        Temp1(iLay)=0.0
        END DO

      DO iLay=1,kProfLayer
        pavg(iLay)=plev(iLay+1)-plev(iLay)
        pavg(iLay)=pavg(iLay)/log(plev(iLay+1)/plev(iLay))
        END DO

c see which set of Mixed Paths the current atmosphere occupies eg 
c set 1 = 1..100, set2= 101..200 etc
c eg if current atmosphere is from MixfilPath 110 to 190, and kProfLayer = 100,
c then we set iMod as 2      idiv(150,100)=1  === 2nd set of mixed paths
c assume each atmosphere has at least 30 layers in it!!!
      iM=idiv(iaRadLayer(30),kProfLayer)+1
      DO iLay=1,kProfLayer
        raProfileTemp(iLay)=raVTemp(iLay+(iM-1)*kProfLayer)
        END DO

      DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
        !map this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
        iL=iL-idiv(iL,kProfLayer)*kProfLayer  
        IF (iL .EQ. 0) THEN
c          iL=idiv(iaRadLayer(iLay),kProfLayer)*kProfLayer
          iL=kProfLayer
          END IF
        rP=plev(iL+1)-10000*delta
        if (rp .LT. plev(kProfLayer+1)) then
          rp = plev(kProfLayer+1)+10000*delta
          end if
        TEMP1(iNumLayer-iLay+1)=FindBottomTemp(rP,raProfileTemp)
c        print *,iLay,iL,rP,TEMP1(iNumLayer-iLay+1)
        END DO
      TEMP1(iNumLayer+1)=rSurfaceTemp

      DO iLay=1,iNumLayer+1
        temp(iLay)=temp1(iLay)
        END DO

c      print *,0,raProfileTemp(iLay),temp1(iNumLayer+1)
c      DO iLay=1,iNumlayer
c        avg=0.5*(temp1(iNumLayer-iLay+1)+temp1(iNumLayer-iLay+2))
c        print *,ilay,raProfileTemp(iLay),temp1(iNumLayer-iLay+1),avg
c        end do

      RETURN
      END
c************************************************************************
c this function quickly computes the background thermal contribution downto
c the layer where the instrument is
      REAL FUNCTION backgnd(raaAbs,waveno,raVTemp,iaaRadLayer,iNumLayer,
     $                      iNpmix,iAtm,rFracTop,j)

      include 'scatter.param'

      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,raVTemp(kMixFilRows),waveno
      INTEGER iAtm,j,iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iNpmix
c     raWaves(j)=waveno

      INTEGER iLay,iaRadLayer(kProfLayer),iaRadLayerTemp(kProfLayer)
      INTEGER iT,iExtraThermal
      REAL rad,r1,r2,rMPTemp,raExtraThermal(kMaxPts),rThetaEff
      REAL rTrans, rEmiss, rAbs

      r1=kPlanck1 
      r2=kPlanck2       

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
      rMPTemp=kTSpace
      rad=exp(r2*waveno/rMPTemp)-1.0 
      rad=r1*((waveno**3))/rad

c go from top of atmosphere to instrument 
      IF (iExtraThermal .LT. 0) THEN
c just have to worry about top layer being fractional 
        rAbs=raaAbs(j,iaRadLayer(iNumLayer))*(1-rFracTop)
        rTrans=exp(-rAbs/rThetaEff)

        rMPTemp=raVTemp(iaRadLayer(iNumLayer))
        rEmiss=exp(r2*waveno/rMPTemp)-1.0  
        rEmiss=r1*((waveno**3))/rEmiss
        rEmiss=(1.0-rTrans)*rEmiss

        rad=rad*rTrans + rEmiss
      ELSE    
c go from top of atmosphere to layer above instrument 
        DO iLay=iT,iNumLayer+1
          rAbs=raaAbs(j,iaRadLayerTemp(iLay))
          rTrans=exp(-rAbs/rThetaEff)
          rMPTemp=raVTemp(iaRadLayerTemp(iNumLayer))
          rEmiss=exp(r2*waveno/rMPTemp)-1.0  
          rEmiss=r1*((waveno**3))/rEmiss
          rEmiss=(1.0-rTrans)*rEmiss
          rad=rad*rTrans + rEmiss
          END DO
c do layer where instrument is .. don't worry about interpolating temp
        DO iLay=iNumLayer,iNumLayer
          rAbs=raaAbs(j,iaRadLayerTemp(iLay))*(1-rFracTop)
          rTrans=exp(-rAbs/rThetaEff)
          rMPTemp=raVTemp(iaRadLayerTemp(iNumLayer))
          rEmiss=exp(r2*waveno/rMPTemp)-1.0  
          rEmiss=r1*((waveno**3))/rEmiss
          rEmiss=(1.0-rTrans)*rEmiss
          rad=rad*rTrans + rEmiss
          END DO
        END IF
 
      backgnd=rad*2.0*kPi

      RETURN
      END
c************************************************************************
c include Evans code here
cccccc      include 'scatter_rtspec.f'   this is compiled separately

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

