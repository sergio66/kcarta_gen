c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with STamnes et al  disort code  *********
c************** Any additional routines are also included here **********
c************************************************************************
c************** All the routines in this file are necessary *************
c************************************************************************

c************************************************************************
c note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
c              disort, layer 1 == TOA, layer kProfLayer = ground
c                      there are nlev = 1 + iNumlayer  levels
c************************************************************************


c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output
      SUBROUTINE doscatter_disort(raWaves,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iDoFlux     = do radiance or flux computation
c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
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
c raNumberDensity = P/RT == number of particles in each layer of atm
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raNumberDensity(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile,iDoFlux
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
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

      CALL interface_disort(raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
c     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        iNp,iaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
c     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
c     $        raLayAngles,raSunAngles,
     $         raLayAngles,iDoFlux,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity) 
 
      RETURN
      END

c************************************************************************
c this subroutine sets up the scattering table info from SSCATMIE.F
      SUBROUTINE SetMieTables_DISORT(           
     $        !!!!!!!!!!!!!!!!!these are the input variables
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable, 
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        !!!!!!!!!!!!!!!!!!these are the output variables
     $    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,
     $    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,
     $    TABPHI2UP, TABPHI2DN,
     $    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, 
     $    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,
     $    iCloudySky, IACLDTOP, IACLDBOT)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c ---------------- inputs needed to read scattering tables -------------------
c this is which atm number is being used, and whether these are binary files
      INTEGER iAtm,iBinaryFile,iNumLayer,iDownward
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
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
c this is just to set everything about clouds relative to TOA layer
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)

c ---------------- outputs from the scattering tables -------------------
c --------------------- produced by Evans Mie code ----------------------
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
      CHARACTER*80 SCATFILE(MAXSCAT)

      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

      INTEGER NSCATTAB, NCLDLAY, NLEV, NABSNU
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      REAL    IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)

      INTEGER iCloudySky        !!!!are there clouds in this atm???
 
c ---------------------------- local variables ----------------------------
      INTEGER iaTable(kMaxClouds*kCloudLayers),iIn,iJ,iReadTable,I
      INTEGER iCloud,iStep
      REAL extinct
      INTEGER LL,II,N,M,iB_atm,iT_Atm,iLayers

      CHARACTER*80 caName
      CHARACTER*1 caScale(MAXSCAT)

      !initialise all scattering info to null

      !!!! need to reference the cloud tops and bottoms wrt TOP layer of
      !!!! defined atm 
      !!!! eg if atm from 971 to 150 mb (plane at 150 mb) ==> 
      !!!!       this occupies kCARTA layers 19-78
      !!!!   if cloud from 248 to 214 mb  ==> 
      !!!!       this occupies kCARTA layers 69 to 72
      !!!! Thus the cloud occupies RTSPEC atmosphere "tau" from 6 to 9
      IF (iDownWard .EQ. 1) THEN
        iB_Atm=iaaRadLayer(iAtm,1) 
        iT_Atm=iaaRadLayer(iAtm,iNumLayer) 
      ELSEIF (iDownWard .EQ. -1) THEN
        iT_Atm=iaaRadLayer(iAtm,1) 
        iB_Atm=iaaRadLayer(iAtm,iNumLayer) 
        END IF

      !!!!!!hwowever we also do fluxes, so even if the atm is defined so it
      !!!!!!is for an uplook instrument, RTSPEC will be called in a downlook
      !!!!!!fashion, and vice versa

      IF (iB_Atm .GT. iT_Atm) THEN
        iCloudySky = iT_Atm
        iT_Atm = iB_Atm
        iB_atm = iCloudySky
        END IF

      iCloudySky = -1  !!!!!!! assume NO clouds associated with this atm

      NSCATTAB=-1000   !!!!!!! total of how many scattering tables (files) 
      DO iIn=1,kMaxClouds*kCloudLayers 
        iaTable(iIn)=-1
        END DO 
      DO iIn=1,MAXSCAT
        ScatFile(iIn)=' ' 
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
          IF (iI .GT. MAXSCAT) THEN
            write(kStdErr,*)'in interface_disort, you have set it up so'
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

        !!check to see if this cloud is to be used with this atm
        DO iJ=1,iaCloudNumAtm(iIn) 
          IF (iaaCloudWhichAtm(iIn,iJ)  .EQ. iAtm) THEN
            iCloudySky = iIn         !!!! set this up
            iaCloudWithThisAtm(iIn) = 1

            !!!!!these are the kCARTA layers 1 to 100 = GND to TOA
            IACLDTOP(iIn)=iaaCloudWhichLayers(iIn,1)+1
            IACLDBOT(iIn)=iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn))

            write(kStdWarn,*)'cloud # ',iIn,' associated with atm # ',iAtm
            write(kStdWarn,*)'cloud is in KCARTA layers ',
     $                        iaCldTop(iIn)-1,' to ',iaCldBot(iIn)

            !!!!!these are the DISORT layers 100 to 1 = GND to TOA
            iaCldbot(iIn) = iT_Atm - iaCldbot(iIn) + 1
            iaCldtop(iIn) = iT_Atm - iaCldtop(iIn) + 1

            iaCldBot(iIn) = iaCldBot(iIn) + 1
            iaCldTop(iIn) = iaCldTop(iIn) + 1

            write(kStdWarn,*)'cloud is in DISORT layers ',
     $                        iaCldTop(iIn)+1,' to ',iaCldBot(iIn)

            END IF
          END DO

       !!check to see which scattering tables to be used with this atm
        DO iJ=1,iaCloudNumLayers(iIn) 
          iI=iaaScatTable(iIn,iJ) 
          IF (iaCloudWithThisAtm(iIn)  .EQ. 1) THEN
            iaScatTable_With_Atm(iI) = 1
            write(kStdWarn,*)'scatter table ',iI,' used with atm # ',iAtm
            END IF
          END DO

        END DO      !!!!!!!!main       DO iIn=1,iNclouds 

C     Only read in scattering tables that are needed for this atm
      iReadTable = 1
      IF (iReadTable .GT. 0) THEN
        IF (iBinaryFile .EQ. 1) THEN
          DO I = 1, NSCATTAB 
            IF (iaScatTable_With_Atm(I) .GT. 0) THEN
              write(kStdWarn,*) 'Reading binary scatter data for table #',I
              CALL READ_SSCATTAB_BINARY(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,
     $          caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), 
     $          NWAVETAB(I), WAVETAB(1,I),
     $          MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
     $          TABPHI1UP(1,I), TABPHI1DN(1,I), 
     $          TABPHI2UP(1,I), TABPHI2DN(1,I))
              IF ((ABS(MUINC(1)-0.2113) .GT. 0.001) .OR.
     $        (ABS(MUINC(2)-0.7887) .GT. 0.001)) THEN
                write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                CALL DoStop
                END IF
              END IF
            ENDDO
        ELSE IF (iBinaryFile .EQ. -1) THEN
          DO I = 1, NSCATTAB 
            IF (iaScatTable_With_Atm(I) .GT. 0) THEN
              write(kStdWarn,*) 'Reading ascii scatter data for table #',I
              CALL READ_SSCATTAB(SCATFILE(I),  !!!!!!MAXTAB, MAXGRID,
     $          caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), 
     $          NWAVETAB(I), WAVETAB(1,I),
     $          MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
     $          TABPHI1UP(1,I), TABPHI1DN(1,I), 
     $          TABPHI2UP(1,I), TABPHI2DN(1,I))
              IF ((ABS(MUINC(1)-0.2113) .GT. 0.001) .OR.
     $          (ABS(MUINC(2)-0.7887) .GT. 0.001)) THEN
                write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                CALL DoStop
                END IF
              END IF
            ENDDO
          END IF    !iBinaryFile .GT. 0
        END IF      !iReadTable  .GT. 0

      !check to see that all MIE tables read in had the same nmuobs used
      !assuming that this atm does have a cloud associated with it
      IF (iCloudySky .GT. 0) THEN
        iLayers = 0
        LL = 0
        DO i=1,nscattab
          IF (iaScatTable_With_Atm(I).GT. 0) THEN
            iLayers = iLayers + nmuobs(i)
            LL = LL + 1  !!keep track of how many scattering tables used
            II = I       !!keep track of which scattering table used with atm
            END IF
          END DO
        IF (int(iLayers*1.0/LL) .NE. nmuobs(II)) THEN
          write (kStdErr,*) iLayers,LL,int(iLayers*1.0/LL),nmuobs(II)
          write (kStdErr,*) 'Some of the Mie Scattering tables had different'
          write (kStdErr,*) 'number of angles used in computing Mie coeffs'
          write (kStdErr,*) 'Please recheck  sscatmie.x and rerun'
          CALL DoStop
          END IF
      ELSE
        write (kStdWarn,*) 'no clouds with atmosphere number ',iAtm,' !!!'
        END IF

c Frank Evans code scales the Mie scattering parameters, so we have to 
c unscale them!!!!!!!!
c          DO I = 1, NSCATTAB 
c            IF (iaScatTable_With_Atm(I).GT. 0) THEN
c              CALL UnScaleMie(
c     $          caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
c     $          ndme(i)*nwavetab(i))
c              END IF
c            END DO              
       
c code from n_rad_jac_scat.f 
       iCloud = -1
       IF (iCloudySky .LT. 0) THEN
          write(kStdWarn,*)'Could not find a cloud for atmosphere #',iAtm
          write(kStdWarn,*)'setting IWP = -100.0'
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
            IF (iaCloudWithThisAtm(i) .EQ. 1) THEN
              ncldlay = ncldlay + iaCloudNumLayers(i)
              write(kStdWarn,*) 'Cloud #, num layers = ',i,iaCloudNumLayers(i)

              write(kStdWarn,*) 'L   iwp  dme  iscattab   kCARTA Layer  : ' 

              DO iStep = 1, iaCloudNumLayers(i)
                iLayers = iLayers + 1
                IWP(iLayers) = raaaCloudParams(i,iStep,1) 
                DME(iLayers) = raaaCloudParams(i,iStep,2) 
                ISCATTAB(iLayers) = iaaScatTable(i,iStep)
                write(kStdWarn,*) iLayers,iwp(iLayers),dme(iLayers),
     $              iscattab(iLayers),iaaCloudWhichLayers(i,iStep)
                END DO
              END IF
            END DO
          END IF

C     Find the levels for top of cloud and observation level
c     remember that these numbers are with respect to the KLAYERS pressure
c     levels and layers
c     these will be reset when passed in and out of GetAbsProfileDISORT
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
c when we call GetAbsProfile, variables icldtop,icldbot,iobs will be reset
c to reflect the rtspec layering
c   TOA    --------------
c             layer 1
c          --------------
c              .....                 DISORT
c          --------------
c        layer iNumLayer-1
c          --------------
c         layer iNumLayer
c   GND --------------------------

      IF (IWP(1) .LE. 0.0) THEN  !we have no cloud; set up fictitious clouds
        IF (iDownWard .GT. 0) THEN    
          !down look instr : set cloud BELOW observer, in kCARTA layer #1
          ICLDTOP = 2
          ICLDBOT = 1
          !down look instr : set cloud BELOW observer, in DISORT layer #2
          ICLDTOP = 2
          ICLDBOT = 3
          IOBS    = iNumLayer     
        ELSE IF (iDownWard .LT. 0) THEN    !up look instr
          !up look instr : set cloud ABOVE observer, in kCARTA layer #iNumLayer
          ICLDTOP = iNumLayer+1
          ICLDBOT = iNumLayer
          !up look instr : set cloud ABOVE observer, in DISORT layer #1
          ICLDTOP = 2
          ICLDBOT = 1
          IOBS    = 1             
          END IF
        END IF

      IF (IWP(1) .GT. 0.0) THEN
        !do not really need icldtop/bot, but just set it up
        ICLDTOP=iaaCloudWhichLayers(iCloudySky,1)+1
        ICLDBOT=iaaCloudWhichLayers(iCloudySky,iaCloudNumLayers(iCloudySky))

        icldbot = iT_Atm - icldbot + 1
        icldtop = iT_Atm - icldtop + 1

        icldbot = icldbot + 1
        icldtop = icldtop + 1

        IF (iDownWard .GT. 0) THEN
          IOBS   = iNumLayer       
        ELSE IF (iDownWard .LT. 0) THEN
          IOBS   = 1
          END IF
        END IF

      iobs=(iNumLayer+1)-iobs+1

 30   FORMAT(I3,' ',A80)

      RETURN
      END

c************************************************************************
c this subtroutine does some initializations for a uplook instr
      SUBROUTINE Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,rSatAngle,iTag,
     $                         raTopIntensity,raSolarBeam,TOA_to_instr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iTag
      REAL raVTemp(kMixFilRows)        !temperature profile
      REAL rFracTop                    !topmost fraction
      REAL raWaves(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs  
      REAL rSatAngle
c output variables
      REAL raTopIntensity(kMaxPts)     !intensity at TOA
      REAL raSolarBeam(kMaxPts)        !intensity diue to solar rad
      REAL TOA_to_instr(kMaxPts)       !abs coeffs

c local variables
      INTEGER iaRadLayer(kProfLayer),iL,iF
      REAL r1,r2,rDummy

      r1=kPlanck1
      r2=kPlanck2

      DO iL=1,kMaxPts
        TOA_to_instr(iL)=0.0
        END DO

c bring incident space radiation down from TOA to instrument
      DO iF=1,kMaxPts
c compute the Plank radiation from space
        raTopIntensity(iF)=exp(r2*raWaves(iF)/kTSpace)-1.0 
        raTopIntensity(iF)=r1*((raWaves(iF))**3)/raTopIntensity(iF)
        END DO 

c set the solar beam intensity ... 
c if the sun is ON and is in the FOV of the instrument, then the BC of 5600K
c   to the TOA is set, else the BC of 2.6K to TOA is set; also fbeam set to 0
c if the sun is ON and is NOT in the FOV of the instrument, then the fbeam is
c   set to solar beam, while BC of TOA is 2.6K is set

      DO iF=1,kMaxPts     !!!!assume sun is NOT ON
        raSolarBeam(iF)=0.0
        END DO

      IF (kSolar .GE. 0) THEN
        IF (abs(abs(rSatAngle)-abs(kSolarAngle)) .GE. 1.0e-3) THEN
          !sun is on, but not in FOV of instr
          !so set fbeam correctly to that of sun
          CALL SolarBeam(kSolar,raSolarBeam,raWaves,iTag) 
          rDummy = abs(cos(kSolarAngle*kPi/180))
c bring solar radiation down from TOA to instrument
          CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,TOA_to_instr) 
          DO iF=1,kMaxPts
            raSolarBeam(iF)=raSolarBeam(iF)*exp(-TOA_to_instr(iF)/rDummy)
            END DO
          END IF

        IF (abs(abs(rSatAngle)-abs(kSolarAngle)) .LE. 1.0e-3) THEN
          !sun is on, and in FOV of instr
          !so set fbeam correctly to 0.0, and the BC of 5600K
          END IF

        END IF     !!IF (kSolar .GE. 0)

      RETURN
      END
c************************************************************************
c this subtroutine does some initializations for a downlook instr
      SUBROUTINE Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,rSatAngle,iTag,
     $                         raTopIntensity,raSolarBeam,TOA_to_instr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iTag
      REAL raVTemp(kMixFilRows)        !temperature profile
      REAL rFracTop                    !topmost fraction
      REAL raWaves(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs  
      REAL rSatAngle
c output variables
      REAL raTopIntensity(kMaxPts)     !intensity at TOA
      REAL raSolarBeam(kMaxPts)        !intensity diue to solar rad
      REAL TOA_to_instr(kMaxPts)       !abs coeffs

c local variables
      INTEGER iaRadLayer(kProfLayer),iL,iF
      REAL r1,r2,rDummy

      r1=kPlanck1
      r2=kPlanck2

c see if there are any layers between TOA and instr; if there are, set
c TOA_to_instr to the cumulative k(TOA to aircraft), else set it to 0
      DO iL=1,kProfLayer
        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
        END DO
      CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,TOA_to_instr) 

c bring incident space radiation down from TOA to instrument
      DO iF=1,kMaxPts
c compute the Plank radiation from space
        raTopIntensity(iF)=exp(r2*raWaves(iF)/kTSpace)-1.0 
        raTopIntensity(iF)=r1*((raWaves(iF))**3)/raTopIntensity(iF)
        END DO 
c this is technically incorrect ... we really should do rad transfer here
      DO iF=1,kMaxpts
        raTopIntensity(iF)=raTopIntensity(iF)*exp(-TOA_to_instr(iF)) 
        END DO
      
      IF (kSolar .GE. 0) THEN
        CALL SolarBeam(kSolar,raSolarBeam,raWaves,iTag) 
        rDummy = abs(cos(kSolarAngle*kPi/180))
c bring solar radiation down from TOA to instrument
        DO iF=1,kMaxPts
         raSolarBeam(iF)=raSolarBeam(iF)*exp(-TOA_to_instr(iF)/rDummy)
          END DO
      ELSE 
        DO iF=1,kMaxPts
          raSolarBeam(iF)=0.0
          END DO
        END IF

        RETURN
        END 
c************************************************************************
c set raCorrelatedK to abs coeffs of layer nearest the ground
      SUBROUTINE SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,
     $                iaaRadLayer,iAtm,iNumLayer, ABSPROF)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

      INTEGER iDownWard            !direction of rad transfer
      REAL raCorrelatedK(kMaxPts)  !to see how the "k distributions" are
      REAL raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      REAL    ABSPROF(MAXNZ,MAXABSNU)      

      INTEGER iI,iJ,iaIndx(kMaxPts),iMethod
      REAL raL2S(kMaxPts),raaWgt(kMaxPts,kProfLayer),raPeak(kMaxPts)

      iMethod = -1

      IF (iMethod .EQ. -1) THEN
        !simply store optical depth of layer closest to gnd
        IF (iDownWard .EQ. 1) THEN
          !! as k(gnd) increases, down look instr only penetrates less into 
          !! atmosphere, so radiance decreases
          DO iI=1,kMaxPts
            raCorrelatedK(iI)=raaAbs(iI,iaaRadLayer(iAtm,1))
            END DO
        ELSE
          !! as k(gnd) increases, uplook instr only penetrates less into 
          !! atmosphere, so radiance increases as you see hot stuff in 
          !!your face
          DO iI=1,kMaxPts
            raCorrelatedK(iI)=raaAbs(iI,iaaRadLayer(iAtm,iNumLayer))
            END DO
          END IF
        END IF

      IF (iMethod .EQ. +1) THEN
        !store optical depth of layer that peaks the weight fcn
        !!!!!!remember ABSPROF(1,:) = TOA, ABSPROF(NLEV-1,:) = GND
        IF (iDownward .EQ. 1) THEN
          !do down look weight fcn
          DO iI=1,kMaxPts
            raL2S(iI)  = 0.0        !!!optical depth to space = 0.0
            iaIndx(iI) = 1
            raPeak(iI) = -1.0e10
            END DO
          DO iI = 1,kMaxPts
            DO iJ = iNumLayer,1,-1
              raaWgt(iI,iJ) = absprof(iNumLayer-iJ+1,iI)
              raaWgt(iI,iJ) = (1-exp(-raaWgt(iI,iJ)))*exp(-raL2S(iI))
              raL2S(iI)     = raL2S(iI) + absprof(iNumLayer-iJ+1,iI)
              END DO
            END DO
        ELSE IF (iDownward .EQ. -1) THEN
          !do up look weight fcn
          DO iI=1,kMaxPts
            raL2S(iI) = 0.0        !!!optical depth to gnd = 0.0
            iaIndx(iI) = 1
            raPeak(iI) = -1.0e10
            END DO
          DO iI = 1,kMaxPts
            DO iJ = 1,iNumLayer
              raaWgt(iI,iJ) = absprof(iNumLayer-iJ+1,iI)
              raaWgt(iI,iJ) = (1-exp(-raaWgt(iI,iJ)))*exp(-raL2S(iI))
              raL2S(iI)     = raL2S(iI) + absprof(iNumLayer-iJ+1,iI)
              END DO
            END DO
          END IF
        !!!having stored the wgt fcns, find where it peaks
        DO iI = 1,kMaxPts
          DO iJ = 1,iNumLayer
            IF (raaWgt(iI,iJ) .GT. raPeak(iI)) THEN
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
      END 

c************************************************************************
c set up the single scatter albedos etc for the clouds
      SUBROUTINE SetUpClouds(nstr, nmuobs, iaCloudWithThisAtm,
     $                iaCldTop,iaCldBot,iaCloudNumLayers,rF,
     $                IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $                TABEXTINCT, TABSSALB, TABASYM, ISCATTAB,
     $                !!!!!!!outputs
     $                extinct,dtauc,ssalb,asym,pmom)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c inputs
      REAL rF                                 !wavenumber
      INTEGER iaCloudWithThisAtm(kMaxClouds)  !is this cloud in this atm?
      INTEGER IACLDTOP(kMaxClouds)            !cloud top layer
      INTEGER IACLDBOT(kMaxClouds)            !cloud bot layer
      INTEGER iaCloudNumLayers(kMaxClouds)    !number of layers cloud occupies
      REAL    IWP(MAXNZ), DME(MAXNZ)          !iwp, particle size
      INTEGER NDME(MAXSCAT), NWAVETAB(MAXSCAT),ISCATTAB(MAXNZ)
      INTEGER nstr                            !number of disort streams
      INTEGER nmuobs                          !number of RTSPEC computed angles
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
c outputs
      REAL extinct                      !extinct coeff from Mie Scatter tables
      DOUBLE PRECISION dtauc(maxcly)    !optical depths; also used as INPUT
      DOUBLE PRECISION ASYM(maxnz)
      DOUBLE PRECISION pmom(0:maxmom,maxcly)  !scattering phase fcn
      DOUBLE PRECISION ssalb(maxcly)   !single scatter albedo of layers

c local variables
      INTEGER LL,M,ICLDTOP,ICLDBOT,N,L,I,nmom
      REAL ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz) 
      DOUBLE PRECISION tauC(kProfLayer),tauCG(kProfLayer) 

      !!!!!!!!!!!!!! ********* CLOUD SCATTERING ************* !!!!!!!!!!!!
      LL=0
      DO M = 1,kMaxClouds
        IF (iaCloudWithThisAtm(M) .GT. 0) THEN
          ICLDTOP = IACLDTOP(M)
          ICLDBOT = IACLDBOT(M)
          DO N = ICLDTOP, ICLDBOT - 1   
            !L is the cloud layer = 1(top)..nlay(bot) in cloud
            L = LL + N-ICLDTOP+1

c            write (kstdwarn,*) 'm,tp,bt,l,iwp(l),dme(l) = ',
c     $ m,icldtop,icldbot,l,iwp(l),dme(l)

            I = ISCATTAB(L)     !!!!!I is the scattering table info number

C      Interpolate to get values of extinction, s.s. albedo, and 
C      phi function values at given obs. mu, waveno, and particle size. 
C      Note: we don't need phi functions for Eddington-only solution
c      This means that while the rtspec code had a choice of 
c            CALL INTERP_SCAT_TABLE2 (rF, DME(L),    versus
c            CALL INTERP_SCAT_TABLE3 (rF, DME(L),     
c      here we only need the simpler first choice as we are not messing 
c      around with the phase functions
            CALL INTERP_SCAT_TABLE2 (rF, DME(L),     
     $                EXTINCT, SSALB_RTSPEC(L), ASYM_RTSPEC(L), 
     $                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), 
     $                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I)) 

c but the indices have to be modified so we put info into the correct layers
c also have to set asymmetry parameter
            TAUC(L) = DBLE(IWP(L)*EXTINCT/1000)
            TAUCG(L) = dtauc(N) + TAUC(L)
            IF (TAUCG(L) .GT. 1.0E-5) THEN 
              SSALB_RTSPEC(L) = SSALB_RTSPEC(L)*TAUC(L)/TAUCG(L)
            ELSE 
              SSALB_RTSPEC(L) = 0.0
              ENDIF 

c so now set the indices correctly and do relevant changes of real --> double
            dtauc(N) = DBLE(TAUCG(L))
            SSALB(N) = DBLE(SSALB_RTSPEC(L))
            IF (ssalb_rtspec(L) .gt. 1.0e-6) THEN
              asym(N)  = DBLE(asym_rtspec(L))
            ELSE
              asym(N)  = DBLE(0.0)
              END IF
  
            !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
            nmom = max(2*nmuobs + 1,nstr)
            nmom = min(maxmom,nmom)
            CALL getmom(3,asym(N),nmom,pmom(0,N))        !!!first arg = 3
            ENDDO                 !DO N=ICLDTOP,ICLDBOT-1
          LL = LL + iaCloudNumLayers(M)
          END IF                  !IF (iaCloudWithThisAtm(M) .GT. 0) THEN
        END DO                    !DO M = 1,kMaxClouds

      RETURN
      END

c************************************************************************
c this subroutine sets up the optical depths, single scatter albedo etc for
c a atmosphere where there is only gas + Rayleigh scattering
      SUBROUTINE SetUpRayleigh(nlev,nstr,nmuobs, rF,raDensity,raThickness,
     $                        dtauc,ssalb,asym,pmom)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'
c inputs
      INTEGER nlev,nstr          !number of levels, number of streams
      INTEGER nmuobs             !number of RTSPEC angles for Mie computations
      REAL rF                    !wavenumber
      REAL raDensity(kProfLayer),raThickness(kProfLayer)
c outputs
      DOUBLE PRECISION dtauc(maxcly)    !optical depths; also used as INPUT
      DOUBLE PRECISION ASYM(maxnz)
      DOUBLE PRECISION pmom(0:maxmom,maxcly)  !scattering phase fcn
      DOUBLE PRECISION ssalb(maxcly)   !single scatter albedo of layers

c local variables
      REAL ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
      DOUBLE PRECISION tauC(kProfLayer),tauCG(kProfLayer)
      REAL Rayleigh
      INTEGER N,nmom

      DO N = 1,NLEV-1   !!!!!!!to include scattering
        !indices have to be modified so we put info into the correct layers
        TAUC(N) = DBLE(rayleigh(rF,raDensity(N),raThickness(N)))
        TAUCG(N) = dtauc(N) + TAUC(N)
        SSALB_RTSPEC(N) = SNGL(TAUC(N)/TAUCG(N))
c so now set the indices correctly and do relevant changes of real --> double
        dtauc(N) = TAUCG(N)
        SSALB(N) = DBLE(SSALB_RTSPEC(N))
        asym(N)  = DBLE(0.0)

        !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
        nmom = max(2*nmuobs + 1,nstr)
        nmom = min(maxmom,nmom)
        CALL getmom(2,asym(N),nmom,pmom(0,N))         !!!first arg = 2
        ENDDO 

      RETURN
      END

c************************************************************************
c this does the final initializations before calling DISORT
      SUBROUTINE FinalInitialization(
     $    !!!!inputs
     $       iDownWard,rSatAngle,rTopIntensity,rSolarBeam,emiss,
     $       rSurfaceTemp,dtauc,dTotalOpticalDepth,iDoFlux,nlev,
     $       iNp,iaOp,
     $    !!!!outputs
     $       usrtau,ntau,utau,usrang,numu,umu,
     $       nphi,phi,fisot,fbeam,umu0,phi0,
     $       ibcnd,lamber,albedo,
     $       btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'
 
c input variables
      INTEGER nlev,iNp,iaOp(kPathsOut)
      DOUBLE PRECISION dtauc(maxcly)    !optical depths; also used as INPUT
      DOUBLE PRECISION dTotalOpticalDepth 
      REAL rSatAngle,rTopIntensity,rSolarBeam,emiss,rSurfaceTemp
      INTEGER iDownWard,iDoFlux
c output variables      
      INTEGER ntau,numu,nphi
      LOGICAL usrtau,usrang
      DOUBLE PRECISION utau(maxulv)     !tau's at which to output results 
      DOUBLE PRECISION umu(maxumu)      !ang's at which to output results 
      DOUBLE PRECISION phi(maxphi)      !azimuthal phi's to output radiance
      DOUBLE PRECISION fisot            !isotropic TOA radiation
      DOUBLE PRECISION fbeam,umu0,phi0  !solar beam

      INTEGER ibcnd
      LOGICAL lamber,plank,onlyfl
      DOUBLE PRECISION albedo,btemp,ttemp,temis,accur
      CHARACTER  header*127             !dumb comment
      LOGICAL prnt(5)                   !prnt(1)=true, print input variables

c local variables   
      INTEGER iI,iJ
      DOUBLE PRECISION d1,d2,dCumulative(maxcly)
     
      d1 = dtauc(1)
      d2 = dTotalOpticalDepth

      IF (iDownward .EQ. 1) THEN
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
        umu(1) = DBLE(ABS(cos(rSatAngle*kPi/180)))
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
        umu(1) = DBLE(-ABS(cos(rSatAngle*kPi/180)))
        nphi = 1          !return intensity at one user zenith angle
        phi(1) = DBLE(0.0)
        END IF

      fisot  = DBLE(rTopIntensity)

      IF (abs(abs(rSatAngle) - abs(kSolarAngle)) .LE. 1e-7) THEN
        kSolarAngle=kSolarAngle+1e-4
        END IF
      fbeam  = DBLE(rSolarBeam)
      umu0   = DBLE(cos(kSolarAngle*kPi/180))
      phi0   = DBLE(0.0)

      ibcnd = 0           !general case
      lamber = .false.    !specify bidir reflectance
      lamber = .true.     !isotropic reflect lower bdry ==> specify albedo
      albedo  = DBLE(1.0-emiss)   !from defns in books,papers

      btemp  = DBLE(rSurfaceTemp)    !!ground temperature
      IF (iDownward .EQ. 1) THEN
        !down look instrument
        ttemp  = DBLE(kTSpace)         !!2.7 Kelvin
      ELSE
        !up look instrument
        IF (kSolar .LT. 0) THEN
          ttemp  = DBLE(kTSpace)         !!2.7 Kelvin or 5600K!!!!!!
        ELSE IF (kSolar .GE. 0) THEN
          IF (abs(abs(kSolarAngle)-abs(rSatAngle)) .LE. 1.0e-3) THEN
            ttemp  = DBLE(5600.0)         !!sun in FOV .. 5600K!!!!!!
          ELSE
            ttemp  = DBLE(kTSpace)        !!sun not in FOV .. 2.7!!!!!!
            END IF
          END IF
        END IF

      temis  = DBLE(1.0)
      plank  = .true.      !emission from layers

      IF (iDoFlux .EQ. -1) THEN
        onlyfl = .false.     !only need intensity as well
        header  = 'scattering computations'
      ELSEIF (iDoFlux .EQ. +1) THEN
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
      END

c************************************************************************
c this is the main interface to DISORT
      SUBROUTINE interface_disort(
        !first the usual kCARTA variables
     $        raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
c     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        iNp,iaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
c     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
c     $        raLayAngles,raSunAngles,
     $        raLayAngles,iDoFlux,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable, 
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity) 

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iDoFlux     = do flux (+1) or radiance (-1)
c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
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
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c iDownward = +1 ==> downward looking instrument
c             -1 ==> upward looking instrument
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raNumberDensity(kProfLayer)
      REAL raLayAngles(kProfLayer)
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL rFracTop,rFracBot
      REAL raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag,iDoFlux
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
      CHARACTER*80 SCATFILE(MAXSCAT)

      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

      INTEGER NSCATTAB, NCLDLAY, NLEV, NABSNU
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      REAL    IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iCloudySky
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)

c ---------------------------------------------------------------------

C   RTSPEC Radiative transfer variables: 
      REAL    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)

c ---------------------------------------------------------------------

c these are direct from DISOTEST.F

      CHARACTER  header*127          !dumb comment
      
      LOGICAL lamber  !true  ==> isotropic reflecting lower bdry
                      !          so need to specify albedo
                      !false ==> bidirectionally reflecting bottom bdry
                      !      ==> need to specify fcn BDREF()
      LOGICAL plank   !use the plank function for local emission
      LOGICAL onlyfl  !true ==> return flux, flux divergence, mean intensitys
                      !falsetrue ==> return flux, flux divergence, mean 
                      !              intensities and intensities
      LOGICAL prnt(5) !prnt(1)=true, print input variables
                      !prnt(2)=true, print fluxes
                      !prnt(3)=true, intensities at user angles and levels
                      !prnt(4)=true, print planar transmitivity and albedo
                      !prnt(5)=true, print phase function moments PMOM
      LOGICAL usrtau  !false ==> rad quantities return at every bdry
                      !true  ==> rad quantities return at NTAU optical depths
                      !          which will be specified in utau(1:ntau)
      LOGICAL usrang  !false ==> rad quantities returned at stream angles
                      !true  ==> rad quantities returned at user angles
                      !          which will be specified in umu(1:numu)

      INTEGER ibcnd            !0 ==> general case beam (fbeam), isotropic
                               !      top illumination (fisot), thermal top 
                               !      emission (temis,ttemp),internal thermal
                               !      sources (temper), reflection at bottom
                               !      (lamber, albedo, bdref), thermal 
                               !      emission from bottom (btemp)
                               !1 ==> return only albedo/trans of entire
                               !      medium vs incident beam angle
      INTEGER nmom             !number of legendre phase polynoms
      INTEGER nlyr             !number of computational layers
      INTEGER nstr             !number of radiation streams
      INTEGER ntau             !associated with LOGICAL usrtau, print results
                               !at this many optical depths
      INTEGER numu             !associated with LOGICAL usrang, specifies how
                               !many polar angles results to be printed (umu)
      INTEGER nphi             !specifies how many azimuth angles results to 
                               !be printed (phi) .. can only be 0 if 
                               !onlyfl = .true. 

      DOUBLE PRECISION accur   !accuracy convergence criterion for azimuth 
                               !(fourier cosine) series .. set between 0-0.01
      DOUBLE PRECISION albedo  !bottom bdry albedo
      DOUBLE PRECISION btemp   !bottom surface temp
      DOUBLE PRECISION dtauc(maxcly)       
                               !optical depths of computational layers
      DOUBLE PRECISION fisot   !intensity of top bdry isotropic illumination
      DOUBLE PRECISION fbeam   !intensity of incident // beam at TOA
      DOUBLE PRECISION phi(maxphi)    
                               !the azimuthal phi's to output radiance      
      DOUBLE PRECISION pmom(0:maxmom,maxcly)  
                               !scattering phase fcn in legendre polynoms
      DOUBLE PRECISION phi0    !solar azimuth
      DOUBLE PRECISION ssalb(maxcly)   
                               !single scatter albedo of computational layers
      DOUBLE PRECISION temper(0:maxcly) !temperature of the levels (0=TOA)
      DOUBLE PRECISION temis            !emissivity of top bdry
      DOUBLE PRECISION ttemp            !temperature of top bdry
      DOUBLE PRECISION wvnmhi, wvnmlo   !bounds within which to do computations
      DOUBLE PRECISION umu(maxumu)      !ang's at which to output results 
      DOUBLE PRECISION umu0         !polar angle cosine of incident solar beam
      DOUBLE PRECISION utau(maxulv)     !tau's at which to output results 

      !!!!!output variables
      DOUBLE PRECISION dfdt(maxulv)  !flux diverge d(net flux)/d(optical depth)
                                     !where 'net flux' includes direct beam
      DOUBLE PRECISION flup(maxulv)  !diffuse up flux
      DOUBLE PRECISION rfldir(maxulv)!direct beam flux (without delta scaling)
      DOUBLE PRECISION rfldn(maxulv) !diffuse down flux = total-direct beam
      DOUBLE PRECISION uavg(maxulv)  !mean intensity (including direct beam)
      DOUBLE PRECISION UU( MAXUMU, MAXULV, MAXPHI )
                                      !intensity if ONLYFL = .false., 0 o/w
      DOUBLE PRECISION albmed(maxumu)!albedo of medium as fcn of cos(umu(*))
                                      !only set if ibcn == 1
      DOUBLE PRECISION trnmed(maxumu)!transmission as fcn of cos(umu(*))
                                      !only set if ibcn == 1
c ---------------------------------------------------------------------
      DOUBLE PRECISION dTotalOpticalDepth
      DOUBLE PRECISION ASYM(maxnz)
c ---------------------------------------------------------------------

c these variables are to get the parameters from Frank Evans Mie Code 
      REAL ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)

c this is to do "correlated k" computations
      REAL raCorrelatedK(kMaxPts)

c more local variables
      INTEGER iaStep(kMaxPts),iDiv,iScatter
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iScat,iIOUN,iF,iFF,iFFMax,iL
      REAL TOA_to_instr(kMaxPts)
      INTEGER iaRadLayer(kProfLayer)
      REAL raSolarBeam(kMaxPts),raTopIntensity(kMaxPts)
      REAL r1,r2,rDummy
      INTEGER iReadTable,iII
      INTEGER iDumper,iRayleigh
      REAL raDensity(kProfLayer)
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle

c these parameters are to step thru some of the 10000 pts
      REAL raaIntenSTEP(kProfLayer,kMaxPts),raIntenSTEP(kMaxPts)
      REAL raaNoscatterSTEP(kProfLayer,kMaxPts),raNoscatterSTEP(kMaxPts)
      REAL rakSTEP(kMaxPts),raFreqSTEP(kMaxPts),radtot
      INTEGER iSTEP,iStepPts

c these variables are to get the parameters from Frank Evans Mie Code 
      REAL extinct
      INTEGER LL,L,I,N,M
      
c this is fro RAYLEIGH
      REAL raThicknessRayleigh(kProfLayer)

c this is to convert from W to mW
      REAL raInten2(kMaxPts)

      rSolarAngle = kSolarAngle

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c  kScatter = 1 works good for both up and down look
c  kScatter = 3 works good for up look
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      r1=kPlanck1
      r2=kPlanck2

      !!! Rayleigh scatter has hardly any effect at infrared wavenumbers
      iRayleigh = -1       !!! -1 : do not want Rayleigh Scattering
                           !!! +1 : do Rayleigh Scattering

      iIOUN=kStdkCarta

      DO iIn=1,kMaxPts
        raIntenSTEP(iIn) = 0.0
        raFreqSTEP(iIn)  = 0.0
        rakSTEP(iIn)     = 0.0
        END DO

      WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',cos(rSatAngle*kPi/180.0),
     $   rSurfaceTemp

      CALL SetMieTables_DISORT(           
     $        !!!!!!!!!!!!!!!!!these are the input variables
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable, 
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        !!!!!!!!!!!!!!!!!!these are the output variables
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,
     $        TABPHI2UP, TABPHI2DN,
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, 
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,
     $        iCloudySky, IACLDTOP, IACLDBOT)

      CALL GetAbsProfileDISORT(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity,
     $      raDensity,raLayerTemp,
     $      iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)

      !!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
      IF (iCloudySky .LT. 0) THEN
        !!!!note that we do not care about background thermal accurately here
        write (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
        DO iii = 1,iNp
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
              END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $           rSatAngle,rSolarAngle,rSurfaceTemp,
     $           raUseEmissivity(iF),raWaves(iF),raInten(iF),iaOp(iii),
     $           iProfileLayers,raPressLevels)
            END DO
          CALL wrtout(iIOUN,caOutName,raWaves,raInten)
          END DO
        GOTO 9876
        END IF

      !!!!!!!!!! if iCloudSky .GT. 0 go thru and do the DISORT stuff !!!!!!!!
      CALL SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,
     $                    iaaRadLayer,iAtm,iNumLayer,ABSPROF)

c set up some things for the instrument
      IF (iDownward .EQ. 1) THEN             !!down look instr
        CALL Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,rSatAngle,iTag,
     $                         raTopIntensity,raSolarBeam,TOA_to_instr)
      ELSE
        CALL Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp,
     $                         rFracTop,raWaves,raaAbs,rSatAngle,iTag,
     $                         raTopIntensity,raSolarBeam,TOA_to_instr)
        END IF

c set the temperature levels, changing to double precision
      DO iL=1,NLEV
        temper(iL-1)=DBLE(temp(iL))
        END DO

c **************************** SPEED UP CODE ***************************
      nstr  = kDis_nstr   ! number of streams used by DISORT (2,4,8,16 ...)
      iStep = kDis_pts    ! number of wavenumber pts to use (1,2,...,10000)
                          ! out of 10000
c **************************** SPEED UP CODE ***************************
      IF (iStep .GT. kMaxPts) THEN
        write(kStdWarn,*) 'Resetting kDis_Pts to kMaxPts'
        iStep = kMaxPts
        END IF

      IF (iStep .LT. 20) THEN
        write(kStdWarn,*) 'Resetting kDis_Pts to 20'
        iStep = 20
        END IF

      !if you want to do 10 pts out of 10000 pts, then you have to do rad
      !transfer on points 1,1001,2001,3001,...10000
      !ie step over 10000/iStep points
      IF (kScatter .NE. 2) THEN
        iStep = iDiv(kMaxPts,iStep)  
        END IF

c set up array of wavenumber indices that we step thru, to do the rad transfer
c (as DISORT is quite slow, we will not use each and every point)
      iScatter = kScatter
      CALL FindWavenumberPoints(iStep,iScatter,raCorrelatedK,iaStep,iFFMax)

      iStepPts = 0
      DO iFF = 1,iFFMax
        IF = iaStep(iFF)
        iStepPts = iStepPts + 1
        DO iL=1,NLEV-1
          dtauc(iL) = DBLE(absprof(iL,iF))
          raTau(iL) = absprof(iL,iF)
          END DO

        wvnmlo = DBLE(raWaves(iF))
        wvnmhi = DBLE(raWaves(iF)+kaFrStep(iTag))

c ++++++++++++++++ this is to include MIE data from RTSPEC +++++++++++++++++

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

c to test no scattering, just replace following doloops with DO N = 1,-1

        IF (iRayleigh .EQ. -1) THEN     !want cloud, not Rayleigh scattering
          CALL SetUpClouds(nstr,nmuobs(1),iaCloudWithThisAtm,
     $                iaCldTop,iaCldBot,iaCloudNumLayers,raWaves(iF),
     $                IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $                TABEXTINCT, TABSSALB, TABASYM, ISCATTAB,
     $                extinct,dtauc,ssalb,asym,pmom)
        ELSE IF (iRayleigh .EQ. +1) THEN   !want Rayleigh, not cloud scattering
          CALL SetUpRayleigh(nlev,nstr,nmuobs(1),raWaves(iF),raDensity,
     $                       raThicknessRayleigh,dtauc,ssalb,asym,pmom)
          END IF

        dTotalOpticalDepth=DBLE(0.0)
        DO iL=1,NLEV-1
cdebug sergio ############### this tests NO scattering
c          raTau(iL) = absprof(iL,iF)
c          dtauc(iL) = DBLE(absprof(iL,iF))
c          asym(iL)  = 0.0
c          ssalb(iL) = 0.0
c          !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
c          nmom = max(2*nmuobs(1) + 1,nstr)
c          nmom = min(maxmom,nmom)
c          CALL getmom(3,asym(iL),nmom,pmom(0,iL))        !!!first arg = 3
cdebug sergio ###############
          dTotalOpticalDepth = dTotalOpticalDepth+dtauc(iL)
          END DO

c ++++++++++++++++ this is to do non scat rad transfer ++++++++++++++++++++
c ++++ use wavenumber at raWaves(iF)
c ++++ save results in raaNoScatterSTEP(iii,iStepPts)
      DO iii=1,iNp
        CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $       rSatAngle,rSolarAngle,rSurfaceTemp,raUseEmissivity(iF),
     $       raWaves(iF),raaNoScatterSTEP(iii,iStepPts),iaOp(iii),
     $       iProfileLayers,raPressLevels)
        END DO

c ++++++++++++++++ final initializations ++++++++++++++++++++++++++++++
c     note we do not need flux here!!!!
      CALL FinalInitialization(
     $       iDownWard,rSatAngle,raTopIntensity(iF),raSolarBeam(iF),
     $       raUseEmissivity(iF),rSurfaceTemp,
     $       dtauc,dTotalOpticalDepth,iDoFlux,nlyr+1,iNp,iaOp,
     $       usrtau,ntau,utau,usrang,numu,umu,
     $       nphi,phi,fisot,fbeam,umu0,phi0,
     $       ibcnd,lamber,albedo,
     $       btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)

c not sending MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM
c      CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
c     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
c     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
c     &                   FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
c     &                   PLANK, ONLYFL, ACCUR, PRNT, HEADER, MAXCLY,
c     &                   MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, RFLDN,
c     &                   FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )
     
      !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
      nmom = max(2*nmuobs(1) + 1,nstr)
      nmom = min(maxmom,nmom)

c      iL=NLEV-1
c      print *, 'cos(rSatAngle),sfctemp = ',cos(rSatAngle*kPi/180.0),
c     $   rSurfaceTemp
c      print *, 'sfc lamber, alb, emiss = ',lamber,albedo,raUseEmissivity(iF)
c      print *, 'nlyr,iL,dtauC = ',nlyr,iL,dtauC(iL)
c      print *,'  btemp,ttemp,temis,plank,onlyfl,accur,prnt,header = '
c      print *,   btemp,ttemp,temis,plank,onlyfl,accur,prnt,header
c      stop

      CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO,
     &                   WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU,
     &                   UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0,
     &                   FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS,
     &                   PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN,
     &                   FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

        DO iii = 1,iNp
          raKStep(iStepPts)     = raCorrelatedK(iF)
          raFreqSTEP(iStepPts)  = raWaves(iF)
          raaIntenSTEP(iii,iStepPts) = SNGL(uu(1,iii,1))
          END DO
 
        ENDDO
c ------------------------ END OF DISORT ---------------------------------

      !interpolate coarse step grid onto fine kCARTA grid   
      DO iii = 1,iNp
        DO iF = 1,iStepPts
          raIntenStep(iF)     = raaIntenStep(iii,iF)
          raNoScatterStep(iF) = raaNoScatterStep(iii,iF)
          END DO
        CALL Interpolator(raFreqStep,rakStep,raIntenStep,raNoScatterSTEP,
     $                 iStepPts,iDownWard,nlev,
     $                 iProfileLayers,raPressLevels,
     $                 raCorrelatedK,raLayerTemp,absprof,raWaves,
     $                 rSatAngle,rSolarAngle,rSurfaceTemp,raUseEmissivity,
     $                 raInten)
        CALL wrtout(iIOUN,caOutName,raWaves,raInten)
        END DO
      
 9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm

 900  FORMAT(8(E12.4,'  '))
 999  FORMAT(I2,'  ',I6,'   ',F9.4,'  ',F8.6,'  ',F8.4,'  ',3(E12.6,'  '))

      RETURN
      END

c************************************************************************
c this subroutine figures out how to do the interpolations
      SUBROUTINE Interpolator(raFreqStep,rakStep,raIntenStep,raNoScatterSTEP,
     $                   iStepPts,iDownWard,nlev,
     $                   iProfileLayers,raPressLevels,
     $                   raCorrelatedK,raLayerTemp,absprof,raWaves,
     $                   rSatAngle,rSolarAngle,rSurfaceTemp,raUseEmissivity,
     $                   raInten)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c inputs
      REAL raPressLevels(kProfLayer+1)
      REAL raFreqStep(kMaxPts),raKstep(kMaxPts)
      REAL raNoScatterSTEP(kMaxPts),raIntenStep(kMaxPts)
      REAL raCorrelatedK(kMaxPts),raWaves(kMaxPts),raLayerTemp(kProfLayer)
      REAL ABSPROF(MAXNZ,MAXABSNU),rSatAngle,rSolarAngle,rSurfaceTemp
      REAL raUseEmissivity(kMaxPts)
      INTEGER iStepPts,nlev,iDownWard,iProfileLayers
c output
      REAL raInten(kMaxPts)
      
c local variables
      REAL DWA(kMaxPts),raTau(kProfLayer),raNoScatterAll(kMaxPts)
      INTEGER Indx(kMaxPts),iF,iL
      REAL radtot,ttorad

      IF (kScatter .EQ. 1) THEN
        !!!interpolate wavenumber points
        IF (iDownWard .EQ. -1) THEN      !!!works better for uplook instr
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
              END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $         rSatAngle,rSolarAngle,rSurfaceTemp,
     $         raUseEmissivity(iF),raWaves(iF),raNoScatterAll(iF),-1,
     $         iProfileLayers,raPressLevels)
            END DO
          CALL INTERP_PLANCK_0(
     $          raFreqStep,raIntenStep,raNoScatterStep,iStepPts,
     $          raWaves,raNoScatterAll,raInten)
        ELSEIF (iDownWard .EQ. 1) THEN   !!!works better for downlook instr
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
              END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $         rSatAngle,rSolarAngle,rSurfaceTemp,
     $         raUseEmissivity(iF),raWaves(iF),raNoScatterAll(iF),-1,
     $         iProfileLayers,raPressLevels)
            END DO
          CALL INTERP_PLANCK_0(
     $          raFreqStep,raIntenStep,raNoScatterStep,iStepPts,
     $          raWaves,raNoScatterAll,raInten)
          END IF
        END IF

      IF ((kScatter .GE. 2) .AND. (iDownWard .EQ. -1)) THEN
        !!!first interpolate in k
        DO iF = 1,kMaxPts
          CALL Linear_Real(raKSTEP,raIntenStep,iStepPts,
     $                     raCorrelatedK(iF),raInten(iF))
          END DO

      ELSEIF ((kScatter .GE. 2) .AND. (iDownWard .EQ. 1)) THEN
        !!!interpolate in k
        !!!do cumulative distr function to see if things might be too spiky
        DO iF=1,kMaxPts
          DO iL=1,NLEV-1
            raTau(iL)  = absprof(iL,iF)
            END DO
          CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $       rSatAngle,rSolarAngle,rSurfaceTemp,
     $       raUseEmissivity(iF),raWaves(iF),raNoScatterAll(iF),-1,
     $       iProfileLayers,raPressLevels)
          END DO
        CALL INTERP_PLANCK_3(
     $          raFreqStep,raKStep,raIntenStep,raNoScatterStep,iStepPts,
     $          raWaves,raNoScatterAll,raCorrelatedK,raInten)

        END IF

      RETURN
      END
c************************************************************************
c this subroutine takes in the input abs coeffs (raaAbsCoeff) where raa(1,:)
c is the lowest layer and raa(kProfLayer,:) is the highest layer .. it then 
c outputs these  abs coeffs into absprof, where absprof(1,:) is the top, and
c absprof(iNumLayer,:) = ground

c remember kCARTA has -----------------
c                      layer 100 = TOA
c                     -----------------
c                    
c                           ...
c
c                     -----------------
c                      layer 1 = GND
c                     -----------------
c for a DOWNLOOK instr, the layering is set as 1,2,3,4,5,......,100
c   rFracTop set for the top layer (100) and rFracBot set for the bot layer (1)
c for a   UPLOOK instr, the layering is set as 100,99,98,........,1
c   rFracTop set for the top layer (100) and rFracBot set for the bot layer (1)
c and this is WHAT IS ALWAYS set prior to the rad transfer computations

c remember DISORT has -----------------
c      and RTSPEC      layer 1 = TOA
c                     -----------------
c                    
c                           ...
c
c                     -----------------
c                      layer 100 = GND
c                     -----------------
c so for an downlook instr, everything is easily coded below (ie simply 
c flip  1 -->  100, 100 --> 1,
c flip  2 -->   99,  99 --> 2,
c flip  3 -->   98,  98 --> 3 etc)


C     sets optical depths for NLEV-1 layers and NABSNU wavenumbers. 
C     The temperature (K) of the profile is also returned. (no height needed)

c same as GetAbsProfile in scatter_rtspec_main,f except that
c 1) calls getverticaltempDISORT
c 2) sets the particle number denisties of the layers
      SUBROUTINE GetAbsProfileDISORT(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity,
     $      raDensity,raLayerTemp,
     $      iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param' 

c these are variables that come in from kcartamain.f 
      REAL raaAbs(kMaxPts,kMixFilRows),raWaves(kMaxPts),rFracTop,rFracBot
      REAL raVTemp(kMixFilRows),rSurfaceTemp
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNpmix
      INTEGER iDownWard,iProfileLayers
      REAL iwp, raNumberDensity(kProfLayer),raLayerTemp(kProfLayer)
c these are variables that we have to set
      INTEGER  NABSNU, NLEV         !!!!!!!!!!!!!MAXNZ, MAXABSNU 
      REAL   TEMP(*), ABSPROF(MAXNZ,*) 
      INTEGER  ICLDTOP,iCLDBOT,IOBS
      REAL raDensity(*)
      REAL raThicknessRayleigh(kProfLayer)

c local variables
      INTEGER iaRadLayer(kProfLayer), iaTemp(kProfLayer), iFr, iL, iLay 
      REAL NU, raVT1(kMixFilRows), InterpTemp

c these are to flip the temperature, abs profiles if instr looks up
      REAL raTemp(kProfLayer+1),raaTemp(kProfLayer,kMaxPts)

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

      DO iFr=1,kMixFilRows 
        raVT1(iFr)=raVTemp(iFr) 
        END DO 

c set the temperature of the bottommost layer correctly 
      IF (iDownWard .EQ. 1) THEN         !downlook instr
c if the bottommost layer is fractional, interpolate!!!!!! 
        iL=iaRadLayer(1) 
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,
     $                       1,iL) 
        write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!! 
c this is hardly going to affect thermal/solar contributions (using this temp  
c instead of temp of full layer at 100 km height!!!!!! 
        iL=iaRadLayer(iNumLayer) 
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,
     $                       -1,iL) 
        write(kStdWarn,*) 'top temp interped to ',raVT1(iL) 
      ELSE IF (iDownWard .EQ. -1) THEN       !uplook instr
c if the bottom layer is fractional, interpolate!!!!!! 
        iL=iaRadLayer(iNumLayer) 
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,
     $                       1,iL) 
        write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the top layer is fractional, interpolate!!!!!! 
        iL=iaRadLayer(1) 
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,
     $                       -1,iL) 
        write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 
        END IF

c now set the kCARTA LAYERS temperatures, used with NoScatterRadTransfer
c recall for DISORT , toa = layer 1   while for kCARTA, toa = 100
c recall for DISORT , gnd = layer 100 while for kCARTA, gnd = 1
c but also that for UPLOOK instr, clear sky kCARTA flips the layering!!!!
      IF (iDownward .EQ. 1) THEN
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
      IF (iDownWard .EQ. -1) THEN
        DO iLay = 1, iNumLayer
          iaTemp(iLay)=iaRadLayer(iLay)
          END DO
        DO iLay = 1, iNumLayer
          iaRadLayer(iNumLayer-iLay+1)=iaTemp(iLay)
          END DO
        END IF

c set the vertical temperatures of the atmosphere 
c instead of temperatures at n layers, it computes temps at (n+1) levels
      CALL SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,
     $                    rSurfaceTemp,iProfileLayers,raPressLevels)

c now set up the abs coeffs
c initialize array to all zeroes
c      print *,'rFracTop,rFracBot = ',rFracTop,rFracBot
      DO iFr=1,kMaxPts
        DO iLay=iNumLayer+1,kProfLayer
          absprof(iLay,iFr)=0.0
          raDensity(iLay)=0.0
          raThicknessRayleigh(iLay)=0.0
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
         raDensity(iNumLayer-iLay+1)=raNumberDensity(iL)
         !!!!!Laythick is from outlayers.param    before Nov 17,2001
         !!!!!Laythick is gotten from the KLAYERS profile after Nov 17, 2001
         !!!!!only need for Rayleigh scattering, which is negligible in infrared
         raThicknessRayleigh(iNumLayer-iLay+1)=raThickness(iL)*nu
         DO iFr=1,kMaxPts 
           !absprof wants level 1 == TOA, level iNumLayer= gnd
           absprof(iNumLayer-iLay+1,iFr)=raaAbs(iFr,iL)*nu
           END DO 
         END DO 

c comment this out on Oct 20, 2000 as is taken care of in interface_disort
c now set iobs
c      iobs=(iNumLayer+1)-iobs+1

c iDownward = +1 ==> downward looking instrument
c             -1 ==> upward looking instrument
c remember there is ONE more level than there are layers
c      icldtop=(iNumLayer+1)-icldtop+1
c      icldbot=(iNumLayer+1)-icldbot+1

      RETURN 
      END
c************************************************************************

c set the vertical temperatures of the atmosphere 
c this sets the temperatures at the pressure level boundaries, using the
c temperatures of the pressure layers that have been supplied by kLayers
c note that temeprature of bottom level is NOT set to rSurfaceTemp
      SUBROUTINE SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,
     $                         rSurfaceTemp,iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c these are variables that come in from kcartamain.f 
      REAL raVTemp(kMixFilRows),rSurfaceTemp,raPressLevels(kProfLayer+1)
      INTEGER iNumLayer,iProfileLayers
c these are variables that we have to set
      REAL    TEMP(*)

c local variables
      INTEGER iaRadLayer(kProfLayer), iL, iLay ,iM, idiv
      REAL FindBottomTemp,Temp1(maxnz)
      REAL pavg(kProfLayer),rP,raProfileTemp(kProfLayer)

      DO iLay=1,MAXNZ
        Temp1(iLay)=0.0
        END DO

c see which set of Mixed Paths the current atmosphere occupies eg 
c set 1 = 1..100, set2= 101..200 etc
c eg if current atmosphere is from MixfilPath 110-190, and kProfLayer = 100,
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
          iL=kProfLayer
          END IF
        rP=raPressLevels(iL+1)-10000*delta
        if (rp .LT. raPressLevels(kProfLayer+1)) then
          rp = raPressLevels(kProfLayer+1)+10000*delta
          end if
        TEMP1(iNumLayer-iLay+1)=FindBottomTemp(rP,raProfileTemp,
     $                                         raPressLevels,iProfileLayers)
        END DO

cccc this is where we differ from RTSPEC
cccc      TEMP1(iNumLayer+1)=rSurfaceTemp
      rP = DISORTsurfPress
      TEMP1(iNumLayer+1)=FindBottomTemp(rP,raProfileTemp,
     $                                         raPressLevels,iProfileLayers)

      DO iLay=1,iNumLayer+1
        temp(iLay)=temp1(iLay)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine is from LBLRTMv5.10

      REAL FUNCTION Rayleigh(fr,wtot,layerthick)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'
      
C     ********** Rayleigh Scattering calculation ********** 
c 
c     The formulation, adopted from MODTRAN_3.5 (using approximation 
c     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization 
c     = 0.0279, output in km-1 for T=273K & P=1 ATM) is as follows: 
C 
c     The rayleigh "absorption" coefficient (scattering in the direct 
c     beam) ray_abs can be defined as 
c 
c         ray_abs = (fr**4/(9.38076E18-1.08426E09*fr**2)) 
c     *        *wmol_tot*conv_cm2mol 
c 
c     where fr is the wavenumber value, wmol_tot is the total 
c     molecular amount in the layer, and conv_cm2mol is the conversion 
c     factor derived by multiplying air density (2.68675E19 mol/cm3) 
c     at 273 K with the number of km per cm (1.e-5 km/cm). 
c 
c     For numerical purposes, the layer amount of all molecules is 
c     calculated above as WTOT, which has been scaled by 1.e-20. We 
c     have taken 1.e20 out of the air density in the denominator 
c     as well. In addition, a factor of 10,000 (cm-1) has been  
c     divided out of fr. Finally, the radiation field is 
c     excluded, so xfr**4 is replaced by xfr**3. When 
c     JRAD=1, the radiation field is put in by multiplying the 
c     absorption by xfr. 
c 
c     Rayleigh scattering in the direct beam is only calculated for 
c     model runs > 3100 cm-1. 
c 
c        Thus the current formulation is 
 
      REAL Losch             !Loschmidt number
      REAL layerthick        !layer thickness
      REAL fr,xrayl,wtot     !fr is the wavenumber
                             !wtot is total number of molecules 
                             !in layer and xrayl is multiplication factor
       
      REAL conv_cm2mol,xfr
      REAL l,l1,a0,a1,a2,a3,y0,y1,y2,y3

      xrayl = 1.0
      conv_cm2mol = 1./(2.68675e-1*1.e5) 

      l=10000.0/fr        !change from cm-1 to um
          
c      xfr = fr/1.e4 
c      y0 = (xfr**4/(9.38076E2-10.8426*xfr**2))*wtot*conv_cm2mol*1e-20
c      y0 = y0*100              !change abs coeff from cm-1 to m-1
c      y0 = y0*layerthick
cc LBLRTM had xfr**3 and not xfr**4
cc      y0 = (xfr**3/(9.38076E2-10.8426*xfr**2))*wtot*conv_cm2mol*XRAYL 


c works quite well in the 300-3000 cm-1 range! use this! 
c this is from SBDART 
c calculate molecular rayleigh scattering coefficient 
c using approximation of shettle et al., 1980 (appl. opt., 2873-4) 
c with the depolarization = 0.0279 instead of 0.035 
c input: fr = frequency in wavenumbers (cm-1) 
c output: y0 = molecular scattering coefficient (km-1) 
c               for temperature = 273 k & pressure = 1 atm. 
      Losch = (101325/1.23e-23/273 * 1e-6)    ! (p0/k T0)
      y1 = fr**4/(9.38076e+18 - 1.08426e+09 * fr ** 2) * (wtot/Losch) 
      y1 = y1*layerthick/1000  !the abs coeff above is in km-1
    
c this is from Radiative Transfer in the Ocean, Stamnes&Thomas pg 73
c      a0=3.9729066
c      a1=4.6547659e-2
c      a2=4.5055995e-4
c      a3=2.3229848e-5
c      y0=0.0
c      y0=a0 + a1/(l*l) + a2/(l*l*l*l) + a3/(l*l*l*l*l*l)
c      y0=1e-28/(l*l*l*l)*y0              !!!cross section in cm2
c      y0 = y0*(500*100)                  !!!ty0pical lay0er ~ 500 m
c      y0=y0*wtot                         !!! cm3 * cm-3 = optical depth

c try this one, from same book
c      a0=1.000292        !!ref index of air
c      !change wavelength to m (1um = 1e-6 m) 
c      !change wtot to m-3 (by multiplying by 1e6)
c      y1=32*(kPi**3)*((a0-1)**2)/(3*l*l*l*l*1e-24*wtot*1e6)
c      y1=y1*250             !!typicl layer ~ 250m   

c try this one, from same book
c      y2=8*kpi/3*((2*kpi/l*1e-6)**4)*(a0*a0-1)/(a0*a0+2)*1e-9*(wtot*1e6)

c try this from Applied Optics
c      y3 = 8*(kPi**3)/3*(wtot*1e6)*((a0*a0-1)/(a0*a0+2)**2)*(l*1e-6)

c try this from Liou
c      !!!need wavelength in um to compute refractive index (eqn 3.70)
c      a0 = 1/(l*l)
c      a0 = 6452.8 + 2949810/(146-a0) + 25540/(41-a0)
c      a0 = a0*1.0e-8 + 1       !!this is the ref index

c      a1 = 0.035               !!this is the anisotropy factor
    
c      !!!now need wavelength in cm
c      l1 = l*1.0e-4
c      y3 = 8*(kPi**3)/3/(wtot)*((a0*a0-1)**2) *(6+3*a1)/(6-7*a1)/(l1**4)*25000

c      !!!need wavelength in cm to compute scattering cross section (eqn 3.72)
c      y3 = 8*(kPi**3)/3/(l1**4)
c      y3 = y3 * ((a0*a0-1)**2)
c      y3 = y3 * (6+3*a1)/(6-7*a1)   
c      Losch = (101325/1.23e-23/273 * 1e-6)    ! (p0/k T0)
c      y3 = y3/Losch/Losch       !!!!!!!ONLY WAY TO MAKE THINGS WORK    
c      y2 = y3   !!!!!!this is eqn 3.72 which is scattering PER molecule
      
c      y3 = y3*wtot*(layerthick*100)  !!each layer ~ 250 m thick ~ 25000 cm
c      print *,fr,wtot,y1,y3

      rayleigh = y1

      RETURN
      END


c************************************************************************
c this subroutine finds out the index of wavenumber points at which to do 
c the radiative transfer
      SUBROUTINE FindWavenumberPoints(iStep,iScatter,raCorrelatedK,
     $                                iaStep,iFFMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iStep               !tells how many points to step thru, or find
      INTEGER iScatter            !1,2 or 3 tells which type to use
      INTEGER iaStep(kMaxPts)     !gives the index of wavenumber points to use
      INTEGER iFFMax              !tells how many of the 10000 pts to use
      REAL raCorrelatedK(kMaxPts) !array of abs coeffs of layer closest to gnd

      INTEGER iF,iL

      DO iF=1,kMaxPts
        iaStep(iF)=0
        END DO

      IF ((iScatter .LT. 1) .OR. (iScatter .GT. 3)) THEN
        write (kStdErr,*) 'For DISORT, need to choose kScatter = 1,2 or 3'
        write (kStdErr,*) 'Please retry'
        CALL DoStop
        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (iScatter .EQ. 1) THEN      
        !!!!! use equal spacing of points; nothing special about chosen points
        !!!!! do radiative transfer on freq index pts 1,1+iStep,1+2*iStep,...
        !!!!!       1+N*iStep,10000
        iL=1
        iF=1
 10     CONTINUE
        iaStep(iL) = iF
        iL=iL+1
        iF=iF+iStep
        IF (iF .LE. kMaxPts) THEN
          GOTO 10
        ELSE
          GOTO 20
          END IF
 20     CONTINUE
        IF (iaStep(iL-1) .NE. kMaxPts) THEN
          !make sure last point radiances are computed for, is the 10000th pt
          iaStep(iL) = kMaxPts
          iFFMax = iL
        ELSE
          iL=iL-1
          iFFMax = iL
          END IF
        END IF            !IF (kScatter .EQ. 1) THEN      
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      IF (iScatter .EQ. 2) THEN      
        !!!!! do radiative transfer on freq index pts a,b,c,...z
        !!!!! where the points are chosen such that they have the smallest
        !!!!! absorption coeffs
        CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
        END IF            !IF (kScatter .EQ. 2) THEN      

      IF (iScatter .EQ. 3) THEN      
        !!!!! do radiative transfer on freq index pts a,b,c,...z
        !!!!! where the points are chosen such that they range from the
        !!!!! smallest to largest absorption coeffs
        CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
        END IF            !IF (kScatter .EQ. 3) THEN      

      RETURN
      END
c************************************************************************
c this subroutine sorts the k values and returns index values
c if iScatter = 2, sort the k values from smallest to largest
c                  return indices of smallest iStep values
c if iScatter = 3, sort the k values from smallest to largest
c                  return indices of smallest to largest k, in steps of iStep
c the outputs are iaStep and iFFMax
      SUBROUTINE FindIndices(iScatter,iStep,arr,iaStep,iFFMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iStep               !tells how many points to step thru, or find
      INTEGER iScatter            !1,2 or 3 tells which type to use
      INTEGER iaStep(kMaxPts)     !gives the index of wavenumber points to use
      INTEGER iFFMax              !tells how many of the 10000 pts to use
      REAL arr(kMaxPts)           !array of abs coeffs of layer closest to gnd

      INTEGER indx(kMaxPts),i,j,IF,IL
      REAL TempArr(kMaxPts)

      CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

c now we have sorted the k so that indx(i) contains the sorted keys info
c      DO i=1,kMaxPts
c        print *,i,indx(i),arr(indx(i))
c        END DO

      IF (iScatter .EQ. 2) THEN  
        !simply give smallest k values for all but last two points.
        !for last point, give largest k value
        !for second last point, give medium k value
        DO i=1,iStep-2
          iaStep(i) = indx(i)
          END DO
        j=kMaxPts-i
        IF (mod(j,2) .eq. 0) THEN
          j=j/2 + i
        ELSE
          j=(j+1)/2 + i
          END IF
c        iaStep(iStep-1) = indx(kMaxPts-100)
        iaStep(iStep-1) = indx(j)
        iaStep(iStep)   = indx(kMaxPts)
        iFFMax = iStep
        END IF

      IF (iScatter .EQ. 3) THEN  !give smallest to largest k values
        iL=1
        iF=1
 10     CONTINUE
        iaStep(iL) = indx(iF)
        iL=iL+1
        iF=iF+iStep
        IF (iF .LE. kMaxPts) THEN
          GOTO 10
        ELSE
          GOTO 20
          END IF
 20     CONTINUE
        IF (iaStep(iL-1) .NE. indx(kMaxPts)) THEN
          !make sure radiances are computed for largest k
          iaStep(iL) = indx(kMaxPts)
          iFFMax = iL
        ELSE
          iL=iL-1
          iFFMax = iL
          END IF
        END IF

c very important .. want to keep integrity of the Planck fcn!!!!!
c so make sure the original 1st and 10000th points are always used
       IF (iScatter .GE. 2) THEN

         !!!!!!make sure (iaStep(i) contains "1")
         j = -1
         iL = 1
 30      CONTINUE
         IF (iaStep(iL) .EQ. 1)  THEN
           j = +1      !!!!!yes, it does contain "1"
         ELSE
           iL = iL + 1 !!!!keep on hoping
           END IF
        IF ((iL .LE. iFFMax) .and. (j .lt. 0)) THEN
          GOTO 30
          END IF
        IF (j .lt. 0) THEN !!!!!!poooey, not found
          i = -1
          iL = 1
 40       CONTINUE
          IF (indx(iL) .EQ. 1) THEN
            i = +1       !!!!!yes, it is found
            iFFMax = iFFMax + 1
            iaStep(iFFMax) = indx(iL)
          ELSE 
            iL = iL + 1
            GOTO 40
            END IF
          END IF

         !!!!!!make sure (iaStep(i) contains "kMaxPts")
         j = -1
         iL = 1
 50      CONTINUE
         IF (iaStep(iL) .EQ. kMaxPts)  THEN
           j = +1      !!!!!yes, it does contain "kMaxPts"
         ELSE
           iL = iL + 1 !!!!keep on hoping
           END IF
        IF ((iL .LE. iFFMax) .and. (j .lt. 0)) THEN
          GOTO 50
          END IF
        IF (j .lt. 0) THEN !!!!!!poooey, not found
          i = -1
          iL = 1
 60       CONTINUE
          IF (indx(iL) .EQ. kMaxPts) THEN
            i = +1       !!!!!yes, it is found
            iFFMax = iFFMax + 1
            iaStep(iFFMax) = indx(iL)
          ELSE 
            iL = iL + 1
            GOTO 60
            END IF
          END IF

c reset iaStep(iFFMax)
          DO i=1,kMaxPts
            TempArr(i)=1.0e10
            END DO
          DO i=1,iFFMax
            TempArr(iaStep(i))=arr(iaStep(i))
            END DO
          CALL NumericalRecipesIndexer(indx,TempArr,kMaxPts)
          DO i=1,iFFMax
            iaStep(i)=indx(i)
            END DO
          END IF
        
      RETURN
      END

c************************************************************************
c this subroutine sorts the k values and returns index values
c if iScatter = 2, sort the k values from smallest to largest
c                  return indices of smallest iStep values
c if iScatter = 3, sort the k values from smallest to largest
c                  return indices of smallest to largest k, in steps of iStep
c the outputs are iaStep and iFFMax
      SUBROUTINE FindIndices_NotUsed(iScatter,iStep,arr,iaStep,iFFMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iStep               !tells how many points to step thru, or find
      INTEGER iScatter            !1,2 or 3 tells which type to use
      INTEGER iaStep(kMaxPts)     !gives the index of wavenumber points to use
      INTEGER iFFMax              !tells how many of the 10000 pts to use
      REAL arr(kMaxPts)           !array of abs coeffs of layer closest to gnd

      INTEGER indx(kMaxPts)
      REAL TempArr(kMaxPts),kmin,kmax,k0,dk

      INTEGER i,j,iF,iL

      CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

c now we have sorted the k so that indx(i) contains the sorted keys info
c      DO i=1,kMaxPts
c        print *,i,indx(i),arr(indx(i))
c        END DO

      IF (iScatter .EQ. 2) THEN  
        !simply give smallest k values for all but last two points.
        !for last point, give largest k value
        !for second last point, give medium k value
        DO i=1,iStep-2
          iaStep(i) = indx(i)
          END DO
        j=kMaxPts-i
        IF (mod(j,2) .eq. 0) THEN
          j=j/2 + i
        ELSE
          j=(j+1)/2 + i
          END IF
c        iaStep(iStep-1) = indx(kMaxPts-100)
        iaStep(iStep-1) = indx(j)
        iaStep(iStep)   = indx(kMaxPts)
        iFFMax = iStep
        END IF

      IF (iScatter .EQ. 3) THEN  
       !smallest to largest k values, equally stepped (kmax-kmin)/(iStep-1)
        kmin = arr(indx(1))
        kmax = arr(indx(kMaxPts))
        dk = (kmax-kmin)/(iStep-1)
        k0 = kmin

        k0 = k0 + dk
        iaStep(1) = indx(1)
        iF=2
        iL=2

 10     CONTINUE
        IF ((arr(indx(iL)) .LT. k0) .AND. (iL .LT. kMaxPts)) THEN
          iL = iL + 1
          GOTO 10
        ELSE
          iaStep(iF) = indx(iL)
          iF = iF + 1
          iL = iL + 1
          k0 = k0 + dk
          IF (k0 .LT. kmax) THEN
            GOTO 10
            END IF
          END IF
        IF (iaStep(iF-1) .NE. indx(kMaxPts)) THEN
          !make sure radiances are computed for largest k
          iaStep(iF) = indx(kMaxPts)
          iFFMax = iF
        ELSE
          iF=iF-1
          iFFMax = iF
          END IF
        END IF

c very important .. want to keep integrity of the Planck fcn!!!!!
c so make sure the original 1st and 10000th points are always used
       IF (iScatter .GE. 2) THEN

         !!!!!!make sure (iaStep(i) contains "1")
         j = -1
         iL = 1
 30      CONTINUE
         IF (iaStep(iL) .EQ. 1)  THEN
           j = +1      !!!!!yes, it does contain "1"
         ELSE
           iL = iL + 1 !!!!keep on hoping
           END IF
        IF (iL .LE. iFFMax) THEN
          GOTO 30
          END IF
        IF (j .lt. 0) THEN !!!!!!poooey, not found
          i = -1
          iL = 1
 40       CONTINUE
          IF (indx(iL) .EQ. 1) THEN
            i = +1       !!!!!yes, it is found
            iFFMax = iFFMax + 1
            iaStep(iFFMax) = indx(iL)
          ELSE 
            iL = iL + 1
            GOTO 40
            END IF
          END IF

         !!!!!!make sure (iaStep(i) contains "kMaxPts")
         j = -1
         iL = 1
 50      CONTINUE
         IF (iaStep(iL) .EQ. kMaxPts)  THEN
           j = +1      !!!!!yes, it does contain "kMaxPts"
         ELSE
           iL = iL + 1 !!!!keep on hoping
           END IF
        IF (iL .LE. iFFMax) THEN
          GOTO 50
          END IF
        IF (j .lt. 0) THEN !!!!!!poooey, not found
          i = -1
          iL = 1
 60       CONTINUE
          IF (indx(iL) .EQ. kMaxPts) THEN
            i = +1       !!!!!yes, it is found
            iFFMax = iFFMax + 1
            iaStep(iFFMax) = indx(iL)
          ELSE 
            iL = iL + 1
            GOTO 60
            END IF
          END IF

c reset iaStep(iFFMax)
          DO i=1,kMaxPts
            TempArr(i)=1.0e10
            END DO
          DO i=1,iFFMax
            TempArr(iaStep(i))=arr(iaStep(i))
            END DO
          CALL NumericalRecipesIndexer(indx,TempArr,kMaxPts)
          DO i=1,iFFMax
            iaStep(i)=indx(i)
            END DO
          END IF
        
      RETURN
      END
c************************************************************************
      SUBROUTINE INTERP_PLANCK_0(XA,YA,CA,N,raFO,raNS,raOut)

c this does a very smart interpolation ahem ahem
c        CALL INTERP_PLANCK_0(
c     $          raFreqStep,raIntenStep,raNoScatterStep,iStepPts,
c     $          raWaves,raNoScatterAll,raInten)

      IMPLICIT NONE

       include '../INCLUDE/kcarta.param'

C linear interpolation
C      -----------------------------------------------------------------
C      XA  : I  : REAL arr : freq array   array(N) in increasing order
C      YA  : I  : REAL arr : intensity y  array(N) from scattering
C      CA  : I  : REAL arr : intensity y  array(N) from clear sky
C      N   : I  : INT      : number of points fom DISORT, in array
c
C      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
C      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
C      raOut : O  : REAL arr : output radiance             array (kMaxPts)
C
       REAL XA(*),YA(*),CA(*)
       REAL raFO(*),raNS(*),raOut(*)
       INTEGER N

C      Local Variables
       INTEGER K,KLO,KHI,iI
       REAL A,B,H,radtot,ttorad
       REAL yn,y0,X,y,w,dx

      DO iI=1,kMaxPts
        X=raFO(iI)
        IF (X .gt. xa(n)) X=xa(n)-1.0e-7
        IF (X .lt. xa(1)) X=xa(1)+1.0e-7

        !Determine between which pair of points X falls (bisect loop)
        KLO=1
        KHI=N

 20     IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
            KHI=K
          ELSE
            KLO=K
            ENDIF
          GOTO 20
          ENDIF
C
        dx = XA(KHI) - x
        H  = XA(KHI) - XA(KLO)

        IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI),X,N
 1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,',X = ',E12.5,', numpts = ',I4,'. Quitting.')
          ENDIF

c METHOD 1
c      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
c      Y=(YA(KHI) + dh) - h*b 

c METHOD 2
c interpolate radiances
C      YA  : I  : REAL arr : intensity y array(N) from scattering
c      CA  : I  : REAL arr : intensity y array(N) from clear sky
c        !!!interpolate in radiance space
c        yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
c        y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
c        B=(yn-y0)/(XA(KHI) - XA(KLO)) !!slope
c        Y=yn - dx*b   
c        Y=Y+raNS(iI) 

        !!!interpolate in temperature space
        !!! get out stuff from clear sky
        yn = radtot(xa(khi),ya(KHI))-radtot(xa(khi),CA(KHI))      
        y0 = radtot(xa(klo),ya(Klo))-radtot(xa(klo),CA(Klo))      
        B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
        Y = yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
        Y = Y + radtot(x,raNS(iI)) !!add on clear sky stuff
        Y = ttorad(x,y)

c interpolate BTs
c      B=(radtot(WA(KHI),YA(KHI))-radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
c      Y=radtot(WA(KHI),YA(KHI))-dx*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
c      Y=ttorad(W,Y)+0*dh              !!change back to rad, add on  dr/dvo dv

c METHOD 3
c      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv

c      print *,w,x,xa(khi),y,ya(khi),(y-ya(khi))/ya(khi)
c      read *,khi
        raOut(iI)=Y
        END DO

      RETURN
      END

c************************************************************************
c this does a linear interpolation, but tries to weight the y axis so that
c information relating to the Planck depenadance on wavenumber, is included
       SUBROUTINE INTERP_PLANCK_1(WA,DWA,XA,YA,N,W,X,Y)

      IMPLICIT NONE

       include '../INCLUDE/kcarta.param'

C linear interpolation
C      -----------------------------------------------------------------
C      WA  : I  : REAL arr : wavenumber w array(N)
C     DWA  : I  : REAL arr : dy/dw so we can weight the interpolations
C      XA  : I  : REAL arr : abs coeff x array(N) in increasing order
C      YA  : I  : REAL arr : intensity y array(N)
C      N   : I  : INT      : number of points in arrays
C      W   : I  : REAL     : wavnumber at which spline is evaluated
C      X   : I  : REAL     : x point at which to evaluate spline
C      Y   : O  : REAL     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),WA(*),DWA(*),X,Y,W
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL A,B,H,dh,dl,radtot,ttorad

c can try METHOD 1
c using linear interp (yh-yh)/(xh-xh) = (yh-y)/(xh-x) ==>
c y = yh - (xh-x)(yh-yl)/(xh-xl)
c modify this to 
c y = YH - (xh-x)(YH-YL)/(xh-xl)
c where dh,dl take into account the variation of planck wrt wavenumber
c rad(vo+dv) = rad(vo) + d(rad)/dv * dv
c hence YH = yh + dr/dw dh and YL = yl + dr/dw dl

c can try METHOD 2
c r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk + dr/dvo dv
C
C      -----------------------------------------------------------------
C
C     Determine between which pair of points X falls (bisect loop)
      KLO=1
      KHI=N

 20   IF ( (KHI - KLO) .GT. 1) THEN
        K=(KHI + KLO)/2
        IF (XA(K) .GT. X) THEN
          KHI=K
        ELSE
          KLO=K
          ENDIF
        GOTO 20
        ENDIF
C
      H=XA(KHI) - XA(KLO)
      dh = W - WA(khi)       !!!wavenumber diff
      dl = WA(klo) - W       !!!wavenumber diff       
      dh = dh * DWA(khi)     !!!first derivative
      dl = dl * DWA(klo)     !!!first derivative

      IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
        ENDIF

c METHOD 1
c      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
c      Y=(YA(KHI) + dh) - h*b 

c METHOD 2
c interpolate radiances
c      B=(YA(KHI) - YA(KLO))/(XA(KHI) - XA(KLO)) !!slope
c      Y=YA(KHI) - h*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
c      Y=Y+dh            !!add on                            + dr/dvo dv
c interpolate BTs
      B=(radtot(WA(KHI),YA(KHI)) - radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
      Y=radtot(WA(KHI),YA(KHI)) - h*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
      Y=ttorad(W,Y)+0*dh                !!change back to rad, add on  dr/dvo dv

c METHOD 3
c      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv

c      print *,w,x,xa(khi),y,ya(khi),(y-ya(khi))/ya(khi)
c      read *,khi

      RETURN
      END

c************************************************************************
      SUBROUTINE INTERP_PLANCK_2(WA,DWA,XA,YA,CA,N,
     $          raFO,raNS,raCK,raOut)

      IMPLICIT NONE

       include '../INCLUDE/kcarta.param'

C linear interpolation
C      -----------------------------------------------------------------
C      WA  : I  : REAL arr : wavenumber w array(N)
C     DWA  : I  : REAL arr : dy/dw        array(N) to weight the interpolations
C      XA  : I  : REAL arr : corr coeff x array(N) in increasing order
C      YA  : I  : REAL arr : intensity y  array(N) from scattering
C      CA  : I  : REAL arr : intensity y  array(N) from clear sky
C      N   : I  : INT      : number of points fom DISORT, in array
c
C      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
C      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
C      raCK  : I  : REAL arr : entire correlated k         array (kMaxPts)
C      raOut : O  : REAL arr : output radiance             array (kMaxPts)
C
       REAL XA(*),YA(*),CA(*),WA(*),DWA(*)
       REAL raFO(*),raNS(*),raCK(*),raOut(*)
       INTEGER N

C      Local Variables
       INTEGER K,KLO,KHI,iI
       REAL A,B,H,dh,dl,radtot,ttorad
       REAL yn,y0,X,y,w,dx

c can try METHOD 1
c using linear interp (yh-yh)/(xh-xh) = (yh-y)/(xh-x) ==>
c y = yh - (xh-x)(yh-yl)/(xh-xl)
c modify this to 
c y = YH - (xh-x)(YH-YL)/(xh-xl)
c where dh,dl take into account the variation of planck wrt wavenumber
c rad(vo+dv) = rad(vo) + d(rad)/dv * dv
c hence YH = yh + dr/dw dh and YL = yl + dr/dw dl

c can try METHOD 2
c r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk + dr/dvo dv
C

      DO iI=1,kMaxPts
        X=raCK(iI)
        W=raFO(iI)

        !Determine between which pair of points X falls (bisect loop)
        KLO=1
        KHI=N

 20     IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
            KHI=K
          ELSE
            KLO=K
            ENDIF
          GOTO 20
          ENDIF
C
        dx = XA(KHI) - x
        H  = XA(KHI) - XA(KLO)
        dh = W - WA(khi)       !!!wavenumber diff
        dl = WA(klo) - W       !!!wavenumber diff       
c        print *,dl,dh,XA(KLO),X,XA(KHI)
        dh = dh * DWA(khi)     !!!first derivative
        dl = dl * DWA(klo)     !!!first derivative

        IF (H .LE. 0.0) THEN
          WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
          ENDIF

c METHOD 1
c      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
c      Y=(YA(KHI) + dh) - h*b 

c METHOD 2
c interpolate radiances
C      YA  : I  : REAL arr : intensity y array(N) from scattering
C      CA  : I  : REAL arr : intensity y array(N) from clear sky
c        !!!interpolate in radiance space
c        yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
c        y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
c        B=(yn-y0)/(XA(KHI) - XA(KLO)) !!slope
c        Y=yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
c        Y=Y+raNS(iI) !!add on clear sky stuff
c        !!add on dr/dvo dv 
c        !!CANNOT DO THIS as we are modelling differences, not actual numbers
c        !!Y=Y+dh

        !!!interpolate in temperature space
        !!! get out stuff from clear sky
        yn = radtot(wa(khi),ya(KHI))-radtot(wa(khi),CA(KHI))      
        y0 = radtot(wa(klo),ya(Klo))-radtot(wa(klo),CA(Klo))      
        B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
        Y = yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
        Y = Y + radtot(w,raNS(iI)) !!add on clear sky stuff
        Y = ttorad(w,y)

c interpolate BTs
c      B=(radtot(WA(KHI),YA(KHI))-radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
c      Y=radtot(WA(KHI),YA(KHI))-dx*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
c      Y=ttorad(W,Y)+0*dh              !!change back to rad, add on  dr/dvo dv

c METHOD 3
c      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv

c      print *,w,x,xa(khi),y,ya(khi),(y-ya(khi))/ya(khi)
c      read *,khi
        raOut(iI)=Y
        END DO

      RETURN
      END

c************************************************************************
c same as PLANCK2 except we do interpolations in K space
      SUBROUTINE INTERP_PLANCK_3(WA,XA,YA,CA,N,raFO,raNS,raCK,raOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

C linear interpolation
C      -----------------------------------------------------------------
C      WA  : I  : REAL arr : wavenumber w array(N)
C      XA  : I  : REAL arr : corr coeff x array(N) in increasing order
C      YA  : I  : REAL arr : intensity y  array(N) from scattering
C      CA  : I  : REAL arr : intensity y  array(N) from clear sky
C      N   : I  : INT      : number of points fom DISORT, in array
c
C      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
C      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
C      raCK  : I  : REAL arr : entire correlated k         array (kMaxPts)
C      raOut : O  : REAL arr : output radiance             array (kMaxPts)
C
      REAL XA(*),YA(*),CA(*),WA(*)
      REAL raFO(*),raNS(*),raCK(*),raOut(*)
      INTEGER N

C     Local Variables
      INTEGER K,KLO,KHI,iI,indx(kMaxPts),iDok,iDoF
      REAL A,B,H,radtot,ttorad
      REAL yn,y0,X,Y1,Y2,YTotal,dx
c these have WA,XA,YA,CA sorted according to wavenumber
      REAL WA_sort(kMaxPts),XA_sort(kMaxPts)
      REAL YA_sort(kMaxPts),CA_sort(kMaxPts)

 1010          FORMAT('ERROR! linear SPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')

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
        IF (iDoK .GT. 0) THEN
          !!!!! -------------  first do dT/dk dk !!!!!!!!!!!!!
          X=raCK(iI)
 
          !Determine between which pair of points X falls (bisect loop)
          KLO=1
          KHI=N
  
 20       IF ( (KHI - KLO) .GT. 1) THEN
            K=(KHI + KLO)/2
            IF (XA(K) .GT. X) THEN
              KHI=K
            ELSE
              KLO=K
              ENDIF
            GOTO 20
            ENDIF

          dx = XA(KHI) - x
          H  = XA(KHI) - XA(KLO)

          IF (H .LE. 0.0) THEN
            WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
            ENDIF

          !!!interpolate in temperature space
          !!! get out stuff from clear sky
          yn = radtot(wa(khi),ya(KHI))-radtot(wa(khi),CA(KHI))      
          y0 = radtot(wa(klo),ya(Klo))-radtot(wa(klo),CA(Klo))      
          B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
          Y1 = yn - dx*b   !!this is T(vo+dv,ko+dk) = T(vo,ko) + dT/dko dk
          END IF

        Y2=0.0
        IF (iDoF .GT. 0) THEN
          !!!!! -------------  then do dT/dv dv !!!!!!!!!!!!!
          X=raFO(iI)
  
          !Determine between which pair of points X falls (bisect loop)
          KLO=1
          KHI=N

 30       IF ( (KHI - KLO) .GT. 1) THEN
            K=(KHI + KLO)/2
            IF (WA_Sort(K) .GT. X) THEN
              KHI=K
            ELSE
              KLO=K
              ENDIF
            GOTO 30
            ENDIF
C
          dx = WA_Sort(KHI) - x
          H  = WA_Sort(KHI) - WA_Sort(KLO)

          IF (H .LE. 0.0) THEN
            WRITE(kStdWarn,1010) KLO,KHI,WA_Sort(KLO),WA_Sort(KHI)
            ENDIF

          !!!interpolate in temperature space
          !!! get out stuff from clear sky
          yn = radtot(wa_sort(khi),ya_sort(KHI)) - 
     $         radtot(wa_sort(khi),ca_sort(KHI))
          y0 = radtot(wa_sort(klo),ya_sort(Klo)) -
     $         radtot(wa_sort(klo),ca_sort(Klo))
          B = (yn-y0)/(wa_sort(KHI) - wa_sort(KLO)) !!slope
          Y2 = yn - dx*b   !!this is T(vo+dv,ko+dk) = T(vo,ko) + dT/dvo dv
          END IF

c remember      iDoK = -1       !!!!!do not interpolate in K, so Y1 == 0.0
c remember      iDoF = +1       !!!!!do interpolate in F    , so Y2 <> 0.0

        YTotal = Y1 + Y2        !!this is the total dT/dk dk + dT/dv dv
        YTotal = YTotal + radtot(raFO(iI),raNS(iI)) !!add on clear sky stuff
        !!!!!now change from BT to radiance
        YTotal = ttorad(raFO(iI),YTotal)
        raOut(iI)=YTotal

        END DO

      RETURN
      END

c************************************************************************
c this subrtouine computes drad/dwavenumber for a Planck blackbody radiance
       SUBROUTINE drad_dv(WA,DWA,YA,N)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

C linear interpolation
C      -----------------------------------------------------------------
C      WA  : I  : REAL arr : wavenumber w array(N)
C     DWA  : O  : REAL arr : dy/dw so we can weight the interpolations
C      YA  : I  : REAL arr : intensity y array(N)
C      N   : I  : INT      : number of points in arrays
C      -----------------------------------------------------------------
C
C     Parameters
      REAL YA(*),WA(*),DWA(*)
      INTEGER N
      INTEGER i
      REAL raT(kMaxPts),r1,r2,v,u,dv,du

c we know r(fo+df) = r(fo) + df dr/df
c Let us assume r(fo) = r(fo,T) ===> first we have to find equivalent temp T
c having found this, we then simply diff dr/df where r = radiance at temp T
c so given radiance r(f), do rad2bt(f,r) --> T
c                         then differentiate d(rPlanck(f,T))/df

       r1=kPlanck1
       r2=kPlanck2

c first convert the intensities to equivalent temps T
c bt = c2 * fr / log(1 + c1 * fr^3 / rd)
       DO i = 1,N
         raT(i) = r2*wa(i)/(log(1+r1*(wa(i)**3)/ya(i)))  
         END DO

c now do dr/df where r=ttorad = planck radiance at temp T, frequency fr
c rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
       DO i = 1,N
         u  = r1*(wa(i)**3)
         du = r1*3*(wa(i)**2) 
         v  = exp(r2*wa(i)/raT(i))-1
         dv = r2/raT(i) * exp(r2*wa(i)/raT(i))
         !!!! now do the derivative
         dwa(i) = (v*du - u*dv)/(v*v)
         END DO

       RETURN
       END

c************************************************************************
c this doe a cumulative K distribution fcn, and then for the first 50%
c of the lowest k values, just does an average of intensity
      SUBROUTINE CumulativeK(raKall,raF,raK,raI,iN)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL raKAll(kMaxPts)        !this is all 10000 pts k(lowest layer)
      REAL raF(kMaxPts)           !these are freqs of pts where DISORT used
      REAL raK(kMaxPts)           !these are ks of pts where DISORT used
      REAL raI(kMaxPts)           !intensities of pts where DISORT used
      INTEGER iN                  !how many points DISORT went over

c note that we will do a sort on k, and then accordingly, average some of the
c points in raK and raI so as to make things smoother

      INTEGER iCnt50,iCnt75,iInCnt50,iInCnt75
      REAL rCnt50,rCnt75,rF,rI,rK,radtot,ttorad
      INTEGER indx(kMaxPts),iI,iJ,iCnt
      REAL NormalizedK(kMaxPts),kmin,kmax,dk
      REAL TempF(kMaxPts),TempI(kMaxPts),TempK(kMaxPts)
      REAL raCumulativeK(kMaxPts),raCumulativeDistr(kMaxPts)

      CALL NumericalRecipesIndexer(indx,raKAll,kMaxPts)
      
      !!set up the normalised vector knorm = ksorted/max(k)
      DO iI=1,kMaxPts
        NormalizedK(iI)=raKall(indx(iI))/raKall(indx(kMaxPts))
        END DO

      kmin=0.0
      kmax=NormalizedK(kMaxPts)   !!!should be 1.0

      iJ = 500
      dk = 1.0/iJ

c compute the cumulative distr function
      iCnt=1
      DO iI=1,iJ
        raCumulativeK(iI)=kmin + iI*dk
 10     CONTINUE
        IF (NormalizedK(iCnt) .LE. raCumulativeK(iI)) THEN
          iCnt = iCnt + 1
          IF (iCnt .LE. kMaxPts) THEN
            GOTO 10
            END IF
          END IF
        iCnt = min(iCnt,kMaxPts)
        raCumulativeDistr(iI)=iCnt*1.0/kMaxPts
        END DO

c now see which index required before 50% of the cumulative function is used
      iCnt50 = 1
      iCnt75 = 1
 30   CONTINUE
      IF (raCumulativeDistr(iCnt50) .LT. 0.5) THEN
        iCnt50 = iCnt50 + 1
        GOTO 30
        END IF
 40   CONTINUE
      IF (raCumulativeDistr(iCnt75) .LT. 0.75) THEN
        iCnt75 = iCnt75 + 1
        GOTO 40
        END IF

      !find values of raK corresponding to the normalised k=0.5,0.75 values
      !remember raK is where DISORT is called, and so these should be 
      !already sorted from smallest to largest
      rCnt50 = raCumulativeK(iCnt50)*raKall(indx(kMaxPts))
      rCnt75 = raCumulativeK(iCnt75)*raKall(indx(kMaxPts))
      iInCnt50 = 1
      iInCnt75 = 1
 50   CONTINUE
      IF (raK(iInCnt50) .LT. rCnt50) THEN
        iInCnt50 = iInCnt50 + 1
        GOTO 50
        END IF
 60   CONTINUE
      IF (raK(iInCnt75) .LT. rCnt75) THEN
        iInCnt75 = iInCnt75 + 1
        GOTO 60
        END IF
c      print *,'normalised cml upto 50% = ',iCnt50,raCumulativeK(iCnt50),rCnt50
c      print *,'normalised cml upto 75% = ',iCnt75,raCumulativeK(iCnt75),rCnt75
c      print *,'input cml upto 50% = ',iInCnt50,raK(iInCnt50)
c      print *,'input cml upto 75% = ',iInCnt75,raK(iInCnt75)

      !since there is so much variation in the lowest radainces, and they 
      !account for so much of raK, let us avg the first few values and 
      !shove the rest around
      IF (raCumulativeK(iCnt50) .LE. 0.1) THEN
        !first initialise the temp arrays
        DO iI=1,iN
          TempF(iI)=raF(iI)
          TempK(iI)=raK(iI)
          TempI(iI)=raI(iI)
          END DO
        !50% of the (k/kmax) are smaller in magnitude that 0.1
        !this will drive the interpolation wrt raKStep nuts!!!!
        !so do an average; 
        !make sure at least FIVE points remain for the interpolations
        IF ((iN - iInCnt50) .LE. 5) THEN
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
        raF(1)=raF(1)
        raK(1)=raK(1)
        raI(1)=raI(1)
        !you have averaged the next few to get something
        raF(2)=rF
        raK(2)=rK
        raI(2)=rI
        !now fill in the larger "k" values
        DO iI=3,iN - iInCnt50 + 2
          raF(iI)=TempF(iInCnt50 + iI - 2)
          rak(iI)=Tempk(iInCnt50 + iI - 2)
          raI(iI)=TempI(iInCnt50 + iI - 2)
          END DO
        IN = iN - iInCnt50 + 2
        END IF

      RETURN
      END
c************************************************************************
