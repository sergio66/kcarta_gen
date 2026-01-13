c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the rtspec routines  **********************
c ***********   this is the LATEST rtspec version  **********************

C     Author: Frank Evans  University of Colorado   May 1999 
C     Co-authors: 
C       Aaron Evans    (multilayer code, spectral segments, double sidebands, 
C                       upward looking geometry, output convolution) 
C       Merritt Deeter (single scattering layer routines) 

c modifications by Sergio De Souza-Machado
c 1) all subroutines now call include 'scatter.param'
c 2) gasRT1 now has surface emissivity passed in, so it can be used
c 3) nlev ---> nlev - 1 whenever we are using optical depths
c    ie there are nlev lavels and so there are nlev-1 layers.
c 4) background thermal has been included!!!!!!!!!!!!!!!!!!!!!!
c    when called, radobs=0.0 if there is no backgnd thermal
c                 radobs=downward thermal radiance at posn of instrument
c                        simply computed using acos(3/5)
               
c************************************************************************
c************************************************************************
      SUBROUTINE COMPUTE_RADIATIVE_TRANSFER (RTMODEL,
     $               MUOBS, IOBS, WAVENO,
     $               NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!     $               MAXTAB, MAXGRID, 
     $               NSCATTAB, MUINC,
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,
     $               TABEXTINCT, TABSSALB, TABASYM,
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN,
     $               RADOBS,
     $               ibdry, TOA_to_instr)  !these are 2 new parameters by SSM

      include 'scatter.param'

C       Calculates the radiance at the desired level (IOBS) and angle (MUOBS)
C     for one wavenumber (WAVENO) for the atmosphere and cloud specified.
C     The temperature and molecular absorption profiles are input
C     (TEMP, TAUGAS) for the NLEV levels.  The multilayer cloud is
C     between levels ICLDTOP and ICLDBOT. Each of the NCLDLAY cloud layers
C     has an ice water path (IWP), median particle diameter (DME), and
C     pointer to the scattering table (ISCATTAB).  The scattering properties
C     are input as tables. Accept mu > 0 (down looking) or mu < 0 (up looking).

C     The observed radiance (RADOBS) is output in W/(m^2 ster cm^-1) units.
C     RTMODEL : S,E,H for single scatter, Eddington, and Hybrid respectively
C     MUOBS : Cosine of observation angle (negative for looking up)
C     IOBS : The level number of observation height
C     WAVENO : Wave number (1/cm) for this monochromatic RT calculation.
C     NLEV : Number of levels in atmosphere
C       Note: Levels/layers numbering start at top of atmopshere
C     TEMP : Temperature array for levels in atmosphere
C     TAUGAS : Optical depth array for gases in layer of atmosphere
C     SFCTEMP : Temperature of surface
C     SFCEMIS : Emissivity of surface
C     NCLDLAY : Number of cloud layers, ie. # of atm layer occupied by cloud
C     ICLDTOP : The level number of cloud top
C     ICLDBOT : The level number of cloud bottom
C     IWP : Ice Water Path (g/m^2) array for cloud layers
C     DME : Median Particle Diameter (um) array for cloud layers
C     ISCATTAB : Scattering table file number array for cloud layers
C       Note: ISCATAB allows multiple scattering table files
C     MAXTAB : Maximum number of elements in scattering table
C     MAXGRID : Maximum number of points on scattring table axis
C     NSCATTAB : Number of scattering tables/files
C     MUINC : Cosines of two incident angles of single scattering solution
C     NMUOBS : Number of viewing cosines in scattering table
C     MUTAB : Viewing cosines making up grid of scattering table
C     NDME : Number of Median Particle Diameters in scattering table
C     DMETAB : Median Particle Diameters making up grid of scattering table
C     NWAVETAB : Number of wavenumbers in scattering table
C     WAVETAB : Wavenumbers making up grid of scattering table
C     TABEXTINCT : Extinction (1/km) matrix, dimensions are wavenumber/Dme
C       Note : Scattering table calculated for IWC = 1 g/m^2, thus 
C              Extinction = TABEXTINCT*IWP/(layer thickness) 
C     TABSSALB : Single Scattering albedo matrix, dimensions are Wavenumber/Dme
C     TABASYM : Asymmetry Parameter matrix, dimensions are Wavenumber/Dme
C     TABPHI : "Phi" function matrix, dimensions are Wavenumber/Dme/Mu
C              "Phi" represents the amount of incident radiation scattered.
C     TABPHI1UP : "Phi" up applied to forward scattering radiation, for MU1
C     TABPHI2UP : "Phi" up applied to forward scattering radiation, for MU2
C     TABPHI1DN : "Phi" down applied to back scattering radiation, for MU1
C     TABPHI2DN : "Phi" down applied to back scattering radiation, for MU2
C     RADOBS : Output of subroutine, Upwelling radiation at observation height

      INTEGER NLEV, IOBS
      INTEGER NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(NCLDLAY)
      REAL    MUOBS, WAVENO, SFCTEMP, SFCEMIS
      REAL    TEMP(NLEV), TAUGAS(NLEV-1)
      REAL    IWP(NCLDLAY), DME(NCLDLAY), IWPTOT
      REAL    RADOBS
      CHARACTER*1 RTMODEL
      INTEGER  NSCATTAB   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAXTAB, MAXGRID
      INTEGER  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
      REAL     MUTAB(MAXGRID,NSCATTAB)
      REAL     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
      REAL     TABASYM(MAXTAB,NSCATTAB)
      REAL     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
      REAL     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
      INTEGER   IBDRY                          !!!for diffusive approx
      REAL      TOA_to_instr                   !!!sum(k) from TOA to instr

C      Local variables
c      INTEGER  MAXNZ
c      PARAMETER (MAXNZ=81)
      INTEGER I, L, N
      REAL    TAUTOT(MAXNZ), COALB(MAXNZ), EXTINCT
      REAL    TAUC(MAXNZ), TAUCG(MAXNZ), SSALB(MAXNZ), ASYM(MAXNZ)
      REAL    PHI1UP(MAXNZ), PHI1DN(MAXNZ)
      REAL    PHI2UP(MAXNZ), PHI2DN(MAXNZ)
      REAL    RADBNDRYUP(2), RADBNDRYDN(2)
      REAL    RADBNDRYUP1(MAXNZ), RADBNDRYDN1(MAXNZ)
      REAL    RADBNDRYUP2(MAXNZ), RADBNDRYDN2(MAXNZ)
      REAL    RAD0UPOBS(MAXNZ), RAD0DNOBS(MAXNZ)
      REAL    FLUXES(3,MAXNZ)
      REAL    FLUXTOP, FLUXBOT, TTOP, TBOT, RADBOTTOM, RADTOP
      REAL    FLUXUP, FLUXDN, FLUXUPSEDD, FLUXDNSEDD

      REAL rDummy1,rDummy2

C     Initialize the optical depths at each level
      DO N = 1, NLEV-1
        TAUTOT(N) = TAUGAS(N)
        COALB(N) = 1.0
        ENDDO

      IWPTOT = 0
      DO N = 1,NCLDLAY
        IWPTOT = IWPTOT + IWP(N)
        ENDDO

C     If clear sky only do GASRT1 and GASRT2
cccccccc this was orig rtspec code       IF (IWPTOT.EQ.0) THEN
      IF ((IWP(1) .LE. 0)) THEN
C       Integrate radiance from surface to observer
C       First go from surface to cloud bottom 
        IF (MUOBS .GT. 0) THEN
          !if looking down, do the nocloud downlook rad trans
          CALL GASRT1_nocloud (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS,
     $        NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT, 
     $        NCLDLAY, RAD0UPOBS, RAD0DNOBS,
     $        ibdry, TOA_to_instr)
        ELSE 
          !this was original Evans code
          CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS, 
     $        NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT,  
     $                   NCLDLAY, RAD0UPOBS, RAD0DNOBS) 
          END IF
        IF(MUOBS.GT.0) THEN 
C         if looking down then go from cloud top to observer
          RADOBS = RAD0UPOBS(1)
          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
     $                ICLDTOP+1, IOBS, MUOBS, RADOBS)
        ELSEIF(MUOBS.LT.0) THEN
C         Then go from cloud bottom to observer
          RADOBS = RAD0DNOBS(NCLDLAY)
          rDummy1=RADOBS
          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
     $                ICLDBOT-1, IOBS, ABS(MUOBS), RADOBS)
          rDummy2=RADOBS
          ENDIF

      ELSE
C       If cloud present do Single-scatter,Eddington or hybrid routine

C       Get the optical properties for the cloud layers
        DO N = ICLDTOP, ICLDBOT - 1
          L = N-ICLDTOP+1
C             Interpolate to get values of extinction, s.s. albedo, and
C             phi function values at given obs. mu, waveno, and particle size.
C             Note: we don't need phi functions for Eddington only solution$   
          I = ISCATTAB(L)

          IF (RTMODEL .EQ. 'E') THEN
            CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L),    
     $                EXTINCT, SSALB(L), ASYM(L),
     $                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I),
     $                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
          ELSE
            CALL INTERP_SCAT_TABLE3 (ABS(MUOBS), WAVENO, DME(L),    
     $                EXTINCT, SSALB(L), ASYM(L),
     $                PHI1UP(L), PHI1DN(L), PHI2UP(L), PHI2DN(L),
     $                NMUOBS(I), MUTAB(1,I), 
     $                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I),
     $                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
     $                TABPHI1UP(1,I), TABPHI1DN(1,I), 
     $                TABPHI2UP(1,I), TABPHI2DN(1,I))   
          ENDIF

C         Compute the optical depth of cloud layer, including gas
          TAUC(L) = IWP(L)*EXTINCT/1000.
          TAUCG(L) = TAUGAS(N) + TAUC(L)
          IF (TAUCG(L) .GT. 1.0E-5) THEN
            SSALB(L) = SSALB(L)*TAUC(L)/TAUCG(L)
          ELSE
            SSALB(L) = 0.0
            ENDIF
          TAUTOT(N) = TAUCG(L)
          COALB(N) = 1.0 - SSALB(L)
          ENDDO

C       Calculate the gas/cloud emission only radiative transfer 
C       for incident radiation on each layer
        CALL GASRT1 (MUINC(1), WAVENO, SFCTEMP, SFCEMIS,
     $               NLEV, TEMP, TAUTOT, COALB,
     $               ICLDTOP, ICLDBOT, NCLDLAY, 
     $               RADBNDRYUP1, RADBNDRYDN1)
        CALL GASRT1 (MUINC(2), WAVENO, SFCTEMP, SFCEMIS, 
     $               NLEV, TEMP, TAUTOT, COALB,
     $               ICLDTOP, ICLDBOT, NCLDLAY, 
     $               RADBNDRYUP2, RADBNDRYDN2)
        CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS,
     $               NLEV, TEMP, TAUTOT, COALB,
     $               ICLDTOP, ICLDBOT, NCLDLAY,  RAD0UPOBS, RAD0DNOBS)

C       Calculate the fluxes for the multilayer boundary conditions
        FLUXTOP = MUINC(1)*.5*RADBNDRYDN1(1) 
     $          + MUINC(2)*.5*RADBNDRYDN2(1)
        FLUXBOT = MUINC(1)*.5*RADBNDRYUP1(NCLDLAY)
     $          + MUINC(2)*.5*RADBNDRYUP2(NCLDLAY)

C       Don't compute multilayer fluxes for a single cloud layer
        IF (NCLDLAY .LE. 1) THEN
          FLUXES(2,1) = FLUXTOP
          FLUXES(1,2) = FLUXBOT
        ELSE
C         For multilayer cloud calculate the multilayer Eddington fluxes
	  CALL EDDRTF (NCLDLAY, TAUCG, SSALB, ASYM,
     $                 TEMP(ICLDTOP), WAVENO, 
     $                 FLUXTOP, FLUXBOT,  FLUXES)
          ENDIF


C       If looking down must loop over cloud layers from bottom to top 
        IF(MUOBS.GT.0) THEN
C         Introduce scalar variable for radiance at bottom of cloud
          RADBOTTOM = RAD0UPOBS(NCLDLAY)
          
C         Loop over the layers of cloud, from bottom to top
          DO L = NCLDLAY,1, -1
            TTOP = TEMP(L-1+ICLDTOP)
            TBOT = TEMP(L+ICLDTOP)
            FLUXUP = FLUXES(1,L+1)
            FLUXDN = FLUXES(2,L)
  
C           Do Single-scattering, Eddington, or Hybrid for each layer 
            IF (RTMODEL .EQ. 'S') THEN
              RADBNDRYUP(1) = RADBNDRYUP1(L)
              RADBNDRYUP(2) = RADBNDRYUP2(L)
              RADBNDRYDN(1) = RADBNDRYDN1(L)
              RADBNDRYDN(2) = RADBNDRYDN2(L)
              CALL SSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC,
     $                    RADBNDRYUP, RADBNDRYDN, RADBOTTOM,
     $                    TAUCG(L), SSALB(L), ASYM(L), 
     $                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L),
     $                    RADOBS)
            ELSE IF (RTMODEL .EQ. 'E') THEN
              CALL EDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, 
     $                      FLUXUP, FLUXDN, RADBOTTOM, 
     $                      TAUCG(L), SSALB(L), ASYM(L), RADOBS)
            ELSE
              FLUXDNSEDD = MUINC(1)*.5*RADBNDRYDN1(L) 
     $                 + MUINC(2)*.5*RADBNDRYDN2(L)
              FLUXUPSEDD = MUINC(1)*.5*RADBNDRYUP1(L) 
     $                 + MUINC(2)*.5*RADBNDRYUP2(L)
              RADBNDRYUP(1) = RADBNDRYUP1(L)
              RADBNDRYUP(2) = RADBNDRYUP2(L)
              RADBNDRYDN(1) = RADBNDRYDN1(L)
              RADBNDRYDN(2) = RADBNDRYDN2(L)
              CALL SEDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC, 
     $                       RADBNDRYUP, RADBNDRYDN, 
     $                       FLUXUPSEDD, FLUXDNSEDD, 
     $                       FLUXUP, FLUXDN, RADBOTTOM,
     $                       TAUCG(L), SSALB(L), ASYM(L),
     $                       PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L),
     $                       RADOBS)
            ENDIF

C           The calculated upwelling radiance at top of layer is used
C           as the incident radiance on bottom of next layer.
            RADBOTTOM = RADOBS
            ENDDO

C           Integrate radiance from cloud top to observer
          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
     $                ICLDTOP, IOBS, MUOBS, RADOBS)

C       If looking up must loop over cloud layers from top to bottom  
        ELSEIF(MUOBS.LT.0) THEN
C         Introduce scalar variable for radiance at top of cloud
          RADTOP = RAD0DNOBS(1)

C         Loop over the layers of cloud, from top to bottom
          DO L = 1,NCLDLAY
            TTOP = TEMP(L-1+ICLDTOP)
            TBOT = TEMP(L+ICLDTOP)
            FLUXUP = FLUXES(1,L+1)
            FLUXDN = FLUXES(2,L)
  
C           Do Single-scattering, Eddington, or Hybrid for each layer 
            IF (RTMODEL .EQ. 'S') THEN
              RADBNDRYUP(1) = RADBNDRYUP1(L)
              RADBNDRYUP(2) = RADBNDRYUP2(L)
              RADBNDRYDN(1) = RADBNDRYDN1(L)
              RADBNDRYDN(2) = RADBNDRYDN2(L)
              CALL SSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, MUINC,
     $                    RADBNDRYDN, RADBNDRYUP, RADTOP,
     $                    TAUCG(L), SSALB(L), ASYM(L), 
     $                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L),
     $                    RADOBS)
            ELSE IF (RTMODEL .EQ. 'E') THEN
              CALL EDDSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, 
     $                      FLUXDN, FLUXUP, RADTOP, 
     $                      TAUCG(L), SSALB(L), ASYM(L), RADOBS)
            ELSE
              FLUXDNSEDD = MUINC(1)*.5*RADBNDRYDN1(L) 
     $                 + MUINC(2)*.5*RADBNDRYDN2(L)
              FLUXUPSEDD = MUINC(1)*.5*RADBNDRYUP1(L) 
     $                 + MUINC(2)*.5*RADBNDRYUP2(L)
              RADBNDRYUP(1) = RADBNDRYUP1(L)
              RADBNDRYUP(2) = RADBNDRYUP2(L)
              RADBNDRYDN(1) = RADBNDRYDN1(L)
              RADBNDRYDN(2) = RADBNDRYDN2(L)
              CALL SEDDSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, MUINC, 
     $                       RADBNDRYDN, RADBNDRYUP, 
     $                       FLUXDNSEDD, FLUXUPSEDD, 
     $                       FLUXDN, FLUXUP, RADTOP,
     $                       TAUCG(L), SSALB(L), ASYM(L),
     $                       PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L),
     $                       RADOBS)

              ENDIF

C           The calculated upwelling radiance at bottom of layer is used
C           as the incident radiance on top of next layer.
            RADTOP = RADOBS
            ENDDO
C         Integrate radiance from cloud top to observer
          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
     $                ICLDBOT, IOBS, ABS(MUOBS), RADOBS)
c          print *,'hello',waveno,radobs
          ENDIF
      
        ENDIF         !if clear sky  or if no cloud
c      print *,waveno,rDummy1,rDummy2,RADOBS
c      read *,l

      RETURN
      END


c************************************************************************      
      SUBROUTINE SSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC,
     $                    IBNDRYUP, IBNDRYDN, I0UPOBS,
     $                    TAUTOT, SSALB, ASYM, 
     $                    PHI1UP, PHI1DN, PHI2UP, PHI2DN,  IUPOBS)

      include 'scatter.param'

C     Calculates upwelling radiance from a single thermally emitting
C     homogeneous layer using the single scattering radiative transfer
C     method of Deeter and Evans (1998).
C    Input Parameters:
C       MUOBS      observation mu (cosine of viewing zenith angle)
C       TTOP       cloud top temperature (K)
C       TBOT       cloud bottom temperature (K)
C       WAVENO     wavenumber (cm^-1)
C       MUINC      cosine zenith angles for the two incident directions
C       IBNDRYUP   two incident radiances upwelling on bottom of layer 
C       IBNDRYDN   two incident radiance downwelling on top of layer
C       I0UPOBS    incident radiance at bottom of layer at observation angle
C       TAUTOT     optical depth of layer
C       SSALB      single scattering albedo of layer
C       ASYM       asymmetry parameter of layer
C       PHI*       single scattering functions at the observation angle
C
C    Output Parameters:
C       IUPOBS     outgoing radiance at top at observation angle
C
C      Input and output radiances are in units of W m^-2 sr^-1 cm and
C    fluxes in units of W m^-2 cm.

      REAL    MUOBS, TTOP, TBOT, WAVENO, MUINC(2)
      REAL    IBNDRYUP(2), IBNDRYDN(2), I0UPOBS, IUPOBS
      REAL    TAUTOT, SSALB, ASYM
      REAL    PHI1UP, PHI1DN, PHI2UP, PHI2DN

      INTEGER IMU
      REAL    BETA, COPLNCK, PLANCK0, PLANCK1
      REAL    EXPBETA, EXPMUOBS, DMUOBS, EXPMAX
      REAL    MU, DMUUP, DMUDN, EXPMU
      REAL    S0UP, S1UP, IUP(2), IDN(2)
      REAL    TERM1, TERM2, TERM3

      EXPMAX = 1.E08

      IF (TAUTOT .LT. 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
      ENDIF

C           Planck function radiances in units of W m^-2 sr^-1 cm
      PLANCK0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TTOP) - 1)
      PLANCK1 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TBOT) - 1)
      BETA = LOG(PLANCK1/PLANCK0)/TAUTOT
      COPLNCK = (1. - SSALB)*PLANCK0
      EXPBETA = PLANCK1/PLANCK0
      EXPMUOBS = EXP(-TAUTOT/MUOBS)
      IF (EXPMUOBS .LT. 1./EXPMAX) EXPMUOBS = 1./EXPMAX
      DMUOBS = 1. - BETA*MUOBS
      IF(ABS(DMUOBS).LT.1.E-4) DMUOBS = SIGN(1.E-4,DMUOBS)

C           Zero order scattering contribution
      S0UP = (1. - EXPBETA*EXPMUOBS)*COPLNCK/DMUOBS

C           IMU (1 or 2) is index for discrete mu's in simple model
      DO IMU = 1, 2
        MU = MUINC(IMU)
        EXPMU = EXP(-TAUTOT/MU)
        IF (EXPMU .LT. 1./EXPMAX) EXPMU = 1./EXPMAX
        DMUUP = 1. - BETA*MU
        DMUDN = 1. + BETA*MU
        IF(ABS(DMUUP).LT.1.E-4) DMUUP = SIGN(1.E-4,DMUUP)
        IF(ABS(DMUDN).LT.1.E-4) DMUDN = SIGN(1.E-4,DMUDN)

        TERM1 = IBNDRYUP(IMU)*EXPMU - COPLNCK*EXPBETA*EXPMU/DMUUP
        TERM2 = (MU/(MUOBS - MU))*(1. - EXPMUOBS/EXPMU)
        TERM3 = COPLNCK*(1. - EXPBETA*EXPMUOBS)/(DMUUP*DMUOBS)
        IUP(IMU) = - (TERM1*TERM2 - TERM3)

        TERM1 = IBNDRYDN(IMU) - COPLNCK/DMUDN
        TERM2 = (MU/(MUOBS + MU))*(1. - EXPMUOBS*EXPMU)
        TERM3 =  COPLNCK*(1. - EXPBETA*EXPMUOBS)/(DMUDN*DMUOBS)
        IDN(IMU) = (TERM1*TERM2 + TERM3)
      ENDDO

C           First order scattering contribution
      S1UP = (SSALB/2.)*(PHI1UP*IUP(1) + PHI2UP*IUP(2)
     $         + PHI1DN*IDN(1) + PHI2DN*IDN(2))

      IUPOBS = I0UPOBS*EXPMUOBS + S0UP + S1UP

      RETURN
      END

c************************************************************************
      SUBROUTINE EDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, 
     $                      FLUXUP, FLUXDN, I0UPOBS, 
     $                      TAUTOT, SSALB, ASYM, IUPOBS)

      include 'scatter.param'

C     Calculates upwelling radiance from a single thermally emitting
C     homogeneous layer using Eddington's second approximation 
C     radiative transfer method.
C   
C    Input Parameters:
C       MUOBS      observation mu (cosine of viewing zenith angle)
C       TTOP       cloud top temperature (K)
C       TBOT       cloud bottom temperature (K)
C       WAVENO     wavenumber (cm^-1)
C       FLUXUP     incident flux on bottom from multilayer Eddington solution
C       FLUXDN     incident flux on top from multilayer Eddington solution
C       I0UPOBS    incident radiance at bottom of layer at observation angle
C       TAUTOT     optical depth of layer
C       SSALB      single scattering albedo of layer
C       ASYM       asymmetry parameter of layer
C
C    Output Parameters:
C       IUPOBS     outgoing radiance at top at observation angle
C
C      Input and output radiances are in units of W m^-2 sr^-1 cm and
C    fluxes in units of W m^-2 cm.

      REAL    MUOBS, TTOP, TBOT
      REAL    WAVENO, FLUXUP,FLUXDN
      REAL    I0UPOBS
      REAL    IUPOBS
      
      
      REAL    SSALB, ASYM, TAUTOT
      REAL    BETA, PLANCK0, PLANCK1, COPLANCK
      REAL    EXPBETA, EXPMUOBS, EXPLAMP, EXPLAMM, EXPMAX
      REAL    T, R, LAMBDA
      REAL    F1UP, F0DN, FP1UP, FP0DN, FH1UP, FH0DN
      REAL    PARTFAC, PARTP, PARTM
      REAL    CFAC, CP, CM, DB, D1, D2, DP, DM, SB, SP, SM

      EXPMAX = 1.E08

      IF (TAUTOT .LT. 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
      ENDIF
    
C           Planck function radiances in units of W m^-2 sr^-1 cm

      PLANCK0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TTOP) - 1)
      PLANCK1 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TBOT) - 1)

C           Compute the Planck function variation coefficient 
      BETA = LOG(PLANCK1/PLANCK0)/TAUTOT
      
C      Find lambda so that it can be compared to beta
C      R, T are 2x2 matrix coefficents, lambda is eigenvalue
      T = (7-SSALB*(4+3*ASYM))*0.25
      R = (1-SSALB*(4-3*ASYM))*0.25
C           solution explodes if r = 0!!
      IF (ABS(R) .LT. 1.E-5) R = 1.E-5  
      LAMBDA = SQRT(T**2-R**2)
      
      IF(ABS(BETA**2-LAMBDA**2).LT.1.E-4) THEN
        BETA = SQRT(LAMBDA**2+1.E-4)
      ENDIF

      COPLANCK = (1-SSALB)*PLANCK0
      EXPBETA = PLANCK1/PLANCK0
      EXPMUOBS = EXP(-TAUTOT/MUOBS)
      IF (EXPMUOBS .LT. 1./EXPMAX) EXPMUOBS = 1./EXPMAX

      EXPLAMP = EXP(LAMBDA*TAUTOT)
      IF (EXPLAMP .GT. EXPMAX) EXPLAMP = EXPMAX
      EXPLAMM = 1.0/EXPLAMP

      F1UP = FLUXUP
      F0DN = FLUXDN

C           Compute particular solution at boundaries
      PARTFAC = COPLANCK/(LAMBDA**2-BETA**2)
      PARTP = PARTFAC*(T-R+BETA)
      PARTM = PARTFAC*(T-R-BETA)
      FP0DN = PARTM
      FP1UP = PARTP*EXPBETA
C           Subtract particular solution off to get homogeneous solution
      FH0DN = F0DN - FP0DN
      FH1UP = F1UP - FP1UP
C           Compute homogeneous coefficients C+ and C-
      CFAC = 1.0/(R**2*EXPLAMP - (LAMBDA-T)**2*EXPLAMM)
      CP = CFAC*(R*FH1UP - (LAMBDA-T)*EXPLAMM*FH0DN)
      CM = CFAC*(R*EXPLAMP*FH0DN - (LAMBDA-T)*FH1UP)
C           Make the constants for the 3 source function terms
      DB = COPLANCK + SSALB*PARTFAC*(2*(T-R) + 3*ASYM*BETA*MUOBS)
      D1 = R-T+LAMBDA
      D2 = 1.5*ASYM*(R+T-LAMBDA)*MUOBS
      DP = CP*SSALB*(D1+D2)
      DM = CM*SSALB*(D1-D2)

C           Finally do the answer
      
      IF(ABS(1-BETA*MUOBS).LE.1.E-5) THEN
        SB = DB*TAUTOT/MUOBS
      ELSE
        SB = DB*(1 - EXPBETA*EXPMUOBS)/(1-BETA*MUOBS)
      ENDIF
      
      IF(ABS(1-LAMBDA*MUOBS).LE.1.E-5) THEN
        SP = DP*TAUTOT/MUOBS
      ELSE
        SP = DP*(1 - EXPLAMP*EXPMUOBS)/(1-LAMBDA*MUOBS)
      ENDIF

      SM = DM*(1 - EXPLAMM*EXPMUOBS)/(1+LAMBDA*MUOBS)

      IUPOBS = I0UPOBS*EXPMUOBS + SB + SP + SM

      RETURN
      END

c************************************************************************
      SUBROUTINE SEDDSCATRT (MUOBS, TTOP, TBOT, WAVENO,
     $                       MUINC, IBNDRYUP, IBNDRYDN, 
     $                       FLUXUPSEDD, FLUXDNSEDD, 
     $                       FLUXUP, FLUXDN, I0UPOBS,   
     $                       TAUTOT, SSALB, ASYM,
     $                       PHI1UP, PHI1DN, PHI2UP, PHI2DN,  IUPOBS)

      include 'scatter.param'

C     Calculates upwelling radiance from a single thermally emitting
C     homogeneous layer using the single scattering/Eddington second 
C     approximation hybrid radiative transfer method of Deeter and Evans (1998).
C
C    Input Parameters:
C       MUOBS      observation mu (cosine of viewing zenith angle)
C       TTOP       cloud top temperature (K)
C       TBOT       cloud bottom temperature (K)
C       WAVENO     wavenumber (cm^-1)
C       MUINC      cosine zenith angles for the two incident directions
C       IBNDRYUP   two incident radiances upwelling on bottom of layer 
C       IBNDRYDN   two incident radiance downwelling on top of layer
C       FLUXUPSEDD incident flux on bottom from single scattering/Eddington
C       FLUXDNSEDD incident flux on top from single scattering/Eddington
C       FLUXUP     incident flux on bottom from multilayer Eddington solution
C       FLUXDN     incident flux on top from multilayer Eddington solution
C       I0UPOBS    incident radiance at bottom of layer at observation angle
C       TAUTOT     optical depth of layer
C       SSALB      single scattering albedo of layer
C       ASYM       asymmetry parameter of layer
C       PHI*       single scattering functions at the observation angle
C
C    Output Parameters:
C       IUPOBS     outgoing radiance at top at observation angle
C
C      Input and output radiances are in units of W m^-2 sr^-1 cm and
C    fluxes in units of W m^-2 cm.

      REAL    MUOBS, TTOP, TBOT
      REAL    WAVENO, TAUTOT
      REAL    IBNDRYUP(2), IBNDRYDN(2)
      REAL    I0UPOBS
      REAL    IUPOBS, FLUXUP, FLUXDN, FLUXUPSEDD, FLUXDNSEDD
      REAL    MUINC(2)

      INTEGER IMU
      REAL    SSALB, ASYM, PHI1UP, PHI1DN, PHI2UP, PHI2DN
      REAL    BETA, COPLANCK, PLANCK0, PLANCK1
      REAL    EXPBETA, EXPMUOBS, DMUOBS, EXPMAX
      REAL    MU, DMUUP, DMUDN, EXPMU
      REAL    S1UP, IUP(2), IDN(2)
      REAL    TERM1, TERM2, TERM3
      REAL    EXPLAMP, EXPLAMM
      REAL    T, R, LAMBDA
      REAL    F1UP, F0DN, FP1UP, FP0DN, FH1UP, FH0DN
      REAL    PARTFAC, PARTP, PARTM
      REAL    CFAC, CP, CM, DB, D1, D2, DP, DM, SB, SP, SM

      REAL    SQR3, EXP3TP, EXP3TM, FAC3TP, FAC3TM, FAC3BETA
      REAL    FACF0DN, FACF1UP, FACG, FACGP, FACGM
      REAL    FACDENOM, DENOM1, DENOM2, FACNUM1, FACNUM2
      REAL    SEDTERM1, SEDTERM2, FACMUBET, FACNUM3, SEDTERM3, S1EDD
     
      IF (TAUTOT .LT. 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
      ENDIF

      EXPMAX = 1.E08
      SQR3 = SQRT(3.0)

C           Planck function radiances in units of W m^-2 sr^-1 cm

      PLANCK0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TTOP) - 1)
      PLANCK1 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TBOT) - 1)
      BETA = LOG(PLANCK1/PLANCK0)/TAUTOT

C      Find lambda so that it can be compared to beta
C      R, T are 2x2 matrix coefficents, lambda is eigenvalue
      T = (7-SSALB*(4+3*ASYM))*0.25
      R = (1-SSALB*(4-3*ASYM))*0.25
C           solution explodes if r = 0!!
      IF (ABS(R) .LT. 1.E-5) R = 1.E-5
      LAMBDA = SQRT(T**2-R**2)

      IF(ABS(BETA**2-LAMBDA**2).LT.1.E-4) THEN
        BETA = SQRT(LAMBDA**2+1.E-4)
      ENDIF

      COPLANCK = (1. - SSALB)*PLANCK0
      EXPBETA = PLANCK1/PLANCK0
      EXPMUOBS = EXP(-TAUTOT/MUOBS)
      IF (EXPMUOBS .LT. 1./EXPMAX) EXPMUOBS = 1./EXPMAX
      DMUOBS = 1. - BETA*MUOBS
      IF(ABS(DMUOBS).LT.1.E-4) DMUOBS = SIGN(1.E-4,DMUOBS)
 
C         Accurate Single Scattering Part
C           IMU (1 or 2) is index for discrete mu's in simple model
      DO IMU = 1, 2
        MU = MUINC(IMU)
        EXPMU = EXP(-TAUTOT/MU)
        IF (EXPMU .LT. 1./EXPMAX) EXPMU = 1./EXPMAX
        DMUUP = 1. - BETA*MU
        DMUDN = 1. + BETA*MU
        IF(ABS(DMUUP).LT.1.E-4) DMUUP = SIGN(1.E-4,DMUUP)
        IF(ABS(DMUDN).LT.1.E-4) DMUDN = SIGN(1.E-4,DMUDN)

        TERM1 = IBNDRYUP(IMU)*EXPMU - COPLANCK*EXPBETA*EXPMU/DMUUP
        TERM2 = (MU/(MUOBS - MU))*(1. - EXPMUOBS/EXPMU)
        TERM3 = COPLANCK*(1. - EXPBETA*EXPMUOBS)/(DMUUP*DMUOBS)
        IUP(IMU) = - (TERM1*TERM2 - TERM3)
 
        TERM1 = IBNDRYDN(IMU) - COPLANCK/DMUDN
        TERM2 = (MU/(MUOBS + MU))*(1. - EXPMUOBS*EXPMU)
        TERM3 =  COPLANCK*(1. - EXPBETA*EXPMUOBS)/(DMUDN*DMUOBS)
        IDN(IMU) = (TERM1*TERM2 + TERM3)
      ENDDO

C           First order scattering contribution
      S1UP = (SSALB/2.)*(PHI1UP*IUP(1) + PHI2UP*IUP(2)
     $         + PHI1DN*IDN(1) + PHI2DN*IDN(2))


C         Eddington part

      EXPLAMP = EXP(LAMBDA*TAUTOT)
      IF (EXPLAMP .GT. EXPMAX) EXPLAMP = EXPMAX
      EXPLAMM = 1.0/EXPLAMP

C           Fluxes sent into subroutine from output of multilayer solution
      F1UP = FLUXUP
      F0DN = FLUXDN

C           Compute particular solution at boundaries
      PARTFAC = COPLANCK/(LAMBDA**2-BETA**2)
      PARTP = PARTFAC*(T-R+BETA)
      PARTM = PARTFAC*(T-R-BETA)
      FP0DN = PARTM
      FP1UP = PARTP*EXPBETA
C           Subtract particular solution off to get homogeneous solution
      FH0DN = F0DN - FP0DN
      FH1UP = F1UP - FP1UP
C           Compute homogeneous coefficients C+ and C-
      CFAC = 1.0/(R**2*EXPLAMP - (LAMBDA-T)**2*EXPLAMM)
      CP = CFAC*(R*FH1UP - (LAMBDA-T)*EXPLAMM*FH0DN)
      CM = CFAC*(R*EXPLAMP*FH0DN - (LAMBDA-T)*FH1UP)
C           Make the constants for the 3 source function terms
      DB = COPLANCK + SSALB*PARTFAC*(2*(T-R) + 3*ASYM*BETA*MUOBS)
      D1 = R-T+LAMBDA
      D2 = 1.5*ASYM*(R+T-LAMBDA)*MUOBS
      DP = CP*SSALB*(D1+D2)
      DM = CM*SSALB*(D1-D2)

C           Finally do the answer

      IF(ABS(1-BETA*MUOBS).LE.1.E-5) THEN
        SB = DB*TAUTOT/MUOBS
      ELSE
        SB = DB*(1 - EXPBETA*EXPMUOBS)/(1-BETA*MUOBS)
      ENDIF
      
      IF(ABS(1-LAMBDA*MUOBS).LE.1.E-5) THEN
        SP = DP*TAUTOT/MUOBS
      ELSE
        SP = DP*(1 - EXPLAMP*EXPMUOBS)/(1-LAMBDA*MUOBS)
      ENDIF

      SM = DM*(1 - EXPLAMM*EXPMUOBS)/(1+LAMBDA*MUOBS)

C           Compute the Eddington single scattering solution
      F0DN = FLUXDNSEDD
      F1UP = FLUXUPSEDD
      EXP3TP = EXP(SQR3*TAUTOT)
      IF (EXP3TP .GT. EXPMAX) EXP3TP = EXPMAX
      EXP3TM = 1./EXP3TP
      FAC3TP = 1 - EXPMUOBS*EXP3TP
      FAC3TM = 1 - EXPMUOBS*EXP3TM
      FAC3BETA = 1./(3. - BETA**2)
      FACF0DN = F0DN - (1.5 - BETA)*COPLANCK*FAC3BETA
      FACF1UP = F1UP - (1.5 + BETA)*COPLANCK*EXPBETA*FAC3BETA
      FACG = 1.5*(2. - SQR3)*ASYM*MUOBS
      FACGP = -1.5 + SQR3 + FACG
      FACGM = -1.5 + SQR3 - FACG
      FACDENOM = EXP3TM*(-7./4. + SQR3)**2 - EXP3TP/16.
      DENOM1 =  FACDENOM*(1 + MUOBS*SQR3)
      DENOM2 = -FACDENOM*(1 - MUOBS*SQR3)
      FACNUM1 = -.25*EXP3TP*FACF0DN + (-7./4.+ SQR3)*FACF1UP
      FACNUM2 = -(-7./4. + SQR3)*EXP3TM*FACF0DN + .25*FACF1UP
      SEDTERM1 = FAC3TM*FACNUM1*FACGM/DENOM1
      SEDTERM2 = FAC3TP*FACNUM2*FACGP/DENOM2
      FACMUBET = 1. - EXPBETA*EXPMUOBS
      FACNUM3 = COPLANCK*3.*(1 + BETA*ASYM*MUOBS)*FAC3BETA
      SEDTERM3 = FACMUBET*FACNUM3/DMUOBS

      S1EDD = SSALB*(SEDTERM1 + SEDTERM2 + SEDTERM3)

C           Zero order scattering contribution
C       S0UP = (1. - EXPBETA*EXPMUOBS)*COPLANCK/DMUOBS

      IUPOBS = I0UPOBS*EXPMUOBS + SB + SP + SM - S1EDD + S1UP

      RETURN
      END


c************************************************************************
       SUBROUTINE EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, 
     $                    ASYMMETRIES, TEMPS, WAVENO, 
     $                    FLUXTOP, FLUXBOT,  FLUXES)

      include 'scatter.param'

C       EDDRTF computes the layer interface fluxes for a multilayer 
C     plane-parallel atmosphere with thermal sources of radiation using the
C     Eddington approximation.  The medium is specified by a number 
C     of homogeneous layers.  The Planck function is assumed linear with 
C     optical depth (slightly inconsistent with the exponential assumption
C     in the single layer routines).  The temperatures, optical depth, 
C     single scattering albedo, and asymmetry parameter are specified 
C     for each layer.  The boundary conditions are the fluxes incident
C     on the top and bottom of the domain. The Eddington fluxes at each 
C     level are returned.
C       The model works by calculating the reflection, transmission, and
C     source terms for each layer from the input properties.  A
C     tri-diagonal matrix solver is then used to compute the fluxes 
C     at each layer interface from the applied boundary conditions.
C
C     Parameters:
C       Input:
C     NLAYER         integer      Number of homogenous layers
C                                  (layers are specified from the top down)
C     OPTDEPTHS      real array   Optical thickness of layers
C     ALBEDOS        real array   Single scattering albedos
C     ASYMMETRIES    real array   Asymmetry parameters
C     TEMPS          real array   Temperatures (K) at layer interfaces
C                                  (e.g. TEMPS(1) is at top of top layer, 
C                                   TEMPS(2) is at bottom of top layer).
C     WAVENO         real         wavenumber (cm^-1, for Planck function)
C
C       Output:
C     FLUXES         real         Eddington fluxes at layer interfaces.
C                                   FLUXES(1,L) is upwelling,
C                                   FLUXES(2,L) is downwelling,
C                                   L=1 is top, L=NUML+1 is bottom

      INTEGER   NLAYER
      REAL      TEMPS(NLAYER+1)
      REAL      OPTDEPTHS(NLAYER)
      REAL      ALBEDOS(NLAYER)
      REAL      ASYMMETRIES(NLAYER)
      REAL      WAVENO
      REAL      FLUXES(3,(NLAYER+1))
 
C               MAXLAY is maximum number of layers
      INTEGER   MAXLAY, MAXN
c      PARAMETER (MAXLAY=200, MAXN=2*MAXLAY+1)
      PARAMETER (MAXLAY=kProfLayer, MAXN=2*MAXLAY+1)
      INTEGER   N, L, I
      REAL*8    DELTAU, G, OMEGA
      REAL*8    LAMBDA, R, T, D, CM, CP, X1, X2
      REAL*8    REFLECT, TRANS, SOURCEP, SOURCEM
      REAL*8    RADP1P, RADP1M, RADP2P, RADP2M
      REAL*8    PI, PLANCK1, PLANCK2, TAU
      REAL*8    EXLP, EXLM, V
      REAL*8    LOWER(MAXN), UPPER(MAXN), DIAG(MAXN), RHS(MAXN)
      REAL      FLUXTOP, FLUXBOT
      PARAMETER (PI=3.1415926535)
 
C               Compute the reflection, transmission, and source
C               coefficients for each layer for the diffuse Eddington
C               two stream problem.


      N = 2*NLAYER+2
      IF (N .GT. MAXN) STOP 'EDDRTF: Exceeded maximum number of layers'

      PLANCK1 = 0.5*KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMPS(1)) - 1)
      TAU = 0.0
      I = 2
      DO L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        IF (DELTAU .LT. 0.0) STOP 'EDDRTF: TAU<0'
C            Special case for zero optical depth
        IF (DELTAU .EQ. 0.0) THEN
          TRANS = 1.0
          REFLECT = 0.0
          SOURCEP = 0.0
          SOURCEM = 0.0
        ELSE
          OMEGA = ALBEDOS(L)
          G = ASYMMETRIES(L)
          R = ( 1.0 - OMEGA*(4.0-3.0*G) )/4.0
          T = ( 7.0 - OMEGA*(4.0+3.0*G) )/4.0
          LAMBDA = SQRT( 3.0*(1.0-OMEGA)*(1.0-OMEGA*G) )
C              Special case for conservative scattering (lambda=0)
          IF (LAMBDA .EQ. 0.0) THEN
            D = 1.0/(1.0+T*DELTAU)
            TRANS = D
            REFLECT = -R*DELTAU*D
          ELSE
            X1 = -R
            X2 = LAMBDA + T
            EXLP = DEXP(MIN(LAMBDA*DELTAU,75.D0))
            EXLM = 1.0/EXLP
            TRANS = 2.*LAMBDA/(X2*EXLP + (LAMBDA-T)*EXLM)
            REFLECT = X1*(EXLP - EXLM) *TRANS /(2.*LAMBDA)
            D = 1.0/(X2**2 *EXLP - X1**2 *EXLM)
          ENDIF
          
C               Calculate thermal source terms
           
            PLANCK2 = 0.5*KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMPS(L+1)) - 1)

            V = 2.0*(PLANCK2-PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
            RADP1P = -V + PLANCK1
            RADP2M =  V + PLANCK2
            RADP2P = -V + PLANCK2
            RADP1M =  V + PLANCK1
            IF (LAMBDA .EQ. 0.0) THEN
              SOURCEP = 0
              SOURCEM = 0
            ELSE
              CP  =  (X1*EXLM*RADP1P - X2*RADP2M) *D
              CM = (-X2*EXLP*RADP1P + X1*RADP2M) *D
              SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
              SOURCEM = X2*CP + X1*CM + RADP1M
            ENDIF
            PLANCK1 = PLANCK2
            FLUXES(3,L) = 0.0
        ENDIF
        DIAG(I) = -REFLECT
        DIAG(I+1) = -REFLECT
        LOWER(I) = 1.0
        LOWER(I+1) = -TRANS
        UPPER(I) = -TRANS
        UPPER(I+1) = 1.0
        RHS(I) = SOURCEM
        RHS(I+1) = SOURCEP
        I = I + 2
      ENDDO
      
C           Setup for and call the tri-diagonal matrix solver
      RHS(1) = FLUXTOP
      DIAG(1) = 0.0
      UPPER(1) = 1.0
c      DIAG(N) = -(1.0-GNDEMIS)
      DIAG(N) = 0.0
      LOWER(N) = 1.0
      RHS(N) = FLUXBOT
      CALL TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
C           Put the fluxes in the output array
     
      I = 1
      DO L = 1, NLAYER+1 
        FLUXES(1,L) = RHS(I)
        FLUXES(2,L) = RHS(I+1)
        I = I + 2
      ENDDO
 
      RETURN
      END
 
c************************************************************************
      SUBROUTINE TRIDIAG (N, LOWER, DIAG, UPPER, RHS)


C       Computes the solution to a tridiagonal system. 
C       N is order of the matrix.  LOWER(2..N) is the subdiagonal,
C       DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the 
C       superdiagonal.  On input RHS is the right hand side, while
C       on output it is the solution vector.  Everything is destroyed.
C       Hacked from Linpack DGTSL.

      include 'scatter.param'

      INTEGER N 
      REAL*8  LOWER(*), DIAG(*), UPPER(*), RHS(*)
      INTEGER K, KB
      REAL*8  T

      IF (N .EQ. 1) THEN
        IF (DIAG(1) .EQ. 0.0) GOTO 990
        RHS(1) = RHS(1)/DIAG(1)
      ENDIF
      LOWER(1) = DIAG(1)
      DIAG(1) = UPPER(1)
      UPPER(1) = 0.0
      UPPER(N) = 0.0
      DO K = 1, N-1
C              Interchange this and next row to the get the largest pivot.
        IF (ABS(LOWER(K+1)) .GE. ABS(LOWER(K))) THEN
          T = LOWER(K+1)
          LOWER(K+1) = LOWER(K)
          LOWER(K) = T
          T = DIAG(K+1)
          DIAG(K+1) = DIAG(K)
          DIAG(K) = T
          T = UPPER(K+1)
          UPPER(K+1) = UPPER(K)
          UPPER(K) = T
          T = RHS(K+1)
          RHS(K+1) = RHS(K)
          RHS(K) = T
        ENDIF
        IF (LOWER(K) .EQ. 0.0) GOTO 990
        T = -LOWER(K+1)/LOWER(K)
        LOWER(K+1) = DIAG(K+1) + T*DIAG(K)
        DIAG(K+1) = UPPER(K+1) + T*UPPER(K)
        UPPER(K+1) = 0.0
        RHS(K+1) = RHS(K+1) + T*RHS(K)
      ENDDO
      IF (LOWER(N) .EQ. 0.0) GOTO 990

C           Back substitute
      RHS(N) = RHS(N)/LOWER(N)
      RHS(N-1) = (RHS(N-1) - DIAG(N-1)*RHS(N))/LOWER(N-1)
      DO KB = 1, N-2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) -DIAG(K)*RHS(K+1) -UPPER(K)*RHS(K+2))/LOWER(K)
      ENDDO
      RETURN

990   CONTINUE
        STOP 'Singular matrix in TRIDIAG'
      END

c************************************************************************
c slight modifications to GASRT1 done by Sergio, to allow backgnd thermal to 
c be computed or not; Planck terms stored in array instead of computing twice
c this computes backgnd thermal using fast diffusivity approx

      SUBROUTINE GASRT1_nocloud (MU, WAVENO, SFCTEMP, SFCEMIS,
     $                NLEV, TEMP, TAU, COALB, ICLDTOP, 
     $                ICLDBOT, NCLDLAY, RAD1UP, RAD0DN,
     $                ibdry, TOA_to_instr)
C       Compute the upwelling radiance at the bottom of each cloud layer
C     (RAD1UP) and the downwelling radiance at the top of each cloud
C     layer (RAD0DN).  Integrates the nonscattering RTE through the domain
C     from top to bottom and then bottom to top.  Only the absorption
C     part of the optical depth in the cloud is considered.  The bottom
C     of the lowest cloud layer is at level ICLDBOT and the top of the
C     highest cloud layer is as level ICLDTOP.  The radiance is computed 
C     at the angle MU and wavenumber WAVENO.  The temperature profile (TEMP),
C     layer optical depth profile (TAU), and single scattering coalbedo
C     (COALB=1-omega) profile are input.  There is specular reflection from
C     the surface yet no multiple scattering effect due of the surface, ie.
C     The gas integrated down from TOA and reflected off surface once with
C     reflectivity equal to one minus emissivity, and this is the initial
C     condition when integrating gas incident on bottom of cloud. 

      include 'scatter.param'

      INTEGER NLEV, ICLDTOP, ICLDBOT, NCLDLAY,ibdry
      REAL    MU, WAVENO, SFCTEMP, SFCEMIS
      REAL    TEMP(NLEV), TAU(NLEV), COALB(NLEV)
      REAL    RAD1UP(NCLDLAY), RAD0DN(NCLDLAY)
      REAL    RAD1UP0, RAD0DN0, TOA_to_instr
      INTEGER I, L
      REAL    PLANCK0, PLANCK1, PLANCKSFC, TAUO, TRANS, DELPLANCK

      REAL raPlanck(kProfLayer+1)

      DO i=1,NLEV
        raPlanck(i)=KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(I)) - 1)
        END DO

C     Loop through layers, starting at top and going down
C     Continue past cloud layer to surface to compute specular
C     reflection.

      IF (kThermal .GE. 0) THEN
       CALL FastBDRYL2GDiffusiveApprox_rtspec(TOA_to_instr, 
     $                MU, WAVENO,NLEV, TEMP, TAU,RAD0DN, RAD0DN0,ibdry)  
        END IF

C     Loop through layers, starting at bottom and going up
c     note because of defn of emissivity used in FastBdry, we multiply
c     rad0dn0 by (1-ems)/pi  
      L = NCLDLAY
      PLANCKSFC = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/SFCTEMP) - 1)

      IF (kThermal .GE. 0) THEN
        RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)/kPi
      ELSE
        RAD1UP0 = SFCEMIS*PLANCKSFC
        END IF

cx      print *,'surf',rad1up0
      IF(ICLDTOP.EQ.(NLEV-1)) RAD1UP(1) = RAD1UP0

c notice if DELPLANCK == 0 then below algorithm reduces to 
c rad=sum( (1-tau(i))B(i) + rad*tau(i)) 
c which is the same as that of the usual kCARTA algorithm
cc        DO iFr=1,kMaxPts
cc          raInten(iFr)=raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc          END DO
cc raaEmission = (1-tau(i))B(i), raaLayTran = tau(i)
c except that there the Planck term B(i) is at the "middle" of the layer,
c while here it is at the top of the layer, which is colder .. hence this 
c overall radiance should be smaller than than of kCARTA!

c so if DELPLANK is non zero, but actually positive, this small correction
c will bring up back to kCARTA radiance levels!!! woohoo
      PLANCK1 = raPlanck(NLEV)
      DO I = NLEV-1, ICLDTOP+1, -1
        PLANCK0 = raPlanck(I)
        TAUO = TAU(I)/MU
        IF (TAUO .LT. 0.001) THEN
          RAD1UP0 = RAD1UP0*(1-TAUO) 
     $              + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I) 
        ELSE
          TRANS = EXP(-TAUO)
          DELPLANCK = (PLANCK1-PLANCK0)/TAUO
          RAD1UP0 = RAD1UP0*TRANS + COALB(I)*( PLANCK0+DELPLANCK
     $                  - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO)) )
        ENDIF
        IF (I .LE. ICLDBOT) THEN
          RAD1UP(L) = RAD1UP0
          L = L - 1
        ENDIF
        PLANCK1 = PLANCK0
cx        print *,i,rad1up0
      ENDDO    

cx      read *,l
      RETURN
      END

c************************************************************************
c slight modifications to GASRT1 done by Sergio, to allow backgnd thermal to 
c be computed or not; Planck terms stored in array instead of computing twice
c this computes backgnd thermal using acos(3/5) at every level

      SUBROUTINE GASRT1_nocloud_works (MU, WAVENO, SFCTEMP, SFCEMIS,
     $                NLEV, TEMP, TAU, COALB, ICLDTOP, 
     $                ICLDBOT, NCLDLAY, RAD1UP, RAD0DN,
     $                TOA_to_instr)
C       Compute the upwelling radiance at the bottom of each cloud layer
C     (RAD1UP) and the downwelling radiance at the top of each cloud
C     layer (RAD0DN).  Integrates the nonscattering RTE through the domain
C     from top to bottom and then bottom to top.  Only the absorption
C     part of the optical depth in the cloud is considered.  The bottom
C     of the lowest cloud layer is at level ICLDBOT and the top of the
C     highest cloud layer is as level ICLDTOP.  The radiance is computed 
C     at the angle MU and wavenumber WAVENO.  The temperature profile (TEMP),
C     layer optical depth profile (TAU), and single scattering coalbedo
C     (COALB=1-omega) profile are input.  There is specular reflection from
C     the surface yet no multiple scattering effect due of the surface, ie.
C     The gas integrated down from TOA and reflected off surface once with
C     reflectivity equal to one minus emissivity, and this is the initial
C     condition when integrating gas incident on bottom of cloud. 

      include 'scatter.param'

      INTEGER NLEV, ICLDTOP, ICLDBOT, NCLDLAY
      REAL    MU, WAVENO, SFCTEMP, SFCEMIS
      REAL    TEMP(NLEV), TAU(NLEV), COALB(NLEV)
      REAL    RAD1UP(NCLDLAY), RAD0DN(NCLDLAY)
      REAL    RAD1UP0, RAD0DN0, TOA_to_instr
      INTEGER I, L, IEND
      REAL    PLANCK0, PLANCK1, PLANCKSFC, TAUO, TRANS, DELPLANCK

      REAL raPlanck(kProfLayer+1),MUDOWN

      DO i=1,NLEV
        raPlanck(i)=KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(I)) - 1)
        END DO

C     Loop through layers, starting at top and going down
C     Continue past cloud layer to surface to compute specular
C     reflection.

      IF (kThermal .GE. 0) THEN
        RAD0DN0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/kTSpace) - 1)
        IF(ICLDBOT.EQ.2) RAD0DN(1) = RAD0DN0

        PLANCK0 = raPlanck(1)

        L = 1
        IF(SFCEMIS.LT.1) IEND = NLEV-1
        IF(SFCEMIS.EQ.1) IEND = ICLDBOT-2

        MUDOWN=3./5.      
        DO I = 1, IEND
C         Planck function radiances in units of W m^-2 sr^-1 cm
          PLANCK1 = raPlanck(I+1)
          TAUO = TAU(I)/MUDOWN
          IF (TAU(I) .LT. 0.001) THEN
            RAD0DN0 = RAD0DN0*(1-TAUO)
     $              + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
          ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAU(I)
            RAD0DN0 = RAD0DN0*TRANS + COALB(I)*( PLANCK1-DELPLANCK
     $                  - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
          ENDIF
          IF (I+1 .GE. ICLDTOP . AND. I .LE. (ICLDBOT-2)) THEN 
            RAD0DN(L) = RAD0DN0
            L = L + 1  
            ENDIF
          PLANCK0 = PLANCK1
          ENDDO    
        END IF

C     Loop through layers, starting at bottom and going up
      L = NCLDLAY
      PLANCKSFC = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/SFCTEMP) - 1)

      IF (kThermal .GE. 0) THEN
        RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)
      ELSE
        RAD1UP0 = SFCEMIS*PLANCKSFC
        END IF

cx      print *,'surf',rad1up0
      IF(ICLDTOP.EQ.(NLEV-1)) RAD1UP(1) = RAD1UP0

c notice if DELPLANCK == 0 then below algorithm reduces to 
c rad=sum( (1-tau(i))B(i) + rad*tau(i)) 
c which is the same as that of the usual kCARTA algorithm
cc        DO iFr=1,kMaxPts
cc          raInten(iFr)=raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc          END DO
cc raaEmission = (1-tau(i))B(i), raaLayTran = tau(i)
c except that there the Planck term B(i) is at the "middle" of the layer,
c while here it is at the top of the layer, which is colder .. hence this 
c overall radiance should be smaller than than of kCARTA!

c so if DELPLANK is non zero, but actually positive, this small correction
c will bring up back to kCARTA radiance levels!!! woohoo
      PLANCK1 = raPlanck(NLEV)
      DO I = NLEV-1, ICLDTOP+1, -1
        PLANCK0 = raPlanck(I)
        TAUO = TAU(I)/MU
        IF (TAUO .LT. 0.001) THEN
          RAD1UP0 = RAD1UP0*(1-TAUO) 
     $              + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I) 
        ELSE
          TRANS = EXP(-TAUO)
          DELPLANCK = (PLANCK1-PLANCK0)/TAUO
          RAD1UP0 = RAD1UP0*TRANS + COALB(I)*( PLANCK0+DELPLANCK
     $                  - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO)) )
        ENDIF
        IF (I .LE. ICLDBOT) THEN
          RAD1UP(L) = RAD1UP0
          L = L - 1
        ENDIF
        PLANCK1 = PLANCK0
cx        print *,i,rad1up0
      ENDDO    

cx      read *,l
      RETURN
      END

c************************************************************************
c slight modifications to GASRT1 done by Sergio so that Planck terms 
c stored in array instead of computing twice

      SUBROUTINE GASRT1 (MU, WAVENO, SFCTEMP, SFCEMIS,
     $                NLEV, TEMP, TAU, COALB, ICLDTOP, 
     $               ICLDBOT, NCLDLAY, RAD1UP, RAD0DN)
C       Compute the upwelling radiance at the bottom of each cloud layer
C     (RAD1UP) and the downwelling radiance at the top of each cloud
C     layer (RAD0DN).  Integrates the nonscattering RTE through the domain
C     from top to bottom and then bottom to top.  Only the absorption
C     part of the optical depth in the cloud is considered.  The bottom
C     of the lowest cloud layer is at level ICLDBOT and the top of the
C     highest cloud layer is as level ICLDTOP.  The radiance is computed 
C     at the angle MU and wavenumber WAVENO.  The temperature profile (TEMP),
C     layer optical depth profile (TAU), and single scattering coalbedo
C     (COALB=1-omega) profile are input.  There is specular reflection from
C     the surface yet no multiple scattering effect due of the surface, ie.
C     The gas integrated down from TOA and reflected off surface once with
C     reflectivity equal to one minus emissivity, and this is the initial
C     condition when integrating gas incident on bottom of cloud. 

      include 'scatter.param'

      INTEGER NLEV, ICLDTOP, ICLDBOT, NCLDLAY
      REAL    MU, WAVENO, SFCTEMP, SFCEMIS
      REAL    TEMP(NLEV), TAU(NLEV), COALB(NLEV)
      REAL    RAD1UP(NCLDLAY), RAD0DN(NCLDLAY)
      REAL    RAD1UP0, RAD0DN0
      INTEGER I, L, IEND
      REAL    PLANCK0, PLANCK1, PLANCKSFC, TAUO, TRANS, DELPLANCK

      REAL raPlanck(kProfLayer+1)

      DO i=1,NLEV
        raPlanck(i)=KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(I)) - 1)
        END DO

C         Loop through layers, starting at top and going down
C      Continue past cloud layer to surface to compute specular
C      reflection.

      RAD0DN0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/kTSpace) - 1)
      IF(ICLDBOT.EQ.2) RAD0DN(1) = RAD0DN0

      PLANCK0 = raPlanck(1)

      L = 1
      IF(SFCEMIS.LT.1) IEND = NLEV-1
      IF(SFCEMIS.EQ.1) IEND = ICLDBOT-2
      
      DO I = 1, IEND
C           Planck function radiances in units of W m^-2 sr^-1 cm
        PLANCK1 = raPlanck(I+1)
        TAUO = TAU(I)/MU
        IF (TAU(I) .LT. 0.001) THEN
          RAD0DN0 = RAD0DN0*(1-TAUO)
     $              + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
        ELSE
          TRANS = EXP(-TAUO)
          DELPLANCK = (PLANCK1-PLANCK0)/TAU(I)
          RAD0DN0 = RAD0DN0*TRANS + COALB(I)*( PLANCK1-DELPLANCK
     $                  - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
        ENDIF
        IF (I+1 .GE. ICLDTOP . AND. I .LE. (ICLDBOT-2)) THEN 
          RAD0DN(L) = RAD0DN0
          L = L + 1  
        ENDIF
        PLANCK0 = PLANCK1
c      print *,i,waveno,RAD0DN0,temp(i),tau(i)
      ENDDO    
c      read *,idummy

C     Loop through layers, starting at bottom and going up
      L = NCLDLAY
      PLANCKSFC = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/SFCTEMP) - 1)

      RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)

      IF(ICLDTOP.EQ.(NLEV-1)) RAD1UP(1) = RAD1UP0

      PLANCK1 = raPlanck(NLEV)
      DO I = NLEV-1, ICLDTOP+1, -1
        PLANCK0 = raPlanck(I)
        TAUO = TAU(I)/MU
        IF (TAUO .LT. 0.001) THEN
          RAD1UP0 = RAD1UP0*(1-TAUO) 
     $              + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I) 
        ELSE
          TRANS = EXP(-TAUO)
          DELPLANCK = (PLANCK1-PLANCK0)/TAUO
          RAD1UP0 = RAD1UP0*TRANS + COALB(I)*( PLANCK0+DELPLANCK
     $                  - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO)) )
        ENDIF
cx        print *,I,TEMP(I),COALB(I),RAD1UP0
        IF (I .LE. ICLDBOT) THEN
          RAD1UP(L) = RAD1UP0
          L = L - 1
        ENDIF

        PLANCK1 = PLANCK0
      ENDDO    

      RETURN
      END

c************************************************************************
      SUBROUTINE GASRT2 (WAVENO, NLEV, TEMP, TAU, ICLD, IOBS, MU, RAD)
C       Compute the upwelling radiance starting at the cloud top (ICLD)
C     and going to the observation level (IOBS) or the downwelling 
C     radiance starting at cloud bottom (ICLD) and going down to 
C     observation level (IOBS) depending on MU. 
C     The radiance is computed at the angle MU and wavenumber WAVENO.
C     The temperature profile (TEMP) and layer optical depth profile (TAU)
C     are input.

      include 'scatter.param'

      INTEGER NLEV, ICLD, IOBS
      REAL    WAVENO, TEMP(NLEV), TAU(NLEV), MU, RAD
      INTEGER I
      REAL    PLANCK0, PLANCK1, TAUO, TRANS, DELPLANCK

C       Loop through layers, if looking down start at bottom and go up
      IF(IOBS.LT.ICLD) THEN
        PLANCK1 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(ICLD)) - 1)
        DO I = ICLD-1, IOBS, -1
          PLANCK0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(I)) - 1)
          TAUO = TAU(I)/MU
          IF (TAUO .LT. 0.001) THEN
            RAD = RAD*(1-TAUO) + TAUO*0.5*(PLANCK0+PLANCK1) 
          ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAUO
            RAD = RAD*TRANS + PLANCK0+DELPLANCK
     $            - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO))
          ENDIF
          PLANCK1 = PLANCK0
        ENDDO    

C       Loop through layers, if looking up start at top and go down

      ELSEIF(IOBS.GT.ICLD) THEN
        PLANCK0 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(ICLD)) - 1)
        DO I = ICLD, IOBS-1
C           Planck function radiances in units of W m^-2 sr^-1 um^-1
          PLANCK1 = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/TEMP(I+1)) - 1)
          TAUO = TAU(I)/MU
          IF (TAU(I) .LT. 0.001) THEN
            RAD = RAD*(1-TAUO)
     $              + TAUO*0.5*(PLANCK0+PLANCK1)
          ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAU(I)
            RAD = RAD*TRANS + ( PLANCK1-DELPLANCK
     $                  - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
          ENDIF
          PLANCK0 = PLANCK1
        ENDDO    
      ENDIF

      RETURN
      END

c************************************************************************
      SUBROUTINE READ_SSCATTAB(SCATFILE,               !!!  MAXTAB, MAXGRID,
     $          NMUOBS, MUTAB, NDME, DMETAB, NWAVE, WAVETAB,
     $          MUINC, TABEXTINCT, TABSSALB, TABASYM,
     $          TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

      include 'scatter.param'

C       Reads in the single scattering table for a number of wavenumbers,
C     particle sizes, and viewing angles.  The scattering properties are
C     computed for a IWC/LWC of 1 g/m^3.  
C       Input parameters:
C     SCATFILE   file name of scattering file
C     MAXTAB     maximum array size for whole table
C     MAXGRID    maximum array size for table grids
C       Output parameters:
C     NMUOBS     number of viewing angle mu grid values
C     MUTAB      viewing angle grid values
C     NDME       number of Dme grid values
C     DMETAB     Dme grid values
C     NWAVE      number of wavenumber grid values
C     WAVETAB    wavenumber grid values
C     MUINC(2)   cosine zenith of two incident angles
C     TABEXTINCT tabulated extinction (km^-1)
C     TABSSALB   tabulated single scattering albedo 
C     TABASYM    tabulated asymmetry parameter
C     TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN  
C                tabulated phase function terms for incident radiance angles

c!!!      INTEGER  MAXTAB, MAXGRID
      INTEGER  NMUOBS, NDME, NWAVE
      REAL     MUTAB(*), DMETAB(*), WAVETAB(*)
      REAL     MUINC(2)
      REAL     TABEXTINCT(*), TABSSALB(*), TABASYM(*) 
      REAL     TABPHI1UP(*), TABPHI1DN(*)
      REAL     TABPHI2UP(*), TABPHI2DN(*)
      CHARACTER*(*) SCATFILE
      INTEGER  IMU, ID, IW, K2, K3

      OPEN (UNIT=2, STATUS='OLD', FILE=SCATFILE)
      READ (2,*) 
      READ (2,*) NMUOBS
      READ (2,*) NDME
      READ (2,*) NWAVE
      IF (MAX(NMUOBS,NDME,NWAVE) .GT. MAXGRID) 
     $    STOP 'READ_SSCATTAB: MAXGRID exceeded'
      IF (NMUOBS*NDME*NWAVE .GT. MAXTAB) 
     $    STOP 'READ_SSCATTAB: MAXTAB exceeded'
      READ (2,*) MUINC(1), MUINC(2)
      READ (2,*) 
      READ (2,*) 
      READ (2,*)
      READ (2,*)
      READ (2,*) 
      READ (2,*) (MUTAB(IMU), IMU = 1, NMUOBS)         
      DO IW = 1, NWAVE
        DO ID = 1, NDME
            K2 = IW-1 + NWAVE*(ID-1)
            K3 = NMUOBS*K2
            READ(2,*)
            READ(2,*)
            READ(2,*) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1),
     $       TABSSALB(K2+1), TABASYM(K2+1)
            READ(2,*) (TABPHI1UP(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI2UP(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI1DN(IMU+K3), IMU = 1, NMUOBS)
            READ(2,*) (TABPHI2DN(IMU+K3), IMU = 1, NMUOBS)
         ENDDO
      ENDDO
      
      CLOSE (UNIT=2)
      RETURN
      END

c************************************************************************
      SUBROUTINE INTERP_SCAT_TABLE2 (WAVENO, DME,
     $                  EXTINCT, SSALB, ASYM,
     $              NDME, DMETAB, NWAVE, WAVETAB,
     $              TABEXTINCT, TABSSALB, TABASYM)
C       Interpolates the scattering properties from the table for
C     a particular wavenumber and particle size.  Does a bilinear 
C     interpolation, but optimized for fixed particle size and slowly 
C     varying wavenumber.  If the DME is the same as last time then we 
C     can just linearly interpolate in wavenumber between stored 
C     scattering values.  If the DME has changed then we linearly 
C     interpolate between the DMETAB grid lines for each of the two 
C     wavenumber grid lines.

      include 'scatter.param'

      REAL     WAVENO, DME
      REAL     EXTINCT, SSALB, ASYM
      INTEGER  NDME, NWAVE
      REAL     DMETAB(NDME), WAVETAB(NWAVE)
      REAL     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME) 
      REAL     TABASYM(NWAVE,NDME)
      INTEGER  IW0, IW1, ID, IL, IU, IM
      LOGICAL  NEWIW
      REAL     FWAV, FDME, FLDME, F
      REAL     OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
      SAVE     IW0, IW1, ID, OLDDME, FDME, FLDME
      SAVE     EXT0,EXT1, ALB0,ALB1, ASYM0,ASYM1
      DATA     IW0/1/, IW1/2/

C         Check that parameter are in range of table
      IF (WAVENO .LT. WAVETAB(1) .OR. WAVENO .GT. WAVETAB(NWAVE))
     .  STOP 'INTERP_SCAT_TABLE: wavenumber out of range.'
      IF (DME .LT. DMETAB(1) .OR. DME .GT. DMETAB(NDME))
     .  STOP 'INTERP_SCAT_TABLE: particle Dme out of range.'

C         See if wavenumber is within last wavenumber grid, otherwise
C           find the grid location and interpolation factor for WAVENO
      NEWIW = .FALSE.
      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (WAVENO .GE. WAVETAB(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE.
      ENDIF

      IF (DME .NE. OLDDME) THEN
C         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (DME .GE. DMETAB(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
      ENDIF

      IF (DME .NE. OLDDME .OR. NEWIW) THEN
C         If not the same Dme or a new wavenumber grid, then 
C           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) 
     $                + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) 
     $                + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)
      ENDIF

C         Linearly interpolate the scattering properties in wavenumber
      FWAV = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
      F = 1-FWAV
      EXTINCT = F*EXT0 + FWAV*EXT1
      SSALB = F*ALB0 + FWAV*ALB1
      ASYM = F*ASYM0 + FWAV*ASYM1

      OLDDME = DME
      RETURN
      END

c************************************************************************
      SUBROUTINE INTERP_SCAT_TABLE3 (MU, WAVENO, DME,
     $                  EXTINCT, SSALB, ASYM,
     $                  PHI1UP, PHI1DN, PHI2UP, PHI2DN,
     $              NMU, MUTAB, NDME, DMETAB, NWAVE, WAVETAB,
     $              TABEXTINCT, TABSSALB, TABASYM,
     $              TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
C       Interpolates the scattering properties from the table for
C     a particular observation angle, wavenumber, and particle size.
C     Does a trilinear interpolation, but optimized for fixed viewing
C     angle and particle size and slowly varying wavenumber.  If the
C     DME and MU are the same as last time then we can just linearly
C     interpolate in wavenumber between stored scattering values.
C     If the DME or MU has changed then we bilinearly interpolate between 
C     the DMETAB and MUTAB grid lines for each of the two wavenumber
C     grid lines.

      include 'scatter.param'
 
      REAL     MU, WAVENO, DME
      REAL     EXTINCT, SSALB, ASYM, PHI1UP, PHI1DN, PHI2UP, PHI2DN
      INTEGER  NMU, NDME, NWAVE
      REAL     MUTAB(NMU), DMETAB(NDME), WAVETAB(NWAVE)
      REAL     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME) 
      REAL     TABASYM(NWAVE,NDME)
      REAL     TABPHI1UP(NMU,NWAVE,NDME), TABPHI1DN(NMU,NWAVE,NDME)
      REAL     TABPHI2UP(NMU,NWAVE,NDME), TABPHI2DN(NMU,NWAVE,NDME)
      INTEGER  IW0, IW1, ID, IMU, IL, IU, IM
      LOGICAL  NEWIW
      REAL     FWAV, FDME, FLDME, FMU, F1, F2, F3, F4
      REAL     OLDMU, OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
      REAL     PHI1UP0, PHI1UP1, PHI1DN0, PHI1DN1
      REAL     PHI2UP0, PHI2UP1, PHI2DN0, PHI2DN1
      SAVE     IW0, IW1, IMU,ID, OLDMU, OLDDME, FDME, FLDME, FMU
      SAVE     EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
      SAVE     PHI1UP0, PHI1UP1, PHI1DN0, PHI1DN1
      SAVE     PHI2UP0, PHI2UP1, PHI2DN0, PHI2DN1
      DATA     IW0/1/, IW1/2/
      
C         Check that parameter are in range of table
      IF (WAVENO .LT. WAVETAB(1) .OR. 
     .  WAVENO .GT. WAVETAB(NWAVE)) THEN 
        PRINT *, WAVENO,WAVETAB(1),WAVETAB(NWAVE)
        STOP 'INTERP_SCAT_TABLE3: wavenumber out of range.'
      ENDIF

      IF (DME .LT. DMETAB(1) .OR. DME .GT. DMETAB(NDME)) THEN
        WRITE(*,*) 'Dme: ', DME, DMETAB(1), DMETAB(NDME)
        STOP 'INTERP_SCAT_TABLE: particle Dme out of range.'
      ENDIF
c      IF (MU .LT. MUTAB(1) .OR. MU .GT. MUTAB(NMU))
c     .  STOP 'INTERP_SCAT_TABLE3: viewing angle mu out of range.'

C         See if wavenumber is within last wavenumber grid, otherwise
C           find the grid location and interpolation factor for WAVENO
      NEWIW = .FALSE.
      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (WAVENO .GE. WAVETAB(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE.
      ENDIF

      IF (MU .NE. OLDMU) THEN
C         Find the grid location and interpolation factor for MU
        IL=1
        IU=NMU
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (MU .GE. MUTAB(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        IMU = MAX(IL,1)
        FMU = (MU-MUTAB(IMU))/(MUTAB(IMU+1)-MUTAB(IMU))
      ENDIF
      IF (DME .NE. OLDDME) THEN
C         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL .GT. 1)
          IM = (IU+IL)/2
          IF (DME .GE. DMETAB(IM)) THEN
            IL = IM
          ELSE
            IU = IM
          ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
      ENDIF

C         If not the same Dme and mu, then bilinearly interpolate things
C           Logarithmically interpolate extinction
      IF (DME .NE. OLDDME .OR. MU .NE. OLDMU .OR. NEWIW) THEN
        F1 = (1-FMU)*(1-FDME)
        F2 = (1-FMU)*FDME
        F3 = FMU*(1-FDME)
        F4 = FMU*FDME
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) 
     $                + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) 
     $                + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)

        PHI1UP0 = F1*TABPHI1UP(IMU,IW0,ID) + F2*TABPHI1UP(IMU,IW0,ID+1)
     $    + F3*TABPHI1UP(IMU+1,IW0,ID) + F4*TABPHI1UP(IMU+1,IW0,ID+1)
        PHI1UP1 = F1*TABPHI1UP(IMU,IW1,ID) + F2*TABPHI1UP(IMU,IW1,ID+1)
     $    + F3*TABPHI1UP(IMU+1,IW1,ID) + F4*TABPHI1UP(IMU+1,IW1,ID+1)

        PHI1DN0 = F1*TABPHI1DN(IMU,IW0,ID) + F2*TABPHI1DN(IMU,IW0,ID+1)
     $    + F3*TABPHI1DN(IMU+1,IW0,ID) + F4*TABPHI1DN(IMU+1,IW0,ID+1)
        PHI1DN1 = F1*TABPHI1DN(IMU,IW1,ID) + F2*TABPHI1DN(IMU,IW1,ID+1)
     $    + F3*TABPHI1DN(IMU+1,IW1,ID) + F4*TABPHI1DN(IMU+1,IW1,ID+1)

        PHI2UP0 = F1*TABPHI2UP(IMU,IW0,ID) + F2*TABPHI2UP(IMU,IW0,ID+1)
     $    + F3*TABPHI2UP(IMU+1,IW0,ID) + F4*TABPHI2UP(IMU+1,IW0,ID+1)
        PHI2UP1 = F1*TABPHI2UP(IMU,IW1,ID) + F2*TABPHI2UP(IMU,IW1,ID+1)
     $    + F3*TABPHI2UP(IMU+1,IW1,ID) + F4*TABPHI2UP(IMU+1,IW1,ID+1)

        PHI2DN0 = F1*TABPHI2DN(IMU,IW0,ID) + F2*TABPHI2DN(IMU,IW0,ID+1)
     $    + F3*TABPHI2DN(IMU+1,IW0,ID) + F4*TABPHI2DN(IMU+1,IW0,ID+1)
        PHI2DN1 = F1*TABPHI2DN(IMU,IW1,ID) + F2*TABPHI2DN(IMU,IW1,ID+1)
     $    + F3*TABPHI2DN(IMU+1,IW1,ID) + F4*TABPHI2DN(IMU+1,IW1,ID+1)
      ENDIF

C         Linearly interpolate the scattering properties in wavenumber
      FWAV = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
      F1 = 1-FWAV
      EXTINCT = F1*EXT0 + FWAV*EXT1
      SSALB = F1*ALB0 + FWAV*ALB1
      ASYM = F1*ASYM0 + FWAV*ASYM1
      PHI1UP = F1*PHI1UP0 + FWAV*PHI1UP1
      PHI1DN = F1*PHI1DN0 + FWAV*PHI1DN1
      PHI2UP = F1*PHI2UP0 + FWAV*PHI2UP1
      PHI2DN = F1*PHI2DN0 + FWAV*PHI2DN1

      OLDDME = DME
      OLDMU = MU
      RETURN
      END

c************************************************************************
      SUBROUTINE INTERPOLATE_PROFILE (N, XO, YO, XN, YN)
C     Interpolates Xo,Yo data to Xn finding Yn

      INTEGER N
      REAL*8    XO(*), XN
      REAL      YO(*), YN
      INTEGER I, IL, IU, IM
      REAL    X

c     Find XO value that is closest to XN
      
      IL=0  
      IU=N+1
10    IF (IU-IL.GT.1) THEN
        IM=(IU+IL)/2
        IF (XN.LT.XO(IM)) THEN
          IU=IM
        ELSE   
          IL=IM
        ENDIF 
      GO TO 10
      ENDIF
      I=IL

C     If XN out of XO range, don't extrapolate, use min/max XO
     
      IF (I .LE. 0) THEN
        I = 1  
        X = 0.0
      ELSE IF (I .GE. N) THEN
        I = N-1
        X = 1.0
      ELSE
        X = (XN-XO(I))/(XO(I+1)-XO(I))
      ENDIF

      YN = YO(I) + (YO(I+1)-YO(I))*X

      RETURN
      END
c************************************************************************
      SUBROUTINE READ_ABS_PROFILE (ABSFILE, BINARYABSFILE,
     $             ABSNU1, ABSNU2, ABSDELNU, 
     $             NABSNU,                       !!!!!!!MAXNZ, MAXABSNU, 
     $             NLEV, HEIGHT, TEMP, ABSPROF)
C       Reads in the gaseous absorption file which contains layer 
C     optical depths for NLEV-1 layers and NABSNU wavenumbers.
C     The height (km) and temperature (K) of the profile are also returned.
C     The file may be in ascii text or binary format.

      include 'scatter.param'

      INTEGER NABSNU, NLEV            !!!!!!!!!!!MAXNZ, MAXABSNU
      LOGICAL BINARYABSFILE
      REAL*8  ABSNU1, ABSNU2, ABSDELNU
      REAL    HEIGHT(*), TEMP(*), ABSPROF(MAXNZ,*)
      CHARACTER ABSFILE*(*)
      INTEGER I, J
      REAL    NU

      IF (BINARYABSFILE) THEN
        OPEN (UNIT=1,FILE=ABSFILE, STATUS='OLD',FORM='UNFORMATTED')
        READ (1) ABSNU1, ABSNU2, ABSDELNU, NABSNU
C        WRITE (*,*) ABSNU1, ABSNU2, ABSDELNU, NABSNU
        IF (NABSNU .GT. MAXABSNU) STOP 'RTSPEC: MAXABSNU exceeded'
        READ (1) NLEV
c        WRITE (*,*) NLEV
        IF (NLEV .GT. MAXNZ) STOP 'RTSPEC: MAXNZ exceeded'
        READ (1) (HEIGHT(I), I=1,NLEV)
c        WRITE (*,'(52(1X,F5.1),:)') (HEIGHT(I), I=1,NLEV)
        READ (1) (TEMP(I), I=1,NLEV)
c        WRITE (*,'(52(1X,F5.1),:)') (TEMP(I), I=1,NLEV)
        DO J = 1, NABSNU
          READ (1) (ABSPROF(I,J), I=1,NLEV-1)
        ENDDO
        CLOSE (1)
      ELSE
        OPEN (UNIT=1, FILE=ABSFILE, STATUS='OLD')
        READ (1,*)
        READ (1,*) ABSNU1, ABSNU2, ABSDELNU, NABSNU

        IF (NABSNU .GT. MAXABSNU) THEN
          PRINT *, 'NABSNU',NABSNU,MAXABSNU
          STOP 'RTSPEC: MAXABSNU exceeded'
        ENDIF
        READ (1,*) NLEV
        IF (NLEV .GT. MAXNZ) STOP 'RTSPEC: MAXNZ exceeded'
        READ (1,*) (HEIGHT(I), I=1,NLEV)
        READ (1,*) (TEMP(I), I=1,NLEV)
        READ (1,*)
        DO J = 1, NABSNU
          READ (1,*) NU, (ABSPROF(I,J), I=1,NLEV-1)
        ENDDO
        CLOSE (1)
      ENDIF
      RETURN
      END

c************************************************************************
c this subroutine does downward thermalrad tansfer from iS to iE
c ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
C THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
c and that raTemp has already been initialized with 2.96k Planck fcn

c this is the same as  SUBROUTINE FastBDRYL2GDiffusiveApprox_rtspec in 
c rad_diff.f, except that it has been modified for rtspec.f applications,
c as 1..iNumLayer ---> iNumLayer .. 1
c 
c for layers   20..1, it uses acos(3/5)
c for layers 100..20, it does t(i-1->0,x1)-t(i->0,x2) 
c    where x1 is calculated at layer i-1, x2 is calculated at layer i

      SUBROUTINE FastBDRYL2GDiffusiveApprox_rtspec(TOA_to_instr,
     $                MU, WAVENO, 
     $                NLEV, TEMP, TAU, RAD0DN, RAD0DN0, ibdry) 
C     Compute the downwelling radiance at the bottom of each cloud layer 
C     using the Fast Diffusivity Approx  
 
      include 'scatter.param' 
 
      INTEGER NLEV
      REAL    MU, WAVENO
      REAL    TEMP(NLEV), TAU(NLEV)
      REAL    RAD0DN(kProfLayer) 
      REAL    RAD0DN0, TOA_to_instr
      INTEGER ibdry
 
      INTEGER iS,iE

c local variables
      INTEGER iLay,iL,iLp1,iBdryP1,iSecondEnd,iCase
      REAL r1,r2,rPlanck,rMPTemp,rFreqAngle,rFreqAngle_m1

c to do the angular integration
      REAL rAngleTr_m1,rAngleTr,rL2G,rL2Gm1
      REAL FindDiffusiveAngleExp,rDiff,rCosDiff

      r1=kPlanck1
      r2=kPlanck2

      iBdryP1=NLev
      iCase=-1

c note that we are not as careful as  FastBDRYL2GDiffusiveApprox in that we
c do not completely fill up atmospehere to include layers above instrument
c (if intrument is not at TOA)    
      iS=nlev-1
      iE=1
c now we have 3 different cases to consider
c CASE A1 : easy -- this is do ENTIRE atmnosphere
c iS=100   iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
c                                     do iB->iE using accurate diffusive approx
c CASE A2 : easy -- this is do instr-gnd
c iS~50    iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
c                                     do iB->iE using accurate diffusive approx
       IF ((iS .GE. iBdry) .AND. (iBdry .GE. iE)) THEN
         iCase=1
         iBdryP1=iBdry+1
         END IF
c CASE B : quite easy -- this is do atmosphere -- instr
c iS=100   iE>iB                  ==> do iS->iE using acos(3/5)
       IF ((iS .GE. iBdry) .AND. (iBdry .LE. iE)) THEN
         iCase=2
         iBdryP1=iE
         END IF
c CASE C : easy -- this is do instr-gnd
c iS~50    iE~1   iB > iS,iE      ==> do iB->iE using accurate diffusive approx
       IF ((iBdry .GE. iS) .AND. (iBdry .GE. iE)) THEN
         iCase=3
         iBdry=iS
         END IF

      IF (iCase .EQ. -1) THEN
        write(kStdErr,*)'In FastBDRYL2GDiffusiveApprox_rtspec, icase = -1'
        CALL DoSTOP
        END IF

c ****** now map 1 .. iNumLayer ------> iNumlayer .. 1      *******
      iBdryP1=nlev-iBdryP1+1
      iBdry=nlev-iBdry+1
      iE=nlev-1
      iS=1
c ******** also note iLm1 ---> iLp1
      
      rDiff=(kThermalAngle*kPi/180.0)
      rCosDiff=cos(rDiff)

      if (TOA_to_instr .LT. 0) THEN
        RAD0DN0=KPLANCK1 *WAVENO**3 
     $            / (EXP(kPlanck2*WAVENO/kTSpace) - 1) 
      else 
        RAD0DN0=TOA_to_instr
        END IF        

c     now just go from TOA to instrument .. assume there are no clouds
c      RAD0DN0=RAD0DN0*exp(-TOA_to_instr)

      RAD0DN(1) = RAD0DN0

c initalize raL2G,raL2Gm1 
      rL2G=0.0
      rL2Gm1=0.0

c calculate rL2Gm1 which is the L2G transmission from layer iS-1 to ground
      DO iLay=2,nlev-1
        iL=iLay
        rL2Gm1=rL2Gm1+tau(iL)
        END DO
c calculate rL2G which is the L2G transmission from layer iS to ground
c and initialise the angles
      rL2G=rL2Gm1+tau(1)

c do top part of atmosphere, where we can use acos(3/5)
      IF ((iCase .EQ. 1)  .OR. (iCase. EQ. 2)) THEN
c go from top of atmosphere to boundary
        DO iLay=iS,iBdryp1
          iL=iLay
          iLp1=iLay+1
          rMPTemp=temp(iL)
c find the diffusive angles for the layer beneath
          rAngleTr_m1=exp(-rL2Gm1/rCosDiff)
          rAngleTr=exp(-rL2G/rCosDiff)
c Planckian emissions
          rPlanck=exp(r2*waveno/rMPTemp)-1.0
          rPlanck=r1*((waveno**3))/rPlanck
          RAD0DN0=RAD0DN0+rPlanck*(rAngleTr_m1-rAngleTr)
          RAD0DN(1) = RAD0DN0          
c get ready for the layer beneath
          rL2G=rL2Gm1
          rL2Gm1=rL2Gm1-tau(iLp1)
          END DO
        END IF

      IF ((iCase .EQ. 1) .OR. (iCase .EQ. 3)) THEN
c go from boundary to ground, or iE
c do bottom part of atmosphere ACCURATELY

        iSecondEnd=nlev-1
        rAngleTr=FindDiffusiveAngleExp(rL2G)
        rFreqAngle=rAngleTr

        DO iLay=iBdry,iSecondEnd
          iL=iLay
          iLp1=iLay+1
          rMPTemp=temp(iL)
c find the diffusive angles for the layer beneath
          rAngleTr_m1=FindDiffusiveAngleExp(rL2Gm1)
          rFreqAngle_m1=rAngleTr_m1
          rAngleTr_m1=exp(-rL2Gm1/cos(rAngleTr_m1))
          rAngleTr=rFreqAngle
          rAngleTr=exp(-rL2G/cos(rAngleTr))
c Planckian emissions
          rPlanck=exp(r2*WAVENO/rMPTemp)-1.0
          rPlanck=r1*((WAVENO**3))/rPlanck
          RAD0DN0=RAD0DN0+rPlanck*(rAngleTr_m1-rAngleTr)
          RAD0DN(1) = RAD0DN0          
c get ready for the layer beneath
          rL2G=rL2Gm1
          rL2Gm1=rL2Gm1-tau(iLp1)
          rFreqAngle=rFreqAngle_m1
          END DO

        END IF

c whether we did gaussian quadrature or diffusive approx, we now need the 2pi
c factor from the azimuthal integration
c however, there is also an average factor of 0.5 ==> overall, we need "pi"
      RAD0DN0=RAD0DN0*kPi

      RETURN
      END

c************************************************************************
c this subroutine checks to see if there are any layers above the instrument
c as they have to be added on to do the solar/backgnd thermal correctly!! 
c same as AddUppermostLayersQ, except it accepts raaAbs as input, and 
c outputs radiance from TOA to instr ---- if instr is at TOA, it outputs -10

      SUBROUTINE  Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,rFracTop, 
     $                                raWaves,raaAbs,raExtra) 
 
      include 'kcarta.param' 
 
c rFracTop tells how much of the upper layer has been used, due to instr posn  
c iaRadLayer = current radiating atmosphere defn : gnd to instrument 
c iNumLayers = number of mixed paths in the defined radiating atmosphere 
c iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument, 
c                  temporarily define atm from GND to TOP of atmosphere 
c iT             = number of layers in this temporary atmosphere 
c iExtra = -1 if no layeres added on, +1 if layers added on 
c raExtra = array initialized to all zeros if instr at TOA
c         = array initialized to sum(k) from TOA to instr if instr inside atm
      INTEGER iNumLayer,iaRadLayer(kProfLayer) 
      REAL raExtra(kMaxPts),rFracTop,raaAbs(kMaxPts,kMixFilRows)
      REAL raVTemp(kMixFilRows),raWaves(kMaxPts)
 
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iExtra 
      INTEGER iI,iFr,iJ

      REAL waveno,rad,k,mudown
 
      iExtra=-1 
 
c check to see the posn of the instrument (defined by layers i1,i2,..iN),  
c relative to physical top of atmosphere, as defined by 100 layers 
      iI=MOD(iaRadLayer(iNumLayer),kProfLayer) 
c if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know 
c that ALL upper layers have been used in the atmosphere defn. 
cwe DO have to check that even if topmaost layer=100, it could still be  
c fractionally weighted due to the posn of instr at top layer being within 
c the layer, not on top of it 

      DO iFr=1,kMaxPts
        raExtra(iFr)=-10.0 
        END DO 
 
      IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .LE. 1.0e-4))THEN 
c current defined atmosphere has all g-100 layers, 100th layer had frac 1.0 
        iExtra=-1 
 
      ELSE IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .GE. 1.0e-4))THEN 
c even though the current defined atmosphere has all g-100 layers,  
c 100th layer had frac 0 < f < 1 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
c        write(kStdWarn,*) 'top most layer is fractional layer. Some' 
c        write(kStdWarn,*) 'portion needed above instrument to calculate' 
c        write(kStdWarn,*) ' thermal/solar' 
 
      ELSE IF ((iI .NE. 0)) THEN 
c current defined atmosphere does not have all g-100 layers 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
c now add on the upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer)=0 
 15     CONTINUE 
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) .NE. 0) THEN 
          iT=iT+1 
          iaRadLayerTemp(iT)=iaRadLayerTemp(iT-1)+1 
c          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT) 
          GO TO 15 
          END IF 
c        write(kStdWarn,*)'added ',iT-iNumLayer,' layers' 
c        write(kStdWarn,*)'above instrument to calculate th/solar/flux' 
        END IF 


ccccccccccccc this is new .. where subroutine differs from AddUpperMostLayers
      if (iExtra .gt. 0) THEN
        MUDOWN=3.0/5.0
        DO iFr=1,kMaxPts
          waveno=raWaves(iFr)
          raExtra(iFr) = KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/kTSpace) - 1)
          END DO
        DO iI=iT,iNumLayer+1,-1
          iJ=iaRadLayerTemp(iI)
          DO iFr=1,kMaxPts
            waveno=raWaves(iFr)
            k=raaAbs(iFr,iJ)
            rad=KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/raVTemp(iJ)) - 1)
            raExtra(iFr)=raExtra(iFr)*exp(-k/MUDOWN)+rad*(1-exp(-k/MUDOWN))
            END DO
          END DO

        DO iI=iNumLayer,iNumLayer
          iJ=iaRadLayerTemp(iI)
          DO iFr=1,kMaxPts
            waveno=raWaves(iFr)
            k=raaAbs(iFr,iJ)*(1-rFracTop)
            rad=KPLANCK1 *WAVENO**3
     $            / (EXP(kPlanck2*WAVENO/raVTemp(iJ)) - 1)
            raExtra(iFr)=raExtra(iFr)*exp(-k/MUDOWN)+rad*(1-exp(-k/MUDOWN))
            END DO
          END DO

        END IF

      RETURN 
      END 

c************************************************************************
c this file reads a binary made from the ASCII sscatmie.x file
      SUBROUTINE READ_SSCATTAB_BINARY(SCATFILE,   !!!  MAXTAB, MAXGRID,
     $          NMUOBS, MUTAB, NDME, DMETAB, NWAVE, WAVETAB,
     $          MUINC, TABEXTINCT, TABSSALB, TABASYM,
     $          TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

      include 'scatter.param'

C       Reads in the single scattering table for a number of wavenumbers,
C     particle sizes, and viewing angles.  The scattering properties are
C     computed for a IWC/LWC of 1 g/m^3.  
C       Input parameters:
C     SCATFILE   file name of scattering file
C     MAXTAB     maximum array size for whole table
C     MAXGRID    maximum array size for table grids
C       Output parameters:
C     NMUOBS     number of viewing angle mu grid values
C     MUTAB      viewing angle grid values
C     NDME       number of Dme grid values
C     DMETAB     Dme grid values
C     NWAVE      number of wavenumber grid values
C     WAVETAB    wavenumber grid values
C     MUINC(2)   cosine zenith of two incident angles
C     TABEXTINCT tabulated extinction (km^-1)
C     TABSSALB   tabulated single scattering albedo 
C     TABASYM    tabulated asymmetry parameter
C     TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN  
C                tabulated phase function terms for incident radiance angles

c!!!      INTEGER  MAXTAB, MAXGRID
      INTEGER  NMUOBS, NDME, NWAVE
      REAL     MUTAB(*), DMETAB(*), WAVETAB(*)
      REAL     MUINC(2)
      REAL     TABEXTINCT(*), TABSSALB(*), TABASYM(*) 
      REAL     TABPHI1UP(*), TABPHI1DN(*)
      REAL     TABPHI2UP(*), TABPHI2DN(*)
      CHARACTER*(*) SCATFILE
      INTEGER  IMU, ID, IW, K2, K3
      INTEGER IERR

      OPEN (UNIT=kTempUnit, STATUS='OLD', FORM='UNFORMATTED', 
     $      FILE=SCATFILE, IOSTAT=IERR)
       IF (IERR .NE. 0) THEN 
          WRITE(kStdErr,1010) IERR, SCATFILE
          CALL DoSTOP 
       ENDIF 
 1010  FORMAT('ERROR! number ',I5,' opening scatter data file:',/,A80) 

      kTempUnitOpen=1
      READ(kTempUnit) NMUOBS
      READ(kTempUnit) NDME
      READ(kTempUnit) NWAVE
      IF (MAX(NMUOBS,NDME,NWAVE) .GT. MAXGRID) 
     $    STOP 'READ_SSCATTAB_BINARY: MAXGRID exceeded'
      IF (NMUOBS*NDME*NWAVE .GT. MAXTAB) 
     $    STOP 'READ_SSCATTAB_BINARY: MAXTAB exceeded'
      READ(kTempUnit) MUINC(1), MUINC(2)
      READ(kTempUnit) (MUTAB(IMU), IMU = 1, NMUOBS)         

      DO IW = 1, NWAVE
        DO ID = 1, NDME
            K2 = IW-1 + NWAVE*(ID-1)
            K3 = NMUOBS*K2
            READ(kTempUnit) DMETAB(ID), WAVETAB(IW), TABEXTINCT(K2+1),
     $       TABSSALB(K2+1), TABASYM(K2+1)
            READ(kTempUnit) (TABPHI1UP(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI2UP(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI1DN(IMU+K3), IMU = 1, NMUOBS)
            READ(kTempUnit) (TABPHI2DN(IMU+K3), IMU = 1, NMUOBS)
         ENDDO
      ENDDO
      
      CLOSE (kTempUnit)
      kTempUnitOpen=-1
      write(kStdWarn,*)'sucess : read in binary scattr data from file = '
      write(kStdWarn,1020) scatfile

 1020 FORMAT(A80)

      RETURN
      END

c************************************************************************
