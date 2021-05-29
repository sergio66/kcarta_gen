! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_rtspec_code

USE basic_common
USE ttorad_common
USE jac_main
!USE jac_pclsam_up
!USE jac_pclsam_down
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_diff_and_quad
USE spline_and_sort_and_common

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the rtspec routines  **********************
! ***********   this is the LATEST rtspec version  **********************

!     Author: Frank Evans  University of Colorado   Feb 2001
!     Co-authors:
!       Aaron Evans    (multilayer code, spectral segments, double sidebands,
!                       upward looking geometry, output convolution)
!       Merritt Deeter (single scattering layer routines)

! modifications by Sergio De Souza-Machado
! 1) all subroutines now call include '../INCLUDE/TempF90/scatterparam.f90'
! 2) gasRT1 now has surface emissivity passed in, so it can be used
! 3) nlev ---> nlev - 1 whenever we are using optical depths
!    ie there are nlev lavels and so there are nlev-1 layers.
! 4) background thermal has been included!!!!!!!!!!!!!!!!!!!!!!
!    when called, radobs=0.0 if there is no backgnd thermal
!                 radobs=downward thermal radiance at posn of instrument
!                        simply computed using acos(3/5)
! 5) for DISORT, we need to unscale the parameters alb,extinct, asym
!    to do this, we need to know information on how they were scaled by
!    sscatmie. read_scattab modified to read in and store this info

! do not "data save" iw0,iw1 in INTERP_SCAT_TABLE2,3 as there are many times
! this subroutine is called, so it is better to start always with 1,2
!************************************************************************
!************************************************************************
! these are the new fixed Feb 2001 routines for use in kCARTA
    SUBROUTINE COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
    MUOBS, IOBS, WAVENO, &
    NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
!!!! .               MAXTAB, MAXGRID,
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
    RADOBS, &
    ibdry, TOA_to_instr,iiDiv)  !3 new parameters by SSM
     
!       Calculates the radiance at the desired level (IOBS) and angle (MUOBS)
!     for one wavenumber (WAVENO) for the atmosphere and cloud specified.
!     The temperature and molecular absorption profiles are input
!     (TEMP, TAUGAS) for the NLEV levels.  The multilayer cloud is
!     between levels ICLDTOP and ICLDBOT. Each of the NCLDLAY cloud layers
!     has an ice water path (IWP), median particle diameter (DME), and
!     pointer to the scattering table (ISCATTAB).  The scattering properties
!     are input as tables. Accepts mu > 0 (down looking) or mu < 0 (up looking).
!     The observed radiance (RADOBS) is output in W/(m^2 ster cm^-1) units.
!     RTMODEL : S,E,H for single scatter, Eddington, and Hybrid respectively
!     MUOBS : Cosine of observation angle (negative for looking up)
!     IOBS : The level number of observation height
!     WAVENO : Wave number (1/cm) for this monochromatic RT calculation.
!     NLEV : Number of levels in atmosphere
!       Note: Levels/layers numbering start at top of atmopshere
!     TEMP : Temperature array for levels in atmosphere
!     TAUGAS : Optical depth array for gases in layer of atmosphere
!     SFCTEMP : Temperature of surface
!     SFCEMIS : Emissivity of surface
!     NCLDLAY : Number of cloud layers, ie. # of atm layer occupied by cloud
!     ICLDTOP : The level number of cloud top
!     ICLDBOT : The level number of cloud bottom
!     IWP : Ice Water Path (g/m^2) array for cloud layers
!     DME : Median Particle Diameter (um) array for cloud layers
!     ISCATTAB : Scattering table file number array for cloud layers
!       Note: ISCATAB allows multiple scattering table files
!     MAXTAB : Maximum number of elements in scattering table
!     MAXGRID : Maximum number of points on scattring table axis
!     NSCATTAB : Number of scattering tables/files
!     MUINC : Cosines of two incident angles of single scattering solution
!     NMUOBS : Number of viewing cosines in scattering table
!     MUTAB : Viewing cosines making up grid of scattering table
!     NDME : Number of Median Particle Diameters in scattering table
!     DMETAB : Median Particle Diameters making up grid of scattering table
!     NWAVETAB : Number of wavenumbers in scattering table
!     WAVETAB : Wavenumbers making up grid of scattering table
!     TABEXTINCT : Extinction (1/km) matrix, dimensions are wavenumber/Dme
!       Note : Scattering table calculated for IWC = 1 g/m^2, thus
!              Extinction = TABEXTINCT*IWP/(layer thickness)
!     TABSSALB : Single Scattering albedo matrix, dimensions are Wavenumber/Dme
!     TABASYM : Asymmetry Parameter matrix, dimensions are Wavenumber/Dme
!     TABPHI : "Phi" function matrix, dimensions are Wavenumber/Dme/Mu
!              "Phi" represents the amount of incident radiation scattered.
!     TABPHI1UP : "Phi" up applied to forward scattering radiation, for MU1
!     TABPHI2UP : "Phi" up applied to forward scattering radiation, for MU2
!     TABPHI1DN : "Phi" down applied to back scattering radiation, for MU1
!     TABPHI2DN : "Phi" down applied to back scattering radiation, for MU2
!     RADOBS : Output of subroutine, Upwelling radiation at observation height

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: NLEV, IOBS,iiDiv
    INTEGER :: NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::    MUOBS, WAVENO, SFCTEMP, SFCEMIS
    REAL ::    TEMP(KPROFLAYER+1), TAUGAS(KPROFLAYER)
    REAL ::    IWP(MAXNZ), DME(MAXNZ), IWPTOT
    REAL ::    RADOBS
    CHARACTER(1) :: RTMODEL
    INTEGER ::  NSCATTAB                     !!!!!!!!!!!!!! MAXTAB, MAXGRID
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
    INTEGER ::   IBDRY                          !!!for diffusive approx
    REAL ::      TOA_to_instr                   !!!sum(k) from TOA to instr

!         Local variables
!      INTEGER  MAXNZ
!      PARAMETER (MAXNZ=150)
! all these were changed from MAXNZ to KPROFLAYER
    INTEGER :: I, L, N
    REAL ::    TAUTOT(kProfLayer), COALB(kProfLayer), EXTINCT
    REAL ::    TAUC(kProfLayer), TAUCG(kProfLayer), SSALB(kProfLayer), &
    ASYM(kProfLayer)
    REAL ::    PHI1UP(kProfLayer), PHI1DN(kProfLayer)
    REAL ::    PHI2UP(kProfLayer), PHI2DN(kProfLayer)
    REAL ::    RADBNDRYUP(2), RADBNDRYDN(2)
    REAL ::    RADBNDRYUP1(kProfLayer), RADBNDRYDN1(kProfLayer)
    REAL ::    RADBNDRYUP2(kProfLayer), RADBNDRYDN2(kProfLayer)
    REAL ::    RAD0UPOBS(kProfLayer), RAD0DNOBS(kProfLayer)
    REAL ::    FLUXES(3,kProfLayer)
    REAL ::    FLUXTOP, FLUXBOT, TTOP, TBOT, RADBOTTOM, RADTOP
    REAL ::    FLUXUP, FLUXDN, FLUXUPSEDD, FLUXDNSEDD
    REAL :: rDummy1,rDummy2
    INTEGER :: Nprime,iDummy

!         Initialize the optical depths at each level
    IF (NLEV > MAXNZ) &
    STOP 'COMPUTE_RADIATIVE_TRANSFER: MAXNZ exceeded'
! this was originally       DO N = 1, NLEV ... wierd!!!!!!!!!!
    DO N = 1, NLEV-1
        TAUTOT(N) = TAUGAS(N)
        COALB(N) = 1.0
    ENDDO

    IWPTOT = 0
    DO N = 1,NCLDLAY
        IWPTOT = IWPTOT + IWP(N)
    ENDDO

! this is heavily modified!!!!!!!!!!!!!!!!!!!!!
! **************** orig version **************************
! cccC          If clear sky only do GASRT1 and GASRT2
! ccc      IF (IWPTOT.EQ.0) THEN
! cccC            Integrate radiance from surface to observer
! cccC           First go from surface to cloud bottom
! ccc        CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS,
! ccc     .        NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT,
! ccc     .                   NCLDLAY, RAD0UPOBS, RAD0DNOBS)
! ccc        IF(MUOBS.GT.0) THEN
! cccC           if looking down then go from cloud top to observer
! ccc          RADOBS = RAD0UPOBS(1)
! ccc          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
! ccc     .                ICLDTOP+1, IOBS, MUOBS, RADOBS)
! ccc        ELSEIF(MUOBS.LT.0) THEN
! cccC           Then go from cloud bottom to observer
! ccc          RADOBS = RAD0DNOBS(NCLDLAY)
! ccc          CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS,
! ccc     .                ICLDBOT-1, IOBS, ABS(MUOBS), RADOBS)
! ccc        ENDIF

!     If clear sky only do GASRT1 and GASRT2
! cccccc this was orig rtspec code       IF (IWPTOT.EQ.0) THEN
    IF ((IWP(1) <= 0)) THEN
    !       Integrate radiance from surface to observer
    !       First go from surface to cloud bottom
        IF (MUOBS >= 0) THEN
        ! f looking down, do the nocloud downlook rad trans
            CALL GASRT1_nocloud (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS, &
            NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT, &
            NCLDLAY, RAD0UPOBS, RAD0DNOBS, &
            ibdry, TOA_to_instr, iiDiv)
        !          CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS,
        !     $        NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT,
        !     $                   NCLDLAY, RAD0UPOBS, RAD0DNOBS, iiDiv)
        ELSE
        ! his was original Evans code
            CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS, &
            NLEV, TEMP, TAUTOT, COALB, ICLDTOP, ICLDBOT, &
            NCLDLAY, RAD0UPOBS, RAD0DNOBS, iiDiv)
        END IF

        IF(MUOBS >= 0) THEN
        !         if looking down then go from cloud top to observer
            RADOBS = RAD0UPOBS(1)
            RADOBS = RAD0UPOBS(NCLDLAY) !!!new
            CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS, &
            ICLDTOP+1, IOBS, MUOBS, RADOBS, iiDiv)
        !          print *,WAVENO, SFCTEMP, SFCEMIS,NLEV,RAD0UPOBS(NCLDLAY),RADOBS
        ELSEIF(MUOBS < 0) THEN
        !         Then go from cloud bottom to observer
            RADOBS = RAD0DNOBS(NCLDLAY)
            rDummy1=RADOBS
            CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS, &
            ICLDBOT-1, IOBS, ABS(MUOBS), RADOBS, iiDiv)
            rDummy2=RADOBS
        ENDIF

    ELSE
    !       cloudy sky case
    !          If cloud present do Single-scatter,Eddington or hybrid routine
    !              Get the optical properties for the cloud layers
        DO N = ICLDTOP, ICLDBOT - 1
            Nprime = N - iiDiv*kProfLayer
            L = N-ICLDTOP+1
        !             Interpolate to get values of extinction, s.s. albedo, and
        !             phi function values at given obs. mu, waveno, and particle size.
        !             Note: we don't need phi functions for Eddington only solution.
            I = ISCATTAB(L)
            IF (RTMODEL == 'E') THEN
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
            ELSE
                CALL INTERP_SCAT_TABLE3 (ABS(MUOBS), WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM(L), &
                PHI1UP(L), PHI1DN(L), PHI2UP(L), PHI2DN(L), &
                NMUOBS(I), MUTAB(1,I), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                TABPHI1UP(1,I), TABPHI1DN(1,I), &
                TABPHI2UP(1,I), TABPHI2DN(1,I))
            ENDIF

        !             Compute the optical depth of cloud layer, including gas
            TAUC(L) = IWP(L)*EXTINCT/1000.
        !!!!!!!!!!!! orig code TAUCG(L) = TAUGAS(N)      + TAUC(L)
        !!!!!!!!!!!! new  code TAUCG(L) = TAUGAS(Nprime) + TAUC(L)
            TAUCG(L) = TAUGAS(Nprime) + TAUC(L)
            IF (TAUCG(L) > 1.0E-5) THEN
                SSALB(L) = SSALB(L)*TAUC(L)/TAUCG(L)
            ELSE
                SSALB(L) = 0.0
            ENDIF
        !!!!!!!! orig code  TAUTOT(N) = TAUCG(L), COALB(N) = 1.0 - SSALB(L)
            TAUTOT(Nprime) = TAUCG(L)
            COALB(Nprime)  = 1.0 - SSALB(L)
        ENDDO

    !              Calculate the gas/cloud emission only radiative transfer
    !                for incident radiation on each layer
        CALL GASRT1 (MUINC(1), WAVENO, SFCTEMP, SFCEMIS, &
        NLEV, TEMP, TAUTOT, COALB, &
        ICLDTOP, ICLDBOT, NCLDLAY, &
        RADBNDRYUP1, RADBNDRYDN1, iiDiv)
        CALL GASRT1 (MUINC(2), WAVENO, SFCTEMP, SFCEMIS, &
        NLEV, TEMP, TAUTOT, COALB, &
        ICLDTOP, ICLDBOT, NCLDLAY, &
        RADBNDRYUP2, RADBNDRYDN2, iiDiv)
        CALL GASRT1 (ABS(MUOBS), WAVENO, SFCTEMP, SFCEMIS, &
        NLEV, TEMP, TAUTOT, COALB, &
        ICLDTOP, ICLDBOT, NCLDLAY, &
        RAD0UPOBS, RAD0DNOBS, iiDiv)
                      
    !        print *,' '
    !        print *,SFCTEMP,SFCEMIS,NLEV,ICLDTOP,ICLDBOT,NCLDLAY,waveno,MUOBS
    !        DO iDummy = 1,NLEV-2
    !          print *,iDummy,NCLDLAY,RAD0UPOBS(iDummy),RAD0DNOBS(iDummy),
    !     $             tautot(NLEV-1-idummy),temp(NLEV-1-iDummy)
    !          END DO

    !              Calculate the fluxes for the multilayer boundary conditions
        FLUXTOP = MUINC(1)*.5*RADBNDRYDN1(1) &
        + MUINC(2)*.5*RADBNDRYDN2(1)
        FLUXBOT = MUINC(1)*.5*RADBNDRYUP1(NCLDLAY) &
        + MUINC(2)*.5*RADBNDRYUP2(NCLDLAY)

    !              Don't compute multilayer fluxes for a single cloud layer
        IF (NCLDLAY <= 1) THEN
            FLUXES(2,1) = FLUXTOP
            FLUXES(1,2) = FLUXBOT
        ELSE
        !              For multilayer cloud calculate the multilayer Eddington fluxes
            CALL EDDRTF (NCLDLAY, TAUCG, SSALB, ASYM, &
            TEMP(ICLDTOP-iiDiv*kProfLayer), WAVENO, &
            FLUXTOP, FLUXBOT,  FLUXES)
        ENDIF

    !        If looking down must loop over cloud layers from bottom to top
        IF(MUOBS >= 0) THEN
        !           Introduce scalar variable for radiance at bottom of cloud
            RADBOTTOM = RAD0UPOBS(NCLDLAY)
                      
        !              Loop over the layers of cloud, from bottom to top

            DO L = NCLDLAY,1, -1
                TTOP = TEMP(L-1 + ICLDTOP - iiDiv*kProfLayer)
                TBOT = TEMP(L   + ICLDTOP - iiDiv*kProfLayer)
                FLUXUP = FLUXES(1,L+1)
                FLUXDN = FLUXES(2,L)
                  
            !                Do Single-scattering, Eddington, or Hybrid for each layer
                IF (RTMODEL == 'S') THEN
                    RADBNDRYUP(1) = RADBNDRYUP1(L)
                    RADBNDRYUP(2) = RADBNDRYUP2(L)
                    RADBNDRYDN(1) = RADBNDRYDN1(L)
                    RADBNDRYDN(2) = RADBNDRYDN2(L)
                    CALL SSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC, &
                    RADBNDRYUP, RADBNDRYDN, RADBOTTOM, &
                    TAUCG(L), SSALB(L), ASYM(L), &
                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L), &
                    RADOBS)
                ELSE IF (RTMODEL == 'E') THEN
                    CALL EDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, &
                    FLUXUP, FLUXDN, RADBOTTOM, &
                    TAUCG(L), SSALB(L), ASYM(L), RADOBS)
                ELSE
                    FLUXDNSEDD = MUINC(1)*.5*RADBNDRYDN1(L) &
                    + MUINC(2)*.5*RADBNDRYDN2(L)
                    FLUXUPSEDD = MUINC(1)*.5*RADBNDRYUP1(L) &
                    + MUINC(2)*.5*RADBNDRYUP2(L)
                    RADBNDRYUP(1) = RADBNDRYUP1(L)
                    RADBNDRYUP(2) = RADBNDRYUP2(L)
                    RADBNDRYDN(1) = RADBNDRYDN1(L)
                    RADBNDRYDN(2) = RADBNDRYDN2(L)
                    CALL SEDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC, &
                    RADBNDRYUP, RADBNDRYDN, &
                    FLUXUPSEDD, FLUXDNSEDD, &
                    FLUXUP, FLUXDN, RADBOTTOM, &
                    TAUCG(L), SSALB(L), ASYM(L), &
                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L), &
                    RADOBS)
                ENDIF

            !             The calculated upwelling radiance at top of layer is used
            !               as the incident radiance on bottom of next layer.
            !            print *,'cld ---> ',RTMODEL
            !            print *,'cld ---> ',L,TTOP,TBOT,RADBOTTOM,TAUCG(L),SSALB(L),RADOBS
                RADBOTTOM = RADOBS
            ENDDO

        !              Integrate radiance from cloud top to observer
            CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS, &
            ICLDTOP, IOBS, MUOBS, RADOBS, iiDiv)
        !         print *,'final ',RADOBS
        !         call dostop

        !        If looking up must loop over cloud layers from top to bottom
        ELSEIF(MUOBS < 0) THEN
        !           Introduce scalar variable for radiance at top of cloud
            RADTOP = RAD0DNOBS(1)

        !              Loop over the layers of cloud, from top to bottom
            DO L = 1,NCLDLAY
                TTOP = TEMP(L-1+ICLDTOP-iiDiv*kProfLayer)
                TBOT = TEMP(L+ICLDTOP-iiDiv*kProfLayer)
                FLUXUP = FLUXES(1,L+1)
                FLUXDN = FLUXES(2,L)
                  
            !                Do Single-scattering, Eddington, or Hybrid for each layer
                IF (RTMODEL == 'S') THEN
                    RADBNDRYUP(1) = RADBNDRYUP1(L)
                    RADBNDRYUP(2) = RADBNDRYUP2(L)
                    RADBNDRYDN(1) = RADBNDRYDN1(L)
                    RADBNDRYDN(2) = RADBNDRYDN2(L)
                    CALL SSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, MUINC, &
                    RADBNDRYDN, RADBNDRYUP, RADTOP, &
                    TAUCG(L), SSALB(L), ASYM(L), &
                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L), &
                    RADOBS)
                ELSE IF (RTMODEL == 'E') THEN
                    CALL EDDSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, &
                    FLUXDN, FLUXUP, RADTOP, &
                    TAUCG(L), SSALB(L), ASYM(L), RADOBS)
                ELSE
                    FLUXDNSEDD = MUINC(1)*.5*RADBNDRYDN1(L) &
                    + MUINC(2)*.5*RADBNDRYDN2(L)
                    FLUXUPSEDD = MUINC(1)*.5*RADBNDRYUP1(L) &
                    + MUINC(2)*.5*RADBNDRYUP2(L)
                    RADBNDRYUP(1) = RADBNDRYUP1(L)
                    RADBNDRYUP(2) = RADBNDRYUP2(L)
                    RADBNDRYDN(1) = RADBNDRYDN1(L)
                    RADBNDRYDN(2) = RADBNDRYDN2(L)
                    CALL SEDDSCATRT (ABS(MUOBS), TBOT, TTOP, WAVENO, MUINC, &
                    RADBNDRYDN, RADBNDRYUP, &
                    FLUXDNSEDD, FLUXUPSEDD, &
                    FLUXDN, FLUXUP, RADTOP, &
                    TAUCG(L), SSALB(L), ASYM(L), &
                    PHI1UP(L),PHI1DN(L),PHI2UP(L),PHI2DN(L), &
                    RADOBS)
                ENDIF

            !             The calculated upwelling radiance at bottom of layer is used
            !               as the incident radiance on top of next layer.
                RADTOP = RADOBS
            ENDDO
        !              Integrate radiance from cloud top to observer
            CALL  GASRT2 (WAVENO, NLEV, TEMP, TAUGAS, &
            ICLDBOT, IOBS, ABS(MUOBS), RADOBS, iiDiv)
        ENDIF
              
    ENDIF

    RETURN
    END SUBROUTINE COMPUTE_RADIATIVE_TRANSFER

!************************************************************************
    SUBROUTINE SSCATRT (MUOBS, TTOP, TBOT, WAVENO, MUINC, &
    IBNDRYUP, IBNDRYDN, I0UPOBS, &
    TAUTOT, SSALB, ASYM, &
    PHI1UP, PHI1DN, PHI2UP, PHI2DN,  IUPOBS)

!     Calculates upwelling radiance from a single thermally emitting
!     homogeneous layer using the single scattering radiative transfer
!     method of Deeter and Evans (1998).
!    Input Parameters:
!       MUOBS      observation mu (cosine of viewing zenith angle)
!       TTOP       cloud top temperature (K)
!       TBOT       cloud bottom temperature (K)
!       WAVENO     wavenumber (cm^-1)
!       MUINC      cosine zenith angles for the two incident directions
!       IBNDRYUP   two incident radiances upwelling on bottom of layer
!       IBNDRYDN   two incident radiance downwelling on top of layer
!       I0UPOBS    incident radiance at bottom of layer at observation angle
!       TAUTOT     optical depth of layer
!       SSALB      single scattering albedo of layer
!       ASYM       asymmetry parameter of layer
!       PHI*       single scattering functions at the observation angle

!    Output Parameters:
!       IUPOBS     outgoing radiance at top at observation angle

!      Input and output radiances are in units of W m^-2 sr^-1 cm and
!    fluxes in units of W m^-2 cm.

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
         
    REAL ::    MUOBS, TTOP, TBOT, WAVENO, MUINC(2)
    REAL ::    IBNDRYUP(2), IBNDRYDN(2), I0UPOBS, IUPOBS
    REAL ::    TAUTOT, SSALB, ASYM
    REAL ::    PHI1UP, PHI1DN, PHI2UP, PHI2DN

    INTEGER :: IMU
    REAL ::    BETA, COPLNCK, PLANCK0, PLANCK1
    REAL ::    EXPBETA, EXPMUOBS, DMUOBS, EXPMAX
    REAL ::    MU, DMUUP, DMUDN, EXPMU
    REAL ::    S0UP, S1UP, IUP(2), IDN(2)
    REAL ::    TERM1, TERM2, TERM3
          
    EXPMAX = 1.E08

    IF (TAUTOT < 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
    ENDIF

!           Planck function radiances in units of W m^-2 sr^-1 cm
    PLANCK0 = ttorad(WAVENO,TTOP)
    PLANCK1 = ttorad(WAVENO,TBOT)
    BETA = LOG(PLANCK1/PLANCK0)/TAUTOT

    COPLNCK = (1. - SSALB)*PLANCK0
    EXPBETA = PLANCK1/PLANCK0
    EXPMUOBS = EXP(-TAUTOT/MUOBS)
    IF (EXPMUOBS < 1./EXPMAX) EXPMUOBS = 1./EXPMAX
    DMUOBS = 1. - BETA*MUOBS
    IF(ABS(DMUOBS) < 1.E-4) DMUOBS = SIGN(1.E-4,DMUOBS)

!c the "exponential in radiance" variation is better than "exponential in T"
!c      print *,'beta,expbeta 0 = ',beta,expbeta
!c      BETA = LOG(TTOP/TBOT)/TAUTOT             !!!!!!
!c      EXPBETA = TTOP/TBOT                      !!!!!!
!c      print *,'beta,expbeta 1 = ',beta,expbeta

!           Zero order scattering contribution
    S0UP = (1. - EXPBETA*EXPMUOBS)*COPLNCK/DMUOBS

!           IMU (1 or 2) is index for discrete mu's in simple model
    DO IMU = 1, 2
        MU = MUINC(IMU)
        EXPMU = EXP(-TAUTOT/MU)
        IF (EXPMU < 1./EXPMAX) EXPMU = 1./EXPMAX
        DMUUP = 1. - BETA*MU
        DMUDN = 1. + BETA*MU
        IF(ABS(DMUUP) < 1.E-4) DMUUP = SIGN(1.E-4,DMUUP)
        IF(ABS(DMUDN) < 1.E-4) DMUDN = SIGN(1.E-4,DMUDN)

        TERM1 = IBNDRYUP(IMU)*EXPMU - COPLNCK*EXPBETA*EXPMU/DMUUP
        TERM2 = (MU/(MUOBS - MU))*(1. - EXPMUOBS/EXPMU)
        TERM3 = COPLNCK*(1. - EXPBETA*EXPMUOBS)/(DMUUP*DMUOBS)
        IUP(IMU) = - (TERM1*TERM2 - TERM3)

        TERM1 = IBNDRYDN(IMU) - COPLNCK/DMUDN
        TERM2 = (MU/(MUOBS + MU))*(1. - EXPMUOBS*EXPMU)
        TERM3 =  COPLNCK*(1. - EXPBETA*EXPMUOBS)/(DMUDN*DMUOBS)
        IDN(IMU) = (TERM1*TERM2 + TERM3)
    ENDDO

!           First order scattering contribution
    S1UP = (SSALB/2.)*(PHI1UP*IUP(1) + PHI2UP*IUP(2) &
    + PHI1DN*IDN(1) + PHI2DN*IDN(2))

    IUPOBS = I0UPOBS*EXPMUOBS + S0UP + S1UP

    RETURN
    END SUBROUTINE SSCATRT

!************************************************************************
    SUBROUTINE EDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, &
    FLUXUP, FLUXDN, I0UPOBS, &
    TAUTOT, SSALB, ASYM, IUPOBS)
!     Calculates upwelling radiance from a single thermally emitting
!     homogeneous layer using Eddington's second approximation
!     radiative transfer method.

!    Input Parameters:
!       MUOBS      observation mu (cosine of viewing zenith angle)
!       TTOP       cloud top temperature (K)
!       TBOT       cloud bottom temperature (K)
!       WAVENO     wavenumber (cm^-1)
!       FLUXUP     incident flux on bottom from multilayer Eddington solution
!       FLUXDN     incident flux on top from multilayer Eddington solution
!       I0UPOBS    incident radiance at bottom of layer at observation angle
!       TAUTOT     optical depth of layer
!       SSALB      single scattering albedo of layer
!       ASYM       asymmetry parameter of layer

!    Output Parameters:
!       IUPOBS     outgoing radiance at top at observation angle

!      Input and output radiances are in units of W m^-2 sr^-1 cm and
!    fluxes in units of W m^-2 cm.

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
          
    REAL ::    MUOBS, TTOP, TBOT
    REAL ::    WAVENO, FLUXUP,FLUXDN
    REAL ::    I0UPOBS
    REAL ::    IUPOBS
          
    REAL ::    SSALB, ASYM, TAUTOT
    REAL ::    BETA, PLANCK0, PLANCK1, COPLANCK
    REAL ::    EXPBETA, EXPMUOBS, EXPLAMP, EXPLAMM, EXPMAX
    REAL ::    T, R, LAMBDA
    REAL ::    F1UP, F0DN, FP1UP, FP0DN, FH1UP, FH0DN
    REAL ::    PARTFAC, PARTP, PARTM
    REAL ::    CFAC, CP, CM, DB, D1, D2, DP, DM, SB, SP, SM

    EXPMAX = 1.E08

    IF (TAUTOT < 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
    ENDIF
        
!           Planck function radiances in units of W m^-2 sr^-1 cm
    PLANCK0 = ttorad(WAVENO,TTOP)
    PLANCK1 = ttorad(WAVENO,TBOT)

!           Compute the Planck function variation coefficient
    BETA = LOG(PLANCK1/PLANCK0)/TAUTOT
    EXPBETA = PLANCK1/PLANCK0

!c the "exponential in radiance" variation is better than "exponential in T"
!c      print *,'beta,expbeta 0 = ',beta,expbeta
!c      BETA = LOG(TTOP/TBOT)/TAUTOT             !!!!!!
!c      EXPBETA = TTOP/TBOT                      !!!!!!
!c      print *,'beta,expbeta 1 = ',beta,expbeta
          
!      Find lambda so that it can be compared to beta
!      R, T are 2x2 matrix coefficents, lambda is eigenvalue
    T = (7-SSALB*(4+3*ASYM))*0.25
    R = (1-SSALB*(4-3*ASYM))*0.25
!           solution explodes if r = 0!!
    IF (ABS(R) < 1.E-5) R = 1.E-5
    LAMBDA = SQRT(T**2-R**2)
          
    IF(ABS(BETA**2-LAMBDA**2) < 1.E-4) THEN
        BETA = BETA/ABS(BETA)*SQRT(LAMBDA**2+1.E-4)
        EXPBETA = EXP(BETA*TAUTOT)
    ENDIF

    COPLANCK = (1-SSALB)*PLANCK0
    EXPMUOBS = EXP(-TAUTOT/MUOBS)
    IF (EXPMUOBS < 1./EXPMAX) EXPMUOBS = 1./EXPMAX

    EXPLAMP = EXP(LAMBDA*TAUTOT)
    IF (EXPLAMP > EXPMAX) EXPLAMP = EXPMAX
    EXPLAMM = 1.0/EXPLAMP

    F1UP = FLUXUP
    F0DN = FLUXDN

!           Compute particular solution at boundaries
    PARTFAC = COPLANCK/(LAMBDA**2-BETA**2)
    PARTP = PARTFAC*(T-R+BETA)
    PARTM = PARTFAC*(T-R-BETA)
    FP0DN = PARTM
    FP1UP = PARTP*EXPBETA
!           Subtract particular solution off to get homogeneous solution
    FH0DN = F0DN - FP0DN
    FH1UP = F1UP - FP1UP
!           Compute homogeneous coefficients C+ and C-
    CFAC = 1.0/(R**2*EXPLAMP - (LAMBDA-T)**2*EXPLAMM)
    CP = CFAC*(R*FH1UP - (LAMBDA-T)*EXPLAMM*FH0DN)
    CM = CFAC*(R*EXPLAMP*FH0DN - (LAMBDA-T)*FH1UP)
!           Make the constants for the 3 source function terms
    DB = COPLANCK + SSALB*PARTFAC*(2*(T-R) + 3*ASYM*BETA*MUOBS)
    D1 = R-T+LAMBDA
    D2 = 1.5*ASYM*(R+T-LAMBDA)*MUOBS
    DP = CP*SSALB*(D1+D2)
    DM = CM*SSALB*(D1-D2)

!           Finally do the answer
          
    IF(ABS(1-BETA*MUOBS) <= 1.E-5) THEN
        SB = DB*TAUTOT/MUOBS
    ELSE
        SB = DB*(1 - EXPBETA*EXPMUOBS)/(1-BETA*MUOBS)
    ENDIF
          
    IF(ABS(1-LAMBDA*MUOBS) <= 1.E-5) THEN
        SP = DP*TAUTOT/MUOBS
    ELSE
        SP = DP*(1 - EXPLAMP*EXPMUOBS)/(1-LAMBDA*MUOBS)
    ENDIF

    SM = DM*(1 - EXPLAMM*EXPMUOBS)/(1+LAMBDA*MUOBS)

    IUPOBS = I0UPOBS*EXPMUOBS + SB + SP + SM

    RETURN
    END SUBROUTINE EDDSCATRT

!************************************************************************
    SUBROUTINE SEDDSCATRT (MUOBS, TTOP, TBOT, WAVENO, &
    MUINC, IBNDRYUP, IBNDRYDN, &
    FLUXUPSEDD, FLUXDNSEDD, &
    FLUXUP, FLUXDN, I0UPOBS, &
    TAUTOT, SSALB, ASYM, &
    PHI1UP, PHI1DN, PHI2UP, PHI2DN,  IUPOBS)
!     Calculates upwelling radiance from a single thermally emitting
!     homogeneous layer using the single scattering/Eddington second
!     approximation hybrid radiative transfer method of Deeter and Evans (1998).

!    Input Parameters:
!       MUOBS      observation mu (cosine of viewing zenith angle)
!       TTOP       cloud top temperature (K)
!       TBOT       cloud bottom temperature (K)
!       WAVENO     wavenumber (cm^-1)
!       MUINC      cosine zenith angles for the two incident directions
!       IBNDRYUP   two incident radiances upwelling on bottom of layer
!       IBNDRYDN   two incident radiance downwelling on top of layer
!       FLUXUPSEDD incident flux on bottom from single scattering/Eddington
!       FLUXDNSEDD incident flux on top from single scattering/Eddington
!       FLUXUP     incident flux on bottom from multilayer Eddington solution
!       FLUXDN     incident flux on top from multilayer Eddington solution
!       I0UPOBS    incident radiance at bottom of layer at observation angle
!       TAUTOT     optical depth of layer
!       SSALB      single scattering albedo of layer
!       ASYM       asymmetry parameter of layer
!       PHI*       single scattering functions at the observation angle

!    Output Parameters:
!       IUPOBS     outgoing radiance at top at observation angle

!      Input and output radiances are in units of W m^-2 sr^-1 cm and
!    fluxes in units of W m^-2 cm.

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    REAL ::    MUOBS, TTOP, TBOT
    REAL ::    WAVENO, TAUTOT
    REAL ::    IBNDRYUP(2), IBNDRYDN(2)
    REAL ::    I0UPOBS
    REAL ::    IUPOBS, FLUXUP, FLUXDN, FLUXUPSEDD, FLUXDNSEDD
    REAL ::    MUINC(2)

    INTEGER :: IMU
    REAL ::    SSALB, ASYM, PHI1UP, PHI1DN, PHI2UP, PHI2DN
    REAL ::    BETA, COPLANCK, PLANCK0, PLANCK1
    REAL ::    EXPBETA, EXPMUOBS, DMUOBS, EXPMAX
    REAL ::    MU, DMUUP, DMUDN, EXPMU
    REAL ::    S1UP, IUP(2), IDN(2)
    REAL ::    TERM1, TERM2, TERM3
    REAL ::    EXPLAMP, EXPLAMM
    REAL ::    T, R, LAMBDA
    REAL ::    F1UP, F0DN, FP1UP, FP0DN, FH1UP, FH0DN
    REAL ::    PARTFAC, PARTP, PARTM
    REAL ::    CFAC, CP, CM, DB, D1, D2, DP, DM, SB, SP, SM

    REAL ::    SQR3, EXP3TP, EXP3TM, FAC3TP, FAC3TM, FAC3BETA
    REAL ::    FACF0DN, FACF1UP, FACG, FACGP, FACGM
    REAL ::    FACDENOM, DENOM1, DENOM2, FACNUM1, FACNUM2
    REAL ::    SEDTERM1, SEDTERM2, FACMUBET, FACNUM3, SEDTERM3, S1EDD
         
    IF (TAUTOT < 1.0E-6) THEN
        IUPOBS = I0UPOBS
        RETURN
    ENDIF

    EXPMAX = 1.E08
    SQR3 = SQRT(3.0)

!           Planck function radiances in units of W m^-2 sr^-1 cm

    PLANCK0 = ttorad(WAVENO,TTOP)
    PLANCK1 = ttorad(WAVENO,TBOT)
    BETA = LOG(PLANCK1/PLANCK0)/TAUTOT
    EXPBETA = PLANCK1/PLANCK0

!c the "exponential in radiance" variation is better than "exponential in T"
!c      print *,'beta,expbeta 0 = ',beta,expbeta
!c      BETA = LOG(TTOP/TBOT)/TAUTOT             !!!!!!
!c      EXPBETA = TTOP/TBOT                      !!!!!!
!c      print *,'beta,expbeta 1 = ',beta,expbeta

!      Find lambda so that it can be compared to beta
!      R, T are 2x2 matrix coefficents, lambda is eigenvalue
    T = (7-SSALB*(4+3*ASYM))*0.25
    R = (1-SSALB*(4-3*ASYM))*0.25
!           solution explodes if r = 0!!
    IF (ABS(R) < 1.E-5) R = 1.E-5
    LAMBDA = SQRT(T**2-R**2)

    IF(ABS(BETA**2-LAMBDA**2) < 1.E-4) THEN
        BETA = BETA/ABS(BETA)*SQRT(LAMBDA**2+1.E-4)
        EXPBETA = EXP(BETA*TAUTOT)
    ENDIF

    COPLANCK = (1. - SSALB)*PLANCK0
    EXPMUOBS = EXP(-TAUTOT/MUOBS)
    IF (EXPMUOBS < 1./EXPMAX) EXPMUOBS = 1./EXPMAX
    DMUOBS = 1. - BETA*MUOBS
    IF(ABS(DMUOBS) < 1.E-4) DMUOBS = SIGN(1.E-4,DMUOBS)
     
!         Accurate Single Scattering Part
!           IMU (1 or 2) is index for discrete mu's in simple model
    DO IMU = 1, 2
        MU = MUINC(IMU)
        EXPMU = EXP(-TAUTOT/MU)
        IF (EXPMU < 1./EXPMAX) EXPMU = 1./EXPMAX
        DMUUP = 1. - BETA*MU
        DMUDN = 1. + BETA*MU
        IF(ABS(DMUUP) < 1.E-4) DMUUP = SIGN(1.E-4,DMUUP)
        IF(ABS(DMUDN) < 1.E-4) DMUDN = SIGN(1.E-4,DMUDN)

        TERM1 = IBNDRYUP(IMU)*EXPMU - COPLANCK*EXPBETA*EXPMU/DMUUP
        TERM2 = (MU/(MUOBS - MU))*(1. - EXPMUOBS/EXPMU)
        TERM3 = COPLANCK*(1. - EXPBETA*EXPMUOBS)/(DMUUP*DMUOBS)
        IUP(IMU) = - (TERM1*TERM2 - TERM3)
         
        TERM1 = IBNDRYDN(IMU) - COPLANCK/DMUDN
        TERM2 = (MU/(MUOBS + MU))*(1. - EXPMUOBS*EXPMU)
        TERM3 =  COPLANCK*(1. - EXPBETA*EXPMUOBS)/(DMUDN*DMUOBS)
        IDN(IMU) = (TERM1*TERM2 + TERM3)
    ENDDO

!           First order scattering contribution
    S1UP = (SSALB/2.)*(PHI1UP*IUP(1) + PHI2UP*IUP(2) &
    + PHI1DN*IDN(1) + PHI2DN*IDN(2))

!         Eddington part

    EXPLAMP = EXP(LAMBDA*TAUTOT)
    IF (EXPLAMP > EXPMAX) EXPLAMP = EXPMAX
    EXPLAMM = 1.0/EXPLAMP

!           Fluxes sent into subroutine from output of multilayer solution
    F1UP = FLUXUP
    F0DN = FLUXDN

!           Compute particular solution at boundaries
    PARTFAC = COPLANCK/(LAMBDA**2-BETA**2)
    PARTP = PARTFAC*(T-R+BETA)
    PARTM = PARTFAC*(T-R-BETA)
    FP0DN = PARTM
    FP1UP = PARTP*EXPBETA
!           Subtract particular solution off to get homogeneous solution
    FH0DN = F0DN - FP0DN
    FH1UP = F1UP - FP1UP
!           Compute homogeneous coefficients C+ and C-
    CFAC = 1.0/(R**2*EXPLAMP - (LAMBDA-T)**2*EXPLAMM)
    CP = CFAC*(R*FH1UP - (LAMBDA-T)*EXPLAMM*FH0DN)
    CM = CFAC*(R*EXPLAMP*FH0DN - (LAMBDA-T)*FH1UP)
!           Make the constants for the 3 source function terms
    DB = COPLANCK + SSALB*PARTFAC*(2*(T-R) + 3*ASYM*BETA*MUOBS)
    D1 = R-T+LAMBDA
    D2 = 1.5*ASYM*(R+T-LAMBDA)*MUOBS
    DP = CP*SSALB*(D1+D2)
    DM = CM*SSALB*(D1-D2)

!           Finally do the answer

    IF(ABS(1-BETA*MUOBS) <= 1.E-5) THEN
        SB = DB*TAUTOT/MUOBS
    ELSE
        SB = DB*(1 - EXPBETA*EXPMUOBS)/(1-BETA*MUOBS)
    ENDIF
          
    IF(ABS(1-LAMBDA*MUOBS) <= 1.E-5) THEN
        SP = DP*TAUTOT/MUOBS
    ELSE
        SP = DP*(1 - EXPLAMP*EXPMUOBS)/(1-LAMBDA*MUOBS)
    ENDIF

    SM = DM*(1 - EXPLAMM*EXPMUOBS)/(1+LAMBDA*MUOBS)

!           Compute the Eddington single scattering solution
    F0DN = FLUXDNSEDD
    F1UP = FLUXUPSEDD
    EXP3TP = EXP(SQR3*TAUTOT)
    IF (EXP3TP > EXPMAX) EXP3TP = EXPMAX
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

!           Zero order scattering contribution
!       S0UP = (1. - EXPBETA*EXPMUOBS)*COPLANCK/DMUOBS

    IUPOBS = I0UPOBS*EXPMUOBS + SB + SP + SM - S1EDD + S1UP

    RETURN
    END SUBROUTINE SEDDSCATRT

!************************************************************************
    SUBROUTINE EDDRTF (NLAYER, OPTDEPTHS, ALBEDOS, &
    ASYMMETRIES, TEMPS, WAVENO, &
    FLUXTOP, FLUXBOT,  FLUXES)
!       EDDRTF computes the layer interface fluxes for a multilayer
!     plane-parallel atmosphere with thermal sources of radiation using the
!     Eddington approximation.  The medium is specified by a number
!     of homogeneous layers.  The Planck function is assumed linear with
!     optical depth (slightly inconsistent with the exponential assumption
!     in the single layer routines).  The temperatures, optical depth,
!     single scattering albedo, and asymmetry parameter are specified
!     for each layer.  The boundary conditions are the fluxes incident
!     on the top and bottom of the domain. The Eddington fluxes at each
!     level are returned.
!       The model works by calculating the reflection, transmission, and
!     source terms for each layer from the input properties.  A
!     tri-diagonal matrix solver is then used to compute the fluxes
!     at each layer interface from the applied boundary conditions.

!     Parameters:
!       Input:
!     NLAYER         integer      Number of homogenous layers
!                                  (layers are specified from the top down)
!     OPTDEPTHS      real array   Optical thickness of layers
!     ALBEDOS        real array   Single scattering albedos
!     ASYMMETRIES    real array   Asymmetry parameters
!     TEMPS          real array   Temperatures (K) at layer interfaces
!                                  (e.g. TEMPS(1) is at top of top layer,
!                                   TEMPS(2) is at bottom of top layer).
!     WAVENO         real         wavenumber (cm^-1, for Planck function)

!       Output:
!     FLUXES         real         Eddington fluxes at layer interfaces.
!                                   FLUXES(1,L) is upwelling,
!                                   FLUXES(2,L) is downwelling,
!                                   L=1 is top, L=NUML+1 is bottom

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER ::   NLAYER
    REAL ::      TEMPS(NLAYER+1)
    REAL ::      OPTDEPTHS(NLAYER)
    REAL ::      ALBEDOS(NLAYER)
    REAL ::      ASYMMETRIES(NLAYER)
    REAL ::      WAVENO
    REAL ::      FLUXES(3,(NLAYER+1))
     
!               MAXLAY is maximum number of layers
    INTEGER ::   MAXLAY, MAXN
    PARAMETER (MAXLAY=200, MAXN=2*MAXLAY+1)
    INTEGER ::   N, L, I
    REAL*8 ::    DELTAU, G, OMEGA
    REAL*8 ::    LAMBDA, R, T, D, CM, CP, A, B, X1, X2
    REAL*8 ::    REFLECT, TRANS, SOURCEP, SOURCEM
    REAL*8 ::    RADP1P, RADP1M, RADP2P, RADP2M
    REAL*8 ::    PI, PLANCK1, PLANCK2, TAU
    REAL*8 ::    EXLP, EXLM, V
    REAL*8 ::    LOWER(MAXN), UPPER(MAXN), DIAG(MAXN), RHS(MAXN)
    REAL ::      FLUXTOP, FLUXBOT
    PARAMETER (PI=3.1415926535)
     
!               Compute the reflection, transmission, and source
!               coefficients for each layer for the diffuse Eddington
!               two stream problem.


    N = 2*NLAYER+2
    IF (N > MAXN) STOP 'EDDRTF: Exceeded maximum number of layers'

    PLANCK1 = 0.5*ttorad(WAVENO,TEMPS(1))
    TAU = 0.0
    I = 2
    DO L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        IF (DELTAU < 0.0) STOP 'EDDRTF: TAU<0'
    !            Special case for zero optical depth
        IF (DELTAU == 0.0) THEN
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
        !              Special case for conservative scattering (lambda=0)
            IF (LAMBDA == 0.0) THEN
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
                      
        !               Calculate thermal source terms

            PLANCK2 = 0.5*ttorad(WAVENO,TEMPS(L+1))

            V = 2.0*(PLANCK2-PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
            RADP1P = -V + PLANCK1
            RADP2M =  V + PLANCK2
            RADP2P = -V + PLANCK2
            RADP1M =  V + PLANCK1
            IF (LAMBDA == 0.0) THEN
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
          
!           Setup for and call the tri-diagonal matrix solver
    RHS(1) = FLUXTOP
    DIAG(1) = 0.0
    UPPER(1) = 1.0
!      DIAG(N) = -(1.0-GNDEMIS)
    DIAG(N) = 0.0
    LOWER(N) = 1.0
    RHS(N) = FLUXBOT
    CALL TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
!           Put the fluxes in the output array
         
    I = 1
    DO L = 1, NLAYER+1
        FLUXES(1,L) = RHS(I)
        FLUXES(2,L) = RHS(I+1)
        I = I + 2
    ENDDO
     
    RETURN
    END SUBROUTINE EDDRTF
     
!************************************************************************
    SUBROUTINE TRIDIAG (N, LOWER, DIAG, UPPER, RHS)
!       Computes the solution to a tridiagonal system.
!       N is order of the matrix.  LOWER(2..N) is the subdiagonal,
!       DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the
!       superdiagonal.  On input RHS is the right hand side, while
!       on output it is the solution vector.  Everything is destroyed.
!       Hacked from Linpack DGTSL.

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: N
    REAL*8 ::  LOWER(*), DIAG(*), UPPER(*), RHS(*)
    INTEGER :: K, KB
    REAL*8 ::  T

    IF (N == 1) THEN
        IF (DIAG(1) == 0.0) GOTO 990
        RHS(1) = RHS(1)/DIAG(1)
    ENDIF
    LOWER(1) = DIAG(1)
    DIAG(1) = UPPER(1)
    UPPER(1) = 0.0
    UPPER(N) = 0.0
    DO K = 1, N-1
    !              Interchange this and next row to the get the largest pivot.
        IF (ABS(LOWER(K+1)) >= ABS(LOWER(K))) THEN
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
        IF (LOWER(K) == 0.0) GOTO 990
        T = -LOWER(K+1)/LOWER(K)
        LOWER(K+1) = DIAG(K+1) + T*DIAG(K)
        DIAG(K+1) = UPPER(K+1) + T*UPPER(K)
        UPPER(K+1) = 0.0
        RHS(K+1) = RHS(K+1) + T*RHS(K)
    ENDDO
    IF (LOWER(N) == 0.0) GOTO 990

!           Back substitute
    RHS(N) = RHS(N)/LOWER(N)
    RHS(N-1) = (RHS(N-1) - DIAG(N-1)*RHS(N))/LOWER(N-1)
    DO KB = 1, N-2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) -DIAG(K)*RHS(K+1) -UPPER(K)*RHS(K+2))/LOWER(K)
    ENDDO
    RETURN

    990 CONTINUE
    STOP 'Singular matrix in TRIDIAG'
    END SUBROUTINE TRIDIAG
     

!************************************************************************
    SUBROUTINE GASRT1 (MU, WAVENO, SFCTEMP, SFCEMIS, &
    NLEV, TEMP, TAU, COALB, ICLDTOP, &
    ICLDBOT, NCLDLAY, RAD1UP, RAD0DN, iiDiv)
!       Compute the upwelling radiance at the bottom of each cloud layer
!     (RAD1UP) and the downwelling radiance at the top of each cloud
!     layer (RAD0DN).  Integrates the nonscattering RTE through the domain
!     from top to bottom and then bottom to top.  Only the absorption
!     part of the optical depth in the cloud is considered.  The bottom
!     of the lowest cloud layer is at level ICLDBOT and the top of the
!     highest cloud layer is as level ICLDTOP.  The radiance is computed
!     at the angle MU and wavenumber WAVENO.  The temperature profile (TEMP),
!     layer optical depth profile (TAU), and single scattering coalbedo
!     (COALB=1-omega) profile are input.  There is specular reflection from
!     the surface yet no multiple scattering effect due of the surface, ie.
!     The gas integrated down from TOA and reflected off surface once with
!     reflectivity equal to one minus emissivity, and this is the initial
!     condition when integrating gas incident on bottom of cloud.

! modify this so that we also do stuff at arccos(3/5)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: NLEV, ICLDTOP, ICLDBOT, NCLDLAY, iiDiv
    REAL ::    MU, WAVENO, SFCTEMP, SFCEMIS
    REAL ::    TEMP(KPROFLAYER+1), TAU(KPROFLAYER), COALB(KPROFLAYER)
    REAL ::    RAD1UP(NCLDLAY), RAD0DN(NCLDLAY)
    REAL ::    RAD1UP0, RAD0DN0
    INTEGER :: I, L, IEND
    REAL ::    PLANCK0, PLANCK1, PLANCKSFC, TAUO, TRANS, DELPLANCK

    REAL :: r1,r2
    REAL :: r35,mu35,rZ,rBooga,rBeta,rTT,rZeta

    INTEGER :: iCldTopA,iCldBotA

    iCldTopA = iCldTop - kProfLayer*iiDiv
    iCldBotA = iCldBot - kProfLayer*iiDiv

!         Loop through layers, starting at top and going down
!      Continue past cloud layer to surface to compute specular
!      reflection.

! ---- at the diffusive angle
    mu35 = 5.0/3.0
    r35 = ttorad(WAVENO,sngl(kTSpace))
    PLANCK0 = ttorad(WAVENO,TEMP(1))
    L = 1
    IF(SFCEMIS < 1) IEND = NLEV-1
    IF(SFCEMIS == 1) IEND = ICLDBOTA-2
          
    DO I = 1, IEND
    !           Planck function radiances in units of W m^-2 sr^-1 cm
        PLANCK1 = ttorad(WAVENO,TEMP(I+1))
        TAUO = TAU(I)/MU35
        IF (TAU(I) < 0.001) THEN
            r35 = r35*(1-TAUO) &
            + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
        ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAUO
            r35 = r35*TRANS + COALB(I)*( PLANCK1-DELPLANCK &
            - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
        ENDIF
        PLANCK0 = PLANCK1
    ENDDO
    r35 = r35 * 2 * kPi

! ---- at the downward view angle
    RAD0DN0 = r1 *WAVENO**3 &
    / (EXP(r2*WAVENO/sngl(kTSpace)) - 1)
!      RAD0DN0 = 0.0

    IF(ICLDBOTA == 2) RAD0DN(1) = RAD0DN0

    PLANCK0 = ttorad(WAVENO,TEMP(1))
    L = 1
    IF(SFCEMIS < 1) IEND = NLEV-1
    IF(SFCEMIS == 1) IEND = ICLDBOTA-2
          
    DO I = 1, IEND
        PLANCK1 = ttorad(WAVENO,TEMP(1+1))
        TAUO = TAU(I)/MU
        IF (TAU(I) < 0.001) THEN
            RAD0DN0 = RAD0DN0*(1-TAUO) &
            + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
        ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAUO
            RAD0DN0 = RAD0DN0*TRANS + COALB(I)*( PLANCK1-DELPLANCK &
            - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
        ENDIF

        IF (I+1 >= ICLDTOPA .AND. I <= (ICLDBOTA-2)) THEN
            RAD0DN(L) = RAD0DN0
            L = L + 1
        ENDIF
        PLANCK0 = PLANCK1
    ENDDO

! ---- at the upward view angle
!          Loop through layers, starting at bottom and going up
    L = NCLDLAY
    PLANCKSFC = ttorad(WAVENO,SFCTEMP)

!      RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)
!      RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)/kPi

    RAD1UP0 = SFCEMIS*PLANCKSFC+R35*(1-SFCEMIS)/kPi

    rZ = RAD1UP0
!      print *,'here',SFCEMIS,PLANCKSFC,R35,RAD1UP0,rZ

    IF(ICLDTOPA == (NLEV-1)) RAD1UP(1) = RAD1UP0

    PLANCK1 = ttorad(WAVENO,TEMP(NLEV))
    DO I = NLEV-1, ICLDTOPA+1, -1

        PLANCK0 = ttorad(WAVENO,TEMP(I))
        TAUO = TAU(I)/MU
        IF (TAUO < 0.001) THEN
            RAD1UP0 = RAD1UP0*(1-TAUO) &
            + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
        ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAUO
            RAD1UP0 = RAD1UP0*TRANS + COALB(I)*( PLANCK0+DELPLANCK &
            - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO)) )
        ENDIF
        rZ = rZ*exp(-TAUO) + (1-exp(-TAUO))*PLANCK1
                
    !        print *,'---->',I,MU,TEMP(I+1),TEMP(I),TAU(I),RAD1UP0,rZ
        IF (I <= ICLDBOTA) THEN
            RAD1UP(L) = RAD1UP0
            L = L - 1
        ENDIF

        PLANCK1 = PLANCK0
    ENDDO

! ---- at the upward view angle, using the twostream code
!          Loop through layers, starting at bottom and going up
    L = NCLDLAY
    PLANCKSFC = r1 *WAVENO**3 &
    / (EXP(r2*WAVENO/SFCTEMP) - 1)

    RAD1UP0 = SFCEMIS*PLANCKSFC+R35*(1-SFCEMIS)/kPi

    rZ = RAD1UP0

    IF(ICLDTOPA == (NLEV-1)) RAD1UP(1) = RAD1UP0

    PLANCK1 = ttorad(WAVENO,TEMP(NLEV))
    DO I = NLEV-1, ICLDTOPA+1, -1
    !        rBooga = log(TEMP(iBeta+1)/TEMP(iBeta))
    !            rbeta = 1/raaAbs(iFr,iL) * rBooga
    !            rTT   = ttorad(raFreq(iFr),TEMP(iBeta))/(1 + rbeta*rCos)
    !            rZeta = (raInten(iFr) - rTT) * exp(-raaAbs(iFr,iL)/rCos)
    !            raInten(iFr) = rZeta + rTT * exp(raaAbs(iFr,iL) * rbeta)

        TAUO = TAU(I)/MU
        rBooga = log(TEMP(I)/TEMP(I+1))
        rBeta = 1/TAU(I) * rBooga
        rTT = ttorad(WAVENO,TEMP(I+1))/(1+rbeta*MU)
        rZeta = (RAD1UP0-rTT) * exp(-TAUO)
        RAD1UP0 = rZeta + rTT * exp(TAU(I)*rBeta)
                
        IF (I <= ICLDBOTA) THEN
            RAD1UP(L) = RAD1UP0
            L = L - 1
        ENDIF
    ENDDO
            
    RETURN
    END SUBROUTINE GASRT1

!************************************************************************
    SUBROUTINE GASRT2 (WAVENO, NLEV, TEMP, TAU, ICLD, IOBS, MU, RAD, iiDiv)
!       Compute the upwelling radiance starting at the cloud top (ICLD)
!     and going to the observation level (IOBS) or the downwelling
!     radiance starting at cloud bottom (ICLD) and going down to
!     observation level (IOBS) depending on MU.
!     The radiance is computed at the angle MU and wavenumber WAVENO.
!     The temperature profile (TEMP) and layer optical depth profile (TAU)
!     are input.

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: NLEV, ICLD, IOBS, iiDiv
    REAL ::    WAVENO, TEMP(KPROFLAYER+1), TAU(KPROFLAYER), MU, RAD
    INTEGER :: I
    REAL ::    PLANCK0, PLANCK1, TAUO, TRANS, DELPLANCK

    REAL :: r35,mu35,rZ,rBooga,rBeta,rTT,rZeta,radrad

!       Loop through layers, if looking down start at bottom and go up

    radrad = rad

    IF(IOBS < ICLD) THEN
        PLANCK1 = ttorad(WAVENO,TEMP(ICLD-iiDiv*kProfLayer))
        DO I = ICLD-1, IOBS, -1
            PLANCK0 = ttorad(WAVENO,TEMP(I-iiDiv*kProfLayer))
            TAUO = TAU(I-iiDiv*kProfLayer)/MU
            IF (TAUO < 0.001) THEN
                RAD = RAD*(1-TAUO) + TAUO*0.5*(PLANCK0+PLANCK1)
            ELSE
                TRANS = EXP(-TAUO)
                DELPLANCK = (PLANCK1-PLANCK0)/TAUO
                RAD = RAD*TRANS + PLANCK0+DELPLANCK &
                - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO))
            ENDIF
            PLANCK1 = PLANCK0
        ENDDO

        DO I = ICLD-1, IOBS, -1
            TAUO = TAU(I)/MU
            rBooga = log(TEMP(I)/TEMP(I+1))
            rBeta = 1/TAU(I) * rBooga
            rTT = ttorad(WAVENO,TEMP(I+1))/(1+rbeta*MU)
            rZeta = (radrad-rTT) * exp(-TAUO)
            radrad = rZeta + rTT * exp(TAU(I)*rBeta)
        ENDDO
    !        print *,'rad, radrad = ',rad,radrad
        rad = radrad

    !       Loop through layers, if looking up start at top and go down
    ELSEIF(IOBS > ICLD) THEN
        PLANCK0 = ttorad(WAVENO,TEMP(ICLD-iiDiv*kProfLayer))
        DO I = ICLD, IOBS-1
        !           Planck function radiances in units of W m^-2 sr^-1 um^-1
            PLANCK1 = ttorad(WAVENO,TEMP(I+1-iiDiv*kProfLayer))
            TAUO = TAU(I-iiDiv*kProfLayer)/MU
            IF (TAU(I-iiDiv*kProfLayer) < 0.001) THEN
                RAD = RAD*(1-TAUO) &
                + TAUO*0.5*(PLANCK0+PLANCK1)
            ELSE
                TRANS = EXP(-TAUO)
                DELPLANCK = (PLANCK1-PLANCK0)/TAUO
                RAD = RAD*TRANS + ( PLANCK1-DELPLANCK &
                - TRANS*(PLANCK1-DELPLANCK*(1.0+TAUO)) )
            ENDIF
            PLANCK0 = PLANCK1
        ENDDO
    ENDIF

    RETURN
    END SUBROUTINE GASRT2

!************************************************************************
!************************************************************************
!************************************************************************
! slight modifications to GASRT1 done by Sergio, to allow backgnd thermal to
! be computed or not; Planck terms stored in array instead of computing twice
! this computes backgnd thermal using fast diffusivity approx
    SUBROUTINE GASRT1_nocloud (MU, WAVENO, SFCTEMP, SFCEMIS, &
    NLEV, TEMP, TAU, COALB, ICLDTOP, &
    ICLDBOT, NCLDLAY, RAD1UP, RAD0DN, &
    ibdry, TOA_to_instr, iiDiv)
!       Compute the upwelling radiance at the bottom of each cloud layer
!     (RAD1UP) and the downwelling radiance at the top of each cloud
!     layer (RAD0DN).  Integrates the nonscattering RTE through the domain
!     from top to bottom and then bottom to top.  Only the absorption
!     part of the optical depth in the cloud is considered.  The bottom
!     of the lowest cloud layer is at level ICLDBOT and the top of the
!     highest cloud layer is as level ICLDTOP.  The radiance is computed
!     at the angle MU and wavenumber WAVENO.  The temperature profile (TEMP),
!     layer optical depth profile (TAU), and single scattering coalbedo
!     (COALB=1-omega) profile are input.  There is specular reflection from
!     the surface yet no multiple scattering effect due of the surface, ie.
!     The gas integrated down from TOA and reflected off surface once with
!     reflectivity equal to one minus emissivity, and this is the initial
!     condition when integrating gas incident on bottom of cloud.

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: NLEV, ICLDTOP, ICLDBOT, NCLDLAY,ibdry, iiDiv
    REAL ::    MU, WAVENO, SFCTEMP, SFCEMIS
    REAL ::    TEMP(KPROFLAYER+1), TAU(KPROFLAYER), COALB(KPROFLAYER)
    REAL ::    RAD1UP(NCLDLAY), RAD0DN(NCLDLAY)
    REAL ::    RAD1UP0, RAD0DN0, TOA_to_instr
    INTEGER :: I, L
    REAL ::    PLANCK0, PLANCK1, PLANCKSFC, TAUO, TRANS, DELPLANCK

    REAL :: raPlanck(kProfLayer+1)

    INTEGER :: iCldTopA,iCldBotA

    iCldTopA = iCldTop - kProfLayer*iiDiv
    iCldBotA = iCldBot - kProfLayer*iiDiv

    DO i=1,NLEV
        raPlanck(i) = ttorad(WAVENO,TEMP(I))
    END DO

!     Loop through layers, starting at top and going down
!     Continue past cloud layer to surface to compute specular
!     reflection.

    IF (kThermal >= 0) THEN
        CALL FastBDRYL2GDiffusive_rts(TOA_to_instr, &
        MU, WAVENO,NLEV, TEMP, TAU,RAD0DN, RAD0DN0,ibdry)
    END IF

!     Loop through layers, starting at bottom and going up
!     note because of defn of emissivity used in FastBdry, we multiply
!     rad0dn0 by (1-ems)/pi
    L = NCLDLAY
    PLANCKSFC = ttorad(WAVENO,SFCTEMP)

    IF (kThermal >= 0) THEN
    !        RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)
        RAD1UP0 = SFCEMIS*PLANCKSFC+RAD0DN0*(1-SFCEMIS)/kPi
    ELSE
        RAD1UP0 = SFCEMIS*PLANCKSFC
    END IF

    IF(ICLDTOPA == (NLEV-1)) RAD1UP(1) = RAD1UP0

! notice if DELPLANCK == 0 then below algorithm reduces to
! rad=sum( (1-tau(i))B(i) + rad*tau(i))
! which is the same as that of the usual kCARTA algorithm
!c        DO iFr=1,kMaxPts
!c          raInten(iFr)=raaEmission(iFr,iLay)+
!c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
!c          END DO
!c raaEmission = (1-tau(i))B(i), raaLayTran = tau(i)
! except that there the Planck term B(i) is at the "middle" of the layer,
! while here it is at the top of the layer, which is colder .. hence this
! overall radiance should be smaller than than of kCARTA!

! so if DELPLANK is non zero, but actually positive, this small correction
! will bring up back to kCARTA radiance levels!!! woohoo
    PLANCK1 = raPlanck(NLEV)
    DO I = NLEV-1, ICLDTOPA+1, -1
        PLANCK0 = raPlanck(I)
        TAUO = TAU(I)/MU
        IF (TAUO < 0.001) THEN
            RAD1UP0 = RAD1UP0*(1-TAUO) &
            + TAUO*0.5*(PLANCK0+PLANCK1)*COALB(I)
        ELSE
            TRANS = EXP(-TAUO)
            DELPLANCK = (PLANCK1-PLANCK0)/TAUO
            RAD1UP0 = RAD1UP0*TRANS + COALB(I)*( PLANCK0+DELPLANCK &
            - TRANS*(PLANCK0+DELPLANCK*(1.0+TAUO)) )
        ENDIF
        IF (I <= ICLDBOTA) THEN
            RAD1UP(L) = RAD1UP0
            L = L - 1
        ENDIF
        PLANCK1 = PLANCK0
    ENDDO

    RETURN
    end SUBROUTINE GASRT1_nocloud

!************************************************************************
END MODULE scatter_rtspec_code
