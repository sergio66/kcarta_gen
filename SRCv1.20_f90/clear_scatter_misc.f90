! Copyright 2006
! University of Maryland Baltimore County
! All Rights Reserved

MODULE clear_scatter_misc

USE basic_common
USE kcoeff_common
USE spline_and_sort_and_common
USE rad_diff_and_quad
USE clear_scatter_basic

IMPLICIT NONE

CONTAINS

! ************************************************************************
! ******************* THESE ARE THE RTSPEC ROUTINES **********************
! ************************************************************************
    SUBROUTINE JACOBIAN_INTERP_SCAT_TABLE2 (WAVENO, DME, &
    dEXTINCT_dr, dSSALB_dr, dASYM_dr, &
    NDME, DMETAB, NWAVE, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM)
!       Interpolates the scattering properties from the table for
!     a particular wavenumber and particle size.  Does a bilinear
!     interpolation, but optimized for fixed particle size and slowly
!     varying wavenumber.  If the DME is the same as last time then we
!     can just linearly interpolate in wavenumber between stored
!     scattering values.  If the DME has changed then we linearly
!     interpolate between the DMETAB grid lines for each of the two
!     wavenumber grid lines.

! also computes the derivatives wrt particle size   d = 2r!!!!

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL ::     WAVENO, DME
    REAL ::     dEXTINCT_dr, dSSALB_dr, dASYM_dr
    INTEGER ::  NDME, NWAVE
    REAL ::     DMETAB(NDME), WAVETAB(NWAVE)
    REAL ::     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME)
    REAL ::     TABASYM(NWAVE,NDME)
    INTEGER ::  IW0, IW1, ID, IL, IU, IM
    LOGICAL ::  NEWIW
    REAL ::     FWAV, FDME, FLDME, F
    REAL ::     OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
    REAL ::     dEXT0, dEXT1, dALB0, dALB1, dASYM0, dASYM1
    DATA     IW0/1/, IW1/2/
    REAL ::     EXTINCT, SSALB, ASYM

    iw0 = 1
    iw1 = 2
    olddme = 0.0

    iw0 = 1
    iw1 = nwave
    olddme = -10.0

!         Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
        write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
        IF (WAVENO < WAVETAB(1)) THEN
            WAVENO = WAVETAB(1)
        ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
            WAVENO = WAVETAB(NWAVE)
        END IF
    ! ALL DoStop
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
    !        write(kStdErr,*) DME,' outside ',DMETAB(1),':',DMETAB(NDME)
    !        write(kStdErr,*) 'INTERP_SCAT_TABLE: particle Dme out of range ... RESET'
        IF (DME < DMETAB(1)) THEN
            DME = DMETAB(1)
        ELSEIF (DME > DMETAB(NDME)) THEN
            DME = DMETAB(NDME)
        END IF
    ! ALL DoStop
    END IF

!     See if wavenumber is within last wavenumber grid, otherwise
!     find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
!      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (WAVENO >= WAVETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE. 
    ENDIF

    IF (DME /= OLDDME) THEN
    !         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (DME >= DMETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
    ENDIF

    IF (DME /= OLDDME .OR. NEWIW) THEN
    !         If not the same Dme or a new wavenumber grid, then
    !           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )

        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)

        dEXT0 = EXT0/DME * LOG(TABEXTINCT(IW0,ID+1)/TABEXTINCT(IW0,ID))/ &
        LOG(DMETAB(ID+1)/DMETAB(ID))
        dEXT1 = EXT1/DME * LOG(TABEXTINCT(IW1,ID+1)/TABEXTINCT(IW1,ID))/ &
        LOG(DMETAB(ID+1)/DMETAB(ID))

        dALB0 = (TABSSALB(IW0,ID+1)-TABSSALB(IW0,ID))/(DMETAB(ID+1)-DMETAB(ID))
        dALB1 = (TABSSALB(IW1,ID+1)-TABSSALB(IW1,ID))/(DMETAB(ID+1)-DMETAB(ID))

        dASYM0 = (TABASYM(IW0,ID+1)-TABASYM(IW0,ID))/(DMETAB(ID+1)-DMETAB(ID))
        dASYM1 = (TABASYM(IW1,ID+1)-TABASYM(IW1,ID))/(DMETAB(ID+1)-DMETAB(ID))

    ENDIF

!     Linearly interpolate the scattering properties in wavenumber
    FWAV    = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
    F       = 1-FWAV
    dEXTINCT_dr = F*dEXT0 + FWAV*dEXT1
    dSSALB_dr   = F*dALB0 + FWAV*dALB1
    dASYM_dr    = F*dASYM0 + FWAV*dASYM1

    OLDDME = DME

    RETURN
    END SUBROUTINE JACOBIAN_INTERP_SCAT_TABLE2

!************************************************************************
    SUBROUTINE INTERP_SCAT_TABLE3 (MU, WAVENO, DME, &
    EXTINCT, SSALB, ASYM, &
    PHI1UP, PHI1DN, PHI2UP, PHI2DN, &
    NMU, MUTAB, NDME, DMETAB, NWAVE, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
!       Interpolates the scattering properties from the table for
!     a particular observation angle, wavenumber, and particle size.
!     Does a trilinear interpolation, but optimized for fixed viewing
!     angle and particle size and slowly varying wavenumber.  If the
!     DME and MU are the same as last time then we can just linearly
!     interpolate in wavenumber between stored scattering values.
!     If the DME or MU has changed then we bilinearly interpolate between
!     the DMETAB and MUTAB grid lines for each of the two wavenumber
!     grid lines.

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
     
    REAL ::     MU, WAVENO, DME
    REAL ::     EXTINCT, SSALB, ASYM, PHI1UP, PHI1DN, PHI2UP, PHI2DN
    INTEGER ::  NMU, NDME, NWAVE
    REAL ::     MUTAB(NMU), DMETAB(NDME), WAVETAB(NWAVE)
    REAL ::     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME)
    REAL ::     TABASYM(NWAVE,NDME)
    REAL ::     TABPHI1UP(NMU,NWAVE,NDME), TABPHI1DN(NMU,NWAVE,NDME)
    REAL ::     TABPHI2UP(NMU,NWAVE,NDME), TABPHI2DN(NMU,NWAVE,NDME)
    INTEGER ::  IW0, IW1, ID, IMU, IL, IU, IM
    LOGICAL ::  NEWIW
    REAL ::     FWAV, FDME, FLDME, FMU, F1, F2, F3, F4
    REAL ::     OLDMU, OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
    REAL ::     PHI1UP0, PHI1UP1, PHI1DN0, PHI1DN1
    REAL ::     PHI2UP0, PHI2UP1, PHI2DN0, PHI2DN1
! sergio do not save iw0,iw1
!      SAVE     IW0, IW1, IMU,ID, OLDMU, OLDDME, FDME, FLDME, FMU
    SAVE     IW0, IW1, IMU,ID, OLDMU, OLDDME, FDME, FLDME, FMU
    SAVE     EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
    SAVE     PHI1UP0, PHI1UP1, PHI1DN0, PHI1DN1
    SAVE     PHI2UP0, PHI2UP1, PHI2DN0, PHI2DN1
    DATA     IW0/1/, IW1/2/
     
    iw0 = 1
    iw1 = 2
         
!         Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. &
    WAVENO > WAVETAB(NWAVE)) THEN
        write(kStdErr,*) WAVENO,WAVETAB(1),WAVETAB(NWAVE)
        write (kStdErr,*) 'INTERP_SCAT_TABLE3: wavenumber out of range.'
        CALL DoStop
    ENDIF

    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
        WRITE(*,*) 'Dme: ', DME, DMETAB(1), DMETAB(NDME)
        STOP 'INTERP_SCAT_TABLE: particle Dme out of range.'
    ENDIF
!      IF (MU .LT. MUTAB(1) .OR. MU .GT. MUTAB(NMU))
!     .  STOP 'INTERP_SCAT_TABLE3: viewing angle mu out of range.'

!         See if wavenumber is within last wavenumber grid, otherwise
!           find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
    IF (WAVENO < WAVETAB(IW0) .OR. WAVENO > WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (WAVENO >= WAVETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE. 
    ENDIF

    IF (MU /= OLDMU) THEN
    !         Find the grid location and interpolation factor for MU
        IL=1
        IU=NMU
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (MU >= MUTAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IMU = MAX(IL,1)
        FMU = (MU-MUTAB(IMU))/(MUTAB(IMU+1)-MUTAB(IMU))
    ENDIF
    IF (DME /= OLDDME) THEN
    !         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (DME >= DMETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
    ENDIF

!         If not the same Dme and mu, then bilinearly interpolate things
!           Logarithmically interpolate extinction
    IF (DME /= OLDDME .OR. MU /= OLDMU .OR. NEWIW) THEN
        F1 = (1-FMU)*(1-FDME)
        F2 = (1-FMU)*FDME
        F3 = FMU*(1-FDME)
        F4 = FMU*FDME
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)

        PHI1UP0 = F1*TABPHI1UP(IMU,IW0,ID) + F2*TABPHI1UP(IMU,IW0,ID+1) &
        + F3*TABPHI1UP(IMU+1,IW0,ID) + F4*TABPHI1UP(IMU+1,IW0,ID+1)
        PHI1UP1 = F1*TABPHI1UP(IMU,IW1,ID) + F2*TABPHI1UP(IMU,IW1,ID+1) &
        + F3*TABPHI1UP(IMU+1,IW1,ID) + F4*TABPHI1UP(IMU+1,IW1,ID+1)

        PHI1DN0 = F1*TABPHI1DN(IMU,IW0,ID) + F2*TABPHI1DN(IMU,IW0,ID+1) &
        + F3*TABPHI1DN(IMU+1,IW0,ID) + F4*TABPHI1DN(IMU+1,IW0,ID+1)
        PHI1DN1 = F1*TABPHI1DN(IMU,IW1,ID) + F2*TABPHI1DN(IMU,IW1,ID+1) &
        + F3*TABPHI1DN(IMU+1,IW1,ID) + F4*TABPHI1DN(IMU+1,IW1,ID+1)

        PHI2UP0 = F1*TABPHI2UP(IMU,IW0,ID) + F2*TABPHI2UP(IMU,IW0,ID+1) &
        + F3*TABPHI2UP(IMU+1,IW0,ID) + F4*TABPHI2UP(IMU+1,IW0,ID+1)
        PHI2UP1 = F1*TABPHI2UP(IMU,IW1,ID) + F2*TABPHI2UP(IMU,IW1,ID+1) &
        + F3*TABPHI2UP(IMU+1,IW1,ID) + F4*TABPHI2UP(IMU+1,IW1,ID+1)

        PHI2DN0 = F1*TABPHI2DN(IMU,IW0,ID) + F2*TABPHI2DN(IMU,IW0,ID+1) &
        + F3*TABPHI2DN(IMU+1,IW0,ID) + F4*TABPHI2DN(IMU+1,IW0,ID+1)
        PHI2DN1 = F1*TABPHI2DN(IMU,IW1,ID) + F2*TABPHI2DN(IMU,IW1,ID+1) &
        + F3*TABPHI2DN(IMU+1,IW1,ID) + F4*TABPHI2DN(IMU+1,IW1,ID+1)
    ENDIF

!         Linearly interpolate the scattering properties in wavenumber
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
    END SUBROUTINE INTERP_SCAT_TABLE3

!************************************************************************
    SUBROUTINE READ_ABS_PROFILE (ABSFILE, BINARYABSFILE, &
    ABSNU1, ABSNU2, ABSDELNU, &
    NABSNU,                        & !!!!!!!MAXNZ, MAXABSNU,
    NLEV, HEIGHT, TEMP, ABSPROF)
!       Reads in the gaseous absorption file which contains layer
!     optical depths for NLEV-1 layers and NABSNU wavenumbers.
!     The height (km) and temperature (K) of the profile are also returned.
!     The file may be in ascii text or binary format.

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    INTEGER :: NABSNU, NLEV            !!!!!!!!!!!MAXNZ, MAXABSNU
    LOGICAL :: BINARYABSFILE
    REAL*8 ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::    HEIGHT(*), TEMP(*), ABSPROF(MAXNZ,*)
    CHARACTER ABSFILE*(*)
    INTEGER :: I, J
    REAL ::    NU

    IF (BINARYABSFILE) THEN
        OPEN (UNIT=1,FILE=ABSFILE, STATUS='OLD',FORM='UNFORMATTED')
        READ (1) ABSNU1, ABSNU2, ABSDELNU, NABSNU
    !        WRITE (*,*) ABSNU1, ABSNU2, ABSDELNU, NABSNU
        IF (NABSNU > MAXABSNU) STOP 'RTSPEC: MAXABSNU exceeded'
        READ (1) NLEV
    !        WRITE (*,*) NLEV
        IF (NLEV > MAXNZ) STOP 'RTSPEC: MAXNZ exceeded'
        READ (1) (HEIGHT(I), I=1,NLEV)
    !        WRITE (*,'(52(1X,F5.1),:)') (HEIGHT(I), I=1,NLEV)
        READ (1) (TEMP(I), I=1,NLEV)
    !        WRITE (*,'(52(1X,F5.1),:)') (TEMP(I), I=1,NLEV)
        DO J = 1, NABSNU
            READ (1) (ABSPROF(I,J), I=1,NLEV-1)
        ENDDO
        CLOSE (1)
    ELSE
        OPEN (UNIT=1, FILE=ABSFILE, STATUS='OLD')
        READ (1,*)
        READ (1,*) ABSNU1, ABSNU2, ABSDELNU, NABSNU

        IF (NABSNU > MAXABSNU) THEN
            write(kStdErr,*) 'NABSNU',NABSNU,MAXABSNU
            write(kStdErr,*) 'RTSPEC: MAXABSNU exceeded'
            CALL DOStop
        ENDIF
        READ (1,*) NLEV
        IF (NLEV > MAXNZ) STOP 'RTSPEC: MAXNZ exceeded'
        READ (1,*) (HEIGHT(I), I=1,NLEV)
        READ (1,*) (TEMP(I), I=1,NLEV)
        READ (1,*)
        DO J = 1, NABSNU
            READ (1,*) NU, (ABSPROF(I,J), I=1,NLEV-1)
        ENDDO
        CLOSE (1)
    ENDIF
    RETURN
    END SUBROUTINE READ_ABS_PROFILE

!************************************************************************
! this subroutine does downward thermalrad tansfer from iS to iE
! ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
! THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
! and that raTemp has already been initialized with kTSpace Planck fcn

! this is the same as  SUBROUTINE FastBDRYL2GDiffusiveApprox_rtspec in
! rad_diff.f, except that it has been modified for rtspec.f applications,
! as 1..iNumLayer ---> iNumLayer .. 1

! for layers   20..1, it uses acos(3/5)
! for layers 100..20, it does t(i-1->0,x1)-t(i->0,x2)
!    where x1 is calculated at layer i-1, x2 is calculated at layer i

    SUBROUTINE FastBDRYL2GDiffusive_rts(TOA_to_instr,MU, WAVENO, &
    NLEV, TEMP, TAU, RAD0DN, RAD0DN0, ibdry)
!     Compute the downwelling radiance at the bottom of each cloud layer
!     using the Fast Diffusivity Approx
     
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
     
    INTEGER :: NLEV
    REAL ::    MU, WAVENO
    REAL ::    TEMP(NLEV), TAU(NLEV)
    REAL ::    RAD0DN(kProfLayer)
    REAL ::    RAD0DN0, TOA_to_instr
    INTEGER :: ibdry
     
    INTEGER :: iS,iE

! local variables
    INTEGER :: iLay,iL,iLp1,iBdryP1,iSecondEnd,iCase
    REAL :: rPlanck,rMPTemp,rFreqAngle,rFreqAngle_m1

! to do the angular integration
    REAL :: rAngleTr_m1,rAngleTr,rL2G,rL2Gm1
    REAL :: rDiff,rCosDiff

    iBdryP1 = NLev
    iCase   = -1

! note that we are not as careful as  FastBDRYL2GDiffusiveApprox in that we
! do not completely fill up atmospehere to include layers above instrument
! (if intrument is not at TOA)
    iS = nlev-1
    iE = 1
! now we have 3 different cases to consider
! CASE A1 : easy -- this is do ENTIRE atmnosphere
! iS=100   iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diffusive approx
! CASE A2 : easy -- this is do instr-gnd
! iS~50    iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diffusive approx
    IF ((iS >= iBdry) .AND. (iBdry >= iE)) THEN
        iCase   = 1
        iBdryP1 = iBdry+1
    END IF
! CASE B : quite easy -- this is do atmosphere -- instr
! iS=100   iE>iB                  ==> do iS->iE using acos(3/5)
    IF ((iS >= iBdry) .AND. (iBdry <= iE)) THEN
        iCase   = 2
        iBdryP1 = iE
    END IF
! CASE C : easy -- this is do instr-gnd
! iS~50    iE~1   iB > iS,iE      ==> do iB->iE using accurate diffusive approx
    IF ((iBdry >= iS) .AND. (iBdry >= iE)) THEN
        iCase = 3
        iBdry = iS
    END IF

    IF (iCase == -1) THEN
        write(kStdErr,*)'In FastBDRYL2GDiffusive_rts, icase = -1'
        CALL DoSTOP
    END IF

! ****** now map 1 .. iNumLayer ------> iNumlayer .. 1      *******
    iBdryP1 = nlev-iBdryP1+1
    iBdry   = nlev-iBdry+1
    iE      = nlev-1
    iS      = 1
! ******** also note iLm1 ---> iLp1
          
    rDiff    = (kThermalAngle*kPi/180.0)
    rCosDiff = cos(rDiff)

    if (TOA_to_instr < 0) THEN
        RAD0DN0 = ttorad(WAVENO,sngl(kTspace))
    else
        RAD0DN0 = TOA_to_instr
    END IF

!     now just go from TOA to instrument .. assume there are no clouds
!      RAD0DN0=RAD0DN0*exp(-TOA_to_instr)

    RAD0DN(1) = RAD0DN0

! initalize raL2G,raL2Gm1
    rL2G   = 0.0
    rL2Gm1 = 0.0

! calculate rL2Gm1 which is the L2G transmission from layer iS-1 to ground
    DO iLay = 2,nlev-1
        iL     = iLay
        rL2Gm1 = rL2Gm1+tau(iL)
    END DO
! calculate rL2G which is the L2G transmission from layer iS to ground
! and initialise the angles
    rL2G = rL2Gm1+tau(1)

! do top part of atmosphere, where we can use acos(3/5)
    IF ((iCase == 1)  .OR. (iCase == 2)) THEN
    ! go from top of atmosphere to boundary
        DO iLay = iS,iBdryp1
            iL      = iLay
            iLp1    = iLay+1
            rMPTemp = temp(iL)
        ! find the diffusive angles for the layer beneath
            rAngleTr_m1 = exp(-rL2Gm1/rCosDiff)
            rAngleTr    = exp(-rL2G/rCosDiff)
        ! Planckian emissions
            rPlanck = ttorad(waveno,rMPTemp)
            RAD0DN0 = RAD0DN0 + rPlanck*(rAngleTr_m1-rAngleTr)
            RAD0DN(1) = RAD0DN0
        ! get ready for the layer beneath
            rL2G   = rL2Gm1
            rL2Gm1 = rL2Gm1-tau(iLp1)
        END DO
    END IF

    IF ((iCase == 1) .OR. (iCase == 3)) THEN
    ! go from boundary to ground, or iE
    ! do bottom part of atmosphere ACCURATELY

        iSecondEnd=nlev-1
        rAngleTr   = FindDiffusiveAngleExp(rL2G)
        rFreqAngle = rAngleTr

        DO iLay = iBdry,iSecondEnd
            iL = iLay
            iLp1    = iLay+1
            rMPTemp = temp(iL)
        ! find the diffusive angles for the layer beneath
            rAngleTr_m1   = FindDiffusiveAngleExp(rL2Gm1)
            rFreqAngle_m1 = rAngleTr_m1
            rAngleTr_m1   = exp(-rL2Gm1/rAngleTr_m1)
            rAngleTr      = rFreqAngle
            rAngleTr      = exp(-rL2G/rAngleTr)
        ! Planckian emissions
            rPlanck = ttorad(waveno,rMPTemp)
            RAD0DN0 = RAD0DN0+rPlanck*(rAngleTr_m1-rAngleTr)
            RAD0DN(1) = RAD0DN0
        ! get ready for the layer beneath
            rL2G       = rL2Gm1
            rL2Gm1     = rL2Gm1-tau(iLp1)
            rFreqAngle = rFreqAngle_m1
        END DO

    END IF

! whether we did gaussian quadrature or diffusive approx, we now need the 2pi
! factor from the azimuthal integration
! however, there is also an average factor of 0.5 ==> overall, we need "pi"
    RAD0DN0 = RAD0DN0*kPi

    RETURN
    end SUBROUTINE FastBDRYL2GDiffusive_rts

!************************************************************************
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
! same as AddUppermostLayersQ, except it accepts raaAbs as input, and
! outputs radiance from TOA to instr ---- if instr is at TOA, it outputs -10

    SUBROUTINE  Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
!                  temporarily define atm from GND to TOP of atmosphere
! iT             = number of layers in this temporary atmosphere
! iExtra = -1 if no layeres added on, +1 if layers added on
! raExtra = array initialized to all zeros if instr at TOA
!         = array initialized to sum(k) from TOA to instr if instr inside atm
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    REAL :: raExtra(kMaxPts),rFracTop,raaAbs(kMaxPts,kMixFilRows)
    REAL :: raVTemp(kMixFilRows),raFreq(kMaxPts)
     
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iExtra
    INTEGER :: iI,iFr,iJ

    REAL :: waveno,rad,k,mudown
     
    iExtra=-1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        waveno=raFreq(iFr)
        raExtra(iFr) = ttorad(WAVENO,sngl(kTSpace))
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4))THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    !        write(kStdWarn,*) 'top most layer is fractional layer. Some'
    !        write(kStdWarn,*) 'portion needed above instrument to calculate'
    !        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
        !          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
    !        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
    !        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF

! ccccccccccc this is new .. where subroutine differs from AddUpperMostLayers
    IF (iExtra > 0) THEN
        MUDOWN=3.0/5.0
        DO iFr=1,kMaxPts
            waveno=raFreq(iFr)
            raExtra(iFr) = ttorad(WAVENO,sngl(kTSpace))
        END DO
        DO iI = iT,iNumLayer+1,-1
            iJ = iaRadLayerTemp(iI)
            DO iFr=1,kMaxPts
                waveno=raFreq(iFr)
                k=raaAbs(iFr,iJ)
                rad = ttorad(WAVENO,raVTemp(iJ))
                raExtra(iFr) = raExtra(iFr)*exp(-k/MUDOWN)+rad*(1-exp(-k/MUDOWN))
            END DO
        END DO

        DO iI = iNumLayer,iNumLayer
            iJ = iaRadLayerTemp(iI)
            DO iFr=1,kMaxPts
                waveno=raFreq(iFr)
                k=raaAbs(iFr,iJ)*(1-rFracTop)
                rad = ttorad(WAVENO,raVTemp(iJ))
                raExtra(iFr) = raExtra(iFr)*exp(-k/MUDOWN)+rad*(1-exp(-k/MUDOWN))
            END DO
        END DO
    ELSE
        write (kStdWarn,*) 'no need to add on any layers from TOA to intr'
    END IF

    RETURN
    END SUBROUTINE 

!************************************************************************
! this subroutine interpolates the flux, based on clear sky rad transfer
! ie we call DISORT or RTSPEC at point spacing = iStep
! so we have to interpolate this slow flux variation onto the fast spectral
! variation which is evident in the monochromatics radiance
    SUBROUTINE InterpolateFlux(raaFlux,iLay,raKC,raFreq,iStep)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input is raaFlux, layer iLay, computed on point spacing = kDis_Pts
! output is raaFlux, layer iLay, computed on point spacing = 1
    REAL :: raaFlux(kMaxPts,kProfLayer+1),raKC(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iLay,iSTep

! local variables
    REAL :: y1,y2,x1,x2,m,c,sf
    INTEGER :: iJ,iI,jNU1,jNU2,iFr,iFr1,iFr2

    jNU1 = 1
    jNU2 = kMaxPts
    iJ = 0

    DO iI = jNU1, jNU2, iStep
        iJ = iJ + 1
        iFr1 = jNU1 + iStep*(iJ-1)
        iFr2 = iFr1 + iStep

        IF (iFr2 > kMaxPts) GOTO 210

    ! ompute the scale factors
        y1 = raaFlux(iFr1,iLay)/raKC(iFr1)
        y2 = raaFlux(iFr2,iLay)/raKC(iFr2)

        x2 = raFreq(iFr2)
        x1 = raFreq(iFr1)

    ! ompute the slope and intercept
        m = (y2 - y1)/(x2 - x1)
        c = y2 - m*x2

    ! ow interpolate raKC linearly, with this scale factor!!!
        DO iFr = iFr1,iFr2
            sf = m*(raFreq(iFr)-x1) + y1
            raaFlux(iFr,iLay) = sf*raKC(iFr)
        END DO
    END DO


    210 CONTINUE
! do the last chunk, upto the 10000th point
    iFr2 = kMaxPts
! ompute the scale factors
    y1 = raaFlux(iFr1,iLay)/raKC(iFr1)
    y2 = raaFlux(iFr2,iLay)/raKC(iFr2)
    x2 = raFreq(iFr2)
    x1 = raFreq(iFr1)
! ompute the slope and intercept
    m = (y2 - y1)/(x2 - x1)
    c = y2 - m*x2
! ow interpolate linearly!!!
    DO iFr = iFr1,iFr2
        sf = m*(raFreq(iFr)-x1) + y1
        raaFlux(iFr,iLay) = sf*raKC(iFr)
    END DO
     
    RETURN
    end SUBROUTINE InterpolateFlux
!************************************************************************
! this subroutine sets up the scattering table info from SSCATMIE.F
    SUBROUTINE SetMieTables_RTSPEC(raFreq, & 
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
    iSergio, &
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)
     
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
     
! iSergio INTEGER that tells if this is RTSPEC or SERGIO's code
    INTEGER :: iSergio
    INTEGER :: iaCldTypes(kMaxClouds)  !! 101 201 or 301 for water, ice, aerosol
    REAL :: raFreq(kMaxPts)
! ---------------- inputs needed to read scattering tables -------------------
! this is which atm number is being used, and whether these are binary files
    INTEGER :: iAtm,iBinaryFile,iNumLayer,iDownward
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this is just to set everything about clouds relative to TOA layer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
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
          
! c      INTEGER  MAXTAB, MAXGRID, MAXSCAT
! c      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
    CHARACTER(120) :: SCATFILE(MAXSCAT)
     
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
     
    INTEGER :: NSCATTAB, NCLDLAY, NLEV, NABSNU
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    REAL ::    IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
     
    INTEGER :: iCloudySky        !!!!are there clouds in this atm???

! ---------------------------- local variables ----------------------------
    INTEGER :: IACLDTOPKCARTA(kMaxClouds), IACLDBOTKCARTA(kMaxClouds)
    INTEGER :: iaTable(kMaxClouds*kCloudLayers),iIn,iJ,iReadTable,I
    INTEGER :: iCloud,iStep
    REAL :: extinct
    INTEGER :: LL,II,N,M,iLayers,iBlah

    INTEGER :: iL,iT,iB,iNumClds,iT_Atm,iB_Atm,iaCldInLayer(kProfLayer)

    REAL ::    raCldLayer(MAXNZ),iwp0(maxnz),dme0(maxnz), rDmePhase,rScat
    REAL ::    dmetab_phase(kProfLayer)
    INTEGER :: indx(MAXNZ),iscattab0(maxnz),iiDiv,i1
     
    CHARACTER(120) :: caName
    CHARACTER(1) :: caScale(MAXSCAT)
     
! nitialise all scattering info to null

    iB = -kProfLayer
    iT = kMixFilRows+1

    iiDiv = 0
    555 CONTINUE
    IF (iiDiv*kProfLayer < iaaRadLayer(iAtm,3)) THEN
        iiDiv = iiDiv + 1
    END IF
    iiDiv = iiDiv - 1

! copied from s_scatter_spectra.f .. all table names etc are unique, so no
! need to make more checks

    iCloudySky = -1        !!!!!!!assume no clouds in sky

!!!! need to reference the cloud tops and bottoms wrt TOP layer of
!!!! defined atm
!!!! eg if atm from 971 to 150 mb (plane at 150 mb) ==>
!!!!       this occupies kCARTA layers 19-78
!!!!   if cloud from 248 to 214 mb  ==>
!!!!       this occupies kCARTA layers 69 to 72
!!!! Thus the cloud occupies RTSPEC atmosphere "tau" from 6 to 9

!!!!!!!!this is all that is needed if only RTSPEC rad transfer were used
    IF (iDownWard == 1) THEN
        iB_Atm = iaaRadLayer(iAtm,1)
        iT_Atm = iaaRadLayer(iAtm,iNumLayer)
    ELSEIF (iDownWard == -1) THEN
        iT_Atm = iaaRadLayer(iAtm,1)
        iB_Atm = iaaRadLayer(iAtm,iNumLayer)
    END IF
!!!!!!hwowever we also do fluxes, so even if the atm is defined so it
!!!!!!is for an uplook instrument, RTSPEC will be called in a downlook
!!!!!!fashion, and vice versa

    IF (iB_Atm > iT_Atm) THEN
        iCloudySky = iT_Atm
        iT_Atm = iB_Atm
        iB_atm = iCloudySky
    END IF
          
! nitialise all scattering info to null

    iCloudySky = -1        !!!!!!!assume no clouds associated with this atm

    iCldTopKcarta = -1
    iCldBotKcarta = kProfLayer+1
    NSCATTAB=-1000
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
        iaCldTopkCarta(iIn) = -1
        iaCldBotkCarta(iIn) = -1
    END DO

    DO iIn=1,iNclouds
        DO iJ=1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            IF (iI > MAXSCAT) THEN
                write(kStdErr,*)'unfortunately, in scatterparam.f90 we have '
                write(kStdErr,*)'MAXSCAT = ',maxscat
                write(kStdErr,*)'please reset and retry'
                CALL DoSTOP
            END IF
            caName=caaaScatTable(iIn,iJ)
            IF (iaTable(iI) < 0) THEN  !nothing associated with this yet
                IF (iI > NSCATTAB) THEN
                    NSCATTAB = iI
                END IF
                iaTable(iI) = 1
                ScatFile(iI) = caName
            END IF
        END DO

    ! check to see if this cloud is to be used with this atm . recall that
    ! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
    ! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    ! INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
        DO iJ=1,iaCloudNumAtm(iIn)
            IF (iaaCloudWhichAtm(iIn,iJ)  == iAtm) THEN
                iCloudySky = iIn         !!!! set this up
                iaCloudWithThisAtm(iIn) = 1
                IACLDTOP(iIn) = iaaCloudWhichLayers(iIn,1)+1
                IACLDBOT(iIn) = iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn))
                iaCldTopkCarta(iIn) = iaCldTop(iIn)     !!! not needed
                iaCldBotkCarta(iIn) = iaCldBot(iIn)     !!! not needed
            ! iCldTopkCarta = iaCldTop(iIn)-1
            ! iCldBotkCarta = iaCldBot(iIn)
                IF (iCldTopkCarta < iaCldTop(iIn)-1) THEN
                    iCldTopkCarta = iaCldTop(iIn)-1
                END IF
                IF (iCldBotkCarta > iaCldBot(iIn)) THEN
                    iCldBotkCarta = iaCldBot(iIn)
                END IF
                write(kStdWarn,*) ' '
                write(kStdWarn,*)'cloud # ',iIn,' associated with atm # ',iAtm
                write(kStdWarn,*)'setmie0 : cloud is in KCARTA layers ', &
                iiDiv*kProfLayer+iaCldTop(iIn)-1,' to ', &
                iiDiv*kProfLayer+iaCldBot(iIn)

            !!!!!these are the RTSPEC layers 100 to 1 = GND to TOA
                iaCldbot(iIn) = iT_Atm - iaCldbot(iIn) + 1
                iaCldtop(iIn) = iT_Atm - iaCldtop(iIn) + 1
                write(kStdWarn,*)'setmie0 : cloud is in RTSPEC layers ', &
                iaCldTop(iIn)+1,' to ',iaCldBot(iIn)

            END IF
        END DO
         
    ! check to see which scattering tables to be used with this atm
        DO iJ=1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            IF (iaCloudWithThisAtm(iIn) == 1) THEN
                iaScatTable_With_Atm(iI) = 1
                write(kStdWarn,*)'scat table ',iI,' for atm,layer # ',iAtm,iJ
            END IF
        END DO
    END DO      !!!!!!!!main       DO iIn=1,iNclouds

!     Only read in scattering tables that are needed for this atm
    iReadTable = 1
          
    IF (iReadTable > 0) THEN
        IF (iBinaryFile == 1) THEN
            DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading binary scatter data for table #',I
                    write(kStdWarn,*) scatfile(I)
                    CALL READ_SSCATTAB_BINARY(SCATFILE(I),   & !!!!!!MAXTAB, MAXGRID,
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))
                    IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR. &
                    (ABS(MUINC(2)-0.7887) > 0.001)) THEN
                        write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                        CALL DoStop
                    END IF
                    IF (iaPhase(I) > 0) THEN
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        ELSE IF (iBinaryFile == -1) THEN
            DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading ascii scatter data for table #',I,SCATFILE(I)
                    CALL READ_SSCATTAB(SCATFILE(I),   & !!!!!!MAXTAB, MAXGRID,
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))

                    IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR. &
                    (ABS(MUINC(2)-0.7887) > 0.001)) THEN
                        write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                        CALL DoStop
                    END IF
                    IF (iaPhase(I) > 0) THEN
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        ELSE IF (iBinaryFile == 0) THEN
            DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading "special" ascii scatter data for table #',I
                    CALL READ_SSCATTAB_SPECIAL(SCATFILE(I), &
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))

                    IF (iaPhase(I) > 0) THEN
                        write(kStdErr,*) 'Right now, incapapable of this silly task!!!'
                        write(kStdErr,*) 'need iaPhase = 0 for iBinaryFIle = 0'
                        CALL DoStop
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        END IF    !iBinaryFile > 0
    END IF      !iReadTable  > 0

! Frank Evans code scales the Mie scattering parameters, so if we are using
! my canned EDDINGTON method, we have to unscale them!!!!!!!!
    IF (iSergio > 0) THEN
        DO I = 1, NSCATTAB
            IF (iaScatTable_With_Atm(I) > 0) THEN
                CALL UnScaleMie( &
                caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                ndme(i)*nwavetab(i))
            END IF
        END DO
    END IF

    DO i=1,MAXNZ
        raCldLayer(I) = +1.0e10
        indx(I)       = -1
        iscattab0(I)  = -1
        dme0(I)       = -1.0
        iwp0(I)       = -1.0
    END DO

    iCloud = -1
    IF (iCloudySky < 0) THEN
    !!!!!this is similar to DISORT interface
        write(kStdWarn,*)'Could not find a cloud for atmosphere #',iAtm
        write(kStdWarn,*)'setting IWP = -100.0'
        iCloud=1    !say cloud number one associated with this atmosphere
        ncldlay=1   !say fictitious cloud occupies one layer
        IWP(1)      = -100.0   !but ensure cloud has NO particles in it!
        DME(1)      = -10.0    !but ensure cloud has NO particles in it!
        ISCATTAB(1) = -1
         
    ELSE
    !!!!!find total number of clouds, and hence total number of layers
    !!!!!this is quite different to DISORT interface, as we have to make
    !!!!!sure that cloud layers (if more than one cloud) are sequential
        NCLDLAY=0
        iLayers = 0  !!!!!  ----- "iLayers" is an important variable -----
        iNumClds = 0
        DO i=1,kMaxClouds
            IF (iaCloudWithThisAtm(i) == 1) THEN
                iNumClds = iNumClds + 1
                iT = i
                ncldlay = ncldlay + iaCloudNumLayers(i)
                write(kStdWarn,*) 'Cloud #, num layers = ',i,iaCloudNumLayers(i)
                write(kStdWarn,*) 'L  KCLay iwp  dme  iscattab             : '
                write(kStdWarn,*) '-----------------------------------------'
                DO iStep = 1, iaCloudNumLayers(i)
                    iLayers = iLayers+1
                    IWP(iLayers) = raaaCloudParams(i,iStep,1)
                    DME(iLayers) = raaaCloudParams(i,iStep,2)
                    ISCATTAB(iLayers) = iaaScatTable(i,iStep)
                    raCldLayer(iLayers) = +1.0 * iaaCloudWhichLayers(i,iStep)
                    write(kStdWarn,*) iLayers,int(raCldLayer(iLayers)),iwp(iLayers),dme(iLayers), &
                    iscattab(iLayers)
                    IWP0(iLayers) = raaaCloudParams(i,iStep,1)
                    DME0(iLayers) = raaaCloudParams(i,iStep,2)
                    ISCATTAB0(iLayers) = iaaScatTable(i,iStep)

                    raCldLayer(iLayers) = +1.0 * &
                    (iT_atm - iaaCloudWhichLayers(i,iStep) + 1)
                !                print *,'cloud , DISORT layer = ',i,raCldLayer(iLayers)
                END DO
            END IF
        END DO

    !!!this is where we totally differ from DISORT, as we have to fill
    !!!in any "in-between" layers with a fictitious empty cloud
        IF (iNumClds == 1) THEN
        ! iT is the highest cloud, and lowest cloud, set from above
        ! so just figure out top and bottom
            iB = iaaCloudWhichLayers(iT,iaCloudNumLayers(iT))
            iT = iaaCloudWhichLayers(iT,1)

            iB = iT_Atm - iB + 1
            iT = iT_Atm - iT + 1
            write (kStdWarn,*) 'RTSPEC cloud layers are from ',iB,' to ',iT

        ELSE

        ! oh boy, have to figure out if the clouds are next to each
        ! other; if not, create a blank cloud in between them
            DO i=1,kProfLayer
            !!!assume all layers are clear sky
                iaCldInLayer(i) = 0
            END DO
                      
            DO i=1,kMaxClouds  !!! see which layers have a cloud in them
                IF (iaCloudWithThisAtm(i) == 1) THEN
                    iT = iaaCloudWhichLayers(i,1)
                    iB = iaaCloudWhichLayers(i,iaCloudNumLayers(i))

                    iB = iT_Atm - iB + 1
                    iT = iT_Atm - iT + 1
                    write (kStdWarn,*) 'cloud # ',i,' : RTSPEC layers are from ', &
                    iB,' to ',iT

                    DO iL = iT,iB
                        iaCldInLayer(iL) = iaCldInLayer(iL) + 1
                    END DO
                END IF
            END DO

            iB = -kProfLayer
            iT = kMixFilRows+1
            DO i=1,kProfLayer
            ! see if more than one cloud per layer
                IF (iaCldInLayer(i) > 1) THEN
                    write(kStdErr,*) 'More than one cloud in kLAYERS layer ',i
                    write(kStdErr,*) 'Please check section SCATTR and retry'
                    CALL DoStop
                END IF
            ! see lowest, highest parts of different clouds simultaneously
                IF (iaCldInLayer(i) == 1) THEN
                    IF (i >= iB) iB = i
                    IF (i <= iT) iT = i
                END IF
            END DO

            write (kStdWarn,*) 'highest/lowest RTSPEC cloud layers = ',iT,iB

        !!!!!! now loop from iB to iT and see if we have to fill in
            DO i = iT,iB
                IF (iaCldInLayer(i) == 0) THEN
                    write(kStdWarn,999) i
                    ncldlay = ncldlay + 1
                    iLayers = iLayers + 1     !!! --- iLayers is VERY IMPORTANT ---
                    IWP(iLayers) = 0.0
                    DME(iLayers) = raaaCloudParams(iCloudySky,1,2)
                    ISCATTAB(iLayers) = iaaScatTable(iCloudySky,1)

                    IWP0(iLayers) = 0.0
                    DME0(iLayers) = raaaCloudParams(iCloudySky,1,2)
                    ISCATTAB0(iLayers) = iaaScatTable(iCloudySky,1)
                    raCldLayer(iLayers) = 1.0 * i

                END IF
            END DO

        !            DO i = 1,iLayers
        !              iI = i
        !              write(kStdWarn,*) 'before ',iI,iwp(iI),dme(iI),iscattab(iI),
        !     $                                    raCldLayer(iI)
        !            END DO

        !!!!!!!!!now sort the layers
            CALL NumericalRecipesIndexer(indx,raCldLayer,iLayers)
            DO i = 1,iLayers
                iwp(i) = iwp0(indx(i))
                dme(i) = dme0(indx(i))
                iscattab(i) = iscattab0(indx(i))
            !              write(kStdWarn,*) 'after ',i,iwp(i),dme(i),iscattab(i),
            !     $                                    raCldLayer(indx(i))
            END DO
        END IF       !!!IF (iNumClds == 1) THEN
    END IF

    999 FORMAT('empty RTSPEC layer ',I3,' found between clouds; set IWP = 0.0')

! not needed
!      WRITE(*,*) 'Observation level (km)'
!      READ (*,*) ZOBS
!      WRITE (*,*) 'Input cloud top height (km)'
!      READ (*,*) ZTOP

!     Find the levels for top of cloud and observation level
!     remember that these numbers are with respect to the KLAYERS pressure
!     levels and layers
!     these will be reset when they are passed in and out of GetAbsProfile
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

! when we call GetAbsProfile, the variables icldtop,icldbot,iobs will be reset
! to reflect the rtspec layering
!   TOA    --------------
!             layer 1
!          --------------
!              .....                 RTSPEC
!          --------------
!        layer iNumLayer-1
!          --------------
!         layer iNumLayer
!   GND --------------------------

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do

    IF (rScat > 0.0) THEN
        ICLDTOP = iT-1
        ICLDBOT = iB

        IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 3)) THEN
            ICLDTOP = icldtop+1
            ICLDBOT = icldbot+1
        END IF

        IF (iDownWard > 0) THEN
            IOBS    = iNumLayer
        ELSE IF (iDownWard < 0) THEN
            IOBS   = 1
        END IF
    END IF

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do
           
    IF (rScat <= 0.0) THEN  !we have no cloud; set up fictitious clouds
    !!!! this is never used by radiative transfer code alone, as a noncloud
    !!!! situation is dealt with by using procedure clearskyradtrans
    !!!! however,  when computing fluxes, this becomes important

        IF (iDownWard > 0) THEN
        ! own look instr : set cloud BELOW observer, in kCARTA layer #1
            ICLDTOP = 2
            ICLDBOT = 1
        ! own look instr : set cloud BELOW observer, RTSPEC layer #iNumLayer
            ICLDTOP = iNumLayer-2
            ICLDBOT = iNumLayer-1
            IOBS    = iNumLayer
        ELSE IF (iDownWard < 0) THEN    !up look instr
        ! p look instr : set cloud ABOVE observer, in kCARTA layer #iNumLayer
            ICLDTOP = iNumLayer+1
            ICLDBOT = iNumLayer
        ! p look instr : set cloud ABOVE observer, in RTSPEC layer #1
            ICLDTOP = 1
            ICLDBOT = 2
            IOBS    = 1
        END IF
    END IF

    RETURN
    end SUBROUTINE SetMieTables_RTSPEC

!************************************************************************
! this subroutine takes in the input abs coeffs (raaAbsCoeff) where raa(1,:)
! is the lowest layer and raa(kProfLayer,:) is the highest layer .. it then
! outputs these  abs coeffs intp absprof, where absprof(1,:) is the top, and
! absprof(iNumLayer,:) = ground

!     sets optical depths for NLEV-1 layers and NABSNU wavenumbers.
!     The temperature (K) of the profile is also returned. (no height needed)

    SUBROUTINE GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
    ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raLayerTemp, &
    iProfileLayers, raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are variables that come in from kcartamain.f
    REAL :: raaAbs(kMaxPts,kMixFilRows),raFreq(kMaxPts),rFracTop,rFracBot
    REAL :: raVTemp(kMixFilRows),rSurfaceTemp
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNpmix
    INTEGER :: iDownWard,iProfileLayers
    REAL :: iwp,raPressLevels(kProfLayer+1),rSurfPress
! these are variables that we have to set
    INTEGER ::  NABSNU, NLEV         !!!!!!!!!!!!!MAXNZ, MAXABSNU
    REAL ::   ABSNU1, ABSNU2, ABSDELNU
    REAL ::   TEMP(*), ABSPROF(MAXNZ,*), raLayerTemp(*)
    INTEGER ::  ICLDTOP,iCLDBOT,IOBS

! local variables
    INTEGER :: iaRadLayer(kProfLayer), iFr, iL, iLay
    REAL :: NU, raVT1(kMixFilRows), InterpTempSurf,rT2

! these are to flip the temperature, abs profiles if instr looks up
    REAL :: raTemp(kProfLayer+1),raaTempAbs(kProfLayer,kMaxPts)

    absnu1=raFreq(1)
    absnu2=raFreq(kMaxPts)
    absdelnu=(absnu2-absnu1)/(kMaxPts-1)
    nabsnu = kMaxPts
    nlev = iNumLayer+1           !this is the number of pressure levels

    DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
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

! set the temperature of the bottommost layer correctly
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO

    IF (iDownWard == 1) THEN         !downlook instr
    ! if the bottommost layer is fractional, interpolate!!!!!!
        iL = iaRadLayer(1)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot, &
        &                        1,iL)
        write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
    !        rT2 = interpTempSurf(iProfileLayers,raPressLevels,raVTemp,rFracBot,
    !     $                       1,iL,rSurfaceTemp,rSurfPress)
    !        write(kStdWarn,*) 'surface bottom temp interped to ',rT2

    ! if the topmost layer is fractional, interpolate!!!!!!
    ! this is hardly going to affect thermal/solar contributions (using this temp
    ! instead of temp of full layer at 100 km height!!!!!!
        iL = iaRadLayer(iNumLayer)
        raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop, &
        -1,iL)
        write(kStdWarn,*) 'top temp interped to ',raVT1(iL)
    ELSE IF (iDownWard == -1) THEN       !uplook instr
    ! if the bottom layer is fractional, interpolate!!!!!!
        iL = iaRadLayer(iNumLayer)
        raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot, &
        &                        1,iL)
        write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
    ! if the top layer is fractional, interpolate!!!!!!
        iL = iaRadLayer(1)
        raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop, &
        -1,iL)
        write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)
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

! set the vertical temperatures of the atmosphere
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer, &
    iDownWard,iProfileLayers,raPressLevels)

! now set up the abs coeffs
! initialize array to all zeroes
    DO iFr=1,kMaxPts
        DO iLay = iNumLayer+1,kProfLayer
            absprof(iLay,iFr) = 0.0
        END DO
    END DO

    IF (iDownWard == 1) THEN        !!!! no problemo, down look instr
        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            nu=1.0
            IF (iLay == 1) THEN
                nu=rFracBot
            ELSE IF (iLay == iNumLayer) THEN
                nu=rFracTop
            END IF
            DO iFr=1,kMaxPts
            ! bsprof wants level 1 == TOA, level iNumLayer= gnd
                absprof(iNumLayer-iLay+1,iFr) = raaAbs(iFr,iL)*nu
            END DO
        END DO
    ELSEIF (iDownWard == -1) THEN        !!!! oopsy, up look instr
        DO iLay=1,iNumLayer
            iL = iaRadLayer(iLay)
            nu=1.0
            IF (iLay == iNumLayer) THEN
                nu=rFracBot
            ELSE IF (iLay == 1) THEN
                nu=rFracTop
            END IF
            DO iFr=1,kMaxPts
            ! bsprof wants level 1 == TOA, level iNumLayer= gnd
                absprof(iNumLayer-iLay+1,iFr) = raaAbs(iFr,iL)*nu
            END DO
        END DO
    END IF

! now set icldtop icldbot, iobs
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
! remember there is ONE more level than there are layers
!      icldtop=(iNumLayer+1)-icldtop+1
!      icldbot=(iNumLayer+1)-icldbot+1
!      iobs=(iNumLayer+1)-iobs+1

! icldtop,icldbot are set in the main calling routine

! 1      if (iDownWard .gt. 0) then
!        iobs=1
!        icldtop = iNumLayer-icldtop+1
!        icldbot = iNumLayer-icldbot+1
!      else if (iDownWard.lt. 0) then
!        iobs = iNumLayer
!        icldtop = iNumLayer-icldtop+1
!        icldbot = iNumLayer-icldbot+1

    IF  (iDownWard == -1) THEN
    ! lip TEMPerature array
    ! emember there is one more level than there are layers
        DO  iLay=1,iNumLayer+1
            raTemp(iLay) = TEMP((iNumLayer+1)-iLay+1)
        END DO
        DO  iLay=1,iNumLayer+1
            TEMP(iLay) = raTemp(iLay)
        END DO
    ! lip absprof array
        DO iFr=1,kMaxPts
            DO  iLay=1,iNumLayer
                raaTempAbs(iLay,iFr) = absprof(iNumLayer-iLay+1,iFr)
            END DO
        END DO
        DO iFr=1,kMaxPts
            DO  iLay=1,iNumLayer
                absprof(iLay,iFr) = raaTempAbs(iLay,iFr)
            END DO
        END DO
    END IF      !!!!!!!if iDownWard == -1

    RETURN
    end SUBROUTINE GetAbsProfileRTSPEC

!************************************************************************
!       this is for scatter_disort
!************************************************************************
! this subroutine calculates the solar beam incident at TOA
    SUBROUTINE SolarBeamDisort(iDoSolar,raSun,raFreq,iTag)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! raSun    = final solar contr
! raFreq  = frequency array
    REAL :: raSun(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iTag,iDoSolar
! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))
     
! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: rCos,raKabs(kMaxPts)
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
    REAL :: rExtraFac

    rExtraFac = 1.0     !! before 01/17/06
    rExtraFac = kPi     !! for a few days after 01/17/06
    rExtraFac = 1.0     !! after 01/27/06
              
    IF (iDoSolar == 0) THEN
    ! se 5700K
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
        ! compute the Plank radiation from the sun
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
    ELSEIF (iDoSolar == 1) THEN
    ! ead in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
    END IF

! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun

! now account for the flux hitting the earth at an angle (use radians)
! instead of rSunAngle, use angle at lowest layer
!      rSunAngle = kSolarAngle
!      rSunAngle = (rSunAngle*kPi/180.0)
!      rCos      = cos(rSunAngle)
           
! kCARTA non scatter would adjust raSun by cos(rSunAngle) * rSolidAngle
! but DISORT does the adjust raSun by cos(rSunAngle), so just use rSolidAngle
! Jan 17 2006 : since this supposed to be solar beam flux, multiply by kPi
!   this should now be consistent with kOmegaSun in pre_defined.param
! Jan 27 2006 : went back to rExtraFac = 1.0, as we adjusted things for
!    kTwoStream and PCLSAM

! 6.785087652174316e-5 = pi(sun diam/sun dist)^2 = pi*((0.6951e9/149.57e9)^2)
!        rOmegaSun = 6.785087652174316e-5      version on Jan 2006

    DO iFr = 1,kMaxPts
        raSun(iFr) = raSun(iFr) * rOmegaSun * rExtraFac
    END DO

    RETURN
    end SUBROUTINE SolarBeamDisort
      
!************************************************************************
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
! same as AddUppermostLayersQ, except it accepts raaAbs as input, and
! outputs radiance from TOA to instr ---- if instr is at TOA, it outputs -10
! same as Find_Radiance_TOA_to_instr (in scatter_rtspec), except it outputs
! the total optical depth bewteen TOA and instrument

    SUBROUTINE Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp,rFracTop, &
    raFreq,raaAbs,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
     
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
!                  temporarily define atm from GND to TOP of atmosphere
! iT             = number of layers in this temporary atmosphere
! iExtra = -1 if no layeres added on, +1 if layers added on
! raExtra = array initialized to all zeros if instr at TOA
!         = array initialized to sum(k) from TOA to instr if instr inside atm
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    REAL :: raExtra(kMaxPts),rFracTop,raaAbs(kMaxPts,kMixFilRows)
    REAL :: raVTemp(kMixFilRows),raFreq(kMaxPts)
     
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iExtra
    INTEGER :: iI,iFr,iJ

    REAL :: waveno,rad,k,mudown
     
    iExtra=-1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4))THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    !        write(kStdWarn,*) 'top most layer is fractional layer. Some'
    !        write(kStdWarn,*) 'portion needed above instrument to calculate'
    !        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
        !          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
    !        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
    !        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF


! c this is new .. where subroutine differs from AddUpperMostLayers
! c this is new .. where subroutine differs from Find_Radiance_TOA_to_Instr

    if (iExtra > 0) THEN
        DO iI = iT,iNumLayer+1,-1
            iJ = iaRadLayerTemp(iI)
            DO iFr=1,kMaxPts
                raExtra(iFr) = raExtra(iFr) +  raaAbs(iFr,iJ)
            END DO
        END DO

        DO iI = iNumLayer,iNumLayer
            iJ = iaRadLayerTemp(iI)
            DO iFr=1,kMaxPts
                raExtra(iFr) = raExtra(iFr) +  raaAbs(iFr,iJ)*(1-rFracTop)
            END DO
        END DO

    END IF

    RETURN
    end SUBROUTINE Find_K_TOA_to_instr

!************************************************************************
! if needed, unscale parameters so DISORT is happy, when compared to RTSPEC
! Frank Evans recommended not to do this, but instead turn DELTAM off in DISORT
! however, leaving DELTAM on does not significantly change things
    SUBROUTINE UnScaleMie(cScale,TABEXTINCT, TABSSALB, TABASYM, N)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL :: tabextinct(*),tabssalb(*),tabasym(*)
    INTEGER :: N
    CHARACTER cScale

    INTEGER :: iI
    REAL :: wn,w0,w1,an,a0_p,a0_m,a0,a1,en,e0,e1,f

! Frank Evans code scales the Mie scattering parameters, so we have to
! unscale them!!!!!!!!
!       Delta function scale the scattering properties in the phase function.
!     The delta function fraction (F) is related to the asymmetry parameter,
!     depending on DELTASCALE: N - F=0, Y - F=g, H - F=g^2, G - F=g^2 w/
!     Gaussian filtering
!       The extinction and single scattering albedo are scaled but the
!     phase function Legendre coefficients are not, because the phase
!     function scaling is done in CALC_PHI.  Instead the scaling fraction
!     is returned in LEGEG(0,*) = 1-F.  The scaled asymmetry parameter is
!     computed and returned.
!       Subroutine modified 12/11/96 to force Gaussian filtering of Legendre
!     coefficients.  Gaussian-width values L0 determined by fitting
!     half-width half-max of forward scattering peak (in terms of mu) of
!     actual phase function to a Gaussian-filtered delta function - phase
!     function with width parameter L0 ;
!       Subroutine FINDL0 determines L0, subroutine FINDHALF determines HWHM
!     of arbitrary phase function, subroutine SCATCALC calculates phase
!     function at specific value of mu
!      REAL EXTINCT, ALBEDO, ASYM
!      CHARACTER DELTASCALE
!      REAL    F, FCTR, L0
     
!      ASYM = LEGEN(1)/3.0
!      IF (DELTASCALE .EQ. 'Y') THEN
!        F = ASYM
!      ELSE IF ((DELTASCALE .EQ. 'H').OR.(DELTASCALE .EQ. 'G')) THEN
!        F = ASYM**2
!      ELSE
!        F = 0.0
!      ENDIF

!   Scale the extinction, single scattering albedo,  asymmetry parameter
!      FCTR = (1 - F*ALBEDO)
!      EXTINCT = EXTINCT*FCTR
!      ALBEDO = ALBEDO*(1-F)/FCTR
!      ASYM = (ASYM-F)/(1-F)
!      IF (DELTASCALE .EQ. 'Y' .OR. DELTASCALE .EQ. 'H') THEN
!        LEGEN(0) = 1-F
!      ENDIF
     
    IF ((cScale == 'n') .OR. (cScale == 'N')) THEN
    !!!nothing scaled, so everything is cool
        GOTO 10
    END IF
    IF ((cScale == 'y') .OR. (cScale == 'Y')) THEN
        write (kStdErr,*) 'Cannot invert the Mie parameters for ''y'' scaling'
        write (kStdErr,*) 'Please rerun sscatmie using n,g or h scaling'
        CALL DoStop
    END IF

! only cases left are "g" or "h" scaling, so go ahead
! x1 is the new scaled stuff that sscatmie put out
! x0 is the original unscaled stuff that DISORT wants
! xn is a check, to see if using "f" we get a0 -> a1, w0 -> w1, e0 ->e1
    DO iI = 1,N
        a1 = tabasym(iI)
        a0 = 1 + 4*a1*a1 - 4*a1            !!!!!this term is always positive
        a0_m = (-1 - sqrt(a0))/(2*(a1-1))  !!!!!a1-1 always < 0
        a0_p = (-1 + sqrt(a0))/(2*(a1-1))  !!!!!a1-1 always < 0
        a0 = min(a0_p,a0_m)
        tabasym(iI) = a0
        f = a0*a0                      !!!this is the scaling parameter

        w1 = tabssalb(iI)
        w0 = w1/(1+f*(w1-1))
        tabssalb(iI) = w0

        e1 = tabextinct(iI)
        e0 = e1/(1-f*w0)
        tabextinct(iI) = e0

    !!!!just a check
    ! n=e0*(1-f*w0)
    ! n=w0*(1-f)/(1-f*w0)
    ! n=(a0-f)/(1-f)
    ! rint *,a0,a1,an,e0,e1,en,w0,w1,wn

    END DO
            
    10 CONTINUE
    RETURN
    end SUBROUTINE UnScaleMie
!************************************************************************
! this integer function finds out WHERE to put the cloud layer WRT kCARTA
! Suppose RTSPEC,DISOR have iNumlayers in Atm, from 1 to iNumLayer
!                      with cloud in layers iC1 ..iCN
! then RTSPEC GND = layer iNumLayer -----> KCARTA layer iaRadLayer(1)
! then RTSPEC TOA = layer 1         -----> KCARTA layer iaRadLayer(iNumLayer)
    INTEGER FUNCTION iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,iL)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    INTEGER :: iAtm,iNumLayer,iL,iUpDown,iS,iE,iI,iiDiv

! his is for downlook instr

    iS = iaaRadLayer(iAtm,1)
    iiDiv = 0
    10 CONTINUE
    IF (iiDiv * kProfLayer < iS) THEN
        iiDiv = iiDiv + 1
        GOTO 10
    END IF
    iiDiv = iiDiv - 1

    iE = iaaRadLayer(iAtm,iNumLayer)
    IF (iS <= iE) THEN
        iUpDown = +1             !!!down look instr
    ELSE
        iUpDown = -1             !!!up look instr
    END IF

    IF (iUpDown > 0) THEN
        iI = iaaRadLayer(iAtm,iNumLayer) - iL + 1
    ELSE
        iI = iaaRadLayer(iAtm,1) - iL + 1
    END IF
    iFindWhereInAtm = iI + (kProfLayer*iiDiv)

    RETURN
    end FUNCTION iFindWhereInAtm
!************************************************************************
! this subroutine copies the relevant parts of raaAbs to raaExtTemp
! right now it hust makes a copy of the atmospheric gases absorbtion
    SUBROUTINE CopyRaaExt_twostream(raaAbs,raaExtTemp,raaScatTemp, &
    raaAsymTemp,iaaRadLayer,iAtm,iNumlayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
    REAL :: raaAbs(kMaxPts,kMixFilRows)        !original, from uncompression
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !temporary copy
    REAL :: raaScatTemp(kMaxPts,kMixFilRows)   !temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !temporary copy
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info

    INTEGER :: iL,iF,iI

    DO iL = 1,iNumLayer
        iI = iaaRadLayer(iAtm,iL)
        DO iF = 1,kMaxPts
            raaExtTemp(iF,iI)  = raaAbs(iF,iI)
            raaScatTemp(iF,iI) = 0.0
            raaAsymTemp(iF,iI) = 0.0
        END DO
    END DO

    RETURN
    end SUBROUTINE CopyRaaExt_twostream

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! it also does the scattering part of the cloud stuff, so that we can
! use the formulation of Pat Arnott of the DRI
    SUBROUTINE AddCloud_twostream(raFreq,raaExtTemp,raaScatTemp, &
    raaAsymTemp,iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaScatTemp(kMaxPts,kMixFilRows)   !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND

! mie scattering tables
    INTEGER :: NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::    IWP(MAXNZ), DME(MAXNZ)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)

! local variables
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do

    IF (rScat > 0.0) THEN
    !!!!first find out where the cloud top is, in kCARTA layering
        N = iCldTop
        iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        N = iCldBot-1
        iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    ! ow get the optical properties for the cloud layers
    ! his is the original code
    !     $      kcarta layers           iaCldTop(iIn)-1,' to ',iaCldBot(iIn)
    !     $      rtspec layers           iaCldTop(iIn)+1,' to ',iaCldBot(iIn)
    ! cloud is in KCARTA layers           45 to           42
    ! cloud is in RTSPEC layers           56 to           59

        DO N = ICLDTOP, ICLDBOT-1
            L  = N-ICLDTOP+1
            I  = ISCATTAB(L)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaExtTemp(iF,iI)
                rAbs   = taugas
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
            !  Compute the optical depth of cloud layer, including gas
                TAUC_L   = IWP(L)*EXTINCT/1000
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L
                raaExtTemp(iF,iI)  = TAUTOT_N
            ! also save the SCAT coeff
                raaScatTemp(iF,iI) = IWP(L)*EXTINCT/1000*SSALB(L)
                SSALB(L) = SSALB(L)*TAUC_L/TAUCG_L
                raaScatTemp(iF,iI) = SSALB(L)*TAUTOT_N
            ! also save the ASYM coeff
                IF (IWP(L) >= 1.0e-5) THEN
                    raaAsymTemp(iF,iI) = ASYM_RTSPEC(L)
                ELSE
                    raaAsymTemp(iF,iI) = 0.0
                END IF
            !        print *,waveno,DME(L),IWP(L),TAUC_L,SSALB(L),ASYM_RTSPEC(L),
            !     $          raaExtTemp(iF,iI),raaScatTemp(iF,iI),raaAsymTemp(iF,iI)
            !            call dostopMesg('AddCloud_twostream $')
            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

! now use the partial fractions
    i1  = iaaRadLayer(iAtm,1)
    i2  = iaaRadLayer(iAtm,iNumLayer)
    iFixHere = -1         !!!do not adjust here, scatter_twostream does it
    iFixHere = +1         !!!do adjust here, scatter_twostream does not
    iFixHere = -1
    IF (iFixHere > 0) THEN
        IF (i1 > i2) THEN
        ! adiation going from eg layer 100 to 1 ==> up look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)  = raaExtTemp(iF,i1) * rFracTop
                raaScatTemp(iF,i1) = raaScatTemp(iF,i1) * rFracTop
                raaExtTemp(iF,i2)  = raaExtTemp(iF,i2) * rFracBot
                raaScatTemp(iF,i2) = raaScatTemp(iF,i2) * rFracBot
            END DO
        ELSEIF (i1 < i2) THEN
        ! adiation going from eg layer 1 to 100 ==> down look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)  = raaExtTemp(iF,i1) * rFracBot
                raaScatTemp(iF,i1) = raaScatTemp(iF,i1) * rFracBot
                raaExtTemp(iF,i2)  = raaExtTemp(iF,i2) * rFracTop
                raaScatTemp(iF,i2) = raaScatTemp(iF,i2) * rFracTop
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE AddCloud_twostream
!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!   *** it also adds on the "backscattered" part for PCLSAM algorithm ***
! this way we have a fast alternative to kTwoStream
    SUBROUTINE AddCloud_pclsam(raFreq, &
    raaExtTemp,raaScatTemp,raaAsymTemp, &
    iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaScatTemp(kMaxPts,kMixFilRows)   !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND

! mie scattering tables
    INTEGER :: NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::    IWP(MAXNZ), DME(MAXNZ)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)

! local variables
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do

    IF (rScat > 0.0) THEN
    !!!!first find out where the cloud top is, in kCARTA layering
        N = iCldTop
        iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        N = iCldBot-1
        iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    ! ow get the optical properties for the cloud layers
        DO iI = 1,kMixFilRows
            DO iF = 1,kMaxPts
                raaAsymTemp(iF,iI) = 0.0
            END DO
        END DO

        DO N = ICLDTOP, ICLDBOT-1
            L  = N-ICLDTOP+1
            I  = ISCATTAB(L)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaExtTemp(iF,iI)
                rAbs   = taugas
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
            !  Compute the optical depth of cloud layer, including gas
                TAUC_L   = IWP(L)*EXTINCT/1000
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L

            ! the SSALB coeff
                rScat    = SSALB(L) * IWP(L)*EXTINCT/1000
                SSALB(L) = SSALB(L)*TAUC_L/TAUCG_L
                raaScatTemp(iF,iI) = SSALB(L)

            ! ---------------> now add on the backscattered part <--------------------
                b = (1.0 - ASYM_RTSPEC(L))/2.0
                TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b))
                raaExtTemp(iF,iI)  = TAUTOT_N
            ! ---------------> now add on the backscattered part <--------------------

                IF (IWP(L) >= 1.0e-5) THEN
                    raaAsymTemp(iF,iI) = ASYM_RTSPEC(L)
                ELSE
                    raaAsymTemp(iF,iI) = 0.0
                END IF

            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

! now use the partial fractions
    i1  = iaaRadLayer(iAtm,1)
    i2  = iaaRadLayer(iAtm,iNumLayer)
    iFixHere = -1         !!!do not adjust here, scatter_twostream does it
    iFixHere = +1         !!!do adjust here, scatter_twostream does not
    iFixHere = -1
    IF (iFixHere > 0) THEN
        IF (i1 > i2) THEN
        ! adiation going from eg layer 100 to 1 ==> up look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)   = raaExtTemp(iF,i1) * rFracTop
                raaExtTemp(iF,i2)   = raaExtTemp(iF,i2) * rFracBot
            ! do not need this since this is a ratio
            !            raaSSAlbTemp(iF,i1) = raaSSAlbTemp(iF,i1) * rFracTop
            !            raaSSAlbTemp(iF,i2) = raaSSAlbTemp(iF,i2) * rFracBot
            END DO
        ELSEIF (i1 < i2) THEN
        ! adiation going from eg layer 1 to 100 ==> down look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)   = raaExtTemp(iF,i1) * rFracBot
                raaExtTemp(iF,i2)   = raaExtTemp(iF,i2) * rFracTop
            ! do not need this since this is a ratio
            !            raaSSAlbTemp(iF,i1) = raaSSAlbTemp(iF,i1) * rFracBot
            !            raaSSAlbTemp(iF,i2) = raaSSAlbTemp(iF,i2) * rFracTop
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE AddCloud_pclsam

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!  1) it also adds on the "backscattered" part for PCLSAM algorithm
! this way we have a fast alternative to kTwoStream
!  2) does the jacobian part for d/d(DME)
! this is for a DOWNLOOK instrument, so we call
!        raaPhaseJacobASYM(iF,iI) = hg2_real_deriv_wrt_g(-mu_sun,mu_sat,ASYM)
    SUBROUTINE AddCloud_pclsam_Jacob_downlook(raFreq,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raaPhaseJacobASYM, &
    iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! mie scattering tables
    INTEGER ::  NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::     IWP(MAXNZ), DME(MAXNZ)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

! local variables
    REAL :: mu_sun,mu_sat
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat
    REAL :: OMEGA, ASYM,tautotal_0

    REAL :: dEXTINCT_dr, dSSALB_dr, dASYM_dr
    REAL :: rW,x1,x2,x3,x4,x5
    REAL :: hg2_real

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do

    IF (rScat > 0.0) THEN
    !!!!first find out where the cloud top is, in kCARTA layering
        N = iCldTop
        iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        N = iCldBot-1
        iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    ! ow get the optical properties for the cloud layers
        DO iI = 1,kProfLayerJac
            DO iF = 1,kMaxPts
                raaPhaseJacobASYM(iF,iI) = 0.0
                raaExtJacobIWP(iF,iI)    = 0.0
                raaSSAlbJacobIWP(iF,iI)  = 0.0
                raaAsymJacobIWP(iF,iI)   = 0.0
                raaExtJacobDME(iF,iI)    = 0.0
                raaSSAlbJacobDME(iF,iI)  = 0.0
                raaAsymJacobDME(iF,iI)   = 0.0
            END DO
        END DO

        DO N = ICLDTOP, ICLDBOT-1
            L  = N-ICLDTOP+1
            I  = ISCATTAB(L)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)
            mu_sat = cos(raLayAngles(iI)*kPi/180)
            mu_sun = cos(raSunAngles(iI)*kPi/180)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaExtTemp(iF,iI)
                rAbs   = taugas
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                CALL JACOBIAN_INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                dEXTINCT_dr, dSSALB_dr, dASYM_dr, &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                OMEGA = SSALB(L)
                ASYM  = ASYM_RTSPEC(L)

            !  Compute the optical depth of cloud layer, including gas
                TAUC_L   = IWP(L)*EXTINCT/1000
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L

            !   the SSALB coeff
                rW       = SSALB(L)
                rScat    = SSALB(L) * IWP(L)*EXTINCT/1000
                SSALB(L) = SSALB(L) * TAUC_L/TAUCG_L
                raaSSAlbTemp(iF,iI) = SSALB(L)

            ! ---------------> now add on the backscattered part <--------------------
                b = (1.0 - ASYM_RTSPEC(L))/2.0
                TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b))
                raaExtTemp(iF,iI)  = TAUTOT_N
            ! ---------------> now add on the backscattered part <--------------------

                IF (IWP(L) >= 1.0e-5) THEN
                    raaAsymTemp(iF,iI) = ASYM_RTSPEC(L)
                ELSE
                    raaAsymTemp(iF,iI) = 0.0
                END IF
            ! -------------------------- now do the jacobians --------------------------
            !! technically we are doing d/d(DME) and not d/d(RME); they are
            ! related by raaXYZJacobRME(iF,iI) = raaXYZJacobDME(iF,iI)

                tautotal_0 = TAUCG_L

            !! --------> d/d(iwp) <---------  !!
                x1 = EXTINCT/1000
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaExtJacobIWP(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2
                 
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaSSAlbJacobIWP(iF,iI) = x2

                raaAsymJacobIWP(iF,iI) = 0.0

            !! --------> d/d(dme) <---------  !!
                x1 = IWP(L)/1000*dEXTINCT_dr
                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                x3 = -1/2*dASYM_dr
                raaExtJacobDME(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2 + &
                TAUCG_L*SSALB(L)*x3

                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                raaSSAlbJacobDME(iF,iI) = x2

                raaAsymJacobDME(iF,iI) = dASYM_dr

            !! --------> d/d(g) <---------  !!
                raaPhaseJacobASYM(iF,iI) = &
                hg2_real_deriv_wrt_g(-mu_sun,mu_sat,ASYM)

            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

! now use the partial fractions????? see last section in
!       SUBROUTINE AddCloud_pclsam( )

    RETURN
    end SUBROUTINE AddCloud_pclsam_Jacob_downlook

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!  1) it also adds on the "backscattered" part for PCLSAM algorithm
! this way we have a fast alternative to kTwoStream
!  2) does the jacobian part for d/d(DME)
! this is for a UPLOOK instrument, so we call
!       raaPhaseJacobASYM(iF,iI) = hg2_real_deriv_wrt_g(-mu_sun,-mu_sat,ASYM)
    SUBROUTINE AddCloud_pclsam_Jacob_uplook(raFreq,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raaPhaseJacobASYM, &
    iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! mie scattering tables
    INTEGER ::  NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::     IWP(MAXNZ), DME(MAXNZ)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

! local variables
    REAL :: mu_sun,mu_sat
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat
    REAL :: OMEGA, ASYM,tautotal_0

    REAL :: dEXTINCT_dr, dSSALB_dr, dASYM_dr
    REAL :: rW,x1,x2,x3,x4,x5
    REAL :: hg2_real

    rScat = 0.0
    do i1 = 1,maxnz
        rScat = rScat + iwp(i1)
    end do

    IF (rScat > 0.0) THEN
    !!!!first find out where the cloud top is, in kCARTA layering
        N = iCldTop
        iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        N = iCldBot-1
        iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    ! ow get the optical properties for the cloud layers
        DO iI = 1,kMixFilRows
            DO iF = 1,kMaxPts
                raaPhaseJacobASYM(iF,iI) = 0.0
                raaExtJacobIWP(iF,iI)    = 0.0
                raaSSAlbJacobIWP(iF,iI)  = 0.0
                raaAsymJacobIWP(iF,iI)   = 0.0
                raaExtJacobDME(iF,iI)    = 0.0
                raaSSAlbJacobDME(iF,iI)  = 0.0
                raaAsymJacobDME(iF,iI)   = 0.0
            END DO
        END DO

        DO N = ICLDTOP, ICLDBOT-1
            L  = N-ICLDTOP+1
            I  = ISCATTAB(L)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)
            mu_sat = cos(raLayAngles(iI)*kPi/180)
            mu_sun = cos(raSunAngles(iI)*kPi/180)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaExtTemp(iF,iI)
                rAbs   = taugas
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                CALL JACOBIAN_INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                dEXTINCT_dr, dSSALB_dr, dASYM_dr, &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                OMEGA = SSALB(L)
                ASYM  = ASYM_RTSPEC(L)

            !  Compute the optical depth of cloud layer, including gas
                TAUC_L   = IWP(L)*EXTINCT/1000
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L

            !   the SSALB coeff
                rW       = SSALB(L)
                rScat    = SSALB(L) * IWP(L)*EXTINCT/1000
                SSALB(L) = SSALB(L) * TAUC_L/TAUCG_L
                raaSSAlbTemp(iF,iI) = SSALB(L)

            ! ---------------> now add on the backscattered part <--------------------
                b = (1.0 - ASYM_RTSPEC(L))/2.0
                TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b))
                raaExtTemp(iF,iI)  = TAUTOT_N
            ! ---------------> now add on the backscattered part <--------------------

                IF (IWP(L) >= 1.0e-5) THEN
                    raaAsymTemp(iF,iI) = ASYM_RTSPEC(L)
                ELSE
                    raaAsymTemp(iF,iI) = 0.0
                END IF
            ! -------------------------- now do the jacobians --------------------------
            !! technically we are doing d/d(DME) and not d/d(RME); they are
            ! related by raaXYZJacobRME(iF,iI) = raaXYZJacobDME(iF,iI)

                tautotal_0 = TAUCG_L

            !! --------> d/d(iwp) <---------  !!
                x1 = EXTINCT/1000
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaExtJacobIWP(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2
                 
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaSSAlbJacobIWP(iF,iI) = x2

                raaAsymJacobIWP(iF,iI) = 0.0

            !! --------> d/d(dme) <---------  !!
                x1 = IWP(L)/1000*dEXTINCT_dr
                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                x3 = -1/2*dASYM_dr
                raaExtJacobDME(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2 + &
                TAUCG_L*SSALB(L)*x3

                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                raaSSAlbJacobDME(iF,iI) = x2

                raaAsymJacobDME(iF,iI) = dASYM_dr

            !! --------> d/d(g) <---------  !!
                raaPhaseJacobASYM(iF,iI) = &
                hg2_real_deriv_wrt_g(-mu_sun,-mu_sat,ASYM)

            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

! now use the partial fractions????? see last section in
!       SUBROUTINE AddCloud_pclsam( )

    RETURN
    end SUBROUTINE AddCloud_pclsam_Jacob_uplook

!************************************************************************
! this subroutine copies the relevant parts of raaAbs to raaExtTemp
    SUBROUTINE CopyRaaAbs(raaAbs,raaExtTemp,iaaRadLayer,iAtm,iNumlayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
    REAL :: raaAbs(kMaxPts,kMixFilRows)        !original, from uncompression
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !temporary copy
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info

    INTEGER :: iL,iF,iI

    DO iL = 1,iNumLayer
        iI = iaaRadLayer(iAtm,iL)
        DO iF = 1,kMaxPts
            raaExtTemp(iF,iI) = raaAbs(iF,iI)
        END DO
    END DO

    RETURN
    end SUBROUTINE CopyRaaAbs

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
    SUBROUTINE AddCloud_absorbonly( &
    raFreq,raaExtTemp,iaaRadLayer,iAtm,iNumlayer, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !temporary copy
    REAL :: raFreq(kMaxPts)                   !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms

! mie scattering tables
    INTEGER :: NCLDLAY, ICLDTOP, ICLDBOT, ISCATTAB(MAXNZ)
    REAL ::    IWP(MAXNZ), DME(MAXNZ)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)

! local variables
    INTEGER :: iL,iF,iI,N,L,I
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ),rScat

    rScat = 0.0
    do iI = 1,maxnz
        rScat = rScat + iwp(iI)
    end do

    IF (rScat > 0) THEN
    !       !Get the optical properties for the cloud layers
        DO N = ICLDTOP, ICLDBOT - 1
            L  = N-ICLDTOP+1
            I      = ISCATTAB(L)
            iL = kProfLayer - N + 1
        !          iI     = iaaRadLayer(iAtm,iL)
            iI = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,iL)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaExtTemp(iF,iI)
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
            !  Compute the optical depth of cloud layer, including gas
            !  differs from typical DISORT/RTSPEC since we only want ABSORPTION
            !  and NOT EXTINCT = ABS + SCATTER
                TAUC_L   = IWP(L)*EXTINCT/1000*(1.0-SSALB(L))
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L
                raaExtTemp(iF,iI) = TAUTOT_N

            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

    RETURN
    end SUBROUTINE AddCloud_absorbonly
!************************************************************************
! this subroutine reads in the phase info associated with the file
    SUBROUTINE READ_PHASE(SCATFILE,raFreq,rDmePhase,ndme,DMETAB_PHASE, &
    raPhasePoints,raComputedPhase)

    IMPLICIT NONE
          
    include '../INCLUDE/scatterparam.f90'

! input : the scattering file name
    CHARACTER(120) :: SCATFILE
    REAL :: raFreq(kMaxPts),rDmePhase      !!!wavenumber and particle size
    REAL :: DMETAB_PHASE(kProfLayer)        !!!partice size info
    INTEGER :: ndme                         !!!number of particles
! output : the phase info associated with the cloud; else use HG
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    CHARACTER(120) :: caPhaseFile
    INTEGER :: iI,iErr,iN,iJ,iS,iJump
    REAL :: rW,rD,rX,rY,rMid,slope

    REAL :: raComputedPhase1(MaxPhase),raComputedPhase2(MaxPhase),rD1,rD2

    rMid = raFreq(kMaxPts/2)

    iI = 1
    10 CONTINUE
    IF (SCATFILE(iI:iI) /= ' ') THEN
        iI = iI + 1
        GOTO 10
    END IF
    caPhaseFile = scatfile
    caPhaseFile(iI:iI+6) = '.phase'

    OPEN (UNIT = kTempUnit, STATUS='OLD', FORM='FORMATTED', &
    FILE=caPhaseFile, IOSTAT = iERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,1010) IERR, SCATFILE
        CALL DoSTOP
    ENDIF
    1010 FORMAT('ERROR! number ',I5,' opening phase scatter data file:',/,A120)
    kTempUnitOpen=1

    read(kTempUnit,*) rW,rD,iN,rX,rY    !!!read the first line

    20 CONTINUE
    IF (rW < rMid) THEN
        DO iJ = 2,iN
        !!!! jump over the current (rW,rD) set of iN  points
            read(kTempUnit,*) rW,rD,iN,rX,rY !!!skip over these lines
        END DO
        read(kTempUnit,*) rW,rD,iN,rX,rY  !!!read the "first" line of (rW,rD')
        GOTO 20
    END IF

!!! now we are at the wavenumber point we wish to use! proceed!

    IF (rDmePhase <= DMETAB_PHASE(1)) THEN
    !!! we are at the data we need (FIRST DIAMETER IS ONE NEEDED!!!)
        raPhasePoints(1)   = rX
        raComputedPhase(1) = rY
        DO iJ = 2,iN
            read(kTempUnit,*) rW,rD,iN,rX,rY
            raPhasePoints(iJ)   = rX
            raComputedPhase(iJ) = rY
        END DO
    ELSEIF (rDmePhase >= DMETAB_PHASE(ndme)) THEN
    !!! skip to the last data set for this wavenumber
    !!! (LAST DIAMETER IS ONE NEEDED!!!)
        DO iS = 1,ndme-1
            DO iJ = 2,iN
                read(kTempUnit,*) rW,rD,iN,rX,rY
            END DO
            read(kTempUnit,*) rW,rD,iN,rX,rY
        END DO
    !!!we are at the data we need
        raPhasePoints(1)   = rX
        raComputedPhase(1) = rY
        DO iJ = 2,iN
            read(kTempUnit,*) rW,rD,iN,rX,rY
            raPhasePoints(iJ)   = rX
            raComputedPhase(iJ) = rY
        END DO

    ELSE
    !!! pah! need to interpolate between particle sizes
        iJump = 1
        30 CONTINUE
        IF (DMETAB_PHASE(iJump) < rDmePhase) THEN
            iJump = iJump + 1
            GOTO 30
        END IF
        iJump = iJump-1
        DO iS = 1,iJump-1
            DO iJ = 2,iN
                read(kTempUnit,*) rW,rD,iN,rX,rY
            END DO
            read(kTempUnit,*) rW,rD,iN,rX,rY
        END DO

    !!!we are before the data we need
        rD1 = rD
        raPhasePoints(1)   = rX
        raComputedPhase1(1) = rY
        DO iJ = 2,iN
            read(kTempUnit,*) rW,rD,iN,rX,rY
            raPhasePoints(iJ)   = rX
            raComputedPhase1(iJ) = rY
        END DO
    !!!we are after the data we need
        read(kTempUnit,*) rW,rD,iN,rX,rY
        rD2 = rD
        raPhasePoints(1)   = rX
        raComputedPhase2(1) = rY
        DO iJ = 2,iN
            read(kTempUnit,*) rW,rD,iN,rX,rY
            raPhasePoints(iJ)   = rX
            raComputedPhase2(iJ) = rY
        END DO
    !!!!now do average
        DO iJ = 1,iN
            slope = (raComputedPhase2(iJ)-raComputedPhase1(iJ))/log(rD2/rD1)
            raComputedPhase(iJ) = raComputedPhase2(iJ)-slope*log(rD2/rDmephase)
        END DO
    END IF

    CLOSE(kTempUnit)
    kTempUnitOpen=-1

!!!now we need to flip things if necessary
    IF (raPhasePoints(1) > raPhasePoints(MaxPhase)) THEN
        DO iJ = 1,MaxPhase/2
            slope = raPhasePoints(iJ)
            raPhasePoints(iJ) = raPhasePoints(MaxPhase-iJ+1)
            raPhasePoints(MaxPhase-iJ+1) = slope

            slope = raComputedPhase(iJ)
            raComputedPhase(iJ) = raComputedPhase(MaxPhase-iJ+1)
            raComputedPhase(MaxPhase-iJ+1) = slope
        END DO
    END IF

!      DO iJ = 1,10
!        print *,iJ,rW,rD1,rDmephase,rD2,raPhasePoints(iJ),
!     $   raComputedPhase1(iJ),raComputedPhase2(iJ),raComputedPhase(iJ)
!      END DO

    RETURN
    END SUBROUTINE READ_PHASE

!************************************************************************
! this function computes the Henyey Greenstein function, assuming the
! cos(phi1-phi2) factor = 1
    DOUBLE PRECISION Function hg2_double(mu1,mu2,g)

    IMPLICIT NONE

    DOUBLE PRECISION :: mu1,mu2,g  !mu1,mu2 are the two angles
!g is the asymmetry
    DOUBLE PRECISION :: normB,mu0,yexact

! normB is normalisation of mu from -1 to 1 and works out to be 2.0
!! we also know that (1/2) integral P(-1,1) = 1
! ormB = 1/sqrt(1+g*g - 2.0*g) - 1/sqrt(1+g*g + 2.0*g)
! ormB = (1-g*g)/g * normB
! ormB = 2.0

!!!compute mu0 = cos ofangle between the two
    mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)
     
    yexact = (1 + g*g - 2.0*g*mu0) * sqrt(1 + g*g - 2.0*g*mu0)
    yexact = (1-g*g)/yexact
! exact = yexact/normB * 2.0
    yexact = yexact

    hg2_double = yexact

    RETURN
    END FUNCTION hg2_double

!************************************************************************
! this function computes the Henyey Greenstein function, assuming the
! cos(phi1-phi2) factor = 1
    REAL Function hg2_real(mu1,mu2,g)

    IMPLICIT NONE

    REAL :: mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry

    REAL :: normB,mu0,yexact

! normB is normalisation of mu from -1 to 1 and works out to be 2.0
!! we also know that (1/2) integral P(-1,1) = 1
! ormB = 1/sqrt(1+g*g - 2.0*g) - 1/sqrt(1+g*g + 2.0*g)
! ormB = (1-g*g)/g * normB
! ormB = 2.0

!!!compute mu0 = cos ofangle between the two
    mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)
     
    yexact = (1 + g*g - 2.0*g*mu0) * sqrt(1 + g*g - 2.0*g*mu0)
    yexact = (1-g*g)/yexact
! exact = yexact/normB * 2.0
    yexact = yexact

    hg2_real = yexact

    RETURN
    end Function hg2_real

!************************************************************************
! this function computes the d/dg for Henyey Greenstein function, assuming the
! cos(phi1-phi2) factor = 1
! see /home/sergio/MATLABCODE/RADTrans/GENERAL_CLOUD/hg2_deriv.m
    REAL Function hg2_real_deriv_wrt_g(mu1,mu2,g)

    IMPLICIT NONE

    REAL :: mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry

    REAL :: nn,normB,dnormB,mu0,y,y1,yder
     
! normB is normalisation of mu from -1 to 1 and works out to be 2.0
! nn    = 1/sqrt(1+g*g - 2.0*g) - 1/sqrt(1+g*g + 2.0*g)
! normB = (1-g*g)/g * nn
! normB = 2.0
! we also know that (1/2) integral P(-1,1) = 1
     
! normB = ((g-2)/((1+g*g -2.0*g)**(3/2)) - (g+2)/((1+g*g +2.0*g)**(3/2)))
! normB = -(1+g*g)/(g*g)*nn - (1-g*g)/g*dnormB
! normB = -2 * dnormB/(normB*normB)
    dnormB = 0
     
    mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2)
     
    y1 = (1-g*g)/((1 + g*g - 2.0*g*mu0)**(3.0/2.0))
!y = y1/normB *2
    y = y1
     
! ow compute the derivative
    yder = (1 + g*g - 2.0*g*mu0)
    yder = (yder*(-2.0*g) - 3*(1-g*g)*(g-mu0))/(yder**(5.0/2.0))
! der = yder/normB * 2
    yder = yder

    hg2_real_deriv_wrt_g = yder

    RETURN
    end Function hg2_real_deriv_wrt_g

!************************************************************************
!     this is to take in arbitrary cloud profiles and add them together
!     using the Sun Shine prescription

! this is from Dave Turner's thesis : refer pg 64 of the thesis
! Studies of the radiative properties of ice and mixed-phase clouds
! Authors: ZHIAN SUN; KEITH P. SHINE
! Source: Quarterly Journal of the Royal Meteorological Society,
! Volume 120, Number 515, January 1994 Part A, pp. 111-137(27)
    SUBROUTINE SetMieTables_RTSPEC_100layer(raFreq, &
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
    iSergio, &
!!!!!!!!!!!!!!!!!! these are the cloud profiles PLUS output
    iaCldTypes,raaKlayersCldAmt,raVTemp, &
!!!!!!!!!!!!!!!!!! these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
    raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)

    implicit none
     
    include '../INCLUDE/scatterparam.f90'

! iSergio INTEGER that tells if this is RTSPEC or SERGIO's code
    INTEGER :: iSergio
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
! ---------------- inputs needed to read scattering tables -------------------
! this is which atm number is being used, and whether these are binary files
    INTEGER :: iAtm,iBinaryFile,iNumLayer,iDownward
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this is just to set everything about clouds relative to TOA layer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
!! for some reason I had used kMaxWater instead of kMaxClouds here??????
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
    INTEGER :: iaCldTypes(kMaxClouds)

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
          
! c      INTEGER  MAXTAB, MAXGRID, MAXSCAT
! c      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
    CHARACTER(120) :: SCATFILE(MAXSCAT)
     
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
     
    INTEGER :: NSCATTAB, NCLDLAY, NLEV, NABSNU
    INTEGER :: ICLDTOP, ICLDBOT, IOBS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
     
    INTEGER :: iCloudySky        !!!!are there clouds in this atm???

! ---------------------------- local variables ----------------------------
    INTEGER :: IACLDTOPKCARTA(kMaxClouds), IACLDBOTKCARTA(kMaxClouds)
    INTEGER :: iaTable(kMaxClouds*kCloudLayers),iIn,iJ,iReadTable,I
    INTEGER :: iCloud,iStep
    REAL :: extinct
    INTEGER :: LL,II,N,M,iLayers,iBlah

    INTEGER :: iL,iT,iB,iNumClds,iT_Atm,iB_Atm,iaCldInLayer(kProfLayer)

    REAL ::    raCldLayer(MAXNZ),iwp0(maxnz),dme0(maxnz), rDmePhase
    REAL ::    dmetab_phase(kProfLayer)
    INTEGER :: indx(MAXNZ),iscattab0(maxnz),iiDiv
     
    CHARACTER(120) :: caName
    CHARACTER(1) :: caScale(MAXSCAT)

    REAL :: rX,rY,raTempLay(kProfLayer),rT
    REAL :: raC(0:3)   !!! KN Liou ice dme parameters

    INTEGER :: iF,iG,iaRadLayer(kProfLayer),iC

    raC(0) = 326.3
    raC(1) = 12.42
    raC(2) = 0.197
    raC(3) = 0.0012

    DO iL = 1,iNumLayer
        iaRadLayer(iL) = iaaRadLayer(iAtm,iL)
        raTempLay(iL)  = raVtemp(iaRadLayer(iL))
    END DO
    DO iL = iNumLayer+1,kProfLayer
        iaRadLayer(iL) = -1
    END DO

! ------------- >> find the scattering params <<---------------
! nitialise all scattering info to null

    iiDiv = 0
    555 CONTINUE
    IF (iiDiv*kProfLayer < iaaRadLayer(iAtm,3)) THEN
        iiDiv = iiDiv + 1
    END IF
    iiDiv = iiDiv - 1

! copied from s_scatter_spectra.f .. all table names etc are unique, so no
! need to make more checks
! Life is easy since cloud profile occupies ALL layers!

    iCloudySky = -1        !!!!!!!assume no clouds in sky

!!!!!!!!this is all that is needed if only RTSPEC rad transfer were used
    IF (iDownWard == 1) THEN
        iB_Atm = iaaRadLayer(iAtm,1)
        iT_Atm = iaaRadLayer(iAtm,iNumLayer)
    ELSEIF (iDownWard == -1) THEN
        iT_Atm = iaaRadLayer(iAtm,1)
        iB_Atm = iaaRadLayer(iAtm,iNumLayer)
    END IF
!!!!!!however we also do fluxes, so even if the atm is defined so it
!!!!!!is for an uplook instrument, RTSPEC will be called in a downlook
!!!!!!fashion, and vice versa

    IF (iB_Atm > iT_Atm) THEN
        iCloudySky = iT_Atm
        iT_Atm = iB_Atm
        iB_atm = iCloudySky
    END IF
! nitialise all scattering info to null

    iCloudySky = -1        !!!!!!!assume no clouds associated with this atm

    iCldTopKcarta = -1
    iCldBotKcarta = kProfLayer+1
    NSCATTAB = -1000
    DO iIn = 1,kMaxClouds*kCloudLayers
        iaTable(iIn) = -1
    END DO
    DO iIn = 1,MAXSCAT
        ScatFile(iIn) = ' '
        iaScatTable_With_Atm(iIn) = -1
    END DO
    DO iIn = 1,kMaxClouds
        iaCloudWithThisAtm(iIn) = -1
        iaCldTop(iIn) = -1
        iaCldBot(iIn) = -1
        iaCldTopkCarta(iIn) = -1
        iaCldBotkCarta(iIn) = -1
    END DO

    DO iIn = 1,iNclouds
        DO iJ = 1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            IF (iI > MAXSCAT) THEN
                write(kStdErr,*)'unfortunately, in scatterparam.f90 we have '
                write(kStdErr,*)'MAXSCAT = ',maxscat
                write(kStdErr,*)'please reset and retry'
                CALL DoSTOP
            END IF
            caName=caaaScatTable(iIn,iJ)
            IF (iaTable(iI) < 0) THEN  !nothing associated with this yet
                IF (iI > NSCATTAB) THEN
                    NSCATTAB = iI
                END IF
                iaTable(iI) = 1
                ScatFile(iI) = caName
            END IF
        END DO

    ! check to see if this cloud is to be used with this atm . recall that
    ! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
    ! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    ! INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
        DO iJ = 1,iaCloudNumAtm(iIn)
            IF (iaaCloudWhichAtm(iIn,iJ)  == iAtm) THEN
                iCloudySky = iIn         !!!! set this up
                iaCloudWithThisAtm(iIn) = 1
                IACLDTOP(iIn) = iaaCloudWhichLayers(iIn,1)+1
                IACLDBOT(iIn) = iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn))
                iaCldTopkCarta(iIn) = iaCldTop(iIn)     !!! not needed
                iaCldBotkCarta(iIn) = iaCldBot(iIn)     !!! not needed
            ! iCldTopkCarta = iaCldTop(iIn)-1
            ! iCldBotkCarta = iaCldBot(iIn)
                IF (iCldTopkCarta < iaCldTop(iIn)-1) THEN
                    iCldTopkCarta = iaCldTop(iIn)-1
                END IF
                IF (iCldBotkCarta > iaCldBot(iIn)) THEN
                    iCldBotkCarta = iaCldBot(iIn)
                END IF
                write(kStdWarn,*)'cloud # ',iIn,' associated with atm # ',iAtm
                write(kStdWarn,*)'setmie1 : cloud is in KCARTA layers ', &
                iiDiv*kProfLayer+iaCldTop(iIn)-1,' to ', &
                iiDiv*kProfLayer+iaCldBot(iIn)
            !!!!!these are the RTSPEC layers 100 to 1 = GND to TOA
                iaCldbot(iIn) = iT_Atm - iaCldbot(iIn) + 1
                iaCldtop(iIn) = iT_Atm - iaCldtop(iIn) + 1
                write(kStdWarn,*)'setmie1 : cloud is in RTSPEC layers ', &
                iaCldTop(iIn)+1,' to ',iaCldBot(iIn)

            END IF
        END DO
         
    ! check to see which scattering tables to be used with this atm
        DO iJ = 1,1
            iI = iaaScatTable(iIn,iJ)
            IF (iaCloudWithThisAtm(iIn) == 1) THEN
                iaScatTable_With_Atm(iI) = 1
                write(kStdWarn,*)'scat table ',iI,' for atm,layer # ',iAtm,iJ
            END IF
        END DO
        DO iJ = 1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            IF (iaCloudWithThisAtm(iIn) == 1) THEN
                iaScatTable_With_Atm(iI) = 1
            END IF
        END DO
    END DO      !!!!!!!!main       DO iIn=1,iNclouds

!     Only read in scattering tables that are needed for this atm
    iReadTable = 1
    IF (iReadTable > 0) THEN
        IF (iBinaryFile == 1) THEN
            DO I = 1, NSCATTAB
                write(kStdWarn,*) ' '
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading binary scatter data for table #',I
                    write(kStdWarn,*) scatfile(I)
                !              print *,'lalala',ka100layerCloudType(I)
                !              IF (ka100layerCloudType(I) .EQ. 201) PRINT *,'usususus 201 ice'
                    CALL READ_SSCATTAB_BINARY(SCATFILE(I),   & !!!!!!MAXTAB, MAXGRID,
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))
                    IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR. &
                    (ABS(MUINC(2)-0.7887) > 0.001)) THEN
                        write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                        CALL DoStop
                    END IF
                    IF (iaPhase(I) > 0) THEN
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        ELSE IF (iBinaryFile == -1) THEN
            DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading ascii scatter data for table #',I
                    CALL READ_SSCATTAB(SCATFILE(I),   & !!!!!!MAXTAB, MAXGRID,
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))

                    IF ((ABS(MUINC(1)-0.2113) > 0.001) .OR. &
                    (ABS(MUINC(2)-0.7887) > 0.001)) THEN
                        write(kStdErr,*) 'RTSPEC: Coded for incident mu=0.2113,0.7887'
                        CALL DoStop
                    END IF
                    IF (iaPhase(I) > 0) THEN
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        ELSE IF (iBinaryFile == 0) THEN
            DO I = 1, NSCATTAB
                IF (iaScatTable_With_Atm(I) > 0) THEN
                    write(kStdWarn,*) 'Reading ascii scatter data for table #',I
                    CALL READ_SSCATTAB_SPECIAL(SCATFILE(I), &
                    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
                    NWAVETAB(I), WAVETAB(1,I), &
                    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                    TABPHI1UP(1,I), TABPHI1DN(1,I), &
                    TABPHI2UP(1,I), TABPHI2DN(1,I))

                    IF (iaPhase(I) > 0) THEN
                        write(kStdErr,*) 'Right now, incapapable of this silly task!!!'
                        write(kStdErr,*) 'need iaPhase = 0 for iBinaryFIle = 0'
                        CALL DoStop
                        DO iBlah = 1,NDME(I)
                            dmetab_phase(iBlah) = DMETAB(iBlah,I)
                        END DO
                        rDmePhase = raaaCloudParams(I,1,2)
                        CALL READ_PHASE(SCATFILE(I),raFreq,rDmePhase,ndme(I),dmetab, &
                        raPhasePoints,raComputedPhase)
                    END IF
                END IF
            ENDDO
        END IF    !iBinaryFile > 0
    END IF      !iReadTable  > 0

! Frank Evans code scales the Mie scattering parameters, so if we are using
! my canned EDDINGTON method, we have to unscale them!!!!!!!!
    IF (iSergio > 0) THEN
        DO I = 1, NSCATTAB
            IF (iaScatTable_With_Atm(I) > 0) THEN
                CALL UnScaleMie( &
                caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
                ndme(i)*nwavetab(i))
            END IF
        END DO
    END IF

    DO i=1,MAXNZ
        raCldLayer(I) = +1.0e10
        indx(I)       = -1
        iscattab0(I)  = -1
        dme0(I)       = -1.0
        iwp0(I)       = -1.0
    END DO

    iCloud = -1
    IF (iCloudySky < 0) THEN
    !!!!!this is similar to DISORT interface
        write(kStdWarn,*)'Could not find a cloud for atmosphere #',iAtm
        write(kStdWarn,*)'setting IWP = -100.0'
        iCloud=1    !say cloud number one associated with this atmosphere
        ncldlay=1   !say fictitious cloud occupies one layer
        raaIWP(1,1)     = -100.0   !but ensure cloud has NO particles in it!
        raaDME(1,1)     = -10.0    !but ensure cloud has NO particles in it!
        iaaSCATTAB(1,1) = -1
         
    ELSE
    !!!!!find total number of clouds, and hence total number of layers
    !!!!!this is quite different to DISORT interface, as we have to make
    !!!!!sure that cloud layers (if more than one cloud) are sequential
        NCLDLAY=0
        iLayers = 0  !!!!!  ----- "iLayers" is an important variable -----
        iNumClds = 0
        DO i=1,kMaxClouds
            IF (iaCloudWithThisAtm(i) == 1) THEN
                iNumClds = iNumClds + 1
                iT = i
                ncldlay = ncldlay + iaCloudNumLayers(i)
                IF (kOuterLoop == 1) THEN
                    write(kStdWarn,*) 'Cloud #, num layers = ',i,iaCloudNumLayers(i)
                    write(kStdWarn,*) 'L  KCLay iscattab  iacType  Ttemp(K)   dme(um)  iwp(g/m2) '
                    write(kStdWarn,*) '----------------------------------------------------------'
                END IF
                DO iStep = 1, iaCloudNumLayers(i)
                    iLayers = iStep
                    raaIWP(iLayers,i)     = raakLayersCldAmt(iaRadLayer(iLayers),i)
                    IF ((ka100layerCloudType(I) >= 200) .AND. (ka100layerCloudType(I) <= 299)) THEN
                        raaDME(iLayers,i)     = raaaCloudParams(i,iStep,2)  !! orig
                    !!! do KN Liou parameters
                        rT = raTempLay(iStep) - 273.15
                        IF (rT < -50.0) rT = -50.0
                        IF (rT > -25.0) rT = -25.0
                        raaDME(iLayers,i) = 0.0
                        DO iC = 0,3
                            raaDME(iLayers,i) = raaDME(iLayers,i) + raC(iC)*(rT ** iC)
                        END DO
                    ELSE
                        raaDME(iLayers,i)     = raaaCloudParams(i,iStep,2)
                    END IF
                    iaaSCATTAB(iLayers,i) = iaaScatTable(i,iStep)
                    raCldLayer(iLayers)   = +1.0 * iaaCloudWhichLayers(i,iStep)
                    IF (kOuterLoop == 1) THEN
                    !                write(kStdWarn,1234) iLayers,int(raCldLayer(iLayers)),raaIWP(iLayers,i),
                    !     $                            raaDME(iLayers,i),iaaScattab(iLayers,i),
                    !     $                            iaCldTypes(i),raTempLay(iStep)
                        write(kStdWarn,1234) iLayers,int(raCldLayer(iLayers)),iaaScattab(iLayers,i), &
                        iaCldTypes(i),raTempLay(iStep),raaDME(iLayers,i),raaIWP(iLayers,i)
                    END IF
                    raCldLayer(iLayers) = +1.0 * &
                    (iT_atm - iaaCloudWhichLayers(i,iStep) + 1)
                END DO
            END IF
        END DO
    END IF
     
    1234 FORMAT(I3,' ',I3,'   ',I3,'         ',I3,' ',F10.4,' ',F10.4,' ',F10.4,' ')

!!!don't worry about "in between empty clouds
    i  = 1
    iT = iaaCloudWhichLayers(i,1)
    iB = iaaCloudWhichLayers(i,iaCloudNumLayers(i))
    write (kStdWarn,*) 'KCARTA (A) cloud layers are from ',iB,' to ',iT
! ow swap the things
    iB = iaaCloudWhichLayers(i,1)
    iT = iaaCloudWhichLayers(i,iaCloudNumLayers(i))
    write (kStdWarn,*) 'KCARTA (B) cloud layers are from ',iB,' to ',iT

    iB = iT_Atm - iB + 1
    iT = iT_Atm - iT + 1
    write (kStdWarn,*) 'RTSPEC cloud layers are from ',iB,' to ',iT

    IF (kWhichScatterCode /= 5) THEN
        ICLDTOP = iT-1
        ICLDBOT = iB
    ELSEIF (kWhichScatterCode == 5) THEN
        ICLDTOP = iT
        ICLDBOT = iB
    END IF
    IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 3)) THEN
        ICLDTOP = icldtop+1
        ICLDBOT = icldbot+1
    END IF

    IF (iDownWard > 0) THEN
        IOBS    = iNumLayer
    ELSE IF (iDownWard < 0) THEN
        IOBS   = 1
    END IF

    write(kStdWarn,*) 'icldtop,icldbot = ',icldtop,icldbot
           
! ------------- >> find the scattering params << --------------

    RETURN
    end SUBROUTINE SetMieTables_RTSPEC_100layer

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!   *** it also adds on the "backscattered" part for PCLSAM algorithm ***
! this way we have a fast alternative to kTwoStream
    SUBROUTINE AddCloud_pclsam_SunShine_100layerclouds(raFreq, &
    raaExtTemp,raaScatTemp,raaAsymTemp, &
    iaaRadLayer,iAtm,iNumlayer,iNclouds, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaScatTemp(kMaxPts,kMixFilRows)   !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    REAL :: raCC(kProfLayer)                   !cloud fraction; default to 1.0
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND
    INTEGER :: iNclouds

! mie scattering tables
    INTEGER :: NCLDLAY, ICLDTOP, ICLDBOT
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)

! local variables
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere,iG
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat

    REAL :: raaEAll(kMaxClouds,kMaxPts,kProfLayer), raE(kMaxPts)
    REAL :: raaWAll(kMaxClouds,kMaxPts,kProfLayer), raW(kMaxPts)
    REAL :: raaGAll(kMaxClouds,kMaxPts,kProfLayer), raG(kMaxPts)
    INTEGER :: iDMEVary
    REAL :: raaGasAbs(kMaxPts,kProfLayer),rX,rY

    iDMEVary = +1    !! if dme varies   across 100 layers, do things slowly
    iDMEVary = -1    !! if dme constant across 100 layers, do things fast

    DO N = 1,kProfLayer
        DO iF = 1,kMaxPts
            raaGasAbs(iF,N) = raaExtTemp(iF,N)
        END DO
    END DO

! --------------------------------->     <----------------------------------
! oddly there is an offset of 1 compared to SUBROUTINE AddCloud_pclsam :
!           iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)  there
!           iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)     here
! --------------------------------->     <----------------------------------

!!!!first find out where the cloud top is, in kCARTA layering
    N = iCldTop
    iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
    N = iCldBot-1
    iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    DO iG = 1,iNclouds
    !! read for "one" cloud layer
        N = iCldTop
        L  = N-ICLDTOP+1
        I  = iaaSCATTAB(L,iG)
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        DO iF = 1,kMaxPts
            waveno = raFreq(iF)
            CALL INTERP_SCAT_TABLE2 (WAVENO, raaDME(L,iG), &
            EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
            NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
            TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
        ! tore optical coeffs of cloud layer
            raE(iF) = EXTINCT/1000
            raW(iF) = SSALB(L)
            raG(iF) = ASYM_RTSPEC(L)
        END DO          !loop over freqs

    !! put info for all cloud layers
        DO N = 1,iNumLayer
            L  = N-ICLDTOP+1
            I  = iaaSCATTAB(L,iG)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        !          print *,iG,N,L,iI,raaIWP(L,iG),raaDME(L,iG),raE(1),raW(1),raG(1)
            DO iF = 1,kMaxPts
            !  Compute the optical depth of cloud layer, including gas
                raaEAll(iG,iF,iI) = raE(iF)*raaIWP(L,iG)*raCC(iI)  !! include raCC????
                raaEAll(iG,iF,iI) = raE(iF)*raaIWP(L,iG)           !! ignore  raCC????
                raaWAll(iG,iF,iI) = raW(iF)
                IF (raaIWP(L,iG) >= 1.0e-10) THEN
                    raaGAll(iG,iF,iI) = raG(iF)
                ELSE
                    raaGAll(iG,iF,iI) = 0.0
                END IF
            END DO    !loop over freqs
        END DO      !loop over cloud layers
    !        print *,'interp : did ',iG,' of ',iNclouds,' clouds'
    END DO        !loop over clouds

! --------------------------->
! combine the different cloud optical properties into one !!!!!!
    DO iL = 1,kProfLayer
        DO iF = 1,kMaxPts
            raaExtTemp(iF,iL) = 0.0
            raaScatTemp(iF,iL) = 0.0
            raaAsymTemp(iF,iL) = 0.0
        END DO
    END DO

! first sum over optical depths to find (weighted) total optical depth
    DO iG = 1,iNclouds
        DO iI = 1,iNumLayer
            iL = iaaRadLayer(iAtm,iI)
            DO iF = 1,kMaxPts
                raaExtTemp(iF,iL) = max(raaExtTemp(iF,iL) + raaEAll(iG,iF,iL),0.0)
            END DO
        END DO
    END DO

! now find the weighted single scattering parameter
    DO iG = 1,iNclouds
        DO iI = 1,iNumLayer
            iL = iaaRadLayer(iAtm,iI)
            DO iF = 1,kMaxPts
                rX = raaWAll(iG,iF,iL)*raaEAll(iG,iF,iL)
                rY = raaExtTemp(iF,iL)
                raaScatTemp(iF,iL) = max(raaScatTemp(iF,iL)+(rX)/(rY+1.0e-16),0.0)
            END DO
        END DO
    END DO

! now find the weighted asymmetry
    DO iG = 1,iNclouds
        DO iI = 1,iNumLayer
            iL = iaaRadLayer(iAtm,iI)
            DO iF = 1,kMaxPts
                rX = raaGAll(iG,iF,iL)*raaWAll(iG,iF,iL)*raaEAll(iG,iF,iL)
                rY = raaExtTemp(iF,iL)*raaScatTemp(iF,iL)
                raaAsymTemp(iF,iL) = raaAsymTemp(iF,iL)+(rX)/(rY+1.0e-16)
            END DO
        END DO
    END DO

! ---->> flip things for the code to work well : flip the layers!!!! <<----
    DO iG = 1,iNumLayer
        N = iG
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        iL = iaaRadLayer(iAtm,iG)
        DO iF = 1,kMaxPts
            raaEall(1,iF,iI) = raaExtTemp(iF,iL)
            raaWall(1,iF,iI) = raaScatTemp(iF,iL)
            raaGall(1,iF,iI) = raaAsymTemp(iF,iL)
        END DO
    END DO
! ---->> flip things for the code to work well : flip the layers!!!! <<----
            
! --------------------------->
! now add GAS abs coeffs + CLOUD abs coeffs
! temporarily save things in raaWAll,raaGAll,raaEall for later flip
    DO iG = 1,iNumLayer
        N = iG
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        DO iF = 1,kMaxPts
            waveno = raFreq(iF)
            taugas = raaGasAbs(iF,iI)
            rAbs   = taugas
        !  Compute the optical depth of cloud layer, including gas
            TAUC_L   = raaEall(1,iF,iI)
            TAUCG_L  = TAUGAS + TAUC_L
            TAUTOT_N = TAUCG_L
        ! the SSALB coeff
            raaScatTemp(iF,iI) = raaWall(1,iF,iI)*TAUC_L/TAUCG_L
        ! now add on the backscattered part
            b = (1.0 - raaGall(1,iF,iI))/2.0
            TAUTOT_N = TAUTOT_N * (1 - raaScatTemp(iF,iI)*(1.0-b))
            raaExtTemp(iF,iI)  = TAUTOT_N
        END DO          !loop over freqs
    END DO        !loop over cloud layers

! now use the partial fractions
    i1  = iaaRadLayer(iAtm,1)
    i2  = iaaRadLayer(iAtm,iNumLayer)
    iFixHere = -1         !!!do not adjust here, scatter_twostream does it
    iFixHere = +1         !!!do adjust here, scatter_twostream does not
    iFixHere = -1
    IF (iFixHere > 0) THEN
        IF (i1 > i2) THEN
        ! adiation going from eg layer 100 to 1 ==> up look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)   = raaExtTemp(iF,i1) * rFracTop
                raaExtTemp(iF,i2)   = raaExtTemp(iF,i2) * rFracBot
            END DO
        ELSEIF (i1 < i2) THEN
        ! adiation going from eg layer 1 to 100 ==> down look instr
            DO iF = 1,kMaxPts
                raaExtTemp(iF,i1)   = raaExtTemp(iF,i1) * rFracBot
                raaExtTemp(iF,i2)   = raaExtTemp(iF,i2) * rFracTop
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE AddCloud_pclsam_SunShine_100layerclouds

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!  1) it also adds on the "backscattered" part for PCLSAM algorithm
! this way we have a fast alternative to kTwoStream
!  2) does the jacobian part for d/d(DME)
! this is for a DOWNLOOK instrument, so we call
!        raaPhaseJacobASYM(iF,iI) = hg2_real_deriv_wrt_g(-mu_sun,mu_sat,ASYM)
    SUBROUTINE AddCloud_pclsam_Jacob_downlook_sunshine( &
    raFreq,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raaPhaseJacobASYM, &
    iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, &
    iNclouds, raaIWP, raaDME, iaaSCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! mie scattering tables
    INTEGER ::  NCLDLAY, ICLDTOP, ICLDBOT
    INTEGER ::  ISCATTAB(MAXNZ),iaaScattab(maxnz,kMaxCLouds),iNclouds, iG
    REAL ::     IWP(MAXNZ), DME(MAXNZ)
    REAL ::     raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

! local variables
    REAL :: mu_sun,mu_sat
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere,iDMEVary
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat
    REAL :: OMEGA, ASYM,tautotal_0

    REAL :: dEXTINCT_dr, dSSALB_dr, dASYM_dr
    REAL :: rW,x1,x2,x3,x4,x5
    REAL :: hg2_real

    INTEGER ::  ISCATTABX(MAXNZ)
    REAL ::     IWPX(MAXNZ), DMEX(MAXNZ)
    REAL :: raaGasAbs(kMaxPts,kProfLayer),rX,rY

    iDMEVary = +1    !! if dme varies   across 100 layers, do things slowly
    iDMEVary = -1    !! if dme constant across 100 layers, do things fast

    DO N = 1,kProfLayer
        DO iF = 1,kMaxPts
            raaGasAbs(iF,N) = raaExtTemp(iF,N)
        END DO
    END DO

! --------------------------------->     <----------------------------------
! oddly there is an offset of 1 compared to SUBROUTINE AddCloud_pclsam :
!           iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N+1)  there
!           iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)     here
! --------------------------------->     <----------------------------------

    IWP(1) = 10.0
!!!!first find out where the cloud top is, in kCARTA layering
    N = iCldTop
    iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
    N = iCldBot-1
    iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

!! ----> set DME as nonzero, assign a cloud scattering table <-----
!! ----> to show thermal struct of atm via d(rad)/d(IWP)     <-----
! o get largest cloud info
    DO iI = 1,kProfLayer
        iscattab(iI)  = -1
        dme(iI)       = -1.0
        iwp(iI)       = -1.0
        iscattabX(iI)  = -1
        dmeX(iI)       = -1.0
        iwpX(iI)       = -1.0
    END DO

    DO N = 1,iNumlayer
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        DO iG = 1,iNclouds
            i2 = 1   !!! set this default
            IF (raaIWP(L,iG) > iwpX(L)) THEN
                i2 = iG
                iwpX(L) = raaIWP(L,iG)
            END IF
        END DO
        iwpX(L) = raaIWP(L,i2)
        dmeX(L) = raaDME(L,i2)
        iscattabX(L) = iaaScattab(L,i2)
    !        print *,'(1)',N,L,iI,i2,iwpX(L),dmeX(L),iscattabX(L)
    END DO
!! ----> set DME as nonzero, assign a cloud scattering table <-----
!! ----> to show thermal struct of atm via d(rad)/d(IWP)     <-----
            
! flip things for the code to work well : flip the layers!!!!
    DO iG = 1,iNumLayer
        N  = iG
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        iL = iaaRadLayer(iAtm,iG)
    !       raaEall(1,iF,iI) = raaExtTemp(iF,iL)
        iwp(iI)      = iwpX(iG)
        dme(iI)      = dmeX(iG)
        iscattab(iI) = iscattabX(iG)
    !        print *,'(2)',N,iG,iI,iL,iwp(iI),dme(iI),iscattab(iI)
    END DO

! ow get the optical properties for the cloud layers
    DO iI = 1,kProfLayerJac
        DO iF = 1,kMaxPts
            raaPhaseJacobASYM(iF,iI) = 0.0
            raaExtJacobIWP(iF,iI)    = 0.0
            raaSSAlbJacobIWP(iF,iI)  = 0.0
            raaAsymJacobIWP(iF,iI)   = 0.0
            raaExtJacobDME(iF,iI)    = 0.0
            raaSSAlbJacobDME(iF,iI)  = 0.0
            raaAsymJacobDME(iF,iI)   = 0.0
        END DO
    END DO

    DO iG = 1,iNumLayer
        N  = iG
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        iL = iaaRadLayer(iAtm,iG)
        I  = ISCATTAB(iI)
        mu_sat = cos(raLayAngles(iI)*kPi/180)
        mu_sun = cos(raSunAngles(iI)*kPi/180)
        iF = 1
    !        print *,'(3)',N,L,iI,iL,raaGasAbs(iF,iL),ISCATTAB(iI),DME(iI),IWP(iI)
        DO iF = 1,kMaxPts
            waveno = raFreq(iF)
            taugas = raaGasAbs(iF,iL)
            rAbs   = taugas
        !  here we only need the simpler first choice as we are not messing
        !  around with the phase functions
            CALL INTERP_SCAT_TABLE2 (WAVENO, DME(iI), &
            EXTINCT, SSALB(iI), ASYM_RTSPEC(iI), &
            NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
            TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

            CALL JACOBIAN_INTERP_SCAT_TABLE2 (WAVENO, DME(iI), &
            dEXTINCT_dr, dSSALB_dr, dASYM_dr, &
            NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
            TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

            OMEGA = SSALB(iI)
            ASYM  = ASYM_RTSPEC(iI)

        !  Compute the optical depth of cloud layer, including gas
            TAUC_L   = IWP(iI)*EXTINCT/1000
            TAUCG_L  = TAUGAS + TAUC_L
            TAUTOT_N = TAUCG_L

        !   the SSALB coeff
            rW        = SSALB(iI)
            rScat     = SSALB(iI) * IWP(iI)*EXTINCT/1000
            SSALB(iI) = SSALB(iI) * TAUC_L/TAUCG_L
            raaSSAlbTemp(iF,iL) = SSALB(iI)

        ! ---------------> now add on the backscattered part <--------------------
            b = (1.0 - ASYM_RTSPEC(iI))/2.0
            TAUTOT_N = TAUTOT_N * (1 - SSALB(iI)*(1.0-b))
            raaExtTemp(iF,iL)  = TAUTOT_N
        ! ---------------> now add on the backscattered part <--------------------

            IF (IWP(iI) >= 1.0e-5) THEN
                raaAsymTemp(iF,iL) = ASYM_RTSPEC(iI)
            ELSE
                raaAsymTemp(iF,iL) = 0.0
            END IF
        ! -------------------------- now do the jacobians --------------------------
        !! technically we are doing d/d(DME) and not d/d(RME); they are
        ! related by raaXYZJacobRME(iF,iI) = raaXYZJacobDME(iF,iI)

            tautotal_0 = TAUCG_L

        !          if (iF .EQ. 1) THEN
        !            print *,iI,iwp(iI),raaGasAbs(1,iL),tautotal_0
        !            end if

        !! --------> d/d(iwp) <---------  !!
            x1 = EXTINCT/1000
            x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
            raaExtJacobIWP(iF,iL) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2
             
            x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
            raaSSAlbJacobIWP(iF,iL) = x2

            raaAsymJacobIWP(iF,iL) = 0.0

        !! --------> d/d(dme) <---------  !!
            x1 = IWP(iI)/1000*dEXTINCT_dr
            x4 = EXTINCT*IWP(iI)/1000/TAUCG_L
            x5 = tautotal_0*SSALB(iI)*dEXTINCT_dr*(1-x4)
            x2 = IWP(iI)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
            x3 = -1/2*dASYM_dr
            raaExtJacobDME(iF,iL) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2 + &
            TAUCG_L*SSALB(iI)*x3

            x4 = EXTINCT*IWP(iI)/1000/TAUCG_L
            x5 = tautotal_0*SSALB(iI)*dEXTINCT_dr*(1-x4)
            x2 = IWP(iI)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
            raaSSAlbJacobDME(iF,iL) = x2

            raaAsymJacobDME(iF,iI) = dASYM_dr

        !! --------> d/d(g) <---------  !!
            raaPhaseJacobASYM(iF,iL) = &
            hg2_real_deriv_wrt_g(-mu_sun,mu_sat,ASYM)

        END DO          !loop over freqs
    END DO        !loop over cloud layers

! now use the partial fractions????? see last section in
!       SUBROUTINE AddCloud_pclsam( )

    RETURN
    end SUBROUTINE AddCloud_pclsam_Jacob_downlook_sunshine

!************************************************************************
! this subroutine adds on the absorptive part of cloud extinction
! basically the same as AddCloud_twostream EXCEPT
!  1) it also adds on the "backscattered" part for PCLSAM algorithm
! this way we have a fast alternative to kTwoStream
!  2) does the jacobian part for d/d(DME)
! this is for a UPLOOK instrument, so we call
!       raaPhaseJacobASYM(iF,iI) = hg2_real_deriv_wrt_g(-mu_sun,-mu_sat,ASYM)
    SUBROUTINE AddCloud_pclsam_Jacob_uplook_sunshine( &
    raFreq,raLayAngles,raSunAngles, &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
    raaPhaseJacobASYM, &
    iaaRadLayer,iAtm,iNumlayer, &
    rFracTop,rFracBot, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    NCLDLAY, ICLDTOP, ICLDBOT, &
    iNclouds, raaIWP, raaDME, iaaSCATTAB, &
    NSCATTAB, MUINC, &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, &
    TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! usual variables
    INTEGER :: iAtm,iNumlayer                  !which atmosphere, num of layers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer) !to get layer info
    REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
    REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
    REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
    REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
    REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
    REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
    REAL :: raFreq(kMaxPts)                    !wavenumber grid
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA    !kcarta cloud top/bottoms
    REAL :: rFracTop,rFracBot                  !layer fractions at TOA,GND
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! mie scattering tables
    INTEGER ::  NCLDLAY, ICLDTOP, ICLDBOT
    INTEGER ::  ISCATTAB(MAXNZ),iaaScattab(maxnz,kMaxCLouds),iNclouds, iG
    REAL ::     IWP(MAXNZ), DME(MAXNZ)
    REAL ::     raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
    INTEGER ::  NSCATTAB
    INTEGER ::  NMUOBS(NSCATTAB), NDME(NSCATTAB), NWAVETAB(NSCATTAB)
    REAL ::     MUTAB(MAXGRID,NSCATTAB)
    REAL ::     DMETAB(MAXGRID,NSCATTAB), WAVETAB(MAXGRID,NSCATTAB)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,NSCATTAB), TABSSALB(MAXTAB,NSCATTAB)
    REAL ::     TABASYM(MAXTAB,NSCATTAB)
    REAL ::     TABPHI1UP(MAXTAB,NSCATTAB), TABPHI1DN(MAXTAB,NSCATTAB)
    REAL ::     TABPHI2UP(MAXTAB,NSCATTAB), TABPHI2DN(MAXTAB,NSCATTAB)
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
    REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

! local variables
    REAL :: mu_sun,mu_sat
    INTEGER :: iL,iF,iI,N,L,I,ikcCldtop,ikcCldbot
    INTEGER :: i1,i2,iFixHere
    REAL :: tauc_L,taucg_L,tautot_n,taugas,waveno,b
    REAL :: extinct,SSALB(MAXNZ), ASYM_RTSPEC(MAXNZ)
    REAL :: dmedme,albedo,asymmetry,rAbs,rAlbedo,rScat
    REAL :: OMEGA, ASYM,tautotal_0

    REAL :: dEXTINCT_dr, dSSALB_dr, dASYM_dr
    REAL :: rW,x1,x2,x3,x4,x5
    REAL :: hg2_real

    REAL :: raaGasAbs(kMaxPts,kProfLayer),rX,rY

    DO N = 1,kProfLayer
        DO iF = 1,kMaxPts
            raaGasAbs(iF,N) = raaExtTemp(iF,N)
        END DO
    END DO

    IWP(1) = 10.0
    IF (IWP(1) > 0.0) THEN
    !!!!first find out where the cloud top is, in kCARTA layering
        N = iCldTop
        iKcCldTop = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        N = iCldBot-1
        iKcCldBot = IFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)

    ! o get largest cloud info
        DO iI = 1,kProfLayer
            iscattab(iI) = -1
            dme(iI)      = -1.0
            iwp(iI)      = -1.0
        END DO
        DO N = 1,iNumlayer
            L  = N-ICLDTOP+1
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
            DO iG = 1,iNclouds
                i2 = 1   !!! set this default
                IF (raaIWP(L,iG) > iwp(L)) THEN
                    i2 = iG
                    iwp(L) = raaIWP(L,iG)
                END IF
            END DO
            iwp(L) = raaIWP(L,i2)
            dme(L) = raaDME(L,i2)
            iscattab(L) = iaaScattab(L,i2)
        END DO

    ! ow get the optical properties for the cloud layers
        DO iI = 1,kMixFilRows
            DO iF = 1,kMaxPts
                raaPhaseJacobASYM(iF,iI) = 0.0
                raaExtJacobIWP(iF,iI)    = 0.0
                raaSSAlbJacobIWP(iF,iI)  = 0.0
                raaAsymJacobIWP(iF,iI)   = 0.0
                raaExtJacobDME(iF,iI)    = 0.0
                raaSSAlbJacobDME(iF,iI)  = 0.0
                raaAsymJacobDME(iF,iI)   = 0.0
            END DO
        END DO

        DO N = 1,iNumLayer
            L  = N-ICLDTOP+1
            I  = ISCATTAB(L)
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
            mu_sat = cos(raLayAngles(iI)*kPi/180)
            mu_sun = cos(raSunAngles(iI)*kPi/180)
            DO iF = 1,kMaxPts
                waveno = raFreq(iF)
                taugas = raaGasAbs(iF,iI)
            !            taugas = raaExtTemp(iF,iI)
                rAbs   = taugas
            !  here we only need the simpler first choice as we are not messing
            !  around with the phase functions
                CALL INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                EXTINCT, SSALB(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                CALL JACOBIAN_INTERP_SCAT_TABLE2 (WAVENO, DME(L), &
                dEXTINCT_dr, dSSALB_dr, dASYM_dr, &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

                OMEGA = SSALB(L)
                ASYM  = ASYM_RTSPEC(L)

            !  Compute the optical depth of cloud layer, including gas
                TAUC_L   = IWP(L)*EXTINCT/1000
                TAUCG_L  = TAUGAS + TAUC_L
                TAUTOT_N = TAUCG_L

            !   the SSALB coeff
                rW       = SSALB(L)
                rScat    = SSALB(L) * IWP(L)*EXTINCT/1000
                SSALB(L) = SSALB(L) * TAUC_L/TAUCG_L
                raaSSAlbTemp(iF,iI) = SSALB(L)

            ! ---------------> now add on the backscattered part <--------------------
                b = (1.0 - ASYM_RTSPEC(L))/2.0
                TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b))
                raaExtTemp(iF,iI)  = TAUTOT_N
            ! ---------------> now add on the backscattered part <--------------------

                IF (IWP(L) >= 1.0e-5) THEN
                    raaAsymTemp(iF,iI) = ASYM_RTSPEC(L)
                ELSE
                    raaAsymTemp(iF,iI) = 0.0
                END IF
            ! -------------------------- now do the jacobians --------------------------
            !! technically we are doing d/d(DME) and not d/d(RME); they are
            ! related by raaXYZJacobRME(iF,iI) = raaXYZJacobDME(iF,iI)

                tautotal_0 = TAUCG_L

            !! --------> d/d(iwp) <---------  !!
                x1 = EXTINCT/1000
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaExtJacobIWP(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2
                 
                x2 = OMEGA*EXTINCT/1000*TAUGAS/(TAUCG_L**2)
                raaSSAlbJacobIWP(iF,iI) = x2

                raaAsymJacobIWP(iF,iI) = 0.0

            !! --------> d/d(dme) <---------  !!
                x1 = IWP(L)/1000*dEXTINCT_dr
                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                x3 = -1/2*dASYM_dr
                raaExtJacobDME(iF,iI) = TAUTOT_N/TAUCG_L*x1 + TAUCG_L*(b-1)*x2 + &
                TAUCG_L*SSALB(L)*x3

                x4 = EXTINCT*IWP(L)/1000/TAUCG_L
                x5 = tautotal_0*SSALB(L)*dEXTINCT_dr*(1-x4)
                x2 = IWP(L)/1000*x5/(TAUCG_L**2) + x4*dSSALB_dr
                raaSSAlbJacobDME(iF,iI) = x2

                raaAsymJacobDME(iF,iI) = dASYM_dr

            !! --------> d/d(g) <---------  !!
                raaPhaseJacobASYM(iF,iI) = &
                hg2_real_deriv_wrt_g(-mu_sun,-mu_sat,ASYM)

            END DO          !loop over freqs
        END DO        !loop over cloud layers
    ENDIF

! now use the partial fractions????? see last section in
!       SUBROUTINE AddCloud_pclsam( )

    RETURN
    end SUBROUTINE AddCloud_pclsam_Jacob_uplook_sunshine

!************************************************************************

END MODULE clear_scatter_misc
