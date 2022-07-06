! Copyright 1997
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

MODULE scatter_disort_aux

USE basic_common
USE ttorad_common
USE kcoeff_common
USE spline_and_sort_and_common
USE rad_diff_and_quad
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_main
USE ttorad_common
USE spline_and_sort_and_common

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subtroutine does some initializations for a uplook instr
    SUBROUTINE Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,rSatAngle,iTag, &
    raTopIntensity,raSolarBeam,TOA_to_instr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input variables
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iTag
    REAL :: raVTemp(kMixFilRows)        !temperature profile
    REAL :: rFracTop                    !topmost fraction
    REAL :: raFreq(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs
    REAL :: rSatAngle
! output variables
    REAL :: raTopIntensity(kMaxPts)     !intensity at TOA
    REAL :: raSolarBeam(kMaxPts)        !intensity diue to solar rad
    REAL :: TOA_to_instr(kMaxPts)       !abs coeffs

! local variables
    INTEGER :: iaRadLayer(kProfLayer),iL,iF
    REAL :: rDummy

    DO iL=1,kMaxPts
      TOA_to_instr(iL) = 0.0
    END DO

! bring incident space radiation down from TOA to instrument
    DO iF=1,kMaxPts
      ! compute the Plank radiation from space
      raTopIntensity(iF) = ttorad(raFreq(iF),sngl(kTSpace))
    END DO

! set the solar beam intensity ...
! if the sun is ON and is in the FOV of the instrument, then the BC of 5700K
!   to the TOA is set, else the BC of 2.6K to TOA is set; also fbeam set to 0
! if the sun is ON and is NOT in the FOV of the instrument, then the fbeam is
!   set to solar beam, while BC of TOA is 2.6K is set

    DO iF=1,kMaxPts     !!!!assume sun is NOT ON
      raSolarBeam(iF) = 0.0
    END DO

    IF (kSolar >= 0) THEN
      IF (abs(abs(rSatAngle)-abs(kSolarAngle)) >= 1.0e-3) THEN
        ! sun is on, but not in FOV of instr
        ! set fbeam correctly to that of sun
        CALL SolarBeamDisort(kSolar,raSolarBeam,raFreq,iTag)
        rDummy = abs(cos(kSolarAngle*kPi/180))
        ! bring solar radiation down from TOA to instrument
        CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
            rFracTop,raFreq,raaAbs,TOA_to_instr)
        DO iF=1,kMaxPts
          raSolarBeam(iF) = raSolarBeam(iF)*exp(-TOA_to_instr(iF)/rDummy)
        END DO
      END IF

      IF (abs(abs(rSatAngle)-abs(kSolarAngle)) <= 1.0e-3) THEN
        ! un is on, and in FOV of instr
        ! o set fbeam correctly to 0.0, and the BC of 5700K
      END IF
    END IF     !!IF (kSolar >= 0)

    RETURN
    end SUBROUTINE Init_UpLook
!************************************************************************
! this subtroutine does some initializations for a downlook instr
    SUBROUTINE Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,rSatAngle,iTag, &
    raTopIntensity,raSolarBeam,TOA_to_instr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input variables
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iTag
    REAL :: raVTemp(kMixFilRows)        !temperature profile
    REAL :: rFracTop                    !topmost fraction
    REAL :: raFreq(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs
    REAL :: rSatAngle
! output variables
    REAL :: raTopIntensity(kMaxPts)     !intensity at TOA
    REAL :: raSolarBeam(kMaxPts)        !intensity diue to solar rad
    REAL :: TOA_to_instr(kMaxPts)       !abs coeffs

! local variables
    INTEGER :: iaRadLayer(kProfLayer),iL,iF
    REAL :: rDummy

! see if there are any layers between TOA and instr; if there are, set
! TOA_to_instr to the cumulative k(TOA to aircraft), else set it to 0
    DO iL=1,kProfLayer
      iaRadLayer(iL) = iaaRadLayer(iAtm,iL)
    END DO
    CALL Find_K_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
      rFracTop,raFreq,raaAbs,TOA_to_instr)

! bring incident space radiation down from TOA to instrument
    DO iF=1,kMaxPts
      raTopIntensity(iF) = ttorad(raFreq(iF),sngl(kTSpace))
    END DO
! this is technically incorrect ... we really should do rad transfer here
    DO iF=1,kMaxpts
      raTopIntensity(iF) = raTopIntensity(iF)*exp(-TOA_to_instr(iF))
    END DO
          
    IF (kSolar >= 0) THEN
      CALL SolarBeamDisort(kSolar,raSolarBeam,raFreq,iTag)
      rDummy = abs(cos(kSolarAngle*kPi/180))
      ! bring solar radiation down from TOA to instrument
      DO iF=1,kMaxPts
        raSolarBeam(iF) = raSolarBeam(iF)*exp(-TOA_to_instr(iF)/rDummy)
      END DO
    ELSE
      DO iF=1,kMaxPts
        raSolarBeam(iF) = 0.0
      END DO
    END IF

    RETURN
    end SUBROUTINE Init_DownLook
!************************************************************************
! set raCorrelatedK to abs coeffs of layer nearest the ground
    SUBROUTINE SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs, &
    iaaRadLayer,iAtm,iNumLayer, ABSPROF)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: iDownWard            !direction of rad transfer
    REAL :: raCorrelatedK(kMaxPts)  !to see how the "k distributions" are
    REAL :: raaAbs(kMaxPts,kMixFilRows) !matrix of abs coeffs
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    REAL ::    ABSPROF(MAXNZ,MAXABSNU)

    INTEGER :: iI,iJ,iaIndx(kMaxPts),iMethod
    REAL :: raL2S(kMaxPts),raaWgt(kMaxPts,kProfLayer),raPeak(kMaxPts)

    iMethod = -1

    IF (iMethod == -1) THEN
      ! store optical depth of layer closest to gnd
      IF (iDownWard == 1) THEN
        !! as k(gnd) increases, down look instr only penetrates less into
        !! atmosphere, so radiance decreases
        DO iI=1,kMaxPts
           raCorrelatedK(iI) = raaAbs(iI,iaaRadLayer(iAtm,1))
        END DO
      ELSE
        !! as k(gnd) increases, uplook instr only penetrates less into
        !! atmosphere, so radiance increases as you see hot stuff in
        !! your face
        DO iI=1,kMaxPts
          raCorrelatedK(iI) = raaAbs(iI,iaaRadLayer(iAtm,iNumLayer))
        END DO
      END IF
    END IF

    IF (iMethod == +1) THEN
      ! store optical depth of layer that peaks the weight fcn
      !!!!!!remember ABSPROF(1,:) = TOA, ABSPROF(NLEV-1,:) = GND
      IF (iDownward == 1) THEN
        ! down look weight fcn
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
     ELSE IF (iDownward == -1) THEN
       ! do up look weight fcn
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
    end SUBROUTINE SetCorrelatedK

!************************************************************************
! set up the single scatter albedos etc for the clouds
    SUBROUTINE SetUpCloudsDISORT(nstr, nmuobs, iaCloudWithThisAtm, &
    iaCldTop,iaCldBot,iaCloudNumLayers,rF, &
    iAtm,iaaRadLayer,iNumLayer, &
    IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM, ISCATTAB, &
    extinct,dtauc,ssalb,asym,pmom,iFF)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! inputs
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iFF
    REAL :: rF                                 !wavenumber
    INTEGER :: iaCloudWithThisAtm(kMaxClouds)  !is this cloud in this atm?
    INTEGER :: IACLDTOP(kMaxClouds)            !cloud top layer
    INTEGER :: IACLDBOT(kMaxClouds)            !cloud bot layer
    INTEGER :: iaCloudNumLayers(kMaxClouds)    !number of layers cloud occupies
    REAL ::    IWP(MAXNZ), DME(MAXNZ)          !iwp, particle size
    INTEGER :: NDME(MAXSCAT), NWAVETAB(MAXSCAT),ISCATTAB(MAXNZ)
    INTEGER :: nstr                            !number of disort streams
    INTEGER :: nmuobs                          !number of RTSPEC computed angles
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
! outputs
    REAL :: extinct                      !extinct coeff from Mie Scatter tables
    DOUBLE PRECISION :: dtauc(maxcly)    !optical depths; also used as INPUT
    DOUBLE PRECISION :: ASYM(maxnz)
    DOUBLE PRECISION :: pmom(0:maxmom,maxcly)  !scattering phase fcn
    DOUBLE PRECISION :: ssalb(maxcly)   !single scatter albedo of layers

! local variables
    INTEGER :: LL,M,ICLDTOP,ICLDBOT,N,L,I,nmom, iiDiv,Nprime
    REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
    DOUBLE PRECISION :: tauC(kProfLayer),tauCG(kProfLayer)
    CHARACTER(1) :: caScale(MAXSCAT)
    INTEGER :: iExtScaling,iSartaTables,iDefault
    REAL :: rE,rW,rG
!>>>>>>>>>>>>>>>>>>>>>>>>>
    iDefault = +1  !! use     the SARTA Tables as is (delta scaled)
    iDefault = -1  !! unscale the SARTA Tables as is (no delta scaled)

    iSartaTables = +1   !! SARTA scattering tables,  delta scaled
    iSartaTables = -1   !! RRTM  does not use delta scaling
    iSartaTables = iaaOverrideDefault(3,7)

    IF ((iSartaTables .NE. iDefault) .AND. (kOuterLoop .EQ. 1) .AND. (iFF .EQ. 1)) THEN
      write(kStdErr,*) 'OH OH KEEP DELTA SCALE for SetUpCloudsDISORT'
      write(kStdErr,'(A,I3,I3,A,I3,I3)') 'SetUpCloudsDISORT iDefault,iExtScaling = ',&
                   iDefault,iSartaTables,' kWhichScatterCode,kScatter = ',kWhichScatterCode,kScatter
    END IF
!>>>>>>>>>>>>>>>>>>>>>>>>>

    iiDiv = 0
    555 CONTINUE
    IF (iiDiv*kProfLayer < iaaRadLayer(iAtm,3)) THEN
      iiDiv = iiDiv + 1
      GOTO 555
    END IF
    iiDiv = iiDiv - 1

!!!!!!!!!!!!!! ********* CLOUD SCATTERING ************* !!!!!!!!!!!!
    LL = 0
    DO M = 1,kMaxClouds
      IF (iaCloudWithThisAtm(M) > 0) THEN
        ICLDTOP = IACLDTOP(M)+1
        ICLDBOT = IACLDBOT(M)+1
        DO N = ICLDTOP, ICLDBOT-1
          Nprime = N-iiDiv*kProfLayer
          ! L is the cloud layer = 1(top)..nlay(bot) in cloud
          L = LL + N-ICLDTOP+1

!          write (kstdErr,'(A,3(I3),2(F12.6))') 'm,tp,bt,l,iwp(l),dme(l) = ',&
!               m,icldtop,icldbot,l,iwp(l),dme(l)

          I = ISCATTAB(L)     !!!!!I is the scattering table info number

          ! Interpolate to get values of extinction, s.s. albedo, and
          ! phi function values at given obs. mu, waveno, and particle size.
          ! Note: we don't need phi functions for Eddington-only solution
          ! This means that while the rtspec code had a choice of
          !       CALL INTERP_SCAT_TABLE2 (rF, DME(L),    versus
          !       CALL INTERP_SCAT_TABLE3 (rF, DME(L),
          ! here we only need the simpler first choice as we are not messing
          ! around with the phase functions
          CALL INTERP_SCAT_TABLE2 (rF, DME(L), &
                EXTINCT, SSALB_RTSPEC(L), ASYM_RTSPEC(L), &
                NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
                TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))

!   this was commmented out till June 2022 then Tang/Yang/Huang/RRTM suggest no use unscaled
!   Frank Evans code scales the Mie scattering parameters, so we have to
!   unscale them!!!!!!!!
          IF (iSartaTables .LT. 0) THEN
            rE = EXTINCT
            rW = SSALB_RTSPEC(L)
            rG = ASYM_RTSPEC(L)
            CALL UnScaleMieOne('G',rE,rW,rG)
            EXTINCT = rE
            SSALB_RTSPEC(L) = rW
            ASYM_RTSPEC(L) = rG

!            IF (iFF .EQ. 1) THEN
!              write(kStdErr,'(A,2(F12.5),4(ES12.5))') &
!              ' DISORT iSartaTables = 1 : UNSCALED rF, rDME, rE rW rG IWP(L)',rF,DME(L),EXTINCT/1000,SSALB_RTSPEC(L),ASYM_RTSPEC(L),IWP(L)
!            END IF

          END IF

          ! but the indices have to be modified so we put info into the correct layers
          ! also have to set asymmetry parameter
          TAUC(L) = DBLE(IWP(L)*EXTINCT/1000.0)
          !!!!!!!!!!!! orig code TAUCG(L) = TAUGAS(N) + TAUC(L)
          TAUCG(L) = dtauc(Nprime) + TAUC(L)

          IF (TAUCG(L) > 1.0D-5) THEN
            !!! rescale SSA to be wc Ecld/ Etotal
            SSALB_RTSPEC(L) = SSALB_RTSPEC(L)*TAUC(L)/TAUCG(L)
          ELSE
            SSALB_RTSPEC(L) = 0.0
          ENDIF

!    write(kStdErr,'(A,3(I4),F12.4,5(ES12.4))') 'OIOIOI',L,ICLDTOP,ICLDBOT,rF,IWP(L),EXTINCT/1000.0,TAUCG(L),SSALB_RTSPE!C(L),asym_rtspec(L)

          ! so now set the indices correctly and do relevant changes of real --> double
          !!!orig code dtauc(N)      = blah, ssalb(N)      = blah
          !!!onew code dtauc(Nprime) = blah, ssalb(Nprime) = blah
          dtauc(Nprime) = DBLE(TAUCG(L))
          SSALB(Nprime) = DBLE(SSALB_RTSPEC(L))
          IF (ssalb_rtspec(L) > 1.0e-6) THEN
            asym(Nprime)  = DBLE(asym_rtspec(L))
          ELSE
            asym(Nprime)  = DBLE(0.0)
          END IF
                  
          !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
          nmom = max(2*nmuobs + 1,nstr)
          nmom = min(maxmom,nmom)
          !! get the phase moments for the HG phase function
          CALL dgetmom(3,asym(Nprime),nmom,pmom(0,Nprime))
        ENDDO                 !DO N=ICLDTOP,ICLDBOT-1
      LL = LL + iaCloudNumLayers(M)
      END IF                  !IF (iaCloudWithThisAtm(M) > 0) THEN
    END DO                    !DO M = 1,kMaxClouds

    RETURN
    end SUBROUTINE SetUpCloudsDISORT

!************************************************************************
! this subroutine sets up the optical depths, single scatter albedo etc for
! a atmosphere where there is only gas + Rayleigh scattering
    SUBROUTINE SetUpRayleigh(nlev,nstr,nmuobs, rF,raDensity,raThickness, &
      dtauc,ssalb,asym,pmom)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
! inputs
    INTEGER :: nlev,nstr          !number of levels, number of streams
    INTEGER :: nmuobs             !number of RTSPEC angles for Mie computations
    REAL :: rF                    !wavenumber
    REAL :: raDensity(kProfLayer),raThickness(kProfLayer)
! outputs
    DOUBLE PRECISION :: dtauc(maxcly)    !optical depths; also used as INPUT
    DOUBLE PRECISION :: ASYM(maxnz)
    DOUBLE PRECISION :: pmom(0:maxmom,maxcly)  !scattering phase fcn
    DOUBLE PRECISION :: ssalb(maxcly)   !single scatter albedo of layers

! local variables
    REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
    DOUBLE PRECISION :: tauC(kProfLayer),tauCG(kProfLayer)
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
      nmom = max(2*nmuobs + 1,nstr)
      nmom = min(maxmom,nmom)
      !! get the phase moments for Rayleigh scattering
      CALL dgetmom(2,asym(N),nmom,pmom(0,N))
    ENDDO

    RETURN
    end SUBROUTINE SetUpRayleigh

!************************************************************************
! this does the final initializations before calling DISORT
    SUBROUTINE FinalInitialization( &
!!!!inputs
    iDownWard,rSatAngle,rTopIntensity,rSolarBeam,emiss, &
    rSurfaceTemp,dtauc,dTotalOpticalDepth,iDoFlux,nlev, &
    iNp,iaOp, &
!!!!outputs
    usrtau,ntau,utau,usrang,numu,umu, &
    nphi,phi,fisot,fbeam,umu0,phi0, &
    ibcnd,lamber,albedo, &
    btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
     
! input variables
    INTEGER :: nlev,iNp,iaOp(kPathsOut)
    DOUBLE PRECISION :: dtauc(maxcly)    !optical depths; also used as INPUT
    DOUBLE PRECISION :: dTotalOpticalDepth
    REAL :: rSatAngle,rTopIntensity,rSolarBeam,emiss,rSurfaceTemp
    INTEGER :: iDownWard,iDoFlux
! output variables
    INTEGER :: ntau,numu,nphi
    LOGICAL :: usrtau,usrang
    DOUBLE PRECISION :: utau(maxulv)     ! tau's at which to output results
    DOUBLE PRECISION :: umu(maxumu)      ! ang's at which to output results
    DOUBLE PRECISION :: phi(maxphi)      ! azimuthal phi's to output radiance
    DOUBLE PRECISION :: fisot            ! isotropic TOA radiation
    DOUBLE PRECISION :: fbeam,umu0,phi0  ! solar beam

    INTEGER :: ibcnd
    LOGICAL :: lamber,plank,onlyfl
    DOUBLE PRECISION :: albedo,btemp,ttemp,temis,accur
    CHARACTER  header*127                ! dumb comment
    LOGICAL :: prnt(5)                   ! print(1) = true, print input variables

! local variables
    INTEGER :: iI,iJ
    DOUBLE PRECISION :: d1,d2,dCumulative(maxcly)
         
    d1 = dtauc(1)
    d2 = dTotalOpticalDepth

    IF (iDownward == 1) THEN
      !down look instrument
      usrtau = .TRUE.   !return intensity at ONE optical depth instead of at all levels
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
        ! we will get out radiance at TOP of layer
        utau(iI) = dCumulative((nlev-1)-iaOp(iI)+1) - d1
      END DO
          
      usrang = .TRUE.   !return intensity at one user polar angle
      numu = 1
      umu(1) = DBLE(ABS(cos(rSatAngle*kPi/180)))
      !nphi = 1          !return intensity at one user zenith angle == 0
      !phi(1) = DBLE(0.0)
      nphi = 1          !return intensity at one user zenith angle == 0
      !phi(1) = DBLE(0.0)

    ELSE
      ! up look instrument
      usrtau = .TRUE.   !return intensity at ONE optical depth instead of at all levels
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
      ! we will get out radiance at BOTTOM of layer; flip iaOp stuff
          utau(iI) = dCumulative(iaOp(iI)) - d1
      END DO

      usrang = .TRUE.   !return intensity at one user polar angle
      numu = 1
      umu(1) = DBLE(-ABS(cos(rSatAngle*kPi/180)))
      !nphi = 1          !return intensity at one user zenith angle == 0
      !phi(1) = DBLE(0.0)
      nphi = 0          !return intensity at one user zenith angle == 0
      !phi(1) = DBLE(0.0)
    END IF

    fisot = 0.0
    fisot = DBLE(rTopIntensity)

    IF (abs(abs(rSatAngle) - abs(kSolarAngle)) <= 1e-7) THEN
        kSolarAngle = kSolarAngle+1e-4
    END IF
    fbeam  = DBLE(rSolarBeam)
    umu0   = DBLE(cos(kSolarAngle*kPi/180))
    phi0   = DBLE(0.0)

    ibcnd = 0                   ! general case
    lamber = .FALSE.            ! specify bidir reflectance
    lamber = .TRUE.             ! isotropic reflect lower bdry ==> specify albedo
    albedo  = DBLE(1.0-emiss)   ! from defns in books,papers

    btemp  = DBLE(rSurfaceTemp)      !!ground temperature
    IF (iDownward == 1) THEN
      !down look instrument
      ttemp  = DBLE(kTSpace)         !!2.7 Kelvin
    ELSE
      !up look instrument
      IF (kSolar < 0) THEN
        ttemp  = DBLE(kTSpace)       !!2.7 Kelvin or 5700K!!!!!!
      ELSE IF (kSolar >= 0) THEN
        IF (abs(abs(kSolarAngle)-abs(rSatAngle)) <= 1.0e-3) THEN
            ttemp  = DBLE(kSunTemp)    !!sun in FOV .. 5700K!!!!!!
        ELSE
            ttemp  = DBLE(kTSpace)     !!sun not in FOV .. 2.7!!!!!!
        END IF
      END IF
    END IF

    temis  = DBLE(1.0)
    plank  = .TRUE.      !emission from layers

    IF (iDoFlux == -1) THEN
      onlyfl = .FALSE.     !only need intensity as well
      header  = 'scattering computations'
    ELSEIF (iDoFlux == +1) THEN
      onlyfl = .TRUE.     !compute flux only
      usrtau = .FALSE.    !compute at every boundary
      ntau = nlev
      numu = 0
      nphi = 0
      header  = 'flux computations'
    END IF

    accur   =  DBLE(0.000)
    accur   =  DBLE(0.000001)
    prnt(1) = .FALSE. 
    prnt(2) = .FALSE. 
    prnt(3) = .FALSE. 
    prnt(4) = .FALSE. 
    prnt(5) = .FALSE. 

    RETURN
    end SUBROUTINE FinalInitialization

!************************************************************************
! this subroutine figures out how to do the interpolations
    SUBROUTINE Interpolator(raFreqStep,rakStep,raIntenStep, &
    raNoScatterSTEP,raaNoScatterStep,iii, &
    iStepPts,iDownWard,nlev, &
    iProfileLayers,raPressLevels, &
    raCorrelatedK,raLayerTemp,absprof,raFreq, &
    rSatAngle,rSolarAngle, &
    rSurfaceTemp,rSurfPress,raUseEmissivity, &
    raInten)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! inputs
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFreqStep(kMaxPts),raKstep(kMaxPts),raIntenStep(kMaxPts)
    REAL :: raNoScatterSTEP(kMaxPts),raaNoScatterStep(kProfLayer,kMaxPts)
    REAL :: raCorrelatedK(kMaxPts),raFreq(kMaxPts),raLayerTemp(kProfLayer)
    REAL :: ABSPROF(MAXNZ,MAXABSNU),rSatAngle,rSolarAngle,rSurfaceTemp
    REAL :: raUseEmissivity(kMaxPts),rSurfPress
    INTEGER :: iStepPts,nlev,iDownWard,iProfileLayers,iii
! output
    REAL :: raInten(kMaxPts)
          
! local variables
    REAL :: DWA(kMaxPts),raTau(kProfLayer),raNoScatterAll(kMaxPts)
    INTEGER :: Indx(kMaxPts),iF,iL

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
        DO iF = 1,kMaxPts
          raNoScatterAll(iF) = raaNoScatterStep(iii,iF)
        END DO
        CALL INTERP_PLANCK_0( &
          raFreqStep,raIntenStep,raNoScatterStep,iStepPts, &
          raFreq,raNoScatterAll,raInten)
      ELSEIF (iDownWard == 1) THEN   !!!works better for downlook instr
        DO iF = 1,kMaxPts
          DO iL = 1,NLEV-1
            raTau(iL)  = absprof(iL,iF)
          END DO
          CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
                rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
                raUseEmissivity(iF),raFreq(iF),raNoScatterAll(iF),-1, &
                iProfileLayers,raPressLevels)
        END DO
      CALL INTERP_PLANCK_0( &
            raFreqStep,raIntenStep,raNoScatterStep,iStepPts, &
            raFreq,raNoScatterAll,raInten)
      END IF
    END IF

    IF ((kScatter >= 2) .AND. (iDownWard == -1)) THEN
      !!!first interpolate in k
      DO iF = 1,kMaxPts
        CALL LINEAR_ONE(raKSTEP,raIntenStep,iStepPts,raCorrelatedK(iF),raInten(iF))
      END DO

    ELSEIF ((kScatter >= 2) .AND. (iDownWard == 1)) THEN
      !!!interpolate in k
      !!!do cumulative distr function to see if things might be too spiky
      DO iF=1,kMaxPts
        DO iL=1,NLEV-1
          raTau(iL)  = absprof(iL,iF)
          END DO
        CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
            rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
            raUseEmissivity(iF),raFreq(iF),raNoScatterAll(iF),-1, &
            iProfileLayers,raPressLevels)
        END DO
        CALL INTERP_PLANCK_3( &
          raFreqStep,raKStep,raIntenStep,raNoScatterStep,iStepPts, &
          raFreq,raNoScatterAll,raCorrelatedK,raInten)
    END IF

    RETURN
    end SUBROUTINE Interpolator
!************************************************************************
! this subroutine takes in the input abs coeffs (raaAbsCoeff) where raa(1,:)
! is the lowest layer and raa(kProfLayer,:) is the highest layer .. it then
! outputs these  abs coeffs into absprof, where absprof(1,:) is the top, and
! absprof(iNumLayer,:) = ground
!
! remember kCARTA has -----------------
!                      layer 100 = TOA
!                     -----------------
!
!                           ...
!
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
    SUBROUTINE GetAbsProfileDISORT(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
    NABSNU, NLEV, TEMP, ABSPROF, &
    ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity, &
    raDensity,raLayerTemp, &
    iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh, &
    rSatAngle,raLayAngles)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! these are variables that come in from kcartamain.f
    REAL :: raaAbs(kMaxPts,kMixFilRows),raFreq(kMaxPts),rFracTop,rFracBot
    REAL :: raVTemp(kMixFilRows),rSurfaceTemp,rSurfPress
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNpmix
    INTEGER :: iDownWard,iProfileLayers
    REAL :: iwp(MAXNZ), raNumberDensity(kProfLayer),raLayerTemp(kProfLayer)
! these are variables that we have to set
    INTEGER ::  NABSNU, NLEV         !!!!!!!!!!!!!MAXNZ, MAXABSNU
    REAL ::   TEMP(*), ABSPROF(MAXNZ,*)
    INTEGER ::  ICLDTOP,iCLDBOT,IOBS
    REAL :: raDensity(*)
    REAL :: raThicknessRayleigh(kProfLayer)
    REAL :: raLayAngles(kProfLayer),rSatAngle

! local variables
    INTEGER :: iaRadLayer(kProfLayer), iaTemp(kProfLayer),iFr, iL, iLay, iiDiv
    REAL :: NU, raVT1(kMixFilRows)

! these are to flip the temperature, abs profiles if instr looks up
    REAL :: raTemp(kProfLayer+1),raaTemp(kProfLayer,kMaxPts),rT2

!!! SUBROUTINE FinalInitialization uses
!!!      usrang = .TRUE.   !return intensity at one user polar angle
!!!      numu = 1
!!!      umu(1) = DBLE(ABS(cos(rSatAngle*kPi/180)))
!!!                Remember, can set iSnell = 0 (iaaOveride(2,7)) and the NO ray tracing due to curvature
!!!                is done. SO then PCLSAM angls will equal DISORT angles at all layers
!!!
!!!                Unfortunately DISORT does it Gaussian streams probably assuming nadir ODs are
!!!                sent in, with plane parallel assumptions
!!! SO THIS BRILLIANT PLAN IS TOAST FOR NOW
!!!                In the future I plan to modify GetAbsProfileDISORT so the layerangles are input. Then
!!!                the NABSNU ODs can be modified so that instead of being nadir ODs and then adjusted by DISORT
!!!                using umu, they will be slightly adjusted with a factor L that changes layer by layer 
!!!                   exp(-k/cos(a+b)) = exp(-k/(adj(L)*cosa)) where a      = sat angle, 
!!!                                                                  adj(L) = ray tracing adjustment factor to OD
!!!                k/cos(satzen + f(L)) = k/cos(satzen) 1/adj(L)
!!!                cos(satzen) adj(L) = cos(satzen + f(L))
!!!                adj(L) = cos(satzen + f(L))/cos(satzen) = cos(LayAngle)/cos(satzen)
!!!
!!!                So instead of just sending     absnu(L) = k(L)       to DISORT
!!!                we will send               absnu_adj(L) = k(L)/adj(L)
!!!                and then DISORT will do    absnu_adj(L)/cos(satzen)
!!!

    nabsnu = kMaxPts
    nlev=iNumLayer+1           !this is the number of pressure levels

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

    iiDiv = 0
    555 CONTINUE
    IF (iiDiv*kProfLayer < iaRadLayer(3)) THEN
      iiDiv = iiDiv + 1
      GOTO 555
    END IF
    iiDiv = iiDiv - 1

    DO iFr=1,kMixFilRows
      raVT1(iFr) = raVTemp(iFr)
    END DO

! set the temperature of the bottommost layer correctly
    IF (iDownWard == 1) THEN         !downlook instr
      ! if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot, &
        &                        1,iL)
      rT2 = InterpTempSurf(iProfileLayers,raPressLevels,raVTemp,rFracBot, &
        &                        1,iL,rSurfaceTemp,rSurfPress)
        write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
      ! if the topmost layer is fractional, interpolate!!!!!!
      ! this is hardly going to affect thermal/solar contributions (using this temp
      ! instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp interped to ',raVT1(iL)
    ELSE IF (iDownWard == -1) THEN       !uplook instr
      ! if the bottom layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
      ! if the top layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
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

! ecause kCARTA flips everything for uplook instrument early on in
! he code, we have to reflip everything for DISORT
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
    CALL SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer, &
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
      ELSEIF (iLay == iNumLayer) THEN
        nu=rFracTop
      END IF
      raDensity(iNumLayer-iLay+1) = raNumberDensity(iL-iiDiv*kProfLayer)
      !!!!!Laythick is from outlayers.param    before Nov 17,2001
      !!!!!Laythick is gotten from the KLAYERS profile after Nov 17, 2001
      !!!!!only need for Rayleigh scattering, which is small in infrared
      raThicknessRayleigh(iNumLayer-iLay+1) = &
      raThickness(iL-iiDiv*kProfLayer)*nu
      DO iFr=1,kMaxPts
       ! absprof wants level 1 == TOA, level iNumLayer= gnd
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
    end SUBROUTINE GetAbsProfileDISORT
!************************************************************************

! set the vertical temperatures of the atmosphere
! this sets the temperatures at the pressure level boundaries, using the
! temperatures of the pressure layers that have been supplied by kLayers
! note that temeprature of bottom level is NOT set to rSurfaceTemp
    SUBROUTINE SetDISORTTemp(TEMP,iaRadLayer,raVTemp,iNumLayer, &
    rSurfaceTemp,iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! these are variables that come in from kcartamain.f
    REAL :: raVTemp(kMixFilRows),rSurfaceTemp,raPressLevels(kProfLayer+1)
    INTEGER :: iNumLayer,iProfileLayers
! these are variables that we have to set
    REAL ::    TEMP(*)

! local variables
    INTEGER :: iaRadLayer(kProfLayer), iL, iLay ,iM
    REAL :: Temp1(maxnz)
    REAL :: pavg(kProfLayer),rP,raProfileTemp(kProfLayer)

    DO iLay = 1,MAXNZ
      Temp1(iLay) = 0.0
    END DO

! see which set of Mixed Paths the current atmosphere occupies eg
! set 1 = 1..100, set2= 101..200 etc
! eg if current atmosphere is from MixfilPath 110-190, and kProfLayer = 100,
! then we set iMod as 2      idiv(150,100) = 1  === 2nd set of mixed paths
! assume each atmosphere has at least 25 layers in it!!!
    iM = idiv(iaRadLayer(25),kProfLayer)+1
    DO iLay = 1,kProfLayer
      raProfileTemp(iLay) = raVTemp(iLay+(iM-1)*kProfLayer)
    END DO

    DO iLay = 1,iNumLayer
      iL = iaRadLayer(iLay)
      ! map this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
      iL = iL-idiv(iL,kProfLayer)*kProfLayer
      IF (iL == 0) THEN
        iL = kProfLayer
      END IF
      rP=raPressLevels(iL+1)-10000*delta
      if (rp < raPressLevels(kProfLayer+1)) then
        rp = raPressLevels(kProfLayer+1)+10000*delta
      end if
      TEMP1(iNumLayer-iLay+1) = FindBottomTemp(rP,raProfileTemp,raPressLevels,iProfileLayers)
    END DO

! cc this is where we differ from RTSPEC
! cc      TEMP1(iNumLayer+1) = rSurfaceTemp
    rP = DISORTsurfPress
    TEMP1(iNumLayer+1) = FindBottomTemp(rP,raProfileTemp,raPressLevels,iProfileLayers)

    DO iLay = 1,iNumLayer+1
      temp(iLay) = temp1(iLay)
    END DO

    RETURN
    end SUBROUTINE SetDISORTTemp

!************************************************************************
! this subroutine is from LBLRTMv5.10

    REAL FUNCTION Rayleigh(fr,wtot,layerthick)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
          
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
    REAL :: layerthick        !layer thickness
    REAL :: fr,xrayl,wtot     !fr is the wavenumber
! tot is total number of molecules
! n layer and xrayl is multiplication factor
           
    REAL :: conv_cm2mol,xfr
    REAL :: l,l1,a0,a1,a2,a3,y0,y1,y2,y3

    xrayl = 1.0
    conv_cm2mol = 1./(2.68675e-1*1.e5)

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
    Losch = (kAtm2mb*100/kBoltzmann/273 * 1e-6)    ! (p0/k T0)
    y1 = fr**4/(9.38076e+18 - 1.08426e+09 * fr ** 2) * (wtot/Losch)
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
    end FUNCTION Rayleigh


!************************************************************************
! this subroutine finds out the index of wavenumber points at which to do
! the radiative transfer
    SUBROUTINE FindWavenumberPoints(iStep,iScatter,raCorrelatedK, &
    iaStep,iFFMax)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iStep               !tells how many points to step thru, or find
    INTEGER :: iScatter            !1,2 or 3 tells which type to use
    INTEGER :: iaStep(kMaxPts)     !gives the index of wavenumber points to use
    INTEGER :: iFFMax              !tells how many of the 10000 pts to use
    REAL :: raCorrelatedK(kMaxPts) !array of abs coeffs of layer closest to gnd

    INTEGER :: iF,iL

    DO iF=1,kMaxPts
      iaStep(iF) = 0
    END DO

    IF ((iScatter < 1) .OR. (iScatter > 3)) THEN
      write (kStdErr,*) 'For DISORT, need to choose kScatter = 1,2 or 3'
      write (kStdErr,*) 'iScatter = ',iScatter
      write (kStdErr,*) 'Please retry'
      CALL DoStop
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    IF (iScatter == 1) THEN
      !!!!! use equal spacing of points; nothing special about chosen points
      !!!!! do radiative transfer on freq index pts 1,1+iStep,1+2*iStep,...
      !!!!!       1+N*iStep,10000
      iL=1
      iF=1
  10 CONTINUE
      iaStep(iL) = iF
      iL=iL+1
      iF=iF+iStep
      IF (iF <= kMaxPts) THEN
        GOTO 10
      ELSE
        GOTO 20
      END IF
      20 CONTINUE
      IF (iaStep(iL-1) /= kMaxPts) THEN
        !make sure last point radiances are computed for, is the 10000th pt
        iaStep(iL) = kMaxPts
        iFFMax = iL
      ELSE
        iL=iL-1
        iFFMax = iL
      END IF
    END IF            !IF (kScatter == 1) THEN
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    IF (iScatter == 2) THEN
      !!!!! do radiative transfer on freq index pts a,b,c,...z
      !!!!! where the points are chosen such that they have the smallest
      !!!!! absorption coeffs
      CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
    END IF            !IF (kScatter == 2) THEN

    IF (iScatter == 3) THEN
      !!!!! do radiative transfer on freq index pts a,b,c,...z
      !!!!! where the points are chosen such that they range from the
      !!!!! smallest to largest absorption coeffs
      CALL FindIndices(iScatter,iStep,raCorrelatedK,iaStep,iFFMax)
    END IF            !IF (kScatter == 3) THEN

    RETURN
    end SUBROUTINE FindWavenumberPoints
!************************************************************************
! this subroutine sorts the k values and returns index values
! if iScatter = 2, sort the k values from smallest to largest
!                  return indices of smallest iStep values
! if iScatter = 3, sort the k values from smallest to largest
!                  return indices of smallest to largest k, in steps of iStep
! the outputs are iaStep and iFFMax
    SUBROUTINE FindIndices(iScatter,iStep,arr,iaStep,iFFMax)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iStep               !tells how many points to step thru, or find
    INTEGER :: iScatter            !1,2 or 3 tells which type to use
    INTEGER :: iaStep(kMaxPts)     !gives the index of wavenumber points to use
    INTEGER :: iFFMax              !tells how many of the 10000 pts to use
    REAL :: arr(kMaxPts)           !array of abs coeffs of layer closest to gnd

    INTEGER :: indx(kMaxPts),i,j,IF,IL
    REAL :: TempArr(kMaxPts)

    CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

! now we have sorted the k so that indx(i) contains the sorted keys info
    IF (iScatter == 2) THEN
      ! imply give smallest k values for all but last two points.
      ! or last point, give largest k value
      ! or second last point, give medium k value
      DO i=1,iStep-2
        iaStep(i) = indx(i)
      END DO
      j = kMaxPts-i
      IF (mod(j,2) == 0) THEN
        j=j/2 + i
      ELSE
        j=(j+1)/2 + i
      END IF
      ! iaStep(iStep-1) = indx(kMaxPts-100)
      iaStep(iStep-1) = indx(j)
      iaStep(iStep)   = indx(kMaxPts)
      iFFMax = iStep
    END IF

    IF (iScatter == 3) THEN  !give smallest to largest k values
      iL=1
      iF=1
      10 CONTINUE
      iaStep(iL) = indx(iF)
      iL=iL+1
      iF=iF+iStep
      IF (iF <= kMaxPts) THEN
        GOTO 10
      ELSE
        GOTO 20
      END IF
      20 CONTINUE
      IF (iaStep(iL-1) /= indx(kMaxPts)) THEN
        ! make sure radiances are computed for largest k
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
      30 CONTINUE
      IF (iaStep(iL) == 1)  THEN
        j = +1      !!!!!yes, it does contain "1"
      ELSE
        iL = iL + 1 !!!!keep on hoping
      END IF
      IF ((iL <= iFFMax) .AND. (j < 0)) THEN
        GOTO 30
      END IF
      IF (j < 0) THEN !!!!!!poooey, not found
        i = -1
        iL = 1
        40 CONTINUE
        IF (indx(iL) == 1) THEN
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
      50 CONTINUE
      IF (iaStep(iL) == kMaxPts)  THEN
        j = +1      !!!!!yes, it does contain "kMaxPts"
      ELSE
        iL = iL + 1 !!!!keep on hoping
      END IF
      IF ((iL <= iFFMax) .AND. (j < 0)) THEN
        GOTO 50
      END IF
      IF (j < 0) THEN !!!!!!poooey, not found
        i = -1
        iL = 1
        60 CONTINUE
        IF (indx(iL) == kMaxPts) THEN
          i = +1       !!!!!yes, it is found
          iFFMax = iFFMax + 1
          iaStep(iFFMax) = indx(iL)
        ELSE
          iL = iL + 1
          GOTO 60
        END IF
      END IF

    ! reset iaStep(iFFMax)
      DO i=1,kMaxPts
        TempArr(i) = 1.0e10
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
    end SUBROUTINE FindIndices

!************************************************************************
! this subroutine sorts the k values and returns index values
! if iScatter = 2, sort the k values from smallest to largest
!                  return indices of smallest iStep values
! if iScatter = 3, sort the k values from smallest to largest
!                  return indices of smallest to largest k, in steps of iStep
! the outputs are iaStep and iFFMax
    SUBROUTINE FindIndices_NotUsed(iScatter,iStep,arr,iaStep,iFFMax)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iStep               !tells how many points to step thru, or find
    INTEGER :: iScatter            !1,2 or 3 tells which type to use
    INTEGER :: iaStep(kMaxPts)     !gives the index of wavenumber points to use
    INTEGER :: iFFMax              !tells how many of the 10000 pts to use
    REAL :: arr(kMaxPts)           !array of abs coeffs of layer closest to gnd

    INTEGER :: indx(kMaxPts)
    REAL :: TempArr(kMaxPts),kmin,kmax,k0,dk

    INTEGER :: i,j,iF,iL

    CALL NumericalRecipesIndexer(indx,arr,kMaxPts)

! now we have sorted the k so that indx(i) contains the sorted keys info
!      DO i=1,kMaxPts
!        print *,i,indx(i),arr(indx(i))
!      END DO

    IF (iScatter == 2) THEN
      ! imply give smallest k values for all but last two points.
      ! or last point, give largest k value
      ! or second last point, give medium k value
      DO i=1,iStep-2
        iaStep(i) = indx(i)
      END DO
      j = kMaxPts-i
      IF (mod(j,2) == 0) THEN
        j=j/2 + i
      ELSE
        j=(j+1)/2 + i
      END IF
      !iaStep(iStep-1) = indx(kMaxPts-100)
      iaStep(iStep-1) = indx(j)
      iaStep(iStep)   = indx(kMaxPts)
      iFFMax = iStep
    END IF

    IF (iScatter == 3) THEN
      ! smallest to largest k values, equally stepped (kmax-kmin)/(iStep-1)
      kmin = arr(indx(1))
      kmax = arr(indx(kMaxPts))
      dk = (kmax-kmin)/(iStep-1)
      k0 = kmin

      k0 = k0 + dk
      iaStep(1) = indx(1)
      iF=2
      iL=2

      10 CONTINUE
      IF ((arr(indx(iL)) < k0) .AND. (iL < kMaxPts)) THEN
        iL = iL + 1
        GOTO 10
      ELSE
        iaStep(iF) = indx(iL)
        iF = iF + 1
        iL = iL + 1
        k0 = k0 + dk
        IF (k0 < kmax) THEN
          GOTO 10
        END IF
      END IF
      IF (iaStep(iF-1) /= indx(kMaxPts)) THEN
        ! make sure radiances are computed for largest k
        iaStep(iF) = indx(kMaxPts)
        iFFMax = iF
      ELSE
        iF=iF-1
        iFFMax = iF
      END IF
    END IF

! very important .. want to keep integrity of the Planck fcn!!!!!
! so make sure the original 1st and 10000th points are always used
    IF (iScatter >= 2) THEN
      !!!!!!make sure (iaStep(i) contains "1")
      j = -1
      iL = 1
      30 CONTINUE
      IF (iaStep(iL) == 1)  THEN
        j = +1      !!!!!yes, it does contain "1"
      ELSE
        iL = iL + 1 !!!!keep on hoping
      END IF
      IF (iL <= iFFMax) THEN
        GOTO 30
      END IF
      IF (j < 0) THEN !!!!!!poooey, not found
        i = -1
        iL = 1
        40 CONTINUE
        IF (indx(iL) == 1) THEN
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
      50 CONTINUE
      IF (iaStep(iL) == kMaxPts)  THEN
        j = +1      !!!!!yes, it does contain "kMaxPts"
      ELSE
        iL = iL + 1 !!!!keep on hoping
      END IF
      IF (iL <= iFFMax) THEN
        GOTO 50
      END IF
      IF (j < 0) THEN !!!!!!poooey, not found
        i = -1
        iL = 1
        60 CONTINUE
        IF (indx(iL) == kMaxPts) THEN
          i = +1       !!!!!yes, it is found
          iFFMax = iFFMax + 1
          iaStep(iFFMax) = indx(iL)
        ELSE
          iL = iL + 1
          GOTO 60
        END IF
      END IF

      ! reset iaStep(iFFMax)
      DO i=1,kMaxPts
        TempArr(i) = 1.0e10
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
    end SUBROUTINE FindIndices_NotUsed
!************************************************************************
    SUBROUTINE INTERP_PLANCK_0(XA,YA,CA,N,raFO,raNS,raOut)

! this does a very smart interpolation ahem ahem
!        CALL INTERP_PLANCK_0(
!     $          raFreqStep,raIntenStep,raNoScatterStep,iStepPts,
!     $          raFreq,raNoScatterAll,raInten)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      XA  : I  : REAL arr : freq array   array(N) in increasing order
!      YA  : I  : REAL arr : intensity y  array(N) from scattering
!      CA  : I  : REAL arr : intensity y  array(N) from clear sky
!      N   : I  : INT      : number of points fom DISORT, in array

!      raFO  : I  : REAL arr : entire wavenumber           array (kMaxPts)
!      raNS  : I  : REAL arr : entire non scatter radiance array (kMaxPts)
!      raOut : O  : REAL arr : output radiance             array (kMaxPts)

    REAL :: XA(*),YA(*),CA(*)
    REAL :: raFO(*),raNS(*),raOut(*)
    INTEGER :: N

!      Local Variables
    INTEGER :: K,KLO,KHI,iI,iMethod,iDefault
    REAL :: A,B,H
    REAL :: yn,y0,X,y,w,dx

    iMethod = 2        !! uses rads; seems to be a good method, as it
!! will not barf if the radiances are 0.0
!! but it does lead to some noticeable problems
    iMethod = 3        !! usesBT(rads); seems to be a good method. needs some
!! tweaks to handle TOA rads for uplook inst (==0)

    DO iI=1,kMaxPts
      X = raFO(iI)
      IF (X > xa(n)) X = xa(n) - 1.0e-7
      IF (X < xa(1)) X = xa(1) + 1.0e-7

      ! determine between which pair of points X falls (bisect loop)
      KLO = 1
      KHI = N

  20  IF ( (KHI - KLO) > 1) THEN
        K = (KHI + KLO)/2
        IF (XA(K) > X) THEN
          KHI = k
        ELSE
          KLO = k
        ENDIF
      GOTO 20
      ENDIF
    
      dx = XA(KHI) - x
      H  = XA(KHI) - XA(KLO)

      IF (H <= 0.0) THEN
        WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI),X,N
  1010 FORMAT('ERROR! linear SPLINT: bad XA array.',/, &
            'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5, &
            ', XA(KHI) = ',E12.5,',X = ',E12.5,', numpts = ',I4,'. Quitting.')
      ENDIF

      IF (iMethod == 1) THEN
        ! ix this, as we need to define dh,dl
        !B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
        !Y=(YA(KHI) + dh) - h*b
        write (kStdErr,*) 'iMethod = 1 invalid in INTERP_PLANCK_0'
        CALL DoStop

      ELSEIF (iMethod == 2) THEN
        ! nterpolate radiances
        ! A  : I  : REAL arr : intensity y array(N) from scattering
        ! A  : I  : REAL arr : intensity y array(N) from clear sky
        !!!interpolate in radiance space
        yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
        y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
        B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
        Y = yn - dx*b
        Y = Y + raNS(iI)

      ELSEIF (iMethod == 3) THEN
        IF (ya(KHI) > 0 .AND. CA(KHI) > 0 .AND.  ya(KLO) > 0 &
             .AND. CA(KLO) > 0) THEN
          ! METHOD 3
          ! interpolate in BT space if rads ~= 0
          !!!interpolate in temperature space
          !!! get out stuff from clear sky
          yn = radtot(xa(khi),ya(KHI))-radtot(xa(khi),CA(KHI))
          y0 = radtot(xa(klo),ya(Klo))-radtot(xa(klo),CA(Klo))
          B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
          Y = yn - dx*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
          Y = Y + radtot(x,raNS(iI)) !!add on clear sky stuff
          Y = ttorad(x,y)
        ELSE
          ! METHOD 3
          ! interpolate in rad space if rads = 0
          yn = ya(KHI)-CA(KHI)      !!! get out stuff from clear sky
          y0 = ya(KLO)-CA(KLO)      !!! get out stuff from clear sky
          B = (yn-y0)/(XA(KHI) - XA(KLO)) !!slope
          Y = yn - dx*b
          Y = Y + raNS(iI)
        END IF

      ELSEIF (iMethod == 4) THEN
        ! METHOD 4 : directly interpolate in BT space
        ! interpolate BTs
        ! fix this
        write (kStdErr,*) 'iMethod = 4 invalid in INTERP_PLANCK_0'
        !B = (radtot(WA(KHI),YA(KHI))-radtot(WA(KLO),YA(KLO)))/
        !     $        (XA(KHI) - XA(KLO))
        !Y=radtot(WA(KHI),YA(KHI))-dx*b !!r(vo+dv,ko+dk)=r(vo,ko)+dr/dko dk
        !Y=ttorad(W,Y)+0*dh             !!change bck to rad, add on dr/dvodv

      ELSEIF (iMethod == 5) THEN
        write (kStdErr,*) 'iMethod = 5 invalid in INTERP_PLANCK_0'
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

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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
    REAL :: XA(*),YA(*),WA(*),DWA(*),X,Y,W
    INTEGER :: N

!      Local Variables
    INTEGER :: K,KLO,KHI
    REAL :: A,B,H,dh,dl

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

 20 IF ( (KHI - KLO) > 1) THEN
      K=(KHI + KLO)/2
      IF (XA(K) > X) THEN
        KHI = k
      ELSE
        KLO = k
      ENDIF
      GOTO 20
    ENDIF

    H = XA(KHI) - XA(KLO)
    dh = W - WA(khi)       !!!wavenumber diff
    dl = WA(klo) - W       !!!wavenumber diff
    dh = dh * DWA(khi)     !!!first derivative
    dl = dl * DWA(klo)     !!!first derivative

    IF (H <= 0.0) THEN
      WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
1010 FORMAT('ERROR! linear SPLINT: bad XA array.',/, &
        'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5, &
        ', XA(KHI) = ',E12.5,'. Quitting.')
    ENDIF

! METHOD 1
!      B=((YA(KHI) + dh) - (YA(KLO) + dl))/(XA(KHI) - XA(KLO)) !!slope
!      Y=(YA(KHI) + dh) - h*b

! METHOD 2
! interpolate radiances
!      B=(YA(KHI) - YA(KLO))/(XA(KHI) - XA(KLO)) !!slope
!      Y=YA(KHI) - h*b   !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
!      Y=Y+dh            !!add on                            + dr/dvo dv
! interpolate BTs
    B = (radtot(WA(KHI),YA(KHI)) - radtot(WA(KLO),YA(KLO)))/(XA(KHI) - XA(KLO))
    Y = radtot(WA(KHI),YA(KHI)) - h*b   !!r(vo+dv,ko+dk) = r(vo,ko) + dr/dko dk
    Y = ttorad(W,Y)+0*dh                !!change back to rad, add on  dr/dvo dv

! METHOD 3
!      Y=YA(KHI) + dh    !!this is r(vo+dv,ko+dk) = r(vo,ko) + dr/dvo dv

    RETURN
    END SUBROUTINE INTERP_PLANCK_1

!************************************************************************
    SUBROUTINE INTERP_PLANCK_2(WA,DWA,XA,YA,CA,N, &
    raFO,raNS,raCK,raOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    REAL :: XA(*),YA(*),CA(*),WA(*),DWA(*)
    REAL :: raFO(*),raNS(*),raCK(*),raOut(*)
    INTEGER :: N

!      Local Variables
    INTEGER :: K,KLO,KHI,iI
    REAL :: A,B,H,dh,dl
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

    ! determine between which pair of points X falls (bisect loop)
      KLO=1
      KHI=N

 20   IF ( (KHI - KLO) > 1) THEN
        K=(KHI + KLO)/2
        IF (XA(K) > X) THEN
          KHI = k
        ELSE
          KLO = k
        ENDIF
        GOTO 20
      ENDIF
      
      dx = XA(KHI) - x
      H  = XA(KHI) - XA(KLO)
      dh = W - WA(khi)       !!!wavenumber diff
      dl = WA(klo) - W       !!!wavenumber diff
      dh = dh * DWA(khi)     !!!first derivative
      dl = dl * DWA(klo)     !!!first derivative
  
      IF (H <= 0.0) THEN
        WRITE(kStdWarn,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010   FORMAT('ERROR! linear SPLINT: bad XA array.',/, &
            'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5, &
            ', XA(KHI) = ',E12.5,'. Quitting.')
      ENDIF

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

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    REAL :: XA(*),YA(*),CA(*),WA(*)
    REAL :: raFO(*),raNS(*),raCK(*),raOut(*)
    INTEGER :: N

!     Local Variables
    INTEGER :: K,KLO,KHI,iI,indx(kMaxPts),iDok,iDoF
    REAL :: A,B,H
    REAL :: yn,y0,X,Y1,Y2,YTotal,dx
! these have WA,XA,YA,CA sorted according to wavenumber
    REAL :: WA_sort(kMaxPts),XA_sort(kMaxPts)
    REAL :: YA_sort(kMaxPts),CA_sort(kMaxPts)

    1010 FORMAT('ERROR! linear SPLINT: bad XA array.',/, &
         'KLO=',I5,', KHI=',I5,', XA(KLO) = ',1PE12.5, &
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
           
        !determine between which pair of points X falls (bisect loop)
        KLO=1
        KHI=N
            
     20 IF ( (KHI - KLO) > 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) > X) THEN
            KHI = k
          ELSE
            KLO = k
          ENDIF
          GOTO 20
        ENDIF

        dx = XA(KHI) - x
        H  = XA(KHI) - XA(KLO)

        IF (H <= 0.0) THEN
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
      IF (iDoF > 0) THEN
        !!!!! -------------  then do dT/dv dv !!!!!!!!!!!!!
        X=raFO(iI)
              
        ! etermine between which pair of points X falls (bisect loop)
        KLO=1
        KHI=N

    30  IF ( (KHI - KLO) > 1) THEN
          K=(KHI + KLO)/2
          IF (WA_Sort(K) > X) THEN
            KHI = k
          ELSE
            KLO = k
          ENDIF
          GOTO 30
      ENDIF
  
      dx = WA_Sort(KHI) - x
      H  = WA_Sort(KHI) - WA_Sort(KLO)

      IF (H <= 0.0) THEN
        WRITE(kStdWarn,1010) KLO,KHI,WA_Sort(KLO),WA_Sort(KHI)
      ENDIF

      !!!interpolate in temperature space
      !!! get out stuff from clear sky
      yn = radtot(wa_sort(khi),ya_sort(KHI)) - &
      radtot(wa_sort(khi),ca_sort(KHI))
      y0 = radtot(wa_sort(klo),ya_sort(Klo)) - &
      radtot(wa_sort(klo),ca_sort(Klo))
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

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! linear interpolation
!      -----------------------------------------------------------------
!      WA  : I  : REAL arr : wavenumber w array(N)
!     DWA  : O  : REAL arr : dy/dw so we can weight the interpolations
!      YA  : I  : REAL arr : intensity y array(N)
!      N   : I  : INT      : number of points in arrays
!      -----------------------------------------------------------------

!     Parameters
    REAL :: YA(*),WA(*),DWA(*)
    INTEGER :: N
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
      raT(i) = r2*wa(i)/(log(1+r1*(wa(i)**3)/ya(i)))
    END DO

! now do dr/df where r=ttorad = planck radiance at temp T, frequency fr
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
    DO i = 1,N
      u  = r1*(wa(i)**3)
      du = r1*3*(wa(i)**2)
      v  = exp(r2*wa(i)/raT(i))-1
      dv = r2/raT(i) * exp(r2*wa(i)/raT(i))
      !!!! now do the derivative
      dwa(i) = (v*du - u*dv)/(v*v)
    END DO

    RETURN
    end SUBROUTINE drad_dv

!************************************************************************
! this doe a cumulative K distribution fcn, and then for the first 50%
! of the lowest k values, just does an average of intensity
    SUBROUTINE CumulativeK(raKall,raF,raK,raI,iN)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: raKAll(kMaxPts)        !this is all 10000 pts k(lowest layer)
    REAL :: raF(kMaxPts)           !these are freqs of pts where DISORT used
    REAL :: raK(kMaxPts)           !these are ks of pts where DISORT used
    REAL :: raI(kMaxPts)           !intensities of pts where DISORT used
    INTEGER :: iN                  !how many points DISORT went over

! note that we will do a sort on k, and then accordingly, average some of the
! points in raK and raI so as to make things smoother

    INTEGER :: iCnt50,iCnt75,iInCnt50,iInCnt75
    REAL :: rCnt50,rCnt75,rF,rI,rK
    INTEGER :: indx(kMaxPts),iI,iJ,iCnt
    REAL :: NormalizedK(kMaxPts),kmin,kmax,dk
    REAL :: TempF(kMaxPts),TempI(kMaxPts),TempK(kMaxPts)
    REAL :: raCumulativeK(kMaxPts),raCumulativeDistr(kMaxPts)

    CALL NumericalRecipesIndexer(indx,raKAll,kMaxPts)
          
! set up the normalised vector knorm = ksorted/max(k)
    DO iI = 1,kMaxPts
        NormalizedK(iI) = raKall(indx(iI))/raKall(indx(kMaxPts))
    END DO

    kmin = 0.0
    kmax = NormalizedK(kMaxPts)   !!!should be 1.0

    iJ = 500
    dk = 1.0/iJ

! compute the cumulative distr function
    iCnt = 1
    DO iI = 1,iJ
      raCumulativeK(iI) = kmin + iI*dk
      10 CONTINUE
      IF (NormalizedK(iCnt) <= raCumulativeK(iI)) THEN
        iCnt = iCnt + 1
        IF (iCnt <= kMaxPts) THEN
          GOTO 10
        END IF
      END IF
      iCnt = min(iCnt,kMaxPts)
      raCumulativeDistr(iI) = iCnt*1.0/kMaxPts
    END DO

! now see which index required before 50% of the cumulative function is used
    iCnt50 = 1
    iCnt75 = 1
    30 CONTINUE
    IF (raCumulativeDistr(iCnt50) < 0.5) THEN
      iCnt50 = iCnt50 + 1
      GOTO 30
    END IF
    40 CONTINUE
    IF (raCumulativeDistr(iCnt75) < 0.75) THEN
      iCnt75 = iCnt75 + 1
      GOTO 40
    END IF

! ind values of raK corresponding to the normalised k=0.5,0.75 values
! emember raK is where DISORT is called, and so these should be
! lready sorted from smallest to largest
    rCnt50 = raCumulativeK(iCnt50)*raKall(indx(kMaxPts))
    rCnt75 = raCumulativeK(iCnt75)*raKall(indx(kMaxPts))
    iInCnt50 = 1
    iInCnt75 = 1
    50 CONTINUE
    IF (raK(iInCnt50) < rCnt50) THEN
      iInCnt50 = iInCnt50 + 1
      GOTO 50
    END IF
    60 CONTINUE
    IF (raK(iInCnt75) < rCnt75) THEN
      iInCnt75 = iInCnt75 + 1
      GOTO 60
    END IF
!      print *,'normalised cml upto 50% = ',iCnt50,raCumulativeK(iCnt50),rCnt50
!      print *,'normalised cml upto 75% = ',iCnt75,raCumulativeK(iCnt75),rCnt75
!      print *,'input cml upto 50% = ',iInCnt50,raK(iInCnt50)
!      print *,'input cml upto 75% = ',iInCnt75,raK(iInCnt75)

! ince there is so much variation in the lowest radainces, and they
! ccount for so much of raK, let us avg the first few values and
! hove the rest around
    IF (raCumulativeK(iCnt50) <= 0.1) THEN
      ! irst initialise the temp arrays
      DO iI = 1,iN
        TempF(iI) = raF(iI)
        TempK(iI) = raK(iI)
        TempI(iI) = raI(iI)
      END DO
    ! 0% of the (k/kmax) are smaller in magnitude that 0.1
    ! this will drive the interpolation wrt raKStep nuts!!!!
    ! so do an average;
    ! make sure at least FIVE points remain for the interpolations
      IF ((iN - iInCnt50) <= 5) THEN
        iInCnt50 = iN-5
      END IF
      rF = 0.0
      rI = 0.0
      rK = 0.0
      DO iI = 2,iInCnt50
        rF = rF+raF(iI)
        rK = rK+raK(iI)
        rI = rI+radtot(raF(iI),raI(iI))      !!!average the temps!!!
      END DO
      rF = rF/(iInCnt50-1)
      rK = rK/(iInCnt50-1)
      rI = rI/(iInCnt50-1)
      rI = ttorad(rF,rI)                      !!!change to radiance!!!
                
      ! now shove this info into the arrays
      ! keep the smallest "k" info!!!!!
      raF(1) = raF(1)
      raK(1) = raK(1)
      raI(1) = raI(1)
      ! you have averaged the next few to get something
      raF(2) = rF
      raK(2) = rK
      raI(2) = rI
      ! now fill in the larger "k" values
      DO iI=3,iN - iInCnt50 + 2
        raF(iI) = TempF(iInCnt50 + iI - 2)
        rak(iI) = Tempk(iInCnt50 + iI - 2)
        raI(iI) = TempI(iInCnt50 + iI - 2)
      END DO
      IN = iN - iInCnt50 + 2
    END IF

    RETURN
    end SUBROUTINE CumulativeK
!************************************************************************

!! see eg /home/sergio/KCARTA/SCATTERCODE/DISORT4.0.99/disort4.0.99/doc/DISORT.txt
      SUBROUTINE InputPrintDebugDisort(raFreq, NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
           USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, &
           IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, BTEMP, TTEMP, TEMIS, &           
           PLANK, WVNMLO, WVNMHI, ONLYFL, ACCUR, PRNT, HEADER, &
           ICLDTOP, ICLDBOT, iaRadLayer,iNumLayer,raPressLevels,rSurfPress)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: ICLDTOP, ICLDBOT, iaRadLayer(kProfLayer),iNumLayer

    REAL :: raFreq(kMaxPts),raPressLevels(kProflayer+1),rSurfPress

    CHARACTER  header*127          !dumb comment
          
    LOGICAL :: lamber  !true  ==> isotropic reflecting lower bdry
!          so need to specify albedo
! alse ==> bidirectionally reflecting bottom bdry
!      ==> need to specify fcn BDREF()
    LOGICAL :: plank   !use the plank function for local emission
    LOGICAL :: onlyfl  !true ==> return flux, flux divergence, mean intensitys
! alsetrue ==> return flux, flux divergence, mean
!              intensities and intensities
    LOGICAL :: prnt(5) !prnt(1) = true, print input variables
! rnt(2) = true, print fluxes
! rnt(3) = true, intensities at user angles and levels
! rnt(4) = true, print planar transmitivity and albedo
! rnt(5) = true, print phase function moments PMOM
    LOGICAL :: usrtau  !false ==> rad quantities return at every bdry
! rue  ==> rad quantities return at NTAU optical depths
!          which will be specified in utau(1:ntau)
    LOGICAL :: usrang  !false ==> rad quantities returned at stream angles
! rue  ==> rad quantities returned at user angles
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
! o nlvl = nlyr + 1
    INTEGER :: nstr             !number of radiation streams
    INTEGER :: ntau             !associated with LOGICAL usrtau, print results
! t this many optical depths
    INTEGER :: numu             !associated with LOGICAL usrang, specifies how
! any polar angles results to be printed (umu)
    INTEGER :: nphi             !specifies how many azimuth angles results to
! e printed (phi) .. can only be 0 if
! nlyfl = .true.

    DOUBLE PRECISION :: accur   !accuracy convergence criterion for azimuth
! fourier cosine) series .. set between 0-0.01
    DOUBLE PRECISION :: albedo  !bottom bdry albedo
    DOUBLE PRECISION :: btemp   !bottom surface temp
    DOUBLE PRECISION :: dtauc(maxcly)
! ptical depths of computational layers
    DOUBLE PRECISION :: fisot   !intensity of top bdry isotropic illumination
    DOUBLE PRECISION :: fbeam   !intensity of incident // beam at TOA
    DOUBLE PRECISION :: phi(maxphi)
! he azimuthal phi's to output radiance
    DOUBLE PRECISION :: pmom(0:maxmom,maxcly)
! cattering phase fcn in legendre polynoms
    DOUBLE PRECISION :: phi0    !solar azimuth
    DOUBLE PRECISION :: ssalb(maxcly)
! ingle scatter albedo of computational layers
    DOUBLE PRECISION :: temper(0:maxcly) !temperature of the levels (0=TOA)
    DOUBLE PRECISION :: temis            !emissivity of top bdry
    DOUBLE PRECISION :: ttemp            !temperature of top bdry
    DOUBLE PRECISION :: wvnmhi, wvnmlo   !bounds within which to do computations
    DOUBLE PRECISION :: umu(maxumu)      !ang's at which to output results
    DOUBLE PRECISION :: umu0         !polar angle cosine of incident solar beam
    DOUBLE PRECISION :: utau(maxulv)     !tau's at which to output results

    INTEGER :: iI,iJ,iK

    write(kStdWarn,*) ' '
    write(kStdWarn,*)           ' DISORT input for raFreq(1) = ',raFreq(1)
    WRITE(kStdWarn,'(A,F12.5)') ' emissivity = ',TEMIS
    write(kStdWarn,'(A,I4)') ' NLYR          Number of computational layers   = ',NLYR
    write(kStdWarn,'(A,I4)') ' NMOM          Number of phase function moments = ',NMOM
    write(kStdWarn,'(A)')    '   IK  IJ  iI      PLEV        TLEV     DTAUC      SSALB      PMOM(0)      PMOM(1)     PMOM(2)      PMOM(3)'
    write(kStdWarn,'(A)')    ' ------------------------------------------------------------------------------------------------------------'
    WRITE(kStdWarn,'(3(I4),2(F12.5))') 0,0,0,0.0,TTEMP
    iI = 1
    iJ = iaRadLayer(iI)
    iK = kProfLayer - iJ + 1
    iI = iI - 1
    iJ = iJ - 1
    iK = iK + 1
    write(kStdWarn,'(I4,I4,I4,2(F12.5))') iK,iJ,iI,raPressLevels(iK+(kProfLayer-iNumLayer)),TEMPER(iI)
    DO iI = 1,NLYR
      IJ = iaRadLayer(iI)
      iK = kProfLayer - iJ + 1
      IF ((iI .EQ. ICLDTOP) .OR. (iI .EQ. ICLDBOT)) THEN
        write(kStdWarn,'(A)')    ' ------------------------------------------------------------------------------------------------------------'
      END IF
      WRITE(kStdWarn,'(3(I4),2(F12.5),6(ES12.5))') & 
       iK,iJ,iI,raPressLevels(iK+(kProfLayer-iNumLayer)),TEMPER(iI),DTAUC(iI),SSALB(iI),PMOM(0,iI),PMOM(1,iI),PMOM(2,iI),PMOM(3,iI)
    END DO
    WRITE(kStdWarn,'(3(I4),2(F12.5))') NLYR+1,NLYR+1,NLYR+1,rSurfPress,BTEMP 
    write(kStdWarn,'(A)')    '   IK  IJ  iI      PLEV        TLEV     DTAUC      SSALB      PMOM(0)      PMOM(1)     PMOM(2)      PMOM(3)'
    write(kStdWarn,'(A)')    ' ------------------------------------------------------------------------------------------------------------'

    IF (LAMBER .EQ. .true.) THEN
      write(kStdWarn,'(A)') 'LAMBER = true = isotropic reflecting boundary'
    ELSE
      write(kStdWarn,'(A)') 'LAMBER = false = bidirectional reflecting boundary'
    END IF
    IF (PLANK .EQ. .true.) THEN
      write(kStdWarn,'(A)') 'PLANK = true = thermal emission'
    ELSE
      write(kStdWarn,'(A)') 'PLANK = false = no thermal emission !!!!!!!!!!!!!!!!!!'
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,I3)') 'NSTR = Number of computational polar angles to be used = streams',NSTR
    
    write(kStdWarn,*) ' '
    IF (USRANG .EQ. .false.) THEN
       write(kStdWarn,'(A)') 'USRANG = false, radiant quantities at gauss quad streams'
    ELSE
      write(kStdWarn,'(A,I3,A)') 'USRANG = true, return rads at ',NUMU,' specified cos(polar angle),polar angle : '
      DO iI = 1,NUMU
        WRITE(kStdWarn,'(I3,2(F12.5))') iI,UMU(iI),ACOS(UMU(iI))*180.0/kPi
      END DO       
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,I3,A)') ' return rads at ',NPHI,' specified azimuth angs (0 ok only if ONLYFL = true)'
    DO iI = 1,NPHI
      WRITE(kStdWarn,'(I3,F12.5)') iI,PHI(iI)
    END DO       

    write(kStdWarn,*) ' '
    IF (USRTAU .EQ. .false.) THEN
       write(kStdWarn,'(A)') 'USRTAU = false, return rads at every layer'
    ELSE
      write(kStdWarn,'(A,I3,A)') 'USRTAU = true, return rads at ',NTAU,' specified layers/ODs : '
      DO iI = 1,NTAU
        WRITE(kStdWarn,'(I3,ES12.5)') iI,UTAU(iI)
      END DO       
    END IF

    write(kStdWarn,*) ' '
    IF (IBCND .EQ. 0) THEN
      write(kStdWarn,'(A)') 'IBCND = 0    GENERAL CASE'
    ELSEIF (IBCND .EQ. 1) THEN
      write(kStdWarn,'(A)') 'IBCND = 1    Return alb/transmissivity (no Planck sources) = Re/Tr = PCRTM at UMU output angles'
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,3(ES12.5))') 'FBEAM intensity,cos(sol0),sol0,phi0',FBEAM,UMU0,acos(UMU0)*180/kPi,PHI0
    write(kStdWarn,'(A,ES12.5)')    'FISOT intensity of top boundary illumination',FISOT

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,2(F12.5),A,ES12.5)')  'WAVENO ',WVNMLO,WVNMHI,' ACCUR',ACCUR
    write(kStdWarn,'(A,L)') 'ONLYFL ',ONLYFL
    write(kStdWarn,'(A,L)') 'PRNT 1 ',PRNT(1)
    write(kStdWarn,'(A,L)') 'PRNT 2 ',PRNT(2)
    write(kStdWarn,'(A,L)') 'PRNT 3 ',PRNT(3)
    write(kStdWarn,'(A,L)') 'PRNT 4 ',PRNT(4)
    write(kStdWarn,'(A,L)') 'PRNT 5 ',PRNT(5)

    write(kStdWarn,*) ' '

    RETURN
    end SUBROUTINE InputPrintDebugDisort

!************************************************************************

END MODULE scatter_disort_aux
