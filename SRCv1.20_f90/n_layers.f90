! Copyright 2014
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:42
 
! University of Maryland Baltimore County
! All Rights Reserved

! this file deals with a simple layers
! does everything in terms of LINEAR interps, not SPLINE
! DO NOT USE        Call r_sort_logspl(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
! DO     USE        Call r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)

! assumes WV is in "wet air", other gases in dry air VMR (klayers unit 12)
! code
!  (a) reads in input text levels profile, gases in whatever units
!      changes the input units to MIXING RATIOS
!  (b) read in upper level VMR profiles for above gases, as well as for "missing gases"
!      adjusts the Ref VMR to agree with input where the input levels end, and then tapers the
!        adjust ratio to 1, with about 4 AIRS levels
!      adjusts the Ref Temp to agree with input where the input levels end, and then tapers the
!        adjust ratio to 0, with about 4 AIRS levels
!   ** the tacking on (or adding new gases) can either be via
!       - using the reference P,PP,T for 100 layers ==> VMR = PP/P
!       - using one of the six AFGL profiles        ==> VMR comes from the ~50 levels
!  (c) if Earth do the adjustment of dry air mix ratios
!  (d) does the sublev integrations

! Useful References
!   http://www.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_1/index.html
!   http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! saved off in ../PDF/wisc_notes_*.pdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! this subroutine reads in a TEST file LEVELS mixing ratio PROFILE
! the text file will have N levs x (G+2) columns of the format

! caStr80   = comment
! numlevs   = number of levels
! Psurf,TSurf,HSurf = surface pressure (mb) and temperature (K) and height (km)
! year lat lon     = year of obs (for co2 MR) and gravity adjust (not yet implemented)
!                     if year < 0, ignore
!                     if year > 0 and kPlanet = 3, then if adding on Gas2, adjust for ppmv
! numgases  = number of gases
! gasIDS    = which gases
! gasUNITS  = eg MR, g/g, RH, same as klayers
!   p1     T1   g1_1  g2_1 .... gG_1
!   p2     T2   g1_2  g2_2 .... gG_2
!   p3     T3   g1_3  g2_3 .... gG_3
!   p4     T4   g1_4  g2_4 .... gG_4

!   pN     TN   g1_N  g2_N .... gG_N

! where p(i) is in mb, T(i) is in Kelvin and N <= kProfLayer*2
! gi_j (j=1--N) are the gas amounts eg dimensionless 0 < VMR < 1, 0 < RH < 100 etc

!        See eg KCARTA/IP_PROFILES/levelsRTP_to_levelstext.m
!             KCARTA/IP_PROFILES/levelsprofile_text1.prf

! see /asl/packages/klayers/Doc/description.txt

! or can read TAPE5, modified TAPE6 from LBLRTM


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! for testing just do SUBROUTINE InputMR_profile(caPfName)    !! for testing

SUBROUTINE InputMR_profile(caPfName,iNumLays,iNumGases,iaG,iLowestLev,  &
    raTout,raAmountOut,raZout, raPout,raaQout,raaPartPressout,  &
    raPbndFinal,raTbndFinal,iZbndFinal)


CHARACTER (LEN=80), INTENT(IN OUT)       :: caPfName
NO TYPE, INTENT(OUT)                     :: iNumLays
NO TYPE, INTENT(IN OUT)                  :: iNumGases
NO TYPE, INTENT(IN OUT)                  :: iaG
NO TYPE, INTENT(OUT)                     :: iLowestLev
REAL, INTENT(IN OUT)                     :: raTout(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raAmountOu
REAL, INTENT(OUT)                        :: raZout(kProfLayer+1)
NO TYPE, INTENT(OUT)                     :: raPout
NO TYPE, INTENT(IN OUT)                  :: raaQout
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: raPbndFina
NO TYPE, INTENT(IN OUT)                  :: raTbndFina
NO TYPE, INTENT(OUT)                     :: iZbndFinal
IMPLICIT NONE

INTEGER :: iplev
INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
INCLUDE '../INCLUDE/airslevelheightsparam.f90'

! input

! output
REAL :: !! in K
REAL :: raAmountOut(kProfLayer)             !! in molecules/cm2
REAL :: !! in meters, notice this is at LEVELS

REAL :: raPout(kProfLayer)                  !! in N/m2 (even though input TAPE5,TAPE6, kRTP=-10, pressures are in mb)
REAL :: raaQout(kProfLayer,kMaxGas)         !! in molecules/cm2
REAL :: raaPartPressout(kProfLayer,kMaxGas) !! in N/m2
INTEGER :: iNumLays,iNumGases,iaG(kMaxGas)
INTEGER :: iLowestLev
REAL :: raPbndFinal(kProfLayer+1)           !! final output pressure boundaries, default = 101 AIRS levels
REAL :: raTbndFinal(kProfLayer+1)           !! final output temperatures at boundaries
INTEGER :: iZbndFinal                       !! number of output plev bdries,
!! should be iNumLays+1 if using LBLRTM bdies, else 101

! local
INTEGER :: iG,iL,iJ,iK,iaGasUnits(kMaxGas)
INTEGER :: iFound,iHighestLay,iNumLevs
REAL :: raP(2*kProfLayer),raT(2*kProfLayer),raaG_VMR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
REAL :: raX(kMaxGas),rX,rP,rT,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
REAL :: rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
INTEGER :: iKCARTADirectionUPorDOWN,iMid,iFloor,iCnt,iOffSet
REAL :: gamma,z,dz,slope,amount,junk
REAL :: rFracBot
REAL :: raPBnd(2*kProfLayer)  !! do we want user defined pressure level boundaries??
INTEGER :: iZbnd              !! do we want user defined pressure level boundaries??

REAL :: raPX(kProfLayer+1),raZX(kProfLayer+1),raTX(kProfLayer+1)
REAL :: raaG_VMRX(kProfLayer+1,kMaxGas),raLayDensityX(kProfLayer+1)
REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)

CALL Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,  &
    rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta)

DO iL = 1,kMaxLayer+1
  raPbndFinal(iL) = PLEV_KCARTADATABASE_AIRS(iL)
  raTbndFinal(iL) = 0.0
END DO
iZbndFinal = kMaxLayer+1

iZbnd = +1   !! default to using the input levels as boundaries when doing levels --> layers integration
iZbnd = -1   !! default to using AIRS 101 levels  as boundaries when doing levels --> layers integration

! >>> read in user supplied prof
! this will be at whatever gas units eg RH g/kg  VMR etc, but then CONVERTED to VMR (gasunit 12)
IF (kRTP == -10) THEN
  CALL ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
      iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,  &
      raP,raT,iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)
  iZbnd = -1
ELSE IF (kRTP == -5) THEN
  CALL ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
      iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases, raP,raT,raAlt,iZbnd,raPBnd,  &
      iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)
ELSE IF (kRTP == -6) THEN
!! edited TAPE6, which contains first few parts of TAPE5, and then just profile info from TAPE6 (mol/cm2)
  CALL ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
      iNumLays,rPSurf,rTSurf,rHSurf,iNumGases, raPX,raTX,raLayDensityX,raZX,  &
      raPressLevelsX,raTPressLevelsX,raAltitudesX,  &
      iaG,iaGasUnits,raaG_VMRX,rPMin,rPMax,rYear,rLat,rLon)
  iZbnd = -1
ELSE
  WRITE(kStdErr,*) 'huh?? reading text levels file or text LBLRTM TAPE5/TAPE6 MODIFIED',kRTP
  CALL DoStop
END IF

IF (kRTP /= -6) THEN
! >>> interp onto the kCARTA database profiles
! >>> also sort the gasIDs into increasing order
  
! find pseudo temperature at this level
! recall pav = (p2-p1)/log(p2/p1)
!        tav = (T2-T1)/log(p2/p1)
  
!! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
!! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
!! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
  
  IF (iZbnd < 0) THEN
! >>>>>>>>>> default AIRS 101 levels
    
! >>>> tack on Standard Profile to User Profile if needed;
! >>>> also add on a few gases (so we have gasIDs [1 2 3 4 5 6 9 12]) if needed
! old and new gases are ALWAYS going be in units = 12 (volume mixing ratio)
! have to be careful when adding in upper level info to old gases, as this could be a mix of units!!!!!!!!!!
! taper the tacked on profile so that within about 4 sublevels, you are just tacking on US Std (or default ref prof)
! can do the CO2 = 370 +yy-2002), CH4 = 1.7 + (yy-2002)*0.01
    CALL Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR,  &
        rPMin,rPMax,rYear)
    
! now check everything already is in VMR (gas unit 12); change to that if needed
! also adjust MR (wet water) if this is planet Earth
    CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)
    
    CALL InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,  &
        PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev)
    
    rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/  &
        (PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) 'iX = iLowestLev ==>'
    WRITE(kStdWarn,*)' PLEV_KCARTADATABASE_AIRS(iX),rPSurf,PLEV_KCARTADATABASE_AIRS(iX+1) = '
    WRITE(kStdWarn,*) PLEV_KCARTADATABASE_AIRS(iLowestLev),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev+1)
    WRITE(kStdWarn,*) 'rfracBot = ',rFracBot
    
! >>> do the integrals using integrals wrt p!!! closer to klayers
    CALL DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,  &
        PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,  &
        raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
    iNumLays = kProflayer-iLowestLev+1  !! this is now LAYER number, max == kProfLayer when iLowestLev = 1
    
    DO iL = 1,101
      raPbndFinal(iL) = raPbndFinal(iL)/100.0    !! convert N/m2 to mb
    END DO
    
  ELSE IF (iZbnd >= 0) THEN
! >>>>>>>>>> user boundaries
    DO iL = 1,iZbnd
      raPbndFinal(iL) = raPbnd(iL)/100.0    !! convert N/m2 to mb
    END DO
    DO iL = iZbnd+1,kProfLayer+1
      raPbndFinal(iL) = 0.005
      raPbndFinal(iL) = raPbnd(iZbnd)/100.0   !! convert N/m2 to mb
    END DO
    iZbndFinal = iZbnd
    
! now check everything already is in VMR (gas unit 12); change to that if needed
! also adjust MR (wet water) if this is planet Earth
    CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)
    
    CALL InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,  &
        PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev,iZbnd,raPBnd)
    
    rFracBot = (rPSurf-raPBnd(iLowestLev+1))/  &
        (raPBnd(iLowestLev)-raPBnd(iLowestLev+1))
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) 'iX = iLowestLev ==>'
    WRITE(kStdWarn,*)' raPBnd(iX),rPSurf,raPBnd(iX+1) = '
    WRITE(kStdWarn,*) raPBnd(iLowestLev),rPSurf,raPbnd(iLowestLev+1)
    WRITE(kStdWarn,*) 'rfracBot = ',rFracBot
    
! >>> do the integrals using integrals wrt p!!! closer to klayers
    CALL DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,  &
        raPBnd,raAlt,iZbnd,rfracBot,  &
        raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
    iNumLevs = iZbnd-iLowestLev+1
    iNumLays = iNumLevs - 1       !! this is now LAYER number, max == iZbnd-1 when iLowestLev = 1
    
  END IF
  
ELSE IF (kRTP == -6) THEN
  CALL DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon,  &
      PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,  &
      raPX,raTX,raLayDensityX,raaG_VMRX,  &
      raPressLevelsX,raTPressLevelsX,raAltitudesX,  &
      raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)
  rfracBot = 1.0
  iNumLevs = iNumLays+1
  iZbnd = iNumLevs
  iOffSet = (kProfLayer+1)-(iNumLays+1)
  DO iL = 1,iZbnd
    raPbndFinal(iL)    = raPressLevelsX(iL)         !! in mb
    raTbndFinal(iL)    = raTPressLevelsX(iL)        !! in K
    raZout(iL+iOffSet) = raAltitudesX(iL)*1000.0    !! reset raZout ... and use  m
  END DO
  iZbndFinal = iNumLevs
END IF

!      do iL = 1,iNumlevs
!        print *,'hmm4',iL,raP(iL),raT(iL),raaG_VMR(iL,1),raPX(iL),raTX(iL),raaG_VMRX(iL,1)
!      end do

RETURN
END SUBROUTINE InputMR_profile

!************************************************************************

SUBROUTINE Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,  &
    rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta)


NO TYPE, INTENT(IN OUT)                  :: iKCARTADir
NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
NO TYPE, INTENT(IN OUT)                  :: DATABASELE
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! output
REAL :: rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
INTEGER :: iKCARTADirectionUPorDOWN
REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1),DATABASELEVHEIGHTS(kMaxLayer+1)

! local
INTEGER :: iL

rPminKCarta = +1E10
rPmaxKCarta = -1E10
rHminKCarta = +1E10
rHmaxKCarta = -1E10
rPmin = +1E10
rPmax = -1E10
iKCARTADirectionUPorDOWN = +1   !! assume increading index == lower pressure ==> going UP
DO iL = 1,kMaxLayer+1
  PLEV_KCARTADATABASE_AIRS(iL) = PLEV_KCARTADATABASE_AIRS(iL) * 100.0  !! change mb --> N/m2
  DATABASELEVHEIGHTS(iL)       = DATABASELEVHEIGHTS(iL) * 1000.0       !! change km --> m
  IF (PLEV_KCARTADATABASE_AIRS(iL) > rPmaxKCarta) rPmaxKCarta = PLEV_KCARTADATABASE_AIRS(iL)
  IF (PLEV_KCARTADATABASE_AIRS(iL) < rPminKCarta) rPminKCarta = PLEV_KCARTADATABASE_AIRS(iL)
  IF (DATABASELEVHEIGHTS(iL) > rHmaxKCarta) rHmaxKCarta = DATABASELEVHEIGHTS(iL)
  IF (DATABASELEVHEIGHTS(iL) < rHminKCarta) rHminKCarta = DATABASELEVHEIGHTS(iL)
END DO
IF ((PLEV_KCARTADATABASE_AIRS(1) > PLEV_KCARTADATABASE_AIRS(2)) .AND.  &
      (DATABASELEVHEIGHTS(1) < DATABASELEVHEIGHTS(2))) THEN
  iKCARTADirectionUPorDOWN = +1   !! increasing index == lower pressure ==> going UP
ELSE IF ((PLEV_KCARTADATABASE_AIRS(1) < PLEV_KCARTADATABASE_AIRS(2)) .AND.  &
      (DATABASELEVHEIGHTS(1) > DATABASELEVHEIGHTS(2))) THEN
  iKCARTADirectionUPorDOWN = -1   !! decreasing index == higher pressure ==> going DOWN
  WRITE (kStdErr,*) 'Bummer : being lazy, wanted this to be going UP'
  CALL DoStop
ELSE
  WRITE (kStdErr,*) 'Inconsistency : database pressures in one dir, heights in another'
  CALL DoStop
END IF

RETURN
END SUBROUTINE Init_n_layers

!************************************************************************
! this subroutine takes raaG_VMR and adds on new gas profiles using raaR100MR
! only to rPmin < rP < rPmax

SUBROUTINE InterpNewGasProfiles_to_InputPlevels(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels,  &
    iNumLevs,rPmin,rPmax,raP,raaG_VMR)


INTEGER, INTENT(IN)                      :: iNumGases0
INTEGER, INTENT(IN OUT)                  :: iaG0(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: raR100Pres
REAL, INTENT(IN)                         :: raaR100MR(kMaxLayer+10,kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iRefLevels
INTEGER, INTENT(IN)                      :: iNumLevs
REAL, INTENT(IN)                         :: rPmin
REAL, INTENT(IN)                         :: rPmax
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(OUT)                        :: raaG_VMR(2*kProfLayer,kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input



REAL :: raR100Press(kMaxLayer+10)
! output


! local
INTEGER :: iG,iL,i1,i2,iJunk
REAL :: raTempMR(kMaxLayer+10),raTempP(kMaxLayer+10),raTempMRX(2*kMaxLayer),rJunk
CHARACTER (LEN=50) :: FMT

INTEGER :: NTOTAL
REAL :: V1S,V2S,DV

! first find pressures that span rPmin < rP < rPmax
IF (raR100Press(iRefLevels) > rPmin*100.0) THEN
  WRITE(kStdErr,*) 'iRefLevels = ',iRefLevels
  WRITE(kStdErr,*) 'oops min pressure in ref database = ',raR100Press(iRefLevels),' N/m2'
  WRITE(kStdErr,*) 'should be SMALLER than min pressure in supplied user profile = ',rPmin*100.0,' N/m2'
  CALL DoStop
END IF

IF (raR100Press(1) < rPmax*100.0) THEN
  WRITE(kStdErr,*) 'iRefLevels = ',iRefLevels
  WRITE(kStdErr,*) 'oops max pressure in ref database = ',raR100Press(1),' N/m2'
  WRITE(kStdErr,*) 'should be GREATER than max pressure in supplied user profile = ',rPmax*100,' N/m2'
  CALL DoStop
END IF

WRITE(kStdWarn,*) ' '

i1 = iRefLevels
10  CONTINUE
IF ((raR100Press(i1)/100.0 < rPmin) .AND. (i1 > 1)) THEN
  i1 = i1 - 1
  GO TO 10
END IF
i1 = MIN(i1 + 1,iReflevels)

i2 = 1
20  CONTINUE
IF ((raR100Press(i2)/100.0 > rPmax) .AND. (i2 < iRefLevels)) THEN
  i2 = i2 + 1
  GO TO 20
END IF
i2 = MAX(1,i2-1)

!      NTOTAL = 100
!      V1S = 1000.0
!      V2S = 1100.0
!      DV = 0.0025
!      Write(kStdWarn,'(/,A,/A,F12.5,/,A,F12.5,/,A,F12.5,/,A,I7)')
!     1    ' Adjusted Limits of Scanned Spectrum:',
!     2    ' V1 = ',V1S,' V2 = ',V2S,' DV = ',DV,' N  =',NTOTAL
!      Write(kStdWarn,'(/A,A,F12.5,A,F12.5,A,F12.5,A,I7)')
!     1    ' Adjusted Limits of Scanned Spectrum:',
!     2    ' V1 = ',V1S,' V2 = ',V2S,' DV = ',DV,' N  =',NTOTAL
!      FMT = '(/A,A,F12.5,A,F12.5,A,F12.5,A,I7)'
!      Write(kStdWarn,FMT)
!     1    ' Adjusted2 Limits of Scanned Spectrum:',
!     2    ' V1x = ',V1S,' V2x = ',V2S,' DVx = ',DV,' Nx  =',NTOTAL

FMT = '(A,F15.7,A,I3,A,F15.7)'
rJunk = raR100Press(1)/100.0
WRITE(kStdWarn,FMT) 'Ref Database max press ',rJunk,' (at level ',1,'); user max press (in mb) = ',rPmax
rJunk = raR100Press(iRefLevels)/100.0
WRITE(kStdWarn,FMT) 'Ref Database min press ',rJunk,' (at level ',iRefLevels,'); user min press (in mb) = ',rPmin
FMT = '(A,I3,I3,A)'
WRITE(kStdWarn,FMT) 'Ref profiles between levels ',i1,i2,' spans rPmin,rPmax'

DO iL = i2,i1
  raTempP(iL-i2+1) = raR100Press(iL)/100.0
END DO

!      print *,(raTempP(iL),iL=1,i1-i2+1)
!      print *,(raP(iL),iL=1,iNumLevs)

!      print *,(raaR100MR(1,iG),iG = 1,iNumGases)
DO iG = iNumGases0 + 1,iNumGases
  DO iL = i2,i1
    raTempMR(iL-i2+1) = raaR100MR(iL,iG)
  END DO
  CALL r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
!      DO iJunk = 1,min(iNumLevs,i1-i2+1)
!        print *,iJunk,raTempP(iJunk),raTempMR(iJunk),raP(iJunk),raTempMRX(iJunk)
!      END DO
  DO iL = 1,iNumLevs
    raaG_VMR(iL,iG) = raTempMRX(iL)
  END DO
END DO

RETURN
END SUBROUTINE InterpNewGasProfiles_to_InputPlevels

!************************************************************************
! this basically loops and reads in ref profile for gas iG, if necessary changes pressures to N/m2
! input gas list (from user)                      = iaG
! output gaslist (preset list = 1 2 3 4 5 6 9 12) = union(user list,preset list)

! called by Tack_on_profile

SUBROUTINE ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,iaPlanetMolecules,iPlanetMolecules,  &
    iaG,iaGasUnits,iNumGases,  &
    raR100Press,raR100Temp,raaR100MR,iRefLevels,laysORlevs)


NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
NO TYPE, INTENT(IN OUT)                  :: iaPlanetMo
NO TYPE, INTENT(IN OUT)                  :: iPlanetMol
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNumGases
NO TYPE, INTENT(IN OUT)                  :: raR100Pres
REAL, INTENT(OUT)                        :: raR100Temp(kMaxLayer+10)
REAL, INTENT(OUT)                        :: raaR100MR(kMaxLayer+10,kMaxGas)
INTEGER, INTENT(OUT)                     :: iRefLevels
NO TYPE, INTENT(IN OUT)                  :: laysORlevs
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input/output

INTEGER :: laysORlevs                                    !! +1 if read 100 US Standard Ref layers, -1 IF READ 50 AFGL levels
INTEGER :: iaPlanetMolecules(kMaxGas),iPlanetMolecules   !! important moecule list
REAL :: PLEV_KCARTADATABASE_AIRS(kmaxLayer+1)
! output
! these are the individual reference profiles, at ~ (kMaxLayer + 10) layers
REAL :: raR100Press(kMaxLayer+10)


! local
REAL :: raRx110Temp(kProfLayer+10),raRx110MR(kProfLayer+10),raRx110Press(kProfLayer+10)

INTEGER :: iL,iJ,iMid,ifloor,iErr,iG,iFound,iNumGases0
REAL :: rX,raR100Amt(kMaxLayer+10),raR100PartPress(kMaxLayer+10)
CHARACTER (LEN=80) :: caRefFName
CHARACTER (LEN=3) :: ca3
INTEGER :: iNumLevsX,iIndex
REAL :: raPPX(kMaxProfLayer),raQX(kMaxProfLayer) !!US Std layer ppress, amt
REAL :: raPX(kMaxProfLayer), raTX(kMaxProfLayer) !!US Std layer press, temp
REAL :: raMMX(kMaxProfLayer)
CHARACTER (LEN=50) :: FMT

!! do a union of iaG and iaPlanetMolecules
iNumGases0 = iNumGases
DO iG = 1,iPlanetMolecules
  iFound = -1
  DO iL = 1,iNumGases
    IF (iaG(iL) == iaPlanetMolecules(iG)) iFound = +1
  END DO
  IF (iFound < 0) THEN
    WRITE(kStdWarn,*)'gasID ',iaPlanetMolecules(iG),' not found in user list, adding .... '
    iNumGases = iNumGases + 1
    iaG(iNumGases) = iaPlanetMolecules(iG)
    iaGasUnits(iNumGases) = 12   !! VMR
  END IF
END DO
IF (iNumGases0 == iNumGases) THEN
  FMT = '(A,I2,A)'
  WRITE(kStdWarn,*) 'you had all ',iPlanetMolecules,' important molecules in the input profile'
ELSE IF (iNumGases0 < iNumGases) THEN
  FMT = '(A,I2,A,I2,A)'
  WRITE(kStdWarn,FMT) 'you had ',iPlanetMolecules,' input gas profiles; now have ',iNumGases,' gas profiles'
END IF

DO iG = 1,iNumGases
  IF ((kPlanet /= 3) .OR. (laysORlevs == +1)) THEN
!! this is old : use only (US) Standard Profile
! read kCARTA kProfLayer reference profile
    CALL FindReferenceName(caRefFName,iaG(iG),-1)
    CALL ReadRefProf(caRefFName,kMaxLayer,raR100Amt,  &
        raR100Temp,raR100Press,raR100PartPress,iErr)
    iRefLevels = kMaxLayer
! change pressures from atm to N/m2
    DO iL = 1,kMaxLayer
      raR100Press(iL)     = raR100Press(iL) * kAtm2mb*100.0
      raR100PartPress(iL) = raR100PartPress(iL) * kAtm2mb*100.0
      raaR100MR(iL,iG)    = raR100PartPress(iL)/raR100Press(iL)
    END DO
    
  ELSE IF ((kPlanet == 3) .AND. (laysORlevs == -1)) THEN
!ca3 = kcaLevsRefProf(1:3)
    ca3 = 'DNE'
    iIndex = INDEX(kcaLevsRefProf,'DNE')
    IF (iIndex > 0) THEN
      WRITE(kStdErr,*) 'oops : do not have a set of Levels Reference Profiles'
      CALL DoStop
    END IF
!! this is new : can use one of the (AFGL) models
    CALL ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iaG(iG),iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
    
!xtest    CALL getAFGL(-1,iaG(iG),raPX,raPPX,raTX,raQX) ! this is what "AddOnAFGLProfile_arblevels" uses
!xtest        DO iL = 1,kMaxProfLayer
!xtest          raMMX(iL) = raPPX(iL)/raPX(iL)
!xtest        END DO
    
    iRefLevels = iNumLevsx
    DO iL = 1,iNumLevsx
      raR100Press(iL)     = raRx110Press(iL)       !! already in N/m2
      raR100Temp(iL)      = raRx110Temp(iL)
      raaR100MR(iL,iG)    = raRx110MR(iL) /1.0E6   !! change from ppmv to MR
      raR100PartPress(iL) = raR100Press(iL) * raaR100MR(iL,iG)
!xtest          IF (iL .LE. 100) print *,iaG(iG),iL,raaR100MR(iL,iG),raMMX(iL),raaR100MR(iL,iG)/raMMX(iL)
    END DO
    DO iL = iNumLevsx+1,kMaxLayer
      raR100Press(iL)     = 0.0
      raaR100MR(iL,iG)    = 0.0
      raR100PartPress(iL) = 0.0
    END DO
  END IF
  
! make sure pressures are decreasing with index ie layers going higher and higher
  IF (raR100Press(1) < raR100Press(2)) THEN
!!! need to swap!!!
    iMid = ifloor(kMaxLayer/2.0)
    DO iL = 1,iMid
      iJ = kMaxLayer-iL+1
      rX = raR100Amt(iJ)
      raR100Amt(iJ) = raR100Amt(iL)
      raR100Amt(iL) = rX
      
      rX = raR100Temp(iJ)
      raR100Temp(iJ) = raR100Temp(iL)
      raR100Temp(iL) = rX
      
      rX = raR100Press(iJ)
      raR100Press(iJ) = raR100Press(iL)
      raR100Press(iL) = rX
      
      rX = raR100PartPress(iJ)
      raR100PartPress(iJ) = raR100PartPress(iL)
      raR100PartPress(iL) = rX
      
      rX = raaR100MR(iJ,iG)
      raaR100MR(iJ,iG) = raaR100MR(iL,iG)
      raaR100MR(iL,iG) = rX
    END DO
  END IF
END DO

RETURN
END SUBROUTINE ReadRefProf_Units_laysORlevs

!************************************************************************
! this subroutine changes the input levels profiles units from VMR to ppmv
!    >>>> and if planet Earth, adjusts mix ratios (wet water) <<<<
! then convers back to VMR

SUBROUTINE adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR)


INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iaGasUnits(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumLevs
INTEGER, INTENT(IN)                      :: iNumGases
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raaG_VMR(2*kProfLayer,kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input


! input/output


! local
INTEGER :: iL,iG,iDry2Wet,iWat,iFound,iQuiet
REAL :: rJunk1,rJunk2
CHARACTER (LEN=40) :: FMT

iFound = -1
iWat = 1
10   CONTINUE
IF (iaG(iWat) == 1) THEN
  iFound = +1
ELSE
  iWat = iWat + 1
  IF (iWat <= iNumGases) THEN
    GO TO 10
  END IF
END IF

IF ((kPlanet == 3) .AND. (iFound == 1)) THEN
  WRITE(kStdWarn,*) 'Planet = 3, EARTH .. gas ',iWat,' corresponds to water'
  IF (iWat /= 1) THEN
    WRITE(kStdErr,*) 'Planet = 3, EARTH .. found water but it should have been g1, not g',iWat
    CALL DoStop
  END IF
ELSE IF ((kPlanet == 3) .AND. (iFound /= 1)) THEN
  WRITE(kStdErr,*) 'Planet = 3, EARTH .. no input gas corresponds to water'
  CallDoStop
END IF

!! change to ppmv
FMT = '(A,/A,/A,/A)'
WRITE(kStdWarn,FMT) 'FINAL PROCESSING : Read in reference gas profiles, text input levels profile etc',  &
    '     (ref. profiles to augment input prof to 0.005 mb and/or add standard panet gases)',  &
    '  converting from VMR to ppmv so we can do wet water fix',  &
    'After this we can do the levels2layers integration'
iQuiet = +1
DO iG = 1,iNumGases
  rJunk1 = raaG_VMR(1,iG)
  IF (iaGasUnits(iG) /= 10) THEN
    CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_VMR,iQuiet)
    rJunk2 = raaG_VMR(1,iG)
    WRITE(kStdWarn,800) iG,iaG(iG),iaGasUnits(iG),rJunk1,rJunk2
  ELSE
    WRITE(kStdWarn,801) iG,iaG(iG),iaGasUnits(iG),rJunk1
  END IF
END DO
800  FORMAT('index=',I3,' gasID=',I3,' InputGasUnits=',I3,' origQ(level1)=',ES12.6,' finalQ_ppmv(level1)=',ES12.6)
801  FORMAT('index=',I3,' gasID=',I3,' InputGasUnits=',I3,' origQ(level1)=',ES12.6,' already in ppmv')

iDry2Wet = -1 !! do not adjust dry air to wet air MR
iDry2Wet = +1 !!        adjust dry air to wet air MR  DEFAULT

! IF (kRTP == -20) iDry2Wet = -1

IF (iDry2Wet < 0) THEN
!! simple change from ppmv (gas units 10) to VMR (which is gas units 12)
  
  WRITE(kStdWarn,*) 'converting back to ppmv from VMR, not doing wet water fix'
  DO iG = 1,iNumGases
    DO iL = 1,iNumLevs
      raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0E6
    END DO
  END DO
  
!! need to make sure gas units code == 10.11.12.20.21 before trying this; see toppmv.f in klayers
ELSE IF ((iDry2Wet > 0) .AND. (kPlanet == 3) .AND. (iaG(1) == 1)) THEN
  
!!! first convert WV mixing ratios
!!! however, klayars assumes the user has provided WET mixing ratios, so NO NEED to do this
!        DO iG = 1,1
!          DO iL = 1,iNumLevs
!            print *,iL,raaG_VMR(iL,iG),1e6*raaG_VMR(iL,iG) / (1.0e6+raaG_VMR(iL,iG))
!            raaG_VMR(iL,iG) =  1.0e6 * raaG_VMR(iL,iG) / (1.0e6 + raaG_VMR(iL,iG))
!          END DO
!        END DO
  
!!! then use WV mixing ratios to fix other gas mixing ratios if they are 10.11.12.20.21
  DO iG = 2,iNumGases
    IF ((iaGasUnits(iG) >= 10) .AND. (iaGasUnits(iG) <= 21)) THEN
      DO iL = 1,iNumLevs
        raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) * (1.0E6 - raaG_VMR(iL,1))/1E6
      END DO
    END IF
  END DO
  
!! finally convert ppmv --> VMR
  WRITE(kStdWarn,*) 'converting back to ppmv from VMR, after doing wet water fix'
  DO iG = 1,iNumGases
    DO iL = 1,iNumLevs
      raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0E6
    END DO
  END DO
  
ELSE
  WRITE(kStdErr,*) 'Need iDry2Wet == +/- 1'
  CALL DoStop
END IF

!! check all MR larger than 0
DO iG = 1,iNumGases
  DO iL = 1,kMaxLayer
    raaG_VMR(iL,iG) = MAX(raaG_VMR(iL,iG),0.0)
  END DO
END DO

RETURN
END SUBROUTINE adjustVMR_Earth

!************************************************************************
! this does the actual change for raaG_VMR(iL,iG)
! klayers/Doc/gas_units_code.txt
! klayers/Doc/toppmv.f
!              ------- 10-19 = volume mixing ratio ---------------------
!         10   parts per million volume mixing ratio (ppmv)
!              Number of gas X molecules per 1E+6 "air" molecules

!         11   parts per billion volume mixing ratio (ppbv)
!              Number of gas X molecules per 1E+9 "air" molecules

!         12   volume mixing ratio (unitless fraction)
!              Number of gas X molecules per "air" molecule

!              ------- 20-29 = mass mixing ratio -----------------------
!         20   mass mixing ratio in (g/kg)
!              Grams of gas X per kilogram of "air"

!         21   mass mixing ratio in (g/g)
!              Grams of gas X per gram of "air"

!              ------- 30-39 = partial pressure ------------------------

!         30   partial pressure in millibars (mb) Note: mb=hPa
!              Pressure of gas X as it exists in the atmosphere.
!              To clarify, this means at the corresponding profile
!              temperature and pressure level total pressure.

!         31   partial pressure in atmospheres (atm)
!              Pressure of gas X as it exists in the atmosphere

!              ------- 40-49 = water vapor humidity units --------------

!         40   relative humidity in (percent)
!              100 times actual vapor pressure divided by saturation
!              vapor pressure

!         41   relative humidity (unitless fraction)
!              Actual vapor pressure divided by saturation vapor
!              pressure

!         42   dew point temperature (Kelvin)
!              Temperaure at which vapor will start to condense out

!         43   dew point temperature (Celcius)
!              Temperaure at which vapor will start to condense out

!              Possible additional units might be
!                 inches or centimenters of water vapor
!                 grams of water vapor


! pIN is in mb

SUBROUTINE changeLVLS_2_ppmv(iGasID,iGasUnits,iNumLevs,iG,PIN,TIN,raaG_VMR,iQuiet)


INTEGER, INTENT(IN)                      :: iGasID
INTEGER, INTENT(IN)                      :: iGasUnits
INTEGER, INTENT(IN)                      :: iNumLevs
INTEGER, INTENT(OUT)                     :: iG
REAL, INTENT(IN)                         :: PIN(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: TIN(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raaG_VMR(2*kProfLayer,kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iQuiet
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input


! input/output


! local
INTEGER :: iL,iCode,NLEV
REAL :: MRIN(2*kProfLayer,kMaxGas),RJUNK,WEXSVP,MDAIR
INTEGER :: iaLocalGasID(kMaxGas)
REAL :: MASSF(kMaxGas)
CHARACTER (LEN=50) :: FMT
CHARACTER (LEN=40) :: caaUnit(50)
CHARACTER (LEN=20) :: cID

INCLUDE '../INCLUDE/gasIDnameparam.f90'
! see cbgids.f in klayers
DATA (iaLocalGasID(iL),iL=1,12) /01,02,03,04,05,06,07,08,09,10,11,12/
DATA (MASSF(iL),iL=1,12) /18.015,44.010,47.9982,44.013,28.011,16.043,31.999,30.006,64.063,46.006,17.031,63.013/

DO iL = 1,50
  caaUnit(iL) = '** Not Known **'
END DO
caaUnit(01) = 'layer amt molecules/cm2'
caaUnit(10) = 'vmr in ppmv'
caaUnit(11) = 'vmr in ppbv'
caaUnit(20) = 'mmr in g/kg'
caaUnit(21) = 'mmr in g/g'
caaUnit(30) = 'part press in mb (not accepted)'
caaUnit(31) = 'part press in atm (not accepted)'
caaUnit(40) = 'RH % (not accepted)'
caaUnit(41) = 'RH fraction (not accepted)'

IF (iGasID > 12) THEN
  MASSF(iGasID) = kAtmMolarMass
END IF

MDAIR = kAtmMolarMass

DO iL = 1,iNumLevs
  MRIN(iL,iG) = raaG_VMR(iL,iG)
END DO

! convert to ppmv (gas unit 10)
iCode = iGasUnits
NLEV  = iNumLevs

IF ((ICODE /= 10) .AND. (iQuiet < 0)) THEN
  FMT = '(I2,A,I2,A,A,A,I2,A,A,A)'
  cID = caGID(iGasID)
  WRITE(kStdWarn,FMT) iG,' : gID=',iGasID,'(',cID,') units code ',ICODE,'(',caaUnit(ICODE),') convert TO ppmv (UNIT code 10)'
END IF

IF ((ICODE == 10) .AND. (iQuiet < 0)) THEN
  FMT = '(I2,A,I2,A,A,A)'
  cID = caGID(iGasID)
  WRITE(kStdWarn,FMT) iG,' : gID=',iGasID,'(',cID,') already in ppmv (unit code 10)'
  
!            parts per billion volume mixing ratio
ELSE IF (ICODE == 11) THEN
!               PPMV = PPBV*1E-3
  DO IL=1,NLEV
    MRIN(IL,IG)=MRIN(IL,IG)*1E-3
    ENDDO
      
!            volume mixing ratio
    ELSE IF (ICODE == 12) THEN
!               PPMV = VMR*1E+6
      DO IL=1,NLEV
        MRIN(IL,IG)=MRIN(IL,IG)*1E+6
        ENDDO
          
!            mass mixing ratio in g/kg
        ELSE IF (ICODE == 20) THEN
!               PPMV = MRgkg*((MDAIR*1E-3)/MASSF)*1E+6
          RJUNK=1E+3*MDAIR/MASSF(IGasID)
          DO IL=1,NLEV
            MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
            ENDDO
              
!            mass mixing ratio in g/g
            ELSE IF (ICODE == 21) THEN
!               PPMV = MRgg*(MDAIR/MASSF)*1E+6
              RJUNK=1E+6*MDAIR/MASSF(IGasID)
              DO IL=1,NLEV
                MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
                ENDDO
                  
!            partial pressure in mb
                ELSE IF (ICODE == 30) THEN
!            PPMV = (PPmb/PIN)*1E+6
                  DO IL=1,NLEV
                    MRIN(IL,IG)=(MRIN(IL,IG)/PIN(IL))*1E+6
                    ENDDO
                      
!            partial pressure in atm
                    ELSE IF (ICODE == 31) THEN
!               PPMV = (PPatm*1013.25/PIN)*1E+6
                      DO IL=1,NLEV
                        MRIN(IL,IG)=(MRIN(IL,IG)*1013.25/PIN(IL))*1E+6
                        ENDDO
                          
!            relative humidy in percent
!            note we need to change PN fron N/m2 to mb, so divide by 100
                        ELSE IF (ICODE == 40 .AND. iGasID == 1) THEN
!            PPMV = (RH%/100)*(SVP/PIN)*1E+6
                          DO IL=1,NLEV
!                print *,iL,MRIN(IL,IG),TIN(iL),PIN(IL)/100.0,WEXSVP( TIN(IL))
                            MRIN(IL,IG)=MRIN(IL,IG)*  &
                                (WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+4
                            ENDDO
!                call dostop
                              
!            relative humidity (fraction)
                            ELSE IF (ICODE == 41 .AND. iGasID == 1) THEN
!               PPMV = RH*(SVP/PIN)*1E+6
                              DO IL=1,NLEV
                                MRIN(IL,IG)=MRIN(IL,IG)*  &
                                    (WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+6
                                ENDDO
                                  
                                ELSE
                                  WRITE(kStdErr,*) 'gasID,icode = ',iGasID,ICODE,' UNKNOWN COMBO to comvert to ppmv!!!'
                                  CALL DoStop
                                END IF
                                
! save
                                DO iL = 1,iNumLevs
                                  raaG_VMR(iL,iG) = MRIN(iL,iG)
                                END DO
                                
                                RETURN
                              END SUBROUTINE changeLVLS_2_ppmv
                              
!************************************************************************
! this subroutine tacks on Standard Profile above what the user has given, to fill in amounts till 0.005 mb
! also if kCARTA has earlier determined there should be more gases than in user input profile (eg suppose
! user only gave WV/O3 but kCARTA also wants other gases), then this is the place additional profiles added in
! also can do the CO2 = 370 ppm +(yy-2002)*2.2
! also can do the CH4 = 1.7 ppm +(yy-2002)*0.01
                              
                              SUBROUTINE Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_VMR,  &
                                  rPMin,rPMax,rYear)
                              
                              
                              NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
                              INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
                              INTEGER, INTENT(IN OUT)                  :: iaGasUnits(kMaxGas)
                              INTEGER, INTENT(IN OUT)                  :: iNumLevs
                              INTEGER, INTENT(IN)                      :: iNumGases
                              REAL, INTENT(OUT)                        :: raP(2*kProfLayer)
                              REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
                              REAL, INTENT(IN OUT)                     :: raaG_VMR(2*kProfLayer,kMaxGas)
                              NO TYPE, INTENT(IN)                      :: rPMin
                              NO TYPE, INTENT(IN)                      :: rPMax
                              REAL, INTENT(IN OUT)                     :: rYear
                              IMPLICIT NONE
                              
                              INCLUDE '../INCLUDE/kcartaparam.f90'
                              
! input var
                              
                              
                              REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! input/output var
                              REAL :: rPmin,rPmax
                              
                              
                              
! local var
!! individual reference profiles, at kMaxLayer layers, in terms of MR = PPress/Press
                              INTEGER :: iaG0(kMaxGas),iNumGases0,i2,iFound,iaGasUnits0(kMaxGas),iIndex
                              REAL :: raaR100MR(kMaxLayer+10,kMaxGas),raR100Temp(kMaxLayer+10),raR100Press(kMaxLayer+10)
                              INTEGER :: iL,iG,iK,iAbove,iMerge,iRefLevels,laysORlevs
                              REAL :: rToffset,rJunk,raOffset(kMaxGas),raJunk(kMaxLayer+10)
                              REAL :: rCO2ppmv,rAdjust,rCH4ppbv
                              CHARACTER (LEN=3) :: ca3
                              CHARACTER (LEN=50) :: FMT
                              INTEGER :: iEarth,iMars,iaEarth(8),iaMars(2),iaPlanetMolecules(kMaxGas),iPlanetMolecules
                              
                              DATA (iaEarth(iG),iG=1,8) /01,02,03,04,05,06,09,12/      !! 8 important Earth atm molecules
                              DATA (iaMars(iG),iG=1,2)  /02,22/                        !! 2 important Mars  atm molecules
                              iEarth = 8
                              iMars = 2
                              DO iG = 1,kMaxGas
                                iaPlanetMolecules(iG) = -1
                              END DO
                              IF (kPlanet == 03) THEN
                                iPlanetMolecules = iEarth
!! earth
                                DO iG = 1,iEarth
                                  iaPlanetMolecules(iG) = iaEarth(iG)
                                END DO
                              ELSE IF (kPlanet == 04) THEN
                                iPlanetMolecules = iMars
!! mars
                                DO iG = 1,iMars
                                  iaPlanetMolecules(iG) = iaMars(iG)
                                END DO
                              ELSE
                                WRITE (kStdErr,*) 'oops right now can only handle Earth or Mars. Sorry, Mork'
                                CALL DoStop
                              END IF
                              
! first read in reference profile, and if needed info for other gases
!     read in additional gases eg from ERA/ECM we get gas 1,3 but we would like gas 1,2,3,4,5,6,8,12
!     read in ref profiles (eg for Gas2 since Mars, Earth, Venus all have CO2 in the atm; get the Standard VMR)
                              iNumGases0 = iNumGases
                              DO iL = 1,iNumGases0
                                iaG0(iL) = iaG(iL)
                                iaGasUnits0(iL) = iaGasUnits(iL)
                              END DO
                              laysORlevs = +1   !! orig, read in reference P,PP and get mix ratio using PP/P fOR each gas
                              laysORlevs = -1   !! new,  read in one of 6 AFGL P/T/ppmv level profiles
!ca3 = kcaLevsRefProf(1:3)
                              ca3 = 'DNE'
                              iIndex = INDEX(kcaLevsRefProf,'DNE')
                              IF ((laysORlevs == -1) .AND. (iIndex > 0)) THEN
                                WRITE(kStdWarn,*) 'oops : do not have a set of Levels Reference Profiles; code is setting laysORlevs = +1'
                                laysORlevs = +1
                              END IF
                              
                              CALL ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,iaPlanetMolecules,iPlanetMolecules,  &
                                  iaG,iaGasUnits,iNumGases,raR100Press,raR100Temp,raaR100MR,iRefLevels,laysORlevs)
                              
                              IF (iNumGases0 /= iNumGases) THEN
                                WRITE(kStdWarn,*)'Orig num of gases, new num of gases = ',iNumGases0,iNumGases
!! new gases are all at VMR (dimensionless fraction) = gasunit 12
!! these new profiles are AFGL profiles from glatm, and NOT at sonde or AIRS 101 levels
!! so they need to be interpolated from (raR100Press,raaR100MR) into raaG_VMR (at sonde P levels between rPmax and rPmin)
                                CALL InterpNewGasProfiles_to_InputPlevels(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels,  &
                                    iNumLevs,rPmin,rPmax,raP,raaG_VMR)
                              END IF
                              
! then find KCARTA levels above which there is NO user supplied info
                              FMT = '(A,F12.5,F12.5,I3)'
                              WRITE(kStdWarn,FMT)' user supplied profile : Pmax(mb),Pmin(mb),numlevs = ',rPMax,rPMin,iNumLevs
                              WRITE(kStdWarn,FMT)' kCARTA Pav Database   : Pmax(mb),Pmin(mb),numlevs = ',  &
                                  raR100Press(1)/100.0,raR100Press(iRefLevels)/100.0,iRefLevels
                              
                              iAbove = kMaxLayer
                              iAbove = iRefLevels
                              10  CONTINUE
                              IF ((raR100Press(iAbove)/100.0 < rPMin) .AND. (iAbove > 1)) THEN
                                iAbove = iAbove - 1
                                GO TO 10
                              END IF
                              iAbove = iAbove + 1
                              IF (iAbove > iRefLevels) THEN
                                FMT = '(A,I3,A,I3)'
                                WRITE(kStdWarn,FMT) 'tacking on ref  profile info from level ',iAbove,' to ',iRefLevels
                                WRITE(kStdWarn,*) 'Hmm ... no need to do this ... looks like user supplied levels profile past 0.005 mb'
                                GO TO 123
                              END IF
                              
                              FMT = '(A,I3,A,I3)'
                              WRITE(kStdWarn,FMT) 'tacking on ref profile info from level ',iAbove,' to ',iRefLevels
                              WRITE(kStdWarn,FMT) 'which should extend user supplied level info from level ',iNumLevs,' TO ',  &
                                  iNumLevs + (iRefLevels-iAbove+1)
                              
! now do linear interpolation from iAbove pressure, down to min(raP)
                              CALL r_sort_loglinear1(raR100Press,raR100Temp,iRefLevels,rPMin*100.0,rJunk,1)
                              rToffset = raT(iNumLevs) - rJunk
                              FMT = '(A,I3,A,I3)'
                              WRITE(kStdWarn,FMT)'tacking on info from kCARTA Pav Dtabase, layers ',iAbove,' to ',iRefLevels
                              WRITE(kStdWarn,*)'ToffSet = ',rToffset
!      iK = -1
!      write(kStdWarn,*) rPMin*100.0,iNumLevs,raT(iNumLevs)
!      DO iL = 1,iRefLevels
!        IF ((raR100Press(iL) .LE. rPMin*100.0) .AND. (iK .EQ. -1)) THEN
!        write(kStdWarn,*) iL,raR100Press(iL)/100,raR100Temp(iL),'*****',raT(iNumLevs),rToffset
!      ELSE
!        write(kStdWarn,*) iL,raR100Press(iL)/100,raR100Temp(iL)
!      END IF
!     END DO
!     call dostop
                              
! also need to figure out the gas multiplier offset
                              DO iG = 1, iNumGases
                                DO iL = 1,iRefLevels
                                  raJunk(iL) = raaR100MR(iL,iG)
                                END DO
                                CALL r_sort_loglinear1(raR100Press,raJunk,iRefLevels,rPMin*100.0,rJunk,1)
!! need to do a CHANGE OF UNITS to ppmv!!!!!
                                raoffset(iG) = raaG_VMR(iNumLevs,iG)/rJunk
                                IF (raoffset(iG) < 0) THEN
!        print *,iG,iaG(iG)
!          iK = -1
!          write(kStdWarn,*) rPMin*100.0,iNumLevs,raaG_VMR(iNumLevs,iG)
!          DO iL = 1,iRefLevels
!            IF ((raR100Press(iL) .LE. rPMin*100.0) .AND. (iK .EQ. -1)) THEN
!            write(kStdWarn,*) iL,raR100Press(iL)/100,raJunk(iL),'*****',raT(iNumLevs),rJunk
!              ELSE
!            write(kStdWarn,*) iL,raR100Press(iL)/100,raJunk(iL)
!          END IF
!          END DO
                                  raoffset(iG) = 1.0
                                END IF
                                WRITE(kStdWarn,*) 'Multiplier for gasID = ',iaG(iG),' units = ',iaGasUnits(iG),' = ',raoffset(iG)
                              END DO
                              
                              iMerge = +1    !! in other words, 3 points away from iMerge, we bump up/down the multiplier TO 1.000
                              DO iL = iAbove,iRefLevels
                                iNumLevs = iNumLevs + 1
                                raP(iNumLevs) = raR100Press(iL)/100.0
                                
!raT(iNumLevs) = raR100Temp(iL) + rToffset
                                IF ((iL-iAbove) > iMerge) THEN
                                  rAdjust = 0.0  !! make the temp adjustment effectively 0
                                ELSE
                                  rJunk = (0 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 0, inside 4 levels; this is slope
                                  rAdjust = rJunk * (iL-iAbove) + rToffset
                                  rAdjust = rToffset  !!!!!!!! ????
                                END IF
                                raT(iNumLevs) = raR100Temp(iL) + rAdjust
                                
                                DO iG = 1,iNumGases
!raaG_VMR(iNumLevs,iG) = raaR100MR(iL,iG) * raOffset(iG)
                                  IF ((iL-iAbove) > iMerge) THEN
                                    rAdjust = 1  !! make the MR adjustment effectively 1
                                  ELSE
                                    rJunk = (1 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 1, inside 4 levels; this is slope
                                    rAdjust = rJunk * (iL-iAbove) + raOffset(iG)
                                    rAdjust = raOffset(iG)   !!!! ????
!            print *,iG,iL,iL-iAbove,raOffset(iG),rAdjust
                                  END IF
                                  raaG_VMR(iNumLevs,iG) = raaR100MR(iL,iG) * rAdjust
                                END DO
                              END DO
                              WRITE(kStdWarn,*) 'Have extended user supplied info to level',iNumLevs
                              
                              123  CONTINUE
!      DO iL = 1,iRefLevels
!        print *,'T profile ',iL,raP(iL),raT(iL)
!      END DO
!      call dostop
                              
! >>>>>>>>>>>>>>>>>>>>>>>>> CO2 adjust; growth rate = 2 ppmv/yr
                              rCO2ppmv = 370.0
                              IF ((rYear > 1970) .AND. (kPlanet == 03)) THEN
!! should be dt = (yy-2002) + mm/12;  rCO2ppmv = 370 + 2*dt;
                                rCO2ppmv = 370.0 + (rYear-2002)*2.0
                                
!! now see if CO2 was originally there
                                iFound = -1
                                i2 = 1
                                20     CONTINUE
                                IF (iaG0(i2) == 2) THEN
                                  iFound = +1
                                ELSE IF (i2 < iNumGases0) THEN
                                  i2 = i2 + 1
                                  GO TO 20
                                END IF
                                
                                IF (iFound > 0) THEN
!! user had CO2 in supplied prof, no need to adjust MR
                                  WRITE(kStdWarn,*) 'User has CO2 in supplied profile, so not adjusting bottom levels'
                                ELSE
!! user did not have CO2 in there, yes need to adjust MR
                                  i2 = 1
                                  iFound = -1
                                  25       CONTINUE
                                  IF (iaG(i2) == 2) THEN
                                    iFound = +1
                                  ELSE IF (i2 < iNumGases) THEN
                                    i2 = i2 + 1
                                    GO TO 25
                                  END IF
                                  IF (iFound < 0) THEN
                                    WRITE(kStdErr,*) 'hmm, expected to find CO2!!! ooops'
                                    CALL DoStop
                                  END IF
                                  rJunk = 0.0
                                  DO iL = 3,7
                                    rJunk = rJunk + raaG_VMR(iL,i2)
                                  END DO
                                  rJunk = rJunk/5.0
                                  WRITE(kStdWarn,*) '  reference CO2 profile that was read in has CO2 MR = ',rJunk*1.0E6,' ppmv'
                                  WRITE(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCO2ppmv,' ppmv'
                                  DO iL = 1,iNumLevs
                                    raaG_VMR(iL,i2) = raaG_VMR(iL,i2) * rCO2ppmv/(rJunk*1.0E6)
                                  END DO
                                END IF
                              END IF
                              
! >>>>>>>>>>>>>>>>>>>>>>>>> CH4 adjust; growth rate = 0.01 ppmb/yr
                              rCH4ppbv = 1.7
                              IF ((rYear > 1970) .AND. (kPlanet == 03)) THEN
                                rCH4ppbv = 1.7 + (rYear-2002)*0.01
                                rCH4ppbv = 1.8   !! hard coded into klayers right now
                                
!! now see if CH4 was originally there
                                iFound = -1
                                i2 = 1
                                30     CONTINUE
                                IF (iaG0(i2) == 6) THEN
                                  iFound = +1
                                ELSE IF (i2 < iNumGases0) THEN
                                  i2 = i2 + 1
                                  GO TO 30
                                END IF
                                
                                IF (iFound > 0) THEN
!! user had CH4 in supplied prof, no need to adjust MR
                                  WRITE(kStdWarn,*) 'User has CH4 in supplied profile, so not adjusting bottom levels'
                                ELSE
!! user did not have CH4 in there, yes need to adjust MR
                                  i2 = 1
                                  iFound = -1
                                  35       CONTINUE
                                  IF (iaG(i2) == 6) THEN
                                    iFound = +1
                                  ELSE IF (i2 < iNumGases) THEN
                                    i2 = i2 + 1
                                    GO TO 35
                                  END IF
                                  IF (iFound < 0) THEN
                                    WRITE(kStdErr,*) 'hmm, expected to find CH4!!! ooops'
                                    CALL DoStop
                                  END IF
                                  rJunk = 0.0
                                  DO iL = 3,7
                                    rJunk = rJunk + raaG_VMR(iL,i2)
                                  END DO
                                  rJunk = rJunk/5.0
                                  WRITE(kStdWarn,*) '  reference CH4 profile that was read in has CH4 MR = ',rJunk*1.0E6,' ppbv'
                                  WRITE(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCH4ppbv,' ppbv'
                                  DO iL = 1,iNumLevs
                                    raaG_VMR(iL,i2) = raaG_VMR(iL,i2) * rCH4ppbv/(rJunk*1.0E6)
                                  END DO
                                END IF
                              END IF
                              
!>>>>>>>>>>>>>>>>>>>>>>>>>
                              
!! check all MR larger than 0
                              DO iG = 1, iNumGases
                                DO iL = 1,kMaxLayer
                                  raaG_VMR(iL,iG) = MAX(raaG_VMR(iL,iG),0.0)
                                END DO
!      print *,iG,iaG(iG),raaG_VMR(1,iG)
                              END DO
                              
                              RETURN
                            END SUBROUTINE Tack_on_profile
                            
!************************************************************************
! this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
                            
                            SUBROUTINE InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,  &
                                PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev)
                            
                            
                            INTEGER, INTENT(IN)                      :: iNumLevs
                            INTEGER, INTENT(IN)                      :: iNumGases
                            INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
                            REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
                            REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
                            REAL, INTENT(IN)                         :: raaG_VMR(2*kProfLayer,kMaxGas)
                            NO TYPE, INTENT(IN OUT)                  :: rPMin
                            NO TYPE, INTENT(IN OUT)                  :: rPMax
                            REAL, INTENT(OUT)                        :: rPSurf
                            REAL, INTENT(IN OUT)                     :: rTSurf
                            NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
                            NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
                            REAL, INTENT(OUT)                        :: raPX(kProfLayer+1)
                            REAL, INTENT(OUT)                        :: raTX(kProfLayer+1)
                            REAL, INTENT(OUT)                        :: raaG_VMRX(kProfLayer+1,kMaxGas)
                            INTEGER, INTENT(OUT)                     :: iLowestLev
                            IMPLICIT NONE
                            
                            INCLUDE '../INCLUDE/kcartaparam.f90'
                            
! input var
                            
                            REAL :: rPmin,rPmax, rPmaxKCarta
                            
                            REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
! output var
                            
                            
                            REAL :: rInJunk,rOutJunk
                            
! local var
                            INTEGER :: iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
                            REAL :: raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_VMRX2(kProfLayer+1,kMaxGas)
                            
                            WRITE(kStdWarn,*) ' ----->>>> using default AIRS 101 levels for layers integration'
                            
! see which kCARTA pressure level is equal to or greater than surface pressure
                            WRITE(kStdWarn,*) '--------------------------------------------------------------------------'
                            WRITE(kStdWarn,*) ' '
                            WRITE(kStdWarn,*) ' >>> prep : klayers within kcarta : interp usergrid onto 101 levels grid >>>'
                            
                            iLowestLev = kMaxLayer
                            10   CONTINUE
                            IF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) < rPSurf) .AND. (iLowestLev > 1)) THEN
                              iLowestLev = iLowestLev - 1
                              GO TO 10
                            ELSE IF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) < rPSurf) .AND. (iLowestLev == 1)) THEN
                              WRITE (kStdErr,*) 'PSurf = ',rPSurf,' and Max Kcarta Database Plev = ',rPmaxKCarta
                              CALL DoStop
                            ELSE IF (PLEV_KCARTADATABASE_AIRS(iLowestLev) >= rPSurf) THEN
                              WRITE(kStdWarn,*)'lowest level = ',iLowestLev
                            END IF
                            
! fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
                            DO iL = 1,kMaxLayer-iLowestLev+1
                              iNumUse = iL
                              raPX(iL) = PLEV_KCARTADATABASE_AIRS(iLowestLev+iL)/100.0  ! input press in mb, so change plev_kcarta TO mb
                            END DO
                            
! now interpolate T and MR onto this grid
                            CALL r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
                            DO iG = 1,iNumGases
                              DO iL = 1,iNumLevs
                                raTemp(iL) = raaG_VMR(iL,iG)
                              END DO
                              CALL r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
                              DO iL = 1,iNumUse
                                raaG_VMRX(iL,iG) = raTempX(iL)
                              END DO
                            END DO
                            
!      DO iL = 1,kMaxLayer-iLowestLev+1
!        write(*,1234) iL,raP(iL),raT(iL),raaG_VMR(iL,1),raPX(iL),raTX(iL),raaG_VMRX(iL,1)
!      END DO
!      call dostop
                            1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
                            
! need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
                            iHigh = iNumUse
                            20   CONTINUE
                            IF (raP(iNumLevs) > raPX(iHigh)) THEN
                              iHigh = iHigh - 1
                              GO TO 20
                            END IF
                            iHigh = iHigh + 1
                            IF (iHigh <= iNumUse) THEN
                              WRITE(kSTdWarn,*) ' '
                              WRITE(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
                              WRITE(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
                              WRITE(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' TO ',iNumUse
                              DO iL = iHigh,iNumUse
                                rInJunk = raPX(iL)
                                CALL r_sort_loglinear1(raP,raT,iNumLevs,rInJunk,rOutJunk,1)
                                raTX(iL) = rOutJunk
                                
                                DO iG = 1,iNumGases
                                  DO iX = 1,iNumLevs
                                    raTemp(iX) = raaG_VMR(iX,iG)
                                  END DO
                                  CALL r_sort_loglinear1(raP,raTemp,iNumLevs,rInJunk,rOutJunk,1)
                                  raaG_VMRX(iL,iG) = rOutJunk
                                END DO
                              END DO
                            END IF
                            
! now sort the gasIDs into increasing order, and do the same for the mixing ratios
                            DO iL = 1,iNumGases
                              raArr(iL) = iaG(iL) * 1.0
                            END DO
                            CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)
                            
! assign temp arrays
                            DO iG = 1,iNumGases
                              iaG2(iG) = iaG(iG)
                            END DO
                            DO iL = 1,iNumUse
                              DO iG = 1,iNumGases
                                raaG_VMRX2(iL,iG) = raaG_VMRX(iL,iG)
                              END DO
                            END DO
                            
! sort according to iaIndex
                            DO iG = 1,iNumGases
                              iaG(iG) = iaG2(iaIndex(iG))
                            END DO
                            DO iL = 1,iNumUse
                              DO iG = 1,iNumGases
                                raaG_VMRX(iL,iG) = raaG_VMRX2(iL,iaIndex(iG))
                              END DO
                            END DO
                            
                            WRITE(kStdWarn,*) '--------------------------------------------------------------------------'
                            
                            RETURN
                          END SUBROUTINE InterpUser2kCARTA
                          
!************************************************************************
! this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
                          
                          SUBROUTINE InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_VMR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,  &
                              PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_VMRX,iLowestLev,iZbnd,raPbnd)
                          
                          
                          INTEGER, INTENT(IN)                      :: iNumLevs
                          INTEGER, INTENT(IN)                      :: iNumGases
                          INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
                          REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
                          REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
                          REAL, INTENT(IN)                         :: raaG_VMR(2*kProfLayer,kMaxGas)
                          NO TYPE, INTENT(IN OUT)                  :: rPMin
                          NO TYPE, INTENT(IN OUT)                  :: rPMax
                          REAL, INTENT(IN OUT)                     :: rPSurf
                          REAL, INTENT(IN OUT)                     :: rTSurf
                          NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
                          NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
                          REAL, INTENT(OUT)                        :: raPX(kProfLayer+1)
                          REAL, INTENT(OUT)                        :: raTX(kProfLayer+1)
                          REAL, INTENT(OUT)                        :: raaG_VMRX(kProfLayer+1,kMaxGas)
                          INTEGER, INTENT(OUT)                     :: iLowestLev
                          INTEGER, INTENT(IN)                      :: iZbnd
                          REAL, INTENT(IN)                         :: raPbnd(2*kProfLayer)
                          IMPLICIT NONE
                          
                          INCLUDE '../INCLUDE/kcartaparam.f90'
                          
! input var
                          
                          REAL :: rPmin,rPmax, rPmaxKCarta
                          
                          REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
                          REAL :: !! user defined boundaries
                          
! output var
                          
                          REAL :: rInJunk,rOutJunk
                          
! local var
                          INTEGER :: iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
                          REAL :: raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_VMRX2(kProfLayer+1,kMaxGas)
                          
                          WRITE(kStdWarn,*) ' ----->>>> using user defined ',iZbnd,' pressure levels for layers integration'
                          
! see which user defined pressure level is equal to or greater than surface pressure
                          WRITE(kStdWarn,*) '--------------------------------------------------------------------------'
                          WRITE(kStdWarn,*) ' '
                          WRITE(kStdWarn,*) ' >>> prep : user defined levels : interp usergrid onto N levels grid >>>'
                          
                          rPmaxKCarta = -1.0
                          DO iL = 1,iZbnd
                            PRINT *,'in InterpUser2UserBnd',iL,raPbnd(iL)
                            IF (raPbnd(iL) > rPmaxKCarta) rPmaxKCarta = raPbnd(iL)
                          END DO
                          WRITE(kStdWarn,*) 'User defined levels for integration : maxP = ',rPmaxKCarta,' Psurf = ',rPSurf,' mb'
                          
!! assume raPbnd(1) > raPbnd(2) > raPbnd(3) > ... > raPbnd(iZbnd) and look
!! for raPbnd greater than raPSurf
                          iLowestLev = iZbnd
                          10   CONTINUE
                          IF ((raPbnd(iLowestLev) < rPSurf) .AND. (iLowestLev > 1)) THEN
                            iLowestLev = iLowestLev - 1
                            GO TO 10
                          ELSE IF ((raPbnd(iLowestLev) < rPSurf) .AND. (iLowestLev == 1)) THEN
                            WRITE (kStdErr,*) 'Error!!!  PSurf = ',rPSurf,' and Max UserLev Plev = ',rPmaxKCarta
                            CALL DoStop
                          ELSE IF (raPbnd(iLowestLev) >= rPSurf) THEN
                            WRITE(kStdWarn,*)'lowest level = ',iLowestLev
                          END IF
                          
! fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
! remember this comes from LBLRTM input profile, so lowest layer will probably CLOSELY or EXACTLY match pSurf
                          DO iL = 1,iZbnd
                            raPX(iL) = raPbnd(iL)/100.0  ! input raPbnd press in N/m2, so change plev_kcarta TO mb
                          END DO
                          iNumUse = iZbnd
                          
! now interpolate T and MR onto this grid
                          CALL r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
                          DO iG = 1,iNumGases
                            DO iL = 1,iNumLevs
                              raTemp(iL) = raaG_VMR(iL,iG)
                            END DO
                            CALL r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
                            DO iL = 1,iNumUse
                              raaG_VMRX(iL,iG) = raTempX(iL)
                            END DO
                          END DO
                          
!      DO iL = 1,iZbnd-iLowestLev+1
!        write(*,1234) iL,raP(iL),raT(iL),raaG_VMR(iL,1),raPX(iL),raTX(iL),raaG_VMRX(iL,1)
!      END DO
!      call dostop
                          1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
                          
! need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
                          iHigh = iNumUse
                          20   CONTINUE
                          IF (raP(iNumLevs) > raPX(iHigh)) THEN
                            PRINT *,'wah iHigh',iNumLevs,iHigh,raP(iNumLevs),raPX(iHigh)
                            iHigh = iHigh - 1
                            GO TO 20
                          END IF
                          
                          iHigh = iHigh + 1
                          IF (iHigh <= iNumUse) THEN
                            WRITE(kSTdWarn,*) ' '
                            WRITE(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
                            WRITE(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
                            WRITE(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' TO ',iNumUse
                            DO iL = iHigh,iNumUse
                              rInJunk = raPX(iL)
                              CALL r_sort_loglinear1(raP,raT,iNumLevs,rInJunk,rOutJunk,1)
                              raTX(iL) = rOutJunk
                              
                              DO iG = 1,iNumGases
                                DO iX = 1,iNumLevs
                                  raTemp(iX) = raaG_VMR(iX,iG)
                                END DO
                                CALL r_sort_loglinear1(raP,raTemp,iNumLevs,rInJunk,rOutJunk,1)
                                raaG_VMRX(iL,iG) = rOutJunk
                              END DO
                            END DO
                          END IF
                          
! now sort the gasIDs into increasing order, and do the same for the mixing ratios
                          DO iL = 1,iNumGases
                            raArr(iL) = iaG(iL) * 1.0
                          END DO
                          CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)
                          
! assign temp arrays
                          DO iG = 1,iNumGases
                            iaG2(iG) = iaG(iG)
                          END DO
                          DO iL = 1,iNumUse
                            DO iG = 1,iNumGases
                              raaG_VMRX2(iL,iG) = raaG_VMRX(iL,iG)
                            END DO
                          END DO
                          
! sort according to iaIndex
                          DO iG = 1,iNumGases
                            iaG(iG) = iaG2(iaIndex(iG))
                          END DO
                          DO iL = 1,iNumUse
                            DO iG = 1,iNumGases
                              raaG_VMRX(iL,iG) = raaG_VMRX2(iL,iaIndex(iG))
                            END DO
                          END DO
                          
                          WRITE(kStdWarn,*) '--------------------------------------------------------------------------'
                          
                          RETURN
                        END SUBROUTINE InterpUser2UserBnd
                        
!************************************************************************
! function to compute vapor pressure wexp
                        
                        REAL FUNCTION WEXSVP(T)
                        
                        
                        REAL, INTENT(IN OUT)                     :: T
                        IMPLICIT NONE
                        
                        
                        DOUBLE PRECISION :: TEMP,LOGP
                        
                        TEMP = DBLE(T)
                        LOGP = -0.58002206D+4 / TEMP + 0.13914993D+1  &
                            - 0.48640239D-1 * TEMP + 0.41764768D-4 * TEMP**2  &
                            - 0.14452093D-7 * TEMP**3 + 0.65459673D+1 * DLOG(TEMP)
                        WEXSVP = SNGL( 1.0D-2 * DEXP( LOGP ) )
                        
!      Hyland & Wexler, 1983, equation for SVP over ice
!           log Pi =  -0.56745359E+4 / T
!                    + 0.63925247E+1
!                    - 0.96778430E-2  T
!                    + 0.62215701E-6  T^2
!                    + 0.20747825E-8  T^3
!                    - 0.94840240E-12 T^4
!                    + 0.41635019E+1 log(T)
!       with T in [K] and Pi in [Pa]
                        
                        RETURN
                      END FUNCTION WEXSVP
                      
!************************************************************************
! copied from klayers/grav.f
!CALL PROTOCOL:
!    GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)
                      
!INPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      LAT     latitude                    degrees
!    REAL      LON     longitude                   degrees
!    REAL      WINDE   wind velecity east          m/s
!    REAL      WINDN   wind velocity north         m/s
!    REAL      Z       altitude                    m
                      
!OUTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      GRAV    Earth gravity               m/s^2
                      
!INPUT/OUTPUT PARAMETERS: none
                      
!RETURN VALUES: none
                      
!PARENT(S):
!    INTLEV
                      
!ROUTINES CALLED: none
                      
!FILES ACCESSED: none
                      
!COMMON BLOCKS: none
                      
!DESCRIPTION:
!    Function to calculate Earth gravity (gravitation plus
!    centripetal acceleration) for points in the atmosphere.
                      
!    It calculates surface gravity using an equation given in
!    "American Institute of Physics Handbook", 1963, page 2-102.
!    This equation is essentially a variant of the International
!    Gravity Formula with an extra term for longitude (which is
!    very minor).
                      
!    Centripetal acceleration is tangental velocity squared over the
!    radius.
                      
!    Gravitation is the gravitational constant times the Earth's mass
!    divided by the square of the radius.
                      
!    Gravity at any point in the atmosphere is simply surface gravity
!    plus the change in gravitation and centripetal acceleration at
!    that point compared to the surface.
                      
!ALGORITHM REFERENCES: see DESCRIPTION
                      
!KNOWN BUGS AND LIMITATIONS: none
                      
!ROUTINE HISTORY:
!    Date     Programmer        Comments
!------------ ----------------- ----------------------------------------
! Mar  8 1995 Scott Hannon/UMBC created
! Jun 23 1995 Scott Hannon      Correct some comments
! 28 Sep 2007 Scott Hannon      Add more centripetel term comments
                      
!END====================================================================
                      
!      =================================================================
                      
                      REAL FUNCTION GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)
!      =================================================================
                      
!-----------------------------------------------------------------------
!      IMPLICIT NONE
!-----------------------------------------------------------------------
                      
                      REAL, INTENT(IN)                         :: Z
                      REAL, INTENT(IN)                         :: WINDE
                      REAL, INTENT(IN OUT)                     :: WINDN
                      REAL, INTENT(IN)                         :: LAT
                      REAL, INTENT(IN)                         :: LON
                      IMPLICIT NONE
                      
!-----------------------------------------------------------------------
!      INCLUDE FILES
!-----------------------------------------------------------------------
!      none
                      
!-----------------------------------------------------------------------
!      EXTERNAL FUNCTIONS
!-----------------------------------------------------------------------
!      none
                      
!-----------------------------------------------------------------------
!      ARGUMENTS
!-----------------------------------------------------------------------
                      
                      
                      
!-----------------------------------------------------------------------
!      LOCAL VARIABLES
!-----------------------------------------------------------------------
                      
                      REAL :: G_SUR, R, COSLT, COSLT2, SINLT2, SIN2LT, COSLON,  &
                          LTRAD, C_SUR, C_Z, RTOT, GRAVZ
                      
!      Constants for (1 - b^2/a^2) with
!      a = 6.378388E+6 m = equatorial radius, and
!      b = 6.356911E+6 m = polar radius.
                      REAL :: B2, ABTERM
                      
!      Constants for normal gravity equation
!      (see "American Institute of Physics Handbook", 1963, pg 2-102)
                      REAL :: G0
                      REAL :: C1, C2, C3
                      
!      Constants for pi/180, 2*pi, and Earth's rotational speed
!      of w=1/86400 rev/s
                      REAL :: PI180,PI2, W
                      
                      DATA B2, ABTERM /4.041031E+13, 6.724285E-3/
                      DATA G0 /9.780455/
                      DATA C1,C2,C3 /5.30157E-3, -5.85E-6, 6.40E-6/
                      DATA PI180,PI2,W /1.7453293E-2, 6.28318531, 1.1574074E-5/
                      
!-----------------------------------------------------------------------
!      SAVE STATEMENTS
!-----------------------------------------------------------------------
!      none
                      
!***********************************************************************
!***********************************************************************
!      EXECUTABLE CODE begins below
!***********************************************************************
!***********************************************************************
                      
!      Calculate longitude term
!      Add offset of 18 degrees, double it, convert to radians, and
!      take the cosine
                      COSLON=COS( PI180*2.0*( LON + 18.0 ) )
                      
!      Calculate the latitude terms
!      Convert Latitude into radians
                      LTRAD=PI180*LAT
!      Calculate sine and cosine terms
                      COSLT = COS(LTRAD)
                      COSLT2 = COSLT**2
                      SINLT2 = ( SIN(LTRAD ) )**2
                      SIN2LT = ( SIN( 2.0*LTRAD ) )**2
                      
!      Calculate the Earth's radius at this latitude
                      R = SQRT( B2/( 1.0 - COSLT2*ABTERM ) )
                      
!      Calculate total distance from Earth's center
                      RTOT = R + Z
                      
!      Calculate gravity at the Earth's surface
                      G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON )
                      
!      Calculate the centripetal term at the Earth's surface
!      Note: the centripetal acceleration due to Earth's rotation
!      is in a direction perpendicular to the Earth's rational axis,
!      but for this gravity function we are only interested in the
!      component parallel to the radial direction.
                      C_SUR = COSLT2*R*(PI2*W)**2
                      
!      Calculate the centripetal term at altitude z (with wind)
                      C_Z = ( ( PI2*RTOT*COSLT*W + WINDE )**2 + (WINDN)**2 )/RTOT
                      
!      Calculate the change in gravitation with altitude
                      GRAVZ=(G_SUR + C_SUR)*(1.0 - R**2/RTOT**2)
                      
                      GRAV_EARTH = G_SUR + (C_SUR - C_Z) - GRAVZ
                      
                      RETURN
                    END FUNCTION GRAV_EARTH
                    
!************************************************************************
! this reads in the reference levels profile
                    
                    SUBROUTINE ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iGasID,iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
                    
                    
                    NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
                    INTEGER, INTENT(IN OUT)                  :: iGasID
                    INTEGER, INTENT(IN OUT)                  :: iNumLevsx
                    NO TYPE, INTENT(IN OUT)                  :: raRx110Pre
                    NO TYPE, INTENT(IN OUT)                  :: raRx110Tem
                    REAL, INTENT(OUT)                        :: raRx110MR(kProfLayer+10)
                    IMPLICIT NONE
                    
                    INCLUDE '../INCLUDE/kcartaparam.f90'
                    
! input
                    
                    REAL :: PLEV_KCARTADATABASE_AIRS(kProfLayer+1)
! output
                    
                    REAL :: raRx110Temp(kProfLayer+10), raRx110Press(kProfLayer+10)
! local
                    REAL :: PLEVx110(kProfLayer+10)
                    REAL :: rLat,raJunk(kProfLayer+10)
                    CHARACTER (LEN=80) :: caPFname,caStr,caComment
                    INTEGER :: iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec
                    
! first extend PLEV_KCARTADATABASE_AIRS a little above its lwest value, so can do interps well
                    DO iI=1,kProfLayer+1
                      PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(iI)
                    END DO
                    rLat = PLEV_KCARTADATABASE_AIRS(kProfLayer+1)/20
                    DO iI=kProfLayer+2,kProfLayer+10
                      PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(kProfLayer+1) - rLat*(iI-(kProfLayer+1))
                    END DO
                    
!      caPFname = '/home/sergio/KCARTA/INCLUDE/glatm_16Aug2010.dat'
                    caPFname = kcaLevsRefProf
                    
                    IF ((kAFGLProf < 1) .OR. (kAFGLProf > 6)) THEN
                      WRITE(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
                      CALL DoStop
                    END IF
                    
                    iErr = 0
                    iIOUN2 = kProfileUnit
                    OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',  &
                        IOSTAT=iErrIO)
                    IF (iErrIO /= 0) THEN
                      iErr=1
                      WRITE(kStdErr,1070) iErrIO, caPfname
                      CALL DoSTOP
                    END IF
                    1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
                    kProfileUnitOpen=1
                    
                    iAFGL = 0
                    iFoundGas = -1
                    
                    10  CONTINUE
                    READ(iIOUN2,1080) caStr
                    IF (caStr(1:1) == '!') THEN
                      GO TO 10   !! keep reading comments
                    ELSE
                      READ (caStr,*) iNumLevsx
                    END IF
                    IF (iNumLevsx > kProfLayer+1) THEN
                      WRITE(kStdErr,*) 'oops iNumLevsx .GT. kProfLayer+1 ',iNumLevsx,kProfLayer+1
                      CALL DoStop
                    END IF
                    
                    20  CONTINUE
                    iAFGL = iAFGL + 1
                    
                    30  CONTINUE
                    READ(iIOUN2,1080) caStr
                    IF (caStr(1:1) == '!') THEN
                      GO TO 30   !! keep reading comments
                    ELSE
                      caComment = caStr
                    END IF
                    
                    40  CONTINUE
                    READ(iIOUN2,1080) caStr
                    IF (caStr(1:1) == '!') THEN
                      GO TO 40   !! keep reading comments
                    ELSE
                      READ (caStr,*) rLat
                    END IF
!      write(kStdWarn,*) 'iAFGL Profile ',iAFGL,' rLat ',rLat, ' : ',caComment
                    
                    READ(iIOUN2,1080) caStr   !! comment
                    READ(iIOUN2,1080) caStr   !! Altitude(km)
                    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)
                    
                    IF (iAFGL == kAFGLProf) THEN
!        write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! Press(mb)
                      READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)
                      
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! Temp(K)
                      READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)
                      
!        write(kStdWarn,*) '   need to find profile for gasID ',iGasID
                    ELSE
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! Press(mb)
                      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                      
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! Temp(K)
                      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                    END IF
                    
                    READ(iIOUN2,1080) caStr   !! comment
                    READ(iIOUN2,1080) caStr   !! density (cm-3)
                    READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)
                    
                    DO iG = 1,7
                      IF (iG == iGasID) THEN
!          write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
                        iFoundGas = +1
                        READ(iIOUN2,1080) caStr   !! comment
                        READ(iIOUN2,1080) caStr   !! gas name
                        READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
                      ELSE
                        READ(iIOUN2,1080) caStr   !! comment
                        READ(iIOUN2,1080) caStr   !! gas name
                        READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                      END IF
                    END DO
                    
                    IF ((iAFGL == kAFGLProf) .AND. (iFoundGas > 0)) THEN
                      GO TO 60   !! found the AFGL prof and found the gas; done!
                    ELSE IF ((iAFGL == 6) .AND. (iFoundGas < 0)) THEN
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! modend
                      GO TO 50   !! found the AFGL prof but not found the gas
                    ELSE IF ((iAFGL < 6) .OR. (iFoundGas < 0)) THEN
!! either did not find the gas or the AFGL prof
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! modend
                      GO TO 20
                    END IF
                    
                    50   CONTINUE
                    READ(iIOUN2,1080) caStr   !! comment
                    READ(iIOUN2,1080) caStr   !! constituent profs
                    READ(iIOUN2,1080) caStr   !! comment
                    READ(iIOUN2,1080) caStr   !! mingas
                    
                    iG = 7
                    55   CONTINUE
                    iG = iG + 1
                    IF (iG == iGasID) THEN
!        write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
                      iFoundGas = +1
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! gas name
                      READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
                    ELSE
                      READ(iIOUN2,1080) caStr   !! comment
                      READ(iIOUN2,1080) caStr   !! gas name
                      READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                    END IF
                    IF ((iG < 28) .AND. (iFoundGas < 0)) THEN
                      GO TO 55
                    END IF
                    
                    IF (iFoundGas < 0) THEN
!! bweh try xsec gases
                      123    CONTINUE
                      READ(iIOUN2,1080) caStr   !! comment
                      IF (caStr(1:1) == '!') THEN
                        GO TO 123
                      ELSE IF (caStr(1:6) == 'DATEND') THEN
                        GO TO 60       !!! oh oh end of file, give up!!!!!
                      ELSE
                        READ (caStr,*) iXsec
                        IF (iXsec == iGasID) THEN
                          iFoundGas = 1
                          READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
                        ELSE
                          READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                          GO TO 123
                        END IF
                      END IF
                    END IF
                    
                    60   CONTINUE
                    CLOSE(iIOUN2)
                    kProfileUnitOpen=-1
                    1080 FORMAT(A80)
                    
                    IF (iFoundGas < 0) THEN
!! finally give up
                      WRITE(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
                      CALL DOStop
                    END IF
                    
!*************************
! finally interp these onto the AIRS pressure levels
                    DO iI = 1,iNumLevsx
                      raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
!        print *,iI,raRx110Press(iI),raRx110Temp(iI),raRx110MR(iI)
                    END DO
                    
                    CALL r_sort_loglinear(raRx110Press,raRx110Temp,iNumLevsx,PLEVx110,raJunk,kProfLayer+10)
                    DO iI = 1,kProfLayer+10
                      raRx110Temp(iI) = raJunk(iI)
                    END DO
                    CALL r_sort_loglinear(raRx110Press,raRx110MR,iNumLevsx,PLEVx110,raJunk,kProfLayer+10)
                    DO iI = 1,kProfLayer+10
                      raRx110MR(iI) = raJunk(iI)
                      raRx110Press(iI) = PLEVx110(iI)
                    END DO
                    iNumLevsx = kProfLayer+10
                    
                    RETURN
                  END SUBROUTINE ReadRefProf_Levels
                  
!************************************************************************
! this subroutine does the klayers integrations, integrals over P
! so this is more similar to kLAYERs, CRTM
! this is newer code; it ignores surface as "lowest" level and only uses the user level info
! this is for AIRS 101 levels
                  
! >>>>>>>> raPout comes in mb
! >>>>>>>> and goes out in N/m2
                  
                  SUBROUTINE DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,  &
                      PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,  &
                      raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
                  
                  
                  REAL, INTENT(IN OUT)                     :: rHSurf
                  REAL, INTENT(IN OUT)                     :: rPSurf
                  REAL, INTENT(IN OUT)                     :: rTSurf
                  INTEGER, INTENT(IN OUT)                  :: iLowestLev
                  INTEGER, INTENT(IN)                      :: iNumGases
                  INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
                  REAL, INTENT(IN)                         :: rLat
                  REAL, INTENT(IN OUT)                     :: rLon
                  NO TYPE, INTENT(IN OUT)                  :: PAVG_KCART
                  NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
                  NO TYPE, INTENT(IN OUT)                  :: DATABASELE
                  REAL, INTENT(IN OUT)                     :: rfracBot
                  REAL, INTENT(IN)                         :: raPX(kProfLayer+1)
                  REAL, INTENT(IN)                         :: raTX(kProfLayer+1)
                  REAL, INTENT(IN OUT)                     :: raaG_VMRX(kProfLayer+1,kMaxGas)
                  REAL, INTENT(OUT)                        :: raPout(kProfLayer)
                  NO TYPE, INTENT(IN OUT)                  :: raAmountOu
                  REAL, INTENT(OUT)                        :: raTout(kProfLayer)
                  REAL, INTENT(OUT)                        :: raZout(kProfLayer+1)
                  REAL, INTENT(OUT)                        :: raaQout(kProfLayer,kMaxGas)
                  NO TYPE, INTENT(IN OUT)                  :: raaPartPre
                  IMPLICIT NONE
                  
                  INCLUDE '../INCLUDE/kcartaparam.f90'
                  
! input
                  REAL :: PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_
                  REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)
                  
                  
! output
                  REAL :: raAmountOut(kProfLayer)
                  REAL :: raaPartPressOut(kProfLayer,kMaxGas)
                  
! local
                  INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
                  REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
                  REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
                  REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
                  REAL :: raXYZPress(kMaxlayer+1+1)
                  REAL :: raXYZTemp(kMaxlayer+1+1)
                  REAL :: raaXYZ_VMR(kMaxlayer+1+1,kMaxGas)
                  REAL :: raXYZ_VMRwater(kMaxlayer+1+1)
                  REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
                  REAL :: zWoo,rConvertQ,grav_earth,wexsvp,rRH
                  INTEGER :: iWoo
                  
                  rConvertQ = 6.022141E23/1E4     ! multiply by this to change moles/m2 to molecules/cm2
                  rConvertQ = kAvog/1000/1.0E4    ! multiply by this to change moles/m2 to molecules/cm2
                  rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor
                  
                  DO iL = 1,kProfLayer
                    DO iG = 1,iNumGases
                      raaQout(iL,iG) = 0.0
                    END DO
                  END DO
                  
                  IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
                    WRITE (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
                    CALL DoStop
                  END IF
                  
                  iNFine = 10
                  iNFine = 200
                  iNFine = 50  !! default
                  
                  WRITE(kStdWarn,*) '  '
                  WRITE(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers integration >>>>>>>>>>>>>>>>>>>'
                  
! >>>>>>>>>>
!! set up temporary arrays/matrices for interping onto sublevels
                  iUpperLev = kMaxLayer-iLowestLev+1
!! set levels 1 .. iNumLevs
                  DO iL = 1,iUpperLev
                    raXYZPress(iL) = raPX(iL)
                    raXYZTemp(iL) = raTX(iL)
                    DO iG = 1,iNumGases
                      raaXYZ_VMR(iL,iG) = raaG_VMRX(iL,iG)
                      IF (iG ==  1) raXYZ_VMRwater(iL) = raaG_VMRX(iL,iG)
                    END DO
                  END DO
! >>>>>>>>>>
                  
! now start integrating
! lowest layer is special case
                  
! http://cimss.ssec.wisc.edu/~paulv/
! http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
!     /home/sergio/KCARTA/DOC/sci_klayers.txt
                  
                  1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
                  
                  z = rHSurf !! need to kludge this
                  zWoo = rHSurf
                  iWoo = -1
                  
                  DO iL = iLowestLev,iLowestLev
                    raPout(iL) = LOG(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
                    raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
                    raPout(iL) = raPout(iL)/100.0    !!! change N/m2 to mb
                    
                    z0 = z
                    iJ = iL - iLowestLev + 1
                    
                    amount = 0.0
                    rTsum  = 0.0
                    rWgt = 0.0
                    
!!! >>>>>>  this is starting out at         rP = PLEV_KCARTADATABASE_AIRS(iL)        rT = rTSurfx  >>>>>>>>
                    rP = PLEV_KCARTADATABASE_AIRS(iL)
                    rP = rP/100.0    !!! change N/m2 to mb
                    CALL r_sort_loglinear1(raXYZPress,raXYZTemp,iUpperLev,rP,rT,1)
                    
!! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = log(a/c)
                    dlnp = LOG(PLEV_KCARTADATABASE_AIRS(iL)) - LOG(PLEV_KCARTADATABASE_AIRS(iL+1))
                    dlnp = dlnp / (iNFine)
                    
!! information for current (sub)level
                    rP_n = rP                             !! current sublev press
                    rT_n = rT                             !! current sublev temp
                    rMR_water_n = raXYZ_VMRwater(1)       !! current water MR
                    CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP,rMR_water_n,1) !! current water MR
                    
!! information for next (sub)level
                    rP_np1 = LOG(rP_n) - dlnp
                    rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                    CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                    CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)  !! next sublev MRw
                    
                    IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
                      iWoo = +1
                      zWoo   = z
                    END IF
                    
                    DO iLoop = 1,iNFine
                      
                      rMolarMass_n   = kAtmMolarMass
                      rMolarMass_np1 = kAtmMolarMass
!!! if Earth, adjust molecular mass for displacement by water
                      IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
                      IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
                      rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
                      rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol
                      
                      density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
                      density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC
                      
                      IF (kPlanet /= 3) THEN
                        gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
                        gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
                      ELSE
                        gravity = grav_earth(z,0.0,0.0,rLat,rLon)
                      END IF
                      
                      Tav   = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
                      SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
                      dz = SpHeatDryAir * Tav / gravity * LOG(rP_n/rP_np1)
                      
                      rTsum = rTsum + Tav * (density_n + density_np1)/2
                      rWgt  = rWgt + (density_n + density_np1)/2
                      
                      Pav = (rP_n-rP_np1)/LOG(rP_n/rP_np1)
                      damount = Pav/Tav/kMGC * dz
!          print *,iLoop,PLEV_KCARTADATABASE_AIRS(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
                      amount = amount + damount
                      
                      DO iG = 1,iNumGases
                        DO iCnt = 1,iUpperLev
                          raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
                        END DO
                        CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1)
                        CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1)
                        raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
                      END DO
                      
!!! now update for next iteration
                      z = z + dz
                      rP_n = rP_np1
                      rT_n = rT_np1
                      rMR_water_n = rMR_water_np1
                      rP = rP_n
                      
                      rP_np1 = LOG(rP_n) - dlnp
                      rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                      CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                      CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                      
                      IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
                        iWoo = +1
                        zWoo   = z
                      END IF
                      
                    END DO
                    
                    raTout(iL)      = rTsum/rWgt                 ! <<< need TAV
                    raAmountOut(iL) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
                    raZout(iL+1)    = z
                    
!        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
!     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
                  END DO
                  
! displace everything  z
                  iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
                  raZout(iL)   = -(zWoo - z0)
                  raZout(iL+1) = +(z-zWoo)
                  WRITE(kStdWarn,*) '  need to displace lowest layer heights'
                  WRITE(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
                  WRITE(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',PLEV_KCARTADATABASE_AIRS(iL),rPSurf,  &
                      PLEV_KCARTADATABASE_AIRS(iL+1)
                  z = z - zWoo
                  
!cc      !no need to adjust this amount for book-keeping, make consistent with eg n_rad_jac_scat.f since we did FULL layer
!cc      raAmountOut(iL) = raAmountOut(iL)/rFracBot
!cc      DO iG = 1,iNumGases
!cc        raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
!cc      END DO
                  
! go to TOA
                  DO iL = iLowestLev + 1,kProfLayer
                    z0 = z
                    
!! compute avg pressure of layer
                    raPout(iL) = LOG(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
                    raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
                    raPout(iL) = raPout(iL)/100.0    !!! change N/m2 to mb
                    
                    iJ = iL - iLowestLev + 1
                    rP = PLEV_KCARTADATABASE_AIRS(iL)
                    rP = rP/100.0    !!! change N/m2 to mb
                    rT = raXYZTemp(iJ)
                    
                    amount = 0.0
                    rTsum  = 0.0
                    rWgt = 0.0
                    
!! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
                    dlnp = LOG(PLEV_KCARTADATABASE_AIRS(iL)) - LOG(PLEV_KCARTADATABASE_AIRS(iL+1))
                    dlnp = dlnp / (iNFine)
                    
!! information for current (sub)level
                    rP_n = rP                              !! current sublev press
                    rT_n = rT                              !! current sublev temp
                    rMR_water_n = raXYZ_VMRwater(iJ)        !! current water MR
                    
!! information for next (sub)level
                    rP_np1 = LOG(rP_n) - dlnp
                    rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                    CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                    CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                    
                    DO iLoop = 1,iNFine
                      
                      rMolarMass_n   = kAtmMolarMass
                      rMolarMass_np1 = kAtmMolarMass
!!! if Earth, adjust molecular mass for displacement by water
                      IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
                      IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
                      rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
                      rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol
                      
                      density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
                      density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC
                      
                      IF (kPlanet /= 3) THEN
                        gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
                        gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
                      ELSE
                        gravity = grav_earth(z,0.0,0.0,rLat,rLon)
                      END IF
                      
                      Tav = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
                      SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
                      dz = SpHeatDryAir * Tav / gravity * LOG(rP_n/rP_np1)
                      
                      rTsum = rTsum + Tav * (density_n + density_np1)/2
                      rWgt  = rWgt + (density_n + density_np1)/2
                      
                      Pav = (rP_n-rP_np1)/LOG(rP_n/rP_np1)
                      damount = Pav/Tav/kMGC * dz
                      amount = amount + damount
                      
                      DO iG = 1,iNumGases
                        DO iCnt = 1,iUpperLev
                          raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
                        END DO
                        CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1)
                        CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1)
                        raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
!            q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
!            q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
!            raaQout(iL,iG) = raaQout(iL,iG) + (q_n + q_np1)/2 * rConvertQ
                      END DO
                      
!!! now update for next iteration
                      z = z + dz
                      rP_n = rP_np1
                      rT_n = rT_np1
                      rMR_water_n = rMR_water_np1
                      rP = rP_n
                      
                      rP_np1 = LOG(rP_n) - dlnp
                      rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                      CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                      CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                      
                    END DO
                    
                    raTout(iL)      = rTsum/rWgt             ! <<< need TAV
                    raAmountOut(iL) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
                    raZout(iL+1)    = z
                    
!        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
!     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
                  END DO
                  
                  WRITE(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
!! raZout is in meters
!! raaQout is in moles/m2
                  DO iL = iLowestLev,kProfLayer
                    rThick = ABS(raZout(iL) - raZout(iL+1)) * 100.0
                    rT = raTout(iL)
                    raPout(iL) = raPout(iL) * 100.0    !! change mb to N/m2
                    DO iG = 1,iNumGases
                      raaPartPressOut(iL,iG) = raaQout(iL,iG) / raAmountOut(iL) * raPout(iL)
                    END DO
                  END DO
                  
! testing
!          !!! this automatically puts partial pressure in ATM, assuming
!          !!! gas amount in kilomoles/cm2, length in cm, T in kelvin
!          !!!note "j"!!!
!      DO iL = iLowestLev,kProfLayer
!        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
!        rT = raTout(iL)
!        DO iG = 1,iNumGases
!          rAmt = raaQout(iL,iG)/kAvog
!          rPP  = rAmt*1.0e9*kMGC*rT / (rThick*kAtm2mb*100.0)
!          IF (iG .EQ. 1) THEN
!            print *,iL,iG,raaPartPressOut(iL,iG)/(kAtm2mb*100).0,rPP,rThick,
!     $             raaPartPressOut(iL,iG)/rPP/(kAtm2mb*100.0),raAmountOut(iL),raaQout(iL,iG)
!          END IF
!        print *,iL,raPout(iL)/100,raTout(iL),raAmountOut(iL),
!     $        raaQout(iL,1),raaQout(iL,2),raaQout(iL,3),raaQout(iL,4),raaQout(iL,5),
!     $        raaG_VMRX(iL-iLowestLev+1,1),raaG_VMRX(iL-iLowestLev+1,2),raaG_VMRX(iL-iLowestLev+1,3),
!     $        raaG_VMRX(iL-iLowestLev+1,4),raaG_VMRX(iL-iLowestLev+1,5)
!        END DO
!      END DO
! testing
                  
                  IF ((iLowestLev-1) >= 1) THEN
                    DO iL = 1,iLowestLev-1
                      raPout(iL) = LOG(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
                      raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
                      raTout(iL) = kStempMin
                      DO iG = 1,iNumGases
                        raaPartPressOut(iL,iG) = 0.0
                        raaQout(iL,iG) = 0.0
                        raaG_VMRX(iL,iG) = 0.0
                      END DO
                    END DO
                  END IF
                  
                  IF (kPlanet == 3) THEN
!! show the RH for gas 1
                    WRITE(kStdWarn,*) ' '
                    WRITE(kStdWarn,*) 'Checking Relative Humidity (integrate over AIRS 101 levels)'
                    WRITE(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Q6'
                    WRITE(kStdWarn,113) '------------------------------------------------------------------------------------------'
                    DO iL = iLowestLev,kProfLayer
                      z    = (raZout(iL) + raZout(iL+1))/2/1000
                      zWoo = (raZout(iL+1) - raZout(iL))/1000
                      rP   = raPout(iL)            !! N/m2
                      rT   = raTout(iL)
                      rPP  = raaPartPressOut(iL,1) !! N/m2
                      rSVP = wexsvp(rT) * 100      !! mb --> N/m2
                      rRH  = rPP/rSVP*100.0
                      IF (rRH <= 100.0) THEN
                        WRITE(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),  &
                            (raaQout(iL,iG),iG=1,6)
                      ELSE
                        WRITE(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),  &
                            (raaQout(iL,iG),iG=1,6),' ***** '
                      END IF
                    END DO
                  END IF
                  111  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '))
                  112  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '),A7)
                  113  FORMAT(A90)
                  
! now find pressure output corresponding to HGT output from LBLRTM
                  IF (kRTP == -20) THEN
                    IF (raRTP_TxtInput(6) > 0) THEN
!!! input level boundaries in km, change to mb
                      DO iL = iLowestLev,kProfLayer
                        raJunk(iL-iLowestLev+1) = raZout(iL)
                        raJunk2(iL-iLowestLev+1) = raPout(iL)
                      END DO
                      CALL r_sort_linear1(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo,1)
                      WRITE(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
                      raRTP_TxtInput(6) = zWoo/100.0  !! mb
                    ELSE IF (raRTP_TxtInput(6) < 0) THEN
!!! input level boundaries in mb
                      raRTP_TxtInput(6) = ABS(raRTP_TxtInput(6)) * 100.0
                      WRITE(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
                    END IF
                  END IF
                  
                  WRITE(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
                  
                  RETURN
                END SUBROUTINE DoIntegrateLevels2Layers_wrtP
                
!************************************************************************
! this subroutine does the klayers integrations, integrals over P
! so this is more similar to kLAYERs, CRTM
! this is newer code; it ignores surface as "lowest" level and only uses the user level info
! this is for USER DEFINED levels
                
! >>>>>>>> raPout comes in mb
! >>>>>>>> and goes out in N/m2
                
                SUBROUTINE DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,  &
                    raPBnd,raAlt,iZbnd,rfracBot,  &
                    raPX,raTX,raaG_VMRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
                
                
                REAL, INTENT(IN OUT)                     :: rHSurf
                REAL, INTENT(IN OUT)                     :: rPSurf
                REAL, INTENT(IN OUT)                     :: rTSurf
                INTEGER, INTENT(IN OUT)                  :: iLowestLev
                INTEGER, INTENT(IN)                      :: iNumGases
                INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
                REAL, INTENT(IN)                         :: rLat
                REAL, INTENT(IN OUT)                     :: rLon
                REAL, INTENT(IN)                         :: raPBnd(2*kProfLayer)
                REAL, INTENT(IN OUT)                     :: raAlt(2*kProfLayer)
                INTEGER, INTENT(IN)                      :: iZbnd
                REAL, INTENT(IN OUT)                     :: rfracBot
                REAL, INTENT(IN)                         :: raPX(kProfLayer+1)
                REAL, INTENT(IN)                         :: raTX(kProfLayer+1)
                REAL, INTENT(IN OUT)                     :: raaG_VMRX(kProfLayer+1,kMaxGas)
                REAL, INTENT(OUT)                        :: raPout(kProfLayer)
                NO TYPE, INTENT(IN OUT)                  :: raAmountOu
                REAL, INTENT(OUT)                        :: raTout(kProfLayer)
                REAL, INTENT(OUT)                        :: raZout(kProfLayer+1)
                REAL, INTENT(OUT)                        :: raaQout(kProfLayer,kMaxGas)
                NO TYPE, INTENT(IN OUT)                  :: raaPartPre
                IMPLICIT NONE
                
                INCLUDE '../INCLUDE/kcartaparam.f90'
                
! input
                
                
                
                
! output
                REAL :: raAmountOut(kProfLayer)
                REAL :: raaPartPressOut(kProfLayer,kMaxGas)
                
! local
                INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
                REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
                REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
                REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
                REAL :: raXYZPress(kMaxlayer+1+1)
                REAL :: raXYZTemp(kMaxlayer+1+1)
                REAL :: raaXYZ_VMR(kMaxlayer+1+1,kMaxGas)
                REAL :: raXYZ_VMRwater(kMaxlayer+1+1)
                REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
                REAL :: zWoo,rConvertQ,grav_earth,wexsvp,rRH
                INTEGER :: iOffset,iWoo
                
                rConvertQ = 6.022141E23/1E4     ! multiply by this to change moles/m2 to molecules/cm2
                rConvertQ = kAvog/1000/1.0E4    ! multiply by this to change moles/m2 to molecules/cm2
                rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor
                
                DO iL = 1,kProfLayer
                  DO iG = 1,iNumGases
                    raaQout(iL,iG) = 0.0
                  END DO
                END DO
                
                IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
                  WRITE (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
                  CALL DoStop
                END IF
                
                iNFine = 10
                iNFine = 200
                iNFine = 50
                
                WRITE(kStdWarn,*) '  '
                WRITE(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
                
! >>>>>>>>>>
!! set up temporary arrays/matrices for interping onto sublevels
                iUpperLev = iZbnd
!! set levels 1 .. iNumLevs
                DO iL = 1,iUpperLev
                  raXYZPress(iL) = raPX(iL)
                  raXYZTemp(iL) = raTX(iL)
                  DO iG = 1,iNumGases
                    raaXYZ_VMR(iL,iG) = raaG_VMRX(iL,iG)
                    IF (iG ==  1) raXYZ_VMRwater(iL) = raaG_VMRX(iL,iG)
                  END DO
!      print *,iNumGases,iL,raPX(iL),raTX(iL),raaG_VMRX(iL,1),raaG_VMRX(iL,7),raaG_VMRX(iL,22)
                END DO
                
! >>>>>>>>>>
                
! now start integrating
! lowest layer is special case
                
! http://cimss.ssec.wisc.edu/~paulv/
! http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
!     /home/sergio/KCARTA/DOC/sci_klayers.txt
                
                1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
                
                z = rHSurf !! need to kludge this
                zWoo = rHSurf
                iWoo = -1
                
!      DO iL=1,iZbnd
!        print *,'sergio1',iL,raPBnd(iL)
!      end do
                
                iOffSet = kProfLayer-(iZbnd-1)
                DO iL = iLowestLev,iLowestLev
                  raPout(iL+iOffset) = LOG(raPBnd(iL)/raPBnd(iL+1))
                  raPout(iL+iOffSet) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL+iOffset)
                  raPout(iL+iOffSet) = raPout(iL+iOffSet)/100.0    !!! change N/m2 to mb
                  
                  z0 = z
                  iJ = iL - iLowestLev + 1
                  
                  amount = 0.0
                  rTsum  = 0.0
                  rWgt = 0.0
                  
!!! >>>>>>  this is starting out at         rP = raPBnd(iL)        rT = rTSurfx  >>>>>>>>
                  rP = raPBnd(iL)
                  rP = rP/100.0    !!! change N/m2 to mb
                  CALL r_sort_loglinear1(raXYZPress,raXYZTemp,iUpperLev,rP,rT,1)
                  
!! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
                  dlnp = LOG(raPBnd(iL)) - LOG(raPBnd(iL+1))
                  dlnp = dlnp / (iNFine)
                  
!! information for current (sub)level
                  rP_n = rP                             !! current sublev press
                  rT_n = rT                             !! current sublev temp
                  rMR_water_n = raXYZ_VMRwater(1)        !! current water MR
                  CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP,rMR_water_n,1) !! current water MR
                  
!! information for next (sub)level
                  rP_np1 = LOG(rP_n) - dlnp
                  rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                  CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                  CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                  
                  IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
                    iWoo = +1
                    zWoo   = z
                  END IF
                  
                  DO iLoop = 1,iNFine
                    
                    rMolarMass_n   = kAtmMolarMass
                    rMolarMass_np1 = kAtmMolarMass
!!! if Earth, adjust molecular mass for displacement by water
                    IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
                    IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
                    rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
                    rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol
                    
                    density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
                    density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC
                    
                    IF (kPlanet /= 3) THEN
                      gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
                      gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
                    ELSE
                      gravity = grav_earth(z,0.0,0.0,rLat,rLon)
                    END IF
                    
                    Tav   = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
                    SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
                    dz = SpHeatDryAir * Tav / gravity * LOG(rP_n/rP_np1)
                    
                    rTsum = rTsum + Tav * (density_n + density_np1)/2
                    rWgt  = rWgt + (density_n + density_np1)/2
                    
                    Pav = (rP_n-rP_np1)/LOG(rP_n/rP_np1)
                    damount = Pav/Tav/kMGC * dz
!          print *,iLoop,raPBnd(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
                    amount = amount + damount
                    
                    DO iG = 1,iNumGases
                      DO iCnt = 1,iUpperLev
                        raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
                      END DO
                      CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1)
                      CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1)
                      raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
                    END DO
                    
!!! now update for next iteration
                    z = z + dz
                    rP_n = rP_np1
                    rT_n = rT_np1
                    rMR_water_n = rMR_water_np1
                    rP = rP_n
                    
                    rP_np1 = LOG(rP_n) - dlnp
                    rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                    CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                    CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                    
                    IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
                      iWoo = +1
                      zWoo   = z
                    END IF
                    
                  END DO
                  
                  raTout(iL+iOffSet)      = rTsum/rWgt                 ! <<< need TAV
                  raAmountOut(iL+iOffSet) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
                  raZout(iL+1+iOffSet)    = z
                  
!        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL+iOffSet),
!     $        raAmountOut(iL+iOffSet),raZout(iL+1+iOffSet),raPout(iL+iOffSet),z/1000
                END DO
                
! displace everything  z
                iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
                raZout(iL+iOffSet)   = -(zWoo - z0)
                raZout(iL+1+iOffSet) = +(z-zWoo)
                WRITE(kStdWarn,*) '  need to displace lowest layer heights'
                WRITE(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
                WRITE(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',raPBnd(iL),rPSurf,  &
                    raPBnd(iL+1)
                z = z - zWoo
                
!cc      !no need to adjust this amount for book-keeping, make consistent with eg n_rad_jac_scat.f since we did FULL layer
!cc      raAmountOut(iL) = raAmountOut(iL)/rFracBot
!cc      DO iG = 1,iNumGases
!cc        raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
!cc      END DO
                
! go to TOA as defined by user
                
                DO iL = iLowestLev + 1,iZbnd-1
                  z0 = z
!! compute avg pressure of layer
                  raPout(iL+iOffSet) = LOG(raPBnd(iL)/raPBnd(iL+1))
                  raPout(iL+iOffSet) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL+iOffset)
                  raPout(iL+iOffSet) = raPout(iL+iOffSet)/100.0    !!! change N/m2 to mb
                  
                  iJ = iL - iLowestLev + 1
                  rP = raPBnd(iL)
                  rP = rP/100.0    !!! change N/m2 to mb
                  rT = raXYZTemp(iJ)
                  
                  amount = 0.0
                  rTsum  = 0.0
                  rWgt = 0.0
                  
!! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
                  dlnp = LOG(raPBnd(iL)) - LOG(raPBnd(iL+1))
                  dlnp = dlnp / (iNFine)
                  
!! information for current (sub)level
                  rP_n = rP                              !! current sublev press
                  rT_n = rT                              !! current sublev temp
                  rMR_water_n = raXYZ_VMRwater(iJ)        !! current water MR
                  
!! information for next (sub)level
                  rP_np1 = LOG(rP_n) - dlnp
                  rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                  CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                  CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                  
                  DO iLoop = 1,iNFine
                    
                    rMolarMass_n   = kAtmMolarMass
                    rMolarMass_np1 = kAtmMolarMass
!!! if Earth, adjust molecular mass for displacement by water
                    IF (kPlanet == 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
                    IF (kPlanet == 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
                    rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
                    rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol
                    
                    density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
                    density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC
                    
                    IF (kPlanet /= 3) THEN
                      gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
                      gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
                    ELSE
                      gravity = grav_earth(z,0.0,0.0,rLat,rLon)
                    END IF
                    
                    Tav = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)
                    SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)
                    dz = SpHeatDryAir * Tav / gravity * LOG(rP_n/rP_np1)
                    
                    rTsum = rTsum + Tav * (density_n + density_np1)/2
                    rWgt  = rWgt + (density_n + density_np1)/2
                    
                    Pav = (rP_n-rP_np1)/LOG(rP_n/rP_np1)
                    damount = Pav/Tav/kMGC * dz
                    amount = amount + damount
                    
                    DO iG = 1,iNumGases
                      DO iCnt = 1,iUpperLev
                        raJunk(iCnt) = raaXYZ_VMR(iCnt,iG)
                      END DO
                      CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1)
                      CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1)
                      raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
!            q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
!            q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
!            raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + (q_n + q_np1)/2 * rConvertQ
                    END DO
                    
!!! now update for next iteration
                    z = z + dz
                    rP_n = rP_np1
                    rT_n = rT_np1
                    rMR_water_n = rMR_water_np1
                    rP = rP_n
                    
                    rP_np1 = LOG(rP_n) - dlnp
                    rP_np1 = EXP(rP_np1)                                                                !! next sublev pressure
                    CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
                    CALL r_sort_loglinear1(raXYZPress,raXYZ_VMRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw
                    
                  END DO
                  
                  raTout(iL+iOffSet)      = rTsum/rWgt             ! <<< need TAV
                  raAmountOut(iL+iOffSet) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
                  raZout(iL+1+iOffSet)    = z
                  
!        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL),
!     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
                END DO
                
                WRITE(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
                DO iL = 1,iOffSet-1
                  raPout(iL) = 0.0
                  raTout(iL) = 0.0
                  raZout(iL+1) = raZout(iZbnd)
                  raAmountOut(iL) = 0.0
                END DO
                
!! raZout is in meters
!! raaQout is in moles/m2
                DO iL = iLowestLev,iZbnd-1
                  rThick = ABS(raZout(iL+iOffSet) - raZout(iL+1+iOffSet)) * 100.0
                  rT = raTout(iL+iOffSet)
                  raPout(iL+iOffSet) = raPout(iL+iOffSet) * 100.0    !! change mb to N/m2
                  DO iG = 1,iNumGases
                    raaPartPressOut(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) / raAmountOut(iL+iOffSet) * raPout(iL+iOffSet)
                  END DO
                END DO
                
! testing
!          !!! this automatically puts partial pressure in ATM, assuming
!          !!! gas amount in kilomoles/cm2, length in cm, T in kelvin
!          !!!note "j"!!!
!      DO iL = iLowestLev,kProfLayer
!        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
!        rT = raTout(iL)
!        DO iG = 1,iNumGases
!          rAmt = raaQout(iL,iG)/kAvog
!          rPP  = rAmt*1.0e9*kMGC*rT / (rThick*kAtm2mb*100.0)
!          IF (iG .EQ. 1) THEN
!            print *,iL,iG,raaPartPressOut(iL,iG)/(kAtm2mb*100).0,rPP,rThick,
!     $             raaPartPressOut(iL,iG)/rPP/(kAtm2mb*100.0),raAmountOut(iL),raaQout(iL,iG)
!          END IF
!        print *,iL,raPout(iL)/100,raTout(iL),raAmountOut(iL),
!     $        raaQout(iL,1),raaQout(iL,2),raaQout(iL,3),raaQout(iL,4),raaQout(iL,5),
!     $        raaG_VMRX(iL-iLowestLev+1,1),raaG_VMRX(iL-iLowestLev+1,2),raaG_VMRX(iL-iLowestLev+1,3),
!     $        raaG_VMRX(iL-iLowestLev+1,4),raaG_VMRX(iL-iLowestLev+1,5)
!        END DO
!      END DO
! testing
                
!      IF ((iLowestLev-1) .GE. 1) THEN
                DO iL = 1,iOffSet-1
                  raPout(iL) = LOG(raPBnd(iL)/raPBnd(iL+1))
                  raPout(iL) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL)
                  raTout(iL) = kStempMin
                  DO iG = 1,iNumGases
                    raaPartPressOut(iL,iG) = 0.0
                    raaQout(iL,iG) = 0.0
                    raaG_VMRX(iL,iG) = 0.0
                  END DO
                END DO
!      END IF
                
                IF (iZbnd+iOffSet <= kProfLayer) THEN
                  DO iL = iZbnd,kProfLayer
                    raPout(iL) = 0.0
                    raTout(iL) = kStempMin
                    DO iG = 1,iNumGases
                      raaPartPressOut(iL,iG) = 0.0
                      raaQout(iL,iG) = 0.0
                      raaG_VMRX(iL,iG) = 0.0
                    END DO
                  END DO
                END IF
                
                IF (kPlanet == 3) THEN
!! show the RH for gas 1
                  WRITE(kStdWarn,*) ' '
                  WRITE(kStdWarn,*) 'Checking Relative Humidity (integrate over user levels)'
                  WRITE(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
                  WRITE(kStdWarn,113) '------------------------------------------------------------------------------------------'
                  DO iL = iLowestLev,iZbnd-1
                    z    = (raZout(iL+iOffset) + raZout(iL+1+iOffset))/2/1000
                    zWoo = (raZout(iL+1+iOffset) - raZout(iL+iOffset))/1000
                    rP   = raPout(iL+iOffset)            !! N/m2
                    rT   = raTout(iL+iOffset)
                    rPP  = raaPartPressOut(iL+iOffset,1) !! N/m2
                    rSVP = wexsvp(rT) * 100      !! mb --> N/m2
                    rRH  = rPP/rSVP*100.0
                    IF (rRH <= 100.0) THEN
                      WRITE(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet),  &
                          (raaQout(iL+iOffSet,iG),iG=1,6)
                    ELSE
                      WRITE(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet),  &
                          (raaQout(iL+iOffSet,iG),iG=1,6),' BAD RH!'
                    END IF
                  END DO
                END IF
                111  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '))
                112  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '),A7)
                113  FORMAT(A90)
                
! now find pressure output corresponding to HGT output from LBLRTM
                IF (kRTP == -20) THEN
                  IF (raRTP_TxtInput(6) > 0) THEN
!!! input level boundaries in km, change to mb
                    DO iL = iLowestLev,kProfLayer
                      raJunk(iL-iLowestLev+1) = raZout(iL)
                      raJunk2(iL-iLowestLev+1) = raPout(iL)
                    END DO
                    CALL r_sort_linear1(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo,1)
                    WRITE(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
                    raRTP_TxtInput(6) = zWoo/100.0  !! mb
                  ELSE IF (raRTP_TxtInput(6) < 0) THEN
!!! input level boundaries in mb
                    raRTP_TxtInput(6) = ABS(raRTP_TxtInput(6)) * 100.0
                    WRITE(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
                  END IF
                END IF
                
                WRITE(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
                
                RETURN
              END SUBROUTINE DoIntegrateLevels2UserLayers_wrtP
              
!************************************************************************
! >>> read in user supplied prof
! this will be at whatever gas units eg RH g/kg  VMR etc, but then CONVERTED to VMR (gasunit 12)
              
              SUBROUTINE  ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
                  iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,raP,raT,  &
                  iaG,iaGasUnits,raaG_VMR,rPMin,rPMax,rYear,rLat,rLon)
              
              
              NO TYPE, INTENT(IN OUT)                  :: caPFname
              NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
              NO TYPE, INTENT(IN OUT)                  :: rHminKCart
              NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
              NO TYPE, INTENT(IN OUT)                  :: rPminKCart
              INTEGER, INTENT(IN OUT)                  :: iNumLevs
              REAL, INTENT(OUT)                        :: rPSurf
              REAL, INTENT(OUT)                        :: rTSurf
              REAL, INTENT(OUT)                        :: rHSurf
              INTEGER, INTENT(IN)                      :: iNumGases
              REAL, INTENT(OUT)                        :: raP(2*kProfLayer)
              REAL, INTENT(OUT)                        :: raT(2*kProfLayer)
              INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
              INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
              REAL, INTENT(OUT)                        :: raaG_VMR(2*kProfLayer,kMaxGas)
              NO TYPE, INTENT(IN OUT)                  :: rPMin
              NO TYPE, INTENT(IN OUT)                  :: rPMax
              REAL, INTENT(IN OUT)                     :: rYear
              REAL, INTENT(IN OUT)                     :: rLat
              REAL, INTENT(IN OUT)                     :: rLon
              IMPLICIT NONE
              
              INCLUDE '../INCLUDE/kcartaparam.f90'
              
! input
              CHARACTER (LEN=80) :: caPfName
              REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
! output
              
              REAL :: rPmin,rPmax, tAugment
              
              
! local var
              INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iReadInOK,iIgnore,iYes
              REAL :: raX(kMaxGas),rX,rP,rT,raP_Nm2(2*kProfLayer),rPRefMin,rPRefMax,rPRefMinAugmented
              CHARACTER (LEN=80) :: caStr
              CHARACTER (LEN=40) :: caaUnit(50)
              
              DO iL = 1,50
                caaUnit(iL) = '** Not Known **'
              END DO
              caaUnit(01) = 'layer amt molecules/cm2'
              caaUnit(10) = 'vmr in ppmv'
              caaUnit(11) = 'vmr in ppbv'
              caaUnit(20) = 'mmr in g/kg'
              caaUnit(21) = 'mmr in g/g'
              caaUnit(30) = 'part press in mb (not accepted)'
              caaUnit(31) = 'part press in atm (not accepted)'
              caaUnit(40) = 'RH % (not accepted)'
              caaUnit(41) = 'RH fraction (not accepted)'
              
              rPRefMax = 1100.0 !! 1100  mb is AIRS RefPress LowerBdry      see SUBR ReadRefProf_Levels
              rPRefMin = 3.0E-3 !! 0.005 mb extended to 0.003 mb = 0.3 N/m2 see SUBR ReadRefProf_Levels
              
              rPmin = +1.0E6
              rPmax = -1.0E+6
              
              WRITE(kStdWarn,*) 'Reading in user supplied TXT LVLS file .....'
              
              iIgnore = -1
              iIOUN2 = kProfileUnit
              OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',  &
                  IOSTAT=iErrIO)
              IF (iErrIO /= 0) THEN
                iErr=1
                WRITE(kStdErr,1070) iErrIO, caPfname
                1070     FORMAT('ERROR! number ',I5,' opening PRFILE path LEVELS profile file'  &
                    ,/,A80)
                CALL DoSTOP
              END IF
              kProfileUnitOpen=1
              
              READ (iIOUN2,5030,ERR=13,END=13) caStr
              
              READ (iIOUN2,*) iNumLevs
              IF (iNumLevs > 2*kProfLayer) THEN
                WRITE(kStdErr,*) 'iNumLevs .GT. 2*kProfLayer',iNumLevs,kProfLayer
                CALL DoStop
              END IF
              
              READ (iIOUN2,*) rPSurf,rTSurf,rHSurf
              IF ((rHSurf > rHmaxKCarta) .OR. (rHSurf < rHminKCarta)) THEN
                WRITE(kStdErr,*) 'need rHmaxKCarta >= rHSurf >= rHminKCarta but have'
                WRITE(kStdErr,*) '(rHmaxKCarta,rHSurf,rHminKCarta) = ',rHmaxKCarta,rHSurf,rHminKCarta
                CALL DoStop
              END IF
              IF ((rPSurf > rPmaxKCarta) .OR. (rPSurf < rPminKCarta)) THEN
                WRITE(kStdErr,*) 'need rPmaxKCarta >= rPSurf >= rPminKCarta but have'
                WRITE(kStdErr,*) '(rPmaxKCarta,rPSurf,rPminKCarta) = ',rPmaxKCarta,rPSurf,rPminKCarta
                CALL DoStop
              END IF
              IF ((rTSurf > kStempMax) .OR. (rTSurf < kStempMin)) THEN
                WRITE(kStdErr,*) 'need kStempMax >= rTSurf >= kStempMin but have'
                WRITE(kStdErr,*) '(kStempMax,rTSurf,kStempMin) = ',kStempMax,rTSurf,kStempMin
                CALL DoStop
              END IF
              
              raRTP_TxtInput(1) = rPSurf
              raRTP_TxtInput(2) = rTSurf
              raRTP_TxtInput(3) = rHSurf
              
              READ (iIOUN2,*) rYear,rLat,rLon
              
              READ (iIOUN2,*) iNumGases
              IF (iNumGases > kMaxGas) THEN
                WRITE(kStdErr,*) 'iNumGases .GT. kMaxGas',iNumGases,kMaxGas
                CALL DoStop
              END IF
              
              READ (iIOUN2,*) (iaG(iJ),iJ=1,iNumGases)
              READ (iIOUN2,*) (iaGasUnits(iJ),iJ=1,iNumGases)
              
              WRITE(kStdWarn,'(A)') '-------------------------------------------------'
              WRITE(kStdWarn,'(A)') 'Index  GasID  UnitsCode  GasUnit'
              WRITE(kStdWarn,'(A)') '-------------------------------------------------'
              DO iG = 1,iNumGases
                WRITE(kStdWarn,'(3(I5),A,A)') iG,iaG(iG),iaGasUnits(iG),'      ',caaUnit(iaGasUnits(iG))
              END DO
              WRITE(kStdWarn,'(A)') '-------------------------------------------------'
              
!! this works fine if levels entered so thatpres(1) > pres(2) > pres(3) ie from ground up
!! may be buggy for reverse (from TOA down)
              WRITE(kStdWarn,*) 'Input levels profile : ',iNumLevs,' levels for ',iNumGases,' gases'
              WRITE(kStdWarn,*) 'PSurf = ',rPSurf,' mb;  TSurf = ',rTSurf,' K '
              iReadInOK = 0
              iIgnore = -1
              DO iL = 1,iNumLevs
                READ (iIOUN2,*) rP,rT,(raX(iJ),iJ=1,iNumGases)
                IF ((rP <= rPRefMax) .AND. (iIgnore < 0)) THEN
                  IF ((rP < rPRefMin) .AND. (rP < rPmin)) iIgnore = +1
                  iReadInOK = iReadInOK + 1
                  raP_Nm2(iReadInOK) = rP * 100.0  !! change from mb to N/m2
                  raP(iReadInOK) = rP
                  raT(iReadInOK) = rT
                  IF (rPmax <= raP(iReadInOK)) rPmax = raP(iReadInOK)
                  IF (rPmin > raP(iReadInOK)) rPmin = raP(iReadInOK)
                  DO iJ = 1,iNumGases
                    raaG_VMR(iReadInOK,iJ) = raX(iJ)
                  END DO
                  WRITE(kStdWarn,5040) iReadInOK,raP(iReadInOK),raT(iReadInOK),raaG_VMR(iReadInOK,1)
                ELSE IF (rP <= rPRefMin) THEN
!! things are OK, we are ignoring pressures for levels above 0.005 mb
                  WRITE(kStdWarn,5041) iL,rP,rT
                ELSE IF ((rP > rPRefMin) .AND. (iIgnore > 0)) THEN
!! things are BAD, this pressure is for levels below 0.005 mb
!! but we already found a pressre for above 0.005 mb, so pinput is UNSORTED??
                  WRITE(kStdWarn,5042) iL,rP,rT
                  CALL DoStop
                END IF
              END DO
              
              13   CONTINUE
              CLOSE(iIOUN2)
              kProfileUnitOpen = -11
              
              IF (iNumLevs /= iReadInOK) THEN
                WRITE(kStdWarn,*) 'Input textfile had ',iNumLevs,' levels of which ',iReadInOK,' lay between ',rPRefMax,rPRefMin,' mb'
              ELSE
                WRITE(kStdWarn,*) 'All ',iNumLevs,' input levels in textfile lay between ',rPRefMax,rPRefMin,' mb'
              END IF
              iNumLevs = iReadInOK
              
              5030 FORMAT(A80)
              5040 FORMAT('iReadInOK, rP (mb), rT (K), raG1 = ',I3,' ',F14.6,' ',F8.3,' ',ES12.5)
              5041 FORMAT('Ignore (too low P) rP (mb), rT (K), raG1 = ',I3,' ',ES12.5,' ',F8.3)
              5042 FORMAT('Bad (press below 0.005 mb (UNSORTED PinLEVS??)) rP (mb), rT (K), raG1 = ',I3,' ',ES12.5,' ',F8.3)
              
              WRITE(kStdWarn,*)'   KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
              WRITE(kStdWarn,*)'   kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
              WRITE(kStdWarn,*)' input file : spres (mb)/Height(km)     = ',rPSurf,rHSurf
              
              rPSurf = rPSurf * 100.0
              rHSurf = rHSurf * 1000.0
              
! make sure pressures are decreasing with index ie layers going higher and higher
              IF (raP(1) < raP(2)) THEN
!!! need to swap!!!
                iMid = ifloor(iNumLevs/2.0)
                DO iL = 1,iMid
                  iJ = iNumLevs-iL+1
                  rX = raP(iJ)
                  raP(iJ) = raP(iL)
                  raP(iL) = rX
                  
                  rX = raP_Nm2(iJ)
                  raP_Nm2(iJ) = raP_Nm2(iL)
                  raP_Nm2(iL) = rX
                  
                  rX = raT(iJ)
                  raT(iJ) = raT(iL)
                  raT(iL) = rX
                  
                  DO iG = 1,iNumGases
                    raX(iG) = raaG_VMR(iJ,iG)
                    raaG_VMR(iJ,iG) = raaG_VMR(iL,iG)
                    raaG_VMR(iL,iG) = raX(iG)
                  END DO
                END DO
              END IF
              
              WRITE(kStdWarn,*) 'input file : Highest altitude (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
              WRITE(kStdWarn,*) ' '
              
              iYes = +1    !!! reset topmost pressure to rPRefMin = 0.005
              iYes = -1    !!! leave topmost pressure level alone, do not reset to rPRefMin = 0.005
!!! actually we augment the refdatabase to 0.2750000 N/m2 = 0.00275
              rPRefMinAugmented = 0.00275
              IF ((rPmin < rPRefMin) .AND. (iYes > 0)) THEN
!!! interpolate gas concentrations last point to rPRefMin
                WRITE(kStdWarn,*) 'highest press point is ',rPmin,' mb higher than ref database point ',rPRefMin,'mb fixing!'
                WRITE(kStdWarn,*) ' '
                DO iG = 1,iNumGases
                  DO iL = 1,iNumLevs
                    raX(iL) = raaG_VMR(iL,iG)
                  END DO
                  CALL r_sort_loglinear1(raP,raX,iNumLevs,rPRefMin,rX,1)
                  raaG_VMR(iNumLevs,iG) = rX
                END DO
                
!!! interpolate T(z) last point to rPRefMin
                CALL r_sort_loglinear1(raP,raT,iNumLevs,rPRefMin,rX,1)
                raT(iNumLevs) = rX
                
!!! set last point to 0.005 mb
                rPmin = rPRefMin
                raP(iNumLevs) = rPRefMin
                raP_Nm2(iNumLevs) = rPRefMin * 100.0
                
              ELSE IF (rPmin < rPRefMinAugmented) THEN
!!! force interpolate gas concentrations last point to rPRefMinAugmented
                WRITE(kStdWarn,*) 'highest press point is ',rPmin,' mb higher than augmented ref dbase pt ',rPRefMinAugmented,'mb fixing!'
                WRITE(kStdWarn,*) ' '
                DO iG = 1,iNumGases
                  DO iL = 1,iNumLevs
                    raX(iL) = raaG_VMR(iL,iG)
                  END DO
                  CALL r_sort_loglinear1(raP,raX,iNumLevs,rPRefMinAugmented,rX,1)
                  raaG_VMR(iNumLevs,iG) = rX
                END DO
                
!!! interpolate T(z) last point to rPRefMin
                CALL r_sort_loglinear1(raP,raT,iNumLevs,rPRefMinAugmented,rX,1)
                raT(iNumLevs) = rX
                
!!! set last point to 0.00275 mb
                rPmin = rPRefMinAugmented
                raP(iNumLevs) = rPRefMinAugmented
                raP_Nm2(iNumLevs) = rPRefMin * 100.0
              END IF
              
              
! now change all units to PPMV (units 10) and thn to VMR (units 12)
              DO iG = 1,iNumGases
                CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP_Nm2,raT,raaG_VMR,-1)
                DO iL = 1,iNumLevs
                  raaG_VMR(iL,iG) =  raaG_VMR(iL,iG) / 1.0E6
                  iaGasUnits(iG) = 12   !! VMR
                END DO
              END DO
              
              RETURN
            END SUBROUTINE  ReadInput_LVL_Profile
            
!************************************************************************
            
            SUBROUTINE ReadRefProf_Levels2(iGasID,raP,iNumLevsIN,raX)
            
            
            INTEGER, INTENT(IN OUT)                  :: iGasID
            REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
            INTEGER, INTENT(IN)                      :: iNumLevsIN
            REAL, INTENT(OUT)                        :: raX(2*kProfLayer)
            IMPLICIT NONE
            INCLUDE '../INCLUDE/kcartaparam.f90'
            
! input
            
            
! output
            
            
! local
            REAL :: raRx110Temp(2*kProfLayer),raRx110MR(2*kProfLayer),raRx110Press(2*kProfLayer)
            REAL :: rLat,raJunk(kProfLayer*2)
            CHARACTER (LEN=80) :: caPFname,caStr,caComment
            INTEGER :: iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec,iNumLevsx
            
            caPFname = kcaLevsRefProf
            
            IF ((kAFGLProf < 1) .OR. (kAFGLProf > 6)) THEN
              WRITE(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
              CALL DoStop
            END IF
            
            iErr = 0
            iIOUN2 = kProfileUnit
            OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',  &
                IOSTAT=iErrIO)
            IF (iErrIO /= 0) THEN
              iErr=1
              WRITE(kStdErr,1070) iErrIO, caPfname
              CALL DoSTOP
            END IF
            1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
            kProfileUnitOpen=1
            
            iAFGL = 0
            iFoundGas = -1
            
            10  CONTINUE
            READ(iIOUN2,1080) caStr
            IF (caStr(1:1) == '!') THEN
              GO TO 10   !! keep reading comments
            ELSE
              READ (caStr,*) iNumLevsx
            END IF
            IF (iNumLevsx > 2*kProfLayer) THEN
              WRITE(kStdErr,*) 'oops iNumLevsx .GT. 2*kProfLayer ',iNumLevsx,kProfLayer*2
              CALL DoStop
            END IF
            
            20  CONTINUE
            iAFGL = iAFGL + 1
            
            30  CONTINUE
            READ(iIOUN2,1080) caStr
            IF (caStr(1:1) == '!') THEN
              GO TO 30   !! keep reading comments
            ELSE
              caComment = caStr
            END IF
            
            40  CONTINUE
            READ(iIOUN2,1080) caStr
            IF (caStr(1:1) == '!') THEN
              GO TO 40   !! keep reading comments
            ELSE
              READ (caStr,*) rLat
            END IF
!      write(kStdWarn,*) 'iAFGL Profile ',iAFGL,' rLat ',rLat, ' : ',caComment
            
            READ(iIOUN2,1080) caStr   !! comment
            READ(iIOUN2,1080) caStr   !! Altitude(km)
            READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)
            
            IF (iAFGL == kAFGLProf) THEN
!        write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! Press(mb)
              READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)
              
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! Temp(K)
              READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)
              
!        write(kStdWarn,*) '   need to find profile for gasID ',iGasID
            ELSE
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! Press(mb)
              READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
              
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! Temp(K)
              READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
            END IF
            
            READ(iIOUN2,1080) caStr   !! comment
            READ(iIOUN2,1080) caStr   !! density (cm-3)
            READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)
            
            DO iG = 1,7
              IF (iG == iGasID) THEN
!          write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
                iFoundGas = +1
                READ(iIOUN2,1080) caStr   !! comment
                READ(iIOUN2,1080) caStr   !! gas name
                READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
              ELSE
                READ(iIOUN2,1080) caStr   !! comment
                READ(iIOUN2,1080) caStr   !! gas name
                READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
              END IF
            END DO
            
            IF ((iAFGL == kAFGLProf) .AND. (iFoundGas > 0)) THEN
              GO TO 60   !! found the AFGL prof and found the gas; done!
            ELSE IF ((iAFGL == 6) .AND. (iFoundGas < 0)) THEN
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! modend
              GO TO 50   !! found the AFGL prof but not found the gas
            ELSE IF ((iAFGL < 6) .OR. (iFoundGas < 0)) THEN
!! either did not find the gas or the AFGL prof
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! modend
              GO TO 20
            END IF
            
            50   CONTINUE
            READ(iIOUN2,1080) caStr   !! comment
            READ(iIOUN2,1080) caStr   !! constituent profs
            READ(iIOUN2,1080) caStr   !! comment
            READ(iIOUN2,1080) caStr   !! mingas
            
            iG = 7
            55   CONTINUE
            iG = iG + 1
            IF (iG == iGasID) THEN
!        write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
              iFoundGas = +1
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! gas name
              READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
            ELSE
              READ(iIOUN2,1080) caStr   !! comment
              READ(iIOUN2,1080) caStr   !! gas name
              READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
            END IF
            IF ((iG < 28) .AND. (iFoundGas < 0)) THEN
              GO TO 55
            END IF
            
            IF (iFoundGas < 0) THEN
!! bweh try xsec gases
              123    CONTINUE
              READ(iIOUN2,1080) caStr   !! comment
              IF (caStr(1:1) == '!') THEN
                GO TO 123
              ELSE IF (caStr(1:6) == 'DATEND') THEN
                GO TO 60       !!! oh oh end of file, give up!!!!!
              ELSE
                READ (caStr,*) iXsec
                IF (iXsec == iGasID) THEN
                  iFoundGas = 1
                  READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
                ELSE
                  READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
                  GO TO 123
                END IF
              END IF
            END IF
            
            60   CONTINUE
            CLOSE(iIOUN2)
            kProfileUnitOpen=-1
            1080 FORMAT(A80)
            
            IF (iFoundGas < 0) THEN
!! finally give up
              WRITE(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
              CALL DOStop
            END IF
            
!*************************
! finally interp these onto the raP pressure levels
            DO iI = 1,iNumLevsx
              raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
!        print *,iI,raRx110Press(iI),raRx110Temp(iI),raRx110MR(iI)
            END DO
            
!      CALL r_sort_loglinear(raRx110Press,raRx110Temp,iNumLevsx,PLEVx110,raJunk,2*kProfLayer)
!      DO iI = 1,2*kProfLayer
!        raRx110Temp(iI) = raJunk(iI)
!      END DO
            CALL r_sort_loglinear(raRx110Press,raRx110MR,iNumLevsx,raP,raJunk,iNumLevsIN)
            DO iI = 1,iNumLevsIN
              raX(iI) = raJunk(iI)
            END DO
!      iNumLevsx = 2*kProfLayer
            
            RETURN
          END SUBROUTINE ReadRefProf_Levels2
!************************************************************************