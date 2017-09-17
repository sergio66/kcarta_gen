! this is a much more detailed LBLRTM TAPE5 reader
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:46
 
! and hmm only seems to work for case mentioned below??? and is not called anywhere else
! see https://svn.ssec.wisc.edu/repos/uwphysret/trunk/mfiles/iasi_radpress_write5.m

SUBROUTINE  ReadInput_DetailedLBLRTM_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
    iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,raP,raT,  &
    iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)


NO TYPE, INTENT(IN OUT)                  :: caPFname
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
INTEGER, INTENT(OUT)                     :: iNumLevs
REAL, INTENT(OUT)                        :: rPSurf
REAL, INTENT(OUT)                        :: rTSurf
REAL, INTENT(OUT)                        :: rHSurf
INTEGER, INTENT(OUT)                     :: iNumGases
REAL, INTENT(OUT)                        :: raP(2*kProfLayer)
REAL, INTENT(OUT)                        :: raT(2*kProfLayer)
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
REAL, INTENT(OUT)                        :: raaG_MR(2*kProfLayer,kMaxGas)
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

REAL :: rPmin,rPmax


! local var
INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
REAL :: raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH
CHARACTER (LEN=80) :: caStr,caStrX,caStrY
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=1) :: c1

rPmin = +1.0E6
rPmax = -1.0E+6
iNXsec = -1

111  FORMAT(A80)

iIOUN2 = kProfileUnit
OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED', IOSTAT=iErrIO)
IF (iErrIO /= 0) THEN
  iErr=1
  WRITE(kStdErr,1070) iErrIO, caPfname
  1070     FORMAT('ERROR! number ',I5,' opening PRFILE path LBLRTM profile file'  &
      ,/,A80)
  CALL DoSTOP
END IF
kProfileUnitOpen=1

! this should read /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
READ (iIOUN2,111,ERR=13,END=13) caStr
READ (iIOUN2,111,ERR=13,END=13) caStr

READ (iIOUN2,*) rF1,rF2
WRITE(kStdWarn,*) 'LBLRTM input indicates start/stop wavenumbers are ',rF1,rF2

READ (iIOUN2,*) rTSurf
WRITE(kStdWarn,*) 'LBLRTM input indicates STEMP is ',rTSurf

READ (iIOUN2,*) (iaJunk(iL),iL = 1,9)
iNumLevs = iaJunk(3)
iNumGases = iaJunk(6)
IF (iNumLevs < 0) THEN
  WRITE(kStdWarn,*) 'looks like boundaries for calculations specified in mb'
  iNumLevs = ABS(iNumLevs)
  iLBROutBdryHorP = -1
ELSE
  iLBROutBdryHorP = +1
  WRITE(kStdWarn,*) 'looks like boundaries for calculations specified in km'
END IF

WRITE(kStdWarn,*) 'TAPE5 implies kCARTA compute ODs at ',iNumLevs,' press/heights, for iNumGases = ',iNumGases
IF (iNumLevs > kProfLayer) THEN
  WRITE(kStdErr,*) 'though irrelevant for kCARTA calcs, it will not be able to read in the pressures/heights'
  WRITE(kStdErr,*) iNumLevs,kProfLayer
  CALL DoStop
END IF
IF (iNumGases > kMaxGas) THEN
  WRITE(kStdErr,*) 'iNumGases .GT. kMaxGas',iNumGases,kMaxGas
  CALL DoStop
END IF

READ (iIOUN2,*) rTophgt,rHSurf,rViewAngle
WRITE(kStdWarn,*) 'LBLRTM input indicates start/stop hgts are ',rHSUrf,rTopHgt,' with view angle ',rViewAngle

READ (iIOUN2,*) (raLBL_Hgts(iJ),iJ=1,iNumLevs)

raRTP_TxtInput(6) = +raLBL_Hgts(iNumLevs)    !! default, assume hgt in km
IF  (iLBROutBdryHorP < 0) raRTP_TxtInput(6) = -raLBL_Hgts(iNumLevs)

READ (iIOUN2,*) iNumLevs
IF (iNumLevs > 2*kProfLayer) THEN
  WRITE(kStdErr,*) 'iNumLevs .GT. 2*kProfLayer',iNumLevs,kProfLayer
  CALL DoStop
END IF

! see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
!       JCHAR = 1-6           - default to value for specified model atmosphere
!              = " ",A         - volume mixing ratio (ppmv)
!              = B             - number density (cm-3)
!              = C             - mass mixing ratio (gm/kg)
!              = D             - mass density (gm m-3)
!              = E             - partial pressure (mb)
!              = F             - dew point temp (K) *H2O only*
!              = G             - dew point temp (C) *H2O only*
!              = H             - relative humidity (percent) *H2O only*
!              = I             - available for user definition
DO iG = 1,iNumGases
  iaG(iG) = iG
  iaGasUnits(iG) = 10   !!! assume hardcoded ppmv
END DO

DO iL = 1,iNumLevs
! READ (iIOUN2,*) rH,rP,rT,caStrX
  READ (iIOUN2,111) caStrY
  READ(caStrY,*) rH,rP,rT
  caStrX = caStrY(31:80)
  caStr30 = caStrX(11:40)
  
  IF (iL == 1) THEN
    DO iJ=1,iNumGases
      c1 = caStr30(iJ:iJ)
      IF (c1 == 'A') iaGasUnits(iJ) = 10
      IF (c1 == ' ') iaGasUnits(iJ) = 10
      IF (c1 == 'B') iaGasUnits(iJ) = -1
      IF (c1 == 'C') iaGasUnits(iJ) = 20
      IF (c1 == 'D') iaGasUnits(iJ) = -1
      IF (c1 == 'E') iaGasUnits(iJ) = -1
      IF (c1 == 'F') iaGasUnits(iJ) = 42
      IF (c1 == 'G') iaGasUnits(iJ) = 43
      IF (c1 == 'H') iaGasUnits(iJ) = 40
      IF (c1 == 'I') iaGasUnits(iJ) = -1
      
      IF (iaGasUnits(iJ) < 0) THEN
        WRITE(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
        CALL DOStop
      ELSE
        WRITE(kStdWarn,*) 'LBLRTM gas number ',iJ,' "units" ',c1,' = kCARTA levels units of ',iaGasUnits(iJ)
      END IF
    END DO
  END IF
  
  READ (iIOUN2,*) (raX(iG),iG=1,iNumGases)
!        raP(iL) = rP * 100.0  !! change from mb to N/m2
  raP(iL) = rP          !! keep in mb
  raT(iL) = rT
  IF (rPmax <= raP(iL)) rPmax = raP(iL)
  IF (rPmin > raP(iL)) rPmin = raP(iL)
  DO iJ = 1,iNumGases
    raaG_MR(iL,iJ) = raX(iJ)
  END DO
  
  IF (iL == 1) THEN
!          rPSurf = raP(1)/100 !! because we redo this below
    rPSurf = raP(1)     !! keep in mb
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
  END IF
END DO

!! now see if there are xsec gases
READ (iIOUN2,111,ERR=13,END=13) caStr
READ(caStr,*) iNXsec
IF (iNXsec > 0) THEN
  READ (iIOUN2,111,ERR=13,END=13) caStr     !!!! xsec names
  CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
  READ (iIOUN2,*) iNumLevsXsec
  IF (iNumLevsXsec > 2*kProfLayer) THEN
    WRITE(kStdErr,*) 'iNumLevsXsec .GT. 2*kProfLayer',iNumLevsXsec,kProfLayer
    CALL DoStop
  END IF
  IF (iNumLevsXsec /= iNumLevs) THEN
    WRITE(kStdErr,*) 'iNumLevsXsec .NE. iNumLevs',iNumLevsXsec,iNumLevs
    CALL DoStop
  END IF
  DO iL = 1,iNumLevsXsec
    READ (iIOUN2,111,ERR=13,END=13) caStrY
    caStr30 = caStrY(16:45)
    IF (iL == 1) THEN
      DO iJ=1,iNXsec
        c1 = caStr30(iJ:iJ)
        IF (c1 == 'A') iaGasUnits(iJ+iNumGases) = 10
        IF (c1 == ' ') iaGasUnits(iJ+iNumGases) = 10
        IF (c1 == 'B') iaGasUnits(iJ+iNumGases) = -1
        IF (c1 == 'C') iaGasUnits(iJ+iNumGases) = 20
        IF (c1 == 'D') iaGasUnits(iJ+iNumGases) = -1
        IF (c1 == 'E') iaGasUnits(iJ+iNumGases) = -1
        IF (c1 == 'F') iaGasUnits(iJ+iNumGases) = 42
        IF (c1 == 'G') iaGasUnits(iJ+iNumGases) = 43
        IF (c1 == 'H') iaGasUnits(iJ+iNumGases) = 40
        IF (c1 == 'I') iaGasUnits(iJ+iNumGases) = -1
        
        IF (iaGasUnits(iJ+iNumGases) < 0) THEN
          WRITE(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
          CALL DOStop
        ELSE
          WRITE(kStdWarn,*) 'LBLRTM xsec number ',iJ,' is of "units" ',c1,' = kCARTA input levels units of ',  &
              iaGasUnits(iJ+iNumGases)
        END IF
      END DO
    END IF
    
    READ (iIOUN2,*) (raX(iG),iG=1,iNXsec)
    DO iJ = 1,iNXsec
      raaG_MR(iL,iJ+iNumGases) = raX(iJ)
    END DO
  END DO
  iNumGases = iNumGases + iNXsec
END IF

13   CONTINUE
CLOSE(iIOUN2)
kProfileUnitOpen = -11

raRTP_TxtInput(1) = rPSurf
raRTP_TxtInput(2) = rTSurf
raRTP_TxtInput(3) = rHSurf    !! km
raRTP_TxtInput(4) = rTophgt   !! km
raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, downwell rad, else upwell rad

rPSurf = rPSurf * 100.0    !! change from mb to N/m2
rHSurf = rHSurf * 1000.0   !! change from km to m

WRITE(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
WRITE(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
WRITE(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf

! make sure pressures are decreasing with index ie layers going higher and higher
IF (raP(1) < raP(2)) THEN
!!! need to swap!!!
  iMid = ifloor(iNumLevs/2.0)
  DO iL = 1,iMid
    iJ = iNumLevs-iL+1
    rX = raP(iJ)
    raP(iJ) = raP(iL)
    raP(iL) = rX
    
    rX = raT(iJ)
    raT(iJ) = raT(iL)
    raT(iL) = rX
    
    DO iG = 1,iNumGases
      raX(iG) = raaG_MR(iJ,iG)
      raaG_MR(iJ,iG) = raaG_MR(iL,iG)
      raaG_MR(iL,iG) = raX(iG)
    END DO
  END DO
END IF
WRITE(kStdWarn,*) 'input file : Highest altitude (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
WRITE(kStdWarn,*) ' '

! now change all units to MR
DO iG = 1,iNumGases
  CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR)
  DO iL = 1,iNumLevs
    raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0E6
    iaGasUnits(iG) = 12
  END DO
END DO

RETURN
END SUBROUTINE  ReadInput_DetailedLBLRTM_Profile

!************************************************************************
! this subroutine does the klayers integrations, integrals over P
! so this is more similar to kLAYERs, CRTM
! this is orig code; it puts surface as "lowest level and adds the user level info above that

SUBROUTINE DoIntegrateLevels2Layers_wrtP_pSurf(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,  &
    PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,  &
    raPX,raTX,raaG_MRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)


REAL, INTENT(IN OUT)                     :: rHSurf
REAL, INTENT(IN)                         :: rPSurf
REAL, INTENT(IN)                         :: rTSurf
INTEGER, INTENT(IN)                      :: iLowestLev
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
REAL, INTENT(IN OUT)                     :: raaG_MRX(kProfLayer+1,kMaxGas)
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
REAL :: raaXYZ_MR(kMaxlayer+1+1,kMaxGas)
REAL :: raXYZ_MRwater(kMaxlayer+1+1)
REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1
REAL :: zWoo,rConvertQ,grav_earth,wexsvp,rRH
INTEGER :: iWoo

rConvertQ = 6.022141E23/1E4     ! multiply by this to change moles/m2 to molecules/cm2
rConvertQ = kAvog/1000/1.0E4    ! multiply by this to change moles/m2 to molecules/cm2

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

! >>>>>>>>>>
!! set up temporary arrays/matrices for interping onto sublevels
!! set first level at surface press, surface temp
raXYZPress(1) = rPSurf
raXYZTemp(1)  = rTSurf
iUpperLev = kMaxLayer-iLowestLev+1
IF ((iUpperLev > kMaxLayer) .OR. (iUpperLev < 1)) THEN
  WRITE (kStdErr,*) '(iUpperLev .GT. kMaxLayer) .OR. (iUpperLev .LT. 1)',iUpperLev,kMaxLayer
  CALL Dostop
END IF
!! do surface by interpolating MR from user supplied/kcarta_database levels, to surface pressure
DO iG = 1,iNumGases
  DO iL = 1,iUpperLev
    raJunk(iL) = raaG_MRX(iL,iG)
  END DO
  CALL r_sort_loglinear(raPX,raJunk,iUpperLev,rPSurf,junk,1)
  raaXYZ_MR(1,iG) = junk
  IF (iG ==  1) raXYZ_MRwater(1) = junk
END DO
! >>>>>>>>>>
!! set levels 2 .. iNumLevs+1
DO iL = 1,iUpperLev
  raXYZPress(iL+1) = raPX(iL)
  raXYZTemp(iL+1) = raTX(iL)
  DO iG = 1,iNumGases
    raaXYZ_MR(iL+1,iG) = raaG_MRX(iL,iG)
    IF (iG ==  1) raXYZ_MRwater(iL+1) = raaG_MRX(iL,iG)
  END DO
END DO
! >>>>>>>>>>

! now start integrating
! lowest layer is special case

! http://cimss.ssec.wisc.edu/~paulv/
! http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
!     /home/sergio/KCARTA/DOC/sci_klayers.txt

z = rHSurf !! need to kludge this
zWoo = rHSurf
iWoo = -1

DO iL = iLowestLev,iLowestLev
  raPout(iL) = LOG(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
  raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
  
  z0 = z
  iJ = iL - iLowestLev + 1
  
  amount = 0.0
  rTsum  = 0.0
  rWgt = 0.0
  
!!! >>>>>>  this is starting out at         rP = PLEV_KCARTADATABASE_AIRS(iL)        rT = rTSurfx  >>>>>>>>
  rP = PLEV_KCARTADATABASE_AIRS(iL)
  CALL r_sort_loglinear(raXYZPress,raXYZTemp,iUpperLev+1,rP,rT,1)
  
  dlnp = LOG(PLEV_KCARTADATABASE_AIRS(iL)) - LOG(PLEV_KCARTADATABASE_AIRS(iL+1))
  dlnp = dlnp / (iNFine)
  
!! information for current (sub)level
  rP_n = rP                             !! current sublev press
  rT_n = rT                             !! current sublev temp
  rMR_water_n = raXYZ_MRwater(1)        !! current water MR
  CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP,rMR_water_n,1) !! current water MR
  
!! information for next (sub)level
  rP_np1 = LOG(rP_n) - dlnp
  rP_np1 = EXP(rP_np1)                                                              !! next sublev pressure
  CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
  CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw
  
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
      DO iCnt = 1,iUpperLev+1
        raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
      END DO
      CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_n,  rMR_n,  1)
      CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_np1,rMR_np1,1)
      raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
    END DO
    
!!! now update for next iteration
    z = z + dz
    rP_n = rP_np1
    rT_n = rT_np1
    rMR_water_n = rMR_water_np1
    rP = rP_n
    
    rP_np1 = LOG(rP_n) - dlnp
    rP_np1 = EXP(rP_np1)                                                              !! next sublev pressure
    CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
    CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw
    
    IF ((rP_n <= rPSurf) .AND. (iWoo < 0)) THEN
      iWoo = +1
      zWoo   = z
    END IF
    
  END DO
  
  raTout(iL)      = rTsum/rWgt                 ! <<< need TAV
  raAmountOut(iL) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
  raZout(iL+1)    = z
  
!        print *,iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
!     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
END DO

! displace everything  z
iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
raZout(iL)   = -(zWoo - z0)
raZout(iL+1) = +(z-zWoo)
WRITE(kStdWarn,*) '  need to displace lowest layer heights'
WRITE(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
WRITE(kStdWarn,*) '  Plowest,pSurf,pHighest (mb) ',PLEV_KCARTADATABASE_AIRS(iL),rPSurf,  &
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
  
  iJ = iL - iLowestLev + 1
  rP = PLEV_KCARTADATABASE_AIRS(iL)
  rT = raXYZTemp(iJ)
  
  amount = 0.0
  rTsum  = 0.0
  rWgt = 0.0
  
  dlnp = LOG(PLEV_KCARTADATABASE_AIRS(iL)) - LOG(PLEV_KCARTADATABASE_AIRS(iL+1))
  dlnp = dlnp / (iNFine)
  
!! information for current (sub)level
  rP_n = rP                              !! current sublev press
  rT_n = rT                              !! current sublev temp
  rMR_water_n = raXYZ_MRwater(iJ)        !! current water MR
  
!! information for next (sub)level
  rP_np1 = LOG(rP_n) - dlnp
  rP_np1 = EXP(rP_np1)                                                              !! next sublev pressure
  CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
  CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw
  
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
      DO iCnt = 1,iUpperLev+1
        raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
      END DO
      CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_n,  rMR_n,  1)
      CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_np1,rMR_np1,1)
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
    rP_np1 = EXP(rP_np1)                                                              !! next sublev pressure
    CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
    CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw
    
  END DO
  
  raTout(iL)      = rTsum/rWgt             ! <<< need TAV
  raAmountOut(iL) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
  raZout(iL+1)    = z
  
!        print *,iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
!     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
END DO

WRITE(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
!! raZout is in meters
!! raaQout is in moles/m2
DO iL = iLowestLev,kProfLayer
  rThick = ABS(raZout(iL) - raZout(iL+1)) * 100.0
  rT = raTout(iL)
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
!     $        raaG_MRX(iL-iLowestLev+1,1),raaG_MRX(iL-iLowestLev+1,2),raaG_MRX(iL-iLowestLev+1,3),
!     $        raaG_MRX(iL-iLowestLev+1,4),raaG_MRX(iL-iLowestLev+1,5)
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
      raaG_MRX(iL,iG) = 0.0
    END DO
  END DO
END IF

IF (kPlanet == 3) THEN
!! show the RH for gas 1
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*) 'Checking Relative Humidity'
  WRITE(kStdWarn,*) ' Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      Q         Q1         Q2        Q3'
  WRITE(kStdWarn,*) '-------------------------------------------------------------------------------------------------------'
  DO iL = iLowestLev,kProfLayer
    z    = (raZout(iL) + raZout(iL+1))/2/1000
    zWoo = (raZout(iL+1) - raZout(iL))/1000
    rP   = raPout(iL)            !! N/m2
    rT   = raTout(iL)
    rPP  = raaPartPressOut(iL,1) !! N/m2
    rSVP = wexsvp(rT) * 100      !! mb --> N/m2
    rRH  = rPP/rSVP*100.0
!          print *,iL,rP/100,rT,rPP/100,rRH
!          IF (rRH .LE. 100.0) THEN
!            write(kStdWarn,111) iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),(raaQout(iL,iG),iG=1,3)
!          ELSE
!            write(kStdWarn,112) iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),(raaQout(iL,iG),iG=1,3),' ***** '
!          END IF
  END DO
END IF
111  FORMAT(I3,' ',1(F9.5,' '),6(F8.4,' '),4(ES9.4,' '))
112  FORMAT(I3,' ',1(F9.5,' '),6(F8.4,' '),4(ES9.4,' '),A7)

RETURN
END SUBROUTINE DoIntegrateLevels2Layers_wrtP_pSurf

!************************************************************************
! this subroutine will take in 100 AIRS layering stuff and interpolate to
! the new arbitrary layering
! WARNING : this assumes that the user has not mucked up KLAYERS layering
!           such that highest Z pressure (lowest pressure) is NOT TOA
!           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
! do the lower atm (usual -1) or upper atm (NLTE +1)

! kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
! so at arbitrary pressure layering, it deals with abs coeffs
! so we do not need raRamt
! but we do need the interpolated temp and partial pressures

SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,  &
    raR100Amt,raR100Temp,raR100Press,raR100PartPress,  &
    raaPress,iGas,iGasID,iNumLayers,  &
    raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)


REAL, INTENT(OUT)                        :: raRAmt(kProfLayer)
REAL, INTENT(OUT)                        :: raRTemp(kProfLayer)
REAL, INTENT(OUT)                        :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(IN)                         :: raR100Amt(kMaxLayer)
REAL, INTENT(IN)                         :: raR100Temp(kMaxLayer)
NO TYPE, INTENT(IN OUT)                  :: raR100Pres
NO TYPE, INTENT(IN OUT)                  :: raR100Part
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN)                      :: iGasID
INTEGER, INTENT(IN)                      :: iNumLayers
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
INTEGER, INTENT(IN OUT)                  :: iError
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

!  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
!  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
! The above equation and 3 data points define the 101 AIRS levels, which
! are in airslevelsparam.f90

! input
! do the lower atm (usual -1) or upper atm (NLTE +1)
INTEGER :: iLowerOrUpper
! these are the individual reference profiles, at kMaxLayer layers

REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the arbitrary profiles stored in matrices

INTEGER :: iSplineType
! these are the kLAYERS pressure levels, layer thick for the current profile
REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
!  output
! these are the individual reference profiles, at kProfLayer layers

REAL :: raRPartPress(kProfLayer)
REAL :: pMax100,pMin100

! local variables
INTEGER :: iI,iNot,kMaxLayer_over_2,iJ
REAL :: raWorkP(kMaxLayer),raXgivenP(kMaxLayer),  &
    raYgivenP(kMaxLayer),raY2P(kMaxLayer)
REAL :: raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2
REAL :: raSortPressLevels(kMaxLayer+1)
REAL :: raSortPressHeights(kMaxLayer+1)

REAL :: raDataBaseThickness(kMaxLayer)
REAL :: DatabaseHeight(kMaxLayer)
REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)
REAL :: DATABASELEV(kMaxLayer+1)

! pressure variables!!!!! ----------------->
! raaPress in atm

CALL databasestuff(iLowerOrUpper,  &
    DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

!     Assign values for interpolation
!     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
rYP1 = 1.0E+16
rYPN = 1.0E+16

!!! find the thickness of the pressure layers in centimeters
!!! hence we have x 1000 (km -> m) x 100 (m -> cm)
DO iI = 1,kMaxLayer
  raDataBaseThickness(iI) = (DATABASELEVHEIGHTS(iI+1) -  &
      DATABASELEVHEIGHTS(iI))*1000.0 * 100.0
END DO

! now store stuff sorted in terms of increasing pressure
kMaxLayer_over_2 = 50
DO iI = 1,kMaxLayer+1
  raSortPressLevels(iI)  = DataBaseLev(101-iI+1)                  !!! mb
  raSortPressHeights(iI) = DataBaseLevHeights(101-iI+1) * 1000.0  !!! m
END DO

!!!this tells how many layers are NOT dumped out by kLAYERS
iNot = kProfLayer-iNumLayers
! now just happily spline everything on!!!!!! for the amts
! recall you need raXgivenP to be increasing
pMax100 = -1.0E10
pMin100 = +1.0E10
DO iI = 1,kMaxLayer
  raXgivenP(iI) = LOG(raR100Press(kMaxLayer-iI+1))
  raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
  raYgivenP(iI) = raR100Amt(kMaxLayer-iI+1)
  IF (raXgivenP(iI) < pMin100) THEN
    pMin100 = raXgivenP(iI)
  END IF
  IF (raXgivenP(iI) > pMax100) THEN
    pMax100 = raXgivenP(iI)
  END IF
END DO
CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
DO iI = 1,iNot
  rxpt = LOG(raaPress(iI,iGas))
  rxpt = raaPress(iI,iGas)
  r = 0.0
  raRAmt(iI) = r
END DO

! these are the layers with info
DO iI = iNot+1,kProfLayer
  rxpt = LOG(raaPress(iI,iGas))
  rxpt = raaPress(iI,iGas)
  IF (iSplineType == +1) THEN
    CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
  ELSE
    CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rXPT,r)
  END IF
  raRAmt(iI) = r
!^^^^^^^^^^^^^ check to see if linear interp will help^^^^^^^^^^^^^^^^^^^
  IF ((r < 0.0)) THEN
    WRITE (kStdWarn,*) 'Making reference profile : negative amt found : '
    WRITE (kStdWarn,*) 'GasID,layer,ref amt,iLorU (-1/+1) : ',iGasID,iI,r,iLowerOrUpper
    CALL rlinear(raXgivenP,raYgivenP,kMaxLayer,rxpt,r,1)
    WRITE (kStdWarn,*) '   trying linear interp : ',iI,r
    raRAmt(iI) = r
    IF ( (r < 0.0). AND. ((rxpt > pMin100) .AND. (rxpt < pMax100))) THEN
!!!things barfed even though pMin100 < rXpt < pMax100 ... bad
      WRITE (kStdErr,*) 'gasID,layer,ref amount(linear interp)= ', iGasID,iI,r
      WRITE (kStdErr,*) 'min(pD),rXpt,max(pD) = ',pMin100,rXpt,pMax100
      CALL DoStop
    ELSE IF ( (r < 0.0). AND.  &
          ((rxpt < pMin100) .OR. (rxpt > pMax100))) THEN
!!!things barfed, pMin100 > rXpt or pMax100 < rXpt ... hmmm
      WRITE (kStdWarn,*) 'gasID,layer,ref amount(linear interp)= ',  &
          iGasID,iI,r
!!!recall raYgiven(iI) is swtiched GND = kMaxLayer,100=TOA
      IF (rxpt < pMin100) r = raYgivenP(1)          !!TOA amt
      IF (rxpt > pMax100) r = raYgivenP(kMaxLayer)  !!GND amt
      WRITE (kStdWarn,*) 'gasID,layer,ref amount (reset to) ',iGasID,iI,r
      WRITE (kStdErr,*)  'gasID,layer,ref amount (reset to) ',iGasID,iI,r
      raRAmt(iI) = r
!do iJ = 1,kMaxLayer
!  print *,'moolah',iGas,iJ,raXgivenP(iJ),raYgivenP(iJ),raaPress(iJ,iGas)
!end do
!call dostop
      
    END IF
  END IF
!^^^^^^^^^^^^^ check to see if linear interp will help^^^^^^^^^^^^^^^^^^^
END DO

!!! raaPress in in ATM, raPressLevels is in MB
!      DO iI = 1,kProfLayer
!       print *,'xaxaxa',iI,raaPress(iI,1)*1013.25,raPressLevels(iI),raaPress(iI,1)*1013.25/raPressLevels(iI)
!      END DO
!      DO iI = 1,kProfLayer
!       print *,'xaxaxa2',iI,raaPress(iI,1),raR100Press(iI),raPressLevels(iI)/1013.25
!      END DO
!      call dostopmesg('xaxaxa$')

! now just happily spline everything on!!!!!! for the temps
DO iI = 1,kMaxLayer
  raXgivenP(iI) = LOG(raR100Press(kMaxLayer-iI+1))
  raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
  raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
END DO
CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
DO iI = 1,iNot
  rxpt = LOG(raaPress(iI,iGas))
  rxpt = raaPress(iI,iGas)
  r = 273.15
  raRTemp(iI) = r
END DO
DO iI = iNot+1,kProfLayer
  rxpt = LOG(raaPress(iI,iGas))
  rxpt = raaPress(iI,iGas)
  IF (iSplineType == +1) THEN
    CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
  ELSE
    CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
  END IF
  raRTemp(iI) = r
END DO

! now just happily spline everything on!!!!!! for the partial pressures
! since only WATER cares about partial pressures, simple interpolation is fine
! for all gases except gasID = 1
!  **************** this is the orig code ********************************
IF (iGasID >= 1) THEN  !!! do for all gases, then redo Water correctly
  DO iI = 1,kMaxLayer
    raXgivenP(iI) = LOG(raR100Press(kMaxLayer-iI+1))
    raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
    raYgivenP(iI) = raR100PartPress(kMaxLayer-iI+1)
  END DO
  CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
  DO iI = 1,iNot
    rxpt = LOG(raaPress(iI,iGas))
    rxpt = raaPress(iI,iGas)
    r = 0.0
    raRPartPress(iI) = r
  END DO
  DO iI = iNot+1,kProfLayer
    rxpt = LOG(raaPress(iI,iGas))
    rxpt = raaPress(iI,iGas)
    IF (iSplineType == +1) THEN
      CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
    ELSE
      CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
    END IF
    IF ((rxpt > 0.0) .AND. (r < 0.0) .AND. (iGasID == 1)) THEN
      WRITE (kStdErr,*) 'In creating water reference profile, negative'
      WRITE (kStdErr,*) 'gas partial pressure found!!!!'
      WRITE (kStdErr,*) 'failure in "orig" test!!!!'
      WRITE (kStdErr,*) 'iSplineType  =',iSplineType
      DO iJ = 1,kMaxLayer
        PRINT *,iGas,iJ,raXgivenP(iJ),raYgivenP(iJ),rxpt
      END DO
      CALL dostop
      r2 = 0.0
      IF (iSplineType == 1) THEN
!!try linear
        CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r2)
        PRINT *,'bad spline for gas .... ',iSplineType,rxpt,r,' with linear',r2
      END IF
      WRITE (kStdErr,*) iGas,'<',iI,'>      <<',rxpt,'>>',r,r2
      WRITE (kStdErr,*) 'looping thru makerefprof spline .... iJ InPP InRefGasAmt GridPP'
      DO iJ = 1,kMaxLayer
        WRITE(kStdErr,*) iJ,raXgivenP(iJ),raYgivenP(iJ),raaPress(iJ,iGas)
      END DO
!!! this is bizarre .. if r2 > 0 and r < 0, just set r = r2 instead
!!! or CALL DoStop
      CALL DoStopMesg('wierd problems in water : MakeRefProf$')
    ELSE IF ((r < 0.0) .AND. (iGasID /= 1)) THEN
      WRITE (kStdWarn,*) 'Warning!!In creating ref profile,negative'
      WRITE (kStdWarn,*) 'gas partial pressure found!!!! Reset to 0'
      WRITE (kStdWarn,*) 'Gas ID, layer, PP = ',iGasID,iI,r
      r = 0.0
    END IF
    raRPartPress(iI) = r
  END DO
END IF
!  **************** this is the new code ********************************
! WATER cares about partial pressures, so be careful!
IF ((iGasID == 1) .AND. (iSplineType > 0)) THEN
  DO iI = 1,iNot
    r = 0.0
    raRPartPress(iI) = r
  END DO
  DO iI = iNot+1,kProfLayer
    r = raRPartPress(iI)
    r0 = raRPartPress(iI)
    CALL PPThruLayers(  &
        iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press,  &
        raDataBaseThickness,raSortPressLevels,raSortPressHeights,  &
        raPressLevels,raRTemp,r)
!-->> Howard Motteler does not bother with this, and so he sets r = r0!!!!
!-->>        print *,iGasID,iI,r,r0,r-r0
!-->>        r = r0
    IF ((r < 0.0) .AND. (iGasID == 1)) THEN
      WRITE (kStdErr,*) 'In creating water reference profile, negative'
      WRITE (kStdErr,*) 'gas partial pressure found!!!!'
      WRITE (kStdErr,*) 'failure in "new" test!!!!'
      WRITE (kStdErr,*) 'iSplineType  =',iSplineType
      CALL DoStop
    ELSE IF ((r < 0.0) .AND. (iGasID /= 1)) THEN
      WRITE (kStdWarn,*) 'Warning!!In creating ref profile,negative'
      WRITE (kStdWarn,*) 'gas partial pressure found!!!! Reset to 0'
      WRITE (kStdWarn,*) 'Gas ID, layer, PP = ',iGasID,iI,r
      r = 0.0
    END IF
    raRPartPress(iI) = r
  END DO
END IF

! simply put in the pressures
DO iI = 1,iNot
!these are "junk"
  raRPress(iI) = raaPress(iNot+1,iGas)
END DO
DO iI = iNot+1,kProfLayer
  raRPress(iI) = raaPress(iI,iGas)
END DO

RETURN
END SUBROUTINE MakeRefProf
!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

!****************
! this is for LAYER TEMPERATURE varying exponentially across layer
! since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
! this VARIES the satellite view angle as it goes through the layers
!   ie does "radiative transfer for instrument"
!****************

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth

SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(  &
    iVaryIN,raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


INTEGER, INTENT(IN OUT)                  :: iVaryIN
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
INTEGER :: iDumpAllUARads
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN,iDownWard
INTEGER :: iCloudLayerTop,iCloudLayerBot

REAL :: TEMP(MAXNZ),ravt2(maxnz)

iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
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

rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
rCos=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

WRITE(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

!cc      IF ((iNLTEStart .LE. kProfLayer) .AND. (iDoSolar .GE. 0)) THEN
!cc        DO iLay = iNumLayer,iNumLayer
!cc          iL = iaRadLayer(iLay)
!cc          IF (iL .NE. kProfLayer) THEN
!cc            write(kStdErr,*) 'NLTE rad code assumes TOA = kProfLayer'
!cc            write(kStdErr,*) 'but you seem to imply aircraft instrument'
!cc            write(kStdErr,*) 'that is NOT at TOA!!!!'
!cc            CALL DoStop
!cc          ELSE
!cc            DO iFr = 1,kMaxPts
!cc              raaLay2Sp(iFr,iL) = raaAbs(iFr,iL)
!cc            END DO
!cc          END IF
!cc        END DO
!cc 777    CONTINUE
!cc        DO iLay = iNumLayer-1,1,-1
!cc          iL = iaRadLayer(iLay)
!cc          DO iFr = 1,kMaxPts
!cc            raaLay2Sp(iFr,iL) = raaLay2Sp(iFr,iL+1)+raaAbs(iFr,iL)
!cc          END DO
!cc        END DO
!cc      END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

iVary = +1
!!!do default stuff; set temperatures at layers
IF (iVary == +1) THEN
  DO iLay=1,kProfLayer
    raVT2(iLay) = raVTemp(iLay)
  END DO
  iL = iaRadLayer(iNumLayer)
  raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
  iL = iaRadLayer(1)
  raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
  raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
iDownward = +1
CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,  &
    iProfileLayers,raPressLevels)
CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)
!      DO iFr = 1,kProflayer
!        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr),
!     $          iNLTEStart,raaPlanckCoeff(1,iFr),raaPlanckCoeff(5001,iFr)
!        end do
!      print *,'stopping here'
!      call dostop

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
DO iLay=1,iNumLayer
  iL = iaRadLayer(iLay)
! first get the Mixed Path temperature for this radiating layer
  rMPTemp = raVT1(iL)
  IF (iL < iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iL >= iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
      raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
!c        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
!c          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
!c          rOmegaSun = kOmegaSun
!c          DO iFr=1,kMaxPts
!c            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
!c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
!c            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) +
!c    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
!c    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
!c            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
!c          END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

DO iFr=1,kMaxPts
  raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  IF ((iVaryIN == 2) .OR. (iVaryIN == 3)) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,  &
        rCos,rFracBot, iVaryIN,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,-1,raInten)
  END IF
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
  IF ((iVaryIN == 2) .OR. (iVaryIN == 3)) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,  &
        rCos,+1.0, iVaryIN,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,-1,raInten)
  END IF
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
DO iLay = iHigh,iHigh
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
  IF (iUpper >= 1) THEN
!!! need to compute stuff at extra layers (100-200 km)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp >= 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      WRITE(kStdWarn,*) 'assume you need to output rad at TOA'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output everything at the top of this'
      WRITE(kStdWarn,*) 'stratosphere'
!do radiative transfer thru this layer, but do not output here
      DO iFr=1,kMaxPts
        raInten(iFr) =  &
            raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
!!! forget about interpolation thru the layers, just dump out the
!!! radiance at the top of startosphere (120-200 km)
      DO iFr=1,iDp
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END DO
    END IF
  END IF
  
  IF (iUpper < 1) THEN
!!! no need to compute stuff at extra layers (100-200 km)
!!! so just do usual stuff
!!! see if this mixed path layer is in the list iaOp to be output
!!! since we might have to do fractions!
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
  END IF
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_V

!************************************************************************
