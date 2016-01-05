c this is a much more detailed LBLRTM TAPE5 reader
c and hmm only seems to work for case mentioned below??? and is not called anywhere else
c see https://svn.ssec.wisc.edu/repos/uwphysret/trunk/mfiles/iasi_radpress_write5.m
      SUBROUTINE  ReadInput_DetailedLBLRTM_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,raP,raT,
     $                          iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)

      IMPLICIT NONE
     
      INCLUDE '../INCLUDE/kcarta.param'

c input
      CHARACTER*80 caPfName
      REAL rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
c output
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas)

c local var
      INTEGER iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
      REAL raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH
      CHARACTER*80 caStr,caStrX,caStrY
      CHARACTER*30 caStr30
      CHARACTER*1  c1

      rPmin = +1.0e6
      rPmax = -1.0e+6
      iNXsec = -1

 111  FORMAT(A80)
 
      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
 1070     FORMAT('ERROR! number ',I5,' opening PRFILE path LBLRTM profile file'
     $            ,/,A80)
          CALL DoSTOP
      ENDIF
      kProfileUnitOpen=1

c this should read /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
      READ (iIOUN2,111,ERR=13,END=13) caStr
      READ (iIOUN2,111,ERR=13,END=13) caStr

      READ (iIOUN2,*) rF1,rF2
      write(kStdWarn,*) 'LBLRTM input indicates start/stop wavenumbers are ',rF1,rF2

      READ (iIOUN2,*) rTSurf
      write(kStdWarn,*) 'LBLRTM input indicates STEMP is ',rTSurf

      READ (iIOUN2,*) (iaJunk(iL),iL = 1,9)
      iNumLevs = iaJunk(3)
      iNumGases = iaJunk(6)
      IF (iNumLevs .LT. 0) THEN
        write(kStdWarn,*) 'looks like boundaries for calculations specified in mb'
        iNumLevs = abs(iNumLevs)
        iLBROutBdryHorP = -1
      ELSE
        iLBROutBdryHorP = +1
        write(kStdWarn,*) 'looks like boundaries for calculations specified in km'
      END IF

      write(kStdWarn,*) 'TAPE5 implies kCARTA compute ODs at ',iNumLevs,' press/heights, for iNumGases = ',iNumGases
      IF (iNumLevs .GT. kProfLayer) THEN
        write(kStdErr,*) 'though irrelevant for kCARTA calcs, it will not be able to read in the pressures/heights'
        write(kStdErr,*) iNumLevs,kProfLayer
        CALL DoStop
      END IF
      IF (iNumGases .GT. kMaxGas) THEN
        write(kStdErr,*) 'iNumGases .GT. kMaxGas',iNumGases,kMaxGas
        CALL DoStop
      END IF

      READ (iIOUN2,*) rTophgt,rHSurf,rViewAngle
      write(kStdWarn,*) 'LBLRTM input indicates start/stop hgts are ',rHSUrf,rTopHgt,' with view angle ',rViewAngle

      READ (iIOUN2,*) (raLBL_Hgts(iJ),iJ=1,iNumLevs)

      raRTP_TxtInput(6) = +raLBL_Hgts(iNumLevs)    !! default, assume hgt in km
      IF  (iLBROutBdryHorP .LT. 0) raRTP_TxtInput(6) = -raLBL_Hgts(iNumLevs)

      READ (iIOUN2,*) iNumLevs
      IF (iNumLevs .GT. 2*kProfLayer) THEN
        write(kStdErr,*) 'iNumLevs .GT. 2*kProfLayer',iNumLevs,kProfLayer
        CALL DoStop
      END IF

c see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
c       JCHAR = 1-6           - default to value for specified model atmosphere
c              = " ",A         - volume mixing ratio (ppmv)
c              = B             - number density (cm-3)
c              = C             - mass mixing ratio (gm/kg)
c              = D             - mass density (gm m-3)
c              = E             - partial pressure (mb)
c              = F             - dew point temp (K) *H2O only*
c              = G             - dew point temp (C) *H2O only*
c              = H             - relative humidity (percent) *H2O only*
c              = I             - available for user definition
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

        IF (iL .EQ. 1) THEN
          DO iJ=1,iNumGases
            c1 = caStr30(iJ:iJ)
            IF (c1 .EQ. 'A') iaGasUnits(iJ) = 10
            IF (c1 .EQ. ' ') iaGasUnits(iJ) = 10
            IF (c1 .EQ. 'B') iaGasUnits(iJ) = -1
            IF (c1 .EQ. 'C') iaGasUnits(iJ) = 20
            IF (c1 .EQ. 'D') iaGasUnits(iJ) = -1
            IF (c1 .EQ. 'E') iaGasUnits(iJ) = -1
            IF (c1 .EQ. 'F') iaGasUnits(iJ) = 42
            IF (c1 .EQ. 'G') iaGasUnits(iJ) = 43             
            IF (c1 .EQ. 'H') iaGasUnits(iJ) = 40
            IF (c1 .EQ. 'I') iaGasUnits(iJ) = -1
            
            IF (iaGasUnits(iJ) .LT. 0) THEN
              write(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
              CALL DOStop
            ELSE
              write(kStdWarn,*) 'LBLRTM gas number ',iJ,' "units" ',c1,' = kCARTA levels units of ',iaGasUnits(iJ)
            END IF
          END DO
        END IF

        READ (iIOUN2,*) (raX(iG),iG=1,iNumGases)
c        raP(iL) = rP * 100.0  !! change from mb to N/m2
        raP(iL) = rP          !! keep in mb
        raT(iL) = rT
        IF (rPmax .LE. raP(iL)) rPmax = raP(iL)
        IF (rPmin .GT. raP(iL)) rPmin = raP(iL)
        DO iJ = 1,iNumGases
          raaG_MR(iL,iJ) = raX(iJ)
        END DO

        IF (iL .EQ. 1) THEN
c          rPSurf = raP(1)/100 !! because we redo this below
          rPSurf = raP(1)     !! keep in mb        
          IF ((rHSurf .GT. rHmaxKCarta) .OR. (rHSurf .LT. rHminKCarta)) THEN
            write(kStdErr,*) 'need rHmaxKCarta >= rHSurf >= rHminKCarta but have'
            write(kStdErr,*) '(rHmaxKCarta,rHSurf,rHminKCarta) = ',rHmaxKCarta,rHSurf,rHminKCarta
            CALL DoStop
          END IF
          IF ((rPSurf .GT. rPmaxKCarta) .OR. (rPSurf .LT. rPminKCarta)) THEN
            write(kStdErr,*) 'need rPmaxKCarta >= rPSurf >= rPminKCarta but have'
            write(kStdErr,*) '(rPmaxKCarta,rPSurf,rPminKCarta) = ',rPmaxKCarta,rPSurf,rPminKCarta
            CALL DoStop
          END IF
          IF ((rTSurf .GT. kStempMax) .OR. (rTSurf .LT. kStempMin)) THEN
            write(kStdErr,*) 'need kStempMax >= rTSurf >= kStempMin but have'
            write(kStdErr,*) '(kStempMax,rTSurf,kStempMin) = ',kStempMax,rTSurf,kStempMin
            CALL DoStop
          END IF
        END IF
      END DO

      !! now see if there are xsec gases
      READ (iIOUN2,111,ERR=13,END=13) caStr
      READ(caStr,*) iNXsec
      IF (iNXsec .GT. 0) THEN
        READ (iIOUN2,111,ERR=13,END=13) caStr     !!!! xsec names 
        CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
        READ (iIOUN2,*) iNumLevsXsec
        IF (iNumLevsXsec .GT. 2*kProfLayer) THEN
          write(kStdErr,*) 'iNumLevsXsec .GT. 2*kProfLayer',iNumLevsXsec,kProfLayer
          CALL DoStop
        END IF
        IF (iNumLevsXsec .NE. iNumLevs) THEN
          write(kStdErr,*) 'iNumLevsXsec .NE. iNumLevs',iNumLevsXsec,iNumLevs
          CALL DoStop
        END IF
        DO iL = 1,iNumLevsXsec
          READ (iIOUN2,111,ERR=13,END=13) caStrY
          caStr30 = caStrY(16:45)
          IF (iL .EQ. 1) THEN
            DO iJ=1,iNXsec
              c1 = caStr30(iJ:iJ)
              IF (c1 .EQ. 'A') iaGasUnits(iJ+iNumGases) = 10
              IF (c1 .EQ. ' ') iaGasUnits(iJ+iNumGases) = 10
              IF (c1 .EQ. 'B') iaGasUnits(iJ+iNumGases) = -1
              IF (c1 .EQ. 'C') iaGasUnits(iJ+iNumGases) = 20
              IF (c1 .EQ. 'D') iaGasUnits(iJ+iNumGases) = -1
              IF (c1 .EQ. 'E') iaGasUnits(iJ+iNumGases) = -1
              IF (c1 .EQ. 'F') iaGasUnits(iJ+iNumGases) = 42
              IF (c1 .EQ. 'G') iaGasUnits(iJ+iNumGases) = 43             
              IF (c1 .EQ. 'H') iaGasUnits(iJ+iNumGases) = 40
              IF (c1 .EQ. 'I') iaGasUnits(iJ+iNumGases) = -1
              
              IF (iaGasUnits(iJ+iNumGases) .LT. 0) THEN
                write(kStdErr,*) 'LBLRTM --> kCarta not set up to deal with gas units ',c1,' for gas number ',iJ
                CALL DOStop
              ELSE
                write(kStdWarn,*) 'LBLRTM xsec number ',iJ,' is of "units" ',c1,' = kCARTA input levels units of ',
     $              iaGasUnits(iJ+iNumGases)
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

      write(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
      write(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
      write(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf
     
c make sure pressures are decreasing with index ie layers going higher and higher
      IF (raP(1) .LT. raP(2)) THEN
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
      write(kStdWarn,*) 'input file : Highest altitude (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
      write(kStdWarn,*) ' '

c now change all units to MR
      DO iG = 1,iNumGases
        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR)
        DO iL = 1,iNumLevs
          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
          iaGasUnits(iG) = 12
        END DO
      END DO         

      RETURN
      END

c************************************************************************
c this subroutine does the klayers integrations, integrals over P
c so this is more similar to kLAYERs, CRTM
c this is orig code; it puts surface as "lowest level and adds the user level info above that
      SUBROUTINE DoIntegrateLevels2Layers_wrtP_pSurf(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,
     $               PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,
     $               raPX,raTX,raaG_MRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL rHSurf,rPSurf,rTSurf,PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      REAL DATABASELEVHEIGHTS(kMaxLayer+1),rfracBot,rLat,rLon
      INTEGER iLowestLev,iNumGases,iaG(kMaxGas)
      REAL raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas)
c output
      REAL raTout(kProfLayer),raAmountOut(kProfLayer),raZout(kProfLayer+1),raPout(kProfLayer)
      REAL raaQout(kProfLayer,kMaxGas),raaPartPressOut(kProfLayer,kMaxGas)

c local
      INTEGER iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
      REAL z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
      REAL rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
      REAL dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
      REAL raXYZPress(kMaxlayer+1+1)
      REAL raXYZTemp(kMaxlayer+1+1)
      REAL raaXYZ_MR(kMaxlayer+1+1,kMaxGas)
      REAL raXYZ_MRwater(kMaxlayer+1+1)
      REAL raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1
      REAL zWoo,rConvertQ,grav_earth,wexsvp,rRH
      INTEGER iWoo

      rConvertQ = 6.022141e23/1e4     ! multiply by this to change moles/m2 to molecules/cm2
      rConvertQ = kAvog/1000/1.0e4    ! multiply by this to change moles/m2 to molecules/cm2

      DO iL = 1,kProfLayer
        DO iG = 1,iNumGases
          raaQout(iL,iG) = 0.0
        END DO
      END DO

      IF ((kPlanet .EQ. 3) .AND. (iaG(1) .NE. 1)) THEN
        write (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
        Call DoStop
      END IF

      iNFine = 10
      iNFine = 200
      iNFine = 50

c >>>>>>>>>>
      !! set up temporary arrays/matrices for interping onto sublevels
      !! set first level at surface press, surface temp
      raXYZPress(1) = rPSurf
      raXYZTemp(1)  = rTSurf
      iUpperLev = kMaxLayer-iLowestLev+1
      IF ((iUpperLev .GT. kMaxLayer) .OR. (iUpperLev .LT. 1)) THEN
        write (kStdErr,*) '(iUpperLev .GT. kMaxLayer) .OR. (iUpperLev .LT. 1)',iUpperLev,kMaxLayer
        Call Dostop
      END IF
      !! do surface by interpolating MR from user supplied/kcarta_database levels, to surface pressure
      DO iG = 1,iNumGases
        DO iL = 1,iUpperLev
          raJunk(iL) = raaG_MRX(iL,iG)
        END DO
        CALL r_sort_loglinear(raPX,raJunk,iUpperLev,rPSurf,junk,1)
        raaXYZ_MR(1,iG) = junk
        IF (iG .EQ.  1) raXYZ_MRwater(1) = junk
      END DO
c >>>>>>>>>>
      !! set levels 2 .. iNumLevs+1
      DO iL = 1,iUpperLev
        raXYZPress(iL+1) = raPX(iL)
        raXYZTemp(iL+1) = raTX(iL)
        DO iG = 1,iNumGases
          raaXYZ_MR(iL+1,iG) = raaG_MRX(iL,iG)
          IF (iG .EQ.  1) raXYZ_MRwater(iL+1) = raaG_MRX(iL,iG)
        END DO
      END DO
c >>>>>>>>>>

c now start integrating
c lowest layer is special case

c http://cimss.ssec.wisc.edu/~paulv/
c http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
c see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
c     /home/sergio/KCARTA/DOC/sci_klayers.txt

      z = rHSurf !! need to kludge this
      zWoo = rHSurf
      iWoo = -1

      DO iL = iLowestLev,iLowestLev
        raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
        raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)

        z0 = z
        iJ = iL - iLowestLev + 1

        amount = 0.0
        rTsum  = 0.0
        rWgt = 0.0

        !!! >>>>>>  this is starting out at         rP = PLEV_KCARTADATABASE_AIRS(iL)        rT = rTSurfx  >>>>>>>>
        rP = PLEV_KCARTADATABASE_AIRS(iL)
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,iUpperLev+1,rP,rT,1)
        
        dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                             !! current sublev press
        rT_n = rT                             !! current sublev temp
        rMR_water_n = raXYZ_MRwater(1)        !! current water MR
        CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP,rMR_water_n,1) !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                              !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw

        IF ((rP_n .LE. rPSurf) .AND. (iWoo .LT. 0)) THEN
          iWoo = +1
          zWoo   = z
        END IF

        DO iLoop = 1,iNFine

          rMolarMass_n   = kAtmMolarMass
          rMolarMass_np1 = kAtmMolarMass
          !!! if Earth, adjust molecular mass for displacement by water
          IF (kPlanet .EQ. 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
          IF (kPlanet .EQ. 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
          rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
          rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

          density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
          density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

          IF (kPlanet .NE. 3) THEN
            gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
            gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
          ELSE
            gravity = grav_earth(z,0.0,0.0,rLat,rLon)
          END IF

          Tav   = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)          
          SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)  
          dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

          rTsum = rTsum + Tav * (density_n + density_np1)/2
          rWgt  = rWgt + (density_n + density_np1)/2
                         
          Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
          damount = Pav/Tav/kMGC * dz 
c          print *,iLoop,PLEV_KCARTADATABASE_AIRS(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
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

          rP_np1 = log(rP_n) - dlnp
          rP_np1 = exp(rP_np1)                                                              !! next sublev pressure
          CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw

          IF ((rP_n .LE. rPSurf) .AND. (iWoo .LT. 0)) THEN
            iWoo = +1
            zWoo   = z
          END IF

        END DO

        raTout(iL)      = rTsum/rWgt                 ! <<< need TAV
        raAmountOut(iL) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
        raZout(iL+1)    = z

c        print *,iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
c     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
      END DO

      ! displace everything  z
      iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
      raZout(iL)   = -(zWoo - z0)
      raZout(iL+1) = +(z-zWoo)
      write(kStdWarn,*) '  need to displace lowest layer heights'
      write(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
      write(kStdWarn,*) '  Plowest,pSurf,pHighest (mb) ',PLEV_KCARTADATABASE_AIRS(iL),rPSurf,
     $                                                   PLEV_KCARTADATABASE_AIRS(iL+1)
      z = z - zWoo

ccc      !no need to adjust this amount for book-keeping, make consistent with eg n_rad_jac_scat.f since we did FULL layer
ccc      raAmountOut(iL) = raAmountOut(iL)/rFracBot 
ccc      DO iG = 1,iNumGases
ccc        raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
ccc      END DO
      
c go to TOA
      DO iL = iLowestLev + 1,kProfLayer
        z0 = z

        !! compute avg pressure of layer
        raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
        raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)

        iJ = iL - iLowestLev + 1
        rP = PLEV_KCARTADATABASE_AIRS(iL)
        rT = raXYZTemp(iJ)
        
        amount = 0.0
        rTsum  = 0.0
        rWgt = 0.0

        dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                              !! current sublev press
        rT_n = rT                              !! current sublev temp
        rMR_water_n = raXYZ_MRwater(iJ)        !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                              !! next sublev pressure
        CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw

        DO iLoop = 1,iNFine

          rMolarMass_n   = kAtmMolarMass
          rMolarMass_np1 = kAtmMolarMass
          !!! if Earth, adjust molecular mass for displacement by water
          IF (kPlanet .EQ. 03) rMolarMass_n   = kAtmMolarMass - rMR_water_n  *(kAtmMolarMass - 18.0)
          IF (kPlanet .EQ. 03) rMolarMass_np1 = kAtmMolarMass - rMR_water_np1*(kAtmMolarMass - 18.0)
          rMolarMass_n   = rMolarMass_n/1000.0      !! change from g/mol to kg/mol
          rMolarMass_np1 = rMolarMass_np1/1000.0    !! change from g/mol to kg/mol

          density_n   = rP_n   / rT_n   * rMolarMass_n   / kMGC
          density_np1 = rP_np1 / rT_np1 * rMolarMass_np1 / kMGC

          IF (kPlanet .NE. 3) THEN
            gravity = kGravity * (1 - 2 * z/(kPlanetRadius * 1000))
            gravity = kGravity/((1 + z/(kPlanetRadius * 1000))**2)
          ELSE
            gravity = grav_earth(z,0.0,0.0,rLat,rLon)
          END IF

          Tav = (rT_n * density_n + rT_np1 * density_np1) / (density_n + density_np1)          
          SpHeatDryAir = kMGC*(1/rMolarMass_n * density_n + 1/rMolarMass_np1 * density_np1) / (density_n + density_np1)         
          dz = SpHeatDryAir * Tav / gravity * log(rP_n/rP_np1)

          rTsum = rTsum + Tav * (density_n + density_np1)/2
          rWgt  = rWgt + (density_n + density_np1)/2
                         
          Pav = (rP_n-rP_np1)/log(rP_n/rP_np1)
          damount = Pav/Tav/kMGC * dz
          amount = amount + damount
          
          DO iG = 1,iNumGases
            DO iCnt = 1,iUpperLev+1
              raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
            END DO
            CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_n,  rMR_n,  1) 
            CALL r_sort_loglinear(raXYZPress,raJunk,iUpperLev+1,rP_np1,rMR_np1,1) 
            raaQout(iL,iG) = raaQout(iL,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
c            q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
c            q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
c            raaQout(iL,iG) = raaQout(iL,iG) + (q_n + q_np1)/2 * rConvertQ
          END DO

          !!! now update for next iteration
          z = z + dz
          rP_n = rP_np1
          rT_n = rT_np1
          rMR_water_n = rMR_water_np1
          rP = rP_n

          rP_np1 = log(rP_n) - dlnp
          rP_np1 = exp(rP_np1)                                                              !! next sublev pressure
          CALL r_sort_loglinear(raXYZPress,raXYZTemp,    iUpperLev+1,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear(raXYZPress,raXYZ_MRwater,iUpperLev+1,rP_np1,rMR_water_np1,1)   !! next sublev MRw

        END DO

        raTout(iL)      = rTsum/rWgt             ! <<< need TAV
        raAmountOut(iL) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
        raZout(iL+1)    = z

c        print *,iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
c     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
      END DO

      write(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
      !! raZout is in meters
      !! raaQout is in moles/m2
      DO iL = iLowestLev,kProfLayer
        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
        rT = raTout(iL)  
        DO iG = 1,iNumGases
          raaPartPressOut(iL,iG) = raaQout(iL,iG) / raAmountOut(iL) * raPout(iL)
        END DO
      END DO

c testing
c          !!! this automatically puts partial pressure in ATM, assuming
c          !!! gas amount in kilomoles/cm2, length in cm, T in kelvin
c          !!!note "j"!!!
c      DO iL = iLowestLev,kProfLayer
c        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
c        rT = raTout(iL)  
c        DO iG = 1,iNumGases
c          rAmt = raaQout(iL,iG)/kAvog
c          rPP  = rAmt*1.0e9*kMGC*rT / (rThick*kAtm2mb*100.0)
c          IF (iG .EQ. 1) THEN
c            print *,iL,iG,raaPartPressOut(iL,iG)/(kAtm2mb*100).0,rPP,rThick,
c     $             raaPartPressOut(iL,iG)/rPP/(kAtm2mb*100.0),raAmountOut(iL),raaQout(iL,iG)
c          END IF
c        print *,iL,raPout(iL)/100,raTout(iL),raAmountOut(iL),
c     $        raaQout(iL,1),raaQout(iL,2),raaQout(iL,3),raaQout(iL,4),raaQout(iL,5),
c     $        raaG_MRX(iL-iLowestLev+1,1),raaG_MRX(iL-iLowestLev+1,2),raaG_MRX(iL-iLowestLev+1,3),
c     $        raaG_MRX(iL-iLowestLev+1,4),raaG_MRX(iL-iLowestLev+1,5)
c        END DO
c      END DO
c testing

      IF ((iLowestLev-1) .GE. 1) THEN
        DO iL = 1,iLowestLev-1
          raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
          raPout(iL) = (PLEV_KCARTADATABASE_AIRS(iL) - PLEV_KCARTADATABASE_AIRS(iL+1))/raPout(iL)
          raTout(iL) = kStempMin
          DO iG = 1,iNumGases
            raaPartPressOut(iL,iG) = 0.0
            raaQout(iL,iG) = 0.0
            raaG_MRX(iL,iG) = 0.0
          END DO
        END DO
      END IF

      IF (kPlanet .EQ. 3) THEN
        !! show the RH for gas 1
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'Checking Relative Humidity'
        write(kStdWarn,*) ' Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      Q         Q1         Q2        Q3'
        write(kStdWarn,*) '-------------------------------------------------------------------------------------------------------'
        DO iL = iLowestLev,kProfLayer
          z    = (raZout(iL) + raZout(iL+1))/2/1000
          zWoo = (raZout(iL+1) - raZout(iL))/1000
          rP   = raPout(iL)            !! N/m2
          rT   = raTout(iL)
          rPP  = raaPartPressOut(iL,1) !! N/m2
          rSVP = wexsvp(rT) * 100      !! mb --> N/m2
          rRH  = rPP/rSVP*100.0
c          print *,iL,rP/100,rT,rPP/100,rRH
c          IF (rRH .LE. 100.0) THEN
c            write(kStdWarn,111) iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),(raaQout(iL,iG),iG=1,3)
c          ELSE
c            write(kStdWarn,112) iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),(raaQout(iL,iG),iG=1,3),' ***** '
c          END IF
        END DO
      END IF
 111  FORMAT(I3,' ',1(F9.5,' '),6(F8.4,' '),4(ES9.4,' '))
 112  FORMAT(I3,' ',1(F9.5,' '),6(F8.4,' '),4(ES9.4,' '),A7)

      RETURN
      END

c************************************************************************
