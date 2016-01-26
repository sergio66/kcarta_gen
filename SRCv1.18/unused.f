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
c this subroutine will take in 100 AIRS layering stuff and interpolate to
c the new arbitrary layering
c WARNING : this assumes that the user has not mucked up KLAYERS layering
c           such that highest Z pressure (lowest pressure) is NOT TOA
c           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
c do the lower atm (usual -1) or upper atm (NLTE +1)

c kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
c so at arbitrary pressure layering, it deals with abs coeffs
c so we do not need raRamt
c but we do need the interpolated temp and partial pressures
      SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $           raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $           raaPress,iGas,iGasID,iNumLayers,
     $           raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
c  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
c The above equation and 3 data points define the 101 AIRS levels, which
c are in airslevels.param

c input
c do the lower atm (usual -1) or upper atm (NLTE +1)
      INTEGER iLowerOrUpper
c these are the individual reference profiles, at kMaxLayer layers
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
c these are the arbitrary profiles stored in matrices
      REAL raaPress(kProfLayer,kGasStore)
      INTEGER iError,iGas,iGasID,iNumLayers,iSplineType
c these are the kLAYERS pressure levels, layer thick for the current profile
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
c  output
c these are the individual reference profiles, at kProfLayer layers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      REAL pMax100,pMin100

c local variables
      INTEGER iI,iNot,kMaxLayer_over_2,iJ
      REAL raWorkP(kMaxLayer),raXgivenP(kMaxLayer),
     $     raYgivenP(kMaxLayer),raY2P(kMaxLayer)
      REAL raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2
      REAL raSortPressLevels(kMaxLayer+1)
      REAL raSortPressHeights(kMaxLayer+1)

      REAL raDataBaseThickness(kMaxLayer)
      REAL DatabaseHeight(kMaxLayer)
      REAL DATABASELEVHEIGHTS(kMaxLayer+1)
      REAL DATABASELEV(kMaxLayer+1)

c pressure variables!!!!! ----------------->
c raaPress in atm

      CALL databasestuff(iLowerOrUpper,
     $                   DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

C     Assign values for interpolation
C     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
      rYP1 = 1.0E+16
      rYPN = 1.0E+16

      !!! find the thickness of the pressure layers in centimeters
      !!! hence we have x 1000 (km -> m) x 100 (m -> cm)
      DO iI = 1,kMaxLayer
        raDataBaseThickness(iI) = (DATABASELEVHEIGHTS(iI+1) - 
     $                            DATABASELEVHEIGHTS(iI))*1000.0 * 100.0
      END DO

c now store stuff sorted in terms of increasing pressure 
      kMaxLayer_over_2 = 50
      DO iI = 1,kMaxLayer+1
        raSortPressLevels(iI)  = DataBaseLev(101-iI+1)                  !!! mb
        raSortPressHeights(iI) = DataBaseLevHeights(101-iI+1) * 1000.0  !!! m
      END DO

      !!!this tells how many layers are NOT dumped out by kLAYERS
      iNot = kProfLayer-iNumLayers
c now just happily spline everything on!!!!!! for the amts
c recall you need raXgivenP to be increasing
      pMax100 = -1.0e10
      pMin100 = +1.0e10
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100Amt(kMaxLayer-iI+1)
        IF (raXgivenP(iI) .LT. pMin100) THEN
          pMin100 = raXgivenP(iI)
        END IF
        IF (raXgivenP(iI) .GT. pMax100) THEN
          pMax100 = raXgivenP(iI)
        END IF
      END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        r = 0.0
        raRAmt(iI) = r
      END DO

c these are the layers with info
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        IF (iSplineType .EQ. +1) THEN
          CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        ELSE
          CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rXPT,r)
        END IF
        raRAmt(iI) = r
c^^^^^^^^^^^^^ check to see if linear interp will help^^^^^^^^^^^^^^^^^^^
        IF ((r .LT. 0.0)) THEN
          write (kStdWarn,*) 'Making reference profile : negative amt found : '
          write (kStdWarn,*) 'GasID,layer,ref amt,iLorU (-1/+1) : ',iGasID,iI,r,iLowerOrUpper
          CALL rlinear(raXgivenP,raYgivenP,kMaxLayer,rxpt,r,1)
          write (kStdWarn,*) '   trying linear interp : ',iI,r
          raRAmt(iI) = r          
          IF ( (r .LT. 0.0). AND. 
     $         ((rxpt .GT. pMin100) .AND. (rxpt .LT. pMax100))) THEN
              !!!things barfed even though pMin100 < rXpt < pMax100 ... bad
            write (kStdErr,*) 'gasID,layer,ref amount(linear interp)= ',
     $ iGasID,iI,r
            write (kStdErr,*) 'min(pD),rXpt,max(pD) = ',pMin100,rXpt,pMax100
            CALL DoStop
          ELSEIF ( (r .LT. 0.0). AND. 
     $          ((rxpt .LT. pMin100) .OR. (rxpt .GT. pMax100))) THEN
              !!!things barfed, pMin100 > rXpt or pMax100 < rXpt ... hmmm
            write (kStdWarn,*) 'gasID,layer,ref amount(linear interp)= ',
     $ iGasID,iI,r
            !!!recall raYgiven(iI) is swtiched GND = kMaxLayer,100=TOA
            IF (rxpt .LT. pMin100) r = raYgivenP(1)          !!TOA amt
            IF (rxpt .GT. pMax100) r = raYgivenP(kMaxLayer)  !!GND amt
            write (kStdWarn,*) 'gasID,layer,ref amount (reset to) ',iGasID,iI,r
            write (kStdErr,*)  'gasID,layer,ref amount (reset to) ',iGasID,iI,r	    
            raRAmt(iI) = r
            !do iJ = 1,kMaxLayer
            !  print *,'moolah',iGas,iJ,raXgivenP(iJ),raYgivenP(iJ),raaPress(iJ,iGas)
            !end do
            !call dostop

          END IF
        END IF
c^^^^^^^^^^^^^ check to see if linear interp will help^^^^^^^^^^^^^^^^^^^
      END DO

      !!! raaPress in in ATM, raPressLevels is in MB
c      DO iI = 1,kProfLayer
c       print *,'xaxaxa',iI,raaPress(iI,1)*1013.25,raPressLevels(iI),raaPress(iI,1)*1013.25/raPressLevels(iI)
c      END DO
c      DO iI = 1,kProfLayer
c       print *,'xaxaxa2',iI,raaPress(iI,1),raR100Press(iI),raPressLevels(iI)/1013.25
c      END DO
c      call dostopmesg('xaxaxa$')
     
c now just happily spline everything on!!!!!! for the temps
      DO iI = 1,kMaxLayer
        raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
        raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
        raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
      END DO
      CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
      DO iI = 1,iNot
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        r = 273.15
        raRTemp(iI) = r
      END DO
      DO iI = iNot+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        IF (iSplineType .EQ. +1) THEN
          CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        ELSE
          CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
        END IF
        raRTemp(iI) = r
      END DO

c now just happily spline everything on!!!!!! for the partial pressures
c since only WATER cares about partial pressures, simple interpolation is fine 
c for all gases except gasID = 1
c  **************** this is the orig code ********************************
      IF (iGasID .GE. 1) THEN  !!! do for all gases, then redo Water correctly
        DO iI = 1,kMaxLayer
          raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
          raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
          raYgivenP(iI) = raR100PartPress(kMaxLayer-iI+1)
        END DO
        CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
        DO iI = 1,iNot
          rxpt = log(raaPress(iI,iGas))
          rxpt = raaPress(iI,iGas)
          r = 0.0
          raRPartPress(iI) = r
        END DO
        DO iI = iNot+1,kProfLayer
          rxpt = log(raaPress(iI,iGas))
          rxpt = raaPress(iI,iGas)
          IF (iSplineType .EQ. +1) THEN
            CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
          ELSE
            CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
          END IF
          IF ((rxpt .GT. 0.0) .AND. (r .LT. 0.0) .AND. (iGasID .EQ. 1)) THEN
            write (kStdErr,*) 'In creating water reference profile, negative'
            write (kStdErr,*) 'gas partial pressure found!!!!'
            write (kStdErr,*) 'failure in "orig" test!!!!'
            write (kStdErr,*) 'iSplineType  =',iSplineType
	    DO iJ = 1,kMaxLayer
	      print *,iGas,iJ,raXgivenP(iJ),raYgivenP(iJ),rxpt
	    END DO
	    call dostop
            r2 = 0.0
            IF (iSplineType .EQ. 1) THEN
              !!try linear
              CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r2)
              print *,'bad spline for gas .... ',iSplineType,rxpt,r,' with linear',r2
            END IF
            write (kStdErr,*) iGas,'<',iI,'>      <<',rxpt,'>>',r,r2
            write (kStdErr,*) 'looping thru makerefprof spline .... iJ InPP InRefGasAmt GridPP'
            DO iJ = 1,kMaxLayer
              write(kStdErr,*) iJ,raXgivenP(iJ),raYgivenP(iJ),raaPress(iJ,iGas)
            END DO   
            !!! this is bizarre .. if r2 > 0 and r < 0, just set r = r2 instead
            !!! or CALL DoStop         
            CALL DoStopMesg('wierd problems in water : MakeRefProf$')
          ELSEIF ((r .LT. 0.0) .AND. (iGasID .NE. 1)) THEN
            write (kStdWarn,*) 'Warning!!In creating ref profile,negative'
            write (kStdWarn,*) 'gas partial pressure found!!!! Reset to 0'
            write (kStdWarn,*) 'Gas ID, layer, PP = ',iGasID,iI,r
            r = 0.0
          END IF
          raRPartPress(iI) = r
        END DO
      END IF
c  **************** this is the new code ********************************
c WATER cares about partial pressures, so be careful!
      IF ((iGasID .EQ. 1) .AND. (iSplineType .GT. 0)) THEN 
        DO iI = 1,iNot
          r = 0.0
          raRPartPress(iI) = r
        END DO
        DO iI = iNot+1,kProfLayer
          r = raRPartPress(iI)
          r0 = raRPartPress(iI)
          CALL PPThruLayers(
     $     iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press,
     $     raDataBaseThickness,raSortPressLevels,raSortPressHeights,
     $     raPressLevels,raRTemp,r)
c-->> Howard Motteler does not bother with this, and so he sets r = r0!!!! 
c-->>        print *,iGasID,iI,r,r0,r-r0
c-->>        r = r0
          IF ((r .LT. 0.0) .AND. (iGasID .EQ. 1)) THEN
            write (kStdErr,*) 'In creating water reference profile, negative'
            write (kStdErr,*) 'gas partial pressure found!!!!'
            write (kStdErr,*) 'failure in "new" test!!!!'
            write (kStdErr,*) 'iSplineType  =',iSplineType
            CALL DoStop
          ELSEIF ((r .LT. 0.0) .AND. (iGasID .NE. 1)) THEN
            write (kStdWarn,*) 'Warning!!In creating ref profile,negative'
            write (kStdWarn,*) 'gas partial pressure found!!!! Reset to 0'
            write (kStdWarn,*) 'Gas ID, layer, PP = ',iGasID,iI,r
            r = 0.0
          END IF
          raRPartPress(iI) = r
        END DO
      END IF

c simply put in the pressures
      DO iI = 1,iNot
        !these are "junk"
        raRPress(iI) = raaPress(iNot+1,iGas)
      END DO
      DO iI = iNot+1,kProfLayer
        raRPress(iI) = raaPress(iI,iGas)
      END DO

      RETURN
      END
c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c****************
c this is for LAYER TEMPERATURE varying exponentially across layer
c since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
c this VARIES the satellite view angle as it goes through the layers
c   ie does "radiative transfer for instrument"
c****************

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(
     $    iVaryIN,raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iDumpAllUARads = dump rads at all layers or only select layers?
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iVaryIN,iDumpAllUARads
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iDownWard
      INTEGER iCloudLayerTop,iCloudLayerBot

      REAL TEMP(MAXNZ),ravt2(maxnz)

      iIOUN = iIOUN_IN

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
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

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

      write(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

ccc      IF ((iNLTEStart .LE. kProfLayer) .AND. (iDoSolar .GE. 0)) THEN
ccc        DO iLay = iNumLayer,iNumLayer
ccc          iL = iaRadLayer(iLay)
ccc          IF (iL .NE. kProfLayer) THEN
ccc            write(kStdErr,*) 'NLTE rad code assumes TOA = kProfLayer'
ccc            write(kStdErr,*) 'but you seem to imply aircraft instrument'
ccc            write(kStdErr,*) 'that is NOT at TOA!!!!'
ccc            CALL DoStop
ccc          ELSE
ccc            DO iFr = 1,kMaxPts
ccc              raaLay2Sp(iFr,iL) = raaAbs(iFr,iL)
ccc            END DO
ccc          END IF
ccc        END DO
ccc 777    CONTINUE
ccc        DO iLay = iNumLayer-1,1,-1
ccc          iL = iaRadLayer(iLay)
ccc          DO iFr = 1,kMaxPts
ccc            raaLay2Sp(iFr,iL) = raaLay2Sp(iFr,iL+1)+raaAbs(iFr,iL)
ccc          END DO
ccc        END DO
ccc      END IF

c note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      iVary = +1
      !!!do default stuff; set temperatures at layers
      IF (iVary .EQ. +1) THEN
        DO iLay=1,kProfLayer
          raVT2(iLay) = raVTemp(iLay)
        END DO
        iL = iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
        iL = iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
      END IF

c set the vertical temperatures of the atmosphere 
c temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
      iDownward = +1
      CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,
     $                   iProfileLayers,raPressLevels)
      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)
c      DO iFr = 1,kProflayer
c        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr),
c     $          iNLTEStart,raaPlanckCoeff(1,iFr),raaPlanckCoeff(5001,iFr)
c        end do
c      print *,'stopping here'
c      call dostop

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
           raaEmission(iFr,iLay)=0.0
         END DO
       END DO
       DO iLay=2,iNumLayer-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)/rCos)
           raaEmission(iFr,iLay)=0.0
         END DO
       END DO
       DO iLay = iNumLayer,iNumLayer
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
           raaEmission(iFr,iLay)=0.0
         END DO
       END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
        ELSEIF (iL .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)	  
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
cc        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
cc          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
cc          rOmegaSun = kOmegaSun
cc          DO iFr=1,kMaxPts
cc            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
cc            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
cc            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) + 
cc    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
cc    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
cc            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
cc          END DO
        END IF
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $         raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $         iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)

      DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        IF ((iVaryIN .EQ. 2) .OR. (iVaryIN .EQ. 3)) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,rFracBot,
     $                      iVaryIN,raInten)
        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,-1,raInten)
        END IF
      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

        IF ((iVaryIN .EQ. 2) .OR. (iVaryIN .EQ. 3)) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,+1.0,
     $                      iVaryIN,raInten)
        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,-1,raInten)
        END IF
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            write(kStdWarn,*) 'assume you need to output rad at TOA'
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
            write(kStdWarn,*) 'and output everything at the top of this'
            write(kStdWarn,*) 'stratosphere'
            !do radiative transfer thru this layer, but do not output here
            DO iFr=1,kMaxPts
              raInten(iFr) = 
     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
            DO iFr=1,iDp
              CALL wrtout(iIOUN,caOutName,raFreq,raInten)
            END DO
          END IF
        END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
           IF (iDp .GT. 0) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
             DO iFr=1,iDp
               CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
               CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
        END IF

      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
