c Copyright 2014
c University of Maryland Baltimore County 
c All Rights Reserved

c this file deals with a simple layers
c does everything in terms of LINEAR interps, not SPLINE
c DO NOT USE        Call r_sort_logspl(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
c DO     USE        Call r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
c
c assumes WV is in "wet air", other gases in dry air MR
c code
c  (a) reads in input text levels profile, gases in whatever units
c      changes the input units to MIXING RATIOS
c  (b) read in upper level MR profiles for above gases, as well as for "missing gases"
c      adjusts the Ref MR to agree with input where the input levels end, and then tapers the
c        adjust ratio to 1, with about 4 AIRS levels
c      adjusts the Ref Temp to agree with input where the input levels end, and then tapers the
c        adjust ratio to 0, with about 4 AIRS levels
c   ** the tacking on (or adding new gases) can either be via
c       - using the reference P,PP,T for 100 layers ==> MR = PP/P
c       - using one of the six AFGL profiles        ==> MR comes from the ~50 levels
c  (c) if Earth do the adjustment of dry air mix ratios
c  (d) does the sublev integrations
c
c Useful References 
c   http://www.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_1/index.html
c   http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
c saved off in ../PDF/wisc_notes_*.pdf
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c this subroutine reads in a TEST file LEVELS mixing ratio PROFILE
c the text file will have N levs x (G+2) columns of the format

c caStr80   = comment
c numlevs   = number of levels
c Psurf,TSurf,HSurf = surface pressure (mb) and temperature (K) and height (km)
c year lat lon     = year of obs (for co2 MR) and gravity adjust (not yet implemented)
c                     if year < 0, ignore
c                     if year > 0 and kPlanet = 3, then if adding on Gas2, adjust for ppmv
c numgases  = number of gases
c gasIDS    = which gases
c gasUNITS  = eg MR, g/g, RH, same as klayers
c   p1     T1   g1_1  g2_1 .... gG_1
c   p2     T2   g1_2  g2_2 .... gG_2
c   p3     T3   g1_3  g2_3 .... gG_3
c   p4     T4   g1_4  g2_4 .... gG_4
c
c   pN     TN   g1_N  g2_N .... gG_N
c
c where p(i) is in mb, T(i) is in Kelvin and N <= kProfLayer*2
c gi_j (j=1--N) are the gas amounts eg dimensionless 0 < MR < 1, RH etc
c
c        See eg KCARTA/IP_PROFILES/levelsRTP_to_levelstext.m
c	       KCARTA/IP_PROFILES/levelsprofile_text1.prf
c
c see /asl/packages/klayers/Doc/description.txt
c
c or can read TAPE5, modified TAPE6 from LBLRTM
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c for testing just do SUBROUTINE InputMR_profile(caPfName)    !! for testing
      SUBROUTINE InputMR_profile(caPfName,iNumLays,iNumGases,iaG,iLowestLev,
     $                           raTout,raAmountOut,raZout,
     $                           raPout,raaQout,raaPartPressout,
     $                           raPbndFinal,raTbndFinal,iZbndFinal)

      IMPLICIT NONE

      INTEGER iplev
      include '../INCLUDE/kcarta.param'
      include '../INCLUDE/KCARTA_database.param'
      include '../INCLUDE/airslevelheights.param'

c input
      CHARACTER*80 caPfName
c output
      REAL raTout(kProfLayer)                  !! in K
      REAL raAmountOut(kProfLayer)             !! in molecules/cm2
      REAL raZout(kProfLayer+1)                !! in meters, notice this is at LEVELS
      REAL raPout(kProfLayer)                  !! in N/m2 (even though input TAPE5,TAPE6, kRTP=-10, pressures are in mb)
      REAL raaQout(kProfLayer,kMaxGas)         !! in molecules/cm2
      REAL raaPartPressout(kProfLayer,kMaxGas) !! in N/m2
      INTEGER iNumLays,iNumGases,iaG(kMaxGas)
      INTEGER iLowestLev
      REAL raPbndFinal(kProfLayer+1)           !! final output pressure boundaries, default = 101 AIRS levels
      REAL raTbndFinal(kProfLayer+1)           !! final output temperatures at boundaries
      INTEGER iZbndFinal                       !! number of output plev bdries,
                                               !! should be iNumLays+1 if using LBLRTM bdies, else 101
      
c local
      INTEGER iG,iL,iJ,iK,iaGasUnits(kMaxGas)
      INTEGER iFound,iHighestLay,iNumLevs
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
      REAL raX(kMaxGas),rX,rP,rT,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
      REAL rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
      INTEGER iKCARTADirectionUPorDOWN,iMid,iFloor,iCnt,iOffSet
      REAL gamma,z,dz,slope,amount,junk
      REAL rFracBot
      REAL raPBnd(2*kProfLayer)  !! do we want user defined pressure level boundaries??
      INTEGER iZbnd              !! do we want user defined pressure level boundaries??

      REAL raPX(kProfLayer+1),raZX(kProfLayer+1),raTX(kProfLayer+1)
      REAL raaG_MRX(kProfLayer+1,kMaxGas),raLayDensityX(kProfLayer+1)
      REAL raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)
      
      CALL Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,
     $                    rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta) 

      DO iL = 1,kMaxLayer+1
        raPbndFinal(iL) = PLEV_KCARTADATABASE_AIRS(iL)
        raTbndFinal(iL) = 0.0
      END DO
      iZbndFinal = kMaxLayer+1

      iZbnd = +1   !! default to using the input levels as boundaries when doing levels --> layers integration
      iZbnd = -1   !! default to using AIRS 101 levels  as boundaries when doing levels --> layers integration
      
c >>> read in user supplied prof
c this will be at whatever gas units eg RH g/kg  VMR etc, but then CONVERTED to MR (gasunit 12)
      IF (kRTP .EQ. -10) THEN
        CALL ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,
     $                          raP,raT,iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)
	iZbnd = -1	
      ELSEIF (kRTP .EQ. -5) THEN
        CALL ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,
     $                          raP,raT,raAlt,iZbnd,raPBnd,
     $                          iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)
      ELSEIF (kRTP .EQ. -6) THEN
        !! edited TAPE6, which contains first few parts of TAPE5, and then just profile info from TAPE6 (mol/cm2)
        CALL ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLays,rPSurf,rTSurf,rHSurf,iNumGases,     
     $                          raPX,raTX,raLayDensityX,raZX,
     $                          raPressLevelsX,raTPressLevelsX,raAltitudesX,
     $                          iaG,iaGasUnits,raaG_MRX,rPMin,rPMax,rYear,rLat,rLon)
        iZbnd = -1
      ELSE
        write(kStdErr,*) 'huh?? reading text levels file or text LBLRTM TAPE5/TAPE6 MODIFIED',kRTP
        Call DoStop
      END IF

      IF (kRTP .NE. -6) THEN
        ! >>> interp onto the kCARTA database profiles
        ! >>> also sort the gasIDs into increasing order

        ! find pseudo temperature at this level
        ! recall pav = (p2-p1)/log(p2/p1)
        !        tav = (T2-T1)/log(p2/p1)
      
        !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
        !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
        !! have read LEVELS profile, so need to tack on necessary info above this profile, and integrate !!
	
        IF (iZbnd .LT. 0) THEN
	  ! >>>>>>>>>> default AIRS 101 levels

          ! >>>> tack on Standard Profile to User Profile; 
          ! >>>> if needed; also add on some more gases if needed, at bottom of profile to TOA
          ! old and new gases are ALWAYS going be in units = 12 (mixing ratio)
          ! have to be careful when adding in upper level info to old gases, as this could be a mix of units!!!!!!!!!!
          ! taper the tacked on profile so that within about 4 sublevels, you are just tacking on US Std (or default ref prof)
          CALL Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_MR,
     $                       rPMin,rPMax,rYear)

          ! now check everything already is in MR (gas unit 12)
          ! also adjust MR if this is planet Earth
          CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_MR)

          CALL InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_MR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,
     $                       PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_MRX,iLowestLev)

          rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/
     $         (PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
          write(kStdWarn,*) ' '
  	  write(kStdWarn,*) 'iX = iLowestLev ==>'
          write(kStdWarn,*)' PLEV_KCARTADATABASE_AIRS(iX),rPSurf,PLEV_KCARTADATABASE_AIRS(iX+1) = '
          write(kStdWarn,*) PLEV_KCARTADATABASE_AIRS(iLowestLev),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev+1)
          write(kStdWarn,*) 'rfracBot = ',rFracBot
  
          ! >>> do the integrals using integrals wrt p!!! closer to klayers
          CALL DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,
     $               PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,
     $               raPX,raTX,raaG_MRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
          iNumLays = kProflayer-iLowestLev+1  !! this is now LAYER number, max == kProfLayer when iLowestLev = 1

          DO iL = 1,101
            raPbndFinal(iL) = raPbndFinal(iL)/100.0    !! convert N/m2 to mb
          END DO

        ELSEIF (iZbnd .GE. 0) THEN
	  ! >>>>>>>>>> user boundaries
          DO iL = 1,iZbnd
            raPbndFinal(iL) = raPbnd(iL)/100.0    !! convert N/m2 to mb
          END DO
	  DO iL = iZbnd+1,kProfLayer+1
	    raPbndFinal(iL) = 0.005
	    raPbndFinal(iL) = raPbnd(iZbnd)/100.0   !! convert N/m2 to mb
	  END DO
          iZbndFinal = iZbnd

          ! now check everything already is in MR (gas unit 12)
          ! also adjust MR if this is planet Earth
          CALL adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_MR)

          CALL InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_MR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,
     $                       PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_MRX,iLowestLev,iZbnd,raPBnd)

          rFracBot = (rPSurf-raPBnd(iLowestLev+1))/
     $         (raPBnd(iLowestLev)-raPBnd(iLowestLev+1))
          write(kStdWarn,*) ' '
  	  write(kStdWarn,*) 'iX = iLowestLev ==>'
          write(kStdWarn,*)' raPBnd(iX),rPSurf,raPBnd(iX+1) = '
          write(kStdWarn,*) raPBnd(iLowestLev),rPSurf,raPbnd(iLowestLev+1)
          write(kStdWarn,*) 'rfracBot = ',rFracBot

          ! >>> do the integrals using integrals wrt p!!! closer to klayers
          CALL DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,
     $               raPBnd,raAlt,iZbnd,rfracBot,
     $               raPX,raTX,raaG_MRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)
          iNumLevs = iZbnd-iLowestLev+1
	  iNumLays = iNumLevs - 1       !! this is now LAYER number, max == iZbnd-1 when iLowestLev = 1

        END IF
           
      ELSEIF (kRTP .EQ. -6) THEN
        CALL DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon,
     $               PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,
     $               raPX,raTX,raLayDensityX,raaG_MRX,
     $               raPressLevelsX,raTPressLevelsX,raAltitudesX,
     $               raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)
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

c      do iL = 1,iNumlevs
c        print *,'hmm4',iL,raP(iL),raT(iL),raaG_MR(iL,1),raPX(iL),raTX(iL),raaG_MRX(iL,1)
c      end do

      RETURN
      END

c************************************************************************
      SUBROUTINE Init_n_layers(iKCARTADirectionUPorDOWN,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,
     $                    rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output
      REAL rPminKCarta,rPmaxKCarta,rPmin,rPmax,rHminKCarta,rHmaxKCarta
      INTEGER iKCARTADirectionUPorDOWN
      REAL PLEV_KCARTADATABASE_AIRS(kMaxLayer+1),DATABASELEVHEIGHTS(kMaxLayer+1)

c local
      INTEGER iL

      rPminKCarta = +1e10
      rPmaxKCarta = -1e10
      rHminKCarta = +1e10
      rHmaxKCarta = -1e10
      rPmin = +1e10
      rPmax = -1e10
      iKCARTADirectionUPorDOWN = +1   !! assume increading index == lower pressure ==> going UP
      DO iL = 1,kMaxLayer+1
        PLEV_KCARTADATABASE_AIRS(iL) = PLEV_KCARTADATABASE_AIRS(iL) * 100.0  !! change mb --> N/m2
        DATABASELEVHEIGHTS(iL)       = DATABASELEVHEIGHTS(iL) * 1000.0       !! change km --> m
        IF (PLEV_KCARTADATABASE_AIRS(iL) .GT. rPmaxKCarta) rPmaxKCarta = PLEV_KCARTADATABASE_AIRS(iL)
        IF (PLEV_KCARTADATABASE_AIRS(iL) .LT. rPminKCarta) rPminKCarta = PLEV_KCARTADATABASE_AIRS(iL)
        IF (DATABASELEVHEIGHTS(iL) .GT. rHmaxKCarta) rHmaxKCarta = DATABASELEVHEIGHTS(iL)
        IF (DATABASELEVHEIGHTS(iL) .LT. rHminKCarta) rHminKCarta = DATABASELEVHEIGHTS(iL)
      END DO
      IF ((PLEV_KCARTADATABASE_AIRS(1) .GT. PLEV_KCARTADATABASE_AIRS(2)) .AND. 
     c    (DATABASELEVHEIGHTS(1) .LT. DATABASELEVHEIGHTS(2))) THEN
        iKCARTADirectionUPorDOWN = +1   !! increasing index == lower pressure ==> going UP
      ELSEIF ((PLEV_KCARTADATABASE_AIRS(1) .LT. PLEV_KCARTADATABASE_AIRS(2)) .AND. 
     c     (DATABASELEVHEIGHTS(1) .GT. DATABASELEVHEIGHTS(2))) THEN
        iKCARTADirectionUPorDOWN = -1   !! decreasing index == higher pressure ==> going DOWN
        write (kStdErr,*) 'Bummer : being lazy, wanted this to be going UP'
        CALL DoStop
      ELSE
        write (kStdErr,*) 'Inconsistency : database pressures in one dir, heights in another'
       CALL DoStop
      END IF

      RETURN
      END

c************************************************************************
c this subroutine takes raaG_MR and adds on new gas profiles using raaR100MR
c only to rPmin < rP < rPmax
      SUBROUTINE AddInNewProfiles(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels,
     $                            iNumLevs,rPmin,rPmax,raP,raaG_MR)

      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/kcarta.param'

c input
      INTEGER iNumGases0,iaG0(kMaxGas),iNumGases,iaG(kMaxGas),iNumLevs,iRefLevels
      REAL rPmin,rPmax
      REAL raP(2*kProfLayer)
      REAL raaR100MR(kMaxLayer+10,kMaxGas),raR100Press(kMaxLayer+10)
c output
      REAL raaG_MR(2*kProfLayer,kMaxGas)

c local
      INTEGER iG,iL,i1,i2
      REAL raTempMR(kMaxLayer+10),raTempP(kMaxLayer+10),raTempMRX(2*kMaxLayer)

c first find pressures that span rPmin < rP < rPmax
      IF (raR100Press(iRefLevels) .GT. rPmin*100.0) THEN
        write(kStdErr,*) 'iRefLevels = ',iRefLevels
        write(kStdErr,*) 'oops min pressure in ref database = ',raR100Press(iRefLevels),' N/m2'
        write(kStdErr,*) 'should be SMALLER than min pressure in supplied user profile = ',rPmin*100.0,' N/m2'
        CALL DoStop
      END IF

      IF (raR100Press(1) .LT. rPmax*100.0) THEN
        write(kStdErr,*) 'iRefLevels = ',iRefLevels
        write(kStdErr,*) 'oops max pressure in ref database = ',raR100Press(1),' N/m2'
        write(kStdErr,*) 'should be GREATER than max pressure in supplied user profile = ',rPmax*100,' N/m2'
        CALL DoStop
      END IF

      write(kStdWarn,*) ' '

      i1 = iRefLevels
  10  CONTINUE
      IF ((raR100Press(i1) .LT. rPmin) .AND. (i1 .GT. 1)) THEN
        i1 = i1 - 1
        GOTO 10
      END IF
      i1 = min(i1 + 1,iReflevels)

      i2 = 1
  20  CONTINUE
      IF ((raR100Press(i2) .GT. rPmax) .AND. (i2 .LT. iRefLevels)) THEN
        i2 = i2 + 1
        GOTO 20
      END IF
      i2 = max(1,i2-1)
      write(kStdWarn,*) 'need to peruse ref profiles between levels ',i1,i2,' to span rPmin,rPmax'
      write(kStdWarn,*) 'Ref Database max press, user max press (in mb) = ',raR100Press(1)/100.0         ,rPmax
      write(kStdWarn,*) 'Ref Database min press, user min press (in mb) = ',raR100Press(iRefLevels)/100.0,rPmin

      DO iL = i2,i1
        raTempP(iL-i2+1) = raR100Press(iL)
      END DO
      DO iG = iNumGases0 + 1,iNumGases
        DO iL = i2,i1
          raTempMR(iL-i2+1) = raaR100MR(iL,iG)
        END DO
        Call r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
        DO iL = 1,iNumLevs
          raaG_MR(iL,iG) = raTempMRX(iL)
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this basically loops and reads in ref profile for gas iG, if necessary changes pressures to N/m2
c input gas list (from user) = iaG
c output gaslist = union(user list,preset list)
      SUBROUTINE ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,
     $                     iaG,iaGasUnits,iNumGases,raR100Press,raR100Temp,raaR100MR,
     $                     iRefLevels,laysORlevs)

      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/kcarta.param'

c these are the individual reference profiles, at kMaxLayer layers
c input/output
      INTEGER iaG(kMaxGas),iaGasUnits(kMaxGas),iNumGases    !! from user supplied list
      INTEGER laysORlevs                                    !! +1 if read 100 US Standard Ref layers, -1 if read 50 AFGL levels
      REAL PLEV_KCARTADATABASE_AIRS(kmaxLayer+1)
c output
      REAL raR100Temp(kMaxLayer+10),raaR100MR(kMaxLayer+10,kMaxGas),raR100Press(kMaxLayer+10)
      REAL raRx110Temp(kProfLayer+10),raRx110MR(kProfLayer+10),raRx110Press(kProfLayer+10)
      INTEGER iRefLevels    !!! iRefLevels ~ 100 + 10 extra levels really high up : see SUBR ReadRefProf_Levels

c local 
      INTEGER iL,iJ,iMid,ifloor,iErr,iG,iaPreset(kMaxGas),iFound
      REAL rX,raR100Amt(kMaxLayer+10),raR100PartPress(kMaxLayer+10)
      CHARACTER*80 caRefFName
      CHARACTER*3 ca3
      INTEGER iEarth,iMars,iaEarth(8),iaMars(2),iPreset,iNumLevsX,iIndex

      DATA (iaEarth(iG),iG=1,8) /01,02,03,04,05,06,09,12/      !! 8 important Earth atm molecules
      DATA (iaMars(iG),iG=1,2)  /02,22/                        !! 2 important Mars  atm molecules
      iEarth = 8
      iMars = 2
       
      DO iG = 1,kMaxGas
        iaPreset(iG) = -1
      END DO
      IF (kPlanet == 03) THEN
        iPreset = iEarth
        !! earth
        DO iG = 1,iEarth
          iaPreset(iG) = iaEarth(iG)
        END DO
      ELSEIF (kPlanet == 04) THEN
        iPreset = iMars
        !! mars
        DO iG = 1,iMars
          iaPreset(iG) = iaMars(iG)
        END DO
      ELSE
        write (kStdErr,*) 'oops right now can only handle Earth or Mars. Sorry, Mork'
        Call DoStop
      END IF
        
      !! now do a union of iaG and iaPreset
      DO iG = 1,iPreset
        iFound = -1
        DO iL = 1,iNumGases
          IF (iaG(iL) .EQ. iaPreset(iG)) iFound = +1
        END DO
        IF (iFound .LT. 0) THEN
          write(kStdWarn,*)'gasID ',iaPreset(iG),' not found in user list, adding .... '
          iNumGases = iNumGases + 1
          iaG(iNumGases) = iaPreset(iG)
          iaGasUnits(iNumGases) = 12   !! VMR
        END IF
      END DO

      DO iG = 1,iNumGases
        IF ((kPlanet .NE. 3) .OR. (laysORlevs .EQ. +1)) THEN
          !! this is old : use only (US) Standard Profile
          ! read kCARTA kProfLayer reference profile
          CALL FindReferenceName(caRefFName,iaG(iG),-1)
          CALL ReadRefProf(caRefFName,kMaxLayer,raR100Amt,
     $         raR100Temp,raR100Press,raR100PartPress,iErr)
          iRefLevels = kMaxLayer
          ! change pressures from atm to N/m2
          DO iL = 1,kMaxLayer
            raR100Press(iL)     = raR100Press(iL) * kAtm2mb*100.0
            raR100PartPress(iL) = raR100PartPress(iL) * kAtm2mb*100.0
            raaR100MR(iL,iG)    = raR100PartPress(iL)/raR100Press(iL)
          END DO

        ELSEIF ((kPlanet .EQ. 3) .AND. (laysORlevs .EQ. -1)) THEN
	  !ca3 = kcaLevsRefProf(1:3)
	  ca3 = 'DNE'
	  iIndex = index(kcaLevsRefProf,'DNE')
          IF (iIndex .GT. 0) THEN
            write(kStdErr,*) 'oops : do not have a set of Levels Reference Profiles'
            CALL DoStop
          END IF
          !! this is new : can use one of the (AFGL) models
          CALL ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iaG(iG),iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
          iRefLevels = iNumLevsx
          DO iL = 1,iNumLevsx
            raR100Press(iL)     = raRx110Press(iL)       !! already in N/m2
            raR100Temp(iL)      = raRx110Temp(iL)
            raaR100MR(iL,iG)    = raRx110MR(iL) /1.0e6   !! change from ppmv to MR
            raR100PartPress(iL) = raR100Press(iL) * raaR100MR(iL,iG)
          END DO
          DO iL = iNumLevsx+1,kMaxLayer
            raR100Press(iL)     = 0.0
            raaR100MR(iL,iG)    = 0.0
            raR100PartPress(iL) = 0.0
          END DO
        END IF

        ! make sure pressures are decreasing with index ie layers going higher and higher
        IF (raR100Press(1) .LT. raR100Press(2)) THEN
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
      END

c************************************************************************
c this subroutine changes the input levels profiles to ppmv to VMR 
c and if planet Earth, adjusts mix ratios
      SUBROUTINE adjustVMR_Earth(iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_MR)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL raP(2*kProfLayer),raT(2*kProfLayer)
c input/output
      REAL raaG_MR(2*kProfLayer,kMaxGas)

c local
      INTEGER iL,iG,iDry2Wet,iWat,iFound

      iFound = -1
      iWat = 1
 10   CONTINUE
      IF (iaG(iWat) .EQ. 1) THEN
        iFound = +1
      ELSE
        iWat = iWat + 1
        IF (iWat .LE. iNumGases) THEN
          GOTO 10
        END IF
      END IF

      IF ((kPlanet .EQ. 3) .AND. (iFound .EQ. 1)) THEN
        write(kStdWarn,*) 'Planet = 3, EARTH .. gas ',iWat,' corresponds to water'
        IF (iWat .NE. 1) THEN
          write(kStdErr,*) 'Planet = 3, EARTH .. found water but it should have been g1, not g',iWat
          CALL DoStop
        END IF
      ELSEIF ((kPlanet .EQ. 3) .AND. (iFound .NE. 1)) THEN
        write(kStdErr,*) 'Planet = 3, EARTH .. no input gas corresponds to water'
         CallDoStop
      END IF

      !! change to ppmv
      DO iG = 1,iNumGases
c        write(kStdWarn,*)iG,iaG(iG),iaGasUnits(iG)
        IF (iaGasUnits(iG) .NE. 10) THEN
          CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR)
        END IF
      END DO

      iDry2Wet = -1 !! do not adjust dry air to wet air MR
      iDry2Wet = +1 !!        adjust dry air to wet air MR
      !! now change from ppmv to VMR (which is gas units 12)

      ! IF (kRTP == -20) iDry2Wet = -1

      IF (iDry2Wet .LT. 0) THEN
        DO iG = 1,iNumGases
          DO iL = 1,iNumLevs
            raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
          END DO
        END DO

      !! need to make sure gas units code == 10.11.12.20.21 before trying this; see toppmv.f in klayers
      ELSEIF ((iDry2Wet .GT. 0) .AND. (kPlanet .EQ. 3) .AND. (iaG(1) .EQ. 1)) THEN

        !!! first convert WV mixing ratios
        !!! however, klayars assumes the user has provided WET mixing ratios, so NO NEED to do this
c        DO iG = 1,1
c          DO iL = 1,iNumLevs
c            print *,iL,raaG_MR(iL,iG),1e6*raaG_MR(iL,iG) / (1.0e6+raaG_MR(iL,iG))
c            raaG_MR(iL,iG) =  1.0e6 * raaG_MR(iL,iG) / (1.0e6 + raaG_MR(iL,iG))
c          END DO
c        END DO

        !!! then use WV mixing ratios to fix other gas mixing ratios if they are 10.11.12.20.21
        DO iG = 2,iNumGases
          IF ((iaGasUnits(iG) .GE. 10) .AND. (iaGasUnits(iG) .LE. 21)) THEN
            DO iL = 1,iNumLevs
              raaG_MR(iL,iG) =  raaG_MR(iL,iG) * (1.0e6 - raaG_MR(iL,1))/1e6
            END DO
          END IF
        END DO

        !! finally convert ppmv --> MR
        DO iG = 1,iNumGases
          DO iL = 1,iNumLevs
            raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
          END DO
        END DO

      ELSE
        write(kStdErr,*) 'Need iDry2Wet == +/- 1'
        Call DoStop
      END IF

      !! check all MR larger than 0
      DO iG = 1,iNumGases
        DO iL = 1,kMaxLayer
          raaG_MR(iL,iG) = max(raaG_MR(iL,iG),0.0)
        END DO
      END DO
 
      RETURN
      END

c************************************************************************
c this does the actual change for raaG_MR(iL,iG)
C klayers/Doc/gas_units_code.txt
C klayers/Doc/toppmv.f
C              ------- 10-19 = volume mixing ratio --------------------- 
C         10   parts per million volume mixing ratio (ppmv)
C              Number of gas X molecules per 1E+6 "air" molecules
C
C         11   parts per billion volume mixing ratio (ppbv)
C              Number of gas X molecules per 1E+9 "air" molecules
C
C         12   volume mixing ratio (unitless fraction)
C              Number of gas X molecules per "air" molecule
C
C              ------- 20-29 = mass mixing ratio -----------------------
C         20   mass mixing ratio in (g/kg)
C              Grams of gas X per kilogram of "air"
C
C         21   mass mixing ratio in (g/g)
C              Grams of gas X per gram of "air"
C
C              ------- 30-39 = partial pressure ------------------------
C
C         30   partial pressure in millibars (mb) Note: mb=hPa
C              Pressure of gas X as it exists in the atmosphere.
C              To clarify, this means at the corresponding profile
C              temperature and pressure level total pressure.
C
C         31   partial pressure in atmospheres (atm)
C              Pressure of gas X as it exists in the atmosphere
C
C              ------- 40-49 = water vapor humidity units --------------
C
C         40   relative humidity in (percent)
C              100 times actual vapor pressure divided by saturation
C              vapor pressure
C
C         41   relative humidity (unitless fraction)
C              Actual vapor pressure divided by saturation vapor
C              pressure
C
C         42   dew point temperature (Kelvin)
C              Temperaure at which vapor will start to condense out
C
C         43   dew point temperature (Celcius)
C              Temperaure at which vapor will start to condense out
C
C              Possible additional units might be
C                 inches or centimenters of water vapor
C                 grams of water vapor

c
c pIN is in mb
c 
      SUBROUTINE changeLVLS_2_ppmv(iGasID,iGasUnits,iNumLevs,iG,PIN,TIN,raaG_MR)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iGasID,iGasUnits,iNumLevs,iG
      REAL PIN(2*kProfLayer),TIN(2*kProfLayer)
c input/output
      REAL raaG_MR(2*kProfLayer,kMaxGas)

c local
      INTEGER iL,iCode,NLEV
      REAL MRIN(2*kProfLayer,kMaxGas),RJUNK,WEXSVP,MDAIR
      INTEGER iaLocalGasID(kMaxGas)
      REAL MASSF(kMaxGas)

c   1 H2O, water
c   2 CO2, carbon dioxide
c   3 O3, ozone
c   4 N2O, nitrous oxide
c   5 CO, carbon monoxide
c   6 CH4, methane
c   7 O2, (diatomic) oxygen
c   8 NO, nitric oxide
c   9 SO2, sulfur dioxide
c  10 NO2
c  11 NH3, ammonia
c  12 HNO3, nitric acid

c see cbgids.f in klayers
      DATA (iaLocalGasID(iL),iL=1,12) /01,02,03,04,05,06,07,08,09,10,11,12/
      DATA (MASSF(iL),iL=1,12) /18.015,44.010,47.9982,44.013,28.011,16.043,31.999,30.006,64.063,46.006,17.031,63.013/

      IF (iGasID .GT. 12) THEN
        MASSF(iGasID) = kAtmMolarMass
      ENDIF

      MDAIR = kAtmMolarMass

      DO iL = 1,iNumLevs
        MRIN(iL,iG) = raaG_MR(iL,iG) 
      END DO

c convert to ppmv (gas unit 10)
      iCode = iGasUnits
      NLEV  = iNumLevs

C            parts per billion volume mixing ratio
             IF (ICODE .EQ. 11) THEN
C               PPMV = PPBV*1E-3
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*1E-3
                ENDDO
C
C            volume mixing ratio
             ELSEIF (ICODE .EQ. 12) THEN
C               PPMV = VMR*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*1E+6
                ENDDO
C
C            mass mixing ratio in g/kg
             ELSEIF (ICODE .EQ. 20) THEN
C               PPMV = MRgkg*((MDAIR*1E-3)/MASSF)*1E+6
                RJUNK=1E+3*MDAIR/MASSF(IGasID)
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
                ENDDO
C
C            mass mixing ratio in g/g
             ELSEIF (ICODE .EQ. 21) THEN
C               PPMV = MRgg*(MDAIR/MASSF)*1E+6
                RJUNK=1E+6*MDAIR/MASSF(IGasID)
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*RJUNK
                ENDDO
C
C            partial pressure in mb
             ELSEIF (ICODE .EQ. 30) THEN
C            PPMV = (PPmb/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=(MRIN(IL,IG)/PIN(IL))*1E+6
                ENDDO
C
C            partial pressure in atm
             ELSEIF (ICODE .EQ. 31) THEN
C               PPMV = (PPatm*1013.25/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=(MRIN(IL,IG)*1013.25/PIN(IL))*1E+6
                ENDDO
C
C            relative humidy in percent
c            note we need to change PN fron N/m2 to mb, so divide by 100
             ELSEIF (ICODE .EQ. 40 .AND. iGasID .EQ. 1) THEN
C            PPMV = (RH%/100)*(SVP/PIN)*1E+6
                DO IL=1,NLEV
c                print *,iL,MRIN(IL,IG),TIN(iL),PIN(IL)/100.0,WEXSVP( TIN(IL))
                   MRIN(IL,IG)=MRIN(IL,IG)*
     $                (WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+4
                ENDDO
c                call dostop
C
C            relative humidity (fraction)
             ELSEIF (ICODE .EQ. 41 .AND. iGasID .EQ. 1) THEN
C               PPMV = RH*(SVP/PIN)*1E+6
                DO IL=1,NLEV
                   MRIN(IL,IG)=MRIN(IL,IG)*
     $                (WEXSVP( TIN(IL) )/(PIN(IL)/100.0))*1E+6
                ENDDO
C
             ENDIF

c save
      DO iL = 1,iNumLevs
        raaG_MR(iL,iG) = MRIN(iL,iG)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine tacks on Standard Profile above what the user has given, to fill in amounts till 0.005 mb
c also if kCARTA has earlier determined there should be more gases than in user input profile (eg suppose
c user only gave WV/O3 but kCARTA also wants other gases), then this is the place additional profiles added in
      SUBROUTINE Tack_on_profile(PLEV_KCARTADATABASE_AIRS,iaG,iaGasUnits,iNumLevs,iNumGases,raP,raT,raaG_MR,
     $                           rPMin,rPMax,rYear)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input var
      INTEGER iNumLevs
      REAL rYear
      REAL PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
c input/output var
      REAL rPmin,rPmax
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas)
      INTEGER iaG(kMaxGas),iaGasUnits(kMaxGas),iNumGases

c local var
      !! individual reference profiles, at kMaxLayer layers, in terms of MR = PPress/Press
      INTEGER iaG0(kMaxGas),iNumGases0,i2,iFound,iaGasUnits0(kMaxGas),iIndex
      REAL raaR100MR(kMaxLayer+10,kMaxGas),raR100Temp(kMaxLayer+10),raR100Press(kMaxLayer+10)
      INTEGER iL,iG,iAbove,iMerge,iRefLevels,laysORlevs
      REAL rToffset,rJunk,raOffset(kMaxGas),raJunk(kMaxLayer+10)
      REAL rCO2ppmv,rAdjust,rCH4ppbv
      CHARACTER*3 ca3

c first read in reference profile, and if needed info for other gases
c     read in additional gases eg from ERA/ECM we get gas 1,3 but we would like gas 1,2,3,4,5,6,8,12
c     read in ref profiles (eg for Gas2 since Mars, Earth, Venus all have CO2 in the atm; get the Standard MR)
      iNumGases0 = iNumGases
      DO iL = 1,iNumGases0
        iaG0(iL) = iaG(iL)
        iaGasUnits0(iL) = iaGasUnits(iL)
      END DO
      laysORlevs = +1   !! orig, read in reference P,PP and get mix ratio using PP/P for each gas
      laysORlevs = -1   !! new,  read in one of 6 AFGL P/T/ppmv level profiles
      !ca3 = kcaLevsRefProf(1:3)
      ca3 = 'DNE'
      iIndex = index(kcaLevsRefProf,'DNE')
      IF ((laysORlevs .EQ. -1) .AND. (iIndex .GT. 0)) THEN
        write(kStdWarn,*) 'oops : do not have a set of Levels Reference Profiles; code is setting laysORlevs = +1'
        laysORlevs = +1
      END IF
      
      Call ReadRefProf_Units_laysORlevs(PLEV_KCARTADATABASE_AIRS,
     $                                  iaG,iaGasUnits,iNumGases,raR100Press,raR100Temp,raaR100MR,iRefLevels,laysORlevs)
      IF (iNumGases0 .NE. iNumGases) THEN
        write(kStdWarn,*)'Orig num of gases, new num of gases = ',iNumGases0,iNumGases
        !! new gases are all at MR (dimensionless fraction) = gasunit 12
        !! need to massage in these new profiles into raaG_MR, between rPmax and rPmin
        CALL AddInNewProfiles(iNumGases0,iaG0,iNumGases,iaG,raR100Press,raaR100MR,iRefLevels,
     $                        iNumLevs,rPmin,rPmax,raP,raaG_MR)
      END IF

c then find KCARTA levels above which there is NO user supplied info
      write(kStdWarn,*)' user supplied profile : Pmax(mb),Pmin(mb),numlevs = ',rPMax,rPMin,iNumLevs
      write(kStdWarn,*)' kCARTA Pav Database   : Pmax(mb),Pmin(mb),numlevs = ',
     $     raR100Press(1)/100.0,raR100Press(iRefLevels)/100.0,iRefLevels
      
      iAbove = kMaxLayer
      iAbove = iRefLevels
  10  CONTINUE
      IF ((raR100Press(iAbove)/100.0 .LT. rPMin) .AND. (iAbove .GT. 1)) THEN
        iAbove = iAbove - 1
        GOTO 10
      END IF
      iAbove = iAbove + 1
      IF (iAbove .GT. iRefLevels) THEN
        write(kStdWarn,*) 'Will be tacking on reference profile info from level ',iAbove,' to ',iRefLevels
        write(kSTdWarn,*) 'Hmm ... no need to do this ... looks like user supplied levels profile past 0.005 mb'
        GOTO 123
      END IF
      
      write(kStdWarn,*) 'Will be tacking on reference profile info from level ',iAbove,' to ',iRefLevels
      write(kStdWarn,*) 'This should extend user supplied level info from level ',iNumLevs,' to ',
     $ iNumLevs + (iRefLevels-iAbove+1) 
      
c now do linear interpolation from iAbove pressure, down to min(raP)
      Call r_sort_loglinear1(raR100Press,raR100Temp,iRefLevels,rPMin*100.0,rJunk,1)
      rToffset = raT(iNumLevs) - rJunk 
      write(kStdWarn,*)'tacking on info from kCARTA Pav Dtabase, layers ',iAbove,' to ',iRefLevels
      write(kStdWarn,*)'ToffSet = ',rToffset

c also need to figure out the gas multiplier offset
      DO iG = 1, iNumGases
        DO iL = 1,iRefLevels
          raJunk(iL) = raaR100MR(iL,iG)
        END DO
        Call r_sort_loglinear1(raR100Press,raJunk,iRefLevels,rPMin*100.0,rJunk,1)
        !! need to do a CHANGE OF UNITS to ppmv!!!!!
        raoffset(iG) = raaG_MR(iNumLevs,iG)/rJunk 
        write(kStdWarn,*) 'Multiplier for gasID = ',iaG(iG),' units = ',iaGasUnits(iG),' = ',raoffset(iG)
      END DO

      iMerge = +1    !! in other words, 3 points away from iMerge, we bump up/down the multiplier to 1.000
      DO iL = iAbove,iRefLevels
        iNumLevs = iNumLevs + 1
        raP(iNumLevs) = raR100Press(iL)

        !raT(iNumLevs) = raR100Temp(iL) + rToffset
        IF ((iL-iAbove) .GT. iMerge) THEN
          rAdjust = 0.0  !! make the temp adjustment effectively 0
        ELSE
          rJunk = (0 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 0, inside 4 levels; this is slope
          rAdjust = rJunk * (iL-iAbove) + rToffset
        END IF
        raT(iNumLevs) = raR100Temp(iL) + rAdjust

        DO iG = 1,iNumGases
          !raaG_MR(iNumLevs,iG) = raaR100MR(iL,iG) * raOffset(iG)
          IF ((iL-iAbove) .GT. iMerge) THEN
            rAdjust = 1  !! make the MR adjustment effectively 1
          ELSE
            rJunk = (1 - raOffset(iG))/(iMerge) !! adjust rOffSet from its given value to 1, inside 4 levels; this is slope
            rAdjust = rJunk * (iL-iAbove) + raOffset(iG)
c            print *,iG,iL,iL-iAbove,raOffset(iG),rAdjust
          END IF
          raaG_MR(iNumLevs,iG) = raaR100MR(iL,iG) * rAdjust
        END DO
      END DO
      write(kStdWarn,*) 'Have extended user supplied info to level',iNumLevs

 123  CONTINUE
c >>>>>>>>>>>>>>>>>>>>>>>>> CO2 adjust; growth rate = 2 ppmv/yr
      rCO2ppmv = 370.0
      IF ((rYear .GT. 1970) .AND. (kPlanet .EQ. 03)) THEN
        !! should be dt = (yy-2002) + mm/12;  rCO2ppmv = 370 + 2*dt;
        rCO2ppmv = 370.0 + (rYear-2002)*2.0

        !! now see if CO2 was originally there
        iFound = -1
        i2 = 1
 20     CONTINUE        
        IF (iaG0(i2) .EQ. 2) THEN
          iFound = +1
        ELSEIF (i2 .LT. iNumGases0) THEN
          i2 = i2 + 1
          GOTO 20
        END IF

        IF (iFound .GT. 0) THEN
          !! user had CO2 in supplied prof, no need to adjust MR
          write(kStdWarn,*) 'User has CO2 in supplied profile, so not adjusting bottom levels'
        ELSE
          !! user did not have CO2 in there, yes need to adjust MR
          i2 = 1
          iFound = -1
 25       CONTINUE        
          IF (iaG(i2) .EQ. 2) THEN
            iFound = +1
          ELSEIF (i2 .LT. iNumGases) THEN
            i2 = i2 + 1
            GOTO 25
          END IF
          IF (iFound .LT. 0) THEN
            write(kStdErr,*) 'hmm, expected to find CO2!!! ooops'
            CALL DoStop
          END IF
          rJunk = 0.0
          DO iL = 3,7
            rJunk = rJunk + raaG_MR(iL,i2)
          END DO
          rJunk = rJunk/5.0
          write(kStdWarn,*) '  reference CO2 profile that was read in has CO2 MR = ',rJunk*1.0e6,' ppmv'
          write(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCO2ppmv,' ppmv'
          DO iL = 1,iNumLevs
            raaG_MR(iL,i2) = raaG_MR(iL,i2) * rCO2ppmv/(rJunk*1.0e6)
          END DO
        END IF
      END IF

c >>>>>>>>>>>>>>>>>>>>>>>>> CH4 adjust; growth rate = 0.01 ppmb/yr
      rCH4ppbv = 1.7
      IF ((rYear .GT. 1970) .AND. (kPlanet .EQ. 03)) THEN
        rCH4ppbv = 1.7 + (rYear-2002)*0.01
        rCH4ppbv = 1.8   !! hard coded into klayers right now

        !! now see if CH4 was originally there
        iFound = -1
        i2 = 1
 30     CONTINUE        
        IF (iaG0(i2) .EQ. 6) THEN
          iFound = +1
        ELSEIF (i2 .LT. iNumGases0) THEN
          i2 = i2 + 1
          GOTO 30
        END IF

        IF (iFound .GT. 0) THEN
          !! user had CH4 in supplied prof, no need to adjust MR
          write(kStdWarn,*) 'User has CH4 in supplied profile, so not adjusting bottom levels'
        ELSE
          !! user did not have CH4 in there, yes need to adjust MR
          i2 = 1
          iFound = -1
 35       CONTINUE        
          IF (iaG(i2) .EQ. 6) THEN
            iFound = +1
          ELSEIF (i2 .LT. iNumGases) THEN
            i2 = i2 + 1
            GOTO 35
          END IF
          IF (iFound .LT. 0) THEN
            write(kStdErr,*) 'hmm, expected to find CH4!!! ooops'
            CALL DoStop
          END IF
          rJunk = 0.0
          DO iL = 3,7
            rJunk = rJunk + raaG_MR(iL,i2)
          END DO
          rJunk = rJunk/5.0
          write(kStdWarn,*) '  reference CH4 profile that was read in has CH4 MR = ',rJunk*1.0e6,' ppbv'
          write(kStdWarn,*) '  this profile will be adjusted to "new" profile with ',rCH4ppbv,' ppbv'
          DO iL = 1,iNumLevs
            raaG_MR(iL,i2) = raaG_MR(iL,i2) * rCH4ppbv/(rJunk*1.0e6)
          END DO
        END IF
      END IF

c>>>>>>>>>>>>>>>>>>>>>>>>>

      !! check all MR larger than 0
      DO iG = 1, iNumGases
        DO iL = 1,kMaxLayer
          raaG_MR(iL,iG) = max(raaG_MR(iL,iG),0.0)
        END DO
      END DO

      RETURN
      END

c************************************************************************
c this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
      SUBROUTINE InterpUser2kCARTA(iNumLevs,iNumGases,iaG,raP,raT,raaG_MR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,
     $                             PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_MRX,iLowestLev)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input var
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rPmaxKCarta
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas)
      REAL PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
c output var
      INTEGER iLowestLev
      REAL raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas)
      REAL rInJunk,rOutJunk

c local var
      INTEGER iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
      REAL raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_MRX2(kProfLayer+1,kMaxGas)

      write(kStdWarn,*) ' ----->>>> using default AIRS 101 levels for layers integration'
      
c see which kCARTA pressure level is equal to or greater than surface pressure
      write(kStdWarn,*) '--------------------------------------------------------------------------'
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' >>> prep : klayers within kcarta : interp usergrid onto 101 levels grid >>>'
      
      iLowestLev = kMaxLayer
 10   CONTINUE
      IF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) .LT. rPSurf) .AND. (iLowestLev .GT. 1)) THEN
        iLowestLev = iLowestLev - 1
        GOTO 10
      ELSEIF ((PLEV_KCARTADATABASE_AIRS(iLowestLev) .LT. rPSurf) .AND. (iLowestLev .EQ. 1)) THEN
        write (kStdErr,*) 'PSurf = ',rPSurf,' and Max Kcarta Database Plev = ',rPmaxKCarta 
        Call DoStop
      ELSEIF (PLEV_KCARTADATABASE_AIRS(iLowestLev) .GE. rPSurf) THEN
        write(kStdWarn,*)'lowest level = ',iLowestLev
      END IF

c fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
      DO iL = 1,kMaxLayer-iLowestLev+1
        iNumUse = iL
        raPX(iL) = PLEV_KCARTADATABASE_AIRS(iLowestLev+iL)/100.0  ! input press in mb, so change plev_kcarta to mb
      END DO
      
c now interpolate T and MR onto this grid
      Call r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raTemp(iL) = raaG_MR(iL,iG)
        END DO
        Call r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
        DO iL = 1,iNumUse
          raaG_MRX(iL,iG) = raTempX(iL)
        END DO
      END DO
      
c      DO iL = 1,kMaxLayer-iLowestLev+1
c        write(*,1234) iL,raP(iL),raT(iL),raaG_MR(iL,1),raPX(iL),raTX(iL),raaG_MRX(iL,1)
c      END DO
c      call dostop
 1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
 
c need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
      iHigh = iNumUse
 20   CONTINUE
      IF (raP(iNumLevs) .GT. raPX(iHigh)) THEN
        iHigh = iHigh - 1
        GOTO 20
      END IF
      iHigh = iHigh + 1
      IF (iHigh .LE. iNumUse) THEN
        write(kSTdWarn,*) ' '
        write(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
        write(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
        write(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' to ',iNumUse
        DO iL = iHigh,iNumUse
          rInJunk = raPX(iL)
          Call r_sort_loglinear1(raP,raT,iNumLevs,rInJunk,rOutJunk,1)
          raTX(iL) = rOutJunk

          DO iG = 1,iNumGases
            DO iX = 1,iNumLevs
              raTemp(iX) = raaG_MR(iX,iG)
            END DO
          Call r_sort_loglinear1(raP,raTemp,iNumLevs,rInJunk,rOutJunk,1)
          raaG_MRX(iL,iG) = rOutJunk
          END DO
        END DO
      END IF

c now sort the gasIDs into increasing order, and do the same for the mixing ratios
      DO iL = 1,iNumGases
        raArr(iL) = iaG(iL) * 1.0
      END DO
      CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)

c assign temp arrays
      DO iG = 1,iNumGases
        iaG2(iG) = iaG(iG)
      END DO
      DO iL = 1,iNumUse
        DO iG = 1,iNumGases
          raaG_MRX2(iL,iG) = raaG_MRX(iL,iG)
        END DO
      END DO

c sort according to iaIndex
      DO iG = 1,iNumGases
        iaG(iG) = iaG2(iaIndex(iG))
      END DO
      DO iL = 1,iNumUse
        DO iG = 1,iNumGases
          raaG_MRX(iL,iG) = raaG_MRX2(iL,iaIndex(iG))
        END DO
      END DO

      write(kStdWarn,*) '--------------------------------------------------------------------------'
      
      RETURN
      END

c************************************************************************
c this subroutine takes in user/tacked on profile, and puts it onto the kCARTA Database levels
      SUBROUTINE InterpUser2UserBnd(iNumLevs,iNumGases,iaG,raP,raT,raaG_MR,rPMin,rPMax,rPSurf,rTSurf,rPmaxKCarta,
     $                             PLEV_KCARTADATABASE_AIRS,raPX,raTX,raaG_MRX,iLowestLev,iZbnd,raPbnd)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input var
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rPmaxKCarta
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas)
      REAL PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      REAL raPbnd(2*kProfLayer)  !! user defined boundaries
      INTEGER iZbnd              !! user defined boundaries
c output var
      INTEGER iLowestLev
      REAL raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas),rInJunk,rOutJunk

c local var
      INTEGER iL,iG,iNumUse,iHigh,iX,iaIndex(kMaxGas),iaG2(kMaxGas)
      REAL raTemp(2*kProfLayer),raTempX(kProfLayer),raArr(kMaxGas),raaG_MRX2(kProfLayer+1,kMaxGas)

      write(kStdWarn,*) ' ----->>>> using user defined ',iZbnd,' pressure levels for layers integration'
      
c see which user defined pressure level is equal to or greater than surface pressure
      write(kStdWarn,*) '--------------------------------------------------------------------------'
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' >>> prep : user defined levels : interp usergrid onto N levels grid >>>'

      rPmaxKCarta = -1.0
      DO iL = 1,iZbnd
        print *,'in InterpUser2UserBnd',iL,raPbnd(iL)
        IF (raPbnd(iL) .GT. rPmaxKCarta) rPmaxKCarta = raPbnd(iL)
      END DO
      write(kStdWarn,*) 'User defined levels for integration : maxP = ',rPmaxKCarta,' Psurf = ',rPSurf,' mb'

      !! assume raPbnd(1) > raPbnd(2) > raPbnd(3) > ... > raPbnd(iZbnd) and look
      !! for raPbnd greater than raPSurf
      iLowestLev = iZbnd
 10   CONTINUE
      IF ((raPbnd(iLowestLev) .LT. rPSurf) .AND. (iLowestLev .GT. 1)) THEN
        iLowestLev = iLowestLev - 1
        GOTO 10
      ELSEIF ((raPbnd(iLowestLev) .LT. rPSurf) .AND. (iLowestLev .EQ. 1)) THEN
        write (kStdErr,*) 'Error!!!  PSurf = ',rPSurf,' and Max UserLev Plev = ',rPmaxKCarta 
        Call DoStop
      ELSEIF (raPbnd(iLowestLev) .GE. rPSurf) THEN
        write(kStdWarn,*)'lowest level = ',iLowestLev
      END IF

c fill in raPressureJunk, which is the pressure grid onto which we need to interpolate T,MR onto
c remember this comes from LBLRTM input profile, so lowest layer will probably CLOSELY or EXACTLY match pSurf
      DO iL = 1,iZbnd
        raPX(iL) = raPbnd(iL)/100.0  ! input raPbnd press in N/m2, so change plev_kcarta to mb
      END DO
      iNumUse = iZbnd
      
c now interpolate T and MR onto this grid
      Call r_sort_loglinear(raP,raT,iNumLevs,raPX,raTX,iNumUse)
      DO iG = 1,iNumGases
        DO iL = 1,iNumLevs
          raTemp(iL) = raaG_MR(iL,iG)
        END DO
        Call r_sort_loglinear(raP,raTemp,iNumLevs,raPX,raTempX,iNumUse)
        DO iL = 1,iNumUse
          raaG_MRX(iL,iG) = raTempX(iL)
        END DO
      END DO
      
c      DO iL = 1,iZbnd-iLowestLev+1
c        write(*,1234) iL,raP(iL),raT(iL),raaG_MR(iL,1),raPX(iL),raTX(iL),raaG_MRX(iL,1)
c      END DO
c      call dostop
 1234 FORMAT(I3,2(' ',F10.4),' ',E10.4,2(' ',F10.4),' ',E10.4)
 
c need to be careful with last point; linear maybe better than spline if raP(iNumLevs) > raPX(iNumUse)
      iHigh = iNumUse
 20   CONTINUE
      IF (raP(iNumLevs) .GT. raPX(iHigh)) THEN
        print *,'wah iHigh',iNumLevs,iHigh,raP(iNumLevs),raPX(iHigh)
        iHigh = iHigh - 1
        GOTO 20
      END IF

      iHigh = iHigh + 1
      IF (iHigh .LE. iNumUse) THEN
        write(kSTdWarn,*) ' '
        write(kStdWarn,*) 'Tacked on profile ends at LOWER pressure than plev_kcartadatabase_airs'
        write(kSTdWarn,*) 'ie,  raP(iNumLevs) > raPX(iNumUse) : ',raP(iNumLevs),raPX(iNumUse)
        write(kSTdWarn,*) 'Replace profile values (done with spline) with linear interp, between',iHigh,' to ',iNumUse
        DO iL = iHigh,iNumUse
          rInJunk = raPX(iL)
          Call r_sort_loglinear1(raP,raT,iNumLevs,rInJunk,rOutJunk,1)
          raTX(iL) = rOutJunk

          DO iG = 1,iNumGases
            DO iX = 1,iNumLevs
              raTemp(iX) = raaG_MR(iX,iG)
            END DO
          Call r_sort_loglinear1(raP,raTemp,iNumLevs,rInJunk,rOutJunk,1)
          raaG_MRX(iL,iG) = rOutJunk
          END DO
        END DO
      END IF

c now sort the gasIDs into increasing order, and do the same for the mixing ratios
      DO iL = 1,iNumGases
        raArr(iL) = iaG(iL) * 1.0
      END DO
      CALL NumericalRecipesIndexer(iaIndex,raArr,iNumGases)

c assign temp arrays
      DO iG = 1,iNumGases
        iaG2(iG) = iaG(iG)
      END DO
      DO iL = 1,iNumUse
        DO iG = 1,iNumGases
          raaG_MRX2(iL,iG) = raaG_MRX(iL,iG)
        END DO
      END DO

c sort according to iaIndex
      DO iG = 1,iNumGases
        iaG(iG) = iaG2(iaIndex(iG))
      END DO
      DO iL = 1,iNumUse
        DO iG = 1,iNumGases
          raaG_MRX(iL,iG) = raaG_MRX2(iL,iaIndex(iG))
        END DO
      END DO

      write(kStdWarn,*) '--------------------------------------------------------------------------'
      
      RETURN
      END

c************************************************************************
c function to compute vapor pressure wexp
       REAL FUNCTION WEXSVP(T)

       IMPLICIT NONE

       REAL T
       DOUBLE PRECISION TEMP,LOGP

       TEMP = DBLE(T)
       LOGP = -0.58002206D+4 / TEMP
     $       + 0.13914993D+1
     $       - 0.48640239D-1 * TEMP
     $       + 0.41764768D-4 * TEMP**2
     $       - 0.14452093D-7 * TEMP**3
     $       + 0.65459673D+1 * DLOG(TEMP)
       WEXSVP = SNGL( 1.0D-2 * DEXP( LOGP ) )
C
C      Hyland & Wexler, 1983, equation for SVP over ice
C           log Pi =  -0.56745359E+4 / T
C                    + 0.63925247E+1
C                    - 0.96778430E-2  T
C                    + 0.62215701E-6  T^2
C                    + 0.20747825E-8  T^3
C                    - 0.94840240E-12 T^4
C                    + 0.41635019E+1 log(T)
C       with T in [K] and Pi in [Pa]
C
       RETURN
       END

c************************************************************************
c copied from klayers/grav.f
!CALL PROTOCOL:
C    GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      LAT     latitude                    degrees          
C    REAL      LON     longitude                   degrees
C    REAL      WINDE   wind velecity east          m/s
C    REAL      WINDN   wind velocity north         m/s
C    REAL      Z       altitude                    m

!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      GRAV    Earth gravity               m/s^2

!INPUT/OUTPUT PARAMETERS: none

!RETURN VALUES: none

!PARENT(S):
C    INTLEV

!ROUTINES CALLED: none

!FILES ACCESSED: none

!COMMON BLOCKS: none

!DESCRIPTION:
C    Function to calculate Earth gravity (gravitation plus
C    centripetal acceleration) for points in the atmosphere.
C
C    It calculates surface gravity using an equation given in
C    "American Institute of Physics Handbook", 1963, page 2-102.
C    This equation is essentially a variant of the International
C    Gravity Formula with an extra term for longitude (which is
C    very minor).
C
C    Centripetal acceleration is tangental velocity squared over the
C    radius.
C
C    Gravitation is the gravitational constant times the Earth's mass
C    divided by the square of the radius.
C
C    Gravity at any point in the atmosphere is simply surface gravity
C    plus the change in gravitation and centripetal acceleration at
C    that point compared to the surface.

!ALGORITHM REFERENCES: see DESCRIPTION

!KNOWN BUGS AND LIMITATIONS: none

!ROUTINE HISTORY:
C    Date     Programmer        Comments
C------------ ----------------- ----------------------------------------
C Mar  8 1995 Scott Hannon/UMBC created
C Jun 23 1995 Scott Hannon      Correct some comments
C 28 Sep 2007 Scott Hannon      Add more centripetel term comments

!END====================================================================

C      =================================================================
       REAL FUNCTION GRAV_EARTH(Z, WINDE, WINDN, LAT, LON)
C      =================================================================
C
C-----------------------------------------------------------------------
C      IMPLICIT NONE
C-----------------------------------------------------------------------
       IMPLICIT NONE

C-----------------------------------------------------------------------
C      INCLUDE FILES
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      EXTERNAL FUNCTIONS
C-----------------------------------------------------------------------
C      none

C-----------------------------------------------------------------------
C      ARGUMENTS
C-----------------------------------------------------------------------

       REAL Z, WINDE, WINDN, LAT, LON

C-----------------------------------------------------------------------
C      LOCAL VARIABLES
C-----------------------------------------------------------------------
C
       REAL G_SUR, R, COSLT, COSLT2, SINLT2, SIN2LT, COSLON,
     $    LTRAD, C_SUR, C_Z, RTOT, GRAVZ
C
C      Constants for (1 - b^2/a^2) with
C      a = 6.378388E+6 m = equatorial radius, and
C      b = 6.356911E+6 m = polar radius.
       REAL B2, ABTERM
C
C      Constants for normal gravity equation
C      (see "American Institute of Physics Handbook", 1963, pg 2-102) 
       REAL G0
       REAL C1, C2, C3
C
C      Constants for pi/180, 2*pi, and Earth's rotational speed
C      of w=1/86400 rev/s
       REAL PI180,PI2, W
C 
       DATA B2, ABTERM /4.041031E+13, 6.724285E-3/
       DATA G0 /9.780455/
       DATA C1,C2,C3 /5.30157E-3, -5.85E-6, 6.40E-6/
       DATA PI180,PI2,W /1.7453293E-2, 6.28318531, 1.1574074E-5/

C-----------------------------------------------------------------------
C      SAVE STATEMENTS
C-----------------------------------------------------------------------
C      none

C***********************************************************************
C***********************************************************************
C      EXECUTABLE CODE begins below
C***********************************************************************
C***********************************************************************
C
C      Calculate longitude term
C      Add offset of 18 degrees, double it, convert to radians, and
C      take the cosine
       COSLON=COS( PI180*2.0*( LON + 18.0 ) )
C
C      Calculate the latitude terms
C      Convert Latitude into radians
       LTRAD=PI180*LAT
C      Calculate sine and cosine terms
       COSLT = COS(LTRAD)
       COSLT2 = COSLT**2
       SINLT2 = ( SIN(LTRAD ) )**2
       SIN2LT = ( SIN( 2.0*LTRAD ) )**2
C
C      Calculate the Earth's radius at this latitude
       R = SQRT( B2/( 1.0 - COSLT2*ABTERM ) )
C
C      Calculate total distance from Earth's center
       RTOT = R + Z
C
C      Calculate gravity at the Earth's surface
       G_SUR = G0*( 1.0 + C1*SINLT2 + C2*SIN2LT + C3*COSLT2*COSLON )
C
C      Calculate the centripetal term at the Earth's surface
C      Note: the centripetal acceleration due to Earth's rotation
C      is in a direction perpendicular to the Earth's rational axis,
C      but for this gravity function we are only interested in the
C      component parallel to the radial direction.
       C_SUR = COSLT2*R*(PI2*W)**2
C
C      Calculate the centripetal term at altitude z (with wind)
       C_Z = ( ( PI2*RTOT*COSLT*W + WINDE )**2 + (WINDN)**2 )/RTOT
C
C      Calculate the change in gravitation with altitude
       GRAVZ=(G_SUR + C_SUR)*(1.0 - R**2/RTOT**2)
C
       GRAV_EARTH = G_SUR + (C_SUR - C_Z) - GRAVZ
C
       RETURN
       END

c************************************************************************
c this reads in the reference levels profile
      SUBROUTINE ReadRefProf_Levels(PLEV_KCARTADATABASE_AIRS,iGasID,iNumLevsx,raRx110Press,raRx110Temp,raRx110MR)
      
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iGasID
      REAL PLEV_KCARTADATABASE_AIRS(kProfLayer+1)
c output
      INTEGER iNumLevsx
      REAL raRx110Temp(kProfLayer+10),raRx110MR(kProfLayer+10),raRx110Press(kProfLayer+10)
c local
      REAL PLEVx110(kProfLayer+10)
      REAL rLat,raJunk(kProfLayer+10)
      CHARACTER*80 caPFname,caStr,caComment
      INTEGER iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec

c first extend PLEV_KCARTADATABASE_AIRS a little above its lwest value, so can do interps well
      DO iI=1,kProfLayer+1
        PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(iI)
      END DO
      rLat = PLEV_KCARTADATABASE_AIRS(kProfLayer+1)/20
      DO iI=kProfLayer+2,kProfLayer+10
        PLEVx110(iI) = PLEV_KCARTADATABASE_AIRS(kProfLayer+1) - rLat*(iI-(kProfLayer+1))
      END DO

c      caPFname = '/home/sergio/KCARTA/INCLUDE/glatm_16Aug2010.dat'
      caPFname = kcaLevsRefProf

      IF ((kAFGLProf .LT. 1) .OR. (kAFGLProf .GT. 6)) THEN
        write(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
        CALL DoStop
      END IF

      iErr = 0
      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
          CALL DoSTOP
      ENDIF
 1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
      kProfileUnitOpen=1

      iAFGL = 0
      iFoundGas = -1
      
  10  CONTINUE
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 10   !! keep reading comments
      ELSE
        READ (caStr,*) iNumLevsx
      END IF
      IF (iNumLevsx .GT. kProfLayer+1) THEN
        write(kStdErr,*) 'oops iNumLevsx .GT. kProfLayer+1 ',iNumLevsx,kProfLayer+1
        CALL DoStop
      END IF

  20  CONTINUE
      iAFGL = iAFGL + 1

  30  CONTINUE  
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 30   !! keep reading comments
      ELSE  
        caComment = caStr          
      END IF

  40  CONTINUE
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 40   !! keep reading comments
      ELSE
        READ (caStr,*) rLat
      END IF
c      write(kStdWarn,*) 'iAFGL Profile ',iAFGL,' rLat ',rLat, ' : ',caComment

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Altitude(km)
      READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

      IF (iAFGL .EQ. kAFGLProf) THEN
c        write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! Press(mb)
        READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)

        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! Temp(K)
        READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)

c        write(kStdWarn,*) '   need to find profile for gasID ',iGasID
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
        IF (iG .EQ. iGasID) THEN
c          write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
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

      IF ((iAFGL .EQ. kAFGLProf) .AND. (iFoundGas .GT. 0)) THEN 
        GOTO 60   !! found the AFGL prof and found the gas; done!
      ELSEIF ((iAFGL .EQ. 6) .AND. (iFoundGas .LT. 0)) THEN 
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! modend
        GOTO 50   !! found the AFGL prof but not found the gas
      ELSEIF ((iAFGL .LT. 6) .OR. (iFoundGas .LT. 0)) THEN 
        !! either did not find the gas or the AFGL prof
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! modend
        GOTO 20
      END IF

 50   CONTINUE 
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! constituent profs
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! mingas

      iG = 7
 55   CONTINUE
      iG = iG + 1
      IF (iG .EQ. iGasID) THEN
c        write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
        iFoundGas = +1
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
      ELSE
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
      END IF
      IF ((iG .LT. 28) .AND. (iFoundGas .LT. 0)) THEN
        GOTO 55
      END IF
      
      IF (iFoundGas .LT. 0) THEN
        !! bweh try xsec gases
 123    CONTINUE
        READ(iIOUN2,1080) caStr   !! comment
        IF (caStr(1:1) .EQ. '!') THEN
          GOTO 123
        ELSEIF (caStr(1:6) .EQ. 'DATEND') THEN
          GOTO 60       !!! oh oh end of file, give up!!!!!              
        ELSE
          READ (caStr,*) iXsec
          IF (iXsec .EQ. iGasID) THEN 
            iFoundGas = 1
            READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)          
          ELSE
            READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)          
            GOTO 123
          END IF
        END IF
      END IF

 60   CONTINUE
      CLOSE(iIOUN2) 
      kProfileUnitOpen=-1
 1080 FORMAT(A80)

      IF (iFoundGas .LT. 0) THEN
        !! finally give up
        write(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
        CALL DOStop
      END IF

c*************************
c finally interp these onto the AIRS pressure levels
      DO iI = 1,iNumLevsx
        raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
c        print *,iI,raRx110Press(iI),raRx110Temp(iI),raRx110MR(iI)
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
      END

c************************************************************************
c this subroutine does the klayers integrations, integrals over P
c so this is more similar to kLAYERs, CRTM
c this is newer code; it ignores surface as "lowest" level and only uses the user level info
c this is for AIRS 101 levels

c >>>>>>>> raPout comes in mb
c >>>>>>>> and goes out in N/m2

      SUBROUTINE DoIntegrateLevels2Layers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,
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
      REAL raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
      REAL zWoo,rConvertQ,grav_earth,wexsvp,rRH
      INTEGER iWoo

      rConvertQ = 6.022141e23/1e4     ! multiply by this to change moles/m2 to molecules/cm2
      rConvertQ = kAvog/1000/1.0e4    ! multiply by this to change moles/m2 to molecules/cm2
      rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor
      
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

      write(kStdWarn,*) '  '
      write(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'      

c >>>>>>>>>>
      !! set up temporary arrays/matrices for interping onto sublevels
      iUpperLev = kMaxLayer-iLowestLev+1
      !! set levels 1 .. iNumLevs
      DO iL = 1,iUpperLev
        raXYZPress(iL) = raPX(iL)
        raXYZTemp(iL) = raTX(iL)
        DO iG = 1,iNumGases
          raaXYZ_MR(iL,iG) = raaG_MRX(iL,iG)
          IF (iG .EQ.  1) raXYZ_MRwater(iL) = raaG_MRX(iL,iG)
        END DO
      END DO
c >>>>>>>>>>

c now start integrating
c lowest layer is special case

c http://cimss.ssec.wisc.edu/~paulv/
c http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
c see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
c     /home/sergio/KCARTA/DOC/sci_klayers.txt

 1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
 
      z = rHSurf !! need to kludge this
      zWoo = rHSurf
      iWoo = -1

      DO iL = iLowestLev,iLowestLev
        raPout(iL) = log(PLEV_KCARTADATABASE_AIRS(iL)/PLEV_KCARTADATABASE_AIRS(iL+1))
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

        !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
        dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))  
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                             !! current sublev press
        rT_n = rT                             !! current sublev temp
        rMR_water_n = raXYZ_MRwater(1)        !! current water MR
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP,rMR_water_n,1) !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

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
            DO iCnt = 1,iUpperLev
              raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
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

          rP_np1 = log(rP_n) - dlnp
          rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
          CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

          IF ((rP_n .LE. rPSurf) .AND. (iWoo .LT. 0)) THEN
            iWoo = +1
            zWoo   = z
          END IF

        END DO

        raTout(iL)      = rTsum/rWgt                 ! <<< need TAV
        raAmountOut(iL) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
        raZout(iL+1)    = z

c        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
c     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
      END DO

      ! displace everything  z
      iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
      raZout(iL)   = -(zWoo - z0)
      raZout(iL+1) = +(z-zWoo)
      write(kStdWarn,*) '  need to displace lowest layer heights'
      write(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
      write(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',PLEV_KCARTADATABASE_AIRS(iL),rPSurf,
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
        raPout(iL) = raPout(iL)/100.0    !!! change N/m2 to mb
	
        iJ = iL - iLowestLev + 1
        rP = PLEV_KCARTADATABASE_AIRS(iL)
	rP = rP/100.0    !!! change N/m2 to mb	
        rT = raXYZTemp(iJ)
        
        amount = 0.0
        rTsum  = 0.0
        rWgt = 0.0

        !! no need to divide by 100 since it cancels   log(a/c)-log(b/c) = log(a/c/b/c) = loag(a/c)
        dlnp = log(PLEV_KCARTADATABASE_AIRS(iL)) - log(PLEV_KCARTADATABASE_AIRS(iL+1))
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                              !! current sublev press
        rT_n = rT                              !! current sublev temp
        rMR_water_n = raXYZ_MRwater(iJ)        !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

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
            DO iCnt = 1,iUpperLev
              raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
            END DO
            CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1) 
            CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1) 
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
          rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
          CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

        END DO

        raTout(iL)      = rTsum/rWgt             ! <<< need TAV
        raAmountOut(iL) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
        raZout(iL+1)    = z

c        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,PLEV_KCARTADATABASE_AIRS(iL+1),rT,raTout(iL),
c     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
      END DO
      
      write(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
      !! raZout is in meters
      !! raaQout is in moles/m2
      DO iL = iLowestLev,kProfLayer
        rThick = abs(raZout(iL) - raZout(iL+1)) * 100.0
        rT = raTout(iL)
	raPout(iL) = raPout(iL) * 100.0    !! change mb to N/m2
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
        write(kStdWarn,*) 'Checking Relative Humidity (integrate over AIRS 101 levels)'
        write(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
        write(kStdWarn,113) '------------------------------------------------------------------------------------------'
        DO iL = iLowestLev,kProfLayer
          z    = (raZout(iL) + raZout(iL+1))/2/1000
          zWoo = (raZout(iL+1) - raZout(iL))/1000
          rP   = raPout(iL)            !! N/m2
          rT   = raTout(iL)
          rPP  = raaPartPressOut(iL,1) !! N/m2
          rSVP = wexsvp(rT) * 100      !! mb --> N/m2
          rRH  = rPP/rSVP*100.0
          IF (rRH .LE. 100.0) THEN
            write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),
     $ (raaQout(iL,iG),iG=1,6)
          ELSE
            write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),
     $ (raaQout(iL,iG),iG=1,6),' ***** '
          END IF
        END DO
      END IF
 111  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '))
 112  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(ES9.3,' '),A7)
 113  FORMAT(A90)

c now find pressure output corresponding to HGT output from LBLRTM
      IF (kRTP .EQ. -20) THEN
        IF (raRTP_TxtInput(6) .GT. 0) THEN
          !!! input level boundaries in km, change to mb
          DO iL = iLowestLev,kProfLayer
            raJunk(iL-iLowestLev+1) = raZout(iL)          
            raJunk2(iL-iLowestLev+1) = raPout(iL)          
          END DO
          CALL r_sort_linear1(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo,1)
          write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
          raRTP_TxtInput(6) = zWoo/100.0  !! mb
        ELSEIF (raRTP_TxtInput(6) .LT. 0) THEN
          !!! input level boundaries in mb     
          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
        END IF
      END IF

      write(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
      
      RETURN
      END

c************************************************************************
c this subroutine does the klayers integrations, integrals over P
c so this is more similar to kLAYERs, CRTM
c this is newer code; it ignores surface as "lowest" level and only uses the user level info
c this is for USER DEFINED levels

c >>>>>>>> raPout comes in mb
c >>>>>>>> and goes out in N/m2

      SUBROUTINE DoIntegrateLevels2UserLayers_wrtP(rHSurf,rPSurf,rTSurf,iLowestLev,iNumGases,iaG,rLat,rLon,
     $               raPBnd,raAlt,iZbnd,rfracBot,
     $               raPX,raTX,raaG_MRX,raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL rHSurf,rPSurf,rTSurf,raPBnd(2*kProfLayer)
      REAL raAlt(2*kProfLayer),rfracBot,rLat,rLon
      INTEGER iLowestLev,iNumGases,iaG(kMaxGas),iZbnd
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
      REAL raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
      REAL zWoo,rConvertQ,grav_earth,wexsvp,rRH
      INTEGER iOffset,iWoo

      rConvertQ = 6.022141e23/1e4     ! multiply by this to change moles/m2 to molecules/cm2
      rConvertQ = kAvog/1000/1.0e4    ! multiply by this to change moles/m2 to molecules/cm2
      rConvertQ = rConvertQ * 100.0   !but if input p is in mb, we need to change to N/m2 by this factor
      
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

      write(kStdWarn,*) '  '
      write(kStdWarn,*) '  >>>>>>>>>>>> doing the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'      

c >>>>>>>>>>
      !! set up temporary arrays/matrices for interping onto sublevels
      iUpperLev = iZbnd
      !! set levels 1 .. iNumLevs
      DO iL = 1,iUpperLev
        raXYZPress(iL) = raPX(iL)
        raXYZTemp(iL) = raTX(iL)
        DO iG = 1,iNumGases
          raaXYZ_MR(iL,iG) = raaG_MRX(iL,iG)
          IF (iG .EQ.  1) raXYZ_MRwater(iL) = raaG_MRX(iL,iG)
        END DO
c      print *,iNumGases,iL,raPX(iL),raTX(iL),raaG_MRX(iL,1),raaG_MRX(iL,7),raaG_MRX(iL,22)	
      END DO

c >>>>>>>>>>

c now start integrating
c lowest layer is special case

c http://cimss.ssec.wisc.edu/~paulv/
c http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
c see /home/sergio/KCARTA/DOC/NotesOnAtmosphericProfilesAndQuantities.pdf
c     /home/sergio/KCARTA/DOC/sci_klayers.txt

 1234 FORMAT(I3,' ',F10.4,' ' ,I3,6(' ',F10.4),1(' ',E10.4),3(' ',F10.4))
 
      z = rHSurf !! need to kludge this
      zWoo = rHSurf
      iWoo = -1

c      DO iL=1,iZbnd
c        print *,'sergio1',iL,raPBnd(iL)
c      end do
      
      iOffSet = kProfLayer-(iZbnd-1)
      DO iL = iLowestLev,iLowestLev
        raPout(iL+iOffset) = log(raPBnd(iL)/raPBnd(iL+1))
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
        dlnp = log(raPBnd(iL)) - log(raPBnd(iL+1))  
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                             !! current sublev press
        rT_n = rT                             !! current sublev temp
        rMR_water_n = raXYZ_MRwater(1)        !! current water MR
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP,rMR_water_n,1) !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

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
c          print *,iLoop,raPBnd(iL),rP,amount*rConvertQ,(amount+damount)*rConvertQ,SpHeatDryAir,kMGC
          amount = amount + damount

          DO iG = 1,iNumGases
            DO iCnt = 1,iUpperLev
              raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
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

          rP_np1 = log(rP_n) - dlnp
          rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
          CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

          IF ((rP_n .LE. rPSurf) .AND. (iWoo .LT. 0)) THEN
            iWoo = +1
            zWoo   = z
          END IF

        END DO

        raTout(iL+iOffSet)      = rTsum/rWgt                 ! <<< need TAV
        raAmountOut(iL+iOffSet) = amount*rConvertQ           ! change moles/m2 to molecules/cm2
        raZout(iL+1+iOffSet)    = z
        
c        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL+iOffSet),
c     $        raAmountOut(iL+iOffSet),raZout(iL+1+iOffSet),raPout(iL+iOffSet),z/1000
      END DO

      ! displace everything  z
      iL = iL - 1    !! recall on exiting the loop, iL is actually incremented by 1
      raZout(iL+iOffSet)   = -(zWoo - z0)
      raZout(iL+1+iOffSet) = +(z-zWoo)
      write(kStdWarn,*) '  need to displace lowest layer heights'
      write(kStdWarn,*) '  Lowest Level (m), rHSurf (m) Upper Levl(m) = ',-(zWoo - z0),z0,+(z-zWoo)
      write(kStdWarn,*) '  plowest,pSurf,pLowest+1 (mb) ',raPBnd(iL),rPSurf,
     $                                                   raPBnd(iL+1)
      z = z - zWoo

ccc      !no need to adjust this amount for book-keeping, make consistent with eg n_rad_jac_scat.f since we did FULL layer
ccc      raAmountOut(iL) = raAmountOut(iL)/rFracBot 
ccc      DO iG = 1,iNumGases
ccc        raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
ccc      END DO
      
c go to TOA as defined by user

      DO iL = iLowestLev + 1,iZbnd-1
        z0 = z
        !! compute avg pressure of layer
        raPout(iL+iOffSet) = log(raPBnd(iL)/raPBnd(iL+1))
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
        dlnp = log(raPBnd(iL)) - log(raPBnd(iL+1))
        dlnp = dlnp / (iNFine)

        !! information for current (sub)level
        rP_n = rP                              !! current sublev press
        rT_n = rT                              !! current sublev temp
        rMR_water_n = raXYZ_MRwater(iJ)        !! current water MR

        !! information for next (sub)level
        rP_np1 = log(rP_n) - dlnp
        rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
        CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
        CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

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
            DO iCnt = 1,iUpperLev
              raJunk(iCnt) = raaXYZ_MR(iCnt,iG)
            END DO
            CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_n,  rMR_n,  1) 
            CALL r_sort_loglinear1(raXYZPress,raJunk,iUpperLev,rP_np1,rMR_np1,1) 
            raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + damount * (rMR_n + rMR_np1)/2 * rConvertQ
c            q_n   = rP_n  /rT_n  /kMGC * dz * rMR_n
c            q_np1 = rP_np1/rT_np1/kMGC * dz * rMR_np1
c            raaQout(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) + (q_n + q_np1)/2 * rConvertQ
          END DO

          !!! now update for next iteration
          z = z + dz
          rP_n = rP_np1
          rT_n = rT_np1
          rMR_water_n = rMR_water_np1
          rP = rP_n

          rP_np1 = log(rP_n) - dlnp
          rP_np1 = exp(rP_np1)                                                                !! next sublev pressure
          CALL r_sort_loglinear1(raXYZPress,raXYZTemp,    iUpperLev,rP_np1,rT_np1,       1)   !! next sublev temp
          CALL r_sort_loglinear1(raXYZPress,raXYZ_MRwater,iUpperLev,rP_np1,rMR_water_np1,1)   !! next sublev MRw

        END DO

        raTout(iL+iOffSet)      = rTsum/rWgt             ! <<< need TAV
        raAmountOut(iL+iOffSet) = amount * rConvertQ     ! change moles/m2 to molecules/cm2
        raZout(iL+1+iOffSet)    = z

c        write(*,1234) iL,(z-z0),iJ,raTX(iJ),rTSurf,rP,raPBnd(iL+1),rT,raTout(iL),
c     $        raAmountOut(iL),raZout(iL+1),raPout(iL),z/1000
      END DO

      write(kStdWarn,*) 'Topmost press level is at z = ',z,' meters'
      DO iL = 1,iOffSet-1
        raPout(iL) = 0.0
	raTout(iL) = 0.0
	raZout(iL+1) = raZout(iZbnd)
	raAmountOut(iL) = 0.0
      END DO

      !! raZout is in meters
      !! raaQout is in moles/m2
      DO iL = iLowestLev,iZbnd-1
        rThick = abs(raZout(iL+iOffSet) - raZout(iL+1+iOffSet)) * 100.0
        rT = raTout(iL+iOffSet)
	raPout(iL+iOffSet) = raPout(iL+iOffSet) * 100.0    !! change mb to N/m2
        DO iG = 1,iNumGases
          raaPartPressOut(iL+iOffSet,iG) = raaQout(iL+iOffSet,iG) / raAmountOut(iL+iOffSet) * raPout(iL+iOffSet)
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

c      IF ((iLowestLev-1) .GE. 1) THEN
        DO iL = 1,iOffSet-1
          raPout(iL) = log(raPBnd(iL)/raPBnd(iL+1))
          raPout(iL) = (raPBnd(iL) - raPBnd(iL+1))/raPout(iL)
          raTout(iL) = kStempMin
          DO iG = 1,iNumGases
            raaPartPressOut(iL,iG) = 0.0
            raaQout(iL,iG) = 0.0
            raaG_MRX(iL,iG) = 0.0
          END DO
        END DO
c      END IF

      IF (iZbnd+iOffSet .LE. kProfLayer) THEN
        DO iL = iZbnd,kProfLayer
          raPout(iL) = 0.0
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
        write(kStdWarn,*) 'Checking Relative Humidity (integrate over user levels)'
        write(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
        write(kStdWarn,113) '------------------------------------------------------------------------------------------'
        DO iL = iLowestLev,iZbnd-1
          z    = (raZout(iL+iOffset) + raZout(iL+1+iOffset))/2/1000
          zWoo = (raZout(iL+1+iOffset) - raZout(iL+iOffset))/1000
          rP   = raPout(iL+iOffset)            !! N/m2
          rT   = raTout(iL+iOffset)
          rPP  = raaPartPressOut(iL+iOffset,1) !! N/m2
          rSVP = wexsvp(rT) * 100      !! mb --> N/m2
          rRH  = rPP/rSVP*100.0
          IF (rRH .LE. 100.0) THEN
            write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet),
     $ (raaQout(iL+iOffSet,iG),iG=1,6)
          ELSE
            write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL+iOffSet),
     $ (raaQout(iL+iOffSet,iG),iG=1,6),' BAD RH!'
          END IF
        END DO
      END IF
 111  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '))
 112  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '),A7)
 113  FORMAT(A90)

c now find pressure output corresponding to HGT output from LBLRTM
      IF (kRTP .EQ. -20) THEN
        IF (raRTP_TxtInput(6) .GT. 0) THEN
          !!! input level boundaries in km, change to mb
          DO iL = iLowestLev,kProfLayer
            raJunk(iL-iLowestLev+1) = raZout(iL)          
            raJunk2(iL-iLowestLev+1) = raPout(iL)          
          END DO
          CALL r_sort_linear1(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo,1)
          write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
          raRTP_TxtInput(6) = zWoo/100.0  !! mb
        ELSEIF (raRTP_TxtInput(6) .LT. 0) THEN
          !!! input level boundaries in mb     
          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
        END IF
      END IF

      write(kStdWarn,*) '  >>>>>>>>>>>> finished the Levels --> Layers intergation >>>>>>>>>>>>>>>>>>>'
      
      RETURN
      END

c************************************************************************
      SUBROUTINE  ReadInput_LVL_Profile(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
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
      INTEGER iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iReadInOK
      REAL raX(kMaxGas),rX,rP,rT,raP_Nm2(2*kProfLayer),rPRefMin,rPRefMax
      CHARACTER*80 caStr


      rPRefMax = 1100.0 !! 1100  mb is AIRS RefPress LowerBdry      see SUBR ReadRefProf_Levels
      rPRefMin = 3.0e-3 !! 0.005 mb extended to 0.003 mb = 0.3 N/m2 see SUBR ReadRefProf_Levels
      
      rPmin = +1.0e6
      rPmax = -1.0e+6

      write(kSTdWarn,*) 'Reading in user supplied TXT LVLS file .....'

      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
 1070     FORMAT('ERROR! number ',I5,' opening PRFILE path LEVELS profile file'
     $            ,/,A80)
          CALL DoSTOP
      ENDIF
      kProfileUnitOpen=1

      READ (iIOUN2,5030,ERR=13,END=13) caStr

      READ (iIOUN2,*) iNumLevs
      IF (iNumLevs .GT. 2*kProfLayer) THEN
        write(kStdErr,*) 'iNumLevs .GT. 2*kProfLayer',iNumLevs,kProfLayer
        CALL DoStop
      END IF

      READ (iIOUN2,*) rPSurf,rTSurf,rHSurf
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

      raRTP_TxtInput(1) = rPSurf
      raRTP_TxtInput(2) = rTSurf
      raRTP_TxtInput(3) = rHSurf
      
      READ (iIOUN2,*) rYear,rLat,rLon

      READ (iIOUN2,*) iNumGases
      IF (iNumGases .GT. kMaxGas) THEN
        write(kStdErr,*) 'iNumGases .GT. kMaxGas',iNumGases,kMaxGas
        CALL DoStop
      END IF
       
      READ (iIOUN2,*) (iaG(iJ),iJ=1,iNumGases)
      READ (iIOUN2,*) (iaGasUnits(iJ),iJ=1,iNumGases)
      
      write(kStdWarn,*) 'Input levels profile : ',iNumLevs,' levels for ',iNumGases,' gases'
      write(kStdWarn,*) 'PSurf = ',rPSurf,' mb;  TSurf = ',rTSurf,' K '
      iReadInOK = 0
      DO iL = 1,iNumLevs
        READ (iIOUN2,*) rP,rT,(raX(iJ),iJ=1,iNumGases)
	IF ((rP .LE. rPRefMax) .AND. (rP .GE. rPRefMin)) THEN
	  iReadInOK = iReadInOK + 1
          raP_Nm2(iReadInOK) = rP * 100.0  !! change from mb to N/m2
          raP(iReadInOK) = rP
          raT(iReadInOK) = rT
          IF (rPmax .LE. raP(iReadInOK)) rPmax = raP(iReadInOK)
          IF (rPmin .GT. raP(iReadInOK)) rPmin = raP(iReadInOK)
          DO iJ = 1,iNumGases
            raaG_MR(iReadInOK,iJ) = raX(iJ)
          END DO
  	  write(kStdWarn,5040) iReadInOK,raP(iReadInOK),raT(iReadInOK),raaG_MR(iReadInOK,1)
	ELSE
  	  write(kStdWarn,5041) iL,rP,rT
	END IF
      END DO
      
 13   CONTINUE
      CLOSE(iIOUN2) 
      kProfileUnitOpen = -11

      IF (iNumLevs .NE. iReadInOK) THEN
        WRITE(kStdWarn,*) 'Input textfile had ',iNumLevs,' levels of which ',iReadInOK,' lay between ',rPRefMax,rPRefMin,' mb'
      ELSE
        WRITE(kStdWarn,*) 'All ',iNumLevs,' input levels in textfile lay between ',rPRefMax,rPRefMin,' mb'
      END IF
      iNumLevs = iReadInOK
      
 5030 FORMAT(A80)
 5040 FORMAT('iReadInOK, rP (mb), rT (K), raG1 = ',I3,' ',E10.5,' ',F8.3,' ',E10.5)
 5041 FORMAT('Bad (too low or too high P) rP (mb), rT (K), raG1 = ',I3,' ',E10.5,' ',F8.3) 
 
      write(kStdWarn,*)'   KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
      write(kStdWarn,*)'   kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
      write(kStdWarn,*)' input file : spres (mb)/Height(km)     = ',rPSurf,rHSurf

      rPSurf = rPSurf * 100.0
      rHSurf = rHSurf * 1000.0

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
        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP_Nm2,raT,raaG_MR)
        DO iL = 1,iNumLevs
          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
          iaGasUnits(iG) = 12
        END DO
      END DO         

      RETURN
      END

c************************************************************************
      SUBROUTINE ReadRefProf_Levels2(iGasID,raP,iNumLevsIN,raX)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c input
      INTEGER iGasID,iNumLevsIN
      REAL raP(2*kProfLayer)
c output
      REAL raX(2*kProfLayer)

c local
      REAL raRx110Temp(2*kProfLayer),raRx110MR(2*kProfLayer),raRx110Press(2*kProfLayer)
      REAL rLat,raJunk(kProfLayer*2)
      CHARACTER*80 caPFname,caStr,caComment
      INTEGER iAFGL,iProf,iIOUN2,iERRIO,iErr,iI,iG,iFoundGas,iXsec,iNumLevsx

      caPFname = kcaLevsRefProf

      IF ((kAFGLProf .LT. 1) .OR. (kAFGLProf .GT. 6)) THEN
        write(kStdErr,*) 'Need 1 <= kAFGLProf <= 6, but kAFGLProf = ',kAFGLProf
        CALL DoStop
      END IF

      iErr = 0
      iIOUN2 = kProfileUnit
      OPEN(UNIT=iIOun2,FILE=caPfname,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
      IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caPfname
          CALL DoSTOP
      ENDIF
 1070 FORMAT('ERROR! number ',I5,' opening GLATM file ',A80)
      kProfileUnitOpen=1

      iAFGL = 0
      iFoundGas = -1
      
  10  CONTINUE
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 10   !! keep reading comments
      ELSE
        READ (caStr,*) iNumLevsx
      END IF
      IF (iNumLevsx .GT. 2*kProfLayer) THEN
        write(kStdErr,*) 'oops iNumLevsx .GT. 2*kProfLayer ',iNumLevsx,kProfLayer*2
        CALL DoStop
      END IF

  20  CONTINUE
      iAFGL = iAFGL + 1

  30  CONTINUE  
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 30   !! keep reading comments
      ELSE  
        caComment = caStr          
      END IF

  40  CONTINUE
      READ(iIOUN2,1080) caStr
      IF (caStr(1:1) .EQ. '!') THEN 
        GOTO 40   !! keep reading comments
      ELSE
        READ (caStr,*) rLat
      END IF
c      write(kStdWarn,*) 'iAFGL Profile ',iAFGL,' rLat ',rLat, ' : ',caComment

      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! Altitude(km)
      READ(iIOUN2,*) (rajunk(iI),iI=1,iNumLevsx)

      IF (iAFGL .EQ. kAFGLProf) THEN
c        write(kStdWarn,*) 'Reading in P/T for kAFGLProf profile ',iAFGL
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! Press(mb)
        READ(iIOUN2,*) (raRx110Press(iI),iI=1,iNumLevsx)

        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! Temp(K)
        READ(iIOUN2,*) (raRx110Temp(iI),iI=1,iNumLevsx)

c        write(kStdWarn,*) '   need to find profile for gasID ',iGasID
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
        IF (iG .EQ. iGasID) THEN
c          write(kStdWarn,*) 'Reading in GasID ',iG,' profile from kAFGLprof profile ',iAFGL
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

      IF ((iAFGL .EQ. kAFGLProf) .AND. (iFoundGas .GT. 0)) THEN 
        GOTO 60   !! found the AFGL prof and found the gas; done!
      ELSEIF ((iAFGL .EQ. 6) .AND. (iFoundGas .LT. 0)) THEN 
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! modend
        GOTO 50   !! found the AFGL prof but not found the gas
      ELSEIF ((iAFGL .LT. 6) .OR. (iFoundGas .LT. 0)) THEN 
        !! either did not find the gas or the AFGL prof
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! modend
        GOTO 20
      END IF

 50   CONTINUE 
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! constituent profs
      READ(iIOUN2,1080) caStr   !! comment
      READ(iIOUN2,1080) caStr   !! mingas

      iG = 7
 55   CONTINUE
      iG = iG + 1
      IF (iG .EQ. iGasID) THEN
c        write(kStdWarn,*) 'Reading in minor absorbing GasID ',iG,' profile from AFGL profile ',iAFGL
        iFoundGas = +1
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)
      ELSE
        READ(iIOUN2,1080) caStr   !! comment
        READ(iIOUN2,1080) caStr   !! gas name
        READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)
      END IF
      IF ((iG .LT. 28) .AND. (iFoundGas .LT. 0)) THEN
        GOTO 55
      END IF
      
      IF (iFoundGas .LT. 0) THEN
        !! bweh try xsec gases
 123    CONTINUE
        READ(iIOUN2,1080) caStr   !! comment
        IF (caStr(1:1) .EQ. '!') THEN
          GOTO 123
        ELSEIF (caStr(1:6) .EQ. 'DATEND') THEN
          GOTO 60       !!! oh oh end of file, give up!!!!!              
        ELSE
          READ (caStr,*) iXsec
          IF (iXsec .EQ. iGasID) THEN 
            iFoundGas = 1
            READ(iIOUN2,*) (raRx110MR(iI),iI=1,iNumLevsx)          
          ELSE
            READ(iIOUN2,*) (raJunk(iI),iI=1,iNumLevsx)          
            GOTO 123
          END IF
        END IF
      END IF

 60   CONTINUE
      CLOSE(iIOUN2) 
      kProfileUnitOpen=-1
 1080 FORMAT(A80)

      IF (iFoundGas .LT. 0) THEN
        !! finally give up
        write(kStdErr,*) 'read 6 AFGL profs, gases 1-7,8-28, and XSEC gases but did not find gas OOPS',iGasID
        CALL DOStop
      END IF

c*************************
c finally interp these onto the raP pressure levels
      DO iI = 1,iNumLevsx
        raRx110Press(iI) = raRx110Press(iI) * 100.0  !! change mb --> N/m2
c        print *,iI,raRx110Press(iI),raRx110Temp(iI),raRx110MR(iI)
      END DO

c      CALL r_sort_loglinear(raRx110Press,raRx110Temp,iNumLevsx,PLEVx110,raJunk,2*kProfLayer)
c      DO iI = 1,2*kProfLayer
c        raRx110Temp(iI) = raJunk(iI)
c      END DO
      CALL r_sort_loglinear(raRx110Press,raRx110MR,iNumLevsx,raP,raJunk,iNumLevsIN)
      DO iI = 1,iNumLevsIN
        raX(iI) = raJunk(iI)
      END DO
c      iNumLevsx = 2*kProfLayer

      RETURN
      END
c************************************************************************
