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
