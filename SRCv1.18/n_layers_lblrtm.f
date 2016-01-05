c Copyright 2014
c University of Maryland Baltimore County 
c All Rights Reserved

c this file does the LBLRTM TAPE5/TAPE6 

c see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
c       http://www.ssec.wisc.edu/~paulv/Fortran90/AtmProfile/Modules.html
c       https://svn.ssec.wisc.edu/repos/uwphysret/trunk/mfiles/compute_F.m
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
c
c      LBLRTM     iaGasUnits(iG) = 12   !!! assume hardcoded VMR
c      LBLRTM     pressures are in mb, which is what klayers wants
c
c      rP(min/max)KCarta is in N/m2 which is x100 what it would be in mb
c
c************************************************************************
c http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
c this one takes the ppmv etc and applies to TOP bdry
c  (ie record 2.1.1 if TOA is level 1 ... GND at level L, this CORRECTLY applies ppmv to (L)
c ALTZ(L-1),  PZ(L-1),  TZ(L-1),  ATLZ(L),  PZ(L),  TZ(L)
c    ALTZ(L-1), PZ(L-1) and TZ(L-1) are only required for the first layer.
c    LBLRTM assumes that these quantites are equal to the top of the previous
c    layer for L > 1.
								      
      SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases,
     $                          raP,raT,raAlt,iZbnd,raPBnd,
     $                          iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)

      IMPLICIT NONE
     
      INCLUDE '../INCLUDE/kcarta.param'
      include '../INCLUDE/KCARTA_database.param'
      
c input
      CHARACTER*80 caPfName
      REAL rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
c output
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
      REAL raPBnd(2*kProfLayer)  !! do we want user defined pressure level boundaries??
      INTEGER iZbnd              !! do we want user defined pressure level boundaries??

c local var
      INTEGER iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP,IPLEV
      REAL rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,raSumCheck(kMaxGas)
      CHARACTER*80 caStr,caStrX,caStrY
      CHARACTER*30 caStr30      
      CHARACTER*1  c1
      CHARACTER*2  c2
      CHARACTER*3  c3
      CHARACTER*4  c4
      CHARACTER*5  c5
      INTEGER iWriteRTP,iNumGasesBAD,iaBadGasProfile(kMaxGas),iReplaceZeroProf
      INTEGER IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS
      INTEGER XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL
      REAL rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,rDVOUT
      REAL raEmiss(3),raRefl(3),rSecnto
      INTEGER ILNFLG,iForm,iBmax
      REAL raX(KMaxGas),rX,rSatZen
      REAL raZbnd(2*kProfLayer)
      INTEGER iDefault,iAIRS101_or_LBL_levels

      iDefault = +1   !! use AIRS101 levels for the integration
      iAIRS101_or_LBL_levels = +1 !! use AIRS101 levels for the integration
      iAIRS101_or_LBL_levels = -1 !! use LBLRTM  levels for the integration

      IF (iDefault .NE. iAIRS101_or_LBL_levels) THEN
        write(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 when doing integration from levels to layers'
	write(kStdErr,*) 'iDefault, iAIRS101_or_LBL_levels = ',iDefault,iAIRS101_or_LBL_levels
      END IF
      
      rPmin = +1.0e6
      rPmax = -1.0e+6
      iNXsec = -1

      iZbnd = -1   !!! assume we want to use default 101 AIRS levels
      DO iL = 1,kProfLayer+1
        raPbnd(iL) = PLEV_KCARTADATABASE_AIRS(iL)
      END DO
      
      write(kSTdWarn,*) 'Reading in LBLRTM TAPE5 .....'

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

      !! read till you get the $ character, RECORD 1.1
 10   CONTINUE
      READ (iIOUN2,111,ERR=13,END=13) caStr
      IF (caStr(1:1) .NE. '$') GOTO 10

      READ (iIOUN2,111,ERR=13,END=13) caStr  
      CALL read_record_1p2(caStr,
     $          IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)
      IF (ICNTNM .EQ. 6) THEN
        READ (iIOUN2,111,ERR=13,END=13) caStr
        CALL read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)	
      END IF
      IF (IEMIT .EQ. 2) THEN
        CALL DOSTOPMesg('oops cannot handle solar file$')
      END IF

      IF ((IHIRAC .GT. 0) .OR. (IAERSL .GT. 0) .OR. (IEMIT .EQ. 1) .OR. (IATM .EQ. 1) .OR. (ILAS .GT. 0)) THEN
        READ (iIOUN2,111,ERR=13,END=13) caStr
        CALL read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)
      END IF

      !! iotflg comes from solar data
      !! IF ((IEMIT .EQ. 1) .OR. (IEMIT .EQ. 2 .AND.  IOTFLG .EQ. 1)) THEN
      IF ((IEMIT .EQ. 1) .OR. (IEMIT .EQ. 2)) THEN      
        READ (iIOUN2,111,ERR=13,END=13) caStr
        CALL read_record_1p4(caStr,rTSurf,raEmiss,raRefl)
      END IF

 111  FORMAT(A80)
      write(kStdWarn,*) 'TAPE 5 has iAtm = ',iAtm
      IF (IATM .EQ. 0) THEN
        !! need to read in profile
	CALL read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)
        CALL read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck,
     $                         rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta,
     $                         iaG,iaGasUnits,iForm,iNumGases)
        IF (iAIRS101_or_LBL_levels .EQ. +1) THEN
  	  iZbnd = -1            !! use default AIRS 101 levels, they should span past the rPSurf
	ELSE
          iZbnd = iNumLevs  !! use the LBLRTM pressure levels	
        END IF	
        IF (iZbnd .GT. 0) THEN
  	  raPbnd(1) = rPSurf
  	  DO iG = 1,iNumLevs
	    raPbnd(iG) = raP(iG) * 100.0 !! change from mb to N/m2, as Psurf is finally in N/m2
c	    raPbnd(iG) = raP(iG)         !! keep in mb
      	  END DO
	END IF
        IF (iXsect .EQ. 1) THEN
	  CALL read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)
	END IF
      ELSE
        rViewAngle = -9999.0
	rTopHgt    = 705000.0
	CALL read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)
	IF (iBmax .GT. 0) THEN
	  iZbnd = +iBmax
	  CALL read_record_3p3b(iIOUN2,iBmax,raZbnd)
	END IF
        CALL read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf,
     $                          rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          raAlt,raP,raT,raSumCheck,
     $                          iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon,
     $                          iZbnd,raZbnd,raPbnd)
      END IF

      write(kStdWarn,*) 'After exiting the equivalent of nolgas and xsecgas in TAPE5 LBLRTM, found',iNumGases,' gases'
      write(kStdWarn,*) '    these are the gasIDs'
      write(kStdWarn,*) (iaG(iJ),iJ=1,iNumGases)
      write(kStdWarn,*) '    these are the gas units'      
      write(kStdWarn,*) (iaGasUnits(iJ),iJ=1,iNumGases)
      write(kStdWarn,*) ' '
      
      !! now see if data was acutally read in
      iNumGasesBAD = 0
      DO iJ = 1,iNumGases
        IF (raSumCheck(iJ) .LT. 1.0e-20) THEN
          write(kStdWarn,*) ' read in TAPE5, following gas had    zero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
          iNumGasesBAD = iNumGasesBAD + 1
          iaBadGasProfile(iNumGasesBAD) = iJ
        ELSE
          write(kStdWarn,*) ' readin TAPE5, following gas had nonzero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
        END IF
      END DO

      !! need to fix the BAD gases
      iDefault = +1
      iReplaceZeroProf = -1    !! assume user knows why there is a ZERO everywhere gas profile
      iReplaceZeroProf = +1    !! assume user wants to replace ZERO everywhere gas profile with climatology

      IF ((iDefault .NE. iReplaceZeroProf) .AND. (iNumGasesBAD .GT. 0)) THEN
        write(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
        write(kStdErr,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf
        write(kStdWarn,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
        write(kStdWarn,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf 
      END IF
      IF ((iNumGasesBAD .GT. 0) .AND. (iReplaceZeroProf .GT. 0)) THEN
        DO iG = 1,iNumGasesBAD
          CALL substitute_tape5_profile_for_climatology(iaBadGasProfile(iG),iNumLevs,raP,raaG_MR)
          iaGasUnits(iaBadGasProfile(iG)) = 10
          raSumCheck(iaBadGasProfile(iG)) = 0.0
          DO iJ = 1,iNumLevs
            raSumCheck(iaBadGasProfile(iG)) = raSumCheck(iaBadGasProfile(iG)) + 
     $                                        raaG_MR(iJ,iaBadGasProfile(iG))
          END DO
          write(kStdWarn,*) 'reset ZERO gas prof for gasID ',iaBadGasProfile(iG),' ppmv sum = ',raSumCheck(iaBadGasProfile(iG))
        END DO
      END IF

 13   CONTINUE
      CLOSE(iIOUN2) 
      kProfileUnitOpen = -1

      raRTP_TxtInput(1) = rPSurf
      raRTP_TxtInput(2) = rTSurf
      raRTP_TxtInput(3) = rHSurf    !! km
      raRTP_TxtInput(4) = rTophgt   !! km
      raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, downwelling radiation to instr, else upwelling rad to instr

      rPSurf = rPSurf * 100.0       !! change from mb to N/m2
      rHSurf = rHSurf * 1000.0      !! change from km to m

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
      write(kStdWarn,*) 'input file : Top alt (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
      write(kStdWarn,*) ' '

      iWriteRTP = +1
      IF (iWriteRTP .GT. 0) THEN
        CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,
     $                  raP,raT,raaG_MR,raAlt)
      END IF
      
c now change all units first to PPMV (unit 10) and then to Volume Mix Ratio (unit 12)
c actually this is a waste of time as we go from VMR (units 12) to PPMV (units 10) back to VMR (unit 12)
      DO iG = 1,iNumGases
        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR)
        iaGasUnits(iG) = 12
        DO iL = 1,iNumLevs
          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6	  
        END DO
      END DO         

      RETURN
      END

c************************************************************************
c this reads in TAPE6 which has the integrated layer amounts
      SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          iNumLays,rPSurf,rTSurf,rHSurf,iNumGases,raPX,raTX,raLayDensityX,raZX,
     $                          iaG,iaGasUnits,raaG_MRX,rPMin,rPMax,rYear,rLat,rLon)

      IMPLICIT NONE
     
      INCLUDE '../INCLUDE/kcarta.param'

c input
      CHARACTER*80 caPfName
      REAL rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
c output
      INTEGER iNumLays,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rHSurf,rYear,rLat,rLon
      REAL raPX(kProfLayer+1),raTX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas),raLayDensityX(kProfLayer+1)
      REAL raZX(kProfLayer+1)

c local var
      INTEGER iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iNumLevs,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
      REAL raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,r1,r2,r3,r4,r5,r6
      CHARACTER*80 caStrX,caStrY
      CHARACTER*120 caStr
      CHARACTER*30 caStr30
      CHARACTER*1  c1
      CHARACTER*2  c2a,c2b
      INTEGER iWriteRTP,iCountGases,iGCnt,iPass
      INTEGER i1,i2,i3,i4,i5,i6

      rPmin = +1.0e6
      rPmax = -1.0e+6
      iNXsec = -1

      write(kSTdWarn,*) 'Reading in modified LBLRTM TAPE6 = HEAD TAPE5 + MID TAPE6.....'
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

c this should read
c /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
      READ (iIOUN2,111,ERR=222,END=222) caStr
      READ (iIOUN2,111,ERR=222,END=222) caStr

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
      IF (iLBROutBdryHorP .EQ. +1) THEN
        write(kStdWarn,*) 'LBLRTM input indicates start/stop hgts are ',rHSUrf,rTopHgt,' with view angle ',rViewAngle
      ELSEIF (iLBROutBdryHorP .EQ. -1) THEN
        write(kStdWarn,*) 'LBLRTM input indicates start/stop press are ',rHSUrf,rTopHgt,' with view angle ',rViewAngle
      END IF

      READ (iIOUN2,*) (raLBL_Hgts(iJ),iJ=1,iNumLevs)

      raRTP_TxtInput(6) = +raLBL_Hgts(iNumLevs)    !! default, assume hgt in km
      IF  (iLBROutBdryHorP .LT. 0) raRTP_TxtInput(6) = -raLBL_Hgts(iNumLevs)

      READ (iIOUN2,*) iNumLays,iNXsec
      write(kSTdWarn,*) 'Reading TAPE6 info for ',iNumLays,' layers and ',iNumGases,' molgases and ',iNXsec,' xsecgases'

      IF (iNumLays .GT. kProfLayer) THEN
        write(kStdErr,*) 'iNumLays .GT. kProfLayer',iNumLays,kProfLayer
        CALL DoStop
      END IF

      !! we now read every gas 1 .. iNumGases  in edited TAPE6 == molecules/cm2
      DO iG = 1,iNumGases
        iaG(iG) = iG
        iaGasUnits(iG) = 1   !!! assume hardcoded molecules/cm2
      END DO

      READ(iIOUN2,111) caStr      
      IF (caStr(1:6) .NE. '0LAYER') THEN
        write(kStdErr,*) 'oops, expecting string : 0LAYER                          P(MB)       T(K)    ALPHL    ALPHD    ALPHV ...'
        write(kStdErr,*) 'but instead got ',caSTr
        CALL DOStop
      END IF
      !READ(iIOUN2,111) caStr      
      !IF (caStr(1:6) .NE. '     ') THEN
      !  write(kStdErr,*) 'expecting blank string but got ',caStr
      !  CALL DOStop
      !END IF

      DO iL = 1,iNumLays
        READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT
        raZX(iL)   = r1
        raZX(iL+1) = r2
c        print *,iL,iNumLays,raZX(iL),rP,rT
      END DO

      iPass = 1
      iCountGases = 0
      iGCnt = min(7,iNumGases)  !! in first pass, read at most 7 gases
      write(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at pass ',iPass
      READ(iIOUN2,111) caStr      !! blank except for 1 at beginning
      READ(iIOUN2,111) caStr      !! eg LBLRTM    14/07/04  22:24:32
      READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
      READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O 
      DO iL = 1,iNumLays
        READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
        raZX(iL)   = r1
        raZX(iL+1) = r2
c        raPX(iL) = rP * 100.0  !! change from mb to N/m2
        raPX(iL) = rP          !! keep in mb
        raTX(iL) = rT
        IF (rPmax .LE. raPX(iL)) rPmax = raPX(iL)
        IF (rPmin .GT. raPX(iL)) rPmin = raPX(iL)
        IF (iL .EQ. 1) THEN
c          rPSurf = raPX(1)/100 !! because we redo this below, change to mb from N/m2
          rPSurf = raPX(1)      !! keep in mb
        END IF
        DO iG=1,iGCnt
          raaG_MRX(iL,iG+iCountGases) = raX(iG)
        END DO
        raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
      END DO
      READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
      READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM 
      iCountGases = iCountGases + iGCnt

 15   CONTINUE
      IF ((iNumGases .EQ. iCountGases) .AND. (iNXsec .EQ. 0)) GOTO 222     !!! done, just close file!!!!
      IF ((iNumGases .EQ. iCountGases) .AND. (iNXsec .GT. 0)) GOTO 23      !!! done with main gases, need to read xsec gases

      ! need to do this if need to read in more gases
      iPass = iPass + 1
      iGCnt = min(8,iNumGases-iCountGases)      
      write(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at pass ',iPass
      READ(iIOUN2,111) caStr      !! eg 1blank
      READ(iIOUN2,111) caStr      !! LBLRTM    14/07/04  22:24:32
      READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
      READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O 
      DO iL = 1,iNumLays
        READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
        raZX(iL)   = r1
        raZX(iL+1) = r2
c        raPX(iL) = rP * 100.0  !! change from mb to N/m2
        raPX(iL) = rP          !! keep in mb
        raTX(iL) = rT
        IF (rPmax .LE. raPX(iL)) rPmax = raPX(iL)
        IF (rPmin .GT. raPX(iL)) rPmin = raPX(iL)
        IF (iL .EQ. 1) THEN
c          rPSurf = raPX(1)/100 !! because we redo this below
          rPSurf = raPX(1)      !! keep in mb
        END IF
        DO iG=1,iGCnt
          raaG_MRX(iL,iG+iCountGases) = raX(iG)
        END DO
        !!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
      END DO
      iCountGases = iCountGases + iGCnt
      READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
      READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM 
      GOTO 15      

 23   CONTINUE    !! finished reading in LAYER amounts for main gases, now read in LEVELS
      DO i1 = 1,iPass
        READ(iIOUN2,111) caStr         !! blank
        READ(iIOUN2,111) caStr      !! ----
        READ(iIOUN2,111) caStr      !!  MIXING RATIOS BY LAYER
        READ(iIOUN2,111) caStr      !!   P(MB)      T(K)   IPATH         H2O           CO2 
        DO i2 = 1,iNumLays
          READ(iIOUN2,111) caStr      !! mix ratios
        END DO
        READ(iIOUN2,111) caStr      !! blank
        READ(iIOUN2,111) caStr      !! LBLRTM    14/07/04  22:24:32
      END DO

      !!! now ready to read in xsec gases, assume you only need 1 pass
      READ(iIOUN2,111) caStr      !! blank
      READ(iIOUN2,111) caStr      !!  *****  CROSS SECTIONS  *****
      READ(iIOUN2,111) caStr      !! MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
      READ(iIOUN2,111) caStr      !! P(MB)      T(K)   IPATH     F11           F14           F22
      caStr = caStr(56:120)
c      print *,caStr
      CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
      IF (iNXsec .GT. 7) THEN
        write(kStdWarn,*) 'oops assumed at most 7 xsec gases, so only 1 pass ... need to modify code for > 7 gases'
        CALL DOStop
      END IF
      DO iL = 1,iNumLays
        READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iNXsec)
        raZX(iL)   = r1
        raZX(iL+1) = r2
c        raPX(iL) = rP * 100.0  !! change from mb to N/m2
        raPX(iL) = rP          !! keep in mb
        raTX(iL) = rT
        IF (rPmax .LE. raPX(iL)) rPmax = raPX(iL)
        IF (rPmin .GT. raPX(iL)) rPmin = raPX(iL)
        IF (iL .EQ. 1) THEN
c          rPSurf = raPX(1)/100 !! because we redo this below
          rPSurf = raPX(1)     !! keep in mb
        END IF
        DO iG=1,iNXsec
          raaG_MRX(iL,iG+iNumGases) = raX(iG)
        END DO
        !!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
      END DO
      
 222  CONTINUE
      CLOSE(iIOUN2) 
      kProfileUnitOpen = -11
      iNumGases = iNumGases + iNXsec

c      DO iG = 1,iNumGases
c        print *,iG,iaG(iG),iaGasUnits(iG)
c      END DO

 111  FORMAT(A120)

      raRTP_TxtInput(1) = rPSurf
      raRTP_TxtInput(2) = rTSurf
      raRTP_TxtInput(3) = rHSurf    !! km
      raRTP_TxtInput(4) = rTophgt   !! km
      raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, then downwell rad, else upwell radn

      rPSurf = rPSurf * 100.0       !! chnage from mb to N/m2
      rHSurf = rHSurf * 1000.0      !! change from km to m

      write(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
      write(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
      write(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf
     
c make sure pressures are decreasing with index ie layers going higher and higher
      IF (raPX(1) .LT. raPX(2)) THEN
        !!! need to swap!!!
        iMid = ifloor(iNumLays/2.0)
        DO iL = 1,iMid
          iJ = iNumLays-iL+1
          rX = raPX(iJ)
          raPX(iJ) = raPX(iL)
          raPX(iL) = rX

          rX = raTX(iJ)
          raTX(iJ) = raTX(iL)
          raTX(iL) = rX

          DO iG = 1,iNumGases
            raX(iG) = raaG_MRX(iJ,iG)
            raaG_MRX(iJ,iG) = raaG_MRX(iL,iG)
            raaG_MRX(iL,iG) = raX(iG)
          END DO
        END DO
      END IF
      write(kStdWarn,*) 'input file : Highest altitude (lowest press) = ',raPX(iNumLays),' mb at level ',iNumLays
      write(kStdWarn,*) ' '

c      iWriteRTP = +1
c      IF (iWriteRTP .GT. 0) THEN
c        CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,
c     $                  raP,raT,raaG_MR)
c      END IF

c now change all units to MR
c      DO iG = 1,iNumGases
c        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR)
c        DO iL = 1,iNumLevs
c          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
c          iaGasUnits(iG) = 12
c        END DO
c      END DO         

      RETURN
      END

c************************************************************************
c parses LBLRTM xsec names and finds gas ids
      SUBROUTINE XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,iLBLTapeType)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      CHARACTER*80 caStr
      INTEGER iNXsec,iLBLTapeType
c input/output
      INTEGER iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)

c local
      INTEGER iaX(kMaxGas),iG,i1,i2
      CHARACTER*10 caStr10
      CHARACTER*15 caStr15

      IF (iLBLTapeType .EQ. -5) THEN
        DO iG = 1,iNXsec
          i1 = 01 + (iG-1)*10
          i2 = 10 + (iG-1)*10
          caStr10 = caStr(i1:i2)
          CALL mapXsecname_to_XsecID(caStr10,iG,iaX)
        END DO

        DO iG = iNumGases+1,iNumGases+iNXsec
          iaG(iG) = iaX(iG-iNumGases)
          iaGasUnits(iG) = 12   !!! assume hardcoded VMR
        END DO

      ELSEIF (iLBLTapeType .EQ. -6) THEN
        DO iG = 1,iNXsec
          i1 = 01 + (iG-1)*15
          i2 = 10 + (iG-1)*15
          caStr15 = caStr(i1:i2)
          caStr10 = caSTR15(1:10)
          CALL mapXsecname_to_XsecID(caStr10,iG,iaX)
        END DO

        DO iG = iNumGases+1,iNumGases+iNXsec
          iaG(iG) = iaX(iG-iNumGases)
          iaGasUnits(iG) = 1   !!! assume hardcoded molecules/cm2
        END DO
      END IF

      RETURN
      END

c************************************************************************
c this maps xsecnames to xsec IDs
      SUBROUTINE mapXsecname_to_XsecID(caStr10,iG,iaX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      CHARACTER*10 caStr10  !! str we need to compare
      INTEGER iG            !! index
c input/output
      INTEGER iaX(kMaxGas)  !! modify xsec ID according to strcmp

c local 
      INTEGER iI,iFound,iChecked
      CHARACTER*10 caaNamesA(13)
      CHARACTER*10 caaNamesB(13)
      CHARACTER*10 caaNamesC(13)
      CHARACTER*10 caX
      
      IF (caStr10(1:1) .NE. ' ') THEN
        caX = caStr10
      ELSE
        CALL adjustleftstr(caStr10,caX)
      END IF

      caaNamesA(1) = 'CCL3F'
        caaNamesB(1) = 'F11'
        caaNamesC(1) = 'CFC11'
      caaNamesA(2) = 'CCL2F2'
        caaNamesB(2) = 'F12'
        caaNamesC(2) = 'CFC12'
      caaNamesA(3) = 'CCLF3'
        caaNamesB(3) = 'F13'
        caaNamesC(3) = 'CFC13'
      caaNamesA(4) = 'CF4'
        caaNamesB(4) = 'F14'
        caaNamesC(4) = 'CFC14'
      caaNamesA(5) = 'CHCL2F'
        caaNamesB(5) = 'CFC21'
        caaNamesC(5) = 'F21'
      caaNamesA(6) = 'CHC2F2'
        caaNamesB(6) = 'CFC22'
        caaNamesC(6) = 'F22'
      caaNamesA(7) = 'C2CL3F3'
        caaNamesB(7) = 'CFC113'
        caaNamesC(7) = 'F113'
      caaNamesA(8) = 'C2CL2F4'
        caaNamesB(8) = 'CFC114'
        caaNamesC(8) = 'F114'
      caaNamesA(9) = 'C2CLF5'
        caaNamesB(9) = 'CFC115'
        caaNamesC(9) = 'F115'
      caaNamesA(10) = 'CCL4'
        caaNamesB(10) = 'CCL4'
        caaNamesC(10) = ' '
      caaNamesA(11) = 'CLONO2'
        caaNamesB(11) = 'CLONO2'
        caaNamesC(11) = 'CLNO3'
      caaNamesA(12) = 'N2O5'
        caaNamesB(12) = 'N2O5'
        caaNamesC(12) = ' '
      caaNamesA(13) = 'HNO4'
        caaNamesB(13) = 'HNO4'
        caaNamesC(13) = ' '

      iaX(iG) = -1
      iFound = -1
      iChecked = 1
 10   CONTINUE
      IF ((caaNamesA(iChecked) .EQ. caX) .OR. (caaNamesB(iChecked) .EQ. caX) .OR. (caaNamesC(iChecked) .EQ. caX)) THEN
        iFound = iChecked
      END IF

      IF ((iFound .LT. 0) .AND. (iChecked .LT. 13)) THEN
        iChecked = iChecked + 1
        GOTO 10
      END IF

      IF (iFound .LT. 0) THEN
        write(kStdErr,*) 'Reading LBLRTM file; could not find xsec match for ',caX
        CALL DOStop
      ELSE
        iaX(iG) = 50 + iChecked
        write(kStdWarn,*) '  LBLRTM xsec gas = ',caX,' corresponds to gasID ',iaX(iG)
      END IF

      RETURN
      END

c                           -------------------------------------------------
c                             Alias(1)           Alias(2)           Alias(3)           Alias(4)
c                            ----------         ----------         ----------         ----------
c                            CLONO2              CLNO3
c                            HNO4
c                            CHCL2F                                 CFC21              F21
c                            CCL4
c                            CCL3F               CFCL3              CFC11              F11
c                            CCL2F2              CF2CL2             CFC12              F12
c                            C2CL2F4             C2F4CL2            CFC114             F114
c                            C2CL3F3             C2F3CL3            CFC113             F113
c                            N2O5
c                            HNO3
c                            CF4                                    CFC14              F14
c                            CHCLF2              CHF2CL             CFC22              F22
c                            CCLF3                                  CFC13              F13
c                            C2CLF5                                 CFC115             F115

c************************************************************************
c this writes the LBLRTM input so you can cut and paste into Matlab, and save RTP file
c now write the klayers stuff
      SUBROUTINE lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,
     $                  raP,raT,raaG_MR,raAlt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iNumLevs,iNumGases,iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL rPmin,rPmax,rPSurf,rTSurf,rHSurf,rF1,rF2
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas)
      REAL raAlt(2*kProfLayer)
      
c local var
      REAL raX(2*kProfLayer)
      INTEGER iL,iG,iLeftjust_lenstr
      CHARACTER*8 ca8
      CHARACTER*6 ca6
      CHARACTER*2 ca2

      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' '

      write(kStdWarn,*) 'h.vcmin = ',rF1,';'
      write(kStdWarn,*) 'h.vcmax = ',rF2,';'
c      write(kStdWarn,*) 'h.pmin = ',rPmin/100.0,';      %% rPmin in N/m2 --> mb'
c      write(kStdWarn,*) 'h.pmax = ',rPmax/100.0,';      %% rPmax in N/m2 -->  mb'
      write(kStdWarn,*) 'h.pmin = ',rPmin,';      %% in mb'
      write(kStdWarn,*) 'h.pmax = ',rPmax,';      %% im mb'
      write(kStdWarn,*) 'h.ngas = ',iNumGases,';'
      IF (iaGasUnits(1) .GT. 1) THEN
        write(kStdWarn,*) 'h.ptype = 0;'
        write(kStdWarn,*) 'h.pfields = 1;'
      ELSEIF (iaGasUnits(1) .EQ. 1) THEN
        write(kStdWarn,*) 'h.ptype = 1;'
        write(kStdWarn,*) 'h.pfields = 1;'
      END IF
      write(kStdWarn,*) 'h.glist = [',(iaG(iG),iG=1,iNumGases),']'';'
      write(kStdWarn,*) 'h.gunit = [',(iaGasUnits(iG),iG=1,iNumGases),']'';'

      write(kStdWarn,*) 'p.nlevs = ',iNumLevs,';'
      write(kStdWarn,*) 'p.spres = ',rPsurf/100.0,';     %% in mb'
      write(kStdWarn,*) 'p.stemp = ',rTSurf,';'
      write(kStdWarn,*) 'p.salti = ',rHSurf,'; %%%% WOWOWOWOWOW'

      write(kStdWarn,*) 'p.satzen = 0.0;'
      write(kStdWarn,*) 'p.scanang = 0.0;'
      write(kStdWarn,*) 'p.solzen = 130.0;'      
      write(kStdWarn,*) 'p.upwell = 1;'
      write(kStdWarn,*) 'p.zobs = 705000.0;'
      write(kStdWarn,*) 'p.nemis = 2;'
      write(kStdWarn,*) 'p.efreq = [600 3000]'';'
      write(kStdWarn,*) 'p.emis = [1.0 1.0]'';'
      write(kStdWarn,*) 'p.rho = [0.0 0.0]'';'      
      
      ca8 = 'p.plevs'
        !CALL write_stringnice(ca8,raP,0.01,8,iNumLevs,-1)   !! raP in N/m2, convert to mb for rtp 
        CALL write_stringnice(ca8,raP,1.00,8,iNumLevs,-1)	 !! raP in mb, keep in mb for rtp
      ca8 = 'p.ptemp'
        CALL write_stringnice(ca8,raT,1.00,8,iNumLevs,+1)
      ca8 = 'p.palt'
        CALL write_stringnice(ca8,raAlt,1000.00,8,iNumLevs,+1)  !!raAlt in km, convert to m for rtp

      DO iG = 1,iNumGases
        CALL int2str(iaG(iG),ca2)
        ca6 = 'p.gas_'
        ca8 = ca6 // ca2
        DO iL = 1,iNumLevs
          raX(iL) = raaG_MR(iL,iG)
        END DO
        CALL write_stringnice(ca8,raX,1.00,8,iNumLevs,-1)
      END DO

      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' '

      RETURN
      END

c************************************************************************
c this simply writes the strings nicely
      SUBROUTINE write_stringnice(ca8,raX,rScale,iNumItemsPerLine,iNumLevs,ForE)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      CHARACTER*8 ca8
      INTEGER iNumLevs,iNumItemsPerLine,ForE
      REAL raX(2*kProfLayer),rScale

c local var
      INTEGER iL,iY,iFull,iPart,i1,i2,iFloor,iLen
      REAL raY(2*kProfLayer)

      iFull = ifloor(iNumLevs*1.0/iNumItemsPerLine)
      
      DO iL = 1,iNumLevs
        raY(iL) = raX(iL) * rScale
      END DO

      IF (iNumItemsPerLine .NE. 8) iNumItemsPerLine = 8
      
      write(kStdWarn,*) ca8,' = [...'
      DO iL = 1,iFull
        i1 = 1 + (iL-1)* iNumItemsPerLine
        i2 = i1 + iNumItemsPerLine - 1
c        print *,iL,i1,i2
        IF (ForE .EQ. +1) THEN
          write(kStdWarn,1108) (raY(iY),iY=i1,i2)
        ELSEIF (ForE .EQ. -1) THEN
          write(kStdWarn,1308) (raY(iY),iY=i1,i2)
	END IF
      END DO
      IF (i2 .LT. iNumLevs) THEN
        i1 = i2 + 1
        i2 = iNumLevs
c        print *,9999,i1,i2
        iLen = i2-i1+1
	IF (iLen .GT. iNumItemsPerLine) THEN
	  write(kSTdErr,*) 'cannot have more than 8 array members per line'
          CALL DoStop
	END IF
	IF (ForE .EQ. +1) THEN
    	  IF (iLen .EQ. 1) THEN
            write(kStdWarn,1101) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 2) THEN
            write(kStdWarn,1102) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 2) THEN
            write(kStdWarn,1102) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 3) THEN
            write(kStdWarn,1103) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 4) THEN
            write(kStdWarn,1104) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 5) THEN
            write(kStdWarn,1105) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 6) THEN
            write(kStdWarn,1106) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 7) THEN
            write(kStdWarn,1107) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 8) THEN
            write(kStdWarn,1108) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 9) THEN
            write(kStdWarn,1109) (raY(iY),iY=i1,i2)
	  ELSE
	    write(kStdErr,*) 'oops so many numbers in the line????'
	    Call DoStop
	  END IF	
	ELSEIF (ForE .EQ. -1) THEN
    	  IF (iLen .EQ. 1) THEN
            write(kStdWarn,1301) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 2) THEN
            write(kStdWarn,1302) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 2) THEN
            write(kStdWarn,1302) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 3) THEN
            write(kStdWarn,1303) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 4) THEN
            write(kStdWarn,1304) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 5) THEN
            write(kStdWarn,1305) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 6) THEN
            write(kStdWarn,1306) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 7) THEN
            write(kStdWarn,1307) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 8) THEN
            write(kStdWarn,1308) (raY(iY),iY=i1,i2)
	  ELSEIF (iLen .EQ. 9) THEN
            write(kStdWarn,1309) (raY(iY),iY=i1,i2)
	  ELSE
	    write(kStdErr,*) 'oops so many numbers in the line????'
	    Call DoStop
	  END IF
	END IF
      END IF
      write(kStdWarn,*) ']'';'

 1101 FORMAT(1(' ',F12.5),' ...')
 1102 FORMAT(2(' ',F12.5),' ...')
 1103 FORMAT(3(' ',F12.5),' ...')
 1104 FORMAT(4(' ',F12.5),' ...')
 1105 FORMAT(5(' ',F12.5),' ...')
 1106 FORMAT(6(' ',F12.5),' ...')
 1107 FORMAT(7(' ',F12.5),' ...')
 1108 FORMAT(8(' ',F12.5),' ...')
 1109 FORMAT(9(' ',F12.5),' ...')

 1301 FORMAT(1(' ',E12.5),' ...')
 1302 FORMAT(2(' ',E12.5),' ...')
 1303 FORMAT(3(' ',E12.5),' ...')
 1304 FORMAT(4(' ',E12.5),' ...')
 1305 FORMAT(5(' ',E12.5),' ...')
 1306 FORMAT(6(' ',E12.5),' ...')
 1307 FORMAT(7(' ',E12.5),' ...')
 1308 FORMAT(8(' ',E12.5),' ...')
 1309 FORMAT(9(' ',E12.5),' ...')
       
      RETURN
      END

c************************************************************************
c this subroutine just takes the LBLRTM mol/cm2 TAPE6 profile and does some simple calcs
      SUBROUTINE DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon,
     $               PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,
     $               raPX,raTX,raLayDensityX,raaG_MRX,
     $               raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL rHSurf,rPSurf,rTSurf,PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      REAL DATABASELEVHEIGHTS(kMaxLayer+1),rfracBot,rLat,rLon
      INTEGER iNumLays,iNumGases,iaG(kMaxGas)
      REAL raPX(kProfLayer+1),raTX(kProfLayer+1),raLayDensityX(kProfLayer+1),raaG_MRX(kProfLayer+1,kMaxGas)
c output
      REAL raTout(kProfLayer),raAmountOut(kProfLayer),raZout(kProfLayer+1),raPout(kProfLayer)
      REAL raaQout(kProfLayer,kMaxGas),raaPartPressOut(kProfLayer,kMaxGas)
      INTEGER iLowestLev

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

      DO iL = 1,kProfLayer
        DO iG = 1,iNumGases
          raaQout(iL,iG) = 0.0
        END DO
      END DO

      IF ((kPlanet .EQ. 3) .AND. (iaG(1) .NE. 1)) THEN
        write (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
        Call DoStop
      END IF

      iLowestLev = kProfLayer - iNumLays + 1

      rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/
     $         (PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
      ! just readjust the layer info read in, and tack on additional gas profiles from refgas
      write(kSTdWarn,*) 'iLowestLev = ',iLowestLev,' rfracBot = ',rFracBot
      write(kStdWarn,*) 'PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev) = ',
     $   PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev)

      DO iL = iLowestLev,kProfLayer
        raZout(iL) = DATABASELEVHEIGHTS(iL)
        raPout(iL) = raPX(iL-iLowestLev+1)
        raTout(iL) = raTX(iL-iLowestLev+1)
        raAmountOut(iL) = raLayDensityX(iL-iLowestLev+1) 
        IF (iL .EQ. iLowestLev) raAmountout(iL) = raAmountout(iL)/rFracBot
        DO iG = 1,iNumGases
          raaQout(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)
          IF (iL .EQ. iLowestLev) raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
          raaPartPressOut(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)/raLayDensityX(iL-iLowestLev+1) * raPX(iL-iLowestLev+1)
        END DO
      END DO

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
        write(kStdWarn,*) 'NOT Checking Relative Humidity'
        write(kStdWarn,*) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      Q         Q1         Q2'
        write(kStdWarn,*) '-------------------------------------------------------------------------------------------------------'
        DO iL = iLowestLev,kProfLayer
          z    = (raZout(iL) + raZout(iL+1))/2/1000
          zWoo = (raZout(iL+1) - raZout(iL))/1000
          rP   = raPout(iL)            !! N/m2
          rT   = raTout(iL)
          rPP  = raaPartPressOut(iL,1) !! N/m2
          rSVP = 0.0
          rRH  = 0.0
          IF (rRH .LE. 100.0) THEN
            write(kStdWarn,111) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),
     $ (raaQout(iL,iG),iG=1,6)
          ELSE
            write(kStdWarn,112) iL-iLowestLev+1,iL,rP/100.0,z,zWoo,rT,rPP/100.0,rSVP/100,rRH,raAmountOut(iL),
     $ (raaQout(iL,iG),iG=1,6),' ***** '
          END IF
        END DO
      END IF
 111  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '))
 112  FORMAT(2(I3,' '),1(F9.4,' '),6(F8.4,' '),7(E9.3,' '),A7)

c now find pressure output corresponding to HGT output from LBLRTM
      IF (kRTP .EQ. -20) THEN
        IF (raRTP_TxtInput(6) .GT. 0) THEN
          !!! input level boundaries in km, change to mb
          DO iL = iLowestLev,kProfLayer
            raJunk(iL-iLowestLev+1) = raZout(iL)          
            raJunk2(iL-iLowestLev+1) = raPout(iL)          
          END DO
          CALL r_sort_linear(raJunk,raJunk2,kProfLayer-iLowestLev+1,raRTP_TxtInput(4)*1000,zWoo,1)
          write(kStdWarn,*)'LBLRTM output height of ',raRTP_TxtInput(4),' km corresponds to ',zWoo,' N/m2'
c          raRTP_TxtInput(6) = zWoo/100.0  !! mb
          raRTP_TxtInput(6) = zWoo        !! keep in mb
        ELSEIF (raRTP_TxtInput(6) .LT. 0) THEN
          !!! input level boundaries in mb     
c          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
c          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) 
          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' mb'
        END IF
      END IF

      RETURN
      END

c************************************************************************
c if TAPE5 has zeros everywhere for a certain gas, this reads in a US Std profile and stuffs it in
      SUBROUTINE substitute_tape5_profile_for_climatology(iG,iNumLevs,raP,raaG_MR)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c input
      INTEGER iG             ! which gasID to put in profile
      INTEGER iNumLevs       ! how many levs in profile
      REAL raP(2*kProfLayer) ! pressure levels
c input/output
      REAL raaG_MR(2*kProfLayer,kMaxGas)

c local vars
      REAL raX(2*kProfLayer)
      INTEGER iL

      CALL ReadRefProf_Levels2(iG,raP,iNumLevs,raX)

      DO iL = 1,iNumLevs
        raaG_MR(iL,iG) = raX(iL)
      END DO

      RETURN
      END

c************************************************************************
c this reads record 1.2 of a LBLRTM TAPE5
      SUBROUTINE read_record_1p2(caStr,
     $              IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)

      IMPLICIT NONE
      
c input
      CHARACTER*80 caStr
c output
      INTEGER IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS

c local      
      CHARACTER*1  c1
      CHARACTER*2  c2
      CHARACTER*3  c3
      CHARACTER*4  c4
      CHARACTER*5  c5

c RECORD 1.2
c      IHIRAC, ILBLF4, ICNTNM, IAERSL,  IEMIT,  ISCAN, IFILTR, IPLOT, ITEST,  IATM,  IMRG,  ILAS,   IOD, IXSECT,  MPTS,  NPTS
c           5,     10,     15,     20,     25,     30,     35,    40,    45,    50, 54-55,    60,    65,     70, 72-75, 77-80
c       4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4
c this should read 
c /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
      !! reads HI=1 F4=1 CN=5 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0
      !! reads HI=1 F4=1 CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0

      c1 = caStr(5:5)
      read(c1,'(I1)') IHIRAC
      c1 = caStr(10:10)
      read(c1,'(I1)') ILBLF4
      c1 = caStr(15:15)
      read(c1,'(I1)') ICNTNM
      c1 = caStr(20:20)
      read(c1,'(I1)') IAERSL
      c1 = caStr(25:25)
      read(c1,'(I1)') IEMIT
      c1 = caStr(30:30)
      read(c1,'(I1)') ISCAN
      c1 = caStr(35:35)
      read(c1,'(I1)') IFILTR
      c1 = caStr(40:40)
      read(c1,'(I1)') IPLOT
      c1 = caStr(45:45)
      read(c1,'(I1)') ITEST
      c1 = caStr(50:50)
      read(c1,'(I1)') IATM
c      c2 = caStr(54:55)
c      read(c1,'(I2)') IMRG      
c      read(c2,2) IMRG
      c1 = caStr(55:55)
      read(c1,'(I1)') IMRG
      c1 = caStr(60:60)
      read(c1,'(I1)') ILAS
      c1 = caStr(65:65)
      read(c1,'(I1)') IOD
      c1 = caStr(70:70)
      read(c1,'(I1)') IXSECT
      c4 = caStr(72:75)
      read(c4,'(I4)') MPTS
      c4 = caStr(77:80)
      read(c4,'(I4)') NPTS

 2    FORMAT(I2)
 4    FORMAT(I4)
 
c      print *,'record 1.2 : ',IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS

      RETURN
      END

c************************************************************************
c this reads record 1.2a of a LBLRTM TAPE5
      SUBROUTINE read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)

      IMPLICIT NONE

      !!   XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL in free format
      !!     XSELF  H2O self broadened continuum absorption multiplicative factor
      !!     XFRGN  H2O foreign broadened continuum absorption multiplicative factor
      !!     XCO2C  CO2 continuum absorption multiplicative factor
      !!     XO3CN  O3 continuum absorption multiplicative factor
      !!     XO2CN  O2 continuum absorption multiplicative factor
      !!     XN2CN  N2 continuum absorption multiplicative factor
      !!      XRAYL Rayleigh extinction multiplicative factor

c input
      CHARACTER*80 caStr
c output
      INTEGER XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

      READ (caStr,*) XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL
      
      RETURN
      END
      
c************************************************************************
c this reads record 1.3 of a LBLRTM TAPE5
      SUBROUTINE read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)

      IMPLICIT NONE
      
c input
      CHARACTER*80 caStr
c output
      REAL rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,rDVOUT
      INTEGER ILNFLG

c local      
      CHARACTER*1  c1
      CHARACTER*10  c10

c RECORD 1.3    (required if IHIRAC > 0; IAERSL > 0; IEMIT = 1; IATM = 1; or ILAS > 0; otherwise omit)
c             V1,     V2,   SAMPLE,   DVSET,  ALFAL0,   AVMASS,   DPTMIN,   DPTFAC,   ILNFLG,     DVOUT
c           1-10,  11-20,    21-30,   31-40,   41-50,    51-60,    61-70,    71-80,     85,      90-100
c          E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3

 10   FORMAT(E10.3)

      rSAMPLE = 4.0
      rALFAL0 = 0.04
      rAVMASS = 36
      rDPTMIN = 0.0002
      rDPTFAC = 0.001
      ILNFLG = 0
      rDVOUT = 0.0
      
      c10 = caStr(01:10)
      read(c10,10) rF1
      c10 = caStr(11:20)
      read(c10,10) rF2

      c10 = caStr(21:30)
      read(c10,10) rSAMPLE
      c10 = caStr(31:40)
      read(c10,10) rDVSET
      c10 = caStr(41:50)
      read(c10,10) rALFAL0
      c10 = caStr(51:60)
      read(c10,10) rAVMASS
      c10 = caStr(61:70)
      read(c10,10) rDPTMIN
      c10 = caStr(71:80)
      read(c10,10) rDPTFAC
c      c1 = caStr(85:85)
c      read(c1,'(I1)') ILFLG
c      c1 = caStr(91:100)
c      read(c10,10) rDVOUT

c      print *,'record 1.3 : ',rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT

      RETURN
      END

c************************************************************************
c this reads record 1.4 of a LBLRTM TAPE5
      SUBROUTINE read_record_1p4(caStr,rTSurf,raEmiss,raRefl)

      IMPLICIT NONE
      
c RECORD 1.4    (required if IEMIT = 1, or both IEMIT=2 and IOTFLG=2; otherwise omit)
c         TBOUND, SREMIS(1), SREMIS(2), SREMIS(3), SRREFL(1), SRREFL(2), SRREFL(3)
c           1-10,     11-20,     21-30,     31-40,     41-50,     51-60,     61-70
c          E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3	  

c input
      CHARACTER*80 caStr
c output
      REAL rTSurf,raEmiss(3),raRefl(3)

c local      
      CHARACTER*10  c10
      INTEGER iI

 10   FORMAT(E10.3)
 
      c10 = caStr(01:10)
      read(c10,10) rTSurf
      DO iI = 1,3
        raEmiss(iI) = -999.0
        c10 = caStr(iI*10+1:iI*10+10)
        read(c10,10) raEmiss(iI)
      END DO
      IF (raEmiss(1) .LT. 0.0) THEN
        !! expecting emiss file
	raEmiss(2) = -1.0
	raEmiss(3) = -1.0
      END IF
      
      DO iI = 4,6
        raRefl(iI-3) = -999.0      
        c10 = caStr(iI*10+1:iI*10+10)
        read(c10,10) raRefl(iI-3)
      END DO
      IF (raRefl(1) .LT. 0.0) THEN
        !! expecting emiss file
	raRefl(2) = -1.0
	raRefl(3) = -1.0
      END IF

c      print *,(raEmiss(iI),iI=1,3)
c      print *,(raRefl(iI),iI=1,3)
c      print *,'record 1.4 : tsurf = ',rTsurf

      RETURN
      END
      
c************************************************************************      
c this reads record 2.1 of a LBLRTM TAPE5
      SUBROUTINE read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)

      IMPLICIT NONE
      
c         IFORM, NLAYRS, NMOL, SECNTO,        ZH1,       ZH2,    ZANGLE
c            2     3-5,   6-10,  11-20,      41-48,     53-60,     66-73
c          1X,I1    I3,    I5,   F10.2,  20X, F8.2,  4X, F8.2,  5X, F8.3

c              IFORM      (0,1) column amount format flag
c      		         = 0  read PAVE, WKL(M,L), WBROADL(L) in F10.4, E10.3, E10.3 formats (default)
c                        = 1  read PAVE, WKL(M,L), WBROADL(L) in E15.7 format
c
c              NLAYRS      number of layers (maximum of 200)
c	      
c                NMOL      value of highest molecule number used (default = 7; maximum of 35)
c		                                 See Table I for molecule numbers.
c              SECNTO      user entered scale factor for the column amount for the layers defined by NLAYRS
c	                                       if positive, looking up
c	                                        if negative, looking down
c          		                          normal value = 1.0
c                 ZH1      observer altitude
c                 ZH2      end point altitude
c              ZANGLE      zenith angle at H1 (degrees)

c input
      INTEGER iIOUN2
c output
      INTEGER iForm,iNumLevs,iNumGases
      REAL rSecnto,rTopHgt,rHSurf,rViewAngle

c local
      CHARACTER*80 caStr      
      CHARACTER*10 c10
      CHARACTER*1  c1
      CHARACTER*3  c3
      CHARACTER*5  c5
      CHARACTER*8  c8      

      REAL rJunk
      INTEGER iI

 8    FORMAT(E8.3)
 10   FORMAT(E10.3)
 111  FORMAT(A80)
 
      READ (iIOUN2,111,ERR=13,END=13) caStr

      c1 = caStr(2:2)
      read(c1,*) iFORM
      c3 = caStr(3:5)
      read(c3,*) iNumLevs
      !! this is really number of layers, so need to increment by 1
      !! iNumLevs = iNumLevs + 1
      c5 = caStr(6:10)
      read(c5,*) iNumGases
      c10 = caStr(11:20)
      read(c10,10) rSecnto   !! scale factor for layer ODS : positive if looking up, -ve if looking down
      
      c8 = caStr(41:48)
      read(c8,8) rTopHgt     !! observer altitude
      c8 = caStr(53:60)
      read(c8,8) rHSurf      !! end point altitude
      c8 = caStr(66:73)
      read(c8,8) rViewAngle  !! zenith angle at observer altitude (rTopHgt) so this is scanang

      IF (rHSurf .GT. rTopHgt) THEN
        rJunk = rTopHgt
	rTopHgt = rHSurf
	rHSurf = rJunk
      END IF
      
c      print *,'record 2.1 : ',iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle
 13   CONTINUE
 
      RETURN
      END
      
c************************************************************************      
c this reads record 2.2 of a LBLRTM TAPE5
      SUBROUTINE read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck,
     $                         rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta,      
     $                         iaG,iaGasUnits,iForm,iNumGases)

      IMPLICIT NONE     
      INCLUDE '../INCLUDE/kcarta.param'
   
c input
      INTEGER iNumLevs,iIOUN2,iForm,iNumGases
      REAL rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta
c output
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),rPSurf,rHSurf,raSumCheck(kMaxGas)
      REAL rTSurf,rPmin,rPmax,raAlt(2*kProfLayer)
      INTEGER iaG(kMaxGas),iaGasUnits(kMaxGas)
      
c local
      INTEGER iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2
      CHARACTER*120 caStrX,caStrY
      CHARACTER*30  caStr30
      CHARACTER*15  c15      
      CHARACTER*10  c10
      CHARACTER*8  c8
      CHARACTER*7  c7
      CHARACTER*3  c3
      CHARACTER*2  c2
      CHARACTER*1  c1
      CHARACTER*80 ca80
      REAL rH,rP,rT,raX(kMaxGas),rJunk,rHSurfJunk

 120  FORMAT(A120)
 15   FORMAT(E15.7) 
 10   FORMAT(F10.4)
 8    FORMAT(F8.3)
 7    FORMAT(F7.2) 
 
      DO iG = 1,iNumGases
        iaG(iG) = iG
        iaGasUnits(iG) = 12   !!! assume hardcoded VMR
      END DO

      DO iJ = 1,iNumGases
        raSumCheck(iJ) = 0.0
      END DO
        
      DO iL = 1,iNumLevs              !!! >>>>>>>>>>>>>>>>>>>> start reading the MOLGAS profiles
        READ (iIOUN2,120) caStrY

        !! first we want to read Pave and Tave
        IF (iForm .EQ. 0) THEN
          c10 = caStrY(01:10)
          read(c10,10) rP
          c10 = caStrY(11:20)
          read(c10,10) rT
          c10 = caStrY(21:30)
          read(c10,10) rJunk
          c3 = caStrY(31:33)
          c2 = caStrY(34:35)
          ca80 = caStrY(37:120)
	ELSEIF (iForm .EQ. 1) THEN
          c15 = caStrY(01:15)
          read(c15,15) rP
          c10 = caStrY(16:25)
          read(c10,10) rT
          c10 = caStrY(26:35)
          read(c10,10) rJunk
          c3 = caStrY(36:38)
          c2 = caStrY(39:40)
          ca80 = caStrY(41:120)
	END IF

	IF (iL. EQ. 1) THEN
          !! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
	  !!                      which are basically
          !!          SurfAlt	Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
          READ (ca80,*) rHSurfJunk,rPSurf,rTSurf,raAlt(iL+1),raP(iL+1),raT(iL+1)   !! p is in mb
	  raAlt(iL) = rHSurfJunk
	  raP(iL)   = rPSurf
	  raT(iL)   = rTSurf
          IF (rPmax .LE. raP(iL)) rPmax = raP(iL)
          IF (rPmin .GT. raP(iL)) rPmin = raP(iL)	    
          IF (rPmax .LE. rPSurf) rPmax = rPSurf
          IF (rPmin .GT. rPSurf) rPmin = rPSurf
	  
	  IF (abs(rHSurfJunk-rHSurf) .GT. 0.001) THEN
	    write(kStdWarn,*) 'resetting rHSurf from first levl info in TAPE5 ',rHSurfJunk,rHSurf
	    rHSurf = rHSurfJunk
	  END IF

          IF (rPmax .LE. raP(iL+1)) rPmax = raP(iL+1)
          IF (rPmin .GT. raP(iL+1)) rPmin = raP(iL+1)	    

          IF (raP(iL+1) .GT. raP(iL)) THEN
	    write(kStdErr,*) 'huh reading LBLRTM TAPE5  iL,raP(iL),raP(iL+1) = ',iL,raP(iL),raP(iL+1)
	    CALL DoStop
	  END IF
	  
c	  print *,iL,raAlt(iL),raP(iL),raT(iL),rHSurfJunk,rPSurf,rTSurf

	ELSE
          !! now read next "BLANK"  and A(z) P(z) and T(z)
	  !!                         which are basically
          !!                         A(z) Pz)  and T(z) for next level, ad continuum
          READ (ca80,*) raAlt(iL+1),raP(iL+1),raT(iL+1)    !! p in mb
          IF (rPmax .LE. raP(iL+1)) rPmax = raP(iL+1)
          IF (rPmin .GT. raP(iL+1)) rPmin = raP(iL+1)
c	  print *,iL,raAlt(iL+1),raP(iL+1),raT(iL+1)	  
	END IF
	
        !! now read MixRatio(z)
        IF (iNumGases .LE. 7) THEN
          READ (iIOUN2,120) caStrY	  
	  read (caStrY,*) (raX(iJ),iJ=1,iNumGases)
	ELSEIF (iNumGases .GT. 7) THEN
	  iGasCntLoop = 0
 100      CONTINUE
          iX1 = iGasCntLoop*7 + 1
	  iX2 = iX1 + 7
          iX2 = min(iX2,iNumGases)
          READ (iIOUN2,120) caStrY	  
	  read (caStrY,*) (raX(iJ),iJ=iX1,iX2)
	  IF (iX2 .LT. iNumGases) THEN
	    iGasCntLoop = iGasCntLoop + 1
	    GOTO 100
	  END IF
        END IF

        DO iJ = 1,iNumGases
          raaG_MR(iL+1,iJ) = raX(iJ)
          raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
        END DO
	IF (iL .EQ. 1) THEN
          DO iJ = 1,iNumGases
            raaG_MR(iL,iJ) = raX(iJ)
            raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
          END DO
        END IF
	
        IF (iL .EQ. 1) THEN
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
	
      END DO       !!!! <<<<<<<<<<<<<<<<<<<< done reading the MOLGAS profiles, loop over levs

      iNumLevs = iNumLevs+1   !!! since we added on surface level
      
      kLBLRTM_toa = raP(iNumLevs)
      write(kStdWarn,*) 'kLBLRTM_toa = ',kLBLRTM_toa,' mb (when used with flux calcs)'
      write(kStdWarn,*) ' '
      
      RETURN
      END

c************************************************************************      
c read the xsect
c this reads record 2.2 of LBLRTM TAPE5
      SUBROUTINE read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)

      IMPLICIT NONE     
      INCLUDE '../INCLUDE/kcarta.param'
   
c input
      INTEGER iIOUN2,iNumLevs,iNumGases,iNXsec
      REAL raP(2*kProfLayer),raT(2*kProfLayer)
c output
      REAL raaG_MR(2*kProfLayer,kMaxGas),raSumCheck(kMaxGas)
      INTEGER iaG(kMaxGas),iaGasUnits(kMaxGas)
      
c local
      INTEGER iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2,iXSBIN,iFormX,iXmol,iNumLevsXsec
      CHARACTER*120 caStr,caStrX,caStrY
      CHARACTER*30 caStr30
      CHARACTER*15  c15      
      CHARACTER*10  c10
      CHARACTER*8  c8
      CHARACTER*7  c7
      CHARACTER*3  c3
      CHARACTER*2  c2
      CHARACTER*1  c1
      CHARACTER*80 ca80
      REAL raPX(2*kProfLayer),raTX(2*kProfLayer),raAltX(2*kProfLayer),raX(kMaxGas),rP,rT,rJunk
      REAL rTSurfX,rHSUrfX,rPSurfX

 111  FORMAT(A80)
 120  FORMAT(A120)
 15   FORMAT(E15.7) 
 10   FORMAT(F10.4)
 8    FORMAT(F8.3)
 7    FORMAT(F7.2) 

      !! now see if there are xsec gases
      READ (iIOUN2,111,ERR=13,END=13) caStr
      READ(caStr,*) iNXsec,iXSBIN
      IF (iNXsec .GT. 0) THEN    !!!! >>>>>>>>>>>>>>> start reading XSCGAS profiles

        DO iG = iNumGases+1,iNumGases+iNXsec
          iaG(iG) = iG
          iaGasUnits(iG) = 12   !!! assume hardcoded VMR
        END DO

        DO iJ = iNumGases+1,iNumGases+iNXsec
          raSumCheck(iJ) = 0.0
        END DO

        READ (iIOUN2,111,ERR=13,END=13) caStr     !!!! xsec names 
        CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
	
        READ (iIOUN2,*) iFormX,iNumLevsXsec,iXmol
        IF (iNumLevsXsec .GT. 2*kProfLayer) THEN
          write(kStdErr,*) 'iNumLevsXsec .GT. 2*kProfLayer',iNumLevsXsec,kProfLayer
          CALL DoStop
        END IF
        IF (iNumLevsXsec .NE. (iNumLevs-1)) THEN
          write(kStdErr,*) 'iNumLevsXsec .NE. iNumLevs',iNumLevsXsec,iNumLevs-1
          CALL DoStop
        END IF

        !! now ready to read in the profiles, same as in molgas above!
        DO iL = 1,iNumLevsXsec
          READ (iIOUN2,111) caStrY

          !! first we want to read Pave and Tave
          IF (iFormX .EQ. 0) THEN
            c10 = caStrY(01:10)
            read(c10,10) rP
            c10 = caStrY(11:20)
            read(c10,10) rT
            c10 = caStrY(21:30)
            read(c10,10) rJunk
            c3 = caStrY(31:33)
            c2 = caStrY(34:35)
           ca80 = caStrY(37:120)
  	  ELSEIF (iFormX .EQ. 1) THEN
            c15 = caStrY(01:15)
            read(c15,15) rP
            c10 = caStrY(16:25)
            read(c10,10) rT
            c10 = caStrY(26:35)
            read(c10,10) rJunk
            c3 = caStrY(36:38)
            c2 = caStrY(39:40)
            ca80 = caStrY(41:120)
          END IF

          IF (iL. EQ. 1) THEN
            !! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
  	    !!                      which are basically
            !!          SurfAlt	Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
            READ (ca80,*) rHSurfX,rPSurfX,rTSurfX,raAltX(iL+1),raPX(iL+1),raTX(iL+1)   !! p is in mb	  
	  ELSE
            READ (ca80,*) raAltX(iL+1),raPX(iL+1),raTX(iL+1)
          END IF

          IF (iNXsec .LE. 7) THEN
            READ (iIOUN2,120) caStrY	  
  	    read (caStrY,*) (raX(iJ),iJ=1,iNXsec)
	  ELSEIF (iNXsec .GT. 7) THEN
	    iGasCntLoop = 0
 100        CONTINUE
            iX1 = iGasCntLoop*7 + 1
	    iX2 = iX1 + 7
            iX2 = min(iX2,iNXsec)
            READ (iIOUN2,120) caStrY	  
	    read (caStrY,*) (raX(iJ),iJ=iX1,iX2)
	    IF (iX2 .LT. iNXsec) THEN
	      iGasCntLoop = iGasCntLoop + 1
	      GOTO 100
	    END IF
          END IF

          IF (iL .eq. 1) THEN
            DO iJ = 1,iNXsec
              raaG_MR(iL,iJ+iNumGases) = raX(iJ)
              raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
            END DO
          END IF 
           
          DO iJ = 1,iNXsec
            raaG_MR(iL+1,iJ+iNumGases) = raX(iJ)
            raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
          END DO

	END DO   !!! loop over levels
      END IF     !!! end reading xsc gas profiles

      iNumGases = iNumGases + iNXsec

 13   CONTINUE
 
      RETURN
      END
c************************************************************************
c this reads record 3.1 of LBLRTM TAPE5
      SUBROUTINE read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)

      IMPLICIT NONE     
      INCLUDE '../INCLUDE/kcarta.param'
	 
c input
      INTEGER iIOUN2
c output
      INTEGER iNumGases,iBmax
      REAL rTopHgt,rSatZen,rHSurf

c local
      INTEGER iJunk,iI,iY
      CHARACTER*80 caStr
      CHARACTER*1 c1
      CHARACTER*2 c2
      CHARACTER*5 c5      
      CHARACTER*10 c10      
      REAL rCO2MIX,rJUNK

c 3.1
c      MODEL,  ITYPE, IBMAX,  NOZERO,  NOPRNT,  NMOL, IPUNCH, IFXTYP,   MUNITS,    RE, HSPACE,  VBAR, CO2MX
c          5,     10,    15,      20,      25,    30,     35,  36-37,    39-40, 41-50,  51-60, 61-70, 71-80
c         I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3, F10.3
      READ(iIOUN2,111) caStr
      
      c5 = caStr(1:5)
      read(c5,*) iJunk
      IF (iJunk .NE. 0) THEN
        write(kStdErr,*) 'in record 3.1 LBLRTM specifies one of its own internal profile OOPS'
        CALL dostop
      END IF

      c5 = caStr(6:10)
      read(c5,*) iJunk
      IF (iJunk .NE. 3) THEN
        write(kStdWarn,*) 'in record 3.1 LBLRTM specifies iType = ',iJunk,' OOPS only want 3 (gnd to space), resetting'
        iJunk = 3
      END IF

      c5 = caStr(11:15)
      read(c5,*) iBmax
      IF (iBmax .EQ. 0) THEN
        write(kStdWarn,*) 'in record 3.1 LBLRTM will use its own boundaries'
      ELSE
        write(kStdWarn,*) 'in record 3.1 LBLRTM will read in user specified boundaries'      
      END IF

      c5 = caStr(16:20)
      read(c5,*) iJunk
      c5 = caStr(21:25)
      read(c5,*) iJunk
      c5 = caStr(26:30)
      read(c5,*) iNumGases
      c5 = caStr(31:35)
      read(c5,*) iJunk

      rCO2MIX = 330.0
      c10 = caStr(71:80)
      iY = -1
      iI = 1
      DO iI = 1,10
        IF (c10(iI:iI) .NE. ' ') iY = +1
      END DO
      IF  (iY .GT. 0) read(c10,*) rCO2mix

c 3.2
c         H1,    H2,   ANGLE,   RANGE,   BETA,   LEN,     HOBS
c       1-10, 11-20,   21-30,   31-40,  41-50, 51-55,    61-70
c      F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3
      
      READ(iIOUN2,111) caStr
      c10 = caStr(1:10)
      read (c10,*) rHSurf
      c10 = caStr(11:20)
      read (c10,*) rTopHgt
      c10 = caStr(21:30)
      read (c10,*) rSatZen

      IF (rHSurf .GT. rTopHgt) THEN
        rJunk = rTopHgt
	rTopHgt = rHSurf
	rHSurf = rJunk
      END IF

 111  FORMAT(A80)
 
      RETURN
      END
c************************************************************************
c this reads in record 3.3b of LBLRTM TAPE5
c altitudes of LBLRTM layer boundaries
      SUBROUTINE read_record_3p3b(iIOUN2,iBmax,raZbnd)

      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c input
      INTEGER iBmax,iIOUN2
c output
      REAL raZbnd(2*kProfLayer)

c local
      INTEGER iI,iCnt,i1,i2,iMax

      iCnt = 1
 20   CONTINUE
      i1 = (iCnt-1)*8 + 1 
      i2 = iCnt*8
      IF (i2 .GT. iBmax) i2 = iBmax
      READ(iIOUN2,*) (raZbnd(iI),iI=i1,i2)
      IF (i2 .LT. iBmax) THEN
        iCnt = iCnt + 1
	GOTO 20
      END IF
 
      RETURN
      END

c************************************************************************      
c this reads in record 3.4 of LBLRTM TAPE5
c user defined profile
      SUBROUTINE read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf,
     $                          rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,
     $                          raAlt,raP,raT,raSumCheck,
     $                          iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon,
     $                          iZbnd,raZbnd,raPbnd)
     
      IMPLICIT NONE
      include '../INCLUDE/kcarta.param'

c input
      INTEGER iNumGases,iIOUN2
c output
      INTEGER iNumLevs
      REAL raSumCheck(kMaxGas),rPSurf,rHSurf,rTSurf
      REAL rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
      INTEGER iaG(kMaxGas),iaGasUnits(kMaxGas)
      REAL rPmin,rPmax,rYear,rLat,rLon
      REAL raP(2*kProfLayer),raT(2*kProfLayer),raaG_MR(2*kProfLayer,kMaxGas),raAlt(2*kProfLayer)
      REAL raZbnd(2*kProfLayer),raPbnd(2*kProfLayer)
      INTEGER iZbnd   !!! are we using default 101 levels, or LBLRTM defined?     
     
c local var
      CHARACTER*5  c5
      INTEGER iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
      REAL raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH
      CHARACTER*80 caStr,caStrX,caStrY
      CHARACTER*30 caStr30
      CHARACTER*1  c1

      DO iJ = 1,iNumGases
        raSumCheck(iJ) = 0.0
      END DO

      rPmin = +1.0e6
      rPmax = -1.0e+6
      iNXsec = -1

c record 3.4
 111  FORMAT(A80)
      read(iIOUN2,111) caStr
      c5 = caStr(1:5)
      read(c5,*) iNumLevs

c record 3.5
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
	IF (iL .EQ. 1) THEN
	  rPSurf = rP
	  rHSurf = rH
	  IF (abs(rTSurf-rT) .GE. 0.001) THEN
	    !! hmm looks like info in first item of Record 3.5 is inconsistent with info in Record 1.4
	    write(kStdWarn,*) 'rT first item of Record 3.5 is inconsistent with TBound info in Record 1.4'
	    write(kSTdWarn,*) rTSurf,rT
	    write(kStdWarn,*) 'artificially set rPSurf to be a little higher than lowest "p" entry'
            rPSurf = rPSurf-0.125
	  END IF
	  !rTSurf = rT
	END IF
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

c        raP(iL) = rP * 100.0  !! change from mb to N/m2
        raP(iL) = rP          !! keep in mb
        raT(iL) = rT
	raAlt(iL) = rH
        IF (rPmax .LE. raP(iL)) rPmax = raP(iL)
        IF (rPmin .GT. raP(iL)) rPmin = raP(iL)
        READ (iIOUN2,*) (raX(iG),iG=1,iNumGases)	
        DO iJ = 1,iNumGases
          raaG_MR(iL,iJ) = raX(iJ)
          raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)	  
        END DO

        IF (iL .EQ. 1) THEN
	   !! rPSurf already set a few lines above, and need to be careful since we want bdry levels set
	   !! so do NOT play with it here
c          rPSurf = raP(1)/100 !! because we redo this below
c          rPSurf = raP(1)     !! keep in mb        
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
            raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)	    
          END DO
        END DO
        iNumGases = iNumGases + iNXsec
      END IF

      IF (iZbnd .GT. 0) THEN
        !!need to convert the user heights to user defined pressure levels
	CALL rspl(raAlt,raP,iNumLevs,raZbnd,raPbnd,iZbnd)
	DO iJ = 1,iZbnd
	  raPbnd(iJ) = raPbnd(iJ) * 100.0  !! change from mb to N/m2, as Psurf is finally in N/m2
	END DO
      END IF
      
 13   CONTINUE
      RETURN
      END
c************************************************************************      
