! Copyright 2014
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:42
 
! University of Maryland Baltimore County
! All Rights Reserved

! this file does the LBLRTM TAPE5/TAPE6

! see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
!       http://www.ssec.wisc.edu/~paulv/Fortran90/AtmProfile/Modules.html
!       https://svn.ssec.wisc.edu/repos/uwphysret/trunk/mfiles/compute_F.m
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

!      LBLRTM     iaGasUnits(iG) = 12   !!! assume hardcoded VMR
!      LBLRTM     pressures are in mb, which is what klayers wants

!      rP(min/max)KCarta is in N/m2 which is x100 what it would be in mb

! ESw.d ENw.d ESw.dEe ENw.dEe output format see http://www.enautica.pt/publico/professores/vfranco/formats.pdf

!************************************************************************
! http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
! this one takes the ppmv etc and applies to TOP bdry
!  (ie record 2.1.1 if TOA is level 1 ... GND at level L, this CORRECTLY applies ppmv to (L)
! ALTZ(L-1),  PZ(L-1),  TZ(L-1),  ATLZ(L),  PZ(L),  TZ(L)
!    ALTZ(L-1), PZ(L-1) and TZ(L-1) are only required for the first layer.
!    LBLRTM assumes that these quantites are equal to the top of the previous
!    layer for L > 1.

SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE5(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
    iNumLevs,rPSurf,rTSurf,rHSurf,iNumGases, raP,raT,raAlt,iZbnd,raPBnd,  &
    iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon)


NO TYPE, INTENT(IN OUT)                  :: caPFname
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
INTEGER, INTENT(IN)                      :: iNumLevs
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(IN)                         :: rTSurf
REAL, INTENT(IN OUT)                     :: rHSurf
INTEGER, INTENT(IN)                      :: iNumGases
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
REAL, INTENT(IN)                         :: raAlt(2*kProfLayer)
INTEGER, INTENT(OUT)                     :: iZbnd
REAL, INTENT(IN OUT)                     :: raPBnd(2*kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
REAL, INTENT(IN OUT)                     :: raaG_MR(2*kProfLayer,kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: rPMin
NO TYPE, INTENT(IN OUT)                  :: rPMax
REAL, INTENT(IN OUT)                     :: rYear
REAL, INTENT(IN OUT)                     :: rLat
REAL, INTENT(IN OUT)                     :: rLon
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INTEGER :: iPLEV
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

! input
CHARACTER (LEN=80) :: caPfName
REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta
! output

REAL :: rPmin,rPmax
REAL :: !
REAL :: raPavg(2*kProfLayer),raTavg(2*KProfLayer)  !!! TAPE5 can have layers average p AND T
REAL :: !! do we want user defined pressure level boundaries??


! local var
INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
REAL :: rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,raSumCheck(kMaxGas)
CHARACTER (LEN=80) :: caStr,caStrX,caStrY
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=2) :: c2
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=4) :: c4
CHARACTER (LEN=5) :: c5
INTEGER :: iWriteRTP,iNumGasesBAD,iaBadGasProfile(kMaxGas),iReplaceZeroProf
INTEGER :: IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS
INTEGER :: XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL
REAL :: rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,rDVOUT
REAL :: raEmiss(3),raRefl(3),rSecnto
INTEGER :: ILNFLG,iForm,iBmax,iOffSet
REAL :: raX(KMaxGas),rX,rSatZen
REAL :: raZbnd(2*kProfLayer)
INTEGER :: iDefault,iAIRS101_or_LBL_levels

iDefault = +1   !! use AIRS101 levels for the integration
iAIRS101_or_LBL_levels = +1 !! use AIRS101 levels for the integration
iAIRS101_or_LBL_levels = -1 !! use LBLRTM  levels for the integration
iAIRS101_or_LBL_levels = iaaOverrideDefault(3,1)
IF (ABS(iAIRS101_or_LBL_levels) > 1) THEN
  WRITE(kStdErr,*) 'invalid iAIRS101_or_LBL_levels ',iAIRS101_or_LBL_levels
  CALL DoStop
END IF
IF (iDefault /= iAIRS101_or_LBL_levels) THEN
  WRITE(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 when doing integration from levels to layers'
  WRITE(kStdErr,*) 'iDefault, iAIRS101_or_LBL_levels = ',iDefault,iAIRS101_or_LBL_levels
END IF

DO iL = 1,2*kProfLayer
  raPavg(iL) = -9999.0
  raTavg(iL) = -9999.0
END DO

rPmin = +1.0E6
rPmax = -1.0E+6
iNXsec = -1

IF (kProfLayer > kMaxLayer) THEN
  WRITE(kStdErr,*) 'oops, want to save LBLRTM temperatures into kLBLRTM_levelT'
  WRITE(kStdErr,*) 'but kProfLayer .GT. kMaLayer+1'
  CALL DoStop
ELSE
  DO iG = 1,kMaxLayer
    kLBLRTM_levelT(iG) = 0.0
    kLBLRTM_layerTavg(iG) = 0.0
  END DO
  kLBLRTM_levelT(kMaxLayer+1) = 0.0
END IF

iZbnd = -1   !!! assume we want to use default 101 AIRS levels
DO iL = 1,kProfLayer+1
  raPbnd(iL) = PLEV_KCARTADATABASE_AIRS(iL)
END DO

WRITE(kSTdWarn,*) 'Reading in LBLRTM TAPE5 .....'

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

!! read till you get the $ character, RECORD 1.1
10   CONTINUE
READ (iIOUN2,111,ERR=13,END=13) caStr
IF (caStr(1:1) /= '$') GO TO 10

READ (iIOUN2,111,ERR=13,END=13) caStr
CALL read_record_1p2(caStr,  &
    IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)
IF (ICNTNM == 6) THEN
  READ (iIOUN2,111,ERR=13,END=13) caStr
  CALL read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)
END IF
IF (IEMIT == 2) THEN
  CALL DOSTOPMesg('oops cannot handle solar file$')
END IF

IF ((IHIRAC > 0) .OR. (IAERSL > 0) .OR. (IEMIT == 1) .OR. (IATM == 1) .OR. (ILAS > 0)) THEN
  READ (iIOUN2,111,ERR=13,END=13) caStr
  CALL read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)
END IF

!! iotflg comes from solar data
!! IF ((IEMIT .EQ. 1) .OR. (IEMIT .EQ. 2 .AND.  IOTFLG .EQ. 1)) THEN
IF ((IEMIT == 1) .OR. (IEMIT == 2)) THEN
  READ (iIOUN2,111,ERR=13,END=13) caStr
  CALL read_record_1p4(caStr,rTSurf,raEmiss,raRefl)
END IF

111  FORMAT(A80)
WRITE(kStdWarn,*) 'TAPE 5 has iAtm = ',iAtm
IF (IATM == 0) THEN
!! need to read in profile
  CALL read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)
  CALL read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raPavg,raTavg,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck,  &
      rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta,  &
      iaG,iaGasUnits,iForm,iNumGases)
  IF (iAIRS101_or_LBL_levels == +1) THEN
    iZbnd = -1            !! use default AIRS 101 levels, they should span past the rPSurf
  ELSE
    iZbnd = iNumLevs  !! use the LBLRTM pressure levels
  END IF
  IF (iZbnd > 0) THEN
    raPbnd(1) = rPSurf
    iOffSet = kProfLayer-(iZbnd-1)
    DO iG = 1,iNumLevs
      IF (iG < iNumLevs) kLBLRTM_layerTavg(iG+iOffset) = raTavg(iG)
      kLBLRTM_levelT(iG+iOffSet) = raT(iG)
      raPbnd(iG) = raP(iG) * 100.0 !! change from mb to N/m2, as Psurf is finally in N/m2
!          raPbnd(iG) = raP(iG)         !! keep in mb
    END DO
  END IF
  IF (iXsect == 1) THEN
    CALL read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)
  END IF
ELSE
  rViewAngle = -9999.0
  rTopHgt    = 705000.0
  CALL read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)
  IF (iBmax > 0) THEN
    iZbnd = +iBmax
    CALL read_record_3p3b(iIOUN2,iBmax,raZbnd)
  END IF
  CALL read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf,  &
      rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
      raAlt,raP,raT,raSumCheck,  &
      iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon, iZbnd,raZbnd,raPbnd)
END IF
!!! TAPE5 can also include scanangle info and output file name info for flux, see Record 6,6.1
!!! Unless you want to read in rCos(flux angles) for the three angles, as yu can set them in nm_radnce, using
!!!   iAtmLoop = 3, raAtmLoop(1,2,3)

WRITE(kStdWarn,*) 'After exiting the equivalent of nolgas and xsecgas in TAPE5 LBLRTM, found',iNumGases,' gases'
WRITE(kStdWarn,*) '    these are the gasIDs'
WRITE(kStdWarn,*) (iaG(iJ),iJ=1,iNumGases)
WRITE(kStdWarn,*) '    these are the gas units'
WRITE(kStdWarn,*) (iaGasUnits(iJ),iJ=1,iNumGases)
WRITE(kStdWarn,*) ' '

!! now see if data was acutally read in
iNumGasesBAD = 0
DO iJ = 1,iNumGases
  IF (raSumCheck(iJ) < 1.0E-20) THEN
    WRITE(kStdWarn,*) ' read in TAPE5, following gas had    zero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
    iNumGasesBAD = iNumGasesBAD + 1
    iaBadGasProfile(iNumGasesBAD) = iJ
  ELSE
    WRITE(kStdWarn,*) ' readin TAPE5, following gas had nonzero profile column sum ',iJ,iaG(iJ),raSumCheck(iJ)
  END IF
END DO

!! need to fix the BAD gases
iDefault = +1
iReplaceZeroProf = -1    !! assume user knows why there is a ZERO everywhere gas profile
iReplaceZeroProf = +1    !! assume user wants to replace ZERO everywhere gas proFILE with climatology
iReplaceZeroProf = iaaOverrideDefault(3,2)
IF (ABS(iReplaceZeroProf) > 1) THEN
  WRITE(kStdErr,*) 'invalid iReplaceZeroProf ',iReplaceZeroProf
  CALL DoStop
END IF
IF ((iDefault /= iReplaceZeroProf) .AND. (iNumGasesBAD > 0)) THEN
  WRITE(kStdErr,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
  WRITE(kStdErr,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf
  WRITE(kStdWarn,*) 'in ReadInput_LBLRTM_ProfileTAPE5 : user wants to replace zero prof with climatology'
  WRITE(kStdWarn,*) 'iDefault = ',iDefault,' iReplaceZeroProf = ',iReplaceZeroProf
END IF

IF ((iNumGasesBAD > 0) .AND. (iReplaceZeroProf > 0)) THEN
  DO iG = 1,iNumGasesBAD
    CALL substitute_tape5_profile_for_climatology(iaBadGasProfile(iG),iNumLevs,raP,raaG_MR)
    iaGasUnits(iaBadGasProfile(iG)) = 10
    raSumCheck(iaBadGasProfile(iG)) = 0.0
    DO iJ = 1,iNumLevs
      raSumCheck(iaBadGasProfile(iG)) = raSumCheck(iaBadGasProfile(iG)) +  &
          raaG_MR(iJ,iaBadGasProfile(iG))
    END DO
    WRITE(kStdWarn,*) 'reset ZERO gas prof for gasID ',iaBadGasProfile(iG),' ppmv sum = ',raSumCheck(iaBadGasProfile(iG))
  END DO
END IF

13   CONTINUE
CLOSE(iIOUN2)
kProfileUnitOpen = -1

raRTP_TxtInput(1) = rPSurf
raRTP_TxtInput(2) = rTSurf
raRTP_TxtInput(3) = rHSurf    !! km
raRTP_TxtInput(4) = rTophgt   !! km
raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, downwelling radiation to instr, ELSE upwelling rad TO instr

rPSurf = rPSurf * 100.0       !! change from mb to N/m2
rHSurf = rHSurf * 1000.0      !! change from km to m

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
WRITE(kStdWarn,*) 'input file : Top alt (lowest press) = ',raP(iNumLevs),' mb at level ',iNumLevs
WRITE(kStdWarn,*) ' '

iWriteRTP = +21   !!!! raP,raT,raaG_MR,raAlt are in levels (g/g or ppmv or whatever) AND we DO have some layer avg info
  IF (iWriteRTP > 0) THEN
    CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,  &
        raP,raT,raaG_MR,raAlt,raPavg,raTavg,iWriteRTP)
  END IF
  
! now change all units first to PPMV (unit 10) and then to Volume Mix Ratio (unit 12)
! actually this is a waste of time as we go from VMR (units 12) to PPMV (units 10) back to VMR (unit 12)
  DO iG = 1,iNumGases
    CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR,+1)
    iaGasUnits(iG) = 12
    DO iL = 1,iNumLevs
      raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0E6
    END DO
  END DO
  
  RETURN
END SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE5

!************************************************************************
! this reads in TAPE6 which has the integrated layer amounts

SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE6(caPFname,rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta,  &
    iNumLays,rPSurf,rTSurf,rHSurf,iNumGases,raPX,raTX,raLayDensityX,raZX,  &
    raPressLevelsX,raTPressLevelsX,raAltitudesX,  &
    iaG,iaGasUnits,raaG_MRX,rPMin,rPMax,rYear,rLat,rLon)


NO TYPE, INTENT(IN OUT)                  :: caPFname
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
INTEGER, INTENT(OUT)                     :: iNumLays
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(OUT)                        :: rTSurf
REAL, INTENT(OUT)                        :: rHSurf
INTEGER, INTENT(OUT)                     :: iNumGases
REAL, INTENT(OUT)                        :: raPX(kProfLayer+1)
REAL, INTENT(OUT)                        :: raTX(kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: raLayDensi
REAL, INTENT(OUT)                        :: raZX(kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: raAltitude
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
REAL, INTENT(OUT)                        :: raaG_MRX(kProfLayer+1,kMaxGas)
NO TYPE, INTENT(OUT)                     :: rPMin
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
!!! these are LAYERS
REAL :: raLayDensityX(kProf

!!! this is LEVELS
REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)

! local var
INTEGER :: iIOUN2,iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iNumLevs,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
REAL :: raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH,r1,r2,r3,r4,r5,r6,rPTOA
CHARACTER (LEN=80) :: caStrX,caStrY
CHARACTER (LEN=120) :: caStr120,caStr
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=2) :: c2a,c2b
CHARACTER (LEN=128) :: cX
INTEGER :: iWriteRTP,iCountGases,iGCnt,iPass,iOffSet
INTEGER :: i1,i2,i3,i4,i5,i6
! for printing to lblrtm2rtp
REAL :: raPY(2*kProfLayer),raTY(2*kProfLayer),raaG_MRY(2*kProfLayer,kMaxGas),raZY(2*kProfLayer)
REAL :: raPavgJunk(2*kProfLayer),raTavgJunk(2*kProfLayer),raNavgJunk(2*kProfLayer),raSumGasAmtTAPE6(2*kProfLayer)
REAL :: raPZ(2*kProfLayer),raTZ(2*kProfLayer)

DO iL = 1,2*kProfLayer
  raSumGasAmtTAPE6(iL) = 0.0
END DO

rPmin = +1.0E6
rPmax = -1.0E+6
iNXsec = -1

WRITE(kSTdWarn,*) 'Reading in modified LBLRTM TAPE6 = 5 line summary of TAPE5 + MID TAPE6.....'
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

! this reads a modified version of TAPE6 so that the following info is there
! >>>>>>>>>>>>>>>>>>> simple header made by YOU
! short description of your yada yada yada
! rF1 rF2
! rTSurf,rPSurf,rPTOA
! iNumLevs iNumGases iNXsec
! rViewAngle
! >>>>>>>>>>>>>>>>> cut from TAPE5
! TAPE5 INFO
! for i = 1 : iNumlayes
!   read string1
!   read string2
! end do
! END TAPE5 INFO next line is start of TAPE6 info
! >>>>>>>>>>>>>>>>> cut from TAPE6
! caStr starting with 0LAYER and then all the subsequent info

!---------------------- read in header ------------------------
WRITE(kStdWarn,*) '   ... reading in 5 line summary of TAPE5/TAPE6 file ....'
READ (iIOUN2,111,ERR=222,END=222) caStr
READ (iIOUN2,*) rF1,rF2
WRITE(kStdWarn,*) 'LBLRTM input indicates start/stop wavenumbers are ',rF1,rF2

READ (iIOUN2,*) rTSurf,rPSurf,rPTOA
WRITE(kStdWarn,*) 'LBLRTM input indicates STEMP/SPRES and TOA_Pres are ',rTSurf,rPSurf,rPTOA

READ (iIOUN2,*) iNumLevs,iNumGases,iNXSec
WRITE(kStdWarn,*) 'TAPE5 implies kCARTA compute ODs at ',iNumLevs,' press/heights, for iNumGases = ',iNumGases
IF (iNumLevs > kProfLayer) THEN
  WRITE(kStdErr,*) 'though irrelevant for kCARTA calcs, it will not be able to read in the pressures/heights'
  WRITE(kStdErr,*) iNumLevs,kProfLayer
  CALL DoStop
END IF
IF (iNumGases+iNXsec > kMaxGas) THEN
  WRITE(kStdErr,*) 'iNumGases+iNXsec .GT. kMaxGas',iNumGases,iNXsec,kMaxGas
  CALL DoStop
END IF

READ (iIOUN2,*) rViewAngle
iLBROutBdryHorP = +1   !! assume start/stop heigt info given

iNumLays = iNumLevs - 1
IF (iNumLays > kProfLayer) THEN
  WRITE(kStdErr,*) 'iNumLays .GT. kProfLayer',iNumLays,kProfLayer
  CALL DoStop
END IF

!---------------------- read in TAPE5 info (for plevs)  -----------
! see subr read_record_2p1p1
WRITE(kSTdWarn,*) ' '
WRITE(kStdWarn,*) '   ... reading in relevant edited TAPE 5 info ....'
READ (iIOUN2,111,ERR=222,END=222) caStr   !!! start TAPE5 info
DO iL = 1,iNumLevs-1
  READ (iIOUN2,111,ERR=222,END=222) caStr120
  caStr = caStr120(1:40)
  READ(caStr,*) raPAvgJunk(iL),raTavgJunk(iL)
  caStr = caStr120(41:120)
  IF (iL. EQ. 1) THEN
    READ (caStr,*) raAltitudesX(1),raPressLevelsX(1),raTPressLevelsX(1),raAltitudesX(2),raPressLevelsX(2),raTPressLevelsX(2)
  ELSE
    READ (caStr,*) raAltitudesX(iL+1),raPressLevelsX(iL+1),raTPressLevelsX(iL+1)
  END IF
  READ (iIOUN2,111,ERR=222,END=222) caStr  !! mixing ratios, not relevant
!! pV = nRT ==> pdz = n/Vdz RT ==> p dz/RT = q = moles/m3, which can be compared to what is in TAPE6!!!!
!! change mb to N/m2 and km to m
  raNAvgJunk(iL) = raPavgJunk(iL)*100.0 * (raAltitudesX(iL+1)-raAltitudesX(iL)) * 1000.0
  raNAvgJunk(iL) = raNAvgJunk(iL) / kMGC /raTavgJunk(iL)    !!! this is in mles/m2
  raNAvgJunk(iL) = raNAvgJunk(iL) * 1E-4 * kAvog/1000.0
END DO

kLBLRTM_TOA = raPressLevelsX(iNumLevs)
WRITE(kStdWarn,*) 'kLBLRTM_toa = ',kLBLRTM_toa,' mb (when used with flux calcs)'
WRITE(kStdWarn,*) ' '

iOffSet = kProfLayer+1 - (iNumLevs)
DO iL = 1,iNumLevs
  kLBLRTM_levelT(iL+iOffSet) = raTPressLevelsX(iL)
END DO

READ (iIOUN2,111,ERR=222,END=222) caStr   !!! end TAPE5 info

!---------------------- read in TAPE6 info ------------------------

! now start reading in the TAPE6 info
! first set gasunits = 1 = molecules/cm2
WRITE(kSTdWarn,*) ' '
WRITE(kStdWarn,*) '   ... reading in relevant edited TAPE 6 info ....'
!! we now read every gas 1 .. iNumGases  in edited TAPE6 == molecules/cm2
DO iG = 1,iNumGases
  iaG(iG) = iG
  iaGasUnits(iG) = 1   !!! assume hardcoded molecules/cm2
END DO

READ(iIOUN2,111) caStr
IF (caStr(1:6) /= '0LAYER') THEN
  WRITE(kStdErr,*) 'oops, expecting string : 0LAYER                          P(MB)       T(K)    ALPHL    ALPHD    ALPHV ...'
  WRITE(kStdErr,*) 'but instead got ',caSTr
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
  IF (iL == 1)        rHSurf = r1
  IF (iL == iNumLays) rTopHgt = r2
!        print *,iL,iNumLays,raZX(iL),rP,rT
END DO

rPmax = rPSurf
rPMin = rPTOA

iPass = 1
iCountGases = 0
iGCnt = MIN(7,iNumGases)  !! in first pass, read at most 7 gases
iGCnt = 8   !! reads first gases, then "all other gases"
!!! might need to re-code this because of OTHER
!!!  H2O           CO2            O3           N2O            CO           CH4            O2            >>> OTHER <<<
!!! because of OTHER need to read at least iNumGases+1  profiles
!!! might need to re-code this because of OTHER
iGCnt = MIN(8,iNumGases+1)   !! reads first gases, then "all other gases"
WRITE(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at FIRST PASS ',iPass
WRITE(kStdWarn,*) '  because the first 7 are regular profile and eight is "OTHER" '
READ(iIOUN2,111) caStr      !! blank except for 1 at beginning
READ(iIOUN2,111) caStr      !! eg LBLRTM    14/07/04  22:24:32
READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O
DO iL = 1,iNumLays
  READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
  DO iG=1,iGCnt
    raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
  END DO
  raZX(iL)   = r1
  raZX(iL+1) = r2
!        raPX(iL) = rP * 100.0  !! change from mb to N/m2
  raPX(iL) = rP          !! keep in mb
  raTX(iL) = rT
  IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
  IF (rPmin > raPX(iL)) rPmin = raPX(iL)
  DO iG=1,iGCnt-1
    raaG_MRX(iL,iG+iCountGases) = raX(iG)
  END DO
!      print *,iL,'read in molgas',raaG_MRX(iL,1),raaG_MRX(iL,2),raaG_MRX(iL,3)
  raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
END DO
READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM
iCountGases = iCountGases + (iGCnt-1)   !!! because of >>>> OTHER <<<<

iOffset = kProfLayer - iNumLays
DO iL = 1,iNumLays
  kLBLRTM_layerTavg(iL+iOffSet) = raTX(iL)
END DO

15   CONTINUE
IF ((iNumGases == iCountGases) .AND. (iNXsec == 0)) GO TO 222     !!! done, just CLOSE FILE!!!!
IF ((iNumGases == iCountGases) .AND. (iNXsec > 0)) GO TO 23      !!! done with main gases, need TO READ xsec gases

! need to do this if need to read in more gases
iPass = iPass + 1
iGCnt = MIN(7,iNumGases+1-iCountGases)
!!! might need to re-code this because of OTHER
!!!  H2O           CO2            O3           N2O            CO           CH4            O2            >>> OTHER <<<
!!! because of OTHER need to read at least iNumGases+1  profiles
!!! might need to re-code this because of OTHER
WRITE(kSTdWarn,*) '  Reading ',iGCnt,' gases from TAPE6 at pass ',iPass
READ(iIOUN2,111) caStr      !! eg 1blank
READ(iIOUN2,111) caStr      !! LBLRTM    14/07/04  22:24:32
READ(iIOUN2,111) caStr      !! eg MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER
READ(iIOUN2,111) caStr      !! eg P(MB)      T(K)   IPATH         H2O
DO iL = 1,iNumLays
  READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iGCnt)
  DO iG=1,iGCnt
    raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
  END DO
  raZX(iL)   = r1
  raZX(iL+1) = r2
!        raPX(iL) = rP * 100.0  !! change from mb to N/m2
  raPX(iL) = rP          !! keep in mb
  raTX(iL) = rT
  IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
  IF (rPmin > raPX(iL)) rPmin = raPX(iL)
  DO iG=1,iGCnt
    raaG_MRX(iL,iG+iCountGases) = raX(iG)
  END DO
!!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
END DO
iCountGases = iCountGases + iGCnt
READ(iIOUN2,111) caStr      !! ACCUMULATED MOLECULAR AMOUNTS FOR TOTAL PATH
READ(iIOUN2,111) caStr      !! eg 0 97  0.000 TO 82.724 KM

GO TO 15

23   CONTINUE    !! finished reading in LAYER amounts for main gases, now read in AUX CUMULATIVE INFO
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
!      print *,caStr
CALL XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,kRTP)
IF (iNXsec > 7) THEN
  WRITE(kStdWarn,*) 'oops assumed at most 7 xsec gases, so only 1 pass ... need to modify code for > 7 gases'
  CALL DOStop
END IF
DO iL = 1,iNumLays
  READ (iIOUN2,*) i1,i2,r1,c2a,r2,c2b,rP,rT,i6,(raX(iG),iG=1,iNXsec)
  DO iG=1,iGCnt
    raSumGasAmtTAPE6(iL) = raSumGasAmtTAPE6(iL) + raX(iG)
  END DO
  raZX(iL)   = r1
  raZX(iL+1) = r2
!        raPX(iL) = rP * 100.0  !! change from mb to N/m2
  raPX(iL) = rP          !! keep in mb
  raTX(iL) = rT
  IF (rPmax <= raPX(iL)) rPmax = raPX(iL)
  IF (rPmin > raPX(iL)) rPmin = raPX(iL)
  DO iG=1,iNXsec
    raaG_MRX(iL,iG+iNumGases) = raX(iG)
  END DO
!      print *,iL,'read in xscgas',raaG_MRX(iL,1),raaG_MRX(iL,2),raaG_MRX(iL,3)
!!!! raLayDensityX(iL) = raX(7)*100.0/20.9    !! turn OXYGEN amount into proxy for AIR, assuming VMR of 0.209
END DO

222  CONTINUE
CLOSE(iIOUN2)
kProfileUnitOpen = -11
iNumGases = iNumGases + iNXsec

111  FORMAT(A120)
112  FORMAT(I3,2(' ',F10.3), 2(' ',F10.3),5(' ',ES12.5),1(' ',F12.5))
113  FORMAT(A122)

WRITE(kStdWarn,*) '... finshed reading in TAPE5/TAPE6 combo ....'
WRITE(kStdWarn,*) ' '
WRITE(kSTdWarn,*) '     rTSurf,rPSurf,rPTOA = ',rTSurf,rPSurf,rPTOA
cX =  &
    ' Index   Plev      Tlev       Pavg       Tavg     WV(amt)       CO2(amt)   O3(amt)    || q=pdz/RT   totalQ(TAPE6)    Ratio'
WRITE(kStdWarn,113) cX
cX =  &
    '--------------------------------------------------------------------------------------||-----------------------------------'

WRITE(kStdWarn,113) cX
DO iL = 1,iNumLays+1
  WRITE(kSTdWarn,112) iL,raPressLevelsX(iL),raTPressLevelsX(iL),raPX(iL),raTX(iL),(raaG_MRX(iL,iJ),iJ=1,3),  &
      raNAvgJunk(iL),raSumGasAmtTAPE6(iL),raNAvgJunk(iL)/raSumGasAmtTAPE6(iL)
END DO
WRITE(kStdWarn,113) cX

raRTP_TxtInput(1) = rPSurf
raRTP_TxtInput(2) = rTSurf
raRTP_TxtInput(3) = rHSurf    !! km
raRTP_TxtInput(4) = rTophgt   !! km
raRTP_TxtInput(5) = rViewAngle  !! if 0 < ang < 90, then downwell rad, else upwell radn

rPSurf = rPSurf * 100.0       !! chnage from mb to N/m2
rHSurf = rHSurf * 1000.0      !! change from km to m

WRITE(kStdWarn,*)'  KCARTA Database : max/min press (mb) = ',rPmaxKCarta/100.0,rPminKCarta/100.0
WRITE(kStdWarn,*)'  kCARTA Database : max/min height (m) = ',rHmaxKCarta,rHminKCarta
WRITE(kStdWarn,*)'input file : spres/sHeight      = ',rPSurf,rHSurf

! make sure pressures are decreasing with index ie layers going higher and higher
IF (raPX(1) < raPX(2)) THEN
!!! need to swap!!!
  PRINT *,'swappy doody doo'
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
WRITE(kStdWarn,*) 'input file : Highest altitude (lowest <press>) = ',raPX(iNumLays),' mb at layer ',iNumLays
WRITE(kStdWarn,*) 'input file : Highest altitude (lowest plev)    = ',rPTOA,' mb at level ',iNumLays+1
WRITE(kStdWarn,*) ' '

iWriteRTP = +1   !!! this is in layers
IF (iWriteRTP > 0) THEN
  DO iL = 1,iNumLevs
    raPZ(IL) = raPX(iL)             !!! average pressure
    raTZ(iL) = raTPressLevelsX(iL)  !!! levels pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    raTY(IL) = raTX(iL)             !!! layers temp
    raZY(IL) = raZX(iL)             !!! layers alt
    raPY(IL) = raPressLevelsX(iL)        !!! always need levels info for p.plevs
    raZY(IL) = raAltitudesX(iL)        !!! levels alt
    DO iG = 1,iNumGases
      raaG_MRY(iL,iG) = raaG_MRX(iL,iG)
    END DO
  END DO
  CALL lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,  &
      raPY,raTY,raaG_MRY,raZY,raPZ,raTZ,iWriteRTP)
END IF

! now change all units to MR
!      DO iG = 1,iNumGases
!        CALL changeLVLS_2_ppmv(iaG(iG),iaGasUnits(iG),iNumLevs,iG,raP,raT,raaG_MR,+1)
!        DO iL = 1,iNumLevs
!          raaG_MR(iL,iG) =  raaG_MR(iL,iG) / 1.0e6
!          iaGasUnits(iG) = 12
!        END DO
!      END DO

RETURN
END SUBROUTINE  ReadInput_LBLRTM_ProfileTAPE6

!************************************************************************
! parses LBLRTM xsec names and finds gas ids

SUBROUTINE XsecNamesLBL(caStr,iaG,iaGasUnits,iNumGases,iNXsec,iLBLTapeType)


CHARACTER (LEN=80), INTENT(IN)           :: caStr
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN)                      :: iNXsec
NO TYPE, INTENT(IN OUT)                  :: iLBLTapeTy
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input

INTEGER :: iLBLTapeType
! input/output


! local
INTEGER :: iaX(kMaxGas),iG,i1,i2
CHARACTER (LEN=10) :: caStr10
CHARACTER (LEN=15) :: caStr15

IF (iLBLTapeType == -5) THEN
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
  
ELSE IF (iLBLTapeType == -6) THEN
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
END SUBROUTINE XsecNamesLBL

!************************************************************************
! this maps xsecnames to xsec IDs

SUBROUTINE mapXsecname_to_XsecID(caStr10,iG,iaX)


CHARACTER (LEN=10), INTENT(IN)           :: caStr10
INTEGER, INTENT(IN OUT)                  :: iG
INTEGER, INTENT(OUT)                     :: iaX(kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input


! input/output
INTEGER :: !! modify xsec ID according to strcmp

! local
INTEGER :: iI,iFound,iChecked
CHARACTER (LEN=10) :: caaNamesA(13)
CHARACTER (LEN=10) :: caaNamesB(13)
CHARACTER (LEN=10) :: caaNamesC(13)
CHARACTER (LEN=10) :: caX

IF (caStr10(1:1) /= ' ') THEN
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
IF ((caaNamesA(iChecked) == caX) .OR. (caaNamesB(iChecked) == caX) .OR. (caaNamesC(iChecked) == caX)) THEN
  iFound = iChecked
END IF

IF ((iFound < 0) .AND. (iChecked < 13)) THEN
  iChecked = iChecked + 1
  GO TO 10
END IF

IF (iFound < 0) THEN
  WRITE(kStdErr,*) 'Reading LBLRTM file; could not find xsec match for ',caX
  CALL DOStop
ELSE
  iaX(iG) = 50 + iChecked
  WRITE(kStdWarn,*) '  LBLRTM xsec gas = ',caX,' corresponds to gasID ',iaX(iG)
END IF

RETURN
END SUBROUTINE mapXsecname_to_XsecID

!                           -------------------------------------------------
!                             Alias(1)           Alias(2)           Alias(3)           Alias(4)
!                            ----------         ----------         ----------         ----------
!                            CLONO2              CLNO3
!                            HNO4
!                            CHCL2F                                 CFC21              F21
!                            CCL4
!                            CCL3F               CFCL3              CFC11              F11
!                            CCL2F2              CF2CL2             CFC12              F12
!                            C2CL2F4             C2F4CL2            CFC114             F114
!                            C2CL3F3             C2F3CL3            CFC113             F113
!                            N2O5
!                            HNO3
!                            CF4                                    CFC14              F14
!                            CHCLF2              CHF2CL             CFC22              F22
!                            CCLF3                                  CFC13              F13
!                            C2CLF5                                 CFC115             F115

!************************************************************************
! this writes the LBLRTM input so you can cut and paste into Matlab, and save RTP file
! now write the klayers stuff

SUBROUTINE lblrtm2rtp(rF1,rF2,rPmin,rPmax,iNumGases,iaG,iaGasUnits,iNumLevs,rPSurf,rTSurf,rHSurf,  &
    raP,raT,raaG_MR,raAlt,raJunkP,raJunkT,iWriteRTP)


REAL, INTENT(IN)                         :: rF1
REAL, INTENT(IN)                         :: rF2
REAL, INTENT(IN)                         :: rPmin
REAL, INTENT(IN)                         :: rPmax
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iaGasUnits(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumLevs
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(IN)                         :: rTSurf
REAL, INTENT(IN)                         :: rHSurf
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
REAL, INTENT(IN)                         :: raaG_MR(2*kProfLayer,kMaxGas)
REAL, INTENT(IN OUT)                     :: raAlt(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raJunkP(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raJunkT(2*kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iWriteRTP
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input

INTEGER :: iWriteRTP    !! +1 for levels, +2 for layers profiles

REAL :: !! this is always LEVELS press
REAL :: !! these could be lay
REAL :: !! these change as iWri
!! from 1 (layers) to 21 (levels)

! local var
REAL :: raX(2*kProfLayer)
INTEGER :: iL,iG,iLeftjust_lenstr
CHARACTER (LEN=8) :: ca8
CHARACTER (LEN=6) :: ca6
CHARACTER (LEN=2) :: ca2

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' '

WRITE(kStdWarn,*) 'h.vcmin = ',rF1,';'
WRITE(kStdWarn,*) 'h.vcmax = ',rF2,';'
!      write(kStdWarn,*) 'h.pmin = ',rPmin/100.0,';      %% rPmin in N/m2 --> mb'
!      write(kStdWarn,*) 'h.pmax = ',rPmax/100.0,';      %% rPmax in N/m2 -->  mb'
WRITE(kStdWarn,*) 'h.pmin = ',rPmin,';      %% in mb'
WRITE(kStdWarn,*) 'h.pmax = ',rPmax,';      %% im mb'
WRITE(kStdWarn,*) 'h.ngas = ',iNumGases,';'
IF (iaGasUnits(1) > 1) THEN
  IF (iWriteRTP <= 1) THEN
    WRITE(kStdErr,*) 'trying to write out LEVELS profile but found inconsistency ...'
    CALL DoStop
  END IF
  WRITE(kStdWarn,*) 'h.ptype = 0;'
  WRITE(kStdWarn,*) 'h.pfields = 1;'
ELSE IF (iaGasUnits(1) == 1) THEN
  IF (iWriteRTP /= 1) THEN
    WRITE(kStdErr,*) 'trying to write out LAYERS profile but found inconsistency ...'
    CALL DoStop
  END IF
  WRITE(kStdWarn,*) 'h.ptype = 1;'
  WRITE(kStdWarn,*) 'h.pfields = 1;'
END IF
WRITE(kStdWarn,*) 'h.glist = [',(iaG(iG),iG=1,iNumGases),']'';'
WRITE(kStdWarn,*) 'h.gunit = [',(iaGasUnits(iG),iG=1,iNumGases),']'';'

WRITE(kStdWarn,*) 'p.nlevs = ',iNumLevs,';'
WRITE(kStdWarn,*) 'p.spres = ',rPsurf/100.0,';  %% in mb'
WRITE(kStdWarn,*) 'p.stemp = ',rTSurf,'; %%%%% if kSurfTemp < 0, o/w use raStemp from nm_radnce'
WRITE(kStdWarn,*) 'p.salti = ',rHSurf,'; %%%% WOWOWOWOWOW'

WRITE(kStdWarn,*) 'p.satzen = 0.0;'
WRITE(kStdWarn,*) 'p.scanang = 0.0;'
WRITE(kStdWarn,*) 'p.solzen = 130.0;'
WRITE(kStdWarn,*) 'p.upwell = 1;'
WRITE(kStdWarn,*) 'p.zobs = 705000.0;'
WRITE(kStdWarn,*) 'p.nemis = 2;'
WRITE(kStdWarn,*) 'p.efreq = [600 3000]'';  %% CHECK THIS WITH LBLRTM EFREQ'
WRITE(kStdWarn,*) 'p.emis = [1.0 1.0]'';    %% AND AGAINST NM_RADNCE EMISFILE AND RAEMISS'
WRITE(kStdWarn,*) 'p.rho = [0.0 0.0]'';     %% AND AGAINST NM_RADNCE RARHO'

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
IF (iWriteRTP == 1) THEN
  WRITE(kStdWarn,*) '% we have read in TAPE6 (layers averages) so have average lay pressures as well'
  ca8 = 'p.plays'
  CALL write_stringnice(ca8,raJunkP,1.00,8,iNumLevs,-1)       !! write out LEVELS T, as you have written out raT which is layers T
  
  WRITE(kStdWarn,*) '% we have read in TAPE6 (layers averages) but also have TAPE5 (levels T info)'
  ca8 = 'p.tlevs'
  CALL write_stringnice(ca8,raJunkT,1.00,8,iNumLevs,-1)       !! write out LEVELS T, as you have written out raT which is layers T
ELSE IF (iWriteRTP == 21) THEN
  WRITE(kStdWarn,*) '% we have read in TAPE5 (level info) but have some lay pressures as well'
  ca8 = 'p.plays'
  CALL write_stringnice(ca8,raJunkP,1.00,8,iNumLevs,-1)
  
  WRITE(kStdWarn,*) '% we have read in TAPE5 (level info) but have some lay temps as well'
  ca8 = 'p.tlays'
  CALL write_stringnice(ca8,raJunkT,1.00,8,iNumLevs,-1)
END IF
WRITE(kStdWarn,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
WRITE(kStdWarn,*) ' '

ca8 = 'p.plevs'
!CALL write_stringnice(ca8,raP,0.01,8,iNumLevs,-1)       !! raP in N/m2, convert to mb for rtp
CALL write_stringnice(ca8,raP,1.00,8,iNumLevs,-1)       !! raP in mb, keep in mb for rtp
ca8 = 'p.ptemp'
CALL write_stringnice(ca8,raT,1.00,8,iNumLevs,+1)
ca8 = 'p.palts'
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

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' % >>>>>>>>>>>>>>>>> LBLRTM --> RTP file cut >>>>>>>>>>>>>>>>>>>'
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' '

RETURN
END SUBROUTINE lblrtm2rtp

!************************************************************************
! this simply writes the strings nicely

SUBROUTINE write_stringnice(ca8,raX,rScale,iNumItemsPerLine,iNumLevs,ForE)


CHARACTER (LEN=8), INTENT(OUT)           :: ca8
REAL, INTENT(IN)                         :: raX(2*kProfLayer)
REAL, INTENT(IN)                         :: rScale
NO TYPE, INTENT(IN OUT)                  :: iNumItemsP
INTEGER, INTENT(IN)                      :: iNumLevs
INTEGER, INTENT(IN OUT)                  :: ForE
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input

INTEGER :: iNumItemsPerLine


! local var
INTEGER :: iL,iY,iFull,iPart,i1,i2,iFloor,iLen
REAL :: raY(2*kProfLayer)

iFull = ifloor(iNumLevs*1.0/iNumItemsPerLine)

DO iL = 1,iNumLevs
  raY(iL) = raX(iL) * rScale
END DO

IF (iNumItemsPerLine /= 8) iNumItemsPerLine = 8

WRITE(kStdWarn,*) ca8,' = [...'
DO iL = 1,iFull
  i1 = 1 + (iL-1)* iNumItemsPerLine
  i2 = i1 + iNumItemsPerLine - 1
!        print *,iL,i1,i2
  IF (ForE == +1) THEN
    WRITE(kStdWarn,1108) (raY(iY),iY=i1,i2)
  ELSE IF (ForE == -1) THEN
    WRITE(kStdWarn,1308) (raY(iY),iY=i1,i2)
  END IF
END DO
IF (i2 < iNumLevs) THEN
  i1 = i2 + 1
  i2 = iNumLevs
!        print *,9999,i1,i2
  iLen = i2-i1+1
  IF (iLen > iNumItemsPerLine) THEN
    WRITE(kSTdErr,*) 'cannot have more than 8 array members per line'
    CALL DoStop
  END IF
  IF (ForE == +1) THEN
    IF (iLen == 1) THEN
      WRITE(kStdWarn,1101) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 2) THEN
      WRITE(kStdWarn,1102) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 2) THEN
      WRITE(kStdWarn,1102) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 3) THEN
      WRITE(kStdWarn,1103) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 4) THEN
      WRITE(kStdWarn,1104) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 5) THEN
      WRITE(kStdWarn,1105) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 6) THEN
      WRITE(kStdWarn,1106) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 7) THEN
      WRITE(kStdWarn,1107) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 8) THEN
      WRITE(kStdWarn,1108) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 9) THEN
      WRITE(kStdWarn,1109) (raY(iY),iY=i1,i2)
    ELSE
      WRITE(kStdErr,*) 'oops so many numbers in the line????'
      CALL DoStop
    END IF
  ELSE IF (ForE == -1) THEN
    IF (iLen == 1) THEN
      WRITE(kStdWarn,1301) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 2) THEN
      WRITE(kStdWarn,1302) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 2) THEN
      WRITE(kStdWarn,1302) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 3) THEN
      WRITE(kStdWarn,1303) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 4) THEN
      WRITE(kStdWarn,1304) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 5) THEN
      WRITE(kStdWarn,1305) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 6) THEN
      WRITE(kStdWarn,1306) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 7) THEN
      WRITE(kStdWarn,1307) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 8) THEN
      WRITE(kStdWarn,1308) (raY(iY),iY=i1,i2)
    ELSE IF (iLen == 9) THEN
      WRITE(kStdWarn,1309) (raY(iY),iY=i1,i2)
    ELSE
      WRITE(kStdErr,*) 'oops so many numbers in the line????'
      CALL DoStop
    END IF
  END IF
END IF
WRITE(kStdWarn,*) ']'';'

1101 FORMAT(1(' ',F12.5),' ...')
1102 FORMAT(2(' ',F12.5),' ...')
1103 FORMAT(3(' ',F12.5),' ...')
1104 FORMAT(4(' ',F12.5),' ...')
1105 FORMAT(5(' ',F12.5),' ...')
1106 FORMAT(6(' ',F12.5),' ...')
1107 FORMAT(7(' ',F12.5),' ...')
1108 FORMAT(8(' ',F12.5),' ...')
1109 FORMAT(9(' ',F12.5),' ...')

1301 FORMAT(1(' ',ES12.5),' ...')
1302 FORMAT(2(' ',ES12.5),' ...')
1303 FORMAT(3(' ',ES12.5),' ...')
1304 FORMAT(4(' ',ES12.5),' ...')
1305 FORMAT(5(' ',ES12.5),' ...')
1306 FORMAT(6(' ',ES12.5),' ...')
1307 FORMAT(7(' ',ES12.5),' ...')
1308 FORMAT(8(' ',ES12.5),' ...')
1309 FORMAT(9(' ',ES12.5),' ...')

RETURN
END SUBROUTINE write_stringnice

!************************************************************************
! this subroutine just takes the LBLRTM mol/cm2 TAPE6 profile and does some simple calcs

SUBROUTINE DoLBLRTMLayers2KCARTALayers(rHSurf,rPSurf,rTSurf,iNumLays,iNumGases,iaG,rLat,rLon,  &
    PAVG_KCARTADATABASE_AIRS,PLEV_KCARTADATABASE_AIRS,DATABASELEVHEIGHTS,rfracBot,  &
    raPX,raTX,raLayDensityX,raaG_MRX,  &
    raPressLevelsX,raTPressLevelsX,raAltitudesX,  &
    raPout,raAmountOut,raTout,raZout,raaQout,raaPartPressOut,iLowestLev)


REAL, INTENT(IN OUT)                     :: rHSurf
REAL, INTENT(OUT)                        :: rPSurf
REAL, INTENT(IN OUT)                     :: rTSurf
INTEGER, INTENT(IN)                      :: iNumLays
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaG(kMaxGas)
REAL, INTENT(IN OUT)                     :: rLat
REAL, INTENT(IN OUT)                     :: rLon
NO TYPE, INTENT(IN OUT)                  :: PAVG_KCART
NO TYPE, INTENT(IN OUT)                  :: PLEV_KCART
NO TYPE, INTENT(IN OUT)                  :: DATABASELE
REAL, INTENT(OUT)                        :: rfracBot
REAL, INTENT(IN)                         :: raPX(kProfLayer+1)
REAL, INTENT(IN)                         :: raTX(kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: raLayDensi
NO TYPE, INTENT(IN OUT)                  :: raaG_MRX
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: raAltitude
REAL, INTENT(OUT)                        :: raPout(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raAmountOu
REAL, INTENT(OUT)                        :: raTout(kProfLayer)
REAL, INTENT(OUT)                        :: raZout(kProfLayer+1)
REAL, INTENT(OUT)                        :: raaQout(kProfLayer,kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
INTEGER, INTENT(OUT)                     :: iLowestLev
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input
REAL :: PAVG_KCARTADATABASE_AIRS(kMaxLayer),PLEV_KCARTADATABASE_
REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)

REAL :: raLayDensityX(kProfLayer+1),raaG_MRX(kProfLayer+1,
! output
REAL :: raAmountOut(kProfLayer)
REAL :: raaPartPressOut(kProfLayer,kMaxGas)
REAL :: raTPressLevelsX(kProfLayer+1),raPressLevelsX(kProfLayer+1),raAltitudesX(kProfLayer+1)


! local
INTEGER :: iL,iCnt,iJ,iNFine,iSmallLoop,iG,iLoop,iUpperLev
REAL :: z,rP,rPTry,rdPWant,rdP,rT,amount,slope,dz,z0,gravity,gamma,junk,rTsum,rWgt,damount
REAL :: rMolarMass_n,rMolarMass_np1,rMR_water_n,rMR_water_np1,rPP,rSVP,rAmt,rThick
REAL :: dlnp,rP_n,rP_np1,rT_n,rT_np1,rTX,rPX,density_n,density_np1,amount_n,amount_np1,Pav,Tav,SpHeatDryAir
REAL :: raXYZPress(kMaxlayer+1+1)
REAL :: raXYZTemp(kMaxlayer+1+1)
REAL :: raaXYZ_MR(kMaxlayer+1+1,kMaxGas)
REAL :: raXYZ_MRwater(kMaxlayer+1+1)
REAL :: raJunk(kMaxLayer),rMR_n,rMR_np1,q_n,q_np1,raJunk2(kMaxLayer)
REAL :: zWoo,rConvertQ,grav_earth,wexsvp,rRH
INTEGER :: iWoo

DO iL = 1,kProfLayer
  DO iG = 1,iNumGases
    raaQout(iL,iG) = 0.0
  END DO
END DO

IF ((kPlanet == 3) .AND. (iaG(1) /= 1)) THEN
  WRITE (kStdErr,*) 'Need GasID(1) = 1 (WV) for Planet Earth'
  CALL DoStop
END IF

iLowestLev = kProfLayer - iNumLays + 1

! d onot need this, as we are using TAPE6 levels
!      rFracBot = (rPSurf-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))/
!     $         (PLEV_KCARTADATABASE_AIRS(iLowestLev)-PLEV_KCARTADATABASE_AIRS(iLowestLev+1))
! just readjust the layer info read in, and tack on additional gas profiles from refgas
rFracBot = 1.0
WRITE(kSTdWarn,*) 'iLowestLev = ',iLowestLev,' rfracBot = ',rFracBot,'( if less than 1, reset to 1)'
WRITE(kStdWarn,*) 'PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev) = ',  &
    PLEV_KCARTADATABASE_AIRS(iLowestLev+1),rPSurf,PLEV_KCARTADATABASE_AIRS(iLowestLev)

DO iL = kProfLayer+1,kProfLayer+1
  raZout(iL) = raAltitudesX(iL-iLowestLev+1)*1000   !! change to meters
END DO
DO iL = iLowestLev,kProfLayer
  raZout(iL) = raAltitudesX(iL-iLowestLev+1)*1000   !! change to meters
  raPout(iL) = raPX(iL-iLowestLev+1) * 100.0        !! change mb to N/m2
  raTout(iL) = raTX(iL-iLowestLev+1)
  raAmountOut(iL) = raLayDensityX(iL-iLowestLev+1)
  IF (iL == iLowestLev) raAmountout(iL) = raAmountout(iL)/rFracBot
  DO iG = 1,iNumGases
    raaQout(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)
    IF (iL == iLowestLev) raaQout(iL,iG) = raaQout(iL,iG)/rFracBot
    raaPartPressOut(iL,iG) = raaG_MRX(iL-iLowestLev+1,iG)/raLayDensityX(iL-iLowestLev+1) * raPX(iL-iLowestLev+1)*100
  END DO
END DO

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

!      DO iL = iLowestLev,kProfLayer+1
!        raPressLevels(iL)  = raPressLevelsX(iL-iLowestLev+1)
!        raTPressLevels(iL) = raTPressLevelsX(iL-iLowestLev+1)
!      END DO

IF (kPlanet == 3) THEN
!! show the RH for gas 1
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*) 'Checking Relative Humidity'
  WRITE(kStdWarn,113) ' iX   Lay  P(mb)    zav(km)  dz(km)   T(K)     PP(mb)  SVP(mb)    RH      QTot      Q1..Qn'
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
!          raRTP_TxtInput(6) = zWoo/100.0  !! mb
    raRTP_TxtInput(6) = zWoo        !! keep in mb
  ELSE IF (raRTP_TxtInput(6) < 0) THEN
!!! input level boundaries in mb
!          raRTP_TxtInput(6) = abs(raRTP_TxtInput(6)) * 100.0
!          write(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' N/m2'
    raRTP_TxtInput(6) = ABS(raRTP_TxtInput(6))
    WRITE(kStdWarn,*)'LBLRTM output pressure of ',raRTP_TxtInput(6),' mb'
  END IF
END IF

RETURN
END SUBROUTINE DoLBLRTMLayers2KCARTALayers

!************************************************************************
! if TAPE5 has zeros everywhere for a certain gas, this reads in a US Std profile and stuffs it in

SUBROUTINE substitute_tape5_profile_for_climatology(iG,iNumLevs,raP,raaG_MR)


INTEGER, INTENT(IN OUT)                  :: iG
INTEGER, INTENT(IN)                      :: iNumLevs
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(OUT)                        :: raaG_MR(2*kProfLayer,kMaxGas)
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input


REAL :: ! pressure levels
! input/output


! local vars
REAL :: raX(2*kProfLayer)
INTEGER :: iL

CALL ReadRefProf_Levels2(iG,raP,iNumLevs,raX)

DO iL = 1,iNumLevs
  raaG_MR(iL,iG) = raX(iL)
END DO

RETURN
END SUBROUTINE substitute_tape5_profile_for_climatolog

!************************************************************************
! this reads record 1.2 of a LBLRTM TAPE5

SUBROUTINE read_record_1p2(caStr,  &
    IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS)


CHARACTER (LEN=80), INTENT(IN)           :: caStr
INTEGER, INTENT(IN OUT)                  :: IHIRAC
INTEGER, INTENT(IN OUT)                  :: ILBLF4
INTEGER, INTENT(IN OUT)                  :: ICNTNM
INTEGER, INTENT(IN OUT)                  :: IAERSL
INTEGER, INTENT(IN OUT)                  :: IEMIT
INTEGER, INTENT(IN OUT)                  :: ISCAN
INTEGER, INTENT(IN OUT)                  :: IFILTR
INTEGER, INTENT(IN OUT)                  :: IPLOT
INTEGER, INTENT(IN OUT)                  :: ITEST
INTEGER, INTENT(IN OUT)                  :: IATM
INTEGER, INTENT(IN OUT)                  :: IMRG
INTEGER, INTENT(IN OUT)                  :: ILAS
INTEGER, INTENT(IN OUT)                  :: IOD
NO TYPE, INTENT(IN OUT)                  :: IXSECT
NO TYPE, INTENT(IN OUT)                  :: MPTS
NO TYPE, INTENT(IN OUT)                  :: NPTS
IMPLICIT NONE

! input

! output


! local
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=2) :: c2
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=4) :: c4
CHARACTER (LEN=5) :: c5

! RECORD 1.2
!      IHIRAC, ILBLF4, ICNTNM, IAERSL,  IEMIT,  ISCAN, IFILTR, IPLOT, ITEST,  IATM,  IMRG,  ILAS,   IOD, IXSECT,  MPTS,  NPTS
!           5,     10,     15,     20,     25,     30,     35,    40,    45,    50, 54-55,    60,    65,     70, 72-75, 77-80
!       4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1,  4X,I1, 4X,I1, 4X,I1, 4X,I1, 3X,A2, 4X,I1, 4X,I1,  4X,I1, 1X,I4, 1X,I4
! this should read
! /home/sergio/IR_NIR_VIS_UV_RTcodes/LBLRTM/LBLRTM12.2/lblrtm/run_examples/run_example_user_defined_upwelling/TAPE5
!! reads HI=1 F4=1 CN=5 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0
!! reads HI=1 F4=1 CN=6 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=0 MG=1 LA=0 OD=1 XS=0    0    0

c1 = caStr(5:5)
READ(c1,'(I1)') IHIRAC
c1 = caStr(10:10)
READ(c1,'(I1)') ILBLF4
c1 = caStr(15:15)
READ(c1,'(I1)') ICNTNM
c1 = caStr(20:20)
READ(c1,'(I1)') IAERSL
c1 = caStr(25:25)
READ(c1,'(I1)') IEMIT
c1 = caStr(30:30)
READ(c1,'(I1)') ISCAN
c1 = caStr(35:35)
READ(c1,'(I1)') IFILTR
c1 = caStr(40:40)
READ(c1,'(I1)') IPLOT
c1 = caStr(45:45)
READ(c1,'(I1)') ITEST
c1 = caStr(50:50)
READ(c1,'(I1)') IATM
!      c2 = caStr(54:55)
!      read(c1,'(I2)') IMRG
!      read(c2,2) IMRG
c1 = caStr(55:55)
READ(c1,'(I1)') IMRG
c1 = caStr(60:60)
READ(c1,'(I1)') ILAS
c1 = caStr(65:65)
READ(c1,'(I1)') IOD
c1 = caStr(70:70)
READ(c1,'(I1)') IXSECT
c4 = caStr(72:75)
READ(c4,'(I4)') MPTS
c4 = caStr(77:80)
READ(c4,'(I4)') NPTS

2    FORMAT(I2)
4    FORMAT(I4)

!      print *,'record 1.2 : ',IHIRAC,ILBLF4,ICNTNM,IAERSL,IEMIT,ISCAN,IFILTR,IPLOT,ITEST,IATM,IMRG,ILAS,IOD,IXSECT,MPTS,NPTS

RETURN
END SUBROUTINE read_record_1p2

!************************************************************************
! this reads record 1.2a of a LBLRTM TAPE5

SUBROUTINE read_record_1p2a(caStr,XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL)


CHARACTER (LEN=80), INTENT(IN OUT)       :: caStr
INTEGER, INTENT(IN OUT)                  :: XSELF
INTEGER, INTENT(IN OUT)                  :: XFRGN
INTEGER, INTENT(IN OUT)                  :: XCO2C
INTEGER, INTENT(IN OUT)                  :: XO3CN
INTEGER, INTENT(IN OUT)                  :: XO2CN
INTEGER, INTENT(IN OUT)                  :: XN2CN
INTEGER, INTENT(IN OUT)                  :: XRAYL
IMPLICIT NONE

!!   XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, XRAYL in free format
!!     XSELF  H2O self broadened continuum absorption multiplicative factor
!!     XFRGN  H2O foreign broadened continuum absorption multiplicative factor
!!     XCO2C  CO2 continuum absorption multiplicative factor
!!     XO3CN  O3 continuum absorption multiplicative factor
!!     XO2CN  O2 continuum absorption multiplicative factor
!!     XN2CN  N2 continuum absorption multiplicative factor
!!      XRAYL Rayleigh extinction multiplicative factor

! input

! output


READ (caStr,*) XSELF,XFRGN,XCO2C,XO3CN,XO2CN,XN2CN,XRAYL

RETURN

END SUBROUTINE read_record_1p2a

!************************************************************************
! this reads record 1.3 of a LBLRTM TAPE5

SUBROUTINE read_record_1p3(caStr,rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT)


CHARACTER (LEN=80), INTENT(IN)           :: caStr
REAL, INTENT(IN OUT)                     :: rF1
REAL, INTENT(IN OUT)                     :: rF2
REAL, INTENT(IN OUT)                     :: rSample
REAL, INTENT(IN OUT)                     :: rDVset
REAL, INTENT(OUT)                        :: rALFAL0
REAL, INTENT(OUT)                        :: rAVMASS
REAL, INTENT(OUT)                        :: rDPTMIN
REAL, INTENT(OUT)                        :: rDPTFAC
INTEGER, INTENT(OUT)                     :: ILNFLG
REAL, INTENT(OUT)                        :: rDVOUT
IMPLICIT NONE

! input

! output



! local
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=10) :: c10

! RECORD 1.3    (required if IHIRAC > 0; IAERSL > 0; IEMIT = 1; IATM = 1; or ILAS > 0; otherwise omit)
!             V1,     V2,   SAMPLE,   DVSET,  ALFAL0,   AVMASS,   DPTMIN,   DPTFAC,   ILNFLG,     DVOUT
!           1-10,  11-20,    21-30,   31-40,   41-50,    51-60,    61-70,    71-80,     85,      90-100
!          E10.3,  E10.3,    E10.3,   E10.3,   E10.3,    E10.3,    E10.3,    E10.3,    4X,I1,  5X,E10.3

10   FORMAT(E10.3)

rSAMPLE = 4.0
rALFAL0 = 0.04
rAVMASS = 36
rDPTMIN = 0.0002
rDPTFAC = 0.001
ILNFLG = 0
rDVOUT = 0.0

c10 = caStr(01:10)
READ(c10,10) rF1
c10 = caStr(11:20)
READ(c10,10) rF2

c10 = caStr(21:30)
READ(c10,10) rSAMPLE
c10 = caStr(31:40)
READ(c10,10) rDVSET
c10 = caStr(41:50)
READ(c10,10) rALFAL0
c10 = caStr(51:60)
READ(c10,10) rAVMASS
c10 = caStr(61:70)
READ(c10,10) rDPTMIN
c10 = caStr(71:80)
READ(c10,10) rDPTFAC
!      c1 = caStr(85:85)
!      read(c1,'(I1)') ILFLG
!      c1 = caStr(91:100)
!      read(c10,10) rDVOUT

!      print *,'record 1.3 : ',rF1,rF2,rSample,rDVset,rALFAL0,rAVMASS,rDPTMIN,rDPTFAC,ILNFLG,rDVOUT

RETURN
END SUBROUTINE read_record_1p3

!************************************************************************
! this reads record 1.4 of a LBLRTM TAPE5

SUBROUTINE read_record_1p4(caStr,rTSurf,raEmiss,raRefl)


CHARACTER (LEN=80), INTENT(IN)           :: caStr
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(OUT)                        :: raEmiss(3)
REAL, INTENT(OUT)                        :: raRefl(3)
IMPLICIT NONE

! RECORD 1.4    (required if IEMIT = 1, or both IEMIT=2 and IOTFLG=2; otherwise omit)
!         TBOUND, SREMIS(1), SREMIS(2), SREMIS(3), SRREFL(1), SRREFL(2), SRREFL(3)
!           1-10,     11-20,     21-30,     31-40,     41-50,     51-60,     61-70
!          E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3,     E10.3

! input

! output


! local
CHARACTER (LEN=10) :: c10
INTEGER :: iI

10   FORMAT(E10.3)

c10 = caStr(01:10)
READ(c10,10) rTSurf
DO iI = 1,3
  raEmiss(iI) = -999.0
  c10 = caStr(iI*10+1:iI*10+10)
  READ(c10,10) raEmiss(iI)
END DO
IF (raEmiss(1) < 0.0) THEN
!! expecting emiss file
  raEmiss(2) = -1.0
  raEmiss(3) = -1.0
END IF

DO iI = 4,6
  raRefl(iI-3) = -999.0
  c10 = caStr(iI*10+1:iI*10+10)
  READ(c10,10) raRefl(iI-3)
END DO
IF (raRefl(1) < 0.0) THEN
!! expecting emiss file
  raRefl(2) = -1.0
  raRefl(3) = -1.0
END IF

!      print *,(raEmiss(iI),iI=1,3)
!      print *,(raRefl(iI),iI=1,3)
!      print *,'record 1.4 : tsurf = ',rTsurf

RETURN
END SUBROUTINE read_record_1p4

!************************************************************************
! this reads record 2.1 of a LBLRTM TAPE5

SUBROUTINE read_record_2p1(iIOUN2,iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN OUT)                  :: iForm
INTEGER, INTENT(IN OUT)                  :: iNumLevs
INTEGER, INTENT(IN OUT)                  :: iNumGases
REAL, INTENT(IN OUT)                     :: rSecnto
REAL, INTENT(IN OUT)                     :: rTopHgt
REAL, INTENT(IN OUT)                     :: rHSurf
REAL, INTENT(IN OUT)                     :: rViewAngle
IMPLICIT NONE

!         IFORM, NLAYRS, NMOL, SECNTO,        ZH1,       ZH2,    ZANGLE
!            2     3-5,   6-10,  11-20,      41-48,     53-60,     66-73
!          1X,I1    I3,    I5,   F10.2,  20X, F8.2,  4X, F8.2,  5X, F8.3

!              IFORM      (0,1) column amount format flag
!                           = 0  read PAVE, WKL(M,L), WBROADL(L) in F10.4, E10.3, E10.3 formats (default)
!                        = 1  read PAVE, WKL(M,L), WBROADL(L) in E15.7 format

!              NLAYRS      number of layers (maximum of 200)

!                NMOL      value of highest molecule number used (default = 7; maximum of 35)
!                                             See Table I for molecule numbers.
!              SECNTO      user entered scale factor for the column amount for the layers defined by NLAYRS
!                                             if positive, looking up
!                                              if negative, looking down
!                                                normal value = 1.0
!                 ZH1      observer altitude
!                 ZH2      end point altitude
!              ZANGLE      zenith angle at H1 (degrees)

! input

! output



! local
CHARACTER (LEN=80) :: caStr
CHARACTER (LEN=10) :: c10
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=5) :: c5
CHARACTER (LEN=8) :: c8

REAL :: rJunk
INTEGER :: iI

8    FORMAT(E8.3)
10   FORMAT(E10.3)
111  FORMAT(A80)

READ (iIOUN2,111,ERR=13,END=13) caStr

c1 = caStr(2:2)
READ(c1,*) iFORM
c3 = caStr(3:5)
READ(c3,*) iNumLevs
!! this is really number of layers, so need to increment by 1
!! iNumLevs = iNumLevs + 1
c5 = caStr(6:10)
READ(c5,*) iNumGases
c10 = caStr(11:20)
READ(c10,10) rSecnto   !! scale factor for layer ODS : positive if looking up, -ve IF looking down

c8 = caStr(41:48)
READ(c8,8) rTopHgt     !! observer altitude
c8 = caStr(53:60)
READ(c8,8) rHSurf      !! end point altitude
c8 = caStr(66:73)
READ(c8,8) rViewAngle  !! zenith angle at observer altitude (rTopHgt) so this is scanang

IF (rHSurf > rTopHgt) THEN
  rJunk = rTopHgt
  rTopHgt = rHSurf
  rHSurf = rJunk
END IF

!      print *,'record 2.1 : ',iForm,iNumLevs,iNumGases,rSecnto,rTopHgt,rHSurf,rViewAngle
13   CONTINUE

RETURN
END SUBROUTINE read_record_2p1

!************************************************************************
! this reads record 2.2 of a LBLRTM TAPE5

SUBROUTINE read_record_2p1p1(iIOUN2,iNumLevs,raP,raT,raPavg,raTavg,raaG_MR,raAlt,rPSurf,rHSurf,rTSurf,rPmin,rPmax,raSumCheck,  &
    rHminKCarta,rHmaxKCarta,rPminKCarta,rPmaxKCarta,  &
    iaG,iaGasUnits,iForm,iNumGases)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN OUT)                  :: iNumLevs
REAL, INTENT(OUT)                        :: raP(2*kProfLayer)
REAL, INTENT(OUT)                        :: raT(2*kProfLayer)
NO TYPE, INTENT(OUT)                     :: raPavg
REAL, INTENT(OUT)                        :: raTavg(2*kProfLayer)
REAL, INTENT(OUT)                        :: raaG_MR(2*kProfLayer,kMaxGas)
REAL, INTENT(OUT)                        :: raAlt(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(OUT)                        :: rHSurf
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(OUT)                        :: rPmin
REAL, INTENT(OUT)                        :: rPmax
NO TYPE, INTENT(OUT)                     :: raSumCheck
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iForm
INTEGER, INTENT(IN)                      :: iNumGases
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input

REAL :: rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta
! output
REAL :: raSumCh
REAL :: raPavg(2*kProfLay


! local
INTEGER :: iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2
CHARACTER (LEN=120) :: caStrX,caStrY
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=15) :: c15
CHARACTER (LEN=10) :: c10
CHARACTER (LEN=8) :: c8
CHARACTER (LEN=7) :: c7
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=2) :: c2
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=80) :: ca80
REAL :: rH,rP,rT,raX(kMaxGas),rJunk,rHSurfJunk

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
  IF (iForm == 0) THEN
    c10 = caStrY(01:10)
    READ(c10,10) rP
    c10 = caStrY(11:20)
    READ(c10,10) rT
    c10 = caStrY(21:30)
    READ(c10,10) rJunk
    c3 = caStrY(31:33)
    c2 = caStrY(34:35)
    ca80 = caStrY(37:120)
  ELSE IF (iForm == 1) THEN
    c15 = caStrY(01:15)
    READ(c15,15) rP
    c10 = caStrY(16:25)
    READ(c10,10) rT
    c10 = caStrY(26:35)
    READ(c10,10) rJunk
    c3 = caStrY(36:38)
    c2 = caStrY(39:40)
    ca80 = caStrY(41:120)
  END IF
  
  raPavg(iL) = rP
  raTavg(iL) = rT
  
  IF (iL. EQ. 1) THEN
!! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
!!                      which are basically
!!          SurfAlt      Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
    READ (ca80,*) rHSurfJunk,rPSurf,rTSurf,raAlt(iL+1),raP(iL+1),raT(iL+1)   !! p is in mb
    raAlt(iL) = rHSurfJunk
    raP(iL)   = rPSurf
    raT(iL)   = rTSurf
    IF (rPmax <= raP(iL)) rPmax = raP(iL)
    IF (rPmin > raP(iL)) rPmin = raP(iL)
    IF (rPmax <= rPSurf) rPmax = rPSurf
    IF (rPmin > rPSurf) rPmin = rPSurf
    
    IF (ABS(rHSurfJunk-rHSurf) > 0.001) THEN
      WRITE(kStdWarn,*) 'resetting rHSurf from first levl info in TAPE5 ',rHSurfJunk,rHSurf
      rHSurf = rHSurfJunk
    END IF
    
    IF (rPmax <= raP(iL+1)) rPmax = raP(iL+1)
    IF (rPmin > raP(iL+1)) rPmin = raP(iL+1)
    
    IF (raP(iL+1) > raP(iL)) THEN
      WRITE(kStdErr,*) 'huh reading LBLRTM TAPE5  iL,raP(iL),raP(iL+1) = ',iL,raP(iL),raP(iL+1)
      CALL DoStop
    END IF
    
!        print *,iL,raAlt(iL),raP(iL),raT(iL),rHSurfJunk,rPSurf,rTSurf
    
  ELSE
!! now read next "BLANK"  and A(z) P(z) and T(z)
!!                         which are basically
!!                         A(z) Pz)  and T(z) for next level, ad continuum
    READ (ca80,*) raAlt(iL+1),raP(iL+1),raT(iL+1)    !! p in mb
    IF (rPmax <= raP(iL+1)) rPmax = raP(iL+1)
    IF (rPmin > raP(iL+1)) rPmin = raP(iL+1)
!        print *,iL,raAlt(iL+1),raP(iL+1),raT(iL+1)
  END IF
  
!! now read MixRatio(z)
  IF (iNumGases <= 7) THEN
    READ (iIOUN2,120) caStrY
    READ (caStrY,*) (raX(iJ),iJ=1,iNumGases)
  ELSE IF (iNumGases > 7) THEN
    iGasCntLoop = 0
    100      CONTINUE
    iX1 = iGasCntLoop*7 + 1
    iX2 = iX1 + 7
    iX2 = MIN(iX2,iNumGases)
    READ (iIOUN2,120) caStrY
    READ (caStrY,*) (raX(iJ),iJ=iX1,iX2)
    IF (iX2 < iNumGases) THEN
      iGasCntLoop = iGasCntLoop + 1
      GO TO 100
    END IF
  END IF
  
  DO iJ = 1,iNumGases
    raaG_MR(iL+1,iJ) = raX(iJ)
    raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
  END DO
  IF (iL == 1) THEN
    DO iJ = 1,iNumGases
      raaG_MR(iL,iJ) = raX(iJ)
      raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
    END DO
  END IF
  
  IF (iL == 1) THEN
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
  
END DO       !!!! <<<<<<<<<<<<<<<<<<<< done reading the MOLGAS profiles, loop over levs

iNumLevs = iNumLevs+1   !!! since we added on surface level

kLBLRTM_toa = raP(iNumLevs)
WRITE(kStdWarn,*) 'kLBLRTM_toa = ',kLBLRTM_toa,' mb (when used with flux calcs)'
WRITE(kStdWarn,*) ' '

RETURN
END SUBROUTINE read_record_2p1p1

!************************************************************************
! read the xsect
! this reads record 2.2 of LBLRTM TAPE5

SUBROUTINE read_record_2p2(iIOUN2,iNumLevs,iNumGases,iNXsec,raP,raT,raaG_MR,raSumCheck,iaG,iaGasUnits)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN OUT)                  :: iNumLevs
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN)                      :: iNXsec
REAL, INTENT(IN OUT)                     :: raP(2*kProfLayer)
REAL, INTENT(IN OUT)                     :: raT(2*kProfLayer)
REAL, INTENT(OUT)                        :: raaG_MR(2*kProfLayer,kMaxGas)
REAL, INTENT(OUT)                        :: raSumCheck(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input


! output



! local
INTEGER :: iYes,iG,iJ,iL,iGasCntLoop,iRemain,iX1,iX2,iXSBIN,iFormX,iXmol,iNumLevsXsec
CHARACTER (LEN=120) :: caStr,caStrX,caStrY
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=15) :: c15
CHARACTER (LEN=10) :: c10
CHARACTER (LEN=8) :: c8
CHARACTER (LEN=7) :: c7
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=2) :: c2
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=80) :: ca80
REAL :: raPX(2*kProfLayer),raTX(2*kProfLayer),raAltX(2*kProfLayer),raX(kMaxGas),rP,rT,rJunk
REAL :: rTSurfX,rHSUrfX,rPSurfX

111  FORMAT(A80)
120  FORMAT(A120)
15   FORMAT(E15.7)
10   FORMAT(F10.4)
8    FORMAT(F8.3)
7    FORMAT(F7.2)

!! now see if there are xsec gases
READ (iIOUN2,111,ERR=13,END=13) caStr
READ(caStr,*) iNXsec,iXSBIN
IF (iNXsec > 0) THEN    !!!! >>>>>>>>>>>>>>> start reading XSCGAS profiles
  
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
  IF (iNumLevsXsec > 2*kProfLayer) THEN
    WRITE(kStdErr,*) 'iNumLevsXsec .GT. 2*kProfLayer',iNumLevsXsec,kProfLayer
    CALL DoStop
  END IF
  IF (iNumLevsXsec /= (iNumLevs-1)) THEN
    WRITE(kStdErr,*) 'iNumLevsXsec .NE. iNumLevs',iNumLevsXsec,iNumLevs-1
    CALL DoStop
  END IF
  
!! now ready to read in the profiles, same as in molgas above!
  DO iL = 1,iNumLevsXsec
    READ (iIOUN2,111) caStrY
    
!! first we want to read Pave and Tave
    IF (iFormX == 0) THEN
      c10 = caStrY(01:10)
      READ(c10,10) rP
      c10 = caStrY(11:20)
      READ(c10,10) rT
      c10 = caStrY(21:30)
      READ(c10,10) rJunk
      c3 = caStrY(31:33)
      c2 = caStrY(34:35)
      ca80 = caStrY(37:120)
    ELSE IF (iFormX == 1) THEN
      c15 = caStrY(01:15)
      READ(c15,15) rP
      c10 = caStrY(16:25)
      READ(c10,10) rT
      c10 = caStrY(26:35)
      READ(c10,10) rJunk
      c3 = caStrY(36:38)
      c2 = caStrY(39:40)
      ca80 = caStrY(41:120)
    END IF
    
    IF (iL. EQ. 1) THEN
!! now read A(z-1)    P(z-1)  T(z-1) and A(z) P(z) and T(z)
!!                      which are basically
!!          SurfAlt      Spres   Stemp      A(z) Pz)  and T(z) in the level above the ground
      READ (ca80,*) rHSurfX,rPSurfX,rTSurfX,raAltX(iL+1),raPX(iL+1),raTX(iL+1)   !! p is in mb
    ELSE
      READ (ca80,*) raAltX(iL+1),raPX(iL+1),raTX(iL+1)
    END IF
    
    IF (iNXsec <= 7) THEN
      READ (iIOUN2,120) caStrY
      READ (caStrY,*) (raX(iJ),iJ=1,iNXsec)
    ELSE IF (iNXsec > 7) THEN
      iGasCntLoop = 0
      100        CONTINUE
      iX1 = iGasCntLoop*7 + 1
      iX2 = iX1 + 7
      iX2 = MIN(iX2,iNXsec)
      READ (iIOUN2,120) caStrY
      READ (caStrY,*) (raX(iJ),iJ=iX1,iX2)
      IF (iX2 < iNXsec) THEN
        iGasCntLoop = iGasCntLoop + 1
        GO TO 100
      END IF
    END IF
    
    IF (iL == 1) THEN
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
END SUBROUTINE read_record_2p2
!************************************************************************
! this reads record 3.1 of LBLRTM TAPE5

SUBROUTINE read_record_3p1_and_3p2(iIOUN2,iNumGases,iBmax,rHSurf,rTopHgt,rSatZen)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iBmax
REAL, INTENT(IN OUT)                     :: rHSurf
REAL, INTENT(IN OUT)                     :: rTopHgt
REAL, INTENT(IN OUT)                     :: rSatZen
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input

! output



! local
INTEGER :: iJunk,iI,iY
CHARACTER (LEN=80) :: caStr
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=2) :: c2
CHARACTER (LEN=5) :: c5
CHARACTER (LEN=10) :: c10
REAL :: rCO2MIX,rJUNK

! 3.1
!      MODEL,  ITYPE, IBMAX,  NOZERO,  NOPRNT,  NMOL, IPUNCH, IFXTYP,   MUNITS,    RE, HSPACE,  VBAR, CO2MX
!          5,     10,    15,      20,      25,    30,     35,  36-37,    39-40, 41-50,  51-60, 61-70, 71-80
!         I5,     I5,    I5,      I5,      I5,    I5,     I5,     I2,   1X, I2, F10.3,  F10.3, F10.3, F10.3
READ(iIOUN2,111) caStr

c5 = caStr(1:5)
READ(c5,*) iJunk
IF (iJunk /= 0) THEN
  WRITE(kStdErr,*) 'in record 3.1 LBLRTM specifies one of its own internal profile OOPS'
  CALL dostop
END IF

c5 = caStr(6:10)
READ(c5,*) iJunk
IF (iJunk /= 3) THEN
  WRITE(kStdWarn,*) 'in record 3.1 LBLRTM specifies iType = ',iJunk,' OOPS only want 3 (gnd to space), resetting'
  iJunk = 3
END IF

c5 = caStr(11:15)
READ(c5,*) iBmax
IF (iBmax == 0) THEN
  WRITE(kStdWarn,*) 'in record 3.1 LBLRTM will use its own boundaries'
ELSE
  WRITE(kStdWarn,*) 'in record 3.1 LBLRTM will read in user specified boundaries'
END IF

c5 = caStr(16:20)
READ(c5,*) iJunk
c5 = caStr(21:25)
READ(c5,*) iJunk
c5 = caStr(26:30)
READ(c5,*) iNumGases
c5 = caStr(31:35)
READ(c5,*) iJunk

rCO2MIX = 330.0
c10 = caStr(71:80)
iY = -1
iI = 1
DO iI = 1,10
  IF (c10(iI:iI) /= ' ') iY = +1
END DO
IF  (iY > 0) READ(c10,*) rCO2mix

! 3.2
!         H1,    H2,   ANGLE,   RANGE,   BETA,   LEN,     HOBS
!       1-10, 11-20,   21-30,   31-40,  41-50, 51-55,    61-70
!      F10.3, F10.3,   F10.3,   F10.3,  F10.3,    I5, 5X,F10.3

READ(iIOUN2,111) caStr
c10 = caStr(1:10)
READ (c10,*) rHSurf
c10 = caStr(11:20)
READ (c10,*) rTopHgt
c10 = caStr(21:30)
READ (c10,*) rSatZen

IF (rHSurf > rTopHgt) THEN
  rJunk = rTopHgt
  rTopHgt = rHSurf
  rHSurf = rJunk
END IF

111  FORMAT(A80)

RETURN
END SUBROUTINE read_record_3p1_and_3p2
!************************************************************************
! this reads in record 3.3b of LBLRTM TAPE5
! altitudes of LBLRTM layer boundaries

SUBROUTINE read_record_3p3b(iIOUN2,iBmax,raZbnd)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN)                      :: iBmax
REAL, INTENT(IN OUT)                     :: raZbnd(2*kProfLayer)
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input

! output


! local
INTEGER :: iI,iCnt,i1,i2,iMax

iCnt = 1
20   CONTINUE
i1 = (iCnt-1)*8 + 1
i2 = iCnt*8
IF (i2 > iBmax) i2 = iBmax
READ(iIOUN2,*) (raZbnd(iI),iI=i1,i2)
IF (i2 < iBmax) THEN
  iCnt = iCnt + 1
  GO TO 20
END IF

RETURN
END SUBROUTINE read_record_3p3b

!************************************************************************
! this reads in record 3.4 of LBLRTM TAPE5
! user defined profile

SUBROUTINE read_record_3p4_and_3p5_and_3p7(iIOUN2,iNumGases,iNumLevs,rPSurf,rHSurf,rTSurf,  &
    rHmaxKCarta,rHminKCarta,rPmaxKCarta,rPminKCarta, raAlt,raP,raT,raSumCheck,  &
    iaG,iaGasUnits,raaG_MR,rPMin,rPMax,rYear,rLat,rLon, iZbnd,raZbnd,raPbnd)


INTEGER, INTENT(IN OUT)                  :: iIOUN2
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN)                      :: iNumLevs
REAL, INTENT(OUT)                        :: rPSurf
REAL, INTENT(OUT)                        :: rHSurf
REAL, INTENT(OUT)                        :: rTSurf
NO TYPE, INTENT(IN OUT)                  :: rHmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rHminKCart
NO TYPE, INTENT(IN OUT)                  :: rPmaxKCart
NO TYPE, INTENT(IN OUT)                  :: rPminKCart
REAL, INTENT(OUT)                        :: raAlt(2*kProfLayer)
REAL, INTENT(OUT)                        :: raP(2*kProfLayer)
REAL, INTENT(OUT)                        :: raT(2*kProfLayer)
REAL, INTENT(OUT)                        :: raSumCheck(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaG(kMaxGas)
INTEGER, INTENT(OUT)                     :: iaGasUnits(kMaxGas)
REAL, INTENT(OUT)                        :: raaG_MR(2*kProfLayer,kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: rPMin
NO TYPE, INTENT(IN OUT)                  :: rPMax
REAL, INTENT(IN OUT)                     :: rYear
REAL, INTENT(IN OUT)                     :: rLat
REAL, INTENT(IN OUT)                     :: rLon
NO TYPE, INTENT(IN)                      :: iZbnd
REAL, INTENT(IN OUT)                     :: raZbnd(2*kProfLayer)
REAL, INTENT(OUT)                        :: raPbnd(2*kProfLayer)
IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

! input

! output


REAL :: rPminKCarta,rPmaxKCarta,rHminKCarta,rHmaxKCarta

REAL :: rPmin,rPmax


INTEGER :: iZbnd   !!! are we using default 101 levels, or LBLRTM defined?

! local var
CHARACTER (LEN=5) :: c5
INTEGER :: iErr,iErrIO,iL,iJ,iG,iMid,ifloor,iaJunk(20),iNumLevsXsec,iNXsec,iLBROutBdryHorP
REAL :: raX(kMaxGas),rX,rP,rT,rF1,rF2,rTophgt,rViewAngle,raLBL_Hgts(kProfLayer),rH
CHARACTER (LEN=80) :: caStr,caStrX,caStrY
CHARACTER (LEN=30) :: caStr30
CHARACTER (LEN=1) :: c1

DO iJ = 1,iNumGases
  raSumCheck(iJ) = 0.0
END DO

rPmin = +1.0E6
rPmax = -1.0E+6
iNXsec = -1

! record 3.4
111  FORMAT(A80)
READ(iIOUN2,111) caStr
c5 = caStr(1:5)
READ(c5,*) iNumLevs

! record 3.5
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
  IF (iL == 1) THEN
    rPSurf = rP
    rHSurf = rH
    IF (ABS(rTSurf-rT) >= 0.001) THEN
!! hmm looks like info in first item of Record 3.5 is inconsistent with info in Record 1.4
      WRITE(kStdWarn,*) 'rT first item of Record 3.5 is inconsistent with TBound info in Record 1.4'
      WRITE(kSTdWarn,*) rTSurf,rT
      WRITE(kStdWarn,*) 'artificially set rPSurf to be a little higher than lowest "p" entry'
      rPSurf = rPSurf-0.125
    END IF
!rTSurf = rT
  END IF
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
  
!        raP(iL) = rP * 100.0  !! change from mb to N/m2
  raP(iL) = rP          !! keep in mb
  raT(iL) = rT
  raAlt(iL) = rH
  IF (rPmax <= raP(iL)) rPmax = raP(iL)
  IF (rPmin > raP(iL)) rPmin = raP(iL)
  READ (iIOUN2,*) (raX(iG),iG=1,iNumGases)
  DO iJ = 1,iNumGases
    raaG_MR(iL,iJ) = raX(iJ)
    raSumCheck(iJ) = raSumCheck(iJ) + raX(iJ)
  END DO
  
  IF (iL == 1) THEN
!! rPSurf already set a few lines above, and need to be careful since we want bdry levels set
!! so do NOT play with it here
!          rPSurf = raP(1)/100 !! because we redo this below
!          rPSurf = raP(1)     !! keep in mb
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
      raSumCheck(iJ+iNumGases) = raSumCheck(iJ+iNumGases) + raX(iJ)
    END DO
  END DO
  iNumGases = iNumGases + iNXsec
END IF

IF (iZbnd > 0) THEN
!!need to convert the user heights to user defined pressure levels
  CALL rspl(raAlt,raP,iNumLevs,raZbnd,raPbnd,iZbnd)
  DO iJ = 1,iZbnd
    raPbnd(iJ) = raPbnd(iJ) * 100.0  !! change from mb to N/m2, as Psurf is finally in N/m2
  END DO
END IF

13   CONTINUE
RETURN
END SUBROUTINE read_record_3p4_and_3p5_and_3p7
!************************************************************************
