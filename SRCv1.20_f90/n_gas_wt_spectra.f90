! Copyright 2000
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:42
 
! University of Maryland Baltimore County
! All Rights Reserved

! this file reads *MOLGAS,*XSCGAS,*FREQ (very easy),*MIXFIL

!************************************************************************
! this subroutine deals with the 'MOLGAS' keyword
! and basically deals with GasID 1--40

! main output parameter is iaMOLgases
! synopsis IN : iNGas     = number of gasIDs to be read in (if -1, all)
!               iaGasesNL = list of gasIDs
!         OUT : iaMOLgases(iI) = set to 1 if gasID was found in this section

SUBROUTINE molgas4(iNGas,iaGasesNL,iaMOLgases)


INTEGER, INTENT(IN OUT)                  :: iNGas
INTEGER, INTENT(IN)                      :: iaGasesNL(kGasComp)
INTEGER, INTENT(OUT)                     :: iaMOLgases(kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iNGas     = number of gases stored in this namelist
! iaGasesNL = gasIDs stored in this namelist
!             the elements contain the gasIDs, or -1 ow
!                so eg if gases1,2,5,101 needed then
!            iaGasesNL(1)=1  iaGasesNL(2)=2  iaGasesNL(3)=5 iaGasesNL(4)=101
! iaMOLgases   = overall gasIDs to be used by kCARTA
!                the elements are set to 0 if gas not required, 1 if required
!                so eg if gases1,2,5,101 needed then
!            iaMOLgases(1)=iaMOLgases(2)=iaMOLgases(5)=iaMOLgases(101)=1
!            all else 0
! MOLGAS or XSCGAS ... array iaInputOrder is built up from this array, each
! time it is returned from a call to either subroutine MOLGAS or xscfil


! local variables
INTEGER :: iTag,iErr,iNotMainGas,MainGas
INTEGER :: iInt,iaTemp(kMaxGas),iC,iCC,iaInDataBase(kGasComp)
INTEGER :: iNgasesCheck
INTEGER :: iCheckCompDataBase

CHARACTER (LEN=7) :: caWord

caWord = '*MOLGAS'

DO iInt = 1,kMaxGas
  iaMOLgases(iInt) = 0
  iaTemp(iInt) = 0
END DO

! allow iNgas = -5/-6 ,-1,or > 0 where
! -114  means allow LBL gases which were in older kcarta versions (upto v1.14) REMOVED
!       which was 1..32,101,102,103
! -5/-6 means allow LBLRTM style gases 1-32
! -1    means allow ALL LBL gases from the new HITRAN208 = 1..42
! > 0   means you choose the gases (using iaGasesNL)

! read the no of gases stored in the file
iErr = -1
IF (iNgas > 0) THEN
  WRITE(kStdWarn,*) 'including user specified lbl gases'
! use the molecular ID's from the namelist file
  DO iInt = 1,iNgas
    iNotMainGas = MainGas(iaGasesNL(iInt))
    IF (iNotMainGas < 0) THEN
      WRITE(kStdErr,*) 'Invalid MOLGAS GasID',iaGasesNL(iInt),' entered'
      WRITE(kStdErr,*) 'Please check *MOLGAS and retry'
      WRITE(kStdErr,*) 'Note : If the GasID printed was -100, you '
      WRITE(kStdErr,*) 'probably entered less gasIDs than you promised '
      CALL DoSTOP
    ELSE
      iaTemp(iInt) = iaGasesNL(iInt)
    END IF
  END DO
  CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)   !! user specified gases
  
ELSE IF (iNgas == -1) THEN
  CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)   !! latest version, uses 42 molgas
  
ELSE IF ((iNGas == -5) .OR. (iNGas == -6)) THEN  !! use 36 gases for LBLRTM
  CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)
  
ELSE
  WRITE(kStdErr,*) 'input user iNXsec = ',iNgas
  WRITE(kStdErr,*) '  valid iNgas : -5/-6 (lblrtm) -1 (all) +N select few'
  CALL DoStop
END IF

RETURN
END SUBROUTINE molgas4

!************************************************************************
! this subroutine deals with the 'XSCGAS' keyword
! and basically deals with GasID 51--63
! synopsis IN : iNXsec    = number of gasIDs to be read in (if -1, all)
!               iaLXsecNL = list of gasIDs
!         OUT : iaXSCgases(iI) = set to 1 if gasID was found in this section

SUBROUTINE xscgas4(iNXsec,iaLXsecNL,iaXSCgases)


INTEGER, INTENT(IN OUT)                  :: iNXsec
INTEGER, INTENT(IN)                      :: iaLXsecNL(kGasXSecHi-kGasXSecLo+1
INTEGER, INTENT(OUT)                     :: iaXSCgases(kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iNXSec     = number of gases stored in this namelist
! iaLXGasesNL = gasIDs stored in this namelist
!             the elements contain the gasIDs, or -1 ow
!                so eg if gases 51,52 needed then
!            iaLXSecNL(1)=51  iaLXSecNL(2)=52
! iaXSCgases   = overall gasIDs to be used by kCARTA
!                the elements are set to 0 if gas not required, 1 if required
!                so eg if gases 51,52 needed then
!            iaXSCgases(51)=iaXSCgases(52)=1
!            all else 0
! MOLGAS or XSCGAS ... array iaInputOrder is built up from this array, each
! time it is returned from a call to either subroutine MOLGAS or xscfil


! local variables
CHARACTER (LEN=7) :: caWord
INTEGER :: iInt,iErr,iC,iCC,iCheckXsecDataBase,iNXsecCheck
INTEGER :: iaTemp(kMaxGas),iTag
INTEGER :: iaInDataBase(kMaxLayer)

caWord = '*XSCGAS'

DO iInt = 1,kMaxGas
  iaXSCgases(iInt) = 0
  iaTemp(iInt) = 0
END DO

! allow iNxsec = -5/-6 ,-1,or > 0 where
! -114  means allow LBL gases which were in older kcarta versions (upto v1.14)
!       which was 51..63 REMOVED
! -5/-6 means allow LBLRTM gases 51..63
! -1    means allow ALL XSC gases from the new HITRAN208 = 51..81
! > 0   means you choose the gases (using iaLXsecNL)

iErr = -1
IF (iNxsec > 0) THEN
  WRITE(kStdWarn,*) 'including user specified xsc gases'
! use the xsec gas ID's from the namelist file
  DO iInt = 1,iNXsec
    IF ((iaLXsecNL(iInt) < kGasXsecLo) .OR.  &
          (iaLXsecNL(iInt) > kGasXsecHi)) THEN
      WRITE(kStdErr,*) 'Invalid XSCGAS GasID',iaLXsecNL(iInt),' entered'
      WRITE(kStdErr,*) 'Please check *XSCGAS and retry'
      WRITE(kStdErr,*) 'Note : If the GasID printed was -100, you '
      WRITE(kStdErr,*) 'probably entered less gasIDs than you promised '
      CALL DoSTOP
    ELSE
      iaTemp(iInt) = iaLXsecNL(iInt)
    END IF
  END DO
  CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! now check the veracity of these gases
ELSE IF (iNxsec == -1) THEN
  CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! add all xsec gases
ELSE IF ((iNxsec == -5) .OR. (iNXsec == -6)) THEN
  CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! add LBLRTM xsec gases
ELSE
  WRITE(kStdErr,*) 'input user iNXsec = ',iNXsec
  WRITE(kStdErr,*) '  valid iNXsec : -5/-6 (lblrtm) -1 (all) +N select few'
  CALL DoStop
END IF

RETURN
END SUBROUTINE xscgas4

!************************************************************************
! this subroutine is for the new spectra
! this reads in the info for the new spectra for ONE gas

SUBROUTINE spectra4(iNumNewGases,  &
    iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks)


NO TYPE, INTENT(IN OUT)                  :: iNumNewGas
INTEGER, INTENT(IN)                      :: iaNewGasID(kGasStore)
INTEGER, INTENT(IN OUT)                  :: iaNewData(kGasStore)
NO TYPE, INTENT(IN OUT)                  :: iaaNewChun
NO TYPE, INTENT(IN OUT)                  :: caaaNewChu
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks

INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaaNewChunks(kGasStore,kNumkCompT)

! local variables
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead,iCount,iErr

caWord = '*SPECTR'
iErr = -1

5030 FORMAT(A130)

iNumLinesRead = 0
13   IF (iNumLinesRead > 0) THEN
  iErr = 1
  WRITE(kStdErr,5010) caWord
  CALL DoSTOP
END IF
5010 FORMAT('Error reading section ',A7)

iNumLinesRead = 1
! read how many gases will have new spectroscopy
iErr = -1
IF ((iNumNewGases < 1) .OR. (iNumNewGases > kGasStore)) THEN
  WRITE(kStdErr,*)'need a valid number of gases in *SPECTRA!!'
  WRITE(kStdErr,*)'must be > 0, < kGasStore!!'
  WRITE(kStdErr,*)'please check and retry!'
  CALL DoSTOP
END IF

DO iCount = 1,iNumNewGases
  IF ((iaNewGasID(iCount) < 1) .OR. (iaNewGasID(iCount) > kMaxGas)) THEN
    WRITE(kStdErr,*) 'iaNewGasID(',iCount,') = ',iaNewGasID(iCount)
    WRITE(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry (subr spectra4)!'
    CALL DoSTOP
  END IF
  IF ((iaNewData(iCount) < 1) .OR. (iaNewData(iCount) > kNumkCompT)) THEN
    WRITE(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry!'
    CALL DoSTOP
  END IF
END DO

! really no more error checking possible, as we have to actually go ahead and
! open the files as necessary

RETURN
END SUBROUTINE spectra4

!************************************************************************
! this subroutine is for nonlte ... does the decision between SLOW or FAST
! this reads in the info for the new spectra for ONE gas

SUBROUTINE nonlte( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
    iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
    caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)


INTEGER, INTENT(IN OUT)                  :: iRTP
INTEGER, INTENT(IN OUT)                  :: iaGasesIn(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNumGases
REAL, INTENT(IN OUT)                     :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raNLTEstre
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iNLTE_Slow
NO TYPE, INTENT(IN OUT)                  :: iaNLTEGasI
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
NO TYPE, INTENT(IN OUT)                  :: caaStrongL
NO TYPE, INTENT(IN OUT)                  :: iaNLTEBand
NO TYPE, INTENT(IN OUT)                  :: raNLTEstar
NO TYPE, INTENT(IN OUT)                  :: caaaNLTEBa
NO TYPE, INTENT(IN OUT)                  :: caaNLTETem
NO TYPE, INTENT(IN OUT)                  :: caaUpperMi
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iAllLayers
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
NO TYPE, INTENT(IN OUT)                  :: raKSolarAn
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iSetBloat
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
NO TYPE, INTENT(IN OUT)                  :: caOutUABlo
NO TYPE, INTENT(IN OUT)                  :: caPlanckUA
NO TYPE, INTENT(IN OUT)                  :: caOutUAfil
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (-1/-2) model
! iSetBloat        tells us if stick to 0.0025cm-1 or bloat up to 0.0005 cm-1
! raNLTEstrength   tells how strongly to add on all the files (default 1.0)
! iNumNLTEGases    tells number of NLTE gases
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines   line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE  tells if the code does upper atm NLTE
! iAllLayersLTE    tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd  tells the code if use weak background lines as well, or not
INTEGER :: iNLTE_SlowORFast

INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
REAL :: raNLTEstrength(kGasStore)
INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
INTEGER :: iaNLTEBands(kGasStore)
REAL :: raNLTEstart(kGasStore)
CHARACTER (LEN=80) :: caaaNLTEBands(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaNLTETemp(kGasStore)

! the solar angle
REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
CHARACTER (LEN=80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
CHARACTER (LEN=80) :: caPlanckUAfile,caOutUAfile

IF (iNLTE_SlowORFast == +1) THEN
  CALL nonlteSLOW_LBL( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
      raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
      iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
      caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
      raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
      iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
      raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
      caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)
  
ELSE IF (iNLTE_SlowORFast == -1) THEN
  CALL nonlteFAST_SARTA( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
      raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
      iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
      caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
      raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
      iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
      raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
      caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)
  
ELSE IF (iNLTE_SlowORFast == -2) THEN
  CALL nonlteFAST_KCOMP( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
      raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
      iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
      caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
      raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
      iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
      raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
      caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)
END IF

RETURN
END SUBROUTINE nonlte

!************************************************************************
! this subroutine is for FAST version of nonlte
! this reads in the info for the new spectra for ONE gas

SUBROUTINE nonlteFAST_SARTA( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
    iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
    caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)


INTEGER, INTENT(IN OUT)                  :: iRTP
INTEGER, INTENT(IN OUT)                  :: iaGasesIn(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNumGases
REAL, INTENT(IN OUT)                     :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raNLTEstre
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iNLTE_Slow
NO TYPE, INTENT(IN OUT)                  :: iaNLTEGasI
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
NO TYPE, INTENT(IN OUT)                  :: caaStrongL
NO TYPE, INTENT(IN OUT)                  :: iaNLTEBand
NO TYPE, INTENT(IN OUT)                  :: raNLTEstar
NO TYPE, INTENT(IN OUT)                  :: caaaNLTEBa
NO TYPE, INTENT(IN OUT)                  :: caaNLTETem
NO TYPE, INTENT(IN OUT)                  :: caaUpperMi
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iAllLayers
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
NO TYPE, INTENT(IN OUT)                  :: raKSolarAn
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(OUT)                     :: iSetBloat
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
NO TYPE, INTENT(IN OUT)                  :: caOutUABlo
NO TYPE, INTENT(IN OUT)                  :: caPlanckUA
NO TYPE, INTENT(IN OUT)                  :: caOutUAfil
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (+1) model
! iSetBloat        tells us if stick to 0.0025cm-1 or bloat up to 0.0005 cm-1
! raNLTEstrength   tells how strongly to add on all the files (default 1.0)
! iNumNLTEGases    tells number of NLTE gases
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines   line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE  tells if the code does upper atm NLTE
! iAllLayersLTE    tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd  tells the code if use weak background lines as well, or not
INTEGER :: iNLTE_SlowORFast

INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
REAL :: raNLTEstrength(kGasStore)
INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
INTEGER :: iaNLTEBands(kGasStore)
REAL :: raNLTEstart(kGasStore)
CHARACTER (LEN=80) :: caaaNLTEBands(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaNLTETemp(kGasStore)

! the solar angle
REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
CHARACTER (LEN=80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
CHARACTER (LEN=80) :: caPlanckUAfile,caOutUAfile

! local vars
INTEGER :: iI,iJ

caPlanckBloatFile = 'blankcaPlanckBloatFile'
caOutBloatFile    = 'blankcaOutBloatFile'
caOutUABloatFile  = 'blankcaOutUABloatFile'
caPlanckUAfile    = 'blankcaPlanckUAfile'
caOutUAfile       = 'blankcaOutUAfile'

IF (iSetBloat > 0) THEN
  iSetBloat = -1
  WRITE(kStdWarn,*) 'Reset iSetBloat to -1 in nonlteFAST_SARTA'
END IF
IF (iDoUpperAtmNLTE > 0) THEN
  iDoUpperAtmNLTE = -1
  WRITE(kStdWarn,*) 'Reset iDoUpperAtmNLTE to -1 in nonlteFAST_SARTA'
END IF

!      DO iI = 1,iNumNLTEGases
!        print *,iaNLTEGasID(iI),raNLTEstrength(iI),raNLTEstart(iI),
!     $           iaNLTEChunks(iI)
!        DO iJ = 1,iaNLTEChunks(iI)
!          print *,iJ,iaaNLTEChunks(iI,iJ)
!          END DO
!        END DO

RETURN
END SUBROUTINE nonlteFAST_SARTA

!************************************************************************
! this subroutine is for GENLN2_LBL nonlte
! this reads in the info for the new spectra for ONE gas

SUBROUTINE nonlteSLOW_LBL( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
    iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
    caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)


INTEGER, INTENT(IN OUT)                  :: iRTP
INTEGER, INTENT(IN)                      :: iaGasesIn(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumGases
REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raNLTEstre
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iNLTE_Slow
NO TYPE, INTENT(IN OUT)                  :: iaNLTEGasI
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
NO TYPE, INTENT(IN OUT)                  :: caaStrongL
NO TYPE, INTENT(IN OUT)                  :: iaNLTEBand
NO TYPE, INTENT(IN OUT)                  :: raNLTEstar
NO TYPE, INTENT(IN OUT)                  :: caaaNLTEBa
NO TYPE, INTENT(IN OUT)                  :: caaNLTETem
NO TYPE, INTENT(IN OUT)                  :: caaUpperMi
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iAllLayers
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
NO TYPE, INTENT(IN OUT)                  :: raKSolarAn
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(OUT)                     :: iSetBloat
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
NO TYPE, INTENT(IN OUT)                  :: caOutUABlo
NO TYPE, INTENT(IN OUT)                  :: caPlanckUA
NO TYPE, INTENT(IN OUT)                  :: caOutUAfil
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (-1/-2) model
! iSetBloat        tells us if stick to 0.0025cm-1 or bloat up to 0.0005 cm-1
! raNLTEstrength   tells how strongly to add on all the files (default 1.0)
! iNumNLTEGases    tells number of NLTE gases
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines   line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE  tells if the code does upper atm NLTE
! iAllLayersLTE    tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd  tells the code if use weak background lines as well, or not
INTEGER :: iNLTE_SlowORFast

INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
REAL :: raNLTEstrength(kGasStore)
INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
INTEGER :: iaNLTEBands(kGasStore)
REAL :: raNLTEstart(kGasStore)
CHARACTER (LEN=80) :: caaaNLTEBands(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaNLTETemp(kGasStore)

! the solar angle
REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
CHARACTER (LEN=80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
CHARACTER (LEN=80) :: caPlanckUAfile,caOutUAfile

! local variables
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead,iErr,iI,iaDumb(kMaxGas),iaGases(kMaxGas)
REAL :: rH,rHN
INTEGER :: iStrongISO,iStrongJL,iStrongJU,iStrongGASID,iLTEIn,OutsideSpectra
INTEGER :: iNLTEStart2350,iJ,iG

INTEGER :: iBand,iGasID,iNum,iISO,iLineMixBand,iInt,iType,iGas,iDoVoigtChi
DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
CHARACTER (LEN=1) :: caJPQR(kHITRAN)

CHARACTER (LEN=80) :: caNONLTETempKC
INTEGER :: iRegr,iMethod,iDefault
REAL :: raTemp(kProfLayer)
CHARACTER (LEN=80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT)
INTEGER :: iJunkNum,iaJunk(kGasStore)

CALL NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands,caaaNLTEBands)

IF (ABS(iSetBloat) /= 1) THEN
  WRITE(kStdErr,*) 'need iSetBloat  = -1,+1',iSetBloat
  CALL DOSTOP
END IF

IF (ABS(iDoUpperAtmNLTE) /= 1) THEN
  WRITE(kStdErr,*) 'need iDoUpperAtmNLTE  = -1,+1',iDoUpperAtmNLTE
  CALL DOSTOP
END IF

!!! how to predict NLTE temperatures before 2004
iMethod = 1      !!! lapack  .. not supported anymore
iMethod = 3      !!! predictor
iMethod = 2      !!! polynom
!!after 2004, just use GRANADA model EXTERNALLY to kCARTA

iDefault = 2
IF (iMethod /= iDefault) THEN
  WRITE(kStdErr,*) 'imethod,idefault in subr nonlte = ',iMethod,iDefault
END IF
IF (iMethod == 3) THEN
  WRITE(kStdErr,*) 'imethod = 3 still not programmed up'
  CALL DoStop
END IF
IF (iMethod == 1) THEN
  WRITE(kStdErr,*) 'imethod = 1 LAPACK has been voided'
  CALL DoStop
END IF

caNONLTETempKC = 'boo_nlte_file'
DO iI = 1,kProfLayer
  raTemp(iI) = raaTemp(iI,1)   !!assume all kinetic temps equal
END DO

IF (iNumNLTEGases > 1) THEN
  WRITE(kStdErr,*) 'kCARTA can only handle ONE NLTE gas ... CO2!!!'
  CALL DoStop
END IF

IF (iNumNLTEGases < 1) THEN
  WRITE(kStdErr,*) 'huh, 0 or -ve number of NLTE gases????'
  CALL DoStop
END IF

iRegr = -1
DO iI = 1,iNumNLTEGases
  IF (caaNLTETemp(iI) == 'nlteguess') THEN
!          IF (iMethod .EQ. +1) THEN
!            !!figure out iRegr and reset it
!            CALL closest_regr_lapack(
!     $         caOutName,raTemp,raPressLevels,iProfileLayers,raKSolarAngle(1),
!     $         iRTP,iRegr,caNONLTETempKC)
    IF (iMethod == +2) THEN
!!leave iRegr as -1
      CALL polynom_nlte(  &
          caOutName,raTemp,raPressLevels,iProfileLayers,raKSolarAngle(1),  &
          iRTP,iRegr,caNONLTETempKC)
!          ELSEIF (iMethod .EQ. +3) THEN
!            !!leave iRegr as -1
!            CALL predict_nlte(
!     $         caOutName,raTemp,raPressLevels,iProfileLayers,raKSolarAngle(1),
!     $         iRTP,iRegr,caNONLTETempKC)
    END IF
    caaNLTETemp(iI) = caNONLTETempKC
  END IF
END DO

caWord = '*NONLTE'
iErr = -1

! we are calling this routine a little early in n_main ... need to reassign
! iaGases correctly
! upto this point eg if gas IDs were 1,3,5,22,51 then
! iaGasesIn(1) = iaGasesIn(3) = iaGasesIn(5) = iaGasesIn(22) = iaGasesIn(51)  =  1
! all other iaGasesIn(iN) = -1
! we have to redo this so that iaGases contains the LIST
! iaGases(1,2,3,4,5) = 1,3,5,22,51  all else -1
DO iInt = 1,kMaxGas
  iaDumb(iInt) = iaGasesIn(iInt)
  iaGases(iInt) = -1
END DO

iType = 1
DO iInt = 1,kMaxGas
  IF (iaDumb(iInt) > 0) THEN
    iaGases(iType)  =  iInt
    iType = iType + 1
  END IF
END DO

5030 FORMAT(A130)

iNumLinesRead = 0
13   IF (iNumLinesRead > 0) THEN
  iErr = 1
  WRITE(kStdErr,5010) caWord
  CALL DoSTOP
END IF
5010 FORMAT('Error reading section ',A7)

DO iI = 1,kGasStore
  IF (ABS(raNLTEstrength(iI)) > 10.0) THEN
    WRITE(kStdWarn,*)'gas',iI, 'multiply NLTE data by',raNLTEstrength(iI)
  END IF
END DO

iNumLinesRead = 1
! read how many gases will have new nonLTE spectroscopy
iErr = -1
IF ((iNumNLTEGases < 1) .OR. (iNumNLTEGases > kGasStore)) THEN
  WRITE(kStdErr,*)'need a valid number of gases in *NONLTE!!'
  WRITE(kStdErr,*)'must be > 0, < kGasStore!!'
  WRITE(kStdErr,*)'please check and retry!'
  CALL DoSTOP
END IF

DO iLTEIn = 1,iNumNLTEGases
  IF ((iaNLTEGasID(iLTEIn) < 1) .OR. (iaNLTEGasID(iLTEIn) > kMaxGas)) THEN
    WRITE(kStdErr,*)'iaNLTEGasID(',iLTEIn,') = ',iaNLTEGasID(iLTEIn)
    WRITE(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry (subr nonlteSLOW_LBL)!'
    CALL DoSTOP
  END IF
  IF ((iaNLTEChunks(iLTEIn) < 1) .OR.  &
        (iaNLTEChunks(iLTEIn) > kNumkCompT)) THEN
    WRITE(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry!'
    CALL DoSTOP
  END IF
END DO

! now figure out above which height the NLTE starts
! iaNLTEStart is for MOST of the bands (ie the weaker bands)
! ********************** this is using the user supplied info *************
DO iLTEIn = 1,iNumNLTEGases
  rH = raNLTEstart(iLTEIn)*1000     !!!NLTE start height in m
  IF (rH <= 0) THEN
!!! easy, use all layers for the NLTE
    iI = kProfLayer - iProfileLayers + 1
    iaNLTEStart(iLTEIn) = iI
    GO TO 24
  ELSE IF (rH > 0) THEN
!!!go ahead and find iStart
    iI = kProfLayer - iProfileLayers + 1
    iI = iI - 1
    iI = MAX(1,iI)
    23       CONTINUE
    IF ((rH >= raLayerHeight(iI)) .AND. (iI. LT. kProfLayer)) THEN
      iI = iI + 1
      GO TO 23
    END IF
    IF ((rH >= raLayerHeight(iI)) .AND. (iI. EQ. kProfLayer)) THEN
      WRITE(kStdErr,*) 'You are starting NLTE too high!!! '
      WRITE(kStdErr,*) 'Max layer height (m) = ',raLayerHeight(iI)
      WRITE(kStdErr,*) 'Your NLTE start (m)  = ',rH
      CALL DoStop
    ELSE
      iaNLTEStart(iLTEIn) = MAX(iI,1)
      WRITE(kStdWarn,*) '(nml file) user defined NLTE starts at lay ',iI
      WRITE(kStdWarn,*) 'for nlte gas number ',iLTEIn,' = gasID ',iaNLTEGasID(iLTEIn)
      WRITE(kStdWarn,*) 'which is layer at ~ ',rH/1000,' km'
    END IF
  END IF
  24     CONTINUE
END DO

! now figure out above which height the NLTE starts, for the strong bands
! iaNLTEStart2350 is for STRONGEST NLTE bands, in 4 um region ahem!
! ************* this is using the D. Edwards/M. Lopez-Puertas info ************
iJunkNum    = -1
rH = 1.0E10             !!!dumb large number, too high!!!!!!
DO iGas = 1,iNumGases
  iLTEIn = OutsideSpectra(iaGases(iGas),iNumNLTEGases,iaNLTEGasID,iJunkNum,iaJunk,2205.0,605.0,2830.0,20)
  IF (iLTEIn > 0) THEN
    rHN = raNLTEstart(iLTEIn)*1000     !!!NLTE start height in m
    rH = MIN(rH,rHN)
!!! now see if we can come up with a better estimate
    iStrongGASID = iaGases(iGas)
    DO iBand = 1,iaNLTEBands(iLTEIn)
!!read in the lineshape parameters for the band
      CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands,  &
          iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift,  &
          daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,  &
          iLineMixBand,iDoVoigtChi)
      iStrongISO   = iISO
      iStrongJU    = nint(daJU(1))
      iStrongJL    = nint(daJL(1))
      CALL FindNLTEHeight(iLTEIn,caaNLTETemp,  &
          iStrongISO,iStrongJL,iStrongJU,iStrongGasID,  &
          iaGases,raaTemp,raaPress,  &
          raPressLevels,raLayerHeight,iProfileLayers,rHN,iBand)
      rH = MIN(rH,rHN)
    END DO
    
!!!now process the results
    IF (rH < 0) THEN
!!! easy, use all layers for the NLTE
      iI = kProfLayer - iProfileLayers + 1
      iNLTEStart2350 = iI
    ELSE
!!!go ahead and find iStart
      iI = kProfLayer - iProfileLayers + 1
      iI = iI - 1
      iI = MAX(1,iI)
      26         CONTINUE
      IF ((rH >= raLayerHeight(iI)) .AND. (iI. LT. kProfLayer)) THEN
        iI = iI + 1
        GO TO 26
      END IF
      IF ((rH >= raLayerHeight(iI)) .AND. (iI. EQ. kProfLayer)) THEN
        WRITE(kStdErr,*) 'You are starting NLTE too high!!! '
        WRITE(kStdErr,*) 'Max layer height (m) = ',raLayerHeight(iI)
        WRITE(kStdErr,*) 'Your NLTE start (m)  = ',rH
        CALL DoStop
      ELSE
        iNLTEStart2350 = iI
        WRITE(kStdWarn,*) 'From reading the NLTE profiles '
        WRITE(kStdWarn,*) 'Strong NLTE starts ',iI,' for gas ',iLTEIn
      END IF
    END IF
  END IF
END DO

WRITE (kStdWarn,*) '   '
WRITE (kStdWarn,*) '   '
WRITE (kStdWarn,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>> '
DO iI = 1,iNumNLTEGases
  iaNLTEStart2350(iI) = MIN(iNLTEStart2350,iaNLTEStart(iI))
  iaNLTEStart(iI)     = iaNLTEStart2350(iI)
  WRITE(kStdWarn,*) 'NLTE gas ',iI,' ( = gid ', iaNLTEGasID(iI), ') : NLTE starts at layer ',iaNLTEStart(iI)
END DO
WRITE (kStdWarn,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>> '

! check to see that if there is more than ONE gas that is in NONLTE, that all
! these nonLTE gases start at the same layer, or the rad transfer code might
! get peeved!!!!
! well ok, i think my code can handle this event

! now do the output file names, if we need to bloat things
IF (iSetBloat > 0) THEN
  DO iG = 1,80
    caOutBloatFile(iG:iG)  = ' '
    caPlanckBloatFile(iG:iG) = ' '
  END DO
  iJ = 80
  30     CONTINUE
  IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
    iJ = iJ - 1
    GO TO 30
  END IF
  DO iG = 1,iJ
    caOutBloatFile(iG:iG) = caOutname(iG:iG)
    caPlanckBloatFile(iG:iG) = caOutname(iG:iG)
  END DO
  iG = iJ + 1
  caOutBloatFile(iG:iG+5)  = '_bloat'
  caPlanckBloatFile(iG:iG+12) = '_bloat_PLANCK'
END IF

IF ((iSetBloat > 0) .AND. (iDoUpperAtmNLTE > 0)) THEN
  DO iG = 1,80
    caOutUABloatFile(iG:iG)  = ' '
!          caUAPlanckBloatFile(iG:iG)  = ' '
  END DO
  iJ = 80
  35     CONTINUE
  IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
    iJ = iJ - 1
    GO TO 35
  END IF
  DO iG = 1,iJ
    caOutUABloatFile(iG:iG) = caOutname(iG:iG)
!          caUAPlanckBloatFile(iG:iG) = caOutname(iG:iG)
  END DO
  iG = iJ + 1
  caOutUABloatFile(iG:iG+9)  = '_UA_bloat'
!        caUAPlanckBloatFile(iG:iG+16) = '_UA_bloat_PLANCK'
END IF

! now do the output file names, if we need to do the UA and dump out all rads
! and/or optical depths
! right now the code only dumps out rads
IF (iDoUpperAtmNLTE > 0) THEN
  DO iG = 1,80
    caPlanckUAFile(iG:iG)  = ' '
    caOutUAFile(iG:iG) = ' '
  END DO
  iJ = 80
  40     CONTINUE
  IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
    iJ = iJ - 1
    GO TO 40
  END IF
  DO iG = 1,iJ
    caOutUAFile(iG:iG) = caOutname(iG:iG)
    caPlanckUAFile(iG:iG) = caOutname(iG:iG)
  END DO
  iG = iJ + 1
  caOutUAFile(iG:iG+2)  = '_UA'
  caPlanckUAFile(iG:iG+9) = '_UA_PLANCK'
  
END IF

! really no more error checking possible, as we have to actually go ahead and
! open the files as necessary

RETURN
END SUBROUTINE nonlteSLOW_LBL

!************************************************************************
! this subroutine is for KCOMPRESSED LBL nonlte
! this reads in the info for the new spectra for ONE gas

SUBROUTINE nonlteFAST_KCOMP( iRTP,iaGasesIn,iNumGases,raaTemp,raaPress,  &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, iaNLTEGasID,iaNLTEChunks,  &
    iaaNLTEChunks,caaStrongLines, iaNLTEBands,raNLTEstart,  &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,  &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,  &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,  &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,  &
    caOutBloatFile,caOutUABloatFile, caPlanckUAfile,caOutUAfile)


INTEGER, INTENT(IN OUT)                  :: iRTP
INTEGER, INTENT(IN)                      :: iaGasesIn(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNumGases
REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raNLTEstre
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iNLTE_Slow
NO TYPE, INTENT(IN OUT)                  :: iaNLTEGasI
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
NO TYPE, INTENT(IN OUT)                  :: caaStrongL
NO TYPE, INTENT(IN OUT)                  :: iaNLTEBand
NO TYPE, INTENT(IN OUT)                  :: raNLTEstar
NO TYPE, INTENT(IN OUT)                  :: caaaNLTEBa
NO TYPE, INTENT(IN OUT)                  :: caaNLTETem
NO TYPE, INTENT(IN OUT)                  :: caaUpperMi
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iaNLTEStar
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iAllLayers
NO TYPE, INTENT(IN OUT)                  :: iUseWeakBa
NO TYPE, INTENT(IN OUT)                  :: raKSolarAn
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(OUT)                     :: iSetBloat
NO TYPE, INTENT(IN OUT)                  :: caPlanckBl
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
NO TYPE, INTENT(IN OUT)                  :: caOutUABlo
NO TYPE, INTENT(IN OUT)                  :: caPlanckUA
NO TYPE, INTENT(IN OUT)                  :: caOutUAfil
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (-1/-2) model
! iSetBloat        tells us if stick to 0.0025cm-1 or bloat up to 0.0005 cm-1
! raNLTEstrength   tells how strongly to add on all the files (default 1.0)
! iNumNLTEGases    tells number of NLTE gases
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines   line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE  tells if the code does upper atm NLTE
! iAllLayersLTE    tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd  tells the code if use weak background lines as well, or not
INTEGER :: iNLTE_SlowORFast

INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
REAL :: raNLTEstrength(kGasStore)
INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
INTEGER :: iaNLTEBands(kGasStore)
REAL :: raNLTEstart(kGasStore)
CHARACTER (LEN=80) :: caaaNLTEBands(kGasStore,kNumkCompT)
CHARACTER (LEN=80) :: caaNLTETemp(kGasStore)

! the solar angle
REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
CHARACTER (LEN=80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
CHARACTER (LEN=80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
CHARACTER (LEN=80) :: caPlanckUAfile,caOutUAfile

! local variables
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead,iErr,iI,iaDumb(kMaxGas),iaGases(kMaxGas)
REAL :: rH,rHN
INTEGER :: iStrongISO,iStrongJL,iStrongJU,iStrongGASID,iLTEIn
INTEGER :: iNLTEStart2350,iJ,iG

INTEGER :: iBand,iGasID,iNum,iISO,iLineMixBand,iInt,iType,iGas,iDoVoigtChi
DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)

CHARACTER (LEN=80) :: caNONLTETempKC
INTEGER :: iRegr,iMethod,iDefault
REAL :: raTemp(kProfLayer)
CHARACTER (LEN=80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT)

CALL NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands,caaaNLTEBands)

IF (ABS(iSetBloat) /= 1) THEN
  WRITE(kStdErr,*) 'need iSetBloat  = -1,+1',iSetBloat
  CALL DOSTOP
END IF

IF (iSetBloat == +1) THEN
  WRITE(kStdErr,*) 'you have specified you want kCompressed NLTE database'
  WRITE(kStdErr,*) 'need iSetBloat  = -1 if iNLTE_SlowORFast = -2'
  CALL DoStop
END IF
caPlanckBloatFile = 'blankcaPlanckBloatFile'
caOutBloatFile    = 'blankcaOutBloatFile'
caOutUABloatFile  = 'blankcaOutUABloatFile'

IF (ABS(iDoUpperAtmNLTE) /= 1) THEN
  WRITE(kStdErr,*) 'need iDoUpperAtmNLTE  = -1,+1',iDoUpperAtmNLTE
  CALL DOSTOP
END IF

caNONLTETempKC = 'boo_nlte_file'
DO iI = 1,kProfLayer
  raTemp(iI) = raaTemp(iI,1)   !!assume all kinetic temps equal
END DO

IF (iNumNLTEGases > 1) THEN
  WRITE(kStdErr,*) 'kCARTA can only handle ONE NLTE gas ... CO2!!!'
  CALL DoStop
END IF

IF (iNumNLTEGases < 1) THEN
  WRITE(kStdErr,*) 'huh, 0 or -ve number of NLTE gases????'
  CALL DoStop
END IF

caWord='*NONLTE'
iErr=-1

! we are calling this routine a little early in n_main ... need to reassign
! iaGases correctly
! upto this point eg if gas IDs were 1,3,5,22,51 then
! iaGasesIn(1)=iaGasesIn(3)=iaGasesIn(5)=iaGasesIn(22)=iaGasesIn(51) = 1
! all other iaGasesIn(iN) = -1
! we have to redo this so that iaGases contains the LIST
! iaGases(1,2,3,4,5) = 1,3,5,22,51  all else -1
DO iInt=1,kMaxGas
  iaDumb(iInt)=iaGasesIn(iInt)
  iaGases(iInt)=-1
END DO

iType=1
DO iInt=1,kMaxGas
  IF (iaDumb(iInt) > 0) THEN
    iaGases(iType) = iInt
    iType = iType + 1
  END IF
END DO

5030 FORMAT(A130)

iNumLinesRead=0
13   IF (iNumLinesRead > 0) THEN
  iErr=1
  WRITE(kStdErr,5010) caWord
  CALL DoSTOP
END IF
5010 FORMAT('Error reading section ',A7)

DO iI = 1,kGasStore
  IF (ABS(raNLTEstrength(iI)) > 10.0) THEN
    WRITE(kStdWarn,*)'gas',iI, 'multiply NLTE data by',raNLTEstrength(iI)
  END IF
END DO

iNumLinesRead=1
! read how many gases will have new nonLTE spectroscopy
iErr=-1
IF ((iNumNLTEGases < 1) .OR. (iNumNLTEGases > kGasStore)) THEN
  WRITE(kStdErr,*)'need a valid number of gases in *NONLTE!!'
  WRITE(kStdErr,*)'must be > 0, < kGasStore!!'
  WRITE(kStdErr,*)'please check and retry!'
  CALL DoSTOP
END IF

DO iLTEIn=1,iNumNLTEGases
  IF ((iaNLTEGasID(iLTEIn) < 1) .OR. (iaNLTEGasID(iLTEIn) > kMaxGas)) THEN
    WRITE(kStdErr,*)'iaNLTEGasID(',iLTEIn,') = ',iaNLTEGasID(iLTEIn)
    WRITE(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry (subr nonlteFAST_KCOMP)!'
    CALL DoSTOP
  END IF
  IF ((iaNLTEChunks(iLTEIn) < 1) .OR.  &
        (iaNLTEChunks(iLTEIn) > kNumkCompT)) THEN
    WRITE(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
    WRITE(kStdErr,*)'please check and retry!'
    CALL DoSTOP
  END IF
END DO

! no need to figure out above which height the NLTE starts
! as this is 30 km (already in kCOMP DATABASE
WRITE(kStdWarn,*) 'kCompressed NLTE database starts NLTE at 30 km '

! now do the output file names, if we need to do the UA and dump out all rads
! and/or optical depths
! right now the code only dumps out rads
IF (iDoUpperAtmNLTE > 0) THEN
  DO iG = 1,80
    caPlanckUAFile(iG:iG)  = ' '
    caOutUAFile(iG:iG) = ' '
  END DO
  iJ = 80
  40     CONTINUE
  IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
    iJ = iJ - 1
    GO TO 40
  END IF
  DO iG = 1,iJ
    caOutUAFile(iG:iG) = caOutname(iG:iG)
    caPlanckUAFile(iG:iG) = caOutname(iG:iG)
  END DO
  iG = iJ + 1
  caOutUAFile(iG:iG+2)  = '_UA'
  caPlanckUAFile(iG:iG+9) = '_UA_PLANCK'
  
END IF

! really no more error checking possible, as we have to actually go ahead and
! open the files as necessary

RETURN
END SUBROUTINE nonlteFAST_KCOMP

!************************************************************************
! this subroutine finds out where the Dave Edwards/ Manuel Lopez Puertas
! profiles starts to differs from the input temperature profile

SUBROUTINE FindNLTEHeight(iLTEIn,caaNLTETemp,  &
    iStrongISO,iStrongJL,iStrongJU,iStrongGasID, iaGases,raaTemp,raaPress,  &
    raPressLevels,raLayerHeight,iProfileLayers,rH,iBand)


INTEGER, INTENT(IN OUT)                  :: iLTEIn
NO TYPE, INTENT(IN OUT)                  :: caaNLTETem
INTEGER, INTENT(IN OUT)                  :: iStrongISO
INTEGER, INTENT(IN)                      :: iStrongJL
INTEGER, INTENT(IN)                      :: iStrongJU
NO TYPE, INTENT(IN OUT)                  :: iStrongGas
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: rH
INTEGER, INTENT(IN OUT)                  :: iBand
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
CHARACTER (LEN=80) :: caaNLTETemp(kGasStore)

! this is pressure levels info
REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
! these are the strongest NLTE band parameters (CO2 at 4 um)
INTEGER :: iStrongGasID
! these are input gas profiles
INTEGER :: iProfileLayers


! input/output variable : rH in meters


! local variables
INTEGER :: iL,iG,raTemp(kProfLayer),iLp1,iLp2
REAL :: rL0,rLp1,rLp2,rDeltaThreshold
CHARACTER (LEN=80) :: caFname                !!! file to read
DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN) !!! quantum numbers
INTEGER :: iaJ_UorL(kHITRAN)           !!! matched upper or lower quant nos
INTEGER :: iGasID, iISO                !!! what to read
INTEGER :: iAllOrSome                  !!! read NLTE (+1) or LTE (-1) stuff?
REAL :: raTTemp1(kNLTEProfLayer),raTPress1(kNLTEProfLayer)
REAL :: raNLTEtemp1(kNLTEProfLayer),raQtips1(kNLTEProfLayer)
INTEGER :: iNumVibLevels

REAL :: raLTE(kProfLayer+1),raNLTE(kProfLayer+1),rH0,rH1,rH2
REAL :: raDeltaUse(kProfLayer+1),raDeltaT(kProfLayer+1)
REAL :: rT_tropopause,raPavg(kProfLayer)
INTEGER :: iT_tropopause,i600,iI,iStart
DOUBLE PRECISION :: dVibCenter
REAL :: raNLTE_STD(kProfLayer),raLTE_STD(kProfLayer)
CHARACTER (LEN=100) :: ca1,ca2

rDeltaThreshold = 0.15   !!need (Tv - Tk) > rDeltaThreshold for NLTE

WRITE(kStdWarn,*) 'User specified start NLTE height = ',rH
iStart = (kProfLayer-iProfileLayers+1)

DO iI = iStart,kProfLayer
  raPavg(iI) = raaPress(iI,1)*kAtm2mb
END DO

!!!!fill up the 1 -> iProfileLayers-1 indices w/ increasing dummy numbers
DO iI = iStart-1,1,-1
  raPavg(iI) = raPavg(iI+1) + 10
END DO

rH0 = rH
!!! set up the atmosphere profile
iG = iStrongGasID

rH1 = 0.0
rH2 = 0.0
DO iL = iStart,kProfLayer
  rH1 = rH1 + raaTemp(iL,1)     !!!! should be water temp
  rH2 = rH2 + raaTemp(iL,iG)    !!!! should be CO2 temp
END DO
rH1 = rH1/iProfileLayers
rH2 = rH2/iProfileLayers
IF ((rH2 - rH1) < 10.0) THEN
!!!! only have CO2 in profile
  iG = 1
END IF

daJU(1) = iStrongJU * 1.0D0
daJL(1) = iStrongJL * 1.0D0
caFName = caaNLTETemp(iLTEIn)

CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL,  &
    iStrongGasID,iStrongISO,+1,  &
    iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTEtemp1,dVibCenter)

WRITE(kStdWarn,*) 'iGASID,ISO,iLSGQ,iUSGQ,vibcntr = ',  &
    iStrongGasID,iStrongISO,iStrongJL,iStrongJU,dVibCenter

!swap so pressures are increasing
IF (raTPress1(1) > raTPress1(iNumVibLevels)) THEN
!!rearrange from small to large
  DO iL = 1,iNumVibLevels
    raQtips1(iL) = raTPress1(iNumVibLevels-iL+1)
  END DO
  DO iL = 1,iNumVibLevels
    raTPress1(iL) = raQTips1(iL)
  END DO
  
  DO iL = 1,iNumVibLevels
    raQtips1(iL) = raTTemp1(iNumVibLevels-iL+1)
  END DO
  DO iL = 1,iNumVibLevels
    raTTemp1(iL) = raQTips1(iL)
  END DO
  
  DO iL = 1,iNumVibLevels
    raQtips1(iL) = raNLTEtemp1(iNumVibLevels-iL+1)
  END DO
  DO iL = 1,iNumVibLevels
    raNLTEtemp1(iL) = raQTips1(iL)
  END DO
END IF

CALL rspl(raTPress1,raTTemp1,iNumVibLevels,raPavg,raLTE,kProfLayer)
CALL rspl(raTPress1,raNLTEtemp1,iNumVibLevels,raPavg,raNLTE,kProfLayer)

DO iL = 1,kProfLayer
  raDeltaT(iL)   = raNLTE(iL)-raLTE(iL)
  raDeltaUse(iL) = raNLTE(iL)-raaTemp(iL,iG)
END DO
iL = kProfLayer+1

WRITE(kStdWarn,*) 'start checking LTE/NLTE profile differences : '

iL = 1
30   CONTINUE
IF (raPressLevels(iL) < 1.0) THEN !!!no info here, skip
  iL = iL + 1
  GO TO 30
END IF
WRITE(kStdWarn,*) 'surface level at  ',iL,' = ',raPressLevels(iL),' mb'

!!!first few layers might have a temperature inversion
40   CONTINUE
IF (raPressLevels(iL) > 600) THEN !!!lower than 4.5 km, skip
  iL = iL + 1
  GO TO 40
END IF
WRITE(kStdWarn,*) '4.5 km roughly at ',iL,' = ',raPressLevels(iL),' mb'
i600 = iL

60   CONTINUE
IF ((raDeltaUse(iL) <= rDeltaThreshold).AND.(iL < kProfLayer)) THEN
!!! delta(profile) - delta(GENLN2 NLTE) < 1 K
!!! profile > GENLN2 LTE
  iL = iL + 1
  GO TO 60
END IF

! check to make sure the next 2 or 3 layers are also in NLTE, ie they also
! have (Tv-Tk) > 0.2
IF (iL <= (kProfLayer-2)) THEN
  rL0 =  raDeltaUse(iL)
  iLp1 = iL + 1
  rLp1 = raDeltaUse(iLp1)
  iLp2 = iL + 2
  rLp2 = raDeltaUse(iLp2)
  IF ((rL0 <= rDeltaThreshold) .OR. (rLp1 <= rDeltaThreshold)  &
        .OR. (rLp2 <= rDeltaThreshold)) THEN
!!! was a false set; the next couple or so layers are still in LTE;
    iL = iL + 1
    GO TO 60
  END IF
END IF

rH = raLayerHeight(iL)
IF (iL >= kProfLayer-2) rH = 45000.0   !!!default to 45 km
WRITE(kStdWarn,*) 'start doing NLTE profile at layer ',  &
    iL,', p = ',raPressLevels(iL),' mb, h = ', rH/1000,' km'

IF (rH > rH0) THEN
!!! oops ... let us start where user specified, instead of
!!! higher up in the atmosphere!
  WRITE(kStdWarn,*) 'User wants NLTE to start at ',rH0/1000,' ,not ',rH/1000,' km'
  rH = rH0
END IF

WRITE(kStdWarn,*) 'Start NLTE at (user supplied VS 2350[T-Tvib]) height = ',rH0/1000,' vs ',rH/1000,' km'

ca1 = 'iI   Pavg    Tk(klayers) | Tk(VibT)          dT  |     Tv           dT '
ca2 = '-------------------------|-----------------------|----------------------'

IF (iBand == 1) THEN
  CALL GetUSSTD_2350(raPavg,iStart,daJL,daJU,iaJ_UorL,raLTE_STD,raNLTE_STD)
  WRITE(kStdWarn,*) ca1
  WRITE(kStdWarn,*) ca2
  DO iI = iStart, kProfLayer
    WRITE(kStdWarn,1234) iI,raPavg(iI),raaTemp(iI,1),'|',  &
        raLTE(iI), raLTE(iI)-raaTemp(iI,1),'|',  &
        raNLTE(iI),raNLTE(iI)-raaTemp(iI,1)
  END DO
END IF

IF (iBand == 1) THEN
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*) ' Comparisons of Tk,tNLTE vs USSTD : '
  ca1 = 'iI   Pavg    Tk(klayers)      TStd         dT  |     Tv      TvSTD      dTv'
  ca2 = '-----------------------------------------------|----------------------------'
  
  WRITE(kStdWarn,*) ca1
  WRITE(kStdWarn,*) ca2
  DO iI = iStart, kProfLayer
    WRITE(kStdWarn,1250) iI,raPavg(iI),raaTemp(iI,1),raLTE_STD(iI),raaTemp(iI,1)-raLTE_STD(iI),'|',  &
        raNLTE(iI), raNLTE_STD(iI),raNLTE(iI)-raNLTE_STD(iI)
  END DO
END IF

1234 FORMAT(I3,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1)
1250 FORMAT(I3,' ',4(F10.5,' '),A1,' ',3(F10.5,' '))

RETURN
END SUBROUTINE FindNLTEHeight

!************************************************************************
! simple routine to get US STD 2350 NLTE temps

SUBROUTINE GetUSSTD_2350(raPavg,iStart,daJL,daJU,iaJ_UorL,raLTE_STD,raNLTE_STD)


REAL, INTENT(IN OUT)                     :: raPavg(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iStart
DOUBLE PRECISION, INTENT(IN OUT)         :: daJL(kHITRAN)
DOUBLE PRECISION, INTENT(IN OUT)         :: daJU(kHITRAN)
INTEGER, INTENT(IN OUT)                  :: iaJ_UorL(kHITRAN)
REAL, INTENT(OUT)                        :: raLTE_STD(kProfLayer)
REAL, INTENT(OUT)                        :: raNLTE_STD(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

CHARACTER (LEN=80) :: caFname                !!! file to read
DOUBLE PRECISION :: !!! quantum numbers
INTEGER :: !!! matched upper or lower quant nos
INTEGER :: iGasID, iISO                !!! what to read
INTEGER :: iAllOrSome                  !!! read NLTE (+1) or LTE (-1) stuff?
REAL :: raTTemp1(kNLTEProfLayer),raTPress1(kNLTEProfLayer)
REAL :: raNLTEtemp1(kNLTEProfLayer),raQtips1(kNLTEProfLayer)
INTEGER :: iNumVibLevels, iStrongGasID,iStrongISO

DOUBLE PRECISION :: dVibCenter
INTEGER :: iI

REAL :: raPjunk(kProfLayer),raTJunk(kProfLayer)
REAL :: raLogTPress1(kNLTEProfLayer),raTX1(kNLTEProfLayer),raJunk(kNLTEProfLayer)
REAL :: rX,rY,raTX2(kNLTEProfLayer)

caFName = 'xnlte_1_1_1_6_sol_0.genln2'
CALL concatCA80(caAuxNLTERefsPath,caFName)

iStrongGasID = 2
iStrongISO = 1
CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL,  &
    iStrongGasID,iStrongISO,+1,  &
    iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTEtemp1,dVibCenter)

DO iI = 1,iNumVibLevels
  raLogTPress1(iI) = LOG(raTPress1(iI))
!        print *,iI,raLogTPress1(iI),raTPress1(iI),raTTemp1(iI),raNLTEtemp1(iI)
END DO

DO iI = 1,iNumVibLevels
  raLogTPress1(iI) = LOG(raTPress1(iNumVibLevels-iI+1))
  raTX1(iI)        = raNLTEtemp1(iNumVibLevels-iI+1)
  raTX2(iI)        = raTTemp1(iNumVibLevels-iI+1)
END DO

!! now interp this onto raPavg
DO iI = iStart,kProfLayer
!        print *,iI,iStart,kProfLayer,iI-iStart+1
  raPJunk(iI-iStart+1) = raPAvg(iI)
  rX = LOG(raPAvg(iI))
  CALL rsplin(raLogTPress1,raTX1,raJunk,iNumVibLevels,rX,rY)
  raNLTE_STD(iI) = rY
  CALL rsplin(raLogTPress1,raTX2,raJunk,iNumVibLevels,rX,rY)
  raLTE_STD(iI) = rY
END DO

RETURN
END SUBROUTINE GetUSSTD_2350

!************************************************************************
! this subroutine maps the bands in nml file caaaNLTEBandsOrig to files

SUBROUTINE NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands,  &
    caaaNLTEBands)


NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEGasI
NO TYPE, INTENT(IN OUT)                  :: iaNLTEBand
NO TYPE, INTENT(IN OUT)                  :: caaaNLTEBa
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars
INTEGER :: iNumNLTEGases,iaNLTEGasID(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NLTE bands bad boys
INTEGER :: iaNLTEBands(kGasStore)
! input/output vars
CHARACTER (LEN=80) :: caaaNLTEBands(kGasStore,kNumkCompT)

! local vars
INTEGER :: iI,iJ,iK,iGasID,iISO,iLSGQ,iUSGQ,iUMBC_LBL
CHARACTER (LEN=80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT),caX,caY
CHARACTER (LEN=4) :: ca4
INTEGER :: iaBandID(kNumkCompT),iCount,iCountX,iFound

DO iJ = 1,kNumkCompT
  iaBandID(iJ) = -1
END DO

DO iI = 1,kGasStore
  DO iJ = 1,kNumkCompT
    caaaNLTEBandsOrig(iI,iJ) = caaaNLTEBands(iI,iJ)
  END DO
END DO

DO iI = 1,iNumNLTEGases
  iCount = 1
  DO iJ = 1,iaNLTEBands(iI)
    caX = caaaNLTEBandsOrig(iI,iJ)
    READ(caX,*) iGasID,iIso,iLSGQ,iUSGQ
    IF (iGasID /= 2) THEN
      WRITE(kStdErr,*) 'Need GasID = 2 in nml variable caaaNLTEBand'
      CALL DoStop
!                                            iso iLSGQ iUSGQ umbc-lbl
!--------------------------------------------------------------------
    ELSE IF ((iLSGQ == 1) .AND. (iUSGQ == 9)) THEN
      IF (iIso == 1) iUMBC_LBL = 2350          !1    1    9    2350
      IF (iIso == 2) iUMBC_LBL = 2351          !2    1    9    2351
      IF (iIso == 3) iUMBC_LBL = 2352          !3    1    9    2352
      IF (iIso == 4) iUMBC_LBL = 2355          !4    1    9    2355
      IF ((iIso > 4) .OR. (iIso < 1)) THEN
        WRITE(kStdErr,*) ' SUBROUTINE NLTEBandMapper   1,9 isotope error'
        WRITE(kStdErr,*) ' isotope = ',iIso
        CALL DoStop
      END IF
      
    ELSE IF ((iLSGQ == 4) .AND. (iUSGQ == 24)) THEN
      IF (iIso == 1) iUMBC_LBL = 2310          !1    4   24    2310
      IF (iIso == 2) iUMBC_LBL = 2311          !2    4   24    2311
      IF ((iIso > 2) .OR. (iIso < 1)) THEN
        WRITE(kStdErr,*) ' SUBROUTINE NLTEBandMapper  4,24 isotope error'
        WRITE(kStdErr,*) ' isotope = ',iIso
        CALL DoStop
      END IF
      
    ELSE IF ((iLSGQ == 2) .AND. (iUSGQ == 16)) THEN
      IF (iIso == 1) iUMBC_LBL = 2320          !1    2   16    2320
      IF (iIso == 2) iUMBC_LBL = 2321          !2    2   16    2321
      IF (iIso == 3) iUMBC_LBL = 2322          !2    2   16    2322
      IF ((iIso > 3) .OR. (iIso < 1)) THEN
        WRITE(kStdErr,*) ' SUBROUTINE NLTEBandMapper  2,16 isotope error'
        WRITE(kStdErr,*) ' isotope = ',iIso
        CALL DoStop
      END IF
      
    ELSE IF ((iLSGQ == 3) .AND. (iUSGQ == 23)) THEN
      IF (iIso == 1) iUMBC_LBL = 2353          !1    3   23    2353
      IF (iIso == 2) iUMBC_LBL = 2253          !2    3   23    2253
      IF ((iIso > 2) .OR. (iIso < 1)) THEN
        WRITE(kStdErr,*) ' SUBROUTINE NLTEBandMapper  3,23 isotope error'
        WRITE(kStdErr,*) ' isotope = ',iIso
        CALL DoStop
      END IF
      
    ELSE IF ((iLSGQ == 5) .AND. (iUSGQ == 25)) THEN
      IF (iIso == 1) iUMBC_LBL = 2354          !1    5   25    2354
      IF (iIso == 2) iUMBC_LBL = 2254          !2    5   25    2254
      IF ((iIso > 2) .OR. (iIso < 1)) THEN
        WRITE(kStdErr,*) ' SUBROUTINE NLTEBandMapper  3,23 isotope error'
        WRITE(kStdErr,*) ' isotope = ',iIso
        CALL DoStop
      END IF
      
!!! weak hot bands that Manuel LopezPuertas suggested to use
    ELSE IF (iIso == 1) THEN
      IF ((iLSGQ == 2) .AND. (iUSGQ == 15)) iUMBC_LBL = 2110
      IF ((iLSGQ == 3) .AND. (iUSGQ == 25)) iUMBC_LBL = 2120
      IF ((iLSGQ == 4) .AND. (iUSGQ == 22)) iUMBC_LBL = 2130
      IF ((iLSGQ == 5) .AND. (iUSGQ == 23)) iUMBC_LBL = 2140
      IF ((iLSGQ == 3) .AND. (iUSGQ == 22)) iUMBC_LBL = 2150
      IF ((iLSGQ == 6) .AND. (iUSGQ == 36)) iUMBC_LBL = 2160
      IF ((iLSGQ == 7) .AND. (iUSGQ == 37)) iUMBC_LBL = 2170
      IF ((iLSGQ == 8) .AND. (iUSGQ == 38)) iUMBC_LBL = 2180
      
    ELSE
      WRITE(kStdErr,*) 'trying to figure out CO2 band ....',caX
      WRITE(kStdErr,*) 'NLTEBandMapper cannot decipher',iIso,iLSGQ,iUSGQ
      CALL DoStop
    END IF
    
!!file name caY are of the form 'blah2350.dat'
    WRITE(ca4,40) iUMBC_LBL
    40       FORMAT(I4)
    DO iK = 1,80
      caY(iK:iK) = ' '
    END DO
    caY(1:4) = ca4
    caY(5:8) = '.dat'
    
    CALL ConcatCA80(caStrongLineParams,caY)
    caaaNLTEBands(iI,iJ) = caY
    WRITE(kStdWarn,*) iI,iJ,'-> NLTE : iISO,iLSGQ,iUSGQ,UMBC-LBL iD: ',  &
        iISO,iLSGQ,iUSGQ,iUMBC_LBL
    
    IF (iCount == 1) THEN
      iaBandID(iCount) = iUMBC_LBL
      iCount = iCount + 1
    ELSE
!!! check to make sure this band_id is not used already
      iFound = -1
      DO iCountX = 1,iCount-1
        IF (iaBandID(iCountX) == iUMBC_LBL) THEN
          iFound = +1
          WRITE(kStdErr,*) 'NLTEBandMapper found ',iCountX,iUMBC_LBL
        END IF
      END DO
      IF (iFound == +1) THEN
        WRITE(kStdErr,*) 'NLTEBandMapper found double count iUMBC_LBL'
        CALL DoStop
      ELSE
        iaBandID(iCount) = iUMBC_LBL
        iCount = iCount + 1
      END IF
    END IF
    
  END DO
END DO

RETURN
END SUBROUTINE NLTEBandMapper

!************************************************************************
! subroutine adds in the molgases

SUBROUTINE add_molgas(iNgas,iaGasesNL,iaMOLgases)


INTEGER, INTENT(IN OUT)                  :: iNgas
INTEGER, INTENT(IN OUT)                  :: iaGasesNL(kGasComp)
INTEGER, INTENT(OUT)                     :: iaMOLgases(kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input/output

! output


! local
INTEGER :: iMax,iC,iCC,iaInDataBase(kGasComp),iTag,iErr,iaTemp(kMaxGas)
INTEGER :: iNotMainGas,iNgasesCheck,iInt,iWhichKC,iYesWV,i101,i102,i103
INTEGER :: iCheckCompDataBase,MainGas,iKLBLRTMgases,iMax_molgas

iWHichKC = iNgas
iKLBLRTMgases = 32

IF (iWhichKC == -1) THEN
!! this is kc1.15 or later
! use ALL gases in the compressed database for v1.15+
!  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9
!         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,28,X,X,31,32,...42
!         PLUS kSelf,kFor,kHeavyWater (101,102,103)
  iMax = kGasComp
  WRITE(kStdWarn,*) 'including all lbl gases from 1 to ', kGasComp
  WRITE(kStdWarn,*) ' plus 101,102,103 for Self,Forn cont, heavy water'
ELSE IF ((iWhichKC == -5) .OR. (iWhichKC == -6)) THEN
!! this is kc1.14 or earlier
!! also LBLRTM style
! use ALL gases in the compressed database for c1.14-
!  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9
!         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,28,X,X,31
!         PLUS kSelf,kFor,kHeavyWater (101,102,103)
!        iMax = 31
!        write(kStdWarn,*) 'including v114- lbl gases from 1 to ',iMax
  iMax = iKLBLRTMgases
  WRITE(kStdWarn,*) 'including LBLRTM gases from 1 to ',iMax
  WRITE(kStdWarn,*) ' plus 101,102,103 for Self,Forn cont, heavy water'
ELSE IF (iWHichKC > 0) THEN
  WRITE(kStdWarn,*) 'using user sepcified list'
ELSE
  WRITE(kStdErr,*) 'need iWhichKC = -1 or -5/-6 or > 0'
  CALL DoStop
END IF

IF (iWhichKC < 0) THEN
  iNgas = 0
  IF (iWhichKC == -1) THEN
    iMax_molgas = kGasComp
  ELSE IF ((iWhichKC == -5) .OR. (iWhichKC == -6)) THEN
    iMax_molgas = iKLBLRTMgases
  END IF
  
  DO iC = 1,kGasComp
    iaInDataBase(iC) = -1
  END DO
  DO iC = 1,iMax_molgas
    iTag = -1
    iCC = iCheckCompDataBase(iC,-100.0,-100.0,iTag,iErr)
    IF (iCC > 0) THEN
      iNgas = iNgas+1
      iaInDataBase(iC) = 1
    END IF
  END DO
  
! now based on which gases were found, reset array iaTemp
  iCC = 1
  DO iC = 1,kGasComp
    IF (iaInDataBase(iC) > 0) THEN
      WRITE(kStdWarn,*)'Including gasID ',iC,' in comp database'
      iaTemp(iCC) = iC
      iCC = iCC+1
    END IF
  END DO
  
  IF (kCKD >= 0) THEN
    DO iC = kNewGasLo,kNewGasHi
      WRITE(kStdWarn,*)'Including gasID ',iC,' in comp database'
      iaTemp(iCC) = iC
      iCC = iCC+1
      iNgas = iNgas+1
    END DO
  END IF
  
!!!heavy water
  iC = kNewGasHi+1
  WRITE(kStdWarn,*)'Including gasID ',iC,' in comp database'
  iaTemp(iCC) = iC
  iCC = iCC+1
  iNgas = iNgas+1
  
ELSE IF (iWhichKC > 0) THEN
  iYesWV = -1
  i101 = -1
  i102 = -1
  i103 = -1
  DO iC = 1,iWhichKC
    iaTemp(iC) = iaGasesNL(iC)
    IF (iaTemp(iC). EQ. 1) iYesWV = +1
    IF (iaTemp(iC). EQ. 101) i101 = +1
    IF (iaTemp(iC). EQ. 102) i102 = +1
    IF (iaTemp(iC). EQ. 103) i103 = +1
  END DO
  IF ((kCKD >= 0) .AND. (iYesWV > 0) .AND. (i101 == -1)) THEN
!! need to add on g101
    iWhichKC = iWhichKC + 1
    iaTemp(iWhichKC) = 101
    iaGasesNL(iWhichKC) = 101
  END IF
  IF ((kCKD >= 0) .AND. (iYesWV > 0) .AND. (i102 == -1)) THEN
!! need to add on g102
    iWhichKC = iWhichKC + 1
    iaTemp(iWhichKC) = 102
    iaGasesNL(iWhichKC) = 102
  END IF
  IF ((kCKD >= 0) .AND. (iYesWV > 0) .AND. (i102 == -1)) THEN
!! need to add on g103
    iWhichKC = iWhichKC + 1
    iaTemp(iWhichKC) = 103
    iaGasesNL(iWhichKC) = 103
  END IF
  iNgas = iWhichKC
END IF

! check the molecular ID's are between 1 and kGasComp , and
! kNewGasLo and kNewGasHi+1
! (should be in the compressed data base)
! if iNgas = -1, everything should have been set correctly above; if user
! entered in the GasIDs him/her self, there could be mistakes
DO iC = 1,iNgas
  iNotMainGas = MainGas(iaTemp(iC))
  IF (iNotMainGas < 0) THEN
! gas does not exist in the compressed base ... stop
    WRITE(kStdErr,777) iC,iaTemp(iC),1,kGasComp
    777       FORMAT('Error in MOLGAS!! iC = ',I2,' found Gas ID   =  ',I2,'(  &
        MOLGAS ID''S SHOULD BE BETWEEN ',I2,' AND ',I2 ,')')
    CALL DoSTOP
  END IF
! check to see if kcomp files do exist
  iTag = -1
  IF (iaTemp(IC) <= kGasComp) THEN
    iCC = iCheckCompDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
    IF (iCC < 0) THEN
      WRITE(kStdWarn,780) iaTemp(iC)
    END IF
  ELSE IF (iaTemp(IC) <= kNewGasHi) THEN
    WRITE(kStdWarn,785) iaTemp(iC)
  ELSE IF (iaTemp(IC) == kNewGasHi+1) THEN
    WRITE(kStdWarn,787) iaTemp(iC)
  END IF
END DO
780  FORMAT('Warning! GasID ',I3,' not in compressed data base')
785  FORMAT('Warning! GasID ',I3,' is a new continuum gas')
787  FORMAT('Warning! GasID ',I3,' is heavy water (gasID  =  1, HITRAN isotope  =  4)')

iNgasesCheck = iNgas

! set the identities of the gases whose abs spectra are in files
! iaMOLgases keeps track of how many times the GasID has been found, so that
! no double counting of the gases (between MOLGAS and XSCGAS) is done
DO iInt = 1,iNgas
  iaMOLgases(iaTemp(iInt)) = iaMOLgases(iaTemp(iInt))+1
END DO

! check to see that the same gas ID has not been entered more than once
DO iInt = 1,iNgas
  IF ((iaMOLgases(iaTemp(iInt)) > 1) .OR.  &
        (iaMOLgases(iaTemp(iInt)) < 0)) THEN
    WRITE(kStdErr,*) 'Gas ID',iaTemp(iInt),' entered more than'
    WRITE(kStdErr,*)'once. Please check *MOLGAS and retry'
    CALL DoSTOP
  END IF
END DO

iNgas = iNgasesCheck
WRITE(kStdWarn,*) 'MOLGAS ... gases stored are ',iNgas
DO iInt = 1,kGasComp
  IF (iaMOLgases(iInt) > 0) THEN
    WRITE(kStdWarn,*) '     going to use compressed database gas ',iInt
  END IF
END DO

DO iInt = kNewGasLo,kNewGasHi
  IF (iaMOLgases(iInt) > 0) THEN
    WRITE(kStdWarn,*) '     going to use new continuum gas ',iInt
  END IF
END DO

iInt = kNewGasHi+1
IF (iaMOLgases(iInt) > 0) THEN
  WRITE(kStdWarn,*) '     going to use new heavy water gas ',iInt
END IF

IF (kCKD >= 0) THEN
  DO iInt = kNewGasLo,kNewGasHi
    IF (iaMOLgases(iInt) <= 0) THEN
      WRITE(kStdWarn,*) 'Cannot have CKD on and gasIDs 101/102 unused'
      WRITE(kStdErr,*) 'Cannot have CKD on and gasIDs 101/102 unused'
      CALL DoSTOP
    END IF
  END DO
END IF

RETURN
END SUBROUTINE add_molgas

!************************************************************************
! subroutine adds in the xsecgases

SUBROUTINE add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)


INTEGER, INTENT(IN OUT)                  :: iNXsec
INTEGER, INTENT(IN)                      :: iaLXsecNL(kGasXSecHi-kGasXSecLo+1
INTEGER, INTENT(OUT)                     :: iaXSCgases(kMaxGas)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input/output



! local
INTEGER :: iInt,iErr,iC,iCC,iCheckXsecDataBase,iNXsecCheck
INTEGER :: iaTemp(kMaxGas),iTag
INTEGER :: iaInDataBase(kMaxLayer),iWhichXSC,iKLBLRTMgases

iWhichXSC = iNXsec

IF (iWhichXSC == -1) THEN
  WRITE(kStdWarn,*) 'including all xsc gases from ',kGasXsecLo, ' to ', kGasXsecHi
! use all gases in the xsec database
!  check to see the following reference profiles exist : 51..82
  iNxsec = 0
  DO iC = kGasXsecLo,kGasXsecHi
    iaInDataBase(iC) = -1
  END DO
  DO iC = kGasXsecLo,kGasXsecHi
    iTag = -1
    iCC = iCheckXsecDataBase(iC,-100.0,-100.0,iTag,iErr)
    IF (iCC > 0) THEN
      iNxsec = iNxsec+1
      iaInDataBase(iC) = 1
    END IF
  END DO
! now based on which gases were found, reset array iaTemp
  iCC = 1
  DO iC = kGasXsecLo,kGasXsecHi
    IF (iaInDataBase(iC) > 0) THEN
      WRITE(kStdWarn,*) 'Including gasID ',iC,' in xsec database'
      iaTemp(iCC) = iC
      iCC = iCC+1
    END IF
  END DO
  
ELSE IF ((iWhichXSC == -5) .OR. (iWhichXSC == -6))THEN
  iKLBLRTMgases = 63
  WRITE(kStdWarn,*) 'including LBLRTM xsc gases from ',kGasXsecLo,' to ',iKLBLRTMgases
! use all gases in the xsec database
!  check to see the following reference profiles exist : 51..63
  iNxsec = 0
  DO iC = kGasXsecLo,iKLBLRTMgases
    iaInDataBase(iC) = -1
  END DO
  DO iC = kGasXsecLo,iKLBLRTMgases
    iTag = -1
    iCC = iCheckXsecDataBase(iC,-100.0,-100.0,iTag,iErr)
    IF (iCC > 0) THEN
      iNxsec = iNxsec+1
      iaInDataBase(iC) = 1
    END IF
  END DO
! now based on which gases were found, reset array iaTemp
  iCC = 1
  DO iC = kGasXsecLo,iKLBLRTMgases
    IF (iaInDataBase(iC) > 0) THEN
      WRITE(kStdWarn,*) 'Including gasID ',iC,' in xsec database'
      iaTemp(iCC) = iC
      iCC = iCC+1
    END IF
  END DO
  
ELSE IF (iWhichXSC > 0) THEN
  DO iC = 1,iNXsec
    iaTemp(iC) = iaLXsecNL(iC)
  END DO
ELSE
  WRITE(kStdErr,*) 'xscgas4 can only process iNXsec = -5/-6,-1, or > 0'
  CALL DOSTOP
END IF

! check the molecular ID's are between kGasXsecLo and kXsecGasHi
! (should be in the cross sec data base)
! if iNXsec = -1, everything should have been set correctly above; if user
! enetered in the GasIDs him/her self, there could be mistakes
DO iC = 1,iNXsec
  IF ((iaTemp(iC) < kGasXsecLo).OR.(iaTemp(iC) > kGasXsecHi)) THEN
! gas does not exist in the xsec data base ... stop
    WRITE(kStdErr,777) iaTemp(iC),kGasXsecLo,kGasXsecHi
    777      FORMAT('Error in XSCGAS!! found Gas ID  = ',I2,'(  &
        XSCGAS ID''S SHOULD BE BETWEEN ',I2,' AND ',I2 ,')')
    CALL DoSTOP
  END IF
! check to see if cross section data files for this gas does exist
  iTag = -1
  iCC = iCheckXsecDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
  IF (iCC < 0) THEN
    WRITE(kStdWarn,780) iaTemp(iC)
    780       FORMAT('Warning! GasID ',I2,' not in xsec data base')
  END IF
END DO

iNXsecCheck=iNXsec
! set the identities of the gases whose molecular cross sections are in file
! iaXSCgases keeps track of how many times the GasID has been found, so that
! no double counting of the gases (between MOLGAS and XSCGAS) is done
DO iInt = 1,iNXsec
  iaXSCgases(iaTemp(iInt)) = iaXSCgases(iaTemp(iInt))+1
END DO

! check to see that the same gas ID has not been entered more than once
DO iInt = 1,iNXsec
  IF ((iaXSCgases(iaTemp(iInt)) > 1) .OR.  &
        (iaXSCgases(iaTemp(iInt)) < 0)) THEN
    WRITE(kStdErr,*) 'Gas ID ',iaTemp(iInt),' entered more than'
    WRITE(kStdErr,*) 'once Please check *XSCGAS and retry'
    CALL DoSTOP
  END IF
END DO

iNXsec = iNXsecCheck
WRITE(kStdWarn,*)'cross section gases stored are ',iNXsec
DO iInt = kGasXsecLo,kGasXsecHi
  IF (iaXSCgases(iInt) > 0) THEN
    WRITE(kStdWarn,*) 'CrossSect database gas ',iInt
  END IF
END DO

RETURN
END SUBROUTINE add_xsecgas
!************************************************************************