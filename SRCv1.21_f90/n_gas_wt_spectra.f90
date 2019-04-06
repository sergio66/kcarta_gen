! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_gas_wt_spectra

USE basic_common
!USE kpredictVT
USE kcoeff_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile
USE s_writefile
USE n_nonlte_common

IMPLICIT NONE

CONTAINS

! this file reads *MOLGAS,*XSCGAS,*FREQ (very easy),*MIXFIL

!************************************************************************
! this subroutine deals with the 'MOLGAS' keyword
! and basically deals with GasID 1--40

! main output parameter is iaMOLgases
! synopsis IN : iNGas     = number of gasIDs to be read in (if -1, all)
!               iaGasesNL = list of gasIDs
!         OUT : iaMOLgases(iI) = set to 1 if gasID was found in this section
    SUBROUTINE molgas4(iNGas,iaGasesNL,iaMOLgases)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iNGas,iaGasesNL(kGasComp),iaMOLgases(kMaxGas)

! local variables
    INTEGER :: iTag,iErr,iNotMainGas
    INTEGER :: iInt,iaTemp(kMaxGas),iC,iCC,iaInDataBase(kGasComp)
    INTEGER :: iNgasesCheck

    CHARACTER(7) :: caWord

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
      write(kStdWarn,*) 'including user specified lbl gases'
      ! use the molecular ID's from the namelist file
      DO iInt = 1,iNgas
        iNotMainGas = MainGas(iaGasesNL(iInt))
        IF (iNotMainGas < 0) THEN
          write(kStdErr,*) 'Invalid MOLGAS GasID',iaGasesNL(iInt),' entered'
          write(kStdErr,*) 'Please check *MOLGAS and retry'
          write(kStdErr,*) 'Note : If the GasID printed was -100, you '
          write(kStdErr,*) 'probably entered less gasIDs than you promised '
          CALL DoSTOP
        ELSE
          iaTemp(iInt) = iaGasesNL(iInt)
        END IF
      END DO
      CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)   !! user specified gases

    ELSEIF (iNgas == -1) THEN
      CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)   !! latest version, uses 42 molgas

    ELSEIF ((iNGas == -5) .OR. (iNGas == -6)) THEN  !! use 36 gases for LBLRTM
      CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)

    ELSE
      write(kStdErr,*) 'input user iNXsec = ',iNgas
      write(kStdErr,*) '  valid iNgas : -5/-6 (lblrtm) -1 (all) +N select few'
      CALL DoStop
    END IF
            
    RETURN
    end SUBROUTINE molgas4

!************************************************************************
! this subroutine deals with the 'XSCGAS' keyword
! and basically deals with GasID 51--63
! synopsis IN : iNXsec    = number of gasIDs to be read in (if -1, all)
!               iaLXsecNL = list of gasIDs
!         OUT : iaXSCgases(iI) = set to 1 if gasID was found in this section
    SUBROUTINE xscgas4(iNXsec,iaLXsecNL,iaXSCgases)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iNXsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1),iaXSCgases(kMaxGas)

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iInt,iErr,iC,iCC,iNXsecCheck
    INTEGER :: iaTemp(kMaxGas),iTag
    INTEGER :: iaInDataBase(kMaxLayer)
    caWord = '*XSCGAS'

    iaXSCgases = 0
    iaTemp     = 0

! allow iNxsec = -5/-6 ,-1,or > 0 where
! -114  means allow LBL gases which were in older kcarta versions (upto v1.14)
!       which was 51..63 REMOVED
! -5/-6 means allow LBLRTM gases 51..63
! -1    means allow ALL XSC gases from the new HITRAN208 = 51..81
! > 0   means you choose the gases (using iaLXsecNL)

    iErr = -1
    IF (iNxsec > 0) THEN
      write(kStdWarn,*) 'including user specified xsc gases'
      ! use the xsec gas ID's from the namelist file
      DO iInt = 1,iNXsec
        IF ((iaLXsecNL(iInt) < kGasXsecLo) .OR. (iaLXsecNL(iInt) > kGasXsecHi)) THEN
          write(kStdErr,*) 'Invalid XSCGAS GasID',iaLXsecNL(iInt),' entered'
          write(kStdErr,*) 'Please check *XSCGAS and retry'
          write(kStdErr,*) 'Note : If the GasID printed was -100, you '
          write(kStdErr,*) 'probably entered less gasIDs than you promised '
          CALL DoSTOP
        ELSE
          iaTemp(iInt) = iaLXsecNL(iInt)
        END IF
      END DO
      CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! now check the veracity of these gases
    ELSEIF (iNxsec == -1) THEN
      CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! add all xsec gases
    ELSEIF ((iNxsec == -5) .OR. (iNXsec == -6)) THEN
      CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)  !! add LBLRTM xsec gases
    ELSE
      write(kStdErr,*) 'input user iNXsec = ',iNXsec
      write(kStdErr,*) '  valid iNXsec : -5/-6 (lblrtm) -1 (all) +N select few'
      CALL DoStop
    END IF

    RETURN
    end SUBROUTINE xscgas4

!************************************************************************
! this subroutine is for the new spectra
! this reads in the info for the new spectra for ONE gas
    SUBROUTINE spectra4(iNumNewGases, &
    iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks
    INTEGER :: iaNewGasID(kGasStore),iaNewData(kGasStore)
    INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iCount,iErr

    caWord = '*SPECTR'
    iErr = -1

    5030 FORMAT(A130)

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7)

    iNumLinesRead = 1
! read how many gases will have new spectroscopy
    iErr = -1
    IF ((iNumNewGases < 1) .OR. (iNumNewGases > kGasStore)) THEN
      write(kStdErr,*)'need a valid number of gases in *SPECTRA!!'
      write(kStdErr,*)'must be > 0, < kGasStore!!'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF

    DO iCount = 1,iNumNewGases
      IF ((iaNewGasID(iCount) < 1) .OR. (iaNewGasID(iCount) > kMaxGas)) THEN
        write(kStdErr,*) 'iaNewGasID(',iCount,') = ',iaNewGasID(iCount)
        write(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
        write(kStdErr,*)'please check and retry (subr spectra4)!'
        CALL DoSTOP
      END IF
      IF ((iaNewData(iCount) < 1) .OR. (iaNewData(iCount) > kNumkCompT)) THEN
        write(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP
      END IF
    END DO
            
! really no more error checking possible, as we have to actually go ahead and
! open the files as necessary

    RETURN
    end SUBROUTINE spectra4

!************************************************************************
! this subroutine is for nonlte ... does the decision between SLOW or FAST
! this reads in the info for the new spectra for ONE gas
    SUBROUTINE nonlte( &
    iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks, &
    iaaNLTEChunks,caaStrongLines, &
    iaNLTEBands,raNLTEstart, &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
    caOutBloatFile,caOutUABloatFile, &
    caPlanckUAfile,caOutUAfile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iRTP,iaGasesIn(kMaxGas),iNumGases,iSetBloat,iNLTE_SlowORFast
    REAL ::    raaTemp(kProfLayer,kGasStore),raaPress(kProfLayer,kGasStore)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    REAL :: raNLTEstrength(kGasStore)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    REAL :: raNLTEstart(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caOutName
! the solar angle
    REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
    REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
    INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile

    IF (iNLTE_SlowORFast == +1) THEN
      CALL nonlteSLOW_LBL( &
        iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
        raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
        iaNLTEGasID,iaNLTEChunks, &
        iaaNLTEChunks,caaStrongLines, &
        iaNLTEBands,raNLTEstart, &
        caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
        raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
        iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
        raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
        caOutBloatFile,caOutUABloatFile, &
        caPlanckUAfile,caOutUAfile)

    ELSEIF (iNLTE_SlowORFast == -1) THEN
      CALL nonlteFAST_SARTA( &
        iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
        raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
        iaNLTEGasID,iaNLTEChunks, &
        iaaNLTEChunks,caaStrongLines, &
        iaNLTEBands,raNLTEstart, &
        caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
        raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
        iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
        raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
        caOutBloatFile,caOutUABloatFile, &
        caPlanckUAfile,caOutUAfile)

    ELSEIF (iNLTE_SlowORFast == -2) THEN
      CALL nonlteFAST_KCOMP( &
        iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
        raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
        iaNLTEGasID,iaNLTEChunks, &
        iaaNLTEChunks,caaStrongLines, &
        iaNLTEBands,raNLTEstart, &
        caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
        raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
        iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
        raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
        caOutBloatFile,caOutUABloatFile, &
        caPlanckUAfile,caOutUAfile)
    END IF

    RETURN
    end SUBROUTINE nonlte

!************************************************************************
! this subroutine is for FAST version of nonlte
! this reads in the info for the new spectra for ONE gas
    SUBROUTINE nonlteFAST_SARTA( &
    iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks, &
    iaaNLTEChunks,caaStrongLines, &
    iaNLTEBands,raNLTEstart, &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
    caOutBloatFile,caOutUABloatFile, &
    caPlanckUAfile,caOutUAfile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iRTP,iaGasesIn(kMaxGas),iNumGases,iSetBloat,iNLTE_SlowORFast
    REAL ::    raaTemp(kProfLayer,kGasStore),raaPress(kProfLayer,kGasStore)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    REAL :: raNLTEstrength(kGasStore)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    REAL :: raNLTEstart(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caOutName
! the solar angle
    REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
    REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
    INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile

! local vars
    INTEGER :: iI,iJ

    caPlanckBloatFile = 'blankcaPlanckBloatFile'
    caOutBloatFile    = 'blankcaOutBloatFile'
    caOutUABloatFile  = 'blankcaOutUABloatFile'
    caPlanckUAfile    = 'blankcaPlanckUAfile'
    caOutUAfile       = 'blankcaOutUAfile'
            
    IF (iSetBloat > 0) THEN
      iSetBloat = -1
      write(kStdWarn,*) 'Reset iSetBloat to -1 in nonlteFAST_SARTA'
    END IF
    IF (iDoUpperAtmNLTE > 0) THEN
      iDoUpperAtmNLTE = -1
      write(kStdWarn,*) 'Reset iDoUpperAtmNLTE to -1 in nonlteFAST_SARTA'
    END IF
           
    RETURN
    end SUBROUTINE nonlteFAST_SARTA

!************************************************************************
! this subroutine is for GENLN2_LBL nonlte
! this reads in the info for the new spectra for ONE gas
    SUBROUTINE nonlteSLOW_LBL( &
    iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks, &
    iaaNLTEChunks,caaStrongLines, &
    iaNLTEBands,raNLTEstart, &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
    caOutBloatFile,caOutUABloatFile, &
    caPlanckUAfile,caOutUAfile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iRTP,iaGasesIn(kMaxGas),iNumGases,iSetBloat,iNLTE_SlowORFast
    REAL ::    raaTemp(kProfLayer,kGasStore),raaPress(kProfLayer,kGasStore)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    REAL :: raNLTEstrength(kGasStore)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    REAL :: raNLTEstart(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caOutName
! the solar angle
    REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
    REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
    INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iErr,iI,iaDumb(kMaxGas),iaGases(kMaxGas)
    REAL :: rH,rHN
    INTEGER :: iStrongISO,iStrongJL,iStrongJU,iStrongGASID,iLTEIn
    INTEGER :: iNLTEStart2350,iJ,iG

    INTEGER :: iBand,iGasID,iNum,iISO,iLineMixBand,iInt,iType,iGas,iDoVoigtChi
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)

    CHARACTER(80) :: caNONLTETempKC
    INTEGER :: iRegr,iMethod,iDefault
    REAL :: raTemp(kProfLayer)
    CHARACTER(80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT)
    INTEGER :: iJunkNum,iaJunk(kGasStore)

    CALL NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands,caaaNLTEBands)

    IF (abs(iSetBloat) /= 1) THEN
      write(kStdErr,*) 'need iSetBloat  = -1,+1',iSetBloat
      CALL DOSTOP
    END IF
            
    IF (abs(iDoUpperAtmNLTE) /= 1) THEN
      write(kStdErr,*) 'need iDoUpperAtmNLTE  = -1,+1',iDoUpperAtmNLTE
      CALL DOSTOP
    END IF
            
!!! how to predict NLTE temperatures before 2004
    iMethod = 1      !!! lapack  .. not supported anymore
    iMethod = 3      !!! predictor
    iMethod = 2      !!! polynom
! after 2004, just use GRANADA model EXTERNALLY to kCARTA

    iDefault = 2
    IF (iMethod /= iDefault) THEN
      write(kStdErr,*) 'imethod,idefault in subr nonlte = ',iMethod,iDefault
    END IF
    IF (iMethod == 3) THEN
      write(kStdErr,*) 'imethod = 3 still not programmed up'
      CALL DoStop
    END IF
    IF (iMethod == 1) THEN
      write(kStdErr,*) 'imethod = 1 LAPACK has been voided'
      CALL DoStop
    END IF

    caNONLTETempKC = 'boo_nlte_file'
    raTemp = raaTemp(:,1)   !!assume all kinetic temps equal

    IF (iNumNLTEGases > 1) THEN
      write(kStdErr,*) 'kCARTA can only handle ONE NLTE gas ... CO2!!!'
      CALL DoStop
    END IF

    IF (iNumNLTEGases < 1) THEN
      write(kStdErr,*) 'huh, 0 or -ve number of NLTE gases????'
      CALL DoStop
    END IF

    iRegr = -1
    DO iI = 1,iNumNLTEGases
      IF (caaNLTETemp(iI) == 'nlteguess') THEN
        IF (iMethod == +2) THEN
          ! leave iRegr as -1
          CALL polynom_nlte( &
              caOutName,raTemp,raPressLevels,iProfileLayers,raKSolarAngle(1), &
              iRTP,iRegr,caNONLTETempKC)
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
    iaDumb = iaGasesIn
    iaGases = -1

    iType = 1
    DO iInt = 1,kMaxGas
      IF (iaDumb(iInt) > 0) THEN
        iaGases(iType)  =  iInt
        iType = iType + 1
      END IF
    END DO

 5030 FORMAT(A130)

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7)

    DO iI = 1,kGasStore
      IF (abs(raNLTEstrength(iI)) > 10.0) THEN
        write(kStdWarn,*)'gas',iI, 'multiply NLTE data by',raNLTEstrength(iI)
      END IF
    END DO

    iNumLinesRead = 1
! read how many gases will have new nonLTE spectroscopy
    iErr = -1
    IF ((iNumNLTEGases < 1) .OR. (iNumNLTEGases > kGasStore)) THEN
      write(kStdErr,*)'need a valid number of gases in *NONLTE!!'
      write(kStdErr,*)'must be > 0, < kGasStore!!'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF

    DO iLTEIn = 1,iNumNLTEGases
      IF ((iaNLTEGasID(iLTEIn) < 1) .OR. (iaNLTEGasID(iLTEIn) > kMaxGas)) THEN
        write(kStdErr,*)'iaNLTEGasID(',iLTEIn,') = ',iaNLTEGasID(iLTEIn)
        write(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
        write(kStdErr,*)'please check and retry (subr nonlteSLOW_LBL)!'
        CALL DoSTOP
      END IF
      IF ((iaNLTEChunks(iLTEIn) < 1) .OR. (iaNLTEChunks(iLTEIn) > kNumkCompT)) THEN
        write(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
        write(kStdErr,*)'please check and retry!'
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
        GOTO 24
      ELSEIF (rH > 0) THEN
        !!!go ahead and find iStart
        iI = kProfLayer - iProfileLayers + 1
        iI = iI - 1
        iI = max(1,iI)
 23     CONTINUE
        IF ((rH >= raLayerHeight(iI)) .AND. (iI < kProfLayer)) THEN
          iI = iI + 1
          GOTO 23
        END IF
        IF ((rH >= raLayerHeight(iI)) .AND. (iI == kProfLayer)) THEN
          write(kStdErr,*) 'You are starting NLTE too high!!! '
          write(kStdErr,*) 'Max layer height (m) = ',raLayerHeight(iI)
          write(kStdErr,*) 'Your NLTE start (m)  = ',rH
          CALL DoStop
        ELSE
          iaNLTEStart(iLTEIn) = max(iI,1)
          write(kStdWarn,*) '(nml file) user defined NLTE starts at lay ',iI
          write(kStdWarn,*) 'for nlte gas number ',iLTEIn,' = gasID ',iaNLTEGasID(iLTEIn)
          write(kStdWarn,*) 'which is layer at ~ ',rH/1000,' km'
        END IF
      END IF
 24   CONTINUE
    END DO
     
! now figure out above which height the NLTE starts, for the strong bands
! iaNLTEStart2350 is for STRONGEST NLTE bands, in 4 um region ahem!
! ************* this is using the D. Edwards/M. Lopez-Puertas info ************
    iJunkNum    = -1
    rH = 1.0e10             !!!dumb large number, too high!!!!!!
    DO iGas = 1,iNumGases
      iLTEIn = OutsideSpectra(iaGases(iGas),iNumNLTEGases,iaNLTEGasID,iJunkNum,iaJunk,2205.0,605.0,2830.0,20)
      IF (iLTEIn > 0) THEN
        rHN = raNLTEstart(iLTEIn)*1000     !!!NLTE start height in m
        rH = min(rH,rHN)
        !!! now see if we can come up with a better estimate
        iStrongGASID = iaGases(iGas)
        DO iBand = 1,iaNLTEBands(iLTEIn)
          ! read in the lineshape parameters for the band
          CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
                iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
                daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR, &
                iLineMixBand,iDoVoigtChi)
          iStrongISO   = iISO
          iStrongJU    = nint(daJU(1))
          iStrongJL    = nint(daJL(1))
          CALL FindNLTEHeight(iLTEIn,caaNLTETemp, &
            iStrongISO,iStrongJL,iStrongJU,iStrongGasID, &
            iaGases,raaTemp,raaPress, &
            raPressLevels,raLayerHeight,iProfileLayers,rHN,iBand)
          rH = min(rH,rHN)
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
          iI = max(1,iI)
 26       CONTINUE
          IF ((rH >= raLayerHeight(iI)) .AND. (iI < kProfLayer)) THEN
            iI = iI + 1
            GOTO 26
          END IF
          IF ((rH >= raLayerHeight(iI)) .AND. (iI == kProfLayer)) THEN
            write(kStdErr,*) 'You are starting NLTE too high!!! '
            write(kStdErr,*) 'Max layer height (m) = ',raLayerHeight(iI)
            write(kStdErr,*) 'Your NLTE start (m)  = ',rH
            CALL DoStop
          ELSE
            iNLTEStart2350 = iI
            write(kStdWarn,*) 'From reading the NLTE profiles '
            write(kStdWarn,*) 'Strong NLTE starts ',iI,' for gas ',iLTEIn
          END IF
        END IF
      END IF
    END DO

    write (kStdWarn,*) '   '
    write (kStdWarn,*) '   '
    write (kStdWarn,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>> '
    DO iI = 1,iNumNLTEGases
      iaNLTEStart2350(iI) = min(iNLTEStart2350,iaNLTEStart(iI))
      iaNLTEStart(iI)     = iaNLTEStart2350(iI)
      write(kStdWarn,*) 'NLTE gas ',iI,' ( = gid ', iaNLTEGasID(iI), ') : NLTE starts at layer ',iaNLTEStart(iI)
    END DO
    write (kStdWarn,*) ' >>>>>>>>>>>>>>>>>>>>>>>>>>> '

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
 30   CONTINUE
      IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
        iJ = iJ - 1
        GOTO 30
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
        ! caUAPlanckBloatFile(iG:iG)  = ' '
      END DO
      iJ = 80
 35   CONTINUE
      IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
        iJ = iJ - 1
        GOTO 35
      END IF
      DO iG = 1,iJ
        caOutUABloatFile(iG:iG) = caOutname(iG:iG)
        !caUAPlanckBloatFile(iG:iG) = caOutname(iG:iG)
      END DO
      iG = iJ + 1
      caOutUABloatFile(iG:iG+9)  = '_UA_bloat'
      ! caUAPlanckBloatFile(iG:iG+16) = '_UA_bloat_PLANCK'
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
 40   CONTINUE
      IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
        iJ = iJ - 1
        GOTO 40
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
    end SUBROUTINE nonlteSLOW_LBL

!************************************************************************
! this subroutine is for KCOMPRESSED LBL nonlte
! this reads in the info for the new spectra for ONE gas
    SUBROUTINE nonlteFAST_KCOMP( &
    iRTP,iaGasesIn,iNumGases,raaTemp,raaPress, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks, &
    iaaNLTEChunks,caaStrongLines, &
    iaNLTEBands,raNLTEstart, &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
    raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
    iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
    caOutBloatFile,caOutUABloatFile, &
    caPlanckUAfile,caOutUAfile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iRTP,iaGasesIn(kMaxGas),iNumGases,iSetBloat,iNLTE_SlowORFast
    REAL ::    raaTemp(kProfLayer,kGasStore),raaPress(kProfLayer,kGasStore)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    REAL :: raNLTEstrength(kGasStore)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells the height at which to start NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaaNONLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    REAL :: raNLTEstart(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caOutName
! the solar angle
    REAL :: raKSolarAngle(kMaxAtm)
! tells the name of the GENLN2 ip file that has mixing ratios!
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
! this is pressure levels info
    REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
    INTEGER :: iProfileLayers

! output variable : converts NLTE start heights to AIRS layers
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iErr,iI,iaDumb(kMaxGas),iaGases(kMaxGas)
    REAL :: rH,rHN
    INTEGER :: iStrongISO,iStrongJL,iStrongJU,iStrongGASID,iLTEIn
    INTEGER :: iNLTEStart2350,iJ,iG

    INTEGER :: iBand,iGasID,iNum,iISO,iLineMixBand,iInt,iType,iGas,iDoVoigtChi
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)

    CHARACTER(80) :: caNONLTETempKC
    INTEGER :: iRegr,iMethod,iDefault
    REAL :: raTemp(kProfLayer)
    CHARACTER(80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT)

    CALL NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands,caaaNLTEBands)

    IF (abs(iSetBloat) /= 1) THEN
      write(kStdErr,*) 'need iSetBloat  = -1,+1',iSetBloat
      CALL DOSTOP
    END IF

    IF (iSetBloat == +1) THEN
     write(kStdErr,*) 'you have specified you want kCompressed NLTE database'
      write(kStdErr,*) 'need iSetBloat  = -1 if iNLTE_SlowORFast = -2'
      CALL DoStop
    END IF
    caPlanckBloatFile = 'blankcaPlanckBloatFile'
    caOutBloatFile    = 'blankcaOutBloatFile'
    caOutUABloatFile  = 'blankcaOutUABloatFile'

    IF (abs(iDoUpperAtmNLTE) /= 1) THEN
      write(kStdErr,*) 'need iDoUpperAtmNLTE  = -1,+1',iDoUpperAtmNLTE
      CALL DOSTOP
    END IF
            
    caNONLTETempKC = 'boo_nlte_file'
    DO iI = 1,kProfLayer
      raTemp(iI) = raaTemp(iI,1)   !!assume all kinetic temps equal
    END DO

    IF (iNumNLTEGases > 1) THEN
      write(kStdErr,*) 'kCARTA can only handle ONE NLTE gas ... CO2!!!'
      CALL DoStop
    END IF

    IF (iNumNLTEGases < 1) THEN
      write(kStdErr,*) 'huh, 0 or -ve number of NLTE gases????'
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
    iaDumb = iaGasesIn
    iaGases = -1

    iType = 1
    DO iInt = 1,kMaxGas
      IF (iaDumb(iInt) > 0) THEN
        iaGases(iType) = iInt
        iType = iType + 1
      END IF
    END DO

 5030 FORMAT(A130)

    iNumLinesRead=0
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7)

    DO iI = 1,kGasStore
      IF (abs(raNLTEstrength(iI)) > 10.0) THEN
        write(kStdWarn,*)'gas',iI, 'multiply NLTE data by',raNLTEstrength(iI)
      END IF
    END DO

    iNumLinesRead=1
! read how many gases will have new nonLTE spectroscopy
    iErr=-1
    IF ((iNumNLTEGases < 1) .OR. (iNumNLTEGases > kGasStore)) THEN
      write(kStdErr,*)'need a valid number of gases in *NONLTE!!'
      write(kStdErr,*)'must be > 0, < kGasStore!!'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF

    DO iLTEIn=1,iNumNLTEGases
      IF ((iaNLTEGasID(iLTEIn) < 1) .OR. (iaNLTEGasID(iLTEIn) > kMaxGas)) THEN
        write(kStdErr,*)'iaNLTEGasID(',iLTEIn,') = ',iaNLTEGasID(iLTEIn)
        write(kStdErr,*)'need a valid gas ID in *SPECTRA!!'
        write(kStdErr,*)'please check and retry (subr nonlteFAST_KCOMP)!'
        CALL DoSTOP
      END IF
      IF ((iaNLTEChunks(iLTEIn) < 1) .OR. (iaNLTEChunks(iLTEIn) > kNumkCompT)) THEN
        write(kStdErr,*)'cannot have so many kCARTA chunks in *SPECTRA!!'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP
      END IF
    END DO

! no need to figure out above which height the NLTE starts
! as this is 30 km (already in kCOMP DATABASE
    write(kStdWarn,*) 'kCompressed NLTE database starts NLTE at 30 km '

! now do the output file names, if we need to do the UA and dump out all rads
! and/or optical depths
! right now the code only dumps out rads
    IF (iDoUpperAtmNLTE > 0) THEN
      DO iG = 1,80
        caPlanckUAFile(iG:iG)  = ' '
        caOutUAFile(iG:iG) = ' '
      END DO
      iJ = 80
 40   CONTINUE
      IF ((caOutName(iJ:iJ) == ' ') .AND. (iJ > 1)) THEN
        iJ = iJ - 1
        GOTO 40
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
    end SUBROUTINE nonlteFAST_KCOMP

!************************************************************************
! this subroutine finds out where the Dave Edwards/ Manuel Lopez Puertas
! profiles starts to differs from the input temperature profile
    SUBROUTINE FindNLTEHeight(iLTEIn,caaNLTETemp, &
    iStrongISO,iStrongJL,iStrongJU,iStrongGasID, &
    iaGases,raaTemp,raaPress, &
    raPressLevels,raLayerHeight,iProfileLayers,rH,iBand)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input parameters
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    INTEGER :: iLTEIn
! this is pressure levels info
    REAL :: raPressLevels(kProfLayer+1),raLayerHeight(kProfLayer)
! these are the strongest NLTE band parameters (CO2 at 4 um)
    INTEGER :: iStrongISO,iStrongJL,iStrongJU,iStrongGasID
! these are input gas profiles
    INTEGER :: iaGases(kMaxGas),iProfileLayers
    REAL ::    raaTemp(kProfLayer,kGasStore),raaPress(kProfLayer,kGasStore)
    INTEGER :: iBand
! input/output variable : rH in meters
    REAL :: rH
     
! local variables
    INTEGER :: iL,iG,raTemp(kProfLayer),iLp1,iLp2
    REAL :: rL0,rLp1,rLp2,rDeltaThreshold
    CHARACTER(80) :: caFname                !!! file to read
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
    CHARACTER(100) :: ca1,ca2

    rDeltaThreshold = 0.15   !!need (Tv - Tk) > rDeltaThreshold for NLTE

    write(kStdWarn,*) 'User specified start NLTE height = ',rH
    iStart = (kProfLayer-iProfileLayers+1)

    raPavg(iStart:kProfLayer) = raaPress(iStart:kProfLayer,1)*kAtm2mb

!!!!fill up the 1 -> iProfileLayers-1 indices w/ increasing dummy numbers
    DO iI = iStart-1,1,-1
      raPavg(iI) = raPavg(iI+1) + 10
    END DO

    rH0 = rH
!!! set up the atmosphere profile
    iG = iStrongGasID
           
    rH1 = 0.0
    rH2 = 0.0
    rH1 = sum(raaTemp(iStart:kProfLayer,1))     !!!! should be water temp
    rH2 = sum(raaTemp(iStart:kProfLayer,iG))    !!!! should be CO2 temp
    rH1 = rH1/iProfileLayers
    rH2 = rH2/iProfileLayers
    IF ((rH2 - rH1) < 10.0) THEN
      !!!! only have CO2 in profile
      iG = 1
    END IF
            
    daJU(1) = iStrongJU * 1.0d0
    daJL(1) = iStrongJL * 1.0d0
    caFName = caaNLTETemp(iLTEIn)

    CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL,iStrongGasID,iStrongISO,+1, &
      iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTEtemp1,dVibCenter)

    write(kStdWarn,*) 'iGASID,ISO,iLSGQ,iUSGQ,vibcntr = ',iStrongGasID,iStrongISO,iStrongJL,iStrongJU,dVibCenter

! wap so pressures are increasing
    IF (raTPress1(1) > raTPress1(iNumVibLevels)) THEN
      ! rearrange from small to large
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

    call spl(raTPress1,raTTemp1,iNumVibLevels,raPavg,raLTE,kProfLayer)
    call spl(raTPress1,raNLTEtemp1,iNumVibLevels,raPavg,raNLTE,kProfLayer)

    raDeltaT   = raNLTE - raLTE
    raDeltaUse(1:kProfLayer) = raNLTE(1:kProfLayer) - raaTemp(1:kProfLayer,iG)
    iL = kProfLayer+1

    write(kStdWarn,*) 'start checking LTE/NLTE profile differences : '

    iL = 1
 30 CONTINUE
    IF (raPressLevels(iL) < 1.0) THEN !!!no info here, skip
      iL = iL + 1
      GOTO 30
    END IF
    write(kStdWarn,*) 'surface level at  ',iL,' = ',raPressLevels(iL),' mb'

!!!first few layers might have a temperature inversion
 40 CONTINUE
    IF (raPressLevels(iL) > 600) THEN !!!lower than 4.5 km, skip
      iL = iL + 1
      GOTO 40
    END IF
    write(kStdWarn,*) '4.5 km roughly at ',iL,' = ',raPressLevels(iL),' mb'
    i600 = iL

 60 CONTINUE
    IF ((raDeltaUse(iL) <= rDeltaThreshold) .AND. (iL < kProfLayer)) THEN
      !!! delta(profile) - delta(GENLN2 NLTE) < 1 K
      !!! profile > GENLN2 LTE
      iL = iL + 1
      GOTO 60
    END IF

! check to make sure the next 2 or 3 layers are also in NLTE, ie they also
! have (Tv-Tk) > 0.2
    IF (iL <= (kProfLayer-2)) THEN
      rL0 =  raDeltaUse(iL)
      iLp1 = iL + 1
      rLp1 = raDeltaUse(iLp1)
      iLp2 = iL + 2
      rLp2 = raDeltaUse(iLp2)
      IF ((rL0 <= rDeltaThreshold) .OR. (rLp1 <= rDeltaThreshold) .OR. (rLp2 <= rDeltaThreshold)) THEN
        !!! was a false set; the next couple or so layers are still in LTE;
        iL = iL + 1
        GOTO 60
      END IF
    END IF

    rH = raLayerHeight(iL)
    IF (iL >= kProfLayer-2) rH = 45000.0   !!!default to 45 km
    write(kStdWarn,*) 'start doing NLTE profile at layer ',iL,', p = ',raPressLevels(iL),' mb, h = ', rH/1000,' km'
          
    IF (rH > rH0) THEN
      !!! oops ... let us start where user specified, instead of
      !!! higher up in the atmosphere!
      write(kStdWarn,*) 'User wants NLTE to start at ',rH0/1000,' ,not ',rH/1000,' km'
      rH = rH0
    END IF

    write(kStdWarn,*) 'Start NLTE at (user supplied VS 2350[T-Tvib]) height = ',rH0/1000,' vs ',rH/1000,' km'

    ca1 = 'iI   Pavg    Tk(klayers) | Tk(VibT)          dT  |     Tv           dT '
    ca2 = '-------------------------|-----------------------|----------------------'

    IF (iBand == 1) THEN
      CALL GetUSSTD_2350(raPavg,iStart,daJL,daJU,iaJ_UorL,raLTE_STD,raNLTE_STD)
      write(kStdWarn,*) ca1
      write(kStdWarn,*) ca2
      DO iI = iStart, kProfLayer
        write(kStdWarn,1234) iI,raPavg(iI),raaTemp(iI,1),'|', &
          raLTE(iI), raLTE(iI)-raaTemp(iI,1),'|', &
          raNLTE(iI),raNLTE(iI)-raaTemp(iI,1)
      END DO
    END iF

    IF (iBand == 1) THEN

      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' Comparisons of Tk,tNLTE vs USSTD : '
      ca1 = 'iI   Pavg    Tk(klayers)      TStd         dT  |     Tv      TvSTD      dTv'
      ca2 = '-----------------------------------------------|----------------------------'

      write(kStdWarn,*) ca1
      write(kStdWarn,*) ca2
      DO iI = iStart, kProfLayer
        write(kStdWarn,1250) iI,raPavg(iI),raaTemp(iI,1),raLTE_STD(iI),raaTemp(iI,1)-raLTE_STD(iI),'|', &
           raNLTE(iI), raNLTE_STD(iI),raNLTE(iI)-raNLTE_STD(iI)
      END DO
    END iF

 1234 FORMAT(I3,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1,' ',2(F10.5,' '),A1)
 1250 FORMAT(I3,' ',4(F10.5,' '),A1,' ',3(F10.5,' '))
     
    RETURN
    end SUBROUTINE FindNLTEHeight

!************************************************************************
! simple routine to get US STD 2350 NLTE temps
    SUBROUTINE GetUSSTD_2350(raPavg,iStart,daJL,daJU,iaJ_UorL,raLTE_STD,raNLTE_STD)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    CHARACTER(80) :: caFname                !!! file to read
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN) !!! quantum numbers
    INTEGER :: iaJ_UorL(kHITRAN)           !!! matched upper or lower quant nos
    INTEGER :: iGasID, iISO                !!! what to read
    INTEGER :: iAllOrSome                  !!! read NLTE (+1) or LTE (-1) stuff?
    REAL :: raTTemp1(kNLTEProfLayer),raTPress1(kNLTEProfLayer)
    REAL :: raNLTEtemp1(kNLTEProfLayer),raQtips1(kNLTEProfLayer)
    INTEGER :: iNumVibLevels,iStart,iStrongGasID,iStrongISO
    REAL :: raPavg(kProfLayer),raNLTE_STD(kProfLayer),raLTE_STD(kProfLayer)
    DOUBLE PRECISION :: dVibCenter
    INTEGER :: iI,iJ

    REAL :: raPjunk(kProfLayer),raTJunk(kProfLayer)
    REAL :: raLogTPress1(kNLTEProfLayer),raTX1(kNLTEProfLayer),raJunk(kNLTEProfLayer)
    REAL :: rX,rY,raTX2(kNLTEProfLayer)
    REAL :: rYP1,rYPN,raY2P1(kNLTEProfLayer),raWorkP1(kNLTEProfLayer), &
                      raY2P2(kNLTEProfLayer),raWorkP2(kNLTEProfLayer)    

    caFName = 'xnlte_1_1_1_6_sol_0.genln2'
    CALL concatCA80(caAuxNLTERefsPath,caFName)

    iStrongGasID = 2
    iStrongISO = 1
    CALL ReadGENLN2_NLTE_Profile(caFName,daJL,daJU,iaJ_UorL, &
    iStrongGasID,iStrongISO,+1, &
    iNumVibLevels,raTPress1,raTTemp1,raQtips1,raNLTEtemp1,dVibCenter)

    DO iI = 1,iNumVibLevels
      raLogTPress1(iI) = log(raTPress1(iI))
      ! print *,iI,raLogTPress1(iI),raTPress1(iI),raTTemp1(iI),raNLTEtemp1(iI)
    END DO

    DO iI = 1,iNumVibLevels
      raLogTPress1(iI) = log(raTPress1(iNumVibLevels-iI+1))
      raTX1(iI)        = raNLTEtemp1(iNumVibLevels-iI+1)
      raTX2(iI)        = raTTemp1(iNumVibLevels-iI+1)
      !write(kStdWarn,*),'GetUSSTD_2350 A',iI,raTPress1(iNumVibLevels-iI+1),raTX2(iI),raTX1(iI)
    END DO

!! now interp this onto raPavg
!     Assign values for interpolation
!     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
    rYP1 = 1.0E+16
    rYPN = 1.0E+16
    CALL sply2(raLogTPress1,raTX1,iNumVibLevels,rYP1,rYPN,raY2P1,raWorkP1)
    CALL sply2(raLogTPress1,raTX2,iNumVibLevels,rYP1,rYPN,raY2P2,raWorkP2)    
    
    DO iI = iStart,kProfLayer
      rX = log(raPAvg(iI))
      CALL splin(raLogTPress1,raTX1,raY2P1,iNumVibLevels,rX,rY)
      raNLTE_STD(iI) = rY
      CALL splin(raLogTPress1,raTX2,raY2P2,iNumVibLevels,rX,rY)
      raLTE_STD(iI) = rY
      IF (iI == iStart) then
        write(kStdWarn,*) 115,raLogTPress1(115),raTX2(115),raTX1(115)
        write(kStdWarn,*) iI,rX,raLTE_STD(iI),raNLTE_STD(iI)
        write(kStdWarn,*) 116,raLogTPress1(116),raTX2(116),raTX1(116)
      end if	
    END DO

    RETURN
    end SUBROUTINE GetUSSTD_2350

!************************************************************************
! this subroutine maps the bands in nml file caaaNLTEBandsOrig to files
    SUBROUTINE NLTEBandMapper(iNumNLTEGases,iaNLTEGasID,iaNLTEBands, &
    caaaNLTEBands)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input vars
    INTEGER :: iNumNLTEGases,iaNLTEGasID(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NLTE bands bad boys
    INTEGER :: iaNLTEBands(kGasStore)
! input/output vars
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)

! local vars
    INTEGER :: iI,iJ,iK,iGasID,iISO,iLSGQ,iUSGQ,iUMBC_LBL
    CHARACTER(80) :: caaaNLTEBandsOrig(kGasStore,kNumkCompT),caX,caY
    CHARACTER(4) :: ca4
    INTEGER :: iaBandID(kNumkCompT),iCount,iCountX,iFound

    iaBandID = -1

    caaaNLTEBandsOrig = caaaNLTEBands

    DO iI = 1,iNumNLTEGases
      iCount = 1
      DO iJ = 1,iaNLTEBands(iI)
        caX = caaaNLTEBandsOrig(iI,iJ)
        read(caX,*) iGasID,iIso,iLSGQ,iUSGQ
        IF (iGasID /= 2) THEN
          write(kStdErr,*) 'Need GasID = 2 in nml variable caaaNLTEBand'
          CALL DoStop

          !                                            iso iLSGQ iUSGQ umbc-lbl
          !--------------------------------------------------------------------
        ELSEIF ((iLSGQ == 1) .AND. (iUSGQ == 9)) THEN
          IF (iIso == 1) iUMBC_LBL = 2350          !1    1    9    2350
          IF (iIso == 2) iUMBC_LBL = 2351          !2    1    9    2351
          IF (iIso == 3) iUMBC_LBL = 2352          !3    1    9    2352
          IF (iIso == 4) iUMBC_LBL = 2355          !4    1    9    2355
          IF ((iIso > 4) .OR. (iIso < 1)) THEN
            write(kStdErr,*) ' SUBROUTINE NLTEBandMapper   1,9 isotope error'
            write(kStdErr,*) ' isotope = ',iIso
            CALL DoStop
          END IF

        ELSEIF ((iLSGQ == 4) .AND. (iUSGQ == 24)) THEN
          IF (iIso == 1) iUMBC_LBL = 2310          !1    4   24    2310
          IF (iIso == 2) iUMBC_LBL = 2311          !2    4   24    2311
          IF ((iIso > 2) .OR. (iIso < 1)) THEN
            write(kStdErr,*) ' SUBROUTINE NLTEBandMapper  4,24 isotope error'
            write(kStdErr,*) ' isotope = ',iIso
            CALL DoStop
          END IF

        ELSEIF ((iLSGQ == 2) .AND. (iUSGQ == 16)) THEN
          IF (iIso == 1) iUMBC_LBL = 2320          !1    2   16    2320
          IF (iIso == 2) iUMBC_LBL = 2321          !2    2   16    2321
          IF (iIso == 3) iUMBC_LBL = 2322          !2    2   16    2322
          IF ((iIso > 3) .OR. (iIso < 1)) THEN
            write(kStdErr,*) ' SUBROUTINE NLTEBandMapper  2,16 isotope error'
            write(kStdErr,*) ' isotope = ',iIso
            CALL DoStop
          END IF
                      
        ELSEIF ((iLSGQ == 3) .AND. (iUSGQ == 23)) THEN
          IF (iIso == 1) iUMBC_LBL = 2353          !1    3   23    2353
          IF (iIso == 2) iUMBC_LBL = 2253          !2    3   23    2253
          IF ((iIso > 2) .OR. (iIso < 1)) THEN
            write(kStdErr,*) ' SUBROUTINE NLTEBandMapper  3,23 isotope error'
            write(kStdErr,*) ' isotope = ',iIso
            CALL DoStop
          END IF

        ELSEIF ((iLSGQ == 5) .AND. (iUSGQ == 25)) THEN
          IF (iIso == 1) iUMBC_LBL = 2354          !1    5   25    2354
          IF (iIso == 2) iUMBC_LBL = 2254          !2    5   25    2254
          IF ((iIso > 2) .OR. (iIso < 1)) THEN
            write(kStdErr,*) ' SUBROUTINE NLTEBandMapper  3,23 isotope error'
            write(kStdErr,*) ' isotope = ',iIso
            CALL DoStop
          END IF

        !!! weak hot bands that Manuel LopezPuertas suggested to use
        ELSEIF (iIso == 1) THEN
          IF ((iLSGQ == 2) .AND. (iUSGQ == 15)) iUMBC_LBL = 2110
          IF ((iLSGQ == 3) .AND. (iUSGQ == 25)) iUMBC_LBL = 2120
          IF ((iLSGQ == 4) .AND. (iUSGQ == 22)) iUMBC_LBL = 2130
          IF ((iLSGQ == 5) .AND. (iUSGQ == 23)) iUMBC_LBL = 2140
          IF ((iLSGQ == 3) .AND. (iUSGQ == 22)) iUMBC_LBL = 2150
          IF ((iLSGQ == 6) .AND. (iUSGQ == 36)) iUMBC_LBL = 2160
          IF ((iLSGQ == 7) .AND. (iUSGQ == 37)) iUMBC_LBL = 2170
          IF ((iLSGQ == 8) .AND. (iUSGQ == 38)) iUMBC_LBL = 2180

        ELSE
          write(kStdErr,*) 'trying to figure out CO2 band ....',caX
          write(kStdErr,*) 'NLTEBandMapper cannot decipher',iIso,iLSGQ,iUSGQ
          CALL DoStop
        ENDIF

        ! file name caY are of the form 'blah2350.dat'
         write(ca4,40) iUMBC_LBL
 40      FORMAT(I4)
         DO iK = 1,80
           caY(iK:iK) = ' '
         END DO
         caY(1:4) = ca4
         caY(5:8) = '.dat'

         CALL ConcatCA80(caStrongLineParams,caY)
         caaaNLTEBands(iI,iJ) = caY
         write(kStdWarn,*) iI,iJ,'-> NLTE : iISO,iLSGQ,iUSGQ,UMBC-LBL iD: ',iISO,iLSGQ,iUSGQ,iUMBC_LBL

         IF (iCount == 1) THEN
           iaBandID(iCount) = iUMBC_LBL
           iCount = iCount + 1
         ELSE
           !!! check to make sure this band_id is not used already
           iFound = -1
           DO iCountX = 1,iCount-1
             IF (iaBandID(iCountX) == iUMBC_LBL) THEN
               iFound = +1
               write(kStdErr,*) 'NLTEBandMapper found ',iCountX,iUMBC_LBL
             END IF
           END DO
           IF (iFound == +1) THEN
             write(kStdErr,*) 'NLTEBandMapper found double count iUMBC_LBL'
             CALL DoStop
           ELSE
             iaBandID(iCount) = iUMBC_LBL
             iCount = iCount + 1
           END IF
         ENDIF
       END DO
    END DO

    RETURN
    end SUBROUTINE NLTEBandMapper

!************************************************************************
! subroutine adds in the molgases
    SUBROUTINE add_molgas(iNgas,iaGasesNL,iaMOLgases)

    IMPLICIT NONE

    INCLUDE '../INCLUDE/kcartaparam.f90'

! input/output
    INTEGER :: iNgas
! output
    INTEGER :: iaGasesNL(kGasComp),iaMOLgases(kMaxGas)

! local
    INTEGER :: iMax,iC,iCC,iaInDataBase(kGasComp),iTag,iErr,iaTemp(kMaxGas)
    INTEGER :: iNotMainGas,iNgasesCheck,iInt,iWhichKC,iYesWV,i101,i102,i103
    INTEGER :: iKLBLRTMgases,iMax_molgas

    iWHichKC = iNgas
    iKLBLRTMgases = 32
          
    IF (iWhichKC == -1) THEN
      !! this is kc1.15 or later
      ! use ALL gases in the compressed database for v1.15+
      !  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9
      !         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,28,X,X,31,32,...42
      !         PLUS kSelf,kFor,kHeavyWater (101,102,103)
      iMax = kGasComp
      write(kStdWarn,*) 'including all lbl gases from 1 to ', kGasComp
      write(kStdWarn,*) ' plus 101,102,103 for Self,Forn cont, heavy water'
    ELSEIF ((iWhichKC == -5) .OR. (iWhichKC == -6)) THEN
      !! this is kc1.14 or earlier
      !! also LBLRTM style
      ! use ALL gases in the compressed database for c1.14-
      !  check to see the following reference profiles exist : 1,2,3,4,5,6,7,8,9
      !         10,11,12,X,X,15,16,X,18,19,20,21,22,23,X,25,26,27,28,X,X,31
      !         PLUS kSelf,kFor,kHeavyWater (101,102,103)
      !        iMax = 31
      !        write(kStdWarn,*) 'including v114- lbl gases from 1 to ',iMax
      iMax = iKLBLRTMgases
      write(kStdWarn,*) 'including LBLRTM gases from 1 to ',iMax
      write(kStdWarn,*) ' plus 101,102,103 for Self,Forn cont, heavy water'
    ELSEIF (iWHichKC > 0) THEN
      write(kStdWarn,*) 'using user sepcified list'
    ELSE
      write(kStdErr,*) 'need iWhichKC = -1 or -5/-6 or > 0'
      CALL DoStop
    END IF

    IF (iWhichKC < 0) THEN
      iNgas = 0
      iF (iWhichKC == -1) THEN
        iMax_molgas = kGasComp
      ELSEIF ((iWhichKC == -5) .OR. (iWhichKC == -6)) THEN
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
          write(kStdWarn,*)'Including gasID ',iC,' in comp database'
          iaTemp(iCC) = iC
          iCC = iCC+1
        END IF
      END DO

      IF (kCKD >= 0) THEN
        DO iC = kNewGasLo,kNewGasHi
          write(kStdWarn,*)'Including gasID ',iC,' in comp database'
          iaTemp(iCC) = iC
          iCC = iCC+1
          iNgas = iNgas+1
        END DO
      END IF

      !!!heavy water
      iC = kNewGasHi+1
      write(kStdWarn,*)'Including gasID ',iC,' in comp database'
      iaTemp(iCC) = iC
      iCC = iCC+1
      iNgas = iNgas+1

    ELSEIF (iWhichKC > 0) THEN
      iYesWV = -1
      i101 = -1
      i102 = -1
      i103 = -1
      DO iC = 1,iWhichKC
        iaTemp(iC) = iaGasesNL(iC)
        IF (iaTemp(iC) == 1) iYesWV = +1
        IF (iaTemp(iC) == 101) i101 = +1
        IF (iaTemp(iC) == 102) i102 = +1
        IF (iaTemp(iC) == 103) i103 = +1
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
 777 FORMAT('Error in MOLGAS!! iC = ',I2,' found Gas ID   =  ',I2,'( &
            MOLGAS ID''s should be between ',I2,' and ',I2 ,')')
    DO iC = 1,iNgas
      iNotMainGas = MainGas(iaTemp(iC))
      IF (iNotMainGas < 0) THEN
        ! gas does not exist in the compressed base ... stop
        WRITE(kStdErr,777) iC,iaTemp(iC),1,kGasComp
        CALL DoSTOP
      END IF
      ! check to see if kcomp files do exist
      iTag = -1
      IF (iaTemp(IC) <= kGasComp) THEN
        iCC = iCheckCompDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
        IF (iCC < 0) THEN
          WRITE(kStdWarn,780) iaTemp(iC)
          END IF
        ELSEIF (iaTemp(IC) <= kNewGasHi) THEN
          WRITE(kStdWarn,785) iaTemp(iC)
        ELSEIF (iaTemp(IC) == kNewGasHi+1) THEN
          WRITE(kStdWarn,787) iaTemp(iC)
        END IF
    END DO
 780 FORMAT('Warning! GasID ',I3,' not in compressed data base')
 785 FORMAT('Warning! GasID ',I3,' is a new continuum gas')
 787 FORMAT('Warning! GasID ',I3,' is heavy water (gasID  =  1, HITRAN isotope  =  4)')
          
    iNgasesCheck = iNgas

! set the identities of the gases whose abs spectra are in files
! iaMOLgases keeps track of how many times the GasID has been found, so that
! no double counting of the gases (between MOLGAS and XSCGAS) is done
    DO iInt = 1,iNgas
      iaMOLgases(iaTemp(iInt)) = iaMOLgases(iaTemp(iInt))+1
    END DO

! check to see that the same gas ID has not been entered more than once
    DO iInt = 1,iNgas
      IF ((iaMOLgases(iaTemp(iInt)) > 1) .OR. (iaMOLgases(iaTemp(iInt)) < 0)) THEN
        write(kStdErr,*) 'Gas ID',iaTemp(iInt),' entered more than'
        write(kStdErr,*)'once. Please check *MOLGAS and retry'
        CALL DoSTOP
      END IF
    END DO

    iNgas = iNgasesCheck
    write(kStdWarn,*) 'MOLGAS ... gases stored are ',iNgas
    DO iInt = 1,kGasComp
      IF (iaMOLgases(iInt) > 0) THEN
        write(kStdWarn,*) '     going to use compressed database gas ',iInt
      END IF
    END DO

    DO iInt = kNewGasLo,kNewGasHi
      IF (iaMOLgases(iInt) > 0) THEN
        write(kStdWarn,*) '     going to use new continuum gas ',iInt
      END IF
    END DO

    iInt = kNewGasHi+1
    IF (iaMOLgases(iInt) > 0) THEN
      write(kStdWarn,*) '     going to use new heavy water gas ',iInt
    END IF

    IF (kCKD >= 0) THEN
      DO iInt = kNewGasLo,kNewGasHi
        IF (iaMOLgases(iInt) <= 0) THEN
          write(kStdWarn,*) 'Cannot have CKD on and gasIDs 101/102 unused'
          write(kStdErr,*) 'Cannot have CKD on and gasIDs 101/102 unused'
          Call DoSTOP
        END IF
      END DO
    END IF

    RETURN
    end SUBROUTINE add_molgas

!************************************************************************
! subroutine adds in the xsecgases
    SUBROUTINE add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)

    IMPLICIT NONE

    INCLUDE '../INCLUDE/kcartaparam.f90'

! input/output
    INTEGER :: iNXsec
    INTEGER :: iaLXsecNL(kGasXSecHi-kGasXSecLo+1),iaXSCgases(kMaxGas)

! local
    INTEGER :: iInt,iErr,iC,iCC,iNXsecCheck
    INTEGER :: iaTemp(kMaxGas),iTag
    INTEGER :: iaInDataBase(kMaxLayer),iWhichXSC,iKLBLRTMgases

    iWhichXSC = iNXsec

    IF (iWhichXSC == -1) THEN
      write(kStdWarn,*) 'including all xsc gases from ',kGasXsecLo, ' to ', kGasXsecHi
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
          write(kStdWarn,*) 'Including gasID ',iC,' in xsec database'
          iaTemp(iCC) = iC
          iCC = iCC+1
        END IF
      END DO

    ELSE IF ((iWhichXSC == -5) .OR. (iWhichXSC == -6))THEN
      iKLBLRTMgases = 63
      write(kStdWarn,*) 'including LBLRTM xsc gases from ',kGasXsecLo,' to ',iKLBLRTMgases
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
          write(kStdWarn,*) 'Including gasID ',iC,' in xsec database'
          iaTemp(iCC) = iC
          iCC = iCC+1
        END IF
      END DO

    ELSEIF (iWhichXSC > 0) THEN
      DO iC = 1,iNXsec
        iaTemp(iC) = iaLXsecNL(iC)
      END DO
    ELSE
      write(kStdErr,*) 'xscgas4 can only process iNXsec = -5/-6,-1, or > 0'
      CALL DOSTOP
    END IF
     
! check the molecular ID's are between kGasXsecLo and kXsecGasHi
! (should be in the cross sec data base)
! if iNXsec = -1, everything should have been set correctly above; if user
! enetered in the GasIDs him/her self, there could be mistakes
  777 FORMAT('Error in XSCGAS!! found Gas ID  = ',I2,'( &
            XSCGAS ID''s should be between ',I2,' and ',I2 ,')')
  780 FORMAT('Warning! GasID ',I2,' not in xsec data base')
    DO iC = 1,iNXsec
      IF ((iaTemp(iC) < kGasXsecLo) .OR. (iaTemp(iC) > kGasXsecHi)) THEN
        ! gas does not exist in the xsec data base ... stop
        WRITE(kStdErr,777) iaTemp(iC),kGasXsecLo,kGasXsecHi
        CALL DoSTOP
      END IF
      ! check to see if cross section data files for this gas does exist
      iTag = -1
      iCC = iCheckXsecDataBase(iaTemp(iC),-100.0,-100.0,iTag,iErr)
      IF (iCC < 0) THEN
        WRITE(kStdWarn,780) iaTemp(iC)
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
      IF ((iaXSCgases(iaTemp(iInt)) > 1) .OR. (iaXSCgases(iaTemp(iInt)) < 0)) THEN
        write(kStdErr,*) 'Gas ID ',iaTemp(iInt),' entered more than'
        write(kStdErr,*) 'once Please check *XSCGAS and retry'
        CALL DoSTOP
      END IF
    END DO

    iNXsec = iNXsecCheck
    write(kStdWarn,*)'cross section gases stored are ',iNXsec
    DO iInt = kGasXsecLo,kGasXsecHi
      IF (iaXSCgases(iInt) > 0) THEN
        write(kStdWarn,*) 'CrossSect database gas ',iInt
      END IF
    END DO

    RETURN
    end SUBROUTINE add_xsecgas
!************************************************************************

! this subroutine checks to see if the gasID is 1-36 or 101-102 or 103
    INTEGER FUNCTION MainGas(iGasID)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iGasID

    INTEGER :: iI,i1,i2
           
    iI = -1        !assume this is a main gas (1-36, or 101-102)

    i1 = -1
    i2 = -1

    IF ((iGasID >= 1) .AND. (iGasID <= kGasComp)) THEN
      i1 = 1          !! main gas
    ELSEIF (iGasID == kNewGasHi+1) THEN
      i1 = 1          !! heavy water
    ELSEIF ((iGasID >= kNewGasLo) .AND. (iGasID <= kNewGasHi)) THEN
      i2 = 1          !! water continuum
    END IF
          
    IF ((i1 == 1) .OR. (i2 == 1)) THEN
      iI = 1
    END IF

    IF ((i2 == 1) .AND. kCKD < 0) THEN
      write(kStdWarn,*) 'Cannot have gases 101,102 with CKD turned off'
      write(kStdErr,*) 'Cannot have gases 101,102 with CKD turned off'
      Call DoSTOP
    END IF

    MainGas = iI

    RETURN
    end FUNCTION MainGas

!************************************************************************
END MODULE n_gas_wt_spectra
