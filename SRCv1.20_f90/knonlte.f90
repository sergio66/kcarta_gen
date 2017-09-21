! Copyright 2015
! University of Maryland Baltimore County
! All Rights Reserved

MODULE knonlte

USE basic_common
USE spline_and_sort
use klineshapes
use kvoigt_cousin
use klinemix
use kreadVTprofiles
use kbloat
use kcousin
use kpredictVT
use kcoeffMAIN
USE kcoeff_basic

IMPLICIT NONE

CONTAINS

! to put in arbitrary CO2 amounts, have to be a little more careful!!!!!!!
!   load in the LTE weak+strong database       (at multiplier = 1)
!   compute weak backgnd LTE, lower atm        (at multiplier = rX)
!   compute strong NLTE lower atm              (at multiplier = rX)
!   compute strong NLTE radiance enhance, LA   (at multiplier = rX)
!   compute weak backgnd LTE, upper atm        (at multiplier = rX)
!   compute strong NLTE upper atm              (at multiplier = rX)
!   compute strong NLTE radiance enhance, UA   (at multiplier = rX)
!   combine and output everything              (at multiplier = rX)
!   divide the OD output by multiplier, so that the raaMix code recorrects
!                                       the OD

! note this is separate from  raNLTEstrength(1) :
!   raNLTEstrength(1) = +1.0 : says use a weight of 1.0 for all the
!                          upper layer optical depths computed for this gas.\\
!                          Together with parameter kCousin_CO2path, this lets
!                          kCARTA do cousin lineshape computations via the
!                          following fudge : if raNLTEstrength(X) = -1.1212 and
!                          iaNLTEGasID(X) = 2, kCARTA not do NLTE but instead
!                          swaps the linemix uncompressed spectra for cousin
!                          uncompressed spectra in the layers above that
!                          specified by raNONLTEstart(X)\\
!                          Otherwise raNLTEstrength is really not used

!************************************************************************
!************************************************************************
!************************************************************************
!  UA pressure levels are from running KCARTADATA/NLTE/UA/read_glatm.m
!************************************************************************
!************************************************************************
!************************************************************************
!* defaults used when producing the NLTE database in Dec 2004
!*                  for NLTE
!*  knonlte.f           iTalk = -1              short summary of NLTE info
!*                                                /band/layer
!*  klineshapes.f       iTalk = -1              short summary of NLTE info
!*                                                /band/layer
!*  kvoigt_cousin.f     iTestGenln2 = -1        use line mix for strong bands
!*  knonlte.f           iNoPressureShiftCO2 = 1 no pressure shift in CO2 lines
!*  klineshapes.f       iQtips_H98 = +1         QTIPS values for partition fcn
!*  kbloat.f            iLinearOrSpline = -1    Use spline to go from 0.0025 to
!*                                                0.0005 spacing, for LTE stuff
!*                                                in BloatCoeffs()
!*
!************************************************************************
!************************************************************************
!        !!! (a) all weak LTE lines have been read in from the compressed
!        !!!     database, computed using a Cousin lineshape
!        !!!     CALL lte_spectra(...)
!        !!!     CALL read_std_optdepths_new(iGasID,raFreq,iStart,
!        !!! (b) these medium LTE lines will have a Cousin lineshape
!        !!!     CALL lte_spectra(...)
!        !!! (c) the strong NLTE lines will have a Cousin lineshape deep in the
!        !!!     4 um band, or a linemix lineshape in the R branchhead
!        !!!     CALL nonlte_spectra(...)
! this means there is a discontinuity at 2355 cm-1, as the strongest sigsig
! band is switched from Cousin to line mix (about a 1 K step)
!************************************************************************
! Ref : Edwards et al : NON LTE Studies, JGR vol 98, pgs 14955 August 1993
! from GENLN4, zpthadj.f ... we are sending in v_upper IUSGQ
! from GENLN4, ntepro.f  ... we are sending in v_upper IUSGQ
! from GENLN4, emissn.f  ... check radiance computation emfac vs abfac
! Ref : KOPRA documentation
!! http://www-imk.fzk.de/imk2/ame/publications/kopra_docu/
! see pg 277,435 of Lopez-Puertas book for the NLTE bands in the 4.3 um
! fundamental 00011 -> 00001    vl,vu = 1,9     2349 cm-1   hit2350
! isotope     00011 -> 00001    vl,vu = 1,9     2283 cm-1   hit2351
!! FH         01111 -> 01101    vl,vu = 2,16    2336 cm-1   hit2320
!! SH         10012 -> 10002    vl,vu = 3,23    2327 cm-1   hit2353
!! hot 1      02211 -> 02201    vl,vu = 4,24    2324 cm-1   hit2310
!! hot 2      10011 -> 10001    vl,vu = 5,25    2326 cm-1   hit2354
! Ref : HITRAN bd_vibs.for tells you the HITRAN spectroscopic notations
! see SUBROUTINE NLTEBandMapper for all this
! so v_lower == 1, v_upper ==  9 for strong SigSig band
!!          00001         00011
!!/asl/data/hitran/HITRAN2k/BD_Vibs.for

! /home/sergio/KCARTA/SRC/NONLTE2/sergio/NONLTE_MATLAB/make_lbl2kc_weakco2.m
! has topts.linemixloop = -1 SO WE NEED TO ADD on all the linemix bnds!
! which I define as 2310,2311,  2320,2321,2322,   2350,2351,2352,2353,2354,2355
! Pretty much bands 2350,2351,2352,2355,  2320,2321,  2310,2311 are in NLTE
! So we need to include the effects of bands 2322, 2353 and 2354!!!!!!!
! These are linemix bands, so just use Cousin for them!

! NOTE (1) :
!      INTEGER iEdwards  !!!**** choose how to do things!!!!!! for the
!                        !!!**** abs coeff and Planck coeff factors!!!!
! this is following Eqn 15,16 of the NLTE paper by Edwards
! Ref : Edwards et al : NON LTE Studies, JGR vol 98, pgs 14955 August 1993
!      iEdwards = +1
! this is following the TES algorithm
!      iEdwards = -1

! places I can see where it could differ from GENLN are in computing downwell
! background radiation, kCARTA does NOT use raaPlanckCoeff!

! MAGIC NUMBERS
!    iOneCmFine = 400                 !!!! 1 cm-1 = 400 pts at 0.0025 cm-1
!    iOffSet = (iWideMeshLoop-1)*2000 !!!! 1-400 on kCARTA =1-2000 high res

! **************************** testing the code ************************
! 1) testing NLTE, only for strongest 4.3 um band (sig sig)
!   (a) in kcartaparam.f90, use kProfLayer = 43
!   (b) use edwards_jan2003.nml (which is in NLTE dir)
!   (c) test middle of band, since he used CO2 continuum lineshape which means
!       that in the wings, he has absorption(continuum) + abs(sigsig*cousin)
!       compared to me, who only has abs(sigsig*cousin); BUT his BT in between
!       lines is BIGGER than mine! (in between lines)
!   (d) can test the default kCARTA spacing
!          dLineStrenMin = -1.0d0     !!!! use all backgnd lines
!          dDeltaFreqNLTE = -1.0d0    !!!use default dkFreqTag
!          kBoxCarUse = 5             !!!use default 5
!        OR the GENLN4 spacing by setting switches at kcartabasic.f
!          dLineStrenMin = 1.0d-25    !!!! use backgnd lines bigger than this
!          dDeltaFreqNLTE = 1.0d-3    !!!use GENLN4
!          kBoxCarUse = 1             !!!use GENLN4

!************************************************************************
! 6/02/03        Include first order linemixing code into the NLTE .. for the
!                strong bands (2350,2351), use linemix plus birnbaum; else use
!                cousin lineshape (for the "other strong bands"); for medium
!                and weak lines, use plain voigt
!                  Subroutine read_lineparameters initially sets this info
!                  Subroutine CousinVsMix        reassigns individual lines
!                                                to linemix, or back to cousin
!                                                (since eg 2350 P branch sux!!)
!                iLineMixBand = -1               plain voigt
!                iLineMixBand = +1               plain voigt times cousin
!                iLineMixBand = +2               plain voigt times linemix
!                The linemix code is slow, and so the code mainly does Cousin
!                only, by resetting iLineMixBand = 1 in sub read_lineparameters
!                Code can also use fudge files
!                   /asl/data/kcarta/KCARTADATA/General/ChiFile/
!                         nonlte_co2_fudge_2380.txt;
!                  see NONLTE2/sergio/SECOND_TEST/save_ratio_k.m for details

!************************************************************************
! this is the main driver!!!!
! see if current gas ID needs nonLTE spectroscopy
    SUBROUTINE NLTEDriver( &
    iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
    iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
    iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
    iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
    iUseWeakBackGnd,raFreq,pProf,iaCont,rSolzen, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
    pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    raVertTemp,iVertTempSet, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    raUpperPress,raUpperPartPress,raUpperTemp, &
    raUpperGasAmt,raUpperNLTETemp, &
    iUpper,iDoUpperAtmNLTE, &
    dLineStrenMin,dDeltaFreqNLTE, &
    caaUpperMixRatio,iNumberUA_NLTEOut, &
    rFreqStart,rFreqEnd,rFileStartFr, &
    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
    caOutUAFile,caOutUABloatFile, &
    iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
    daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
    daaPlanckCoeffBloat, &
    daaUpperPlanckCoeff, &
    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    iDoDQ,daaDT,daaDQ, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
! rCO2mult  = tells how much input .rtp file wants the CO2 to be multiplied by
!             remeber this is as compared to kCO2ppmv
    REAL :: rCO2mult
! iaNLTEBands     tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
! iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
    INTEGER :: iaNLTEBands(kGasStore)
    INTEGER :: iGas,iaGases(kMaxGas),iSplineType,iUseWeakBackGnd
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    INTEGER :: iNumNLTEGases,iSetBloat
    INTEGER :: iTag,iActualTag,iProfileLayers,iL_low,iL_high
    INTEGER :: iVertTempSet,iaCont(kMaxGas)
    INTEGER :: iUpper,iNumberUA_NLTEOut
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iNLTE_SlowORFast
    REAL :: rFileStartFr
! caaaNLTEBands     tells the name of the files containing the line parameters
! caaNLTETemp       tells the name of the files containing the nonLTE temps
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaStrongLines(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    CHARACTER(80) :: caOutUAFile,caOutUABloatFile
! pProf is the avg layer pressure
    REAL :: pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
    REAL :: raLayerHeight(kProfLayer),rSolzen
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
    REAL :: raNLTEstrength(kGasStore),raFreq(kMaxPts)
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    REAL :: raVertTemp(kProfLayer)
! stuff for the upper levels
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    REAL :: rFreqStart,rFreqEnd
! strengths etc
    DOUBLE PRECISION :: dLineStrenMin,dDeltaFreqNLTE

! output params
! iLTEIn            is just a variable that tells us if nonLTE to be used
    INTEGER :: iLTEIn,iWhichChunk,iFunnyCousin
! tells the nonscattering code, at which layer to start NONLTE rad transfer
    INTEGER :: iNLTEStart
! daaGasAbCoeff has the uncompressed gas absorption coeff
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer),raaCO2_LTE(kMaxPts,kProfLayer)
! daaNLTEGasAbCoeff has the nonLTE gas absorption coeff
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
! this has the bloated stuff
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
! this is ABOVE the 0.005 mb kcarta TOA
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
! these are for Matlab style kCOmp Corner Weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer)
    INTEGER :: iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer)
    REAL ::    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer)
    REAL ::    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer)
    INTEGER :: iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer)
    REAL ::    raQ21(kProfLayer),raQ22(kProfLayer)

    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iDoDQ,iCount,iRefLayer

    INTEGER :: iL

    iLTEIn = -1
          
    IF (iNLTE_SlowORFast == +1) THEN
        write(kStdWarn,*) 'doing slow LBL NLTE ...'
        CALL NLTE_SLOW_LBL( &
        iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
        iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
        iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
        iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
        iUseWeakBackGnd,raFreq,pProf, &
        iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
        pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
        pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
        raRAmt,raRTemp,raRPress,raRPartPress, &
        raVertTemp,iVertTempSet, &
        raTAmt,raTTemp,raTPress,raTPartPress, &
        raUpperPress,raUpperPartPress,raUpperTemp, &
        raUpperGasAmt,raUpperNLTETemp, &
        iUpper,iDoUpperAtmNLTE, &
        dLineStrenMin,dDeltaFreqNLTE, &
        caaUpperMixRatio,iNumberUA_NLTEOut, &
        rFreqStart,rFreqEnd,rFileStartFr, &
        iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
        caOutUAFile,caOutUABloatFile, &
        iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
        daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
        daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
        daFreqBloat, &
        daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
        daaPlanckCoeffBloat, &
        daaUpperPlanckCoeff, &
        daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
        daaUpperPlanckCoeffBloat, &
        daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat)
    ELSEIF (iNLTE_SlowORFast == -2) THEN
        write(kStdWarn,*) 'doing fast kcompressed NLTE ...'
        CALL NLTE_Fast_Compressed( &
        iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
        iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
        iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
        iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
        iUseWeakBackGnd,raFreq,pProf,iaCont,rSolzen, &
        iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
        pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
        pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
        raRAmt,raRTemp,raRPress,raRPartPress, &
        raVertTemp,iVertTempSet, &
        raTAmt,raTTemp,raTPress,raTPartPress, &
        raUpperPress,raUpperPartPress,raUpperTemp, &
        raUpperGasAmt,raUpperNLTETemp, &
        iUpper,iDoUpperAtmNLTE, &
        dLineStrenMin,dDeltaFreqNLTE, &
        caaUpperMixRatio,iNumberUA_NLTEOut, &
        rFreqStart,rFreqEnd,rFileStartFr, &
        iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
        caOutUAFile,caOutUABloatFile, &
        iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
        daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
        daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
        daFreqBloat, &
        daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
        daaPlanckCoeffBloat, &
        daaUpperPlanckCoeff, &
        daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
        daaUpperPlanckCoeffBloat, &
        daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
        iDoDQ,daaDT,daaDQ, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    ELSE
        write(kStdErr,*) 'hmm NLTE Driver neediNLTE_SlowORFast = +1,-2'
        write(kStdErr,*) 'not ',iNLTE_SlowORFast
        CALL DoStop
    ENDIF

    RETURN
    end SUBROUTINE NLTEDriver

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! fast kCOMP calcs for NLTE ODs and Planck Coeffs
    SUBROUTINE NLTE_Fast_Compressed( &
    iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
    iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
    iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
    iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
    iUseWeakBackGnd,raFreq,pProf,iaCont,rSolzen, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
    pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    raVertTemp,iVertTempSet, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    raUpperPress,raUpperPartPress,raUpperTemp, &
    raUpperGasAmt,raUpperNLTETemp, &
    iUpper,iDoUpperAtmNLTE, &
    dLineStrenMin,dDeltaFreqNLTE, &
    caaUpperMixRatio,iNumberUA_NLTEOut, &
    rFreqStart,rFreqEnd,rFileStartFr, &
    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
    caOutUAFile,caOutUABloatFile, &
    iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
    daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
    daaPlanckCoeffBloat, &
    daaUpperPlanckCoeff, &
    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    iDoDQ,daaDT,daaDQ, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
! rCO2mult  = tells how much input .rtp file wants the CO2 to be multiplied by
!             remeber this is as compared to kCO2ppmv
    REAL :: rCO2mult
! iaNLTEBands     tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
! iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
    INTEGER :: iaNLTEBands(kGasStore)
    INTEGER :: iGas,iaGases(kMaxGas),iSplineType,iUseWeakBackGnd
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    INTEGER :: iNumNLTEGases,iSetBloat
    INTEGER :: iTag,iActualTag,iProfileLayers,iL_low,iL_high
    INTEGER :: iVertTempSet,iaCont(kMaxGas)
    INTEGER :: iUpper,iNumberUA_NLTEOut
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iNLTE_SlowORFast
    REAL :: rFileStartFr,rSolzen
! caaaNLTEBands     tells the name of the files containing the line parameters
! caaNLTETemp       tells the name of the files containing the nonLTE temps
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaStrongLines(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    CHARACTER(80) :: caOutUAFile,caOutUABloatFile
! pProf is the avg layer pressure
    REAL :: pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
    REAL :: raNLTEstrength(kGasStore),raFreq(kMaxPts)
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    REAL :: raVertTemp(kProfLayer)
! stuff for the upper levels
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    REAL :: rFreqStart,rFreqEnd
! strengths etc
    DOUBLE PRECISION :: dLineStrenMin,dDeltaFreqNLTE

! output params
! iLTEIn            is just a variable that tells us if nonLTE to be used
    INTEGER :: iLTEIn,iWhichChunk,iFunnyCousin
! tells the nonscattering code, at which layer to start NONLTE rad transfer
    INTEGER :: iNLTEStart
! daaGasAbCoeff has the uncompressed gas absorption coeff
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer),raaCO2_LTE(kMaxPts,kProfLayer)
! daaNLTEGasAbCoeff has the nonLTE gas absorption coeff
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
! this has the bloated stuff
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
! this is ABOVE the 0.005 mb kcarta TOA
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

    REAL :: raUpperPartPress_Std(kProfLayer)

! these are for Matlab style kCOmp Corner Weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer)
    INTEGER :: iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer)
    REAL ::    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer)
    REAL ::    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer)
    INTEGER :: iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer)
    REAL ::    raQ21(kProfLayer),raQ22(kProfLayer)

    DOUBLE PRECISION :: daaWeakOptDepth(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDTp(kMaxPtsJac,kProfLayerJac)   !! planck coeffs d/dT nonzero
    DOUBLE PRECISION :: daaDQp(kMaxPtsJac,kProfLayerJac)   !! planck coeffs d/dq    zero

    INTEGER :: iDoDQ,iZeroPlanck

    INTEGER :: iErr,iCount,iDoNLTE
    INTEGER :: i_NLTEFile_TYPE,iL,iFr
    REAL :: kFrStep,rSolzenX,rSolzenY,rJunk

    INTEGER :: iBand,iType,iUpperStd_Num,iNumMixRatioLevs,iIOUN
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer),rX,rY
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    REAL :: raInterpTempUA(kProfLayer),raUpperPress1013(kProfLayer)
    CHARACTER(80) :: caOutName
    REAL :: raX(kMaxPts),rMult,rMult0
    REAL :: raRPressX(kMaxLayer),raRPPressX(kMaxLayer)
    REAL :: raRAmtx(kMaxLayer),raRTempx(kMaxLayer)

! junk to read in HITRAN
    INTEGER :: iNum,iISO,iLineMixBand,iGasID,iDoVoigtChi
    REAL :: raNLTEtemp(kProfLayer),raVibQFT(kProfLayer),raLTETemp(kProfLayer)
    DOUBLE PRECISION :: daPlanck(kMaxPts)
    DOUBLE PRECISION :: daK(kMaxPts)
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daChi(kMaxPts),daChiBloat(kBloatPts)
    DOUBLE PRECISION :: dVibCenter      !!!from D. Edwards NLTE files
    INTEGER :: iStartUse,iaJ_UorL(kHITRAN)
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)
              
    INTEGER :: iJunkNum,iaJunk(kGasStore)

! see if current gas ID needs nonLTE spectroscopy
    iLTEIn      = -1
    iWhichChunk = -1
    iDoNLTE     = -1
    iJunkNum    = -1

    iLTEIn = OutsideSpectra(iaGases(iGas),iNumNLTEGases,iaNLTEGasID,iJunkNum,iaJunk,raFreq(1),605.0,2830.0,20)

    IF (iLTEIn > 0) THEN
        CALL LowerAtmNLTERefs(raRPressX,raRPPressX,raRTempx,raRAmtx)
    !        DO iL = 1, kProfLayer
    !          print *,iL,raRTemp(iL)-raRTempx(iL),raRAmt(iL)/raRAmtX(iL)
    !        END DO
    !      call dostop
    END IF

    iBand = 2350
    iBand = 1
    IF ((iDoUpperAtmNLTE > 0) .AND. (iLTEIn > 0)) THEN
    !! read in alt/press/temp/upper mix ratio from caaUpperMixRatio
    !!  into raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs
    !!
    !! this is DIFFERENT from the UA NLTE compressed database!
    !! which is read in and saved into raUpper*_Std(iI=1,iiUpperStd_Num)
        CALL MixRatio(iaGases(iGas),rCO2mult,iLTEIn,caaUpperMixRatio, &
        raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
        raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    !! this reads in the GENLN2 style actual NLTE vib temp profiles
    !! for the diffferent CO2 bands, as stored in caaNLTETemp
    !! results stored raUpper*(iI=1,iUpper)
        CALL read_upperatm_lte_temperature( &
        iaGases(iGas),iNLTEStart,iLTEin,iBand,caaNLTETemp, &
        raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        pProf,raPresslevels,raLayerHeight,raThickness, &
        iUpper,raUpperTemp,raUpperGasAmt,raUpperPress,raUpperPartPress, &
        raUpperPressLevels,raUpperThickness)

    !        DO iL = 1,iNumMixRatioLevs
    !          print *,iL,raUpper_Pres(iL),raUpper_MixRatio(iL)
    !        END DO
    !        print *,' '
    !        DO iL = 1,iUpperStd_Num
    !          print *,iL,raUpperPress_std(iL),raUpperMixRatio_std(iL)
    !        END DO
    !        print *,' '
    !        DO iL = 1,iUpper
    !          print *,iL,raUpperPress(iL),raUpperTemp(iL)
    !        END DO
    !        CALL DOSTOP
         
        IF (iDoUpperAtmNLTE > 0) THEN
            rMult0 = 0.8      !!        modify the Planck coeffs
            rMult0 = 0.7275   !!        modify the Planck coeffs
            rMult0 = 0.775    !!        modify the Planck coeffs
            rMult0 = 1.0000   !! do not modify the Planck coeffs
            rMult0 = 0.75     !!        modify the Planck coeffs

        ! need to interp1 (raUpperPress,raUpperTemp) onto (raUpperPress_Std)
            DO iL = 1,iUpper
                raUpperPress1013(iL) = raUpperPress(iL)*1013.25
            END DO

            iZeroPlanck = iUpperStd_Num + 1
            IF (raUpperPress_Std(iUpperStd_Num) < raUpperPress1013(iUpper))THEN
            !! need to turn off the Planck function else we may get bad results
                iZeroPlanck = 1
                10 CONTINUE
                IF (raUpperPress_Std(iZeroPlanck) > raUpperPress1013(iUpper) .AND. &
                iZeroPlanck < iUpperStd_Num) THEN
                    iZeroPlanck = iZeroPlanck + 1
                    GOTO 10
                END IF
            END IF

            IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
                write(kStdWarn,*) 'wil reduce UA Planckcoeffs in successive layers by factor of ',rMult0
                write(kStdWarn,*) 'User Profile : n,minP(TOA_UA) = ',iUpper,raUpperPress1013(iUpper)
                write(kStdWarn,*) 'US Std Prof  : n,minP(TOA_UA) = ',iUpperStd_Num,raUpperPress_Std(iUpperStd_Num)
                write(kStdWarn,*) ' '
                IF (iZeroPlanck < (iUpperStd_Num + 1)) THEN
                    write(kStdWarn,*) 'will zero out UA Planck coeffs from layer ',iZeroPlanck
                END IF
            END IF

            CALL logrspl(raUpperPress1013,raUpperTemp,iUpper,raUpperPress_Std,raInterpTempUA,iUpperStd_Num)
            iUpper = iUpperStd_Num
        !! thus we have mapped raUpperPress1013,raUpperTemp,iUpper  ----->>>>
        !!                     raUpperPress_Std,raInterpTempUA,iUpperStd_Num

            IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
                write(kStdWarn,*) 'UPPER ATM Profile for kCompressed (use default MixRatios in UA)'
                write(kStdWarn,*) ' iL     StdPress    StdTemp   StdCO2Amt  UserTemp  User-Std'
                write(kStdWarn,*) '          mb          K       kmol/cm2      K        K'
                write(kStdWarn,*) '------------------------------------------------------------'
                DO iL = 1,iUpperStd_Num
                    write(kStdWarn,*) iL,raUpperPress_Std(iL),raUpperTemp_Std(iL),raUpperCO2Amt_Std(iL), &
                    raInterpTempUA(iL),raInterpTempUA(iL)-raUpperTemp_Std(iL)
                    IF (iL+1 == iZeroPlanck) THEN
                        write(kStdWarn,*) '------------------------------------------------------------'
                    END IF
                END DO
            END IF
        END IF

        IF (kNLTEOutUAOpen == -1) THEN
        !!! open the output file etc
            iType = +1
            CALL OpenUAFile(iType,iUpper,caOutUAFile,rFreqStart,rFreqEnd, &
            iDumpAllUASpectra,iDumpAllUARads, &
            iFileIDLo,iFileIDHi,iTag, &
            iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
            IF (iDumpAllUARads > 0) THEN
                iNumberUA_NLTEOut = 1 + iUpper
            ELSEIF (iDumpAllUARads <= 0) THEN
                iNumberUA_NLTEOut = 1 + 1
            END IF
        END IF

    END IF

    IF (rSolzen >= 00.0 .AND. rSolzen < 40.0) THEN
        rSolzenX = 00.0
        rSolzenY = 40.0
    ELSEIF (rSolzen >= 40.0 .AND. rSolzen < 60.0) THEN
        rSolzenX = 40.0
        rSolzenY = 60.0
    ELSEIF (rSolzen >= 60.0 .AND. rSolzen < 80.0) THEN
        rSolzenX = 60.0
        rSolzenY = 80.0
    ELSEIF (rSolzen >= 80.0 .AND. rSolzen < 85.0) THEN
        rSolzenX = 80.0
        rSolzenY = 85.0
    ELSEIF (rSolzen >= 85.0 .AND. rSolzen <= 90.0) THEN
        rSolzenX = 85.0
        rSolzenY = 90.0
    END IF

    kFrStep = kaFrStep(iTag)
    iNLTEStart = 1            !! since we have the compressed ODs and planck
!! modifiers all the way thru the atm

    IF (iLTEIn > 0) THEN
    ! read in the NON LTE temperatures amd Vibrational Partition Fcns
        iBand = 1
        iGasID = iaGases(iGas)
        CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
        iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
        daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,iLineMixBand,iDoVoigtChi)

        CALL read_nonlte_temperature(iGasID,iISO,iLTEin,iBand,caaNLTETemp, &
        pProf,raPressLevels,raLayerHeight,raThickness,iProfileLayers, &
        raTPress,raTPartPress,raTTemp,raTAmt,daJL,daJU, &
        iaJ_UorL,raLTETemp,raNLTETemp,raVibQFT,iAllLayersLTE,dVibCenter)

        iWhichChunk = NewDataChunk(iLTEIn,iaNLTEChunks,iaaNLTEChunks,rFileStartFr)

        IF ((iWhichChunk > 0) .AND. (iaGases(iGas) == 2)) THEN
        ! uncompress lower atm ODs
            i_NLTEFile_TYPE = 100 + nint(rSolzenX)
            CALL compressedNLTE(iaGases(iGas),rFileStartFr,iTag,iActualTag, &
            kProfLayer,iL_low,iL_high, &
            raTAmt,raRAmtx,raTTemp,raRTempx, &
            iErr,iDoDQ,pProf,iProfileLayers, &
            daaDQ,daaDT,daaGasAbCoeff,iSplineType,i_NLTEFile_TYPE, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)


        ! uncompress lower atm planck
            i_NLTEFile_TYPE = 200 + nint(rSolzenX)
            CALL compressedNLTE(iaGases(iGas),rFileStartFr,iTag,iActualTag, &
            kProfLayer,iL_low,iL_high, &
            raTAmt,raRAmtx,raTTemp,raRTempx, &
            iErr,iDoDQ,pProf,iProfileLayers, &
            daaDQp,daaDTp,daaPlanckCoeff,iSplineType,i_NLTEFile_TYPE, &
            iaP1,iaP2,raP1,raP2, &
            iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
            iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
            iaQ11,iaQ12,raQ11,raQ12, &
            iaQ21,iaQ22,raQ21,raQ22)


            IF (kJacobian >= 0) THEN
            ! T dependance in both abs coeffs and planck
            ! q dependance only in abs coeffs
                DO iL = 1,kProfLayer
                    DO iFr = 1,kMaxPts
                        daaDT(iFr,iL) = daaDT(iFr,iL) + daaDTp(iFr,iL)
                    END DO
                END DO
            END IF

        !!!now do the stuff ABOVE the kCARTA TOA (above 0.005 mb)
            IF (iDoUpperAtmNLTE > 0) THEN
                write(kStdWarn,*) ' '
                write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>----------------------->>>>>>>>>>>>>>>'

                write(kStdWarn,*) 'Doing stratosphere NLTE abs coeff ....'
                write(kStdWarn,*) 'upper atm NLTE abs coeffs ...'

                DO iFr = 1,iUpper
                !! these are pressures in mb
                    pProfNLTE_upatm(iFr)      = raUpperPress(iFr)*kAtm2mb
                    pProfNLTE_upatm(iFr)      = raUpperPress_Std(iFr)*kAtm2mb
                    raUpperPartPress_Std(iFr) = raUpperPress_Std(iFr) * 385 * 1e-6
                END DO

                CALL read_std_optdepths_upper_UA( &
                iaGases(iGas),iSplineType,raUpperPressLevels, &
                pProfNLTE_upatm,iUpper, &
                rFileStartFr,iTag,iActualTag, &
                iaGases(iGas),raFreq, &
                raUpperCO2Amt_Std,raInterpTempUA,raUpperPress_Std, &
                raUpperPartPress_Std, &
                iUpper,daaWeakOptDepth)

            ! uncompress upper atm ODs ... no jacs possible ..
            ! so use daaDQp,daaDTp
                i_NLTEFile_TYPE = 300 + nint(rSolzenX)
                CALL compressedNLTE(iaGases(iGas),rFileStartFr,iTag,iActualTag, &
                kProfLayer,1,iUpperStd_Num, &
                raUpperCO2Amt_Std,raUpperCO2Amt_Std,raInterpTempUA,raUpperTemp_Std, &
                iErr,iDoDQ,raUpperPress_Std,iUpperStd_Num, &
                daaDQp,daaDTp,daaUpperNLTEGasAbCoeff,iSplineType,i_NLTEFile_TYPE, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

            ! uncompress upper atm planck ... no jacs possible
                i_NLTEFile_TYPE = 400 + nint(rSolzenX)
                CALL compressedNLTE(iaGases(iGas),rFileStartFr,iTag,iActualTag, &
                kProfLayer,1,iUpperStd_Num, &
                raUpperCO2Amt_Std,raUpperCO2Amt_Std,raInterpTempUA,raUpperTemp_Std, &
                iErr,iDoDQ,raUpperPress_Std,iUpperStd_Num, &
                daaDQp,daaDTp,daaUpperPlanckCoeff,iSplineType,i_NLTEFile_TYPE, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

                rMult = 1.0
                DO iL = 1,iZeroPlanck-1
                !              print *,'mwah',iL,daaWeakOptDepth(1,iL),daaUpperNLTEGasAbCoeff(1,iL),daaUpperPlanckCoeff(1,iL)
                    DO iFr = 1,kMaxPts
                        daaUpperPlanckCoeff(iFr,iL) = daaUpperPlanckCoeff(iFr,iL) * rMult
                    END DO
                    rMult = rMult * rMult0
                END DO

                DO iL = iZeroPlanck,kProfLayer
                    DO iFr = 1,kMaxPts
                        daaUpperPlanckCoeff(iFr,iL) = 0.0d0
                    END DO
                END DO

                DO iL = 1,kProfLayer
                    DO iFr = 1,kMaxPts
                        daaUpperSumNLTEGasAbCoeff(iFr,iL) = daaUpperNLTEGasAbCoeff(iFr,iL)
                    END DO
                END DO

                IF (iDumpAllUASpectra > 0) THEN
                ! dump out the UA spectra, iUpper paths worth of them!
                    CALL wrtout_head_uafile(caOutUAFile, &
                    raFreq(1),raFreq(kMaxPts),raFreq,iTag,1,iUpper)
                    caOutName = 'DumDum'
                    iIOUN = kNLTEOutUA
                    DO iL = 1,iUpper
                        rjunk = sngl(daaUpperSumNLTEGasAbCoeff(1,iL))
                        write(kStdWarn,*) 'dump out UA ODs : gid, layer, amt, OD = ', &
                        iaGases(iGas),iL,raUpperCO2Amt_Std(iL),rJunk
                        DO iFr = 1,kMaxPts
                            raX(iFr) = sngl(daaUpperSumNLTEGasAbCoeff(iFr,iL))
                        END DO
                        CALL wrtout(iIOUN,caOutName,raFreq,raX)
                    END DO
                END IF

            !             do iL = 1,kProfLayer
            !               print *,'place A',iaGases(iGas),iL,daaGasAbCoeff(1,iL),daaPlanckCoeff(1,iL),
            !     $                    daaUpperNLTEGasAbCoeff(1,iL),daaUpperPlanckCoeff(1,iL)
            !             end do

            END IF
        END IF
    END IF

    RETURN
    end SUBROUTINE NLTE_Fast_Compressed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! SLOW LBL calcs for NLTE ODs and Planck Coeffs
    SUBROUTINE NLTE_SLOW_LBL( &
    iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
    iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
    iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
    iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
    iUseWeakBackGnd,raFreq,pProf, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
    pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
    raRAmt,raRTemp,raRPress,raRPartPress, &
    raVertTemp,iVertTempSet, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    raUpperPress,raUpperPartPress,raUpperTemp, &
    raUpperGasAmt,raUpperNLTETemp, &
    iUpper,iDoUpperAtmNLTE, &
    dLineStrenMin,dDeltaFreqNLTE, &
    caaUpperMixRatio,iNumberUA_NLTEOut, &
    rFreqStart,rFreqEnd,rFileStartFr, &
    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
    caOutUAFile,caOutUABloatFile, &
    iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
    daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    daFreqBloat, &
    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
    daaPlanckCoeffBloat, &
    daaUpperPlanckCoeff, &
    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
! rCO2mult  = tells how much input .rtp file wants the CO2 to be multiplied by
!             remeber this is as compared to kCO2ppmv
    REAL :: rCO2mult
! iaNLTEBands     tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
! iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
    INTEGER :: iaNLTEBands(kGasStore)
    INTEGER :: iGas,iaGases(kMaxGas),iSplineType,iUseWeakBackGnd
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    INTEGER :: iNumNLTEGases,iSetBloat
    INTEGER :: iTag,iActualTag,iProfileLayers,iL_low,iL_high
    INTEGER :: iVertTempSet
    INTEGER :: iUpper,iNumberUA_NLTEOut
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iNLTE_SlowORFast
    REAL :: rFileStartFr
! caaaNLTEBands     tells the name of the files containing the line parameters
! caaNLTETemp       tells the name of the files containing the nonLTE temps
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaStrongLines(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    CHARACTER(80) :: caOutUAFile,caOutUABloatFile
! pProf is the avg layer pressure
    REAL :: pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
    REAL :: raNLTEstrength(kGasStore),raFreq(kMaxPts)
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    REAL :: raVertTemp(kProfLayer)
! stuff for the upper levels
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    REAL :: rFreqStart,rFreqEnd
! strengths etc
    DOUBLE PRECISION :: dLineStrenMin,dDeltaFreqNLTE

! output params
! iLTEIn            is just a variable that tells us if nonLTE to be used
    INTEGER :: iLTEIn,iWhichChunk,iFunnyCousin
! tells the nonscattering code, at which layer to start NONLTE rad transfer
    INTEGER :: iNLTEStart
! daaGasAbCoeff has the uncompressed gas absorption coeff
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer),raaCO2_LTE(kMaxPts,kProfLayer)
! daaNLTEGasAbCoeff has the nonLTE gas absorption coeff
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
! this has the bloated stuff
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
! this is ABOVE the 0.005 mb kcarta TOA
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local vars, mainly functions called
!     iDoVoigtChi = -1  !! do not multiply by chi function in voigt_chi; all
!                      !! necesary things will be done by calling
!                      !! Co2_4um_fudge_nlte_fast
!     iDoVoigtChi = +1  !! do multiply by chi function in voigt_chi
    INTEGER :: iDoVoigtChi
    INTEGER :: iL,iFr,iDoNLTE

    INTEGER :: iJunkNum,iaJunk(kGasStore)

! see if current gas ID needs nonLTE spectroscopy
    iLTEIn      = -1
    iWhichChunk = -1
    iDoNLTE     = -1
    iJunkNum    = -1

!! set up nlte strengths etc
    CALL SetUpNLTEStrengths(dLineStrenMin,dDeltaFreqNLTE,iDoVoigtChi)

    iLTEIn = OutsideSpectra(iaGases(iGas),iNumNLTEGases,iaNLTEGasID,iJunkNum,iaJunk,raFreq(1),605.0,2830.0,20)
    IF (iLTEIn > 0) THEN
        iWhichChunk = NewDataChunk(iLTEIn,iaNLTEChunks,iaaNLTEChunks,rFileStartFr)

        IF ((iWhichChunk > 0) .AND. (iaGases(iGas) == 2) &
         .AND. (abs(raNLTEstrength(iLTEIn)-1.0) <= 0.01))  THEN
            raNLTEstrength(iLTEIn) = rCO2mult
        !! change the coeffs to coeff*rCO2mult, because all the coeffs
        !! that are computed, have this factor
            iDoNLTE = +1
            DO iL = 1,kProfLayer
                DO iFr = 1,kMaxPts
                    daaGasAbCoeff(iFr,iL) = daaGasAbCoeff(iFr,iL) * rCO2mult
                END DO
            END DO
        END IF

        IF ((iWhichChunk > 0) .AND. (iaGases(iGas) == 2) &
         .AND. (raNLTEstrength(iLTEIn) < 0))  THEN
            IF (iSetBloat > 0) THEN
                write(kStdErr,*) 'Cannot bloat up H92 database calculations'
                write(kStdErr,*) 'used with GENLN2 (Cousin)'
                CALL DoStop
            END IF
        !!! do a computation where you change the linemix spectra to
        !!! cousin spectra, using old kCARTA database
            write(kStdWarn,*) 'Replacing kCARTA database optical depths '
            write(kStdWarn,*) 'with GENLN2 cousin optical depths'
            CALL CousinContribution(iaGases(iGas), &
            rFileStartFr,iTag,iActualTag,iProfileLayers,iL_low,iL_high, &
            raTAmt,raRAmt,raTTemp,raRTemp,pProf, &
            raNLTEstrength(iLTEIn),iaNLTEStart(iLTEin), &
            iUpper,daaGasAbCoeff,iSplineType)
            iFunnyCousin = +1

        ELSEIF (iWhichChunk > 0) THEN
        !!!first compute LTE background optical depths
            IF (iUseWeakBackGnd > 0) THEN
                write(kStdWarn,*) 'Replacing kCARTA database optical depths with backgnd'
                write(kStdWarn,*) 'LTE abs coeffs in necessary upper layers ....'
            ! use pProfNLTE or pProf??????
                CALL lte_spectra(iTag,iActualTag,iLTEin,iaNLTEStart(iLTEin),raFreq, &
                iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
                pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
                raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
                iProfileLayers,iSplineType, &
                raVertTemp,iVertTempSet,rFileStartFr, &
                raTAmt,raTTemp,raTPress,raTPartPress,iaGases(iGas), &
                daaNLTEGasAbCoeff,daaPlanckCoeff,daaSumNLTEGasAbCoeff, &
                iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
                dLineStrenMin,dDeltaFreqNLTE, &
                iDoVoigtChi,raNLTEstrength(iLTEIn))
                CALL AddNLTECoeffs(daaGasAbCoeff,daaNLTEGasAbCoeff, &
                iaNLTEStart(iLTEin),-1)
            ELSE
                CALL AddNLTECoeffs(daaGasAbCoeff,daaNLTEGasAbCoeff, &
                iaNLTEStart(iLTEin),0)
            END IF

            IF (iSetBloat > 0) THEN
            !!!! initialize the bloated matrices to LTE stuff
            ! use pProfNLTE or pProf??????
                CALL BloatCoeffsDriver(iTag,iActualTag, &
                rFileStartFr,raFreq,daaGasAbCoeff, &
                raaRestOfLTEGases,raaCO2_LTE, &
                daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
                daFreqBloat,iaNLTEStart(iLTEin),raNLTEstrength(iLTEIn), &
                daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
                daaPlanckCoeffBloat, &
                iGas,iaGases(iGas),iL_low,iL_high,iUseWeakBackGnd, &
                raRAmt,raRTemp,raRPress,raRPartPress, &
                pProfNLTE,iProfileLayers, &
                raTAmt,raTTemp,raTPress,raTPartPress,iSplineType)
            END IF
        !!!!compute the line shapes for the LineMix bands
        !!!!some of which are in LTE, some in NLTE
        !!!!also find contribution for numerator in beta (planck coeff)
            write(kStdWarn,*) 'Adding in NLTE upper layers abs coeff ....'
        ! use pProfNLTE or pProf??????
            CALL nonlte_spectra(iTag,iActualTag,iLTEin, &
            iaNLTEStart(iLTEin),iaNLTEStart2350(iLTEin),raFreq, &
            iaNLTEBands,caaaNLTEBands,caaNLTETemp, &
            pProfNLTE,raPressLevels,raLayerHeight, &
            raThickness,iProfileLayers, &
            raTAmt,raTTemp,raTPress,raTPartPress, &
            daaNLTEGasAbCoeff,daaPlanckCoeff,daaSumNLTEGasAbCoeff, &
            iDoUpperAtmNLTE,iAllLayersLTE,dDeltaFreqNLTE, &
            raNLTEstrength(iLTEIn),iDoVoigtChi, &
            iSetBloat,daaSumNLTEGasAbCoeffBloat,daaNLTEGasAbCoeffBloat, &
            daaPlanckCoeffBloat,daFreqBloat)

        !!!update at which layer NONLTE starts
            iNLTEStart = min(iNLTEStart,iaNLTEStart(iLTEin))
            iNLTEStart = min(iNLTEStart,iaNLTEStart2350(iLTEin))

            CALL AddNLTECoeffs(daaGasAbCoeff,daaNLTEGasAbCoeff, &
            iNLTEStart,-1)

        !!!now do the stuff ABOVE the kCARTA TOA (above 0.005 mb)
            IF (iDoUpperAtmNLTE > 0) THEN
                write(kStdWarn,*) ' '
                write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>----------------------->>>>>>>>>>>>>>>'
                write(kStdWarn,*) 'Doing stratosphere NLTE abs coeff ....'
                write(kStdWarn,*) 'upper atm NLTE abs coeffs ...'
                IF (iUseWeakBackGnd > 0) THEN
                    write(kStdWarn,*) 'backgrnd upper atm LTE abs coeffs ...'
                    CALL lte_spectra_upper( &
                    iTag,iActualTag, &
                    rFileStartFr,iaNLTEStart(iLTEIn),rCO2mult,raFreq, &
                    iaGases(iGas),iLTEIn,iaNLTEBands,caaStrongLines, &
                    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
                    pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
                    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
                    raUpperPress,raUpperPartPress,raUpperTemp, &
                    raUpperGasAmt,raUpperNLTETemp, &
                    daaUpperPlanckCoeff, &
                    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
                    daaUpperPlanckCoeffBloat, &
                    daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                    iUpper,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
                    dLineStrenMin,dDeltaFreqNLTE, &
                    iDoVoigtChi,raNLTEstrength(iLTEIn),iSetBloat,daFreqBloat, &
                    iGas,iSplinetype,iNumberUA_NLTEOut, &
                    rFreqStart,rFreqEnd, &
                    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
                    iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks, &
                    caOutUAFile,caOutUABloatFile)
                !! have already put abs coeffs into daaUpperNLTEGasAbCoeff,
                !! and ignore other gases, so just forget about AddNLTECoeffs
                END IF
                CALL UpHigh_Nonlte_Spectra(iTag,iActualTag,iaNLTEStart(iLTEIn), &
                raFreq,iaGases(iGas),iLTEIn,iaNLTEBands, &
                caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,rCO2mult, &
                pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
                pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
                raUpperPress,raUpperPartPress,raUpperTemp, &
                raUpperGasAmt,raUpperNLTETemp, &
                daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff, &
                daaUpperSumNLTEGasAbCoeff, &
                iUpper,iDoUpperAtmNLTE,iAllLayersLTE, &
                dDeltaFreqNLTE,iDoVoigtChi,raNLTEstrength(iLTEIn), &
                iSetBloat,daFreqBloat, &
                daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                daaUpperNLTEGasAbCoeffBloat, &
                iGas,iSplinetype,iNumberUA_NLTEOut, &
                rFreqStart,rFreqEnd, &
                iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
                iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks, &
                caOutUAFile,caOutUABloatFile)
            END IF   !IF (iDoUpperAtmNLTE > 0) THEN
        END IF     !IF ((iWhichChunk > 0) .AND. 
    !      (raNLTEstrength(iLTEIn) .LT. 0)
    END IF       !IF (iLTEIn > 0) THEN

    IF ((iWhichChunk > 0) .AND. (iaGases(iGas) == 2) &
     .AND. (iDONLTE > 0))  THEN
        raNLTEstrength(iLTEIn) = rCO2mult
    !! change the coeffs to coeff*rCO2mult, because all the coeffs
    !! that are computed, have this factor
        DO iL = 1,kProfLayer
            DO iFr = 1,kMaxPts
                daaGasAbCoeff(iFr,iL) = daaGasAbCoeff(iFr,iL)/rCO2mult
            END DO
        END DO
    END IF

    IF (iLTEIn > 0) THEN
        IF ((iWhichChunk > 0) .AND. (iaGases(iGas) == 2) &
         .AND. (raNLTEstrength(iLTEIn) < 0))  THEN
            iLTEIn = -1
            write (kStdWarn,*) 'using the linemix CO2 database from GENLN2'
        ELSEIF ((iWhichChunk < 0) .OR. (raNLTEstrength(iLTEIn) < 0)) THEN
            iLTEIn = -1
            write (kStdWarn,*) 'this is NLTE gas, but is NOT a NLTE chunk!'
        END IF
    END IF

    RETURN
    END SUBROUTINE NLTE_SLOW_LBL

!************************************************************************
! this sets up some params
    SUBROUTINE SetUpNLTEStrengths(dLineStrenMin,dDeltaFreqNLTE,iDoVoigtChi)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! dLineStrenMin is which min line strength to use
    DOUBLE PRECISION :: dLineStrenMin    !!to prune the database
    DOUBLE PRECISION :: dDeltaFreqNLTE       !!dF to use (default = 0.0025 cm-1)
!     iDoVoigtChi = +1  !! do multiply by chi function in voigt_chi
    INTEGER :: iDoVoigtChi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dLineStrenMin = -1.0d0     !!! use all backgnd lines
    dDeltaFreqNLTE = -1.0d0    !!! use default dkFreqTag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dLineStrenMin = 1.0d-23    !!! use backgnd lines bigger than this
    dDeltaFreqNLTE = 1.0d-3    !!! use GENLN2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dLineStrenMin = -1.0d0     !!! use all backgnd lines
    dDeltaFreqNLTE = -1.0d0    !!! use default dkFreqTag

    iDoVoigtChi = +1           !!! multiply by chi fcn in voigt_chi

    RETURN
    end SUBROUTINE SetUpNLTEStrengths

!************************************************************************
! this subroutine checks to see if we need to zero the daaPlanckCoeff for
! the particular chunk
    SUBROUTINE ZeroPlanckCoeff(iNumLayers, &
    rStartBlock,iTag,iDoUpperAtmNLTE, &
    daaPlanckCoeff,daaSumNLTEGasAbCoeff,daaNLTEGasAbCoeff, &
    daaUpperPlanckCoeff, &
    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
    iChunk_DoNLTE,iNumGases,iaGases, &
    iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,rFileStartFr, &
    iSetBloat,daaSumNLTEGasAbCoeffBloat, &
    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,daFreqBloat, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat, &
    rFrLow,rFrHigh,iTotalStuff,iFileIDLo,iFileIDHi, &
    raaRestOfLTEGases,raaCO2_LTE,caOutBloatFile,caPlanckBloatFile)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iNumLayers                 !number of layers in atm #1
    REAL ::    rStartBlock                !start freq of current 10000 pt chunk
    INTEGER :: iNumNLTEGases              !number of nonlte gases to consider
    INTEGER :: iNLTE_SlowORFast           !use slow (+1) or fast (-1/-2) model
    INTEGER :: iaNLTEGasID(kGasStore)     !gasIDs that are not in LTE
    INTEGER :: iaNLTEChunks(kGasStore)    !for above gases, how many data sets
! ot in LTE
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    REAL ::    rFileSTartFr               !current start chunk
    INTEGER :: iNumGases                  !number of gases
    INTEGER :: iaGases(kMaxGas)           !which gas IDs
    INTEGER :: iDoUpperAtmNLTE            !should we do upper stratosphere
    INTEGER :: iSetBloat                  !do we bloat to 0.0005 cm-1
    INTEGER :: iTag                       !what is the spacing we need?
    CHARACTER(80) :: caOutBloatFile        !file to dump high res info to
    CHARACTER(80) :: caPlanckBloatFile     !file to dump high res info to
    REAL :: rFrLow,rFrHigh                !start and stop freqs
    INTEGER :: iTotalStuff                !total number of outputs
    INTEGER :: iFileIDLo,iFileIDHi        !which 0.0025 kchunks to expect
! output parameters
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbcoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)

    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

    REAL :: raaRestOfLTEGases(kMaxPts,kProfLayer)
    REAL ::  raaCO2_LTE(kMaxPts,kProfLayer)

    INTEGER :: iChunk_DoNLTE,iUsualLTEGas
          
! local variables
    INTEGER :: iLTEIn,iWhichChunk,iFr,iL,iGas,iJump,iType,iTypeUA
    INTEGER :: iIOUN,iFileErr
    DOUBLE PRECISION :: dTemp1,dfFine
    INTEGER :: iJunkNum,iaJunk(kGasStore)

    iChunk_DoNLTE = -1 !!assume everything IS in LTE
    iJunkNum    = -1

    write(kStdWarn,*) 'Checking NLTE Gases using OutsideSpectra ...'
    DO iGas = 1,iNumGases
        iWhichChunk = -1
        iLTEIn = OutsideSpectra(iaGases(iGas),iNumNLTEGases,iaNLTEGasID,iJunkNum,iaJunk,2205.0,605.0,2830.0,20)
        IF ((iLTEIn > 0) .AND. (kSolarAngle >= 0 .AND. kSolarAngle <= 90)) THEN
            iWhichChunk = NewDataChunk(iLTEIn,iaNLTEChunks,iaaNLTEChunks,rFileStartFr)
            IF ((iWhichChunk > 0) .AND. (iNLTE_SlowORFast == +1)) THEN
                iChunk_DoNLTE = +1           !!!need to do NLTE slowly
                GOTO 10
            ELSEIF ((iWhichChunk > 0) .AND. (iNLTE_SlowORFast == -1)) THEN
                iChunk_DoNLTE = +2           !!!need to do NLTE with SARTA model
                GOTO 10
            ELSEIF ((iWhichChunk > 0) .AND. (iNLTE_SlowORFast == -2)) THEN
                iChunk_DoNLTE = +3          !!!need to do NLTE with compressed database
                GOTO 10
            END IF
        END IF
    END DO

    10 CONTINUE

    IF ((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) THEN
        DO iL = 1,kProfLayer
            DO iFr = 1,kMaxPts
            !!! do the kCARTA 100 AIRS layers
                daaSumNLTEGasAbCoeff(iFr,iL) = 0.0d0
                daaNLTEGasAbCoeff(iFr,iL)    = 0.0d0
                daaPlanckCoeff(iFr,iL)       = 0.0d0
            END DO
        END DO
    END IF

    IF (((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) .AND. &
    (iDoUpperAtmNLTE == 1)) THEN
        DO iL = 1,kProfLayer
            DO iFr = 1,kMaxPts
            !!! do the upper atm layers
                daaUpperSumNLTEGasAbCoeff(iFr,iL) = 0.0d0
                daaUpperNLTEGasAbCoeff(iFr,iL)    = 0.0d0
                daaUpperPlanckCoeff(iFr,iL)       = 0.0d0
            END DO
        END DO
    END IF

! set up the high res wavenumber array
    IF ((iChunk_DoNLTE == 1) .AND. (iSetBloat == 1)) THEN
        iJump = iFloor(kBoxCarUse*1.0/2)
        dfFine = kaFineFrStep(iTag)
        dTemp1 = rStartBlock*1.0d0
        DO iFr = 1,kMaxPts*kBoxCarUse-iJump
        ! do the bulk of the points
            daFreqBloat(iFr+iJump) = dTemp1 + (iFr-1)*dfFine
        END DO
        DO iFr = 1,iJump
        ! do the end points
            daFreqBloat(iFr) = dTemp1 - (3-iFr)*dfFine
        END DO
        DO iL = 1,kProfLayer
            DO iFr = 1,kBloatPts
            !!! do the kCARTA 100 AIRS layers
                daaSumNLTEGasAbCoeffBloat(iFr,iL) = 0.0d0
                daaNLTEGasAbCoeffBloat(iFr,iL)    = 0.0d0
                daaPlanckCoeffBloat(iFr,iL)       = 0.0d0
                daaUpperPlanckCoeffBloat(iFr,iL)       = 0.0d0
                daaUpperSumNLTEGasAbCoeffBloat(iFr,iL) = 0.0d0
                daaUpperNLTEGasAbCoeffBloat(iFr,iL)    = 0.0d0
            END DO
        END DO

        DO iL = 1,kProfLayer
            DO iFr = 1,kMaxPts
            !!! do the kCARTA 100 AIRS layers
                raaRestOfLTEGases(iFr,iL) = 0.0
                raaCO2_LTE(iFr,iL)        = 0.0
            END DO
        END DO
    END IF

!!! BLOATING can only be done by iNLTE_SlowORFast == +1 or equivalently
!!!                              iChunk_DoNLTE    == +1
    IF ((iChunk_DoNLTE == 1) .AND. (kBloatOutOpen < 0)) THEN
        IF (iSetBloat == 1) THEN
            iType = +1
            CALL OpenOutputBloatFile( &
            iType,iNumLayers,caOutBloatFile,rFrLow,rFrHigh, &
            iFileIDLo,iFileIDHi,iTag,iTotalStuff, &
            iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        END IF
    END IF

!!! kPlanckOut is set by parameter kFlux in nm_params section
!!! if set to 0, then subroutine in s_writefile turns on planck writing
    IF ((iChunk_DoNLTE == 1) .AND. (kBloatPlanckOpen < 0) &
     .AND. (kPlanckOut == 0)) THEN
        IF (iSetBloat == 1) THEN
            iType = -1
            CALL OpenOutputBloatFile( &
            iType,iNumLayers,caPlanckBloatFile,rFrLow,rFrHigh, &
            iFileIDLo,iFileIDHi,iTag,iTotalStuff, &
            iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        END IF
    END IF

    RETURN
    end SUBROUTINE ZeroPlanckCoeff

!************************************************************************
! this simply adds on two double matrices, from layers iSTart -> kProflayer
! store results in daaSave
! if iAdd = +1, then you add the abs coeffs together in relevant layers
! if iAdd = -1, then you replace daa with daaAdd in relevant layers
! if iAdd =  0, then you replace daa with 0.0 in relevant layers

! note that here we are adding CUMULATIVE things, so we can check to ensure
! kAbs > 0!!
    SUBROUTINE AddNLTECoeffs(daaUpdate,daaAdd,iStart,iAdd)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
    DOUBLE PRECISION :: daaAdd(kMaxPts,kProfLayer)
    INTEGER :: iStart,iAdd
! input/output
    DOUBLE PRECISION :: daaUpdate(kMaxPts,kProfLayer)

! local vars
    REAL :: rX
    INTEGER :: iFr,iL

    rX = 1.0   !!since now ALL the subroutines use raNLTEstrength(iLTEIn)

    IF (iAdd > 0) THEN
    ! ou are adding daaUpdate + daaAdd in the relevant layers
    ! his could be eg you already have the weak background from run7co2,
    ! nd now you are adding on the nonLTE contribution
        DO iL = iStart,kProfLayer
            DO iFr = 1,kMaxPts
                daaUpdate(iFr,iL) = max(daaUpdate(iFr,iL)+daaAdd(iFr,iL)*rX,0.0d0)
            END DO
        END DO
    END IF

    IF (iAdd < 0) THEN
    ! ou are replacing daa with daaAdd in the relevant layers
    ! his could be eg you have uncompressed the kCARTA database, but now
    ! ou need to know the background abs coeffs in LTE
        DO iL = iStart,kProfLayer
            DO iFr = 1,kMaxPts
                daaUpdate(iFr,iL) = max(daaAdd(iFr,iL)*rX,0.0d0)
            END DO
        END DO
    END IF

    IF (iAdd == 0) THEN
    ! ou are replacing daaUpdate with 0.0 in the relevant layers
        DO iL = iStart,kProfLayer
            DO iFr = 1,kMaxPts
                daaUpdate(iFr,iL) = 0.0d0
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE AddNLTECoeffs

!************************************************************************
! this computes the cumulative multiplication modification to Planck, which is
! just 1.000000000000000
    SUBROUTINE SetPlanckCoeff_Cousin(iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iNLTEStart
! output params
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iL,iFr

    DO iL = 1,kProfLayer
        DO iFr = 1,kMaxPts
            raaPlanckCoeff(iFr,iL) = 1.00000
        END DO
    END DO

    RETURN
    end SUBROUTINE SetPlanckCoeff_Cousin

!************************************************************************
! this subroutine changes double --> real for the KCARTA atmosphere
! this computes the cumulative multiplication modification to Planck
! this is for the usual 100 kCARTA database layers
    SUBROUTINE SetPlanckCoeff(iChunk_DoNLTE,iNLTEStart,iAtm,iaaRadLayer, &
    daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
    raaSumAbcoeff,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
! cumulative nonlte gas optical depths
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer) !!planck mod, so far
    REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)             !!mixed path abscoeff
    INTEGER :: iAtm                                        !!which atmosphere
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)             !!set of mixed paths
    INTEGER :: iNLTEStart                                  !!where nonLTE starts
    INTEGER :: iChunk_DoNLTE                               !!if 3, then came from compr. files
! output params
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    DOUBLE PRECISION :: dTemp,dTemp1
    INTEGER :: iFr,iL,iM,iA

! initialize to 1.0!!!!
    DO iL = 1,kProfLayer
        DO iFr = 1,kMaxPts
            raaPlanckCoeff(iFr,iL) = 1.0
        END DO
    END DO

    iA = 1    !!!asuume current atmospher uses mixed path layers 1-100
    iM = iaaRadLayer(iAtm,1)
          
! ind which set of Mixed Paths current atm uses eg 1-100, 101-200 etc
    DO iL = 1,kProfLayer
        IF (iM <= kProfLayer*iL) THEN
            iA = iL
            GOTO 10
        END IF
    END DO
    10 CONTINUE
          
! recall raaSumAbCoeff(iFr,iM) is the CUMULATIVE abscoeff=sum(LTE) + sum(NLTE)
! where sum(LTE) = all gases that are in LTE, plus CO2 in LTE plus CO2 in NLTE
!       sum(NLTE) = CO2 in NLTE
! thus we need to break down Eqn 20 of the KOPRA paper appropriately
!                            Eqn 16 of the JGR D. Edwards paper appropriately
! by finding numerator of (Eqn 20 kopra or Eqn 16, JGR); denom = raaSumAbCoeff
    IF (iChunk_DoNLTE == 1) THEN
    !! actually computed this using LBL
        DO iL = iNLTEStart,kProfLayer
            iM = iL + (iA-1)*kProfLayer
            DO iFr = 1,kMaxPts
            !!!raaSumAbCoeff = raaSumLTE + raaSumNLTE
                dTemp = raaSumAbCoeff(iFr,iM)*1.0d0-daaSumNLTEGasAbCoeff(iFr,iL)
            ! dTemp is thus the LTE component of numerator (beta = 1.0);
            ! add on the NLTE component (where beta is some other number)
                dTemp = max(dTemp + daaPlanckCoeff(iFr,iL),0.0d0)
            ! normalise by raaSumAbCoeff
                dTemp  = dTemp + dDeltaNLTE
                dTemp1 = raaSumAbCoeff(iFr,iM)*1.0d0 + dDeltaNLTE
                raaPlanckCoeff(iFr,iL) = sngl(dTemp/dTemp1)
            END DO
        END DO
    ELSEIF (iChunk_DoNLTE == 3) THEN
    !! simply used compressed tables
        DO iL = iNLTEStart,kProfLayer
            iM = iL + (iA-1)*kProfLayer
            DO iFr = 1,kMaxPts
                raaPlanckCoeff(iFr,iL) = max(sngl(daaPlanckCoeff(iFr,iL)),0.0)
            END DO
        END DO
    END IF
       
    RETURN
    end SUBROUTINE SetPlanckCoeff

!************************************************************************
! this subroutine changes double --> real for the UPPER stratosphere
! this computes the cumulative multiplication modification to Planck
    SUBROUTINE SetUpperPlanckCoeff( &
    iChunk_DoNLTE,iUpper,daaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff, &
    raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat,iSetBloat)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iUpper,iSetBloat,iChunk_DoNLTE
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raFreq(kMaxPts)
! output parameters
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
! input/output parameters
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
    DOUBLE PRECISION :: dTemp
    INTEGER :: iFr,iL,iMethod

    IF (iChunk_DoNLTE == +3) THEN
        iMethod = -3    !!! we get necessary stuff from kCompressed NLTE files
    ELSE
        iMethod = -1    !!! assume only the NLTE components of CO2 added in
        iMethod = +1    !!! assume NLTE components of CO2, plus (LTE CO2, others)
    END IF

! initialize to 1.0!!!!
    DO iL = 1,kProfLayer
        DO iFr = 1,kMaxPts
            raaUpperPlanckCoeff(iFr,iL) = 1.0
        END DO
    END DO
      
    IF (iMethod == -3) THEN
    !! we get stuff direct from the kCompressed Files
        DO iL = 1,iUpper
            DO iFr = 1,kMaxPts
                raaUpperSumNLTEGasAbCoeff(iFr,iL) = sngl(daaUpperNLTEGasAbCoeff(iFr,iL))
                raaUpperPlanckCoeff(iFr,iL)       = sngl(daaUpperPlanckCoeff(iFr,iL))
            END DO
        END DO
              
    ! recall raaSumAbCoeff(iFr,iM) is the CUMULATIVE abscoeff=sum(LTE) + sum(NLTE)
    ! where sum(LTE) = all gases that are in LTE, plus CO2 in LTE plus CO2 in NLTE
    !       sum(NLTE) = CO2 in NLTE
    ! thus we need to break down Eqn 20 of the KOPRA paper appropriately
    !                            Eqn 16 of the JGR D. Edwards paper appropriately
    ! by finding numerator of (Eqn 20 kopra or Eqn 16, JGR); denom = raaSumAbCoeff

    ELSEIF (iMethod == +1) THEN
    !! have NLTE(CO2) plus (LTE CO2 and other gases)
    ! regular 0.0025 cm-1
        DO iL = 1,iUpper
            DO iFr = 1,kMaxPts
            !!!raaSumAbCoeff = raaSumLTE + raaSumNLTE
                dTemp = daaUpperSumNLTEGasAbCoeff(iFr,iL) - &
                daaUpperNLTEGasAbCoeff(iFr,iL)
            ! dTemp is thus the LTE component of numerator (beta = 1.0);
            ! add on the NLTE component (where beta is some other number)
                dTemp = max(dTemp + daaUpperPlanckCoeff(iFr,iL),0.0d0) + dDeltaNLTE
            ! normalise by raaSumAbCoeff
                dTemp = dTemp/(daaUpperSumNLTEGasAbCoeff(iFr,iL) + dDeltaNLTE)
                raaUpperPlanckCoeff(iFr,iL)    = real(dTemp)
                dTemp = daaUpperSumNLTEGasAbCoeff(iFr,iL)
                raaUpperSumNLTEGasAbCoeff(iFr,iL) = real(dTemp)
            END DO
        END DO

    ! bloated 0.0005 cm-1
        IF (iSetBloat > 0) THEN
            DO iL = 1,iUpper
                DO iFr = 1,kBloatPts
                !!!raaSumAbCoeff = raaSumLTE + raaSumNLTE
                    dTemp = daaUpperSumNLTEGasAbCoeffBloat(iFr,iL) - &
                    daaUpperNLTEGasAbCoeffBloat(iFr,iL)
                ! dTemp is thus the LTE component of numerator (beta = 1.0);
                ! add on the NLTE component (where beta is some other number)
                    dTemp = max(dTemp + daaUpperPlanckCoeffBloat(iFr,iL),0.0d0) &
                    + dDeltaNLTE
                ! normalise by raaSumAbCoeff
                    dTemp = dTemp/(daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)+dDeltaNLTE)
                    daaUpperPlanckCoeffBloat(iFr,iL)    = real(dTemp)
                !              dTemp = daaUpperSumNLTEGasAbCoeff(iFr,iL)
                !              raaUpperNLTEGasAbCoeff(iFr,iL) = real(dTemp)
                END DO
            END DO
        END IF
                   
    ELSEIF (iMethod == -1) THEN
    ! don't worry about the LTE parts of CO2 and other gases
    ! regular 0.0025 cm-1
        DO iL = 1,iUpper
            DO iFr = 1,kMaxPts
                dTemp = max(daaUpperPlanckCoeff(iFr,iL),0.0d0)
                dTemp = dTemp + dDeltaNLTE
            ! normalise by raaSumAbCoeff
                dTemp = dTemp/(daaUpperSumNLTEGasAbCoeff(iFr,iL) + dDeltaNLTE)
                raaUpperPlanckCoeff(iFr,iL)   = sngl(dTemp)
                dTemp = daaUpperSumNLTEGasAbCoeff(iFr,iL) + dDeltaNLTE
                raaUpperSumNLTEGasAbCoeff(iFr,iL)= sngl(dTemp)
            END DO
        END DO

    ! bloated 0.0005 cm-1
        IF (iSetBloat > 0) THEN
            DO iL = 1,iUpper
                DO iFr = 1,kBloatPts
                    dTemp = max(daaUpperPlanckCoeffBloat(iFr,iL),0.0d0)
                    dTemp = dTemp + dDeltaNLTE
                ! normalise by raaSumAbCoeff
                    dTemp = dTemp/(daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)+dDeltaNLTE)
                    daaUpperPlanckCoeffBloat(iFr,iL)   = dTemp
                !             dTemp = daaUpperSumNLTEGasAbCoeff(iFr,iL)
                !             raaUpperNLTEGasAbCoeff(iFr,iL)= real(dTemp)
                END DO
            END DO
        END IF
    END IF

    RETURN
    end SUBROUTINE SetUpperPlanckCoeff

!************************************************************************
! This subroutine computes the LTE abs coeffs for the background lines,
!   for stuff in the KCARTA atmosphere
! The weakest lines are simple : just read in precomputed US Standard optical
!   depths, and adjust as necessary to current pressure layering and
!   gas amounts (see read_std_optdepths)
! The strong lines, for which run7co2.m does linemixing, are a little more
!   complicated, as the code uses the current profile to compute voigt optical
!   depths * birnbaum of these strong lines. These lines are in LTE so use the
!   current temperature profile

    SUBROUTINE lte_spectra(iTag,iActualTag,iLTEin,iStart,raFreq, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
    pProf,raPressLevels,raLayerHeight,raThickness, &
    raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
    iProfileLayers,iSplineType, &
    raVertTemp,iVertTempSet,rFileStartFr, &
    raTAmt,raTTemp,raTPress,raTPartPress,iGasID, &
    daaNLTEGasAbCoeff,daaPlanckCoeff,daaSumNLTEGasAbCoeff, &
    iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    dLineStrenMin,dDeltafreqNLTE,iDoVoigtChi,rNLTEstrength)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
! iTag      tells the current frequency step
! iLTEin    tells which set of data to use
! iStart    tells which layer to start NLTE calcs
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! caaStrongLines     line param files associated with strong lines, in LTE
! raTTemp etc     tell you the LTE gas profile
! raFreq          tells you the wavenumber array
! pProf           tells you the avg pressure in the AIRS layers
! dDeltaFreq are for default (0.0025 spacing, 5 point boxcar) or not
! dLineStrenMin is which min line strength to use
    DOUBLE PRECISION :: dDeltafreqNLTE,dLineStrenMin
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    REAL :: pProf(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    INTEGER :: iaNLTEBands(kGasStore),iLTEIn,iSTart,iTag,iActualTag,iGASID
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaStrongLines(kGasStore)
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raFreq(kMaxPts),rNLTEstrength
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! these are the individual reference profiles
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! this sets the vertical temp profile, to be checked for all gases
    REAL :: raVertTemp(kProfLayer),rFileStartFr
    INTEGER :: iVertTempSet
    INTEGER :: iL_low,iL_high,iSplineType
    INTEGER :: iProfileLayers   !number of layers in RTP file
! output parameter
! daaNLTEGasAbCoeffCoeff tells how to modify the abs coeff
! daaPlanckCoeff         tells how to modify the Planck function
! daaSumNLTEGasAbCoeff   tells running sum of NLTE GasAbCoeff
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbcoeff(kMaxPts,kProfLayer)

! local variables
    CHARACTER(80) :: caStrong
    INTEGER :: iFr,iL,iBand,iNum,iISO,iLineMixBand,iLowerOrUpper
    DOUBLE PRECISION :: daK(kMaxPts),dK
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN),daIso(kHITRAN)
    DOUBLE PRECISION :: dLTE
    DOUBLE PRECISION :: daChi(kMaxPts),daChiBloat(kBloatPts)
    DOUBLE PRECISION :: daaWeakOptDepth(kMaxPts,kProfLayer)
    INTEGER :: iDefault,iUnCompressType

! //////////////----> first compute lineshape of stronger lines, at LTE
    caStrong = caaStrongLines(iLTEIn)

!!! read in line params for the strong lines that are NOT in NLTE; so
!!! the default assumption is that they are in LTE
    CALL read_stronglineLTE_lineparameters( &
    iGasID,caaaNLTEBands,iaNLTEBands,iLTEin,dLineStrenMin,caStrong, &
    iNum,daIso,daElower,daLineCenter,daJL,daJU,daPshift, &
    daStren296,daW_For,daW_self,daW_temp)

    IF (iGasID == 2) THEN
    !!! use cousin for these bands in LTE,
    !!! don't get fancy with linemix/birnbaum right now, as we assume
    !!! (a) all weak LTE lines have been read in from the compressed
    !!!     database, computed using a Cousin lineshape
    !!! (b) these medium LTE lines will have a Cousin lineshape
    !!! (c) the strong NLTE lines will have a Cousin lineshape deep in the
    !!!     4 um band, or a linemix lineshape in the R branchhead
        iLineMixBand = +1
    ELSE
        iLineMixBand = -1
    END IF
! compute lineshapes, using lineparameters, at layers blah..kProflayer
    IF (iUseWeakBackGnd == 1) THEN
        write(kStdWarn,*) 'Adding in LTE contribution of strong lines ...'
        DO iL = iStart,kProfLayer   !!!loop over layers
            dLTE  = raTTemp(iL)*1.0d0
            IF (iNum > 0) THEN
                write (kStdWarn,*) 'LBL for Strong BackGnd (LTE) Lines : iL,T,Q = ',iL,sngl(dLTE),raTamt(iL)
                CALL compute_lte_spectra_fast(iTag,iActualTag, &
                daK,raFreq,iGasID,iNum,daIso, &
                daElower,daLineCenter,daJL,daJU,daPshift, &
                daStren296,daW_For,daW_self,daW_temp,dLTE,iLineMixBand, &
                iL,raTAmt,raTTemp,raTPress,raTPartPress, &
                dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength)
            ELSE
                write(kStdWarn,*) '   no Strong Bands are in LTE!!! -----> ',iL
                DO iFr = 1,kMaxPts
                    daK(iFr) = 0.0d0
                END DO
            END IF
            DO iFr = 1,kMaxPts
                dK = daK(iFr) * dble(rNLTEStrength)
                daaPlanckCoeff(iFr,iL)       = dK
                daaNLTEGasAbCoeff(iFr,iL)    = dK
                daaSumNLTEGasAbCoeff(iFr,iL) = dK
            END DO    !loop over layers
        END DO
    END IF        !IF (iUseWeakBackGnd == 1) THEN

    write(kStdWarn,*) ' '

! //////////////----> then read in optical depths of the weaker lines
! compute lineshapes, using lineparameters, at layers blah..kProflayer
    IF (iUseWeakBackGnd == 1) THEN
    !!! read in optical depths for the weak LTE lines
        write(kStdWarn,*) 'Adding in LTE contribution of weak lines ...'
    !--> this set computed at the US Standard profile plus 10 T
    !--> offsets, so it is literally a compressed database
        iUnCompressType = -2  !!! usual AIRS layers weak CO2 backgnd
        iLowerOrUpper = -1    !!! usual kLAYERS layers
        CALL read_std_optdepths_new(iGasID,raFreq,iStart, &
        raTAmt,raTTemp,raTPress,raTPartPress, &
        raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
        iProfileLayers,iSplineType,pProf, &
        raVertTemp,iVertTempSet,rFileStartFr, &
        iTag,iActualTag,iUnCompressType, &
        iLowerOrUpper,daaWeakOptDepth)
        DO iL = iStart,kProfLayer   !!!loop over layers
            DO iFr = 1,kMaxPts
                dK = daaWeakOptDepth(iFr,iL) * dble(rNLTEStrength)
                daaPlanckCoeff(iFr,iL)       = daaPlanckCoeff(iFr,iL)       + dK
                daaNLTEGasAbCoeff(iFr,iL)    = daaNLTEGasAbCoeff(iFr,iL)    + dK
                daaSumNLTEGasAbCoeff(iFr,iL) = daaSumNLTEGasAbCoeff(iFr,iL) + dK
            END DO    !loop over layers
        END DO
    END IF        !IF (iUseWeakBackGnd == 1) THEN

!**** therefore, at the end of this routine
!**** daaNLTEGasAbCoeff = daaSumNLTEGasAbCoeff
!**** and in essence, daaNLTEGasAbCoeff = daaPlanckCoeff also

    RETURN
    end SUBROUTINE lte_spectra

!************************************************************************
! This subroutine computes the LTE abs coeffs for the background lines,
!   for stuff "above" the KCARTA atmosphere ie in the stratospphere
! The weakest lines are simple : just read in precomputed US Standard optical
!   depths, use the uppermost layer (lowest pressure) and scale optical depths
!   accordingly
! The strong lines, for which run7co2.m does linemixing, are a little more
!   complicated, as the code uses the current profile to compute voigt optical
!   depths * birnbaum of these strong lines. These lines are in LTE so use the
!   current temperature profile
    SUBROUTINE lte_spectra_upper(iTag,iActualTag, &
    rFileStartFr,iNLTEStart,rCO2mult, &
    raFreq,iGasID,iLTEIn,iaNLTEBands,caaStrongLines, &
    caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
    pProf,raPressLevels,raLayerHeight,raThickness, &
    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
    raUpperPress,raUpperPartPress,raUpperTemp, &
    raUpperGasAmt,raUpperNLTETemp, &
    daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff, &
    daaUpperSumNLTEGasAbCoeff, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat, &
    iUpper,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
    dLineStrenMin,dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength, &
    iSetBloat,daFreqBloat, &
    iGas,iSplineType,iNumberUA_NLTEOut, &
    rFrLow,rFrHigh, &
    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
    iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks, &
    caOutUAfile,caOutUABloatFile)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTEStart tells where the non LTE atmosphere starts
! iTag      tells the current frequency step
! iLTEin    tells which set of data to use
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! raFreq        tells you the wavenumber array
! pProf         tells you the avg pressure in the layers (in mbar)
! dDeltaFreq,kBoxCarUse are for default (0.0025 spacing, 5 point boxcar) or not
! caaStrongLines      line param files associated with strong lines, in LTE
! rCO2mult  = tells how much input .rtp file wants the CO2 to be multiplied by
!             remeber this is as compared to kCO2ppmv
    REAL :: rCO2mult
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    INTEGER ::      iGas,iSplineType  !!!needed to uncompress database
    CHARACTER(80) :: caaStrongLines(kGasStore)
    DOUBLE PRECISION :: dDeltafreqNLTE,dLineStrenMin
    INTEGER :: iUSeWeakBackGnd,iSetBloat
    REAL :: pProf(kProfLayer),raPressLevels(kProfLayer+1)     !!in mb
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer) !!in meters
    INTEGER :: iaNLTEBands(kGasStore),iLTEIn,iGasID,iTag,iActualTag,iNLTEStart
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    REAL :: raFreq(kMaxPts),rNLTEstrength,rFileStartFr
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! these are just for dumping stuff to the UA file
    REAL :: rFrLow,rFrHigh
    INTEGER :: iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi
    INTEGER :: iNumNLTEGases,iaNLTEChunks(kGasStore)
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE
    CHARACTER(80) :: caOutUAFile,caOutUABloatFile

! output parameters
! iUpper    tells which layer to end NLTE calcs
! daaPlanckCoeff        tells how to modify the Planck function
! while the rest are the P,PP,T(lte),T(NLTE),Q for the upper atmosphere
    INTEGER :: iUpper,iNumberUA_NLTEOut
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    REAL :: pProfNLTE_upatm(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)

! local variables
    CHARACTER(80) :: caStrong
    INTEGER :: iFr,iL,iBand,iNum,iISO,iLineMixBand,iDefault,iType
    REAL :: raLTETemp(kProfLayer),raNLTEtemp(kProfLayer)
    DOUBLE PRECISION :: daK(kMaxPts),dK
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN),daIso(kHITRAN)
    DOUBLE PRECISION :: dLTE
    DOUBLE PRECISION :: daChi(kMaxPts),daChiBloat(kBloatPts)
    DOUBLE PRECISION :: daaWeakOptDepth(kMaxPts,kProfLayer)
    INTEGER :: iChi,iNumMixRatioLevs,iDo
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer),rX,rY
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! his reads in mix ratios from GND to 0.005 mb to 0.00005 mb
!! read in alt/press/temp/upper mix ratio from caaUpperMixRatio
!!  into raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs
!!
!! this is DIFFERENT from the UA NLTE compressed database!
!! which is saved into raUpper*_Std
    CALL MixRatio(iGasID,rCO2mult,iLTEIn,caaUpperMixRatio, &
    raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
    raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
    raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    CALL read_upperatm_lte_temperature( &
    iGasID,iNLTEStart,iLTEin,iBand,caaNLTETemp, &
    raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
    pProf,raPresslevels,raLayerHeight,raThickness, &
    iUpper,raUpperTemp,raUpperGasAmt,raUpperPress,raUpperPartPress, &
    raUpperPressLevels,raUpperThickness)

    IF (kNLTEOutUAOpen == -1) THEN
    !!! open the output file etc
        iType = +1
        CALL OpenUAFile(iType,iUpper,caOutUAFile,rFrLow,rFrHigh, &
        iDumpAllUASpectra,iDumpAllUARads, &
        iFileIDLo,iFileIDHi,iTag, &
        iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        IF (iDumpAllUARads > 0) THEN
            iNumberUA_NLTEOut = 1 + iUpper
        ELSEIF (iDumpAllUARads <= 0) THEN
            iNumberUA_NLTEOut = 1 + 1
        END IF

        IF ((iSetBloat > 0) .AND. (kBloatNLTEOutUAOpen == -1)) THEN
            CALL OpenBloatUAFile( &
            iType,iUpper,caOutUABloatFile,rFrLow,rFrHigh, &
            iDumpAllUARads, &
            iFileIDLo,iFileIDHi,iTag, &
            iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        END IF
    END IF

    IF   (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'Upper atmosphere pressure levels'
        write(kStdWarn,*) '  iL       P(iL)         P(iL+1)           Pav'
        write(kStdWarn,*) '-----------------------------------------------'
    END IF

    DO iL = 1,iUpper
        rX = raUpperPressLevels(iL) - raUpperPressLevels(iL+1)
        rY = log(raUpperPressLevels(iL)/raUpperPressLevels(iL+1))
        pProfNLTE_upatm(iL) = rX/rY
        IF   (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1)) THEN
            write(kStdWarn,*) iL,raUpperPressLevels(iL),raUpperPressLevels(iL+1), &
            pProfNLTE_upatm(iL)
        END IF
    END DO

    write(kStdWarn,*) '-----------------------------------------------'

! ---------------> first compute lineshape of stronger lines, at LTE
    caStrong = caaStrongLines(iLTEIn)

!!! read in line params for the strong lines that are NOT in NLTE; so
!!! the default assumption is that they are in LTE
    CALL read_stronglineLTE_lineparameters( &
    iGasID,caaaNLTEBands,iaNLTEBands,iLTEin,dLineStrenMin,caStrong, &
    iNum,daIso,daElower,daLineCenter,daJL,daJU,daPshift, &
    daStren296,daW_For,daW_self,daW_temp)

    IF (iGasID == 2) THEN
    !!!use cousin for these strong bands in LTE,
    !!!don't get fancy with linemix/birnbaum right now
        iLineMixBand = +1
    ELSE
        iLineMixBand = -1
    END IF
! compute lineshapes, using lineparameters, at layers 1 .. iUpper
    iDo = iUseWeakBackGnd
    IF (iDo > 0) THEN
        write(kStdWarn,*) 'Adding in LTE contribution of strong lines ...'
        DO iL = 1,iUpper   !!!loop over layers
            IF (iNum > 0) THEN
                dLTE  = raUpperTemp(iL)*1.0d0
                write (kStdWarn,*) 'LBL for UA Strong BackGnd (LTE) Lines : iL,T,Q = ',iL,sngl(dLTE),raUpperGasAmt(iL)
                CALL compute_lte_spectra_fast(iTag,iActualTag, &
                daK,raFreq,iGasID,iNum,daIso, &
                daElower,daLineCenter,daJL,daJU,daPshift, &
                daStren296,daW_For,daW_self,daW_temp,dLTE,iLineMixBand, &
                iL,raUpperGasAmt,raUpperTemp,raUpperPress,raUpperPartPress, &
                dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength)
            ELSE
                write(kStdWarn,*) '   no Strong Bands are in LTE!!! -----> ',iL
                DO iFr = 1,kMaxPts
                    daK(iFr) = 0.0d0
                END DO
            END IF
            DO iFr = 1,kMaxPts
                dK = daK(iFr) * dble(rNLTEStrength)
                daaUpperNLTEGasAbCoeff(iFr,iL)    = 0.0  !!! currently doing LTE
                daaUpperSumNLTEGasAbCoeff(iFr,iL) = dK   !!! but accumulate opt dep
                daaUpperPlanckCoeff(iFr,iL)       = dK
            END DO
        END DO    !loop over layers
    END IF

    write(kStdWarn,*) ' '

! ---------------> then read in optical depths of the weaker lines
! ---------------> which were computed at US standard profile, and readjust to
! ---------------> the current layering structure
! ---------------> note we use the UA kCARTA database

! compute lineshapes, using lineparameters, at layers blah..kProflayer
    iDo = iUseWeakBackGnd

    IF (iDo > 0) THEN
    !!! read in optical depths for the weak LTE lines
        write(kStdWarn,*) 'Adding in LTE contribution of UA weak lines ...'
        CALL read_std_optdepths_upper_UA( &
        iGas,iSplineType,raUpperPressLevels, &
        pProfNLTE_upatm,iUpper, &
        rFileStartFr,iTag,iActualTag, &
        iGasID,raFreq, &
        raUpperGasAmt,raUpperTemp,raUpperPress,raUpperPartPress, &
        iUpper,daaWeakOptDepth)
        DO iL = 1,iUpper   !!!loop over layers
            DO iFr = 1,kMaxPts
                dK = daaWeakOptDepth(iFr,iL) * dble(rNLTEStrength)
                daaUpperNLTEGasAbCoeff(iFr,iL) = daaUpperNLTEGasAbCoeff(iFr,iL)+dK
                daaUpperSumNLTEGasAbCoeff(iFr,iL) = &
                daaUpperSumNLTEGasAbCoeff(iFr,iL)+dK
                daaUpperPlanckCoeff(iFr,iL)    = daaUpperPlanckCoeff(iFr,iL) + dK
            END DO    !loop over layers
        END DO
    END IF

    IF (iSetBloat > 0) THEN
        CALL bloatUAstuff(daaUpperNLTEGasAbCoeff,daaUpperPlanckCoeff,iUpper, &
        raFreq,daFreqBloat, &
        daaUpperPlanckCoeffBloat,daaUpperNLTEGasAbCoeffBloat, &
        daaUpperSumNLTEGasAbCoeffBloat)
    END IF

    RETURN
    end SUBROUTINE lte_spectra_upper

!************************************************************************
! This subroutine computes the NLTE abs coeffs based on the lineparameters
! and the LTE parameters, for stuff in the KCARTA atmosphere
    SUBROUTINE nonlte_spectra(iTag,iActualTag,iLTEin, &
    iStart,iStart2350,raFreq, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp, &
    pProf,raPressLevels,raLayerHeight,raThickness,iProfileLayers, &
    raTAmt,raTTemp,raTPress,raTPartPress, &
    daaNLTEGasAbCoeff,daaPlanckCoeff,daaSumNLTEGasAbCoeff, &
    iDoUpperAtmNLTE,iAllLayersLTE, &
    dDeltafreqNLTE,rNLTEstrength,iDoVoigtChi, &
    iSetBloat,daaSumNLTEGasAbCoeffBloat,daaNLTEGasAbCoeffBloat, &
    daaPlanckCoeffBloat,daFreqBloat)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
! iTag       tells the current frequency step
! iLTEin     tells which set of data to use
! iStart     tells which layer to start NLTE calcs
! iStart2350 tells which layer to start NLTE calcs for strongest 4um CO2 band
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! raTTemp etc   tell you the LTE gas profile
! raFreq        tells you the wavenumber array
! pProf         tells you the avg pressure in the AIRS layers
! dDeltaFreq,kBoxcarUse are for default (0.0025 spacing, 5 point boxcar) or not
! the daaXBloat are matrices, at high resolution, done if iSetBloat > 0
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: dDeltafreqNLTE
    INTEGER :: iSetBloat,iProfileLayers
    REAL :: pProf(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer)
    INTEGER :: iaNLTEBands(kGasStore),iLTEIn,iSTart,iStart2350,iTag,iActualTag
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raFreq(kMaxPts),rNLTEstrength
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! output parameter
! daaNLTEGasAbCoeffCoeff  tells how to modify the abs coeff
! daaPlanckCoeff         tells how to modify the Planck function
! daaSumNLTEGasAbCoeff keeps track of the total NLTE gas abs coeff (JU == 9)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE

! local variables
    INTEGER :: iFr,iL,iBand,iGasID,iNum,iISO,iLineMixBand,iUpdateSumNLTE,iPrintTalk
    REAL :: raNLTEtemp(kProfLayer),raVibQFT(kProfLayer),raLTETemp(kProfLayer)
    DOUBLE PRECISION :: daPlanck(kMaxPts)
    DOUBLE PRECISION :: daK(kMaxPts)
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daChi(kMaxPts),daChiBloat(kBloatPts)
    DOUBLE PRECISION :: dVibCenter      !!!from D. Edwards NLTE files
    INTEGER :: iStartUse,iaJ_UorL(kHITRAN)
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)
          
! this is for bloating
    DOUBLE PRECISION :: daKBloat(kBloatPts),daPBloat(kBloatPts),d1,d2,d3

    DO iL = 1,kProfLayer
        raVibQFT(iL) = 1.0
    END DO

    DO iBand = 1,iaNLTEBands(iLTEIn)
    ! read in the lineshape parameters for the band
        CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
        iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
        daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,iLineMixBand,iDoVoigtChi)

    ! read in the NON LTE temperatures amd Vibrational Partition Fcns
        CALL read_nonlte_temperature(iGasID,iISO,iLTEin,iBand,caaNLTETemp, &
        pProf,raPressLevels,raLayerHeight,raThickness,iProfileLayers, &
        raTPress,raTPartPress,raTTemp,raTAmt,daJL,daJU, &
        iaJ_UorL,raLTETemp,raNLTETemp,raVibQFT,iAllLayersLTE,dVibCenter)

        IF ((iGasID == 2) .AND. (iISO == 1) .AND. &
        (nint(daJU(1)) == 9)) THEN
            iStartUse = iStart2350
            write(kStdWarn,*) ' .. strongest nlte band, istart = ',iStartUse
        ELSE
            iStartUse = iStart
            write(kStdWarn,*) ' .. weaker nlte band, istart = ',iStartUse
        ENDIF

    ! compute lineshapes, using lineparameters, at layers blah..kProflayer
        DO iL = iStartUse,kProfLayer
        ! loop over layers

            IF ((iBand == 1 .AND. iL == iStartUse) .AND. &
            (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1))) THEN
                iPrintTalk = +1
            ELSE
                iPrintTalk = -1
            END IF

            CALL compute_nlte_spectra_planck_fast(iBand, &
            iPrintTalk,iTag,iActualTag,daK,daPlanck,raFreq, &
            iGasID,iNum,iISO,daElower,daLineCenter, &
            daJL,daJU,iaJ_UorL,daPshift, &
            daStren296,daW_For,daW_self,daW_temp, &
            daJLowerQuantumRot,caJPQR,dVibCenter,iLineMixBand, &
            iL,raLTETemp,raNLTETemp,raVibQFT,iUpdateSumNLTE, &
            raTAmt,raTTemp,raTPress,raTPartPress,iAllLayersLTE, &
            dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength, &
            iSetBloat,daKBloat,daPBloat,daFreqBloat)

            DO iFr = 1,kMaxPts
                daaSumNLTEGasAbCoeff(iFr,iL)= daK(iFr)+daaSumNLTEGasAbCoeff(iFr,iL)
                daaNLTEGasAbCoeff(iFr,iL)   = daK(iFr) + daaNLTEGasAbCoeff(iFr,iL)
                daaPlanckCoeff(iFr,iL)      = daPlanck(iFr)+daaPlanckCoeff(iFr,iL)
            END DO

            IF (iSetBloat > 0) THEN
                DO iFr = 1,kBloatPts
                    daaSumNLTEGasAbCoeffBloat(iFr,iL) = daKBloat(iFr) + &
                    daaSumNLTEGasAbCoeffBloat(iFr,iL)
                    daaNLTEGasAbCoeffBloat(iFr,iL) = daKBloat(iFr) + &
                    daaNLTEGasAbCoeffBloat(iFr,iL)
                    daaPlanckCoeffBloat(iFr,iL)    = daPBloat(iFr) + &
                    daaPlanckCoeffBloat(iFr,iL)
                END DO
            END IF
                       
        END DO    !loop over layers iL = iStart,kProfLayer
    END DO      !loop over bands

    RETURN
    end SUBROUTINE nonlte_spectra

!************************************************************************
! this subroutine does the stuff high up in the STRATOSphere!
! assumes that the only contribution is from the VERY strong NLTE bands
! and nothing else from the weak bands or other gases!
! reads in some silly stuff for CO2
    SUBROUTINE UpHigh_Nonlte_Spectra(iTag,iActualTag, &
    iNLTEStart,raFreq,iGasID,iLTEIn, &
    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,rCO2mult, &
    pProf,raPressLevels,raLayerHeight,raThickness, &
    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
    raUpperPress,raUpperPartPress,raUpperTemp, &
    raUpperGasAmt,raUpperNLTETemp, &
    daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
    iUpper,iDoUpperAtmNLTE,iAllLayersLTE, &
    dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength, &
    iSetBloat,daFreqBloat, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat, &
    iGas,iSplineType,iNumberUA_NLTEOut, &
    rFrLow,rFrHigh,iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
    iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks, &
    caOutUAfile,caOutUABloatFile)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
! iNLTEStart tells where the non LTE atmosphere starts
! iTag      tells the current frequency step
! iLTEin    tells which set of data to use
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp   tells the name of the files containing the nonLTE temps
! raTTemp etc     tell you the LTE gas profile
! raFreq         tells you the wavenumber array
! pProf           tells you the avg pressure in the layers (in mbar)
! dDeltaFreq,kBoxCarUse are for default (0.0025 spacing, 5 point boxcar) or not
! rCO2mult  = tells how much input .rtp file wants the CO2 to be multiplied by
!             remeber this is as compared to kCO2ppmv
    REAL :: rCO2mult
    DOUBLE PRECISION :: dDeltafreqNLTE
    INTEGER :: iSetBloat
    REAL :: pProf(kProfLayer),raPressLevels(kProfLayer+1)     !!in mb
    REAL :: raLayerHeight(kProfLayer),raThickness(kProfLayer) !!in meters
    INTEGER :: iaNLTEBands(kGasStore),iLTEIn,iGasID,iTag,iActualTag,iNLTEStart
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    REAL :: raFreq(kMaxPts),rNLTEstrength
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
! these are just for dumping stuff to the UA file
    REAL :: rFrLow,rFrHigh
    INTEGER :: iDumpAllUARads,iDumpAllUASpectra,iFileIDLo,iFileIDHi
    INTEGER :: iNumNLTEGases,iaNLTEChunks(kGasStore)
    INTEGER :: iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE
    INTEGER :: iGas,iSplineType,iNumberUA_NLTEOut
    CHARACTER(80) :: caOutUAFile,caOutUABloatFile

! output parameters
! iUpper    tells which layer to end NLTE calcs
! daaPlanckCoeff        tells how to modify the Planck function
! while the rest are the P,PP,T(lte),T(NLTE),Q for the upper atmosphere
    INTEGER :: iUpper
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    REAL :: pProfNLTE_upatm(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
    INTEGER :: iFr,iL,iBand,iNum,iISO,iLineMixBand,iUpdateSumNLTE,iType
    REAL :: raVibQFT(kProfLayer)
    DOUBLE PRECISION :: daPlanck(kMaxPts)
    DOUBLE PRECISION :: daK(kMaxPts)
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daChi(kMaxPts),daChiBloat(kBloatPts)
    DOUBLE PRECISION :: dVibCenter      !from Dave Edwards NLTE profiles
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)
    INTEGER :: iChi,iNumMixRatioLevs,iaJ_UorL(kHITRAN)
    REAL :: raUpper_Pres(2*kProfLayer),raUpper_MixRatio(2*kProfLayer)

! this is for bloating
    DOUBLE PRECISION :: daKBloat(kBloatPts),daPBloat(kBloatPts)

! this is for dumping stuff
    INTEGER :: iIOUN,iPrintTalk
    CHARACTER(80) :: caOutName
    REAL :: raX(kMaxPts)

! this is mainly junk
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer),rJunk
    INTEGER :: iUpperStd_Num

    DO iL = 1,kProfLayer
        raVibQFT(iL) = 1.0
    END DO

! the file really should have been opened in SUBROUTINE nonlte_spectra
    IF (kNLTEOutUAOpen == -1) THEN
    !!! open the output file etc
        iType = +1
        CALL OpenUAFile(iType,iUpper,caOutUAFile,rFrLow,rFrHigh, &
        iDumpAllUASpectra,iDumpAllUARads, &
        iFileIDLo,iFileIDHi,iTag, &
        iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        IF (iDumpAllUARads > 0) THEN
            iNumberUA_NLTEOut = 1 + iUpper
        ELSEIF (iDumpAllUARads <= 0) THEN
            iNumberUA_NLTEOut = 1 + 1
        END IF

        IF ((iSetBloat > 0) .AND. (kBloatNLTEOutUAOpen == -1)) THEN
            CALL OpenBloatUAFile( &
            iType,iUpper,caOutUABloatFile,rFrLow,rFrHigh, &
            iDumpAllUARads, &
            iFileIDLo,iFileIDHi,iTag, &
            iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)
        END IF
    END IF

    IF (iGasID /= 2) THEN
        write(kStdWarn,*) 'Currently kCARTA has no info above 0.005 mb'
        write(kStdWarn,*) 'for GASID = ',iGasID
        CALL DoStop
    ELSEIF ((iGasID == 2) .AND. (iDoUpperAtmNLTE > 0)) THEN
    !! read in the mixing ratios (from files supplied by eg glatm.dat
    !! read in alt/press/temp/upper mix ratio from caaUpperMixRatio
    !!  into raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs
    !!
    !! this is DIFFERENT from the UA NLTE compressed database!
    !! which is saved into raUpper*_Std
        CALL MixRatio(iGasID,rCO2mult,iLTEIn,caaUpperMixRatio, &
        raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
        raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
        raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

        DO iBand = 1,iaNLTEBands(iLTEIn)
        ! do a ***dummy*** read to get the iL,iU gas quantum numbers
            CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
            iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
            daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,iLineMixBand,iDoVoigtChi)

        ! read in the upper atm NLTE temperatures
            CALL read_upperatm_nonlte_temperature( &
            iGasID,iISO,iNLTEStart,iLTEin,iBand,caaNLTETemp, &
            raUpper_Pres,raUpper_MixRatio,iNumMixRatioLevs, &
            pProf,raPresslevels,raLayerHeight,raThickness, &
            iUpper,raUpperTemp,raUpperGasAmt,raUpperPress,raUpperPartPress, &
            daJL,daJU,iaJ_UorL,raUpperNLTETemp,raVibQFT, &
            iAllLayersLTE,dVibCenter, &
            raUpperPressLevels,raUpperThickness, &
            raUpperPress_Std,raUpperMixRatio_Std,raUpperDZ_Std, &
            raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

            IF ((iUpper >= 1) .AND. (iDoUpperAtmNLTE > 0)) THEN
            ! read in the lineshape parameters
                CALL read_lineparameters(iLTEin,iBand,caaaNLTEBands, &
                iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,daPshift, &
                daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR,iLineMixBand,iDoVoigtChi)
            ! compute lineshapes, using lineparameters,
            ! at layers kProflayer+1 .. kProfLayer+iUpper
                DO iL = 1,iUpper

                    IF  ((iBand == 1 .AND. iL == 1) .AND. &
                    (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
                    (abs(kLongOrShort) <= 1))) THEN
                        iPrintTalk = +1
                    ELSE
                        iPrintTalk = -1
                    END IF

                    CALL compute_nlte_spectra_planck_fast(iBand, &
                    iPrintTalk,iTag,iActualTag,daK,daPlanck,raFreq, &
                    iGasID,iNum,iISO,daElower,daLineCenter, &
                    daJL,daJU,iaJ_UorL,daPshift, &
                    daStren296,daW_For,daW_self,daW_temp,daJLowerQuantumRot,caJPQR, &
                    dVibCenter,iLineMixBand, &
                    iL,raUpperTemp,raUpperNLTETemp,raVibQFT,iUpdateSumNLTE, &
                    raUpperGasAmt,raUpperTemp,raUpperPress,raUpperPartPress, &
                    iAllLayersLTE,dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength, &
                    iSetBloat,daKBloat,daPBloat,daFreqBloat)

                    DO iFr = 1,kMaxPts
                        daaUpperNLTEGasAbCoeff(iFr,iL) = &
                        daK(iFr) + daaUpperNLTEGasAbCoeff(iFr,iL)
                        daaUpperSumNLTEGasAbCoeff(iFr,iL) = &
                        daK(iFr) + daaUpperSumNLTEGasAbCoeff(iFr,iL)
                        daaUpperPlanckCoeff(iFr,iL)   = &
                        daPlanck(iFr) + daaUpperPlanckCoeff(iFr,iL)
                    END DO

                    IF (iSetBloat > 0) THEN
                        DO iFr = 1,kBloatPts
                            daaUpperSumNLTEGasAbCoeffBloat(iFr,iL) = daKBloat(iFr) + &
                            daaUpperSumNLTEGasAbCoeffBloat(iFr,iL)
                            daaUpperNLTEGasAbCoeffBloat(iFr,iL) = daKBloat(iFr) + &
                            daaUpperNLTEGasAbCoeffBloat(iFr,iL)
                            daaUpperPlanckCoeffBloat(iFr,iL) = daPBloat(iFr) + &
                            daaUpperPlanckCoeffBloat(iFr,iL)
                        END DO
                    END IF

                END DO    !loop over layers
            END IF      !if iUpper >= 1
        END DO
    END IF

    IF (iDumpAllUASpectra > 0) THEN
    ! dump out the UA spectra, iUpper paths worth of them!
        CALL wrtout_head_uafile(caOutUAFile, &
        raFreq(1),raFreq(kMaxPts),raFreq,iTag,1,iUpper)
        caOutName = 'DumDum'
        iIOUN = kNLTEOutUA
        DO iL = 1,iUpper
            rjunk = sngl(daaUpperSumNLTEGasAbCoeff(1,iL))
            write(kStdWarn,*) 'dump out UA ODs : gid, layer, amt, OD = ',iGasID,iL,raUpperGasAmt(iL),rJunk
            DO iFr = 1,kMaxPts
                raX(iFr) = sngl(daaUpperSumNLTEGasAbCoeff(iFr,iL))
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raX)
        END DO
    END IF

    RETURN
    end SUBROUTINE UpHigh_Nonlte_Spectra

!************************************************************************
! this subroutine does some inits for the running meshes
    SUBROUTINE InitRunningMesh(iTag,raFreq,dDeltaFreqNLTE, &
    dXNear,dXMedium,dXCoarse,iaClose,df,dfFine, &
    iOneCmFine,iFineMeshBoxPts,iOneCmMedium,iMediumMeshBoxPts, &
    iOneCmCoarse,iCoarseMeshBoxPts,iNWide,daK,daPlanck,daFreq, &
    iSetBloat,daKBloat,daPBloat)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iTag,iSetBloat
    REAL :: raFreq(kMaxPts)
    DOUBLE PRECISION :: dDeltaFreqNLTE
! output params
    DOUBLE PRECISION :: daFreq(kMaxPts),daK(kMaxPts),daPlanck(kMaxPts)
    DOUBLE PRECISION :: dXNear,dXMedium,dXCoarse,dfFine,f0,df
    DOUBLE PRECISION :: daFreqFineMesh(kMaxPtsBox),daFreqOutMesh(kMaxPts)
    DOUBLE PRECISION :: daKBloat(kBloatPts),daPBloat(kBloatPts)
    INTEGER :: iaClose(kHITRAN),iNWide
    INTEGER :: iOneCmFine,iOneCmMedium,iOneCmCoarse
    INTEGER :: iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts

! local vars
    INTEGER :: iFr,iJump

!!!!see run7.m for defn of xnear,xmed,xfar ... these are bins
    dXNear = 1.0d0
    dXMedium  = 2.0d0
    dXCoarse  = 1.0d2

    DO iFr = 1,kHITRAN
        iaClose(iFr)  = +1
    END DO

    IF (dDeltaFreqNLTE < 0) THEN
        df = kaFrStep(iTag)
        dfFine = kaFineFrStep(iTag)
    ELSE
        df = dDeltaFreqNLTE
        dfFine = df/(kBoxCarUse*1.0d0)
    END IF

!!!for the fine mesh
    iJump = iFloor(kBoxCarUse*1.0/2.0)
    iOneCmFine = nint(1.0d0/df)
    iFineMeshBoxPts = iOneCmFine*kBoxCarUse !!!compute at upto 5 times finer
    iFineMeshBoxPts = iFineMeshBoxPts+iJump !!!need 2 points of left of chunk
    iFineMeshBoxPts = iFineMeshBoxPts - 1   !!!but since go to eg 2381.9975
!!!instead of 2382.0000, need 1
!!!less point on the right,
!!!instead of 2 more points!

!!!for the medium mesh
    iOneCmMedium = nint(1.0d0/kaMediumFrStep(iTag))
    iMediumMeshBoxPts = iOneCmMedium + 1

!!!for the coarse mesh
    iOneCmCoarse = nint(1.0d0/kaCoarseFrStep(iTag))
    iCoarseMeshBoxPts = iOneCmCoarse + 1

! split 25 cm-1 into 25 chunks of 1 cm-1 width each
    iNWide = nint(raFreq(kMaxPts)-raFreq(1))

    IF (dDeltaFreqNLTE < 0.0d0) THEN
        df = 1.0d0 * kaFrStep(iTag)      !!!!use default
        iOneCmFine = 400                 !!!! 1 cm-1 = 400 pts at 0.0025 cm-1
        iOneCmFine = nint(1.0d0/df)
    ELSE
        df = dDeltaFreqNLTE
        iOneCmFine = nint(1.0d0/df)
    ! OneCmFine = 1000
    END IF

    f0 = raFreq(1)*1.0d0
    DO iFr = 1,kMaxPts
        daK(iFr)      = 0.0d0      !!!!!! ---------> initialise <------
        daPlanck(iFr) = 0.0d0      !!!!!! ---------> initialise <------
        daFreq(iFr)   = f0 + (iFr-1)*df
    END DO

    IF (iSetBloat > 0) THEN
        DO iFr = 1,kBloatPts
            daKBloat(iFr) = 0.0d0      !!!!!! ---------> initialise <------
            daPBloat(iFr) = 0.0d0      !!!!!! ---------> initialise <------
        END DO
    END IF

    RETURN
    end SUBROUTINE InitRunningMesh

!************************************************************************
! this sets up the running mesh!!!!
    SUBROUTINE SetRunningMesh(iWideMeshloop,daFreq,dXNear,dXMedium,dXCoarse, &
    iNum,daLineShift,i1,i2,iTag,iActualTag,dfFine, &
    iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts, &
    iDoFine,iDoMedium,iDoCoarse, &
    dX1,dX2,d1,d2,d3,d4,d5,d6,iCounter10,iCounter25,iCounter100, &
    iaClose,daFreqOutMesh,daMedium,daCoarse, &
    daFreqFineMesh,daTempClose,daTempMedium,daTempCoarse, &
    iISO,dJU,daFudgeF,daFudgeM,daFudgeC,iDoVoigtChi)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iDoFine,iDoMedium,iDoCoarse !!do we do fine,med,coarse?
    INTEGER :: iWideMeshLoop,iNum,i1,i2,iTag,iActualTag,iISO
    INTEGER :: iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts
    DOUBLE PRECISION :: dXNear,dXMedium,dXCoarse,dfFine,dJU
    DOUBLE PRECISION :: daFreq(kMaxPts),daLineShift(kHITRAN)
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! output params
    INTEGER :: iaClose(kHITRAN)
    INTEGER :: iCounter10,iCounter25,iCounter100
    DOUBLE PRECISION :: d1,d2,d3,d4,d5,d6,dX1,dX2
    DOUBLE PRECISION :: daTempClose(kMaxPtsBox)
    DOUBLE PRECISION :: daTempMedium(kMaxPtsBox),daTempCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: daFreqFineMesh(kMaxPtsBox),daFreqOutMesh(kMaxPts)
    DOUBLE PRECISION :: daMedium(kMaxPtsBox),daCoarse(kMaxPtsBox)
! fudge factors, if necessary
    DOUBLE PRECISION :: daFudgeF(kMaxPtsBox)
    DOUBLE PRECISION :: daFudgeM(kMaxPtsBox),daFudgeC(kMaxPtsBox)

! local variables
    INTEGER :: iLines,iFr,iN,iM
    DOUBLE PRECISION :: daFine(kMaxPtsBox)

    IF (iDoFine < 0) THEN
        write(kStdErr,*) 'Wow .. you have iDoFine = -1!!!'
        CALL DoStop
    END IF

    IF (iDoMedium < 0) THEN
        write(kStdWarn,*) 'Wow .. you have iDoMed = -1!!!'
    END IF

    IF (iDoCoarse < 0) THEN
        write(kStdWarn,*) 'Wow .. you have iDoCoarse = -1!!!'
    END IF

! this is always done, as we ALWAYS do the fine stuff
    dX1 = daFreq(1) + (iWideMeshloop - 1) * 1.0d0
    dX2 = daFreq(1) + (iWideMeshloop)     * 1.0d0

!!!find the lines that are within +/- dXNear cm-1 from dX1,dX2
    d1 = dX1 - dXNear
    d2 = dX2 + dXnear
     
    IF (iDoMedium > 0) THEN   !!!run7lbl style
    !!!find the lines that are slightly further away from dX1,dX2
        d3 = dX1 - dXMedium
        d4 = dX2 + dXMedium

    !!!find the lines that are even further away from dX1,dX2
        d5 = dX1 - dXCoarse
        d6 = dX2 + dXCoarse

    ELSE   !!!GENLN2 style, no medium strength bins
    !!!find the lines that are slightly further away from dX1,dX2
        d3 = dX1 + dXMedium
        d4 = dX2 - dXMedium

    !!!find the lines that are even further away from dX1,dX2
        d5 = dX1 - dXCoarse
        d6 = dX2 + dXCoarse
    END IF

    iCounter10 = 0
    iCounter25 = 0
    iCounter100 = 0
    DO iLines = 1,iNum
        IF ((daLineShift(iLines) >= d1) .AND. &
        (daLineShift(iLines) <= d2)) THEN
        !!!this line center lies within +/- dXNear of current 1cm-1 chunk
        !!! these are the CLOSE lines
            iaClose(iLines) = +1
            iCounter10 = iCounter10 + 1
        ELSEIF ((daLineShift(iLines) >= d3) .AND. &
            (daLineShift(iLines) < d1)) THEN
        !!!this line center lies just outside current kCARTA chunk
        !!! these are the MEDIUM lines
            iaClose(iLines) = 0
            iCounter25 = iCounter25 + 1
        ELSEIF ((daLineShift(iLines) > d2) .AND. &
            (daLineShift(iLines) <= d4)) THEN
        !!!this line center lies just outside current kCARTA chunk
        !!! these are the MEDIUM lines
            iaClose(iLines) = 0
            iCounter25 = iCounter25 + 1
        ELSE
        !!!this line center lies outside the current kCARTA chunk
        !!! these are the FAR lines
            iaClose(iLines) = -1
            iCounter100 = iCounter100 + 1
        END IF
    END DO

    iFr = 3
    daFreqFineMesh(iFr)   = nint(dX1)*1.0d0
    daFreqFineMesh(iFr-1) = nint(dX1)*1.0d0 - dfFine
    daFreqFineMesh(iFr-2) = nint(dX1)*1.0d0 - 2.0*dfFine

!!!! set the fine, medium and coarse grids up
    DO iFr = 1,iFineMeshBoxPts
        daFreqFineMesh(iFr) = nint(dX1)*1.0d0 + (iFr-3)*dfFine
        daFine(iFr)         = daFreqFineMesh(iFr)
        daTempClose(iFr)    = 0.0d0
    END DO
    DO iFr = 1,iMediumMeshBoxPts
        daMedium(iFr)       = nint(dX1)*1.0d0 +(iFr-1)*kaMediumFrStep(iTag)
        daTempMedium(iFr)   = 0.0d0
    END DO
    DO iFr = 1,iCoarseMeshBoxPts
        daCoarse(iFr)       = nint(dX1)*1.0d0 +(iFr-1)*kaCoarseFrStep(iTag)
        daTempCoarse(iFr)   = 0.0d0
    END DO
!!!set the current output mesh
    DO iFr = i1,i2
        daFreqOutMesh(iFr-i1+1) = daFreq(iFr)
    END DO

    IF ((nint(sngl(dJU)) == 9) .AND. (iISO == 1) .AND. &
    (iDoVoigtChi > 0)) THEN
    ! oh oh check to see if we need to adjust the linemix to bring it
    ! approximately equal to UMBC-LBL
    ! this parameter set in kcartamain/kcartabasic
        iN = kMaxPts
        IF ((daFreq(1) <= 2505.0d0) .AND. (daFreq(iN) >= 2355.0d0)) THEN
            CALL co2_4um_nlte_fudge(daFine,   daFudgeF, iFineMeshBoxPts)
            CALL co2_4um_nlte_fudge(daMedium, daFudgeM, iMediumMeshBoxPts)
            CALL co2_4um_nlte_fudge(daCoarse, daFudgeC, iCoarseMeshBoxPts)
        END IF
    END IF
      
    RETURN
    end SUBROUTINE SetRunningMesh

!************************************************************************
! this subroutine computes the lineshapes AND
! it also computes the contribution to the numerator, for the Planck coeff
! this is done for ONE band, for ONE layer
    SUBROUTINE compute_nlte_spectra_planck_fast( &
    iBand,iPrintTalk,iTag,iActualTag,daK,daPlanck,raFreq, &
    iGasID,iNum,iISO,daElower,daLineCenter,daJL,daJU,iaJ_UorL,daPshift, &
    daStren296,daW_For,daW_self,daW_temp, &
    daJLowerQuantumRot,caJPQR,dVibCenter,iLineMixBand, &
    iL,raLTETemp,raNLTEtemp,raVibQFT,iUpdateSumNLTE, &
    raTAmt,raTTemp,raTPress,raTPartPress,iAllLayersLTE, &
    dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength, &
    iSetBloat,daKBloat,daPBloat,daFreqBloat)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: daKBloat(kBloatPts),daPBloat(kBloatPts)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    INTEGER :: iPrintTalk  !!!if iTalk = +1 print all sorts of stuff
!!!if iTalk = -1 print header only if iPrintTalk = 1
    INTEGER :: iBand       !!! this is band# x out of N bands to do
    INTEGER :: iAllLayersLTE   !!! -1 upper layers in NLTE, +1 all layers LTE
    INTEGER :: iGasID,iNum,iISO,iL,iTag,iActualTag,iLineMixBand
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: daJLowerQuantumRot(kHITRAN)
    CHARACTER(1) ::      caJPQR(kHITRAN)
    DOUBLE PRECISION :: dVibCenter
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raFreq(kMaxPts),rNLTEstrength
    REAL :: raLTETemp(kProfLayer),raNLTEtemp(kProfLayer),raVibQFT(kProfLayer)
    DOUBLE PRECISION :: dDeltaFreqNLTE
    INTEGER :: iSetBloat,iaJ_UorL(kHITRAN)
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! output params ... the arrays are looped over ALL lines in band
    DOUBLE PRECISION :: daK(kMaxPts),daPlanck(kMaxPts)
    INTEGER :: iUpdateSumNLTE  !!! -1 if band is in LTE, +1 if band is in NLTE

! local variables
    INTEGER :: iFr,iLines,iI,iJ
    DOUBLE PRECISION :: dP,dPP,dGasAmt,dNLTE,dLTE,dLTEr1r2
    DOUBLE PRECISION :: daQtipsFactor(kHITRAN)
    DOUBLE PRECISION :: daRTop(kHITRAN),daRBot(kHITRAN)
    DOUBLE PRECISION :: daCFactor(kHITRAN),daKFactor(kHITRAN)
    DOUBLE PRECISION :: v0,dPartitionFcn,dMass
    DOUBLE PRECISION :: daBroad(kHITRAN),daLineShift(kHITRAN)
    DOUBLE PRECISION :: daStrenNLTE(kHITRAN)
    DOUBLE PRECISION :: dAlpha,dBeta,factor
! this is to spline on the far off lines!!!!
    INTEGER :: iDoCoarse,iDoMedium,iDoFine,iTrue
    INTEGER :: iCounter10,iCounter25,iCounter100
    DOUBLE PRECISION :: daTemp(kMaxPts),daFreq(kMaxPts),daTempP(kMaxPts)
    DOUBLE PRECISION :: d1,d2,d3,d4,d5,d6,dX1,dX2
    DOUBLE PRECISION :: daMedium(kMaxPtsBox),yaMedium(kMaxPtsBox)
    DOUBLE PRECISION :: daCoarse(kMaxPtsBox),yaCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: daTempClose(kMaxPtsBox),daTempClose1(kMaxPtsBox)
    DOUBLE PRECISION :: daTempMedium(kMaxPtsBox),daTempCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: daPlanckClose(kMaxPtsBox),daPlanckClose1(kMaxPtsBox)
    DOUBLE PRECISION :: daPlanckMedium(kMaxPtsBox),daPlanckCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: dXNear,dXMedium,dXCoarse,dfFine,f0,df
    DOUBLE PRECISION :: daFreqFineMesh(kMaxPtsBox),daFreqOutMesh(kMaxPts)
    INTEGER :: iCount,iaClose(kHITRAN),iNWide,iWideMeshLoop,i1,i2
    INTEGER :: iOneCmFine,iOneCmMedium,iOneCmCoarse
    INTEGER :: iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts
    DOUBLE PRECISION :: daFudgeF(kMaxPtsBox)
    DOUBLE PRECISION :: daFudgeM(kMaxPtsBox),daFudgeC(kMaxPtsBox)

    INTEGER :: iNum1,iTooFar,iTooFarX,iLine,iaLineMix(kHITRAN)
    INTEGER :: iDefault,iNoPressureShiftCO2
    DOUBLE PRECISION :: daYmix(kHitran)
    DOUBLE PRECISION :: dT,daChi(kMaxPtsBox),dVibQFT
    CHARACTER cCousOrBirn

! for birnbaum
    DOUBLE PRECISION :: chiBirn(kMaxPts),xBirn(kMaxPts)
    INTEGER :: iNptsBirn

! if we want to dump out line parameters !!!! ugh!!! makes warning.msg LONG!!!!
    INTEGER :: iTalk,iOneLine

    iDoFine = +1     !!!go ahead and loop over the "near" lines default
    iDoFine = -1     !!!ignore the "near" lines

    iDoMedium = +1     !!!go ahead and loop over the "med" lines default
    iDoMedium = -1     !!!ignore the "med" lines

    iDoCoarse = +1     !!!go ahead and loop over the "far" lines default
    iDoCoarse = -1     !!!ignore the "far" lines yeah!

    iDoFine = +1
    iDoMedium = +1
    iDoCoarse = +1

    iDefault = -1       !!!!just summarize things        DEFAULT
    iTalk = +1          !!!!oh no!!! give out all the info!!!
    iTalk = -1          !!!!just summarize things        DEFAULT
    IF (iTalk /= iDefault) THEN
        print *,'iTalk,iDefault = ',iTalk,iDefault
    END IF

    iDefault            = +1 !default for CO2 in the strong bands (linemix)
    iNoPressureShiftCO2 = -1 !allow p shifts for CO2 strong linemix lines
    iNoPressureShiftCO2 = +1 !default for CO2 in the strong bands (linemix)

    iNoPressureShiftCO2 = +1

    IF (iNoPressureShiftCO2 /= iDefault) THEN
        print *,'iNoPressureShiftCO2,iDefault = ',iNoPressureShiftCO2,iDefault
    END IF

!!!these are from kLAYERS
    dP      = raTPress(iL)*1.0d0
    dPP     = raTPartPress(iL)*1.0d0
    dGasAmt = raTAmt(iL)*1.0d0
    dLTE    = raTtemp(iL)*1.0d0         !!KLAYERS local kinetic temperature

!!! these are from the VT files
    IF (iAllLayersLTE == -1) THEN
        dLTEr1r2 = raLTEtemp(iL)*1.0d0     !!local kinetic temp, from VT file
        dNLTE    = raNLTEtemp(iL)*1.0d0    !!band vib temp, from VT file
    ELSEIF (iAllLayersLTE == +1) THEN
    ! before March 3, 2005
        dLTEr1r2 = raTtemp(iL)*1.0d0       !!KLAYERS local kinetic temperature
        dNLTE    = raTtemp(iL)*1.0d0       !!KLAYERS local kinetic temperature
    ! after March 3, 2005
        dLTEr1r2 = raLTEtemp(iL)*1.0d0     !!local kinetic temp, from VT file
        dNLTE    = raLTEtemp(iL)*1.0d0     !!local kinetic temp, from VT file
    END IF
    dVibQFT = raVibQFT(iL)*1.0d0

    IF (iGasID /= 2) THEN
    ! find the shifted line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines) + dP*daPshift(iLines)
        END DO
    ELSEIF ((iGasID == 2) .AND. (iNoPressureShiftCO2 == -1)) THEN
    ! find the shifted line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines) + dP*daPshift(iLines)
        END DO
    ELSEIF ((iGasID == 2) .AND. (iNoPressureShiftCO2 == +1)) THEN
    ! no pressure shifts for line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines)
        END DO
    END IF

    CALL InitRunningMesh(iTag,raFreq,dDeltaFreqNLTE, &
    dXNear,dXMedium,dXCoarse,iaClose,df,dfFine, &
    iOneCmFine,iFineMeshBoxPts,iOneCmMedium,iMediumMeshBoxPts, &
    iOneCmCoarse,iCoarseMeshBoxPts,iNWide,daK,daPlanck,daFreq, &
    iSetBloat,daKBloat,daPBloat)

    iUpdateSumNLTE = -1   !!! assume the band is in LTE
    IF (abs(dLTE - dNLTE) >= 5.0d-2) THEN
        iUpdateSumNLTE = +1
    END IF

    DO iFr = 1,kHITRAN
        daRTop(iFr)   = 1.0d0
        daRBot(iFr)   = 1.0d0
        iaClose(iFr)  = +1
    END DO

! find the qfcn for this GAS,ISOTOPE combination at LTE
! dPartition fcn = q296/qT
! this band has all molecules with same mass, isotope number
! so dPartitionFcn applies equally to all the lines!
    CALL qfcn(dLTE,iGasID,iISO,dPartitionFcn,dMass)

!!!find the broadening coefficient
    CALL Broad(iGasID,iNum,daLineShift, &
    daW_For,daW_Self,daW_temp,dP,dPP,dLTE,daBroad)

!!! check to see if we use Cousin or LineMix/Birnbaum or Nothing
    IF (iLineMixBand == 2) THEN
    !!!! check if we need Birn+Linemix or Cousin
    !!!! individual lines are done either with first order linemix, or with
    !!!! cousin. The reason is that eg P branch for 2350 band, gives sucky
    !!!! results if linemixing is used, while the R branch works fine!
        CALL CousinVsMix(iaLineMix,iLineMixBand, &
        int(daJL(1)),int(daJU(1)),iISO,daFreq, &
        daLineShift,daStren296,dVibCenter,iNum)
        CALL birnbaum_coarse(chiBirn,xBirn,iNptsBirn,dLTE,tau2_birn(dP,dPP))
    !      print *,'doing birnbaum_coarse ',iLineMixBand,int(daJL(1)),int(daJU(1)),iISO,daFreq(1)
    ELSE
        DO i1 = 1,iNum
            iaLineMix(i1) = iLineMixBand  !!stick to Cousin or nothing
        END DO
    END IF

! find the line strengths at *** LTE *****
! also find QTIPS modifier due to vibration partition function at NLTE
    CALL Strengths_nlte_vibpartfcn( &
    iL,iNum,iGasID,iISO, &
    dGasAmt,dPartitionFcn,dLTE,dNLTE,dVibQFT, &
    daLineShift,daELower,daStren296,daStrenNLTE,daQtipsFactor)

! find the ratios of upper and lower level populations, at
! the local temperature, and the non local vibrational temperature
! this is so that we can find the modification to the Planck function
    CALL NLTEPopulationRatios( &
    iL,iNum,dP,dLTE,dLTEr1r2,dNLTE,dVibCenter, &
    daElower,daLineShift,iaJ_UorL, &
    daRTop,daRBot,daCFactor,daKFactor)

! divide lines into FINE,MEDIUM,COARSE meshes
! lines to use. take [f3 f4] as the limits of the current wide mesh
! (a) fine lines are always within  [f3-xnear  f4+xnear]     xnear ~ 1 cm-1
!     ie these lines are always within +/- 1 cm-1 of current wide mesh
! (b) medium lines are always within
!             [f3-xmed f3-xnear] U [f4+xnear f4+xmed]      xmed ~ 2 cm-1
!     ie these lines are always 1-2 cm-1 away from current wide mesh
! (c) far lines are always within
!             [f3-xfar f3-xmed] U [f4+xmed f4+xfar]        xfar ~ 25 cm-1
!     ie these lines are always 2-25 cm-1 away from current wide mesh

! used to be 2317.2
    iOneLine = -1
    IF (iOneLine > 0) THEN
        print *,'iOneline = 2324.97'
        DO iLines = 1,iNum
            IF (dabs(daLineShift(iLines)-2324.97) > 0.5d0) THEN
                daStrenNLTE(iLines) = 0.0d0
            END IF
        END DO
    END IF

    i1 = 1
    i2 = iOneCmFine
    DO iWideMeshLoop = 1,iNWide
        CALL SetRunningMesh(iWideMeshloop,daFreq,dXNear,dXMedium,dXCoarse, &
        iNum,daLineShift,i1,i2,iTag,iActualTag,dfFine, &
        iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts, &
        iDoFine,iDoMedium,iDoCoarse, &
        dX1,dX2,d1,d2,d3,d4,d5,d6,iCounter10,iCounter25,iCounter100, &
        iaClose,daFreqOutMesh,daMedium,daCoarse, &
        daFreqFineMesh,daTempClose,daTempMedium,daTempCoarse, &
        iISO,daJU(1),daFudgeF,daFudgeM,daFudgeC,iDoVoigtChi)

        DO iFr = 1,iFineMeshBoxPts
            daPlanckClose(iFr) = 0.0d0
        END DO
        DO iFr = 1,iMediumMeshBoxPts
            daPlanckMedium(iFr) = 0.0d0
        END DO
        DO iFr = 1,iCoarseMeshBoxPts
            daPlanckCoarse(iFr) = 0.0d0
        END DO

    ! now loop over the close lines to get the optical depths!
        IF (iDoFine > 0) THEN
            iTooFarX = +10000
            iTrue = -1
            DO iLines = 1,iNum
                IF (iaClose(iLines) > 0) THEN
                ! do the "close" lines
                    iTrue = +1
                    CALL AlphaBeta_Factors( &
                    iLines,iUpdateSUMNLTE,rNLTEstrength,daRTop,daRBot, &
                    daCFactor,daKFactor,dAlpha,dBeta)
                    CALL voigt_chi( &
                    daFreqFineMesh,iFineMeshBoxPts,daLineShift(iLines), &
                    dLTE,dMass,daBroad(iLines),iaLineMix(iLines),dP,dPP, &
                    daJL(iLines),daJU(iLines),daJLowerQuantumRot(iLines),iISO,dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,daTempClose1,daFudgeF,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iFineMeshBoxPts
                        daTempClose(iFr) = daTempClose(iFr) + &
                        dAlpha*daStrenNLTE(iLines)*daTempClose1(iFr)
                        daPlanckClose(iFr) = daPlanckClose(iFr) + &
                        dBeta*daStrenNLTE(iLines)*daTempClose1(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                print *,' compute_nlte_spectra_planck_fast FINE lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                print *,' compute_nlte_spectra_planck_fast FINE lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!! now do the boxcar
            IF (iTrue > 0) THEN
                CALL DoBoxCar(daTempClose,iFineMeshBoxPts,i1,i2,daK)
                CALL DoBoxCar(daPlanckClose,iFineMeshBoxPts,i1,i2,daPlanck)
            END IF

            IF (iSetBloat > 0) THEN
                CALL AccumulateBloat(daFreqFineMesh,daTempClose,iFineMeshBoxPts, &
                daMedium      ,daTempMedium,iCount, &
                iWideMeshLoop,i1,i2,1,daKBloat)
                CALL AccumulateBloat(daFreqFineMesh,daPlanckClose,iFineMeshBoxPts, &
                daMedium      ,daPlanckMedium,iCount, &
                iWideMeshLoop,i1,i2,1,daPBloat)
            END IF
        END IF

        IF (iDoMedium > 0) THEN
        ! ah well; bite the bullet and do the med lines, even though they
        ! probably contribute nuthin!
            iTooFarX = +10000
            iCount = iMediumMeshBoxPts
            DO iLines = 1,iNum
                IF (iaClose(iLines) == 0) THEN
                ! do the "med" lines
                    CALL AlphaBeta_Factors( &
                    iLines,iUpdateSUMNLTE,rNLTEstrength,daRTop,daRBot, &
                    daCFactor,daKFactor,dAlpha,dBeta)
                    CALL voigt_chi( &
                    daMedium,iCount,daLineShift(iLines),dLTE,dMass, &
                    daBroad(iLines),iaLineMix(iLines),dP,dPP, &
                    daJL(iLines),daJU(iLines),daJLowerQuantumRot(iLines), &
                    iISO,dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,yaMedium,daFudgeM,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iCount
                        daTempMedium(iFr) = daTempMedium(iFr) + &
                        dAlpha*daStrenNLTE(iLines)*yaMedium(iFr)
                        daPlanckMedium(iFr) = daPlanckMedium(iFr) + &
                        dBeta*daStrenNLTE(iLines)*yaMedium(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                print *,' compute_nlte_spectra_planck_fast MEDIUM lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                print *,' compute_nlte_spectra_planck_fast MEDIUM lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!interpolate the med lines sum onto the near lines
            CALL dspl(daMedium,daTempMedium,iCount,daFreqOutMesh,daTemp,i2-i1+1)
            CALL dspl(daMedium,daPlanckMedium,iCount,daFreqOutMesh, &
            daTempP,i2-i1+1)
        !!!add on the med lines to the near lines
            DO iFr = i1,i2
                daK(iFr)      = daK(iFr)      + daTemp(iFr-i1+1)
                daPlanck(iFr) = daPlanck(iFr) + daTempP(iFr-i1+1)
            END DO
            IF (iSetBloat > 0) THEN
                CALL AccumulateBloat(daFreqFineMesh,daTempClose,iFineMeshBoxPts, &
                daMedium      ,daTempMedium,iCount, &
                iWideMeshLoop,i1,i2,-1,daKBloat)
                CALL AccumulateBloat(daFreqFineMesh,daPlanckClose,iFineMeshBoxPts, &
                daMedium      ,daPlanckMedium,iCount, &
                iWideMeshLoop,i1,i2,-1,daPBloat)
            END IF
        END IF

        IF (iDoCoarse > 0) THEN
        ! ah well; bite the bullet and do the far lines, even though they
        ! probably contribute nuthin!
            iTooFarX = +10000
            iCount = iCoarseMeshBoxPts
            DO iLines = 1,iNum
                IF (iaClose(iLines) < 0) THEN
                ! do the "far" lines
                    CALL AlphaBeta_Factors( &
                    iLines,iUpdateSUMNLTE,rNLTEstrength,daRTop,daRBot, &
                    daCFactor,daKFactor,dAlpha,dBeta)
                    CALL voigt_chi( &
                    daCoarse,iCount,daLineShift(iLines),dLTE,dMass, &
                    daBroad(iLines),iaLineMix(iLines), &
                    dP,dPP,daJL(iLines),daJU(iLines),daJLowerQuantumRot(iLines),iISO,dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,yaCoarse,daFudgeC,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iCount
                        daTempCoarse(iFr) = daTempCoarse(iFr) + &
                        dAlpha*daStrenNLTE(iLines)*yaCoarse(iFr)
                        daPlanckCoarse(iFr) = daPlanckCoarse(iFr) + &
                        dBeta*daStrenNLTE(iLines)*yaCoarse(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                print *,' compute_nlte_spectra_planck_fast COARSE lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                print *,' compute_nlte_spectra_planck_fast COARSE lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!interpolate the far lines sum onto the near lines
            CALL dspl(daCoarse,daTempCoarse,iCount,daFreqOutMesh,daTemp,i2-i1+1)
            CALL dspl(daCoarse,daPlanckCoarse,iCount,daFreqOutMesh, &
            daTempP,i2-i1+1)
        !!!add on the far lines to the near lines
            DO iFr = i1,i2
                daK(iFr)      = daK(iFr)      + daTemp(iFr-i1+1)
                daPlanck(iFr) = daPlanck(iFr) + daTempP(iFr-i1+1)
            END DO
            IF (iSetBloat > 0) THEN
                CALL AccumulateBloat(daFreqFineMesh,daTempClose,iFineMeshBoxPts, &
                daCoarse      ,daTempCoarse,iCount, &
                iWideMeshLoop,i1,i2,-1,daKBloat)
                CALL AccumulateBloat(daFreqFineMesh,daPlanckClose,iFineMeshBoxPts, &
                daCoarse      ,daPlanckCoarse,iCount, &
                iWideMeshLoop,i1,i2,-1,daPBloat)
            END IF
        END IF

        i1 = i1 + iOneCmFine
        i2 = i2 + iOneCmFine
    ENDDO

! now make sure everything > 0
    DO iFr = 1,kMaxPts
        daK(iFr)      = max(0.0d0,daK(iFr))
        daPlanck(iFr) = max(0.0d0,daPlanck(iFr))
    END DO

    IF (iSetBloat > 0) THEN
        DO iFr = 1,kBloatPts
            daKBloat(iFr)      = max(0.0d0,daKBloat(iFr))
            daPBloat(iFr)      = max(0.0d0,daPBloat(iFr))
        END DO
    END IF

    IF ((iBand == 1) .AND. &
    ( ((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1))) THEN
        iTalk = +1
    ELSE
        iTalk = -1
    END IF

    IF (iPrintTalk > 0) THEN
        write(kStdWarn,*)'doing LBL for Strong NLTE Lines :'
        write(kStdWarn,*)'lay Band   Line C/B  v0         r1         r2            Kfac    1-Cfac/Kfac'
        write(kStdWarn,*)'   center  indx      '
        write(kStdWarn,*)'------------------------------------------------------------------------------'
        DO iLines = 1,iNum
            IF (iaLineMix(iLines) == 1) THEN        !cousin
                cCousOrBirn = 'C'
            ELSEIF (iaLineMix(iLines) == 2) THEN    !linemix * birn
                cCousOrBirn = 'B'
            ELSEIF (iaLineMix(iLines) < 0) THEN    !nothing ... plain voigt
                cCousOrBirn = 'V'
            ELSE
                write(kStdErr,*) 'error in ialineMix',iLines,iaLineMix(iLines)
                CALL DoStop
            END IF
            IF (iLines == 1) THEN
                write(kStdWarn,1010) iL,dVibCenter,iLines,cCousOrBirn,daLineShift(iLines), &
                daRBot(iLines),daRTop(iLines),daKFactor(iLines), &
                &          1.0d0 - daCFactor(iLines)/daKFactor(iLines)
            END IF
        END DO
    END IF
          
    1010 FORMAT(I3,' ',D12.6,' ',I4,' ',A1,' ',1(D12.6,' '),5(D12.6,' '))

    RETURN
    end SUBROUTINE compute_nlte_spectra_planck_fast

!************************************************************************
! this subroutine computes the background lineshapes
! loops over ALL lines >= strmin, using LTE
! does it like run7.m where it divides the 25 cm-1 chunk into 1 cm-1 chunks
    SUBROUTINE compute_lte_spectra_fast( &
    iTag,iActualTag,daK,raFreq,iGasID,iNum,daIso, &
    daElower,daLineCenter,daJL,daJU,daPshift, &
    daStren296,daW_For,daW_self,daW_temp,dVibCenter,iLineMixBand, &
    iL,raTAmt,raTTemp,raTPress,raTPartPress, &
    dDeltaFreqNLTE,iDoVoigtChi,rNLTEstrength)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    INTEGER :: iGasID,iNum,iL,iTag,iActualTag,iLineMixBand
    DOUBLE PRECISION :: daElower(kHITRAN),daLineCenter(kHITRAN),daISO(kHITRAN)
    DOUBLE PRECISION :: daJL(kHITRAN),daJU(kHITRAN)
    DOUBLE PRECISION :: daPshift(kHITRAN),daStren296(kHITRAN),daW_for(kHITRAN)
    DOUBLE PRECISION :: daW_self(kHITRAN),daW_temp(kHITRAN)
    DOUBLE PRECISION :: dVibCenter,dDeltaFreqNLTE
    REAL :: raTTemp(kProfLayer),raTAmt(kProfLayer),raTPress(kProfLayer)
    REAL :: raTPartPress(kProfLayer)
    REAL :: raFreq(kMaxPts),rNLTEstrength
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! output params
    DOUBLE PRECISION :: daK(kMaxPts)

! local variables
    INTEGER :: iFr,iLines,iI,iJ,iTooFar,iTooFarX
    DOUBLE PRECISION :: dP,dPP,dGasAmt,dLTE
    DOUBLE PRECISION :: v0,dPartitionFcn,dMass
    DOUBLE PRECISION :: daBroad(kHITRAN),daLineShift(kHITRAN)
    DOUBLE PRECISION :: daStrenLTE(kHITRAN)
    DOUBLE PRECISION :: daPartitionFcn(kHITRAN),daMass(kHITRAN)
    DOUBLE PRECISION :: daQtipsFactor(kHITRAN)

    DOUBLE PRECISION :: daKBloat(kBloatPts),daPBloat(kBloatPts)
    INTEGER :: iZZZBloat

! this is to spline on the far off lines!!!!
    INTEGER :: iDoCoarse,iDoMedium,iDoFine,iTrue,iISO
    INTEGER :: iCounter10,iCounter25,iCounter100
    DOUBLE PRECISION :: daTemp(kMaxPts),daFreq(kMaxPts)
    DOUBLE PRECISION :: d1,d2,d3,d4,d5,d6,dX1,dX2
    DOUBLE PRECISION :: daPlanck(kMaxPts)   !!not really needed here
    DOUBLE PRECISION :: daMedium(kMaxPtsBox),yaMedium(kMaxPtsBox)
    DOUBLE PRECISION :: daCoarse(kMaxPtsBox),yaCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: daTempClose(kMaxPtsBox),daTempClose1(kMaxPtsBox)
    DOUBLE PRECISION :: daTempMedium(kMaxPtsBox),daTempCoarse(kMaxPtsBox)
    DOUBLE PRECISION :: dXNear,dXMedium,dXCoarse,dfFine,f0,df
    DOUBLE PRECISION :: daFreqFineMesh(kMaxPtsBox),daFreqOutMesh(kMaxPts)
    INTEGER :: iCount,iaClose(kHITRAN),iNWide,iWideMeshLoop,i1,i2
    INTEGER :: iOneCmFine,iOneCmMedium,iOneCmCoarse
    INTEGER :: iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts
    DOUBLE PRECISION :: daFudgeF(kMaxPtsBox)
    DOUBLE PRECISION :: daFudgeM(kMaxPtsBox),daFudgeC(kMaxPtsBox)

! for birnbaum
    DOUBLE PRECISION :: chiBirn(kMaxPts),xBirn(kMaxPts)
    INTEGER :: iNptsBirn,iDefault,iNoPressureShiftCO2,iOneLine

    iZZZBloat = -1

    iDoFine = +1     !!!go ahead and loop over the "near" lines default
    iDoFine = -1     !!!ignore the "near" lines

    iDoMedium = +1     !!!go ahead and loop over the "med" lines default
    iDoMedium = -1     !!!ignore the "med" lines

    iDoCoarse = +1     !!!go ahead and loop over the "far" lines default
    iDoCoarse = -1     !!!ignore the "far" lines yeah!

    iDoFine = +1
    iDoMedium = +1
    iDoCoarse = +1

    iDefault            = +1 !default for CO2 in the strong bands (linemix)
    iNoPressureShiftCO2 = -1 !allow p shifts for CO2 strong linemix lines
    iNoPressureShiftCO2 = +1 !default for CO2 in the strong bands (linemix)

    iNoPressureShiftCO2 = +1

    IF (iNoPressureShiftCO2 /= iDefault) THEN
        print *,'iNoPressureShiftCO2,iDefault = ',iNoPressureShiftCO2,iDefault
    END IF

    dP  = raTPress(iL)*1.0d0
    dPP = raTPartPress(iL)*1.0d0
    dGasAmt  = raTAmt(iL)*1.0d0
    dLTE     = raTtemp(iL)*1.0d0             !!local kinetic temperature

    IF (iGasID /= 2) THEN
    ! find the shifted line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines) + dP*daPshift(iLines)
        END DO
    ELSEIF ((iGasID == 2) .AND. (iNoPressureShiftCO2 == -1)) THEN
    ! find the shifted line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines) + dP*daPshift(iLines)
        END DO
    ELSEIF ((iGasID == 2) .AND. (iNoPressureShiftCO2 == +1)) THEN
    ! no pressure shifts for line centers
        DO iLines = 1, iNum
            daLineShift(iLines) = daLineCenter(iLines)
        END DO
    END IF

    CALL InitRunningMesh(iTag,raFreq,dDeltaFreqNLTE, &
    dXNear,dXMedium,dXCoarse,iaClose,df,dfFine, &
    iOneCmFine,iFineMeshBoxPts,iOneCmMedium,iMediumMeshBoxPts, &
    iOneCmCoarse,iCoarseMeshBoxPts,iNWide,daK,daPlanck,daFreq, &
    iZZZBloat,daKBloat,daPBloat)

    dP  = raTPress(iL)*1.0d0
    dPP = raTPartPress(iL)*1.0d0
    dGasAmt  = raTAmt(iL)*1.0d0
    dLTE     = raTtemp(iL)*1.0d0        !!lKLAYERS local kinetic temperature

! find the shifted line centers and loop over lines, for this layer
    DO iLines = 1, iNum
    ! find the qfcn for this GAS,ISOTOPE combination at LTE
        CALL qfcn(dLTE,iGasID,int(daISO(iLines)),dPartitionFcn,daMass(iLines))
        daPartitionFcn(iLines) = dPartitionFcn
    END DO

! find the line strengths at *** LTE *****
! daQtipsFactor(iN) = 1.0
    CALL Strengths_lte_only(iNum,iGasID,dGasAmt,daPartitionFcn,dLTE, &
    daLineShift,daELower,daStren296,daStrenLTE,daQtipsFactor)

    iOneLine = -1
    IF (iOneLine > 0) THEN
        print *,'iOneline = 2324.97'
        DO iLines = 1,iNum
            IF (dabs(daLineShift(iLines)-2324.97) > 0.5d0) THEN
                daStrenLTE(iLines) = 0.0d0
            END IF
        END DO
    END IF

!!!find the broadening coefficients
    CALL Broad(iGasID,iNum,daLineShift, &
    daW_For,daW_Self,daW_temp,dP,dPP,dLTE,daBroad)

! divide lines into FINE,MEDIUM,COARSE meshes
! lines to use. take [f3 f4] as the limits of the current wide mesh
! (a) fine lines are always within  [f3-xnear  f4+xnear]     xnear ~ 1 cm-1
!     ie these lines are always within +/- 1 cm-1 of current wide mesh
! (b) medium lines are always within
!             [f3-xmed f3-xnear] U [f4+xnear f4+xmed]      xmed ~ 2 cm-1
!     ie these lines are always 1-2 cm-1 away from current wide mesh
! (c) far lines are always within
!             [f3-xfar f3-xmed] U [f4+xmed f4+xfar]        xfar ~ 25 cm-1
!     ie these lines are always 2-25 cm-1 away from current wide mesh

    i1 = 1
    i2 = iOneCmFine
    DO iWideMeshLoop = 1,iNWide
        CALL SetRunningMesh(iWideMeshloop,daFreq,dXNear,dXMedium,dXCoarse, &
        iNum,daLineShift,i1,i2,iTag,iActualTag,dfFine, &
        iFineMeshBoxPts,iMediumMeshBoxPts,iCoarseMeshBoxPts, &
        iDoFine,iDoMedium,iDoCoarse, &
        dX1,dX2,d1,d2,d3,d4,d5,d6,iCounter10,iCounter25,iCounter100, &
        iaClose,daFreqOutMesh,daMedium,daCoarse, &
        daFreqFineMesh,daTempClose,daTempMedium,daTempCoarse, &
        int(daISO(iLines)),dAJU(1),daFudgeF,daFudgeM,daFudgeC,iDoVoigtChi)

    ! now loop over the close lines to get the optical depths!
        IF (iDoFine > 0) THEN
            iTrue = -1
            iTooFarX = +10000
            DO iLines = 1,iNum
                IF (iaClose(iLines) > 0) THEN
                ! do the "close" lines
                    iTrue = +1
                    CALL voigt_chi( &
                    daFreqFineMesh,iFineMeshBoxPts,daLineShift(iLines), &
                    dLTE,daMass(iLines),daBroad(iLines),iLineMixBand, &
                    dP,dPP,daJL(iLines),daJU(iLines),-1.0d0,int(daISO(iLines)),dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,daTempClose1,daFudgeF,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iFineMeshBoxPts
                        daTempClose(iFr) = daTempClose(iFr) + &
                        daStrenLTE(iLines)*daTempClose1(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast FINE lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast FINE lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!! now do the boxcar
            IF (iTrue > 0) THEN
                CALL DoBoxCar(daTempClose,iFineMeshBoxPts,i1,i2,daK)
            END IF
        END IF

        IF (iDoMedium > 0) THEN
        ! ah well; bite the bullet and do the med lines, even though they
        ! probably contribute nuthin!
            iTooFarX = +10000
            iCount = iMediumMeshBoxPts
            DO iLines = 1,iNum
                IF (iaClose(iLines) == 0) THEN
                ! do the "med" lines
                    CALL voigt_chi( &
                    daMedium,iCount,daLineShift(iLines),dLTE, &
                    daMass(iLines),daBroad(iLines),iLineMixBand,dP,dPP, &
                    daJL(iLines),daJU(iLines),-1.0d0,int(daISO(iLines)),dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,yaMedium,daFudgeM,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iCount
                        daTempMedium(iFr) = daTempMedium(iFr) + &
                        daStrenLTE(iLines)*yaMedium(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast MEDIUM lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast MEDIUM lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!interpolate the med lines sum onto the near lines
            CALL dspl(daMedium,daTempMedium,iCount,daFreqOutMesh,daTemp,i2-i1+1)
        !!!add on the med lines to the near lines
            DO iFr = i1,i2
                daK(iFr) = daK(iFr) + daTemp(iFr-i1+1)
            END DO
        END IF

        IF (iDoCoarse > 0) THEN
        ! ah well; bite the bullet and do the far lines, even though they
        ! probably contribute nuthin!
            iTooFarX = +10000
            iCount = iCoarseMeshBoxPts
            DO iLines = 1,iNum
                IF (iaClose(iLines) < 0) THEN
                ! do the "far" lines
                    CALL voigt_chi( &
                    daCoarse,iCount,daLineShift(iLines),dLTE, &
                    daMass(iLines),daBroad(iLines),iLineMixBand,dP,dPP, &
                    daJL(iLines),daJU(iLines),-1.0d0,int(daISO(iLines)),dVibCenter, &
                    xBirn,chiBirn,iNptsBirn,yaCoarse,daFudgeC,iDoVoigtChi,iTooFar)
                    IF ((iTooFar > 0) .AND. (iTooFar < iTooFarX)) THEN
                        iTooFarX = iTooFar
                    END IF
                    DO iFr = 1,iCount
                        daTempCoarse(iFr) = daTempCoarse(iFr) + &
                        daStrenLTE(iLines)*yaCoarse(iFr)
                    END DO
                END IF
            END DO
            IF (iTooFarX < 100) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast COARSE lines dVibCenter ---- bad lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            ELSEIF ((iTooFarX >= 100) .AND. (iTooFarX /= 10000)) THEN
                iISO = int(daISO(1))
                print *,' compute_lte_spectra_fast COARSE lines dVibCenter ---- OK lines Jrotatation(iTooFar) = ',iTooFarX
                print *,'   iNumLines  : vibctr,JL,JU,iISO = ',iNum,dVibCenter,int(daJL(1)),int(daJU(1)),iISO
            END IF
        !!!interpolate the far lines sum onto the near lines
            CALL dspl(daCoarse,daTempCoarse,iCount,daFreqOutMesh,daTemp,i2-i1+1)
        !!!add on the far lines to the near lines
            DO iFr = i1,i2
                daK(iFr) = daK(iFr) + daTemp(iFr-i1+1)
            END DO
        END IF

        i1 = i1 + iOneCmFine
        i2 = i2 + iOneCmFine
    ENDDO

! now make sure everything > 0
    DO iFr = 1,kMaxPts
        daK(iFr) = max(0.0d0,daK(iFr))
    END DO

    RETURN
    end SUBROUTINE compute_lte_spectra_fast

!************************************************************************
! allows a NLTE atmosphere, does SLOW ACCURATE GENLN2 calc

! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles, &
    iTag,raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer),rCO2MixRatio
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2,rOmegaSun
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    INTEGER :: iCloudLayerTop,iCloudLayerBot
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,troplayer
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm), &
        raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar

    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO

! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'slow nlte radmodel ...'
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh=iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay)  = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
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
    IF (iNLTEStart > kProfLayer) THEN
        iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
        write (kStdErr,*) 'Normal rad transfer .... no NLTE'
        write (kStdErr,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
        write (kStdErr,*) 'should be calling rad_trans_SAT_LOOK_DOWN'
        CALL DoStop
    ELSE
        iLay = 1
        987 CONTINUE
        iL=iaRadLayer(iLay)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF ((iLModKprofLayer < iNLTEStart) .AND. (iLay < iNumLayer)) THEN
            iLay = iLay + 1
            GOTO 987
        END IF
        iSTopNormalRadTransfer = iLay
        write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
    END IF

    DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF (iLModKprofLayer < iNLTEStart) THEN
        ! ormal, no LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iLModKprofLayer >= iNLTEStart) THEN
        ! ew; LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    DO iFr=1,kMaxPts
        raInten(iFr)=raSurface(iFr)*raUseEmissivity(iFr)+ &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
        raSun(iFr)*raSunRefl(iFr)
    END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    DO iLay=iHigh,iHigh
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

        IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp >= 1) THEN

                write(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer'
                write(kStdWarn,*) 'This is the top of the usual AIRS atmosphere'
                write(kStdWarn,*) '   you have iDoUpperATM > 0'
                write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
                write(kStdWarn,*) 'and output stuff into the blah_UA file'
                write(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of'
                write(kStdWarn,*) 'stratosphere into both this and the UA file'

            ! o radiative transfer thru this layer
                DO iFr=1,kMaxPts
                    raInten(iFr) = &
                    raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
                END DO

            ! ow do complete rad transfer thru upper part of atmosphere
                CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)), &
                iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of stratosphere (120-200 km)

                write(kStdWarn,*) 'finally outputting radiances at TOTAL Complete TOA into'
                write(kStdWarn,*) 'usual binary file (iLay = ',iLay,')'

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
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        END IF

    !c no need to do radiative transfer thru this layer
    !c        DO iFr=1,kMaxPts
    !c          raInten(iFr)=raaEmission(iFr,iLay)+
    !c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
    !c          END DO

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_SLOW

!************************************************************************
! allows a NLTE atmosphere, does FAST SARTA NLTE calc

! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FAST(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles, &
    iTag,raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2,rOmegaSun
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    INTEGER :: iCloudLayerTop,iCloudLayerBot
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,troplayer
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    REAL :: SUNCOS ! solar zenith angle cosine at surface
    REAL ::  SCOS1 ! solar zenith angle cosine at layer1
    REAL ::  VSEC1 ! satellite view zenith angle secant at layer1

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm), &
        raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'fast nlte radmodel ...'
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh=iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay)  = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
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
    IF (iNLTEStart <= kProfLayer) THEN
        write (kStdWarn,*) 'this routine expects SARTA NLTE',iNLTEStart
        CALL DoStop
    ELSE
    !        write (kStdWarn,*) 'Normal rad transfer .... '
    !        write (kStdWarn,*) 'adding on SARTA NLTE if DAYTIME'
        iLay = 1
        iSTopNormalRadTransfer = kProfLayer + 1
    END IF

    DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    iSpecular = +1    !some specular refl, plus diffuse
    iSpecular = -1    !no   specular refl, only diffuse
    IF (iSpecular > 0) THEN
        write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN_NLTE_FAST'
        CALL loadspecular(raFreq,raSpecularRefl)
        DO iFr=1,kMaxPts
        ! aSpecularRefl(iFr) = 0.0272   !!! smooth water
            raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+ &
            raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
            raSun(iFr)*(raSpecularRefl(iFr) + raSunRefl(iFr))
        END DO
    ELSE
        DO iFr=1,kMaxPts
            raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+ &
            raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
            raSun(iFr)*raSunRefl(iFr)
        END DO
    END IF

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    DO iLay=iHigh,iHigh
    ! o rad transfer to TOA (80 km)
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+ &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

        suncos = raSunAngles(iaRadLayer(1))           !! at surface
        scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
        vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

        suncos = cos(suncos*kPi/180.0)
        scos1  = cos(scos1*kPi/180.0)
        vsec1  = 1/cos(vsec1*kPi/180.0)

        IF (iDoSolar >= 0) THEN
        !         do iFr = 1,iNumlayer
        !           write(kStdWarn,*) iFr,iaRadLayer(iFr),raSunAngles(iaRadLayer(iFr))
        !           end do
            write(kStdWarn,*)'day .... add SARTA_NLTE at solangle for chunk ',raSunAngles(iaRadLayer(1)),raFreq(1)
            CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
            iaRadLayer,iNumlayer,raInten,rCO2MixRatio)
        ELSE
            write(kStdWarn,*)'nighttime ... do not need SARTA_NLTE'
        END IF

        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FAST

!************************************************************************
! allows a NLTE atmosphere, does USUAL FAST computation using COMPPRESSED NLTE tables

! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles, &
    iTag,raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer),rCO2MixRatio
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2,rOmegaSun
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    INTEGER :: iCloudLayerTop,iCloudLayerBot
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,troplayer
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    iIOUN = iIOUN_IN

    rThermalRefl = 1.0/kPi
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm), &
        raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar

    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO

! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'slow nlte radmodel ...'
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh=iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay)  = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
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
    IF (iNLTEStart > kProfLayer) THEN
        iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
        write (kStdErr,*) 'Normal rad transfer .... no NLTE'
        write (kStdErr,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
        write (kStdErr,*) 'should be calling rad_trans_SAT_LOOK_DOWN'
        CALL DoStop
    ELSE
        iLay = 1
        987 CONTINUE
        iL=iaRadLayer(iLay)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF ((iLModKprofLayer < iNLTEStart) .AND. (iLay < iNumLayer)) THEN
            iLay = iLay + 1
            GOTO 987
        END IF
        iSTopNormalRadTransfer = iLay
        write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
        write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
    END IF

    DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        IF (iLModKprofLayer == 0) THEN
            iLModKprofLayer = kProfLayer
        END IF
        IF (iLModKprofLayer < iNLTEStart) THEN
        ! ormal, no LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iLModKprofLayer >= iNLTEStart) THEN
        ! ew; LTE emission stuff
            DO iFr=1,kMaxPts
                rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
                raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    DO iFr=1,kMaxPts
        raInten(iFr)=raSurface(iFr)*raUseEmissivity(iFr)+ &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
        raSun(iFr)*raSunRefl(iFr)
    END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
            raInten(iFr)=raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    DO iLay=iHigh,iHigh
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

        IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp >= 1) THEN

                write(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer'
                write(kStdWarn,*) 'This is the top of the usual AIRS atmosphere'
                write(kStdWarn,*) '   you have iDoUpperATM > 0'
                write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
                write(kStdWarn,*) 'and output stuff into the blah_UA file'
                write(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of'
                write(kStdWarn,*) 'stratosphere into both this and the UA file'

            ! o radiative transfer thru this layer
                DO iFr=1,kMaxPts
                    raInten(iFr) = &
                    raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
                END DO

            ! ow do complete rad transfer thru upper part of atmosphere
                CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)), &
                iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of stratosphere (120-200 km)

                write(kStdWarn,*) 'finally outputting radiances at TOTAL Complete TOA into'
                write(kStdWarn,*) 'usual binary file (iLay = ',iLay,')'

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
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        END IF

    !c no need to do radiative transfer thru this layer
    !c        DO iFr=1,kMaxPts
    !c          raInten(iFr)=raaEmission(iFr,iLay)+
    !c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
    !c          END DO

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR

!************************************************************************
! this subroutine very quickly does the radiative transfer
! since the optical depths are soooooooooo small, use double precision
    SUBROUTINE UpperAtmRadTrans(raInten,raFreq,rSatAngle, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
!   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!   upper atm abs coeff and planck coeff
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!   input wavevector and integer stating which layer to stop rad transfer at
    REAL :: raFreq(kMaxPts),rSatAngle
    INTEGER :: iUpper
! do we want to do upper atm NLTE computations?
    INTEGER :: iDoUpperAtmNLTE
! do we dump all or some rads?
    INTEGER :: iDumpAllUARads
! input/output pararameter
!   this contains the radiance incident at kCARTA TOA (0.005 mb)
!   it will finally contain the radiance exiting from TOP of UPPER ATM
    REAL :: raInten(kMaxPts)

! local variables
    INTEGER :: iFr,iL,iIOUN
    REAL :: rEmission,rTrans,rMu,raInten0(kMaxPts)
    DOUBLE PRECISION :: daInten(kMaxPts),dTrans,dEmission
    CHARACTER(80) :: caOutName

    caOutName = 'DumDum'
    iIOUN = kNLTEOutUA
      
    IF (iDoUpperAtmNLTE <= 0) THEN
        write (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
        CALL DoStop
    ELSE
        write(kStdWarn,*) 'Doing UA (NLTE) radtransfer at 0.0025 cm-1 '
    END IF

! compute radiance intensity thru NEW uppermost layers of atm
    DO iFr = 1,kMaxPts
        raInten0(iFr) = raInten(iFr)
        daInten(iFr)  = dble(raInten(iFr))
    END DO

    iL = 0
    IF (kNLTEOutUAOpen > 0) THEN
        write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
    ! always dump out the 0.005 mb TOA radiance if the UA file is open
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF

    rMu = cos(rSatAngle*kPi/180.0)

    DO iL = 1,iUpper - 1

        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = (raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = (raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            (ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans

            raInten(iFr) = sngl(daInten(iFr))
        END DO

        IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
            write(kStdWarn,*) 'dumping out UA rads at iL = ',iL
        ! dump out the radiance at this HIGH pressure level
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    DO iL = iUpper,iUpper
        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = dble(raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = dble(raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            dble(ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans
            raInten(iFr) = sngl(daInten(iFr))

        END DO

        IF (kNLTEOutUAOpen > 0) THEN
        ! always dump out the 0.000025 mb TOA radiance if the UA file is open
            write(kStdWarn,*) 'dumping out 0.000025 mb UA rads iL = ',iL
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE UpperAtmRadTrans

!************************************************************************
! this is based on subroutine "othergases"
! this subroutine calls the routines to read in the k-compressed data
! iGasID tells which gas type, rFileStartFr identifies the frequency range
! have to send in ref temp, amount and profile temp,amount
    SUBROUTINE compressedNLTE(iGasID,rFileStartFr,iTag,iActualTag, &
    iProfLayer,iL,iU, &
    raPAmt,raRAmt,raPTemp,raRTemp,iErr,iDoDQ,pProf,iProfileLayers, &
    daaDQ,daaDT,daaAbsCoeff,iSplineType,i_NLTEFile_TYPE, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! iGasID     = GASID
! rFileStartFr    = current k-comp block of 25 cm-1 that is being processed
! iTag       = current k-comp block of 25 cm-1 that is being processed
! iProfLayer = number of layers in profile === kProfLayer
! iL,iU      = min/max layer number (=1,kMaxlayer)
! daaAbs     = final uncompressed abs coefficient for gas iGasID
! iErr       = errors (mainly associated with file I/O, could be associated
!              with incorrect number of layers in compresse database etc)
! raP/RAmt   = arrays containing actual/reference gas amounts
! raP/RPart  = arrays containing actual/reference gas partial pressures
! raP/RTemp  = arrays containing actual/reference gas temperatures
! daaDT      = analytic Jacobian wrt temperature
! daaDQ      = analytic Jacobian wrt amount
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! i_NLTEFile_TYPE     > 100 for LA ods, > 200 for LA planckcoeff
!                     > 300 for UA ods, > 300 for UA planckcoeff
    INTEGER :: iGasID,iErr,iProfLayer,iL,iU,iDoDQ,iTag,iActualTag
    INTEGER :: iProfileLayers,iSplineType,i_NLTEFile_TYPE
    REAL :: raPAmt(kProfLayer),raRAmt(kProfLayer),pProf(kProfLayer)
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! these are for Matlab style kCOmp Corner Weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer)
    INTEGER :: iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer)
    REAL ::    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer)
    REAL ::    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer)
    INTEGER :: iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer)
    REAL ::    raQ21(kProfLayer),raQ22(kProfLayer)

! local variables associated with uncompressing the database
    CHARACTER(120) :: caFName
    INTEGER :: iIOUN,iFileGasID,iNpts,iNLay,iKtype,iNk,iKm,iKn,iUm,iUn
    INTEGER :: iT0,iaTsort(kMaxTemp),iUAoriLA
    DOUBLE PRECISION :: dSfreq,dFStep,daToffset(kMaxTemp)
    DOUBLE PRECISION :: daaaKX(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaUX(kMaxPts,kMaxK)
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge,WhichGasPosn
    INTEGER :: iDefault,iLowerOrUpper,iJ

! iUAirLA = -1 for usual lower atm (so can do jacs), +1 for extra upper atm (so no jacs)
    IF (i_NLTEFile_TYPE < 300) THEN
        iUAoriLA = -1   !! lower atm
    ELSEIF (i_NLTEFile_TYPE >= 300) THEN
        iUAoriLA = +1   !! upper atm
    END IF

    iIOUN = kCompUnit
    CALL CompFileName(i_NLTEFile_TYPE,iGasID,rFileStartFr,iTag,iActualTag,caFName)
    CALL rdcomp(caFName,iIOUN,iFileGasID,dSfreq,dFStep,iNPts,iNLay, &
    iKtype,iNk,iKm,iKn,iUm,iUn,daToffset,iT0,iaTsort, &
    daaaKX,daaUX)

! check that the file has the data for the correct gas
    IF (iFileGasID /= iGasID) THEN
        iErr=1
        WRITE(kStdErr,1000) caFName,iFileGasID,iGasID
        1000 FORMAT('Error! file : ',/,A120,/, &
        'contains data for GasID ',I3,' not desired GasID ',I3)
        CALL DoSTOP
    END IF

! check that the data file has the right number of layers
    IF (iNLay /= kMaxLayer) THEN
        iErr=1
        WRITE(kStdWarn,1010) caFName,iNLay,kMaxLayer
        1010 FORMAT('WARNING NLTE FastComp file : ',/,A120,/, &
        'contains data for ',i3,' layers but kMaxLayer = ',I3)
    !        CALL DoSTOP
    END IF

! interpolate compressed data in temperature, to get abs coeff matrix
    IF ((kJacobian >= 0) .AND. (iUAoriLA == -1)) THEN
    !! can do jacs for usual 100 LA layers, with default CO2 profile
        CALL GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
        daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlinetype)
    ELSEIF ((kJacobian >= 0) .AND. (iUAoriLA == +1)) THEN
    !! cannot do jacs for UA, so just do this
        iLowerOrUpper = +2
        CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
        pProf,iProfileLayers,iSplineType,iLowerOrUpper)
    ELSEIF ((kJacobian < 0) .AND. (iUAoriLA == -1)) THEN
    !! no jacs, usual LA
        iLowerOrUpper = -1
        iLowerOrUpper = +3
    !        CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,
    !     $         raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
    !     $         pProf,iProfileLayers,iSplineType,iLowerOrUpper)
        CALL x2GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
        pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    ELSEIF ((kJacobian < 0) .AND. (iUAoriLA == +1)) THEN
    !! no jacs, UA
        iLowerOrUpper = +2
        CALL GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
        pProf,iProfileLayers,iSplineType,iLowerOrUpper)
    END IF

! because of the iKtype=1,2 possibility, do any necessary jacobians calcs HERE!
    IF (kJacobian >= 0) THEN
        IF (iDoDQ > 0) THEN
            IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                CALL FinalAmtDeriv(daaDQ,iKtype)
            END IF
        END IF
        IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
        (kActualJacs == 32) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
            CALL FinalTempDeriv(iKtype,daaAbsCoeff,daaDT,raPAmt)
        END IF
    END IF

! convert absorption coefficient correctly if necessary
    IF (iKtype == 2) THEN
        CALL RaisePower(daaAbsCoeff)
    END IF

! now compute
!    optical depth = actual gas amount * abs coeff
!    planck coeff  = ref    gas amount * abs coeff

    IF ((i_NLTEFile_TYPE >= 100) .AND. (i_NLTEFile_TYPE < 200)) THEN
        CALL AmtScale(daaAbsCoeff,raPAmt)
    ELSEIF ((i_NLTEFile_TYPE >= 300) .AND. (i_NLTEFile_TYPE < 400)) THEN
        CALL AmtScale(daaAbsCoeff,raPAmt)
    ELSEIF ((i_NLTEFile_TYPE >= 200) .AND. (i_NLTEFile_TYPE < 300)) THEN
        CALL AmtScale(daaAbsCoeff,raRAmt)
    ELSEIF ((i_NLTEFile_TYPE >= 400) .AND. (i_NLTEFile_TYPE < 500)) THEN
        CALL AmtScale(daaAbsCoeff,raRAmt)
    END IF

!      do iJ = 1,kProfLayer
!        print *,i_NLTEFile_TYPE,iJ,raPAmt(iJ),raRAmt(iJ),raPAmt(iJ)/raRAmt(iJ)
!      end do

! NOTE I have not added the co2 chi functions here as they *should* be in
! the kCompressed Database produced by running kCARTA with these chi fcns on

    RETURN
    end SUBROUTINE compressedNLTE

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE knonlte
