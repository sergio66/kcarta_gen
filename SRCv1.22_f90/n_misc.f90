! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

!! use must precede IMPLICIT and INCLUDE statements

MODULE n_misc

USE ieee_arithmetic
USE basic_common
USE freqfile
USE spline_and_sort_and_common
USE s_misc
USE n_rad_jac_scat

IMPLICIT NONE

CONTAINS

!************************************************************************
     logical function is_badnum(x)

     IMPLICIT NONE

     include '../INCLUDE/TempF90/kcartaparam.f90'

     real :: x
     logical :: test

     test = .false.

     if (x > huge(x)) then
       write(kStdWarn,*) 'x is huge ',x
       write(kStdErr,*)  'x is huge ',x
       test = .true.
     end if
     if (abs(x) > huge(x)) then
       write(kStdWarn,*) '+/-x is huge ',x
       write(kStdErr,*)  '+/-x is huge ',x
       test = .true.
     end if
       
    is_badnum = test
    end FUNCTION is_badnum

!************************************************************************
     logical function is_goodnum(x)

     logical :: test1,test2,test
     real :: x

     test2 = ieee_is_nan(x)
     test1 = is_badnum(x)

     test = .true.
     if ((test1 .EQV. .true.) .OR. (test2 .EQV. .true.)) test = .false.

     is_goodnum = test
  
     end FUNCTION is_goodnum

!************************************************************************
! set the default parameter values, for those that are not set in *PARAM
! read the parameter file to set parameters that have optional values
    SUBROUTINE SetDefaultParams

! NOTE !!!! also double check subroutine EXTRAPAR in strings2.f
! NOTE !!!! also double check subroutine SETDEFAULTPARAMS in misc.f

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iI,iJ

! acos(3/5) * 180/kPi
    kThermalAngle = 53.1313010235415598
          
! set default values here
! kLayer2Sp
!     -2     Layer transmittance            t(i)=exp(-k(i))
!     -1     Layer Optical depth            k(i)
!      1     Layer-to-Space Optical Depth   k2s(i)=sum(j=i,n)(k(j))
!      2     Layer-to-Space transmittance   t2s(i)=sum(j=i,n)exp(-k(j))
!      3     Layer-to-Ground Optical Depth  k2s(i)=sum(j=1,i)(k(j))
!      4     Layer-to-Ground transmittance  t2s(i)=sum(j=1,i)exp(-k(j))
    kLayer2Sp    = -1

! kCKD      == -1,00,21,23,24 sets which continuum version calculation to do
!              -1 : no continuum
! ----------------------------------------------------------------------------
! AER-CKD      01 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD version of Dec 2002
!       AIRS   02 : version 02    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Aug 2003
!       AIRS   03 : version 03    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Jan 2004
!       AIRS   04 : version 04    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!       AIRS   05 : version 05    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     (so far it looks like CKD v4)
!                                     On top of that it puts Dave Tobin's cs,cf
!                                     between 1300-1800 cm-1
!       AIRS   06 : version 06    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     extremely similar to CKD 4, except have
!                                     fixed the problem at 850 cm-1 which is
!                                     really HNO3 and not a continuum problem
! In asl:/carrot/s1/strow/Tobin/Tobin_radish/NISTdata2/New/, use
! cs0_tsl0_hr.mat and cf0_tsl0_hr.mat.  makecons.m with flag=7
! looks like it should load things in correctly.
! ----------------------------------------------------------------------------
! AER-STD      25 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Sep 2011
! ----------------------------------------------------------------------------
! AER-STD      27 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Feb 2016
! AER-STD      32 : self, foreign   is by CKD and modified by Mlawer/Alvarado
!                                 ... this is the MT_CKD 3.2 version of Feb 2017
! ----------------------------------------------------------------------------
! old versions of AER CKD
! AER-CKD      00 : version 00
! AER-CKD      21 : version 2.1
! AER-CKD      23 : version 2.3
! AER-CKD      24 : version 2.4
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! these were our interim versions (made by us, basically MT_CKD 1)
! no longer in use!

!    RAL       12 : self = version 2.4, foreign = blend of 0.0 and 2.4
!                                       from 0-1575, 1625-3000, use v2.4
!                                       from 1575-1625 linearly blend v2.4,0.0
!                                       using a triangle
!    RAL       13 : self = version 2.3, foreign = dave tobin's v2.3
!    RAL       90 : self = version 2.4, foreign = blend of 2.4,Tobin's thesis
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use Tobins thesis
!    RAL       50 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst50.m used linear tempr interpolation
!    RAL       51 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst51.m used power tempr interpolation
!                                       Need to be careful ... I put in Dave
!                                       Tobin's fixes for CS in the 600-1100
!                                       cm-1 range; see mst51_tobin.m in
!                                       the SPECTRA directory (~0.9*ckd24)
!    RAL       52 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst52.m used power tempr interpolation
!                                       and says CF(296)=CF(243); so use the
!                                       foreign broadened 243 data to get CS243
!    RAL       55 : self, foreign       same as 51 above, but uses
!                                a) CS fudge factor of 0.87 for v <= 1137 cm-1
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       56 : self, foreign       same as 51 above, but uses
!                                a) John Taylor fudge between 600-1200 cm-1
!                                   instead of Dave Tobin fudge
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                   these are better than the above fudges
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       60 : self, foreign     is a hybrid of 51,55 done by scott
!                                     so that bias errors are reduced, as
!                                     a function of water column amount.
!                                     eventually will include tuning coeffs
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25 32]  ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

    kCKD = 25

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
    kGasTemp = -1

! kLongOrShort == whether the user wants to save ALL header info (+1) or
!                 just a portion (-1) or output just the basic results (0)
    kLongOrShort = 1

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs
    kActualJacs = -1    !!! default

! kActualJacsT = -1, kActualJacsB = -1
! if we set kActualJacs = 100, then we can also set the layers of the column
! that we want to perturb; default = -1/-1 means all layers (kActualJacs = 100)
! else kActualJacsABCXYZ sets these as kActualJacsT=ABC, kActualJacsB=XYZ
    kActualJacsT = -1
    kActualJacsB = -1

! kJacobOutput == -1 if we output d(radiance)/dq,d(radiance)/dT
!                  0 if we output d(radiance)/dq * q, d(radiance)/dT
!                  1 if we output d(BT)/dq * q, d(BT)/dT
    kJacobOutput = 1

! kFlux == -1 if we do not want flux/OLR computations or output NLTE
!                                                        Planck modifiers
!       ==  1 if we want DNWELL flux at every layer      units = mW m-2
!       ==  2 if we want flux computations : output      units = kelvin day-1
!       ==  3 if we want UPWELL flux at every layer      units = mW m-2
!       ==  4 if we want outgoing OLR only at TOA,       units = mW m-2
!       ==  5 if we want outgoing OLR at TOA, ILR at GND units = mW m-2
!       ==  6 if we want DNWELL and UPWELL flux at each layer units = mW m-2
!   --> ==  0 if we want to output NLTE Planck modifiers
    kFlux = -1

! kPlanckOut == -1 if we do not want to output Planck modifiers, 1 if we do
    kPlanckOut = -1

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
! kSurfTemp = -1 == want to use user supplied surface temp in *RADNCE
!              1 == want to use user supplied surface temp in *RADNCE as
!                     an offset to pressure interpolated temperature
! so in above RTP cases if kSurfTemp < 0 : use raTSurf(iI)
!                          kSurfTemp > 0 : use raTSurf(iI) + InterpedTemp
    kSurfTemp = -1

! kRTP = -5  : read LBLRTM style LAYERS profile (edited TAPE 6); set atm from namelist
! kRTP = -6  : read LBLRTM style LEVELS profile (       TAPE 5); set atm from namelist
! kRTP = -10 : read TEXT style LEVELS   profile; set atm from namelist
! kRTP = -2  : read GENLN4 style LAYERS profile; set atm from namelist
! kRTP = -1  : read old style kLAYERS   profile; set atm from namelist
! kRTP =  0  : read RTP style kLAYERS   profile; set atm from namelist
! kRTP = +1  : read RTP style kLAYERS   profile; set atm from RTP file
! kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
    kRTP = 1

! kTempJac == -2 if we only want d/dT(planck) in temperature jacobian
!             -1 if we only want d/dT((1-tau)(tau->sp)) in temp jacobian
!              0 if we want complete d/dT(planck) + d/dT(tau) in temp jac
    kTempJac = 0

! the following cannot be controlled by the user using *PARAMS
! all the radiance parameters and the Jacobian parameter
! kJacobian ==  1 if analytic Jacobians are to be computed
!              -1 if we do the standard Genln2 computations w/o Jacobians
    kSolar        = 1     !turn on solar
    kSolarAngle   = 0.0   !solar angle
    kSolarRefl    = -1.0  !use (1-ems)/pi
    kThermal      = 0     !use fast diffusive approx
    kThermalAngle = -1.0  !use acos(3/5) in upper layers
    kThermalJacob = 1     !use thermal backgnd in Jacobians

    kJacobian = -1         !do not do Jacobians
    kWhichScatterCode = -1 !do not do scattering computations
    kScatter  = -1         !do not do scattering computations

    k100layerCloud = -1   !assume rtp file does NOT have 100 layer cloud

! 2016
! allow nm_params to define defaults
!   GENERAL      iaaDefaults(1,:) = iSARTAChi   iSPlineType iCO2Chi  iMatlabORf77   iWhichScatterCode iMethod
!   RT/BackTherm iaaDefaults(2,:) = iGaussPts iGaussQuad  iSnell      iInterpType  iWhichRT  (kTemperVary set in mn_radnce)
!   Misc         iaaDefaults(3,:) = mix table --> klayers, dump debug info to text etc
! not done
!   NLTE         iaaDefaults(4,:) = iCurrent    iTalk       iTestGenln2  iNoPressureShiftCO2 iQtips_H98
!                             iLinearOrSpline iDoCO2Continuum iMethod
!   TAPE5/6     iaaDefaults(4,:) = iReplaceZeroProf iAIRS101_or_LBL_levels IPLEV iAddLBLRTM
!      INTEGER iaaOverrideDefault(4,10)
!      COMMON/comBlockDefault/iaaOverrideDefault
    iaaOverrideDefault = -9999
          
! GENERAL
    caaTextOverrideDefault  = 'notset'
          
    iaaOverrideDefault(1,1) = -1    !!! iSARTAChi = -1  for no tuning, see kcartabasic/kcartamain/kcartajpl
                                    !!!                 kcartaparallel and finally used in kcoeffMAIN.f
    iaaOverrideDefault(1,2) = +1    !!! iSplinetype = +1 for SUBR iSetSplineType in kcartamisc.f
    iaaOverrideDefault(1,3) = +2    !!! iCO2Chi = +2     for SUBR multiply_co2_chi_functions in kcoeffMAIN.f
    iaaOverrideDefault(1,4) = +1    !!! iMatlabORf77 = +1  use Maltab style uncompression,  kcoeffMAIN.f
    iaaOverrideDefault(1,5) = +5    !!! iWHichScatterCode = 5 for PCLSAM in rtp_interface.f
    iaaOverrideDefault(1,6) = +1    !!! iReadP = 1 when assuming GENLN2 style profile in n_pth_mix.f
    iaaOverrideDefault(1,7) = -1    !!! iLogOrLinear = -1 when interp scat tables SUBR INTERP_SCAT_TABLE2
                                    !!!   in clear_scatter_misc.f
    iaaOverrideDefault(1,8) = -1    !!! -1 to keep h.vcmin/h.vcmax as read in from RTPfile, +1 to override with rf1,rf2
    iaaOverrideDefault(1,9) =  0    !!! iCO2WVCA = +2 will turn it on for SUBR add_co2_wv_n2_continuum in kcoeffMAIN.f/kcoeff_common.f90
! new for cloudy jacs
    iaaOverrideDefault(1,10) =  -1  !!! output only final jacobian (ie wieghted sum); if you set +1, it will output all jacobians
                                    !!!   so the file can be 5 times larger than needed (ice/water/both/clear/total) HUGE! SAD!

! RadTrans
! see n_main.f90
!    kTemperVary  = -1          !assume const-in-tau temperature variation  till Oct 2019
!    kTemperVary  = +43         !assume linear-in-tau temperature variation after Oct 2019
    iaaOverrideDefault(2,1) = kTemperVary !!! kTemperVary .... can be reset in nm_radnce, and then subr SetkTemperVary
    iaaOverrideDefault(2,2) = +3    !!! THIS IS LBLRTM STYLE iGaussPts = 3 for flux and downwell gauss quad
!!!   see SUBR IntegrateOverAngles_LinearInTau in rad_quad.f
    iaaOverrideDefault(2,3) = 0     !!! SUBR BackGndThermal in rad_diff.f
!!! iDothermal = kThermal; if iDoThermal = -1, no backgnd thermal computed
!!!                                      =  0, backgndthermal with diffusive approx << DEFAULT >>
!!!                                            --->>> control further with iaaOverrideDefault(2,4) <<<---
!!!                                     = +10 backgndthermal with diffusive approx [NEW]
!!!                                       constant angle (typicall acos(3/5)) in all layers, but refl th = (1-e) [Nick Nalli] NEW
!!!                                            --->>> control further with iaaOverrideDefault(2,4) <<<---
!!!                                      = +1, use integration over angles, const-in-tau  layer T
!!!                                            --->>> control further with iaaOverrideDefault(2,5) <<<---
!!!                                      = +2, use integration over angles, linear-in-tau layer T
!!!   this is the main routine, called by all downwelling RT routines in rad_main.
!!!   all of them have -1 for iDoAcos35
!!!     calls IntegrateOverAngles_LinearInTau (iDoThermal = 2)
!!!     calls IntegrateOverAngles             (iDoThermal = 1) -- also can set iaaOverrideDefault(2,5)
!!!     calls DoDiffusivityApprox             (iDoThermal = 0) << DEFAULT >>
    iaaOverrideDefault(2,4) = -1    !!! SUBR radnce4RTP in n_rtp.f90
!!!   raKThermalAngle(iC) = iaaOverrideDefault(2,4) in rtp_interface.f which then becomes kThermalAngle
!!!     = -1, fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 (acos(3/5))
!!!
!!!   this sets  kSetThermalAngle = -1 for acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!!!                               = +1 for constant angle (typically acos(3/5)) in all layers
!!!                               = -2     same as -1, except linear-in-tau T variation
!!!                               = +2     LBLRTM style 3 exponential gauss quad, not yet implemented
!!! SUBR DoDiffusivityApprox in rad_diff.f uses this info
!!!   iDiffMethod = kSetThermalAngle
!!!     = -1 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!!!          differs from iaaOverrideDefault(2,5) = 0 since here, upper layers use acos(3/5) lower layers are accurate
!!!                                                   while there, use layer-varying accurate acos(rDiffusive)
!!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 = acos(3/5) << DEFAULT >>
!!!           this can be controlled by kThermalAngle, either in nm_params for kRTP = +1
!!!                                                   or rakThermalAngle() for kRTP = 0,-1
!!!           so in nm_params : set iaaOverride(2,4) = 1, kThermalAngle = 50.0 and that works!!!
!!!     = -2 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers, linear in tau T
!!!     = +2 diffusive background using LBLRTM style 3 exponetial gauss quad, not yet implemented
    iaaOverrideDefault(2,5) = 0     !!! SUBR IntegrateOverAngles in rad_quad.f, called by SUBR BackGndThermal
!!!   iGaussQuad =    -1 for integrate using newton quad 0:90/20:90 (VERY SLOW)
!!!                    0 for accurate diffusivity                   (AT ALL LAYERS << DEFAULT >>)
!!!                      so this differs from iaaOverrideDefault(2,3) = 0,iaaOverrideDefault(2,4) = -1
!!!                      where acos(3/5) is used in upper ayers, and accurate diffusive angle in lower layers
!!!                   +1 for gausslegendre w(i) at theta(i)         (QUITE SLOW)
    iaaOverrideDefault(2,6) = +1    
!!! iUsualUpwell = +1 for upwell RT with surface term, << DEFAULT >>
!!!                -1 with no surface,
!!!                -2 to only dump downwelling backgnd
!!!   see SUBR find_radiances in rad_main.f
    iaaOverrideDefault(2,7) = -1    
!!! iUseSnell = -1 for No  Snell law raytrace YES layer curvature effects, similar to SARTA <<<default >>>
!!!           = +1 for Yes Snell law raytrace YES layer curvature effects
!!!           = 0  for No  Snell law raytrace NO  layer curvature effects  PLANE PARALLEL
!!!   see SUBR FindLayerAngles in rad_angles.f
    iaaOverrideDefault(2,8) = +1    !!! iInterpType = +1 to turn (pav,Tav) into (plevs,Tlevs), only used if kTemperVary = 43

!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022
!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022
!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022
!!! >>>> --- >>>    see SUBR Get_Temp_Plevs in n_pth_mix.f
!!! >>>> --- >>>    iaaOverrideDefault(2,9) = -1    !!! >>>> --->>> iLBLRTM_highres = -1 do not estimate/fix problems because use 0.0025 cm-1, when kTemperVary = 43 << DEFAULT>>
!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022
!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022
!!! >>>> --- >>>    use of iaaOerrideDefault(2,9) for HighRes LayerInTemp adjustments has been retired July 2022

! see scatter_pclsam_code.f90, this is the adjust chou routine
    iaaOverrideDefault(2,9) = 35    !!!! this is now used for ChouAdjust * 100
    iaaOverrideDefault(2,9) = -3    !!!! this is now used for ChouAdjust = parametrized routine

!!!   see SUBR rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMISS in rad_main.f
!!!   0 for ABS clouds, 2 for RTPSEC, 3 for DISORT
!!!    iaaOverrideDefault(2,10) = 5    !!! kWhichScatterCode =  5 for PCLSAM (Default)                 TILL  MAY 2021
!!!    iaaOverrideDefault(2,10) = 5    !!! kWhichScatterCode = -1/+1 for adjusting Planet/Sun distance AFTER MAY 2021
    iaaOverrideDefault(2,10) = -1   !!! if kPlanet == 3, then adjust Planet/Sun dist (-1 no, default, +1 yes)
      
! n_layers_lblrtm.f and n_pth_mix.f  TAPE5/6 
    iaaOverrideDefault(3,1) = -1    !!! iAIRS101_or_LBL_levels use LBLRTM, not AIRS 101 levels, for integration
    iaaOverrideDefault(3,2) = +1    !!! iReplaceZeroProf = +1 to add in profiles TAPE5 does not have
    iaaOverrideDefault(3,3) = -1    !!! iAddLBLRTM = -1 when gas profile missing from TAPE5/6, do not add it in
! plus EXTRA THINGS
    iaaOverrideDefault(3,4) = -1    !!! dump out surface terms (downwell therm, emiss, rho) : -1 = NO (default), +1 = yes
    !iaaOverrideDefault(3,5) = +1   !!! if there are clouds in rtp, keep them (TwoSlab); if this is -1, switch to clear sky
    iaaOverrideDefault(3,5) = 0     !!! 0 means do scatter/no scatter as set in rtp (default), works with iForceScatterCalc_EvenIfNoCld
                                    !!! +1 = force scatter (PCLSAM/DISORT), -1 = force clear                         
    iaaOverrideDefault(3,6) = -9999 !!! if -9999 ignore this switch
                                    !!! else switch between -1 (no sun) 0 (ttorad(v,SunTemp) and +1 (tables)
                                    !!!   remember solar tables are at 0.0025 cm-1 so will have problems for 605-905 cm-1
                                    !!!   if iaKSolar = +1, so switch to iaKSolar = 0
    iaaOverrideDefault(3,7) = +1    !!! SARTA scat tables are delta scaled, leave them scaled <DEFAULT>, -1 is to un delta scale
    iaaOverrideDefault(3,8) = -1    !!! the pressure levels in "klayers.op.rtp,klayers.rp.rtp" and "../INCLUDE/TempF90/KCARTA_databaseparam.f90" DISagree
                                    !!! default = -1
                                    !!! if set to -1, then for example compiled assuming 
                                    !!!            Earth AIRS STD  100 and rtp from /home/chepplew/gitLib/klayersV205/Src_rtpV201/oco2
                                    !!! if set to +1, then for example compiled assuming 
                                    !!!   see eg Makefile_tar_objs_data_f90_datafix for 
                                    !!!            Earth AIRS STD  100 and rtp from /asl/packages/klayersV205/BinV201/klayers_airs_wetwater
                                    !!!            Earth AIRS PBL  100 and rtp from /home/chepplew/gitLib/klayersV205/Src_rtpV201/oco2
                                    !!!            Mars  CIRAS STD 100 and rtp from /home/sergio//SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/klayersV205_Mars


    RETURN
    end SUBROUTINE SetDefaultParams

!************************************************************************
    SUBROUTINE Check_iaaOverrideDefault

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iI,iJ,iTemp

!   iaaOverrideDefault(1:4,1:10)

    !************************************************************************    
    iI = 1
    iJ = 0

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,1)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1) .AND. (iTemp /= 2) .AND. (itemp /= 3)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
       'iSartaChi : need iaaOverrideDefault(',iI,',',iJ,') = -1,+1,2,3 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,2)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1) .AND. (abs(iTemp) /= 2)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
       'iSplineType : need iaaOverrideDefault(',iI,',',iJ,') = +/- 1,2 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,3)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (iTemp < 0 .OR. iTemp > 2) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
        'iCO2Chi : need iaaOverrideDefault(',iI,',',iJ,') = 0,1,2 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,4)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) .NE. 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
        'iMatlabORf77 : need iaaOverrideDefault(',iI,',',iJ,') = -1,+1 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,5)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (iTemp < 0 .OR.  iTemp > 6) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
        'iWhichScatterCode : need iaaOverrideDefault(',iI,',',iJ,') = 0..6 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,6)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
       'iReadRP : need iaaOverrideDefault(',iI,',',iJ,') = +/- 1 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,7)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
       'iLogOrLinear : need iaaOverrideDefault(',iI,',',iJ,') = +/- 1 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,8)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
        'rfmin/rfmax : need iaaOverrideDefault(',iI,',',iJ,') = +/- 1 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,9)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (iTemp .NE. 0 .and. iTemp .NE. 2 .and. iTemp .NE. 4 .and. iTemp .NE. 6) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') &
        'CO2 N2 WV CIA : need iaaOverrideDefault(',iI,',',iJ,') = 0,2,4,6,8 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(1,10)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') &
        'size of cloudy jac file : all or just weighted sum  : need iaaOverrideDefault(',iI,',',iJ,') = +/- 1 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    !************************************************************************    
    iI = 2
    iJ = 0

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,1)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1 .AND. iTemp /= 2 .AND. iTemp /= 42 .AND. iTemp /= 43) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'kTemperVary : need iaaOverrideDefault(',iI,',',iJ,') = -1,1,2,42,43 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,2)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (iTemp /= 1 .AND. iTemp /= 2 .AND. iTemp /= 3 .AND. iTemp /= 4) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iNumGaussPts : need iaaOverrideDefault(',iI,',',iJ,') = 1..4 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,3)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((iTemp < -1 .AND. iTemp > 2) .AND. (iTemp /= 10)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iDoThermalBackGnd : need iaaOverrideDefault(',iI,',',iJ,') = -1,0,1,2,10 not ', &
         iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,4)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1) .AND. (abs(iTemp) /= 2)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iraKThermalAngle : need iaaOverrideDefault(',iI,',',iJ,') = +/-1,2 not ', &
      iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,5)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) > 0) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iGaussQuad : need iaaOverrideDefault(',iI,',',iJ,') = -1,0,+1 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,6)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1) .AND. (abs(iTemp) /= 2)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iUpwell : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 , -2 and not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,7)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) > 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iSnell : need iaaOverrideDefault(',iI,',',iJ,') = -1,0,+1 not ',  &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,8)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iInterpType : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

!!! >>>> --- >>> 
!!! >>>> --- >>>     iJ = iJ+1
!!! >>>> --- >>>     iTemp = iaaOverrideDefault(2,9)
!!! >>>> --- >>>     iTemp = iaaOverrideDefault(iI,iJ)
!!! >>>> --- >>>     IF (abs(iTemp) /= 1) THEN
!!! >>>> --- >>>       write(kStdErr,'(A,I2,A,I2,A,I2)') 'iEstimateHighRes : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ', &
!!! >>>> --- >>>         iaaOverrideDefault(iI,iJ)
!!! >>>> --- >>>       CALL DoStop
!!! >>>> --- >>>     END IF
!!! >>>> --- >>>     IF (iTemp > 0) THEN
!!! >>>> --- >>>       write(kStdErr,'(A)') 'Plan to retire iEstimateHighRes iaaOverrideDefault(2,9)'
!!! >>>> --- >>>       CALL DoStop
!!! >>>> --- >>>     END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,9)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (iTemp .LT. -3 .OR. iTemp .GT. 100) THEN
      write(kStdErr,'(A,I2,A,I2,A,I8)') & 
        'iAdjustChou : need iaaOverrideDefault(',iI,',',iJ,') between 0 and 100, or -1,-2,-3, not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

!    iJ = iJ+1
!    iTemp = iaaOverrideDefault(2,10)
!    iTemp = iaaOverrideDefault(iI,iJ)
!odd looks like iaaOverrideDefault(1,5)
!    IF (iTemp < 0 .OR. iTemp > 5) THEN
!      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iWhichScatterCode : need iaaOverrideDefault(',iI,',',iJ,') = 0..5 not ',iaaOverrideDefault(iI,iJ)
!      CALL DoStop
!    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(2,10)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iAdjust Sun-Planet Distance : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ', &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    !************************************************************************    
    iI = 3
    iJ = 0

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,1)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iAIRS101_or_LBL_levels : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ',  &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,2)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iReplaceXeroProf for LBLRTM : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ',  &
      iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,3)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') 'iAddLBLRTM missing profile : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ',  &
      iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,4)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF (abs(iTemp) /= 1) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') & 
        'Dump out surface terms (upwell radiation) : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ',  &
        iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,5)
    iTemp = iaaOverrideDefault(iI,iJ)
    !IF (abs(iTemp) /= 1) THEN
    !  write(kStdErr,'(A,I2,A,I2,A,I2)') & 
    !    'Keep rtp cloud info (+1) or switch to clear only (-1)  : need iaaOverrideDefault(',iI,',',iJ,') = +/-1 not ', &
    !    iaaOverrideDefault(iI,iJ)
    !  CALL DoStop
    !END IF
    IF ((abs(iTemp) > 1)) THEN
      write(kStdErr,'(A,I2)') &
        'Can figure clear or scattering as in rtp file (0, default), force scattering calcs even if clear (+1), & 
         force clear calc even if clouds (-1) .. not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,6)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1) .AND. (iTemp .NE. 0) .AND. (iTemp .NE. -9999)) THEN
      write(kStdErr,'(A,I2,A,I2,A,I2)') &
        'Keep sun as -9999 or -1,0,+1 : need iaaOverrideDefault(',iI,',',iJ,') = +/-1,0,-9999 not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,7)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1)) THEN
      write(kStdErr,'(A,I2)') &
        'Keep scattering tables DeltaScaled (+1, default) or unscale them (-1) .. not ',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    iJ = iJ+1
    iTemp = iaaOverrideDefault(3,8)
    iTemp = iaaOverrideDefault(iI,iJ)
    IF ((abs(iTemp) /= 1)) THEN
      write(kStdErr,'(A,I2)') &
        'use ../INCLUDE/TempF90/KCARTA_databaseparam.f90 : PLEV_KCARTADATABASE_AIRS consistent & 
        with klayers grid (+1/default) or not (-1)',iaaOverrideDefault(iI,iJ)
      CALL DoStop
    END IF

    RETURN
    end SUBROUTINE Check_iaaOverrideDefault

!************************************************************************
! now check parameters in *PARAM
! this subroutine checks parameters in *PARAMS ... abort if they
! do not make sense ..
    SUBROUTINE CheckParams

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: i0,iT,iB,iJ,iGah,iConstOrVary
    CHARACTER(9) :: iIOUN9

    write(kStdWarn,*) 'checking parameters (from *PARAMS) .... '
    write(kStdWarn,*) '  '

    IF ((iabs(kLayer2Sp) /= 1) .AND. (iabs(kLayer2Sp) /= 2) .AND. (kLayer2Sp /= 3) .AND. (kLayer2Sp /= 4)) THEN
      write(kStdErr,*) 'In *PARAMS, need kLayer2Sp = +/-1,+/-2,+3,+4'
      write(kStdErr,*) 'kLayer2Sp == do layer-to-space calc or not'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF
    IF (iabs(kGasTemp) /= 1) THEN
      write(kStdErr,*) 'In *PARAMS, program needs kGasTemp = +/- 1'
      write(kStdErr,*) 'kGasTemp = use CO2 temperature profile or not'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

! CKD    releases are 0,21,23,24
! MT_CKD releases are 1,25
! our modifications are 2,3,4,5 (derived from 1)
!                       51,55,60 (derived from analysing RAL data)
! and various 12,13,50,52,56 which might be gotten rid of eventually
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25] [] ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

    IF ((kCKD /= -1) &
!!! the std CKD pre-2002 versions
     .AND. (kCKD /= 0) .AND. (kCKD /= 21) .AND. (kCKD /= 23) .AND. (kCKD /= 24) &
!!! these are MT_CKD1 and research versions from AIRS data
     .AND. (kCKD /= 1) .AND. (kCKD /= 4) .AND. (kCKD /= 6) &
!!! these are MT_CKD25 owards more recent versions
     .AND. (kCKD /= 25) .AND. (kCKD /= 27) .AND. (kCKD /= 32)) &
    THEN
      write(kStdErr,*) 'In *PARAMS, need kCKD = [-1] for no continuum OR'
      write(kStdErr,*) '                 CKD    versions 0,21,23 or 24'
      write(kStdErr,*) '              MT_CKD    versions 1,  [4,6]'
      write(kStdErr,*) '              MT_CKD    versions 25  [   ]'
      write(kStdErr,*) '              MT_CKD    versions 27  [   ]'
      write(kStdErr,*) '       (latest AER versions =  1, released Dec 2002)'
      write(kStdErr,*) '       (latest AER versions = 25, released Dec 2010)'
      write(kStdErr,*) '       (latest AER versions = 27, released Feb 2016)'
      write(kStdErr,*) '       (latest AER versions = 32, released Feb 2017)'	
      write(kStdErr,*) '           [ are UMBC modifications of CKD or MT-CKD] '
      write(kStdErr,*) 'kCKD is water continuum calculation version'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

    IF (iabs(kLongOrShort) > 2) THEN
      write(kStdErr,*) 'In *PARAMS, program needs kLongOrShort = -2,-1,0,+1,+2'
      write(kStdErr,*) 'kLongOrShort = complete header info (+1) or not (-1),  long warning.msg'
      write(kStdErr,*) 'kLongOrShort = complete header info (+2) or not (-2), short warning.msg'
      write(kStdErr,*) '               or file containing results only (0)'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs

    IF ((kActualJacs /= -1)  .AND. (kActualJacs /= 20)   .AND. &
      (kActualJacs /= 30)  .AND. (kActualJacs /= 40)   .AND. &
      (kActualJacs /= 50)  .AND. (kActualJacs /= 100)  .AND. &
      (kActualJacs /= -2)  .AND. (kActualJacs /= 32)   .AND. &
      (kActualJacs /= 102) .AND. (kActualJacs < 102)) THEN
      write(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
      write(kStdErr,*) '  or 100 or 100ABCXYZ'
      write(kStdErr,*) 'OR -2, 32,102 or 102ABCXYZ'
      write(kStdErr,*) 'kActualJacs = 100(2) ==> column gas/temp jacs '
      write(kStdErr,*) '   for all layers, plus stemp jac'
      write(kStdErr,*) 'kActualJacs = 100(2)ABCXYZ ==> column gas/temp jacs '
      write(kStdErr,*) '   for layers ABC to XYZ, plus stemp jac'
      write(kStdErr,*) 'kActualJacs = actually compute all profile jacs (-1)'
      write(kStdErr,*) '              or Q(z),T(z),W(z) jacs (0s elsewhere)'
      write(kStdErr,*) ' '
      write(kStdErr,*) 'You have set this as ',kActualJacs
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF
    iT = -1
    iB = -1
    iGah = -1

    IF (kActualJacs > 102) THEN
      IF (kActualJacs < 100000001) THEN
        write(kStdErr,*) 'too few characters in kActualJacs : ',kActualJacs
        write(kStdErr,*) 'check it !'
        CALL DoStop
      END IF
      IF (kActualJacs > 102999999) THEN
        write(kStdErr,*) 'too many characters in kActualJacs : ',kActualJacs
        write(kStdErr,*) 'max possible value is 102999999'
        CALL DoStop
      END IF

      i0 = kActualJacs
      !! i0 better be 9 characters long
      write(iIOUN9,99) kActualJacs
      iJ = 9
 10   CONTINUE
      IF ((iIOUN9(iJ:iJ) /= ' ') .AND. (iJ >= 1)) THEN
        !write(kStdErr,*) 'good',iJ,iIOUN9(iJ:iJ),kActualJacs
        iJ = iJ - 1
        GOTO 10
      ELSEIF (iJ > 0) THEN
        iGah = iJ
            write(kStdErr,*) 'space in  kActualJacs at ',iJ,iIOUN9(iJ:iJ)
      END IF
      IF (iGah > 0) THEN
        write(kStdErr,*) 9-iGah,' chars in kActualJacs = ',kActualJacs
        write(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
        write(kStdErr,*) '  or 100 or 100ABCXYZ'
        write(kStdErr,*) '  or 102 or 102ABCXYZ'
        write(kStdErr,*) 'need 9 characters 10E ABC XYZ .. try again'
        CALL DoStop
      END IF

      write(kStdWarn,*) 'kActualJacs passed test ... '
      IF (kActualJacs <= 102000000) THEN
        iT = i0 - 100000000
        iT = int(iT/1000)
        iB = i0 - int(i0/1000)*1000
        kActualJacs = 100
      ELSE
        i0 = i0 - 2000000
        iT = i0 - 100000000
        iT = int(iT/1000)
        iB = i0 - int(i0/1000)*1000
        kActualJacs = 102
      END IF

      IF (iT < iB) THEN
        iJ = iT
        iT = iB
        iB = iJ
      END IF
      IF (iT > kProfLayer) THEN
        write(kStdWarn,*) 'IT = ',iT,' greater than kProfLayer = ',kProfLayer
        write(kStdWarn,*) 'resetting iT = kProfLayer'
        iT = kProfLayer
      END IF
      IF (iB > kProfLayer) THEN
        write(kStdWarn,*) 'IB = ',iB,' greater than kProfLayer = ',kProfLayer
        write(kStdWarn,*) 'resetting iB = 1'
        iB = 1
      END IF
      kActualJacsT = iT
      kActualJacsB = iB
    END IF
 99 FORMAT(I9)

    IF ((iabs(kJacobOutput) /= 1) .AND. (kJacobOutput /= 0) .AND. &
      (kJacobOutput /= 2))  THEN
      write(kStdErr,*) 'kJacobOutput = ',kJacobOutput
      write(kStdErr,*) 'In *PARAMS, need kJacobOutput =-1,0,+1,+2'
      write(kStdErr,*) 'kJacobOutput = format to output Jacobians'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

    IF ((kFlux < -1) .OR. (kFlux >  7)) THEN
      write(kStdErr,*) 'In *PARAMS, program needs kFlux = -1 OR 1,2,3,4,5,6,7'
      write(kStdErr,*) 'where kFlux = do/do not compute fluxes'
      write(kStdErr,*) 'OR         program needs kFlux = -1,+1,2,3,4,5,6,7 OR 0'
      write(kStdErr,*) 'where kFlux = do not/do  output NLTE Planck'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
    IF ((abs(kSurfTemp)-1.0) >= 1e-5) THEN
      write(kStdErr,*) 'In *PARAMS, program needs kSurfTemp = +/-1.0'
      write(kStdErr,*) 'where kSurfTemp tells the program how to use'
      write(kStdErr,*) 'the surftemperatures in *RADNCE'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

!!!kRTP = -6 : read LBLRTM       LAYERS profile; set atm from namelist
!!!kRTP = -5 : read LBLRTM       LEVELS profile; set atm from namelist
!!!kRTP = -10 : read LEVELS            profile; set atm from namelist
!!!kRTP = -2  : read GENLN2 style LAYERS profile; set atm from namelist
!!!kRTP = -1  : read old    style LAYERS profile; set atm from namelist
!!!kRTP =  0  : read RTP   style kLAYERS profile; set atm from namelist
!!!kRTP = +1  : read RTP   style kLAYERS profile; set atm from RTP file
!!!kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
    IF ((kRTP /= -5) .AND. (kRTP /= -6) .AND. (kRTP /= +2) .AND. &
    (kRTP /= -10) .AND. ((kRTP < -2) .OR. (kRTP > 1))) THEN
      write(kStdErr,*) 'Need to set RTP = -10,-6,-5,-2,-1,0,+1,+2'
      write(kStdErr,*) 'Please reset kRTP and retry'
      CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == 1)) THEN
      write(kStdErr,*) 'Cannot read surface temperature info from RTP file'
      write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
      write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
      CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == 2)) THEN
      write(kStdErr,*) 'Cannot read surface temperature info from JPL/NOAA input'
      write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
      write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
      CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. ((kRTP == -5) .OR. (kRTP == -6))) THEN
      write(kStdErr,*) 'Will read surface temperature info from LBLRTM file'
      write(kStdErr,*) 'and ask kCARTA to add on raTSurf offset from nm_radnces!!!'
      !        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
      !        CALL DoSTOP
    END IF

    IF ((kSurfTemp > 0) .AND. (kRTP == -10)) THEN
      write(kStdErr,*) 'Cannot read surface temperature info from LEVELS TXT file'
      write(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
      write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
      CALL DoSTOP
    END IF
     
    IF ((kTempJac < -2) .OR. (kTempJac > 0)) THEN
      write(kStdErr,*) 'In *PARAMS, program needs kTempJac=-2,-1,0'
      write(kStdErr,*) 'where kTempJac = use Planck or tau or both '
      write(kStdErr,*) 'when doing d/dT'
      write(kStdErr,*) 'Please reset and retry'
      CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE CheckParams

!************************************************************************
! this sets the temperature variation after nm_radnce is read in
    SUBROUTINE SetkTemperVary(iTemperVary)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input param
    INTEGER :: iTemperVary      !!! from namelist file

! local var
    INTEGER :: iConstOrVary

! this is TEMPERATURE variation in layer
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileUPWELL_LINEAR_IN_TAU
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileDNWELL_LINEAR_IN_TAU
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
!      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      kTemperVary = +1     !!!temperature in layer varies
!      kTemperVary = +2     !!!temperature in layer varies linearly, simple
!      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
!      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM,
!                           !!!  no O(tau) approx, very similar to kTemperVary=4
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM,
!                           !!!  debugged for small O(tau), used with EliMlawer 12/2015
!      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!!  x/6 as x-->0 compared to kTemperVary = +42 *****
!      IF (kFlux .LE. 0) THEN
!        kTemperVary = -1     !!! temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      ELSE
!        kTemperVary = +43    !!! temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!! x/6 as x-->0 compared to kTemperVary = +42 ****
!      END IF

    kTemperVary = iTemperVary
          
    iConstOrVary = -1   !! if do flux, do linear vary T with tau
    iConstOrVary = +1   !! if only RaDTrans, then do constant T in layer, default SARTA/kCARTA for RT only
          
    IF (kFlux <= 0) THEN
      IF (iConstOrVary > 0) THEN
        kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
        iaaOverrideDefault(2,1) = -1
        write(kStdWarn,*) 'kFlux <= 0 so set kTemperVary,iaaOverrideDefault(2,1) to -1'
        write(kStdErr,*)  'kFlux <= 0 so set kTemperVary,iaaOverrideDefault(2,1) to -1'
      ELSEIF (iConstOrVary < 0) THEN
        kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
        iaaOverrideDefault(2,1) = +43
        iaaOverrideDefault(2,4) = +2
        !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
        write(kStdWarn,*) 'kFlux < 0 but set kTemperVary, iaaOverrideDefault(2,1) to 43, iaaOverrideDefault(2,4) to +2'
        write(kStdErr,*)  'kFlux < 0 but set kTemperVary, iaaOverrideDefault(2,1) to 43, iaaOverrideDefault(2,4) to +2'
      END IF
    ELSEIF (kFlux > 0) THEN
      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
      iaaOverrideDefault(2,1) = +43
      iaaOverrideDefault(2,4) = +2
      !!!  x/6 as x-->0 compared to kTemperVary = +42 ****
      write(kStdWarn,'(A)') 'kFlux > 0 so set kTemperVary,iaaOverrideDefault(2,1) to 43; iaaOverrideDefault(2,4) to +2'
      write(kStdErr,'(A)')  'kFlux > 0 so set kTemperVary,iaaOverrideDefault(2,1) to 43; iaaOverrideDefault(2,4) to +2'
    END IF

!!! new, do what the user wishes!!!
    IF ((kFlux <= 0) .AND. (iTemperVary > 0)) THEN
      kTemperVary = +43
    END IF
          
!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
! kTemperVary = +43
! IF (iTemperVary .NE. kTemperVary) THEN
!  write(kStdWarn,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  write(kStdErr,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  kTemperVary = iTemperVary
! END IF
!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
                
    IF (kTemperVary == -1) THEN
      write(kStdWarn,*) 'kTemperVary = -1     !layer temp constant (SARTA DEFAULT)'
    ELSEIF (kTemperVary == +1) THEN
      write(kStdWarn,*) 'kTemperVary = +1     !layer temp varies'
    ELSEIF (kTemperVary == +2) THEN
      write(kStdWarn,*) 'kTemperVary = +2     !layer temp varies linearly, simple v2'
    ELSEIF (kTemperVary == +3) THEN
      write(kStdWarn,*) 'kTemperVary = +3     !layer temp varies linearly, ala LBLRTM v3'
    ELSEIF (kTemperVary == +4) THEN
      write(kStdWarn,*) 'kTemperVary = +4     !layer temp varies linearly, ala LBLRTM v4 O(tau^2)'
    ELSEIF (kTemperVary == +41) THEN
      write(kStdWarn,*) 'kTemperVary = +41    !layer temp varies linearly, ala LBLRTM v4 (Pade)'
    ELSEIF (kTemperVary == +42) THEN
      write(kStdWarn,*) 'kTemperVary = +42    !layer temp varies linearly, ala LBLRTM v4 O(tau)'
    ELSEIF (kTemperVary == +43) THEN
      write(kStdWarn,*) 'kTemperVary = +43    !layer temp varies linearly, ala LBLRTM v4 O(tau) -> tau/6 NEW DEFAULT'
    ELSE
      write(kStdErr,*)'kTemperVary = ',kTemperVary,'unknown option'
      CALL DoStop
    END IF

    iaaOverrideDefault(2,1) = kTemperVary
          
    RETURN
    end SUBROUTINE SetkTemperVary

!************************************************************************
! this function does the solar viewing angle conversion, by Sergio
! given solar zenith angle == angle wrt nadir at Earth surface,
! can compute solar angle at height h : from satzen(h=0) to scanang(h=H)

! ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
! SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

! directly compare to saconv.m
! directly compare to saconv.m
! directly compare to saconv.m

    REAL FUNCTION SACONV_SUN(LSZA, SURFALT, ALT )
                            
    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'
           
! input param
    REAL :: LSZA         !! solar/satellite zenith angle at local surface (which is not necessarily 0)
    REAL :: SURFALT      !! surface altitude ABOVE EARTH surface (km) (eg [6400 Km + ] aircraft at 05 km)
    REAL :: ALT          !! layer altitude   ABOVE EARTH surface (km) (eg [6400 Km + ] layer    at 15 km)

    REAL :: rX,rPi,rEarth

!       rPi = 3.1415927
!       rEarth = 6370.0
           
    rPi = kPi
    rEarth = kPlanetRadius
           
!! p = Snell's law = n(i) r(i) sin(theta(i)) = constant =>
!!     (R+h(i)) sin (theta(i)) = constant (assume n(i) = 1)
!! so Re sin(theta_earth) = (Re+h(i))sin (theta(i))

    rX = (SURFALT + rEarth)/(ALT + rEarth) * sin(LSZA * rPi/180.0)
!! SURFALT < ALT ==> ratio in front of sin(LSZA) < 1
!!   ==> local angles get "smaller" the higher you go

    if (rX > 1.00000) rX = 1.00000000
    if (rX < 0.00000) rX = 0.00000000

    SACONV_SUN = ASIN(rX) * 180/rPi

    RETURN
    END FUNCTION SACONV_SUN

!************************************************************************

! ORIG_SACONV_SUN assumes surfalt == 0       y = ORIG_SACONV_SUN( SZA,          ALT )
! SACONV_SUN      assumes nonzero surfalt    y = SACONV_SUN     (LSZA, SURFALT, ALT )

! directly compare to saconv.m
! directly compare to saconv.m
! directly compare to saconv.m

! this function does the surface --> arb height viewing angle conversion, by Scott
! copied from saconv.f in SARTA package, modified to
!   have input lay height in km
!   have output angle in degrees
! Function to convert the surface solar/satellite zenith angle SZA into the
! local solar/satellite angle at altitude ALT (sunang/scanang)
!    from surface satzen(h=0) to scanang(h=H) at arb altitude H

! NPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SZA     Surface Zenith Angle        degrees


! UTPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL fun  SACONV  local (height H) angle    degrees

!    Function to convert the Zenith Angle SZA at the Earth's
!    surface into a the local angle at altitude ALT.
!    The local angle generally varies slightly with altitude
!    due to the curvature of the Earth and its atmosphere.
!    The effect is largest at the maximum zenith angle, and
!    disappears as the zenith angle approaches 0 degrees.
!    Currently this function only considers the geometry of the
!    situation, and no refractive effects are included.

!    The layers of the atmosphere may be considered as concentric
!    rings with some average altitude. A ray traced thru these rings
!    at any viewing angle other than nadir will have a slightly
!    different angle (relative to the outward radial at the point
!    of intersection) in each ring.

!    The local angle may be calculated (using The Law of
!    Sines) if we know:
!       The solar/satellire zenith angle, SZA, at the Earth's surface (ALT=0)
!       The layer altitude, ALT.
!       The radius of the Earth, RE.

!    The solution uses the law of sines and sin(180 - x) = sin(x)
    REAL FUNCTION ORIG_SACONV_SUN( SZA, ALT )
     
    implicit none
    include '../INCLUDE/TempF90/kcartaparam.f90'
           
    real :: sza,alt

! local variables
    real :: saconv,conv, re, ra

!      ------------------
!      Assign some values
!      ------------------
!      CONV = pi/180 = degrees to radians conversion factor
      CONV=1.7453292E-02
      CONV = 0.017453292519943

!      RE = radius of the Earth (in km)
!      RE=6.37E+03
      RE = kPlanetRadius

!      RA = radius of the point to calc the angle at (in km)
!      Note: layer altitude already in kilometers
      RA = rE + ALT

!      -----------------
!      Do the conversion
!      -----------------

      SACONV=ASIN( (RE/RA) * SIN(CONV*SZA) )
           
! change back to degrees
    ORIG_SACONV_SUN = SACONV/conv

    RETURN
    END FUNCTION ORIG_SACONV_SUN
           
!************************************************************************
!    Function to convert the AIRS satellite viewing angle into the
!    local path angle (      Viewing Angle CONVersion   )
! for downward looking instrument : from scanang(h=H) to satzen(h=0)

! directly compare to vaconv.m
! directly compare to vaconv.m
! directly compare to vaconv.m
! can Matlab debug version, but WARNING its input arguments are switched thusly
!   [zang]=vaconv( sva, salt, alt );
! cd /asl/matlab2012/science
! >> alt = [0:10:80]
! >> [zang]=vaconv(50, 705*1000, alt'*1000 )

! uses law of sines

! NPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SVA     Satellite viewing angle     degrees
!    REAL      SALT    Satellite altitude          km

    REAL FUNCTION VACONV( SVA, ALT, SALT )

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'
           
    REAL :: SVA, ALT, SALT

!      LOCAL VARIABLES
    REAL :: CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink

!      CONV = pi/180 = degrees to radians conversion factor
    CONV = 1.7453292E-02
    CONV = 0.017453292519943

    theta = sva*conv

!      RE = radius of the Earth (in km)
!       RE=6.37E+03
    RE = kPlanetRadius

!      RA = radius of the point to calc the angle at (in km)
    RA = rE + ALT

!      RS = radius of the satellite orbit (in km)
    RS = rE + SALT


!      -----------------
!      Do the conversion
!      -----------------

    rBoink =  (RS/RA) * SIN(CONV*SVA)
    IF (abs(rBoink) > 1.0) rBoink = 1.00000000
           
    RTEMP=ASIN(rBoink)

! change back to degrees

!   write(kSTdWarn,*) 'miaow',RE,ALT,SALT,sva,RTEMP/conv
    VACONV = RTEMP/conv
            
    RETURN
    END FUNCTION VACONV

!************************************************************************
!    Function to convert the AIRS satellite viewing angle into the
!    local path angle (      Viewing Angle CONVersion   ) WITH RefIndex
! for downward looking instrument

! directly compare to vaconv.m, but this uses Snell
! directly compare to vaconv.m, but this uses Snell
! directly compare to vaconv.m, but this uses Snell
! can Matlab debug version, but WARNING its input arguments are switched thusly
!   [zang]=vaconv( sva, salt, alt );
! cd /asl/matlab2012/science
! >> alt = [0:10:80]
! >> [zang]=vaconv(50, 705*1000, alt'*1000 )

! uses law of sines


! NPUT PARAMETERS:
!    type      name    purpose                     units
!    --------  ------  --------------------------  ---------------------
!    REAL      ALT     Average layer altitude      km
!    REAL      SVA     Satellite viewing angle     degrees
!    REAL      SALT    Satellite altitude          km

    REAL FUNCTION VACONV_Snell( SVA, ALT, SALT, rNumberDensity )

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: SVA, ALT, SALT, rNumberDensity

!      LOCAL VARIABLES
    REAL :: CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink
    REAL :: lambda,n0,gamma,ref_ind,mu
    REAL :: rLosch             !Loschmidt number

    rLosch = (kAtm2mb*100/kBoltzmann/273 * 1e-6)    ! (p0/k T0)

!      CONV = pi/180 = degrees to radians conversion factor
    CONV = 1.7453292E-02
    CONV = 0.017453292519943
    theta = sva*conv

!      RE = radius of the Earth (in km)
    RE = 6367.512

! http://mintaka.sdsu.edu/GF/explain/atmos_refr/horizon.html
! http://www.ess.uci.edu/~cmclinden/link/xx/node45.html for TONS of detail
!      RA = radius of the point to calc the angle at (in km)
    RA = rE + ALT

!      RS = radius of the satellite orbit (in km)
    RS = rE + SALT

! ref ind of air at 10 um
    lambda = 10
    n0 = 64.328 + 29498.1/(146-1/lambda/lambda) + 255.4/(41 - 1/lambda/lambda)
    n0 = 1 + n0/1e6
    gamma = (n0*n0-1)/(n0*n0+2)

    IF (ALT >= 200) THEN
      ref_ind = 1.0
    ELSE
      mu = rNumberDensity/rLosch * gamma
      ref_ind = sqrt((2*mu+1)/(1-mu))
    END IF


!      -----------------
!      Do the conversion
!      -----------------

    rBoink =  abs((RS/RA) * SIN(CONV*SVA)/ref_ind)
    if (abs(rBoink) > 1.0) rBoink = 1.00000000

    RTEMP = ASIN(rBoink)
           
! change back to degrees
    VACONV_SNELL = RTEMP/conv
            
    RETURN
    end FUNCTION VACONV_Snell

!************************************************************************

END MODULE n_misc
