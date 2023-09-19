%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

system_slurm_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% see eg set_convolver.m
%  iInstr = 1;   % AIRS only
%  iInstr = 2;   % IASI only
%  iInstr = 4;   % CRIS all hi/CHIRP/CrIS lo
%  iInstr = 14;  % AIRS + CRIS all hi/CHIRP/CrIS lo
%  iInstr = 124; % AIRS + IASI + CRIS all hi/CHIRP/CrIS lo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
use_this_rtp = '/home/sergio/MATLABCODE/RANDOM_LARRABEE/bdry_layer_wv.rtp';
use_this_rtp = '/home/sergio/kcartaV118/WORK/test_sergio_fixzobs_clear.rtp';

%% to do cloudy calcs, 34 latbins, so Howard. Chris. Larrabee can test airs->cris, iasi->cris
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS/sartaANDkcarta_op_cloud.rtp';
use_this_rtp = 'sartaANDkcarta_op_cloud.rtp';

%% for JPL George calcs, george aumann one day worth of AIRS data
use_this_rtp = '/home/sergio/MATLABCODE/AUMANN/KCARTA_CLOUDS_May2015/loc.op.rtp';   %% default SARTA/klayers CO2 profile
use_this_rtp = '/home/sergio/PCRTM_XIANGLEI/RUN/pcrtm_co2profile.op.rtp';           %% pcrtm profile upto 1 mb
use_this_rtp = '/home/sergio/PCRTM_XIANGLEI/RUN/pcrtm_co2profile_V2.op.rtp';        %% pcrtm profile upto 0.005 mb

%% for EUMETSAT
use_this_rtp = 'EumetSat_profile_Sept2015/sonde_iasi1.op.rtp';
use_this_rtp = 'EumetSat_profile_Sept2015/ecmwf_iasi1.op.rtp';

%% for flux
use_this_rtp = '/asl/s1/sergio/feb2002_raw_op_airs.rad.constemiss.rtp';
use_this_rtp = 'RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/test.op.rtp';          %% usual 100 AIRS layers
use_this_rtp = 'RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/test_1kmMIPAS.op.rtp'; %% MIPAS 1 km layers

%% clear asky calcs, 7377 profs
use_this_rtp = '/home/sergio/PCRTM_XIANGLEI/RUN/junk_ECM7377_notunclr.rp.rtp';

%% testing LBLRTM vs KVCARTA radiative transfer on 3795 profiles
use_this_rtp = 'JUNK/era_airibrad_day001_random_2015_01_01_tropical_ocean_night.op.rtp';

%% BRDF; check it is from test_brdf1.m
use_this_rtp = '/home/sergio/MATLABCODE/BRDF/junk1.op.rtp';

%% cloudy calcs PCLSAM 100 layers
use_this_rtp_ip = '/asl/s1/sergio/home/PCRTM_XIANGLEI/RUN//forITOVS_ECM_100layercloud.ip.rtp';
use_this_rtp = '/asl/s1/sergio/home/PCRTM_XIANGLEI/RUN//forITOVS_ECM_100layercloud.op.rtp';
iSigmaIASI = -1;   %% use cc = 1     at each level : this is a suggestion to use cc() = 1 at each level then r = tcc * rcloud + (1-tcc) * rclear
iSigmaIASI = +1;   %% use 0 < cc < 1 at each level : this is what Guido/Giuliano are using
iSigmaIASI = 2 ;   %% MRO ala Marco Matricardi

%% maddy CRTM vs kCARTA
use_this_rtp = '../../../IP_PROFILES/eric_maddy.op.rtp';
use_this_rtp = '../../../IP_PROFILES/eric_maddy2.op.rtp';

%% cloudy calcs PCLSAM 2slab
use_this_rtp = '/asl/s1/sergio/home/PCRTM_XIANGLEI/RUN/forITOVS_ECM_2slab.op.rtp';           %% twoslab; cumsum = 999, not so good
use_this_rtp = '/asl/s1/sergio/home/PCRTM_XIANGLEI/RUN/for_kcarta_241_pcrtm_clouds7.op.rtp'; %% twoslab; cumsum = -1, much better

%% testing new sarta
use_this_rtp = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_rp.so2.latlon.const_emiss.rtp';
use_this_rtp = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_rp.so2.latlon.const_emiss_0.8.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm/xyzregr49_1013_400ppm.op.rtp';
use_this_rtp = 'xyzregr49_1013_400ppm.op.rtp';
use_this_rtp = 'cris_hr_ecmwf_csarta_g2s4_d20160120_clear_night_ocean.op.rtp';
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/SAF704_regrprofs_25May2016_xmb_2235g4.op.rtp';
use_this_rtp = 'SAF704_regrprofs_25May2016_xmb_2235g4.op.rtp';

%% for antonia check
use_this_rtp = '/home/sergio/SARTA_CLOUDY/SARTA_CRIS_HiRes/src/test.op.rtp';
%use_this_rtp = '/home/sergio/SARTA_CLOUDY/SARTA_CRIS_HiRes/src/test_820km.op.rtp';

use_this_rtp = '/asl/s1/sergio/TIGR_ALL/tigr1_600_1200.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/AIRS_TUNING/AIRS_TUNING_RTP/pert_tune_tunmlt_jan18_2016_6000tropicalprofiles.op.rtp';
use_this_rtp = 'pert_tune_tunmlt_jan18_2016_6000tropicalprofiles.op.rtp';

%% Larrabee, AIRS STM Mar 2016 and Oct 2017
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49.op.rtp';                     %% scanang = 22
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_400ppm.op.rtp';              %% scanang = 22
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_385ppm.op.rtp';              %% scanang = 22
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_385ppm.2235chans.op.rtp';    %% scanang = 22

use_this_rtp = 'JUNK_HITRAN_385ppmCO2/crisg4_oct16_newjunk49.rp.rtp';
use_this_rtp = 'JUNK_HITRAN_400ppmCO2/crisg4_oct16_newjunk49.rp.rtp';
use_this_rtp = 'JUNK_HITRAN_370ppmCO2_PERTURBS_forITOVS/crisg4_oct16_newjunk49.rp.rtp';

%%% this is the SARTA CrIS HIRes test!!!!!!!
use_this_rtp = '/home/strow/Work/Rta/sarta/test/rtp_drivers/SAF_6angs_704profs_1013mb_unitemis.rtp';
  use_this_rtp = 'RTP/SAF_6angs_704profs_1013mb_unitemis.rtp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% chris is trying to get SARTA-CRIS hires going

%% for background thermal
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm/regr49_1013_400ppm.op.rtp';
use_this_rtp = 'regr49_1013_400ppm.op.rtp';

use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_2surfaces.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_2surfaces.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_1013mb_unitemis.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_1013mb_seaemis.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_1080mb_seaemis.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_1080mb_unitemis.rtp';
use_this_rtp = 'save_SAF_6angs_704profs_1013mb_unitemis.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/save_SAF_6angs_704profs_1013mb_seaemis.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/save_SAF_6angs_704profs_1013mb_unitemis.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/save_SAF_6angs_704profs_1013mb_0p8emis.rtp';
use_this_rtp = 'regr_rtp_6angs_49profs_1013mb_0p8emis.rtp';

%% testing CO2
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_420ppm_2235g4.op.rtp';
use_this_rtp = 'regr49_1013_420ppm_2235g4.op.rtp';
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_410ppm.op.rtp';
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_420ppm.op.rtp';
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_400ppm.op.rtp';

use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1013_2235.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49/regr49_1013_385ppm.op.rtp';
use_this_rtp = 'regr49_1013_385ppm.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA//ECMWF_SAF_137Profiles/save_SAF_704_profiles_25-May-2016_xmb_zeroemiss.op.rtp';
use_this_rtp = 'save_SAF_704_profiles_25-May-2016_xmb_zeroemiss.op.rtp';
use_this_rtp = '/home/chepplew/projects/sarta/cris_hr/regr49_1100_400ppm_2235_e1rh0.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA_SARTA/RUN_KCARTA/testperturb.rp.rtp';

use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/CO_Tuning/cris_mean_hires_a2v4_ref_2018_clear_desc_ocean.op.rtp';
use_this_rtp = 'RTP/cris_mean_hires_a2v4_ref_2018_clear_desc_ocean.op.rtp';

use_this_rtp = 'RTP/barnet_climate.op.rtp';

%% for HITRAN 2018 meeting
use_this_rtp = 'RTP/era_airxbcal_2017_day031_clear_ocean_night.op.rtp';
use_this_rtp = 'RTP/era_retr_era_airxbcal_day031_clear_4.z4sergio_clear_retrieve_allfov.rtp';

%% chris NH3 tests
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Mar2018_NH3/stdNH3_1100mb_op_400ppm.rtp';
use_this_rtp = 'RTP/stdNH3_1100mb_op_400ppm.rtp';
use_this_rtp = 'RTP/stdNH3_1100mb_op_400ppm_zeroemiss.rtp';
use_this_rtp = 'RTP/stdNH3_1100mb_op_400ppm_unitemiss.rtp';
%% strow NH3 test
use_this_rtp = '/home/strow/Work/Rta/cris_testin_nh3.rtp';
use_this_rtp = 'RTP/cris_testin_nh3_x40.rtp';

% chris is testing his new sarta for AIRS
%cp /home/sergio/MATLABCODE/REGR_PROFILES_SARTA_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Sept2018_AIRS2645/regr49_1013_400ppm_unitemiss.op.rtp .
use_this_rtp = 'regr49_1013_400ppm_unitemiss.op.rtp';
use_this_rtp = 'regr49_1100_400ppm_unitemiss.op.rtp';

% chris is testing SNOs
%{
addpath /home/sergio/MATLABCODE/matlib/clouds/sarta
use_this_rtp = '/home/chepplew/data/sno/rtp/sno_sample.rtp';
[h,ha,p,pa] = rtpread(use_this_rtp);
run_sarta.clear = -1;
run_sarta.cloud = -1;
[prof,orig_slabs] = driver_sarta_cloud_rtp(h,ha,p,pa,run_sarta)
h.pfields = 4 + 1;
h.vcmin = 605;
h.vcmax = 2830;
prof.upwell = ones(size(prof.stemp));
prof.zobs = ones(size(prof.stemp)) * 705000;%% no problem, scanag = satzen = 0
rtpwrite('RTP/sno_sample.ip.rtp',h,ha,prof,pa);
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';
klayerser = ['!' klayers ' fin=RTP/sno_sample.ip.rtp fout=RTP/sno_sample.op.rtp'];
%klayerser = ['!' klayers ' fin=/home/chepplew/data/sno/rtp/sno_sample.rtp fout=RTP/sno_sample.op.rtp'];
eval(klayerser)
%}
use_this_rtp = 'RTP/sno_sample.op.rtp';

%cp /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/Aux_jacs_AIRS_August2018_2002_2017/LATS40_avg/Desc/latbin1_40.rp.rtp RTP/.
use_this_rtp = 'RTP/latbin1_40.rp.rtp';

%% this is looking at ENA and TWP sondes, clear sky AIRS obs ... 36 profiles
%% /home/MATLABCODE/CRODGERS_FAST_CLOUD/SONDE_VALIDATION/KCARTA_comparisons
use_this_rtp = 'RTP/ena_twp385.rp.rtp';
use_this_rtp = 'RTP/ena_twp.rp.rtp';

%% testing trapezoid jacs!!!
%use_this_rtp = '/home/sergio/KCARTA/TEST/Test_v120_vs_v121_TRAPEZOIDJACS/trapzpert.rp.rtp';

% trying to compare newest sarta vs kcarta for tropical clear
use_this_rtp = 'RTP/make_365_anom_tropical_profs_for_sartaVSkcarta_radscomparison.op.rtp';
use_this_rtp = 'RTP/make_365_anom_smidlat_profs_for_sartaVSkcarta_radscomparison.op.rtp';
use_this_rtp = 'RTP/make_365_anom_spolar_profs_for_sartaVSkcarta_radscomparison.op.rtp';

%% hmm why does new SARTA clear crash with 7377 set from oerge? (2009/03/01)
use_this_rtp = 'RTP/junkclr.op.rtp';

%% this is for testing various kcarta databases against eg GEISA or HITRAN LM and uncertainty
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_emiss0.985.op.rtp';          %% scanang = 22, sea emiss
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49.op.rtp';                     %% scanang = 22, sea emiss

%% checking flux with cloud and without cloud against RRTM ... works out pretty well for cloud and excellent for clear
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/Clouds_RadiativeKernels/favgPert.op.rtp';
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/Clouds_RadiativeKernels/favgPert_clear.op.rtp';
use_this_rtp = 'RTP/favgPert_clear.op.rtp';
use_this_rtp = 'RTP/favgPert.op.rtp';

%% for IASI NG
use_this_rtp = '/asl/s1/sergio/rtp/rtp_airicrad_v6/2018/06/29/cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.rtp';
use_this_rtp = 'RTP/cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp';
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_0deg.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/Add_2_FFTs/p240random_2018_06_29_clrsky.op.rtp';
use_this_rtp = 'RTP/p240random_2018_06_29_clrsky.op.rtp';

%% smoke
use_this_rtp = 'RTP/australiasmoke_20120_01_04_g040.op.rtp';

%% checking jacobians for AMT 2020 stability paper
%% /home/sergio/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/NewJacs_forAMT2020
use_this_rtp = 'RTP/pbeforeavg.op.rtp';
use_this_rtp = 'RTP/pafteravg.op.rtp';

%% to do the TWP tile .... run off single footprint clear code, 
%% moo = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/LookAtTimeSeries/RTP_PROFSV2/Cld/CRODGERS_FAST_CLOUD_Retrievals/V5Jacs//jacsV5_retr_gran_lonbin_67_latbin_35_JOB_2515.mat');
%% rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/tile_lonbin_67_latbin_35_JOB_2515_clear.rtp',moo.hoem,[],moo.poem,[]);
%% can try llittle cloudy
%% moo = load('/home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/LookAtTimeSeries/RTP_PROFSV2/Cld/CRODGERS_FAST_CLOUD_Retrievals/V4Jacs//jacsV4_ret>> ran_lonbin_67_latbin_35_JOB_2515.mat');
use_this_rtp = 'RTP/tile_lonbin_67_latbin_35_JOB_2515_clear.rtp';

%% test nicknalli emis, symbolic link to (Summer 2020)
%% /home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junk.rp.rtp
use_this_rtp = 'RTP/testnalli_masudaemis.rtp';    %% masuda emis
use_this_rtp = 'RTP/testnalliemis_resetTWV.rtp';  %% new emis, reset ST, CO2,CH4, WV
use_this_rtp = 'RTP/testnalliemis0.rtp';           %% new emis, June 2020
use_this_rtp = 'RTP/testnalliemis1.rtp';           %% new emis, Aug 2020
use_this_rtp = 'RTP/pnalli_aug30_2020_masuda19_65profs.op.rtp'; %% smaller subset, 65 profs
use_this_rtp = 'RTP/pnalli_aug30_2020_nalli30_65profs.op.rtp';  %% smaller subset, 65 profs
%% test nicknalli emis, symbolic link to (Nov2020)
%% /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/BRDF_EMISSIVITY_NALLI/TestNalliEMiss/
use_this_rtp = 'RTP/nalli_resetTWV_era_Nov24_2020.rp.rtp';  %% new emis
use_this_rtp = 'RTP/nalli_era_Nov24_2020.rp.rtp';  %% new emis
use_this_rtp = 'RTP/masuda_resetTWV_era_Nov24_2020.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/masuda_era_Nov24_2020.rp.rtp'; %% masuda emis
%% /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/BRDF_EMISSIVITY_NALLI/TestNalliEMiss/
use_this_rtp = 'RTP/masuda_resetTWV_era_Feb23_2021.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/masuda_era_Feb23_2021.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/nalli_resetTWV_era_Feb23_2021.rp.rtp';  %% new emis
use_this_rtp = 'RTP/nalli_era_Feb23_2021.rp.rtp';  %% new emis
%% /umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/BRDF_EMISSIVITY_NALLI/TestNalliEMiss/
%% test nicknalli emis, symbolic link to (Mar2021)
use_this_rtp = 'RTP/nalli_era_Mar17_2021.rp.rtp';  %% new emis 
use_this_rtp = 'RTP/masuda_era_Mar17_2021.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/nalli_resetTWV_era_Mar17_2021.rp.rtp';  %% new emis
use_this_rtp = 'RTP/masuda_resetTWV_era_Mar17_2021.rp.rtp'; %% masuda emis
%% test nicknalli emis, symbolic link to (Apr2021),5000 profiles
use_this_rtp = 'RTP/nalli_era_Mar17_2021.rp.rtp';  %% new emis 
use_this_rtp = 'RTP/masuda_era_Mar17_2021.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/nalli_resetTWV_era_Mar17_2021.rp.rtp';  %% new emis
use_this_rtp = 'RTP/masuda_resetTWV_era_Mar17_2021.rp.rtp'; %% masuda emis
%% test nicknalli emis, symbolic link to (Mar2021)
use_this_rtp = 'RTP/nalli_resetTWV_era_Apr23_2021.rp.rtp';  %% new emis
use_this_rtp = 'RTP/nalli_era_Apr23_2021.rp.rtp';  %% new emis 
use_this_rtp = 'RTP/masuda_era_Apr23_2021.rp.rtp'; %% masuda emis
use_this_rtp = 'RTP/masuda_resetTWV_era_Apr23_2021.rp.rtp'; %% masuda emis

%% testing H2016_LBMRTMco2 vs H2008_UMBC_CO2
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49.op.rtp';

%% Chris H. perturbed profiles Jan 2019
%% /home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_unitemis_nadir_pert.rtp
use_this_rtp = 'RTP/r49_1100_400p_unitemis_nadir_pert.rtp';
%% /home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_unitemis_8angs_pert.rtp
use_this_rtp = 'RTP/r49_1100_400p_unitemis_8angs_pert.rtp';
%% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_ep75_7angs_unpert.rtp
use_this_rtp = 'RTP/r49_1013_400p_ep75_7angs_unpert.rtp';
%% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_e1_7angs_unpert.rtp
use_this_rtp = 'RTP/r49_1013_400p_e1_7angs_unpert.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_unitemis_seaemis_7angs_night.rtp
use_this_rtp = 'RTP/r49_1013_400p_unitemis_seaemis_7angs_night.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1100_98lev_400p_unitemis_seaemis_7angs_night.rtp
use_this_rtp = 'RTP/r49_1100_98lev_400p_unitemis_seaemis_7angs_night.rtp';
% home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_seaemis_8angs_pert_v2.rtp
use_this_rtp = 'RTP/r49_1100_400p_seaemis_8angs_pert_v2.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_seaemis_8angs_pert_v2.rtp
use_this_rtp = 'RTP/r49_1100_400p_seaemis_8angs_pert_v2.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert.rtp
use_this_rtp = 'RTP/r49_1013_400p_8angs_sfc_pert.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert_v2.rtp
use_this_rtp = 'RTP/r49_1013_400p_8angs_sfc_pert_v2.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert_815zobs.rtp
use_this_rtp = 'RTP/r49_1013_400p_8angs_sfc_pert_815zobs.rtp';
%%%%
%%% for CHIRP
use_this_rtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert.rtp';
% /home/chepplew/data/sarta/prod_2019/generic/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp
  use_this_rtp = 'RTP/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2020/generic/r49_1100_400p_unitemis_8angs_gas_pert_v1.rtp';
  use_this_rtp = 'RTP/r49_1100_400p_unitemis_8angs_gas_pert_v1.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_8angs_sfc_pert_v2.rtp';
  use_this_rtp = 'RTP/r49_1013_400p_8angs_sfc_pert_v2.rtp';

use_this_rtp = '/home/sergio/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES_AIRS_IASI_CRIS/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = 'RTP/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = '/asl/val/airs_2009/ARMTWP/Phase1/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
% chris files
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/armtwp_nearby_semiclear_jan04_resetT_nte_co2.rtp';       % 1
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/tairs_382co2_105o3_104n2o_070co_101ch4_050hno3_rad.rtp'; % 2
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/voem_start.rtp';                                         % 3
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp3_start.rtp';                                         % 4
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp2_start.rtp';                                         % 5
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp1_start.rtp';                                         % 6
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp3_start.rtp';                                         % 7
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp2_start.rtp';                                         % 8
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp1_start.rtp';                                         % 9
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/minn_start.rtp';                                         % 10
use_this_rtp = '/home/chepplew/data/sarta_tuning/frm_scott/tiasi_382co2_105o3_104n2o_101ch4_050hno3_rad.rtp_1';             % 11
%
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp2_start.rtp';                                         % 8
use_this_rtp = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp2_start.rtp';                                         % 5
use_this_rtp = '/home/chepplew/data/sarta/validation/2018d152G013_levs.rtp';
% chris h tests of sarta
use_this_rtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_seaemis_8angs_pert_v2.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1100_400p_unitemis_8angs_pert.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2019/generic/r49_1013_400p_unitemis_7angs_night.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/r49_1013_400p_unitemis_7angs_night.rtp';
use_this_rtp = 'RTP/r49_1013_400p_unitemis_7angs_night.rtp';

%% larrabee wants col gas jacs (actually stemp only)
use_this_rtp = 'RTP/latbin20_45angles.op.rtp';
use_this_rtp = 'RTP/pert_testKC.op.rtp';   %% from ~/MATLABCODE/oem_pkg_run/AIRS_new_clear_scan_August2019/NewStuff_forAMT2020_Response2Reviewers/NewJacs/
use_this_rtp = 'RTP/latbin1_40.op.rtp';
use_this_rtp = 'RTP/latbin1_40.op_400ppm.rtp';

%% testing ecRad vs RRTM flux
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/testRRTM_ECRAD.rtp';

%% allsky trends for tiles : see
%%   /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/make_summary_latbin_files_txt.m
%%   these start in Jan 1 of every year
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019.rtp';
use_this_rtp = 'RTP/summary_17years_all_lat_all_lon_2002_2019.rtp';
use_this_rtp = 'RTP/summary_17years_all_lat_all_lon_2002_2019_palts.rtp';
%%
%% /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/call_save_split_apart_rtp_howard_bins_startSept2002.m
%%    these start Sept 1, 2002 and go on
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp';
use_this_rtp = 'RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002.rtp';
use_this_rtp = 'RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp';
%% this is latbin 44, lon 12
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/LatBin44/summary_latbin_44_lonbin_12.rtp'; %% 388 steps
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/LatBin32/summary_latbin_32_lonbin_12.rtp'; %% 388 steps
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/LatBin64/summary_latbin_64_lonbin_12.rtp'; %% 388 steps
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin01/all12monthavg_T_WV_grid_latbin_01_lonbin_12.rtp'; %% 25 T x WV grids
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin64/all12monthavg_T_WV_grid_latbin_64_lonbin_12.rtp'; %% 25 T x WV grids
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin32/all12monthavg_T_WV_grid_latbin_32_lonbin_12.rtp'; %% 25 T x WV grids
%% alllsky profiles, averaged over longitudes
use_this_rtp = '/~/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR_zonalavg/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/simulate64binsERA5_32.rp.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/simulate64binsERA5_32.op.rtp'; %% 16416 profiles, no need to do ANY jacs
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/SimulateTimeSeries/simulate64binsERA5_15.op.rtp'; %% 16416 profiles, no need to do ANY jacs

use_this_rtp = '/home/chepplew/data/sarta/validation/sng_2020_subs_for_kcarta.rtp';
use_this_rtp = '/home/chepplew/data/sarta/validation/sng_2020_subs_for_kcarta_v2.rtp';
use_this_rtp = 'RTP/ecmwf_airicrad_day268_2021clear.rtp';

use_this_rtp = '/home/chepplew/data/scratch/iasi_35987367.kla_1';
use_this_rtp = '/home/sergio/KCARTA/IP_PROFILES/junk49_400ppm.op.rtp';              %% scanang = 22
use_this_rtp = '/home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/make_plot_PCTS.op.rtp';

iInstr = 1;
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_nocld_USSTD.op.rtp';          %% testing DISORT
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_watercld_USSTD.op.rtp';       %% testing DISORT
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_tangtesticecld_USSTD.op.rtp'; %% testing DISORT
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_icecld_USSTD.op.rtp';         %% testing DISORT
use_this_rtp = '/home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/parametrize_pclsam_disort_AFGL_1.rp.rtp';
use_this_rtp = '/home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/parametrize_pclsam_disort_AFGL_1_49.rp.rtp';
use_this_rtp = '/home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/parametrize_pclsam_disort_arbCZT_AFGL_1_49.rp.rtp';
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_iceNwatercld_USSTD_7cngwat.op.rtp'; %% testing PCLSAM with Chou adjustment 
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_watercld_USSTD_7cngwat.op.rtp';     %% testing PCLSAM with Chou adjustment
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_icecld_USSTD_7cngwat.op.rtp';       %% testing PCLSAM with Chou adjustment
use_this_rtp = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_watercld_USSTD_7cngwat.op.rtp';     %% testing PCLSAM with Chou adjustment
use_this_rtp = '/asl/s1/sergio/forITOVS_May2019/KCARTA_PCRTM_Tang//pcrtm_11.01.064_7chans_15min.op.rtp';
use_this_rtp = '/asl/s1/sergio/forITOVS_May2019/KCARTA_PCRTM_Tang//pcrtm_11.01.055_7chans_15min.op.rtp';
use_this_rtp = '/asl/s1/sergio/forITOVS_May2019/KCARTA_PCRTM_Tang//pcrtm_10.31.215_7chans_15min.op.rtp';

iInstr = 124;
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA_SARTA/ECMWF_SAF_137Profiles/save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/save_SAF_704_profiles_29-Apr-2016_1100mb_400ppmv_unitemis.op.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/test_3profs.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/test_9profiles.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/tigr_cris_fsr.rtp';
use_this_rtp = '/home/chepplew/data/sarta/prod_2022/generic/test_pert_gas1.rtp';

iInstr = 124;
use_this_rtp = '/home/chepplew/data/Sergio/airs_2018d259.op.rtp';    %% comparing AIRS CriS IASI clear, 30000 fovs
use_this_rtp = '/home/chepplew/data/Sergio/iasi2_20180916.op.rtp_1'; %% comparing AIRS CriS IASI clear, 34000 fovs
use_this_rtp = '/home/chepplew/data/Sergio/j1_20190916.op.rtp';      %% comparing AIRS CriS IASI clear, 40500 fovs

iIntr = 1;
use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_12years_all_lat_all_lon_2002_2014_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_07years_all_lat_all_lon_2012_2019_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
use_this_rtp = '/asl/s1/sergio/MakeAvgObsStats2002_2020_startSept2002_v3/TimeSeries/ERA5/Tile_Center12months/DESC/2012/FixedNAN/all4608_era5_full12months_Qcumulative09.rtp'; % BT1231cld quants = [0 0.03 0.05 0.1 0.2 0.5 0.8 0.9 0.95 0.97 1.0];

clear iInstr iDoConvolve
iInstr = 4;
use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/ChangeJPSSTilt/junk2A.op.rtp';                                     iInstr = 4;
use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/ChangeJPSSTilt/junk.op.rtp';                                       iInstr = 4;
use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/ChangeJPSSTilt/junk3.op.rtp';                                      iInstr = 4;
use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/ChangeJPSSTilt/junk4.op.rtp';                                      iInstr = 4;
use_this_rtp = '//umbc/xfs2/strow/asl/s1/sergio/home/MATLABCODE/QUICKTASKS_TELECON/DaveTobin_AMS2023/dtobin_ams2023.rp.rtp';  iInstr = 4;

clear iInstr iDoConvolve
use_this_rtp = '/home/chepplew/data/scratch/mktemp_fnnPyDLB_cris.op.rtp'; iInstr = 0; iDoConvolve = -1;  %% only kcarta outputs, no convolve

clear iInstr iDoConvolve
iInstr = 124; use_this_rtp = '/home/sergio/KCARTA/WORK/wierd_jpss2.op.rtp'; 

clear iInstr iDoConvolve
iInstr = 2; iDoConvolve = +1;
use_this_rtp = '/home/chepplew/data/sarta/validation/iasi2_2020d081g099_op.rtp';
use_this_rtp = '/home/chepplew/data/sarta/validation/iasi_1_2016d246_clear_sea_gn173_dep108.rtp_1';
use_this_rtp = '/home/chepplew/data/sarta/validation/iasi_1_2016d246_clear_sea_gn173_dep0.rtp_1';
use_this_rtp = '/home/chepplew/data/sarta/validation/hdo/hdo_test_prf_frm_jpl.rtp_1';
use_this_rtp = '/home/chepplew/data/scratch/mktemp_N2dzVSCD_airs_001_kl.rtp';
use_this_rtp = '/home/chepplew/data/scratch/mktemp_N2dzVSCD_airs_001_kl.rtp_1';

clear iInstr iDoConvolve
iInstr = 4; iDoConvolve = +1;
use_this_rtp = '/home/sergio/MATLABCODE/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm/xyzregr49_1013_400ppm.op.rtp';
use_this_rtp = '/home/sergio/MATLABCODE_Git/REGR_PROFILES_SARTA/RUN_KCARTA/REGR49_400ppm_H2016_Feb2020_AIRS2834_CHIRP/regr49_1013_400ppm_unitemiss.op.rtp';

clear iInstr iDoConvolve
iInstr = 1; iDoConvolve = +1;
use_this_rtp = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac/newdayx_1_100_12150.op.rtp';
use_this_rtp = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac/cloudy_airs_l1c_ecm_sarta_baum_ice.2018.06.29.086_cumsum_-1.op.rtp';
use_this_rtp = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/srcF77_jac/newdayx_clr.op.rtp';
use_this_rtp = '/home/chepplew/data/scratch/mktemp_gLJIFbf4_kl2.op.rtp';
use_this_rtp = '/home/chepplew/data/scratch/mktemp_xwCLnR3G_airs_l1c__op.rtp';

clear iInstr iDoConvolve
iInstr = 1; iDoConvolve = +1;
%% ERA5 desc avg
%% >>>>>>>>>> see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_makeavgprofile_ERA5_monthly_desc_or_asc.m
%% >>>>>>>>>> see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_makeavgprofile_ERA5_monthly_desc_or_asc.m
use_this_rtp = 'RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp';
use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp'; %% actually has clouds
use_this_rtp = 'RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_12years_all_lat_all_lon_2002_2014_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_07years_all_lat_all_lon_2012_2019_monthlyERA5.rp.rtp';
use_this_rtp = 'RTP/summary_atm_N_cld_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';  %% clouds
%%%
use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';            %% clear, raw,  see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp';       %% clear, pert, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m
%%% 
%% >>>>>>>>>> see /home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/driver_makeavgprofile_ERA5_monthly_desc_or_asc.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% moved to clust_do_kcarta_driver.m
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
%% use_this_rtp0 = use_this_rtp;
%% set_convolver
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
%%
%% see eg set_convolver.m
%%  iInstr = 1;   % AIRS only
%%  iInstr = 2;   % IASI only
%%  iInstr = 4;   % CRIS all hi/CHIRP/CrIS lo
%%  iInstr = 14;  % AIRS + CRIS all hi/CHIRP/CrIS lo
%%  iInstr = 124; % AIRS + IASI + CRIS all hi/CHIRP/CrIS lo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
