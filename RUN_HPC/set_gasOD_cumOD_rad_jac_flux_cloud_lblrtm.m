%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = 605; f2 = 2830; mm = -1;   %% dummies, needed for now

iDoRad = +0;             %% do usual rad
iDoRad = +10;            %% do usual rad, test NNalli Emis : backgnd thermal uses p.satzen for kThermalAngle, not acos(3/5)
iDoRad = +1;   gg = 3;   %% do individual gas OD, need to set the .nml template switches yourself (eg cumOD, cumTrans)
iDoRad = +2;             %% do cumulative gas OD
iDoRad = +3;   gg = 4;   %% do rads/jacs/fluxes and if needed, jacobian for this gas

iDoJac = +100; %% do column Jacobians, dnlook from TOA
iDoJac = -100; %% do column Jacobians, uplook from GND
iDoJac = +1;   %% do Jacobians
iDoJac = -1;   %% do rads and mebbe ODs

%% see definitions of kFlux
iDoFlux = +6;  %% do up/down fluxes
iDoFlux = +2;  %% do heating rates
iDoFlux = -1;  %% no fluxes
iDoFlux = +5;  %% do ILR/OLR only

iDoCloud = +1;   %% yes TwoSlabclouds    %% bkcarta.x ALWAYS TURNS CLOUDS OFF
iDoCloud = +100; %% yes 100 layer clouds %% bkcarta.x ALWAYS TURNS CLOUDS OFF
iDoCloud = -1; %% no clouds

iDoLBLRTM = +9999; %% use kcarta_LBLRTM
iDoLBLRTM = -1; %% use our optical depths
iDoLBLRTM = +1; %% use LBLRTM optical depths for CO2
iDoLBLRTM = +2; %% use LBLRTM optical depths for CO2, CH4
iDoLBLRTM = +3; %% use LBLRTM optical depths for CO2, O3, CH4
iDoLBLRTM = +4; %% use LBLRTM optical depths for CO2, O3, CH4, CO
iDoLBLRTM = +5; %% use LBLRTM optical depths for CO2, O3, CH4, CO, N2
iDoLBLRTM = +6; %% use LBLRTM optical depths for CO2, O3, CH4, CO, N2, O2

iDoLBLRTM = -1; %% use our optical depths
iDoLBLRTM = +3; %% use LBLRTM optical depths for CO2, O3, CH4
iDoLBLRTM = +2; %% use LBLRTM optical depths for CO2, CH4

iDo_rt_1vs43 = 43;  %% use LINEAR in tau RT (ala LBLRTM)
iDo_rt_1vs43 = -1;  %% use CONST  in tau RT (ala SARTA)

iHITRAN = 2008;
iHITRAN = 2012;
iHITRAN = 2016;   

iKCKD =  1;    %% MT CKD  1, about 2003
iKCKD =  6;    %% MT CKD  6, Sergio/Scott mod to MT CKD 1, about 2008
iKCKD = 25;    %% MT CKD 25, about 2015
iKCKD = 32;    %% MT CKD 32, about 2018
iKCKD = 43;    %% MT CKD 32, about 2026

iArb_RADatPLEV = -1;  %% dump rad at TOA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%                      THESE ARE DEFAULT            %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_iKCKD_HITRAN_iDoLBLRTM_iDoRad_iDoCloud_iDoJac %% default is CKD43,HITRAN2024,rad and no jacs .... override below

iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids *** USE THIS FOR WV,O3,T ****
iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +100; gg = 2456;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,G4,G5,G6,G51,G52,T(z),ST

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
iKCKD =  32; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs

iArb_RADatPLEV = 250.00;   %% if you want to dump rads out at eg 250 mb instead of default
iArb_RADatPLEV = -1;       %% if you want to dump rads out at default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_kcarta_exec_iHITRAN 

if iDoLBLRTM == 9999
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_lblrtm12.4';
  disp('>>>>>>>>>>>> hope you chose CKD 25 !!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
end

iDoDefault = +1;
iDoDefault = 2008;
iDoDefault = -1;

if iDoDefault == 2008
  iDoRad = +3;       %% do rads/jacs/fluxes
  iDoJac = -1;       %% do rads and mebbe ODs
  iDoFlux = -1;      %% no fluxes
  iDoCloud = -1;     %% no clouds
  iDoLBLRTM = -1;    %% use our optical depths
  iDo_rt_1vs43 = -1; %% const in tau radiative transfer
  kcartaexec = '/home/sergio/KCARTA/BIN/Oct10_2016/bkcarta_H2008.x';
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta_H2008.x';
  disp('>>>>>>>>>>>> doing H2008 kcarta --- hope you chose CKD 1 or CKD 6 !!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')  
elseif iDoDefault == 1
  iDoRad = +3;       %% do rads/jacs/fluxes
  iDoJac = -1;       %% do rads and mebbe ODs
  iDoFlux = -1;      %% no fluxes
  iDoCloud = -1;     %% no clouds
  iDoLBLRTM = -1;    %% use our optical depths
  iDo_rt_1vs43 = -1; %% const in tau radiative transfer
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x';  
end
if iDoRad == 20
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20';
end

%%kcartaexec    = '/home/sergio/KCARTA/BIN/kcarta.x90';
%%kcartaexec    = '/home/sergio/KCARTA/BIN/kcarta.x90_400ppmv_H16';
%%kcartaexec    = '/home/sergio/KCARTA/BIN/bkcarta.x_f90_120_400ppmv_H16';

strIceCloud   = '/asl/s1/sergio/CLOUDS_MIEDATA/CIRRS_PYANG_MODIS_CERES/COARSE_RRTM/kcarta_200_3000_pingyang_modisL2.dat';
%% <<<<<<<<<< this is what SARTA uses >>>>>>>>>>
strWaterCloud   = '/asl/s1/sergio/CLOUDS_MIEDATA/WATER250/water_405_2905_250';   %% this is what SARTA uses
strIceCloud     = '/asl/s1/sergio/CLOUDS_MIEDATA/CIRRUS_BRYANBAUM/v2013/ice_yangbaum_GHM_333_2980_forkcarta';
%% <<<<<<<<<< this is what SARTA uses >>>>>>>>>>
%strWaterCloud = '/asl/s1/sergio/CLOUDS_MIEDATA/WATER_290_GAMMADIST/waterALL_105_3005_290';              %% testing FIR
%strWaterCloud = '/asl/s1/sergio/CLOUDS_MIEDATA/Pengwang_Zhai/WATER/water_pengwang_600_3000_forkcarta';  %% testing Pengwang Zhai SOS

caComment = [date ' ' kcartaexec ' iDoRad=' num2str(iDoRad) ' iDoLBLRTM=' num2str(iDoLBLRTM) ' iDo_rt_1vs43=' num2str(iDo_rt_1vs43,'%02d')];
caComment = [caComment ' iDoCloud=' num2str(iDoCloud)];
ooh = strfind(caComment,'/');
if length(ooh) > 0
  caComment = strrep(caComment, '/', '\/');
end
if length(caComment) > 160
  disp('warning : set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m is truncating caComment to 160 chars')
  caComment = caComment(1:160);
end
