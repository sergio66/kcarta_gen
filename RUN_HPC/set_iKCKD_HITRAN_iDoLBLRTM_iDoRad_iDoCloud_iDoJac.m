%iDoLBLRTM = +3; %% use LBLRTM optical depths for CO2, O3, CH4
%iDoLBLRTM = +9999; %% use kcarta_LBLRTM
%iDo_rt_1vs43 = 43;  %% use LINEAR in tau RT (ala LBLRTM)

%iDoCloud = +100; %% yes 100 layer clouds %% bkcarta.x ALWAYS TURNS CLOUDS OFF
%iDoCloud = +1;   %% yes 2slab clouds     %% bkcarta.x ALWAYS TURNS CLOUDS OFF

%iDoJac = +1;   %% do Jacobians

%iDoRad = 0;     %% usual rad
%iDoLBLRTM = +3; %% use LBLRTM optical depths for CO2, O3, CH4
%iDo_rt_1vs43 = 43;  %% use LINEAR in tau RT (ala LBLRTM)

%iDoLBLRTM = +8; %% use 8 ext ODs (1,103,3,4,5,6,9,12) when doing uncertainty originally 2017
%iDoLBLRTM = +7; %% use 8 ext ODs (1,103,3,4,5,6,9     when doing uncertainty again      Apr2018,Oct2019
%iDoLBLRTM = -1; %% use our optical depths

%iDoJac = +1;   %% do Jacobians
%iDoRad = +3;   gg = 5;   %% do rads/jacs/fluxes and if needed, jacobian for this gas

iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = 7;  uncstr = 'Rn';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = 7;  uncstr = 'P+';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = 7;  uncstr = 'B+';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = 7;  uncstr = 'S+';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = 7;  uncstr = 'W+';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iHITRAN = 2016; iDoLBLRTM = 7;  uncstr = 'Rn';  %% W+,S+,B+,P+,Rn, link to template_Qrad_HITRANunc.nml
iHITRAN = 2008; iDoLBLRTM = -1;                %% use our optical depths
iHITRAN = 2008; iDoLBLRTM = 2;                 %% use LBLRTM ODs
iHITRAN = 2012; iDoLBLRTM = -1;                %% use our optical depths
iHITRAN = 2012; iDoLBLRTM = 2;                 %% use LBLRTM ODs
iHITRAN = 2016; iDoLBLRTM = -1;                %% use our optical depths
iHITRAN = 2016; iDoLBLRTM = 2;                 %% use LBLRTM ODs

iKCKD =  6; iHITRAN = 2008; iDoLBLRTM = -1; %% use our optical depths
iKCKD =  6; iHITRAN = 2012; iDoLBLRTM = -1; %% use our optical depths
iKCKD =  6; iHITRAN = 2016; iDoLBLRTM = -1; %% use our optical depths

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 7;  %% the uncertainties, link to template_Qrad_HITRANunc.nml
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 0; %% use UMBC ODs
iKCKD =  1;  iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs
iKCKD =  6;  iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs
iKCKD =  25; iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs
iKCKD =  32; iHITRAN = 2012; iDoLBLRTM = 2; %% use LBLRTM ODs, H12
iKCKD =  32; iHITRAN = 2015; iDoLBLRTM = 2; %% use LBLRTM ODs, G15
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs   ************************************* BEST DEFAULT

%iDoJac = +100; %% do column Jacobians, please set things correctly in template_Qcoljacobian.nml

%iKCKD =  25; iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs  %% for NOAA 2018 meeting
%iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; %% use LBLRTM ODs  %% for NOAA 2018 meeting
%iKCKD =  25; iHITRAN = 2012; iDoLBLRTM = 1; %% use LBLRTM ODs  %% for NOAA 2018 meeting
%iKCKD =  25; iHITRAN = 2016; iDoLBLRTM = 1; %% use LBLRTM ODs  %% for NOAA 2018 meeting and for testing quite a few of the CO2 versions
%iKCKD =  25; iHITRAN = 2015; iDoLBLRTM = 1; %% use LBLRTM ODs  %% for NOAA 2018 meeting + GEISA 2015
%% turn off XSEC for GEISA 2015 / HITRAN 2016 tests for NAOO 2018 meeting
%% iNxsec = -1 --> iNxsec = 0 and iKCKD =  25; iHITRAN = 2016/2015/2012; iDoLBLRTM = 1; %% use LBLRTM ODs for CO2 only

%% turn off/on CO2/WV continuum function : go into template_Qrad.nml and explicitly do this
%               iaaOverrideDefault(1,9) = +2 : do WV/CO2 continuum
%               iaaOverrideDefault(1,9) = +4 : do WV/N2  continuum	       
%               iaaOverrideDefault(1,9) = +6 : do WV/CO2 + WV/N2  continuum	       
%iaaOverride(1,9) = 0     %% default, all off

%iKCKD =  25; iHITRAN = 2016; iDoLBLRTM = 2; iDoCloud = +1; %% use LBLRTM ODs, do clouds for SNOs etc

%uncstr = 'nothing';
%iHITRAN = 2016; iKCKD = 32; iDoLBLRTM = -1;  uncstr = 'XYZ';  %% the uncertainties, link to template_Qrad_HITRANunc.nml
uncstr = 'nothing';

%iDoRad = +1;   gg = 2;   %% do individual gas OD, need to set the .nml template switches yourself (eg cumOD, cumTrans)
%iDoFlux = +5;  %% do GND ILR/ [tropopase toa] OLR only GOOD
%iDoFlux = +7;  %% do [GND tropopase toa] ILR OLR
%iDoRad = +3;   gg = 2;  iDoJac = +1;  %% do rads/jacs/fluxes and if needed, jacobian for this gas
%iDoRad = +3;   gg = 2;  iDoJac = -1;  %% do rads/jacs/fluxes

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 3;     iDoJac = +1; iDoCloud = -1; %% use LBLRTM ODs   ************************************* BEST DEFAULT

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1; iDoCloud = +1; %% use LBLRTM ODs   ************************************* BEST DEFAULT

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2;     iDoRad = 3; iDoJac = +100; iDoCloud = -1; %% use LBLRTM ODs   ***************************** col jac BEST DEFAULT
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 5912;  iDoRad = 3; iDoJac = +1;   iDoCloud = -1; %% use LBLRTM ODs   ************************************* BEST DEFAULT
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2346;  iDoRad = 3; iDoJac = +1;   iDoCloud = -1; %% use LBLRTM ODs   ************************************* BEST DEFAULT
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2456;  iDoRad = 3; iDoJac = +1;   iDoCloud = -1; %% use LBLRTM ODs   ************************************* BEST DEFAULT
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = -1; %% use LBLRTM ODs   ************************************* BEST DEFAULT

%%%%%%% note when I do this,I have changed convolver so only AIRS 2834 chans are done, else waste time convolving!!!  so reset set_convolver.m after this!!!! %%%%%%%%%%
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2346;  iDoRad = 3; iDoJac = +100; iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,3,4,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2456;  iDoRad = 3; iDoJac = +100; iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,4,5,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* CLR JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2346;  iDoRad = 3; iDoJac = +100; iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* COL CLD JACS eg for 64x72 grids : does G2,3,4,6,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 2456;  iDoRad = 3; iDoJac = +100; iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* COL CLD JACS eg for 64x72 grids : does G2,4,5,6,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* CLD JACS 97 layers eg for 64x72 grids

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; iDoRad = 0;  iDoCloud = +1; iDoJac = -1;               %% use LBLRTM ODs   ************************************* BEST DEFAULT, but does "special"
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; iDoRad = 0;  iDoCloud = -1; iDoJac = -1;               %% use LBLRTM ODs   ************************************* BEST DEFAULT, but does "special"
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +100; gg = 2346;  iDoCloud = -1; %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,3,4,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +100; gg = 2456;  iDoCloud = -1; %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,4,5,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +1;   gg = 1001;  iDoCloud = -1; %% clrsky, use LBLRTM ODs   ************************************* CLR JACS 97 layers eg for 64x72 grids

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 0;  iDoCloud = +1; iDoJac = -1;               %% use LBLRTM ODs   ************************************* BEST DEFAULT, but does "special"
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 0;  iDoCloud = -1; iDoJac = -1;               %% use LBLRTM ODs   ************************************* BEST DEFAULT, but does "special"
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +100; gg = 2346; iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,3,4,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +100; gg = 2456; iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* COL CLR JACS eg for 64x72 grids : does G2,4,5,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoJac = +1;   gg = 1001; iDoCloud = -1;  %% clrsky, use LBLRTM ODs   ************************************* CLR JACS 97 layers eg for 64x72 grids

%%%%%%% note when I do this,I have changed convolver so only AIRS 2834 chans are done, else waste time convolving!!!  so reset set_convolver.m after this!!!! %%%%%%%%%%

iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 10; iDoCloud = -1; iDoJac = -1;              %% use LBLRTM ODs   ************************************* BEST DEFAULT, Nalli Emiss, no jacs
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% use LBLRTM ODs   ************************************* BEST DEFAULT               no jacs
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = +1; iDoJac = -1;              %% allsky, use LBLRTM ODs  ****************************** BEST DEFAULT cldsky        no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky        no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = +1; iDoJac = -1;              %% allsky, use LBLRTM ODs  ****************************** BEST DEFAULT cldsky        no jacs
iKCKD =  06; iHITRAN = 2012; iDoLBLRTM = -1; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use UMBC CO2/CH4 ODs  ************************ BEST DEFAULT clrsky        no jacs
iKCKD =  06; iHITRAN = 2008; iDoLBLRTM = -1; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use UMBC CO2/CH4 ODs  ************************ BEST DEFAULT clrsky        no jacs
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky        no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = +1; gg = 1001;   %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = +1; gg = 2346;   %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = +1; gg = 2456;   %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = +1; gg = 5912;   %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky        no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = -1; iDoJac = -1;              %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky        no jacs
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoJac = +100; gg = 2346;iDoCloud = -1;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,3,4,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoJac = +100; gg = 2456;iDoCloud = -1;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,4,5,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2016; iDoLBLRTM = 2;  iDoRad = 3;  iDoJac = +1;   gg = 2;  iDoCloud = -1;   %% use LBLRTM ODs   ************************************* BEST DEFAULT
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2;  iDoRad = 3;  iDoCloud = +1; iDoJac = -1;              %% allsky, use LBLRTM ODs  ****************************** BEST DEFAULT allsky  PCLSAM no jacs
%              iHITRAN = 2012;
%              iHITRAN = 2008;

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 20; iDoCloud = +1; iDoJac = -1;               %% allsky, use LBLRTM ODs, DISORT *********************** DISORT                      no jacs

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 5912;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids this does G 5,9,12
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 2346;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids this doess G 2,3,4,6,51,52
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 2456;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids this doess G 2,4,5,6,51,52

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* CLD JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; gg = 2346;  iDoRad = 3; iDoJac = +100; iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* COL CLD JACS eg for 64x72 grids : does G2,3,4,6,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; gg = 2456;  iDoRad = 3; iDoJac = +100; iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* COL CLD JACS eg for 64x72 grids : does G2,4,5,6,T(z),ST

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids *** USE THIS FOR WV,O3,T ****
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = +1;  %% allsky, use LBLRTM ODs  ************************************* CLD JACS 97 layers eg for 64x72 grids

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids *** USE THIS FOR WV,O3,T ****
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +100; gg = 2346;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,3,4,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +100; gg = 2456;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,4,5,6,51,52,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +100; gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G1,101,101,103,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -100; gg = 1001;  %% clrsky, use LBLRTM ODs  ******* UPLOOK *************** COL CLR JACS eg for 64x72 grids : does G1,101,101,103,T(z),ST

iDoFlux = +5;  %% OLR, ILR and trop
iDoFlux = -1;  %% no fluxes
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -100; gg = 1003;  %% clrsky, use LBLRTM ODs  ******* UPLOOK *************** COL CLR JACS eg for 64x72 grids : does G1,G2,G3,G4,G5,G6,101,101,103,T(z),ST
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;   gg = 1003;  %% clrsky, use LBLRTM ODs  ******* UPLOOK *************** RADS ONLY    eg for 64x72 grids : does G1,G2,G3,G4,G5,G6,101,101,103,T(z),ST

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; gg = 1001;  iDoRad = 3; iDoJac = +1;   iDoCloud = +1;  %% allsky, use LBLRTM ODs   ************************************* CLD JACS 97 layers eg for 64x72 grids
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = +1; iDoJac = -1;               %% allsky, use LBLRTM ODs  ****************************** BEST DEFAULT allsky  PCLSAM no jacs
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 20; iDoCloud = +1; iDoJac = -1;               %% allsky, use LBLRTM ODs, DISORT *********************** DISORT                      no jacs

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 2346;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids this doess G 2,4,5,6,51,52
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 2456;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids this doess G 2,3,4,6,51,52
iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids *** USE THIS FOR WV,O3,T ****

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +1;   gg = 1001;  %% clrsky, use LBLRTM ODs  ****************************** CLR JACS 97 layers eg for 64x72 grids *** USE THIS FOR WV,O3,T ****
iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = +100; gg = 2456;  %% clrsky, use LBLRTM ODs  ****************************** COL CLR JACS eg for 64x72 grids : does G2,G4,G5,G6,G51,G52,T(z),ST

iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
iKCKD =  32; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
iKCKD =  43; iHITRAN = 2024; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
