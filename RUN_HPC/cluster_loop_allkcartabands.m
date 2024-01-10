%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

set_rtp
set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm

%kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv_H16';
%kcartaexec = '/home/sergio/KCARTA/BIN/Oct26_2016/kcarta.x_lblrtm12.4';
%kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv_H12';

disp('>>>>>')
fprintf(1,'kcartaexec   = %s \n',kcartaexec);
fprintf(1,'f1,f2        = %4i %4i \n',f1,f2);
fprintf(1,'iDoRad       = %2i \n',iDoRad);
if iDoJac > 0
  fprintf(1,'  doing jacs for %3i \n',gg)
else
  disp('  no jacs')
end
fprintf(1,'iDoFlux      = %2i \n',iDoFlux);
fprintf(1,'iDoCloud     = %2i \n',iDoCloud);
fprintf(1,'iDoLBLRTM    = %2i \n',iDoLBLRTM);
fprintf(1,'iDo_rt_1vs43 = %2i \n',iDo_rt_1vs43);
fprintf(1,'iHITRAN      = %2i \n',iHITRAN);
fprintf(1,'iKCKD        = %2i \n',iKCKD);
disp('>>>>>')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
use_this_rtp0 = use_this_rtp;
%set_convolver
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% everywhere you find '/' in use_this_rtp, replace it with '\/'
ooh = strfind(use_this_rtp,'/');
if length(ooh) > 0
  use_this_rtp = strrep(use_this_rtp, '/', '\/');
end

%% everywhere you find '/' in strIceCloud, replace it with '\/'
ooh = strfind(strIceCloud,'/');
if length(ooh) > 0
  strIceCloud = strrep(strIceCloud, '/', '\/');
end

%% everywhere you find '/' in strWaterCloud, replace it with '\/'
ooh = strfind(strWaterCloud,'/');
if length(ooh) > 0
  strWaterCloud = strrep(strWaterCloud, '/', '\/');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%      DATA kaMinFr        /   15.0000,   30.0000,   50.0000,   80.0000,
%c                    140.0000,  300.0000,  500.0000,  605.0000,
%c                    2830.0000,  3550.0000,  5550.0000,  8250.0000,
%c                    12000.0000,  25000.0000  /
%      DATA kaMaxFr        /   30.0000,   50.0000,   80.0000,  150.0000,
%c                    310.0000,  510.0000,  605.0000,  2830.0000,
%c                    3580.0000,  5650.0000,  8400.0000,  12250.0000,
%c                    25000.0000,  44000.0000  /

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 1;
end
% JOB = 49;
% JOB = 1;
% JOB = 269
% JOB = 705;
% JOB = 16
% JOB = 20

[hjunk,hajunk,pjunk,pajunk] = rtpread(use_this_rtp0);
mmw = mmwater_rtp(hjunk,pjunk);
emiss900 = interp1(pjunk.efreq(1:pjunk.nemis(JOB),JOB),pjunk.emis(1:pjunk.nemis(JOB),JOB),900);
surf_properties = [pjunk.stemp(JOB) mmw(JOB) emiss900 pjunk.rlat(JOB) pjunk.rlon(JOB)];

ctype = -9999;
ctype2 = -9999;
iCldORClr = -1;
if isfield(pjunk,'ctype')
  ctype  = pjunk.ctype(JOB);
  ctype2 = pjunk.ctype2(JOB);
  if pjunk.ctype(JOB) > 0 | pjunk.ctype2(JOB) > 0
    iCldORClr = +1;
  else
    iCldORClr = -1;
  end
end

strWaterCloud = '/asl/s1/sergio/CLOUDS_MIEDATA/WATER250/water_405_2905_250';
  strWaterCloud = strrep(strWaterCloud, '/', '\/');
strIceCloud   = '/asl/s1/sergio/CLOUDS_MIEDATA/CIRRUS_BRYANBAUM/v2013/ice_yangbaum_GHM_333_2980_forkcarta';
  strIceCloud = strrep(strIceCloud, '/', '\/');
iNclouds = 2;
iCloudType1 = ctype;
iCloudType1 = 101;
if iCloudType1 == 101
  strCloudType1 = strWaterCloud;
elseif iCloudType1 == 201
  strCloudType1 = strIceCloud;
end
iCloudType2 = ctype2;
iCloudType2 = 201;
if iCloudType2 == 101
  strCloudType2 = strWaterCloud;
elseif iCloudType2 == 201
  strCloudType2 = strIceCloud;
end
sedderC = ' ';
sedderC = [sedderC ' -e "s/iCloudType1/'     num2str(ctype) '/g"  -e "s/strCloudType1/'  strCloudType1 '/g"'];
sedderC = [sedderC ' -e "s/iCloudType2/'    num2str(ctype2) '/g"  -e "s/strCloudType2/'  strCloudType2 '/g"'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kcbands = [15 30 50 80 140 300 500 605 805 2830];

tic
iiBin = JOB;
fluxall = [];
jall = [];
dall = [];
wall = [];

%%% for iBand = length(kcbands)-1 : length(kcbands)-1

for iBand = 1 : length(kcbands)-1
  f1 = kcbands(iBand);
  f2 = kcbands(iBand+1);

  fprintf(1,'%3i iBand  fr1,fr2 = %8.2f %8.2f \n',iBand,f1,f2);

  sedder = ['!sed -e "s/FF1/' num2str(f1) '/g"  -e "s/FF2/' num2str(f2) '/g" '];
  sedder = [sedder ' -e "s/MMM/'     num2str(iiBin) '/g"'];
  sedder = [sedder ' -e "s/TTT/'     num2str(iDo_rt_1vs43) '/g"'];
  sedder = [sedder ' -e "s/XYZXYZ/'  use_this_rtp   '/g"'];  %%%<<<-- to sed rtpfname
  sedder = [sedder ' -e "s/CKDCKD/'  num2str(iKCKD) '/g"'];
  sedder = [sedder ' -e "s/DOLBLRTM/' num2str(iDoLBLRTM) '/g"'];
  if iDoJac == -100
    sedder = [sedder ' -e "s/DOOUTPRESSLEVEL/' num2str(1200) '/g"'];
  elseif iDoJac == +100
    sedder = [sedder ' -e "s/DOOUTPRESSLEVEL/' num2str(0) '/g"'];
  end

  outnml  = ['kcband' num2str(iBand) '.nml_' num2str(iiBin)];
  outname = ['JUNK/kcband' num2str(iBand) '.dat_' num2str(iiBin)];  

  if iDoJac ~= -1
    outjacname = ['JUNK/kcband' num2str(iBand) '.jac_' num2str(iiBin)];  
    outjacnameCOL = [outjacname '_COL'];
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  if exist(outname)
    rmer = ['!/bin/rm ' outname ' ' outjacnameCOL];
    eval(rmer);
  end    
  if iDoJac == -1
    loop_all_kcartabands_radsNFlux
  elseif iDoJac == +1
    disp('not very wise, huge 100 layer jac files will be concatenated together, try column jacs')
    error('loop_all_kcartabands_radsNFlux_100layerjac, see loop_all_kcartabands_radsNFlux_coljac')
  elseif abs(iDoJac) == +100
    loop_all_kcartabands_radsNFlux_coljac
    rmer = ['!/bin/rm ' outname ' ' outjacnameCOL]; eval(rmer);
  end    
end

toc

fprintf(1,'iDoRad = %2i iDoFlux = %2i \n',iDoRad,iDoFlux)

if iDoRad == 3  & iDoFlux < 0
  saver = ['save JUNK/allkcbands_prof' num2str(JOB) '.mat wall dall kcbands'];
elseif iDoRad == 3  & iDoFlux == 5
  saver = ['save JUNK/allkcbands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];
  saver = ['save JUNK/allkc5bands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];  %% originally was called allkcbands which now is for rad calcs only
elseif iDoRad == 3  & iDoFlux == 7
  saver = ['save JUNK/all7kcbands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];
end
if abs(iDoJac) == +100
  %% 9 column jacs : [gid(1 2 3 6 101 102 103)     T ST]  --> 6 column jacs [gid(1 = (1 101 102 103), 2 3 6)     T ST]
  jallx = zeros(length(jall),6);
  jallx(:,1)       = sum(jall(:,[1 5 6 7]),2);
  jallx(:,[2 3 4]) = jall(:,[2 3 4]);
  jallx(:,[5 6])   = jall(:,[8 9]);
  jac_comment = ['gid 1,2,3,6  T,ST   where gid1 = sum of 1,101,102,103'];
  jall = jallx;
  saver = [saver ' jall jac_comment'];
end
surf_prop_comment = 'stemp mmw emiss(900cm-1)  rlat rlon';
saver = [saver ' surf_prop_comment surf_properties'];

eval(saver)

%{
for Oct 2018  AIRS STM
set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
stemp = p.stemp(iiBin);
figure(2); plot(wall,dall,wall,ttorad(wall,p.stemp(iiBin)));
save /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/Clouds_Flux_Rates/kcarta_flux_trop.mat wall dall stemp
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
