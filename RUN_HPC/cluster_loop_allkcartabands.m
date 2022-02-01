%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools

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
% JOB = 49;
% JOB = 1;
% JOB = 269
% JOB = 705;
% JOB = 16
% JOB = 20

[hjunk,hajunk,pjunk,pajunk] = rtpread(use_this_rtp0);
ctype  = pjunk.ctype(JOB);
ctype2 = pjunk.ctype2(JOB);

if pjunk.ctype(JOB) > 0 | pjunk.ctype2(JOB) > 0
  iCldORClr = +1;
else
  iCldORClr = -1;
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
  
  outnml  = ['kcband' num2str(iBand) '.nml_' num2str(iiBin)];
  outname = ['JUNK/kcband' num2str(iBand) '.dat_' num2str(iiBin)];  

  if iDoRad == 3  & iDoFlux < 0
    if iCldORClr == -1 | iBand < 7
      %% can only do clear
      sedder = [sedder ' template_Qrad_allbands.nml  > ' outnml];
    elseif iCldORClr == +1 & iBand >= 7
      sedder = [sedder ' ' sedderC ' '];
      sedder = [sedder ' template_Qrad_2cloud_allbands.nml  > ' outnml];
    end
    eval(sedder)
    kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
    if exist(outname)
      rmer = ['!/bin/rm ' outname];
      eval(rmer);
    end    
    eval(kcartaer)
  
    if iCldORClr == -1  | iBand < 7
      %% can only do clear
      [d,w] = readkcstd(outname);
      dall = [dall; d];
      wall = [wall w];
    elseif iCldORClr == +1 & iBand >= 7
      %cldparamsname = [outname '_CLD'];
      %cunk = load(cldparamsname);
      %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

      [d,w] = readkcstd(outname);
      [mmjunk,nnjunk] = size(d);
      dall = [dall; d(:,nnjunk)];
      wall = [wall w];
    end

  elseif iDoRad == 3 & iDoFlux == 5  
    if iCldORClr == -1 | iBand < 7
      %% can only do clear
      sedder = [sedder ' template_Qflux5_allbands.nml  > ' outnml];
    elseif iCldORClr == +1 & iBand >= 7
      sedder = [sedder ' ' sedderC ' '];
      sedder = [sedder ' template_Qflux5_2cloud_allbands.nml  > ' outnml];
    end
    eval(sedder)
    kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
    if exist(outname)
      rmer = ['!/bin/rm ' outname];
      eval(rmer);
    end    
    eval(kcartaer)

    if iCldORClr == -1 | iBand < 7
      %% can only do clear
      [d,w] = readkcstd(outname);
      dall = [dall; d];
      wall = [wall w];

      [fx,w] = readkcflux([outname '_OLR3']);
      fluxall = [fluxall; fx];

    elseif iCldORClr == +1 & iBand >= 7  
      %cldparamsname = [outname '_CLD'];
      %cunk = load(cldparamsname);
      %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

      [d,w] = readkcstd(outname);
      [mmjunk,nnjunk] = size(d);
      dall = [dall; d(:,nnjunk)];
      wall = [wall w];

      [fx,w] = readkcflux([outname '_OLR3']);
      fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
    end

  elseif iDoRad == 3 & iDoFlux == 7  
    if iCldORClr == -1 | iBand < 7
      %% can only do clear
      sedder = [sedder ' template_Qflux7_allbands.nml  > ' outnml];
    elseif iCldOrClr == +1 &  iBand >= 7
      sedder = [sedder ' ' sedderC ' '];
      sedder = [sedder ' template_Qflux7_2cloud_allbands.nml  > ' outnml];
    end
    eval(sedder)
    kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
    if exist(outname)
      rmer = ['!/bin/rm ' outname];
      eval(rmer);
    end    
    eval(kcartaer)

    if iCldORClr == -1 | iBand < 7
      %% can only do clear
      [d,w] = readkcstd(outname);
      dall = [dall; d];
      wall = [wall w];

      [fx,w] = readkcflux([outname '_IOLR3']);
      fluxall = [fluxall; fx];

    elseif iCldORClr == +1 & iBand >= 7  
      %cldparamsname = [outname '_CLD'];
      %cunk = load(cldparamsname);
      %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

      [d,w] = readkcstd(outname);
      [mmjunk,nnjunk] = size(d);
      dall = [dall; d(:,nnjunk)];
      wall = [wall w];

      [fx,w] = readkcflux([outname '_IOLR3']);
      fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
    end
  
  else
    fprintf(1,'iDoRad = %3i iDoFlux = %3i \n',iDoRad,iDoFlux)
    error('huh?? unknown combo of iDoRad/iDoFlux!!!')
  end
  
  figure(2); plot(wall,rad2bt(wall,dall)); pause(0.1)
  rmer = ['!/bin/rm ' outnml ' status' num2str(iiBin)]; eval(rmer)

end
toc

if iDoRad == 3  & iDoFlux < 0
  saver = ['save JUNK/allkcbands_prof' num2str(JOB) '.mat wall dall kcbands'];
elseif iDoRad == 3  & iDoFlux == 5
  saver = ['save JUNK/allkcbands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];
  saver = ['save JUNK/allkc5bands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];  %% originally was called allkcbands which now is for rad calcs only
elseif iDoRad == 3  & iDoFlux == 7
  saver = ['save JUNK/all7kcbands_prof' num2str(JOB) '.mat wall dall fluxall kcbands'];
end
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
