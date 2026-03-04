%% % /bin/rm slurm* JUNK/rad.dat*; ; sbatch --array=G1-G2 sergio_matlab_jobB.sbatch

%% need to modify template_QXYZ.nml CORRECTLY for the rtp file to process!

%addpath /home/sergio/MATLABCODE
addpath /home/sergio/git/matlabcode
addpath /home/sergio/KCARTA/MATLAB

% system_slurm_stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iInstr = 1; iDoConvolve = 1;
use_this_rtp = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1013_400ppm_unitemiss.op.rtp';

kcartaexec20 = '/home/sergio/KCARTA/BIN/kcarta.x90_v1.22_400ppmv_H20';                %% H20, v1.22, UMBC CO2
kcartaexec24 = '/home/sergio/KCARTA/BIN/kcarta.x90_v1.22_400ppmv_H24';                %% H24, v1.22, UMBC CO2

%% everywhere you find '/' in use_this_rtp, replace it with '\/'
ooh = strfind(use_this_rtp,'/');
if length(ooh) > 0
  use_this_rtp = strrep(use_this_rtp, '/', '\/');
end
MMM = 12; %% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/TEST_H2020_H2024_CKD32_CKD43/compare_the_H16_H20_H24_hitranversions_with_kcarta.m which is copied into this here directory, pardner

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% there are 71 gases to get rid of, see the lists of gasIDs using "compare_l2s_H16_H20.m" 
JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
if length(JOB) == 0
  JOB = 1;
  JOB = 72;
end

load('H16_20_24_gidlist.mat');
glist0 = gidlist.h24;
if JOB > length(glist0) + 1
  error('JOB > length(glist0) + 1')
elseif JOB == length(glist0) + 1
  disp('JOB == length(glist0) + 1     so will do all gases, wgt 1 and the filename will have 0')
end

if JOB <= length(glist0)
  gid = glist0(JOB);
  WGT  = 0.0;  
else
  gid = 25;
  WGT  = 1.0;
end  
fprintf(1,'processing JOB %5i == gasID removed %5i \n',JOB,gid);

CKDCKD = 32;
FF1    = 605;
FF2    = 2830;
XYZXYZ = use_this_rtp;
GGG    = gid;

outfilejunk = ['JUNK/individual_prof_convolved_kcartaH2020_H2024_' num2str(gid) '.mat'];
if JOB == length(glist0) + 1
  outfilejunk = ['JUNK/individual_prof_convolved_kcartaH2020_H2024_' num2str(0) '.mat'];
end
junkdir = dir(outfilejunk);

if length(junkdir) == 0 
  outnml   = ['nml_GGG_' num2str(GGG) '.nml'];
  rad20out = ['JUNK/radH20_' num2str(gid) '.mat'];
  rad24out = ['JUNK/radH24_' num2str(gid) '.mat'];  
  if JOB == length(glist0) + 1
    %% this is the control, all gases
    rad20out = ['JUNK/radH20_' num2str(0) '.mat'];
    rad24out = ['JUNK/radH24_' num2str(0) '.mat'];    
  end
  
  sedder = ['!sed ' ];
  sedder = [sedder ' -e "s/CKDCKD/'     num2str(CKDCKD) '/g"'];
  sedder = [sedder ' -e "s/FF1/'        num2str(FF1) '/g"'];
  sedder = [sedder ' -e "s/FF2/'        num2str(FF2) '/g"'];
  sedder = [sedder ' -e "s/MMM/'        num2str(MMM) '/g"'];
  
  sedder = [sedder ' -e "s/GGG/'        num2str(GGG) '/g"'];
  sedder = [sedder ' -e "s/WGT/'        num2str(WGT) '/g"'];
  
  sedder = [sedder ' -e "s/DOLBLRTM/'   num2str(2)   '/g"'];
  sedder = [sedder ' -e "s/XYZXYZ/'     XYZXYZ       '/g"'];      
  sedder = [sedder ' template_Qrad0.nml  > ' outnml];
  eval(sedder)
  
  kcartaer = ['!' kcartaexec20 ' ' outnml ' ' rad20out];
  eval(kcartaer)
  [rad,w] = readkcstd(rad20out);
  [fc,qc20] = quickconvolve(w,rad,0.25,0.25);

  kcartaer = ['!' kcartaexec24 ' ' outnml ' ' rad24out];
  eval(kcartaer)
  [rad,w] = readkcstd(rad24out);
  [fc,qc24] = quickconvolve(w,rad,0.25,0.25);

  rmer = ['!/bin/rm ' rad20out ' ' rad24out ' ' outnml];
  eval(rmer);
  
  saver = ['save ' outfilejunk ' fc qc20 qc24'];
  eval(saver);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now run loop_analyze_H2020_H2024.m')
