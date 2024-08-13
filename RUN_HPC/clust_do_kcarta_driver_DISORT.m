function [] = clust_do_kcarta_driver_DISORT(iProf)
%% % /bin/rm slurm* JUNK/rad.dat*; ; sbatch --array=G1-G2 sergio_matlab_jobB.sbatch 09

%% need to modify template_Qradcloud_2cloud_DISORT.nml CORRECTLY for the rtp file to process!
%% read in the individual chunks and put together entire spectrum using read_disort_chunks

%% this is basically same as clust_do_kcarta_driver_DISORT_wholespectrum and has been retired
%% this is basically same as clust_do_kcarta_driver_DISORT_wholespectrum and has been retired
%% this is basically same as clust_do_kcarta_driver_DISORT_wholespectrum and has been retired

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_rtp
set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm

disp('>>>>>')
disp(' clust_do_kcarta_driver_DISORT')
fprintf(1,'kcartaexec   = %s \n',kcartaexec);
fprintf(1,'f1,f2 WILL BE RESET = %4i %4i \n',f1,f2);
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

disp('DISORT runs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
use_this_rtp0 = use_this_rtp;
set_convolver
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

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% PERTAINS TO NOTHING
if length(JOB) == 0
  JOB = 12;  %% do the 905-930 cm-1 chunk
  JOB = -99999999;
end
JOB = -99999999;
disp('JOB is irrelevant in this clust_do_kcarta_driver_DISORT.m : does either 605-2830 cm-1 or 310-510 cm-1') 
disp('JOB is irrelevant in this clust_do_kcarta_driver_DISORT.m : does either 605-2830 cm-1 or 310-510 cm-1') 
disp('JOB is irrelevant in this clust_do_kcarta_driver_DISORT.m : does either 605-2830 cm-1 or 310-510 cm-1') 

disp('try clust_do_kcarta_driver_DISORT_individualprofile_clusterchunks.m ... you ive it a profile and JOB is the chunk number')
disp('try clust_do_kcarta_driver_DISORT_individualprofile_clusterchunks.m ... you ive it a profile and JOB is the chunk number')
disp('try clust_do_kcarta_driver_DISORT_individualprofile_clusterchunks.m ... you ive it a profile and JOB is the chunk number')
%error('this is kinda retired ... use clust_do_kcarta_driver_DISORT_wholespectrum.m since there JOB ======= iProf')

if nargin == 0
  disp('no arguments, iProf == iProf = 1')
  iProf = 1;  
end

iIRorFIR = -1;
iIRorFIR = +1;
if iIRorFIR == +1
  %% 89 chunks
  f1 = 605;
  f2 = 2830;
  iKCKD = iKCKD;
elseif iIRorFIR == -1
  %% 20 chunks
  f1 = 310;
  f2 = 510;
  iKCKD = 1;
end
fprintf(1,'f1,f2 RESET TO %4i %4i \n',f1,f2);

%%% fprintf(1,'processing kCARTA freq chunk JOB %5i profile %5i iDoRad = %3i \n',JOB,iProf,iDoRad);
iiBin = iProf;   %% THIS IS THE BIN OR PROFILE in the RTPFIL, used in sed_1_2_cloudfiles.m
fprintf(1,'processing kCARTA freq chunk %5i-%5i profile %5i iDoRad = %3i \n',f1,f2,iProf,iDoRad);
iDISORT = +1;
iDISopt = 09;

do_kcarta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp('NOW READ IN THE SEPARATE FREQ CHUNKS "read_disort_chunks" and then YOU call parts of do_convolve');
% outnameCLD = [outname '_CLD'];
% rmerCLD = ['!/bin/rm ' outnameCLD];
% eval(rmerCLD);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'kcarta exitcode = %3i  << iDoConvolve = %3i iInstr = %3i >> iDoRad = %3i \n',exitcode,iDoConvolve,iInstr,iDoRad);
%  do_convolve(iInstr,iiBin);
%  do_convolve_jac(iInstr,iiBin);

if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == -1) & exitcode == 0
  do_convolve(iInstr,iiBin,iDoRad);
%if iDoConvolve > 0 & iDoRad == 3 & exitcode == 0
%  do_convolve(iInstr,iiBin);
%end
elseif iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == 1 | abs(iDoJac) == 100) & exitcode == 0
  do_convolve(iInstr,iiBin);
  fprintf(1,'jacobian gasID gg = %4i iDoJac = %4i iDoCLoud = %4i \n',gg,iDoJac,iDoCloud)
  do_convolve_jac(gg,iInstr,iiBin,iDoJac,iDoCloud);
end

