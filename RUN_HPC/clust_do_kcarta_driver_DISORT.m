function [] = clust_do_kcarta_driver_DISORT(iiBin)
%% % /bin/rm slurm* JUNK/rad.dat*; ; sbatch --array=G1-G2 sergio_matlab_jobB.sbatch

%% need to modify template_Qradcloud_2cloud_DISORT.nml CORRECTLY for the rtp file to process!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_rtp
set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm

disp('>>>>>')
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

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));   %% PERTAINS TO WAVENUMBER CHUNK, NOT PROFILE
%JOB = 46    %% can be between 1-89 for the 89 kCARTA chunks

if nargin == 0
  iiBin = 1;
end

f1 = 605 + (JOB-1)*25;
f2 = f1 + 25;
fprintf(1,'f1,f2 RESET TO %4i %4i \n',f1,f2);

fprintf(1,'processing JOB %5i == same profile %5i \n',JOB,iiBin);
iDISORT = +1;
do_kcarta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('NOW READ IN THE SEPARATE FREQ CHUNKS "read_disort_chunks" and then YOU call parts of do_convolve');
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(1,'kcarta exitcode = %3i  << iDoConvolve = %3i iInstr = %3i >> iDoRad = %3i \n',exitcode,iDoConvolve,iInstr,iDoRad);
%  do_convolve(iInstr,iiBin);
%  do_convolve_jac(iInstr,iiBin);

if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac <= 0) & exitcode == 0
  do_convolve(iInstr,iiBin,iDoRad);
%if iDoConvolve > 0 & iDoRad == 3 & exitcode == 0
%  do_convolve(iInstr,iiBin);
%end
elseif iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == 1 | iDoJac == 100) & exitcode == 0
  do_convolve(iInstr,iiBin);
  fprintf(1,'jacobian gasID gg = %4i iDoJac = %4i iDoCLoud = %4i \n',gg,iDoJac,iDoCloud)
  do_convolve_jac(gg,iInstr,iiBin,iDoJac,iDoCloud);
end

