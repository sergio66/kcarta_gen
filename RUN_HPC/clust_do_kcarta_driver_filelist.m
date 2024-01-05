%% % sbatch --array=G1-G2 sergio_matlab_jobB.sbatch

%% need to modify template_Qrad.nml CORRECTLY for the rtp file to process!

addpath /home/sergio/MATLABCODE
system_slurm_stats

%%%%%%%%%%%%%%%%%%%%%%%%%

list_of_profiles_to_process = 'JUNK/notdone.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%

set_rtp
set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm
kcartaexec

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

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
thelist = load(list_of_profiles_to_process);
iiBin = thelist(JOB);
fprintf(1,'processing JOB %5i --> profile %5i \n',JOB,iiBin);
do_kcarta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% if iDoConvolve > 0 & (iDoRad == 3 | iDoRad == 10) & exitcode == 0
if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == -1) & exitcode == 0
  do_convolve(iInstr,iiBin);
end
if iDoConvolve > 0 & (iDoRad == 3 | iDoRad == 10) & (iDoJac == 1 | abs(iDoJac) == 100) & exitcode == 0
  fprintf(1,'jacobian gasID gg = %4i iDoJac = %4i iDoCLoud = %4i \n',gg,iDoJac,iDoCloud)
  do_convolve_jac(gg,iInstr,iiBin,iDoJac,iDoCloud);
end
