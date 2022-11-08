%% % /bin/rm slurm* JUNK/rad.dat*; ; sbatch --array=G1-G2 sergio_matlab_jobB.sbatch

%% need to modify template_QXYZ.nml CORRECTLY for the rtp file to process!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set_rtp
set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm
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
fprintf(1,'<<< iaChunkSize  = %4i >>> \n',iaChunkSize);
disp('>>>>>')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT TOUCH THESE LAST TWO LINES. EDIT set_convolver as needed
use_this_rtp0 = use_this_rtp;
set_convolver
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
%JOB = 49;
%JOB = 2222;
%JOB = 308

% iaChunkSize = 30; %% assume 30 jobs will take 30 x 1.5 = 45 mins, set in sbatch script
% obviously if iChunkSize == 1 then each node processes only ONE profile

iiBinAll = (1:iaChunkSize) + (JOB-1)*iaChunkSize;
for ix = 1 : length(iiBinAll)
  iiBin = iiBinAll(ix);
  fprintf(1,'JOB = %4i index %2i of %2i : processing %5i \n',JOB,ix,length(iiBinAll),iiBin);
  outfilejunk = ['JUNK/individual_prof_convolved_kcarta_crisHI_' num2str(iiBin) '.mat'];
  outfilejunk = ['JUNK/individual_prof_convolved_kcarta_*_' num2str(iiBin) '.mat'];
  junkdir = dir(outfilejunk);

  %if ~exist(outfilejunk)
  if length(junkdir) == 0
    do_kcarta

    fprintf(1,'done kcarta run ... iDoConvolve = %2i iDoRad = %2i iDoJac = %3i \n',iDoConvolve,iDoRad,iDoJac)
    if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & exitcode == 0
      do_convolve(iInstr,iiBin);
    end
    if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == 1 | iDoJac == 100) & exitcode == 0
      fprintf(1,'jacobian gasID gg = %4i iDoJac = %4i iDoCLoud = %4i \n',gg,iDoJac,iDoCloud)
      do_convolve_jac(gg,iInstr,iiBin,iDoJac,iDoCloud);
    end

  else
    fprintf(1,'%s already exists, next!!! \n',outfilejunk)
  end

end
