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
disp('>>>>>')

use_this_rtp00  = use_this_rtp;
strIceCloud00   = strIceCloud;
strWaterCloud00 = strWaterCloud;

iInstr0 = 2; %% IASI
iInstr0 = 1; %% AIRS

iaH = [2008 2012 2015 2016];
iaH = [2016.3 2017];  %% 2017 is high res w 3 gases, while 2016.3 is usual res with 3 gases

iaH = [2017];         %% 2017 is high res w 3 gases, while 2016.3 is usual res with 3 gases
iaH = [2016.3];       %% 2016.3 is usual res with 3 gases, while 2017 is high res w 3 gases

for iiLoop = 1 : length(iaH)

  use_this_rtp  = use_this_rtp00;
  strIceCloud   = strIceCloud00;
  strWaterCloud = strWaterCloud00;

  iHITRAN = iaH(iiLoop);
  set_kcarta_exec_iHITRAN

  fprintf(1,'iHITRAN = %4i kcartaexec   = %s \n',iHITRAN,kcartaexec);

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
  
  JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
  % JOB = 49;
  % JOB = 1;
  % JOB = 269
  % JOB = 705;
  % JOB = 16
  % JOB = 2000
  % JOB = 218

  iiBin = JOB;
  rmer = ['!/bin/rm JUNK/rad.dat' num2str(iiBin) ' JUNK/jac.dat' num2str(iiBin) ' JUNK/jac.dat' num2str(iiBin) '_COL'];
  eval(rmer)

  fprintf(1,'processing JOB %5i == same profile %5i \n',JOB,iiBin);

  iInstr = iInstr0; 

  if iInstr == 1
    fout = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
  elseif iInstr == 2
    fout = ['JUNK/individual_prof_convolved_kcarta_iasi_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
  elseif iInstr == 4 
    fout = ['JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
  elseif iInstr == 14 | iInstr == 24 | iInstr == 124
    fout = ['JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
  end

%fout
%exist(fout,'file')
%error('kjgs')

  if ~exist(fout,'file')
    %iDo_rt_1vs43 = -1;  %% this does not make a difference unless you set in template_reset_nm_spectra.nml
    %iDo_rt_1vs43 = +43; %% this does not make a difference unless you set in template_reset_nm_spectra.nml
    do_kcarta
  
    iInstr = iInstr0;  
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
  
    if iInstr == 1
      mver = ['!/bin/mv JUNK/individual_prof_convolved_kcarta_airs_' num2str(iiBin) '.mat'];
      mver = [mver ' JUNK/individual_prof_convolved_kcarta_airs_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
    elseif iInstr == 2
      mver = ['!/bin/mv JUNK/individual_prof_convolved_kcarta_iasi_' num2str(iiBin) '.mat'];
      mver = [mver ' JUNK/individual_prof_convolved_kcarta_iasi_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
    elseif iInstr == 4 
      mver = ['!/bin/mv JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(iiBin) '.mat'];
      mver = [mver ' JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
    elseif iInstr == 14 | iInstr == 24 | iInstr == 124
      mver = ['!/bin/mv JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(iiBin) '.mat'];
      mver = [mver ' JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(iiBin) '_HITorGEI' num2str(iHITRAN) '.mat'];
    end
    eval(mver)
  
  else
    fprintf(1,'%s already exists \n',fout);
  end
  
  %disp('ret to continue'); pause
end
  
