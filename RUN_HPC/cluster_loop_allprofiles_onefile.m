%% % /bin/rm slurm* JUNK/rad.dat*; ; sbatch --array=G1-G2 sergio_matlab_jobB.sbatch

%% need to modify template_QXYZ.nml CORRECTLY for the rtp file to process!

%% this is for the rtp files served by eg
%%  /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_starts_Sept2002/call_save_split_apart_rtp_howard_bins_startSept2002_16daysavgs.m
%% and then checked by eg
%%  /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/TILES_TILES_TILES_MakeAvgCldProfs2002_2020/Code_starts_Sept2002/driver_jacs_T_WV_grid_allmonths_VERSUS_allsteps.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
% JOB = 49;
% JOB = 1;
% JOB = 269
% JOB = 705;
% JOB = 16
% JOB = 2000

theX = mod(JOB,72);             theXstr = num2str(theX,'%02d');
theY = floor((JOB-0.1)/72) + 1; theYstr = num2str(theY,'%02d');

%{
sedder for namelist = !sed -e "s/FF1/605/g"  -e "s/FF2/2830/g"  -e "s/MMM/124/g" -e "s/TTT/-1/g" -e "s/XYZXYZ/\/home\/sergio\/KCARTA\/WORK\/RUN_TARA\/GENERIC_RADSnJACS_MANYPROFILES\/RTP\/armtwp_nearby_semiclear_jan04_resetT_nte.rt
p/g" -e "s/CKDCKD/32/g" -e "s/DOLBLRTM/2/g" -e "s/COMMENTCOMMENT/04-Jul-2021 \/home\/sergio\/KCARTA\/BIN\/kcarta.x_f90_122_385ppmv_H08_CO2_UMBC iDoRad=3 iDoLBLRTM=2 iDo_rt_1vs43=-1 iDoC/g" template_Qrad.nml  > run_nml124
kcartaer = !time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_385ppmv_H08_CO2_UMBC run_nml124 JUNK/rad.dat124; echo $? >& status124
%}

%% 25 T x WV grids
use_this_rtp = ['/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin' theYstr '/all12monthavg_T_WV_grid_latbin_' theYstr '_lonbin_' theXstr '.rtp']; 

dirout = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Oct2020_startSept2002_trendsonly/AllDemJacsAnomaly/T_WV_Grid/'];
  dirout = [dirout theYstr '/' theXstr '/'];
  if ~exist(dirout)
    mker = ['!mkdir -p ' dirout];
    eval(mker);
  end

loop_allprofiles_onefile = +1;

for ttBin = 1 : 25   %% number of timesteps in file
  use_this_rtp = ['/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin' theYstr '/all12monthavg_T_WV_grid_latbin_' theYstr '_lonbin_' theXstr '.rtp']; %% else the / keeps growing and the \/ keeps growing

  fxstat   = [dirout '/status' num2str(ttBin)];
  fxnml    = [dirout '/run_nml' num2str(ttBin)];
  fxrad    = [dirout '/rad.dat' num2str(ttBin)];
  fxjac    = [dirout '/jac.dat' num2str(ttBin)];
  fxcoljac = [dirout '/jac.dat' num2str(ttBin) '_COL'];

  fxconvrad = [dirout '/individual_prof_convolved_kcarta_airs_' num2str(ttBin) '.mat'];
  fxconvjac = [dirout '/individual_prof_convolved_kcarta_airs_' num2str(ttBin) '_jac.mat'];
  fxconvcol = [dirout '/individual_prof_convolved_kcarta_airs_' num2str(ttBin) '_coljac.mat'];

  if ~exist(fxconvjac)
    rmerx = ['!/bin/rm ' fxnml ' ' fxrad ' ' fxjac ' ' fxcoljac ' ' fxstat]; eval(rmerx);
    
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
        
    iiBin = JOB;
    iiBin = -1;
    iiBin = ttBin;
    fprintf(1,'processing JOB %5i ==> XX = %2i YY = %2i ttBin = %3i \n',JOB,theX,theY,ttBin);
    do_kcarta
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf(1,'kcarta exitcode = %3i  << iDoConvolve = %3i iInstr = %3i >> iDoRad = %3i \n',exitcode,iDoConvolve,iInstr,iDoRad);
    %  do_convolve(iInstr,iiBin);
    %  do_convolve_jac(iInstr,iiBin);
    
    if iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac <= 0) & exitcode == 0
      do_convolve(iInstr,iiBin,iDoRad,dirout);
    %if iDoConvolve > 0 & iDoRad == 3 & exitcode == 0
    %  do_convolve(iInstr,iiBin);
    %end
    elseif iDoConvolve > 0 & (iDoRad == 0 | iDoRad == 3 | iDoRad == 10) & (iDoJac == 1 | iDoJac == 100) & exitcode == 0
      do_convolve(iInstr,iiBin,iDoRad,dirout);
      fprintf(1,'jacobian gasID gg = %4i iDoJac = %4i iDoCLoud = %4i \n',gg,iDoJac,iDoCloud)
      do_convolve_jac(gg,iInstr,iiBin,iDoJac,iDoCloud,dirout);
    end
  
    rmerx = ['!/bin/rm ' fxnml ' ' fxrad ' ' fxjac ' ' fxcoljac ' ' fxstat]; 
    rmerx = ['!/bin/rm '       ' ' fxrad ' ' fxjac ' ' fxcoljac ' ' fxstat]; 
    eval(rmerx);
    disp(' ')
    disp('<<<<<<<<<<<<<<<<<<<<<    >>>>>>>>>>>>>>>>>>>>>>> ')
    disp(' ')

  %disp('ret')
  %pause
  end   %% if ~exist(fxjac)
end     %% for ttBin
