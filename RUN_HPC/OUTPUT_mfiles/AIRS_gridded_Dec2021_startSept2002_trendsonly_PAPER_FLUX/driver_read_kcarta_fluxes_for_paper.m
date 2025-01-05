%% this is /home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m

%{
cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES
  set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m : iKCKD =  32; iHITRAN = 2020; iDoLBLRTM = 2; iDoRad = 3;  iDoCloud = -1; iDoJac = -1;               %% clrsky, use LBLRTM ODs  ****************************** BEST DEFAULT clrsky         no jacs
  
  rad0 : 
    set_rtp.m : use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';            %% clear, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
    rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p cpu2021                              --array=1-4608 sergio_matlab_jobB.sbatch 6
    watch 'ls -lt /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof* | wc -l'

    when done
      mv JUNK/allkcbands_prof* JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Raw/.

  pert : 
    set_rtp.m : use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp';            %% clear, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
    rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p cpu2021                              --array=1-4608 sergio_matlab_jobB.sbatch 6
    watch 'ls -lt /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof* | wc -l'

    when done
      mv JUNK/allkcbands_prof* JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Pert/.

%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_rtp.m
%% use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp';            %% clear, raw,  see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/read_fileMean17years.m
%% use_this_rtp = 'RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp';       %% clear, pert, see /home/sergio/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m
%%%

iMakePertProf = -1;
if iMakePertProf > 0
  pert_kcarta_fluxes_for_paper
  rtpwrite('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5_pert.rp.rtp',h,ha,pert,pa);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iRawOrPert = +1;
iaFound = [];
for ii = 1 : 4608
  fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/allkcbands_prof' num2str(ii) '.mat'];
  if iRawOrPert > 0
    fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Raw/allkcbands_prof' num2str(ii) '.mat'];
  elseif iRawOrPert < 0
    fin = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/Pert/allkcbands_prof' num2str(ii) '.mat'];
  end
  if exist(fin)
    iaFound(ii) = +1;
    loader = ['a = load(''' fin ''');'];
    eval(loader);
    allrads.dall(:,ii) = a.dall;
  else
    iaFound(ii) = 0;
  end
end
allrads.wall = a.wall;

moo = find(iaFound > 0);
fprintf(1,'found %4i; max = %4i \n',sum(iaFound),max(moo))

moo = find(allrads.wall >= 605.00-1e-05);
[fc,qc] = convolve_airs(allrads.wall(moo),allrads.dall(moo,:),1:2834); 
moo = find(allrads.wall <= fc(1));
[fgc,qgc] = quickconvolve(allrads.wall(moo),allrads.dall(moo,:),0.5,0.5); 

oo = find(iaFound==1); oo = oo(end); 
plot(a.wall,a.dall,[fgc; fc],[qgc(:,oo); qc(:,oo)])

comment = 'see /home/MATLABCODE/oem_pkg_run/AIRS_gridded_STM_May2021_trendsonlyCLR/driver_read_kcarta_fluxes_for_paper.m'
if iRawOrPert > 0
  save /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/raw.mat fc qc fgc qgc comment
else
  save /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly_PAPER_FLUX/pert.mat fc qc fgc qgc comment
end
