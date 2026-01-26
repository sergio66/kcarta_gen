#!/bin/bash

# after launchung job watch squeue -u sergio
#
# or check_progress.sc
#
# run this with   sbatch -p 2018 --array=1-49 --output='testslurm' sergio_matlab_chip.sbatch 1
# run this with   sbatch -p 2018 --array=1-49 sergio_matlab_chip.sbatch 1
# run this with   sbatch --exclude=cnode203,cnode267 --array=1-49 sergio_matlab_chip.sbatch 1
#
# run this with   /bin/rm batchout*.* slurm* JUNK/individual_prof_convolved_kcarta* JUNK/rad.dat* JUNK/jac.dat*; sbatch -p 2018 --array=1-49 sergio_matlab_chip.sbatch 1
# run this with   /bin/rm batchout*.* slurm* JUNK/individual_prof_convolved_kcarta* JUNK/rad.dat*; sbatch -p 2018 --array=1-Nmax%128 sergio_matlab_chip.sbatch 1
# run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p 2018 --array=1-Nmax%128 sergio_matlab_chip.sbatch 1
# run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p 2018 --exclude=cnode203,cnode267 --array=1-Nmax%128 sergio_matlab_chip.sbatch 1
# run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p 2021                              --array=1-Nmax%128 sergio_matlab_chip.sbatch 1

########################################################################
## for DISORT : submit_many_disort_runs.m      can launch jobs for  sbatch -p 2021  --array=1111 sergio_matlab_chip.sbatch 9 profileles A to B, does whole 605-2830 cm-1
#   run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;    sbatch -p 2021  --array=1111 sergio_matlab_chip.sbatch 9 789   runs DISORT all chunks, profile 789
#                                           yes 11111 is there cos JOB is not needed, wierd
#                                           submit_many_disort_runs.m   can launch multiple profiles then gather up results using read_disort_chunks
## can do the same thing with this ...!!!! so I've kinda retired option 9 (clust_do_kcarta_driver_DISORT.m)
#   run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p 2021      --array=1-Nmax%128 sergio_matlab_chip.sbatch 10      
#     to run DISORT for 89 chunks,straighforward readkcstd for all profiles
## can put ONE profile (eg 789) onto 89 prpocessors, then use read_disort_chunks.m to gather
#   run this with   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc; sbatch -p 2021      --array=1-89 sergio_matlab_chip.sbatch 12 789
#   this needs a call to read_disort_chunks after it runs
## 
## so to run just ONE profile PQR can do one or the other of these
## rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;
##   sbatch -p 2021  --array=1111         sergio_matlab_chip.sbatch 9 PQR      ONE profile, VEDY VEDY SLOW, one processor
##   sbatch -p 2021  --array=PQR          sergio_matlab_chip.sbatch 10         ONE profile, VEDY VEDY SLOW, one processor
##   sbatch -p 2021  --array=1-89         sergio_matlab_chip.sbatch 11 PQR     <<< ONE profile, VERY FAST, 89 processors >>>
##
## or to run multiple profiles, different chunks 1-89 
#   this needs a call to read_disort_chunks after it runs
#   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;  sbatch -p 2021  --array=1-Nprofs/iaChunksize sergio_matlab_chip.sbatch 12 78   << make sure you do Nprofs/iaChunkSize!!!! 
#                                                                       << run DISORT for freqchunk 78, profiles (1:iaChunkSize) + (JOB-1)*iaChunkSize
#                                        submit_many_disort_run_loops.m   can launch multiple profiles then gather up results using read_disort_chunks
# >>>>>>>>>>>>>>>>>>>>>>>>> ########################################################################
# if putting ONE job per node (ie no loop)
#    run this with   /bin/rm batchout*.out batchout*.err JUNK/rad.dat*; sbatch --array=1-588%128 sergio_matlab_chip.sbatch 1
#    where the %128 limits to 128 jobs at a go, to limit simultaneous resource need
# if putting MANY jobs per node (ie loop)
#   run this with   /bin/rm batchout*.out batchout*.err JUNK/rad.dat*; sbatch --array=1-Nprofs/iaChunksize sergio_matlab_chip.sbatch 2 << make sure you do Nprofs/iaChunkSize!!!!
#                   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;  sbatch --partition=2018   --array=1-Nprofs/iaChunksize sergio_matlab_chip.sbatch 2       << make sure you do Nprofs/iaChunkSize!!!!
#                   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;  sbatch --partition=2021  --array=1-Nprofs/iaChunksize sergio_matlab_chip.sbatch 2       << make sure you do Nprofs/iaChunkSize!!!!
#                   rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;  sbatch -p 2021           --array=1-Nprofs/iaChunksize sergio_matlab_chip.sbatch 12 78   << make sure you do Nprofs/iaChunkSize!!!! 
#                                                                                                << run DISORT for freqchunk 78, profiles (1:iaChunkSize) + (JOB-1)*iaChunkSize
#                                        submit_many_disort_run_loops.m   can launch multiple profiles then gather up results using read_disort_chunks
# >>>>>>>>>>>>>>>>>>>>>>>>> ########################################################################

#  Name of the job:
#SBATCH --job-name=KCARTA_DRIVER

# requeue; if held type scontrol release JOBID
#SBATCH --requeue

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#SBATCH --account=pi_strow

#SBATCH --cluster=chip-cpu

########################################################################

#  Specify the qos and run time (format:  dd-hh:mm:ss)

## 100 layer clouds and TwoSlab clouds jacs rads, loop 
##SBATCH --qos=medium
##SBATCH --time=3:59:00

# DISORT rads, using option 10 where we loop
##SBATCH --qos=medium
##SBATCH --time=1:59:00

# clear, or TwoSlab clouds, or fluxes
##SBATCH --qos=short
##SBATCH --time=0:59:00

#########################

#  Specify the desired partition develop/batch/prod
## can shorthand this as -p high_mem on the command line
##SBATCH --partition=batch
##SBATCH --partition=2018
##SBATCH --partition=2021
##SBATCH --partition=2024

#### to make sure we have exclusive privilege, two slab jacs
##SBATCH --qos=pi_strow
##SBATCH --partition=pi_strow
##SBATCH --time=3:59:00

#### we run on 2024 but run risk of someone pre-empting, clear sky
#SBATCH --qos=shared
#SBATCH --partition=2024
#SBATCH --time=0:59:00
########################################################################

##  This is in MB, very aggressive but I have been running outta memory
##SBATCH --mem-per-cpu=24000
##  This is in MB, less aggressive
##SBATCH --mem-per-cpu=12000
##  This is in MB, very lean, don;t go below this as java starts crying, good for almost everything
##SBATCH --mem-per-cpu=4000
##  This is in MB, less aggressive
##SBATCH --mem-per-cpu=8000
##  This is in MB, less aggressive
#SBATCH --mem-per-cpu=16000

########################################################################

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

##  Specify the job array (format:  start-stop:step)
## DO NOT USE srun matlab -nodisplay -r "clust_do_kcarta_driver; exit"   but instead matlab -nodisplay -r "clust_do_kcarta_driver; exit" 

if [[ "$1" -eq "" ]]; then
  # this is individual clear/cld runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is individual clear/cld runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver; exit"
elif [[ "$1" -eq "2" ]]; then
  # this is for looped clear runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "iaChunkSize = 30; clust_do_kcarta_driver_loop; exit"
elif [[ "$1" -eq "3" ]]; then
  # this is for looped allsky runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "iaChunkSize = 15; clust_do_kcarta_driver_loop; exit"
elif [[ "$1" -eq "4" ]]; then
  # this is for clear/cloud filelist
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_filelist; exit"
elif [[ "$1" -eq "5" ]]; then
  # this is for clear/cloud looped filelist
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "iaChunkSize = 15; clust_do_kcarta_driver_filelist_loop; exit"
elif [[ "$1" -eq "6" ]]; then
  # this is for all bands
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "cluster_loop_allkcartabands; exit"
elif [[ "$1" -eq "7" ]]; then
  # this is for all profiles inside a rtp file, but can do many rtp files eg
  # use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin32/all12monthavg_T_WV_grid_latbin_32_lonbin_12.rtp'; %% 25 T x WV grids
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "cluster_loop_allprofiles_onefile; exit"
elif [[ "$1" -eq "8" ]]; then
  # this loops over H2008,H2012,H2016
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_allHITRAN; exit"

######################### DISORT #########################
######################### DISORT #########################
######################### DISORT #########################

elif [[ "$1" -eq "9" ]]; then
  # this does entire spectrum
  # look at submit_many_disort_runs.m  --- this is RETIRED --- use "clust_do_kcarta_driver_DISORT_wholespectrum" instead of this and/or "submit_many_disort_runs.m"
  if [[ $# -eq 1 ]]; then
    ## this does profile number 1
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_DISORT(1); exit"
  elif [[ $# -eq 2 ]]; then
    ## this does profile number ($2)
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_DISORT($2); exit"
  fi
elif [[ "$1" -eq "10" ]]; then
  # this is individual entire spectrum DISORT runs, 90+ minutes each, this is slow but straightforward
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_DISORT_wholespectrum; exit"
elif [[ "$1" -eq "11" ]]; then
  ## this needs a call to read_disort_chunks after it runs
  ## this needs a call to read_disort_chunks after it runs
  ## this needs a call to read_disort_chunks after it runs
  # this loops over 25 cm-1 chunks for DISORT, if you only have a handful of profiles this is fast since it breaks profiles/freqchunks for separate processors
  if [[ $# -eq 1 ]]; then
    ## this does profile number 1
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_DISORT_individualprofile_clusterchunks(1); exit"
  elif [[ $# -eq 2 ]]; then
    ## this does profile number ($2)
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_do_kcarta_driver_DISORT_individualprofile_clusterchunks($2); exit"
  fi
elif [[ "$1" -eq "12" ]]; then
  ## this needs a call to read_disort_chunks after it runs
  ## this needs a call to read_disort_chunks after it runs
  ## this needs a call to read_disort_chunks after it runs
  # this loops over 25 cm-1 chunks for DISORT, this is complicated
  if [[ $# -eq 1 ]]; then
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "iaChunkSize = 50; clust_do_kcarta_driver_DISORT_loop(1,iaChunkSize); exit"
  elif [[ $# -eq 2 ]]; then
    /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "iaChunkSize = 50; clust_do_kcarta_driver_DISORT_loop($2,iaChunkSize); exit"
  fi
fi
