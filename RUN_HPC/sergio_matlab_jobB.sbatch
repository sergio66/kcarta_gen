#!/bin/bash

# after launchung job watch squeue -u sergio
#
# or check_progress.sc
#
# run this with   sbatch -p high_mem --array=1-49 --output='testslurm' sergio_matlab_jobB.sbatch 1
# run this with   sbatch -p high_mem --array=1-49 sergio_matlab_jobB.sbatch 1
# run this with   sbatch --exclude=cnode203,cnode267 --array=1-49 sergio_matlab_jobB.sbatch 1
#
# run this with   /bin/rm batchout*.* slurm* JUNK/individual_prof_convolved_kcarta* JUNK/rad.dat* JUNK/jac.dat*; sbatch -p high_mem --array=1-49 sergio_matlab_jobB.sbatch 1
# run this with   /bin/rm batchout*.* slurm* JUNK/individual_prof_convolved_kcarta* JUNK/rad.dat*; sbatch -p high_mem --array=1-Nmax%128 sergio_matlab_jobB.sbatch 1
# run this with   rmer_slurm_JUNKrad.sc; sbatch -p high_mem --array=1-Nmax%128 sergio_matlab_jobB.sbatch 1
# run this with   rmer_slurm_JUNKrad.sc; sbatch -p high_mem --exclude=cnode203,cnode267 --array=1-Nmax%128 sergio_matlab_jobB.sbatch 1
# run this with   rmer_slurm_JUNKrad.sc; sbatch -p cpu2021                              --array=1-Nmax%128 sergio_matlab_jobB.sbatch 1
#
# >>>>>>>>>>>>>>>>>>>>>>>>>
# if putting ONE job per node (ie no loop)
# run this with   /bin/rm batchout*.out batchout*.err JUNK/rad.dat*; sbatch --array=1-588%128 sergio_matlab_jobB.sbatch 1
#    where the %128 limits to 128 jobs at a go, to limit simultaneous resource need
# if putting MANY jobs per node (ie loop)
# run this with   /bin/rm batchout*.out batchout*.err JUNK/rad.dat*; sbatch --array=1-Nprofs/iaChunksize sergio_matlab_jobB.sbatch 2 << make sure you do Nprofs/iaChunkSize!!!!
#                 rmer_slurm_JUNKrad.sc;  sbatch --partition=high_mem --array=1-Nprofs/iaChunksize sergio_matlab_jobB.sbatch 2 << make sure you do Nprofs/iaChunkSize!!!!
#                 rmer_slurm_JUNKrad.sc;  sbatch --partition=cpu2021  --array=1-Nprofs/iaChunksize sergio_matlab_jobB.sbatch 2 << make sure you do Nprofs/iaChunkSize!!!!
# >>>>>>>>>>>>>>>>>>>>>>>>>

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

#  Specify the desired partition develop/batch/prod
##SBATCH --partition=batch
## can shorthand this as -p high_mem on the command line
#SBATCH --partition=high_mem

########################################################################
#  Specify the qos and run time (format:  dd-hh:mm:ss)

## 100 layer clouds
##SBATCH --qos=medium+
##SBATCH --time=3:59:00

# clear, or TwoSlab clouds, or fluxes
##SBATCH --qos=short+
##SBATCH --time=0:59:00

# TwoSlab clouds jacs rads, loop 
#SBATCH --qos=medium+
#SBATCH --time=3:59:00

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

## defualt error output is to (by default slurm-{job id}_{array #}.out)
##  here we help spread the pain, into a directory
##SBATCH -o /home/sbuczko1/logs/sbatch/run_airibrad_rand_rtp-%A_%a.out
##SBATCH -e /home/sbuczko1/logs/sbatch/run_airibrad_rand_rtp-%A_%a.err
#######  here we help spread the pain, or into the dir where the call was initiated ... do this for fluu error reporting
##SBATCH -o batchout-%A_%a.out
##SBATCH -e batchout-%A_%a.err

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
  matlab -nodisplay -r "clust_do_kcarta_driver; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is individual clear/cld runs
  matlab -nodisplay -r "clust_do_kcarta_driver; exit"
elif [[ "$1" -eq "2" ]]; then
  # this is for looped clear runs
  matlab -nodisplay -r "iaChunkSize = 30; clust_do_kcarta_driver_loop; exit"
elif [[ "$1" -eq "3" ]]; then
  # this is for looped allsky runs
  matlab -nodisplay -r "iaChunkSize = 15; clust_do_kcarta_driver_loop; exit"
elif [[ "$1" -eq "4" ]]; then
  # this is for clear/cloud filelist
  matlab -nodisplay -r "clust_do_kcarta_driver_filelist; exit"
elif [[ "$1" -eq "5" ]]; then
  # this is for clear/cloud looped filelist
  matlab -nodisplay -r "iaChunkSize = 15; clust_do_kcarta_driver_filelist_loop; exit"
elif [[ "$1" -eq "6" ]]; then
  # this is for all bands
  matlab -nodisplay -r "cluster_loop_allkcartabands; exit"
elif [[ "$1" -eq "7" ]]; then
  # this is for all profiles inside a rtp file, but can do many rtp files eg
  # use_this_rtp = '/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/16dayAvgLatBin32/all12monthavg_T_WV_grid_latbin_32_lonbin_12.rtp'; %% 25 T x WV grids
  matlab -nodisplay -r "cluster_loop_allprofiles_onefile; exit"
elif [[ "$1" -eq "8" ]]; then
  # this loops over H2008,H2012,H2016
  matlab -nodisplay -r "clust_do_kcarta_driver_allHITRAN; exit"
fi
