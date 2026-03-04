#!/bin/bash

# after launchung job watch squeue -u sergio
#
# or check_progress.sc
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# run this with   /bin/rm slurm*; sbatch --array=1-72 --output='/dev/null/' sergio_matlab_chip.sbatch 1
# run this with   /bin/rm slurm*; sbatch --array=1-72                       sergio_matlab_chip.sbatch 1
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
# >>>>>>>>>>>>>>>>>>>>>>>>> ########################################################################

#  Name of the job:
#SBATCH --job-name=CHANGED_GASES

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

#### to make sure we have exclusive privilege, two slab jacs
##SBATCH --qos=pi_strow
##SBATCH --partition=pi_strow
##SBATCH --time=0:59:00

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
  # this is individual clear runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_run_kcarta_H2020_H2024; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is individual clear runs
  /usr/ebuild/installs/software/MATLAB/2023b/bin/matlab -nodisplay -r "clust_run_kcarta_H2020_H2024; exit"    
fi
