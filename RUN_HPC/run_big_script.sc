## instead of singleton paul suggest jobs 1-100 doing a 1-100%16 will limit it to just 16 jobs going at once.
## so eg sbatch --array=1-588%64 sergio_matlab_jobB.sbatch

/bin/rm slurm*.out JUNK/rad.dat* 
sbatch --array=1-50   sergio_matlab_jobB.sbatch 
sbatch --array=51-100 -d singleton sergio_matlab_jobB.sbatch 
sbatch --array=101-150 -d singleton sergio_matlab_jobB.sbatch 
sbatch --array=151-200 -d singleton sergio_matlab_jobB.sbatch 
sbatch --array=201-250 -d singleton sergio_matlab_jobB.sbatch
sbatch --array=251-300 -d singleton sergio_matlab_jobB.sbatch 
sbatch --array=301-350 -d singleton sergio_matlab_jobB.sbatch
sbatch --array=351-400 -d singleton sergio_matlab_jobB.sbatch 
sbatch --array=401-450 -d singleton sergio_matlab_jobB.sbatch
sbatch --array=451-500 -d singleton sergio_matlab_jobB.sbatch
sbatch --array=501-550 -d singleton sergio_matlab_jobB.sbatch
sbatch --array=551-600 -d singleton sergio_matlab_jobB.sbatch
