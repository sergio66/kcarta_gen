disp('make sure you have run      rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;      before')
disp('make sure you have run      rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;      before')
disp('make sure you have run      rmer_slurm_JUNKconv.sc; rmer_slurm_JUNKrad.sc;      before')

JMAX = 660;

echoer = ['!echo $PATH'];            eval(echoer)
echoer = ['!echo $LD_LIBRARY_PATH']; eval(echoer)

for JOBB = 1 : JMAX
  clear JOB
  clust_do_kcarta_driver
end  
