addpath /asl/matlib/h4tools

disp('make sure you run rmer_slurm_JUNKrad.sc before hand')
disp('make sure you run rmer_slurm_JUNKrad.sc before hand')
disp('make sure you run rmer_slurm_JUNKrad.sc before hand')

disp('ret to continue'); pause

set_rtp
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %4i profiles \n',use_this_rtp,length(p.stemp))

iaX = input('Enter profiles to launch ... be nice and launch at most 100 - 250 profiles (100 * 89 chunks = 8900 processors) : [iS:iE] : ')
iS = iaX(1);
iE = iaX(2);
for ii = iS: iE
  launcher = ['!sbatch -p high_mem                              --array=1-89       sergio_matlab_jobB.sbatch 9 ' num2str(ii)];
  fprintf(1,'%s \n',launcher);
  eval(launcher)
  pause(0.1)
end
disp('wait a while (about an hour to do 100 profiles), then read them in with read_disort_chunks.m')
