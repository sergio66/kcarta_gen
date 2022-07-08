addpath /asl/matlib/h4tools

disp('make sure you run rmer_slurm_JUNKrad.sc before hand')
disp('make sure you run rmer_slurm_JUNKrad.sc before hand')
disp('make sure you run rmer_slurm_JUNKrad.sc before hand')

disp('ret to continue'); pause

set_rtp
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %4i profiles \n',use_this_rtp,length(p.stemp))
disp('  so for example if iaChunkSize == 50 and there are 5000 profiles')
disp('  then we need to launch 5000/50 = 100 jobs [01:10] [11:20] ... [95:100]  and this code loops over freq chunks 1-89')

iaX = input('Enter profiles to launch ... be nice and launch at most 10 proffilechunks (10 profilechunks * 89 chunks = 890 processors) : [iS:iE] : ')
iS = iaX(1);
iE = iaX(2);
for ii = 1:89
  launcher = ['!sbatch -p high_mem --array=' num2str(iS) '-' num2str(iE) ' sergio_matlab_jobB.sbatch 10 ' num2str(ii)];
  fprintf(1,'%s \n',launcher);
  eval(launcher)
  pause(0.1)
end
disp('wait a while (about an hour to do 100 profiles), then read them in with read_disort_chunks.m')
