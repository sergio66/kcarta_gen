%addpath /asl/packages/rtp_prod2/util/
addpath /home/sergio/MATLABCODE/matlib/rtp_prod2/util/

junk = clock; junk(6) = round(junk(6)); fprintf(1,'JOB started %4i/%2i/%2i %2i:%2i:%2i \n',junk); clear junk

% grab the slurm array index for this process
slurmindex = str2num(getenv('SLURM_ARRAY_TASK_ID'));   % 0-19999

% collect some system parameters to log
[~, hostname] = system('hostname');
slurm_job_id = getenv('SLURM_JOB_ID');
slurm_array_job_id = getenv('SLURM_ARRAY_JOB_ID');
fprintf(1, '*** Hostname: %s\tJobID: %s\tArray JobID: %s\n', hostname, slurm_job_id, slurm_array_job_id);
       
slurm_job_partition = getenv('SLURM_JOB_PARTITION');
slurm_restart_count = getenv('SLURM_RESTART_COUNT');
fprintf(1, '*** Partition: %s\tRestart Count: %s\n', slurm_job_partition, slurm_restart_count);
	      
slurm_submit_host = getenv('SLURM_SUBMIT_HOST');
slurm_submit_dir = getenv('SLURM_SUBMIT_DIR');
fprintf(1, '*** Submit host: %s\tSubmit dir: %s\n', slurm_submit_host, slurm_submit_dir);
		     
[sID, sTempPath] = genscratchpath();
fprintf(1, '*** Temp path: %s\tTemp sID: %s\n', sTempPath, sID);
fprintf(1, '*** Task run start %s\n', char(datetime('now')));

disp('space available on /scratch')
sizer = ['!df -h /scratch'];
eval(sizer)

if length(slurm_job_id) > 0
  tempscratchdir = sprintf('/scratch/%d', str2num(slurm_job_id));
else
  tempscratchdir = sprintf('/scratch/%d', str2num(slurm_job_id));
  tempscratchdir = sprintf('/umbc/lustre/strow/sergio/scratch/%d', str2num(slurm_job_id));
end

ee = exist(tempscratchdir, 'dir');
fprintf(1,'existence of tempscratchdir %s = %3i \n',tempscratchdir,ee);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%% this is from Howard
%function ctmp = ccsds_tmpfile

  jobid = str2num(getenv('SLURM_JOB_ID'));         % job ID
  jarid = str2num(getenv('SLURM_ARRAY_TASK_ID'));  % job array ID
  procid = str2num(getenv('SLURM_PROCID'));        % relative process ID
  if isempty(jobid) || isempty(jarid) || isempty(procid)  || exist(sprintf('/scratch/%d', jobid), 'dir') == 0
    fprintf(1, 'warning: using current directory for CCSDS temp file\n')
    rng('shuffle');
    ctmp = sprintf('ccsds_%05d.tmp', randi(99999));
  else
    ctmp = sprintf('/scratch/%d/ccsds_%03d_%03d.tmp', jobid, jarid, procid);
  end
%}
