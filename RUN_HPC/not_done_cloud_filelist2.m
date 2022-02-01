iMax = 72*64;
dirout = ['JUNK/AIRS_gridded_Oct2020_trendsonly/'];
dirout = ['JUNK/AIRS_gridded_Oct2020_startSept2002_trendsonly/AllDemJacs/'];

iNotFound = 0;
iSmallFile = 0;
iFound = 0;
iaFound = zeros(1,iMax);

fid = fopen('JUNK/notdone.txt','w');

for ii = 1 : iMax
  fname = [dirout '/individual_prof_convolved_kcarta_airs_100.mat'];

  fname = [dirout '/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];
  fnamex = [dirout '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];

  MinSize = 6000000; %% cloud jacobians!!!!  iSetCoud = 9999
  %MinSize = 9000000; %% cloud jacobians!!!!  iSetCloud = -1

  if exist(fname) & exist(fnamex)
    thedir = dir(fname);
    if thedir.bytes < MinSize
      %rmer = ['!/bin/rm ' fname ' ' fnamex];
      fprintf(1,'Small file %6i %12i \n',ii,thedir.bytes)
      % eval(rmer)
      iSmallFile = iSmallFile + 1;
      fprintf(fid,'%6i \n',ii)
    else
      iFound = iFound + 1;
      iaFound(ii) = 1;
    end
  elseif ~exist(fname) & exist(fnamex)
    iNotFound = iNotFound + 1;
    fprintf(fid,'%6i \n',ii);
    rmer = ['!/bin/rm ' fnamex];eval(rmer)
    fprintf(1,' notfound A %6i \n',ii);
  elseif exist(fname) & ~exist(fnamex)
    iNotFound = iNotFound + 1;
    fprintf(fid,'%6i \n',ii);
    rmer = ['!/bin/rm ' fname];eval(rmer)
    fprintf(1,' notfound B %6i \n',ii);
  elseif ~exist(fname) & ~exist(fnamex)
    iNotFound = iNotFound + 1;
    fprintf(fid,'%6i \n',ii);
    fprintf(1,' notfound both %6i \n',ii);
  end
end
fclose(fid);

fprintf(1,'of expected %6i files, found %6i, small = %6i, not found = %6i \n',iMax,iFound,iSmallFile,iNotFound)

disp(' ')
disp('if many files remain to be done can run eg     sbatch --array=1-(iNotFound+iSmallFile)/iChunkSize sergio_matlab_jobB.sbatch 5')
disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[badrow,badcol] = find(iaFound == 0);
baddy = unique(badcol)'
if length(baddy) > 0
  %for bb = 1 : length(baddy)
  %  str = ['sbatch -exclude=cnode203,cnode204,cnode260,cnode267 --array=' num2str(baddy(bb)) ' sergio_matlab_jobB.sbatch'];
  %  fprintf(fid,'%s \n',str);
  %end

  fid = fopen('notdone_filelist.sc','w');
  fprintf(1,'found that %4i of %4i timesteps did not finish : see badanom.sc \n',length(baddy),72*64)
  str = ['sbatch --account=pi_strow --exclude=cnode[204,225,267] --array='];
  str = ['sbatch --account=pi_strow  --array='];
  fprintf(fid,'%s',str);
  iX = nice_output(fid,baddy);   %% now put in continuous strips
  fprintf(1,'length(badanom) = %4i Num Continuous Strips = %4i \n',length(baddy),iX)
  str = [' sergio_matlab_jobB.sbatch 1'];
  %% going to call individual processors to do individual profiles
  fprintf(fid,'%s \n',str);
  fclose(fid);
  disp(' ')
  disp('  >>> or can run notdone_filelist.sc to run individual profiles not done')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
