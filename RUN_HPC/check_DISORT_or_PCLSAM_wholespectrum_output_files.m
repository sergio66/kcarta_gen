set_rtp
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %5i profiles \n',use_this_rtp,length(p.stemp));

iDorP = input('Enter DISORT (-1) or Chou/PCLSAM (+1, default) : ');
if length(iDorP) == 0
  iDorP = +1;
end

%iNum = input('Enter number of files to look for ');
iNum = length(p.stemp);

iaFound = zeros(1,iNum);
for ii = 1 : iNum
  if iDorP < 0
    fname = ['JUNK/rad.dat' num2str(ii) '_0605'];
  else
    fname = ['JUNK/rad.dat' num2str(ii)];
  end

  if exist(fname)
    thedir = dir(fname);
    if thedir.bytes < 3566468
      iaFound(ii) = -1;
    else
      iaFound(ii) = +1;
    end
  end
end
plot(iaFound);
made  = find(iaFound == 1);     fprintf(1,'files yes made = %5i of %5i \n',length(made),iNum)
notmade  = find(iaFound == 0);  fprintf(1,'files not made = %5i of %5i \n',length(notmade),iNum)
toosmall = find(iaFound == -1); fprintf(1,'files too small = %5i of %5i \n',length(toosmall),iNum)

iRemove = input('remove the files that are too small? (-1 default /+1) : ');
if length(iRemove) == 0
  iRemove = -1;
end
if iRemove > 0
  for ii = 1 : length(toosmall)
    if iDorP < 0
      fname = ['JUNK/rad.dat' num2str(toosmall(ii)) '_0605'];
    else
      fname = ['JUNK/rad.dat' num2str(toosmall(ii))];
    end
    if exist(fname)
      lser = ['!ls -lt ' fname];
      eval(lser);
      rmer = ['!/bin/rm ' fname];
      eval(rmer);
    end
    if iDorP < 0
      fname = ['JUNK/rad.dat' num2str(toosmall(ii)) '_0605_CLD'];
    else
      fname = ['JUNK/rad.dat' num2str(toosmall(ii)) '_CLD'];
    end
    if exist(fname)
      lser = ['!ls -lt ' fname];
      eval(lser);
      rmer = ['!/bin/rm ' fname];
      eval(rmer);
    end
  end
end

baddy = find(iaFound <= 0);
if length(baddy) > 0
  excludelist = ' ';
  excludelist = 'cnode013,cnode018,cnode019';

  fid = fopen('missing_disort.sc','w');

  %for ii = 1 : length(baddy)
  %  str = ['sbatch -p high_mem  --exclude=' excludelist ' --array=' num2str(baddy(ii)) ' sergio_matlab_jobB.sbatch 10'];
  %  fprintf(fid,'%s \n',str);
  %end

  fprintf(1,'found that %5i of %5i DISORT whole spectra did not finish : see missing_disort.sc \n',length(baddy),iNum)
  str = ['sbatch -p high_mem --exclude=' excludelist ' --array='];
  fprintf(fid,'%s',str);
  iX = nice_output(fid,baddy);   %% now put in continuous strips
  fprintf(1,'length(badanom) = %4i Num Continuous Strips = %4i \n',length(baddy),iX)
  if iDorP < 0
    str = [' sergio_matlab_jobB.sbatch 10'];
  else
    str = [' sergio_matlab_jobB.sbatch 1'];
  end
  fprintf(fid,'%s \n',str);

  fclose(fid);
  disp('look at missing_disort.sc : if too many processes you can edit it and split it between eg cpu2021 and high_mem');
end

