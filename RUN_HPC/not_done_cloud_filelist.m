disp(' run not_done_cloud_filelist.m first ')
disp(' run not_done_cloud_filelist.m first ')
disp('RET to continue'); pause

iS =  0001; iE = 2500;
iS =  5001; iE = 7377;

iaArray = iS : iE;
iaX     = zeros(size(iaArray));

thedir = dir('JUNK/rad.dat*');
for ii = 1 : length(thedir)
  zname = thedir(ii).name;
  if zname(end:end) == 'D'
    %% this is rad.datX_CLD, ignore
  else
    zname = str2num(zname(8:end));
    boo = find(iaArray == zname);
    if length(boo) > 0
      iaX(boo) = 1;
    end
  end
end

bad = find(iaX == 0);
bad = iaArray(bad)

if length(bad) > 0
  fid = fopen('notdone.txt','w');
  for ii = 1 : length(bad)
    fprintf(fid,'%4i \n',bad(ii));
  end
  fclose(fid)
end
