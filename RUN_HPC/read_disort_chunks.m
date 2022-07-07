addpath /home/sergio/KCARTA/MATLAB
wdis = zeros(89*10000,1);

iaProf = input('Enter rtp profile number list (one for only one, 2 for (Start: Stop) : ');

if length(iaProf) == 1
  iaList = iaProf;
elseif length(iaProf) == 2
  iaList = iaProf(1):iaProf(2);
end

raddis = zeros(89*10000,length(iaList));

for iPP = 1 : length(iaList)
  iProf = iaList(iPP);
  iCnt = 0;
  iaFound = zeros(1,89);
  for ii = 1 : 89
    ind = (1:10000) + (ii-1)*10000;
    f1 = 605 + (ii-1)*25;
    fname = ['JUNK/rad.dat' num2str(iProf) '_' num2str(f1,'%04d')];
    %fprintf(1,'%4i %s %3i \n',f1,fname,exist(fname));
    if exist(fname)
      iCnt = iCnt + 1;
      iaFound(ii) = 1;
      [dx,wx] = readkcstd(fname);
      wdis(ind) = wx;
      raddis(ind,iPP) = dx;
      %plot(wx,rad2bt(wx,dx)); pause
    end
  end
  
  fprintf(1,'read in %2i of 89 DISORT chunks into [wdis,raddis] \n',iCnt)
  if iCnt > 0
    plot(wdis,rad2bt(wdis,raddis));
  end
  if iCnt < 89
    disp('missing chunks')
    boo = find(iaFound == 0)
  end
end
