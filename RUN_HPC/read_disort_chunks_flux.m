addpath /home/sergio/KCARTA/MATLAB
wdis = zeros(89*10000,1);
fluxdis = zeros(89*10000,1);

iProf = input('Enter rtp profile number : ');

iCnt = 0;
for ii = 1 : 89
  ind = (1:10000) + (ii-1)*10000;
  f1 = 605 + (ii-1)*25;
  fname = ['JUNK/rad.dat' num2str(iProf) '_' num2str(f1,'%04d') '_OLR'];
  %fprintf(1,'%4i %s %3i \n',f1,fname,exist(fname));
  if exist(fname)
    iCnt = iCnt + 1;
    [dx,wx] = readkcflux(fname);
    wdis(ind) = wx;
    fluxdis(ind) = dx;
    %plot(wx,rad2bt(wx,dx)); pause
  end
end

fprintf(1,'read in %2i of 89 DISORT chunks into [wdis,fluxdis] \n',iCnt)
if iCnt > 0
  plot(wdis,fluxdis);
end
