%% see ~/KCARTA/INCLUDE/pre_defined.param
%     DATA kaMinFr        /   15.0000,   30.0000,   50.0000,   80.0000,
%     c                    140.0000,  300.0000,  500.0000,  805.0000,
%     c                    2830.0000,  3550.0000,  5550.0000,  8250.0000,
%     c                    12000.0000,  25000.0000  /
%      DATA kaMaxFr        /   30.0000,   50.0000,   80.0000,  150.0000,
%     c                    310.0000,  510.0000,  880.0000,  2830.0000,
%     c                    3580.0000,  5650.0000,  8400.0000,  12250.0000,
%     c                    25000.0000,  44000.0000  /

f1 = [15 30 50 80 140 300 500 605 2830 3550 5550 8250 12000 25000];
f2 = f1(2:end); f2(length(f1)) = 44000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iHIT = 2020;
iHIT = 2016;
if iHIT == 2016
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res';
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res';
elseif iHIT == 2020
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20_orig605_805res';
end
outnml = 'junk.nml';

iaDo = 1 : length(f1);
iaDo = [8];

iMakeData = -1;
iMakeData = +1;
if iMakeData > 0
  for iii = 1 : length(iaDo)
    ii = iaDo(iii);
    sedder = ['!sed -e "s/FF1/' num2str(f1(ii)) '/g"  -e "s/FF2/' num2str(f2(ii)) '/g" '];
    sedder = [sedder ' tempate_quickuse_l2s_kcVERS.nml  > ' outnml];
    eval(sedder);
  
    outname = ['../L2SComparisons/l2s_kc122_H' num2str(iHIT-2000) '_' num2str(f1(ii)) '_' num2str(f2(ii)) '.dat'];
    if exist(outname)
      rmer = ['!rm ' outname];
      eval(rmer);
    end
  
    fprintf(1,'%2i %s \n',ii,outname);
    kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(ii)];
    eval(kcartaer);
  
    rmer = ['!rm status' num2str(ii) ' junk.nml'];
    eval(rmer)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
iReadData = +1;
if iReadData > 0
  clf
  dall = [];
  wall = [];
  d6qc = [];
  fqc = [];

  for iii = 1 : length(iaDo)
    ii = iaDo(iii);
  
    outname = ['../L2SComparisons/l2s_kc122_H16_'                     num2str(f1(ii)) '_' num2str(f2(ii)) '.dat'];
    outname = ['../L2SComparisons/l2s_kc122_H' num2str(iHIT-2000) '_' num2str(f1(ii)) '_' num2str(f2(ii)) '.dat'];
    if exist(outname)
      [d,w,caVersion, detail] = readkcstd_detail(outname);
      iaGasID = detail.iaGasID;

      saver = ['save ../L2SComparisons/l2s_kc122_H' num2str(iHIT-2000) '_' num2str(f1(ii)) '_' num2str(f2(ii)) '.mat d w iaGasID'];
      eval(saver);
      %[d,w] = readkcstd(outname);
      ptspacing(ii) = mean(diff(w));
      meanw(ii)     = mean(w);

      wall = [wall w];
      dall = [dall; d(:,1:6)];
     
%% AVIRIS-NG has 5 nm spec from 360-2500 nm = 425 chaans      
      l1 = 10000/max(w);
      l2 = 10000/min(w);
      lav = 10000/mean(w);
      %% want spacing of 5 nm
      lavp005 = lav+0.005;
      wavp005 = mean(w) - 10000/lavp005;

%      if max(w) < 3000
%        [fc,qc] = quickconvolve(w,d(:,1:6),100,100);
%      else 
%        [fc,qc] = quickconvolve(w,d(:,1:6),5000,5000);
%      end
      [fc,qc] = quickconvolve(w,d(:,1:6),wavp005,wavp005);
      meandwfc(ii) = mean(diff(fc));

      fqc = [fqc fc];
      d6qc = [d6qc; qc];

      figure(1); semilogy(w,d(:,1),'b',w,d(:,2),'g',w,d(:,3),'y',w,d(:,4),'c',w,d(:,5),'m',w,d(:,6),'r')
      hl = legend(num2str((1:6)'),'location','best','fontsize',10);
      %disp('ret to continue'); pause
      pause(0.1);
      hold on
    else
      fprintf(1,'%2i %s dne \n',ii,outname)
    end

  end

  [ptspacing; meandwfc]

  hold off
  semilogy(10000./wall,dall);       hl = legend(num2str((1:6)'),'location','best','fontsize',10);
    xlabel('Wavelength(um)'); 

  moo = find(10000./wall >= 0.36 & 10000./wall <= 2.5); %% AVIRIS-NG has 425 channels in this region
  moo = find(10000./fqc >= 0.36 & 10000./fqc <= 2.5); %% AVIRIS-NG has 425 channels in this region
  whos moo wall

  plot(10000./wall,exp(-sum(dall')));

  dnew = dall .* (ones(length(dall),1) * [1 1 1 1 1 1]);
  dnew = dall .* (ones(length(dall),1) * [1 1 1 1.2 1 1]);
  plot(10000./wall,exp(-sum(dall')),'b.-',10000./wall,exp(-sum(dnew')),'r',10000./fqc,exp(-sum(d6qc')),'k.-'); 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
