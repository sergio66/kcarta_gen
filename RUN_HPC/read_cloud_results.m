addpath /home/sergio/KCARTA/MATLAB
i2 = input('Enter number of files to read from JUNK : ');

for ii = 1 : i2
  fname = ['JUNK/rad.dat' num2str(ii)];
  [d,w] = readkcstd_smart(fname);
  [mm,nn] = size(d);
  d = d(:,nn);
  data(:,ii) = d;
end

addpath /home/sergio/MATLABCODE
[fc,rc] = convolve_airs(w,data,1:2378);
tc = rad2bt(fc,rc);

%% [h,ha,p,pa] = rtpread('sartaANDkcarta_op_cloud.rtp');
%% plot(p.stemp,tc(1291,:),'o')