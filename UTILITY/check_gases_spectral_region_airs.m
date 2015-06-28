%% this is to look for gases in spectral regions

if ~exist('kcdata') & ~exist('kcfreq')
  addpath ../MATLAB
  addpath /home/sergio/MATLABCODE
  [kcdata,kcfreq] = readkcstd('l2s_kc116_72gases.dat');  
  [fc,qc] = convolve_airs(kcfreq,kcdata,1:2378);
  [fc,qcsumdata] = convolve_airs(kcfreq,sum(kcdata'),1:2378);
end

ff = input('enter AIRS channel center freq you want to check : ');
woo = find(fc >= ff-2.5 & fc <= ff+2.5); whos woo
for ii = 1 : 72
  %% kcdata already contains ONE unit of gas ii; add on 9 more units to see if there is a change in total L2S
  [fc,qcx] = convolve_airs(kcfreq,sum(kcdata') + 9*kcdata(:,ii)',1:2378);
  plot(fc(woo),exp(-qc(woo,ii)),'bo-',fc(woo),exp(-qcsumdata(woo)),'ko-',...
       fc(woo),exp(-10*qc(woo,ii)),'ro-',fc(woo),exp(-qcx(woo)),'k--','linewidth',2)
  ax = axis; line([ff ff],[ax(3) ax(4)],'color','green','linewidth',3)
  hl = legend('gas(ii)','total us std','gas(ii)x10','total us std + gas(ii)x10','location','southwest');
  set(hl,'fontsize',10);
  title(num2str(ii)); grid; pause;
end
