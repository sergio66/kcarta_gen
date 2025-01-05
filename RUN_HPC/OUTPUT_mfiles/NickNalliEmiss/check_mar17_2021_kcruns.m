addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlab2012/rtptoolsV201/

clear all; clf; clc; 
load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

dirMR  = 'Mar17_2021/MasudaEmis_ResetTWV_LMCO2OD';


dirMR  = 'Mar17_2021/NalliEmis_OrigERA';   %% WOW, moved original calcs to Mar17_2021/NalliEmis_OrigERA_ORIG, now FINE
dirMR  = 'Mar17_2021/NalliEmis_ResetTWV';  %% FINE
dirMR  = 'Mar17_2021/MasudaEmis_ResetTWV'; %% PWOBLEM, oops I used Nalli eiss setting for kcarta .... gives 0.1 K diffs ... mebbe did not wait long enuff, so moved original calcs, now FINE
dirMR  = 'Mar17_2021/MasudaEmis_OrigERA';  %% FINE

dirMR
dirMRx = 'Mar17_2021/JUNK';

for ii = 1 : 186
  fname = [dirMR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radMR(:,ii) = x.rKc;

  fname = [dirMRx '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radMRx(:,ii) = x.rKc;
end

fKc  = x.fKc(ichan);
tMR   = rad2bt(x.fKc(ichan),radMR(ichan,:));
tMRx  = rad2bt(x.fKc(ichan),radMRx(ichan,:));

%figure(1)
%plot(fKc,nanmean(tobs'-tMRx'),'b',fKc,nanmean(tobs'-tMR'),'r');
%title('Bias'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

figure(1)
plot(fKc,nanmean(tMR'-tMRx'),'b.-',fKc,nanstd(tMR'-tMRx'),'r');
title('Bias'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(fKc) max(fKc)])

