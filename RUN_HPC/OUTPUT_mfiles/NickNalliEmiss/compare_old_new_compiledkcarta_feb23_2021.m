addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlab2012/rtptoolsV201/

clear all; clc; figure(1); clf; 
load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

dirMR  = 'Feb23_2021/MasudaEmis_ResetTWV';
dirMRx = 'Feb23_2021/MasudaEmis_ResetTWV_TestNewCompiledKCARTA';

dirMRx = 'Feb23_2021/JUNK/';
dirMR  = 'Feb23_2021/MasudaEmis_ResetTWV';
dirMR  = 'Feb23_2021/MasudaEmis_OrigERA';
dirMR  = 'Feb23_2021/NalliEmis_ResetTWV';  %% WHOA, looks like I used  iDoRad = 3; in set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m
dirMR  = 'Feb23_2021/NalliEmis_OrigERA';   %% WHOA, looks like I used  iDoRad = 3; in set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m
dirMR  = 'Feb23_2021/NalliEmis_ResetTWV';  %% fixed whew needed iDoRad = 10; in set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m
dirMR  = 'Feb23_2021/NalliEmis_OrigERA';   %% fixed whew needed iDoRad = 10; in set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m

dirMR

for ii = 1 : 186
  fname = [dirMR '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat']; %% for Masuda, Old Nalli (in the _OLD_Oops_forgot_to_set_iDoRad\=10)
  fname = [dirMR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];           %% for redone Nalli
  x = load(fname);
  radMR(:,ii) = x.rKc;

  fname = [dirMRx '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];  
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
plot(fKc,nanmean(tMR'-tMRx'),'b',fKc,nanstd(tMR'-tMRx'),'r');
title('Bias and StdDev'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(fKc) max(fKc)])

