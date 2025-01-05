addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlab2012/rtptoolsV201/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

[h,ha,pMR,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_resetTWV_era_Apr23_2021.rp.rtp');
[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,pMR,1:length(pMR.stemp),2);

dirMR  = 'Apr23_2021/MasudaEmis_ResetTWV';
dirNew  = '../';

iMax = 168;
iMax = 186;
iMax = length(pMR.stemp);

iCount = 0;
for ii = 1 : iMax
  if mod(ii,1000) == 0
    fprintf(1,'x')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  fnameMR = [dirMR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  fnameNew = [dirNew '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];

  if exist(fnameMR) & exist(fnameNew)
    iCount = iCount + 1;
    iaFound(iCount) = ii;
    x = load(fnameMR);
    radMR(:,iCount) = x.rKc;
    x = load(fnameNew);
    radNew(:,iCount) = x.rKc;
  end

end
fprintf(1,'\n')

iMax = iCount;
[length(pMR.stemp) iMax]

fKc  = x.fKc(ichan);
tobs = rad2bt(h.vchan,pMR.robs1(:,iaFound));
tMR  = rad2bt(x.fKc(ichan),radMR(ichan,:));
tNew = rad2bt(x.fKc(ichan),radNew(ichan,:));

figure(1)
plot(fKc,nanmean(tobs'-tMR'),'b',fKc,nanmean(tobs'-tNew'),'r',...
     fKc,nanstd(tobs'-tMR'),'c',fKc,nanstd(tobs'-tNew'),'m'); grid
hl = legend('121 mean','122 mean','121 std','122 std','location','best');
title('Bias'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

figure(2)
plot(fKc,nanmean(tNew'-tMR'),'b',fKc,nanstd(tNew'-tMR'),'r'); grid
hl = legend('122-121 mean','122-121 std','location','best');
title('Bias'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

