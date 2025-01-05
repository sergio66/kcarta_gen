addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlab2012/rtptoolsV201/

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

[h,ha,pM0,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_era_Feb23_2021.rp.rtp');
[h,ha,pMR,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_resetTWV_era_Feb23_2021.rp.rtp');
[h,ha,pN0,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_era_Feb23_2021.rp.rtp');
[h,ha,pNR,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_resetTWV_era_Feb23_2021.rp.rtp');

dirM0 = 'Feb23_2021/MasudaEmis_OrigERA';
dirMR = 'Feb23_2021/MasudaEmis_ResetTWV';
dirN0 = 'Feb23_2021/NalliEmis_OrigERA';
dirNR = 'Feb23_2021/NalliEmis_ResetTWV';

for ii = 1 : 186
  fname = [dirM0 '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  x = load(fname);
  radM0(:,ii) = x.rKc;

  fname = [dirMR '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  x = load(fname);
  radMR(:,ii) = x.rKc;

  fname = [dirN0 '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat']; %% old, done wrong
  fname = [dirN0 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];           %% fixed 
  x = load(fname);
  radN0(:,ii) = x.rKc;

  fname = [dirNR '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  fname = [dirNR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radNR(:,ii) = x.rKc;
end

fKc  = x.fKc(ichan);
tobs = rad2bt(h.vchan,pM0.robs1);
tM0  = rad2bt(x.fKc(ichan),radM0(ichan,:));
tMR  = rad2bt(x.fKc(ichan),radMR(ichan,:));
tN0  = rad2bt(x.fKc(ichan),radN0(ichan,:));
tNR  = rad2bt(x.fKc(ichan),radNR(ichan,:));

figure(1)
plot(fKc,nanmean(tobs'-tM0'),'b',fKc,nanmean(tobs'-tN0'),'r',...
     fKc,nanmean(tobs'-tMR'),'c',fKc,nanmean(tobs'-tNR'),'m'); grid
hl = legend('Masuda ERA','Nalli ERA','Masuda ResetTWV','Nalli ResetTWV','location','best');
title('Bias'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

figure(2)
plot(fKc,nanstd(tobs'-tM0'),'b',fKc,nanstd(tobs'-tN0'),'r',...
     fKc,nanstd(tobs'-tMR'),'c',fKc,nanstd(tobs'-tNR'),'m'); grid
hl = legend('Masuda ERA','Nalli ERA','Masuda ResetTWV','Nalli ResetTWV','location','best');
title('Std dev'); xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

mmW_ERA   = mmwater_rtp(h,pM0);
mmW_Reset = mmwater_rtp(h,pMR);

results.stemp_ERA = pM0.stemp;
results.mmw_ERA   = mmW_ERA;
results.stemp_ResetERA = pMR.stemp;
results.mmw_ResetERA   = mmW_Reset;
results.solzen = pM0.solzen;
results.satzen = pM0.satzen;
results.scanang = pM0.scanang;
results.wspeed = pM0.wspeed;
results.rlon = pM0.rlon;
results.rlat = pM0.rlat;
results.masuda_efreq = pM0.efreq(1:30,:);
results.masuda_emis  = pM0.emis(1:30,:);
results.masuda_rho   = pM0.rho(1:30,:);
results.nalliE_efreq  = pN0.efreq(1:30,:);
results.nalliE_emis   = pN0.emis(1:30,:);
results.nalliE_rho    = pN0.rho(1:30,:);
results.bad = find(pN0.nemis ~= 30);
results.vhan = h.vchan;
results.robs = pM0.robs1;
results.rcalc_masuda_ERA      = radM0(ichan,:);
results.rcalc_masuda_ResetERA = radMR(ichan,:);
results.rcalc_nalli_ERA       = radN0(ichan,:);
results.rcalc_nalli_ResetERA  = radNR(ichan,:);

figure(3); plot(results.masuda_efreq,results.masuda_emis,'b',results.nalliE_efreq,results.nalliE_emis,'r')

% save nalli_Feb23_2021.mat results
