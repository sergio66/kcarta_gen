addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlab2012/rtptoolsV201/
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

[h,ha,pM0,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_era_Apr23_2021.rp.rtp');
[h,ha,pMR,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_resetTWV_era_Apr23_2021.rp.rtp');
[h,ha,pN0,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_era_Apr23_2021.rp.rtp');
[h,ha,pNR,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_resetTWV_era_Apr23_2021.rp.rtp');
[length(pM0.stemp)  length(pMR.stemp) length(pN0.stemp)  length(pNR.stemp)]

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75,ppmvSURF] = layers2ppmv(h,pMR,1:length(pMR.stemp),2);

dirM0  = 'Apr23_2021/MasudaEmis_OrigERA';
dirMR  = 'Apr23_2021/MasudaEmis_ResetTWV';
dirN0  = 'Apr23_2021/NalliEmis_OrigERA';
dirNR  = 'Apr23_2021/NalliEmis_ResetTWV';
dirMRx = 'Apr23_2021/MasudaEmis_ResetTWV_LMCO2OD/';

iMax = 168;
iMax = 186;
iMax = length(pM0.stemp);

for ii = 1 : iMax
  if mod(ii,1000) == 0
    fprintf(1,'x')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  fname = [dirM0 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radM0(:,ii) = x.rKc;

  fname = [dirMR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radMR(:,ii) = x.rKc;

  fname = [dirN0 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radN0(:,ii) = x.rKc;

  fname = [dirNR '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radNR(:,ii) = x.rKc;

  fname = [dirMRx '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radMRx(:,ii) = x.rKc;
end
fprintf(1,'\n')

fKc  = x.fKc(ichan);
tobs = rad2bt(h.vchan,pM0.robs1(:,1:iMax));
tM0  = rad2bt(x.fKc(ichan),radM0(ichan,:));
tMR  = rad2bt(x.fKc(ichan),radMR(ichan,:));
tN0  = rad2bt(x.fKc(ichan),radN0(ichan,:));
tNR  = rad2bt(x.fKc(ichan),radNR(ichan,:));
tMRx = rad2bt(x.fKc(ichan),radMRx(ichan,:));

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

figure(3); 
i800  = find(fKc >= 821,1);
i1231 = find(fKc >= 1231.1,1);
plot(pM0.stemp,tobs(i800,:)-tM0(i800,:),'b.',pMR.stemp,tobs(i800,:)-tMR(i800,:),'c.',...
     pN0.stemp,tobs(i800,:)-tN0(i800,:),'r.',pNR.stemp,tobs(i800,:)-tNR(i800,:),'m.')
addpath /home/sergio/MATLABCODE/SHOWSTATS
[n,nx,ny,nmeanM0,nstdM0] = myhist2d(pM0.stemp,tobs(i800,:)-tM0(i800,:),270:5:310,-4:0.25:+4);
[n,nx,ny,nmeanMR,nstdMR] = myhist2d(pMR.stemp,tobs(i800,:)-tMR(i800,:),270:5:310,-4:0.25:+4);
[n,nx,ny,nmeanN0,nstdN0] = myhist2d(pN0.stemp,tobs(i800,:)-tN0(i800,:),270:5:310,-4:0.25:+4);
[n,nx,ny,nmeanNR,nstdNR] = myhist2d(pNR.stemp,tobs(i800,:)-tNR(i800,:),270:5:310,-4:0.25:+4);
plot(270:5:310,nmeanM0,'b.-',270:5:310,nmeanMR,'c.-',270:5:310,nmeanN0,'r.-',270:5:310,nmeanNR,'m.-'); grid
  hl = legend('Masuda 0','Masuda Reset','Nalli0','Nalli Reset','location','best','fontsize',10); title('820 cm-1 bias'); xlabel('stemp');

figure(4);
addpath /home/sergio/MATLABCODE/SHOWSTATS
[n,nx,ny,nmeanM0,nstdM0] = myhist2d(pM0.satzen,tobs(i800,:)-tM0(i800,:),0:2.5:50,-4:0.25:+4);
[n,nx,ny,nmeanMR,nstdMR] = myhist2d(pMR.satzen,tobs(i800,:)-tMR(i800,:),0:2.5:50,-4:0.25:+4);
[n,nx,ny,nmeanN0,nstdN0] = myhist2d(pN0.satzen,tobs(i800,:)-tN0(i800,:),0:2.5:50,-4:0.25:+4);
[n,nx,ny,nmeanNR,nstdNR] = myhist2d(pNR.satzen,tobs(i800,:)-tNR(i800,:),0:2.5:50,-4:0.25:+4);
plot(0:2.5:50,nmeanM0,'b.-',0:2.5:50,nmeanMR,'c.-',0:2.5:50,nmeanN0,'r.-',0:2.5:50,nmeanNR,'m.-'); grid
  hl = legend('Masuda 0','Masuda Reset','Nalli0','Nalli Reset','location','best','fontsize',10); title('820 cm-1 bias'); xlabel('satzen')

figure(5);
addpath /home/sergio/MATLABCODE/SHOWSTATS
[n,nx,ny,nmeanM0,nstdM0] = myhist2d(pM0.wspeed,tobs(i800,:)-tM0(i800,:),0:1:20,-4:0.25:+4);
[n,nx,ny,nmeanMR,nstdMR] = myhist2d(pMR.wspeed,tobs(i800,:)-tMR(i800,:),0:1:20,-4:0.25:+4);
[n,nx,ny,nmeanN0,nstdN0] = myhist2d(pN0.wspeed,tobs(i800,:)-tN0(i800,:),0:1:20,-4:0.25:+4);
[n,nx,ny,nmeanNR,nstdNR] = myhist2d(pNR.wspeed,tobs(i800,:)-tNR(i800,:),0:1:20,-4:0.25:+4);
plot(0:1:20,nmeanM0,'b.-',0:1:20,nmeanMR,'c.-',0:1:20,nmeanN0,'r.-',0:1:20,nmeanNR,'m.-'); grid
  hl = legend('Masuda 0','Masuda Reset','Nalli0','Nalli Reset','location','best','fontsize',10); title('820 cm-1 bias'); xlabel('wspeed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
plot(fKc,nanmean(tobs'-tMR'),'b',fKc,nanmean(tobs'-tMRx'),'r',...
     fKc,nanstd(tobs'-tMR'),'c--',fKc,nanstd(tobs'-tMRx'),'m--','linewidth',2); grid
plotaxis2;
hl = legend('LBLRTM','LM','location','best','fontsize',10);
xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])
%axis([2190 2440 -2 +2])
%axis([650 850   -2 +2])

boo.fKc = fKc;
boo.lblrtm_bias = nanmean(tobs'-tMR');
boo.lblrtm_std  = nanstd(tobs'-tMR');
boo.lm_bias = nanmean(tobs'-tMRx');
boo.lm_std  = nanstd(tobs'-tMRx');
%save('lblrtm_vs_lm.mat','-struct','boo');

plot(boo.fKc,boo.lblrtm_bias,'b',boo.fKc,boo.lm_bias,'r',...
     boo.fKc,boo.lblrtm_std,'c--',boo.fKc,boo.lm_std,'m--','linewidth',2); grid
plotaxis2;
hl = legend('LBLRTM','LM','location','best','fontsize',10);
xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tropics = find(abs(pMR.rlat) < 30);
midlat  = find(abs(pMR.rlat) >= 30 & abs(pMR.rlat) < 60);
polar   = find(abs(pMR.rlat) >= 60);

ind = tropics; whos ind; mean(pMR.stemp(ind))
  booI.fKc = fKc;
  booI.obs    = nanmean(tobs(:,ind)');
  booI.lblrtm = nanmean(tMR(:,ind)');
  booI.lm     = nanmean(tMRx(:,ind)');
  booI.lblrtm_bias = nanmean(tobs(:,ind)'-tMR(:,ind)');
  booI.lblrtm_std  = nanstd(tobs(:,ind)'-tMR(:,ind)');
  booI.lm_bias = nanmean(tobs(:,ind)'-tMRx(:,ind)');
  booI.lm_std  = nanstd(tobs(:,ind)'-tMRx(:,ind)');

figure(7)
  plot(booI.fKc,booI.obs,'k',booI.fKc,booI.lblrtm,'b',booI.fKc,booI.lm,'r','linewidth',2); grid
  plotaxis2;
  hl = legend('Obs','LBLRTM','LM','location','best','fontsize',10);
  xlabel('wavenumber'); ylabel('BT(K)'); xlim([min(h.vchan) max(h.vchan)])

figure(8)
  plot(booI.fKc,booI.lblrtm_bias,'b',booI.fKc,booI.lm_bias,'r',...
       booI.fKc,booI.lblrtm_std,'c--',booI.fKc,booI.lm_std,'m--','linewidth',2); grid
  plotaxis2;
  hl = legend('LBLRTM','LM','location','best','fontsize',10);
  xlabel('wavenumber'); ylabel('Obs-Cal BT(K)'); xlim([min(h.vchan) max(h.vchan)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

figure(4); plot(results.masuda_efreq,results.masuda_emis,'b',results.nalliE_efreq,results.nalliE_emis,'r')

% save nalli_Apr23_2021.mat results
