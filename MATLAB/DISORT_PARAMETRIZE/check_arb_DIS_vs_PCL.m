addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

iP = 1;  %% proff1
iPstr = ['Prof' num2str(iP)];

%%% see eg make_profiles_arb49.m
ctype         = [101 201];
water_cpsize0 = 11:3:29;  
ice_cpsize0   = 20:15:120;
cngwat0       = [0 1e-2 1e-1 1 10 30 100 300];
cprtop0       = 150:100:950;
scanang0      = 0:10:50;
%%% see eg make_profiles_arb49.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iVers = input('Enter files to check (1/default) 9072 profiles strictly uses cngwat.cpszie from N-D tables or (-1) 6700 randomly perturbed values from N-D tables : ');
if length(iVers) == 0
  iVers = +1;
end

if iVers == +1
  frtp   = '/home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/parametrize_pclsam_disort_AFGL_1_49.rp.rtp';
  fdis   = [iPstr '/ChouTestLookupTableArbProfiles/conv4440DISORT.mat'];
  ftest  = [iPstr '/ChouTestLookupTableArbProfiles/conv4440PCLSAM_BEST.mat'];  %% variable factor
  ftest0 = [iPstr '/ChouTestLookupTableArbProfiles/conv4440PCLSAM_NOADJ.mat'];
else
  frtp   = '/home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/parametrize_pclsam_disort_arbCZT_AFGL_1_49.rp.rtp';
  fdis   = [iPstr '/ChouTestLookupTableArbProfiles/conv4440_random_pert_DISORT.mat'];
  ftest  = [iPstr '/ChouTestLookupTableArbProfiles/conv4440_random_pert_PCLSAM_BEST.mat'];  %% variable factor
  ftest0 = [iPstr '/ChouTestLookupTableArbProfiles/conv4440_random_pert_PCLSAM_NOADJ.mat'];
  ftest0 = [iPstr '/ChouTestLookupTableArbProfiles/conv4440_random_pert_PCLSAM_BEST_VaryMatr.mat'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dis   = load(fdis);
pcl   = load(ftest);
pcl0   = load(ftest0);
[h,ha,p,pa] = rtpread(frtp);

i1231 = find(dis.fc >= 1231,1);

figure(1)
mmw = mmwater_rtp(h,p);
scatter(p.stemp,mmw,50,abs(p.plat),'filled'); colormap jet; colorbar; xlabel('stemp');ylabel('mmw'); title('colorbar - lat')

disBT = rad2bt(dis.fc,dis.qc);
pclBT = rad2bt(pcl.fc,pcl.qc);
pcl0BT = rad2bt(pcl.fc,pcl0.qc);

figure(1)
plot(dis.fc,nanmean(disBT'-pcl0BT'),'b',dis.fc,nanstd(disBT'-pcl0BT'),'c--',dis.fc,nanmean(disBT'-pclBT'),'r',dis.fc,nanstd(disBT'-pclBT'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best');; title('ALL')
disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iWater      = find(p.ctype == 101 & p.ctype2 == -9999 & p.cprtop >= 400);
iIce        = find(p.ctype == 201 & p.ctype2 == -9999 & p.cprtop <= 600);
iIceNWater  = find(p.ctype == 101 & p.ctype2 == 201);

iPolar = find(p.stemp <= 273);
iMid = find(p.stemp > 273 & p.stemp <= 293);
iHot = find(p.stemp > 293);

iPolar = intersect(iPolar,union(iWater,iIce));
iMid = intersect(iMid,union(iWater,iIce));
iHot = intersect(iHot,union(iWater,iIce));

figure(2)
plot(dis.fc,nanmean(disBT(:,iPolar)'-pcl0BT(:,iPolar)'),'b',dis.fc,nanstd(disBT(:,iPolar)'-pcl0BT(:,iPolar)'),'c--',dis.fc,nanmean(disBT(:,iPolar)'-pclBT(:,iPolar)'),'r',dis.fc,nanstd(disBT(:,iPolar)'-pclBT(:,iPolar)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Polar')

figure(3)
plot(dis.fc,nanmean(disBT(:,iMid)'-pcl0BT(:,iMid)'),'b',dis.fc,nanstd(disBT(:,iMid)'-pcl0BT(:,iMid)'),'c--',dis.fc,nanmean(disBT(:,iMid)'-pclBT(:,iMid)'),'r',dis.fc,nanstd(disBT(:,iMid)'-pclBT(:,iMid)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Mid')

figure(4)
plot(dis.fc,nanmean(disBT(:,iHot)'-pcl0BT(:,iHot)'),'b',dis.fc,nanstd(disBT(:,iHot)'-pcl0BT(:,iHot)'),'c--',dis.fc,nanmean(disBT(:,iHot)'-pclBT(:,iHot)'),'r',dis.fc,nanstd(disBT(:,iHot)'-pclBT(:,iHot)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Hot')

figure(5); 
plot(p.stemp(iPolar),disBT(i1231,iPolar),'bo',p.stemp(iMid),disBT(i1231,iMid),'gx',p.stemp(iHot),disBT(i1231,iHot),'r.',p.stemp(iPolar),p.stemp(iPolar),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iPolar),p.stemp(iPolar)-disBT(i1231,iPolar),'bo',p.stemp(iMid),p.stemp(iMid)-disBT(i1231,iMid),'gx',p.stemp(iHot),p.stemp(iHot)-disBT(i1231,iHot),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iPolar),pclBT(i1231,iPolar)-disBT(i1231,iPolar),'bo',p.stemp(iMid),pclBT(i1231,iMid)-disBT(i1231,iMid),'gx',p.stemp(iHot),pclBT(i1231,iHot)-disBT(i1231,iHot),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iPolar),pclBT(i1231,iPolar)-disBT(i1231,iPolar),'bo',disBT(i1231,iMid),pclBT(i1231,iMid)-disBT(i1231,iMid),'gx',disBT(i1231,iHot),pclBT(i1231,iHot)-disBT(i1231,iHot),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
plot(dis.fc,nanmean(disBT(:,iWater)'-pcl0BT(:,iWater)'),'b',dis.fc,nanstd(disBT(:,iWater)'-pcl0BT(:,iWater)'),'c--',dis.fc,nanmean(disBT(:,iWater)'-pclBT(:,iWater)'),'r',dis.fc,nanstd(disBT(:,iWater)'-pclBT(:,iWater)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld')

figure(3)
plot(dis.fc,nanmean(disBT(:,iIce)'-pcl0BT(:,iIce)'),'b',dis.fc,nanstd(disBT(:,iIce)'-pcl0BT(:,iIce)'),'c--',dis.fc,nanmean(disBT(:,iIce)'-pclBT(:,iIce)'),'r',dis.fc,nanstd(disBT(:,iIce)'-pclBT(:,iIce)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cloud')

figure(4)
plot(dis.fc,nanmean(disBT(:,iIceNWater)'-pcl0BT(:,iIceNWater)'),'b',dis.fc,nanstd(disBT(:,iIceNWater)'-pcl0BT(:,iIceNWater)'),'c--',dis.fc,nanmean(disBT(:,iIceNWater)'-pclBT(:,iIceNWater)'),'r',dis.fc,nanstd(disBT(:,iIceNWater)'-pclBT(:,iIceNWater)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('IceNWater Cld ')

figure(5); 
plot(p.stemp(iWater),disBT(i1231,iWater),'bo',p.stemp(iIce),disBT(i1231,iIce),'gx',p.stemp(iIceNWater),disBT(i1231,iIceNWater),'r.',p.stemp(iWater),p.stemp(iWater),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iWater),p.stemp(iWater)-disBT(i1231,iWater),'bo',p.stemp(iIce),p.stemp(iIce)-disBT(i1231,iIce),'gx',p.stemp(iIceNWater),p.stemp(iIceNWater)-disBT(i1231,iIceNWater),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iWater),pclBT(i1231,iWater)-disBT(i1231,iWater),'bo',p.stemp(iIce),pclBT(i1231,iIce)-disBT(i1231,iIce),'gx',p.stemp(iIceNWater),pclBT(i1231,iIceNWater)-disBT(i1231,iIceNWater),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iWater),pclBT(i1231,iWater)-disBT(i1231,iWater),'bo',disBT(i1231,iIce),pclBT(i1231,iIce)-disBT(i1231,iIce),'gx',disBT(i1231,iIceNWater),pclBT(i1231,iIceNWater)-disBT(i1231,iIceNWater),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(dis.fc,nanmean(disBT(:,iWater)'-pcl0BT(:,iWater)'),'b',dis.fc,nanstd(disBT(:,iWater)'-pcl0BT(:,iWater)'),'c--',dis.fc,nanmean(disBT(:,iWater)'-pclBT(:,iWater)'),'r',dis.fc,nanstd(disBT(:,iWater)'-pclBT(:,iWater)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld')

szeWater   = unique(p.cpsize(iWater));
szeWater   = water_cpsize0;
wS1 = find(p.cpsize(iWater) <= szeWater(2));
wS2 = find(p.cpsize(iWater) > szeWater(2) & p.cpsize(iWater) <= szeWater(4));
wS3 = find(p.cpsize(iWater) > szeWater(4));

figure(2)
plot(dis.fc,nanmean(disBT(:,iWater(wS1))'-pcl0BT(:,iWater(wS1))'),'b',dis.fc,nanstd(disBT(:,iWater(wS1))'-pcl0BT(:,iWater(wS1))'),'c--',dis.fc,nanmean(disBT(:,iWater(wS1))'-pclBT(:,iWater(wS1))'),'r',dis.fc,nanstd(disBT(:,iWater(wS1))'-pclBT(:,iWater(wS1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld SmallDme')

figure(3)
plot(dis.fc,nanmean(disBT(:,iWater(wS2))'-pcl0BT(:,iWater(wS2))'),'b',dis.fc,nanstd(disBT(:,iWater(wS2))'-pcl0BT(:,iWater(wS2))'),'c--',dis.fc,nanmean(disBT(:,iWater(wS2))'-pclBT(:,iWater(wS2))'),'r',dis.fc,nanstd(disBT(:,iWater(wS2))'-pclBT(:,iWater(wS2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld MediumDME')

figure(4)
plot(dis.fc,nanmean(disBT(:,iWater(wS3))'-pcl0BT(:,iWater(wS3))'),'b',dis.fc,nanstd(disBT(:,iWater(wS3))'-pcl0BT(:,iWater(wS3))'),'c--',dis.fc,nanmean(disBT(:,iWater(wS3))'-pclBT(:,iWater(wS3))'),'r',dis.fc,nanstd(disBT(:,iWater(wS3))'-pclBT(:,iWater(wS3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld LargeDME ')

figure(5); 
plot(p.stemp(iWater(wS1)),disBT(i1231,iWater(wS1)),'bo',p.stemp(iWater(wS2)),disBT(i1231,iWater(wS2)),'gx',p.stemp(iWater(wS3)),disBT(i1231,iWater(wS3)),'r.',p.stemp(iWater(wS1)),p.stemp(iWater(wS1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iWater(wS1)),p.stemp(iWater(wS1))-disBT(i1231,iWater(wS1)),'bo',p.stemp(iWater(wS2)),p.stemp(iWater(wS2))-disBT(i1231,iWater(wS2)),'gx',p.stemp(iWater(wS3)),p.stemp(iWater(wS3))-disBT(i1231,iWater(wS3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iWater(wS1)),pclBT(i1231,iWater(wS1))-disBT(i1231,iWater(wS1)),'bo',p.stemp(iWater(wS2)),pclBT(i1231,iWater(wS2))-disBT(i1231,iWater(wS2)),'gx',p.stemp(iWater(wS3)),pclBT(i1231,iWater(wS3))-disBT(i1231,iWater(wS3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iWater(wS1)),pclBT(i1231,iWater(wS1))-disBT(i1231,iWater(wS1)),'bo',disBT(i1231,iWater(wS2)),pclBT(i1231,iWater(wS2))-disBT(i1231,iWater(wS2)),'gx',disBT(i1231,iWater(wS3)),pclBT(i1231,iWater(wS3))-disBT(i1231,iWater(wS3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topWater   = unique(p.cprtop(iWater));
topWater   = cprtop0;
wT1 = find(p.cprtop(iWater) <= topWater(3));
wT2 = find(p.cprtop(iWater) > topWater(3) & p.cprtop(iWater) <= topWater(6));
wT3 = find(p.cprtop(iWater) > topWater(6));

figure(2)
plot(dis.fc,nanmean(disBT(:,iWater(wT1))'-pcl0BT(:,iWater(wT1))'),'b',dis.fc,nanstd(disBT(:,iWater(wT1))'-pcl0BT(:,iWater(wT1))'),'c--',dis.fc,nanmean(disBT(:,iWater(wT1))'-pclBT(:,iWater(wT1))'),'r',dis.fc,nanstd(disBT(:,iWater(wT1))'-pclBT(:,iWater(wT1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld HighAlt')

figure(3)
plot(dis.fc,nanmean(disBT(:,iWater(wT2))'-pcl0BT(:,iWater(wT2))'),'b',dis.fc,nanstd(disBT(:,iWater(wT2))'-pcl0BT(:,iWater(wT2))'),'c--',dis.fc,nanmean(disBT(:,iWater(wT2))'-pclBT(:,iWater(wT2))'),'r',dis.fc,nanstd(disBT(:,iWater(wT2))'-pclBT(:,iWater(wT2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld MediumAlt')

figure(4)
plot(dis.fc,nanmean(disBT(:,iWater(wT3))'-pcl0BT(:,iWater(wT3))'),'b',dis.fc,nanstd(disBT(:,iWater(wT3))'-pcl0BT(:,iWater(wT3))'),'c--',dis.fc,nanmean(disBT(:,iWater(wT3))'-pclBT(:,iWater(wT3))'),'r',dis.fc,nanstd(disBT(:,iWater(wT3))'-pclBT(:,iWater(wT3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld LowAlt')

figure(5); 
plot(p.stemp(iWater(wT1)),disBT(i1231,iWater(wT1)),'bo',p.stemp(iWater(wT2)),disBT(i1231,iWater(wT2)),'gx',p.stemp(iWater(wT3)),disBT(i1231,iWater(wT3)),'r.',p.stemp(iWater(wT1)),p.stemp(iWater(wT1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iWater(wT1)),p.stemp(iWater(wT1))-disBT(i1231,iWater(wT1)),'bo',p.stemp(iWater(wT2)),p.stemp(iWater(wT2))-disBT(i1231,iWater(wT2)),'gx',p.stemp(iWater(wT3)),p.stemp(iWater(wT3))-disBT(i1231,iWater(wT3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iWater(wT1)),pclBT(i1231,iWater(wT1))-disBT(i1231,iWater(wT1)),'bo',p.stemp(iWater(wT2)),pclBT(i1231,iWater(wT2))-disBT(i1231,iWater(wT2)),'gx',p.stemp(iWater(wT3)),pclBT(i1231,iWater(wT3))-disBT(i1231,iWater(wT3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iWater(wT1)),pclBT(i1231,iWater(wT1))-disBT(i1231,iWater(wT1)),'bo',disBT(i1231,iWater(wT2)),pclBT(i1231,iWater(wT2))-disBT(i1231,iWater(wT2)),'gx',disBT(i1231,iWater(wT3)),pclBT(i1231,iWater(wT3))-disBT(i1231,iWater(wT3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amtWater   = unique(p.cngwat(iWater));
amtWater   = cngwat0;
wA1 = find(p.cngwat(iWater) <= amtWater(3));
wA2 = find(p.cngwat(iWater) > amtWater(3) & p.cngwat(iWater) <= amtWater(5));
wA3 = find(p.cngwat(iWater) > amtWater(5));

figure(2)
plot(dis.fc,nanmean(disBT(:,iWater(wA1))'-pcl0BT(:,iWater(wA1))'),'b',dis.fc,nanstd(disBT(:,iWater(wA1))'-pcl0BT(:,iWater(wA1))'),'c--',dis.fc,nanmean(disBT(:,iWater(wA1))'-pclBT(:,iWater(wA1))'),'r',dis.fc,nanstd(disBT(:,iWater(wA1))'-pclBT(:,iWater(wA1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld LowOD')

figure(3)
plot(dis.fc,nanmean(disBT(:,iWater(wA2))'-pcl0BT(:,iWater(wA2))'),'b',dis.fc,nanstd(disBT(:,iWater(wA2))'-pcl0BT(:,iWater(wA2))'),'c--',dis.fc,nanmean(disBT(:,iWater(wA2))'-pclBT(:,iWater(wA2))'),'r',dis.fc,nanstd(disBT(:,iWater(wA2))'-pclBT(:,iWater(wA2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld MediumOD')

figure(4)
plot(dis.fc,nanmean(disBT(:,iWater(wA3))'-pcl0BT(:,iWater(wA3))'),'b',dis.fc,nanstd(disBT(:,iWater(wA3))'-pcl0BT(:,iWater(wA3))'),'c--',dis.fc,nanmean(disBT(:,iWater(wA3))'-pclBT(:,iWater(wA3))'),'r',dis.fc,nanstd(disBT(:,iWater(wA3))'-pclBT(:,iWater(wA3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld HighOD')

figure(5); 
plot(p.stemp(iWater(wA1)),disBT(i1231,iWater(wA1)),'bo',p.stemp(iWater(wA2)),disBT(i1231,iWater(wA2)),'gx',p.stemp(iWater(wA3)),disBT(i1231,iWater(wA3)),'r.',p.stemp(iWater(wA1)),p.stemp(iWater(wA1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iWater(wA1)),p.stemp(iWater(wA1))-disBT(i1231,iWater(wA1)),'bo',p.stemp(iWater(wA2)),p.stemp(iWater(wA2))-disBT(i1231,iWater(wA2)),'gx',p.stemp(iWater(wA3)),p.stemp(iWater(wA3))-disBT(i1231,iWater(wA3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iWater(wA1)),pclBT(i1231,iWater(wA1))-disBT(i1231,iWater(wA1)),'bo',p.stemp(iWater(wA2)),pclBT(i1231,iWater(wA2))-disBT(i1231,iWater(wA2)),'gx',p.stemp(iWater(wA3)),pclBT(i1231,iWater(wA3))-disBT(i1231,iWater(wA3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iWater(wA1)),pclBT(i1231,iWater(wA1))-disBT(i1231,iWater(wA1)),'bo',disBT(i1231,iWater(wA2)),pclBT(i1231,iWater(wA2))-disBT(i1231,iWater(wA2)),'gx',disBT(i1231,iWater(wA3)),pclBT(i1231,iWater(wA3))-disBT(i1231,iWater(wA3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angleWater = unique(p.scanang(iWater));
angleWater = scanang0;  
wX1 = find(p.scanang(iWater) <= 10);
wX2 = find(p.scanang(iWater) > 10 & p.scanang(iWater) <= 30);
wX3 = find(p.scanang(iWater) > 30);

figure(2)
plot(dis.fc,nanmean(disBT(:,iWater(wX1))'-pcl0BT(:,iWater(wX1))'),'b',dis.fc,nanstd(disBT(:,iWater(wX1))'-pcl0BT(:,iWater(wX1))'),'c--',dis.fc,nanmean(disBT(:,iWater(wX1))'-pclBT(:,iWater(wX1))'),'r',dis.fc,nanstd(disBT(:,iWater(wX1))'-pclBT(:,iWater(wX1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld LowAngle')

figure(3)
plot(dis.fc,nanmean(disBT(:,iWater(wX2))'-pcl0BT(:,iWater(wX2))'),'b',dis.fc,nanstd(disBT(:,iWater(wX2))'-pcl0BT(:,iWater(wX2))'),'c--',dis.fc,nanmean(disBT(:,iWater(wX2))'-pclBT(:,iWater(wX2))'),'r',dis.fc,nanstd(disBT(:,iWater(wX2))'-pclBT(:,iWater(wX2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld MediumAngle')

figure(4)
plot(dis.fc,nanmean(disBT(:,iWater(wX3))'-pcl0BT(:,iWater(wX3))'),'b',dis.fc,nanstd(disBT(:,iWater(wX3))'-pcl0BT(:,iWater(wX3))'),'c--',dis.fc,nanmean(disBT(:,iWater(wX3))'-pclBT(:,iWater(wX3))'),'r',dis.fc,nanstd(disBT(:,iWater(wX3))'-pclBT(:,iWater(wX3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Water Cld HighAngle')

figure(5); 
plot(p.stemp(iWater(wX1)),disBT(i1231,iWater(wX1)),'bo',p.stemp(iWater(wX2)),disBT(i1231,iWater(wX2)),'gx',p.stemp(iWater(wX3)),disBT(i1231,iWater(wX3)),'r.',p.stemp(iWater(wX1)),p.stemp(iWater(wX1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iWater(wX1)),p.stemp(iWater(wX1))-disBT(i1231,iWater(wX1)),'bo',p.stemp(iWater(wX2)),p.stemp(iWater(wX2))-disBT(i1231,iWater(wX2)),'gx',p.stemp(iWater(wX3)),p.stemp(iWater(wX3))-disBT(i1231,iWater(wX3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iWater(wX1)),pclBT(i1231,iWater(wX1))-disBT(i1231,iWater(wX1)),'bo',p.stemp(iWater(wX2)),pclBT(i1231,iWater(wX2))-disBT(i1231,iWater(wX2)),'gx',p.stemp(iWater(wX3)),pclBT(i1231,iWater(wX3))-disBT(i1231,iWater(wX3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iWater(wX1)),pclBT(i1231,iWater(wX1))-disBT(i1231,iWater(wX1)),'bo',disBT(i1231,iWater(wX2)),pclBT(i1231,iWater(wX2))-disBT(i1231,iWater(wX2)),'gx',disBT(i1231,iWater(wX3)),pclBT(i1231,iWater(wX3))-disBT(i1231,iWater(wX3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(dis.fc,nanmean(disBT(:,iIce)'-pcl0BT(:,iIce)'),'b',dis.fc,nanstd(disBT(:,iIce)'-pcl0BT(:,iIce)'),'c--',dis.fc,nanmean(disBT(:,iIce)'-pclBT(:,iIce)'),'r',dis.fc,nanstd(disBT(:,iIce)'-pclBT(:,iIce)'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld')

szeIce   = unique(p.cpsize(iIce));
szeIce   = ice_cpsize0;
iS1 = find(p.cpsize(iIce) <= szeIce(2));
iS2 = find(p.cpsize(iIce) > szeIce(2) & p.cpsize(iIce) <= szeIce(4));
iS3 = find(p.cpsize(iIce) > szeIce(4));

figure(2)
plot(dis.fc,nanmean(disBT(:,iIce(iS1))'-pcl0BT(:,iIce(iS1))'),'b',dis.fc,nanstd(disBT(:,iIce(iS1))'-pcl0BT(:,iIce(iS1))'),'c--',dis.fc,nanmean(disBT(:,iIce(iS1))'-pclBT(:,iIce(iS1))'),'r',dis.fc,nanstd(disBT(:,iIce(iS1))'-pclBT(:,iIce(iS1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld SmallDme')

figure(3)
plot(dis.fc,nanmean(disBT(:,iIce(iS2))'-pcl0BT(:,iIce(iS2))'),'b',dis.fc,nanstd(disBT(:,iIce(iS2))'-pcl0BT(:,iIce(iS2))'),'c--',dis.fc,nanmean(disBT(:,iIce(iS2))'-pclBT(:,iIce(iS2))'),'r',dis.fc,nanstd(disBT(:,iIce(iS2))'-pclBT(:,iIce(iS2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld MediumDME')

figure(4)
plot(dis.fc,nanmean(disBT(:,iIce(iS3))'-pcl0BT(:,iIce(iS3))'),'b',dis.fc,nanstd(disBT(:,iIce(iS3))'-pcl0BT(:,iIce(iS3))'),'c--',dis.fc,nanmean(disBT(:,iIce(iS3))'-pclBT(:,iIce(iS3))'),'r',dis.fc,nanstd(disBT(:,iIce(iS3))'-pclBT(:,iIce(iS3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld LargeDME ')

figure(5); 
plot(p.stemp(iIce(iS1)),disBT(i1231,iIce(iS1)),'bo',p.stemp(iIce(iS2)),disBT(i1231,iIce(iS2)),'gx',p.stemp(iIce(iS3)),disBT(i1231,iIce(iS3)),'r.',p.stemp(iIce(iS1)),p.stemp(iIce(iS1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iIce(iS1)),p.stemp(iIce(iS1))-disBT(i1231,iIce(iS1)),'bo',p.stemp(iIce(iS2)),p.stemp(iIce(iS2))-disBT(i1231,iIce(iS2)),'gx',p.stemp(iIce(iS3)),p.stemp(iIce(iS3))-disBT(i1231,iIce(iS3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iIce(iS1)),pclBT(i1231,iIce(iS1))-disBT(i1231,iIce(iS1)),'bo',p.stemp(iIce(iS2)),pclBT(i1231,iIce(iS2))-disBT(i1231,iIce(iS2)),'gx',p.stemp(iIce(iS3)),pclBT(i1231,iIce(iS3))-disBT(i1231,iIce(iS3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iIce(iS1)),pclBT(i1231,iIce(iS1))-disBT(i1231,iIce(iS1)),'bo',disBT(i1231,iIce(iS2)),pclBT(i1231,iIce(iS2))-disBT(i1231,iIce(iS2)),'gx',disBT(i1231,iIce(iS3)),pclBT(i1231,iIce(iS3))-disBT(i1231,iIce(iS3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

topIce   = unique(p.cprtop(iIce));
topIce   = cprtop0;
iT1 = find(p.cprtop(iIce) <= topIce(3));
iT2 = find(p.cprtop(iIce) > topIce(3) & p.cprtop(iIce) <= topIce(6));
iT3 = find(p.cprtop(iIce) > topIce(6));

figure(2)
plot(dis.fc,nanmean(disBT(:,iIce(iT1))'-pcl0BT(:,iIce(iT1))'),'b',dis.fc,nanstd(disBT(:,iIce(iT1))'-pcl0BT(:,iIce(iT1))'),'c--',dis.fc,nanmean(disBT(:,iIce(iT1))'-pclBT(:,iIce(iT1))'),'r',dis.fc,nanstd(disBT(:,iIce(iT1))'-pclBT(:,iIce(iT1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld HighAlt')

figure(3)
plot(dis.fc,nanmean(disBT(:,iIce(iT2))'-pcl0BT(:,iIce(iT2))'),'b',dis.fc,nanstd(disBT(:,iIce(iT2))'-pcl0BT(:,iIce(iT2))'),'c--',dis.fc,nanmean(disBT(:,iIce(iT2))'-pclBT(:,iIce(iT2))'),'r',dis.fc,nanstd(disBT(:,iIce(iT2))'-pclBT(:,iIce(iT2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld MediumAlt')

figure(4)
plot(dis.fc,nanmean(disBT(:,iIce(iT3))'-pcl0BT(:,iIce(iT3))'),'b',dis.fc,nanstd(disBT(:,iIce(iT3))'-pcl0BT(:,iIce(iT3))'),'c--',dis.fc,nanmean(disBT(:,iIce(iT3))'-pclBT(:,iIce(iT3))'),'r',dis.fc,nanstd(disBT(:,iIce(iT3))'-pclBT(:,iIce(iT3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld LowAlt')

figure(5); 
plot(p.stemp(iIce(iT1)),disBT(i1231,iIce(iT1)),'bo',p.stemp(iIce(iT2)),disBT(i1231,iIce(iT2)),'gx',p.stemp(iIce(iT3)),disBT(i1231,iIce(iT3)),'r.',p.stemp(iIce(iT1)),p.stemp(iIce(iT1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iIce(iT1)),p.stemp(iIce(iT1))-disBT(i1231,iIce(iT1)),'bo',p.stemp(iIce(iT2)),p.stemp(iIce(iT2))-disBT(i1231,iIce(iT2)),'gx',p.stemp(iIce(iT3)),p.stemp(iIce(iT3))-disBT(i1231,iIce(iT3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iIce(iT1)),pclBT(i1231,iIce(iT1))-disBT(i1231,iIce(iT1)),'bo',p.stemp(iIce(iT2)),pclBT(i1231,iIce(iT2))-disBT(i1231,iIce(iT2)),'gx',p.stemp(iIce(iT3)),pclBT(i1231,iIce(iT3))-disBT(i1231,iIce(iT3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iIce(iT1)),pclBT(i1231,iIce(iT1))-disBT(i1231,iIce(iT1)),'bo',disBT(i1231,iIce(iT2)),pclBT(i1231,iIce(iT2))-disBT(i1231,iIce(iT2)),'gx',disBT(i1231,iIce(iT3)),pclBT(i1231,iIce(iT3))-disBT(i1231,iIce(iT3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amtIce   = unique(p.cngwat(iIce));
amtIce   = cngwat0;
iA1 = find(p.cngwat(iIce) <= amtIce(3));
iA2 = find(p.cngwat(iIce) > amtIce(3) & p.cngwat(iIce) <= amtIce(5));
iA3 = find(p.cngwat(iIce) > amtIce(5));

figure(2)
plot(dis.fc,nanmean(disBT(:,iIce(iA1))'-pcl0BT(:,iIce(iA1))'),'b',dis.fc,nanstd(disBT(:,iIce(iA1))'-pcl0BT(:,iIce(iA1))'),'c--',dis.fc,nanmean(disBT(:,iIce(iA1))'-pclBT(:,iIce(iA1))'),'r',dis.fc,nanstd(disBT(:,iIce(iA1))'-pclBT(:,iIce(iA1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld LowOD')

figure(3)
plot(dis.fc,nanmean(disBT(:,iIce(iA2))'-pcl0BT(:,iIce(iA2))'),'b',dis.fc,nanstd(disBT(:,iIce(iA2))'-pcl0BT(:,iIce(iA2))'),'c--',dis.fc,nanmean(disBT(:,iIce(iA2))'-pclBT(:,iIce(iA2))'),'r',dis.fc,nanstd(disBT(:,iIce(iA2))'-pclBT(:,iIce(iA2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld MediumOD')

figure(4)
plot(dis.fc,nanmean(disBT(:,iIce(iA3))'-pcl0BT(:,iIce(iA3))'),'b',dis.fc,nanstd(disBT(:,iIce(iA3))'-pcl0BT(:,iIce(iA3))'),'c--',dis.fc,nanmean(disBT(:,iIce(iA3))'-pclBT(:,iIce(iA3))'),'r',dis.fc,nanstd(disBT(:,iIce(iA3))'-pclBT(:,iIce(iA3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld HighOD')

figure(5); 
plot(p.stemp(iIce(iA1)),disBT(i1231,iIce(iA1)),'bo',p.stemp(iIce(iA2)),disBT(i1231,iIce(iA2)),'gx',p.stemp(iIce(iA3)),disBT(i1231,iIce(iA3)),'r.',p.stemp(iIce(iA1)),p.stemp(iIce(iA1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iIce(iA1)),p.stemp(iIce(iA1))-disBT(i1231,iIce(iA1)),'bo',p.stemp(iIce(iA2)),p.stemp(iIce(iA2))-disBT(i1231,iIce(iA2)),'gx',p.stemp(iIce(iA3)),p.stemp(iIce(iA3))-disBT(i1231,iIce(iA3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iIce(iA1)),pclBT(i1231,iIce(iA1))-disBT(i1231,iIce(iA1)),'bo',p.stemp(iIce(iA2)),pclBT(i1231,iIce(iA2))-disBT(i1231,iIce(iA2)),'gx',p.stemp(iIce(iA3)),pclBT(i1231,iIce(iA3))-disBT(i1231,iIce(iA3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iIce(iA1)),pclBT(i1231,iIce(iA1))-disBT(i1231,iIce(iA1)),'bo',disBT(i1231,iIce(iA2)),pclBT(i1231,iIce(iA2))-disBT(i1231,iIce(iA2)),'gx',disBT(i1231,iIce(iA3)),pclBT(i1231,iIce(iA3))-disBT(i1231,iIce(iA3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angleIce = unique(p.scanang(iIce));
angleIce = scanang0;
iX1 = find(p.scanang(iIce) <= 10);
iX2 = find(p.scanang(iIce) > 10 & p.scanang(iIce) <= 30);
iX3 = find(p.scanang(iIce) > 30);

figure(2)
plot(dis.fc,nanmean(disBT(:,iIce(iX1))'-pcl0BT(:,iIce(iX1))'),'b',dis.fc,nanstd(disBT(:,iIce(iX1))'-pcl0BT(:,iIce(iX1))'),'c--',dis.fc,nanmean(disBT(:,iIce(iX1))'-pclBT(:,iIce(iX1))'),'r',dis.fc,nanstd(disBT(:,iIce(iX1))'-pclBT(:,iIce(iX1))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld LowAngle')

figure(3)
plot(dis.fc,nanmean(disBT(:,iIce(iX2))'-pcl0BT(:,iIce(iX2))'),'b',dis.fc,nanstd(disBT(:,iIce(iX2))'-pcl0BT(:,iIce(iX2))'),'c--',dis.fc,nanmean(disBT(:,iIce(iX2))'-pclBT(:,iIce(iX2))'),'r',dis.fc,nanstd(disBT(:,iIce(iX2))'-pclBT(:,iIce(iX2))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld MediumAngle')

figure(4)
plot(dis.fc,nanmean(disBT(:,iIce(iX3))'-pcl0BT(:,iIce(iX3))'),'b',dis.fc,nanstd(disBT(:,iIce(iX3))'-pcl0BT(:,iIce(iX3))'),'c--',dis.fc,nanmean(disBT(:,iIce(iX3))'-pclBT(:,iIce(iX3))'),'r',dis.fc,nanstd(disBT(:,iIce(iX3))'-pclBT(:,iIce(iX3))'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','New bias','New std','fontsize',8,'location','best'); title('Ice Cld HighAngle')

figure(5); 
plot(p.stemp(iIce(iX1)),disBT(i1231,iIce(iX1)),'bo',p.stemp(iIce(iX2)),disBT(i1231,iIce(iX2)),'gx',p.stemp(iIce(iX3)),disBT(i1231,iIce(iX3)),'r.',p.stemp(iIce(iX1)),p.stemp(iIce(iX1)),'k-')
  xlabel('stemp'); ylabel('BT1231 DISORT');
plot(p.stemp(iIce(iX1)),p.stemp(iIce(iX1))-disBT(i1231,iIce(iX1)),'bo',p.stemp(iIce(iX2)),p.stemp(iIce(iX2))-disBT(i1231,iIce(iX2)),'gx',p.stemp(iIce(iX3)),p.stemp(iIce(iX3))-disBT(i1231,iIce(iX3)),'r.')
  xlabel('stemp'); ylabel('stemp-BT1231 DISORT');
plot(p.stemp(iIce(iX1)),pclBT(i1231,iIce(iX1))-disBT(i1231,iIce(iX1)),'bo',p.stemp(iIce(iX2)),pclBT(i1231,iIce(iX2))-disBT(i1231,iIce(iX2)),'gx',p.stemp(iIce(iX3)),pclBT(i1231,iIce(iX3))-disBT(i1231,iIce(iX3)),'r.')
  xlabel('stemp'); ylabel('BT1231 PCLSAM-DISORT');
plot(disBT(i1231,iIce(iX1)),pclBT(i1231,iIce(iX1))-disBT(i1231,iIce(iX1)),'bo',disBT(i1231,iIce(iX2)),pclBT(i1231,iIce(iX2))-disBT(i1231,iIce(iX2)),'gx',disBT(i1231,iIce(iX3)),pclBT(i1231,iIce(iX3))-disBT(i1231,iIce(iX3)),'r.')
  xlabel('BT1231 DISORT'); ylabel('BT1231 PCLSAM-DISORT');

disp('ret to continue'); pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
