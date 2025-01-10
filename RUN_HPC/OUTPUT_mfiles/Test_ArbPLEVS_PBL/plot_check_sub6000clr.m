%% see /asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/make_clear_rtp_039.m, copied here

%% done by Chris??????? or Sergio?????
[hC,ha,psarta_clr_rad_airs,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_sartaH2020_rp.rtp'); %% H2020
tsarta_clr_rad_airs0 = rad2bt(hC.vchan,psarta_clr_rad_airs.rcalc);
if ~exist('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_sartaH2020_chrisexec_rp.rtp')
  cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev fin=sub6000clr_airs_g12_sartaH2020_rp.rtp fout=sub6000clr_airs_g12_sartaH2020_chrisexec_rp.rtp'];
  eval(sartaer);
  cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/Test_ArbPLEVS_PBL/
end
[hC,ha,psarta_clr_rad_airsX,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_sartaH2020_chrisexec_rp.rtp'); %% H2020
tsarta_clr_rad_airsX = rad2bt(hC.vchan,psarta_clr_rad_airsX.rcalc);

tsarta_clr_rad_airs = tsarta_clr_rad_airsX;

%% done by Sergio
[hS,ha,psarta_clr_rad_airs2,pa] = rtpread('/home/chepplew/projects/klayers_wrk/sergio_sub6000clr_airs_g12_sar_may19.rtp');                                                     %% H2016 WIERD CHANS
[hS,ha,psarta_clr_rad_airs2,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_rp_airs_l1c_2834_cloudy_may19.rtp');  %% H2016
[hS,ha,psarta_clr_rad_airs2,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_rp.rtp');                             %% H2016
tsarta_clr_rad_airs2 = rad2bt(hS.vchan,psarta_clr_rad_airs2.rcalc);

scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.spres); colorbar
scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.satzen); colorbar
scatter_coast(psarta_clr_rad_airs.rlon,psarta_clr_rad_airs.rlat,10,psarta_clr_rad_airs.stemp-tsarta_clr_rad_airs(1520,:)); colorbar

h = hS;

plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tsarta_clr_rad_airs2'),h.vchan,nanstd(tsarta_clr_rad_airs'-tsarta_clr_rad_airs2'))
title('SARTA : H2020 - H2016'); legend('bias','std dev','location','best');
disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nprof = length(psarta_clr_rad_airs.stemp);
p = psarta_clr_rad_airs;
  landfrac = p.landfrac;
  day   = find(p.solzen < 90);
  night = find(p.solzen >= 90);

  day   = find(p.solzen < 90 & landfrac == 0);
  night = find(p.solzen >= 90 & landfrac == 0);
  %night = find(p.solzen >= 135 & landfrac == 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');
ichan = h2645.h.ichan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a12 = load('6555_profiles/AIRS_G12/airs_G12_summary.mat');
p12 = load('6555_profiles/PBL_G12/pbl_G12_summary.mat');
a80 = load('6555_profiles/AIRS_G80/airs_G80_summary.mat');
p80 = load('6555_profiles/PBL_G80/pbl_G80_summary.mat');

fKc = a12.fKc;
tpbl12 = rad2bt(fKc(ichan),p12.pbl_kcarta(ichan,:));
tpbl80 = rad2bt(fKc(ichan),p80.pbl_kcarta(ichan,:));
tairs12 = rad2bt(fKc(ichan),a12.airs_kcarta(ichan,:));
tairs80 = rad2bt(fKc(ichan),a80.airs_kcarta(ichan,:));

figure(1);
plot(fKc(ichan),nanmean(tairs12'-tairs80'),'b.-',fKc(ichan),nanmean(tairs12'-tpbl80'),'r',fKc(ichan),nanmean(tairs12'-tpbl12'),'k')
plotaxis2; legend('airs12-airs80','airs12-pbl80','airs12-pbl12','location','best','fontsize',10)
title('Bias');

figure(2)
plot(fKc(ichan),nanstd(tairs12'-tairs80'),'b.-',fKc(ichan),nanstd(tairs12'-tpbl80'),'r',fKc(ichan),nanstd(tairs12'-tpbl12'),'k')
plotaxis2; legend('airs12-airs80','airs12-pbl80','airs12-pbl12','location','best','fontsize',10)
title('St Dev')

disp('ret to continue'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Y,I] = sort(h.vchan);
h.vchanSort = Y;

[Y,I1,I2] = intersect(h.vchan,fKc(ichan));
for ii = 1 : length(ichan)
  junk = abs(h.vchan-fKc(ichan(ii)));
  I(ii) = find(junk == min(junk),1);
  chandiff(ii) = junk(I(ii));
end
h.vchanSort = h.vchan(I);

plot(h.vchanSort,h.vchanSort-fKc(ichan))
plot(h.vchanSort,chandiff); title('ChanDiff (sorted) \newline kCARTA-SARTA')
pause(1);

plot(h.vchan,nanmean(tsarta_clr_rad_airs'),'b.-',h.vchanSort,nanmean(tsarta_clr_rad_airs(I,:)'),'r',...
     fKc(ichan),nanmean(tairs12'),'k')
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,night)'),'b.-',h.vchanSort,nanmean(tsarta_clr_rad_airs(I,night)'),'r',...
     fKc(ichan),nanmean(tairs12(:,night)'),'k')
pause(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tsarta_clr_rad_airs = tsarta_clr_rad_airs(I,:);

%tpbl12 = tpbl12(I,:);
%tpbl80 = tpbl80(I,:);
%tairs12 = tairs12(I,:);
%tairs80 = tairs80(I,:);

tpbl = tpbl80;
tairs = tairs80;
figure(1); 
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
  title('G80 PBL');

figure(2);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
  title('G80 AIRS')

figure(3);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',...
     h.vchanSort,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('G80')

figure(4); 
plot(h.vchanSort,nanmean(tairs'-tpbl'),'k',h.vchanSort,nanstd(tairs'-tpbl'),'g'); plotaxis2;
  title('G80 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tpbl = tpbl12;
tairs = tairs12;
figure(1); 
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
  title('G12 PBL');

figure(2);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
  title('G12 AIRS')

figure(3);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',...
     h.vchanSort,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('G12')

figure(4); 
plot(h.vchanSort,nanmean(tairs'-tpbl'),'k',h.vchanSort,nanstd(tairs'-tpbl'),'g'); plotaxis2;
  title('G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ix = day;

figure(1); 
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c'); plotaxis2;
  title('D G12 PBL');
 ylabel('SARTA-KCARTA (b) mean (c) std')

figure(2);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
  title('D G12 AIRS')
 ylabel('SARTA-KCARTA (r) mean (m) std')

figure(3);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c',...
     h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('D G12')

figure(4); 
plot(h.vchanSort,nanmean(tairs(:,ix)'-tpbl(:,ix)'),'k',h.vchanSort,nanstd(tairs(:,ix)'-tpbl(:,ix)'),'g'); plotaxis2;
  title('D G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ix = night;

figure(1); 
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c'); plotaxis2;
  title('N G12 PBL');

figure(2);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
  title('N G12 AIRS')

figure(3);
plot(h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c',...
     h.vchanSort,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchanSort,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('N G12')

figure(4); 
plot(h.vchanSort,nanmean(tairs(:,ix)'-tpbl(:,ix)'),'k',h.vchanSort,nanstd(tairs(:,ix)'-tpbl(:,ix)'),'g'); plotaxis2;
  title('N G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
raLat = -90:5:+90; 
for ii = 1 : length(raLat)-1
  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1));
  thenum(ii) = length(oo);
end
plot(meanvaluebin(raLat),thenum)
[~,~,pjunk,~] = rtpread('/home/chepplew/data/scratch/sub6000clr_pbl_g12_op.rtp');
plot(pjunk.rlon-p.rlon)

i2326 = find(h.vchanSort >= 2326,1);
figure(5); clf; scatter_coast(p.rlon,p.rlat,10,tsarta_clr_rad_airs(i2326,:) - tairs(i2326,:)); colormap jet
title('BT2326 Sarta Cal - kCarta Cal')

figure(6); clf
for ii = 1 : length(raLat)-1
  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1));
  four_um_bias(ii) = -nanmean(tsarta_clr_rad_airs(i2326,oo) - tairs(i2326,oo));
  four_um_solzen(ii) = nanmean(p.solzen(oo));

  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1) & p.landfrac == 0 & p.solzen < 90);
  nltebias(ii) = -nanmean(tsarta_clr_rad_airs(i2326,oo) - tairs(i2326,oo));
  nlte_solzen(ii) = nanmean(p.solzen(oo));

  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1) & p.landfrac == 0 & p.solzen > 90);
  ltebias(ii) = -nanmean(tsarta_clr_rad_airs(i2326,oo) - tairs(i2326,oo));
  lte_solzen(ii) = nanmean(p.solzen(oo));
end

plot(meanvaluebin(raLat),nltebias,'r',meanvaluebin(raLat),ltebias,'b',meanvaluebin(raLat),four_um_bias,'k','linewidth',2)
  xlabel('Latitude'); ylabel('KCARTA-SARTA'); legend('day/ocean','night/ocean','D/N/L/O=ALL','location','best','fontsize',10)
plot(meanvaluebin(raLat),nlte_solzen,'r',meanvaluebin(raLat),lte_solzen,'b',meanvaluebin(raLat),four_um_solzen,'k','linewidth',2)
  xlabel('Latitude'); ylabel('SOLZEN'); legend('day/ocean','night/ocean','D/N/L/O=ALL','location','best','fontsize',10)

plot(p.solzen,tairs(i2326,:)-tsarta_clr_rad_airs(i2326,:),'.')
xlabel('Solzen'); ylabel('KCARTA-SARTA'); title('BTD at 2326 cm-1')
plot(p.solzen,tairs(i2326,:),'rx',p.solzen,tsarta_clr_rad_airs(i2326,:),'b.')
xlabel('Solzen'); ylabel('KCARTA(r) SARTA(b)'); title('BT at 2326 cm-1')


i1039 = find(h.vchanSort >= 1039,1);
figure(5); clf; scatter_coast(p.rlon,p.rlat,10,tsarta_clr_rad_airs(i1039,:) - tairs(i1039,:)); colormap jet
title('BT1039 Sarta Cal - kCarta Cal')
for ii = 1 : length(raLat)-1
  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1));
  o3bias(ii) = nanmean(tsarta_clr_rad_airs(i1039,oo) - tairs(i1039,oo));
end
plot(meanvaluebin(raLat),o3bias)

plot(meanvaluebin(raLat),o3bias,meanvaluebin(raLat),nltebias); ylabel('SARTA-KCARTA')
plotaxis2;
  legend('ozone','nlte','location','best','fontsize',10); 
xlabel('Latitude'); ylabel('SARTA-KCARTA bias [K]');


%{
mmw = mmwater_rtp(h,p);
o3du = dobson_rtp(h,p);
plot(p.rlat,o3du,'.'); xlabel('Latitude'); ylabel('O3 dobson units');
%}

%{
[h49,~,p49,~] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/regr49_airs.op.rtp');
[h49,~,p49x,~] = rtpread('/asl/s1/sergio/RTP_pin_feb2002/pin_feb2002_sea_airsnadir_op.so2.latlon.rtp');
[h49,~,p49,~] = rtpread('/asl/packages/klayersV205/Data/adafgl_16Aug2010_op.rtp');
[h49,~,p49,~] = rtpread('/asl/packages/klayersV205/Test_rtpV201/pin_feb2002_sea_airsnadir_op.rtp');
p49.rlat = p49x.rlat;
mmw49 = mmwater_rtp(h49,p49);
o3du49 = dobson_rtp(h49,p49);
plot(p.rlat,o3du,'b.',p49.rlat,o3du49,'r.'); xlabel('Latitude'); ylabel('O3 dobson units');
plot(p.stemp,o3du,'b.',p49.stemp,o3du49,'r.'); xlabel('SKT [K]'); ylabel('O3 dobson units');
plot(p.stemp,mmw,'b.',p49.stemp,mmw49,'r.'); xlabel('SKT [K]'); ylabel('Col WV');
loglog(p.gas_3,p.plevs,'b',p49.gas_3,p49.plevs,'r')

loglog(p.gas_3,p.plevs,'b',p49.gas_3,p49.plevs,'r');
xlim([1e13 1e18]);title('b=ECM  r=Regr')
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

[hxx,ha,pS_a12,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sergio_sub6000clr_airs.op.rtp');
[hxx,ha,pS_p12,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sergio_sub6000clr_pbl.op.rtp');
[hxx,ha,pS_a80,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sergio_sub6000clr_airs_g80.op.rtp');
[hxx,ha,pS_p80,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sergio_sub6000clr_pbl_g80.op.rtp');

[hxx,ha,pC_a12,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_rp.rtp');
[hxx,ha,pC_p12,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_pbl_g12_op.rtp');
[hxx,ha,pC_a80,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g80_op.rtp');
[hxx,ha,pC_p80,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_pbl_g80_op.rtp');

pS=pS_a12; pC=pC_a12; compare_two_structures(pS,pC);
pS=pS_a80; pC=pC_a80; compare_two_structures(pS,pC);

pS=pS_p12; pC=pC_p12; compare_two_structures(pS,pC);
pS=pS_p80; pC=pC_p80; compare_two_structures(pS,pC);

l2s = load('/home/sergio/KCARTA/L2SComparisons/l2s_kc122_H20_605_2830.mat');
l2s_gid = l2s.iaGasID;
[~,l2s] = convolve_airs(double(l2s.w),double(l2s.d),double(ichan));

%for ig = 1 : 72; 
for ig = 3; 
  plot(h.vchanSort,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchanSort,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',h.vchanSort,exp(-l2s(:,ig))); 
  plotaxis2; title(num2str(l2s_gid(ig)));
  pause
end
%}
