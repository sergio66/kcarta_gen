%% see /asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/make_clear_rtp_039.m, copied here

[h,ha,psarta_clr_rad_airs,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_rp.rtp');
tsarta_clr_rad_airs = rad2bt(h.vchan,psarta_clr_rad_airs.rcalc);

scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.spres); colorbar
scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.satzen); colorbar
scatter_coast(psarta_clr_rad_airs.rlon,psarta_clr_rad_airs.rlat,10,psarta_clr_rad_airs.stemp-tsarta_clr_rad_airs(1520,:)); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nprof = length(psarta_clr_rad_airs.stemp);
p = psarta_clr_rad_airs;
  landfrac = p.landfrac;
  day   = find(p.solzen < 90);
  night = find(p.solzen >= 90);

  day   = find(p.solzen < 90 & landfrac == 0);
  night = find(p.solzen >= 90 & landfrac == 0);

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

plot(fKc(ichan),nanmean(tairs12'-tairs80'),'b.-',fKc(ichan),nanmean(tairs12'-tpbl80'),'r',fKc(ichan),nanmean(tairs12'-tpbl12'),'k')
plotaxis2; legend('airs12-airs80','airs12-pbl80','airs12-pbl12','location','best','fontsize',10)

plot(fKc(ichan),nanstd(tairs12'-tairs80'),'b.-',fKc(ichan),nanstd(tairs12'-tpbl80'),'r',fKc(ichan),nanstd(tairs12'-tpbl12'),'k')
plotaxis2; legend('airs12-airs80','airs12-pbl80','airs12-pbl12','location','best','fontsize',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tpbl = tpbl80;
tairs = tairs80;
figure(1); 
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
  title('G80 PBL');

figure(2);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
  title('G80 AIRS')

figure(3);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('G80')

figure(4); 
plot(h.vchan,nanmean(tairs'-tpbl'),'k',h.vchan,nanstd(tairs'-tpbl'),'g'); plotaxis2;
  title('G80 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tpbl = tpbl12;
tairs = tairs12;
figure(1); 
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
  title('G12 PBL');

figure(2);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
  title('G12 AIRS')

figure(3);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('G12')

figure(4); 
plot(h.vchan,nanmean(tairs'-tpbl'),'k',h.vchan,nanstd(tairs'-tpbl'),'g'); plotaxis2;
  title('G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ix = day;

figure(1); 
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c'); plotaxis2;
  title('D G12 PBL');

figure(2);
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
  title('D G12 AIRS')

figure(3);
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('D G12')

figure(4); 
plot(h.vchan,nanmean(tairs(:,ix)'-tpbl(:,ix)'),'k',h.vchan,nanstd(tairs(:,ix)'-tpbl(:,ix)'),'g'); plotaxis2;
  title('D G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ix = night;

figure(1); 
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c'); plotaxis2;
  title('N G12 PBL');

figure(2);
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
  title('N G12 AIRS')

figure(3);
plot(h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tpbl(:,ix)'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,ix)'-tairs(:,ix)'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
title('N G12')

figure(4); 
plot(h.vchan,nanmean(tairs(:,ix)'-tpbl(:,ix)'),'k',h.vchan,nanstd(tairs(:,ix)'-tpbl(:,ix)'),'g'); plotaxis2;
  title('N G12 : AIRS100 - PBL');

disp('ret to continue'); pause; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i2326 = find(h.vchan >= 2626,1);
figure(5); clf; scatter_coast(p.rlon,p.rlat,10,tsarta_clr_rad_airs(i2326,:) - tairs(i2326,:)); colormap jet
title('BT2326 Sarta Cal - kCarta Cal')

raLat = -90:5:+90; 
for ii = 1 : length(raLat)-1
  oo = find(p.rlat > raLat(ii) & p.rlat <= raLat(ii+1));
  thenum(ii) = length(oo);
end
plot(meanvaluebin(raLat),thenum)
[~,~,pjunk,~] = rtpread('/home/chepplew/data/scratch/sub6000clr_pbl_g12_op.rtp');
plot(pjunk.rlon-p.rlon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
l2s = load('/home/sergio/KCARTA/L2SComparisons/l2s_kc122_H20_605_2830.mat');
l2s_gid = l2s.iaGasID;
[~,l2s] = convolve_airs(double(l2s.w),double(l2s.d),double(ichan));

%for ig = 1 : 72; 
for ig = 3; 
  plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',h.vchan,exp(-l2s(:,ig))); 
  plotaxis2; title(num2str(l2s_gid(ig)));
  pause
end
%}
