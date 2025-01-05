sarta = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';

sartaer = ['!' sarta ' fin=summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_zonalavg64lats.rtp fout=junk.rad.rtp listp=32 listj=-1']
sartaer = ['!' sarta ' fin=junk.op.rtp fout=junk.rad.rtp listp=32 listj=-1']

addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools

[h,ha,p,pa] = rtpread('junk.op.rtp');

[w,jacT,iaProf,iaNumLay] = readsarta_jac('junk.rad.rtp_jacTZ',100);
[w,jac1,iaProf,iaNumLay] = readsarta_jac('junk.rad.rtp_jacG1',1);
[w,jac3,iaProf,iaNumLay] = readsarta_jac('junk.rad.rtp_jacG3',3);
[w,jacWGT,iaProf,iaNumLay] = readsarta_jac('junk.rad.rtp_WGTFCN',200);

figure(1); clf; pcolor(w,1:096,squeeze(jacT(1,:,1:096))');   title('Tz SARTA');   shading flat; set(gca,'ydir','reverse'); colorbar; colormap jet; caxis([-0.05 +0.15])
figure(2); clf; pcolor(w,1:096,squeeze(jac1(1,:,1:096))');   title('WVz SARTA');  shading flat; set(gca,'ydir','reverse'); colorbar; colormap jet; caxis([-0.8 +0.1])
figure(3); clf; pcolor(w,1:096,squeeze(jac3(1,:,1:096))');   title('O3z SARTA');  shading flat; set(gca,'ydir','reverse'); colorbar; colormap jet; caxis([-0.8 +0.1])
figure(4); clf; pcolor(w,1:096,squeeze(jacWGT(1,:,1:096))'); title('WGTz SARTA'); shading flat; set(gca,'ydir','reverse'); colorbar; colormap jet; caxis([0 +0.15])

kc = load('individual_prof_convolved_kcarta_airs_32_jac_clr.mat');
kc = load('individual_prof_convolved_kcarta_airs_32_jac.mat');
figure(5); ind = (1:96) + (3-1)*96; clf; pcolor(kc.fKc(1:2378),1:096,kc.rKc(1:2378,ind)'); title('Tz   KCARTA'); shading flat; colorbar; colormap jet; caxis([-0.05 +0.15])
figure(6); ind = (1:96) + (1-1)*96; clf; pcolor(kc.fKc(1:2378),1:096,kc.rKc(1:2378,ind)'); title('WVz  KCARTA'); shading flat; colorbar; colormap jet; caxis([-0.8 +0.1])
figure(7); ind = (1:96) + (2-1)*96; clf; pcolor(kc.fKc(1:2378),1:096,kc.rKc(1:2378,ind)'); title('O3z  KCARTA'); shading flat; colorbar; colormap jet; caxis([-0.8 +0.1])
figure(8); ind = (1:96) + (4-1)*96; clf; pcolor(kc.fKc(1:2378),1:096,kc.rKc(1:2378,ind)'); title('WGTz KCARTA'); shading flat; colorbar; colormap jet; caxis([0 0.15])
