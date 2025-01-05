clear all

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /home/sergio/MATLABCODE/SHOWSTATS

nedt = instr_chans('airs',2);

file{1} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/armtwp_nearby_semiclear_jan04_resetT_nte_co2.rtp';
file{2} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/tairs_382co2_105o3_104n2o_070co_101ch4_050hno3_rad.rtp';
file{3} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/voem_start.rtp';
file{4} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp3_start.rtp';
file{5} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp2_start.rtp';
file{6} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/twp1_start.rtp';
file{7} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp3_start.rtp';
file{8} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp2_start.rtp';
file{9} = '/home/chepplew/data/sarta_tuning/sorted/from_scott/sgp1_start.rtp';
file{10}= '/home/chepplew/data/sarta_tuning/sorted/from_scott/minn_start.rtp';
file{11}= '/home/chepplew/data/sarta_tuning/frm_scott/tiasi_382co2_105o3_104n2o_101ch4_050hno3_rad.rtp_1';

topts.sarta08 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';         %% HITRAN 2008 spectroscopy OLD DEFAULT SCOTT, better results?
topts.sarta08 = '/asl/packages/sartaV108/BinV201/sarta_airs_PGEv6_postNov2003';         %% HITRAN 2008 spectroscopy OLD DEFAULT SCOTT, better results?
topts.sarta08 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_iceaggr_waterdrop_desertdust_slabcloud_hg3'; %% no tuning?
topts.sarta08 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140'; %% no tuning?

topts.sarta16 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v2';                                        %% HITRAN 2016 spectroscopy NEW DEFAULT CHRIS, odder results? FAST, wrong data location
topts.sarta16 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';                                        %% HITRAN 2016 spectroscopy NEW DEFAULT CHRIS, odder results? FAST, correct data location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iFile = input('Enter file AIRS (1:10) IASI 11 : ');
[h,ha,p,pa] = rtpread(file{iFile});
fprintf(1,'%s : length(p.stemp) = %5i .. running sarta .. \n',file{iFile},length(p.stemp));

if iFile >= 3 & iFile <= 10
  sartaer = ['!time ' topts.sarta08 ' fin=' file{iFile} ' fout=junk.rp.rtp >& ugh']; eval(sartaer); [h08,ha08,p08,pa08] = rtpread('junk.rp.rtp'); st08 = rad2bt(h08.vchan,p08.rcalc);
  sartaer = ['!time ' topts.sarta16 ' fin=' file{iFile} ' fout=junk.rp.rtp >& ugh']; eval(sartaer); [h16,ha16,p16,pa16] = rtpread('junk.rp.rtp'); st16 = rad2bt(h16.vchan,p16.rcalc);
elseif iFile >= 11
  disp('no sarta runs because this is IASI')
else
  disp('no sarta runs because the rtp files have issues')
end

if iFile <= 10
  rad = zeros(length(p.stemp),2834,4);
else
  rad = zeros(length(p.stemp),4231,4);
end

iaaFound = zeros(length(p.stemp),4);
for ii = 1 : length(p.stemp)
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  if iFile <= 10
    strI = 'airs';
    a.rKc = zeros(2834,1);
    b.rKc = zeros(2834,1);
    c.rKc = zeros(2834,1);
    d.rKc = zeros(2834,1);
  else
    strI = 'iasi';
    a.rKcIasi = zeros(8461,1);
    b.rKcIasi = zeros(8461,1);
    c.rKcIasi = zeros(8461,1);
    d.rKcIasi = zeros(8461,1);
    a.rKcIasi = zeros(4231,1);
    b.rKcIasi = zeros(4231,1);
    c.rKcIasi = zeros(4231,1);
    d.rKcIasi = zeros(4231,1);
  end
  fx = ['FILE' num2str(iFile) '/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2016.mat']; 
  if exist(fx)
    a = load(fx);
    iaaFound(ii,1) = 1;
  end
  fx = ['FILE' num2str(iFile) '/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2015.mat']; 
  if exist(fx)
    b = load(fx);
    iaaFound(ii,2) = 1;
  end
  fx = ['FILE' num2str(iFile) '/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2012.mat']; 
  if exist(fx)
    c = load(fx);
    iaaFound(ii,3) = 1;
  end
  fx = ['FILE' num2str(iFile) '/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2008.mat']; 
  if exist(fx)
    d = load(fx);
    iaaFound(ii,4) = 1;
  end
  
  if iFile == 5 | iFile == 8
    fx = ['FILE' num2str(iFile) '_HighRes/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2017.mat']; 
    fx = ['FILE' num2str(iFile) '_HighRes_linear_in_tau/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2017.mat']; 
    if exist(fx)
      e = load(fx);
      iaaFound(ii,5) = 1;
      rad_high(ii,:) = e.rKc;
    end

    fx = ['FILE' num2str(iFile) '_usualres_3gases/individual_prof_convolved_kcarta_' strI '_' num2str(ii) '_HITorGEI2016.3.mat']; 
    if exist(fx)
      f = load(fx);
      iaaFound(ii,6) = 1;
      rad_low(ii,:) = f.rKc;
    end

  end

  if iFile <= 10
    rad(ii,:,:) = [a.rKc b.rKc c.rKc d.rKc];
  else
    rad(ii,:,:) = [a.rKcIasi(1:4231) b.rKcIasi(1:4231) c.rKcIasi(1:4231) d.rKcIasi(1:4231)];
  end
end
fprintf(1,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iFile <= 10
  g = dogoodchan;
  nchan = 2378;
else
  g = 1 : 4231;
  nchan = 4231;
end

boo = sum(iaaFound(:,1:4),2); allfour = find(boo == 4); fprintf(1,'found %5i of %5i had all 4 kcarta runs \n',length(allfour),length(p.stemp))
if length(allfour) < length(p.stemp)
  disp('bad profiles')
  bad = find(boo < 4)
end

tobs = rad2bt(h.vchan,p.robs1(:,allfour)); tobs = real(tobs);
t16  = rad2bt(h.vchan,squeeze(rad(allfour,1:nchan,1))');
t15  = rad2bt(h.vchan,squeeze(rad(allfour,1:nchan,2))');
t12  = rad2bt(h.vchan,squeeze(rad(allfour,1:nchan,3))');
t08  = rad2bt(h.vchan,squeeze(rad(allfour,1:nchan,4))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iFile == 5 | iFile == 8
  t16_low   = rad2bt(h.vchan,rad_low(:,1:nchan)');
  t16_high  = rad2bt(h.vchan,rad_high(:,1:nchan)');
  figure(1); clf
  plot(h.vchan(g),nanmean(tobs(g,:)'-t16_high(g,:)'),'r.-',h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'b','linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('High res H16','Usual res H16','location','best','fontsize',10); 
    title('mean(obs-kCARTA) : High vs USUAL res H16')
    axis([640 840 -1 +1])
    axis([640 740 -1 +1])

  figure(2); clf
  plot(h.vchan(g),nanmean(tobs(g,:)'-t16_high(g,:)') - nanmean(tobs(g,:)'-t16(g,:)'),'b','linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('High res - Usual res H16','location','best','fontsize',10); 
   title('mean (High vs USUAL ALL gas res H16)')
    axis([640 840 -1 +1])
    axis([640 740 -1 +1])
  plot(h.vchan(g),nanmean(t16_high(g,:)'-t16(g,:)'),'b',...
       h.vchan(g),nanstd(t16_high(g,:)'-t16(g,:)'),'c',...
         'linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('mean','std','location','best','fontsize',10); title('High vs USUAL ALL gas res H16')
    axis([640 840 -1 +1])
    axis([640 740 -1 +1])

  figure(3); clf
  ix = find(h.vchan >= 640 & h.vchan <= 740); ix = intersect(ix,g);
  Y = tobs(ix,:); Y1 = tobs(ix,:)-t16(ix,:); Y2 = tobs(ix,:)-t16_high(ix,:);
  Y = Y(:);       Y1 = Y1(:);                Y2 = Y2(:); 
  [~,I] = sort(Y);
  plot(Y(I),Y1(I),'b.-',Y(I),Y2(I),'r.-','linewidth',2);
    plotaxis2;
    hl = legend('H16','H16 High','location','best','fontsize',10); grid on
    xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias']);
    xlim([200 300]);  ylim([-2 +2])
  
  [n,nx,ny,m1,s1] = myhist2d(Y,Y1,180:2:280,-5:0.1:+5);
  [n,nx,ny,m2,s2] = myhist2d(Y,Y2,180:2:280,-5:0.1:+5);
  plot(180:2:280,m1,'b.-',180:2:280,m2,'r.-','linewidth',2);
  errorbar((180:2:280)-0.2,m1,s1,'color','b','linewidth',2); hold on
  errorbar((180:2:280)-0.1,m2,s2,'color','r','linewidth',2); hold on
    plotaxis2;
    hl = legend('H16','H16 High','location','best','fontsize',10); grid on
    xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias,std']);
    xlim([180 280]);  ylim([-1 +1])

  figure(4); clf
  plot(h.vchan(g),tobs(g,1)-t16_high(g,1),'rx-',h.vchan(g),tobs(g,1)-t16_low(g,1),'b','linewidth',2);  
  plotaxis2;
  title('Profile 1 (r) bias high res (r) bias usual res');
  xlim([640 840])
  ylim([-1 +1])
  
  disp('finished high vs usual res ... ret to continue');pause;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%
if iFile == 5 | iFile == 8
  t16_low   = rad2bt(h.vchan,rad_low(:,1:nchan)');
  t16_high  = rad2bt(h.vchan,rad_high(:,1:nchan)');

  figure(1); clf
  plot(h.vchan(g),nanmean(tobs(g,:)'-t16_high(g,:)'),'r.-',h.vchan(g),nanmean(tobs(g,:)'-t16_low(g,:)'),'b','linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('High res H16','Usual res H16 3gas','location','best','fontsize',10); 
    title('mean(Obs-kCARTA) : High vs USUAL 3gas res H16')
    axis([640 840 -1 +1])
    axis([640 740 -1 +1])

  figure(2); clf
  plot(h.vchan(g),nanmean(tobs(g,:)'-t16_high(g,:)') - nanmean(tobs(g,:)'-t16_low(g,:)'),'b',...
        h.vchan(g),nanstd((tobs(g,:)'-t16_high(g,:)')) - ((tobs(g,:)'-t16_low(g,:)')),'c',...
       'linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('High res - Usual res H16 3gas','location','best','fontsize',10); title('High vs USUAL 3 gas res H16')
  plot(h.vchan(g),nanmean(t16_high(g,:)' - t16_low(g,:)'),'b',...
       h.vchan(g),nanstd(t16_high(g,:)' - t16_low(g,:)'),'c',...
       'linewidth',2);
    plotaxis2; xlim([640 840]); hl = legend('mean','std','location','best','fontsize',10); 
    title('High vs USUAL 3 gas res H16')
    axis([640 840 -1 +1])
    axis([640 740 -1 +1])

  figure(3); clf
  ix = find(h.vchan >= 640 & h.vchan <= 740); ix = intersect(ix,g);
  Y = tobs(ix,:); Y1 = tobs(ix,:)-t16_low(ix,:); Y2 = tobs(ix,:)-t16_high(ix,:);
  Y = Y(:);       Y1 = Y1(:);                Y2 = Y2(:); 
  [~,I] = sort(Y);
  plot(Y(I),Y1(I),'b.-',Y(I),Y2(I),'r.-','linewidth',2);
    plotaxis2;
    hl = legend('H16','H16 High','location','best','fontsize',10); grid on
    xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias']);
    xlim([200 300]);  ylim([-2 +2])
  
  [n,nx,ny,m1,s1] = myhist2d(Y,Y1,180:2:280,-5:0.1:+5);
  [n,nx,ny,m2,s2] = myhist2d(Y,Y2,180:2:280,-5:0.1:+5);
  plot(180:2:280,m1,'b.-',180:2:280,m2,'r.-','linewidth',2);
  errorbar((180:2:280)-0.2,m1,s1,'color','b','linewidth',2); hold on
  errorbar((180:2:280)-0.1,m2,s2,'color','r','linewidth',2); hold on
    plotaxis2;
    hl = legend('H16 3gas','H16 High','location','best','fontsize',10); grid on
    xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias,std']);
    xlim([180 280]);  ylim([-1 +1])

  figure(4); clf
  plot(h.vchan(g),t16_high(g,1) - t16_low(g,1),'bx-',h.vchan(g),nanmean(t16_high(g,:)'-t16_low(g,:)'),'r','linewidth',2);  
  title('high res - usual res (b) Prof 1 (r) all ');
  plotaxis2;
  xlim([640 840])

  disp('finished high vs usual 3gas res ... ret to continue');pause;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf
plot(h.vchan(g),nanmean(tobs(g,:)'-t08(g,:)'),'k.-',h.vchan(g),nanmean(tobs(g,:)'-t12(g,:)'),'b',h.vchan(g),nanmean(tobs(g,:)'-t15(g,:)'),'g',h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'r',...
     h.vchan(g),nanstd(tobs(g,:)'-t08(g,:)'),'k--',h.vchan(g),nanstd(tobs(g,:)'-t12(g,:)'),'b--',h.vchan(g),nanstd(tobs(g,:)'-t15(g,:)'),'g--',h.vchan(g),nanstd(tobs(g,:)'-t16(g,:)'),'r--',...
     h.vchan(g),+nedt(g),'k',h.vchan(g),-nedt(g),'k','linewidth',2)
plotaxis2;
xlim([640 1640]);  ylim([-2 +2])
yyaxis right
  plot(h.vchan(g),nanmean(tobs(g,:)'),'kx-','linewidth',0.25)
  ax = gca;
  ax.YAxis(1).Color = 'k';
  ax.YAxis(2).Color = 'k';
hl = legend('H08','H12','G15','H16','<BT>','location','best','fontsize',10); grid on
xlabel('Wavenumber \nu cm^{-1}'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' (solid) bias (dash) std']);

figure(1); clf
plot(h.vchan(g),nanmean(tobs(g,:)'-t08(g,:)'),'k.-',h.vchan(g),nanmean(tobs(g,:)'-t12(g,:)'),'b',h.vchan(g),nanmean(tobs(g,:)'-t15(g,:)'),'g',h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'r',...
     h.vchan(g),nanstd(tobs(g,:)'-t08(g,:)'),'k--',h.vchan(g),nanstd(tobs(g,:)'-t12(g,:)'),'b--',h.vchan(g),nanstd(tobs(g,:)'-t15(g,:)'),'g--',h.vchan(g),nanstd(tobs(g,:)'-t16(g,:)'),'r--',...
     h.vchan(g),+nedt(g),'k',h.vchan(g),-nedt(g),'k','linewidth',2)
plotaxis2;
xlim([640 1640]);  ylim([-2 +2])
hl = legend('H08','H12','G15','H16','location','best','fontsize',10); grid on
xlabel('Wavenumber \nu cm^{-1}'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' (solid) bias (dash) std']);

figure(2)
ix = find(h.vchan >= 1320 & h.vchan <= 1700); ix = intersect(ix,g);
Y = tobs(ix,:); Y1 = tobs(ix,:)-t08(ix,:); Y2 = tobs(ix,:)-t12(ix,:); Y3 = tobs(ix,:)-t15(ix,:); Y4 = tobs(ix,:)-t16(ix,:);
Y = Y(:); Y1 = Y1(:); Y2 = Y2(:); Y3 = Y3(:); Y4 = Y4(:); 
[~,I] = sort(Y);
plot(Y(I),Y1(I),'k.-',Y(I),Y2(I),'b.-',Y(I),Y3(I),'g.-',Y(I),Y4(I),'r.-','linewidth',2);
  plotaxis2;
  hl = legend('H08','H12','G15','H16','location','best','fontsize',10); grid on
  xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias']);
  xlim([200 300]);  ylim([-2 +2])

[n,nx,ny,m1,s1] = myhist2d(Y,Y1,200:5:300,-5:0.2:+5);
[n,nx,ny,m2,s2] = myhist2d(Y,Y2,200:5:300,-5:0.2:+5);
[n,nx,ny,m3,s3] = myhist2d(Y,Y3,200:5:300,-5:0.2:+5);
[n,nx,ny,m4,s4] = myhist2d(Y,Y4,200:5:300,-5:0.2:+5);
plot(200:5:300,m1,'k.-',200:5:300,m2,'b.-',200:5:300,m3,'g.-',200:5:300,m4,'r.-','linewidth',2);
errorbar((200:5:300)-0.2,m1,s1,'color','k','linewidth',2); hold on
errorbar((200:5:300)-0.1,m2,s2,'color','b','linewidth',2); hold on
errorbar((200:5:300)+0.1,m3,s3,'color','g','linewidth',2); hold on
errorbar((200:5:300)+0.2,m4,s4,'color','r','linewidth',2); hold off
  plotaxis2;
  hl = legend('H08','H12','G15','H16','location','best','fontsize',10); grid on
  xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias,std']);
  xlim([200 300]);  ylim([-2 +2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp('ret to continue to SARTA/kCARTA comparisons'); pause
figure(1); clf
plot(h.vchan(g),nanmean(tobs(g,:)'-t08(g,:)'),'k.-',h.vchan(g),nanmean(tobs(g,:)'-st08(g,allfour)'),'b',h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'g',h.vchan(g),nanmean(tobs(g,:)'-st16(g,allfour)'),'r',...
     h.vchan(g),nanstd(tobs(g,:)'-t08(g,:)'),'k--',h.vchan(g),nanstd(tobs(g,:)'-st08(g,allfour)'),'b--',h.vchan(g),nanstd(tobs(g,:)'-t16(g,:)'),'g--',h.vchan(g),nanstd(tobs(g,:)'-st16(g,allfour)'),'r--',...
     h.vchan(g),+nedt(g),'k',h.vchan(g),-nedt(g),'k','linewidth',2)
plotaxis2;
xlim([640 1640]);  ylim([-2 +2])
hl = legend('kCARTA H08','Sarta H08','kCARTA H16','Sarta H16''location','best','fontsize',10); grid on
xlabel('Wavenumber \nu cm^{-1}'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' (solid) bias (dash) std']);

figure(2); clf
ix = find(h.vchan >= 1320 & h.vchan <= 1700); ix = intersect(ix,g);
Y = tobs(ix,:); Y1 = tobs(ix,:)-t08(ix,:); Y2 = tobs(ix,:)-st08(ix,allfour); Y3 = tobs(ix,:)-t16(ix,:); Y4 = tobs(ix,:)-st16(ix,allfour);
Y = Y(:); Y1 = Y1(:); Y2 = Y2(:); Y3 = Y3(:); Y4 = Y4(:); 
[~,I] = sort(Y);
plot(Y(I),Y1(I),'k.-',Y(I),Y2(I),'b.-',Y(I),Y3(I),'g.-',Y(I),Y4(I),'r.-','linewidth',2);
  plotaxis2;
  hl = legend('H08','Sarta H08','H16','Sarta H16','location','best','fontsize',10); grid on
  xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias']);
  xlim([200 300]);  ylim([-2 +2])

[n,nx,ny,m1,s1] = myhist2d(Y,Y1,200:5:300,-5:0.2:+5);
[n,nx,ny,m2,s2] = myhist2d(Y,Y2,200:5:300,-5:0.2:+5);
[n,nx,ny,m3,s3] = myhist2d(Y,Y3,200:5:300,-5:0.2:+5);
[n,nx,ny,m4,s4] = myhist2d(Y,Y4,200:5:300,-5:0.2:+5);
plot(200:5:300,m1,'k.-',200:5:300,m2,'b.-',200:5:300,m3,'g.-',200:5:300,m4,'r.-','linewidth',2);
errorbar((200:5:300)-0.2,m1,s1,'color','k','linewidth',2); hold on
errorbar((200:5:300)-0.1,m2,s2,'color','b','linewidth',2); hold on
errorbar((200:5:300)+0.1,m3,s3,'color','g','linewidth',2); hold on
errorbar((200:5:300)+0.2,m4,s4,'color','r','linewidth',2); hold off
  plotaxis2;
  hl = legend('H08','Sarta H08','H16','Sarta H16','location','best','fontsize',10); grid on
  xlabel('BT obs'); ylabel('BT(K)'); title(['File ' num2str(iFile) ' bias,std']);
  xlim([200 300]);  ylim([-2 +2])

