addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD

%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_rtp.m
use_this_rtp = '/home/sergio/kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES_AIRS_IASI_CRIS/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = 'RTP/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = '/asl/val/airs_2009/ARMTWP/Phase1/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
use_this_rtp = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';

[h,ha,p,pa] = rtpread(use_this_rtp);
tobs = rad2bt(h.vchan,p.robs1);

for ii = 1 : 180
  fname = ['H2016/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radH16(:,ii) = x.rKc;

  fname = ['H2012/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radH12(:,ii) = x.rKc;

  fname = ['H2008/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  radH08(:,ii) = x.rKc;
end

fc = x.fKc;
t16 = rad2bt(fc,radH16);
t12 = rad2bt(fc,radH12);
t08 = rad2bt(fc,radH08);

g = dogoodchan;

figure(1)
plot(h.vchan(g),nanmean(tobs(g,:)'-t08(g,:)'),'b',h.vchan(g),nanmean(tobs(g,:)'-t12(g,:)'),'g',...
     h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'r','linewidth',2)
  plotaxis2;
  hl = legend('H08','H12','H16','location','best','fontsize',10);
  xlim([640 1640])
ylabel('Bias (K)'); xlabel('Wavenumber cm-1');

figure(2)
plot(h.vchan(g),nanstd(tobs(g,:)'-t08(g,:)'),'b',h.vchan(g),nanstd(tobs(g,:)'-t12(g,:)'),'g',...
     h.vchan(g),nanstd(tobs(g,:)'-t16(g,:)'),'r','linewidth',2)
  plotaxis2;
  hl = legend('H08','H12','H16','location','best','fontsize',10);
  xlim([640 1640]); ylim([0 1.5]);
ylabel('Std (K)'); xlabel('Wavenumber cm-1'); 

figure(1); xlim([1240 1640])
figure(2); xlim([1240 1640])
addpath /home/sergio/MATLABCODE/SHOWSTATS
figure(3); oo = find(h.vchan >= 1240 & h.vchan <= 1640); oo = intersect(h.ichan(oo),g);
  plot(tobs(oo,:),tobs(oo,:)-t08(oo,:),'.',tobs(oo,:),tobs(oo,:)-t12(oo,:),'g.',tobs(oo,:),tobs(oo,:)-t16(oo,:),'r.')
  ylabel('BTobs-BTcal'); xlabel('BTobs'); title('1240-1640 cm-1')

  [n,nx,ny,m08,s08] = myhist2d(tobs(oo,:),tobs(oo,:)-t08(oo,:),200:5:300,-5:0.25:+5);
  [n,nx,ny,m12,s12] = myhist2d(tobs(oo,:),tobs(oo,:)-t12(oo,:),200:5:300,-5:0.25:+5);
  [n,nx,ny,m16,s16] = myhist2d(tobs(oo,:),tobs(oo,:)-t16(oo,:),200:5:300,-5:0.25:+5);
  plot(200:5:300,m08,'b',200:5:300,m12,'g',200:5:300,m16,'r','linewidth',2)
  ylabel('BTobs-BTcal'); xlabel('BTobs'); title('1240-1640 cm-1'); plotaxis2;
  hl = legend('H08','H12','H16','location','best','fontsize',10);
