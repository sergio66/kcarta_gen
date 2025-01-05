addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/KCARTA/MATLAB

set(0,'DefaultLegendAutoUpdate','off')

fin = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/armtwp_nearby_semiclear_jan04_resetT_nte.rtp';
fCO2 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/armtwp_nearby_semiclear_jan04_resetT_nte_addCO2_CH4.rtp';

if ~exist(fCO2)
  [h0,ha,p0,pa] = rtpread(fin);

  %% now need to add in CO2 and CH4
  addpath //home/sergio/MATLABCODE/SARTA_Tuning/Scott_ResetCode/Sergio_resets
  [h,ha,p,pa] = reset_co2ppm(fin,fCO2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hh,hha,pp,ppa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/armtwp_nearby_semiclear_jan04_resetT_nte_oldrtime.rtp');
get_attr(hha)

[h,ha,p,pa] = rtpread(fCO2);
[yy,mm,dd,hh] = tai2utcSergio(p.rtime+offset1958_to_1993);
tuning = load('/asl/data/sarta_database/Data_jan04untun/Coef/tunmlt_jan04deliv.txt');
plot(tuning(:,2),tuning(:,[3:5]),'linewidth',2); xlim([640 1640])
hl = legend('fixed','WV lines','WV con','location','best');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 180
  istr = num2str(ii);
  f08 = ['H08/individual_prof_convolved_kcarta_crisHI_crisMED_' istr '.mat'];
  f12 = ['H12/individual_prof_convolved_kcarta_crisHI_crisMED_' istr '.mat'];
  f16 = ['H16/individual_prof_convolved_kcarta_crisHI_crisMED_' istr '.mat'];

  a = load(f08);
  rKc08(:,ii) = a.rKc;

  a = load(f12);
  rKc12(:,ii) = a.rKc;

  a = load(f16);
  rKc16(:,ii) = a.rKc;

end
tobs = rad2bt(h.vchan,p.robs1);
tcal = rad2bt(h.vchan,p.rcalc);
t08 = rad2bt(a.fKc,rKc08); t08 = t08(1:2378,:);
t12 = rad2bt(a.fKc,rKc12); t12 = t12(1:2378,:);
t16 = rad2bt(a.fKc,rKc16); t16 = t16(1:2378,:);

addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
g = dogoodchan;
plot(h.vchan(g),nanmean(tobs(g,:)'-tcal(g,:)'),'k',h.vchan(g),nanmean(tobs(g,:)'-t08(g,:)'),'b',...
     h.vchan(g),nanmean(tobs(g,:)'-t12(g,:)'),'g',h.vchan(g),nanmean(tobs(g,:)'-t16(g,:)'),'r')
xlim([640 1640]); ylim([-1 +1])
plotaxis2;
hl = legend('sarta1.06','H08','H12','H16','location','best');
hold on;
plot(h.vchan(g),nanmean(tobs(g,:)'-tcal(g,:)'),'k','linewidth',2);
hold off

[ppmvLAY,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(h,p,1:length(p.stemp),2);
