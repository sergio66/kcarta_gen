use_this_rtp = 'RTP/r49_1013_400p_8x2x2_2834_pbl.rtp';     %% 49 profiles x 2 surfaces x 8 view  angles * 2 sol angles   TEST1
use_this_rtp = 'RTP/r49_1013_400p_8x2x2_2834_airslay.rtp'; %% 49 profiles x 2 surfaces x 8 view  angles * 2 sol angles   TEST2 
[h,ha,p,pa] = rtpread(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' use_this_rtp]);

if ~exist('TEST2_Jan29_2025/sartaSergioH2020.rtp.rtp')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/
  sarta = ['/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020'];
  sartaer = ['!time ' sarta ' fin=' use_this_rtp  ' fout=JUNK/Test_ArbPLEVS_PBL/TEST2_Jan29_2025/sartaSergioH2020.rtp.rtp'];
  eval(sartaer)
  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 1568
  fin = ['TEST1_Jan29_2025/individual_prof_convolved_kcarta_airs_iasi_crisHI_crisMED_' num2str(ii) '.mat'];
  a = load(fin);
  pbl100(ii,:) = a.rKc;
end

for ii = 1 : 1568
  fin = ['TEST2_Jan29_2025/individual_prof_convolved_kcarta_airs_iasi_crisHI_crisMED_' num2str(ii) '.mat'];
  a = load(fin);
  airs100(ii,:) = a.rKc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hSARTA100,~,pSARTA100,~] = rtpread('TEST2_Jan29_2025/sartaSergioH2020.rtp.rtp');
tSARTA100 = rad2bt(hSARTA100.vchan,pSARTA100.rcalc);

fKc = a.fKc;
tpbl100  = rad2bt(fKc,pbl100');
tairs100 = rad2bt(fKc,airs100');

oo = find(abs(p.satzen) < 50);
[Y,I] = sort(fKc);
plot(fKc-hSARTA100.vchan)

plot(fKc(I),nanmean(tairs100(I,:)'-tSARTA100(I,:)'),fKc(I),nanstd(tairs100(I,:)'-tSARTA100(I,:)'))
xlabel('Wavenumber [cm-1]'); ylabel('AIRS100 : KCARTA-SARTA \newline \delta BT [K]');
legend('bias','stddev','location','best');
set(gca,'fontsize',12)
%% sergioprint('/home/sergio/PAPERS/REPORTS/ForStrow/RTA/rta_report_2025/diff_airs100_pbl100');

plot(fKc(I),nanmean(tairs100(I,oo)'-tSARTA100(I,oo)'),fKc(I),nanstd(tairs100(I,oo)'-tSARTA100(I,oo)'))
xlabel('Wavenumber [cm-1]'); ylabel('AIRS100 satzen < 50: KCARTA-SARTA \newline \delta BT [K]');
legend('bias','stddev','location','best');
set(gca,'fontsize',12)

plot(1:1568,tairs100(1291,:),1:1568,tpbl100(1291,:),1:1568,tSARTA100(1291,:))
plot(1:1568,tairs100(1291,:)-tpbl100(1291,:),1:1568,tairs100(1291,:)-tSARTA100(1291,:),1:1568,tairs100(1291,:)-p.stemp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


plot(fKc(I),nanmean(tairs100(I,:)'-tpbl100(I,:)'),fKc(I),nanstd(tairs100(I,:)'-tpbl100(I,:)'))
xlabel('Wavenumber [cm-1]'); ylabel('AIRS100 - PBL100 \newline \delta BT [K]');
legend('bias','stddev','location','best');
set(gca,'fontsize',12)
%% sergioprint('/home/sergio/PAPERS/REPORTS/ForStrow/RTA/rta_report_2025/diff_airs100_pbl100');

night = find(p.solzen > 90);
day = find(p.solzen <= 90);

plot(fKc(I),nanmean(tairs100(I,day)'-tpbl100(I,day)'),fKc(I),nanstd(tairs100(I,day)'-tpbl100(I,day)'))
xlabel('Wavenumber [cm-1]'); ylabel('AIRS100 - PBL100 \newline \delta BT [K]');
legend('bias','stddev','location','best');
set(gca,'fontsize',12)
%% sergioprint('/home/sergio/PAPERS/REPORTS/ForStrow/RTA/rta_report_2025/diff_airs100_pbl100_day');

plot(fKc(I),nanmean(tairs100(I,night)'-tpbl100(I,night)'),fKc(I),nanstd(tairs100(I,night)'-tpbl100(I,night)'))
xlabel('Wavenumber [cm-1]'); ylabel('AIRS100 - PBL100 \newline \delta BT [K]');
legend('bias','stddev','location','best');
set(gca,'fontsize',12)
%% sergioprint('/home/sergio/PAPERS/REPORTS/ForStrow/RTA/rta_report_2025/diff_airs100_pbl100_night');
