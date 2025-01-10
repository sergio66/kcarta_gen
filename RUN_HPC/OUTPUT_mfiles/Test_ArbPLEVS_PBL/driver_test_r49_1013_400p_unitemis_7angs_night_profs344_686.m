%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AIRS SARTA
if ~exist('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night.rp.rtp')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/

  [hjunk,hajunk,pjunk,pajunk] = rtpread('RTP/r49_1013_400p_unitemis_7angs_night.rtp');
  pjunk.rho = 0 * pjunk.rho;
  rtpwrite('RTP/r49_1013_400p_unitemis_7angs_night.rtp',hjunk,hajunk,pjunk,pajunk);

  toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';                                                           %% HITRAN 2016 spectroscopy with jacobians 
  toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';                        %% HITRAN 2016 spectroscopy with jacobians 
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod fin=RTP/r49_1013_400p_unitemis_7angs_night.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_H16.rp.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19 fin=RTP/r49_1013_400p_unitemis_7angs_night.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_H16.rp.rtp'];
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod fin=RTP/r49_1013_400p_unitemis_7angs_night.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_H16_sergio.rp.rtp'];
  eval(sartaer)

  %% H2020
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=RTP/r49_1013_400p_unitemis_7angs_night.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night.rp.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev fin=RTP/r49_1013_400p_unitemis_7angs_night.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night.rp.rtp'];
  eval(sartaer)
  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CRIS SARTA
if ~exist('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night.rtp')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/

  [hjunk,hajunk,pjunk,pajunk] = rtpread('RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp');
  pjunk.rho = 0 * pjunk.rho;
  rtpwrite('RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp',hjunk,hajunk,pjunk,pajunk);

  %% H2020 and H2016
  % -rwxrwxr-x 1 sergio pi_strow 1502832 Jan  7 16:55 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new
  % -rwxrwxr-x 1 sergio pi_strow 1501144 Jan  6 20:47 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new
  % -rwxrwxr-x 1 sergio pi_strow 1851216 Apr 15  2024 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new_apr15_2024

  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new_apr15_2024 '];
  sartaer = [sartaer ' fin=RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_cris_H16_sergio.rp.rtp'];
  eval(sartaer)

  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new'];
  sartaer = [sartaer ' fin=RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_cris_H20_sergio.rp.rtp'];
  eval(sartaer)

  sartaer = ['!/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev'];
  sartaer = [sartaer ' fin=RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp fout=RTP/r49_1013_400p_unitemis_7angs_night_cris_H20_sergio_chrisexec.rp.rtp'];
  eval(sartaer)
  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('reading in 343 UNIT EMISS files')
for ii = 1 : 343
  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night/'];
  %fname = [fname 'individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  fname = [fname 'individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  a = load(fname);
    fairs = a.fKc;
    orig_airs_kc(:,ii) = a.rKc;
    fcris0 = a.hi_fcris;
    orig_cris_kc(:,ii) = a.hi_rcris_all;

  %fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night_SergioTest_Jan2025/'];
  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night_SergioTest_Jan2025_FixRho/'];
  fname = [fname 'individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    a = load(fname);
    new_airs_kc(:,ii) = a.rKc;

  %fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night_SergioTest_Jan2025/'];
  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night_SergioTest_Jan2025_FixRho/'];
  fname = [fname 'individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
    a = load(fname);
    fcris1 = a.hi_fcris;
    new_cris_kc(:,ii) = a.hi_rcris_all;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wahOld = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night/individual_prof_convolved_kcarta_crisHI_crisMED_105.mat');
wahNew = load('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/r49_1013_400p_unitemis_7angs_night_SergioTest_Jan2025/individual_prof_convolved_kcarta_airs_105.mat');

[h,ha,p0,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49/RTP/r49_1013_400p_unitemis_7angs_night.rtp');
[h,ha,p1,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night.rtp');

[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_H16.rp.rtp');
[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_H16_sergio.rp.rtp');
[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night.rp.rtp');

[sum(p0.stemp-p1.stemp) sum(p0.stemp-ps.stemp)]

airs_told = rad2bt(fairs,orig_airs_kc);
airs_tnew = rad2bt(fairs,new_airs_kc);

[Y,I] = sort(fairs);
plot(fairs(I),nanmean(airs_told(I,:)'-airs_tnew(I,:)'),fairs(I),nanstd(airs_told(I,:)'-airs_tnew(I,:)'))
title('old KC test - new KC test')

ind = 1 : 2378;
ts = rad2bt(fairs(ind),ps.rcalc);
plot(fairs(ind),nanmean(airs_told(ind,:)'-airs_tnew(ind,:)'),'b',fairs(ind),nanstd(airs_told(ind,:)'-airs_tnew(ind,:)'),'c',...
     fairs(ind),nanmean(airs_told(ind,:)'-ts(ind,:)'),'r',fairs(ind),nanstd(airs_told(ind,:)'-ts(ind,:)'),'m')
title('AIRS');
legend('ChrisSaved-kCARTAnew bias','ChrisSaved-kCARTAnew stddev','ChrisSaved-SARTAnew bias','ChrisSaved-SARTAnew stddev','location','best','fontsize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue to Cris'); pause

[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_cris_H16_sergio.rp.rtp');
[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_cris_H20_sergio.rp.rtp');
[hs,~,ps,~] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_cris_H20_sergio_chrisexec.rp.rtp');

[sum(p0.stemp-p1.stemp) sum(p0.stemp-ps.stemp)]

cris_told = rad2bt(fcris0,orig_cris_kc);
cris_tnew = rad2bt(fcris1,new_cris_kc);

[Y,I0,I1] = intersect(fcris0,fcris1);
plot(fcris0(I0),nanmean(cris_told(I0,:)'-cris_tnew(I1,:)'),fcris0(I0),nanstd(cris_told(I0,:)'-cris_tnew(I1,:)'))
title('old KC test - new KC test')
ylim([-1 +1]*0.2)

ind = I0;
ts = rad2bt(hs.vchan,ps.rcalc);
[Y,I3,I4] = intersect(fcris0,hs.vchan);
plot(fcris0(I0),nanmean(cris_told(I0,:)'-cris_tnew(I1,:)'),fcris0(I0),nanstd(cris_told(I0,:)'-cris_tnew(I1,:)'),...
     fcris0(I3),nanmean(cris_told(I3,:)'-ts(I4,:)'),'r',fcris0(I3),nanstd(cris_told(I3,:)'-ts(I4,:)'),'m')
ylim([-1 +1])
ylim([-0.1 +0.5])

title('CrIS HiRes');
legend('ChrisSaved-kCARTAnew bias','ChrisSaved-kCARTAnew stddev','ChrisSaved-SARTAnew bias','ChrisSaved-SARTAnew stddev','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to compare to CHris Data'); pause;

% BTW - the kCARTA calcs I used to do the CrIS HR comparison are at:
% /home/chepplew/data/sarta/prod_2022/generic//kcarta/REGR49/r49_1013_98lev_400p_unitemis_seaemis_7angs_night/
% and this rtp set:
% /home/chepplew/data/sarta/prod_2022/generic/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp
% (option to choose unit emisivity or sea). (edited) 

[hjunk,~,pjunk,~] = rtpread('/home/chepplew/data/sarta/prod_2022/generic/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp');
[hjunk,pjunk] = subset_rtp(hjunk,pjunk,[],[],(1:343) + 343);
compare_two_structures(ps,pjunk);
plot(ps.efreq,ps.emis,'b',pjunk.efreq,pjunk.emis,'r'); title('emissivity (b) sergio (r) chris')
plot(ps.efreq,ps.rho,'b',pjunk.efreq,pjunk.rho,'r'); title('reflectivity (b) sergio (r) chris')
