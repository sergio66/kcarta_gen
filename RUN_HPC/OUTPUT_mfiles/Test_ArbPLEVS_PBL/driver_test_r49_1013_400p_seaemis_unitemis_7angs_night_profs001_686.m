%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AIRS SARTA
f1A = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
f2A = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_airs';
f1A = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
f2A = '../sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_airs';
if ~exist(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/' f2A '_H2020.rtp'])
  disp('making AIRS files')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest

  toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';                                                           %% HITRAN 2016 spectroscopy with jacobians 
  toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';                        %% HITRAN 2016 spectroscopy with jacobians 
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod                                      fin=' f1A ' fout=' f2A '_H2016.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19                                           fin=' f1A ' fout=' f2A '_H2016.rtp'];
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod fin=' f1A ' fout=' f2A '_H2016.rtp'];
  %eval(sartaer)

  %% H2020
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=' f1A ' fout=' f2A '_H2020.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev                                   fin=' f1A ' fout=' f2A '_H2020.rtp'];
  eval(sartaer)
  sartaer
  eval(['ls -lt ' f2A '_H2020.rtp'])

  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% AIRS DAY SARTA
f1A  = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
f1AD = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_day.rtp';
f2AD = '../sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_airs_day';
if ~exist(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/' f2AD '_H2020.rtp'])
  disp('making AIRS DAYfiles')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest

  [hjunk,hajunk,pjunk,pajunk] = rtpread(f1A);
  pjunk.solzen = 0 * pjunk.solzen;
  rtpwrite(f1AD,hjunk,hajunk,pjunk,pajunk);

  toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';                                                           %% HITRAN 2016 spectroscopy with jacobians 
  toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';                        %% HITRAN 2016 spectroscopy with jacobians 
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod                                      fin=' f1AD ' fout=' f2AD '_H2016.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19                                           fin=' f1AD ' fout=' f2AD '_H2016.rtp'];
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod fin=' f1AD ' fout=' f2AD '_H2016.rtp'];
  %eval(sartaer)

  %% H2020
  sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=' f1AD ' fout=' f2AD '_H2020.rtp'];
  sartaer = ['!/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_p2022jul22_dev                                   fin=' f1AD ' fout=' f2AD '_H2020.rtp'];
  eval(sartaer)
  sartaer
  eval(['ls -lt ' f2AD '_H2020.rtp'])

  [hjunk,hajunk,pjunk,pajunk] = rtpread([f2A '_H2020.rtp']);
  [hjunk,hajunk,pjunk2,pajunk] = rtpread([f2AD '_H2020.rtp']);
  compare_two_structures(pjunk,pjunk2)
  plot(1:686,pjunk.solzen,'b',1:686,pjunk2.solzen,'r.-')
  tjunk  = rad2bt(hjunk.vchan,pjunk.rcalc);  
  tjunk2 = rad2bt(hjunk.vchan,pjunk2.rcalc);  
  plot(hjunk.vchan,nanmean(tjunk2'-tjunk'))
  title('SARTA : Solzen=0 - Solzen=150')

  error(';ajf;ajf')
  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CRIS SARTA
f1C = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris.rtp'
f2C = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris';
f1C = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris.rtp';
f2C = '../sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_cris';
if ~exist(['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/' f2C '_H2020.rtp'])
  disp('making CrIS files')
  cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest

  [hjunk0,hajunk0,pjunk0,pajunk0] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/r49_1013_400p_unitemis_7angs_night_crisHiRes.rtp');
  [hjunk,hajunk,pjunk,pajunk] = rtpread(f1A);
  hjunk = hjunk0;
  rtpwrite(f1C,hjunk,hajunk,pjunk,pajunk);

  %% H2020 and H2016
  % -rwxrwxr-x 1 sergio pi_strow 1502832 Jan  7 16:55 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new
  % -rwxrwxr-x 1 sergio pi_strow 1501144 Jan  6 20:47 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new
  % -rwxrwxr-x 1 sergio pi_strow 1851216 Apr 15  2024 /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new_apr15_2024

  %sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_dec17_iceGHMbaum_wdrop_ddust_sc_hg3_new_apr15_2024 '];
  %sartaer = [sartaer ' fin=' f1C ' fout=' f2C '_H2016.rtp'];    
  %eval(sartaer)

  %sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_crisg4_hires_jan25_H2020_iceGHMbaum_wdrop_ddust_sc_hg3_new'];
  %sartaer = [sartaer ' fin=' f1C ' fout=' f2C '_H2016.rtp'];
  %eval(sartaer)

  sartaer = ['!/home/chepplew/gitLib/sarta/bin/cris_hrg4_p2022jul22_dev'];
  sartaer = [sartaer ' fin=' f1C ' fout=' f2C '_H2020.rtp'];
  eval(sartaer)
  cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('reading in 2 x 343 KCARTAfiles')
dir0 = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest/';

for ii = 1 : 2*343
  if mod(ii,100) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/' dir0 '/'];
  %fname = [fname 'individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  fname = [fname 'individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  a = load(fname);
    fairs = a.fKc;
    orig_airs_kc(:,ii) = a.rKc;
    fcris0 = a.hi_fcris;
    orig_cris_kc(:,ii) = a.hi_rcris_all;

  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/'];
  fname = [fname 'individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    a = load(fname);
    new_airs_kc(:,ii) = a.rKc;

  fname = ['/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/'];
  fname = [fname 'individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
    a = load(fname);
    fcris1 = a.hi_fcris;
    new_cris_kc(:,ii) = a.hi_rcris_all;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1A = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest/r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
f2A = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_airs';
f1A = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night.rtp';
f2A = '../sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_airs';
cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest
  [h,ha,p0,pa] = rtpread(f1A);
  [hs,~,ps,~] = rtpread([f2A '_H2020.rtp']);
cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
[sum(p0.stemp-ps.stemp)]

airs_told = rad2bt(fairs,orig_airs_kc);
airs_tnew = rad2bt(fairs,new_airs_kc);

[Y,I] = sort(fairs);
plot(fairs(I),nanmean(airs_told(I,:)'-airs_tnew(I,:)'),fairs(I),nanstd(airs_told(I,:)'-airs_tnew(I,:)'))
title('old KC test - new KC test')

ind = 1 : 2378;
ts = rad2bt(fairs(ind),ps.rcalc);
plot(fairs(ind),nanmean(airs_told(ind,:)'-airs_tnew(ind,:)'),'b',fairs(ind),nanstd(airs_told(ind,:)'-airs_tnew(ind,:)'),'c',...
     fairs(ind),nanmean(airs_told(ind,:)'-ts(ind,:)'),'r',fairs(ind),nanstd(airs_told(ind,:)'-ts(ind,:)'),'m')
title('AIRS Sea+Unit emissi');
legend('ChrisSaved-kCARTAnew bias','ChrisSaved-kCARTAnew stddev','ChrisSaved-SARTAnew bias','ChrisSaved-SARTAnew stddev','location','best','fontsize',10);
disp('ret to continue'); pause

ix = 1; ix = (1:343) + (ix-1)*343;
plot(fairs(ind),nanmean(airs_told(ind,ix)'-airs_tnew(ind,ix)'),'b',fairs(ind),nanstd(airs_told(ind,ix)'-airs_tnew(ind,ix)'),'c',...
     fairs(ind),nanmean(airs_told(ind,ix)'-ts(ind,ix)'),'r',fairs(ind),nanstd(airs_told(ind,ix)'-ts(ind,ix)'),'m')
title('AIRS Sea Emiss');
legend('ChrisSaved-kCARTAnew bias','ChrisSaved-kCARTAnew stddev','ChrisSaved-SARTAnew bias','ChrisSaved-SARTAnew stddev','location','best','fontsize',10);
disp('ret to continue'); pause

ix = 2; ix = (1:343) + (ix-1)*343;
plot(fairs(ind),nanmean(airs_told(ind,ix)'-airs_tnew(ind,ix)'),'b',fairs(ind),nanstd(airs_told(ind,ix)'-airs_tnew(ind,ix)'),'c',...
     fairs(ind),nanmean(airs_told(ind,ix)'-ts(ind,ix)'),'r',fairs(ind),nanstd(airs_told(ind,ix)'-ts(ind,ix)'),'m')
title('AIRS Unit Emiss');
legend('ChrisSaved-kCARTAnew bias','ChrisSaved-kCARTAnew stddev','ChrisSaved-SARTAnew bias','ChrisSaved-SARTAnew stddev','location','best','fontsize',10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to continue to Cris'); pause

f1C = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris.rtp'
f2C = 'JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris';
f1C = 'r49_1013_98lev_400p_unitemis_seaemis_7angs_night_cris.rtp';
f2C = '../sarta_prod_2022_REGR49_BestTest_SergioJan2025/r49_cris';

cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/ChrisH/sarta_prod_2022_REGR49_BestTest
  [h,ha,p0,pa] = rtpread(f1C);
  [hs,~,ps,~] = rtpread([f2C '_H2020.rtp']);
cd /umbc/xfs2/strow/asl/s1/sergio/KCARTA_RUNTARA/JUNK/Test_ArbPLEVS_PBL
[sum(p0.stemp-ps.stemp)]

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

mmw = mmwater_rtp(h,p0);
o3du = dobson_rtp(h,p0);
