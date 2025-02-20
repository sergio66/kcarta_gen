%% see plot_check_sub6000clr.m

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

i2326 = find(hC.vchan >= 2326,1);
scatter(psarta_clr_rad_airs.solzen,tsarta_clr_rad_airsX(i2326,:)-tsarta_clr_rad_airs0(i2326,:),10,psarta_clr_rad_airs.rlat,'filled'); colorbar
xlabel('Solzen'); ylabel('BT 2326 : Chis-Sergio sarta'); title('Colorbar = Latitude'); colormap jet

ugh = find(tsarta_clr_rad_airsX(i2326,:)-tsarta_clr_rad_airs0(i2326,:) > 7 & psarta_clr_rad_airs.solzen < 30 & psarta_clr_rad_airs.rlat < 30)
scatter(psarta_clr_rad_airs.solzen(ugh),tsarta_clr_rad_airsX(i2326,ugh)-tsarta_clr_rad_airs0(i2326,ugh),10,psarta_clr_rad_airs.rlat(ugh),'filled'); colorbar

[hjunk,pjunk] = subset_rtp(hC,psarta_clr_rad_airs,[],hC.ichan(i2326),ugh(1));
[hjunk,pjunk] = cat_rtp(hjunk,pjunk,hjunk,pjunk);
pjunk.solzen(1) = 150;

rtpwrite('//home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/TEST_RTP/bad_nlte_sartajac_H2020.rtp',hjunk,ha,pjunk,pa)

cd /home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/TEST_RTP/
sartaer = ['!../bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=TEST_RTP/bad_nlte_sartajac_H2020.rtp fout=TEST_RTP/bad_nlte_sartajac_H2020_sarta.rtp'];
eval(sartaer)
cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/Test_ArbPLEVS_PBL/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% AHA should now work, for some UNKNOWN and UNCOMMUNICATED reason, I found incFTC_airs_jan25_H2020_pclsam.f -> incFTC_airs_l1c_p2022jul22_dev.f has
%%%  PARAMETER(COFNTE = .FALSE.) instead of .TRUE.
cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/
sartaer = ['!/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=sub6000clr_airs_g12_sartaH2020_rp.rtp fout=sub6000clr_airs_g12_sartaH2020_rp_NEW.rtp'];
eval(sartaer)
cd /home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/Test_ArbPLEVS_PBL/
