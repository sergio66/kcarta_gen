addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/PLOTTER

[h,ha,p,pa] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/cloudy_airs_l1c_ecm_sarta_baum_ice.2011.03.11.039_cumsum_-1.rtp');
scatter_coast(p.rlon,p.rlat,50,p.stemp)

p.cfrac   = 0 * p.cfrac;
p.cfrac2  = 0 * p.cfrac2;
p.cfrac12 = 0 * p.cfrac12;
p.cngwat  = 0 * p.cngwat;
p.cngwat2 = 0 * p.cngwat2;
p.ctype  = -999 * ones(size(p.ctype));
p.ctype2 = -999 * ones(size(p.ctype2));

p.rcalc = p.sarta_rclearcalc; p = rmfield(p,'sarta_rclearcalc');
rtpwrite('/asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/cloudy_airs_l1c_ecm_sarta_baum_ice.2011.03.11.039_cumsum_-1_noclds.ip.rtp',h,ha,p,pa);

cd /asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/

klayers = ['/asl/packages/klayersV205/BinV201/klayers_airs'];
klayerser = ['!' klayers ' fin=cloudy_airs_l1c_ecm_sarta_baum_ice.2011.03.11.039_cumsum_-1_noclds.ip.rtp '];
klayerser = [klayerser   ' fout=cloudy_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_airs.op.rtp']; 
eval(klayerser)

%% see /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/clustbatch_redo_stemp_wv_cloud_filelist.m
toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';                                          
              %% HITRAN 2016 spectroscopy NEW DEFAULT CHRIS, odder results? SLOW IGNORE THIS TOOOOOO SLOW
toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';   
              %% HITRAN 2008 spectroscopy OLD DEFAULT SCOTT, better results, with WV line tuning!!!!
toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v2';
              %% HITRAN 2016 spectroscopy NEW DEFAULT CHRIS, odder results? FAST, wrong data location
toptsSARTA.sarta_airs = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';       
              %% HITRAN 2016 spectroscopy NEW DEFAULT CHRIS, odder results? FAST, correct data location
toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod_debug_save';
             %% HITRAN 2020 spectroscopy with jacobians nanananana
toptsSARTA.sarta_airs = '/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_may19_prod';
              %% HITRAN 2020 spectroscopy with jacobians nanananana
sarta = toptsSARTA.sarta_airs;
sartaer = ['!' sarta ' fin=cloudy_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_airs.op.rtp']; 
sartaer = [sartaer   ' fout=cloudy_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_airs.rp.rtp']; 
eval(sartaer)

klayers = ['/home/chepplew/gitLib/klayersV205/BinV201/klayers_pbl_wetwater_test'];
klayerser = ['!' klayers ' fin=cloudy_airs_l1c_ecm_sarta_baum_ice.2011.03.11.039_cumsum_-1_noclds.ip.rtp '];
klayerser = [klayerser   '  fout=cloudy_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_pbl.op.rtp']; 
eval(klayerser)
