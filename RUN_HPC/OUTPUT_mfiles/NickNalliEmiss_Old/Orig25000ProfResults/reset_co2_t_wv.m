addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools  % for rtpread & rtpwrite
addpath /asl/matlib/rtptools % for bits2pfields & pfields2bits
addpath /asl/matlib/aslutil  % for rad2bt

addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/SHOWSTATS
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/QUICKTASKS_TELECON/SARTA_Tuning/Scott_ResetCode/Sergio_resets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iFix = -1; %% orig
iFix = +1; %% new
if iFix < 0
  %% I think this is what we SHOULD use since Scott developed reset code using Masuda
  fin = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis.rtp';
  fout = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis_resetTWV.rtp';

  %% after analyzing fin, and resetting CO2/Tsurf/WV ... save the params to nicknalli
  finx  = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';
  foutx = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';
else
  %% hmm maybe we should do it this way after all since the reset code is agnostic about emissivity ....
  fin  = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';
  fout = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';

  %% after analyzing fin, and resetting CO2/Tsurf/WV ... save the params to nicknalli
  finx  = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';
  foutx = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';
end

fip = mktempS('.ip.rtp');
fop = mktempS('.op.rtp');
frp = mktempS('.rp.rtp');

%sarta1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_tra';                               %%% WIERD
sarta1 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v2';
sarta2 = '/home/sergio/SARTA_CLOUDY/BinV201/sarta_apr08_m140x_iceGHMbaum_waterdrop_desertdust_slabcloud_hg3';
sarta3 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';
sarta4 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_may19_prod';                                          %%% WIERD

sarta5 = '/home/sergio/SARTA_CLOUDY/LatestSarta/sarta/bin/airs_l1c_2834_cloudy_may19_prod_tuning';  %% VARIABLE TUNING

sarta = sarta5;

klayers = '/asl/packages/klayers/Bin/klayers_airs';
klayers = '/asl/packages/klayersV205/BinV201/klayers_airs';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hin,ha,pin,pa] = rtpread(fin);
if ~isfield(pin,'rcalc')
  if isfield(pin,'rclr')
    disp('found p.rclr')
    pin.rcalc = pin.rclr;
  end
end
BTobs = rad2bt(hin.vchan,pin.robs1);
BTcal0 = rad2bt(hin.vchan,pin.rcalc);

tic
iDoCO2 = +1;
if iDoCO2 > 0
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
  disp('running reset co2 ....')
  [hco2,ha,pco2,pa] = reset_co2ppm(fin);
  toc
else
  disp('skip reset co2')
  hco2 = hin;
  pco2 = pin;
end
rtpwrite(fip,hco2,ha,pco2,pa);

iDoCH4 = +1;
if iDoCH4 > 0
  disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
  disp('running reset ch4 ....')
  hch4   = hco2;
  [hch4,pch4] = reset_ch4ppm(hco2,pco2);
  toc
else
  disp('skip reset ch4')
  hch4 = hco2;
  pch4 = pco2;
end
rtpwrite(fip,hch4,ha,pch4,pa);

disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
[hout,ha,pout,pa,deltastruct] = reset_Tsurf_water(fip,[],sarta);
BTcalF = rad2bt(hout.vchan,pout.rcalc);
toc

figure(1); clf
plot(hout.vchan,nanmean(BTobs'-BTcal0'),'b',hout.vchan,nanmean(BTobs'-BTcalF'),'r')
plotaxis2; title('Bias before and after resetST/WV')
figure(2); clf
plot(hout.vchan,nanstd(BTobs'-BTcal0'),'b',hout.vchan,nanstd(BTobs'-BTcalF'),'r')
plotaxis2; title('Std before and after resetST/WV')

figure(3); clf
plot(0:0.001:2,histc(deltastruct.wmult,0:0.001:2),'b',-4:0.1:+4,histc(deltastruct.dstemp,-4:0.1:+4),'r')
hl= legend('WV mult','ST offset','location','best');

rtpwrite(fout,hout,ha,pout,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(fin,finx)
  %% now fix nalli
  %% now read in finx and reset the stemp,WV
  [hx,hax,px,pax] = rtpread(finx);
  px.gas_1 = pout.gas_1;
  px.gas_2 = pout.gas_2;
  px.gas_6 = pout.gas_6;
  px.stemp = pout.stemp;
  rtpwrite(foutx,hx,hax,px,pax);
end
toc
