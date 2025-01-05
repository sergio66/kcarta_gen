addpath /asl/matlib/h4tools

fM = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis.rtp';    %% masuda emis
fN = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';           %% new emis
fN = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';  %% new emis, reset ST, CO2,CH4, WV

[h,ha,pM,pa] = rtpread(fM);
[h,ha,pN,pa] = rtpread(fN);

[sum(pM.gas_1(:)-pN.gas_1(:)) sum(pM.gas_2(:)-pN.gas_2(:)) sum(pM.gas_3(:)-pN.gas_3(:))]
[sum(pM.ptemp(:)-pN.ptemp(:)) sum(pM.stemp(:)-pN.stemp(:))]
[sum(pM.satzen(:)-pN.satzen(:)) sum(pM.wspeed(:)-pN.wspeed(:)) sum(pM.rlat(:)-pN.rlat(:)) sum(pM.rlon(:)-pN.rlon(:))]

pNick = struct;
pNick.satzen = pM.satzen;
pNick.wspeed = pM.wspeed;
pNick.stemp = pM.stemp;
pNick.rlon = pN.rlon;
pNick.rlat = pN.rlat;
pNick.masuda.nemis = ones(size(pN.rlon))*19;
pNick.masuda.efreq = pM.efreq(1:19,:);
pNick.masuda.emis  = pM.emis(1:19,:);
pNick.nalli.nemis = ones(size(pN.rlon))*30;
pNick.nalli.efreq = pN.efreq(1:30,:);
pNick.nalli.emis  = pN.emis(1:30,:);

comment1 = 'satzen  windspeed (m/s)  SST   rlon  rlat';
comment2 = 'masuda : number of points (19) at wavenumber (efreq) and value (emis)';
comment3 = 'nalli  : number of points (19) at wavenumber (efreq) and value (emis)';
save /asl/ftp/pub/sergio/nicknalli.mat pNick comment*
