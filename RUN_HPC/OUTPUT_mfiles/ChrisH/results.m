addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

[h,ha,p,pa] = rtpread('/home/chepplew/data/sarta/validation/sng_2020_subs_for_kcarta.rtp');    %% Vers1
[h,ha,p,pa] = rtpread('/home/chepplew/data/sarta/validation/sng_2020_subs_for_kcarta_v2.rtp'); %% Vers2

tObs = rad2bt(h.vchan,p.robs1);

for ii = 1 : length(p.stemp)
  fname = ['Vers2/H2016/Clear/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = load(fname);
  raaClear2016(:,ii) = loader.rKc;;

  fname = ['Vers2/H2016/Cloud/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = load(fname);
  raaCloud2016(:,ii) = loader.rKc;;

  fname = ['Vers2/H2020/Clear/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  loader = load(fname);
  raaClear2020(:,ii) = loader.rKc;;
end

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

raaClear2016 = raaClear2016(ichan,:);
raaCloud2016 = raaCloud2016(ichan,:);
raaClear2020 = raaClear2020(ichan,:);

fKc = loader.fKc(ichan);

taaClear2016 = rad2bt(fKc,raaClear2016);
taaClear2020 = rad2bt(fKc,raaClear2020);
taaCloud2016 = rad2bt(fKc,raaCloud2016);

plot(fKc,nanmean(tObs'-taaClear2020'),fKc,nanmean(tObs'-taaClear2016'),fKc,nanmean(tObs'-taaCloud2016'))
plot(fKc,nanmean(taaClear2016'-taaClear2020'),fKc,nanmean(taaClear2016'-taaCloud2016'))
plot(fKc,nanmean(taaClear2016'),fKc,nanmean(taaClear2020'),fKc,nanmean(taaCloud2016'))

