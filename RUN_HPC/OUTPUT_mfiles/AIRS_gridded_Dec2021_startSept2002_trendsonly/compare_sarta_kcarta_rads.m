addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil

[h,ha,p,pa] = rtpread('/home/sergio/MATLABCODE/oem_pkg_run/FIND_NWP_MODEL_TRENDS/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

raaRad = zeros(2834,4608);
for ii = 1 : 4608
  fname = ['AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fname);
  raaRad(:,ii) = x.rKc;
end

raaRad = raaRad(ichan,:);
fKc = x.fKc(ichan);
tS  = rad2bt(fKc,p.rcalc);
tKC = rad2bt(fKc,raaRad);

plot(fKc,nanmean(tKC'-tS'),'b',fKc,nanstd(tKC'-tS'),'c')
xlim([640 1640]); hl = legend('mean','std','location','best');
  xlabel('wavenumber cm-1'); ylabel('\delta BT(K)'); title('kCARTA-SARTA');
addpath /home/sergio/MATLABCODE
plotaxis2;
