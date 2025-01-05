addpath /asl/matlib/aslutil

thedir = ['../'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  rad_airs_1013_22deg(:,ii) = a.rKc;
end

fairs = a.fKc;
new_bt_airs_1013_22deg = rad2bt(fairs,rad_airs_1013_22deg);

if exist('kcartaexec_oct08_2024.mat')
  kcartaexec_oct08_2024 = load('kcartaexec_oct08_2024.mat');
  boo = kcartaexec_oct08_2024.new_bt_airs_1013_22deg - new_bt_airs_1013_22deg;
  mboo = nanmean(boo(:));
  sboo = nanstd(boo(:));
  fprintf(1,'mean diff with kcartaexec_oct08_2024 = %8.6f +/- %8.6f \n',mboo,sboo);
  plot(fairs,nanmean(boo,2),fairs,nanstd(boo,[],2),'c');
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
title('KCARTA-KCARTA Oct08/2024 : 22 deg'); disp('ret to continue'); pause
  
end

[h,ha,p22,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/regr49_airs_1013_22deg.rp.rtp');
[h,ha,p00,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/regr49_airs_1013.rp.rtp');
t22 = rad2bt(h.vchan,p22.rcalc);
t00 = rad2bt(h.vchan,p00.rcalc);
plot(fairs(1:2378),nanmean(new_bt_airs_1013_22deg(1:2378,:)' - t22'),fairs(1:2378),nanstd(new_bt_airs_1013_22deg(1:2378,:)' - t22'),'c');
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
title('KCARTA-SARTA : 22 deg'); disp('ret to continue'); pause

plot(fairs(1:2378),nanmean(new_bt_airs_1013_22deg(1:2378,:)' - t22'),fairs(1:2378),nanstd(new_bt_airs_1013_22deg(1:2378,:)' - t22'),'c',...
     fairs(1:2378),new_bt_airs_1013_22deg(1:2378,1) - t22(:,1),'r'); plotaxis2; 
legend('bias','std dev','TRP only','location','best','fontsize',10);
title('KCARTA-SARTA : 22 deg'); disp('ret to continue'); pause

