addpath /asl/matlib/aslutil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thedir = ['ORIGCODE/Regr49_airs_1013'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  orig_rad_airs_1013_00deg(:,ii) = a.rKc;
end

thedir = ['ORIGCODE/Regr49_airs_1013_22deg'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  orig_rad_airs_1013_22deg(:,ii) = a.rKc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

thedir = ['ORIGCODE/Regr49_pbl_1013'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  orig_rad_pbl_1013_00deg(:,ii) = a.rKc;
end

thedir = ['ORIGCODE/Regr49_pbl_1013_22deg'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  orig_rad_pbl_1013_22deg(:,ii) = a.rKc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thedir = ['NEWCODE/Regr49_airs_1013'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  new_rad_airs_1013_00deg(:,ii) = a.rKc;
end

thedir = ['NEWCODE/Regr49_airs_1013_22deg'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  new_rad_airs_1013_22deg(:,ii) = a.rKc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

thedir = ['NEWCODE/Regr49_pbl_1013'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  new_rad_pbl_1013_00deg(:,ii) = a.rKc;
end

thedir = ['NEWCODE/Regr49_pbl_1013_22deg'];
for ii = 1 : 49
  fname = [thedir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  a = load(fname);
  new_rad_pbl_1013_22deg(:,ii) = a.rKc;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fairs = a.fKc;

orig_bt_airs_1013_00deg = rad2bt(fairs,orig_rad_airs_1013_00deg);
orig_bt_airs_1013_22deg = rad2bt(fairs,orig_rad_airs_1013_22deg);
orig_bt_pbl_1013_00deg  = rad2bt(fairs,orig_rad_pbl_1013_00deg);
orig_bt_pbl_1013_22deg  = rad2bt(fairs,orig_rad_pbl_1013_22deg);

new_bt_airs_1013_00deg = rad2bt(fairs,new_rad_airs_1013_00deg);
new_bt_airs_1013_22deg = rad2bt(fairs,new_rad_airs_1013_22deg);
new_bt_pbl_1013_00deg  = rad2bt(fairs,new_rad_pbl_1013_00deg);
new_bt_pbl_1013_22deg  = rad2bt(fairs,new_rad_pbl_1013_22deg);

%%%%%%%%%%%%%%%%%%%%%%%%%

%% these should be zero or close to zero
plot(fairs,nanmean(orig_bt_airs_1013_00deg'-new_bt_airs_1013_00deg'),'b',...
     fairs,nanstd(orig_bt_airs_1013_00deg'-new_bt_airs_1013_00deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('00 deg : AIRS 100 layers : orig - new'); disp('ret to continue'); pause

plot(fairs,nanmean(orig_bt_airs_1013_22deg'-new_bt_airs_1013_22deg'),'b',...
     fairs,nanstd(orig_bt_airs_1013_22deg'-new_bt_airs_1013_22deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('22 deg : AIRS 100 layers : orig - new'); disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% these should have issues, as orig code did arb PLEVS reference profiles wrongly
plot(fairs,nanmean(orig_bt_pbl_1013_00deg'-new_bt_pbl_1013_00deg'),'b',...
     fairs,nanstd(orig_bt_pbl_1013_00deg'-new_bt_pbl_1013_00deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('00 deg : PBL 100 layers : orig - new'); disp('ret to continue'); pause

plot(fairs,nanmean(orig_bt_pbl_1013_22deg'-new_bt_pbl_1013_22deg'),'b',...
     fairs,nanstd(orig_bt_pbl_1013_22deg'-new_bt_pbl_1013_22deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('22 deg : PBL 100 layers : orig - new'); disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

%% these should be very good, since code is fixed
%%       so now we are just showing layering diffs

plot(fairs,nanmean(new_bt_airs_1013_00deg'-new_bt_pbl_1013_00deg'),'b',...
     fairs,nanstd(new_bt_airs_1013_00deg'-new_bt_pbl_1013_00deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('00 deg : AIRS-PBL 100 layers : new - new'); disp('ret to continue'); pause

plot(fairs,nanmean(new_bt_airs_1013_22deg'-new_bt_pbl_1013_22deg'),'b',...
     fairs,nanstd(new_bt_airs_1013_22deg'-new_bt_pbl_1013_22deg'),'c')
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
  title('22 deg : AIRS-PBL 100 layers : new - new'); disp('ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%

[h,ha,p22,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/regr49_airs_1013_22deg.rp.rtp');
[h,ha,p00,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/regr49_airs_1013.rp.rtp');
t22 = rad2bt(h.vchan,p22.rcalc);
t00 = rad2bt(h.vchan,p00.rcalc);
plot(fairs(1:2378),nanmean(new_bt_airs_1013_22deg(1:2378,:)' - t22'),fairs(1:2378),nanstd(new_bt_airs_1013_22deg(1:2378,:)' - t22'),'c');
  plotaxis2; legend('bias','std dev','location','best','fontsize',10);
title('KCARTA-SARTA : 22 deg'); disp('ret to continue'); pause

plot(fairs(1:2378),new_bt_airs_1013_22deg(1:2378,1) - t22(:,1)); plotaxis2; 
title('Tropical : KCARTA-SARTA : 22 deg'); disp('ret to continue'); pause
