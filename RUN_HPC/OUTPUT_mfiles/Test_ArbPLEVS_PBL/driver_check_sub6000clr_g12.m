%% see /asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/make_clear_rtp_039.m, copied here

[h,ha,psarta_clr_rad_airs,pa] = rtpread('/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_rp.rtp');
tsarta_clr_rad_airs = rad2bt(h.vchan,psarta_clr_rad_airs.rcalc);

scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.spres); colorbar
scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.satzen); colorbar
scatter_coast(psarta_clr_rad_airs.rlon,psarta_clr_rad_airs.rlat,10,psarta_clr_rad_airs.stemp-tsarta_clr_rad_airs(1520,:)); colorbar

p = psarta_clr_rad_airs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nprof = length(psarta_clr_rad_airs.stemp);

if ~exist('6555_profiles/PBL_G12/pbl_G12_summary.mat')
  use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_pbl_g12_op.rtp';
  h2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');

  pbl_kcarta = nan(2834,nprof);
  pbl_found = zeros(1,nprof);
  for ii = 1 : nprof
    if mod(ii,1000) == 0
      fprintf(1,'+')
    elseif mod(ii,100) == 0
      fprintf(1,'.')
    end
    fin = ['6555_profiles/PBL_G12/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    if exist(fin)
      x = load(fin);
      pbl_kcarta(:,ii) = x.rKc;
      pbl_found(ii) = +1;
    end
  end
  fKc = x.fKc;
  fprintf(1,'found %5i of the PBL files \n',sum(pbl_found))
  save 6555_profiles/PBL_G12/pbl_G12_summary.mat pbl_found pbl_kcarta fKc h2645

  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),'b.',psarta_clr_rad_airs.stemp,rad2bt(1231,pbl_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231');
  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:) - rad2bt(1231,pbl_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231 sarta-PBL kcarta');
else
  load('6555_profiles/PBL_G12/pbl_G12_summary.mat');
end

ichan = h2645.h.ichan;
tpbl = rad2bt(fKc(ichan),pbl_kcarta(ichan,:));
day   = find(p.solzen < 90);
night = find(p.solzen >= 90);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'k',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'g',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,day)'-tpbl(:,day)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,day)'-tpbl(:,day)'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,night)'-tpbl(:,night)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,night)'-tpbl(:,night)'),'m'); plotaxis2;
  title('PBL : thicker layers at TOA \newline are hurting O3 and NLTE','fontsize',10)
  hl = legend('bias all','std dev all','bias day','std dev day','bias night','std dev night','location','best','fontsize',10);
  title('PBL G12')

disp('ret to continue'); pause

xlim([645 2780])
xlim([645 1680])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('6555_profiles/AIRS_G12/airs_G12_summary.mat')
  use_this_rtp = '/home/sergio/MATLABCODE/QUICKTASKS_TELECON/PBL_BillIrion_klayers/PROFILES/sub6000clr_airs_g12_op.rtp';
  h2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');

  airs_kcarta = nan(2834,nprof);
  airs_found = zeros(1,nprof);
  for ii = 1 : nprof
    if mod(ii,1000) == 0
      fprintf(1,'+')
    elseif mod(ii,100) == 0
      fprintf(1,'.')
    end
    fin = ['6555_profiles/AIRS_G12/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    if exist(fin)
      x = load(fin);
      airs_kcarta(:,ii) = x.rKc;
      airs_found(ii) = +1;
    end
  end
  fKc = x.fKc;
  fprintf(1,'found %5i of the AIRS files \n',sum(airs_found))
  save 6555_profiles/AIRS_G12/airs_G12_summary.mat airs_found airs_kcarta fKc h2645
else
  load('6555_profiles/AIRS_G12/airs_G12_summary.mat');
end
plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),'b.',psarta_clr_rad_airs.stemp,rad2bt(1231,airs_kcarta(1291,:)),'r.')
  xlabel('ECM stemp'); ylabel('BT1231');
plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:) - rad2bt(1231,airs_kcarta(1291,:)),'r.')
  xlabel('ECM stemp'); ylabel('BT1231 sarta-AIRS kcarta');

ichan = h2645.h.ichan;
tairs = rad2bt(fKc(ichan),airs_kcarta(ichan,:));
day   = find(p.solzen < 90);
night = find(p.solzen >= 90);
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'k',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'g',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,day)'-tairs(:,day)'),'b',h.vchan,nanstd(tsarta_clr_rad_airs(:,day)'-tairs(:,day)'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs(:,night)'-tairs(:,night)'),'r',h.vchan,nanstd(tsarta_clr_rad_airs(:,night)'-tairs(:,night)'),'m'); plotaxis2;
  title('AIRS : thicker layers at TOA \newline are hurting O3 and NLTE','fontsize',10)
  hl = legend('bias all','std dev all','bias day','std dev day','bias night','std dev night','location','best','fontsize',10);
  title('AIRS G12')

disp('ret to continue'); pause

xlim([645 2780])
xlim([645 1680])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_check_sub6000clr

