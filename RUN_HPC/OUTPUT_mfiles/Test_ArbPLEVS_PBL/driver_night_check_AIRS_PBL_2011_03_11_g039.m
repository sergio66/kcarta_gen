%% see /asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/make_clear_rtp_039.m, copied here

[h,ha,psarta_clr_rad_airs,pa] = rtpread('/asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/night_clear_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_airs.rp.rtp');
tsarta_clr_rad_airs = rad2bt(h.vchan,psarta_clr_rad_airs.rcalc);

scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.spres); colorbar
scatter(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),10,psarta_clr_rad_airs.satzen); colorbar
scatter_coast(psarta_clr_rad_airs.rlon,psarta_clr_rad_airs.rlat,10,psarta_clr_rad_airs.stemp-tsarta_clr_rad_airs(1520,:)); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('NIGHT_PBL_2011_03_11_g039/night_pbl_summary.mat')
  use_this_rtp = '/umbc/xfs2/strow/asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/night_clear_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_pbl.op.rtp';
  h2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');

  pbl_kcarta = nan(2834,1215);
  pbl_found = zeros(1,1215);
  for ii = 1 : 1215
    if mod(ii,1000) == 0
      fprintf(1,'+')
    elseif mod(ii,100) == 0
      fprintf(1,'.')
    end
    fin = ['NIGHT_PBL_2011_03_11_g039/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    if exist(fin)
      x = load(fin);
      pbl_kcarta(:,ii) = x.rKc;
      pbl_found(ii) = +1;
    end
  end
  fKc = x.fKc;
  fprintf(1,'found %5i of the PBL files \n',sum(pbl_found))
  save NIGHT_PBL_2011_03_11_g039/night_pbl_summary.mat pbl_found pbl_kcarta fKc h2645

  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),'b.',psarta_clr_rad_airs.stemp,rad2bt(1231,pbl_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231');
  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:) - rad2bt(1231,pbl_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231 sarta-PBL kcarta');

  ichan = h2645.h.ichan;
  tpbl = rad2bt(fKc(ichan),pbl_kcarta(ichan,:));
  plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
    title('PBL : thicker layers at TOA \newline are hurting O3','fontsize',10)
  xlim([645 2780])
  xlim([645 1680])
else
  load('NIGHT_PBL_2011_03_11_g039/night_pbl_summary.mat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('NIGHT_AIRS_2011_03_11_g039/night_airs_summary.mat')
  use_this_rtp = '/umbc/xfs2/strow/asl/s1/sergio/rtp/rtp_airicrad_v6/2011/03/11/night_clear_airs_l1c_test_2011.03.11.039_cumsum_-1_noclds_airs.op.rtp';
  h2645 = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/h2645structure.mat');

  airs_kcarta = nan(2834,1215);
  airs_found = zeros(1,1215);
  for ii = 1 : 1215
    if mod(ii,1000) == 0
      fprintf(1,'+')
    elseif mod(ii,100) == 0
      fprintf(1,'.')
    end
    fin = ['NIGHT_AIRS_2011_03_11_g039/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    if exist(fin)
      x = load(fin);
      airs_kcarta(:,ii) = x.rKc;
      airs_found(ii) = +1;
    end
  end
  fKc = x.fKc;
  fprintf(1,'found %5i of the AIRS files \n',sum(airs_found))
  save NIGHT_AIRS_2011_03_11_g039/night_airs_summary.mat airs_found airs_kcarta fKc h2645

  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:),'b.',psarta_clr_rad_airs.stemp,rad2bt(1231,airs_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231');
  plot(psarta_clr_rad_airs.stemp,tsarta_clr_rad_airs(1520,:) - rad2bt(1231,airs_kcarta(1291,:)),'r.')
    xlabel('ECM stemp'); ylabel('BT1231 sarta-AIRS kcarta');

  ichan = h2645.h.ichan;
  tairs = rad2bt(fKc(ichan),airs_kcarta(ichan,:));
  plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
    title('AIRS100 : still have wierd O3; no NLTE','fontsize',10)
  xlim([645 2780])
  xlim([645 1680])
else
  load('NIGHT_AIRS_2011_03_11_g039/night_airs_summary.mat');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ichan = h2645.h.ichan;

tpbl = rad2bt(fKc(ichan),pbl_kcarta(ichan,:));
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c'); plotaxis2;
  title('PBL : thicker layers at TOA \newline are hurting O3','fontsize',10)
disp('ret to continue'); pause; 

tairs = rad2bt(fKc(ichan),airs_kcarta(ichan,:));
plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
  title('AIRS100 : still have wierd O3','fontsize',10)
disp('ret to continue'); pause; 

plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',...
     h.vchan,nanmean(tsarta_clr_rad_airs'-tairs'),'r',h.vchan,nanstd(tsarta_clr_rad_airs'-tairs'),'m'); plotaxis2;
legend('PBL bias','PBL stddev','AIRS100 bias','AIRS100 stddev','location','best','fontsize',10)
disp('ret to continue'); pause; 

plot(h.vchan,nanmean(tairs'-tpbl'),'k',h.vchan,nanstd(tairs'-tpbl'),'g'); plotaxis2;
  title('AIRS100 - PBL');
disp('ret to continue'); pause; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
l2s = load('/home/sergio/KCARTA/L2SComparisons/l2s_kc122_H20_605_2830.mat');
l2s_gid = l2s.iaGasID;
[~,l2s] = convolve_airs(double(l2s.w),double(l2s.d),double(ichan));

%for ig = 1 : 72; 
for ig = 3; 
  plot(h.vchan,nanmean(tsarta_clr_rad_airs'-tpbl'),'b',h.vchan,nanstd(tsarta_clr_rad_airs'-tpbl'),'c',h.vchan,exp(-l2s(:,ig))); 
  plotaxis2; title(num2str(l2s_gid(ig)));
  pause
end
%}
