addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/h4tools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/rtptools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/aslutil
addpath /home/sergio/git/matlabcode/JPL_DUST_Nov2014/GEOPHYSICAL/VER_AUG2014  %% subset_rtp_allcloudfields.m

if ~exist('pycalcs')
  load convolved_kcarta_vs_sarta.mat
end  

%{
% latest and greated overrwrite of sarta clear jacs from convolved_kcarta_vs_sarta.mat
% see driver_test_the_jacs.m
% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp LISTJ=-1 LISTP=1 JACUNIT=0
[fsarta,tjacsarta]  = readsarta_jac('/home/sergio/git/kcarta_gen/TEST_JACS/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacTZ',100); tjacsarta  = squeeze(tjacsarta);  tjacsarta = tjacsarta(:,1:97);   tjacsarta = fliplr(tjacsarta); %% 98 is surface temp
[fsarta,q1jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/TEST_JACS/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG1',1);   q1jacsarta = squeeze(q1jacsarta); q1jacsarta = q1jacsarta(:,1:97); q1jacsarta = fliplr(q1jacsarta);
[fsarta,q2jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/TEST_JACS/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG2',2);   q2jacsarta = squeeze(q2jacsarta); q2jacsarta = q2jacsarta(:,1:97); q2jacsarta = fliplr(q2jacsarta);
%% then after jacs do
% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp
[hy,ha,py,pa] = rtpread('/home/sergio/git/kcarta_gen/TEST_JACS/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp');
%}

figure(1); plot(fc,rad2bt(fc,rc22),'b',fc,rad2bt(fc,rc18),'c',hy.vchan,rad2bt(hy.vchan,py.rcalc(:,1)),'r'); xlim([645 1645])
disp('ret to continue'); pause

isarta = find(fsarta >= 1231,1); ikcarta = find(fc >= 1231,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('1231 cm-1 Tz jacs');
%figure(2); plot(xq1c22(ikcarta,1:97),1:97,'bx-',xq1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta WV jac has continuum, kcarta only WV lines');
figure(2); plot(q1c22(ikcarta,1:97),1:97,'bx-',q1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('1231 cm-1 sarta and kcarta WV jac have continuum');
disp('showed 1231 cm-1 ret to continue'); pause

isarta = find(fsarta >= 1419,1); ikcarta = find(fc >= 1419,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('1419 cm-1 Tz jacs'); ylim([20 60])
%figure(2); plot(xq1c22(ikcarta,1:97),1:97,'bx-',xq1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta WV jac has continuum, kcarta only WV lines'); ylim([20 60])
figure(2); plot(q1c22(ikcarta,1:97),1:97,'bx-',q1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('1419 cm-1 sarta and kcarta WV jac have continuum'); ylim([20 60])
disp('showed 1419 cm-1 ret to continue'); pause

isarta = find(fsarta >= 721,1); ikcarta = find(fc >= 721,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('721 cm-1 Tz jacs');
figure(2); plot(q2c22(ikcarta,1:97),1:97,'bx-',q2c18(ikcarta,1:97),1:97,'c',q2jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('721 cm-1 CO2 jac'); ylim([20 60])
disp('showed 721 cm-1 ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pycalcs = rad2bt(hy.vchan,py.rcalc);
rcalcs  = py.rcalc;
py.plays = plevs2plays(py.plevs);

% iLay = 2;
% figure(1); plot(fc,tc22(:,iLay),'b',fc,tc18(:,iLay),'c',fsarta,tjacsarta(:,iLay),'r');       legend('kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
% %figure(2); plot(fc,xq1c22(:,iLay),'b',fc,xq1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('kc22','kc18','sarta','location','best'); xlim([650 1650]); title('WV lines (kc) and lines+cont (sarta) jac');
% figure(2); plot(fc,q1c22(:,iLay),'b',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r');    legend('kc22','kc18','sarta','location','best'); xlim([650 1650]); title('WV + cont lines (kc and sarta)');
% figure(3); plot(fc,q2c22(:,iLay),'b',fc,q2c18(:,iLay),'c',fsarta,q2jacsarta(:,iLay),'r');    legend('kc22','kc18','sarta','location','best'); xlim([650 1050]); title('CO2 jac');
% disp('iLay = 2 ret to continue'); pause

kJacobOutput = -1
if kJacobOutput == 0
  disp('the jacobians are in dBT/dX ie everything converted to BT(K)')
  dT = 0.10;
  iLay = 1;
  figure(1); plot(hy.vchan,(pycalcs(:,7)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 1');
  figure(2); plot(hy.vchan,(pycalcs(:,4)-pycalcs(:,1))/dq4,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 1');
  disp('iLay = 1 ret to continue'); pause
  
  iLay = 2;
  figure(1); plot(hy.vchan,(pycalcs(:,5)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 2');
  figure(2); plot(hy.vchan,(pycalcs(:,2)-pycalcs(:,1))/dq5,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 2');
  disp('iLay = 2 ret to continue'); pause
  
  iLay = 7;
  figure(1); plot(hy.vchan,(pycalcs(:,6)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 7');
  figure(2); plot(hy.vchan,(pycalcs(:,3)-pycalcs(:,1))/dq10,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 7');
  disp('iLay = 7 ret to continue'); pause

elseif kJacobOutput == -1
  disp('the jacobians are in drad/dX ie everything still in radiances')
  dT = 0.10;
  iLay = 1;
  figure(1); plot(hy.vchan,(rcalcs(:,7)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 1');
  figure(2); plot(hy.vchan,(rcalcs(:,4)-rcalcs(:,1))/dq4,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 1');
  disp('iLay = 1 ret to continue'); pause
  
  iLay = 2;
  figure(1); plot(hy.vchan,(rcalcs(:,5)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 2');
  figure(2); plot(hy.vchan,(rcalcs(:,2)-rcalcs(:,1))/dq5,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 2');
  disp('iLay = 2 ret to continue'); pause
  
  iLay = 7;
  figure(1); plot(hy.vchan,(rcalcs(:,6)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac, Lay = 7');
  figure(2); plot(hy.vchan,(rcalcs(:,3)-rcalcs(:,1))/dq10,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVjac, Lay = 7');
  disp('iLay = 7 ret to continue'); pause
end

%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta
%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta
%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iTile = +1;
if iTile > 0

  %% tile plots
  %want a 2x3 tiled layout
  figure(4); clf
  ta = tiledlayout(1,2);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  isarta1 = find(fsarta >= 1231,1); ikcarta1 = find(fc >= 1231,1);
  isarta2 = find(fsarta >= 1419,1); ikcarta2 = find(fc >= 1419,1);
  theplays = py.plays(1:97,1); theplays = flipud(theplays);
  
  tafov(1) = nexttile;
  semilogy(tc22(ikcarta1,1:96),theplays(1:96),'bx-',tjacsarta(isarta1,1:96),theplays(1:96),'rx-',tc22(ikcarta2,1:96),theplays(1:96),'cx-',tjacsarta(isarta2,1:96),theplays(1:96),'m');
  legend('kc 1231','sarta 1231','kcarta 1419','sarta 1419','location','best'); title('T(p) jacs'); set(gca,'ydir','reverse'); ylim([10 1000])
  xlabel('d(rad)/dT'); ylabel('Pressure [mb]')
  
  tafov(2) = nexttile;
  plot(q1c22(ikcarta1,1:96)/1e5,theplays(1:96),'bx-',q1jacsarta(isarta1,1:96)/1e5,theplays(1:96),'rx-',q1c22(ikcarta2,1:96)/1e8,theplays(1:96),'cx-',q1jacsarta(isarta2,1:96)/1e8,theplays(1:96),'m');
  legend('kc 1231','sarta 1231','kcarta 1419','sarta 1419','location','east'); title('normalized WV(p) jacs'); set(gca,'ydir','reverse'); ylim([100 1000])
  xlabel('d(rad)/dQ [normalized]'); ylabel('Pressure [mb]')
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(5); clf
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  iLay = 1;
  fprintf(1,'pressure for layer %2i = %8.3f \n',iLay,theplays(iLay))  
  tafov(1) = nexttile;  
    plot(hy.vchan,(rcalcs(:,7)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-',  fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]); title('1000 mb layer');
    ylabel('d(rad)/dT')
  tafov(2) = nexttile;
    plot(hy.vchan,(rcalcs(:,4)-rcalcs(:,1))/dq4,'k',fc,q1c22(:,iLay),'bx-',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]);
    ylabel('d(rad)/dQ'); xlabel('Wavenumber [cm-1]')
  
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  tafov(1).XTickLabel = '';
  tafov(1).XLabel.String = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(6); clf
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];
  
  iLay = 2;
  fprintf(1,'pressure for layer %2i = %8.3f \n',iLay,theplays(iLay))  
  tafov(1) = nexttile;  
    plot(hy.vchan,(rcalcs(:,5)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-',   fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]); title('970 mb layer');
    ylabel('d(rad)/dT')  
  tafov(2) = nexttile;    
    plot(hy.vchan,(rcalcs(:,2)-rcalcs(:,1))/dq5,'k',fc,q1c22(:,iLay),'bx-',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]); 
    ylabel('d(rad)/dQ'); xlabel('Wavenumber [cm-1]')
  
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  tafov(1).XTickLabel = '';
  tafov(1).XLabel.String = [];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(7); clf
  ta = tiledlayout(2,1);
  ta.OuterPosition = [0.0375 0.0375 0.925 0.925];

  iLay = 7;
  fprintf(1,'pressure for layer %2i = %8.3f \n',iLay,theplays(iLay))
  tafov(1) = nexttile;  
    plot(hy.vchan,(rcalcs(:,6)-rcalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-',   fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]); title('840 mb layer');
    ylabel('d(rad)/dT')  
  tafov(2) = nexttile;    
    plot(hy.vchan,(rcalcs(:,3)-rcalcs(:,1))/dq10,'k',fc,q1c22(:,iLay),'bx-',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','sarta','location','best');    xlim([650 1650]); 
    ylabel('d(rad)/dQ'); xlabel('Wavenumber [cm-1]')
  
  ta.Padding = 'compact';
  ta.TileSpacing = 'compact';
  tafov(1).XTickLabel = '';
  tafov(1).XLabel.String = [];

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dirReport = '/home/sergio/git/matlabcode/QUICKTASKS_TELECON/SuddenStratWarming_SSW/RTA_Report/';
  figure(4); sergioprintfig([dirReport '/t_and_q_vertical_profile_jacs_1231_1419_channels']);
  figure(5); sergioprintfig([dirReport '/spectral_t_q_jacs_1000_mb_layer']);
  figure(6); sergioprintfig([dirReport '/spectral_t_q_jacs_0970_mb_layer']);
  figure(7); sergioprintfig([dirReport '/spectral_t_q_jacs_0840_mb_layer']);  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end  
