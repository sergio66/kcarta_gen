addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/h4tools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/rtptools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/aslutil
addpath /home/sergio/git/matlabcode/JPL_DUST_Nov2014/GEOPHYSICAL/VER_AUG2014  %% subset_rtp_allcloudfields.m

if ~exist('pycalcs')
  load convolved_kcarta_vs_sarta.mat
end  

figure(1); plot(fc,rad2bt(fc,rc22),'b',fc,rad2bt(fc,rc18),'c',hy.vchan,rad2bt(hy.vchan,py.rcalc(:,1)),'r'); xlim([645 1645])
disp('ret to continue'); pause

isarta = find(fsarta >= 1231,1); ikcarta = find(fc >= 1231,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('Tz jacs');
%figure(2); plot(xq1c22(ikcarta,1:97),1:97,'bx-',xq1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta WV jac has continuum, kcarta only WV lines');
figure(2); plot(q1c22(ikcarta,1:97),1:97,'bx-',q1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta and kcarta WV jac have continuum');
disp('showed 1231 cm-1 ret to continue'); pause

isarta = find(fsarta >= 1419,1); ikcarta = find(fc >= 1419,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('Tz jacs'); ylim([20 60])
%figure(2); plot(xq1c22(ikcarta,1:97),1:97,'bx-',xq1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta WV jac has continuum, kcarta only WV lines'); ylim([20 60])
figure(2); plot(q1c22(ikcarta,1:97),1:97,'bx-',q1c18(ikcarta,1:97),1:97,'c',q1jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('sarta and kcarta WV jac have continuum'); ylim([20 60])
disp('showed 1419 cm-1 ret to continue'); pause

isarta = find(fsarta >= 721,1); ikcarta = find(fc >= 721,1);
figure(1); plot(tc22(ikcarta,1:97),1:97,'bx-',tc18(ikcarta,1:97),1:97,'c',tjacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('Tz jacs');
figure(2); plot(q2c22(ikcarta,1:97),1:97,'bx-',q2c18(ikcarta,1:97),1:97,'c',q2jacsarta(isarta,1:97),1:97,'r'); legend('kc22','kc18','sarta','location','best'); title('CO2 jac'); ylim([20 60])
disp('showed 721 cm-1 ret to continue'); pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pycalcs = rad2bt(hy.vchan,py.rcalc);

iLay = 2;
figure(1); plot(fc,tc22(:,iLay),'b',fc,tc18(:,iLay),'c',fsarta,tjacsarta(:,iLay),'r');       legend('kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
%figure(2); plot(fc,xq1c22(:,iLay),'b',fc,xq1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('kc22','kc18','sarta','location','best'); xlim([650 1650]); title('WV lines (kc) and lines+cont (sarta) jac');
figure(2); plot(fc,q1c22(:,iLay),'b',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r');    legend('kc22','kc18','sarta','location','best'); xlim([650 1650]); title('WV + cont lines (kc and sarta)');
figure(3); plot(fc,q2c22(:,iLay),'b',fc,q2c18(:,iLay),'c',fsarta,q2jacsarta(:,iLay),'r');    legend('kc22','kc18','sarta','location','best'); xlim([650 1050]); title('CO2 jac');
disp('iLay - 2 ret to continue'); pause

dT = 0.10;
iLay = 1;
figure(1); plot(hy.vchan,(pycalcs(:,7)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
figure(2); plot(hy.vchan,(pycalcs(:,4)-pycalcs(:,1))/dq4,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
disp('iLay - 1 ret to continue'); pause

iLay = 2;
figure(1); plot(hy.vchan,(pycalcs(:,5)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
figure(2); plot(hy.vchan,(pycalcs(:,2)-pycalcs(:,1))/dq5,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
disp('iLay - 2 ret to continue'); pause

iLay = 7;
figure(1); plot(hy.vchan,(pycalcs(:,6)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
figure(2); plot(hy.vchan,(pycalcs(:,3)-pycalcs(:,1))/dq10,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
disp('iLay - 7 ret to continue'); pause

%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta
%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta
%% save convolved_kcarta_vs_sarta.mat fc *c22* *c18* dq* dT hy py fsarta tjacsarta q1jacsarta q2jacsarta
