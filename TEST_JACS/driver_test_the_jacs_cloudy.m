addpath /home/sergio/git/matlabcode/matlibSergio/matlab2012/h4tools
addpath /home/sergio/git/matlabcode/matlibSergio/matlab2012/rtptools
addpath /home/sergio/git/matlabcode/matlibSergio/matlab2012/aslutil
addpath /home/sergio/git/matlabcode/TIME


%%% CLEAR CLEAR CLEAR : see driver_test_the_jacs.m
[h,ha,p,pa] = rtpread('ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp');
[yy,mm,dd,hh] = tai2utcSergio(p.rtime(1)); fprintf(1,'%4i/%02i/%02i at %8.5f hrs \n',yy,mm,dd,hh);      %%% ---> 2023/4/2
tclr = rad2bt(h.vchan,p.rcalc);

%% we got this from the /home/sergio/asl/asl/rtp/airs/airs_l1c_v675/clear/2023/ecmwf_airicrad_day092_clear.rtp

and subset for ones with spres == 1000 mb (so we know we would have fractional bottom layer) which turned out to be FOV 4758, and make emiss = 1, rho = 0

%% then just run sarta without jacs to make rtp file for all 7 profiles (so you can check finite difference jacs)
%% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now add in a modicum of clouds to test PCLSAM jacobians
p.ctype  = ones(size(p.stemp)) * 201;
p.cprtop = ones(size(p.stemp)) * 202;
p.cprbot = ones(size(p.stemp)) * 299;
p.cpsize = ones(size(p.stemp)) * 50;
p.cngwat = ones(size(p.stemp)) * 10;
p.cfrac  = ones(size(p.stemp)) * 0.25;

p.ctype2  = ones(size(p.stemp)) * 101;
p.cprtop2 = ones(size(p.stemp)) * 805;
p.cprbot2 = ones(size(p.stemp)) * 900;
p.cpsize2 = ones(size(p.stemp)) * 20;
p.cngwat2 = ones(size(p.stemp)) * 10;
p.cfrac2  = ones(size(p.stemp)) * 0.25;
p.cfrac12  = ones(size(p.stemp)) * 0.125;

%{
%% to test easy clear sky jacs
p.cfrac = p.cfrac * 0;   p.cfrac2 = p.cfrac2 * 0;   p.cfrac12 = p.cfrac12 * 0;
p.cngwat = p.cngwat * 0; p.cngwat2 = p.cngwat2 * 0; 
%}

rtpwrite('ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.op.rtp',h,ha,p,pa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% to do SARTA jacs, with output radiance usunits = rad use listp=1 listj=-1 jacunit=0   (proflie 1, all jacs, output units)
%% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.sarta.rtp listp=1 listj=-1 jacunit=0

addpath ../MATLAB

[fsarta,tjacsarta]  = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacTZ',100); tjacsarta  = squeeze(tjacsarta);  tjacsarta = tjacsarta(:,1:97);   tjacsarta = fliplr(tjacsarta); %% 98 is surface temp
[fsarta,q1jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG1',1);   q1jacsarta = squeeze(q1jacsarta); q1jacsarta = q1jacsarta(:,1:97); q1jacsarta = fliplr(q1jacsarta);
[fsarta,q2jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG2',2);   q2jacsarta = squeeze(q2jacsarta); q2jacsarta = q2jacsarta(:,1:97); q2jacsarta = fliplr(q2jacsarta); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% then just run sarta without jacs to make rtp file for all 7 profiles (so you can check finite difference jacs)
%% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.sarta.rtp

[hy,hay,py,pay] = rtpread('ecmwf_airicrad_day092_cloud_unitemiss_profile_4758.sarta.rtp');
pycalcs = rad2bt(hy.vchan,py.rcalc);
plot(h.vchan,pycalcs - tclr)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dT = 0.10;
% iLay = 1;
% figure(1); plot(hy.vchan,(pycalcs(:,7)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
% figure(2); plot(hy.vchan,(pycalcs(:,4)-pycalcs(:,1))/dq4,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
% 
% iLay = 2;
% figure(1); plot(hy.vchan,(pycalcs(:,5)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
% figure(2); plot(hy.vchan,(pycalcs(:,2)-pycalcs(:,1))/dq5,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
% 
% iLay = 7;
% figure(1); plot(hy.vchan,(pycalcs(:,6)-pycalcs(:,1))/dT,'k',fc,tc22(:,iLay),'bx-', fc,tc18(:,iLay),'c', fsarta,tjacsarta(:,iLay),'r');  legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('Tjac');
% figure(2); plot(hy.vchan,(pycalcs(:,3)-pycalcs(:,1))/dq10,'k',fc,q1c22(:,iLay),'bx-',fc,q1c18(:,iLay),'c',fsarta,q1jacsarta(:,iLay),'r'); legend('finite diff','kc22','kc18','sarta','location','best');    xlim([650 1650]); title('WVTjac');
% 

pnew = py;
nx = 100;
dq5  = pnew.gas_1(nx+1-5,2)  - pnew.gas_1(nx+1-5,1);     dq5  = dq5/6.02214076e26;      % 7.898596e-08     profile 2, dT would be profile 5
dq10 = pnew.gas_1(nx+1-10,3) - pnew.gas_1(nx+1-10,1);    dq10 = dq10/6.02214076e26;     % 4.978980e-08     profile 3, dT would be profile 6
dq4  = pnew.gas_1(nx+1-4,4)  - pnew.gas_1(nx+1-4,1);     dq4  = dq4/6.02214076e26;      % 8.585113e-08     profile 4, dT would be profile 7

f = hy.vchan;

iLay = 2;
%% if jacunit = 0
figure(1); plot(f,(pycalcs(:,5)-pycalcs(:,1))/(dT),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(pycalcs(:,2)-pycalcs(:,1))/(dq5),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])
%% if jacunit = 1
figure(1); plot(f,(rad2bt(f,pycalcs(:,5))-rad2bt(f,pycalcs(:,1)) )/(dT),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(rad2bt(f,pycalcs(:,2))-rad2bt(f,pycalcs(:,1)) )/(dq5),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])

iLay = 7;
%% if jacunit = 0
figure(1); plot(f,(pycalcs(:,6)-pycalcs(:,1))/(dT),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(pycalcs(:,3)-pycalcs(:,1))/(dq10),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])
%% if jacunit = 1
figure(1); plot(f,(rad2bt(f,pycalcs(:,6))-rad2bt(f,pycalcs(:,1)) )/(dT),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(rad2bt(f,pycalcs(:,3))-rad2bt(f,pycalcs(:,1)) )/(dq10),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])

iLay = 1;
%% if jacunit = 0
figure(1); plot(f,(pycalcs(:,7)-pycalcs(:,1))/(dT/2),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(pycalcs(:,4)-pycalcs(:,1))/(dq4/2),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])
%% if jacunit = 1
figure(1); plot(f,(rad2bt(f,pycalcs(:,7))-rad2bt(f,pycalcs(:,1)) )/(dT/2),'k',fsarta,tjacsarta(:,iLay),'r'); title(['T jac' num2str(iLay,'%02i')])
figure(2); plot(f,(rad2bt(f,pycalcs(:,4))-rad2bt(f,pycalcs(:,1)) )/(dq4/2),'k',fsarta,q1jacsarta(:,iLay),'r'); title(['Q jac' num2str(iLay,'%02i')])

