addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/h4tools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/rtptools
addpath ../../../WorkDirDec2025/matlabcode/matlibSergio/matlib/aslutil
addpath /home/sergio/git/matlabcode/JPL_DUST_Nov2014/GEOPHYSICAL/VER_AUG2014  %% subset_rtp_allcloudfields.m

[h,ha,p,pa] = rtpread('ecmwf_airicrad_day092_clear_unitemiss.op.rtp');
[hx,px] = subset_rtp_allcloudfields(h,p,[],[],4758);
[hnew,pnew] = replicate_rtp_headprof(hx,px,1,7);

nlevs = pnew.nlevs(1);  %%% 98, so nlays = 97; recall kProfLayer = 100

%%   RTP             |          INTERNAL KCARTA                    |   OUTPUT JAC
%% ------------------|---------------------------------------------|-----------
%% rtp layer 1 == 0.005 mb -->  kcarta 100   kProfLayer-1 + 1      |      97
%% rtp layer 2 == 0.015 mb -->  kcarta 99    kProfLayer-2 + 1      |      96
%% rtp layer 3 == 0.105 mb -->  kcarta 98    kProfLayer-3 + 1      |      95
%% ...               |                                             |
%% rtp layer 88      |      --> kcarta 13    kProfLayer-i + 1      |      10
%% rtp layer 89      |      --> kcarta 12    kProfLayer-i + 1      |       9 
%% rtp layer 90      |      --> kcarta 11    kProfLayer-i + 1      |       8
%% rtp layer 91      |      --> kcarta <10>  kProfLayer-i + 1      |       7     PERTURB this iJax = 10
%% rtp layer 92      |      --> kcarta 09    kProfLayer-i + 1      |       6
%% rtp layer 93      |      --> kcarta 08    kProfLayer-i + 1      |       5
%% rtp layer 94      |      --> kcarta 07    kProfLayer-i + 1      |       4
%% rtp layer 95      |      --> kcarta 06    kProfLayer-i + 1      |       3
%% rtp layer 96      |      --> kcarta <05>  kProfLayer-i + 1      |       2     PERTURB this iJax = 5
%% rtp layer 97      |      --> kcarta 04    kProfLayer-nlays + 1  |       1
%% ------------------|---------------------------------------------|-------------
%% rtp layer 98      |      --> kcarta 03    kProfLayer-i + 1      |      DNE
%% rtp layer 99      |      --> kcarta 02    kProfLayer-i + 1      |      DNE
%% rtp layer 100     |      --> kcarta 01    kProfLayer-i + 1      |      DNE
%% rtp layer 101     |      --> kcarta XX    kProfLayer-i + 1      |      DNE 
%% ------------------|--------------------------------------------------------

%% perturb water kcarta layers 5,10,5:10

nx = 100;  %% start with this
nx+1-4     %% 97  surface layer
nx+1-5     %% 96  two layers above surface
nx+1-10    %% 91  seven layers above surface

pnew.gas_1(nx+1-5,2)    = pnew.gas_1(nx+1-5,2) * 1.01;
pnew.gas_1(nx+1-10,3)   = pnew.gas_1(nx+1-10,3) * 1.01;
%pnew.gas_1(nx+1-(5:10),4) = pnew.gas_1(nx+1-(5:10),4) * 1.01;
pnew.gas_1(nx+1-4,4)    = pnew.gas_1(nx+1-4,4) * 1.01;

%% perturb T layers 5,10,5:10
ls -lt ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtppnew.ptemp(nx+1-5,5)    = pnew.ptemp(nx+1-5,5) + 0.1;
pnew.ptemp(nx+1-10,6)   = pnew.ptemp(nx+1-10,6) + 0.1;
%pnew.ptemp(nx+1-(5:10),7) = pnew.ptemp(nx+1-(5:10),7) + 0.1;
pnew.ptemp(nx+1-4,7) = pnew.ptemp(nx+1-4,7) + 0.1; 

q4 = pnew.gas_1(nx+1-4,1);   q4 = q4/6.02214076e26;
q5 = pnew.gas_1(nx+1-5,1);   q5 = q5/6.02214076e26;
q10 = pnew.gas_1(nx+1-10,1); q10 = q10/6.02214076e26;
dq5 = pnew.gas_1(nx+1-5,2)   - pnew.gas_1(nx+1-5,1);     dq5  = dq5/6.02214076e26;      % 7.898596e-08
dq10 = pnew.gas_1(nx+1-10,3) - pnew.gas_1(nx+1-10,1);    dq10 = dq10/6.02214076e26;     % 4.978980e-08
dq4 = pnew.gas_1(nx+1-4,4)   - pnew.gas_1(nx+1-4,1);     dq4  = dq4/6.02214076e26;      % 8.585113e-08

fprintf(1,'iL = %3i     q4   dq4  = %8.6e   %8.6e \n',nx+1-04,q4,dq4)
fprintf(1,'iL = %3i     q5   dq5  = %8.6e   %8.6e \n',nx+1-05,q5,dq5)
fprintf(1,'iL = %3i     q10  dq10 = %8.6e   %8.6e \n',nx+1-10,q10,dq10)

plot(pnew.ptemp(:,5:7)-pnew.ptemp(:,1),1:101,'o-');   set(gca,'ydir','reverse'); ylim([85 100]); grid on
plot(pnew.gas_1(:,2:4)./pnew.gas_1(:,1),1:101,'o-'); set(gca,'ydir','reverse'); ylim([85 100]); grid on

rtpwrite('ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp',hnew,ha,pnew,pa);

%% to do SARTA jacs, wit output radiance usunits = rad use listp=1 listj=-1 jacunit=0   (proflie 1, all jacs, output units)
%% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp listp=1 listj=-1 jacunit=0
%% 
%% then just run sarta without jacs to make rtp file for all 7 profiles (so you can check finite difference jacs)
%% /home/sergio/git/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/JACvers/bin/jac_airs_l1c_2834_cloudy_jan25_H2020 fin=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp fout=ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp

error('osjsjfs')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% using /bin/rm junk.*; ../BIN/kcarta.x90_v1.22     quickuseF90_v1.22.nml junk.dat junk.jac 

if ~exist('pnew')
  [hnew,ha,pnew,pa] = rtpread('ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp');
end

%% raw calcs
[r22,w]= readkcstd('junk.dat');
[j22,w]= readkcjac('junk.jac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% going into kCARTA and recompiling each time with iPertT = +/- 1  iPertQ = +/-1

%% %% and then eventually
%% %% go into kcarta and perturn T using iJax = 5
%% [r22_T,w]= readkcstd('junk.dat');
%% dT = 0.01;
%% %% go into kcarta and perturb Q using iJax = 5
%% [r22_Q,w]= readkcstd('junk.dat');
%% dq5 = 7.8986e-08; %% or whatever
%% 
%% %% and then eventually
%% %% go into kcarta and perturn T using iJax = 10
%% [r22_T10,w]= readkcstd('junk.dat');
%% dT = 0.01;
%% %% go into kcarta and perturb Q using iJax = 10
%% [r22_Q1-,w]= readkcstd('junk.dat');
%% dq10 = 4.9790e-08; %% or whatever
%% 
%% figure(1); plot(w,j22(:,2),'b.-',w,(r22_Q-r22)/dq5,'r'); xlim([970 980]); title('v1.22 Q')
%% figure(2); plot(w,j22(:,2+97),'b.-',w,(r22_T-r22)/dT,'r'); xlim([970 980]); title('v1.22 T')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% using /bin/rm junk.*; ../BIN/kcarta.x90_v1.22     quickuseF90_v1.22.nml junk.dat junk.jac -XYZ   where XYZ = 1 .. 7 for the above 7 profiles you will first get

dT = 0.1;

%%% also directly from the above 7 profile RTP file, iRTPCommandLine = 2,5 for Q,T lay5
[r22x_Q5,w]= readkcstd('junk.dat');  %% -XYZ = 2
[r22x_T5,w]= readkcstd('junk.dat');  %% -XYZ = 5
dq5 = 7.8986e-08;
dq5x = dq5;

%%% also directly from the above 7 profile RTP file, iRTPCommandLine = 3,6 for Q,T lay 10
[r22x_Q10,w]= readkcstd('junk.dat');  %% -XYZ = 3
[r22x_T10,w]= readkcstd('junk.dat');  %% -XYZ = 6
dq10 = 4.9790e-08;
dq10x = dq10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% using /bin/rm junk.*; ../BIN/kcarta.x_v1.18     quickuseF77_v1.18.nml junk.dat junk.jac 
[r18,w]= readkcstd('junk.dat');
[j18,w]= readkcjac('junk.jac');

%% using /bin/rm junk.*; ../BIN/kcarta.x_v1.10     quickuseF77_v1.10.nml junk.dat junk.jac 
[r10,w]= readkcstd('junk.dat');
[j10,w]= readkcjac('junk.jac');

%% using /bin/rm junk.*; ../BIN/kcarta.x_v1.07     quickuseF77_v1.07.nml junk.dat junk.jac 
[r07,w]= readkcstd('junk.dat');
[j07,w]= readkcjac('junk.jac');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(w,rad2bt(w,r07),'k',w,rad2bt(w,r10),'g',w,rad2bt(w,r18),'b',w,rad2bt(w,r22),'r'); legend('1.07','1.10','1.18','1.22','location','best')
plot(w,rad2bt(w,r22) - rad2bt(w,r07),'k',w,rad2bt(w,r22) - rad2bt(w,r10),'g',w,rad2bt(w,r22) - rad2bt(w,r18),'b');
  legend('1.07','1.10','1.18','location','best')

addpath /home/sergio/git/matlabcode
[fc,qc] = quickconvolve(w,[r07 r10 r18 r22],1,1); tc = rad2bt(fc,qc);

plot(fc,tc(:,4)-tc)
  legend('1.07','1.10','1.18','location','best')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = 2;

figure(1); plot(w,j07(:,ind+5),'k.-',w,j10(:,ind+5),'g.-',w,j18(:,ind+5),'b.-',w,j22(:,ind+5),'r.-',w,(r22x_Q10-r22)/dq10x,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('Q10')
figure(1); plot(w,j07(:,97+ind+5),'k.-',w,j10(:,97+ind+5),'g.-',w,j18(:,97+ind+5),'b.-',w,j22(:,97+ind+5),'r.-',w,(r22x_T10-r22)/dT,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('T10')
  
figure(1); plot(w,j07(:,ind+0),'k.-',w,j10(:,ind+0),'g.-',w,j18(:,ind+0),'b.-',w,j22(:,ind+0),'r.-',w,(r22x_Q5-r22)/dq5x,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('Q05')
figure(1); plot(w,j07(:,97+ind+0),'k.-',w,j10(:,97+ind+0),'g.-',w,j18(:,97+ind+0),'b.-',w,j22(:,97+ind+0),'r.-',w,(r22x_T5-r22)/dT,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('T05')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fc,jc07] = quickconvolve(w,j07,1,1);
[fc,jc10] = quickconvolve(w,j10,1,1);
[fc,jc18] = quickconvolve(w,j18,1,1);
[fc,jc22] = quickconvolve(w,j22,1,1);

[fc,deltaR_Q10] = quickconvolve(w,r22x_Q10-r22,1,1);
[fc,deltaR_T10] = quickconvolve(w,r22x_T10-r22,1,1);
[fc,deltaR_Q5]  = quickconvolve(w,r22x_Q5-r22,1,1);
[fc,deltaR_T5]  = quickconvolve(w,r22x_T5-r22,1,1);

ind = 2;

figure(1); plot(fc,jc07(:,ind+5),'k.-',fc,jc10(:,ind+5),'g.-',fc,jc18(:,ind+5),'b.-',fc,jc22(:,ind+5),'r.-',fc,deltaR_Q10/dq10x,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('Q10')
figure(1); plot(fc,jc07(:,97+ind+5),'k.-',fc,jc10(:,97+ind+5),'g.-',fc,jc18(:,97+ind+5),'b.-',fc,jc22(:,97+ind+5),'r.-',fc,deltaR_T10/dT,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('T10')

figure(1); plot(fc,jc07(:,ind+0),'k.-',fc,jc10(:,ind+0),'g.-',fc,jc18(:,ind+0),'b.-',fc,jc22(:,ind+0),'r.-',fc,deltaR_Q5/dq5x,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('Q5')
figure(1); plot(fc,jc07(:,97+ind+0),'k.-',fc,jc10(:,97+ind+0),'g.-',fc,jc18(:,97+ind+0),'b.-',fc,jc22(:,97+ind+0),'r.-',fc,deltaR_T5/dT,'c');
  legend('Jac 1.07','Jac 1.10','Jac 1.18','Jac 1.22','kcarta finite diff','location','best');
  title('T5')


figure(1); plot(fc,jc22(:,ind+5),'r.-',fc,deltaR_Q10/dq10x,'m',fc,jc22(:,ind+0),'b.-',fc,deltaR_Q5/dq5x,'c');
  legend('J22 lay10','finite lay10','J22 lay5','finite lay5','location','best'); title('Q jac')
figure(1); plot(fc,jc22(:,97+ind+5),'r.-',fc,deltaR_T10/dT,'m',fc,jc22(:,97+ind+0),'b.-',fc,deltaR_T5/dT,'c');
  legend('J22 lay10','finite lay10','J22 lay5','finite lay5','location','best'); title('T jac')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% from WORK directory 
%%%% /bin/cp -a quickuse*.nml driver_test_the_jacs.m test_jac.sc text_ecmwf_airicrad_day092_clear_unitemiss_4758.txt ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp ../TEST_JACS/
%%%% /bin/cp -a quickuse*.nml driver_test_the_jacs.m test_jac.sc text_ecmwf_airicrad_day092_clear_unitemiss_4758.txt ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp ../TEST_JACS/
%%%% /bin/cp -a quickuse*.nml driver_test_the_jacs.m test_jac.sc text_ecmwf_airicrad_day092_clear_unitemiss_4758.txt ecmwf_airicrad_day092_clear_unitemiss_profile_4758.op.rtp ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp ../TEST_JACS/

%% when you have SARTA jacs (2645 chans) and kcarta run from 605-2830 cm-1
%[fairs,tjacairs] = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacTZ',100);              %% 1-97 = atm, 98 = stemp
%plot(radsOut.freqAllChunks,sum(jacsOut.tjacAllChunks,2),'b',w22,sum(j10(:,indT),2),'g',w22,sum(j22(:,indT),2),'r',fairs,sum(tjacairs(:,1:97),2),'kx-');

[fsarta,tjacsarta]  = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacTZ',100); tjacsarta  = squeeze(tjacsarta);  tjacsarta = tjacsarta(:,1:97);   tjacsarta = fliplr(tjacsarta); %% 98 is surface temp
[fsarta,q1jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG1',1);   q1jacsarta = squeeze(q1jacsarta); q1jacsarta = q1jacsarta(:,1:97); q1jacsarta = fliplr(q1jacsarta);
[fsarta,q2jacsarta] = readsarta_jac('/home/sergio/git/kcarta_gen/WORK/ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp_jacG2',2);   q2jacsarta = squeeze(q2jacsarta); q2jacsarta = q2jacsarta(:,1:97); q2jacsarta = fliplr(q2jacsarta); 
[hy,hay,py,pay] = rtpread('ecmwf_airicrad_day092_clear_unitemiss_profile_4758.sarta.rtp');

%% run kcarta 1.18, only G1, G2
[r18,w18] = readkcstd('junk.dat');
[j18,w18] = readkcjac('junk.jac');
xj18_G1 = j18(:,1:97); j18_G2 = j18(:,(1:97)+1*97); j18_T = j18(:,(1:97)+2*97);
[fc,tc18]  = quickconvolve(w18,j18_T,1,1);
[fc,xq1c18] = quickconvolve(w18,xj18_G1,1,1);
[fc,q2c18] = quickconvolve(w18,j18_G2,1,1);
[fc,rc18]  = quickconvolve(w18,r18,1,1);
%% run kcarta 1.18, with G1,G101,G102,G103
[r18,w18] = readkcstd('junk.dat');
[j18,w18] = readkcjac('junk.jac');
j18_G1 = j18(:,(1:97)+0*97) + j18(:,(1:97)+1*97) + j18(:,(1:97)+2*97) + j18(:,(1:97)+3*97);
[fc,q1c18] = quickconvolve(w18,j18_G1,1,1);

%% run kcarta 1.22, only G1, G2
[r22,w22] = readkcstd('junk.dat');
[j22,w22] = readkcjac('junk.jac');
xj22_G1 = j22(:,1:97); j22_G2 = j22(:,(1:97)+1*97); j22_T = j22(:,(1:97)+2*97);
[fc,tc22]  = quickconvolve(w22,j22_T,1,1);
[fc,xq1c22] = quickconvolve(w22,xj22_G1,1,1);
[fc,q2c22] = quickconvolve(w22,j22_G2,1,1);
[fc,rc22]  = quickconvolve(w22,r22,1,1);
%% run kcarta 1.22, with G1,G101,G102,G103
[r22,w22] = readkcstd('junk.dat');
[j22,w22] = readkcjac('junk.jac');
j22_G1 = j22(:,(1:97)+0*97) + j22(:,(1:97)+1*97) + j22(:,(1:97)+2*97) + j22(:,(1:97)+3*97);
[fc,q1c22] = quickconvolve(w22,j22_G1,1,1);

plot_the_results_clear
