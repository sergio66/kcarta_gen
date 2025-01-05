addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

clear all

%==> /asl/rta/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcononly.txt <==
% AIRS coefficient tuning file, Scott Hannon <hannon@umbc.edu>
% AIRS RTA "m140" April 2008
% Created: 16 May 2008, Scott Hannon - water continuum tuning from
%    interpolation of "dec05deliv"; all other tunings are ones
%
% Note: channels with ID > 2378 are fake channels
%  ID   freq      fixed   w lines   w con    ozone     CO       CH4    nonLTE
%----  --------  -------  -------  -------  -------  -------  -------  -------
tuning = load('/asl/rta/sarta_database/Data_AIRS_apr08/Coef/tunmlt_wcon_nte.txt');
tuning = load('/asl/rta/sarta_database/Data_AIRS_apr08/Coef//tunmlt_17jun2011_ch4.txt');
[Y,I] = sort(tuning(:,2)); plot(tuning(I,2),tuning(I,4));  xlim([1200 1700])

%% see /asl/rta/sarta_database/Data_AIRS_may19/Coef/refprof_trace400

gas6 = load('/home/sergio/SPECTRA/IPFILES/std_gx6x_6');  %% need molecules/cm^2
gas6(:,5) = gas6(:,5)*6.023e26;

fin = '/home/chepplew/data/sarta_tuning/sorted/from_scott/armtwp_nearby_semiclear_jan04_resetT_nte_co2.rtp';
[h0,ha,p0,pa] = rtpread(fin);
vchan = h0.vchan;
h0 = rmfield(h0,'vchan');
  [ppmLAY,ppmAVG] = layers2ppmv(h0,p0,1:180,2);
  plot(p0.co2ppm,ppmAVG,'.')
h0.ngas = 5; %%% YEAH SOMETHING SCREWY HERE ..... h0.glist(6) = 6; h0.gunit(6) = 1; p0.gas_6 = flipud(gas6(:,5))*ones(1,180);

iTry = +1;
if iTry == -1
  %% try 1
  fout0 = 'raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';                                                                                                            rtpwrite(fout0,h0,ha,p0,pa);
  fout1 = 'set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';        h=h0;p=p0; h.gunit(2) = 10; p=rmfield(p,'gas_2');                                       rtpwrite(fout1,h,ha,p,pa);
  fout2 = 'set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';  h=h0;p=p0; h.gunit(2) = 10; p=rmfield(p,'gas_2'); p.co2ppm = ones(size(p.co2ppm))*400;  rtpwrite(fout2,h,ha,p,pa);
  fout3 = 'set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';          h=h0;p=p0; p.gas_2 = p.gas_2 * 430/370;                                                 rtpwrite(fout3,h,ha,p,pa);
elseif iTry == +1
  %% try 2
  fout0 = 'raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';                                                                                                                         rtpwrite(fout0,h0,ha,p0,pa);
  fout1 = 'set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';        h=h0;p=p0; h.ngas=4; h.glist=h.glist([1 3 4 5]); h.gunit=h.gunit([1 3 4 5]); p=rmfield(p,'gas_2');   rtpwrite(fout1,h,ha,p,pa);
  fout2 = 'set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';  h=h0;p=p0; h.ngas=4; h.glist=h.glist([1 3 4 5]); h.gunit=h.gunit([1 3 4 5]); p=rmfield(p,'gas_2'); p.co2ppm = ones(size(p.co2ppm))*400;  
                                                                                            rtpwrite(fout2,h,ha,p,pa);
  fout3 = 'set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp';          h=h0;p=p0; p.gas_2 = p.gas_2 * 430/370;                                                 rtpwrite(fout3,h,ha,p,pa);
end

sarta = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3';         %% see eg /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/clustbatch_redo_stemp_wv_cloud_filelist.m
%sarta = '/asl/packages/sartaV108/BinV201/sarta_airs_PGEv6_preNov2003';
%sarta = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';

disp(' ')
fprintf(1,'sarta = %s \n',sarta)
disp(' ')

rmer = ['!/bin/rm *.rp.rtp']; eval(rmer);
sartaer = ['!time ' sarta ' fin=raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp                   fout=raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp'];                   eval(sartaer);
sartaer = ['!time ' sarta ' fin=set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp       fout=set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp'];       eval(sartaer);
sartaer = ['!time ' sarta ' fin=set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp fout=set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp']; eval(sartaer);
sartaer = ['!time ' sarta ' fin=set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.op.rtp         fout=set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp'];         eval(sartaer);

[hf0,ha,pf0,pa] = rtpread('raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp');
[hf0,ha,pf1,pa] = rtpread('set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp');
[hf0,ha,pf2,pa] = rtpread('set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp');
[hf0,ha,pf3,pa] = rtpread('set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp');

addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/CLOUD
g = dogoodchan;
h.vchan = vchan;
tobs = rad2bt(vchan,p0.robs1);
t0 = rad2bt(vchan,pf0.rcalc);
t1 = rad2bt(vchan,pf1.rcalc);
t2 = rad2bt(vchan,pf2.rcalc);
t3 = rad2bt(vchan,pf3.rcalc);

figure(1)
plot(vchan(g),nanmean(tobs(g,:)'-t0(g,:)'),'kx-',vchan(g),nanmean(tobs(g,:)'-t1(g,:)'),'b.-',vchan(g),nanmean(tobs(g,:)'-t2(g,:)'),'g',vchan(g),nanmean(tobs(g,:)'-t3(g,:)'),'r'); xlim([640 1640]);
  hl = legend('raw gas2=369ppm','gunit=10,co2ppm=370','gunit=10,co2ppm=400','gas2=430','location','best','fontsize',8); grid
plot(vchan(g),tobs(g,1)-t0(g,1),'kx-',vchan(g),tobs(g,1)-t1(g,1),'b.-',vchan(g),tobs(g,1)-t2(g,1),'g',vchan(g),tobs(g,1)-t3(g,1),'r'); xlim([640 1640]);
  hl = legend('raw gas2=369ppm','gunit=10,co2ppm=370','gunit=10,co2ppm=400','gas2=430','location','best','fontsize',8); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kcarta = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res';

do_kcartaer

addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
[r0,w] = readkcstd('rad0.dat');
[r1,w] = readkcstd('rad1.dat');
[r2,w] = readkcstd('rad2.dat');
[r3,w] = readkcstd('rad3.dat');
[fc,qc] = convolve_airs(w,[r0 r1 r2 r3],1:2378); tc = rad2bt(fc,qc);

figure(2)
plot(vchan(g),nanmean(tobs(g,:)'-tc(g,1)'),'kx-',vchan(g),nanmean(tobs(g,:)'-tc(g,2)'),'b.-',vchan(g),nanmean(tobs(g,:)'-tc(g,3)'),'g',vchan(g),nanmean(tobs(g,:)'-tc(g,4)'),'r'); xlim([640 1640]);
  hl = legend('raw gas2=369ppm','gunit=10,co2ppm=370','gunit=10,co2ppm=400','gas2=430','location','best','fontsize',8); grid
plot(vchan(g),tobs(g,1)-tc(g,1),'kx-',vchan(g),tobs(g,1)-tc(g,2),'b.-',vchan(g),tobs(g,1)-tc(g,3),'g',vchan(g),tobs(g,1)-tc(g,4),'r'); xlim([640 1640]);
  hl = legend('raw gas2=369ppm','gunit=10,co2ppm=370','gunit=10,co2ppm=400','gas2=430','location','best','fontsize',8); grid

figure(1); axis([640 1640 -4 +4])
figure(2); axis([640 1640 -4 +4])

figure(3)
plot(vchan(g),t0(g,1)-tc(g,1),'kx-',vchan(g),t1(g,1)-tc(g,2),'b.-',vchan(g),t2(g,1)-tc(g,3),'g',vchan(g),t3(g,1)-tc(g,4),'r'); xlim([640 1640]); title('SARTA-kCARTA')
  hl = legend('raw gas2=369ppm','gunit=10,co2ppm=370','gunit=10,co2ppm=400','gas2=430','location','best','fontsize',8); grid
figure(3); axis([640 1640 -0.5 +0.5])
