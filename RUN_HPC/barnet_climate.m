clear all

addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools

[h0,ha,p0,pa] = rtpread('RTP/stdNH3_1100mb_op_400ppm.rtp');
h0.nchan = 2378;
h0.ichan = (1:2378)';
[hnew,pnew] = replicate_rtp_headprof(h0,p0,49,9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iDo = -1;
if iDo > 0
[ppmvCO2,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(hnew,pnew,1,2);
[ppmvCH4,ppmvAVG,ppmvMAX,pavgLAY,tavgLAY,ppmv500,ppmv75] = layers2ppmv(hnew,pnew,1,6);

%% 40, 80, 160, 280, 370, 400, 420, 560, and 1120 
ii = 1; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 040/400;
ii = 2; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 080/400;
ii = 3; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 160/400;
ii = 4; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 280/400;
ii = 5; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 370/400;
ii = 6; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 400/400;
ii = 7; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 420/400;
ii = 8; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 560/400;
ii = 9; pnew.gas_2(:,ii) = pnew.gas_2(:,ii) * 1120/400; 

rtpwrite('RTP/barnet_climate.op.rtp',hnew,ha,pnew,pa);

disp('sbatch --array=1-9 sergio_matlab_jobB.sbatch')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/KCARTA/MATLAB

iDo = +1;
if iDo > 0

clear all
%  iasi-like, meaning 0.5 cm-1 resolution 0.25 cm-1 step gaussian.

[hnew,ha,pnew,pa] = rtpread('RTP/barnet_climate.op.rtp');

[d1,w] = readkcstd('JUNK/rad.dat1');
[d2,w] = readkcstd('JUNK/rad.dat2');
[d3,w] = readkcstd('JUNK/rad.dat3');
[d4,w] = readkcstd('JUNK/rad.dat4');
[d5,w] = readkcstd('JUNK/rad.dat5');
[d6,w] = readkcstd('JUNK/rad.dat6');
[d7,w] = readkcstd('JUNK/rad.dat7');
[d8,w] = readkcstd('JUNK/rad.dat8');
[d9,w] = readkcstd('JUNK/rad.dat9');

xppm = [040 080 160 280 370 400 420 560 1120];

t = rad2bt(w,[d1 d2 d3 d4 d5 d6 d7 d8 d9]);
plot(w,t-t(:,6)*ones(1,9))
[fc,qc] = quickconvolve(w,[d1 d2 d3 d4 d5 d6 d7 d8 d9],0.5,0.25); tc = rad2bt(fc,qc);
figure(1)
plot(fc,tc)
hl = legend('040','080','160','280','370','400','420','560','1120'); set(hl,'fontsize',10)

figure(2)
plot(fc,tc-tc(:,6)*ones(1,9),'linewidth',2)
hl = legend('040','080','160','280','370','400','420','560','1120'); set(hl,'fontsize',10)
axis([605 2830 -20 +10]); ylabel('\Delta BT'); xlabel('wavnumber cm-1')
grid
title('x ppm - 400 ppm')

figure(3)
f0 = 760;
f0 = 735;
f0 = 750;
aha = find(fc >= f0,1);
[Y,I] = sort(xppm);
plot(xppm(I),pnew.stemp(1)-tc(aha,I),'o-'); grid; xlabel('CO2 ppm'); ylabel(['stemp - BT at ' num2str(f0) ' cm-1']);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iDo = -1;
if iDo > 0

fid = fopen('barnet.txt','w');
junk = [0 xppm];
fprintf(fid,'%10.5f %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e \n',junk);
junk = [fc' qc];
fprintf(fid,'%10.5f %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e %10.7e \n',junk');
fclose(fid);

wah = load('barnet.txt');
wahCO2 = wah(1,:);
len = length(wah);
wah = wah(2:len,:);
wahf = wah(:,1);
wahr = wah(:,2:10);
waht = rad2bt(wahf,wahr);
f0 = 750;
aha = find(wahf >= f0,1);
plot(xppm,pnew.stemp(1)-waht(aha,:),'o-');
grid; xlabel('CO2 ppm'); ylabel(['stemp - BT at ' num2str(f0) ' cm-1']);

end