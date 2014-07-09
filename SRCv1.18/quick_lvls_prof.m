figure(1); clf
figure(2); clf
figure(3); clf
figure(4); clf
figure(5); clf

xstartup
format short e
iProf = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is from /home/sergio/klayersV205/Data/glatm.dat

USStd_P = [...
 1.013E+03, 8.988E+02, 7.950E+02, 7.012E+02, 6.166E+02,   ...
 5.405E+02, 4.722E+02, 4.111E+02, 3.565E+02, 3.080E+02,   ...       
 2.650E+02, 2.270E+02, 1.940E+02, 1.658E+02, 1.417E+02,   ...       
 1.211E+02, 1.035E+02, 8.850E+01, 7.565E+01, 6.467E+01,   ...       
 5.529E+01, 4.729E+01, 4.047E+01, 3.467E+01, 2.972E+01,   ...       
 2.549E+01, 1.743E+01, 1.197E+01, 8.010E+00, 5.746E+00,   ...       
 4.150E+00, 2.871E+00, 2.060E+00, 1.491E+00, 1.090E+00,   ...       
 7.978E-01, 4.250E-01, 2.190E-01, 1.090E-01, 5.220E-02,   ...       
 2.400E-02, 1.050E-02, 4.460E-03, 1.840E-03, 7.600E-04,   ...       
 3.200E-04, 1.450E-04, 7.100E-05, 4.010E-05, 2.540E-05]           ;
                                                                
USStd_T = [...
    288.20,    281.70,    275.20,    268.70,    262.20,   ...
    255.70,    249.20,    242.70,    236.20,    229.70,   ...       
    223.30,    216.80,    216.70,    216.70,    216.70,   ...       
    216.70,    216.70,    216.70,    216.70,    216.70,   ...       
    216.70,    217.60,    218.60,    219.60,    220.60,   ...       
    221.60,    224.00,    226.50,    230.00,    236.50,   ...       
    242.90,    250.40,    257.30,    264.20,    270.60,   ...       
    270.70,    260.80,    247.00,    233.30,    219.60,   ...       
    208.40,    198.60,    188.90,    186.90,    188.40,   ...       
    195.10,    208.80,    240.00,    300.00,    360.00];           



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

airslevels = load('/home/sergio/MATLABCODE/airslevels.dat');

[h,ha,p,pa] = rtpread('/asl/s1/sergio/junk.ip.rtp');

data = [p.plevs(1:p.nlevs(iProf),iProf)      p.ptemp(1:p.nlevs(iProf),iProf)];
data = [data p.gas_1(1:p.nlevs(iProf),iProf) p.gas_2(1:p.nlevs(iProf),iProf) p.gas_3(1:p.nlevs(iProf),iProf)];
data = [data p.gas_5(1:p.nlevs(iProf),iProf) p.gas_6(1:p.nlevs(iProf),iProf) p.gas_9(1:p.nlevs(iProf),iProf)];

boo = find(data(:,1) > 1.0);
boo = find(data(:,1) > 5e-3);
data = data(boo,:);

fid = fopen('levels_prof1.txt','w');
str = 'test : prof1 from Regr49';
fprintf(fid,'''%s ''\n',str);
fprintf(fid,'%3i \n',length(boo));
fprintf(fid,'%8.4f %8.4f %8.4f \n',p.spres(1),p.stemp(1),0.0);
fprintf(fid,'%8.4f %8.4f %8.4f      \n',2008.0,p.plat(1),0.0);
fprintf(fid,'%3i \n',h.ngas);
fprintf(fid,'%3i %3i %3i %3i %3i %3i \n',h.glist);
fprintf(fid,'%3i %3i %3i %3i %3i %3i \n',h.gunit);

fprintf(fid,'%8.6e %8.5f %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e \n',data');

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h2,ha2,p2,pa2] = rtpread('/asl/s1/sergio/junk.op.rtp');
pavgN  = p2.plevs(1:100)-p2.plevs(2:101);;
pavgD  = log(p2.plevs(1:100)./p2.plevs(2:101));
pavg = pavgN./pavgD;
ix = 001:100; 
dataOut = [pavg' p2.ptemp(ix,iProf) p2.gas_1(ix,iProf) p2.gas_2(ix,iProf) p2.gas_3(ix,iProf) p2.gas_4(ix,iProf) p2.gas_5(ix,iProf)];

%% now load in outputt from kcarta. and look at 
%      print *,iL,PAVG_KCARTADATABASE_AIRS(iL),raTout(iL),raAmountOut(iL),raaQout(iL,1),raaQout(iL,2),
%     $        raaG_MRX(iL-iLowestLev+1,1),raaG_MRX(iL-iLowestLev+1,2)
ix = 90:101; ix = 90 : 97;flipud([p2.ptemp(ix,iProf) p2.gas_1(ix,iProf) p2.gas_2(ix,iProf) p2.gas_3(ix,iProf)])

%% this confirms that I have interpolated temps fron USER to KCARTA levels quite well!!!!
kc2 = load('ugh2');
semilogy(data(:,2),data(:,1),'bo-',kc2(:,3),kc2(:,2)/100,'rx-','linewidth',2); grid
  set(gca,'ydir','reverse')
ax = axis;
for ii = 3 : 4 : length(airslevels)
  line([170 300],[airslevels(ii) airslevels(ii)],'color','r','linestyle','--');
end
axis([170 300 0.001 1020]);

%kcx = load('ugh');
kcx = load('/home/sergio/SPECTRA/IPFILES/std_co2');

kc = load('ugh1');

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
%semilogy(data(:,2),data(:,1),dataOut(:,2),dataOut(:,1),'r',kc(:,3),kc(:,2),'k',kcx(:,4),kcx(:,6),'g'); grid
semilogy(data(:,2),data(:,1),'bo-',dataOut(:,2),dataOut(:,1),'r',kc(:,3),kc(:,2),'k',kcx(:,4),kcx(:,2)*1013.25,'g',...
         USStd_T,USStd_P,'c'); grid
%semilogy(data(:,2),data(:,1),dataOut(:,2),dataOut(:,1),'r.-',kc(:,3),kc(:,2),'ko-'); grid
  set(gca,'ydir','reverse')
axis([170 300 0.001 1020]);

semilogy(flipud(dataOut(1:97,2)) - kc(:,3),kc(:,2),'o-'); grid
  set(gca,'ydir','reverse'); 
axis([-0.5 +0.5 0 1000]);
title('klayer - kcarta layers')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
loglog(dataOut(:,3),dataOut(:,1),'ro-',kc(:,5),kc(:,2),'k+-'); grid
  set(gca,'ydir','reverse'); title('WV')
  axis([1e16 1e21 0 1000]);
pause(1)

loglog(dataOut(:,4),dataOut(:,1),'ro-',kc(:,6),kc(:,2),'k+-'); grid
  set(gca,'ydir','reverse'); title('CO2')
  axis([1e16 1e21 0 1000]);
pause(1)

loglog(dataOut(:,5),dataOut(:,1),'ro-',kc(:,7),kc(:,2),'k+-'); grid
  set(gca,'ydir','reverse'); title('O3')
  axis([1e13 1e18 0 1000]);
pause(1)

loglog(dataOut(:,6),dataOut(:,1),'ro-',kc(:,8),kc(:,2),'k+-'); grid
  set(gca,'ydir','reverse'); title('N2O')
  axis([1e11 1e18 0 1000]);
pause(1)

loglog(dataOut(:,7),dataOut(:,1),'ro-',kc(:,9),kc(:,2),'k+-'); grid
  set(gca,'ydir','reverse'); title('CO')
  axis([1e14 1e17 0 1000]);
pause(1)

semilogy(flipud(dataOut(1:97,3:7))./kc(:,5:9),kc(:,2)); grid
semilogy(flipud(dataOut(1:97,3:5))./kc(:,5:7),kc(:,2),'linewidth',2); grid
  set(gca,'ydir','reverse'); 
axis([0.95 1.15 0 1000])
title('ratio klayers/ kcarta layers')

figure(3)
aN = flipud(dataOut(1:97,3:7))-kc(:,5:9);
aD = flipud(dataOut(1:97,3:7));
aN = flipud(dataOut(1:97,3:5))-kc(:,5:7);
aD = flipud(dataOut(1:97,3:5));
semilogy(aN./aD * 100,kc(:,2),'linewidth',2); grid
  set(gca,'ydir','reverse'); 
axis([-5 +5 0  1000])
title('percent diff (klayers - kcarta layers)/klayers')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath ../MATLAB
[d0,w0] = readkcstd('../WORK/junk0.dat');
[dL50,wL50] = readkcstd('../WORK/junkLVL.dat');    %% iNFine = 50
[dL10,wL10] = readkcstd('../WORK/junkLVL10.dat');    %% iNFine = 10
[dL200,wL200] = readkcstd('../WORK/junkLVL200.dat');    %% iNFine = 200

[fc,qc] = quickconvolve(w0,[d0 dL50 dL10],0.5,0.5);
[fc,qc] = quickconvolve(w0,[d0 dL50 dL200],0.5,0.5);

figure(4)
plot(w0,rad2bt(w0,d0)-rad2bt(w0,dL50),w0,rad2bt(w0,d0)-rad2bt(w0,dL10),'r')

figure(5) 
plot(fc,rad2bt(fc,qc(:,1))-rad2bt(fc,qc(:,2)),'b',...
     fc,rad2bt(fc,qc(:,1))-rad2bt(fc,qc(:,3)),'r',...
     fc,rad2bt(fc,qc(:,2))-rad2bt(fc,qc(:,3)),'k')
