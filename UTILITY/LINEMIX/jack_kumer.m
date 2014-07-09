cd /carrot/s1/sergio/AIRSCO2/KCARTA_TRIES

[all_lte,wlte] = readkcstd('all_lte.dat');    %%lte computation, all gases
[co2_lte,wlte] = readkcstd('co2_lte.dat');    %%lte computation, only CO2

%%% 10 secs
%[d2205_all,w2205] = readkcstd('all2205.dat');      all gases used, CO2 in NLTE
%                                                   cousin lineshape
%[d2205_lte,w2205] = readkcstd('cousin_2205.dat');  only CO2, all bands in LTE
%                                                   cousin lineshape
%[d2205_cou,w2205] = readkcstd('cousin_2205.dat');  only CO2, some band in NLTE
%                                                   cousin lineshape
%[d2205_cou,w2205] = readkcstd('newcous_2205.dat'); same as above; older vers
%[d2205_new,w2205] = readkcstd('co2_2205new.dat');  linemixing, new vers
%[d2205,    w2205] = readkcstd('co2_2205.dat');     linemixing, old vers

%%% 10 secs
[d2205_all,w2205] = readkcstd('all2205.dat');
[d2205_lte,w2205] = readkcstd('cousin_2205.dat');
[d2205_cou,w2205] = readkcstd('cousin_2205.dat');
[d2205_cou,w2205] = readkcstd('newcous_2205.dat');
[d2205_new,w2205] = readkcstd('co2_2205new.dat');
[d2205,    w2205] = readkcstd('co2_2205.dat');

% 186 min
[d2230_all,w2230] = readkcstd('all2230.dat');
[d2230_lte,w2230] = readkcstd('couslte_2230.dat');
[d2230_cou,w2230] = readkcstd('cousin_2230.dat');
[d2230_cou,w2230] = readkcstd('newcous_2230.dat');
[d2230_new,w2230] = readkcstd('co2_2230new.dat');
[d2230,    w2230] = readkcstd('co2_2230.dat');

% 273 min
[d2255_all,w2255] = readkcstd('all2255.dat');
[d2255_lte,w2255] = readkcstd('couslte_2255.dat');
[d2255_cou,w2255] = readkcstd('cousin_2255.dat');
[d2255_cou,w2255] = readkcstd('newcous_2255.dat');
[d2255_new,w2255] = readkcstd('co2_2255new.dat');
[d2255,    w2255] = readkcstd('co2_2255.dat');

% 292 min
[d2280_all,w2280] = readkcstd('all2280.dat');
[d2280_lte,w2280] = readkcstd('couslte_2280.dat');
[d2280_cou,w2280] = readkcstd('cousin_2280.dat');
[d2280_cou,w2280] = readkcstd('newcous_2280.dat');
[d2280_new,w2280] = readkcstd('co2_2280new.dat');
[d2280,    w2280] = readkcstd('co2_2280.dat');

% 263 min
[d2305_all,w2305] = readkcstd('all2305.dat');
[d2305_lte,w2305] = readkcstd('couslte_2305.dat');
[d2305_cou,w2305] = readkcstd('cousin_2305.dat');
[d2305_cou,w2305] = readkcstd('newcous_2305.dat');
[d2305_new,w2305] = readkcstd('co2_2305new.dat');
[d2305,    w2305] = readkcstd('co2_2305.dat');

% 154 min
[d2330_all,w2330] = readkcstd('all2330.dat');
[d2330_lte,w2330] = readkcstd('couslte_2330.dat');
[d2330_cou,w2330] = readkcstd('cousin_2330.dat');
[d2330_cou,w2330] = readkcstd('newcous_2330.dat');
[d2330_new,w2330] = readkcstd('co2_2330new.dat');
[d2330,    w2330] = readkcstd('co2_2330.dat');

% 65 min
[d2355_all,w2355] = readkcstd('all2355.dat');
[d2355_lte,w2355] = readkcstd('couslte_2355.dat');
[d2355_cou,w2355] = readkcstd('cousin_2355.dat');
[d2355_cou,w2355] = readkcstd('newcous_2355.dat');
[d2355_new,w2355] = readkcstd('co2_2355new.dat');
[d2355,    w2355] = readkcstd('co2_2355.dat');

% 21 min
[d2380_all,w2380] = readkcstd('all2380.dat');
[d2380_lte,w2380] = readkcstd('couslte_2380.dat');
[d2380_cou,w2380] = readkcstd('cousin_2380.dat');
[d2380_cou,w2380] = readkcstd('newcous_2380.dat');
%[d2380_new,w2380] = readkcstd('co2_2380new.dat');
[d2380_new,w2380] = readkcstd('cousin_2380.dat');
[d2380,    w2380] = readkcstd('co2_2380.dat');

% 13 min
[d2405_all,w2405] = readkcstd('all2405.dat');
[d2405_lte,w2405] = readkcstd('couslte_2405.dat');
[d2405_cou,w2405] = readkcstd('cousin_2405.dat');
[d2405_cou,w2405] = readkcstd('newcous_2405.dat');
[d2405_new,w2405] = readkcstd('co2_2405new.dat');
[d2405,    w2405] = readkcstd('co2_2405.dat');

% 17 min
[d2430_all,w2430] = readkcstd('all2430.dat');
[d2430_lte,w2430] = readkcstd('couslte_2430.dat');
[d2430_cou,w2430] = readkcstd('cousin_2430.dat');
[d2430_cou,w2430] = readkcstd('newcous_2430.dat');
[d2430_new,w2430] = readkcstd('co2_2430new.dat');
[d2430,    w2430] = readkcstd('co2_2430.dat');

% 19 min
%[d2455_all,w2455] = readkcstd('all2455.dat');
%[d2455_lte,w2455] = readkcstd('couslte_2455.dat');
%[d2455_cou,w2455] = readkcstd('cousin_2455.dat');
%[d2455_cou,w2455] = readkcstd('newcous_2455.dat');
%[d2455_new,w2455] = readkcstd('co2_2455new.dat');
%[d2455,    w2455] = readkcstd('co2_2455.dat');

wall = [w2205;     w2230;     w2255;     w2280;     w2305;     w2330];
dall = [d2205_all; d2230_all; d2255_all; d2280_all; d2305_all; d2330_all];
dclte= [d2205_lte; d2230_lte; d2255_lte; d2280_lte; d2305_lte; d2330_lte];
dcou = [d2205_cou; d2230_cou; d2255_cou; d2280_cou; d2305_cou; d2330_cou];
dnew = [d2205_new; d2230_new; d2255_new; d2280_new; d2305_new; d2330_new];
dold = [d2205;     d2230;     d2255;     d2280;     d2305;     d2330];

wall = [wall;  w2355;     w2380;     w2405;     w2430];
dall = [dall;  d2355_all; d2380_all; d2405_all; d2430_all];
dclte= [dclte; d2355_lte; d2380_lte; d2405_lte; d2430_lte];
dcou = [dcou;  d2355_cou; d2380_cou; d2405_cou; d2430_cou];
dnew = [dnew;  d2355_new; d2380_new; d2405_new; d2430_new];
dold = [dold;  d2355;     d2380;     d2405;     d2430];

plot(wall,rad2bt(wall,dold),wall,rad2bt(wall,dnew),...
     wall,rad2bt(wall,dcou),wlte,rad2bt(wlte,co2_lte))
plot(wall,rad2bt(wall,dnew)-rad2bt(wall,co2_lte),...
     wall,rad2bt(wall,dold)-rad2bt(wall,co2_lte),...
     wall,rad2bt(wall,dcou)-rad2bt(wall,co2_lte))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error('oops')

fname = '/carrot/s1/sergio/AIRSCO2/NONLTE/GRANULE_14july2002/';
fname = [fname 'uniform_clear029.rtp'];
[h,ha,p,pa] = rtpread(fname);
fobs = h.vchan; 
tobs  = rad2bt(fobs,p.robs1); mo = mean(tobs');
tcalc = rad2bt(fobs,p.rcalc); mc = mean(tcalc');

edw_lte  = load('/carrot/s1/sergio/AIRSCO2/DEDWARDS/cut_lte.dat');
edw_nlte = load('/carrot/s1/sergio/AIRSCO2/DEDWARDS/cut_nlte.dat');

aa = 2;     %%%AIRS res FWHM
bb = 1;     %%%AIRS res spacing

ed_f = edw_lte(:,1);
wall = 2205:0.0025:2455-0.0025;
ed_lte  = spline(ed_f,edw_lte(:,2),wall);
ed_nlte = spline(ed_f,edw_nlte(:,2),wall);
[fce,qce] = quickconvolve(wall,[ed_lte ed_nlte]*1000,aa,bb);

%[fc,qc] =quickconvolve(wall,[co2_lte dclte dold dnew dcou],aa,bb);
[fc,qc] =quickconvolve(wall,[co2_lte all_lte dall dclte dold dnew dcou],aa,bb);
plot(fc,rad2bt(fc,qc(3,:)')-rad2bt(fc,qc(2,:)),...
     fce,rad2bt(fce,qce(2,:)')-rad2bt(fce,qce(1,:)'),fobs,mo - mc,'.-',...
     fce,rad2bt(fce,qce(1,:)')-rad2bt(fc,qc(2,:)'),'-.');
axis([2200 2450 -5 15]); grid; legend('kcarta','genln2','obs','G-K(LTE)',0);

plot(fce,rad2bt(fce,qce'),fc,rad2bt(fc,qc([2 3],:)'))
axis([2200 2450 200 300]); grid; legend('edw L','edw NL','kc L','kc NL',0); 

all_nlte = dall;           %%%% nonLTE at high resolution
kc_conv_lte  = qc(2,:)';   %%%% LTE convolved; FWHM = 0.3; spacing = 0.3/5
kc_conv_nlte = qc(3,:)';   %%%% NLTE convolved; FWHM = 0.3; spacing = 0.3/5
save /carrot/s1/sergio/kumer_high.mat wall all_lte all_nlte
save /carrot/s1/sergio/kumer_low.mat  fc   kc_conv_lte  kc_conv_nlte

aa = 0.3;        %%%high res FWHM
bb = aa/5;       %%%high res spacing
aa = 2;     %%%AIRS res FWHM
bb = 1;     %%%AIRS res spacing

%% kopra4 ==> start at 80, kopra3 ==> start at 90 ONLY SigSig in NLTE
%% kopra5 ==> start at 90, SigSig,PiPi,DeltDelt in NLTE
[kopraLTE,w]  = readkcstd('kopraLTE_3.dat');
[kopraLTEa,wa]  = readkcstd('kopraLTE_3a.dat');

[kopraNLTE,w] = readkcstd('kopraNLTE_5.dat');
[kopraNLTEa,wa] = readkcstd('kopraNLTE_5a.dat');

%[kopraNLTEb,wb] = readkcstd('newtryb.dat');
%[kopraNLTE,w] = readkcstd('newtry.dat');
%[kopraNLTEa,wa] = readkcstd('newtrya.dat');
%w  = [w; wb];
%kopraNLTE = [kopraNLTE; kopraNLTEb];

[kopraNLTEa,wa] = readkcstd('newtry2_2205_2255.dat');
[kopraNLTEb,wb] = readkcstd('newtry2_2255_2280.dat');
[kopraNLTEc,wc] = readkcstd('newtry2_2280_2305.dat');
[kopraNLTEd,wd] = readkcstd('newtry2_2305_2330.dat');
[kopraNLTEe,we] = readkcstd('newtry2_2330_2355.dat');
[kopraNLTEf,wf] = readkcstd('newtry2_2355_2455.dat');
w  = [wa; wb; wc; wd; we; wf];
kopraNLTE = [kopraNLTEa; kopraNLTEb; kopraNLTEc; ... 
             kopraNLTEd; kopraNLTEe; kopraNLTEf];

ww = [w]; 
%kopraNLTE_all = [kopraNLTEa; kopraNLTE]; 
kopraNLTE_all = [kopraNLTE]; 
kopraLTE_all = [kopraLTEa; kopraLTE];
[fce,qce] = quickconvolve(wall,[ed_lte; ed_nlte]*1000,aa,bb);
[fk,kopra] = quickconvolve(ww,[kopraLTE_all kopraNLTE_all],aa,bb);
plot(fce,rad2bt(fce,qce(2,:)')-rad2bt(fce,qce(1,:)'),fobs,mo - mc,'.-',...
      fk,rad2bt(fk,kopra(2,:)')-rad2bt(fk,kopra(1,:)'))
axis([2340 2455 -10 10])

fA = 2350; fB = 2420;
fA = 2240; fB = 2280;

ii = find(ww >= fA & ww <= fB);
jack_mono = [ww(ii)  kopraLTE_all(ii)  kopraNLTE_all(ii)]; whos jack_mono
fid=fopen('jack_monoA.txt','w');
fprintf(fid,'%12.10e %12.10e %12.10e \n',jack_mono');
fclose(fid);

aa = 0.3;        %%%high res FWHM
bb = aa/5;       %%%high res spacing
[fk,kopra] = quickconvolve(ww,[kopraLTE_all kopraNLTE_all],aa,bb); 
ii = find(fk >= fA & fk <= fB);
jack_high = [fk(ii);  kopra(1,ii);  kopra(2,ii)]'; whos jack_high
fid=fopen('jack_highA.txt','w');
fprintf(fid,'%12.10e %12.10e %12.10e \n',jack_high');
fclose(fid);

aa = 2;     %%%AIRS res FWHM
bb = 1;     %%%AIRS res spacing
[fk,kopra] = quickconvolve(ww,[kopraLTE_all kopraNLTE_all],aa,bb); 
ii = find(fk >= fA & fk <= fB);
jack_low = [fk(ii);  kopra(1,ii);  kopra(2,ii)]'; whos jack_low
fid=fopen('jack_lowA.txt','w');
fprintf(fid,'%12.10e %12.10e %12.10e \n',jack_low');
fclose(fid);

aa = 1;     %%% med res FWHM
bb = 1;     %%% AIRS res spacing
aa = 1;     %%% med res FWHM
bb = 0.5;   %%% J.Kumer requested this
[fk,kopra] = quickconvolve(ww,[kopraLTE_all kopraNLTE_all],aa,bb); 
ii = find(fk >= fA & fk <= fB);
jack_med = [fk(ii);  kopra(1,ii);  kopra(2,ii)]'; whos jack_med
fid=fopen('jack_medA.txt','w');
fprintf(fid,'%12.10e %12.10e %12.10e \n',jack_med');
fclose(fid);




