%% data from http://rredc.nrel.gov/solar/spectra/am0/

aVIS = load('E490_00a_AM0_short.csv');
aNIR = load('E490_00a_AM0_SHORT_IR.csv');
whos

a(:,1) = [aVIS(:,1); aNIR(:,1)];
a(:,2) = [aVIS(:,2); aNIR(:,2)];

xstartup
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solidangle = 6.785087652174316e-5;
sb = 5.6704e-8;

xx = 10 : 10 : 200000;
yy = ttorad(xx,5800);
%% integral over (0,2pi) dphi  (p/pi/2) cos (theta) d(cos theta) = pi
[trapz(xx,yy)/1000*pi sb*(5800)^4]    %% pi * int ttorad(x,T) = sig T^4
yyxx = ttorad(xx,5800) * solidangle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 500 : 10 : 100000; y = ttorad(x,5800);
[lout,rout] = rads_wnum2lamda(x,y);           
rout = rout * 1000 * solidangle;  %% change from W /m2/um/sr to mW /m2/um/sr and account for solidangle

%%% >>>>>> this gives 1360 W/m2 <<<<<<<<<<<
[trapz(x,y)/1000*solidangle trapz(lout,rout)/1000 trapz(a(:,1),a(:,2)*1000)/1000]
%%% >>>>>> this gives 1360 W/m2 <<<<<<<<<<<

figure(1);
plot(a(:,1),a(:,2)*1000,lout,rout)
axis([0 10 0 2500000])
xlabel('wavelength um'); ylabel('mW/m2/um/sr')

[boo,woo] = rads_lamda2wnum(a(:,1),a(:,2));
%%% >>>>>> this gives 1360 W/m2 <<<<<<<<<<<
[trapz(lout,rout)/1000 trapz(a(:,1),a(:,2)*1000)/1000 trapz(boo,woo)/1000]
%%% >>>>>> this gives 1360 W/m2 <<<<<<<<<<<

figure(2)
semilogy(boo,woo,x,y*solidangle)
xlabel('wavenumber cm-1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
this is from KCARTA/INCLUDE/pre_defined.param
      DATA kaMinFr        /   15.0000,   30.0000,   50.0000,   80.0000,
     $                    140.0000,  300.0000,  500.0000,  605.0000,
     $                    2830.0000,  3550.0000,  5550.0000,  8250.0000,
     $                    12000.0000,  25000.0000  /
      DATA kaMaxFr        /   30.0000,   50.0000,   80.0000,  150.0000,
     $                    310.0000,  510.0000,  605.0000,  2830.0000,
     $                    3580.0000,  5650.0000,  8400.0000,  12250.0000,
     $                    25000.0000,  44000.0000  /
      DATA kaFrStep       /  5.0000e-05,  1.0000e-04,  1.5000e-04,  2.5000e-04,
     $                    5.0000e-04,  1.0000e-03,  1.5000e-03,  2.5000e-03,
     $                    2.5000e-03,  1.0000e-02,  1.5000e-02,  2.5000e-02,
     $                    5.0000e-02,  1.0000e-01  /
%}

kaMinFr = [15.0000,   30.0000,   50.0000,   80.0000, ...
           140.0000,  300.0000,  500.0000,  605.0000, ...
           2830.0000,  3550.0000,  5550.0000,  8250.0000, ...
           12000.0000,  25000.0000];

kaMaxFr = [30.0000,   50.0000,   80.0000,  150.0000, ...
           310.0000,  510.0000,  605.0000,  2830.0000, ...
           3580.0000,  5650.0000,  8400.0000,  12250.0000, ...
           25000.0000,  44000.0000  ];

kaFrStep = [5.0000e-05,  1.0000e-04,  1.5000e-04,  2.5000e-04, ...
            5.0000e-04,  1.0000e-03,  1.5000e-03,  2.5000e-03, ...
            2.5000e-03,  1.0000e-02,  1.5000e-02,  2.5000e-02, ...
            5.0000e-02,  1.0000e-01  ];

figure(3); clf

iCnt = 0;
fkc = [];
rkc = [];
iNeed = 0;
for ii = 1 : length(kaMinFr)
  df = kaFrStep(ii);
  npts = (kaMaxFr(ii) - kaMinFr(ii))/df;
  nchunk = npts/10000;
  fprintf(1,'Start Freq    Df    NCHUNK = %8.6f  %8.6f  %3i \n',kaMinFr(ii),df*10000,nchunk)
  for jj = 1 : nchunk
    iCnt = iCnt + 1;
    ystart(iCnt) = kaMinFr(ii) + (jj-1)*df*10000;
    yend(iCnt)   = kaMinFr(ii) + (jj-0)*df*10000;

    if ystart(iCnt) >= 310
      fname = ['rad_solar' num2str(ystart(iCnt)) '.dat'];
      if exist(fname)
        fid = fopen(fname, 'r');
        if fid == -1
          error(sprintf('can not open %s', fname));
        end
    
        % header = fstart,fstop, df
        j = fread(fid, 1, 'integer*4');
        t = fread(fid, 3, 'real*8');
        j = fread(fid, 1, 'integer*4');
        [fstart, fend, df] = mdeal(t);

        freqs = ystart(iCnt) : df : yend(iCnt)-df; 
        freqs = fstart : df : fend; 
  
        j = fread(fid, 1, 'integer*4');
        temp = fread(fid, j/8, 'real*8');
        j = fread(fid, 1, 'integer*4');

        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
        fkc = [fkc freqs];
        rkc = [rkc; temp];
        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
        %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle

        %figure(3); plot(freqs,temp/(4*pi),'k'); hold on
        figure(3); plot(freqs,temp*1000*solidangle,'k.-'); hold on
      else
        fprintf(1,'hmmm %s does not exist??? \n',fname)
      end
    end
  end
end

figure(3);
plot(boo,woo,'r',x,y*solidangle,'g')
hold off
xlabel('wavenumber cm-1'); ylabel('TOA radiance mW/m2/sr/cm-1')

trapz(x,y)
trapz(x,y)/1000*pi
[trapz(x,y)/1000*pi trapz(xx,yy)/1000*pi sb*(5800)^4 trapz(a(:,1),a(:,2))]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% from sarta, see rdsun.f and incFTC.f
FNSUN ='/asl/data/sarta_database/Data_AIRS_apr08/Solar/solar_m140x.txt';
FNSUN = 'solar_m140x.txt';
sarta = load(FNSUN);  %% effectively ttorad(nu,5600) in W/m2/sr/cm-1  so need to do the solidangle conversion,

 figure(4); plot(x,y*solidangle,fkc,rkc*1000*solidangle,sarta(:,2),sarta(:,3)*solidangle*1000)
figure(4); plot(x,y*solidangle,fkc,rkc*1000*solidangle,sarta(:,2),sarta(:,3)*solidangle*1000,'.'); 
  axis([600 3000 0 20])
xlabel('wavenumber cm-1'); ylabel('TOA radiance mW/m2/sr/cm-1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/
figure(5)
[sergioL,sergioR] = rads_wnum2lamda(fkc,rkc);
plot(sergioL,sergioR*1000*solidangle,'k.',a(:,1),a(:,2),'r');
xlabel('wavelength um'); ylabel('TOA radiance W/m2/sr/um')
axis([0 4 0 2500])