%% data from http://rredc.nrel.gov/solar/spectra/am0/

aVIS = load('E490_00a_AM0_short.csv');
aNIR = load('E490_00a_AM0_SHORT_IR.csv');
whos

a(:,1) = [aVIS(:,1); aNIR(:,1)];
a(:,2) = [aVIS(:,2); aNIR(:,2)];

xstartup
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS/

x = 500 : 100 : 100000; y = ttorad(x,5900);
[lout,rout] = rads_wnum2lamda(x,y);           

solidangle = 6.785087652174316e-5;

figure(1);
plot(a(:,1),a(:,2),lout,rout*solidangle)
axis([0 10 0 2500])
xlabel('wavelength um');

[boo,woo] = rads_lamda2wnum(a(:,1),a(:,2));

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

iCnt = 0;
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
%%%%%%%%% if ystart(iCnt) > 22000    %%% orig plan was to only do UV
    if ystart(iCnt) >= 2830          %%% but have decided to update NIR onwards, esp since there were many gaps
      iNeed = iNeed + 1;

      fprintf(1,'looks like we need this solar file %5i \n',ystart(iCnt))
      freqs = ystart(iCnt) : df : yend(iCnt)-df;
      whos freqs
      solarrads = interp1(boo,woo,freqs); figure(3); plot(freqs,solarrads,'o-'); 
      title(num2str(iNeed)); pause(0.1)

      fname = ['rad_solar' num2str(ystart(iCnt)) '.dat'];
      fid=fopen(fname,'w','ieee-le');
      fs=freqs(1);
      fe=freqs(end);
      fprintf(1,'%s \n',fname);
      fprintf(1,'  fs,fe,length(freqs) = %10.5f %10.5f %6i \n \n',fs,fe,length(freqs)')
      
      % header = fstart,fstep
      filemark= 8 + 8 + 8;
      fwrite(fid,filemark,'integer*4');
      fwrite(fid,[fs,fe,df],'real*8');
      fwrite(fid,filemark,'integer*4');

      %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
      %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
      %Save the solar radout array
      filemark= 8 * length(freqs);
      fwrite(fid,filemark,'integer*4');
      fwrite(fid,solarrads/(1000 * solidangle),'real*8');
      fwrite(fid,filemark,'integer*4');
      %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle
      %%% KCARTA solar datafiles are in W/m2/sr/cm-1 and then kCARTA internally multiplies by 1000 * solidangle

      fclose(fid);

    end
  end
end

! ls -lt rad_solar2????.dat
!pwd
! ls -lt rad_solar1[0-9][0-9][0-9][0-9].dat

figure(3); plot(ystart,'o-'); grid
  line([0 iCnt],[0605 0605]);
  line([0 iCnt],[2830 2830]);