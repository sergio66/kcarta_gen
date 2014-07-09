%% careful : we might need to replace upper 15 layers!

iStartNLTE = input('Enter how many upper KCARTA layers to replace! ');

load /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2155.mat
fra = fr;
k2310a = k2310;
k2320a = k2320;
k2350a = k2350;
k2351a = k2351;
k2352a = k2352;

load /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2255.mat
frb = fr;
k2310b = k2310;
k2320b = k2320;
k2350b = k2350;
k2351b = k2351;
k2352b = k2352;

load /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat
frc = fr;
k2310c = k2310;
k2320c = k2320;
k2350c = k2350;
k2351c = k2351;
k2352c = k2352;

clear fr k2310 k2320 k2350 k2351 k2352

fr    = [fra    frb    frc];
k2310 = [k2310a k2310b k2310c];
k2320 = [k2320a k2320b k2320c];
k2350 = [k2350a k2350b k2350c];
k2351 = [k2351a k2351b k2351c];
k2352 = [k2352a k2352b k2352c];

%%%% now choose which ones i want to add together!
%%%% /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_*.dat
clear k00
len = length(fr);

iOrig  = 100-iStartNLTE;
iaOrig = 1 : iOrig;
iaNew  = iOrig+1 : 100;

k00(iaOrig,1:len) = zeros(iOrig,len);
k00(iaNew,1:len) = k2350 + k2351 + k2352;

%%then load all the CO2 spectra, from kCARTA database, from 2155 to 2455
[backgnd,w] = readkcstd('/taro/s1/sergio/AIRSCO2/NONLTE/co2_spectra.dat');
k00 = backgnd' - k00;
k00(iaOrig,1:len) = zeros(iOrig,len);

whos fr k00
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/k23505152.mat fr k00
plot(fr,k00(iOrig+1,:),fr,k2350(1,:)); pause(1)

cd /home/sergio/KCARTA/UTILITY
umbclbl_2_kcarta
%Enter start chunk freq (>=  605 cm-1) : 2155
%Enter stop  chunk freq (<= 2805 cm-1) : 2430
%Enter input (.mat) filename : '/taro/s1/sergio/AIRSCO2/NONLTE/CO2/k23505152'
%Enter output prefix (.dat) filename : '/taro/s1/sergio/AIRSCO2/NONLTE/CO2'
%Enter gasID : 2

cd /home/sergio/KCARTA/SRCv1.11/NONLTE