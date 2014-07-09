cd /home/sergio/SPECTRA
fa = 2155; fb = 2455; 

f1 = 2355;
f2 = 2455;
fname = '/home/sergio/SPECTRA/IPFILES/std_co2_last5';

clear topts;
topts.birn = 'n';
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2350;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2350nodoc] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355nodoc.mat

clear topts;
topts.birn = 'n';
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2351;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2351nodoc] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355nodoc.mat

clear topts;
topts.birn = 'n';
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2352;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2352nodoc] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355nodoc.mat

clear topts;
topts.birn = 'n';
topts.PQRallowed = [-13 +13];                   %%pipi  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2320;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2320nodoc] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355nodoc.mat

clear topts;
topts.birn = 'n';
topts.PQRallowed = [-12 +12];                  %%deltdelt
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2310;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2310nodoc] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355nodoc.mat

cd /home/sergio/KCARTA/SRCv1.11/NONLTE/
