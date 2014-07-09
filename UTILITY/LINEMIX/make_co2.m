cd /home/sergio/SPECTRA
fa = 2155; fb = 2455; 

f1 = 2355;
f2 = 2455;
%%% assuming NLTE starts at 60 km (upper 5 kCARTA layers)
fname = '/home/sergio/SPECTRA/IPFILES/std_co2_last5';

%%% looking carefully at Dave Edwards file /SRCv1.11/NONLTE/midlat_day.dat
%%% we see that NLTE starts in the upper 10 kCARTA layers; for fun
%%% just replace the upper 15 kCARTA layers
%%% ie even near the tropopause!!!!!!
fname = '/home/sergio/SPECTRA/IPFILES/std_co2_last15';

clear topts;
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2350;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2350] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat

clear topts;
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2351;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2351] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat

clear topts;
topts.PQRallowed = [-11 +11];                  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2352;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2352] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat

clear topts;
topts.PQRallowed = [-13 +13];                   %%pipi  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2320;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2320] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat

clear topts;
topts.PQRallowed = [-12 +12];                  %%deltdelt
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.band   = 2310;                          
topts.LVF = 'V';
topts.mainloop = -1; 
[fr,k2310] = run7co2(2,f1,f2,fname,topts);
save /taro/s1/sergio/AIRSCO2/NONLTE/CO2/co2_2355.mat

cd /home/sergio/KCARTA/SRCv1.11/NONLTE/
