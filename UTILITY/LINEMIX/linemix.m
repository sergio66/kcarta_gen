cd /home/sergio/SPECTRA
fa = 2155; fb = 2455; 

f1 = 2380;
f2 = 2405;

% need to edit this file to put in the correct Temperature!!!!
% afetr this, use linemix_combine.m 
fname = '/home/sergio/SPECTRA/IPFILES/co2one_lowp';

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [-11]; topts.band   = 2350; 
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [+11]; topts.band   = 2350;                          
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [-11]; topts.band   = 2351;                          
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [+11]; topts.band   = 2351; 
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [-11]; topts.band   = 2352;                          
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [+11]; topts.band   = 2352;                          
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [-13]; topts.band   = 2320; %%pipi  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [+13]; topts.band   = 2320; %%pipi  
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [-12]; topts.band   = 2310; %%deltdelt
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/SPECTRA; clear topts;
topts.PQRallowed = [+12]; topts.band   = 2310; %%deltdelt
topts.HITRAN = '/asl/data/hitran/h2k.oldiso';
topts.LVF = 'V'; topts.mainloop = -1; topts.birn = 'n';
[fr,k] = run7co2(2,f1,f2,fname,topts);

cd /home/sergio/KCARTA/SRCv1.11/NONLTE

