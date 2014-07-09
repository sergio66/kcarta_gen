cd /home/sergio/SPECTRA

%%% need to modify run7co2 so it dumps out lineORIG into 
%%% save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines.mat  lineORIG  or
%%% save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines2.mat lineORIG  or
%%% save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines3.mat lineORIG  or
%%% save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines4.mat lineORIG  or
%%% save /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/otherlines5.mat lineORIG  
 
fname = 'IPFILES/co2one';
clear topts
topts.band = [2310      2320      2350 2351 2352     ];             %% X
topts.band = [2310 2311 2320 2321 2350 2351 2352     ];             %% 2
topts.band = [2310 2311 2320 2321 2350 2351 2352            2355];  %% 3
topts.band = [2310      2320      2350 2351      2353 2354];        %% 4
topts.band = [                    2350                    ];        %% 5 
run7co2(2,2205,2505,fname,topts);

%%then just cut and paste appropriate stuff from KCARTA/UTILITY/lineparameters
cd /home/sergio/KCARTA/SRCv1.11/NONLTE/
cd /carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/
%topts.band = [2310      2320      2350 2351 2352  ];     load otherlines.mat
%topts.band = [2310 2311 2320 2321 2350 2351 2352   ];    load otherlines2.mat
%topts.band = [2310 2311 2320 2321 2350 2351 2352 2355];  load otherlines3.mat
%topts.band = [2310 2320 2350 2351 2353 2354];            load otherlines4.mat
topts.band = [          2350               ];            load otherlines5.mat
lines = lineORIG; clear lineORIG

idgas = 2;

semilogy(lines.wnum,lines.stren,'.'); grid 
fname = '/carrot/s1/sergio/AIRSCO2/CO2_BANDS_PARAM/weakco2lines5.dat';
fid=fopen(fname,'w','ieee-le');

% Write header info
% header1 = gasid, npts, isotope

isotope = lines.iso;
stren   = lines.stren;
elower  = lines.els;
freq    = lines.wnum;
j_lower = lines.ilsgq;
j_upper = lines.iusgq;
p_shift = lines.tsp;
w       = lines.abroad;
w_s     = lines.sbroad;
w_temp  = lines.abcoef;

npts = length(stren);
filemark= 4 + 4;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[idgas,npts],'integer*4');
fwrite(fid,filemark,'integer*4');

%Save the needed vectors 
filemark= 8 * npts;

fwrite(fid,filemark,'integer*4');
fwrite(fid,isotope,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,elower,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,freq,'real*8');
fwrite(fid,filemark,'integer*4');

lowerj = j_lower * 1.0;
fwrite(fid,filemark,'integer*4');
fwrite(fid,lowerj,'real*8');
fwrite(fid,filemark,'integer*4');

upperj = j_upper * 1.0;
fwrite(fid,filemark,'integer*4');
fwrite(fid,upperj,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,p_shift,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,stren,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w_s,'real*8');
fwrite(fid,filemark,'integer*4');

fwrite(fid,filemark,'integer*4');
fwrite(fid,w_temp,'real*8');
fwrite(fid,filemark,'integer*4');

fclose(fid);
