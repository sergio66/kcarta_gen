%%% this reads in CO2 370 ppmv stuff, and finally outputs CO2 378 ppmv stuff
%%% to update this, simply change 378 to your new number

dirx = '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN/fbin/etc.ieee-le/378ppmv/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end
dirx = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/378ppmv/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end

dirP0 = '/asl/data/kcarta/KCARTADATA/RefProf.For.v107up_CO2ppmv370/';
dirx  = '/asl/data/kcarta/KCARTADATA/RefProf.For.v107up_CO2ppmv378/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  cper = ['!/bin/cp ' dirP0 '* ' dirx '.'];
  eval(cper)
  end

inCO2  = 370;
outCO2 = 378;

gid = 2;
ratio = outCO2/inCO2;
cdir    = '/carrot/s1/sergio/CO2ppmv378/MATFILES/';
dtype = 'ieee-le';

dachunk = 605:25:2830;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----> basic stuff');
fdirIN  = '/asl/data/kcarta/v24.ieee-le/co2.ieee-le/';
fdirNEW = '/carrot/s1/sergio/CO2ppmv378/co2.ieee-le/';

rmer = ['!/bin/rm ' cdir '*.mat']; eval(rmer);
rmer = ['!/bin/rm ' fdirNEW '*.dat']; eval(rmer);

for ii = 1 : length(dachunk)
  vchunk = dachunk(ii);
  fname = sprintf('%s/r%d_g%d.dat', fdirIN, vchunk, gid);
  inputname = fname; 
  ee = exist(inputname,'file');
 
  if ee > 0
    % inputs
    %   gid    -  gas ID
    %   vchunk -  wavenumber start of 25 1/cm chunk
    %   cdir   -  directory for matlab output files
    %   fdir   -  directory for fortran source files
    %   dtype  -  fortran data type
    for2matconvert(gid, vchunk, cdir, fdirIN, ratio, dtype)

    % inputs
    %   gid    -  gas ID
    %   vchunk -  wavenumber start of 25 1/cm chunk
    %   cdir   -  directory for matlab mat source files
    %   fdir   -  directory for fortran output files
    %   dtype  -  output data type (for matlab open)
    mat2for(gid, vchunk, cdir, fdirNEW, dtype)

    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----> background upper stuff');
fdirIN  = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/370ppmv/'
fdirNEW = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/378ppmv/'

rmer = ['!/bin/rm ' cdir '*.mat']; eval(rmer);
rmer = ['!/bin/rm ' fdirNEW '*.dat']; eval(rmer);

for ii = 1 : length(dachunk)
  vchunk = dachunk(ii);
  fname = sprintf('%s/r%d_g%d.dat', fdirIN, vchunk, gid);
  inputname = fname;
  ee = exist(inputname,'file');
  if ee > 0
    for2matconvert(gid, vchunk, cdir, fdirIN, ratio, dtype)
    mat2for(gid, vchunk, cdir, fdirNEW, dtype)
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----> background usual stuff');
fdirIN  = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN/fbin/etc.ieee-le/370ppmv/'
fdirNEW = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN/fbin/etc.ieee-le/378ppmv/'

rmer = ['!/bin/rm ' cdir '*.mat']; eval(rmer);
rmer = ['!/bin/rm ' fdirNEW '*.dat']; eval(rmer);

for ii = 1 : length(dachunk)
  vchunk = dachunk(ii);
  fname = sprintf('%s/r%d_g%d.dat', fdirIN, vchunk, gid);
  inputname = fname;
  ee = exist(inputname,'file');
  if ee > 0
    for2matconvert(gid, vchunk, cdir, fdirIN, ratio, dtype)
    mat2for(gid, vchunk, cdir, fdirNEW, dtype)
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% also gave to fix the refprof files

error('ack');

%%% MAIN ONE ----------------------------------------------------
cd /asl/data/kcarta/KCARTADATA/RefProf.For.v107up_CO2ppmv378/
haha = load('refgas2');

haha(:,3) = haha(:,3)*ratio;
haha(:,5) = haha(:,5)*ratio;

fid = fopen('refgas2','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% BACKGND : USUAL LAYERS  --------------------------------------------------
cd /asl/data/kcarta/KCARTADATA/NLTE/USUALLAYERS
haha = load('refgas2_370ppmv');

haha(:,3) = haha(:,3)*ratio;
haha(:,5) = haha(:,5)*ratio;

fid = fopen('refgas2_378ppmv','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% BACKGND : UA  --------------------------------------------------
cd /asl/data/kcarta/KCARTADATA/NLTE/UA
haha = load('refgas2_370ppmv');

haha(:,3) = haha(:,3)*ratio;
haha(:,5) = haha(:,5)*ratio;

fid = fopen('refgas2_378ppmv','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% make the links
cd /asl/data/kcarta/KCARTADATA/General
lnser = ['!ln -s ../NLTE/USUALLAYERS/refgas2_378ppmv refgas2Back_378ppmv'];
eval(lnser)
lnser = ['!ln -s ../NLTE/UA/refgas2_378ppmv refgas2BackUA_378ppmv'];
eval(lnser)
cper = ['!cp ../RefProf.For.v107up_CO2ppmv378/refgas2 refgas2_378ppmv'];
eval(cper)