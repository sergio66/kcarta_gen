%%% this reads in CO2 370 ppmv stuff, and finally outputs CO2 400 ppmv stuff
%%% to update this, simply change 400 to your new number

dirx = '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN/fbin/etc.ieee-le/400ppmv/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end
dirx = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/400ppmv/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end

%%% ****** be careful about this; the new ref profs are for July2010 ******
dirP0 = '/asl/data/kcarta/KCARTADATA/RefProf_July2010.For.v115up_CO2ppmv385/';
dirPX = '/asl/data/kcarta/KCARTADATA/RefProf_July2010.For.v115up_CO2ppmv400/';
ratioPX = 400/385;
ee = exist(dirPX,'dir');
if ee == 0
  mker = ['!mkdir ' dirPX]; eval(mker);
  end
%% copy the files over, so that you always have a fresh refgas copy so that 
%% you can use the following later, without worrying about whether you are
%% mutiplyinrefgas2 by ratioPX again and again and again and again ....
%%   haha(:,3) = haha(:,3)*ratioPX;
%%   haha(:,5) = haha(:,5)*ratioPX;
cper = ['!/bin/cp ' dirP0 '* ' dirPX '.'];
eval(cper)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inCO2  = 370;
outCO2 = 400;

gid = 2;
ratio = outCO2/inCO2;
cdir    = '/asl/s1/sergio/CO2ppmv400/MATFILES/';
dtype = 'ieee-le';

dachunk = 605:25:2830;

disp('----> basic stuff');
fdirIN  = '/asl/data/kcarta/v24.ieee-le/co2.ieee-le/';        %% 370 ppmv
fdirNEW = '/asl/s1/sergio/CO2ppmv400/co2.ieee-le/';    %% 400 ppmv
fdirASL = '/asl/data/kcarta/UMBC_CO2_H1998.ieee-le/CO2ppmv400.ieee-le/';

rmer = ['!/bin/rm ' cdir '*.mat']; eval(rmer);
rmer = ['!/bin/rm ' fdirNEW '*.dat']; eval(rmer);

dirx = '/asl/s1/sergio/CO2ppmv400/';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end
dirx = '/asl/s1/sergio/CO2ppmv400/MATFILES';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end
dirx = '/asl/s1/sergio/CO2ppmv400/co2.ieee-le';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end
dirx = /asl/data/kcarta/UMBC_CO2_H1998.ieee-le/CO2ppmv400.ieee-le';
ee = exist(dirx,'dir');
if ee == 0
  mker = ['!mkdir ' dirx]; eval(mker);
  end

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

%% now copy from fdirNEW to fdirASL, so this can be used in kcarta.param
cper = ['!/bin/cp ' fdirNEW '* ' fdirASL '.']; eval(cper)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('----> background upper stuff');
fdirIN  = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/370ppmv/'
fdirNEW = ...
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN_UPPER/fbin/etc.ieee-le/400ppmv/'

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
  '/asl/data/kcarta/BACKGND_11SEPT2005_COUSIN/fbin/etc.ieee-le/400ppmv/'

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
%% also need to fix the refprof files
%% if you;re not careful and keep running this code, you will keep increasing
%% co2 amount by ratioPX or by ratio .......
error('ack')

%%% MAIN ONE -------------------- made from 385 ppmv ----------------------
cder = ['cd ' dirPX]; eval(cder);
haha = load('refgas2');

haha(:,3) = haha(:,3)*ratioPX;
haha(:,5) = haha(:,5)*ratioPX;

fid = fopen('refgas2','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% BACKGND : USUAL LAYERS  --------------------- made from 370 ppmv ---------
cd /asl/data/kcarta/KCARTADATA/NLTE/USUALLAYERS
haha = load('refgas2_370ppmv');

haha(:,3) = haha(:,3)*ratio;
haha(:,5) = haha(:,5)*ratio;

fid = fopen('refgas2_400ppmv','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% BACKGND : UA  --------------------- made from 370 ppmv ---------
cd /asl/data/kcarta/KCARTADATA/NLTE/UA
haha = load('refgas2_370ppmv');

haha(:,3) = haha(:,3)*ratio;
haha(:,5) = haha(:,5)*ratio;

fid = fopen('refgas2_400ppmv','w');
fprintf(fid,'%2i %11.4e %11.4e %8.3f %11.4e \n',haha');
fclose(fid)

%%% make the links
cd /asl/data/kcarta/KCARTADATA/General
lnser = ['!ln -s ../NLTE/USUALLAYERS/refgas2_400ppmv refgas2Back_400ppmv'];
eval(lnser)
lnser = ['!ln -s ../NLTE/UA/refgas2_400ppmv refgas2BackUA_400ppmv'];
eval(lnser)
cper = ...
  ['!cp ../RefProf_July2010.For.v115up_CO2ppmv400/refgas2 refgas2_400ppmv'];
eval(cper)

cd /home/sergio/KCARTA/UTILITY