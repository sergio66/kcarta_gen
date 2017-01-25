function [h,ha,p,pa] = levelstext_to_levelsRTP(fTXT,satzenORscanang,satheight);

%% kCARTA can take in a NWP/sonde levels profile in text
%% and process it internally without going through klayers.
%% this function erps in the same file and produces a levels RTP file
%%
%% input fTXT : text file containing levels profile in format (kRTP = -10) the kcarta expects
%%       satzenORscanang : if +ve you sent in satzen, if negative you sent in scanang
%%       satheight : optional; if not supplied set to 705000 m (AIRS height in meters)
%%
%% also see levelsRTP_to_levelstext.m

addpath /home/sergio/MATLABCODE/TIME

fid = fopen(fTXT,'r');

S = fgets(fid); %% comment
S = fgets(fid); %% nlevs
  nlevs = str2num(S);
  p.nlevs = nlevs;
S = fgets(fid); %% spres,stemp,salti
  junk = str2num(S);
  p.spres = junk(1);
  p.stemp = junk(2);
  p.salti = junk(3);
S = fgets(fid); %% year, lat, lon
  junk = str2num(S);
  p.rtime = (junk(1)-1958)*365*24*60*60;
  p.rlat = junk(2);
  p.rlon = junk(3);
  p.plat = p.rlat;
  p.plon = p.rlon;
S = fgets(fid); %% numgases
  ngas = str2num(S);
S = fgets(fid); %% glist
  glist = str2num(S);
  if length(glist) ~= ngas
    error('need same number of gases as stated in ngas')
  end
S = fgets(fid); %% gunit
  gunit = str2num(S);
  if length(gunit) ~= ngas
    error('need same number of gunits as stated in ngas')
  end

p.satheight = 705000;
p.upwell = +1;
p.solzen  = 150;
if nargin == 2
  p.zobs = 705000;
else
  p.zobs = satheight;
end

if abs(satzenORscanang) < eps
  p.satzen = 0.0;
  p.scanang = 0.0;
elseif  satzenORscanang > 0
  p.satzen = satzenORscanang;
  p.scanang  = saconv(p.satzen,p.zobs);
elseif  satzenORscanang < 0
  p.scanang = abs(satzenORscanang);
  p.satzen  = vaconv(p.scanang,p.zobs,0);
end  

fprintf(1,'nlevs = %3i ngas = %3i \n',nlevs,ngas);

for ii = 1 : nlevs  
  S = fgets(fid); 
  junk = str2num(S);
  if length(junk) ~= ngas + 2
    error('need correct number of columns for p,T,gas profile')
  end
  p.plevs(ii) = junk(1);
  p.ptemp(ii) = junk(2);
  junk = junk(3:end);
  for gg = 1 : ngas
    x = junk(gg);
    str = ['p.gas_' num2str(glist(gg)) '(' num2str(ii) ') = x;'];
    eval(str);
  end  
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pa = {{'profiles','rtime','seconds since 1958'}};
ha = {{'header','hdf file',fTXT}};

h.ngas  = ngas;
h.glist = glist';
h.gunit = gunit';
h.vcmin = 605;
h.vcmax = 2830;
h.ptype = 0;    %% levels
h.pfields = 1;  %% profile only
h.pmin = min(p.plevs);
h.pmax = max(p.plevs);

p.nemis = 2;
p.efreq = [500  3000]';
p.emis  = [0.98 0.98]';
p.rho   = (1-p.emis)/pi;

p.ptemp = p.ptemp';
p.plevs = p.plevs';
for gg = 1 : ngas
  x = junk(gg);
  str = ['p.gas_' num2str(glist(gg)) ' = p.gas_' num2str(glist(gg)) ''';' ];
  eval(str);
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%