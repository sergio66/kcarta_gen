
function [plev, temp, gasid, gasmx, lat] = gproread(fname)

% function [plev, temp, gasid, gasmx, lat] = gproread(fname)
%
% gproread reads a genln2/klayers format "user" (i.e. input) 
% level profile
% 
% input
%   fname - profile file name
%
% outputs
%   plev  - n-vector of pressures, in mb
%   temp  - n-vector of temperatures, in K
%   gasid - k-vector of HITRAN gas IDs
%   gasmx - n x k array of mixing ratios, in ppmv
%   lat   - profile latitude
%
% H. Motteler, 11 Apr 00
% last modified 5 Feb 01

[fid, msg] = fopen(fname, 'r');
if fid == -1
  error(['error opening ', fname]);
end

k = 0;
tline = fgetl(fid);
while tline(1) == '!'
  tline = fgetl(fid);
  k = k + 1;
end

frewind(fid);

for i = 1:k
  fgetl(fid);
end

% read first data line
porh = fscanf(fid, '%s', 1);
nlev = fscanf(fid, '%g', 1);
ngas = fscanf(fid, '%g', 1);

% read second data line
gasid = fscanf(fid, '%g', ngas);

% The following seems to be necessary to flush a dangling
% <EOL>, for the subsequent fgetl() to work.  The possibly
% more logical fscanf(fid, '\n') doesn't work.

tline = fgetl(fid);
if ~isempty(deblank(tline))
   error('gproread(): hoplessly confused')
end

line = fgetl(fid);
A = sscanf(line, '%g');
lat = A(1);

hflag = porh(2) == 'H';

pdata = fscanf(fid, '%g', [3+hflag+ngas, nlev])';

plev = pdata(:, hflag+1);
temp = pdata(:, hflag+2);
gasmx = pdata(:, hflag+4:hflag+4+ngas-1);

fclose(fid);

