function gprowrite(fname, plev, temp, gasid, gasmx, lat)

% function gprowrite(fname, plev, temp, gasid, gasmx, lat)
%
% gprowrite writes a genln2/klayers format "user" (i.e. input) 
% level profile
%
% inputs
%   fname - output profile file name
%   plev  - n-vector of pressures, in mb
%   temp  - n-vector of temperatures, in K
%   gasid - k-vector of HITRAN gas IDs
%   gasmx - n x k array of mixing ratios, in ppmv
%   lat   - profile latitude, in degrees
%
% output
%   a genln2 user profile

[fid, msg] = fopen(fname, 'w');
if fid == -1
  error(msg)
end

nlev = length(plev);
ngas = length(gasid);

fprintf(fid, '! output from writepro.m\n');

fprintf(fid, '''P''  %g  %g\n', nlev, ngas);

fprintf(fid, '%g ', gasid);
fprintf(fid, '\n');

fprintf(fid, '%g  0.0\n', lat);

for i = 1:nlev

  fprintf(fid, '%9.3e %7.3f -1.0 ', plev(i), temp(i));

  fprintf(fid, '%9.3e ', gasmx(i, :));
  fprintf(fid, '\n');

end

fclose(fid);

