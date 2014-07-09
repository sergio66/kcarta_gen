
function gpro2rtp(gdir, gpro, rtpf)

% NAME
%
%   gpro2rtp -- translate genln2 profiles to an RTP profile set
% 
% SYNOPSIS
%
%   gpro2rtp(gdir, gpro, rtpf) 
%
% INPUTS
%
%   gdir - GENLN2 profile directory (prefix for filenames)
%   gpro - GENLN2 profile filename, or cell array of filenames
%   rtpf - RTP output filename
%
% OUTPUT
%
%   an RTP file containing the listed GENLN2 profiles
%
%   gpro2rtp fills out the RTP header fields glist, gunit, pmin,
%   pmax and the profile fields plat, psurf, nlevs, plevs, ptemp,
%   and gamnt
%

if strcmp(class(gpro), 'char')
  gpro = {gpro};
end

nprof = length(gpro);

% -----------------------------------------------
% read the GENLN profiles into a structure array
% -----------------------------------------------
clear prof
pmax = 0;
pmin = 10000;
for i = 1 : nprof

  % read a GENLN2 format profile
  pname = gpro{i};
  pfile = [gdir, pname];
  [plevs, ptemp, glist, gamnt, plat] = gproread(pfile);

  pmax = max([plevs(:); pmax]);
  pmin = min([plevs(:); pmin]);

  if i == 1
    glist1 = glist;
  elseif glist1 ~= glist
    i, glist, glist1
    error('all profiles in an RTP file must all have the same gas set')
  end

  prof(i).plat = plat;
  prof(i).spres = max(plevs);
  prof(i).nlevs = length(plevs);
  prof(i).plevs = plevs(:);
  prof(i).ptemp = ptemp(:);
  prof(i).gamnt = gamnt;

end

head.ptype = 0;
head.pfields = 1;
head.glist = glist;
head.gunit = ones(length(head.glist), 1) * 10; % 10 is PPMV
head.pmin = pmin;
head.pmax = pmax;

% profile attributes
pattr = { {'plevs', 'units', 'millibars'}, ...
          {'ptemp', 'units', 'kelvins'}};

% header attributes
hattr = { {'header', 'title', 'gpro2rtp.m translation of GENLN2 profiles'}, ...
          {'header', 'date', date}, ...
	  {'header', 'gas units', 'PPMV'}};

% write the profile in RTP format
rtpwrite2(rtpf, head, hattr, prof, pattr);

