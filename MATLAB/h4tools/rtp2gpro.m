
function rtp2gpro(rtpf, gdir, gpre)

% NAME
%
%   rtp2gpro -- translate an RTP profile set to genln2 profiles
% 
% SYNOPSIS
%
%   rtp2gpro(rtpf, gdir, gpre)
%
% INPUTS
%
%   rtpf - RTP output filename
%   gdir - optional directory for GENLN2 profiles
%   gpre - optional alternate prefix for output names
%
% OUTPUT
%
%   files of the form <gdir>/p<i>, one file per profile
%

if nargin < 3
  gpre = 'p';
end
if nargin < 2
  gdir = './';
end

% read the RTP file
[head, hattr, prof, pattr] = rtpread2(rtpf);

glist = head.glist;

for i = 1 : length(prof)

  pfile = [gdir, gpre, num2str(i)];

  gprowrite(pfile, ...
            prof(i).plevs, ...
	    prof(i).ptemp, ...
	    head.glist,    ...
	    prof(i).gamnt, ...
	    prof(i).plat);
end

