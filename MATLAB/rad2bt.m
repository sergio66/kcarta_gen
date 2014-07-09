function bt = rad2bt(fr, rad);

% function bt = rad2bt(fr, rad);
%
% Translate radiances to brightness temperatures
%
% Input: 
%    fr - vector of wavenumber (cm-1)
%    rad - vector or column-order array of radiances (mW/m2 per cm-1 per strad)
%
% Output
%    bt - vector or column-order array of brightness temp's
%

% H. Motteler, 8/27/98
% Update: 24 Jun 02 Scott Hannon - use milliwatts instead of Watts; adjust
%    Planck and Boltzmann to use NIST/CODATA98 values.  Also renamed "rd"
%    to "rad" to match name used in "bt2rad.m".
% Update: 27 May 2009, Breno Imbiriba - modify dimension checks to correct
%    inconsistency when f is [1 x 1].
% Note: constants are NIST/CODATA98 values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
c2 = (h*c/k) * 100;

% return bt = c2 * fr / log(1 + c1 * fr^3 / rad)

[m1 n1]=size(fr);
[m2 n2]=size(rad);
if(n1~=1)
  fr=fr';
  [m1 n1]=size(fr);
end
if(m2~=m1)
  if(n2==m1)
    rad=rad';
  else
    error(['Wrong argument size']);
  end
end
n=n2;

bt = c2 * (fr * ones(1,n)) ./ log(1 + c1 * (fr.^3 * ones(1,n)) ./ rad);

