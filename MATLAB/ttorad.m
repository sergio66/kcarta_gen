function rad = bt2rad(fr, bt);

% function rad = bt2rad(fr, bt);
%
% Translate brightness temperatures to radiances
%
% Input: 
%    fr  - vector of wavenumber (cm-1)
%    bt  - vector or column-order array of brightness temperature (Kelvin)
%
% Output:
%    rad - vector or column-order array of radiance (mW/m2 per cm-1 per strad)
%
    
% H. Motteler, 8/27/98
% Update: 24 Jun 02 Scott Hannon - use milliwatts instead of Watts; adjust
%    Planck and Boltzmann to use NIST/CODATA98 values
% Note: constants are NIST/CODATA98 values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1


% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
c2 = (h*c/k) * 100;


% return rad = c1 * fr^3 / (exp(c2*fr/bt) - 1)

fr = fr(:);
[m,n] = size(bt);
if m == 1
  bt = bt';
  n = 1;
end

rad  = c1 * (fr.^3 * ones(1,n)) ./ (exp((c2 * fr * ones(1,n)) ./ bt) - 1);

%%% end of function %%%%
