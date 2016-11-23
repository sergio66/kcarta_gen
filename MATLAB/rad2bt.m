%
% NAME
%   rad2bt - translate radiance to brightness temperature
%
% SYNOPSIS
%   bt = rad2bt(fr, rad);
%
% INPUTS
%   fr  - n-vector of wavenumbers, cm-1
%   rad - n x k array of radiances, mW/m2 per strad
%
% OUTPUT
%   bt  - n x k array of brightness temps, K
%
% DISCUSSION
%   radiance units are milliwatts per square meter per steradian.
%
%   fr can be a row or column vector, and rad a vector or array.
%   When the length of fr equals the number of rows of rad, it is
%   applied along columns.  If fr is scalar it is applied to all
%   elements of rad, and if rad is scalar it is applied to all
%   elements of fr.  If rad is a row vector with the same length 
%   as fr, then for this one case fr is applied along the row.
%
% AUTHOR
%   H. Motteler
%

function bt = rad2bt(fr, rad);

% Constants; values from NIST (CODATA98)
c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1

% Compute radiation constants c1 and c2
c1 = 2*h*c*c * 1e+11;
c2 = (h*c/k) * 100;

% the scalar calculation, for reference
% bt = c2 * fr / log(1 + c1 * fr^3 / rad)

% set up the data
d1 = size(fr);              % fr dimension list
fr = fr(:);                 % make fr a column vector
k = length(fr);             % fr vector length
d2 = size(rad);             % rad dimension list
[m,n] = size(rad);          % rad size as a 2d array
rad = reshape(rad, m, n);   % treat rad as a 2d array
j = n;                      % copy fr across j columns

% fr and rad special cases 
if k ~= m                   % length of fr != rows of rad
  j = 1;                    % keep fr as a single column
  if k > 1 && m==1 && n==1  % fr is a vector and rad a scalar
    d2 = d1;                % reshape rad at the end to match fr
  elseif k == n && m == 1   % rad is a row vec with same len as fr
    rad = rad(:);           % reshape rad as a column, to match fr
  elseif k ~= 1             % if k == 1 then fr is scalar, otherwise...
    error('the length of fr must equal the number of rows of rad')
  end
end

% do the vectorized calculation
bt = c2 * (fr * ones(1,j)) ./ log(1 + c1 * (fr.^3 * ones(1,j)) ./ rad);

% restore rad original shape
bt = reshape(bt, d2);

