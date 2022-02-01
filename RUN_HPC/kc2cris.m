%% copied from /asl/packages/airs_decon/source
%
% NAME
%   kc2cris - convolve kcarta to CrIS channel radiances
%
% SYNOPSIS
%   [rad2, frq2] = kc2cris(user, rad1, frq1, opts)
%
% INPUTS
%   user  - CrIS user grid struct
%   rad1  - kcarta radiances, m x n array
%   frq1  - kcarta frequencies, m-vector
%   opts  - optional parameters
%
% user fields
%   v1  - user grid low
%   v2  - user grid high
%   vr  - out-of-band rolloff
%   dv  - user grid dv
%
% opts fields
%   pL  - passband low (default user.v1)
%   pH  - passband high (default user.v2)
%   rL  - LHS rolloff width (default user.vr)
%   rH  - RHS rolloff width (default user.vr)
%   ng  - number of guard channels (default 0)
%
% OUTPUTS
%   rad2  - CrIS radiances, k x n array
%   frq2  - CrIS frequency grid, k-vector
%
% DISCUSSION
%   kc2cris and kc2iasi are similar internally--see finterp.pdf 
%   from the airs_decon or ccast documentation for the relevant
%   derivations.  rad1 is an m x n array of kcarta radiances and
%   can quickly reach memory limits for large n.
%
% HM, 22 Oct 2014
%

function [rad2, frq2] = kc2cris(user, rad1, frq1, opts)

% set defaults
pL = user.v1;  % passband low
pH = user.v2;  % passband high
rL = user.vr;  % LHS rolloff width
rH = user.vr;  % RHS rolloff width
ng = 0;        % no guard channels

% check opts fields
if nargin == 4
  if isfield(opts, 'pL'), pL = opts.pL; end
  if isfield(opts, 'pH'), pH = opts.pH; end
  if isfield(opts, 'rL'), rL = opts.rL; end
  if isfield(opts, 'rH'), rH = opts.rH; end
  if isfield(opts, 'ng'), ng = opts.ng; end
end

% check that array sizes match
frq1 = frq1(:);
[m, nobs] = size(rad1);
if m ~= length(frq1)
  error('rad1 and frq1 do not conform')
end

%-----------------------------------
% set up interferometric parameters
%-----------------------------------

dv1 = 0.0025;   % kcarta dv
dv2 = user.dv;  % user grid dv

% check input frequency spacing
if abs(dv1 - (frq1(2) - frq1(1))) > 1e-10
  error('input frequency spacing not 0.0025 1/cm')
end

% get rational approx to dv1/dv2
[m1, m2] = rat(dv1/dv2);
% if ~isclose(m1/m2, dv1/dv2, 4)
%   error('no rational approximation for dv1 / dv2')
% end

% get the tranform sizes
vb = pH + rH;
for k = 4 : 24
  if m2 * 2^k * dv1 >= vb, break, end
end
N1 = m2 * 2^k;
N2 = m1 * 2^k;

% get (and check) dx
dx1 = 1 / (2*dv1*N1);
dx2 = 1 / (2*dv2*N2);
% if ~isclose(dx1, dx2, 4)
%   error('dx1 and dx2 are different')
% end
dx = dx1;

% fprintf(1, 'kc2cris: N1 = %7d, N2 = %5d, dx = %6.3e\n', N1, N2, dx);

%-------------------------------
% take kcarta to CrIS radiances
%-------------------------------

% set kcarta radiance passband to the user grid
rad1 = bandpass(frq1, rad1, pL, pH, rL, rH);

% embed kcarta radiance in a 0 to Vmax grid
ftmp = (0:N1)' * dv1;
rtmp = zeros(N1+1, nobs);
[ix, jx] = seq_match(ftmp, frq1);
rtmp(ix, :) = rad1(jx, :);

% radiance to interferogram
igm1 = real(ifft([rtmp; flipud(rtmp(2:N1, :))]));
igm1 = igm1(1:N1+1, :);

% apply a time-domain apodization
% dtmp = (0:N2)' * dx;
% apod = gaussapod(dtmp, 2) * ones(1, nobs);
% igm1(1:N2+1, :) = igm1(1:N2+1, :) .* apod;

% interferogram to radiance
rad2 = real(fft([igm1(1:N2+1,:); flipud(igm1(2:N2,:))]));
frq2 = (0:N2)' * dv2;

% return the user grid plus any guard channels
dg = ng * dv2;
ix = find(user.v1 - dg <= frq2 & frq2 <= user.v2 + dg);
rad2 = rad2(ix, :);
frq2 = frq2(ix);

