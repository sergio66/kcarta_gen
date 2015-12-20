function flux_levels = compute_flux_from_4radstreams(w,rads)

% SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE
% use following angles and weights, see subr FindGauss2old
%  daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
%      = 81.96604706186704     65.39188287931897     43.68425288798865     19.45630246759402
%  daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0
%
% rads = 89000 * (4 * N) levels
% for each level, do
% flux_levels(:,ii) = (daW(1)*d2x(:,1)+daW(2)*d2x(:,2)+daW(3)*d2x(:,3)+daW(4)*d2x(:,4)) * 2*pi
% >> total flux = sum(flux_levels,1) * dw/1000 ;
%
% the nml file has
%  nm_params
%     !! testing no flux
%     kFlux = -1
%  nm_radnce
%   iAtmLoop = 3
%    raAtmLoop(1) = 81.9660
%    raAtmLoop(2) = 65.3919
%    raAtmLoop(3) = 43.6843
%    raAtmLoop(4) = 19.4563
%    raAtmLoop(5) = -10.00
%  nm_output
%    !dump out radiance at all levels
%    iaPrinter(1)   =    3
%    iaGPMPAtm(1)   =   -1
%    iaNp(1)        =   -1

daX = [0.1397599 0.4164096 0.7231570 0.9428958];  %% compute up and downwell rads at thiese angles
daW = [0.0311810 0.1298475 0.2034646 0.1355069];  %% these are the weights

dw = mean(diff(w));

[m,n] = size(rads);

Nlays = floor(n/4);
if abs(Nlays - n/4) > eps
  error('hmm, looks like input number of layers is not divisible by 4')
else
  fprintf(1,'number of layers = %3i \n',Nlays)
end

%% integral over dphi = integral (2pi) --> pi (because there is a factor of 2 in the theta integral)
%% but then the above 4 angles are only over the postivie angles, also need to do negative angles
%% so factor of 2 pi
flux_levels = zeros(m,Nlays);
for jj = 1 : 4
  ind = (1:Nlays) + (jj-1)*Nlays;
  flux_levels = flux_levels + daW(jj)*rads(:,ind)*2*pi;
end

plot(sum(flux_levels,1) * dw/1000,1:Nlays,'linewidth',2); title('flux W/m2 at each level'); grid
xlabel('W/m2'); ylabel('level number')  