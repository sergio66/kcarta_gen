function [theta1,theta2,ref_ind] = scanang2satzen(thetaIN,H_sat,p)

% http://mintaka.sdsu.edu/GF/explain/atmos_refr/horizon.html
% http://www.ess.uci.edu/~cmclinden/link/xx/node45.html for TONS of detail

% input p       = local pressures (mb) at which to find angles
%       thetaIN = satellite scanang (ie wrt nadir) in degrees
%       H_sat   = satellite altitude (in km)
% output
%  theta1 = from Scott's vaconv (no refractive index change)
%  theta2 = from Sergio's vaconv (refractive index change)
%       n = refractive index

R = 6367.512; %% earth radius

p = sort(p);   %% make sure p is from TOP (lowest pressure) to BOTTOM (largest pressure)
h = [[H_sat-100: -100 : 100] (p2h(p)/1000)'];
h = [(p2h(p)/1000)'];

max(p2h(p)/1000)

for ii = 1 : length(h)
  x = (R + H_sat)/(R + h(ii)) * sin(thetaIN * pi/180);
  x = min(x,1);
  theta1(ii) = asin(x) * 180/pi;
end

%% 

p0 = 1013.25;
lambda = 10;
n0 = 64.328 + 29498.1/(146-1/lambda/lambda) + 255.4/(41 - 1/lambda/lambda);
n0 = 1 + n0/1e6;
gamma = (n0*n0-1)/(n0*n0+2);

ref_ind = 1.0;   %% in deep space

for ii = 1 : length(h)
  if h >= 90
    ref_ind(ii) = 1;
  else
    mu = p(ii)/p0*gamma;
    ref_ind(ii) = sqrt((2*mu+1)/(1-mu));
  end
  x = (R + H_sat)/(R + h(ii)) * sin(thetaIN * pi/180) / ref_ind(ii);
  x = min(x,1);
  theta2(ii) = asin(x) * 180/pi;
end

figure(1);
  plot(1:length(h),theta1,'bo-',1:length(h),theta2,'r',1:length(h),ones(size(h))*thetaIN,'ko-')
figure(2);
  plot(1:length(h),theta1-theta2)
