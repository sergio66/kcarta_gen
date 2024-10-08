function ht = p2h(pin)
%this function takes pressure in mb, and converts to height in m

%load /asl/matlab/kcarta/airsheights.dat
%load /asl/matlab/kcarta/airslevels.dat
load /home/sergio/MATLABCODE/airsheights.dat
load /home/sergio/MATLABCODE/airslevels.dat

h = airsheights;
p = airslevels;
for ii=1:100
  pavg(ii) = (p(ii+1)-p(ii))/log(p(ii+1)/p(ii));
end
%plot(h,pavg)

%ht = interp1(pavg,h,pin);
ht = interp1(pavg,h,pin,[],'extrap');

% boo = find(~isfinite(ht) | ht > 7.5e4);
% ht(boo) = 7.5e4;
% if ((isnan(ht)) | ht > 7.05e4)
%   ht = 8.09e4;
% end

hmax = 90e3;  %% 90 km
boo = find(~isfinite(ht) | ht > hmax);
ht(boo) = hmax;

