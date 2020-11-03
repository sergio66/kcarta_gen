function [rad,jac,w] = readkcPCLSAM_coljac(radfile,jacfile,kLongOrShort);

addpath /home/sergio/KCARTA/MATLAB

if nargin == 2
  kLongOrShort = 0;
end

if kLongOrShort == 0
  [rad,w]  = readkcstd(radfile);
else
  [rad,w]  = readkcBasic(radfile);
end

[xjac,w] = readkcBasic(jacfile);

[mm, nn]  = size(rad);  %% assuming TOA output, this will give number of subpxels+1
[mmj,nnj] = size(xjac);  %% now have to add these up nicely

subpixels = nn-1;   %% this is number of subpixels

numjacs = nnj/nn;   %% this is number of Q(1),Q(2) .. T,stemp col jacs

if mm ~= mmj
  error('oops number of wavnumber points in radiance/jac are different')
end

fprintf(1,'there are %2i subpixels (from rads) and hence %2i coljacs \n',subpixels,numjacs)

ind0 = (1:numjacs);
jac = zeros(mm,numjacs);
for ii = 1 : subpixels
  ind = (ii-1)*numjacs + ind0;
  jac = jac + xjac(:,ind);
end


