function [idchan, freq, coef] = rd_nte_le(fname,numcoeff);

% function [idchan, freq, coef] = rd_nte(fname,numcoeff);
%
% Reads in binary FORTRAN non-LTE data file created by program "wrt_nte.m".
% same as /home/hannon/Fit_deltaR_nonLTE/Src/rd_nte.m and expects one of
% Scott Hannon's big-endian files
%
% Input:
%    fname = {string} name of binary FORTRAN data file to read
%    numcoeff [optional, default 7] : Scott's NLTE model had 6 coeffs,
%                                     later became 7 
%
% Output:
%    idchan   = [nchan x 1] channel ID
%    freq     = [nchan x 1] center frequency
%    coef     = [nchan x 6] non-LTE coeffcients
%

% Created: 15 March 2005, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1
  numcoeff = 7;
end

% Expected value of ifm = 4*(1 + 1 + 6) = 32
% This is 4 bytes each for (idchan, freq, & 6 coefs)
ifm_exp=32;

% Expected value of ifm = 4*(1 + 1 + 7) = 36
% This is 4 bytes each for (idchan, freq, & 7 coefs)
ifm_exp=36;

ifm_exp = 4*(1+1+numcoeff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code below should not require modifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = dir(fname);
if d.bytes == 0
  fname
  error('oh oh file DNE')
end

junk = d.bytes/(ifm_exp + 8);
nchan = round(junk);
if (abs(junk - nchan) > 0.0001)
   error('Unexpected file size')
end


% Open output file
fid=fopen(fname,'r','ieee-le');


% Dimension output arrays
idchan=zeros(nchan,1);
freq=zeros(nchan,1);
coef=zeros(nchan,6);


% Loop over the channels
for ic=1:nchan

   % Read FORTRAN start-of-record marker
   [ifm,count]=fread(fid,1,'integer*4');
   if (count == 0)
      ic
      disp('The FORTRAN data file contains fewer channels than expected')
   end
   if (ifm ~= ifm_exp)
      count
      ifm
      ifm_exp
      error('FORTRAN start-of-record marker is unexpected size')
   end

   % Read data for this channel
   idchan(ic)=fread(fid,1,'integer*4');
   freq(ic)=fread(fid,1,'real*4');
   coef(ic,:)=fread(fid,[1,6],'real*4');

   % Read FORTRAN end-of-record marker
   ifm=fread(fid,1,'integer*4');
   if (ifm ~= ifm_exp)
      ifm
      ifm_exp
      error('FORTRAN end-of-record marker is unexpected size')
   end

end % end of loop over bands

fclose(fid);

%%% end of function %%%
