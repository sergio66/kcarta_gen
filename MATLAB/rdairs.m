% this script reads in the data file which is the result of the
% kCARTA data convolved with the specified function
%
% note the only gases whose d/dq are output are those whose ID <= 6
% useful matrices are : raaAmt  : gas amounts  
%                      raaTemp : gas temperatures
%                      raCenterFreq = SRF center frequencies
%                      raaConData   = convolved output
%
% if there are N layers and G gases, then there will be (N+1)G+3
% Jacobians that are output : N for each of the G gases, N for the
% temperature and one each for Jacobians wrt sureface temp, surface
% emissivity and background thermal
%
% FORTRAN to MATLAB conversion of the binary file FILENAME
% automatically adds  .datCON to the FILENAME

%binary data in form
%iS,iE,iStore  iS=start index,iE=end index of AIRS (1 <= iS <= iE <= 2372)
%              iStore = number of spectra to expect 
%                              (1 if iChoice=2 in readatmos.f
%                              ny if iChoice=2 in readatmos.f)
%raCenterFreq(iJ)iJ=iS,iE       AIRS center channel freq
%followed by data, for iK=1,iStore
%raConvolved(iK,iJ)iJ=iS,iE

%infix=input('Enter INPUT binary file name (w/o .datCON extension) ','s');
%extn='.datCON';
%filename=[infix extn];
%fid=fopen([filename],'r','ieee-be');

filename = input('Enter INPUT binary file name ','s');
fid = fopen([filename],'r','native');

%first read iS,iE,iConSaved
fmjunk1   = fread(fid,1,'integer*4');
iS        = fread(fid,1,'integer*4');
iE        = fread(fid,1,'integer*4');
iConSaved = fread(fid,1,'integer*4');
fmjunk1   = fread(fid,1,'integer*4');

%now read CenterFreqs
iBlSz   = iE-iS+1;
raCFreq = zeros(iBlSz,1);
fmjunk1 = fread(fid,1,'integer*4');
ra      = fread(fid,iBlSz,'real*4');
fmjunk1 = fread(fid,1,'integer*4');
raCFreq = ra;

%read in the convolved data
raaCData = zeros(iBlSz,iConSaved);
for ii = 1:iConSaved
  fmjunk1        = fread(fid,1,'integer*4');
  ra             = fread(fid,iBlSz,'real*4');
  fmjunk1        = fread(fid,1,'integer*4');
  raaCData(:,ii) = ra;
  end 
fclose(fid);

raaCData = raaCData';
raCFreq  = raCFreq';

clear ans extn fid filename fmjunk1 iBlSz iConSaved iE iS ii infix ra;
