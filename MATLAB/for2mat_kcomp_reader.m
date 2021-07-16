function [gid,fr,kcomp,B] = for2mat_reader(filename,dfchunk);

% function for2mat(gid, vchunk, cdir, fdir, dtype)
%
% produce .mat files from old-format kcarta fortran data
%
% inputs
%   filename - ieee-le compressed data file
%   dfchunk  - optional  chunksize for 10000 pts, assumes 25cm-1
% outputs
%   gid    -  gas ID
%   fr     -  wavenumbers
%   kcomp  - compressed OD
%   B      - basis vectors

if nargin == 1
  dfchunk = 25;
end

% default output data type
dtype = 'ieee-le';

fname = filename;

fid=fopen(fname, 'r' ,dtype);
if fid == -1
  error('can''t open input file')
end

ktype = 2;
npts  = 10000;
fstep = 0.0025;
nlay  = 100;
% [npts,nd]=size(B);

% temperature offsets
ntemp = 11;
toff  = -50.0:10.0:50.0;

% read header info
filemark  = 4 + 8 + 8 + 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4;
filemark2 = fread(fid, 1, 'integer*4');
if filemark ~= filemark2, error('header read format error'), end
gid = fread(fid, 1, 'integer*4');
  fprintf(1,'reading in info for gid = %3i \n',gid);

tmp = fread(fid, 2, 'real*8');
vchunk = tmp(1); fstep2 = tmp(2);
  fprintf(1,'reading in info for chunk = %8.3f \n',vchunk);

tmp = fread(fid, 8, 'integer*4');
npts = tmp(1); nd=tmp(4); ntemp2 = tmp(5);
if ntemp ~= ntemp2, error('ntemp mismatch'), end
filemark2 = fread(fid, 1, 'integer*4');

% read the temperature offsets
filemark= 8 * ntemp;
filemark2 = fread(fid, 1, 'integer*4');
if filemark ~= filemark2, error('temperature offset read format error'), end
toff2 = fread(fid,ntemp,'real*8');
filemark2 = fread(fid, 1, 'integer*4');

% read the compressed absorptions
if gid > 1
  % do anything but water
  kcomp = zeros(nd,nlay,ntemp);
  filemark = 8 * ntemp * nlay;
  for i = 1:nd
    filemark2 = fread(fid, 1, 'integer*4');
    if filemark ~= filemark2, error('absorption read format error'), end
    kcomp(i,:,:) = fread(fid, [nlay,ntemp], 'real*8');
    filemark2 = fread(fid, 1, 'integer*4');
  end
else
  % do 5 partial pressures, for water
  kcomp = zeros(nd,nlay,ntemp,5);
  for pi = 1 : 5
    filemark = 8 * ntemp * nlay;
    for i = 1:nd
      filemark2 = fread(fid, 1, 'integer*4');
      if filemark ~= filemark2, error('absorption read format error'), end
      kcomp(i,:,:,pi) = fread(fid, [nlay,ntemp], 'real*8');
      filemark2 = fread(fid, 1, 'integer*4');
    end
  end
end

% read the basis set 
B = zeros(npts, nd);
filemark= 8 * npts;
for i = 1:nd
   filemark2 = fread(fid, 1, 'integer*4');
   if filemark ~= filemark2, error('basis read format error'), end
   B(:,i) = fread(fid, npts, 'real*8');
   filemark2 = fread(fid, 1, 'integer*4');
end

fclose(fid);

% save the matlab variables:
%
%   fr        1 x 10000           frequency scale
%   gid       1 x 1               gas ID
%   B     10000 x <d>             basis
%   kcomp   <d> x 100 x 11 x <p>  tabulated absorptions
%

fr = vchunk + (0:9999)*dfchunk/10000;

whos fr kcomp B
