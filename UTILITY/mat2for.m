function mat2for(gid, vchunk, cdir, fdir, dtype)

% function mat2for(gid, vchunk, cdir, fdir, dtype)
%
% produce old-format fortran data for kcarta
%
% inputs
%   gid    -  gas ID
%   vchunk -  wavenumber start of 25 1/cm chunk
%   cdir   -  directory for matlab mat source files
%   fdir   -  directory for fortran output files
%   dtype  -  output data type (for matlab open)
%

% default output data type
if nargin < 5
 dtype = 'ieee-be';
end

% default dir's for fortran compressed data 
if nargin < 4
  if gid == 1
    fdir = '/home/motteler/abstab/fbin/h2o.ieee-be';
  else
    fdir = '/home/motteler/abstab/fbin/etc.ieee-be';
  end
end

% default dir's for matlab compressed data
if nargin < 3
 if gid == 1
   cdir = '/home/motteler/abstab/kcomp.h2o'; 
 elseif 1 < gid  & gid <= 50
   cdir = '/home/motteler/abstab/kcomp';
 else
   cdir = '/home/motteler/abstab/kcomp.xsec';
 end
end

% load the matlab variables:
%
%   fr        1 x 10000           frequency scale
%   gid       1 x 1               gas ID
%   B     10000 x <d>             basis
%   kcomp   <d> x 100 x 11 x <p>  tabulated absorptions
%
eval(sprintf('load %s/cg%dv%d.mat', cdir, gid, vchunk));

% 
% variable names from old load
% 
%   GAS_NUM            1x1              8  double array
%   MAXU               1x1              8  double array
%   START_FREQ         1x1              8  double array
%   bt_ref             3x10000     240000  double array
%   k_comp            30x1100      264000  double array
%   k_max_errors       1x3             24  double array
%   k_std_errors       1x3             24  double array
%   u              10000x30       2400000  double array

fname = sprintf('r%d_g%d.dat', vchunk, gid);

ktype=2;

npts=10000;
fstep=0.0025;

nlay = 100;
[npts,nd]=size(B);

% temperature offsets
ntemp=11;
toff=-50.0:10.0:50.0;

zaza = [fdir,'/',fname];
ee = exist(zaza,'file');
iDo = +1;
if ee ~= 0
  fprintf(1,'%s already exists \n',zaza);
  disp('try again!')
  iDo = -1;
  end

if iDo > 0
  fprintf(1,'writing to %s \n',zaza);

  fid=fopen([fdir,'/',fname], 'w' ,dtype);
  if fid == -1
    error('can''t open output file')
  end

  % Write header info
  filemark= 4 + 8 + 8 + 4 + 4 + 4 + 4 + 4 + 4 + 4 + 4;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,gid,'integer*4');
  fwrite(fid,[vchunk,fstep],'real*8');
  fwrite(fid,[npts,nlay,ktype,nd,ntemp,nlay,npts,nd],'integer*4');
  fwrite(fid,filemark,'integer*4');

  % Write the temperature offsets
  filemark= 8 * ntemp;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,toff,'real*8');
  fwrite(fid,filemark,'integer*4');

  % write the compressed absorptions
  if gid > 1
    % do anything but water
    filemark= 8 * ntemp * nlay;
    for i = 1:nd
      fwrite(fid,filemark,'integer*4');
      fwrite(fid, kcomp(i,:), 'real*8');
      fwrite(fid,filemark,'integer*4');
    end
  else
    % do 5 partial pressures, for water
    for pi = 1 : 5
      filemark= 8 * ntemp * nlay;
      for i = 1:nd
        fwrite(fid,filemark,'integer*4');
        fwrite(fid, kcomp(i,:,:,pi), 'real*8');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end

  % write the basis set 
  filemark= 8 * npts;
  for i = 1:nd
    fwrite(fid,filemark,'integer*4');
    fwrite(fid,B(:,i),'real*8');
    fwrite(fid,filemark,'integer*4');
  end

  fclose(fid);
  end

