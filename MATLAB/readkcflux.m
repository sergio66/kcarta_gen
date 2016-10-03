function [data, wnums, plevs, hgt, hr] = readkcflux(kfile, plevsIN, dfile)

% function [data, wnums, plevs, hgt, hr] = readkcflux(kfile, plevsIN, dfile)
%
% readkcflux is a simple reader & unchunker for kcarta flux output files
%
% INPUTS
%
%   kfile    - kcarta jacobian output file
%   plevsIN  - input pressure levels (should be same grid as where flux was computed) (optional)
%   dfile    - readkcflux output file (optional)
%
% OUTPUTS
%
%   data   - a w by n array of data from kcarta
%   wnums  - a w by 1 vector of data wavenumbers
%   plevs  - output pressure (mb) .. useful if plevsIN was not sent in
%   hgt    - output height   (km) .. 
%   hr    - optional heating rate, if the input file is *_ALL, which
%          - has up/down fluxes from which heating rate can be derived
%
% If the input parameter dfile is specified, then the data array
% is written to file dfile, and the return values [data, wnums]
% are left unassigned.
%
% The w by n data array has one row for each output wavenumber and
% one column for each "output data block" row.  In practice, this
% means the data array columns are a concatenation of whatever was 
% asked for in the kcarta driver file.
% 
% LIMITATIONS
% 

% Large data sets (e.g., mixed path sets of more than one or two
% "chunks") take up too much space to be returned in memory; for 
% such cases an output file should always be specified.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hr = [];

[fin,msg] = fopen(kfile, 'r');
if fin == -1
  error(['error opening input file\n', msg]);
  end
fid=fin;                    %<------------- my modification
%%%%%%%%%%%%%%%%%%%%%%
% READ HEADER RECORDS
%%%%%%%%%%%%%%%%%%%%%%

% MAIN HEADER

% comment
flen    = fread(fin, 1, 'integer*4');
comment = fread(fin, 120, 'char');
comment = setstr(comment');
flen    = fread(fin, 1, 'integer*4');

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
iNumLayers = nlayer;              %<------------- my modification

% start, stop frequency
flen = fread(fin, 1, 'integer*4');
fmin = fread(fin, 1, 'real*4');
fmax = fread(fin, 1, 'real*4');
flen = fread(fin, 1, 'integer*4');

% low and high chunk index
flen      = fread(fin, 1, 'integer*4');
lowchunk  = fread(fin, 1, 'integer*4');
highchunk = fread(fin, 1, 'integer*4');
flen      = fread(fin, 1, 'integer*4');

% INFO ABOUT how many types of fluxes (only one : net=up-down)
% number of fluxes
flen   = fread(fin, 1, 'integer*4');
nimportant = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

%number of atmospheres
flen   = fread(fin, 1, 'integer*4');
natmos = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

%number of layers per atmosphere
flen   = fread(fin, 1, 'integer*4');
numlay = fread(fin, natmos, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

%%%%%%%%%%%%%%%%%%%%%%
% DATA RECORD SUMMARY
%%%%%%%%%%%%%%%%%%%%%%

nchunk = highchunk-lowchunk+1;         %number of 10000 pt chunks

flen       = fread(fin, 1, 'integer*4');
iTotal     = fread(fin, 1, 'integer*4'); % total number of slabs
iNumAtm    = fread(fin, 1, 'integer*4'); % number of atmospheres
iImportant = fread(fin, 1, 'integer*4'); % number of fluxes (==1)
flen       = fread(fin, 1, 'integer*4');

flen1       = fread(fin, 1, 'integer*4');
iaNumLayers = fread(fin, iNumAtm, 'integer*4'); % rows in each ODB
flen2       = fread(fin, 1, 'integer*4');
nODBrows    = iaNumLayers;

% sanity checks
if flen1 ~= flen2
  error('Fortran records out of phase!');
end

if nchunk ~= iTotal
  error('readkcflux: chunk counts do not match!');
end

if iNumAtm ~= natmos
  error('readkcflux: number of atmos do not match!');
end


nODBs=natmos;                        %fluxes computed for each atmos
fprintf(2, 'readkcflux: %d chunks, %d ODBs(num of atmos), %d total rows\n', ...
        nchunk, nODBs, sum(nODBrows));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output file or array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noutrow = nchunk * 10000;       % number of output rows
noutcol = sum(nODBrows);        % number of output columns
ODBrowdat = zeros(10000, 1);    % ODB row data buffer

if nargin == 1 | nargin == 2
 
  % initialize output arrays
  data = zeros(noutrow, noutcol);
  wnums = zeros(noutrow, 1);

else

  % pre-extend the output file
  totbytes = noutrow * noutcol * 4;
  fprintf(2, 'readkcflux: creating %d x %d element output array\n', ...
          noutrow, noutcol);

  [fout,msg] = fopen(dfile, 'w');
  if fout == -1
    error(['error creating output file', msg]);
  end
  z = zeros(noutrow,1);
  for i = 1:noutcol
    if fwrite(fout, z, 'float') ~= noutrow
      error(['fwrite in extend failed, i=',num2str(i)]);
    end
  end
  fclose(fout);
  clear z

  % open the extended file in "update" mode
  [fout,msg] = fopen(dfile, 'r+');
  if fout == -1
    error(['error opening output file', msg]);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read "output data blocks"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for chunk = 1:nchunk
  cumODBrow = 1;
  for atmospheres=1:iNumAtm   %%%%%go atm by atm

    % read ODB header
    flen    = fread(fin, 1, 'integer*4');
    type    = fread(fin, 1, 'integer*4'); % ODB type
    subtype = fread(fin, 1, 'integer*4'); % ODB subtype
    nrow    = fread(fin, 1, 'integer*4'); % number of rows, this ODB
    flen    = fread(fin, 1, 'integer*4');

    flen    = fread(fin, 1, 'integer*4');
    ncols   = fread(fin, 1, 'integer*4'); % spectral pts in a chunk
    frlow   = fread(fin, 1, 'real*4');    % low freq of chunk
    frhigh  = fread(fin, 1, 'real*4');    % high freq of chunk
    frinc   = fread(fin, 1, 'real*4');    % frequency increment
    flen    = fread(fin, 1, 'integer*4');

    if chunk == 1
      fprintf(2, 'readkcflux: ODB atm# = %d, subtype = %d, rows = %d\n', ...
              type, subtype, nrow);
    end

    ODBrowmax=numlay(atmospheres);

    % sanity check 
    if nrow ~= ODBrowmax
      fprintf(2, 'readkcflux: WARNING -- bad ODB row count,  ');
      fprintf(2, 'header = %-4d ODB = %-4d\n', numlay(atmospheres), nrow);
      end

    %read in data
    for ODBrow = 1:ODBrowmax
      flen1 = fread(fin, 1, 'integer*4');
      [ODBrowdat, count] = fread(fin, [10000,1], 'real*4');
      if count ~= 10000
        error(['fread failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
      end
      flen2 = fread(fin, 1, 'integer*4');
       if nargin == 1 | nargin == 2
        % write ODBrowdat to the data array
        outrows = (chunk-1) * 10000 + 1 : chunk * 10000;
        outcol = cumODBrow;
        data(outrows, outcol) = ODBrowdat;
        wnums(outrows) = frlow + (0:9999)*frinc;
      else
        % write ODBrowdat to the output file
        outpos = ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        if fseek(fout, outpos, -1) ~= 0
          error(['fseek failed,odb=',num2str(odb),' chunk=',num2str(chunk)]);
          end
        if fwrite(fout, ODBrowdat, 'real*4') ~= 10000
          error(['fwrite failed,odb=',num2str(odb),'chunk=',num2str(chunk)]);
          end     
        end
      cumODBrow = cumODBrow + 1;
      end % loop on current ODB rows
    end % loop on atmospheres
  end % loop on chunks

fclose(fin);

if nargin == 3 
  fclose (fout);
end

if iNumAtm > 1
  datax = data;
  clear data
  for aa = 1 : iNumAtm
    len = nrow;
    inds = (aa-1)*len + (1:len);
    data (aa,:,:) = datax(:,inds);
  end
  clear datax
  [aa,mm,nn] = size(data);
else
  [mm,nn] = size(data);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataX = data;

if nargin == 1
  plevs = load('airslevels.dat');  %% from GND to TOA, decreasing p; 101 levels
  if length(strfind(kfile,'_ALL')) > 0    
    plevs = plevs(101-nn/2+1:end);
  else
    plevs = plevs(101-nn+1:end);
  end
else
  plevs = plevsIN;
  if plevs(20) < plevs(30)
    plevs = flip(plevs);
  end
  if length(strfind(kfile,'_ALL')) > 0      
    plevs = plevs(length(plevs)-nn/2+1:end);
  else
    plevs = plevs(length(plevs)-nn+1:end);    
  end
end
hgt   = p2h(plevs)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ilen = 1 : length(wnums);
dv  = 0.0025;
dvx = dv;
dvx = mean(diff(wnums));
if (max(diff(wnums)) - min(diff(wnums))) > 0.0001
  error('oops! cannot figure out point spacing!');
  end

iPower = -9:+3;
remain = dvx ./ (10.^iPower);
boo = find(remain < 10,1);
dvxnew = round(10*remain(boo))/10 * 10^(iPower(boo));

if (abs(dvxnew/dvx) - 1) < 1e-5
  wnums = wnums(1) + (ilen-1)*dvxnew;
else
  disp('warning ... in readkcflux.m ... could not figure out wnum spacing!!!')
  wnums = floor(wnums(1)) + (ilen-1)*dvx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for aa = 1 : iNumAtm
  if iNumAtm > 1
    data = squeeze(dataX(aa,:,:));
  else
    data = dataX;
  end
  plot_fluxes
  plot_heatingrates
  if iNumAtm > 1
    fprintf(1,'showing heating rate/fluxes for atmosphere %2i of %2i \n',aa,iNumAtm);
    disp('ret'); pause
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = dataX;