function [data, wnums] = readkcBasic(kfile)

% function [data, wnums] = readkcBasic(kfile)
%
% readkc is a simple reader & unchunker for kcarta output files
%
% INPUTS
%
%   kfile  - kcarta output file
%   dfile  - readkc output file (optional)   ..... NOT YET IMPLEMENTED
%
% OUTPUTS
%
%   data   - a w by n array of data from kcarta
%   wnums  - a w by 1 vector of data wavenumbers
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

[fin,msg] = fopen(kfile, 'r');
if fin == -1
  error(['error opening input file\n', msg]);
  end
fid=fin;                    %<------------- my modification
%%%%%%%%%%%%%%%%%%%%%%
% READ HEADER RECORDS
%%%%%%%%%%%%%%%%%%%%%%

% MAIN HEADER

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
kProfLayer=nlayer;              %<------------- my modification

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

% read regular freq step size (v(2)-v(1))
flen  = fread(fin, 1, 'real*4');
dv0   = fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

% read block freq step size of each regular res kCARTA 10000 pt chunk 
% (v(10001) - v(1))
flen    = fread(fin, 1, 'real*4');
dblock0 = fread(fin, 1, 'real*4');
flen   = fread(fin, 1, 'real*4');

% read total number of iTotal outputs per regular 10000 pt chunk
flen  = fread(fin, 1, 'integer*4');
iOut0 = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

fprintf(1,'number of layers = %3i \n',nlayer);
fprintf(1,'start, stop freqs  of regular kCARTA = %12.6f %12.6f \n',fmin,fmax);
fprintf(1,'start, stop chunks of regular kCARTA = %12.6f %12.6f \n',lowchunk,highchunk);
fprintf(1,'regular point spacing, span per 10000 pts = %12.6f %12.6f \n',dv0,dblock0);
fprintf(1,'regular number of output = %4i \n',iOut0);
fprintf(1,'***************************************\n');

% read iType
flen  = fread(fin, 1, 'integer*4');
iType = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

% read iNumberofLayers in atm 1
flen        = fread(fin, 1, 'integer*4');
iNumLayers1 = fread(fin, 1, 'integer*4');
flen        = fread(fin, 1, 'integer*4');

% read boxcar points 
flen  = fread(fin, 1, 'integer*4');
iBox  = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

% read fine freq step size (v(2)-v(1))
flen  = fread(fin, 1, 'integer*4');
dv    = fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'integer*4');

% read actual blocks that kCARTA does NLTE for
% low and high chunk index
flen         = fread(fin, 1, 'integer*4');
lowchunk1    = fread(fin, 1, 'integer*4');
highchunk1   = fread(fin, 1, 'integer*4');
iNLTEChunks1 = fread(fin, 1, 'integer*4');
flen         = fread(fin, 1, 'integer*4');

% read actual start/stop freqs that kCARTA bloats
flen  = fread(fin, 1, 'real*4');
fmin1 = fread(fin, 1, 'real*4');
fmax1 = fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

% read block freq step size of each high res kCARTA 10000 pt chunk 
flen  = fread(fin, 1, 'real*4');
dblock= fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

% read total number of iTotal outputs per fine res 50000 -> 5*10000 pt chunk
flen  = fread(fin, 1, 'integer*4');
iOut  = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if (iType > 0)
  disp('reading in bloated regular binary file');
elseif (iType < 0)
  disp('reading in bloated planck binary file');
elseif (iType == 0) then
  error('oops unsupported type while in bloated binary file');
  end

fprintf(1,'number of boxcar pts = %3i \n',iBox);
fprintf(1,'start, stop freq of bloated kCARTA = %12.6f %12.6f \n',fmin1,fmax1);
fprintf(1,'start, stop, number chunks of bloated kCARTA = %4i %4i %4i \n',lowchunk1,highchunk1,iNLTEChunks1);
fprintf(1,'bloated point spacing, span per 50000 pts = %12.6f %12.6f \n',dv,dblock);
if (iType > 0)
  fprintf(1,'regular*nbox number of output = %4i \n',iOut);
  iOut0 = iOut0;   %%unchanged!
elseif (iType < 0)
  fprintf(1,'regular*nbox number of output = %4i \n',iNumLayers1*iBox);
  iOut0 = iNumLayers1;   %%changed!
elseif (iType == 0) 
  error('oops unsupported type while in bloated binary file');
  end

fprintf(1,'***************************************\n');

%%%%%%%%%%%%%%%%%%%%% have read header info ... create arrays/matrices %%%%%
iChunks = highchunk - lowchunk + 1;
iChunks = iNLTEChunks1;
fprintf(1,'Number of bloated kCarta chunks     = %4i \n',iChunks);
fprintf(1,'Number of outputs per chunk         = %4i \n',iOut0);

fmin1 = fmin1 - dv*floor(iBox/2);
ii=(1:10000*iChunks*iBox) - 1;
wnums = fmin1 + dv*ii;

data = zeros(10000*iChunks*iBox,iOut0);
for iC = 1 : iChunks
  fprintf(1,'  reading in chunk number %4i .... \n',iC);

  % number of params 
  flen    = fread(fin, 1, 'integer*4');
  kMaxPts = fread(fin, 1, 'integer*4');
  f1      = fread(fin, 1, 'real*4');
  f2      = fread(fin, 1, 'real*4');
  df      = fread(fin, 1, 'real*4');
  flen    = fread(fin, 1, 'integer*4');

  flen       = fread(fin, 1, 'integer*4');
  kBloatPts  = fread(fin, 1, 'integer*4');
  df         = fread(fin, 1, 'real*8');
  kBoxCarUse = fread(fin, 1, 'integer*4');
  flen       = fread(fin, 1, 'integer*4');

  flen    = fread(fin, 1, 'integer*4');
  daStart  = fread(fin, kBoxCarUse, 'real*8');
  flen    = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');
  daEnd    = fread(fin, kBoxCarUse, 'real*8');
  flen    = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');
  daStartBox  = fread(fin, kBoxCarUse, 'real*8');
  flen    = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');
  daEndBox  = fread(fin, kBoxCarUse, 'real*8');
  flen    = fread(fin, 1, 'integer*4');
  
  for jj = 1 : iOut0
    for iB = 1 : iBox
      flen  = fread(fin, 1, 'integer*4');
      ra    = fread(fin, 10000, 'real*4');
      flen  = fread(fin, 1, 'integer*4');
      index = (1:10000) + (iB-1)*10000;
      index = index + (iC-1)*10000*iBox;
      data(index,jj) = ra;
      end
    end
  end

fclose(fin);

