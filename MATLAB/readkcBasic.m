function [data, wnums, caVersion] = readkcBasic(kfile)

% function [data, wnums, caVersion] = readkcBasic(kfile)
%
% readkcBasic is a simple reader & unchunker for kcarta output files
% and is also used for reading kCARTA clear sky column jacs 
%
% use readkcPCLSAM_coljac forPCLSAM 2slab column jacs; 
% that wrapper calls this
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
%   caVersion - descriptive  string set in kcarta.param at compile time,
%               CKD vers and nml comment (in nm_outout)
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

% version number
flen    = fread(fin, 1, 'integer*4');
version = fread(fin, 80, 'char');
caVersion.include_param = setstr(version');
version = caVersion.include_param;
flen    = fread(fin, 1, 'integer*4');

version_number=str2num(version(2:5));

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
iNumLayers=nlayer;              %<------------- my modification

% number of params 
flen    = fread(fin, 1, 'integer*4');
nparams = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');

flen    = fread(fin, 1, 'integer*4');
rparams = fread(fin, nparams, 'real*4');
flen    = fread(fin, 1, 'integer*4');
caVersion.ckd = rparams(2);

% comment
flen    = fread(fin, 1, 'integer*4');
if  (version_number >= 1.22)
  comment = fread(fin, 160, 'char');
elseif  (version_number >= 1.18)
  comment = fread(fin, 120, 'char');
else
  comment = fread(fin, 80, 'char');
end  
comment = setstr(comment');
flen    = fread(fin, 1, 'integer*4');
caVersion.comment = comment;

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

% read iLongOrShort
flen  = fread(fin, 1, 'integer*4');
htype = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if htype ~= 0
  fprintf(1,'this reader cannot work for kLongOrShort = -1,+1!! \n');
  error('Use readkcstd.m instead!!!');
end

% read freq step size (v(2)-v(1))
flen  = fread(fin, 1, 'real*4');
dv    = fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

% read block freq step size of each kCARTA 10000 pt chunk (v(10001) - v(1))
flen  = fread(fin, 1, 'real*4');
dblock= fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

% read total number of outputs per 10000 pt chunk
flen  = fread(fin, 1, 'integer*4');
iOut  = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

%%%%%%%%%%%%%%%%%%%%% have read header info ... create arrays/matrices %%%%%
iChunks = highchunk - lowchunk + 1;

fprintf(1,'Number of kCarta chunks     = %4i \n',iChunks);
fprintf(1,'Number of outputs per chunk = %4i \n',iOut);

ii=(1:10000*iChunks) - 1;
wnums = fmin + dv*ii;

data = zeros(10000*iChunks,iOut);
for ii = 1 : iChunks
  fprintf(1,'  reading in chunk number %4i .... \n',ii);
  index = (1:10000) + (ii-1)*10000;
  for jj = 1 : iOut
    flen  = fread(fin, 1, 'integer*4');
    ra    = fread(fin, 10000, 'real*4');
    flen  = fread(fin, 1, 'integer*4');
    data(index,jj) = ra;
    end
  end

fclose(fin);

dv  = 0.0025;
dvx = dv;
dvx = mean(diff(wnums));
if (max(diff(wnums)) - min(diff(wnums))) > 0.0001
  error('oops! cannot figure out point spacing!');
  end
ilen = 1 : length(wnums);
wnums = floor(wnums(1)) + (ilen-1)*dvx;
