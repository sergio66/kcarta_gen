function [data, wnums] = readkcUA(kfile)

% function [data, wnums] = readkcUA(kfile)
%
% readkc is a simple reader & unchunker for kcarta UA output files
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
flen = fread(fin, 1, 'integer*4');
iOut = fread(fin, 1, 'integer*4');
flen = fread(fin, 1, 'integer*4');

% read type (+1 for regular, -1 for Planck)
flen  = fread(fin, 1, 'integer*4');
iType = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

fprintf(1,'number of layers = %3i \n',nlayer);
fprintf(1,'start, stop freqs  of regular kCARTA = %12.6f %12.6f \n',fmin,fmax);
fprintf(1,'start, stop chunks of regular kCARTA = %12.6f %12.6f \n',lowchunk,highchunk);
fprintf(1,'regular point spacing, span per 10000 pts = %12.6f %12.6f \n',dv0,dblock0);
fprintf(1,'regular number of output = %4i \n',iOut);
fprintf(1,'***************************************\n');

% read actual blocks that kCARTA does NLTE for
% low and high chunk index
flen        = fread(fin, 1, 'integer*4');
lowchunk    = fread(fin, 1, 'integer*4');
highchunk   = fread(fin, 1, 'integer*4');
iNLTEChunks = fread(fin, 1, 'integer*4');
flen        = fread(fin, 1, 'integer*4');

% read actual start/stop freqs that kCARTA outputs/bloats
flen  = fread(fin, 1, 'real*4');
fmin1 = fread(fin, 1, 'real*4');
fmax1 = fread(fin, 1, 'real*4');
flen  = fread(fin, 1, 'real*4');

if (iType > 0)
  disp('reading in regular UA binary file');
elseif (iType < 0)
  disp('reading in planck UA binary file');
elseif (iType == 0) then
  error('oops unsupported type while in bloated binary file');
  end

fprintf(1,'***************************************\n');

%%%%%%%%%%%%%%%%%%%%% have read header info ... create arrays/matrices %%%%%
iChunks = highchunk - lowchunk + 1;
iChunks = iNLTEChunks;
fprintf(1,'Number of     kCarta chunks     = %4i \n',iChunks);
fprintf(1,'Number of outputs per chunk     = %4i \n',iOut);

ii=(1:10000*iChunks) - 1; 
wnums = fmin1 + dv0*ii; 

fprintf(1,'reading %3i levels from each chunk, of which 1st and Nth are where UA starts/ends \n',iOut)
 
data = zeros(10000*iChunks,iOut); 
for ii = 1 : iChunks 
  fprintf(1,'  reading in chunk number %4i .... \n',ii); 
  index = (1:10000) + (ii-1)*10000; 

  if iOut == 17
    % only dumping out rads

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
   
    %fprintf(1,'%3i %3i %3i \n',type,subtype,nrow);
    %fprintf(1,'%6i %8.2f %8.2f %8.6e \n',ncols,frlow,frhigh,frinc);

    for jj = 1 : iOut
      flen  = fread(fin, 1, 'integer*4'); 
      ra    = fread(fin, 10000, 'real*4'); 
      flen  = fread(fin, 1, 'integer*4'); 
      data(index,jj) = ra; 
    end 
  end
%{
  else
    % first dumping out ODs

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
   
    %fprintf(1,'%3i %3i %3i \n',type,subtype,nrow);
    %fprintf(1,'%6i %8.2f %8.2f %8.6e \n',ncols,frlow,frhigh,frinc);

    booX = iOut-2;
    for jj = 1 : iOut-2
      flen  = fread(fin, 1, 'integer*4');
      ra    = fread(fin, 10000, 'real*4'); 
      flen  = fread(fin, 1, 'integer*4'); 
      data(index,jj) = ra; 
      %if jj <= 10
      %  fprintf(1,'%3i %8.6e \n',jj,ra(1));
      %end
    end 
  
    % then dumping out rads

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
   
    %fprintf(1,'%3i %3i %3i \n',type,subtype,nrow);
    %fprintf(1,'%6i %8.2f %8.2f %8.6e \n',ncols,frlow,frhigh,frinc);

    for jj = 1 : 2
      flen  = fread(fin, 1, 'integer*4'); 
      ra    = fread(fin, 10000, 'real*4'); 
      flen  = fread(fin, 1, 'integer*4');
      data(index,jj+booX) = ra; 
    end 
  end
%}


end 
 
fclose(fin);

