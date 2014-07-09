function [data, wnums, detail] = readkcstd_detail(kfile, dfile)

% function [data, wnums] = readkcstd(kfile, dfile)
%
% readkc is a simple reader & unchunker for kcarta output files
%
% INPUTS
%
%   kfile  - kcarta output file
%   dfile  - readkc output file (optional)
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

% this is the same as /asl/matlab/aslutil/readkc.m except it can read long 
% headers as well

% H. Motteler, 8/24/98
% Updated for version 1.03 by Scott Hannon, 30 July 1999
% Updated to read in "long version" kcarta files by Sergio Machado
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
version = setstr(version');
flen    = fread(fin, 1, 'integer*4');

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

% comment
flen    = fread(fin, 1, 'integer*4');
comment = fread(fin, 80, 'char');
comment = setstr(comment');
flen    = fread(fin, 1, 'integer*4');

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

if htype == 0
  fprintf(1,'this reader cannot work for kLongOrShort = 0!! \n');
  error('Use readkcBasic.m instead!!!');
  end

%%%%
% This block added by Scott for version 1.03+
version_number=str2num(version(2:5));
if ((version_number >= 1.03) & (version_number <= 1.05)) 
  % read in M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
  flen = fread(fin,1,'integer*4');
  iaS  = fread(fin,6,'integer*4');
  flen = fread(fin,1,'integer*4');

  % read in pressure levels
  flen = fread(fin,1,'integer*4');
  detail.raP  = fread(fin,nlayer+1,'real*4');
  flen = fread(fin,1,'integer*4');

% Larry M. found bug in pProf writing out one too many elements
elseif ((version_number >= 1.06) & (version_number <= 1.08))
  % read in M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
  flen = fread(fin,1,'integer*4');
  iaS  = fread(fin,6,'integer*4');
  flen = fread(fin,1,'integer*4');

  % read in pressure levels
  flen = fread(fin,1,'integer*4');
  raP  = fread(fin,nlayer,'real*4');
  flen = fread(fin,1,'integer*4');

elseif ((version_number >= 1.09))
  % read in pressure levels
  flen = fread(fin,1,'integer*4');
  raP  = fread(fin,nlayer,'real*4');
  flen = fread(fin,1,'integer*4');
  end

%%%%%% GAS PATH HEADER

if htype < 0                 %%%%%%%%%%%%%%%%%%% short version 
  % number of paths to be output
  flen    = fread(fin, 1, 'integer*4');
  ngasout = fread(fin, 1, 'integer*4');
  flen    = fread(fin, 1, 'integer*4');

  if (ngasout > 0)
    flen     = fread(fin, 1, 'integer*4');
    gaspaths = fread(fin, ngasout, 'integer*4');
    flen     = fread(fin, 1, 'integer*4');
    end
  end

if htype > 0     
  %%%%%%%% long version ..... copied direct from readgaspaths.m

  %read iNumPaths=iNumGases*iNumLayers 
  flen       = fread(fid,1,'integer*4'); 
  iNumPaths  = fread(fid,1,'integer*4'); 
  flen       = fread(fid,1,'integer*4'); 
  ngasout    = iNumPaths;        %<---------------------- my modification
   
  iNumGases = iNumPaths/iNumLayers;   %this is number of gases 
  iaGasID   = zeros(1,iNumGases); 
  raaAmt    = zeros(iNumLayers,iNumGases); 
  raaTemp   = zeros(iNumLayers,iNumGases); 
 
  %now read the PathNum GasID Temperature Amount for each path in the profile 
  for ii=1:iNumGases 
    for jj=1:iNumLayers 
      flen   = fread(fid,1,'integer*4'); 
      iPath  = fread(fid,1,'integer*4'); 
      iGasID = fread(fid,1,'integer*4'); 
      rTemp  = fread(fid,1,'real*4'); 
      rAmt   = fread(fid,1,'real*4'); 
      flen   = fread(fid,1,'integer*4'); 
      raaAmt(jj,ii)  = rAmt; 
      raaTemp(jj,ii) = rTemp; 
      end  
    iaGasID(ii)=iGasID; 
    end 
  detail.raaAmt    = raaAmt;
  detail.raaTemp   = raaTemp;
  detail.iaGasID   = iaGasID;

  %read in Number of Paths to be output each time 
  flen         = fread(fid,1,'integer*4'); 
  iNumPathsOut = fread(fid,1,'integer*4'); 
  flen         = fread(fid,1,'integer*4'); 
  if ((iNumPathsOut < 0) | (iNumPathsOut > iNumPaths)) 
    error('Mistake in number of paths to be output!!!') 
    end 
  ia=0; 
  clear ia; 
  if (iNumPathsOut > 0) 
    %read in paths 
    iaOutPaths = zeros(1,iNumPathsOut); 
    flen       = fread(fid,1,'integer*4'); 
    ia         = fread(fid,iNumPathsOut,'integer*4'); 
    iaOutPaths(1,:)=ia'; 
    clear ia 
    flen       = fread(fid,1,'integer*4'); 
    gaspaths   = iaOutPaths;       %<------------ my modification
    end 
  end        %if htype > 0 

% MIXED PATH HEADER
if htype < 0                    %%%%%%%%%%%%%%%%%%%%short version
  flen      = fread(fin, 1, 'integer*4');
  nmixpaths = fread(fin, 1, 'integer*4');
  flen      = fread(fin, 1, 'integer*4');

  if (nmixpaths > 0)

    flen    = fread(fin, 1, 'integer*4');
    nmixout = fread(fin, 1, 'integer*4');
    flen    = fread(fin, 1, 'integer*4');

    if (nmixout > 0)
      flen     = fread(fin, 1, 'integer*4');
      mixpaths = fread(fin, nmixout, 'integer*4');
      flen     = fread(fin, 1, 'integer*4');
      end
    end
  end

if htype > 0     
  %%%% long version ... copied direct from readmixedpaths.m
  %read in number of mixed paths in *MIXFIL 
  flen      = fread(fid,1,'integer*4'); 
  iNpmix    = fread(fid,1,'integer*4'); 
  flen      = fread(fid,1,'integer*4'); 
  nmixpaths = iNpmix;       %<------------ my modification

  %read in #of lines that have info in *MIXFIL 
  flen           = fread(fid,1,'integer*4'); 
  iMixFileLines  = fread(fid,1,'integer*4'); 
  flen           =fread(fid,1,'integer*4'); 
 
  if (iNpmix > 0) 
    caaMixPathInfo = zeros(iMixFileLines,130); 
    %read in the mixing table 
    for ii=1:iMixFileLines 
      flen           = fread(fid,1,'integer*4'); 
      caMixPathInfo  = fread(fid,130,'char'); 
      flen           = fread(fid,1,'integer*4'); 
      end 
 
    %read in mixed path temperatures 
    raMixTemp = zeros(1,iNpmix); 
    flen      = fread(fid,1,'integer*4'); 
    ra        = fread(fid,iNpmix,'real*4'); 
    flen      = fread(fid,1,'integer*4'); 
    raMixTemp(1,:)=ra(1:length(ra))'; 
    clear ra  
 
    %read in Number of mixed paths to be output each time 
    flen            = fread(fid,1,'integer*4'); 
    iNumMixPathsOut = fread(fid,1,'integer*4'); 
    flen            = fread(fid,1,'integer*4'); 
    nmixout         = iNumMixPathsOut;       %<------------ my modification

    if ((iNumMixPathsOut < 0) | (iNumMixPathsOut > iNpmix)) 
      error('Mistake in number of mixed paths to be output!!!') 
      end 
    iaOutMixPaths=0; 
    clear iaOutMixPaths; 
    if (iNumMixPathsOut > 0) 
      %read in paths 
      iaOutMixPaths = zeros(1,iNumMixPathsOut); 
      flen          = fread(fid,1,'integer*4'); 
      ia            = fread(fid,iNumMixPathsOut,'integer*4'); 
      flen          = fread(fid,1,'integer*4'); 
      %fprintf(1,'%3i\n',iaOutMixPaths); 
      iaOutMixPaths(1,:)=ia'; 
      clear ia   
      mixpaths      =iaOutMixPaths;       %<------------ my modification
      end 
    end 
  end                %if LS > 0 

% ATMOSPHERE HEADER

% number of atmospheres in *RADFIL
flen   = fread(fin, 1, 'integer*4');
natmos = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% max number of emissivity points per atm
flen  = fread(fin, 1, 'integer*4');
nemis = fread(fin, 1, 'integer*4');
flen  = fread(fin, 1, 'integer*4');

if (natmos > 0)

  for i = 1:natmos

    % read in atmosphere number, # of mixed paths in atmosphere
    flen = fread(fin, 1, 'integer*4');
    iAtm = fread(fin, 1, 'integer*4');
    iNumPathsAtmos = fread(fin, 1, 'integer*4');
    flen = fread(fin, 1, 'integer*4');

    if (iNumPathsAtmos > 0)

      %read in paths making up atmosphere
      flen = fread(fin, 1, 'integer*4');
      ia   = fread(fin, iNumPathsAtmos, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      flen   = fread(fin, 1, 'integer*4');
      detail.rTbdy(i)  = fread(fin, 1, 'real*4');
      detail.rTinit(i) = fread(fin, 1, 'real*4');
      detail.rSat(i)   = fread(fin, 1, 'real*4');
      detail.rHgt(i)   = fread(fin, 1, 'real*4');
      flen(i)   = fread(fin, 1, 'integer*4');

       
      flen           = fread(fin, 1, 'integer*4');
      detail.ikSolar(i)        = fread(fin, 1, 'integer*4');
      detail.rkSolarAngle(i)   = fread(fin, 1, 'real*4');
      detail.rkSolarRefl(i)    = fread(fin, 1, 'real*4');
      detail.ikThermal(i)      = fread(fin, 1, 'integer*4');
      detail.rkThermalAngle(i) = fread(fin, 1, 'real*4');
      detail.ikThermalJacob(i) = fread(fin, 1, 'integer*4');
      flen           = fread(fin, 1, 'integer*4');

      flen = fread(fin, 1, 'integer*4');
      iEms = fread(fin, 1, 'integer*4');
      flen = fread(fin, 1, 'integer*4');

      %now read the freq, ems data points
      for jj=1:iEms
        flen = fread(fin, 1, 'integer*4');
        rff  = fread(fin, 1, 'real*4');
        ree  = fread(fin, 1, 'real*4');
        detail.efreq(i,jj) = rff;
        detail.emis(i,jj)  = ree;
        flen = fread(fin, 1, 'integer*4');
        end 
     
      %now read in layers to be output
      flen          = fread(fin, 1, 'integer*4');
      iNumLayersOut = fread(fin, 1, 'integer*4');
      flen          = fread(fin, 1, 'integer*4');

      if (iNumLayersOut > 0)
        % read in mixed paths to be output
        flen = fread(fin, 1, 'integer*4');
        ia   = fread(fin, iNumLayersOut, 'integer*4');
        flen = fread(fin, 1, 'integer*4');
        end

      if (iNumLayersOut > 0)
        % read in pressures at which radiances to be output
        flen = fread(fin, 1, 'integer*4');
        ra   = fread(fin, iNumLayersOut, 'real*4');
        flen = fread(fin, 1, 'integer*4');
        end
      end
    end
  end

%%%%%%%%%%%%%%%%%%%%%%
% DATA RECORD SUMMARY
%%%%%%%%%%%%%%%%%%%%%%

flen   = fread(fin, 1, 'integer*4');
nchunk = fread(fin, 1, 'integer*4'); % total number of chunks
nODBs  = fread(fin, 1, 'integer*4'); % number of (ODBs) "output data blocks"
flen   = fread(fin, 1, 'integer*4');

flen1    = fread(fin, 1, 'integer*4');
nODBrows = fread(fin, nODBs, 'integer*4'); % rows in each ODB
flen2    = fread(fin, 1, 'integer*4');

% sanity checks
if flen1 ~= flen2
  error('Fortran records out of phase!');
end

if nchunk ~= 1+highchunk-lowchunk
  error('readkc: chunk counts do not match!');
end

fprintf(2, 'readkc: %d chunks, %d ODBs, %d total rows\n', ...
        nchunk, nODBs, sum(nODBrows));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize output file or array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noutrow = nchunk * 10000;       % number of output rows
noutcol = sum(nODBrows);        % number of output columns
ODBrowdat = zeros(10000, 1);    % ODB row data buffer

if nargin == 1
 
  % initialize output arrays
  data = zeros(noutrow, noutcol);
  wnums = zeros(noutrow, 1);

else

  % pre-extend the output file
  totbytes = noutrow * noutcol * 4;
  fprintf(2, 'readkc: creating %d x %d element output array\n', ...
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

  for odb = 1:nODBs

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
      fprintf(2, 'readkc: ODB type = %d, subtype = %d, rows = %d\n', ...
                type, subtype, nrow);
    end

    % sanity check
    if nrow ~= nODBrows(odb)
      fprintf(2, 'readkc: WARNING -- bad ODB row count,  ');
      fprintf(2, 'header = %-4d ODB = %-4d\n', nODBrows(odb), nrow);
    end

    % read a row of ODB data
    for ODBrow = 1:nODBrows(odb)
      flen1 = fread(fin, 1, 'integer*4');
      [ODBrowdat, count] = fread(fin, [10000,1], 'real*4');
      if count ~= 10000
        error(['fread failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
      end
      flen2 = fread(fin, 1, 'integer*4');

      if nargin == 1

        % write ODBrowdat to the data array
        outrows = (chunk-1) * 10000 + 1 : chunk * 10000;
        outcol = cumODBrow;
        data(outrows, outcol) = ODBrowdat;
        wnums(outrows) = frlow + (0:9999)*frinc;

      else

        % write ODBrowdat to the output file
        outpos = ((cumODBrow-1)*nchunk*10000 + (chunk-1)*10000)*4;
        if fseek(fout, outpos, -1) ~= 0
          error(['fseek failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end
        if fwrite(fout, ODBrowdat, 'real*4') ~= 10000
          error(['fwrite failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
        end     

      end

      cumODBrow = cumODBrow + 1;

    end % loop on current ODB rows

  end % loop on ODBs

end % loop on chunks

fclose(fin);

if nargin == 2 
  fclose (fout);
end

dv  = 0.0025;
dvx = dv;
dvx = mean(diff(wnums));
if (max(diff(wnums)) - min(diff(wnums))) > 0.0001
  error('oops! cannot figure out point spacing!');
  end
ilen = 1 : length(wnums);
wnums = floor(wnums(1)) + (ilen-1)*dvx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
detail.version = version;
detail.iNumLayers = iNumLayers;
detail.nparams = nparams;
detail.rparams = rparams;
detail.comment = comment;
detail.fmin = fmin;
detail.fmax = fmax;
detail.lowchunk  = lowchunk;
detail.highchunk = highchunk;
detail.htype     = htype;
