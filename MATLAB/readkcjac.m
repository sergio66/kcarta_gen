function [data, wnums, natmos, numlay, ngases] = readkcjac(kfile, dfile)

% function [data, wnums, natmos, numlay, ngases] = readkcjac(kfile, dfile)
%
% readkcjac is a simple reader & unchunker for kcarta jacobian output files
%
% INPUTS
%
%   kfile  - kcarta jacobian output file
%   dfile  - readkcjac output file (optional)
%
% OUTPUTS
%
%   data   - a w by n array of data from kcarta
%   wnums  - a w by 1 vector of data wavenumbers
%   natmos - integer telling you how many atmospheres are in the file
%   numlay - integer array telling you how many layers each atmosphere has
%   ngases - integer telling you how many gases there are

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
iNumLayers=nlayer;              %<------------- my modification

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

% INFO ABOUT GASES WHOSE JACOBIAN IS COMPUTED
% number of gases
flen   = fread(fin, 1, 'integer*4');
ngases = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% number of atm
flen   = fread(fin, 1, 'integer*4');
natmos = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

%number of layers per atmosphere
flen   = fread(fin, 1, 'integer*4');
numlay = fread(fin, natmos, 'integer*4');
flen   = fread(fin, 1, 'integer*4');

% PROFILE INFO ABOUT GASES WHOSE JACOBIAN IS COMPUTED
iNumGases=ngases;
%now read the PathNum GasID Temperature Amount for each path in the profile 
for ii=1:iNumGases 
  for jj=1:iNumLayers 
    fmjunk1=fread(fid,1,'integer*4'); 
    iGasID=fread(fid,1,'integer*4');
    rTemp=fread(fid,1,'real*4'); 
    rAmt=fread(fid,1,'real*4'); 
    fmjunk1=fread(fid,1,'integer*4'); 
    raaAmt(jj,ii)=rAmt; 
    raaTemp(jj,ii)=rTemp; 
    end  
  iaGasID(ii)=iGasID; 
  end 
 
%%%%%%%%%%%%%%%%%%%%%%
% DATA RECORD SUMMARY
%%%%%%%%%%%%%%%%%%%%%%

nchunk=highchunk-lowchunk+1;         %number of 10000 pt chunks

flen   = fread(fin, 1, 'integer*4');
iTotal = fread(fin, 1, 'integer*4'); % total number of chunks
iNumAtm= fread(fin, 1, 'integer*4'); % number of atmospheres
iGasDQ = fread(fin, 1, 'integer*4'); % number of gases whose d/dq is done
flen   = fread(fin, 1, 'integer*4');

flen1       = fread(fin, 1, 'integer*4');
iaNumLayers = fread(fin, iNumAtm, 'integer*4'); % rows in each ODB
flen2       = fread(fin, 1, 'integer*4');

% sanity checks
if flen1 ~= flen2
  error('Fortran records out of phase!');
end

if nchunk ~= iTotal;
  error('readkcjac: chunk counts do not match!');
end

if ngases ~= iGasDQ;
  error('readkcjac: chunk counts do not match!');
end

blah  = ngases+1+1;                  %ngases d/dq and d/dT,weight fcns
for jj=1:natmos       
  nODBrows(jj)=blah*numlay(jj)+4; 
  end

nODBs=natmos;                        %jacobians computed for each atmos
nOBDs=natmos*(blah+1);               %in blah, also have 4 d/d(surface params)
fprintf(2, 'readkcjac: %d chunks, %d ODBs(num of atmos), %d total rows\n', ...
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
  fprintf(2, 'readkcjac: creating %d x %d element output array\n', ...
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

    for jacobItem = 1:iGasDQ+3    %%%%read the gas d/dqs + d/Dt + wgt + surface
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
        fprintf(2, 'readkcjac: ODB atm# = %d, subtype = %d, rows = %d\n', ...
                type, subtype, nrow);
      end

      %%%%%% remember for gas d/dq, or temp d/dT, or wgt fcn we have 
      %%%%%% to read in numlay(atmospheres) sets of data  
      %%%%%% but for surface d/dsurface we only read in 4 sets of data
      if jacobItem <= (iGasDQ+2)
        ODBrowmax=numlay(atmospheres);
      else
        ODBrowmax=4;
        end

      % sanity check 
      if nrow ~= ODBrowmax
        fprintf(2, 'readkcjac: WARNING -- bad ODB row count,  ');
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
            error(['fseek failed,odb=',num2str(odb),' chunk=',num2str(chunk)]);
            end
          if fwrite(fout, ODBrowdat, 'real*4') ~= 10000
            error(['fwrite failed,odb=',num2str(odb),'chunk=',num2str(chunk)]);
            end     
          end
        cumODBrow = cumODBrow + 1;
        end % loop on current ODB rows
      end % loop on jacobItem

    end % loop on atmospheres

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
