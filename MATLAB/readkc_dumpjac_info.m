function [data,wnums,iNumLayer] = readkc_dumpjac_info(kfile)

%% see SUBROUTINE DumpJacobianInfo in s_writefile.f90 and jac_down.f90

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

% start, stop frequency
flen = fread(fin, 1, 'integer*4');
fmin = fread(fin, 1, 'real*4');
fmax = fread(fin, 1, 'real*4');
flen = fread(fin, 1, 'integer*4');

% number of layers
flen   = fread(fin, 1, 'integer*4');
nlayer = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
iNumLayer = nlayer;              %<------------- my modification

% which type of data (gas jac, temp jac, raaRad etc
flen   = fread(fin, 1, 'integer*4');
iWhichData = fread(fin, 1, 'integer*4');
flen   = fread(fin, 1, 'integer*4');
iNumLayers=nlayer;              %<------------- my modification


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read "output data blocks"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iProfLayer = 100;

for odb = 1:iProfLayer
  flen1 = fread(fin, 1, 'integer*4');
  [odbdat, count] = fread(fin, [10000,1], 'real*4');
  if count ~= 10000
    error(['fread failed, odb=',num2str(odb),' chunk=',num2str(chunk)]);
  end
  flen2 = fread(fin, 1, 'integer*4');
  data(odb,:) = odbdat;   
end % loop on layers

fclose(fin);

dvx  = 0.0025;
ind = 1:10000;
wnums = fmin + (ind-1)*dvx;
data = data';
