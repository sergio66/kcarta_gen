
% A demo using h4sdwrite() to save SRF's
%
% NOTES:
%
% 1. all data and attributes are loaded from the .mat file
%
% 2. SDS attribute assignment fails in Matlab 5.3/R11, but 
%    works OK in R12
%
% 3. don't casually change the SDS array or attribute names used 
%    below; this is what the outside world sees as defining the 
%    dataset, and users may want to use these names in their SRF 
%    readers
%
% 4. don't accidentally transpose any of the arrays or vectors;
%    as with the array names, users are likely to expect whatever
%    they first see, and may code things with these assumptions.
%    It should be OK to change the size of dimensions, e.g., from
%    471 to 500 tabulation points
% 
% 5. as far as possible, the data type that is passed to h4sdwrite
%    is what is saved in the HDF file; see htype.m for the list of 
%    SD types that matlab can write.  Be careful with changing data
%    types; many users will have these hard-coded in their readers
%
% 6. array names, dimension order and sizes, and data types can 
%    all be easily double checked in the resulting HDF file, with
%    utilities like HDFLook
%
% 7. HDF stores data in row order, while Matlab stores it in column
%    order.  When we represent the respective data formats with the
%    dimension lists flipped, as below, no explicit transpose is
%    needed to go from one format to the other--in the binary read of 
%    HDF data, in the test code below, the HDF data is read in and 
%    compared directly with the Matlab data; no transpose is done
%    before the comparison.
%
%       Matlab column-order view of the SRF data
%
%          chanid     1 x 2378
%          freq       1 x 2378
%          fwgrid   471 x 1   
%          srfval   471 x 2378
%          width      1 x 2378
%
%      HDFLook row-order view of the SRF data
%
%         chanid    2378 x 1
%         freq      2378 x 1
%         fwgrid       1 x 471 
%         srfval    2378 x 471
%         width     2378 x 1
%
% H. Motteler, 28 Nov 00
%

% set filenames
hfile = '../srfhdf/srftablesV10.hdf';
mfile = '../srfhdf/srftablesV10.mat';

% get the data
load(mfile)

% convert to more compact data types; only the channel centers
% (i.e., freq) need to be saved in the default double precision
chanid = int16(chanid);
fwgrid = single(fwgrid);
srfval = single(srfval);
width = single(width);

% set file attributes
fattr = { ...
  {'author',  author}, ...
  {'version', version}, ...
  {'comment', comment} };

% set array names, values, attributes, and dimension info
slist = { ...
  { 'chanid', chanid, {{'units', chanid_units}}, {'one', 'nchan'}}, ...
  { 'freq',   freq,   {{'units', freq_units}},   {'one', 'nchan'}}, ...
  { 'fwgrid', fwgrid, {{'units', fwgrid_units}}, {'npts', 'one'}}, ...
  { 'srfval', srfval, {{'units', srfval_units}}, {'npts', 'nchan'}},...
  { 'width',  width,  {{'units', width_units}},  {'one', 'nchan'}} };

% call h4sdwrite to write the HDF file
h4sdwrite(hfile, slist, fattr);

% call h4sdread to read the data back, for a check
[slist1, fattr1] = h4sdread(hfile);

% compare original and re-read data
[isequal(slist, slist1), isequal(fattr, fattr1)]

% dump HDF SRF as binary data and check results again
cmd = sprintf('!hdp dumpsds -b -n srfval %s > srfval.bin', hfile);
eval(cmd);
fid = fopen('srfval.bin', 'r');
% [srf1, count] = fread(fid, inf, 'double');
[srf1, count] = fread(fid, inf, 'single');
srf1 = reshape(srf1, 471, 2378);
isequal(srfval, srf1)

% clean up
! rm srfval.bin 2> /dev/null

