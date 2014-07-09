
function h4sdwrite(hfile, slist, fattr)

% NAME
%
%   h4sdwrite -- write a list of matlab arrays as HDF4 SDs
%
% SYNOPSIS
%
%   h4sdwrite(hfile, slist, fattr)
%
% INPUTS
%
%   hfile  - hdf file name
%   slist  - cell list of numeric arrays and attributes
%   fattr  - optional cell listof general attributes
%
% OUTPUT
%
%   an HDF file with an SDS for each array in slist
%
% INPUT FORMAT
%
%   the variable "slist" is a cell list of one or more elements
%   where each element has the form {name, data}, or optionally
%   {name, data, attrs} or {name, data, attrs, dims}, with
%
%       name    - ascii name for data array
%       data    - numeric data array (can be multidimensional)
%       attrs   - optional cell list of {name, value} attributes
%       dims    - optional cell list of dimension names
%
% EXAMPLE 1
%
%   % set array names and values
%   slist = { {'A1', A1}, {'A2', A2} };
%
%   % write the HDF file "Adata"
%   h4sdwrite('Adata', slist);
%
% EXAMPLE 2
% 
%   % set file attributes
%   fattr = { ...
%     {'author',  'Scott Hannon'}, ...
%     {'version', '1.0'}, ...
%     {'comment', 'freqgrid = fwgrid*width + freq' } };
% 
%   % set array names, values, single attributes and dimension info
%   slist = { ...
%     { 'chanid', chanid, {{'units', 'integer'}}, {'one',  'nchan'}}, ...
%     { 'freq',   freq,   {{'units', '1/cm'}},    {'one',  'nchan'}}, ...
%     { 'fwgrid', fwgrid, {{'units', 'FWHM'}},    {'npts', 'one'}}, ...
%     { 'srfval', srfval, {{'units', 'peak=1'}},  {'npts', 'nchan'}},...
%     { 'width',  width,  {{'units', '1/cm'}},    {'one',  'nchan'}}  };
%
%   % write the HDF file "srfV10"
%   h4sdwrite('srfV10', slist, fattr);
%
%
% H. Motteler, 22 Nov 00
%

% create a new HDF SD file
sd_id = hdfsd('start', hfile, 'create');
if sd_id == -1
  error('HDF sdstart failed')
end

% set any file attributes
if nargin == 3
  for i = 1 : length(fattr)
    status = hdfsd('setattr', sd_id, fattr{i}{1}, fattr{i}{2});
    if status == -1
      error('HDF sdsetattr failed')
    end
  end
end

% write the individual SD's out
for i = 1 : length(slist)

  % cases for optional arg's 
  if length(slist{i}) == 4
    sds_id = ...
      mat2sdsid(sd_id, slist{i}{1}, slist{i}{2}, slist{i}{3}, slist{i}{4});
  elseif length(slist{i}) == 3
    sds_id = ...
      mat2sdsid(sd_id, slist{i}{1}, slist{i}{2}, slist{i}{3});
  elseif length(slist{i}) == 2
    sds_id = ...
      mat2sdsid(sd_id, slist{i}{1}, slist{i}{2});
  else
    error('bad "slist" value')
  end

  % close the current SDS
  status = hdfsd('endaccess',sds_id);
  if status == -1
    error('HDF sdendaccess failed')
  end
end

% close the file
status = hdfsd('end',sd_id);
if status == -1
  error('HDF sdend failed')
end

