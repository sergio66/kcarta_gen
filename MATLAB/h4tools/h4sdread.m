
function [alist, fattr] = h4sdread(hfile)

% NAME
%
%   h4sdread -- read HDF4 SDs as a list of matlab arrays
%
% SYNOPSIS
%
%   [alist, fattr] = h4sdread(hfile)
%
% INPUTS
%
%   hfile  - hdf file name
%
% OUTPUT
%
%   alist  - cell list of numeric arrays and attributes
%   fattr  - optional cell array of general attributes
%
% OUTPUT FORMAT
%
%   the variable "alist" is a cell list of one or more elements
%   where each element has the form {name, data}, or optionally
%   {name, data, attrs} or {name, data, attrs, dims}, with
%
%       name    - ascii name for data array
%       data    - numeric data array (can be multidimensional)
%       attrs   - optional cell list of {name, value} attributes
%       dims    - optional cell list of dimension names
%
%
% H. Motteler, 22 Nov 00
%

access_mode = 'read';
sd_id = hdfsd('start', hfile, access_mode);
if sd_id == -1
  error('HDF sdstart failed')
end

% general info about file contents
[ndatasets, nglobal_attr, status] = hdfsd('fileinfo',sd_id);
if status == -1
  error('HDF sdfileinfo failed')
end

% read any file attributes
fattr = {};
for attr_index = 0 : nglobal_attr - 1

  [name,data_type,count,status] = hdfsd('attrinfo',sd_id,attr_index);
  if status == -1
    error('HDF sdattrinfo failed')
  end

  [data,status] = hdfsd('readattr', sd_id, attr_index);
  if status == -1
    error('HDF sdreadattr failed')
  end

  fattr{attr_index + 1} = {name, data};
end

% loop on SDS indices and read SDS's
alist = {};
for sds_index = 0 : ndatasets-1

  % get sds ID from index
  sds_id = hdfsd('select',sd_id, sds_index);

  [name, A, attrs, dims] = sdsid2mat(sd_id, sds_id);

  % close the current SDS
  status = hdfsd('endaccess',sds_id);
  if status == -1
    error('HDF sdendaccess failed')
  end

  alist{sds_index + 1} = {name, A, attrs, dims};

end

% close the file
status = hdfsd('end',sd_id);
if sd_id == -1
  error('HDF sdend failed')
end

