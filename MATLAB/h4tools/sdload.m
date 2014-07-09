%
% NAME
%
%   sdload -- read HDF4 SD's and NetCDF files as Matlab arrays
%
% SYNOPSIS
%
%   [dat, atr, dim] = sdload(hfile)
%
% INPUT
%
%   hfile  - HDF SD or NetCDF file name
%
% OUTPUT
%
%   dat    - structure whose fields are arrays, one for each SDS
%   atr    - optional structure of file and variable attributes
%   dim    - optional structure of variable dimension attributes
%
% ATTRIBUTES
% 
%   The attribute structure atr may contain both file and variable
%   attributes.  File attributes are fields of atr with names that do
%   not match any variable name.  The attribute values are character
%   arrays.
% 
%   If dat.df is a field of the dat structure, i.e. a variable to be
%   saved, then atr.df is an optional corresponding variable attribute.
%   Variables can have more than one attribute, so variable attributes
%   are cell arrays of the form { {name, value}, {name, value}, ... }
%   with "name" the attribute name, and "value" the value.  Both name
%   and value should be character arrays.
% 
% DIMENSIONS
% 
%   The optional HDF4 dimensions are similar to variable attributes, but
%   don't have explicit names.
% 
%   If dat.df is a field of the dat structure, i.e. a variable to be
%   saved, then dim.df is an optional corresponding dimension attribute.
%   These are cell arrays of the form {name, name, ...}, with one name
%   for each dimension of dat.df.  The names should be character arrays.
% 
% BUGS
% 
%   Not all HDF SDS names can be represented as Matlab variables; 
%   to work around this, any character in an SDS name that's not a
%   letter or digit is replaced with an underscore.  If that is a 
%   problem, the lower lever h4sdread can be used instead of sdload.
%
%   Compatiblilty of the HDF SD's and NetCDF formats is limited by 
%   the libraries used.  sdload.m should read both HDF and NetCDF 
%   without any problems, but the output of sdsave.m is HDF, and may 
%   not be readable by NetCDF tools unless they were compiled with 
%   the HDF libraries.
% 
% H. Motteler, 31 Oct 02
%

function [dat, atr, dim] = sdload(hfile)

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

  % copy file attrs to atr structure
  vname = mkvar(name);
  eval(sprintf('atr.%s = data;', vname));
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

  % copy data, attrs, and dims to structs
  vname = mkvar(name);
  eval(sprintf('dat.%s = A;', vname));
  eval(sprintf('atr.%s = attrs;', vname));
  eval(sprintf('dim.%s = dims;', vname));
end

% close the file
status = hdfsd('end',sd_id);
if sd_id == -1
  error('HDF sdend failed')
end

%
% function s2 = mkvar(s1);
%
% replace everything in s1 that's not a letter or digit with '_'
%
function s2 = mkvar(s1);

vok = isletter(s1) | ('0' <= s1 & s1 <= '9');
s2 = s1;
s2(~vok) = '_';
if s2(1) == '_'
  s2(1) = 'x';
end

