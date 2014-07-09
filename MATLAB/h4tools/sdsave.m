%
% NAME
%
%   sdsave -- write matlab numeric arrays as HDF4 SDs
%
% SYNOPSIS
%
%   sdsave(hfile, dat, atr, dim)
%
% INPUTS
%
%   hfile  - hdf file name
%   dat    - structure of arrays of numeric data
%   atr    - optional structure of file and variable attributes
%   dim    - optional structure of variable dimension attributes
%
% OUTPUT
%
%   an HDF file with an SDS for each field in dat
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
%   Variable can have more than one attribute, so variable attributes
%   are cell arrays of the form {{name, value}, {name, value}, ...}
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
%   Dimension names should have a single size associated with them.
%   If you name a dimension "latitude" for the first dimension of
%   an m x n variable V1, and then try to name the first dimension
%   of a j x k variable V2 "latitude", for j ~= m, you will get an
%   error.
%
% BUGS
% 
%   Compatiblilty of the HDF SD's and NetCDF formats is limited by 
%   the libraries used.  sdload.m should read both HDF and NetCDF 
%   without any problems, but the output of sdsave.m is HDF, and may 
%   not be readable by NetCDF tools unless they were compiled with 
%   the HDF libraries.
% 
%   Matlab variables are inherently at least 2-dimensional, so a
%   load from a non-matlab source followed by a save will result
%   in the HDF libraries adding a "fakeDim<n>" dimension name to 
%   scalar and 1-D vector values.
%
% H. Motteler, 31 Oct 02
%

function sdsave(hfile, dat, atr, dim)

% create a new HDF SD file
sd_id = hdfsd('start', hfile, 'create');
if sd_id == -1
  error('HDF sdstart failed')
end

% set any file attributes
if nargin >= 3
  afields = fieldnames(atr);
  for j = 1 : length(afields)
    aname = afields{j};
    % only do general attributes here
    if ~isfield(dat, aname)
      status = hdfsd('setattr', sd_id, aname, getfield(atr,  aname));
      if status == -1
        error('HDF sdsetattr failed')
      end
    end
  end
end

% default is no dimension data
if nargin < 4
  dim = {};
end

% default is no attribute data
if nargin < 3
  atr = {};
end

% make sure we have a minimal set of parameters
if nargin < 2
  error('sdsave needs an HDF filename and a structure of data to save\n');
end

% loop on "dat" structure fields
dfields = fieldnames(dat);
for j = 1 : length(dfields)

  % current field name
  dname = dfields{j};

  % copy current field's data
  dval1 = getfield(dat, dname);

  % make sure dval1 is numeric data
  % if ~isnumeric(dval1)
  %   fprintf(2, 'sdsave warning: %s is not numeric data, not saved\n', dname);
  %   continue
  % end    

  % check for an attribute field
  if isfield(atr, dname)
    aval1 = getfield(atr, dname);
    if ~isvarattr(aval1)
      aval1 = {};
      fprintf(2, 'sdsave warning: bad attribute format, field %s\n', dname);
    end
  else
    aval1 = {};
  end

  % check for a dimension field
  if isfield(dim, dname)
    mval1 = getfield(dim, dname);
    if ~isvardim(mval1)
      mval1 = {};
      fprintf(2, 'sdsave warning: bad dimension format, field %s\n', dname);
    end
  else
    mval1 = {};
  end

  % create an SD for this field
  sds_id = mat2sdsid(sd_id, dname, dval1, aval1, mval1);

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

% isvarattr tests if its argument is a cell array of the form
% {{name, value}, {name, value}, ...}, where "name" and "value"
% are character arrays
%
function b = isvarattr(d)
b = 1;
if ~iscell(d)
  fprintf(2, 'sdsave: variable attribute list must be a cell list\n');
  fprintf(2, '        of the form {{name, value}, {name, value}, ...}\n');
  b = 0;
  return
end
for i = 1 : length(d)
  if ~(iscell(d{i}) & length(d{i}) == 2)
    fprintf(2, 'sdsave: variable attribute list must be a cell list\n');
    fprintf(2, '        of the form {{name, value}, {name, value}, ...}\n');
    b = 0;
    return
  end
  if ~(ischar(d{i}{1}) & ischar(d{i}{2}))
    fprintf(2, 'sdsave: attribute names and values must be char strings\n')
    b = 0;
    return
  end
end

% isvardim tests if its argument is a cell array of the form
% {name, name, ...}, where "name" is a character array
%
function b = isvardim(d)
b = 1;
if ~iscell(d)
  fprintf(2, 'sdsave: variable dimension list must be a cell array\n');
  b = 0;
  return
end
for i = 1 : length(d)
  if ~ischar(d{i})
    fprintf(2, 'sdsave: variable dimensions must be character strings\n')
    b = 0;
    return
  end
end

