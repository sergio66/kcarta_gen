
function h4vswrite(hfile, vlist)

% NAME
%
%   h4vswrite -- write a structure of arrays as an HDF VS 
%
% SYNOPSIS
%
%   h4vswrite(hfile, vlist)
%
% INPUTS
%
%   hfile  - hdf file name
%   vlist  - cell list of structures arrays and attributes
%
% OUTPUT
%
%   an HDF file with a vdata for each matlab structure array
%
% INPUT DATA FORMAT
%
%   The input vlist is a cell array with one element for each 
%   vdata to be written.  Each element of vlist has the form 
%
%                {vname, vstru, vattr, vclass},   
%   where
%
%    - vname is the vdata name
%
%    - vstru is a structure whose fields are fsize-by-nrec arrays.
%      fsize is the field size and nrec is the the number of records; 
%      nrec should be the same for all fields.
%
%    - vattr is a cell list of attribute specifications.  Each attribute 
%      specification is a list of three strings of the form
%
%                      {fname, aname, aval}
%   
%      where fname is either a field name or for general attributes
%      the vdata name, aname is the attribute name, and aval the 
%      attribute value.
%
%    - vclass is an optional vdata class name; if this is not
%      specified it defaults to "structure array"
%
% NOTES
%
%   vstru is a structure of arrays rather than an array of
%   structures, because manipulating the former is much faster
%   in Matlab, for large arrays; see stransp1.m and stransp2.m 
%   to go between the two formats.
%   
%
% H. Motteler, 9 July 01
%

% default class for structure arrays
sclass = 'structure array';

% create a new hdf file
access = 'create';
file_id = hdfh('open', hfile, access, 0);
if file_id == -1
  error('HDF hopen failed');
end

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% loop on vdatas in the vlist
for i = 1 : length(vlist)

  % write the i-th vdata
  switch length(vlist{i})

    case 2   
      % data with no attributes or vdata class
      mat2vsfid(file_id, vlist{i}{1}, vlist{i}{2}, {}, sclass)

    case 3
      % data with attributes but no vdata class 
      mat2vsfid(file_id, vlist{i}{1}, vlist{i}{2}, vlist{i}{3}, sclass);

    case 4
      % data with attributes and vdata class
      mat2vsfid(file_id, vlist{i}{1}, vlist{i}{2}, vlist{i}{3}, vlist{i}{4});

    otherwise
      error('wrong number of elements in vlist');
  end
end

% end vgroup interface access
status = hdfv('end',file_id);
if status == -1
  error('HDF vend failed')
end

% close the HDF file
status = hdfh('close',file_id);
if status == -1
  error('HDF hclose failed')
end

