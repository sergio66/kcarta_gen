
function vlist = h4vsread(hfile)

% NAME
%
%   h4vsread -- read an HDF VS as a structure of arrays
%
% SYNOPSIS
%
%   vlist = h4vsread(hfile) 
%
% INPUTS
%
%   hfile  - hdf file name
%
% OUTPUT
%
%   vlist  - cell list of structures arrays and attributes
%
% OUTPUT DATA FORMAT
%
%   The input vlist is a cell array with one element for each 
%   vdata to be written.  Each element of vlist has the form 
%
%                 {vname, vstru, vattr, vclass},   
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
%    - vclass is an optional vdata class name
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

% open the HDF file
file_id = hdfh('open', hfile, 'read', 0);
if file_id == -1
  error('HDF hopen failed');
end

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% get number of top-level vdatas
[vrefs, count] = hdfvs('lone', file_id, 1000);

% loop on groups
j = 0;
vlist = {};
for i = 1 : length(vrefs)

  % get the id, class, and name of the next vdata
  vdata_id = hdfvs('attach', file_id, vrefs(i), 'r');
  if vdata_id == -1
    error('HDF vsattach failed');
  end

  [class_name,status] = hdfvs('getclass',vdata_id);
  if status == -1
    error('HDF vsgetclass failed');
  end

  [vname, status] = hdfvs('getname', vdata_id);
  if status == -1
    error('HDF vsgetname failed');
  end

  isattr = hdfvs('isattr',vdata_id);

  status = hdfvs('detach', vdata_id);
  if status == -1
    error('HDF vsdetach failed');
  end

  if ~isattr

    % read the vdata
    [vstru, vattr] = vsfid2mat(file_id, vname);

    % add results to output cell array
    j = j + 1;
    vlist{j}{1} = vname;
    vlist{j}{2} = vstru;
    vlist{j}{3} = vattr;
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

