
function vdata_id = mat2vsfid(file_id, vname, vstru, vattr, vclass)

% NAME
%
%   mat2vsfid  - write a matlab structure to an HDF VS
%
% SYNOPSIS
%
%   vdata_id = mat2vsfid(file_id, vname, vstru, vattr, vclass) 
%
% INPUTS
%
%   file_id   - HDF file ID from hopen
%   vname     - name for the new vdata
%   vstru     - structure of arrays, the main data
%   vattr     - cell array of attributes, optional
%   vclass    - vdata class, char string, optional
%
% OUTPUT
%
%   vdata_id  - the vdata ID, after a detach
%
% INPUT DATA FORMATS
%
%   The input vstru is a structure whose fields are fsize-by-nrec
%   arrays.  fsize is the field size and nrec the number of records;
%   nrec should be the same for all fields.
%
%   The input vattr is a cell list of attribute specifications.
%   Each attribute specification is a list of three strings, of 
%   the form
%                      {fname, aname, aval}
%   
%   where fname is either a field name or for general attributes the
%   vdata name, aname is the attribute name, and aval the attribute
%   value.
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

% set defaults
if nargin < 5
  vclass = 'struct array';
end
if nargin < 4
  vattr = {};
end

% ------------------
% create a new vdata
% ------------------

access = 'w';
vdata_ref = -1; % flag to create
vdata_id = hdfvs('attach', file_id, vdata_ref, access);
if vdata_id == -1
  error('HDF vsattach failed');
end

% give it a name and class
status = hdfvs('setname', vdata_id, vname);
if status == -1
  error('HDF vssetname failed');
end

status = hdfvs('setclass', vdata_id, vclass);
if status == -1
  error('HDF vssetclass failed');
end


% -----------------------------------------------------
% get structure field names and define the vdata fields
% -----------------------------------------------------

dfields = fieldnames(vstru);
for j = 1 : length(dfields)

  % get field name and value
  fname = dfields{j};

  % find type and size of fields the first record
  eval(sprintf('ftype = class(vstru.%s);', fname));
  ftype = htype(ftype);

  eval(sprintf('[fsize, nrec] = size(vstru.%s);', fname));

  status = hdfvs('fdefine', vdata_id, fname, ftype, fsize);
  if status == -1
    error('HDF vsfdefine failed')
  end

  % build the field list
  if j == 1
    fieldlist = fname;
  else
    fieldlist = [fieldlist,',',fname];
  end
end

status = hdfvs('setfields', vdata_id, fieldlist);
if status == -1
  error('HDF vssetfields failed')
end

%----------------------
% write any attributes
%----------------------

if ~isempty(vattr)
  for j = 1 : length(vattr)

    aset = vattr{j};
    fname = aset{1};     % vdata or field name
    aname = aset{2};	 % attribute name
    aval = aset{3};	 % attribute value

    if strcmp(fname, vname)
      % if fname == vname, we have a "general" attribute spec 
      field_index = -1;
    else 
      % we have a field attribute specification
      % see if fname really is a field name
      [field_index, status] = hdfvs('findex', vdata_id, fname);
      if status == -1
         error(sprintf('%s is not a field name', fname))
      end
    end

    % set the field attribute
    status = hdfvs('setattr', vdata_id, field_index, aname, aval);
    if status == -1
       error('vssetattr, set vdata field attribute failed')
    end
  end
end  

% ----------------------------------------------------------
% translate the input structure to a cell array for vswrite
% ----------------------------------------------------------

% dcell = {};
% for j = 1 : length(dfields)
% 
%   % get field name and value
%   fname = dfields{j};
% 
%   dcell{j} = eval(sprintf('[vstru(:).%s]', fname));
% end

dcell = struct2cell(vstru)';

clear vstru

% ----------------
% write the vdata
% ----------------

status = hdfvs('write', vdata_id, dcell);
if status == -1
  error('HDF vswrite failed')
end

clear dcell

% detach vdata_id

status = hdfvs('detach', vdata_id);
if status == -1
  error('HDF vsdetach failed')
end

