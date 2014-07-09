
function [vstru, vattr] = vsfid2mat(file_id, vname);

% NAME
%
%   vsfid2mat -- read an HDF VS into a matlab structure
%
% SYNOPSIS
%
%   [vstru, vattr] = vsfid2mat(file_id, vname)
%
% INPUTS
%
%   file_id   - HDF file ID from hopen
%   vname     - name for the new vdata
%
% OUTPUT
%
%   vstru     - structure of arrays, the main data
%   vattr     - cell array of attributes, optional
%
% OUTPUT DATA FORMATS
%
%   The output vstru is a structure whose fields are fsize-by-nrec
%   arrays.  fsize is the field size and nrec the number of records;
%   nrec should be the same for all fields.
%
%   The output vattr is a cell list of attribute specifications.
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

% --------------------------
% set-up and read the vdata
% --------------------------

% attach to the vdata "vname"
vdata_ref = hdfvs('find', file_id, vname);
if vdata_ref == -1
  error('HDF vsfind failed');
end
access = 'r';
vdata_id = hdfvs('attach', file_id, vdata_ref, access);
if vdata_id == -1
  error('HDF vsattach failed');
end

% get info about the vdata
[n, interlace, fields, nbytes, vdata_name, status] = ...
				hdfvs('inquire', vdata_id);
if status == -1
  error('HDF vsinquire failed');
end

% set the fields to be read
status = hdfvs('setfields', vdata_id, fields);
if status == -1
  error('HDF vssetfields failed');
end

% read all n vdata records
[data, count] = hdfvs('read', vdata_id, n);     % about 60% of runtime
if count ~= n
  error('HDF vsread failed');
end

% -----------------------------------------
% copy the vdata cell array to a structure
% -----------------------------------------

% separate the field strings
cind = findstr(',', fields);
cind = [0, cind, length(fields)+1];
field_list = {};
vfield_list = {};
for i = 1 : length(cind) - 1
  a = cind(i)+1;
  b = cind(i+1)-1;
  field_list{i} = fields(a:b);
  vfield_list{i} = mkvar(fields(a:b));
end

vstru = cell2struct(data, vfield_list, 1);

% ---------------------
% read any attributes
% ---------------------

vattr = {};
nattrs = hdfvs('nattrs',vdata_id);

k = 1; % total found attribute count

for j = 0 : length(field_list)  % loop on fields

  if j == 0
    % special case for "general" vdata attr's
    fname = vname;
    field_index = -1; 
  else
    % vdata field attributes
    fname = field_list{j};
     [field_index, status] = hdfvs('findex', vdata_id, fname);
    if status == -1
      error('HDF vsfindex failed')
    end
  end

  nfattrs = hdfvs('fnattrs', vdata_id, field_index);

  for attr_index = 0 : nfattrs - 1

    [aname, data_type, count, nbytes, status] = hdfvs('attrinfo',...
                    vdata_id, field_index, attr_index);
    if status == -1
      error('HDF vsattrinfo failed')
    end

    [aval, status] = hdfvs('getattr', vdata_id, field_index, attr_index);
    if status == -1
      error('HDF vsgetattr failed')
    end

    vattr{k} = {fname, aname, aval};
    k = k + 1;

  end
end

if nattrs ~= k - 1
  error('attribute count mismatch')
end

% ------------------------
% detach the vs interface
% ------------------------

status = hdfvs('detach', vdata_id);
if status == -1
  error('HDF vsdetach failed')
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

