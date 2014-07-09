
function h4sgwrite(hfile, dstr, astr)

% NAME
%
%   h4sgwrite -- write a matlab structure array as a vgroups of SD's
%
% SYNOPSIS
%
%   h4sgwrite(hfile, dstr, astr)
%
% INPUTS
%
%   hfile  - hdf file name
%   dstr   - 1-d array of structures containing data to write
%   astr   - structure of global attributes 
%
% OUTPUT
%
%   an HDF file with a vgroup of SDS's for each structure index
%
% BUGS
%
%   - still "in progress", passes some basic tests
%   - HDF seems to have a hard limit of 5000 SDS's total
%   - reopen and write to an existing big file is slow
%   - needs the capability to update single fields
%   - probably some performance tuning would help, maybe
%     calculate num_dds_block from size of input struct
%
%
% H. Motteler, 14 Nov 00
%

% check if hfile already exists
if hdfh('ishdf', hfile)
  access = 'write';
else
  access = 'create';
end

% open new or existing hfile
file_id = hdfh('open', hfile, access, 4000);
if file_id == -1
  error('HDF hopen failed');
else
  fprintf(2, 'opened %s mode %s\n', hfile, access);
end

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% initialize the SD interface
sd_id = hdfsd('start', hfile, 'write');
if sd_id == -1
  error('HDF sdstart failed');
end

% assign optional top level SD attributes

if nargin == 3
  afields = fieldnames(astr);
  for j = 1 : length(afields)
    aname = afields{j};
    eval(sprintf('aval = astr.%s;', aname));
    status = hdfsd('setattr', sd_id, aname, aval);
    if status == -1
      error('HDF sdsetattr failed')
    end
  end
end


% loop on structure indices / vgroups

for i = 1 : length(dstr)

  % create a new vgroup
  access = 'w';
  vgroup_ref = -1; % flag to create
  vgroup_id = hdfv('attach', file_id, vgroup_ref, access);
  if vgroup_id == -1
    error('HDF vattach failed')
  end

  % assign the group the name "record i"
  sname = sprintf('record %d', i);
  status = hdfv('setname', vgroup_id, sname);
  if status == -1
    error('HDF vsetname failed')
  end

  % assign the group the class "srecord"
  sclass = 'srecord';
  status = hdfv('setclass', vgroup_id, sclass);
  if status == -1
    error('HDF vsetclass failed')
  end

  % loop on structure fields
  dfields = fieldnames(dstr);
  for j = 1 : length(dfields)

    % get field name and value
    dname = dfields{j};
    eval(sprintf('dval = dstr(%d).%s;', i, dname));

    % check for non-empty field
    if ~isempty(dval)
  
      % create an SDS for this field
      sds_id = mat2sdsid(sd_id, dname, dval);

      % add the SDS's to the vgroups
      sds_ref = hdfsd('idtoref', sds_id);
      tag = 720;
      status = hdfv('addtagref', vgroup_id, tag, sds_ref);
      if status == -1
        error('HDF vaddtagref failed')
      end

      % end access to the current SDS
      status = hdfsd('endaccess',sds_id);
      if status == -1
        error('HDF sdendaccess failed')
      end

    end % dval for non-empty field
  end % loop on structure fields

  % end access to the current vgroup
  status = hdfv('detach', vgroup_id);
  if status == -1
    error('HDF vdetach failed')
  end

end % loop on structure indices / vgroups


% end SD interface access
status = hdfsd('end',sd_id);
if status == -1
  error('HDF sdend failed')
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


