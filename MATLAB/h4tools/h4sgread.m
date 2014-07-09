
function [dstr, astr] = h4sgread(hfile)

% NAME
%
%   h4sgread -- read vgroups of SDS's to a matlab structure array
%
% SYNOPSIS
%
%   [dstr, astr] = h4sgread(hfile)
%
% INPUTS
%
%   hfile  - hdf file name
%
% OUTPUTS
%
%   dstr  - structure of data arrays
%   astr  - structure of global attributes
%
% BUGS
%
%   - still "in progress", passes basic tests
%   - astr is not set if no attr's in file
%   - what is the extra vgroup at the top level?  
%     testing for the "srecords" class is a hack
%
%
% H. Motteler, 14 Nov 00
%


% open the HDF file
file_id = hdfh('open', hfile, 'read', 0);
if file_id == -1
  error('HDF hopen failed');
end

% get number of SDS's in this file
num = hdfh('number',file_id,720);

% initialize the V interface
status = hdfv('start',file_id);
if status == -1
  error('HDF vstart failed');
end

% initialize the SD interface
sd_id = hdfsd('start', hfile, 'read');
if sd_id == -1
  error('HDF sdstart failed');
end

% general info about SD contents
[ndatasets,nglobal_attr,status] = hdfsd('fileinfo',sd_id);
if status == -1
  error('HDF sdfileinfo failed');
end


% read global attributes

% loop on global attributes
for j = 0 : nglobal_attr - 1
  [aname,data_type,count,status] = hdfsd('attrinfo', sd_id, j);
  if status == -1
    error('HDF sdattrinfo failed');
  end

  [adata,status] = hdfsd('readattr', sd_id, j);
  if status == -1
    error('HDF sdreadattr failed');
  end

  % copy attribute name and value
  % astr.<aname> = adata;
  eval(sprintf('astr.%s = adata;', aname));

end


% step thru vgroups, read SDS's as structure fields

% get number of top-level vgroups
[vrefs, count] = hdfv('lone', file_id, 1);
[vrefs, count] = hdfv('lone', file_id, count);

% loop on vgroups / structure records
k = 0; % structure index
for i = 1 : length(vrefs)

  vgroup_id = hdfv('attach', file_id, vrefs(i), 'r');
  if vgroup_id == -1
    error('HDF vattach failed')
  end

  [vgroup_name, status] = hdfv('getname', vgroup_id);
  if status == -1
    error('HDF vgetname failed')
  end

  [class_name, status] = hdfv('getclass', vgroup_id);
  if status == -1
    error('HDF vgetclass failed')
  end

  % see if we have an "srecord" of SDS's for our structure
  if strcmp(class_name, 'srecord')

    k = k + 1; % increment structure index

    % count = hdfv('ntagrefs',vgroup_id) 

    [tag2, refs2, count] = hdfv('gettagrefs',vgroup_id, 1000);

    % loop on contents of the vgroup
    for j = 1 : length(refs2)
      % test for SD tag
      if tag2(j) == 720
        sds_index = hdfsd('reftoindex', sd_id, refs2(j));
        sds_id = hdfsd('select', sd_id, sds_index);
        [sname,rank,dimsizes,data_type,nattrs,status] = ...
				          hdfsd('getinfo',sds_id);
        if status == -1
          error('HDF getinfo failed')
        end

	% read the data, using values from the getinfo call
	start = zeros(1,rank);
	stride = ones(1,rank);
	edge = dimsizes;
	[sdata,status] = hdfsd('readdata',sds_id,start,stride,edge);
        if status == -1
          error('HDF readdata failed')
        end

        status = hdfsd('endaccess',sds_id);
        if status == -1
          error('HDF sdendaccess failed')
        end

	% add the array to the structure to be returned
	% dstr(k).<sname> = sdata;
	eval(sprintf('dstr(%d).%s = sdata;', k, sname));

      end % test for SD tag
    end % loop on vgroup elements
  end % test for vgroup of class 'srecord'

  % end access to the current vgroup
  status = hdfv('detach', vgroup_id);
  if status == -1
    error('HDF vdetach failed')
  end

end % loop on vgroups


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

