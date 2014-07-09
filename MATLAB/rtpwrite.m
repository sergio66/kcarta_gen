
function rtpwrite(hfile, head, hattr, prof, pattr)

% NAME
%
%   rtpwrite -- write a set of RTP Version 2.01 profiles
%
% SYNOPSIS
%
%   rtpwrite(hfile, head, hattr, prof, pattr)
%
% INPUTS
% 
%   hfile  - filename for the HDF RTP file
%   head   - structure for the header data 
%   hattr  - header general and field attributes
%   prof   - structure for profile field arrays
%   pattr  - profile general and field attributes
%
% OUTPUT
%
%   an HDF file containing the relevant profile and header data
%
% STRUCTURE FORMATS
%
%   The fields of the profile structure are fsize-by-nrec arrays,
%   where fsize is the field size and nrec the number of records.
%   The fields of the header are fsize-by-one arrays, where fsize 
%   is the header field size.
%
% HDF ATTRIBUTE FORMAT
%
%   hattr and pattr are cell lists of attribute specifications.
%   Each attribute specification is a cell list of three strings of
%   the form
%                      {fname, aname, aval}
%   
%   where fname is either a field name or for general attributes,
%   or "header" or "profiles" for general attributes.  Here is an
%   example of setting profile attributes:
% 
%   pattr = ...
%     { {'profiles', 'comment',  '1761 TIGR Profiles'}, ...
%       {'plevs',    'comment',  'standard AIRS pressure levels'}, ...
%       {'plevs',    'units',    'millibars'}, ...
%       {'gas_1',    'comment',  'water'} }, ...
%       {'gas_1',    'units',    'PPMV'} }
%
%   Fields can have more than one attribute, and do not have 
%   to have any.  Header attributes can be set in a similar way.
%   If no  attributes are to be set, the attribute array should 
%   be set to the empty array.
%
% NOTES
%
%   As a convention, general attributes for the file as a whole
%   should be set in hattr, as "header" attributes.
%
%   Profiles are a structure of arrays rather than an array of
%   structures, because manipulating the former is much faster
%   in Matlab, for large arrays; see stransp1.m and stransp2.m 
%   to go between the two formats.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H. Motteler, 9 July 01;
% Updated for RTP version 2.01, 24 Oct 2008 S.Hannon
% Update: 20 Oct 2011, S.Hannon - write calflag & pnote as uint8

%disp('  Using updated RTPWRITE')

bad = -9999;  % JPL bad value

% --------------------------------------------
% do sanity check of the supplied field sizes
% --------------------------------------------

hfields = fieldnames(head);
for j = 1 : length(hfields);
  fname = hfields{j};
  eval(sprintf('[m,n] = size(head.%s);', fname));
  if n ~= 1
    error('header fields must be column vectors');
  end
end

pfields = fieldnames(prof);
fname = pfields{1};
eval(sprintf('[m,nprof] = size(prof.%s);', fname));
for j = 2 : length(pfields);
  fname = pfields{j};
  eval(sprintf('[m,n] = size(prof.%s);', fname));
  if n ~= nprof
    error('profile structure fields must all have the same number of columns');
  end
end
  
% -------------------------------------
% squeeze data to (mostly) 32-bit types
% -------------------------------------

% HEAD
for i = 1 : length(hfields)
  fname = hfields{i};
  switch fname

    % set the int32 fields
    case {'ptype', 'pfields', 'ngas', 'glist', 'gunit', ...
          'instid', 'pltfid', 'nchan', 'ichan', 'iudef', 'itype'}
      eval(sprintf('head.%s = int32(head.%s);', fname, fname));

    % assume everything else should be a float32
    otherwise
      eval(sprintf('head.%s = single(head.%s);', fname, fname));
  end
end

% PROF
for i = 1 : length(pfields);
  fname = pfields{i};
  switch fname

    % set the int32 fields
    case {'nemis', 'landtype', 'nlevs', 'clrflag', ...
          'ctype', 'ctype2', 'upwell', 'robsqual', ...
          'findex', 'atrack', 'xtrack', 'ifov', 'iudef', 'itype'}
      eval(sprintf('prof.%s = int32(prof.%s);', fname, fname));

    % set the double fields
    % (TAI should be a double)
    case {'ptime', 'rtime'}
      eval(sprintf('prof.%s = double(prof.%s);', fname, fname));

    % set the byte-array fields
    case {'pnote', 'calflag'}
      eval(sprintf('prof.%s = uint8(prof.%s);', fname, fname));

    % assume everything else should be a float32
    otherwise
      eval(sprintf('prof.%s = single(prof.%s);', fname, fname));
  end
end

% -------------------------------
% write profiles as HDF 4 vdatas
% -------------------------------

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

% write the header and profiles
mat2vsfid(file_id, 'header', head, hattr);
mat2vsfid(file_id, 'profiles', prof, pattr);

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

