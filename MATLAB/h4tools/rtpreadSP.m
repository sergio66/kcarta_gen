
function [head, hattr, prof, pattr] = rtpread(hfile)

% NAME
%
%   rtpread -- read a set of RTP profiles
%
% SYNOPSIS
%
%   [head, hattr, prof, pattr] = rtpread(hfile)
%
% INPUTS
% 
%   hfile  - filename for the HDF RTP file
%
% OUTPUT
%
%   head   - structure for the header data 
%   hattr  - header general and field attributes
%   prof   - structure for profile field arrays
%   pattr  - profile general and field attributes
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
%   rtpread returns a structure whose fields are arrays where
%   the columns are the individual profiles.  It does not build
%   a 2D constituent array gamnt.  rtpread2 returns RTP data as
%   a structure array, with one structure record per profile, 
%   and builds the gamnt array from the individual constituent 
%   fields.
%
% H. Motteler, 12 Jun 01
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

% read the header and profiles
[head, hattr] = vsfid2mat(file_id, 'header');
[prof, pattr] = vsfid2mat(file_id, 'profiles');

% move data to 64 bit types, as Matlab does not support 
% arithmetic directly on int32 or float32 types

% expand header fields to double
% hfields = fieldnames(head);
% for i = 1 : length(hfields)
%   fname = hfields{i};
%   eval(sprintf('head.%s = double(head.%s);', fname, fname));
% end

% expand profile fields to double
% pfields = fieldnames(prof);
% for j = 1 : length(pfields)
%   fname = pfields{j};
%   switch fname
% 
%     % calflag and pnote should be strings
%     % if these are already strings, "char" doesn't have any effect
%     case {'calflag', 'pnote'}
%       eval(sprintf('prof.%s = char(prof.%s);', fname, fname));
% 
%     % return everything else as a double
%     otherwise
%       eval(sprintf('prof.%s = double(prof.%s);', fname, fname));
%   end
% end

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

