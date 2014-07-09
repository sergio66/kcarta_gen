
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
%   prof   - structure array of profile data
%   pattr  - profile general and field attributes
%
% STRUCTURE FORMATS
%
%   Profiles are returned as an array of structures.  Fields
%   of the profiles are returned as fsize-by-one column vectors, 
%   where fsize is the field size.  An exception is the pseudo-
%   field gamnt, which is an nlevs-by-ngas array.  The fields of 
%   the header are also fsize-by-one arrays, where fsize is the 
%   header field size.
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
%   rtpread2 is relatively slow, due to performance limitations 
%   of arrays of structures in Matlab, and also due to building 
%   the "gamnt" pseudo-field.  rtpread is faster; but it returns
%   a single structure whose fields are fsize-by-nrec arrays, and
%   does not build the gamnt array.
%
%   As a convention, general attributes for the file as a whole
%   should be set in hattr, as "header" attributes.
%
%
% H. Motteler, 12 Jun 01
%

% call the lower-level RTP reader
[head, hattr, prof, pattr] = rtpread(hfile);

[mlevs,n] = size(prof.plevs);

% transpose the profile data
prof = stransp1(prof);

% if there is no glist, just return prof as-is
if ~isfield(head, 'glist')
  return
end

% build an inverse gas map
k = max(double(head.glist));
ginv = zeros(k,1);
for i = 1 : double(head.ngas)
  ginv(double(head.glist(i))) = i;
end

% bundle gas_<i> fields together to make the gamnt array
pfields = fieldnames(prof);
for i = 1 : length(prof)
  prof(i).gamnt = zeros(mlevs, head.ngas);
  for j = 1 : length(pfields)
    fname = pfields{j};
    if strncmp(fname, 'gas_', 4) & ~strcmp(fname, 'gas_xx')
      gid = str2num(fname(5:length(fname)));
      k = ginv(gid);  % k is index of gid, in glist
      eval(sprintf('prof(i).gamnt(:,k) = prof(i).gas_%d;', gid))
    end
  end
end

