
function rtpwrite2(hfile, head, hattr, prof, pattr)

% NAME
%
%   rtpwrite2 -- write a set of RTP profiles
%
% SYNOPSIS
%
%   rtpwrite2(hfile, head, hattr, prof, pattr)
%
% INPUTS
% 
%   hfile  - filename for the HDF RTP file
%   head   - structure for the header data 
%   hattr  - header general and field attributes
%   prof   - structure array of profile data
%   pattr  - profile general and field attributes
%
% OUTPUT
%
%   an HDF file containing the relevant profile and header data
%
% STRUCTURE FORMATS
%
%   Profiles are represented as an array of structures.  Fields should
%   generally be fsize-by-one column vectors, where fsize is the field
%   size.  An exception is the pseudo-field gamnt, which is an nlevs-by-
%   ngas array.  The pnote string can be a row.
%   
%   The profile surface fields rho, emis, mwemis, mwstb, the atmospheric
%   fields plevs, plays, palts, ptemp, and the comment field pnote can
%   vary in size from profile to profile.  To write these out as HDF
%   vdata records, these fields are filled out to their longest lengths
%   as found in prof, and the corresponding RTP size field is added, if 
%   it is not already present.
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
%   rtpwrite2 is relatively slow, due to performance limitations
%   of arrays of structures in Matlab, and also due to unpacking
%   the "gamnt" pseudo-field and filling variable-sized fields.
%
% H. Motteler, 12 Jun 01
%

bad = -9999;  % JPL bad value

% find the max sizes for profile arrays
nprof = length(prof);

% get the max number of IR reflectance points
mrho = 0;
if isfield(prof, 'rho')
  if isfield(prof, 'nrho')
    % use the supplied nrho
    for i = 1 : nprof
      mrho = max(mrho, prof(i).nrho);
    end
  else
    for i = 1 : nprof
      prof(i).nrho = length(prof(i).rho);
      mrho = max(mrho, prof(i).nrho);
    end
  end
end
rhofill = ones(mrho,1) * bad;

% get the max number of IR emiss points
memis = 0;
if isfield(prof, 'emis')
  if isfield(prof, 'nemis')
    % use the supplied nemis
    for i = 1 : nprof
      memis = max(memis, prof(i).nemis);
    end
  else
    % add an nemis field
    for i = 1 : nprof
      prof(i).nemis = length(prof(i).emis);
      memis = max(memis, prof(i).nemis);
    end
  end
end
emisfill = ones(memis,1) * bad;

% get the max number of MW emis points
mwmemis = 0;
if isfield(prof, 'mwemis')
  if isfield(prof, 'mwnemis')
    % use the supplied mwnemis
    for i = 1 : nprof
      mwmemis = max(mwmemis, prof(i).mwnemis);
    end
  else
    % add a mwnemis field
    for i = 1 : nprof
      prof(i).mwnemis = length(prof(i).mwemis);
      mwmemis = max(mwmemis, prof(i).mwnemis);
    end
  end
end
mwemisfill = ones(mwmemis,1) * bad;

% get the maximum number of MW surf Tb points
mwmstb = 0;
if isfield(prof, 'mwstb')
  if isfield(prof, 'mwnstb')
    % use the supplied mwnstb
    for i = 1 : nprof
      mwmstb = max(mwmstb, prof(i).mwnstb);
    end
  else
    % add a mwnstb field
    for i = 1 : nprof
      prof(i).mwnstb = length(prof(i).mwstb);
      mwmstb = max(mwmstb, prof(i).mwnstb);
    end
  end
end
mwstbfill = ones(mwmstb,1) * bad;

% get the max number of levels
mlevs = 0;
if isfield(prof, 'plevs')
  if isfield(prof, 'nlevs')
    % use the supplied nlevs
    for i = 1 : nprof
      mlevs = max(mlevs, prof(i).nlevs);
    end
  else
    % add an nlevs field
    for i = 1 : nprof
      prof(i).nlevs = length(prof(i).plevs);
      mlevs = max(mlevs, prof(i).nlevs);
    end
  end
end
levfill = ones(mlevs,1) * bad;

% get length of the longest pnote string
mpnote = 0;
if isfield(prof, 'pnote')
  for i = 1 : nprof
    mpnote = max(mpnote, length(prof(i).pnote));
  end
end
pnotefill = char(ones(mpnote,1)*double(' '));

% loop on profiles, fill variable-sized fields to max val's
for i = 1 : nprof
    
  % loop on profile fields
  pfields = fieldnames(prof);
  for j = 1 : length(pfields);
    fname = pfields{j};
    switch fname

      % fill rho to mrho
      case {'rho '}
        t = rhofill;
        k = prof(i).nrho;
        t(1:k) = prof(i).rho(1:k);
        prof(i).rho = t;

      % fill emis to memis
      case {'emis'}
        t = emisfill;
        k = prof(i).nemis;
        t(1:k) = prof(i).emis(1:k);
        prof(i).emis = t;

      % fill mwemis to mwmemis
      case {'mwemis'}
        t = mwemisfill;
        k = prof(i).mwnemis;
        t(1:k) = prof(i).mwemis(1:k);
        prof(i).mwemis = t;

      % fill mwstb to mwmstb
      case {'mwstb'}
        t = mwstbfill;
        k = prof(i).mwstb;
        t(1:k) = prof(i).mwstb(1:k);
        prof(i).mwstb = t;

      % fill profile fields to mlevs
      case {'plevs','plays','palts','ptemp'}
        t = levfill;
        k = prof(i).nlevs;
        eval(sprintf('t(1:k) = prof(i).%s(1:k);', fname));
        eval(sprintf('prof(i).%s = t;', fname));
    
      % save gamnt as gas_<n> fields and fill to mlevs
      case {'gamnt'}
        k = prof(i).nlevs;
        for gind = 1:length(head.glist)
          t = levfill; 
	  t(1:k) = prof(i).gamnt(1:k,gind);
          eval(sprintf('prof(i).gas_%d = t;', head.glist(gind)));
        end
        prof(i).gamnt = bad;  % don't pass the gamnt array on

      % fill pnote to mpnote
      case {'pnote'}
        t = pnotefill;
        k = length(prof(i).pnote);
        t(1:k) = prof(i).pnote;
        prof(i).pnote = t;

    end
  end
end

% set head.ngas
if isfield(head, 'glist')
  if ~isfield(head, 'ngas')
    head.ngas = length(head.glist);
  end
else
  head.ngas = 0;
end

% set head.nchan
if isfield(head, 'vchan')
  if ~isfield(head, 'nchan')
    head.nchan = length(head.vchan);
  end
elseif isfield(head, 'ichan')
  if ~isfield(head, 'nchan')
    head.nchan = length(head.ichan);
  end
else
  head.nchan = 0;
end

% set head.mwnchan
if isfield(head, 'mwfchan')
  if ~isfield(head, 'mwnchan')
    head.mwnchan = length(head.mwfchan);
  end
else
  head.mwnchan = 0;
end

% transpose the profile data
prof = stransp2(prof);

% call rtpwrite to write the data
rtpwrite(hfile, head, hattr, prof, pattr);

