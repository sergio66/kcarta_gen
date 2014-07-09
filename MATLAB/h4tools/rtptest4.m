
% save Scott's 48 fitting profiles in RTP format
% 
% 

% The fitting profiles, in GENLN input (level) format,
% with filenames of the form "myp1", "myp2", etc.
profdir = '/home/motteler/asl/profiles/fitpro_july00/';

% nprof = 48;
nprof = 48;

% hdf output file
hfile = 'test4.hdf';

% -----------------------------------------
% read the profiles into a structure array
% -----------------------------------------

clear prof
pmax = 0;
pmin = 10000;
for i = 1 : nprof

  % read a GENLN2 format profile
  pname = sprintf('myp%d', i);
  pfile = [profdir, pname];
  [plevs, ptemp, glist, gamnt, plat] = gproread(pfile);

  pmax = max([plevs(:); pmax]);
  pmin = min([plevs(:); pmin]);

  if i == 1
    glist1 = glist;
  elseif glist1 ~= glist
    i, glist, glist1
    error('HDF vdata RTP profiles must all have the same gas set')
  end

  % save the profile values as matlab structure fields
  prof(i).plat = plat;
  prof(i).nlevs = length(plevs);
  prof(i).plevs = plevs(:);
  prof(i).ptemp = ptemp(:);
  prof(i).gamnt = gamnt;

end

head.glist = glist;
head.pmin = pmin;
head.pmax = pmax;

%---------------------
% set demo attributes
%---------------------

pattr = { {'plat', 'units', 'millibars'}, ...
          {'ptemp', 'units', 'kelvins'} };

k = length(pattr);
for j = 1 : length(glist)
  gid = glist(j);
  gname = sprintf('gas_%d', gid);
  pattr{j + k -1} = {gname, 'units', 'PPMV'};
end

% set sample header attributes
hattr = { {'header', 'author', 'S. Hannon; RTP version, H. Motteler'}, ...
          {'header', 'date', '4 Dec 00'}, ...
          {'header', 'comment', '48 fitting profiles in RTP format'} };

% write the test data
rtpwrite2(hfile, head, hattr, prof, pattr);

% read it back in
[head2, hattr2, prof2, pattr2] = rtpread2(hfile);

% compare
isequal(hattr, hattr2)
isequal(pattr, pattr2)

structcmp(head, head2)

structcmp(prof, prof2)

return

% the following assumes ftest2 copies copies its input to ftest2.hdf
! cp test4.hdf rtp
! cd rtp && ftest2 > test.log 2>&1
! tail rtp/test.log

[head3, hattr3, prof3, pattr3] = rtpread2('rtp/ftest2.hdf');


