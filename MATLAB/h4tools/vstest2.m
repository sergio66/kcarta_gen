
%
% test of h4vswrite and h4vsread
%
% note: to dump headers: hdp dumpvd -h vstest2.hdf
%

% generate some test data
% ------------------------

% set hdf file and vs names
vname1 = 'test vs #1';
vname2 = 'test vs #2';
hfile = 'vstest2.hdf';
! rm vstest2.hdf 2> /dev/null

% generate test structure 1
for j = 1:10
  s1(j).field1 = rand(4,1);
  s1(j).field2 = rand(7,1);
  s1(j).field3 = rand(3,1);
end
s1 = stransp2(s1);

% generate test structure 2
for j = 1:10
  s2(j).field1 = rand(5,1);
  s2(j).field2 = rand(3,1);
  s2(j).field3 = rand(4,1);
  s2(j).field4 = rand(1,1);
end
s2 = stransp2(s2);

% generate some attributes
a1 = {{vname1, 'gen attr 1', 'gen attr val 1'}, ...
      {vname1, 'gen attr 2', 'gen attr val 2'}, ...
      {'field2', 'field 2 attr 1', 'field 2 attr 1 val'}, ...
      {'field2', 'field 2 attr 2', 'field 2 attr 2 val'}, ...
      {'field3', 'field 3 attr', 'field 3 attr val'}};

a2 = {{vname2, 'gen attr 1', 'gen attr val 1'}, ...
      {vname2, 'gen attr 2', 'gen attr val 2'}, ...
      {'field2', 'field 2 attr 1 set2', 'field 2 attr 1 set 2 val'}, ...
      {'field2', 'field 2 attr 2 set2', 'field 2 attr 2 set 2 val'}, ...
      {'field3', 'field 3 attr', 'field 3 attr val'}};

save vstest2 vname1 vname2 hfile s1 s2 a1 a2


% write as HDF data
% -------------------

fprintf(1, 'calling h4vswrite()...\n');
h4vswrite(hfile, { {vname1, s1, a1}, {vname2, s2, a2}})


% read the data back in 
% ----------------------

clear all
load vstest2

fprintf(1, 'calling h4vsread()...\n');
vlist = h4vsread(hfile);

vname1r = vlist{1}{1};
s1r = vlist{1}{2};
a1r = vlist{1}{3};

vname2r = vlist{2}{1};
s2r = vlist{2}{2};
a2r = vlist{2}{3};


% compare results
% ----------------

[ isequal(s1, s1r);
  isequal(s2, s2r);
  isequal(a1, a1r);
  isequal(a2, a2r) ]'

% clean up
! rm vstest2.hdf vstest2.mat 2> /dev/null

