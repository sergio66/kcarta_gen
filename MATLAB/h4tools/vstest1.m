
%
% test of mat2vsfid and vsfid2mat, file attributes
%
% note: to dump headers: hdp dumpvd -h vstest.hdf

% generate some test data

% set hdf file and vs names
vname1 = 'test vs #1';
vname2 = 'test vs #2';
hfile = 'vstest1.hdf';
! rm vstest1.hdf 2> /dev/null

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

save vstest vname1 vname2 hfile s1 s2 a1 a2

% test write of HDF data
% -----------------------

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

% [refs,count] = hdfv('lone',file_id,100)

fprintf(1, 'calling mat2vsfid()...\n');
mat2vsfid(file_id, vname1, s1, a1);

fprintf(1, 'calling mat2vsfid()...\n');
mat2vsfid(file_id, vname2, s2, a2);

% [refs,count] = hdfv('lone',file_id,100)

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


% read the data back in 
% ----------------------

clear all
load vstest

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

fprintf(1, 'calling vsfid2mat()...\n');
[sr1, ar1] = vsfid2mat(file_id, vname1);

fprintf(1, 'calling vsfid2mat()...\n');
[sr2, ar2] = vsfid2mat(file_id, vname2);

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


% compare results
% ----------------

[ isequal(s1, sr1);
  isequal(s2, sr2);
  isequal(a1, ar1);
  isequal(a2, ar2) ]'

% clean up
! rm vstest1.hdf 2> /dev/null

