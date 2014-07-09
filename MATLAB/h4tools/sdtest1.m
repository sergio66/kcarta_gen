
%
% test of h4sgread.m and h4sgwrite.m
%
% NOTE: max number of SDS's seems to be 5000!


% generate some test data
% ------------------------

hfile = 'sdtest1.hdf';
! rm sdtest1.hdf 2> /dev/null

a.attr1 = 'A val';
a.attr2 = 'B val';

clear s

% for j = 1:1250
for j = 1:200
  s(j).field1 = rand(100,2);
  s(j).field2 = rand(100,2,2);
  s(j).field3 = rand(200,1);
  s(j).field4 = rand(100,1);
end

save sdtest1 hfile a s


% write the data
% ---------------

fprintf(1, 'calling h4sgwrite()...\n');
h4sgwrite(hfile, s, a);


% read the data
% --------------

clear all
load sdtest1

fprintf(1, 'calling h4sgread()...\n');
[dstr, astr] = h4sgread(hfile);

fprintf(1, 'comparing results...\n');
isequal(s, dstr)

! rm sdtest1.hdf sdtest1.mat 2> /dev/null

