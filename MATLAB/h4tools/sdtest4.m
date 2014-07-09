
% test new versions of sdload() and sdsave()
%
% NOTE: it appears that arrays can't have the same
% dimension names!
%

hfile = 'sdtest3.hdf';

A1 = rand(3,2,4);
A2 = rand(2,3);
A3 = rand(3,1);

A1attrs = {{'attr11', 'val11'}, {'attr12', 'val12'}};
A1dims = {'X', 'Y', 'Z'};

A2attrs = {{'attr21', 'val21'}};
A2dims = {'dim1', 'dim2'};

A3attrs = {{'attr31', 'val31'}};
% A3dims = {'dim1', 'dim2'};  % this gives an error
A3dims = {'dim1b', 'dim2b'};  % this works OK

alist = {{'A1', A1, A1attrs, A1dims}, ...
         {'A2', A2, A2attrs, A2dims}, ...
	 {'A3', A3, A3attrs, A3dims} };

fattr = {{'comment', 'a test'}, {'date', date}};

% test sdload with h4sdwrite
% h4sdwrite(hfile, alist, fattr);

% test sdload with sdsave
dat2.A1=A1;      dat2.A2=A2;      dat2.A3=A3; 
atr2.A1=A1attrs; atr2.A2=A2attrs; atr2.A3=A3attrs; 
dim2.A1=A1dims;  dim2.A2=A2dims;  dim2.A3=A3dims;
eval(sprintf('atr2.%s=''%s'';', fattr{1}{1}, fattr{1}{2}))
eval(sprintf('atr2.%s=''%s'';', fattr{2}{1}, fattr{2}{2}))
sdsave(hfile, dat2, atr2, dim2);

% read the data with sdload
[dat, atr, dim] = sdload(hfile);

% check the results
[isequal(A1, dat.A1), isequal(A2, dat.A2), isequal(A3, dat.A3)]
[isequal(A1attrs,atr.A1), isequal(A2attrs,atr.A2), isequal(A3attrs, atr.A3)]
[isequal(fattr{1}{2}, atr.comment), ...
 isequal(fattr{2}{2}, atr.date)]
[isequal(A1dims,dim.A1), isequal(A2dims,dim.A2), isequal(A3dims, dim.A3)]

