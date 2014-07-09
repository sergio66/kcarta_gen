
% test new versions of h4sdwrite() and h4sdread()
%
% a useful reminder from the user's guide:
% 
%   Note that if, for example, the creation of a second dimension
%   named "Altitude" is attempted and the size of the dimension is
%   different from the existing dimension named "Altitude", an error
%   condition will be generated.
% 

hfile = 'sdtest3.hdf';

A1 = rand(3,2,4);
A2 = rand(2,3);
A3 = rand(2,3);  

A1attrs = {{'attr11', 'val11'}, {'attr12', 'val12'}};
A1dims = {'X', 'Y', 'Z'};

A2attrs = {{'attr21', 'val21'}};
A2dims = {'dim1', 'dim2'};

A3attrs = {{'attr31', 'val31'}};
A3dims = {'dim1', 'dim2'};

alist = {{'A1', A1, A1attrs, A1dims}, ...
         {'A2', A2, A2attrs, A2dims}, ...
	 {'A3', A3, A3attrs, A3dims} };

fattr = {{'comment', 'a test'}, {'date', date}};

h4sdwrite(hfile, alist, fattr);

[alist1, fattr1] = h4sdread(hfile);

[isequal(alist, alist1), isequal(fattr, fattr1)]

