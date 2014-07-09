
function s2 = stransp1(s1)

% function s2 = stransp1(s1)
%
% transpost a structure of arrays into an array of structures
%
% example: if s1 is
% 
%   a: [4x6]
%   b: [3x6]
%   c: [2x6]
% 
% then s2 is an array of 6 structures, each with fields
% 
%   a: [4x1]
%   b: [3x1]
%   c: [2x1]
%
% The second dimension of the fields of s1 must all be the same.
%

flist = fieldnames(s1);

% get second dimension of the first s1 field
[m,n] = size(getfield(s1, flist{1}));

for i = 1 : n  % loop on records
  for j = 1 : length(flist)  % loop on fields

     field = flist{j};

     % want s2(i).field = s1.field(:,i)

     eval(sprintf('s2(%d).%s = s1.%s(:,%d);', i, field, field, i));
  end
end

