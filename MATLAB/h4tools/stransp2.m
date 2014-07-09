
function s1 = stransp2(s2)

% function s1 = stransp2(s2)
%
% transpost an array of structures into a structure of arrays
%
% example: if s2 is an array of 6 structures, each with fields
% 
%   a: [4x1]
%   b: [3x1]
%   c: [2x1]
%
% then s1 is a structure
% 
%   a: [4x6]
%   b: [3x6]
%   c: [2x6]
% 
% The second dimension of the fields of s2 must all be 1.
%

flist = fieldnames(s2);

% get number of records
n = length(s2);

% check the second dimension of the first field of s1(1)
[m,n2] = size(getfield(s2(1), flist{1}));
if n2 ~= 1
  error('the second dimension of the fields of s2 must all be 1');
end

for i = 1 : n  % loop on records
  for j = 1 : length(flist)  % loop on fields

     field = flist{j};

     % want s1.field(:,i) = s2(i).field

     eval(sprintf('s1.%s(:,%d) = s2(%d).%s;', field, i, i, field));
  end
end

