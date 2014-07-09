
function [a1,a2,a3,a4,a5,a6] = mdeal(a);

if length(a) ~= nargout
  error('input and output lengths must match')
end

a1 = a(1);
if nargout == 1, return, end
a2 = a(2);
if nargout == 2, return, end
a3 = a(3);
if nargout == 3, return, end
a4 = a(4);
if nargout == 4, return, end
a5 = a(5);
if nargout == 5, return, end
a6 = a(6);
if nargout == 6, return, end

