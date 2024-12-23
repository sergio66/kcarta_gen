function [cline] = readnextline(fid, comchar);

% function [cline] = readnextline(fid, comchar);
%
% Read next line not starting with comchar from input file.
%
% Input:
%    fid = [1 x 1] integer file unit ID
%    comchar = [1 x 1] character
%
% Output:
%    cline = [string] character string
%

% Created: 13 July 2010, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skip = 1;
while (skip == 1)
   cline = fgetl(fid);
   if (cline(1) == comchar)
     skip = 1;
   else
     skip = 0;
   end
end

%%% end of program %%%
