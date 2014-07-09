pmin = 800; pmax = 1100;
npts = 101;
dp = (pmax-pmin)/(npts-1);
klayers_Plev = pmin:dp:pmax;
klayers_Plev = fliplr(klayers_Plev);
whos klayers_Plev
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the new klayers from 0-83 km, but in 120 pressure levels
fid = fopen('cbplev_test.f','w');
string = 'c see make_cbplev_test.m';
fprintf(fid,'%s \n',string);
string = 'c 100 pressure levels for using kLAYERS from 0-83 km';
fprintf(fid,'%s \n',string);
fprintf(fid,'%s \n',string);

string = 'C      ------------------------------------------------------------';
fprintf(fid,'%s \n',string);
string = 'C PLEV with the AIRS layer boundary pressure levels (in mb)';
fprintf(fid,'%s \n',string);
string = 'C      ------------------------------------------------------------';
fprintf(fid,'%s \n',string);
string = '       DATA  (PLEV(I), I = 101, 1, -1 )';
fprintf(fid,'%s \n',string);
fprintf(fid,'     $ / \n');

PlevUpFLIP = fliplr(klayers_Plev);
for jj = 1 : 25
  ind = (1:4) + (jj-1)*4;
  data = PlevUpFLIP(ind);
  fprintf(fid,'     $ %12.8e,  %12.8e, %12.8e, %12.8e, \n',data); 
  end

ind = 101;
data = PlevUpFLIP(ind);
fprintf(fid,'     $ %12.8e/ \n',data);

string = '      END';
fprintf(fid,'%s \n',string);

fclose(fid);

