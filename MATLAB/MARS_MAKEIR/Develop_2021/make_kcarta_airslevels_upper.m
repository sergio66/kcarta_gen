XYZ = fliplr(plevsMarsUpper);

fid = fopen('/home/sergio/KCARTA/INCLUDE/MARS_database_params_2021/airslevels_upper.param_mars','w');

str = 'c see /home/sergio/HITRAN2UMBCLBL/MARS_MAKEIR/Develop_2021/produce_kcarta_paramfiles_2021.m';
fprintf(fid,'%s \n\n',str);

str = 'c 101 pressure levels for the DEFINITION of the kCARTA MARS database';
fprintf(fid,'%s \n',str);
str = 'c Sergio Machado UMBC May 2021';
fprintf(fid,'%s \n\n',str);
str = 'c this is for the upper atm ';
fprintf(fid,'%s \n',str);

str = ' ';
fprintf(fid,'%s \n',str);
str = '       INTEGER iXPlanet6';
fprintf(fid,'%s \n',str);
str = '       DATA iXPlanet6 /04/';
fprintf(fid,'%s \n',str);
str = ' ';
fprintf(fid,'%s \n',str);

str = '       REAL DATABASELEV(101)';
fprintf(fid,'%s \n',str);
str = '       INTEGER IPLEV';
fprintf(fid,'%s \n\n',str);

str = 'C      -----------------------------------------------------------------';
fprintf(fid,'%s \n',str);
str = 'C      PLEV with the AIRS layer boundary pressure levels (in mb)';
fprintf(fid,'%s \n',str);
str = 'C      -----------------------------------------------------------------';
fprintf(fid,'%s \n',str);
str = '       DATA  ( DATABASELEV(IPLEV), IPLEV = 101, 1, -1 )';
fprintf(fid,'%s \n',str);
str = '     $ /';
fprintf(fid,'%s \n',str);
for ii = 1:25
  index = (1:4) + (ii-1)*4;
  fprintf(fid,'     $ %10.6e, %10.6e, %10.6e, %10.6e, \n',XYZ(index));
end
index = 101;
fprintf(fid,'     $ %10.6e \n',XYZ(index));
str = '     $ /';
fprintf(fid,'%s \n\n',str);

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
