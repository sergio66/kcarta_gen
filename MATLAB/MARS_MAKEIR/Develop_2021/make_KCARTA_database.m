fid = fopen('/home/sergio/KCARTA/INCLUDE/MARS_database_params_2021/KCARTA_database.param_mars','w');

str = 'c see /home/sergio/HITRAN2UMBCLBL/MARS_MAKEIR/Develop_2021/produce_kcarta_paramfiles_2021.m';
fprintf(fid,'%s \n',str);
str = 'c Sergio Machado UMBC May 2021';
fprintf(fid,'%s \n\n',str);

str = ' ';
fprintf(fid,'%s \n',str);
str = '       INTEGER iXPlanet0';
fprintf(fid,'%s \n',str);
str = '       DATA iXPlanet0 /04/';
fprintf(fid,'%s \n',str);
str = ' ';
fprintf(fid,'%s \n',str);

str = '      REAL PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)  ';
fprintf(fid,'%s \n',str);
str = '      REAL PAVG_KCARTADATABASE_AIRS(kMaxLayer)    ';
fprintf(fid,'%s \n\n',str);

str = '      DATA  (PLEV_KCARTADATABASE_AIRS(IPLEV),IPLEV=     1,  101,   1)';
fprintf(fid,'%s \n',str);
str = '     $ /';
fprintf(fid,'%s \n',str);
for ii = 1:25
  index = (1:4) + (ii-1)*4;
  fprintf(fid,'     $ %10.4f, %10.4f, %10.4f, %10.4f, \n',Y(index));
end
index = 101;
fprintf(fid,'     $ %10.4f \n',Y(index));
str = '     $ /';
fprintf(fid,'%s \n\n',str);

str = '      DATA  (PAVG_KCARTADATABASE_AIRS(IPLEV),IPLEV=     1,  100,   1)';
fprintf(fid,'%s \n',str);
str = '     $ /';
fprintf(fid,'%s \n',str);
for ii = 1:25
  index = (1:4) + (ii-1)*4;
  if ii ~= 25
    fprintf(fid,'     $ %10.4f, %10.4f, %10.4f, %10.4f, \n',Y(index));
  else
    fprintf(fid,'     $ %10.4f, %10.4f, %10.4f, %10.4f \n',Yav(index));
  end
end
str = '     $ /';
fprintf(fid,'%s \n\n',str);

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
