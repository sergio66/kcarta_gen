cper = ['!/bin/cp -a ../../SRCv1.18/*.f .'];
eval(cper)

thedir = dir('*.f');

for ii = 1 : length(thedir)
  fnameIN  =  thedir(ii).name;
  fprintf(1,'%3i %s \n',ii,fnameIN)
end
disp('ret to continue'); pause

iiX = 11;
iiX = 1 : length(thedir);
for iii = 1 : length(iiX)
  ii = iiX(iii);
  fnameIN  =  thedir(ii).name;
  xfnameIN =  ['x' thedir(ii).name];    
  fnameOUT1 = ['TEMPF90/' fnameIN '90'];
  fnameOUT2 = ['../' fnameIN '90'];

  %% first get rid of tabs
  tabber = ['!sed -e "s/\t/      /g" ' fnameIN ' > ' xfnameIN];
  eval(tabber)
  mver = ['!/bin/mv ' xfnameIN ' ' fnameIN];
  eval(mver)

%{
  %% OLD WAY : not working so well
  %% then shorten the lines
  fid = fopen('junkIN','w');
  fprintf(fid,'%s\n',fnameIN);
  fclose(fid);
  shorter = ['!f77_nchar_to_120char < junkIN'];
  eval(shorter)
  mver = ['!/bin/mv ' xfnameIN ' ' fnameIN];
  eval(mver)

  %% run this code to convert f77 file to f90 file
  fprintf(1,'processing %3i of %3i %s \n',ii,length(thedir),fnameIN);
  f77tof90er = ['!../../F77toF90/f2f90convert/f2f90 ' fnameIN(1:length(fnameIN)-2) ' 2 10 t f >& ugh'];
  eval(f77tof90er)
%}

  %% NEW WAY
  %% run this code to convert f77 file to f90 file
  fprintf(1,'processing %3i of %3i %s \n',ii,length(thedir),fnameIN);
  fid = fopen('junkIN','w');
  fprintf(fid,'%s\n',fnameIN);
  fclose(fid);  
  f77tof90er = ['!../../F77toF90/TO_F90/to_f90 < junkIN >& ugh'];
  eval(f77tof90er)

  %% move the f90 file to temporary dir
  mver = ['!mv ' fnameIN(1:length(fnameIN)-2) '.f90 ' fnameOUT1];
  eval(mver)

  %% sed the file in the temp dir to change include file names from *.parm to *param.param
  sedder = ['!sed -e "s/kcarta.param/kcartaparam.f90/g"  -e "s/scatter.param/scatterparam.f90/g" ' ];
  sedder = [sedder ' -e "s/airsheights.param/airsheightsparam.f90/g" '];
  sedder = [sedder ' -e "s/airslevels.param/airslevelsparam.f90/g" '];
  sedder = [sedder ' -e "s/airslevelheights.param/airslevelheightsparam.f90/g" '];  
  sedder = [sedder ' -e "s/airsheights_uppers.param/airsheights_uppersparam.f90/g" '];
  sedder = [sedder ' -e "s/airslevels_upper.param/airslevels_upperparam.f90/g" '];
  sedder = [sedder ' -e "s/airslevelheights_upper.param/airslevelheights_upperparam.f90/g" '];  
  sedder = [sedder ' -e "s/kcarta_jpl.param/kcarta_jplparam.f90/g" '];
  sedder = [sedder ' -e "s/KCARTA_database.param/KCARTA_databaseparam.f90/g" '];
  sedder = [sedder ' -e "s/gasIDname.param/gasIDnameparam.f90/g" '];        
  sedder = [sedder ' -e "s/include''../include ''../g" ' fnameOUT1 ' > ' fnameOUT2];
  eval(sedder)
  disp(' ')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

