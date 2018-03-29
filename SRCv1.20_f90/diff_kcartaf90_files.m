function [nfiles,commonfiles] = diff_kcartaf90_files(thedir1,thedir2);

clc

%{
thedir1 = '/home/sergio/KCARTA/SRCv1.20_f90/';
thedir2 = '/home/sergio/KCARTA/SRCv1.20_f90/WORKS_Feb28_2018_Verified/';
%}

thefiles = dir([thedir1 '/*.f90']);

nfiles = length(thefiles);

fprintf(1,'checking %3i kCARTA f90 files in %s versus those in %s \n',nfiles,thedir1,thedir2);
disp('ret to continue'); pause

if exist('diff_kcarta_f90.txt')
  rmer = ['!/bin/rm diff_kcarta_f90.txt'];
  eval(rmer);
end

echoer = ['!echo "diffing ' thedir1 ' and ' thedir2 ' "  > diff_kcarta_results.txt '];
eval(echoer);

commonfiles = 0;
for ii = 1 : nfiles
  fprintf(1,'diffing  %s vs %s \n',thefiles(ii).name,thedir2);

  echoer = ['!echo " ========>>>>>>>> diffing ' thefiles(ii).name ' " > junk.txt '];
  eval(echoer);
  catter = ['!cat junk.txt >> diff_kcarta_results.txt'];
  eval(catter);
  
  fname = [thedir1 '/' thefiles(ii).name];
  ee = exist([thedir2 '/' thefiles(ii).name]);
  if ee > 0
    commonfiles = commonfiles + 1;
    differ = ['!diff ' fname ' ' thedir2 '/. >& ugh'];
    eval(differ);
    lser = dir('ugh');
    if lser.bytes > 0
      disp(' ')
      fprintf(1,'>>>>>>>>>>> %s is different \n',thefiles(ii).name)
      morer = ['!more ugh']; eval(morer)
      catter = ['!cat ugh >> diff_kcarta_results.txt'];
      eval(catter);      
      %disp('ret to continue'); pause
      pause(0.1)
      disp(' ')
    end
  else
    fprintf(1,'%s is in DIR1 but not in DIR2 \n',thefiles(ii).name);
    disp(' ')
  end
end

eval(['!/bin/rm ugh junk.txt']);
