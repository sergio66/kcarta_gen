%%% tests everything except rtp_interface.f and *unused.f

clc
thedir1 = '/home/sergio/KCARTA/SRCv1.16/';

%% had some issues with rtp file when checking in layer thicknesses
thedir2 = '/home/sergio/KCARTA/SRCv1.17/WORKS_March31_2012/';

%% these two are same
thedir2 = '/home/sergio/KCARTA/SRCv1.17/WORKS_April30_2012/';
thedir2 = '/home/sergio/KCARTA/SRCv1.16/WORKS_April30_2012/';
thedir2 = '/home/sergio/KCARTA/SRCv1.17/WORKS_May01_2012/';

thedir2 = '/home/sergio/KCARTA/SRCv1.18/';

thefiles = dir([thedir1 '/*.f']);

nfiles = length(thefiles);

fprintf(1,'checking %3i kCARTA f77 files versus those in %s \n',nfiles,thedir2);

for ii = 1 : nfiles
  fprintf(1,'diffing  %s vs %s \n',thefiles(ii).name,thedir2);
  fname = [thedir1 thefiles(ii).name];
  ee = exist([thedir2 thefiles(ii).name]);
  if ee > 0
    differ = ['!diff ' fname ' ' thedir2 '/. >& ugh'];
    eval(differ);
    lser = dir('ugh');
    if lser.bytes > 0
      disp(' ')
      fprintf(1,'>>>>>>>>>>> %s is different \n',thefiles(ii).name)
      morer = ['!more ugh']; eval(morer)
      disp('ret to continue'); pause
      disp(' ')
    end
  else
    fprintf(1,'%s is in DIR1 but not in DIR2 \n',thefiles(ii).name);
    disp(' ')
  end
end
