%%% tests everything except rtp_interface.f and *unused.f

g77test = '/usr/bin/g77 -fcase-lower -ffixed-line-length-120 -c ';
thedir = dir('/home/sergio/kcartaV114/SRCv1.10/*.f');
nfiles = length(thedir);
for ii = 1 : nfiles
  fname = thedir(ii).name;
  if (length(strfind('rtp_interface.f',fname)) == 0 & ...
      length(strfind('unused.f',fname)) == 0)
    fprintf(1,'g77 testing %s \n',fname);
    g77er = ['!' g77test ' ' fname]; eval(g77er);
    disp('     tested file ... hit return')
    pause(1)
    end
  end
