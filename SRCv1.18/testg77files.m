%%% tests everything except rtp_interface.f and *unused.f

g77test = '/usr/bin/g77 -fcase-lower -ffixed-line-length-120 -c ';
thedir = dir('/home/sergio/kcartaV118/SRCv1.18/*.f');
nfiles = length(thedir)

for ii = 1:nfiles
  %% pretend every month has 31 days and we start on 01/01/2000
  lala = datevec(thedir(ii).date,31);
  yy(ii) = lala(1); mm(ii) = lala(2); dd(ii) = lala(3);
  hh(ii) = lala(4); nn(ii) = lala(5); ss(ii) = lala(4);
  indextime(ii) = (yy(ii)-2000)*372*24*60*60 + (mm(ii)-1)*31*24*60*60 + (dd(ii)-1)*24*60*60 + ...
                  (hh(ii)-1)*60*60 + (nn(ii)-1)*60 + ss(ii);
end

[Y,I] = sort(indextime,'descend');
for ii = 1 : nfiles
  fname = thedir(I(ii)).name;
  if (length(strfind('rtp_interface.f',fname)) == 0 & ...
      length(strfind('unused.f',fname)) == 0 & ...
      length(strfind('rtpdemo.f',fname)) == 0)
    fprintf(1,'g77 testing %s \n',fname);
    g77er = ['!' g77test ' ' fname]; eval(g77er);
    disp('     tested file ... hit return')
    pause
    end
  end
