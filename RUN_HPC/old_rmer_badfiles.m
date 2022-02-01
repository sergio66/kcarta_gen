%% 1 * 89 * 10000 pts clr sky rad file                          is about 3566244 bytes
%% 3 * 89 * 10000 pts 1 cld, 1 clr, 1 wgt sum rad file          is about 10566244 bytes
%% 5 * 89 * 10000 pts 2 cld, 1 cld12, 1 clr, 1 wgt sum rad file is about 17566244 bytes

homedir = pwd;
if ~strcmp(homedir(end-11:end),'MANYPROFILES')
  homedir
  error('should be in kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES')
end

error('try rmer_badfiles2.m ???')

iCnt = 0;
thedir = dir(['JUNK/rad.dat*']);
for ii = 1 : length(thedir);
  sz = thedir(ii).bytes;
  if sz < 3000000
    iCnt = iCnt + 1;
    fname = ['JUNK/' thedir(ii).name];
    fprintf(1,'%s is of size %8i \n',fname,sz);
    rmerx = ['!/bin/rm ' fname];
    eval(rmerx);
  end
end
fprintf(1,'%6i out of %6i seem bad \n',iCnt,length(thedir))
