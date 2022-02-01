%% 1 * 89 * 10000 pts clr sky rad file                          is about 3566244 bytes
%% 3 * 89 * 10000 pts 1 cld, 1 clr, 1 wgt sum rad file          is about 10566244 bytes
%% 5 * 89 * 10000 pts 2 cld, 1 cld12, 1 clr, 1 wgt sum rad file is about 17566244 bytes

homedir = pwd;
if ~strcmp(homedir(end-11:end),'MANYPROFILES')
  homedir
  error('should be in kcartaV118/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES')
end

iOne = 3566244;
iOne = 3560000;

iCnt1 = 0;
iCnt3 = 0;
iCnt5 = 0;

iYes = input('have you moved the rad.datXYZ_CLD files away??? (-1/+1) : ');
if iYes < 0
  mver = ['!/bin/mv JUNK/*_CLD JUNK/CLD/.'];
  eval(mver);
end

thedir = dir(['JUNK/rad.dat*']);
for ii = 1 : length(thedir);
  sz = thedir(ii).bytes;
  if sz > 3.1*iOne &  sz < 5*iOne
    iCnt5 = iCnt5 + 1;
    fname = ['JUNK/' thedir(ii).name];
    fprintf(1,'bad5 %s is of size %8i \n',fname,sz);
    rmerx = ['!/bin/rm ' fname];
    eval(rmerx);    
  elseif sz > 1.1*iOne &  sz < 3*iOne
    iCnt3 = iCnt3 + 1;
    fname = ['JUNK/' thedir(ii).name];
    fprintf(1,'bad3 %s is of size %8i \n',fname,sz);
    rmerx = ['!/bin/rm ' fname];
    eval(rmerx);
  elseif sz > 0 &  sz < iOne
    iCnt1 = iCnt1 + 1;
    fname = ['JUNK/' thedir(ii).name];
    fprintf(1,'bad1 %s is of size %8i \n',fname,sz);
    rmerx = ['!/bin/rm ' fname];
    eval(rmerx);    
  end
end
fprintf(1,'%6i,%6i,%6i out of %6i seem bad \n',iCnt1,iCnt3,iCnt5,length(thedir))
