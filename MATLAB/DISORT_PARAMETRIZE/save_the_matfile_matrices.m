fprintf(1,'ICE   scanang x cprtop x cngwat x cpsize = %4i x %4i x %4i x %4i \n',length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),  length(unique(p.cpsize(booI))))
fprintf(1,'WATER scanang x cprtop x cngwat x cpsize = %4i x %4i x %4i x %4i \n',length(unique(p.scanang(booW))),length(unique(p.cprtop(booW))),length(unique(p.cngwat(booW))),  length(unique(p.cpsize(booW))))
%% make a N dim matrix, one for Ice and one for Water,      Freq x scanang x cprtop x cngwat x cpsize

%% see ~/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/make_profiles.m
%% for cc = 1 : length(ctype)
%%   for dd = 1 : length(water_cpsize)
%%     for qq = 1 : length(cngwat)
%%       for tt = 1 : length(cprtop)
%%         for ss = 1 : length(scanang)

if iVers == 0
  %% Howard says this should work, with matr(a,b,c,d) having a as innermost loop, dummy!!!!
  matrI  = reshape(iaChouFac(bestI),length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI))));
  matrW  = reshape(iaChouFac(bestW),length(unique(p.cpsize(booW))),length(unique(p.cngwat(booW))),length(unique(p.cprtop(booW))),length(unique(p.scanang(booW))));
  indexI = reshape(booI,length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI))));
  indexW = reshape(booW,length(unique(p.cpsize(booW))),length(unique(p.cngwat(booW))),length(unique(p.cprtop(booW))),length(unique(p.scanang(booW))));
elseif iVers == 1
  %% ultimately I am sending in matr(cpsize,cngwat,cprtop,scanang) ie this transpoe etc works, but Howard says matr(a,b,c,d) with d as inner most loop is just plain wierd
  %% BUT NOW I SEE THE LOOPS in make_profiles.m SO I THINK THIS IS RIGHT !!!!! SCANANG is innermost loop, so is first index!!! then for a unknown but CORRECT reason I permute so scanang becomes innermost
  matrI  = reshape(iaChouFac(bestI),length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); matrI = permute(matrI,[ 4 3 2 1]);
  matrW  = reshape(iaChouFac(bestW),length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); matrW = permute(matrW,[ 4 3 2 1]); 
  indexI = reshape(booI,length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); indexI = permute(indexI,[ 4 3 2 1]);
  indexW = reshape(booW,length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); indexW = permute(indexW,[ 4 3 2 1]);
end

cpsizeI = unique(p.cpsize(booI));
cpsizeW = unique(p.cpsize(booW));
cprtop = unique(p.cprtop(booI));
cngwat = unique(p.cngwat(booI));
scanang = unique(p.scanang(booI));
comment = 'see reshape(bestI,length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI)))) and make_profiles.m';
if iReg == 1  
  if iVers == 1
    generic_name = 'generic_605_1655_I_W.mat';
  elseif iVers == 0
    generic_name = 'generic_605_1655_I_W_vers0.mat';
  end
else
  fooA = num2str(floor(f(iaInd(1))),'%04d');
  fooB = num2str(ceil(f(iaInd(end))),'%04d');
  if iVers == 1
    generic_name = ['generic_' fooA '_' fooB '_I_W.mat'];
  elseif iVers == 0
    generic_name = ['generic_' fooA '_' fooB '_I_W_vers0.mat'];
  end
end
saver = ['save ' generic_name ' matrI matrW cpsizeI cpsizeW cprtop cngwat scanang comment indexI indexW iVers'];
eval(saver)

%% [matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_605_1655_I_W.mat',      'generic_605_1655_I_W.bin',      605,2830,iVers=1);  YAY
%% [matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_605_1655_I_W_vers0.mat','generic_605_1655_I_W_vers0.bin',605,2830,iVers=0);  BOO

fprintf(1,'saved %s \n',generic_name)
