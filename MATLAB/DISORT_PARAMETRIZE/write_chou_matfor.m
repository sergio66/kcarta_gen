function [matrII,matrWW,indexII,indexWW] = write_chou_matfor(fnameIN,fnameOUT,f1,f2);

%% see /home/sergio/HITRAN2UMBCLBL/FORTRAN/mat2for/mat2for.m
%% example write_chou_matfor('generic_605_1655_I_W.mat','generic_605_1655_I_W.bin',605,2830);

dtype = 'ieee-le';

load(fnameIN); % this contains matrI matrW cpsizeI cpsizeW cprtop cngwat scanang comment for f1 to f2

if exist(fnameOUT)
  fnameOUT
  error('file already exists')
end

fourOReight = 8; %% since we are writing out doubles
fourOReight = 4; %% since we are writing out reals

fid = fopen(fnameOUT,'w',dtype);

%% write freq header info
filemark = fourOReight + fourOReight;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[f1 f2],'real*4');
fwrite(fid,filemark,'integer*4');

%% write lenght of arrays header info
filemark = 4 + 4 + 4 + 4;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[length(cpsizeI) length(cngwat) length(cprtop) length(scanang)],'integer*4');
fwrite(fid,filemark,'integer*4');

%% write the ice particle sizes (um)
filemark = fourOReight * length(cpsizeI);
fwrite(fid,filemark,'integer*4');
fwrite(fid,cpsizeI,'real*4');
fwrite(fid,filemark,'integer*4');

%% write the water particle sizes (um)
filemark = fourOReight * length(cpsizeW);
fwrite(fid,filemark,'integer*4');
fwrite(fid,cpsizeW,'real*4');
fwrite(fid,filemark,'integer*4');

%% write the cloud amounts g/m2
filemark = fourOReight * length(cngwat);
fwrite(fid,filemark,'integer*4');
fwrite(fid,cngwat,'real*4');
fwrite(fid,filemark,'integer*4');

%% write the cloud top mb
filemark = fourOReight * length(cprtop);
fwrite(fid,filemark,'integer*4');
fwrite(fid,cprtop,'real*4');
fwrite(fid,filemark,'integer*4');

%% write the scanang
filemark = fourOReight * length(scanang);
fwrite(fid,filemark,'integer*4');
fwrite(fid,scanang,'real*4');
fwrite(fid,filemark,'integer*4');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iVers = 0; %% cpsize cngwat cprtop scanang
iVers = 1; %% scanang cpsize cngwat cprtop 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iVers == 0
  %% now loop and save the matrixI hehehehe [sze cng cpr ang]
  matrII = matrI;
  filemark = fourOReight * length(scanang);
  for dd = 1 : length(cpsizeI)
    for qq = 1 : length(cngwat)
      for tt = 1 : length(cprtop)
        data = matrII(dd,qq,tt,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end
  
  %% now loop and save the matrixW hehehehe [sze cng cpr ang]
  matrWW = matrW;
  filemark = fourOReight * length(scanang); 
  for dd = 1 : length(cpsizeW)
    for qq = 1 : length(cngwat)
      for tt = 1 : length(cprtop)
        data = matrWW(dd,qq,tt,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end

  indexII = indexI;
  %% now loop and save the indexixI hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeI)
      for qq = 1 : length(cngwat)
        data = indexII(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end
  
  indexWW = indexW;
  %% now loop and save the indexixW hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeW)
      for qq = 1 : length(cngwat)
        data = indexWW(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iVers == 1
  matrII = permute(matrI,[4 1 2 3]);
  %% now loop and save the matrixI hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeI)
      for qq = 1 : length(cngwat)
        data = matrII(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end
  
  %% now loop and save the matrixW hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  matrWW = permute(matrW,[4 1 2 3]);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeW)
      for qq = 1 : length(cngwat)
        data = matrWW(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end

  %%% tis is not really needed, just doing it for checking
  indexII = permute(indexI,[4 1 2 3]);
  %% now loop and save the indexixI hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeI)
      for qq = 1 : length(cngwat)
        data = indexII(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end
  
  indexWW = permute(indexW,[4 1 2 3]);
  %% now loop and save the indexixW hehehehe  [ang sze cng cpr]
  filemark = fourOReight * length(cprtop);
  for ss = 1 : length(scanang)
    for dd = 1 : length(cpsizeW)
      for qq = 1 : length(cngwat)
        data = indexWW(ss,dd,qq,:);
        fwrite(fid,filemark,'integer*4');
        fwrite(fid,data,'real*4');
        fwrite(fid,filemark,'integer*4');
      end
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% to debug KCARTA  FUNCTION read_chou_scale_parametrized
%[matrII(2,3,4,5) indexII(2,3,4,5)]
%[matrWW(2,1,3,4) indexWW(2,1,3,4)]
%[matrWW(6,1,7,2) indexWW(6,1,7,2)]
%[matrII(6,1,7,2) indexII(6,1,7,2)]
[matrWW(5,1,6,5) indexWW(5,1,6,5)]
  
fclose(fid);

