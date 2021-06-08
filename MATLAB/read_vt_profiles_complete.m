function [p,T,QV,Tnlte,bandcenter,iso,vibinfo] = read_vt_profiles_complete(fname);

%% same as "read_vt_profiles" except it also gives vibrational state info : LEVEL, IDAFGL, ENERGRY etc
%% /home/sergio/KCARTA/NONLTE_PRODUCTION/VT_48PROFILES_120_400ppmv_v121_H16_NLTEH16_Apr2021_NewNLTEProfiles_SAVETHISGOOD/sergio_merge/vt1_s0.prf
%! /home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_D
%! The vib temps are given for
%! the following number of CO2 bands : iCount =                           47
%! vib. states at following number of pressure levels : iNV =            103
%!   IL    MOL   IDMOL  ISO   IDISO  LEVEL   IDAFGL  ENERGY(cm-1)
%!    1 CO2       2    626      1     1101      2     667.38000
%!    2 CO2       2    626      1    10002      3    1285.40906
%!    3 CO2       2    626      1     2201      4    1335.13196
%!    4 CO2       2    626      1    10001      5    1388.18506
%!    5 CO2       2    626      1    11102      6    1932.46997
%
% from TAPE4 examples
%------ CO2 VIBRATIONAL STATE DATA (NUM(I),ISO(I),ID(I),EE(I),NDG(I),I=1,N)
%    1     1 '00001'      0.000     1
%    2     1 '10002'   1285.410     1

%% [p,T,QV,Tnlte,bandcenter,iso,vibinfo] = read_vt_profiles_complete(fname);
%% reads in the NLTE profiles produced/given by Manul Lopez Puertas
%% these files can be ingested by kCARTA or GENLN2

%% input
%%   fname eg '/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/...
%%             VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt5_s0.prf';
%%         eg '/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/
%%             nlte_1_1_1_3_sol_0.genln2';
%%
%% output
%%     p   = n x 1 array of pressures
%%     T   = n x 1 array of kinetic temps (LTE)
%%    QV   = n x 6 matrix of vib partition functions (for 6 isotopes)
%% Tnlte   = n x 30 matrix of nlte vibrational temps
%%   bc    = 1 x 30 matrix of band centers
%% vibinfo = structure of ISO, IUSGQ,energy 

p = [];
T = [];
QV  = [];
Tnlte = [];

ee = exist(fname);
if ee == 0
  fname
  error('file not found')
end

fid = fopen(fname,'r');

%% read header comments
GoOn = 1;
iCnt = 0;

while GoOn > 0
  strAll = fgetl(fid);
  str1   = strAll(1);
  str56   = strAll(5:6);
  if strcmp(str1,'!') & ~strcmp(str56,'IL')
    GoOn = 1;
    iCnt = iCnt + 1;
  else
    GoOn = -1;
  end
end

%% read bands info
GoOn = 1;
iCntB = 0;

while GoOn > 0
  strAll = fgetl(fid);
  str1   = strAll(1);
  str56   = strAll(5:6);
  if ~strcmp(str1,'*') & ~strcmp(str56,'  ')
    GoOn = 1;
    iCntB = iCntB + 1;
    vibinfo.iso3(iCntB)    = str2num(strAll(23:25));
    vibinfo.iso(iCntB)     = str2num(strAll(32:32));    
    vibinfo.iusgq(iCntB,:) =        (strAll(37:41));    
    vibinfo.hitID(iCntB)   = str2num(strAll(47:48));    
    vibinfo.energy(iCntB)  = str2num(strAll(53:62));
  else
    GoOn = -1;
  end
end

strAll = fgetl(fid);   %% this should be *** 1 47 1 103

iDebug = -1;

if iDebug > 0
  fprintf(1,'have read in  %3i lines of comments starting with ! \n',iCnt+1)
  fprintf(1,'First useful string = %s \n',strAll);
end
[unk,junkI1,numBands,junk2,numLevs] = strread(strAll,'%s%d%d%d%d', 'delimiter',' ');

strAll = fgetl(fid);
[iGasID,numIso] = strread(strAll,'%d%d', 'delimiter',' ');

%if iDebug > 0
fprintf(1,'there are %3i bands, %3i levels, %3i isos in %s \n',numBands,numLevs,numIso,fname);
%end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% so now we know how many blocks of 5 to read, and lines containing less 
%% than 5 elements 
iNumFullBlocks = floor(numLevs/5);
iLeftover = numLevs - 5*iNumFullBlocks;
%%%%%%%%%%%%%%%%%%%%%%%%%

%% get p and T
p = get_vt_blocks(iNumFullBlocks,iLeftover,fid);
strAll = fgetl(fid);
T = get_vt_blocks(iNumFullBlocks,iLeftover,fid);

%if strfind(fname,'/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/')
%  %% this is one of the US STD Toffsets I have made
%  semilogy(T,p); set(gca,'ydir','reverse');
%end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the QV

strAll = fgetl(fid);
for ii = 1 : numIso
  strAll = fgetl(fid);
  QV(:,ii) = get_vt_blocks(iNumFullBlocks,iLeftover,fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% reading Tnlte for  1 bandcenter = 2349.1433 
%% reading Tnlte for  2 bandcenter = 2283.4880 
%% reading Tnlte for  3 bandcenter = 2332.1130 
%% reading Tnlte for  4 bandcenter = 2340.0140 
%% reading Tnlte for  5 bandcenter = 3004.0122 
%% reading Tnlte for  6 bandcenter = 2920.2390 
%% reading Tnlte for  7 bandcenter = 2982.1120 
%% reading Tnlte for  8 bandcenter = 3659.2730 
%% reading Tnlte for  9 bandcenter = 3557.3122 
%% reading Tnlte for 10 bandcenter = 3612.8420 
%% reading Tnlte for 11 bandcenter = 3527.7380 
%% reading Tnlte for 12 bandcenter = 3714.7830 
%% reading Tnlte for 13 bandcenter = 3632.9100 
%% reading Tnlte for 14 bandcenter = 2797.1360 
%% reading Tnlte for 15 bandcenter = 3714.7830 
%% reading Tnlte for 16 bandcenter = 3612.8420 
%% reading Tnlte for 17 bandcenter = 4247.7050 
%% reading Tnlte for 18 bandcenter = 4314.9130 
%% reading Tnlte for 19 bandcenter = 4390.6290 
%% reading Tnlte for 20 bandcenter = 667.3801 
%% reading Tnlte for 21 bandcenter = 1285.4090 
%% reading Tnlte for 22 bandcenter = 1335.1320 
%% reading Tnlte for 23 bandcenter = 1388.1850 
%% reading Tnlte for 24 bandcenter = 1932.4702 
%% reading Tnlte for 25 bandcenter = 2003.2463 
%% reading Tnlte for 26 bandcenter = 2076.8560 
%% reading Tnlte for 27 bandcenter = 648.4784 
%% reading Tnlte for 28 bandcenter = 1265.8282 
%% reading Tnlte for 29 bandcenter = 1297.2640 
%% reading Tnlte for 30 bandcenter = 1370.0630 
%% reading Tnlte for 31 bandcenter = 662.3734 
%% reading Tnlte for 32 bandcenter = 1259.4260 
%% reading Tnlte for 33 bandcenter = 1325.1410 
%% reading Tnlte for 34 bandcenter = 1365.8440 
%% reading Tnlte for 35 bandcenter = 1901.7370 
%% reading Tnlte for 36 bandcenter = 664.7289 
%% reading Tnlte for 37 bandcenter = 1272.2870 
%% reading Tnlte for 38 bandcenter = 1329.8430 
%% reading Tnlte for 39 bandcenter = 1376.0274 
%% reading Tnlte for 40 bandcenter = 1916.6950 

bandCntPlot = 2283;  %% weaker    Band
bandCntPlot = 2332;  %% weaker    Band
bandCntPlot = 2340;  %% weaker    Band
bandCntPlot = 2350;  %% strongest Band

%% get the vibtemps
strAll = fgetl(fid);
for ii = 1 : numBands
  strAll = fgetl(fid);
  xarr = strread(strAll);
  %xarr = xarr(1:end-1);
  bandcenter(ii) = xarr(5);
  iso(ii)        = xarr(3);
  if iDebug > 0
    fprintf(1,'reading Tnlte for %2i/%2i bandcenter = %8.4f \n',ii,numBands,bandcenter(ii))
  end
  Tnlte(:,ii) = get_vt_blocks(iNumFullBlocks,iLeftover,fid);
  if (abs(bandcenter(ii)-bandCntPlot) < 1)
    iB = ii;
    if iDebug > 0
      fprintf(1,'%3i %3i %3i \n',xarr(2:4))
    end
    if strfind(fname,'/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/')
      %% this is one of the "US STD Toffsets" I have made
      semilogy(T,p,'cs-',Tnlte(:,ii),p,'ms-'); set(gca,'ydir','reverse');
    else
      semilogy(T,p,'bo-',Tnlte(:,ii),p,'ro-'); set(gca,'ydir','reverse');
    end
  end
end
fclose(fid);

%whos p T QV Tnlte bandcenter

%QV = QV(:,1);
%Tnlte = Tnlte(:,iB);
