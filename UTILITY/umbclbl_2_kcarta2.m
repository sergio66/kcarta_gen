% this file takes data from run6, and outputs it in a style that
% kCARTA would accept, in the *SPECTRA section

%% this template file assumes each 25cm-1 chunk from run6+, is stored in a 
%% separate file; so the code loops over loading in each .mat file and saving 
%% the data to an appropriate kCARTA .dat file
%% note : it expects the data to be saved as "fr" and "k00" in each .mat file

%The header info contains the following integers on one line : 
%$idgas,npts,nlay$. These are the gasID, number of wavenumber points (should 
%equal kMaxPts=10000) and number of layers (which should equal kProfLayer). 
%The next line in the header contains two reals : $sfreq, fstep$ which are the 
%start frequency and  wavenumber step respectively. These numbers should 
%correspond to the corresponding kCompressed file the data is replacing :\\
%\medskip
%{\sf 
%\ttab idgas npts nlay\\
%\ttab sfreq fstep\\
%\ttab }

%After this, the actual data should be stored in layer form, as double 
%precision variable : \\
%\medskip
%{\sf 
%\ttab daAbsLayer1(J),J=1,kMaxPts)\\
%\ttab daAbsLayer2(J),J=1,kMaxPts)\\
%\ttab daAbsLayer3(J),J=1,kMaxPts)\\
%\ttab ...\\
%\ttab daAbsLayerN(J),J=1,kMaxPts)\\
%\ttab }
%where as usual, layer 1 is the ground (bottommost) layer and layer kProfLayer
%is the highest layer.

f0 = 605;
fend = 2805;
while f0 < fend
  f0
  loader = ...
    ['load /taro/s1/sergio/HITRAN2000WATER/CAMEX1/water' num2str(f0) '.mat'];
  eval([loader]);

  sfreq = fr(1);
  sfreq = floor(sfreq);
  fname = ...
     ['/taro/s1/sergio/HITRAN2000WATER/CAMEX1/gas1_' num2str(sfreq) '.dat'];
  fid=fopen(fname,'w','ieee-le');

  idgas = 1;
  npts = 10000;
  nlay = 100;

  % Write header info
  % header1 = gasid, npts, nlay
  filemark= 4 + 4 + 4;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,[idgas,npts,nlay],'integer*4');
  fwrite(fid,filemark,'integer*4');

  % header2 = sfreq,fstep
  filemark= 8 + 8;
  sfreq = fr(1);
  df=0.0025;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,[sfreq,df],'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save the k matrices sequentially
  filemark= 8 * npts;
  for i = 1:nlay
    fwrite(fid,filemark,'integer*4');
    eval(['k' int2str(i) '=k00(i,:);'])
    fwrite(fid,eval(['k' int2str(i)]),'real*8');
    fwrite(fid,filemark,'integer*4');
    end

  fclose(fid);

  f0 = f0 + 25;
  end
