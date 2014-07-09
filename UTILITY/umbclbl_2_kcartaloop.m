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

%fstart   = input('Enter start chunk freq (>=  605 cm-1) : ');
%fend     = input('Enter stop  chunk freq (<= 2805 cm-1) : ');
%fnamein  = input('Enter input (.mat) filename : ');
%fnameout = input('Enter output prefix (.dat) filename : ');
%idgas    = input('Enter gasID : ');

%%% note fnamein,fnameout do not have "blah_" format
% info.fstart = 605; info.fend = 780; info.idgas = 2;
% info.fnamein = 'blah'; info.fnameout = 'kc_std_umbc12';
fstart   = info.fstart;
fend     = info.fend;
fnamein  = info.fnamein;
fnameout = info.fnameout;
idgas    = info.idgas;
clear k*

f0 = fstart;
  loader = ['load ' fnamein num2str(f0) '.mat'];
  eval([loader]);
  whos    %%check for fr (1x10000) and k00 (100 x 10000)
%error('ooh');

clf;
f0 = fstart;
while f0 <= fend
  loader = ['load ' fnamein num2str(f0) '.mat'];
  eval([loader]);

  %%this comes after checking "error (ooh)" above
  fr = w;
  %k00 = dumbc12_shift;
  k00 = d2002;
  plot(fr,k00(1,:)); hold on

  fprintf(1,'f0 = %8.3f    fr(1) = %8.3f \n',f0,fr(1));

  sfreq = fr(1);
  fname = [fnameout '_gas' num2str(idgas) '_' num2str(sfreq) '.dat'];
  fid=fopen(fname,'w','ieee-le');

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
  df = mean(diff(fr));
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
