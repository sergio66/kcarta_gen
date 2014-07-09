% this file takes data from run6+, and outputs it in a style that
% kCARTA would accept, in the *SPECTRA section

%% this template file assumes all necessary 25cm-1 chunks from run6+, are
%% stored in ONE file; so the code loads this master .mat file and loops over
%% saving the data to an appropriate kCARTA .dat file
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

fstart = input('Enter start chunk freq (>=  605 cm-1) : ');
fend   = input('Enter stop  chunk freq (<= 2805 cm-1) : ');
fnamein  = input('Enter input (.mat) filename : ');
fnameout = input('Enter output prefix (.dat) filename : ');
idgas    = input('Enter gasID : ');

loader = ['load ' fnamein '.mat'];
eval([loader]);

aa = exist('fr'); 
if aa ~= 1
  error('need to have your "fr" array in the .mat file');
  end
aa = exist('k00'); 
if aa ~= 1
  error('need to have your "k00" matrix in the .mat file');
  end

[mm,nn] = size(k00);
if mm > nn
  k00 = k00';
  end

[mm,nn] = size(k00);
if mm > 100
  error('hmm, edit file becuz kCARTA usually cannot handle > 100 layers');
  end
if mm < 100
  k000 = zeros(100,nn);
  k000(100-mm+1:100,:) = k00;
  k00 = k000; clear k000;
  end

f0 = fstart;
while f0 <= fend

  ii1 = find(abs(fr - f0) <= 0.0025/2);
  ii  = ii1 : ii1+9999;

  sfreq = fr(ii(1));
  sfreq = floor(sfreq);
  fname = [fnameout '_' num2str(sfreq) '.dat'];
  fprintf(1,'fname is %s \n',fname);
  fid=fopen(fname,'w','ieee-le');

  fprintf(1,'  start index = %5i  freq(start index) = %8.6f \n',ii(1),sfreq);
  %%%idgas = 1;
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
  sfreq = fr(ii(1));
  df=0.0025;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,[sfreq,df],'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save the k matrices sequentially
  filemark= 8 * npts;
  for i = 1:nlay
    fwrite(fid,filemark,'integer*4');
    eval(['k' int2str(i) '=k00(i,ii);'])
    fwrite(fid,eval(['k' int2str(i)]),'real*8');
    fwrite(fid,filemark,'integer*4');
    end

  fclose(fid);

  f0 = f0 + 25;
  end
