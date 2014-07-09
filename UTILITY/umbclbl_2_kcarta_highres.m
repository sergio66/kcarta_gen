%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% almost same as umbclbl_2_kcarta but outputs stuff at 0.0005 cm-1 res
%% so you HAVE to be CAREFUL with the start and stop points 
%% remember instead of starting out at eg 2305.0000, we start at
%%            2305.00 - 0.0005 - 0.0005 = 2304.9990 
%% so need to add on 2 extra points at the beginning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  %% remember instead of starting out at eg 2305.0000, we start at
  %%            2305.00 - 0.0005 - 0.0005 = 2304.9990 
  %% so need to add on 2 extra points at the beginning

fr0(1) = fstart(1) - 2*0.0005;
fr0(2) = fstart(1) - 1*0.0005;
fr0(3:length(fr)+2) = fr;
fr = fr0; clear fr0;

[m,n] = size(k00);
k000 = zeros(m,n+2);
k000(:,3:n+2) = k00;
k00 = k000; clear k000;

f0 = fstart; 

while f0 <= fend

  ii1 = find(abs(fr - (f0 - 2*0.0005)) <= 0.0005/2);

  ii  = ii1 : ii1+49999;

  fprintf(1,'fr(A) = %20.12f \n',fr(ii(1)))
  fprintf(1,'fr(B) = %20.12f \n',fr(ii(length(ii))))

  sfreq = fr(ii(3));
  sfreq = floor(sfreq);
  fname = [fnameout '_' num2str(sfreq) '.dat'];
  fprintf(1,'fname is %s \n',fname);
  fid=fopen(fname,'w','ieee-le');

  fprintf(1,'  start index = %5i  freq(start index) = %8.6f \n',ii(3),sfreq);
  %%%idgas = 1;
  npts = 50000;
  nlay = 100;

  % Write header info
  % header1 = gasid, npts, nlay
  filemark= 4 + 4 + 4;
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,[idgas,npts,nlay],'integer*4');
  fwrite(fid,filemark,'integer*4');

  % header2 = sfreq,fstep
  filemark= 8 + 8;
  sfreq = fr(ii(3));
  df = 0.0005;
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
