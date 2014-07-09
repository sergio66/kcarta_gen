% this file takes data from run6, and outputs it in a style that
% kCARTA would accept, in the *SPECTRA section

%% copied from umbclbl_2_kcarta.m .. used for variable spacing

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

%%% put gas2 to trick things

f1 = 4500;
f2 = 22000;
nbox = 5; 
pointsPerChunk = 10000; 

fend = f2;
f0 = f1;
topts = runXtopts_params_smart(f0);  
dv = topts.ffin*nbox*pointsPerChunk; 

gases = [1 2 3 4 5 6 7];

gasid1 = 1;
gasid3 = 3;

dvx = -1;
f0 = f1;
while f0 < fend
  topts = runXtopts_params_smart(f0); 
  dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = f0 + dv; 
  if abs(dv-dvx) > 1
    disp('----------------------------------------------------------')
    end
  dvx = dv;
  fprintf(1,'%10.2f  %10.6f  %10.6f  %10.6f\n',f0,fmax,topts.ffin,dv);
  f0 = f0 + dv;
  end

wall  = [];
dallA = [];  %% absorption from HITRAN lines
dallC = [];  %% CKD continuum absoprtion
dallO = [];  %% O3 chappius absoprtion
dallX = [];  %% sum total absorption

cd /home/sergio/SPECTRA
f0 = f1;
while f0 < fend
  topts = runXtopts_params_smart(f0); 
  dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = f0 + dv; 

  toptswc.ffin = topts.ffin;
  toptswc.CKD  = 1;
  toptswc.nbox = 5;

  toptswc.ffin = topts.ffin * 5;
  toptswc.CKD  = 1;
  toptswc.nbox = 1;

  fnameIN = ['/carrot/s1/sergio/RUN8_VISDATABASE/trp' num2str(f0) '.mat'];

  if (exist(fnameIN))
    %%% water continuum
    clear w dcon

    fip = ['IPFILES/trp_g' num2str(gasid1)];
    figure(1)
    [w,dcon] = run7watercontinuum(gasid1,f0,fmax,fip,toptswc);  
    
    loader = ['load ' fnameIN];
    eval([loader]);

    %%% chappius
    clear dozone

    fip = ['IPFILES/trp_g' num2str(gasid)];
    qamt = load(fip);
    qamt = qamt(:,5)*6.023e26;
    figure(1)
    for ll = 1 : 100
      dozone(:,ll) = chappius_ozone_vis(wvis0,qamt(ll));
      end
    dozone = dozone';

    fr = wvis0;
    k00 = dvis0 + dcon + dozone;

    figure(2);
    wall = [wall fr];
    dallA = [dallA sum(dvis0)];
    dallC = [dallC sum(dcon)];
    dallO = [dallO sum(dozone)];
    dallX = [dallX sum(k00)];

    plot(10000./wall,exp(-dallA),10000./wall,exp(-dallC),...
         10000./wall,exp(-dallX),'r'); pause(0.1)

    sfreq = fr(1);
    sfreq = floor(sfreq);
    fname = ...
       ['/carrot/s1/sergio/RUN8_VISDATABASE/wo_gas2_' num2str(sfreq) '.dat'];
    fid=fopen(fname,'w','ieee-le');

    idgas = 2;
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
    fprintf(1,'%s %10.2f  %10.6f  %10.6f\n',fnameIN,sfreq,df,dv);

  else
    fprintf(1,'%s NOT FOUND \n',fnameIN);
    end

  f0 = f0 + dv;
  end

modis_channels = [470 550 659 865 1240 1640 2130]*1e-3;
modis_channels = 10000./modis_channels;
figure(1);clf
plot(wall,exp(-dallA),wall,exp(-dallC),...
     wall,exp(-dallX),'r',...
     modis_channels,ones(size(modis_channels)),'kx'); 

figure(2);
plot(10000./wall,exp(-dallA),10000./wall,exp(-dallC),...
     10000./wall,exp(-dallX),'r',...
     10000./modis_channels,ones(size(modis_channels)),'kx'); 

