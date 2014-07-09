cd /home/sergio/SPECTRA

%%% this dumps out the weak backgnd lines, that do not use linemix
topts.linemixloop = -1;

gasID = 2;
profname = '/home/sergio/SPECTRA/IPFILES/std_co2';
iFr = 2205:25:2505;
%%iFr = 2305:25:2505;
co2 = load(profname);
dir = '/carrot/s1/sergio/AIRSCO2/LTE_WEAK_STD_OPTDEPTH/gas2_';

for ii = 1 : length(iFr)
%for ii = 1 : 1
  clear fr k00

  fmin = iFr(ii);
  fmax = fmin + 25;
  [fr,k00]=run7co2(gasID,fmin,fmax,profname,topts);

  sfreq = fr(1);
  sfreq = floor(sfreq);
  fname = ...
     [dir num2str(sfreq) '.dat'];
  fid=fopen(fname,'w','ieee-le');

  idgas = gasID;
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
  df = (fmax-fmin)/length(fr);
  fwrite(fid,filemark,'integer*4');
  fwrite(fid,[sfreq,df],'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save p,pp,t,q
  filemark= 8 * nlay;
  fwrite(fid,filemark,'integer*4');
  eval(['p=co2(:,2);'])
  fwrite(fid,p,'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save p,pp,t,q
  filemark= 8 * nlay;
  fwrite(fid,filemark,'integer*4');
  eval(['pp=co2(:,3);'])
  fwrite(fid,pp,'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save p,pp,t,q
  filemark= 8 * nlay;
  fwrite(fid,filemark,'integer*4');
  eval(['t=co2(:,4);'])
  fwrite(fid,t,'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save p,pp,t,q
  filemark= 8 * nlay;
  fwrite(fid,filemark,'integer*4');
  eval(['q=co2(:,5);'])
  fwrite(fid,q,'real*8');
  fwrite(fid,filemark,'integer*4');

  %Save the optical depth
  filemark= 8 * npts;
  for i = 1:nlay
    fwrite(fid,filemark,'integer*4');
    eval(['k' int2str(i) '=k00(i,:);'])
    fwrite(fid,eval(['k' int2str(i)]),'real*8');
    fwrite(fid,filemark,'integer*4');
    end

  fclose(fid);

  end 
