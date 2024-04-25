function [w,d,iaProf,iaNumLay] = readsarta_jac(fname,iGID,rThreshold)

if nargin == 2
  %% comparing performance of /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/RODGERS/RODGERS21
  %%   when using finite diff jacs vs analystic jacs, the former is
  %%   better, since the jacs whose magnitude < 1e-4 are set to 0

  %{
  %% in MATLABCODE/CRODGERS_FAST_CLOUD after restart21A with ChrisSarta or SergioSarta, see compare_finitediff_analytic_jacresults_HALO_2019_04_25.m
  chris = load('chris_AIRS_HALO_2019_04_25_sarta_jac.mat');
  sergio = load('sergio_AIRS_HALO_2019_04_25_sarta_jac.mat');

  figure(100);
  iii = 10; semilogy(1:594,sergio.jacA(iii,:),'b.-',1:594,chris.jacA(iii,:),'rx-',1:594,jacA(iii,:),'g.-'); grid; boo=chris.jacA(iii,:); boo(boo < eps) = nan; nanmin(boo)
  ans = 1.5259e-04
  %}

  rThreshold = 1.0e-5;
  rThreshold = 1.0e-3;
  rThreshold = 4.0e-4;
end

%% see s_writefile.f90 line 2516
%%   readsarta_jac.m   : output is (numprof,numchan,numlay,'double')  [natural]
%%   readsarta_jacv2.m : output is (numlay,numchan,numprof,'single')  [CRODGERS_FAST_CLOUD/RODGERS/JAC_CODE/jacX_analytic.m]

%{
if iGID == 100
  disp('expecting ST +  temperature jacs')
elseif iGID == 200
  disp('expecting WGT FCNS')
elseif iGID == 300
  disp('expecting CLD JACS')
elseif length(intersect(iGID,[1 2 3 4 5 6 9 12])) == 1
  disp('expecting gas jacs')
else
  error('iGID = [1 2 3 4 5 6 9 12] for for [WV,OZ,...HNO3]     or  [100,200,300] for  [TZ,WGTFCN,CLD] jacs')
end
%}

if length(intersect(iGID,[100 200 300 1 2 3 4 5 6 9 12])) == 0
  error('iGID = [1 2 3 4 5 6 9 12] for for [WV,OZ,...HNO3]     or  [100,200,300] for  [TZ,WGTFCN,CLD] jacs')
end

if ~exist(fname)
  str = ['file "' fname '" DNE'];
  error(str)
end

[fin,msg] = fopen(fname,'r','ieee-be');  %% SARTA COMPILED with byte swap to IEEE_BE
if fin == -1
  error(['error opening input file ', fin,msg]);
end

flen    = fread(fin, 1, 'integer*4');
numprof = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');

flen    = fread(fin, 1, 'integer*4');
numchan = fread(fin, 1, 'integer*4');
flen    = fread(fin, 1, 'integer*4');

%fprintf(1,'expecting %5i profiles with %4i channels \n',numprof,numchan)

w = [];
if iGID == 100
  d = zeros(numprof,numchan,101);
elseif iGID == 300
  d = zeros(numprof,numchan,7);  %% cfrac1,amt1,sze1,cfrac2,amt2,sze2,cfrac12
  d = zeros(numprof,numchan,11); %% cfrac1,amt1,sze1,top1,bot2,cfrac2,amt2,sze2,top2,bot2,cfrac12
  d = zeros(numprof,numchan,12); %% cfrac1,amt1,sze1,top1,bot2,cfrac2,amt2,sze2,top2,bot2,cfrac12,stemp
else
  d = zeros(numprof,numchan,100);
end

flen = fread(fin, 1, 'integer*4');
w = fread(fin,numchan,'real*4');
flen = fread(fin, 1, 'integer*4');

%disp('looking at basic memory usage in readsarta_jac.m, d=zeros(numprof,numchan,X) where X=100 for profile or eg 12 for clds')
%whos d w
%monitor_memory_whos
%disp('looking at basic memory usage in readsarta_jac.m')

for iP = 1 : numprof
  flen = fread(fin, 1, 'integer*4');
  iaProf(iP)   = fread(fin, 1, 'integer*4');
  iaNumLay(iP) = fread(fin, 1, 'integer*4');
  iNC          = fread(fin, 1, 'integer*4');
  iWhich       = fread(fin, 1, 'integer*4');
  flen = fread(fin, 1, 'integer*4');
  if iNC ~= numchan
    [iNC numchan]
    error('iNC (inside file) and numchan (head of file) are different!!!')
  end
  if iWhich ~= iGID
    [iWhich iGID]
    error('iWhich (in file)  and iGID (user supplied) are different!!!')
  end

  for iL = 1 : iaNumLay(iP)
    flen = fread(fin, 1, 'integer*4');
    junk = fread(fin,numchan,'real*4');
    d(iP,:,iL) = junk;
    flen = fread(fin, 1, 'integer*4');
  end

  %if mod(iP,1000) == 0
  %  fprintf(1,'profile iP = %6i loop iL = %3i of %3i readsarta_jac.m, mem usuage  = %4i \n',iP,iL,iaNumLay(iP)); 
  %  whos d junk
  %  monitor_memory_whos
  %end

end

fclose(fin);

d(abs(d) < rThreshold) = 0.0;
