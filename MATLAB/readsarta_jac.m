function [w,d,iaProf,iaNumLay] = readsarta_jac(fname,iGID)

%% see s_writefile.f90 line 2516

if iGID == 100
  disp('expecting ST +  temperature jacs')
elseif iGID == 200
  disp('expecting WGT FCNS')
elseif length(intersect(iGID,[1 2 3 4 5 6 9 12])) == 1
  disp('expecting gas jacs')
else
  error('iGID = [1 2 3 4 5 6 9 12]  [100,200] for WV,OZ,WGTFCN,TZ jacs')
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

fprintf(1,'expecting %5i profiles with %4i channels \n',numprof,numchan)

w = [];
if iGID == 100
  d = zeros(numprof,numchan,101);
else
  d = zeros(numprof,numchan,100);
end

flen = fread(fin, 1, 'integer*4');
w = fread(fin,numchan,'real*4');
flen = fread(fin, 1, 'integer*4');

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
end

fclose(fin);
