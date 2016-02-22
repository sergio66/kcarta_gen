sedder = [sedder0 ' template_twoslabLOOP.nml  >  ' outnmlx];
eval(sedder)

fout = ['JUNK/kctwoslab_chunk_' num2str(f1) '_prof_' num2str(iRTP) '.dat'];
fout = mktemp(fout);       %% oops when we create this temp name, we make an empty file
rmer = ['!/bin/rm ' fout]; eval(rmer); %% rm the empty file before trying to run kCARTA

ughN = ['JUNK/ugh' num2str(iRTP)];
kcer = ['!' kcX ' ' outnmlx ' ' fout ' >& ' ughN];
fprintf(1,'running kcarta : %s \n',kcer);
eval(kcer);

[r2s,w] = readkcstd(fout);
fluxout = [fout '_OLR'];
[flux2s,w] = readkcflux(fluxout);   %% chose kFLux .EQ. 4

rmer = ['!/bin/rm ' fout ' ' fluxout ' ' tmp_rtp ' ' ughN];
eval(rmer);
