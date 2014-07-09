function [w,r] = readsolar(filename);

%% reads the double precision solar files used by kCARTA

fin = fopen(filename,'r','ieee-le');

flen = fread(fin, 1, 'integer*4');
fmin = fread(fin, 1, 'real*8');
fmax = fread(fin, 1, 'real*8');
df   = fread(fin, 1, 'real*8');
flen = fread(fin, 1, 'integer*4');

flen  = fread(fin, 1, 'integer*4');
r     = fread(fin, 10000, 'real*8');
flen  = fread(fin, 1, 'integer*4');

w = (1:10000)-1;
w = w*df + fmin;

fclose(fin);