%writes solar spectra to a file so that kCARTAv1.04+ can read it in

load /salsify/users/hannon/Solar_data/Chris_Barnet/fine_solar.mat

[mm,nn]=size(fout);
if (mm > nn)
  npts=mm;
else
  npts=nn;
  end

fid=fopen('rad_solar.dat','w','ieee-be');

ind=find((fout >= 605.0) & (fout <= 2830-0.0025));

fs=fout(ind(1))
fe=fout(ind(length(ind)))
length(ind)

df=abs(fout(length(fout))-fout(1));
df=df/(length(fout)-1)

% header = fstart,fstep
filemark= 8 + 8 + 8;
fwrite(fid,filemark,'integer*4');
fwrite(fid,[fs,fe,df],'real*8');
fwrite(fid,filemark,'integer*4');

%Save the solar radout array
filemark= 8 * length(ind);
fwrite(fid,filemark,'integer*4');
fwrite(fid,rout(ind),'real*8');
fwrite(fid,filemark,'integer*4');

fclose(fid);
