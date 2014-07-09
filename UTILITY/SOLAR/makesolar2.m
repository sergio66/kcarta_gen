% this is from /salsify/users/sergio/KCARTA/MATLAB/GENERAL
% writes solar spectra to many files so that kCARTAv1.04+ can read it in

load /salsify/users/hannon/Solar_data/Chris_Barnet/fine_solar.mat

[mm,nn]=size(fout);
if (mm > nn)
  npts=mm;
else
  npts=nn;
  end

blah=(2805-605)/25 + 1;

number=605;
ind=find((fout >= number) & (fout <= number+25-0.0025));

for ii=1:blah
  number=605+(ii-1)*25;
  fname=['rad_solar' num2str(number) '.dat'];

  %%%%%%%%%%%%%%%%%%%  check out ieee-le or ieee-be %%%%%%%%%%%%%%%%%%%%%
  fid=fopen(fname,'w','ieee-le');
  %%%%%%%%%%%%%%%%%%%  check out ieee-le or ieee-be %%%%%%%%%%%%%%%%%%%%%

  fs=fout(ind(1));
  fe=fout(ind(length(ind)));

  df=abs(fout(length(fout))-fout(1));
  df=df/(length(fout)-1);

  fname
  fprintf(1,'fs,fe,length(ind) = %10.5f %10.5f %6i \n \n',fs,fe,length(ind)')

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

  ind=ind+10000;

  end
