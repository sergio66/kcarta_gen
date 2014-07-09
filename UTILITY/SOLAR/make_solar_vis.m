% this is from /salsify/users/sergio/KCARTA/MATLAB/GENERAL
% writes solar spectra to many files so that kCARTAv1.04+ can read it in

load_solar_vis_data

nbox = 5; 
pointsPerChunk = 10000; 

wnVIS1 = 3000; 
wnVIS2 = 6050; 

wnVIS1 = 6050; 
wnVIS2 = 22000; 

fmin = wnVIS1;  
topts = runXtopts_params_smart(fmin);  
dv = topts.ffin*nbox*pointsPerChunk; 

ind = 1 : 10000;

wx = fliplr(wx);
rx = fliplr(rx);
while fmin <= wnVIS2 

  topts = runXtopts_params_smart(fmin); 
  dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = fmin + dv; 
  wvis0 = []; 

  fout = (ind-1)*topts.ffin*nbox + fmin;
  rout = interp1(wx,rx,fout)/6.785087652174316e-5/1000; 
      %%kCARTA takes care of both of these
  fname=['rad_solar' num2str(fmin) '.dat'];

  %%%%%%%%%%%%%%%%%%%  check out ieee-le or ieee-be %%%%%%%%%%%%%%%%%%%%%
  fid=fopen(fname,'w','ieee-le');
  %%%%%%%%%%%%%%%%%%%  check out ieee-le or ieee-be %%%%%%%%%%%%%%%%%%%%%

  fs=fout(ind(1));
  fe=fout(ind(length(ind)));

  df=abs(fout(length(fout))-fout(1));
  df=df/(length(fout)-1);

  fprintf(1,'%s \n',fname);
  fprintf(1,'  fs,fe,length(ind) = %10.5f %10.5f %6i \n \n',fs,fe,length(ind)')

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

  fmin = fmin + dv;
  end

