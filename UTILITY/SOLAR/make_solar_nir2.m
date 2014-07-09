% this is from /salsify/users/sergio/KCARTA/MATLAB/GENERAL
% writes solar spectra to many files so that kCARTAv1.04+ can read it in

load_solar_vis_data
%%% or can load /home/sergio/KCARTADATA/General/SOLAR_ATMOS/sergio_solar.mat

nbox = 5; 
pointsPerChunk = 10000; 

wnVIS1 = 2830; 
wnVIS2 = 3330; 

fmin = wnVIS1;

fuse = 2000;  
topts = runXtopts_params_smart(fuse);  
dv = topts.ffin*nbox*pointsPerChunk; 

ind = 1 : 10000;

wx = fliplr(wx);
rx = fliplr(rx);
while fmin <= wnVIS2 

  %topts = runXtopts_params_smart(fmin); 
  %dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = fmin + dv; 

  wvis0 = []; 

  fout = (ind-1)*topts.ffin*nbox + fmin;
  rout = interp1(wx,rx,fout)/6.785087652174316e-5/1000; 

  rjunk = ttorad(fout,6100);  %% KCARTA worries about the 6.785087652174316e-5 
  if fmin >= 4250 & fmax <= 4450
    iLoadJack = input('load Jack Kumers data instead of SBDART? ');
    iX = input(' all lines included (+1) or lines removed (-1)? ');
    if iLoadJack > 0
      if iX > 0
        dada = load('/home/sergio/KUMER/lockheed_ascii_4250_4450_solarA.txt');
        rout1 = interp1(dada(:,1),dada(:,2),fout,[],'extrap');
        rout1 = rout1/6.785087652174316e-5;
      else
        dada = load('/home/sergio/KUMER/lockheed_ascii_4250_4450_solarB.txt');
        rout1 = interp1(dada(:,1),dada(:,3),fout,[],'extrap');
        rout1 = rout1/6.785087652174316e-5;
        end
      plot(fout,rout,fout,rjunk/1000,fout,rout1); 
      title('SBDART  ttorad(w,6100) kumer');
      pause
      rout = rout1;
      end
    end

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

