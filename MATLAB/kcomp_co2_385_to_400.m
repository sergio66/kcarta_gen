%% this is in /asl/data/kcarta/H2012.ieee-le/IR605/lblrtm12.4/etc.ieee-le/CO2_400ppmv
%% and in /home/sergio/HITRAN2UMBCLBL/FORTRAN/mat2for
%% and in /home/sergio/KCARTA/MATLAB

addpath /home/sergio/HITRAN2UMBCLBL/FORTRAN/for2mat
addpath /home/sergio/HITRAN2UMBCLBL/FORTRAN/mat2for

inputCO2ppm = 385;
outputCO2ppm = 400;

thedir = dir(['../CO2_385ppmv/*_g2.dat']);
for ii = 1 : length(thedir)
  filename = ['../CO2_385ppmv/' thedir(ii).name];
  fout = [thedir(ii).name];  
  [gid0,fr0,kcomp0,B0] = for2mat_kcomp_reader(filename);
  kcompX = kcomp0 * (400/385)^(1/4);
  kcompX = kcomp0 * (outputCO2ppm/inputCO2ppm)^(1/4);
  dtype = 'ieee-le';
  for2for_gasmult(gid0,fr0,kcompX,B0,1.000,fout,dtype);
  [gid1,fr1,kcomp1,B1] = for2mat_kcomp_reader(fout);

  od0 = (B0*squeeze(kcomp0(:,:,6))); od0 = od0 .^ 4;
  od1 = (B1*squeeze(kcomp1(:,:,6))); od1 = od1 .^ 4;

  fprintf(1,'gid = %2i vchunk = %4i change OD using ratio %3i/%3i = %8.6f \n',gid0,fr0(1),outputCO2ppm,inputCO2ppm);
  figure(1); clf; semilogy(fr0,sum(od0'),'b',fr1,sum(od1'),'r'); title(num2str(fr0(1)))
  figure(2); clf; plot(fr0,sum(od1') ./ sum(od0'),'r'); title(num2str(fr0(1)))
  pause(0.1)
  
end
