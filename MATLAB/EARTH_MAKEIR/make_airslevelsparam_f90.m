addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

[h,ha,p,pa] = rtpread('klayersV205_Data/Data/adafgl_16Aug2010_op.rtp');

fid = fopen('airslevels.param_earth','w');
fprintf(fid,'! 101 pressure levels for the DEFINITION of the kCARTA database \n');
fprintf(fid,'  \n');
fprintf(fid,'      INTEGER iXPlanet2 \n');
fprintf(fid,'      DATA iXPlanet2 /03/ \n');
fprintf(fid,'  \n');
fprintf(fid,'       REAL DATABASELEV(kMaxLayer+1) \n');
fprintf(fid,'       INTEGER IPLEV \n');
fprintf(fid,'  \n');
fprintf(fid,'!       REAL delta  \n');
fprintf(fid,'!       PARAMETER(delta=1.0e-9)    !just an "epsilon" needed when finding \n');
fprintf(fid,'!                                  !pressure layer bndries \n');
fprintf(fid,'  \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      from /home/sergio/KCARTA/MATLAB/EARTH_MAKEIR/make_airslevelsparam_f90.m \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plevs = p.plevs(:,6);
plevs = sort(plevs,'ascend');

fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      Load up p.plev, write them out  US Standard Profile (in mb) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');

fprintf(fid,'       DATA  (DatabaseLEV(IPLEV), IPLEV = 101,52,-1 ) & \n');

ii = 1;
  ind = (1:5)+(ii-1)*5;
  data = plevs(ind);
  fprintf(fid,'       /%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, & \n',data);

for ii = 2 : 10
  ind = (1:5)+(ii-1)*5;
  data = plevs(ind);
  if ii < 10
    fprintf(fid,'        %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, & \n',data);
  else
    fprintf(fid,'        %12.4f, %12.4f, %12.4f, %12.4f, %12.4f / \n',data);
  end
end

fprintf(fid,'       DATA  (DatabaseLEV(IPLEV), IPLEV = 51,1,-1 ) & \n');
ii = 11;
  ind = (1:5)+(ii-1)*5;
  data = plevs(ind);
  fprintf(fid,'       /%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, & \n',data);

for ii = 12 : 20
  ind = (1:5)+(ii-1)*5;
  data = plevs(ind);
  fprintf(fid,'        %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, & \n',data);
end

ii = 21;
  ind = (1:5)+(ii-1)*5;
   ind = ind(ind <= 101);
  data = plevs(ind);
  fprintf(fid,'        %12.4f / \n',data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);

