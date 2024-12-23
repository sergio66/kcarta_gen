addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

[h,ha,p,pa] = rtpread('klayersV205_Data/Data/adafgl_16Aug2010_op.rtp');

fid = fopen('airslevelheights.param_earth','w');
fprintf(fid,'! from KCARTADATA/RefProf.For.v107up/mypref.op.new \n');
fprintf(fid,'! the heights of the pressure levels for the kCARTA DATABASE \n');

fprintf(fid,'! 101 level heights for the DEFINITION of the kCARTA database \n');
fprintf(fid,'  \n');
fprintf(fid,'      INTEGER iXPlanet3 \n');
fprintf(fid,'      DATA iXPlanet3 /03/ \n');
fprintf(fid,'  \n');
fprintf(fid,'       REAL DatabaseHeight(kMaxLayer+1) \n');
fprintf(fid,'       INTEGER IPLEV \n');
fprintf(fid,'  \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      from /home/sergio/KCARTA/MATLAB/EARTH_MAKEIR/make_airslevelheightsparam_f90.m \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plevsHGT = p.palts(:,6)/1000;
plevsHGT = sort(plevsHGT,'descend');

fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      Load up p.palts, write them out  US Standard Profile (in km) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'       DATA  (DatabaseLEVHEIGHTS(IPLEV), IPLEV = 101,1,-1 ) & \n');

ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = plevsHGT(ind);
  fprintf(fid,'       /%12.5f, %12.5f, %12.5f, %12.5f, & \n',data);

for ii = 2 : 25
  ind = (1:4)+(ii-1)*4;
  data = plevsHGT(ind);
  fprintf(fid,'        %12.5f, %12.5f, %12.5f, %12.5f, & \n',data);
end
ii = 26;
  ind = (1:4)+(ii-1)*4;
   ind = ind(ind <= 101);
  data = plevsHGT(ind);
  fprintf(fid,'        %12.5f / \n',data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);

