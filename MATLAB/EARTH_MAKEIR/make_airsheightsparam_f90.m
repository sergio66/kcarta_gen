addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

[h,ha,p,pa] = rtpread('klayersV205_Data/Data/adafgl_16Aug2010_op.rtp');

fid = fopen('airsheights.param_earth','w');
fprintf(fid,'! 100 layer heights for the DEFINITION of the kCARTA database \n');
fprintf(fid,'  \n');
fprintf(fid,'      INTEGER iXPlanet3 \n');
fprintf(fid,'      DATA iXPlanet3 /03/ \n');
fprintf(fid,'  \n');
fprintf(fid,'       REAL DatabaseHeight(kMaxLayer) \n');
fprintf(fid,'       INTEGER IPHEIGHT \n');
fprintf(fid,'  \n');
fprintf(fid,'! note that the program expects heights in km, so these numbers will have \n');
fprintf(fid,'! to be read in and divided by 1000.0 \n');
fprintf(fid,'  \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      from /home/sergio/KCARTA/MATLAB/EARTH_MAKEIR/make_airsheightsparam_f90.m \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = make_rtp_plays(p);
playsHGT = interp1(log(p.plevs(:,6)),p.palts(:,6),log(p.plays(1:100,6)),[],'extrap');
playsHGT = interp1(log(p.plevs(:,6)),p.palts(:,6),log(p.plays(1:100,6)),[],'extrap');
playsHGT = sort(playsHGT,'ascend');

fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      Load up p.palts, get height at avg pressure from the AIRS layer heights US Standard Profile (in m) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'       DATA  (DatabaseHEIGHT(IPHEIGHT), IPHEIGHT = 1,100,1 ) & \n');

ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = playsHGT(ind);
  fprintf(fid,'       /%13.7e, %13.7e, %13.7e, %13.7e, & \n',data);

for ii = 2 : 24
  ind = (1:4)+(ii-1)*4;
  data = playsHGT(ind);
  fprintf(fid,'        %13.7e, %13.7e, %13.7e, %13.7e, & \n',data);
end
ii = 25;
  ind = (1:4)+(ii-1)*4;
  data = playsHGT(ind);
  fprintf(fid,'        %13.7e, %13.7e, %13.7e, %13.7e / \n',data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);

