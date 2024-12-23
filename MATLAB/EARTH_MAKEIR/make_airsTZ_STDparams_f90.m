addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

[h,ha,p,pa] = rtpread('klayersV205_Data/Data/adafgl_16Aug2010_op.rtp');

fid = fopen('airsTZ_STD.param_earth','w');
fprintf(fid,'! 100 layer heights and Temperatures and GasAmounts for the DEFINITION of the kCARTA database \n');
fprintf(fid,'  \n');
fprintf(fid,'      INTEGER iXPlanet8 \n');
fprintf(fid,'      DATA iXPlanet8 /08/ \n');
fprintf(fid,'  \n');
fprintf(fid,'       REAL DatabaseTZ(kMaxLayer) ! Kelvin\n');
fprintf(fid,'       REAL DatabaseQZ(kMaxLayer) ! molecules/cm2\n');
fprintf(fid,'       INTEGER IPLAY \n');
fprintf(fid,'  \n');
fprintf(fid,'! note that the program expects T(z) in K and gas amounts in moles/cm2 \n');
fprintf(fid,'  \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      from /home/sergio/KCARTA/MATLAB/EARTH_MAKEIR/make_airsTZ_STDparams_f90.m \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 8.31; %% J/mol/K

p = make_rtp_plays(p);
playsHGT = p.plays(1:100,6);
[playsHGT,I] = sort(playsHGT,'ascend');

TlaysHGT = p.ptemp(1:100,6);
TlaysHGT = TlaysHGT(I);

dZlaysHGT = abs(p.palts(1:100,6) - p.palts(2:101,6));
dZlaysHGT = dZlaysHGT(I);

QlaysHGT = playsHGT*100 .* dZlaysHGT ./R ./TlaysHGT; %% moles/m2 * m = moles/m2
MRlaysHGT = QlaysHGT/1e4 * 6.023e23;                 %% molecules/cm2
QlaysHGT = QlaysHGT/1e4;                             %% moles/cm2
QlaysHGT = QlaysHGT/1e3;                             %% kilomoles/cm2


fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      Use p.ptemp,p.salti t get T(z) Q(z) at avg pressure from the AIRS layer heights US Standard Profile (in K and kmol/cm2) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'!      This is TZ in K \n')
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'       DATA  (DatabaseTZ(IPLAY), IPLAY = 1,100,1 ) & \n');
ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = TlaysHGT(ind);
  fprintf(fid,'       /%13.7f, %13.7f, %13.7f, %13.7f, & \n',data);

for ii = 2 : 24
  ind = (1:4)+(ii-1)*4;
  data = TlaysHGT(ind);
  fprintf(fid,'        %13.7f, %13.7f, %13.7f, %13.7f, & \n',data);
end
ii = 25;
  ind = (1:4)+(ii-1)*4;
  data = TlaysHGT(ind);
  fprintf(fid,'        %13.7f, %13.7f, %13.7f, %13.7f / \n',data);

fprintf(fid,' \n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
fprintf(fid,' \n');

fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'!      This is QZ (all gases) in kmol/cm2 \n')
fprintf(fid,'!      -----------------------------------------------------------------\n');

fprintf(fid,'       DATA  (DatabaseQZ(IPLAY), IPLAY = 1,100,1 ) & \n');
ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = QlaysHGT(ind);
  fprintf(fid,'       /%13.7e, %13.7e, %13.7e, %13.7e, & \n',data);

for ii = 2 : 24
  ind = (1:4)+(ii-1)*4;
  data = QlaysHGT(ind);
  fprintf(fid,'        %13.7e, %13.7e, %13.7e, %13.7e, & \n',data);
end
ii = 25;
  ind = (1:4)+(ii-1)*4;
  data = QlaysHGT(ind);
  fprintf(fid,'        %13.7e, %13.7e, %13.7e, %13.7e / \n',data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now check CO2 ppmv
printarray(p.gas_1(1:100,6)./MRlaysHGT * 1e6,'WV PPMV')
printarray(p.gas_2(1:100,6)./MRlaysHGT * 1e6,'CO2 PPMV')
printarray(p.gas_4(1:100,6)./MRlaysHGT * 1e6,'N2O PPMV')
printarray(p.gas_6(1:100,6)./MRlaysHGT * 1e6,'CH4 PPMV')

disp('QlaysHGT./MRlaysHGT = 0.166666e-26 = 1/(Na*1000) = 0.166666e-26 = 1.66666e-27')
