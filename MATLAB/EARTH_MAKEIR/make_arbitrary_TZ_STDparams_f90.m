addpath /home/sergio/git/matlabcode/matlibSergio/matlib2025/h4tools/
addpath /home/sergio/git/matlabcode

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% lrwxrwxrwx 1 sergio pi_sergio   32 Dec 19  2025 arbitrary_TZ_STDparam_earth -> arbitrary_TZ_STD.param_earth_PBL

% See ~/git/kcarta_gen/SRCv1.22_f90/Revisions_arbitrary_plevs.txt
%   07/25/26        ~/git/sarta_scatter_rtp_klayers_sergio/KLAYERS_CHepplew_PBL/klayersV205/Grid/ is linked to
%                   ~/git/kcarta_gen/MATLAB/EARTH_MAKEIR/make_arbitrary_TZ_STDparams_f90.m

iWhichNewType = 1;   %% USSTD for 100 L2 layers --> OCO2 PBL layers (125 m thick at surface)
iWhichNewType = 50;  %% USSTD for 100 L2 layers --> aircraft at 1 mb   = 50 km
iWhichNewType = 20;  %% USSTD for 100 L2 layers --> aircraft at 57 mb  = 20 km
iWhichNewType = 12;  %% USSTD for 100 L2 layers --> aircraft at 250 mb = 12 km

dir0 = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/';

if iWhichNewType == 1
  %% <<<<<<<<<< have to define the OCO2 PBL file here, and define the output name here >>>>>>>>>>>>>>>>>>
  fin = '/home/chepplew/projects/klayers_wrk/regr49_pbl.op.rtp';    iSTD = 49;
  fin = [dir0 /regr49_pbl_probably370ppm.op.rtp'];                  iSTD = 49;  
  fout = 'arbitrary_TZ_STD.param_earth_PBL';
  comment = '! ARBITRARY Irion/Hepplewhite PBL layer heights and Temperatures and GasAmounts';
  %% <<<<<<<<<< have to define the OCO2 PBL file here, and define the output name here >>>>>>>>>>>>>>>>>>
elseif iWhichNewType == 12
  %% <<<<<<<<<< have to define the 12 km (250 mb) aircraft  file here, and define the output name here >>>>>>>>>>>>>>>>>>
  iSTD = 49;    
  fin = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1100_with_co2_400ppm_9gases_unitemiss_aircraft_12km.op.rtp';
  fout = 'aircraft_12km_250mb_TZ_STD.param_earth';
  comment = '! ARBITRARY pend = 250 mb (12 km) layer heights and Temperatures and GasAmounts';
  %% <<<<<<<<<< have to define the 12 km (250 mb) aircraft  file here, and define the output name here >>>>>>>>>>>>>>>>>>
else
  error('unknow iWhichNewType')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fin)
  error('cannot find your definition rtp file here');
end
if exist(fout)
  fprintf(1,'fout %s exists \n',fout)
  error('fout exists')
end

[h0,ha0,p0,pa0] = rtpread('klayersV205_Data/Data/adafgl_16Aug2010_op.rtp');
[h, ha, p, pa ] = rtpread(fin);
fid = fopen(fout,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'%s \n', comment);
fprintf(fid,'  \n');
fprintf(fid,'      INTEGER iXPlanet9 \n');
fprintf(fid,'      DATA iXPlanet8 /09/ \n');
fprintf(fid,'  \n');
fprintf(fid,'       REAL ARBDatabaseTZ(kMaxLayer)           ! Kelvin \n');
fprintf(fid,'       REAL ARBDatabaseQZ(kMaxLayer)           ! molecules/cm2 \n');
fprintf(fid,'       REAL ARBDATABASELEVHEIGHTS(kMaxLayer+1) ! km \n');
fprintf(fid,'       REAL ARBDATABASEPLEVS(kMaxLayer+1)      ! mb \n');
%%%%% fprintf(fid,'      INTEGER IPLAY \n');
fprintf(fid,'  \n');
fprintf(fid,'! note that the program expects T(z) in K and gas amounts in moles/cm2 \n');
fprintf(fid,'  \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      from /home/sergio/KCARTA/MATLAB/EARTH_MAKEIR/make_arbitrary_TZ_STDparams_f90.m \n');
fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'  \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = 8.31; %% J/mol/K

p = make_rtp_plays(p);

plevs    = p.plevs(1:101,iSTD);

plevsHGT = p.palts(1:101,iSTD)/1000;

playsHGT = p.plays(1:100,iSTD);
[playsHGT,I] = sort(playsHGT,'ascend');

TlaysHGT = p.ptemp(1:100,iSTD);
TlaysHGT = TlaysHGT(I);

dZlaysHGT = abs(p.palts(1:100,iSTD) - p.palts(2:101,iSTD));
dZlaysHGT = dZlaysHGT(I);

QlaysHGT = playsHGT*100 .* dZlaysHGT ./R ./TlaysHGT; %% moles/m2 * m = moles/m2
MRlaysHGT = QlaysHGT/1e4 * 6.023e23;                 %% molecules/cm2
QlaysHGT = QlaysHGT/1e4;                             %% moles/cm2
QlaysHGT = QlaysHGT/1e3;                             %% kilomoles/cm2


fprintf(fid,'!      ----------------------------------------------------------------- \n');
fprintf(fid,'!      Use p.ptemp,p.salti t get T(z) Q(z) at avg pressure from the AIRS layer heights US Standard Profile (in K and kmol/cm2) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'!      This is TZ in K \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'       DATA  (ARBDatabaseTZ(IPLAY), IPLAY = 1,100,1 ) & \n');
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
fprintf(fid,'!      This is QZ (all gases) in kmol/cm2 \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');

fprintf(fid,'       DATA  (ARBDatabaseQZ(IPLAY), IPLAY = 1,100,1 ) & \n');
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

fprintf(fid,' \n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
fprintf(fid,' \n');

fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'!      This is PLEVS with the ARBITRARY layer boundary pressure heights (in mb) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');

fprintf(fid,'       DATA  (ARBDATABASEPLEVS(IPLAY), IPLAY = 1,101,1 ) & \n');
ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = plevs(ind);
  fprintf(fid,'       /%13.7f, %13.7f, %13.7f, %13.7f, & \n',data);

for ii = 2 : 25
  ind = (1:4)+(ii-1)*4;
  data = plevs(ind);
  fprintf(fid,'        %13.7f, %13.7f, %13.7f, %13.7f, & \n',data);
end
ii = 26;
  ind = (1:4)+(ii-1)*4;
  ind = ind(1);
  data = plevs(ind);
  fprintf(fid,'        %13.7f / \n',data);

fprintf(fid,' \n');
fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n');
fprintf(fid,' \n');

fprintf(fid,'!      -----------------------------------------------------------------\n');
fprintf(fid,'!      This is PALTS with the ARBITRARY layer boundary pressure heights (in km) \n');
fprintf(fid,'!      -----------------------------------------------------------------\n');

fprintf(fid,'       DATA  (ARBDATABASELEVHEIGHTS(IPLAY), IPLAY = 1,101,1 ) & \n');
ii = 1;
  ind = (1:4)+(ii-1)*4;
  data = plevsHGT(ind);
  fprintf(fid,'       /%13.7f, %13.7f, %13.7f, %13.7f, & \n',data);

for ii = 2 : 25
  ind = (1:4)+(ii-1)*4;
  data = plevsHGT(ind);
  fprintf(fid,'        %13.7f, %13.7f, %13.7f, %13.7f, & \n',data);
end
ii = 26;
  ind = (1:4)+(ii-1)*4;
  ind = ind(1);
  data = plevsHGT(ind);
  fprintf(fid,'        %13.7f / \n',data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now check CO2 ppmv
printarray(p.gas_1(1:100,iSTD)./MRlaysHGT * 1e6,'WV PPMV')
printarray(p.gas_2(1:100,iSTD)./MRlaysHGT * 1e6,'CO2 PPMV')
printarray(p.gas_4(1:100,iSTD)./MRlaysHGT * 1e6,'N2O PPMV')
printarray(p.gas_6(1:100,iSTD)./MRlaysHGT * 1e6,'CH4 PPMV')

disp('QlaysHGT./MRlaysHGT = 0.166666e-26 = 1/(Na*1000) = 0.166666e-26 = 1.66666e-27')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,' >>>> now symbolically linking    arbitrary_TZ_STDparam_earth    to    %s \n',fout);

rmer = ['!/bin/rm arbitrary_TZ_STDparam_earth'];        eval(rmer);
lner = ['!ln -s ' fout ' arbitrary_TZ_STDparam_earth']; eval(lner)   

disp(' ')
disp('now see Makefile_v122_Intelf90 --> Makefile_tar_objs_data_f90 --> Makefile_tar_objs_data_f90_datafix');
disp('that has cd ../INCLUDE; cd EARTH_database_params; ./lner_EARTH_database_params.sc')
disp('that latter has you going back to SRCv1.22_f90 to run "cp_param_files_to_f90_v122.sc" ')
disp('here we run cp_param_files_to_f90_v122_arb.sc')
disp(' ')
disp('now re-Make the kcarta code')

runner = ['!cp_param_files_to_f90_v122_arb.sc']; eval(runner);
