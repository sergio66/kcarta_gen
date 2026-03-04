%{
%% from SRCv1.22_F90
make -f makefile 385_H12_CO2_UMBC_default_f90_short
make -f makefile 400_H16_default_f90_short
make -f makefile 400_H20_default_f90_short

%% from WORK
rm junk12.dat junk16.dat junk20.dat
time ../../BIN/kcarta.x90_v1.22_385ppmv_H12_CO2_UMBC quickuseF90.nml junk12.dat
time ../../BIN/kcarta.x90_v1.22_400ppmv_H16 quickuseF90.nml          junk16.dat
time ../../BIN/kcarta.x90_v1.22_400ppmv_H20 quickuseF90.nml          junk20.dat
time ../../BIN/kcarta.x90_v1.22_400ppmv_H24 quickuseF90.nml          junk24.dat
%}

addpath ../../MATLAB
[d12,w] = readkcstd('junk12.dat'); t12 = rad2bt(w,d12);
[d16,w] = readkcstd('junk16.dat'); t16 = rad2bt(w,d16);
[d20,w] = readkcstd('junk20.dat'); t20 = rad2bt(w,d20);
[d24,w] = readkcstd('junk24.dat'); t24 = rad2bt(w,d24);

figure(1);
plot(w,[t12 t16 t20 t24]); legend('H12','H16','H20','H24');
plot(w,t24-t12,'b',w,t24-t16,'g',w,t24-t20,'r'); legend('H24-H12','H24-H16','H24-H20','location','best'); title('monochromatic kcarta')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/git/matlabcode
[fc,qc] = quickconvolve(w,[d12 d16 d20 d24],1,1); tc = rad2bt(fc,qc); tc12 = tc(:,1); tc16 = tc(:,2); tc20 = tc(:,3); tc24 = tc(:,4);

figure(2); 
plot(fc,[tc12 tc16 tc20 tc24]); legend('H12','H16','H20','H24');
plot(fc,tc24-tc12,'b',fc,tc24-tc16,'g',fc,tc24-tc20,'r'); legend('H24-H12','H24-H16','H24-H20','location','best'); title('gaussian FWHM = 1 cm-1')
