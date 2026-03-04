%{
%% from SRCv1.22_F90
make -f makefile 385_H12_CO2_UMBC_default_f90_short
make -f makefile 400_H16_default_f90_short
make -f makefile 400_H20_default_f90_short

%% from WORK
rm junk12.dat junk16.dat junk20.dat junk12.dat_CLD junk16.dat_CLD junk20.dat_CLD
time ../../BIN/kcarta.x90_v1.22_385ppmv_H12_CO2_UMBC quickuse_cldF90.nml junk12.dat
time ../../BIN/kcarta.x90_v1.22_400ppmv_H16 quickuse_cldF90.nml          junk16.dat
time ../../BIN/kcarta.x90_v1.22_400ppmv_H20 quickuse_cldF90.nml          junk20.dat
%}

addpath ../../MATLAB
[d12,w] = readkcstd('junk12.dat'); t12 = rad2bt(w,d12);
[d16,w] = readkcstd('junk16.dat'); t16 = rad2bt(w,d16);
[d20,w] = readkcstd('junk20.dat'); t20 = rad2bt(w,d20);

addpath /umbc/rs/pi_sergio/WorkDirDec2025/matlabcode
[dx12,w] = readkcstd_twoslabclds('junk12.dat','junk12.dat_CLD'); tx12 = rad2bt(w,dx12);
[dx16,w] = readkcstd_twoslabclds('junk16.dat','junk16.dat_CLD'); tx16 = rad2bt(w,dx16);
[dx20,w] = readkcstd_twoslabclds('junk20.dat','junk20.dat_CLD'); tx20 = rad2bt(w,dx20);

[sum(d12(:,5)-dx12) sum(d16(:,5)-dx16) sum(d20(:,5)-dx20)]
[sum(abs(d12(:,5)-dx12)) sum(abs(d16(:,5)-dx16)) sum(abs(d20(:,5)-dx20))]

i1231 = find(w >= 1231,1);
figure(1); plot(1:6,[t12(i1231,:) tx12(i1231)],'o-')

figure(1);
plot(w,[t12(:,5) t16(:,5) t20(:,5)]); legend('H12','H16','H20');
plot(w,t20(:,5),'r',w,t20(:,4),'b'); legend('H20cld','H20clr');

plot(w,[t12(:,5) t16(:,5) t20(:,5)],w,[t12(:,4) t16(:,4) t20(:,4)],'--'); legend('H12cld','H16cld','H20cld','H12clr','H16clr','H20clr');
plot(w,t20(:,5)-t12(:,5),w,t20(:,5)-t16(:,5)); legend('H20-H12','H20-H16','location','best'); title('monochromatic kcarta')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/git/matlabcode
[fc,qc_clr] = quickconvolve(w,[d12(:,4) d16(:,4) d20(:,4)],1,1); tc_clr = rad2bt(fc,qc_clr); tc12_clr = tc_clr(:,1); tc16_clr = tc_clr(:,2); tc20_clr = tc_clr(:,3);
[fc,qc_cld] = quickconvolve(w,[d12(:,5) d16(:,5) d20(:,5)],1,1); tc_cld = rad2bt(fc,qc_cld); tc12_cld = tc_cld(:,1); tc16_cld = tc_cld(:,2); tc20_cld = tc_cld(:,3);

figure(2); 
plot(fc,[tc12_clr tc16_clr tc20_clr]); legend('H12','H16','H20');
plot(fc,tc20_clr-tc12_clr,fc,tc20_clr-tc16_clr); legend('H20-H12','H20-H16','location','best'); title('gaussian FWHM = 1 cm-1 CLR')

figure(2); 
plot(fc,[tc12_cld tc16_cld tc20_cld]); legend('H12','H16','H20');
plot(fc,tc20_cld-tc12_cld,fc,tc20_cld-tc16_cld); legend('H20-H12','H20-H16','location','best'); title('gaussian FWHM = 1 cm-1 CLD')

figure(2); 
plot(fc,tc12_clr-tc12_cld,fc,tc16_clr-tc16_cld,fc,tc20_clr-tc20_cld); legend('H12','H16','H20','location','best'); title('gaussian FWHM = 1 cm-1 CLR-CLD')
