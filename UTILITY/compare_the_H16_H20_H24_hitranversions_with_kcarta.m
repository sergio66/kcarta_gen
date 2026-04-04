%% this is copied into /home/sergio/KCARTA/UTILITY
%% it came from running /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/sergio_matlab_chip.sbatch with
%%   set_rtp.m    iInstr = 1; iDoConvolve = 1;
%%                use_this_rtp = '/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1013_400ppm_unitemiss.op.rtp';
%% for KCARTA H16,H20,H24 (see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_kcarta_exec_iHITRAN.m)
%% and then saving the convolved AIRS SRF output intho the subdirectories H2020_CKD32/ H2024_CKD32/ H2024_CKD43/
%%
%%    kcartaexec16 = '/home/sergio/KCARTA/BIN/kcarta.x90_v1.22_400ppmv_H16';                %% H15, v1.22, UMBC CO2
%%    kcartaexec20 = '/home/sergio/KCARTA/BIN/kcarta.x90_v1.22_400ppmv_H20';                %% H20, v1.22, UMBC CO2
%%    kcartaexec24 = '/home/sergio/KCARTA/BIN/kcarta.x90_v1.22_400ppmv_H24';                %% H24, v1.22, UMBC CO2
%%
%% the results show after convolution that if you use same CKD3.2, then  H20 and H24 are basically identical except at 1113 cm-1 and 1247 cm1 ... an it depends on SKT/MMW
%%   so probably water
%%   Will analyze by looping over removing gases, using /home/sergio/KCARTA/UTILITY/clust_run_kcarta_H2020_H2024.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/git/matlabcode
addpath /home/sergio/git/matlabcode/PLOTTER

dir0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/TEST_H2020_H2024_CKD32_CKD43/';
dirH24_43N2O2_CO2CH4 = [dir0 'H2024_CKD43_newN2O2_newCO2CH4/'];
dirH24_43N2O2        = [dir0 'H2024_CKD43_newN2O2/'];
dirH24_43 = [dir0 'H2024_CKD43/'];
dirH24_32 = [dir0 'H2024_CKD32/'];
dirH20_32 = [dir0 'H2020_CKD32/'];

for ii = 1 : 49
  fin = [dirH24_43N2O2_CO2CH4 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fin);
  d2024_43B(ii,:) = x.rKc;

  fin = [dirH24_43N2O2 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fin);
  d2024_43A(ii,:) = x.rKc;

  fin = [dirH24_43 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fin);
  d2024_43(ii,:) = x.rKc;

  fin = [dirH24_32 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fin);
  d2024_32(ii,:) = x.rKc;

  fin = [dirH20_32 '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  x = load(fin);
  d2020_32(ii,:) = x.rKc;

end

fKc = x.fKc;

t2024_43B = rad2bt(fKc,d2024_43B');
t2024_43A = rad2bt(fKc,d2024_43A');
t2024_43  = rad2bt(fKc,d2024_43');
t2024_32  = rad2bt(fKc,d2024_32');
t2020_32  = rad2bt(fKc,d2020_32');

figure(1);
plot(fKc,mean(t2024_43'-t2024_32'),'b',fKc,mean(t2024_43'-t2020_32'),'g',fKc,mean(t2024_32'-t2020_32'),'r')
xlim([645 1645]); ylabel('Mean')
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H20 CKD32','H24 CKD32 - H20 CKD32','location','best');

figure(2);
plot(fKc,std(t2024_43'-t2024_32'),'b',fKc,std(t2024_43'-t2020_32'),'g',fKc,std(t2024_32'-t2020_32'),'r')
xlim([645 1645]); ylabel('Std');
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H20 CKD32','H24 CKD32 - H20 CKD32','location','best');

figure(3)
i900 = find(fKc >= 900,1);
plot(t2024_43(i900,:),t2024_43(i900,:)-t2024_32(i900,:),'b.',t2024_43(i900,:),t2020_32(i900,:)-t2020_32(i900,:),'gx',t2024_43(i900,:),t2024_43(i900,:)-t2020_32(i900,:),'ro')
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H20 CKD32','H24 CKD32 - H20 CKD32','location','best');
xlabel('BT900 H24 CKD43')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);
plot(fKc,mean(t2024_43'-t2024_32'),'b',fKc,mean(t2024_43'-t2024_43A'),'g',fKc,mean(t2024_43'-t2024_43B'),'r')
xlim([645 1645]); ylabel('Mean')
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H24 CKD43 N2O2','H24 CKD43 - H24 CKD43 N2O2 CO2CH4','location','best');

figure(5);
plot(fKc,std(t2024_43'-t2024_32'),'b',fKc,std(t2024_43'-t2024_43A'),'g',fKc,std(t2024_43'-t2024_43B'),'r')
xlim([645 1645]); ylabel('Std')
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H24 CKD43 N2O2','H24 CKD43 - H24 CKD43 N2O2 CO2CH4','location','best');

figure(6)
i900 = find(fKc >= 900,1);
plot(t2024_43(i900,:),t2024_43(i900,:)-t2024_32(i900,:),'b.',t2024_43(i900,:),t2024_43(i900,:)-t2024_43A(i900,:),'gx',t2024_43(i900,:),t2024_43(i900,:)-t2024_43B(i900,:),'ro')
legend('H24 CKD43 - H24 CKD32','H24 CKD43 - H24 CKD43 N2O2','H24 CKD43 - H24 CKD43 N2O2 CO2CH4','location','best');
xlabel('BT900 H24 CKD43')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/matlibSergio/matlib2025/h4tools
[h,ha,p,pa] = rtpread('/home/sergio/git/matlabcode/REGR_PROFILES_SARTA/REGR49_PROFILES_for_kCARTA_breakouts_for_SARTA/regr49_1013_400ppm_unitemiss.op.rtp');
mmw = mmwater_rtp(h,p);

figure(7);
plot(mmw,t2024_43(i900,:)-t2024_32(i900,:),'b.',mmw,t2024_32(i900,:)-t2020_32(i900,:),'gx',mmw,t2024_43(i900,:)-t2020_32(i900,:),'ro')
xlabel('mmw'); ylabel('BTD (K)')
ax = axis;
line([mmw(01) mmw(01)],[ax(3) ax(4)],'color','r');
line([mmw(49) mmw(49)],[ax(3) ax(4)],'color','b');
plotaxis2;
legend('H24 CKD43 - H24 CKD32','H24 CKD32 - H20 CKD32','H24 CKD43 - H20 CKD32','location','best');
text(mmw(01)+1,1,'TRP','color','r')
text(mmw(49)+1,1,'STD','color','b')
set(gca,'fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1111 = find(fKc >= 1113.55,1);
i1273 = find(fKc >= 1273,1);

figure(8);
plot(1:49,t2024_32(i1111,:)-t2020_32(i1111,:),'b.',1:49,t2024_32(i1273,:)-t2020_32(i1273,:),'rx');
legend('1111 cm-1','1272 cm-1','location','best');
title('H24 CKD32 = H20 CKD 32')

btd = t2024_32(i1111,:)-t2020_32(i1111,:);
worst = find(btd == max(abs(btd)));
  fprintf(1,'suggest you loop through gases in profile %2i to find bump at 1113.5 cm-1 \n',worst)
  fprintf(1,' which has stemp = %8.3f K mmw = %8.3f mm \n',p.stemp(worst),mmw(worst))
leastbad = find(btd == min(abs(btd)));
  fprintf(1,'suggest you loop through gases in profile %2i to find bump at 1113.5 cm-1 \n',leastbad)
  fprintf(1,' which has stemp = %8.3f K mmw = %8.3f mm \n',p.stemp(leastbad),mmw(leastbad))
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fig 18 of Mlawer, E. J., Mascio, J., Turner, D. D., Payne, V. H., Flynn, C. J., & Pincus, R.(2024).
%% A more transparent infrared window. Journal of Geophysical Research: Atmospheres, 129,
%% e2024JD041366. https://doi.org/10.1029/2024JD041366

[Y,I] = sort(fKc);
figure(9);
ind = [1 2 3 4 5 49];
plot(fKc(I),t2024_43(I,ind)-t2024_32(I,ind),'linewidth',2)
ylabel('H24 CKD43 - H24 CKD32'); xlabel('Wavenumer cm-1')
legend('TRP','MLS','MLW','SAS','SAW','STD','location','best')
xlim([700 1300]);
set(gca,'fontsize',14)

mycolors = [
    1 0 0;   % Red
    0 0 1;   % Blue
    0.5 0.5 0.5; % Gray
    1 0.5 0; % Orange
    1 0.7 1; % Pink
    0 0 0    % Black
];

% Set the color order for the current axes (gca) or figure
colororder(gca, mycolors); 

