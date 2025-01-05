addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD

clc; disp('checking the profiles')

[h,ha,pM0_Feb23,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_era_Feb23_2021.rp.rtp');
[h,ha,pMR_Feb23,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_resetTWV_era_Feb23_2021.rp.rtp');
[h,ha,pN0_Feb23,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_era_Feb23_2021.rp.rtp');
[h,ha,pNR_Feb23,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_resetTWV_era_Feb23_2021.rp.rtp');

[h,ha,pM0_Mar17,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_era_Mar17_2021.rp.rtp');
[h,ha,pMR_Mar17,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/masuda_resetTWV_era_Mar17_2021.rp.rtp');
[h,ha,pN0_Mar17,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_era_Mar17_2021.rp.rtp');
[h,ha,pNR_Mar17,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/nalli_resetTWV_era_Mar17_2021.rp.rtp');

compare_profs(pM0_Feb23,pM0_Mar17,-1,'Masuda 0'); disp('ret to continue'); pause; 
compare_profs(pMR_Feb23,pMR_Mar17,-1,'Masuda R'); disp('ret to continue'); pause; 
compare_profs(pN0_Feb23,pN0_Mar17,-1,'Nalli  0'); disp('ret to continue'); pause; 
compare_profs(pNR_Feb23,pNR_Mar17,-1,'Nalli  R'); disp('ret to continue'); pause; 

mmw = mmwater_rtp(h,pNR_Mar17);
scatter_coast(pNR_Mar17.rlon,pNR_Mar17.rlat,10,pNR_Mar17.satzen); colormap jet
oo = find(abs(pNR_Mar17.rlat) <= 30 & pNR_Mar17.satzen > 30);
oo = find(abs(pNR_Mar17.rlat) <= 30 & pNR_Mar17.satzen > 30 & pNR_Mar17.stemp > 295 & mmw > 30);
oo = oo(1)    %% 16   
              %% now run kCARTA and save the nml files, rad files so you can dump out what is needed
%{
edit as needed
-rw-rw-r-- 1 sergio pi_strow 17572 Mar 25 17:56 set_rtp.m
-rw-rw-r-- 1 sergio pi_strow 15359 Mar 25 17:56 set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm.m

then run kcarta
 1017  25/03/21 17:55:27 mv /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/run_nml16 Mar17_2021/DUMP_SURFACE/run_nml16_Masuda_ResetTWV
 1018  25/03/21 17:56:20 mv /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/rad.dat16 Mar17_2021/DUMP_SURFACE/rad.dat16_Masuda_ResetTWV


 1021  25/03/21 17:58:49 mv /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/rad.dat16 Mar17_2021/DUMP_SURFACE/rad.dat16_Nalli_ResetTWV
 1022  25/03/21 17:59:00 mv /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/run_nml16 Mar17_2021/DUMP_SURFACE/run_nml16_Nalli_ResetTWV

then edit the .nml files for
3/15/21         iaaOverrideDefault(3,4) = -1 default
                iaaOverrideDefault(3,4) = +1 allows rad_main debug surface radiation printout to StdErr

      raInten = raSurface*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 5',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if

and rerun to get the dumps
cd Mar17_2021/DUMP_SURFACE
  run_kc.sc

grep -i 'DBGSURF' out_masudaemis.txt > out_masudaemis2.txt
grep -i 'DBGSURF' out_nalliemis.txt > out_nalliemis2.txt

then edit those two files to get ri of first 10 columns and read into matlab

surface_nalli = load('Mar17_2021/DUMP_SURFACE/out_nalliemis2.txt');
surface_masuda = load('Mar17_2021/DUMP_SURFACE/out_masudaemis2.txt');
fr = surface_nalli(:,1);
ii=2; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('Emissivity')
ii=3; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('Rho Solar')
ii=4; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('1 or 1/pi refl fac')
ii=5; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('ttorad(v,Ts)')
ii=6; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('BackGnd Thermal')
ii=7; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('Solar at surface')
ii=8; plot(fr,surface_masuda(:,ii),'b.-',fr,surface_nalli(:,ii),'r'); title('total I\_0 at surface')
ii=8; plot(fr,surface_masuda(:,ii)-surface_masuda(:,5),'b.-',fr,surface_nalli(:,ii)-surface_nalli(:,5),'r'); title('Sol+ReflTh surface')
  hl = legend('masuda','nalli','location','best'); 
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; disp('checking the kcarta calcs')

%% I had messed up the masuda reset calcs here
a1 = load('nalli_Mar17_2021_orig.mat');
a0 = load('nalli_Mar17_2021.mat');

%%% the Nalli ERA and Reset for a0 was totally messed
a1 = load('nalli_Feb23_2021.mat');
a0 = load('nalli_Feb23_2021_oops_messedupNalliEmis_Kcarta_settings.mat');

%% evey thing = 0 shows that I did not change things from Mar22 to current version, good
a1 = load('nalli_Mar17_2021_savedMar22_looksgood.mat');
a0 = load('nalli_Mar17_2021.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%
%% shows that a1 (Feb 23) is now quite similar to a0 (Mar 17) : good,the small diffs due to relefctivity rho being tweaked
a1 = load('nalli_Feb23_2021.mat');
a0 = load('nalli_Mar17_2021.mat');

%% shows that fioles on ftp site are now correct
a1 = load('/asl/ftp/pub/sergio/NalliEmiss/nalli_Feb23_2021.mat');
a0 = load('/asl/ftp/pub/sergio/NalliEmiss/nalli_Mar17_2021.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%

sum(sum(a0.results.rcalc_nalli_ResetERA-a1.results.rcalc_nalli_ResetERA))
sum(sum(a0.results.rcalc_nalli_ERA-a1.results.rcalc_nalli_ERA))
sum(sum(a0.results.rcalc_masuda_ResetERA-a1.results.rcalc_masuda_ResetERA))
sum(sum(a0.results.rcalc_masuda_ERA-a1.results.rcalc_masuda_ERA))
