#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_code.f90
#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_main.f90
#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_flux.f90

make -f Makefile_v122_Intelf90 scat
ls -lt ../BIN/kcarta.x90
mv ../BIN/kcarta.x90 /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20
ls -lt /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20
