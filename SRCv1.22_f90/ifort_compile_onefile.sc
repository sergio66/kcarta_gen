echo " "
echo "to start"
echo "  (1) basic : only compile and link newly edited files, type               ifort_compile_onefile.sc"
echo "  (2) from scratch, except keep SAME ../INCLUDE files,  type   make clean; ifort_compile_onefile.sc"
echo "  (3) completely anew, including NEW ../INCLUDE files,  type   make clean; make"
echo " "

#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_code.f90
#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_main.f90
#ifort    -c -fpp -implicitnone -msse2 -extend-source 132 -O2 -names lowercase -assume nounderscore -heap-arrays -diag-disable 8291 -mcmodel=medium -shared-intel   scatter_pclsam_flux.f90

make -f Makefile_v122_Intelf90 scat
ls -lt ../BIN/kcarta.x90
mv ../BIN/kcarta.x90 /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20
ls -lt /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20
