echo " "
echo "to start"
echo "  (1) basic : only compile and link newly edited files, type               gfortran_compile_onefile.sc"
echo "  (2) from scratch, except keep SAME ../INCLUDE files,  type   make clean; gfortran_compile_onefile.sc"
echo "  (3) completely anew, including NEW ../INCLUDE files,  type   make clean; make"
echo " "

#echo "in the disort codes (scatter_disort_main.f90, scatter_disort_flux.f90, scatter_disort_code.f90) look for "
#echo "        print *,'dag, gfortran cannot compile these and says they should all be double' "
#echo "        call dostop "
#echo " and scatter_disort_code.f90 : fix for ifort (I have not figured out the automatic -fdble for gfortran) "
#echo " change all these back to single "
#echo "    DBDREF( DBLE(WVNMLO), DBLE(WVNMHI), DBLE(CMU(IQ)), DBLE(CMU(JQ)), & "
#echo "                    DBLE(PI*GMU(K)) ) * COS( MAZIM*PI*GMU( K ) ) "
#echo " "

echo "in kcartamain.f90 look for "
echo "        print *,'dag, gfortran cannot compile these and says they should all be double' "
echo "        c      use ifport             ! for getenv "
echo " and incomment for ifort (I have not figured out the automatic -fdble for gfortran) "
echo " see https://community.intel.com/t5/Intel-Fortran-Compiler/Errors-using-GETPID/m-p/1094827"
echo " "

make -f Makefile_v122_gfortran  scat
ls -lt ../BIN/kcarta.x90

#mv ../BIN/kcarta.x90 /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20
#ls -lt /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20

mv ../BIN/kcarta.x90 ../BIN/kcarta.x_f90_122_400ppmv_H20
ls -lt ../BIN/kcarta.x_f90_122_400ppmv_H20
