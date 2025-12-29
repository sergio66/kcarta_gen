echo "making kCARTA for H08, H12, G15, H16"

make
make -f makefile 400_H20_default_f90_orig605_805res
make -f makefile 400_H16_default_f90_orig605_805res
make -f makefile 400_G15_NLTEH16_f90_orig605_805res

make -f makefile 385_H12_CO2_UMBC_default_f90
make -f makefile 385_H12_CO2_UMBC_default_f90_X     ## Dec 2025, when we "lost" the disks

make -f makefile 385_H08_CO2_UMBC_default_f90

########################################################################

make -f makefile 400_H20_default_f90
make -f makefile 400_H16_default_f90

########################################################################

make -f makefile 400_H12p8_highres605_1205
make -f makefile 400_H12p8_veryhighres605_1205
