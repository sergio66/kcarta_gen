echo "making kCARTA for H08, H12, G15, H16"

make -f makefile 385_H12_CO2_UMBC_default_f90
make -f makefile 385_H08_CO2_UMBC_default_f90
make -f makefile 400_H16_default_f90_orig605_805res
make -f makefile 400_G15_NLTEH16_f90_orig605_805res
