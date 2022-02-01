rm hres.dat vhres.dat usualres.dat

### these are coded with LNLRTM12.8 default, so no need to "switch" .. points to same data!
time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_highres605_1205_12p8       quickuse_high_veryhigh_res.nml hres.dat
mv quickuse_warning.msg quickuse_warning_veryhighres.msg

### these are coded with LNLRTM12.8 default, so no need to "switch" .. points to same data!
time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_veryhighres605_1205_12p8   quickuse_high_veryhigh_res.nml vhres.dat
mv quickuse_warning.msg quickuse_warning_highres.msg

## here we do not switch to alternate datanase (ie usual res is computed using UMBC CO2 ods
#time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res quickuse_high_veryhigh_res.nml                usualres.dat
### but we HAVE to switch from kCARTA default UMBC to LBLRTM12.8
time /home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res quickuse_high_veryhigh_res_yesaltdatabase.nml usualres.dat
mv quickuse_warning.msg quickuse_warning_usualres.msg
