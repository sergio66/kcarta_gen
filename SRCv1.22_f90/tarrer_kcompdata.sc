#!/bin/bash
KCOMPDATA_H2016='/asl/ftp/pub/packages/water_hdo_etc_H2016_IR_v1.ieee-le.tar'
KCOMPDATA_LBLRTM12p8='/asl/ftp/pub/packages/co2_ch4_lblrtm12p8.tar'

touch $KCOMPDATA_H2016
#rm ${KCOMPDATA_H2016} 2> /dev/null || true;	\
echo "tar -cvf ${KCOMPDATA_H2016}			\
  /asl/rta/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/*.dat	   \
  /asl/rta/kcarta/H2016.ieee-le/IR605/etc.ieee-le/*.dat    \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*.dat   \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*32.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*25.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor1.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor1.bin" 
tar -cf ${KCOMPDATA_H2016}			\
  /asl/rta/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/*.dat	   \
  /asl/rta/kcarta/H2016.ieee-le/IR605/etc.ieee-le/*.dat    \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*.dat   \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CT-N2*  \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*32.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*25.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor1.bin \
  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/CKDFor1.bin 

touch $KCOMPDATA_LBLRTM12p8
#rm ${KCOMPDATA_LBLRTM12p8} 2> /dev/null || true;	\
echo "tar -cvf ${KCOMPDATA_LBLRTM12p8}			\
  /asl/rta/kcarta_sergio/KCDATA/General/ChiFile            \
  /asl/rta/kcarta_sergio/UMBC_CO2_H1998.ieee-le/CO2ppmv400.ieee-le        \
  /asl/rta/kcarta_sergio/KCDATA/NLTE/SARTA_COEFS/nonLTE7_m150.le.dat      \
  /asl/data/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/CO2_400ppmv \
  /asl/data/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/*g6.dat"
tar -cf ${KCOMPDATA_LBLRTM12p8}			\
  /asl/rta/kcarta_sergio/KCDATA/General/ChiFile            \
  /asl/rta/kcarta_sergio/UMBC_CO2_H1998.ieee-le/CO2ppmv400.ieee-le        \
  /asl/rta/kcarta_sergio/KCDATA/NLTE/SARTA_COEFS/nonLTE7_m150.le.dat      \
  /asl/data/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/CO2_400ppmv \
  /asl/data/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/*g6.dat
#chmod 664 ${KCOMPDATA_LBLRTM12p8}

# just use MT CKD3.2 its the latest, and no need for the compressed CKD files ....
#  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/*25.bin \ 
#  /asl/rta/kcarta_sergio/KCDATA/General/CKDieee_le/Compr/*/*   \
### OH OH need to make /asl/rta/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/r1005_g103.dat
### OH OH need to make /asl/rta/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/r1005_g103.dat
### OH OH need to make /asl/rta/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/r1005_g103.dat

### these are NLTE, only giving /asl/rta/kcarta_sergio/KCDATA/NLTE/SARTA_COEFS/nonLTE7_m150.le.dat
#   		 /asl/rta/kcarta/v10.ieee-le/etc.ieee-le \
#                /asl/rta/kcarta_sergio/KCDATA/General/refgas2_400ppmv \
#                /asl/rta/kcarta_sergio/KCDATA/General/refgas2Back_400ppmv \
#                /asl/rta/kcarta_sergio/KCDATA/General/refgas2BackUA_400ppmv \
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/UA/ \
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/LA/ \
#                /asl/data/kcarta_sergio/KCDATA/NLTE/AIRSCO2/CO2_H16/co2_allpqr_
#                /asl/rta/kcarta_sergio/BACKGND_COUSIN/fbin/etc.ieee-le/400ppmv_H16/
#                /asl/rta/kcarta_sergio/BACKGND_COUSIN_UPPER/fbin/etc.ieee-le/400ppmv_H16/
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/AIRSCO2/LINEMIX_H16_400/
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/UA/us_std_prof_ua_400ppmv
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/LA/us_std_prof_la_400ppmv
#                /asl/rta/kcarta_sergio/KCDATA/NLTE/SARTA_COEFS/nonLTE7_m150.le.dat
