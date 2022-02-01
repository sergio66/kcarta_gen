# [sergio@taki-usr2 GENERIC_RADSnJACS_MANYPROFILES]$ ls -lt /home/sergio/KCARTA/BIN/kcarta.x_f90_120_400ppmv_H16
# -rwxrwxr-x 1 sergio pi_strow 4513720 Feb 14 05:19 /home/sergio/KCARTA/BIN/kcarta.x_f90_120_400ppmv_H16
# [sergio@taki-usr2 GENERIC_RADSnJACS_MANYPROFILES]$ ls -lt /home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res
# -rwxrwxr-x 1 sergio pi_strow 4739344 Oct 26 13:20 /home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res

time /home/sergio/KCARTA/BIN/kcarta.x_f90_120_400ppmv_H16 xrun_nml7 v120.dat
mv warning.msg warning120.msg

time /home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res xrun_nml7 v121.dat
mv warning.msg warning121.msg
