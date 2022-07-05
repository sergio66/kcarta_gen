if iHITRAN == 2008
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x.114';
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x.115';
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x.116';
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x.116_hartmann';
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta_H2008.x';
  kcartaexec = '/home/sergio/KCARTA/BIN/Oct10_2016/bkcarta_H2008.x';                  %% OLD
  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x_385ppmv_H08';                       %% NEWER
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H08_orig605_805res'; %% H08, v1.21, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_385ppmv_H08_CO2_UMBC';       %% H08, v1.22, UMBC CO2

elseif iHITRAN == 2012  
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv';                            %% can do cloudy calcs with this, H2012, v1.18
%  kcartaexec = '/home/sergio/KCARTA/BIN/bkcarta.x_10gasjac_10MP';                    %% do CLR calcs, H2012; can do 10 jacs, 10 sets of MP DNE
%  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv';                       %% the whole kit and caboodle. H12, v1.20
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv_test_Apr15_2016Absoft';
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv_test_Oct01_2016Ifort';
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv';                            %% latest gen v1.18 code, allows raAltComprDirsScale
                                                                                      %%   is UMBC CO2 messed up???
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv';                        %% latest gen v1.20 code, allows raAltComprDirsScale
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H12';                %% H12, v1.21, allows raAltComprDirsScale, set to one thread test
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H12_orig605_805res'; %% H12, v1.21, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_385ppmv_H12_CO2_UMBC';       %% H12, v1.22, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1

elseif iHITRAN == 2015
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv_G15';                    %% G15, v1.20, allows raAltComprDirsScale
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_G15';                %% G15, v1.21, allows raAltComprDirsScale, set to one thread test
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_G15';                %% G15, v1.22, allows raAltComprDirsScale, set to one thread test

elseif iHITRAN == 2016 | iHITRAN == 2016.3
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_400ppmv_H16';                        %% H16, v1.18, allows raAltComprDirsScale
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv_H16';                    %% H16, v1.20, allows raAltComprDirsScale
  %kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x90';                                 %% H16, v1.20, allows raAltComprDirsScale, set to one thread test

  %%%%%%%%%%
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16';                %% H16, v1.21, allows raAltComprDirsScale, set to one thread test .. for 500-880 only cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_120_400ppmv_H16';                %% H16, v1.20, should straight default to 0.0025 cm-1, 605-2830 cm-1 across board, kTemperVary = -1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res'; %% H16, v1.21, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res'; %% H16, v1.22, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1, sun-earth dist modification
  %%%%%%%%%%

elseif iHITRAN == 2017
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_highres605_1205_12p8';
elseif iHITRAN == 2018
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_veryhighres605_1205_12p8_error';

elseif iHITRAN == 2020
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20_orig605_805res';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if iDoFlux > 0 | iDoCloud > 0
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x';
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_400ppmv';    %% the whole kit and caboodle  
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16'; %% H16, v1.21, allows raAltComprDirsScale, set to one thread test .. for 500-880,805-2830 cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res'; %% H16, v1.21, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res'; %% H16, v1.22, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1, sun-earth dist modification
  kcartaexec = '/home/sergio/KCARTA/BIN/kcarta.x_f90_122_400ppmv_H20';                %% H20, v1.22, allows raAltComprDirsScale, set to one thread test .. for 605-2830 cm-1, sun-earth dist modification
end

