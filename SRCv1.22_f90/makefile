#echo "IMPORTANT -- make sure you ALWAYS run cp_param_files_to_f90_v122.sc AFTER doing the datafixXYZ_H20AB"

$(info   )
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info IMPORTANT (1) -- try to build on oldest nodes rather than taki2 so it runs on every node)
$(info IMPORTANT (2) -- make sure you ALWAYS run cp_param_files_to_f90_v122.sc AFTER doing the datafixXYZ_H20AB, if you do things by hand)
$(info IMPORTANT (3) -- since this is 500-880,805-2830 cm-1 res use cp_param_files_to_f90_v122.sc and not cp_param_files_to_f90_v120.sc)
$(info IMPORTANT (4) -- if you want to use eg g77 change the make -f Makefile_v122_Intelf90 scat to eg  make -f Makefile_v120_gnu scat)
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info   )

400_H20_default_f90_short:
	make -f Makefile_v122_Intelf90 datafix400_v122_H2020
	echo "this datafix400_v122_H2020 looks suspiciously similar to datafix400_v122_H2020_orig605_805res"
	cp_param_files_to_f90_v122.sc
	make -f Makefile_v122_Intelf90 scat
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_H20

400_H20_default_f90:
	make -f Makefile_v122_Intelf90 datafix400_v122_H2020
	echo "this datafix400_v122_H2020 looks suspiciously similar to datafix400_v122_H2020_orig605_805res"
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_400ppmv_H20
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_400ppmv_H20
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_H20

#########################################################################

400_H16_default_f90:
	make -f Makefile_v122_Intelf90 datafix400_v122_H2016
	echo "this datafix400_v122_H2016 looks suspiciously similar to datafix400_v122_H2016_orig605_805res"
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_400ppmv_H16
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_400ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_H16

mars_H16_default_f90:
	make -f Makefile_v122_Intelf90 mars_datafix400_v122_H2016
	cp_param_files_to_f90_mars.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta_mars.x_f90_122_400ppmv_H16
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT_mars.x_f90_122_400ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta_mars.x_f90_122_400ppmv_H16

385_H16_default_f90:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2016
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H16
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H16

400_H12p8_highres605_1205:
	make -f Makefile_v122_Intelf90 datafix400_v122_H2016_highres
	cp_param_files_to_f90HIGHRES.sc
	echo "CO2 is from LBLRTM12.8"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90 ../BIN/bkcarta.x_f90_122_highres605_1205_12p8
	mv ../BIN/kcarta.x90 ../BIN/kcarta.x_f90_122_highres605_1205_12p8	

400_H12p8_veryhighres605_1205:
	make -f Makefile_v122_Intelf90 datafix400_v122_H2016_veryhighres
	cp_param_files_to_f90VERYHIGHRES.sc
	echo "CO2 is from LBLRTM12.8"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90 ../BIN/bkcarta.x_f90_122_veryhighres605_1205_12p8
	mv ../BIN/kcarta.x90 ../BIN/kcarta.x_f90_122_veryhighres605_1205_12p8	

########################################################################
## these uses ../INCLUDE/pre_defined_orig605_805res.param which is
##   500-605 cm-1 at 0.0015 cm-1 res, 605-2830 cm-1 at 0.0025 cm-1

400_H20_default_f90_orig605_805res:
	make -f Makefile_v122_Intelf90 datafix400_v122origres_H2020
	echo "this datafix400_v122_H2020 looks suspiciously similar to datafix400_v122_H2020_orig605_805res"
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_400ppmv_H20_orig605_805res
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_400ppmv_H20_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_H20_orig605_805res

400_H16_default_f90_orig605_805res:
	make -f Makefile_v122_Intelf90 datafix400_v122origres_H2016
	echo "this datafix400_v122_H2016 looks suspiciously similar to datafix400_v122_H2016_orig605_805res"
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_400ppmv_H16_orig605_805res
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_400ppmv_H16_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_H16_orig605_805res

400_G15_NLTEH16_f90_orig605_805res:
	make -f Makefile_v122_Intelf90 datafix400_v122_G2015_NLTEH2016
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_400ppmv_G15
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_400ppmv_G15
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_400ppmv_G15

#########################
## these are orig res, and 385 pmv only
## Dec 2025, when we "lost" the disks
385_H12_CO2_UMBC_default_f90:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2012_CO2_UMBC
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H12_CO2_UMBC
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H12_CO2_UMBC
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H12_CO2_UMBC

## these are orig res, and 385 pmv only
## Dec 2025, when we "lost" the disks
385_H12_CO2_UMBC_default_f90_X_noneed:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2012_CO2_UMBC_X
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H12_CO2_UMBC_X
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H12_CO2_UMBC_X
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H12_CO2_UMBC_X
#########################


## these are orig res, and 385 pmv only
385_H08_CO2_UMBC_default_f90:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2008_CO2_UMBC
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H08_CO2_UMBC
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H08_CO2_UMBC
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H08_CO2_UMBC

## these are orig res, and 385 pmv only
x_NOexist_385_H04_CO2_UMBC_default_f90:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2004_CO2_UMBC
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H04_CO2_UMBC
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H04_CO2_UMBC
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H04_CO2_UMBC

## these are orig res, and 385 pmv only
x_NOexist_385_H2k_CO2_UMBC_default_f90:
	make -f Makefile_v122_Intelf90 datafix385_v122_H2k_CO2_UMBC
	cp_param_files_to_f90_v122.sc
#	make -f Makefile_v122_Intelf90 basic
#	make -f Makefile_v122_Intelf90 basicTXT	
	make -f Makefile_v122_Intelf90 scat
#	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_122_385ppmv_H2k_CO2_UMBC
#	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_122_385ppmv_H2k_CO2_UMBC
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_122_385ppmv_H2k_CO2_UMBC

clean: 
	rm *.optrpt *genmod.* *.o *.mod *~ 2> /dev/null || true
########################################################################
