#echo "IMPORTANT -- make sure you ALWAYS run cp_param_files_to_f90.sc AFTER doing the datafixXYZ_H20AB"

$(info   )
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info IMPORTANT (1) -- try to build on oldest nodes rather than taki2 so it runs on every node)
$(info IMPORTANT (2) -- make sure you ALWAYS run cp_param_files_to_f90.sc AFTER doing the datafixXYZ_H20AB, if you do things by hand)
$(info IMPORTANT (3) -- since this is 500-880,805-2830 cm-1 res use cp_param_files_to_f90.sc and not cp_param_files_to_f90_v120.sc)
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info   )

400_H16_default_f90:
	make -f Makefile_v121_Intelf90 datafix400_v121_H2016
	cp_param_files_to_f90.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H16
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H16

385_H16_default_f90:
	make -f Makefile_v121_Intelf90 datafix385_v121_H2016
	cp_param_files_to_f90.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_385ppmv_H16
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_385ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_385ppmv_H16

400_H12p8_highres605_1205:
	make -f Makefile_v121_Intelf90 datafix400_v121_H2016_highres
	cp_param_files_to_f90HIGHRES.sc
	echo "CO2 is from LBLRTM12.8"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90 ../BIN/bkcarta.x_highres605_1205_12p8
	mv ../BIN/kcarta.x90 ../BIN/kcarta.x_highres605_1205_12p8	

400_H12p8_veryhighres605_1205:
	make -f Makefile_v121_Intelf90 datafix400_v121_H2016_veryhighres
	cp_param_files_to_f90VERYHIGHRES.sc
	echo "CO2 is from LBLRTM12.8"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90 ../BIN/bkcarta.x_veryhighres605_1205_12p8
	mv ../BIN/kcarta.x90 ../BIN/kcarta.x_veryhighres605_1205_12p8	

########################################################################
## these uses ../INCLUDE/pre_defined_orig605_805res.param which is
##   500-605 cm-1 at 0.0015 cm-1 res, 605-2830 cm-1 at 0.0025 cm-1
400_H16_default_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v121origres_H2016
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H16_orig605_805res
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H16_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H16_orig605_805res

400_H16_default_f90_orig605_805res_v120:
	make -f Makefile_v121_Intelf90 datafix400_v120origres_H2016
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H16_orig605_805res
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H16_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H16_orig605_805res

400_H12_default_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120origres_H2012
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H12_orig605_805res
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H12_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H12_orig605_805res

400_H08_default_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120origres_H2008
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H08_orig605_805res
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H08_orig605_805res
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H08_orig605_805res

#########################

400_H16_g110_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120_H2016_g110
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H16_g110
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H16_g110
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H16_g110

400_G15_NLTEH16_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120_G2015_NLTEH2016
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_G15
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_G15
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_G15

400_H16_NLTEH08_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120_H2016_NLTEH2008
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H16_NLTE08
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H16_NLTE08
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H16_NLTE08

400_H16_NLTEH04_f90_orig605_805res:
	make -f Makefile_v121_Intelf90 datafix400_v120_H2016_NLTEH2004
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v121_Intelf90 basic
	make -f Makefile_v121_Intelf90 basicTXT	
	make -f Makefile_v121_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_121_400ppmv_H16_NLTE04
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_121_400ppmv_H16_NLTE04
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_121_400ppmv_H16_NLTE04

########################################################################
### these are basically "old" SRCv1.18
385_H12_default:
	make -f Makefile datafix385_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_385ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_385ppmv

385_H12_lblrtm12.4:
	make -f Makefile datafix385_H2012_LBLRTM
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_lblrtm12.4
	mv ../BIN/kcarta.x ../BIN/kcarta.x_lblrtm12.4	

385_H12_highres605_1205:
	/bin/rm ../INCLUDE/kcarta.param
	ln -s ../INCLUDE/kcarta_v118HIGHRES_IR.param ../INCLUDE/kcarta.param
	echo "MAKE SURE rtp_interface.f f1,f2 are 605, 1205 (current limits of HIGHRES database) !!!"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_highres605_1205
	mv ../BIN/kcarta.x ../BIN/kcarta.x_highres605_1205	

400_H12_default:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv

clean: 
	rm *.optrpt *genmod.* *.o *.mod *~ 2> /dev/null || true
########################################################################
