#echo "IMPORTANT -- make sure you ALWAYS run cp_param_files_to_f90.sc AFTER doing the datafixXYZ_H20AB"

$(info   )
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info IMPORTANT (1) -- try to build on oldest nodes rather than taki2 so it runs on every node)
$(info IMPORTANT (2) -- make sure you ALWAYS run cp_param_files_to_f90.sc AFTER doing the datafixXYZ_H20AB, if you do things by hand)
$(info IMPORTANT (3) -- since this is CONSTANT 605-2830 cm-1 res use cp_param_files_to_f90_v120.sc and not cp_param_files_to_f90.sc)
$(info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>)
$(info   )

400_H16_default_f90:
	make -f Makefile_v120_Intelf90 datafix400_v120_H2016
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT	
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H16
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H16
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H16

400_H16_g110_f90:
	make -f Makefile_v120_Intelf90 datafix400_v120_H2016_g110
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT	
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H16_g110
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H16_g110
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H16_g110

400_G15_NLTEH16_f90:
	make -f Makefile_v120_Intelf90 datafix400_G2015_NLTEH2016
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT	
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_G15
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_G15
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_G15

400_H16_NLTEH08_f90:
	make -f Makefile_v120_Intelf90 datafix400_v120_H2016_NLTEH2008
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT	
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H16_NLTE08
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H16_NLTE08
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H16_NLTE08

400_H16_NLTEH04_f90:
	make -f Makefile_v120_Intelf90 datafix400_v120_H2016_NLTEH2004
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT	
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H16_NLTE04
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H16_NLTE04
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H16_NLTE04

400_H12_default_f90:
	make -f Makefile_v120_Intelf90 datafix400_v120_H2012
	cp_param_files_to_f90_v120.sc
	make -f Makefile_v120_Intelf90 basic
	make -f Makefile_v120_Intelf90 basicTXT		
	make -f Makefile_v120_Intelf90 scat
	mv ../BIN/bkcarta.x90        ../BIN/bkcarta.x_f90_120_400ppmv_H12
	mv ../BIN/bkcartaTXT.x90     ../BIN/bkcartaTXT.x_f90_120_400ppmv_H12
	mv ../BIN/kcarta.x90         ../BIN/kcarta.x_f90_120_400ppmv_H12

########################################################################
### these are basically "old" SRCv1.18 and may not exist anymore
385_H12_default:
	make -f Makefile datafix385_v120_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_385ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_385ppmv

385_H12_lblrtm12.4:
	make -f Makefile datafix385_v120_H2012_LBLRTM
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
	make -f Makefile datafix400_v120_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv

clean: 
	rm *.optrpt *genmod.* *.o *.mod *~ 2> /dev/null || true

########################################################################
