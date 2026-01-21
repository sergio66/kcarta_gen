400_H16_default:
	make -f Makefile_tar_objs_data datafix400_H2016
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_f77_v1.18_400ppmv_H16_NLTEH16_120NLTElayers_new
	mv ../BIN/kcarta.x_v1.18  ../BIN/kcarta.x_f77_v1.18_400ppmv_H16_NLTEH16_120NLTElayers_new
#	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_400ppmv_H16
#	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_400ppmv_H16

400_H16_NLTEH2008:
	make -f Makefile_tar_objs_data datafix400_H2016_NLTEH2008
	make -f Makefile
	make -f Makefile scat
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_400ppmv_H16
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_400ppmv_H16			  

400_H12_default:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_400ppmv
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_400ppmv

400_H12_default_test:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_400ppmv_test
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_400ppmv_test

########################################################################
## this is nothing special, it basically is SRCv1.18 with the strutures dot
## replaced by percent sign 
400_H12_default:
	make -f Makefile_v118_Intelf90 datafix400_H2012
	make -f Makefile_v118_Intelf90
	make -f Makefile_v118_Intelf90 scat
	mv ../BIN/bkcarta.x_v1.18        ../BIN/bkcarta.x_v1.18_400ppmv
	mv ../BIN/bkcartaTXT.x_v1.18     ../BIN/bkcartaTXT.x_v1.18_400ppmv
	mv ../BIN/kcarta.x_v1.18         ../BIN/kcarta.x_v1.18_400ppmv

########################################################################

#########################
## these are orig res, and 385 pmv only
## Dec 2025, when we "lost" the disks
385_H12_default:
	echo "making 385_H12_default"
	make -f Makefile datafix385_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_385ppmv
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_385ppmv
#########################

385_H08_default:
	make -f Makefile datafix385_H2008
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_385ppmv_H08
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_385ppmv_H08

385_H12_lblrtm12.4:
	make -f Makefile datafix385_H2012_LBLRTM
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_lblrtm12.4
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_lblrtm12.4	

385_H12_highres605_1205:
	/bin/rm ../INCLUDE/kcarta.param
	ln -s ../INCLUDE/kcarta_v118HIGHRES_IR.param ../INCLUDE/kcarta.param
	echo "MAKE SURE rtp_interface.f f1,f2 are 605, 1205 (current limits of HIGHRES database) !!!"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x_v1.18 ../BIN/bkcarta.x_v1.18_highres605_1205
	mv ../BIN/kcarta.x_v1.18 ../BIN/kcarta.x_v1.18_highres605_1205	

clean: 
	rm *.o *~ 2> /dev/null || true

########################################################################

