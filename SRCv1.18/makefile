400_H16_default:
	make -f Makefile_tar_objs_data datafix400_H2016
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv_H16
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv_H16

400_H16_NLTEH2008:
	make -f Makefile_tar_objs_data datafix400_H2016_NLTEH2008
	make -f Makefile
	make -f Makefile scat
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv_H16
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv_H16			  

400_H12_default:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv

400_H12_default_test:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv_test
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv_test

########################################################################
## this is nothing special, it basically is SRCv1.18 with the strutures dot
## replaced by percent sign 
400_H12_default_f90:
	make -f Makefile_v118_Intelf90 datafix400_H2012
	make -f Makefile_v118_Intelf90
	make -f Makefile_v118_Intelf90 scat
	mv ../BIN/bkcarta.x        ../BIN/bkcarta.x_f90_400ppmv
	mv ../BIN/bkcartaTXT.x     ../BIN/bkcartaTXT.x_f90_400ppmv
	mv ../BIN/kcarta.x         ../BIN/kcarta.x_f90_400ppmv

########################################################################

385_H12_default:
	make -f Makefile datafix385_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_385ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_385ppmv

385_H08_default:
	make -f Makefile datafix385_H2008
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_385ppmv_H08
	mv ../BIN/kcarta.x ../BIN/kcarta.x_385ppmv_H08

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

clean: 
	rm *.o *~ 2> /dev/null || true

########################################################################

