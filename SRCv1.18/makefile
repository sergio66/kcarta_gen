385_default:
	make -f Makefile datafix385_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_385ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_385ppmv

385_lblrtm12.4:
	make -f Makefile datafix385_H2012_LBLRTM
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_lblrtm12.4
	mv ../BIN/kcarta.x ../BIN/kcarta.x_lblrtm12.4	

385_highres605_1205:
	/bin/rm ../INCLUDE/kcarta.param
	ln -s ../INCLUDE/kcarta_v118HIGHRES_IR.param ../INCLUDE/kcarta.param
	echo "MAKE SURE rtp_interface.f f1,f2 are 605, 1205 (current limits of HIGHRES database) !!!"
	read -p "Press [Enter] key to continue ..."
	make -f Makefile clean	
	make -f Makefile
	make -f Makefile scat	
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_highres605_1205
	mv ../BIN/kcarta.x ../BIN/kcarta.x_highres605_1205	

########################################################################

400_default:
	make -f Makefile datafix400_H2012
	make -f Makefile
	make -f Makefile scat		
	mv ../BIN/bkcarta.x ../BIN/bkcarta.x_400ppmv
	mv ../BIN/kcarta.x ../BIN/kcarta.x_400ppmv

400_default_f90:
	make -f Makefile_v118_Intelf90 datafix400_H2012
	/bin/cp ../INCLUDE/scatter.param  ../INCLUDE/scatter.f; fixcon ../INCLUDE/scatter.f ../INCLUDE/scatter.f90
	sed -i 's/.param/.f90/g' ../INCLUDE/scatter.f90 ../INCLUDE/scatter.f90
	/bin/cp ../INCLUDE/kcarta.param  ../INCLUDE/kcarta.f;   fixcon ../INCLUDE/kcarta.f ../INCLUDE/kcarta.f90
	sed -i 's/.param/.f90/g' ../INCLUDE/kcarta.f90 ../INCLUDE/kcarta.f90	
	/bin/cp ../INCLUDE/pre_defined.param  ../INCLUDE/pre_defined.f;   fixcon ../INCLUDE/pre_defined.f ../INCLUDE/pre_defined.f90
	/bin/cp ../INCLUDE/post_defined.param  ../INCLUDE/post_defined.f;   fixcon ../INCLUDE/post_defined.f ../INCLUDE/post_defined.f90			
	make -f Makefile_v118_Intelf90
	make -f Makefile_v118_Intelf90 scat
	make -f Makefile_v118_Intelf90 parallel	
	mv ../BIN/bkcarta.x        ../BIN/bkcarta.x_f90_400ppmv
	mv ../BIN/bkcartaTXT.x     ../BIN/bkcartaTXT.x_f90_400ppmv
	mv ../BIN/kcarta.x         ../BIN/kcarta.x_f90_400ppmv
	mv ../BIN/kcartaparallel.x ../BIN/kcartaparallel.x_f90_400ppmv	
