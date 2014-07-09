#
# makefile to build kcarta.x
#
# "make" or "make install" will compile kcarta.x and move it to
# its home in the ../BIN directory
# 
# "make kcarta.x" will compile kcarta.x and leave it in the current
# directory
#
# Example compiler flags and libraries are included for several
# compilers and architectures.  Initial testing should be done with
# array checking on; When things are working OK, you can leave out
# the array checking and crank up the optimization.

SHELL = /bin/sh

# SGI Fortran
# ------------
# SGI compiler options
# -u  : turn off implicit typing of variables
# -g  : generate debugging information (turns off optimization)
# -C  : do run time subscript range checking
# -w0 : inform about unused variables 
# -O3 : heavy optimization
# -64 : 64-bit objects (libraries must match)
#
#F77=/usr/bin/f77
## FLAGS = -static -C -col120
#FLAGS = -static -O2 -64 -col120 -C
#UFLAG = -u
#LIBS  = -lblas


# Linux with Absoft Fortran
# --------------------------
# Absoft compiler options
# -f	fold all names to lower case
# -N109 fold all names to upper case
# -N15  so that LAPACK lib can be included; adds trailing underscores
# -W	wide source file
# -w	suppress warning messages (absoft is more fussy than SGI or g77)
# -A	alignment warnings
# -C	check array bounds
# -O    some optimizations
# -N3   add record info to unformatted files
# -s    static allocation
#
###### this was old Absoft  F77=/usr/bin/f77 ... howard has changed .bashrc
#F77=/usr/local/absoft/bin/f77
F77=/asl/opt/absoft/absoft10.0/bin/af77
## FLAGS = -w -W -C -N3 -s -N109
FLAGS = -w -W -O -N3 -s -f -N15 -A
LIBS = -lU77 -L ../LIB/blas.ref -lblas

# Linux with g77
# ---------------
# F77=/usr/bin/g77
# FLAGS= -fcase-lower -ffixed-line-length-120
# LIBS = -L ../LIB -lblas

# SunOS FORTRAN
# -------------
# for compiling on air2.ssec.wisc.edu:Solaris 2.5.1 on a Sun Ultra-2
# -e     wide sources
# -fast  plausible optimizations
# -w     ignore warning
# -g     link w/ debugging info (don't use w/ optimization)
# 
#F77 = /opt/SUNWspro/bin/f77
#### FLAGS = -c -g -u -O0  don't seem to need -static
#### FLAGS = -e -fast -w 
#FLAGS = -e -fast -w
#LIBS = -lblas

# SUN Solaris Fortran
# -------------------
# SUN compiler options
# -g  : generate debugging information (turns off optimization)
# -O3 : heavy optimization
# -xarch: Sparc chip family
#       v8plus for UltraSPARC or UltraSPARC II  (sun4u)
#       v8 for SuperSPARC or SuperSPARC II or HyperSPARC (sun4m or Sun4d)
#       v8a for microSPARC or microSPARC II  (sun4m)
#       v7 for earlier SPARC chips  (sun4c or sun4m for Sun4-630/670/690)
# -dalign to be compatible with the Sun Performance library
# -xlic_lib=sunperf Link option for selecting the Sun Performance library
#F77=f77
#FLAGS = -O3 -dalign -xarch=v8plus
#LIBS  = -xlic_lib=sunperf

###### want to change this to /asl/www/pub/kcarta
# The following variables are used by "make dist" and "make ftp", to
# bundle up a kcarta distribution.  KTAR is a path and filename for
# the tarfile created by "make dist", and KFTP is the machine and
# path to the FTP archive.
KTAR  = /salsify/scratch4/kcarta107.tar
KTARBUG  = /salsify/scratch4/kcarta_bugfix107.tar
KFTP  = kale:/usr/local/ftp/pub/kcarta

PARAMS = kcarta.param Makefile

PARAMS_SCAT = kcarta.param scatter.param gauss.param \
              scatter_rtspec.f  Makefile

OBJSSTRING = s_writefile.o s_misc.o n_main.o n_gas_wt_spectra.o n_pth_mix.o \
             n_rad_jac_scat.o   n_output.o 
 
OBJSMAIN=kcartamain.o       kcartamisc.o         

OBSSVD = kcoeffMAIN.o       kcoeffSPL.o    kcoeffSPLJAC.o

######## this is if you want to use CKD as defined by the original people
######## with the water lineshape = lorentz
######## this code has the smooth 10cm-1 self interpolation for CKD2.3+
OBJSCON_WB = calcon.o h2oft0_wb.o h2ost0_wb.o h2ost1_wb.o \
          calconCKD2.4.o calconCKD2.4.data.o

######## this is if you want to use CKD as defined by GENLN2
######## with the water lineshape = lorentz-lorentz(25)
OBJSCON = calcon.o h2oft0.o h2ost0.o h2ost1.o \
          calconCKD2.4.o calconCKD2.4.data.o

OBJSXSEC = calxsc.o calq.o 

OBSRAD = rad_main.o rad_diff.o  rad_quad.o rad_misc.o rad_flux.o

OBSJAC = jac_main.o jac_up.o jac_down.o

OBSSCAT = scatter_main.o scatter_rtspec.o scatter_flux.o

OBJS = $(OBJSMAIN) $(OBSSVD) $(OBJSSTRING) $(OBSRAD) $(OBSJAC) \
       $(OBSSCAT) $(OBJSXSEC) $(OBJSCON)
       
install: kcarta.x
	mv kcarta.x ../BIN/kcarta.x

kcarta.x : $(OBJS)
	$(F77) -o kcarta.x $(FLAGS) $(OBJS) $(LIBS) 

kcartamain.o : kcartamain.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcartamain.f

kcartamisc.o : kcartamisc.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcartamisc.f

kcoeffMAIN.o : kcoeffMAIN.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffMAIN.f

kcoeffSPL.o : kcoeffSPL.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffSPL.f

kcoeffSPLJAC.o : kcoeffSPLJAC.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffSPLJAC.f

rad_main.o : rad_main.f $(PARAMS) 
	$(F77) -c $(FLAGS) $(UFLAG) rad_main.f

rad_misc.o : rad_misc.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) rad_misc.f

rad_diff.o : rad_diff.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) rad_diff.f

rad_quad.o : rad_quad.f $(PARAMS) gauss.param
	$(F77) -c $(FLAGS) $(UFLAG) rad_quad.f

rad_flux.o : rad_flux.f $(PARAMS) gauss.param
	$(F77) -c $(FLAGS) $(UFLAG) rad_flux.f

scatter_main.o : scatter_main.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_main.f

scatter_rtspec.o : scatter_rtspec.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_rtspec.f

scatter_flux.o : scatter_flux.f $(PARAMS_SCAT)
	$(F77) -c $(FLAGS) $(UFLAG) scatter_flux.f

jac_main.o : jac_main.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_main.f

jac_up.o : jac_up.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_up.f

jac_down.o : jac_down.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_down.f

s_writefile.o : s_writefile.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) s_writefile.f

s_misc.o : s_misc.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) s_misc.f

n_main.o : n_main.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) n_main.f

n_pth_mix.o : n_pth_mix.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) n_pth_mix.f

n_gas_wt_spectra.o : n_gas_wt_spectra.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) n_gas_wt_spectra.f

n_rad_jac_scat.o : n_rad_jac_scat.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) n_rad_jac_scat.f

n_output.o : n_output.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) n_output.f

#do the CON AND XSEC stuff
calcon.o: calcon.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) calcon.f

calconCKD2.4.o: calconCKD2.4.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) calconCKD2.4.f

calconCKD2.4.data.o: calconCKD2.4.data.f $(PARAMS)
	$(F77) -c $(FLAGS) calconCKD2.4.data.f

####this is for the correct CKD defn
h2oft0_wb.o: h2oft0_wb.f $(PARAMS)
	$(F77) -c $(FLAGS) h2oft0_wb.f

h2ost0_wb.o: h2ost0_wb.f $(PARAMS)
	$(F77) -c $(FLAGS) h2ost0_wb.f

h2ost1_wb.o: h2ost1_wb.f $(PARAMS)
	$(F77) -c $(FLAGS) h2ost1_wb.f

#####this is for the GENLN2 versions
h2oft0.o: h2oft0.f $(PARAMS)
	$(F77) -c $(FLAGS) h2oft0.f

h2ost0.o: h2ost0.f $(PARAMS)
	$(F77) -c $(FLAGS) h2ost0.f

h2ost1.o: h2ost1.f $(PARAMS)
	$(F77) -c $(FLAGS) h2ost1.f

#### this is if kCompressed Database dne for gasID 51-63
calq.o: calq.f $(PARAMS) 
	$(F77) -c $(FLAGS) calq.f 
 
calxsc.o: calxsc.f $(PARAMS) 
	$(F77) -c $(FLAGS) $(UFLAG) calxsc.f 

clean: 
	rm *.o *~ 2> /dev/null || true

dist:
	cd ../.. ;				\
	rm $(KTAR) 2> /dev/null || true;	\
	tar -cf $(KTAR)				\
		KCARTA/EVANS_COLORADO/*.f       \
		KCARTA/EVANS_COLORADO/Makefile  \
		KCARTA/BIN			\
		KCARTA/COMPARISON/COMPv1.02     \
		KCARTA/DATA			\
		KCARTA/DOC			\
		KCARTA/KLAYERS			\
		KCARTA/LIB			\
		KCARTA/MATLAB			\
		KCARTA/MATLAB103		\
		KCARTA/OUTPUT			\
		KCARTA/README.1ST		\
		KCARTA/RUN			\
		KCARTA/SCRIPTS			\
		KCARTA/SRC			\
		KCARTA/SRCv1.07			\
		KCARTA/UTILITY103		\
		KCARTA/UTILITY;			\
	chmod 664 $(KTAR)

bugfix:
	cd ../.. ;				\
	rm $(KTARBUG) 2> /dev/null || true;	\
	tar -cf $(KTARBUG)			\
		KCARTA/SRCv1.07			\
		KCARTA/KLAYERS			\
		KCARTA/MATLAB103		\
		KCARTA/SCRIPTS			\
		KCARTA/EVANS_COLORADO/*.f       \
		KCARTA/EVANS_COLORADO/Makefile  \
		KCARTA/UTILITY103;		\
	chmod 664 $(KTARBUG)

ftp: $(KTAR)
	scp $(KTAR) $(KFTP)

