#
# makefile to build kcarta.x
#
# "make" or "make scat" will compile kcarta.x and move it to
# its home in the ../BIN directory
# 
# "make kcarta.x" will compile kcarta.x and leave it in the current
# directory
#
# Example compiler flags and libraries are included for several
# compilers and architectures.  Initial testing should be done with
# array checking on; When things are working OK, you can leave out
# the array checking and crank up the optimization.

# Things to do before typing "Make"
# Remeber all .param and include files are in ../INCLUDE (kcarta.param, 
#    scatter.param, pre_defined.param, gauss.param) while for rtpdefs.f, we 
#    have a Makefile directive for the include file (see RTP section below)

# ==================================
# 1) Type "make clean" to clear up everything
# 2) Set the correct paths, parameters in ../INCLUDE/kcarta.param
#    (eg see the IEEE_BE vs IEEE_LE note below)
# 3) Set the correct HDF (HDF_INSTALL,HDFLIB) and RTP (RTPLIB) flags
#    Also set up the correct path to the RTP include file rtpdefs.f 
#    (RTPINC and RTPDEFS_FILE)
# 4) Set up the options for your compiler (SGI,Linux,Sun,...)
#    This could include setting up the correct path to the blas library ...
#    see LIBS for this
# 5) You might need to set up the RTP library package. Go to where it is
#    (eg ../../RTP/Src) Check that the Makefile options there are set up for 
#    your compiler, type "make clean" and then "make"
# 6) You might need to set up the BLAS library package. Go to where it is
#    (eg ../LIB/blas.ref) Check that the Makefile options there are set up 
#    for your compiler, type "make clean" and then "make"
#          >>>>>> next two are optional <<<<<<<<<<
# 7) Set up the SCATTER_ALL parameter below, which tells the linker which
#    scattering routines to use in the complete kcarta package (kcarta.x)
#         >>>>>> make the executable <<<<<<<
# 8) if you want basic    kcarta package, type "make" or "make basic"
#    if you want complete kcarta (jacobs+all scattering), type "make scat"

# Other things to do before getting package working
# ==================================================
# 1) You need to set up KALAYERS eventually. Go to where it is 
#    (eg ../../KLAYERS/Src/)  Check that the Makefile options 
#    there are set up for your compiler. You will also need to set the 
#    correct HDF (HDF_INSTALL,HDFLIB) and RTP (RTPLIB) flags
#    Type "make clean" and then "make"
# 2) Go to ../UTILITY  Check that the Makefile options there are set up for 
#    your compiler. Note that SGI F77 and LINUX g77 support fseeks, which
#    make the readers fast. If you are stuck with Absoft F77 then you cannot
#    use fseeks, and need to go to the ReEADwoFSEEK directory!
#      If you use SGI F77, simply type "make" and things should be fine
#      If you use LINUX g77, simply type "make" and things will crap out as
#        g77 defines fseek as a subroutine, not a function. So you have to go 
#        edit readkcBasic.f and readbinary.f, turning off the call as a 
#        function, and turning it on as a subroutine
#      If you use Absoft F77 or some other compiler that does not support 
#        fseeks, go to READwoFSEEK and type "make clean", "make" after
#        appropriately editing the Makefile there
#    Note that this only compiles the BASIC reader .. if you want the more
#    advanced reader, and/or the jacobian and flux readers, you need to
#    type "make execs"
# 3) Go to ../SCRIPTS and type
#       basic.sc BASIC/USStandardProf_NEW dodo.dat
#    If this works, looks like everything is cool as the script
#      changes a point profile to a RTP file (using MATLAB)
#      sends this profile thru klayers_rtp and outputs a layer averaged 
#      profile sends this profile, plus appropriate driver namelist file, 
#      thru kcarta sends the resulting output thru the f77 readkcBasic reader
#      (the resulting output eg dod.dat and other work files are in ../WORK)

###################### NOTE : IEEE_BE vs IEEE_LE ###########################
## There are a number of binary files that are supplied with this code.   ##
## We supply both endian versions (big endian and little endian).         ##
## However, it is the users responsibility to ensure that the code picks  ##
## up and reads the correct version of the files. The main ones are the   ##
## entire kCARTA database. In addition, there are xsecdata, CKD (water    ##
## continuum) and SOLAR radiance datafiles that need to be used correctly ##
##                                                                        ##
## To read the main kCARTA database, ensure that the correct kWaterPath   ##
## and kCompPath parameters are set in kcarta.param, before compilation   ##
## For example, we could have                                             ##
##      PARAMETER (kWaterPath = '/sergio/kcarta/h2o.ieee-be/')            ##
##      PARAMETER (kCompPath = '/sergio/kcarta/etc.ieee-be/')             ##
##      PARAMETER (kCKDPath       = '../DATA/General/CKDieee_be/')        ##
##      PARAMETER (kSolarPath     = '../DATA/General/SOLAR_be/')          ##
##                                                                        ##
## To make doubly and triply sure things have worked, do the following    ##
## To read the correct SOLAR data, CKD data and xsecdata, go to the       ##
## ../DATA/General subdirectory. If necessary, delete the supplied        ##
## symbolic links and make the new required ones eg :                     ##
##   rm SOLAR                                                             ##
##   ln -s SOLARieee_Xe SOLAR                                             ##
##   rm CKD                                                               ##
##   ln -s CKDieee_Xe CKD                                                 ##
##   rm xsecdata.dat                                                      ##
##   ln -s xsecdata.dat.XXXXX xsecdata.dat                                ##
############################################################################

############################################################################
##          set up the appropriate RTP and HDF paths and libraries        ##
##          current options are SGI and LINUX                             ##
############################################################################
# ---------
# HDF flags
# ---------
# HDF_INSTALL : set this to point to the HDF location
# HDFLIB : hdf lib, relative to install path above
#
# sgi
# HDF_INSTALL = /salsify/users/motteler/hdf
# HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz
#
# linux
HDF_INSTALL = /usr/local/hdf
HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz

# ---------
# RTP flags
# ---------
# RTPLIB : point this to the lib you created in RTP/Src_rtp
#
# sgi
# RTPLIB = -L/asl/packages/rtpV105/lib -lrtp
# RTPLIB = -L/salsify/users/motteler/rtp/lib -lrtp
#
# linux
RTPLIB = -L/asl/packages/rtpV105/lib -lrtp

# ---------------
# RTP include dir
# ---------------
# this is where the RTP package files are stored
RTPINC = -I/asl/packages/rtpV105/include
RTPDEFS_FILE = /asl/packages/rtpV105/include/rtpdefs.f

SHELL = /bin/sh

############################################################################
## set up the appropriate compiler and its options                        ##
## NOTE 1) only SGI F77 and Absoft have been tested with the autodoubling ##
##         required by DISORT                                             ##
##      2) g77 does not support Fortan structures, so it cannot be used   ##
##         when compiling the RTP related code                            ##
##      3) If you ignore the above two limitations, and can figure out if ##
##         your compiler supports autodoubling and F77 structures, then   ##
##         we have a set of options for SunOS and SUN Solaris F77         ##
############################################################################

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
#LIBS  = -lblas -lcomplib.sgimath
#FDOUBLE = -r8

# Linux with Absoft Fortran
# --------------------------
# Absoft compiler options
# -f	fold all names to lower case
# -N109 fold all names to upper case
# -W	wide source file
# -w	suppress warning messages (absoft is more fussy than SGI or g77)
# -A	alignment warnings
# -C	check array bounds
# -O    some optimizations
# -N3   add record info to unformatted files
# -s    static allocation
# -N2   force intrinsic double functions
# -N113 force double precision 
# -N114 force untyped variables as warnings
# -g    symbol and line info included, for debugging purposes
F77=f77 
UFLAG = -N114 
#FLAGS = -v -g -W -A -C -N3 -s -N109 
FLAGS = -w -W -C -O2 -A -N3 -s -f    -N15 
FLAGS_RTP = -w -W -O -C -N3 -s 
#LIBS = -lU77 -L../LIB/blas.ref -lblas -L/usr/local/absoft/lib -llapack 
LIBS = -lU77 -L/usr/local/absoft/lib -llapack -lblas 
FDOUBLE = -N2 -N113 


#Linux with PDF Fortran
#----------------------
#F77 = pgf77
#UFLAG = -Mdclchk
##FLAGS = -Mextend -Mbounds -Mupcase -O2
#FLAGS = -Mextend -Mnoupcase -O2
#FLAGS_RTP = -Mextend -Mnoupcase -O2
#LIBS = -L ../LIB/blas.ref -lblas   
#FDOUBLE = -R8

# Linux with g77
# ---------------
#F77=/usr/bin/g77
#FLAGS= -fcase-lower -ffixed-line-length-120
#LIBS = -L ../LIB -lblas

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
# -e     wide sources
# -O3 : heavy optimization
# -w     ignore warning
# -xarch: Sparc chip family
#       v8plus for UltraSPARC or UltraSPARC II  (sun4u)
#       v8 for SuperSPARC or SuperSPARC II or HyperSPARC (sun4m or Sun4d)
#       v8a for microSPARC or microSPARC II  (sun4m)
#       v7 for earlier SPARC chips  (sun4c or sun4m for Sun4-630/670/690)
# -dalign to be compatible with the Sun Performance library
# -xlic_lib=sunperf Link option for selecting the Sun Performance library
#F77=f77
#FLAGS = -O3 -dalign -xarch=v8plus -e -w
#LIBS  = -xlic_lib=sunperf

###### want to change this to /asl/www/pub/kcarta
# The following variables are used by "make dist" and "make ftp", to
# bundle up a kcarta distribution.  KTAR is a path and filename for
# the tarfile created by "make dist", and KFTP is the machine and
# path to the FTP archive.

KTAR     = /salsify/scratch4/kcarta110.tar
KTARBUG  = /salsify/scratch4/kcarta_bugfix110.tar
KTARBUGGIVEOUT  = /salsify/scratch4/kcarta_bugfix_giveout110.tar
KGIVEOUT = /salsify/scratch4/kcartagiveout110.tar
KDATA    = /salsify/scratch4/kcartaDATA.tar
KMATLAB  = /salsify/scratch4/kcartaMATLAB.tar
KFTP     = kale:/usr/local/ftp/pub/kcarta

KTAR     = /asl/ftp/pub/packages/kcarta110.tar
KTARBUG  = /asl/ftp/pub/packages/kcarta_bugfix110.tar
KTARBUGGIVEOUT  = /asl/ftp/pub/packages/kcarta_bugfix_giveout110.tar
KGIVEOUT = /asl/ftp/pub/packages/kcartagiveout110.tar
KDATA    = /asl/ftp/pub/packages/kcartaDATA.tar
KMATLAB  = /asl/ftp/pub/packages/kcartaMATLAB.tar
KFTP     = /asl/ftp/pub/packages/kcarta

KTAR     = /asl/packages/kcarta110.tar
KTARBUG  = /asl/packages/kcarta_bugfix110.tar
KTARBUGGIVEOUT  = /asl/packages/kcarta_bugfix_giveout110.tar
KGIVEOUT = /asl/packages/kcartagiveout110.tar
KDATA    = /asl/packages/kcartaDATA.tar
KMATLAB  = /asl/packages/kcartaMATLAB.tar
KCWRAP   = /asl/packages/kcwrap.tar
KFTP     = /asl/packages/kcarta

PARAMS = ../INCLUDE/kcarta.param ../INCLUDE/pre_defined.param Makefile

PARAMS_RTP = ../INCLUDE/kcarta.param ../INCLUDE/pre_defined.param Makefile

###### OBJSMAIN       is for all features : jacobians, scatter, flux
###### OBJSMAIN_SCAT  is for some of kCARTA features : jacobians, ONE scatter
######    have to edit file kcartascat.f and 
###### OBJSMAIN_BASIC is for none of kCARTA features

OBJSMAIN       = kcartamain.o       kcartamisc.o         

OBJSMAIN_SCAT  = kcartascat.o       kcartamisc.o         

OBJSMAIN_BASIC = kcartabasic.o      kcartamisc.o         

OBJSSTRING     = s_writefile.o       s_misc.o         n_main.o      \
                 n_gas_wt_spectra.o  n_pth_mix.o      n_rad_jac_scat.o  \
                 n_output.o          rtp_interface.o
 
OBSSVD         = kcoeffMAIN.o       kcoeffSPL.o    kcoeffSPLJAC.o

###########################################################################
### note that we eventually use lookup tables generated from Matlab     ###
###                   for the continuum                                 ###
### depending on the database used, set up the OBJSCON to indicate      ###
### which water continnum type to use. We generated our latest version  ###
### of our kCompressed database assuming the CKD lineshape, so this     ###
### pretty much meas use OBJSCON_WB!!!                                  ###

####### since we are using new kCARTA database, use this version
####### to use CKD as originally defined (wb === without basement term!!!)
####### with the water lineshape = lorentz-lorentz(25)
####### this code has the smooth 10cm-1 self interpolation for CKD2.3+
OBJSCON_WB = calcon.o h2oft0_wb.o h2ost0_wb.o h2ost1_wb.o \
          calconCKD2.4.o calconCKD2.4.data.o

####### if using the 1996-1999 kCARTA database, use this version
####### to use CKD as used by older kCARTA versions (defined by GENLN2)
####### with the water lineshape = lorentz
OBJSCON_OLD_GENLN2 = calcon.o h2oft0.o h2ost0.o h2ost1.o \
          calconCKD2.4.o calconCKD2.4.data.o

## choose one of the above two
OBJSCON = $(OBJSCON_WB)
#######################################################################

OBJSXSEC     = calxsc.o calq.o 

OBSRAD       = rad_main.o rad_diff.o  rad_quad.o rad_misc.o rad_flux.o

OBSRAD_BASIC = rad_main.o rad_diff.o  rad_quad.o rad_misc.o

OBSJAC       = jac_main.o jac_up.o jac_down.o

#######################################################################
###   this is to set up your favorite scattering routine package    ###
###   RTSPEC : fast but does not include solar beam                 ###
###   DISORT : slow, computes the whole kit and caboodle            ###
###   SERGIO : a trial at using eddington to compute Jacobians      ###

PARAMS_SCAT = ../INCLUDE/kcarta.param ../INCLUDE/scatter.param \
	      ../INCLUDE/pre_defined.param ../INCLUDE/gauss.param Makefile

OBSSCAT_ALL    = scatter_all_misc.o

OBSSCAT_TWOSTR = scatter_twostream_main.o scatter_twostream_code_guts.o \
		 scatter_twostream_code.o scatter_twostream_flux.o 
 
OBSSCAT_RTSPEC = scatter_rtspec_main.o scatter_rtspec_code.o \
                 scatter_rtspec_flux.o

OBSSCAT_DISORT = scatter_disort_main.o scatter_disort_code.o \
                 scatter_disort_flux.o scatter_disort_misc.o

OBSSCAT_DIS = $(OBSSCAT_DISORT) $(OBSSCAT_ALL) 
OBSSCAT_RTS = $(OBSSCAT_RTSPEC) $(OBSSCAT_ALL) 
OBSSCAT_TWO = $(OBSSCAT_TWOSTR) $(OBSSCAT_ALL) 
OBSSCAT3 = $(OBSSCAT_DISORT) $(OBSSCAT_RTSPEC) $(OBSSCAT_TWOSTR) \
           $(OBSSCAT_ALL) 

############################################################################
## set up the appropriate scattering pakages                              ##
## NOTE 1) kcarta.x can include one or more of DISORT (slow, with solar), ##
##         RTSPEC (fast, no solar), TWOSTREAM (fast, with solar, less accurate)
##      2) skcarta.x can only inlude one of the above three               ##
##      3) bkcarta only has clear sky radiances, fluxes, jacobians        ##
############################################################################

### for COMPLETE kcarta package (kcarta.x) choose one of the 5 options
### SCATTER_ALL = OBSSCAT3
### default is to turn on ALL scattering packages (DISORT,RTSPEC,TWOSTREAM)
### in addition, you get clear sky radiance, jacobian, flux code compiled in
SCATTER_ALL = $(OBSSCAT3)   
OBJS = $(OBJSMAIN) $(OBSSVD) $(OBSRAD) $(OBSJAC) \
       $(SCATTER_ALL) $(OBJSXSEC) $(OBJSCON) $(OBJSSTRING) 

### for BASIC clear sky package (bkcarta.x), here we are!!!!!!!!
OBJS_BASIC = $(OBJSMAIN_BASIC) $(OBSSVD) $(OBSRAD_BASIC) \
       $(OBJSXSEC) $(OBJSCON) $(OBSSCAT_ALL) $(OBSJAC) $(OBJSSTRING) 

########################################################################
################## MAIN PART OF MAKEFILE ###############################
########################################################################
# executables live in ../BIN
EXECS = basic scat
EXECS = basic
all : $(EXECS)

kcartaparam: 
	rm ../INCLUDE/kcarta.param; \
      cd ../INCLUDE; \
      ln -s kcarta_v110.param kcarta.param; \
      cd ../SRCv1.10
     

############## the two executables : 
##### (1) complete kcarta (memory hog) with jacobians, all scattering codes
scat: kcarta.x
	mv kcarta.x ../BIN/kcarta.x

kcarta.x : $(OBJS)
	$(F77) -o kcarta.x $(FLAGS) $(OBJS) $(LIBS) $(RTPLIB) $(HDFLIB)

kcartamain.o : kcartamain.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcartamain.f

##### (2) basic kcarta (easy on memory) with no jacobs, fluxes or scattering
basic: bkcarta.x
	mv bkcarta.x ../BIN/bkcarta.x

bkcarta.x : $(OBJS_BASIC)
	$(F77) -o bkcarta.x $(FLAGS) $(OBJS_BASIC) $(LIBS) $(RTPLIB) $(HDFLIB)

kcartabasic.o : kcartabasic.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcartabasic.f

################

kcartamisc.o : kcartamisc.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcartamisc.f

### kCompressed Database uncompression routines
kcoeffMAIN.o : kcoeffMAIN.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffMAIN.f

kcoeffSPL.o : kcoeffSPL.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffSPL.f

kcoeffSPLJAC.o : kcoeffSPLJAC.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) kcoeffSPLJAC.f

### Non scattering radtransfer  code (NON SCATTERING ONLY)
rad_main.o : rad_main.f $(PARAMS) 
	$(F77) -c $(FLAGS) $(UFLAG) rad_main.f

rad_misc.o : rad_misc.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) rad_misc.f

rad_diff.o : rad_diff.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) rad_diff.f

rad_quad.o : rad_quad.f $(PARAMS) ../INCLUDE/gauss.param
	$(F77) -c $(FLAGS) $(UFLAG) rad_quad.f

rad_flux.o : rad_flux.f $(PARAMS) ../INCLUDE/gauss.param
	$(F77) -c $(FLAGS) $(UFLAG) rad_flux.f

#### many of the routines in here are common to all scattering codes 
scatter_all_misc.o : scatter_all_misc.f  $(PARAMS_SCAT)  
	$(F77) -c $(FLAGS) $(UFLAG) scatter_all_misc.f 

### this is Sergio Machado's twostream (+solar) radtransfer code 
scatter_twostream_main.o : scatter_twostream_main.f  $(PARAMS_SCAT)  
	$(F77) -c $(FLAGS) $(UFLAG) scatter_twostream_main.f 
 
scatter_twostream_code.o : scatter_twostream_code.f  $(PARAMS_SCAT)  
	$(F77) -c $(FLAGS) $(UFLAG) scatter_twostream_code.f 
 
scatter_twostream_code_guts.o : scatter_twostream_code_guts.f  $(PARAMS_SCAT)  
	$(F77) -c $(FLAGS) $(UFLAG) $(FDOUBLE) scatter_twostream_code_guts.f 
#	$(F77) -c $(FLAGS) $(UFLAG) scatter_twostream_code_guts.f 
 
scatter_twostream_flux.o : scatter_twostream_flux.f  $(PARAMS_SCAT)  
	$(F77) -c $(FLAGS) $(UFLAG) scatter_twostream_flux.f 
 
#### this is Frank Evans RTSPEC
scatter_rtspec_main.o : scatter_rtspec_main.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_rtspec_main.f

scatter_rtspec_code.o : scatter_rtspec_code.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_rtspec_code.f

scatter_rtspec_flux.o : scatter_rtspec_flux.f $(PARAMS_SCAT)
	$(F77) -c $(FLAGS) $(UFLAG) scatter_rtspec_flux.f

#### this is Kurt Stamnes DISORT
scatter_disort_main.o : scatter_disort_main.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_disort_main.f

scatter_disort_flux.o : scatter_disort_flux.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(UFLAG) scatter_disort_flux.f

scatter_disort_code.o : scatter_disort_code.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(FDOUBLE) $(UFLAG) scatter_disort_code.f

scatter_disort_misc.o : scatter_disort_misc.f  $(PARAMS_SCAT) 
	$(F77) -c $(FLAGS) $(FDOUBLE) $(UFLAG) scatter_disort_misc.f

### Jacobian code (NON SCATTERING ONLY)
jac_main.o : jac_main.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_main.f

jac_up.o : jac_up.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_up.f

jac_down.o : jac_down.f $(PARAMS)
	$(F77) -c $(FLAGS) $(UFLAG) jac_down.f

### misc and I/O routines
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

rtp_interface.o : rtp_interface.f $(PARAMS_RTP) $(RTPDEFS_FILE)
	$(F77) -c $(RTPINC) $(FLAGS) $(UFLAG) rtp_interface.f

### CON AND XSEC stuff
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

#### this is if kCompressed Database DOES NOT EXIST for gasID 51-63
calq.o: calq.f $(PARAMS) 
	$(F77) -c $(FLAGS) calq.f 
 
calxsc.o: calxsc.f $(PARAMS) 
	$(F77) -c $(FLAGS) $(UFLAG) calxsc.f 

clean: 
	rm *.o *~ 2> /dev/null || true

################### DISTRIBUTIONS ########################################
#### 1) to make the distributions, go to /home/sergio and look for kcartaV1xx
#### 2) since all the KCARTA directories are under kcartaV1xx, rename this to
####    kcartaV110 eg            mv kcartaV111 kcartaV110
#### 3) rm old sybolic link      rm KCARTA
#### 4) make new symbolic link   ln -s  kcartaV110 KCARTA
#### 5) cd to this dir              cd KCARTA 
#### 6) rm old symbolic link SRC    rm SRC 
#### 7) make new symbolic link SRC  ln -s SRCv1.10 SRC 
#### 8) cd to this dir              cd SRC
#### 9) make whatever distribution you want! eg make dist

#### (A) : this has ALL the source code 
####       clear sky rad transfer, clear sky jacobians, scatter, nlte and
####       the nlte profiles
dist:
	cd ../.. ;				\
	rm $(KTAR) 2> /dev/null || true;	\
	tar -cvf $(KTAR)			\
		kcartaV110/LIB                  \
		kcartaV110/INCLUDE 		\
		kcartaV110/BIN			\
		kcartaV110/SCATTERCODE          \
		kcartaV110/TEST                 \
		kcartaV110/DOC			\
		kcartaV110/MATLAB		\
		kcartaV110/README.1ST		\
		kcartaV110/SCRIPTS		\
		kcartaV110/SRC			\
		kcartaV110/SRCv1.10		\
		kcartaV110/UTILITY;		\
	chmod 664 $(KTAR)

bugfix:
	cd ../.. ;				\
	rm $(KTARBUG) 2> /dev/null || true;	\
	tar -cf $(KTARBUG)			\
		kcartaV110/README.1ST		\
		kcartaV110/INCLUDE		\
		kcartaV110/DOC/kcarta1.10*	\
		kcartaV110/SRCv1.10/Makefile	\
		kcartaV110/SRCv1.10/*.f;	\
	chmod 664 $(KTARBUG)

##############################
#### (B) this is the MATLAB directory, for reading in binary kCARTA output
matlab:
	cd ../.. ;				\
	rm $(KMATLAB) 2> /dev/null || true;     \
	tar -cf $(KMATLAB)			\
		kcartaV110/MATLAB;		\
	chmod 664 $(KMATLAB)
		        
data:
	cd ../.. ;				\
	rm $(KDATA) 2> /dev/null || true;	\
	tar -cf $(KDATA)			\
		KCARTADATA;			\
	chmod 664 $(KDATA)

##############################
### (C) this is the very useful wrapper that Howard Motteler wrote
###     the main version would be in /home/motteler/radtrans/kcwrap/kcwrap
###     some copies would exist in the KCARTA/SCRIPTS directory
kcwrap:
	cd ../.. ;				\
	rm $(KCWRAP) 2> /dev/null || true;	\
	tar -cf $(KCWRAP)			\
		/asl/matlab/rtptools		\
		/asl/matlab/h4tools		\
		/asl/data/airs/srf/srftablesV10.hdf \
		kcartaV110/SCRIPTS;		\
	chmod 664 $(KCWRAP)

ftp: $(KTAR)
	scp $(KTAR) $(KFTP)

##############################
#### (D) these are the basic files we give out : 
####     clear sky radtransfer and jacobians, as well as nlte source code

#these are the files for the basic kCARTA package
#kc1 = c*.f h*.f j*.f kcarta*.f kcoeff*.f n*.f \
#      r*.f s_*.f  INCLUDE/*.param Makefile
giveout:
	cd ../.. ;				\
	rm $(KGIVEOUT) 2> /dev/null || true;	\
	tar -cf $(KGIVEOUT)			\
                kcartaV110/LIB                  \
                kcartaV110/INCLUDE              \
                kcartaV110/BIN                  \
                kcartaV110/TEST                 \
                kcartaV110/DOC                  \
                kcartaV110/MATLAB               \
                kcartaV110/README.1ST           \
                kcartaV110/SCRIPTS              \
                kcartaV110/INCLUDE/*.param      \
                kcartaV110/SRCv1.10/*.param     \
                kcartaV110/SRCv1.10/c*.f        \
                kcartaV110/SRCv1.10/h*.f        \
                kcartaV110/SRCv1.10/j*.f        \
                kcartaV110/SRCv1.10/kcarta*.f   \
                kcartaV110/SRCv1.10/kcoeff*.f   \
                kcartaV110/SRCv1.10/n*.f        \
                kcartaV110/SRCv1.10/r*.f        \
                kcartaV110/SRCv1.10/s_*.f       \
                kcartaV110/SRCv1.10/scatter_all_misc.f       \
                kcartaV110/SRCv1.10/Makefile    \
		kcartaV110/UTILITY;		\
	chmod 664 $(KGIVEOUT)

giveout_bugfix:
	cd ../.. ;				\
	rm $(KTARBUGGIVEOUT) 2> /dev/null || true;	\
	tar -cf $(KTARBUGGIVEOUT)			\
		kcartaV110/README.1ST		\
		kcartaV110/INCLUDE		\
		kcartaV110/DOC/kcarta1.10*	\
		kcartaV110/SRCv1.10/Makefile    \
                kcartaV110/SRCv1.10/c*.f        \
                kcartaV110/SRCv1.10/h*.f        \
                kcartaV110/SRCv1.10/j*.f        \
                kcartaV110/SRCv1.10/k*.f   \
                kcartaV110/SRCv1.10/n*.f        \
                kcartaV110/SRCv1.10/r*.f        \
                kcartaV110/SRCv1.10/s_*.f       \
                kcartaV110/SRCv1.10/scatter_all_misc.f       \
                kcartaV110/SRCv1.10/Makefile;    \
	chmod 664 $(KTARBUGGIVEOUT)

################### DISTRIBUTIONS ########################################
