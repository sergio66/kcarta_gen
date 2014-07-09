#
# makefile to build utility programs
# see ../SRC/Makefile for more info on compiler options
# note : if you want to use this, make sure your compiler supports fseek
#        without reservations. According to the manuals, Absoft supports
#        fseek/ftell, but I could never get things to work correctly!
#        so I use g77 or PDF Fortran for the readers.
# OTOH, the RTP files have structures in them, which g77 does not support, so
#       makeRTPfile.f has to be compiled with Absoft

### for Linux with g77 and other (eg Absoft, PDF) compilers, the RTP stuff is 
### compiled with Absoft or PDF (hence we have RTPF77,RTPFLAGS,RTPLIBS etc),
### while readers are compiled with g77 (hence we have F77,FLAGS,LIBS etc)

### for SGI, you can compile everything with f77 if you want; if so, you do 
### need to change the calling for fseek() from subroutine to function

### worked fine with RTPv103 ...

SHELL = /bin/sh

# RTP flags
# ---------
# RTPLIB : point this to the lib you created in RTP
#
# sgi
# RTPLIB = -L/asl/packages/rtpV105/lib -lrtp
#
# linux
RTPLIB = -L/asl/packages/rtpV105/lib -lrtp

# ---------------
# RTP include dir
# ---------------
# this is where the RTP package files are stored
RTPINC = -I/asl/packages/rtpV105/include
RTPDEFS_FILE = /asl/packages/rtpV105/include/rtpdefs.f

# ---------
# HDF flags
# ---------
# HDF_INSTALL : set this to point to the HDF location
# HDFLIB : hdf library, relative to install path above
# 
# linux
HDF_INSTALL = /home/motteler/asl/hdf/HDF4.1r3/NewHDF
HDF_INSTALL = /usr/local/hdf
HDF_INSTALL = /asl/opt
HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz
#
# sgi
#HDF_INSTALL = /salsify/users/motteler/hdftest
#HDFLIB   = -L$(HDF_INSTALL)/lib -lmfhdf -ldf -ljpeg -lz

## ------------------------ COMPILER TYPES --------------------------

## ------------------------ SGI --------------------------
# SGI Fortran ... supports fseek and structures, so can be used to compile
# all files
# SGI compiler options
# -u  : turn off implicit typing of variables
# -g  : generate debugging information (turns off optimization)
# -C  : do run time subscript range checking
# -w0 : inform about unused variables 
# -O3 : heavy optimization
# -64 : 64-bit objects (libraries must match)
#F77   = /usr/bin/f77
#FLAGS = -static -C -u -col120
#LIBS = 
#RTPF77   = /usr/bin/f77
#RTPFLAGS = -static -C -u -col120
#RTPLIBS = 

## ------------------------ Linux g77 --------------------------
# Linux with g77 .. supports fseek but NO structiures, and so CANNOT be
# used to compile makeRTPfile.x
#F77   = /usr/bin/g77
#FLAGS = -fcase-lower -ffixed-line-length-120 
#F77   = /asl/opt/bin/g95
#FLAGS = -ffixed-line-length-132 -ftrace=full 
RTPF77 = pgf77
#LIBS = 

## ------------------------ Linux Absoft f77 --------------------------
# Linux with Absoft Fortran ... does not support fseek the way we want!
# and so it SHOULD NOT be used to compile the F77 readers. However it 
# supports structures, and so it can be used to compile makeRTPfile.x
# Absoft compiler options
# -f    fold all names to lower case
# -N109 fold all names to upper case
# -W    wide source file
# -w    suppress warning messages (absoft is more fussy than SGI or g77)
# -A    alignment warnings
# -C    check array bounds
# -O    some optimizations
# -N3   add record info to unformatted files
# -s    static allocation
# -N2   force intrinsic double functions
# -N113 force double precision 
# -N114 force untyped variables as warnings
#F77   = /asl/opt/absoft/absoft10.0/bin/af77
#FLAGS = -w -W -s -C -N109 -N114 -N3
#RTPF77   = f77
#RTPF77   = /asl/opt/absoft/absoft10.0/bin/af77
#RTPFLAGS = -w -W -s -C -N109 -N114 -N3
#RTPLIBS = -lU77

## ------------------------ Linux Portland Group f77 -------------------
#Linux with PDF Fortran supports fseek and structures, so can be used to 
#compile all files
F77 = pgf77
FLAGS = -Mextend -Mbounds -Mnoupcase -O2
LIBS = 
RTPF77 = pgf77
RTPFLAGS = -Mextend -Mbounds -Mnoupcase -O2
RTPLIBS = 

## ------------------------ SUN f77 --------------------------
# SunOS fortran supports fseek and structures, and so can be used to 
# compile all files
# -------------
# for compiling on air2.ssec.wisc.edu:Solaris 2.5.1 on a Sun Ultra-2
# -e     wide sources
# -fast  plausible optimizations
# -w     ignore warning
# -g     link w/ debugging info (don't use w/ optimization)
### FLAGS = -c -g -u -O0 -static
### FLAGS = -e -fast -w -static
#F77 = /opt/SUNWspro/bin/f77
#FLAGS = -e -fast -w
#LIBS = 
#RTPF77 = /opt/SUNWspro/bin/f77
#RTPFLAGS = -e -fast -w
#RTPLIBS = 

# EXECS defines the set of all programs to be compiled
# BASIC defines the set of core programs to be compiled
BASIC = readkcarta.x readkcBasic.x compdatabase.x makeRTPfile.x
EXECS = readkcarta.x readjacob.x compdatabase.x readflux.x readkcBasic.x \
        makeRTPfile.x

PARAMS = convolve.param ../INCLUDE/kcarta.param

# executables live in ../BIN
#

basic: $(BASIC) Makefile
	mv $(BASIC) ../BIN

all: $(EXECS) Makefile
	mv $(EXECS) ../BIN

# executables
#

##RTP creator
makeRTPfile.x: makeRTPfile.o 
	$(RTPF77) -o makeRTPfile.x $(RTPFLAGS) makeRTPfile.o \
	$(RTPLIB) $(RTPLIBS) $(HDFLIB) 

readkcarta.x : readkcarta.o readbinary.o
	$(F77) -o readkcarta.x $(FLAGS) readkcarta.o readbinary.o  $(RTPLIBS)

readkcBasic.x : readkcBasic.o readbinary.o
	$(F77) -o readkcBasic.x $(FLAGS) readkcBasic.o readbinary.o $(RTPLIBS)

readjacob.x : readjacob.o readbinaryjacob.o
	$(F77) -o readjacob.x $(FLAGS) readjacob.o readbinaryjacob.o $(RTPLIBS)

readflux.x : readflux.o readbinaryjacob.o
	$(F77) -o readflux.x $(FLAGS) readflux.o readbinaryjacob.o $(RTPLIBS)

compdatabase.x: compdatabase.f
	$(F77) -o compdatabase.x $(FLAGS) compdatabase.f

makeinp.x: makeinp.f
	$(F77) -o makeinp.x $(FLAGS) makeinp.f

# object files
makeRTPfile.o : makeRTPfile.f $(RTPDEFS_FILE) 
	$(RTPF77) -c makeRTPfile.f $(RTPINC) $(RTPFLAGS) makeRTPfile.f

readbinary.o : readbinary.f $(PARAMS)
	$(F77) -c $(FLAGS) readbinary.f

readbinaryjacob.o : readbinaryjacob.f $(PARAMS)  
	$(F77) -c $(FLAGS) readbinaryjacob.f

responsefcns11.o : responsefcns11.f $(PARAMS)
	$(F77) -c $(FLAGS) responsefcns11.f

readkcarta.o : readkcarta.f $(PARAMS)
	$(F77) -c $(FLAGS) readkcarta.f

readkcBasic.o : readkcBasic.f $(PARAMS)
	$(F77) -c $(FLAGS) readkcBasic.f

readjacob.o : readjacob.f $(PARAMS)
	$(F77) -c $(FLAGS) readjacob.f

readflux.o : readflux.f $(PARAMS)
	$(F77) -c $(FLAGS) readflux.f

#not needed for SRCv1.02+
#makeprofile.x: makeprofile.f
#	$(F77) -o makeprofile.x $(FLAGS) makeprofile.f
#readjacob2.x : readjacob2.o readbinaryjacob.o
#	$(F77) -o readjacob2.x $(FLAGS) readjacob2.o readbinaryjacob.o
#sturner_offset.o : sturner_offset.f 
#	$(F77) -c $(FLAGS) sturner_offset.f
#offset.x : sturner_offset.o
#	$(F77) -o offset.x $(FLAGS) sturner_offset.o 

clean : 
	rm *.o *~ 2> /dev/null || true

