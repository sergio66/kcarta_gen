#
# makefile to build utility programs
#
# see ../SRC/Makefile for more info on compiler options

SHELL = /bin/sh

# SGI Fortran
# ------------
#F77=/usr/bin/f77
#FLAGS= -static -C -u -col120
#LIBS = 

# Linux with Absoft Fortran
# --------------------------
F77=/usr/bin/f77
F77 = f77
FLAGS = -f -w -W -C -N3
FLAGS = -w -W -O -N3 -s -N109
LIBS = -lU77

# Linux with g77
# ---------------
#F77   = /usr/bin/g77
#FLAGS = -ffixed-line-length-120
#LIBS = 

# SunOS fortran
# -------------
# for compiling on air2.ssec.wisc.edu:Solaris 2.5.1 on a Sun Ultra-2
# -e     wide sources
# -fast  plausible optimizations
# -w     ignore warning
# -g     link w/ debugging info (don't use w/ optimization)
#F77 = /opt/SUNWspro/bin/f77
### FLAGS = -c -g -u -O0 -static
### FLAGS = -e -fast -w -static
#FLAGS = -e -fast -w
#LIBS = 

# EXECS defines the set of programs to be compiled
EXECS = readkcarta.x readjacob.x readkcBasic.x
BASIC = readkcBasic.x

# executables live in ../BIN
#
basic: $(BASIC)
	mv $(BASIC) ../../BIN

all: $(EXECS)
	mv $(EXECS) ../../BIN

# executables
#
readkcarta.x : readkcarta2.o readbinary2.o
	$(F77) $(FLAGS) -o readkcarta.x readkcarta2.o readbinary2.o $(LIBS)

readkcBasic.x : readkcBasic2.o readbinary2.o
	$(F77) $(FLAGS) -o readkcBasic.x readkcBasic2.o readbinary2.o $(LIBS)

readjacob.x : readjacob2.o readbinaryjacob2.o
	$(F77) $(FLAGS) -o readjacob.x readjacob2.o readbinaryjacob2.o $(LIBS)

# object files
#
readbinary2.o : readbinary2.f convolve.param
	$(F77) -c $(FLAGS) readbinary2.f

readbinaryjacob2.o : readbinaryjacob2.f  convolve.param
	$(F77) -c $(FLAGS) readbinaryjacob2.f

readkcarta2.o : readkcarta2.f convolve.param
	$(F77) -c $(FLAGS) readkcarta2.f

readkcBasic2.o : readkcBasic2.f convolve.param
	$(F77) -c $(FLAGS) readkcBasic2.f

readjacob2.o : readjacob2.f convolve.param
	$(F77) -c $(FLAGS) readjacob2.f

clean : 
	rm *.o *~ 2> /dev/null || true

