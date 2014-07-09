# copied from /asl/packages/rtpV105/src
# 
# FSEEK library makefile

# set FSEEKHOME to point to the local FSEEK installation
FSEEK_LOCAL_HOME = ..

# ------------------
# C compiler options
# -------------------

CFLAGS = -O
CC = cc

# ------------
# make targets
# -------------

all: libfseek_local.a install

# -------------
# C/Fortran API
# -------------

FORSRC = fseek_ftell.c 
FOROBS = fseek_ftell.o 

$(FOROBS): $(FORSRC) 
	   $(CC) $(CFLAGS) -c $(FORSRC) 


# -----------
# FSEEK library
# -----------

libfseek_local.a: $(FOROBS) 
	  ar -rcs libfseek_local.a $(FOROBS)

# -------
# install
# -------

install:libfseek_local.a
	cp -a libfseek_local.a $(FSEEK_LOCAL_HOME)/lib

# -------
# cleanup
# -------

clean:
	rm *.o librtp.a 2> /dev/null || true







