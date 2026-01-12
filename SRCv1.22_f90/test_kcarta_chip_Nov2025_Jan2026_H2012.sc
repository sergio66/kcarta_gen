## this is Nov 2025 - Jan 2026 when our disks were down
## and we really only have the H2012 kCARTA databases on/asl/data/kcarta

## everything else was on xfs2 which is TWITCHING?DEAD

make clean; make -f makefile 385_H12_CO2_UMBC_default_f90
make -f Makefile_v122_Intelf90 scat
ifort_compile_onefile.sc
