## this is Nov 2025 - Jan 2026 when our disks were down
## and we really only have the H2012 kCARTA databases on/asl/data/kcarta

## everything else was on xfs2 which is TWITCHING?DEAD

make clean; make -f makefile 385_H12_default
make -f Makefile_v118_Intel scat
ifort_compile_onefile.sc
