## run with    run90.sc >& ugh 

## see /home/motteler/shome/rtp/Makefile, copied to /home/sergio/git/rtp/rtpV221
# get HDF libs and Intel compilers
module load HDF/4.2.14-GCCcore-8.3.0
module load intel-compilers/2021.4.0

# old stuff, libraries probably broken
#ifort -c makeRTPfile.f90 -I/home/sergio/git/rtp/rtpV201/include -extend-source 132 -check bounds -O2

echo "now ifort the code and link"
#ifort -c makeRTPfile.f90 -I/home/motteler/shome/rtp/include     -extend-source 132 -check bounds -O2
#ifort -o makeRTPfile.x -extend-source 132 -check bounds -O2 makeRTPfile.o \
#   -L/home/motteler/shome/rtp/lib -lrtp \
#   -L/usr/ebuild/software/HDF/4.2.14-GCCcore-8.3.0/lib -ldf
ifort -c makeRTPfile.f90 -I/home/sergio/git/rtp/rtpV221/include -extend-source 132 -check bounds -O2
ifort -o makeRTPfile.x -extend-source 132 -check bounds -O2 makeRTPfile.o \
   -L/home/sergio/git/rtp/rtpV221/lib -lrtp \
   -L/usr/ebuild/software/HDF/4.2.14-GCCcore-8.3.0/lib -ldf

mv makeRTPfile.x ../BIN/makeRTPfile.x

#########################

echo "now run the code"
rm junk.ip.rtp junk.op.rtp
../BIN/makeRTPfile.x fin=day_py4cats.atmX fout=junk.ip.rtp

echo "now look at junk.ip.rtp"
/home/sergio/git/rtp/rtpV201/utils/rtpdump -h -p junk.ip.rtp

#/asl/packages/klayersV205/BinV201/klayers_airs fin=junk.ip.rtp fout=junk.op.rtp
/home/sergio/SARTA_CLOUDY_RTP_KLAYERS_NLEVELS/klayersV205_140levs/BinV201/klayers_airs_wetwater140levs fin=junk.ip.rtp fout=junk.op.rtp
