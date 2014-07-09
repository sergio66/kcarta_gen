/bin/rm readbinary.f
/bin/rm readbinaryjacob.f
/bin/rm readflux.f
/bin/rm readjacob.f
/bin/rm readkcarta.f
/bin/rm readkcBasic.f

echo "linking .f files to the .f.others code ...."

ln -s readbinary.f.others       readbinary.f
ln -s readbinaryjacob.f.others  readbinaryjacob.f
ln -s readflux.f.others         readflux.f
ln -s readjacob.f.others        readjacob.f
ln -s readkcarta.f.others       readkcarta.f
ln -s readkcBasic.f.others      readkcBasic.f

