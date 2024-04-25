########################################################################
echo "linking .f files to the .f.g77 code ...."

/bin/rm readbinary.f
/bin/rm readbinaryjacob.f
/bin/rm readflux.f
/bin/rm readjacob.f
/bin/rm readkcarta.f
/bin/rm readkcBasic.f

ln -s readbinary.f.g77       readbinary.f
ln -s readbinaryjacob.f.g77  readbinaryjacob.f
ln -s readflux.f.g77         readflux.f
ln -s readjacob.f.g77        readjacob.f
ln -s readkcarta.f.g77       readkcarta.f
ln -s readkcBasic.f.g77      readkcBasic.f

########################################################################
echo "linking .f90 files to the .f.g77 code ...."

/bin/rm readbinary.f90
/bin/rm readbinaryjacob.f90
/bin/rm readflux.f90
/bin/rm readjacob.f90
/bin/rm readkcarta.f90
/bin/rm readkcBasic.f90

ln -s readbinary.f.g77       readbinary.f90
ln -s readbinaryjacob.f.g77  readbinaryjacob.f90
ln -s readflux.f.g77         readflux.f90
ln -s readjacob.f.g77        readjacob.f90
ln -s readkcarta.f.g77       readkcarta.f90
ln -s readkcBasic.f.g77      readkcBasic.f90

