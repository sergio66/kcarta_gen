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
