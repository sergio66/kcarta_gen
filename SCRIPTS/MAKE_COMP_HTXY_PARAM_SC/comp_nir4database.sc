:
# this makes comp.param and xsec.param, which contains lists of the 
# available compressed data files. Make sure you look at INCLUDE/kcarta.param
# and set up the correct paths to kWaterPath,kCompPath and 
# kXsecParamFile,kCompParamFile

### WARNING (1) : 
### look at kcarta.param and set correct paths to kWaterPath,kCompPath
 
shome=`pwd`

#first do the water
cd /spinach/s6/sergio/RUN8_NIRDATABASE/NIR8250_12250/fbin/h2o.ieee-le  ##kWaterPath
ls -1  >& $shome/waterdatabase_nir8250_12250
cd $shome

#then do CO2
#cd $shome
#cd /spinach/s6/sergio/RUN8_NIRDATABASE/NIR8250_12250/fbin/etc.ieee-le  ##kCO2Path
#ls -1  >& $shome/co2database_nir8250_12250
#cd $shome

#then do the rest of the gases
cd $shome
cd /spinach/s6/sergio/RUN8_NIRDATABASE/NIR8250_12250/fbin/etc.ieee-le  ##kCompPath
ls -1  >& $shome/othersdatabase_nir8250_12250
cd $shome

## but now have to get rid of gas2 within othersdatabase1
## sed '/_g2.dat/ d' othersdatabase1 > othersdatabase  

#####now put the two files together and process them!!!!!!!!
cat waterdatabase_nir8250_12250 othersdatabase_nir8250_12250 > compdatabase_nir8250_12250

if [ -r compdatabase ]
then
  rm compdatabase
fi
ln -s compdatabase_nir8250_12250 compdatabase

######### run code compdatabase.x to produce comp.param
cp ../BIN/compdatabase.x .
if [ -r comp.param ]
then
  rm comp.param
fi
compdatabase.x > dope

mv comp.param comp_nir8250_12250.param
echo " ---->>>> results are stored in comp_nir8250_12250.param <<<<----"

rm compdatabase_nir8250_12250 waterdatabase_nir8250_12250 comp1.param
rm othersdatabase_nir8250_12250 compdatabase.x dope compdatabase




