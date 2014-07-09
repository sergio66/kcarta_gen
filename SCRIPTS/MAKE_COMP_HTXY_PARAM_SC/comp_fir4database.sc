:
# this makes comp.param and xsec.param, which contains lists of the 
# available compressed data files. Make sure you look at INCLUDE/kcarta.param
# and set up the correct paths to kWaterPath,kCompPath and 
# kXsecParamFile,kCompParamFile

### WARNING (1) : 
### look at kcarta.param and set correct paths to kWaterPath,kCompPath
 
shome=`pwd`

#first do the water
cd /asl/data/kcarta/H2012.ieee-le/FIR80/h2o.ieee-le  ##kWaterPath
ls -1  >& $shome/waterdatabase_fir80_150
cd $shome

#then do CO2
#cd $shome
#cd /spinach/s6/sergio/RUN8_NIRDATABASE/FIR80/fbin/etc.ieee-le  ##kCO2Path
#ls -1  >& $shome/co2database_fir80_150
#cd $shome

#then do the rest of the gases
cd $shome
cd /asl/data/kcarta/H2012.ieee-le/FIR80/etc.ieee-le  ##kCompPath
ls -1  >& $shome/othersdatabase_fir80_150
cd $shome

## but now have to get rid of gas2 within othersdatabase1
## sed '/_g2.dat/ d' othersdatabase1 > othersdatabase  

#####now put the two files together and process them!!!!!!!!
cat waterdatabase_fir80_150 othersdatabase_fir80_150 > compdatabase_fir80_150

if [ -r compdatabase ]
then
  rm compdatabase
fi
ln -s compdatabase_fir80_150 compdatabase

######### run code compdatabase.x to produce comp.param
cp ../../BIN/compdatabase.x .
if [ -r comp.param ]
then
  rm comp.param
fi
compdatabase.x > dope

mv comp.param PARAM_TEMP/comp_fir80_150.param
echo " ---->>>> results are stored in PARAM_TEMP/comp_fir80_150.param <<<<----"

rm compdatabase_fir80_150 waterdatabase_fir80_150 comp1.param
rm othersdatabase_fir80_150 compdatabase.x dope compdatabase




