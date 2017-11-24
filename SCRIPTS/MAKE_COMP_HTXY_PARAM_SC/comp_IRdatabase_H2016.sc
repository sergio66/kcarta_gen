:
# this makes comp.param and xsec.param, which contains lists of the 
# available compressed data files. Make sure you look at INCLUDE/kcarta.param
# and set up the correct paths to kWaterPath,kCompPath and 
# kXsecParamFile,kCompParamFile

### WARNING (1) : 
### look at kcarta.param and set correct paths to kWaterPath,kCompPath
 
shome=`pwd`

#first do the water
cd /asl/data/kcarta/H2016.ieee-le/IR605/hdo.ieee-le/                     ##kWaterPath
ls -1 *g1.dat >& $shome/waterdatabase_ir605_2830
cd $shome

#then do CO2   COMMENT THIS OUT AS NEEDED
cd $shome
cd /asl/data/kcarta_sergio/UMBC_CO2_H1998.ieee-le/CO2ppmv400.ieee-le/    ##CO2Path
ls -1  >& $shome/co2database_ir605_2830
cd $shome

#then do the rest of the gases
cd $shome
cd /asl/data/kcarta/H2016.ieee-le/IR605/etc.ieee-le/                     ##kCompPath
ls -1  >& $shome/othersdatabase_ir605_2830
cd $shome

## but now have to get rid of gas2 within othersdatabase1
## sed '/_g2.dat/ d' othersdatabase1 > othersdatabase  

#####now put the two files together and process them!!!!!!!!
cat waterdatabase_ir605_2830 othersdatabase_ir605_2830 > compdatabase_ir605_2830

#####now put the three files together and process them!!!!!!!!  COMMENT THIS OUT AS NEEDED
cat waterdatabase_ir605_2830 othersdatabase_ir605_2830 co2database_ir605_2830 > compdatabase_ir605_2830

if [ -r compdatabase ]
then
  rm compdatabase
fi
ln -s compdatabase_ir605_2830 compdatabase

######### run code compdatabase.x to produce comp.param
cp ../../UTILITY/compdatabase.x .
if [ -r comp.param ]
then
  rm comp.param
fi
compdatabase.x > dope

mv comp.param ABCcomp_ir605_2830.param
echo " ---->>>> temp results are stored in ABCcomp_ir605_2830.param <<<<----"

## now we have to quit double counting gases in LBL vs XSC
## now we have decided gases 30/81 : gas 30 is too weak
## now we have decided gases 35/61 : gas 35 is too weak
## now we have decided gases 41/80 : gas 41 is too weak
## now we have decided gases 42/54 : gas 42 is too weak, except at 631 cm-1
sed -e '/ 30  /d' -e '/ 35  /d' -e '/ 41  /d' -e '/42     1180/d' ABCcomp_ir605_2830.param > PARAM_HT2016/comp_ir605_2830.param
echo " ---->>>> deleted gases 30 35 41 42 <<<<<<<<<<<<<---------"
echo " ---->>>> FINAL results are stored in PARAM_HT2016/comp_ir605_2830.param <<<<----"

rm compdatabase_ir605_2830 waterdatabase_ir605_2830 comp1.param
rm othersdatabase_ir605_2830 compdatabase.x dope compdatabase ABCcomp_ir605_2830.param
rm co2database_ir605_2830
