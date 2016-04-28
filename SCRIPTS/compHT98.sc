:
# this makes comp.param and xsec.param, which contains lists of the 
# available compressed data files. Make sure you look at INCLUDE/kcarta.param
# and set up the correct paths to kWaterPath,kCompPath and 
# kXsecParamFile,kCompParamFile

### WARNING (1) : 
### look at kcarta.param and set correct paths to kWaterPath,kCompPath
 
shome=`pwd`

#first do the water
#cd /asl/data/kcarta/v20.ieee-le/h2o.ieee-le                   ##path to kWaterPath
cd /asl/data/kcarta/H2012.ieee-le/IR605/lblrtm2/h2o.ieee-le/  ##path to kWaterPath
ls -1  >& $shome/waterdatabase
cd $shome

#then do CO2, if it is separate --- else just stick to the empty file (from "rm" and "touch")
doCO2=1
doCO2=-1
if [ -r $shome/co2database ]
then
  rm $shome/co2database
fi
touch $shome/co2database
if [ $doCO2 -gt 0 ]
then
  cd $shome
  cd /asl/data/kcarta/v24.ieee-le/co2.ieee-le    ##path to kCO2Path
  ls -1  >& $shome/co2database
  cd $shome
fi

#then do the rest of the gases
cd $shome
#cd /asl/data/kcarta/v20.ieee-le/etc.ieee-le                   ##path to kCompPath
cd /asl/data/kcarta/H2012.ieee-le/IR605/lblrtm2/etc.ieee-le/  ##path to kCompPath
ls -1  >& $shome/othersdatabase1
cd $shome

## but if we have separate CO2 database, have to get rid of gas2 within othersdatabase1
if [ "$doCO2" -gt "0" ]
then
  sed '/_g2.dat/ d' $shome/othersdatabase1 > $shome/othersdatabase
else
  cp $shome/othersdatabase1 $shome/othersdatabase
fi

#####now put the two files together and process them!!!!!!!!
cat waterdatabase co2database othersdatabase > compdatabase
cp ../UTILITY/compdatabase.x .

if [ -r comp.param ]
then
  rm comp.param
fi

compdatabase.x > dope
rm comp1.param compdatabase waterdatabase othersdatabase compdatabase.x dope
rm co2database othersdatabase1
mv comp.param comp0.param

########break up the files into comp.param and xsec.param
#awk '($1 <= 50)' comp0.param > comp107.param
#awk '($1 >= 51)' comp0.param > xsec107.param
awk '($1 <= 50)' comp0.param > compNEW.param
awk '($1 >= 51)' comp0.param > xsecNEW.param

### WARNING (2) : 
### look at kcarta.param and set correct paths to kXsecParamFile,kCompParamFile

# copy the files to where they are needed; assume that it will pretty much
# be in /asl/data/kcarta/KCARTADATA/General/ for time eternal
# save as files that reflect the database used eg H1998, H2000  
# the names are semi-logical eg H1996  means "built with HITRAN 1996"
#                            eg H1998  means "built with HITRAN 1998"
#                            eg HT1998 means "built with HITRAN 1998 + Toth"
# edit as needed
#
# orig stuff, 1998
#cp comp107.param                              ../SRC/compHT1998.param
#cp comp107.param /asl/data/kcarta/KCARTADATA/General/compHT1998.param
#cp xsec107.param                              ../SRC/xsecHT1998.param
#cp xsec107.param /asl/data/kcarta/KCARTADATA/General/xsecHT1998.param
# 
# lblrtm based IR database
cat compNEW.param xsecNEW.param > /home/sergio/KCARTA/SCRIPTS/MAKE_COMP_HTXY_PARAM_SC/PARAM_TEMP/testH2012_lblrtm

rm comp0.param
#rm comp107.param xsec107.param

