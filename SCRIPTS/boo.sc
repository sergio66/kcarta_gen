:
shome=`pwd` 
#then do CO2 
cd $shome 
cd /asl/data/kcarta/v24.ieee-le/co2.ieee-le    ##path to kCO2Path 
ls -1  >& $shome/co2database 
cd $shome 

#then do the rest of the gases 
cd $shome 
cd /asl/data/kcarta/v20.ieee-le/etc.ieee-le    ##path to kCompPath 
ls -1  >& $shome/othersdatabase1 
cd $shome 
 
## but now have to get rid of gas2 within othersdatabase1 
sed '/_g2.dat/ d' othersdatabase1 > othersdatabase 
