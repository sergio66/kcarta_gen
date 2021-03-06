:
# this makes comp.param and xsec.param, which contains lists of the 
# available compressed data files. Make sure you look at INCLUDE/kcarta.param
# and set up the correct paths to kWaterPath,kCompPath and 
# kXsecParamFile,kCompParamFile

### WARNING (1) : 
### look at kcarta.param and set correct paths to kWaterPath,kCompPath
 
shome=`pwd`

#first do the water
cd /asl/data/kcarta/v10.ieee-le/h2o.ieee-le    ##path to kWaterPath
ls -1  >& $shome/waterdatabase
cd $shome

#then do the rest of the gases
cd $shome
cd /asl/data/kcarta/v10.ieee-le/etc.ieee-le    ##path to kCompPath
ls -1  >& $shome/othersdatabase
cd $shome

#####now put the two files together and process them!!!!!!!!
cat waterdatabase othersdatabase > compdatabase
cp ../BIN/compdatabase.x .

if [ -r comp.param ]
then
  rm comp.param
fi

compdatabase.x > dope
rm comp1.param compdatabase waterdatabase othersdatabase compdatabase.x dope
mv comp.param comp0.param

########break up the files into comp.param and xsec.param
awk '($1 <= 50)' comp0.param > comp107.param
awk '($1 >= 51)' comp0.param > xsec107.param

### WARNING (2) : 
### look at kcarta.param and set correct paths to kXsecParamFile,kCompParamFile

# copy the files to where they are needed; assume that it will pretty much
# be in /asl/data/kcarta/KCARTADATA/General/ for time eternal
# save as files that reflect the database used eg H1998, H2000  
# the names are semi-logical eg H1996  means "built with HITRAN 1996"
#                            eg H1998  means "built with HITRAN 1998"
#                            eg HT1998 means "built with HITRAN 1998 + Toth"
cp comp107.param                              ../SRC/compH1992.param
cp comp107.param /asl/data/kcarta/KCARTADATA/General/compH1992.param
cp xsec107.param                              ../SRC/xsecH1992.param
cp xsec107.param /asl/data/kcarta/KCARTADATA/General/xsecH1992.param
rm comp0.param comp107.param xsec107.param

