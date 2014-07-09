:
:
### run this script with basic.sc inRTPfile outkcBINfile
### this script takes a RTP sonde point profile for klayers to use; after 
### that KCARTA is run; after that the kcarta readers are called to produce 
### a simple binary file that can easily be read in

### SYNOPSIS : basic.sc inRTPfile outkcBINfile
### 
### name of script                        = basic.sc
### path/name of input RTP point profile  = inRTPfile
### path/name of output KCARTA binary     = ../WORK/outkcBINfile
###
### other files created (in ../WORK subdir) are
###    klayers.op.rtp  RTP layer file from KLAYERS
###    basic110.nml    namelist file to drive KCARTA (from SCRIPTS/BASIC)
###    warning.msg     error or chirpy messages that come from kCARTA
###    outfile.bin     file that come out of readkcBasic.x

### assumes the following subdir structure
### ---- KLAYERS --- Bin
### ---- KCARTA ---- WORK
###         |------- SCRIPTS
###         |------- BIN
###         |------- UTILITY

### you can change things if necessary

SCRIPTS="/home/sergio/KCARTA/SCRIPTS"     #this is the SCRIPTS directory
WORK="/home/sergio/KCARTA/WORK"           #this is the WORK directory
BIN_KC="/home/sergio/KCARTA/BIN"          #this is the KCARTA binary directory
BIN_KL="/home/sergio/KLAYERS/Bin"         #this is the KLAYERS binary subdir

clear

################ first make sure files and directories are OK ################
####check to see if the WORK subdir exits

if [ -d ${WORK} ]
then
#  echo "../WORK directory exists .. no need to create"
  echo " "
else 
  echo "########################################"
  echo "ERROR : no WORK subdirectory ... exiting"
  echo "########################################"
  exit
fi

######check to see if input profile for klayers.x $1 exists
if [ -r "$1" ]
then
  echo " "
else
  echo "####################################################"
  echo "kLAYERS input file " $1 "does not exist .... exiting"
  echo "####################################################"
  exit
fi

######check to see if kcarta.x output file $2 exists
if [ -r "${WORK}/$2" ]
then
  echo "##################################################################"
  echo "kCARTA output file " ${WORK}/$2 "already exists! Cannot overwrite!"
  echo "##################################################################"
  exit
fi

######check to see if readkcBasic.x output file exists
if [ -r ${WORK}/outfile.bin ]

then
  echo "##############################################################"
  echo "readkcBasic output file " ${WORK} "outfile.bin already exists!"
  echo "##############################################################"
  echo "Cannot overwrite!"
  exit
fi

################ then make sure executables exist ################
#### make sure that the klayers executable exists
if [ -x ${BIN_KL}/klayers_airs ]
then
  echo "  "
else 
  echo "########################################################"
  echo "ERROR : KLAYERS executable does not exist in " ${BIN_KL}
  echo "go to the KLAYERS/Src and compile
  echo "########################################################"
  exit
fi

###### make sure that the KCARTA code exists
if [ -x ${BIN_KC}/bkcarta.x ]
then
  echo "  "
else 
  echo "####################################################"
  echo "ERROR : bkcarta.x does not exist in " ${BIN_KC}
  echo "go to the KCARTA/SRC and compile
  echo "####################################################"
  exit
fi

###### make sure that the reader exists
if [ -x ${BIN_KC}/readkcBasic.x ]
then
  echo "  "
else 
  echo "####################################################"
  echo "ERROR : readkcBasic.x does not exist in " ${BIN_KC}
  echo "go to the KCARTA/UTILITY and compile
  echo "####################################################"
  exit
fi

######## run the RTP point profile thru klayers.x ################
cp $1 ${WORK}/klayers.ip.rtp
cd ${WORK}
cp ${SCRIPTS}/BASIC/klayers110 .
SEDSTR='s/\//\\\//g'
REALOUT1=`echo ${BIN_KL} | sed $SEDSTR`
REALOUT2=`echo ${WORK}   | sed $SEDSTR`
sed -e "s/LAY/$REALOUT1/g" -e "s/WRK/$REALOUT2/g" klayers110 > klayers110run

cp klayers110run ${BIN_KL}/.
/bin/rm klayers110run klayers110

cd ${BIN_KL}
chmod +x klayers110run
more klayers110run
echo "running sonde rtp profile "$1" thru --> klayers_airs <-- "
klayers110run
/bin/rm klayers110run

########### now run the basic.nml file ####################
echo "     "
echo "     "
echo "     "
cd ${WORK}

/bin/rm warning.msg  2> /dev/null !! true
/bin/cp ${SCRIPTS}/BASIC/basic110 .
/bin/cp ${SCRIPTS}/BASIC/basic110.nml .

SEDSTR='s/\//\\\//g'
REALOUT=`echo $2 | sed $SEDSTR`
sed -e "s/outkcfile/$REALOUT/g" basic110 > kcartabasic110
echo "running --> kcarta <-- to produce outputfile "$2
chmod +x kcartabasic110
kcartabasic110
echo "display last few lines of warning.msg, which tells if kCARTA worked!!"
tail -3 warning.msg
/bin/rm basic110 kcartabasic110 klayers.ip.rtp klayers.sonde.* makeRTPin

########### now read in the file and reshape the matrix ####################
echo "     "
echo "     "
echo "     "

/bin/rm gnuplotfile.txt  2> /dev/null !! true
/bin/rm outfile.bin  2> /dev/null !! true
/bin/cp ${SCRIPTS}/BASIC/basic3 .

SEDSTR='s/\//\\\//g'
REALREAD=`echo $2 | sed $SEDSTR`
sed -e "s/outkcfile/$REALREAD/g" basic3 > readkcBasic3

echo "reshaping file to binary file --> outfile.bin <-- "
chmod +x readkcBasic3
readkcBasic3
/bin/rm basic3 readkcBasic3


######### if you want to display text file with gnuplot
if [ -r /usr/bin/gnuplot ]
then
  /bin/rm gnuplotfile.txt  2> /dev/null !! true
  /bin/cp ${SCRIPTS}/BASIC/basic3txt .

  SEDSTR='s/\//\\\//g'
  REALREAD=`echo $2 | sed $SEDSTR`
  sed -e "s/outkcfile/$REALREAD/g" basic3txt > readkcBasic3txt

  echo "reshaping first column of binary file --> gnuplotfile.txt <-- "
  chmod +x readkcBasic3txt
  readkcBasic3txt
  gnuplot ${SCRIPTS}/BASIC/gnuplotter
  /bin/rm basic3txt readkcBasic3txt gnuplotfile.txt
fi
########### come back to SCRIPTS ########################
cd ${SCRIPTS}
