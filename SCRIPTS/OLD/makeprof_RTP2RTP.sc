:
### run this script with makeprof.sc inRTPfile outRTPfile
### this script takes a input level RTP file and runs it thru klayers

### SYNOPSIS : basic.sc inRTPfile outRTPfile
### 
### name of script                           = basic.sc
### path/name of input RTP point profile     = inRTPfile
### path/name of RTP layer file from KLAYERS = outRTPfile (resides in ${WORK})

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


################ then make sure executables exist ################
#### make sure that the klayers executable exists
if [ -x ${BIN_KL}/klayers_airs ]
then
  echo "  "
else 
  echo "########################################################"
  echo "ERROR : KLAYERS executable does not exist in " ${BIN_KL}
  echo "go to the KLAYERS/Src and compile                       "
  echo "########################################################"
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

### go to WORK and rename klayers.op.rtp appropriately
cd ${WORK}
/bin/rm klayers.ip.rtp 
/bin/mv klayers.op.rtp $2 

########### come back to SCRIPTS ########################
cd ${SCRIPTS}

