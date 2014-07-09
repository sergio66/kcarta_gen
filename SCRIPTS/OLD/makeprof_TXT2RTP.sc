:
### run this script with makeprof.sc inTEXTprofile outRTPfile
### this script uses f77 to take a text sonde point profile, convert it
### to a RTP file for klayers to use; final output is a KLAYERS RTP file

### SYNOPSIS : basic.sc inTEXTprofile outRTPfile
### 
### name of script                           = basic.sc
### path/name of input text point profile    = inprofile
### path/name of RTP layer file from KLAYERS = outRTPfile 

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
HOME=pwd

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

###### make sure that the text point profile --> RTP file converter exists
if [ -x ${BIN_KC}/makeRTPfile.x ]
then
  echo "  "
else 
  echo "####################################################"
  echo "ERROR : makeRTPfile.x does not exist in " ${BIN_KC}
  echo "go to the KCARTA/UTILITY and compile                "
  echo "####################################################"
  exit
fi

######## first run the text point profile thru makeRTPfile.x ################
cp $1 ${WORK}/klayers.sonde.ip
cd ${WORK}
echo "running point profile "$1" thru --> makeRTPfile.x <-- "

echo "klayers.sonde.ip" > makeRTPin
${BIN_KC}/makeRTPfile.x < makeRTPin
/bin/mv makeRTPfile.ip.rtp klayers.ip.rtp

######## then run the RTP point profile thru klayers.x ################
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
/bin/rm klayers.ip.rtp klayers.sonde.ip makeRTPin
#/bin/rm klayers.sonde.ip makeRTPin
/bin/mv klayers.op.rtp $2 

########### come back to SCRIPTS ########################
cd ${SCRIPTS}
