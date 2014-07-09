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

SCRIPTS="/home/sergio/KCARTA/SCRIPTS"     #the SCRIPTS directory
#BIN_KC="/home/sergio/KCARTA/BIN"          #the KCARTA binary directory
#BIN_KL="/home/sergio/KLAYERS/Bin"         #the KLAYERS binary subdir
BIN_KC="/home/sergio/KCARTA/BIN"          #the KCARTA binary directory
BIN_KL="/asl/packages/klayers/Bin/"       #the KLAYERS binary subdir
HOME=$(pwd)                               #the HOME dir, where you are

clear

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
cp $1 ${HOME}/klayers.sonde.ip
cd ${HOME}
echo "running point profile "$1" thru --> makeRTPfile.x <-- "

echo "klayers.sonde.ip" > makeRTPin
${BIN_KC}/makeRTPfile.x < makeRTPin

pwd
/bin/mv makeRTPfile.ip.rtp klayers.ip.rtp
ls -lt klayers.ip.rtp 
echo "ran point profile "$1" thru --> makeRTPfile.x <-- "
echo " "
echo " "
echo " "
echo " "
echo " "
echo " "

######## then run the RTP point profile thru klayers.x ################
cd ${HOME}
cp ${SCRIPTS}/BASIC/klayers110 .
SEDSTR='s/\//\\\//g'
REALOUT1=`echo ${BIN_KL} | sed $SEDSTR`
REALOUT2=`echo ${HOME}   | sed $SEDSTR`
sed -e "s/LAY/$REALOUT1/g" -e "s/WRK/$REALOUT2/g" klayers110 > klayers110run

chmod +x klayers110run
more klayers110run
echo "running sonde rtp profile "$1" thru --> klayers_airs <-- "
./klayers110run
/bin/rm klayers110run klayers110

### go to HOME and rename klayers.op.rtp appropriately
cd ${HOME}
/bin/rm klayers.ip.rtp klayers.sonde.ip makeRTPin
#/bin/rm klayers.sonde.ip makeRTPin
/bin/mv klayers.op.rtp $2 

########### come back to SCRIPTS ########################
cd ${SCRIPTS}

########### come back to HOME ########################
cd ${HOME}
