clear

########################################################################
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "files in this dir made by /home/sergio/XYZXYZ/blah.m"

echo " "
echo "Entering lner_EARTH_database_params.sc, which symbolically links kCARTA pressurelevels, heights etc"
echo " "
echo "MAKE SURE YOU ARE in the EARTH_database_para_subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the EARTH_database_para_subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the EARTH_database_para_subdir!!!!!!!!!! at this initial start point!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."

########################################################################
echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd ../TempF90
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the EARTH_database_para_subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the EARTH_database_para_subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the EARTH_database_para_subdir!!!!!!!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
#/bin/rm airsheights_upper.param airsheights.param airslevelheights_upper.param airslevels_upper.param
#/bin/rm KCARTA_database.param airslevelheights.param airslevels.param
if [[ -L airsheights_upper.param ]] ;      then rm airsheights_upper.param ; fi
if [[ -L airsheights.param ]] ;            then rm airsheights.param ; fi
if [[ -L airslevelheights_upper.param ]] ; then rm airslevelheights_upper.param ; fi
if [[ -L airslevels_upper.param ]] ;       then rm airslevels_upper.param ; fi
if [[ -L KCARTA_database.param ]] ;        then rm KCARTA_database.param ; fi
if [[ -L airslevelheights.param ]] ;       then rm airslevelheights.param ; fi
if [[ -L airslevels.param ]] ;             then rm airslevels.param ; fi
if [[ -L gasIDnameparam.param ]] ;         then rm gasIDnameparam.param ; fi

ln -s ../EARTH_database_params/airsheights_upper.param_earth      airsheights_upper.param       
ln -s ../EARTH_database_params/airsheights.param_earth            airsheights.param
ln -s ../EARTH_database_params/airslevelheights_upper.param_earth airslevelheights_upper.param
ln -s ../EARTH_database_params/airslevels_upper.param_earth       airslevels_upper.param
ln -s ../EARTH_database_params/KCARTA_database.param_earth        KCARTA_database.param
ln -s ../EARTH_database_params/airslevelheights.param_earth       airslevelheights.param
ln -s ../EARTH_database_params/airslevels.param_earth             airslevels.param
ln -s ../EARTH_database_params/gasIDname.param_earth              gasIDnameparam.param

########################################################################
echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd Earth
echo "now gone up one directory level, and then to TempF90/Earth ..."
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Earth directory, and NOT NOT NOT NOT in the EARTH_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Earth directory, and NOT NOT NOT NOT in the EARTH_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Earth directory, and NOT NOT NOT NOT in the EARTH_database_param_ subdir!!!!!!!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
#/bin/rm airsheights_upper.param airsheights.param airslevelheights_upper.param airslevels_upper.param
#/bin/rm KCARTA_database.param airslevelheights.param airslevels.param
if [[ -L airsheights_upper.param ]] ;      then rm airsheights_upper.param ; fi
if [[ -L airsheights.param ]] ;            then rm airsheights.param ; fi
if [[ -L airslevelheights_upper.param ]] ; then rm airslevelheights_upper.param ; fi
if [[ -L airslevels_upper.param ]] ;       then rm airslevels_upper.param ; fi
if [[ -L KCARTA_database.param ]] ;        then rm KCARTA_database.param ; fi
if [[ -L airslevelheights.param ]] ;       then rm airslevelheights.param ; fi
if [[ -L airslevels.param ]] ;             then rm airslevels.param ; fi
if [[ -L gasIDnameparam.param ]] ;         then rm gasIDnameparam.param ; fi

ln -s ../../EARTH_database_params/airsheights_upper.param_earth      airsheights_upper.param       
ln -s ../../EARTH_database_params/airsheights.param_earth            airsheights.param
ln -s ../../EARTH_database_params/airslevelheights_upper.param_earth airslevelheights_upper.param
ln -s ../../EARTH_database_params/airslevels_upper.param_earth       airslevels_upper.param
ln -s ../../EARTH_database_params/KCARTA_database.param_earth        KCARTA_database.param
ln -s ../../EARTH_database_params/airslevelheights.param_earth       airslevelheights.param
ln -s ../../EARTH_database_params/airslevels.param_earth             airslevels.param
ln -s ../../EARTH_database_params/gasIDname_earth.param              gasIDnameparam.param

########################################################################

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd ../../EARTH_database_params
echo " "
echo " current dir = "
pwd
echo "Done"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
