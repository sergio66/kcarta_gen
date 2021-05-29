clear

########################################################################
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "files in this dir made by /home/sergio/HITRAN2UMBCLBL/MARS_MAKEIR/Develop_2021/produce_kcarta_paramfiles_2021.m"

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo "Starting lner_MARS_database_params.sc, which symbolically links kCARTA pressurelevels, heights etc"
echo " "
echo "MAKE SURE YOU ARE in the MARS_database_param_ subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the MARS_database_param_ subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the MARS_database_param_ subdir!!!!!!!!!! at this initial start point!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."

########################################################################
echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd ../TempF90
echo "now gone up one directory level, and then to TempF90 ..."
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90 directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
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

ln -s ../MARS_database_params/airsheights_upper.param_mars      airsheights_upper.param       
ln -s ../MARS_database_params/airsheights.param_mars            airsheights.param
ln -s ../MARS_database_params/airslevelheights_upper.param_mars airslevelheights_upper.param
ln -s ../MARS_database_params/airslevels_upper.param_mars       airslevels_upper.param
ln -s ../MARS_database_params/KCARTA_database.param_mars        KCARTA_database.param
ln -s ../MARS_database_params/airslevelheights.param_mars       airslevelheights.param
ln -s ../MARS_database_params/airslevels.param_mars             airslevels.param
ln -s ../MARS_database_params/gasIDname.param_mars              gasIDnameparam.param

########################################################################
echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd Mars
echo "now gone up one directory level, and then to TempF90/Mars ..."
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Mars directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Mars directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE/TempF90/Mars directory, and NOT NOT NOT NOT in the MARS_database_param_ subdir!!!!!!!!!!"
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

ln -s ../../MARS_database_params/airsheights_upper.param_mars      airsheights_upper.param       
ln -s ../../MARS_database_params/airsheights.param_mars            airsheights.param
ln -s ../../MARS_database_params/airslevelheights_upper.param_mars airslevelheights_upper.param
ln -s ../../MARS_database_params/airslevels_upper.param_mars       airslevels_upper.param
ln -s ../../MARS_database_params/KCARTA_database.param_mars        KCARTA_database.param
ln -s ../../MARS_database_params/airslevelheights.param_mars       airslevelheights.param
ln -s ../../MARS_database_params/airslevels.param_mars             airslevels.param
ln -s ../../MARS_database_params/gasIDname.param_mars              gasIDnameparam.param

########################################################################

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd ../../MARS_database_params
echo " "
echo " current dir = "
pwd
echo "Done"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
