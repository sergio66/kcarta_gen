#clear

#rm ../INCLUDE/kcarta.param

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
if [[ -L airsTZ_STD.param ]] ;             then rm airsTZ_STD.param ; fi
if [[ -L gasIDnameparam.param ]] ;         then rm gasIDnameparam.param ; fi

echo "linking airsheights_upper.param"
ln -s ../EARTH_database_params/airsheights_upper.param_earth      airsheights_upper.param       
echo "linking airsheight.param"
ln -s ../EARTH_database_params/airsheights.param_earth            airsheights.param
echo "linking airslevlheights_upper.param"
ln -s ../EARTH_database_params/airslevelheights_upper.param_earth airslevelheights_upper.param
echo "linking airslevels_upper.param"
ln -s ../EARTH_database_params/airslevels_upper.param_earth       airslevels_upper.param
echo "linking KCARTA_database.param"
ln -s ../EARTH_database_params/KCARTA_database.param_earth        KCARTA_database.param
echo "linking airslevelheights.param"
ln -s ../EARTH_database_params/airslevelheights.param_earth       airslevelheights.param
echo "linking airslevels.param"
ln -s ../EARTH_database_params/airslevels.param_earth             airslevels.param
echo "linking airsTZ_STD.param"
ln -s ../EARTH_database_params/airsTZ_STD.param_earth             airsTZ_STD.param
echo "linking gasIDnameparam.param"
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
if [[ -L airsTZ_STD.param ]] ;             then rm airsTZ_STD.param ; fi
if [[ -L gasIDnameparam.param ]] ;         then rm gasIDnameparam.param ; fi

echo "linking airsheights_upper.param"
ln -s ../../EARTH_database_params/airsheights_upper.param_earth      airsheights_upper.param       
echo "linking airsheight.param"
ln -s ../../EARTH_database_params/airsheights.param_earth            airsheights.param
echo "linking airslevlheights_upper.param"
ln -s ../../EARTH_database_params/airslevelheights_upper.param_earth airslevelheights_upper.param
echo "linking airslevels_upper.param"
ln -s ../../EARTH_database_params/airslevels_upper.param_earth       airslevels_upper.param
echo "linking KCARTA_database.param"
ln -s ../../EARTH_database_params/KCARTA_database.param_earth        KCARTA_database.param
echo "linking airslevelheights.param"
ln -s ../../EARTH_database_params/airslevelheights.param_earth       airslevelheights.param
echo "linking airslevels.param"
ln -s ../../EARTH_database_params/airslevels.param_earth             airslevels.param
echo "linking airsTZ_STD.param"
ln -s ../../EARTH_database_params/airsTZ_STD.param_earth             airsTZ_STD.param
echo "linking gasIDnameparam.param"
ln -s ../../EARTH_database_params/gasIDname.param_earth              gasIDnameparam.param

########################################################################

echo " "
echo ">>>>>>>>>>>>>>>>>>>>>>>>>"
cd ../../EARTH_database_params
echo " "
echo " current dir = "
pwd
echo "   Now going to run cp_param_files_to_f90_v122.sc ... which translates the .param_earth files to .f90 files"
echo "   Now going to run cp_param_files_to_f90_v122.sc ... which translates the .param_earth files to .f90 files"
echo "   Now going to run cp_param_files_to_f90_v122.sc ... which translates the .param_earth files to .f90 files"
echo " "
echo "Done"
echo ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
read -p "Press [Enter] key to continue ..."
