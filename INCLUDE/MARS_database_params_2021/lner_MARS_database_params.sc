echo " "
echo "Entering lner_MARS_database_params.sc, which symbolically links kCARTA pressurelevels, heights etc"
echo " "
echo "MAKE SURE YOU ARE in the MARS_ subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the MARS_ subdir!!!!!!!!!! at this initial start point!!!!"
echo "MAKE SURE YOU ARE in the MARS_ subdir!!!!!!!!!! at this initial start point!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."
cd ../

echo "MAKE SURE YOU ARE in the INCLUDE directory, and NOT NOT NOT NOT in the MARS_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE directory, and NOT NOT NOT NOT in the MARS_ subdir!!!!!!!!!!"
echo "MAKE SURE YOU ARE in the INCLUDE directory, and NOT NOT NOT NOT in the MARS_ subdir!!!!!!!!!!"
echo " "
echo " current dir = "
pwd
read -p "Press [Enter] key to continue ..."

/bin/rm airsheights_upper.param airsheights.param airslevelheights_upper.param airslevels_upper.param
/bin/rm KCARTA_database.param airslevelheights.param airslevels.param

ln -s MARS_database_params/airsheights_upper.param_mars      airsheights_upper.param       
ln -s MARS_database_params/airsheights.param_mars            airsheights.param
ln -s MARS_database_params/airslevelheights_upper.param_mars airslevelheights_upper.param
ln -s MARS_database_params/airslevels_upper.param_mars       airslevels_upper.param
ln -s MARS_database_params/KCARTA_database.param_mars        KCARTA_database.param
ln -s MARS_database_params/airslevelheights.param_mars       airslevelheights.param
ln -s MARS_database_params/airslevels.param_mars             airslevels.param

cd MARS_database_params
echo " "
echo " current dir = "
pwd
echo "Done"
