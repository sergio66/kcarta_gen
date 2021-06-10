/bin/rm ../INCLUDE/TempF90/*.f90 ../INCLUDE/TempF90/*.param

/bin/cp ../INCLUDE/scatter.param  ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f; fixcon ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90
sed -e "s/.param/param.f90/g" ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90 > ../INCLUDE/TempF90/scatterparam.f90; rm ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90 ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f

/bin/cp ../INCLUDE/kcarta.param  ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f;   fixcon ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90
sed -e "s/pre_defined.param/pre_definedparam_mars.f90/g" -e "s/post_defined.param/post_definedparam.f90/g" ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90 > ../INCLUDE/TempF90/kcartaparam.f90; rm ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90	../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f

#########################
#### pre_defined_mars.param copied from pre_defined_orig605_805res.param
/bin/cp ../INCLUDE/pre_defined_mars.param        ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_mars.f;            fixcon ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_mars.f            ../INCLUDE/TempF90/pre_definedparam.f90;             rm ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_mars.f

## post defined same for Earth and Mars
/bin/cp ../INCLUDE/post_defined.param            ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f;                fixcon ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f                ../INCLUDE/TempF90/post_definedparam.f90;                 rm ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f

#########################

/bin/cp ../INCLUDE/TempF90/Mars/airsheights.param             ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_mars.f;            fixcon ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_mars.f            ../INCLUDE/TempF90/airsheightsparam.f90;             rm ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_mars.f

/bin/cp ../INCLUDE/TempF90/Mars/airslevels.param              ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_mars.f;             fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_mars.f             ../INCLUDE/TempF90/airslevelsparam.f90;              rm ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_mars.f

/bin/cp ../INCLUDE/TempF90/Mars/airslevelheights.param        ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_mars.f;       fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_mars.f       ../INCLUDE/TempF90/airslevelheightsparam.f90;        rm ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_mars.f

/bin/cp ../INCLUDE/TempF90/Mars/airsheights_upper.param       ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_mars.f;      fixcon ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_mars.f      ../INCLUDE/TempF90/airsheights_upperparam.f90;       rm ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_mars.f

/bin/cp ../INCLUDE/TempF90/Mars/airslevels_upper.param        ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_mars.f;       fixcon ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_mars.f       ../INCLUDE/TempF90/airslevels_upperparam.f90;        rm ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_mars.f

/bin/cp ../INCLUDE/TempF90/Mars/airslevelheights_upper.param  ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_mars.f; fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_mars.f ../INCLUDE/TempF90/airslevelheights_upperparam.f90;  rm ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_mars.f

###
/bin/cp ../INCLUDE/MARS_database_params_2021/KCARTA_database.param_mars   ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam_mars.f;  fixcon ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam_mars.f        ../INCLUDE/TempF90/KCARTA_databaseparam.f90;         rm ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam_mars.f

/bin/cp ../INCLUDE/MARS_database_params_2021/gasIDname.param_mars         ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam_mars.f;        fixcon ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam_mars.f              ../INCLUDE/TempF90/gasIDnameparam.f90;               rm ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam_mars.f

/bin/cp ../INCLUDE/gauss.param                   ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f;                  fixcon ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f                  ../INCLUDE/TempF90/gaussparam.f90;                   rm ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f

########################################################################
echo "showing files in ../INCLUDE/TempF90/ .. should be todays date"
ls -lt ../INCLUDE/TempF90/*.f90 ../INCLUDE/TempF90/*.param
read -p "Press [Enter] key to continue ..."
