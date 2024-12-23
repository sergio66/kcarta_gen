pwd
/bin/rm ../INCLUDE/TempF90/*.f90 ../INCLUDE/TempF90/*.param

/bin/cp ../INCLUDE/scatter.param  ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f; fixcon ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90
sed -e "s/.param/param.f90/g" ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90 > ../INCLUDE/TempF90/scatterparam.f90; rm ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f90 ../INCLUDE/TempF90/JunkTempDir/xscatterparam.f

/bin/cp ../INCLUDE/kcarta.param  ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f;   fixcon ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90
sed -e "s/pre_defined.param/pre_definedparam.f90/g" -e "s/post_defined.param/post_definedparam.f90/g" ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90 > ../INCLUDE/TempF90/kcartaparam.f90; rm ../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f90	../INCLUDE/TempF90/JunkTempDir/xkcartaparam.f

#########################
#/bin/cp ../INCLUDE/pre_definedHIGHRES_IR.param   ../INCLUDE/TempF90/JunkTempDir/pre_definedHIGHRES_IRparam.f;  fixcon ../INCLUDE/TempF90/JunkTempDir/pre_definedHIGHRES_IRparam.f  ../INCLUDE/TempF90/pre_definedHIGHRES_IRparam.f90; rm ../INCLUDE/TempF90/JunkTempDir/pre_definedHIGHRES_IRparam.f

/bin/cp ../INCLUDE/pre_defined_orig605_805res.param  ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_orig605_805res.f;   fixcon ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_orig605_805res.f  ../INCLUDE/TempF90/pre_definedparam_orig605_805res.f90; rm ../INCLUDE/TempF90/JunkTempDir/pre_definedparam_orig605_805res.f

/bin/cp ../INCLUDE/pre_defined.param             ../INCLUDE/TempF90/JunkTempDir/pre_definedparam.f;            fixcon ../INCLUDE/TempF90/JunkTempDir/pre_definedparam.f            ../INCLUDE/TempF90/pre_definedparam.f90;           rm ../INCLUDE/TempF90/JunkTempDir/pre_definedparam.f

/bin/cp ../INCLUDE/post_defined.param            ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f;           fixcon ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f           ../INCLUDE/TempF90/post_definedparam.f90;           rm ../INCLUDE/TempF90/JunkTempDir/post_definedparam.f

#########################

/bin/cp ../INCLUDE/TempF90/Earth/airsTZ_STD.param             ../INCLUDE/TempF90/JunkTempDir/airsTZ_STDparam_earth.f;            fixcon ../INCLUDE/TempF90/JunkTempDir/airsTZ_STDparam_earth.f            ../INCLUDE/TempF90/airsTZ_STDparam.f90;             rm ../INCLUDE/TempF90/JunkTempDir/airsTZ_STDparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airsheights.param             ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_earth.f;            fixcon ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_earth.f            ../INCLUDE/TempF90/airsheightsparam.f90;             rm ../INCLUDE/TempF90/JunkTempDir/airsheightsparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airslevels.param              ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_earth.f;             fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_earth.f             ../INCLUDE/TempF90/airslevelsparam.f90;              rm ../INCLUDE/TempF90/JunkTempDir/airslevelsparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airslevelheights.param        ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_earth.f;       fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_earth.f       ../INCLUDE/TempF90/airslevelheightsparam.f90;        rm ../INCLUDE/TempF90/JunkTempDir/airslevelheightsparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airsheights_upper.param       ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_earth.f;      fixcon ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_earth.f      ../INCLUDE/TempF90/airsheights_upperparam.f90;       rm ../INCLUDE/TempF90/JunkTempDir/airsheights_upperparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airslevels_upper.param        ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_earth.f;       fixcon ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_earth.f       ../INCLUDE/TempF90/airslevels_upperparam.f90;        rm ../INCLUDE/TempF90/JunkTempDir/airslevels_upperparam_earth.f

/bin/cp ../INCLUDE/TempF90/Earth/airslevelheights_upper.param  ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_earth.f; fixcon ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_earth.f ../INCLUDE/TempF90/airslevelheights_upperparam.f90;  rm ../INCLUDE/TempF90/JunkTempDir/airslevelheights_upperparam_earth.f

#########################
/bin/cp ../INCLUDE/EARTH_database_params/KCARTA_database.param_earth   ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam.f;        fixcon ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam.f        ../INCLUDE/TempF90/KCARTA_databaseparam.f90;     rm ../INCLUDE/TempF90/JunkTempDir/KCARTA_databaseparam.f

/bin/cp ../INCLUDE/EARTH_database_params/gasIDname.param_earth               ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam.f;              fixcon ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam.f              ../INCLUDE/TempF90/gasIDnameparam.f90;              rm ../INCLUDE/TempF90/JunkTempDir/gasIDnameparam.f

#/bin/cp ../INCLUDE/gauss.param                   ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f;                  fixcon ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f                  ../INCLUDE/TempF90/gaussparam.f90;                   rm ../INCLUDE/TempF90/JunkTempDir/gaussparam_mars.f
/bin/cp ../INCLUDE/gauss.param                   ../INCLUDE/TempF90/JunkTempDir/gaussparam_earth.f;                  fixcon ../INCLUDE/TempF90/JunkTempDir/gaussparam_earth.f                  ../INCLUDE/TempF90/gaussparam.f90;                   rm ../INCLUDE/TempF90/JunkTempDir/gaussparam_earth.f

########################################################################
echo "showing files in ../INCLUDE/TempF90/ .. should be todays date"
ls -lt ../INCLUDE/TempF90/*.f90 ../INCLUDE/TempF90/*.param
read -p "Press [Enter] key to continue ..."
