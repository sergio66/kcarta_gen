/bin/cp ../INCLUDE/scatter.param  ../INCLUDE/xscatterparam.f; fixcon ../INCLUDE/xscatterparam.f ../INCLUDE/xscatterparam.f90
sed -e "s/.param/param.f90/g" ../INCLUDE/xscatterparam.f90 > ../INCLUDE/scatterparam.f90; rm ../INCLUDE/xscatterparam.f90 ../INCLUDE/xscatterparam.f

/bin/cp ../INCLUDE/kcarta.param  ../INCLUDE/xkcartaparam.f;   fixcon ../INCLUDE/xkcartaparam.f ../INCLUDE/xkcartaparam.f90
sed -e "s/pre_defined.param/pre_definedparam_mars.f90/g" -e "s/post_defined.param/post_definedparam_mars.f90/g" ../INCLUDE/xkcartaparam.f90 > ../INCLUDE/kcartaparam.f90; rm ../INCLUDE/xkcartaparam.f90	../INCLUDE/xkcartaparam.f

/bin/cp ../INCLUDE/pre_defined_mars.param        ../INCLUDE/pre_definedparam_mars.f;            fixcon ../INCLUDE/pre_definedparam_mars.f            ../INCLUDE/pre_definedparam_mars.f90;             rm ../INCLUDE/pre_definedparam_mars.f

/bin/cp ../INCLUDE/pre_defined_mars.param        ../INCLUDE/pre_definedparam_mars.f;            fixcon ../INCLUDE/pre_definedparam_mars.f            ../INCLUDE/pre_definedparam_mars.f90;             rm ../INCLUDE/pre_definedparam_mars.f

/bin/cp ../INCLUDE/post_defined_mars.param       ../INCLUDE/post_definedparam_mars.f;           fixcon ../INCLUDE/post_definedparam_mars.f           ../INCLUDE/post_definedparam_mars.f90;            rm ../INCLUDE/post_definedparam_mars.f

/bin/cp ../INCLUDE/airsheights.param             ../INCLUDE/airsheightsparam_mars.f;            fixcon ../INCLUDE/airsheightsparam_mars.f            ../INCLUDE/airsheightsparam_mars.f90;             rm ../INCLUDE/airsheightsparam_mars.f

/bin/cp ../INCLUDE/airslevels.param              ../INCLUDE/airslevelsparam_mars.f;             fixcon ../INCLUDE/airslevelsparam_mars.f             ../INCLUDE/airslevelsparam_mars.f90;              rm ../INCLUDE/airslevelsparam_mars.f

/bin/cp ../INCLUDE/airslevelheights.param        ../INCLUDE/airslevelheightsparam_mars.f;       fixcon ../INCLUDE/airslevelheightsparam_mars.f       ../INCLUDE/airslevelheightsparam_mars.f90;        rm ../INCLUDE/airslevelheightsparam_mars.f

/bin/cp ../INCLUDE/airsheights_upper.param       ../INCLUDE/airsheights_upperparam_mars.f;      fixcon ../INCLUDE/airsheights_upperparam_mars.f      ../INCLUDE/airsheights_upperparam_mars.f90;       rm ../INCLUDE/airsheights_upperparam_mars.f

/bin/cp ../INCLUDE/airslevels_upper.param        ../INCLUDE/airslevels_upperparam_mars.f;       fixcon ../INCLUDE/airslevels_upperparam_mars.f       ../INCLUDE/airslevels_upperparam_mars.f90;        rm ../INCLUDE/airslevels_upperparam_mars.f

/bin/cp ../INCLUDE/airslevelheights_upper.param  ../INCLUDE/airslevelheights_upperparam_mars.f; fixcon ../INCLUDE/airslevelheights_upperparam_mars.f ../INCLUDE/airslevelheights_upperparam_mars.f90;  rm ../INCLUDE/airslevelheights_upperparam_mars.f

/bin/cp ../INCLUDE/KCARTA_database.param         ../INCLUDE/KCARTA_databaseparam_mars.f;        fixcon ../INCLUDE/KCARTA_databaseparam_mars.f        ../INCLUDE/KCARTA_databaseparam_mars.f90;         rm ../INCLUDE/KCARTA_databaseparam_mars.f

/bin/cp ../INCLUDE/gasIDname.param               ../INCLUDE/gasIDnameparam_mars.f;              fixcon ../INCLUDE/gasIDnameparam_mars.f              ../INCLUDE/gasIDnameparam_mars.f90;               rm ../INCLUDE/gasIDnameparam_mars.f

/bin/cp ../INCLUDE/gauss.param                   ../INCLUDE/gaussparam_mars.f;                  fixcon ../INCLUDE/gaussparam_mars.f                  ../INCLUDE/gaussparam_mars.f90;                   rm ../INCLUDE/gaussparam_mars.f





