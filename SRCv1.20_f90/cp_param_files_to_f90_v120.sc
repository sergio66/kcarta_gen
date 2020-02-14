/bin/cp ../INCLUDE/scatter.param  ../INCLUDE/xscatterparam.f; fixcon ../INCLUDE/xscatterparam.f ../INCLUDE/xscatterparam.f90
sed -e "s/.param/param.f90/g" ../INCLUDE/xscatterparam.f90 > ../INCLUDE/scatterparam.f90; rm ../INCLUDE/xscatterparam.f90 ../INCLUDE/xscatterparam.f

/bin/cp ../INCLUDE/kcarta.param  ../INCLUDE/xkcartaparam.f;   fixcon ../INCLUDE/xkcartaparam.f ../INCLUDE/xkcartaparam.f90
sed -e "s/pre_defined_orig605_805res.param/pre_definedparam.f90/g" -e "s/post_defined.param/post_definedparam.f90/g" ../INCLUDE/xkcartaparam.f90 > ../INCLUDE/kcartaparam.f90; rm ../INCLUDE/xkcartaparam.f90	../INCLUDE/xkcartaparam.f

#/bin/cp ../INCLUDE/pre_definedHIGHRES_IR.param   ../INCLUDE/pre_definedHIGHRES_IRparam.f;  fixcon ../INCLUDE/pre_definedHIGHRES_IRparam.f  ../INCLUDE/pre_definedHIGHRES_IRparam.f90; rm ../INCLUDE/pre_definedHIGHRES_IRparam.f

/bin/cp ../INCLUDE/pre_defined_orig605_805res.param             ../INCLUDE/pre_definedparam.f;            fixcon ../INCLUDE/pre_definedparam.f            ../INCLUDE/pre_definedparam.f90;           rm ../INCLUDE/pre_definedparam.f

/bin/cp ../INCLUDE/post_defined.param            ../INCLUDE/post_definedparam.f;           fixcon ../INCLUDE/post_definedparam.f           ../INCLUDE/post_definedparam.f90;           rm ../INCLUDE/post_definedparam.f

/bin/cp ../INCLUDE/airsheights.param             ../INCLUDE/airsheightsparam.f;            fixcon ../INCLUDE/airsheightsparam.f            ../INCLUDE/airsheightsparam.f90;            rm ../INCLUDE/airsheightsparam.f

/bin/cp ../INCLUDE/airslevels.param              ../INCLUDE/airslevelsparam.f;             fixcon ../INCLUDE/airslevelsparam.f             ../INCLUDE/airslevelsparam.f90;             rm ../INCLUDE/airslevelsparam.f

/bin/cp ../INCLUDE/airslevelheights.param        ../INCLUDE/airslevelheightsparam.f;       fixcon ../INCLUDE/airslevelheightsparam.f       ../INCLUDE/airslevelheightsparam.f90;       rm ../INCLUDE/airslevelheightsparam.f

/bin/cp ../INCLUDE/airsheights_upper.param       ../INCLUDE/airsheights_upperparam.f;      fixcon ../INCLUDE/airsheights_upperparam.f      ../INCLUDE/airsheights_upperparam.f90;      rm ../INCLUDE/airsheights_upperparam.f

/bin/cp ../INCLUDE/airslevels_upper.param        ../INCLUDE/airslevels_upperparam.f;       fixcon ../INCLUDE/airslevels_upperparam.f        ../INCLUDE/airslevels_upperparam.f90;      rm ../INCLUDE/airslevels_upperparam.f

/bin/cp ../INCLUDE/airslevelheights_upper.param  ../INCLUDE/airslevelheights_upperparam.f; fixcon ../INCLUDE/airslevelheights_upperparam.f ../INCLUDE/airslevelheights_upperparam.f90; rm ../INCLUDE/airslevelheights_upperparam.f

/bin/cp ../INCLUDE/kcarta_jpl.param              ../INCLUDE/kcarta_jplparam.f;             fixcon ../INCLUDE/kcarta_jplparam.f             ../INCLUDE/kcarta_jplparam.f90;             rm ../INCLUDE/kcarta_jplparam.f

/bin/cp ../INCLUDE/KCARTA_database.param         ../INCLUDE/KCARTA_databaseparam.f;        fixcon ../INCLUDE/KCARTA_databaseparam.f        ../INCLUDE/KCARTA_databaseparam.f90;        rm ../INCLUDE/KCARTA_databaseparam.f

/bin/cp ../INCLUDE/gasIDname.param               ../INCLUDE/gasIDnameparam.f;              fixcon ../INCLUDE/gasIDnameparam.f              ../INCLUDE/gasIDnameparam.f90;              rm ../INCLUDE/gasIDnameparam.f

/bin/cp ../INCLUDE/gauss.param                   ../INCLUDE/gaussparam.f;                  fixcon ../INCLUDE/gaussparam.f                  ../INCLUDE/gaussparam.f90;                  rm ../INCLUDE/gaussparam.f





