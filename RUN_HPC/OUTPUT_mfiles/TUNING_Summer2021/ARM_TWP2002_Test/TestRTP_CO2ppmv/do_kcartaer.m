rmer = ['!/bin/rm rad0.dat rad1.dat rad2.dat rad3.dat warning.msg']; eval(rmer);

use_this_rtp = 'raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp'; 
  sedder = ['!sed -e "s/XYZXYZ/'  use_this_rtp   '/g" template_Qrad.nml  >  outnml0'];  eval(sedder); kcartaer = ['!time ' kcarta  ' outnml0  rad0.dat'];                   eval(kcartaer);
  mver = ['!/bin/mv warning.msgMMM warning.msg0']; eval(mver); 
  grepper = ['!grep -in ''CO2 profile in rtp file is'' warning.msg0']; eval(grepper)
  grepper = ['!grep -in ''Planet Earth : CO2 colavg ppm'' warning.msg0']; eval(grepper)
disp('ret to continue'); pause

use_this_rtp = 'set_gunit10_raw_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp';
  sedder = ['!sed -e "s/XYZXYZ/'  use_this_rtp   '/g" template_Qrad.nml  >  outnml1'];  eval(sedder); kcartaer = ['!time ' kcarta  ' outnml1  rad1.dat'];                   eval(kcartaer);
  mver = ['!/bin/mv warning.msgMMM warning.msg1']; eval(mver);
  grepper = ['!grep -in ''CO2 profile in rtp file is'' warning.msg1']; eval(grepper)
  grepper = ['!grep -in ''Planet Earth : CO2 colavg ppm'' warning.msg1']; eval(grepper)
disp('ret to continue'); pause

use_this_rtp = 'set_gunit10_co2ppm400_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp';
  sedder = ['!sed -e "s/XYZXYZ/'  use_this_rtp   '/g" template_Qrad.nml  >  outnml2'];  eval(sedder); kcartaer = ['!time ' kcarta  ' outnml2  rad2.dat'];                   eval(kcartaer);
  mver = ['!/bin/mv warning.msgMMM warning.msg2']; eval(mver); 
  grepper = ['!grep -in ''CO2 profile in rtp file is'' warning.msg2']; eval(grepper)
  grepper = ['!grep -in ''Planet Earth : CO2 colavg ppm'' warning.msg2']; eval(grepper)
disp('ret to continue'); pause

use_this_rtp = 'set_gas_2_430_armtwp_nearby_semiclear_jan04_resetT_nte_co2.rp.rtp';
  sedder = ['!sed -e "s/XYZXYZ/'  use_this_rtp   '/g" template_Qrad.nml  >  outnml3'];  eval(sedder); kcartaer = ['!time ' kcarta  ' outnml3  rad3.dat'];                   eval(kcartaer);
  mver = ['!/bin/mv warning.msgMMM warning.msg3']; eval(mver); 
  grepper = ['!grep -in ''CO2 profile in rtp file is'' warning.msg3']; eval(grepper)
  grepper = ['!grep -in ''Planet Earth : CO2 colavg ppm'' warning.msg3']; eval(grepper)
disp('ret to continue'); pause
