#!/bin/sh
#
# this is a simple demo, normally this would be run with a loop
### if you have runaway processes, do
### killall bkcarta.x.CHINEW; killall kcwrap2; killall kcwrap_CKD51

# no graphic environment
DISPLAY=""
export DISPLAY

clustcmd -s -n 17 32 'ps aux | grep sergio' | more
clustcmd -s -n 17 32 killall -9 matlab
clustcmd -s -n 17 32 killall -9 /usr/local/matlabR14/bin/glnx86/MATLAB
clustcmd -s -n 17 32 killall -9 klayers_airs
clustcmd -s -n 17 32 killall -9 sarta_dec05_iceaggr_waterdrop_volz_1_Mar08
clustcmd -s -n 17 32 killall -9 loop_nodes_tropical_dcc_library.sc

#ssh -x airs1  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs2  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs3  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs4  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs5  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs6  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs7  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs8  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs9  "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs10 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs11 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs12 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs13 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs14 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs15 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"
#ssh -x airs16 "killall /usr/local/matlabR14/bin/glnx86/MATLAB"

ssh -x airs17 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs18 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs19 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs20 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs21 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs22 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs23 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs24 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs25 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs26 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs27 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs28 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs29 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs30 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs31 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"
ssh -x airs32 "killall /usr/local/matlabR2006b/bin/glnx86/MATLAB"

ssh -x airs17 "killall do_lte_profiles.sc"
ssh -x airs18 "killall do_lte_profiles.sc"
ssh -x airs19 "killall do_lte_profiles.sc"
ssh -x airs20 "killall do_lte_profiles.sc"
ssh -x airs21 "killall do_lte_profiles.sc"
ssh -x airs22 "killall do_lte_profiles.sc"
ssh -x airs23 "killall do_lte_profiles.sc"
ssh -x airs24 "killall do_lte_profiles.sc"
ssh -x airs25 "killall do_lte_profiles.sc"
ssh -x airs26 "killall do_lte_profiles.sc"
ssh -x airs27 "killall do_lte_profiles.sc"
ssh -x airs28 "killall do_lte_profiles.sc"
ssh -x airs29 "killall do_lte_profiles.sc"
ssh -x airs30 "killall do_lte_profiles.sc"
ssh -x airs31 "killall do_lte_profiles.sc"
ssh -x airs32 "killall do_lte_profiles.sc"
