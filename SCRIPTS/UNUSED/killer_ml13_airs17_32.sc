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
clustcmd -s -n 17 32 killall -9 klayers

#ssh -x airs1  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs2  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs3  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs4  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs5  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs6  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs7  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs8  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs9  "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs10 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs11 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs12 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs13 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs14 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs15 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
#ssh -x airs16 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs17 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs18 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs19 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs20 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs21 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs22 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs23 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs24 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs25 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs26 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs27 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs28 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs29 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs30 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs31 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
ssh -x airs32 "killall /usr/local/matlabR13.1/bin/glnx86/matlab"
