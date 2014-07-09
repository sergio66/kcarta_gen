#!/bin/sh
#
# this is a simple demo, normally this would be run with a loop
### if you have runaway processes, do
### killall bkcarta.x.CHINEW; killall kcwrap2; killall kcwrap_CKD51

# no graphic environment
DISPLAY=""
export DISPLAY

clustcmd -s -n 1 16 'ps aux | grep sergio' | more
clustcmd -s -n 1 16 killall -9 matlab
clustcmd -s -n 1 16 killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB
clustcmd -s -n 1 16 killall -9 /usr/local/matlabR14.2/bin/glnx86/matlab_helper
#clustcmd -s -n 1 16 killall -9 klayers

ssh -Y airs1  "killall -9 bkcartaNLTE120.x"
ssh -Y airs2  "killall -9 bkcartaNLTE120.x"
ssh -Y airs3  "killall -9 bkcartaNLTE120.x"
ssh -Y airs4  "killall -9 bkcartaNLTE120.x"
ssh -Y airs5  "killall -9 bkcartaNLTE120.x"
ssh -Y airs6  "killall -9 bkcartaNLTE120.x"
ssh -Y airs7  "killall -9 bkcartaNLTE120.x"
ssh -Y airs8  "killall -9 bkcartaNLTE120.x"
ssh -Y airs9  "killall -9 bkcartaNLTE120.x"
ssh -Y airs10 "killall -9 bkcartaNLTE120.x"
ssh -Y airs11 "killall -9 bkcartaNLTE120.x"
ssh -Y airs12 "killall -9 bkcartaNLTE120.x"
ssh -Y airs13 "killall -9 bkcartaNLTE120.x"
ssh -Y airs14 "killall -9 bkcartaNLTE120.x"
ssh -Y airs15 "killall -9 bkcartaNLTE120.x"
ssh -Y airs16 "killall -9 bkcartaNLTE120.x"

ssh -Y airs1  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs2  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs3  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs4  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs5  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs6  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs7  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs8  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs9  "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs10 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs11 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs12 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs13 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs14 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs15 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"
ssh -Y airs16 "killall -9 /usr/local/matlabR14.2/bin/glnx86/MATLAB"

