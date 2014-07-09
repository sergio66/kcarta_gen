#!/bin/sh
#
# this is a simple demo, normally this would be run with a loop
### if you have runaway processes, do
### killall bkcarta.x.CHINEW; killall kcwrap2; killall kcwrap_CKD51

# no graphic environment
DISPLAY=""
export DISPLAY

ssh -x airs1  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs2  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs3  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs4  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs5  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs6  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs7  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs8  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs9  "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs10 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs11 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs12 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs13 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs14 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs15 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs16 "cd /airs/s1/sergio; /bin/rm k*;"
ssh -x airs17 "cd /airs/s1/sergio; /bin/rm k*;"
