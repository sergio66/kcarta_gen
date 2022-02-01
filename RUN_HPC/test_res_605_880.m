addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
clear all
[rv,wv] = readkcstd('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/vhres.dat');
[rh,wh] = readkcstd('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/hres.dat');
[ru,wu] = readkcstd('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/usualres.dat');

plot(wu,rad2bt(wu,ru),'b.-',wh,rad2bt(wh,rh),'gx-',wv,rad2bt(wv,rv),'r+-')
  disp('ret to continue'); pause
xlim([700 750])
  disp('ret to continue'); pause
xlim([660 670])
  disp('ret to continue'); pause
xlim([667 668])
  disp('ret to continue'); pause
xlim([667.3 667.4])
  disp('ret to continue'); pause

[fu,qu] = convolve_airs(wu,ru,1:500);
[fh,qh] = convolve_airs(wh,rh,1:500);
[fv,qv] = convolve_airs(wv,rv,1:500);

figure(1); plot(fh,rad2bt(fh,qh),'r.-',fh,rad2bt(fv,qv),'gx-',fu,rad2bt(fu,qu),'b+-')
figure(2); plot(fh,rad2bt(fh,qh)-rad2bt(fu,qu),'g',fh,rad2bt(fh,qh)-rad2bt(fu,qu),'r'); 
    axis([640 840 -0.5 +0.5]); title('highres-usualres'); plotaxis2;
