%{
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:44 generic_2380_2830_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:44 generic_2280_2380_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:44 generic_2130_2280_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:44 generic_2005_2130_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:43 generic_1655_2005_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:42 generic_1380_1655_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:42 generic_1280_1380_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:41 generic_1105_1280_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:41 generic_0980_1105_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:40 generic_0805_0980_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:40 generic_0667_0670_I_W.mat XXX
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:39 generic_0607_0805_I_W.mat
-rw-rw-r-- 1 sergio pi_strow 12K Jul 24 15:38 generic_605_1655_I_W.mat
%}

%[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_605_1655_I_W.mat',      'generic_605_1655_I_W.bin',      605,2830,iVers=1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_0607_0805_I_W.mat','generic_0605_0805_I_W.bin',0605,0805,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_0805_0980_I_W.mat','generic_0805_0980_I_W.bin',0805,0980,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_0980_1105_I_W.mat','generic_0980_1105_I_W.bin',0980,1105,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_1105_1280_I_W.mat','generic_1105_1280_I_W.bin',1105,1280,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_1280_1380_I_W.mat','generic_1280_1380_I_W.bin',1280,1380,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_1380_1655_I_W.mat','generic_1380_1655_I_W.bin',1380,1655,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_1655_2005_I_W.mat','generic_1655_2005_I_W.bin',1655,2005,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_2005_2130_I_W.mat','generic_2005_2130_I_W.bin',2005,2130,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_2130_2280_I_W.mat','generic_2130_2280_I_W.bin',2130,2280,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_2280_2380_I_W.mat','generic_2280_2380_I_W.bin',2280,2380,1);  %% YAY
[matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_2380_2830_I_W.mat','generic_2380_2830_I_W.bin',2380,2830,1);  %% YAY
