clear all

addpath /home/sergio/HITRAN2UMBCLBL/FORTRAN/for2mat/

fvhigh = '/asl/data/kcarta/H2016.ieee-le/IR605/etc.ieee-le/lblrtm0.0002/';           %% SHOULD BE 385 ppm
fhigh  = '/asl/data/kcarta/H2016.ieee-le/IR605/etc.ieee-le/lblrtm0.0005/';           %% SHOULD BE 385 ppm
flow   = '/asl/rta/kcarta/H2016.ieee-le/IR605/lblrtm12.8/etc.ieee-le/CO2_400ppmv/';  %% SHOULD BE 400 ppm

[vgid,vfr,vkcomp,vB] = for2mat_kcomp_reader([fvhigh '/y791_g2.dat'],0.0002*10000);
[hgid,hfr,hkcomp,hB] = for2mat_kcomp_reader([fhigh  '/x790_g2.dat'],0.0005*10000);
[lgid,lfr,lkcomp,lB] = for2mat_kcomp_reader([flow   '/r780_g2.dat']);

[vgid,vfr,vkcomp,vB] = for2mat_kcomp_reader([fvhigh '/y719_g2.dat'],0.0002*10000);
[hgid,hfr,hkcomp,hB] = for2mat_kcomp_reader([fhigh  '/x720_g2.dat'],0.0005*10000);
[lgid,lfr,lkcomp,lB] = for2mat_kcomp_reader([flow   '/r705_g2.dat']);

%% see eg /home/sergio/MATLABCODE/KCMIX2/PACKAGE_UPnDOWNLOOK_2014_NLTE/opticaldepths.m
vk = vB * squeeze(vkcomp(:,:,6)); vk = vk.^4;
hk = hB * squeeze(hkcomp(:,:,6)); hk = hk.^4;
lk = lB * squeeze(lkcomp(:,:,6)); lk = lk.^4;

semilogy(lfr,sum(lk,2),'b.-',hfr,sum(hk,2),'gx-',vfr,sum(vk,2),'r')
semilogy(lfr,sum(lk,2),'b.-',hfr,sum(hk,2)*400/385,'g',vfr,sum(vk,2)*400/385,'r') %% remember kh is OD for 385 ppm, so make them equal
