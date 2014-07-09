%function raaRayleigh = rayleigh(raFreq,raVT1)

profile = load('/home/sergio/SPECTRA/IPFILES/std_co2');
raVT1 = profile(:,4);
raFreq = 0.25:0.01:1.0;
raFreq = 10000./raFreq;

fit1 =  9.38076e+18;
fit2 = -1.08426e+09;

pzero = 1013.25;
tzero = 273.15;

% rayleigh scattering coefficient (1/km) 

raPressLevels = load('/home/sergio/MATLABCODE/airslevels.dat');

for iLay = 1:100
  p_nz(iLay) = raPressLevels(iLay) - raPressLevels(iLay+1);
  p_nz(iLay) = p_nz(iLay)/log(raPressLevels(iLay)/raPressLevels(iLay+1));
  t_nz(iLay) = raVT1(iLay);
  dz(iLay)   = abs(p2h(raPressLevels(iLay)) - p2h(raPressLevels(iLay+1)))/1000;
  end

v = raFreq;
sig = v.^4./(fit1+fit2*v.^2);
for iLay = 1:100
  % assume 5km scale ht.
  raaRayleigh(:,iLay)=sig*(p_nz(iLay)/pzero)/(t_nz(iLay)/tzero)*dz(iLay);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%see ../DOC/PDF/RayleighScatteringForAtmos.pdf
%Wavelength Atmospheric models
%(nm) 15'N, mean 450N, winter 450N, summer 60'W, winter 60'N, summer
Frolich_Shaw = [...
260.0 2.18505 2.18313 2.18522 2.18196 2.17904;
280.0 1.57963 1.57832 1.57979 1.57748 1.57533;
300.0 1.17273 1.17180 1.17286 1.17118 1.16957;
320.0 0.89022 0.88953 0.89032 0.88907 0.88783;
340.0 0.68867 0.68816 0.68876 0.68780 0.68683;
360.0 0.54153 0.54114 0.54161 0.54086 0.54010;
380.0 0.43196 0.43166 0.43202 0.43143 0.43082;
400.0 0.34894 0.34869 0.34899 0.34851 0.34802;
420.0 0.28505 0.28486 0.28509 0.28471 0.28430;
440.0 0.23522 0.23506 0.23525 0.23494 0.23460;
460.0 0.19587 0.19574 0.19590 0.19564 0.19536;
480.0 0.16445 0.16434 0.16448 0.16426 0.16402;
500.0 0.13911 0.13902 0.13913 0.13895 0.13875;
550.0 0.09424 0.09418 0.09426 0.09413 0.09400;
600.0 0.06613 0.06609 0.06614 0.06606 0.06596;
650.0 0.04778 0.04775 0.04779 0.04773 0.04766;
700.0 0.03539 0.03537 0.03540 0.03535 0.03530;
750.0 0.02678 0.02676 0.02678 0.02675 0.02671;
800.0 0.02063 0.02062 0.02064 0.02061 0.02058;
850.0 0.01616 0.01615 0.01616 0.01614 0.01611;
900.0 0.01283 0.01282 0.01283 0.01282 0.01280;
1000.0 0.00840 0.00839 0.00840 0.00839 0.00838;
1100.0 0.00572 0.00572 0.00572 0.00572 0.00571;
1200.0 0.00404 0.00403 0.00404 0.00403 0.00402;
1300.0 0.00293 0.00292 0.00293 0.00292 0.00292;
1400.0 0.00217 0.00217 0.00217 0.00217 0.00217;
1500.0 0.00165 0.00165 0.00165 0.00165 0.00164;
];
%Elevation: 0.00 km.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot(10000./raFreq,sum(raaRayleigh(:,4:100)'),'ko-','linewidth',2); hold on
plot(Frolich_Shaw(:,1)/1000,Frolich_Shaw(:,2:6)); hold off

figure(2)
semilogy(10000./raFreq,sum(raaRayleigh(:,4:100)'),'ko-','linewidth',2); hold on
semilogy(Frolich_Shaw(:,1)/1000,Frolich_Shaw(:,2:6)); hold off
