function [kcartajacdump] = readkc_all_dumpjac(raFreq,iaGasID)
%function [raaRad,raaTau,raaOneMinusTau,raaLay2Gnd,raaGeneral,raaAllDT,raaaAllDQ] = readkc_all_dumpjac(raFreq,iaGasID)

%% compare to the kcmix in ~/git/kcarta/JACDOWN/jac_downlook.m
%% compare to the kcmix in ~/git/kcarta/JACDOWN/jac_downlook.m
%% compare to the kcmix in ~/git/kcarta/JACDOWN/jac_downlook.m

wnum2str = num2str(raFreq(1),'%04d');

strGas  = num2str(iaGasID(1),'%03d');
thestr0 = ['dump_jacobian_' wnum2str];

for ig = 1 : length(iaGasID)
  boo = [thestr0 '_gasID_' num2str(iaGasID(ig),'%03d') '.dat'];
  [junk,w] = readkc_dumpjac_info(boo);
  raaaAllDQ(ig,:,:) = junk;
end

w = 1:10000;
dv = 0.0025;
kcartajacdump.w = raFreq(1)+ (w-1)*dv;

kcartajacdump.w;

kcartajacdump.jacTG          = readkc_dumpjac_info([thestr0 '_temp.dat']);
kcartajacdump.jacQG          = readkc_dumpjac_info([thestr0 '_gasID_' strGas '.dat']);

kcartajacdump.rawjacTG       = readkc_dumpjac_info([thestr0 '_rawD_DT_temp.dat']);
kcartajacdump.rawjacQG       = readkc_dumpjac_info([thestr0 '_rawD_DQ_gasID_' strGas '.dat']);

kcartajacdump.raaLay2Gnd     = readkc_dumpjac_info([thestr0 '_raaLay2Gnd.dat']);
kcartajacdump.raaLay2Sp      = readkc_dumpjac_info([thestr0 '_raaLay2Sp.dat']);

kcartajacdump.raaOneMinusTau = readkc_dumpjac_info([thestr0 '_raaOneMinusTau.dat']);
kcartajacdump.raaTau         = readkc_dumpjac_info([thestr0 '_raaTau.dat']);

kcartajacdump.raaRad         = readkc_dumpjac_info([thestr0 '_raaRad.dat']);
kcartajacdump.raaRadDT       = readkc_dumpjac_info([thestr0 '_raaRadDT.dat']);

kcartajacdump.raaGeneral     = readkc_dumpjac_info([thestr0 '_raaGeneral.dat']);
kcartajacdump.absc           = readkc_dumpjac_info([thestr0 '_raaAbsCoeff.dat']);

[~,~,iNumLayer] = readkc_dumpjac_info([thestr0 '_raaGeneral.dat']);;
kcartajacdump.iNumLayer = iNumLayer;
