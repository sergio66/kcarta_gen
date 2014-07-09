path(path,'/home/sergio/SPECTRA/FORTRANLINUX');

clf;
T = [150 175 200 225 250 275 300 325 350]; T=T';

%% do the +ve freqs 
freq = 0.0;
w_tot = 0.001;
tau2 = 0.005;                   %%%vmr = 363 ppmv??
tau2 = 4.769324235354260E-003;  %%%vmr = 370 ppmv which is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = 0 : 1: 250;

chip150 = birnbaumWORKS(f,freq,w_tot,150,tau2);
chip175 = birnbaumWORKS(f,freq,w_tot,175,tau2);
chip200 = birnbaumWORKS(f,freq,w_tot,200,tau2);
chip225 = birnbaumWORKS(f,freq,w_tot,225,tau2);
chip250 = birnbaumWORKS(f,freq,w_tot,250,tau2);
chip275 = birnbaumWORKS(f,freq,w_tot,275,tau2);
chip300 = birnbaumWORKS(f,freq,w_tot,300,tau2);
chip325 = birnbaumWORKS(f,freq,w_tot,325,tau2);
chip350 = birnbaumWORKS(f,freq,w_tot,350,tau2);

chip=[chip150;chip175;chip200;chip225;chip250;chip275;chip300;chip325;chip350];
plot(f,chip);
axis tight

yp150 = polyfit(f,chip150,3);
yp175 = polyfit(f,chip175,3);
yp200 = polyfit(f,chip200,3);
yp225 = polyfit(f,chip225,3);
yp250 = polyfit(f,chip250,3);
yp275 = polyfit(f,chip275,3);
yp300 = polyfit(f,chip300,3);
yp325 = polyfit(f,chip325,3);
yp350 = polyfit(f,chip350,3);
yp = [yp150;yp175;yp200;yp225;yp250;yp275;yp300;yp325;yp350];

semilogy(T,abs(yp));
zp1 = polyfit(T,yp(:,1),3);
zp2 = polyfit(T,yp(:,2),3);
zp3 = polyfit(T,yp(:,3),3);
zp4 = polyfit(T,yp(:,4),3);
ZP = [zp1; zp2; zp3; zp4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = -250 : 1: 0;

chim150 = birnbaumWORKS(f,freq,w_tot,150,tau2);
chim175 = birnbaumWORKS(f,freq,w_tot,175,tau2);
chim200 = birnbaumWORKS(f,freq,w_tot,200,tau2);
chim225 = birnbaumWORKS(f,freq,w_tot,225,tau2);
chim250 = birnbaumWORKS(f,freq,w_tot,250,tau2);
chim275 = birnbaumWORKS(f,freq,w_tot,275,tau2);
chim300 = birnbaumWORKS(f,freq,w_tot,300,tau2);
chim325 = birnbaumWORKS(f,freq,w_tot,325,tau2);
chim350 = birnbaumWORKS(f,freq,w_tot,350,tau2);

chim=[chim150;chim175;chim200;chim225;chim250;chim275;chim300;chim325;chim350];
plot(f,chim);
axis tight

ym150 = polyfit(f,chim150,3);
ym175 = polyfit(f,chim175,3);
ym200 = polyfit(f,chim200,3);
ym225 = polyfit(f,chim225,3);
ym250 = polyfit(f,chim250,3);
ym275 = polyfit(f,chim275,3);
ym300 = polyfit(f,chim300,3);
ym325 = polyfit(f,chim325,3);
ym350 = polyfit(f,chim350,3);
ym = [ym150;ym175;ym200;ym225;ym250;ym275;ym300;ym325;ym350];

semilogy(T,abs(ym));
zm1 = polyfit(T,ym(:,1),3);
zm2 = polyfit(T,ym(:,2),3);
zm3 = polyfit(T,ym(:,3),3);
zm4 = polyfit(T,ym(:,4),3);
ZM = [zm1; zm2; zm3; zm4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% this is to test the coeffs %%%%%%%%%%%%%%

f = -250:1:250;
iip = find(f >= 0); fp = f(iip);
iim = find(f <  0); fm = f(iim);

T = 240;
zm = polyval(zm1,T)*(fm.^3) + polyval(zm2,T)*(fm.^2) + ...
    polyval(zm3,T)*(fm) + polyval(zm4,T);
zp = polyval(zp1,T)*(fp.^3) + polyval(zp2,T)*(fp.^2) + ...
     polyval(zp3,T)*(fp) + polyval(zp4,T);
zfit = [zm zp];

zcorrect = birnbaumWORKS(f,freq,w_tot,T,tau2);

subplot(211); plot(f,zfit,f,zcorrect,'r')
subplot(212); plot(f,(zfit-zcorrect)./zcorrect);


