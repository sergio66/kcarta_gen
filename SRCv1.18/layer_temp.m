function [yusual,yusualfine,ylinear,yexp,yLBL,tau] = layer_temp(nu,Tin,TBot,TTop,tau0,theta)

% [yusual,yusualfine,ylinear,yexp,yLBL,tau] = layer_temp(nu,Tin,TBot,TTop,tau0,theta)
%
% input
%   nu        = wavenumber 
%   Tin       = what the equivalent BT incident at bottom of layer corresponds to
%   Tbot,Ttop = temperatures at bottom/top of the layer
%   tau0      = total OD of layer
%   theta     = view angle (degrees)
% output
%   yusual     = what is in SARTA (using layer temp = (Ttop + Tbot)/2)
%   yusualfine = divvying up the layer into 100 chunks, doing SARTA RT on each little chunk
%   ylinear    = RT using simple linear-in-tau variation
%   yexp       = RT using exponential-in-tau variation
%   yLBL       = RT using LBLRTM fancy linear in tau variation
%   tau        = tau0 is divided into 100 chunks, spaced tau = 0 : dt : tau0
%
% example : dz = 0.25; [a,b,c,d,e,tau] = layer_temp(805,305,300,300-lapse*dz,50,22);
%           dz = 0.25; [a,b,c,d,e,tau] = layer_temp(805,302,300,300-lapse*dz,0.5,22);

addpath /asl/matlib/aslutil

dt = tau0/101;
tau = 0 : dt : tau0;

mu = cos(theta*pi/180);

rIn  = ttorad(nu,Tin);
rBot = ttorad(nu,TBot);
rTop = ttorad(nu,TTop);

rMid = ttorad(nu,(TTop+TBot)/2);
yusual = rIn.*exp(-tau/mu) + rMid.*(1-exp(-tau/mu));

tz = interp1([0 tau0],[TBot TTop],tau);
%plot(tau,tz); error('oo')
yusualfine(1) = rIn;
for ii = 2 : length(tz)
  yusualfine(ii) = yusualfine(ii-1) * exp(-dt/mu) + ttorad(nu,tz(ii))*(1-exp(-dt/mu));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yusual0 = rIn.*exp(-tau/mu) + rBot.*(1-exp(-tau/mu));

ylinear = zeros(size(yusual));
gamma = (rTop-rBot)/tau0;
ylinear = gamma*mu*(tau/mu-1) + gamma*mu*exp(-tau/mu);
ylinear = ylinear + yusual0;

yexp    = zeros(size(yusual));
beta    = 1/tau0*log(rTop/rBot);
yexp    = rIn.*exp(-tau/mu) + rBot/(1+beta*mu).*(exp(beta*tau)-exp(-tau/mu));

yLBL    = zeros(size(yusual));
f       = 1 - 2*(1./tau0 - exp(-tau0)./(1-exp(-tau0)));
Beff    = rMid + (rTop-rMid).*f;
yLBL    = rIn.*exp(-tau/mu) + Beff .*(1-exp(-tau/mu));

plot(tau,rad2bt(nu,yusual),'bo-',tau,rad2bt(nu,yusualfine),'cd-',...
     tau,rad2bt(nu,ylinear),'g',tau,rad2bt(nu,yexp),'rs-',...
     tau,rad2bt(nu,yLBL),'k');
hl = legend('SARTA/kCARTA','SARTA/kCARTA in 100','linear','exp','LBLRTM','location','northeast');
set(hl,'fontsize',10); grid
xlabel('tau'); ylabel('equaivalent BT (K)')

tusual = rad2bt(nu,yusual); tusual = tusual(end);
tusualfine = rad2bt(nu,yusualfine); tusualfine = tusualfine(end);
tlinear = rad2bt(nu,ylinear); tlinear = tlinear(end);
texp = rad2bt(nu,yexp); texp = texp(end);
tLBL = rad2bt(nu,yLBL); tLBL = tLBL(end);

[tusual tusualfine tlinear texp tLBL]