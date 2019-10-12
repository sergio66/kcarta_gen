%% see JGR 1992 v97
%% Line by Line Calculations of Atmospheric Fluxes and Cooling Rates' Application to Water Vapor
%% SHEPARDA. CLOUGH, MICHAEL J. IACONO, AND JEAN-LUC MONCET

%% Eqn 9

addpath /home/sergio/KCARTA/MATLAB

junk = input('Enter TL TU : ');
TL = junk(1); TU = junk(2);  

tau = 0.0001 : 0.01 : 10; T = exp(-tau);

v = 651.95;  %% deep in the CO2 15 um
BU = ttorad(v,TU);
BL = ttorad(v,TL);

I0 = BU * (1-T) - (BL-BU)*T + (BL-BU)./tau.*(1-T);        %% eqn 9 EXACT

Bav = (BU+BL)/2;
I1 = BU*(1-T) - 2*(Bav-BU)*T + 2 * (Bav-BU)./tau.*(1-T);  %% eqn 11

I2 = (1-T).*(BU + 2*(Bav-BU)*(1./tau - T./(1-T)));        %% eqn 13,14

a = 0.2;
I3 = (1-T).*(Bav + (a*tau)*BU)./(1+a*tau);                %% eqn 15, PADE USED IN LBLRTM

a = 0.193; b = 0.013;
I4 = (1-T).*(Bav + (a*tau + b*tau.*tau)*BU)./(1+a*tau+b*tau.*tau);   %% eqn 16

subplot(211); plot(tau,[I0; I1; I2; I3; I4]);
subplot(212); plot(tau,I1./I0,tau,I2./I0,tau,I2./I0,tau,I3./I0)

clf
plot(tau,I0,'kx-',tau,I1,'g',tau,I2,'r',tau,I3,'b.-',tau,I4,'c');
hl = legend('I0=exact','I1','I2','I3=pade','I4','location','best');
