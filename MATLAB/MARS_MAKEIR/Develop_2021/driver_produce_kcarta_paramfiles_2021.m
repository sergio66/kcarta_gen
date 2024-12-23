%% fit ax^2 + bx + c to 
%%   p(001) = 7.000 mb
%%   p(041) = 3.500 mb
%%   p(101) = 0.0005 mb
%% and when you evaluate the polynomial between 1--101, make sure all pressures are positive!!!
%% may need to adjust second element to achieve this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = polyfit([1 41 101],[1100 500 0.005],2); pGnd = 1013;  %% EARTH
P = polyfit([1 31 101],[6.75 3.50 0.001],2); pGnd = 6.31;   %% MARS
P = polyfit([1 29 101],[6.75 3.50 0.0001],2); pGnd = 6.31;   %% MARS
Y = polyval(P,1:101);
figure(5); plot(Y,1:101,'o-'); grid; ax = axis; 
  line([pGnd pGnd],[1 101],'color','k'); text(5,80,'Surface');
  line([0.5 0.5],[1 101],'color','r'); text(1,80,'Tropopause');

bad = find(Y < 0);
if length(bad) > 0
  error('found negative pressures!')
end
title('MARS Std Atm')

mars = load('mars.txt');
plevsMars   = Y;
  [max(plevsMars) min(plevsMars)]
heightsMars = interp1(log(mars(:,2)),mars(:,1),log(plevsMars),[],'extrap');
  [max(heightsMars) min(heightsMars)]
tempMars    = interp1(log(mars(:,2)),mars(:,3),log(plevsMars),[],'extrap');

figure(4)
  semilogy(mars(:,1),mars(:,2),'bo-',heightsMars,plevsMars,'rx-'); set(gca,'ydir','reverse');
  line([-5 +100],[pGnd pGnd],'color','k')
  line([-5 +100],[0.0001 0.0001],'color','k')
  xlabel('height(km)');   ylabel('plevs(mb)');
  grid
figure(3)
  semilogy(mars(:,3),mars(:,2),'bo-',tempMars,plevsMars,'rx-'); set(gca,'ydir','reverse');
  line([100 300],[pGnd pGnd],'color','k')
  line([100 300],[0.0001 0.0001],'color','k')
  xlabel('T(K)');   ylabel('plevs(mb)');
  grid

heightsMars2014 = heightsMars;
plevsMars2014 = plevsMars;
tempMars2014 = tempMars;

disp('this was 2014 : ret to continue');
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% this is from mars_sergio_rtp.m

h = 0:250:120000;
p = zeros(size(h));
T = zeros(size(h));

%% T in C, h in meters, p in kPa
oo = find(h <= 7000);
 T(oo) = -31 - 0.000998 * h (oo);       %% -dT/dh = 0.998 K/km
 p(oo) = .699 * exp(-0.00009 * h(oo));

oo = find(h > 7000);
  T(oo) = -23.4 - 0.00222 * h(oo) + 1e-8*h(oo).^2;      %% -dT/dh = 2.2 K/km, I did the quadratic modification 
  p(oo) = .699 * exp(-0.00009 * h(oo));

T = T + 273;
n = p ./ (0.1921 * T);   %% density using eqn of state

p = p*1000/100;  %% KPa to Pa = N/m2 to mb
T = T;           %% in K
n = n/1e6;       %% per cm3

hp = 80; matr = [1 1 1; hp*hp hp 1; 10201 101 1]; yvec = [13 0.3 0.001].^(2/7);  ABC = matr\yvec'     %% for Mars
figure(4); A = ABC(1); B = ABC(2); C = ABC(3); ii = 1:101; P = (A*ii.^2 + B*ii + C).^(7/2); semilogy(ii,P,'o-')

TP = interp1(log(p),T,log(P),[],'extrap'); figure(3); semilogy(T,p,TP,P); set(gca,'ydir','reverse');
HP = interp1(log(p),h,log(P),[],'extrap'); figure(4); semilogy(h,p,HP,P); set(gca,'ydir','reverse');
NP = interp1(log(p),n,log(P),[],'extrap'); figure(2); semilogy(n,p,NP,P); set(gca,'ydir','reverse');

Y = P;
plevsMars   = P;
  [max(plevsMars) min(plevsMars)]
heightsMars = HP/1000;
  [max(heightsMars) min(heightsMars)]
tempMars    = TP;

figure(4)
  semilogy(heightsMars2014,plevsMars2014,'bo-',heightsMars,plevsMars,'rx-'); set(gca,'ydir','reverse');
  line([-5 +100],[pGnd pGnd],'color','k')
  line([-5 +100],[0.0001 0.0001],'color','k')
  xlabel('height(km)');   ylabel('plevs(mb)');
  grid
figure(3)
  semilogy(tempMars2014,plevsMars2014,'bo-',tempMars,plevsMars,'rx-'); set(gca,'ydir','reverse');
  line([100 300],[pGnd pGnd],'color','k')
  line([100 300],[0.0001 0.0001],'color','k')
  xlabel('T(K)');   ylabel('plevs(mb)');
  grid

heightsMars2021 = heightsMars;
plevsMars2021 = plevsMars;
tempMars2021 = tempMars;

disp('this is new 2021 : ret to continue');
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YN = Y(1:end-1)-Y(2:end);
YD = log(Y(1:end-1)./Y(2:end));
Yav = YN ./ YD;

  make_KCARTA_database

  make_kcarta_airslevelheights
  make_kcarta_airslevels
  make_kcarta_airsheights

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% not really interested in UA NLTE, but kCARTA needs these placeholders

%% fit ax^2 + bx + c to 
%%   p(001) = 1e-2 mb
%%   p(081) = 1e-3 mb
%%   p(101) = 1e-51 mb
%% and when you evaluate the polynomial between 1--101, make sure all pressures are positive!!!
%% may need to adjust second element to achieve this
Pu = polyfit([1 81 101],[1e-3 1e-4 1e-6],2); 
Yu = polyval(Pu,1:101);
plevsMarsUpper   = Yu; 
heightsMarsUpper = interp1(log(mars(:,2)),mars(:,1),log(plevsMarsUpper),[],'extrap');;

Yu = plevsMarsUpper;
YN = Yu(1:end-1)-Yu(2:end);
YD = log(Yu(1:end-1)./Yu(2:end));
YavUpper = YN ./ YD;

figure(5)
plot(Yu,1:101)
bad = find(Yu < 0);
if length(bad) > 0
  error('found negative pressures!')
end
title('MARS Upper Atm')
disp('ret to continue');
pause

  make_kcarta_airslevelheights_upper
  make_kcarta_airslevels_upper
  make_kcarta_airsheights_upper

