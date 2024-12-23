%% fit ax^2 + bx + c to 
%%   p(001) = 7.000 mb
%%   p(041) = 3.500 mb
%%   p(101) = 0.0005 mb
%% and when you evaluate the polynomial between 1--101, make sure all pressures are positive!!!
%% may need to adjust second element to achieve this

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
  grid
figure(3)
  semilogy(mars(:,3),mars(:,2),'bo-',tempMars,plevsMars,'rx-'); set(gca,'ydir','reverse');
  line([100 300],[pGnd pGnd],'color','k')
  line([100 300],[0.0001 0.0001],'color','k')
  grid

disp('ret to continue');
pause

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

