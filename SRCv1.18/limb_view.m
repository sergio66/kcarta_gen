H = 705;
R = 6400;

theta = 0 : 1 : 89;
theta = theta - 0.000001;
theta = theta * pi/180;

a = 1 + 1./(tan(theta) .^ 2);
b = -2*(R+H)./tan(theta);
c = (R+H)^2 - R^2;

plus = -b + sqrt( b.*b - 4 * a .* c);  plus = plus./(2*a);
minus = -b - sqrt( b.*b - 4 * a .* c); minus = minus./(2*a);

oo = find(abs(imag(plus)) > eps | real(plus) < 0); plus(oo) = NaN;
oo = find(abs(imag(minus)) > eps | real(minus) < 0); minus(oo) = NaN;

figure(1);
  plot(theta*180/pi,plus,'bo-',theta*180/pi,minus,'rx-')
grid

thetamax = (H+2*R)*H/((R+H)^2);
thetamax = acos(sqrt(thetamax))*180/pi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma_min = R/(R+H); gamma_min = asin(gamma_min)*180/pi

figure(2)
x = 0 : 100 : 8000;
theta = 0.001; theta = theta * pi/180; m = -1/tan(theta); line01 = m * x + (R+H);
theta = 30; theta = theta * pi/180; m = -1/tan(theta); line30 = m * x + (R+H);
theta = 50; theta = theta * pi/180; m = -1/tan(theta); line50 = m * x + (R+H);
theta = 64; theta = theta * pi/180; m = -1/tan(theta); line64 = m * x + (R+H);
theta = 80; theta = theta * pi/180; m = -1/tan(theta); line80 = m * x + (R+H);
plot(x,[line01; line30; line50; line64; line80]);
hold on
plot(x,line64,'c','linewidth',4);
hold off
viscircles([0 0],6400);
axis equal
axis([-10000 +10000 -10000 +10000]); grid on
axis([-10000 +10000 -00000 +10000]); grid on