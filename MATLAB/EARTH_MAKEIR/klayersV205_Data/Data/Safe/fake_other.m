% Generate a complete gas profile from just one lower trop and one
% mid-strat mixing ratios
%

% Created: 09 July 2010, Scott Hannon - based on fake_cfc.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gasid = input('enter gas ID number: ');
mrlowtrop = input('enter lower trop mixing ratio (ppm) : ');
mrmidstrat = input('enter mid-strat mixing ratio (ppm) : ');

p=[1.013E+03, 8.988E+02, 7.950E+02, 7.012E+02, 6.166E+02, ...
   5.405E+02, 4.722E+02, 4.111E+02, 3.565E+02, 3.080E+02, ...
   2.650E+02, 2.270E+02, 1.940E+02, 1.658E+02, 1.417E+02, ...
   1.211E+02, 1.035E+02, 8.850E+01, 7.565E+01, 6.467E+01, ...
   5.529E+01, 4.729E+01, 4.047E+01, 3.467E+01, 2.972E+01, ...
   2.549E+01, 1.743E+01, 1.197E+01, 8.010E+00, 5.746E+00, ...
   4.150E+00, 2.871E+00, 2.060E+00, 1.491E+00, 1.090E+00, ...
   7.978E-01, 4.250E-01, 2.190E-01, 1.090E-01, 5.220E-02, ...
   2.400E-02, 1.050E-02, 4.460E-03, 1.840E-03, 7.600E-04, ...
   3.200E-04, 1.450E-04, 7.100E-05, 4.010E-05, 2.540E-05];
%
mr = zeros(1,50);


% constant from surface to 300 mb
ind = find(p >= 300);
mr(ind) = mrlowtrop;
mr300 = mrlowtrop;

% constant from 30 mb to TOA
ind = find(p <= 30);
mr(ind) = mrmidstrat;


% inverse drop off from 300 to 30 mb
ind = find(p > 30 & p < 300);
mr(ind) = (1./p(ind) - 1./30)*(mrlowtrop - mrmidstrat)/(1./300 - 1./30) + ...
   mrmidstrat;


% Write formated data to screen
fprintf(1,' \n',[]);
fprintf(1,'%3i\n',gasid);
for ii=1:9
   jj = (ii-1)*5 + (1:5);
   fprintf(1,'%10.3e, %10.3e, %10.3e, %10.3e, %10.3e,\n',mr(jj));
end
ii=10;
jj = (ii-1)*5 + (1:5);
fprintf(1,'%10.3e, %10.3e, %10.3e, %10.3e, %10.3e\n',mr(jj));
fprintf(1,' \n',[]);


%%% uncomment for testing
plot(mr,p,'.-'),grid
%%%

%%% end oif program %%%
