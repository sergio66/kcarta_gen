%[d,w] = readkcstd('l2s_kc115.dat');
[d,w] = readkcstd('l2s_kc116_72gases.dat');

iaGasID_main = [[1:29] [31 32 33 36 37 38 39 40 42]];
iaGasID_xsec = [51:81];

iaGasID = [iaGasID_main iaGasID_xsec 101 102 103];

whos

for ii = 1 : length(iaGasID)
  figure(1); plot(w,exp(-d(:,ii))); axis([605 2830 0 1]); grid
  title(num2str(iaGasID(ii)));

  boo = find(w >= 800 & w <= 1300);
  figure(2); plot(w,exp(-d(:,ii))); axis([800 1300 min(exp(-d(boo,ii)))-eps 1]); grid
  title(num2str(iaGasID(ii)));

  pause(1)
end

ii = input('Enter gas to examine??? ');
while ii > 0
  [Y,iA,iB] = intersect(iaGasID,ii);
  if length(Y) == 1
    figure(1); plot(w,exp(-d(:,iA))); axis([605 2830 0 1]); grid
    title(num2str(iaGasID(iA)));

    boo = find(w >= 800 & w <= 1300);
    figure(2); plot(w,exp(-d(:,iA))); axis([800 1300 min(exp(-d(boo,iA)))-eps 1]); grid
    title(num2str(iaGasID(iA)));
  else
    disp('oops cannot find that gas in this list ...')
  end
  ii = input('Enter gas to examine??? ');
end
