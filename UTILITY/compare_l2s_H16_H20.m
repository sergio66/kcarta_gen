clear all
h16 = load('../L2SComparisons/l2s_kc122_H16_605_2830.mat');
h20 = load('../L2SComparisons/l2s_kc122_H20_605_2830.mat');

for ii = 1 : 72; 
  boo = find(h20.d(:,ii) > eps & h16.d(:,ii) > eps);
  boo1 = find(h16.d(:,ii) > eps);
  boo2 = find(h20.d(:,ii) > eps);

  figure(1); semilogy(h20.w(boo),h20.d(boo,ii)); title(num2str(h20.iaGasID(ii)));

  figure(2); plot(h20.w(boo),h20.d(boo,ii)./h16.d(boo,ii)); 
    title(num2str(h20.iaGasID(ii))); ylim([0 2]); plotaxis2; 
    ylabel('H2020/H2016');

  figure(3); dr = 0.75 : 0.0001 : 1.25; 
    semilogy(dr,histc(h20.d(boo,ii)./h16.d(boo,ii),dr)/length(boo)); 
    title(num2str(h20.iaGasID(ii))); plotaxis2;
    ylabel('hist(H2020/H2016)');

  junk = [length(boo1)/10000 length(boo2)/10000 mean(h20.d(boo,ii)./h16.d(boo,ii)) std(h20.d(boo,ii)./h16.d(boo,ii))];
  %fprintf(1,'gid %2i numchunks H2016:H2020 = %5.2f %5.2f ratio H2020/H2016 = %8.6f +/- %8.6f \n',h20.iaGasID(ii),junk);
  junk(1:2) = ceil(junk(1:2));
  fprintf(1,'gid %2i numchunks H2016:H2020 = %2i %2i ratio H2020/H2016 = %8.6f +/- %8.6f \n',h20.iaGasID(ii),junk);

  disp('ret to continue'); pause
  pause(0.5); 
end
