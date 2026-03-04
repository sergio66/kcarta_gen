clear all

h16 = load('../L2SComparisons/l2s_kc122_H16_605_2830.mat');
h20 = load('../L2SComparisons/l2s_kc122_H20_605_2830.mat');
h24 = load('../L2SComparisons/l2s_kc122_H24_605_2830.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('convolving')
iY = input('convolve before seeing results (-1/+1 default) : ');
if length(iY) == 0
  iY = +1;
end

addpath /home/sergio/git/matlabcode
if iY > 0
  [h16.w,h16.d] = quickconvolve(h16.w,h16.d,0.25,0.25);
  [h20.w,h20.d] = quickconvolve(h20.w,h20.d,0.25,0.25);
  [h24.w,h24.d] = quickconvolve(h24.w,h24.d,0.25,0.25);
end
w = h16.w;

comment = 'see /home/sergio/KCARTA/UTILITY/compare_l2s_H16_H20.m';
gidlist.h16 = h16.iaGasID;
gidlist.h20 = h20.iaGasID;
gidlist.h24 = h24.iaGasID;
save H16_20_24_gidlist.mat gidlist comment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('looking for gases in H16,H20 not in H24')
[YI,i20,I16] = intersect(h20.iaGasID,h16.iaGasID);
[YJ,J20,J24] = intersect(h20.iaGasID,h24.iaGasID);
[YK,K16,K24] = intersect(h16.iaGasID,h24.iaGasID);

if length(h24.iaGasID) ~= length(h20.iaGasID)
  missing = setdiff(h20.iaGasID,h24.iaGasID);
  fprintf(1,'hmm, one less gas in H2024???? = %2i but looks like the ODs are zero anyway???? \n',missing)
  junk = find(h20.iaGasID == missing);
  plot(h16.w,h16.d(:,junk),'b.-',h20.w,h20.d(:,junk),'r'); legend('H16','H20')

  disp('    this file /home/sergio/KCARTA/SCRIPTS/MAKE_COMP_HTXY_PARAM_SC/comp_IRdatabase_H2024.sc')
  disp('        says quit double counting gases 30/81 35/61 41/80 42/54')
  disp('         ---->>>> deleted gases 30 35 41 42 <<<<<<<<<<<<<---------')

  title(['H24 is missing gasID ' num2str(missing)])
  disp('ret to continue'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('off we go .....')
figure(4); clf
figure(5); clf

od16 = zeros(size(h24.w));
od20 = zeros(size(h24.w));
od24 = zeros(size(h24.w));

ia24 = 1 : length(h24.iaGasID);
ia24 = [1 length(h24.iaGasID)];
for ii24 = 1 : length(ia24)
  i24 = ia24(ii24);

  i16 = find(h16.iaGasID == h24.iaGasID(i24));
  i20 = find(h20.iaGasID == h24.iaGasID(i24));  
  
  boo = find(h20.d(:,i20) > eps & h16.d(:,i16) > eps & h24.d(:,i24) > eps);
  boo1 = find(h16.d(:,i16) > eps);
  boo2 = find(h20.d(:,i20) > eps);
  boo3 = find(h24.d(:,i24) > eps);  

  figure(1); semilogy(w(boo),h20.d(boo,i20)); title(num2str(h20.iaGasID(i20)));
  if iY > 0
    figure(1); semilogy(w(boo),h16.d(boo,i16),'b',w(boo),h20.d(boo,i20),'g',w(boo),h24.d(boo,i24),'r'); title(num2str(h20.iaGasID(i20)));
    legend('h16','h20','h24')
  end
  
  figure(2); plot(w(boo),h16.d(boo,i16)./h20.d(boo,i20),'b.-',w(boo),h24.d(boo,i24)./h20.d(boo,i20),'r'); 
    title(num2str(h20.iaGasID(i20))); ylim([0 2]); plotaxis2; 
    legend('H2016/H2020','H2024/H2020','location','best');

  figure(3); dr = 0.75 : 0.0001 : 1.25; 
    semilogy(dr,histc(h16.d(boo,i16)./h20.d(boo,i20),dr)/length(boo),'b.-',dr,histc(h24.d(boo,i24)./h20.d(boo,i20),dr)/length(boo),'r')
    title(num2str(h20.iaGasID(i20))); plotaxis2;
    legend('hist(H2020/H2016)','hist(H2024/H2016)','location','best');

  figure(4);
  if i24 <= 9
    semilogy(w,h24.d(:,i24)); hold on;
    xlim([605 1650])
    ylim([1e-10 1e5])
  end

  od16 = od16 + h16.d(:,i16);
  od20 = od20 + h20.d(:,i20);
  od24 = od24 + h24.d(:,i24);
  
  tr16 = exp(-od16);
  tr20 = exp(-od20);
  tr24 = exp(-od24);
  
  bt16 = rad2bt(w,ttorad(w,290).*exp(-od16));
  bt20 = rad2bt(w,ttorad(w,290).*exp(-od20));
  bt24 = rad2bt(w,ttorad(w,290).*exp(-od24));

  figure(5);
    plot(w,tr16,'g',w,tr20,'b',w,tr24,'r')
    xlim([1100 1150])    
  figure(6);
    plot(w,bt16,'g',w,bt20,'b',w,bt24,'r')
    xlim([1100 1150])    
  end


  junk = [length(boo1)/10000 length(boo2)/10000 mean(h16.d(boo,i16)./h20.d(boo,i20)) std(h16.d(boo,i16)./h20.d(boo,i20))];
  %fprintf(1,'gid %2i numchunks H2016:H2020 = %5.2f %5.2f ratio H2016/H2022 = %8.6f +/- %8.6f \n',h20.iaGasID(i20),junk);
  junk(1:2) = ceil(junk(1:2));
  fprintf(1,'gid %2i numchunks H2016:H2020 = %2i %2i ratio H2016/H2020 = %8.6f +/- %8.6f \n',h20.iaGasID(i20),junk);

  junk = [length(boo1)/10000 length(boo2)/10000 mean(h24.d(boo,i24)./h20.d(boo,i20)) std(h24.d(boo,i24)./h20.d(boo,i20))];
  %fprintf(1,'gid %2i numchunks H2024:H2020 = %5.2f %5.2f ratio H2024/H2022 = %8.6f +/- %8.6f \n',h20.iaGasID(i20),junk);
  junk(1:2) = ceil(junk(1:2));
  fprintf(1,'gid %2i numchunks H2024:H2020 = %2i %2i ratio H2024/H2020 = %8.6f +/- %8.6f \n',h20.iaGasID(i20),junk);

  disp('ret to continue'); pause
  
  pause(0.5); 
end
