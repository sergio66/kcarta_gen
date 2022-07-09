set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'need to read in %5i profiles \n',length(p.stemp));

rad = zeros(89*10000,length(p.stemp));

iDorC = input('Enter DISORT (-1) or Chou/PCLASM (+1) : ');

for ii = 1 : length(p.stemp)
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  if iDorC == -1
    fname = ['JUNK/rad.dat' num2str(ii) '_0605'];
  else
    fname = ['JUNK/rad.dat' num2str(ii)];
  end
  [d,w] = readkcstd(fname);
  if iDorC == +1
    [mm,nn] = size(d);
    d = d(:,nn);
  end
  rad(:,ii) = d;
end
fprintf(1,'\n');
disp('now save w rad into the filename you need ...')
%{
eg 
save -v7.3 /home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/Prof1/Chou0.00/kcartarad.mat w rad
radsing = single(rad);
save -v7.3 /home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/Prof1/Chou0.00/kcartaradsing.mat w radsing
%}

i1231 = find(w >= 1231,1);
plot(p.stemp,rad2bt(1231,rad(i1231,:))); ylabel('BT1231'); xlabel('stemp');
plot(rad2bt(1231,rad(i1231,:))); ylabel('BT1231'); 
semilogx(p.cngwat,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('cngwat');
semilogx(p.cprtop,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('cprtop');
plot(p.ctype,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('ctype');
