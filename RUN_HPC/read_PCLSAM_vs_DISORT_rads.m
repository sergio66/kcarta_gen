set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'need to read in %5i profiles \n',length(p.stemp));

rad = zeros(89*10000,length(p.stemp));

iDorP = input('Enter DISORT (-1) or Chou/PCLSAM (+1, default) : ');
if length(iDorP) == 0
  iDorP = +1;
end

badfile = [];;
for ii = 1 : length(p.stemp)
  if mod(ii,1000) == 0
    fprintf(1,'+')
  elseif mod(ii,100) == 0
    fprintf(1,'.')
  end

  if iDorP == -1
    fname = ['JUNK/rad.dat' num2str(ii) '_0605'];
  else
    fname = ['JUNK/rad.dat' num2str(ii)];
  end
  try 
    [d,w] = readkcstd(fname);
    if iDorP == +1
      [mm,nn] = size(d);
      d = d(:,nn);
    end
    rad(:,ii) = d;
  catch
    fprintf(1,'problem with file %s \n',fname)
    badfile = [badfile ii];
  end
end
fprintf(1,'\n');
disp('now save w rad into the filename you need ...')
%{
eg 
save -v7.3 /home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/Prof1/Chou0.00/kcartarad.mat w rad
radsing = single(rad);
wsing   = single(w);
save -v7.3 /home/sergio/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/Prof1/Chou0.00/kcartaradsing.mat wsing radsing
%}

i1231 = find(w >= 1231,1);
plot(p.stemp,rad2bt(1231,rad(i1231,:))); ylabel('BT1231'); xlabel('stemp');
plot(rad2bt(1231,rad(i1231,:))); ylabel('BT1231'); 
semilogx(p.cngwat,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('cngwat');
semilogx(p.cprtop,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('cprtop');
plot(p.ctype,rad2bt(1231,rad(i1231,:)),'.'); ylabel('BT1231'); xlabel('ctype');
