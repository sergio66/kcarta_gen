xstartup
for ii = 1 : 49
  fname = ['JUNK/rad.dat' num2str(ii)];
  [d,w] = readkcstd(fname);
  dnew(ii,:) = d;

  fname = ['JUNK_OLDCKD6/rad.dat' num2str(ii)];
  [d,w] = readkcstd(fname);
  dold(ii,:) = d;
end

told = rad2bt(w,dold');
tnew = rad2bt(w,dnew');
plot(w,mean(told'-tnew'),w,std(told'-tnew')); grid

ylabel('ckd6 - ckd25 = umbc - mtckd25')
hl = legend('bias','std');