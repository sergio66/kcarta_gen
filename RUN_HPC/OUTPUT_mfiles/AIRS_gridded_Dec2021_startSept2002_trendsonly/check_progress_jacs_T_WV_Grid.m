lser = ['!ls -lt AllDemJacsAnomaly/T_WV_Grid/*/*/ind*jac.mat | wc -l']; eval(lser)

for ii = 1 : 72
  for jj = 1 : 64
    istr = num2str(ii,'%02d');
    jstr = num2str(jj,'%02d');
    thedir = dir(['AllDemJacsAnomaly/T_WV_Grid/' jstr '/' istr '/ind*jac.mat']);
    thecount(ii,jj) = length(thedir);
  end
end

jett = jet(26); jett(1,:) = 1;
figure(1); pcolor(1:72,1:64,thecount'); shading flat; colorbar; colormap(jett); title('count of jac files made so far'); xlabel('LonBin'); ylabel('LatBin');
figure(2); plot(1:64,sum(thecount),'o-'); xlabel('Latbin');
figure(2); plot(1:64,sum(thecount)/(72*25),'o-'); xlabel('Latbin Done (N/(72*25))'); 
  xticks([0:4:64]); grid on
fprintf(1,'done %8i so far out of %8i or %8.5f percent \n',sum(sum(thecount)),25*64*72,sum(sum(thecount))/(25*64*72)*100)
