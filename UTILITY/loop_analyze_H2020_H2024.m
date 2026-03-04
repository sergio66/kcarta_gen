load('H16_20_24_gidlist.mat');
glist0 = gidlist.h24;

disp(' ')
disp('if you do    load ../L2SComparisons/l2s_kc122_H24_605_2830.mat')
disp('             semilogy(w,d(:,1),w,d(:,71))')
disp('             you will clearly see that KCARTA can differentiate between G1 and G103')
disp('This not self evident if wantlist = [1 103] and iDiff = +1 (showing the DBT between H24 and H20)')
disp('  But is self evident if wantlist = [1 103] and iDiff = -1 (showing the rawBT of H24 and H20)')
disp(' ')
disp('Therefore g1,g103 in H20 and H24 are essentially identical after convolution')
disp(' ')

xlim1 = [1050 1300];
xlim2 = [2500 2800];
xlim3 = [0650 0850];
xlim4 = [1250 1650];

gg = 0;
ii = length(glist0) + 1;
fname = ['JUNK/individual_prof_convolved_kcartaH2020_H2024_' num2str(gg) '.mat'];
if exist(fname)
  iaDone(ii) = 1;
  iaGasID(ii) = gg;
  loader = ['x = load(''' fname ''');'];
  eval(loader);
  diffBT0 = rad2bt(x.fc,x.qc24)-rad2bt(x.fc,x.qc20);
  allBT0  = rad2bt(x.fc,x.qc24);
  figure(1)
    plot(x.fc,diffBT0); xlim(xlim1)
    title(['ALL GASES ' num2str(gg)]);
    xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
  figure(2)
    plot(x.fc,diffBT0); xlim(xlim2)
    title(['ALL GASES ' num2str(gg)]);
    xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
  figure(3)
    plot(x.fc,diffBT0); xlim(xlim3)
    title(['ALL GASES ' num2str(gg)]);
    xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
  figure(4)
    plot(x.fc,diffBT0); xlim(xlim4)
    title(['ALL GASES ' num2str(gg)]);
    xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
  pause;
end

iaDone  = zeros(1,length(glist0));
iaGasID = nan(1,length(glist0));

%%%%%%%%%%%%%%%%%%%%%%%%%

iaDo = glist0;
wantlist = [1 103];         [Y,iaDo] = intersect(glist0,wantlist);
wantlist = [1 4 56 72 103]; [Y,iaDo] = intersect(glist0,wantlist);   %% [4 56 72] are the gases that show differences

%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : length(iaDo)
  gg = glist0(iaDo(ii));
  fname = ['JUNK/individual_prof_convolved_kcartaH2020_H2024_' num2str(gg) '.mat'];
  if exist(fname)
    iaDone(ii) = 1;
    iaGasID(ii) = gg;
    loader = ['x = load(''' fname ''');'];
    eval(loader);
    
    iDiff = -1;  %% BT24,BT20,all
    iDiff = +1;  %% DBT = BT24-BT20     DEFAULT
    
    if iDiff > 0    
      figure(1);
        plot(x.fc,diffBT0,'k',x.fc,rad2bt(x.fc,x.qc24)-rad2bt(x.fc,x.qc20),'r'); xlim(xlim1)
        title(num2str(gg));
        legend('all gases','for this gas','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(2);
        plot(x.fc,diffBT0,'k',x.fc,rad2bt(x.fc,x.qc24)-rad2bt(x.fc,x.qc20),'r'); xlim(xlim2)
        title(num2str(gg));
        legend('all gases','for this gas','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(3);
        plot(x.fc,diffBT0,'k',x.fc,rad2bt(x.fc,x.qc24)-rad2bt(x.fc,x.qc20),'r'); xlim(xlim3)
        title(num2str(gg));
        legend('all gases','for this gas','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(4);
        plot(x.fc,diffBT0,'k',x.fc,rad2bt(x.fc,x.qc24)-rad2bt(x.fc,x.qc20),'r'); xlim(xlim4)
        title(num2str(gg));
        legend('all gases','for this gas','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
    else
      figure(1);
        plot(x.fc,allBT0,'k.-',x.fc,rad2bt(x.fc,x.qc20),'b',x.fc,rad2bt(x.fc,x.qc24),'r'); xlim(xlim1)
        title(num2str(gg));
        legend('all gases','H20','H24','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(2);
        plot(x.fc,allBT0,'k.-',x.fc,rad2bt(x.fc,x.qc20),'b',x.fc,rad2bt(x.fc,x.qc24),'r'); xlim(xlim2)
        title(num2str(gg));
        legend('all gases','H20','H24','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(3);
        plot(x.fc,allBT0,'k.-',x.fc,rad2bt(x.fc,x.qc20),'b',x.fc,rad2bt(x.fc,x.qc24),'r'); xlim(xlim3)
        title(num2str(gg));
        legend('all gases','H20','H24','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')    
      figure(4);
        plot(x.fc,allBT0,'k.-',x.fc,rad2bt(x.fc,x.qc20),'b',x.fc,rad2bt(x.fc,x.qc24),'r'); xlim(xlim4)
        title(num2str(gg));
        legend('all gases','H20','H24','location','best');
        xlabel('Wavenumber cm-1'); ylabel('BTD (K)')
    end
    pause;
  end
end
