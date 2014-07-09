figure(1); clf
sol = 0;
str = ['s' num2str(sol)];
%% manuel
for ii = 1 : 48
  fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
  fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt' num2str(ii) '_' str '.prf'];
  [p(:,ii),T(:,ii),QV(:,:,ii),Tnlte(:,:,ii),bandcenter(:,ii),iso(:,ii)] = read_vt_profiles(fname);  
  hold on
end
hold off; title('MANUEL'); grid

for ii = 2 : 48
  junk = T(:,ii);
  newjunk = interp1(log10(p(:,ii)),junk,log10(p(:,1)),[],'extrap');
  T(:,ii) = newjunk;
end

[mm,nn,oo] = size(QV);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(QV(:,jj,ii));
    newjunk = interp1(log10(p(:,ii)),junk,log10(p(:,1)),[],'extrap');
    QV(:,jj,ii) = newjunk;
  end
end

[mm,nn,oo] = size(Tnlte);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(Tnlte(:,jj,ii));
    newjunk = interp1(log10(p(:,ii)),junk,log10(p(:,1)),[],'extrap');
    Tnlte(:,jj,ii) = newjunk;
  end
end

for ii = 2 : 48
  p(:,ii) = p(:,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
%% sergio
for ii = 1 : 48
  fname = ['/strowdata1/shared/sergio/kcartaV118/SRCv1.16/NONLTE/sergio/AIRSDATA_NLTE_VIBTEMPS_Puertas/'];
  fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/REGR49/0/'];
  fname = [fname 'regr49_nlte_1_1_1_' num2str(ii) '_sol_0.genln2'];
  [ps(:,ii),Ts(:,ii),QVs(:,:,ii),Tnltes(:,:,ii),bandcenters(:,ii),isos(:,ii)] = read_vt_profiles(fname);  
  hold on
end
hold off; title('SERGIO'); grid

for ii = 2 : 48
  junk = Ts(:,ii);
  newjunk = interp1(log10(ps(:,ii)),junk,log10(ps(:,1)),[],'extrap');
  Ts(:,ii) = newjunk;
end

[mm,nn,oo] = size(QVs);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(QVs(:,jj,ii));
    newjunk = interp1(log10(ps(:,ii)),junk,log10(ps(:,1)),[],'extrap');
    QVs(:,jj,ii) = newjunk;
  end
end

[mm,nn,oo] = size(Tnltes);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(Tnltes(:,jj,ii));
    newjunk = interp1(log10(ps(:,ii)),junk,log10(ps(:,1)),[],'extrap');
    Tnltes(:,jj,ii) = newjunk;
  end
end

for ii = 2 : 48
  ps(:,ii) = ps(:,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3); clf

thebandx = [2.3491    2.2835    2.3321    2.3400]*1000;
%
% see top of knonlte.f
%
% !   IL    MOL   IDMOL  ISO   IDISO  LEVEL   IDAFGL  ENERGY(cm-1)
% !    8 CO2       2    626      1       11      9    2349.14300
% !   10 CO2       2    626      1     1111     16    3004.01200
% !   11 CO2       2    626      1    10012     23    3612.84200
% !   12 CO2       2    626      1     2211     24    3659.27300
% !   13 CO2       2    626      1    10011     25    3714.78300
% !   18 CO2       2    636      2       11      9    2283.48800
thebandx = [2349.14300 2283.48800 3004.01200 3612.84200 3659.27300 3714.78300];
theisox  = [1          2          1          1          1          1         ];
for jj = 1 : length(thebandx)
  boo  = find((abs(bandcenter(:,1)-thebandx(jj)) < 0.1) & (iso(:,1) == theisox(jj)),1);
  boos = find((abs(bandcenters(:,1)-thebandx(jj)) < 0.1) & (isos(:,1) == theisox(jj)),1);
  cumok = [];

  for ii = 1 : 48
    Tss(:,ii) = interp1(log10(ps(:,ii)),Ts(:,ii),log10(p(:,ii)),[],'extrap');
    Tnltess(:,ii) = interp1(log10(ps(:,ii)),squeeze(Tnltes(:,boos,ii)),log10(p(:,ii)),[],'extrap');
  end

  poo = squeeze(Tnlte(:,boo,:));
  yuck = nanstd(Tnltess'-squeeze(Tnlte(:,boo,:))');
  muck = nanmean(Tnltess'-squeeze(Tnlte(:,boo,:))');

  for kk = 1:length(p(:,1))
    ok  = find( ((Tnltess(kk,:) - poo(kk,:)) >= (muck(kk) - 3*yuck(kk)) ) & ...
                ((Tnltess(kk,:) - poo(kk,:)) <= (muck(kk) + 3*yuck(kk)) ));

    if kk == 1
      cumok = ok;
    else
      cumok = intersect(cumok,ok);
    end
    plot(1:48,abs(Tnltess(kk,:)-poo(kk,:)),'ro',1:48,yuck(kk)*ones(1,48),'k',...
                                                1:48,yuck(kk)*ones(1,48)*3,'gx'); 
    grid; title(num2str(kk))
  end
    
  themeanDT(:,jj) = nanmean(Tnltess(:,cumok)'-squeeze(Tnlte(:,boo,cumok))');

  semilogy(Tss-T,p,'b',Tnltess-squeeze(Tnlte(:,boo,:)),p,'r',themeanDT(:,jj),p(:,1),'ko-'); 
    set(gca,'ydir','reverse'); grid  
  line([-70 +70],[5e-3 5e-3],'color','k','linewidth',2)
  title([num2str(thebandx(jj))  ' cm-1 : sergio - manuel (b) LTE (r) NLTE'])
  axis([-10 +10 1e-4 1000])
  ret
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4); clf
addpath  /home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs
%% now reverse engineer files
for ii = 1 : 48
  fname = ['/strowdata1/shared/sergio/kcartaV118/SRCv1.16/NONLTE/sergio/AIRSDATA_NLTE_VIBTEMPS_Puertas/'];
  fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/REGR49/0/'];
  fnameOUT = [fname 'NEWregr49_nlte_1_1_1_' num2str(ii) '_sol_0.genln2'];
  fname = [fname 'regr49_nlte_1_1_1_' num2str(ii) '_sol_0.genln2'];
  [ps2(:,ii),Ts2(:,ii),QVs2(:,:,ii),Tnltes2(:,:,ii),bandcenters2(:,ii),isos2(:,ii)] = ...
     adjust_vt_profiles(fname,fnameOUT,thebandx,theisox,themeanDT,mean(p'));    
  hold on
end

for ii = 2 : 48
  junk = Ts2(:,ii);
  newjunk = interp1(log10(ps2(:,ii)),junk,log10(ps2(:,1)),[],'extrap');
  Ts2(:,ii) = newjunk;
end

[mm,nn,oo] = size(QVs2);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(QVs2(:,jj,ii));
    newjunk = interp1(log10(ps2(:,ii)),junk,log10(ps2(:,1)),[],'extrap');
    QVs2(:,jj,ii) = newjunk;
  end
end

[mm,nn,oo] = size(Tnltes2);
for ii = 2 : oo
  for jj = 1 : nn
    junk = squeeze(Tnltes2(:,jj,ii));
    newjunk = interp1(log10(ps2(:,ii)),junk,log10(ps2(:,1)),[],'extrap');
    Tnltes2(:,jj,ii) = newjunk;
  end
end

for ii = 2 : 48
  ps2(:,ii) = ps2(:,1);
end

figure(5); clf
for jj = 1 : length(thebandx)
  boo  = find((abs(bandcenter(:,1)-thebandx(jj)) < 0.1) & (iso(:,1) == theisox(jj)),1);
  boos = find((abs(bandcenters(:,1)-thebandx(jj)) < 0.1) & (isos(:,1) == theisox(jj)),1);
  boos2 = find((abs(bandcenters2(:,1)-thebandx(jj)) < 0.1) & (isos2(:,1) == theisox(jj)),1);
  cumok = [1:48];
  cumok = [];

  for ii = 1 : 48
    Tss(:,ii) = interp1(log10(ps(:,ii)),Ts(:,ii),log10(p(:,ii)),[],'extrap');
    Tnltess(:,ii)  = interp1(log10(ps(:,ii)),squeeze(Tnltes(:,boos,ii)),log10(p(:,ii)),[],'extrap');
    Tnltess2(:,ii) = interp1(log10(ps2(:,ii)),squeeze(Tnltes2(:,boos2,ii)),log10(p(:,ii)),[],'extrap');
  end

  poo = squeeze(Tnlte(:,boo,:));
  yuck = nanstd(Tnltess'-squeeze(Tnlte(:,boo,:))');
  muck = nanmean(Tnltess'-squeeze(Tnlte(:,boo,:))');

  for kk = 1:length(p(:,1))
    ok  = find( ((Tnltess(kk,:) - poo(kk,:)) >= (muck(kk) - 3*yuck(kk)) ) & ...
                ((Tnltess(kk,:) - poo(kk,:)) <= (muck(kk) + 3*yuck(kk)) ));

    if kk == 1
      cumok = ok;
    else
      cumok = intersect(cumok,ok);
    end
    plot(1:48,abs(Tnltess(kk,:)-poo(kk,:)),'ro',1:48,yuck(kk)*ones(1,48),'k',...
                                                1:48,yuck(kk)*ones(1,48)*3,'gx'); 
    grid; title(num2str(kk))
  end
    
  themeanDT2(:,jj) = nanmean(Tnltess2(:,cumok)'-squeeze(Tnlte(:,boo,cumok))');

  semilogy(Tss-T,p,'b',Tnltess-squeeze(Tnlte(:,boo,:)),p,'r',themeanDT(:,jj),p(:,1),'ko-',...
                      Tnltess2-squeeze(Tnlte(:,boo,:)),p,'m',themeanDT2(:,jj),p(:,1),'gx-')
    set(gca,'ydir','reverse'); grid  
  line([-70 +70],[5e-3 5e-3],'color','k','linewidth',2)
  title([num2str(thebandx(jj))  ' cm-1 : sergio - manuel (b) LTE (r) NLTE'])
  axis([-10 +10 1e-4 1000])
  ret
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

comment = 'see /home/sergio/KCARTA/MATLAB/compare_vt01_48_profiles_sergioVSmanuel.m';
thep = mean(p');
save /home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/adjustNLTE_VT.mat comment thebandx  theisox themeanDT thep