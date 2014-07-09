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


