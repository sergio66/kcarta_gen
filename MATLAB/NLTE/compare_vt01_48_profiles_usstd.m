figure(1); clf
sol = 0;
str = ['s' num2str(sol)];
for ii = 1 : 48
  fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
  fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt' num2str(ii) '_' str '.prf'];
  [p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);  
  p48(:,ii)       = p';
  T48(:,ii)       = T'; 
  QV48(:,:,ii)    = QV;
  Tnlte48(:,:,ii) = Tnlte;
  band48(:,ii)    = bandcenter';
  hold on
end

ii = 1;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_1_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 2;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 3;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

hold off; title('OLD'); grid

%{
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_1_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold off; title('NEW'); grid
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% now interp everything onto pStd
for ii = 1 : 48
  px = p48(:,ii);
  tx = T48(:,ii);
  ty = interp1(log10(px),tx,log10(pStd(:,2)),[],'extrap');
  Tnew48(:,ii) = ty;

  Q = squeeze(QV48(:,:,ii));
  Qy = interp1(log10(px),Q,log10(pStd(:,2)),[],'extrap');
  QVnew48(:,:,ii) = Qy;

  tx = squeeze(Tnlte48(:,:,ii));
  ty = interp1(log10(px),tx,log10(pStd(:,2)),[],'extrap');
  Tnltenew48(:,:,ii) = ty;
end

figure(2); 
  semilogy(Tnew48,p,'b'); set(gca,'ydir','reverse')
  semilogy(Tnew48 - TStd(:,2)*ones(1,48),p,'b'); 
    set(gca,'ydir','reverse');  grid
    line([-50 -50],[log(1000) log(1e-4)],'color','k')
    line([+50 +50],[log(1000) log(1e-4)],'color','k')

fprintf(1,'bandcenter(1,1) for USStd = %812.4f \n',bandStd(1,1));
ix48 = find(band48(:,1) == bandStd(1,1));
  semilogy(Tnew48 - TStd(:,2)*ones(1,48),p,'b',...
           squeeze(Tnltenew48(:,ix48,:)) - squeeze(TnlteStd(:,1,2))*ones(1,48),p,'r'); 
  set(gca,'ydir','reverse');     grid

pStdALL = pStd(:,2)*ones(1,48);
figure(3)
 aa = Tnew48 - TStd(:,2)*ones(1,48);
 bb = squeeze(Tnltenew48(:,ix48,:)) - squeeze(TnlteStd(:,1,2))*ones(1,48);
 cc = log10(pStd(:,2))*ones(1,48);
 scatter(aa(:),bb(:),20,cc(:),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');

%% this is how bb varies with aa in the 1100 - 0.005 LA usal kCARTA atm
iLA = find(pStd(:,2) >= 5e-3);
iProf=1; scatter(bb(iLA,iProf),aa(iLA,iProf),30,log10(pStd(iLA,2)),'filled'); grid
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press)'); colorbar
 lala = find(pStdALL > 0.005); lala = lala(:);
 scatter(aa(lala),bb(lala),20,cc(lala),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press) : ALL : 1100 - 0.005 mb'); colorbar; ret

%% this is how bb varies with aa in the 0.005 - 0.00005 mb UA usal kCARTA atm
iLA = find(pStd(:,2) <= 5e-3);
iProf=1; scatter(bb(iLA,iProf),aa(iLA,iProf),30,log10(pStd(iLA,2)),'filled'); grid
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press)'); colorbar
 lala = find(pStdALL < 0.005); lala = lala(:);
 scatter(aa(lala),bb(lala),20,cc(lala),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press) : UA : 0.005 - 0.0005mb'); colorbar; ret

%% this is how bb varies with aa in the transition region
lala = find(log10(pStd(:,2)) < 0.4 & log10(pStd(:,2)) >= -1.25);
iProf=1; scatter(bb(lala,iProf),aa(lala,iProf),30,log10(pStd(lala,2)),'filled'); grid
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press)'); colorbar
 lala = find(pStdALL < 10^0.4 & pStdALL > 10^-1.25); lala = lala(:);
 scatter(aa(lala),bb(lala),20,cc(lala),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press) : TRANSITION'); colorbar; ret

%% this is how bb varies with aa upto 10^0.4 mb
lala = find(log10(pStd(:,2)) >0.4);
iProf=1; scatter(bb(lala,iProf),aa(lala,iProf),30,log10(pStd(lala,2)),'filled'); grid
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press)'); colorbar
 lala = find(pStdALL > 10^0.4); lala = lala(:);
 scatter(aa(lala),bb(lala),20,cc(lala),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press) : 1100 - 5 mb'); colorbar; ret

%% this is how bb varies with aa at low pressures
lala = find(log10(pStd(:,2)) < -1.25 & pStd(:,2) > 0.005);
iProf=1; scatter(bb(lala,iProf),aa(lala,iProf),30,log10(pStd(lala,2)),'filled'); grid
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press)'); colorbar
 lala = find(log10(pStdALL) < -1.25 & pStdALL > 0.005); lala = lala(:);
 scatter(aa(lala),bb(lala),20,cc(lala),'filled')
grid; colorbar; xlabel('LTE - LTE(std)'); ylabel('NLTE - NLTE(std)');
xlabel('NLTE(profile)-NLTE(USStd)'); xlabel('LTE(profile)-LTE(USStd)');
title('colorbar = log10(press) : 5 mb - 0.005 mb'); colorbar; ret

figure(4)
  errorbar_x(nanmean(aa'),log10(pStd(:,2)),2*nanstd(aa'),'bo-'); hold on
  errorbar_x(nanmean(bb'),log10(pStd(:,2)),2*nanstd(bb'),'ro-'); hold off
  set(gca,'ydir','reverse'); grid
  line([-40 +40],[log10(5e-3) log10(5e-3)],'color','k')

%% this shows how std dev of (TLTE-TSTD) varies with p
%% this shows how std dev of (TNLTE-TNLTE-STD) varies with p
plot(log10(pStd(:,2)),2*nanstd(aa'),'bo-',log10(pStd(:,2)),2*nanstd(bb'),'ro-')
 line([log10(5e-3) log10(5e-3)],[0 40],'color','k'); grid

%% notice how NLTE becomes almost a constant value ie TNLTE - TNLTE-STD --> 0 at small p
boo = zeros(size(pStd(:,2)));
boo(log10(pStd(:,2)) >= 0.4) = 25;

boo(log10(pStd(:,2)) < 0.4 & log10(pStd(:,2)) >= -1.25) = 15;
x1=log10(0.4); y1=25; x2=log10(-1.25); y2=2.5; slope = (y1-y2)/(x1-x2);
x1=0.4; y1=25; x2=-1.25; y2=2.5; slope = (y1-y2)/(x1-x2);
lala = find(log10(pStd(:,2)) < 0.4 & log10(pStd(:,2)) >= -1.25);
xlala = log10(pStd(lala,2));
ylala = 25-slope*(x1-xlala);
boo(lala) = ylala;

boo(log10(pStd(:,2)) < -1.25) = 2.5;
plot(log10(pStd(:,2)),2*nanstd(aa'),'bo-',log10(pStd(:,2)),2*nanstd(bb'),'ro-',...
     log10(pStd(:,2)),boo,'ko-')
 line([log10(5e-3) log10(5e-3)],[0 40],'color','k'); grid

