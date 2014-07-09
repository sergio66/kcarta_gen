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

hold off; title('ORIG'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2); clf
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
fname = [fname 'xnlte_1_1_1_1_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 2;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 3;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

hold off; title('X NLTE'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3); clf
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
fname = [fname 'ynlte_1_1_1_1_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 2;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'ynlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

ii = 3;
fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'ynlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);
  pStd(:,ii)       = p';
  TStd(:,ii)       = T'; 
  QVStd(:,:,ii)    = QV;
  TnlteStd(:,:,ii) = Tnlte;
  bandStd(:,ii)    = bandcenter';

hold off; title('Y NLTE'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
