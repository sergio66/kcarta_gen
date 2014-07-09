figure(1); clf
sol = 0;
str = ['s' num2str(sol)];
for ii = 1 : 48
  fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
  fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt' num2str(ii) '_' str '.prf'];
  [p(:,ii),T(:,ii),QV(:,ii),Tnlte(:,ii),bandcenter(:,ii)] = read_vt_profiles(fname);  
  hold on
end

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_1_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold off; title('OLD'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf
sol = 0;
str = ['s' num2str(sol)];
for ii = 1 : 48
  fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
  fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt' num2str(ii) '_' str '.prf'];
  [p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);  
  hold on
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
