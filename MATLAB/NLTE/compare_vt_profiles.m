figure(1)
fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
%fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt5_s0.prf'];
fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt1_s0.prf'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold on

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_3_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'nlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold off; title('OLD'); grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
fname = ['/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/'];
%fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt5_s0.prf'];
%fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt1_s0.prf'];
fname = [fname 'VT_48PROFILES_120_400ppmv/sergio_mergeMAIN/vt10_s0.prf'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold on

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_3_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_6_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

fname = ['/home/sergio/HITRAN2UMBCLBL/MAKEIR_4umNLTE/NLTEProfs/IPFILES/0/'];
fname = [fname 'xnlte_1_1_1_11_sol_0.genln2'];
[p,T,QV,Tnlte,bandcenter] = read_vt_profiles(fname);

hold off; title('NEW'); grid
