%x = load('JUNK/xconvolved_kcarta_airs_put_together.mat');
umbc = load('JUNK/xconvolved_kcarta_airs_put_together_CKD25_umbcOD.mat');
lbl2 = load('JUNK/xconvolved_kcarta_airs_put_together_CKD25_lblrtmOD_varyT.mat');
lbl = load('JUNK/xconvolved_kcarta_airs_put_together_CKD25_lblrtmOD.mat');
umbc8 = load('JUNK/xconvolved_kcarta_airs_put_together_CKD6_umbcOD_H2008.mat');
[h,ha,p,pa] = rtpread('pert_tune_tunmlt_jan18_2016_6000tropicalprofiles.op.rtp');

g = dogoodchan;

tobs = rad2bt(umbc.fKc,p.robs1);
%tx = rad2bt(umbc.fKc,x.dKcall);
tl2 = rad2bt(umbc.fKc,lbl2.dKcall);
tl = rad2bt(umbc.fKc,lbl.dKcall);
tu = rad2bt(umbc.fKc,umbc.dKcall);
tu8 = rad2bt(umbc.fKc,umbc8.dKcall);

%% x and umbc should be the same
%[sum(sum(tx-tl)) sum(sum(tx-tl2)) sum(sum(tx-tu)) sum(sum(tx-tu8))]
tx
plot(umbc.fKc,nanmean(tu'-tl'),umbc.fKc,nanmean(tu'-tl2'),umbc.fKc,nanmean(tu'-tu8'),'linewidth',2)
axis([600 800 -0.5 +0.5]); grid
hl = legend('UMBC-LBL','UMBC-LBL,varyT','UMBC-UMBC H2008');
axis([640 800 -0.5 +0.5]); grid
axis([640 800 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([1200 1400 -0.5 +0.5]); grid
axis([1200 1400 -0.5 +0.5]); grid

plot(umbc.fKc,nanmean(tu8'-tl'),umbc.fKc,nanmean(tu8'-tl2'),umbc.fKc,nanmean(tu8'-tu'),'linewidth',2)
plot(umbc.fKc,nanmean(tu8'-tl'),umbc.fKc,nanmean(tu8'-tl2'),umbc.fKc,nanmean(tu8'-tu'),'linewidth',2)
hl = legend('UMBCH08-LBL','UMBCH08-LBLvaryT','UMBCH08-UMBCH12');
axis([1200 1400 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([640 800 -0.5 +0.5]); grid
axis([640 800 -0.5 +0.5]); grid
dbt1 = nanmean(tu8'-tl');
dbt2 = nanmean(tu8'-tl2');
dbt3 = nanmean(tu8'-tu');
comment = 'dbt1=tu8-tl, dbt2=tu8-tl2, dbt3=tu8-tu12 where tu8=umbcH08,tu12=umbH12,tl=lblrtm,tl2=lblrtm with linear-in-tau';
comment2 = 'see umbc_vs_lblrtm_6000tropical_jan18_2016.m';
save dbt_umbc_lblrtm.mat dbt1 dbt2 dbt3 comment comment2

plot(umbc.fKc,nanmean(tu8'-tl'),umbc.fKc,nanmean(tu8'-tl'),umbc.fKc,nanmean(tu8'-tu'),'linewidth',2)
plot(umbc.fKc(g),nanmean(tobs(g,:)'-tl(g,:)'),'b',umbc.fKc(g),nanmean(tobs(g,:)'-tl2(g,:)'),'g',umbc.fKc(g),nanmean(tobs(g,:)'-tu(g,:)'),'r+-',...
     umbc.fKc(g),nanmean(tobs(g,:)'-tu8(g,:)'),'k','linewidth',2)
hl = legend('tobs-LBL','tobs-LBLvaryT','tobs-UMBCH12','tobs-UMBCH08');
axis([1200 1400 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([2200 2420 -0.5 +0.5]); grid
axis([640 800 -0.5 +0.5]); grid
axis([640 800 -1.5 +1.5]); grid
