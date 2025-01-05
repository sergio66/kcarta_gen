addpath /asl/matlib/h4tools
addpath /asl/matlib/science
addpath /asl/matlib/aslutil

[h,~,pN,~] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/TestNalliEMiss/pnalli_aug30_2020_nalli30_65profs.op.rtp');
[h,~,pM,~] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/TestNalliEMiss/pnalli_aug30_2020_masuda19_65profs.op.rtp');

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat

for ii = 1 : 60
  fN = ['Nalli30/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  loader = ['x = load(''' fN ''');'];
  eval(loader)
  rcalc_nalli(:,ii) = x.rKc;

  fM = ['Masuda19/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  loader = ['x = load(''' fM ''');'];
  eval(loader)
  rcalc_masuda(:,ii) = x.rKc;
end

fKc     = x.fKc(ichan);
tobs    = rad2bt(h.vchan,pN.robs1(:,1:60));
tmasuda = rad2bt(x.fKc(ichan),rcalc_masuda(ichan,:));
tnalli = rad2bt(x.fKc(ichan),rcalc_nalli(ichan,:));
plot(fKc,nanmean(tobs'-tmasuda'),'b',fKc,nanmean(tobs'-tnalli'),'r',...
     fKc,nanstd(tobs'-tmasuda'),'c',fKc,nanstd(tobs'-tnalli'),'m')
grid
hl = legend('Masuda Bias','Nalli Bias','Masuda Std','Nalli Std','location','best');

%{
rcalc_masudax = rcalc_masuda(ichan,:);
rcalc_nallix  = rcalc_nalli(ichan,:);
save /asl/ftp/pub/sergio/NalliEmiss/nalli_sept18_2020.mat pN pM rcalc_nallix rcalc_masudax fKc
%}
