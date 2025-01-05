%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrERA5.m
%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Aug2022_startSept2002_endAug2014_trendsonly/clust_put_together_jacs_clrERA5.m

addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

%%%%%%%% ORIG CODE %%%%%%%
%% recall log10(X) = log(X)/log(10)
%% kcarta does Q d(BT)/dQ == dBT/d(logQ) = dBT/dQ/Q = Q d(BT)/dQ
%% but if we change to log10 then dBT/dlog10(Q) = dBT/d(log(Q)/log(10)) = log10 dBT/dlog(Q) = log10 Q dBT/dQ = log10 dBT/dlog(Q)
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%% disp('WARNING here we use log10 jacs ie multiply gas jacs by log_e_(10) ie jac --> jac * log(10) = 2.3026')
%%%%%%%% ORIG CODE %%%%%%%

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 32

miaow = load('sarta_chans_for_l1c.mat');
ind2834to2645 = miaow.ichan;

iOldORNew = +1;   %% the 17 year ERA-I
iOldORNew = +2;   %% the 19 year ERA5 2002/09-2021/08
iOldORNew = +12;  %% the 12 year ERA5 2002/09-2014/08
iOldORNew = +07;  %% the 07 year ERA5 2012/05-2019/04
iOldORNew = +20;  %% the 20 year ERA5 2002/09-2022/08

if iOldORNew < 0
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  %%% NEED TO REDO
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin' num2str(JOB,'%02i') '.mat'];
elseif iOldORNew == 1
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_newSARTA_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_newSARTA_' num2str(JOB,'%02i') '.mat'];
elseif iOldORNew == 2
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_Dec2021_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_Dec2021_' num2str(JOB,'%02i') '.mat'];
elseif iOldORNew == 12
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjac12yearLatBin' num2str(JOB,'%02i') '.mat'];  
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_12yr_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_12yr_' num2str(JOB,'%02i') '.mat'];
elseif iOldORNew == 07
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjac12yearLatBin' num2str(JOB,'%02i') '.mat'];  
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_07yr_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_07yr_' num2str(JOB,'%02i') '.mat'];
elseif iOldORNew == 20
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjac20yearLatBin' num2str(JOB,'%02i') '.mat'];  
  SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];  
  foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjacLatBin_kCARTA_ERA5_20yr_' num2str(JOB,'%02i') '.mat'];
  foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020_startSept2002/Retrieval/LatBin65/SubsetJacLatbin/kcarta_clr_subjac_nostruct_LatBin_kCARTA_ERA5_20yr_' num2str(JOB,'%02i') '.mat'];
else
  error('unknown iOldORNew')
end

if exist(foutsubjac)
%  foutsubjac
%  error('foutsubjac exists')
end
if exist(foutsubjac2)
%  foutsubjac2
%  error('foutsubjac2 exists')
end

sarta = load(SARTAjac);

if iOldORNew == 0
  [h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020_startSept2002/summary_17years_all_lat_all_lon_2002_2019.rtp');  %% already has palts
elseif iOldORNew == 1
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_17years_all_lat_all_lon_2002_2019_palts_startSept2002_CLEAR.rtp');
elseif iOldORNew == 2
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_19years_all_lat_all_lon_2002_2021_monthlyERA5.rp.rtp');
elseif iOldORNew == 12
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_12years_all_lat_all_lon_2002_2014_monthlyERA5.rp.rtp');
elseif iOldORNew == 07
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_07years_all_lat_all_lon_2012_2019_monthlyERA5.rp.rtp');
elseif iOldORNew == 20
  [h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/summary_20years_all_lat_all_lon_2002_2022_monthlyERA5.rp.rtp');
else
  error('unknown iOldORNew')
end

plot(1:72,p.stemp(sarta.subjac.indices),'b-',1:72,sarta.subjac.stemp,'r.-')
pause(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

factor_log10 = log(10); %% this changes gas jac scaling to log10
factor_log10 = 1.0;     %% this keeps   gas jac scaling to loge

ind_lat_junk = JOB;
ind_subset_junk = (1:72) + (ind_lat_junk-1)*72;

for lon = 1 : 72
  if mod(lon,10) == 0
    fprintf(1,'o');
  else
    fprintf(1,'.')
  end

  ind_subset_junk_ii = ind_subset_junk(lon);
  frad0 = ['AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '.mat'];
  fz    = ['AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_jac.mat'];
  fcol  = ['AllDemJacsClr/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_coljac.mat'];
  iaIndices(lon) = ind_subset_junk_ii;

  arad0 = load(frad0);
  az    = load(fz);
  acol  = load(fcol); acol.rKc(:,1:6) = acol.rKc(:,1:6) * factor_log10;
  
  trad0 = rad2bt(az.fKc,arad0.rKc);
  %tcol  = rad2bt(az.fKc,acol.rKc);

  clear aout
  aout.fKc = az.fKc;

  %% kcarta nml = 2 4 5 6 51 52 T ST
  aout.jac(:,1) = acol.rKc(:,1);     %%% CO2
  aout.jac(:,2) = acol.rKc(:,2);     %%% N2O
  aout.jac(:,3) = acol.rKc(:,4);     %%% CH4
  aout.jac(:,4) = acol.rKc(:,5);     %%% CFC11
  aout.jac(:,5) = acol.rKc(:,6);     %%% CFC12
  aout.jac(:,6) = acol.rKc(:,8);     %%% ST

%  aout.jac(:,1) = (tcol(:,1)-trad0)*factor_log10;     %%% CO2
%  aout.jac(:,2) = (tcol(:,3)-trad0)*factor_log10;     %%% N2O
%  aout.jac(:,3) = (tcol(:,4)-trad0)*factor_log10;     %%% CH4
%  aout.jac(:,4) = tcol(:,3) * 0; %%% cng1
%  aout.jac(:,5) = tcol(:,3) * 0; %%% cng2
%  aout.jac(:,6) = tcol(:,6)-trad0;     %%% ST
  
  [~,numlays] = size(az.rKc);
  %% numlays = (numlays-6)/4;  %% see ~/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/template_Q2346_51_52jacobian.nml
  numlays = (numlays-4)/4;  %% see ~/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/template_Q2346_51_52jacobian.nml
  
  figure(1); plot(aout.fKc,aout.jac); 
  
  ind = (1:numlays);
  
  ix = ind + numlays*0; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*factor_log10;  %% WV
    [an(1) an(end) ix(1) ix(end)]; 
    wvind = fliplr(an);
    figure(2); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacWV(1:numlays,:,lon)),1),'g'); title('WV')
     hl = legend('sum(Kc(1:97))','sarta','location','best','fontsize',10);

  ix = ind + numlays*2; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix);  %% T
    [an(1) an(end) ix(1) ix(end)];
    tzind = fliplr(an);
    figure(3); plot(aout.fKc,sum(az.rKc(:,ix),2),'r',aout.fKc,acol.rKc(:,7),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacT(1:numlays,:,lon)),1),'g'); title('T')
%    figure(3); plot(aout.fKc,sum(az.rKc(:,ix),2),'r',aout.fKc,(tcol(:,7)-trad0),'b.-',...
%                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacT(1:numlays,:,lon)),1),'g'); title('T')
     hl = legend('sum(Kc(1:97))','col Kc','sarta','location','best','fontsize',10);

  figure(4); clf
  %% except we do not do col O3 in COL CLR jacs
  ix = ind + numlays*1; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*factor_log10;  %% O3
    [an(1) an(end) ix(1) ix(end)];
    o3ind = fliplr(an);
%{
    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',aout.fKc,acol.rKc(:,2),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')
%    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*factor_log10,'r',aout.fKc,(tcol(:,2)-trad0)*factor_log10,'b.-',...
%                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')
     hl = legend('sum(Kc(1:97))','col Kc','sarta','location','best','fontsize',10);
%}

  ixx = numlays*4 + 1;  %% surface temp
  figure(1); plot(aout.fKc,aout.jac,aout.fKc,az.rKc(:,ixx),'b');   
  figure(5); plot(aout.fKc,az.rKc(:,ixx),'ro-',aout.fKc,aout.jac(:,6),'bx-',...
                  aout.fKc(ind2834to2645),squeeze(sarta.subjac.jacST(:,lon)),'g.');  title('SurfT')
     hl = legend('Kc','for retr','sarta','location','best','fontsize',10);
  figure(6); plot(aout.fKc,aout.jac(:,1),'b',...
                  aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacCO2z(1:numlays,:,lon)),1),'g'); title('CO2')
     hl = legend('for retr','sarta','location','best','fontsize',10);

  aout.fKc = aout.fKc(ind2834to2645);
  aout.jac = aout.jac(ind2834to2645,:);

  kcarta.subjac.coljacCO2(:,lon)   = aout.jac(:,1);
  kcarta.subjac.coljacN2O(:,lon)   = aout.jac(:,2);
  %kcarta.subjac.coljacCO(:,lon)   = aout.jac(:,X);
  kcarta.subjac.coljacCH4(:,lon)   = aout.jac(:,3);
  kcarta.subjac.coljacCFC11(:,lon) = aout.jac(:,4);
  kcarta.subjac.coljacCFC12(:,lon) = aout.jac(:,5);
%  kcarta.subjac.jacCld1(:,lon)   = sarta.subjac.jacCld1(:,lon);
%  kcarta.subjac.jacCld2(:,lon)   = sarta.subjac.jacCld2(:,lon);
%  kcarta.subjac.jacSze1(:,lon)   = sarta.subjac.jacSze1(:,lon);
%  kcarta.subjac.jacSze2(:,lon)   = sarta.subjac.jacSze2(:,lon);
%  kcarta.subjac.jacTop1(:,lon)   = sarta.subjac.jacTop1(:,lon);
%  kcarta.subjac.jacTop2(:,lon)   = sarta.subjac.jacTop2(:,lon);
  kcarta.subjac.jacST(:,lon)     = aout.jac(:,6);
  if lon == 1
    kcarta.subjac.jacT = zeros(100,2645,72);
    kcarta.subjac.jacWV = zeros(100,2645,72);
    kcarta.subjac.jacO3 = zeros(100,2645,72);
    kcarta.subjac.jacCO2z = zeros(100,2645,72);

    kcarta.subjac.jacT(1:numlays,:,lon)  = aout.jac(:,tzind)';
    kcarta.subjac.jacWV(1:numlays,:,lon) = aout.jac(:,wvind)';
    kcarta.subjac.jacO3(1:numlays,:,lon) = aout.jac(:,o3ind)';
    kcarta.subjac.jacCO2z(1,:,lon)       = aout.jac(:,1);
  else
    kcarta.subjac.jacT(1:numlays,:,lon)  = aout.jac(:,tzind)';
    kcarta.subjac.jacWV(1:numlays,:,lon) = aout.jac(:,wvind)';
    kcarta.subjac.jacO3(1:numlays,:,lon) = aout.jac(:,o3ind)';
    kcarta.subjac.jacCO2z(1,:,lon)       = aout.jac(:,1);
  end

  plot(mean(kcarta.subjac.coljacN2O')); pause(0.1);
  %disp('ret'); pause
  pause(0.1)
end
fprintf(1,'\n')

addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS
i500mb = find(p.plevs(:,3000) >= 500,1);
ppmv2 = layers2ppmv(h,p,1:length(p.stemp),2);
ppmv4 = layers2ppmv(h,p,1:length(p.stemp),4);
ppmv6 = layers2ppmv(h,p,1:length(p.stemp),6);

% kcarta.subjac.rlon  = sarta.subjac.rlon;
% kcarta.subjac.rlat  = sarta.subjac.rlat;
% kcarta.subjac.stemp = sarta.subjac.stemp;
% kcarta.subjac.nlevs = sarta.subjac.nlevs;
% kcarta.subjac.ppmv2 = sarta.subjac.ppmv2;
% kcarta.subjac.ppmv4 = sarta.subjac.ppmv4;
% kcarta.subjac.ppmv6 = sarta.subjac.ppmv6;

kcarta.subjac.rlon  = p.rlon(ind_subset_junk);
kcarta.subjac.rlat  = p.rlat(ind_subset_junk);
kcarta.subjac.stemp = p.stemp(ind_subset_junk);
kcarta.subjac.nlevs = p.nlevs(ind_subset_junk);
kcarta.subjac.ppmv2 = ppmv2(i500mb,ind_subset_junk);
kcarta.subjac.ppmv4 = ppmv4(i500mb,ind_subset_junk);
kcarta.subjac.ppmv6 = ppmv6(i500mb,ind_subset_junk);
kcarta.subjac.indices = iaIndices;

sum(sarta.subjac.indices - iaIndices)
  
subjac = kcarta.subjac;

kcarta.subjac.comment{1} = 'see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Oct2020_startSept2002_trendsonly/clust_put_together_jacs_clr.m     for ERAI';
kcarta.subjac.comment{2} = 'see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Dec2021_startSept2002_trendsonly/clust_put_together_jacs_clrERA5.m for ERA5';
kcarta.subjac.comment{3} = 'see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Mar2023_startSept2002_endAug2022_trendsonlyclust_put_together_jacs_clrERA5.m for ERA5';

subjac = kcarta.subjac;
clear sarta kcarta

if ~exist(foutsubjac)
  saver = ['save ' foutsubjac ' subjac'];
  eval(saver)

  save(foutsubjac2,'-struct', 'subjac');
else
  fprintf(1,'%s already exists \n',foutsubjac)
end

disp('finished')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
