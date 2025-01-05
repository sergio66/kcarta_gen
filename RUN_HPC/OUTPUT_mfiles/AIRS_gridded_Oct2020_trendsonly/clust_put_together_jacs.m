addpath /asl/matlib/h4tools

JOB = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%JOB = 32

miaow = load('sarta_chans_for_l1c.mat');
ind2834to2645 = miaow.ichan;

SARTAjac = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/subjacLatBin' num2str(JOB,'%02i') '.mat'];
sarta = load(SARTAjac);

foutsubjac  = ['kcarta_subjacLatBin' num2str(JOB,'%02i') '.mat'];
foutsubjac  = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/kcarta_subjacLatBin' num2str(JOB,'%02i') '.mat'];
foutsubjac2 = ['/asl/s1/sergio/rtp/MakeAvgProfs2002_2020/Retrieval/LatBin65/SubsetJacLatbin/kcarta_subjac_nostruct_LatBin' num2str(JOB,'%02i') '.mat'];

[h,ha,p,pa] = rtpread('/asl/s1/sergio/MakeAvgProfs2002_2020/summary_17years_all_lat_all_lon_2002_2019_palts.rtp');
plot(1:72,p.stemp(sarta.subjac.indices),'b-',1:72,sarta.subjac.stemp,'r.-')
pause(1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind_lat_junk = JOB;
ind_subset_junk = (1:72) + (ind_lat_junk-1)*72;

for lon = 1 : 72
  if mod(lon,10) == 0
    fprintf(1,'o');
  else
    fprintf(1,'.')
  end

  ind_subset_junk_ii = ind_subset_junk(lon);
  fz   = ['AllDemJacs/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_jac.mat'];
  fcol = ['AllDemJacs/individual_prof_convolved_kcarta_airs_' num2str(ind_subset_junk_ii) '_coljac.mat'];
  iaIndices(lon) = ind_subset_junk_ii;

  az = load(fz);
  acol = load(fcol); acol.rKc(:,1:4) = acol.rKc(:,1:4) * log(10);
  
  clear aout
  aout.fKc = az.fKc;
  aout.jac(:,1) = acol.rKc(:,1);     %%% CO2
  aout.jac(:,2) = acol.rKc(:,3);     %%% N2O
  aout.jac(:,3) = acol.rKc(:,4);     %%% CH4
  aout.jac(:,4) = acol.rKc(:,3) * 0; %%% cng1
  aout.jac(:,5) = acol.rKc(:,3) * 0; %%% cng2
  aout.jac(:,6) = acol.rKc(:,6);     %%% ST
  
  [~,numlays] = size(az.rKc);
  numlays = (numlays-4)/4;
  
  figure(1); plot(aout.fKc,aout.jac); 
  
  ind = (1:numlays);
  
  ix = ind + numlays*0; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*log(10);  %% WV
    [an(1) an(end) ix(1) ix(end)]; 
    wvind = fliplr(an);
    figure(2); plot(aout.fKc,sum(az.rKc(:,ix),2)*log(10),'r',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacWV(1:numlays,:,lon)),1),'g'); title('WV')

  ix = ind + numlays*2; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix);  %% T
    [an(1) an(end) ix(1) ix(end)];
    tzind = fliplr(an);
    figure(3); plot(aout.fKc,sum(az.rKc(:,ix),2),'r',aout.fKc,acol.rKc(:,5),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacT(1:numlays,:,lon)),1),'g'); title('T')

  ix = ind + numlays*1; [~,an] = size(aout.jac); an = an + (1:numlays); an = fliplr(an); aout.jac(:,an) = az.rKc(:,ix)*log(10);  %% O3
    [an(1) an(end) ix(1) ix(end)];
    o3ind = fliplr(an);
    figure(4); plot(aout.fKc,sum(az.rKc(:,ix),2)*log(10),'r',aout.fKc,acol.rKc(:,2),'b.-',...
                    aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacO3(1:numlays,:,lon)),1),'g'); title('O3')

  ixx = numlays*4 + 1;  %% surface temp
  figure(1); plot(aout.fKc,aout.jac,aout.fKc,az.rKc(:,ixx));   
  figure(5); plot(aout.fKc,az.rKc(:,ixx),'r',aout.fKc,aout.jac(:,6),'b',...
                  aout.fKc(ind2834to2645),squeeze(sarta.subjac.jacST(:,lon)),'g.');  title('SurfT')
  figure(6); plot(aout.fKc,aout.jac(:,1),'b',...
                  aout.fKc(ind2834to2645),sum(squeeze(sarta.subjac.jacCO2z(1:numlays,:,lon)),1),'g'); title('CO2')

  aout.fKc = aout.fKc(ind2834to2645);
  aout.jac = aout.jac(ind2834to2645,:);

  kcarta.subjac.coljacCH4(:,lon) = aout.jac(:,3);
  %kcarta.subjac.coljacCO(:,lon) = aout.jac(:,X);
  kcarta.subjac.coljacN2O(:,lon) = aout.jac(:,2);
  kcarta.subjac.jacCld1(:,lon)   = sarta.subjac.jacCld1(:,lon);
  kcarta.subjac.jacCld2(:,lon)   = sarta.subjac.jacCld2(:,lon);
  kcarta.subjac.jacSze1(:,lon)   = sarta.subjac.jacSze1(:,lon);
  kcarta.subjac.jacSze2(:,lon)   = sarta.subjac.jacSze2(:,lon);
  kcarta.subjac.jacTop1(:,lon)   = sarta.subjac.jacTop1(:,lon);
  kcarta.subjac.jacTop2(:,lon)   = sarta.subjac.jacTop2(:,lon);
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

  %disp('ret'); pause
  pause(0.1)
end
fprintf(1,'\n')

kcarta.subjac.rlon = sarta.subjac.rlon;
kcarta.subjac.rlat = sarta.subjac.rlat;
kcarta.subjac.stemp = sarta.subjac.stemp;
kcarta.subjac.nlevs = sarta.subjac.nlevs;
kcarta.subjac.ppmv2 = sarta.subjac.ppmv2;
kcarta.subjac.ppmv4 = sarta.subjac.ppmv4;
kcarta.subjac.ppmv6 = sarta.subjac.ppmv6;
kcarta.subjac.indices = iaIndices;

sum(sarta.subjac.indices - iaIndices)
  
subjac = kcarta.subjac;

kcarta.subjac.comment = 'see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/JUNK/AIRS_gridded_Oct2020_trendsonly/clust_put_together_jacs.m';

subjac = kcarta.subjac;
clear sarta kcarta

saver = ['save ' foutsubjac ' subjac'];
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
