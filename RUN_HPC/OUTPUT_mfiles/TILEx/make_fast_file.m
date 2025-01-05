addpath /home/sergio/KCARTA/MATLAB

l1c = load('/home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat');
ii = 1;
  f0   = ['individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  f0 = load(f0);
  l1c.ichan([1520 2600])            %% should be 1291 and 2333
  f0.fKc(l1c.ichan([1520 2600]))    %% should be 1231 and 2616 cm-1
  plot(f0.fKc(l1c.ichan))
  plot(f0.fKc(l1c.ichan),rad2bt(f0.fKc(l1c.ichan),f0.rKc(l1c.ichan)))

for ii = 1 : 400
  if mod(ii,100) == 0
    fprintf(1,'+')
  elseif mod(ii,10) == 0
    fprintf(1,'.')
  end

  f0   = ['individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  fcol = ['individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat'];
  fw   = ['individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];   %% has comnibes 1,101,102,103 to be g1, so the jacs are WV,O3,T,WGT and 4 surf

  f0 = load(f0);
  fcol = load(fcol);
  fw = load(fw);

  t0   = rad2bt(f0.fKc,f0.rKc);
  tcol = rad2bt(f0.fKc,fcol.rKc);   %% 2 4 5 6 51 52

  logx = log(1+0.001);; %% yeah that is what kCARTA/INCLUDE/pre_defined.params says is the size of the col mult
  jacCO2(ii,:) = (tcol(l1c.ichan,1)-t0(l1c.ichan))/logx;
  jacN2O(ii,:) = (tcol(l1c.ichan,2)-t0(l1c.ichan))/logx;
  jacCO(ii,:)  = (tcol(l1c.ichan,3)-t0(l1c.ichan))/logx;
  jacCH4(ii,:) = (tcol(l1c.ichan,4)-t0(l1c.ichan))/logx;
  jacCFC11(ii,:) = (tcol(l1c.ichan,5)-t0(l1c.ichan))/logx;
  jacCFC12(ii,:) = (tcol(l1c.ichan,6)-t0(l1c.ichan))/logx;

  ind = (1:97)+0*97; jacWV(ii,:,:) = fw.rKc(l1c.ichan,ind);
  ind = (1:97)+1*97; jacO3(ii,:,:) = fw.rKc(l1c.ichan,ind);
  ind = (1:97)+2*97; jacT(ii,:,:)  = fw.rKc(l1c.ichan,ind);
  ind = (1:97)+3*97; 
    ind = ind(end); ind = ind+1; jacST(ii,:) = fw.rKc(l1c.ichan,ind);
end

fprintf(1,'\n');

%{
fairs = f0.fKc(l1c.ichan);
iairs = l1c.ichan;
comment = 'see eg /home/sergio/MATLABCODE/oem_pkg_run_sergio_AuxJacs/MakeAvgCldProfs2002_2020/LookAtTimeSeries/RTP_PROFSV2/Cld/CRODGERS_FAST_CLOUD_Retrievals/V4Jacs//jacsV4_retr_gran_lonbin_67_latbin_35_JOB_2515.mat';
save kcartajacs_clearsky_lonbin_67_latbin_35_JOB_2515.mat comment fairs iairs jacWV jacO3 jacT jacST jacCO2 jacN2O jacCH4
%}

