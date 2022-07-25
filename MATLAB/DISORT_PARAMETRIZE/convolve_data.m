addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB

iC = 1;  %% AIRS
iC = 0;  %% quickconvolve(w,rad,1,1)

iP = 1;  %% proff1
iPstr = ['Prof' num2str(iP)];

if iC == 0
  if ~exist([iPstr '/DISORT/conv4440.mat'])
    disp('loading and convolving disort radiances to QuickConvolve')
    clear w rad
    load([iPstr '/DISORT/disortrad.mat']);
    [fc,qc] = quickconvolve(w,rad,0.5,0.5);
    saver = ['save ' iPstr '/DISORT/conv4440.mat fc qc'];
    eval(saver)
  else
    disp('already have gaussian disort')
  end

  ix = [[0:5:95] [100 : 5 : 135]];
  for ii = 1 : length(ix)
    switch ii
      case 1
       str = '0.00/';
      case 2
       str = '0.05/';
      case 3
       str = '0.10/';
      case 4
       str = '0.15/';
      case 5
       str = '0.20/';
      case 6
       str = '0.25/';
      case 7
       str = '0.30/';
      case 8
       str = '0.35/';
      case 9
       str = '0.40/';
      case 10
       str = '0.45/';
      case 11
       str = '0.50/';
      case 12
       str = '0.55/';
      case 13
       str = '0.60/';
      case 14
       str = '0.65/';
      case 15
       str = '0.70/';
      case 16
       str = '0.75/';
      case 17
       str = '0.80/';
      case 18
       str = '0.85/';
      case 19
       str = '0.90/';
      case 20
       str = '0.95/';
      %%% these are testing the parametrization
      case 21
       str = 'TestLookupTable/';  %% tropical only
      case 22,23,24,25,26,27,28
       str = 'TestLookupTableArbProfiles/'; %% all 49 regression, DISORT, PCLSAM TANG, PCLSAM ONLY
    end
    %str = [str '/'];
    fprintf(1,'%2i str = %s \n',ii,str);
    if ii <= 21
      fileIN  = 'kcartarad.mat';
      fileOUT = 'conv4440.mat'; 

    elseif ii == 22
      fileIN  = 'kcartaradDISORT.mat';
      fileOUT = 'conv4440DISORT.mat'; 
    elseif ii == 23
      fileIN  = 'kcartaradPCLSAM_BEST.mat';
      fileOUT = 'conv4440PCLSAM_BEST.mat'; 
    elseif ii == 24
      fileIN  = 'kcartaradPCLSAM_NOADJ.mat';
      fileOUT = 'conv4440PCLSAM_NOADJ.mat'; 

    elseif ii == 25
      fileIN  = 'kcartarad_random_pert_DISORT.mat';
      fileOUT = 'conv4440_random_pert_DISORT.mat'; 
    elseif ii == 26
      fileIN  = 'kcartarad_random_pert_PCLSAM_BEST.mat';
      fileOUT = 'conv4440_random_pert_PCLSAM_BEST.mat'; 
    elseif ii == 27
      fileIN  = 'kcartarad_random_pert_PCLSAM_NOADJ.mat';
      fileOUT = 'conv4440_random_pert_PCLSAM_NOADJ.mat'; 
    elseif ii == 28
      fileIN  = 'kcartarad_random_pert_PCLSAM_BEST_VaryMatr.mat';
      fileOUT = 'conv4440_random_pert_PCLSAM_BEST_VaryMatr.mat'; 
    end

    fxin  = [iPstr '/Chou' str fileIN];
    fxout = [iPstr '/Chou' str fileOUT];
    if ~exist(fxout) & exist(fxin)
      fprintf(1,'loading and convolving PCLSAM radiances in %s to QuickConvolve %s \n',fxin,fxout);
      clear w rad
      load(fxin);
      [fc,qc] = quickconvolve(w,rad,0.5,0.5);
      saver = ['save ' fxout ' fc qc'];
      eval(saver)
    else
      disp('already have gaussian conv')
    end
  end

elseif iC == 1
  if ~exist([iPstr '/DISORT/airs2834.mat'])
    disp('loading and convolving disort radiances to AIRS')
    clear w rad
    load([iPstr '/DISORT/disortrad.mat']);
    [fc,qc] = convolve_airs(w,rad,1:2834);
    saver = ['save ' iPstr '/DISORT/airs2834.mat fc qc'];
    eval(saver)
  else
    disp('already have AIRS disort')
  end

  ix = [[0:5:95] [100 : 5 : 135]];
  for ii = 1 : length(ix)
    switch ii
      case 1
       str = '0.00/';
      case 2
       str = '0.05/';
      case 3
       str = '0.10/';
      case 4
       str = '0.15/';
      case 5
       str = '0.20/';
      case 6
       str = '0.25/';
      case 7
       str = '0.30/';
      case 8
       str = '0.35/';
      case 9
       str = '0.40/';
      case 10
       str = '0.45/';
      case 11
       str = '0.50/';
      case 12
       str = '0.55/';
      case 13
       str = '0.60/';
      case 14
       str = '0.65/';
      case 15
       str = '0.70/';
      case 16
       str = '0.75/';
      case 17
       str = '0.80/';
      case 18
       str = '0.85/';
      case 19
       str = '0.90/';
      case 20
       str = '0.95/';
      %%% these are testing the parametrization
      case 21
       str = 'TestLookupTable/';  %% tropical only
      case 22,23,24,25,26,27
       str = 'TestLookupTableArbProfiles/'; %% all 49 regression, DISORT, PCLSAM TANG, PCLSAM ONLY
    end
    %str = [str '/'];
    fprintf(1,'%2i str = %s \n',ii,str);
    if ii <= 21
      fileIN  = 'kcartarad.mat';
      fileOUT = 'airs2834.mat'; 

    elseif ii == 22
      fileIN  = 'kcartaradDISORT.mat';
      fileOUT = 'airs2834DISORT.mat'; 
    elseif ii == 23
      fileIN  = 'kcartaradPCLSAM_BEST.mat';
      fileOUT = 'airs2834PCLSAM_BEST.mat'; 
    elseif ii == 24
      fileIN  = 'kcartaradPCLSAM_NOADJ.mat';
      fileOUT = 'airs2834PCLSAM_NOADJ.mat'; 

    elseif ii == 25
      fileIN  = 'kcartarad_random_pert_DISORT.mat';
      fileOUT = 'airs2834_random_pert_DISORT.mat'; 
    elseif ii == 26
      fileIN  = 'kcartarad_random_pert_PCLSAM_BEST.mat';
      fileOUT = 'airs2834_random_pert_PCLSAM_BEST.mat'; 
    elseif ii == 27
      fileIN  = 'kcartarad_random_pert_PCLSAM_NOADJ.mat';
      fileOUT = 'airs2834_random_pert_PCLSAM_NOADJ.mat'; 
    elseif ii == 28
      fileIN  = 'kcartarad_random_pert_PCLSAM_BEST_VaryMatr.mat';
      fileOUT = 'airs2834_random_pert_PCLSAM_BEST_VaryMatr.mat'; 

    end

    fxin  = [iPstr '/Chou' str fileIN];
    fxout = [iPstr '/Chou' str fileOUT];
    if ~exist(fxout) & exist(fxin)
      fprintf(1,'loading and convolving PCLSAM radiances in %s to AIRSConvolve %s \n',fxin,fxout);
      clear w rad
      load(fxin);
      [fc,qc] = convolve_airs(w,rad,1:2834);
      saver = ['save ' fxout ' fc qc'];
      eval(saver)
    else
      disp('already have AIRS conv')
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
fdis    = [iPstr '/DISORT/conv4440.mat'];
ftest   = [iPstr '/ChouTestLookupTable/conv4440.mat'];
ftest0  = [iPstr '/Chou0.00/conv4440.mat'];
ftest30 = [iPstr '/Chou0.30/conv4440.mat'];
dis   = load(fdis);
pcl   = load(ftest);
pcl0  = load(ftest0);
pcl30 = load(ftest30);

disBT   = rad2bt(dis.fc,dis.qc);
pclBT   = rad2bt(pcl.fc,pcl.qc);
pcl0BT  = rad2bt(pcl.fc,pcl0.qc);
pcl30BT = rad2bt(pcl.fc,pcl30.qc);
plot(dis.fc,nanmean(disBT'-pcl0BT'),'b',dis.fc,nanstd(disBT'-pcl0BT'),'c--',dis.fc,nanmean(disBT'-pcl30BT'),'g',dis.fc,nanstd(disBT'-pcl30BT'),'y--',dis.fc,nanmean(disBT'-pclBT'),'r',dis.fc,nanstd(disBT'-pclBT'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','Tang bias','Tang std','New bias','New std','fontsize',8,'location','best');
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
check_arb_DIS_vs_PCL
%}
