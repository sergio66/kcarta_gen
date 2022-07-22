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
    disp('already have AIRS disort')
  end

  ix = [0:5:100];
  for ii = 1 : length(ix)
    switch ii
      case 1
       str = '0.00';
      case 2
       str = '0.05';
      case 3
       str = '0.10';
      case 4
       str = '0.15';
      case 5
       str = '0.20';
      case 6
       str = '0.25';
      case 7
       str = '0.30';
      case 8
       str = '0.35';
      case 9
       str = '0.40';
      case 10
       str = '0.45';
      case 11
       str = '0.50';
      case 12
       str = '0.55';
      case 13
       str = '0.60';
      case 14
       str = '0.65';
      case 15
       str = '0.70';
      case 16
       str = '0.75';
      case 17
       str = '0.80';
      case 18
       str = '0.85';
      case 19
       str = '0.90';
      case 20
       str = '0.95';
      case 21
       str = 'TestLookupTable';
    end
    fprintf(1,'str = %s \n',str);
    if ~exist([iPstr '/Chou' str '/conv4440.mat']) & exist([iPstr '/Chou' str '/kcartarad.mat'])
      disp('loading and convolving PCLSAM radiances to QuickConvolve')
      clear w rad
      load([iPstr '/Chou' str '/kcartarad.mat']);
      [fc,qc] = quickconvolve(w,rad,0.5,0.5);
      saver = ['save ' iPstr '/Chou' str '/conv4440.mat fc qc'];
      eval(saver)
    else
      disp('already have AIRS PCLSAM')
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

  ix = [0:5:100];
  for ii = 1 : length(ix)
    switch ii
      case 1
       str = '0.00';
      case 2
       str = '0.05';
      case 3
       str = '0.10';
      case 4
       str = '0.15';
      case 5
       str = '0.20';
      case 6
       str = '0.25';
      case 7
       str = '0.30';
      case 8
       str = '0.35';
      case 9
       str = '0.40';
      case 10
       str = '0.45';
      case 11
       str = '0.50';
      case 12
       str = '0.55';
      case 13
       str = '0.60';
      case 14
       str = '0.65';
      case 15
       str = '0.70';
      case 16
       str = '0.75';
      case 17
       str = '0.80';
      case 18
       str = '0.85';
      case 19
       str = '0.90';
      case 20
       str = '0.95';
      case 21
       str = 'TestLookupTable';
    end
    fprintf(1,'str = %s \n',str);
    if ~exist([iPstr '/Chou' str '/airs2834.mat']) & exist([iPstr '/Chou' str '/kcartarad.mat'])
      disp('loading and convolving PCLSAM radiances to AIRS')
      clear w rad
      load([iPstr '/Chou' str '/kcartarad.mat']);
      [fc,qc] = convolve_airs(w,rad,1:2834);
      saver = ['save ' iPstr '/Chou' str '/airs2834.mat fc qc'];
      eval(saver)
    else
      disp('already have AIRS PCLSAM')
    end
  end

end

%{
fdis  = [iPstr '/DISORT/conv4440.mat'];
ftest = [iPstr '/ChouTestLookupTable/conv4440.mat'];
ftest0 = [iPstr '/Chou0.00/conv4440.mat'];
ftest30 = [iPstr '/Chou0.30/conv4440.mat'];
dis = load(fdis);
pcl = load(ftest);
pcl0 = load(ftest0);
pcl30 = load(ftest30);

disBT = rad2bt(dis.fc,dis.qc);
pclBT = rad2bt(pcl.fc,pcl.qc);
pcl0BT = rad2bt(pcl.fc,pcl0.qc);
pcl30BT = rad2bt(pcl.fc,pcl30.qc);
plot(dis.fc,nanmean(disBT'-pcl0BT'),'b',dis.fc,nanstd(disBT'-pcl0BT'),'c--',dis.fc,nanmean(disBT'-pcl30BT'),'g',dis.fc,nanstd(disBT'-pcl30BT'),'y--',dis.fc,nanmean(disBT'-pclBT'),'r',dis.fc,nanstd(disBT'-pclBT'),'m--')
axis([605 1650 -1 +1])
plotaxis2;
ylabel('DISORT - PCLSAM (K)'); xlabel('Wavenumber cm-1')
hl = legend('Orig bias','Orig std','Tang bias','Tang std','New bias','New std','fontsize',8,'location','best');
%}
