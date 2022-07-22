addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

fairs = instr_chans;
nedt  = instr_chans('airs',2);
g = dogoodchan;
g2 = find(nedt < 1);
g2 = intersect(g2,g);
plot(fairs(g2),nedt(g2))

iP = 1;  %% proff1
iPstr = ['Prof' num2str(iP)];

[h,ha,p,pa] = rtpread(['parametrize_pclsam_disort_AFGL_' num2str(iP) '.rp.rtp']);

iC = 1;  %% AIRS
iC = 0;  %% quickconvolve(w,rad,1,1)

if iC == 1
  disp('loading AIRS convolutions')
  x = load([iPstr '/DISORT/airs2834.mat']);   disort = x.qc; f = x.fc;
  x = load([iPstr '/Chou0.00/airs2834.mat']); a0p00 = x.qc;
  x = load([iPstr '/Chou0.05/airs2834.mat']); a0p05 = x.qc;
  x = load([iPstr '/Chou0.10/airs2834.mat']); a0p10 = x.qc;
  x = load([iPstr '/Chou0.15/airs2834.mat']); a0p15 = x.qc;
  x = load([iPstr '/Chou0.20/airs2834.mat']); a0p20 = x.qc;
  x = load([iPstr '/Chou0.25/airs2834.mat']); a0p25 = x.qc;
  x = load([iPstr '/Chou0.30/airs2834.mat']); a0p30 = x.qc;
  x = load([iPstr '/Chou0.35/airs2834.mat']); a0p35 = x.qc;
  x = load([iPstr '/Chou0.40/airs2834.mat']); a0p40 = x.qc;
  x = load([iPstr '/Chou0.45/airs2834.mat']); a0p45 = x.qc;
  x = load([iPstr '/Chou0.50/airs2834.mat']); a0p50 = x.qc;
  x = load([iPstr '/Chou0.55/airs2834.mat']); a0p55 = x.qc;
  x = load([iPstr '/Chou0.60/airs2834.mat']); a0p60 = x.qc;
  x = load([iPstr '/Chou0.65/airs2834.mat']); a0p65 = x.qc;
  x = load([iPstr '/Chou0.70/airs2834.mat']); a0p70 = x.qc;
  x = load([iPstr '/Chou0.75/airs2834.mat']); a0p75 = x.qc;
  x = load([iPstr '/Chou0.80/airs2834.mat']); a0p80 = x.qc;
  x = load([iPstr '/Chou0.85/airs2834.mat']); a0p85 = x.qc;
  x = load([iPstr '/Chou0.90/airs2834.mat']); a0p90 = x.qc;
  clear x
else
  disp('loading QuickConvolve convolutions')
  x = load([iPstr '/DISORT/conv4440.mat']);   disort = x.qc; f = x.fc;
  x = load([iPstr '/Chou0.00/conv4440.mat']); a0p00 = x.qc;
  x = load([iPstr '/Chou0.05/conv4440.mat']); a0p05 = x.qc;
  x = load([iPstr '/Chou0.10/conv4440.mat']); a0p10 = x.qc;
  x = load([iPstr '/Chou0.15/conv4440.mat']); a0p15 = x.qc;
  x = load([iPstr '/Chou0.20/conv4440.mat']); a0p20 = x.qc;
  x = load([iPstr '/Chou0.25/conv4440.mat']); a0p25 = x.qc;
  x = load([iPstr '/Chou0.30/conv4440.mat']); a0p30 = x.qc;
  x = load([iPstr '/Chou0.35/conv4440.mat']); a0p35 = x.qc;
  x = load([iPstr '/Chou0.40/conv4440.mat']); a0p40 = x.qc;
  x = load([iPstr '/Chou0.45/conv4440.mat']); a0p45 = x.qc;
  x = load([iPstr '/Chou0.50/conv4440.mat']); a0p50 = x.qc;
  x = load([iPstr '/Chou0.55/conv4440.mat']); a0p55 = x.qc;
  x = load([iPstr '/Chou0.60/conv4440.mat']); a0p60 = x.qc;
  x = load([iPstr '/Chou0.65/conv4440.mat']); a0p65 = x.qc;
  x = load([iPstr '/Chou0.70/conv4440.mat']); a0p70 = x.qc;
  x = load([iPstr '/Chou0.75/conv4440.mat']); a0p75 = x.qc;
  x = load([iPstr '/Chou0.80/conv4440.mat']); a0p80 = x.qc;
  x = load([iPstr '/Chou0.85/conv4440.mat']); a0p85 = x.qc;
  x = load([iPstr '/Chou0.90/conv4440.mat']); a0p90 = x.qc;
  clear x
end

tdisort = rad2bt(f,disort);
t0p00 = rad2bt(f,a0p00);
t0p05 = rad2bt(f,a0p05);
t0p10 = rad2bt(f,a0p10);
t0p15 = rad2bt(f,a0p15);
t0p20 = rad2bt(f,a0p20);
t0p25 = rad2bt(f,a0p25);
t0p30 = rad2bt(f,a0p30);
t0p35 = rad2bt(f,a0p35);
t0p40 = rad2bt(f,a0p40);
t0p45 = rad2bt(f,a0p45);
t0p50 = rad2bt(f,a0p50);
t0p55 = rad2bt(f,a0p55);
t0p60 = rad2bt(f,a0p60);
t0p65 = rad2bt(f,a0p65);
t0p70 = rad2bt(f,a0p70);
t0p75 = rad2bt(f,a0p75);
t0p80 = rad2bt(f,a0p80);
t0p85 = rad2bt(f,a0p85);
t0p90 = rad2bt(f,a0p90);

iaChouFac = 0 : 0.05 : 0.90;

for iReg = 1:14
  if iReg == 1
    iaInd = find(f <= 1655);
    disp('f <= 1655')
  elseif iReg == 2
    iaInd = find(f <= 805);
    disp('f <= 805')
  elseif iReg == 3
    iaInd = find(f <= 667 | (f >= 670 & f <= 805));
    disp('f <= 805, avoid 667')
  elseif iReg == 4
    iaInd = find(f >= 667 & f <= 670);
    disp('f = 667')  
  elseif iReg == 5
    iaInd = find(f > 805 & f <= 980);
    disp('805 < f <= 980')
  elseif iReg == 6
    iaInd = find(f > 980 & f <= 1105);
    disp('980 < f <= 1105')
  elseif iReg == 7
    iaInd = find(f > 1105 & f <= 1280);
    disp('1105 < f <= 1280')
  elseif iReg == 8
    iaInd = find(f > 1280 & f <= 1380);
    disp('1280 < f <= 1380')
  elseif iReg == 9
    iaInd = find(f > 1380 & f < 1655);
    disp('1380 < f <= 1655')
  elseif iReg == 10
    iaInd = find(f > 1655 & f < 2005);
    disp('1655 < f <= 2005')
  elseif iReg == 11
    iaInd = find(f > 2005 & f < 2130);
    disp('2005 < f <= 2130')
  elseif iReg == 12
    iaInd = find(f > 2130 & f < 2280);
    disp('2130 < f <= 2280')
  elseif iReg == 13
    iaInd = find(f > 2280 & f < 2380);
    disp('2280 < f <= 2380')
  elseif iReg == 14
    iaInd = find(f > 2380 & f < 2830);
    disp('2380 < f <= 2830')
  end

  chisqr00 = (t0p00(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr00 = chisqr00.^2; chisqr00 = sum(chisqr00,1)/length(f(iaInd));
  chisqr05 = (t0p05(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr05 = chisqr05.^2; chisqr05 = sum(chisqr05,1)/length(f(iaInd));
  chisqr10 = (t0p10(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr10 = chisqr10.^2; chisqr10 = sum(chisqr10,1)/length(f(iaInd));
  chisqr15 = (t0p15(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr15 = chisqr15.^2; chisqr15 = sum(chisqr15,1)/length(f(iaInd));
  chisqr20 = (t0p20(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr20 = chisqr20.^2; chisqr20 = sum(chisqr20,1)/length(f(iaInd));
  chisqr25 = (t0p25(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr25 = chisqr25.^2; chisqr25 = sum(chisqr25,1)/length(f(iaInd));
  chisqr30 = (t0p30(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr30 = chisqr30.^2; chisqr30 = sum(chisqr30,1)/length(f(iaInd));
  chisqr35 = (t0p35(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr35 = chisqr35.^2; chisqr35 = sum(chisqr35,1)/length(f(iaInd));
  chisqr40 = (t0p40(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr40 = chisqr40.^2; chisqr40 = sum(chisqr40,1)/length(f(iaInd));
  chisqr45 = (t0p45(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr45 = chisqr45.^2; chisqr45 = sum(chisqr45,1)/length(f(iaInd));
  chisqr50 = (t0p50(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr50 = chisqr50.^2; chisqr50 = sum(chisqr50,1)/length(f(iaInd));
  chisqr55 = (t0p55(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr55 = chisqr55.^2; chisqr55 = sum(chisqr55,1)/length(f(iaInd));
  chisqr60 = (t0p60(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr60 = chisqr60.^2; chisqr60 = sum(chisqr60,1)/length(f(iaInd));
  chisqr65 = (t0p65(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr65 = chisqr65.^2; chisqr65 = sum(chisqr65,1)/length(f(iaInd));
  chisqr70 = (t0p70(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr70 = chisqr70.^2; chisqr70 = sum(chisqr70,1)/length(f(iaInd));
  chisqr75 = (t0p75(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr75 = chisqr75.^2; chisqr75 = sum(chisqr75,1)/length(f(iaInd));
  chisqr80 = (t0p80(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr80 = chisqr80.^2; chisqr80 = sum(chisqr80,1)/length(f(iaInd));
  chisqr85 = (t0p85(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr85 = chisqr85.^2; chisqr85 = sum(chisqr85,1)/length(f(iaInd));
  chisqr90 = (t0p90(iaInd,:)-tdisort(iaInd,:))./tdisort(iaInd,:); chisqr90 = chisqr90.^2; chisqr90 = sum(chisqr90,1)/length(f(iaInd));
  chisqr = [chisqr00; chisqr05; chisqr10; chisqr15; chisqr20; chisqr25; chisqr30; chisqr35; chisqr40; chisqr45; chisqr50; chisqr55; chisqr60; chisqr65; chisqr70; chisqr75; chisqr80; chisqr85; chisqr90];
  clear tbest
  for ii = 1 : 6048
    junk = chisqr(:,ii);
    meanchisqr(ii) = nanmean(junk);
    stdchisqr(ii) = nanstd(junk);
    junk = find(junk == min(junk),1);
    best(ii) = junk;
    miaow = ['bestspectra = t0p' num2str((junk-1)*5,'%02d') ';'];
    eval(miaow);
    tbest(:,ii) = bestspectra(iaInd,ii);
  end

  booI = find(p.cngwat > 0 & p.ctype == 201); 
  chisqr00 = (t0p00(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr00 = chisqr00.^2; chisqr00 = sum(chisqr00,1)/length(f(iaInd));
  chisqr05 = (t0p05(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr05 = chisqr05.^2; chisqr05 = sum(chisqr05,1)/length(f(iaInd));
  chisqr10 = (t0p10(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr10 = chisqr10.^2; chisqr10 = sum(chisqr10,1)/length(f(iaInd));
  chisqr15 = (t0p15(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr15 = chisqr15.^2; chisqr15 = sum(chisqr15,1)/length(f(iaInd));
  chisqr20 = (t0p20(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr20 = chisqr20.^2; chisqr20 = sum(chisqr20,1)/length(f(iaInd));
  chisqr25 = (t0p25(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr25 = chisqr25.^2; chisqr25 = sum(chisqr25,1)/length(f(iaInd));
  chisqr30 = (t0p30(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr30 = chisqr30.^2; chisqr30 = sum(chisqr30,1)/length(f(iaInd));
  chisqr35 = (t0p35(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr35 = chisqr35.^2; chisqr35 = sum(chisqr35,1)/length(f(iaInd));
  chisqr40 = (t0p40(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr40 = chisqr40.^2; chisqr40 = sum(chisqr40,1)/length(f(iaInd));
  chisqr45 = (t0p45(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr45 = chisqr45.^2; chisqr45 = sum(chisqr45,1)/length(f(iaInd));
  chisqr50 = (t0p50(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr50 = chisqr50.^2; chisqr50 = sum(chisqr50,1)/length(f(iaInd));
  chisqr55 = (t0p55(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr55 = chisqr55.^2; chisqr55 = sum(chisqr55,1)/length(f(iaInd));
  chisqr60 = (t0p60(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr60 = chisqr60.^2; chisqr60 = sum(chisqr60,1)/length(f(iaInd));
  chisqr65 = (t0p65(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr65 = chisqr65.^2; chisqr65 = sum(chisqr65,1)/length(f(iaInd));
  chisqr70 = (t0p70(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr70 = chisqr70.^2; chisqr70 = sum(chisqr70,1)/length(f(iaInd));
  chisqr75 = (t0p75(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr75 = chisqr75.^2; chisqr75 = sum(chisqr75,1)/length(f(iaInd));
  chisqr80 = (t0p80(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr80 = chisqr80.^2; chisqr80 = sum(chisqr80,1)/length(f(iaInd));
  chisqr85 = (t0p85(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr85 = chisqr85.^2; chisqr85 = sum(chisqr85,1)/length(f(iaInd));
  chisqr90 = (t0p90(iaInd,booI)-tdisort(iaInd,booI))./tdisort(iaInd,booI); chisqr90 = chisqr90.^2; chisqr90 = sum(chisqr90,1)/length(f(iaInd));
  chisqr = [chisqr00; chisqr05; chisqr10; chisqr15; chisqr20; chisqr25; chisqr30; chisqr35; chisqr40; chisqr45; chisqr50; chisqr55; chisqr60; chisqr65; chisqr70; chisqr75; chisqr80; chisqr85; chisqr90];
  clear tbestI
  for iii = 1 : length(booI)
    ii = booI(iii);
    junk = chisqr(:,iii);
    meanchisqrI(iii) = nanmean(junk);
    stdchisqrI(iii) = nanstd(junk);
    junk = find(junk == min(junk),1);
    bestI(iii) = junk;
    miaow = ['bestspectra = t0p' num2str((junk-1)*5,'%02d') ';'];
    eval(miaow);
    tbestI(:,iii) = bestspectra(iaInd,ii);
  end

  booW = find(p.cngwat > 0 & p.ctype == 101); 
  chisqr00 = (t0p00(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr00 = chisqr00.^2; chisqr00 = sum(chisqr00,1)/length(f(iaInd));
  chisqr05 = (t0p05(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr05 = chisqr05.^2; chisqr05 = sum(chisqr05,1)/length(f(iaInd));
  chisqr10 = (t0p10(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr10 = chisqr10.^2; chisqr10 = sum(chisqr10,1)/length(f(iaInd));
  chisqr15 = (t0p15(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr15 = chisqr15.^2; chisqr15 = sum(chisqr15,1)/length(f(iaInd));
  chisqr20 = (t0p20(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr20 = chisqr20.^2; chisqr20 = sum(chisqr20,1)/length(f(iaInd));
  chisqr25 = (t0p25(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr25 = chisqr25.^2; chisqr25 = sum(chisqr25,1)/length(f(iaInd));
  chisqr30 = (t0p30(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr30 = chisqr30.^2; chisqr30 = sum(chisqr30,1)/length(f(iaInd));
  chisqr35 = (t0p35(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr35 = chisqr35.^2; chisqr35 = sum(chisqr35,1)/length(f(iaInd));
  chisqr40 = (t0p40(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr40 = chisqr40.^2; chisqr40 = sum(chisqr40,1)/length(f(iaInd));
  chisqr45 = (t0p45(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr45 = chisqr45.^2; chisqr45 = sum(chisqr45,1)/length(f(iaInd));
  chisqr50 = (t0p50(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr50 = chisqr50.^2; chisqr50 = sum(chisqr50,1)/length(f(iaInd));
  chisqr55 = (t0p55(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr55 = chisqr55.^2; chisqr55 = sum(chisqr55,1)/length(f(iaInd));
  chisqr60 = (t0p60(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr60 = chisqr60.^2; chisqr60 = sum(chisqr60,1)/length(f(iaInd));
  chisqr65 = (t0p65(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr65 = chisqr65.^2; chisqr65 = sum(chisqr65,1)/length(f(iaInd));
  chisqr70 = (t0p70(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr70 = chisqr70.^2; chisqr70 = sum(chisqr70,1)/length(f(iaInd));
  chisqr75 = (t0p75(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr75 = chisqr75.^2; chisqr75 = sum(chisqr75,1)/length(f(iaInd)); 
  chisqr80 = (t0p80(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr80 = chisqr80.^2; chisqr80 = sum(chisqr80,1)/length(f(iaInd));
  chisqr85 = (t0p85(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr85 = chisqr85.^2; chisqr85 = sum(chisqr85,1)/length(f(iaInd));
  chisqr90 = (t0p90(iaInd,booW)-tdisort(iaInd,booW))./tdisort(iaInd,booW); chisqr90 = chisqr90.^2; chisqr90 = sum(chisqr90,1)/length(f(iaInd));
  chisqr = [chisqr00; chisqr05; chisqr10; chisqr15; chisqr20; chisqr25; chisqr30; chisqr35; chisqr40; chisqr45; chisqr50; chisqr55; chisqr60; chisqr65; chisqr70; chisqr75; chisqr80; chisqr85; chisqr90];
  clear tbestW
  for iii = 1 : length(booW)
    ii = booW(iii);
    junk = chisqr(:,iii);
    meanchisqrW(iii) = nanmean(junk);
    stdchisqrW(iii) = nanstd(junk);
    junk = find(junk == min(junk),1);
    bestW(iii) = junk;
    miaow = ['bestspectra = t0p' num2str((junk-1)*5,'%02d') ';'];
    eval(miaow);
    tbestW(:,iii) = bestspectra(iaInd,ii);
  end

  figure(9)  
  plot(f(iaInd),nanmean(t0p00(iaInd,:)-tdisort(iaInd,:),2),'b',f(iaInd),nanmean(tbest-tdisort(iaInd,:),2),'r',...
       f(iaInd),nanstd(t0p00(iaInd,:)-tdisort(iaInd,:),0,2),'c--',f(iaInd),nanstd(tbest-tdisort(iaInd,:),0,2),'m--',fairs(g2),nedt(g2),'k'); plotaxis2; title('BOTH')
  xlim([min(f(iaInd)) max(f(iaInd))])
  xlabel('Wavenumber cm-1'); ylabel('BT(K)')
  figure(10)  
  plot(f(iaInd),nanmean(t0p00(iaInd,booI)-tdisort(iaInd,booI),2),'b',f(iaInd),nanmean(tbestI-tdisort(iaInd,booI),2),'r',...
       f(iaInd),nanstd(t0p00(iaInd,booI)-tdisort(iaInd,booI),0,2),'c--',f(iaInd),nanstd(tbestI-tdisort(iaInd,booI),0,2),'m--',fairs(g2),nedt(g2),'k'); plotaxis2; title('ICE')
  xlim([min(f(iaInd)) max(f(iaInd))])
  xlabel('Wavenumber cm-1'); ylabel('BT(K)')
  figure(11)  
  plot(f(iaInd),nanmean(t0p00(iaInd,booW)-tdisort(iaInd,booW),2),'b',f(iaInd),nanmean(tbestW-tdisort(iaInd,booW),2),'r',...
       f(iaInd),nanstd(t0p00(iaInd,booW)-tdisort(iaInd,booW),0,2),'c--',f(iaInd),nanstd(tbestW-tdisort(iaInd,booW),0,2),'m--',fairs(g2),nedt(g2),'k'); plotaxis2; title('WATER')
  xlim([min(f(iaInd)) max(f(iaInd))])
  xlabel('Wavenumber cm-1'); ylabel('BT(K)')

  if iReg == 1
    fprintf(1,'ICE   scanang x cprtop x cngwat x cpsize = %4i x %4i x %4i x %4i \n',length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),  length(unique(p.cpsize(booI))))
    fprintf(1,'WATER scanang x cprtop x cngwat x cpsize = %4i x %4i x %4i x %4i \n',length(unique(p.scanang(booW))),length(unique(p.cprtop(booW))),length(unique(p.cngwat(booW))),  length(unique(p.cpsize(booW))))
    %% make a N dim matrix, one for Ice and one for Water,      Freq x scanang x cprtop x cngwat x cpsize
    %% see ~/KCARTA/TEST/DISORT_vs_PCLSAM/PARAMETRIZE/make_profiles.m
    matrI  = reshape(iaChouFac(bestI),length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI))));
    matrW  = reshape(iaChouFac(bestW),length(unique(p.cpsize(booW))),length(unique(p.cngwat(booW))),length(unique(p.cprtop(booW))),length(unique(p.scanang(booW))));
    indexI = reshape(booI,length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI))));
    indexW = reshape(booW,length(unique(p.cpsize(booW))),length(unique(p.cngwat(booW))),length(unique(p.cprtop(booW))),length(unique(p.scanang(booW))));

    matrI  = reshape(iaChouFac(bestI),length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); matrI = permute(matrI,[ 4 3 2 1]);
    matrW  = reshape(iaChouFac(bestW),length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); matrW = permute(matrW,[ 4 3 2 1]); 
    indexI = reshape(booI,length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); indexI = permute(indexI,[ 4 3 2 1]);
    indexW = reshape(booW,length(unique(p.scanang(booI))),length(unique(p.cprtop(booI))),length(unique(p.cngwat(booI))),length(unique(p.cpsize(booI)))); indexW = permute(indexW,[ 4 3 2 1]);

    cpsizeI = unique(p.cpsize(booI));
    cpsizeW = unique(p.cpsize(booW));
    cprtop = unique(p.cprtop(booI));
    cngwat = unique(p.cngwat(booI));
    scanang = unique(p.scanang(booI));
    comment = 'see reshape(bestI,length(unique(p.cpsize(booI))),length(unique(p.cngwat(booI))),length(unique(p.cprtop(booI))),length(unique(p.scanang(booI)))) and make_profiles.m';
    save generic_605_1655_I_W.mat matrI matrW cpsizeI cpsizeW cprtop cngwat scanang comment indexI indexW
    %% [matrII,matrWW,indexII,indexWW] = write_chou_matfor('generic_605_1655_I_W.mat','generic_605_1655_I_W.bin',605,2830);
  end

  jett = jet(64); jett = jett(32-24:32+24,:);
  
  figure(1); clf
  boo = find(p.cngwat > 0 & p.ctype == 201); 
  scatter(p.cprtop(boo),log10(p.cngwat(boo)),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),log10(p.cngwat(boo)),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('log(cngwat)'); 
  title('Ice'); colormap(jett)
  
  figure(2); clf
  boo = find(p.cngwat > 0 & p.ctype == 101); 
  scatter(p.cprtop(boo),log10(p.cngwat(boo)),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),log10(p.cngwat(boo)),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('log(cngwat)'); 
  title('Water'); colormap(jett)
  
  figure(3); clf
  boo = find(p.cngwat > 0 & p.ctype == 201); 
  scatter(p.cprtop(boo),p.cpsize(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),p.cpsize(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('cpsize'); 
  title('Ice'); colormap(jett)
  
  figure(4); clf
  boo = find(p.cngwat > 0 & p.ctype == 101); 
  scatter(p.cprtop(boo),p.cpsize(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),p.cpsize(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('cpsize'); 
  title('Water'); colormap(jett)
  
  figure(5); clf
  boo = find(p.cngwat > 0 & p.ctype == 201); 
  scatter(p.cprtop(boo),p.scanang(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),p.scanang(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('scanang'); 
  title('Ice'); colormap(jett)
  
  figure(6); clf
  boo = find(p.cngwat > 0 & p.ctype == 101); 
  scatter(p.cprtop(boo),p.scanang(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter(p.cprtop(boo),p.scanang(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; xlabel('cprtop'); ylabel('scanang'); 
  title('Water'); colormap(jett)
  
  addpath /home/sergio/MATLABCODE/PLOTTER/matlab-axis-label-alignment-master
  figure(7); clf
  set(gca, 'dataaspectratio', [1 1 0.5], 'projection', 'perspective', 'box', 'on')
  boo = find(p.cngwat > 0 & p.ctype == 201); 
  scatter3(p.cprtop(boo),log10(p.cngwat(boo)),p.cpsize(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter3(p.cprtop(boo),log10(p.cngwat(boo)),p.cpsize(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar;  colormap jet; 
  xlh = xlabel('cprtop','fontsize',16); ylabel('log(cngwat)','fontsize',16); zlabel('cpsize','fontsize',16); %xtickangle(45)
  title('Ice'); colormap(jett)
  h = rotate3d;
  set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
  set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
  set(gcf, 'ResizeFcn', @align_axislabel)
  align_axislabel([], gca)
  xlh.Position = [30 -1 0];
  
  figure(8); clf
  boo = find(p.cngwat > 0 & p.ctype == 101); 
  scatter3(p.cprtop(boo),log10(p.cngwat(boo)),p.cpsize(boo),50,log10(iaChouFac(best(boo))),'filled'); 
  scatter3(p.cprtop(boo),log10(p.cngwat(boo)),p.cpsize(boo),50,(iaChouFac(best(boo))),'filled'); 
  caxis([min(iaChouFac) max(iaChouFac)]); colorbar; colormap jet; 
  xlh = xlabel('cprtop','fontsize',16); ylabel('log(cngwat)','fontsize',16); zlabel('cpsize','fontsize',16); xtickangle(45)
  title('Water'); colormap(jett)
  h = rotate3d;
  set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
  set(h, 'ActionPostCallback', 'set(gcf,''windowbuttonmotionfcn'','''')')
  set(gcf, 'ResizeFcn', @align_axislabel)
  align_axislabel([], gca)
  xlh.Position = [30 -1 0];

  figure(12)
  %% if you plot DISORT vs PCLSAM with correction CHou -- 0, profile 390 is "biggest bias" at 900 cm-1; this is a water cloud, booW(336) = 390, iaChouFac(bestW(336)) = 0.4
  %% i900 = find(f >= 900,1); plot(tdisort(i900,:)-t0p00(i900,:)); iWorst = find(tdisort(i900,:)-t0p00(i900,:) == min(tdisort(i900,:)-t0p00(i900,:)))           == 390
  %% i1231 = find(f >= 1231,1); plot(tdisort(i1231,:)-t0p00(i1231,:)); iWorst = find(tdisort(i1231,:)-t0p00(i1231,:) == min(tdisort(i1231,:)-t0p00(i1231,:)))   == 408
  %% print_cloud_params(h,p,390)      %% TURNS OUT TO BE WATER CLOUD, so does 408        iaChouFac(bestW(find(booW == 390))) =  0.4,  iaChouFac(bestW(find(booW == 408))) =  0.3
  %% find(booW == 390)
  iWorst = 390;
  iWorst = 408;
  if iWorst == 380
   subplot(211); plot(f,tdisort(:,iWorst),'b',f,t0p00(:,iWorst),'k',f,t0p40(:,iWorst),'r')
    subplot(212); plot(f,tdisort(:,iWorst)-t0p00(:,iWorst),'b',f,tdisort(:,iWorst)-t0p40(:,iWorst),'r')
  elseif iWorst == 408
   subplot(211); plot(f,tdisort(:,iWorst),'b',f,t0p00(:,iWorst),'k',f,t0p30(:,iWorst),'r')
    subplot(212); plot(f,tdisort(:,iWorst)-t0p00(:,iWorst),'b',f,tdisort(:,iWorst)-t0p30(:,iWorst),'r')
  end
  disp('ret to continue'); pause
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bluesky = find(p.cngwat <= 0.01);
bluesky = find(p.cngwat == 0);

both = 1 : length(p.stemp);  %% hig and lo clouds says 0.2 is best
both = find((p.ctype == 201 & p.cprtop <= 500) | (p.ctype == 101 & p.cprtop > 500)); %% logical, high ice/low water cloud says 0.3 is best
%both = find((p.ctype == 201 & p.cprtop <= 500) | (p.ctype == 101 & p.cprtop > 500)); %% upside down, low ice cloud says 0.1 is best

ice = find(p.ctype == 201);
ice = find(p.ctype == 201 & p.cprtop <= 500); %% logical,     high ice cloud says 0.3 is best
%ice = find(p.ctype == 201 & p.cprtop > 500);  %% upside down, low ice cloud says 0.1 is best

water = find(p.ctype == 101);
water = find(p.ctype == 101 & p.cprtop > 500);  %% logical, low water cloud says 0.1 is best
%water = find(p.ctype == 101 & p.cprtop <= 500); %% upside down, high water cloud says 0.3 is best

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(f(1:2378),nanmean(tdisort(1:2378,both)'-t0p00(1:2378,both)'),'bx-',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p00(1:2378,both)'),'bx--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p05(1:2378,both)'),'gx-',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p05(1:2378,both)'),'gx--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p10(1:2378,both)'),'rx-',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p10(1:2378,both)'),'rx--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p15(1:2378,both)'),'kx-',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p15(1:2378,both)'),'kx--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p20(1:2378,both)'),'c',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p20(1:2378,both)'),'c--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p25(1:2378,both)'),'m',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p25(1:2378,both)'),'m--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p30(1:2378,both)'),'y',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p30(1:2378,both)'),'y--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p35(1:2378,both)'),'b',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p35(1:2378,both)'),'b--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p40(1:2378,both)'),'g',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p40(1:2378,both)'),'g--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p45(1:2378,both)'),'r',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p45(1:2378,both)'),'r--',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p50(1:2378,both)'),'k',f(1:2378),nanstd(tdisort(1:2378,both)'-t0p50(1:2378,both)'),'k--')
     grid; axis([640 1640 -1 +1]); title('Both');

figure(2)
plot(f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p00(1:2378,ice)'),'bx-',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p00(1:2378,ice)'),'bx--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p05(1:2378,ice)'),'gx-',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p05(1:2378,ice)'),'gx--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p10(1:2378,ice)'),'rx-',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p10(1:2378,ice)'),'rx--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p15(1:2378,ice)'),'kx-',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p15(1:2378,ice)'),'kx--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p20(1:2378,ice)'),'c',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p20(1:2378,ice)'),'c--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p25(1:2378,ice)'),'m',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p25(1:2378,ice)'),'m--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p30(1:2378,ice)'),'y',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p30(1:2378,ice)'),'y--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p35(1:2378,ice)'),'b',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p35(1:2378,ice)'),'b--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p40(1:2378,ice)'),'g',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p40(1:2378,ice)'),'g--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p45(1:2378,ice)'),'r',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p45(1:2378,ice)'),'r--',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p50(1:2378,ice)'),'k',f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p50(1:2378,ice)'),'k--')
grid; axis([640 1640 -1 +1]); title('Ice')

figure(3)
plot(f(1:2378),nanmean(tdisort(1:2378,water)'-t0p00(1:2378,water)'),'bx-',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p00(1:2378,water)'),'bx--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p05(1:2378,water)'),'gx-',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p05(1:2378,water)'),'gx--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p10(1:2378,water)'),'rx-',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p10(1:2378,water)'),'rx--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p15(1:2378,water)'),'kx-',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p15(1:2378,water)'),'kx--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p20(1:2378,water)'),'c',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p20(1:2378,water)'),'c--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p25(1:2378,water)'),'m',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p25(1:2378,water)'),'m--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p30(1:2378,water)'),'y',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p30(1:2378,water)'),'y--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p35(1:2378,water)'),'b',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p35(1:2378,water)'),'b--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p40(1:2378,water)'),'g',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p40(1:2378,water)'),'g--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p45(1:2378,water)'),'r',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p45(1:2378,water)'),'r--',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p50(1:2378,water)'),'k',f(1:2378),nanstd(tdisort(1:2378,water)'-t0p50(1:2378,water)'),'k--')
grid; axis([640 1640 -1 +1]); title('Water')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

both = 1 : length(p.stemp);
both = find((p.ctype == 201 & p.cprtop <= 500) | (p.ctype == 101 & p.cprtop > 500));
figure(1); 
subplot(211);
plot(f(1:2378),nanmean(tdisort(1:2378,both)'-t0p00(1:2378,both)'),'bx-',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p05(1:2378,both)'),'gx-',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p10(1:2378,both)'),'rx-',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p15(1:2378,both)'),'kx-',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p20(1:2378,both)'),'c',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p25(1:2378,both)'),'m',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p30(1:2378,both)'),'y',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p35(1:2378,both)'),'b',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p40(1:2378,both)'),'g',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p45(1:2378,both)'),'r',...
     f(1:2378),nanmean(tdisort(1:2378,both)'-t0p50(1:2378,both)'),'k')
grid; axis([640 1640 -1 +1]); title('Both'); ylabel('Mean');
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);
subplot(212);
plot(f(1:2378),nanstd(tdisort(1:2378,both)'-t0p00(1:2378,both)'),'bx-',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p05(1:2378,both)'),'gx-',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p10(1:2378,both)'),'rx-',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p15(1:2378,both)'),'kx-',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p20(1:2378,both)'),'c',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p25(1:2378,both)'),'m',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p30(1:2378,both)'),'y',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p35(1:2378,both)'),'b',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p40(1:2378,both)'),'g',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p45(1:2378,both)'),'r',...
     f(1:2378),nanstd(tdisort(1:2378,both)'-t0p50(1:2378,both)'),'k')
grid; axis([640 1640 0 +1]); ylabel('Std')
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);

figure(2); 
subplot(211);
plot(f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p00(1:2378,ice)'),'bx-',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p05(1:2378,ice)'),'gx-',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p10(1:2378,ice)'),'rx-',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p15(1:2378,ice)'),'kx-',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p20(1:2378,ice)'),'c',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p25(1:2378,ice)'),'m',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p30(1:2378,ice)'),'y',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p35(1:2378,ice)'),'b',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p40(1:2378,ice)'),'g',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p45(1:2378,ice)'),'r',...
     f(1:2378),nanmean(tdisort(1:2378,ice)'-t0p50(1:2378,ice)'),'k')
grid; axis([640 1640 -1 +1]); title('Ice'); ylabel('Mean');
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);
subplot(212);
plot(f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p00(1:2378,ice)'),'bx-',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p05(1:2378,ice)'),'gx-',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p10(1:2378,ice)'),'rx-',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p15(1:2378,ice)'),'kx-',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p20(1:2378,ice)'),'c',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p25(1:2378,ice)'),'m',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p30(1:2378,ice)'),'y',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p35(1:2378,ice)'),'b',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p40(1:2378,ice)'),'g',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p45(1:2378,ice)'),'r',...
     f(1:2378),nanstd(tdisort(1:2378,ice)'-t0p50(1:2378,ice)'),'k')
grid; axis([640 1640 0 +1]); ylabel('Std')
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);

figure(3); 
subplot(211);
plot(f(1:2378),nanmean(tdisort(1:2378,water)'-t0p00(1:2378,water)'),'bx-',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p05(1:2378,water)'),'gx-',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p10(1:2378,water)'),'rx-',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p15(1:2378,water)'),'kx-',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p20(1:2378,water)'),'c',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p25(1:2378,water)'),'m',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p30(1:2378,water)'),'y',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p35(1:2378,water)'),'b',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p40(1:2378,water)'),'g',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p45(1:2378,water)'),'r',...
     f(1:2378),nanmean(tdisort(1:2378,water)'-t0p50(1:2378,water)'),'k')
grid; axis([640 1640 -1 +1]); title('Water'); ylabel('Mean');
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);
subplot(212);
plot(f(1:2378),nanstd(tdisort(1:2378,water)'-t0p00(1:2378,water)'),'bx-',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p05(1:2378,water)'),'gx-',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p10(1:2378,water)'),'rx-',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p15(1:2378,water)'),'kx-',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p20(1:2378,water)'),'c',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p25(1:2378,water)'),'m',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p30(1:2378,water)'),'y',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p35(1:2378,water)'),'b',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p40(1:2378,water)'),'g',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p45(1:2378,water)'),'r',...
     f(1:2378),nanstd(tdisort(1:2378,water)'-t0p50(1:2378,water)'),'k')
grid; axis([640 1640 0 +1]); ylabel('Std')
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);

figure(4)
subplot(211);
plot(f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p00(1:2378,bluesky)'),'bx-',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p05(1:2378,bluesky)'),'gx-',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p10(1:2378,bluesky)'),'rx-',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p15(1:2378,bluesky)'),'kx-',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p20(1:2378,bluesky)'),'c',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p25(1:2378,bluesky)'),'m',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p30(1:2378,bluesky)'),'y',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p35(1:2378,bluesky)'),'b',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p40(1:2378,bluesky)'),'g',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p45(1:2378,bluesky)'),'r',...
     f(1:2378),nanmean(tdisort(1:2378,bluesky)'-t0p50(1:2378,bluesky)'),'k')
grid; axis([640 1640 -0.1 +0.1]); title('Bluesky'); ylabel('Mean');
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);
subplot(212);
plot(f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p00(1:2378,bluesky)'),'bx-',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p05(1:2378,bluesky)'),'gx-',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p10(1:2378,bluesky)'),'rx-',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p15(1:2378,bluesky)'),'kx-',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p20(1:2378,bluesky)'),'c',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p25(1:2378,bluesky)'),'m',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p30(1:2378,bluesky)'),'y',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p35(1:2378,bluesky)'),'b',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p40(1:2378,bluesky)'),'g',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p45(1:2378,bluesky)'),'r',...
     f(1:2378),nanstd(tdisort(1:2378,bluesky)'-t0p50(1:2378,bluesky)'),'k')
grid; axis([640 1640 0 +0.2]); ylabel('Std')
plotaxis2; hl = legend(num2str(iaChouFac'),'location','best','fontsize',8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
whos both ice water bluesky
