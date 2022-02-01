addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD

clear all

clear_list = 1 : 7377;
clear_list = load('/home/sergio/PCRTM_XIANGLEI/RUN/clear7377.txt');
[h,ha,p,pa] = rtpread('/home/sergio/PCRTM_XIANGLEI/RUN/junk_ECM7377_notunclr.rp.rtp');

fairs = instr_chans;

load /asl/s1/sergio/data/ITOVS_2016/forITOVS_ECM.mat
%         rad_allsky: [2378x7377 double]
%         rad_clrsky: [2378x7377 double]
%         sarta_clear: [2378x7377 single]
%         sarta_cloud: [2378x7377 single]
%         rcalc: [2378x7377 double]
sum(sum(p0.rcalc-p0.sarta_clear))
sum(sum(p0.rcalc-p0.sarta_cloud))
sum(sum(p0.rcalc-p0.rad_clrsky))
sum(sum(p0.rcalc-p0.rad_allsky))
sum(sum(p.rcalc-p0.sarta_clear))

clear dall
iCnt = 0;
for iix = 1 : length(clear_list)
  ii = clear_list(iix);
  fname = ['JUNK/rad.dat' num2str(ii)];
  if exist(fname)
    [d,w] = readkcstd(fname);
    iCnt = iCnt + 1;
    thefile(iCnt) = ii;    
  else
    d = zeros(1,89000);
  end
  dall(iCnt,:) = d;
  if mod(iCnt,100) == 0
    fprintf(1,'%4i out of %4i clear .. found %4i so far \n',ii,length(clear_list),iCnt)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off', 'MATLAB:imagesci:hdf:removalWarningHDFSD')
[wc,rc] = convolve_airs(w,dall,1:2378);

tKC = rad2bt(wc,rc);
tS  = rad2bt(wc,p.rcalc(:,clear_list));
%% tS  = rad2bt(wc,p0.rcalc(:,clear_list));   %% sum(sum(p0.rcalc-p0.rad_allsky)) <<<<<<<<
tAIRS = real(rad2bt(wc,p.robs1(:,clear_list)));

[mean(tAIRS(1291,:)-tKC(1291,:)) std(tAIRS(1291,:)-tKC(1291,:)) mean(tAIRS(1291,:)-tS(1291,:)) std(tAIRS(1291,:)-tS(1291,:)) mean(tS(1291,:)-tKC(1291,:)) std(tS(1291,:)-tKC(1291,:))]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = dogoodchan;
figure(1); clf
boon = find(p.solzen(clear_list) > 90); %% night
plot(wc(g),nanmean(tS(g,boon)'-tAIRS(g,boon)'),'b',wc(g),nanmean(tKC(g,boon)'-tAIRS(g,boon)'),'r',wc(g),nanmean(tKC(g,boon)'-tS(g,boon)'),'k',...
     wc(g),nanstd(tS(g,boon)'-tAIRS(g,boon)'),'b--',wc(g),nanstd(tKC(g,boon)'-tAIRS(g,boon)'),'r--',wc(g),nanstd(tKC(g,boon)'-tS(g,boon)'),'k--','linewidth',2)
axis([650 2780 -2 +2]); grid
hl = legend('SARTA-AIRS','KCARTA-AIRS','kCARTA-SARTA'); set(hl,'fontsize',10);
title('night')

figure(2); clf
bood = find(p.solzen(clear_list) < 90); %% day
plot(wc(g),nanmean(tS(g,bood)'-tAIRS(g,bood)'),'b',wc(g),nanmean(tKC(g,bood)'-tAIRS(g,bood)'),'r',wc(g),nanmean(tKC(g,bood)'-tS(g,bood)'),'k',...
     wc(g),nanstd(tS(g,bood)'-tAIRS(g,bood)'),'b--',wc(g),nanstd(tKC(g,bood)'-tAIRS(g,bood)'),'r--',wc(g),nanstd(tKC(g,bood)'-tS(g,bood)'),'k--','linewidth',2)
axis([650 2780 -2 +2]); grid
hl = legend('SARTA-AIRS','KCARTA-AIRS','kCARTA-SARTA'); set(hl,'fontsize',10);
title('day')

figure(3); clf
dnboo = 1 : iCnt;
plot(wc(g),nanmean(tS(g,dnboo)'-tAIRS(g,dnboo)'),'b',wc(g),nanmean(tKC(g,dnboo)'-tAIRS(g,dnboo)'),'r',wc(g),nanmean(tKC(g,dnboo)'-tS(g,dnboo)'),'k',...
     wc(g),nanstd(tS(g,dnboo)'-tAIRS(g,dnboo)'),'b--',wc(g),nanstd(tKC(g,dnboo)'-tAIRS(g,dnboo)'),'r--',wc(g),nanstd(tKC(g,dnboo)'-tS(g,dnboo)'),'k--','linewidth',2)
axis([650 2780 -2 +2]); grid
hl = legend('SARTA-AIRS','KCARTA-AIRS','kCARTA-SARTA'); set(hl,'fontsize',10);
title('all')


%{
% saver = ['save JUNK/sarta_kcarta7377_H2012.mat wc rc clear_list'];
saver = ['save JUNK/sarta_kcarta7377_H2012_lblrtm.mat wc rc clear_list'];
% saver = ['save JUNK/sarta_kcarta7377_H2008.mat wc rc clear_list'];
eval(saver)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% george simple test

%bt1231 = btemp(1231,airs_obs(1291,:))';
%bt1227 = btemp(1227,airs_obs(1284,:))';
bt1231 = rad2bt(fairs(1291),p.robs1(1291,:));
bt1227 = rad2bt(fairs(1284),p.robs1(1284,:));

%% This test assumes an emissivity of 0.98 and 30 degree scan angle.
q3 = bt1231-bt1227;  sst1231r5 = bt1231+0.2806+1.2008.*q3+0.2962.*q3.^2+1.0489/cos(30/57.3);

% find reasonable clear night data
vf1n=find(abs(sst1231r5-double(p.stemp))<1 & p.solzen>90);
% find reasonable clear day data
vf1d=find(abs(sst1231r5-double(p.stemp))<1 & p.solzen<90);

commonN = intersect(vf1n,clear_list(boon));
commonD = intersect(vf1d,clear_list(bood));
whos vf1* boo* common*


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
