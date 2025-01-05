addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE/
addpath /home/sergio/MATLABCODE/TIME
addpath /home/sergio/MATLABCODE/PLOTTER
addpath /home/sergio/MATLABCODE/SHOWSTATS

if ~exist('all_the_fovs.mat')
  iCnt = 0;
  addpath /home/sergio/KCARTA/MATLAB
  for ii = 1:1:26000
    f0 = ['Masuda/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
    f1 = ['Default/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
    f2 = ['ScanangV2_yesresetSTWV/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat'];
  
    if exist(f0) & exist(f1) & exist(f2)
      iCnt = iCnt + 1;
      iaFind(iCnt) = ii;
      a0 = load(f0);
      a1 = load(f1);
      a2 = load(f2);
      raaMasuda(:,iCnt)  = a0.rKc;
      raaDefault(:,iCnt) = a1.rKc;
      raaScanang(:,iCnt) = a2.rKc;
    end
  
    if mod(iCnt,100) == 0
      fprintf(1,'iCnt = %5i read %5i kcarta run \n',iCnt,ii)
      fairs = a1.fKc;
      t0 = rad2bt(fairs,raaMasuda);
      t1 = rad2bt(fairs,raaDefault);
      t2 = rad2bt(fairs,raaScanang);
      plot(fairs,nanmean(t1'-t2'),'r',fairs,nanstd(t1'-t2'),'m',fairs,nanmean(t1'-t0'),'b',fairs,nanstd(t1'-t0'),'c'); 
      xlim([640 1440]); ylim([-0.25 +0.25])
      pause(0.1)
    end
  end
  
  save all_the_fovs.mat fairs raaDefault raaScanang raaMasuda iaFind
else
  disp('loading the kCARTA data');
  load all_the_fovs.mat
end

load /home/sergio/MATLABCODE/CRODGERS_FAST_CLOUD/sarta_chans_for_l1c.mat
t0 = rad2bt(fairs(ichan),raaMasuda(ichan,:));
t1 = rad2bt(fairs(ichan),raaDefault(ichan,:));
t2 = rad2bt(fairs(ichan),raaScanang(ichan,:));

addpath /asl/matlib/h4tools
iX = +1;  %% nalli emis
iX = -1;  %% masuda emis

iS = -1; %% scott sarta
iS = +1; %% chris sarta

if iX > 0
  %if iS < 0
    [hhh,hhha,pppS,pppa] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junk.rp.rtp');
    load indices_orig_rtp.mat
    disp('scott sarta, new emis')
  %elseif iS > 0
    [hhh,hhha,pppC,pppa] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junk.rp.rtp_chrisHsarta');
    load indices_orig_rtp.mat
    disp('chris sarta, new emis')
  %end
elseif iX < 0
  %if iS < 0
    [hhh,hhha,pppS,pppa] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junkMasuda.rp.rtp');
    disp('scott sarta, masuda emis')
  %elseif iS > 0
    [hhh,hhha,pppC,pppa] = rtpread('/home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junkMasuda.rp.rtp_chrisHsarta');
    disp('chris sarta, masuda emis')
  %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tobsS = rad2bt(hhh.vchan,pppS.robs1);
tcalS = rad2bt(hhh.vchan,pppS.rcalc);
tobsC = rad2bt(hhh.vchan,pppC.robs1);
tcalC = rad2bt(hhh.vchan,pppC.rcalc);

nansum(nansum(tobsS-tobsC))
figure(1)
plot(hhh.vchan,nanmean(tcalC'-tcalS'),'k'); title('Chris-Scott')

figure(2);
plot(hhh.vchan,nanmean(tobsS'-tcalS'),'b',hhh.vchan,nanmean(tobsS'-tcalC'),'r',hhh.vchan,nanmean(tcalC'-tcalS'),'k')
hl = legend('Obs-Scott2008','Obs-Chris2018','Chris2018-Scott2008','location','best');

tobs = tobsC;
tcal = tcalC;

whos tobs* tcal* raa*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(0,'DefaultLegendAutoUpdate','off')
addpath /home/sergio/MATLABCODE/PLOTTER

plot(hhh.vchan,nanmean(tobs'-tcal'),'k',hhh.vchan,nanmean(tobs'-t0'),'g',...
     hhh.vchan,nanmean(tobs'-t1'),'b',hhh.vchan,nanmean(tobs'-t2'),'r')
title('Bias Obs-Cal'); xlabel('Wavenumber cm-1');
hl = legend('SARTA','kCARTA Masuda','kCARTA Emis Default','kCARTA Emis Scanang','location','best');
set(hl,'fontsize',10)

plotaxis2;
xlim([640 1440]); ylim([-0.25 +0.25])
xlim([640 1440]); ylim([-2 +2])

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp');
[yy,mm,dd,hh] = tai2utcSergio(p.rtime);

x22 = 1 : length(p.satzen);
x22 = find(abs(p.satzen-22) < 0.5);
x22 = find(p.rlat > -50); %% avoid S. Pole wince this is July 2019 so maybe ice

figure(2)
scatter_coast(p.rlon(x22),p.rlat(x22),10,p.stemp(x22))
colormap jet

figure(3);
i820 = find(hhh.vchan >= 820,1);
i900 = find(hhh.vchan >= 900,1);
i961 = find(hhh.vchan >= 961,1);
i1231 = find(hhh.vchan >= 1231,1);

[X,I] = sort(p.stemp(x22));
plot(p.stemp(x22(I)),tobs(i820,x22(I))-t2(i820,x22(I)),'r',p.stemp(x22(I)),tobs(i900,x22(I))-t2(i900,x22(I)),'b',...
     p.stemp(x22(I)),tobs(i961,x22(I))-t2(i961,x22(I)),'g',p.stemp(x22(I)),tobs(i1231,x22(I))-t2(i1231,x22(I)),'k')
hl = legend('820','900','961','1231'); ylabel('Bias (K)'); xlabel('stemp')
title('Nick Nalli Model')

plot(p.stemp(x22(I)),tobs(i820,x22(I))-t0(i820,x22(I)),'r',p.stemp(x22(I)),tobs(i900,x22(I))-t0(i900,x22(I)),'b',...
     p.stemp(x22(I)),tobs(i961,x22(I))-t0(i961,x22(I)),'g',p.stemp(x22(I)),tobs(i1231,x22(I))-t0(i1231,x22(I)),'k')
hl = legend('820','900','961','1231'); ylabel('Bias (K)'); xlabel('stemp')
title('Masuda Model')

plot(p.stemp(x22(I)),tobs(i820,x22(I))-t2(i820,x22(I)),'r',p.stemp(x22(I)),tobs(i820,x22(I))-t0(i820,x22(I)),'b')
hl = legend('Nick Nalli','Masuda'); ylabel('820 cm-1 Bias (K)'); xlabel('stemp')

%%%%%%%%%%%%%%%%%%%%%%%%%
[n,nx,ny,nickmean,nickstd] = myhist2d(p.stemp(x22(I)),tobs(i820,x22(I))-t2(i820,x22(I)),270:1:320,-4:0.1:+4);
[n,nx,ny,masudamean,masudastd] = myhist2d(p.stemp(x22(I)),tobs(i820,x22(I))-t0(i820,x22(I)),270:1:320,-4:0.1:+4);
errorbar((270:1:320)-0.25,nickmean,nickstd,'color','r','linewidth',2); hold on
errorbar((270:1:320)+0.25,masudamean,masudastd,'color','b','linewidth',2); hold off
plotaxis2; 
hl = legend('Nick Nalli','Masuda'); ylabel('820 cm-1 Bias (K)'); xlabel('stemp')

subplot(211); plot((270:1:320),nickmean,'r',(270:1:320),masudamean,'b'); ylabel('Bias');
title('820 cm-1 stats (b) Masuda (r) NNalli'); 
plotaxis2;
subplot(212); plot((270:1:320),nickstd,'r',(270:1:320),masudastd,'b'); ylabel('Std'); xlabel('Stemp');

%%%%%%%%%%%%%%%%%%%%%%%%%
[n,nx,ny,nickmean,nickstd] = myhist2d(p.stemp(x22(I)),tobs(i900,x22(I))-t2(i900,x22(I)),270:1:320,-4:0.1:+4);
[n,nx,ny,masudamean,masudastd] = myhist2d(p.stemp(x22(I)),tobs(i900,x22(I))-t0(i900,x22(I)),270:1:320,-4:0.1:+4);
errorbar((270:1:320)-0.25,nickmean,nickstd,'color','r','linewidth',2); hold on
errorbar((270:1:320)+0.25,masudamean,masudastd,'color','b','linewidth',2); hold off
plotaxis2; 
hl = legend('Nick Nalli','Masuda'); ylabel('900 cm-1 Bias (K)'); xlabel('stemp')

subplot(211); plot((270:1:320),nickmean,'r',(270:1:320),masudamean,'b'); ylabel('Bias');
title('900 cm-1 stats (b) Masuda (r) NNalli'); 
plotaxis2;
subplot(212); plot((270:1:320),nickstd,'r',(270:1:320),masudastd,'b'); ylabel('Std'); xlabel('Stemp');

clf
dbt = 270 : 320; plot(dbt,histc(p.stemp,dbt))

i856 = find(hhh.vchan >= 856.375,1);
i863 = find(hhh.vchan >= 863.125,1);
i939 = find(hhh.vchan >= 939.250,1);
i962 = find(hhh.vchan >= 962.500,1);
ddbt0 = tobs([i820 i856 i863 i900 i939 i962],:) - t0([i820 i856 i863 i900 i939 i962],:);
ddbt2 = tobs([i820 i856 i863 i900 i939 i962],:) - t2([i820 i856 i863 i900 i939 i962],:);
plot(p.stemp,ddbt0(1,:)-ddbt0(6,:),'b.',p.stemp,ddbt2(1,:)-ddbt0(6,:),'r.')
[n,nx,ny,nmean0,nstd0] = myhist2d(p.stemp(x22),ddbt0(1,x22)-ddbt0(6,x22),270:1:320,-4:0.1:+4);
[n,nx,ny,nmean2,nstd2] = myhist2d(p.stemp(x22),ddbt2(1,x22)-ddbt2(6,x22),270:1:320,-4:0.1:+4);
errorbar(270:1:320,nmean2,nstd2,'color','r'); hold on
errorbar(270:1:320,nmean0,nstd0,'color','b'); hold off
title('double diff : obs-cal : 820-960 cm-1'); hl = legend('Nick Nalli','Masuda');
plotaxis2; 
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(hhh.vchan,nanmean(tobs'-tcal'),'k',hhh.vchan,nanmean(tobs'-t0'),'g')
title('Bias Obs-Cal'); xlabel('Wavenumber cm-1');
hl = legend('SARTA Masuda','kCARTA Masuda','location','best');
set(hl,'fontsize',10)

addpath /home/sergio/MATLABCODE/PLOTTER
plotaxis2;
xlim([640 1440]); ylim([-0.25 +0.25])
xlim([640 1440]); ylim([-2 +2])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(hhh.vchan,nanmean(tobs'-tcal') - nanmean(tobs'-t0'),'g')
title('Double DIff Bias Obs-Cal : SARTA-kCARTA '); xlabel('Wavenumber cm-1');
plotaxis2;
xlim([640 1440]); ylim([-1 +1])

