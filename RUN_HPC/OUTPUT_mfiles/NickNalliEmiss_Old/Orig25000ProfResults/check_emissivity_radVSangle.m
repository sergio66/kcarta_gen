addpath /home/sergio/KCARTA/MATLAB
addpath /home/sergio/MATLABCODE

%% I think this is what we SHOULD use since Scott developed reset code using Masuda
finM = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis.rtp';
foutM = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis_resetTWV.rtp';

%% hmm maybe we should do it this way after all since the reset code is agnostic about emissivity ....
finN  = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';
foutN = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';

[h,ha,pM,pa] = rtpread(finM);
[h,ha,pN,pa] = rtpread(finN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,w] = readkcstd('testup.dat');
[fc,qc] = quickconvolve(w,d,0.5,0.5);
theta = 0:10:60;
theta = [[0:10:80] 85];
figure(1); plot(fc,rad2bt(fc,qc)); legend(num2str(theta'));

figure(2); tc = rad2bt(fc,qc); i900 = find(fc >= 900,1); plot(theta,tc(i900,:),'o-');
ax = axis;
acos35 =  acos(3/5)*180/pi;
yy = interp1(theta,tc(i900,:),acos35);
line([0 60],[yy yy],'color','r'); line([acos35 acos35],[ax(3) ax(4)],'color','r');
meansatzen =  mean(pM.satzen);
yy = interp1(theta,tc(i900,:),meansatzen);
line([0 60],[yy yy],'color','b'); line([meansatzen meansatzen],[ax(3) ax(4)],'color','b');

acos35 =  acos(3/5)*180/pi;
yy35 = interp1(theta,tc(i900,:),acos35);
meansatzen =  mean(pM.satzen);
yys = interp1(theta,tc(i900,:),meansatzen);
figure(2); tc = rad2bt(fc,qc); i900 = find(fc >= 900,1); plot(theta,tc(i900,:),'ko-');
hold on; plot(acos35,yy35,'rx',meansatzen,yys,'bo','linewidth',2,'markersize',10); hold off
ylabel('BT900 cm-1'); xlabel('satzen'); title('Uplook Instr')

acos35 =  acos(3/5)*180/pi;
yy35 = interp1(theta,qc(i900,:),acos35);
meansatzen =  mean(pM.satzen);
yys = interp1(theta,qc(i900,:),meansatzen);
figure(2); tc = rad2bt(fc,qc); i900 = find(fc >= 900,1); plot(theta,qc(i900,:),'ko-');
hold on; plot(acos35,yy35,'rx',meansatzen,yys,'bo','linewidth',2,'markersize',10); hold off
ylabel('Rad 900 cm-1'); xlabel('satzen'); title('Uplook Instr')

plot(theta,sin(2*theta*pi/180))
%% integral [mu d(mu) r(mu)] = integral [cos(x) d(cos(x)) r(x) dx] = integral [0.5 sin(2x) r(x) dx]
integrand = sin(2*theta*pi/180).*qc(i900,:);  
yy = interp1(theta,integrand,[meansatzen acos35]);
plot(theta,sin(2*theta*pi/180).*qc(i900,:),'b',[meansatzen acos35],yy,'rx')
xlabel('theta'); ylabel('integrand'); grid; title('integrand = [mu d(mu) r(mu)]')

disp('ret to continue to ecompare KCARTA calcs : masuda-nalli vs angles'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load all_the_fovsV1_correctnnalli_codekcarta.mat

i900 = find(fairs >= 900,1);
plot(pN.satzen,raaScanang(i900,:)-raaMasuda(i900,:),'.')
[n,nx,ny,nMmean,nMstd] = myhist2d(pN.satzen,raaScanang(i900,:)-raaMasuda(i900,:),0:5:70,-1:0.1:1); %% 900 cm-1
plot(0:5:70,nMmean,'o-'); xlabel('satzen'); ylabel('kCARTA dRad'); title('NNalli-Masuda calc at 900 cm-1')
plotaxis2;
disp('ret to continue to emiss plots'); pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fN = pN.efreq(1:30,1);
fM = pM.efreq(1:19,1);

intersect(fN,fM)

addpath /home/sergio/MATLABCODE/SHOWSTATS

ii = 5; %% 900 cm-1, typical nalli > masud
ii = 4; %% 860 cm-1, typical nalli < masuda
figure(1)
dx = 0 : 2 : 60;
[n,nx,ny,nMmean,nMstd] = myhist2d(pM.satzen,pM.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
[n,nx,ny,nNmean,nNstd] = myhist2d(pN.satzen,pN.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
plot(dx,nMmean,'b.-',dx,nNmean,'rx-'); 
  hl = legend('Masuda','Nalli','location','best'); title([num2str(fN(ii)) ' emissivity']); xlabel('satzen')

figure(2)
dx = 0 : 1 : 20;
[n,nx,ny,nMmean,nMstd] = myhist2d(pM.wspeed,pM.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
[n,nx,ny,nNmean,nNstd] = myhist2d(pN.wspeed,pN.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
plot(dx,nMmean,'b.-',dx,nNmean,'rx-'); 
  hl = legend('Masuda','Nalli','location','best'); title([num2str(fN(ii)) ' emissivity']); xlabel('wspeed')

for ii = 1 : 12
  figure(1)
  dx = 0 : 2 : 60;
  [n,nx,ny,nMmean,nMstd] = myhist2d(pM.satzen,pM.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
  [n,nx,ny,nNmean,nNstd] = myhist2d(pN.satzen,pN.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
  plot(dx,nMmean,'b.-',dx,nNmean,'rx-'); 
  hl = legend('Masuda','Nalli','location','best'); title([num2str(fN(ii)) ' emissivity']); xlabel('satzen')
  
  figure(2)
  dx = 0 : 1 : 20;
  [n,nx,ny,nMmean,nMstd] = myhist2d(pM.wspeed,pM.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
  [n,nx,ny,nNmean,nNstd] = myhist2d(pN.wspeed,pN.emis(ii,:),dx,0.8:0.01:1.0); %% 900 cm-1
  plot(dx,nMmean,'b.-',dx,nNmean,'rx-'); 
  hl = legend('Masuda','Nalli','location','best'); title([num2str(fN(ii)) ' emissivity']); xlabel('wspeed')

  pause(0.1);
end
  
