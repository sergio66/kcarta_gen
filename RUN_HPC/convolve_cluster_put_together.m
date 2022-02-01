clear all

error('kinda junky code, does not have latest CRIS-lo,med,hi or latest IASI')

addpath /asl/matlib/rtptools
addpath /asl/matlib/h4tools
addpath /home/sergio/KCARTA/MATLAB

set_gasOD_cumOD_rad_jac_flux_cloud_lblrtm;

set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %4i profiless \n',use_this_rtp,length(p.stemp))

%iMax = input('Enter number of profiles : ')
iMax = length(p.stemp);

airs_convolve_file_numchans  %% gives latest clist/sfile
iWhichInstr = input('Enter (1) AIRS (4) CRIS HI (14) AIRS=CRIS HI : ');
if iWhichInstr == 1
  strinstr = 'airs';
  leninstr = 35000;
  chinstr = 2378;
  chinstr = length(clist);
elseif iWhichInstr == 4
  strinstr = 'crisHI';
  leninstr = 21500;
  chinstr = 2235;  
elseif iWhichInstr == 14
  strinstr = 'AIRS_crisHI';
  leninstr = 21500;
  chinstr = 2378;  
end

for ii = 1 : iMax
  if mod(ii,100) == 0
    fprintf(1,'looking for convolved mat file %5i of %5i \n',ii,iMax)
  end
  if exist(['JUNK/individual_prof_convolved_kcarta_' strinstr '_' num2str(ii) '.mat'])
    thedir = dir(['JUNK/individual_prof_convolved_kcarta_' strinstr '_' num2str(ii) '.mat']);
    if thedir.bytes > leninstr
      %% correct size
      iaYes(ii) = 2;
    else
      %% not completely done
      iaYes(ii) = 1;
    end
  else
    iaYes(ii) = -1;
  end
end

good = 1 : iMax;
small = find(iaYes == 1);
good = setdiff(good,small);
bad = find(iaYes < 0);
good = setdiff(good,bad);
ohoh = union(small,bad);
if length(ohoh) > 0
  whos ohoh
  fid = fopen('JUNK/notdone.txt','w');
  for ii = 1 : length(ohoh)
    fprintf(fid,'%4i \n',ohoh(ii));
  end
  fclose(fid);
end

good = 1 : iMax;
small = find(iaYes == 1);
good = setdiff(good,small);
if length(small) > 0
  small
  disp('above files incomplete')
  iYes = input('continue (+1)Y (-1)N : ');
  if iYes < 0
    error('stopping')
  end
end

bad = find(iaYes < 0);
good = setdiff(good,bad);
if length(bad) > 0
  bad
  disp('above files not found')
  iYes = input('continue (+1)Y (-1)N : ');
  if iYes < 0
    error('stopping')
  end
end

fprintf(1,'reading in %5i of %5i \n',length(good),length(p.stemp))
if iWhichInstr ~= 14
  dKcall = zeros(chinstr,iMax);
else
  dKcall.airs = zeros(2378,iMax);
  dKcall.crisHi = zeros(2235,iMax);  
end

for ii = 1 : iMax
  if length(intersect(ii,good)) == 1
    thefile = ['JUNK/individual_prof_convolved_kcarta_' strinstr '_' num2str(ii) '.mat'];
    loader = [' a = load(''' thefile ''');'];
    eval(loader)
    wah = a.rKc;
    [mmsize,nnsize] = size(wah);
    if nnsize == 1
      if iWhichInstr == 1
        dKcall(:,ii) = a.rKc';
      elseif iWhichInstr == 4      
        dKcall(:,ii) = a.rcris_all';
      elseif iWhichInstr == 14
        dKcall.airs(:,ii) = a.rKc';    
        dKcall.crisHi(:,ii) = a.rcris_all';
      end
    else
      if ii == 1
        if iWhichInstr ~= 14
          dKcall = zeros(iMax,2378,nnsize);
        else
          dKcall.airs = zeros(iMax,2378,nnsize);
          dKcall.crisHi = zeros(iMax,2235,nnsize);  
        end
      end     
      if iWhichInstr == 1
        dKcall(ii,:,:) = a.rKc;
      elseif iWhichInstr == 4      
        dKcall(ii,:,:) = a.rcris_all;
      elseif iWhichInstr == 14
        dKcall.airs(ii,:,:) = a.rKc;    
        dKcall.crisHi(ii,:,:) = a.rcris_all;
      end
    end    
  end
  if mod(ii,100) == 0
    fprintf(1,'    reading in %4i of %4i \n',ii,iMax);
  end
end

whos good

goodlist = good;
if iWhichInstr == 1
  fKc = a.fKc;
elseif iWhichInstr == 4
  fKc = a.fcris;
elseif iWhichInstr == 14
  fKc.airs = a.fKc;
  fKc.crisHi = a.fcris;
end  
saver = ['save JUNK/xconvolved_kcarta_' strinstr '_put_together.mat goodlist use_this_rtp fKc dKcall kcartaexec'];
eval(saver)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/sergio/MATLABCODE/CLOUD/
if iWhichInstr == 1
  g = dogoodchan;
else
  g = 1 : chinstr;
end

if ~isfield(h,'vchan')
  if iWhichInstr == 1
    h.vchan = instr_chans;
  elseif iWhichInstr == 4
    h.vchan = instr_chans('crishi');
  end
end

if iWhichInstr == 14
  error('too lazy to code this for AIRS n CrIShi')
end

tobs = real(rad2bt(h.vchan,p.robs1));
tkCarta = real(rad2bt(h.vchan,dKcall));

figure(1); clf
if isfield(p,'rcalc')
  tSarta  = rad2bt(h.vchan,p.rcalc);
  plot(h.vchan(g),nanmean(tobs(g,good)'-tSarta(g,good)'),'b', h.vchan(g),nanstd(tobs(g,good)'-tSarta(g,good)'),'c',...
       h.vchan(g),nanmean(tobs(g,good)'-tkCarta(g,good)'),'r',h.vchan(g),nanstd(tobs(g,good)'-tkCarta(g,good)'),'m')
else
  plot(h.vchan(g),nanmean(tobs(g,good)'-tkCarta(g,good)'),'r',h.vchan(g),nanstd(tobs(g,good)'-tkCarta(g,good)'),'m')
end  

figure(2); clf
plot(tobs(1291,good), tobs(1291,good)-tkCarta(1291,good), '.'); xlabel('BT1231 obs'); ylabel('BT1231 obs-cal')
plot(tobs(1291,good), tkCarta(1291,good), '.', tobs(1291,:),tobs(1291,:)); xlabel('BT1231 obs'); ylabel('BT1231 cal')

figure(3); clf
addpath /home/sergio/MATLABCODE/SHOWSTATS/
[n x y nmean nstd] = myhist2d(tobs(1291,good),tobs(1291,good)-tkCarta(1291,good), 200 : 5 : 300, -10 : 2 : +10, 1, 1);
pcolor(x,y,log10(n)); shading flat; colormap jet; hold on; errorbar(200:5:300,nmean,nstd,'color','k','linewidth',2); hold off
jett = jet; jett(1,:) = 1; colormap(jett); colorbar
xlabel('BT1231 obs'); ylabel('BT1231 obs-kcarta')
mesh(x,y,n); shading flat; colormap jet; hold on; errorbar(200:5:300,nmean,nstd,'color','k','linewidth',2); hold off
axis([200 300 -10 +10 0 100])

figure(4); clf
dbt = 200 : 1 : 320;
trp = intersect(good,find(abs(p.rlat) <= 30));
plot(dbt,hist(tobs(1291,good),dbt),'b',dbt,hist(tkCarta(1291,good),dbt),'r',...
     dbt,hist(tobs(1291,trp),dbt),'c',dbt,hist(tkCarta(1291,trp),dbt),'m',...
     'linewidth',2);
ylabel('hist'); xlabel('BT1231'); hl=legend('tobs all','tKc all','tobs trp','tKc trp','location','northeast');
set(hl,'fontsize',10)

figure(4); clf
dbt = -20 : 1 : 120;
trp = intersect(good,find(abs(p.rlat) <= 30 & p.landfrac == 0));
trp = intersect(good,find(abs(p.rlat) <= 30));
plot(dbt,hist(p.stemp(good)-tobs(1291,good),dbt),'b',dbt,hist(p.stemp(good)-tkCarta(1291,good),dbt),'r',...
     dbt,hist(p.stemp(trp)-tobs(1291,trp),dbt),'c',dbt,hist(p.stemp(trp)-tkCarta(1291,trp),dbt),'m',...
     'linewidth',2);
%semilogy(dbt,hist(p.stemp(good)-tobs(1291,good),dbt),'b',dbt,hist(p.stemp(good)-tkCarta(1291,good),dbt),'r',...
%     dbt,hist(p.stemp(trp)-tobs(1291,trp),dbt),'c',dbt,hist(p.stemp(trp)-tkCarta(1291,trp),dbt),'m',...
%     'linewidth',2);
ylabel('hist'); xlabel('stemp-BT1231'); hl=legend('tobs all','tKc all','tobs trp ocean','tKc trp ocean','location','northeast');
set(hl,'fontsize',10)

trp = intersect(good,find(abs(p.rlat) <= 30));
plot(dbt,hist(p.stemp(trp)-tobs(1291,trp),dbt),'k',dbt,hist(p.stemp(trp)-tkCarta(1291,trp),dbt),'c',...
     'linewidth',2);
ylabel('hist'); xlabel('stemp-BT1231'); hl=legend('tobs trp','tKc trp','location','northeast'); grid
set(hl,'fontsize',10)
