%[sergio@c24-52 GENERIC_RADSnJACS_MANYPROFILES]$ ls -lt JUNK/individual_prof_convolved_kcarta_airs_2223*.mat
%-rw-rw-r-- 1 sergio pi_sergio 8457689 May 28 23:27 JUNK/individual_prof_convolved_kcarta_airs_2223_jac.mat
%-rw-rw-r-- 1 sergio pi_sergio   43622 May 28 23:27 JUNK/individual_prof_convolved_kcarta_airs_2223.mat

N = 4608;

iRegOrCol_Jac = +1; %% full 97 layer jacs for G_N,T,WGT
iRegOrCol_Jac = -1; %% col jacs G_N,T,ST

addpath /home/sergio/git/matlabcode

for ii = 1 : N
  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% these are the .mat files already made
  frad = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
  if iRegOrCol_Jac > 0
    fjac = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];
  else
    fjac = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat'];
  end
    
  thedir = dir(frad);
  if length(thedir) > 0
    iaSizeRad(ii) = thedir.bytes;
  else
    iaSizeRad(ii) = 0;
  end
  
  thedir = dir(fjac);
  if length(thedir) > 0  
    iaSizeJac(ii) = thedir.bytes;
  else
    iaSizeJac(ii) = 0;    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  %% these rad.dat* or jac.dat* files exist so something has maybe gone wrong
  badKCrad = ['JUNK/rad.dat' num2str(ii)];
  badKCjac = ['JUNK/jac.dat' num2str(ii)];  
    
  thedir = dir(badKCrad);
  if length(thedir) > 0
    iaSizeKCRad(ii) = thedir.bytes;
  else
    iaSizeKCRad(ii) = 0;
  end
  
  thedir = dir(badKCjac);
  if length(thedir) > 0  
    iaSizeKCJac(ii) = thedir.bytes;
  else
    iaSizeKCJac(ii) = 0;    
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  if mod(ii,100) == 0
    plot(1:ii,iaSizeRad/(max(iaSizeRad)+1),'b',1:ii,iaSizeJac/(max(iaSizeJac)+1),'r')
    pause(0.1);
  end
end  

figure(1);     plot(1:N,iaSizeRad/(max(iaSizeRad)+1),'b',1:ii,iaSizeJac/(max(iaSizeJac)+1),'r')
figure(2);     plot(1:N,iaSizeKCRad/(max(iaSizeKCRad)+1),'b',1:ii,iaSizeKCJac/(max(iaSizeKCJac)+1),'r')

badRad = find(iaSizeRad == 0);
badJac = find(iaSizeJac == 0);
badBoth = find(iaSizeRad == 0 | iaSizeJac == 0);
badOne   = find(iaSizeRad > 0  & iaSizeJac == 0);  %% if this happened, then the rad conv file should be deleted so we try again
badOneX  = find(iaSizeRad == 0  & iaSizeJac > 0);  %% this should not exist since the rad convolutions are done before jac convolutions
whos badBoth badRad badJac badOne badOneX

if length(intersect(badOne,badBoth)) > 0  %%%% <<<< these are the ind*rad that should be deleted >>>
  disp('delete these mat files where only RAD conv was done')
  printarray(intersect(badOne,badBoth))
  iaX = intersect(badOne,badBoth);
  for iii = 1 : length(iaX)
    ii = iaX(iii);
    frad = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat'];
    if iRegOrCol_Jac > 0
      fjac = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat'];
    else
      fjac = ['JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat'];
    end
    lser = ['!ls -lt  ' frad ' ' fjac]; eval(lser);
    rmer = ['!/bin/rm ' frad ' ' fjac]; eval(rmer);    
  end
end

%badKCRad = find(iaSizeKCRad == 0);
%badKCJac = find(iaSizeKCJac == 0);
%badKCBoth = find(iaSizeKCRad == 0 | iaSizeNBDJac == 0);
%whos badKCBoth badKCRad badKCJac

baddy = badBoth;
if length(baddy) > 0
  write_out_jobsnotdone_for_cluster(baddy,1:N,-1);
else
  disp('YAY looks like all JOBS are done .. but you need to check for corrupted files (incomplete)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
TURNS OUT KCARTA fails on a number of profiles, since the emissivity of first 5 freqs hinges is NaN over land (Africa, Asia, the Americas)
BUT LUCKILY ONLY first 5 hingepoimts, which are almost all high alt channels

%% fix bad NAN for 4608 profiles
[h,ha,p,pa] = rtpread('RTP/summary_23years_all_lat_all_lon_2002_2025_monthlyERA5.op.rtp');
figure(3); plot(p.efreq(:,baddy),p.emis(:,baddy))
[badX,badY] = find(isnan(p.emis));
unique(badY) - baddy'
figure(4); clf; scatter_coast(p.rlon(baddy),p.rlat(baddy),50,p.stemp(baddy)); colormap jet; %% all over land - Africa, Asia, N/S America
unique(badX)  %%% first five channels [1 2 3 4 5]
p.efreq(unique(badX),baddy(1))   %% 645   673   701   729   757    so these are almost  all high alt channels

%% plot(p.stemp(baddy(1)) - rad2bt(h.vchan,p.rcalc(:,baddy(1))))   BAH no sarta
p.efreq(6,1)    %%% 785

pnew = p;
for iii = 1 : length(baddy)
  ii = baddy(iii);
  pnew.emis(1:5,ii) = p.emis(6,ii) * ones(5,1);
  pnew.rho(1:5,ii)  = p.rho(6,ii) * ones(5,1);
end
[badX2,badY2] = find(isnan(pnew.emis)); whos badX2 badY2  YAY

%}
