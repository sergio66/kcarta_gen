function [fUP,fDN,heatRate] = do_flux_630_2805_linear_in_tau_const_angle(w,d,p,iVers,radsum);

% SUBROUTINE flux_moment_slowloopLinearVaryT in rad_flux.f
% SUBROUTINE FindGauss2(nn,daX,daW)          in rad_misc.f

% everything needs to be ordered so layer 1 = TOA and layerN/levelN+1 = GND
%
% input
%  w = wavenumber
%  p = one profile structure, contains p.nlevs and
%      p.plevs ordered from TOA to GND (ie leve1 = TOA, level N is GND)
%      p.tlevs
%      p.plays
%      p.tlays
%      also spres and stemp
%  d = matrix of ODs, size (num wavenumbers, numlays), with columns numbered 1(TOA) to N(gnd)
%  iVers= -1 (constant T),4 (linear tau,O(tau^2)) 41 (Pade) 42 (linear in tau, O(tau))
%  [radsum] = [optional] radsum fluxes/heating rates, Nlevels x 6 rows [level,p,Fup,Fdn,NetF,Heat]
%
% also look at /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB/flux_rrtm_padeVSkvcarta.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(p,'nlevs')
  error('need levels info');
end
if ~isfield(p,'plevs')
  error('need plevels info');
end

nlevs = p.nlevs;

plevsx = p.plevs(1:p.nlevs);
if plevsx(1) < plevsx(2)
  plevsx = flipud(plevsx);
end

mn_sergio = p.nlevs;
if nargin == 5
  [mn_radsum,rows_radsum] = size(radsum);
  if mn_radsum ~= p.nlevs
    error('p.nlevs is different than size of radsum input')
  end
end
mn = mn_sergio;

daX = [0.1397599 0.4164096 0.7231570 0.9428958 1.00000]; %% we really need DOUBLE these since we also need -daX
daW = [0.0311810 0.1298475 0.2034646 0.1355069 0.00000]; %% or sum(daW) = 0.5, we need 1.0 
						      
daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0 

daX = 1.5; daX = 1./daX;
daW = 0.5;

daX = [4.70941630 1.69338507 1.09719858]; daX = 1./daX;
daW = [0.0698269799 0.2292411064 0.2009319137];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radsUP = zeros(length(daX),p.nlevs,length(w));
radsDN = zeros(length(daX),p.nlevs,length(w));
fUP   = zeros(p.nlevs,length(w));
fDN = zeros(p.nlevs,length(w));

for ii = 1 : p.nlevs
  %% upwelling flux
  laykc = mn-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;
  fUP(ii,:) = 0;
  %fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))  
end

[p.plevs(p.nlevs-1) p.spres p.plevs(p.nlevs)];
frac = (p.spres - p.plevs(p.nlevs-1))/(p.plevs(p.nlevs)-p.plevs(p.nlevs-1));

%%%%%%%%%%%%%%%%%%%%%%%%%

%% assume surface emiss = 1
if ~isfield(p,'plays')
  playsN = p.plevs(1:end-1)-p.plevs(2:end);
  playsD = log(p.plevs(1:end-1)./p.plevs(2:end));
  plays = playsN./playsD;
else
  plays = p.plays;
end

if ~isfield(p,'tlevs')
  vt1x = interp1(log(1013.25*plays(1:p.nlevs-1)),p.ptemp(1:p.nlevs-1),log(1013.25*p.plevs(1:p.nlevs)),[],'extrap'); %% this is linear
  vt1 = interp1(log(plays(1:p.nlevs-1)),p.ptemp(1:p.nlevs-1),log(p.plevs(1:p.nlevs)),'spline','extrap');  %% Get_Temp_Plevs default is spline
  vt1(1:3) = vt1x(1:3); %% use linear inter for first (TOA), as that if what f77 kcarta does
else
  vt1 = p.tlevs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : p.nlevs
  %% upwelling flux
  play  = (p.nlevs-1) - (ii-1)+1;
  laykc = play;
  fUP(ii,:) = 0;
  if ii == 1
    fprintf(1,'>> ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.tlays(play),9999)
  else
    fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.tlays(play),d(1,nlevs-ii+1))
  end
  for jj = 1 : length(daX)
    if ii == 1
      %% surface term
      radsUP(jj,ii,:) = ttorad(w,p.stemp);
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);
    elseif ii == 2
      %% through first partial layer
      layfrac = frac;
      od = d(:,nlevs-ii+1)*layfrac;
      mu = daX(jj);
      junk = find_rad_upwell(w,squeeze(radsUP(jj,ii-1,:)),od,mu,play,p.tlays,vt1,p.nlevs,iVers);
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    else
      %% all full layers
      layfrac = 1;
      od = d(:,nlevs-ii+1)*layfrac;
      mu = daX(jj);      
      junk = find_rad_upwell(w,squeeze(radsUP(jj,ii-1,:)),od,mu,play,p.tlays,vt1,p.nlevs,iVers);      
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    end
  end
end

fUP = fUP * 2 * pi;    %% pi is integral over azimuth, while "2" is because sum(daW) = 0.5 ie need angles on either side of nadir

%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')

fDN = zeros(size(fUP));

for ii = p.nlevs : -1 : 1
  %% downwelling flux
  laykc = (nlevs-1) - ii +1;
  play  = laykc;
  fDN(ii,:) = 0;
  if ii < p.nlevs 
    fprintf(1,'ii- laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.tlays(play),d(1,laykc))
  else
    fprintf(1,'ii- laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,101,101,9999,0);
  end
  for jj = 1 : length(daX)
    if ii == p.nlevs
      radsDN(jj,ii,:) = ttorad(w,2.6);
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);
    elseif ii == 1
      layfrac = frac;    
      od = d(:,laykc)*layfrac;
      mu = daX(jj);
      junk = find_rad_dnwell(w,squeeze(radsUP(jj,ii+1,:)),od,mu,play,p.tlays,vt1,p.nlevs,iVers);      
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    else
      layfrac = 1;
      od = d(:,laykc)*layfrac;
      mu = daX(jj);
      junk = find_rad_dnwell(w,squeeze(radsDN(jj,ii+1,:)),od,mu,play,p.tlays,vt1,p.nlevs,iVers);      
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    end
  end
end

fDN = fDN * 2 * pi;    %% pi is integral over azimuth, while "2" is because sum(daW) = 0.5 ie need angles on either side of nadir

delta_w = mean(diff(w));

sum_fDN = sum(fDN,2)*delta_w/1000;
sum_fUP = sum(fUP,2)*delta_w/1000;
figure(1); clf
semilogy(sum_fDN,plevsx,'b.-',sum_fUP,plevsx,'r.-'); grid
  hl = legend('DN','UP');
  set(gca,'ydir','reverse')

if nargin == 5
  %sum(plevsx-radsum(:,2))
  %[sum_fDN radsum(:,4)]
  
  figure(1); clf
  semilogy(sum_fDN,plevsx,'b.-',sum_fUP,plevsx,'r.-',radsum(:,4),radsum(:,2),'c.-',radsum(:,3),radsum(:,2),'m.-'); grid
  hl = legend('OD DN','OD UP','RADSUM DN','RADSUM UP');
  set(gca,'ydir','reverse')

  figure(2); clf
  semilogy(sum_fDN - flipud(radsum(:,4)),plevsx,'b.-',sum_fUP - flipud(radsum(:,3)),plevsx,'m.-'); grid
  hl = legend('OD-RADSUM DN','OD-RADSUM UP');
  set(gca,'ydir','reverse')
  axis([-0.25 +0.25 0 1000]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
net_ml = sum_fUP  - sum_fDN;

% radsum.f
% HEATFAC = 1.0E-7*(GRAV * SECDY)/(CPDAIR)
GRAV = 9.80665E+02;  CPDAIR = 1.00464; SECDY = 8.64E+04;
HEATFAC = 8.4391;
HEATFAC = 1.0E-7*(GRAV * SECDY)/(CPDAIR);

heat_ml = diff(net_ml') ./ diff(plevsx') * HEATFAC;
heatRate = heat_ml;

hgt = p2h(plevsx)/1000;
%hgt = 0.5*(hgt(2:end)+hgt(1:end-1));
figure(3); clf

if nargin == 5
  plot([heat_ml 0],hgt,'b.-',flipud(radsum(:,6)),hgt,'r.-',[heat_ml 0] - (flipud(radsum(:,6)))',hgt,'k.-','linewidth',2);
  hl = legend('matlab flux->hr','RADSUM','matlab-RADSUM','location','southwest'); grid
else
  plot([heat_ml 0],hgt,'b.-','linewidth',2);
  hl = legend('matlab flux->hr','location','southwest'); grid
end
