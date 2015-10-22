% SUBROUTINE flux_moment_slowloopLinearVaryT in rad_flux.f
% SUBROUTINE FindGauss2(nn,daX,daW)          in rad_misc.f

plevsx = p.plevs(1:p.nlevs);
plevsx = flipud(plevsx);

%% see ~/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB//driver_rrtm_no_xsec_HP.m
%% comment = 'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET FLUX       HEATING RATE';
rrtm_up_flux = sum(squeeze(rrtm_results.heating_rate_info(1,4:16,:,3)));
rrtm_dn_flux = sum(squeeze(rrtm_results.heating_rate_info(1,4:16,:,4)));

for ii = 1 : p.nlevs
  %% upwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;
  fUP(ii,:) = 0;
  fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))  
end
%disp('ret'); pause
pause(1)

daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0 
	  
radsUP = zeros(4,p.nlevs,length(w));
radsDN = zeros(4,p.nlevs,length(w));
fUP   = zeros(p.nlevs,length(w));
fDN = zeros(p.nlevs,length(w));

%[p.plevs(p.nlevs-1) p.spres p.plevs(p.nlevs)]
frac = (p.spres - p.plevs(p.nlevs-1))/(p.plevs(p.nlevs)-p.plevs(p.nlevs-1))
%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : p.nlevs
  %% upwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;
  fUP(ii,:) = 0;
  if ii == 1
    fprintf(1,'>> ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))
  else
    fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))
  end
  for jj = 1 : 4
    if ii == 1
      radsUP(jj,ii,:) = ttorad(w,p.stemp);
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);
    elseif ii == 2
      layfrac = frac;
      od = d(:,laykc)*layfrac;
      mu = daX(jj);
      planckrad = ttorad(w,p.ptemp(play));
      junk = squeeze(radsUP(jj,ii-1,:)).*exp(-od/mu) + planckrad.*(1-exp(-od/mu));
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    else
      layfrac = 1;
      od = d(:,laykc)*layfrac;
      mu = daX(jj);      
      planckrad = ttorad(w,p.ptemp(play));
      junk = squeeze(radsUP(jj,ii-1,:)).*exp(-od/mu) + planckrad.*(1-exp(-od/mu));
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    end
  end
end

figure(1);
kcUP = flux(:,1:86);
sum_kcUP = sum(kcUP,1)*0.0025/1000;
sum_fUP  = sum(fUP,2)*0.0025/1000*2*pi;
semilogy(sum_kcUP,plevsx,'b.-',sum_fUP,plevsx,'r.-',rrtm_up_flux,rrtm_results.plevs,'k.-'); grid
  hl = legend('kCARTA','matlab','RRTM');
  set(gca,'ydir','reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')

for ii = p.nlevs : -1 : 1
  %% downwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;  
  fDN(ii,:) = 0;
  fprintf(1,'ii- laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))    
  for jj = 1 : 4
    if ii == p.nlevs
      radsDN(jj,ii,:) = ttorad(w,2.6);
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);
    elseif ii == 1
      layfrac = frac;    
      od = d(:,laykc)*layfrac;
      mu = daX(jj);
      planckrad = ttorad(w,p.ptemp(play));
      junk = squeeze(radsDN(jj,ii+1,:)).*exp(-od/mu) + planckrad.*(1-exp(-od/mu));
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    else
      layfrac = 1;
      od = d(:,laykc)*layfrac;
      mu = daX(jj);      
      planckrad = ttorad(w,p.ptemp(play));
      junk = squeeze(radsDN(jj,ii+1,:)).*exp(-od/mu) + planckrad.*(1-exp(-od/mu));
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    end
  end
end

figure(2)
kcDN = flux(:,86+(1:86));
sum_kcDN = sum(kcDN,1)*0.0025/1000;
sum_fDN  = sum(fDN,2)*0.0025/1000*2*pi;
semilogy(sum_kcDN,plevsx,'b.-',sum_fDN,plevsx,'r.-',rrtm_dn_flux,rrtm_results.plevs,'k.-'); grid
  hl = legend('kCARTA','matlab','RRTM');
  set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plevs = p.plevs(1:p.nlevs-1);
plevs = p.plevs(2:p.nlevs);
plevs = flipud(plevs);

net_kc = sum_kcUP - sum_kcDN;
net_ml = sum_fUP  - sum_fDN;

heat_kc = net_kc(1:85); heat_kc = diff(heat_kc) ./ diff(plevs') * 8.4391; 
heat_ml = net_ml(1:85); heat_ml = diff(heat_ml') ./ diff(plevs') * 8.4391;

hgt = p2h(plevs)/1000;
hgt = 0.5*(hgt(2:end)+hgt(1:end-1));
figure(3);
plot(heat_kc,hgt,'bx-',...
    sum(hr,1)*0.0025,p2h(plevs)/1000,'c',...
    heat_ml,hgt,'ro-',...
    rrtm_results.kc605to2830_hr,p2h(rrtm_results.plevs)/1000,'k',...
    'linewidth',2);
hl = legend('kCARTA flux->hr','raw kCARTA HR','matlab flux->hr','RRTM','location','southwest'); grid


quick_rrtm_surf_flux   = ttorad(630:0.0025:3250,p.stemp)*pi/1000*0.0025; 
quick_kcarta_surf_flux = ttorad(605:0.0025:2830,p.stemp)*pi/1000*0.0025;
junk = [max(rrtm_up_flux) sum(quick_rrtm_surf_flux) max(sum_kcUP) sum(quick_kcarta_surf_flux)];
disp('SURFACE UP FLUXES')
disp('RRTM : 630:3250           KC : 605:2830');
disp('actual  estimate      actual    estimate')
fprintf(1,'%8.6f  %8.6f  %8.6f  %8.6f \n',junk);