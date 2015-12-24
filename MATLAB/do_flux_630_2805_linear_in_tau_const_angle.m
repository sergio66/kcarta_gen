function [fUP,fDN,radsUP,radsDN] = do_flux_630_2805_linear_in_tau_const_angle(w,d,p,iVers);

% SUBROUTINE flux_moment_slowloopLinearVaryT in rad_flux.f
% SUBROUTINE FindGauss2(nn,daX,daW)          in rad_misc.f

% input
%  p = one profile structure
%  d = ODs
%  rrtm = rrtm fluxes/heating rates
%  iVers= -1 (constant T),4 (linear tau,O(tau^2)) 41 (Pade) 42 (linear in tau, O(tau))
%
% also look at /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB/flux_rrtm_padeVSkvcarta.m

%{
example usage :
take an   rtp profile, have p.upwell = +1, dump out ODS, fluxes, 4 upwelling rads, read them in
take same rtp profile, have p.upwell = +2, dump out      fluxes, 4 dnwelling rads, read them in

[h,ha,p,pa] = rtpread('test_night_unit_emiss_uplook.rtp');  %% p.upwell = 2
[h,ha,p,pa] = rtpread('test_night_unit_emiss.rtp');         %% p.upwell = 1

[ods,w]     = readkcstd('junky.dat');   %% dump out ODs at all layers
[radKCDN,w] = readkcstd('junky.dat');   %% dump out DNwell rads at all levels
[radKCUP,w] = readkcstd('junky.dat');   %% dump out UPwell rads at all levels

[rad0,w]    = readkcstd('junky.dat');       %% dump out TOA rad
[flux0,w]   = readkcflux('junky.dat_ALL');  %% dump out up/dn fluxes

[fUP,fDN,radsUP,radsDN] = do_flux_630_2805_linear_in_tau_const_angle(w,ods,p,42);

>> whos rads*
  Name             Size                      Bytes  Class     Attributes

  radKCDN      30000x490                 117600000  double   % dnwell rads, f77
  radKCUP      30000x490                 117600000  double   % upwell rads, f77
  radsDN          4x99x30000             95040000  double   % dnwell rads, this code
  radsUP          4x99x30000             95040000  double   % upwell rads, this code


nlays = 98;
nlays = 85;

nlevs = nlays+1;

figure(1); clf; 
%% downwelling radiation
rs = radKCDN(:,(1:nlays)+0*nlays);    %% from f77 kcarta
ru = squeeze(radsDN(1,:,:));  %% from this matlab code
lay = nlays-6:nlays-1; plot(w,rad2bt(w,rs(:,lay))-rad2bt(w,ru(nlevs-lay,:)')); title('downwell')

figure(2); clf
%% upwelling radiation
rs = radKCUP(:,(1:nlays)+0*nlays);
ru = squeeze(radsUP(1,:,:));  %% from this matlab code
lay = nlays-6:nlays-1; plot(w,rad2bt(w,rs(:,lay))-rad2bt(w,ru(lay+1,:)'),'.-'); title('upwell')

figure(1); clf; plot(w,flux0(:,(1:nlevs)+0*nlevs)-fUP');
figure(2); clf; plot(w,flux0(:,(1:nlevs)+1*nlevs)-fDN');

dxx = 0.0025/1000;

palts = p.palts(1:p.nlevs); palts = flipud(palts)/1000;

figure(1); clf
plot(sum(flux0(:,(1:nlevs)+0*nlevs),1)*dxx,palts,'bx-',sum(fUP,2)*dxx,palts,'r',...
     -sum(flux0(:,(1:nlevs)+1*nlevs),1)*dxx,palts,'cx-',-sum(fDN,2)*dxx,palts,'m','linewidth',2)
grid on
hl = legend('f77 kcarta up','matlab up','f77 kcarta dn','matlab dn','location','northwest');
set(hl,'fontsize',10)
xlabel('Flux W/m2'); ylabel('hgt km')

figure(2); clf
plot(sum(flux0(:,(1:nlevs)+0*nlevs),1)*dxx - sum(fUP,2)'*dxx,palts,'b',...
     sum(flux0(:,(1:nlevs)+1*nlevs),1)*dxx - sum(fDN,2)'*dxx,palts,'r','linewidth',2)
grid on
hl = legend('f77 kcarta - matlab up','f77 kcarta - matlab dn','location','northwest');
set(hl,'fontsize',10)
xlabel('\Delta Flux W/m2'); ylabel('hgt km')

%%% for Chris H, RTM flux
%dlist = dir('/asl/s1/chepplew/data/rfm/od_lay_zen_equ_15um_a_*.asc');
dlist = dir('/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/FromChrisH_Dec2015/od_lay_zen_equ_15um_a_*.asc');
nfils = length(dlist);
fprintf(1,'Found %d files\n',nfils);
clear nu fz fzmn;
fzall = [];
for ifn = 1:nfils
  fin = ['/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/FromChrisH_Dec2015/' dlist(ifn).name];
  FH = fopen(fin,'r');
  for i=1:3
     hdr1{i} = fgetl(FH);
  end

  junk   = fgetl(FH);                           % field specs
  params = textscan(junk,'%f %f %f %f %s');
  npnts  = params{1}(1);
  nu     = [params{2}(1):params{3}(1):params{4}(1)];     % {2} wn start, {3} wn step, {4} wn end.
  fz     = fscanf(FH,'%e',[Inf]);
  fzall = [fzall fz];
  fclose(FH);
  fzmn(ifn) = sum(fz) * params{3};
  %clear fz junk;
end
whos fzmn

alts = [0:1:89];
figure(1);clf;semilogx(fzmn,alts,'.-');grid on;
xlabel('Opt Depth');ylabel('Height km');title('Layer OD 630-700wn zen.mix.ctm.h12');
%saveas(gcf,'./figs/lay_od_zen_630_700wn_zen_mix_ctm_h12.png','png')

plot(1:100,ods(1,:),1:90,fzall(1,:))
semilogy(1:85,ods(1,16:100),'bo-',1:90,fzall(1,:),'r','linewidth',2); hl=legend('kcarta','rrtm');
figure(1); semilogy(w,sum(ods,2),nu,sum(fzall,2),'linewidth',2); hl=legend('kcarta','rrtm');
figure(2); allKC = interp1(w,sum(ods,2),nu); plot(nu,sum(fzall,2) ./ allKC'); title('RRTM/KCARTA')

wah = sum(ods,1); plot(wah(12:100),p.palts); set(gca,'ydir','reverse')

%}

plevsx = p.plevs(1:p.nlevs);
plevsx = flipud(plevsx);

%% see ~/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB//driver_rrtm_no_xsec_HP.m
%% comment = 'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET FLUX       HEATING RATE';
%rrtm_up_flux = sum(squeeze(rrtm_results.heating_rate_info(1,4:16,:,3)));
%rrtm_dn_flux = sum(squeeze(rrtm_results.heating_rate_info(1,4:16,:,4)));

for ii = 1 : p.nlevs
  %% upwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;
  fUP(ii,:) = 0;
  %fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6f \n',ii,laykc,play,p.ptemp(play),d(300,laykc))  
end
%disp('ret'); pause
pause(1)

daX = [0.1397599 0.4164096 0.7231570 0.9428958 1.00000]; %% we really need DOUBLE these since we also need -daX
daW = [0.0311810 0.1298475 0.2034646 0.1355069 0.00000]; %% or sum(daW) = 0.5, we need 1.0 

daX = [4.70941630 1.69338507 1.09719858]; daX = 1./daX;
daW = [0.0698269799 0.2292411064 0.2009319137];
						      
daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0 

woo = find(w >= 630-0.0025/2);
woo = woo';
radsUP = zeros(length(daX),p.nlevs,length(woo));
radsDN = zeros(length(daX),p.nlevs,length(woo));
fUP   = zeros(p.nlevs,length(woo));
fDN = zeros(p.nlevs,length(woo));

%[p.plevs(p.nlevs-1) p.spres p.plevs(p.nlevs)]
frac = (p.spres - p.plevs(p.nlevs-1))/(p.plevs(p.nlevs)-p.plevs(p.nlevs-1))
%%%%%%%%%%%%%%%%%%%%%%%%%

%% assume surface emiss = 1
playsN = p.plevs(1:end-1)-p.plevs(2:end);
playsD = log(p.plevs(1:end-1)./p.plevs(2:end));
plays = playsN./playsD;

vt1x = interp1(log(1013.25*plays(1:p.nlevs-1)),p.ptemp(1:p.nlevs-1),log(1013.25*p.plevs(1:p.nlevs)),[],'extrap'); %% this is linear
vt1 = interp1(log(plays(1:p.nlevs-1)),p.ptemp(1:p.nlevs-1),log(p.plevs(1:p.nlevs)),'spline','extrap');  %% Get_Temp_Plevs default is spline
vt1(1:3) = vt1x(1:3); %% use linear inter for first (TOA), as that if what f77 kcarta does

%ugh2 = load('ugh2');
%zaza = ugh2(:,3);
%zaza = [zaza' [204.2082 190.539]];
%plot(ugh2(:,3)-flipud(vt1(3:99)))
%plot(zaza-fliplr((vt1')))
%keyboard

for ii = 1 : p.nlevs
  %% upwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1;
  play  = (p.nlevs-1) - (ii-1)+1;
  fUP(ii,:) = 0;
  if ii == 1
    fprintf(1,'>> ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.ptemp(play),d(300,laykc))
  else
    fprintf(1,'ii+ laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.ptemp(play),d(300,laykc))
  end
  for jj = 1 : length(daX)
    if ii == 1
      %% surface term
      radsUP(jj,ii,:) = ttorad(w(woo),p.stemp);
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);
    elseif ii == 2
      %% through first partial layer
      layfrac = frac;
      od = d(woo,laykc)*layfrac;
      mu = daX(jj);
      junk = find_rad_upwell(w(woo),squeeze(radsUP(jj,ii-1,:)),od,mu,play,p.ptemp,vt1,p.nlevs,iVers);
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    else
      %% all full layers
      layfrac = 1;
      od = d(woo,laykc)*layfrac;
      mu = daX(jj);      
      junk = find_rad_upwell(w(woo),squeeze(radsUP(jj,ii-1,:)),od,mu,play,p.ptemp,vt1,p.nlevs,iVers);      
      radsUP(jj,ii,:) = junk;
      fUP(ii,:) = fUP(ii,:) + squeeze(radsUP(jj,ii,:))'*daW(jj);      
    end
  end
end

fUP = fUP * 2 * pi;    %% pi is integral over azimuth, while "2" is because sum(daW) = 0.5 ie need angles on either side of nadir

fDN = zeros(size(fUP));

%{
figure(1);
kcUP = flux(woo,1:86);
sum_kcUP = sum(kcUP,1)*0.0025/1000;
sum_fUP  = sum(fUP,2)*0.0025/1000*2*pi;
semilogy(sum_kcUP,plevsx,'b.-',sum_fUP,plevsx,'r.-',rrtm_up_flux,rrtm_results.plevs,'k.-'); grid
  hl = legend('kCARTA','matlab','RRTM');
  set(gca,'ydir','reverse')
%}
  
%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')

for ii = p.nlevs : -1 : 1
  %% downwelling flux
  laykc = 101-(p.nlevs-1) + (ii-1)-1+1;
  play  = (p.nlevs-1) - (ii-1);  
  fDN(ii,:) = 0;
  if ii < p.nlevs 
    fprintf(1,'ii- laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,laykc,play,p.ptemp(play),d(300,laykc))
  else
    fprintf(1,'ii- laykc play t(z) od(z) = %3i %3i %3i %8.6f %8.6e \n',ii,101,101,9999,0);
  end
  for jj = 1 : length(daX)
    if ii == p.nlevs
      radsDN(jj,ii,:) = ttorad(w(woo),2.6);
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);
    elseif ii == 1
      layfrac = frac;    
      od = d(woo,laykc)*layfrac;
      mu = daX(jj);
      junk = find_rad_dnwell(w(woo),squeeze(radsUP(jj,ii+1,:)),od,mu,play,p.ptemp,vt1,p.nlevs,iVers);      
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    else
      layfrac = 1;
      od = d(woo,laykc)*layfrac;
      mu = daX(jj);
      junk = find_rad_dnwell(w(woo),squeeze(radsDN(jj,ii+1,:)),od,mu,play,p.ptemp,vt1,p.nlevs,iVers);      
      radsDN(jj,ii,:) = junk;
      fDN(ii,:) = fDN(ii,:) + squeeze(radsDN(jj,ii,:))'*daW(jj);      
    end
  end
end

fDN = fDN * 2 * pi;    %% pi is integral over azimuth, while "2" is because sum(daW) = 0.5 ie need angles on either side of nadir

disp('returning here')
return

figure(2)
kcDN = flux(woo,86+(1:86));
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

% rrtm.f
% HEATFAC = 1.0E-7*(GRAV * SECDY)/(CPDAIR)
GRAV = 9.80665E+02;  CPDAIR = 1.00464; SECDY = 8.64E+04;
HEATFAC = 8.4391;
HEATFAC = 1.0E-7*(GRAV * SECDY)/(CPDAIR);

heat_kc = net_kc(1:85); heat_kc = diff(heat_kc) ./ diff(plevs') * HEATFAC; 
heat_ml = net_ml(1:85); heat_ml = diff(heat_ml') ./ diff(plevs') * HEATFAC;

hgt = p2h(plevs)/1000;
hgt = 0.5*(hgt(2:end)+hgt(1:end-1));
figure(3); clf
plot(heat_kc,hgt,'bx-',...
    sum(hr(woo,:),1)*0.0025,p2h(plevs)/1000,'c',...
    heat_ml,hgt,'ro-',...
    rrtm_results.kc605to2830_hr,p2h(rrtm_results.plevs)/1000,'k',...
    'linewidth',2);
hl = legend('kCARTA flux->hr','raw kCARTA HR','matlab flux->hr','RRTM','location','southwest'); grid

quick_rrtm_surf_flux   = ttorad(630:0.0025:3250,p.stemp)*pi/1000*0.0025; 
quick_kcarta_surf_flux = ttorad(630:0.0025:2830,p.stemp)*pi/1000*0.0025;
junk = [max(rrtm_up_flux) sum(quick_rrtm_surf_flux) max(sum_kcUP) sum(quick_kcarta_surf_flux)];
disp('SURFACE UP FLUXES')
disp('RRTM : 630:3250           KC : 605:2830');
disp('actual  estimate      actual    estimate')
fprintf(1,'%8.6f  %8.6f  %8.6f  %8.6f \n',junk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hgtx = p2h(plevsx)/1000;

figure(4);
kcUP = flux(woo,1:86);
sum_kcUP = sum(kcUP,1)*0.0025/1000;
sum_fUP  = sum(fUP,2)*0.0025/1000*2*pi;
junk = interp1(rrtm_results.plevs,rrtm_up_flux,plevsx)';
%semilogy(sum_kcUP-junk,plevsx,'b.-',sum_fUP-junk',plevsx,'r.-');
plot(sum_kcUP-junk,hgtx,'b.-',sum_fUP-junk',hgtx,'r.-'); 

hold on

figure(4)
kcDN = flux(woo,86+(1:86));
sum_kcDN = sum(kcDN,1)*0.0025/1000;
sum_fDN  = sum(fDN,2)*0.0025/1000*2*pi;
junk = interp1(rrtm_results.plevs,rrtm_dn_flux,plevsx)';
%semilogy(sum_kcDN-junk,plevsx,'c.-',sum_fDN-junk',plevsx,'m.-'); grid;
plot(sum_kcDN-junk,hgtx,'c.-',sum_fDN-junk',hgtx,'m.-'); grid;

  hl = legend('up kCARTA-RRTM','up matlab-RRTM','dn kCARTA-RRTM','dn matlab-RRTM');
  %set(gca,'ydir','reverse');
  grid on

hold off

%%%%%%%%%%%%%%%%%%%%%%%%%
