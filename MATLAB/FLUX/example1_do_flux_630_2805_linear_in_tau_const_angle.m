% example usage :
% take an   rtp profile, have p.upwell = +1, dump out ODS, fluxes, 4 upwelling rads, read them in
% take same rtp profile, have p.upwell = +2, dump out      fluxes, 4 dnwelling rads, read them in

[h,ha,p,pa] = rtpread('test_night_unit_emiss_uplook.rtp');  %% p.upwell = 2
[h,ha,p,pa] = rtpread('test_night_unit_emiss.rtp');         %% p.upwell = 1

[ods,w]     = readkcstd('junky.dat');   %% dump out ODs at all 101 layers
[radKCDN,w] = readkcstd('junky.dat');   %% dump out DNwell rads at all 101 levels
[radKCUP,w] = readkcstd('junky.dat');   %% dump out UPwell rads at all 101 levels

[rad0,w]    = readkcstd('junky.dat');       %% dump out TOA rad
[flux0,w]   = readkcflux('junky.dat_ALL');  %% dump out up/dn fluxes

[fUP,fDN,radsUP,radsDN] = do_flux_630_2805_linear_in_tau_const_angle(w,ods,p,42);

plot(w,sum(ods,2),w,sum(ods101,2));
plot(w,sum(ods,2)./sum(ods101,2)); title('f77 kcarta 1 km vs 101 layers')

nlays = 98;
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

