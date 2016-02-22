% example usage :
% take an   rtp profile, have p.upwell = +1, dump out ODS, fluxes, 4 upwelling rads, read them in
% take same rtp profile, have p.upwell = +2, dump out      fluxes, 4 dnwelling rads, read them in

[h,ha,p,pa] = rtpread('test_night_unit_emiss_uplook.rtp');  %% p.upwell = 2
[h,ha,p,pa] = rtpread('test_night_unit_emiss.rtp');         %% p.upwell = 1

[ods101,w]     = readkcstd('junky.dat');   %% dump out ODs at all usual AIRS 101 layers  <<<<<<<<<

[ods,w]     = readkcstd('junky.dat');   %% dump out ODs at all 1 km layers
[radKCDN,w] = readkcstd('junky.dat');   %% dump out DNwell rads at all 1 km levels
[radKCUP,w] = readkcstd('junky.dat');   %% dump out UPwell rads at all 1 km levels

[rad0,w]    = readkcstd('junky.dat');       %% dump out TOA rad
[flux0,w]   = readkcflux('junky.dat_ALL');  %% dump out up/dn fluxes

[fUP,fDN,radsUP,radsDN] = do_flux_630_2805_linear_in_tau_const_angle(w,ods,p,42);
[fUP5,fDN5,radsU5P,radsDN5] = do_flux_630_2805_linear_in_tau_const_angle(w,ods*5,p,42);

plot(w,sum(ods,2),w,sum(ods101,2));
plot(w,sum(ods,2)./sum(ods101,2)); title('f77 kcarta 1 km vs 101 layers')

%nlays = 98;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

[h,ha,p2,pa] = rtpread('/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/test.newklayers.ip.rtp');
px.ptemp = flipud(p2.ptemp(1:90));
px.plevs = flipud(p2.plevs(1:90));
px.nlevs = 90;
px.spres = p2.spres;
px.stemp = p2.stemp;
[xfUP,xfDN,xradsUP,xradsDN] = do_flux_630_2805_linear_in_tau_const_angle(nu,fzall,px,42);
plot(sum(fUP,2)*mean(diff(nu))/1000,1:90,-sum(fDN,2)*mean(diff(nu))/1000,1:90,'linewidth',2); grid; title('Quick Matlab calc')
xlabel('Flux W/m2'); ylabel('level number')

plot(w,squeeze(radsUP(4,86,:)),nu,squeeze(xradsUP(4,90,:)))
plot(w,rad2bt(w,squeeze(radsUP(4,86,:))),nu,rad2bt(nu,squeeze(xradsUP(4,90,:))))

figure(1); clf
plot(sum(fUP,2)*mean(diff(w))/1000,1:86,'bo-',sum(xfUP,2)*mean(diff(nu))/1000,1:90,'rx-',...
     sum(fDN,2)*mean(diff(w))/1000,1:86,'co-',sum(xfDN,2)*mean(diff(nu))/1000,1:90,'mx-',...
     'linewidth',2); grid
hl = legend('sergio upwell','rtm upwell','sergio dnwell','rtm dnwell');

figure(2); clf
plot(sum(fUP5,2)*mean(diff(w))/1000,1:86,'bo-',sum(xfUP,2)*mean(diff(nu))/1000,1:90,'rx-',...
     sum(fDN5,2)*mean(diff(w))/1000,1:86,'co-',sum(xfDN,2)*mean(diff(nu))/1000,1:90,'mx-',...
     'linewidth',2); grid
hl = legend('sergioX5 upwell','rtm upwell','sergioX5 dnwell','rtm dnwell');

figure(3); clf
pkc = diff(p.plevs(1:p.nlevs)); pkc = flipud(pkc);
pX  = [diff(px.plevs)];         pX  = flipud(pX);
net1 = sum(fUP,2) - sum(fDN,2);   net1 = diff(net1)./pkc * mean(diff(w))/1000 * -8.4391;
net5 = sum(fUP5,2) - sum(fDN5,2); net5 = diff(net5)./pkc * mean(diff(w))/1000 * -8.4391;
netR = sum(xfUP,2) - sum(xfDN,2); netR = diff(netR)./pX   * mean(diff(nu))/1000 * -8.4391;
plot(net1,1:85,'b',net5,1:85,'g',netR,1:89,'r','linewidth',2); title('pseudo heating rate')
hl = legend('kcX1','kcX5','RTM','location','southwest'); grid on; 
axis([-20 +5 0 100]);

