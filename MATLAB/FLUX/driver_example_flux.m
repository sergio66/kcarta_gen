addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/science
addpath /asl/matlib/aslutil
addpath /home/sergio/KCARTA/MATLAB

cd /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES
iFov    = 1;
fop     = 'RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/test_1kmMIPAS.op.rtp';
od_file   = 'RFM_FLUX/ods.dat';

%% constant layer temp
flux_file = 'RFM_FLUX/rad.dat1_ALL_1kmMIPAS';
hr_file   = 'RFM_FLUX/rad.dat1_HEAT_1kmMIPAS';

%% linear in tau layer temp
flux_file = 'RFM_FLUX/rad.dat1_ALL_1kmMIPAS_varytemp';
hr_file   = 'RFM_FLUX/rad.dat1_HEAT_1kmMIPAS_varytemp';

[h,ha,p,pa] = rtpread(fop);
[h,p] = subset_rtp(h,p,[],[],iFov);

[d,w] = readkcstd(od_file);

disp('reading kcarta fluxes')
[flux,wf] = readkcflux(flux_file);
disp('reading kcarta hr')
[hr,wf] = readkcflux(hr_file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('doing RRTM')

addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM
addpath /home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB
fip = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/test.newklayers.ip.rtp';

rrtm_results = driver_rrtm_no_xsec_HP(fip,2,[]);

[hx,hax,px,pax] = rtpread(fip);
pxx = px;
pxx.gas_1 = interp1(log10(px.plevs),px.gas_1,log10(p.plevs));
pxx.gas_2 = interp1(log10(px.plevs),px.gas_2,log10(p.plevs));
pxx.gas_3 = interp1(log10(px.plevs),px.gas_3,log10(p.plevs));
pxx.gas_4 = interp1(log10(px.plevs),px.gas_4,log10(p.plevs));
pxx.gas_5 = interp1(log10(px.plevs),px.gas_5,log10(p.plevs));
pxx.gas_6 = interp1(log10(px.plevs),px.gas_6,log10(p.plevs));
pxx.ptemp = interp1(log10(px.plevs),px.ptemp,log10(p.plevs));
pxx.plevs = p.plevs;
pxx.nlevs = px.nlevs;

pxx.gas_1 = pxx.gas_1(1:85);
pxx.gas_2 = pxx.gas_2(1:85);
pxx.gas_3 = pxx.gas_3(1:85);
pxx.gas_4 = pxx.gas_4(1:85);
pxx.gas_5 = pxx.gas_5(1:85);
pxx.gas_6 = pxx.gas_6(1:85);
pxx.ptemp = pxx.ptemp(1:85);
pxx.plevs = pxx.plevs(1:85);
pxx.nlevs = 85;

fipx = '/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/AnuDudhia_RFM/xtest.newklayers.ip.rtp';
rtpwrite(fipx,hx,hax,pxx,pax);
rrtm_resultsx = driver_rrtm_no_xsec_HP(fipx,2,[],'/home/sergio/klayersV205/BinV201/klayers_airs_wetwater_1km');

figure(3);
semilogy(rrtm_results.kc605to2830_hr,rrtm_results.plevs,'bo-',...
         rrtm_resultsx.kc605to2830_hr,rrtm_resultsx.plevs,'r','linewidth',2)
hl = legend('AIRS 101 levels','1 kmlevels'); grid
set(gca,'ydir','reverse');
axis([-10 +2 0 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd /home/sergio/KCARTA/MATLAB/

%% now do up/down welling fluxes
% do_flux
% do_flux_630_2805
do_flux_arb

