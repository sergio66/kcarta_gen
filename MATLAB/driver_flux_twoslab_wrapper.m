function [flux,fluxdata] = driver_flux_twoslab_wrapper(hIN,pIN,iRTP);

%% this generic code takes in an arbitrary rtp file (and rtp profile within it),
%% and then breaks it apart so that it can do a flux calc by summing over the clear,
%% cloudy and overlap cloudy streams

%% really only at 100 cm-1 do I start to kick in different scattering tables!!!!

chunkstart = [15 30 50  80 140 300 500  605 2830]; 
chunkstop  = [30 50 80 140 300 500 605 2830 3580]; 

kc0 = '/home/sergio/KCARTA/BIN/kcarta_clear.x';
kc1 = '/home/sergio/KCARTA/BIN/kcarta.x';

%% first subset the rtp profile we are interested in
[h0,p0] = subset_rtp_clouds(hIN,pIN,[],[],iRTP);

if length(intersect(hIN.glist,[201 202 203])) > 0
  error('this seems to be the profile for 100 layer clouds ... expecting TwoSlab')
end

rand('seed', sum(1000*clock))

outnmlx = ['JUNK/outnml_' num2str(iRTP) '.nml'];
outnmlx = mktemp(outnmlx);

if hIN.ptype == 0
  error('this seems to be LEVELS profile');
end

%% then do checks
if p0.cfrac <= eps & p0.cfrac2 <= eps
  p0.cfrac = 0;
  p0.cngwat = 0;
  p0.cfrac2 = 0;
  p0.cngwat2 = 0;
  p0.cfrac12 = 0;
  iType0 = 0;   %% no clouds
elseif p0.cfrac > eps & p0.cfrac2 <= eps
  p0.cfrac2 = 0;
  p0.cngwat2 = 0;
  p0.cfrac12 = 0;
  iType0 = 1;   %% cloud 1
elseif p0.cfrac <= eps & p0.cfrac2 > eps
  p0.cfrac = 0;
  p0.cngwat = 0;
  p0.cfrac12 = 0;
  iType0 = 2;   %% cloud 2
else
   iType0 = 12; %% clouds 1 and 2
end

c11 = p0.cfrac  - p0.cfrac12;
c22 = p0.cfrac2 - p0.cfrac12;
fclr = 1 - [c11 + c22] - p0.cfrac12;

flux = [];

ff        = [];
clearflux = [];
cloud1flux = [];
cloud2flux = [];
cloud12flux = [];

['ctype  ' num2str(p0.ctype)  ' ' num2str(p0.ctype2)]
['cngwat ' num2str(p0.cngwat) ' ' num2str(p0.cngwat2)]
['cpsize ' num2str(p0.cpsize) ' ' num2str(p0.cpsize2)]
['cprtop ' num2str(p0.cprtop) ' ' num2str(p0.cprtop2)]
['cprbot ' num2str(p0.cprbot) ' ' num2str(p0.cprbot2)]
['cfrac  ' num2str(p0.cfrac)  ' ' num2str(p0.cfrac2) ' '  num2str(p0.cfrac12)]

tmp_rtp = mktemp('thetmp_rtp.rtp');
tmp_rtpX = tmp_rtp(6:end);   %% skip the '/tmp/'

tmp_rtp = ['JUNK/temprtp_' num2str(iRTP) 'try.rtp'];
tmp_rtpX = tmp_rtp(6:end);   %%% skip the JUNK/;

fprintf(1,'tmp_rtp  = %s \n',tmp_rtp)
fprintf(1,'tmp_rtpX = %s \n',tmp_rtpX)

tic
for cc = 1 : length(chunkstart)-1
  f1 = chunkstart(cc);
  f2 = chunkstop(cc);

  if f2 <= 140
    iType = 0;   %% do not have scattering tables for very small wavenumber, so just do clear sky
  else
    iType = iType0;  %% have scattering tables, so go ahead and do clear/cloudy calcs
  end

  %% this defines the cloud scattering tables
  if f1 < 605
    iEFGH = 105;
    iABCD = 1905;
  else
    iEFGH = 405;
    iABCD = 3005;
  end
    
  sedder0 = ['!sed -e "s/FF1/' num2str(f1) '/g"  -e "s/FF2/' num2str(f2) '/g" '];
  sedder0 = [sedder0 ' -e "s/XYZ/'    num2str(iRTP) '/g"'];
  sedder0 = [sedder0 ' -e "s/ABCD/'    num2str(iABCD) '/g"'];
  sedder0 = [sedder0 ' -e "s/EFGH/'    num2str(iEFGH) '/g"'];
  sedder0 = [sedder0 ' -e "s/FNAME/'   tmp_rtpX '/g"'];

  hx = h0;
  px = p0;

  hx.vcmin = f1;
  hx.vcmax = f2;

  ff = [ff (f1+f2)/2];

  if iType == 0
    %% NO clouds, just do clear sky fluxes
    fprintf(1,'  ---> iType = 0, << doing clearsky flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    rtpwrite(tmp_rtp,hx,[],px,[]);
    kcX = kc0; flux_twoslab_wrapper
    clearflux   = [clearflux nansum(flux2s)*mean(diff(w))];
    if iType0 == 0
      cloud1flux  = [cloud1flux  0];
      cloud2flux  = [cloud2flux  0];
      cloud12flux = [cloud12flux 0];
    else
      cloud1flux   = [cloud1flux nansum(flux2s)*mean(diff(w))];
      cloud2flux   = [cloud2flux nansum(flux2s)*mean(diff(w))];
      cloud12flux  = [cloud12flux nansum(flux2s)*mean(diff(w))];
    end

  elseif iType == 1
    %%  cloud 1, just do clear sky fluxes and cloud1 fluxes
    fprintf(1,'  ---> iType = 1, << doing clearsky flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac   = 0;
    pxX.cfrac12 = 0;
    pxX.cngwat  = 0;    
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc0; flux_twoslab_wrapper
    clearflux = [clearflux nansum(flux2s)*mean(diff(w))];

    fprintf(1,'  ---> iType = 1, << doing cloud1 flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc1; flux_twoslab_wrapper    
    cloud1flux = [cloud1flux nansum(flux2s)*mean(diff(w))];

    cloud2flux  = [cloud2flux  0];
    cloud12flux = [cloud12flux 0];

  elseif iType == 2
    %%  cloud 2, just do clear sky fluxes and cloud1 fluxes
    fprintf(1,'  ---> iType = 2, << doing clearsky flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac2  = 0;
    pxX.cfrac12 = 0;
    pxX.cngwat2 = 0;    
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc0; flux_twoslab_wrapper
    clearflux = [clearflux nansum(flux2s)*mean(diff(w))];

    fprintf(1,'  ---> iType = 2, << doing cloud1 flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac12 = 0;
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc1; flux_twoslab_wrapper    
    cloud2flux = [cloud2flux nansum(flux2s)*mean(diff(w))];

    cloud1flux  = [cloud1flux  0];
    cloud12flux = [cloud12flux 0];

  elseif iType == 12
    %%  cloud 1 and 2, do clear sky fluxes and cloud1,cloud2, cloud12 fluxes
    fprintf(1,'  ---> iType = 12, << doing clearsky flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac   = 0;
    pxX.cngwat  = 0;    
    pxX.cfrac2  = 0;
    pxX.cngwat2 = 0;    
    pxX.cfrac12 = 0;
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc0; flux_twoslab_wrapper
    clearflux = [clearflux nansum(flux2s)*mean(diff(w))];

    fprintf(1,'  ---> iType = 12, << doing cloud1 flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac2   = 0;
    pxX.cngwat2  = 0;
    pxX.cfrac12  = 0;    
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc1; flux_twoslab_wrapper    
    cloud1flux = [cloud1flux nansum(flux2s)*mean(diff(w))];

    fprintf(1,'  ---> iType = 12, << doing cloud2 flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    pxX.cfrac   = 0;
    pxX.cngwat  = 0;    
    pxX.cfrac12 = 0;
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc1; flux_twoslab_wrapper    
    cloud2flux = [cloud2flux nansum(flux2s)*mean(diff(w))];

    fprintf(1,'  ---> iType = 12, << doing cloud12 flux at [f1,f2] = %4i %4i >> \n',f1,f2)
    pxX = px;
    rtpwrite(tmp_rtp,hx,[],pxX,[]);
    kcX = kc1; flux_twoslab_wrapper    
    cloud12flux = [cloud12flux nansum(flux2s)*mean(diff(w))];

  end

  flux = clearflux/1000*fclr + cloud1flux/1000*c11 + cloud2flux/1000*c22 + ...
         cloud12flux/1000*p0.cfrac12;

  plot(ff,clearflux/1000,ff,cloud1flux/1000,ff,cloud2flux/1000,ff,cloud12flux/1000); hold on
  plot(ff,flux,'k','linewidth',2); hold off
  pause(1)

  disp(' ')

end
toc

rmer = ['!/bin/rm ' outnmlx];
eval(rmer);

flux = nansum(clearflux)   * fclr + ...
       nansum(cloud1flux)  * c11 + ...
       nansum(cloud2flux)  * c22 + ...
       nansum(cloud12flux) * p0.cfrac12;
flux = flux/1000;   %% change from mW/m2 to W/m2

fluxdata.ff          = ff;
fluxdata.clearflux   = clearflux/1000;
fluxdata.cloud1flux  = cloud1flux/1000;
fluxdata.cloud2flux  = cloud2flux/1000;
fluxdata.cloud12flux = cloud12flux/1000;

fluxdata.c1  = p0.cfrac;
fluxdata.c11 = c11;
fluxdata.c2  = p0.cfrac2;
fluxdata.c22 = c22;
fluxdata.c12  = p0.cfrac12;
fluxdata.fclr = fclr;
