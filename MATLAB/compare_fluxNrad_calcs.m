function compare_fluxNrad_calcs(w,rad,flux,UD)

%% compares RT calcs at all levels (except at start boundary) vs flux at all levels
%% M = number of layers, so M+1 = number of levels
%% so flux at all levels = (M+1) levels, and we provide up/down, so 2(M+1)
%% to do flux from rads, we need the rads, at four angles, so 4*M
%%
%%   w    = wavenumber  (1 x Npts)
%%   rads = radiances (N x (M*4))
%%   flux = fluxes    (N * 2*(M+1)), first half are upwelling, last half are downwelling
%%   UD   = direction (1 for upwell,-1 for dnwell)

[m1,n1] = size(w);
[nRad,mRad] = size(rad);
[nFlux,mFlux] = size(flux);

if n1 ~= nRad
  error('different number of wavenumber points in supplied rad')
elseif n1 ~= nFlux
  error('different number of wavenumber points in supplied flux')
end

if (mRad/4 ~= mFlux/2-1)
  error('different number of layers in supplied rad/flux')
end

dw = mean(diff(w));

figure(1); clf
flux_levels = compute_flux_from_4radstreams(w,rad);
quick_flux_from_rads = sum(flux_levels,1) * dw/1000;   %% flux as computed from 4 rads

M = mRad/4;  %% number of layers
if UD == -1
  %% looking at downwelling flux
  levs = (1:M+1);
  levs = levs + (M+1);
  subset_flux = flux(:,levs);
  subset_flux = sum(subset_flux,1)*dw/1000;
  subset_flux = fliplr(subset_flux);    %% need to flip dir
  subset_flux = subset_flux(2:end);     %% we only need from levels 2 : end
elseif UD == +1
  %% looking at upwelling flux
  levs = (1:M+1);
  subset_flux = flux(:,levs);
  subset_flux = sum(subset_flux,1)*dw/1000;
  subset_flux = subset_flux(2:end);     %% we only need from levels 2 : end
end
figure(1); clf; plot(quick_flux_from_rads,1:M,'bx-',subset_flux,1:M,'r','linewidth',2); xlabel('flux W/m2'); ylabel('level');
  hl = legend('from 4 kcarta rads','from kcarta flux','location','east'); set(hl,'fontsize',10);; grid
  if UD == -1
    set(gca,'ydir','reverse')
  end
figure(2); clf; plot(quick_flux_from_rads - subset_flux,1:M,'linewidth',2); grid
  xlabel('\delta flux W/m2'); ylabel('level');
  if UD == -1
    set(gca,'ydir','reverse')
  end
