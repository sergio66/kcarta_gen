function compare_fluxNrad_calcss(w,rad,flux,UD)

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

%{
> flux_levels = compute_flux_from_4radstreams(w,d2x);
number of levels =  97
Out of memory. Type HELP MEMORY for your options.

Error in compute_flux_from_4radstreams (line 51)
  flux_levels = flux_levels + daW(jj)*rads(:,ind)*2*pi;

>> whos
  Name                  Size                     Bytes  Class     Attributes

  ans                   1x1                          8  double
    d20              890000x97                 690640000  double
      d2x              890000x388               2762560000  double
        flux1            890000x196               1395520000  double
	  flux_levels      890000x97                 690640000  double
	    w                     1x890000               7120000  double

>> clear flux_levels
>> flux_levels = compute_flux_from_4radstreams(w,d2x);
number of levels =  97
>> sergio = sum(flux_levels,1) * dw/1000;
Undefined function or variable 'dw'.

>> dw=0.0025;
>> sergio = sum(flux_levels,1) * dw/1000;
>> plot(sergio-bah(2:86),flipud(plevsIN(1:85)),'rx-'); grid
Undefined function or variable 'bah'.

>> bah = flux1(:,98:196); bah = sum(bah,1)*dw/1000;
>> plot(sergio-bah(2:97))
Error using -
Matrix dimensions must agree.

>> whos sergio bah
  Name        Size            Bytes  Class     Attributes

  bah         1x99              792  double
    sergio      1x97              776  double

>> bah = flux1(:,99:196); bah = sum(bah,1)*dw/1000;
>> plot(sergio-bah(2:97))
Error using -
Matrix dimensions must agree.

>> whos sergio bah
  Name        Size            Bytes  Class     Attributes

  bah         1x98              784  double
    sergio      1x97              776  double

>> plot(sergio-bah(2:98))
>> plot(sergio)
>> plot(sergio,1:97)
>> plot(sergio,1:97,bah(2:98),1:97)
>> plot(sergio,1:97,flipud(bah(2:98)),1:97)
>>
>> plot(sergio,1:97,fliplr(bah(2:98)),1:97)
>> plot(sergio,1:97,fliplr(bah(1:97)),1:97)
>> plot(sergio-fliplr(bah(1:97)),1:97)

%}
