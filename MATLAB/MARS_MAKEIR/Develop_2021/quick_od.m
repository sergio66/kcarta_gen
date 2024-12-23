mars = load('mars.txt');

%% pretend this is a levels profile, need to check
plot(mars(:,2),1:90); xlabel('p (mb)'); ylabel('level number');

plot(mars(:,2),mars(:,1)); xlabel('p (mb)'); ylabel('height(km)');
%% goes from 6.4 mb to 3.2 mb in height from 0 to 7.5 km
%% p = p0 exp(-h/H)
%% log(p/p0) = -h/H ==> H = h / -log(p/p0) = 7.5/ -log(3.2/6.4) = 10.82 km

%% so roughly going from 6.4 to 7 mb would be a height change of 
%% h = H * -log(p/p0) = 10.82 * log(7/6.4) = 0.96 km

semilogy(mars(:,3),mars(:,2)); set(gca,'ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% pV = nRT ==> q = n/V L = p/RT L

L = diff(mars(:,1)) * 1000;

px = mars(:,2);
pN = px(1:end-1)-px(2:end);
pD = log(px(1:end-1)./px(2:end));
pav = pN ./ pD;

Tav = interp1(log(mars(:,2)),mars(:,3),log(pav));

semilogy(mars(:,3),mars(:,2),'b',Tav,pav,'r'); set(gca,'ydir','reverse')

R = 8.31;   %% gas constant

%% 1013 mb = 101325 N/m2
density = pav*101325/1013 .* L /R ./Tav;    %%% moles/m2
semilogy(density,pav); set(gca,'ydir','reverse')

q = density/10000 * 6.023e23;               %% molecules/cm2
semilogy(q,pav); set(gca,'ydir','reverse')

%  [index    press  part.press   temp   col.density ] =
%  [integer  real   real         real   real        ]  =
%  [nounits  atm    atm          K      kmoles/cm2  ]

MR = 0.95;   %% CO2 mixing ratio
data = [(1:length(q))' pav/1013.25 pav/1013.25*MR Tav  density/10000/1000*MR];
fid = fopen('mars_co2_amt','w');
fprintf(fid,'%3i %12.8e %12.8e %12.8f %12.8e \n',data');
fclose(fid);

addpath /home/sergio/SPECTRA
%[outwave,out_array] = run8(2,605,855,'./mars_co2_amt');
%error('oioi')
%[outwave,out_array] = driver_run8_parallel(2,605,855,'./mars_co2_amt',25);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OhOh, probably need finer res than 0.0025 cm-1

doppler_widths_wavenumber

error('OhOh');

%% run these separately on different processors

for wn = 605 : 25 : 830
  outname = ['mars_co2_ods_' num2str(wn) '.mat'];
  ee = exist(outname);
  if ee == 0
    [outwave,out_array] = run8(2,wn,wn+25,'./mars_co2_amt');
    saver = ['save ' outname '  outwave out_array data'];
    eval(saver)
  end
end

for wn = 830 : -25 : 605
  outname = ['mars_co2_ods_' num2str(wn) '.mat'];
  ee = exist(outname);
  if ee == 0
    [outwave,out_array] = run8(2,wn,wn+25,'./mars_co2_amt');
    saver = ['save ' outname '  outwave out_array data'];
    eval(saver)
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
