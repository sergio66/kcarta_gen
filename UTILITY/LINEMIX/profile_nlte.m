%this script produces the NLTE profile
edw=load('/home/sergio/KCARTA/SRCv1.11/NONLTE/midlat_day.dat');

ed_press = edw(02:25,:); ed_press = ed_press(:);
ed_lte   = edw(27:50,:); ed_lte   = ed_lte(:);
ed_nlte  = edw(52:75,:); ed_nlte  = ed_nlte(:);

co2 = load('/home/sergio/SPECTRA/IPFILES/std_co2');
[y,ii]=sort(ed_press); 

nlte = interp1(ed_press(ii),ed_nlte(ii),co2(:,2)*1013);

semilogx(ed_press(ii),ed_lte(ii),ed_press(ii),ed_nlte(ii),...
       co2(:,2)*1013,co2(:,4),co2(:,2)*1013,nlte,'.-',co2(:,2)*1013,nlte,'+')

ii = find(isnan(nlte)); nlte(ii) = co2(ii,4); ii

newstuff = co2; newstuff(:,4) = nlte;
fid = fopen('/taro/s1/sergio/AIRSCO2/NONLTE/hit2350_day_profile3','w');
fprintf(fid,'%3i  %8.6e %8.6e  %8.6f  %8.6e \n', newstuff');
fclose(fid);