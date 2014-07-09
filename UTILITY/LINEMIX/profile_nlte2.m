%this script produces the NLTE profile
heights = load('/home/sergio/MATLAB/airsheights.dat');

co2 = load('/home/sergio/SPECTRA/IPFILES/std_co2');
co2 = load('/home/sergio/SPECTRA/IPFILES/mls_co2');
[y,ii]=sort(ed_press); 
[h,ha,p,pa] = rtpread('/home/sergio/AIRSPRODUCTS/FIRSTLIGHT/edward.op.rtp');
theprof  = interp1(p.plevs(1:99),p.ptemp(1:99),co2(:,2)*1013);

edw=load('/home/sergio/KCARTA/SRCv1.11/NONLTE/midlat_day2.dat');
ed_press = edw(02:25,:); ed_press = ed_press(:);
ed_lte      = edw(27:50,:); ed_lte      = ed_lte(:);
ed_nlte2350 = edw(52:75,:);   ed_nlte2350 = ed_nlte2350(:);
ed_nlte2351 = edw(77:100,:);  ed_nlte2351 = ed_nlte2351(:);
ed_nlte2352 = edw(102:125,:); ed_nlte2352 = ed_nlte2352(:);
lte      = interp1(ed_press(ii),ed_lte(ii),co2(:,2)*1013);
nlte2350 = interp1(ed_press(ii),ed_nlte2350(ii),co2(:,2)*1013);
nlte2351 = interp1(ed_press(ii),ed_nlte2351(ii),co2(:,2)*1013);
nlte2352 = interp1(ed_press(ii),ed_nlte2352(ii),co2(:,2)*1013);

papa = [nlte2350 nlte2351 nlte2352];
semilogx(ed_press(ii),ed_lte(ii),...
         co2(:,2)*1013,co2(:,4),'.-',co2(:,2)*1013,papa)
semilogx(ed_press(ii),ed_lte(ii),ed_press(ii),ed_nlte2350(ii),...
 co2(:,2)*1013,co2(:,4),co2(:,2)*1013,nlte2350,'.-',co2(:,2)*1013,nlte2350,'+')

semilogx(co2(:,2)*1013,co2(:,4),'.-',...
         co2(:,2)*1013,(co2(:,4)-lte)./co2(:,4),'.-g',...
         co2(:,2)*1013,(co2(:,4)-papa(:,1))./co2(:,4),'.-r')

semilogx(co2(:,2)*1013,(lte-papa(:,1))./lte,'.-r')
plot(heights/1000,theprof-lte,'.-',heights/1000,theprof-papa(:,1),'.-');
grid
semilogx(co2(:,2)*1013,(co2(:,4)-lte)./co2(:,4),'.-g',...
         co2(:,2)*1013,(co2(:,4)-papa(:,1))./co2(:,4),'.-r')
grid

ii = find(isnan(nlte)); nlte(ii) = co2(ii,4); ii

newstuff = co2; newstuff(:,4) = nlte;
fid = fopen('/taro/s1/sergio/AIRSCO2/NONLTE/hit2350_day_profile3','w');
fprintf(fid,'%3i  %8.6e %8.6e  %8.6f  %8.6e \n', newstuff');
fclose(fid);