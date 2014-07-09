f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/150K_band2350_linemixR.dat';
d150 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/175K_band2350_linemixR.dat';
d175 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/200K_band2350_linemixR.dat';
d200 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/225K_band2350_linemixR.dat';
d225 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/250K_band2350_linemixR.dat';
d250 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/275K_band2350_linemixR.dat';
d275 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/300K_band2350_linemixR.dat';
d300 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/325K_band2350_linemixR.dat';
d325 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/350K_band2350_linemixR.dat';
d350 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/375K_band2350_linemixR.dat';
d375 = load(f);
f = '/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/400K_band2350_linemixR.dat';
d400 = load(f);

pa = [d150(:,3) d175(:,3) d200(:,3) d225(:,3) d250(:,3) d275(:,3) ...
      d300(:,3) d325(:,3) d350(:,3) d375(:,3) d400(:,3)];
plot(d400(:,2),pa);
grid

T = 150:25:400;
ii = 54; plot(T,pa(ii,:),'.-'); fprintf(1,'%8.6f \n',d400(ii,2));