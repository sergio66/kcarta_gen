iG = 2;

band = [2310 2320 2350 2351 2352];
iIso = [1 1 1 2 3];
jlow =  [ 1  1  1  1  1];
jhigh = [24 16  9  9  9]; 
for ii = 150:10:400
  for jj = 1 : 5

    fprintf(1,' %4i  %5i \n',ii,band(jj));

    %read(iIOUN,*) iGasIDy,iNum,iIso,jLow,jHigh%


    d1 = ['/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/' num2str(ii) 'K_band'];
    d1 = [d1 num2str(band(jj)) '_linemixP.dat'];

    d2 = ['/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/' num2str(ii) 'K_band'];
    d2 = [d2 num2str(band(jj)) '_linemixR.dat'];

    d3 = ['/home/sergio/KCARTA/SRCv1.11/NONLTE/LINEMIX/' num2str(ii) 'K_band'];
    d3 = [d3 num2str(band(jj)) '_linemix.dat'];

    d11 = load(d1);
    d22 = load(d2);

    d33 = [d11; d22]';
    iNum = length(d33);
    fid = fopen(d3,'w'); 
    fprintf(fid,'%2i %5i %2i %2i %2i\n',iG,iNum,iIso(jj),jlow(jj),jhigh(jj));

    fprintf(fid,'%4i %8.6e   %8.6e \n',d33); 
    fclose(fid); 

    end
  end


  