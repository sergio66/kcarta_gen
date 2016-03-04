if length(strfind(kfile,'_ALL')) > 0
  ix = input(' this looks like up/dn fluxes at all levels : plot the fluxes and heating rates?? (-1/+1) : ');
  if ix > 0
    [mm,nn] = size(data);
    dw = mean(diff(wnums));
    nnx = 1:nn/2;

    fprintf(1,'   w(1) w(end) dw = %8.6f %8.6f %8.6f wn \n',wnums(1),wnums(end),dw)

    upflux = sum(data(:,1:nn/2))*dw/1000;
    dnflux = sum(data(:,nn/2+(1:nn/2)))*dw/1000;
      fprintf(1,'     up/dn flux at GND = %8.6f   %8.6f W/m2 \n',upflux(1),dnflux(1))
      fprintf(1,'     up/dn flux at TOA = %8.6f   %8.6f W/m2 \n',upflux(end),dnflux(end))      
    netxkc = upflux - dnflux;  %% upwell - downwell flux at each level

    n = -1;
    if n < 0
      plevsx = plevs;
    else
      plevsx = round(plevs*10^n)/10^n;
    end

    figure(1);
      plot(upflux,hgt,'ro-',dnflux,hgt,'bo-'); hold on; 
      plot(netxkc,hgt,'k',diff(netxkc) ./  diff(plevsx') * 8.4391,hgt(1:end-1),'g','linewidth',2); hold off; grid
      title(' (r) upwell flux (b) dnwell flux W/m2 \newline (k) net=up-dn W/m2 (g) 10*flux div=dnet/dz  W/m2/layer',...
            'fontsize',10);
      ylabel('hgt (km)')

    figure(2);
      plot(upflux,hgt,'ro-',-dnflux,hgt,'bo-'); grid
      title(' (r) upwell flux (b) dnwell flux W/m2')
      ylabel('hgt (km)')

    figure(3); 
      %n = input('   enter number of decimal points for HR quickcalc (1,2,3 ... or -1 for all) : ');
      n = -1;
      if n < 0
        plevsx = plevs;
      else
        plevsx = round(plevs*10^n)/10^n;
      end
      htxkc  = diff(netxkc) ./ diff(plevsx') * 8.4391;
      plot(htxkc,hgt(1:end-1),'r'); grid on; title('cooling rate K/day')
      ylabel('hgt (km)')

      hr(:,1) = htxkc;
      hr(:,2) = hgt(1:end-1);
  end

  ix = input(' subset the wavnumbers and then plot the fluxes and heating rates?? (-1/+1) : ');
  if ix > 0
    [mm,nn] = size(data);
    dw = mean(diff(wnums));
    nnx = 1:nn/2;

    fprintf(1,'w(1) w(end) dw = %8.6f %8.6f %8.6f wn \n',wnums(1),wnums(end),dw)
    raX = input('   Enter START/STOP wavenumbers : [a b] : ');
    if raX(1) < wnums(1)
      raX(1) = wnums(1);
    end
    if raX(2) > wnums(end)
      raX(2) = wnums(end);
    end
    iaX = find(wnums >= raX(1) & wnums < raX(2));
    
    upflux = sum(data(iaX,1:nn/2))*dw/1000;
    dnflux = sum(data(iaX,nn/2+(1:nn/2)))*dw/1000;
      fprintf(1,'     up/dn flux at GND = %8.6f   %8.6f W/m2 \n',upflux(1),dnflux(1))
      fprintf(1,'     up/dn flux at TOA = %8.6f   %8.6f W/m2 \n',upflux(end),dnflux(end))          
    netxkc = upflux - dnflux;  %% upwell - downwell flux at each level

    figure(1);
      plot(upflux,hgt,'ro-',dnflux,hgt,'bo-'); hold on; 
      plot(netxkc,hgt,'k',diff(netxkc)*10,hgt(1:end-1),'g','linewidth',2); hold off; grid
      title(' (b) upwell flux (r) dnwell flux W/m2 \newline (k) net=up-dn W/m2 (g) 10*flux div=dnet/dz  W/m2/layer',...
            'fontsize',10);
      ylabel('hgt (km)')

    figure(2);
      plot(upflux,hgt,'ro-',-dnflux,hgt,'bo-'); grid
      title(' (r) upwell flux (b) dnwell flux W/m2')
      ylabel('hgt (km)')

    figure(3); 
      n = input('   enter number of decimal points for HR quickcalc (1,2,3 ... or -1 for all) : ');
      if n < 0
        plevsx = plevs;
      else
        plevsx = round(plevs*10^n)/10^n;
      end
      htxkc  = diff(netxkc) ./ diff(plevsx') * 8.4391;
      plot(htxkc,hgt(1:end-1),'r'); grid on; title('cooling rate K/day')
      ylabel('hgt (km)')

      hr(:,1) = htxkc;
      hr(:,2) = hgt(1:end-1);
  end

end
