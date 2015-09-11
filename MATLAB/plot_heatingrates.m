if length(strfind(kfile,'_HEAT')) > 0
  ix = input(' this looks like Heating Rates at all levels : plot the heating rates?? (-1/+1) : ');
  if ix > 0
    [mm,nn] = size(data);
    dw = mean(diff(wnums));

    fprintf(1,'   w(1) w(end) dw = %8.6f %8.6f %8.6f wn \n',wnums(1),wnums(end),dw)

    hr = sum(data)*dw;
    figure(1);
      plot(hr,hgt,'ro-','linewidth',2); grid
      title('heating rate')
      ylabel('hgt (km)')
  end

  ix = input(' subset the wavnumbers and then plot the fluxes and heating rates?? (-1/+1) : ');
  if ix > 0
    [mm,nn] = size(data);
    dw = mean(diff(wnums));

    fprintf(1,'w(1) w(end) dw = %8.6f %8.6f %8.6f wn \n',wnums(1),wnums(end),dw)
    raX = input('   Enter START/STOP wavenumbers : [a b] : ');
    if raX(1) < wnums(1)
      raX(1) = wnums(1);
    end
    if raX(2) > wnums(end)
      raX(2) = wnums(end);
    end
    iaX = find(wnums >= raX(1) & wnums < raX(2));
    
    hr = sum(data(iaX,:))*dw;
    figure(1);
      plot(hr,hgt,'ro-','linewidth',2); grid
      title('heating rate')
      ylabel('hgt (km)')
  end
end
