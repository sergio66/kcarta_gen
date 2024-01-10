if iDoRad == 3  & iDoFlux < 0
  if (iCldORClr == -1 | iBand < 7)
    %% can only do clear
    sedder = [sedder ' template_Qrad_allbands.nml  > ' outnml];
  elseif iCldORClr == +1 & iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qrad_2cloud_allbands.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
  if exist(outname)
    rmer = ['!/bin/rm ' outname];
    eval(rmer);
  end    
  eval(kcartaer)

  if iCldORClr == -1  | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    dall = [dall; d];
    wall = [wall w];
  elseif iCldORClr == +1 & iBand >= 7
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [mmjunk,nnjunk] = size(d);
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];
  end

elseif iDoRad == 3 & iDoFlux == 5  
  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    sedder = [sedder ' template_Qflux5_allbands.nml  > ' outnml];
  elseif iCldORClr == +1 & iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qflux5_2cloud_allbands.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
  if exist(outname)
    rmer = ['!/bin/rm ' outname];
    eval(rmer);
  end    
  eval(kcartaer)

  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    dall = [dall; d];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_OLR3']);
    fluxall = [fluxall; fx];

  elseif iCldORClr == +1 & iBand >= 7  
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [mmjunk,nnjunk] = size(d);
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_OLR3']);
    fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
  end

elseif iDoRad == 3 & iDoFlux == 7  
  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    sedder = [sedder ' template_Qflux7_allbands.nml  > ' outnml];
  elseif iCldOrClr == +1 &  iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qflux7_2cloud_allbands.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& status' num2str(iiBin)];
  if exist(outname)
    rmer = ['!/bin/rm ' outname];
    eval(rmer);
  end    
  eval(kcartaer)

  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    dall = [dall; d];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_IOLR3']);
    fluxall = [fluxall; fx];

  elseif iCldORClr == +1 & iBand >= 7  
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [mmjunk,nnjunk] = size(d);
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_IOLR3']);
    fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
  end

else
  fprintf(1,'iDoRad = %3i iDoFlux = %3i \n',iDoRad,iDoFlux)
  error('huh?? unknown combo of iDoRad/iDoFlux!!!')
end

figure(2); plot(wall,rad2bt(wall,dall)); pause(0.1)
rmer = ['!/bin/rm ' outnml ' status' num2str(iiBin)]; eval(rmer)
