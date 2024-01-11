if iDoRad == 3  & iDoFlux < 0
  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    sedder = [sedder ' template_Qrad_allbands_coljac.nml  > ' outnml];
  elseif iCldORClr == +1 & iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qrad_2cloud_allbands_coljac.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outjacname '; echo $? >& status' num2str(iiBin)];
  eval(kcartaer)

  if iCldORClr == -1  | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    jall = [jall; j];
    dall = [dall; d];
    wall = [wall w];
  elseif iCldORClr == +1 & iBand >= 7
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    [mmjunk,nnjunk] = size(d);
    jall = [jall; j];
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];
  end

elseif iDoRad == 3 & iDoFlux == 5  
  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    sedder = [sedder ' template_Qflux5_allbands_coljac.nml  > ' outnml];
  elseif iCldORClr == +1 & iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qflux5_2cloud_allbands_coljac.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outjacname '; echo $? >& status' num2str(iiBin)];
  eval(kcartaer)

  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    jall = [jall; j];
    dall = [dall; d];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_OLR3']);
    fluxall = [fluxall; fx];
    rmer = ['!/bin/rm ' outname '_OLR3']; eval(rmer);

  elseif iCldORClr == +1 & iBand >= 7  
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    [mmjunk,nnjunk] = size(d);
    jall = [jall; j];
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_OLR3']);
    fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
    rmer = ['!/bin/rm ' outname '_OLR3']; eval(rmer);
  end

elseif iDoRad == 3 & iDoFlux == 7  
  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    sedder = [sedder ' template_Qflux7_allbands_coljac.nml  > ' outnml];
  elseif iCldORClr == +1 &  iBand >= 7
    sedder = [sedder ' ' sedderC ' '];
    sedder = [sedder ' template_Qflux7_2cloud_allbands_coljac.nml  > ' outnml];
  end
  eval(sedder)
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outjacname '; echo $? >& status' num2str(iiBin)];
  eval(kcartaer)

  if iCldORClr == -1 | iBand < 7
    %% can only do clear
    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    jall = [jall; j];
    dall = [dall; d];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_IOLR3']);
    fluxall = [fluxall; fx];
    rmer = ['!/bin/rm ' outname '_IOLR3']; eval(rmer);

  elseif iCldORClr == +1 & iBand >= 7  
    %cldparamsname = [outname '_CLD'];
    %cunk = load(cldparamsname);
    %[djunk,w] = readkcstd_twoslabclds(outname,cldparamsname);

    [d,w] = readkcstd(outname);
    [j,w] = readkcBasic(outjacnameCOL);
    [mmjunk,nnjunk] = size(d);
    jall = [jall; j];
    dall = [dall; d(:,nnjunk)];
    wall = [wall w];

    [fx,w] = readkcflux([outname '_IOLR3']);
    fluxall = [fluxall; squeeze(fx(nnjunk,:,:))];
    rmer = ['!/bin/rm ' outname '_IOLR3']; eval(rmer);
  end

else
  fprintf(1,'iDoRad = %3i iDoFlux = %3i \n',iDoRad,iDoFlux)
  error('huh?? unknown combo of iDoRad/iDoFlux!!!')
end

figure(1); plot(wall,rad2bt(wall,dall)); pause(0.1)
rmer = ['!/bin/rm ' outnml ' status' num2str(iiBin)]; eval(rmer)
