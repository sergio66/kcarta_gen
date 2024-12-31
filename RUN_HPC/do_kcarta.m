if ~exist('loop_allprofiles_onefile')
  %% generic usual stuff
  outname    = ['JUNK/rad.dat' num2str(iiBin)];
  outnamejac = ['JUNK/jac.dat'  num2str(iiBin)];
  outnml     = ['run_nml'  num2str(iiBin)];
  if exist('iDISopt')
    outnml     = ['run_nml'  num2str(iiBin) '_iDISopt_' num2str(iDISopt,'%02i')];
  end
  outstat    = ['status'  num2str(iiBin)];
else
  outname    = fxrad;
  outnamejac = fxjac;
  outnml     = fxnml;
  if exist('iDISopt')
    outnml     = [outnml '_iDISopt_' num2str(iDISopt,'%02i')];
  end
  outstat    = fxstat;
end

if ~exist('%iIRorFIR')
  iIRorFIR  = -1;
  iIRorFIR = +1;
end
if ~exist('iDISopt')
   iDISopt = -1;
end
if iDISopt ~= 11 & iDISopt ~= 12
  if iIRorFIR == +1
    %% 89 chunks
    f1 = 605;
    f2 = 2830;
    iKCKD = iKCKD;
  elseif iIRorFIR == -1
    %% 20 chunks
    f1 = 310;
    f2 = 510;
    iKCKD = 1;
  end
  fprintf(1,'f1,f2 RESET TO %4i %4i \n',f1,f2);
else
  fprintf(1,'DISORT chunking : f1,f2 REMAINS at %4i %4i \n',f1,f2);
end

sedder = ['!sed -e "s/FF1/' num2str(f1) '/g"  -e "s/FF2/' num2str(f2) '/g" '];
if ~exist('loop_allprofiles_onefile')
  %% generic usual stuff
  sedder = [sedder ' -e "s/MMM/'     num2str(iiBin) '/g"'];
else
  sedder = [sedder ' -e "s/MMM/'     num2str(ttBin) '/g"'];
end
sedder = [sedder ' -e "s/TTT/'     num2str(iDo_rt_1vs43,'%02d') '/g"'];
sedder = [sedder ' -e "s/XYZXYZ/'  use_this_rtp   '/g"'];  %%%<<<-- to sed rtpfname
sedder = [sedder ' -e "s/CKDCKD/'  num2str(iKCKD) '/g"'];  

fprintf(1,'outname = %s \n outnml = %s \n',outname,outnml)

if length(uncstr) == 2
  fprintf(1,'uncstr = %s \n',uncstr)
  sedder = [sedder ' -e "s/UNCUNC/'  uncstr '/g"'];
end

if iDoLBLRTM > 0
  %% do three gases : CO2, O3 and CH4 using LBLRTM ods  
  fprintf(1,' >>>>>>>>> iDoLBLRTM = %2i ... using other ods \n',iDoLBLRTM);
  %sedder = [sedder ' -e "s/DOLBLRTM/3/g"'];
  sedder = [sedder ' -e "s/DOLBLRTM/' num2str(iDoLBLRTM) '/g"'];  
else
  sedder = [sedder ' -e "s/DOLBLRTM/-1/g"'];
end
sedder = [sedder ' -e "s/COMMENTCOMMENT/' caComment '/g"'];

if iDoRad == 0  %% special  eg in nm_radnce, we have iAtmLoop = 3 etc
  disp('do_kcarta.m here 0')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];
  sedder = [sedder ' template_Qspecial.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];

elseif iDoRad == 10  %% rads, using NNalli Emis
  disp('do_kcarta.m here 10 Nalli')
  addpath /asl/matlib/h4tools
  [h,ha,p,pa] = rtpread(use_this_rtp0);
  if abs(p.satzen(iiBin)) > 90
    error('do_kcarta.m : p.satzen > 90')
  end
  sedder = [sedder ' -e "s/SATZENSATZEN/'    num2str(p.satzen(iiBin)) '/g"']; 
  sedder = [sedder ' template_Qrad0_NalliEmiss.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoRad == 20  %% rads, using DISORT
  disp('do_kcarta.m here 20 DISORT')
  outstat  = [outstat  '_' num2str(f1,'%04d')];
  outnml   = [outnml  '_' num2str(f1,'%04d')];
  outname  = [outname '_' num2str(f1,'%04d') '_iDISopt_' num2str(iDISopt,'%02i')];
  type1nml = 'template_Qradcloud_2cloud_DISORT.nml';
  type2nml = 'template_Qradcloud_2cloud_DISORT.nml';
  sed_1_2_cloudfiles
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoRad == 1  %% individual gas OD
  disp('do_kcarta.m here 1')
  sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];
  sedder = [sedder ' template_QgasOD.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoRad == 2  %% cumulative mixed path ODs
  disp('do_kcarta.m here 2')
  sedder = [sedder ' template_QsumMP.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud < 0 & gg ~= 1001 & gg ~= 2346 & gg ~= 2456 & gg ~= 5912
  disp('do_kcarta.m here A1 (clr sky rads and jacs, not WV1,101,102,103)')
  sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Qjacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud < 0 & gg == 1001
  disp('do_kcarta.m here A2 (clr sky rads and jacs for WV1,101,102,103)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_QWVjacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud < 0 & gg == 2346
  disp('do_kcarta.m here A3 (clr sky rads and jacs for GID 2,3,4,6,51,52)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Q2346_51_52jacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud < 0 & gg == 2456
  disp('do_kcarta.m here A3X (clr sky rads and jacs for GID 2,4,5,6,51,52)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Q2456_51_52jacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud < 0 & gg == 5912
  disp('do_kcarta.m here A3 (clr sky rads and jacs for GID 5,9,11,12,61,103)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Q5_9_11_12_61_103_jacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud < 0 & gg ~= 1001 & gg ~= 2346 & gg ~= 5912 & gg ~= 2456
  disp('do_kcarta.m here A4 (clr sky rads and col jacs, not WV1,101,102,103)')
  sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Qcoljacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud < 0 & gg == 1001
  disp('do_kcarta.m here A4 (clr sky rads and col jacs, WV1,101,102,103)')
  %%% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_QcolWV_jacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
elseif iDoJac == -100 & iDoFlux < 0 & iDoCloud < 0 & gg == 1001
  disp('do_kcarta.m here A4 (UPLOOK clr sky rads and col jacs, WV1,101,102,103)')
  %%% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_QcolWV_uplook_jacobian.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
elseif iDoJac == -100 & iDoFlux < 0 & iDoCloud < 0 & gg == 1003
  disp('do_kcarta.m here A4 (UPLOOK clr sky rads and col jacs, WV1,101,102,103 and G3)')
  %%% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_QcolWV_uplook_jacobian3.nml  > ' outnml];      
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud < 0 & gg == 2346
  disp('do_kcarta.m here A5 (clr sky rads and col jacs, for G2,3,4,6,51,52)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Qcol2346_51_52jacobian.nml  > ' outnml];  %currently g 2,4,5,6,51,52    
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud < 0 & gg == 2456
  disp('do_kcarta.m here A5 (clr sky rads and col jacs, for G2,4,5,6,51,52)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Qcol2456_51_52jacobian.nml  > ' outnml];  %currently g 2,4,5,6,51,52    
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud < 0 & gg == 5912
  disp('do_kcarta.m here A5 (clr sky rads and col jacs, for G5,9,11,12,61,103)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  sedder = [sedder ' template_Qcol5_9_11_12_61_103_jacobian.nml template_Qcol2346_51_52jacobian.nml  > ' outnml];  %currently g 5,9,11,12,61,103
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
%%%%%%%%%%%%%%%%%%%%%%%%%
elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud == 1 & gg == 1001
  disp('do_kcarta.m here AA1 (allsky rads and jacs for WV1,101,102,103)')
  type1nml = 'template_QWVcldjacobian_1cloud.nml';
  type2nml = 'template_QWVcldjacobian_2cloud.nml';
  sed_1_2_cloudfiles
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 1 & iDoFlux < 0 & iDoCloud == 1 & gg == 2346
  disp('do_kcarta.m here AA2 (allsky rads and jacs for G!I 2,3,4,6)')
  type1nml = 'template_Q2346cldjacobian_1cloud.nml';
  type2nml = 'template_Q2346cldjacobian_2cloud.nml';
  sed_1_2_cloudfiles
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];

elseif iDoJac == 100 & iDoFlux < 0 & iDoCloud == 1 & gg == 2346
  disp('do_kcarta.m here A5 (allsky rads and col jacs, for G2,3,4,6,51,52)')
  %% sedder = [sedder ' -e "s/GGG/'    num2str(gg) '/g"'];   %% this gives gasID for jacobian
  type1nml = 'template_Qcol2346cldjacobian_1cloud.nml';
  type2nml = 'template_Qcol2346cldjacobian_2cloud.nml';
  type1nml = 'template_Qcol2346_51_52_cldjacobian_1cloud.nml';
  type2nml = 'template_Qcol2346_51_52_cldjacobian_2cloud.nml';
  sed_1_2_cloudfiles
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname ' ' outnamejac '; echo $? >& ' outstat];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif iDoJac < 0 & iDoFlux < 0 & iDoCloud < 0 & (iHITRAN ~= 2016.3 & iHITRAN ~= 2017 & iHITRAN ~= 2018)
  disp('do_kcarta.m here B (clr sky rad, usual res)')
  sedder = [sedder ' template_Qrad.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux < 0 & iDoCloud < 0 & (iHITRAN == 2017 | iHITRAN == 2018)
  disp('do_kcarta.m here B (clr sky rad, HIGH or VHIGH res, 3 GASES ONLY)')
  sedder = [sedder ' template_Qrad_high_veryhigh_res.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux < 0 & iDoCloud < 0 & (iHITRAN == 2016.3)
  disp('do_kcarta.m here B (clr sky rad, USUAL RES, 3 GASES ONLY)')
  sedder = [sedder ' template_Qrad_high_veryhigh_res_yesaltdatabase.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];

%%%%%%%%%%%%%%%%%%%%%%%%%  
elseif iDoJac < 0 & iDoFlux == 2 & iDoCloud < 0
  disp('do_kcarta.m here HEATING RATES C2')
  sedder = [sedder ' template_Qflux2.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux == 5 & iDoCloud < 0
  disp('do_kcarta.m here UP/DOWN ILR/OLR FLUX C5')
  sedder = [sedder ' template_Qflux5.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux == 6 & iDoCloud < 0
  disp('do_kcarta.m here UP/DOWN FLUX C6')
  sedder = [sedder ' template_Qflux6.nml  > ' outnml];
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux < 0 & iDoCloud == 1
  disp('do_kcarta.m here DoCld')
  type1nml = 'template_Qradcloud_1cloud.nml';
  type2nml = 'template_Qradcloud_2cloud.nml';
  sed_1_2_cloudfiles  
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];
  
elseif iDoJac < 0 & iDoFlux < 0 & iDoCloud == 100
  disp('do_kcarta.m here DoCld100')
  type1nml = 'template_Qradcloud100_1cloud.nml';
  type2nml = 'template_Qradcloud100_2cloud.nml';
  sed_1_2_cloudfiles  
  write_cloudfrac_profileMMM_100
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname];  
  kcartaer = ['!time ' kcartaexec ' ' outnml ' ' outname '; echo $? >& ' outstat];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
else
  disp(' ')
  fprintf(1,'iDoJac  iDoFlux iDoCLoud = %4i %4i %4i \n',iDoJac,iDoFlux,iDoCloud)
  error('hmmmmm did not find specs to go grab a template nml file')
end

fprintf(1,'sedder str   = %s \n',sedder)
fprintf(1,'kcartaer str = %s \n',kcartaer)

eval(sedder)
eval(kcartaer)

%junk = ['status' num2str(iiBin)];
junk = outstat;
loader = ['exitcode = load(''' junk ''');'];
eval(loader);
fprintf(1,'iiBin,exitcode = %6i %2i \n',iiBin,exitcode)

rmer = ['!/bin/rm ' outnml ' ' outstat];  
rmer = ['!/bin/rm ' outnml ' status' num2str(iiBin)];
fprintf(1,'%s \n',rmer);
eval(rmer)

rmer = ['!/bin/rm ' outstat];
fprintf(1,'%s \n',rmer);
eval(rmer);

if iDoCloud == 100
  rmer = ['!/bin/rm ' outcloud100]; eval(rmer);
end

