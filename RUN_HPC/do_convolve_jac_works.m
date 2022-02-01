function [w,d] = do_convolve_jac(gg,iInstr,ii,iDoJac,iDoCloud);

if nargin == 4
  iDoCloud == 0;
end

addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

addpath /asl/matlib/fconv
addpath /asl/packages/ccast/source/

addpath /home/sergio/MATLABCODE/FCONV/
addpath /home/sergio/MATLABCODE/FFTCONV/
%addpath /home/sergio/Backup_asl_matlab_Feb2013
addpath /asl/matlab2012/sconv
%addpath /asl/matlab2012/fconv

set_rtp;
% iInstr = input('enter (1) AIRS (2) IASI+IASING (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

if iDoJac == 1
  [dall,w,natmos, numlay, ngases] = readkcjac(['JUNK/jac.dat' num2str(ii)]);
  fjac = ['JUNK/jac.dat' num2str(ii)];
else
  [dall,w] = readkcBasic(['JUNK/jac.dat' num2str(ii) '_COL']);
  fjac = ['JUNK/jac.dat' num2str(ii) '_COL'];
end

fprintf(1,' >>>> gg = %4i iDoCloud = %4i \n',gg,iDoCloud)
if gg ~= 1001 & gg ~= 2346
  fprintf(1,'doing jacobians for gasID %3i \n',gg)
elseif gg == 1001 & iDoCloud <= 0
  fprintf(1,'doing clear sky jacobians for gasID %3i ==> WV 1+101+102+103 ... add together and then O3 \n',1001)
  %% remember do g1,101,102,103  and T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  nnlay = (nn-4)/(4+1+1+1);
  jjuse = 1:nnlay;
  dall1001   = dall(:,jjuse+0*nnlay) + dall(:,jjuse+2*nnlay) + dall(:,jjuse+3*nnlay) + dall(:,jjuse+4*nnlay);
  dallO3     = dall(:,jjuse+1*nnlay);
  dallT      = dall(:,jjuse+5*nnlay);
  dallWgtFcn = dall(:,jjuse+6*nnlay);
  dallSurf   = max(jjuse+6*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
  %whos dall*
  dall = [dall1001 dallO3 dallT dallWgtFcn dallSurf];
  ngases = 1;
elseif gg == 1001 & iDoCloud > 0
  fprintf(1,'doing all sky jacobians for gasID %3i ==> WV 1+101+102+103 ... add together and then O3 \n',1001)
  %% remember do g1,101,102,103  and 201,202 and T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  nnlay = (nn-4)/(4+2+1+1+1);
  jjuse = 1:nnlay;
  dall1001   = dall(:,jjuse+0*nnlay) + dall(:,jjuse+2*nnlay) + dall(:,jjuse+3*nnlay) + dall(:,jjuse+4*nnlay);
  dallO3     = dall(:,jjuse+1*nnlay);
  dallT      = dall(:,jjuse+7*nnlay);
  dallWgtFcn = dall(:,jjuse+8*nnlay);
  dallSurf   = max(jjuse+8*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
  %whos dall*
  dall = [dall1001 dallO3 dallT dallWgtFcn dallSurf];
  ngases = 1;
end

[mm,nn] = size(dall);
if iDoJac == 1
  junk = [mm nn natmos numlay ngases];
  fprintf(1,'jacobian is of size %8i x %3i for %2i atmospheres each with %3i layers and %3i gases \n',junk);
else
  junk = [mm nn];
  fprintf(1,'col jacobian is of size %8i x %3i  \n',junk);
end
d = dall;

if iInstr == 1
  clist = 1:2378;
  sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf'; 
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf'; 
  airs_convolve_file_numchans  %% gives latest clist/sfile
  [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
  whos dall fKc rKc

  if iDoJac == 1
    saver = ['save JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat fKc rKc sfile  use_this_rtp '];
  elseif iDoJac == 100
    saver = ['save JUNK/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat fKc rKc sfile  use_this_rtp '];
  end
  eval(saver)
  
elseif iInstr == 2
  %fiasi = instr_chans('iasi')';
  %/asl/matlab/fftconv/s2fconvkc.m

  iMax = nn;   
  mod20 = ceil(iMax/20);

  for ix = 1 : mod20
    iLay = (1:20) + (ix-1)*20;
    fprintf(1,'chunk %4i of %4i \n',ix,mod20);
    if max(iLay) > iMax
      iLay = iLay(1):iMax;
    end
    if length(iLay) > 1
      %[rch,wch]=xfconvkc_serg_iasi(dall(:,iLay),w,'iasi12992','gauss',6);
      %rKcIasi = interp1(wch',rch',fiasi);
      [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall(:,iLay));
    else
      %[rch,wch]=xfconvkc_serg_iasi(dall(:,iLay),w,'iasi12992','gauss',6);
      %rKcIasi = interp1(wch',rch',fiasi);
      [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall(:,iLay));
    end
    rKcIasi_all(:,iLay) = rKcIasi;
    rKcIasiNG_all(:,iLay) = rKcIasiNG;    
  end
  if iDoJac == 1  
    saver = ['save JUNK/individual_prof_convolved_kcarta_iasi_' num2str(ii) '_jac.mat fiasi* rKcIasi*_all use_this_rtp'];
  else
    saver = ['save JUNK/individual_prof_convolved_kcarta_iasi_' num2str(ii) '_coljac.mat fiasi* rKcIasi*_all use_this_rtp'];
  end
  eval(saver)

elseif iInstr == 3
  error('too lazy to do CrIS LO')
elseif iInstr == 4 | iInstr == 14 | iInstr == 124

  iMax = nn;   
  mod20 = ceil(iMax/20);
  
  toptsHI.user_res = 'hires';  %% Full Spectral Res = FSR
  [hi_fcris,hi_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsHI);

  toptsMED.user_res = 'midres'; %% CHIRP
  [med_fcris,med_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsMED);

  toptsLO.user_res = 'lowres'; %% Nominal Spectral Res = NSR
  [lo_fcris,lo_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsLO);

  %solzen = p.solzen;
  %scanang = p.scanang;
  %satzen = p.satzen;
  %stemp = p.stemp;
  if iDoJac == 1    
    saver = ['save JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_jac.mat *rcris_all *fcris use_this_rtp'];
  else
    saver = ['save JUNK/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_coljac.mat *rcris_all *fcris use_this_rtp'];
  end
  
  if iInstr == 14 | iInstr == 124
    disp('also doing AIRS')
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf';
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
    if iDoJac == 1      
      saver = ['save JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_jac.mat *rcris_all *fcris use_this_rtp'];
      saver = [saver ' fKc rKc sfile'];
    else
      saver = ['save JUNK/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_coljac.mat *rcris_all *fcris use_this_rtp'];
      saver = [saver ' fKc rKc sfile'];
    end
  end

  if iInstr == 24 | iInstr == 124
    disp('also doing IASI')
    %fiasi = instr_chans('iasi')';
    %/asl/matlab/fftconv/s2fconvkc.m
    %[rch,wch]=xfconvkc_serg_iasi(dall,w,'iasi12992','gauss',6);
    %rKcIasi = interp1(wch',rch',fiasi);

    [mm,nn] = size(dall);
    iMax = nn;   
    mod20 = ceil(iMax/20);

    for ix = 1 : mod20
      iLay = (1:20) + (ix-1)*20;
      fprintf(1,'chunk %4i of %4i \n',ix,mod20);
      if max(iLay) > iMax
        iLay = iLay(1):iMax;
      end
      if length(iLay) > 1
        %[rch,wch]=xfconvkc_serg_iasi(dall(:,iLay),w,'iasi12992','gauss',6);
        %rKcIasi = interp1(wch',rch',fiasi);
        [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall(:,iLay));
      else
        %[rch,wch]=xfconvkc_serg_iasi(dall(:,iLay),w,'iasi12992','gauss',6);
        %rKcIasi = interp1(wch',rch',fiasi);
        [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall(:,iLay));
      end
      rKcIasi_all(:,iLay) = rKcIasi;
      rKcIasiNG_all(:,iLay) = rKcIasiNG;    
    end

    clear rKcIasi rKcIasiNG
    saver = [saver ' fiasi* rKcIasi*'];
  end

  eval(saver);
  
end
fprintf(1,'in do_convolve_jac.m iInstr = %4i saver = %s \n',iInstr,saver)
