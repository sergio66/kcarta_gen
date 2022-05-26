function [w,d] = do_convolve_jac(gg,iInstr,ii,iDoJac,iDoCloud,outdir);

if nargin == 4
  iDoCloud = 0;
  outdir = 'JUNK/';
elseif nargin == 5
  outdir = 'JUNK/';
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

use_this_rtp = 'junk.rp.rtp';
if strcmp(outdir,'/JUNK')
  set_rtp;
end

% iInstr = input('enter (1) AIRS (2) IASI+IASING (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

%% ~/KCARTA/SRCv1.21_f90/Revisions.txt
%6/25/19        Changed the column jac mult    from 0.1 to 0.001
%               Changed the column temp offset from 1.0 to 0.01

if iDoJac == 1
  fjac = ['JUNK/jac.dat' num2str(ii)];
  [dall,w,natmos, numlay, ngases] = readkcjac([outdir '/jac.dat' num2str(ii)]);

elseif iDoJac == 100 & iDoCloud <= 0
  %% g2,4,5,6,51,52,T(z),ST
  fjac = [outdir '/jac.dat' num2str(ii) '_COL'];
  [radall,w] = readkcstd([outdir '/rad.dat' num2str(ii)]);    
  [dall,w]   = readkcBasic([outdir '/jac.dat' num2str(ii) '_COL']);  

  [mmm,nnnj] = size(dall);
  [mmm,nnnr] = size(radall);
  radall = radall(:,nnnr) * ones(1,nnnj);
  dall = rad2bt(w,dall)-rad2bt(w,radall);

  %% trace gas jacs
  iii = 1 : nnnj-2
  dall(:,iii) = dall(:,iii)/log(1.001);

  %% T(z) and ST jac
  iii = nnnj-1 : nnnj
  dall(:,iii) = dall(:,iii)/0.01;

elseif iDoJac == 100 & iDoCloud == 1
  %% g2,4,5,6,T(z),ST
  fjac = [outdir '/jac.dat' num2str(ii) '_COL'];
  [radall,dall,w] = readkcPCLSAM_coljac([outdir '/rad.dat' num2str(ii)],[outdir '/jac.dat' num2str(ii) '_COL']);

  [mmm,nnnj] = size(dall);
  [mmm,nnnr] = size(radall);
  radall = radall(:,nnnr) * ones(1,nnnj);
  dall = rad2bt(w,dall)-rad2bt(w,radall);

  %% trace gas jacs
  iii = 1 : nnnj-2
  dall(:,iii) = dall(:,iii)/log(1.001);

  %% T(z) and ST jac
  iii = nnnj-1 : nnnj
  dall(:,iii) = dall(:,iii)/0.01;

else
  error('huh? Need iDoJac == 1 or 100')
end

whos w dall

fprintf(1,' >>>> gg = %4i iDoCloud = %4i \n',gg,iDoCloud)
if gg ~= 1001 & gg ~= 2346 & gg ~= 5912 & iDoJac == 1
  fprintf(1,'doing jacobians for gasID %3i \n',gg)
elseif gg == 1001 & iDoJac == 1
  if iDoCloud <= 0
    fprintf(1,'doing clear sky jacobians for gasID %3i ==> WV 1+101+102+103 ... add together and then O3 \n',1001)
  elseif iDoCloud == 1
    fprintf(1,'doing allsky jacobians for gasID %3i ==> WV 1+101+102+103 ... add together and then O3 \n',1001)
  end
  %% remember do g1,101,102,103  and T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  if iDoCloud <= 0
    nnlay = (nn-4)/(4+1+1+1);
  else
    nnlay = (nn-4)/(4+1+1+1+1+1);
  end
  fprintf(1,'   do_convolve_jac.m : gg=%4i iDoCloud=%2i     iDoJac=%2i mm,nn=%10i %4i     nnlay=%8.6f \n',gg,iDoCloud,iDoJac,mm,nn,nnlay);
  jjuse = 1:nnlay;
  dall1001   = dall(:,jjuse+0*nnlay) + dall(:,jjuse+2*nnlay) + dall(:,jjuse+3*nnlay) + dall(:,jjuse+4*nnlay);
  dallO3     = dall(:,jjuse+1*nnlay);
  dallT      = dall(:,jjuse+5*nnlay);
  dallWgtFcn = dall(:,jjuse+6*nnlay);
  dallSurf   = max(jjuse+6*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
  %whos dall*
  dall = [dall1001 dallO3 dallT dallWgtFcn dallSurf];
  ngases = 2;

  addpath /home/sergio/MATLABCODE
  [fc,qc] = quickconvolve(w,dall1001,1,1);
  figure(1); plot(fc,nansum(qc')); title('WV coljac'); pause(0.1);

elseif gg == 2346 & iDoJac == 1
  if iDoCloud <= 0
    fprintf(1,'doing clear sky jacobians for gasID %3i ==> 2 3 4 6 51 52 \n',2346)
  elseif iDoCloud == 1
    fprintf(1,'doing clear sky jacobians for gasID %3i ==> 2 3 4 6 51 52 \n',2346)
  end
  %% remember also do T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  if iDoCloud <= 0
    nnlay = (nn-4)/(8);
  else
    nnlay = (nn-4)/(10); %% should this also be 8?
  end
  fprintf(1,'   do_convolve_jac.m : gg=%4i iDoCloud=%2i     iDoJac=%2i mm,nn=%10i %4i     nnlay=%8.6f \n',gg,iDoCloud,iDoJac,mm,nn,nnlay);
  jjuse = 1:nnlay;
  dall2      = dall(:,jjuse+0*nnlay);
  dall3      = dall(:,jjuse+1*nnlay);
  dall4      = dall(:,jjuse+2*nnlay);
  dall6      = dall(:,jjuse+3*nnlay);
  dall51     = dall(:,jjuse+4*nnlay);
  dall52     = dall(:,jjuse+5*nnlay);
  dallT      = dall(:,jjuse+6*nnlay);
  dallWgtFcn = dall(:,jjuse+7*nnlay);
  dallSurf   = max(jjuse+8*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
  %whos dall*
  dall = [dall2 dall3 dall4 dall6 dall51 dall52 dallT dallWgtFcn dallSurf];
  ngases = 6;

  addpath /home/sergio/MATLABCODE
  [fc,qc] = quickconvolve(w,dall2,1,1);
  figure(1); plot(fc,nansum(qc')); title('CO2 coljac'); pause(0.1);

elseif gg == 5912 & iDoJac == 1
  if iDoCloud <= 0
    fprintf(1,'doing clear sky jacobians for gasID %3i ==> 5 9 11 12 61 103 \n',5912)
  elseif iDoCloud == 1
    fprintf(1,'doing clear sky jacobians for gasID %3i ==> 5 9 11 12 61 103 \n',5912)
  end
  %% remember also do T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  if iDoCloud <= 0
    nnlay = (nn-4)/(8);
  else
    nnlay = (nn-4)/(10); %% should this also be 8?
  end
  fprintf(1,'   do_convolve_jac.m : gg=%4i iDoCloud=%2i     iDoJac=%2i mm,nn=%10i %4i     nnlay=%8.6f \n',gg,iDoCloud,iDoJac,mm,nn,nnlay);
  jjuse = 1:nnlay;
  dall5      = dall(:,jjuse+0*nnlay);
  dall9      = dall(:,jjuse+1*nnlay);
  dall11      = dall(:,jjuse+2*nnlay);
  dall12     = dall(:,jjuse+3*nnlay);
  dall61     = dall(:,jjuse+4*nnlay);
  dall103    = dall(:,jjuse+5*nnlay);
  dallT      = dall(:,jjuse+6*nnlay);
  dallWgtFcn = dall(:,jjuse+7*nnlay);
  dallSurf   = max(jjuse+8*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
  %whos dall*
  dall = [dall5 dall9 dall11 dall12 dall61 dall103 dallT dallWgtFcn dallSurf];
  ngases = 6;

  addpath /home/sergio/MATLABCODE
  [fc,qc] = quickconvolve(w,dall5,1,1);
  figure(1); plot(fc,nansum(qc')); title('CO coljac'); pause(0.1);

%%elseif gg == 1001 & iDoCloud > 0
%%  fprintf(1,'doing all sky jacobians for gasID %3i ==> WV 1+101+102+103 ... add together and then O3 \n',1001)
%%  %% remember do g1,101,102,103  and 201,202 and T and WgtFcn  1:nn  and then 4 surface
%%  [mm,nn] = size(dall);
%%  nnlay = (nn-4)/(4+2+1+1+1);
%%  jjuse = 1:nnlay;
%%  dall1001   = dall(:,jjuse+0*nnlay) + dall(:,jjuse+2*nnlay) + dall(:,jjuse+3*nnlay) + dall(:,jjuse+4*nnlay);
%%  dallO3     = dall(:,jjuse+1*nnlay);
%%  dallT      = dall(:,jjuse+7*nnlay);
%%  dallWgtFcn = dall(:,jjuse+8*nnlay);
%%  dallSurf   = max(jjuse+8*nnlay); dallSurf = dallSurf+1:nn; dallSurf = dall(:,dallSurf);
%%  %whos dall*
%%  dall = [dall1001 dallO3 dallT dallWgtFcn dallSurf];
%%  ngases = 2;

elseif gg == 2346 & iDoJac == 100
  if iDoCloud <= 0
    fprintf(1,'doing clear sky COL jacobians for gasID %3i ==> 2 3 4 6 51 52 \n',2346)
  elseif iDoCloud == 1
    fprintf(1,'doing clear sky COL jacobians for gasID %3i ==> 2 3 4 6 \n',2346)
  end
  %% remember also do T and WgtFcn  1:nn  and then 4 surface
  [mm,nn] = size(dall);
  nnlay = (nn-4)/(6);
  fprintf(1,'   do_convolve_jac.m : COL JAC gg=%4i iDoCloud=%2i     iDoJac=%2i mm,nn=%10i %4i     nnlay=%8.6f \n',gg,iDoCloud,iDoJac,mm,nn,nnlay);
  jjuse = 1:nnlay;
  dall2      = dall(:,1);
  dall3      = dall(:,2);
  dall4      = dall(:,3);
  dall6      = dall(:,4);

  if iDoCloud > 0
    dallT      = dall(:,5);
    dallST     = dall(:,6);
    dall = [dall2 dall3 dall4 dall6 dallT dallST];
    ngases = 4;
  elseif iDoCloud < 0
    dall51     = dall(:,5);
    dall52     = dall(:,6);
    dallT      = dall(:,7);
    dallST     = dall(:,8);
    dall = [dall2 dall3 dall4 dall6 dall51 dall52 dallT dallST];
    ngases = 6;
  end

  addpath /home/sergio/MATLABCODE
  [fc,qc] = quickconvolve(w,dall2,1,1);
  figure(1); plot(fc,qc); title('CO2 coljac');

else
  error('kjsklsjglkjgskjgs')
end

whos w dall

[mm,nn] = size(dall);
if iDoJac == 1
  junk = [mm nn natmos numlay ngases];
  fprintf(1,'jacobian is of size %8i x %3i for %2i atmospheres each with %3i layers and %3i gases \n',junk);
elseif iDoJac == 100
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
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '_jac.mat fKc rKc sfile  use_this_rtp '];
  elseif iDoJac == 100
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '_coljac.mat fKc rKc sfile  use_this_rtp '];
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
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_iasi_' num2str(ii) '_jac.mat fiasi* rKcIasi*_all use_this_rtp'];
  else
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_iasi_' num2str(ii) '_coljac.mat fiasi* rKcIasi*_all use_this_rtp'];
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
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_jac.mat *rcris_all *fcris use_this_rtp'];
  else
    saver = ['save ' outdir '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_coljac.mat *rcris_all *fcris use_this_rtp'];
  end
  
  if iInstr == 14 | iInstr == 124
    disp('also doing AIRS')
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf';
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
    if iDoJac == 1      
      saver = ['save ' outdir '/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_jac.mat *rcris_all *fcris use_this_rtp'];
      saver = [saver ' fKc rKc sfile'];
    else
      saver = ['save ' outdir '/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '_coljac.mat *rcris_all *fcris use_this_rtp'];
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