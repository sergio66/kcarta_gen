function [] = do_convolve_all(w,dall,iInstr);

addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

addpath /asl/matlib/fconv
addpath /asl/packages/ccast/source/

addpath /home/sergio/MATLABCODE/FCONV/
addpath /home/sergio/MATLABCODE/FFTCONV/
addpath /home/sergio/Backup_asl_matlab_Feb2013
addpath /asl/matlab2012/sconv
%addpath /asl/matlab2012/fconv

%[dall,w,natmos, numlay, ngases] = readkcjac(['JUNK/jac.dat' num2str(ii)]);

%set_rtp;
% iInstr = input('enter (1) AIRS (2) IASI+IASING (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

[mm,nn] = size(dall);
junk = [mm nn];
fprintf(1,'jacobian is of size %8i x %3i\n',junk);

if iInstr == 1
  clist = 1:2378;
  sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf'; 
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf'; 
  airs_convolve_file_numchans  %% gives latest clist/sfile
  [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
  whos dall fKc rKc
  
  saver = ['save JUNK/convolve_all.mat fKc rKc sfile '];
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
    rKcIasi_all(:,iLay)   = rKcIasi;
    rKcIasiNG_all(:,iLay) = rKcIasiNG;    
  end
  
  clear rKcIasi rKcIasiNG
  saver = ['save JUNK/convolve_all.mat fiasi* rKcIasi*_all'];
  eval(saver)

elseif iInstr == 3
  error('too lazy to do CrIS LO')
elseif iInstr == 4 | iInstr == 14

  iMax = nn;   
  mod20 = ceil(iMax/20);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%
  toptsHI.user_res = 'hires';
  [hi_fcris,hi_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsHI);

  toptsMED.user_res = 'midres';
  [med_fcris,med_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsMED);

  toptsLO.user_res = 'lowres';
  [lo_fcris,lo_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsLO);

  %%%%%%%%%%%%%%%%%%%%%%%%%

  %solzen = p.solzen;
  %scanang = p.scanang;
  %satzen = p.satzen;
  %stemp = p.stemp;
  saver = ['save JUNK/convolve_all.mat *rcris_all *fcris'];

  if iInstr == 14
    disp('also doing AIRS')
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf';
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
    saver = ['save JUNK/convolve_all.mat rcris_all fcris'];
    saver = [saver ' fKc rKc sfile'];
  end
  eval(saver);
  
end

  

