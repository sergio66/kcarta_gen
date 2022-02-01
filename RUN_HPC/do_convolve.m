function [w,d] = do_convolve(iInstr,ii,iDoRad,outdir);

if nargin == 2
  iDoRad = 3;  %% assume radiances so you need last one ......
  outdir = 'JUNK/';
elseif nargin == 3
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

[djunk,w,caVers] = readkcstd_smart([outdir '/rad.dat' num2str(ii)]);
[mm,nn] = size(djunk);  %% if cloudy calc, could have 890000 x 5 rads
fprintf(1,'read in %s and found size of kcdata = %6i x %6i \n',[outdir '/rad.dat' num2str(ii)],mm,nn)

if iDoRad == 3
  dall = djunk(:,nn);  %%% >>>>>>>> huh, why have this line???? Becuz : if cloudy calc, could have 890000 x 5 rads and only want the last one
else  
  dall = djunk;  %% convolve all of them!
end

% iInstr = input('enter (1) AIRS (2) IASI+IASI NG (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

solzen = -1;
scanang = -1;
satzen = -1;
stemp = -1;
use_this_rtp = 'junk.rp.rtp';
if strcmp(outdir,'/JUNK')
  set_rtp;
end
good = 1;

if iInstr == 1
  clist = 1:2378;
  sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf'; 
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
  airs_convolve_file_numchans  %% gives latest clist/sfile
  [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
  whos dall fKc rKc
  
  %solzen = p.solzen;
  %scanang = p.scanang;
  %satzen = p.satzen;
  %stemp = p.stemp;
  saver = ['save ' outdir '/individual_prof_convolved_kcarta_airs_' num2str(ii) '.mat fKc rKc sfile  use_this_rtp solzen scanang satzen stemp caVers good'];
  eval(saver)
  
elseif iInstr == 2
  %fiasi = instr_chans('iasi')';
  %/asl/matlab/fftconv/s2fconvkc.m

  %[rch,wch]=xfconvkc_serg_iasi(dall,w,'iasi12992','gauss',6);
  %rKcIasi = interp1(wch',rch',fiasi);
  [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall);

  saver = ['save ' outdir '/individual_prof_convolved_kcarta_iasi_' num2str(ii) '.mat fiasi* rKcIasi* use_this_rtp solzen scanang satzen stemp caVers good'];
  eval(saver)
  
elseif iInstr == 3
  error('too lazy to do CrIS LO')
  
elseif iInstr == 4 | iInstr == 14 | iInstr == 124

  iMax = 1;   %% assumes ONLY ONE RADIANCE, no ODS
  mod20 = ceil(iMax/20);

  toptsHI.user_res = 'hires';
  [hi_fcris,hi_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsHI);

  toptsMED.user_res = 'midres';
  [med_fcris,med_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsMED);

  toptsLO.user_res = 'lowres';
  [lo_fcris,lo_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsLO);

  %solzen = p.solzen;
  %scanang = p.scanang;
  %satzen = p.satzen;
  %stemp = p.stemp;
  saver = ['save ' outdir '/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '.mat *rcris_all *fcris use_this_rtp solzen scanang satzen stemp caVers good'];

  if iInstr == 14 | iInstr == 124
    disp('also doing AIRS')
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf';
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
    %saver = ['save ' outdir '/individual_prof_convolved_kcarta_AIRS_crisHI_crisMED_' num2str(ii) '.mat *rcris_all *fcris use_this_rtp solzen scanang satzen stemp caVers good'];
    saver = [saver ' fKc rKc sfile'];
  end

  if iInstr == 24 | iInstr == 124
    disp('also doing IASI')

    %fiasi = instr_chans('iasi')';
    %/asl/matlab/fftconv/s2fconvkc.m

    %[rch,wch]=xfconvkc_serg_iasi(dall,w,'iasi12992','gauss',6);
    %rKcIasi = interp1(wch',rch',fiasi);
    [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall);
    saver = [saver ' fiasi* rKcIasi*'];
  end
  
  eval(saver);
  
end

fprintf(1,'saver = %s ------- using filename %s \n',saver,use_this_rtp)  

