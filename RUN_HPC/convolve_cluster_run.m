addpath /home/sergio/MATLABCODE
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

clear all

set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %4i profiless \n',use_this_rtp,length(p.stemp))

%iMax = input('Enter number of profiles : ')
iMax = length(p.stemp);

for ii = 1 : iMax
  if exist(['JUNK/rad.dat' num2str(ii)])
    thedir = dir(['JUNK/rad.dat' num2str(ii)]);
    if thedir.bytes > 3.5e6
      %% correct size
      iaYes(ii) = 2;
    else
      %% not completely done
      iaYes(ii) = 1;
    end
  else
    iaYes(ii) = -1;
  end
end

good = 1 : iMax;
small = find(iaYes == 1);
good = setdiff(good,small);
if length(small) > 0
  small
  disp('above files incomplete')
  iYes = input('continue (+1)Y (-1)N : ');
  if iYes < 0
    error('stopping')
  end
end

bad = find(iaYes < 0);
good = setdiff(good,bad);
if length(bad) > 0
  bad
  disp('above files not found')
  iYes = input('continue (+1)Y (-1)N : ');
  if iYes < 0
    error('stopping')
  end
end

dall = zeros(iMax,890000);

for ii = 1 : iMax
  if intersect(ii,good)
    [djunk,w,caVers] = readkcstd_smart(['JUNK/rad.dat' num2str(ii)]);
    [mm,nn] = size(djunk);  %% if cloudy calc, could have 890000 x 5 rads
    d = djunk(:,nn);    
  else
    d = zeros(1,890000);
  end
  dall(ii,:) = d;
  if mod(ii,100) == 0
    fprintf(1,'    reading in %4i of %4i \n',ii,iMax);
  end
end

iInstr = input('enter (1) AIRS (2) IASI (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

if iInstr == 1
  clist = 1:2378;
  sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf'; 
  sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf'; 
  airs_convolve_file_numchans  %% gives latest clist/sfile
  [fKc,rKc] = convolve_airs(w,dall,clist,sfile);

  solzen = p.solzen;
  scanang = p.scanang;
  satzen = p.satzen;
  stemp = p.stemp;
  save JUNK/xconvolved_kcarta_airs.mat fKc rKc sfile  use_this_rtp solzen scanang satzen stemp caVers good
elseif iInstr == 2
  error('too lazy to do IASI')
elseif iInstr == 3
  error('too lazy to do CrIS LO')
elseif iInstr == 4 | iInstr == 14

  mod20 = ceil(iMax/20);

  toptsHI.user_res = 'hires';
  [hi_fcris,hi_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsHI);

  toptsMED.user_res = 'midres';
  [med_fcris,med_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsMED);

  toptsLO.user_res = 'lowres';
  [lo_fcris,lo_rcris_all] = convolve_cris_all_chooseres(w,dall,toptsLO);

  solzen = p.solzen;
  scanang = p.scanang;
  satzen = p.satzen;
  stemp = p.stemp;
  saver = ['save JUNK/xconvolved_kcarta_crisHI.mat *rcris_all *fcris use_this_rtp solzen scanang satzen stemp caVers good'];

  if iInstr == 14
    disp('also doing AIRS')
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf';
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf';
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,rKc] = convolve_airs(w,dall,clist,sfile);
    saver = ['save JUNK/xconvolved_kcarta_AIRS_crisHI.mat *rcris_all *fcris use_this_rtp solzen scanang satzen stemp caVers good'];
    saver = [saver ' fKc rKc sfile'];
  end
  eval(saver);
  
end

  

