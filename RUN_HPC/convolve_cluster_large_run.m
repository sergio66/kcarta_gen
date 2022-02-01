addpath /home/sergio/MATLABCODE
addpath /home/sergio/MATLABCODE/CLOUD
addpath /home/sergio/KCARTA/MATLAB
addpath /asl/matlib/h4tools/

% /asl/matlab/fftconv/s2fconvkc.m
% addpath /home/sergio/MATLABCODE/FCONV/
% addpath /home/sergio/MATLABCODE/FFTCONV/

addpath /asl/matlib/fconv
addpath /asl/packages/ccast/source/

addpath /home/sergio/MATLABCODE/FCONV/
addpath /home/sergio/MATLABCODE/FFTCONV/
addpath /home/sergio/Backup_asl_matlab_Feb2013
addpath /asl/matlab2012/sconv
%addpath /asl/matlab2012/fconv

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

clear all

set_rtp;
[h,ha,p,pa] = rtpread(use_this_rtp);
fprintf(1,'%s has %4i profiless \n',use_this_rtp,length(p.stemp))

%iMax = input('Enter number of profiles : ')
iMax = length(p.stemp);

figure(1); clf
iaYes = ones(1,iMax) * -1;
thedir = dir(['JUNK/rad.dat*']);
for ii = 1 : length(thedir)
  zname = thedir(ii).name;
  whichone = str2num(zname(8:end));
  whichsize = thedir(ii).bytes;
  if whichsize >= 3.5e6
    iaYes(whichone) = 2;
  elseif whichsize < 3.5e6
    iaYes(whichone) = 1;
  end
end
plot(iaYes,'.-'); title('iaYes')
pause(0.1);

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

iInstr = input('enter (1) AIRS (2) IASI (3) CrIS lo (4) CrIS hi (14) AIRS and CRIS Hi : ');

fiasi = instr_chans('iasi')';
rKc = []; rKciasi = []; rcris_all = [];

%% grab enough memory for 500 convolutions at a go
iNumChunk = ceil(iMax/500);
for iLoop = 1 : iNumChunk
  iaInd = (1:500) + (iLoop-1)*500;
  if max(iaInd) > iMax
    iaInd = iaInd(1):iMax;
  end
  iLen = length(iaInd);
  dall = zeros(iLen,890000);

  for iii = 1 : iLen
    ii = iaInd(iii);
    if intersect(ii,good)
      [djunk,w,caVers] = readkcstd_smart(['JUNK/rad.dat' num2str(ii)]);
      [mm,nn] = size(djunk);  %% if cloudy calc, could have 890000 x 5 rads
      d = djunk(:,nn);
    else
      d = zeros(1,890000);
    end
    dall(iii,:) = d;
    if mod(iii,100) == 0
      fprintf(1,'    reading in %4i of %4i \n',ii,iMax);
    end
  end

  if iInstr == 1
    clist = 1:2378;
    sfile = '/asl/matlib/srftest/srftables_m140f_withfake_mar08.hdf'; 
    sfile = '/asl/matlab2012/srftest/srftables_m140f_withfake_mar08.hdf'; 
    airs_convolve_file_numchans  %% gives latest clist/sfile
    [fKc,xrKc] = convolve_airs(w,dall,clist,sfile);
    rKc = [rKc xrKc];
    
    solzen = p.solzen;
    scanang = p.scanang;
    satzen = p.satzen;
    stemp = p.stemp;
    saver = ['save JUNK/xconvolved_kcarta_airs.mat fKc rKc sfile  use_this_rtp solzen scanang satzen stemp caVers good'];
    
      if isfield(p,'rcalc') & isfield(p,'robs1')
        figure(1); plot(p.robs1(2333,iaInd),p.rcalc(2333,iaInd),'b',p.robs1(2333,iaInd),xrKc(2333,:),'r',[0 1],[0 1],'k','linewidth',2)
                   xlabel('robs1(2616 cm-1');  ylabel('rcalc(2616 cm-1'); title('(b) orig sarta (r) new kcarta')
        bto = rad2bt(2616,p.robs1(2333,iaInd));
        btc = rad2bt(2616,p.rcalc(2333,iaInd));
       kbtc = rad2bt(2616,xrKc(2333,:));	
	dbt = -5 : 0.2 : +5; dn0 = hist(bto-btc,dbt); dnk = hist(bto-kbtc,dbt);
	figure(2); plot(dbt,dn0,'b',dbt,dnk,'r','linewidth',2);
                   xlabel('2616 : btObs-btCalc');  title('(b) orig sarta (r) new kcarta')
		   grid
        %disp('ret to continue'); pause
	pause(0.1)
      elseif ~isfield(p,'rcalc') & isfield(p,'robs1')
        figure(1); plot(p.robs1(2333,iaInd),xrKc(2333,:),'r',[0 1],[0 1],'k','linewidth',2)
                   xlabel('robs1(2616 cm-1');  ylabel('rcalc(2616 cm-1'); title('(b) orig sarta (r) new kcarta')
        bto = rad2bt(2616,p.robs1(2333,iaInd));
       kbtc = rad2bt(2616,xrKc(2333,:));	
	dbt = -5 : 0.2 : +5; dnk = hist(bto-kbtc,dbt);
	figure(2); plot(dbt,dnk,'r','linewidth',2);
                   xlabel('2616 : btObs-btCalc');  title('(r) new kcarta')
		   grid
        %disp('ret to continue'); pause
	pause(0.1)
      end

  elseif iInstr == 2
    %/asl/matlab/fftconv/s2fconvkc.m
    %error('should be using kc2iasi')
    %[rch,wch]=xfconvkc_serg_iasi(dall,w,'iasi12992','gauss',6);
    %xrKcIasi = interp1(wch',rch',fiasi);
    %rKcIasi = [rKcIasi xrKcIasi];

    [fiasi,rKcIasi,fiasiNG,rKcIasiNG] = iasi_convolve(w,dall);

    solzen = p.solzen;
    scanang = p.scanang;
    satzen = p.satzen;
    stemp = p.stemp;
    saver = ['save JUNK/xconvolved_kcarta_iasi.mat fiasi rKcIasi use_this_rtp solzen scanang satzen stemp caVers good'];

  elseif iInstr == 3
    error('too lazy to do CrIS LO')
    
  elseif iInstr == 4 | iInstr == 14

    mod20 = ceil(iLen/20);
  
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
      [fKc,xrKc] = convolve_airs(w,dall,clist,sfile);
      rKc = [rKc xrKc];
      if isfield(p,'rcalc') & isfield(p,'robs1')
        figure(1); plot(p.robs1(2333,iaInd),p.rcalc(2333,iaInd),'b',p.robs1(2333,iaInd),xrKc(2333,:),'r',[0 1],[0 1],'k','linewidth',2)
                   xlabel('robs1(2616 cm-1');  ylabel('rcalc(2616 cm-1'); title('(b) orig sarta (r) new kcarta')
        bto = rad2bt(2616,p.robs1(2333,iaInd));
        btc = rad2bt(2616,p.rcalc(2333,iaInd));
       kbtc = rad2bt(2616,xrKc(2333,:));	
	dbt = -5 : 0.2 : +5; dn0 = hist(bto-btc,dbt); dnk = hist(bto-kbtc,dbt);
	figure(2); plot(dbt,dn0,'b',dbt,dnk,'r','linewidth',2);
                   xlabel('2616 : btObs-btCalc');  title('(b) orig sarta (r) new kcarta')
		   grid
        %disp('ret to continue'); pause
	pause(0.1)
      elseif ~isfield(p,'rcalc') & isfield(p,'robs1')
        figure(1); plot(p.robs1(2333,iaInd),xrKc(2333,:),'r',[0 1],[0 1],'k','linewidth',2)
                   xlabel('robs1(2616 cm-1');  ylabel('rcalc(2616 cm-1'); title('(b) orig sarta (r) new kcarta')
        bto = rad2bt(2616,p.robs1(2333,iaInd));
       kbtc = rad2bt(2616,xrKc(2333,:));	
	dbt = -5 : 0.2 : +5; dnk = hist(bto-kbtc,dbt);
	figure(2); plot(dbt,dnk,'r','linewidth',2);
                   xlabel('2616 : btObs-btCalc');  title('(r) new kcarta')
		   grid
        %disp('ret to continue'); pause
	pause(0.1)
      end
      saver = ['save JUNK/xconvolved_kcarta_AIRS_crisHI.mat *rcris_all *fcris use_this_rtp solzen scanang satzen stemp caVers good'];
      saver = [saver ' fKc rKc sfile'];
    end
    whos rKc rcris_all
  else
    error('not allowed iInstr')
  end
  eval(saver);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if isfield(p,'rcalc') & isfield(p,'robs1')
        figure(1); plot(p.robs1(2333,good),p.rcalc(2333,good),'b',p.robs1(2333,good),rKc(2333,good),'r',[0 1],[0 1],'k','linewidth',2)
                   xlabel('robs1(2616 cm-1');  ylabel('rcalc(2616 cm-1'); title('(b) orig sarta (r) new kcarta')
        bto = rad2bt(2616,p.robs1(2333,good));
        btc = rad2bt(2616,p.rcalc(2333,good));
       kbtc = rad2bt(2616,rKc(2333,good));	
	dbt = -5 : 0.2 : +5; dn0 = hist(bto-btc,dbt); dnk = hist(bto-kbtc,dbt);
	figure(2); plot(dbt,dn0,'b',dbt,dnk,'r','linewidth',2);
                   xlabel('2616 : btObs-btCalc');  title('(b) orig sarta (r) new kcarta')
		   grid
        %disp('ret to continue'); pause
	pause(0.1)

        g = dogoodchan;
        bto = rad2bt(h.vchan(g),p.robs1(g,good));
        btc = rad2bt(h.vchan(g),p.rcalc(g,good));
        wah = nanstd(bto'-btc'); wah = find(wah < 2);
        g = g(wah);
	
        bto = rad2bt(h.vchan(g),p.robs1(g,good));
        btc = rad2bt(h.vchan(g),p.rcalc(g,good));
       kbtc = rad2bt(h.vchan(g),rKc(g,good));
	 
       figure(3); clf
       plot(h.vchan(g),nanmean(bto'-btc'),'b',h.vchan(g),nanstd(bto'-btc'),'c',...
            h.vchan(g),nanmean(bto'-kbtc'),'r',h.vchan(g),nanstd(bto'-kbtc'),'m')
	     axis([650 2750 -2 +2]); grid
      end

