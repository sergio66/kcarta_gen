addpath /home/sergio/KCARTA/MATLAB
wdis = zeros(89*10000,1);

iDISopt = input('Enter iDISopt (11/Default) is just one profile, 12 is multiple profiles : ');
if length(iDISopt) == 0
  iDISopt = 11;
end

iaProf = input('Enter rtp profile number list ([defulat = nothing, becomes profile 1], one number for only a sigle profile, two numbers for (Start: Stop) : ');
if length(iaProf) == 0
  iaProf = 1;
end

if length(iaProf) == 1
  iaList = iaProf;
elseif length(iaProf) == 2
  iaList = iaProf(1):iaProf(2);
end

raddis = nan(89*10000,length(iaList));

iaaFound = zeros(length(iaList),89);
for iPP = 1 : length(iaList)
  iProf = iaList(iPP);
  iCnt = 0;
  for ii = 1 : 89
    ind = (1:10000) + (ii-1)*10000;
    f1 = 605 + (ii-1)*25;
    fname = ['JUNK/rad.dat' num2str(iProf) '_' num2str(f1,'%04d') '_iDISopt_' num2str(iDISopt,'%02i') ];
    %fprintf(1,'%4i %s %3i \n',f1,fname,exist(fname));
    if exist(fname)
      thedir = dir(fname);
      if thedir.bytes > 40000
        iCnt = iCnt + 1;
        iaaFound(iPP,ii) = 1;
        [dx,wx] = readkcstd(fname);
        wdis(ind) = wx;
        raddis(ind,iPP) = dx;
        %plot(wx,rad2bt(wx,dx)); pause
      end
    end
  end
  
  fprintf(1,'profile %4i read in %2i of 89 DISORT chunks into [wdis,raddis] \n',iProf,iCnt)
  if iCnt > 0 & length(iaList) == 1
    plot(wdis,rad2bt(wdis,raddis));
  end
  if iCnt < 89 & length(iaList) == 1
    disp('missing chunks')
    boo = find(iaFound == 0);
  end
end
figure(1); plot(wdis,nanmean(rad2bt(wdis,raddis),2));
%sum(iaaFound(:))/89/length(iaList)
if length(iaList) > 1
  figure(2); imagesc(1:89,iaList,iaaFound); colorbar; ylabel('Profile #'); xlabel('kCARTA chunk')
end

if sum(iaaFound(:))/89/length(iaList) == 1
  fprintf(1,'found all %5i files %5i\n',89*length(iaList),sum(iaaFound(:)))
elseif sum(iaaFound(:))/89/length(iaList) < 1
  fprintf(1,'need to find %5i files but only found %5i \n',89*length(iaList),sum(iaaFound(:)))
  iDelete = input('Find and delete bad files (-1/+1) : ');
  if iDelete > 0
    fid = fopen('missing_disort.sc','w');
    [mm,nn] = find(iaaFound == 0);

    for ii = 1 : length(mm)
      iProf = mm(ii) + (iaList(1)-1);
      f1 = 605 + (nn(ii)-1)*25;
      fname = ['JUNK/rad.dat' num2str(iProf) '_' num2str(f1,'%04d')];
      if exist(fname)
        thedir = dir(fname);
        thedir = thedir.bytes;
      else
        thedir = 0;
      end       
      fprintf(1,'%4i %s %4i \n',ii,fname,thedir);

      launcher = ['sbatch -p high_mem    --array=' num2str(nn(ii)) ' sergio_matlab_jobB.sbatch 11 ' num2str(iProf)];
      fprintf(fid,'%s \n',launcher);
      if exist(fname)
        rmer = ['!/bin/rm ' fname];
        eval(rmer);
      end

      fnameCLD = ['JUNK/rad.dat' num2str(iProf) '_' num2str(f1,'%04d') '_CLD'];
      if exist(fnameCLD)
        rmer = ['!/bin/rm ' fnameCLD];
        eval(rmer);
      end

    end    %% loop over missing files    
    fclose(fid);
    disp('can look at missing_disort.sc')
  end      %% if iDelete > 0
end        %%  start loop of deleting missing files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iSave = input('Save the data??? (-1/+1) : ');
if iSave > 0
  %badprofx  = mm + (iaList(1)-1);
  %badchunkx = nn;
  saver = ['save disort_profs_iDISopt_' num2str(iDISopt,'%02i') '_' num2str(iaList(1),'%05d') '_' num2str(iaList(end),'%05d') '.mat wdis raddis iaList'];
  eval(saver)
end
