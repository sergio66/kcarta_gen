%%%%% ####### seq 1 30000 > filelist.hpc

iMax = 10000;
iMax = input('Enter how many profiles you need to process : ');
iWhich = input('Enter (-1) for all files from 1 : N or (+1) to check all files between 1: N : ');

fid = fopen('filelist_hpc','w');

if iWhich < 0
  for ii = 1 : iMax
    fprintf(fid,'%5i \n',ii);
  end

 else
  for ii = 1 : iMax
    filename = ['JUNK/individual_prof_convolved_kcarta_crisHI_' num2str(ii) '.mat'];
    if ~exist(filename)
      fprintf(fid,'%5i \n',ii);
    end
  end
end

fclose(fid);
