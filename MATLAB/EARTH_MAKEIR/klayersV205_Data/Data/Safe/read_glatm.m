% Program read_glatm.m
%
% Read glatm.dat file into matlab
%

% Created: 13 July 2010, Scott Hannon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edit this section as needed

% Comment character
comchar='!';

% Expected number of levels
nlev = 50;

% Expected number of models
nmod = 6;

% Max allowed gas ID
maxID = 80;

% Expected ID of model gases
modID=1:7;

% Surface pressure and altitude
% For kcarta database we want 1100 mb at -689 meters
spres = 1100;
salti = -689;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare matlab arrays

head.pfields=1; % profile
head.ptype=0;   % levels
head.ngas = 0;
head.glist = zeros(maxID,1);
head.gunit = 10*ones(maxID,1); % ppmv
%
prof.plat = zeros(1,nmod);
prof.nlevs = nlev*ones(1,nmod);
prof.plevs = zeros(nlev,nmod);
prof.ptemp = zeros(nlev,nmod);
prof.pnote = char(zeros(80,nmod));
for ii=1:length(modID)
   ig = modID(ii);
   eval(['prof.gas_' int2str(ig) '=zeros(nlev,nmod);'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read file

% Open file
fid = fopen('glatm_9July2010.dat','r');
%fid = fopen('glatm.dat','r');

% Read number of levels = 50
cline = readnextline(fid,comchar);
if (~strcmp(cline(1:2),'50'));
   cline
   error('unexpected number of levels; hardcoded for 50')
end

% Loop over the 6 models
for ii=1:6
   % Read model name
   cline = readnextline(fid,comchar);
   prof.pnote(1:length(cline),ii) = cline;

   % Read model plat
   cline = readnextline(fid,comchar);
    prof.plat(ii) = str2num(cline);

   % Read model altitude (ignore)
   for jj=1:10
      cline = readnextline(fid,comchar);
   end

   % Read model plevs
   for jj=1:10
      cline = readnextline(fid,comchar);
      kk = (jj-1)*5 + (1:5);
      prof.plevs(kk,ii) = str2num(cline);
   end

   % Read model ptemp
   for jj=1:10
      cline = readnextline(fid,comchar);
      kk = (jj-1)*5 + (1:5);
      prof.ptemp(kk,ii) = str2num(cline);
   end

   % Read model density (ignore)
   for jj=1:10
      cline = readnextline(fid,comchar);
   end

   for kk=1:7
      % Read gas ID and profile
      cline = readnextline(fid,comchar);
      ig = str2num(cline(1:2));
      gstr = ['prof.gas_' int2str(ig) '(kk,ii)=str2num(cline);'];
      for jj=1:10
         cline = readnextline(fid,comchar);
         kk = (jj-1)*5 + (1:5);
         eval(gstr);
      end
      if (ii == 1)
         head.ngas = round(head.ngas + 1); % exact integer
         head.glist(head.ngas) = ig;
      end
   end

   % Read modend string
   cline = readnextline(fid,comchar);

end


% Read minigas string
cline = readnextline(fid,comchar);

% Read all remaining gases
domore = 1;
while (domore == 1)
  % Read gas ID or DATAEND
   cline = readnextline(fid,comchar);
   if (strcmp(cline(1),'D'))
      domore = 0;
   else
      ig = str2num(cline(1:2));
      gstr = ['prof.gas_' int2str(ig) '=zeros(nlev,nmod);'];
      eval(gstr);
      gstr = ['prof.gas_' int2str(ig) '(kk,:)=str2num(cline)''*ones(1,nmod);'];
      for jj=1:10
         cline = readnextline(fid,comchar);
         kk = (jj-1)*5 + (1:5);
         eval(gstr);
      end
      head.ngas = round(head.ngas + 1); % exact integer
      head.glist(head.ngas) = ig;
   end
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up
head.glist=head.glist(1:head.ngas);
head.gunit=head.gunit(1:head.ngas);

prof.spres = spres*ones(1,nmod);
prof.salti = salti*ones(1,nmod);

clear ans cline comchar domore fid gstr ig ii jj kk maxID modID nlev nmod
clear stemp salti

disp('done');

%%% end of program %%%
