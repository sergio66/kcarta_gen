function op_rtp_to_txt(filename, profnum, head, prof);

% function op_rtp_to_txt(filename, profnum, head, prof);
%
% Take a KLAYERS "op" type output file in RTP format and dump
% it to a text file in the format used by KCARTA v1.07.
%
% Input:
%    filename : (string) name of text file to create
%    profnum  : (1 x 1) index of desired profile in "prof"
%    head     : (RTP head structure)
%    prof     : (RTP prof structure)
%
% Output: none (except for the text file)
%

% Created 30 August 2001, Scott Hannon
% Update: 11 Nov 2009, S.Hannon - calc plays from plevs
% Update: 13 Nov 2009, S.Hannon - bug fix: add fclose(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nlevs=prof.nlevs(profnum);
nlays=nlevs-1;
nlaytot=100*head.ngas;
plevs = prof.plevs(:,profnum);
i2 = 2:nlevs;
i1 = 1:nlays;
plays = ( plevs(i2) - plevs(i1) )./log( plevs(i2)./plevs(i1) );

% Open output file
fid=fopen(filename,'w');

% Write total number of layers
fprintf(fid,'% 5d\n', nlaytot);

% Loop over the gases 
for igas=1:head.ngas
   eval(['amount=prof.gas_' int2str(head.glist(igas)) '(1:nlays,profnum);']);
   temp=prof.ptemp(1:nlays,profnum);
   zlo=prof.palts(1:nlays,profnum);
   zhi=prof.palts(2:nlevs,profnum);
   zmean=(zlo + zhi)/2000; % mean altitude in km
   dz=abs(zhi - zlo);
   % Calc partial pressure
   pp = amount .* temp ./ ( 2.6867775E+19 * 273.15 * 100*dz );
   % Convert amount from molecules/cm^2 to kilomoles/cm^2
   amount=amount/6.02214199E+26;
   for ilev=100:-1:nlays+1
      fprintf(fid, ...
      '% 3d % 11.4e % 8.3f % 8.3f % 11.4e % 11.4e % 11.4e % 8.3f\n',...
      head.glist(igas), 0, 273.15, 0, ...
      0, 0, 0, 0);
   end
   for ilev=nlays:-1:1
      fprintf(fid, ...
      '% 3d % 11.4e % 8.3f % 8.3f % 11.4e % 11.4e % 11.4e % 8.3f\n',...
      head.glist(igas), amount(ilev), temp(ilev), 0, ...
      plays(ilev)/1013.25, 0, pp(ilev), zmean(ilev));
   end
end

junk = fclose(fid);

%%% end of function %%%
