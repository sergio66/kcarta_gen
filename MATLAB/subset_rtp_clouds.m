function [head, prof]=subset_rtp_clds(headin, profin, glist, clist, plist);

% function [head, prof]=subset_rtp(headin, profin, glist, clist, plist);
%
% Subsets an RTP head & prof structure.
%
% Input:
%    headin = RTP "head" structure
%    profin = RTP "prof" structure
%    glist = (1 x ngas) list of gas IDs to be retained
%    clist = (1 x nchan) list of channel IDs to be retained
%    plist = (1 x nprof) list of profile indices to be retained
%
% Note: if g/c/plist=[], all elements are retained
% 
% Output:
%    head = subsetted "head" structure
%    prof = subsetted "prof" structure
%
% Warning: does not subset non-standard RTP variables! (not even gamnt)
% Note: assumes all profin fields are dimensioned [<whatever> x nprof]
%
% Copied from /asl/matlab/rtptools/subset_rtp.m
% with ciwc, clwc, cc added in


% Created: 13 September 2001, Scott Hannon
% Last update: 1 February 2002, Scott Hannon - add new rtpV103 vars
% Update: 25 June 2002, Scott Hannon - add new rtpV105 vars
% Fix: 22 Oct 2002 Scott Hannon - "calflg" corrected to calflag.
% Update: 15 July 2005, Scott Hannon - re-write checks to make it work
%    with pullchans RTP files. 
% Update: 14 Nov 2008, S.Hannon - update for rtpV201 (removed MW and a
%    few other fields and add some new ones)
% Update: 05 Dec 2008, S.Hannon - add missing clrflag.
% Update: 07 May 2009, S.Hannon - add missing cstemp2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%
% Check headin & profin
%%%%%%%%%%%%%%%%%%%%%%%

% Check for required header fields
if (~isfield(headin,'ngas'))
   disp('Error: input head lacks ngas');
   return
end
ngasin=headin.ngas;
%
if (~isfield(headin,'nchan'))
   disp('Error: input head lacks nchan');
   return
end
nchanin=headin.nchan;
%
% Note: assumes all the 2nd dimension of all profin fields is nprofin
names = fieldnames(profin);
eval( ['d = size(profin.' names{1} ');' ] )
nprofin = d(2);

% Check gas info
if (length(glist) > 0)
   if (~isfield(headin,'glist'))
      disp('Error: input head lacks glist');
      return
   end
   if (~isfield(headin,'gunit'))
      disp('Error: input head lacks gunit');
      return
   end
   %
   ngas=length(glist);
   [c,indg,ib]=intersect(headin.glist,glist);
   if (length(indg) ~= ngas)
      disp('Error: input structures do not contain all gases in glist');
      return
   end
else
   indg=1:ngasin;
   ngas=ngasin;
end

% Check channel info
if (length(clist) > 0)
   if (~isfield(headin,'ichan'))
      disp('Error: input head lacks ichan');
      return
   end
   %
   nchan=length(clist);
   [c,indc,ib]=intersect(headin.ichan,clist);
   if (length(indc) ~= nchan)
      disp('Error: input structures do not contain all channels in clist');
      return
   end
   if (~isfield(profin,'robs1'))
      if (~isfield(profin,'rcalc'))
         disp('Error: no robs1 or rcalc fields for clist to subset');
         return
      end
   end
else
   indc=1:nchanin;
   nchan=nchanin;
end

% Check profile info
if (length(plist) > 0)
   nprof = length(plist);
   % Note: glist is not necessarily sorted
   junk = max( abs( plist - round(plist) ) );
   if (junk > 0)
      disp('Error: plist contains a non-integer index')
      return
   end
   if (min(plist) < 1)
      disp('Error: plist contains a negative or zero index')
      return
   end
   if (max(plist) > nprofin)
      disp('Error: plist contains an index larger than nprofin')
      return
   end
   indp=plist;
else
   indp=1:nprofin;
   nprof=nprofin;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output "head" structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (isfield(headin,'ptype'))
   head.ptype=headin.ptype;
end
if (isfield(headin,'pfields'))
   head.pfields=headin.pfields;
else
   head.pfields=1;
end

if (isfield(headin,'pmin'))
   head.pmin=headin.pmin;
end
if (isfield(headin,'pmax'))
   head.pmax=headin.pmax;
end

head.ngas=ngas; % already tested for existance
if (isfield(headin,'glist'))
   head.glist=headin.glist(indg);
end
if (isfield(headin,'gunit'))
   head.gunit=headin.gunit(indg);
end

if (isfield(headin,'pltfid'))
   head.pltfid=headin.pltfid;
end
if (isfield(headin,'instid'))
   head.instid=headin.instid;
end

head.nchan=nchan; % already tested for existance
if (isfield(headin,'ichan'))
   head.ichan=headin.ichan(indc);
end
if (isfield(headin,'vchan'))
   head.vchan=headin.vchan(indc);
end
if (isfield(headin,'vcmin'))
   head.vcmin=headin.vcmin;
   if (nchanin ~= nchan)
      disp('You should verify head.vcmin is correct for the new channel set')
   end
end
if (isfield(headin,'vcmax'))
   head.vcmax=headin.vcmax;
   if (nchanin ~= nchan)
      disp('You should verify head.vcmax is correct for the new channel set')
   end
end

if (isfield(headin,'iudef'))
   head.iudef=headin.iudef;
end
if (isfield(headin,'itype'))
   head.itype=headin.itype;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output "prof" structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Profile location
if (isfield(profin,'plat'))
   prof.plat=profin.plat(indp);
end
if (isfield(profin,'plon'))
   prof.plon=profin.plon(indp);
end
if (isfield(profin,'ptime'))
   prof.ptime=profin.ptime(indp);
end


% Surface
if (isfield(profin,'stemp'))
   prof.stemp=profin.stemp(indp);
end
if (isfield(profin,'salti'))
   prof.salti=profin.salti(indp);
end
if (isfield(profin,'spres'))
   prof.spres=profin.spres(indp);
end
if (isfield(profin,'landfrac'))
   prof.landfrac=profin.landfrac(indp);
end
if (isfield(profin,'landtype'))
   prof.landtype=profin.landtype(indp);
end
%
if (isfield(profin,'wspeed'))
   prof.wspeed=profin.wspeed(indp);
end
%
if (isfield(profin,'nemis'))
   prof.nemis=profin.nemis(indp);
   if (isfield(profin,'efreq'))
      prof.efreq=profin.efreq(:,indp);
   end
   if (isfield(profin,'emis'))
      prof.emis=profin.emis(:,indp);
   end
   if (isfield(profin,'rho'))
      prof.rho=profin.rho(:,indp);
   end
end


% Profiles
if (isfield(profin,'nlevs'))
   prof.nlevs=profin.nlevs(indp);
end
if (isfield(profin,'plevs'))
   prof.plevs=profin.plevs(:,indp);
end
if (isfield(profin,'palts'))
   prof.palts=profin.palts(:,indp);
end
if (isfield(profin,'ptemp'))
   prof.ptemp=profin.ptemp(:,indp);
end
for ig=1:ngas
   gstr=int2str(head.glist(ig));
   if (isfield(profin,['gas_' gstr]))
      eval(['prof.gas_' gstr '=profin.gas_' gstr '(:,indp);']);
   else
      disp(['Warning: expected profin field gas_' gstr ' does not exist'])
      return
   end
end
%
if (isfield(profin,'gtotal'))
   prof.gtotal=profin.gtotal(indg,indp);
end
if (isfield(profin,'gxover'))
   prof.gxover=profin.gxover(indg,indp);
end
if (isfield(profin,'txover'))
   prof.txover=profin.txover(indp);
end
if (isfield(profin,'co2ppm'))
   prof.co2ppm=profin.co2ppm(indp);
end

% Clear flag
if (isfield(profin,'clrflag'))
   prof.clrflag=profin.clrflag(indp);
end

% Clouds
if (isfield(profin,'ctype'))
   prof.ctype=profin.ctype(indp);
end
if (isfield(profin,'cfrac'))
   prof.cfrac=profin.cfrac(indp);
end
if (isfield(profin,'cemis'))
   prof.cemis=profin.cemis(:,indp);
end
if (isfield(profin,'crho'))
  prof.crho=profin.crho(:,indp);
end
if (isfield(profin,'cprtop'))
   prof.cprtop=profin.cprtop(indp);
end
if (isfield(profin,'cprbot'))
   prof.cprbot=profin.cprbot(indp);
end
if (isfield(profin,'cngwat'))
   prof.cngwat=profin.cngwat(indp);
end
if (isfield(profin,'cpsize'))
   prof.cpsize=profin.cpsize(indp);
end
if (isfield(profin,'cstemp'))
   prof.cstemp=profin.cstemp(indp);
end
%
if (isfield(profin,'ctype2'))
   prof.ctype2=profin.ctype2(indp);
end
if (isfield(profin,'cfrac2'))
   prof.cfrac2=profin.cfrac2(indp);
end
if (isfield(profin,'cemis2'))
   prof.cemis2=profin.cemis2(:,indp);
end
if (isfield(profin,'crho2'))
  prof.crho2=profin.crho2(:,indp);
end
if (isfield(profin,'cprtop2'))
   prof.cprtop2=profin.cprtop2(indp);
end
if (isfield(profin,'cprbot2'))
   prof.cprbot2=profin.cprbot2(indp);
end
if (isfield(profin,'cngwat2'))
   prof.cngwat2=profin.cngwat2(indp);
end
if (isfield(profin,'cpsize2'))
   prof.cpsize2=profin.cpsize2(indp);
end
if (isfield(profin,'cstemp2'))
   prof.cstemp2=profin.cstemp2(indp);
end
if (isfield(profin,'cfrac12'))
   prof.cfrac12=profin.cfrac12(indp);
end
if (isfield(profin,'cc'))
   prof.cc=profin.cc(:,indp);
end
if (isfield(profin,'ciwc'))
   prof.ciwc=profin.ciwc(:,indp);
end
if (isfield(profin,'clwc'))
   prof.clwc=profin.clwc(:,indp);
end


% Radiance viewing parameters
if (isfield(profin,'pobs'))
   prof.pobs=profin.pobs(indp);
end
if (isfield(profin,'zobs'))
   prof.zobs=profin.zobs(indp);
end
if (isfield(profin,'upwell'))
   prof.upwell=profin.upwell(indp);
end
if (isfield(profin,'scanang'))
   prof.scanang=profin.scanang(indp);
end
if (isfield(profin,'satzen'))
   prof.satzen=profin.satzen(indp);
end
if (isfield(profin,'satazi'))
   prof.satazi=profin.satazi(indp);
end
%
if (isfield(profin,'solzen'))
   prof.solzen=profin.solzen(indp);
end
if (isfield(profin,'solazi'))
   prof.solazi=profin.solazi(indp);
end
if (isfield(profin,'sundist'))
   prof.sundist=profin.sundist(indp);
end
if (isfield(profin,'glint'))
   prof.glint=profin.glint(indp);
end
%
if (isfield(profin,'rlat'))
   prof.rlat=profin.rlat(indp);
end
if (isfield(profin,'rlon'))
   prof.rlon=profin.rlon(indp);
end
if (isfield(profin,'rfill'))
   prof.rfill=profin.rfill(indp);
end
if (isfield(profin,'rtime'))
   prof.rtime=profin.rtime(indp);
end
if (isfield(profin,'findex'))
   prof.findex=profin.findex(indp);
end
if (isfield(profin,'atrack'))
   prof.atrack=profin.atrack(indp);
end
if (isfield(profin,'xtrack'))
   prof.xtrack=profin.xtrack(indp);
end
if (isfield(profin,'ifov'))
   prof.ifov=profin.ifov(indp);
end


% Radiance
if (isfield(profin,'robs1'))
   prof.robs1=profin.robs1(indc,indp);
end
if (isfield(profin,'calflag'))
   prof.calflag=profin.calflag(indc,indp);
end
if (isfield(profin,'robsqual'))
   prof.robsqual=profin.robsqual(indp);
end
if (isfield(profin,'freqcal'))
   prof.freqcal=profin.freqcal(indp);
end
%
if (isfield(profin,'rcalc'))
   prof.rcalc=profin.rcalc(indc,indp);
end


% User Defined data
if (isfield(profin,'pnote'))
   prof.pnote=profin.pnote(:,indp);
end
if (isfield(profin,'udef'))
   prof.udef=profin.udef(:,indp);
end
if (isfield(profin,'iudef'))
   prof.iudef=profin.iudef(:,indp);
end
if (isfield(profin,'itype'))
   prof.itype=profin.itype(indp);
end

%%% end of file %%%
