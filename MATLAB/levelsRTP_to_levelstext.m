function [] = levelsRTP_to_levelstext(iProfile,fRTP,fTXT)

%{
test
levelsRTP_to_levelstext(1,'junk49.ip.rtp');
%}

%% also see levelstext_to_levelsRTP.m

if nargin < 2
  error('need at least two inputs (iProfile,fRTP)')
elseif nargin == 2
  fTXT = ['levelsprofile_text' num2str(iProfile) '.prf'];
end

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools

[h,ha,p,pa] = rtpread(fRTP);
[h,p] = subset_rtp(h,p,[],[],iProfile);

fid = fopen(fTXT,'w');

caStr = ['see /home/sergio/KCARTA/MATLAB/levelsRTP_to_levelstext.m rtp=' fRTP ' iProfile=' num2str(iProfile)];
if length(caStr) > 80
  caStr = caStr(1:80);
end
fprintf(fid,'%s\n',caStr);

fprintf(fid,'%3i \n',p.nlevs);
fprintf(fid,'%8.3f %8.3f %8.3f \n',p.spres,p.stemp,p.salti);
if intersect(h.glist,2)
  %% already have CO2
  fprintf(fid,'%4i   %8.3f %8.3f \n',-9999,0.0,0.0);
else
  xstr = date; xstr(end-3:end)
  fprintf(fid,'%s   %8.3f %8.3f \n',xstr,0.0,0.0);
end

fprintf(fid,'%3i \n',h.ngas);

for ii = 1 : h.ngas-1
  fprintf(fid,'%3i ',h.glist(ii));
end
ii = h.ngas;
fprintf(fid,'%3i \n',h.glist(ii));

for ii = 1 : h.ngas-1
  fprintf(fid,'%3i ',h.gunit(ii));
end
ii = h.ngas;
fprintf(fid,'%3i \n',h.gunit(ii));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ll = 1 : p.nlevs
  fprintf(fid,'%8.6e %8.2f ',p.plevs(ll),p.ptemp(ll));
  for ii = 1 : h.ngas
    gid = h.glist(ii);
    str = ['X = p.gas_' num2str(gid) ';'];
    eval(str);
    fprintf(fid,'%8.6e ',X(ll));
  end
  fprintf(fid,' \n');
end

fclose(fid);

