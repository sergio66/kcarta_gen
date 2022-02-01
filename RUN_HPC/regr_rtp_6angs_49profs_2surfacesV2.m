addpath /asl/matlib/h4tools
addpath /asl/matlib/aslutil
addpath /asl/matlib/rtptools
addpath /home/sergio/MATLABCODE

inang = [0 32.8244 44.8285 53.4704 59.8336 65.0428];
inang = [0 32.8244 44.8285 53.4704 59.8336 60.1100];
inang = [0 32.8244 44.8285 53.4704 55.0000 57.7300];

rSatHgt = 705; % AIRS
rSatHgt = 841; % CrIS

%% go from satzen (at gnd) to scanang (at satellite) USE THIS
[inang2scanang] = saconv(inang , [rSatHgt rSatHgt rSatHgt rSatHgt rSatHgt rSatHgt]*1000)

%% go from scanang (at satellite) to satzen (at gnd)
[inang2satzen] = vaconv(inang , [rSatHgt rSatHgt rSatHgt rSatHgt rSatHgt rSatHgt]*1000, [0 0 0 0 0 0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1 = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1013_400ppm.op.rtp';
file2 = '/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr49_1100_400ppm.op.rtp';

[h,ha,p1,pa] = rtpread(file1);
[h,ha,p2,pa] = rtpread(file2);

scanang = zeros(1,6);
satzen  = zeros(1,6);
for ii = 1 : 6
  junk = p1;
  if ii == 1
    p1out = junk;
  else
    junk.scanang = ones(size(junk.stemp)) * inang2scanang(ii);
    junk.satzen = ones(size(junk.stemp)) * inang(ii);
    [h,p1out] = cat_rtp(h,p1out,h,junk);
  end
end
p1out.zobs = ones(size(p1out.zobs)) * rSatHgt * 1000;
rtpwrite('/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_1013mb_seaemis.rtp',h,ha,p1out,pa);
fprintf(1,'spres, emis  %8.6f %8.6f \n',[p1out.spres(1) p1out.emis(1,1)])

p1out.emis = ones(size(p1out.emis));
p1out.rho  = zeros(size(p1out.emis));  
rtpwrite('/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_1013mb_unitemis.rtp',h,ha,p1out,pa);
fprintf(1,'spres, emis  %8.6f %8.6f \n',[p1out.spres(1) p1out.emis(1,1)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1 : 6
  junk = p2;
  if ii == 1
    p2out = junk;
  else
    junk.scanang = ones(size(junk.stemp)) * inang2scanang(ii);
    junk.satzen = ones(size(junk.stemp)) * inang(ii);
    [h,p2out] = cat_rtp(h,p2out,h,junk);
  end
end
p2out.zobs = ones(size(p2out.zobs)) * rSatHgt * 1000;
p2out.spres = 1080 * ones(size(p2out.spres));
rtpwrite('/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_1080mb_seaemis.rtp',h,ha,p2out,pa);
fprintf(1,'spres, emis  %8.6f %8.6f \n',[p2out.spres(1) p2out.emis(1,1)])

p2out.emis = ones(size(p2out.emis));
p2out.rho  = zeros(size(p2out.emis));  
rtpwrite('/home/sergio/MATLABCODE/REGR_PROFILES/RUN_KCARTA/REGR49_400ppm/regr_rtp_6angs_49profs_1080mb_unitemis.rtp',h,ha,p2out,pa);
fprintf(1,'spres, emis  %8.6f %8.6f \n',[p2out.spres(1) p2out.emis(1,1)])

