addpath /home/sergio/MATLABCODE
addpath /asl/matlib/h4tools

ctype        = [101 201];
water_cpsize = 11:3:29;
ice_cpsize   = 20:15:120;
cngwat       = [0 1e-2 1e-1 1 10 30 100 300];
cprtop       = 150:100:950;
scanang      = 0:10:50;

if length(water_cpsize) ~= length(ice_cpsize)
  error('water_cpsize ~= ice_cpsize');
end

number = length(ctype)*length(water_cpsize)*length(cngwat)*length(cprtop)*length(scanang);

[hall,ha,pall,pa] = rtpread('/home/sergio/KCARTA/IP_PROFILES/junk49_400ppm.op.rtp');

iRTP = input('enter replicate which of 49 profs? ');
[hnew,pnew] = replicate_rtp_headprof(hall,pall,iRTP,number);

ii = 0;
for cc = 1 : length(ctype)
  for dd = 1 : length(water_cpsize)
    for qq = 1 : length(cngwat)
      for tt = 1 : length(cprtop)
        for ss = 1 : length(scanang)
          ii = ii + 1;
          pnew.ctype(ii) = ctype(cc);
          if cc == 1
            pnew.cpsize(ii) = water_cpsize(dd);
          else
            pnew.cpsize(ii) = ice_cpsize(dd);
          end
          pnew.cngwat(ii) = cngwat(qq);
          pnew.cprtop(ii) = cprtop(tt);
          pnew.scanang(ii) = scanang(ss);
        end
      end
    end
  end
end

pnew.cfrac = ones(size(pnew.cfrac));
pnew.satzen  = vaconv(pnew.scanang,pnew.zobs,zeros(size(pnew.zobs)))
pnew.cprbot = pnew.cprtop + 50;

foutOP = ['parametrize_pclsam_disort_AFGL_' num2str(iRTP) '.op.rtp'];
foutRP = ['parametrize_pclsam_disort_AFGL_' num2str(iRTP) '.rp.rtp'];
rtpwrite(foutOP,hnew,ha,pnew,pa)

sartaer = ['!time /home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3 fin=' foutOP ' fout=' foutRP];
eval(sartaer);

[hnew,ha,pnew,pa] = rtpread(foutRP);
tnew = rad2bt(hnew.vchan,pnew.rcalc);
plot(pnew.stemp-tnew(1291,:)); title('stemp - BT1231 calc');
