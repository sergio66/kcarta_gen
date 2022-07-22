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

%iRTP = input('enter replicate which of 49 profs? ');
iRTP = 1;
[hnew,pnew1] = replicate_rtp_headprof(hall,pall,iRTP,number);

ii = 0;
for cc = 1 : length(ctype)
  for dd = 1 : length(water_cpsize)
    for qq = 1 : length(cngwat)
      for tt = 1 : length(cprtop)
        for ss = 1 : length(scanang)
          ii = ii + 1;
          pnew1.ctype(ii) = ctype(cc);
          if cc == 1
            pnew1.cpsize(ii) = water_cpsize(dd);
          else
            pnew1.cpsize(ii) = ice_cpsize(dd);
          end
          pnew1.cngwat(ii) = cngwat(qq);
          pnew1.cprtop(ii) = cprtop(tt);
          pnew1.scanang(ii) = scanang(ss);
        end
      end
    end
  end
end

pnew1.cfrac = ones(size(pnew1.cfrac));
pnew1.satzen  = vaconv(pnew1.scanang,pnew1.zobs,zeros(size(pnew1.zobs)));
pnew1.cprbot = pnew1.cprtop + 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now combine ICE n WATER

number = length(ctype)*length(water_cpsize)*length(cngwat)*length(cprtop)*length(scanang);
number = length(water_cpsize)*length(cngwat)*length(cprtop)*length(scanang);

iRTP = 1;
[hnew,pnew2] = replicate_rtp_headprof(hall,pall,iRTP,number);

ii = 0;
for dd = 1 : length(water_cpsize)
  for qq = 1 : length(cngwat)
    for tt = 1 : length(cprtop)
      for ss = 1 : length(scanang)
        top   = cprtop(tt);
        top2  = cprtop(length(cprtop)-tt+1);
        bot   = top + 50;
        bot2  = top2 + 50;

        ii = ii + 1;
        pnew2.ctype(ii) = 101;
        pnew2.ctype2(ii) = 201;
        pnew2.cpsize(ii)  = water_cpsize(dd);
        pnew2.cpsize2(ii) = ice_cpsize(length(water_cpsize)-dd+1);
        pnew2.cngwat(ii)  = cngwat(qq);
        pnew2.cngwat2(ii) = cngwat(length(cngwat)-qq+1);
        pnew2.scanang(ii) = scanang(ss);

        if top == top2
          pnew2.cprtop(ii)  = top - 25;
          pnew2.cprtop2(ii) = top + 25;
        elseif bot < top2
          pnew2.cprtop(ii)  = top;
          pnew2.cprtop2(ii) = top2;
        elseif bot > top2
          pnew2.cprtop(ii)  = top2;
          pnew2.cprtop2(ii) = top;
        end
      end
    end
  end
end

pnew2.cfrac   = ones(size(pnew2.cfrac));
pnew2.cfrac2  = ones(size(pnew2.cfrac));
pnew2.cfrac12 = ones(size(pnew2.cfrac));
pnew2.satzen  = vaconv(pnew2.scanang,pnew2.zobs,zeros(size(pnew2.zobs)));
pnew2.cprbot  = pnew2.cprtop + 50;
pnew2.cprbot2 = pnew2.cprtop2 + 50;

bad = find(pnew2.cprbot > pnew2.cprtop2);
if length(bad) > 0
  whos bad
  error('ksjg;kjggs bad')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /asl/matlib/rtptools
[hnew,pnew3] = cat_rtp(hnew,pnew1,hnew,pnew2);

%% now put in the profiles 1 -- 49
iMod = length(pnew3.stemp)/49;

pnewF = pnew3;
for ii = 1 : 49
  ind = (1 : 49 : length(pnew3.stemp)+100) + (ii-1);
  boo = find(ind <= length(pnew3.stemp));
  ind = ind(boo);
  pnewF.stemp(ind) = pall.stemp(ii);
  pnewF.plat(ind)  = pall.plat(ii);
  pnewF.plon(ind)  = pall.plon(ii);
  pnewF.plevs(:,ind)  = pall.plevs(:,ii) * ones(1,length(ind));
  pnewF.palts(:,ind)  = pall.palts(:,ii) * ones(1,length(ind));
  pnewF.ptemp(:,ind)  = pall.ptemp(:,ii) * ones(1,length(ind));
  pnewF.gas_1(:,ind)  = pall.gas_1(:,ii) * ones(1,length(ind));
  pnewF.gas_2(:,ind)  = pall.gas_2(:,ii) * ones(1,length(ind));
  pnewF.gas_3(:,ind)  = pall.gas_3(:,ii) * ones(1,length(ind));
  pnewF.gas_4(:,ind)  = pall.gas_4(:,ii) * ones(1,length(ind));
  pnewF.gas_5(:,ind)  = pall.gas_5(:,ii) * ones(1,length(ind));
  pnewF.gas_6(:,ind)  = pall.gas_6(:,ii) * ones(1,length(ind));
  pnewF.gas_9(:,ind)  = pall.gas_9(:,ii) * ones(1,length(ind));
  pnewF.gas_12(:,ind) = pall.gas_12(:,ii) * ones(1,length(ind));
end
mmw = mmwater_rtp(hnew,pnewF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

foutOP = ['parametrize_pclsam_disort_AFGL_1_49.op.rtp'];
foutRP = ['parametrize_pclsam_disort_AFGL_1_49.rp.rtp'];
rtpwrite(foutOP,hnew,ha,pnewF,pa)

sartaer = ['!time /home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19_prod_v3 fin=' foutOP ' fout=' foutRP];
eval(sartaer);

[hnew,ha,pnew,pa] = rtpread(foutRP);
tnew = rad2bt(hnew.vchan,pnew.rcalc);
plot(pnew.stemp-tnew(1291,:)); title('stemp - BT1231 calc');
