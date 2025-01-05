addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE

sarta1 = '/asl/packages/sartaV108/BinV201/sarta_apr08_m140_wcon_nte';
sarta2 = '/home/chepplew/gitLib/sarta/bin/airs_l1c_2834_cloudy_may19';

junk = load('individual_prof_convolved_kcarta_crisHI_crisMED_16.mat');

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op.rtp');
h.ichan = (1:2834)';
h.nchan = 2834;
h.vchan = junk.fKc;
p.robs1 = zeros(2834,40);
p.rcalc = zeros(2834,40);

rtpwrite('sartajac.op.rtp',h,ha,p,pa);
eval(['!' sarta1 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pjunk,pax] = rtpread('sartajac.rp.rtp');
rads1_40 = pjunk.rcalc;

rtpwrite('sartajac.op.rtp',h,ha,p,pa);
eval(['!' sarta2 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pjunk,pax] = rtpread('sartajac.rp.rtp');
rads2_40 = pjunk.rcalc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnew] = replicate_rtp_headprof(h,p,20,98);
for ii = 1 :97
  pnew.gas_1(ii,ii) = pnew.gas_1(ii,ii) * 1.1;
end
rtpwrite('sartajac.op.rtp',hnew,ha,pnew,pa);

eval(['!' sarta1 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxWV1,pax] = rtpread('sartajac.rp.rtp');
eval(['!' sarta2 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxWV2,pax] = rtpread('sartajac.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnew] = replicate_rtp_headprof(h,p,20,98);
for ii = 1 :97
  pnew.gas_3(ii,ii) = pnew.gas_3(ii,ii) * 1.1;
end
rtpwrite('sartajac.op.rtp',hnew,ha,pnew,pa);

eval(['!' sarta1 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxO31,pax] = rtpread('sartajac.rp.rtp');
eval(['!' sarta2 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxO32,pax] = rtpread('sartajac.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[hnew,pnew] = replicate_rtp_headprof(h,p,20,99);
for ii = 1 :97
  pnew.ptemp(ii,ii) = pnew.ptemp(ii,ii) + 1;
end
pnew.stemp(98) = pnew.stemp(98) + 1;
rtpwrite('sartajac.op.rtp',hnew,ha,pnew,pa);

eval(['!' sarta1 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxTz1,pax] = rtpread('sartajac.rp.rtp');
eval(['!' sarta2 ' fin=sartajac.op.rtp fout=sartajac.rp.rtp']);
[hx,hax,pxTz2,pax] = rtpread('sartajac.rp.rtp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tjunk = rad2bt(hx.vchan,pxTz1.rcalc);
Tz1jac = tjunk(:,1:97) - tjunk(:,99)*ones(1,97);
ST1jac = tjunk(:,98) - tjunk(:,99);

tjunk = rad2bt(hx.vchan,pxTz2.rcalc);
Tz2jac = tjunk(:,1:97) - tjunk(:,99)*ones(1,97);
ST2jac = tjunk(:,98) - tjunk(:,99);

tjunk = rad2bt(hx.vchan,pxWV1.rcalc);
WV1jac = tjunk(:,1:97) - tjunk(:,98)*ones(1,97); WV1jac = WV1jac/log(1+0.1);

tjunk = rad2bt(hx.vchan,pxWV2.rcalc);
WV2jac = tjunk(:,1:97) - tjunk(:,98)*ones(1,97); WV2jac = WV2jac/log(1+0.1);

tjunk = rad2bt(hx.vchan,pxO31.rcalc);
rads1 = tjunk(:,98);
O31jac = tjunk(:,1:97) - tjunk(:,98)*ones(1,97); O31jac = O31jac/log(1+0.1);

tjunk = rad2bt(hx.vchan,pxO32.rcalc);
rads2 = tjunk(:,98);
O32jac = tjunk(:,1:97) - tjunk(:,98)*ones(1,97); O32jac = O32jac/log(1+0.1);

freqs2834 = hx.vchan;
save sarta_finitediff_jacs_prof20.mat freqs2834 rads* *jac sarta1 sarta2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); clf
  subplot(211); plot(freqs2834,ST1jac); title('STEMP'); 
  subplot(212); plot(freqs2834,ST1jac-ST2jac);

figure(2); clf
  subplot(211); plot(freqs2834,sum(Tz1jac')); title('Tz'); 
  subplot(212); plot(freqs2834,sum(Tz1jac'-Tz2jac'));

figure(3); clf
  subplot(211); plot(freqs2834,sum(WV1jac')); title('WV'); 
  subplot(212); plot(freqs2834,sum(WV1jac'-WV2jac'));

figure(4); clf
  subplot(211); plot(freqs2834,sum(O31jac')); title('O3'); 
  subplot(212); plot(freqs2834,sum(O31jac'-O32jac'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('ret to see percent change'); pause

figure(1); clf
  subplot(211); plot(freqs2834,ST1jac); title('STEMP'); 
  subplot(212); plot(freqs2834,100*(ST1jac-ST2jac)./ST1jac);
    ylim([-20 +20])

figure(2); clf
  subplot(211); plot(freqs2834,sum(Tz1jac')); title('Tz'); 
  subplot(212); plot(freqs2834,100*sum(Tz1jac'-Tz2jac')./sum(Tz1jac'));
    ylim([-20 +20])

figure(3); clf
  subplot(211); plot(freqs2834,sum(WV1jac')); title('WV'); 
  subplot(212); plot(freqs2834,100*sum(WV1jac'-WV2jac')./sum(WV1jac'));
    ylim([-20 +20])

figure(4); clf
  subplot(211); plot(freqs2834,sum(O31jac')); title('O3'); 
  subplot(212); plot(freqs2834,100*sum(O31jac'-O32jac')./sum(O31jac'));
    ylim([-20 +20])
