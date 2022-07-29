addpath /home/sergio/KCARTA/MATLAB

[~,~,ice,~] = rtpread('/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_icecld_USSTD_7cngwat.op.rtp');
[~,~,water,~] = rtpread('/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/rrtm_lw_TAMU_v0/run_TangJAS2018/regr49_watercld_USSTD_7cngwat.op.rtp');

for ii = 1 :7
  fin = ['Prof1/GUIDO/FarIR/DISORT/ICE/rad.dat' num2str(ii) '_0310'];
  [disort_ice(:,ii),w] = readkcstd(fin);

  fin = ['Prof1/GUIDO/FarIR/DISORT/WATER/rad.dat' num2str(ii) '_0310'];
  [disort_water(:,ii),w] = readkcstd(fin);

  fin = ['Prof1/GUIDO/FarIR/CHOU/ICE/rad.dat' num2str(ii)];
  [junk,w] = readkcstd(fin);
  [mm,nn] = size(junk);
  chou_ice(:,ii) = junk(:,nn);

  fin = ['Prof1/GUIDO/FarIR/CHOU/WATER/rad.dat' num2str(ii)];
  [junk,w] = readkcstd(fin);
  [mm,nn] = size(junk);
  chou_water(:,ii) = junk(:,nn);

  fin = ['Prof1/GUIDO/FarIR/TANG/ICE/rad.dat' num2str(ii)];
  [junk,w] = readkcstd(fin);
  [mm,nn] = size(junk);
  tang_ice(:,ii) = junk(:,nn);

  fin = ['Prof1/GUIDO/FarIR/TANG/WATER/rad.dat' num2str(ii)];
  [junk,w] = readkcstd(fin);
  [mm,nn] = size(junk);
  tang_water(:,ii) = junk(:,nn);
end

[fc,qcDI] = quickconvolve(w,disort_ice,1,1);   tcDI = rad2bt(fc,qcDI);
[fc,qcDW] = quickconvolve(w,disort_water,1,1); tcDW = rad2bt(fc,qcDW);
[fc,qcCI] = quickconvolve(w,chou_ice,1,1);     tcCI = rad2bt(fc,qcCI);
[fc,qcCW] = quickconvolve(w,chou_water,1,1);   tcCW = rad2bt(fc,qcCW);
[fc,qcTI] = quickconvolve(w,tang_ice,1,1);     tcTI = rad2bt(fc,qcTI);
[fc,qcTW] = quickconvolve(w,tang_water,1,1);   tcTW = rad2bt(fc,qcTW);

plot(fc,nanmean(tcDI'),fc,nanmean(tcDW'))
plot(fc,nanmean(tcDI'-tcCI'),'b',fc,nanmean(tcDW'-tcCW'),'r',fc,nanmean(tcDI'-tcTI'),'c',fc,nanmean(tcDW'-tcTW'),'m'); 
hl = legend('Chou ICE','Chou WATER','Tang ICE','Tang WATER','location','best'); ylabel('DISORT-CHOU BT(K)');


figure(1)
  plot(fc,tcDI,'linewidth',2); hl = legend(num2str(ice.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT BT(K)'); title('Ice Cloud 200 mb, 50 um deff')
figure(2)
  plot(fc,tcDW,'linewidth',2); hl = legend(num2str(water.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT BT(K)'); title('Water Cloud 800 mb, 20 um deff')

figure(1)
  plot(fc,tcDI-tcCI,'linewidth',2); hl = legend(num2str(ice.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT-CHOU BT(K)'); title('Ice Cloud 200 mb, 50 um deff')
figure(2)
  plot(fc,tcDW-tcCW,'linewidth',2); hl = legend(num2str(water.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT-CHOU BT(K)'); title('Water Cloud 800 mb, 20 um deff')

figure(1)
  plot(fc,tcDI-tcTI,'linewidth',2); hl = legend(num2str(ice.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT-TANG BT(K)'); title('Ice Cloud 200 mb, 50 um deff')
figure(2)
  plot(fc,tcDW-tcTW,'linewidth',2); hl = legend(num2str(water.cngwat'),'location','best','fontsize',10); xlabel('Wavenumber cm-1'); ylabel('DISORT-TANG BT(K)'); title('Water Cloud 800 mb, 20 um deff')
