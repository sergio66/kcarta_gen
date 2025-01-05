addpath /asl/matlib/h4tools
addpath /home/sergio/MATLABCODE/CONVERT_GAS_UNITS

[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op.rtp');
[h,ha,p,pa] = rtpread('/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/latbin1_40.op_400ppm.rtp');

ppmv2 = layers2ppmv(h,p,1:40,2); %% 400 ppm
ppmv4 = layers2ppmv(h,p,1:40,4); %% 0.32 ppm
ppmv6 = layers2ppmv(h,p,1:40,6); %% 1.8 ppm

%% plot(x0.hi_fcris,x.hi_rcris_all(:,1))  CO2
%% plot(x0.hi_fcris,x.hi_rcris_all(:,2))  N2O
%% plot(x0.hi_fcris,x.hi_rcris_all(:,3))  CO
%% plot(x0.hi_fcris,x.hi_rcris_all(:,4))  CH4
%% plot(x0.hi_fcris,x.hi_rcris_all(:,5))  CFC11
%% plot(x0.hi_fcris,x.hi_rcris_all(:,6))  CFC12
%% plot(x0.hi_fcris,x.hi_rcris_all(:,7))  Tz
%% plot(x0.hi_fcris,x.hi_rcris_all(:,8))  Stemp

iX = 1;
iX = 2;
for ii = 1 : 40
  if iX == 1
    fjac = ['New_hamm_ng=2/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_coljac.mat'];
  else
    fjac = ['New_sinc_ng=2/individual_prof_convolved_kcarta_crisHI_crisMED_' num2str(ii) '_coljac.mat'];
  end
  x = load(fjac);
  jacHI(ii,:,:) = x.hi_rcris_all(:,[1 2 4 8]);
  jacLO(ii,:,:) = x.lo_rcris_all(:,[1 2 4 8]);
  jacMED(ii,:,:) = x.med_rcris_all(:,[1 2 4 8]);
end

jacHI(:,:,1) = jacHI(:,:,1) * 2.2/400;
jacHI(:,:,2) = jacHI(:,:,2) * 1/320;
jacHI(:,:,3) = jacHI(:,:,3) * 5/1800;

jacMED(:,:,1) = jacMED(:,:,1) * 2.2/400;
jacMED(:,:,2) = jacMED(:,:,2) * 1/320;
jacMED(:,:,3) = jacMED(:,:,3) * 5/1800;

jacLO(:,:,1) = jacLO(:,:,1) * 2.2/400;
jacLO(:,:,2) = jacLO(:,:,2) * 1/320;
jacLO(:,:,3) = jacLO(:,:,3) * 5/1800;

fLO = x.lo_fcris;
fHI = x.hi_fcris;
fMED = x.med_fcris;

clear fjac
comment = 'normalization 2.2/400 for CO2, 1/320 for N2O and 5/1800 for CH4; 1.0 for stemp';
if iX == 1
  save hamm_ng2_coljac.mat jac* f* comment
elseif iX == 2
  save sinc_ng2_coljac.mat jac* f* comment
end
