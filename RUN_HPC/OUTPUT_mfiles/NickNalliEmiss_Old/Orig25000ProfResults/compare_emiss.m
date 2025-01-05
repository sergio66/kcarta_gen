%% test nicknalli emis, symbolic link to
%% /home/sergio/MATLABCODE/BRDF_EMISSIVITY_NALLI/junk.rp.rtp

%% see /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/set_rtp.m

use_this_rtp0 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalli_masudaemis.rtp';    %% masuda emis
use_this_rtp1 = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis.rtp';           %% new emis
use_this_rtp1x = '/home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/RTP/testnalliemis_resetTWV.rtp';  %% new emis, reset ST, CO2,CH4, WV

if ~exist('p1x')
  [h,ha,p0,pa] = rtpread(use_this_rtp0);
  [h,ha,p1,pa] = rtpread(use_this_rtp1);
  [h,ha,p1x,pa] = rtpread(use_this_rtp1x);
end

disp('p0 = Masuda  p1 = Nalli p1x = Nalli+RESET')
fprintf(1,'f0 = %s \n',use_this_rtp0);
fprintf(1,'f1 = %s \n',use_this_rtp1);
fprintf(1,'f1x = %s \n',use_this_rtp1x);

fprintf(1,'stemp : p0-p1 p0-p1x p1-p1x %8.6e %8.6e %8.6e \n',sum(p0.stemp-p1.stemp),sum(p0.stemp-p1x.stemp),sum(p1.stemp-p1x.stemp))
fprintf(1,'wv    : p0-p1 p0-p1x p1-p1x %8.6e %8.6e %8.6e \n',sum(p0.gas_1(:)-p1.gas_1(:)),sum(p0.gas_1(:)-p1x.gas_1(:)),sum(p1.gas_1(:)-p1x.gas_1(:)))
fprintf(1,'efreq : p0-p1 p0-p1x p1-p1x %8.6e %8.6e %8.6e \n',sum(p0.efreq(:)-p1.efreq(:)),sum(p0.efreq(:)-p1x.efreq(:)),sum(p1.efreq(:)-p1x.efreq(:)))
fprintf(1,'emis  : p0-p1 p0-p1x p1-p1x %8.6e %8.6e %8.6e \n',nansum(p0.emis(:)-p1.emis(:)),nansum(p0.emis(:)-p1x.emis(:)),nansum(p1.emis(:)-p1x.emis(:)))
