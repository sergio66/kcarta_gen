%%% input hpcfname

%% need to modify template_Qrad.nml CORRECTLY for the rtp file to process!

f1 = 605; f2 = 2830; mm = -1;   %% dummies, needed for now
use_this_rtp = '/home/sergio/MATLABCODE/RANDOM_LARRABEE/bdry_layer_wv.rtp'

iDoJac = +1;   %% do Jacobians
iDoJac = -1;   %% do rads

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% everywhere you find '/' in use_this_rtp, replace it with '\/'

ooh = strfind(use_this_rtp,'/');
if length(ooh) > 0
  use_this_rtp = strrep(use_this_rtp, '/', '\/');
end

iiBin = JOB;
hpcfname
do_kcarta

