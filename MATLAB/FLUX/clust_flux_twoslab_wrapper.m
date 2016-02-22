% local running to test
% clustcmd -L clust_flux_twoslab_wrapper.m 001:12150
%
% otherwise when happy
% clustcmd -q medium -n 128 clust_flux_twoslab_wrapper.m 001:12150
%
% or
% clustcmd -q long_contrib -n 128 clust_flux_twoslab_wrapper.m 001:12150

addpath /asl/matlib/h4tools
addpath /asl/matlib/rtptools
addpath /asl/matlib/aslutil

iRTP = 7609;
iRTP = JOB;

fnameIN = '/home/sergio/MATLABCODE/RTPMAKE/MFILES_2011_03_11/junk039.op.rtp';

outname = ['JUNK/RESULTS/flux_' num2str(iRTP,'%05d') '.mat'];

ee = exist(outname);
if ee == 0
  [hIN,haIN,pIN,paIN] = rtpread(fnameIN);
  [flux,fluxdata] = driver_flux_twoslab_wrapper(hIN,pIN,iRTP);
  rlat = pIN.rlat(iRTP); 
  rlon = pIN.rlon(iRTP);
  solzen = pIN.solzen(iRTP);
  saver = ['save ' outname ' flux fluxdata rlat rlon solzen'];
  eval(saver)
else
  fprintf(1,'%s already exists \n',outname);
end
