%% oops had factor of 1000 too low gas amounts in orig OD files, so need to fix

error('already fixed!!')

addpath /home/sergio/MATLABCODE/
addpath /asl/matlib/aslutil/

figure(1); clf;
for wn = 605 : 25 : 830
  outname = ['mars_co2_ods_' num2str(wn) '.mat'];
  fprintf(1,'doing %3i chunk \n',wn)
  loader = ['load ' outname];
  eval(loader);
  data(:,5) = data(:,5) * 1000;
  out_array = out_array * 1000;
  saver = ['save ' outname ' data out_array outwave'];
  eval(saver);
end
