xstartup
fin0  = '/asl/s1/sergio/pin_feb2002_sea_airsnadir_ip.so2.rtp';
fin  = 'junk.ip.rtp';
fout = 'junk.op.rtp';

cper = ['!/bin/cp ' fin0 ' ' fin]; eval(cper)
[h,ha,p,pa] = rtpread(fin);
p.spres(:) = 1100;
rtpwrite(fin,h,ha,p,pa);

for iProf = 2 : 6
  klayers = '/asl/packages/klayers/Bin/klayers_airs ';
  %% note we use nwant=-1 and mnafgl=iX  so as to use the correct filler (else default = 6 = USSTD)
  klayerser = ['!' klayers ' fin=' fin ' fout=' fout ' nwant=-1 mnafgl=' num2str(iProf-1)];
  eval(klayerser);

  [h,ha,p,pa] = rtpread(fout);
  [hx,px]     = subset_rtp(h,p,[],[],iProf-1);
  for iGas = 1:81
    fSTD = ['/asl/data/kcarta/KCARTADATA/RefProf_July2010.For.v115up_CO2ppmv385/refgas' num2str(iGas)];
    ee = exist(fSTD);
    if ee > 0
      fNEW = ['/asl/data/kcarta/KCARTADATA/RefProf_July2010.For.v115up_CO2ppmv385/'];
      fNEW = [fNEW num2str(iProf) '/afglgas' num2str(iGas)];
      cper = ['!/bin/cp ' fSTD ' ' fNEW];
      eval(cper);
      fprintf(1,'copied gas %3i from refgas \n',iGas);

      fieldx = ['gas_' num2str(iGas)];
      if isfield(p,fieldx)
        gaser = ['gas = px.' fieldx ';'];
        eval(gaser)
        gas = gas/6.023e26;
        fid = fopen(fNEW,'w');
        str = ['!AFGL Profile' num2str(iProf) ' made from /home/sergio/KCARTA/UTILITY/AFGLprofs.m'];
        fprintf(fid,'%s \n',str);
        str = ['! GasID ' num2str(iGas)];
        fprintf(fid,'%s \n',str);
        str = ['! lay    pressure    part pres    temp      amount'];
        fprintf(fid,'%s \n',str);
        str = ['!----  -----------  -----------  -------  -----------'];
        fprintf(fid,'%s \n',str);
        dz = [diff(px.palts); 200];   %% in meters
        num_moles = (px.plays*100 .* abs(dz))/8.31./px.ptemp;  %% this is moles per sq meter
        num_moles = num_moles/10000;                           %% this is moles per sq cm
        num_moles = num_moles/1000;                            %% this is kilomoles per sq cm
        mix_ratio = gas ./num_moles;
        array = [ones(101,1)*iGas px.plays/1013.25 px.plays/1013.25 .* mix_ratio  px.ptemp gas];
        array = array(1:100,:);
        array = flipud(array);
        fprintf(fid,'%3i %10.4e %10.4e %10.4f %10.4e \n',array');
        fclose(fid);
        fprintf(1,'  made gas %3i %s \n',iGas,fNEW)
%if iGas == 2
%  keyboard
%end

      end
    end
  end
end

rmer = ['!/bin/rm ' fin ' ' fout]; eval(rmer);

%% final profs look like
%{
!AFGL Profile2 made from /home/sergio/KCARTA/UTILITY/AFGLprofs.m 
! GasID 6 
! lay    pressure    part pres    temp      amount 
!----  -----------  -----------  -------  ----------- 
  6 1.0712e+00 1.8660e-06   303.3398 1.8081e-09 
  6 1.0427e+00 1.8190e-06   301.9182 1.7852e-09 
  6 1.0146e+00 1.7727e-06   300.4778 1.7620e-09 
  6 9.8687e-01 1.7269e-06   299.0181 1.7387e-09 
  6 9.5955e-01 1.6817e-06   297.5387 1.7152e-09 
  6 9.3264e-01 1.6371e-06   296.0393 1.6916e-09 
%}