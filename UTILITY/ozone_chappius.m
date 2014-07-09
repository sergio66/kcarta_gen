%%% put gas2 to trick things

f1 = 4500;
f2 = 22000;
nbox = 5; 
pointsPerChunk = 10000; 

fend = f2;
f0 = f1;
topts = runXtopts_params_smart(f0);  
dv = topts.ffin*nbox*pointsPerChunk; 

gases = [1 2 3 4 5 6 7];

gasid = 3;

dvx = -1;
f0 = f1;
while f0 < fend
  topts = runXtopts_params_smart(f0); 
  dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = f0 + dv; 
  if abs(dv-dvx) > 1
    disp('----------------------------------------------------------')
    end
  dvx = dv;
  fprintf(1,'%10.2f  %10.6f  %10.6f  %10.6f\n',f0,fmax,topts.ffin,dv);
  f0 = f0 + dv;
  end

wall  = [];
dallA = [];
dallC = [];
dallX = [];

cd /home/sergio/SPECTRA
f0 = f1;
while f0 < fend
  topts = runXtopts_params_smart(f0); 
  dv = topts.ffin*nbox*pointsPerChunk; 
  fmax = f0 + dv; 

  fnameIN = ['/carrot/s1/sergio/RUN8_VISDATABASE/trp' num2str(f0) '.mat'];

  if (exist(fnameIN))
    loader = ['load ' fnameIN];
    eval([loader]);

    clear dozone

    fip = ['IPFILES/trp_g' num2str(gasid)];
    qamt = load(fip);
    qamt = qamt(:,5)*6.023e26;
    figure(1)
    for ll = 1 : 100
      dozone(:,ll) = chappius_ozone_vis(wvis0,qamt(ll));
      end
    dozone = dozone';

    fr = wvis0;
    k00 = dvis0 + dozone;

    figure(2);
    wall = [wall fr];
    dallA = [dallA sum(dvis0)];
    dallC = [dallC sum(dozone)];
    dallX = [dallX sum(k00)];

    plot(10000./wall,exp(-dallA),10000./wall,exp(-dallC),...
         10000./wall,exp(-dallX),'r'); pause(0.1)

  else
    fprintf(1,'%s NOT FOUND \n',fnameIN);
    end

  f0 = f0 + dv;
  end

modis_channels = [470 550 659 865 1240 1640 2130]*1e-3;
modis_channels = 10000./modis_channels;
figure(1);clf
plot(wall,exp(-dallA),wall,exp(-dallC),...
     wall,exp(-dallX),'r',...
     modis_channels,ones(size(modis_channels)),'kx'); 

figure(2);
plot(10000./wall,exp(-dallA),10000./wall,exp(-dallC),...
     10000./wall,exp(-dallX),'r',...
     10000./modis_channels,ones(size(modis_channels)),'kx'); 

figure(3); plot(10000./wall,exp(-dallC))
axis([0.5 1.0 0 1]); grid