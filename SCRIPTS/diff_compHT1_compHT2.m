addpath  /home/sergio/MATLABCODE/PLOTMISC
addpath  /home/sergio/SPECTRA
addpath  /home/sergio/SPECTRA/READ_XSEC

%%%%% if you have 2 monitors
%{
set(0, 'DefaultFigureRendererMode', 'manual')
set(0,'DefaultFigureRenderer','zbuffer')
%}
%%%%% if you have 2 monitors

%% function haha = shade2(fig,X0,Y0,W,H,color,transperancy)

file1 = '/asl/data/kcarta/KCARTADATA/General/compHT2010.param';

file2 = '/asl/data/kcarta/KCARTADATA/General/compHT2012.param';
file2 = '/home/sergio/KCARTA/SCRIPTS/MAKE_COMP_HTXY_PARAM_SC/PARAM_TEMP/testH2012.param';

data1 = load(file1);
data2 = load(file2);

H0 = 0.25;
H0 = 0.10;
H0 = 0.05;

iG = input('Enter gas list : (-1) for all, [A B] for selection, (+N) for specific gas : ');
if iG(1) < 0
  glist = 1:81;
elseif iG(1) > 0 & length(iG) == 1
  glist = iG; 
elseif iG(1) > 0 & length(iG) > 1
  glist = iG(1) : iG(2); 
end

for ggx = 1 : length(glist)
  gg = glist(ggx);

  if gg < 51
    lines = show_vis_ir_lines_wavenumber(2012,7,gg);
    %% map ln(min) --> 0, ln(max) --> 1
    damin = log(min(lines.stren));
    damax = log(max(lines.stren));
    map = (log(lines.stren)-damin)/(damax-damin);
    % figure(2); clf; semilogy(lines.wnum,lines.stren*1e20);
    % figure(3); clf; plot(lines.wnum,map);
  end

  figure(1); clf
  woo1 = find(data1(:,1) == gg);
  bwop1 = data1(woo1,:);

  woo2 = find(data2(:,1) == gg);
  bwop2 = data2(woo2,:);

  for ii = 1 : length(woo1)
    X0 = bwop1(ii,2);
    X1 = bwop1(ii,3);
    W  = X1-X0;

    Y0 = -H0 - H0/2;
    H  = H0;
    haha = shade2(1,X0,Y0,W,H,'b',0.75);

    %%%
    X0 = bwop1(ii,2);
    X1 = bwop1(ii,3);
    W  = X1-X0;

    Y0 = +2*H0 - H0/2;
    H  = H0;
    haha = shade2(1,X0,Y0,W,H,'b',0.10);

  end

  for ii = 1 : length(woo2)
    X0 = bwop2(ii,2);
    X1 = bwop2(ii,3);
    W  = X1-X0;

    Y0 = +H0 - H0/2;
    H  = H0;
    haha = shade2(1,X0,Y0,W,H,'r',0.75);

    %%%  
    X0 = bwop2(ii,2);
    X1 = bwop2(ii,3);
    W  = X1-X0;

    Y0 = -2*H0 - H0/2;
    H  = H0;
    haha = shade2(1,X0,Y0,W,H,'r',0.10);

  end

  grid on
  title(['Gas ' num2str(gg) ' blue = old, red = new']);

  if gg < 51
    figure(1); hold on; plot(lines.wnum,(map-0.5)*0.04,'k.','markersize',0.25); hold off
  else
    figure(1); hold on; plot_bands(gg); hold off
  end

  if ggx < length(glist)
    disp('<RET> to continue'); pause
  end

end

