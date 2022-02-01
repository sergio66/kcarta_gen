%% /home/sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/write_cloudfrac_profileMMM_100.m
%% see do_kcarta.m
%%      sed_1_2_cloudfiles.m

addpath /asl/matlib/h4tools

outcloud100 = ['cloudfrac_profile' num2str(iiBin) '_100'];

pjunk = hdfread(use_this_rtp0,'profiles','Fields',{'nlevs','cfrac','cfrac2','cfrac12','ctype','ctype2','plevs','ptemp','gas_201','gas_202','cprtop','cprtop2','cprbot','cprbot2'});
  pop.nlevs = pjunk{1};
  pop.cfrac = pjunk{2};
  pop.cfrac2 = pjunk{3};
  pop.cfrac12 = pjunk{4};
  pop.ctype  = pjunk{5};
  pop.ctype2 = pjunk{6};
  pop.plevs = pjunk{7};  
  pop.ptemp = pjunk{8};
  pop.gas_201 = pjunk{9};
  pop.gas_202 = pjunk{10};  
  pop.cprtop = pjunk{11};
  pop.cprtop2 = pjunk{12};
  pop.cprbot = pjunk{13};
  pop.cprbot2 = pjunk{14};
clear pjunk

pjunk = hdfread(use_this_rtp_ip,'profiles','Fields',{'tcc','cc','ciwc','clwc','plevs'});
  pip.tcc   = pjunk{1};
  pip.cc    = pjunk{2};
  pip.ciwc  = pjunk{3};
  pip.clwc  = pjunk{4};  
  pip.plevs = pjunk{5};
clear pjunk

%[hy,hay,pop,pay] = rtpread(use_this_rtp0);
%[hx,hax,pip,pax] = rtpread(use_this_rtp_ip);
fid = fopen(outcloud100,'w');
  str = ['% p1.cfrac + p1.cfrac2 - p1.cfrac12 = p1.tcc'];
  fprintf(fid,'%s \n',str);
  str = ['% iN/iY nlays=nlevs-1    tcc   cfrac   cfrac2  cfrac12  junk1'];
  fprintf(fid,'%s \n',str);
  [length(pop.nlevs) iiBin]
  
  fprintf(fid,'%2i  %2i %8.6f %8.6f %8.6f %8.6f %8.6f \n',iSigmaIASI,pop.nlevs(iiBin)-1,pip.tcc(iiBin),pop.cfrac(iiBin),pop.cfrac2(iiBin),pop.cfrac12(iiBin),0.0);

  playsN = pop.plevs(1:100,:) - pop.plevs(2:101,:);
  playsD = log(pop.plevs(1:100,:) ./ pop.plevs(2:101,:));
  plays = playsN ./ playsD;
  plays(101,:) = plays(100,:)+1;
  
  str = ['% individual layers : laynum  cc      gas_201(W)  gas_202(I) gas_203(A) ptemp plevs'];
  str = ['% individual layers : laynum  plays   cc  gas_201(W)  gas_202(I) gas_203(A)  ptemp'];  
  fprintf(fid,'%s \n',str);

  cc = interp1(log10(pip.plevs(:,iiBin)),pip.cc(:,iiBin),log10(plays(:,iiBin)),[],'extrap');
  cc = cc(1:100);
%  matrix = [(100:-1:1)' cc  pop.gas_201(1:100,iiBin) pop.gas_202(1:100,iiBin) pop.gas_201(1:100,iiBin)*0 pop.ptemp(1:100,iiBin) plays(1:100,iiBin)];
  matrix = [(100:-1:1)' plays(1:100,iiBin) cc  pop.gas_201(1:100,iiBin) pop.gas_202(1:100,iiBin) pop.gas_201(1:100,iiBin)*0 pop.ptemp(1:100,iiBin)];
  if plays(50) < plays(51)
    %% recall kCARTA has p(1) = 1100 mb, p(2) = 1070 mb ... p(101) = 0.005 mb so want the pressures DECREASING from GND to TOA
    matrix = flipud(matrix);
  end
  matrix(isnan(matrix)) = 0;
  fprintf(fid,'%3i  %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f \n',matrix');

fclose(fid);

[pop.ctype(iiBin) pop.ctype2(iiBin)]
[pip.tcc(iiBin)  pop.cfrac(iiBin)  pop.cfrac2(iiBin)  pop.cfrac12(iiBin) ]

pop.plays = plays;
figure(1); clf
  plot(pip.cc(:,iiBin),pip.plevs(:,iiBin),'kd-',cc,pop.plays(1:100,iiBin),'gd-',...
       pip.clwc(:,iiBin)*1e5,pip.plevs(:,iiBin),'bs-',pop.gas_201(:,iiBin),pop.plays(:,iiBin),'rd-','linewidth',2)
  ax = axis; line([ax(1) ax(2)],[pop.cprtop2(iiBin) pop.cprtop2(iiBin)],'color','k'); line([ax(1) ax(2)],[pop.cprbot2(iiBin) pop.cprbot2(iiBin)],'color','k'); axis(ax);        
  title('CTYPE = 101 (H2O) (k) cc levels (g) cc lays (b) clwc levels (r) gas201 lays')
  set(gca,'ydir','reverse')
figure(2); clf
  plot(pip.cc(:,iiBin),pip.plevs(:,iiBin),'kd-',cc,pop.plays(1:100,iiBin),'gd-',...
       pip.ciwc(:,iiBin)*1e5,pip.plevs(:,iiBin),'bs-',pop.gas_202(:,iiBin),pop.plays(:,iiBin),'rd-','linewidth',2)
  ax = axis; line([ax(1) ax(2)],[pop.cprtop(iiBin) pop.cprtop(iiBin)],'color','k'); line([ax(1) ax(2)],[pop.cprbot(iiBin) pop.cprbot(iiBin)],'color','k'); axis(ax);        
  title('CTYPE = 201 (ICE) (k) cc levels (g) cc lays (b) ciwc levels (r) gas202 lays')
  set(gca,'ydir','reverse')

figure(1); set(gca,'fontsize',10)
figure(2); set(gca,'fontsize',10)