g1515 = 'l2s_kc120_G15_xsecG15.dat';  %% GEISA 2015
f1616 = 'l2s_kc120_H16.dat';          %% molgas = H2016, xscgas = H2016
f1612 = 'l2s_kc120_H16_xsecH12.dat';  %% molgas = H2016, xscgas = H2012
f1212 = 'l2s_kc120_H12.dat';          %% molgas = H2012, xscgas = H2012

addpath /home/sergio/KCARTA/MATLAB

[d1515,w] = readkcstd(g1515);
[d1616,w] = readkcstd(f1616);
[d1612,w] = readkcstd(f1612);
[d1212,w] = readkcstd(f1212);

gaslist = load('/home/sergio/KCARTA/SCRIPTS/MAKE_COMP_HTXY_PARAM_SC/PARAM_HT2016/comp_ir605_2830.param');
gaslist = unique(gaslist(:,1));
gaslist = [gaslist; [101 102 103]'];

whos gaslist d* w

[mm,nn] = size(d1616);
[mm,nn15] = size(d1515);

iX = input('Choose gasID to examine (-1 to exit) : ');
while iX > 0
  ii = intersect(iX,gaslist);
  ii1515 = ii;
  if iX == 101
    ii1515 = nn15-2;
  elseif iX == 102
    ii1515 = nn15-1;
  elseif iX == 103
    ii1515 = nn15;
  end
  
  if length(ii) == 1
    ii = find(gaslist == iX);
    figure(1); semilogy(w,d1212(:,ii),'b',w,d1612(:,ii),'g',w,d1616(:,ii),'r',w,d1515(:,ii1515),'k'); grid; title(num2str(gaslist(ii)));
      hl = legend('1212','1612','1616','G1515');   set(hl,'fontsize',10)
    figure(2); plot(w,d1612(:,ii) - d1212(:,ii),'g',w,d1616(:,ii) - d1212(:,ii),'r'); grid; title(num2str(gaslist(ii)));
      hl = legend('1612-1212','1616-1212');   set(hl,'fontsize',10)  
    figure(3); plot(w,d1612(:,ii)./d1212(:,ii),'g',w,d1616(:,ii)./d1212(:,ii),'r'); grid; title(num2str(gaslist(ii)));
      hl = legend('1612/1212','1616/1212');   set(hl,'fontsize',10)
    figure(4); plot(w,d1515(:,ii1515)./d1616(:,ii),'g',w,d1212(:,ii)./d1616(:,ii),'r'); grid; title(num2str(gaslist(ii)));
      hl = legend('G1515/H1616','H1212/H1616');   set(hl,'fontsize',10)      
      ax = axis;
      if ax(4) > 2
        ax(4) = 2;
      end
      if ax(3) < 0.5
        ax(3) = 0.5;
      end
      axis(ax);
  else
    fprintf(1,'%3i not part of H2016 gaslist \n',iX);
  end
  iX = input('Choose gasID to examine (-1 to exit) : ');  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('now showing all gases')
for ii = 1 : nn
  figure(1); semilogy(w,d1212(:,ii),'b',w,d1612(:,ii),'b',w,d1616(:,ii),'r'); grid; title(num2str(gaslist(ii)));
    hl = legend('1212','1612','1616');   set(hl,'fontsize',10)
  figure(2); plot(w,d1612(:,ii) - d1212(:,ii),'g',w,d1616(:,ii) - d1212(:,ii),'r'); grid; title(num2str(gaslist(ii)));
    hl = legend('1612-1212','1616-1212');   set(hl,'fontsize',10)  
  figure(3); plot(w,d1612(:,ii)./d1212(:,ii),'g',w,d1616(:,ii)./d1212(:,ii),'r'); grid; title(num2str(gaslist(ii)));
    hl = legend('1612/1212','1616/1212');   set(hl,'fontsize',10)
    ax = axis;
    if ax(4) > 2
      ax(4) = 2;
    end
    if ax(3) < 0.5
      ax(3) = 0.5;
    end
    axis(ax);
    
  fprintf(1,'%2i of %2i : gid = %3i \n',ii,nn,gaslist(ii));
  %disp('ret'); pause;
  pause(1)
end  