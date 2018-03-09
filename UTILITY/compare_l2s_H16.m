f1616 = 'l2s_kc120_H16.dat';          %% molgas = H2016, xscgas = H2016
f1612 = 'l2s_kc120_H16_xsecH12.dat';  %% molgas = H2016, xscgas = H2012
f1212 = 'l2s_kc120_H12.dat';          %% molgas = H2012, xscgas = H2012

addpath /home/sergio/KCARTA/MATLAB

[d1616,w] = readkcstd(f1616);
[d1612,w] = readkcstd(f1612);
[d1212,w] = readkcstd(f1212);

gaslist = load('/home/sergio/KCARTA/SCRIPTS/MAKE_COMP_HTXY_PARAM_SC/PARAM_HT2016/comp_ir605_2830.param');
gaslist = unique(gaslist(:,1));
gaslist = [gaslist; [101 102 103]'];

whos gaslist d* w

[mm,nn] = size(d1616);

iX = input('Choose gasID to examine (-1 to exit) : ');
while iX > 0
  ii = intersect(iX,gaslist);
  if length(ii) == 1
    ii = find(gaslist == iX);
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