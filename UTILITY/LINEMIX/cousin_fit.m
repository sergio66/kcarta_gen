path(path,'/home/sergio/SPECTRA/');

clf;
T = [150 175 200 225 250 275 300 325 350]; T=T';

pt = 1.0; ps = 3.6e-4; pf = pt - ps;
%% do the +ve freqs 
freq = 0.0;
w_tot = 0.005;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% f = 0 --> 0.5 ==> chi = 1.0

for ii = 1 : 9
 if (ii == 1) f = 0.5:0.1:3; end;
 if (ii == 2) f = 3:0.1:5;   end;
 if (ii == 3) f = 5:0.1:9;   end;
 if (ii == 4) f = 9:0.1:20;  end;
 if (ii == 5) f = 20:0.1:22; end;
 if (ii == 6) f = 22:0.1:23; end;
 if (ii == 7) f = 23:0.1:25; end;
 if (ii == 8) f = 25:1:50;   end;
 if (ii == 9) f = 50:1:300;  end;

  if (f(1) < 1)
    f1 = 1;
  else
    f1 = round(f(1));
    end
  fname = ['/home/sergio/KCARTA/SRCv1.11/NONLTE/cousin_' num2str(f1) '.txt']

  chip150 = cousin1(f,freq,w_tot,150,pf,ps);
  chip175 = cousin1(f,freq,w_tot,175,pf,ps);
  chip200 = cousin1(f,freq,w_tot,200,pf,ps);
  chip225 = cousin1(f,freq,w_tot,225,pf,ps);
  chip250 = cousin1(f,freq,w_tot,250,pf,ps);
  chip275 = cousin1(f,freq,w_tot,275,pf,ps);
  chip300 = cousin1(f,freq,w_tot,300,pf,ps);
  chip325 = cousin1(f,freq,w_tot,325,pf,ps);
  chip350 = cousin1(f,freq,w_tot,350,pf,ps);

  chip=[chip150;chip175;chip200;chip225;chip250;...
        chip275;chip300;chip325;chip350];
  plot(f,chip);
  axis tight

  if (f(1) < 5)
    Tpower = 2;
  elseif (f(1) < 50)
    Tpower = 3;
  else
    Tpower = 4;
    end
  yp150 = polyfit(f,chip150,Tpower); 
  yp175 = polyfit(f,chip175,Tpower);
  yp200 = polyfit(f,chip200,Tpower);
  yp225 = polyfit(f,chip225,Tpower);
  yp250 = polyfit(f,chip250,Tpower);
  yp275 = polyfit(f,chip275,Tpower);
  yp300 = polyfit(f,chip300,Tpower);
  yp325 = polyfit(f,chip325,Tpower);
  yp350 = polyfit(f,chip350,Tpower);
  yp = [yp150;yp175;yp200;yp225;yp250;yp275;yp300;yp325;yp350];
  plot(f,chip150,f,polyval(yp150,f),'r'); pause(0.5)
  plot(f,chip175,f,polyval(yp175,f),'r'); pause(0.5)
  plot(f,chip200,f,polyval(yp200,f),'r'); pause(0.5)
  plot(f,chip225,f,polyval(yp225,f),'r'); pause(0.5)
  plot(f,chip250,f,polyval(yp250,f),'r'); pause(0.5)
  plot(f,chip275,f,polyval(yp275,f),'r'); pause(0.5)
  plot(f,chip300,f,polyval(yp300,f),'r'); pause(0.5)
  plot(f,chip325,f,polyval(yp325,f),'r'); pause(0.5)
  plot(f,chip350,f,polyval(yp350,f),'r'); pause(0.5)

  npower = 3;
  semilogy(T,abs(yp));
  if (f(1) < 5)
    zp1=polyfit(T,yp(:,1),npower);
    zp2=polyfit(T,yp(:,2),npower);
    zp3=polyfit(T,yp(:,3),npower);
    ZP = [zp1;  zp2; zp3];
    ZM = [zp1; -zp2; zp3];
  elseif (f(1) < 50)
    zp1 = polyfit(T,yp(:,1),npower);
    zp2 = polyfit(T,yp(:,2),npower);
    zp3 = polyfit(T,yp(:,3),npower);
    zp4 = polyfit(T,yp(:,4),npower);
    ZP = [zp1; zp2; zp3; zp4];
    ZM = [zp1; -zp2; zp3; -zp4];
  else
    zp1 = polyfit(T,yp(:,1),npower);
    zp2 = polyfit(T,yp(:,2),npower);
    zp3 = polyfit(T,yp(:,3),npower);
    zp4 = polyfit(T,yp(:,4),npower);
    zp5 = polyfit(T,yp(:,5),npower);
    ZP = [zp1; zp2; zp3; zp4; zp5];
    ZM = [zp1; -zp2; zp3; -zp4; zp5];
    end
  fid = fopen(fname,'w');
  fprintf(fid,'%12.8e %12.8e %12.8e %12.8e \n',ZP');
  fclose(fid);

  end

error('oops');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cousin is symmetric, so c(f-f0) = c(f0-f)
%ZM = [zp1; -zp2; zp3; -zp4];
%%%%%%%%%%%%%%%%%%%%%%%%%%%% this is to test the coeffs %%%%%%%%%%%%%%

F = 1 : 0.1: 300;

ii0 = find(F >= 0.0  & F < 0.5);   f0 = F(ii0); z0 = ones(size(f0));
ii1 = find(F >= 0.5  & F < 3.0);   f1 = F(ii1); z1 = ones(size(f1));
ii2 = find(F >= 3.0  & F < 5.0);   f2 = F(ii2); z2 = ones(size(f2));
ii3 = find(F >= 5.0  & F < 9.0);   f3 = F(ii3); z3 = ones(size(f3));
ii4 = find(F >= 9.0  & F < 20.0);  f4 = F(ii4); z4 = ones(size(f4));
ii5 = find(F >= 20.0 & F < 22.0);  f5 = F(ii5); z5 = ones(size(f5));
ii6 = find(F >= 22.0 & F < 23.0);  f6 = F(ii6); z6 = ones(size(f6));
ii7 = find(F >= 23.0 & F < 25.0);  f7 = F(ii7); z7 = ones(size(f7));
ii8 = find(F >= 25.0 & F < 50.0);  f8 = F(ii8); z8 = ones(size(f8));
ii9 = find(F >= 50.0);             f9 = F(ii9); z9 = ones(size(f9));

T = 180;
ffit = [f0 f1 f2 f3 f4 f5 f6 f7 f8 f9];

zfit = z0;

for ii = 1 : 9
 if (ii == 1) f = 0.5:0.1:3; fp = f1; end;
 if (ii == 2) f = 3:0.1:5;   fp = f2; end;
 if (ii == 3) f = 5:0.1:9;   fp = f3; end;
 if (ii == 4) f = 9:0.1:20;  fp = f4; end;
 if (ii == 5) f = 20:0.1:22; fp = f5; end;
 if (ii == 6) f = 22:0.1:23; fp = f6; end;
 if (ii == 7) f = 23:0.1:25; fp = f7; end;
 if (ii == 8) f = 25:1:50;   fp = f8; end;
 if (ii == 9) f = 50:1:300;  fp = f9; end;

  if (f(1) < 1)
    f1 = 1;
  else
    f1 = round(f(1));
    end
  fname = ['/home/sergio/KCARTA/SRCv1.11/NONLTE/cousin_' num2str(f1) '.txt'];

  zz = load(fname);
  if ii <= 2
    zp1 = zz(1,:);
    zp2 = zz(2,:);
    zp3 = zz(3,:);
    zp = polyval(zp1,T)*(fp.^2) + polyval(zp2,T)*(fp.^1) + polyval(zp3,T);
  elseif (ii <= 8)
    zp1 = zz(1,:);
    zp2 = zz(2,:);
    zp3 = zz(3,:);
    zp4 = zz(4,:);
    zp = polyval(zp1,T)*(fp.^3) + polyval(zp2,T)*(fp.^2) + ...
       polyval(zp3,T)*(fp) + polyval(zp4,T);
  else
    zp1 = zz(1,:);
    zp2 = zz(2,:);
    zp3 = zz(3,:);
    zp4 = zz(4,:);
    zp5 = zz(5,:);
    zp = polyval(zp1,T)*(fp.^4) + polyval(zp2,T)*(fp.^3) + ...
       polyval(zp3,T)*(fp.^2) + polyval(zp4,T)*(fp) + polyval(zp5,T);
    end

  zfit = [zfit zp];
  end

zcorrect = cousin1(F,freq,w_tot,T,pf,ps);
%plot(ffit,zfit,F,zcorrect,'r')
%whos F zfit ffit zcorrect 
%setdiff(F,ffit)
subplot(211); plot(ffit,zfit,F,zcorrect,'r')
subplot(212); plot(F,(zfit-zcorrect)./zcorrect);


