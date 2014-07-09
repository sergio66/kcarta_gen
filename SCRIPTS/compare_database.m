database_old = '/home/sergio/KCARTADATA/General/compHT2007_extended0909.param_fromASL';
dold1 = load(database_old);
database_old = '/home/sergio/KCARTADATA/General/xsecHT2007_extended0909.param_fromASL';
dold2 = load(database_old);

dold = [dold1; dold2];

database_new = '/home/sergio/KCARTA/SCRIPTS/comp_ir605_2830.param';
dnew = load(database_new);

woof = find(dold(:,4) == 20); dold = dold(woof,:);
woof = find(dnew(:,4) == 20); dnew = dnew(woof,:);

gf = 605 : 25 : 2830;
gf = 605 : 25 : 2805;

g1 = 1; g2 = 81;
g1 = 1; g2 = 81;

for ii = g1 : g2
  boo1 = find(dold(:,1) == ii);
  boo2 = find(dnew(:,1) == ii);
  
  figure(1); clf

  if (length(boo1) > 0 & length(boo2) > 0)
    gasfound = find(dold(:,1) == ii); 
    if length(gasfound) > 0
      doldgas = dold(gasfound,:);
      iFound = -1;
      for jj = 1 : length(gf)
        olddatabase(jj) = -1;
        for kk = 1 : length(doldgas(:,1))
          b1 = doldgas(kk,2); b2 = doldgas(kk,3);
          if gf(jj) >= b1 & gf(jj) < b2
            olddatabase(jj) = +1;
            end
          end
        end
      end

    gasfound = find(dnew(:,1) == ii); 
    if length(gasfound) > 0
      dnewgas = dnew(gasfound,:);
      iFound = -1;
      for jj = 1 : length(gf)
        newdatabase(jj) = -1;
        for kk = 1 : length(dnewgas(:,1))
          b1 = dnewgas(kk,2); b2 = dnewgas(kk,3);
          if gf(jj) >= b1 & gf(jj) < b2
            newdatabase(jj) = +1;
            end
          end
        end
      end

    plot(gf,olddatabase,'bo',gf,newdatabase+0.1,'rx')
    title(num2str(ii))
    pause(2)
    end
  end

