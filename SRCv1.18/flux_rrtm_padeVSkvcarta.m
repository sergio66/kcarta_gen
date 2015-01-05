info = load('/home/sergio/IR_NIR_VIS_UV_RTcodes/RRTM/v3.3/rrtm_lw/MATLAB/INPUT_RRTM_quickdump');

pinfo = info(:,2);
tm1 = info(:,4);
tav = info(:,2);
tp1 = info(:,6);

wn = 500 : 1: 3000;
OD = -3:1:2; OD = 10.^OD;

a = 0.193;
b = 0.013;

for ii = 1 : length(OD)
  for jj = 1 : length(tm1)
    tau = OD(ii);
    bm1 = ttorad(wn,tm1(jj));
    bp1 = ttorad(wn,tp1(jj));
    bav = 0.5 * (bm1 + bp1);

    num   = bav + (a*tau + b*tau*tau)*bp1;
    denom = 1 + a*tau + b*tau*tau;
    rrtm_pade(ii,jj,:) = num./denom;

    trans = exp(-tau);
    kcarta_exact(ii,jj,:) = bp1*(1-trans) - (bm1-bp1)*trans + (bm1-bp1)/tau*(1-trans);
  end
end

for ii = 1 : length(OD)
  figure(6); plot(wn,squeeze(kcarta_exact(ii,:,:)));
    title([' OD = ' num2str(OD(ii))])
  figure(7); plot(wn,squeeze(kcarta_exact(ii,:,:) ./ rrtm_pade(ii,:,:)));
    title([' OD = ' num2str(OD(ii))])
  disp('ret to continue');
  pause
end;