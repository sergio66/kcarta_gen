function r = find_rad_dnwell(w,r0,od,mu,iLayNum,playstemp,plevstemp,iNlevs,iVers)

%% does rt through one layer, dnwelling
%% input
%%   w = wavenumber
%%  r0 = initial rad at bottom of layer
%%  od = layer OD
%%  mu = cos(theta)
%%  iLayNum = layer number
%%  playstemp = array containing layer temperatures
%%  plevstemp = array containing level temperatures
%%  iNlevs    = number of levels in profile
%%  iVers     = which version of RT to do
%%              -1 : contant layer temp
%%               4  : linear in tau, small OD accurate to O(od^2)
%%               42 : linear in tau, small OD accurate to O(od)
%%               41 : linear in tau, pade
%% output
%%  f = final rad

iBeta = mod(iLayNum,iNlevs-1);
if iBeta == 0
  iBeta = iNlevs-1;
elseif iLayNum == iNlevs
  iBeta = iNlevs-1;
end

fprintf(1,'%3i %3i %8.6f %8.6f \n',iLayNum,iBeta,plevstemp(iBeta+1),plevstemp(iBeta))

if iVers == -1
  planckrad = ttorad(w,playstemp(iLayNum));
  r = r0.*exp(-od/mu) + planckrad .* (1 - exp(-od/mu));
elseif iVers == 4
  %% note these are "flipped" wrt f77 kcarta as there the lowest layesr == 1, here lowest layer = 100
  raIntenP  = ttorad(w,plevstemp(iBeta+1));
  raIntenP1 = ttorad(w,plevstemp(iBeta));
  raIntenAvg = ttorad(w,playstemp(iBeta));
  rAbs = od/mu;
  rZeta = 2*(raIntenAvg-raIntenP);
  boo = find(rAbs > 0.1);
  if length(boo) > 0
    rTrans(boo) = exp(-rAbs(boo));
    rFcn(boo) = (1-rTrans(boo)').*(raIntenP1(boo) + rZeta(boo)./rAbs(boo)) - rTrans(boo)' .* rZeta(boo);
  end
  boo = find(rAbs <= 0.1);
  if length(boo) > 0  
    rTrans(boo) = 1 - rAbs(boo) + 0.5*rAbs(boo).^2;
    rZeta2(boo) = rZeta(boo).*(rAbs(boo)/2-(rAbs(boo).^2)/3+(rAbs(boo).^3)/6);
    rFcn(boo) = (1-rTrans(boo)').*raIntenP1(boo) + rZeta2(boo);
  end
  r = r0'.*rTrans + rFcn;
elseif iVers == 41
  error('ooo')
elseif iVers == 42
  %% note these are "flipped" wrt f77 kcarta as there the lowest layesr == 1, here lowest layer = 100
  raIntenP  = ttorad(w,plevstemp(iBeta+1));
  raIntenP1 = ttorad(w,plevstemp(iBeta));
  raIntenAvg = ttorad(w,playstemp(iBeta));
  rAbs = od/mu;
  rZeta = 2*(raIntenAvg-raIntenP);
  boo = find(rAbs > 0.05);
  if length(boo) > 0  
    rTrans(boo) = exp(-rAbs(boo));
    rFcn(boo) = (1-rTrans(boo)').*(raIntenP(boo) + rZeta(boo)./rAbs(boo)) - rTrans(boo)' .* rZeta(boo);
  end
  boo = find(rAbs <= 0.05);
  if length(boo) > 0  
    rTrans(boo) = 1 - rAbs(boo);
    rFcn(boo) = rAbs(boo).*raIntenP(boo) + rZeta(boo).*(1-rAbs(boo)/2) - rTrans(boo)' .* rZeta(boo);
  end
  r = r0'.*rTrans + rFcn;
else
  error('unknown iVers')
end
