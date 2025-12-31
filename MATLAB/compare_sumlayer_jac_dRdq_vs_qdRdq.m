%% this is raw radiance, no layer perturbations. So jacobian is accurate

%% can have dR/dq (kJacobOoutput = -1) or q dR/dq (kJacobOoutput = 0)
[r0,w] = readkcstd('junk.dat');
[j0,w] = readkcjac('junk.jac');
j0gas = j0(:,1:97);

[f0,q0] = quickconvolve(w,r0,1,1);
[fj,qjac] = quickconvolve(w,j0gas,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kJacobOoutput == -1
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  
  %% CO2 jac, 1 layers total
  %%    iJax = 6
  %%    iGasJac = 2
  %%    iJax2 = iJax+0    !!! can alter this to do dQ for many adjacent layers
  
  [r2,w] = readkcstd('junk.dat');
  [f2,q2] = quickconvolve(w,r2,1,1);
  
  figure(1); dq2 = 3.76770E-09; jfinite2 = (r2-r0)/dq2; plot(w,jfinite2,'b.-',w,j0gas(:,3)); legend('finite diff','analytic')
  figure(2); jxfinite2 = (q2-q0)/dq2; plot(f0,jxfinite2,'b.-',f0,qjac(:,3)); legend('finite diff','analytic')
  figure(1); xlim([763 764]); ax = axis; figure(2); axis(ax)
  
  figure(1); dq2 = 3.76770E-09; jfinite2 = (r2-r0)/dq2; plot(w,jfinite2,'b.-',w,j0gas(:,3)*1.04); legend('finite diff','analytic * 1.04 since scale factor for iNewGas')
  figure(2); jxfinite2 = (q2-q0)/dq2; plot(f0,jxfinite2,'b.-',f0,qjac(:,3)*1.04); legend('finite diff','analytic * 1.04 since scale factor for iNewGas')
  
  figure(2); jxfinite2 = (q2-q0)/dq2; plot(f0,jxfinite2./qjac(:,3)); legend('finite diff/analytic')
  ylim([0 2]); grid on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  
  %% now perturb eg 5 layers
  %% subroutine Set_Ref_Current_Profs in kcartamisc.f90
  %%
  %% WV jac, 5+1 = 6 layers total
  %%    iJax = 52
  %%    iGasJac = 1
  %%    iJax2 = iJax+5    !!! can alter this to do dQ for many adjacent layers
  %%
  
  [r2,w] = readkcstd('junk.dat');
  [f2,q2] = quickconvolve(w,r2,1,1);
  
  %% this is the text that is displayed from subroutine Set_Ref_Current_Profs, from jac (or radiating atmospere) layers 49:54
  Q = [1.39081E-09 1.13613E-09 9.83971E-10 8.42834E-10 7.61239E-10 7.18113E-10];
    dQ = Q/100;                                                                     %% when I put in constant dumb dQ/Q perturbation for kJacobOoutput == -1
    dQ = [1.39081E-11 1.39081E-11 1.39081E-11 1.39081E-11 1.39081E-11 1.39081E-11]; %% when I put in smart constant dQ  perturbation for kJacobOoutput == -1
      
  wahjac = qjac(:,49:54); dQjac = ones(length(f2),1)*dQ;  wah2jac = sum(wahjac.*dQjac,2); whos wahjac dQjac wah2jac
  dq2 = sum(dQ);

  figure(1); jfinite2 = (r2-r0)/dq2; plot(w,jfinite2,'b.-',w,sum(j0gas(:,49:54),2)/6);      legend('finite diff','analytic')
  figure(2); jxfinite2 = (q2-q0)/dq2; plot(f0,jxfinite2,'b.-',f0,sum(qjac(:,49:54),2)/6)    %% this shows raw finite diff jac, and incorrectly summed jac
  figure(2); jxfinite2xx = (q2-q0);   plot(f0,jxfinite2xx,'b.-',f0,wah2jac,'g')             %% this shows delta finite diff radiance, and correct summation of jac x dq
  
end if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if kJacobOoutput == 0
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  %% if you only have one radiating layer eg gas = 2, iJax = 6,iJax2 = 6
  
  %% CH4 jac, 1 layers total
  %%    iJax = 52
  %%    iGasJac = 6
  %%    iJax2 = iJax+0    !!! can alter this to do dQ for many adjacent layers
  
  [r2,w] = readkcstd('junk.dat');
  [f2,q2] = quickconvolve(w,r2,1,1);

  figure(1); dq2 = 5.14660E-12; jfinite2 = (r2-r0)/dq2; plot(w,jfinite2,'b.-',w,j0gas(:,49)); legend('finite diff','analytic')
  figure(2); jxfinite2 = (q2-q0)/dq2; plot(f0,jxfinite2,'b.-',f0,qjac(:,49)); legend('finite diff','analytic')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  %% if you have eg five radiating layers eg gas = 1, iJax = 52,iJax2 = 57
  
  %% now perturb eg 5 layers
  %% subroutine Set_Ref_Current_Profs in kcartamisc.f90
  %%
  %% CH4 jac, 5+1 = 6 layers total
  %%    iJax = 52
  %%    iGasJac = 6
  %%    iJax2 = iJax+5    !!! can alter this to do dQ for many adjacent layers
  %%
  
  [r2,w] = readkcstd('junk.dat');
  [f2,q2] = quickconvolve(w,r2,1,1);
  
  %% this is the text that is displayed from subroutine Set_Ref_Current_Profs, from jac (or radiating atmospere) layers 49:54
  Q = [5.14660E-10 4.91849E-10 4.69499E-10 4.47608E-10 4.25757E-10 4.04108E-10]; deltaQ = Q/100
  %% rewrite dQ
  dQ = deltaQ./Q;  %% [0.01 0.01 0.01 0.01 0.01 0.01]
  
  wahjac = qjac(:,49:54); dQjac = ones(length(f2),1)*dQ;  wah2jac = sum(wahjac.*dQjac,2); whos wahjac dQjac wah2jac

  figure(1); dq2 = 0.01; jfinite2 = (r2-r0)/log(1+dq2); plot(w,jfinite2,'b.-',w,sum(j0gas(:,49:54),2)); legend('finite diff','analytic')
  figure(2); jxfinite2 = (q2-q0)/log(1+dq2);            plot(f0,jxfinite2,'b.-',f0,sum(qjac(:,49:54),2))                  %% this shows raw finite diff jac, and correctly summed jac since this is dR/dlinQ
  
end if
