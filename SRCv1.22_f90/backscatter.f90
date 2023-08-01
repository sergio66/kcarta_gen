! backscatter coeffs, see Maestr/Martinazzo Journal of Quantitative Spectroscopy & Radiative Transfer 271 (2021) 107739
! Assessment of the accuracy of scaling methods for radiance simulations at far and mid infrared wavelengths
! Michele Martinazzoa, Davide Magurnoa, William Cossicha, Carmine Seriob, Guido Masiellob, Tiziano Maestri

! ---------------> now add on the backscattered part <--------------------
IF ((iScaling .EQ. +1) .OR. (iScaling .EQ. +4)) THEN
  !! this is SIMILARITY SCALING, same as iScaling == -10
  !! Tang 2018, Eqn 15    dt' = dt (1-f) = dt(1-w(1-b))
  b = (1.0 + ASYM_RTSPEC(L))/2.0               

  !! CHOU SCALING TERM with similarity scaling adjustment
  !! Tang 2018, Eqn 14    dt' = dt(1-w/2(1+g))
  TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*b) 

! ELSEIF (iScaling .EQ. -10) THEN  
!   !! this is SIMILARITY SCALING, same as iScaling == +1
!   !! Chou 1999, Eqn 11 with i=1 (ignore other terms), till May 2022 incorect?
!   b = (1.0 - ASYM_RTSPEC(L))/2.0               
!   b = 1.0 - (1.0 + ASYM_RTSPEC(L))/2.0              
!
!   !! CHOU SCALING TERM
!   !! Chou 1999, Eqn 12,13 dt' = dt (1-f) = dt(1-w(1-b))
!   !! Tang 2018, Eqn 14    dt' = dt (1-f) = dt(1-w(1-b))
!   TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b)) 

ELSEIF ((iScaling .EQ. +2) .OR. (iScaling .EQ. +5)) THEN
  !! the "correct" version is 
  !!  b = 1.0 - (0.5 + 0.3738*ASYM_RTSPEC(L)+0.0076*(ASYM_RTSPEC(L)**2) + 0.1186*(ASYM_RTSPEC(L)**3))   
  !! or b ~~ (1-0.5) - 0.3738*ASYM_RTSPEC(L)
  !!      ~~ 1/2 - (0.3738+0.1262)x +0.1262x
  !!      ~~ (1 - x)/2 + 0.1262x so pretty close to what I used till May 2022

  !! Chou 1999, Eqn 11 with i=1 (ignore other terms), after June 2022          
  b = 1.0 - (0.5 + 0.3738*ASYM_RTSPEC(L)+0.0076*(ASYM_RTSPEC(L)**2) + 0.1186*(ASYM_RTSPEC(L)**3))   

  !! CHOU SCALING TERM, and CHOU with SCALING ADJUSTMENT
  !! Chou 1999, Eqn 12,13 dt' = dt (1-f) = dt(1-w(1-b))
  !! Tang 2018, Eqn 14    dt' = dt (1-f) = dt(1-w(1-b))
  TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b)) 

ELSEIF ((iScaling .EQ. +3) .AND. ((CTYPE .GE. 100) .AND. (CTYPE .LE. 199))) THEN
  !! !!!! uniboWAT scaling, Table 3 of Maertri/Martinazzo paper, JSQRT 2021
  b = 1.0 - (0.5 + 0.2884*ASYM_RTSPEC(L) + 0.5545*(ASYM_RTSPEC(L)**2) - 0.3429*(ASYM_RTSPEC(L)**3))   
  TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b)) 

ELSEIF ((iScaling .EQ. +3) .AND. ((CTYPE .GE. 200) .AND. (CTYPE .LE. 299))) THEN
  !! !!!! uniboICE scaling, Table 3 of Maertri/Martinazzo paper, JSQRT 2021
  b = 1.0 - (0.5 + 0.4452*ASYM_RTSPEC(L) - 0.3189*(ASYM_RTSPEC(L)**2) + 0.3737*(ASYM_RTSPEC(L)**3))   
  TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b)) 

ELSEIF ((iScaling .EQ. +3) .AND. (CTYPE .GE. 300)) THEN
  !! Chou 1999, Eqn 11 with i=1 (ignore other terms), after June 2022, use for AEROSOL
  b = 1.0 - (0.5 + 0.3738*ASYM_RTSPEC(L)+0.0076*(ASYM_RTSPEC(L)**2) + 0.1186*(ASYM_RTSPEC(L)**3))   

  !! CHOU SCALING TERM, and CHOU with SCALING ADJUSTMENT
  !! Chou 1999, Eqn 12,13 dt' = dt (1-f) = dt(1-w(1-b))
  !! Tang 2018, Eqn 14    dt' = dt (1-f) = dt(1-w(1-b))
  TAUTOT_N = TAUTOT_N * (1 - SSALB(L)*(1.0-b)) 

ELSE
  write(kStdErr,*) 'iScaling,CTYPE = ',iScaling,CTYPE
  CALL DoStop
END IF
