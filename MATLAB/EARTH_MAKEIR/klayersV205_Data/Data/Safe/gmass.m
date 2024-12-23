% Program gmass
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atomic mass from American Institute Of Physics Handbook (1963) Section7-9
H = 1.00797;
C = 12.01115;
N = 14.0067;
O = 15.9994;
F = 18.9984;
S = 32.064;
Cl = 35.453;
Br = 79.909;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=H+C+O+O+H;
fprintf(1,'gas ID=32 HCOOH mass=%9.4f\n',m);
disp(' ')

m=H+O+O;
fprintf(1,'gas ID=33 HO2 mass=%9.4f\n',m);
disp(' ')

m=O;
fprintf(1,'gas ID=34 O mass=%9.4f\n',m);
disp(' ')

m=Cl+O+N+O*2;
fprintf(1,'gas ID=35 ClONO2 mass=%9.4f\n',m);
disp(' ')

m=N+O;
fprintf(1,'gas ID=36 NO+ mass=%9.4f\n',m);
disp(' ')

m=H+O+Br;
fprintf(1,'gas ID=37 HOBr mass=%9.4f\n',m);
disp(' ')

m=C*2+H*4;
fprintf(1,'gas ID=38 C2H4 mass=%9.4f\n',m);
disp(' ')

m=C+H*3+O+H;
fprintf(1,'gas ID=39 CH3OH mass=%9.4f\n',m);
disp(' ')

m=C+H*3+Br;
fprintf(1,'gas ID=40 CH3Br mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+N;
fprintf(1,'gas ID=41 CH3CN mass=%9.4f\n',m);
disp(' ')

m=C+F*4;
fprintf(1,'gas ID=42 CF4 mass=%9.4f\n',m);
disp(' ')

disp('%%%%%%%%%%%%%%%% xsec %%%%%%%%%%%%%%%%')

m=C*2+F*6;
fprintf(1,'gas ID=64 C2F6 mass=%9.4f\n',m);
disp(' ')

m=C+H+Cl*2+C+F*3;
fprintf(1,'gas ID=65 CHCl2CF3 mass=%9.4f\n',m);
disp(' ')

m=C+H+Cl+F+C+F*3;
fprintf(1,'gas ID=66 CHClFCF3 mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+Cl*2+F;
fprintf(1,'gas ID=67 CH3CCl2F mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+Cl+F*2;
fprintf(1,'gas ID=68 CH3CClF2 mass=%9.4f\n',m);
disp(' ')

m=C+H+Cl*2+C+F*2+C+F*3;
fprintf(1,'gas ID=69 CHCl2CF2CF3 mass=%9.4f\n',m);
disp(' ')

m=C+Cl+F*2+C+F*2+C+H+Cl+F;
fprintf(1,'gas ID=70 CClF2CF2CHClF mass=%9.4f\n',m);
disp(' ')

m=C+H*2+F*2;
fprintf(1,'gas ID=71 CH2F2 mass=%9.4f\n',m);
disp(' ')

m=C+F+H*2+C+F*3;
fprintf(1,'gas ID=72 CFH2CF3 mass=%9.4f\n',m);
disp(' ')

m=C+F*3+C+H*3;
fprintf(1,'gas ID=73 CF3CH3 mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+H+F*2;
fprintf(1,'gas ID=74 CH3CHF2 mass=%9.4f\n',m);
disp(' ')

m=C*6+H*6;
fprintf(1,'gas ID=75 C6H6 mass=%9.4f\n',m);
disp(' ')

m=C+H+F*2+C+F*3;
fprintf(1,'gas ID=76 CHF2CF3 mass=%9.4f\n',m);
disp(' ')

m=C+H+F*2+C+H+F*2;
fprintf(1,'gas ID=77 CHF2CHF2 mass=%9.4f\n',m);
disp(' ')

m=S+F*5+C+F*3;
fprintf(1,'gas ID=78 SF5CF3 mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+O+O+O+N+O*2;
fprintf(1,'gas ID=79 CH3C(O)OONO2 mass=%9.4f\n',m);
disp(' ')

m=C+H*3+C+N;
fprintf(1,'gas ID=80 CH3CN mass=%9.4f\n',m);
disp(' ')


%%% end of program %%%

