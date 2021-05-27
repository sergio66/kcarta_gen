grep -in '(iFr' *.f90 >& ughX1
sed '/kbloat.f90/d'               ughX1 > ughX2
sed '/klinemix.f90/d'             ughX2 > ughX3
sed '/knonlte.f90/d'              ughX3 > ughX4
sed '/kvoigt_cousin.f90/d'        ughX4 > ughX5
sed '/scatter_graycld_code.f90/d' ughX5 > ughX6
sed '/scatter_rtspec_code.f90/d'  ughX6 > ughX7
sed '/scatter_rtspec_flux.f90/d'  ughX7 > ughX8

mv ughX8 ugh
rm ughX*