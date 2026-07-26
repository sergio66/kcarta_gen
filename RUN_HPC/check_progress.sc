date; ls -lt JUNK/ind*jac.mat | wc -l; squeue -u sergio  | wc -l

while true; do ls JUNK/individual_prof_convolved_kcarta_* | wc -l ; sleep 15; done
