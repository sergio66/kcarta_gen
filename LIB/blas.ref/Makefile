# ***** Linux with Absoft Fortran *****
#
# compiler options
# -f	fold all names to lower case
# -W	wide 
# -w	suppress warning messages (absoft is more fussy than g77)
# -A	alignment warnings
# -C	check array bounds
# -O    optimizations
# -N3   add record info to unformatted files
# -N109 fold all names to upper case
# -N15  so that LAPACK lib can be included; adds trailing underscores
#
#F77 = /usr/bin/f77
#F77 = f77
#F77 = /asl/opt/absoft/absoft10.0/bin/af77
#FFLAGS = -w -W -C -N3 -N109
#FFLAGS = -w -W -O -N3 -N109
FFLAGS = -w -W -O3 -N3 -f -N15

#Linux with PGF Fortran
#----------------------
#F77 = pgf77
####FFLAGS = -Mextend -Mbounds -Mupcase -O2
####FFLAGS = -Mextend -Mnoupcase -O3
#FFLAGS = -Mextend -Mnoupcase -O3 -Mbounds

#Linux with ifort
#----------------------
#F77 = /cm/shared/apps/intel/composer_xe/2015.1.133/bin/intel64/ifort
F77 = ifort
FFLAGS = -names lowercase -extend-source 132 -check bounds -msse2 -nothreads -noreentrancy

#linux with gfortran
#-------------------------
#F77 = /usr/ebuild/software/GCCcore/11.2.0/bin/gfortran
#FFLAGS = -ffixed-line-length-none -O3 -fbounds-check

libblas.a:
	$(F77) -c $(FFLAGS) *.f
	ar -rcs libblas.a *.o

clean:
	rm libblas.a 2> /dev/null || true

