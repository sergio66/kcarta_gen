

       kCARTA:  Setting  up  the


                    distribution



     Sergio De Souza-Machado, L. Larrabee Strow,
          Howard Motteler, and Scott Hannon



University of Maryland Baltimore County, Baltimore, MD 21250 USA


                         Copyright 1997
             University of Maryland Baltimore County
                       All Rights Reserved
                      v0.96 March 18, 1998
                                1


Contents


1  Introduction                                                     3


2  Directories                                                       3
   2.1 Files in SRC subdirectory  . . . . . . . . . . . . . . . . . . . . . .  *
 * 4
   2.2 Files in MATLAB subdirectory  . . . . . . . . . . . . . . . . . . .   5
   2.3 Files in UTILITY subdirectory  . . . . . . . . . . . . . . . . . . .   6
   2.4 Files in COMPARISON subdirectory   . . . . . . . . . . . . . . .   6
   2.5 Files in SCRIPTS subdirectory  . . . . . . . . . . . . . . . . . . .   7


3  How to generate data for testing                                 8


4  Settings of parameter files for tests                              9
                                  2


1   Introduction


This is a beta version of the kCARTA atmospheric radiative transfer code
produced by the Atmospheric Spectroscopy Laboratory at the University
of Maryland Baltimore County (UMBC). The Atmospheric Spectroscopy
Laboratory is part of the Joint Center for Earth Systems Technology
(JCET) which is a joint center between UMBC and NASA/Goddard Space
Flight Center. Please contact L. Strow (strow@umbc.edu) or S. DeSouza-
Machado (sergio@umbc.edu) with any questions.
   More detailed information on using kCARTA, including specifics of
driver files and theory of operation, can be found in kcarta.ps;  this
"readme" currently contains some additional information about utility
programs, files, and testing the kCARTA distribution, but it may not be
entirely up-to-date, and will eventually be subsumed by kcarta.ps, and
by the top-level README.1ST, file, which has more information about
compiling and linking the programs.
   Caveats:  As stated above, this is a beta version, so use it at your
own risk.  We would appreciate hearing about problems in running or
using these programs. The documentation is in a very rough state, and
we will be working on it in the future.  The NetCDF routines have not
been as well tested as the rest of the distribution. There are some small
errors on integer wavenumber boundaries that will be fixed in the fu-
ture. These errors are inconsequential if you are convolving with a re-
alistic instrument function.  The absorption coefficients used to gener-
ate the compressed database are mostly based on HITRAN92, and were
computed with GENLN2V3.0 (David Edwards, NCAR, edwards@ucar.edu).
The source for this code may undergo extensive revision in the future,
so please don't count on future versions of kCARTA being compatible.
   If kCARTA is used in any reports or publications, we would appreci-
ate a reference.  We have included, in the DOC directory, a preprint of
an article submitted to JQSRT on the kCARTA algorithm.  If you need
additional references (SPIE paper) please contact L. Strow.
2   Directories


These directory/subdirectories contain the k-compressed code, tests and
so on.  Many of the subdirectories have additional "readme" files, de-
scribing further what is contained in the files. A brief synopsis of each
subdirectory follows.
                                  3


BIN: executable files from running "make" in /SRC; executable files from
     running "make" in /UTILITY


COMPARISON:    results of some test runs, to compare against


DOC:  documentation, as well as some papers


DATA:  contains the following subdirectories


     CompDataBase:    k-compressed (RestOfGases) database

     WaterDataBase:  k-compressed (Water) database

     RefProf: reference profiles

     TestProf: regression profiles

     Template:  template input files

     General: cross sectional gas database, template files


OUTPUT:   subdirectory where output from kcarta.sc goes to


RUN:  run the executable kcarta.x from here


SCRIPTS:  kcarta.sc produces output to compare against the contents of
     the COMPARISON directory; comp.sc produces summary of avail-
     able k-compressed files


SRC: main source code


UTILITY:  auxiliary Fortran source files_readers, comp.f, etc.


MATLAB:   MatLab files to read in the binary output from kcarta.x



2.1  Files in SRC subdirectory


The main source code is stored in this subdirectory, as are some parame-
ter files (kcarta.param, airslevels.param, airheights.param). The supplied
Makefile produces the executable kcarta.x and installs it in the BIN direc-
tory.


kcarta.param: parameter declarations


kcartamain.f: main driver file


strings3.f:necessary string processing routines
                                  4


strings.f:additional string processing routines


stringinput.f:additional string processing routines


rdprof.f:file I/O routines


kcoeff*.f:k-compressed UNPACK routines


radiance.f:forward model routines


thermalbackgnd.f:  supplementary radiance routines


jacob.f:jacobian code


misc.f:miscellaneous subroutines


netcdf.f:NetCDF routines


The programs calcon.f, calq.f, calcxsc.f, h2oft0.f, and h2ost1.f do the
continuum and cross-section calculations.



2.2  Files in MATLAB subdirectory


The main driver files are


readkcarta.m: reads in the binary output file from kcarta.x. Allows user
     to save all spectra for a particular path, or ALL the data for a par-
     ticular 25 cm 1 set. This file calls in the subroutine *.m files that
     are present in the subdirectory.


readjacob.m:  reads in the jacobian binary output file from kcarta.x. Al-
     lows user to save all spectra for a particular path, or ALL the data
     for a particular 25 cm 1 set. This file calls in the subroutine *.m
     files that are present in the subdirectory.


rdairs.m: this file reads in the binary data file produced by running the
     output of kcarta.x through readkcarta.f (a FORTRAN reader)


loadcdf.m: reads in kCARTA output that has been stored in NetCDF for-
     mat



                                  5


2.3  Files in UTILITY subdirectory


The supplied Makefile compiles these files and stores the executables in
the BIN subdirectory


readkcarta.f:reads in the binary output file from kcarta.x, strips the
     header information and saves raw data in a binary file (that can
     easily be read in by eg MATLAB rdairs.m or separate convolution
     routines)


readjacob .f:reads in the binary jacobian output file from kcarta.x, strips
     the header information and saves raw data in a binary file (that can
     easily be read in by eg MATLAB rdairs.m or separate convolution
     routines)


compdatabase.f:  is used to create the parameter file comp.param, which
     lists the compressed files in the database


makeprofile.f: allows the user to specify which regression profile to use
     for testing (1..48) and then creates a file. It then creates a profile
     file that *PTHFIL uses, containing the path profiles of 36 gases(1..27
     and 51..63). Almost all of the 36 gases have the same pressure/part
     pressure/amounts as the reference profile for the gas.  However,
     the temperature profile comes from the specified test profile.  In
     addition, water (gasID 1) and ozone(gasID 3) amounts are vari-
     able, being set from the specified test profile.  Also, CO2 (gasID
     2) amount = reference amount*(363/330).


makeinp.f: is a file that takes the template KCARTA input file, and pro-
     duces a file tailored to the user requirements : frequency, profile
     and output filename.



2.4  Files in COMPARISON subdirectory


Six frequency ranges were used for testing purposes, as follows

    iFreqIndex = 1   r1=755.0   to   r2=780.0
    iFreqIndex = 2  r1=1005.0   to   r2=1055.0
    iFreqIndex = 3  r1=1230.0   to   r2=1255.0
    iFreqIndex = 4  r1=1430.0   to   r2=1505.0
    iFreqIndex = 5  r1=1530.0   to   r2=1605.0
    iFreqIndex = 6  r1=2355.0   to   r2=2405.0
                                  6


   From makeprofile.f and makeinp.f above, followed by a run of kcarta.x,
some examples of the file names are


testprof7:profiles for 36 gases, using test profile 7


testa7.ip:genln2 input file that uses profile 7,freqrange 1


testa7.dat:resultant kcarta.x binary output file


testprof14: profiles for 36 gases, using test profile 14


teste14.ip:genln2 input file that uses profile 14,freqrange 5


teste14.dat:resultant kcarta.x binary output file


   Note that the present kcarta.sc file deletes the testprof*, test*.ip files
to save disk space; the main information in this files will be stored in the
headers of the binary output files.



2.5  Files in SCRIPTS subdirectory


kcarta.sc:is a script file that ran various programs, including kcarta.x
     to create the data in ../COMPARISON. To be able to run this script,
     you may need to type "chmod +x kcarta.sc" at the UNIX prompt


comp.sc: is a script file that uses compdatabase.x to go through the
     available k-compressed files in DATA/CompDataBase and DATA/WaterDataBase,
     eventually storing the summary result in SRC/comp.param, which
     is a file needed by kcarta.x.  To make this script executable, type
     "chmod +x comp.sc" at the UNIX prompt. Each time you add/remove
     files from the compressed database, this script should be run to
     update comp.param and thus let kcarta.x run merrily.


diffemall.sc:is a script file that runs the data files in OUTPUT thru our
     ../BIN/readkcarta.x reader, saving the output as text files.  These
     text files are then individually diff'ed against the text files in COM-
     PARISON, so the user can see if things worked as expected (there
     might be slight roundoff errors in the radiances).
                                  7


3   How to generate data for testing


  1. Run "make" in /SRC - this creates executables kcarta.x Assuming
     the compilation is successful, the user can type "make clean" at
     the UNIX prompt to remove the *.o files Run "make" in /UTILITY -
     this creates executables makeinput.x, makeprofile.x, readkcarta.x,
     readjacob.x, compdatabase.x All executables are stored in BIN after
     compiling


  2. At this point, one can either starts fooling around with some of
     the auxiliary executables, or assume all is well with the world and
     start running the script /SCRIPTS/kcarta.sc (to make kcarta.sc an
     executable, you might have to type "chmod +x kcarta.sc" at the
     UNIX prompt, when in the /SCRIPTS subdirectory.  Also, it would
     behoove the user to check that the subdir variable in kcarta.sc
     is set to OUTPUT, and that NO files are present in the OUTPUT
     subdirectory).


  3. Assuming all is well, kcarta.sc will call makeprofile.x and makeinp.x
     to create a profile for *PTHFIL, atmosphere information for *RAD-
     NCE and an input driver file. kcarta.x is then called, with the output
     for this run going to files that are saved in /OUTPUT. This is re-
     peated for a set of test profiles/frequency combinations. Through-
     out this, temporary files inputfile and outfile are created in /SRC,
     as is another temporary file /SCRIPTS/header.head ... these are all
     deleted at the end of the script run.


  4. The user can now check to see if he/she has the same results as the
     files supplied by us in COMPARISON. To do this, run diffemall.sc in
     the SCRIPT subdirectory.

     As mentioned above, what this script file does is runs readkcarta.x
     in the BIN subdirectory through the data files you produced (which
     should be in the OUTPUT subdirectory). The resulting text output
     is in a 2 colummn format, the first column being the wavenumber
     and the second column being the actual data (in this case, a radi-
     ance)



     Then the UNIX "diff" command is used to compare the two rele-
     vant text files (one in the /OUTPUT subdirectory, the other in the
     /COMPARISON directory).  Due to roundoff error, there might be
                                  8


     some slight differences (if the difference are large, then something
     is very wrong!!)


  5. Instead of resaving the data in a text format, MATLAB v 5.0 can be
     run to directly read in the data files.  The driver .m file is called
     readkcarta. It prompts the user for the name of the binary file to
     be read in w/o the .dat extension (e.g. COMPARISON/testa1). The
     user can then


      (a) save ALL data for a particular 25 cm 1 wavenumber region.
          This means all the path/MP/radiance spectra that were saved
          for that region, will be stored in raaData. The data can be plot-
          ted as plot(raFreq,raaData(ii,:)) where ii is the path/MP/radiance
          the user wishes to see

      (b) save ENTIRE data, across the entire frequency range, for one of
          the paths/MPs/layers that are saved .. this goes into raEntire,
          and can then be plotted as plot(raFreq,raEntire)


  6. The results from e.g. OUTPUT/testa1 can then be read in similarly,
     and compared.
4   Settings of parameter files for tests


For the kcarta.sc test run, atmospheres are built up for a downward
looking instrument. The effects of background thermal radiation as well
as effects of solar radiation, are booth turned off.  The atmosphere is
defined between the maximum and minimum AIRS pressure levels. The
surface emissivity is found in emissivty.dat, while the surface tempera-
ture is the temperature of the lowest layer. The settings of the files used
in running the script are as follows
!            **************************************************
!            *                PROGRAM KCARTA : TEMPLATE FILE      *
!            *                ------------------------------      *
!            * compressed data file base for gases 1-28         *
!            * continuum data file calculation for gas 1,22    *
!            * cross section data file base for gases 51-63    *
!            *                                                          *
!            * version 0.4   Sergio DeSouza-Machado   12/20/97   *



                                  9


!            **************************************************
!the FLAG* are to substitute the correct testing parameters (also in
!the radiance file) FLAG1 for freqs, FLAG2 for output filename,
! and FLAG3 for the profile file
!-------------------------------------------------------------------
!this is for FRQNCY -- this section is always altered by makeinp.x,
!with the start/stop frequencies being modified
*FRQNCY
!FLAG1
 1660.0 1670.0
!********************************************************************
!this is for MOLGAS --- all gases in the compressed database included
*MOLGAS
-1
!********************************************************************
!this is for XSCGAS - all gases in xsecdata   base included.
*XSCGAS
'../DATA/General/xsecdata.dat'
-1
!********************************************************************
!this is for PRFILE - this section is changed by makeprofile.x ... it
! gives the name of the new profile (created from the test regression
! profiles 1-40)
*PRFILE
!FLAG3
'../OUTPUT/testprof1'
!********************************************************************
!this is for WEIGHT - all gases are included, with a weight of 1.0
*WEIGHT
1
1 -1 1.0 -1
!********************************************************************
!this is for RADNCE - this is always modified by makeinp.f
! start pressure = 1200, stop pressure = 0.0 ==> downward looking inst
! 2.96 k = temperature of cold dark space
! 303.34 = current surface temp (always reset by makeinp.f)
! 0.0 ==> satellite is nadir looking
! emissivity.dat contains emissivities
*RADNCE
1



                                 10


!FLAG4
  1   1200.0 0.0   2.96    303.34000    0.0 -1.0
  -1 -1.0 -1 -1 1.0 1
'../OUTPUT/emissivity.dat'
!********************************************************************
!this is for Jacobians - no Jacobians!!
!*JACOBN
!1
!1
!********************************************************************
!this is for output - the file name is always changfed by makeinp.f
!output radiance is at 0.0 mb === top of atmosphere
*OUTPUT
!this is the   comment .. same as the file name for now ....
!FLAG2
'../OUTPUT/testa1.dat'
!this is the ouput file name OPFNAM
!FLAG2
'../OUTPUT/testa1.dat'
!this sets iDat
3
1 1
0.0
!********************************************************************
!end of input
*ENDINP
!********************************************************************
   See kcarta.ps for examples and more detailed information on the for-
mat of driver files.

                                 11
