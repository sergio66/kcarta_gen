The files in this subdirectory read the data output from kCARTAv0.4+
Some have the capacity to further manipulate the data by doing convolutions
this is for kCARTA v1.03+

1)readkcarta.f          This is the same as readatmos.m, and can read the 
                         basic outout binary file in two options 
                         (a) read ALL data for a 25 cm-1 block
                         (b) read ENTIRE data for one of the output spectra
                         Output is in the form of a long text file
                         The data can have been output in long(+1) or 
                         short(-1) form
2)readjacob.f            This is the same as readjacob.m, and can read the
                         basic output JACOBIAN file in one option : read
                         ALL data for a 25 cm-1 block

3)makeinp.f              this program is used by kcarta.sc to create a
                         copy of the template file 
                         DATA/Template/kcarta_template, with the required
                         frequency start/stop points, surface temperature
                         and output file name

4)makeprofile.f          this program is used by kcarta.sc to create a
                         profile file from the supplied regression
                         profiles in DATA/TestProf

5)compdatabase.f         is a program that the user should run (using
                         comp.sc) whenever there is an update to the 
                         files in the compressed database. This program
                         produces parameter file comp.param
 to use this code, go to the subdir containing the database
 eg cd /salsify/production3/motteler/sergio/kcarta/h2o.ieee-be
 do ls -1 *.dat > & /salsify/users/sergio/KCARTA/UTILITY/compdatabase
 then run this program
 you will have to do it separately for the water database and the "rest
 of gases" database; combine the two and then separate out the minor or
 cross section gases. eg put the results into 2 files, 
 '../DATA/General/xsec107.param' and '../DATA/General/comp107.param'

convolve.param is a parameter file used by some of the above programs

The rest are utility files that are linked in by some of the above files. 
/UTILv0.97 contains files for our upcoming heralded hoped for release,
which will allow kCARTA to have an extended database (~255-2805)
