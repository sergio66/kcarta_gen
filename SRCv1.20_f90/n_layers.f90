! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_layers

USE basic_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile
USE n_mr_common

IMPLICIT NONE

CONTAINS

! this file deals with a simple layers
! does everything in terms of LINEAR interps, not SPLINE
! DO NOT USE        Call r_sort_logspl(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)
! DO     USE        Call r_sort_loglinear(raTempP,raTempMR,i1-i2+1,raP,raTempMRX,iNumLevs)

! assumes WV is in "wet air", other gases in dry air VMR (klayers unit 12)
! code
!  (a) reads in input text levels profile, gases in whatever units
!      changes the input units to MIXING RATIOS
!  (b) read in upper level VMR profiles for above gases, as well as for "missing gases"
!      adjusts the Ref VMR to agree with input where the input levels end, and then tapers the
!        adjust ratio to 1, with about 4 AIRS levels
!      adjusts the Ref Temp to agree with input where the input levels end, and then tapers the
!        adjust ratio to 0, with about 4 AIRS levels
!   ** the tacking on (or adding new gases) can either be via
!       - using the reference P,PP,T for 100 layers ==> VMR = PP/P
!       - using one of the six AFGL profiles        ==> VMR comes from the ~50 levels
!  (c) if Earth do the adjustment of dry air mix ratios
!  (d) does the sublev integrations

! Useful References
!   http://www.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_1/index.html
!   http://cimss.ssec.wisc.edu/~paulv/Fortran90/Profile_Utility/profile_units_conversion_2/index.html
! saved off in ../PDF/wisc_notes_*.pdf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! this subroutine reads in a TEST file LEVELS mixing ratio PROFILE
! the text file will have N levs x (G+2) columns of the format

! caStr80   = comment
! numlevs   = number of levels
! Psurf,TSurf,HSurf = surface pressure (mb) and temperature (K) and height (km)
! year lat lon     = year of obs (for co2 MR) and gravity adjust (not yet implemented)
!                     if year < 0, ignore
!                     if year > 0 and kPlanet = 3, then if adding on Gas2, adjust for ppmv
! numgases  = number of gases
! gasIDS    = which gases
! gasUNITS  = eg MR, g/g, RH, same as klayers
!   p1     T1   g1_1  g2_1 .... gG_1
!   p2     T2   g1_2  g2_2 .... gG_2
!   p3     T3   g1_3  g2_3 .... gG_3
!   p4     T4   g1_4  g2_4 .... gG_4

!   pN     TN   g1_N  g2_N .... gG_N

! where p(i) is in mb, T(i) is in Kelvin and N <= kProfLayer*2
! gi_j (j=1--N) are the gas amounts eg dimensionless 0 < VMR < 1, 0 < RH < 100 etc

!        See eg KCARTA/IP_PROFILES/levelsRTP_to_levelstext.m
!             KCARTA/IP_PROFILES/levelsprofile_text1.prf

! see /asl/packages/klayers/Doc/description.txt

! or can read TAPE5, modified TAPE6 from LBLRTM


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE n_layers
