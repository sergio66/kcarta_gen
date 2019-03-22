! ************************************************************************
! ******************* THESE ARE THE DISORT ROUTINES **********************
! ************************************************************************
! ccccccc cp BDREF.f /salsify/users/sergio/KCARTA/SRCv1.08/
! RCS version control information:
! $Header: BDREF.f,v 2.1 2000/03/27 21:40:51 laszlo Exp $
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    REAL FUNCTION  BDREF( WVNMLO, WVNMHI, MU, MUP, DPHI )

    include '../INCLUDE/scatterparam.f90'

!      Supplies surface bi-directional reflectivity.


!      This is only a "stub" version. The user must replace this
!      by his/her own BDREF function.


!      NOTE 1: Bidirectional reflectivity in DISORT is defined
!              by Eq. 39 in STWL.
!      NOTE 2: Both MU and MU0 (cosines of reflection and incidence
!              angles) are positive.

!  INPUT:

!    WVNMLO : Lower wavenumber (inv cm) of spectral interval

!    WVNMHI : Upper wavenumber (inv cm) of spectral interval

!    MU     : Cosine of angle of reflection (positive)

!    MUP    : Cosine of angle of incidence (positive)

!    DPHI   : Difference of azimuth angles of incidence and reflection
!                (radians)


!   Called by- DREF, SURFAC

! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

    REAL ::      DPHI, MU, MUP, WVNMHI, WVNMLO
!     ..
!     .. External Subroutines ..

    EXTERNAL  ERRMSG
!     ..

    REAL :: temp

!      WRITE ( *, '(//,7(1X,A,/))' )
!     &  'To use a bidirectionally reflecting lower boundary you must',
!     &  'replace file BDREF.f with your own file. In that file, you ',
!     &  'should supply the bidirectional reflectivity, as a function ',
!     &  'of the cosine of angle of reflection, the cosine of angle ',
!     &  'of incidence, and the difference of azimuth angles of ',
!     &  'incidence and reflection. See DISORT.doc for more information',
!     &  'and subroutine BDREF in file DISOTEST.f for an example.'
!      CALL ERRMSG( 'BDREF--Please supply a surface BDRF model', .TRUE. )
!      BDREF = 0.0

! this is just 1 = e + b

    BDREF = DISORTraBiDirRefl(DISORTiI)

    RETURN
    END FUNCTION 

!************************************************************************
    SUBROUTINE  GETMOM( IPHAS, GG, NMOM, PMOM )

!        Calculate phase function Legendre expansion coefficients
!        in various special cases


!       INPUT: IPHAS   Phase function options
!                      1 : Isotropic
!                      2 : Rayleigh
!                      3 : Henyey-Greenstein with asymmetry factor GG
!                      4 : Haze L as specified by Garcia/Siewert
!                      5 : Cloud C.1 as specified by Garcia/Siewert

!              GG      Asymmetry factor for Henyey-Greenstein case

!              NMOM    Index of highest Legendre coefficient needed
!                        ( = number of streams 'NSTR'  chosen
!                         for the discrete ordinate method)

!      OUTPUT: PMOM(K)  Legendre expansion coefficients (K=0 to NMOM)
!                         (be sure to dimension '0:maxval' in calling
!                          program)

!      Reference:  Garcia, R. and C. Siewert, 1985: Benchmark Results
!                     in Radiative Transfer, Transp. Theory and Stat.
!                     Physics 14, 437-484, Tables 10 And 17
! ------------------------------------------------------------------

!     .. Scalar Arguments ..

    INTEGER ::   IPHAS, NMOM
    REAL ::      GG
!     ..
!     .. Array Arguments ..

    REAL ::      PMOM( 0:NMOM )
!     ..
!     .. Local Scalars ..

    INTEGER ::   K
!     ..
!     .. Local Arrays ..

    REAL ::      CLDMOM( 299 ), HAZELM( 82 )
!     ..
!     .. External Subroutines ..

    EXTERNAL  ERRMSG
!     ..
!     .. Intrinsic Functions ..

    INTRINSIC MIN
!     ..

    DATA HAZELM /  2.41260, 3.23047, 3.37296, 3.23150, 2.89350, &
    &                2.49594, 2.11361, 1.74812, 1.44692, 1.17714, &
    &                0.96643, 0.78237, 0.64114, 0.51966, 0.42563, &
    &                0.34688, 0.28351, 0.23317, 0.18963, 0.15788, &
    &                0.12739, 0.10762, 0.08597, 0.07381, 0.05828, &
    &                0.05089, 0.03971, 0.03524, 0.02720, 0.02451, &
    &                0.01874, 0.01711, 0.01298, 0.01198, 0.00904, &
    &                0.00841, 0.00634, 0.00592, 0.00446, 0.00418, &
    &                0.00316, 0.00296, 0.00225, 0.00210, 0.00160, &
    &                0.00150, 0.00115, 0.00107, 0.00082, 0.00077, &
    &                0.00059, 0.00055, 0.00043, 0.00040, 0.00031, &
    &                0.00029, 0.00023, 0.00021, 0.00017, 0.00015, &
    &                0.00012, 0.00011, 0.00009, 0.00008, 0.00006, &
    &                0.00006, 0.00005, 0.00004, 0.00004, 0.00003, &
    &                0.00003, 3*0.00002, 8*0.00001 /

    DATA  ( CLDMOM(K), K = 1, 159 ) / &
    &   2.544,  3.883,  4.568,  5.235,  5.887,  6.457,  7.177,  7.859, &
    &   8.494,  9.286,  9.856, 10.615, 11.229, 11.851, 12.503, 13.058, &
    &  13.626, 14.209, 14.660, 15.231, 15.641, 16.126, 16.539, 16.934, &
    &  17.325, 17.673, 17.999, 18.329, 18.588, 18.885, 19.103, 19.345, &
    &  19.537, 19.721, 19.884, 20.024, 20.145, 20.251, 20.330, 20.401, &
    &  20.444, 20.477, 20.489, 20.483, 20.467, 20.427, 20.382, 20.310, &
    &  20.236, 20.136, 20.036, 19.909, 19.785, 19.632, 19.486, 19.311, &
    &  19.145, 18.949, 18.764, 18.551, 18.348, 18.119, 17.901, 17.659, &
    &  17.428, 17.174, 16.931, 16.668, 16.415, 16.144, 15.883, 15.606, &
    &  15.338, 15.058, 14.784, 14.501, 14.225, 13.941, 13.662, 13.378, &
    &  13.098, 12.816, 12.536, 12.257, 11.978, 11.703, 11.427, 11.156, &
    &  10.884, 10.618, 10.350, 10.090,  9.827,  9.574,  9.318,  9.072, &
    &   8.822, 8.584, 8.340, 8.110, 7.874, 7.652, 7.424, 7.211, 6.990, &
    &   6.785, 6.573, 6.377, 6.173, 5.986, 5.790, 5.612, 5.424, 5.255, &
    &   5.075, 4.915, 4.744, 4.592, 4.429, 4.285, 4.130, 3.994, 3.847, &
    &   3.719, 3.580, 3.459, 3.327, 3.214, 3.090, 2.983, 2.866, 2.766, &
    &   2.656, 2.562, 2.459, 2.372, 2.274, 2.193, 2.102, 2.025, 1.940, &
    &   1.869, 1.790, 1.723, 1.649, 1.588, 1.518, 1.461, 1.397, 1.344, &
    &   1.284, 1.235, 1.179, 1.134, 1.082, 1.040, 0.992, 0.954, 0.909 /
    DATA  ( CLDMOM(K), K = 160, 299 ) / &
    &   0.873, 0.832, 0.799, 0.762, 0.731, 0.696, 0.668, 0.636, 0.610, &
    &   0.581, 0.557, 0.530, 0.508, 0.483, 0.463, 0.440, 0.422, 0.401, &
    &   0.384, 0.364, 0.349, 0.331, 0.317, 0.301, 0.288, 0.273, 0.262, &
    &   0.248, 0.238, 0.225, 0.215, 0.204, 0.195, 0.185, 0.177, 0.167, &
    &   0.160, 0.151, 0.145, 0.137, 0.131, 0.124, 0.118, 0.112, 0.107, &
    &   0.101, 0.097, 0.091, 0.087, 0.082, 0.079, 0.074, 0.071, 0.067, &
    &   0.064, 0.060, 0.057, 0.054, 0.052, 0.049, 0.047, 0.044, 0.042, &
    &   0.039, 0.038, 0.035, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026, &
    &   0.024, 0.023, 0.022, 0.021, 0.020, 0.018, 0.018, 0.017, 0.016, &
    &   0.015, 0.014, 0.013, 0.013, 0.012, 0.011, 0.011, 0.010, 0.009, &
    &   0.009, 3*0.008, 2*0.007, 3*0.006, 4*0.005, 4*0.004, 6*0.003, &
    &   9*0.002, 18*0.001 /


    IF ( IPHAS < 1 .OR. IPHAS > 5 ) &
    CALL ERRMSG( 'GETMOM--bad input variable IPHAS', .TRUE. )

    IF ( IPHAS == 3 .AND. (GG <= -1.0 .OR. GG >= 1.0) ) &
    CALL ERRMSG( 'GETMOM--bad input variable GG', .TRUE. )

    IF ( NMOM < 2 ) &
    CALL ERRMSG( 'GETMOM--bad input variable NMOM', .TRUE. )


    PMOM(0) = 1.0
    DO  10  K = 1, NMOM
        PMOM(K) = 0.0
    10 END DO


    IF ( IPHAS == 2 )  THEN
    !                                       ** Rayleigh phase function
        PMOM(2) = 0.1

    ELSE IF ( IPHAS == 3 ) THEN
    !                                       ** Henyey-Greenstein phase fcn
        DO  20  K = 1, NMOM
            PMOM(K) = GG**K
        20 END DO

    ELSE IF ( IPHAS == 4 ) THEN
    !                                        ** Haze-L phase function
        DO  30  K = 1, MIN(82,NMOM)
            PMOM(K) = HAZELM(K) / ( 2*K+1 )
        30 END DO

    ELSE IF ( IPHAS == 5 ) THEN
    !                                        ** Cloud C.1 phase function
        DO  40  K = 1, MIN(298,NMOM)
            PMOM(K) = CLDMOM(K) / ( 2*K+1 )
        40 END DO

    END IF

    RETURN
    END SUBROUTINE 

! ************************************************************************
