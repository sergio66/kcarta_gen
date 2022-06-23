! Copyright 2015
! University of Maryland Baltimore County
! All Rights Resered

! this file deals with solar insolation
! https://data.giss.nasa.gov/modelE/ar5plots/solar.html
! Dr Gary Russell, https://www.giss.nasa.gov/staff/grussell.html

MODULE solar_insolation

! fixcon SOLAR_INSOLATION/SRLOCATSergio_fixdaylatlon.FOR solar_insolation.f90
!
! this is same as SRLOCAT.FOR except I have changed output format statements and replaced
! "always daytime" and "always nighttime" with 0:00-24:00 and 0:00 to 0.00
! so eg do "solinsSergio_looplat.x > junk" to get daily 2012 solar insolation at all lat
! 13 columns  eg
!         If (JDATE >= 1)  Write (6,941) RLAT,RLON,JYEAR,MONTH,JDATE,
!     *                   IDAWNH,IDAWNM, IDUSKH,IDUSKM, SRINC,COSZS,
!     *                   COSZT,SUNDIS
!
! rlat rlon 2012 01 01     7 18    16 41     160.34   0.360  0.98  1.03
!
! see plot_solar_insolation.m
!

IMPLICIT NONE

CONTAINS

!************************************************************************
      SUBROUTINE GetSolarInsolation(JYEAR,JMON,JDATE,snglLAT,snglLON,snglSUNDIS)

      IMPLICIT NONE
      include '../INCLUDE/TempF90/kcartaparam.f90'

! input
      INTEGER JYEAR,JMON,JDATE   !! input
      REAL    SnglLat,SnglLon
! output
      REAL    snglSUNDIS


      REAL*8    RLat,RLon,ECCEN, OBLIQ, OMEGVP, DATE, DAY

!      Implicit Real*8 (A-H,O-Z)
      REAL*8 TWOPI,EDAYzY,RSUNd,REFRAC
      REAL*8 SIND, COSD, SUNDIS, SUNLON, SUNLAT, EQTIME, COSZT, COSZS
      REAL*8  RSmEzM, DUSK, DAWN, SRINC
      
      INTEGER IDAWNH, IDUSKH, IDAWNM, IDUSKM, YEAR,MONTH
      INTEGER DAY2MO(12),ITZONE
      Parameter (TWOPI=6.283185307179586477d0, EDAYzY=365.2425d0, &
                 RSUNd =.267d0, &  !  mean radius of Sun in degrees 
                 REFRAC=.583d0)  !  effective Sun disk increase
      Integer*4 DAYzMO(12), DATMAX
      Character ARG*80, CTZONE(-12:12)*33
!****
      Data DAY2MO /31,28,31, 30,31,30, 31,31,30, 31,30,31/
      Data CTZONE /'Yankee Standard Time', &
                    'Phoenix Island Standard Time', &
                    'Hawaiian Standard Time', &
                    'Alaska Standard Time', &
                    'Pacific Standard Time', &
                    'Mountain Standard Time', &
                    'Central Standard Time', &
                    'Eastern Standard Time', &
                    'Atlantic Standard Time', &
                    'Brazilian Standard Time', &
                    'Fernando de Noronha Standard Time', &
                    'Jan Mayen Standard Time', &
                    'Greenwich Mean Time', &
                    'Central Europe Time', &
                    'Eastern Europe Time', &
                    'Moscow Standard Time', &
                    'Afghan Standard Time', &
                    'Indian Standard Time', &
                    'Bangladesh Standard Time', &
                    'Vietnam Standard Time', &
                    'Australian Western Standard Time', &
                    'Japanese Standard Time', &
                    'Australian Eastern Standard Time', &
                    'New Caledonia Standard Time', &
                    'Fiji Standard Time'/
!****
!      NARGS = IArgC()
!      If (NARGS <= 0)  GoTo 800
!      RLAT   =  40.78
!      RLAT   = 0.0
!      RLON   = -73.97

! testing nick nalli emiss
!      JYEAR  = 2019
!      JMON = 7
!      JDATE  = 2

!      JYEAR  = 2019
!      JMON = 1
!      JDATE  = 1


      RLAT = snglLat*1.0D0
      RLON = snglLon*1.0D0

!**** Determine time zone
  200 ITZONE = Nint (RLON/15)

!**** Determine orbital parameters
      YEAR = JYEAR
      Call ORBPAR (YEAR,ECCEN,OBLIQ,OMEGVP)

      MONTH=JMON

!**** Loop over days
        DATMAX = DAYzMO(MONTH)
        If (MONTH == 2 .and. QLEAPY(JYEAR))  DATMAX = 29

        DATE = JDATE-1 + .5 - RLON/360
        Call YMDtoD (JYEAR,MONTH,DATE,   DAY)
        Call ORBIT  (ECCEN,OBLIQ,OMEGVP, DAY, &
                    SIND,COSD,SUNDIS,SUNLON,SUNLAT,EQTIME)
        Call COSZIJ (RLAT,SIND,COSD, COSZT,COSZS)
        RSmEzM = (REFRAC + RSUNd/SUNDIS) * TWOPI/360
        Call SUNSET (RLAT,SIND,COSD,RSmEzM, DUSK)
        If (DUSK >=  999999)  GoTo 500
        If (DUSK <= -999999)  GoTo 600
!****
!**** Daylight and nightime at this location on this day
!****
        DAWN   = (-DUSK-EQTIME)*24/TWOPI + 12 - RLON/15 + ITZONE
        DUSK   = ( DUSK-EQTIME)*24/TWOPI + 12 - RLON/15 + ITZONE
        IDAWNH = DAWN
        IDUSKH = DUSK
        IDAWNM = Nint ((DAWN-IDAWNH)*60)
        IDUSKM = Nint ((DUSK-IDUSKH)*60)
        if (IDAWNM .LT. 0) IDAWNM = 0
        if (IDUSKM .LT. 0) IDUSKM = 0
        SRINC  = 1367*COSZT / SUNDIS**2
        If (JDATE >= 1)  Write (kStdWarn,941) RLAT,RLON,JYEAR,MONTH,JDATE, &
                         IDAWNH,IDAWNM, IDUSKH,IDUSKM, SRINC,COSZS, &
                         COSZT,SUNDIS
        GoTo 700
!****
!****   Daylight at all times at this location on this day
!****
  500   SRINC = 1367*COSZT / SUNDIS**2
        If (JDATE >= 1)  Write (kStdWarn,941) RLAT,RLON,JYEAR,MONTH,JDATE, &
                         0,0,24,0,SRINC,COSZS,COSZT,SUNDIS
        GoTo 700
!****
!**** Nightime at all times at this location on this day
!****
  600   SRINC = 1367*COSZT / SUNDIS**2
        If (JDATE >= 1)  Write (kStdWarn,941) RLAT,RLON,JYEAR,MONTH,JDATE, &
                         0,0,0,0,SRINC,COSZS,COSZT,SUNDIS
  700   Continue
  699  Continue

       snglSUNDIS = real(SUNDIS)   !!! the output I need!!!!!
       GoTo 999
!****

  810 Write(kStdWarn,*) 'Unable to decipher latitude: ',Trim(ARG)
      Stop 810
  811 Write(kStdWarn,*) 'Specified latitude is out or range: ',RLAT
      Stop 811
  812 Write(kStdWarn,*) 'Unable to decipher longitude: ',Trim(ARG)
      Stop 812
  813 Write(kStdWarn,*) 'Specified longitude is out or range: ',RLON
      Stop 813
  814 Write(kStdWarn,*) 'Unable to decipher year: ',Trim(ARG)
      Stop 814
  816 Write(kStdWarn,*) 'Unable to decipher month: ',Trim(ARG)
      Stop 816
  817 Write(kStdWarn,*) 'Specified month is out or range: ',MONTH
      Stop 817
!****
  900 Format ('% Content-type: text/plain' )
  930 Format ('% Daily Insolation Parameters assuming 1367 W/m2 flux' / &
              '% Latitude:',F8.3,5X,'Longitude:',F8.3 / &
      '% Time Zone: ',A / &
      '% (Longitudes',F8.1,'  to',F7.1,')                                  '/ &
      '%                                           SRINC     COSZS     COSZT  SUNDIS' / &
      '%                                                    Sunlight    Time   ' / &
      '%                                                    Weighted   Weighted' / &
      '%                                          Daily     Cosine     Cosine ' / &
      '%                                          Average      of        of    ' / &
      '%                                          Sunlight   Zenith    Zenith  ' / &
      '%      RLAT    RLON    Date   Sunrise   Sunset    (W/m2)    Angle    Angle  ' / &
      '%      ----    ----   -------   ------    ------    -----  ------  ------')
!  941 Format (I4,2('/',I2.2), 2(I6,':',I2.2), F11.2, F8.3)
!  942 Format (        I10.2 , 2(I6,':',I2.2), F11.2, F8.3)
  941 Format (F8.3,' ',F8.3,' ',I4,2(' ',I2.2), 2(I6,' ',I2.2), F11.2, 3(' ',F8.3))
  942 Format (        I10.2 , 2(I6,' ',I2.2), F11.2, F8.3)
  951 Format (I4,2('/',I2.2),'   Always Day Light', F11.2, F8.3)
  952 Format (        I10.2 ,'   Always Day Light', F11.2, F8.3)
  961 Format (I4,2('/',I2.2),'   Always Nightime ', F11.2, F8.3)
  962 Format (        I10.2 ,'   Always Nightime ', F11.2, F8.3)
  999 End

!************************************************************************
      Subroutine ORBPAR (YEAR, ECCEN,OBLIQ,OMEGVP)
!****
!**** ORBPAR calculates the three orbital parameters as a function of
!**** YEAR.  The source of these calculations is: Andre L. Berger,
!**** 1978, "Long-Term Variations of Daily Insolation and Quaternary
!**** Climatic Changes", JAS, v.35, p.2362.  Also useful is: Andre L.
!**** Berger, May 1978, "A Simple Algorithm to Compute Long Term
!**** Variations of Daily Insolation", published by Institut
!**** D'Astronomie de Geophysique, Universite Catholique de Louvain,
!**** Louvain-la Neuve, No. 18.
!****
!**** Tables and equations refer to the first reference (JAS).  The
!**** corresponding table or equation in the second reference is
!**** enclosed in parentheses.  The coefficients used in this
!**** subroutine are slightly more precise than those used in either
!**** of the references.  The generated orbital parameters are precise
!**** within plus or minus 1000000 years from present.
!****
!**** Input:  YEAR   = years A.D. are positive, B.C. are negative
!**** Output: ECCEN  = eccentricity of orbital ellipse
!****         OBLIQ  = latitude of Tropic of Cancer in radians
!****         OMEGVP = longitude of perihelion =
!****                = spatial angle from vernal equinox to perihelion
!****                  in radians with sun as angle vertex
!****

!      Implicit Real*8 (A-H,O-Z)
      INTEGER YEAR,I
      REAL*8 TWOPI, PIz180, ECCEN, OBLIQ, OMEGVP, YEAR1950, SUMC, ARG, YM1950
      REAL*8 OBLIQD, ESINPI, ECOSPI, PIE, FSINFD, PSI
      
      Parameter (TWOPI=6.283185307179586477d0, PIz180=TWOPI/360d0)
      Real*8 TABLE1(3,47),TABLE4(3,19),TABLE5(3,78)
!**** Table 1 (2).  Obliquity relative to mean ecliptic of date: OBLIQD
      Data TABLE1/-2462.2214466d0, 31.609974d0, 251.9025d0, &
                   -857.3232075d0, 32.620504d0, 280.8325d0, &
                   -629.3231835d0, 24.172203d0, 128.3057d0, &
                   -414.2804924d0, 31.983787d0, 292.7252d0, &
                   -311.7632587d0, 44.828336d0,  15.3747d0, &
                    308.9408604d0, 30.973257d0, 263.7951d0, &
                   -162.5533601d0, 43.668246d0, 308.4258d0, &
                   -116.1077911d0, 32.246691d0, 240.0099d0, &
                    101.1189923d0, 30.599444d0, 222.9725d0, &
                    -67.6856209d0, 42.681324d0, 268.7809d0, &
                     24.9079067d0, 43.836462d0, 316.7998d0, &
                     22.5811241d0, 47.439436d0, 319.6024d0, &
                    -21.1648355d0, 63.219948d0, 143.8050d0, &
                    -15.6549876d0, 64.230478d0, 172.7351d0, &
                     15.3936813d0,  1.010530d0,  28.9300d0, &
                     14.6660938d0,  7.437771d0, 123.5968d0, &
                    -11.7273029d0, 55.782177d0,  20.2082d0, &
                     10.2742696d0,   .373813d0,  40.8226d0, &
                      6.4914588d0, 13.218362d0, 123.4722d0, &
                      5.8539148d0, 62.583231d0, 155.6977d0, &
                     -5.4872205d0, 63.593761d0, 184.6277d0, &
                     -5.4290191d0, 76.438310d0, 267.2772d0, &
                      5.1609570d0, 45.815258d0,  55.0196d0, &
                      5.0786314d0,  8.448301d0, 152.5268d0, &
                     -4.0735782d0, 56.792707d0,  49.1382d0, &
                      3.7227167d0, 49.747842d0, 204.6609d0, &
                      3.3971932d0, 12.058272d0,  56.5233d0, &
                     -2.8347004d0, 75.278220d0, 200.3284d0, &
                     -2.6550721d0, 65.241008d0, 201.6651d0, &
                     -2.5717867d0, 64.604291d0, 213.5577d0, &
                     -2.4712188d0,  1.647247d0,  17.0374d0, &
                      2.4625410d0,  7.811584d0, 164.4194d0, &
                      2.2464112d0, 12.207832d0,  94.5422d0, &
                     -2.0755511d0, 63.856665d0, 131.9124d0, &
                     -1.9713669d0, 56.155990d0,  61.0309d0, &
                     -1.8813061d0, 77.448840d0, 296.2073d0, &
                     -1.8468785d0,  6.801054d0, 135.4894d0, &
                      1.8186742d0, 62.209418d0, 114.8750d0, &
                      1.7601888d0, 20.656133d0, 247.0691d0, &
                     -1.5428851d0, 48.344406d0, 256.6114d0, &
                      1.4738838d0, 55.145460d0,  32.1008d0, &
                     -1.4593669d0, 69.000539d0, 143.6804d0, &
                      1.4192259d0, 11.071350d0,  16.8784d0, &
                     -1.1818980d0, 74.291298d0, 160.6835d0, &
                      1.1756474d0, 11.047742d0,  27.5932d0, &
                     -1.1316126d0,  0.636717d0, 348.1074d0, &
                      1.0896928d0, 12.844549d0,  82.6496d0/
!**** Table 4 (1).  Fundamental elements of the ecliptic: ECCEN sin(pi)
      Data TABLE4/ .01860798d0,  4.207205d0,  28.620089d0, &
                   .01627522d0,  7.346091d0, 193.788772d0, &
                  -.01300660d0, 17.857263d0, 308.307024d0, &
                   .00988829d0, 17.220546d0, 320.199637d0, &
                  -.00336700d0, 16.846733d0, 279.376984d0, &
                   .00333077d0,  5.199079d0,  87.195000d0, &
                  -.00235400d0, 18.231076d0, 349.129677d0, &
                   .00140015d0, 26.216758d0, 128.443387d0, &
                   .00100700d0,  6.359169d0, 154.143880d0, &
                   .00085700d0, 16.210016d0, 291.269597d0, &
                   .00064990d0,  3.065181d0, 114.860583d0, &
                   .00059900d0, 16.583829d0, 332.092251d0, &
                   .00037800d0, 18.493980d0, 296.414411d0, &
                  -.00033700d0,  6.190953d0, 145.769910d0, &
                   .00027600d0, 18.867793d0, 337.237063d0, &
                   .00018200d0, 17.425567d0, 152.092288d0, &
                  -.00017400d0,  6.186001d0, 126.839891d0, &
                  -.00012400d0, 18.417441d0, 210.667199d0, &
                   .00001250d0,  0.667863d0,  72.108838d0/
!**** Table 5 (3).  General precession in longitude: psi
      Data TABLE5/ 7391.0225890d0, 31.609974d0, 251.9025d0, &
                   2555.1526947d0, 32.620504d0, 280.8325d0, &
                   2022.7629188d0, 24.172203d0, 128.3057d0, &
                  -1973.6517951d0,  0.636717d0, 348.1074d0, &
                   1240.2321818d0, 31.983787d0, 292.7252d0, &
                    953.8679112d0,  3.138886d0, 165.1686d0, &
                   -931.7537108d0, 30.973257d0, 263.7951d0, &
                    872.3795383d0, 44.828336d0,  15.3747d0, &
                    606.3544732d0,  0.991874d0,  58.5749d0, &
                   -496.0274038d0,  0.373813d0,  40.8226d0, &
                    456.9608039d0, 43.668246d0, 308.4258d0, &
                    346.9462320d0, 32.246691d0, 240.0099d0, &
                   -305.8412902d0, 30.599444d0, 222.9725d0, &
                    249.6173246d0,  2.147012d0, 106.5937d0, &
                   -199.1027200d0, 10.511172d0, 114.5182d0, &
                    191.0560889d0, 42.681324d0, 268.7809d0, &
                   -175.2936572d0, 13.650058d0, 279.6869d0, &
                    165.9068833d0,  0.986922d0,  39.6448d0, &
                    161.1285917d0,  9.874455d0, 126.4108d0, &
                    139.7878093d0, 13.013341d0, 291.5795d0, &
                   -133.5228399d0,  0.262904d0, 307.2848d0, &
                    117.0673811d0,  0.004952d0,  18.9300d0, &
                    104.6907281d0,  1.142024d0, 273.7596d0, &
                     95.3227476d0, 63.219948d0, 143.8050d0, &
                     86.7824524d0,  0.205021d0, 191.8927d0, &
                     86.0857729d0,  2.151964d0, 125.5237d0, &
                     70.5893698d0, 64.230478d0, 172.7351d0, &
                    -69.9719343d0, 43.836462d0, 316.7998d0, &
                    -62.5817473d0, 47.439436d0, 319.6024d0, &
                     61.5450059d0,  1.384343d0,  69.7526d0, &
                    -57.9364011d0,  7.437771d0, 123.5968d0, &
                     57.1899832d0, 18.829299d0, 217.6432d0, &
                    -57.0236109d0,  9.500642d0,  85.5882d0, &
                    -54.2119253d0,  0.431696d0, 156.2147d0, &
                     53.2834147d0,  1.160090d0,  66.9489d0, &
                     52.1223575d0, 55.782177d0,  20.2082d0, &
                    -49.0059908d0, 12.639528d0, 250.7568d0, &
                    -48.3118757d0,  1.155138d0,  48.0188d0, &
                    -45.4191685d0,  0.168216d0,   8.3739d0, &
                    -42.2357920d0,  1.647247d0,  17.0374d0, &
                    -34.7971099d0, 10.884985d0, 155.3409d0, &
                     34.4623613d0,  5.610937d0,  94.1709d0, &
                    -33.8356643d0, 12.658184d0, 221.1120d0, &
                     33.6689362d0,  1.010530d0,  28.9300d0, &
                    -31.2521586d0,  1.983748d0, 117.1498d0, &
                    -30.8798701d0, 14.023871d0, 320.5095d0, &
                     28.4640769d0,  0.560178d0, 262.3602d0, &
                    -27.1960802d0,  1.273434d0, 336.2148d0, &
                     27.0860736d0, 12.021467d0, 233.0046d0, &
                    -26.3437456d0, 62.583231d0, 155.6977d0, &
                     24.7253740d0, 63.593761d0, 184.6277d0, &
                     24.6732126d0, 76.438310d0, 267.2772d0, &
                     24.4272733d0,  4.280910d0,  78.9281d0, &
                     24.0127327d0, 13.218362d0, 123.4722d0, &
                     21.7150294d0, 17.818769d0, 188.7132d0, &
                    -21.5375347d0,  8.359495d0, 180.1364d0, &
                     18.1148363d0, 56.792707d0,  49.1382d0, &
                    -16.9603104d0,  8.448301d0, 152.5268d0, &
                    -16.1765215d0,  1.978796d0,  98.2198d0, &
                     15.5567653d0,  8.863925d0,  97.4808d0, &
                     15.4846529d0,  0.186365d0, 221.5376d0, &
                     15.2150632d0,  8.996212d0, 168.2438d0, &
                     14.5047426d0,  6.771027d0, 161.1199d0, &
                    -14.3873316d0, 45.815258d0,  55.0196d0, &
                     13.1351419d0, 12.002811d0, 262.6495d0, &
                     12.8776311d0, 75.278220d0, 200.3284d0, &
                     11.9867234d0, 65.241008d0, 201.6651d0, &
                     11.9385578d0, 18.870667d0, 294.6547d0, &
                     11.7030822d0, 22.009553d0,  99.8233d0, &
                     11.6018181d0, 64.604291d0, 213.5577d0, &
                    -11.2617293d0, 11.498094d0, 154.1631d0, &
                    -10.4664199d0,  0.578834d0, 232.7153d0, &
                     10.4333970d0,  9.237738d0, 138.3034d0, &
                    -10.2377466d0, 49.747842d0, 204.6609d0, &
                     10.1934446d0,  2.147012d0, 106.5938d0, &
                    -10.1280191d0,  1.196895d0, 250.4676d0, &
                     10.0289441d0,  2.133898d0, 332.3345d0, &
                    -10.0034259d0,  0.173168d0,  27.3039d0/
!****
      YM1950 = YEAR-1950
!****
!**** Obliquity from Table 1 (2):
!****   OBLIQ# = 23.320556 (degrees)             Equation 5.5 (15)
!****   OBLIQD = OBLIQ# + sum[A cos(ft+delta)]   Equation 1 (5)
!****
      SUMC = 0
      Do I=1,47
        ARG    = PIz180*(YM1950*TABLE1(2,I)/3600+TABLE1(3,I))
        SUMC   = SUMC + TABLE1(1,I)*Cos(ARG)
      END DO
      OBLIQD = 23.320556d0 + SUMC/3600
      OBLIQ  = OBLIQD*PIz180
!****
!**** Eccentricity from Table 4 (1):
!****   ECCEN sin(pi) = sum[M sin(gt+beta)]           Equation 4 (1)
!****   ECCEN cos(pi) = sum[M cos(gt+beta)]           Equation 4 (1)
!****   ECCEN = ECCEN sqrt[sin(pi)^2 + cos(pi)^2]
!****
      ESINPI = 0
      ECOSPI = 0
      Do I=1,19
        ARG    = PIz180*(YM1950*TABLE4(2,I)/3600+TABLE4(3,I))
        ESINPI = ESINPI + TABLE4(1,I)*Sin(ARG)
        ECOSPI = ECOSPI + TABLE4(1,I)*Cos(ARG)
      END DO
      ECCEN  = Sqrt (ESINPI*ESINPI+ECOSPI*ECOSPI)
!****
!**** Perihelion from Equation 4,6,7 (9) and Table 4,5 (1,3):
!****   PSI# = 50.439273 (seconds of degree)         Equation 7.5 (16)
!****   ZETA =  3.392506 (degrees)                   Equation 7.5 (17)
!****   PSI = PSI# t + ZETA + sum[F sin(ft+delta)]   Equation 7 (9)
!****   PIE = atan[ECCEN sin(pi) / ECCEN cos(pi)]
!****   OMEGVP = PIE + PSI + 3.14159                 Equation 6 (4.5)
!****
      PIE    = ATan2(ESINPI,ECOSPI)
      FSINFD = 0
      DO I=1,78
        ARG    = PIz180*(YM1950*TABLE5(2,I)/3600+TABLE5(3,I))
        FSINFD = FSINFD + TABLE5(1,I)*Sin(ARG)
      END DO
      PSI    = PIz180*(3.392506d0+(YM1950*50.439273d0+FSINFD)/3600)
      OMEGVP = Modulo (PIE+PSI+.5*TWOPI, TWOPI)
!****
      Return
      End

!************************************************************************
      Subroutine COSZIJ (RLAT,SIND,COSD, COSZT,COSZS)
!****
!**** COSZIJ calculates the daily average cosine of the zenith angle
!**** weighted by time and weighted by sunlight.
!****
!**** Input: RLAT = latitude (degrees)
!****   SIND,COSD = sine and cosine of the declination angle
!****
!**** Output: COSZT = sum(cosZ*dT) / sum(dT)
!****         COSZS = sum(cosZ*cosZ*dT) / sum(cosZ*dT)
!****
!**** Intern: DAWN = time of DAWN (temporal radians) at mean local time
!****         DUSK = time of DUSK (temporal radians) at mean local time
!****
!      Implicit Real*8 (A-H,O-Z)
      REAL*8 TWOPI,RLAT,SIND,COSD, COSZT,COSZS
      REAL*8 SINJ, COSJ, SJSD, CJSD, CDUSK, DUSK, CJCD, SDUSK, DAWN, SDAWN, S2DAWN
      REAL*8 ECOSZ, QCOSZ, S2DUSK

      Parameter (TWOPI=6.283185307179586477d0)
!****
      SINJ = Sin(TWOPI*RLAT/360)
      COSJ = Cos(TWOPI*RLAT/360)
      SJSD = SINJ*SIND
      CJCD = COSJ*COSD
      If (SJSD+CJCD <= 0)  GoTo 20
      If (SJSD-CJCD >= 0)  GoTo 10
!**** Compute DAWN and DUSK (at local time) and their sines
      CDUSK = - SJSD/CJCD
      DUSK  = ACos(CDUSK)
      SDUSK = SQRT(CJCD*CJCD-SJSD*SJSD) / CJCD
      S2DUSK= 2*SDUSK*CDUSK
      DAWN  = - DUSK
      SDAWN = - SDUSK
      S2DAWN= - S2DUSK
!**** Nightime at initial and final times with daylight in between
      ECOSZ = SJSD*(DUSK-DAWN) + CJCD*(SDUSK-SDAWN)
      QCOSZ = SJSD*ECOSZ + CJCD*(SJSD*(SDUSK-SDAWN) + &
              .5*CJCD*(DUSK-DAWN + .5*(S2DUSK-S2DAWN)))
      COSZT = ECOSZ/TWOPI
      COSZS = QCOSZ/ECOSZ
      Return
!**** Constant daylight at this latitude
   10 DAWN  = -999999
      DUSK  = -999999
      ECOSZ = SJSD*TWOPI
      QCOSZ = SJSD*ECOSZ + .5*CJCD*CJCD*TWOPI
      COSZT = SJSD  !  = ECOSZ/TWOPI
      COSZS = QCOSZ/ECOSZ
      Return
!**** Constant nightime at this latitude
   20 DAWN  = 999999
      DUSK  = 999999
      COSZT = 0
      COSZS = 0
      Return
      End

!************************************************************************
      Subroutine SUNSET (RLAT,SIND,COSD,RSmEzM, DUSK)
!****
!**** Input: RLAT = latitude (degrees)
!****   SIND,COSD = sine and cosine of the declination angle
!****      RSmEzM = (Sun Radius - Earth Radius) / (distance to Sun)
!****
!**** Output: DUSK = time of DUSK (temporal radians) at mean local time
!****
!      Implicit Real*8 (A-H,O-Z)
      REAL*8 TWOPI,RLAT,SIND,COSD,RSmEzM,DUSK
      REAL*8 SINJ,COSJ, SJSD,CJCD, CDUSK

      Parameter (TWOPI=6.283185307179586477d0)
!****
      SINJ = Sin (TWOPI*RLAT/360)
      COSJ = Cos (TWOPI*RLAT/360)
      SJSD = SINJ*SIND
      CJCD = COSJ*COSD
      If (SJSD+RSmEzM+CJCD <= 0)  GoTo 20
      If (SJSD+RSmEzM-CJCD >= 0)  GoTo 10
!**** Compute DUSK (at local time)
      CDUSK = - (SJSD + RSmEzM) / CJCD  !  cosine of DUSK
      DUSK  = ACos (CDUSK)
      Return
!**** Constant daylight at this latitude
   10 DUSK = 999999
      Return
!**** Constant nightime at this latitude
   20 DUSK = -999999
      Return
      End

!************************************************************************
      Subroutine ORBIT (ECCEN,OBLIQ,OMEGVP, DAY, &
                        SIND,COSD,SUNDIS,SUNLON,SUNLAT,EQTIME)
!****
!**** ORBIT receives orbital parameters and time of year, and returns
!**** distance from Sun, declination angle, and Sun's overhead position.
!**** Reference for following calculations is:  V.M.Blanco and
!**** S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
!**** 135 - 151.  Existence of Moon and heavenly bodies other than
!**** Earth and Sun are ignored.  Earth is assumed to be spherical.
!****
!**** Program author: Gary L. Russell 2008/09/22
!**** Angles, longitude and latitude are measured in radians.
!****
!**** Input: ECCEN  = eccentricity of the orbital ellipse
!****        OBLIQ  = latitude of Tropic of Cancer
!****        OMEGVP = longitude of perihelion (sometimes Pi is added) =
!****               = spatial angle from vernal equinox to perihelion
!****                 with Sun as angle vertex
!****        DAY    = days measured since 2000 January 1, hour 0
!****
!**** Constants: EDAYzY = tropical year = Earth days per year = 365.2425
!****            VE2000 = days from 2000 January 1, hour 0 until vernal
!****                     equinox of year 2000 = 31 + 29 + 19 + 7.5/24
!****
!**** Intermediate quantities:
!****    BSEMI = semi minor axis in units of semi major axis
!****   PERIHE = perihelion in days since 2000 January 1, hour 0
!****            in its annual revolution about Sun
!****       TA = true anomaly = spatial angle from perihelion to
!****            current location with Sun as angle vertex
!****       EA = eccentric anomaly = spatial angle measured along
!****            eccentric circle (that circumscribes Earth's orbit)
!****            from perihelion to point above (or below) Earth's
!****            absisca (where absisca is directed from center of
!****            eccentric circle to perihelion)
!****       MA = mean anomaly = temporal angle from perihelion to
!****            current time in units of 2*Pi per tropical year
!****   TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP
!****   EAofVE = EA(VE) = eccentric anomaly of vernal equinox
!****   MAofVE = MA(VE) = mean anomaly of vernal equinox
!****   SLNORO = longitude of Sun in Earth's nonrotating reference frame
!****   VEQLON = longitude of Greenwich Meridion in Earth's nonrotating
!****            reference frame at vernal equinox
!****   ROTATE = change in longitude in Earth's nonrotating reference
!****            frame from point's location on vernal equinox to its
!****            current location where point is fixed on rotating Earth
!****   SLMEAN = longitude of fictitious mean Sun in Earth's rotating
!****            reference frame (normal longitude and latitude)
!****
!**** Output: SIND = sine of declination angle = sin(SUNLAT)
!****         COSD = cosine of the declination angle = cos(SUNLAT)
!****       SUNDIS = distance to Sun in units of semi major axis
!****       SUNLON = longitude of point on Earth directly beneath Sun
!****       SUNLAT = latitude of point on Earth directly beneath Sun
!****       EQTIME = Equation of Time =
!****              = longitude of fictitious mean Sun minus SUNLON
!****
!**** From the above reference:
!**** (4-54): [1 - ECCEN*cos(EA)]*[1 + ECCEN*cos(TA)] = (1 - ECCEN^2)
!**** (4-55): tan(TA/2) = sqrt[(1+ECCEN)/(1-ECCEN)]*tan(EA/2)
!**** Yield:  tan(EA) = sin(TA)*sqrt(1-ECCEN^2) / [cos(TA) + ECCEN]
!****    or:  tan(TA) = sin(EA)*sqrt(1-ECCEN^2) / [cos(EA) - ECCEN]
!****
!     Use C540, Only: EDAYzY,VE2000,TWOPI
!      Implicit  Real*8 (A-H,M-Z)
      REAL*8 TWOPI,EDAYzY,VE2000
      REAL*8 ECCEN,OBLIQ,OMEGVP, DAY
      REAL*8 SIND,COSD,SUNDIS,SUNLON,SUNLAT,EQTIME
      REAL*8 BSEMI, TAofVE, EAofVE, MAofVE, MA, EA, dEA, TA
      REAL*8 SUNX, SUNY, SLNORO, VEQLON, ROTATE, SLMEAN

      Parameter (TWOPI=6.283185307179586477d0, &
                 EDAYzY=365.2425d0, VE2000=79.3125)
!****
!**** Determine EAofVE from geometry: tan(EA) = b*sin(TA) / [e+cos(TA)]
!**** Determine MAofVE from Kepler's equation: MA = EA - e*sin(EA)
!**** Determine MA knowing time from vernal equinox to current day
!****
      BSEMI  = Sqrt (1 - ECCEN*ECCEN)
      TAofVE = - OMEGVP
      EAofVE = ATan2 (BSEMI*Sin(TAofVE), ECCEN+Cos(TAofVE))
      MAofVE = EAofVE - ECCEN*Sin(EAofVE)
!     PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
      MA     = Modulo (TWOPI*(DAY-VE2000)/EDAYzY + MAofVE, TWOPI)
!****
!**** Numerically invert Kepler's equation: MA = EA - e*sin(EA)
!****
      EA  = MA + ECCEN*(Sin(MA) + ECCEN*Sin(2*MA)/2)
   10 dEA = (MA - EA + ECCEN*Sin(EA)) / (1 - ECCEN*Cos(EA))
      EA  = EA + dEA
      If (Abs(dEA) > 1d-10)  GoTo 10
!****
!**** Calculate distance to Sun and true anomaly
!****
      SUNDIS = 1 - ECCEN*Cos(EA)
      TA     = ATan2 (BSEMI*Sin(EA), Cos(EA)-ECCEN)
!****
!**** Change reference frame to be nonrotating reference frame, angles
!**** fixed according to stars, with Earth at center and positive x
!**** axis be ray from Earth to Sun were Earth at vernal equinox, and
!**** x-y plane be Earth's equatorial plane.  Distance from current Sun
!**** to this x axis is SUNDIS sin(TA-TAofVE).  At vernal equinox, Sun
!**** is located at (SUNDIS,0,0).  At other times, Sun is located at:
!****
!**** SUN = (SUNDIS cos(TA-TAofVE),
!****        SUNDIS sin(TA-TAofVE) cos(OBLIQ),
!****        SUNDIS sin(TA-TAofVE) sin(OBLIQ))
!****
      SIND   = Sin(TA-TAofVE) * Sin(OBLIQ)
      COSD   = Sqrt (1 - SIND*SIND)
      SUNX   = Cos(TA-TAofVE)
      SUNY   = Sin(TA-TAofVE) * Cos(OBLIQ)
      SLNORO = ATan2 (SUNY,SUNX)
!****
!**** Determine Sun location in Earth's rotating reference frame
!**** (normal longitude and latitude)
!****
      VEQLON = TWOPI*VE2000 - TWOPI/2 + MAofVE - TAofVE  !  modulo 2*Pi
      ROTATE = TWOPI*(DAY-VE2000)*(EDAYzY+1)/EDAYzY
      SUNLON = Modulo (SLNORO-ROTATE-VEQLON, TWOPI)
               If (SUNLON > TWOPI/2)  SUNLON = SUNLON - TWOPI
      SUNLAT = ASin (Sin(TA-TAofVE)*Sin(OBLIQ))
!****
!**** Determine longitude of fictitious mean Sun
!**** Calculate Equation of Time
!****
      SLMEAN = TWOPI/2 - TWOPI*(DAY-Floor(DAY))
      EQTIME = Modulo (SLMEAN-SUNLON, TWOPI)
               If (EQTIME > TWOPI/2)  EQTIME = EQTIME - TWOPI
!****
      Return
      End

!************************************************************************
      LOGICAL Function QLEAPY (JYEAR)
!****
!**** Determine whether the given JYEAR is a Leap Year or not.
!****
      INTEGER JYEAR
      Logical JUNKLEAPY

      JUNKLEAPY = .False.
      JUNKLEAPY = JUNKLEAPY .or.  Mod(JYEAR,  4) == 0
      JUNKLEAPY = JUNKLEAPY .and. Mod(JYEAR,100) /= 0
      JUNKLEAPY = JUNKLEAPY .or.  Mod(JYEAR,400) == 0
      
      QLEAPY = JUNKLEAPY

      Return
      End

!************************************************************************
      Subroutine YMDtoD (JYEAR,JMONTH,DATE, DAY)
!****
!**** For a given JYEAR (A.D.), JMONTH and DATE (between 0 and 31),
!**** calculate number of DAYs measured from 2000 January 1, hour 0.
!****
!      Implicit Real*8 (A-H,O-Z)
      INTEGER JYEAR,JMONTH,JDAY4C,JDAY1C,JDAY4Y,JDAY1Y
      REAL*8   DATE, DAY
      Integer*4 JDSUMN(12),JDSUML(12)
      INTEGER N4CENT, IYR4C, N1CENT, IYR1C, N4YEAR, IYR4Y, N1YEAR

      Parameter (JDAY4C = 365*400 + 97, & !  number of days in 4 centuries
                 JDAY1C = 365*100 + 24, & !  number of days in 1 century
                 JDAY4Y = 365*  4 +  1, & !  number of days in 4 years
                 JDAY1Y = 365)            !  number of days in 1 year
      Data JDSUMN /0,31,59, 90,120,151, 181,212,243, 273,304,334/
      Data JDSUML /0,31,60, 91,121,152, 182,213,244, 274,305,335/
!****
      N4CENT = Floor ((JYEAR-2000)/400d0)
      IYR4C  = JYEAR-2000 - N4CENT*400
      N1CENT = IYR4C/100
      IYR1C  = IYR4C - N1CENT*100
      N4YEAR = IYR1C/4
      IYR4Y  = IYR1C - N4YEAR*4
      N1YEAR = IYR4Y
      DAY    = N4CENT*JDAY4C
      If (N1CENT > 0)  GoTo 10
!**** First of every fourth century: 16??, 20??, 24??, etc.
      DAY = DAY + N4YEAR*JDAY4Y
      If (N1YEAR > 0)  GoTo 200
      GoTo 100
!**** Second to fourth of every fourth century: 21??, 22??, 23??, etc.
   10 DAY = DAY + JDAY1C+1 + (N1CENT-1)*JDAY1C
      If (N4YEAR > 0)  GoTo 20
!**** First four years of every second to fourth century when there is
!**** no leap year during the four years: 2100-2103, 2200-2203, etc.
      DAY = DAY + N1YEAR*JDAY1Y
      GoTo 210
!**** Subsequent four years of every second to fourth century when
!**** there is a leap year: 2104-2107, 2108-2111 ... 2204-2207, etc.
   20 DAY = DAY + JDAY4Y-1 + (N4YEAR-1)*JDAY4Y
      If (N1YEAR > 0)  GoTo 200
!****
!**** Current year is a leap year
!****
  100 DAY = DAY + JDSUML(JMONTH) + DATE
      Return
!****
!**** Current year is not a leap frog year
!****
  200 DAY = DAY + JDAY1Y+1 + (N1YEAR-1)*JDAY1Y
  210 DAY = DAY + JDSUMN(JMONTH) + DATE
      Return
      End

!************************************************************************
      Subroutine DtoYMD (DAY, JYEAR,JMONTH,DATE)
!****
!**** For a given DAY measured from 2000 January 1, hour 0, determine
!**** the JYEAR (A.D.), JMONTH and DATE (between 0 and 31).
!****
      Implicit  Real*8 (A-H,O-Z)
      Integer*4 JDSUMN(12),JDSUML(12)
      Integer*8 N4CENT, N1CENT, N4YEAR, N1YEAR
      REAL*8 DAY,DATE, DAY4C, DAY1C, DAY4Y, DAY1Y
      INTEGER JYEAR,JMONTH,JDAY4C,JDAY1C,JDAY4Y,JDAY1Y, M
      Parameter (JDAY4C = 365*400 + 97, & !  number of days in 4 centuries
                 JDAY1C = 365*100 + 24, & !  number of days in 1 century
                 JDAY4Y = 365*  4 +  1, & !  number of days in 4 years
                 JDAY1Y = 365)            !  number of days in 1 year
      Data JDSUMN /0,31,59, 90,120,151, 181,212,243, 273,304,334/
      Data JDSUML /0,31,60, 91,121,152, 182,213,244, 274,305,335/
!****
      N4CENT = Floor (DAY/JDAY4C)
      DAY4C  = DAY - N4CENT*JDAY4C
      N1CENT = (DAY4C-1)/JDAY1C
      If (N1CENT > 0)  GoTo 10
!**** First of every fourth century: 16??, 20??, 24??, etc.
      DAY1C  = DAY4C
      N4YEAR = DAY1C/JDAY4Y
      DAY4Y  = DAY1C - N4YEAR*JDAY4Y
      N1YEAR = (DAY4Y-1)/JDAY1Y
      If (N1YEAR > 0)  GoTo 200
      GoTo 100
!**** Second to fourth of every fourth century: 21??, 22??, 23??, etc.
   10 DAY1C  = DAY4C - N1CENT*JDAY1C - 1
      N4YEAR = (DAY1C+1)/JDAY4Y
      If (N4YEAR > 0)  GoTo 20
!**** First four years of every second to fourth century when there is
!**** no leap year: 2100-2103, 2200-2203, 2300-2303, etc.
      DAY4Y  = DAY1C
      N1YEAR = DAY4Y/JDAY1Y
      DAY1Y  = DAY4Y - N1YEAR*JDAY1Y
      GoTo 210
!**** Subsequent four years of every second to fourth century when
!**** there is a leap year: 2104-2107, 2108-2111 ... 2204-2207, etc.
   20 DAY4Y  = DAY1C - N4YEAR*JDAY4Y + 1
      N1YEAR = (DAY4Y-1)/JDAY1Y
      If (N1YEAR > 0)  GoTo 200
!****
!**** Current year is a leap frog year
!****
  100 DAY1Y = DAY4Y
      Do 120 M=1,11
  120 If (DAY1Y < JDSUML(M+1))  GoTo 130
!     M=12
  130 JYEAR  = 2000 + N4CENT*400 + N1CENT*100 + N4YEAR*4 + N1YEAR
      JMONTH = M
      DATE   = DAY1Y - JDSUML(M)
      Return
!****
!**** Current year is not a leap frog year
!****
  200 DAY1Y  = DAY4Y - N1YEAR*JDAY1Y - 1
  210 Do 220 M=1,11
  220 If (DAY1Y < JDSUMN(M+1))  GoTo 230
!     M=12
  230 JYEAR  = 2000 + N4CENT*400 + N1CENT*100 + N4YEAR*4 + N1YEAR
      JMONTH = M
      DATE   = DAY1Y - JDSUMN(M)
      Return
      End

!************************************************************************
      Function VERNAL (JYEAR)
!****
!**** For a given year, VERNAL calculates an approximate time of vernal
!**** equinox in days measured from 2000 January 1, hour 0.
!****
!**** VERNAL assumes that vernal equinoxes from one year to the next
!**** are separated by exactly 365.2425 days, a tropical year
!**** [Explanatory Supplement to The Astronomical Ephemeris].  If the
!**** tropical year is 365.2422 days, as indicated by other references,
!**** then the time of the vernal equinox will be off by 2.88 hours in
!**** 400 years.
!****
!**** Time of vernal equinox for year 2000 A.D. is March 20, 7:36 GMT
!**** [NASA Reference Publication 1349, Oct. 1994].  VERNAL assumes
!**** that vernal equinox for year 2000 will be on March 20, 7:30, or
!**** 79.3125 days from 2000 January 1, hour 0.  Vernal equinoxes for
!**** other years returned by VERNAL are also measured in days from
!**** 2000 January 1, hour 0.  79.3125 = 31 + 29 + 19 + 7.5/24.
!****
      Implicit  Real*8 (A-H,O-Z)
      REAL*8 EDAYzY, VE2000, VERNAL
      INTEGER JYEAR
      Parameter (EDAYzY=365.2425d0, VE2000=79.3125)
      VERNAL = VE2000 + (JYEAR-2000)*EDAYzY
      Return
      End
!************************************************************************

END MODULE solar_insolation
