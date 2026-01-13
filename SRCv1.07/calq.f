       SUBROUTINE CALQ (IDXSEC, TPROF, Q, idQT, rDT)

C      Calculate total internal partition function for XSEC molecule
C      Also allows for simple computation of d/dT if kJacobian=kSpline=1 ...
c           this is enabled by setting idQT=+1 at calling time

       INTEGER NXSEC,NTMP,ID1
       PARAMETER (NXSEC=13, NTMP=18, ID1=51)
C
C      Arguements
       INTEGER IDXSEC
       REAL TPROF, Q
c these next two allow for computation of d/dt if idQT=1
       REAL rDT
       INTEGER idQT
C
C      Local Variables
       REAL TREF, TMP(NTMP), VIB(NXSEC,NTMP), QROT, VIBT, QVIB
       INTEGER IX, JREF, J, JLO, JHI
       
       REAL dqrot,dqvib

C      -----------------------------------------------------------------
C      Data block
C      ----------
C
       DATA TREF /296.0/
       DATA JREF /14/
C
C      Temperatures for vibration partition function interpolation
       DATA (TMP(J),J=1,18)/
     $ 170.0,    180.0,    190.0,    200.0,    210.0,    220.0,
     $ 230.0,    240.0,    250.0,    260.0,    270.0,    280.0,
     $ 290.0,    296.0,    300.0,    310.0,    320.0,    330.0/
C
C      Data for XSEC vibrational partition functions
C
C      IDXSEC=51: CFCL3 (F11)
       DATA (VIB(51 + 1 - ID1, J),J=1,18)/
     $ 1.51031,  1.61072,  1.72276,  1.84728,  1.98524,  2.13770,
     $ 2.30584,  2.49093,  2.69438,  2.91772,  3.16262,  3.43089,
     $ 3.72449,  3.91368,  4.04554,  4.39633,  4.77933,  5.19722/
C
C      IDXSEC=52: CF2Cl2 (F12)
       DATA (VIB(52 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=53: CClF3 (F13)
       DATA (VIB(53 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=54: CF4 (F14)
       DATA (VIB(54 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=55: CHCl2F (F21)
       DATA (VIB(55 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=56: CHClF2 (F22)
       DATA (VIB(56 + 1 - ID1, J),J=1,18)/
     $ 1.08491,  1.10448,  1.12629,  1.15033,  1.17665,  1.20528,
     $ 1.23625,  1.26961,  1.30541,  1.34372,  1.38461,  1.42816,
     $ 1.47446,  1.50359,  1.52359,  1.57567,  1.63079,  1.68909/
C
C      IDXSEC=57: C2Cl3F3 (F113)
       DATA (VIB(57 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=58: C2Cl2F2 (F114)
       DATA (VIB(58 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=59: C2ClF5 (F115)
       DATA (VIB(59 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      IDXSEC=60: CCl4
       DATA (VIB(60 + 1 - ID1, J),J=1,18)/
     $ 1.82149,  1.98259,  2.16373,  2.36682,  2.59396,  2.84752,
     $ 3.13009,  3.44455,  3.79404,  4.18202,  4.61229,  5.08897,
     $ 5.61659,  5.95964,  6.20006,  6.84474,  7.55644,  8.34148/
C
C      IDXSEC=61: ClONO2
       DATA (VIB(61 + 1 - ID1, J),J=1,18)/
     $ 1.80702,  1.91774,  2.03739,  2.16660,  2.30604,  2.45642,
     $ 2.61853,  2.79320,  2.98129,  3.18376,  3.40161,  3.63590,
     $ 3.88777,  4.04783,  4.15843,  4.44916,  4.76131,  5.09634/
C
C      IDXSEC=62: N2O5
       DATA (VIB(62 + 1 - ID1, J),J=1,18)/
     $ 1.13198,  1.16234,  1.19659,  1.23494,  1.27764,  1.32494,
     $ 1.37716,  1.43460,  1.49765,  1.56668,  1.64214,  1.72450,
     $ 1.81429,  1.87196,  1.91207,  2.01846,  2.13414,  2.25985/
C
C      IDXSEC=63: HNO4
       DATA (VIB(63 + 1 - ID1, J),J=1,18)/
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,
     $ 1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000/
C
C      -----------------------------------------------------------------
C
C      Set array index for XSEC molecule ID
       IX=1 + IDXSEC - ID1
C
C      Calculate the rotational partition function
       QROT = (TREF/TPROF)**1.5
C
C      Calculate the vibrational partition function (interpolate
C      for temperature)
C
       JLO=1
       JHI=2
C      Find which temperatures to use
 20    IF (TPROF .GT. TMP(JHI) .AND. JHI .LT. NTMP) THEN
          JLO=JHI
          JHI=JHI + 1
          GOTO 20
       ENDIF
C
       VIBT = VIB(IX, JLO) + ( TPROF - TMP(JLO) )*
     $    ( VIB(IX, JHI) - VIB(IX, JLO) )/( TMP(JHI) - TMP(JLO) )
C
       QVIB = VIB(IX, JREF)/VIBT
C
       Q = QROT*QVIB

       IF (idQT .EQ. 1) THEN
c if kJacobian == 0,1, and kSPline == 1, then exactly calculate d/dT from data
         dqrot=(-1.5)*QROT/(TPROF)
         dqvib= ( VIB(IX, JHI) - VIB(IX, JLO) )/( TMP(JHI) - TMP(JLO) )
         dqvib=dqvib*(-1.0*QVIB/VIBT)
         rDT=dqrot*QVIB+dqvib*QROT
       ELSE
         rDT=0.0
         END IF

       RETURN
       END            
