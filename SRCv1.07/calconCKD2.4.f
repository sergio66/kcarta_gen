c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c this is CKDv2.4, obtained from LBLRTMv5.10
c http://metosrv2.umd.edu/~bobe/LBLRTM/

c modified to make it kCARTA compatible by Sergio De Souza-Machado

c note that since we want to do jacobians, the code has been rewritten 
c so that it finds h20s96,h20s260,f296 at the coarse spacing, and then
c immediately linearly interpolates this onto the finer kCARTA grid
c this means subroutines SL296,SL260,F296 are slightly modified from the
c LBLRTM module, as they have to output the coarse grid

c************************************************************************
       SUBROUTINE CALCON24( IDGAS, NFREQ, raFreq, FSTEP, NLAY, raT, raP, 
     $    raPartP, raAMNT, raaCON, raadDQ, raadDT, iDoDQ,
     $    rSelfMult,rForMult)

c      this is CKD v 2.4

C     SUBROUTINE CONTNM CONTAINS THE CONTINUUM DATA 
C     WHICH IS INTERPOLATED INTO THE ARRAY ABSRB    

      include 'kcarta.param'

c      IMPLICIT NONE

c ----------  these arguments are from calcon.f -------------------
C COMMON BLOCKS CONTAINING CONTINUUM VALUES 

      INTEGER NPTS1,NPTS2,NPTF,NPTFCO2,NPN2RT0,NPN2RT1
      INTEGER NPN2_F,NPO2_F
      INTEGER NPTO3CHAP,NPTO3HH0,NPTO3HH1,NPTO3HH2,NPTO3HUV
      REAL V1S1,V2S1,DVS1
      REAL V1S2,V2S2,DVS2
      REAL V1F,V2F,DVF
      REAL V1FCO2,V2FCO2,DVFCO2
      REAL V1N2RT0,V2N2RT0,DVFN2RT0
      REAL V1N2RT1,V2N2RT1,DVFN2RT1
      REAL V1N2_F,V2N2_F,DVFN2_F
      REAL V1O2_F,V2O2_F,DVFO2_F
      REAL V1O3CHAP,V2O3CHAP,DVO3CHAP
      REAL V1O3HH0,V2O3HH0,DVO3HH0
      REAL V1O3HH1,V2O3HH1,DVO3HH1
      REAL V1O3HH2,V2O3HH2,DVO3HH2
      REAL V1O3HUV,V2O3HUV,DVO3HUV

      REAL raSSFREQ(kLBL) 
      REAL H2OS96(2003),H2OS60(2003),H2OF(2003),CO2F(1003)
      REAL N296(73),N220(73),xn(118),xnt(118)
      REAL X(3150),Y(3150),Z(3150)
      REAL o3HH0(2687),o3HH1(2687),o3HH2(2687),o3HUV(133),o2F(206)

      COMMON /SH2O/ V1S1,V2S1,DVS1,NPTS1,H2OS96
      COMMON /S260/ V1S2,V2S2,DVS2,NPTS2,H2OS60
      COMMON /FH2O/ V1F,V2F,DVF,NPTF,H2OF
      COMMON /FCO2/ V1FCO2,V2FCO2,DVFCO2,NPTFCO2,CO2F
      COMMON /N2RT0/V1N2RT0,V2N2RT0,DVFN2RT0,NPN2RT0,N296
      COMMON /N2RT1/V1N2RT1,V2N2RT1,DVFN2RT1,NPN2RT1,N220
      COMMON /n2_f/ V1N2_F,V2N2_F,DVFN2_F,NPN2_F,xn,xnt
      COMMON /O3CHAP/ V1O3CHAP,V2O3CHAP,DVO3CHAP,NPTO3CHAP,X,Y,Z
      COMMON /O3HH0/ V1O3HH0,V2O3HH0,DVO3HH0,NPTO3HH0,o3HH0
      COMMON /O3HH1/ V1O3HH1,V2O3HH1,DVO3HH1,NPTO3HH1,o3HH1
      COMMON /O3HH2/ V1O3HH2,V2O3HH2,DVO3HH2,NPTO3HH2,o3HH2
      COMMON /O3HUV/ V1O3HUV,V2O3HUV,DVO3HUV,NPTO3HUV,o3HUV
      COMMON /o2_f  /V1O2_F,V2O2_F,DVFO2_F,NPO2_F,o2F 

C Arguements
       INTEGER IDGAS, NFREQ, NLAY, iDoDQ
       REAL raFreq(kMaxPts), FSTEP, raT(kProfLayer), raP(kProfLayer),
     $    raPartP(kProfLayer),raAmnt(kProfLayer)
       REAL raaCON(kMaxPts,kProfLayer)
       REAL raadDQ(kMaxPtsJac,kProfLayerJac)
       REAL raadDT(kMaxPtsJac,kProfLayerJac)
       REAL rSelfMult,rForMult

C Local variables
       REAL raXAMT(kProfLayer),raC11(kProfLayer),
     $    raTFAC(kProfLayer), A1, A2, A3

C Variables for new water continuum adjustment terms
       INTEGER JFAC
       REAL SFAC, FSCAL, XFAC(0:50)

c used in calculating d/dq, d/dT
      REAL ra1,ra2,ra3,rQ,rF

c general terms used in self, foreign
      INTEGER iL

c ---------------------------------------------------------------

       REAL vs2,vf2,vf6,vf4

       REAL V0S1, V0S2, V0S3, HWSQ1, HWSQ2, HWSQ3, 
     $    BETAS1, BETAS3, 
     $    FACTRS1, FACTRS2, FACTRS3

       REAL V0F2, HWSQF2, V0F3

       REAL rSH2O, VS4, V0F1, HWSQF1, BETAF1, FACTRF1,
     $    V0F1a, HWSQF1a, HWSQF3, BETAF1a,  FACTRF1a, BETAF2, 
     $    FACTRF2, BETAF3, FACTRF3
      
      REAL V1ABS,V2ABS,DVABS,NPTABS,ABSRB(kLBL)
      REAL V1C,V2C,DVC
      REAL raContSelf(kMaxPts),raContFor(kMaxPts)
      REAL raSH2OT0(kMaxPts),raSH2OT1(kMaxPts),raFH2O(kMaxPts) 
      INTEGER NPTC

      COMMON /ABSORB/V1ABS,V2ABS,DVABS,NPTABS,ABSRB
                                                   
      REAL XSELF,XFRGN
 
      REAL SH2OT0(kLBL),SH2OT1(kLBL),FH2O(kLBL)
cc     *          CN2T0(kLBL),FCO2(kLBL),CT1(kLBL),CT2(kLBL) 
cc      REAL CCH0(5150),CCH1(5150),CCH2(5150)

cc      REAL ABSBSV(kLBL), XLOSMT
                                                             
cc      EQUIVALENCE (C0,SH2OT0,CN2T0,FCO2) , (C1,SH2OT1,CT1),             
cc     *            (C2,FH2O,CT2)                                         
             
      REAL P0, T0                                                
      DATA P0 / 1013. /,T0 / 296. /                                     
cc      DATA XLOSMT / 2.68675E+19 /                                       

      REAL TS,AVOG
      DATA AVOG/6.022045E+26/
      DATA TS/296.0/

      REAL V1
      INTEGER iI,iJ

      REAL SFACL,SFACH,FL

c     These are self-continuum modification factors from 700-1200 cm-1
      DATA (XFAC(iI),iI=0,50)/
     1    1.00000,1.01792,1.03767,1.05749,1.07730,1.09708,
     2    1.10489,1.11268,1.12047,1.12822,1.13597,1.14367,
     3    1.15135,1.15904,1.16669,1.17431,1.18786,1.20134,
     4    1.21479,1.22821,1.24158,1.26580,1.28991,1.28295,
     5    1.27600,1.26896,1.25550,1.24213,1.22879,1.21560,
     6    1.20230,1.18162,1.16112,1.14063,1.12016,1.10195,
     7    1.09207,1.08622,1.08105,1.07765,1.07398,1.06620,
     8    1.05791,1.04905,1.03976,1.02981,1.00985,1.00000,
     9    1.00000,1.00000,1.00000/
                                                             
C     ASSIGN SCCS VERSION NUMBER TO MODULE 
C
C     Continuum calculation flags:
C     ---------------------------
C     ICNTNM Value      Self     Foreign    Rayleigh     Others
C           0            no        no          no          no
C           1            yes       yes         yes         yes
C           2            no        yes         yes         yes
C           3            yes       no          yes         yes
C           4            no        no          yes         yes
C           5            yes       yes         no          yes
C           6   READ IN XSELF, XFRGN, XCO2C, XO3CN, XO2CN, XN2CN, 
C               and XRAYL in Record 1.2a
C


C=======================================================================

C               ********    WATER VAPOR   ********                      

C H2O (water) continuum 
      IF (IDGAS .NE. 1)THEN 
        write(kStdErr,*)'in CKD2.4, need gasID = 1' 
        CALL DoSTOP 
        END IF 

c initialize a few common variables here
       DVABS = 1.         
       DVABS = FStep
       V1ABS = INT(raFreq(1))    
       IF (V1.LT.0.) V1ABS = V1ABS-1.
       V1ABS = V1ABS-3.*DVABS        
       V2ABS = INT(raFreq(kMaxPts)+3.*DVABS+0.5)  
       NPTABS = (V2ABS-V1ABS)/DVABS+1.5
       XSELF=1.0                                  
       XFRGN=1.0                                  
C=======================================================================
C                             SELF


      !this is to put it into Genln2/kCARTA units
      DO iL=1,kProfLayer                     !loop over layers
        raTFAC(iL) = (raT(iL)-T0)/(260.-T0)              
        raXamt(iL)=raAmnt(iL)*avog
        raC11(iL)=0.5*kPlanck2/raT(iL)
        END DO

      IF ((raFreq(kMaxPts).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN
         ! for self, foreign broadening the v1c,v2c,dvc,nptc parameters
         ! are all the same (sl296, sl260,frn296)
         CALL SL296(raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT0)              
         CALL SL260(raSSFREQ,V1C,V2C,DVC,NPTC,SH2OT1)            
         !now linearly interp or spline them onto the finer grid 

         CALL XLINEAR(raSSFreq,SH2OT0,nptc,raFreq,raSH2OT0,nfreq) 
         CALL XLINEAR(raSSFreq,SH2OT1,nptc,raFreq,raSH2OT1,nfreq) 
         END IF

C     Only calculate if V2 > -20. cm-1 and V1 <  20000. cm-1
      IF ((raFreq(kMaxPts).gt.-20.0).AND.(raFreq(1).lt.20000.)) THEN

        DO iL=1,kProfLayer                     !loop over layers

C--------------------------------------------------------------------
C     *****    Continuum Correction Patches    ********
C--------------------------------------------------------------------
c                   SELF
          V0S1 = 0.
          HWSQ1 = 100.**2
          BETAS1 = 1.E-04
c         FACTRS1 = 0.3                      ! CKD2.2 value
          FACTRS1 = 0.688
 
          V0S2 = 1050.
          HWSQ2 = 200.**2
          FACTRS2 = -0.2333
 
          V0S3 = 1310.
          HWSQ3 = 120.**2
          BETAS3 = 5.E-06
          FACTRS3 = -0.15
C--------------------------------------------------------------------
c         Loop calculating self continuum optical depth
          DO 20 iJ = 1, kMaxPts
            RF = raFreq(iJ)
            rSH2O = 0.                            
            IF (raSH2OT0(iJ).GT.0.) THEN            
              rSH2O=raSH2OT0(iJ)*(raSH2OT1(iJ)/raSH2OT0(iJ))**raTFAC(iL) 
      
              SFAC = 1.
              IF (RF.GE.700. .AND.  RF.LE.1200.) THEN 
                JFAC = (RF-700.)/10. + 0.00001
                SFAC = XFAC(JFAC)

ccc this is new, copied from CKD2.1 make interpolation smoother across 10cm-1
                JFAC=INT( (RF - 700.0)/10.0 + 0.0001)
                SFACL=XFAC(JFAC)
                SFACH=XFAC(JFAC+1)
                FL=700.0 + FLOAT(JFAC*10)
                SFAC=SFACL + (RF - FL)*(SFACH - SFACL)/10.0

                ENDIF
      
C       ---------------------------------------------------------
C         Correction to self continuum (1 SEPT 85); factor of    
C         0.78 at 1000 and  .......
                                  
              VS2 = (RF-V0S1)**2
              VS4 = VS2*VS2
              SFAC = SFAC*(1.+FACTRS1*(HWSQ1/(RF**2+(BETAS1*VS4)+HWSQ1)))  
 
              VS2 = (RF-V0S2)**2
              SFAC = SFAC*(1.+FACTRS2*(HWSQ2/(VS2+HWSQ2)))
    
              VS2 = (RF-V0S3)**2
              VS4 = VS2*VS2
              SFAC = SFAC*(1.+FACTRS3*(HWSQ3/(VS2+(BETAS3*VS4)+HWSQ3))) 
                                            
              rSH2O = SFAC * rSH2O
C    ---------------------------------------------------------
              ENDIF                 !  IF (SH2OT0(iJ).GT.0.) THEN            
 
c            raContSelf(iJ) = (rSH2O*RH2O)*XSELF
            raContSelf(iJ) = (raPartP(iL)*rSH2O)*XSELF
                                  
 20         CONTINUE                                           

C=======================================================================
c                     FOREIGN

          V0F1 = 350.
          HWSQF1 = 200.**2
          BETAF1 = 5.e-09 
          FACTRF1 = -0.7
 
          V0F1a = 630.
          HWSQF1a = 65.**2
          BETAF1a = 2.e-08 
          FACTRF1a = +0.75
 
          V0F2 =1130.
          HWSQF2 = 330.**2
          BETAF2 = 8.E-11
          FACTRF2 = -0.97
 
          V0F3 = 1975.
          HWSQF3 = 250.**2
          BETAF3 = 5.E-06
          FACTRF3 = -0.65
 
c        ------------------------------------------------------------

          CALL FRN296(raSSFREQ,V1C,V2C,DVC,NPTC,FH2O)         
          CALL XLINEAR(raSSFreq,FH2O,nptc,raFreq,raFH2O,nfreq) 

          DO 24 iJ = 1, kMaxPts
            RF = raFreq(iJ)
   
C           CORRECTION TO FOREIGN CONTINUUM
   
            VF2 = (RF-V0F1)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = (1.+FACTRF1*(HWSQF1/(VF2+(BETAF1*VF6)+HWSQF1)))
 
            VF2 = (RF-V0F1a)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF1a*(HWSQF1a/(VF2+(BETAF1a*VF6)+HWSQF1a)))
 
            VF2 = (RF-V0F2)**2
            VF6 = VF2 * VF2 * VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF2*(HWSQF2/(VF2+(BETAF2*VF6)+HWSQF2)))
 
            VF2 = (RF-V0F3)**2
            VF4 = VF2*VF2
            FSCAL = FSCAL* 
     *              (1.+FACTRF3*(HWSQF3/(VF2+BETAF3*VF4+HWSQF3)))
      
            raFH2O(iJ)=raFH2O(iJ)*FSCAL
      
c            raContFor(iJ) = (raFH2O(iJ)*RFRGN)*XFRGN
            raContFor(iJ) = (raFH2O(iJ)*(raP(iL)-raPartP(iL)))*XFRGN
 
 24         CONTINUE                                                  
 
C           ------------------------------------------------------------
c       Interpolate to total optical depth grid
c       If spectral range covers the CKD_2.3 1 cm-1 H2O continuum
c       as well, then stop interpolation at 2200 cm-1.

c          npts_lo = 1
c          if (v1abs .lt. 2200.) then
c             npts_lo =  npts_hi + 1
c             v1c = v1c+npts_hi
c             endif

C         The factor of 1.e-20 is handled this way to avoid underflows
          DO iJ=1,kMaxPts     !was NPTC
            RF = raFreq(iJ)   !was V1C+DVC*FLOAT(iJ-1)
            a1=RF*raXamt(iL)*tanh(raC11(iL)*RF)
            a2=TS/raT(iL)
            a3=1.0e-20*(rSelfMult*raContFor(iJ) + rForMult*raContSelf(iJ))
            raaCon(iJ,iL)=a1*a2*a3

c see if we need to do d/dT and d/dq 
            IF ((kJacobian .GT. 0).AND.(iDoDQ .GT. 0)) THEN  
   
              ra1=EXP( raC11(iL)*RF)+EXP(-raC11(iL)*RF) 
              ra1=1.0/(ra1*ra1) 
              ra1=-ra1*2.0*RF*RF*kPlanck2*raXamt(iL) 
              ra1=ra1/(raT(iL)*raT(iL)) 
 
              ra2=-TS/(raT(iL)*raT(iL)) 
   
              ra3=((raSH2OT1(iJ)/raSH2OT0(iJ))**(raTfac(iL))) 
              ra3=ra3*alog(raSH2OT1(iJ)/raSH2OT0(iJ)) 
              ra3=ra3*raSH2OT0(iJ) 
              ra3=ra3*(-2.7777777e-2)*1.0e-20*raPartP(iL)*rSelfMult
 
              rQ=A1/raAMNT(iL) 
              raadDQ(iJ,iL)=rQ*A2*A3 
   
              raadDT(iJ,iL)=ra1*A2*A3+A1*ra2*A3+A1*A2*ra3 
              END IF 
 
c see if we only need to do d/dT 
            IF ((kJacobian .GE. 0).AND.(iDoDQ .LT. 0)) THEN  
   
              ra1=EXP( raC11(iL)*RF)+ EXP(-raC11(iL)*RF) 
              ra1=1.0/(ra1*ra1) 
              ra1=-ra1*2.0*RF*RF*kPlanck2*raXamt(iL) 
              ra1=ra1/(raT(iL)*raT(iL)) 
 
              ra2=-TS/(raT(iL)*raT(iL)) 
   
              ra3=((raSH2OT1(iJ)/raSH2OT0(iJ))**(raTfac(iL))) 
              ra3=ra3*alog(raSH2OT1(iJ)/raSH2OT0(iJ)) 
              ra3=ra3*raSH2OT0(iJ) 
              ra3=ra3*(-2.7777777e-2)*1.0e-20*raPartP(iL)*rSelfMult
   
              raadDT(iJ,iL)=ra1*A2*A3+A1*ra2*A3+A1*A2*ra3 
   
              END IF  !IF ((kJacobian .GE. 0).AND.(iDoDQ .LT. 0)) THEN  
 
            END DO             !end loop over freqs

          END DO               !end loop over layers
        ENDIF

      RETURN
      END
c************************************************************************

      SUBROUTINE SL296(raF,V1C,V2C,DVC,NPTC,C)                            
            
      include 'kcarta.param'
                          
      REAL V1C,V2C,DVC
      REAL C(*),raF(*)                                        
      INTEGER NPTC
    
      REAL V1ABS,V2ABS,DVABS,ABSRB(kLBL)             
      REAL V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /SH2O/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                        
      V1C = V1ABS-DVC                                                  
      V2C = V2ABS+DVC         
                                          
      I1 = (V1C-V1S)/DVS                                               
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       

      DO 10 J = 1, NPTC                                                
         raF(J)=v1c + (j-1)*dvc
         I = I1+J                                                      
         C(J) = 0.                                                     
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                         
         C(J) = S(I)                                                   
   10 CONTINUE                                                         

      RETURN                                                           
                                          
      END                                                              
 
c************************************************************************
      SUBROUTINE SL260 (raF,V1C,V2C,DVC,NPTC,C)                              

      include 'kcarta.param'

      REAL V1C,V2C,DVC
      REAL C(*),raF(*)                                      
      INTEGER NPTC
    
      REAL V1ABS,V2ABS,DVABS,ABSRB(kLBL)             
      REAL V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /S260/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC
         raF(J)=v1c + (j-1)*dvc                       
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.LT.1).OR.(I.GT.NPTS)) GO TO 10                           
         C(J) = S(I)                                                     
   10 CONTINUE                                                           

      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
      SUBROUTINE FRN296 (raF,V1C,V2C,DVC,NPTC,C)                             
                                       
      include 'kcarta.param'
     
      REAL V1C,V2C,DVC
      REAL C(*),raF(*)                                     
      INTEGER NPTC
    
      REAL V1ABS,V2ABS,DVABS,ABSRB(kLBL)             
      REAL V1S,V2S,DVS,S(2003) 
      INTEGER NPTABS,NPTS

      COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB
      COMMON /FH2O/ V1S,V2S,DVS,NPTS,S                           
      
      INTEGER I1,I2,I,J

      DVC = DVS                                                          
      V1C = V1ABS-DVC                                                    
      V2C = V2ABS+DVC                                                    
                                            
      I1 = (V1C-V1S)/DVS                                                 
      IF (V1C.LT.V1S) I1 = -1 
        
      V1C = V1S+DVS*FLOAT(I1)        
      I2 = (V2C-V1S)/DVS             
      NPTC = I2-I1+3                 
      IF (NPTC.GT.NPTS) NPTC=NPTS+1
      V2C = V1C+DVS*FLOAT(NPTC-1)       
 
      DO 10 J = 1, NPTC
         raF(J)=v1c + (j-1)*dvc
         I = I1+J                                                        
         C(J) = 0.                                                       
         IF ((I.GE.1).AND.(I.LE.NPTS)) THEN                              
            C(J) = S(I)                                                  
         ENDIF                                                           
   10 CONTINUE                                                           

      RETURN                                                             
                                            
      END                                                                
 
c************************************************************************
