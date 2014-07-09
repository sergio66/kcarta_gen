c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c*******************   this is for the JACOBIANS ************************
c************************************************************************
c*****************  NEW  (SPLINE) UNCOMPRESSIONS ROUTINES ***************
c************************************************************************
c this subroutine does the spline interpolation across the temperatures,
c followed by the "un"compression to yield the absorption coeff matrix
c this is (MAINLY!!) for gases other than water even though d/dT(water) 
c uses this subroutine
      SUBROUTINE GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $  raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  daaDQ,daaDT,iDoDQ,raQAirs,iGasID)

      include 'kcarta.param'

c raQAirs   = ratio of kLAYERS/KAIRS reference gas amounts  
c daaAbsCoeff = abs coeff matrix, after the temperature interpolations
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx      = the 11 compressed matrices that will be interpolated
c daaUx       = the uncompression matrix
c raP/Rtemp   = actual/reference temperature profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = number of layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
c daaDT       = analytic Jacobian wrt Temperature
c daaDQ       = analytic Jacobian wrt amount
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),raQAirs(kProfLayer)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer),
     $    daaUx(kMaxPts,kMaxK)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iDoDQ,iGasID

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB
c for the jacobian
      DOUBLE PRECISION daaT(kMaxK,kProfLayer)

      INTEGER iI,iL,iM,iN,iK

      CALL SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx, 
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $   raQAirs,iGasID,daaT,1) 

c multiply daaUx with daaKpro to get daaAbsCoeff
ccccc user supplied info
c this is the assembly language matrix multiplication
c  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,
     $       iLDB,dBeta,daaAbsCoeff,iLDC)
       
      IF (kJacobian .GT. 0) THEN
        CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                       iLDA,iLDB,iLDC,dbeta,iUm,iUn)
        CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,
     $         iLDB,dBeta,daaDT,iLDC)
        END IF

      IF  ((kJacobian .GT. 0)  .AND. (iDoDQ .GT. 0)) THEN
        DO  iI=1,kMaxPtsJac
          DO iL=1,kProfLayerJac
            daaDQ(iI,iL)=daaAbsCoeff(iI,iL)
            END DO
          END DO 
        END IF

      RETURN
      END
c************************************************************************
c this subroutine does the interpolation of a compressed matrix, in temperature
      SUBROUTINE SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx,
     $        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,iGasID,daaT,iActuallyDoDT)

      include 'kcarta.param'
      include 'NewRefProfiles/outincLAY.param' 

c iDoDT       = do we actually do d/dT (as this is called 5 useless times
c                                        by GetAbsCoeffWater)
c daakPro     = final temperature interpolated matrix
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx      = the 11 compressed matrices that will be interpolated in temp
c raP/Rtemp   = actual/reference temperature profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = number of layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),raQAirs(kProfLayer)
      DOUBLE PRECISION daaKpro(kMaxk,kProfLayer),
     $    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer)
c for the jacobian
      DOUBLE PRECISION daaT(kMaxK,kProfLayer)
      DOUBLE PRECISION FirstDeriv
      INTEGER iActuallyDoDT

      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID

C     for interpolating daaaKX in temperature 
      DOUBLE PRECISION daWork(kMaxTemp),dYP1,dYPN,dXPT, 
     $    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp) 

C     for interpolating daaaKX in pressure 
      DOUBLE PRECISION daWorkP(kMaxLayer), 
     $    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer) 
 
C this is the actual matrix that will be created after interpolation in  
c pressure, and then interpolated in temperature 
      DOUBLE PRECISION daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d 
      INTEGER iI,iJ,iK,iM,loffset,loffset2 
 
C     Assign values for interpolation 
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points 
      dYP1=1.0E+16 
      dYPN=1.0E+16 
 
      DO iI=1,iKm 
        daXgiven(iI)=daToffset(iaTsort(iI)) 
        ENDDO 

      IF ((MSubLayer .EQ. 1) .AND. (MThickLayer .EQ. 1)) THEN 
c       user has not changed layering from AIRS, so just go ahead 
c       and set daaaKxNew=daaaKx 
        DO iI=1,iNk 
          DO iJ=1,iKm 
            DO iK=1,kProfLayer 
              daaaKxNew(iI,iJ,iK)=daaaKx(iI,iJ,iK) 
              END DO 
            END DO 
          END DO 
      ELSE 
c  even if kProfLayer =========== kMaxLayer, allow for possibility that 
c  user has changed layering, so we have to do PRESSURE interpolation 
c  of daaaKx onto daaaKnNew 
c  first do the BOTTOM set of layers, which match the AIRS layers 
        DO iI=1,iNk 
          DO iJ=1,iKm 
            DO iK=1,(M1000mb-1) 
              daaaKxNew(iI,iJ,iK)=daaaKx(iI,iJ,iK) 
              END DO 
            END DO 
          END DO 
c  then do the MIDDLE set of layers, which match the AIRS layers 
        loffset=0 
        loffset=loffset + (M1000mb-1) 
        loffset=loffset + (M100mb-M1000mb)*MSubLayer 
        loffset=loffset + 1 
        DO iI=1,iNk 
          DO iJ=1,iKm 
            DO iK=M100mb,M50mb-1 
              daaaKxNew(iI,iJ,loffset+(iK-M100mb))=daaaKx(iI,iJ,iK) 
              END DO 
            END DO 
          END DO 
c  then do the TOP set of layers, which match the AIRS layers 
        loffset=0 
        loffset=loffset + (M1000mb-1) 
        loffset=loffset + (M100mb-M1000mb)*MSubLayer 
        loffset=loffset + (M50mb-M100mb) 
        loffset=loffset + INT((M10mb-M50mb)/(MThickLayer*1.0)) 
        loffset=loffset + 1 
        DO iI=1,iNk 
          DO iJ=1,iKm 
            DO iK=M10mb,kMaxLayer 
              daaaKxNew(iI,iJ,loffset+(iK-M10mb))=daaaKx(iI,iJ,iK) 
              END DO 
            END DO 
          END DO 
c do the SUBDIVIDED/THICKENED set of layers ..need AIRS layers interpolations 
        DO iI=1,iNk 
          DO iJ=1,iKm 
            !first set up the daXgivenP,daYgivenP arrays 
            DO iK=1,kMaxLayer    
              !notice how daXgiven is initialised with increasing pressure 
              !do interpolate in log(pressure) 
              daXgivenP(iK)=log(pavgAIRS(kMaxLayer-iK+1)) 
              daYgivenP(iK)=daaaKx(iI,iJ,kMaxLayer-iK+1) 
              END DO 
            CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN, 
     $                    daY2P,daWorkP) 
c  do the SUBDIVIDED set of layers ..need AIRS layers interpolations 
            IF (MSubLayer .GT. 1) THEN 
              loffset=M1000mb 
              DO iK=M1000mb,M100mb-1 
                DO iM=1,MSubLayer 
                  dxpt=log(pProf(loffset)) 
                  CALL dsplin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d) 
                  IF (iGASID .gt. 2) THEN 
                    daaaKxNew(iI,iJ,loffset)=d*raQAirs(loffset) 
                  ELSE 
                    daaaKxNew(iI,iJ,loffset)=d*(raQAirs(loffset)**0.25) 
                    END IF 
                  loffset=loffset+1 
                  END DO 
                END DO 
            ELSE   
              DO iK=M1000mb,M100mb-1 
                daaaKxNew(iI,iJ,iK)=daaaKx(iI,iJ,iK) 
                END DO 
              END IF 
c  do the THICKENED set of layers ..need AIRS layers interpolations 
            IF (MThickLayer .GT. 1) THEN 
              loffset2=0 
              loffset2=loffset2 + (M1000mb-1) 
              loffset2=loffset2 + (M100mb-M1000mb)*MSubLayer 
              loffset2=loffset2 + (M50mb-M100mb) 
              loffset2=loffset2 + 1 
              DO iK=M50mb,M10mb-1 
                dxpt=log(pProf(loffset2)) 
                CALL dsplin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d) 
                IF (iGASID .gt. 2) THEN 
                  daaaKxNew(iI,iJ,loffset2)=d*raQAirs(loffset2) 
                ELSE 
                  daaaKxNew(iI,iJ,loffset2)=d*(raQAirs(loffset2)**0.25) 
                  END IF 
                loffset2=loffset2+1 
                END DO 
            ELSE   
              DO iK=M50mb,M10mb-1 
                daaaKxNew(iI,iJ,iK)=daaaKx(iI,iJ,iK) 
                END DO 
              END IF 
            END DO 
          END DO 
        END IF 

C     now do the spline Interpolation of the K vectors in TEMPERATURE 
      DO iI=1,iNk            !Loop over the K vectors 
        DO iK=1,kProfLayer   !Loop over the layers 
          dXPT=raPtemp(iK) - raRTemp(iK) 
          DO iJ=1,iKm        !Interpolate KX across iKm for the profile temp 
            daYgiven(iJ)=daaaKxNew(iI, iaTsort(iJ), iK) 
            ENDDO 
          CALL DSPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork) 
c subtract the ref temp from the profile temperature 
          dXPT=raPtemp(iK) - raRTemp(iK) 
          CALL DSPLIN(daXgiven,daYgiven,daY2,iKm,dXPT,daaKpro(iI,iK)) 
          IF ((kJacobian .GT. 0) .AND. (iActuallydoDT .GT. 0)) THEN
            daaT(iI,iK)=FirstDeriv(dXPT,daXgiven,daYgiven,daY2,iKm)
            END IF
          ENDDO 
        ENDDO 

      RETURN
      END 
c************************************************************************
c interpolate the compressed data in temperature AND water amount
      SUBROUTINE GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset,
     $  daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,
     $  raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  daaDQ,daaDT,iDoDQ,raQAirs,iGasID)

      include 'kcarta.param'

c raQAirs   = ratio of kLAYERS/KAIRS reference gas amounts  
c daaAbsCoeff = abs coeff matrix, after the temperature interpolations
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx1..5  = the 11 compressed matrices that will be interpolated in temp
c  after which the water partial pressure * 0.01,1.0,3.3,6.7,10.0 interpolation
c daaUx       = the uncompression matrix
c raP/Rtemp   = actual/reference temperature profiles
c raP/Rpart   = actual/reference partial pressure profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = number of layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
c daaDQ      = analytic Jacobian wrt water amount
c daaDT      = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),raQAirs(kProfLayer)
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK),
     $  daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer),
     $  daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer),
     $  daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer),
     $  daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer),
     $  daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iDODQ
      INTEGER iGasID
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB

c this is to calculate d/dq,d/dT
      DOUBLE PRECISION daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION daaQ(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer),
     $                 daaK_from_Temp(kMaxPts,kProfLayer)

      INTEGER iM,iN,iK

c first interpolate in each of the five water offset matrices, in temperature
c and in pressure to obtain the five current temp profile water matrices
c hence daaAn = the n th matrix, interpolated in layer temperature T
c note that daaT is just a dummy role here!!!!!!!!!!!!!!!!!!
      CALL SplineTempInterpolateJAC(daaA1,daToffset,daaaKx1,
     $        raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,1,daaT,-1)
      CALL SplineTempInterpolateJAC(daaA2,daToffset,daaaKx2,
     $        raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,1,daaT,-1)
      CALL SplineTempInterpolateJAC(daaA3,daToffset,daaaKx3,
     $        raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,1,daaT,-1)
      CALL SplineTempInterpolateJAC(daaA4,daToffset,daaaKx4,
     $        raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,1,daaT,-1)
      CALL SplineTempInterpolateJAC(daaA5,daToffset,daaaKx5,
     $        raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $        raQAirs,1,daaT,-1)

c then interpolate the five matrices in water partial pressure to get the 
c Compressed Absorption Matrix for water
c   daaKpro will have the daaAn interpolated in water amount
      CALL WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5,
     $    raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ)
c multiply daaUx with daaKpro to get daaAbsCoeff
ccccc user supplied info
c this is the assembly language matrix multiplication
c Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,
     $       iLDB,dBeta,daaAbsCoeff,iLDC)

      IF (kJacobian .GT. 0) THEN !do temperature jacobians
        CALL WaterAmount_TempJAC(daaaTemp,daaaKx1,daaaKx2,
     $         daaaKx3,daaaKx4,daaaKx5,raPPart,raRPart,               
     $         iNk,iKm,iKn,iUm,iUn)
c recall this subroutine does the uncompression of daaDT!!!
c as the last parameter (iDoDQ)=1, daaDQ is unchanged from its initial value
        CALL GetAbsCoeffJAC(daak_from_Temp,daToffset,daaaTemp,daaUx,
     $         raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $         daaDQ,daaDT,1,raQAirs,1) 
        IF ((kJacobian .GT. 0) .AND. (iDoDQ .GT. 0)) THEN
          CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)
          CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaQ,
     $       iLDB,dBeta,daaDQ,iLDC)
          END IF
        END IF

      RETURN
      END
c************************************************************************
c this interpolate the five matrices in partial press to get the Compressed
c Absorption Matrix for water, to give the final relevant daaKpro
      SUBROUTINE WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5,
     $  raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ)

      include 'kcarta.param'

c daaA1..5   = the matrices that will be interpolated in partial pressure
c raP/Rpart  = actual/rerefence water partial pressures
c daaKpro    = final resulting matrix for water amounts q
c daaQ     = final resulting matrix for water amounts q+dq
c see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      DOUBLE PRECISION daaA1(kMaxK,kProfLayer),daaA2(kMaxK,kProfLayer),
     $                 daaA3(kMaxK,kProfLayer),daaA4(kMaxK,kProfLayer),
     $                 daaA5(kMaxK,kProfLayer),daaKPro(kMaxK,kProfLayer)
      DOUBLE PRECISION daaQ(kMaxK,kProfLayer)
      INTEGER iNk,iKm,iKn,iUm,iUn

C     for interpolating
      DOUBLE PRECISION daWork(kMaxWater),dYP1,dYPN,dXPT,
     $    daXgiven(kMaxWater),daYgiven(kMaxWater),daY2(kMaxWater)
      DOUBLE PRECISION FirstDeriv
      INTEGER iI,iL

C     Assign some values for interpolation of K vectors
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

C     Do the spline Interpolation of the K vectors
C     Loop over the K singular vectors
      DO iI=1,iNk
C       Loop over the layers
        DO iL=1,kProfLayer
C         Interpolate across kMaxWater for the profile amount
          daXgiven(1)=0.1*raRPart(iL)
          daXgiven(2)=1.0*raRPart(iL)
          daXgiven(3)=3.3*raRPart(iL)
          daXgiven(4)=6.7*raRPart(iL)
          daXgiven(5)=10.0*raRPart(iL)

          daYgiven(1)=daaA1(iI,iL)
          daYgiven(2)=daaA2(iI,iL)
          daYgiven(3)=daaA3(iI,iL)
          daYgiven(4)=daaA4(iI,iL)
          daYgiven(5)=daaA5(iI,iL)          
          CALL DSPLY2(daXgiven,daYgiven,kMaxWater,dYP1,dYPN,daY2,daWork)

c directly take the Part Press amount in the actual profile as the X point
          dXPT=raPPart(iL)
          IF (dXPT .LT. daXgiven(1)) THEN
            dXPT=daXgiven(1)
            END IF
          IF (dXPT .GT. daXgiven(5)) THEN
            dXPT=daXgiven(5)
            END IF
          CALL DSPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,
     $                daaKpro(iI,iL))

          IF (kJacobian .GT. 0) THEN
            daaQ(iI,iL)=FirstDeriv(dXPT,daXgiven,daYgiven,daY2,
     $                             kMaxWater)
            END IF

          ENDDO
        ENDDO

      RETURN
      END 
c************************************************************************
c this subroutine performs interpolations in water amount, to produce 11
c matrices that will finally be interpolated in temperature, in 
c subroutine GetAbsCoeffJAC, to find d/dT

c !!!!!!! notice this is SIMILAR to splinetempinterp -- there we 
c         spline interpolate in pressure, and then spline interp in temp
c    daaaKx(kMaxK,kMaxTemp,kMaxLayer) --> daaaKx(kMaxK,kMaxTemp,kProfLayer) -->
c                                         daaKPRO(kMaxK,kProfLayer)
c here we take profiles raPPart(kProLayer),raRPart(kProfLayer) and then
c take matrices daaaKx1,daaaKx2,...,daaaKx5(kMaxK,kMaxTemp,kMaxLayer)
c and integrate the sublayers to produce daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
c this matrix will then be used to find d/dT
      SUBROUTINE WaterAmount_TempJAC(daaaTemp,daaaKx1,daaaKx2,
     $         daaaKx3,daaaKx4,daaaKx5,raPPart,raRPart,               
     $         iNk,iKm,iKn,iUm,iUn)

      include 'kcarta.param'
      include 'NewRefProfiles/outincLAY.param'

c daaaTemp    = final temperature interpolated matrix
c daaaKx1..5  = the 11 compressed matrices that will be interpolated in amount
c raP/RPart   = actual/reference amount profiles
c iNlay       = number of layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      DOUBLE PRECISION daaaTemp(kMaxk,kMaxTemp,kMaxLayer),  !!!NOT kProfLayer
     $    daaaKx1(kMaxK,kMaxTemp,kMaxLayer),
     $    daaaKx2(kMaxK,kMaxTemp,kMaxLayer),
     $    daaaKx3(kMaxK,kMaxTemp,kMaxLayer),
     $    daaaKx4(kMaxK,kMaxTemp,kMaxLayer),
     $    daaaKx5(kMaxK,kMaxTemp,kMaxLayer)
      INTEGER iNk,iKm,iKn,iUm,iUn

C     for interpolating KX
      DOUBLE PRECISION daWork(kMaxTemp),dYP1,dYPN,dXPT,
     $    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)
      
      INTEGER iI,iJ,iL,iLMin,iLMax,loffset,iM

C     Assign some values for interpolation of K vectors
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

c--------------------
c loop over the BOTTOM layers
      iLMin=1
      iLMax=M1000mb-1
      loffset=0
      DO iL=iLMin,iLMax
        daXgiven(1)=0.1*raRPart(iL)
        daXgiven(2)=1.0*raRPart(iL)
        daXgiven(3)=3.3*raRPart(iL)
        daXgiven(4)=6.6*raRPart(iL)
        daXgiven(5)=10.0*raRPart(iL)
        dXPT=raPPart(iL)

c loop over the  iNk singular vectors
        DO iI=1,iNk
c loop over iKm temperature matrices
          DO iJ=1,iKm
            daYgiven(1)=daaaKx1(iI,iJ,iL)
            daYgiven(2)=daaaKx2(iI,iJ,iL)
            daYgiven(3)=daaaKx3(iI,iJ,iL)
            daYgiven(4)=daaaKx4(iI,iJ,iL)
            daYgiven(5)=daaaKx5(iI,iJ,iL)

C  Interpolate KX across the 5 temperature matrices for the profile amount
            CALL DSPLY2(daXgiven,daYgiven,kMaxWater,
     $                  dYP1,dYPN,daY2,daWork)
            CALL DSPLIN(daXgiven,daYgiven,daY2,kMaxWater,dXPT,
     $                  daaaTemp(iI,iJ,iL))
            END DO
          ENDDO
        ENDDO
c--------------------
c loop over the MIDDLE layers
      iLMin=M1000mb
      iLMax=M100mb-1

      DO iL=iLMin,iLMax
        daXgiven(1)=0.0
        daXgiven(2)=0.0
        daXgiven(3)=0.0
        daXgiven(4)=0.0
        daXgiven(5)=0.0
        dXPT=0.0

        loffset=(iL-M1000mb)*MSubLayer + M1000mb 
        DO iM=0,MSubLayer-1
          daXgiven(1)=daXgiven(1)+0.1*raRPart(iM+loffset)
          daXgiven(2)=daXgiven(2)+1.0*raRPart(iM+loffset)
          daXgiven(3)=daXgiven(3)+3.3*raRPart(iM+loffset)
          daXgiven(4)=daXgiven(4)+6.6*raRPart(iM+loffset)
          daXgiven(5)=daXgiven(5)+10.0*raRPart(iM+loffset)
          dXPT=dXPT+raPPart(iM+loffset)
          END DO

c loop over the  iNk singular vectors
        DO iI=1,iNk
c loop over iKm temperature matrices ...... NOTICE NO OFFSET!!!!!!
          DO iJ=1,iKm
            daYgiven(1)=daaaKx1(iI,iJ,iL)
            daYgiven(2)=daaaKx2(iI,iJ,iL)
            daYgiven(3)=daaaKx3(iI,iJ,iL)
            daYgiven(4)=daaaKx4(iI,iJ,iL)
            daYgiven(5)=daaaKx5(iI,iJ,iL)

C  Interpolate KX across the 5 temperature matrices for the profile amount
            CALL DSPLY2(daXgiven,daYgiven,kMaxWater,
     $                  dYP1,dYPN,daY2,daWork)
            CALL DSPLIN(daXgiven,daYgiven,daY2,kMaxWater,dXPT,
     $                  daaaTemp(iI,iJ,iL))
            END DO
          ENDDO
        ENDDO
c--------------------
c loop over the TOP layers
      iLMin=M100mb
      iLMax=kMaxLayer
      loffset=(M100mb-M1000mb)*MSubLayer + (M1000mb-1) + 1
      !notice if MSubLayer=1,loffset=M100mb
      DO iL=iLMin,iLMax
        daXgiven(1)=0.1*raRPart((iL-iLMin)+loffset)
        daXgiven(2)=1.0*raRPart((iL-iLMin)+loffset)
        daXgiven(3)=3.3*raRPart((iL-iLMin)+loffset)
        daXgiven(4)=6.6*raRPart((iL-iLMin)+loffset)
        daXgiven(5)=10.0*raRPart((iL-iLMin)+loffset)
        dXPT=raPPart((iL-iLMin)+loffset)

c loop over the  iNk singular vectors
        DO iI=1,iNk
c loop over iKm temperature matrices ...... NOTICE NO OFFSET!!!!!!
          DO iJ=1,iKm
            daYgiven(1)=daaaKx1(iI,iJ,iL)
            daYgiven(2)=daaaKx2(iI,iJ,iL)
            daYgiven(3)=daaaKx3(iI,iJ,iL)
            daYgiven(4)=daaaKx4(iI,iJ,iL)
            daYgiven(5)=daaaKx5(iI,iJ,iL)

C  Interpolate KX across the 5 temperature matrices for the profile amount
            CALL DSPLY2(daXgiven,daYgiven,kMaxWater,
     $                  dYP1,dYPN,daY2,daWork)
            CALL DSPLIN(daXgiven,daYgiven,daY2,kMaxWater,dXPT,
     $                  daaaTemp(iI,iJ,iL))
            END DO
          ENDDO
        ENDDO
c--------------------

      RETURN
      END

c************************************************************************
c this subroutine calculates the analytic approximation to the first 
c derivative, using spline approximations
      DOUBLE PRECISION FUNCTION FirstDeriv(x,daXgiven,daYgiven,daY2,iKm)

      include 'kcarta.param'

c daXgiven   == x ordinates (known)
c daYgiven   == y(x)        (known)
c daY2       == d2y/dx2     (known)
c x          == value where we need to find dy/dx at
c iKm        == number of set data points daXgiven(1..iKm)
      DOUBLE PRECISION daXgiven(kMaxTemp),daYgiven(kMaxTemp),
     $                 daY2(kMaxTemp),x
      INTEGER iKm

c local variables
      INTEGER iK
      DOUBLE PRECISION a,b,dAns
      INTEGER KLO,KHI

C     Determine between which pair of points X falls (bisect loop)
      KLO=1
      KHI=iKm
 20   IF ( (KHI - KLO) .GT. 1) THEN
         iK=(KHI + KLO)/2
         IF (daXgiven(iK) .GT. X) THEN
            KHI=iK
         ELSE
            KLO=iK
         ENDIF
         GOTO 20
      ENDIF

      dAns=0.0

      IF (iK .LT. iKm) THEN
        a=(daXgiven(iK+1)-x)/(daXgiven(iK+1)-daXgiven(iK))
        b=1-a

        dAns=(daYgiven(iK+1)-daYgiven(iK))/(daXgiven(iK+1)-daXgiven(iK))
     $     -(3.0*a*a-1.0)/6.0*(daXgiven(iK+1)-daXgiven(iK))*daY2(iK)
     $     +(3.0*b*b-1.0)/6.0*(daXgiven(iK+1)-daXgiven(iK))*daY2(iK+1)
        END IF

      FirstDeriv=dAns

      RETURN
      END

c************************************************************************
