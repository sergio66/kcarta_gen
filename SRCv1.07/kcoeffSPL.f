c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c this subroutine initialises the parameters for call to dGEMM
      SUBROUTINE InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)

      include 'kcarta.param'

      INTEGER iLDA,iLDB,iLDC,iM,iN,iK,iUm,iUn
      DOUBLE PRECISION dAlpha,dBeta
      CHARACTER*1 cTRANSA,cTRANSB

      cTRANSA='N'
      cTRANSB='N'
      iM=iUm
      iN=kProfLayer      !used to be iN=iKn
      iK=iUn
      dAlpha=1.0
      iLDA=kMaxPts
      iLDB=kMaxK
      dBeta=0.0
      iLDC=kMaxPts

      RETURN
      END
c************************************************************************
c this calls subroutine to do pressure, temperature interpolations, then
c does the uncompression
c this is for gases other than water
      SUBROUTINE GetAbsCoeffOLD(daaAbsCoeff,daToffset,daaaKx,daaUx,
     $  raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  raQAirs,iGasID)

      include 'kcarta.param'
      include 'NewRefProfiles/outincLAY.param'

c daaAbsCoeff = abs coeff matrix, after the temperature interpolations
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx      = the 11 compressed matrices that will be interpolated
c daaUx       = the uncompression matrix
c raP/Rtemp   = actual/reference temperature profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = AIRS number of layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNLay x iKm) === iNk x (100 x 11) --------> careful!!
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx
c raQAirs   = ratio of kLAYERS/KAIRS reference gas amounts  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk

      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),raQAirs(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer),
     $    daaUx(kMaxPts,kMaxK)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC,iM,iN,iK
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB

      CALL SplineTempInterpolateOLD(daaKpro,daToffset,daaaKx,raPTemp,
     $      raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,iGasID)

c multiply daaUx with daaKpro to get daaAbsCoeff
ccccc user supplied info
c this is the assembly language matrix multiplication
c  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,
     $       iLDB,dBeta,daaAbsCoeff,iLDC)

      RETURN
      END
c************************************************************************
c interpolate the compressed data in temperature AND water amount
      SUBROUTINE GetAbsCoeffWaterOLD(daaAbsCoeff,daToffset,
     $  daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,
     $  raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,
     $  raQAirs,iGasID)

      include 'kcarta.param'

c daaAbsCoeff = abs coeff matrix, after the temperature interpolations
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx1..5  = the 11 compressed matrices that will be interpolated in temp
c  after which the water partial pressure * 0.01,1.0,3.3,6.7,10.0 interpolation
c daaUx       = the uncompression matrix
c raP/Rtemp   = actual/reference temperature profiles
c raP/Rpart   = actual/reference partial pressure profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = number of AIRS layers (=kMaxLayer)
c iNk         = number of singular vectors
c iKm         = number of temperature matrices (=11)
c iKn         = number of layers in k-comp matrices 
c   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
c iUm         = number of freq points in daaUx
c iUn         = number of singular vectors in daaUx  
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
c raQAirs   = ratio of kLAYERS/KAIRS reference gas amounts  
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),raQAirs(kProfLayer)
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK),
     $  daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer),
     $  daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer),
     $  daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer),
     $  daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer),
     $  daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB

      INTEGER iM,iN,iK

c first interpolate in each of the five water offset matrices, in temperature
c to obtain the five current temp profile water matrices
c hence daaAn = the n th matrix, interpolated in layer temperature T
      CALL SplineTempInterpolateOLD(daaA1,daToffset,daaaKx1,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,1)
      CALL SplineTempInterpolateOLD(daaA2,daToffset,daaaKx2,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,1)
      CALL SplineTempInterpolateOLD(daaA3,daToffset,daaaKx3,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,1)
      CALL SplineTempInterpolateOLD(daaA4,daToffset,daaaKx4,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,1)
      CALL SplineTempInterpolateOLD(daaA5,daToffset,daaaKx5,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,1)

c then interpolate the five matrices in water partial pressure to get the 
c Compressed Absorption Matrix for water
c daaKpro will have the daaAn interpolated in water amount
      CALL WaterAmountInterpolateOLD(daaA1,daaA2,daaA3,daaA4,daaA5,
     $     raQAirs,raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn)

c multiply daaUx with daaKpro to get daaAbsCoeff
ccccc user supplied info
c this is the assembly language matrix multiplication
c Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,
     $       iLDB,dBeta,daaAbsCoeff,iLDC)

      RETURN
      END
c************************************************************************
c this subroutine does the interpolation of a compressed matrix
c first it does a pressure interpolation
c   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
c then temperature interpolation
c   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
      SUBROUTINE SplineTempInterpolateOLD(daaKpro,daToffset,daaaKx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,raQAirs,iGasID)

      include 'kcarta.param'
      include 'NewRefProfiles/outincLAY.param'

c daakPro     = final temperature interpolated matrix
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx      = the 11 compressed matrices that will be interpolated in temp
c raP/Rtemp   = actual/reference temperature profiles
c iaTsort     = integer indices of temperature offsets
c iNlay       = number of AIRS layers (=kMaxLayer)
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
          DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
            daYgiven(iJ)=daaaKxNew(iI, iaTsort(iJ), iK)
            ENDDO
          CALL DSPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
c subtract the ref temp from the profile temperature
          dXPT=raPtemp(iK) - raRTemp(iK)
          CALL DSPLIN(daXgiven,daYgiven,daY2,iKm,dXPT,daaKpro(iI,iK))
          ENDDO
        ENDDO

      RETURN
      END 
c************************************************************************
c this interpolate the five matrices in partial press to get the Compressed
c Absorption Matrix for water, to give the final relevant daaKpro
      SUBROUTINE WaterAmountInterpolateOLD(
     $             daaA1,daaA2,daaA3,daaA4,daaA5,raQAirs,
     $             raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn)

      include 'kcarta.param'

c daaA1..5   = the matrices that will be interpolated in partial pressure
c daaT1..5   = the (T+dT)matrices that will be interpolated in part pressure
c raP/Rpart  = actual/rerefence water partial pressures
c daaKpro    = final resulting matrix for water amounts q
c daaQ     = final resulting matrix for water amounts q+dq
c see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
      REAL raPPart(kProfLayer),raRPart(kProfLayer),raQAirs(kProfLayer)
      DOUBLE PRECISION daaA1(kMaxK,kProfLayer),daaA2(kMaxK,kProfLayer),
     $                 daaA3(kMaxK,kProfLayer),daaA4(kMaxK,kProfLayer),
     $                 daaA5(kMaxK,kProfLayer),daaKPro(kMaxK,kProfLayer)
      INTEGER iNk,iKm,iKn,iUm,iUn

C     for interpolating
      DOUBLE PRECISION daWork(kMaxWater),dYP1,dYPN,dXPT,
     $    daXgiven(kMaxWater),daYgiven(kMaxWater),daY2(kMaxWater),d

      INTEGER iI,iK

C     Assign some values for interpolation of K vectors
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

C     Do the spline Interpolation of the K vectors
C     Loop over the K singular vectors
      DO iI=1,iNk
C       Loop over the layers
        DO iK=1,kProfLayer
C         Interpolate across kMaxWater for the profile amount
          daXgiven(1)=0.1*raRPart(iK)
          daXgiven(2)=1.0*raRPart(iK)
          daXgiven(3)=3.3*raRPart(iK)
          daXgiven(4)=6.7*raRPart(iK)
          daXgiven(5)=10.0*raRPart(iK)

          daYgiven(1)=daaA1(iI,iK)
          daYgiven(2)=daaA2(iI,iK)
          daYgiven(3)=daaA3(iI,iK)
          daYgiven(4)=daaA4(iI,iK)
          daYgiven(5)=daaA5(iI,iK)          
          CALL DSPLY2(daXgiven,daYgiven,kMaxWater,dYP1,dYPN,daY2,daWork)

c directly take the Part Press amount in the actual profile as the X point
          dXPT=raPPart(iK)
          IF (dXPT .LT. daXgiven(1)) THEN
            dXPT=daXgiven(1)
            END IF
          IF (dXPT .GT. daXgiven(5)) THEN
            dXPT=daXgiven(5)
            END IF
          CALL DSPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,d)
          daaKpro(iI,iK)=d

          ENDDO
        ENDDO

      RETURN
      END 

c************************************************************************
