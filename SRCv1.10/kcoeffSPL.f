c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c this subroutine initialises the parameters for call to dGEMM
      SUBROUTINE InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                     iLDA,iLDB,iLDC,dbeta,iUm,iUn)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
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
c   daaUx=iUm x iUn  = 10000 x iUn
c   WITH iUn === iNk
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer),
     $    daaUx(kMaxPts,kMaxK)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC,iM,iN,iK
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB

      CALL SplineTempInterpolateOLD(daaKpro,daToffset,daaaKx,raPTemp,
     $      raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers)

c multiply daaUx with daaKpro to get daaAbsCoeff
ccccc user supplied info
c this is the assembly language matrix multiplication
c  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,
     $                    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,
     $       iLDB,dBeta,daaAbsCoeff,iLDC)

      RETURN
      END
c************************************************************************
c interpolate the compressed data in temperature AND water amount
      SUBROUTINE GetAbsCoeffWaterOLD(daaAbsCoeff,daToffset,
     $  daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,
     $  raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $  pProf,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c pProf       = actual layers from kLAYERS avg pressure
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
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
      REAL raPPart(kProfLayer),raRPart(kProfLayer)
      DOUBLE PRECISION daaAbsCoeff(kMaxPts,kProfLayer),
     $    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK),
     $  daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer),
     $  daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer),
     $  daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer),
     $  daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer),
     $  daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers

C     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER iLDA,iLDB,iLDC
      DOUBLE PRECISION dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER*1 cTRANSA,cTRANSB

      INTEGER iM,iN,iK

c first interpolate in each of the five water offset matrices, in temperature
c to obtain the five current temp profile water matrices
c hence daaAn = the n th matrix, interpolated in layer temperature T
      CALL SplineTempInterpolateOLD(daaA1,daToffset,daaaKx1,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,
     $   pProf,iProfileLayers)
      CALL SplineTempInterpolateOLD(daaA2,daToffset,daaaKx2,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,
     $   pProf,iProfileLayers)
      CALL SplineTempInterpolateOLD(daaA3,daToffset,daaaKx3,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,
     $   pProf,iProfileLayers)
      CALL SplineTempInterpolateOLD(daaA4,daToffset,daaaKx4,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,
     $   pProf,iProfileLayers)
      CALL SplineTempInterpolateOLD(daaA5,daToffset,daaaKx5,
     $   raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,
     $   pProf,iProfileLayers)

c then interpolate the five matrices in water partial pressure to get the 
c Compressed Absorption Matrix for water
c daaKpro will have the daaAn interpolated in water amount
      CALL WaterAmountInterpolateOLD(daaA1,daaA2,daaA3,daaA4,daaA5,
     $                  raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,
     $                  pProf,iProfileLayers)

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
c note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
c first it does a pressure interpolation
c   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
c then temperature interpolation
c   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
      SUBROUTINE SplineTempInterpolateOLD(daaKpro,daToffset,daaaKx,
     $      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,
     $      pProf,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      INTEGER iplev 
      include '../INCLUDE/KCARTA_database.param'

c pProf       = actual layers from kLAYERS avg pressure
c daakPro     = final temperature interpolated matrix, in ABS COEFF
c daToffset   = temperature offsets that the compressed matrices were made at
c daaaKx      = the 11 compressed matrices that will be interpolated in temp
c               NOTICE THESE HAVE 
c                  (OPTICAL DEPTHS)^1/4 = (ref amt * abs coeff)^(1/4)
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
      REAL raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
      DOUBLE PRECISION daaKpro(kMaxk,kProfLayer),
     $    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer)
      INTEGER iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers

C     for interpolating daaaKX in temperature
      DOUBLE PRECISION daWork(kMaxTemp),dYP1,dYPN,dXPT,
     $    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)

C     for interpolating daaaKX in pressure
      DOUBLE PRECISION daWorkP(kMaxLayer),
     $    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)

C this is the actual matrix that will be created after interpolation in 
c pressure, and then interpolated in temperature
      DOUBLE PRECISION daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
      INTEGER iI,iJ,iK,iM,loffset,loffset2,iLowest

c these are to read in the original 100 layers AIRS ref profile for the gas
c these are to read in the new kProfLayers ref profile for the gas
      CHARACTER*80 caFName
      INTEGER iE
      REAL raOrig100A(kMaxLayer),raOrig100T(kMaxLayer) 
      REAL raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer) 

      !read in the orig 100 layer prof
      write (kStdWarn,*) 'Temperature interpolation (5 times for water, once 
     $  for other gases)'
      write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
      write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
      CALL FindReferenceName(caFName,iGasID,-1)  
      CALL ReadRefProf(caFName,kProfileUnit,kMaxLayer,raOrig100A,raOrig100T,
     $         raOrig100P,raOrig100PP,iE)

C     Assign values for interpolation
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

      iLowest = kProfLayer - iProfileLayers + 1

      DO iI=1,iKm
        daXgiven(iI)=daToffset(iaTsort(iI))
        ENDDO

c  even if kProfLayer =========== kMaxLayer, allow for possibility that
c  user has changed layering, so we have to do PRESSURE interpolation
c  of daaaKx onto daaaKnNew

      DO iI=1,iNk
        DO iJ=1,iKm
          !first set up the daXgivenP,daYgivenP arrays
          DO iK=1,kMaxLayer   
            iE = kMaxLayer-iK+1
            !notice how daXgiven is initialised with increasing pressure
            !do interpolate in log(pressure)
            daXgivenP(iK)=log(pavg_kcartadatabase_AIRS(iE))*1.0d0
            !notice how doYgiven is normalised to daaaKx/(100layer ref amount)
            !this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
            !is being changed to raw (abs coeff)^(1/4)
            daYgivenP(iK)=daaaKx(iI,iJ,iE)/(raOrig100A(iE)**0.25)
            END DO
          CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN,
     $                    daY2P,daWorkP)

c  do the new set of layers ..need AIRS layers interpolations
          DO iK=iLowest,kProfLayer
            dxpt=log(pProf(iK))*1.0d0
            CALL dsplin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d)
            daaaKxNew(iI,iJ,iK)=d
            END DO
          END DO
        END DO

C     now do the spline Interpolation of the K vectors in TEMPERATURE
      DO iI=1,iNk                         !Loop over the K vectors
        DO iK=iLowest,kProfLayer   !Loop over the layers
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
     $             daaA1,daaA2,daaA3,daaA4,daaA5,
     $             raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,
     $             pProf,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaA1..5   = the matrices that will be interpolated in partial pressure
c raP/Rpart  = actual/rerefence water partial pressures
c daaKpro    = final resulting matrix for water amounts q
c daaQ     = final resulting matrix for water amounts q+dq
c see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
      REAL raPPart(kProfLayer),raRPart(kProfLayer),pProf(KprofLayer)
      DOUBLE PRECISION daaA1(kMaxK,kProfLayer),daaA2(kMaxK,kProfLayer),
     $                 daaA3(kMaxK,kProfLayer),daaA4(kMaxK,kProfLayer),
     $                 daaA5(kMaxK,kProfLayer),daaKPro(kMaxK,kProfLayer)
      INTEGER iNk,iKm,iKn,iUm,iUn,iProfileLayers

C     for interpolating
      DOUBLE PRECISION daWork(kMaxWater),dYP1,dYPN,dXPT,
     $    daXgiven(kMaxWater),daYgiven(kMaxWater),daY2(kMaxWater),d

      INTEGER iI,iK,iLowest

      iLowest = kProfLayer - iProfileLayers + 1

C     Assign some values for interpolation of K vectors
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

C     Do the spline Interpolation of the K vectors
C     Loop over the K singular vectors
      DO iI=1,iNk
C       Loop over the layers
        DO iK=iLowest,kProfLayer
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
