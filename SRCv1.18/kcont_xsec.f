c Copyright 2011
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c********* this file has the continuum/old xsec routines ****************
c************************************************************************
c this adds on the continuum : 
      SUBROUTINE AddContinuum(iGasID,iTag,iActualTag,rFileStartFr,
     $                   iRefLayer,iProfileLayers,raFreq,raAmt,raTemp,
     $                   kFrStep,raPress,raPartPress,
     $                   iL,iU,daaTemp,daaDQ,daaDT,iDoDQ,iSPlineType,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,pProf,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iGasID    = iaGasID(iCount) = gas ID of current gas
c iRefLayer = number of layers in the reference profiles (=kProfLayer)
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
c daaTemp   = matrix containing the uncompressed k-spectra
c raFreq    = wavenumber array
c iErr      = errors (mainly associated with file I/O)
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c kFrStep  = kFreqStep
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)
      INTEGER iGasID,iL,iU,iErr,iRefLayer,iDoDQ
      INTEGER iTag,iActualTag,iSplineType,iProfileLayers
      REAL raFreq(kMaxPts),kFrStep,pProf(kProfLayer),rFileStartFr
      REAL raAmt(kProfLayer),raTemp(kProfLayer)
      REAL raPartPress(kProfLayer),raPress(kProfLayer)
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c the Matlab weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer),
     $        iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer),
     $        raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer),
     $        raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer),
     $        iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer),
     $        raQ21(kProfLayer),raQ22(kProfLayer)

c iMethod   = 0 for old style, 1 for new style
c local variables
      INTEGER iLay,iFr,iIOUN,iLoc
      REAL rFStep,rMin,rConvFac
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iMin,iMax,iErr1

      INTEGER iDefault,iUseMethod,iMethod

      iDefault = +1      !!! new method, uses kCompressed chunks
      iDefault = +0      !!! original method, uses 2901 (fr x T) lookup table

      iDefault = +0      !!! original method, uses 2901 (fr x T) lookup table
      iMethod = 0        !!! assume we need OLD method

      IF ((abs(iSplineType) .EQ. 2) .AND. (iActualTag .EQ. 20)) THEN
        !! warning so far only CKD1,6,25 has these tables ....
        iMethod = 1        !! use compressed data for 605-2830 cm-1
        iDefault = +1      !! new method, uses kCompressed chunks

        !decided to switch back to look up tables, as Matlab version uses this
        iMethod = 0        !! use look up tables
        iDefault = 0       !! use look up tables
      END IF

      IF (iDefault .NE. iMethod) THEN
        write(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
      END IF

      CALL ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)

      IF (iMethod .EQ. 0) THEN
        ! reads the oldstyle 2901 point lookup table which is fcn(fr,T)
        ! old style, been using this for years
        write (kStdWarn,*) ' WV cont : using Lookup Table (fr x T) format' 
        rConvFac = 1.0   !!! scale factor = 1 in the ORIG tables
        CALL f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ,
     $                         rFStep,raFreq,kFrStep,
     $                         raTemp,daaDQ,daaDT,daaCon)
      ELSEIF (iMethod .EQ. 1) THEN
        ! reads the Compressed Style lookup fcn(fr,T)
        ! new style, since Sept 2011
        write (kStdWarn,*) ' WV cont : using kCompressed Database format'
        rConvFac = 1.0e-23    !! SELF/FORN continuum coeffs in compressed 
                              !! database were remultiplied by 1.0e+23
        CALL xCKDgases(iGasID,rFileStartFr,iTag,iActualTag,
     $                 iRefLayer,iL,iU,
     $                 raAmt,raRAmt,raTemp,raRTemp,
     $                 iErr1,iDoDQ,pProf,iProfileLayers,
     $                 daaDQ,daaDT,daaCon,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
      END IF

      CALL aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt,
     $             raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c does some simple inits
      SUBROUTINE ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL raAmt(kProfLayer),kFrStep
c output
      INTEGER iMin,iMax,iLay,iFr
      REAL rFStep
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      
c local
      INTEGER iN

c everywhere we have         DO iLay=1,kProfLayer
c  we replace with           DO iLay=iMin,iMax 
c      iMin = 1
c      iMax = kProfLayer
      iMax = -1
      DO iN = 1,kProfLayer
        IF ((raAmt(iN) .GT. 0) .AND. (iN .GT. iMax)) THEN
          iMax = iN
          END IF
        END DO
      iMin = kProfLayer + 1
      DO iN = kProfLayer,1,-1
        IF ((raAmt(iN) .GT. 0) .AND. (iN .LT. iMin)) THEN
          iMin = iN
        END IF
      END DO

c initialize raaCon to all zeros
      DO iFr=1,kMaxPts
        DO iLay=iMin,iMax
          daaCon(iFr,iLay)=0.0
        END DO
      END DO

      iLay = kProfLayer
      iFr = kMaxPts
      rFStep = kFrStep

      RETURN
      END

c************************************************************************
      SUBROUTINE f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ,
     $                               rFStep,raFreq,kFrStep,
     $                               raTemp,daaDQ,daaDT,daaCon)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input  
      INTEGER iGasID,iTag,iDoDQ
      REAL raFreq(kMaxPts),rFStep,kFrStep
      REAL raTemp(kProfLayer)
c output
c   daaDQ     = analytic Jacobian wrt gas amount
c   daaDT     = analytic Jacobian wrt temperature
c   daaCon    = continuum coeffs (multiplied by 1)
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)       
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)

c local vars
      INTEGER iLay,iFr,iIOUN,iLoc,iDefault,iUseMethod,iMethod,iErr
      REAL rMin
      DOUBLE PRECISION der1,der2,doop
      REAL rEC,rECCount

c these are to read in the binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)
      INTEGER iCKD,iM,iN

      INTEGER iMin,iMax
      REAL a1,a2
      CHARACTER*120 caFName

c figure out the filename
      CALL CKDFileName(caFName,iGasID,iTag)
      
      iIOUN = kTempUnit 
c open file and load in data
      OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='UNFORMATTED',
     $    IOSTAT=iErr)
      IF (iErr .NE. 0) THEN
        WRITE(kStdErr,1080) iErr, caFname
        WRITE(kStdWarn,1080) iErr, caFname
 1080   FORMAT('ERROR! number ',I5,' opening CKD binary file:',/,A120)
        CALL DoSTOP
      ENDIF

      kTempUnitOpen=1
      READ(iIOUN) d1,d2,df       !read start/stop freq, df
      IF (d1 .GT. raFreq(1)) THEN
        write(kStdErr,*) 'CKD file has freqs that start too high!'
        CALL DoStop
      END IF
      IF (d2 .LT. raFreq(kMaxPts)) THEN
        write(kStdErr,*) 'CKD file has freqs end start too low!'
        CALL DoStop
      END IF

c&&&&&&&&
      READ(iIOUN) iLoc,iCKD,iM,iN
      !read line type, CKD vers, # of temp,freq pts
      IF (iLoc .LT. 0) THEN
        write(kStdErr,*) iLoc    
        write(kStdErr,*) 'need correct continuum for local lineshape !'
        CALL DoStop
      END IF
      IF (iCKD .NE. kCKD) THEN
        write(kStdErr,*) iCKD,kCKD
        write(kStdErr,*) 'CKD versions do not match!'
        CALL DoStop
      END IF
      IF (iM .NE. kTempCKD) THEN
        write(kStdErr,*) iM,kTempCKD
        write(kStdErr,*) 'Need more Temp offsets in CKD binary file'
        CALL DoStop
      END IF
      IF (iN .NE. kFreqCKD) THEN
        write(kStdErr,*) iN,kFreqCKD
        write(kStdErr,*) 'Read in ',iN, ' points from CKD binary file'
        write(kStdErr,*) 'Always need ',kFreqCKD, ' points in the file'
        CALL DoStop
      END IF

      READ(iIOUN) (daTemprt(iLay),iLay=1,iM)   !read the temps

      DO iLay=1,iM
        !for current temp, read the continuum data for each freq
        READ(iIOUN) (daaCKD(iLay,iFr),iFr=1,iN) 
      END DO

      CLOSE(iIOUN)
      kTempUnitOpen=-1
c&&&&&&&&
    
      iMethod = +1       !!! quadratic
      iMethod = +2       !!! linear, T dependance in Cs,Cf
      iMethod = +3       !!! linear, T dependance in Cs only

      iDefault = +1      !!! used by Scott in making Fast Models
      iDefault = +3      !!! been used in KCARTA since 2005

      IF (iDefault .NE. iMethod) THEN
        write(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
      END IF

c interpolate
      IF (iMethod .EQ. 1) THEN
        ! this is what was used in SRCv1.07, which is what Scott Hannon used
        ! for the Fast Models
        CALL ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,
     $                 kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

      ELSEIF (iMethod .EQ. 2) THEN
        ! this is new and improved version, which should be quicker than Quad
        ! THIS ONE has T dependnce both for self and foreign
        CALL ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,
     $                 kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

      ELSEIF (iMethod .EQ. 3) THEN
        ! this is new and improved version, which should be quicker than Quad
        ! THIS ONE has T dependance ONLY for self
        CALL ComputeCKD_Linear_March2002(d1,d2,df,
     $                   daTemprt,daaCKD,raFreq,raTemp,
     $                   kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)
      END IF

      RETURN
      END

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c this subroutine completes the CKD calcs, and d/dT,d/dQ if needed
      SUBROUTINE aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt,
     $             raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iMin,Imax tell over what layers to do the CKD calcs (where gasamt nonzero!)
c input
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)       !! ckd coeffs
      INTEGER iMin,iMax,iGasID,iDoDQ
      REAL raAmt(kProfLayer),raTemp(kProfLayer)
      REAL raPartPress(kProfLayer),raPress(kProfLayer)
      REAL raFreq(kMaxPts),kFrStep,rConvFac
c output
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer)      !! ckd ODs
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)

c local variables
      INTEGER iLay,iFr,iErr
      REAL rFStep,rMin
      DOUBLE PRECISION der1,der2,doop
      REAL rEC,rECCount
      REAL raMult(kProfLayer)
      REAL a1,a2

      rMin=1.0e10 
      rECCount=0.0 
      DO iLay=iMin,iMax
        DO iFr=1,kMaxPts 
          IF (daaCon(iFr,iLay) .lt. 0.0) THEN
c            write(kStdWarn,*) 'continuum < 0 for iGasID,ilay,iFr',
c     $                         iGasID,iLay,iFr,daaCon(iFr,iLay) 
            daaCon(iFr,iLay) = 0.0
            rECCount = rECCount + 1.0 
            IF (daaCon(iFr,iLay) .lt. rMin) THEN
              rMin = daaCon(iFr,iLay) 
            END IF 
          END IF
        END DO
      END DO
      IF (rECCount .GT. 0) THEN 
        write(kStdWarn,*) ' >>> WARNING !!! Continuum for iGasID ',iGasID,' < 0 for ',int(rECCount),' points'
      END IF

c add the continuum abs coeffs to the k-compressed coeffs raaCon = daaTemp
c keeping track of whether the added values are all greater than 0.0
      IF (iGasID .EQ. kNewGasLo) THEN      !CSelf
        DO iLay =iMin,iMax
          raMult(iLay)=raAmt(iLay)*kAvog*raPartPress(iLay)
        END DO
      ELSEIF (iGasID .EQ. kNewGasHi) THEN      !CFor
        DO iLay =iMin,iMax
          raMult(iLay)=raAmt(iLay)*kAvog*(raPress(iLay)-raPartPress(iLay))
        END DO
      END IF

      DO iLay =iMin,iMax
        a1 = raMult(iLay)*296.0/raTemp(iLay)
        a2 = kPlanck2/2/raTemp(iLay)
        DO iFr=1,kMaxPts
          daaTemp(iFr,iLay) = rConvFac * daaCon(iFr,iLay)*a1*raFreq(iFr)*
     $                               tanh(a2*raFreq(iFr))
        END DO
      END DO

      IF (rECCount .gt. 0.5) THEN
        rEC=rEC/rECCount
        write(kStdWarn,*)'Error in CKD data!!! Some values negative!!!'
        write(kStdWarn,*)'and reset to 0.0'
        write(kStdWarn,*) rECCount,' values have avg value ',rEC
        write(kStdWarn,*) 'min negative value in 10000*100 = ',rMin
        IF (abs(rMin) .GT. 1.0e-7) THEN
          iErr=1
          write(kStdErr,*)'Error in CKD data!!! Some values negative!!!'
          CALL DoSTOP
        END IF
      END IF          

      IF (iDoDQ .GE. -1) THEN
        IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 20)) THEN

          !!the gas amount jacobians
          IF (iGasID .EQ. 101) THEN   !!!there is a factor of 2
            DO iLay=iMin,iMax
              DO iFr=1,kMaxPts
                daaDQ(iFr,iLay) = 2 * daaTemp(iFr,iLay)/raAmt(iLay) 
              END DO
            END DO
          ELSEIF (iGasID .EQ. 102) THEN
            DO iLay=iMin,iMax
              DO iFr=1,kMaxPts
                daaDQ(iFr,iLay) = daaTemp(iFr,iLay)/raAmt(iLay) 
              END DO
            END DO
          END IF
        END IF
 
        IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 30). OR. 
     $         (kActualJacs .EQ. 32) .OR. 
     $         (kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102)) THEN
          !!this is the temperature jacobian
          DO iLay=iMin,iMax
            DO iFr=1,kMaxPts
              doop = tanh(kPlanck2*raFreq(iFr)/2/raTemp(iLay))
              der1 = -(296/raTemp(iLay)**2)* doop
              der1 = der1-(296*kPlanck2*raFreq(iFr)/2/(raTemp(iLay)**3))/
     $           (cosh(kPlanck2*raFreq(iFr)/2/raTemp(iLay))**2)
              der1 = der1*daaCon(iFr,iLay)

              der2 = (296/raTemp(iLay))* doop
              der2 = der2*daaDT(iFr,iLay)

              daaDT(iFr,iLay) = raMult(iLay)*raFreq(iFr)*(der1+der2)*rConvFac
            END DO
          END DO
        END IF

      END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, linear in temperature
c does temperature interpolation for SELF and for foreign
c as this is the new CKD UMBC thingies
      SUBROUTINE ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep   = kFreqStep
c iMin,iMax = layers over which to loop
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID,iMin,iMax
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2
c this is for CKD1 temeperature dependant multipliers from 780-980 cm-1
      INTEGER iNpts,iMult
      DOUBLE PRECISION psT(kMaxLayer),psK(kMaxLayer)
      DOUBLE PRECISION pfT(kMaxLayer),pfK(kMaxLayer)
      DOUBLE PRECISION xNum(kMaxLayer)
      DOUBLE PRECISION daaCKD1Mult(kMaxPts,kProfLayer)

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c      iMult = -1
c      IF ((kCKD .EQ. 1) .AND. 
c     $  (raFreq(1) .LT. 979.0) .AND. (raFreq(kMaxPts) .GE. 780)) THEN
c        iMult = +1
c        write(kStdWarn,*) 'CKD 1 Temperature correction for chunk  ',raFreq(1)
c        CALL GetCKD1Mult(iGasID,raTemp,raFreq,
c     $                   iNpts,xNum,psT,psK,pfT,pfK,daaCKD1Mult)
c        END IF

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
      END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=iMin,iMax
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
        END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)

              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
 
            END DO
          END DO
        END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)

              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!
  
            END DO
          END DO
        END IF
      END IF

      IF (iMult .GT. 0) THEN
        DO iLay = iMin,iMax
          DO iFr = 1,kMaxPts
            daaCon(iFr,iLay) = daaCon(iFr,iLay) * daaCKD1Mult(iFr,iLay)
          END DO
        END DO
      END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, linear in temperature
c only does temperature interpolation for SELF, and not for foreign
c as this was the CKD 0,2.1,2.3,2.4 thingies
      SUBROUTINE ComputeCKD_Linear_March2002(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep  = kFreqStep
c iMin,iMax = layers over which to loop
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID,iMin,iMax
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
      END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=iMin,iMax
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
      END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find out line that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
  
            END DO
          END DO
        END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)

              a = (y3-y2)/(x3-x2)
              b = y3-a*x3

              x = raTemp(iLay)
              daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
              daaDT(iFr,iLay)  = a           !this is temp jacobian!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
              daaDT(iFr,iLay)=0.0d0  !this is temp jacobian!
  
            END DO
          END DO
        END IF
      END IF

      RETURN
      END
c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses interpolations : linear in freq, quadratic in temperature
      SUBROUTINE ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep   = kFreqStep
c iMin,iMax = layers over which to loop
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iGasID,iMin,iMax
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iLay,iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp
      DOUBLE PRECISION a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
      END DO   

c for each layer temp, find where the nearest low CKD tempr grid point is
      DO iLay=iMin,iMax
        iaTempIndex(iLay)=1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
      END DO   

      IF (iDoDQ .LT. -1) THEN         !no need to do temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y1=daaCKD(iL-1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL-1,iF+1)-daaCKD(iL-1,iF))/df
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
              x1 = daTemprt(iL-1)
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)
              z1=y1-y3
              z2=y2-y3
              t1=x1*x1-x3*x3
              t2=x1-x3
              t3=x2*x2-x3*x3
              t4=x2-x3
              b=(z1*t3-z2*t1)/(t3*t2-t1*t4)
              a=(z2-b*t4)/t3
              c=y3-a*x3*x3-b*x3

              x=raTemp(iLay)
              daaCon(iFr,iLay)=a*x*x + b*x + c     !this is temp dependance!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
  
            END DO
          END DO
        END IF

      ELSE          !!!!!!!!!do the temp jacobians
        IF (iGasID .EQ. kNewGasLo) THEN
          !self continuum has temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
 
              y1=daaCKD(iL-1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL-1,iF+1)-daaCKD(iL-1,iF))/df
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
              y3=daaCKD(iL+1,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL+1,iF+1)-daaCKD(iL+1,iF))/df
          
              !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
              x1 = daTemprt(iL-1)
              x2 = daTemprt(iL)
              x3 = daTemprt(iL+1)
              z1=y1-y3
              z2=y2-y3
              t1=x1*x1-x3*x3
              t2=x1-x3
              t3=x2*x2-x3*x3
              t4=x2-x3
              b=(z1*t3-z2*t1)/(t3*t2-t1*t4)
              a=(z2-b*t4)/t3
              c=y3-a*x3*x3-b*x3

              x=raTemp(iLay)
              daaCon(iFr,iLay)=a*x*x + b*x + c     !this is temp dependance!
              daaDT(iFr,iLay) =2*a*x + b           !this is temp jacobian!

            END DO
          END DO

        ELSEIF (iGasID .EQ. kNewGasHi) THEN
          !foreign continuum has no temp dependance
          DO iLay=iMin,iMax
            iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
            DO iFr=1,kMaxPts
              iF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
 
              y2=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
          
              daaCon(iFr,iLay)=y2    !this is temp dependance!
              daaDT(iFr,iLay)=0.0d0  !this is temp jacobian!
  
            END DO
          END DO
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the CKD coeffs in the data file in temp,freq
c this uses spline interpolations ........ blooody slow
      SUBROUTINE ComputeCKD_SlowSpline(d1,d2,df,daTemprt,daaCKD,
     $                      raFreq,raTemp,kFrStep,
     $                      daaCon,iDoDQ,daaDQ,daaDT,iMin,iMax)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c daaCon    = output continuum coeffs (normalised)
c daaDQ     = analytic Jacobian wrt gas amount
c daaDT     = analytic Jacobian wrt temperature
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
c raFreq    = wavenumber array
c kFrStep   = kFreqStep
c iMin,iMax = layers over which to loop
      DOUBLE PRECISION daaCon(kMaxPts,kProfLayer)
      INTEGER iDoDQ,iMin,iMax
      REAL raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
c these were from the CKD binary file
      DOUBLE PRECISION d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

c local variables
      INTEGER iFr,iL,iF,iFloor
      INTEGER iaFrIndex(kMaxPts)
      DOUBLE PRECISION daFrDelta(kMaxPts),dTemp,daC(kTempCKD)
      DOUBLE PRECISION dyp1,dypn,daY2(kTempCKD),daWork(kTempCKD)

C     Assign values for interpolation
C     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
      dYP1=1.0E+16
      dYPN=1.0E+16

      dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

c for each freq point in raFreq, find where nearest low CKD freq grid point is
      DO iFr=1,kMaxPts
        iaFrIndex(iFr)=1 + iFloor(real((raFreq(iFr)-d1)/df))
        daFrDelta(iFr)=raFreq(iFr)*1d0-(d1+(iaFrIndex(iFr)-1)*df)
      END DO   

      DO iFr=1,kMaxPts
        iF=iaFrIndex(iFr)   !index of closest freq lower than raFreq(iFr)
        DO iL=1,kTempCKD
          daC(iL)=daaCKD(iL,iF) + 
     $          daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
        END DO
        CALL DSPLY2(daTemprt,daC,kTempCKD,dYP1,dYPN,daY2,daWork) 
        DO iL=iMin,iMax
          CALL DSPLIN(daTemprt,daC,daY2,kTempCKD,dble(raTemp(iL)),daaCon(iFr,iL)) 
        END DO
      END DO

      RETURN
      END

c************************************************************************
c************************************************************************
c************************************************************************
c**  this is the old code (cross sections and continuum gateway calls) **
c************************************************************************
c************************************************************************
c************************************************************************
c this subroutine is the gateway call to XSEC/calxsc 
c which is now REDUNDANT (stored in CONTINUUM_BLOCKDATA_AND_OLDXSEC)
c compute the contribution of Gases 29-63 (if present) 
      SUBROUTINE CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep, 
     $      daaTemp,raVTemp,iVTSet,raFreq,iErr,caXsecName, 
     $      raTAmt,raTTemp,raTPress,raTPart,iaCont, 
     $      daaDQ,daaDT,iDoDQ) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
c iCount    = which of the iNumGases is being processed 
c iGasID    = iaGasID(iCount) = gas ID of current gas 
c iRefLayer = number of layers in the reference profiles (=kProfLayer) 
c iL,iU     = min/max layer number for each gas profile (=1,kProfLayer) 
c iaCont    = whether or not to do continuum calculation .. iaCont(iCount) 
c caXecF    = file name of cross section data 
c daaTemp   = matrix containing the uncompressed k-spectra 
c raVtemp   = vertical temperature profile for the Mixed paths 
c iVTSet    = has the vertical temp been set, to check current temp profile 
c raFreq    = wavenumber array 
c daaDQ     = analytic jacobian wrt gas amount 
c daaDT     = analytic jacobian wrt temperature 
c iErr      = errors (mainly associated with file I/O) 
c iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs 
c kFrStep  = kFreqStep 
      INTEGER iaCont(kMaxGas),iDoDQ 
      CHARACTER*80 caXsecName 
      DOUBLE PRECISION daaTemp(kMaxPts,kProfLayer) 
      INTEGER iGasID,iL,iU,iErr,iRefLayer,iCount,iVTSet 
      REAL raTPress(kProfLayer),raTPart(kProfLayer),kFrStep 
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer) 
      REAL raFreq(kMaxPts),raVTemp(kProfLayer)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac) 
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac) 
 
c local variables 
      REAL raXamnt(kProfLayer),raAdjustProf(kProfLayer),rMin 
      REAL raaXsec(kMaxPts,kProfLayer),rFStep,rCheckTemp 
      INTEGER iNmol, iaMolid(MXXMOL),iLay,iFr,iWhichUsed
 
c used in d/dq, d/dT if kSpline = 1 
      REAL raaSplineDQ(kMaxPtsJac,kProfLayerJac) 
      REAL raaSplineDT(kMaxPtsJac,kProfLayerJac) 
      INTEGER idQT 
 
      REAL AVOG 
      DATA AVOG/6.022045E+26/ 
 
c to see if all the abs coefficients are less than zero ...tau=exp(-abs) <=1 
      REAL rEC,rECCount 
 
c first check to see if exact calculations of d/dq,d/dT should be enabled 
      IF (kJacobian .GT. 0) THEN  
        idQT=1 
      ELSE 
        idQT=-1 
      END IF 

      iLay = kProfLayer 
      iFr = kMaxPts 
      rFStep = kFrStep 
 
c only calculate the current gas cross-section 
      DO iLay=1,MXXMOL 
        iaMolid(iLay)=0 
      END DO 
      iNmol=1 
      iaMolid(1)=iGasID 
 
      iLay = kProfLayer 
 
      IF (iErr .LE. 0) THEN 
c set the vertical temperature profile if iCount=1 (first profile read) 
        IF ((iVTSet .LT. 0) .AND. (iGasID .LE. kGasComp)) THEN  
          write(kStdWarn,*) 'Setting vertical temp profile ...' 
          DO iLay=1,kProfLayer 
            raVTemp(iLay)=raTTemp(iLay) 
          END DO 
        END IF 
c if previous profiles have been read in, check to make sure the  
c temperature profiles are the same!!!! 
        IF ((iVTSet .GT. 0) .AND. (iGasID .LE. kGasComp)) THEN  
          write(kStdWarn,*) 'Checking the vertical temp profile ...' 
          DO iLay=1,kProfLayer 
            rCheckTemp=raTTemp(iLay)-raVTemp(iLay) 
            IF (abs(rCheckTemp) .GE. 1.0e-3) THEN 
              write(kStdWarn,*) 'Warning!!Temp profiles do not match!!!' 
              write(kStdWarn,*) 'Gas#,layer, gastemp, vertical temp = ' 
              write(kStdWarn,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay) 
              write(kStdErr,*) 'Warning!!Temp profiles do not match!!!' 
              write(kStdErr,*) 'Gas#,layer, gastemp, vertical temp = ' 
              write(kStdErr,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay) 
            END IF 
          END DO 
        END IF 
      END IF 
 
      IF (iErr .LT. 0) THEN 
c convert amount from k.moles/cm^2 to molecules/cm^2 
        DO iLay=1,kProfLayer 
          raXamnt(iLay)=raTAmt(iLay)*AVOG 
          raAdjustProf(iLay)=raTTemp(iLay) 
        ENDDO 
 
c make sure that the relevant area in the cross section matrix is zeroed 
        DO iLay=1,kProfLayer 
          DO iFr=1,kMaxPts 
            raaXsec(iFr,iLay)=0.0 
          END DO 
        END DO 
 
c if we need to do jacobian calculations using splines, 
c initialize the matrices here  
        IF (kJacobian .GT. 0) THEN 
          IF (iDoDQ .GT. 0) THEN 
            DO iLay=1,kProfLayerJac 
              DO iFr=1,kMaxPtsJac 
                raaSplineDQ(iFr,iLay)=0.0 
                raaSplineDT(iFr,iLay)=0.0 
              END DO 
            END DO 
          ELSE IF (iDoDQ .LT. 0) THEN 
            DO iLay=1,kProfLayerJac 
              DO iFr=1,kMaxPtsJac 
                raaSplineDT(iFr,iLay)=0.0 
              END DO 
            END DO 
          END IF   
        END IF   
 
 
c this call calculates ONLY the cross sections  
c            OR 
c as idQT=1 this call calculates the cross sections, AND d/dq,d/dT!!!, saving 
c them in raaSplineDQ ,raaSplineDT respectively 

        iLay = kProfLayer 
        iFr = kMaxPts 
        rFStep = kFrStep 

        write(kStdErr,*) 'CALXSC no longer supported!!!!'
        CALL DoStop
cccccc        CALL CALXSC(caXsecName,iFr,raFreq,rFStep, 
cccccc     $      iLay, raAdjustProf,raXamnt,iNmol,iaMolid, 
cccccc     $      raaXsec,iWhichUsed, 
cccccc     $      idQT, raaSplineDQ ,raaSplineDT, iGasID, iDoDQ) 
 
c since CALXSC computed the cross-section and saved it in a 3d matrix, 
c pull out the relevant cross section data into a 2d matrix and add it to  
c the cumulative abscoeff matrix daaTemp 
c first check whether the CALXSC routine actually found a matching  
c gas+wavenumber region 
        rMin=1.0e10 
        rEC=0.0 
        rECCount=0.0 
        IF ((iWhichUsed .GT. 0) .AND. (iWhichUsed .LE. iNmol)) THEN 
          WRITE(kStdWarn,1000) iCount,iGasID 
 1000          FORMAT('adding on cross section for GAS(',I2,') = ',I2) 
          DO iLay=1,kProfLayer 
            DO iFr=1,kMaxPts 
              IF (raaXsec(iFr,iLay) .LT. 0.0) THEN 
c flag this error, as the cross section values should all be > 0 
                iErr=1 
                rEC=rEC+raaXsec(iFr,iLay) 
                rECCount=rECCount+1 
                IF (raaXsec(iFr,iLay) .LT. rMin) THEN 
                  rMin = raaXsec(iFr,iLay) 
                END IF 
                raaXsec(iFr,iLay)=0.0 
              END IF 
              daaTemp(iFr,iLay)=raaXsec(iFr,iLay) 
            END DO 
          END DO 
          ELSE 
            WRITE(kStdWarn,1010) iCount,iGasID 
 1010              FORMAT ('no XSEC contribution due to GAS(',I2,') = ',I2) 
          END IF 
        IF (rECCount .gt. 0.5) THEN 
          rEC=rEC/rECCount 
          write(kStdWarn,*) 'Error in XSEC data!!! Some values negative!' 
          write(kStdWarn,*) 'and reset to 0.0' 
          write(kStdWarn,*) 'rECCount values have avg value ',rEC 
          write(kStdWarn,*) 'min negative value in 10000*100 = ',rMin 
          IF (abs(rMin) .GT. 1.0e-7) THEN 
            iErr=1 
            CALL DoSTOP 
            END IF 
          END IF           
c end main if statement 
        END IF 
 
c do the inclusion of the exact derivatives here 
      IF (kJacobian .GT. 0) THEN 
c exact calculations have already been performed in calXSC,calq ...  
c so apart from the avog factor in raaDQ, no need to do much here!!! 
        IF (iDoDQ .GT. 0) THEN           
          DO iLay=1,kProfLayerJac 
            DO iFr=1,kMaxPtsJac 
              daaDQ(iFr,iLay)=raaSplineDQ(iFr,iLay)*avog 
              daaDT(iFr,iLay)=raaSplineDT(iFr,iLay) 
            END DO 
          END DO 
        ELSE IF (iDoDQ .LT. 0) THEN 
          DO iLay=1,kProfLayerJac 
            DO iFr=1,kMaxPtsJac 
              daaDT(iFr,iLay)=raaSplineDT(iFr,iLay) 
            END DO 
          END DO 
        END IF 
      END IF 
         
      RETURN 
      END 
 
c************************************************************************ 
c this gets the temperature dependant multipliers for CKD1
      SUBROUTINE GetCKD1Mult(iGasID,raTemp,raFreq,
     $                       iNpts,dxAERI,psT,psK,pfT,pfK,
     $                       daaCKD1Mult)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params
      INTEGER iGasID
      REAL raTemp(kProfLayer),raFreq(kMaxPts)
c output params 
c these are read from the data file
      INTEGER iNpts
      DOUBLE PRECISION psT(kMaxLayer),psK(kMaxLayer)
      DOUBLE PRECISION pfT(kMaxLayer),pfK(kMaxLayer)
      DOUBLE PRECISION dxAERI(kMaxLayer)
c this are the temperature corrected coeffs
      DOUBLE PRECISION daaCKD1Mult(kMaxPts,kProfLayer)

c local vars
      CHARACTER*80 FNAM
      INTEGER iL,iFr,iN,iIOUN,iErr
      DOUBLE PRECISION daXT(kMaxLayer),daXK(kMaxLayer)
      DOUBLE PRECISION dT,daX(kMaxLayer),daM(kMaxPts),x,y,daFreq(kMaxPts)

      FNAM = '/home/sergio/KCARTA/INCLUDE/ckd1mult.dat'

      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
          WRITE(kStdErr,1010) IERR, FNAM
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
        ENDIF
      kTempUnitOpen=1

      iNpts = 0
      IF (iGasID .EQ. kNewGasLo) THEN
 10     CONTINUE
        iNpts = iNpts + 1
        READ(iIOUN,*,ERR=20) dxAERI(iNpts),daXT(iNpts),daXK(iNpts),x,y
c        print *,iNpts,dxAERI(iNpts),daXT(iNpts),daXK(iNpts)
        GOTO 10
      ELSEIF (iGasID .EQ. kNewGasHi) THEN
 11     CONTINUE
        iNpts = iNpts + 1
        READ(iIOUN,*,ERR=20) dxAERI(iNpts),x,y,daXT(iNpts),daXK(iNpts)
c        print *,iNpts,dxAERI(iNpts),daXT(iNpts),daXK(iNpts)
        GOTO 11
      END IF

 20   CONTINUE
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      DO iFr = 1,kMaxPts
        daFreq(iFr) = raFreq(iFr) * 1.0d0
      END DO

      DO iL = 1,kProfLayer
        dT = raTemp(iL)*1.0d0
        IF (dT .GE. 296.0d0) THEN  !!! no need to adjust things
          DO iFr = 1,kMaxPts
            daaCKD1Mult(iFr,iL) = 1.0d0
          END DO
        ELSE
          DO iFr = 1,iNpts
            !!!! find the temperature dependance
            daX(iFr) = daXT(iFr)*dT + daXK(iFr)
          END DO
          CALL dspl(dxAERI,daX,iNpts,daFreq,daM,kMaxPts)
          DO iFr = 1,kMaxPts
            daaCKD1Mult(iFr,iL) = daM(iFr)
          END DO
        END IF
c        print *,iL,iGasID,iNpts,kProfLayer,dT,raFreq(1),daaCKD1Mult(1,iL)
      END DO
                    
      RETURN 
      END 
 
c************************************************************************ 
