! Copyright 2011
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcont_xsec

USE basic_common
USE spline_and_sort_and_common
USE kcoeff_FAST
USE s_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
!********* this file has the continuum/old xsec routines ****************
!************************************************************************
! this adds on the continuum :
    SUBROUTINE AddContinuum(iGasID,iTag,iActualTag,rFileStartFr, &
    iRefLayer,iProfileLayers,raFreq,raAmt,raTemp, &
    kFrStep,raPress,raPartPress, &
    iL,iU,daaTemp,daaDQ,daaDT,iDoDQ,iSPlineType, &
    raRAmt,raRTemp,raRPress,raRPartPress,pProf, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
! iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
! daaTemp   = matrix containing the uncompressed k-spectra
! raFreq    = wavenumber array
! iErr      = errors (mainly associated with file I/O)
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! kFrStep  = kFreqStep
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    INTEGER :: iGasID,iL,iU,iErr,iRefLayer,iDoDQ
    INTEGER :: iTag,iActualTag,iSplineType,iProfileLayers
    REAL :: raFreq(kMaxPts),kFrStep,pProf(kProfLayer),rFileStartFr
    REAL :: raAmt(kProfLayer),raTemp(kProfLayer)
    REAL :: raPartPress(kProfLayer),raPress(kProfLayer)
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! the Matlab weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

! iMethod   = 0 for old style, 1 for new style
! local variables
    INTEGER :: iLay,iFr,iIOUN,iLoc
    REAL :: rFStep,rMin,rConvFac
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    INTEGER :: iMin,iMax,iErr1

    INTEGER :: iDefault,iUseMethod,iMethod

    iDefault = +1      !!! new method, uses kCompressed chunks
    iDefault = +0      !!! original method, uses 2901 (fr x T) lookup table

    iDefault = +0      !!! original method, uses 2901 (fr x T) lookup table
    iMethod = 0        !!! assume we need OLD method

    IF ((abs(iSplineType) == 2) .AND. (iActualTag == 20)) THEN
      !! warning so far only CKD1,6,25 has these tables ....
      iMethod = 1        !! use compressed data for 605-2830 cm-1
      iDefault = +1      !! new method, uses kCompressed chunks

      !decided to switch back to look up tables, as Matlab version uses this
      iMethod = 0        !! use look up tables
      iDefault = 0       !! use look up tables
    END IF

    IF (iDefault /= iMethod) THEN
      write(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
    END IF

    CALL ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)

    IF (iMethod == 0) THEN
      ! reads the oldstyle 2901 point lookup table which is fcn(fr,T)
      ! old style, been using this for years
      write (kStdWarn,*) ' WV cont : using Lookup Table (fr x T) format'
      rConvFac = 1.0   !!! scale factor = 1 in the ORIG tables
      CALL f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ, &
        rFStep,raFreq,kFrStep, &
        raTemp,daaDQ,daaDT,daaCon)
    ELSEIF (iMethod == 1) THEN
      ! reads the Compressed Style lookup fcn(fr,T)
      ! new style, since Sept 2011
      write (kStdWarn,*) ' WV cont : using kCompressed Database format'
      rConvFac = 1.0e-23    !! SELF/FORN continuum coeffs in compressed
                              !! database were remultiplied by 1.0e+23
      CALL xCKDgases(iGasID,rFileStartFr,iTag,iActualTag, &
        iRefLayer,iL,iU, &
        raAmt,raRAmt,raTemp,raRTemp, &
        iErr1,iDoDQ,pProf,iProfileLayers, &
        daaDQ,daaDT,daaCon,iSplineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
        iaQ21,iaQ22,raQ21,raQ22)
    END IF

    CALL aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt, &
                raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)

    RETURN
    end SUBROUTINE AddContinuum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! does some simple inits
    SUBROUTINE ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
    REAL :: raAmt(kProfLayer),kFrStep
! output
    INTEGER :: iMin,iMax,iLay,iFr
    REAL :: rFStep
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
          
! local
    INTEGER :: iN

! everywhere we have         DO iLay=1,kProfLayer
!  we replace with           DO iLay=iMin,iMax
!      iMin = 1
!      iMax = kProfLayer
    iMax = -1
    DO iN = 1,kProfLayer
      IF ((raAmt(iN) > 0) .AND. (iN > iMax)) THEN
        iMax = iN
      END IF
    END DO
    iMin = kProfLayer + 1
    DO iN = kProfLayer,1,-1
      IF ((raAmt(iN) > 0) .AND. (iN < iMin)) THEN
        iMin = iN
      END IF
    END DO

! initialize raaCon to all zeros
    daaCon = 0.0

    iLay = kProfLayer
    iFr = kMaxPts
    rFStep = kFrStep

    RETURN
    end SUBROUTINE ckd_init

!************************************************************************
    SUBROUTINE f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ, &
    rFStep,raFreq,kFrStep, &
    raTemp,daaDQ,daaDT,daaCon)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
    INTEGER :: iGasID,iTag,iDoDQ
    REAL :: raFreq(kMaxPts),rFStep,kFrStep
    REAL :: raTemp(kProfLayer)
! output
!   daaDQ     = analytic Jacobian wrt gas amount
!   daaDT     = analytic Jacobian wrt temperature
!   daaCon    = continuum coeffs (multiplied by 1)
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

! local vars
    INTEGER :: iLay,iFr,iIOUN,iLoc,iDefault,iUseMethod,iMethod,iErr
    REAL :: rMin
    DOUBLE PRECISION :: der1,der2,doop
    REAL :: rEC,rECCount

! these are to read in the binary file
    DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)
    INTEGER :: iCKD,iM,iN

    INTEGER :: iMin,iMax
    REAL :: a1,a2
    CHARACTER(120) :: caFName

! figure out the filename
    CALL CKDFileName(caFName,iGasID,iTag)

 1080 FORMAT('ERROR! number ',I5,' opening CKD binary file:',/,A120)          
    iIOUN = kTempUnit
! open file and load in data
    OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='UNFORMATTED', &
      IOSTAT=iErr)
    IF (iErr /= 0) THEN
      WRITE(kStdErr,1080) iErr, caFname
      WRITE(kStdWarn,1080) iErr, caFname
      CALL DoSTOP
    ENDIF

    kTempUnitOpen=1
    READ(iIOUN) d1,d2,df       !read start/stop freq, df
    IF (d1 > raFreq(1)) THEN
      write(kStdErr,*) 'CKD file has freqs that start too high!'
      CALL DoStop
    END IF
    IF (d2 < raFreq(kMaxPts)) THEN
      write(kStdErr,*) 'CKD file has freqs end start too low!'
      CALL DoStop
    END IF

!&&&&&&&&
    READ(iIOUN) iLoc,iCKD,iM,iN
! ead line type, CKD vers, # of temp,freq pts
    IF (iLoc < 0) THEN
      write(kStdErr,*) iLoc
      write(kStdErr,*) 'need correct continuum for local lineshape !'
      CALL DoStop
    END IF
    IF (iCKD /= kCKD) THEN
      write(kStdErr,*) iCKD,kCKD
      write(kStdErr,*) 'CKD versions do not match!'
      CALL DoStop
    END IF
    IF (iM /= kTempCKD) THEN
      write(kStdErr,*) iM,kTempCKD
      write(kStdErr,*) 'Need more Temp offsets in CKD binary file'
      CALL DoStop
    END IF
    IF (iN /= kFreqCKD) THEN
      write(kStdErr,*) iN,kFreqCKD
      write(kStdErr,*) 'Read in ',iN, ' points from CKD binary file'
      write(kStdErr,*) 'Always need ',kFreqCKD, ' points in the file'
      CALL DoStop
    END IF

    READ(iIOUN) (daTemprt(iLay),iLay=1,iM)   !read the temps

    DO iLay=1,iM
      ! or current temp, read the continuum data for each freq
      READ(iIOUN) (daaCKD(iLay,iFr),iFr=1,iN)
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen=-1
!&&&&&&&&
        
    iMethod = +1       !!! quadratic
    iMethod = +2       !!! linear, T dependance in Cs,Cf
    iMethod = +3       !!! linear, T dependance in Cs only

    iDefault = +1      !!! used by Scott in making Fast Models
    iDefault = +3      !!! been used in KCARTA since 2005

    IF (iDefault /= iMethod) THEN
      write(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
    END IF

! interpolate
    IF (iMethod == 1) THEN
      ! this is what was used in SRCv1.07, which is what Scott Hannon used
      ! for the Fast Models
      CALL ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp, &
        kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

    ELSEIF (iMethod == 2) THEN
      ! this is new and improved version, which should be quicker than Quad
      ! THIS ONE has T dependnce both for self and foreign
      CALL ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp, &
        kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

    ELSEIF (iMethod == 3) THEN
      ! this is new and improved version, which should be quicker than Quad
      ! THIS ONE has T dependance ONLY for self
      CALL ComputeCKD_Linear_March2002(d1,d2,df, &
        daTemprt,daaCKD,raFreq,raTemp, &
        kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)
    END IF

    RETURN
    end SUBROUTINE f_by_T_CKDfilereader

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! this subroutine completes the CKD calcs, and d/dT,d/dQ if needed
    SUBROUTINE aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt, &
    raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iMin,Imax tell over what layers to do the CKD calcs (where gasamt nonzero!)
! input
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)       !! ckd coeffs
    INTEGER :: iMin,iMax,iGasID,iDoDQ
    REAL :: raAmt(kProfLayer),raTemp(kProfLayer)
    REAL :: raPartPress(kProfLayer),raPress(kProfLayer)
    REAL :: raFreq(kMaxPts),kFrStep,rConvFac
! output
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)      !! ckd ODs
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iLay,iFr,iErr
    REAL :: rFStep,rMin
    DOUBLE PRECISION :: dader1(kMaxPts),dader2(kMaxPts),dadoop(kMaxPts)
    DOUBLE PRECISION :: daaJunk(kMaxPts,kProfLayer),dDiff
    REAL :: rEC,rECCount
    REAL :: raMult(kProfLayer)
    REAL :: a1,a2
    REAL :: raa2DAmt(kMaxPts,kProfLayer),raOne(kMaxPts)

    raOne = 1.0
    DO iLay = 1,kProfLayer
      raa2DAmt(:,kProfLayer) = raOne * raAmt(iLay)
    END DO

! raAmt(iLay) is in kmoles/cm2   pressures are in atm, CKD coeffs are in (molecules/cm2 1/cm-1)-1 1/atm
!      DO iLay =  iMin,iMax
!        print *,iLay,raAmt(iLay),raPartPress(iLay),raPress(iLay)
!      END DO
!      CALL DoStop
! so raAmt*kAvog*pressure is in molecules/cm2*atm
! so CX * v tanh(hcv/2KT) * raAmt*kAvog*pressure = [(molecules/cm2 1/cm-1)-1 1/atm] [cm-1] [molecules/cm2*atm] = OD no units

    rMin = 1.0e10
    rECCount = 0.0
    daaJunk = daaCon
    WHERE (daaCon < 0)
      daaCon = 0.0d0
    END WHERE

! add the continuum abs coeffs to the k-compressed coeffs raaCon = daaTemp
! keeping track of whether the added values are all greater than 0.0
    IF (iGasID == kNewGasLo) THEN      !CSelf
      raMult(iMin:iMax)=raAmt(iMin:iMax)*kAvog*raPartPress(iMin:iMax)
    ELSEIF (iGasID == kNewGasHi) THEN      !CFor
      raMult(iMin:iMax)=raAmt(iMin:iMax)*kAvog*(raPress(iMin:iMax)-raPartPress(iMin:iMax))
    END IF

    DO iLay =iMin,iMax
      a1 = raMult(iLay)*296.0/raTemp(iLay)
      a2 = kPlanck2/2/raTemp(iLay)
      daaTemp(:,iLay) = rConvFac * daaCon(:,iLay)*a1*raFreq*tanh(a2*raFreq)
    END DO
    
    IF (iDoDQ >= -1) THEN
      IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
        ! the gas amount jacobians
        IF (iGasID == 101) THEN   !!!there is a factor of 2
          daaDQ(:,iMin:iMax) = 2 * daaTemp(:,iMin:iMax)/raa2DAmt(:,iMin:iMax)
        ELSEIF (iGasID == 102) THEN
          daaDQ(:,iMin:iMax) = daaTemp(:,iMin:iMax)/raa2DAmt(:,iMin:iMax)
        END IF
      END IF
         
      IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
        (kActualJacs == 32) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
        ! this is the temperature jacobian
        DO iLay=iMin,iMax
          dadoop = tanh(kPlanck2*raFreq/2/raTemp(iLay))
          dader1 = -(296/raTemp(iLay)**2)*dadoop
          dader1 = dader1-(296*kPlanck2*raFreq/2/(raTemp(iLay)**3))/ &
                   (cosh(kPlanck2*raFreq/2/raTemp(iLay))**2)
          dader1 = dader1*daaCon(:,iLay)

          dader2 = (296/raTemp(iLay))* dadoop
          dader2 = dader2*daaDT(:,iLay)
          daaDT(:,iLay) = raMult(iLay)*raFreq*(dader1+dader2)*rConvFac
        END DO
      END IF

    END IF

    RETURN
    end SUBROUTINE aux_ckd

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, linear in temperature
! does temperature interpolation for SELF and for foreign
! as this is the new CKD UMBC thingies
    SUBROUTINE ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD, &
    raFreq,raTemp,kFrStep, &
    daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    INTEGER :: iDoDQ,iGasID,iMin,iMax
    REAL :: raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! these were from the CKD binary file
    DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

! local variables
    INTEGER :: iLay,iFr,iL,iF
    INTEGER :: iaFrIndex(kMaxPts),iaF(kMaxPts),iaTempIndex(kProfLayer)
    DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp
    DOUBLE PRECISION :: daA(kMaxPts),daB(kMaxPts)
    DOUBLE PRECISION :: c,x1,x2,x3,y1,day2(kMaxPts),day3(kMaxPts),t1,t2,t3,t4,x,z1,z2
! this is for CKD1 temeperature dependant multipliers from 780-980 cm-1
    INTEGER :: iNpts,iMult
    DOUBLE PRECISION :: psT(kMaxLayer),psK(kMaxLayer)
    DOUBLE PRECISION :: pfT(kMaxLayer),pfK(kMaxLayer)
    DOUBLE PRECISION :: xNum(kMaxLayer)
    DOUBLE PRECISION :: daaCKD1Mult(kMaxPts,kProfLayer)

    dTemp = daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

    iMult = -1

! for each freq point in raFreq, find where nearest low CKD freq grid point is
    iaFrIndex = 1 + int(floor(real((raFreq-d1)/df)))
    iaF = iaFrIndex
    daFrDelta = raFreq*1.0d0-(d1+(iaFrIndex-1)*df)

! for each layer temp, find where the nearest low CKD tempr grid point is
    DO iLay = iMin,iMax
      iaTempIndex(iLay) = 1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
    END DO

    IF (iDoDQ < -1) THEN         !no need to do temp jacobians
      IF (iGasID == kNewGasLo) THEN
        ! self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find out line that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!

        END DO

      ELSEIF (iGasID == kNewGasHi) THEN
        ! foreign continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
          !find out line that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!
                    
        END DO
      END IF

   ELSE   !!!!!!!!!do the temp jacobians
      IF (iGasID == kNewGasLo) THEN
        !self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!
          daaDT(:,iLay)  = daA             !this is temp jacobian!

        END DO

      ELSEIF (iGasID == kNewGasHi) THEN
        !foreign continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
          !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!
          daaDT(:,iLay)  = daA             !this is temp jacobian!
                    
        END DO
      END IF
    END IF

    IF (iMult > 0) THEN
      daaCon(:,iMin:iMax) = daaCon(:,iMin:iMax) * daaCKD1Mult(:,iMin:iMax)
    END IF

    RETURN
    end SUBROUTINE ComputeCKD_Linear

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, linear in temperature
! only does temperature interpolation for SELF, and not for foreign
! as this was the CKD 0,2.1,2.3,2.4 thingies
    SUBROUTINE ComputeCKD_Linear_March2002(d1,d2,df,daTemprt,daaCKD, &
    raFreq,raTemp,kFrStep, &
    daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep  = kFreqStep
! iMin,iMax = layers over which to loop
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    INTEGER :: iDoDQ,iGasID,iMin,iMax
    REAL :: raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! these were from the CKD binary file
    DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

! local variables
    INTEGER :: iLay,iFr,iL,iF
    INTEGER :: iaFrIndex(kMaxPts),iaF(kMaxPts),iaTempIndex(kProfLayer)
    DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp
    DOUBLE PRECISION :: daA(kMaxPts),daB(kMaxPts)
    DOUBLE PRECISION :: c,x1,x2,x3,y1,day2(kMaxPts),day3(kMaxPts),t1,t2,t3,t4,x,z1,z2

    dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
    iaFrIndex = 1 + int(floor(real((raFreq-d1)/df)))
    iaF = iaFrIndex
    daFrDelta = raFreq*1.0d0-(d1+(iaFrIndex-1)*df)

! for each layer temp, find where the nearest low CKD tempr grid point is
    DO iLay = iMin,iMax
      iaTempIndex(iLay) = 1 + iFloor(real((raTemp(iLay)-daTemprt(1))/dTemp))
    END DO

    IF (iDoDQ < -1) THEN         !no need to do temp jacobians
      IF (iGasID == kNewGasLo) THEN
        !self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find out line that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!

         END DO

      ELSEIF (iGasID == kNewGasHi) THEN
        !foreign continuum has no temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
                     
          day2 = daaCKD(iL,iaF) + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df                              
          daaCon(:,iLay) = day2    !this is temp dependance!
        END DO
      END IF

    ELSE          !!!!!!!!!do the temp jacobians
      IF (iGasID == kNewGasLo) THEN
        !self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find quadratic that goes thru the 2 "x(j)" points to give "y(j)"
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          day2 = daaCKD(iL,iaF)   + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          day3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df                              
          daA = (day3-day2)/(x3-x2)
          daB = day3-daA*x3
          daaCon(:,iLay) = daA*x + daB     !this is temp dependance!
          daaDT(:,iLay)  = daA             !this is temp jacobian!

        END DO

      ELSEIF (iGasID == kNewGasHi) THEN
        !foreign continuum has no temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
                     
          day2 = daaCKD(iL,iaF) + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df                              
          daaCon(:,iLay) = day2    !this is temp dependance!
          daaDT(:,iLay) = 0.0d0  !this is temp jacobian! 
        END DO
      END IF
    END IF

    RETURN
    end SUBROUTINE ComputeCKD_Linear_March2002

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, quadratic in temperature
    SUBROUTINE ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD, &
    raFreq,raTemp,kFrStep, &
    daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    INTEGER :: iDoDQ,iGasID,iMin,iMax
    REAL :: raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! these were from the CKD binary file
    DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

! local variables
    INTEGER :: iLay,iFr,iL,iF
    INTEGER :: iaFrIndex(kMaxPts),iaTempIndex(kProfLayer),iaF(kMaxPts)
    DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp,daDeltaF(kMaxPts)
    DOUBLE PRECISION :: x1,x2,x3,y1(kMaxPts),y2(kMaxPts),y3(kMaxPts)
    DOUBLE PRECISION :: a(kMaxPts),b(kMaxPts),c(kMaxPts),t1,t2,t3,t4
    DOUBLE PRECISION :: x,z1(kMaxPts),z2(kMaxPts)
    REAL :: raDeltaF(kMaxPts)

    dTemp = daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
    daDeltaF = (raFreq*1.0d0-d1)/df
    raDeltaF = real(daDeltaF)
    iaFrIndex = 1 + int(floor(raDeltaF))
    iaF = iaFrIndex
    daFrDelta = raFreq*1.0d0-(d1+(iaFrIndex-1)*df)

! for each layer temp, find where the nearest low CKD tempr grid point is
    iaTempIndex(iMin:iMax) = 1 + int(floor(real((raTemp(iMin:iMax)-daTemprt(1))/dTemp)))

    IF (iDoDQ < -1) THEN         !no need to do temp jacobians
      IF (iGasID == kNewGasLo) THEN
        ! self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
          x1 = daTemprt(iL-1)
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          y1 = daaCKD(iL-1,iaF) + daFrDelta*(daaCKD(iL-1,iaF+1)-daaCKD(iL-1,iaF))/df
          y2 = daaCKD(iL,iaF) +   daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          y3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df
                            
          z1 = y1-y3
          z2 = y2-y3
          t1 = x1*x1-x3*x3
          t2 = x1-x3
          t3 = x2*x2-x3*x3
          t4 = x2-x3
          b = (z1*t3-z2*t1)/(t3*t2-t1*t4)
          a = (z2-b*t4)/t3
          c = y3-a*x3*x3-b*x3

          daaCon(:,iLay) = a*x*x + b*x + c     !this is temp dependance!
        END DO

      ELSEIF (iGasID == kNewGasHi) THEN
        ! foreign continuum has no temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)

          y2 = daaCKD(iL,iaF) + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
                              
          daaCon(:,iLay) = y2    !this is temp dependance!
                    
        END DO
      END IF
    ELSE          !!!!!!!!!do the temp jacobians
      IF (iGasID == kNewGasLo) THEN
        !self continuum has temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          iL = iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
          !find quadratic that goes thru the 3 "x(j)" points to give "y(j)"
          x1 = daTemprt(iL-1)
          x2 = daTemprt(iL)
          x3 = daTemprt(iL+1)
          x = raTemp(iLay)

          y1 = daaCKD(iL-1,iaF) + daFrDelta*(daaCKD(iL-1,iaF+1)-daaCKD(iL-1,iaF))/df
          y2 = daaCKD(iL,iaF) +   daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
          y3 = daaCKD(iL+1,iaF) + daFrDelta*(daaCKD(iL+1,iaF+1)-daaCKD(iL+1,iaF))/df
                            
          z1 = y1-y3
          z2 = y2-y3
          t1 = x1*x1-x3*x3
          t2 = x1-x3
          t3 = x2*x2-x3*x3
          t4 = x2-x3
          b = (z1*t3-z2*t1)/(t3*t2-t1*t4)
          a = (z2-b*t4)/t3
          c = y3-a*x3*x3-b*x3

          daaCon(:,iLay) = a*x*x + b*x + c     !this is temp dependance!
          daaDT(:,iLay)  = 2*a*x + b           !this is temp jacobian!

        END DO
      ELSEIF (iGasID == kNewGasHi) THEN
        !f oreign continuum has no temp dependance
        DO iLay = iMin,iMax
          iL = iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)

          y2 = daaCKD(iL,iaF) + daFrDelta*(daaCKD(iL,iaF+1)-daaCKD(iL,iaF))/df
                              
          daaCon(:,iLay) = y2    !this is temp dependance!
          daaDT(:,iLay) = 0.0d0  !this is temp jacobian!
                    
        END DO
       END IF
    END IF

    RETURN
    end SUBROUTINE ComputeCKD_Quadratic

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses spline interpolations ........ blooody slow
    SUBROUTINE ComputeCKD_SlowSpline(d1,d2,df,daTemprt,daaCKD, &
    raFreq,raTemp,kFrStep, &
    daaCon,iDoDQ,daaDQ,daaDT,iMin,iMax)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop
    DOUBLE PRECISION :: daaCon(kMaxPts,kProfLayer)
    INTEGER :: iDoDQ,iMin,iMax
    REAL :: raFreq(kMaxPts),raTemp(kProfLayer),kFrStep
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
! these were from the CKD binary file
    DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)

! local variables
    INTEGER :: iFr,iL,iF
    INTEGER :: iaFrIndex(kMaxPts)
    DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp,daC(kTempCKD)
    DOUBLE PRECISION :: dyp1,dypn,daY2(kTempCKD),daWork(kTempCKD)

!     Assign values for interpolation
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1 = 1.0E+16
    dYPN = 1.0E+16

    dTemp = daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
    iaFrIndex = 1 + int(floor(real((raFreq-d1)/df)))
    daFrDelta = raFreq*1.0d0-(d1+(iaFrIndex-1)*df)

    DO iFr=1,kMaxPts
      iF=iaFrIndex(iFr)   !index of closest freq lower than raFreq(iFr)
      DO iL=1,kTempCKD
        daC(iL)=daaCKD(iL,iF) + daFrDelta(iFr)*(daaCKD(iL,iF+1)-daaCKD(iL,iF))/df
      END DO
      CALL DSPLY2(daTemprt,daC,kTempCKD,dYP1,dYPN,daY2,daWork)
      DO iL=iMin,iMax
        CALL DSPLIN(daTemprt,daC,daY2,kTempCKD,dble(raTemp(iL)),daaCon(iFr,iL))
      END DO
    END DO

    RETURN
    end SUBROUTINE ComputeCKD_SlowSpline

!************************************************************************
!************************************************************************
!************************************************************************
!**  this is the old code (cross sections and continuum gateway calls) **
!************************************************************************
!************************************************************************
!************************************************************************
! this subroutine is the gateway call to XSEC/calxsc
! which is now REDUNDANT (stored in CONTINUUM_BLOCKDATA_AND_OLDXSEC)
! compute the contribution of Gases 29-63 (if present)
    SUBROUTINE CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep, &
    daaTemp,raVTemp,iVTSet,raFreq,iErr,caXsecName, &
    raTAmt,raTTemp,raTPress,raTPart,iaCont, &
    daaDQ,daaDT,iDoDQ)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! iCount    = which of the iNumGases is being processed
! iGasID    = iaGasID(iCount) = gas ID of current gas
! iRefLayer = number of layers in the reference profiles (=kProfLayer)
! iL,iU     = min/max layer number for each gas profile (=1,kProfLayer)
! iaCont    = whether or not to do continuum calculation .. iaCont(iCount)
! caXecF    = file name of cross section data
! daaTemp   = matrix containing the uncompressed k-spectra
! raVtemp   = vertical temperature profile for the Mixed paths
! iVTSet    = has the vertical temp been set, to check current temp profile
! raFreq    = wavenumber array
! daaDQ     = analytic jacobian wrt gas amount
! daaDT     = analytic jacobian wrt temperature
! iErr      = errors (mainly associated with file I/O)
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! kFrStep  = kFreqStep
    INTEGER :: iaCont(kMaxGas),iDoDQ
    CHARACTER(80) :: caXsecName
    DOUBLE PRECISION :: daaTemp(kMaxPts,kProfLayer)
    INTEGER :: iGasID,iL,iU,iErr,iRefLayer,iCount,iVTSet
    REAL :: raTPress(kProfLayer),raTPart(kProfLayer),kFrStep
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
     
! local variables
    REAL :: raXamnt(kProfLayer),raAdjustProf(kProfLayer),rMin
    REAL :: raaXsec(kMaxPts,kProfLayer),rFStep,rCheckTemp
    INTEGER :: iNmol, iaMolid(MXXMOL),iLay,iFr,iWhichUsed
     
! used in d/dq, d/dT if kSpline = 1
    REAL :: raaSplineDQ(kMaxPtsJac,kProfLayerJac)
    REAL :: raaSplineDT(kMaxPtsJac,kProfLayerJac)
    INTEGER :: idQT
     
    REAL :: AVOG
    DATA AVOG/6.022045E+26/
     
! to see if all the abs coefficients are less than zero ...tau=exp(-abs) <=1
    REAL :: rEC,rECCount
     
! first check to see if exact calculations of d/dq,d/dT should be enabled
    IF (kJacobian > 0) THEN
        idQT=1
    ELSE
        idQT=-1
    END IF

    iLay = kProfLayer
    iFr = kMaxPts
    rFStep = kFrStep
     
! only calculate the current gas cross-section
    iaMolid = 0
    iNmol=1
    iaMolid(1)=iGasID
     
    iLay = kProfLayer
     
    IF (iErr <= 0) THEN
      ! set the vertical temperature profile if iCount=1 (first profile read)
      IF ((iVTSet < 0) .AND. (iGasID <= kGasComp)) THEN
        write(kStdWarn,*) 'Setting vertical temp profile ...'
        raVTemp = raTTemp
      END IF
      ! if previous profiles have been read in, check to make sure the
      ! temperature profiles are the same!!!!
      IF ((iVTSet > 0) .AND. (iGasID <= kGasComp)) THEN
        write(kStdWarn,*) 'Checking the vertical temp profile ...'
        DO iLay=1,kProfLayer
          rCheckTemp=raTTemp(iLay)-raVTemp(iLay)
          IF (abs(rCheckTemp) >= 1.0e-3) THEN
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
     
    IF (iErr < 0) THEN
      ! convert amount from k.moles/cm^2 to molecules/cm^2
      raXamnt = raTAmt*AVOG
      raAdjustProf = raTTemp
         
      ! make sure that the relevant area in the cross section matrix is zeroed
      raaXsec = 0.0
         
      ! if we need to do jacobian calculations using splines,
      ! initialize the matrices here
      IF (kJacobian > 0) THEN
        IF (iDoDQ > 0) THEN
          raaSplineDQ = 0.0
          raaSplineDT = 0.0
        ELSE IF (iDoDQ < 0) THEN
           raaSplineDT = 0.0
        END IF
      END IF
         
         
      ! this call calculates ONLY the cross sections
      !            OR
      ! as idQT=1 this call calculates the cross sections, AND d/dq,d/dT!!!, saving
      ! them in raaSplineDQ ,raaSplineDT respectively

      iLay = kProfLayer
      iFr = kMaxPts
      rFStep = kFrStep

      write(kStdErr,*) 'CALXSC no longer supported!!!!'
      CALL DoStop
      ! cccc        CALL CALXSC(caXsecName,iFr,raFreq,rFStep,
      ! cccc     $      iLay, raAdjustProf,raXamnt,iNmol,iaMolid,
      ! cccc     $      raaXsec,iWhichUsed,
      ! cccc     $      idQT, raaSplineDQ ,raaSplineDT, iGasID, iDoDQ)
         
      ! since CALXSC computed the cross-section and saved it in a 3d matrix,
      ! pull out the relevant cross section data into a 2d matrix and add it to
      ! the cumulative abscoeff matrix daaTemp
      ! first check whether the CALXSC routine actually found a matching
      ! gas+wavenumber region
      rMin=1.0e10
      rEC=0.0
      rECCount=0.0
      IF ((iWhichUsed > 0) .AND. (iWhichUsed <= iNmol)) THEN
        WRITE(kStdWarn,1000) iCount,iGasID
 1000   FORMAT('adding on cross section for GAS(',I2,') = ',I2)
        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            IF (raaXsec(iFr,iLay) < 0.0) THEN
              ! flag this error, as the cross section values should all be > 0
              iErr=1
              rEC=rEC+raaXsec(iFr,iLay)
              rECCount=rECCount+1
              IF (raaXsec(iFr,iLay) < rMin) THEN
                  rMin = raaXsec(iFr,iLay)
              END IF
              raaXsec(iFr,iLay)=0.0
            END IF
            daaTemp(iFr,iLay)=raaXsec(iFr,iLay)
          END DO
        END DO
      ELSE
        WRITE(kStdWarn,1010) iCount,iGasID
 1010   FORMAT ('no XSEC contribution due to GAS(',I2,') = ',I2)
      END IF
      IF (rECCount > 0.5) THEN
        rEC=rEC/rECCount
        write(kStdWarn,*) 'Error in XSEC data!!! Some values negative!'
        write(kStdWarn,*) 'and reset to 0.0'
        write(kStdWarn,*) 'rECCount values have avg value ',rEC
        write(kStdWarn,*) 'min negative value in 10000*100 = ',rMin
        IF (abs(rMin) > 1.0e-7) THEN
          iErr=1
          CALL DoSTOP
        END IF
      END IF
      ! end main if statement
    END IF
     
! do the inclusion of the exact derivatives here
    IF (kJacobian > 0) THEN
      ! exact calculations have already been performed in calXSC,calq ...
      ! so apart from the avog factor in raaDQ, no need to do much here!!!
      IF (iDoDQ > 0) THEN
        daaDQ = raaSplineDQ * avog
        daaDT = raaSplineDT
      ELSE IF (iDoDQ < 0) THEN
        daaDT = raaSplineDT
      END IF
    END IF
             
    RETURN
    end SUBROUTINE CrossSectionOLD
     

!************************************************************************
END MODULE kcont_xsec
