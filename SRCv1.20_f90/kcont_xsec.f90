! Copyright 2011
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:41
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!********* this file has the continuum/old xsec routines ****************
!************************************************************************
! this adds on the continuum :

SUBROUTINE AddContinuum(iGasID,iTag,iActualTag,rFileStartFr,  &
    iRefLayer,iProfileLayers,raFreq,raAmt,raTemp, kFrStep,raPress,raPartPress,  &
    iL,iU,daaTemp,daaDQ,daaDT,iDoDQ,iSPlineType,  &
    raRAmt,raRTemp,raRPress,raRPartPress,pProf, iaP1,iaP2,raP1,raP2,  &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)


INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
NO TYPE, INTENT(IN OUT)                  :: rFileStart
INTEGER, INTENT(IN OUT)                  :: iRefLayer
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: kFrStep
REAL, INTENT(IN OUT)                     :: raPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPartPres
INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN OUT)                  :: iU
DOUBLE PRECISION, INTENT(IN OUT)         :: daaTemp(kMaxPts,kProfLayer)
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iDoDQ
NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
REAL, INTENT(IN OUT)                     :: raRAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

INTEGER :: iErr
INTEGER :: iSplineType,iProfileLayers
REAL :: rFileStartFr

REAL :: raPartPress(kProfLayer)

REAL :: raRPartPress(kProfLayer)


! the Matlab weights








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

IF ((ABS(iSplineType) == 2) .AND. (iActualTag == 20)) THEN
!! warning so far only CKD1,6,25 has these tables ....
  iMethod = 1        !! use compressed data for 605-2830 cm-1
  iDefault = +1      !! new method, uses kCompressed chunks
  
!decided to switch back to look up tables, as Matlab version uses this
  iMethod = 0        !! use look up tables
  iDefault = 0       !! use look up tables
END IF

IF (iDefault /= iMethod) THEN
  WRITE(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
END IF

CALL ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)

IF (iMethod == 0) THEN
! reads the oldstyle 2901 point lookup table which is fcn(fr,T)
! old style, been using this for years
  WRITE (kStdWarn,*) ' WV cont : using Lookup Table (fr x T) format'
  rConvFac = 1.0   !!! scale factor = 1 in the ORIG tables
  CALL f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ,  &
      rFStep,raFreq,kFrStep, raTemp,daaDQ,daaDT,daaCon)
ELSE IF (iMethod == 1) THEN
! reads the Compressed Style lookup fcn(fr,T)
! new style, since Sept 2011
  WRITE (kStdWarn,*) ' WV cont : using kCompressed Database format'
  rConvFac = 1.0E-23    !! SELF/FORN continuum coeffs in compressed
!! database were remultiplied by 1.0e+23
  CALL xCKDgases(iGasID,rFileStartFr,iTag,iActualTag, iRefLayer,iL,iU,  &
      raAmt,raRAmt,raTemp,raRTemp, iErr1,iDoDQ,pProf,iProfileLayers,  &
      daaDQ,daaDT,daaCon,iSplineType, iaP1,iaP2,raP1,raP2,  &
      iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
      iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
      iaQ21,iaQ22,raQ21,raQ22)
END IF

CALL aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt,  &
    raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)

RETURN
END SUBROUTINE AddContinuum

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! does some simple inits

SUBROUTINE ckd_init(raAmt,iMin,iMax,daaCon,iLay,iFr,rFStep,kFrStep)


REAL, INTENT(IN OUT)                     :: raAmt(kProfLayer)
INTEGER, INTENT(OUT)                     :: iMin
INTEGER, INTENT(OUT)                     :: iMax
DOUBLE PRECISION, INTENT(OUT)            :: daaCon(kMaxPts,kProfLayer)
INTEGER, INTENT(OUT)                     :: iLay
INTEGER, INTENT(OUT)                     :: iFr
REAL, INTENT(OUT)                        :: rFStep
REAL, INTENT(IN)                         :: kFrStep
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input

! output




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
DO iFr=1,kMaxPts
  DO iLay=iMin,iMax
    daaCon(iFr,iLay)=0.0
  END DO
END DO

iLay = kProfLayer
iFr = kMaxPts
rFStep = kFrStep

RETURN
END SUBROUTINE ckd_init

!************************************************************************

SUBROUTINE f_by_T_CKDfilereader(iGasID,iTag,iMin,iMax,iDoDQ,  &
    rFStep,raFreq,kFrStep, raTemp,daaDQ,daaDT,daaCon)


INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iMin
INTEGER, INTENT(IN OUT)                  :: iMax
INTEGER, INTENT(IN OUT)                  :: iDoDQ
REAL, INTENT(IN OUT)                     :: rFStep
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: kFrStep
REAL, INTENT(IN OUT)                     :: raTemp(kProfLayer)
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaCon(kMaxPts,kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input



! output
!   daaDQ     = analytic Jacobian wrt gas amount
!   daaDT     = analytic Jacobian wrt temperature
!   daaCon    = continuum coeffs (multiplied by 1)




! local vars
INTEGER :: iLay,iFr,iIOUN,iLoc,iDefault,iUseMethod,iMethod,iErr
REAL :: rMin
DOUBLE PRECISION :: der1,der2,doop
REAL :: REC,rECCount

! these are to read in the binary file
DOUBLE PRECISION :: d1,d2,df,daaCKD(kTempCKD,kFreqCKD),daTemprt(kTempCKD)
INTEGER :: iCKD,iM,iN


REAL :: a1,a2
CHARACTER (LEN=120) :: caFName

! figure out the filename
CALL CKDFileName(caFName,iGasID,iTag)

iIOUN = kTempUnit
! open file and load in data
OPEN(UNIT=iIOUN,FILE=caFname,STATUS='OLD',FORM='UNFORMATTED', IOSTAT=iErr)
IF (iErr /= 0) THEN
  WRITE(kStdErr,1080) iErr, caFname
  WRITE(kStdWarn,1080) iErr, caFname
  1080   FORMAT('ERROR! number ',I5,' opening CKD binary file:',/,A120)
  CALL DoSTOP
END IF

kTempUnitOpen=1
READ(iIOUN) d1,d2,df       !read start/stop freq, df
IF (d1 > raFreq(1)) THEN
  WRITE(kStdErr,*) 'CKD file has freqs that start too high!'
  CALL DoStop
END IF
IF (d2 < raFreq(kMaxPts)) THEN
  WRITE(kStdErr,*) 'CKD file has freqs end start too low!'
  CALL DoStop
END IF

!&&&&&&&&
READ(iIOUN) iLoc,iCKD,iM,iN
!read line type, CKD vers, # of temp,freq pts
IF (iLoc < 0) THEN
  WRITE(kStdErr,*) iLoc
  WRITE(kStdErr,*) 'need correct continuum for local lineshape !'
  CALL DoStop
END IF
IF (iCKD /= kCKD) THEN
  WRITE(kStdErr,*) iCKD,kCKD
  WRITE(kStdErr,*) 'CKD versions do not match!'
  CALL DoStop
END IF
IF (iM /= kTempCKD) THEN
  WRITE(kStdErr,*) iM,kTempCKD
  WRITE(kStdErr,*) 'Need more Temp offsets in CKD binary file'
  CALL DoStop
END IF
IF (iN /= kFreqCKD) THEN
  WRITE(kStdErr,*) iN,kFreqCKD
  WRITE(kStdErr,*) 'Read in ',iN, ' points from CKD binary file'
  WRITE(kStdErr,*) 'Always need ',kFreqCKD, ' points in the file'
  CALL DoStop
END IF

READ(iIOUN) (daTemprt(iLay),iLay=1,iM)   !read the temps

DO iLay=1,iM
!for current temp, read the continuum data for each freq
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
  WRITE(kStdErr,*) 'adding on CKD : iDefault,iMethod = ',iDefault,iMethod
END IF

! interpolate
IF (iMethod == 1) THEN
! this is what was used in SRCv1.07, which is what Scott Hannon used
! for the Fast Models
  CALL ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,  &
      kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)
  
ELSE IF (iMethod == 2) THEN
! this is new and improved version, which should be quicker than Quad
! THIS ONE has T dependnce both for self and foreign
  CALL ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,raFreq,raTemp,  &
      kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)
  
ELSE IF (iMethod == 3) THEN
! this is new and improved version, which should be quicker than Quad
! THIS ONE has T dependance ONLY for self
  CALL ComputeCKD_Linear_March2002(d1,d2,df, daTemprt,daaCKD,raFreq,raTemp,  &
      kFrStep,daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)
END IF

RETURN
END SUBROUTINE f_by_T_CKDfilereader

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! this subroutine completes the CKD calcs, and d/dT,d/dQ if needed

SUBROUTINE aux_ckd(iMin,iMax,iGasID,raPress,raPartPress,raTemp,raAmt,  &
    raFreq,kFrStep,rConvFac,daaCon,daaTemp,daaDQ,daaDT,iDoDQ)


INTEGER, INTENT(IN)                      :: iMin
INTEGER, INTENT(IN)                      :: iMax
INTEGER, INTENT(IN OUT)                  :: iGasID
REAL, INTENT(IN OUT)                     :: raPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPartPres
REAL, INTENT(IN)                         :: raTemp(kProfLayer)
REAL, INTENT(IN)                         :: raAmt(kProfLayer)
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: kFrStep
REAL, INTENT(IN)                         :: rConvFac
DOUBLE PRECISION, INTENT(OUT)            :: daaCon(kMaxPts,kProfLayer)
DOUBLE PRECISION, INTENT(OUT)            :: daaTemp(kMaxPts,kProfLayer)
DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iDoDQ
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iMin,Imax tell over what layers to do the CKD calcs (where gasamt nonzero!)
! input
DOUBLE PRECISION :: !! ckd coeffs


REAL :: raPartPress(kProfLayer)

! output
DOUBLE PRECISION :: !! ckd ODs



! local variables
INTEGER :: iLay,iFr,iErr
REAL :: rFStep,rMin
DOUBLE PRECISION :: der1,der2,doop
REAL :: REC,rECCount
REAL :: raMult(kProfLayer)
REAL :: a1,a2

rMin=1.0E10
rECCount=0.0
DO iLay=iMin,iMax
  DO iFr=1,kMaxPts
    IF (daaCon(iFr,iLay) < 0.0) THEN
!            write(kStdWarn,*) 'continuum < 0 for iGasID,ilay,iFr',
!     $                         iGasID,iLay,iFr,daaCon(iFr,iLay)
      daaCon(iFr,iLay) = 0.0
      rECCount = rECCount + 1.0
      IF (daaCon(iFr,iLay) < rMin) THEN
        rMin = daaCon(iFr,iLay)
      END IF
    END IF
  END DO
END DO
IF (rECCount > 0) THEN
  WRITE(kStdWarn,*) ' >>> WARNING !!! Continuum for iGasID ',iGasID,' < 0 for ',INT(rECCount),' points'
END IF

! add the continuum abs coeffs to the k-compressed coeffs raaCon = daaTemp
! keeping track of whether the added values are all greater than 0.0
IF (iGasID == kNewGasLo) THEN      !CSelf
  DO iLay =iMin,iMax
    raMult(iLay)=raAmt(iLay)*kAvog*raPartPress(iLay)
  END DO
ELSE IF (iGasID == kNewGasHi) THEN      !CFor
  DO iLay =iMin,iMax
    raMult(iLay)=raAmt(iLay)*kAvog*(raPress(iLay)-raPartPress(iLay))
  END DO
END IF

DO iLay =iMin,iMax
  a1 = raMult(iLay)*296.0/raTemp(iLay)
  a2 = kPlanck2/2/raTemp(iLay)
  DO iFr=1,kMaxPts
    daaTemp(iFr,iLay) = rConvFac * daaCon(iFr,iLay)*a1*raFreq(iFr)*  &
        TANH(a2*raFreq(iFr))
  END DO
END DO

IF (rECCount > 0.5) THEN
  REC=REC/rECCount
  WRITE(kStdWarn,*)'Error in CKD data!!! Some values negative!!!'
  WRITE(kStdWarn,*)'and reset to 0.0'
  WRITE(kStdWarn,*) rECCount,' values have avg value ',REC
  WRITE(kStdWarn,*) 'min negative value in 10000*100 = ',rMin
  IF (ABS(rMin) > 1.0E-7) THEN
    iErr=1
    WRITE(kStdErr,*)'Error in CKD data!!! Some values negative!!!'
    CALL DoSTOP
  END IF
END IF

IF (iDoDQ >= -1) THEN
  IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
    
!!the gas amount jacobians
    IF (iGasID == 101) THEN   !!!there is a factor of 2
      DO iLay=iMin,iMax
        DO iFr=1,kMaxPts
          daaDQ(iFr,iLay) = 2 * daaTemp(iFr,iLay)/raAmt(iLay)
        END DO
      END DO
    ELSE IF (iGasID == 102) THEN
      DO iLay=iMin,iMax
        DO iFr=1,kMaxPts
          daaDQ(iFr,iLay) = daaTemp(iFr,iLay)/raAmt(iLay)
        END DO
      END DO
    END IF
  END IF
  
  IF ((kActualJacs == -1) .OR. (kActualJacs == 30). OR.  &
        (kActualJacs == 32) .OR.  &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
!!this is the temperature jacobian
    DO iLay=iMin,iMax
      DO iFr=1,kMaxPts
        doop = TANH(kPlanck2*raFreq(iFr)/2/raTemp(iLay))
        der1 = -(296/raTemp(iLay)**2)* doop
        der1 = der1-(296*kPlanck2*raFreq(iFr)/2/(raTemp(iLay)**3))/  &
            (COSH(kPlanck2*raFreq(iFr)/2/raTemp(iLay))**2)
        der1 = der1*daaCon(iFr,iLay)
        
        der2 = (296/raTemp(iLay))* doop
        der2 = der2*daaDT(iFr,iLay)
        
        daaDT(iFr,iLay) = raMult(iLay)*raFreq(iFr)*(der1+der2)*rConvFac
      END DO
    END DO
  END IF
  
END IF

RETURN
END SUBROUTINE aux_ckd

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, linear in temperature
! does temperature interpolation for SELF and for foreign
! as this is the new CKD UMBC thingies

SUBROUTINE ComputeCKD_Linear(d1,d2,df,daTemprt,daaCKD,  &
    raFreq,raTemp,kFrStep, daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)


DOUBLE PRECISION, INTENT(IN OUT)         :: d1
DOUBLE PRECISION, INTENT(IN OUT)         :: d2
DOUBLE PRECISION, INTENT(IN OUT)         :: df
DOUBLE PRECISION, INTENT(IN)             :: daTemprt(kTempCKD)
DOUBLE PRECISION, INTENT(IN)             :: daaCKD(kTempCKD,kFreqCKD)
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: kFrStep
DOUBLE PRECISION, INTENT(OUT)            :: daaCon(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iDoDQ
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(OUT)            :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iMin
INTEGER, INTENT(IN)                      :: iMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop





! these were from the CKD binary file


! local variables
INTEGER :: iLay,iFr,iL,IF,iFloor
INTEGER :: iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp
DOUBLE PRECISION :: a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2
! this is for CKD1 temeperature dependant multipliers from 780-980 cm-1
INTEGER :: iNpts,iMult
DOUBLE PRECISION :: psT(kMaxLayer),psK(kMaxLayer)
DOUBLE PRECISION :: pfT(kMaxLayer),pfK(kMaxLayer)
DOUBLE PRECISION :: xNum(kMaxLayer)
DOUBLE PRECISION :: daaCKD1Mult(kMaxPts,kProfLayer)

dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

!      iMult = -1
!      IF ((kCKD .EQ. 1) .AND.
!     $  (raFreq(1) .LT. 979.0) .AND. (raFreq(kMaxPts) .GE. 780)) THEN
!        iMult = +1
!        write(kStdWarn,*) 'CKD 1 Temperature correction for chunk  ',raFreq(1)
!        CALL GetCKD1Mult(iGasID,raTemp,raFreq,
!     $                   iNpts,xNum,psT,psK,pfT,pfK,daaCKD1Mult)
!        END IF

! for each freq point in raFreq, find where nearest low CKD freq grid point is
DO iFr=1,kMaxPts
  iaFrIndex(iFr)=1 + iFloor(REAL((raFreq(iFr)-d1)/df))
  daFrDelta(iFr)=raFreq(iFr)*1D0-(d1+(iaFrIndex(iFr)-1)*df)
END DO

! for each layer temp, find where the nearest low CKD tempr grid point is
DO iLay=iMin,iMax
  iaTempIndex(iLay)=1 + iFloor(REAL((raTemp(iLay)-daTemprt(1))/dTemp))
END DO

IF (iDoDQ < -1) THEN         !no need to do temp jacobians
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
!find out line that goes thru the 2 "x(j)" points to give "y(j)"
        x2 = daTemprt(iL)
        x3 = daTemprt(iL+1)
        
        a = (y3-y2)/(x3-x2)
        b = y3-a*x3
        
        x = raTemp(iLay)
        daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
        
      END DO
    END DO
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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

IF (iMult > 0) THEN
  DO iLay = iMin,iMax
    DO iFr = 1,kMaxPts
      daaCon(iFr,iLay) = daaCon(iFr,iLay) * daaCKD1Mult(iFr,iLay)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE ComputeCKD_Linear

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, linear in temperature
! only does temperature interpolation for SELF, and not for foreign
! as this was the CKD 0,2.1,2.3,2.4 thingies

SUBROUTINE ComputeCKD_Linear_March2002(d1,d2,df,daTemprt,daaCKD,  &
    raFreq,raTemp,kFrStep, daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)


DOUBLE PRECISION, INTENT(IN OUT)         :: d1
DOUBLE PRECISION, INTENT(IN OUT)         :: d2
DOUBLE PRECISION, INTENT(IN OUT)         :: df
DOUBLE PRECISION, INTENT(IN)             :: daTemprt(kTempCKD)
DOUBLE PRECISION, INTENT(IN)             :: daaCKD(kTempCKD,kFreqCKD)
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: kFrStep
DOUBLE PRECISION, INTENT(OUT)            :: daaCon(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iDoDQ
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(OUT)            :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iMin
INTEGER, INTENT(IN)                      :: iMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep  = kFreqStep
! iMin,iMax = layers over which to loop





! these were from the CKD binary file


! local variables
INTEGER :: iLay,iFr,iL,IF,iFloor
INTEGER :: iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp
DOUBLE PRECISION :: a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
DO iFr=1,kMaxPts
  iaFrIndex(iFr)=1 + iFloor(REAL((raFreq(iFr)-d1)/df))
  daFrDelta(iFr)=raFreq(iFr)*1D0-(d1+(iaFrIndex(iFr)-1)*df)
END DO

! for each layer temp, find where the nearest low CKD tempr grid point is
DO iLay=iMin,iMax
  iaTempIndex(iLay)=1 + iFloor(REAL((raTemp(iLay)-daTemprt(1))/dTemp))
END DO

IF (iDoDQ < -1) THEN         !no need to do temp jacobians
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
!find out line that goes thru the 2 "x(j)" points to give "y(j)"
        x2 = daTemprt(iL)
        x3 = daTemprt(iL+1)
        
        a = (y3-y2)/(x3-x2)
        b = y3-a*x3
        
        x = raTemp(iLay)
        daaCon(iFr,iLay) = a*x + b     !this is temp dependance!
        
      END DO
    END DO
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        
        daaCon(iFr,iLay)=y2    !this is temp dependance!
        
      END DO
    END DO
  END IF
  
ELSE          !!!!!!!!!do the temp jacobians
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        
        daaCon(iFr,iLay)=y2    !this is temp dependance!
        daaDT(iFr,iLay)=0.0D0  !this is temp jacobian!
        
      END DO
    END DO
  END IF
END IF

RETURN
END SUBROUTINE ComputeCKD_Linear_March2002
!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses interpolations : linear in freq, quadratic in temperature

SUBROUTINE ComputeCKD_Quadratic(d1,d2,df,daTemprt,daaCKD,  &
    raFreq,raTemp,kFrStep, daaCon,iDoDQ,daaDQ,daaDT,iGasID,iMin,iMax)


DOUBLE PRECISION, INTENT(IN OUT)         :: d1
DOUBLE PRECISION, INTENT(IN OUT)         :: d2
DOUBLE PRECISION, INTENT(IN OUT)         :: df
DOUBLE PRECISION, INTENT(IN)             :: daTemprt(kTempCKD)
DOUBLE PRECISION, INTENT(IN)             :: daaCKD(kTempCKD,kFreqCKD)
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: kFrStep
DOUBLE PRECISION, INTENT(OUT)            :: daaCon(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iDoDQ
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(OUT)            :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iMin
INTEGER, INTENT(IN)                      :: iMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop





! these were from the CKD binary file


! local variables
INTEGER :: iLay,iFr,iL,IF,iFloor
INTEGER :: iaFrIndex(kMaxPts),iaTempIndex(kProfLayer)
DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp
DOUBLE PRECISION :: a,b,c,x1,x2,x3,y1,y2,y3,t1,t2,t3,t4,x,z1,z2

dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
DO iFr=1,kMaxPts
  iaFrIndex(iFr)=1 + iFloor(REAL((raFreq(iFr)-d1)/df))
  daFrDelta(iFr)=raFreq(iFr)*1D0-(d1+(iaFrIndex(iFr)-1)*df)
END DO

! for each layer temp, find where the nearest low CKD tempr grid point is
DO iLay=iMin,iMax
  iaTempIndex(iLay)=1 + iFloor(REAL((raTemp(iLay)-daTemprt(1))/dTemp))
END DO

IF (iDoDQ < -1) THEN         !no need to do temp jacobians
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y1=daaCKD(iL-1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL-1,IF+1)-daaCKD(iL-1,IF))/df
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        
        daaCon(iFr,iLay)=y2    !this is temp dependance!
        
      END DO
    END DO
  END IF
  
ELSE          !!!!!!!!!do the temp jacobians
  IF (iGasID == kNewGasLo) THEN
!self continuum has temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)!closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq index lower than raFreq(iFr)
        
        y1=daaCKD(iL-1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL-1,IF+1)-daaCKD(iL-1,IF))/df
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        y3=daaCKD(iL+1,IF) +  &
            daFrDelta(iFr)*(daaCKD(iL+1,IF+1)-daaCKD(iL+1,IF))/df
        
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
    
  ELSE IF (iGasID == kNewGasHi) THEN
!foreign continuum has no temp dependance
    DO iLay=iMin,iMax
      iL=iaTempIndex(iLay)  !closest temp index lower than raTemp(iLay)
      DO iFr=1,kMaxPts
        IF=iaFrIndex(iFr)   !closest freq indexlower than raFreq(iFr)
        
        y2=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
        
        daaCon(iFr,iLay)=y2    !this is temp dependance!
        daaDT(iFr,iLay)=0.0D0  !this is temp jacobian!
        
      END DO
    END DO
  END IF
END IF

RETURN
END SUBROUTINE ComputeCKD_Quadratic

!************************************************************************
! this subroutine computes the CKD coeffs in the data file in temp,freq
! this uses spline interpolations ........ blooody slow

SUBROUTINE ComputeCKD_SlowSpline(d1,d2,df,daTemprt,daaCKD,  &
    raFreq,raTemp,kFrStep, daaCon,iDoDQ,daaDQ,daaDT,iMin,iMax)


DOUBLE PRECISION, INTENT(IN OUT)         :: d1
DOUBLE PRECISION, INTENT(IN OUT)         :: d2
DOUBLE PRECISION, INTENT(IN OUT)         :: df
DOUBLE PRECISION, INTENT(IN)             :: daTemprt(kTempCKD)
DOUBLE PRECISION, INTENT(IN)             :: daaCKD(kTempCKD,kFreqCKD)
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: kFrStep
DOUBLE PRECISION, INTENT(IN OUT)         :: daaCon(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iDoDQ
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN)                      :: iMin
INTEGER, INTENT(IN)                      :: iMax
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! daaCon    = output continuum coeffs (normalised)
! daaDQ     = analytic Jacobian wrt gas amount
! daaDT     = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
! raFreq    = wavenumber array
! kFrStep   = kFreqStep
! iMin,iMax = layers over which to loop





! these were from the CKD binary file


! local variables
INTEGER :: iFr,iL,IF,iFloor
INTEGER :: iaFrIndex(kMaxPts)
DOUBLE PRECISION :: daFrDelta(kMaxPts),dTemp,daC(kTempCKD)
DOUBLE PRECISION :: dyp1,dypn,daY2(kTempCKD),daWork(kTempCKD)

!     Assign values for interpolation
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
dYP1=1.0E+16
dYPN=1.0E+16

dTemp=daTemprt(2)-daTemprt(1)       !temperature spacing in CKD file

! for each freq point in raFreq, find where nearest low CKD freq grid point is
DO iFr=1,kMaxPts
  iaFrIndex(iFr)=1 + iFloor(REAL((raFreq(iFr)-d1)/df))
  daFrDelta(iFr)=raFreq(iFr)*1D0-(d1+(iaFrIndex(iFr)-1)*df)
END DO

DO iFr=1,kMaxPts
  IF=iaFrIndex(iFr)   !index of closest freq lower than raFreq(iFr)
  DO iL=1,kTempCKD
    daC(iL)=daaCKD(iL,IF) + daFrDelta(iFr)*(daaCKD(iL,IF+1)-daaCKD(iL,IF))/df
  END DO
  CALL DSPLY2(daTemprt,daC,kTempCKD,dYP1,dYPN,daY2,daWork)
  DO iL=iMin,iMax
    CALL DSPLIN(daTemprt,daC,daY2,kTempCKD,DBLE(raTemp(iL)),daaCon(iFr,iL))
  END DO
END DO

RETURN
END SUBROUTINE ComputeCKD_SlowSpline

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

SUBROUTINE CrossSectionOLD(iCount,iGasID,iRefLayer,iL,iU,kFrStep,  &
    daaTemp,raVTemp,iVTSet,raFreq,iErr,caXsecName,  &
    raTAmt,raTTemp,raTPress,raTPart,iaCont, daaDQ,daaDT,iDoDQ)


INTEGER, INTENT(IN OUT)                  :: iCount
INTEGER, INTENT(IN)                      :: iGasID
INTEGER, INTENT(IN OUT)                  :: iRefLayer
INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN OUT)                  :: iU
REAL, INTENT(IN)                         :: kFrStep
DOUBLE PRECISION, INTENT(OUT)            :: daaTemp(kMaxPts,kProfLayer)
REAL, INTENT(OUT)                        :: raVTemp(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iVTSet
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(OUT)                     :: iErr
CHARACTER (LEN=80), INTENT(IN OUT)       :: caXsecName
REAL, INTENT(IN)                         :: raTAmt(kProfLayer)
REAL, INTENT(IN)                         :: raTTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTPress(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTPart(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaCont(kMaxGas)
DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(OUT)            :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iDoDQ
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: REC,rECCount

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
DO iLay=1,MXXMOL
  iaMolid(iLay)=0
END DO
iNmol=1
iaMolid(1)=iGasID

iLay = kProfLayer

IF (iErr <= 0) THEN
! set the vertical temperature profile if iCount=1 (first profile read)
  IF ((iVTSet < 0) .AND. (iGasID <= kGasComp)) THEN
    WRITE(kStdWarn,*) 'Setting vertical temp profile ...'
    DO iLay=1,kProfLayer
      raVTemp(iLay)=raTTemp(iLay)
    END DO
  END IF
! if previous profiles have been read in, check to make sure the
! temperature profiles are the same!!!!
  IF ((iVTSet > 0) .AND. (iGasID <= kGasComp)) THEN
    WRITE(kStdWarn,*) 'Checking the vertical temp profile ...'
    DO iLay=1,kProfLayer
      rCheckTemp=raTTemp(iLay)-raVTemp(iLay)
      IF (ABS(rCheckTemp) >= 1.0E-3) THEN
        WRITE(kStdWarn,*) 'Warning!!Temp profiles do not match!!!'
        WRITE(kStdWarn,*) 'Gas#,layer, gastemp, vertical temp = '
        WRITE(kStdWarn,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay)
        WRITE(kStdErr,*) 'Warning!!Temp profiles do not match!!!'
        WRITE(kStdErr,*) 'Gas#,layer, gastemp, vertical temp = '
        WRITE(kStdErr,*) iCount,iLay,raTTemp(iLay),raVTemp(iLay)
      END IF
    END DO
  END IF
END IF

IF (iErr < 0) THEN
! convert amount from k.moles/cm^2 to molecules/cm^2
  DO iLay=1,kProfLayer
    raXamnt(iLay)=raTAmt(iLay)*AVOG
    raAdjustProf(iLay)=raTTemp(iLay)
    ENDDO
      
! make sure that the relevant area in the cross section matrix is zeroed
      DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
          raaXsec(iFr,iLay)=0.0
        END DO
      END DO
      
! if we need to do jacobian calculations using splines,
! initialize the matrices here
      IF (kJacobian > 0) THEN
        IF (iDoDQ > 0) THEN
          DO iLay=1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
              raaSplineDQ(iFr,iLay)=0.0
              raaSplineDT(iFr,iLay)=0.0
            END DO
          END DO
        ELSE IF (iDoDQ < 0) THEN
          DO iLay=1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
              raaSplineDT(iFr,iLay)=0.0
            END DO
          END DO
        END IF
      END IF
      
      
! this call calculates ONLY the cross sections
!            OR
! as idQT=1 this call calculates the cross sections, AND d/dq,d/dT!!!, saving
! them in raaSplineDQ ,raaSplineDT respectively
      
      iLay = kProfLayer
      iFr = kMaxPts
      rFStep = kFrStep
      
      WRITE(kStdErr,*) 'CALXSC no longer supported!!!!'
      CALL DoStop
!ccccc        CALL CALXSC(caXsecName,iFr,raFreq,rFStep,
!ccccc     $      iLay, raAdjustProf,raXamnt,iNmol,iaMolid,
!ccccc     $      raaXsec,iWhichUsed,
!ccccc     $      idQT, raaSplineDQ ,raaSplineDT, iGasID, iDoDQ)
      
! since CALXSC computed the cross-section and saved it in a 3d matrix,
! pull out the relevant cross section data into a 2d matrix and add it to
! the cumulative abscoeff matrix daaTemp
! first check whether the CALXSC routine actually found a matching
! gas+wavenumber region
      rMin=1.0E10
      REC=0.0
      rECCount=0.0
      IF ((iWhichUsed > 0) .AND. (iWhichUsed <= iNmol)) THEN
        WRITE(kStdWarn,1000) iCount,iGasID
        1000          FORMAT('adding on cross section for GAS(',I2,') = ',I2)
        DO iLay=1,kProfLayer
          DO iFr=1,kMaxPts
            IF (raaXsec(iFr,iLay) < 0.0) THEN
! flag this error, as the cross section values should all be > 0
              iErr=1
              REC=REC+raaXsec(iFr,iLay)
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
        1010              FORMAT ('no XSEC contribution due to GAS(',I2,') = ',I2)
      END IF
      IF (rECCount > 0.5) THEN
        REC=REC/rECCount
        WRITE(kStdWarn,*) 'Error in XSEC data!!! Some values negative!'
        WRITE(kStdWarn,*) 'and reset to 0.0'
        WRITE(kStdWarn,*) 'rECCount values have avg value ',REC
        WRITE(kStdWarn,*) 'min negative value in 10000*100 = ',rMin
        IF (ABS(rMin) > 1.0E-7) THEN
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
        DO iLay=1,kProfLayerJac
          DO iFr=1,kMaxPtsJac
            daaDQ(iFr,iLay)=raaSplineDQ(iFr,iLay)*avog
            daaDT(iFr,iLay)=raaSplineDT(iFr,iLay)
          END DO
        END DO
      ELSE IF (iDoDQ < 0) THEN
        DO iLay=1,kProfLayerJac
          DO iFr=1,kMaxPtsJac
            daaDT(iFr,iLay)=raaSplineDT(iFr,iLay)
          END DO
        END DO
      END IF
    END IF
    
    RETURN
  END SUBROUTINE CrossSectionOLD
  
!************************************************************************
! this gets the temperature dependant multipliers for CKD1
  
  SUBROUTINE GetCKD1Mult(iGasID,raTemp,raFreq,  &
      iNpts,dxAERI,psT,psK,pfT,pfK, daaCKD1Mult)
  
  
  INTEGER, INTENT(IN OUT)                  :: iGasID
  REAL, INTENT(IN)                         :: raTemp(kProfLayer)
  REAL, INTENT(IN)                         :: raFreq(kMaxPts)
  INTEGER, INTENT(OUT)                     :: iNpts
  DOUBLE PRECISION, INTENT(IN)             :: dxAERI(kMaxLayer)
  DOUBLE PRECISION, INTENT(IN OUT)         :: psT(kMaxLayer)
  DOUBLE PRECISION, INTENT(IN OUT)         :: psK(kMaxLayer)
  DOUBLE PRECISION, INTENT(IN OUT)         :: pfT(kMaxLayer)
  DOUBLE PRECISION, INTENT(IN OUT)         :: pfK(kMaxLayer)
  NO TYPE, INTENT(IN OUT)                  :: daaCKD1Mul
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  
! input params
  
  
! output params
! these are read from the data file
  
  
  
  
! this are the temperature corrected coeffs
  DOUBLE PRECISION :: daaCKD1Mult(kMaxPts,kProfLayer)
  
! local vars
  CHARACTER (LEN=80) :: FNAM
  INTEGER :: iL,iFr,iN,iIOUN,iErr
  DOUBLE PRECISION :: daXT(kMaxLayer),daXK(kMaxLayer)
  DOUBLE PRECISION :: dT,daX(kMaxLayer),daM(kMaxPts),x,y,daFreq(kMaxPts)
  
  FNAM = '/home/sergio/KCARTA/INCLUDE/ckd1mult.dat'
  
  iIOUN = kTempUnit
  OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='FORMATTED', IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE(kStdErr,1010) IERR, FNAM
    1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
    CALL DoSTOP
  END IF
  kTempUnitOpen=1
  
  iNpts = 0
  IF (iGasID == kNewGasLo) THEN
    10     CONTINUE
    iNpts = iNpts + 1
    READ(iIOUN,*,ERR=20) dxAERI(iNpts),daXT(iNpts),daXK(iNpts),x,y
!        print *,iNpts,dxAERI(iNpts),daXT(iNpts),daXK(iNpts)
    GO TO 10
  ELSE IF (iGasID == kNewGasHi) THEN
    11     CONTINUE
    iNpts = iNpts + 1
    READ(iIOUN,*,ERR=20) dxAERI(iNpts),x,y,daXT(iNpts),daXK(iNpts)
!        print *,iNpts,dxAERI(iNpts),daXT(iNpts),daXK(iNpts)
    GO TO 11
  END IF
  
  20   CONTINUE
  CLOSE(iIOUN)
  kTempUnitOpen=-1
  
  DO iFr = 1,kMaxPts
    daFreq(iFr) = raFreq(iFr) * 1.0D0
  END DO
  
  DO iL = 1,kProfLayer
    dT = raTemp(iL)*1.0D0
    IF (dT >= 296.0D0) THEN  !!! no need to adjust things
      DO iFr = 1,kMaxPts
        daaCKD1Mult(iFr,iL) = 1.0D0
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
!        print *,iL,iGasID,iNpts,kProfLayer,dT,raFreq(1),daaCKD1Mult(1,iL)
  END DO
  
  RETURN
END SUBROUTINE GetCKD1Mult

!************************************************************************
