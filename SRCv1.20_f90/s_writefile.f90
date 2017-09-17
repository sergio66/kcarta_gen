! Copyright 1997
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:44
 
! University of Maryland Baltimore County
! All Rights Reserved

!***********************************************************************
!**************** MAIN BINARY OUTPUT SUBROUTINES ARE HERE **************
!***********************************************************************
! this subroutine writes out the main header for each results section
! this is a major change from before
! note the program does not keep on opening and closing the file

SUBROUTINE wrtout_head(iIOUN,caOutName,rFrLow,rFrHigh,rDelta,  &
    iMainType,iSubMainType,iNumberOut)


INTEGER, INTENT(IN OUT)                  :: iIOUN
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
REAL, INTENT(OUT)                        :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(OUT)                        :: rDelta
INTEGER, INTENT(IN)                      :: iMainType
NO TYPE, INTENT(IN OUT)                  :: iSubMainTy
INTEGER, INTENT(OUT)                     :: iNumberOut
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIOUN      = unit file number
! caOutName  = binary output file name
! rFrLow     = lowest wavenumber to be output
! rFrHigh    = highest wavenumber to be output
! rDelta     = point spacing
! iMainType  = 1 for path spectra                   iAtmNumber for Jacobians
!              2 for MP spectra                     iAtmNumber for fluxes
!              3 for radiances
! iSubMainType = kLayer2Sp (-2,-1,+1,+2,+3,+4) for path     GAS ID (1..63) for GasJac
!              = kLayer2Sp (-2,-1,+1,+2,+3,+4) for MP       0           for Temp Jac
!              = iAtmNumber for radiances                 -10         for wgt fcn
!                                                         -20         for surf Jac
!              = +1 for upward flux, -1 for downward flux
! iNumberOut   = number of the relevant spectra to look for



INTEGER :: iSubMainType

IF (ABS((rFrHigh-rFrLow)-(rDelta*kMaxPts)) >= 2*rDelta) THEN
  WRITE(kStdErr,*) 'Wow! need (rFrHigh-rFrLow)-(rDelta*kMaxPts)) < 2*rDelta'
  WRITE(kStdErr,*) 'rFrHigh,rFrLow,rDelta,kMaxPts = '
  WRITE(kStdErr,*) rFrHigh,rFrLow,rDelta,kMaxPts
  CALL DoSTOP
END IF
WRITE(kStdWarn,*) 'dump out : ',iNumberOut,' for iMain=',iMainType

IF (kLongOrShort /= 0) THEN
  WRITE(iIOUN) iMainType,iSubMainType,iNumberOut
  WRITE(iIOUN) kMaxPts,rFrLow,rFrHigh,rDelta
END IF

RETURN
END SUBROUTINE wrtout_head
!************************************************************************
! this subroutine writes out the results
! this is a major change from before
! note the program does not keep on opening and closing the file

SUBROUTINE wrtout(iIOUN,caOutName,raFreq,raInten)


INTEGER, INTENT(IN OUT)                  :: iIOUN
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIOUN      = unit file number
! raInten    = array containing computed data (radiances,spectra,jacobians ..)
! raFreq    = array containing wavenumbers
! caOutName  = binary output file name





! local variables
INTEGER :: iInt

WRITE(iIOUN) (raInten(iInt),iInt=1,kMaxPts)

RETURN
END SUBROUTINE wrtout

!************************************************************************
!     these are to dump out Planck multipliers (for NLTE)
!************************************************************************
! this subroutine dumps out the Planck Modifiers for UA

SUBROUTINE  DumpPlanckUA(iAtm,iUpper,caPlanckFile,  &
    raFreq,rDelta,raaUpperPlanckCoeff)


INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iUpper
NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input variables

REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)

CHARACTER (LEN=80) :: caPlanckFile

! local variables
INTEGER :: iFr,iL,iLay,iIOUN,iBeta
REAL :: raTemp(kMaxPts)

iIOUN = kStdPlanck
CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,1,iUpper)

WRITE(kStdWarn,*) 'subroutine UADumpPlanck is dumping out UAPlanck modifiers'
WRITE(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iUpper,' UA layers'

DO iL = 1,iUpper
  iLay = iL
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raaUpperPlanckCoeff(iFr,iLay)
  END DO
  CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
END DO

RETURN
END SUBROUTINE  DumpPlanckUA

!************************************************************************
! this subroutine dumps out the Planck Modifiers

SUBROUTINE  DumpPlanck(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,  &
    raFreq,rDelta,raaPlanckCoeff)


INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input variables
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

CHARACTER (LEN=80) :: caPlanckFile

! local variables
INTEGER :: iFr,iL,iLay,iIOUN,iBeta
REAL :: raTemp(kMaxPts)

iIOUN = kStdPlanck
CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,1,iaNumLayer(iAtm))

WRITE(kStdWarn,*) 'subroutine DumpPlanck is dumping out Planck modifiers'
WRITE(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iaNumLayer(iAtm),' layers'

DO iL = 1,iaNumLayer(iAtm)
  iLay = iaaRadLayer(iAtm,iL)
  iBeta = MOD(iLay,kProfLayer)
  IF (iBeta == 0) THEN
    iBeta = kProfLayer
  END IF
  iLay = iBeta
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raaPlanckCoeff(iFr,iLay)
  END DO
  CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
END DO

RETURN
END SUBROUTINE  DumpPlanck

!************************************************************************
! this subroutine dumps out the Planck Modifiers, set at 1.0 (ie no NLTE yet
! in this chunk!!!!)

SUBROUTINE  DumpPlanckOne(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,  &
    raFreq,rDelta,raaPlanckCoeff)


INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input variables
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

CHARACTER (LEN=80) :: caPlanckFile

! local variables
INTEGER :: iFr,iL,iLay,iIOUN,iBeta
REAL :: raTemp(kMaxPts)

WRITE(kStdWarn,*) 'subroutine DumpPlanckOne is dumping out Planck modifiers = 1'
WRITE(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iaNumLayer(iAtm),' layers'

iIOUN = kStdPlanck
CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,1,iaNumLayer(iAtm))

DO iL = 1,iaNumLayer(iAtm)
  iLay = iaaRadLayer(iAtm,iL)
  iBeta = MOD(iLay,kProfLayer)
  IF (iBeta == 0) THEN
    iBeta = kProfLayer
  END IF
  iLay = iBeta
  DO iFr = 1,kMaxPts
    raTemp(iFr) = 1.0
  END DO
  CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
END DO

RETURN
END SUBROUTINE  DumpPlanckOne

!************************************************************************
! this subroutine dumps out the Planck Modifiers for UA as ones

SUBROUTINE  DumpPlanckUAOne(iAtm,iUpper,caPlanckFile,  &
    raFreq,rDelta,raaUpperPlanckCoeff)


INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iUpper
NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input variables

REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)

CHARACTER (LEN=80) :: caPlanckFile

! local variables
INTEGER :: iFr,iL,iLay,iIOUN,iBeta
REAL :: raTemp(kMaxPts)

iIOUN = kStdPlanck
CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,1,iUpper)

WRITE(kStdWarn,*) 'subroutine UADumpPlanckOne is dumping out UAPlanck modifiers = 1'
WRITE(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iUpper,' UA layers'

DO iL = 1,iUpper
  iLay = iL
  DO iFr = 1,kMaxPts
    raTemp(iFr) = 1.0
  END DO
  CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
END DO

RETURN
END SUBROUTINE  DumpPlanckUAOne

!************************************************************************
!      these are optical depths (individual gas or cumulative ODs)
!************************************************************************

! this function checks the appropriate set of paths to be output, to make sure
! that each path/MP/layer is being output only once! (else the bookkkeeping
! in readatmos.m becomes messed up)

INTEGER FUNCTION CheckDoubleEntry(iaaOp,iRow,iNumLay)


INTEGER, INTENT(IN)                      :: iaaOp(kMaxPrint,kPathsOut)
INTEGER, INTENT(IN OUT)                  :: iRow
INTEGER, INTENT(IN)                      :: iNumLay
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
! iaaOp       = list of paths/MP/layers to be output for each print option
! iRow        = which row to be checked
! iNumLay     = number of elements to be checked in the iRow th row of iaaOp


INTEGER :: iAns,iJ

! assume everything OK
iAns=1

! remember the entries in the row have already been sorted in ascending order,
! and so all that has to be done is to check that adjacent elements are
! different from each other
DO iJ=1,iNumLay-1
  IF (iaaOp(iRow,iJ) == iaaOp(iRow,iJ+1)) THEN
    WRITE(kStdErr,*)'Outtype',iRow,' has',iaaOp(iRow,iJ)
    WRITE(kStdErr,*)'entered more than once'
    iAns=-1
  END IF
END DO

CheckDoubleEntry=iAns

RETURN
END FUNCTION CheckDoubleEntry

!************************************************************************
! this function goes thru the list of layers to be output iaaOp(iJ,1..iEnd)
! to see if they do exist in atmosphere profile #iI, which is in
! iaaRadLayer(iI,1..iNL)

INTEGER FUNCTION OutputLayerInProfile(iaaRadLayer,iI,iNL, iaaOp,iJ,iEnd)


NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iI
INTEGER, INTENT(IN)                      :: iNL
INTEGER, INTENT(IN OUT)                  :: iaaOp(kMaxPrint,kPathsOut)
INTEGER, INTENT(IN OUT)                  :: iJ
INTEGER, INTENT(IN OUT)                  :: iEnd
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iaaOp       = list of layers to be output
! iJ          = print output option number being considered
! iEnd        = numbver of layers in print option #iJ
! iaaRadLayer = for each atmosphere, list of layers used to buid up the atm
! iI          = atmosphere being considered
! iNL         = number of layers in the iI th atmosphere

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


! local variables
INTEGER :: iaTempRad(kMixFilRows)
INTEGER :: iOK,iK
INTEGER :: DoOutputLayer

! naively assume everything OK
iOK=1

! now actually check that iaaOp(iJ,:) is in iaaRadLayer(iI,1..NL)
! first store the relevant iI'th atmosphere of iaaRadLayer(iI,1..NL)
DO iK=1,iNL
  iaTempRad(iK)=iaaRadLayer(iI,iK)
END DO
! and then sort it
CALL DoSort(iaTempRad,iNL)

! now check that the elements of iaaOp to be output, iaaOp(iJ,1..iEnd) are
! all in iaTempRad
iK=0
20   iK=iK+1
IF ((iK <= iEnd) .AND. (iOk > 0)) THEN
  iOk=DoOutputLayer(iaaOp(iJ,iK),iNL,iaTempRad)
  IF (iOK < 0) THEN
    WRITE(kStdErr,*) 'output layer#',iaaOp(iJ,iK),' not found in atm#',iI
  END IF
  GO TO 20
END IF

OutputLayerInProfile=iOK

RETURN
END FUNCTION OutputLayerInProfile

!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now output
! the individual GAS PATH spectra, according to what kLayer2Sp is set to
! kLayer2Sp = +4 : gas transmittances Layer to ground sum(j=i,1) exp(-k(j))
! kLayer2Sp = +3 : gas optical depth  Layer to ground sum(j=i,1) k(j)
! kLayer2Sp = +2 : gas transmittances Layer to space  sum(j=i,N) exp(-k(j))
! kLayer2Sp = +1 : gas optical depth  Layer to space  sum(j=i,N) (k(j))
! kLayer2Sp = -1 : gas optical depth                  k(i)
! kLayer2Sp = -2 : gas transmittances                 exp(-k(j))
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_trans_path(raFreq,rFrLow,rFrHigh, raaGasAbs,iPrinter,  &
    raTAmt,raTTemp,raTPress,raTPartPress, caOutName,  &
    iFileID, iaPath,iNp,iaOp)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iPrinter
REAL, INTENT(IN OUT)                     :: raTAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPartPre
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iaPath(kProfLayer)
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaGasAbs  = single gas abs coeffs
! iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to
! iaPath     = list of the paths corresponding to the current gas
! iNp        = total number of paths to be output
! iaOp       = list of the paths to be output





REAL :: raTPartPress(kProfLayer)

! local variables
INTEGER :: iInt,iDiv,iDp,iStart,iPath,iLay,DoOutputLayer,iIOUN
REAL :: raL2S(kMaxPts)

iIOUN = kStdkCarta

iStart=iDiv(iaPath(1),kProfLayer)

! write spectra to unformatted file
! if iPrinter=1 then have to check for valid paths
DO iLay=1,kProfLayer
! check to see if this path should be output
  iPath=iStart*kProfLayer + iLay
!        IF (iPath .NE. iaPath(iLay)) THEN
!          write(kStdErr,*) 'iPath .NE. iaPath(iLay)'
!          Call DoSTOP
!        END IF
  iDp=DoOutputLayer(iPath,iNp,iaOp)
  IF (iDp > 0) THEN
    IF (kLayer2Sp == -1) THEN
      WRITE(kStdWarn,*)'output GAS OD : iPath,P,T,A,OD = ',  &
          iPath,raTPress(iLay)*kAtm2mb,raTTemp(iLay),raTAmt(iLay),raaGasAbs(1,iLay)
      DO iInt=1,kMaxPts
        raL2S(iInt)=raaGasAbs(iInt,iLay)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    ELSE IF (kLayer2Sp == -2) THEN
      WRITE(kStdWarn,*)'outputting GAS layer trans at iPath = ',iPath
      DO iInt=1,kMaxPts
        raL2S(iInt)=EXP(-raaGasAbs(iInt,iLay))
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    ELSE IF (kLayer2Sp == 1) THEN
      WRITE(kStdWarn,*)'outputting GAS layer2sp OD at iPath = ',iPath
      CALL GasOptLayer2Space(raaGasAbs,raL2S,iLay)
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    ELSE IF (kLayer2Sp == 2) THEN
      WRITE(kStdWarn,*)'outputting GAS layer2sp trans at iPath = ',iPath
      CALL GasTranLayer2Space(raaGasAbs,raL2S,iLay)
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    ELSE IF (kLayer2Sp == 3) THEN
      WRITE(kStdWarn,*)'outputting GAS layer2gnd OD at iPath = ',iPath
      CALL GasOptLayer2Ground(raaGasAbs,raL2S,iLay)
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    ELSE IF (kLayer2Sp == 4) THEN
      WRITE(kStdWarn,*)'outputting GAS layer2gnd trans at iPath = ',iPath
      CALL GasTranLayer2Ground(raaGasAbs,raL2S,iLay)
      CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
    END IF
  END IF
END DO

!      WRITE (6,1001)caOutName
! 1001 FORMAT('Successfully saved unformatted PATH results to ',/,A80)

RETURN
END SUBROUTINE out_trans_path

!************************************************************************
! this is the generic call to MP output

SUBROUTINE DoOutputMixedPaths( iFound,iPrinter,caOutName,  &
    raFreq,rFreqStart,rFreqEnd, raaSumAbCoeff,  &
    iNpmix,iFileID,iNp,iaOp)


INTEGER, INTENT(IN OUT)                  :: iFound
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFreqStart
REAL, INTENT(IN OUT)                     :: rFreqEnd
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
INTEGER, INTENT(IN)                      :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input params


REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)


! local vars
INTEGER :: iDummy,iIpmix,DoOutputLayer

IF ((iFound > 0) .AND. (iPrinter == 2)) THEN
  
! we have a list of mixed paths to output!!!
  IF (iNp < 0) THEN
    WRITE(kStdWarn,*) 'OUTPUTTING ALL MIXED PATHS ... '
  END IF
  
  IF ((kLayer2Sp == 4) .AND. (iNp < 0)) THEN
! this indicates we want L2G transmittances for ALL mixed paths!!!
    WRITE(kStdWarn,*)'     outputting ALL L2G transmittances'
    CALL out_FASTL2Gtrans_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE IF ((kLayer2Sp == 3) .AND. (iNp < 0)) THEN
    WRITE(kStdWarn,*) '     outputting ALL L2G optical depths'
! this indicates we want L2G ods for ALL mixed paths!!!
    CALL out_FASTL2Goptdp_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE IF ((kLayer2Sp == 2) .AND. (iNp < 0)) THEN
! this indicates we want L2S transmittances for ALL mixed paths!!!
    WRITE(kStdWarn,*)'     outputting ALL L2S transmittances'
    CALL out_FASTL2Strans_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE IF ((kLayer2Sp == 1) .AND. (iNp < 0)) THEN
    WRITE(kStdWarn,*) '     outputting ALL L2S optical depths'
! this indicates we want L2S ods for ALL mixed paths!!!
    CALL out_FASTL2Soptdp_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE IF ((kLayer2Sp == -1) .AND. (iNp < 0)) THEN
! this indicates we want L2S transmittances for ALL mixed paths!!!
    WRITE(kStdWarn,*) '    outputting ALL layer optical depths'
    CALL out_FASToptdp_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE IF ((kLayer2Sp == -2) .AND. (iNp < 0)) THEN
! this indicates we want L2S transmittances for ALL mixed paths!!!
    WRITE(kStdWarn,*) '    outputting ALL layer transmittances'
    CALL out_FASTtrans_MP(raFreq,rFreqStart,rFreqEnd,  &
        raaSumAbCoeff,iPrinter, caOutName,  &
        iNpmix,iFileID)
  ELSE
! now loop over the mixed paths, outputting whichever ones are necessary
    DO iIpmix = 1,iNpmix
! see if mixed path iIpmix is set from *OUTPUT
      iDummy = -1
      iDummy = DoOutputLayer(iIpmix,iNp,iaOp)
! if the printing option=2 and iIpmix has been found in the list of
! paths to be output then go ahead and print the relevant (iIpmix th) row of
! SUM abs spectra
      IF ((iPrinter == 2) .AND. (iDummy > 0)) THEN
        CALL out_trans_MP(raFreq,rFreqStart,rFreqEnd,  &
            raaSumAbCoeff,iPrinter, caOutName,  &
            iIpmix,iNpmix,iFileID)
      END IF
! go to next MIXED PATH set by incrementing iIpmix
    END DO
  END IF
! end if (iFound==1)
END IF

RETURN
END SUBROUTINE DoOutputMixedPaths

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH transmittance spectra
! kLayer2Sp = 4  : MP transmittances Layer to ground sum(j=1,i) exp(-k(j))
! kLayer2Sp = 3  : MP optical depth Layer to ground  sum(j=1,i) (k(j))
! kLayer2Sp = 2  : MP transmittances Layer to space sum(j=i,N) exp(-k(j))
! kLayer2Sp = 1  : MP optical depth Layer to space  sum(j=i,N) (k(j))
! kLayer2Sp = -1 : MP optical depth                 k(i
! kLayer2Sp = -2 : MP transmittances                exp(-k(j)
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_trans_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iIpmix,iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iIpmix
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iIpmix     = current mixed path number
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iInt,iPath,iIOUN
REAL :: raL2S(kMaxPts)

iIOUN = kStdkCarta

! write spectra to unformatted file
! this is for mixed paths
iPath=iIpmix
WRITE(kStdWarn,*)'output MIXED path optical depths at iIpmix = ',iIpmix
IF (kLayer2Sp == -1) THEN
  DO iInt=1,kMaxPts
    raL2S(iInt)=raaSumAbs(iInt,iIpmix)
  END DO
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
ELSE IF (kLayer2Sp == -2) THEN
  DO iInt=1,kMaxPts
    raL2S(iInt)=EXP(-raaSumAbs(iInt,iIpmix))
  END DO
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
ELSE IF (kLayer2Sp == 1) THEN
  CALL MP2SpaceOptDp(raaSumAbs,raL2S,iIpmix,iNpmix)
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
ELSE IF (kLayer2Sp == 2) THEN
  CALL MP2SpaceTrans(raaSumAbs,raL2S,iIpmix,iNpmix)
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
ELSE IF (kLayer2Sp == 3) THEN
  CALL MP2GroundOptDp(raaSumAbs,raL2S,iIpmix,iNpmix)
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
ELSE IF (kLayer2Sp == 4) THEN
  CALL MP2GroundTrans(raaSumAbs,raL2S,iIpmix,iNpmix)
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
END IF

!      WRITE (6,1001)caOutName
! 1001 FORMAT('Successfully saved unformatted MP results to ',/,A80)

RETURN
END SUBROUTINE out_trans_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH L2Stransmittance spectra FAST. This assumes each atmos
! is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
! has less than 100 layers (i.e. iNpmix is not a multiple of 100),
! then the "upper layers" are filled with spectra of 0, so their
! transmittance is 1

! this assumes kLayer2Sp = 2, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASTL2Gtrans_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iPath,iIOUN,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax,iI
REAL :: raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

iIOUN = kStdkCarta

! first figure out how many 100 atmosphere layer blocks there are
iAtmBlocks=1
60   CONTINUE
IF (kProfLayer*iAtmBlocks < iNpMix) THEN
  iAtmBlocks=iAtmBlocks+1
  GO TO 60
END IF
iAmax=iAtmBlocks

! make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
! saying the last "atmosphere" has less than set number of layers
iWarn=MOD(iNpmix,kProfLayer)
IF (iWarn /= 0) THEN
  iAmax=iAmax-1
END IF

! calculate L2S and write spectra to unformatted file
iPath=0

! first do the "complete" atmospheres
DO iA=1,iAmax
! put the current atmosphere into raaTempArray
  DO iLay=1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
! compute the L2G optical depths
  DO iLay = 2,kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,kProfLayer
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=EXP(-raaTempArray(iFr,iLay))
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END DO

! then do the last, "incomplete" atmospheres
! put the current atmosphere into raaTempArray
IF (iWarn /= 0) THEN
  iA=iAmax+1
  DO iLay=1,iWarn
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
  DO iLay=iWarn+1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=0.0
    END DO
  END DO
! compute the L2G optical depth
  DO iLay=2,iWarn
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,iWarn
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=EXP(-raaTempArray(iFr,iLay))
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END IF

RETURN
END SUBROUTINE out_FASTL2Gtrans_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH L2S optical depth FAST. This assumes that in each atmos
! is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
! is found to have less than 100 layers (i.e. iNpmix is not a multiple of 100),
! then the "upper layers" are filled with spectr of 0, so their transmittance
! is 1

! this assumes kLayer2Sp = 1, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASTL2Goptdp_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iPath,iIOUN,iI,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax
REAL :: raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

iIOUN = kStdkCarta

! first figure out how many kProfLayer atmosphere layer blocks there are
iAtmBlocks=1
60   CONTINUE
IF (kProfLayer*iAtmBlocks < iNpMix) THEN
  iAtmBlocks=iAtmBlocks+1
  GO TO 60
END IF
iAmax=iAtmBlocks

! make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
! saying the last "atmosphere" has less than kProfLayer layers
iWarn=MOD(iNpmix,kProfLayer)
IF (iWarn /= 0) THEN
  iAmax=iAmax-1
END IF

! calculate L2S and write spectra to unformatted file
iPath=0

! first do the "complete" atmospheres
DO iA=1,iAmax
! put the current atmosphere into raaTempArray
  DO iLay=1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
! compute the L2S optical depths
  DO iLay = 2,kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,kProfLayer
    iPath=iPath+1
    DO iFr=1,kMaxPts
      raL2S(iFr)=raaTempArray(iFr,iLay)
    END DO
    WRITE(kStdWarn,*) 'output MIXED path spectra  at ',iPath
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END DO

! then do the last, "incomplete" atmospheres
! put the current atmosphere into raaTempArray
IF (iWarn /= 0) THEN
  iA=iAmax+1
  DO iLay=1,iWarn
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
  DO iLay=iWarn+1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=0.0
    END DO
  END DO
! compute the L2S optical depths
  DO iLay=iWarn-1,1,-1
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=2,iWarn
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=raaTempArray(iFr,iLay)
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END IF

RETURN
END SUBROUTINE out_FASTL2Goptdp_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH L2Stransmittance spectra FAST. This assumes each atmos
! is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
! has less than 100 layers (i.e. iNpmix is not a multiple of 100),
! then the "upper layers" are filled with spectra of 0, so their
! transmittance is 1

! this assumes kLayer2Sp = 2, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASTL2Strans_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iPath,iIOUN,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax,iI
REAL :: raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

iIOUN = kStdkCarta

! first figure out how many 100 atmosphere layer blocks there are
iAtmBlocks=1
60   CONTINUE
IF (kProfLayer*iAtmBlocks < iNpMix) THEN
  iAtmBlocks=iAtmBlocks+1
  GO TO 60
END IF
iAmax=iAtmBlocks

! make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
! saying the last "atmosphere" has less than set number of layers
iWarn=MOD(iNpmix,kProfLayer)
IF (iWarn /= 0) THEN
  iAmax=iAmax-1
END IF

! calculate L2S and write spectra to unformatted file
iPath=0

! first do the "complete" atmospheres
DO iA=1,iAmax
! put the current atmosphere into raaTempArray
  DO iLay=1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
! compute the L2S optical depths
  DO iLay = kProfLayer-1,1,-1
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,kProfLayer
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=EXP(-raaTempArray(iFr,iLay))
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END DO

! then do the last, "incomplete" atmospheres
! put the current atmosphere into raaTempArray
IF (iWarn /= 0) THEN
  iA=iAmax+1
  DO iLay=1,iWarn
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
  DO iLay=iWarn+1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=0.0
    END DO
  END DO
! compute the L2S optical depth
  DO iLay=iWarn-1,1,-1
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,iWarn
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=EXP(-raaTempArray(iFr,iLay))
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END IF

RETURN
END SUBROUTINE out_FASTL2Strans_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH L2S optical depth FAST. This assumes that in each atmos
! is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
! is found to have less than 100 layers (i.e. iNpmix is not a multiple of 100),
! then the "upper layers" are filled with spectr of 0, so their transmittance
! is 1

! this assumes kLayer2Sp = 1, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASTL2Soptdp_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iPath,iIOUN,iI,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax
REAL :: raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

iIOUN = kStdkCarta

! first figure out how many kProfLayer atmosphere layer blocks there are
iAtmBlocks=1
60   CONTINUE
IF (kProfLayer*iAtmBlocks < iNpMix) THEN
  iAtmBlocks=iAtmBlocks+1
  GO TO 60
END IF
iAmax=iAtmBlocks

! make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
! saying the last "atmosphere" has less than kProfLayer layers
iWarn=MOD(iNpmix,kProfLayer)
IF (iWarn /= 0) THEN
  iAmax=iAmax-1
END IF

! calculate L2S and write spectra to unformatted file
iPath=0

! first do the "complete" atmospheres
DO iA=1,iAmax
! put the current atmosphere into raaTempArray
  DO iLay=1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
! compute the L2S optical depths
  DO iLay = kProfLayer-1,1,-1
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,kProfLayer
    iPath=iPath+1
    DO iFr=1,kMaxPts
      raL2S(iFr)=raaTempArray(iFr,iLay)
    END DO
    WRITE(kStdWarn,*) 'output MIXED path spectra  at ',iPath
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END DO

! then do the last, "incomplete" atmospheres
! put the current atmosphere into raaTempArray
IF (iWarn /= 0) THEN
  iA=iAmax+1
  DO iLay=1,iWarn
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
    END DO
  END DO
  DO iLay=iWarn+1,kProfLayer
    iI=iLay+(iA-1)*kProfLayer
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=0.0
    END DO
  END DO
! compute the L2S optical depths
  DO iLay=iWarn-1,1,-1
    DO iFr=1,kMaxPts
      raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+ raaTempArray(iFr,iLay+1)
    END DO
  END DO
! print them out
  DO iLay=1,iWarn
    iPath=iPath+1
    WRITE(kStdWarn,*) 'output MIXED path spectra at ',iPath
    DO iFr=1,kMaxPts
      raL2S(iFr)=raaTempArray(iFr,iLay)
    END DO
    CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
  END DO
END IF

RETURN
END SUBROUTINE out_FASTL2Soptdp_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH optical depth FAST.

! this assumes kLayer2Sp = -1, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASToptdp_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iIOUN,iI,iFr
REAL :: raL2S(kMaxPts)

iIOUN = kStdkCarta

DO iI=1,iNpMix
  DO iFr=1,kMaxPts
    raL2S(iFr)=raaSumAbs(iFr,iI)
  END DO
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
END DO

RETURN
END SUBROUTINE out_FASToptdp_MP

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now output
! the MIXED PATH transmittance spectra FAST.

! this assumes kLayer2Sp = -2, iOp = -1
! check to see if we want the raw spectra, or kLayer2Space

SUBROUTINE out_FASTtrans_MP(raFreq,rFrLow,rFrHigh, raaSumAbs,iPrinter,  &
    caOutName, iNpmix,iFileID)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN)                         :: raaSumAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPrinter
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = array containin all the frequencies in the current 25 cm-1 block
! rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1,
!                  these need not correspond to 1,10000)
! raaSumAbs  = mixed path sum of the abs coeffs
! iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
! iNpmix     = number of mixed paths
! iFileID       = which of the 25 cm-1 k-comp files is being processed
! caOutName  = name of binary file that output goes to





! local variables
INTEGER :: iI,iFr,iIOUN
REAL :: raL2S(kMaxPts)

iIOUN = kStdkCarta

DO iI=1,iNpMix
  DO iFr=1,kMaxPts
    raL2S(iFr)=EXP(-raaSumAbs(iFr,iI))
  END DO
  CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
END DO

RETURN
END SUBROUTINE out_FASTtrans_MP

!************************************************************************
! this subroutine does the layer to space optical depths by doing
! the calculation between frequency indices from layer iLay-->100
! the result is stored in raL2S
! this is mainly used by the OpticalLayer2Space printout routines in rdprof.f

SUBROUTINE GasOptLayer2Space(raaGasAbs,raL2S,iLay)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kProfLayer)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iLay
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay      = current layer number
! raaGasAbs = matrix containing the raw abs coeffs
! raL2S     = array containing the layer-to-space results



INTEGER :: iI,iJ

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=iLay,kProfLayer
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

RETURN
END SUBROUTINE GasOptLayer2Space

!************************************************************************
! this subroutine does the layer to space transmittances for gas iG by doing
! the calculation between frequency indices from layer iLay-->100
! the result is stored in raL2S
! this is mainly used by the OpticalLayer2Space printout routines in rdprof.f

SUBROUTINE GasTranLayer2Space(raaGasAbs,raL2S,iLay)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kProfLayer)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iLay
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay      = current layer number
! raaGasAbs = matrix containing the raw abs coeffs
! raL2S     = array containing the layer-to-space results



INTEGER :: iI,iJ

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=iLay,kProfLayer
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

DO iI=1,kMaxPts
  raL2S(iI)=EXP(-raL2S(iI))
END DO

RETURN
END SUBROUTINE GasTranLayer2Space

!************************************************************************
! this subroutine does the layer to ground optical depths by doing
! the calculation between frequency indices from layer iLay-->1
! the result is stored in raL2S
! this is mainly used by the OpticalLayer2Space printout routines in rdprof.f

SUBROUTINE GasOptLayer2Ground(raaGasAbs,raL2S,iLay)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kProfLayer)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iLay
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay      = current layer number
! raaGasAbs = matrix containing the raw abs coeffs
! raL2S     = array containing the layer-to-space results



INTEGER :: iI,iJ

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=1,iLay
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

RETURN
END SUBROUTINE GasOptLayer2Ground

!************************************************************************
! this subroutine does the layer to ground transmittances for gas iG by doing
! the calculation between frequency indices from layer 1-->iLay
! the result is stored in raL2G
! this is mainly used by the OpticalLayer2Space printout routines in rdprof.f

SUBROUTINE GasTranLayer2Ground(raaGasAbs,raL2G,iLay)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kProfLayer)
REAL, INTENT(OUT)                        :: raL2G(kMaxPts)
INTEGER, INTENT(IN)                      :: iLay
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay      = current layer number
! raaGasAbs = matrix containing the raw abs coeffs
! raL2S     = array containing the layer-to-space results



INTEGER :: iI,iJ

DO iI=1,kMaxPts
  raL2G(iI)=0.0
END DO

DO iJ=1,iLay
  DO iI=1,kMaxPts
    raL2G(iI)=raL2G(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

DO iI=1,kMaxPts
  raL2G(iI)=EXP(-raL2G(iI))
END DO

RETURN
END SUBROUTINE GasTranLayer2Ground

!************************************************************************
! this subroutine does the MP to space transmission coeffs by doing
! the calculation between frequency indices from layer iLay-->iLM
! where iLm=nearest 100 index above iLay, or iNpmix
! the result is stored in raL2S

SUBROUTINE MP2SpaceTrans(raaGasAbs,raL2S,iIpmix,iNpmix)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iIpmix
INTEGER, INTENT(IN)                      :: iNpmix
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIpMix    = current mixed path number
! iNpMix    = total number of mixed paths
! raaGasAbs = matrix containing the raw mixed path abs coeffs
! raL2S     = array containing the layer-to-space results



! local variables
INTEGER :: iI,iJ,iMax,iFound,iUpper,idiv

! first find the index where the (space) "upper" limit is at
iFound=-1
iI=0
iJ=1
iMax=idiv(iNpmix,kProfLayer)
iUpper=MOD(iNpmix,kProfLayer)
IF (iUpper /= 0) THEN
  iMax=iMax+1
END IF

15   CONTINUE
IF ((iFound < 0) .AND. (iJ <= iMax)) THEN
  IF ((kProfLayer*iI <= iIpmix) .AND. (kProfLayer*iJ >= iIpmix)) THEN
    iFound=1
    iUpper = kProfLayer*iJ
    IF (iUpper > iNpmix) THEN
      iUpper=iNpmix
    END IF
  ELSE
    iI=iI+1
    iJ=iJ+1
  END IF
  GO TO 15
END IF

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=iIpmix,iUpper
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

DO iI=1,kMaxPts
  raL2S(iI)=EXP(-raL2S(iI))
END DO

RETURN
END SUBROUTINE MP2SpaceTrans

!************************************************************************
! this subroutine does the MP to space optical depths by doing
! the calculation between frequency indices from layer iLay-->iLM
! where iLm=nearest 100 index above iLay, or iNpmix
! the result is stored in raL2S

SUBROUTINE MP2SpaceOptDp(raaGasAbs,raL2S,iIpmix,iNpmix)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iIpmix
INTEGER, INTENT(IN)                      :: iNpmix
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIpMix    = current mixed path number
! iNpMix    = total number of mixed paths
! raaGasAbs = matrix containing the raw mixed path abs coeffs
! raL2S     = array containing the layer-to-space results



! local variables
INTEGER :: iI,iJ,iMax,iFound,iUpper,idiv

! first find the index where the (space) "upper" limit is at
iFound=-1
iI=0
iJ=1
iMax=idiv(iNpmix,kProfLayer)
iUpper=MOD(iNpmix,kProfLayer)
IF (iUpper /= 0) THEN
  iMax=iMax+1
END IF

15   CONTINUE
IF ((iFound < 0) .AND. (iJ <= iMax)) THEN
  IF ((kProfLayer*iI <= iIpmix) .AND. (kProfLayer*iJ >= iIpmix)) THEN
    iFound=1
    iUpper = kProfLayer*iJ
    IF (iUpper > iNpmix) THEN
      iUpper=iNpmix
    END IF
  ELSE
    iI=iI+1
    iJ=iJ+1
  END IF
  GO TO 15
END IF

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=iIpmix,iUpper
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

RETURN
END SUBROUTINE MP2SpaceOptDp

!************************************************************************
! this subroutine does the MP to space transmission coeffs by doing
! the calculation between frequency indices from layer iLay-->iLM
! where iLm=nearest 100 index above iLay, or iNpmix
! the result is stored in raL2S

SUBROUTINE MP2GroundTrans(raaGasAbs,raL2S,iIpmix,iNpmix)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iIpmix
INTEGER, INTENT(IN)                      :: iNpmix
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIpMix    = current mixed path number
! iNpMix    = total number of mixed paths
! raaGasAbs = matrix containing the raw mixed path abs coeffs
! raL2S     = array containing the layer-to-space results



! local variables
INTEGER :: iI,iJ,iMax,iFound,iUpper,idiv

! first find the index where the (space) "upper" limit is at
iFound=-1
iI=0
iJ=1
iMax=idiv(iNpmix,kProfLayer)
iUpper=MOD(iNpmix,kProfLayer)
IF (iUpper /= 0) THEN
  iMax=iMax+1
END IF

15   CONTINUE
IF ((iFound < 0) .AND. (iJ <= iMax)) THEN
  IF ((kProfLayer*iI <= iIpmix) .AND. (kProfLayer*iJ >= iIpmix)) THEN
    iFound=1
    iUpper = kProfLayer*iJ
    IF (iUpper > iNpmix) THEN
      iUpper=iNpmix
    END IF
  ELSE
    iI=iI+1
    iJ=iJ+1
  END IF
  GO TO 15
END IF

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=1,iIpmix
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

DO iI=1,kMaxPts
  raL2S(iI)=EXP(-raL2S(iI))
END DO

RETURN
END SUBROUTINE MP2GroundTrans

!************************************************************************
! this subroutine does the MP to space optical depths by doing
! the calculation between frequency indices from layer iLay-->iLM
! where iLm=nearest 100 index above iLay, or iNpmix
! the result is stored in raL2S

SUBROUTINE MP2GroundOptDp(raaGasAbs,raL2S,iIpmix,iNpmix)


REAL, INTENT(IN)                         :: raaGasAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raL2S(kMaxPts)
INTEGER, INTENT(IN)                      :: iIpmix
INTEGER, INTENT(IN)                      :: iNpmix
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iIpMix    = current mixed path number
! iNpMix    = total number of mixed paths
! raaGasAbs = matrix containing the raw mixed path abs coeffs
! raL2S     = array containing the layer-to-space results



! local variables
INTEGER :: iI,iJ,iMax,iFound,iUpper,idiv

! first find the index where the (space) "upper" limit is at
iFound=-1
iI=0
iJ=1
iMax=idiv(iNpmix,kProfLayer)
iUpper=MOD(iNpmix,kProfLayer)
IF (iUpper /= 0) THEN
  iMax=iMax+1
END IF

15   CONTINUE
IF ((iFound < 0) .AND. (iJ <= iMax)) THEN
  IF ((kProfLayer*iI <= iIpmix) .AND. (kProfLayer*iJ >= iIpmix)) THEN
    iFound=1
    iUpper = kProfLayer*iJ
    IF (iUpper > iNpmix) THEN
      iUpper=iNpmix
    END IF
  ELSE
    iI=iI+1
    iJ=iJ+1
  END IF
  GO TO 15
END IF

DO iI=1,kMaxPts
  raL2S(iI)=0.0
END DO

DO iJ=1,iIpmix
  DO iI=1,kMaxPts
    raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
  END DO
END DO

RETURN
END SUBROUTINE MP2GroundOptDp

!************************************************************************
! this subroutine appends _VT_iRegr to end of caOutName to get caVTName
! this is when we use regression to try to predict the best NLTE temps

SUBROUTINE VTName_regr(iRegr,caOutName,caVTFile)


INTEGER, INTENT(IN OUT)                  :: iRegr
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
CHARACTER (LEN=80), INTENT(OUT)          :: caVTFile
IMPLICIT NONE




CHARACTER (LEN=2) :: caString
INTEGER :: iInt,iInt1

DO iInt=1,80
  caVTFile(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caOutName(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caVTFile(1:iInt)=caOutName(1:iInt)
caVTFile(iInt+1:iInt+3)='_VT'

! now process iRegr so that we end up with a right padded string
! eg 2 ---> '2 ', 12 ---> '12' etc
WRITE(caString,15) iRegr
15   FORMAT(I2)
! this is right justified ... change to left justified
iInt=1
16   CONTINUE
IF (caString(iInt:iInt) == ' ') THEN
  iInt=iInt+1
  GO TO 16
END IF

iInt1 = 80
21   CONTINUE
IF ((caVTFile(iInt1:iInt1) == ' ') .AND. (iInt1 >= 1)) THEN
  iInt1 = iInt1-1
  GO TO 21
END IF
iInt1 = iInt1 + 1

caVTFile(iInt1:iInt1+(2-iInt))=caString(iInt:2)

RETURN
END SUBROUTINE VTName_regr

!************************************************************************
! this subroutine appends _RTP_irtp to end of caOutName to get caVTName

SUBROUTINE VTName_rtp(iRtp,caOutName,caVTFile)


INTEGER, INTENT(IN OUT)                  :: iRtp
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
CHARACTER (LEN=80), INTENT(OUT)          :: caVTFile
IMPLICIT NONE




CHARACTER (LEN=5) :: caString
INTEGER :: iInt,iInt1

DO iInt=1,80
  caVTFile(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caOutName(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caVTFile(1:iInt)=caOutName(1:iInt)
caVTFile(iInt+1:iInt+5)='_RTP_'

! now process iRtp so that we end up with a right padded string
! eg 2 ---> '2 ', 12 ---> '12' etc
WRITE(caString,15) iRtp
15   FORMAT(I5)
! this is right justified ... change to left justified
iInt=1
16   CONTINUE
IF (caString(iInt:iInt) == ' ') THEN
  iInt=iInt+1
  GO TO 16
END IF

iInt1 = 80
21   CONTINUE
IF ((caVTFile(iInt1:iInt1) == ' ') .AND. (iInt1 >= 1)) THEN
  iInt1 = iInt1-1
  GO TO 21
END IF
iInt1 = iInt1 + 1

caVTFile(iInt1:iInt1+(5-iInt))=caString(iInt:5)

RETURN
END SUBROUTINE VTName_rtp

!************************************************************************
! this subroutine appends _VT_iRTP to end of caOutName to get caVTName
! so iRegr is NOT used!
! this is when we use a polynom fit, or a predictor fit, to estimate NLTE temps
! this bloody thing DOES NOT WORK!!!!
! this bloody thing DOES NOT WORK!!!!
! this bloody thing DOES NOT WORK!!!!

SUBROUTINE VTName(iRTP,caOutName,caVTFile)


INTEGER, INTENT(IN OUT)                  :: iRTP
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
CHARACTER (LEN=80), INTENT(OUT)          :: caVTFile
IMPLICIT NONE

! input params


! output params


CHARACTER (LEN=5) :: caString
INTEGER :: iInt,iInt1,iI,iJ

DO iInt=1,80
  caVTFile(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caOutName(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caVTFile(1:iInt)        = caOutName(1:iInt)
caVTFile(iInt+1:iInt+5) = '_RTP_'

! now process iRTP so that we end up with a right padded string
! eg 2 ---> '2 ', 12 ---> '12' etc
! assume at most there are 99999 profiles in one RTP file
caString = '     '
WRITE(caString,15) iRTP
15   FORMAT(I5)
! this is right justified ... change to left justified
iInt=1
16   CONTINUE
IF (caString(iInt:iInt) == ' ') THEN
  iInt=iInt+1
  GO TO 16
END IF

iInt1 = 80
21   CONTINUE
IF ((caVTFile(iInt1:iInt1) == ' ') .AND. (iInt1 >= 1)) THEN
  iInt1 = iInt1-1
  GO TO 21
END IF
iInt1 = iInt1 + 1

!      caVTFile(iInt1:iInt1+(5-iInt)) = caString(iInt:5)
iJ = 0
DO iI = iInt,5
  caVTFile(iInt1+iJ:iInt1+iJ) = caString(iI:iI)
  iJ = iJ + 1
END DO

RETURN
END SUBROUTINE VTName

!************************************************************************
! this subroutine appends _COL to end of caOutName to get "new" caJacobFileName

SUBROUTINE jacob2Name(caJacobFile2,caJacobFile)


NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

CHARACTER (LEN=80) :: caJacobFile,caJacobFile2

INTEGER :: iInt

DO iInt=1,80
  caJacobFile2(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caJacobFile(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caJacobFile2(1:iInt) = caJacobFile(1:iInt)
caJacobFile2(iInt+1:iInt+4) = '_COL'

RETURN
END SUBROUTINE jacob2Name

!************************************************************************
! this subroutine appends _FLUX to end of caOutName to get caFluxName

SUBROUTINE FluxName(caFluxFile,caOutName)


CHARACTER (LEN=80), INTENT(OUT)          :: caFluxFile
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'



INTEGER :: iInt

DO iInt=1,80
  caFluxFile(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caOutName(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caFluxFile(1:iInt)=caOutName(1:iInt)

! before Jun 2013
!      IF (kFLux .LE. 2) THEN
!        caFluxFile(iInt+1:iInt+5)='_FLUX'
!      ELSEIF (kFLux .LE. 5) THEN
!        caFluxFile(iInt+1:iInt+4)='_OLR'
!      ELSE
!        write(kStdErr,*) 'Unknown option for kFlux = ',kFlux
!        Call DoStop
!      END IF

! after Jun 2013
WRITE(kStdWarn,*) 'kFlux = ',kFlux
IF (kFLux == 1) THEN
  caFluxFile(iInt+1:iInt+5)='_DOWN'   !! down welling flux at bottom of 100 layers
ELSE IF (kFLux == 2) THEN
  caFluxFile(iInt+1:iInt+5)='_HEAT'   !! up-down heating rate at all 100 layers
ELSE IF (kFLux == 3) THEN
  caFluxFile(iInt+1:iInt+3)='_UP'     !! up welling flux at to[ of 100 layers
ELSE IF (kFLux == 4) THEN
  caFluxFile(iInt+1:iInt+4)='_OLR'    !! upwelling flux at TOA
ELSE IF (kFLux == 5) THEN
  caFluxFile(iInt+1:iInt+5)='_OLR3'   !! upwelling flux at TOA, tropopause,dnwell at gnd
ELSE IF (kFLux == 6) THEN
  caFluxFile(iInt+1:iInt+4)='_ALL'    !! up,down welling flux at all tops/bottoms of each layer
ELSE
  WRITE(kStdErr,*) 'Unknown option for kFlux = ',kFlux
  CALL DoStop
END IF

RETURN
END SUBROUTINE FluxName

!************************************************************************
! this subroutine appends _PLANCK to end of caOutName to get caPlanckName

SUBROUTINE PlanckName(caPlanckFile,caOutName)


NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
IMPLICIT NONE

CHARACTER (LEN=80) :: caPlanckFile

INTEGER :: iInt

DO iInt=1,80
  caPlanckFile(iInt:iInt)=' '
END DO

iInt=80
11   CONTINUE
IF ((caOutName(iInt:iInt) == ' ') .AND. (iInt >= 1)) THEN
  iInt=iInt-1
  GO TO 11
END IF

caPlanckFile(1:iInt)=caOutName(1:iInt)
caPlanckFile(iInt+1:iInt+7)='_PLANCK'

RETURN
END SUBROUTINE PlanckName

!************************************************************************
! this subroutine writes out the header for the summary text file "header.head"
! this subroutine writes out the header for the binary output file ...
! any other writes to this file will be appended to the end

! as the header is written out, the subroutine checks that print option 1 has
! been specified not more than once, print option 3 specified not more than
! once and print option 7 not more than once for the ith atmosphere

! kLongOrShort enables shorted output files to be written out (-1,+1),
! or just the basic results file (0)

SUBROUTINE PrepareOutput(caDriver,caOutName,caJacobFile,caJacobFile2,  &
    caFluxFile,caPlanckFile,iOutFileName,iNumNLTEGases,  &
    rFrLow,rFrHigh,iFileIDLo,iFileIDHi,caComment,  &
    iNumGases,iaGases,raaAmt,raaTemp,raaPress,raaPartPress,  &
    raaRAmt,       raaRPartPress, raPressLevels,iProfileLayers,  &
    iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVT,  &
    iNatm,iNatm2,iaNumLayers,iaaRadLayer,  &
    raTSpace,raTSurf,raSatAngle,raSatHeight, raaaSetEmissivity,iaSetEms,  &
    iOutTypes,iaPrinter,iaAtmPr,iaNp,iaaOp,raaUserPress, iJacob,iaJacob,  &
    iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,  &
    rakThermalAngle,iakThermalJacob,iaOutNumbers, iTotal,iTotalStuff,  &
    iDoUpperAtmNLTE,iDumpAllUASpectra,iDumpAllUARads)


CHARACTER (LEN=80), INTENT(IN OUT)       :: caDriver
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
CHARACTER (LEN=80), INTENT(OUT)          :: caFluxFile
NO TYPE, INTENT(IN OUT)                  :: caPlanckFi
NO TYPE, INTENT(IN OUT)                  :: iOutFileNa
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
REAL, INTENT(IN)                         :: rFrLow
REAL, INTENT(IN)                         :: rFrHigh
INTEGER, INTENT(IN)                      :: iFileIDLo
INTEGER, INTENT(IN)                      :: iFileIDHi
CHARACTER (LEN=120), INTENT(IN OUT)      :: caComment
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
REAL, INTENT(IN)                         :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
REAL, INTENT(IN OUT)                     :: raaRAmt(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaRPartPr
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
INTEGER, INTENT(IN)                      :: iNpmix
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: caaMixFile
NO TYPE, INTENT(IN OUT)                  :: iMixFileLi
REAL, INTENT(IN OUT)                     :: raMixVT(kMixFilRows)
INTEGER, INTENT(IN)                      :: iNatm
INTEGER, INTENT(IN OUT)                  :: iNatm2
NO TYPE, INTENT(IN OUT)                  :: iaNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
INTEGER, INTENT(IN)                      :: iaSetEms(kMaxAtm)
INTEGER, INTENT(IN)                      :: iOutTypes
INTEGER, INTENT(IN)                      :: iaPrinter(kMaxPrint)
INTEGER, INTENT(IN OUT)                  :: iaAtmPr(kMaxPrint)
INTEGER, INTENT(IN)                      :: iaNp(kMaxPrint)
INTEGER, INTENT(IN OUT)                  :: iaaOp(kMaxPrint,kPathsOut)
NO TYPE, INTENT(IN OUT)                  :: raaUserPre
INTEGER, INTENT(IN)                      :: iJacob
INTEGER, INTENT(IN)                      :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iakSolar(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
INTEGER, INTENT(IN OUT)                  :: iakThermal(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakThermal
NO TYPE, INTENT(IN OUT)                  :: iakThermal
NO TYPE, INTENT(IN OUT)                  :: iaOutNumbe
INTEGER, INTENT(IN OUT)                  :: iTotal
NO TYPE, INTENT(IN OUT)                  :: iTotalStuf
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raPressLevels, iProfileLayers = actual number of pressure layers
! iOutFileName = does caOutName exist, or is stuff dumped to screen? for fluxes
! caComment   = comment that the user put into *RADNCE
! iJacob      = number of d/dq we will output (>=1)
! iaJacob     = list of GasID's whose  d/dq we will output
! caDriver     = root name of input file (not really needed here)
! caOutName   = name of output binary file, specified in *OUTPUT
! caDummyFile = name of output binary file, to save Dummy stuff
! rFrLow      = lower wavenumber bound
! rFrHigh     = upper wavenumber bound
! iFileIDLo      = from rFrLow, lower file index (1-kNumkFile)
! iFileIDHi      = from rFrHigh, upper file index (1-kNumkFile)
! iNumGases   = number of gases read in from *GASFIL + *XSCFIL
! iaGases     = integer array containig the GAS ID's
! raaAmt      = array containing amount profiles for each gas
! raaTemp     = array containing temp profiles for each gas
! raaPartPress= array containing partial pressure profiles for each gas
! iNpMix      = number of mixed paths
! raaMix      = complete mixing table
! caMixFName  = name of file containing character MP info (from *MIXFIL)
! iMixFileLines = number of lines containig character info for mixtable
! iNatm       = number of atmospheres read in from *RADIL
! iNatm2      = number of atmospheres read in from *OUTPUT
! iaNumLayers = number of layers in each atmosphere
! iaaRadLayer = for each atmosphere, a list of the layers making it up
! raTSpace    = for each atmosphere, the atmosphere backgnd temperature
! raTSurf     = for each atmosphere, the surface temperature
! raSatAngle  = for each atmosphere, the satellite view angle
! raSatHeight = for each atmosphere, the satellite height
! iOutTypes   = number of output options found in *OUTPUT
! iaPrinter   = for each output option, what to be printed (1,2 or 3)
! iaAtmPr     = number of layers to be printed for each atmosphere
! iaNp        = number of paths/MP/layers to be  output for each print option
! iaaOp       = list of paths/MP/layers to be output for each print option
! raaUserPress= list of output pressures for each print option
! iaSetEms    = number of emissivity data points per atmosphere
! raaaSetEms  = emissivity data
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! rakSolarRefl   =solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaOutNumbers = how many of each print option there are to output
! iTotal = how many of the kcomp files are gonna be unchunked
! iDumpAllUARads = do we dump rads for all layers (-1) or a specific number?
REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
REAL :: rakSolarRefl(kMaxAtm)
INTEGER :: iaOutNumbers(kMaxPrint),iOutFileName
INTEGER :: iakThermalJacob(kMaxAtm)


CHARACTER (LEN=130) :: caaMixFileLines(kProfLayer)

INTEGER :: iMixFileLines
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raaUserPress(kMaxPrint,kProfLayer)

REAL :: raSatHeight(kMaxAtm),raPressLevels(kProfLayer+1)

REAL :: raaPartPress(kProfLayer,kGasStore)
REAL :: raaRPartPress(kProfLayer,kGasStore)

INTEGER :: iaNumLayers(kMaxAtm)
REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
INTEGER :: iProfileLayers
INTEGER :: iDoUpperAtmNLTE,iNumNLTEGases,iDumpAllUARads,iDumpAllUASpectra

INTEGER :: iIOUN,iIOUN1,iIOUN2,iIOUN_JAC2,iI,iJ,iK,iFileErr,iEnd,iP,iOk
INTEGER :: iOutputOptionNum,iNumLay,DoGasJacob
CHARACTER (LEN=80) :: caJacobFile,caJacobFile2, caPlanckFile

INTEGER :: CheckDoubleEntry,iImportant
INTEGER :: iNatmJac,iaLayerJac(kMaxAtm),iIOUN_Flux,iIOUN_Planck,iIOUN_Cloud
REAL :: raParams(kMaxUserSet),raPActualAvg(kProfLayer),rP
CHARACTER (LEN=4) :: caStrJunk(7)
CHARACTER (LEN=80) :: caFCloudName
REAL :: raSumTotalGasAmt(kMaxGas)

!this is for kLongOrShort = 0
INTEGER :: iTag,iTotalStuff
INTEGER :: iCo2,iaCO2path(kProfLayer),iaLayerFlux(kMaxAtm)
INTEGER :: iaJunkFlux(2*kProfLayer)

CHARACTER (LEN=120) :: caStr1,caStr2,caStr3

INCLUDE '../INCLUDE/gasIDnameparam.f90'

DO iI = 1, kMaxGas
  raSumTotalGasAmt(iI) = 0.0
END DO

iDumpAllUASpectra = -1 !!assume no need to dump out CO2 UA spectra
IF (iDoUpperAtmNLTE > 0) THEN
!!need to see if we are dumping out CO2 paths
  iCO2 = -1
  DO iTag = 1,iNumGases
    IF (iaGases(iTag) == 2) THEN
      iCO2 = iTag
    END IF
  END DO
  IF (iCO2 > 0) THEN
    DO iTag = 1,kProfLayer
!!! these are the CO2 paths
      iaCO2path(iTag) = iTag + (iCO2-1)*kProfLayer
    END DO
!!check to see if any of them are being output
    DO iI = 1,iOutTypes
      IF (iaPrinter(iI) == 1) THEN !!! these are paths to be output
        iEnd = iaNp(iI)
        DO iTag = 1,kProfLayer
          DO iP = 1,iEnd
            IF (iaaOp(iI,iP) == iaCO2path(iTag)) THEN
              iDumpAllUASpectra = +1
            END IF
          END DO
        END DO
      END IF
    END DO
  ELSE
    WRITE(kStdErr,*) 'oh oh, need to do NLTE and not using CO2????'
    CALL DoStop
  END IF
END IF

iDumpAllUARads = 0     !!assume dump out rads only for specific levels
!!if it remains at 0 thru the subroutine, then it
!!is set to -1 at end of routine
!!if it ever increases beyond +1, kCARTA barfs

!compute avg layer pressure, to be output
DO iI = 1,kProfLayer
  raPActualAvg(iI) = 0.0
END DO
iOK = kProfLayer - iProfileLayers + 1
DO iI = iOK,kProfLayer
  rP = raPressLevels(iI) - raPresslevels(iI+1)
  raPActualAvg(iI) = rP/LOG(raPressLevels(iI)/raPresslevels(iI+1))
END DO

WRITE(kStdWarn,*)'Preparing header info for the output files ..'

DO iImportant=1,kMaxPrint
  iaOutNumbers(iImportant)=0
END DO

800  FORMAT(A1)
801  FORMAT(A80)

iIOUN =  kStdWarn
iIOUN1 = kStdkCarta

iTotalStuff = 0 !total num of spectra,mixed spectra, radiances to output

! set up raParams; recall all are integers, so change to real!!!
raParams(1)  = kLayer2Sp    * 1.0
raParams(2)  = kCKD         * 1.0
raParams(3)  = kGasTemp     * 1.0
raParams(4)  = kLongOrShort * 1.0
raParams(5)  = kJacobOutput * 1.0
raParams(6)  = kFlux        * 1.0
raParams(7)  = kSurfTemp    * 1.0
raParams(8)  = kTempJac     * 1.0
raParams(9)  = kRTP         * 1.0
raParams(10) = kActualJacs  * 1.0
raParams(11) = -98765.0  !!dummy variables
raParams(12) = -98765.0  !!dummy variables

WRITE(kStdWarn,*) 'kCKD = ',kCKD

! check to make sure that the user first wanted to specify paths, then mixed
! paths, then radiances; else STOP ie we have to find paPrinter(iI)=1,2,3 or
! in that order eg 1,1,1,1 or 2,3,3,3,3  but not 3,1,1 etc
iP=-1
DO iI=1,iOutTypes
  IF (iaPrinter(iI) >= iP) THEN
    iP = iaPrinter(iI) !ok; ia(printer(iI) is going UP or staying same
  ELSE !oh oh iaPrinter(iI) is coming down!!! readers cannot handle this
    WRITE(kStdErr,*)'in *OUTPUT need to specify paths, then mixed'
    WRITE(kStdErr,*)'   paths, and then radiances ie order is '
    WRITE(kStdErr,*)'   important. Please rewrite *OUTPUT section'
    CALL DoSTOP
  END IF
END DO

! if kLongOrShort = 1, output all info, and summarize in kStdWarn
! if kLongOrShort = 0, output ONLY DATA, but summarize in kStdWarn
!                      we need to figure out iTotalStuf ==> do this loop
! if kLongOrShort = -1, output shortened version of binary file

! first open summary text file as a fresh file to be written to
IF (kLongOrShort >= 0) THEN
  
  WRITE(iIOUN,*) 'SUMMARY OF INPUT DRIVER FILE ... '
  WRITE(iIOUN,*) '     '
  WRITE(iIOUN,*) '     '
  
! first output general info -----------------------------------------
  WRITE(iIOUN,*) 'GENERAL INFORMATION : '
  WRITE(iIOUN,*) caVersion
  WRITE(iIOUN,*) 'Number of layers = ',kProfLayer
  WRITE(iIOUN,*) 'Number of parameters in *PARAMS = ',kMaxUserSet
  WRITE(iIOUN,*) (raParams(iI),iI=1,kMaxUserSet)
  WRITE(iIOUN,*) caComment
  WRITE(iIOUN,*) 'Freq endpts = ',rFrLow,rFrHigh
  WRITE(iIOUN,*) 'k-comp file endpts = ',iFileIDLo,iFileIDHi
  WRITE(iIOUN,*) 'Long/short output file = ',kLongOrShort
! v1.04 to v1.08 had these next 2 lines; now remove them
!        WRITE(iIOUN,*) 'M1000mb,M100mb,MSub = ',M1000mb,M100mb,MSubLayer
!        WRITE(iIOUN,*) 'M50mb,M10mb,MThick = ',M50mb,M10mb,MThickLayer
  WRITE(iIOUN,*) 'average layer pressures ...'
  WRITE(iIOUN,*) (raPActualAvg(iI),iI=1,kProfLayer)
  WRITE(iIOUN,*) (iaaOverrideDefault(1,iI),iI=1,10)    !! general settings
  WRITE(iIOUN,*) (iaaOverrideDefault(2,iI),iI=1,10)    !! radtrans settings
  WRITE(iIOUN,*) (iaaOverrideDefault(3,iI),iI=1,10)    !! iLBLRTM settings
  WRITE(iIOUN,*) '***********************************************'
  
! then output path ID stuff ------------------------------------------
  WRITE(iIOUN,*) 'PATH ID INFO'
  WRITE(iIOUN,*) 'iNumPaths = ',iNumGases*kProfLayer
  WRITE(iIOUN,*) 'Num Layers in Profile = ',iProfileLayers
  WRITE(iIOUN,*) 'So start showing info from layer ',kProfLayer-iProfileLayers+1
  
  caStr1 = '  Path# GID  Press      PartP        PPMV        Temp         Amnt    ||      RefPP     RefAmt  Amt/RAmt'
  caStr2 = '             (atm)      (atm)                     (K)     (kmole/cm2) ||       (atm)  (kmol/cm2)        '
  caStr2 = '             (mb )      (mb )                     (K)  (molecule/cm2) ||       (mb ) (molecule/cm2)          '
  caStr3 = '----------------------------------------------------------------------||--------------------------------'
  DO iI=1,iNumGases
    WRITE(iIOUN,234) iI,iaGases(iI),caGID(iaGases(iI))
    WRITE(iIOUN,7169) caStr1
    WRITE(iIOUN,7169) caStr2
    WRITE(iIOUN,7169) caStr3
    DO iJ=kProfLayer-iProfileLayers+1,kProfLayer
      iP=(iI-1)*kProfLayer+iJ
      raSumTotalGasAmt(iaGases(iI)) = raSumTotalGasAmt(iaGases(iI)) + raaAmt(iJ,iI)
! this is writing pressures in mb
      WRITE(iIOUN,7170)iP,iaGases(iI),raPActualAvg(iJ),raaPartPress(iJ,iI)*kAtm2mb,  &
          raaPartPress(iJ,iI)*kAtm2mb/raPActualAvg(iJ)*1.0E6,raaTemp(iJ,iI),raaAmt(iJ,iI)*kAvog,'||',  &
          raaRPartPress(iJ,iI)*kAtm2mb,raaRAmt(iJ,iI)*kAvog,raaAmt(iJ,iI)/raaRAmt(iJ,iI)
! this is writing pressures in atm
!            WRITE(iIOUN,7171)iP,iaGases(iI),raPActualAvg(iJ)/kAtm2mb,raaPartPress(iJ,iI),
!     $                   raaPartPress(iJ,iI)/(raPActualAvg(iJ)/kAtm2mb)*1.0e6,raaTemp(iJ,iI),raaAmt(iJ,iI),'||',
!     $                   raaRPartPress(iJ,iI),raaRAmt(iJ,iI),raaAmt(iJ,iI)/raaRAmt(iJ,iI)
    END DO
    WRITE(kStdWarn,*) '++++++++++++++++++++++++++++++++'
  END DO
  234    FORMAT ('index = ',I3,' gas HITRAN ID = ',I3,' molecule = ',A20)
  
! write sum
  WRITE(iIOUN,*) ' iI  iGasID  Name     total molecules/cm2'
  WRITE(iIOUN,*) '--------------------------------------------------'
  DO iI = 1,iNumGases
    WRITE(iIOUN,7172) iI,iaGases(iI),caGID(iaGases(iI)),raSumTotalGasAmt(iaGases(iI))*kAvog
  END DO
  WRITE(iIOUN,*) '--------------------------------------------------'
  
! then output list of paths to be output
  WRITE(iIOUN,*) 'list of paths to be output ...'
  iP=0
  DO iI=1,iOutTypes
    IF (iaPrinter(iI) == 1) THEN
      iP=iP+1
      iOutputOptionNum=iI
      IF (iaNp(iI) > iNumGases*kProfLayer) THEN
        WRITE(kStdErr,*)'Cannot have more paths to be output than '
        WRITE(kStdErr,*)'iNumGases*kProfLayer'
        WRITE(kStdErr,*)iI,iaNp(iI)
        CALL DoSTOP
      END IF
      IF (iaNp(iI) > 0) THEN
        iEnd=iaNp(iI)
        iTotalStuff = iTotalStuff + iEnd
        WRITE(iIOUN,*) 'Number of Paths to be output=',iEnd
        WRITE(iIOUN,*) (iaaOp(iI,iJ),iJ=1,iEnd)
      ELSE
        iEnd = kProfLayer*iNumGases
        iTotalStuff = iTotalStuff + iEnd
        WRITE(iIOUN,*) 'Number of Paths to be output=',iEnd
        WRITE(iIOUN,*) (iaaOp(iI,iJ),iJ=1,iEnd)
      END IF
      iaOutNumbers(iI)=iEnd
      iNumLay=iEnd
    END IF
  END DO
  IF (iP > 1) THEN
    WRITE(kStdErr,*)'Found more than 1 set of options in *OUTPUT where'
    WRITE(kStdErr,*)'iDat set to 1 (outputting single paths)'
    WRITE(kStdErr,*)'... exiting program!'
    CALL DoSTOP
  END IF
  IF (iP == 0) THEN
! inform user no single set of paths to be output
    WRITE(iIOUN,*) iP
  END IF
  IF (iP == 1) THEN
    iOK=CheckDoubleEntry(iaaOp,iOutputOptionNum,iNumLay)
    IF (iOK < 0) THEN
      WRITE(kStdErr,*)'On checking list of paths to output, have found '
      WRITE(kStdErr,*)'some to be output more than once. messing up'
      WRITE(kStdErr,*)'everything in readatmos.m .. please try again'
      CALL DoSTOP
    END IF
  END IF
  
! then output mixed paths -----------------------------------------------
  WRITE(iIOUN,*) '***********************************************'
  WRITE(iIOUN,*) 'MIXED PATHS'
  WRITE(iIOUN,*) 'Number of mixed paths = ',iNpmix
  WRITE(iIOUN,*) 'Number of uncommented lines in mixtable=', iMixFileLines
  IF (iNpMix > 0) THEN
    DO iI=1,iMixFileLines
      WRITE(iIOUN,*) caaMixFileLines(iI)
    END DO
    
! then output list of mixed paths temperatures
    WRITE(iIOUN,*) 'mixed path temperatures ....'
    WRITE(iIOUN,*) (raMixVT(iI),iI=1,iNpmix)
! then output list of mixed paths to be printed
    WRITE(iIOUN,*) 'list of mixed paths to be output ...'
    iP=0
    DO iI=1,iOutTypes
      IF (iaPrinter(iI) == 2) THEN
        iP=iP+1
        iOutputOptionNum=iI
        IF (iaNp(iI) > iNpmix) THEN
          WRITE(kStdErr,*)'Error! Cannot print out more mixed paths'
          WRITE(kStdErr,*)'than iIpmix'
          CALL DoSTOP
        END IF
        IF (iaNp(iI) > 0) THEN
          iEnd=iaNp(iI)
          iTotalStuff = iTotalStuff + iEnd
          WRITE(iIOUN,*) '# mixed paths to be output=',iEnd
          WRITE(iIOUN,*) (iaaOp(iI,iJ),iJ=1,iEnd)
        ELSE
          iEnd=iNpmix
          iTotalStuff = iTotalStuff + iEnd
          WRITE(iIOUN,*) '# mixed paths to be output=',iEnd
          WRITE(iIOUN,*) (iaaOp(iI,iJ),iJ=1,iEnd)
        END IF
        iNumLay=iEnd
        iaOutNumbers(iI)=iEnd
      END IF
    END DO
    IF (iP > 1) THEN
      WRITE(kStdErr,*)'Found more than 1 options set in *OUTPUT where'
      WRITE(kStdErr,*)'iDat set to 2 (outputting mixed paths)'
      WRITE(kStdErr,*)'... exiting program!'
      CALL DoSTOP
    END IF
    IF (iP == 0) THEN
! inform user no mixed set of paths to be output
      WRITE(iIOUN,*) iP
    END IF
    IF (iP == 1) THEN
      iOK=CheckDoubleEntry(iaaOp,iOutputOptionNum,iNumLay)
      IF (iOK < 0) THEN
        WRITE(kStdErr,*)'On checking list of MIXED paths to output,found'
        WRITE(kStdErr,*)'some that will be output more than once'
        WRITE(kStdErr,*)'Please check *OUTPUT'
        CALL DoSTOP
      END IF
    END IF
    
  END IF
  
! now print the atmosphere information -------------------------------------
  WRITE(iIOUN,*) '***********************************************'
  WRITE(iIOUN,*) 'ATMOSPHERES'
  WRITE(iIOUN,*) 'Numb of atmospheres read from *RADFIL = ',iNatm
  WRITE(iIOUN,*) 'Max numb of emiss. pts = ',kEmsRegions
  IF (iNatm < iNatm2) THEN
    WRITE(kStdErr,*)'RADFIL had information for ',iNatm,' atmospheres'
    WRITE(kStdErr,*)'OUTPUT says there is info for ',iNatm2,' atms!'
    WRITE(kStdErr,*)'please recheck and retry!!!'
    CALL DoSTOP
  END IF
  
  DO iI=1,iNatm
    IF (iaNumLayers(iI) > iNpmix) THEN
      WRITE(kStdErr,*)'atm#',iI,' has too many layers!! (> ',iNpmix,')'
      CALL DoSTOP
    END IF
    IF (iaNumLayers(iI) < 0) THEN
      WRITE(kStdErr,*)'atm#',iI,' has 0 or less layers'
      CALL DoSTOP
    END IF
    WRITE(iIOUN,*) 'atm#',iI,' has',iaNumLayers(iI),' MP layers :'
    WRITE(iIOUN,*) (iaaRadLayer(iI,iJ),iJ=1,iaNumLayers(iI))
    WRITE(iIOUN,*) 'TSpace,TSurf,SatAngle,SatHeight = ',  &
        raTSpace(iI),raTSurf(iI),raSatAngle(iI),raSatHeight(iI)
    WRITE(iIOUN,*) 'kSolar,kSolarAngle,kSolarRefl = ',  &
        iakSolar(iI),rakSolarAngle(iI),rakSolarRefl(iI)
    WRITE(iIOUN,*) 'kThermal,kThermalAngle,kThermalJacob = ',  &
        iakThermal(iI),rakThermalAngle(iI),iakThermalJacob(iI)
! output the emissivities for this atmosphere, first freq pt then ems
    WRITE(iIOUN,*) 'The surface freqs/emissivities are : '
    WRITE(iIOUN,*) iaSetEms(iI)
    DO iJ=1,iaSetEms(iI)
      WRITE(iIOUN,*)raaaSetEmissivity(iI,iJ,1), raaaSetEmissivity(iI,iJ,2)
    END DO
! now output the list of radiances to be printed, for this atmosphere
    iP=0
    DO iJ=1,iOutTypes
      IF ((iaPrinter(iJ) == 3).AND.(iaAtmPr(iJ) == iI)) THEN
        iP=iP+1
        iOutputOptionNum=iJ
!              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .GT. 1)) THEN
!                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
!                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
!                CALL DoSTOP
!              END IF
!              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .LE. -1)) THEN
!                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
!                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
!                CALL DoSTOP
!              END IF
        IF (iaNp(iJ) <= -1) THEN  !!!dumping out rads for ALL layers
!                IF ((iNumNLTEGases .GT. 0) .AND. (iDoUpperAtmNLTE .GT. 0))THEN
!                  iDumpAllUARads = iDumpAllUARads + 1
!                  write(kStdWarn,*) iJ,iaNp(iJ)
!                  write(kStdWarn,*) 'looks like we will dump US rads'
!                END IF
          IF ((iDumpAllUARads > 1) .AND.(iDoUpperAtmNLTE <= 0))THEN
            WRITE(kStdErr,*) 'this is too complicated!!!'
            WRITE(kStdErr,*) 'NLTE code assumes ONE radiating atm'
            CALL DoStop
          END IF
        END IF
        IF (iaNp(iJ) > iaNumLayers(iI)) THEN
          WRITE(kStdErr,*)'Atm #',iI,'too many radiances to be printed!'
          CALL DoSTOP
        END IF
        IF (iaNp(iJ) == 0) THEN
          WRITE(kStdErr,*) 'Atm #',iI,' has 0 radiances to be printed!'
          CALL DoSTOP
        END IF
        IF (iaNp(iJ) > 0) THEN
          iEnd=iaNp(iJ)
          iaOutNumbers(iJ)=iEnd
          iTotalStuff = iTotalStuff + iEnd
          WRITE(iIOUN,*) '# of radiances to be printed=',iEnd
          WRITE(iIOUN,*) (iaaOp(iJ,iK),iK=1,iEnd)
          WRITE(iIOUN,*) (raaUserPress(iJ,iK),iK=1,iEnd)
          iOk=1
          DO iK=1,iEnd
            IF ((iaaOp(iJ,iK) > iaNumLayers(iI)) .OR.  &
                  (iaaOp(iJ,iK) <= 0))  THEN
              iOk=-1
            END IF
          END DO
          IF (iOk < 0) THEN
            WRITE(kStdErr,*)'Atm# ',iI,' has invalid radiance to output!'
            CALL DoSTOP
          END IF
        ELSE
          iEnd=iaNumLayers(iI)
          iaOutNumbers(iJ)=iEnd
          iTotalStuff = iTotalStuff + iEnd
          WRITE(iIOUN,*) '# of radiances to be printed=',iEnd
!note how we are flipping between iI and iJ here
!and instead of outting iaaOp, we are outputting iaaRadLayer
          WRITE(iIOUN,*) (iaaRadLayer(iI,iK),iK=1,iEnd)
          WRITE(iIOUN,*) (raaUserPress(iJ,iK),iK=1,iEnd)
        END IF
        iNumLay=iEnd
      END IF
    END DO
    IF (iP == 0) THEN
! inform user no mixed paths to be output for this atmosphere
      WRITE(iIOUN,*) iP
    END IF
    IF (iP > 1) THEN
      WRITE(kStdErr,*)'Have found more than 1 set of options in *OUTPUT '
      WRITE(kStdErr,*)'where iDat set to 3 (outputting radiances)'
      WRITE(kStdErr,*)'for atmosphere #',iI,'... exiting program!'
      CALL DoSTOP
    END IF
    
  END DO
!cccc        CLOSE(iIOUN)
WRITE(iIOUN,*) 'END SUMMARY OF INPUT DRIVER FILE ... '
WRITE(iIOUN,*) '     '
WRITE(iIOUN,*) '     '

! end if kLongOrShort > 0
END IF
!--------------- OUTPUT BINARY FILE --------------------------------

iTag = -1
DO iI = 1,kW
  IF ((rFrLow >= kaMinFr(iI)) .AND. (rFrHigh <= kaMaxFr(iI))) THEN
    iTag = iI
  END IF
END DO
IF ((iTag < 1) .OR. (iTag > kW)) THEN
  WRITE (kStdErr,*) 'Could not find kaTag for FreqStart,FreqStop!!'
  CALL DoStop
END IF

IF ((kRTP >= 0) .AND. (kWhichScatterCode == 5)) THEN
!! write out cloud info, straight from RTP file and after manipulation
  DO iI = 1,80
    caFCloudName(iI:iI) = ' '
  END DO
  caFCloudName = caOutName
  iJ = 80
  DO WHILE ((caFCloudName(iJ:iJ) == ' ') .AND. (iJ. GT. 0))
    iJ = iJ -1
  END DO
  iJ = iJ + 1
  caFCloudName(iJ:iJ+3) = '_CLD'
  
  caStrJunk(1) = 'typ '
  caStrJunk(2) = 'ctop'
  caStrJunk(3) = 'cbot'
  caStrJunk(4) = 'cng '
  caStrJunk(5) = 'csz '
  caStrJunk(6) = 'frac'
  caStrJunk(7) = 'fr12'
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*)'  after expand_scatter'
  WRITE(kStdWarn,*)'    Cloud1  (before/after)        Cloud2 (before/after)'
  WRITE(kStdWarn,*)'-------------------------------------------------------'
  DO iI=1,7
    WRITE(kStdWarn,*) caStrJunk(iI),raaRTPCloudParams0(1,iI),raaRTPCloudParamsF(1,iI),'<>',  &
        raaRTPCloudParams0(2,iI),raaRTPCloudParamsF(2,iI)
  END DO
  
  OPEN(UNIT=iIOUN_Cloud,FILE=caFCloudName,FORM='FORMATTED',STATUS='NEW',  &
      IOSTAT=iFileErr)
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,103) iFileErr, caFCloudName
    WRITE(kStdErr,*)'make sure the cloud dump file does not exist!'
    CALL DoSTOP
  END IF
  kTempUnitOpen = 1
  IF (k100layerCloud == +1) THEN
    WRITE(iIOUN_Cloud,*) '% 100 layer cloud, but dumping out what was in rtp file for TWO SLAB CLOUD'
  ELSE
    WRITE(iIOUN_Cloud,*) '% TWO SLAB CLOUD'
  END IF
  WRITE(iIOUN_Cloud,*) '% rows 1-7 are ctype(101/201/301=W/I/A),cprtop(mb),cprbot(mb)'
  WRITE(iIOUN_Cloud,*) '% cngwat(g/m2),cpsize(um),cfrac and cfrac12'
  WRITE(iIOUN_Cloud,*) '% cols 1-4 are CLOUD 1 (old/new) and CLOUD 2 (old/new)'
  DO iI = 1, 7
    WRITE(iIOUN_Cloud,*) '% ',caStrJunk(iI),raaRTPCloudParams0(1,iI),raaRTPCloudParamsF(1,iI),raaRTPCloudParams0(2,iI),  &
        raaRTPCloudParamsF(2,iI)
  END DO
  DO iI = 1, 7
    WRITE(iIOUN_Cloud,*) raaRTPCloudParams0(1,iI),raaRTPCloudParamsF(1,iI),raaRTPCloudParams0(2,iI),raaRTPCloudParamsF(2,iI)
  END DO
  CLOSE(iIOUN_Cloud)
  kTempUnitOpen = -1
END IF

! next open unformatted file as a fresh file to be written to
IF (iIOUN1 /= 6) THEN
  OPEN(UNIT=iIOUN1,FILE=caOutName,STATUS='NEW',  &
      FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,103) iFileErr, caOutName
    WRITE(kStdErr,*)'make sure the file does not exist!'
    CALL DoSTOP
  END IF
END IF
103  FORMAT('ERROR! number ',I5,' opening binary data file : ',/,A80)

kStdkCartaOpen=1
WRITE(kStdWarn,*) 'Opened following file for general binary output : '
WRITE(kStdWarn,*) caOutName

IF (kLongOrShort == 0) THEN
  WRITE(iIOUN1) caVersion               !!kcarta version number
  WRITE(iIOUN1) kProfLayer
  WRITE(iIOUN1) kMaxUserSet
  WRITE(iIOUN1) (raParams(iI),iI=1,kMaxUserSet)
  WRITE(iIOUN1) caComment
  WRITE(iIOUN1) rFrLow,rFrHigh
  WRITE(iIOUN1) iFileIDLo,iFileIDHi
  WRITE(iIOUN1) kLongOrShort
!!!this is new stuff, to tell reader how many chunks to read
  WRITE(iIOUN1) kaFrStep(iTag)          !!freq step size
  WRITE(iIOUN1) kaBlSize(iTag)          !!10000 point freq block size
  WRITE(iIOUN1) iTotalStuff             !!number of outputs
  
  WRITE(kStdWarn,*) 'doing the SHORT BASIC version : ......'
  WRITE(kStdWarn,*) '  '
  WRITE(kStdWarn,*) caVersion               !!kcarta version number
  WRITE(kStdWarn,*) rFrLow,rFrHigh          !!start and stop wavenumber
  WRITE(kStdWarn,*) kaFrStep(iTag)          !!freq step size
  WRITE(kStdWarn,*) kaBlSize(iTag)          !!10000 point freq block size
  WRITE(kStdWarn,*) iImportant              !!number of outputs
  
  GO TO 9999
END IF

!    or if kLongOrShort .EQ. -1,+1 then do the standard binary file
! first output general info -----------------------------------------
WRITE(iIOUN1) caVersion
WRITE(iIOUN1) kProfLayer
WRITE(iIOUN1) kMaxUserSet
WRITE(iIOUN1) (raParams(iI),iI=1,kMaxUserSet)
WRITE(iIOUN1) caComment
WRITE(iIOUN1) rFrLow,rFrHigh
WRITE(iIOUN1) iFileIDLo,iFileIDHi
WRITE(iIOUN1) kLongOrShort
! v1.04 to v1.08 had the next line; now remove it
!      WRITE(iIOUN1) M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
WRITE(iIOUN1) (raPactualAvg(iI),iI=1,kProfLayer)
WRITE(iIOUN1) (iaaOverrideDefault(1,iI),iI=1,10)    !! general settings
WRITE(iIOUN1) (iaaOverrideDefault(2,iI),iI=1,10)    !! radtrans settings
WRITE(iIOUN1) (iaaOverrideDefault(3,iI),iI=1,10)    !! iLBLRTM settings

! then output path ID stuff ------------------------------------------
IF (kLongOrShort > 0) THEN
  WRITE(iIOUN1) iNumGases*kProfLayer
  DO iI=1,iNumGases
    DO iJ=1,kProfLayer
      iP=(iI-1)*kProfLayer+iJ
      WRITE(iIOUN1)iP,iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
    END DO
  END DO
END IF
! then output list of paths to be output
iP=0
DO iI=1,iOutTypes
  IF (iaPrinter(iI) == 1) THEN
    iP=iP+1
    IF (iaNp(iI) > iNumGases*kProfLayer) THEN
      WRITE(kStdErr,*)'Error! Cannot have more paths to be output than '
      WRITE(kStdErr,*)'iNumGases*kProfLayer'
      CALL DoSTOP
    END IF
    IF (iaNp(iI) > 0) THEN
      iEnd=iaNp(iI)
      WRITE(iIOUN1) iEnd
      WRITE(iIOUN1) (iaaOp(iI,iJ),iJ=1,iEnd)
    ELSE
      iEnd = kProfLayer*iNumGases
      WRITE(iIOUN1) iEnd
      WRITE(iIOUN1) (iaaOp(iI,iJ),iJ=1,iEnd)
    END IF
    iaOutNumbers(iI)=iEnd
  END IF
END DO
IF (iP > 1) THEN
  WRITE(kStdErr,*)'Have found more than 1 set of options in *OUTPUT'
  WRITE(kStdErr,*)'where iDat set to 1 (outputting single paths)'
  WRITE(kStdErr,*)'... exiting program!'
  CLOSE(iIOUN1)
  CALL DoSTOP
END IF
IF (iP == 0) THEN
! inform user no single set of paths to be output
  WRITE(iIOUN1) iP
END IF

! then output mixed paths -----------------------------------------------
! note this important CHANGE !!!
WRITE(iIOUN1) iNpmix
IF (kLongOrShort > 0) THEN
  WRITE(iIOUN1) iMixFileLines
  IF (iNpMix > 0) THEN
    DO iI=1,iMixFileLines
      WRITE(iIOUN1) caaMixFileLines(iI)
    END DO
! then output list of mixed paths temperatures
    WRITE(iIOUN1) (raMixVT(iI),iI=1,iNpmix)
  END IF
END IF

! this is EXTRA if statement, since above one was commented out
IF (iNpMix > 0) THEN
! then output list of mixed paths to be printed
  iP=0
  DO iI=1,iOutTypes
    IF (iaPrinter(iI) == 2) THEN
      iP=iP+1
      IF (iaNp(iI) > iNpmix) THEN
        WRITE(kStdErr,*)'Error! Cannot print out more mixed paths'
        WRITE(kStdErr,*)'than iIpmix'
        CALL DoSTOP
      END IF
      IF (iaNp(iI) > 0) THEN
        iEnd=iaNp(iI)
        WRITE(iIOUN1) iEnd
        WRITE(iIOUN1) (iaaOp(iI,iJ),iJ=1,iEnd)
      ELSE
        iEnd=iNpmix
        WRITE(iIOUN1) iEnd
        WRITE(iIOUN1) (iaaOp(iI,iJ),iJ=1,iEnd)
      END IF
      iaOutNumbers(iI)=iEnd
    END IF
  END DO
  
  IF (iP > 1) THEN
    WRITE(kStdErr,*)'Have found more than 1 options set in *OUTPUT where'
    WRITE(kStdErr,*)'iDat set to 3 (outputting mixed paths)'
    WRITE(kStdErr,*)'... exiting program!'
    CLOSE(iIOUN1)
    CALL DoSTOP
  END IF
  
  IF (iP == 0) THEN
! inform user no mixed set of paths to be output
    WRITE(iIOUN1) iP
  END IF
END IF

! now print the atmosphere information -------------------------------------
WRITE(iIOUN1) iNatm
WRITE(iIOUN1) kEmsRegions
DO iI=1,iNatm
  IF (iaNumLayers(iI) > iNpmix) THEN
    WRITE(kStdErr,*)'atm#',iI,' has too many layers!! (> ',iNpmix,')'
    CALL DoSTOP
  END IF
  IF (iaNumLayers(iI) < 0) THEN
    WRITE(kStdErr,*)'atm#',iI,' has 0 layers'
    CALL DoSTOP
  END IF
  WRITE(iIOUN1) iI,iaNumLayers(iI)
  WRITE(iIOUN1) (iaaRadLayer(iI,iJ),iJ=1,iaNumLayers(iI))
  WRITE(iIOUN1) raTSpace(iI),raTSurf(iI),raSatAngle(iI), raSatHeight(iI)
  WRITE(iIOUN1) iakSolar(iI),rakSolarAngle(iI),rakSolarRefl(iI),  &
      iakThermal(iI),rakThermalAngle(iI),iakThermalJacob(iI)
  
! output the emissivities for this atmosphere, first freq pt then ems
  WRITE(iIOUN1) iaSetEms(iI)
  DO iJ=1,iaSetEms(iI)
    WRITE(iIOUN1)raaaSetEmissivity(iI,iJ,1), raaaSetEmissivity(iI,iJ,2)
  END DO
! now output the list of mixed paths to be printed, for this atmosphere
  iP=0
  DO iJ=1,iOutTypes
    IF ((iaPrinter(iJ) == 3).AND.(iaAtmPr(iJ) == iI)) THEN
      iP=iP+1
      IF (iaNp(iJ) > iNpmix) THEN
!            IF (iaNp(iJ) .GT. iaNumLayers(iI)) THEN
        WRITE(kStdErr,*)'Atm #',iI,' has too many layers to be printed!'
        CALL DoSTOP
      END IF
      IF (iaNp(iJ) == 0) THEN
        WRITE(kStdErr,*)'Atm #',iI,' has 0 layers to be printed!'
        CALL DoSTOP
      END IF
      IF (iaNp(iJ) > 0) THEN
        iEnd=iaNp(iJ)
        iaOutNumbers(iJ)=iEnd
        WRITE(iIOUN1) iEnd
        WRITE(iIOUN1) (iaaOp(iJ,iK),iK=1,iEnd)
        WRITE(iIOUN1) (raaUserPress(iJ,iK),iK=1,iEnd)
        iOk=1
        DO iK=1,iEnd
          IF ((iaaOp(iJ,iK) > iaNumLayers(iI)) .OR.  &
                (iaaOp(iJ,iK) <= 0))  THEN
            iOk=-1
          END IF
        END DO
        IF (iOk < 0) THEN
          WRITE(kStdErr,*)'Atm# ',iI,' has invalid layer to be output!'
          CALL DoSTOP
        END IF
      ELSE !!this is the case when ALL levels need to be printed
        iEnd=iaNumLayers(iI)
        WRITE(iIOUN1) iEnd
        iaOutNumbers(iJ)=iEnd
!note how we are flipping between iI and iJ here
!and instead of outting iaaOp, we are outputting iaaRadLayer
        WRITE(iIOUN1) (iaaRadLayer(iI,iK),iK=1,iEnd)
        WRITE(iIOUN1) (raaUserPress(iJ,iK),iK=1,iEnd)
        
        IF ((iNumNLTEGases > 0) .AND. (iDoUpperAtmNLTE > 0)) THEN
          iDumpAllUARads = iDumpAllUARads + 1
          WRITE(kStdWarn,*) iJ,iaNp(iJ)
          WRITE(kStdWarn,*) 'looks like we will dump UA rads'
        END IF
        IF ((iDumpAllUARads > 1) .AND.(iDoUpperAtmNLTE <= 0)) THEN
          WRITE(kStdErr,*) 'this is too complicated!!!'
          WRITE(kStdErr,*) 'NLTE code assumes ONE radiating atmosphere'
          CALL DoStop
        END IF
        
      END IF
    END IF
  END DO
  IF (iP == 0) THEN
! inform user no mixed paths to be output for this atmosphere
    WRITE(iIOUN1) iP
  END IF
  IF (iP > 1) THEN
    WRITE(kStdErr,*)'Have found more than 1 set of options in *OUTPUT '
    WRITE(kStdErr,*)'where iDat set to 3 (outputting radiances)'
    WRITE(kStdErr,*)'for atmosphere #',iI,'... exiting program!'
    CLOSE(iIOUN1)
    CALL DoSTOP
  END IF
END DO

! having gotten this far, this tells the reader how many things to expect
! total num kcomp files, total number of output options
! for each output option, how many prints to expect
WRITE(iIOUN1) iTotal,iOutTypes
WRITE(iIOUN1) (iaOutNumbers(iI),iI=1,iOutTypes)

!new do not open and close!!
!      IF (iIOUN1 .NE. 6) THEN
!        CLOSE(iIOUN1)
!      END IF

4000 FORMAT(A130)

!--------------- COL JACOBIAN BINARY FILE --------------------------------
IF (kJacobian >= 0 .AND. kActualJacs >= 100) THEN
  WRITE(kStdWarn,*) 'opening column,stemp (small) jacobian output file'
  kStdJacob2 = kStdJacob2KK
  iIOUN_JAC2 = kStdJacob2
  
  IF (iOutFileName < 0) THEN  !std kcarta jacob2 output dumped to screen
    caJacobFile2 = 'columnjacob.dat'
  ELSE !std kcarta output dumped to file, so find flux file name
    CALL jacob2Name(caJacobFile2,caJacobFile)
  END IF
  
  OPEN(UNIT=iIOUN_JAC2,FILE=caJacobFile2,STATUS='NEW',  &
      FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,403) iFileErr, caJacobFile2
    WRITE(kStdErr,*)'make sure the file does not exist!'
    CALL DoSTOP
  END IF
  
  403  FORMAT('ERROR! number ',I5,' opening COL JACOBIAN binary file :  &
      ',/,A80)
  
  kStdJacob2Open=1
  WRITE(kStdWarn,*) 'Opened file for col jacobian binary output : '
  WRITE(kStdWarn,*) caJacobFile2
  
! now essentially dump out header that resmebles kLongOrShort == 0
  WRITE(iIOUN_Jac2) caVersion               !!kcarta version number
  WRITE(iIOUN_Jac2) kProfLayer
  WRITE(iIOUN_Jac2) kMaxUserSet
  WRITE(iIOUN_Jac2) (raParams(iI),iI=1,kMaxUserSet)
  WRITE(iIOUN_Jac2) caComment
  WRITE(iIOUN_Jac2) rFrLow,rFrHigh
  WRITE(iIOUN_Jac2) iFileIDLo,iFileIDHi
  WRITE(iIOUN_Jac2) 0            !! should be kLongOrShort, but use 0
!!!this is new stuff, to tell reader how many chunks to read
  WRITE(iIOUN_Jac2) kaFrStep(iTag)          !!freq step size
  WRITE(iIOUN_Jac2) kaBlSize(iTag)          !!10000 point freq block size
  iImportant=iJacob+1+1                     !!r(gasQ'),r(temp'),r(stemp')
  WRITE(iIOUN_Jac2) iImportant              !!number of outputs
  
END IF
!--------------- JACOBIAN BINARY FILE --------------------------------
IF (kJacobian >= 0 .AND. kActualJacs < 100) THEN
  WRITE(kStdWarn,*) 'opening regular (large) jacobian output file (X(z))'
  WRITE(kStdWarn,*) 'kJacobOutput = ',kJacobOutput
  iIOUN2 = kStdJacob
  
  IF (iIOUN2 /= 6) THEN
! open unformatted file as a fresh file to be written to
    OPEN(UNIT=iIOUN2,FILE=caJacobFile,STATUS='NEW',  &
        FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
    IF (iFileErr /= 0) THEN
      WRITE(kStdErr,203) iFileErr, caJacobFile
      WRITE(kStdErr,*)'make sure the file does not exist!'
      CALL DoSTOP
    END IF
  END IF
  
  203  FORMAT('ERROR! number ',I5,' opening JACOBIAN binary file : ',/,A80)
  
  kStdJacobOpen=1
  WRITE(kStdWarn,*) 'Opened following file for jacobian binary output : '
  WRITE(kStdWarn,*) caJacobFile
  
! write general header information
  WRITE(iIOUN2) caVersion
  WRITE(iIOUN2) caComment
  WRITE(iIOUN2) kProfLayer
  WRITE(iIOUN2) rFrLow,rFrHigh
  WRITE(iIOUN2) iFileIDLo,iFileIDHi
  
! write some specific info, which will be used for each 25 cm-1 chunk,
! for each atmosphere
  
! first figure out how many gasid's to output
  iImportant=iJacob
  WRITE(iIOUN2) iImportant
  WRITE(kStdWarn,*) ' iImportant gases = ',iImportant
  
! then figure out, of the atmospheres that have been read in, which actually
! have a radiance and hence jacobian calculation associated with them
! assume all error checking done in section above (when blah.dat is created)
  iNatmJac=0
  DO iI=1,iNatm
! now output the list of mixed paths to be printed, for this atmosphere
    DO iJ=1,iOutTypes
      IF ((iaPrinter(iJ) == 3).AND.(iaAtmPr(iJ) == iI)) THEN
        iNatmJac=iNatmJac+1
        IF ((iaNp(iJ) > 0) .OR. (iaNp(iJ) == -1)) THEN
! set the number of layers in this atmosphere
          iaLayerJac(iNatmJac)=iaNumLayers(iI)
        END IF
      END IF
    END DO
  END DO
  
  WRITE(iIOUN2) iNatmJac
  WRITE(iIOUN2) (iaLayerJac(iI),iI=1,iNatmJac)
  
  WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
  WRITE(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)
  
  iI=1
  205    CONTINUE
  iP=DoGasJacob(iaGases(iI),iaJacob,iJacob)
  IF ((iI <= iNumGases) .AND. (iP > 0)) THEN
    DO iJ=1,kProfLayer
      WRITE(iIOUN2)iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
    END DO
    iI=iI+1
    GO TO 205
  END IF
  IF (iI < iNumGases) THEN
    iI=iI+1
    GO TO 205
  END IF
  
!!check to see if we are outputting d/d(DME), d/d(IWP)
  IF ((kJacobian > 0) .AND. (kScatter > 0)) THEN
    DO iJ=1,kProfLayer
      WRITE(iIOUN2) 201,-300.0,01.0
    END DO
    DO iJ=1,kProfLayer
      WRITE(iIOUN2) 201,-300.0,01.0
    END DO
  END IF
  
! recall that at the end, we also compute d/dSurface_temp,d/dSurface_Emis
! and d Thermal BackGnd/d(Surface_Emis)
! and d SolarBackGnd/d(Surface_Emis)
  
! having gotten this far, this tells the reader how many things to expect
! total num kcomp files, total number of output options=iNatm,number of gases
! for each atmosphere, how many layers
  WRITE(iIOUN2) iTotal,iNatmJac,iImportant
  WRITE(iIOUN2) (iaLayerJac(iI),iI=1,iNatmJac)
  
END IF

!--------------- FLUX BINARY FILE -----------------------------------------
!----- the opening header info is almost same as that for Jacobian files ---
IF (kFlux > 0) THEN
  
  iIOUN_Flux = kStdFlux
  
  IF (iOutFileName < 0) THEN  !std kcarta output dumped to screen
    caFluxFile='flux.dat'
  ELSE !std kcarta output dumped to file, so find flux file name
    CALL FluxName(caFluxFile,caOutName)
  END IF
  
! open unformatted file as a fresh file to be written to
  OPEN(UNIT=iIOUN_Flux,FILE=caFluxFile,STATUS='NEW',  &
      FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,303) iFileErr, caFluxFile
    WRITE(kStdErr,*)'make sure the file does not exist!'
    CALL DoSTOP
  END IF
  
  303  FORMAT('ERROR! number ',I5,' opening FLUX binary file : ',/,A80)
  
  kStdFluxOpen=1
  WRITE(kStdWarn,*) 'Opened following file for flux ouput : '
  WRITE(kStdWarn,*) caFluxFile
  
! write general header information
  WRITE(iIOUN_Flux) caComment
!        IF ((kFlux .NE. 6) .OR. (kFlux .NE. 1,2,3)) THEN
!          WRITE(iIOUN_Flux) kProfLayer
!        ELSE
!          WRITE(iIOUN_Flux) kProfLayer+1
!        END IF
  WRITE(iIOUN_Flux) kProfLayer
  WRITE(iIOUN_Flux) rFrLow,rFrHigh
  WRITE(iIOUN_Flux) iFileIDLo,iFileIDHi
!cccccc this is from Jacobians
!cc first figure out how many gasid's to output
!cc        iImportant=iJacob
!cc        WRITE(iIOUN_Flux) iImportant
!cc        write(kStdWarn,*) ' iImportant gases = ',iImportant
! figure out how many types of fluxes to output (duh!!!!!!!)
  iImportant = 1
  WRITE(iIOUN_Flux) iImportant
  WRITE(kStdWarn,*) 'Up-Down Fluxes = ',iImportant
  
! then figure out, of the atmospheres that have been read in, which actually
! have a radiance and hence jacobian calculation associated with them
! assume all error checking done in section above (when blah.dat is created)
  iNatmJac = 0
  DO iI=1,iNatm
! now output the list of mixed paths to be printed, for this atmosphere
    DO iJ=1,iOutTypes
      IF ((iaPrinter(iJ) == 3).AND.(iaAtmPr(iJ) == iI)) THEN
        iNatmJac = iNatmJac+1
        IF (iaNp(iJ) /= 0) THEN   !!ie it is -1 or a positive number
!               IF (iaNp(iJ) .GT. 0) THEN
! set the number of layers in this atmosphere
          iaLayerJac(iNatmJac) = iaNumLayers(iI)
        END IF
      END IF
    END DO
  END DO
  
  IF (kFLux <= 3)  THEN   !!! outputting at all levels
    DO iI = 1,iNatmJac
      iaJunkFlux(iI) = 1 * (iaLayerJac(iI)+1)
    END DO
    
    WRITE(iIOUN_Flux) iNatmJac
    WRITE(iIOUN_Flux) (iaJunkFlux(iI),iI=1,iNatmJac)
    
    WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
    WRITE(kStdWarn,*) (iaJunkFlux(iI),iI=1,iNatmJac)
    
! how many things to expect :
! total num kcomp files, total number of output options=iNatm,number of fluxes
! for each atmosphere, dump OLR AND ILR
    WRITE(iIOUN_Flux) iTotal,iNatmJac,iImportant
    WRITE(iIOUN_Flux) (iaJunkFLux(iI),iI=1,iNatmJac)
    
  ELSE IF (kFLux == 4) THEN   !!! outputting only OLR at TOA
    DO iI = 1,iNatmJac
      iaLayerFlux(iI) = 1
    END DO
    
    WRITE(iIOUN_Flux) iNatmJac
    WRITE(iIOUN_Flux) (iaLayerFlux(iI),iI=1,iNatmJac)
    
    WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
    WRITE(kStdWarn,*) (iaLayerFlux(iI),iI=1,iNatmJac)
    
! no need to dump out gas amounts, temperatures; so now tell the reader how
! many things to expect :
! total num kcomp files, total number of output options=iNatm,number of fluxes
! for each atmosphere, how many layers
    WRITE(iIOUN_Flux) iTotal,iNatmJac,iImportant
    WRITE(iIOUN_Flux) (iaLayerFLux(iI),iI=1,iNatmJac)
    
  ELSE IF (kFLux == 5) THEN   !!! outputting only OLR/ILR at TOA/GND/TRPause
    DO iI = 1,iNatmJac
      iaLayerFlux(iI) = 3
    END DO
    
    WRITE(iIOUN_Flux) iNatmJac
    WRITE(iIOUN_Flux) (iaLayerFlux(iI),iI=1,iNatmJac)
    
    WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
    WRITE(kStdWarn,*) (iaLayerFlux(iI),iI=1,iNatmJac)
    
! no need to dump out gas amounts, temperatures; so now tell the reader how
! many things to expect :
! total num kcomp files, total number of output options=iNatm,number of fluxes
! for each atmosphere, dump OLR AND ILR
    WRITE(iIOUN_Flux) iTotal,iNatmJac,iImportant
    WRITE(iIOUN_Flux) (iaLayerFLux(iI),iI=1,iNatmJac)
    
  ELSE IF (kFLux == 6) THEN   !!! outputing OLR/ILR at each layer, plus GND and TOA
    DO iI = 1,iNatmJac
      iaJunkFlux(iI) = 2 * (iaLayerJac(iI)+1)
    END DO
    
    WRITE(iIOUN_Flux) iNatmJac
    WRITE(iIOUN_Flux) (iaJunkFlux(iI),iI=1,iNatmJac)
    
    WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
    WRITE(kStdWarn,*) (iaJunkFlux(iI),iI=1,iNatmJac)
    
! how many things to expect :
! total num kcomp files, total number of output options=iNatm,number of fluxes
! for each atmosphere, dump OLR AND ILR
    WRITE(iIOUN_Flux) iTotal,iNatmJac,iImportant
    WRITE(iIOUN_Flux) (iaJunkFLux(iI),iI=1,iNatmJac)
  END IF
  
END IF

!--------------- PLANCK BINARY FILE -----------------------------------------
!----- the opening header info is almost same as that for Jacobian files ---
!----- control is thru kFlux (-1,+1,2,3,4 for off, 0 for on)
kPlanckOut = kFlux
IF ((kPlanckOut == 0) .AND. (iNumNLTEGases > 0)) THEN
  
  iIOUN_Planck = kStdPlanck
  
  IF (iOutFileName < 0) THEN  !std kcarta output dumped to screen
    caPlanckFile='planck.dat'
  ELSE !std kcarta output dumped to file, so find planck file name
    CALL PlanckName(caPlanckFile,caOutName)
  END IF
  
! open unformatted file as a fresh file to be written to
  OPEN(UNIT=iIOUN_Planck,FILE=caPlanckFile,STATUS='NEW',  &
      FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
  IF (iFileErr /= 0) THEN
    WRITE(kStdErr,304) iFileErr, caPlanckFile
    WRITE(kStdErr,*)'make sure the file does not exist!'
    CALL DoSTOP
  END IF
  
  304    FORMAT('ERROR! number ',I5,' opening PLANCK binary file : ',/,A80)
  
  kStdPlanckOpen=1
  WRITE(kStdWarn,*) 'Opened following file for planck output :'
  WRITE(kStdWarn,*) caPlanckFile
  
! write general header information
  WRITE(iIOUN_Planck) caComment
  WRITE(iIOUN_Planck) kProfLayer
  WRITE(iIOUN_Planck) rFrLow,rFrHigh
  WRITE(iIOUN_Planck) iFileIDLo,iFileIDHi
!cccccc this is from Jacobians
!cc first figure out how many gasid's to output
!cc        iImportant=iJacob
!cc        WRITE(iIOUN_Planck) iImportant
!cc        write(kStdWarn,*) ' iImportant gases = ',iImportant
! figure out how many types of plancks to output (duh!!!!!!!)
  iImportant=1
  WRITE(iIOUN_Planck) iImportant
  
! then figure out, of the atmospheres that have been read in, which actually
! have a radiance calculation associated with them
! assume all error checking done in section above (when blah.dat is created)
!        print *,iNatm
!        print *,(iaNumLayers(iI),iI=1,iNatm)
!        print *,iOutTypes
!        print *,(iaPrinter(iI),iI=1,iOutTypes)
!        DO iJ=1,iOutTypes
!          print *,iJ,iaPrinter(iJ),iaAtmPr(iJ),iaNp(iJ)
!        END DO
!        stop
  
  iNatmJac=0
  DO iI=1,iNatm
! now output the list of mixed paths to be printed, for this atmosphere
    DO iJ=1,iOutTypes
      IF ((iaPrinter(iJ) == 3) .AND. (iaAtmPr(iJ) == iI)) THEN
        iNatmJac=iNatmJac+1
!              IF (iaNp(iJ) .GT. 0) THEN
        IF (iaNp(iJ) /= 0) THEN   !!ie it is -1 or a positive number
! set the number of layers in this atmosphere
          iaLayerJac(iNatmJac)=iaNumLayers(iI)
        END IF
      END IF
    END DO
  END DO
  
  WRITE(iIOUN_Planck) iNatmJac
  WRITE(iIOUN_Planck) (iaLayerJac(iI),iI=1,iNatmJac)
  
  WRITE(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
  WRITE(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)
  
! no need to dump out gas amounts, temperatures; so now tell the reader how
! many things to expect :
! total num kcomp files, total num of output options=iNatm,number of planckes
! for each atmosphere, how many layers
  WRITE(iIOUN_Planck) iTotal,iNatmJac,iImportant
  WRITE(iIOUN_Planck) (iaLayerJac(iI),iI=1,iNatmJac)
  
END IF

9999 CONTINUE     !would have jumped here if kLongShort = 0

IF (iDumpAllUARads <= 0) THEN
  iDumpAllUARads = -1
END IF

! f77 format specifiers : the 1P is to shift things over by 1 decimal point
!                          but unfortunately it moves everything, so have to negate with 0P, SAFER TO use ES instead
! http://docs.oracle.com/cd/E19957-01/805-4939/z40007437a2e/index.html
! ../DOC/F77FormatSpecifiers.pdf
7169 FORMAT(A104)
! 7170 FORMAT(I4,' ',I4,' ',1(F11.5,' '),2(1P E11.4,' '),0PF11.5,'  ',(1P E11.4),A2,2(1P E11.4,' '),0PF11.7)
7170 FORMAT(I4,' ',I4,' ',1(F11.5,' '),2(ES11.4,' '),F11.5,'  ',(ES11.4),A2,2(ES11.4,' '),F11.7)
7171 FORMAT(I4,' ',I4,' ',3(1P E11.5,' '),0PF11.5,'  ',1P E11.5,A2,2(' ',E11.5),1(' ',E11.3))
7172 FORMAT(I3,' ',I3,' ',A20,' ',ES11.5)

RETURN
END SUBROUTINE PrepareOutput

!************************************************************************
! this subroutine opens the output UA file

! very similar to OpenBloatFile

! only thing is, it dumps out (iUpper + 1) number of things
! if the user asks for all CO2 opt depths to be dumped out, then also the
!    UA opt depth spectra are dumped out!!!!!!!!!
! if only specific rads are being dumped out in the regular AIRS atm (0-80 km)
!    then TWO rads are output in this file, one at the regular TOA (0.005 mb)
!    and the other at the new TOA (0.000025 mb)
!    this is if iDumpAllUARads == -1
! if rads at all layers are being dumped out in the regular AIRS atm (0-80 km)
!    then rads at all layers are also output in this file, plus one at the
!    regular TOA (0.005 mb)
!    this is if iDumpAllUARads == +1

SUBROUTINE OpenUAFile(iType,iUpperLayers,caUAFile,rFrLow,rFrHigh,  &
    iDumpAllUASpectra,iDumpAllUARads, iFileIDLo,iFileIDHi,iTag,  &
    iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)


INTEGER, INTENT(IN OUT)                  :: iType
NO TYPE, INTENT(IN OUT)                  :: iUpperLaye
CHARACTER (LEN=80), INTENT(OUT)          :: caUAFile
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(OUT)                     :: iFileIDLo
INTEGER, INTENT(OUT)                     :: iFileIDHi
INTEGER, INTENT(OUT)                     :: iTag
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'



INTEGER :: iDumpAllUARads,iDumpAllUASpectra

INTEGER :: iUpperLayers !!number of layers in atm #1

INTEGER :: iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

! local vars
INTEGER :: iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor,iTotalStuff
REAL :: r1,r2,r3

iTotalStuff = 0

IF (iDumpAllUASpectra > 0) THEN
  iTotalStuff = iUpperLayers          !!!dump out all the UA opt depths
END IF

IF (iDumpAllUARads > 0) THEN
!!!rads at 0.005 mb and at all layers upto 0.000025 mb
  iTotalStuff = iTotalStuff + (1 + iUpperLayers)
ELSE
!!! default = rads at 0.005 mb and at 0.000025 mb
  iTotalStuff = iTotalStuff + (1 + 1)
END IF

! find number of chunks to dump
iNLTEChunks = -1
i1 = +123456
i2 = -123456
DO iI = 1,iNumNLTEGases
  IF (iaNLTEChunks(iI) >= iNLTEChunks) THEN
    iNLTEChunks = iaNLTEChunks(iI)
  END IF
  DO iJ = 1,iaNLTEChunks(iI)
    IF (iaaNLTEChunks(iI,iJ) <= i1) THEN
      i1 = iaaNLTEChunks(iI,iJ)
    END IF
    IF (iaaNLTEChunks(iI,iJ) >= i2) THEN
      i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
    END IF
  END DO
END DO
r1 = i1*1.0
r2 = i2*1.0
iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1

!!!now check to see we actually start from r1, while running kCARTA!
IF (r1 < rFrLow) THEN
  r1 = rFrLow
  i1 = iFloor(r1)
  iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
END IF

!!!now check to see we actually do go to r2, while running kCARTA!
!      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
!        r2 = rFrHigh-kaBlSize(iTag)
IF (r2 > rFrHigh) THEN
  r2 = rFrHigh
  i2 = iFloor(r2)
  iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
END IF

IF (iType > 0) THEN
  iIoun1 = kNLTEOutUA
ELSE IF (iType < 0) THEN
  iIoun1 = kStdPlanckUA
ELSE
  WRITE(kStdErr,*) 'Unknown print option for UA files'
  CALL DoStop
END IF

! open unformatted file as a fresh file to be written to
OPEN(UNIT=iIoun1,FILE=caUAFile,STATUS='NEW',  &
    FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
IF (iFileErr /= 0) THEN
  WRITE(kStdErr,304) iFileErr,iIOUN1,caUAFile
  WRITE(kStdErr,*)'make sure the file does not exist!'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) 'Printing UA FILE summary stats at open : '
WRITE(kStdWarn,*) caUAFile
WRITE(kStdWarn,*) 'Regular freqs etc .... '
WRITE(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
WRITE(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
WRITE(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
WRITE(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff

WRITE(kStdWarn,*) i1,i2,iNLTEChunks   !!which chunks we are dealing with
!!and how many to expect
WRITE(kStdWarn,*) r1,r2               !!start/stop freqs for NLTE chunks

304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file :  &
    ',/,A80)

IF (iType > 0) THEN
  kNLTEOutUAOpen = 1
  WRITE(kStdWarn,*) 'Opened following file for UA general output :'
  WRITE(kStdWarn,*) caUAFile
ELSE IF (iType < 0) THEN
  kStdPlanckUAOpen = 1
  WRITE(kStdWarn,*) 'Opened following file for UA planck output :'
  WRITE(kStdWarn,*) caUAFile
END IF

WRITE(iIOUN1) kProfLayer
WRITE(iIOUN1) rFrLow,rFrHigh
WRITE(iIOUN1) iFileIDLo,iFileIDHi     !!user should expect x5 of these!
WRITE(iIOUN1) kaFrStep(iTag)          !!regular freq step size
WRITE(iIOUN1) kaBlSize(iTag)          !!regular 10k point freq block size
WRITE(iIOUN1) iTotalStuff             !!regular number of outputs
WRITE(iIOUN1) iType                   !!regular or Planck file

WRITE(iIOUN1) i1,i2,iNLTEChunks       !!which chunks we are dealing with
!!and how many to expect
WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks


RETURN
END SUBROUTINE OpenUAFile

!************************************************************************
! this subroutine opens the output UA file, for bloated version.
! HAS NOT BEEN TESTED AT ALL

! very similar to OpenBloatFile

! only thing is, it dumps out (iUpper + 1) number of things
! if only specific rads are being dumped out in the regular AIRS atm (0-80 km)
!    then TWO rads are output in this file, one at the regular TOA (0.005 mb)
!    and the other at the new TOA (0.000025 mb)
!    this is if iDumpAllUARads == -1
! if rads at all layers are being dumped out in the regular AIRS atm (0-80 km)
!    then rads at all layers are also output in this file, plus one at the
!    regular TOA (0.005 mb)
!    this is if iDumpAllUARads == +1

SUBROUTINE OpenBloatUAFile(  &
    iType,iUpperLayers,caOutUABloatFile,rFrLow,rFrHigh, iDumpAllUARads,  &
    iFileIDLo,iFileIDHi,iTag, iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)


INTEGER, INTENT(OUT)                     :: iType
NO TYPE, INTENT(IN OUT)                  :: iUpperLaye
NO TYPE, INTENT(IN OUT)                  :: caOutUABlo
REAL, INTENT(IN OUT)                     :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(OUT)                     :: iFileIDLo
INTEGER, INTENT(OUT)                     :: iFileIDHi
INTEGER, INTENT(OUT)                     :: iTag
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iaNLTEChun
NO TYPE, INTENT(IN OUT)                  :: iaaNLTEChu
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

CHARACTER (LEN=80) :: caOutUABloatFile

INTEGER :: iDumpAllUARads

INTEGER :: iUpperLayers !!number of layers in atm #1

INTEGER :: iaNLTEChunks(kGasStore)
INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

! local vars
INTEGER :: iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor,iTotalStuff
REAL :: r1,r2,r3

iTotalStuff = 2         !!! default = rads at 0.005 mb and at 0.000025 mb
IF (iDumpAllUARads > 0) THEN
  iTotalStuff = 1 + iUpperLayers      !!!rads at 0.005 mb and at all
!!!layers upto 0.000025 mb
END IF

! find number of chunks to dump
iNLTEChunks = -1
i1 = +123456
i2 = -123456
DO iI = 1,iNumNLTEGases
  IF (iaNLTEChunks(iI) >= iNLTEChunks) THEN
    iNLTEChunks = iaNLTEChunks(iI)
  END IF
  DO iJ = 1,iaNLTEChunks(iI)
    IF (iaaNLTEChunks(iI,iJ) <= i1) THEN
      i1 = iaaNLTEChunks(iI,iJ)
    END IF
    IF (iaaNLTEChunks(iI,iJ) >= i2) THEN
      i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
    END IF
  END DO
END DO
r1 = i1*1.0
r2 = i2*1.0
iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1

!!!now check to see we actually start from r1, while running kCARTA!
IF (r1 < rFrLow) THEN
  r1 = rFrLow
  i1 = iFloor(r1)
  iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
END IF

!!!now check to see we actually do go to r2, while running kCARTA!
!      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
!        r2 = rFrHigh-kaBlSize(iTag)
IF (r2 > rFrHigh) THEN
  r2 = rFrHigh
  i2 = iFloor(r2)
  iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
END IF

IF (iType > 0) THEN
  iIoun1 = kBloatNLTEOutUA
ELSE IF (iType < 0) THEN
  iIoun1 = kBloatPlanckUA
ELSE
  WRITE(kStdErr,*) 'Unknown print option for bloat UA files'
  CALL DoStop
END IF

! open unformatted file as a fresh file to be written to
OPEN(UNIT=iIoun1,FILE=caOutUABloatFile,STATUS='NEW',  &
    FORM='UNFORMATTED',IOSTAT=iFileErr)
! if file error, inform user and stop program
IF (iFileErr /= 0) THEN
  WRITE(kStdErr,304) iFileErr,iIOUN1,caOutUABloatFile
  WRITE(kStdErr,*)'make sure the file does not exist!'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) 'Printing Bloat UA FILE summary stats at open : '
WRITE(kStdWarn,*) caOutUABloatFile
WRITE(kStdWarn,*) 'Regular freqs etc .... '
WRITE(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
WRITE(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
WRITE(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
WRITE(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff
WRITE(kStdWarn,*) 'Bloated  freqs etc .... '
WRITE(kStdWarn,*) '  iType = = ',iType
WRITE(kStdWarn,*) '  iNumlayers = ',iUpperLayers
WRITE(kStdWarn,*) '  kBoxCarUse = ',kBoxCarUse
WRITE(kStdWarn,*) '  df_fine = ',kaFineFrStep(iTag)
WRITE(kStdWarn,*) 'start,stop high res integers : ',i1,i2
WRITE(kStdWarn,*) 'number of highres chunks = ',iNLTEChunks
WRITE(kStdWarn,*) 'start,stop high res doubles : ',r1,r2
WRITE(kStdWarn,*) ' '

304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file :  &
    ',/,A80)

IF (iType > 0) THEN
  kBloatNLTEOutUAOpen = 1
  WRITE(kStdWarn,*) 'Opened following file for UA bloat general output :'
  WRITE(kStdWarn,*) caOutUABloatFile
ELSE IF (iType < 0) THEN
  WRITE(kStdErr,*) 'ooer cannot open bloated UA planck files'
  CALL DoStop
  kBloatPlanckUAOpen = 1
  WRITE(kStdWarn,*) 'Opened following file for UA bloat planck output :'
!        write(kStdWarn,*) caOutUABloatPlanckFile
END IF

WRITE(iIOUN1) kProfLayer
WRITE(iIOUN1) rFrLow,rFrHigh
WRITE(iIOUN1) iFileIDLo,iFileIDHi     !!user should expect x5 of these!
WRITE(iIOUN1) kaFrStep(iTag)          !!regular freq step size
WRITE(iIOUN1) kaBlSize(iTag)          !!regular 10k point freq block size
WRITE(iIOUN1) iTotalStuff             !!regular number of outputs
!!!!!now do the fine stuff
WRITE(iIOUN1) iType                   !!tells reader if reg +1 or plck -1
WRITE(iIOUN1) iUpperLayers            !!tells #of layers in atm 1
WRITE(iIOUN1) kBoxCarUse              !!number of boxcar pts
WRITE(iIOUN1) SNGL(kaFineFrStep(iTag))!!fine freq step size
WRITE(iIOUN1) i1,i2,iNLTEChunks       !!which chunks we are dealing with
!!and how many to expect
WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks
WRITE(iIOUN1) kaBlSize(iTag)/kBoxCarUse
!!fine 10000 point freq block size
WRITE(iIOUN1) iTotalStuff*kBoxCarUse  !!fine number of outputs

RETURN
END SUBROUTINE OpenBloatUAFile

!************************************************************************
! this subroutine writes out the main header for each results section
! this is a major change from before
! note the program does not keep on opening and closing the file

SUBROUTINE wrtout_head_uafile(caOutUAFile, rFrLow,rFrHigh,raFreq,iTag,  &
    iPathORRad,iNumberNLTEOut)


NO TYPE, INTENT(IN OUT)                  :: caOutUAFil
REAL, INTENT(OUT)                        :: rFrLow
REAL, INTENT(IN OUT)                     :: rFrHigh
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(OUT)                     :: iTag
INTEGER, INTENT(IN OUT)                  :: iPathORRad
NO TYPE, INTENT(IN OUT)                  :: iNumberNLT
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! caOutUAFile  = binary output file name
! rFrLow     = lowest wavenumber to be output
! rFrHigh    = highest wavenumber to be output
! rDelta     = point spacing
! iMainType  = 1 for path spectra                   iAtmNumber for Jacobians
!              2 for MP spectra                     iAtmNumber for fluxes
!              3 for radiances
! iSubMainType = kLayer2Sp (-2,-1,1,2,3,4) for path     GAS ID (1..63) for GasJac
!              = kLayer2Sp (-2,-1,1,2,3,4) for MP       0           for Temp Jac
!              = iAtmNumber for radiances             -10         for wgt fcn
!                                                     -20         for surf Jac
!              = +1 for upward flux, -1 for downward flux
! iNumberOut   = number of the relevant spectra to look for
! iPathOrRad   = +1 for CO2 ua path, +3 for ua rad

CHARACTER (LEN=80) :: caOutUAFile

INTEGER :: iNumberNLTEOut

REAL :: rDelta
INTEGER :: iMainType,iSubMainType,iIOUN

rDelta = kaFrStep(iTag)
IF (ABS((rFrHigh-rFrLow)-(rDelta*kMaxPts)) >= 2*rDelta) THEN
  WRITE(kStdErr,*) 'Wow! need (rFrHigh-rFrLow)-(rDelta*kMaxPts)) < rDelta'
  WRITE(kStdErr,*) 'rFrHigh,rFrLow,rDelta,kMaxPts,iTag = '
  WRITE(kStdErr,*) rFrHigh,rFrLow,rDelta,kMaxPts,iTag
  CALL DoSTOP
END IF
WRITE(kStdWarn,*) 'dump out : ',iNumberNLTEOut,' for iMain=',iMainType

iIOUN        = kNLTEOutUA

IF (iPathORRad == +1) THEN
  iMainType    = 1          !! right now these are CO2 opt depths
  iSubMainType = 1          !! right now assume only 1 atm
ELSE IF (iPathORRad == +3) THEN
  iMainType    = 3          !! right now these are rads
  iSubMainType = 1          !! right now assume only 1 atm
END IF

IF (kLongOrShort /= 0) THEN
  WRITE(iIOUN) iMainType,iSubMainType,iNumberNLTEOut
  WRITE(iIOUN) kMaxPts,rFrLow,rFrHigh,rDelta
END IF

RETURN
END SUBROUTINE wrtout_head_uafile
!************************************************************************
