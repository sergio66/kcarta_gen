c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c***********************************************************************
c**************** MAIN BINARY OUTPUT SUBROUTINES ARE HERE **************
c***********************************************************************
c this subroutine writes out the main header for each results section
c this is a major change from before
c note the program does not keep on opening and closing the file
      SUBROUTINE wrtout_head(iIOUN,caOutName,rFrLow,rFrHigh,rDelta,
     $                       iMainType,iSubMainType,iNumberOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iIOUN      = unit file number
c caOutName  = binary output file name
c rFrLow     = lowest wavenumber to be output   
c rFrHigh    = highest wavenumber to be output  
c rDelta     = point spacing
c iMainType  = 1 for path spectra                   iAtmNumber for Jacobians
c              2 for MP spectra                     iAtmNumber for fluxes
c              3 for radiances
c iSubMainType = kLayer2Sp (-2,-1,1,2) for path     GAS ID (1..63) for GasJac
c              = kLayer2Sp (-2,-1,1,2) for MP       0           for Temp Jac
c              = iAtmNumber for radiances           -10         for wgt fcn
c                                                   -20         for surf Jac
c              = +1 for upward flux, -1 for downward flux
c iNumberOut   = number of the relevant spectra to look for

      CHARACTER*80 caOutName
      REAL rFrLow,rFrHigh,rDelta
      INTEGER iMainType,iSubMainType,iNumberOut,iIOUN

      IF (abs((rFrHigh-rFrLow)-(rDelta*kMaxPts)).GE. 2*rDelta) THEN
         write(kStdErr,*) 'Wow! need (rFrHigh-rFrLow)-(rDelta*kMaxPts))
     $ < 2*rDelta'
         write(kStdErr,*) 'rFrHigh,rFrLow,rDelta,kMaxPts = '
         write(kStdErr,*) rFrHigh,rFrLow,rDelta,kMaxPts
         CALL DoSTOP
       END IF
      write(kStdWarn,*) 'dump out : ',iNumberOut,' for iMain=',iMainType

      IF (kLongOrShort .NE. 0) THEN
        WRITE(iIOUN) iMainType,iSubMainType,iNumberOut
        WRITE(iIOUN) kMaxPts,rFrLow,rFrHigh,rDelta
      END IF
      
      RETURN
      END
c************************************************************************
c this subroutine writes out the results
c this is a major change from before
c note the program does not keep on opening and closing the file
      SUBROUTINE wrtout(iIOUN,caOutName,raFreq,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iIOUN      = unit file number
c raInten    = array containing computed data (radiances,spectra,jacobians ..)
c raFreq    = array containing wavenumbers
c caOutName  = binary output file name

      INTEGER iIOUN
      REAL raFreq(kMaxPts),raInten(kMaxPts)
      CHARACTER*80 caOutName

c local variables
      INTEGER iInt

      WRITE(iIOUN) (raInten(iInt),iInt=1,kMaxPts)

      RETURN
      END

c************************************************************************
c     these are to dump out Planck multipliers (for NLTE)
c************************************************************************
c this subroutine dumps out the Planck Modifiers for UA
      SUBROUTINE  DumpPlanckUA(iAtm,iUpper,caPlanckFile,
     $                       raFreq,rDelta,raaUpperPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iAtm,iUpper
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      REAL raFreq(kMaxPts),rDelta !!array and spacing
      CHARACTER*80 caPlanckFile

c local variables
      INTEGER iFr,iL,iLay,iIOUN,iBeta
      REAL raTemp(kMaxPts)

      iIOUN = kStdPlanck
      CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),
     $                 rDelta,iAtm,1,iUpper)

      write(kStdWarn,*) 'subroutine UADumpPlanck is dumping out UAPlanck modifiers'
      write(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iUpper,' UA layers'

      DO iL = 1,iUpper
        iLay = iL
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaUpperPlanckCoeff(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine dumps out the Planck Modifiers
      SUBROUTINE  DumpPlanck(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,
     $                       raFreq,rDelta,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iAtm
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raFreq(kMaxPts),rDelta !!array and spacing
      CHARACTER*80 caPlanckFile

c local variables
      INTEGER iFr,iL,iLay,iIOUN,iBeta
      REAL raTemp(kMaxPts)

      iIOUN = kStdPlanck
      CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),
     $                 rDelta,iAtm,1,iaNumLayer(iAtm))

      write(kStdWarn,*) 'subroutine DumpPlanck is dumping out Planck modifiers'
      write(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iaNumLayer(iAtm),' layers'

      DO iL = 1,iaNumLayer(iAtm)
        iLay = iaaRadLayer(iAtm,iL)
        iBeta = MOD(iLay,kProfLayer) 
        IF (iBeta .EQ. 0) THEN 
          iBeta = kProfLayer 
        END IF 
        iLay = iBeta
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaPlanckCoeff(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine dumps out the Planck Modifiers, set at 1.0 (ie no NLTE yet
c in this chunk!!!!)
      SUBROUTINE  DumpPlanckOne(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,
     $                          raFreq,rDelta,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iAtm
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raFreq(kMaxPts),rDelta !!array and spacing
      CHARACTER*80 caPlanckFile

c local variables
      INTEGER iFr,iL,iLay,iIOUN,iBeta
      REAL raTemp(kMaxPts)

      write(kStdWarn,*) 'subroutine DumpPlanckOne is dumping out Planck modifiers = 1'
      write(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iaNumLayer(iAtm),' layers'

      iIOUN = kStdPlanck
      CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),
     $                 rDelta,iAtm,1,iaNumLayer(iAtm))

      DO iL = 1,iaNumLayer(iAtm)
        iLay = iaaRadLayer(iAtm,iL)
        iBeta = MOD(iLay,kProfLayer) 
        IF (iBeta .EQ. 0) THEN 
          iBeta = kProfLayer 
        END IF 
        iLay = iBeta
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 1.0
        END DO
        CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine dumps out the Planck Modifiers for UA as ones
      SUBROUTINE  DumpPlanckUAOne(iAtm,iUpper,caPlanckFile,
     $                       raFreq,rDelta,raaUpperPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iAtm,iUpper
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      REAL raFreq(kMaxPts),rDelta !!array and spacing
      CHARACTER*80 caPlanckFile

c local variables
      INTEGER iFr,iL,iLay,iIOUN,iBeta
      REAL raTemp(kMaxPts)

      iIOUN = kStdPlanck
      CALL wrtout_head(iIOUN,caPlanckFile,raFreq(1),raFreq(kMaxPts),
     $                 rDelta,iAtm,1,iUpper)

      write(kStdWarn,*) 'subroutine UADumpPlanckOne is dumping out UAPlanck modifiers = 1'
      write(kStdWarn,*) '  for Atm # ',iAtm,' which has ',iUpper,' UA layers'

      DO iL = 1,iUpper
        iLay = iL
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 1.0
        END DO
        CALL wrtout(iIOUN,caPlanckFile,raFreq,raTemp)
      END DO

      RETURN
      END

c************************************************************************
c      these are optical depths (individual gas or cumulative ODs)
c************************************************************************

c this function checks the appropriate set of paths to be output, to make sure
c that each path/MP/layer is being output only once! (else the bookkkeeping 
c in readatmos.m becomes messed up)
      INTEGER FUNCTION CheckDoubleEntry(iaaOp,iRow,iNumLay)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c iaaOp       = list of paths/MP/layers to be output for each print option
c iRow        = which row to be checked
c iNumLay     = number of elements to be checked in the iRow th row of iaaOp
      INTEGER iaaOp(kMaxPrint,kPathsOut),iRow,iNumLay

      INTEGER iAns,iJ

c assume everything OK
      iAns=1

c remember the entries in the row have already been sorted in ascending order,
c and so all that has to be done is to check that adjacent elements are 
c different from each other
      DO iJ=1,iNumLay-1
        IF (iaaOp(iRow,iJ) .EQ. iaaOp(iRow,iJ+1)) THEN
          write(kStdErr,*)'Outtype',iRow,' has',iaaOp(iRow,iJ)
          write(kStdErr,*)'entered more than once'
          iAns=-1
        END IF
      END DO

      CheckDoubleEntry=iAns
      
      RETURN
      END

c************************************************************************
c this function goes thru the list of layers to be output iaaOp(iJ,1..iEnd)
c to see if they do exist in atmosphere profile #iI, which is in
c iaaRadLayer(iI,1..iNL)
      INTEGER FUNCTION OutputLayerInProfile(iaaRadLayer,iI,iNL,
     $                                       iaaOp,iJ,iEnd)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iaaOp       = list of layers to be output
c iJ          = print output option number being considered
c iEnd        = numbver of layers in print option #iJ
c iaaRadLayer = for each atmosphere, list of layers used to buid up the atm
c iI          = atmosphere being considered
c iNL         = number of layers in the iI th atmosphere
      INTEGER iaaOp(kMaxPrint,kPathsOut)
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iI,iNL,iJ,iEnd

c local variables
      INTEGER iaTempRad(kMixFilRows)
      INTEGER iOK,iK
      INTEGER DoOutputLayer

c naively assume everything OK
      iOK=1

c now actually check that iaaOp(iJ,:) is in iaaRadLayer(iI,1..NL)
c first store the relevant iI'th atmosphere of iaaRadLayer(iI,1..NL)
      DO iK=1,iNL
        iaTempRad(iK)=iaaRadLayer(iI,iK)
      END DO
c and then sort it
      CALL DoSort(iaTempRad,iNL)

c now check that the elements of iaaOp to be output, iaaOp(iJ,1..iEnd) are
c all in iaTempRad
      iK=0
 20   iK=iK+1
      IF ((iK .LE. iEnd) .AND. (iOk .GT. 0)) THEN
        iOk=DoOutputLayer(iaaOp(iJ,iK),iNL,iaTempRad)
        IF (iOK .LT. 0) THEN
          write(kStdErr,*) 'output layer#',iaaOp(iJ,iK),' not found 
     $ in atm#',iI
        END IF
        GO TO 20
      END IF        

      OutputLayerInProfile=iOK

      RETURN
      END

c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now output 
c the individual GAS PATH spectra, according to what kLayer2Sp is set to
c kLayer2Sp = 2  : gas transmittances Layer to space sum(j=i,N) exp(-k(j))
c kLayer2Sp = 1  : gas optical depth Layer to space  sum(j=i,N) (k(j))
c kLayer2Sp = -1 : gas optical depth                 k(i)
c kLayer2Sp = -2 : gas transmittances                exp(-k(j))
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_trans_path(raFreq,rFrLow,rFrHigh,
     $           raaGasAbs,iPrinter,
     $           raTAmt,raTTemp,raTPress,raTPartPress,
     $           caOutName,
     $           iFileID,
     $           iaPath,iNp,iaOp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaGasAbs  = single gas abs coeffs
c iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
c iaPath     = list of the paths corresponding to the current gas 
c iNp        = total number of paths to be output
c iaOp       = list of the paths to be output
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaGasAbs(kMaxPts,kProfLayer)
      INTEGER iPrinter,iFileID
      INTEGER iNp,iaOp(kPathsOut),iaPath(kProfLayer)
      CHARACTER*80 caOutName
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer),raTPress(kProfLayer),raTPartPress(kProfLayer)

c local variables
      INTEGER iInt,iDiv,iDp,iStart,iPath,iLay,DoOutputLayer,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN = kStdkCarta

      iStart=iDiv(iaPath(1),kProfLayer)

c write spectra to unformatted file
c if iPrinter=1 then have to check for valid paths
      DO iLay=1,kProfLayer
c check to see if this path should be output
        iPath=iStart*kProfLayer + iLay
c        IF (iPath .NE. iaPath(iLay)) THEN
c          write(kStdErr,*) 'iPath .NE. iaPath(iLay)' 
c          Call DoSTOP
c        END IF
        iDp=DoOutputLayer(iPath,iNp,iaOp)
        IF (iDp .GT. 0) THEN
          IF (kLayer2Sp .EQ. -1) THEN
            write(kStdWarn,*)'output GAS OD : iPath,P,T,A,OD = ',
     $ iPath,raTPress(iLay)*kAtm2mb,raTTemp(iLay),raTAmt(iLay),raaGasAbs(1,iLay)
            DO iInt=1,kMaxPts
              raL2S(iInt)=raaGasAbs(iInt,iLay)
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
          ELSE IF (kLayer2Sp .EQ. -2) THEN
            write(kStdWarn,*)'outputting GAS layer trans at iPath = ',iPath
            DO iInt=1,kMaxPts
              raL2S(iInt)=exp(-raaGasAbs(iInt,iLay))
            END DO
            CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
          ELSE IF (kLayer2Sp .EQ. 1) THEN
            write(kStdWarn,*)'outputting GAS layer2sp OD at iPath = ',iPath
            CALL GasOptLayer2Space(raaGasAbs,raL2S,iLay)
            CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
          ELSE IF (kLayer2Sp .EQ. 2) THEN
            write(kStdWarn,*)'outputting GAS layer2sp trans at iPath = ',iPath
            CALL GasTranLayer2Space(raaGasAbs,raL2S,iLay)
            CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
          END IF  
        END IF
      END DO

c      WRITE (6,1001)caOutName
c 1001 FORMAT('Successfully saved unformatted PATH results to ',/,A80)
     
      RETURN
      END

c************************************************************************
c this is the generic call to MP output
      SUBROUTINE DoOutputMixedPaths(
     $                       iFound,iPrinter,caOutName,
     $                       raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,
     $                       iNpmix,iFileID,iNp,iaOp)

      IMPLICIT NONE

      INCLUDE '../INCLUDE/kcarta.param'

c input params
      INTEGER iFound,iPrinter,iNpmix,iFileID,iNp,iaOp(kPathsOut)
      REAL raFreq(kMaxPts),rFreqStart,rFreqEnd
      REAL raaSumAbCoeff(kMaxPts,kMixFilRows) 
      CHARACTER*80 caOutName

c local vars 
      INTEGER iDummy,iIpmix,DoOutputLayer

      IF ((iFound .GT. 0) .AND. (iPrinter .EQ. 2)) THEN

c we have a list of mixed paths to output!!!
        IF (iNp .LT. 0) THEN
          write(kStdWarn,*) 'OUTPUTTING ALL MIXED PATHS ... '
        END IF

        IF ((kLayer2Sp .EQ. 2) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
          write(kStdWarn,*)'     outputting ALL L2S transmittances'    
          CALL out_FASTL2Strans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
        ELSE IF ((kLayer2Sp .EQ. 1) .AND. (iNp .LT. 0)) THEN
          write(kStdWarn,*) '     outputting ALL L2S optical depths'    
c this indicates we want L2S transmittances for ALL mixed paths!!!
          CALL out_FASTL2Soptdp_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
        ELSE IF ((kLayer2Sp .EQ. -1) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
          write(kStdWarn,*) '    outputting ALL layer optical depths'    
          CALL out_FASToptdp_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
        ELSE IF ((kLayer2Sp .EQ. -2) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
          write(kStdWarn,*) '    outputting ALL layer transmittances'
          CALL out_FASTtrans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
        ELSE
c now loop over the mixed paths, outputting whichever ones are necessary
          DO iIpmix = 1,iNpmix
c see if mixed path iIpmix is set from *OUTPUT
            iDummy = -1
            iDummy = DoOutputLayer(iIpmix,iNp,iaOp)
c if the printing option=2 and iIpmix has been found in the list of 
c paths to be output then go ahead and print the relevant (iIpmix th) row of
c SUM abs spectra
            IF ((iPrinter .EQ. 2) .AND. (iDummy .GT. 0)) THEN
              CALL out_trans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iIpmix,iNpmix,iFileID)
            END IF
c go to next MIXED PATH set by incrementing iIpmix
          END DO
        END IF
c end if (iFound==1)
      END IF

      RETURN
      END

c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH transmittance spectra
c kLayer2Sp = 2  : MP transmittances Layer to space sum(j=i,N) exp(-k(j))
c kLayer2Sp = 1  : MP optical depth Layer to space  sum(j=i,N) (k(j))
c kLayer2Sp = -1 : MP optical depth                 k(i
c kLayer2Sp = -2 : MP transmittances                exp(-k(j)
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_trans_MP(raFreq,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iIpmix,iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iIpmix     = current mixed path number
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iIpmix,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iInt,iPath,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN = kStdkCarta

c write spectra to unformatted file
c this is for mixed paths
      iPath=iIpmix
      write(kStdWarn,*)'output MIXED path optical depths at iIpmix = ',iIpmix
      IF (kLayer2Sp .EQ. -1) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt)=raaSumAbs(iInt,iIpmix)
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      ELSE IF (kLayer2Sp .EQ. -2) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt)=exp(-raaSumAbs(iInt,iIpmix))
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      ELSE IF (kLayer2Sp .EQ. 1) THEN
        CALL MP2SpaceOptDp(raaSumAbs,raL2S,iIpmix,iNpmix)
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      ELSE IF (kLayer2Sp .EQ. 2) THEN
        CALL MP2SpaceTrans(raaSumAbs,raL2S,iIpmix,iNpmix)
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      END IF  

c      WRITE (6,1001)caOutName
c 1001 FORMAT('Successfully saved unformatted MP results to ',/,A80)
     
      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH L2Stransmittance spectra FAST. This assumes each atmos
c is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
c has less than 100 layers (i.e. iNpmix is not a multiple of 100),
c then the "upper layers" are filled with spectra of 0, so their 
c transmittance is 1

c this assumes kLayer2Sp = 2, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASTL2Strans_MP(raFreq,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iPath,iIOUN,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax,iI
      REAL raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

      iIOUN = kStdkCarta

c first figure out how many 100 atmosphere layer blocks there are
      iAtmBlocks=1
 60   CONTINUE
      IF (kProfLayer*iAtmBlocks .LT. iNpMix) THEN
        iAtmBlocks=iAtmBlocks+1
        GO TO 60
      END IF
      iAmax=iAtmBlocks
        
c make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
c saying the last "atmosphere" has less than set number of layers
      iWarn=MOD(iNpmix,kProfLayer)
      IF (iWarn .NE. 0) THEN
        iAmax=iAmax-1
      END IF
        
c calculate L2S and write spectra to unformatted file
      iPath=0

c first do the "complete" atmospheres
      DO iA=1,iAmax
c put the current atmosphere into raaTempArray
        DO iLay=1,kProfLayer
          iI=iLay+(iA-1)*kProfLayer
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
          END DO
        END DO
c compute the L2S optical depths
        DO iLay = kProfLayer-1,1,-1
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+
     $                             raaTempArray(iFr,iLay+1)
          END DO
        END DO
c print them out
        DO iLay=1,kProfLayer
          iPath=iPath+1
          write(kStdWarn,*) 'output MIXED path spectra at ',iPath
          DO iFr=1,kMaxPts
            raL2S(iFr)=exp(-raaTempArray(iFr,iLay))
          END DO
          CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
        END DO
      END DO

c then do the last, "incomplete" atmospheres
c put the current atmosphere into raaTempArray
      IF (iWarn .NE. 0) THEN
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
c compute the L2S optical depth
        DO iLay=iWarn-1,1,-1
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+
     $                               raaTempArray(iFr,iLay+1)
          END DO
        END DO
c print them out
        DO iLay=1,iWarn
          iPath=iPath+1
          write(kStdWarn,*) 'output MIXED path spectra at ',iPath
          DO iFr=1,kMaxPts
            raL2S(iFr)=exp(-raaTempArray(iFr,iLay))
          END DO
          CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
        END DO
      END IF

      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH L2S optical depth FAST. This assumes that in each atmos
c is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
c is found to have less than 100 layers (i.e. iNpmix is not a multiple of 100),
c then the "upper layers" are filled with spectr of 0, so their transmittance
c is 1

c this assumes kLayer2Sp = 1, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASTL2Soptdp_MP(raFreq,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iPath,iIOUN,iI,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax
      REAL raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

      iIOUN = kStdkCarta

c first figure out how many kProfLayer atmosphere layer blocks there are
      iAtmBlocks=1
 60   CONTINUE
      IF (kProfLayer*iAtmBlocks .LT. iNpMix) THEN
        iAtmBlocks=iAtmBlocks+1
        GO TO 60
      END IF
      iAmax=iAtmBlocks
        
c make sure that iNpmix can be exactly divided into kProfLayer, else set a flag
c saying the last "atmosphere" has less than kProfLayer layers
      iWarn=MOD(iNpmix,kProfLayer)
      IF (iWarn .NE. 0) THEN
        iAmax=iAmax-1
      END IF
        
c calculate L2S and write spectra to unformatted file
      iPath=0

c first do the "complete" atmospheres
      DO iA=1,iAmax
c put the current atmosphere into raaTempArray
        DO iLay=1,kProfLayer
          iI=iLay+(iA-1)*kProfLayer
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaSumAbs(iFr,iI)
          END DO
        END DO
c compute the L2S optical depths
        DO iLay = kProfLayer-1,1,-1
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+
     $                             raaTempArray(iFr,iLay+1)
          END DO
        END DO
c print them out
        DO iLay=1,kProfLayer
          iPath=iPath+1
          DO iFr=1,kMaxPts
            raL2S(iFr)=raaTempArray(iFr,iLay)
          END DO
          write(kStdWarn,*) 'output MIXED path spectra  at ',iPath
          CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
        END DO
      END DO

c then do the last, "incomplete" atmospheres
c put the current atmosphere into raaTempArray
      IF (iWarn .NE. 0) THEN
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
c compute the L2S optical depths
        DO iLay=iWarn-1,1,-1
          DO iFr=1,kMaxPts
            raaTempArray(iFr,iLay)=raaTempArray(iFr,iLay)+
     $                               raaTempArray(iFr,iLay+1)
          END DO
        END DO
c print them out
        DO iLay=1,iWarn
          iPath=iPath+1
          write(kStdWarn,*) 'output MIXED path spectra at ',iPath
          DO iFr=1,kMaxPts
            raL2S(iFr)=raaTempArray(iFr,iLay)
          END DO
          CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
        END DO
      END IF

      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH optical depth FAST. 

c this assumes kLayer2Sp = -1, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASToptdp_MP(raFreq,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iIOUN,iI,iFr
      REAL raL2S(kMaxPts)

      iIOUN = kStdkCarta

      DO iI=1,iNpMix
        DO iFr=1,kMaxPts
          raL2S(iFr)=raaSumAbs(iFr,iI)
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      END DO

      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH transmittance spectra FAST. 

c this assumes kLayer2Sp = -2, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASTtrans_MP(raFreq,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raFreq(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iI,iFr,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN = kStdkCarta

      DO iI=1,iNpMix
        DO iFr=1,kMaxPts
          raL2S(iFr)=exp(-raaSumAbs(iFr,iI))
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raL2S)
      END DO

      RETURN
      END

c************************************************************************
c this subroutine does the layer to space optical depths by doing
c the calculation between frequency indices from layer iLay-->100
c the result is stored in raL2S
c this is mainly used by the OpticalLayer2Space printout routines in rdprof.f
      SUBROUTINE GasOptLayer2Space(raaGasAbs,raL2S,iLay)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLay      = current layer number
c raaGasAbs = matrix containing the raw abs coeffs
c raL2S     = array containing the layer-to-space results 
      INTEGER iLay
      REAL raL2S(kMaxPts),raaGasAbs(kMaxPts,kProfLayer)

      INTEGER iI,iJ

      DO iI=1,kMaxPts
        raL2S(iI)=0.0
      END DO

      DO iJ=iLay,kProfLayer
        DO iI=1,kMaxPts
          raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
        END DO
      END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the layer to space transmittances for gas iG by doing
c the calculation between frequency indices from layer iLay-->100
c the result is stored in raL2S
c this is mainly used by the OpticalLayer2Space printout routines in rdprof.f
      SUBROUTINE GasTranLayer2Space(raaGasAbs,raL2S,iLay)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLay      = current layer number
c raaGasAbs = matrix containing the raw abs coeffs
c raL2S     = array containing the layer-to-space results 
      INTEGER iLay
      REAL raL2S(kMaxPts),raaGasAbs(kMaxPts,kProfLayer)

      INTEGER iI,iJ

      DO iI=1,kMaxPts
        raL2S(iI)=0.0
      END DO

      DO iJ=iLay,kProfLayer
        DO iI=1,kMaxPts
          raL2S(iI)=raL2S(iI)+raaGasAbs(iI,iJ)
        END DO
      END DO

      DO iI=1,kMaxPts
        raL2S(iI)=exp(-raL2S(iI))
      END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the MP to space transmission coeffs by doing
c the calculation between frequency indices from layer iLay-->iLM
c where iLm=nearest 100 index above iLay, or iNpmix 
c the result is stored in raL2S
      SUBROUTINE MP2SpaceTrans(raaGasAbs,raL2S,iIpmix,iNpmix)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iIpMix    = current mixed path number
c iNpMix    = total number of mixed paths
c raaGasAbs = matrix containing the raw mixed path abs coeffs
c raL2S     = array containing the layer-to-space results 
      INTEGER iIpmix,iNpmix
      REAL raL2S(kMaxPts),raaGasAbs(kMaxPts,kMixFilRows)

c local variables
      INTEGER iI,iJ,iMax,iFound,iUpper,idiv

c first find the index where the (space) "upper" limit is at
      iFound=-1
      iI=0
      iJ=1
      iMax=idiv(iNpmix,kProfLayer)
      iUpper=MOD(iNpmix,kProfLayer)
      IF (iUpper .NE. 0) THEN
        iMax=iMax+1
      END IF

 15   CONTINUE
      IF ((iFound .LT. 0) .AND. (iJ .LE. iMax)) THEN
        IF ((kProfLayer*iI .LE. iIpmix) .AND. 
     $      (kProfLayer*iJ .GE. iIpmix)) THEN
          iFound=1
          iUpper = kProfLayer*iJ
          IF (iUpper .GT. iNpmix) THEN
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
        raL2S(iI)=exp(-raL2S(iI))
      END DO

      RETURN
      END 

c************************************************************************
c this subroutine does the MP to space optical depths by doing
c the calculation between frequency indices from layer iLay-->iLM
c where iLm=nearest 100 index above iLay, or iNpmix 
c the result is stored in raL2S
      SUBROUTINE MP2SpaceOptDp(raaGasAbs,raL2S,iIpmix,iNpmix)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iIpMix    = current mixed path number
c iNpMix    = total number of mixed paths
c raaGasAbs = matrix containing the raw mixed path abs coeffs
c raL2S     = array containing the layer-to-space results 
      INTEGER iIpmix,iNpmix
      REAL raL2S(kMaxPts),raaGasAbs(kMaxPts,kMixFilRows)

c local variables
      INTEGER iI,iJ,iMax,iFound,iUpper,idiv

c first find the index where the (space) "upper" limit is at
      iFound=-1
      iI=0
      iJ=1
      iMax=idiv(iNpmix,kProfLayer)
      iUpper=MOD(iNpmix,kProfLayer)
      IF (iUpper .NE. 0) THEN
        iMax=iMax+1
      END IF

 15   CONTINUE
      IF ((iFound .LT. 0) .AND. (iJ .LE. iMax)) THEN
        IF ((kProfLayer*iI .LE. iIpmix) .AND. 
     $      (kProfLayer*iJ .GE. iIpmix)) THEN
          iFound=1
          iUpper = kProfLayer*iJ
          IF (iUpper .GT. iNpmix) THEN
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
      END 

c************************************************************************
c this subroutine appends _VT_iRegr to end of caOutName to get caVTName
c this is when we use regression to try to predict the best NLTE temps
      SUBROUTINE VTName_regr(iRegr,caOutName,caVTFile)

      IMPLICIT NONE

      CHARACTER*80 caOutName,caVTFile
      INTEGER iRegr
 
      CHARACTER*2 caString
      INTEGER iInt,iInt1

      DO iInt=1,80 
        caVTFile(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caOutName(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caVTFile(1:iInt)=caOutName(1:iInt)
      caVTFile(iInt+1:iInt+3)='_VT'

c now process iRegr so that we end up with a right padded string 
c eg 2 ---> '2 ', 12 ---> '12' etc 
      WRITE(caString,15) iRegr
 15   FORMAT(I2) 
c this is right justified ... change to left justified 
      iInt=1 
 16   continue 
      IF (caString(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 16 
      END IF 

      iInt1 = 80 
 21   CONTINUE 
      IF ((caVTFile(iInt1:iInt1) .EQ. ' ') .AND. (iInt1 .GE. 1)) THEN 
        iInt1 = iInt1-1 
        GO TO 21 
      END IF 
      iInt1 = iInt1 + 1

      caVTFile(iInt1:iInt1+(2-iInt))=caString(iInt:2) 
      
      RETURN
      END

c************************************************************************
c this subroutine appends _RTP_irtp to end of caOutName to get caVTName
      SUBROUTINE VTName_rtp(iRtp,caOutName,caVTFile)

      IMPLICIT NONE

      CHARACTER*80 caOutName,caVTFile
      INTEGER iRtp
 
      CHARACTER*5 caString
      INTEGER iInt,iInt1

      DO iInt=1,80 
        caVTFile(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caOutName(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caVTFile(1:iInt)=caOutName(1:iInt)
      caVTFile(iInt+1:iInt+5)='_RTP_'

c now process iRtp so that we end up with a right padded string 
c eg 2 ---> '2 ', 12 ---> '12' etc 
      WRITE(caString,15) iRtp
 15   FORMAT(I5) 
c this is right justified ... change to left justified 
      iInt=1 
 16   continue 
      IF (caString(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 16 
      END IF 

      iInt1 = 80 
 21   CONTINUE 
      IF ((caVTFile(iInt1:iInt1) .EQ. ' ') .AND. (iInt1 .GE. 1)) THEN 
        iInt1 = iInt1-1 
        GO TO 21 
      END IF 
      iInt1 = iInt1 + 1

      caVTFile(iInt1:iInt1+(5-iInt))=caString(iInt:5) 
      
      RETURN
      END

c************************************************************************
c this subroutine appends _VT_iRTP to end of caOutName to get caVTName
c so iRegr is NOT used!
c this is when we use a polynom fit, or a predictor fit, to estimate NLTE temps
c this bloody thing DOES NOT WORK!!!!
c this bloody thing DOES NOT WORK!!!!
c this bloody thing DOES NOT WORK!!!!
      SUBROUTINE VTName(iRTP,caOutName,caVTFile)

      IMPLICIT NONE

c input params
      CHARACTER*80 caOutName
      INTEGER iRTP
c output params
      CHARACTER*80 caVTFile
 
      CHARACTER*5 caString
      INTEGER iInt,iInt1,iI,iJ

      DO iInt=1,80 
        caVTFile(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caOutName(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caVTFile(1:iInt)        = caOutName(1:iInt)
      caVTFile(iInt+1:iInt+5) = '_RTP_'

c now process iRTP so that we end up with a right padded string 
c eg 2 ---> '2 ', 12 ---> '12' etc 
c assume at most there are 99999 profiles in one RTP file
      caString = '     '
      WRITE(caString,15) iRTP
 15   FORMAT(I5) 
c this is right justified ... change to left justified 
      iInt=1 
 16   continue 
      IF (caString(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 16 
      END IF 

      iInt1 = 80 
 21   CONTINUE 
      IF ((caVTFile(iInt1:iInt1) .EQ. ' ') .AND. (iInt1 .GE. 1)) THEN 
        iInt1 = iInt1-1
        GO TO 21 
      END IF 
      iInt1 = iInt1 + 1

c      caVTFile(iInt1:iInt1+(5-iInt)) = caString(iInt:5) 
      iJ = 0
      DO iI = iInt,5
        caVTFile(iInt1+iJ:iInt1+iJ) = caString(iI:iI) 
        iJ = iJ + 1
      END DO

      RETURN
      END

c************************************************************************
c this subroutine appends _COL to end of caOutName to get "new" caJacobFileName
      SUBROUTINE jacob2Name(caJacobFile2,caJacobFile)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 caJacobFile,caJacobFile2
 
      INTEGER iInt

      DO iInt=1,80 
        caJacobFile2(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caJacobFile(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caJacobFile2(1:iInt) = caJacobFile(1:iInt)
      caJacobFile2(iInt+1:iInt+4) = '_COL'

      RETURN
      END

c************************************************************************
c this subroutine appends _FLUX to end of caOutName to get caFluxName
      SUBROUTINE FluxName(caFluxFile,caOutName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 caOutName,caFluxFile
 
      INTEGER iInt

      DO iInt=1,80 
        caFluxFile(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caOutName(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caFluxFile(1:iInt)=caOutName(1:iInt)

c before Jun 2013
c      IF (kFLux .LE. 2) THEN
c        caFluxFile(iInt+1:iInt+5)='_FLUX'
c      ELSEIF (kFLux .LE. 5) THEN
c        caFluxFile(iInt+1:iInt+4)='_OLR'
c      ELSE
c        write(kStdErr,*) 'Unknown option for kFlux = ',kFlux
c        Call DoStop
c      END IF

c after Jun 2013
      write(kStdWarn,*) 'kFlux = ',kFlux
      IF (kFLux .EQ. 1) THEN
        caFluxFile(iInt+1:iInt+5)='_DOWN'   !! down welling flux at bottom of 100 layers
      ELSEIF (kFLux .EQ. 2) THEN
        caFluxFile(iInt+1:iInt+5)='_HEAT'   !! up-down heating rate at all 100 layers
      ELSEIF (kFLux .EQ. 3) THEN
        caFluxFile(iInt+1:iInt+3)='_UP'     !! up welling flux at to[ of 100 layers
      ELSEIF (kFLux .EQ. 4) THEN
        caFluxFile(iInt+1:iInt+4)='_OLR'    !! upwelling flux at TOA
      ELSEIF (kFLux .EQ. 5) THEN
        caFluxFile(iInt+1:iInt+5)='_OLR3'   !! upwelling flux at TOA, tropopause,dnwell at gnd
      ELSEIF (kFLux .EQ. 6) THEN
        caFluxFile(iInt+1:iInt+4)='_ALL'    !! up,down welling flux at all tops/bottoms of each layer
      ELSE
        write(kStdErr,*) 'Unknown option for kFlux = ',kFlux
        Call DoStop
      END IF

      RETURN
      END

c************************************************************************
c this subroutine appends _PLANCK to end of caOutName to get caPlanckName
      SUBROUTINE PlanckName(caPlanckFile,caOutName)

      IMPLICIT NONE

      CHARACTER*80 caOutName,caPlanckFile
 
      INTEGER iInt

      DO iInt=1,80 
        caPlanckFile(iInt:iInt)=' ' 
      END DO 

      iInt=80 
 11   CONTINUE 
      IF ((caOutName(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
      END IF 

      caPlanckFile(1:iInt)=caOutName(1:iInt)
      caPlanckFile(iInt+1:iInt+7)='_PLANCK'
      
      RETURN
      END

c************************************************************************
c this subroutine writes out the header for the summary text file "header.head"
c this subroutine writes out the header for the binary output file ...
c any other writes to this file will be appended to the end

c as the header is written out, the subroutine checks that print option 1 has
c been specified not more than once, print option 3 specified not more than 
c once and print option 7 not more than once for the ith atmosphere

c kLongOrShort enables shorted output files to be written out (-1,+1), 
c or just the basic results file (0)
      SUBROUTINE PrepareOutput(caDriver,caOutName,caJacobFile,caJacobFile2,
     $      caFluxFile,caPlanckFile,iOutFileName,iNumNLTEGases,
     $      rFrLow,rFrHigh,iFileIDLo,iFileIDHi,caComment,
     $      iNumGases,iaGases,raaAmt,raaTemp,raaPress,raaPartPress,
     $                        raaRAmt,       raaRPartPress,
     $      raPressLevels,iProfileLayers,
     $      iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVT,
     $      iNatm,iNatm2,iaNumLayers,iaaRadLayer,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      raaaSetEmissivity,iaSetEms,
     $      iOutTypes,iaPrinter,iaAtmPr,iaNp,iaaOp,raaUserPress,
     $      iJacob,iaJacob,
     $      iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $      rakThermalAngle,iakThermalJacob,iaOutNumbers,
     $      iTotal,iTotalStuff,
     $      iDoUpperAtmNLTE,iDumpAllUASpectra,iDumpAllUARads)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      
c raPressLevels, iProfileLayers = actual number of pressure layers
c iOutFileName = does caOutName exist, or is stuff dumped to screen? for fluxes
c caComment   = comment that the user put into *RADNCE
c iJacob      = number of d/dq we will output (>=1)
c iaJacob     = list of GasID's whose  d/dq we will output
c caDriver     = root name of input file (not really needed here)
c caOutName   = name of output binary file, specified in *OUTPUT
c caDummyFile = name of output binary file, to save Dummy stuff
c rFrLow      = lower wavenumber bound
c rFrHigh     = upper wavenumber bound
c iFileIDLo      = from rFrLow, lower file index (1-kNumkFile)
c iFileIDHi      = from rFrHigh, upper file index (1-kNumkFile)
c iNumGases   = number of gases read in from *GASFIL + *XSCFIL
c iaGases     = integer array containig the GAS ID's
c raaAmt      = array containing amount profiles for each gas
c raaTemp     = array containing temp profiles for each gas
c raaPartPress= array containing partial pressure profiles for each gas
c iNpMix      = number of mixed paths
c raaMix      = complete mixing table
c caMixFName  = name of file containing character MP info (from *MIXFIL)
c iMixFileLines = number of lines containig character info for mixtable
c iNatm       = number of atmospheres read in from *RADIL
c iNatm2      = number of atmospheres read in from *OUTPUT
c iaNumLayers = number of layers in each atmosphere
c iaaRadLayer = for each atmosphere, a list of the layers making it up
c raTSpace    = for each atmosphere, the atmosphere backgnd temperature
c raTSurf     = for each atmosphere, the surface temperature
c raSatAngle  = for each atmosphere, the satellite view angle
c raSatHeight = for each atmosphere, the satellite height
c iOutTypes   = number of output options found in *OUTPUT
c iaPrinter   = for each output option, what to be printed (1,2 or 3)
c iaAtmPr     = number of layers to be printed for each atmosphere
c iaNp        = number of paths/MP/layers to be  output for each print option
c iaaOp       = list of paths/MP/layers to be output for each print option
c raaUserPress= list of output pressures for each print option
c iaSetEms    = number of emissivity data points per atmosphere 
c raaaSetEms  = emissivity data 
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c rakSolarRefl   =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaOutNumbers = how many of each print option there are to output
c iTotal = how many of the kcomp files are gonna be unchunked
c iDumpAllUARads = do we dump rads for all layers (-1) or a specific number?
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaOutNumbers(kMaxPrint),iOutFileName
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      CHARACTER*80 caDriver,caOutName,caComment
      CHARACTER*130 caaMixFileLines(kProfLayer)
      INTEGER iaPrinter(kMaxPrint),iaAtmPr(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iOutTypes,iMixFileLines
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iFileIDLo,iFileIDHi
      INTEGER iNumGases,iaGases(kMaxGas),iNpmix,iTotal
      REAL raMixVT(kMixFilRows),raaUserPress(kMaxPrint,kProfLayer)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raPressLevels(kProfLayer+1)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raaPartPress(kProfLayer,kGasStore)
      REAL raaRAmt(kProfLayer,kGasStore),raaRPartPress(kProfLayer,kGasStore)      
      REAL raaMix(kMixFilRows,kGasStore),rFrLow,rFrHigh
      INTEGER iNatm,iNatm2,iaNumLayers(kMaxAtm),iJacob,iaJacob(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iProfileLayers
      INTEGER iDoUpperAtmNLTE,iNumNLTEGases,iDumpAllUARads,iDumpAllUASpectra

      INTEGER iIOUN,iIOUN1,iIOUN2,iIOUN_JAC2,iI,iJ,iK,iFileErr,iEnd,iP,iOk
      INTEGER iOutputOptionNum,iNumLay,DoGasJacob
      CHARACTER*80 caJacobFile,caJacobFile2,caFluxFile,caPlanckFile

      INTEGER CheckDoubleEntry,iImportant
      INTEGER iNatmJac,iaLayerJac(kMaxAtm),iIOUN_Flux,iIOUN_Planck,iIOUN_Cloud
      REAL raParams(kMaxUserSet),raPActualAvg(kProfLayer),rP
      CHARACTER*4 caStrJunk(7)
      CHARACTER*80 caFCloudName
      REAL raSumTotalGasAmt(kMaxGas)
      
      !this is for kLongOrShort = 0
      INTEGER iTag,iTotalStuff      
      INTEGER iCo2,iaCO2path(kProfLayer),iaLayerFlux(kMaxAtm)
      INTEGER iaJunkFlux(2*kProfLayer)

      CHARACTER*120 caStr1,caStr2
      
      include '../INCLUDE/gasIDname.param'
      
      DO iI = 1, kMaxGas
        raSumTotalGasAmt(iI) = 0.0
      END DO
      
      iDumpAllUASpectra = -1 !!assume no need to dump out CO2 UA spectra
      IF (iDoUpperAtmNLTE .GT. 0) THEN
        !!need to see if we are dumping out CO2 paths
        iCO2 = -1
        DO iTag = 1,iNumGases
          IF (iaGases(iTag) .EQ. 2) THEN
            iCO2 = iTag
          END IF
        END DO
        IF (iCO2 .GT. 0) THEN
          DO iTag = 1,kProfLayer
            !!! these are the CO2 paths
            iaCO2path(iTag) = iTag + (iCO2-1)*kProfLayer
          END DO
          !!check to see if any of them are being output
          DO iI = 1,iOutTypes
            IF (iaPrinter(iI) .EQ. 1) THEN !!! these are paths to be output
              iEnd = iaNp(iI)
              DO iTag = 1,kProfLayer
                DO iP = 1,iEnd
                  IF (iaaOp(iI,iP) .EQ. iaCO2path(iTag)) THEN
                    iDumpAllUASpectra = +1
                  END IF
                END DO
              END DO
            END IF
          END DO
        ELSE
          write(kStdErr,*) 'oh oh, need to do NLTE and not using CO2????'
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
        raPActualAvg(iI) = rP/log(raPressLevels(iI)/raPresslevels(iI+1))
      END DO
      
      write(kStdWarn,*)'Preparing header info for the output files ..'

      DO iImportant=1,kMaxPrint
        iaOutNumbers(iImportant)=0
      END DO
       
 800  FORMAT(A1)
 801  FORMAT(A80)

      iIOUN =  kStdWarn
      iIOUN1 = kStdkCarta

      iTotalStuff = 0 !total num of spectra,mixed spectra, radiances to output

c set up raParams; recall all are integers, so change to real!!!
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

      write(kStdWarn,*) 'kCKD = ',kCKD

c check to make sure that the user first wanted to specify paths, then mixed 
c paths, then radiances; else STOP ie we have to find paPrinter(iI)=1,2,3 or
c in that order eg 1,1,1,1 or 2,3,3,3,3  but not 3,1,1 etc
      iP=-1
      DO iI=1,iOutTypes
        IF (iaPrinter(iI) .GE. iP) THEN
          iP = iaPrinter(iI) !ok; ia(printer(iI) is going UP or staying same
        ELSE !oh oh iaPrinter(iI) is coming down!!! readers cannot handle this
          write(kStdErr,*)'in *OUTPUT need to specify paths, then mixed'
          write(kStdErr,*)'   paths, and then radiances ie order is '
          write(kStdErr,*)'   important. Please rewrite *OUTPUT section'
          CALL DoSTOP
        END IF
      END DO

c if kLongOrShort = 1, output all info, and summarize in kStdWarn
c if kLongOrShort = 0, output ONLY DATA, but summarize in kStdWarn
c                      we need to figure out iTotalStuf ==> do this loop
c if kLongOrShort = -1, output shortened version of binary file

c first open summary text file as a fresh file to be written to
      IF (kLongOrShort .GE. 0) THEN

        WRITE(iIOUN,*) 'SUMMARY OF INPUT DRIVER FILE ... '
        WRITE(iIOUN,*) '     '
        WRITE(iIOUN,*) '     '

c first output general info -----------------------------------------
        WRITE(iIOUN,*) 'GENERAL INFORMATION : '
        WRITE(iIOUN,*) caVersion
        WRITE(iIOUN,*) 'Number of layers = ',kProfLayer
        WRITE(iIOUN,*) 'Number of parameters in *PARAMS = ',kMaxUserSet
        WRITE(iIOUN,*) (raParams(iI),iI=1,kMaxUserSet)
        WRITE(iIOUN,*) caComment
        WRITE(iIOUN,*) 'Freq endpts = ',rFrLow,rFrHigh
        WRITE(iIOUN,*) 'k-comp file endpts = ',iFileIDLo,iFileIDHi
        WRITE(iIOUN,*) 'Long/short output file = ',kLongOrShort
c v1.04 to v1.08 had these next 2 lines; now remove them
c        WRITE(iIOUN,*) 'M1000mb,M100mb,MSub = ',M1000mb,M100mb,MSubLayer
c        WRITE(iIOUN,*) 'M50mb,M10mb,MThick = ',M50mb,M10mb,MThickLayer
        WRITE(iIOUN,*) 'average layer pressures ...'
        WRITE(iIOUN,*) (raPActualAvg(iI),iI=1,kProfLayer)
        WRITE(iIOUN,*) '***********************************************'

c then output path ID stuff ------------------------------------------
        WRITE(iIOUN,*) 'PATH ID INFO'
        WRITE(iIOUN,*) 'iNumPaths = ',iNumGases*kProfLayer
        WRITE(iIOUN,*) 'Num Layers in Profile = ',iProfileLayers
        WRITE(iIOUN,*) 'So start showing info from layer ',kProfLayer-iProfileLayers+1

        caStr1 = '  Path# GID  Press      PartP        PPMV        Temp         Amnt    ||      RefPP     RefAmt  Amt/RAmt'
	caStr2 = '----------------------------------------------------------------------||--------------------------------'	  
        DO iI=1,iNumGases
	  write(iIOUN,234) iI,iaGases(iI),caGID(iaGases(iI))
          WRITE(iIOUN,7170) caStr1
          WRITE(iIOUN,7170) caStr2 
          DO iJ=kProfLayer-iProfileLayers+1,kProfLayer
            iP=(iI-1)*kProfLayer+iJ
	    raSumTotalGasAmt(iaGases(iI)) = raSumTotalGasAmt(iaGases(iI)) + raaAmt(iJ,iI)
            WRITE(iIOUN,7171)iP,iaGases(iI),raPActualAvg(iJ)/kAtm2mb,raaPartPress(iJ,iI),
     $                   raaPartPress(iJ,iI)/(raPActualAvg(iJ)/kAtm2mb)*1.0e6,raaTemp(iJ,iI),raaAmt(iJ,iI),'||',
     $                   raaRPartPress(iJ,iI),raaRAmt(iJ,iI),raaAmt(iJ,iI)/raaRAmt(iJ,iI)
c earlier code was giving raaRPartPress(iJ,iI)/raaPartPress(iJ,iI) = 385/360 for CO2, tape6 test from eli     
c            WRITE(iIOUN,7171)iP,iaGases(iI),raPActualAvg(iJ)/kAtm2mb,raaPartPress(iJ,iI),
c     $                   raaPartPress(iJ,iI)/(raPActualAvg(iJ)/kAtm2mb)*100.0,raaTemp(iJ,iI),raaAmt(iJ,iI),'||',
c     $                   raaRPartPress(iJ,iI),raaRAmt(iJ,iI),raaRPartPress(iJ,iI)/raaPartPress(iJ,iI)
          END DO
	  write(kStdWarn,*) '++++++++++++++++++++++++++++++++'
        END DO
 234    FORMAT ('index = ',I3,' gas HITRAN ID = ',I3,' molecule = ',A20)
 
c write sum
        write(iIOUN,*) ' iI  iGasID  Name     total molecules/cm2'
	write(iIOUN,*) '--------------------------------------------------'
        DO iI = 1,iNumGases
	  write(iIOUN,7172) iI,iaGases(iI),caGID(iaGases(iI)),raSumTotalGasAmt(iaGases(iI))*kAvog
	END DO
	write(iIOUN,*) '--------------------------------------------------'
 
c then output list of paths to be output
        WRITE(iIOUN,*) 'list of paths to be output ...'
        iP=0
        DO iI=1,iOutTypes
          IF (iaPrinter(iI) .EQ. 1) THEN
            iP=iP+1
            iOutputOptionNum=iI
            IF (iaNp(iI) .GT. iNumGases*kProfLayer) THEN
              write(kStdErr,*)'Cannot have more paths to be output than '
              write(kStdErr,*)'iNumGases*kProfLayer'
              write(kStdErr,*)iI,iaNp(iI)
              CALL DoSTOP 
            END IF
            IF (iaNp(iI) .GT. 0) THEN
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
        IF (iP .GT. 1) THEN
          write(kStdErr,*)'Found more than 1 set of options in *OUTPUT where'
          write(kStdErr,*)'iDat set to 1 (outputting single paths)'
          write(kStdErr,*)'... exiting program!'
          CALL DoSTOP 
        END IF
        IF (iP .EQ. 0) THEN
c inform user no single set of paths to be output
          WRITE(iIOUN,*) iP
        END IF
        IF (iP .EQ. 1) THEN
          iOK=CheckDoubleEntry(iaaOp,iOutputOptionNum,iNumLay)
          IF (iOK .LT. 0) THEN
           write(kStdErr,*)'On checking list of paths to output, have found '
           write(kStdErr,*)'some to be output more than once. messing up'
           write(kStdErr,*)'everything in readatmos.m .. please try again'
           CALL DoSTOP
         END IF
       END IF
          
c then output mixed paths -----------------------------------------------
        WRITE(iIOUN,*) '***********************************************'
        WRITE(iIOUN,*) 'MIXED PATHS'
        WRITE(iIOUN,*) 'Number of mixed paths = ',iNpmix
        WRITE(iIOUN,*) 'Number of uncommented lines in mixtable=',
     $                iMixFileLines
        IF (iNpMix .GT. 0) THEN
          DO iI=1,iMixFileLines
            WRITE(iIOUN,*) caaMixFileLines(iI)
          END DO

c then output list of mixed paths temperatures
          WRITE(iIOUN,*) 'mixed path temperatures ....'
          WRITE(iIOUN,*) (raMixVT(iI),iI=1,iNpmix)
c then output list of mixed paths to be printed
          WRITE(iIOUN,*) 'list of mixed paths to be output ...'
          iP=0
          DO iI=1,iOutTypes
            IF (iaPrinter(iI) .EQ. 2) THEN
              iP=iP+1
              iOutputOptionNum=iI
              IF (iaNp(iI) .GT. iNpmix) THEN
                write(kStdErr,*)'Error! Cannot print out more mixed paths'
                write(kStdErr,*)'than iIpmix'
                CALL DoSTOP
              END IF
              IF (iaNp(iI) .GT. 0) THEN
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
          IF (iP .GT. 1) THEN
            write(kStdErr,*)'Found more than 1 options set in *OUTPUT where'
            write(kStdErr,*)'iDat set to 2 (outputting mixed paths)'
            write(kStdErr,*)'... exiting program!'
            CALL DoSTOP 
          END IF
          IF (iP .EQ. 0) THEN
c inform user no mixed set of paths to be output
            WRITE(iIOUN,*) iP
          END IF
          IF (iP .EQ. 1) THEN
            iOK=CheckDoubleEntry(iaaOp,iOutputOptionNum,iNumLay)
            IF (iOK .LT. 0) THEN
              write(kStdErr,*)'On checking list of MIXED paths to output,found'
              write(kStdErr,*)'some that will be output more than once'
              write(kStdErr,*)'Please check *OUTPUT'
              CALL DoSTOP
            END IF
          END IF

      END IF

c now print the atmosphere information -------------------------------------
        WRITE(iIOUN,*) '***********************************************'
        WRITE(iIOUN,*) 'ATMOSPHERES'
        WRITE(iIOUN,*) 'Numb of atmospheres read from *RADFIL = ',iNatm
        WRITE(iIOUN,*) 'Max numb of emiss. pts = ',kEmsRegions
        IF (iNatm .LT. iNatm2) THEN
          write(kStdErr,*)'RADFIL had information for ',iNatm,' atmospheres'
          write(kStdErr,*)'OUTPUT says there is info for ',iNatm2,' atms!'
          write(kStdErr,*)'please recheck and retry!!!'
          CALL DoSTOP
        END IF

        DO iI=1,iNatm
          IF (iaNumLayers(iI) .GT. iNpmix) THEN
            write(kStdErr,*)'atm#',iI,' has too many layers!! (> ',iNpmix,')'
            CALL DoSTOP
          END IF
          IF (iaNumLayers(iI) .LT. 0) THEN
            write(kStdErr,*)'atm#',iI,' has 0 or less layers'
            CALL DoSTOP
          END IF
          WRITE(iIOUN,*) 'atm#',iI,' has',iaNumLayers(iI),' MP layers :'
          WRITE(iIOUN,*) (iaaRadLayer(iI,iJ),iJ=1,iaNumLayers(iI))
          WRITE(iIOUN,*) 'TSpace,TSurf,SatAngle,SatHeight = ',
     $        raTSpace(iI),raTSurf(iI),raSatAngle(iI),raSatHeight(iI)
          WRITE(iIOUN,*) 'kSolar,kSolarAngle,kSolarRefl = ',
     $     iakSolar(iI),rakSolarAngle(iI),rakSolarRefl(iI)
          WRITE(iIOUN,*) 'kThermal,kThermalAngle,kThermalJacob = ',
     $         iakThermal(iI),rakThermalAngle(iI),iakThermalJacob(iI)
c output the emissivities for this atmosphere, first freq pt then ems
          WRITE(iIOUN,*) 'The surface freqs/emissivities are : '
          WRITE(iIOUN,*) iaSetEms(iI)
          DO iJ=1,iaSetEms(iI)
            WRITE(iIOUN,*)raaaSetEmissivity(iI,iJ,1),
     $                    raaaSetEmissivity(iI,iJ,2)
          END DO
c now output the list of radiances to be printed, for this atmosphere
          iP=0
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
              iP=iP+1
              iOutputOptionNum=iJ
c              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .GT. 1)) THEN
c                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
c                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
c                CALL DoSTOP
c              END IF
c              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .LE. -1)) THEN
c                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
c                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
c                CALL DoSTOP
c              END IF
              IF (iaNp(iJ) .LE. -1) THEN  !!!dumping out rads for ALL layers
c                IF ((iNumNLTEGases .GT. 0) .AND. (iDoUpperAtmNLTE .GT. 0))THEN
c                  iDumpAllUARads = iDumpAllUARads + 1
c                  write(kStdWarn,*) iJ,iaNp(iJ)
c                  write(kStdWarn,*) 'looks like we will dump US rads'
c                END IF                 
                IF ((iDumpAllUARads .GT. 1) .AND.(iDoUpperAtmNLTE .LE. 0))THEN
                  write(kStdErr,*) 'this is too complicated!!!'
                  write(kStdErr,*) 'NLTE code assumes ONE radiating atm'
                  CALL DoStop
                END IF                 
              END IF
              IF (iaNp(iJ) .GT. iaNumLayers(iI)) THEN
                write(kStdErr,*)'Atm #',iI,'too many radiances to be printed!'
                CALL DoSTOP
              END IF
              IF (iaNp(iJ) .EQ. 0) THEN
                write(kStdErr,*) 'Atm #',iI,' has 0 radiances to be printed!'
                CALL DoSTOP
              END IF
              IF (iaNp(iJ) .GT. 0) THEN
                iEnd=iaNp(iJ)
                iaOutNumbers(iJ)=iEnd
                iTotalStuff = iTotalStuff + iEnd
                WRITE(iIOUN,*) '# of radiances to be printed=',iEnd
                WRITE(iIOUN,*) (iaaOp(iJ,iK),iK=1,iEnd)        
                WRITE(iIOUN,*) (raaUserPress(iJ,iK),iK=1,iEnd)        
                iOk=1
                DO iK=1,iEnd
                  IF ((iaaOp(iJ,iK) .GT. iaNumLayers(iI)) .OR.
     $                 (iaaOp(iJ,iK) .LE. 0))  THEN
                    iOk=-1
                  END IF
                END DO
                IF (iOk .LT. 0) THEN
                  write(kStdErr,*)'Atm# ',iI,' has invalid radiance to output!'
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
          IF (iP .EQ. 0) THEN
c inform user no mixed paths to be output for this atmosphere
            WRITE(iIOUN,*) iP
          END IF
          IF (iP .GT. 1) THEN
            write(kStdErr,*)'Have found more than 1 set of options in *OUTPUT '
            write(kStdErr,*)'where iDat set to 3 (outputting radiances)'
            write(kStdErr,*)'for atmosphere #',iI,'... exiting program!'
            CALL DoSTOP 
          END IF
 
        END DO
ccccc        CLOSE(iIOUN)
        WRITE(iIOUN,*) 'END SUMMARY OF INPUT DRIVER FILE ... '
        WRITE(iIOUN,*) '     '
        WRITE(iIOUN,*) '     '

c end if kLongOrShort > 0
      END IF
c--------------- OUTPUT BINARY FILE --------------------------------

      iTag = -1
      DO iI = 1,kW
        IF ((rFrLow .ge. kaMinFr(iI)) .AND. (rFrHigh .le. kaMaxFr(iI))) THEN
          iTag = iI
        END IF
      END DO        
      IF ((iTag .LT. 1) .OR. (iTag .GT. kW)) THEN
        write (kStdErr,*) 'Could not find kaTag for FreqStart,FreqStop!!'
        CALL DoStop
      END IF

      IF ((kRTP .GE. 0) .AND. (kWhichScatterCode .EQ. 5)) THEN
        !! write out cloud info, straight from RTP file and after manipulation
	DO iI = 1,80
	  caFCloudName(iI:iI) = ' '
	END DO
        caFCloudName = caOutName
	iJ = 80
	DO WHILE ((caFCloudName(iJ:iJ) .EQ. ' ') .AND. (iJ. GT. 0))
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
        write(kStdWarn,*) ' '
        write(kStdWarn,*)'  after expand_scatter'
        write(kStdWarn,*)'    Cloud1  (before/after)        Cloud2 (before/after)'
        write(kStdWarn,*)'-------------------------------------------------------'            
        DO iI=1,7
          write(kStdWarn,*) caStrJunk(iI),raaRTPCloudParams0(1,iI),raaRTPCloudParamsF(1,iI),'<>',
     $raaRTPCloudParams0(2,iI),raaRTPCloudParamsF(2,iI)
        END DO

        OPEN(UNIT=iIOUN_Cloud,FILE=caFCloudName,FORM='FORMATTED',STATUS='NEW',
     $       IOSTAT=iFileErr)
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,103) iFileErr, caFCloudName
          write(kStdErr,*)'make sure the cloud dump file does not exist!' 
          CALL DoSTOP
        END IF
	kTempUnitOpen = 1
	write(iIOUN_Cloud,*) '% rows 1-7 are ctype(101/201/301=W/I/A),cprtop(mb),cprbot(mb),cngwat(g/m2),cpsize(um),cfrac and cfrac12'
	write(iIOUN_Cloud,*) '% cols 1-4 are CLOUD 1 (old/new) and CLOUD 2 (old/new)'
        DO iI=1,7
          write(iIOUN_Cloud,*) raaRTPCloudParams0(1,iI),raaRTPCloudParamsF(1,iI),raaRTPCloudParams0(2,iI),raaRTPCloudParamsF(2,iI)
        END DO
        CLOSE(iIOUN_Cloud)
	kTempUnitOpen = -1
      END IF
c      write(kStdWarn,*) 's_writefile debug jpl clouds sarta vs kcarta'
c      call dostop	
      
c next open unformatted file as a fresh file to be written to
      IF (iIOUN1 .NE. 6) THEN
        OPEN(UNIT=iIOUN1,FILE=caOutName,STATUS='NEW',
     $    FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,103) iFileErr, caOutName
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
        END IF
      END IF
 103  FORMAT('ERROR! number ',I5,' opening binary data file : ',/,A80)

      kStdkCartaOpen=1
      write(kStdWarn,*) 'Opened following file for general binary output : '
      write(kStdWarn,*) caOutName
       
      IF (kLongOrShort .EQ. 0) THEN
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

        GOTO 9999
      END IF
        
c    or if kLongOrShort .EQ. -1,+1 then do the standard binary file
c first output general info -----------------------------------------
      WRITE(iIOUN1) caVersion
      WRITE(iIOUN1) kProfLayer
      WRITE(iIOUN1) kMaxUserSet
      WRITE(iIOUN1) (raParams(iI),iI=1,kMaxUserSet)
      WRITE(iIOUN1) caComment
      WRITE(iIOUN1) rFrLow,rFrHigh
      WRITE(iIOUN1) iFileIDLo,iFileIDHi
      WRITE(iIOUN1) kLongOrShort
c v1.04 to v1.08 had the next line; now remove it
c      WRITE(iIOUN1) M1000mb,M100mb,MSubLayer,M50mb,M10mb,MThickLayer
      WRITE(iIOUN1) (raPactualAvg(iI),iI=1,kProfLayer)        

c then output path ID stuff ------------------------------------------
      IF (kLongOrShort .GT. 0) THEN
        WRITE(iIOUN1) iNumGases*kProfLayer
        DO iI=1,iNumGases
          DO iJ=1,kProfLayer
            iP=(iI-1)*kProfLayer+iJ
            WRITE(iIOUN1)iP,iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
          END DO
        END DO
      END IF
c then output list of paths to be output
      iP=0
      DO iI=1,iOutTypes
        IF (iaPrinter(iI) .EQ. 1) THEN
          iP=iP+1
          IF (iaNp(iI) .GT. iNumGases*kProfLayer) THEN
            write(kStdErr,*)'Error! Cannot have more paths to be output than '
            write(kStdErr,*)'iNumGases*kProfLayer'
            CALL DoSTOP 
          END IF
          IF (iaNp(iI) .GT. 0) THEN
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
      IF (iP .GT. 1) THEN
        write(kStdErr,*)'Have found more than 1 set of options in *OUTPUT'
        write(kStdErr,*)'where iDat set to 1 (outputting single paths)'
        write(kStdErr,*)'... exiting program!'
        CLOSE(iIOUN1)
        CALL DoSTOP 
      END IF
      IF (iP .EQ. 0) THEN
c inform user no single set of paths to be output
        WRITE(iIOUN1) iP
      END IF

c then output mixed paths -----------------------------------------------
c note this important CHANGE !!!
      WRITE(iIOUN1) iNpmix
      IF (kLongOrShort .GT. 0) THEN 
        WRITE(iIOUN1) iMixFileLines
        IF (iNpMix .GT. 0) THEN
          DO iI=1,iMixFileLines
            WRITE(iIOUN1) caaMixFileLines(iI)
          END DO
c then output list of mixed paths temperatures
          WRITE(iIOUN1) (raMixVT(iI),iI=1,iNpmix)
        END IF
      END IF

c this is EXTRA if statement, since above one was commented out
      IF (iNpMix .GT. 0) THEN
c then output list of mixed paths to be printed
        iP=0
        DO iI=1,iOutTypes
          IF (iaPrinter(iI) .EQ. 2) THEN
            iP=iP+1
            IF (iaNp(iI) .GT. iNpmix) THEN
              write(kStdErr,*)'Error! Cannot print out more mixed paths'
              write(kStdErr,*)'than iIpmix'
              CALL DoSTOP
            END IF
            IF (iaNp(iI) .GT. 0) THEN
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

        IF (iP .GT. 1) THEN
          write(kStdErr,*)'Have found more than 1 options set in *OUTPUT where'
          write(kStdErr,*)'iDat set to 3 (outputting mixed paths)'
          write(kStdErr,*)'... exiting program!'
          CLOSE(iIOUN1)
          CALL DoSTOP 
        END IF

        IF (iP .EQ. 0) THEN
c inform user no mixed set of paths to be output
          WRITE(iIOUN1) iP
        END IF
      END IF

c now print the atmosphere information -------------------------------------
      WRITE(iIOUN1) iNatm
      WRITE(iIOUN1) kEmsRegions
      DO iI=1,iNatm
        IF (iaNumLayers(iI) .GT. iNpmix) THEN
          write(kStdErr,*)'atm#',iI,' has too many layers!! (> ',iNpmix,')'
          CALL DoSTOP
        END IF
        IF (iaNumLayers(iI) .LT. 0) THEN
          write(kStdErr,*)'atm#',iI,' has 0 layers'
          CALL DoSTOP
        END IF
        WRITE(iIOUN1) iI,iaNumLayers(iI)
        WRITE(iIOUN1) (iaaRadLayer(iI,iJ),iJ=1,iaNumLayers(iI))
        WRITE(iIOUN1) raTSpace(iI),raTSurf(iI),raSatAngle(iI),
     $                                         raSatHeight(iI) 
        WRITE(iIOUN1) iakSolar(iI),rakSolarAngle(iI),rakSolarRefl(iI),
     $     iakThermal(iI),rakThermalAngle(iI),iakThermalJacob(iI)

c output the emissivities for this atmosphere, first freq pt then ems
        WRITE(iIOUN1) iaSetEms(iI)
        DO iJ=1,iaSetEms(iI)
          WRITE(iIOUN1)raaaSetEmissivity(iI,iJ,1),
     $                    raaaSetEmissivity(iI,iJ,2)
        END DO
c now output the list of mixed paths to be printed, for this atmosphere
        iP=0
        DO iJ=1,iOutTypes
          IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
            iP=iP+1
            IF (iaNp(iJ) .GT. iNpmix) THEN
c            IF (iaNp(iJ) .GT. iaNumLayers(iI)) THEN
              write(kStdErr,*)'Atm #',iI,' has too many layers to be printed!'
              CALL DoSTOP
            END IF
            IF (iaNp(iJ) .EQ. 0) THEN
              write(kStdErr,*)'Atm #',iI,' has 0 layers to be printed!'
              CALL DoSTOP
            END IF
            IF (iaNp(iJ) .GT. 0) THEN
              iEnd=iaNp(iJ)
              iaOutNumbers(iJ)=iEnd
              WRITE(iIOUN1) iEnd
              WRITE(iIOUN1) (iaaOp(iJ,iK),iK=1,iEnd)        
              WRITE(iIOUN1) (raaUserPress(iJ,iK),iK=1,iEnd)        
              iOk=1
              DO iK=1,iEnd
                IF ((iaaOp(iJ,iK) .GT. iaNumLayers(iI)) .OR.
     $                 (iaaOp(iJ,iK) .LE. 0))  THEN
                  iOk=-1
                END IF
              END DO
              IF (iOk .LT. 0) THEN
                write(kStdErr,*)'Atm# ',iI,' has invalid layer to be output!'
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

              IF ((iNumNLTEGases .GT. 0) .AND. (iDoUpperAtmNLTE .GT. 0)) THEN
                iDumpAllUARads = iDumpAllUARads + 1
                write(kStdWarn,*) iJ,iaNp(iJ)
                write(kStdWarn,*) 'looks like we will dump UA rads'
              END IF                 
              IF ((iDumpAllUARads .GT. 1) .AND.(iDoUpperAtmNLTE .LE. 0)) THEN
                write(kStdErr,*) 'this is too complicated!!!'
                write(kStdErr,*) 'NLTE code assumes ONE radiating atmosphere'
                CALL DoStop
              END IF                 

            END IF
          END IF
        END DO
        IF (iP .EQ. 0) THEN
c inform user no mixed paths to be output for this atmosphere
          WRITE(iIOUN1) iP
        END IF
        IF (iP .GT. 1) THEN
          write(kStdErr,*)'Have found more than 1 set of options in *OUTPUT '
          write(kStdErr,*)'where iDat set to 3 (outputting radiances)'
          write(kStdErr,*)'for atmosphere #',iI,'... exiting program!'
          CLOSE(iIOUN1)
          CALL DoSTOP 
        END IF
      END DO

c having gotten this far, this tells the reader how many things to expect
c total num kcomp files, total number of output options
c for each output option, how many prints to expect
      write(iIOUN1) iTotal,iOutTypes   
      write(iIOUN1) (iaOutNumbers(iI),iI=1,iOutTypes)

cnew do not open and close!!
c      IF (iIOUN1 .NE. 6) THEN
c        CLOSE(iIOUN1)
c      END IF

 4000 FORMAT(A130)

c--------------- COL JACOBIAN BINARY FILE --------------------------------
      IF (kJacobian .GE. 0 .AND. kActualJacs .GE. 100) THEN     
        write(kStdWarn,*) 'opening column,stemp (small) jacobian output file'
        kStdJacob2 = kStdJacob2KK
        iIOUN_JAC2 = kStdJacob2

        IF (iOutFileName .LT. 0) THEN  !std kcarta jacob2 output dumped to screen
          caJacobFile2 = 'columnjacob.dat'
        ELSE !std kcarta output dumped to file, so find flux file name
          CALL jacob2Name(caJacobFile2,caJacobFile)
        END IF

        OPEN(UNIT=iIOUN_JAC2,FILE=caJacobFile2,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,403) iFileErr, caJacobFile2
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
        END IF

 403  FORMAT('ERROR! number ',I5,' opening COL JACOBIAN binary file : 
     $      ',/,A80)

        kStdJacob2Open=1
        write(kStdWarn,*) 'Opened file for col jacobian binary output : '
        write(kStdWarn,*) caJacobFile2

c now essentially dump out header that resmebles kLongOrShort == 0
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
c--------------- JACOBIAN BINARY FILE --------------------------------
      IF (kJacobian .GE. 0 .AND. kActualJacs .LT. 100) THEN     
        write(kStdWarn,*) 'opening regular (large) jacobian output file (X(z))'
        write(kStdWarn,*) 'kJacobOutput = ',kJacobOutput
        iIOUN2 = kStdJacob

        IF (iIOUN2 .NE. 6) THEN
c open unformatted file as a fresh file to be written to
          OPEN(UNIT=iIOUN2,FILE=caJacobFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
          IF (iFileErr .NE. 0) THEN
            WRITE(kStdErr,203) iFileErr, caJacobFile
            write(kStdErr,*)'make sure the file does not exist!' 
            CALL DoSTOP
          END IF
        END IF

 203  FORMAT('ERROR! number ',I5,' opening JACOBIAN binary file : 
     $      ',/,A80)

        kStdJacobOpen=1
        write(kStdWarn,*) 'Opened following file for jacobian binary output : '
        write(kStdWarn,*) caJacobFile

c write general header information
        WRITE(iIOUN2) caComment
        WRITE(iIOUN2) kProfLayer
        WRITE(iIOUN2) rFrLow,rFrHigh
        WRITE(iIOUN2) iFileIDLo,iFileIDHi

c write some specific info, which will be used for each 25 cm-1 chunk, 
c for each atmosphere
       
c first figure out how many gasid's to output
        iImportant=iJacob
        WRITE(iIOUN2) iImportant
        write(kStdWarn,*) ' iImportant gases = ',iImportant

c then figure out, of the atmospheres that have been read in, which actually 
c have a radiance and hence jacobian calculation associated with them
c assume all error checking done in section above (when blah.dat is created)
        iNatmJac=0
        DO iI=1,iNatm
c now output the list of mixed paths to be printed, for this atmosphere
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
              iNatmJac=iNatmJac+1
              IF ((iaNp(iJ) .GT. 0) .OR. (iaNp(iJ) .EQ. -1)) THEN
c set the number of layers in this atmosphere
                iaLayerJac(iNatmJac)=iaNumLayers(iI)
              END IF
            END IF
          END DO
        END DO

        WRITE(iIOUN2) iNatmJac
        WRITE(iIOUN2) (iaLayerJac(iI),iI=1,iNatmJac)

        write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
        write(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)

        iI=1
 205    CONTINUE
        iP=DoGasJacob(iaGases(iI),iaJacob,iJacob)
        IF ((iI .LE. iNumGases) .AND. (iP .GT. 0)) THEN
          DO iJ=1,kProfLayer
            WRITE(iIOUN2)iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
          END DO
          iI=iI+1
          GOTO 205
        END IF
        IF (iI .LT. iNumGases) THEN
          iI=iI+1
          GOTO 205
        END IF

        !!check to see if we are outputting d/d(DME), d/d(IWP)
        IF ((kJacobian .GT. 0) .AND. (kScatter .GT. 0)) THEN  
          DO iJ=1,kProfLayer
            WRITE(iIOUN2) 201,-300.0,01.0
          END DO
          DO iJ=1,kProfLayer
            WRITE(iIOUN2) 201,-300.0,01.0
          END DO
        END IF

c recall that at the end, we also compute d/dSurface_temp,d/dSurface_Emis
c and d Thermal BackGnd/d(Surface_Emis)
c and d SolarBackGnd/d(Surface_Emis)

c having gotten this far, this tells the reader how many things to expect
c total num kcomp files, total number of output options=iNatm,number of gases
c for each atmosphere, how many layers
        write(iIOUN2) iTotal,iNatmJac,iImportant
        write(iIOUN2) (iaLayerJac(iI),iI=1,iNatmJac)               

      END IF

c--------------- FLUX BINARY FILE -----------------------------------------
c----- the opening header info is almost same as that for Jacobian files ---
      IF (kFlux .GT. 0) THEN     

        iIOUN_Flux = kStdFlux

        IF (iOutFileName .LT. 0) THEN  !std kcarta output dumped to screen
          caFluxFile='flux.dat'
        ELSE !std kcarta output dumped to file, so find flux file name
          CALL FluxName(caFluxFile,caOutName)
        END IF

c open unformatted file as a fresh file to be written to
        OPEN(UNIT=iIOUN_Flux,FILE=caFluxFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,303) iFileErr, caFluxFile
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
        END IF

 303  FORMAT('ERROR! number ',I5,' opening FLUX binary file : 
     $      ',/,A80)

        kStdFluxOpen=1
        write(kStdWarn,*) 'Opened following file for flux ouput : '
        write(kStdWarn,*) caFluxFile

c write general header information
        WRITE(iIOUN_Flux) caComment
c        IF ((kFlux .NE. 6) .OR. (kFlux .NE. 1,2,3)) THEN
c          WRITE(iIOUN_Flux) kProfLayer
c        ELSE
c          WRITE(iIOUN_Flux) kProfLayer+1
c        END IF
        WRITE(iIOUN_Flux) kProfLayer
        WRITE(iIOUN_Flux) rFrLow,rFrHigh
        WRITE(iIOUN_Flux) iFileIDLo,iFileIDHi
ccccccc this is from Jacobians 
ccc first figure out how many gasid's to output
ccc        iImportant=iJacob
ccc        WRITE(iIOUN_Flux) iImportant
ccc        write(kStdWarn,*) ' iImportant gases = ',iImportant
c figure out how many types of fluxes to output (duh!!!!!!!)
        iImportant = 1
        WRITE(iIOUN_Flux) iImportant
        write(kStdWarn,*) 'Up-Down Fluxes = ',iImportant

c then figure out, of the atmospheres that have been read in, which actually 
c have a radiance and hence jacobian calculation associated with them
c assume all error checking done in section above (when blah.dat is created)
        iNatmJac = 0
        DO iI=1,iNatm
c now output the list of mixed paths to be printed, for this atmosphere
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
              iNatmJac = iNatmJac+1
              IF (iaNp(iJ) .NE. 0) THEN   !!ie it is -1 or a positive number
c               IF (iaNp(iJ) .GT. 0) THEN
c set the number of layers in this atmosphere
                iaLayerJac(iNatmJac) = iaNumLayers(iI)
              END IF
            END IF
          END DO
        END DO

        IF (kFLux .LE. 3)  THEN   !!! outputting at all levels
          DO iI = 1,iNatmJac
            iaJunkFlux(iI) = 1 * (iaLayerJac(iI)+1)
          END DO

          WRITE(iIOUN_Flux) iNatmJac
          WRITE(iIOUN_Flux) (iaJunkFlux(iI),iI=1,iNatmJac)

          write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
          write(kStdWarn,*) (iaJunkFlux(iI),iI=1,iNatmJac)

c how many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, dump OLR AND ILR
          write(iIOUN_Flux) iTotal,iNatmJac,iImportant
          write(iIOUN_Flux) (iaJunkFLux(iI),iI=1,iNatmJac)               

        ELSEIF (kFLux .EQ. 4) THEN   !!! outputting only OLR at TOA
          DO iI = 1,iNatmJac
            iaLayerFlux(iI) = 1
          END DO

          WRITE(iIOUN_Flux) iNatmJac
          WRITE(iIOUN_Flux) (iaLayerFlux(iI),iI=1,iNatmJac)

          write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
          write(kStdWarn,*) (iaLayerFlux(iI),iI=1,iNatmJac)

c no need to dump out gas amounts, temperatures; so now tell the reader how 
c many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, how many layers
          write(iIOUN_Flux) iTotal,iNatmJac,iImportant
          write(iIOUN_Flux) (iaLayerFLux(iI),iI=1,iNatmJac)               

        ELSEIF (kFLux .EQ. 5) THEN   !!! outputting only OLR/ILR at TOA/GND/TRPause
          DO iI = 1,iNatmJac
            iaLayerFlux(iI) = 3
          END DO

          WRITE(iIOUN_Flux) iNatmJac
          WRITE(iIOUN_Flux) (iaLayerFlux(iI),iI=1,iNatmJac)

          write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
          write(kStdWarn,*) (iaLayerFlux(iI),iI=1,iNatmJac)

c no need to dump out gas amounts, temperatures; so now tell the reader how 
c many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, dump OLR AND ILR
          write(iIOUN_Flux) iTotal,iNatmJac,iImportant
          write(iIOUN_Flux) (iaLayerFLux(iI),iI=1,iNatmJac)               

        ELSEIF (kFLux .EQ. 6) THEN   !!! outputing OLR/ILR at each layer, plus GND and TOA
          DO iI = 1,iNatmJac
            iaJunkFlux(iI) = 2 * (iaLayerJac(iI)+1)
          END DO

          WRITE(iIOUN_Flux) iNatmJac
          WRITE(iIOUN_Flux) (iaJunkFlux(iI),iI=1,iNatmJac)

          write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
          write(kStdWarn,*) (iaJunkFlux(iI),iI=1,iNatmJac)

c how many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, dump OLR AND ILR
          write(iIOUN_Flux) iTotal,iNatmJac,iImportant
          write(iIOUN_Flux) (iaJunkFLux(iI),iI=1,iNatmJac)               
        END IF

      END IF

c--------------- PLANCK BINARY FILE -----------------------------------------
c----- the opening header info is almost same as that for Jacobian files ---
c----- control is thru kFlux (-1,+1,2,3,4 for off, 0 for on)
      kPlanckOut = kFlux
      IF ((kPlanckOut .EQ. 0) .AND. (iNumNLTEGases .GT.0)) THEN     

        iIOUN_Planck = kStdPlanck

        IF (iOutFileName .LT. 0) THEN  !std kcarta output dumped to screen
          caPlanckFile='planck.dat'
        ELSE !std kcarta output dumped to file, so find planck file name
          CALL PlanckName(caPlanckFile,caOutName)
        END IF

c open unformatted file as a fresh file to be written to
        OPEN(UNIT=iIOUN_Planck,FILE=caPlanckFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,304) iFileErr, caPlanckFile
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
        END IF

 304    FORMAT('ERROR! number ',I5,' opening PLANCK binary file : 
     $      ',/,A80)

        kStdPlanckOpen=1
        write(kStdWarn,*) 'Opened following file for planck output :'
        write(kStdWarn,*) caPlanckFile

c write general header information
        WRITE(iIOUN_Planck) caComment
        WRITE(iIOUN_Planck) kProfLayer
        WRITE(iIOUN_Planck) rFrLow,rFrHigh
        WRITE(iIOUN_Planck) iFileIDLo,iFileIDHi
ccccccc this is from Jacobians 
ccc first figure out how many gasid's to output
ccc        iImportant=iJacob
ccc        WRITE(iIOUN_Planck) iImportant
ccc        write(kStdWarn,*) ' iImportant gases = ',iImportant
c figure out how many types of plancks to output (duh!!!!!!!)
        iImportant=1
        WRITE(iIOUN_Planck) iImportant

c then figure out, of the atmospheres that have been read in, which actually 
c have a radiance calculation associated with them
c assume all error checking done in section above (when blah.dat is created)
c        print *,iNatm
c        print *,(iaNumLayers(iI),iI=1,iNatm)
c        print *,iOutTypes
c        print *,(iaPrinter(iI),iI=1,iOutTypes)
c        DO iJ=1,iOutTypes
c          print *,iJ,iaPrinter(iJ),iaAtmPr(iJ),iaNp(iJ)
c        END DO
c        stop
        
        iNatmJac=0
        DO iI=1,iNatm    
c now output the list of mixed paths to be printed, for this atmosphere
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3) .AND. (iaAtmPr(iJ).EQ.iI)) THEN
              iNatmJac=iNatmJac+1
c              IF (iaNp(iJ) .GT. 0) THEN
              IF (iaNp(iJ) .NE. 0) THEN   !!ie it is -1 or a positive number
c set the number of layers in this atmosphere
                iaLayerJac(iNatmJac)=iaNumLayers(iI)
              END IF
            END IF
          END DO
        END DO

        WRITE(iIOUN_Planck) iNatmJac
        WRITE(iIOUN_Planck) (iaLayerJac(iI),iI=1,iNatmJac)

        write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
        write(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)

c no need to dump out gas amounts, temperatures; so now tell the reader how 
c many things to expect : 
c total num kcomp files, total num of output options=iNatm,number of planckes
c for each atmosphere, how many layers
        write(iIOUN_Planck) iTotal,iNatmJac,iImportant
        write(iIOUN_Planck) (iaLayerJac(iI),iI=1,iNatmJac)               

      END IF

 9999 CONTINUE     !would have jumped here if kLongShort = 0

      IF (iDumpAllUARads .LE. 0) THEN
        iDumpAllUARads = -1
      END IF

c f77 format specifiers : the 1P is to shift things over by 1 decimal point
c http://docs.oracle.com/cd/E19957-01/805-4939/z40007437a2e/index.html
c ../DOC/F77FormatSpecifiers.pdf
 7170 FORMAT(A104)
 7171 FORMAT(I4,' ',I4,' ',3(1P E11.5,' '),0PF11.5,'  ',1P E11.5,A2,2(' ',E11.5),1(' ',E11.3))  
 7172 FORMAT(I3,' ',I3,' ',A20,' ',ES11.5)

      RETURN
      END

c************************************************************************
c this subroutine opens the output UA file

c very similar to OpenBloatFile

c only thing is, it dumps out (iUpper + 1) number of things
c if the user asks for all CO2 opt depths to be dumped out, then also the 
c    UA opt depth spectra are dumped out!!!!!!!!!
c if only specific rads are being dumped out in the regular AIRS atm (0-80 km)
c    then TWO rads are output in this file, one at the regular TOA (0.005 mb)
c    and the other at the new TOA (0.000025 mb)
c    this is if iDumpAllUARads == -1
c if rads at all layers are being dumped out in the regular AIRS atm (0-80 km)
c    then rads at all layers are also output in this file, plus one at the 
c    regular TOA (0.005 mb)
c    this is if iDumpAllUARads == +1

      SUBROUTINE OpenUAFile(iType,iUpperLayers,caUAFile,rFrLow,rFrHigh,
     $                      iDumpAllUASpectra,iDumpAllUARads,
     $                      iFileIDLo,iFileIDHi,iTag,
     $                      iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 caUAFile
      REAL rFrLow,rFrHigh
      INTEGER iFileIDLo,iFileIDHi,iTag,iDumpAllUARads,iDumpAllUASpectra

      INTEGER iUpperLayers !!number of layers in atm #1
      INTEGER iType      !!+1 if regular output; -1 if planck output
      INTEGER iaNLTEChunks(kGasStore)
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

c local vars
      INTEGER iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor,iTotalStuff
      REAL    r1,r2,r3

      iTotalStuff = 0

      IF (iDumpAllUASpectra .GT. 0) THEN
        iTotalStuff = iUpperLayers          !!!dump out all the UA opt depths
      END IF

      IF (iDumpAllUARads .GT. 0) THEN
        !!!rads at 0.005 mb and at all layers upto 0.000025 mb
        iTotalStuff = iTotalStuff + (1 + iUpperLayers)
      ELSE
        !!! default = rads at 0.005 mb and at 0.000025 mb
        iTotalStuff = iTotalStuff + (1 + 1)
      END IF

c find number of chunks to dump
      iNLTEChunks = -1
      i1 = +123456
      i2 = -123456
      DO iI = 1,iNumNLTEGases
        IF (iaNLTEChunks(iI) .GE. iNLTEChunks) THEN
          iNLTEChunks = iaNLTEChunks(iI)
        END IF
        DO iJ = 1,iaNLTEChunks(iI)
          IF (iaaNLTEChunks(iI,iJ) .LE. i1) THEN
            i1 = iaaNLTEChunks(iI,iJ)
          END IF
          IF (iaaNLTEChunks(iI,iJ) .GE. i2) THEN
            i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
          END IF
        END DO
      END DO
      r1 = i1*1.0
      r2 = i2*1.0 
      iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
 
      !!!now check to see we actually start from r1, while running kCARTA!
      IF (r1 .LT. rFrLow) THEN
        r1 = rFrLow
        i1 = iFloor(r1)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
      END IF

      !!!now check to see we actually do go to r2, while running kCARTA!
c      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
c        r2 = rFrHigh-kaBlSize(iTag)
      IF (r2 .GT. rFrHigh) THEN
        r2 = rFrHigh
        i2 = iFloor(r2)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
      END IF

      IF (iType .GT. 0) THEN
        iIoun1 = kNLTEOutUA
      ELSEIF (iType .LT. 0) THEN
        iIoun1 = kStdPlanckUA
      ELSE
        write(kStdErr,*) 'Unknown print option for UA files'
        CALL DoStop
      END IF

c open unformatted file as a fresh file to be written to
      OPEN(UNIT=iIoun1,FILE=caUAFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
      IF (iFileErr .NE. 0) THEN
        WRITE(kStdErr,304) iFileErr,iIOUN1,caUAFile
        write(kStdErr,*)'make sure the file does not exist!' 
        CALL DoSTOP
      END IF

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Printing UA FILE summary stats at open : '
      write(kStdWarn,*) caUAFile
      write(kStdWarn,*) 'Regular freqs etc .... '
      write(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
      write(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
      write(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
      write(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff

      WRITE(kStdWarn,*) i1,i2,iNLTEChunks   !!which chunks we are dealing with
                                            !!and how many to expect
      WRITE(kStdWarn,*) r1,r2               !!start/stop freqs for NLTE chunks 

 304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file : 
     $      ',/,A80)

      IF (iType .GT. 0) THEN
        kNLTEOutUAOpen = 1
        write(kStdWarn,*) 'Opened following file for UA general output :'
        write(kStdWarn,*) caUAFile
      ELSEIF (iType .LT. 0) THEN
        kStdPlanckUAOpen = 1
        write(kStdWarn,*) 'Opened following file for UA planck output :'
        write(kStdWarn,*) caUAFile
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
      END

c************************************************************************
c this subroutine opens the output UA file, for bloated version.
c HAS NOT BEEN TESTED AT ALL

c very similar to OpenBloatFile

c only thing is, it dumps out (iUpper + 1) number of things
c if only specific rads are being dumped out in the regular AIRS atm (0-80 km)
c    then TWO rads are output in this file, one at the regular TOA (0.005 mb)
c    and the other at the new TOA (0.000025 mb)
c    this is if iDumpAllUARads == -1
c if rads at all layers are being dumped out in the regular AIRS atm (0-80 km)
c    then rads at all layers are also output in this file, plus one at the 
c    regular TOA (0.005 mb)
c    this is if iDumpAllUARads == +1

      SUBROUTINE OpenBloatUAFile(
     $                      iType,iUpperLayers,caOutUABloatFile,rFrLow,rFrHigh,
     $                      iDumpAllUARads,
     $                      iFileIDLo,iFileIDHi,iTag,
     $                      iNumNLTEGases,iaNLTEChunks,iaaNLTEChunks)

      IMPLICIT NONE     

      include '../INCLUDE/kcarta.param'

      CHARACTER*80 caOutUABloatFile
      REAL rFrLow,rFrHigh
      INTEGER iFileIDLo,iFileIDHi,iTag,iDumpAllUARads

      INTEGER iUpperLayers !!number of layers in atm #1
      INTEGER iType      !!+1 if regular output; -1 if planck output
      INTEGER iaNLTEChunks(kGasStore)
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)

c local vars
      INTEGER iIOUN1,iFileErr,i1,i2,iI,iJ,iNLTEChunks,iFloor,iTotalStuff
      REAL    r1,r2,r3

      iTotalStuff = 2         !!! default = rads at 0.005 mb and at 0.000025 mb
      IF (iDumpAllUARads .GT. 0) THEN
        iTotalStuff = 1 + iUpperLayers      !!!rads at 0.005 mb and at all
                                            !!!layers upto 0.000025 mb
      END IF

c find number of chunks to dump
      iNLTEChunks = -1
      i1 = +123456
      i2 = -123456
      DO iI = 1,iNumNLTEGases
        IF (iaNLTEChunks(iI) .GE. iNLTEChunks) THEN
          iNLTEChunks = iaNLTEChunks(iI)
        END IF
        DO iJ = 1,iaNLTEChunks(iI)
          IF (iaaNLTEChunks(iI,iJ) .LE. i1) THEN
            i1 = iaaNLTEChunks(iI,iJ)
          END IF
          IF (iaaNLTEChunks(iI,iJ) .GE. i2) THEN
            i2 = iaaNLTEChunks(iI,iJ) +  kaBlSize(iTag)
          END IF
        END DO
      END DO
      r1 = i1*1.0
      r2 = i2*1.0 
      iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
 
      !!!now check to see we actually start from r1, while running kCARTA!
      IF (r1 .LT. rFrLow) THEN
        r1 = rFrLow
        i1 = iFloor(r1)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
      END IF

      !!!now check to see we actually do go to r2, while running kCARTA!
c      IF (r2 .GE. rFrHigh-kaFrStep(iTag)) THEN
c        r2 = rFrHigh-kaBlSize(iTag)
      IF (r2 .GT. rFrHigh) THEN
        r2 = rFrHigh
        i2 = iFloor(r2)
        iNLTEChunks = iFloor((r2-kaFrStep(iTag)-r1)/kaBlSize(iTag))+1
      END IF

      IF (iType .GT. 0) THEN
        iIoun1 = kBloatNLTEOutUA
      ELSEIF (iType .LT. 0) THEN
        iIoun1 = kBloatPlanckUA
      ELSE
        write(kStdErr,*) 'Unknown print option for bloat UA files'
        CALL DoStop
      END IF

c open unformatted file as a fresh file to be written to
      OPEN(UNIT=iIoun1,FILE=caOutUABloatFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
      IF (iFileErr .NE. 0) THEN
        WRITE(kStdErr,304) iFileErr,iIOUN1,caOutUABloatFile
        write(kStdErr,*)'make sure the file does not exist!' 
        CALL DoSTOP
      END IF

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Printing Bloat UA FILE summary stats at open : '
      write(kStdWarn,*) caOutUABloatFile
      write(kStdWarn,*) 'Regular freqs etc .... '
      write(kStdWarn,*) '  rFrLow,rFrHigh = ',rFrLow,rFrHigh
      write(kStdWarn,*) '  iFileIDLo,iFileIDHi = ',iFileIDLo,iFileIDHi
      write(kStdWarn,*) '  df, chunksize = ',kaFrStep(iTag),kaBlSize(iTag)
      write(kStdWarn,*) '  iTag,iTotalStuff = ',iTag,iTotalStuff
      write(kStdWarn,*) 'Bloated  freqs etc .... '
      write(kStdWarn,*) '  iType = = ',iType
      write(kStdWarn,*) '  iNumlayers = ',iUpperLayers
      write(kStdWarn,*) '  kBoxCarUse = ',kBoxCarUse
      write(kStdWarn,*) '  df_fine = ',kaFineFrStep(iTag)
      write(kStdWarn,*) 'start,stop high res integers : ',i1,i2
      write(kStdWarn,*) 'number of highres chunks = ',iNLTEChunks
      write(kStdWarn,*) 'start,stop high res doubles : ',r1,r2
      write(kStdWarn,*) ' '

 304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening BLOATED binary file : 
     $      ',/,A80)

      IF (iType .GT. 0) THEN
        kBloatNLTEOutUAOpen = 1
        write(kStdWarn,*) 'Opened following file for UA bloat general output :'
        write(kStdWarn,*) caOutUABloatFile
      ELSEIF (iType .LT. 0) THEN
        write(kStdErr,*) 'ooer cannot open bloated UA planck files'
        CALL DoStop
        kBloatPlanckUAOpen = 1
        write(kStdWarn,*) 'Opened following file for UA bloat planck output :'
c        write(kStdWarn,*) caOutUABloatPlanckFile
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
      WRITE(iIOUN1) sngl(kaFineFrStep(iTag))!!fine freq step size
      WRITE(iIOUN1) i1,i2,iNLTEChunks       !!which chunks we are dealing with
                                            !!and how many to expect
      WRITE(iIOUN1) r1,r2                   !!start/stop freqs for NLTE chunks 
      WRITE(iIOUN1) kaBlSize(iTag)/kBoxCarUse
                                            !!fine 10000 point freq block size
      WRITE(iIOUN1) iTotalStuff*kBoxCarUse  !!fine number of outputs

      RETURN
      END

c************************************************************************
c this subroutine writes out the main header for each results section
c this is a major change from before
c note the program does not keep on opening and closing the file
      SUBROUTINE wrtout_head_uafile(caOutUAFile,
     $                        rFrLow,rFrHigh,raFreq,iTag,
     $                        iPathORRad,iNumberNLTEOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caOutUAFile  = binary output file name
c rFrLow     = lowest wavenumber to be output   
c rFrHigh    = highest wavenumber to be output  
c rDelta     = point spacing
c iMainType  = 1 for path spectra                   iAtmNumber for Jacobians
c              2 for MP spectra                     iAtmNumber for fluxes
c              3 for radiances
c iSubMainType = kLayer2Sp (-2,-1,1,2) for path     GAS ID (1..63) for GasJac
c              = kLayer2Sp (-2,-1,1,2) for MP       0           for Temp Jac
c              = iAtmNumber for radiances           -10         for wgt fcn
c                                                   -20         for surf Jac
c              = +1 for upward flux, -1 for downward flux
c iNumberOut   = number of the relevant spectra to look for
c iPathOrRad   = +1 for CO2 ua path, +3 for ua rad

      CHARACTER*80 caOutUAFile
      REAL rFrLow,rFrHigh,raFreq(kMaxPts)
      INTEGER iTag,iNumberNLTEOut,iPathORRad

      REAL rDelta
      INTEGER iMainType,iSubMainType,iIOUN

      rDelta = kaFrStep(iTag)
      IF (abs((rFrHigh-rFrLow)-(rDelta*kMaxPts)).GE. 2*rDelta) THEN
         write(kStdErr,*) 'Wow! need (rFrHigh-rFrLow)-(rDelta*kMaxPts))
     $ < rDelta'
         write(kStdErr,*) 'rFrHigh,rFrLow,rDelta,kMaxPts,iTag = '
         write(kStdErr,*) rFrHigh,rFrLow,rDelta,kMaxPts,iTag
         CALL DoSTOP
       END IF
      write(kStdWarn,*) 'dump out : ',iNumberNLTEOut,' for iMain=',iMainType

      iIOUN        = kNLTEOutUA

      IF (iPathORRad .EQ. +1) THEN
        iMainType    = 1          !! right now these are CO2 opt depths
        iSubMainType = 1          !! right now assume only 1 atm
      ELSEIF (iPathORRad .EQ. +3) THEN
        iMainType    = 3          !! right now these are rads
        iSubMainType = 1          !! right now assume only 1 atm
      END IF

      IF (kLongOrShort .NE. 0) THEN
        WRITE(iIOUN) iMainType,iSubMainType,iNumberNLTEOut
        WRITE(iIOUN) kMaxPts,rFrLow,rFrHigh,rDelta
      END IF

      RETURN
      END
c************************************************************************
