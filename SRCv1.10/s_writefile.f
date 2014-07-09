c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c***********************************************************************
c**************** MAIN BINARY OUTPUT SUBROUTINES ARE HERE **************
c***********************************************************************
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
c              = kLayer2Sp (-2,-1,1,2) for MP       0              for Temp Jac
c              = iAtmNumber for radiances           -10            for wgt fcn
c                                                   -20            for surf Jac
c              = +1 for upward flux, -1 for downward flux
c iNumberOut   = number of the relevant spectra to look for

      CHARACTER*80 caOutName
      REAL rFrLow,rFrHigh,rDelta
      INTEGER iMainType,iSubMainType,iNumberOut,iIOUN

      IF (abs((rFrHigh-rFrLow)-(rDelta*kMaxPts)).GE. 2*rDelta) THEN
         write(kStdErr,*) 'Wow! need (rFrHigh-rFrLow)-(rDelta*kMaxPts))
     $ < rDelta'
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
      SUBROUTINE wrtout(iIOUN,caOutName,raWaves,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iIOUN      = unit file number
c raInten    = array containing computed data (radiances,spectra,jacobians ..)
c raWaves    = array containing wavenumbers
c caOutName  = binary output file name

      INTEGER iIOUN
      REAL raWaves(kMaxPts),raInten(kMaxPts)
      CHARACTER*80 caOutName

c local variables
      INTEGER iInt

      WRITE(iIOUN) (raInten(iInt),iInt=1,kMaxPts)

c      WRITE (6,1001)caOutName
c 1001 FORMAT('Successfully saved unformatted results to file ',/,A80)
     
      RETURN
      END
c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the individual GAS PATH spectra, according to what kLayer2Sp is set to
c kLayer2Sp = 2  : gas transmittances Layer to space sum(j=i,N) exp(-k(j))
c kLayer2Sp = 1  : gas optical depth Layer to space  sum(j=i,N) (k(j))
c kLayer2Sp = -1 : gas optical depth                 k(i
c kLayer2Sp = -2 : gas transmittances                exp(-k(j)
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_trans_path(raWaves,rFrLow,rFrHigh,
     $           raaGasAbs,iPrinter,
     $           caOutName,
     $           iFileID,
     $           iaPath,iNp,iaOp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaGasAbs  = single gas abs coeffs
c iPrinter   = 1,2 or 3 ... will be 1 if this routine is called
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
c iaPath     = list of the paths corresponding to the current gas 
c iNp        = total number of paths to be output
c iaOp       = list of the paths to be output
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaGasAbs(kMaxPts,kProfLayer)
      INTEGER iPrinter,iFileID
      INTEGER iNp,iaOp(kPathsOut),iaPath(kProfLayer)
      CHARACTER*80 caOutName

c local variables
      INTEGER iInt,iDiv,iDp,iStart,iPath,iLay,DoOutputLayer,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN=kStdkCarta

      iStart=iDiv(iaPath(1),kProfLayer)

c write spectra to unformatted file
c if iPrinter=1 then have to check for valid paths
      DO iLay=1,kProfLayer
c check to see if this path should be output
        iPath=iStart*kProfLayer + iLay
c        IF (iPath .NE. iaPath(iLay)) THEN
c          write(kStdErr,*) 'iPath .NE. iaPath(iLay)' 
c          Call DoSTOP
c          END IF
        iDp=DoOutputLayer(iPath,iNp,iaOp)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*)'outputting GAS abs coeffs at iPath = ',iPath
          IF (kLayer2Sp .EQ. -1) THEN
            DO iInt=1,kMaxPts
              raL2S(iInt)=raaGasAbs(iInt,iLay)
              END DO
            CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
          ELSE IF (kLayer2Sp .EQ. -2) THEN
            DO iInt=1,kMaxPts
              raL2S(iInt)=exp(-raaGasAbs(iInt,iLay))
              END DO
            CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
          ELSE IF (kLayer2Sp .EQ. 1) THEN
            CALL GasOptLayer2Space(raaGasAbs,raL2S,iLay)
            CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
          ELSE IF (kLayer2Sp .EQ. 2) THEN
            CALL GasTranLayer2Space(raaGasAbs,raL2S,iLay)
            CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
            END IF  
          END IF
        END DO

c      WRITE (6,1001)caOutName
c 1001 FORMAT('Successfully saved unformatted PATH results to ',/,A80)
     
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
      SUBROUTINE out_trans_MP(raWaves,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iIpmix,iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iIpmix     = current mixed path number
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iIpmix,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iInt,iPath,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN=kStdkCarta

c write spectra to unformatted file
c this is for mixed paths
      iPath=iIpmix
      write(kStdWarn,*)'output MIXED path Xsects at iIpmix = ',iIpmix
      IF (kLayer2Sp .EQ. -1) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt)=raaSumAbs(iInt,iIpmix)
          END DO
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
      ELSE IF (kLayer2Sp .EQ. -2) THEN
        DO iInt=1,kMaxPts
          raL2S(iInt)=exp(-raaSumAbs(iInt,iIpmix))
          END DO
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
      ELSE IF (kLayer2Sp .EQ. 1) THEN
        CALL MP2SpaceOptDp(raaSumAbs,raL2S,iIpmix,iNpmix)
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
      ELSE IF (kLayer2Sp .EQ. 2) THEN
        CALL MP2SpaceTrans(raaSumAbs,raL2S,iIpmix,iNpmix)
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
        END IF  

c      WRITE (6,1001)caOutName
c 1001 FORMAT('Successfully saved unformatted MP results to ',/,A80)
     
      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH L2Stransmittance spectra FAST. This assumes that in each atmos
c is in layers of 100, and so can do the L2S quickly. if one of "atmospheres"
c is found to have less than 100 layers (i.e. iNpmix is not a multiple of 100),
c then the "upper layers" are filled with spectr of 0, so their transmittance
c is 1

c this assumes kLayer2Sp = 2, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASTL2Strans_MP(raWaves,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iPath,iIOUN,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax,iI
      REAL raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

      iIOUN=kStdkCarta

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
        DO iLay=kProfLayer-1,1,-1
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
          CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
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
          CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
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
      SUBROUTINE out_FASTL2Soptdp_MP(raWaves,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iPath,iIOUN,iI,iAtmBlocks,iWarn,iA,iFr,iLay,iAmax
      REAL raL2S(kMaxPts),raaTempArray(kMaxPts,kProfLayer)

      iIOUN=kStdkCarta

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
        DO iLay=kProfLayer-1,1,-1
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
          CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
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
          CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH optical depth FAST. 

c this assumes kLayer2Sp = -1, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASToptdp_MP(raWaves,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iIOUN,iI,iFr
      REAL raL2S(kMaxPts)

      iIOUN=kStdkCarta

      DO iI=1,iNpMix
        DO iFr=1,kMaxPts
          raL2S(iFr)=raaSumAbs(iFr,iI)
          END DO
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
        END DO

      RETURN
      END

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now output 
c the MIXED PATH transmittance spectra FAST. 

c this assumes kLayer2Sp = -2, iOp = -1
c check to see if we want the raw spectra, or kLayer2Space
      SUBROUTINE out_FASTtrans_MP(raWaves,rFrLow,rFrHigh,
     $           raaSumAbs,iPrinter,
     $           caOutName,
     $           iNpmix,iFileID)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raWaves    = array containin all the frequencies in the current 25 cm-1 block
c rFrLow,rFrHigh = freq start/stop points for the output (if iChunkSize=1, 
c                  these need not correspond to 1,10000)
c raaSumAbs  = mixed path sum of the abs coeffs
c iPrinter   = 1,2 or 3 ... will be 2 if this routine is called
c iNpmix     = number of mixed paths
c iFileID       = which of the 25 cm-1 k-comp files is being processed
c caOutName  = name of binary file that output goes to
      REAL raWaves(kMaxPts),rFrLow,rFrHigh
      REAL raaSumAbs(kMaxPts,kMixFilRows)
      INTEGER iPrinter,iNpmix,iFileID
      CHARACTER*80 caOutName

c local variables
      INTEGER iI,iFr,iIOUN
      REAL raL2S(kMaxPts)

      iIOUN=kStdkCarta

      DO iI=1,iNpMix
        DO iFr=1,kMaxPts
          raL2S(iFr)=exp(-raaSumAbs(iFr,iI))
          END DO
        CALL wrtout(iIOUN,caOutName,raWaves,raL2S)
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
          iUpper=kProfLayer*iJ
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
          iUpper=kProfLayer*iJ
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
c this subroutine appends _FLUX to end of caOutName to get caFluxName
      SUBROUTINE FluxName(caFluxFile,caOutName)

      IMPLICIT NONE

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
      caFluxFile(iInt+1:iInt+5)='_FLUX'
      
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
      SUBROUTINE PrepareOutput(caFRoot,caOutName,caJacobFile,iOutFileName,
     $      rFrLow,rFrHigh,iFileIDLo,iFileIDHi,caComment,
     $      iNumGases,iaGases,raaAmt,raaTemp,raPressLevels,iProfileLayers,
     $      iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVT,
     $      iNatm,iNatm2,iaNumLayers,iaaRadLayer,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      raaaSetEmissivity,iaSetEms,
     $      iOutTypes,iaPrinter,iaAtmPr,iaNp,iaaOp,raaUserPress,
     $      iJacob,iaJacob,
     $      iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $      rakThermalAngle,iakThermalJacob,iaOutNumbers,iTotal)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels, iProfileLayers = actual number of pressure layers
c iOutFileName = does caOutName exist, or is stuff dumped to screen? for fluxes
c caComment   = comment that the user put into *RADNCE
c iJacob      = number of d/dq we will output (>=1)
c iaJacob     = list of GasID's whose  d/dq we will output
c caFRoot     = root name of input file (not really needed here)
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
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaOutNumbers(kMaxPrint),iOutFileName
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      CHARACTER*80 caFRoot,caOutName,caComment
      CHARACTER*130 caaMixFileLines(kProfLayer)
      INTEGER iaPrinter(kMaxPrint),iaAtmPr(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iOutTypes,iMixFileLines
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iFileIDLo,iFileIDHi
      INTEGER iNumGases,iaGases(kMaxGas),iNpmix,iTotal
      REAL raMixVT(kMixFilRows),raaUserPress(kMaxPrint,kProfLayer)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raPressLevels(kProfLayer+1)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaMix(kMixFilRows,kGasStore),rFrLow,rFrHigh
      INTEGER iNatm,iNatm2,iaNumLayers(kMaxAtm),iJacob,iaJacob(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iProfileLayers

      INTEGER iIOUN,iIOUN1,iIOUN2,iI,iJ,iK,iFileErr,iEnd,iP,iOk
      INTEGER iOutputOptionNum,iNumLay,DoGasJacob
      CHARACTER*80 caJacobFile,caFluxFile

      INTEGER CheckDoubleEntry,iImportant
      INTEGER iNatmJac,iaLayerJac(kMaxAtm),iIOUN_Flux
      REAL raParams(kMaxUserSet),raPActualAvg(kProfLayer),rP

      !this is for kLongOrShort = 0
      INTEGER iTag,iTotalStuff      

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
      raParams(7)  = kTempJac     * 1.0
      raParams(8)  = kSurfTemp    * 1.0
      raParams(9)  = kRTP         * 1.0
      raParams(10) = -98765.0  !!dummy variables
      raParams(11) = -98765.0  !!dummy variables
      raParams(12) = -98765.0  !!dummy variables 

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

c if kLongOrShort = 1, output all info, and summarize in kStdWarn (yawwwwwwn!)
c if kLongOrShort = 0, output ONLY DATA, but summarize in kStdWarn (yawwwwwwn!)
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
        WRITE(iIOUN,*) (raPActualAvg(iI),iI=1,kProfLayer)
        WRITE(iIOUN,*) '***********************************************'

c then output path ID stuff ------------------------------------------
        WRITE(iIOUN,*) 'PATH ID INFO'
        WRITE(iIOUN,*) 'iNumPaths = ',iNumGases*kProfLayer
        WRITE(iIOUN,*) 'Path#      GasID       Temperature   Amount'
        WRITE(iIOUN,*) '-------------------------------------------'
        DO iI=1,iNumGases
          DO iJ=1,kProfLayer
            iP=(iI-1)*kProfLayer+iJ
            WRITE(iIOUN,*)iP,iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
            END DO
          END DO
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
              iEnd=kProfLayer*iNumGases
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
c                END IF
c              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .LE. -1)) THEN
c                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
c                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
c                CALL DoSTOP
c                END IF
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

c next open unformatted file as a fresh file to be written to
      IF (iIOUN1 .NE. 6) THEN
        OPEN(UNIT=iIOUN1,FILE=caOutName,STATUS='NEW',
     $    FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,103) iFileErr, caOutName
 103      FORMAT('ERROR! number ',I5,' opening binary data file : 
     $    ',/,A80)
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
          END IF
        END IF
 
      kStdkCartaOpen=1
       
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
        WRITE(kStdWarn,*) iTotalStuff             !!number of outputs

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
            iEnd=kProfLayer*iNumGases
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
            ELSE 
              iEnd=iaNumLayers(iI)
              WRITE(iIOUN1) iEnd
              iaOutNumbers(iJ)=iEnd
              !note how we are flipping between iI and iJ here
              !and instead of outting iaaOp, we are outputting iaaRadLayer
              WRITE(iIOUN1) (iaaRadLayer(iI,iK),iK=1,iEnd)        
              WRITE(iIOUN1) (raaUserPress(iJ,iK),iK=1,iEnd)        
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
c        END IF

 4000 FORMAT(A130)

c--------------- JACOBIAN BINARY FILE --------------------------------
      IF (kJacobian .GE. 0) THEN     

        iIOUN2 = kStdJacob

        IF (iIOUN2 .NE. 6) THEN
c open unformatted file as a fresh file to be written to
          OPEN(UNIT=iIOUN2,FILE=caJacobFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
          IF (iFileErr .NE. 0) THEN
            WRITE(kStdErr,203) iFileErr, caJacobFile
 203        FORMAT('ERROR! number ',I5,' opening JACOBIAN binary file : 
     $      ',/,A80)
            write(kStdErr,*)'make sure the file does not exist!' 
            CALL DoSTOP
            END IF
          END IF

        kStdJacobOpen=1

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
              IF (iaNp(iJ) .GT. 0) THEN
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
 303      FORMAT('ERROR! number ',I5,' opening FLUX binary file : 
     $      ',/,A80)
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
          END IF

        kStdFluxOpen=1

c write general header information
        WRITE(iIOUN_Flux) caComment
        WRITE(iIOUN_Flux) kProfLayer
        WRITE(iIOUN_Flux) rFrLow,rFrHigh
        WRITE(iIOUN_Flux) iFileIDLo,iFileIDHi
ccccccc this is from Jacobians 
ccc first figure out how many gasid's to output
ccc        iImportant=iJacob
ccc        WRITE(iIOUN_Flux) iImportant
ccc        write(kStdWarn,*) ' iImportant gases = ',iImportant
c figure out how many types of fluxes to output (duh!!!!!!!)
        iImportant=1
        WRITE(iIOUN_Flux) iImportant
        write(kStdWarn,*) 'Up-Down Fluxes = ',iImportant

c then figure out, of the atmospheres that have been read in, which actually 
c have a radiance and hence jacobian calculation associated with them
c assume all error checking done in section above (when blah.dat is created)
        iNatmJac=0
        DO iI=1,iNatm
c now output the list of mixed paths to be printed, for this atmosphere
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
              iNatmJac=iNatmJac+1
              IF (iaNp(iJ) .GT. 0) THEN
c set the number of layers in this atmosphere
                iaLayerJac(iNatmJac)=iaNumLayers(iI)
                END IF
              END IF
            END DO
          END DO

        WRITE(iIOUN_Flux) iNatmJac
        WRITE(iIOUN_Flux) (iaLayerJac(iI),iI=1,iNatmJac)

        write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
        write(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)

c no need to dump out gas amounts, temperatures; so now tell the reader how 
c many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, how many layers
        write(iIOUN_Flux) iTotal,iNatmJac,iImportant
        write(iIOUN_Flux) (iaLayerJac(iI),iI=1,iNatmJac)               

        END IF

 9999 CONTINUE     !would have jumped here if kLongShort = 0

      RETURN
      END

c************************************************************************

c this is same as PrepareOutput() above, and was used prior to Dec 2001.
c however, there was a slight bug in that if kLongOrShort = 0, then 
c iaOutNumbers is not set, and you could get wierd stuff in warning.msg like
c      dump out :   0  for iMain=  3          instead of
c      dump out :   1  for iMain=  3 

c this subroutine writes out the header for the summary text file "header.head"
c this subroutine writes out the header for the binary output file ...
c any other writes to this file will be appended to the end

c as the header is written out, the subroutine checks that print option 1 has
c been specified not more than once, print option 3 specified not more than 
c once and print option 7 not more than once for the ith atmosphere

c kLongOrShort enables shorted output files to be written out (-1,+1), 
c or just the basic results file (0)
      SUBROUTINE PrepareOutputWORKS(caFRoot,caOutName,caJacobFile,iOutFileName,
     $      rFrLow,rFrHigh,iFileIDLo,iFileIDHi,caComment,
     $      iNumGases,iaGases,raaAmt,raaTemp,raPressLevels,iProfileLayers,
     $      iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVT,
     $      iNatm,iNatm2,iaNumLayers,iaaRadLayer,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      raaaSetEmissivity,iaSetEms,
     $      iOutTypes,iaPrinter,iaAtmPr,iaNp,iaaOp,raaUserPress,
     $      iJacob,iaJacob,
     $      iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $      rakThermalAngle,iakThermalJacob,iaOutNumbers,iTotal)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels, iProfileLayers = actual number of pressure layers
c iOutFileName = does caOutName exist, or is stuff dumped to screen? for fluxes
c caComment   = comment that the user put into *RADNCE
c iJacob      = number of d/dq we will output (>=1)
c iaJacob     = list of GasID's whose  d/dq we will output
c caFRoot     = root name of input file (not really needed here)
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
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaOutNumbers(kMaxPrint),iOutFileName
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      CHARACTER*80 caFRoot,caOutName,caComment
      CHARACTER*130 caaMixFileLines(kProfLayer)
      INTEGER iaPrinter(kMaxPrint),iaAtmPr(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iOutTypes,iMixFileLines
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iFileIDLo,iFileIDHi
      INTEGER iNumGases,iaGases(kMaxGas),iNpmix,iTotal
      REAL raMixVT(kMixFilRows),raaUserPress(kMaxPrint,kProfLayer)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raPressLevels(kProfLayer+1)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaMix(kMixFilRows,kGasStore),rFrLow,rFrHigh
      INTEGER iNatm,iNatm2,iaNumLayers(kMaxAtm),iJacob,iaJacob(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iProfileLayers

      INTEGER iIOUN,iIOUN1,iIOUN2,iI,iJ,iK,iFileErr,iEnd,iP,iOk
      INTEGER iOutputOptionNum,iNumLay,DoGasJacob
      CHARACTER*80 caJacobFile,caFluxFile

      INTEGER CheckDoubleEntry,iImportant
      INTEGER iNatmJac,iaLayerJac(kMaxAtm),iIOUN_Flux
      REAL raParams(kMaxUserSet),raPActualAvg(kProfLayer),rP

      !this is for kLongOrShort = 0
      INTEGER iTag,iTotalStuff      

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
      raParams(7)  = kTempJac     * 1.0
      raParams(8)  = kSurfTemp    * 1.0
      raParams(9)  = kRTP         * 1.0
      raParams(10) = -98765.0  !!dummy variables
      raParams(11) = -98765.0  !!dummy variables
      raParams(12) = -98765.0  !!dummy variables 

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

c if kLongOrShort = 1, output all info, and summarize in kStdWarn (yawwwwwwn!)
c if kLongOrShort = 0, output ONLY DATA, but summarize in kStdWarn (yawwwwwwn!)
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
        WRITE(iIOUN,*) (raPActualAvg(iI),iI=1,kProfLayer)
        WRITE(iIOUN,*) '***********************************************'

c then output path ID stuff ------------------------------------------
        WRITE(iIOUN,*) 'PATH ID INFO'
        WRITE(iIOUN,*) 'iNumPaths = ',iNumGases*kProfLayer
        WRITE(iIOUN,*) 'Path#      GasID       Temperature   Amount'
        WRITE(iIOUN,*) '-------------------------------------------'
        DO iI=1,iNumGases
          DO iJ=1,kProfLayer
            iP=(iI-1)*kProfLayer+iJ
            WRITE(iIOUN,*)iP,iaGases(iI),raaTemp(iJ,iI),raaAmt(iJ,iI)
            END DO
          END DO
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
              iEnd=kProfLayer*iNumGases
              iTotalStuff = iTotalStuff + iEnd
              WRITE(iIOUN,*) 'Number of Paths to be output=',iEnd
              WRITE(iIOUN,*) (iaaOp(iI,iJ),iJ=1,iEnd)        
              END IF
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
c                END IF
c              IF ((kScatter .GT. 0) .AND. (iaNp(iJ) .LE. -1)) THEN
c                write(kStdErr,*)'Atm #',iI,' too many radiances to print!'
c                write(kStdErr,*)'Scattering included ==> only do TOA radiance'
c                CALL DoSTOP
c                END IF
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

c next open unformatted file as a fresh file to be written to
      IF (iIOUN1 .NE. 6) THEN
        OPEN(UNIT=iIOUN1,FILE=caOutName,STATUS='NEW',
     $    FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
        IF (iFileErr .NE. 0) THEN
          WRITE(kStdErr,103) iFileErr, caOutName
 103      FORMAT('ERROR! number ',I5,' opening binary data file : 
     $    ',/,A80)
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
          END IF
        END IF
 
      kStdkCartaOpen=1
       
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
        WRITE(kStdWarn,*) iTotalStuff             !!number of outputs

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
            iEnd=kProfLayer*iNumGases
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
            ELSE 
              iEnd=iaNumLayers(iI)
              WRITE(iIOUN1) iEnd
              iaOutNumbers(iJ)=iEnd
              !note how we are flipping between iI and iJ here
              !and instead of outting iaaOp, we are outputting iaaRadLayer
              WRITE(iIOUN1) (iaaRadLayer(iI,iK),iK=1,iEnd)        
              WRITE(iIOUN1) (raaUserPress(iJ,iK),iK=1,iEnd)        
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
c        END IF

 4000 FORMAT(A130)

c--------------- JACOBIAN BINARY FILE --------------------------------
      IF (kJacobian .GE. 0) THEN     

        iIOUN2 = kStdJacob

        IF (iIOUN2 .NE. 6) THEN
c open unformatted file as a fresh file to be written to
          OPEN(UNIT=iIOUN2,FILE=caJacobFile,STATUS='NEW',
     $      FORM='UNFORMATTED',IOSTAT=iFileErr)
c if file error, inform user and stop program
          IF (iFileErr .NE. 0) THEN
            WRITE(kStdErr,203) iFileErr, caJacobFile
 203        FORMAT('ERROR! number ',I5,' opening JACOBIAN binary file : 
     $      ',/,A80)
            write(kStdErr,*)'make sure the file does not exist!' 
            CALL DoSTOP
            END IF
          END IF

        kStdJacobOpen=1

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
              IF (iaNp(iJ) .GT. 0) THEN
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
 303      FORMAT('ERROR! number ',I5,' opening FLUX binary file : 
     $      ',/,A80)
          write(kStdErr,*)'make sure the file does not exist!' 
          CALL DoSTOP
          END IF

        kStdFluxOpen=1

c write general header information
        WRITE(iIOUN_Flux) caComment
        WRITE(iIOUN_Flux) kProfLayer
        WRITE(iIOUN_Flux) rFrLow,rFrHigh
        WRITE(iIOUN_Flux) iFileIDLo,iFileIDHi
ccccccc this is from Jacobians 
ccc first figure out how many gasid's to output
ccc        iImportant=iJacob
ccc        WRITE(iIOUN_Flux) iImportant
ccc        write(kStdWarn,*) ' iImportant gases = ',iImportant
c figure out how many types of fluxes to output (duh!!!!!!!)
        iImportant=1
        WRITE(iIOUN_Flux) iImportant
        write(kStdWarn,*) 'Up-Down Fluxes = ',iImportant

c then figure out, of the atmospheres that have been read in, which actually 
c have a radiance and hence jacobian calculation associated with them
c assume all error checking done in section above (when blah.dat is created)
        iNatmJac=0
        DO iI=1,iNatm
c now output the list of mixed paths to be printed, for this atmosphere
          DO iJ=1,iOutTypes
            IF ((iaPrinter(iJ) .EQ. 3).AND.(iaAtmPr(iJ).EQ.iI)) THEN
              iNatmJac=iNatmJac+1
              IF (iaNp(iJ) .GT. 0) THEN
c set the number of layers in this atmosphere
                iaLayerJac(iNatmJac)=iaNumLayers(iI)
                END IF
              END IF
            END DO
          END DO

        WRITE(iIOUN_Flux) iNatmJac
        WRITE(iIOUN_Flux) (iaLayerJac(iI),iI=1,iNatmJac)

        write(kStdWarn,*)'had',iNatmJac,' out of',iNatm,' atm to output'
        write(kStdWarn,*) (iaLayerJac(iI),iI=1,iNatmJac)

c no need to dump out gas amounts, temperatures; so now tell the reader how 
c many things to expect : 
c total num kcomp files, total number of output options=iNatm,number of fluxes
c for each atmosphere, how many layers
        write(iIOUN_Flux) iTotal,iNatmJac,iImportant
        write(iIOUN_Flux) (iaLayerJac(iI),iI=1,iNatmJac)               

        END IF

 9999 CONTINUE     !would have jumped here if kLongShort = 0

      RETURN
      END

c************************************************************************
