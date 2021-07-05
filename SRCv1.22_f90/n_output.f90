! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_output

USE basic_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile

IMPLICIT NONE

CONTAINS

! this is the OUTPUT section
!************************************************************************
! note if SCATTER is on (use RTSPEC or DISORT), then the user can only output
! radiances at instrument height (ie -1 option, or more than one pressure level
! to output radiances at, is not allowed)
! hopefully this will be changed in the not so far future
!************************************************************************
! are ALWAYS defined wrt to a FULL layer (the pressure level boundaries of
! which are defined in the KLAYERS files)

! the following two paragraphs are for upgoing radiation (down look instr)

! eg     if top AIRS layer is from 10 mb to 0.005 mb
!    and if pressure where aircraft flies is  5.0 mb
! then only fraction rF=(10-5)/(10-0.005) ~ 1/2 = rFracTop of the top
! layer is used ....
!    --------------

!    //////////////          use this LOWER portion
!    --------------c
! so if the Print Option for this atm is
! 3
! 1 1
! 0
! implying that one radiance, at aircraft height, is to be output, then
! raaOp(iAtm,1) = 0.5 (=rFracTop)
! if the user had said he wanted a radiance output at pressure 7.5 mb :
! 3
! 1 1
! 7.5
! then raaOp(iAtm,1) = (10-7)/(10-0.005) ~ 0.25

! note that we have to be very careful about the bottommost layer
! eg     if bottom AIRS layer is from 1000 mb to 900 mb
!    and if surface pressure  is 950 mb
! then only fraction rF=(950-900)/(1000-900) ~ 1/2 = rFracBot of the bottom
! layer is used
!    --------------
!    //////////////          use this UPPER portion
!    //////////////

!    --------------
! so if the Print Option for this atm is
! 3
! 1 1
! 925.0
! implying that one radiance, at 925 mb, is to be output, then
! raaOp(iAtm,1) = (925-900)/(1000-9000) = 0.25

! also note that if atmosphere is defined between eg 1000 to 0.005 mb, and the
! user says to output radiance at 1000 mb, then this means that the radiation
! at the surface is output :
! raRad(vu) = raEms(vu) * ttorad(Tsurf,vu) + thermal(vu) + solar(vu)

! also note if user specifies radiances to be output at EACH layer
! 3
! 1 -1
! then for each layer, the radiation at uppermost part of the layer is output
! ie bottom layer : surface to top of layer           --> output rad
!      next layer : bottom of layer to top of layer   --> output rad
!        ....          .......                      ...
!       top layer : bottom of layer to aircraft posn  --> output rad

! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


! SAME CONSIDERATIONS APPLY TO UPWARD LOOKING INSTRUMENT .. have to be very
! careful about topmost layer

! also note that if atmosphere is defined between eg 0.005 to 10000 mb, and
! user says to output radiance at 0.005 mb, then this means that the radiation
! at the very TopOfAtmosphere is output :
! raRad(vu) = ttorad(2.73,vu) + solar(vu)

! also note if user specifies radiances to be output at EACH layer
! 3
! 1 -1
! then for each layer, the radiation at bottommost part of the layer is output
! ie    top layer : TOA  to bottom of layer           --> output rad
!      next layer : top of layer to bottom of layer   --> output rad
!        ....          .......                      ...
!    bottom layer : top of layer to instr posn        --> output rad

! this is for the FUTURE         Sergio Dec 8, 2000
! Note : if kScatterCode = 0     (clear sky)  then the above is always true
!      : if kScatterCode = -1,+1 (cloudy sky) then the above can be false
!         - if the "-1" option is used for pressure levels, then everything
!           is ok ie radiances at ALL pressure levels used
!         - if the "-1" option is not used for pressure levels, then pressures
!           are rounded so that they coincide as closely as possible to the
!           pressure boundaries, so that
!           DISORT and RTSPEC can easily output radiances at layer boundaries

!************************************************************************

! this subroutine deals with the 'OUTPUT' keyword
! relevant output is the printing switch, iPrinter
! if iDat = 3 then read       iAtm==atmospher to be output,
!                             iNp==no of mixed paths to be output
! ELSE if iDat =1,2 then read iWhich==which gas or miexd set to be output
!                             iNp==no of mixed pathes to be output

! option 1 : gas paths
! if user wants to output spectra for ALL gases, ALL 100 layers, file says
!      -1 -1
! else if user wants to output spectra for GasID iG, all 100 layers, file says
!      iG  -1
! else if user wants to output spectra for GasID iG, iN layers, file says
!      iG  iN
!      followed by list of iN layers
! else if user  wants to output spctra for ALL gases, iN layers, file says
!      -1  iN
!      followed by list of iN layers

! option 2 : mixed paths
! if user wants to output spectra for ALL MP sets, ALL 100 layers, file says
!      -1 -1
! else if user wants to output spectra for set iS, all 100 layers, file says
!      iS  -1
! else if user wants to output spectra for set iS, iN layers, file says
!      iS  iN
!      followed by list of iN layers
! else if user  wants to output spctra for ALL sets, iN layers, file says
!      -1  iN
!      followed by list of iN layers

! option 3 : radiances
! if user wants to output radiances for ALL atm, ALL layers, input file says
!       -1 -1
! else if wants to output radiances for Nth atmosphere, all layers, file says
!        N  -1
! else if wants to output radiances for ALL atmosphere, at L presss, file says
!       -1  L
!       followed by list of L pressures
! hence always read in 2 integers

    SUBROUTINE output4(iaPrinter,iaGPMPAtm,iaNp, &
    iaaOp,raaOp,raaUserPress, &
    iaNumLayer,iaaRadLayer,iNatm,iNatm2, &
    iOutTypes,iNumGases,iNpmix, &
    raFracTop,raFracBot,raaPrBdry,raPressLevels, &
    iaGases,caComment)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! caComment  = 120 character user given comment
! iaGases    = array containing info about GAS ID's fouund in MOL/XSCGAS
! iaPrinter  = array containing the print options (1,2 or 3)
! iaGPMPAtm    = array containing which atmospheres to be printed (if print
!              option ii=1 or 2, then iaGPMPAtm(ii)=0)
! iaNp       = array with  the # of paths to be printed for each print option
! iaaOp      = matrix with list of paths/MP/radiances output for each option
! raaOp      = matrix with fractional list of where radiances should be output
! raaUserPress= matrix with user specified pressures
! iaNumLayer = number of layers in each atmosphere, from *RADFIL
! iaaRadLayer= list of layers in each atmosphere, from *RADFIL
! iOutTypes  = total number of outputs per each 25cm-1 run
! iNumGases  = total number of gases read in from GASFIL + XSCFIL
! iNpmix     = total number of ixed paths read in from MIXFIL
! iErr       = flags an error (mainly file I/O)
! iNatm      = number of atmospheres that RADFIL had in it
! iNatm2     = number of atmospheres that OUTPUT claims there should be
! raFracTop/Bot = how the top/bottom layers of mixing table have been affected
!                 by the mixing table
! raaPrBdry     = pressure boundaries of the atmospheres
    INTEGER :: iNatm,iaGases(kMaxGas)
    INTEGER :: iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iaNp(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iNumGases,iOutTypes,iNatm2
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNpmix
    REAL :: raaOp(kMaxPrint,kPathsOut),raaUserPress(kMaxPrint,kProfLayer)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
    REAL :: raPressLevels(kProfLayer+1)    !!!!actual pressures of levels
    CHARACTER(7) :: caWord
    CHARACTER(160) :: caComment

    INTEGER :: iPrinter,iAtm,iNp,iaOp(kPathsOut)
    INTEGER :: iNumLinesRead,iI,iJ,iUpDown,iList,iListType
    REAL :: raTemp(kProfLayer)

! these are temporarily used
    INTEGER :: iaGasID(kMaxGas),iGS,iMPSets,iErr,iOutTypesMax
    INTEGER :: iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
    INTEGER :: iaaOpT(kMaxPrT,kPathsOut),iOut2
    REAL :: raaOpT(kMaxPrT,kPathsOut),raaUserPressT(kMaxPrT,kProfLayer)

    INTEGER :: iaPrinter1(kMaxPrint),iaGPMPAtm1(kMaxPrint)
    INTEGER :: iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint),iMaxMinLBLRTMOp
    REAL :: raaOp1(kMaxPrint,kPathsOut),rMaxMinLBLRTMOp

    IF ((kRTP == -10) .OR. (kRTP == -5) .OR. (kRTP == -6)) THEN
      write (kStdWarn,*) 'Need to reset some output params as they came from text LVLS/LBLRTM code'
      write(kStdWarn,*) 'raaOp(1,1) = ',raaOp(1,1),' --> ',raRTP_TxtInput(6),' mb'
      rMaxMinLBLRTMOp = +1.0e10
      iMaxMinLBLRTMOp = 1
      DO iI=1,kMaxPrint
        DO iJ=1,kProfLayer
          IF (raaOp(iI,iJ) < rMaxMinLBLRTMOp) THEN
            iMaxMinLBLRTMOp = iJ
            rMaxMinLBLRTMOp = raaOp(iI,iJ)
          END IF
        END DO
      raaOp(iI,iMaxMinLBLRTMOp) = raRTP_TxtInput(6)
      END DO
    END IF

    iaPrinter1 = iaPrinter
    iaGPMPAtm1 = iaGPMPAtm
    iaNP1      = iaNP

    iaaOp1 = iaaOp
    raaOp1 = raaOp

    caWord = '*OUTPUT'
    iErr = -1

    iOutTypes = 1
    iNatm2    = -1

    iaGasID = -100

    iJ=0
    DO iI=1,kMaxGas
      IF (iaGases(iI) > 0) THEN
        iJ = iJ+1
        iaGasID(iJ) = iI
      END IF
    END DO
    IF (iNumGases /= iJ) THEN
      write(kStdErr,*) 'discrepancy in number of gases stored'
      write(kStdErr,*) 'iNumGases, iGasStored = ',iNumGases,iJ
      CALL DoSTOP
    END IF

    iNumLinesRead=0
 13 IF (iNumLinesRead > 0) THEN
      iErr=1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
    5010 FORMAT('Error reading section ',A7,' of main file')

    iNumLinesRead=1

! figure out as a first guess how many printing options there have to be
    iOutTypesMax=0
    DO iErr=1,kMaxPrint
      IF (iaPrinter(iErr) > 0) THEN
        iOutTypesMax=iOutTypesMax+1
      END IF
    END DO
          
    iOutTypes=1
    DO iOut2=1,iOutTypesMax
      iPrinter=iaPrinter1(iOut2)
      iaPrinterT(iOutTypes)=iPrinter
      IF ((iPrinter /= 1) .AND. (iPrinter /= 2) .AND. (iPrinter /= 3))THEN
        write(kStdErr,*) 'iPrinterOption = ',iPrinter
        write(kStdErr,*)'in *OUTPUT .. Need iPrinterOption= 1,2 or 3'
        CALL DoSTOP
      END IF

      ! now read which atmosphere, gaspath or MP set to print and do check ...
      IF (iPrinter == 3) THEN
        iGs = -1000
        iAtm = iaGPMPAtm1(iOut2)
        iNp = iaNp1(iOut2)
        CALL Read3(iAtm,iNp,iNpmix,iNatm,iNatm2,iListType, &
            iaGPMPAtmT,iaNpT,iOutTypes,iaPrinterT,iPrinter)
      ELSE IF (iPrinter == 2) THEN
        iAtm = 0                 !make sure this is set to 0
        iGs = iaGPMPAtm1(iOut2)
        iNp = iaNp1(iOut2)
        !***warning*** iNp can be reset from -1 to another number by Read2
        ! see subroutine ReCheck_iNp
        ! if initially -1, it is left as -1 if ALL mixed paths to be output
        CALL Read2(iGS,iNp,iNpmix,iMPSets,iAtm, &
            iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,iNumGases)
      ELSE IF (iPrinter == 1) THEN
        iAtm = 0                 !make sure this is set to 0
        iGs = iaGPMPAtm1(iOut2)
        iNp = iaNp1(iOut2)
        !***warning*** iNp can be reset from -1 to another number by Read1
        ! see subroutine ReCheck_iNp
        ! if initially -1, it is left as -1 if ALL paths to be output
        CALL Read1(iGS,iNp,iNpmix,iMPSets,iAtm, &
            iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter, &
            iNumGases,iaGasID)
      END IF

      ! if iNp=0 then nothing to be printed --- abort program
      IF (iNp == 0) THEN
        write(kStdErr,*) 'in *OUTPUT!! You have chosen a valid printing'
        write(kStdErr,*) 'option (1,2 or 3) but do not specify any '
        write(kStdErr,*) 'paths/mixed paths/radiances to be printed !!!'
        CALL DoSTOP
      END IF

      ! now consider the case where iNp < 0 but we only set ONE single option
      ! iAtm = 0, iNp < 0 ==> we want to print out a set of mixed paths or
      !                                          a set of gas optical depths
      ! iAtm > 0, iNp < 0 ==> want to print out radiances at all levels of an
      !                       atmosphere
      ! &&&&&&&&&& BEGIN CASE 1
      IF ((iAtm >= 0)  .AND. (iNp < 0)) THEN
        IF (iPrinter == 1) THEN
          ! have to output ALL the single paths, corresponding to all KProfLayers for
          ! each of the iNumGases gases present
          DO iI = 1,kProfLayer*iNumGases
            iaOp(iI) = iI
            iaaOpT(iOutTypes,iI) = iI
          END DO
        ELSE IF (iPrinter == 2) THEN
          ! have to output ALL the mixed paths, corresponding to all 1..iNpMix
          DO iI = 1,iNpmix
            iaOp(iI)=iI
            iaaOpT(iOutTypes,iI)=iI
          END DO
          ! if iNp < 0 and iPrinter = 3, then set iaOp to "ALL" values, as follows
        ELSE IF (iPrinter == 3) THEN !set pressure fraction=1.0
          ! this subroutine sets iaOp and iaaOpT
          CALL AllLayersOutputPress(iaaRadLayer,iAtm,iaNumLayer, &
                iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT, &
                raaPrBdry,raFracTop,raFracBot,raPressLevels)
        END IF
      END IF       !IF ((iAtm >= 0)  .AND. (iNp < 0)) THEN
      ! &&&&&&&&&& END CASE 1

      ! &&&&&&&&&& BEGIN CASE 2
      ! now consider the case where iNp > 0
      ! if iNp > 0 then look for list of layers/mixed paths/radiances
      ! &&&&&&&&&& BEGIN CASE 2A
      IF ((iAtm >= 0)  .AND. (iNp > 0)) THEN
        IF (iPrinter == 3) THEN
          ! can easily set things for the ONE atmosphere
          raTemp(1:iNp)=raaOp1(iOut2,1:iNp)
          IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
            iUpDown = 1            !upward travelling radiation
          ELSE
            iUpDown = -1           !downward travelling radiation
          END IF
          CALL DoSortReal(raTemp,iNp,-iUpDown)
          CALL PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iAtm, &
                iOutTypes,raTemp,iNp,raaPrBdry,raaUserPressT,raPressLevels)
          iaaOpT(iOutTypes,1:iNp)=iaOp(1:iNp)
          ! &&&&&&&&&& END CASE 2A

          ! for iPrinter = 1,2 (path or MP)
          ! iList     = set how many items found in list
          ! iListType = 1,2,3,4    1 if -1  -1
          !                        2 if -1  iL
          !                        3 if iGS -1
          !                        4 if iGS iL

        ELSE IF (iPrinter /= 3) THEN
          ! &&&&&&&&&& BEGIN CASE 2B
          IF  ((iListType == 2) .OR. (iListType == 4)) THEN
            IF  (iListType == 2) THEN
              iaOp(1:iList)=iaaOp1(iOut2,1:iList)
              CALL DoSort(iaOp,iList)
              !having read in path/mixed path integers for ALL gas/MP sets
              !now for all gas/MP sets, set the relevant paths/mixed paths
              CALL DoSetAll(iaOp,iList,iNp,iPrinter,iNumGases,iMPSets)
            ELSE IF  (iListType == 4) THEN
              iaOp(1:iNp)=iaaOp1(iOut2,1:iNp)
              CALL DoSort(iaOp,iNp)
              !having read in path/mixed path integers for the gas/MP set
              !now for all gas/MP sets, set the relevant paths/mixed paths
              CALL DoSetSome(iaOp,iList,iNp,iPrinter,iNumGases,iaGasID,iMPSets,iGS)
            END IF
         ! &&&&&&&&&& END CASE 2B

       ELSE IF ((iPrinter == 2) .AND. (iListType == 3)) THEN
         ! &&&&&&&&&& BEGIN CASE 2C
         ! have to set all MPs for a particular MP set
         DO iI = 1,kProfLayer
           iaOp(iI) = iI+(iGS-1)*kProfLayer
         END DO
         ! &&&&&&&&&& END CASE 2C

       ELSE IF ((iPrinter == 1) .AND. (iListType == 3)) THEN  !no list
         ! &&&&&&&&&& BEGIN CASE 2D
         ! have to set all paths for a particular gas
         ! see where in the list gas iGS lies
         iI = 1
 111     CONTINUE
         IF (iaGasID(iI) /= iGS) THEN
           iI = iI+1
           GO TO 111
         END IF
         iGS = iI
         DO iI=1,kProfLayer
           iaOp(iI)=iI+(iGS-1)*kProfLayer
         END DO
       END IF       !if iListype = 2,4 or 3
       ! &&&&&&&&&& END CASE 2D

       END IF         !ELSE IF (iPrinter /= 3) THEN
      ! &&&&&&&&&& END CASE 2
    END IF           !      IF ((iNp > 0) .AND. (iAtm >= 0)) THEN

    ! whew!!! we have now parsed in Print Option (1,2 or 3), the line following
    ! it (iAtm/iGS  iNp) and if necessary, the next line that has the list of iNp
    ! paths/mixed paths/pressures ... the list is stored in iaOpT
    ! if iAtm > 0, everything is almost set .. all we have to do is make sure this
    ! option can be included in iaPrinter(kMaxPrint) etc
    IF (iAtm >= 0) THEN
      CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iAtm,iNpmix, &
            iNumGases,iNp,iaNumLayer)
      ! if there is a duplicate, do not update iOutTypes, else increment
      ! OutTypes by 1
      ! if we find there are more than kMaxPrint options being stored,
      ! then HALT!!!
      IF (iPrinter /= 3) THEN  !thus iAtm was set to 0 (no radiance)
        CALL CheckIaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases, &
                iNp,iaNpT,iaPrinterT)
      END IF
      IF ((iPrinter == 3) .AND. (iAtm > 0)) THEN
        CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases, &
                iNp,iaNpT,iaPrinterT,iNatm,iAtm,iaGPMPAtmT, &
                raaOpT,raaUserPressT,iaNumLayer,iUpDown)
      END IF
    END IF           !if iAtm >= 0)

    ! &&&&&&&&&& BEGIN CASE 3
    ! if iAtm < 0 and iNp > 0, we now have to set the same printing output paths
    ! for all the atmospheres
    IF (iAtm < 0) THEN
      IF (iPrinter /= 3) THEN
        write(kStdErr,*)'Something gone wrong in *OUTPUT'
        CALL DoSTOP
      END IF
      ! &&&&&&&&&& BEGIN CASE 3A
      IF (iNp > 0) THEN
        raTemp(1:iNp) = raaOp1(iOut2,1:iNp)
        DO iJ = 1,iNatm
          IF (iaaRadLayer(iJ,1) < iaaRadLayer(iJ,2)) THEN
            iUpDown = 1            !upward travelling radiation
          ELSE
            iUpDown = -1           !downward travelling radiation
          END IF

          iaGPMPAtmT(iOutTypes) = iJ
          iaNpT(iOutTypes)      = iNp
          iaPrinterT(iOutTypes) = iPrinter
                                
          CALL DoSortReal(raTemp,iNp,-iUpDown)
          CALL PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iJ, &
            iOutTypes,raTemp,iNp,raaPrBdry,raaUserPressT,raPressLevels)
          DO iI=1,iNp
            iaaOpT(iOutTypes,iI)=iaOp(iI)
          END DO
          CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iJ,iNpmix, &
            iNumGases,iNp,iaNumLayer)
          CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases, &
            iNp,iaNpT,iaPrinterT,iNatm,iJ,iaGPMPAtmT, &
            raaOpT,raaUserPressT,iaNumLayer,iUpDown)
        END DO
      END IF   !iAtm < 0, iNp > 0
      ! &&&&&&&&&& END CASE 3A

      ! &&&&&&&&&& BEGIN CASE 3B
      IF (iNp < 0) THEN
        DO iJ = 1,iNatm
          IF (iaaRadLayer(iJ,1) < iaaRadLayer(iJ,2)) THEN
            iUpDown = 1            !upward travelling radiation
          ELSE
            iUpDown = -1           !downward travelling radiation
          END IF

          iaGPMPAtmT(iOutTypes) = iJ
          iaNpT(iOutTypes)      = iNp
          iaPrinterT(iOutTypes) = iPrinter
                        
          CALL AllLayersOutputPress(iaaRadLayer,iJ,iaNumLayer, &
            iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT, &
            raaPrBdry,raFracTop,raFracBot,raPressLevels)

          CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iJ,iNpmix, &
            iNumGases,iNp,iaNumLayer)
          CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases, &
            iNp,iaNpT,iaPrinterT,iNatm,iJ,iaGPMPAtmT, &
            raaOpT,raaUserPressT,iaNumLayer,iUpDown)
        END DO
      END IF   !iAtm < 0, iNp < 0
      ! &&&&&&&&&& END CASE 3B
        END IF     !overall iAtm < 0 if loop
    ! &&&&&&&&&& END CASE 3
    END DO       !do loop over iOut2

    iOutTypes=iOutTypes-1   !set it to the correct number actually found

! now set the elements of the temporary arrays/matrices to those that are
! sent out to the main program
    CALL SetActualPrintOptions(iOutTypes,iaPrinter,iaPrinterT, &
      iaGPMPAtm,iaGPMPAtmT,iaNp,iaNpT,iaaOp,iaaOpT,raaOp,raaOpT, &
      raaUserPress,raaUserPressT,iaNumlayer)

    RETURN
    end SUBROUTINE output4

!************************************************************************
! this subroutine "explodes" the input list for ONE mixed paths set/Gas ID set
    SUBROUTINE DoSetSome(iaOp,iList,iNp,iPrinter,iNumGases,iaGasID, &
    iMPSets,iGS)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iaGasID    = array containig GAS ID's
! iGS        = GAS ID or MP set number
! iList      = number of paths/mixed paths set in original iaOp
! iNp        = final number of paths/mixed paths set in iaOp for print option
! iaOp       = matrix with list of paths/MP/radiances output for each option
! iNumGases  = total number of gases read in from GASFIL + XSCFIL
! iMPSETs    = total number of mixed paths sets from *WEIGHT
! iPrinter   = printing option (1 or 2)
    INTEGER :: iPrinter,iNp,iaOp(kPathsOut),iMPSets,iNumGases,iList
    INTEGER :: iaGasID(kMaxGas),iGS
          
    INTEGER :: iI,iD,iaOpTemp(kPathsOut)

! make a working copy of iaOp
    iaOpTemp = iaOp

    IF (iPrinter == 2) THEN
      iD = iGS
      IF (iD > iMPSets) THEN
        write (kStdErr,*) 'oops! mixed path set iGS > those set in *WEIGHT'
        CALL DoSTOP
      END IF
    ELSE
      iI = 1
 11   CONTINUE
      IF (iGS /= iaGasID(iI)) THEN
        iI = iI+1
        GO TO 11
      END IF
      IF (iI <= iNumGases) THEN
        iD = iI
      ELSE
        write(kStdErr,*)'could not find gas ID iGS in list iaGASID'
        CALL DoSTOP
      END IF
    END IF

! heck that the input list read in from file has numbers <= kProfLayer
    DO iI = 1,iList
      IF ((iaOp(iI) > kProfLayer) .OR. (iaOp(iI) < 1)) THEN
        IF (iPrinter == 1) THEN
          write(kStdErr,*) 'path number : ',iaOp(iI)
          write(kStdErr,*)'Invalid path in list found in *OUTPUT'
        ELSEIF (iPrinter == 2) THEN
          write(kStdErr,*) 'mixed path number : ',iaOp(iI)
          write(kStdErr,*) 'Invalid mixedpath in list found in *OUTPUT'
        END IF
      CALL DoSTOP
      END IF
    END DO

! now that everything seems OK, go ahead and replicate list
! thus eg if iGS = 4 th gas, and the list contains 3 items 15 34 72 then we
! need the following final list 315 334 372
    DO iI = 1,iList
      iaOpTemp(iI) = iaOpTemp(iI)+(iD-1)*kProfLayer
    END DO

! reset iaOp
    iaOp = iaOpTemp

    RETURN
    end SUBROUTINE DoSetSome
!************************************************************************
! this subroutine "explodes" the input list for ALL mixed paths sets/Gas IDs
    SUBROUTINE DoSetAll(iaOp,iList,iNp,iPrinter,iNumGases,iMPSets)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iList      = number of paths/mixed paths set in original iaOp
! iNp        = final number of paths/mixed paths set in iaOp for print option
! iaOp       = matrix with list of paths/MP/radiances output for each option
! iNumGases  = total number of gases read in from GASFIL + XSCFIL
! iMPSETs    = total number of mixed paths sets from *WEIGHT
! iPrinter   = printing option (1 or 2)
    INTEGER :: iPrinter,iNp,iaOp(kPathsOut),iMPSets,iNumGases,iList
          
    INTEGER :: iI,iD,iaOpTemp(kPathsOut)

! make a working copy of iaOp
    iaOpTemp = iaOp

    IF (iPrinter == 1) THEN
      IF (iNp /= iList*iNumGases) THEN
        write(kStdErr,*)'bad *OUTPUT for paths output option -1 iN'
        write(kStdErr,*) '(iNp /= iList*iNumGases)'
        CALL DoSTOP
      END IF
      !check that the input list read in from file has numbers <= kProfLayer
      DO iI = 1,iList
        IF ((iaOp(iI) > kProfLayer) .OR. (iaOp(iI) < 1)) THEN
          write(kStdErr,*) 'path number : ',iaOp(iI)
          write(kStdErr,*) 'Invalid path in list found in *OUTPUT'
          write(kStdErr,*)'((iaOp(iI)>kProfLayer) or (iaOp(iI) <1)) '
          CALL DoSTOP
        END IF
      END DO
      !now that everything seems OK, go ahead and replicate list
      ! thus eg is there are 4 gases, and the list contains 3 items 15 34 72 then we
      ! need follwoing final list 15 34 72  115 134 172  215 234 272  315 334 372
      DO iD = 1,iNumGases
        DO iI = 1,iList
          iaOpTemp(iI+(iD-1)*iList) = iaOpTemp(iI)+(iD-1)*kProfLayer
        END DO
      END DO
    END IF

    IF (iPrinter == 2) THEN
      IF (iNp /= iList*iMPSets) THEN
        write(kStdErr,*)'bad *OUTPUT for MP output option -1 iN'
        write(kStdErr,*)'(iNp /= iList*iMPSets)'
        CALL DoSTOP
      END IF
      !check that the input list read in from file has numbers <= kProfLayer
      DO iI = 1,iList
        IF ((iaOp(iI) > kProfLayer) .OR. (iaOp(iI) < 1)) THEN
          write(kStdErr,*) 'Invalid MP in list found in *OUTPUT'
          write(kStdErr,*) '(iaOp(iI) > kProfLayer)or(iaOp(iI) < 1)'
          write(kStdErr,*) 'iI,iaOp(iI) = ',iI,iaOp(iI)
          CALL DoSTOP
        END IF
      END DO

      !now that everything seems OK, go ahead and replicate list
      !thus eg is there are 4 gases, and the list contains 3 items 15 34 72 then we
      !need follwoing final list 15 34 72  115 134 172  215 234 272  315 334 372
      DO iD = 1,iMPSets
        DO iI = 1,iList
          iaOpTemp(iI+(iD-1)*iList) = iaOpTemp(iI)+(iD-1)*kProfLayer
        END DO
      END DO
    END IF

! reset iaOp
    iaOp = iaOpTemp

    RETURN
    end SUBROUTINE DoSetAll
!************************************************************************
! this sets the actual print option variables from the temporary ones
    SUBROUTINE SetActualPrintOptions(iOutTypes,iaPrinter,iaPrinterT, &
    iaGPMPAtm,iaGPMPAtmT,iaNp,iaNpT,iaaOp,iaaOpT,raaOp,raaOpT, &
    raaUserPress,raaUserPressT,iaNumLayer)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! iaPrinterT = array containing the print options (1,2 or 3)
! iaGPMPAtmT   = array containing which atmospheres to be printed (if print
!              option ii=1 or 2, then iaGPMPAtm(ii)=0)
! iaNpT      = array with  the # of paths to be printed for each print option
! iOutTypes  = total number of outputs per each 25cm-1 run
    INTEGER :: iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
    INTEGER :: iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iaNp(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaaOpT(kMaxPrT,kPathsOut)
    REAL :: raaOp(kMaxPrint,kPathsOut),raaOpT(kMaxPrT,kPathsOut)
    REAL :: raaUserPress(kMaxPrint,kProfLayer)
    REAL :: raaUserPressT(kMaxPrT,kProfLayer)
    INTEGER :: iOutTypes,iaNumLayer(kMaxAtm)

    INTEGER :: iI,iJ,iM

    DO iI = 1,iOutTypes
      iaPrinter(iI) = iaPrinterT(iI)
      iaGPMPAtm(iI) = iaGPMPAtmT(iI)
      iaNp(iI) = iaNpT(iI)
    END DO

    DO iI = 1,iOutTypes
      DO iJ = 1,kPathsOut
        iaaOp(iI,iJ) = iaaOpT(iI,iJ)
      END DO
    END DO

    DO iI = 1,iOutTypes
      DO iJ = 1,kProfLayer
        raaOp(iI,iJ) = raaOpT(iI,iJ)
        raaUserPress(iI,iJ) = raaUserPressT(iI,iJ)
      END DO
    END DO

    write(kStdWarn,*) 'namelist file specified ',iOutTypes,' total outputs '
    IF (iOutTypes <= 0) THEN
      write(kStdErr,*) 'need to specify >= 1 output for kCARTA to run!'
      CALL DoStop
    END IF

    write(kStdWarn,*)'Output_number  OutOption(1/2/3)   Atm number'
    write(kStdWarn,*)'Num of path/MP/rads to output (if -1, list all)'
    write(kStdWarn,*)'Item      Type(1=path,2=MP,3=rad)      List'
    write(kStdWarn,*)'----------------------------------------------'
    DO ii = 1,iouttypes
      write(kStdWarn,30) ii,iaPrinter(iI),iaGPMPAtm(iI)
      write(kStdWarn,*) (iaaOp(ii,ij),ij=1,iaNp(ii))
      IF (iaPrinter(iI) == 3) THEN
        IF (iaNp(iI) == -1) THEN
          iM = iaNumLayer(iaGPMPAtm(iI))
        ELSE
          iM = iaNp(iI)
        END IF
        write(kStdWarn,*) 'list of pressures (partial fractions)'
        DO ij = 1,iM
          write(kStdWarn,*) raaUserPress(ii,ij),'(',raaOp(ii,ij),')'
        END DO
      END IF
      write(kStdWarn,*) ' '
    END DO
 30 FORMAT('     ',3(I5,'           '))

    RETURN
    end SUBROUTINE SetActualPrintOptions

!************************************************************************
! this subroutine parses in the line after finding print option3
    SUBROUTINE Read3(iAtm,iNp,iNpmix,iNatm,iNatm2,iListType, &
    iaGPMPAtmT,iaNpT,iOutTypes,iaPrinterT,iPrinter)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! INPUT is basically iAtm,iNp which are atm # and num radiances respectively
! iaPrinterT = array containing the print options (1,2 or 3)
! iaGPMPAtmT   = array containing which atmospheres to be printed (if print
!              option ii=1 or 2, then iaGPMPAtm(ii)=0)
! iaNpT      = array with  the # of paths to be printed for each print option
! iOutTypes  = total number of outputs per each 25cm-1 run
! iNpmix     = total number of ixed paths read in from MIXFIL
! iNatm      = number of atmospheres that RADFIL had in it
! iNatm2     = number of atmospheres that OUTPUT claims there should be
! iListType  = what type printing type (all or some layers)
    INTEGER :: iNatm,iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
    INTEGER :: iOutTypes,iNatm2,iNpmix,iListType
    INTEGER :: iPrinter,iAtm,iNp

    CHARACTER(7) :: caWord

    caWord = '*OUTPUT'

    IF (iNpmix < 0) THEN
      write(kStdErr,*)'Have not found *WEIGHTS ... cannot process'
      write(kStdErr,*)'make sure *WEIGHT is set before *OUTPUT'
      CALL DoSTOP
    END IF

    IF (iNatm == 0) THEN
      write(kStdErr,*)'Have not found *RADNCE ... cannot process'
      write(kStdErr,*)'make sure *RADNCE is before *OUTPUT'
      CALL DoSTOP
    END IF
     
    IF (iAtm == 0) THEN
      write(kStdErr,*) 'Cannot output spectra for atmosphere # 0 !!'
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF
    IF (iAtm > iNatm) THEN
      write(kStdErr,*) 'atmosphere # ',iAtm,' undefined in *RADNCE'
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF

    IF (iAtm > iNatm2) THEN
      iNatm2 = iAtm
    END IF
    IF (iNp < 0) THEN
      iListType = 1
    ELSE
      iListType = 2
    END IF
     
! do things one at a time!!!! ie even if iAtm < 0, think of ONE atmosphere!
    IF (iAtm > 0) THEN
      iaGPMPAtmT(iOutTypes) = iAtm
    ELSE
      iaGPMPAtmT(iOutTypes) = 1
    END IF

    iaNpT(iOutTypes) = iNp
    iaPrinterT(iOutTypes) = iPrinter
     
    RETURN
    end SUBROUTINE Read3

!************************************************************************
! this subroutine parses in the line after finding print option2
    SUBROUTINE Read2(iGS,iNp,iNpmix,iMPSets,iAtm, &
    iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,iNumGases)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! input is basically iGas and iNp which are MPset# and num MPs to be output
! iMPSets    = number of mixed path sets found in *WEIGHT
! iGS        = which GAS ID/mixed path set
! iaGPMPAtmT   = array containing which atmospheres to be printed (if print
!              option ii=1 or 2, then iaGPMPAtm(ii)=0)
! iaNpT      = array with  the # of paths to be printed for each print option
! iOutTypes  = total number of outputs per each 25cm-1 run
! iNpmix     = total number of ixed paths read in from MIXFIL
! iNatm      = number of atmospheres that RADFIL had in it
! iNatm2     = number of atmospheres that OUTPUT claims there should be
! iList      = how many paths/mixed paths to be output (-1 or > 0)
! iListType  = what type printing type (all or some layers)
! iAtm       = dummy for this (set to 0)
    INTEGER :: iAtm,iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT),iGS,iMPSets
    INTEGER :: iOutTypes,iNpmix,iListType,iList,iPrinter,iNp,iNumGases

    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead

    caWord = '*OUTPUT'

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
    5010 FORMAT('Error reading section ',A7,' of main file')

    iAtm = 0
    IF (iNpmix < 0) THEN
      write(kStdErr,*)'Have not found mixing table ... cannot process'
      write(kStdErr,*)'make sure *WEIGHT is set before *OUTPUT'
      CALL DoSTOP
    END IF
    iMPSets = iFloor(iNpmix*1.0/(kProfLayer))

    IF (iGS == 0) THEN
      write(kStdErr,*)'Cannot output spectra for mixed path set # 0 !'
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF
    IF (iGS > iMPSets) THEN
      write(kStdErr,*)'MP set # ',iGS,' undefined in *WEIGHT'
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF

! now recheck iNp and reset if necessary
    CALL ReCheck_iNp(iList,iListType,iGS,iNp,iMPSets,iNumGases,2)
    iaGPMPAtmT(iOutTypes) = iAtm     !=0 as NO radiances to be output
    iaNpT(iOutTypes) = iNp

    RETURN
    end SUBROUTINE Read2

!************************************************************************
! this subroutine parses in the line after finding print option1
    SUBROUTINE Read1(iGS,iNp,iNpmix,iMPSets,iAtm, &
    iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter, &
    iNumGases,iaGasID)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! input is basically iGas and iNp which are Gas# and numPaths to be output
! iMPSets    = number of mixed path sets found in *WEIGHT
! iGS        = which GAS ID/mixed path set
! iaGPMPAtmT   = array containing which atmospheres to be printed (if print
!              option ii=1 or 2, then iaGPMPAtm(ii)=0)
! iaNpT      = array with  the # of paths to be printed for each print option
! iOutTypes  = total number of outputs per each 25cm-1 run
! iNpmix     = total number of ixed paths read in from MIXFIL
! iNatm      = number of atmospheres that RADFIL had in it
! iNatm2     = number of atmospheres that OUTPUT claims there should be
! iList      = how many paths/mixed paths to be output (-1 or > 0)
! iListType  = what type printing type (all or some layers)
! iAtm       = dummy for this (set to 0)
    INTEGER :: iAtm,iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT),iGS,iMPSets
    INTEGER :: iOutTypes,iNpmix,iListType,iList,iPrinter,iNp
    INTEGER :: iNumgases,iaGasID(kMaxGas)

    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iI

    caWord = '*OUTPUT'

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7,' of main file')

    iNumLinesRead = 1

    iAtm = 0

    IF (iGS == 0) THEN
      write(kStdErr,*)'Cannot output spectra for GAS ID  = 0 !!'
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF
    IF (iGS > kMaxGas) THEN
      write(kStdErr,*)'Gas ID ',iGS,' undefined '
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF
    IF ((DoOutputLayer(iGS,iNumGases,iaGasID) < 0) .AND. (iGS > 0)) THEN
      write(kStdErr,*)'List of Gas IDs from MOLGAS/XSCGAS = '
      write(kStdErr,*)(iaGasID(iI),iI=1,iNumGases)
      write(kStdErr,*)'Gas ID ',iGS,' not set from *MOLGAS,*XSCGAS '
      write(kStdErr,*)'check *OUTPUT!!'
      CALL DoSTOP
    END IF

! now recheck iNp and reset if necessary
    CALL ReCheck_iNp(iList,iListType,iGS,iNp,iMPSets,iNumGases,1)
    iaGPMPAtmT(iOutTypes) = iAtm     ! = 0 as NO radiances to be output
    iaNpT(iOutTypes) = iNp

    RETURN
    end SUBROUTINE Read1

!************************************************************************
! this subroutine does the checking of iNp .. if necessary, it resets it
    SUBROUTINE ReCheck_iNp(iList,iListType,iGS,iNp, &
    iMPSets,iNumGases,iPrinter)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iList     = set how many items found in list
! iListType = 1,2,3,4    1 if -1  -1
!                        2 if -1  iL
!                        3 if iGS -1
!                        4 if iGS iL
! iGS       = GAS ID os Mixed path set
! iNp       = from inout file, how many gas paths/mixed paths to be output

    INTEGER :: iList,iListType,iGS,iNp,iMPSets,iNumGases,iPrinter

    INTEGER :: iI

    IF (iPrinter == 1) THEN
      iI = iNumGases
    ELSE
      iI = iMPSets
    END IF

    iList = 0                !set how many items to be found in list
    IF ((iGS == -1) .AND. (iNp == -1)) THEN
      iList = iNp
      iListType = 1
      iNp  =  -1             !for ALL gases/MP sets, output all layers
    ELSE IF ((iGS == -1) .AND. (iNp /= -1)) THEN
      iList = iNp
      iListType = 2
      iNp  =  iNp*iI         !for ALL gases/MP sets, output some layers
    ELSE IF ((iGS > 0) .AND. (iNp == -1)) THEN
      iList = iNp
      iListType = 3
      iNp  =  kProfLayer      !for one gas/MP set, output all layers
    ELSE IF ((iGS > 0) .AND. (iNp /= -1)) THEN
      iList = iNp
      iListType = 4
      iNp  =  iNp            !for one gas/MP set, output iNp layers
    END IF

    RETURN
    end SUBROUTINE ReCheck_iNp

!************************************************************************
! this subroutine takes care of the case iPrinter = 3 (radiance), iAtm > 0,
! (and so the user has a specific atmosphere in mind), iNp = -1
    SUBROUTINE AllLayersOutputPress(iaaRadLayer,iAtm,iaNumLayer, &
    iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT, &
    raaPrBdry,raFracTop,raFracBot,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! iaaRadLayer= list of layers in each radiating atmosphere
! iaaOpT     = matrix with list of paths/MP/radiances output for each option
! raaOpT     = matrix with fractional list of radiances output
! raaUserPressT= matrix with user specified pressures
! iaNumLayer = number of layers in each atmosphere, from *RADFIL
! iOutTypes  = total number of outputs per each 25cm-1 run
! raFracTop/Bot = how the top/bottom layers of mixing table have been affected
!                 by the mixing table
! raaPrBdry     = pressure boundaries of the atmospheres
! iAtm          = which atmoshere we are considering
! iaOp          = temporary part of iaaOpT
    REAL :: raPressLevels(kProfLayer+1)    !!!!actual pressures of levels
    INTEGER :: iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iAtm,iaOp(kPathsOut)
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raaOpT(kMaxPrT,kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)

    INTEGER :: iUpDown,iI,iLay

    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
      iUpDown = 1        !radiation going up
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,2)) THEN
      iUpDown = -1       !radiation going down
    END IF
! the following loop just does the number of layers in the relevant atmosphere
! taking care of the case iAtm > 0,iNp < 0
    DO iI = 2,iaNumLayer(iAtm)-1    !do complete layers
      iaOp(iI) = iI
      iaaOpT(iOutTypes,iI) = iI
      raaOpT(iOutTypes,iI) = 1.0  !define complete layer (fraction = 1.0)
      iLay = MP2Lay(iaaRadLayer(iAtm,iI))
      IF (iUpDown > 0) THEN
        raaUserPressT(iOutTypes,iI) = raPressLevels(iLay+1)
      ELSE
        raaUserPressT(iOutTypes,iI) = raPressLevels(iLay)
      END IF
    END DO
    IF (iUpDown > 0) THEN
      !special : bottom
      iI       = 1
      iaOp(iI) = iI
      iaaOpT(iOutTypes,iI) = iI
      iLay     = MP2Lay(iaaRadLayer(iAtm,iI))
      raaUserPressT(iOutTypes,iI) = raPressLevels(iLay+1)-delta
      !since we flip defns of partial layers for bottom ==> fraction <= 1.0
      raaOpT(iOutTypes,iI) = 1.0 !was incorrectly raFracBot(iAtm)
      raaOpT(iOutTypes,iI) = raFracBot(iAtm)

      !special : top
      iI       = iaNumLayer(iAtm)
      iaOp(iI) = iI
      iaaOpT(iOutTypes,iI) = iI
      raaUserPressT(iOutTypes,iI) = raaPrBdry(iAtm,2)
      raaOpT(iOutTypes,iI)        = raFracTop(iAtm)

    ELSE IF (iUpDown < 0) THEN
      !special : top
      iI       = 1
      iaOp(iI) = iI
      iaaOpT(iOutTypes,iI) = iI
      iLay     = MP2Lay(iaaRadLayer(iAtm,iI))
      raaUserPressT(iOutTypes,iI) = raPressLevels(iLay)
      !since we flip defns of partial layers for top ==> fraction = 1.0
      raaOpT(iOutTypes,iI) = 1.0 !was incorrectly raFracTop(iAtm)
      raaOpT(iOutTypes,iI) = raFracTop(iAtm)

      !special : bottom
      iI       = iaNumLayer(iAtm)
      iaOp(iI) = iI
      iaaOpT(iOutTypes,iI) = iI
      raaUserPressT(iOutTypes,iI) = raaPrBdry(iAtm,2)
      raaOpT(iOutTypes,iI) = raFracBot(iAtm)
    END IF

    RETURN
    end SUBROUTINE AllLayersOutputPress
!************************************************************************
! this subroutine checks to see iaOp is kosher
    SUBROUTINE SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iAtm,iNpmix, &
    iNumGases,iNp,iaNumLayer)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! iPrinter   = which is the current printing option
! iNp        = curreent number of paths/mixed paths/radiances to print
! iNpmix     = number of mixed paths
! iNumGases  = number of gases
! iaaOpT     = matrix with list of paths/MP/radiances output for each option
! iaNumLayer = number of layers in each atmosphere, from *RADFIL
! iOutTypes  = total number of outputs per each 25cm-1 run
! iAtm          = which atmoshere we are considering
! iaOp          = temporary part of iaaOpT
    INTEGER :: iaaOpT(kMaxPrT,kPathsOut),iaOp(kPathsOut),iOutTypes
    INTEGER :: iPrinter,iAtm,iNpmix,iNumGases,iNp,iaNumLayer(kMaxAtm)

    INTEGER :: iI

    DO iI = 1,iNp
      iaaOpT(iOutTypes,iI) = iaOp(iI)
      IF ((iPrinter == 3) .AND. (iAtm > 0)) THEN
        !all numbers here should be between 1 and iaNumLayer(relevant atm)
        IF ((iaOp(iI) < 1) .OR. (iaOp(iI) > iaNumlayer(iAtm))) THEN
          write(kStdErr,*)'Printing option selected = ',iPrinter
          write(kStdErr,*)'layer # ',iI,' = ',iaOp(iI)
          write(kStdErr,*)'Error in reading *output section'
          write(kStdErr,*)'Radiating MP must be between 1 and ',iaNumLayer(iAtm)
          CALL DoSTOP
        END IF
      END IF
      IF (iPrinter == 2) THEN
        ! all numbers here should be between 1 and iNpmix
        IF ((iaOp(iI) < 1) .OR. (iaOp(iI) > iNpmix)) THEN
          write(kStdErr,*)'Printing option selected = ',iPrinter
          write(kStdErr,*)'mixed path # ',iI,' = ',iaOp(iI)
          write(kStdErr,*)'Error in reading *output section'
          write(kStdErr,*)'Mixed Path number must be between 1 and ',iNpmix
          CALL DoSTOP
        END IF
      END IF
      IF ((iPrinter == 1)) THEN
        ! all numbers here should be between 1 and kProfLayer*iNumGases
        IF ((iaOp(iI) < 1) .OR. (iaOp(iI) > kProfLayer*iNumGases)) THEN
          write(kStdErr,*)'Printing option selected = ',iPrinter
          write(kStdErr,*)'# of gases in GASFIL+XSCFIL = ',iNumGases
          write(kStdErr,*)'path # ',iI,' = ',iaOp(iI)
          write(kStdErr,*)'Error in reading *output section'
          write(kStdErr,*)'Path number must be between 1 and ',kProfLayer*iNumGases
          CALL DoSTOP
        END IF
      END IF
    END DO

    RETURN
    end SUBROUTINE SetIaaOpT
!************************************************************************
! this subroutine checks to see if iPrinter=3 has already been stored
! if so then it supplements the previous relevant vector of iaaOpT,
! raaOpT,raaUserPressT with the new information
! this is only to check ONE atmosphere being added on
    SUBROUTINE CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases, &
    iNp,iaNpT,iaPrinterT,iNatm,iAtm,iaGPMPAtmT, &
    raaOpT,raaUserPressT,iaNumLayer,iUD)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! iUD        = is radiation travelling up or down
! iaNumLayer = number of radiating layers in the defined atmospheres
! iNatm      = number of atmospheres read in from *RADNCE
! iAtm       = which atmosphere(s) this print option is concerned with
! iaGPMPAtmT   = array containig atmopsphere printing options
! raaOpT     = matrix list of fractions to output radiances at
! raaOpT     = matrix list of pressures to output radiances at
! iPrinter   = which is the current printing option
! iNp        = curreent number of paths/mixed paths/radiances to print
! iaNpT      = array containig number of paths/mixed paths/radiances to print
! iaPrinterT = array containig printing options
! iNpmix     = number of mixed paths
! iNumGases  = number of gases
! iaaOpT     = matrix with list of paths/MP/radiances output for each option
! iOutTypes  = total number of outputs per each 25cm-1 run
    INTEGER :: iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iaPrinterT(kMaxPrT)
    INTEGER :: iPrinter,iNpmix,iNumGases,iNp,iaNpT(kMaxPrT),iUD
    INTEGER :: iAtm,iNatm,iaGPMPAtmT(kMaxPrT),iaNumLayer(kMaxAtm)
    REAL :: raaOpT(kMaxPrT,kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)
        
    INTEGER :: iI,iJ,iJ1,iFound,iaOp1(kPathsOut),iaOp2(kPathsOut),i12
    REAL :: ra1(kProfLayer),ra2(kProfLayer),raP1(kProfLayer),raP2(kProfLayer)

! heck to see if any of the previous print options  == iPrinter,iAtm
    iFound = -1
    DO iI = 1,iOutTypes-1
      IF ((iaPrinterT(iI) == iPrinter) .AND. (iaGPMPAtmT(iI) == iAtm))THEN
        iFound = iI
      END IF
    END DO

    IF (iFound > 0) THEN    !have found a match!!! have to merge together
      !irst check to see if one or the other option wants ALL paths output
      IF (iaNpT(iFound) == -1) THEN
        iaNpT(iFound) = -1      !old option included ALL layers
        iJ1=iFound
      END IF
      IF (iaNpT(iOutTypes) == -1) THEN
        iaNpT(iFound) = -1      !new option includes ALL layers
        iJ1=iOutTypes
      END IF
      IF (iaNpT(iFound) == -1) THEN
        iJ = iaNumLayer(iAtm)
        DO iI = 1,iJ
          iaaOpT(iFound,iI) = iI
          raaOpT(iFound,iI) = raaOpT(iJ1,iI)
          raaUserPressT(iFound,iI) = raaUserPressT(iJ1,iI)
        END DO
      END IF

      IF ((iaNpT(iOutTypes) /= -1) .AND. (iaNpT(iFound) /= -1)) THEN
        !have to slug it out! ah well!!!!
        DO iI = 1,iaNpT(iFound)
          iaOp1(iI) = iaaOpT(iFound,iI)
          ra1(iI) = raaOpT(iFound,iI)
          raP1(iI) = raaUserPressT(iFound,iI)
        END DO
        DO iI = 1,iaNpT(iOutTypes)
          iaOp2(iI) = iaaOpT(iOutTypes,iI)
          ra2(iI) = raaOpT(iOutTypes,iI)
          raP2(iI) = raaUserPressT(iOutTypes,iI)
        END DO
        CALL MergeRads(iaOp1,iaOp2,ra1,ra2,raP1,raP2, &
            iaNpT(iFound),iaNpT(iOutTypes),i12,iUD)
        !now update raaOpT,raaUserPressT,iaaOpT,iaNpT
        iaNpT(iFound) = i12
        DO iI = 1,i12
          iaaOpT(iFound,iI) = iaOp1(iI)
          raaOpT(iFound,iI) = ra1(iI)
          raaUserPressT(iFound,iI) = raP1(iI)
        END DO
      END IF
    END IF

    IF (iFound < 0) THEN !we really have found a new print option
      iOutTypes = iOutTypes+1
    END IF

    IF (iOutTypes > kMaxPrT) THEN
      write(kStdErr,*)'iOutTypes,kMaxPrT = ',iOutTypes,kMaxPrT
      write(kStdErr,*)'Cannot store so many printing options!!!'
      write(kStdErr,*)'Either increase kMaxPrint in kcartaparam.f90, or '
      write(kStdErr,*)'decrease # of print options found in *OUTPUT'
      CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE CheckRaaOpT

!************************************************************************
! this subroutine checks to see if iPrinter=1 or 2 has already been stored
! if so then it supplements the previous relevant vector of iaaOpT with the
! new information
    SUBROUTINE CheckIaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix, &
    iNumGases,iNp,iaNpT,iaPrinterT)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! iPrinter   = which is the current printing option
! iNp        = curreent number of paths/mixed paths/radiances to print
! iaNpT      = array containig number of paths/mixed paths/radiances to print
! iaPrinterT = array containig printing options
! iNpmix     = number of mixed paths
! iNumGases  = number of gases
! iaaOpT      = matrix with list of paths/MP/radiances output for each option
! iOutTypes  = total number of outputs per each 25cm-1 run
    INTEGER :: iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iaPrinterT(kMaxPrT)
    INTEGER :: iPrinter,iNpmix,iNumGases,iNp,iaNpT(kMaxPrT)

    INTEGER :: iI,iJ,iFound,iaOp1(kPathsOut),iaOp2(kPathsOut)

! heck to see if any of the previous print options  == iPrinter
    iFound = -1
    DO iI = 1,iOutTypes-1
      IF (iaPrinterT(iI) == iPrinter) THEN
        iFound  =  iI
      END IF
    END DO

    IF (iFound > 0) THEN    !have found a match!!! have to merge together
      !first check to see if one or the other option wants ALL paths output
      IF (iaNpT(iFound) == -1) THEN
        iaNpT(iFound) = -1      !old option included ALL paths/mixed paths
      END IF
      IF (iaNpT(iOutTypes) == -1) THEN
        iaNpT(iFound) = -1      !new option includes ALL paths/mixed paths
      END IF
      IF (iaNpT(iFound) == -1) THEN
        IF (iPrinter == 1) THEN
          iJ = iNumGases*kProfLayer
        ELSE
          iJ = iNpMix
        END IF
        DO iI = 1,iJ
          iaaOpT(iFound,iI) = iI
        END DO
      END IF

      IF ((iaNpT(iOutTypes) /= -1) .AND. (iaNpT(iFound) /= -1)) THEN
        !have to slug it out! ah well!!!!
        DO iI = 1,iaNpT(iFound)
          iaOp1(iI) = iaaOpT(iFound,iI)
        END DO
        DO iI = 1,iaNpT(iOutTypes)
          iaOp2(iI) = iaaOpT(iOutTypes,iI)
        END DO
        CALL MergeLists(iaOp1,iaOp2,iaNpT(iFound),iaNpT(iOutTypes),iI)
        !now update iaaOpT,iaNpT
        iaNpT(iFound) = iI
        DO iI = 1,iaNpT(iFound)
          iaaOpT(iFound,iI) = iaOp1(iI)
        END DO
      END IF
    END IF

    IF (iFound < 0) THEN !we really have found a new print option
      iOutTypes = iOutTypes+1
    END IF

    IF (iOutTypes > kMaxPrT) THEN
      write(kStdErr,*)'iOutTypes,kMaxPrT = ',iOutTypes,kMaxPrT
      write(kStdErr,*)'Cannot store so many printing options!!!'
      write(kStdErr,*)'Either increase kMaxPrint in kcartaparam.f90, or'
      write(kStdErr,*)'decrease # of print options found in *OUTPUT'
      CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE CheckIaaOpT
!************************************************************************
! this subroutine takes two lists, merges them together. If duplicate items
! are found in the second list, they are discarded. The new list is then
! sorted and stored in list1. The total length of the new list1 is iSum
    SUBROUTINE MergeLists(iaOp1,iaOp2,i1,i2,iSum)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iaOp1(kPathsOut),iaOp2(kPathsOut),i1,i2,iSum

    INTEGER :: iI,iaRes(kPathsOut),iOK

! irst put in elements of iaOp1 into iaRes
    DO iI = 1,i1
      iaRes(iI) = iaOp1(iI)
    END DO
    iSum = i1

! ow put in elements of iaOp2 into iaRes, ensuring no duplicates
    DO iI = 1,i2
      iOK = DoOutputLayer(iaOp2(iI),iSum,iaRes)
      IF (iOK < 0) THEN   !iaOp2(iI) not found .. add to list
        iSum = iSum+1
        iaRes(iSum) = iaOp2(iI)
        CALL dosort(iaRes,iSum)
      ELSE
        write(kStdWarn,*) 'following is duplicate : ',iaOp2(iI)
      END IF
    END DO

! ort iaRes
    CALL dosort(iaRes,iSum)

! pdate iaOp1
    DO iI = 1,iSum
      iaOp1(iI) = iaRes(iI)
    END DO

    RETURN
    end SUBROUTINE MergeLists

!************************************************************************
! this subroutine takes 6 lists, merges them together. If duplicate items
! are found in the second list, they are discarded. The new list is then
! sorted and stored in list1. The total length of the new list1 is iSum
! the merging is done based elememts in pressure list raP1,raP2
    SUBROUTINE MergeRads(ia1,ia2,ra1,ra2,raP1,raP2,i1,i2,iSum,iUD)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! i1,i2     = lenght of the lits blah1,blah2
! iSum      = final length of combined list .. stored in ia1,ra1,raP1
! iUD       = direction of radiation travel
! ra1,ra2   = arrays containing layer fractions
! raP1,raP2 = arrays containing pressure levels
! ia1,ia2   = arrays containing layers
    INTEGER :: ia1(kPathsOut),ia2(kPathsOut),i1,i2,iSum,iUD
    REAL :: raP1(kProfLayer),raP2(kProfLayer),ra1(kProfLayer),ra2(kProfLayer)

    INTEGER :: iI,iaRes(kPathsOut),iNew
    REAL :: raRes(kProfLayer),raPRes(kProfLayer)

! irst put in elements of ia1 into iaRes
    DO iI = 1,i1
      iaRes(iI) = ia1(iI)
      raRes(iI) = ra1(iI)
      raPRes(iI) = raP1(iI)
    END DO
    iSum = i1

! ow put in elements of ia2 into iaRes, ensuring no duplicates
    DO iI = 1,i2
      CALL DoMergePress(raP2(iI),raPRes,ia2(iI),iaRes,ra2(iI),raRes,iSum,iUD,iNew)
      IF (iNew < 0) THEN   !raP2(iI) already present in array ..
        write(kStdWarn,*) 'following is duplicate : ',raP2(iI)
      END IF
    END DO

    CALL DoSortPress(iaRes,raRes,raPRes,iSum,-iUD)

!update ia1
    DO iI = 1,iSum
      ia1(iI) = iaRes(iI)
      ra1(iI) = raRes(iI)
      raP1(iI) = raPRes(iI)
    END DO

    RETURN
    end SUBROUTINE MergeRads

!************************************************************************
! this subroutine checks to see if we insert pressure rP into array raPRes
! if so then it also inserts some other elements into the respective arrays
    SUBROUTINE DoMergePress(rP,raPRes,i2,iaRes,r2,raRes,iSum, &
    iUD,iNew)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNew   = is this a new pressure?If so iNew = +1
! iUD    = direction of radiation travel
! rP     = see if this is a new pressure
! raPRes = this is the array we check to see if rP exists inside
! i2     = this is the pressure layer assocoated with rP
! iaRes  = this is array containing the pressure layer numbers to be output
! r2     = this is the fraction associated with rP
! raRes  = this is the array containing the fractions of the press. layers
!          to be output
! iSum   = how far down the arrays to search
    REAL :: rP,r2
    INTEGER :: iaRes(kPathsOut),i2,iNew,iSum,iUD
    REAL :: raRes(kProfLayer),raPRes(kProfLayer)

    INTEGER :: iI

! recall the pressures are stored in in(de)creasing order, depending on iUD
! but also note that as subroutine DoOutptuRadianec checks to see if the LAYER
! being processed has an output pressure level, and it goes thru the WHOLE
! list while doing the check, we do not REALLY have to set everything in
! increasing or decreasing order (ie iUD not really needed)
    iNew = 1   !assume new pressure
    iI = 1
 10 CONTINUE
    IF (abs(raPRes(iI)-rP) <= delta) THEN !this pressure was there before
      iNew = -1
    END IF
    IF ((iI < iSum) .AND. (iNew > 0)) THEN
      ! till not found pressure and have some more elements to check in array
      iI = iI+1
      GO TO 10
    END IF

    IF (iNew > 0) THEN     !pressure not found .. add it on
      IF (iSum < kProfLayer) THEN !we add on relevant elements to arrays
        iSum = iSum+1
        raPRes(iSum) = rP
        iaRes(iSum) = i2
        raRes(iSum) = r2
      ELSE IF (iSum == kProfLayer) THEN !we cannot add it on
        write(kStdErr,*)'cannot merge press',rP,'into previous list'
        CALL DoSTOP
      END IF
    END IF

    RETURN
    end SUBROUTINE DoMergePress

!************************************************************************
! this subroutine converts the pressures found in *OUTPUT to which
! layer this is in the radiating atmosphere
! since DISORT and RTSPEC would like to output stuff at layer boundaries, we
! have to round the pressures so that they coincide with layer boundaries
    SUBROUTINE PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iAtm, &
    iOutTypes,raPress,iNp,raaPrBdry,raaUserPressT,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: kMaxPrT
    PARAMETER (kMaxPrT=kMaxPrint+1)

! raaOpT     = matrix with fractional list of layers where radiances output
! iaNumLayer = number of layers in each atmosphere, from *RADFIL
! iaaRadLayer= list of layers in each atmosphere, from *RADFIL
! raaPrBdry  = matrix containing start/stop pressures for each atmosphere
! iAtm       = which atm
! iOutTypes  = which output option this is
! raaUserPressT= which pressures to output the sutff at

! raPress     = array containing pressure levels where radiances to be output
! iNp        = number of pressures where radiances are to be output
! iaOp       = this subroutine figures out which layers are to be output
    REAL :: raPressLevels(kProfLayer+1)    !!!!actual pressures of levels
    INTEGER :: iaOp(kPathsOut),iAtm,iNp,iOutTypes
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raaOpT(kMaxPrT,kProfLayer),raaPrBdry(kMaxAtm,2)
    REAL :: raPress(kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)

    REAL :: rP,rF,raPressDiff(kProfLayer+1),rPressDiffMin,rPScatter,rPScatter1
    INTEGER :: iI,i1,i2,iOne,iFound,iDirection,iW,iOffSet
    INTEGER :: iStartLay,iSilly,iSillySave

! check the radiation direction
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
      iDirection = 1     !going up
    ELSE
      iDirection = -1    !going down
    END IF

! find which set of 100 mixed paths the layers for this atmosphere are from
! 1-100 = set 1, 101-200 = set 2, 201-300 = set 3 etc
! note : the reason i use iaaRadLayer(iAtm,2) is as follows : suppose we have
! uplooking instr, defined between mixed paths 100 --> 1  ... if we used
!         iW == iFloor((iaaRadLayer(iAtm,1)*1.0)/(kProfLayer*1.0))+1
!            == iFloor(100.0/100.0)+1=2  instead of 1

    iW = iFloor((iaaRadLayer(iAtm,2)*1.0)/(kProfLayer*1.0))+1
    iOffSet = (iW-1)*kProfLayer
          
! set which radiating pressure layer atmosphere starts at, to compare to iI
    iStartLay = MP2Lay(iaaRadLayer(iAtm,1))

    DO iI = 1,iNp

      ! see if this pressure lies within the start/stop defined in raaPrBdry
      ! if not, then redefine
      rP = raPress(iI)
      IF (iDirection > 0) THEN
        IF (rP < 0) THEN
          rP = 0.0
          rP = raaPrBdry(iAtm,2)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to TOA',rP
        END IF
        IF (rP > raaPrBdry(iAtm,1)) THEN
          rP = raaPrBdry(iAtm,1)-delta*1.001
          rP = raaPrBdry(iAtm,1)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to GND',rP
        END IF
        IF (rP < raaPrBdry(iAtm,2)) THEN
          rP = raaPrBdry(iAtm,2)+delta*1.001
          rP = raaPrBdry(iAtm,2)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to TOA',rP
        END IF
      END IF
      IF (iDirection < 0) THEN
        IF (rP < 0) THEN
          rP = 1200.0
          rP = raaPrBdry(iAtm,2)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to GND',rP
        END IF
        IF (rP < raaPrBdry(iAtm,1)) THEN
          rP = raaPrBdry(iAtm,1)+delta*1.001
          rP = raaPrBdry(iAtm,1)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to TOA',rP
        END IF
        IF (rP > raaPrBdry(iAtm,2)) THEN
          rP = raaPrBdry(iAtm,2)-delta*1.001
          rP = raaPrBdry(iAtm,2)
          write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
          write(kStdWarn,*) raPress(iI),' to GND',rP
        END IF
      END IF

      ! if DISORT or RTSPEC, as a first cut, find out which pressure level rP
      ! lies closest to. we really HAVE to be careful about up and downlook stuff
      IF ((kWhichScatterCode /= 0) .AND. &
         ((kWhichScatterCode /= 1) .OR. (kWhichScatterCode /= 4) .OR. &
         (kWhichScatterCode /= 5))) THEN
        rPressDiffMin = +1.0e10
        DO i1 = 1,kProfLayer + 1
          raPressDiff(i1) = abs(rP - raPressLevels(i1))
        END DO
        DO i1 = 1,kProfLayer + 1
          IF (raPressDiff(i1) < rPressDiffMin) THEN
            rPressDiffMin =  raPressDiff(i1)
            i2 = i1
          END IF
          END DO
          rPScatter  = rP

          iSilly = 1
 111      CONTINUE
          IF ((rP <= raPressLevels(iSilly)) .AND. &
            (rP >= raPressLevels(iSilly+1))) THEN
            GOTO 121
          ELSE
            iSilly = iSilly + 1
            GOTO 111
          END IF
 121      CONTINUE
          IF (iDirection < 0) THEN
            rPScatter1 = raPressLevels(iSilly) !uplook, so reset to layer bottom
            write(kStdWarn,*) 'low level of output layer is ',rPScatter1, 'mb'
          ELSEIF (iDirection > 0) THEN
            rPScatter1 = raPressLevels(iSilly+1) !downlook, reset to layer top
            write(kStdWarn,*) 'top level of output layer is ',rPScatter1, 'mb'
          END IF

          IF ((i2 /= iSilly) .OR. (abs(rPscatter-rPScatter1) >= 1.0e-3)) THEN
            rPScatter = rPScatter1
          END IF

          ! having saved this in variable rPscatter, make sure you are still w/in atm!!!
          IF (iDirection < 0) THEN  !!! you probably want P ~ 1200 mb
            !!!! raaPrBdry(iAtm,1) = TOA, raaPrBdry(iAtm,2) = GND
            !!!! check that things are ok for downward dir to instr at GND
            IF (rPScatter > raaPrBdry(iAtm,2)) THEN
              rPScatter = raaPrBdry(iAtm,2)
            END IF
            IF (rPScatter < raaPrBdry(iAtm,1)) THEN
              rPScatter = raaPrBdry(iAtm,1)
            END IF
          END IF
          IF (iDirection > 0) THEN
            !!!! raaPrBdry(iAtm,1) = GND, raaPrBdry(iAtm,2) = TOA
            !!!! check that things are ok for upward direction to intr at TOA
            IF (rPScatter > raaPrBdry(iAtm,1)) THEN
              rPScatter = raaPrBdry(iAtm,1)
            END IF
            IF (rPScatter < raaPrBdry(iAtm,2)) THEN
              rPScatter = raaPrBdry(iAtm,2)
            END IF
          END IF
          rP = rPScatter
          write (kStdWarn,*) 'for RTSPEC/DISORT, reset pressure to ',rP
        END IF
                     
        ! now find out the pressure layer this pressure lies within
        iFound = -1
        i1     = 1
        i2     = 2
        10 CONTINUE
        IF ((rP <= raPressLevels(i1)) .AND. (rP > raPressLevels(i2))) THEN
          iFound = 1
          iOne = i1+(iW-1)*kProfLayer  !set this correctly as Mixed Path number
        END IF
        IF ((iFound < 0) .AND. (i1 < kProfLayer)) THEN
          i1 = i1+1
          i2 = i2+1
          GO TO 10
        END IF
        IF ((iFound < 0)) THEN
          IF (abs(rP-raPressLevels(kProfLayer+1)) <= delta) THEN
            i1     = kProfLayer
            iFound = 1
            iOne = i1+(iW-1)*kProfLayer  !set this correctly as Mix Path number
          ELSE
            write(kStdErr,*) 'atm#',iAtm,' could not find pressure ',rP
            write(kStdErr,*) 'within KLAYERS pressure levels. Check'
            write(kStdErr,*) '*RADNCE and *OUTPUT sections'
            CALL DoSTOP
          END IF
        END IF

    ! now see which radiating layer this corresponds to in defined atmosphere
    ! eg if defined atmospghere was 3 950.0  0.15
    ! then 950.0 is within 6th pressure layer, 0.15 within 96th pressure layer
    ! and the 3 implies this is the third group of 100 layer weigths
    ! ==> this atmosphere consists of the 91 mixed paths
    !       206,207,208, ... 294,295,296
    ! with mixed paths 206 and 296 having fractional weights
    
    ! so now if the pressure at which we want radiance to be output is 151.2664,
    ! this is the 50th AIRS pressure level, but it would be the (250-206+1)=45th
    ! radiating layer in the defined atmosphere. Ditto upward looking instrument

    ! we know rP has been checked to lie within the pressure boundaries of
    ! the current atm, so if this is bottom or topmost layer, cannot exceed
    ! rFracBot or rFracTop respectively
        IF (iDirection > 0) THEN   !radiation goes up
          iaOp(iI) = iOne-iaaRadLayer(iAtm,1)+1
          !want BOTTOM fraction of layer
          rF = (raPressLevels(i1)-rP)/(raPressLevels(i1)-raPressLevels(i2))
          IF (iOne == (iaaRadLayer(iAtm,1))) THEN
            !want TOP fraction of  lowest layer
            !print *,'aha : bottom layer for down look instr!!!!!!'
            rF = 1-rF
          END IF
          IF (abs(rF-1.00) <= delta) THEN
            rF = 1.00
          END IF
          raaUserPressT(iOutTypes,iI) = rP
          raaOpT(iOutTypes,iI) = rF
          write(kStdWarn,*) 'iDirection,rP,rF = ',iDirection,rP,rF

        ELSE IF (iDirection < 0) THEN !radiation goes up
          iaOp(iI) = iaaRadLayer(iAtm,1)-iOne+1
          !want BOTTOM fraction of layer
          rF = (rP-raPressLevels(i2))/(raPressLevels(i1)-raPressLevels(i2))
          IF (iOne == (iaaRadLayer(iAtm,1))) THEN
            !want BOTTOM  fraction of highest layer
            !print *,'aha : top layer for up look instr!!!!!!'
            rF = 1-rF
          END IF
          IF (abs(rF-1.00) <= delta) THEN
            rF = 1.00
          END IF
          raaUserPressT(iOutTypes,iI) = rP
          raaOpT(iOutTypes,iI) = rF
          write(kStdWarn,*) 'iDirection,rP,rF = ',iDirection,rP,rF
        END IF

        IF ((iaOp(iI) > iaNumLayer(iAtm)) .OR. (iaOp(iI) < 1))THEN
          write(kStdErr,*) 'iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,End),iNumLayer,iOne,='
          write(kStdErr,*) iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iaNumlayer(iAtm)),iaNumlayer(iAtm),iOne
          write(kStdErr,*) 'iI,iaOp(iI) = ',iI,iaOp(iI)
          write(kStdErr,*) 'Pressure ',rP,'invalid for atm# ',iAtm, 'see press levels below : '
	  DO i2 = 1,kProfLayer
	    write(kStdErr,*) i2,raPresslevels(i2)
	  END DO
          CALL DoSTOP	    
        END IF

        IF (raaOpT(iOutTypes,iI) < 0.0) THEN
          raaOpT(iOutTypes,iI) = 0.0
        END IF
        IF (raaOpT(iOutTypes,iI) > 1.0) THEN
          raaOpT(iOutTypes,iI) = 1.0
        END IF
         
    END DO

    RETURN
    end SUBROUTINE PressTOLayers

!************************************************************************
END MODULE n_output
