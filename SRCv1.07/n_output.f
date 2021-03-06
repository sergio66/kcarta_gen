c Copyright 2000
c University of Maryland Baltimore County 
c All Rights Reserved

c this is the OUTPUT section
c are ALWAYS defined wrt to a FULL layer (the pressure level boundaries of 
c which are defined in NewRefProfiles/outincLAY.param)

c the following two paragraphs are for upgoing radiation (down look instr)

c eg     if top AIRS layer is from 10 mb to 0.005 mb 
c    and if pressure where aircraft flies is  5.0 mb
c then only fraction rF=(10-5)/(10-0.005) ~ 1/2 = rFracTop of the top 
c layer is used .... 
c    --------------
c
c    //////////////          use this LOWER portion
c    --------------c
c so if the Print Option for this atm is
c 3
c 1 1
c 0
c implying that one radiance, at aircraft height, is to be output, then
c raaOp(iAtm,1) = 0.5 (=rFracTop)
c if the user had said he wanted a radiance output at pressure 7.5 mb : 
c 3
c 1 1
c 7
c then raaOp(iAtm,1) = (10-7)/(10-0.005) ~ 0.25

c note that we have to be very careful about the bottommost layer
c eg     if bottom AIRS layer is from 1000 mb to 900 mb 
c    and if surface pressure  is 950 mb
c then only fraction rF=(950-900)/(1000-900) ~ 1/2 = rFracBot of the bottom
c layer is used
c    --------------
c    //////////////          use this UPPER portion
c    //////////////
c
c    --------------
c so if the Print Option for this atm is
c 3
c 1 1
c 925
c implying that one radiance, at 925 mb, is to be output, then
c raaOp(iAtm,1) = (925-900)/(1000-9000) = 0.25

c also note that if atmosphere is defined between eg 1000 to 0.005 mb, and the
c user says to output radiance at 1000 mb, then this means that the radiation
c at the surface is output : 
c raRad(vu) = raEms(vu) * ttorad(Tsurf,vu) + thermal(vu) + solar(vu)

c also note if user specifies radiances to be output at EACH layer 
c 3
c 1 -1
c then for each layer, the radiation at uppermost part of the layer is output
c ie bottom layer : surface to top of layer           --> output rad
c      next layer : bottom of layer to top of layer   --> output rad 
c        ....          .......                      ...
c       top layer : bottom of layer to aircraft posn  --> output rad 

c @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


c SAME CONSIDERATIONS APPLY TO UPWARD LOOKING INSTRUMENT .. have to be very
c careful about topmost layer

c also note that if atmosphere is defined between eg 0.005 to 10000 mb, and 
c user says to output radiance at 0.005 mb, then this means that the radiation
c at the very TopOfAtmosphere is output : 
c raRad(vu) = ttorad(2.73,vu) + solar(vu)

c also note if user specifies radiances to be output at EACH layer 
c 3
c 1 -1
c then for each layer, the radiation at bottommost part of the layer is output
c ie    top layer : TOA  to bottom of layer           --> output rad
c      next layer : top of layer to bottom of layer   --> output rad 
c        ....          .......                      ...
c    bottom layer : top of layer to instr posn        --> output rad 

c************************************************************************

c this subroutine deals with the 'OUTPUT' keyword
c relevant output is the printing switch, iPrinter
c if iDat = 3 then read       iAtm==atmospher to be output,
c                             iNp==no of mixed paths to be output
c ELSE if iDat =1,2 then read iWhich==which gas or miexd set to be output
c                             iNp==no of mixed pathes to be output

c option 1 : gas paths
c if user wants to output spectra for ALL gases, ALL 100 layers, file says
c      -1 -1
c else if user wants to output spectra for GasID iG, all 100 layers, file says
c      iG  -1
c else if user wants to output spectra for GasID iG, iN layers, file says
c      iG  iN
c      followed by list of iN layers
c else if user  wants to output spctra for ALL gases, iN layers, file says
c      -1  iN
c      followed by list of iN layers

c option 2 : mixed paths
c if user wants to output spectra for ALL MP sets, ALL 100 layers, file says
c      -1 -1
c else if user wants to output spectra for set iS, all 100 layers, file says
c      iS  -1
c else if user wants to output spectra for set iS, iN layers, file says
c      iS  iN
c      followed by list of iN layers
c else if user  wants to output spctra for ALL sets, iN layers, file says
c      -1  iN
c      followed by list of iN layers

c option 3 : radiances
c if user wants to output radiances for ALL atm, ALL layers, input file says 
c       -1 -1
c else if wants to output radiances for Nth atmosphere, all layers, file says
c        N  -1
c else if wants to output radiances for ALL atmosphere, at L presss, file says
c       -1  L
c       followed by list of L pressures
c hence always read in 2 integers

      SUBROUTINE output4(iaPrinter,iaGPMPAtm,iaNp,
     $           iaaOp,raaOp,raaUserPress,
     $           iaNumLayer,iaaRadLayer,iNatm,iNatm2,
     $           iOutTypes,iNumGases,iNpmix,
     $           raFracTop,raFracBot,raaPrBdry,
     $           iaGases,caComment)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c caComment  = 80 character user given comment
c iaGases    = array containing info about GAS ID's fouund in MOL/XSCGAS
c iaPrinter  = array containing the print options (1,2 or 3)
c iaGPMPAtm    = array containing which atmospheres to be printed (if print
c              option ii=1 or 2, then iaGPMPAtm(ii)=0)
c iaNp       = array with  the # of paths to be printed for each print option
c iaaOp      = matrix with list of paths/MP/radiances output for each option
c raaOp      = matrix with fractional list of where radiances should be output
c raaUserPress= matrix with user specified pressures
c iaNumLayer = number of layers in each atmosphere, from *RADFIL
c iaaRadLayer= list of layers in each atmosphere, from *RADFIL
c iOutTypes  = total number of outputs per each 25cm-1 run
c iNumGases  = total number of gases read in from GASFIL + XSCFIL
c iNpmix     = total number of ixed paths read in from MIXFIL
c iErr       = flags an error (mainly file I/O)
c iNatm      = number of atmospheres that RADFIL had in it
c iNatm2     = number of atmospheres that OUTPUT claims there should be
c raFracTop/Bot = how the top/bottom layers of mixing table have been affected
c                 by the mixing table
c raaPrBdry     = pressure boundaries of the atmospheres
      INTEGER iNatm,iaGases(kMaxGas)
      INTEGER iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iNumGases,iOutTypes,iNatm2
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNpmix
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      CHARACTER*7 caWord
      CHARACTER*80 caComment

      INTEGER iPrinter,iAtm,iNp,iaOp(kPathsOut)
      INTEGER iNumLinesRead,iI,iJ,iUpDown,iList,iListType
      REAL raTemp(kProfLayer)

c these are temporarily used
      INTEGER iaGasID(kMaxGas),iGS,iMPSets,iErr,iOutTypesMax
      INTEGER iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
      INTEGER iaaOpT(kMaxPrT,kPathsOut),iOut2
      REAL raaOpT(kMaxPrT,kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)

      INTEGER iaPrinter1(kMaxPrint),iaGPMPAtm1(kMaxPrint)
      INTEGER iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
      REAL raaOp1(kMaxPrint,kProfLayer)

      DO iI=1,kMaxPrint
        iaPrinter1(iI)=iaPrinter(iI)
        iaGPMPAtm1(iI)=iaGPMPAtm(iI)
        iaNP1(iI)=iaNP(iI)
        END DO
      DO iI=1,kMaxPrint
        DO iJ=1,kPathsOut
          iaaOp1(iI,iJ)=iaaOp(iI,iJ)
          END DO
        END DO
      DO iI=1,kMaxPrint
        DO iJ=1,kProfLayer
          raaOp1(iI,iJ)=raaOp(iI,iJ)
          END DO
        END DO

      caWord='*OUTPUT'
      iErr=-1

      iOutTypes=1
      iNatm2=-1

      DO iI=1,kMaxGas
        iaGasID(iI)=-100
        END DO

      iJ=0
      DO iI=1,kMaxGas
        IF (iaGases(iI) .GT. 0) THEN
          iJ=iJ+1
          iaGasID(iJ)=iI
          END IF
        END DO
      IF (iNumGases .NE. iJ) THEN
        write(kStdErr,*) 'discrepancy in number of gases stored'
        write(kStdErr,*) 'iNumGases, iGasStored = ',iNumGases,iJ
        CALL DoSTOP
        END IF

c      print *,(iaPrinter(iJ),iJ=1,kMaxPrint)
c      print *,(iaGPMPAtm(iJ),iJ=1,kMaxPrint)
c      print *,(iaNp(iJ),iJ=1,kMaxPrint)

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error reading section ',A7,' of main file')

      iNumLinesRead=1
c don't have to read in the comment!
      caComment=caComment

c figure out as a first guess how many printing options there have to be
      iOutTypesMax=0
      DO iErr=1,kMaxPrint
        IF (iaPrinter(iErr) .gt. 0) THEN
          iOutTypesMax=iOutTypesMax+1
          END IF
        END DO
      
      iOutTypes=1
      DO iOut2=1,iOutTypesMax
        iPrinter=iaPrinter1(iOut2)
        iaPrinterT(iOutTypes)=iPrinter
        IF ((iPrinter.NE.1).AND.(iPrinter.NE.2).AND.(iPrinter.NE.3))THEN
          write(kStdErr,*) 'iPrinterOption = ',iPrinter
          write(kStdErr,*)'in *OUTPUT .. Need iPrinterOption= 1,2 or 3'
          CALL DoSTOP
          END IF

c now read which atmosphere, gaspath or MP set to print and do check ...
        IF (iPrinter .EQ. 3) THEN 
          iGs=-1000
          iAtm=iaGPMPAtm1(iOut2)
          iNp=iaNp1(iOut2)
          CALL Read3(iAtm,iNp,iNpmix,iNatm,iNatm2,iListType,
     $               iaGPMPAtmT,iaNpT,iOutTypes,iaPrinterT,iPrinter)
        ELSE IF (iPrinter .EQ. 2) THEN
          iAtm=0                 !make sure this is set to 0
          iGs=iaGPMPAtm1(iOut2)
          iNp=iaNp1(iOut2)
          !***warning*** iNp can be reset from -1 to another number by Read2
          !see subroutine ReCheck_iNp
          !if initially -1, it is left as -1 if ALL mixed paths to be output
          CALL Read2(iGS,iNp,iNpmix,iMPSets,iAtm,
     $      iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,iNumGases)
        ELSE IF (iPrinter .EQ. 1) THEN 
          iAtm=0                 !make sure this is set to 0
          iGs=iaGPMPAtm1(iOut2)
          iNp=iaNp1(iOut2)
          !***warning*** iNp can be reset from -1 to another number by Read1
          !see subroutine ReCheck_iNp
          !if initially -1, it is left as -1 if ALL paths to be output
          CALL Read1(iGS,iNp,iNpmix,iMPSets,iAtm,
     $      iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,
     $      iNumGases,iaGasID)
          END IF

c if iNp=0 then nothing to be printed --- abort program
        IF (iNp .EQ. 0) THEN
          write(kStdErr,*) 'in *OUTPUT!! You have chosen a valid printing'
          write(kStdErr,*) 'option (1,2 or 3) but do not specify any '
          write(kStdErr,*) 'paths/mixed paths/radiances to be printed !!!'
          CALL DoSTOP
          END IF

c now consider the case where iNp < 0 but we only set ONE single option 
c iAtm = 0, iNp < 0 ==> we want to print out a set of mixed paths or 
c                                          a set of gas optical depths
c iAtm > 0, iNp < 0 ==> want to print out radiances at all levels of an 
c                       atmosphere
c &&&&&&&&&& BEGIN CASE 1
        IF ((iAtm .GE. 0)  .AND. (iNp .LT. 0)) THEN
          IF (iPrinter .EQ. 1) THEN
c have to output ALL the single paths, corresponding to all KProfLayers for 
c each of the iNumGases gases present
            DO iI=1,kProfLayer*iNumGases
              iaOp(iI)=iI
              iaaOpT(iOutTypes,iI)=iI
              END DO
          ELSE IF (iPrinter .EQ. 2) THEN
c have to output ALL the mixed paths, corresponding to all 1..iNpMix
            DO iI=1,iNpmix
              iaOp(iI)=iI
              iaaOpT(iOutTypes,iI)=iI
              END DO
c if iNp < 0 and iPrinter = 3, then set iaOp to "ALL" values, as follows
          ELSE IF (iPrinter .EQ. 3) THEN !set pressure fraction=1.0
            !this subroutine sets iaOp and iaaOpT
            CALL AllLayersOutputPress(iaaRadLayer,iAtm,iaNumLayer,
     $                     iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT,
     $                     raaPrBdry,raFracTop,raFracBot)
            END IF
          END IF       !IF ((iAtm .GE. 0)  .AND. (iNp .LT. 0)) THEN
c &&&&&&&&&& END CASE 1

c &&&&&&&&&& BEGIN CASE 2
c now consider the case where iNp > 0
c if iNp > 0 then look for list of layers/mixed paths/radiances
c &&&&&&&&&& BEGIN CASE 2A
        IF ((iAtm .GE. 0)  .AND. (iNp .GT. 0)) THEN
          IF (iPrinter .EQ. 3) THEN
c can easily set things for the ONE atmosphere
            DO iI=1,iNp
              raTemp(iI)=raaOp1(iOut2,iI)
              END DO
            IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN
              iUpDown=1            !upward travelling radiation
            ELSE
              iUpDown=-1           !downward travelling radiation
              END IF
            CALL DoSortReal(raTemp,iNp,-iUpDown)
            CALL PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iAtm,
     $                 iOutTypes,raTemp,iNp,raaPrBdry,raaUserPressT)
            DO iI=1,iNp
              iaaOpT(iOutTypes,iI)=iaOp(iI)
              END DO
c &&&&&&&&&& END CASE 2A

c for iPrinter = 1,2 (path or MP)
c iList     = set how many items found in list
c iListType = 1,2,3,4    1 if -1  -1
c                        2 if -1  iL
c                        3 if iGS -1
c                        4 if iGS iL

          ELSE IF (iPrinter .NE. 3) THEN
c &&&&&&&&&& BEGIN CASE 2B
            IF  ((iListType .EQ. 2) .OR. (iListType .EQ. 4)) THEN 
              IF  (iListType .EQ. 2) THEN 
                DO iI=1,iList
                  iaOp(iI)=iaaOp1(iOut2,iI)
                  END DO
                CALL DoSort(iaOp,iList)
                !having read in path/mixed path integers for ALL gas/MP sets
                !now for all gas/MP sets, set the relevant paths/mixed paths
                CALL DoSetAll(iaOp,iList,iNp,iPrinter,iNumGases,iMPSets)
              ELSE IF  (iListType .EQ. 4) THEN 
                DO iI=1,iNp
                  iaOp(iI)=iaaOp1(iOut2,iI)
                  END DO
                CALL DoSort(iaOp,iNp)
                !having read in path/mixed path integers for the gas/MP set
                !now for all gas/MP sets, set the relevant paths/mixed paths
                CALL DoSetSome(iaOp,iList,iNp,iPrinter,iNumGases,iaGasID,
     $                       iMPSets,iGS)
                END IF
c &&&&&&&&&& END CASE 2B

c &&&&&&&&&& BEGIN CASE 2C
            ELSE IF ((iPrinter .EQ. 2) .AND. (iListType .EQ. 3)) THEN  
              !have to set all MPs for a particular MP set
              DO iI=1,kProfLayer
                iaOp(iI)=iI+(iGS-1)*kProfLayer
                END DO
c &&&&&&&&&& END CASE 2C

c &&&&&&&&&& BEGIN CASE 2D
            ELSE IF ((iPrinter .EQ. 1) .AND. (iListType .EQ. 3)) THEN  !no list
              !have to set all paths for a particular gas
              !see where in the list gas iGS lies
              iI=1
 111          CONTINUE
              IF (iaGasID(iI) .NE. iGS) THEN
                iI=iI+1
                GO TO 111
                END IF
              iGS=iI
              DO iI=1,kProfLayer
                iaOp(iI)=iI+(iGS-1)*kProfLayer
                END DO
              END IF       !if iListype = 2,4 or 3
c &&&&&&&&&& END CASE 2D
            END IF         !ELSE IF (iPrinter .NE. 3) THEN 
c &&&&&&&&&& END CASE 2
          END IF           !      IF ((iNp .GT. 0) .AND. (iAtm .GE. 0)) THEN 

c whew!!! we have now parsed in Print Option (1,2 or 3), the line following
c it (iAtm/iGS  iNp) and if necessary, the next line that has the list of iNp
c paths/mixed paths/pressures ... the list is stored in iaOpT
c if iAtm > 0, everything is almost set .. all we have to do is make sure this
c option can be included in iaPrinter(kMaxPrint) etc
        IF (iAtm .GE. 0) THEN
          CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iAtm,iNpmix,
     $              iNumGases,iNp,iaNumLayer)
            !if there is a duplicate, do not update iOutTypes, else increment
            !iOutTypes by 1
            !if we find there are more than kMaxPrint options being stored,
            !then HALT!!!
          IF (iPrinter .NE. 3) THEN  !thus iAtm was set to 0 (no radiance)
            CALL CheckIaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases,
     $                   iNp,iaNpT,iaPrinterT)
            END IF 
          IF ((iPrinter .EQ. 3) .AND. (iAtm .GT. 0)) THEN
            CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases,
     $                   iNp,iaNpT,iaPrinterT,iNatm,iAtm,iaGPMPAtmT,
     $                   raaOpT,raaUserPressT,iaNumLayer,iUpDown)
            END IF
          END IF           !if iAtm .GE. 0)

c &&&&&&&&&& BEGIN CASE 3
c if iAtm < 0 and iNp > 0, we now have to set the same printing output paths 
c for all the atmospheres
        IF (iAtm .LT. 0) THEN
          IF (iPrinter .NE. 3) THEN
            write(kStdErr,*)'Something gone wrong in *OUTPUT'
            CALL DoSTOP
            END IF
c &&&&&&&&&& BEGIN CASE 3A
          IF (iNp .GT. 0) THEN 
            DO iJ=1,iNp
              raTemp(iJ)=raaOp1(iOut2,iJ)
              END DO
            DO iJ=1,iNatm
              IF (iaaRadLayer(iJ,1) .LT. iaaRadLayer(iJ,2)) THEN
                iUpDown=1            !upward travelling radiation
              ELSE
                iUpDown=-1           !downward travelling radiation
                END IF

              iaGPMPAtmT(iOutTypes)=iJ
              iaNpT(iOutTypes)=iNp 
              iaPrinterT(iOutTypes)=iPrinter
            
              CALL DoSortReal(raTemp,iNp,-iUpDown)
              CALL PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iJ,
     $         iOutTypes,raTemp,iNp,raaPrBdry,raaUserPressT)
              DO iI=1,iNp
                iaaOpT(iOutTypes,iI)=iaOp(iI)
                END DO
              CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iJ,iNpmix,
     $                 iNumGases,iNp,iaNumLayer)
              CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases,
     $                     iNp,iaNpT,iaPrinterT,iNatm,iJ,iaGPMPAtmT,
     $                     raaOpT,raaUserPressT,iaNumLayer,iUpDown)
              END DO
            END IF   !iAtm < 0, iNp > 0
c &&&&&&&&&& END CASE 3A

c &&&&&&&&&& BEGIN CASE 3B
          IF (iNp .LT. 0) THEN
            DO iJ=1,iNatm
              IF (iaaRadLayer(iJ,1) .LT. iaaRadLayer(iJ,2)) THEN
                iUpDown=1            !upward travelling radiation
              ELSE
                iUpDown=-1           !downward travelling radiation
                END IF

              iaGPMPAtmT(iOutTypes)=iJ
              iaNpT(iOutTypes)=iNp 
              iaPrinterT(iOutTypes)=iPrinter
            
              CALL AllLayersOutputPress(iaaRadLayer,iJ,iaNumLayer,
     $                     iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT,
     $                     raaPrBdry,raFracTop,raFracBot)

              CALL SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iJ,iNpmix,
     $              iNumGases,iNp,iaNumLayer)
              CALL CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases,
     $                   iNp,iaNpT,iaPrinterT,iNatm,iJ,iaGPMPAtmT,
     $                   raaOpT,raaUserPressT,iaNumLayer,iUpDown)
              END DO
            END IF   !iAtm < 0, iNp < 0
c &&&&&&&&&& END CASE 3B
          END IF     !overall iAtm < 0 if loop
c &&&&&&&&&& END CASE 3
        END DO       !do loop over iOut2

      iOutTypes=iOutTypes-1   !set it to the correct number actually found

c now set the elements of the temporary arrays/matrices to those that are 
c sent out to the main program
      CALL SetActualPrintOptions(iOutTypes,iaPrinter,iaPrinterT,
     $         iaGPMPAtm,iaGPMPAtmT,iaNp,iaNpT,iaaOp,iaaOpT,raaOp,raaOpT,
     $         raaUserPress,raaUserPressT,iaNumlayer)

      RETURN
      END

c************************************************************************
c this subroutine "explodes" the input list for ONE mixed paths set/Gas ID set
      SUBROUTINE DoSetSome(iaOp,iList,iNp,iPrinter,iNumGases,iaGasID,
     $                     iMPSets,iGS)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c iaGasID    = array containig GAS ID's
c iGS        = GAS ID or MP set number
c iList      = number of paths/mixed paths set in original iaOp
c iNp        = final number of paths/mixed paths set in iaOp for print option
c iaOp       = matrix with list of paths/MP/radiances output for each option
c iNumGases  = total number of gases read in from GASFIL + XSCFIL
c iMPSETs    = total number of mixed paths sets from *WEIGHT
c iPrinter   = printing option (1 or 2)
      INTEGER iPrinter,iNp,iaOp(kPathsOut),iMPSets,iNumGases,iList
      INTEGER iaGasID(kMaxGas),iGS
      
      INTEGER iI,iD,iaOpTemp(kPathsOut)

c make a working copy of iaOp
      DO iI=1,kPathsOut
        iaOpTemp(iI)=iaOp(iI)
        END DO

      IF (iPrinter .EQ. 2) THEN
        iD=iGS
        IF (iD .GT. iMPSets) THEN
          write (kStdErr,*) 'oops! mixed path set iGS > those set in *WEIGHT'
          CALL DoSTOP
          END IF
      ELSE
        iI=1
 11     CONTINUE
        IF (iGS .NE. iaGasID(iI)) THEN
          iI=iI+1
          GO TO 11
          END IF
        IF (iI .LE. iNumGases) THEN
          iD=iI
        ELSE
          write(kStdErr,*)'could not find gas ID iGS in list iaGASID'
          CALL DoSTOP
          END IF
        END IF

      !check that the input list read in from file has numbers <= kProfLayer
      DO iI=1,iList
        IF ((iaOp(iI) .GT. kProfLayer) .OR. (iaOp(iI) .LT. 1)) THEN
          IF (iPrinter .EQ. 1) THEN
            write(kStdErr,*) iaOp(iI)
            write(kStdErr,*)'Invalid path in list found in *OUTPUT'
          ELSEIF (iPrinter .EQ. 2) THEN
            write(kStdErr,*) iaOp(iI)
            write(kStdErr,*)'Invalid mixedpath in list found in *OUTPUT'
            END IF
          CALL DoSTOP
          END IF
        END DO
      !now that everything seems OK, go ahead and replicate list
c thus eg if iGS = 4 th gas, and the list contains 3 items 15 34 72 then we
c need the following final list 315 334 372
      DO iI=1,iList
        iaOpTemp(iI)=iaOpTemp(iI)+(iD-1)*kProfLayer
        END DO

c reset iaOp
      DO iI=1,kPathsOut
        iaOp(iI)=iaOpTemp(iI)
        END DO

      RETURN
      END
c************************************************************************
c this subroutine "explodes" the input list for ALL mixed paths sets/Gas IDs
      SUBROUTINE DoSetAll(iaOp,iList,iNp,iPrinter,iNumGases,iMPSets)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c iList      = number of paths/mixed paths set in original iaOp
c iNp        = final number of paths/mixed paths set in iaOp for print option
c iaOp       = matrix with list of paths/MP/radiances output for each option
c iNumGases  = total number of gases read in from GASFIL + XSCFIL
c iMPSETs    = total number of mixed paths sets from *WEIGHT
c iPrinter   = printing option (1 or 2)
      INTEGER iPrinter,iNp,iaOp(kPathsOut),iMPSets,iNumGases,iList
      
      INTEGER iI,iD,iaOpTemp(kPathsOut)

c make a working copy of iaOp
      DO iI=1,kPathsOut
        iaOpTemp(iI)=iaOp(iI)
        END DO

      IF (iPrinter .EQ. 1) THEN
        IF (iNp .NE. iList*iNumGases) THEN
          write(kStdErr,*)'bad *OUTPUT for paths output option -1 iN'
          write(kStdErr,*) '(iNp .NE. iList*iNumGases)'
          CALL DoSTOP
          END IF
        !check that the input list read in from file has numbers <= kProfLayer
        DO iI=1,iList
          IF ((iaOp(iI) .GT. kProfLayer) .OR. (iaOp(iI) .LT. 1)) THEN
            write(kStdErr,*) 'Invalid path in list found in *OUTPUT'
            write(kStdErr,*)'((iaOp(iI)>kProfLayer) or (iaOp(iI) <1)) '
            CALL DoSTOP
            END IF
          END DO
        !now that everything seems OK, go ahead and replicate list
c thus eg is there are 4 gases, and the list contains 3 items 15 34 72 then we
c need follwoing final list 15 34 72  115 134 172  215 234 272  315 334 372
        DO iD=1,iNumGases
          DO iI=1,iList
            iaOpTemp(iI+(iD-1)*iList)=iaOpTemp(iI)+(iD-1)*kProfLayer
            END DO
          END DO
        END IF

      IF (iPrinter .EQ. 2) THEN
        IF (iNp .NE. iList*iMPSets) THEN
          write(kStdErr,*)'bad *OUTPUT for MP output option -1 iN'
          write(kStdErr,*)'(iNp .NE. iList*iMPSets)'
          CALL DoSTOP
          END IF
        !check that the input list read in from file has numbers <= kProfLayer
        DO iI=1,iList
          IF ((iaOp(iI) .GT. kProfLayer) .OR. (iaOp(iI) .LT. 1)) THEN
            write(kStdErr,*) 'Invalid MP in list found in *OUTPUT'
            write(kStdErr,*) '(iaOp(iI) > kProfLayer)or(iaOp(iI) < 1)'
            write(kStdErr,*) 'iI,iaOp(iI) = ',iI,iaOp(iI)
            CALL DoSTOP
            END IF
          END DO
        !now that everything seems OK, go ahead and replicate list
c thus eg is there are 4 gases, and the list contains 3 items 15 34 72 then we
c need follwoing final list 15 34 72  115 134 172  215 234 272  315 334 372
        DO iD=1,iMPSets
          DO iI=1,iList
            iaOpTemp(iI+(iD-1)*iList)=iaOpTemp(iI)+(iD-1)*kProfLayer
            END DO
          END DO
        END IF

c reset iaOp
      DO iI=1,kPathsOut
        iaOp(iI)=iaOpTemp(iI)
        END DO

      RETURN
      END
c************************************************************************
c this sets the actual print option variables from the temporary ones
      SUBROUTINE SetActualPrintOptions(iOutTypes,iaPrinter,iaPrinterT,
     $         iaGPMPAtm,iaGPMPAtmT,iaNp,iaNpT,iaaOp,iaaOpT,raaOp,raaOpT,
     $         raaUserPress,raaUserPressT,iaNumLayer)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c iaPrinterT = array containing the print options (1,2 or 3)
c iaGPMPAtmT   = array containing which atmospheres to be printed (if print
c              option ii=1 or 2, then iaGPMPAtm(ii)=0)
c iaNpT      = array with  the # of paths to be printed for each print option
c iOutTypes  = total number of outputs per each 25cm-1 run
      INTEGER iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
      INTEGER iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaaOpT(kMaxPrT,kPathsOut)
      REAL raaOp(kMaxPrint,kPathsOut),raaOpT(kMaxPrT,kPathsOut)
      REAL raaUserPress(kMaxPrint,kProfLayer)
      REAL raaUserPressT(kMaxPrT,kProfLayer)
      INTEGER iOutTypes,iaNumLayer(kMaxAtm)

      INTEGER iI,iJ,iM

      DO iI=1,iOutTypes
        iaPrinter(iI)=iaPrinterT(iI)
        iaGPMPAtm(iI)=iaGPMPAtmT(iI)
        iaNp(iI)=iaNpT(iI)
        END DO

      DO iI=1,iOutTypes
        DO iJ=1,kPathsOut
          iaaOp(iI,iJ)=iaaOpT(iI,iJ)
          END DO
        END DO

      DO iI=1,iOutTypes
        DO iJ=1,kProfLayer
          raaOp(iI,iJ)=raaOpT(iI,iJ)
          raaUserPress(iI,iJ)=raaUserPressT(iI,iJ)
          END DO
        END DO

      write(kStdWarn,*)'File specified ',iOutTypes,' number of outputs '
      write(kStdWarn,*)'Output_number  OutOption(1/2/3)   Atm number'
      write(kStdWarn,*)'Num of path/MP/rads to output (if -1, list all)'
      write(kStdWarn,*)'Item      Type(1=path,2=MP,3=rad)      List'
      write(kStdWarn,*)'----------------------------------------------'
      DO ii=1,iouttypes
        write(kStdWarn,*) ii,iaPrinter(iI),iaGPMPAtm(iI)
        write(kStdWarn,*) (iaaOp(ii,ij),ij=1,iaNp(ii)) 
        IF (iaPrinter(iI) .EQ. 3) THEN
          IF (iaNp(iI) .EQ. -1) THEN
            iM=iaNumLayer(iaGPMPAtm(iI))
          ELSE
            iM=iaNp(iI)
            END IF
          write(kStdWarn,*) 'list of pressures (partial fractions)'
          DO ij=1,iM
            write(kStdWarn,*)raaUserPress(ii,ij),'(',raaOp(ii,ij),')'
            END DO
          END IF
        write(kStdWarn,*) ' '
        END DO

      RETURN
      END

c************************************************************************
c this subroutine parses in the line after finding print option3
      SUBROUTINE Read3(iAtm,iNp,iNpmix,iNatm,iNatm2,iListType,
     $               iaGPMPAtmT,iaNpT,iOutTypes,iaPrinterT,iPrinter)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c INPUT is basically iAtm,iNp which are atm # and num radiances respectively
c iaPrinterT = array containing the print options (1,2 or 3)
c iaGPMPAtmT   = array containing which atmospheres to be printed (if print
c              option ii=1 or 2, then iaGPMPAtm(ii)=0)
c iaNpT      = array with  the # of paths to be printed for each print option
c iOutTypes  = total number of outputs per each 25cm-1 run
c iNpmix     = total number of ixed paths read in from MIXFIL
c iNatm      = number of atmospheres that RADFIL had in it
c iNatm2     = number of atmospheres that OUTPUT claims there should be
c iListType  = what type printing type (all or some layers)
      INTEGER iNatm,iaPrinterT(kMaxPrT),iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT)
      INTEGER iOutTypes,iNatm2,iNpmix,iListType
      INTEGER iPrinter,iAtm,iNp

      CHARACTER*7 caWord

      caWord='*OUTPUT'

      IF (iNpmix .LT. 0) THEN 
        write(kStdErr,*)'Have not found *WEIGHTS ... cannot process' 
        write(kStdErr,*)'make sure *WEIGHT is set before *OUTPUT' 
        CALL DoSTOP 
        END IF 
      IF (iNatm .EQ. 0) THEN 
        write(kStdErr,*)'Have not found *RADNCE ... cannot process' 
        write(kStdErr,*)'make sure *RADNCE is before *OUTPUT' 
        CALL DoSTOP 
        END IF 
 
      IF (iAtm .EQ. 0) THEN 
        write(kStdErr,*) 'Cannot output spectra for atmosphere # 0 !!' 
        write(kStdErr,*)'check *OUTPUT!!' 
        CALL DoSTOP 
        END IF 
      IF (iAtm .GT. iNatm) THEN 
        write(kStdErr,*) 'atmosphere # ',iAtm,' undefined in *RADNCE' 
        write(kStdErr,*)'check *OUTPUT!!' 
        CALL DoSTOP 
        END IF 

      IF (iAtm .GT. iNatm2) THEN 
        iNatm2=iAtm 
        END IF 
      IF (iNp .LT. 0) THEN 
        iListType=1 
      ELSE 
        iListType=2 
        END IF 
 
c do things one at a time!!!! ie even if iAtm < 0, think of ONE atmosphere!
      IF (iAtm .GT. 0) THEN 
        iaGPMPAtmT(iOutTypes)=iAtm 
      ELSE 
        iaGPMPAtmT(iOutTypes)=1 
        END IF  

      iaNpT(iOutTypes)=iNp  
      iaPrinterT(iOutTypes)=iPrinter 
 
      RETURN
      END

c************************************************************************
c this subroutine parses in the line after finding print option2
      SUBROUTINE Read2(iGS,iNp,iNpmix,iMPSets,iAtm,
     $      iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,iNumGases)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c input is basically iGas and iNp which are MPset# and num MPs to be output
c iMPSets    = number of mixed path sets found in *WEIGHT
c iGS        = which GAS ID/mixed path set
c iaGPMPAtmT   = array containing which atmospheres to be printed (if print
c              option ii=1 or 2, then iaGPMPAtm(ii)=0)
c iaNpT      = array with  the # of paths to be printed for each print option
c iOutTypes  = total number of outputs per each 25cm-1 run
c iNpmix     = total number of ixed paths read in from MIXFIL
c iNatm      = number of atmospheres that RADFIL had in it
c iNatm2     = number of atmospheres that OUTPUT claims there should be
c iList      = how many paths/mixed paths to be output (-1 or > 0)
c iListType  = what type printing type (all or some layers)
c iAtm       = dummy for this (set to 0)
      INTEGER iAtm,iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT),iGS,iMPSets
      INTEGER iOutTypes,iNpmix,iListType,iList,iPrinter,iNp,iNumGases

      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iFloor

      caWord='*OUTPUT'

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error reading section ',A7,' of main file')

      iAtm=0 
      IF (iNpmix .LT. 0) THEN
        write(kStdErr,*)'Have not found mixing table ... cannot process'
        write(kStdErr,*)'make sure *WEIGHT is set before *OUTPUT'
        CALL DoSTOP
        END IF
      iMPSets=iFloor(iNpmix*1.0/(kProfLayer))

      IF (iGS .EQ. 0) THEN
        write(kStdErr,*)'Cannot output spectra for mixed path set # 0 !'
        write(kStdErr,*)'check *OUTPUT!!'
        CALL DoSTOP
        END IF
      IF (iGS .GT. iMPSets) THEN
        write(kStdErr,*)'MP set # ',iGS,' undefined in *WEIGHT'
        write(kStdErr,*)'check *OUTPUT!!'
        CALL DoSTOP
        END IF

c now recheck iNp and reset if necessary
      CALL ReCheck_iNp(iList,iListType,iGS,iNp,iMPSets,iNumGases,2)
      iaGPMPAtmT(iOutTypes)=iAtm     !=0 as NO radiances to be output 
      iaNpT(iOutTypes)=iNp 

      RETURN
      END

c************************************************************************
c this subroutine parses in the line after finding print option1
      SUBROUTINE Read1(iGS,iNp,iNpmix,iMPSets,iAtm,
     $      iList,iListType,iaGPMPAtmT,iaNpT,iOutTypes,iPrinter,
     $      iNumGases,iaGasID)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c input is basically iGas and iNp which are Gas# and numPaths to be output
c iMPSets    = number of mixed path sets found in *WEIGHT
c iGS        = which GAS ID/mixed path set
c iaGPMPAtmT   = array containing which atmospheres to be printed (if print
c              option ii=1 or 2, then iaGPMPAtm(ii)=0)
c iaNpT      = array with  the # of paths to be printed for each print option
c iOutTypes  = total number of outputs per each 25cm-1 run
c iNpmix     = total number of ixed paths read in from MIXFIL
c iNatm      = number of atmospheres that RADFIL had in it
c iNatm2     = number of atmospheres that OUTPUT claims there should be
c iList      = how many paths/mixed paths to be output (-1 or > 0)
c iListType  = what type printing type (all or some layers)
c iAtm       = dummy for this (set to 0)
      INTEGER iAtm,iaGPMPAtmT(kMaxPrT),iaNpT(kMaxPrT),iGS,iMPSets
      INTEGER iOutTypes,iNpmix,iListType,iList,iPrinter,iNp
      INTEGER iNumgases,iaGasID(kMaxGas)

      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iI,DoOutputLayer

      caWord='*OUTPUT'

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error reading section ',A7,' of main file')

      iNumLinesRead=1

      iAtm=0

      IF (iGS .EQ. 0) THEN
        write(kStdErr,*)'Cannot output spectra for GAS ID  = 0 !!'
        write(kStdErr,*)'check *OUTPUT!!'
        CALL DoSTOP
        END IF
      IF (iGS .GT. kMaxGas) THEN
        write(kStdErr,*)'Gas ID ',iGS,' undefined '
        write(kStdErr,*)'check *OUTPUT!!'
        CALL DoSTOP
        END IF
      IF ((DoOutputLayer(iGS,iNumGases,iaGasID) .LT. 0) .AND.
     $        (iGS .GT. 0)) THEN
        write(kStdErr,*)'List of Gas IDs from MOLGAS/XSCGAS = '
        write(kStdErr,*)(iaGasID(iI),iI=1,iNumGases)
        write(kStdErr,*)'Gas ID ',iGS,' not set from *MOLGAS,*XSCGAS '
        write(kStdErr,*)'check *OUTPUT!!'
        CALL DoSTOP
        END IF

c now recheck iNp and reset if necessary
      CALL ReCheck_iNp(iList,iListType,iGS,iNp,iMPSets,iNumGases,1)
      iaGPMPAtmT(iOutTypes)=iAtm     !=0 as NO radiances to be output 
      iaNpT(iOutTypes)=iNp 

      RETURN
      END

c************************************************************************
c this subroutine does the checking of iNp .. if necessary, it resets it
      SUBROUTINE ReCheck_iNp(iList,iListType,iGS,iNp,
     $                       iMPSets,iNumGases,iPrinter)

      include 'kcarta.param'

c iList     = set how many items found in list
c iListType = 1,2,3,4    1 if -1  -1
c                        2 if -1  iL
c                        3 if iGS -1
c                        4 if iGS iL
c iGS       = GAS ID os Mixed path set
c iNp       = from inout file, how many gas paths/mixed paths to be output

      INTEGER iList,iListType,iGS,iNp,iMPSets,iNumGases,iPrinter

      INTEGER iI

      IF (iPrinter .EQ. 1) THEN
        iI=iNumGases
      ELSE        
        iI=iMPSets
        END IF

      iList=0                !set how many items to be found in list
      IF ((iGS .EQ. -1) .AND. (iNp .EQ. -1)) THEN
        iList=iNp
        iListType=1
        iNp = -1             !for ALL gases/MP sets, output all layers
      ELSE IF ((iGS .EQ. -1) .AND. (iNp .NE. -1)) THEN
        iList=iNp
        iListType=2
        iNp = iNp*iI         !for ALL gases/MP sets, output some layers
      ELSE IF ((iGS .GT. 0) .AND. (iNp .EQ. -1)) THEN
        iList=iNp
        iListType=3
        iNp = kProfLayer      !for one gas/MP set, output all layers
      ELSE IF ((iGS .GT. 0) .AND. (iNp .NE. -1)) THEN
        iList=iNp
        iListType=4
        iNp = iNp            !for one gas/MP set, output iNp layers
        END IF

      RETURN
      END

c************************************************************************
c this subroutine takes care of the case iPrinter = 3 (radiance), iAtm > 0,
c (and so the user has a specific atmosphere in mind), iNp = -1 
      SUBROUTINE AllLayersOutputPress(iaaRadLayer,iAtm,iaNumLayer,
     $                     iaOp,iaaOpT,raaOpT,iOutTypes,raaUserPressT,
     $                     raaPrBdry,raFracTop,raFracBot)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c iaaRadLayer= list of layers in each radiating atmosphere
c iaaOpT     = matrix with list of paths/MP/radiances output for each option
c raaOpT     = matrix with fractional list of radiances output
c raaUserPressT= matrix with user specified pressures
c iaNumLayer = number of layers in each atmosphere, from *RADFIL
c iOutTypes  = total number of outputs per each 25cm-1 run
c raFracTop/Bot = how the top/bottom layers of mixing table have been affected
c                 by the mixing table
c raaPrBdry     = pressure boundaries of the atmospheres
c iAtm          = which atmoshere we are considering
c iaOp          = temporary part of iaaOpT
      INTEGER iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iAtm,iaOp(kPathsOut)
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raaOpT(kMaxPrT,kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)

      INTEGER iUpDown,iI,iLay,MP2Lay

      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN
        iUpDown = 1        !radiation going up
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,2)) THEN
        iUpDown = -1       !radiation going down
        END IF
c the following loop just does the number of layers in the relevant atmosphere
c taking care of the case iAtm > 0,iNp < 0
      DO iI=2,iaNumLayer(iAtm)-1    !do complete layers
        iaOp(iI)=iI
        iaaOpT(iOutTypes,iI)=iI
        raaOpT(iOutTypes,iI)=1.0  !define complete layer (fraction=1.0)
        iLay=MP2Lay(iaaRadLayer(iAtm,iI))
        IF (iUpDown .GT. 0) THEN
          raaUserPressT(iOutTypes,iI)=plev(iLay+1)
        ELSE
          raaUserPressT(iOutTypes,iI)=plev(iLay)
          END IF
        END DO
        IF (iUpDown .GT. 0) THEN   
          !special : bottom 
          iI=1
          iaOp(iI)=iI
          iaaOpT(iOutTypes,iI)=iI
          iLay=MP2Lay(iaaRadLayer(iAtm,iI))
          raaUserPressT(iOutTypes,iI)=plev(iLay+1)-delta
          !since we flip defns of partial layers for bottom ==> fraction = 0.0
          raaOpT(iOutTypes,iI)=0.0 !was incorrectly raFracBot(iAtm)

          !special : top
          iI=iaNumLayer(iAtm)
          iaOp(iI)=iI
          iaaOpT(iOutTypes,iI)=iI
          raaUserPressT(iOutTypes,iI)=raaPrBdry(iAtm,2)
          raaOpT(iOutTypes,iI)=raFracTop(iAtm)

        ELSE IF (iUpDown .LT. 0) THEN
          !special : top
          iI=1
          iaOp(iI)=iI
          iaaOpT(iOutTypes,iI)=iI
          iLay=MP2Lay(iaaRadLayer(iAtm,iI))
          raaUserPressT(iOutTypes,iI)=plev(iLay)
          !since we flip defns of partial layers for top ==> fraction = 0.0
          raaOpT(iOutTypes,iI)=0.0 !was incorrectly raFracTop(iAtm)

          !special : bottom 
          iI=iaNumLayer(iAtm)
          iaOp(iI)=iI
          iaaOpT(iOutTypes,iI)=iI
          raaUserPressT(iOutTypes,iI)=raaPrBdry(iAtm,2)
          raaOpT(iOutTypes,iI)=raFracBot(iAtm)
         END IF

      RETURN
      END
c************************************************************************
c this subroutine checks to see iaOp is kosher
      SUBROUTINE SetIaaOpT(iaaOpT,iaOp,iOutTypes,iPrinter,iAtm,iNpmix,
     $                    iNumGases,iNp,iaNumLayer)
 
      include 'kcarta.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c iPrinter   = which is the current printing option
c iNp        = curreent number of paths/mixed paths/radiances to print
c iNpmix     = number of mixed paths
c iNumGases  = number of gases
c iaaOpT     = matrix with list of paths/MP/radiances output for each option
c iaNumLayer = number of layers in each atmosphere, from *RADFIL
c iOutTypes  = total number of outputs per each 25cm-1 run
c iAtm          = which atmoshere we are considering
c iaOp          = temporary part of iaaOpT
      INTEGER iaaOpT(kMaxPrT,kPathsOut),iaOp(kPathsOut),iOutTypes
      INTEGER iPrinter,iAtm,iNpmix,iNumGases,iNp,iaNumLayer(kMaxAtm)

      INTEGER iI

      DO iI=1,iNp
        iaaOpT(iOutTypes,iI)=iaOp(iI)
        IF ((iPrinter .EQ. 3) .AND. (iAtm .GT. 0)) THEN
c all numbers here should be between 1 and iaNumLayer(relevant atm)
          IF ((iaOp(iI) .LT. 1).OR.(iaOp(iI) .GT. iaNumlayer(iAtm))) 
     $      THEN  
            write(kStdErr,*)'Printing option selected = ',iPrinter
            write(kStdErr,*)'layer # ',iI,' = ',iaOp(iI)
            write(kStdErr,*)'Error in reading *output section'
            write(kStdErr,*)'Radiating MP must be between 1 and ',
     $                 iaNumLayer(iAtm)
            CALL DoSTOP
            END IF
          END IF
        IF (iPrinter .EQ. 2) THEN
c all numbers here should be between 1 and iNpmix
          IF ((iaOp(iI) .LT. 1).OR.(iaOp(iI).GT.iNpmix)) THEN  
            write(kStdErr,*)'Printing option selected = ',iPrinter
            write(kStdErr,*)'mixed path # ',iI,' = ',iaOp(iI)
            write(kStdErr,*)'Error in reading *output section'
            write(kStdErr,*)'Mixed Path number must be between 1 and ',
     $                 iNpmix
            CALL DoSTOP
            END IF
          END IF
        IF ((iPrinter .EQ. 1)) THEN
c all numbers here should be between 1 and kProfLayer*iNumGases
          IF ((iaOp(iI) .LT. 1).OR.(iaOp(iI).GT.kProfLayer*iNumGases))
     $      THEN  
            write(kStdErr,*)'Printing option selected = ',iPrinter
            write(kStdErr,*)'# of gases in GASFIL+XSCFIL = ',iNumGases
            write(kStdErr,*)'path # ',iI,' = ',iaOp(iI)
            write(kStdErr,*)'Error in reading *output section'
            write(kStdErr,*)'Path number must be between 1 and ',
     $                 kProfLayer*iNumGases
            CALL DoSTOP
            END IF
          END IF
        END DO

      RETURN
      END
c************************************************************************
c this subroutine checks to see if iPrinter=3 has already been stored
c if so then it supplements the previous relevant vector of iaaOpT,
c raaOpT,raaUserPressT with the new information
c this is only to check ONE atmosphere being added on
      SUBROUTINE CheckRaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,iNumGases,
     $                   iNp,iaNpT,iaPrinterT,iNatm,iAtm,iaGPMPAtmT,
     $                   raaOpT,raaUserPressT,iaNumLayer,iUD)

      include 'kcarta.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c iUD        = is radiation travelling up or down
c iaNumLayer = number of radiating layers in the defined atmospheres
c iNatm      = number of atmospheres read in from *RADNCE
c iAtm       = which atmosphere(s) this print option is concerned with
c iaGPMPAtmT   = array containig atmopsphere printing options
c raaOpT     = matrix list of fractions to output radiances at
c raaOpT     = matrix list of pressures to output radiances at
c iPrinter   = which is the current printing option
c iNp        = curreent number of paths/mixed paths/radiances to print
c iaNpT      = array containig number of paths/mixed paths/radiances to print
c iaPrinterT = array containig printing options
c iNpmix     = number of mixed paths
c iNumGases  = number of gases
c iaaOpT     = matrix with list of paths/MP/radiances output for each option
c iOutTypes  = total number of outputs per each 25cm-1 run
      INTEGER iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iaPrinterT(kMaxPrT)
      INTEGER iPrinter,iNpmix,iNumGases,iNp,iaNpT(kMaxPrT),iUD
      INTEGER iAtm,iNatm,iaGPMPAtmT(kMaxPrT),iaNumLayer(kMaxAtm)
      REAL raaOpT(kMaxPrT,kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)
    
      INTEGER iI,iJ,iJ1,iFound,iaOp1(kPathsOut),iaOp2(kPathsOut),i12
      REAL ra1(kProfLayer),ra2(kProfLayer),raP1(kProfLayer),raP2(kProfLayer)

      !check to see if any of the previous print options  == iPrinter,iAtm
      iFound = -1
      DO iI=1,iOutTypes-1
        IF ((iaPrinterT(iI).EQ.iPrinter).AND.(iaGPMPAtmT(iI).EQ.iAtm))THEN
          iFound = iI
          END IF
        END DO

      IF (iFound .GT. 0) THEN    !have found a match!!! have to merge together
        !first check to see if one or the other option wants ALL paths output
        IF (iaNpT(iFound) .EQ. -1) THEN
          iaNpT(iFound) = -1      !old option included ALL layers
          iJ1=iFound
          END IF
        IF (iaNpT(iOutTypes) .EQ. -1) THEN
          iaNpT(iFound) = -1      !new option includes ALL layers
          iJ1=iOutTypes
          END IF                  
        IF (iaNpT(iFound) .EQ. -1) THEN
          iJ=iaNumLayer(iAtm)
          DO iI=1,iJ
            iaaOpT(iFound,iI)=iI
            raaOpT(iFound,iI)=raaOpT(iJ1,iI)
            raaUserPressT(iFound,iI)=raaUserPressT(iJ1,iI)
            END DO
          END IF

        IF ((iaNpT(iOutTypes) .NE. -1).AND.(iaNpT(iFound) .NE. -1)) THEN
          !have to slug it out! ah well!!!!
          DO iI=1,iaNpT(iFound)
            iaOp1(iI)=iaaOpT(iFound,iI)
            ra1(iI)=raaOpT(iFound,iI)
            raP1(iI)=raaUserPressT(iFound,iI)
            END DO
          DO iI=1,iaNpT(iOutTypes)
            iaOp2(iI)=iaaOpT(iOutTypes,iI)
            ra2(iI)=raaOpT(iOutTypes,iI)
            raP2(iI)=raaUserPressT(iOutTypes,iI)
            END DO
          CALL MergeRads(iaOp1,iaOp2,ra1,ra2,raP1,raP2,
     $                   iaNpT(iFound),iaNpT(iOutTypes),i12,iUD)
          !now update raaOpT,raaUserPressT,iaaOpT,iaNpT
          iaNpT(iFound)=i12
          DO iI=1,i12
            iaaOpT(iFound,iI)=iaOp1(iI)
            raaOpT(iFound,iI)=ra1(iI)
            raaUserPressT(iFound,iI)=raP1(iI)
            END DO
          END IF
        END IF

      IF (iFound .LT. 0) THEN !we really have found a new print option
        iOutTypes=iOutTypes+1
        END IF

      IF (iOutTypes .GT. kMaxPrT) THEN
        write(kStdErr,*)'Cannot store so many printing options!!!' 
        write(kStdErr,*)'Either increase kMaxPrint in kcarta.param, or '
        write(kStdErr,*)'decrease # of print options found in *OUTPUT'
        CALL DoSTOP
        END IF

      RETURN
      END

c************************************************************************
c this subroutine checks to see if iPrinter=1 or 2 has already been stored
c if so then it supplements the previous relevant vector of iaaOpT with the
c new information
      SUBROUTINE CheckIaaOpT(iaaOpT,iOutTypes,iPrinter,iNpmix,
     $                      iNumGases,iNp,iaNpT,iaPrinterT)

      include 'kcarta.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c iPrinter   = which is the current printing option
c iNp        = curreent number of paths/mixed paths/radiances to print
c iaNpT      = array containig number of paths/mixed paths/radiances to print
c iaPrinterT = array containig printing options
c iNpmix     = number of mixed paths
c iNumGases  = number of gases
c iaaOpT      = matrix with list of paths/MP/radiances output for each option
c iOutTypes  = total number of outputs per each 25cm-1 run
      INTEGER iaaOpT(kMaxPrT,kPathsOut),iOutTypes,iaPrinterT(kMaxPrT)
      INTEGER iPrinter,iNpmix,iNumGases,iNp,iaNpT(kMaxPrT)

      INTEGER iI,iJ,iFound,iaOp1(kPathsOut),iaOp2(kPathsOut)

      !check to see if any of the previous print options  == iPrinter
      iFound = -1
      DO iI=1,iOutTypes-1
        IF (iaPrinterT(iI) .EQ. iPrinter) THEN
          iFound = iI
          END IF
        END DO

      IF (iFound .GT. 0) THEN    !have found a match!!! have to merge together
        !first check to see if one or the other option wants ALL paths output
        IF (iaNpT(iFound) .EQ. -1) THEN
          iaNpT(iFound) = -1      !old option included ALL paths/mixed paths
          END IF
        IF (iaNpT(iOutTypes) .EQ. -1) THEN
          iaNpT(iFound) = -1      !new option includes ALL paths/mixed paths
          END IF                  
        IF (iaNpT(iFound) .EQ. -1) THEN
          IF (iPrinter .EQ. 1) THEN
            iJ=iNumGases*kProfLayer
          ELSE
            iJ=iNpMix
            END IF
          DO iI=1,iJ
            iaaOpT(iFound,iI)=iI
            END DO
          END IF

        IF ((iaNpT(iOutTypes) .NE. -1).AND.(iaNpT(iFound) .NE. -1)) THEN
          !have to slug it out! ah well!!!!
          DO iI=1,iaNpT(iFound)
            iaOp1(iI)=iaaOpT(iFound,iI)
            END DO
          DO iI=1,iaNpT(iOutTypes)
            iaOp2(iI)=iaaOpT(iOutTypes,iI)
            END DO
          CALL MergeLists(iaOp1,iaOp2,iaNpT(iFound),iaNpT(iOutTypes),iI)
          !now update iaaOpT,iaNpT
          iaNpT(iFound)=iI
          DO iI=1,iaNpT(iFound)
            iaaOpT(iFound,iI)=iaOp1(iI)
            END DO
          END IF
        END IF

      IF (iFound .LT. 0) THEN !we really have found a new print option
        iOutTypes=iOutTypes+1
        END IF

      IF (iOutTypes .GT. kMaxPrT) THEN
        write(kStdErr,*)'Cannot store so many printing options!!!' 
        write(kStdErr,*)'Either increase kMaxPrint in kcarta.param, or'
        write(kStdErr,*)'decrease # of print options found in *OUTPUT'
        CALL DoSTOP
        END IF

      RETURN
      END
c************************************************************************
c this subroutine takes two lists, merges them together. If duplicate items 
c are found in the second list, they are discarded. The new list is then 
c sorted and stored in list1. The total length of the new list1 is iSum
      SUBROUTINE MergeLists(iaOp1,iaOp2,i1,i2,iSum)

      include 'kcarta.param'

      INTEGER iaOp1(kPathsOut),iaOp2(kPathsOut),i1,i2,iSum

      INTEGER iI,iaRes(kPathsOut),DoOutputLayer,iOK

      !first put in elements of iaOp1 into iaRes
      DO iI=1,i1
        iaRes(iI)=iaOp1(iI)
        END DO
      iSum=i1

      !now put in elements of iaOp2 into iaRes, ensuring no duplicates
      DO iI=1,i2
        iOK=DoOutputLayer(iaOp2(iI),iSum,iaRes)
        IF (iOK .LT. 0) THEN   !iaOp2(iI) not found .. add to list
          iSum=iSum+1
          iaRes(iSum)=iaOp2(iI)
          CALL dosort(iaRes,iSum)
        ELSE
          write(kStdWarn,*) 'following is duplicate : ',iaOp2(iI)
          END IF
        END DO

      !sort iaRes
      CALL dosort(iaRes,iSum)

      !update iaOp1      
      DO iI=1,iSum
        iaOp1(iI)=iaRes(iI)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine takes 6 lists, merges them together. If duplicate items 
c are found in the second list, they are discarded. The new list is then 
c sorted and stored in list1. The total length of the new list1 is iSum
c the merging is done based elememts in pressure list raP1,raP2
      SUBROUTINE MergeRads(ia1,ia2,ra1,ra2,raP1,raP2,i1,i2,iSum,iUD)

      include 'kcarta.param'

c i1,i2     = lenght of the lits blah1,blah2
c iSum      = final length of combined list .. stored in ia1,ra1,raP1
c iUD       = direction of radiation travel
c ra1,ra2   = arrays containing layer fractions
c raP1,raP2 = arrays containing pressure levels
c ia1,ia2   = arrays containing layers
      INTEGER ia1(kPathsOut),ia2(kPathsOut),i1,i2,iSum,iUD
      REAL raP1(kProfLayer),raP2(kProfLayer),ra1(kProfLayer),ra2(kProfLayer)

      INTEGER iI,iaRes(kPathsOut),iNew
      REAL raRes(kProfLayer),raPRes(kProfLayer)

      !first put in elements of ia1 into iaRes
      DO iI=1,i1
        iaRes(iI)=ia1(iI)
        raRes(iI)=ra1(iI)
        raPRes(iI)=raP1(iI)
        END DO
      iSum=i1

      !now put in elements of ia2 into iaRes, ensuring no duplicates
      DO iI=1,i2
        CALL DoMergePress(raP2(iI),raPRes,ia2(iI),iaRes,ra2(iI),raRes,
     $           iSum,iUD,iNew)
        IF (iNew .LT. 0) THEN   !raP2(iI) already present in array .. 
          write(kStdWarn,*) 'following is duplicate : ',raP2(iI)
          END IF
        END DO

      CALL DoSortPress(iaRes,raRes,raPRes,iSum,-iUD)

      !update ia1      
      DO iI=1,iSum
        ia1(iI)=iaRes(iI)
        ra1(iI)=raRes(iI)
        raP1(iI)=raPRes(iI)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine checks to see if we insert pressure rP into array raPRes
c if so then it also inserts some other elements into the respective arrays
      SUBROUTINE DoMergePress(rP,raPRes,i2,iaRes,r2,raRes,iSum,
     $                              iUD,iNew)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c iNew   = is this a new pressure?If so iNew = +1
c iUD    = direction of radiation travel
c rP     = see if this is a new pressure
c raPRes = this is the array we check to see if rP exists inside
c i2     = this is the pressure layer assocoated with rP
c iaRes  = this is array containing the pressure layer numbers to be output
c r2     = this is the fraction associated with rP
c raRes  = this is the array containing the fractions of the press. layers
c          to be output
c iSum   = how far down the arrays to search
      REAL rP,r2
      INTEGER iaRes(kPathsOut),i2,iNew,iSum,iUD
      REAL raRes(kProfLayer),raPRes(kProfLayer)

      INTEGER iI

c recall the pressures are stored in in(de)creasing order, depending on iUD
c but also note that as subroutine DoOutptuRadianec checks to see if the LAYER
c being processed has an output pressure level, and it goes thru the WHOLE
c list while doing the check, we do not REALLY have to set everything in
c increasing or decreasing order (ie iUD not really needed)
      iNew=1   !assume new pressure
      iI=1
 10   CONTINUE
      IF (abs(raPRes(iI)-rP) .LE. delta) THEN !this pressure was there before
        iNew=-1
        END IF
      IF ((iI .LT. iSum) .AND. (iNew .GT. 0)) THEN 
        !still not found pressure and have some more elements to check in array
        iI=iI+1
        GO TO 10
        END IF

      IF (iNew .GT. 0) THEN     !pressure not found .. add it on
        IF (iSum .LT. kProfLayer) THEN !we add on relevant elements to arrays
          iSum=iSum+1
          raPRes(iSum)=rP
          iaRes(iSum)=i2
          raRes(iSum)=r2
        ELSE IF (iSum .EQ. kProfLayer) THEN !we cannot add it on
          write(kStdErr,*)'cannot merge press',rP,'into previous list'
          CALL DoSTOP
          END IF
        END IF

      RETURN
      END

c************************************************************************
c this subroutine converts the pressures found in *OUTPUT to which 
c layer this is in the radiating atmosphere
      SUBROUTINE PressTOLayers(raaOpT,iaOp,iaNumLayer,iaaRadLayer,iAtm,
     $         iOutTypes,raPress,iNp,raaPrBdry,raaUserPressT)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      INTEGER kMaxPrT
      PARAMETER (kMaxPrT=kMaxPrint+1)

c raaOpT     = matrix with fractional list of layers where radiances output
c iaNumLayer = number of layers in each atmosphere, from *RADFIL
c iaaRadLayer= list of layers in each atmosphere, from *RADFIL
c raaPrBdry  = matrix containing start/stop pressures for each atmosphere
c iAtm       = which atm
c iOutTypes  = which output option this is
c raaUserPressT= which pressures to output the sutff at

c raPress     = array containing pressure levels where radiances to be output
c iNp        = number of pressures where radiances are to be output

c iaOp       = this subroutine figures out which layers are to be output

      INTEGER iaOp(kPathsOut),iAtm,iNp,iOutTypes
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raaOpT(kMaxPrT,kProfLayer),raaPrBdry(kMaxAtm,2)
      REAL raPress(kProfLayer),raaUserPressT(kMaxPrT,kProfLayer)

      REAL rP,rF
      INTEGER iI,i1,i2,iOne,iFound,iDirection,iFloor,iW,iOffSet
      INTEGER iStartLay,MP2Lay

c check the radiation direction
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN
        iDirection=1     !going up
      ELSE
        iDirection=-1    !going down
        END IF

c find which set of 100 mixed paths the layers for this atmosphere are from
c 1-100 = set 1, 101-200 = set 2, 201-300 = set 3 etc
c note : the reason i use iaaRadLayer(iAtm,2) is as follows : suppose we have
c uplooking instr, defined between mixed paths 100 --> 1  ... if we used
c         iW == iFloor((iaaRadLayer(iAtm,1)*1.0)/(kProfLayer*1.0))+1 
c            == iFloor(100.0/100.0)+1=2  instead of 1

      iW=iFloor((iaaRadLayer(iAtm,2)*1.0)/(kProfLayer*1.0))+1
      iOffSet=(iW-1)*kProfLayer
      
c set which radiating pressure layer atmosphere starts at, to compare to iI
      iStartLay=MP2Lay(iaaRadLayer(iAtm,1))

      DO iI=1,iNp

c see if this pressure lies within the start/stop defined in raaPrBdry
c if not, then redefine
        rP=raPress(iI)
        IF (iDirection .GT. 0) THEN
          IF (rP .GT. raaPrBdry(iAtm,1)) THEN
            rP = raaPrBdry(iAtm,1)-delta*1.001
            rP = raaPrBdry(iAtm,1)
            write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
            write(kStdWarn,*) raPress(iI),' to ',rP
            END IF
          IF (rP .LT. raaPrBdry(iAtm,2)) THEN
            rP = raaPrBdry(iAtm,2)+delta*1.001
            rP = raaPrBdry(iAtm,2)
            write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
            write(kStdWarn,*) raPress(iI),' to ',rP
            END IF
          END IF
        IF (iDirection .LT. 0) THEN
          IF (rP .LT. raaPrBdry(iAtm,1)) THEN
            rP = raaPrBdry(iAtm,1)+delta*1.001
            rP = raaPrBdry(iAtm,1)
            write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
            write(kStdWarn,*) raPress(iI),' to ',rP
            END IF
          IF (rP .GT. raaPrBdry(iAtm,2)) THEN
            rP = raaPrBdry(iAtm,2)-delta*1.001
            rP = raaPrBdry(iAtm,2)
            write(kStdWarn,*) 'For atm ',iAtm,' reset output press from'
            write(kStdWarn,*) raPress(iI),' to ',rP
            END IF
          END IF

c now find out the pressure layer this pressure lies within
        iFound=-1
        i1=1
        i2=2
 10     CONTINUE
        IF ((rP .LE.plev(i1)).AND.(rP .GT.plev(i2))) THEN
          iFound=1
          iOne=i1+(iW-1)*kProfLayer  !set this correctly as Mixed Path number
          END IF
        IF ((iFound .LT. 0) .AND. (i1 .LT. kProfLayer)) THEN
          i1=i1+1
          i2=i2+1
          GO TO 10
          END IF
        IF ((iFound .LT. 0)) THEN
          IF (abs(rP-plev(kProfLayer+1)) .LE. delta) THEN
            i1=kProfLayer
            iFound=1
            iOne=i1+(iW-1)*kProfLayer  !set this correctly as Mix Path number
          ELSE
            write(kStdErr,*) 'atm#',iAtm,' could not find pressure ',rP
            write(kStdErr,*) 'within KLAYERS pressure levels. Check'
            write(kStdErr,*) '*RADNCE and *OUTPUT sections'
            CALL DoSTOP
            END IF
          END IF

c now see which radiating layer this corresponds to in defined atmosphere
c eg if defined atmospghere was 3 950.0  0.15
c then 950.0 is within 6th pressure layer, 0.15 within 96th pressure layer
c and the 3 implies this is the third group of 100 layer weigths 
c ==> this atmosphere consists of the 91 mixed paths 
c       206,207,208, ... 294,295,296
c with mixed paths 206 and 296 having fractional weights
c
c so now if the pressure at which we want radiance to be output is 151.2664,
c this is the 50th AIRS pressure level, but it would be the (250-206+1)=45th
c radiating layer in the defined atmosphere. Ditto upward looking instrument

c we know rP has been checked to lie within the pressure boundaries of 
c the current atm, so if this is bottom or topmost layer, cannot exceed 
c rFracBot or rFracTop respectively
        IF (iDirection .GT. 0) THEN   !radiation goes up
          iaOp(iI)=iOne-iaaRadLayer(iAtm,1)+1
          !want BOTTOM fraction of layer
          rF=(plev(i1)-rP)/(plev(i1)-plev(i2))
          IF (iOne .EQ. (iaaRadLayer(iAtm,1))) THEN  
            !want TOP fraction of  lowest layer
c            print *,'aha : bottom layer for down look instr!!!!!!'
            rF=1-rF
            END IF
          IF (abs(rF-1.00) .LE. delta) THEN
            rF=1.00
            END IF
          raaUserPressT(iOutTypes,iI)=rP
          raaOpT(iOutTypes,iI)=rF
          write(kStdWarn,*) iDirection,rP,rF

        ELSE IF (iDirection .LT. 0) THEN !radiation goes up
          iaOp(iI)=iaaRadLayer(iAtm,1)-iOne+1
          !want BOTTOM fraction of layer
          rF=(rP-plev(i2))/(plev(i1)-plev(i2))
          IF (iOne .EQ. (iaaRadLayer(iAtm,1))) THEN  
            !want BOTTOM  fraction of highest layer
c            print *,'aha : top layer for up look instr!!!!!!'
            rF=1-rF
            END IF
          IF (abs(rF-1.00) .LE. delta) THEN
            rF=1.00
            END IF
          raaUserPressT(iOutTypes,iI)=rP
          raaOpT(iOutTypes,iI)=rF
          write(kStdWarn,*) iDirection,rP,rF
          END IF

        IF ((iaOp(iI) .GT. iaNumLayer(iAtm)).OR.(iaOp(iI) .LT. 1))THEN
          write(kStdErr,*) 'iaaRadLayer(iAtm,1),iOne,=',
     $iaaRadLayer(iAtm,1),iOne
          write(kStdErr,*) 'iI,iaOp(iI) = ',iI,iaOp(iI)
          write(kStdErr,*) 'Pressure ',rP,'invalid for atm# ',iAtm
          CALL DoSTOP
          END IF

        IF (raaOpT(iOutTypes,iI) .LT. 0.0) THEN
          raaOpT(iOutTypes,iI) = 0.0
          END IF
        IF (raaOpT(iOutTypes,iI) .GT. 1.0) THEN
          raaOpT(iOutTypes,iI) = 1.0
          END IF
 
        END DO

      RETURN
      END

c************************************************************************
