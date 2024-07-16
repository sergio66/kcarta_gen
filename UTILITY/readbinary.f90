!*********** THIS IS THE LONG+SHORT VERSION *****************************
! i.e. complete or partial data header saved !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!************************************************************************
!************* READ IN THE k-compressed files ***************************
!************************************************************************

! this reads in the main header information
      SUBROUTINE readmainheaderLONG(iIOUN,rFr1,rFr2,iLorS)

      IMPLICIT NONE 

      include 'convolve.param'

! iLorS      = which version data header was saved in
! iIOUN      = file number
! iSetLow    = integer index (between 1 and 89) of kcomp files start point
! iSetHigh   = integer index (between 1 and 89) of kcomp files end point
! rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
      INTEGER iIOUN,iLorS
      REAL rFr1,rFr2
      
! iCKD       = water continuum version
! caComment  = comment assocaiated with data run
! iNumLayers = number of layers in each gas profile
      CHARACTER*80 caVers
      CHARACTER*80 caComment
      INTEGER iCKD,iI,iNumLayers,iSetLow,iSetHigh,ia6(6),versNum,iaJunk(10)
      INTEGER iKMaxUserSet,char2num,iNumParams
      CHARACTER*3 ca3
      REAL rVers,rTen,pProf(kProfLayer+1),rr,raParams(kMaxUserSet)

      read(iIOUN) caVers
!      print*,'123456789012345678901234567890123456789012345678901234567890123456789012'
!      print *,'0         1         2         3         4         5         6         7         '
!      print *,caVers
! few lines later vcaVers is parsed for eg 1.03, 1.22 etc

      read(iIOUN) iNumLayers

      read(iIOUN) iKmaxUserSet
      IF (iKmaxUserSet .NE. kMaxUserSet) THEN
        print *,'Whoops! wrong reader version!! have different number'
        print *,'of parameters in *PARAMS',iKmaxUserSet,kMaxUserSet
        STOP
        END IF
      read(iIOUN) (raParams(iI),iI=1,iKmaxUserSet)

! read comment
      read(iIOUN) caComment

! read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh

!       print *,iNumLayers,iKmaxUserSet,rFr1,rFr2,iSetLow,iSetHigh

! read how data was saved
      read(iIOUN) iLorS
      if (iLorS .EQ. 0) THEN
        print *,'iLorS = 0, use readkcBasic.f'
        STOP
      END IF

! process STRING vx.xx into REAL
      ca3(1:1)=caVers(2:2)
      ca3(2:3)=caVers(4:5)
      rVers=0.0
      rTen=1.0
      DO iI=1,3
        rr=char2num(ca3,iI)
        rVers=rVers+rr/rTen
        rTen=rTen*10.0
        END DO
      print *,'kCARTA version = ',rVers
      if ((rVers .GE. 1.03) .and. (rVers .LE. 1.05)) THEN
        print *,'reading in MsubLayer, MThickLayer info ....'
        read(iIOUN) (ia6(iI),iI=1,6)
        read(iIOUN) (pProf(iI),iI=1,kProfLayer+1)   !!wrong hehehe
      elseif ((rVers .GE. 1.06) .and. (rVers .LE. 1.08)) THEN
        print *,'reading in MsubLayer, MThickLayer info ....'
        read(iIOUN) (ia6(iI),iI=1,6)
        read(iIOUN) (pProf(iI),iI=1,kProfLayer)   !!fixed by Larry McMillin
      elseif (rVers .GE. 1.09) THEN
        read(iIOUN) (pProf(iI),iI=1,kProfLayer)   !!fixed by Larry McMillin
        if (rVers .GE. 1.18) THEN
          !!! read iaaParams
          read(iIOUN) (iaJunk(iI),iI=1,10)
          read(iIOUN) (iaJunk(iI),iI=1,10)
          read(iIOUN) (iaJunk(iI),iI=1,10)
        end if
      end if

!      print *,'read in mainheader .......'

      RETURN
      END
!************************************************************************
! this reads in the main header information for iLorS = 0 (readkcBasic.f)
      SUBROUTINE readmainheaderSHORT(iIOUN,rFr1,rFr2,rDelta,iSetLow,iSetHigh, &
                      iTotalStuff,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut)

      IMPLICIT NONE 

      include 'convolve.param'
! iSetLow,iSetHigh = Chunk start, stop index
! rDelta     = point spacing
! iIOUN      = file number
! rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
! iChunks    = number of 10000 pr chunks that kCARTA deemed to process
! iTotal     = number of outputs per 10000 points
      INTEGER iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut
      INTEGER iIOUN,iLorS,iChunks,iTotalStuff,iSetLow,iSetHigh
      REAL rFr1,rFr2,raParams(kMaxUserSet),rDelta
      CHARACTER*80 caVers

! local variables
! iSetLow    = integer index (between 1 and 89) of kcomp files start point
! iSetHigh   = integer index (between 1 and 89) of kcomp files end point
! iLorS      = which version data header was saved in
! iCKD       = water continuum version
! caComment  = comment assocaiated with data run
! iNumLayers = number of layers in each gas profile
      CHARACTER*80 caComment
      INTEGER iCKD,iI,iNumLayers,ia6(6)
      INTEGER iKMaxUserSet,char2num
      CHARACTER*3 ca3
      REAL rVers,rTen,pProf(kProfLayer+1),rr,rBlocks

      read(iIOUN) caVers
!      print*,'123456789012345678901234567890123456789012345678901234567890123456789012'
!      print *,'0         1         2         3         4         5         6         7         '
!      print *,caVers
! few lines later vcaVers is parsed for eg 1.03, 1.22 etc

      read(iIOUN) iNumLayers

      read(iIOUN) iKmaxUserSet
      IF (iKmaxUserSet .NE. kMaxUserSet) THEN
        print *,'Whoops! wrong reader version!! have different number'
        print *,'of parameters in *PARAMS'
        STOP
        END IF
      read(iIOUN) (raParams(iI),iI=1,iKmaxUserSet)

! read comment
      read(iIOUN) caComment

! read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh
! read how data was saved
      read(iIOUN) iLorS
      IF (iLorS .NE. 0) THEN
        print *,'Sorry .. use readkcarta.f to read in this file .... '
        STOP
      END IF

! read frequency step size
      read(iIOUN) rdelta
! read chunk size = rDelta*10000
      read(iIOUN) rBlocks
! read how many output options saved per 10000 point chunk
      read(iIOUN) iTotalStuff
      read(iIOUN) iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut

      iChunks = iSetHigh - iSetLow + 1   !!!number of 10000 pt chunks to read

!      print *,'read in mainheader .......'

      RETURN
      END

!************************************************************************
! this combines Long Short
      SUBROUTINE readmainheaderLorS(iIOUN,rFr1,rFr2,iLorS, &
                      iChunks,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut, &
                      rDelta,iSetLow,iSetHigh)

      IMPLICIT NONE 

      include 'convolve.param'
! iSetLow,iSetHigh = Chunk start, stop index
! rDelta     = point spacing
! iIOUN      = file number
! rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
! iChunks    = number of 10000 pr chunks that kCARTA deemed to process
! iTotal     = number of outputs per 10000 points
      INTEGER iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut
      INTEGER iIOUN,iLorS,iChunks,iTotalStuff,iSetLow,iSetHigh
      REAL rFr1,rFr2,raParams(kMaxUserSet),rDelta
      CHARACTER*80 caVers

! local variables
! iSetLow    = integer index (between 1 and 89) of kcomp files start point
! iSetHigh   = integer index (between 1 and 89) of kcomp files end point
! iLorS      = which version data header was saved in
! iCKD       = water continuum version
! caComment  = comment assocaiated with data run
! iNumLayers = number of layers in each gas profile
      CHARACTER*80 caComment
      INTEGER iCKD,iI,iNumLayers,ia6(6),iaJunk(10)
      INTEGER iKMaxUserSet,char2num
      CHARACTER*3 ca3
      REAL rVers,rTen,pProf(kProfLayer+1),rr,rBlocks

      iChunks = -1
      iTotalStuff = -9999
      iNumGasPathsOut = -9999
      iNumMixPathsOut = -9999
      iNumRadsOut = -9999
      rDelta = -9999.0
      iSetLow = -9999
      iSetHigh = -9999

      read(iIOUN) caVers
!      print*,'123456789012345678901234567890123456789012345678901234567890123456789012'
!      print *,'0         1         2         3         4         5         6         7         '
!      print *,caVers
! few lines later vcaVers is parsed for eg 1.03, 1.22 etc

      read(iIOUN) iNumLayers

      read(iIOUN) iKmaxUserSet
      IF (iKmaxUserSet .NE. kMaxUserSet) THEN
        print *,'Whoops! wrong reader version!! have different number'
        print *,'of parameters in *PARAMS'
        STOP
        END IF
      read(iIOUN) (raParams(iI),iI=1,iKmaxUserSet)

! read comment
      read(iIOUN) caComment

! read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh
! read how data was saved
      read(iIOUN) iLorS

      IF (iLorS .EQ. 0) THEN
        print *,'this is SHORT version'

! read frequency step size
        read(iIOUN) rdelta
! read chunk size = rDelta*10000
        read(iIOUN) rBlocks
! read how many output options saved per 10000 point chunk
        read(iIOUN) iTotalStuff
        read(iIOUN) iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut

        iChunks = iSetHigh - iSetLow + 1   !!!number of 10000 pt chunks to read

      ELSE
! process STRING vx.xx into REAL
        print *,'this is LONG version'
        ca3(1:1)=caVers(2:2)
        ca3(2:3)=caVers(4:5)
        rVers=0.0
        rTen=1.0
        DO iI=1,3
          rr=char2num(ca3,iI)
          rVers=rVers+rr/rTen
          rTen=rTen*10.0
          END DO
        print *,'kCARTA version = ',rVers
        if ((rVers .GE. 1.03) .and. (rVers .LE. 1.05)) THEN
          print *,'reading in MsubLayer, MThickLayer info ....'
          read(iIOUN) (ia6(iI),iI=1,6)
          read(iIOUN) (pProf(iI),iI=1,kProfLayer+1)   !!wrong hehehe
        elseif ((rVers .GE. 1.06) .and. (rVers .LE. 1.08)) THEN
          print *,'reading in MsubLayer, MThickLayer info ....'
          read(iIOUN) (ia6(iI),iI=1,6)
          read(iIOUN) (pProf(iI),iI=1,kProfLayer)   !!fixed by Larry McMillin
        elseif (rVers .GE. 1.09) THEN
          read(iIOUN) (pProf(iI),iI=1,kProfLayer)   !!fixed by Larry McMillin
          if (rVers .GE. 1.18) THEN
            !!! read iaaParams
            read(iIOUN) (iaJunk(iI),iI=1,10)
            read(iIOUN) (iaJunk(iI),iI=1,10)
            read(iIOUN) (iaJunk(iI),iI=1,10)
          end if
        end if
      END IF
       
      RETURN
      END

!************************************************************************

! this function takes a character in a string, and gives Roman numeral
      INTEGER FUNCTION char2num(ca3,iI)

      IMPLICIT NONE 
       
      CHARACTER*3 ca3
      INTEGER iI
      
      CHARACTER*1 c1
      INTEGER iJ,iFound,iK

      iK=ICHAR('0')
      iFound = -1
      c1=ca3(iI:iI)
      
      iJ=32
 10   CONTINUE
      IF ((char(iJ) .NE. c1) .AND. (iJ .LT. 128)) THEN
        iJ=iJ+1
        GOTO 10
        END IF

      IF (char(iJ) .EQ. c1) THEN
        iFound = 1
        END IF

      IF (iFound .LT. 0) THEN
        print *,'cannot find ASCII numeral for iI = ',iI,' in ',ca3
        STOP
        END IF

      char2num=iJ-iK

      RETURN
      END
        
!************************************************************************
! this function reads in the individual gas paths, and list of 
! paths to be output
      SUBROUTINE readgaspaths(iIOUN,iNumPathsOut,iLorS)

      IMPLICIT NONE 

      include 'convolve.param'

! iIOUN          = file ioUNIT number
! iNumPathsOut   = number of single gas paths to be output
! iLorS          = Long or Short header format
      INTEGER iIOUN,iNumPathsOUt,iLorS

! iaOutPaths     = array containing the paths to be output
! raaAmt,raaTemp = matrices containing the gas amts,temperatures
! iaGasID        = array containing ID's of gases (kcomp and xsec)
! iNumLayers     = number of layers in each gas profile
! iNumPaths      = total number of gas paths == iNumGases*iNumLayers
      INTEGER iNumLayers,iNumPaths,iaGasId(kGasStore)
      REAL rAmt,rTemp
      INTEGER iPath,iGasID,iNumGases,ii,jj,iaOutPaths(kPathsOut)
      REAL raaTemp(kProfLayer,kGasStore),raaAmt(kProfLayer,kGasStore)

      iNumLayers=kProfLayer
! read iNumPaths=iNumGases*iNumLayers
! -------- optional due to iLorS ----------------
      IF (iLorS .GT. 0) THEN
        read(iIOUN) iNumPaths
        iNumGases=iNumPaths/iNumLayers

! now read the PathNum GasID Temperature Amount for each path in the profile
        DO ii=1,iNumGases
          DO jj=1,iNumLayers
            read(iIOUN) iPath,iGasID,rTemp,rAmt
            raaAmt(jj,ii)=rAmt
            raaTemp(jj,ii)=rTemp
            END DO 
         iaGasID(ii)=iGasID
         END DO
       END IF
! -------- end optional due to iLorS ----------------

! read in Number of Paths to be output each time
       read(iIOUN) iNumPathsOut
       IF (iLorS .LT.  0) THEN
         iNumPaths=iNumPathsOut
         END IF

!       print *,'output the following paths',iNumPathsOut
       IF ((iNumPathsOut .LT. 0).OR.(iNumPathsOut .GT. iNumPaths)) THEN
          print *,'Mistake in number of paths to be output!!!'
          STOP
          END IF
      IF (iNumPathsOut .GT.  0) THEN
! read in paths
        read(iIOUN) (iaOutPaths(ii), ii=1,iNumPathsOut)
        END IF

!       print *,'the total number of gas paths is ',iNumPaths
!       print *,'which gives num of gases',iNumGases
!       print *,'Read in the Path nums, gasIDs, amts, temps'

      RETURN
      END

!************************************************************************
! now read in the mixed path info
      SUBROUTINE readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS)

      IMPLICIT NONE 

      include 'convolve.param'

! iIOUN          = file IOUNIT number
! iNpmix         = number of mixed paths
! iNumMixPathsOut= number of mixed paths to be output
! iLorS          = long or short data header
      INTEGER iIOUN,iNpmix,iNumMixPathsOut,iLorS

! local variables
! iaOutMixPaths  = integer array containing list of paths to be output
! iMixFileLines  = number of lines containing character description of MP's
! caaMixedPathInfo == character description of mixed path info
! raMixTemp      = array containing the Mixed Path Temperatures
      INTEGER iaOutMixPaths(kPathsOut),iMixFileLines
      REAL raMixTemp(kMixFilRows)
      CHARACTER*130 caaMixedPathInfo(kMixFilRows)
      INTEGER ii,kk
      CHARACTER*130 caMixedPathInfo

! read in number of mixed paths in *MIXFIL
      read(iIOUN) iNpmix

! -------- optional due to iLorS ----------------
      IF (iLorS .GT. 0) THEN
! read in #of lines that have info in *MIXFIL
        read(iIOUN) iMixFileLines

        IF  (iNpmix .GT. 0) THEN
! read in the mixing table
          DO ii=1,iMixFileLines
            read(iIOUN) caMixedPathInfo
!           caaMixedPathInfo(ii,1:130)=caMixedPathInfo(1:130)
            END DO

! read in mixed path temperatures
          read(iIOUN) (raMixTemp(ii), ii=1,iNpmix)
          END IF
        END IF
! -------- end optional due to iLorS ----------------

! read in Number of mixed paths to be output each time
      IF (iNpmix .GT. 0) THEN
        read(iIOUN) iNumMixPathsOut
        IF ((iNumMixPathsOut.LT.0).OR.(iNumMixPathsOut.GT.iNpmix)) THEN
            print *,'Mistake in number of mixed paths to be output!!!'
            STOP
            END IF

        IF (iNumMixPathsOut .GT. 0) THEN
! read in paths
          read(iIOUN) (iaOutMixPaths(ii),ii=1,iNumMixPathsOut)
          END IF

        END IF

!      print *,'Number of Mixed Paths = ',iNpmix
!      print *,'Number of lines in MIXFIL to read = ',iMixFileLines

      RETURN
      END

!************************************************************************  
! this subroutine reads in the atmosphere information
      SUBROUTINE readatmospheres(iIOUN,iNpmix,iNatm,iaNumLayersInAtm,iaNumLayersOut)

      IMPLICIT NONE 

      include 'convolve.param'

! iIOUN         = file unit number
! iNpmix        = number of mixed paths (from readmixedpaths)
! iNatm         = number of atmospheres
! iaNumLayersInAtm = integer array with number of layers in each atmosphere
! iaNumLayersOut= integer array with number of layers out per atmosphere
      INTEGER iaNumLayersOut(kMaxAtm),iaNumLayersInAtm(kMaxAtm)
      INTEGER iIOUN,iNpmix,iNatm

! local variables
! raaOutPress   = matrix containing pressures where radiances are output
! iaAtmNum      = integer array with atmospher number (1..iNatm)
! iaaRadLayer   = for each atmosphere, this lists which layers build up the atm
! iaaOutRadPaths= for each atmosphere, list which layers are to be output
! raTbdy,raTinit,raSat = for each atmosphere, these arrays contain atm temp,
!                          surface temp and satellite viewing angle
! raaaEms      =emissivity
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! rakSolarRefl   =solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off      
! raheight = height of instrument above earth
      INTEGER iNumPathsAtmos,iAtm,iNumLayersOut,ii,jj,kk
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),raHeight(kMaxAtm)
      INTEGER iakThermal(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raaaEms(kMaxAtm,kEmsRegions,2)
      REAL raTbdy(kMaxAtm),raTinit(kMaxAtm),raSat(kMaxAtm)
      REAL raaOutPress(kMaxAtm,kProfLayer)
      INTEGER iaAtmNum(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iaEms(kMaxAtm)
      INTEGER iaaOutRadPaths(kMaxAtm,kPathsOut)

!      print *,' reading in the atmospheres info .... '

! read in number of atmospheres in *RADFIL
      read(iIOUN) iNatm
! read in kEmsRegions
      read(iIOUN) ii                       
      IF (ii .NE. kEmsRegions) THEN
        print *,'please recompile code ... kEmsRegions incorrect'
        STOP
        END IF

!      print *,'There are atmospheres in RADFIL ',iNatm

      IF (iNatm .GT. 0) THEN
        DO ii=1,iNatm
         iaAtmNum(ii)=ii
         iaNumLayersinAtm(ii)=0
         iaNumLayersOut(ii)=0
         iaEms(ii)=0
         END DO

       DO ii=1,iNatm
         DO kk=1,iNpmix
           iaaOutRadPaths(ii,kk)=0
           END DO
         END DO

       DO ii=1,iNatm
         DO kk=1,kProfLayer
           raaOutPress(ii,kk)=0.0
           END DO
         END DO

       DO ii=1,iNatm
!         print *,'reading in info for atm num ',ii

! read in atmosphere number, # of mixed paths in atmosphere
         read(iIOUN) iAtm,iNumPathsAtmos
         iaNumLayersInAtm(ii)=iNumPathsAtmos
!         print *,'Num of Paths in Atmosphere = ',iNumPathsAtmos

!         IF((iNumPathsAtmos.LT.0).OR.(iNumPathsAtmos.GT.iNpmix)) THEN
         IF((iNumPathsAtmos.LT.0).OR.(iNumPathsAtmos.GT.kProfLayer)) THEN
           print *,'Mistake in num of mixed paths used in atmosphere !'
           STOP
           END IF
        IF (iNumPathsAtmos .GT.  0) THEN
!          print *,'Reading the Mixed Path numbers ...'
! read in paths making up atmosphere
          read(iIOUN) (iaaRadLayer(ii,jj),jj=1,iNumPathsAtmos) 
          read(iIOUN) raTbdy(ii),raTinit(ii),raSat(ii),raHeight(ii)
          read(iIOUN) iakSolar(iI),rakSolarAngle(iI),rakSolarRefl(iI),iakThermal(iI),rakThermalAngle(iI),iakThermalJacob(iI)
          read(iIOUN) iaEms(ii)
          DO jj=1,iaEms(ii)
            read(iIOUN) raaaEms(ii,jj,1),raaaEms(ii,jj,2)
            END DO

! now read in layers to be output
          read(iIOUN) iNumLayersOut
          iaNumLayersOut(ii)=iNumLayersOut
!          print *,'This atm has layers to be output',iNumLayersOut

      IF((iNumLayersOut.LT.0).OR.(iNumLayersOut.GT.iNumPathsAtmos))THEN 
        print *,'Mistake in number of layers to be output for Atm#!',ii
        STOP
        END IF
          IF (iNumLayersOut .GT. 0) THEN
! read in mixed paths to be output
            read(iIOUN) (iaaOutRadPaths(ii,jj),jj=1,iNumLayersOut)
            read(iIOUN) (raaOutPress(ii,jj),jj=1,iNumLayersOut)
!            print *,(iaaOutRadPaths(ii,jj),jj=1,iNumLayersOut)
            END IF
          END IF

          END DO
        END IF

      RETURN
      END
!************************************************************************
! this prints out a summary of the saved paths/MP/layers
      SUBROUTINE printsummary(iaOutPaths,iaOutMixPaths,iNatm,iaNumLayersOut,iaaOutRadPaths,ii1,ii2)

      IMPLICIT NONE 

      include 'convolve.param'
! iaOutPaths     = integer array with list of gas   paths to be output
! iaOutMixPaths  = integer array with list of mixed paths to be output
! iNatm          = number of atmospheres
! iaNumLayersOut = for each atmosphere, how many layers to be output
! iaaOutRadPaths = for each atmosphere, list of layers to be output
! ii1            =iNumPathsOut
! ii2            =iNumMixPathsOut
      INTEGER iaOutPaths(kPathsOut),iaOutMixPaths(kPathsOut),iNatm
      INTEGER iaNumLayersOut(kMaxAtm)
      INTEGER iaaOutRadPaths(kMaxAtm,kPathsOut)
      INTEGER ii1,ii2

! local variables
      INTEGER ii,iit,ind(kPathsOut),jj,kk

      print *,'summarizing output info ....'
      print *,'option#       path/MP/layer'
      print *,' ---------------------------'

      DO ii=1,ii1 
        ind(ii)=ii
        END DO
 
      print *,'num of paths to be output = ',ii1
      DO ii=1,ii1
        print *,ii,iaOutPaths(ii)
        END DO

      DO ii=1,ii2
        ind(ii)=ii1+ii
        END DO  
      print *,'Num of mixed paths to be output : ',ii2
      DO ii=1,ii2
        print *,ii,iaOutMixPaths(ii)
        END DO

      print *,'num of atmosphers in RADFIL = ',iNatm
      DO ii=1,iNatm
        iit=iaNumLayersOut(ii)
        kk=0
        DO jj=1,ii-1 
          kk=kk+iaNumLayersOut(jj)
          END DO
        DO jj=1,iit
          ind(jj)=ii1+ii2+kk+jj
          END DO
        print *,'Atm ',ii,' has ',iit,' layers to be output '
        DO jj=1,iaNumLayersOut(ii)
          print *,ind(jj),iaaOutRadPaths(ii,jj)
          END DO
        END DO

      RETURN
      END

!************************************************************************
! this reads in a "block" of data
      SUBROUTINE readblock(iFrLow,iFrHigh,rL,rH,ra,iIOUN,iBlSz)

      IMPLICIT NONE 

      include 'convolve.param'

! iFrLow,iFrHigh = integer index of output start/stop point (1--10000)
! rL,rH          = freq start/stop points, corresponding to above two indices
! iBlSz          = block size == number of points to be read in (400 or 10000)
! iIOUN          = file IOUNIT number
! ra             = real array containing the data that is read in
      INTEGER iFrLow,iFrHigh,iIOUN,iBlSz,ii
      REAL rL,rH,ra(kMaxPts)

! first read the header
      read(iIOUN) iFrLow,iFrHigh,rL,rH

! now read the data
      read(iIOUN) (ra(ii),ii=1,iBlSz)

      RETURN
      END
!************************************************************************
! this subroutine reads in the data
      SUBROUTINE readdata(iIOUN,iDataSize,raEntire,iDataPoints,iTotal)

      IMPLICIT NONE 

      include 'convolve.param'

! iDataSize=number of bytes from start of data in this block, 
!           to same start of data
!           in next kcomp block (so we can directly jump there)
      INTEGER iIOUN,iDataSize,iDataPoints,iTotal
      REAL raEntire(kMaxEntire)

      REAL raTemp(kMaxPts)
      INTEGER iI,iJ,iDummy

      iDataPoints=0
      DO iI=1,iTotal      
        read(iIOUN) (raTemp(iJ),iJ=1,kMaxPts)

        DO iJ=1,kMaxPts
          iDataPoints=iDataPoints+1
          raEntire(iDataPoints)=raTemp(iJ)
          END DO

        IF (iI .NE. iTotal) THEN
          CALL fseek_local(iIOUN,idataSize,1)       !!!!if using g77
          END IF

        END DO

      RETURN
      END

!************************************************************************
! this subroutine reads in the data for iLorS == 0
      SUBROUTINE ReadData2S(iIOUN,iJ,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut,raFreq,raaEntire,rDelta,rFrLow)

      IMPLICIT NONE 

      include 'convolve.param'

! iDataSize=number of bytes from start of data in this block, 
!           to same start of data
!           in next kcomp block (so we can directly jump there)
!!! typically iTotalChunks = 89 (chunks) 
      INTEGER iIOUN,iJ,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut    
      REAL raaEntire(kMaxEntire,100),raFreq(kMaxEntire),rDelta,rFrLow

      INTEGER iMainType,iSubMainType,iNumberOut,ikMaxPts

      REAL raTemp(kMaxPts)
      INTEGER iI,iCnt,iX,iDummy
      INTEGER iaIndX(kMaxPts),iaInd(kMaxPts)

      iaIndX = (/(iI,iI=1,kMaxPts)/)
      iaInd = iaIndX + (iJ-1)*10000

      raFreq(iaInd) = rFrLow + (iaIndX-1)*rDelta
      rFrLow = rFrLow + rDelta*kMaxPts

      iCnt = 0

      IF (iNumGasPathsOut .GT. 0) THEN
        DO iI=1,iNumGasPathsOut
          iCnt = iCnt + 1
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      IF (iNumMixPathsOut .GT. 0) THEN
        DO iI=1,iNumMixPathsOut
          iCnt = iCnt + 1
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      IF (iNumRadsOut .GT. 0) THEN
        DO iI=1,iNumRadsOut
          iCnt = iCnt + 1 
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      RETURN
      END

!************************************************************************
! this subroutine reads in the data for iLorS > 0
      SUBROUTINE ReadData2L(iIOUN,iJ,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut,raFreq,raaEntire)

      IMPLICIT NONE 

      include 'convolve.param'

! iDataSize=number of bytes from start of data in this block, 
!           to same start of data
!           in next kcomp block (so we can directly jump there)
!!! typically iTotalChunks = 89 (chunks) 
      INTEGER iIOUN,iJ,iNumGasPathsOut,iNumMixPathsOut,iNumRadsOut    
      REAL raaEntire(kMaxEntire,100),raFreq(kMaxEntire)

      INTEGER iMainType,iSubMainType,iNumberOut,ikMaxPts
      REAL rFrLow,rFrHigh,rDelta

      REAL raTemp(kMaxPts)
      INTEGER iI,iCnt,iX,iDummy
      INTEGER iaIndX(kMaxPts),iaInd(kMaxPts)

      iaIndX = (/(iI,iI=1,kMaxPts)/)
      iaInd = iaIndX + (iJ-1)*10000

      iCnt = 0

      IF (iNumGasPathsOut .GT. 0) THEN
        READ(iIOUN) iMainType,iSubMainType,iNumberOut
        READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta
!       write(*,'(5(I10),3(F20.12))') iJ,iNumGasPathsOut,iMainType,iSubMainType,iNumberOut,rFrLow,rFrHigh,rDelta
        IF (iJ .EQ. 1) THEN
          write(*,'(A,2(I10),1(F20.12))') 'GASOD ',iJ,iNumGasPathsOut,rFrLow
        END IF
        raFreq(iaInd) = rFrLow + (iaIndX-1)*rDelta
        DO iI=1,iNumGasPathsOut
          iCnt = iCnt + 1
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      IF (iNumMixPathsOut .GT. 0) THEN
        READ(iIOUN) iMainType,iSubMainType,iNumberOut
        READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta
!       write(*,'(5(I10),3(F20.12))') iJ,iNumMixPathsOut,iMainType,iSubMainType,iNumberOut,rFrLow,rFrHigh,rDelta
        IF (iJ .EQ. 1) THEN
          write(*,'(A,2(I10),1(F20.12))') 'MIXOD ',iJ,iNumMixPathsOut,rFrLow
        END IF
        raFreq(iaInd) = rFrLow + (iaIndX-1)*rDelta
        DO iI=1,iNumMixPathsOut
          iCnt = iCnt + 1
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      IF (iNumRadsOut .GT. 0) THEN
        DO iI=1,iNumRadsOut
          iCnt = iCnt + 1 
          READ(iIOUN) iMainType,iSubMainType,iNumberOut
          READ(iIOUN) ikMaxPts,rFrLow,rFrHigh,rDelta
!          write(*,'(5(I10),3(F20.12))') iJ,iNumRadsOut,iMainType,iSubMainType,iNumberOut,rFrLow,rFrHigh,rDelta
          IF ((iI .EQ. 1) .AND. (iJ .EQ. 1)) THEN
            write(*,'(A,2(I10),1(F12.5),A,I3)') 'RADS  ',iJ,iNumRadsOut,rFrLow,' need to read ',iNumRadsOut
          END IF
          raFreq(iaInd) = rFrLow + (iaIndX-1)*rDelta
          read(iIOUN) (raTemp(iX),iX=1,kMaxPts)
          raaEntire(iaInd,iCnt) = raTemp
        END DO
      END IF

      IF (iJ .GT. 1) THEN
        write(*,'(A,I3,A,F12.5)') 'chunk iJ = ',iJ,' has start freq rFrLow = ',raFreq(1)
      END IF

      RETURN
      END

!************************************************************************
! this subroutine does the command line stuff
      SUBROUTINE DoCommandLine(caInName,caOutName,iTextOrBinary,iWhich)

      IMPLICIT NONE 

      CHARACTER*80 caInName   !!!Input binary file, written out by kcarta
      CHARACTER*80 caOutName  !!!Output file, written out by this code
                              !!!if only 2 args, this is binary file
                              !!!if      3 args, this is text file
      INTEGER iTextOrBinary   !!! -1 if text   file (3 arguments)
                              !!! +1 if binary file (2 arguments)
      INTEGER iWhich          !!! if text file, which one to output 
 
      CHARACTER*80 caWhich 
      INTEGER iDummy,iError 
 
!      print *,'Enter INPUT binary file name '
!      read 1000,caInName
!      print *,'Do you want TEXT (-1) or BINARY (+1) output?'
!      read *,iTextOrBinary
!      IF (iTextOrBinary .EQ. -1) THEN
!        print *,'input file has ',iTotalStuff,' outputs'
!        print *,'which one do you want output???'
!        read *,iWhich
!        END IF
!      print *,'Enter OUTPUT file name '
!      read 1000,caOutName

! this is the number of args
      INTEGER iarg! 

      !use command line stuff 
      iDummy = iargc() 
 
      IF (iDummy .EQ. 2) THEN
        iWhich = -1
        iTextOrBinary = +1          !!!!binary file
        END IF
      IF (iDummy .EQ. 3) THEN
        iWhich = -1
        iTextOrBinary = -1          !!!!text file
        END IF
      
      IF (iDummy .LT. 2) THEN  
        print *, 'less than two arguments in command line is NOT allowed' 
        STOP
        END IF 

      IF (iDummy .GT. 3) THEN  
        print *, 'more than three arguments in command line is NOT allowed' 
        STOP
        END IF 

      DO iError=1,iDummy 
        IF (iError .EQ. 1) THEN 
          CALL getarg(1,caInName) 
          END IF 
 
        IF (iError .EQ. 2) THEN 
          CALL getarg(2,caOutName) 
          END IF 
 
        IF (iError .EQ. 3) THEN 
          CALL getarg(3,caWhich) 
          CALL GetNumber(caWhich,iWhich)
          END IF
         END DO


      RETURN
      END
!************************************************************************
      SUBROUTINE GetNumber(caName,iW)

      CHARACTER*80 caName
      INTEGER iW

      CHARACTER*80 caTempName
      INTEGER iR,iL,iInt,iFactor,iI
      CHARACTER c

! find the "right" length of the input root name
      iR=len(caName)
 11   continue
      IF (caName(iR:iR) .eq. ' ') THEN
        iR=iR-1
        GO TO 11      
        END IF

! find the "left" length of the input root name
      iL=1
 12   continue
      IF (caName(iL:iL) .eq. ' ') THEN
        iL=iL+1
        GO TO 12      
        END IF

! thus the entered word exists between iL:iR
      IF (iL .GT. iR) THEN
        print *, 'Fatal Error! Invalid (blank) string in file !'
        STOP
        END IF

! now rearrange the string, so that it is right padded with blanks
! this is basically equivalent to  ADJUSTR(caName)
      DO iInt=1,80
        caTempName(iInt:iInt)=' '
        END DO
      caTempName(1:iR-iL+1)=caName(iL:iR)
      caName(1:80)=caTempName(1:80)

      iW = 0
      iFactor = 1
      DO iInt = iR-iL,1,-1
        iFactor = 10*iFactor
        END DO

      DO iInt = 1,iR-iL+1
        ! = caName(iInt:iInt)        
        READ(c,13) iI
        iW = iW + iFactor*iI
        iFactor = iFactor/10
        END DO

 13   FORMAT(I1)

      RETURN
      END


!************************************************************************
! this is the fseek function, as defined by PGI F77

! fseek
! Position file at offset.
! Synopsis
!
! integer function fseek(lu, offset, from)
! integer lu
! integer offset
! integer from
!
! fseek repositions a file connected to logical unit lu. offset is an 
! offset in bytes relative to the position specified by "from" = 0,1,2 :
!
!   0    beginning of the file 
!   1    current position 
!   2    end of the file

! If successful, the value returned by fseek will be zero; otherwise, 
! it's a system error code. 

! USE IFPORT
! https://ahamodel.uib.no/intel/GUID-118CDDD3-D29C-4691-89BC-AFB077AA55CD.html

      SUBROUTINE fseek_local(iIOUN,iStep,iFrom)

#ifdef IFORT
    USE IFPORT
#endif

      IMPLICIT NONE 

      include 'convolve.param'

      INTEGER iIOUN,iStep
      INTEGER iDebug,iI,iLoop
!      INTEGER fseek,ftell   !! absoft and pgf treat these as fcns, g77 as subr
      INTEGER*8 iX8
      CHARACTER cC
      REAL raX(kMaxPts*100),rX

      integer(4) istat4, offset4, ipos4,iIOUN4
      integer iFrom,iX,iFileWhere,my_pos

      iFileWhere = ftell(iIOUN)
      write(*,'(A,I4,I8,I2,I12)') 'just entered fseek_local ; iIOUN,iStep,iFrom .. iFileWhere = ',iIOUN,iStep,iFrom,iFileWhere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this is to test if things work
! stopped working on Absift v10.0
! seems to work on PGI F77
      iIOUN4 = iIOUN

      iDebug = -1
      IF (iDebug  .GT. 0) THEN

        offset4 = 0
        ipos4 = 1

        istat4 = fseek(iIOUN4,offset4,ipos4)
!        call fseek(iIOUN,offset,ipos,istat)
        write(*,'(A,I3,I10,F12.5)') 'istat (if 0 thats good) = ',istat4
        iFileWhere = ftell(iIOUN)
        read(iIOUN) rX
        write(*,'(A,I3,I10,F12.5)') '--------- A1   << DEBUG >> iIOUN,iFileWhere,rX = ',iIOUN,iFileWhere,rX

        offset4 = offset4 + 10
        istat4 = fseek(iIOUN4,offset4,ipos4)
!        call fseek(iIOUN,offset,ipos,istat)
        write(*,'(A,I3,I10,F12.5)') 'istat (if 0 thats good) = ',istat4
        iFileWhere = ftell(iIOUN)
        read(iIOUN) rX
        write(*,'(A,I3,I10,F12.5)') '--------- A2   << DEBUG >> iIOUN,iFileWhere,rX = ',iIOUN,iFileWhere,rX

        offset4 = offset4 + 20
        istat4 = fseek(iIOUN4,offset4,ipos4)
!        call fseek(iIOUN,offset,ipos,istat)
        write(*,'(A,I3,I10,F12.5)') 'istat (if 0 thats good) = ',istat4
        iFileWhere = ftell(iIOUN)
        read(iIOUN) rX
        write(*,'(A,I3,I10,F12.5)') '--------- A3   << DEBUG >> iIOUN,iFileWhere,rX = ',iIOUN,iFileWhere,rX

      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!      CALL fseek(iIOUN,iStep,iFrom,iX)   !!if built in g77 works

      IF (iFrom .EQ. 2) THEN
        print *,'oops cannot handle fseek(iU,iStep,2)'
        STOP
      END IF

      iFileWhere = ftell(iIOUN)
      print *,'before possible rewind ... iFileWhere,iFrom = ',iFileWhere,iFrom
      IF (iFrom .EQ. 0) THEN
        rewind(iIOUN)
      END IF
      iFileWhere = ftell(iIOUN)
      print *,'after possible rewind ... iFileWhere,iFrom = ',iFileWhere,iFrom
      
!      inquire(iIOUN, POS=my_pos)
!      print *,'my_pos = ',my_pos

!      IF (iFrom .EQ. 1) THEN
!        DO iLoop = 1,iStep
!          iFileWhere = ftell(iIOUN)
!          print *,'----------loop >>>>>>><<<<<<',iIOUN,iLoop,iStep,iFileWhere,kMaxPts
!          read(iIOUN) (raX(iI),iI=1,kMaxPts)
!          print *,iLoop,raX(1),raX(kMaxPts)
!        END DO
!      END IF

      IF (iFrom .EQ. 1) THEN
        DO iLoop = 1,1
          iFileWhere = ftell(iIOUN)
          print *,'----------loop >>>>>>><<<<<<',iIOUN,iLoop,iStep,iFileWhere,kMaxPts*iStep
          read(iIOUN) (raX(iI),iI=1,kMaxPts*iStep)
          print *,iLoop,raX(1),raX(kMaxPts*iStep)
        END DO
      END IF

      RETURN
      END

!************************************************************************
