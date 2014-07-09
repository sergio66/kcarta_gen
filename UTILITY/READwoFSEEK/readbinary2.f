c*********** THIS IS THE LONG+SHORT VERSION *****************************
c i.e. complete or partial data header saved !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c************************************************************************
c************* READ IN THE k-compressed files ***************************
c************************************************************************
c this reads in the main header information
      SUBROUTINE readmainheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh,iLorS)

      include 'convolve.param'

c iLorS      = which version data header was saved in
c iIOUN      = file number
c iSetLow    = integer index (between 1 and 89) of kcomp files start point
c iSetHigh   = integer index (between 1 and 89) of kcomp files end point
c rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
      INTEGER iIOUN,iLorS
      REAL rFr1,rFr2,raParams(kMaxUserSet)
      CHARACTER*80 caVers

c iCKD       = water continuum version
c caComment  = comment assocaiated with data run
c iNumLayers = number of layers in each gas profile
      CHARACTER*80 caComment
      INTEGER iCKD,iI,iNumLayers,iSetLow,iSetHigh,ia6(6)
      INTEGER iKMaxUserSet,char2num
      CHARACTER*3 ca3
      REAL rVers,rTen,pProf(kProfLayer+1),rr

      read(iIOUN) caVers
      read(iIOUN) iNumLayers

      read(iIOUN) iKmaxUserSet
      IF (iKmaxUserSet .NE. kMaxUserSet) THEN
        print *,'Whoops! wrong reader version!! have different number'
        print *,'of parameters in *PARAMS'
        STOP
        END IF
      read(iIOUN) (raParams(iI),iI=1,iKmaxUserSet)

c read comment
      read(iIOUN) caComment
c read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh
c read how data was saved
      read(iIOUN) iLorS
      IF (iLorS .EQ. 0) THEN
        print *,'Sorry .. use readkcbasic.f to read in this file .... '
        STOP
        END IF

c process STRING vx.xx into REAL
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
        end if 
      print *,'read in mainheader .......' 

      RETURN
      END
c************************************************************************
c this reads in the main header information for iLorS = 0 (readkcBasic.f)
      SUBROUTINE readmainheaderSHORT(iIOUN,rFr1,rFr2,rDelta,iSetLow,iSetHigh,
     $                               iTotalStuff)

      include 'convolve.param'
c iSetLow,iSetHigh = Chunk start, stop index
c rDelta     = point spacing
c iIOUN      = file number
c rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
c iChunks    = number of 10000 pr chunks that kCARTA deemed to process
c iTotal     = number of outputs per 10000 points
      INTEGER iIOUN,iLorS,iChunks,iTotalStuff,iSetLow,iSetHigh
      REAL rFr1,rFr2,raParams(kMaxUserSet),rDelta
      CHARACTER*80 caVers

c local variables
c iSetLow    = integer index (between 1 and 89) of kcomp files start point
c iSetHigh   = integer index (between 1 and 89) of kcomp files end point
c iLorS      = which version data header was saved in
c iCKD       = water continuum version
c caComment  = comment assocaiated with data run
c iNumLayers = number of layers in each gas profile
      CHARACTER*80 caComment
      INTEGER iCKD,iI,iNumLayers,ia6(6)
      INTEGER iKMaxUserSet,char2num
      CHARACTER*3 ca3
      REAL rVers,rTen,pProf(kProfLayer+1),rr,rBlocks

      read(iIOUN) caVers
      read(iIOUN) iNumLayers

      read(iIOUN) iKmaxUserSet
      IF (iKmaxUserSet .NE. kMaxUserSet) THEN
        print *,'Whoops! wrong reader version!! have different number'
        print *,'of parameters in *PARAMS'
        STOP
        END IF
      read(iIOUN) (raParams(iI),iI=1,iKmaxUserSet)

c read comment
      read(iIOUN) caComment
c read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh
c read how data was saved
      read(iIOUN) iLorS
      IF (iLorS .NE. 0) THEN
        print *,'Sorry .. use readkcarta.f to read in this file .... '
        STOP
        END IF

c read frequency step size
      read(iIOUN) rdelta
c read chunk size = rDelta*10000
      read(iIOUN) rBlocks
c read how many output options saved per 10000 point chunk
      read(iIOUN) iTotalStuff

      iChunks = iSetHigh - iSetLow + 1   !!!number of 10000 pt chunks to read

c      print *,'read in mainheader .......'

      RETURN
      END
c************************************************************************
c this function takes a character in a string, and gives Roman numeral
      INTEGER FUNCTION char2num(ca3,iI)
       
      CHARACTER*3 ca3
      INTEGER iI
      
      CHARACTER*1 c1
      INTEGER iJ,iFound,iK

      iK=ICHAR('0')
      iFound = -1
      c1=ca3(iI:iI)
      
      iJ=32
 10      CONTINUE
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
        
c************************************************************************
c this function reads in the individual gas paths, and list of  
c paths to be output 
      SUBROUTINE readgaspaths(iIOUN,iNumPathsOut,iLorS) 
 
      include 'convolve.param' 
 
c iIOUN          = file ioUNIT number 
c iNumPathsOut   = number of single gas paths to be output 
c iLorS          = Long or Short header format 
      INTEGER iIOUN,iNumPathsOUt,iLorS 
 
c iaOutPaths     = array containing the paths to be output 
c raaAmt,raaTemp = matrices containing the gas amts,temperatures 
c iaGasID        = array containing ID's of gases (kcomp and xsec) 
c iNumLayers     = number of layers in each gas profile 
c iNumPaths      = total number of gas paths == iNumGases*iNumLayers 
      INTEGER iNumLayers,iNumPaths,iaGasId(kMaxGas) 
      REAL rAmt,rTemp 
      INTEGER iPath,iGasID,iNumGases,ii,jj,iaOutPaths(kPathsOut) 
      REAL raaTemp(kProfLayer,kMaxGas),raaAmt(kProfLayer,kMaxGas) 
 
      iNumLayers=kProfLayer 
c read iNumPaths=iNumGases*iNumLayers 
c -------- optional due to iLorS ---------------- 
      IF (iLorS .GT. 0) THEN 
        read(iIOUN) iNumPaths 
        iNumGases=iNumPaths/iNumLayers 
 
c now read the PathNum GasID Temperature Amount for each path in the profile 
        DO ii=1,iNumGases 
          DO jj=1,iNumLayers 
            read(iIOUN) iPath,iGasID,rTemp,rAmt 
            raaAmt(jj,ii)=rAmt 
            raaTemp(jj,ii)=rTemp 
            END DO  
         iaGasID(ii)=iGasID 
         END DO 
       END IF 
c -------- end optional due to iLorS ---------------- 
 
c read in Number of Paths to be output each time 
       read(iIOUN) iNumPathsOut 
       IF (iLorS .LT.  0) THEN 
         iNumPaths=iNumPathsOut 
         END IF 
 
c       print *,'output the following paths',iNumPathsOut 
       IF ((iNumPathsOut .LT. 0).OR.(iNumPathsOut .GT. iNumPaths)) THEN 
          print *,'Mistake in number of paths to be output!!!' 
          STOP 
          END IF 
      IF (iNumPathsOut .GT.  0) THEN 
c read in paths 
        read(iIOUN) (iaOutPaths(ii), ii=1,iNumPathsOut) 
        END IF 
 
c       print *,'the total number of gas paths is ',iNumPaths 
c       print *,'which gives num of gases',iNumGases 
c       print *,'Read in the Path nums, gasIDs, amts, temps' 
 
      RETURN 
      END 

c************************************************************************
c now read in the mixed path info 
      SUBROUTINE readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS) 
 
      include 'convolve.param' 
 
c iIOUN          = file IOUNIT number 
c iNpmix         = number of mixed paths 
c iNumMixPathsOut= number of mixed paths to be output 
c iLorS          = long or short data header 
      INTEGER iIOUN,iNpmix,iNumMixPathsOut,iLorS 
 
c local variables 
c iaOutMixPaths  = integer array containing list of paths to be output 
c iMixFileLines  = number of lines containing character description of MP's 
c caaMixedPathInfo == character description of mixed path info 
c raMixTemp      = array containing the Mixed Path Temperatures 
      INTEGER iaOutMixPaths(kPathsOut),iMixFileLines 
      REAL raMixTemp(kMixFilRows) 
      CHARACTER*130 caaMixedPathInfo(kMixFilRows) 
      INTEGER ii,kk 
      CHARACTER*130 caMixedPathInfo 
 
c read in number of mixed paths in *MIXFIL 
      read(iIOUN) iNpmix 
 
c -------- optional due to iLorS ---------------- 
      IF (iLorS .GT. 0) THEN 
c read in #of lines that have info in *MIXFIL 
        read(iIOUN) iMixFileLines 
 
        IF  (iNpmix .GT. 0) THEN 
c read in the mixing table 
          DO ii=1,iMixFileLines 
            read(iIOUN) caMixedPathInfo 
c           caaMixedPathInfo(ii,1:130)=caMixedPathInfo(1:130) 
            END DO 
 
c read in mixed path temperatures 
          read(iIOUN) (raMixTemp(ii), ii=1,iNpmix) 
          END IF 
        END IF 
c -------- end optional due to iLorS ---------------- 
 
c read in Number of mixed paths to be output each time 
      IF (iNpmix .GT. 0) THEN 
        read(iIOUN) iNumMixPathsOut 
        IF ((iNumMixPathsOut.LT.0).OR.(iNumMixPathsOut.GT.iNpmix)) THEN 
            print *,'Mistake in number of mixed paths to be output!!!' 
            STOP 
            END IF 
 
        IF (iNumMixPathsOut .GT. 0) THEN 
c read in paths 
          read(iIOUN) (iaOutMixPaths(ii),ii=1,iNumMixPathsOut) 
          END IF 
 
        END IF 
 
c      print *,'Number of Mixed Paths = ',iNpmix 
c      print *,'Number of lines in MIXFIL to read = ',iMixFileLines 
 
      RETURN 
      END 
c************************************************************************
c this subroutine reads in the atmosphere information 
      SUBROUTINE readatmospheres(iIOUN,iNpmix,iNatm, 
     $           iaNumLayersInAtm,iaNumLayersOut) 
 
      include 'convolve.param' 
 
c iIOUN         = file unit number 
c iNpmix        = number of mixed paths (from readmixedpaths) 
c iNatm         = number of atmospheres 
c iaNumLayersInAtm = integer array with number of layers in each atmosphere 
c iaNumLayersOut= integer array with number of layers out per atmosphere 
      INTEGER iaNumLayersOut(kMaxAtm),iaNumLayersInAtm(kMaxAtm) 
      INTEGER iIOUN,iNpmix,iNatm 
 
c local variables 
c raaOutPress   = matrix containing pressures where radiances are output 
c iaAtmNum      = integer array with atmospher number (1..iNatm) 
c iaaRadLayer   = for each atmosphere, this lists which layers build up the atm 
c iaaOutRadPaths= for each atmosphere, list which layers are to be output 
c raTbdy,raTinit,raSat = for each atmosphere, these arrays contain atm temp, 
c                          surface temp and satellite viewing angle 
c raaaEms      =emissivity 
c rakSolarAngle = solar angles for the atmospheres 
c rakThermalAngle=thermal diffusive angle 
c rakSolarRefl   =solar reflectance 
c iakthermal,iaksolar = turn on/off solar and thermal 
c iakthermaljacob=turn thermal jacobians on/off       
c raheight = height of instrument above earth 
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
 
c      print *,' reading in the atmospheres info .... ' 
 
c read in number of atmospheres in *RADFIL 
      read(iIOUN) iNatm 
c read in kEmsRegions 
      read(iIOUN) ii
      IF (ii .GT. kEmsRegions) THEN 
        print *,'please recompile code ... kmsRegions incorrect' 
        STOP 
        END IF 
 
c      print *,'There are atmospheres in RADFIL ',iNatm 
 
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
c         print *,'reading in info for atm num ',ii 
 
c read in atmosphere number, # of mixed paths in atmosphere 
         read(iIOUN) iAtm,iNumPathsAtmos 
         iaNumLayersInAtm(ii)=iNumPathsAtmos 
c         print *,'Num of Paths in Atmosphere = ',iNumPathsAtmos 

c         IF((iNumPathsAtmos.LT.0).OR.(iNumPathsAtmos.GT.iNpmix)) THEN 
         IF((iNumPathsAtmos.LT.0).OR.(iNumPathsAtmos.GT.kProfLayer)) THEN 
           print *,'Mistake in num of mixed paths used in atmosphere !' 
           STOP 
           END IF 
        IF (iNumPathsAtmos .GT.  0) THEN 
c          print *,'Reading the Mixed Path numbers ...' 
c read in paths making up atmosphere 
          read(iIOUN) (iaaRadLayer(ii,jj),jj=1,iNumPathsAtmos)  
          read(iIOUN) raTbdy(ii),raTinit(ii),raSat(ii),raHeight(ii) 
          read(iIOUN) iakSolar(iI),rakSolarAngle(iI), 
     $       rakSolarRefl(iI),iakThermal(iI), 
     $       rakThermalAngle(iI),iakThermalJacob(iI) 
          read(iIOUN) iaEms(ii) 
          DO jj=1,iaEms(ii) 
            read(iIOUN) raaaEms(ii,jj,1),raaaEms(ii,jj,2) 
            END DO 
 
c now read in layers to be output 
          read(iIOUN) iNumLayersOut 
          iaNumLayersOut(ii)=iNumLayersOut 
c          print *,'This atm has layers to be output',iNumLayersOut 
 
      IF((iNumLayersOut.LT.0).OR.(iNumLayersOut.GT.iNumPathsAtmos))THEN  
        print *,'Mistake in number of layers to be output for Atm#!',ii 
        STOP 
        END IF 
          IF (iNumLayersOut .GT. 0) THEN 
c read in mixed paths to be output 
            read(iIOUN) (iaaOutRadPaths(ii,jj),jj=1,iNumLayersOut) 
            read(iIOUN) (raaOutPress(ii,jj),jj=1,iNumLayersOut) 
c            print *,(iaaOutRadPaths(ii,jj),jj=1,iNumLayersOut) 
            END IF 
          END IF 
 
          END DO 
        END IF 
 
      RETURN 
      END 

c************************************************************************
c this prints out a summary of the saved paths/MP/layers 
      SUBROUTINE printsummary(iaOutPaths,iaOutMixPaths,iNatm, 
     $                     iaNumLayersOut,iaaOutRadPaths,ii1,ii2) 
 
      include 'convolve.param' 
c iaOutPaths     = integer array with list of gas   paths to be output 
c iaOutMixPaths  = integer array with list of mixed paths to be output 
c iNatm          = number of atmospheres 
c iaNumLayersOut = for each atmosphere, how many layers to be output 
c iaaOutRadPaths = for each atmosphere, list of layers to be output 
c ii1            =iNumPathsOut 
c ii2            =iNumMixPathsOut 
      INTEGER iaOutPaths(kPathsOut),iaOutMixPaths(kPathsOut),iNatm 
      INTEGER iaNumLayersOut(kMaxAtm) 
      INTEGER iaaOutRadPaths(kMaxAtm,kPathsOut) 
      INTEGER ii1,ii2 
 
c local variables 
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

c************************************************************************
c this subroutine reads in the data
      SUBROUTINE readdata(iSetLow,iMaxSet,iBlSz,
     $           iIOUN,iStore,iNumPathsOut,iNumMixPathsOut,
     $           iaNumLayersOut,iNatm,iLoop, raEntire,raFreq,iFreqPts,
     $           iTotal,iOutNumber,iaNumOut,iOutTypeLoop,iLpLp)

      include 'convolve.param'

c iSetLow    = lowest  k-comp block number (1-89  corresponds to 605-2805)
c iMaxSet    = biggest k-comp block number (1-89  corresponds to 605-2805)
c iBlSz      = number of data points per chunk
c iIOUN      = file IOUNIT number
c iStore     = total number of paths/MP/layers output per each 245 cm-1 block
c iNumPathsOut    = number of paths output each time
c iNumMixPathsOut = number of mixed paths output each time
c iaNumLayersOut  = number of layers output each time, for each atmosphere
c iNatm      = number of atmospheres
c iLoop      = which of the iStore paths/MP/layers to be entirely saved
c raEntire   = all the data for one of the chosen paths/MP/layers is saved here
c raFreq     = frequency array for the saved data
c iFreqPts   = total number of points saved
      INTEGER iMaxSet,iSetLow,iIOUN,iStore
      INTEGER iSetHigh,iNatm,iBlSz,iLoop
      INTEGER iNumPathsOut,iNumMixPathsOut,iaNumLayersOut(kMaxAtm)
      REAL raEntire(kMaxEntire),raFreq(kMaxEntire)
      INTEGER  iTotal,iOutNumber,iaNumOut(kMaxPrint),iOutTypeLoop,iLpLp      

c local variables
      REAL raEntire1(kMaxEntire),raFreq1(kMaxEntire),ra(kMaxPts)
      INTEGER ii,jj,kk,indd(kMaxEntire),RRR,LLL,iFreqPts,iNumPtsOut
      INTEGER iFrLow,iFrHigh,iNbl,iEnd,iSSeeTT,iSMT,ikMaxPts
      REAL rL,rH,rDelta
      INTEGER iWhich25block,kkmax,kkremain

      iSetHigh=iMaxSet

      iWhich25block=0

 111  CONTINUE
      iWhich25block=iWhich25block+1
c first go to the offset we need to start with
      ii=0
      DO jj=1,iOutTypeLoop
c first read the header
        read(iIOUN) iNbl,iSMT,iSSeeTT
        read(iIOUN) ikMaxPts,rL,rH,rDelta

        IF (jj .EQ. iOutTypeLoop) THEN
          kkmax=iLpLp-1
          kkremain=iaNumOut(jj)-1
        ELSE
          kkmax=iaNumOut(jj)
          kkremain=0
          END IF

        DO kk=1,kkmax
c then read the data
          ii=ii+1
          CALL  readblockDATAONLY(ra,iIOUN,iBlSz)
          END DO
        END DO

      IF (ii .NE. (iLoop-1)) THEN
        print *,'Something wrong in readdata',ii,iLoop
        STOP
        END IF
   
c NOW proceed to read and save the current data into relevant part of raEntire
      CALL read_data(iIOUN,iSetLow,iSetHigh,iStore,iBlSz,
     $       raEntire,raFreq,rL,rH,rDelta,iSetLow+(iWhich25block-1))

c now go to end of current 25 cm-1 set of data
      DO kk=kkmax+1,kkremain
        CALL  readblockDATAONLY(ra,iIOUN,iBlSz)
        END DO

      DO jj=iOutTypeLoop+1,iOutNumber
c first read the header
        read(iIOUN) iNbl,iSMT,iSSeeTT
        read(iIOUN) ikMaxPts,rL,rH,rDelta

        kkmax=iaNumOut(jj)
        kkremain=0
        DO kk=1,kkmax
          CALL  readblockDATAONLY(ra,iIOUN,iBlSz)
          END DO
        END DO
               
c now loop over remaining 25 cm-1 chunks
      IF (iWhich25block  .LT. (iSetHigh-iSetLow+1)) THEN
        GOTO 111
        END IF

      iFreqPts=kMaxPts*(iSetHigh-iSetLow+1)
      close(iIOUN)

      RETURN
      END

c************************************************************************
c this reads in a "block" of data
      SUBROUTINE readblockDATAONLY(ra,iIOUN,iBlSz)

      include 'convolve.param'

c iBlSz          = block size == number of points to be read in (400 or 10000)
c iIOUN          = file IOUNIT number
c ra             = real array containing the data that is read in
      INTEGER iIOUN,iBlSz,ii
      REAL ra(kMaxPts)

c now read the data
      read(iIOUN) (ra(ii),ii=1,iBlSz)

      RETURN
      END

c************************************************************************
c this reads in the entire frequency range for ONE of the paths/MPs/layers
      SUBROUTINE read_data(iIOUN,iSetLow,iSetHigh,iStore,iBlSz,
     $          raEntire,raFreq,rL,rH,rDelta,ii)

      include 'convolve.param'

c ii         = which of the 25 cm-1 chunks we are reading in currently
c iSetLow    = lowest  k-comp block number (1-89  corresponds to 605-2805)
c iSetHigh   = biggest k-comp block number (1-89  corresponds to 605-2805)
c iStore     = how many paths/MP/layers have been saved for each 25 cm-1 block
c iBlSz      = number of data points per chunk
c iIOUN      = file IOUNIT number
c raEntire   = save entire data for iLoop path/MP/layer
c raFreq     = frequency array for the saved data
c rL,rH          = freq start/stop points, corresponding to above two indices
      INTEGER iSetLow,iSetHigh,iBlSz,iIOUN,iStore,ii
      INTEGER iNbl,iSSeeTT,iSMT,ikMaxPts
      REAL raEntire(kMaxEntire),raFreq(kMaxEntire),rL,rH,rDelta

c local variables
      REAL ra(kMaxPts)
      INTEGER jj,kk,mm,ind(kMaxPts),iFrLow,iFrHigh     

      DO mm=1,kMaxPts
        ind(mm)=kMaxPts*(ii-iSetLow)+mm
        END DO
      CALL readblockDATAONLY(ra,iIOUN,iBlSz)
      DO mm=1,kMaxPts
        raEntire(ind(mm))=ra(mm)
        raFreq(ind(mm))=rL+rDelta*(mm-1)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does the command line stuff
      SUBROUTINE DoCommandLine(caInName,caOutName,iTextOrBinary,iWhich)

      CHARACTER*80 caInName   !!!Input binary file, written out by kcarta
      CHARACTER*80 caOutName  !!!Output file, written out by this code
                              !!!if only 2 args, this is binary file
                              !!!if      3 args, this is text file
      INTEGER iTextOrBinary   !!! -1 if text   file (3 arguments)
                              !!! +1 if binary file (2 arguments)
      INTEGER iWhich          !!! if text file, which one to output

      CHARACTER*80 caWhich
      INTEGER iDummy,iError

c      print *,'Enter INPUT binary file name '
c      read 1000,caInName
c      print *,'Do you want TEXT (-1) or BINARY (+1) output?'
c      read *,iTextOrBinary
c      IF (iTextOrBinary .EQ. -1) THEN
c        print *,'input file has ',iTotalStuff,' outputs'
c        print *,'which one do you want output???'
c        read *,iWhich
c        END IF
c      print *,'Enter OUTPUT file name '
c      read 1000,caOutName


c this is the number of args
      INTEGER iargc 

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
c************************************************************************
      SUBROUTINE GetNumber(caName,iW)

      CHARACTER*80 caName
      INTEGER iW

      CHARACTER*80 caTempName
      INTEGER iR,iL,iInt,iFactor,iI
      CHARACTER c

c find the "right" length of the input root name
      iR=len(caName)
 11   continue
      IF (caName(iR:iR) .eq. ' ') THEN
        iR=iR-1
        GO TO 11      
        END IF

c find the "left" length of the input root name
      iL=1
 12   continue
      IF (caName(iL:iL) .eq. ' ') THEN
        iL=iL+1
        GO TO 12      
        END IF

c thus the entered word exists between iL:iR
      IF (iL .GT. iR) THEN
        print *, 'Fatal Error! Invalid (blank) string in file !'
        STOP
        END IF

c now rearrange the string, so that it is right padded with blanks
c this is basically equivalent to  ADJUSTR(caName)
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
        c = caName(iInt:iInt)        
        READ(c,13) iI
        iW = iW + iFactor*iI
        iFactor = iFactor/10
        END DO

 13   FORMAT(I1)

      RETURN
      END


c************************************************************************
