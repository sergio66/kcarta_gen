c************************ READS IN JACOBIANS ****************************
c************************************************************************
c************* READ IN THE k-compressed files ***************************
c************************************************************************
c this reads in the main header information 
      SUBROUTINE readmainjacobheader(iIOUN,rFr1,rFr2,iSetLow,iSetHigh) 
 
c iIOUN      = file number 
c rFR1,rFr2  = frequency start/stop pts (in cm-1) 
      INTEGER iIOUN 
      REAL rFr1,rFr2  
 
c caComment  = comment assocaiated with data run 
c iSetLow    = integer index (between 1 and 89) of kcomp files start point 
c iSetHigh   = integer index (between 1 and 89) of kcomp files end point 
c iNumLayers = number of layers in each gas profile 
      INTEGER iNumLayers,iSetLow,iSetHigh 
      CHARACTER*80 caComment 
 
c read comment 
      read(iIOUN) caComment 
 
c read number of layers 
      read(iIOUN) iNumLayers 
 
c read start,stop frequency  iSetLo,iSetHi 
      read(iIOUN) rFr1,rFr2  
      read(iIOUN) iSetLow,iSetHigh 
 
      RETURN 
      END 
c************************************************************************ 
c this reads in the additional Jacobian info
      SUBROUTINE readjacobinfo(iIOUN,iNumgases,iNatm,iaNumLayers)

      include 'convolve.param'

c iIOUN       = file number
c iNumgases   = number of gases
c iNatm       = number of atmosphers
c iaNumLayers = number of layers per atmosphere
      INTEGER iIOUN,iNumGases,iNatm,iaNumLayers(kMaxAtm)

      INTEGER ii

      read(iIOUN) iNumGases
      read(iIOUN) iNatm
      read(iIOUN) (iaNumLayers(ii),ii=1,iNatm)
       
      RETURN
      END

c************************************************************************
c this function reads in the individual gas profiles from the jacobian file
      SUBROUTINE readgasjacob(iIOUN,iNumGasesDQ,iNumLayers)

      include 'convolve.param'

c iIOUN          = file ioUNIT number
c iNumLayers     = number of layers in each gas profile
c iaGasID        = array containing ID's of gases for Jacobian calcs
c iNumGasesDQ    = number of gases whose d/dq is output
c raaAmt,raaTemp = matrices containing the gas amts,temperatures
      INTEGER iIOUN,iNumLayers,iaGasId(kGasStore),iNumGasesDQ
      REAL raaTemp(kProfLayer,kGasStore),raaAmt(kProfLayer,kGasStore)

c local variables
      REAL rAmt,rTemp
      INTEGER iPath,iGasID,ii,jj

c now read the PathNum GasID Temperature Amount for each path in the profile
        DO ii=1,iNumGasesDQ
          DO jj=1,iNumLayers
            read(iIOUN) iGasID,rTemp,rAmt
            raaAmt(jj,ii)=rAmt
            raaTemp(jj,ii)=rTemp
            END DO 
         iaGasID(ii)=iGasID
         END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the data
      SUBROUTINE readdata(iSetLow,iMaxSet,iBlSz,
     $           iIOUN,iStore,iNatm,iLoop, raEntire,raFreq,iFreqPts,
     $           iTotal,iOutNumber,iaNumOut,iOutTypeLoop,iLpLp)

      include 'convolve.param'

c iSetLow    = lowest  k-comp block number (1-89  corresponds to 605-2805)
c iMaxSet    = biggest k-comp block number (1-89  corresponds to 605-2805)
c iBlSz      = number of data points per chunk
c iIOUN      = file IOUNIT number
c iStore     = total number of paths/MP/layers output per each 245 cm-1 block
c iNatm      = number of atmospheres
c iLoop      = which of the iStore paths/MP/layers to be entirely saved
c raEntire   = all the data for one of the chosen paths/MP/layers is saved here
c raFreq     = frequency array for the saved data
c iFreqPts   = total number of points saved
      INTEGER iMaxSet,iSetLow,iIOUN,iStore
      INTEGER iSetHigh,iNatm,iBlSz,iLoop
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

