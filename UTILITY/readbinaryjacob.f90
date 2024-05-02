!************************ READS IN JACOBIANS ****************************
!************************************************************************
!************* READ IN THE k-compressed files ***************************
!************************************************************************
! this reads in the main header information
      SUBROUTINE readmainjacobheader(iIOUN,rFr1,rFr2,iNumLayers)

! iIOUN      = file number
! rFR1,rFr2  = frequency srtart/stop pts (in cm-1)
      INTEGER iIOUN,iNumLayers
      REAL rFr1,rFr2 

! caComment  = comment assocaiated with data run
! iSetLow    = integer index (between 1 and 89) of kcomp files start point
! iSetHigh   = integer index (between 1 and 89) of kcomp files end point
! iNumLayers = number of layers in each gas profile
      INTEGER iSetLow,iSetHigh
      CHARACTER*80 caComment

! read comment
      read(iIOUN) caComment

! read number of layers
      read(iIOUN) iNumLayers

! read start,stop frequency  iSetLo,iSetHi
      read(iIOUN) rFr1,rFr2 
      read(iIOUN) iSetLow,iSetHigh

      RETURN
      END
!************************************************************************
! this reads in the additional Jacobian info
      SUBROUTINE readjacobinfo(iIOUN,iNumgases,iNatm,iaNumLayers)

      include 'convolve.param'

! iIOUN       = file number
! iNumgases   = number of gases
! iNatm       = number of atmosphers
! iaNumLayers = number of layers per atmosphere
      INTEGER iIOUN,iNumGases,iNatm,iaNumLayers(kMaxAtm)

      INTEGER ii

      read(iIOUN) iNumGases
      read(iIOUN) iNatm
      read(iIOUN) (iaNumLayers(ii),ii=1,iNatm)
       
      RETURN
      END

!************************************************************************
! this function reads in the individual gas profiles from the jacobian file
      SUBROUTINE readgasjacob(iIOUN,iNumGasesDQ,iNumLayers)

      include 'convolve.param'

! iIOUN          = file ioUNIT number
! iNumLayers     = number of layers in each gas profile
! iNumGasesDQ    = number of gases whose d/dq is output
      INTEGER iIOUN,iNumLayers,iNumGasesDQ

! local variables
! iaGasID        = array containing ID's of gases for Jacobian calcs
! raaAmt,raaTemp = matrices containing the gas amts,temperatures
      REAL rAmt,rTemp
      REAL raaTemp(kProfLayer,kGasStore),raaAmt(kProfLayer,kGasStore)
      INTEGER iPath,iGasID,ii,jj,iaGasId(kGasStore)

      
! now read the PathNum GasID Temperature Amount for each path in the profile
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

!************************************************************************
! this subroutine reads in the data
      SUBROUTINE readdata(iIOUN,iDataSize,raEntire,iDataPoints,iTotal)

      include 'convolve.param'

! iDataSize=number of bytes from start of data in this block, to same start 
!           of data in next kcomp block (so we can directly jump there)
      INTEGER iIOUN,iDataSize,iDataPoints,iTotal
      REAL raEntire(kMaxEntire)

      REAL raTemp(kMaxPts)
      INTEGER iI,iJ,iDummy
      INTEGER FTELL

      iDataPoints=0
      DO iI=1,iTotal      
        read(iIOUN) (raTemp(iJ),iJ=1,kMaxPts)

        DO iJ=1,kMaxPts
          iDataPoints=iDataPoints+1
          raEntire(iDataPoints)=raTemp(iJ)
          END DO

        IF (iI .NE. iTotal) THEN
          CALL fseek(iIOUN,idatasize,1)         !!!!if using g77
          END IF

        END DO

!      print *,(raTemp(iJ),iJ=9991,10000) 
      !this is the last piece of borsh and nonsense

      RETURN
      END
!************************************************************************
