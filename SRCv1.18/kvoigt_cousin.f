c Copyright 2002  
c University of Maryland Baltimore County  
c All Rights Reserved 

c this has the voigt and cousin stuff
c************************************************************************
c this simple routine does a simple boxcar ala the Matlab/GENLN2 code ...
c ie simple boxcar integration of daTempClose from i = 1..iFineMeshBoxPts
c but puts things smartly into daK from i1 to i2 (ie keeps track of the 
c moving mesh)
      SUBROUTINE DoBoxCar(daTempClose,iFineMeshBoxPts,i1,i2,daK)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters
      DOUBLE PRECISION daTempClose(kMaxPtsBox)    !!!input
      INTEGER iBoxCarUSe                          !!!width of boxcar
      INTEGER iFineMeshBoxPts                     !!!number of relevant pts
      INTEGER i1,i2                               !!!store results within here
c output parameter
      DOUBLE PRECISION daK(kMaxPts)

c local variables
      INTEGER iI,iJ,iIm2,iIm1,iIp1,iIp2,iSpline
      DOUBLE PRECISION daSpline(kMaxPts)

      IF (kBoxCarUse .EQ. 5) THEN
        !do a 5 point average
        iI = -2
        iSpline = 0
        DO iJ = 1,iFineMeshBoxPts-2
          iI   = iI + 5 
          iIm2 = iI - 2 
          iIm1 = iI - 1 
          iIp1 = iI + 1 
          iIp2 = iI + 2 
          iSpline = iSpline + 1
          daSpline(iSpline) = 
     $           (daTempClose(iIm2) + daTempClose(iIm1) + daTempClose(iI) +
     $            daTempClose(iIp1) + daTempClose(iIp2))/5.0d0
        END DO
      ELSEIF (kBoxCarUse .EQ. 3) THEN
        !do a 3 point average
        iI = -1
        iSpline = 0
        DO iJ = 1,iFineMeshBoxPts-1
          iI   = iI + 2 
          iIm1 = iI - 1 
          iIp1 = iI + 1 
          iSpline = iSpline + 1
          daSpline(iSpline) = 
     $        (daTempClose(iIm1) + daTempClose(iI) + daTempClose(iIp1))/3.0d0
        END DO
      ELSEIF (kBoxCarUse .EQ. 1) THEN
        !do a 1 point average .... do nothing!!!!
        DO iJ = 1,iFineMeshBoxPts
          daSpline(iJ) = daTempClose(iJ)
        END DO
      END IF

      DO iJ = i1,i2
        daK(iJ) = max(daSpline(iJ-i1+1),0.0d0)
      END DO

      RETURN
      END

c************************************************************************
c do we stick to linemix, or switch to Cousin?
c for testing purposes with Genln2, we have parameter iTestGenln2 = +1
c where we swtich EVERYTHING to Cousin; else the default is set at 
c iTestGenln2 = -1 so we use linemixing for the strongest bands
      SUBROUTINE CousinVsMix(iaLineMix,iLineMixBand,iJL,iJU,iISO,daFreq,
     $                       daLineShift,daStren296,dVibCenter,iNum)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters
      INTEGER iLineMixBand                  !do we use linemixing coeffs
      DOUBLE PRECISION daFreq(kMaxPts)      !the wavevector
      DOUBLE PRECISION daLineShift(kHITRAN) !line centers
      DOUBLE PRECISION daStren296(kHITRAN)  !line strengths at 296 K
      INTEGER iNum                          !number of lines in band
      DOUBLE PRECISION dVibCenter           !vib center of band
      INTEGER iISO                          !isotope
      INTEGER iJL                           !lower quantum number      
      INTEGER iJU                           !upper quantum number
c output parameter
      INTEGER iaLineMix(kHITRAN)            !individual lines linemix coeffs
                                            !set to +2 for linemix, +1 for cous
c local vars
      INTEGER iOut,iI,iStrong,iaWeak(kHITRAN),iSwitched
      INTEGER iTestGenln2,iDefault
      DOUBLE PRECISION d1,d2,dLineShift
      DOUBLE PRECISION dfL_Pbranch,dfR_Pbranch,dfL_Rbranch,dfR_Rbranch
      DOUBLE PRECISION dSmin,dSmax,dDelta,dThreshold

      iSwitched   = -1

      iDefault    = -1       !!default; use line mix for strong bands
      iTestGenln2 = +1       !!test with GENLN2, so use Cousin everywhere
      iTestGenln2 = -1       !!default; use line mix for strong bands

      IF (iTestGenln2 .NE. iDefault) THEN
        print *,'iTestGenln2,iDefault = ',iTestGenln2,iDefault
      END IF
      
      dsMin = +1.0d36
      dsMax = -1.0d36
      dThreshold = 1.0d0/2.0d0    !!!make sure this is  0 <= threshold <= 1
      dThreshold = 0.0d0

      DO iI = 1,iNum
        iaWeak(iI) = +1      !!!!assume line is "weak"
        IF (daStren296(iI) .LE. dSmin) dSmin = daStren296(iI)
        IF (daStren296(iI) .GE. dSmax) dSmax = daStren296(iI)
      END DO
      dDelta = log10(dSmax)-log10(dSmin)
      dThreshold = log10(dSmin) + dThreshold*dDelta
 
      DO iI = 1,iNum
        IF (log10(daStren296(iI)) .GT. dThreshold) iaWeak(iI) = -1
      END DO

      d1 = daFreq(1)
      d2 = daFreq(kMaxPts)

      !!!stick to linemix close to/inside the 2350,2351 and 2310,2320 bands
      !!!else switch to Cousin

      DO iI = 1,iNum

        iOut = 2   
        iStrong = -1

        dLineShift = daLineShift(iI)

        IF ((iLineMixBand .EQ. 2) .AND. (iJU .EQ. 9) .AND. (iISO .EQ. 1)) THEN
          !!strongest sigsig band  ... do linemix here, epecially for R branch
          !!P branch from 2150 to 2387
          dfL_Pbranch = 2155.0d0
          dfR_Pbranch = 2380.0d0
          !!R branch from 2328 to 2436 cm-1
          dfL_Rbranch = 2330.0d0
          dfR_Rbranch = 2430.0d0
          iStrong = +1
        ELSEIF ((iLineMixBand .EQ. 2).AND.(iJU .EQ. 9).AND.(iISO .EQ. 2)) THEN
          !!next strongest sigsig band (isotope)
          !!P branch from 2167.036826 to 2296.920256
          dfL_Pbranch = 2180.0d0
          dfR_Pbranch = 2280.0d0
          !!R branch from 2269.261676 to 2345.936568
          dfL_Rbranch = 2280.0d0
          dfR_Rbranch = 2355.0d0
          iStrong = +1
          IF (iaWeak(iI) .GT. 0) iStrong = -1   !!!weak line in band
        ELSEIF ((iLineMixBand .EQ.2).AND.(iJU .EQ. 24).AND.(iISO .EQ. 1)) THEN
          !!strongest deltdelt band
          !!P branch from 2227.255209 to 2321.797519
          dfL_Pbranch = 2230.0d0
          dfR_Pbranch = 2330.0d0
          !!R branch from 2326.429307 to 2369.735017
          dfL_Rbranch = 2330.0d0
          dfR_Rbranch = 2380.0d0
          iStrong = +1
          IF (iaWeak(iI) .GT. 0) iStrong = -1   !!!weak line in band
        ELSEIF ((iLineMixBand.EQ.2).AND.(iJU .EQ. 16).AND.(iISO .EQ. 1)) THEN
          !!strongest pipi band
          !!P branch from 2227.569529 2335.086137
          dfL_Pbranch = 2230.0d0
          dfR_Pbranch = 2330.0d0
          !!R branch from 2338.151553 2381.80012
          dfL_Rbranch = 2330.0d0
          dfR_Rbranch = 2380.0d0
          iStrong = +1
          IF (iaWeak(iI) .GT. 0) iStrong = -1   !!!weak line in band
        END IF

        IF ((iStrong .EQ. 1) .AND. (iTestGenln2 .EQ. 1)) THEN
          iSwitched = +1
          iStrong = -1
        END IF

        IF ((iLineMixBand .EQ. 2) .AND. (iStrong .GT. 0)) THEN
          !!current 10000 point chunk is way outside the band
          IF ((d1 .LE. dfL_Pbranch).AND.(dLineShift .LE. dVibCenter)) THEN
            iOut = 1             !!!!switch to Cousin 
          ELSEIF ((d2 .GE. dfR_Pbranch) .AND. (dLineShift .LE. dVibCenter))THEN
            iOut = 1             !!!!switch to Cousin
          ELSEIF ((d1 .LE. dfL_Rbranch) .AND. (dLineShift .GE. dVibCenter))THEN
            iOut = 1             !!!!switch to Cousin 
          ELSEIF ((d2 .GE. dfR_Rbranch) .AND. (dLineShift .GE. dVibCenter))THEN
            iOut = 1             !!!!switch to Cousin 
          END IF
          
          !! for some odd reason, do all the P branch lines for 2350 sigsig, 
          !! using the Cousin lineshape
          !! this is because at large pressure I notice FULLMIX << FIRSTORDER
          IF ((dLineShift .LE. dVibCenter)) THEN
            iOut = 1             !!!!switch to Cousin 
          END IF
            
        ELSEIF ((iLineMixBand .EQ. 2) .AND. (iStrong .LT. 0)) THEN
          !! current line is one of the weaker lines in the band OR
          !! weaker bands in NLTE ... just use Cousin
          iOut = 1
        END IF
          
        iaLineMix(iI) = iOut
      END DO

      IF (iSwitched .EQ. +1) THEN
        write(kStdWarn,*) 'testing against GENL2; reset from linemix to cousin'
      END IF
       
      RETURN
      END

c************************************************************************
c this subroutine calls the voigt function at whatever resolution
      SUBROUTINE voigt_chi(daFreq,iFreqPts,dLineShift,
     $                 dLTE,dMass,dBroad,iLineMix,dP,dPP,
     $                 dJL,dJU,dJLowerQuantumRot,iISO,dVibCenter,
     $                 xBirn,chiBirn,iNptsBirn,daLineshape,daFudge,iDoVoigtChi,iTooFar)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters
      INTEGER iFreqPts   !number of relevant points in daFreq 
      INTEGER iLineMix   !do we use linemixing coeffs or cousin
      DOUBLE PRECISION daFreq(kMaxPts)                !the wavevector
      DOUBLE PRECISION dP,dPP                         !layer pressure
      DOUBLE PRECISION dLineShift,dMass,dBroad        !line parameters
      DOUBLE PRECISION dLTE,dJL,dJU,dJLowerQuantumRot
      DOUBLE PRECISION dVibCenter                     !vib center of band
      INTEGER iISO
      DOUBLE PRECISION daFudge(kMaxPtsBox)  !fudge so that we mimic UMBC-LBL
      INTEGER iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
c this is for birnbaum, if needed
      INTEGER iNptsBirn
      DOUBLE PRECISION xBirn(kMaxPts),chiBirn(kMaxPts)
c output parameters
      DOUBLE PRECISION daLineshape(kMaxPtsBox)    ! lineshape, at high res
      INTEGER iToofar                             ! is the input linecenter too far from linecenters
                                                  ! found in linemix file (-1 No, +x yes == in(dJLowerQuantumRot))
                                                  ! if larger than 100, high J so weak line, so dont worry
c local parameters
      DOUBLE PRECISION daImag(kMaxPtsBox)                   !for linemixing calcs
      DOUBLE PRECISION daHighFreqWavenumbers(kMaxPtsBox)    !for computing lineshape, at high res      
      DOUBLE PRECISION daChi(kMaxPtsBox),dT
      DOUBLE PRECISION df,f0,tau2_birn,tau2
      DOUBLE PRECISION daYmix(kHITRAN),daYmixALL(kHITRAN),dY
      INTEGER iFr,iN,iLine,iNum

c      iDoVoigtChi = -1  !! do not multiply by chi function in voigt_chi; all 
c                        !!   necesary things will be done by calling 
c                        !!   Co2_4um_fudge_nlte_fast
c      iDoVoigtChi = +1  !! do multiply by chi function in voigt_chi

      iTooFar = -1
      iN = iFreqPts

      !!! this is NEW Jan 2017
      !!! this is NEW Jan 2017 !!!!
      iN = iFreqPts*5
      DO iFr = 1,kMaxPtsBox
        daHighFreqWavenumbers(iFr) = dble(daFreq(1)) + (iFr-1)*(daFreq(2)-daFreq(1))/(kBoxCarUse * 1.0d0)
      END DO
      !!! compute the voigt lineshape at high resolution, used to be daFreq instead of daHighFreqWavenumbers before Jan 2017
      !!! CALL DoVoigt(daLineshape,daImag,daFreq,dLineShift,dLTE,dMass,dBroad,iN,dP,dPP)
      write(kStdErr,*)  ' >>>>>>>>>>> check this voigt stuff before making NLTE production run!!!! <<<<<<<<<<<<<'
      write(kStdWarn,*) ' >>>>>>>>>>> check this voigt stuff before making NLTE production run!!!! <<<<<<<<<<<<<'            
      CALL DoVoigt(daLineshape,daImag,daHighFreqWavenumbers,dLineShift,dLTE,dMass,dBroad,iN,dP,dPP)
      !!! this is NEW Jan 2017
      !!! this is NEW Jan 2017

      IF (iLineMix .EQ. 1) THEN    !this does cousin
        dT = dLTE
        CALL Cousin(daChi,daHighFreqWavenumbers,dLineShift,dBroad,dT,dP,dPP,1,iN)
        DO iFr = 1,iN
          daLineshape(iFr) = daLineshape(iFr) * daChi(iFr)
        END DO
      END IF

      IF (iLineMix .EQ. 2) THEN   !this does linemix
        dT = dLTE
c compare linemix vs linemixALL	
c        CALL linemix(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)
c        CALL linemixALL(dLineShift,daYmixALL,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)
c	do iN = 1,120
c          print *,iN,dLineShift,daYmix(iN),daYmixALL(iN),daYmixALL(iN)/daYmix(iN)
c	end do
c	call dostop
        CALL linemix(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)

ccc        tau2 = tau2_birn(dP,dPP)
ccc        CALL birnbaum(daChi,daFreq,dLineShift,dBroad,dT,tau2,iN)
c	print *,'doing birnbaum_interp ',dJL,dJU,iISO,daFreq(1)	
        CALL birnbaum_interp(daChi,daHighFreqWavenumbers,dLineShift,iN,chiBirn,xBirn,iNptsBirn)
	
        dY = daYmix(iLine)*dP
        DO iFr = 1,iN
          daLineshape(iFr) = (daLineshape(iFr)+daImag(iFr)*dY)*daChi(iFr)
        END DO
        IF ((dJU .EQ. 9) .AND. (iISO .EQ. 1) .AND. (iDoVoigtChi .GT. 0)) THEN
          !!oh oh check to see if we need to adjust the linemix to bring it
          !!approximnately equal to UMBC-LBL
          IF ((daFreq(1) .LE. 2505.0d0) .AND. (daFreq(iN) .GE. 2355.0d0)) THEN
            DO iFr = 1,iN
              daLineshape(iFr) = daLineshape(iFr)*daFudge(iFr)
            END DO
          END IF
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this is the fudger to bring strongest R branch linemix in line with UMBC-LBL
c need to interpolate the "chi" functions onto dafreq
      SUBROUTINE co2_4um_nlte_fudge(daFreq,daFudge,iN)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters
      INTEGER iN             !number of points
      DOUBLE PRECISION daFreq(kMaxPtsBox)      !the wavevector
c input/output parameters
      DOUBLE PRECISION daFudge(kMaxPtsBox)     !fudge factor interped to daFreq

c local variables
      INTEGER iFileStartFr  !!its ok to keep this as integer
      INTEGER iIOUN,iERR,iChi,iFr,iChunk
      CHARACTER*120 fname
      DOUBLE PRECISION daF(kMaxPts),daChi(kMaxPts)

      iChi = -1
      iFileStartFr = int(daFreq(1)+0.25)   !!! the +0.25 ensures it is in the
                                           !!! correct chunk

      IF ((iFileStartFR .GE. 2355) .AND. (iFileSTartFr .LE. 2505)) THEN
        IF ((iFileStartFR .GE. 2355) .AND. (iFileSTartFr .LT. 2380)) THEN 
          iChunk = 2355
          FNAME = 'nonlte_co2_fudge_2355.txt'
          iChi = +1
        ELSEIF ((iFileStartFR .GE. 2380) .AND. (iFileSTartFr .LT. 2405)) THEN 
          iChunk = 2380
          FNAME = 'nonlte_co2_fudge_2380.txt'
          iChi = +1
        ELSEIF ((iFileStartFR .GE. 2405) .AND. (iFileSTartFr .LT. 2430)) THEN 
          iChunk = 2405
          FNAME = 'nonlte_co2_fudge_2405.txt'
          iChi = +1
        ELSEIF ((iFileStartFR .GE. 2430) .AND. (iFileSTartFr .LT. 2455)) THEN 
          iChunk = 2430
          FNAME = 'nonlte_co2_fudge_2430.txt'
          iChi = +1
        ELSEIF ((iFileStartFR .GE. 2455) .AND. (iFileSTartFr .LT. 2480)) THEN 
          iChunk = 2455
          FNAME = 'nonlte_co2_fudge_2455.txt'
          iChi = +1
        ELSEIF ((iFileStartFR .GE. 2480) .AND. (iFileSTartFr .LT. 2505)) THEN 
          iChunk = 2480
          FNAME = 'nonlte_co2_fudge_2480.txt'
          iChi = +1
        END IF
      END IF

      IF (iChi .GT. 0) THEN
        CALL FindChiFileName(fname)
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          WRITE(kStdErr,*) 'In subroutine co2_4um_nlte_fudge'
          WRITE(kStdErr,1010) IERR, FNAME
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
        ENDIF
        kTempUnitOpen = 1
        READ(iIOUN,*) (daF(iFr),daChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen = -1

        CALL dspl(daF,daChi,kMaxPts,daFreq,daFudge,iN)
      END IF

      RETURN
      END

c************************************************************************
c this is the subroutine that actually calls the voigt fcn
c kinda dumb to call this instead of calling vhh2RI direct, but in a past time
c this subroutine used to call cousin or birnbaum/linem code as well
      SUBROUTINE DoVoigt(daReal1,daImag1,daFreq1,dLineShift,dLTE,
     $                   dMass,dBroad,iN,dP,dPP)

      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 

c input parameters
      INTEGER iN                                       !how many points
      DOUBLE PRECISION daFreq1(kMaxPtsBox)             !the wavevector
      DOUBLE PRECISION dLineShift,dMass,dBroad         !line parameters
      DOUBLE PRECISION dP,dPP                          !layer pressure
      DOUBLE PRECISION dLTE
c output parameters
      DOUBLE PRECISION daReal1(kMaxPtsBox)             !lineshape (real part)
      DOUBLE PRECISION daImag1(kMaxPtsBox)             !lineshape (imag part)

c local variables
      INTEGER iFr

      CALL vhh2RI(daReal1,daImag1,daFreq1,dLineShift,dLTE,dMass,dBroad,iN)

cc for testing : turn off lineshape +/- 25 cm away!
cc this is bad!
cc      DO iFr = 1,iN
cc        IF (abs(dLineShift - daFreq1(iFr)) .GT. 25.0d0) THEN
cc          daReal1(iFr) = 0.0d0
cc          daImag(iFr) = 0.0d0
cc        END IF
cc      END DO
cc this is bad!
cc for testing : turn off lineshape +/- 25 cm away!

      RETURN
      END

c************************************************************************
c this is a quick approx to tanh
      DOUBLE PRECISION FUNCTION tanhsergio(x)

      DOUBLE PRECISION x
      DOUBLE PRECISION y

      y = 1.0D0 !default value

      IF (abs(x) .LT. 6.5d0) y = tanh(x)

      tanhsergio = y

      RETURN
      END

c************************************************************************
c this is the /home/sergio/SPECTRA/FORTRANLINUX/vhh2RI.f code  
c (used by run7.m) 
c this is a Voigt-VanHuber function
      subroutine vhh2RI(wr,wi,v,v0,Temp,m,brd,n_in)
c this is the CURRENT GENLN2 routine 
c this is the same as voivec2.f in ../Genln2 

c subroutine voigt1(wr,wi,0,Temp,m,brd,n_in)
c gets the real and imag parts of the fcn
c v    = frequency array
c v0   = center freq
c T    = temperature
c m    = molecular mass (amu)
c brd  = broadening

      include '../INCLUDE/kcarta.param' 

      DOUBLE PRECISION wr(kMaxPtsBox),v(kMaxPtsBox),v0,Temp,m,brd
      DOUBLE PRECISION wi(kMaxPtsBox)
      DOUBLE PRECISION tanhsergio

      integer n_in

c PROGRAM        VOIVEC     SUBROUTINE          
c
c PURPOSE        COMPUTE COMPLEX PROBABILITY FUNCTION
c
c VERSION        3.X   D.P. EDWARDS   28/05/92
c
c DESCRIPTION    THIS ROUTINE CALCULATES THE COMPLEX PROBABILITY
c                FUNCTION USING A VECTORIZED VERSION OF THE
c                HUMLICEK JQSRT V27 437 1982 PAPER.
c                THE CALCULATION IS PERFORMED FOR THE ARRAY OF X,Y
c                PAIRS FOR A GIVEN LINE OVER THE FINE MESH POINTS
c                OF THE CURRENT WIDE MESH. 
c
c ARGUMENTS      NUM    I*4 I/P NUMBER OF FINE GRID INTERVALS 
c                NL     I*4 I/P FREQUENCY BDY TO START CALCULATION
c                NH     I*4 I/P FREQUENCY BDY TO STOP CALCULATION
c                NSTP   I*4 I/P FREQUENCY BDY STEP
c                 FF     DP  I/P FINE WAVENUMBER GRID [cm-1]
c                 WNUM   DP  I/P WAVENUMBER OF LINE CENTRE [cm-1] 
c                 REPWID R*4 I/P SQRT(ln2)/(DOPPLER WIDTH [cm-1])
c                 Y      R*4 I/P VOIGT Y PARAMETER OF LINE
c                 VT     COM O/P COMPLEX PROBABILITY FUNCTION

c local variables
       DOUBLE PRECISION factor,c2
       integer ii,iDoVanHuber
       complex*16 U,T,VTP(kMaxPtsBox),VTM(kMaxPtsBox)

c local variables

      DOUBLE PRECISION k,c_light,amu,mass,r2,alpha_doppler,g0
      DOUBLE PRECISION repwid,Y,X,S1V,S2V,fact,c2_0

      IF (n_in .gt. kMaxPtsBox) THEN
        write(kStdErr,*)'in vhh2RI.f, n_in .gt. kMaxPtsBox',n_in,kMaxPtsBox
        Call DoStop
      END IF

c  SORT THE (X,Y) PAIRS INTO THE 4 REGIONS OF THE HUMLIcEK 
c  EXPRESSIONS OF THE VOIGT LINE PROFILE.
c

c do the doppler widths first 
      k       = kBoltzmann
      c_light = 2.99792458d8         !ms-1
      amu     = 1.6605402d-27            !nucleon mass/kg
      mass    = m                       !change to kg

c alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
c r2=2*log(2)*k/amu
      r2=11526.218d0
      alpha_doppler=v0/c_light*sqrt(r2*Temp/mass)
      repwid=0.8325546d0/alpha_doppler

c do the g0 factor 
c  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895 
c g0=sqrt(log(2)/pi)/alpha_doppler
      g0= 0.83255461 * 0.5641895 / alpha_doppler


c------------------------------  do the v-vi part ---------------------------
c define arrays for the new Voigt fcn 
      Y=brd/alpha_doppler*0.83255461

       do 20 II=1,n_in
         X =  (V(II) - V0)*REPWID
         S1V = abs(X) + Y
         S2V = (0.195*abs(X)) - 0.176
         T = dcmplx(Y,-X)
c
c  FOR REGION 1 OF HUMLIcEK
c
         if (S1V .ge. 15.0) then
           VTP(II) = T*0.5641896/(0.5+(T*T))
c
c  REGION 2 OF HUMLIcEK
c
         elseif (S1V .ge. 5.5) then
           U = T*T
           VTP(II) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
c
c  REGION 3 OF HUMLIcEK
c
         elseif (Y .ge. S2V) then
           VTP(II) =
     1     (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+
     2     T*.5642236))))/
     3     (16.4955+T*(38.82363+T*(39.27121+
     4     T*(21.69274+T*(6.699398+T)))))
c
c  REGION 4 OF HUMLIcEK
c
         else
           U = T*T
           VTP(II)=cdexp(U)-T*(36183.31-U*(3321.9905-
     1     U*(1540.787-U*(219.0313-U*
     2     (35.76683-U*(1.320522-U*.56419))))))/
     3     (32066.6-U*(24322.84-U*
     4     (9022.228-U*(2186.181-U*(364.2191-
     5     U*(61.57037-U*(1.841439-U)))))))
       endif
c
       U = dcmplx(1.0d0, -1.0d0)
c
   20  continue
c
c

c------------------------------  do the v+vi part ---------------------------
c define arrays for the new Voigt fcn 
      Y=brd/alpha_doppler*0.83255461

       do 30 II=1,n_in
         X =  (V(II) + V0)*REPWID
         S1V = abs(X) + Y
         S2V = (0.195*abs(X)) - 0.176
         T = dcmplx(Y,-X)
c
c  FOR REGION 1 OF HUMLIcEK
c
         if (S1V .ge. 15.0) then
           VTM(II) = T*0.5641896/(0.5+(T*T))
c
c  REGION 2 OF HUMLIcEK
c
         elseif (S1V .ge. 5.5) then
           U = T*T
           VTM(II) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
c
c  REGION 3 OF HUMLIcEK
c
         elseif (Y .ge. S2V) then
           VTM(II) =
     1     (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+
     2     T*.5642236))))/
     3     (16.4955+T*(38.82363+T*(39.27121+
     4     T*(21.69274+T*(6.699398+T)))))
c
c  REGION 4 OF HUMLIcEK
c
         else
           U = T*T
           VTM(II)=cdexp(U)-T*(36183.31-U*(3321.9905-
     1     U*(1540.787-U*(219.0313-U*
     2     (35.76683-U*(1.320522-U*.56419))))))/
     3     (32066.6-U*(24322.84-U*
     4     (9022.228-U*(2186.181-U*(364.2191-
     5     U*(61.57037-U*(1.841439-U)))))))
       endif
c
       U = dcmplx(1.0d0, -1.0d0)
c
   30  continue
c

c --------------- now do the VHH part --------------------------------
c from the voigt subroutine we need some adjustment factors
cC
cC  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
cC       REPWID = 0.83255461/DOPWID
c       H0 = REPWID*0.5641895*STRPAR   ---> STRPAR is multiplied outside
c       Y = WIDPAR*REPWID
cC
cC  CALCULATE THE COMPLEX PROBABILITY FUNCTION
cC
c       CALL VOIVEC(NUM,NL,NH,NSTP,FF,WNUM,REPWID,Y,VT)
cC 
cC  COMPUTE ABSORPTION DUE TO VOIGT LINE SHAPE
cC
c       DO 10 IP=NL,NH,NSTP
c         XABS(IP) = H0*REAL(VT(IP))
c 10        CONTINUE       
cC

      c2=1.4387863d0        !K/ cm-1  from Genln2 manual 
 
      c2_0=0.5d0*c2*v0/Temp 
      c2_0=v0*tanhsergio(c2_0) 
      fact=repwid*0.5641895d0/c2_0 
 
      c2=0.5*c2/Temp 

      iDoVanHuber = -1         !!!          do plain voigt 
      iDoVanHuber =  0         !!!          do Lorentz
      iDoVanHuber = +1         !!! default, do VanHuber

      iDoVanHuber = +1

      IF (iDoVanHuber .NE. +1) THEN
        print *,'doing plain voigt or lorentz, NOT vhh',iDoVanHuber
      END IF

      IF (iDoVanHuber .GT. 0) THEN 
        !!! do VanHuber
        DO  ii = 1,n_in 
          factor  = fact*v(ii)*tanhsergio(c2*v(ii)) 
          vtp(ii) = factor*(vtp(ii)+vtm(ii)) 
          wr(ii)  = dreal(vtp(ii)) 
          wi(ii)  = dimag(vtp(ii)) 
        END DO 

      ELSEIF (iDoVanHuber .LT. 0) THEN 
        !!!! just do voigt; bad at lower wavenumbers!
        DO  ii = 1,n_in 
          factor  = fact*v(ii)
          wr(ii) = factor*dreal(vtp(ii)) 
          wi(ii) = factor*dimag(vtp(ii)) 
        END DO

      ELSEIF (iDoVanHuber .EQ. 0) THEN 
        !!! just do lorentz; bad at lower wavenumbers, and high altitudes
        !!! 1/pi = 0.318...
        factor = 0.31830988d0
        DO ii=1,n_in 
          vtp(ii) = brd/(brd*brd + (v(ii)-v0)**2)
          wr(ii) = factor*dreal(vtp(ii)) 
          wi(ii) = 0.0d0
        END DO
      END IF

      return
      end

c************************************************************************
