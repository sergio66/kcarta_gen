! Copyright 2002
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kvoigt_cousin

USE basic_common
USE kcousin
USE klinemix
USE s_misc

IMPLICIT NONE

CONTAINS

! this has the voigt and cousin stuff
!************************************************************************
! this simple routine does a simple boxcar ala the Matlab/GENLN2 code ...
! ie simple boxcar integration of daTempClose from i = 1..iFineMeshBoxPts
! but puts things smartly into daK from i1 to i2 (ie keeps track of the
! moving mesh)
    SUBROUTINE DoBoxCar(daTempClose,iFineMeshBoxPts,i1,i2,daK)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    DOUBLE PRECISION :: daTempClose(kMaxPtsBox)    !!!input
    INTEGER :: iBoxCarUSe                          !!!width of boxcar
    INTEGER :: iFineMeshBoxPts                     !!!number of relevant pts
    INTEGER :: i1,i2                               !!!store results within here
! output parameter
    DOUBLE PRECISION :: daK(kMaxPts)

! local variables
    INTEGER :: iI,iJ,iIm2,iIm1,iIp1,iIp2,iSpline
    DOUBLE PRECISION :: daSpline(kMaxPts)

    IF (kBoxCarUse == 5) THEN
    ! o a 5 point average
        iI = -2
        iSpline = 0
        DO iJ = 1,iFineMeshBoxPts-2
            iI   = iI + 5
            iIm2 = iI - 2
            iIm1 = iI - 1
            iIp1 = iI + 1
            iIp2 = iI + 2
            iSpline = iSpline + 1
            daSpline(iSpline) = &
            (daTempClose(iIm2) + daTempClose(iIm1) + daTempClose(iI) + &
            daTempClose(iIp1) + daTempClose(iIp2))/5.0d0
        END DO
    ELSEIF (kBoxCarUse == 3) THEN
    ! o a 3 point average
        iI = -1
        iSpline = 0
        DO iJ = 1,iFineMeshBoxPts-1
            iI   = iI + 2
            iIm1 = iI - 1
            iIp1 = iI + 1
            iSpline = iSpline + 1
            daSpline(iSpline) = &
            (daTempClose(iIm1) + daTempClose(iI) + daTempClose(iIp1))/3.0d0
        END DO
    ELSEIF (kBoxCarUse == 1) THEN
    ! o a 1 point average .... do nothing!!!!
        DO iJ = 1,iFineMeshBoxPts
            daSpline(iJ) = daTempClose(iJ)
        END DO
    END IF

    DO iJ = i1,i2
        daK(iJ) = max(daSpline(iJ-i1+1),0.0d0)
    END DO

    RETURN
    end SUBROUTINE DoBoxCar

!************************************************************************
! do we stick to linemix, or switch to Cousin?
! for testing purposes with Genln2, we have parameter iTestGenln2 = +1
! where we swtich EVERYTHING to Cousin; else the default is set at
! iTestGenln2 = -1 so we use linemixing for the strongest bands
    SUBROUTINE CousinVsMix(iaLineMix,iLineMixBand,iJL,iJU,iISO,daFreq, &
    daLineShift,daStren296,dVibCenter,iNum)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iLineMixBand                  !do we use linemixing coeffs
    DOUBLE PRECISION :: daFreq(kMaxPts)      !the wavevector
    DOUBLE PRECISION :: daLineShift(kHITRAN) !line centers
    DOUBLE PRECISION :: daStren296(kHITRAN)  !line strengths at 296 K
    INTEGER :: iNum                          !number of lines in band
    DOUBLE PRECISION :: dVibCenter           !vib center of band
    INTEGER :: iISO                          !isotope
    INTEGER :: iJL                           !lower quantum number
    INTEGER :: iJU                           !upper quantum number
! output parameter
    INTEGER :: iaLineMix(kHITRAN)            !individual lines linemix coeffs
! et to +2 for linemix, +1 for cous
! local vars
    INTEGER :: iOut,iI,iStrong,iaWeak(kHITRAN),iSwitched
    INTEGER :: iTestGenln2,iDefault
    DOUBLE PRECISION :: d1,d2,dLineShift
    DOUBLE PRECISION :: dfL_Pbranch,dfR_Pbranch,dfL_Rbranch,dfR_Rbranch
    DOUBLE PRECISION :: dSmin,dSmax,dDelta,dThreshold

    iSwitched   = -1

    iDefault    = -1       !!default; use line mix for strong bands
    iTestGenln2 = +1       !!test with GENLN2, so use Cousin everywhere
    iTestGenln2 = -1       !!default; use line mix for strong bands

    IF (iTestGenln2 /= iDefault) THEN
        print *,'iTestGenln2,iDefault = ',iTestGenln2,iDefault
    END IF
          
    dsMin = +1.0d36
    dsMax = -1.0d36
    dThreshold = 1.0d0/2.0d0    !!!make sure this is  0 <= threshold <= 1
    dThreshold = 0.0d0

    DO iI = 1,iNum
        iaWeak(iI) = +1      !!!!assume line is "weak"
        IF (daStren296(iI) <= dSmin) dSmin = daStren296(iI)
        IF (daStren296(iI) >= dSmax) dSmax = daStren296(iI)
    END DO
    dDelta = log10(dSmax)-log10(dSmin)
    dThreshold = log10(dSmin) + dThreshold*dDelta
     
    DO iI = 1,iNum
        IF (log10(daStren296(iI)) > dThreshold) iaWeak(iI) = -1
    END DO

    d1 = daFreq(1)
    d2 = daFreq(kMaxPts)

!!!stick to linemix close to/inside the 2350,2351 and 2310,2320 bands
!!!else switch to Cousin

    DO iI = 1,iNum

        iOut = 2
        iStrong = -1

        dLineShift = daLineShift(iI)

        IF ((iLineMixBand == 2) .AND. (iJU == 9) .AND. (iISO == 1)) THEN
        ! strongest sigsig band  ... do linemix here, epecially for R branch
        ! P branch from 2150 to 2387
            dfL_Pbranch = 2155.0d0
            dfR_Pbranch = 2380.0d0
        ! R branch from 2328 to 2436 cm-1
            dfL_Rbranch = 2330.0d0
            dfR_Rbranch = 2430.0d0
            iStrong = +1
        ELSEIF ((iLineMixBand == 2) .AND. (iJU == 9) .AND. (iISO == 2)) THEN
        ! next strongest sigsig band (isotope)
        ! P branch from 2167.036826 to 2296.920256
            dfL_Pbranch = 2180.0d0
            dfR_Pbranch = 2280.0d0
        ! R branch from 2269.261676 to 2345.936568
            dfL_Rbranch = 2280.0d0
            dfR_Rbranch = 2355.0d0
            iStrong = +1
            IF (iaWeak(iI) > 0) iStrong = -1   !!!weak line in band
        ELSEIF ((iLineMixBand == 2) .AND. (iJU == 24) .AND. (iISO == 1)) THEN
        ! strongest deltdelt band
        ! P branch from 2227.255209 to 2321.797519
            dfL_Pbranch = 2230.0d0
            dfR_Pbranch = 2330.0d0
        ! R branch from 2326.429307 to 2369.735017
            dfL_Rbranch = 2330.0d0
            dfR_Rbranch = 2380.0d0
            iStrong = +1
            IF (iaWeak(iI) > 0) iStrong = -1   !!!weak line in band
        ELSEIF ((iLineMixBand == 2) .AND. (iJU == 16) .AND. (iISO == 1)) THEN
        ! strongest pipi band
        ! P branch from 2227.569529 2335.086137
            dfL_Pbranch = 2230.0d0
            dfR_Pbranch = 2330.0d0
        ! R branch from 2338.151553 2381.80012
            dfL_Rbranch = 2330.0d0
            dfR_Rbranch = 2380.0d0
            iStrong = +1
            IF (iaWeak(iI) > 0) iStrong = -1   !!!weak line in band
        END IF

        IF ((iStrong == 1) .AND. (iTestGenln2 == 1)) THEN
            iSwitched = +1
            iStrong = -1
        END IF

        IF ((iLineMixBand == 2) .AND. (iStrong > 0)) THEN
        ! current 10000 point chunk is way outside the band
            IF ((d1 <= dfL_Pbranch) .AND. (dLineShift <= dVibCenter)) THEN
                iOut = 1             !!!!switch to Cousin
            ELSEIF ((d2 >= dfR_Pbranch) .AND. (dLineShift <= dVibCenter))THEN
                iOut = 1             !!!!switch to Cousin
            ELSEIF ((d1 <= dfL_Rbranch) .AND. (dLineShift >= dVibCenter))THEN
                iOut = 1             !!!!switch to Cousin
            ELSEIF ((d2 >= dfR_Rbranch) .AND. (dLineShift >= dVibCenter))THEN
                iOut = 1             !!!!switch to Cousin
            END IF
                      
        !! for some odd reason, do all the P branch lines for 2350 sigsig,
        !! using the Cousin lineshape
        !! this is because at large pressure I notice FULLMIX << FIRSTORDER
            IF ((dLineShift <= dVibCenter)) THEN
                iOut = 1             !!!!switch to Cousin
            END IF
                        
        ELSEIF ((iLineMixBand == 2) .AND. (iStrong < 0)) THEN
        !! current line is one of the weaker lines in the band OR
        !! weaker bands in NLTE ... just use Cousin
            iOut = 1
        END IF
                  
        iaLineMix(iI) = iOut
    END DO

    IF (iSwitched == +1) THEN
        write(kStdWarn,*) 'testing against GENL2; reset from linemix to cousin'
    END IF
           
    RETURN
    end SUBROUTINE CousinVsMix

!************************************************************************
! this subroutine calls the voigt function at whatever resolution
    SUBROUTINE voigt_chi(daFreq,iFreqPts,dLineShift, &
    dLTE,dMass,dBroad,iLineMix,dP,dPP, &
    dJL,dJU,dJLowerQuantumRot,iISO,dVibCenter, &
    xBirn,chiBirn,iNptsBirn,daLineshape,daFudge,iDoVoigtChi,iTooFar)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iFreqPts   !number of relevant points in daFreq
    INTEGER :: iLineMix   !do we use linemixing coeffs or cousin
    DOUBLE PRECISION :: daFreq(kMaxPts)                !the wavevector
    DOUBLE PRECISION :: dP,dPP                         !layer pressure
    DOUBLE PRECISION :: dLineShift,dMass,dBroad        !line parameters
    DOUBLE PRECISION :: dLTE,dJL,dJU,dJLowerQuantumRot
    DOUBLE PRECISION :: dVibCenter                     !vib center of band
    INTEGER :: iISO
    DOUBLE PRECISION :: daFudge(kMaxPtsBox)  !fudge so that we mimic UMBC-LBL
    INTEGER :: iDoVoigtChi   !!set this THE SAME in SetRunningMesh,voigt_chi
! this is for birnbaum, if needed
    INTEGER :: iNptsBirn
    DOUBLE PRECISION :: xBirn(kMaxPts),chiBirn(kMaxPts)
! output parameters
    DOUBLE PRECISION :: daLineshape(kMaxPtsBox)    ! lineshape, at high res
    INTEGER :: iToofar                             ! is the input linecenter too far from linecenters
! found in linemix file (-1 No, +x yes == in(dJLowerQuantumRot))
! if larger than 100, high J so weak line, so dont worry
! local parameters
    DOUBLE PRECISION :: daImag(kMaxPtsBox)                   !for linemixing calcs
    DOUBLE PRECISION :: daHighFreqWavenumbers(kMaxPtsBox)    !for computing lineshape, at high res
    DOUBLE PRECISION :: daChi(kMaxPtsBox),dT
    DOUBLE PRECISION :: df,f0,tau2_birn,tau2
    DOUBLE PRECISION :: daYmix(kHITRAN),daYmixALL(kHITRAN),dY
    INTEGER :: iFr,iN,iLine,iNum

!      iDoVoigtChi = -1  !! do not multiply by chi function in voigt_chi; all
!                        !!   necesary things will be done by calling
!                        !!   Co2_4um_fudge_nlte_fast
!      iDoVoigtChi = +1  !! do multiply by chi function in voigt_chi

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

!!!    write(kStdErr,*)  ' >>>>>>>>>>> check this voigt stuff before making NLTE production run!!!! <<<<<<<<<<<<<'
!!!    write(kStdWarn,*) ' >>>>>>>>>>> check this voigt stuff before making NLTE production run!!!! <<<<<<<<<<<<<'
    CALL DoVoigt(daLineshape,daImag,daHighFreqWavenumbers,dLineShift,dLTE,dMass,dBroad,iN,dP,dPP)
!!! this is NEW Jan 2017
!!! this is NEW Jan 2017

    IF (iLineMix == 1) THEN    !this does cousin
        dT = dLTE
        CALL Cousin(daChi,daHighFreqWavenumbers,dLineShift,dBroad,dT,dP,dPP,1,iN)
        DO iFr = 1,iN
            daLineshape(iFr) = daLineshape(iFr) * daChi(iFr)
        END DO
    END IF

    IF (iLineMix == 2) THEN   !this does linemix
      dT = dLTE
      iN = iFreqPts	
    ! compare linemix vs linemixALL
    !        CALL linemix(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)
    !        CALL linemixALL(dLineShift,daYmixALL,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)
    !      do iN = 1,120
    !          print *,iN,dLineShift,daYmix(iN),daYmixALL(iN),daYmixALL(iN)/daYmix(iN)
    !      end do
    !      call dostop
        CALL linemix(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dT,dP)

    ! c        tau2 = tau2_birn(dP,dPP)
    ! c        CALL birnbaum(daChi,daFreq,dLineShift,dBroad,dT,tau2,iN)
    !      print *,'doing birnbaum_interp ',dJL,dJU,iISO,daFreq(1)
        CALL birnbaum_interp(daChi,daHighFreqWavenumbers,dLineShift,iN,chiBirn,xBirn,iNptsBirn)
              
        dY = daYmix(iLine)*dP
        DO iFr = 1,iN
            daLineshape(iFr) = (daLineshape(iFr)+daImag(iFr)*dY)*daChi(iFr)
        END DO
        IF ((dJU == 9) .AND. (iISO == 1) .AND. (iDoVoigtChi > 0)) THEN
        ! oh oh check to see if we need to adjust the linemix to bring it
        ! approximnately equal to UMBC-LBL
            IF ((daFreq(1) <= 2505.0d0) .AND. (daFreq(iN) >= 2355.0d0)) THEN
                DO iFr = 1,kMaxPts
                    daLineshape(iFr) = daLineshape(iFr)*daFudge(iFr)
                END DO
            END IF
        END IF
    END IF

    RETURN
    end SUBROUTINE voigt_chi

!************************************************************************
! this is the fudger to bring strongest R branch linemix in line with UMBC-LBL
! need to interpolate the "chi" functions onto dafreq
    SUBROUTINE co2_4um_nlte_fudge(daFreq,daFudge,iN)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iN             !number of points
    DOUBLE PRECISION :: daFreq(kMaxPtsBox)      !the wavevector
! input/output parameters
    DOUBLE PRECISION :: daFudge(kMaxPtsBox)     !fudge factor interped to daFreq

! local variables
    INTEGER :: iFileStartFr  !!its ok to keep this as integer
    INTEGER :: iIOUN,iERR,iChi,iFr,iChunk
    CHARACTER(120) :: fname
    DOUBLE PRECISION :: daF(kMaxPts),daChi(kMaxPts)

    iChi = -1
    iFileStartFr = int(daFreq(1)+0.25)   !!! the +0.25 ensures it is in the
!!! correct chunk

    IF ((iFileStartFR >= 2355) .AND. (iFileSTartFr <= 2505)) THEN
        IF ((iFileStartFR >= 2355) .AND. (iFileSTartFr < 2380)) THEN
            iChunk = 2355
            FNAME = 'nonlte_co2_fudge_2355.txt'
            iChi = +1
        ELSEIF ((iFileStartFR >= 2380) .AND. (iFileSTartFr < 2405)) THEN
            iChunk = 2380
            FNAME = 'nonlte_co2_fudge_2380.txt'
            iChi = +1
        ELSEIF ((iFileStartFR >= 2405) .AND. (iFileSTartFr < 2430)) THEN
            iChunk = 2405
            FNAME = 'nonlte_co2_fudge_2405.txt'
            iChi = +1
        ELSEIF ((iFileStartFR >= 2430) .AND. (iFileSTartFr < 2455)) THEN
            iChunk = 2430
            FNAME = 'nonlte_co2_fudge_2430.txt'
            iChi = +1
        ELSEIF ((iFileStartFR >= 2455) .AND. (iFileSTartFr < 2480)) THEN
            iChunk = 2455
            FNAME = 'nonlte_co2_fudge_2455.txt'
            iChi = +1
        ELSEIF ((iFileStartFR >= 2480) .AND. (iFileSTartFr < 2505)) THEN
            iChunk = 2480
            FNAME = 'nonlte_co2_fudge_2480.txt'
            iChi = +1
        END IF
    END IF

    IF (iChi > 0) THEN
        CALL FindChiFileName(fname)
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
        IF (IERR /= 0) THEN
            WRITE(kStdErr,*) 'In subroutine co2_4um_nlte_fudge'
            WRITE(kStdErr,1010) IERR, FNAME
            1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
            CALL DoSTOP
        ENDIF
        kTempUnitOpen = 1
        READ(iIOUN,*) (daF(iFr),daChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen = -1

        CALL dspl(daF,daChi,kMaxPts,daFreq,daFudge,iN)
    END IF

    RETURN
    end SUBROUTINE co2_4um_nlte_fudge

!************************************************************************
! this is the subroutine that actually calls the voigt fcn
! kinda dumb to call this instead of calling vhh2RI direct, but in a past time
! this subroutine used to call cousin or birnbaum/linem code as well
    SUBROUTINE DoVoigt(daReal1,daImag1,daFreq1,dLineShift,dLTE, &
    dMass,dBroad,iN,dP,dPP)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
    INTEGER :: iN                                       !how many points
    DOUBLE PRECISION :: daFreq1(kMaxPtsBox)             !the wavevector
    DOUBLE PRECISION :: dLineShift,dMass,dBroad         !line parameters
    DOUBLE PRECISION :: dP,dPP                          !layer pressure
    DOUBLE PRECISION :: dLTE
! output parameters
    DOUBLE PRECISION :: daReal1(kMaxPtsBox)             !lineshape (real part)
    DOUBLE PRECISION :: daImag1(kMaxPtsBox)             !lineshape (imag part)

! local variables
    INTEGER :: iFr

    CALL vhh2RI(daReal1,daImag1,daFreq1,dLineShift,dLTE,dMass,dBroad,iN)

!c for testing : turn off lineshape +/- 25 cm away!
!c this is bad!
!c      DO iFr = 1,iN
!c        IF (abs(dLineShift - daFreq1(iFr)) .GT. 25.0d0) THEN
!c          daReal1(iFr) = 0.0d0
!c          daImag(iFr) = 0.0d0
!c        END IF
!c      END DO
!c this is bad!
!c for testing : turn off lineshape +/- 25 cm away!

    RETURN
    end SUBROUTINE DoVoigt

!************************************************************************
! this is a quick approx to tanh
    DOUBLE PRECISION FUNCTION tanhsergio(x)

    DOUBLE PRECISION :: x
    DOUBLE PRECISION :: y

    y = 1.0D0 !default value

    IF (abs(x) < 6.5d0) y = tanh(x)

    tanhsergio = y

    RETURN
    END FUNCTION tanhsergio

!************************************************************************
! this is the /home/sergio/SPECTRA/FORTRANLINUX/vhh2RI.f code
! (used by run7.m)
! this is a Voigt-VanHuber function
    subroutine vhh2RI(wr,wi,v,v0,Temp,m,brd,n_in)
! this is the CURRENT GENLN2 routine
! this is the same as voivec2.f in ../Genln2

! subroutine voigt1(wr,wi,0,Temp,m,brd,n_in)
! gets the real and imag parts of the fcn
! v    = frequency array
! v0   = center freq
! T    = temperature
! m    = molecular mass (amu)
! brd  = broadening

    include '../INCLUDE/kcartaparam.f90'

    DOUBLE PRECISION :: wr(kMaxPtsBox),v(kMaxPtsBox),v0,Temp,m,brd
    DOUBLE PRECISION :: wi(kMaxPtsBox)

    integer :: n_in

! PROGRAM        VOIVEC     SUBROUTINE

! PURPOSE        COMPUTE COMPLEX PROBABILITY FUNCTION

! VERSION        3.X   D.P. EDWARDS   28/05/92

! DESCRIPTION    THIS ROUTINE CALCULATES THE COMPLEX PROBABILITY
!                FUNCTION USING A VECTORIZED VERSION OF THE
!                HUMLICEK JQSRT V27 437 1982 PAPER.
!                THE CALCULATION IS PERFORMED FOR THE ARRAY OF X,Y
!                PAIRS FOR A GIVEN LINE OVER THE FINE MESH POINTS
!                OF THE CURRENT WIDE MESH.

! ARGUMENTS      NUM    I*4 I/P NUMBER OF FINE GRID INTERVALS
!                NL     I*4 I/P FREQUENCY BDY TO START CALCULATION
!                NH     I*4 I/P FREQUENCY BDY TO STOP CALCULATION
!                NSTP   I*4 I/P FREQUENCY BDY STEP
!                 FF     DP  I/P FINE WAVENUMBER GRID [cm-1]
!                 WNUM   DP  I/P WAVENUMBER OF LINE CENTRE [cm-1]
!                 REPWID R*4 I/P SQRT(ln2)/(DOPPLER WIDTH [cm-1])
!                 Y      R*4 I/P VOIGT Y PARAMETER OF LINE
!                 VT     COM O/P COMPLEX PROBABILITY FUNCTION

! local variables
    DOUBLE PRECISION :: factor,c2
    integer :: ii,iDoVanHuber
    complex*16 U,T,VTP(kMaxPtsBox),VTM(kMaxPtsBox)

! local variables

    DOUBLE PRECISION :: k,c_light,amu,mass,r2,alpha_doppler,g0
    DOUBLE PRECISION :: repwid,Y,X,S1V,S2V,fact,c2_0

    IF (n_in > kMaxPtsBox) THEN
        write(kStdErr,*)'in vhh2RI.f, n_in > kMaxPtsBox',n_in,kMaxPtsBox
        Call DoStop
    END IF

!  SORT THE (X,Y) PAIRS INTO THE 4 REGIONS OF THE HUMLIcEK
!  EXPRESSIONS OF THE VOIGT LINE PROFILE.


! do the doppler widths first
    k       = kBoltzmann
    c_light = 2.99792458d8         !ms-1
    amu     = 1.6605402d-27            !nucleon mass/kg
    mass    = m                       !change to kg

! alpha_doppler=v0*sqrt(2*log(2)*k*T/mass/c_light/c_light)
! r2=2*log(2)*k/amu
    r2=11526.218d0
    alpha_doppler=v0/c_light*sqrt(r2*Temp/mass)
    repwid=0.8325546d0/alpha_doppler

! do the g0 factor
!  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
! g0=sqrt(log(2)/pi)/alpha_doppler
    g0= 0.83255461 * 0.5641895 / alpha_doppler


!------------------------------  do the v-vi part ---------------------------
! define arrays for the new Voigt fcn
    Y=brd/alpha_doppler*0.83255461

    do 20 II=1,n_in
        X =  (V(II) - V0)*REPWID
        S1V = abs(X) + Y
        S2V = (0.195*abs(X)) - 0.176
        T = dcmplx(Y,-X)
    
    !  FOR REGION 1 OF HUMLIcEK
    
        if (S1V >= 15.0) then
            VTP(II) = T*0.5641896/(0.5+(T*T))
        
        !  REGION 2 OF HUMLIcEK
        
        elseif (S1V >= 5.5) then
            U = T*T
            VTP(II) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
        
        !  REGION 3 OF HUMLIcEK
        
        elseif (Y >= S2V) then
            VTP(II) = &
            (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+ &
            T*.5642236))))/ &
            (16.4955+T*(38.82363+T*(39.27121+ &
            T*(21.69274+T*(6.699398+T)))))
        
        !  REGION 4 OF HUMLIcEK
        
        else
            U = T*T
            VTP(II)=cdexp(U)-T*(36183.31-U*(3321.9905- &
            U*(1540.787-U*(219.0313-U* &
            (35.76683-U*(1.320522-U*.56419))))))/ &
            (32066.6-U*(24322.84-U* &
            (9022.228-U*(2186.181-U*(364.2191- &
            U*(61.57037-U*(1.841439-U)))))))
        endif
    
        U = dcmplx(1.0d0, -1.0d0)
    
    20 END DO



!------------------------------  do the v+vi part ---------------------------
! define arrays for the new Voigt fcn
    Y=brd/alpha_doppler*0.83255461

    do 30 II=1,n_in
        X =  (V(II) + V0)*REPWID
        S1V = abs(X) + Y
        S2V = (0.195*abs(X)) - 0.176
        T = dcmplx(Y,-X)
    
    !  FOR REGION 1 OF HUMLIcEK
    
        if (S1V >= 15.0) then
            VTM(II) = T*0.5641896/(0.5+(T*T))
        
        !  REGION 2 OF HUMLIcEK
        
        elseif (S1V >= 5.5) then
            U = T*T
            VTM(II) = T*(1.410474 + U*.5641896)/(.75 + U*(3.+U))
        
        !  REGION 3 OF HUMLIcEK
        
        elseif (Y >= S2V) then
            VTM(II) = &
            (16.4955+T*(20.20933+T*(11.96482+T*(3.778987+ &
            T*.5642236))))/ &
            (16.4955+T*(38.82363+T*(39.27121+ &
            T*(21.69274+T*(6.699398+T)))))
        
        !  REGION 4 OF HUMLIcEK
        
        else
            U = T*T
            VTM(II)=cdexp(U)-T*(36183.31-U*(3321.9905- &
            U*(1540.787-U*(219.0313-U* &
            (35.76683-U*(1.320522-U*.56419))))))/ &
            (32066.6-U*(24322.84-U* &
            (9022.228-U*(2186.181-U*(364.2191- &
            U*(61.57037-U*(1.841439-U)))))))
        endif
    
        U = dcmplx(1.0d0, -1.0d0)
    
    30 END DO


! --------------- now do the VHH part --------------------------------
! from the voigt subroutine we need some adjustment factors
!C
!C  SQRT(ln 2) = 0.83255461, 1/SQRT(PI) = 0.5641895
!C       REPWID = 0.83255461/DOPWID
!       H0 = REPWID*0.5641895*STRPAR   ---> STRPAR is multiplied outside
!       Y = WIDPAR*REPWID
!C
!C  CALCULATE THE COMPLEX PROBABILITY FUNCTION
!C
!       CALL VOIVEC(NUM,NL,NH,NSTP,FF,WNUM,REPWID,Y,VT)
!C
!C  COMPUTE ABSORPTION DUE TO VOIGT LINE SHAPE
!C
!       DO 10 IP=NL,NH,NSTP
!         XABS(IP) = H0*REAL(VT(IP))
! 10        CONTINUE
!C

    c2=1.4387863d0        !K/ cm-1  from Genln2 manual
     
    c2_0=0.5d0*c2*v0/Temp
    c2_0=v0*tanhsergio(c2_0)
    fact=repwid*0.5641895d0/c2_0
     
    c2=0.5*c2/Temp

    iDoVanHuber = -1         !!!          do plain voigt
    iDoVanHuber =  0         !!!          do Lorentz
    iDoVanHuber = +1         !!! default, do VanHuber

    iDoVanHuber = +1

    IF (iDoVanHuber /= +1) THEN
        print *,'doing plain voigt or lorentz, NOT vhh',iDoVanHuber
    END IF

    IF (iDoVanHuber > 0) THEN
    !!! do VanHuber
        DO  ii = 1,n_in
            factor  = fact*v(ii)*tanhsergio(c2*v(ii))
            vtp(ii) = factor*(vtp(ii)+vtm(ii))
            wr(ii)  = dreal(vtp(ii))
            wi(ii)  = dimag(vtp(ii))
        END DO

    ELSEIF (iDoVanHuber < 0) THEN
    !!!! just do voigt; bad at lower wavenumbers!
        DO  ii = 1,n_in
            factor  = fact*v(ii)
            wr(ii) = factor*dreal(vtp(ii))
            wi(ii) = factor*dimag(vtp(ii))
        END DO

    ELSEIF (iDoVanHuber == 0) THEN
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
    end subroutine vhh2RI

!************************************************************************
END MODULE kvoigt_cousin
