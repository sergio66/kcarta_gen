! Copyright 2003
! University of Maryland Baltimore County
! All Rights Reserved

! this has the linemixing code
!************************************************************************
! this  subroutine finds the linemix files corresponding to dLTE, dJU
! not only does linemix coeffs depend on T, but they weakly depend on the ppmv!
! so watch out for this later on!
    SUBROUTINE FindLineMixFiles(dLTE,dJU,iISO,caFName1,caFName2,dT1,dT2)

    include '../INCLUDE/kcartaparam.f90'

! input parameters
    DOUBLE PRECISION :: dLTE,dJU      !local kinetic temp, upper state quantum
    INTEGER :: iISO                   !isotope
! output parameters
    CHARACTER(120) :: caFName1, caFName2
    DOUBLE PRECISION :: dT1,dT2
        
! local variables
    INTEGER :: iI,iJ,iFound
    CHARACTER(30) :: caFName,caT1,caT2
    CHARACTER(3) :: ca3
    DOUBLE PRECISION :: dT

    DO iI = 1,120
        caFName1(iI:iI) = ' '
        caFName2(iI:iI) = ' '
    END DO

    DO iI = 1,30
        caT1(iI:iI) = ' '
        caT2(iI:iI) = ' '
    END DO

    caFName1 = caLineMixDir
    caFName2 = caLineMixDir

    IF ((iISO == 1) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFName = 'band2350_linemix.dat'
    ELSEIF ((iISO == 2) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFName = 'band2351_linemix.dat'
    ELSEIF ((iISO == 3) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFName = 'band2352_linemix.dat'
    ELSEIF ((iISO == 1) .AND. (nint(sngl(dJU)) == 16)) THEN
        caFName = 'band2320_linemix.dat'
    ELSEIF ((iISO == 1) .AND. (nint(sngl(dJU)) == 24)) THEN
        caFName = 'band2310_linemix.dat'
    ELSE
        write(kStdErr,*) 'Unidentified CO2 isotope for linemixing'
        write(kStdErr,*) 'iISO,dJU = ',iISO,dJU
        CALL DoStop
    END IF

    iFound = -1
    IF (dLTE < 150.0d0+0.001d0) THEN
        caT1 = '150K_'
        caT2 = '160K_'
        dT1 = 150.0d0
        dT2 = 160.0d0
        iFound = +1
    ELSE IF (dLTE >= 400.0d0-0.001d0) THEN
        caT1 = '390K_'
        caT2 = '400K_'
        dT1 = 390.0d0
        dT2 = 400.0d0
        iFound = +1
    ELSE
        iI = 1
        11 CONTINUE
        dT1 = 150.0d0 + (iI-1)*10.0d0
        dT2 = 150.0d0 + (iI)*10.0d0
        IF ((dLTE >= dT1) .AND. (dLTE < dT2) .AND. (iFound < 0)) THEN
            iFound = +1
        ELSE
            iI = iI + 1
            GOTO 11
        END IF
    END IF

    IF (iFound < 0) THEN
        write (kStdErr,*) 'Could not make bounding linemix filenames!'
        Call DoStop
    END IF

! his now assumes iFound .GT. 0
    12 FORMAT(I3)
    WRITE(ca3,12) nint(dT1)
    caT1(1:3) = ca3(1:3)
    caT1(4:5) = 'K_'
    WRITE(ca3,12) nint(dT2)
    caT2(1:3) = ca3(1:3)
    caT2(4:5) = 'K_'

! put everything together
    iI = 1
    10 CONTINUE
    IF ((caFname1(iI:iI) /= ' ') .AND. (iI < 120)) THEN
        iI = iI + 1
        GOTO 10
    END IF
    iJ = 1
    15 CONTINUE
    IF ((caT1(iJ:iJ) /= ' ') .AND. (iJ < 30)) THEN
        iJ = iJ + 1
        GOTO 15
    END IF
    iJ = iJ-2

    caFName1(iI:iI+iJ) = caT1(1:iJ+1)
    caFName2(iI:iI+iJ) = caT2(1:iJ+1)

    iI = 1
    20 CONTINUE
    IF ((caFname1(iI:iI) /= ' ') .AND. (iI < 120)) THEN
        iI = iI + 1
        GOTO 20
    END IF
    25 CONTINUE
    IF ((caFName(iJ:iJ) /= ' ') .AND. (iJ < 30)) THEN
        iJ = iJ + 1
        GOTO 25
    END IF
    iJ = iJ-2

    caFName1(iI:iI+iJ) = caFName(1:iJ+1)
    caFName2(iI:iI+iJ) = caFName(1:iJ+1)
      
    RETURN
    end SUBROUTINE FindLineMixFiles

!************************************************************************
! this subroutine reads in the LINEMIX parameters, from TWO small files
    SUBROUTINE linemix(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dLTE,dP)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: dLineShift   !line center
    DOUBLE PRECISION :: dLTE         !local kinetic temperature
    DOUBLE PRECISION :: dP           !layer pressure
    INTEGER :: iISO                  !isotope type
    DOUBLE PRECISION :: dJL,dJU            !!lower and upper quantum vibration numbers
    DOUBLE PRECISION :: dJLowerQuantumRot  !!lower           quantum rotation number
! output parameters
    DOUBLE PRECISION :: daYmix(kHITRAN)   !line mix coeffs
    INTEGER :: iLine                      !which one corresponds to dLineShift
    INTEGER :: iNum                       !number of Ymix lines for this band
    INTEGER :: iTooFar                    !is dLineShift too far from lines in files?
          
! local variables
    INTEGER :: iI,iErr,iIOUN,jLow,jHigh,iGasID,iISOTOPE,iFound,iJL,iJU,iJLowQR
    DOUBLE PRECISION :: daF1(kHITRAN),daF2(kHITRAN)
    DOUBLE PRECISION :: daY1(kHITRAN),daY2(kHITRAN)
    CHARACTER(120) :: caFname1,caFName2
    CHARACTER(30) ::  caX
    DOUBLE PRECISION :: dDiff0,dDiff1,dT1,dT2,dBlah
    CHARACTER(1) :: caPorR_1(kHITRAN),caPorR_2(kHITRAN)
          
    iTooFar = -1         !!assume this line DOES exist in linemix files
    iJL = int(dJL)
    iJU = int(dJU)
    iJLowQR = int(dJLowerQuantumRot)
          
    CALL FindLineMixFiles(dLTE,dJU,iISO,caFName1,caFName2,dT1,dT2)
          
! open first file
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName1,FORM='formatted',STATUS='OLD', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine linemix, error reading'
        write (kStdErr,*) 'file (number 1) that has linemix parameters'
        WRITE(kStdErr,1070) iErr, caFName1
        CALL DoSTOP
    END IF
    kTempUnitOpen = 1

    read(iIOUN,*) iGasID,iNum,iISOTOPE,jLow,jHigh
    IF ((iGasID == 2) .AND. (iISOTOPE == iISO) .AND. &
    (jHigh == nint(dJU))) THEN
        dDiff0 = 1.0d0
    ELSE
        write(kStdErr,*) 'linemix reader : other bands'
        Call DoStop
    END IF

    IF (iNum > kHITRAN) THEN
        write(kStdErr,*) 'File has ',iNum,' line parameters for gas ',iGasID
        write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
        write(kStdErr,*) 'Please check kHITRAN in kcartaparam.f90 and fix'
        CALL DoStop
    END IF

    DO iI = 1,iNum
    !        read(iIOUN,*) jLow,daF1(iI),daY1(iI)
        read(iIOUN,*) jLow,daF1(iI),daY1(iI),caPorR_1(iI)
    END DO

    close (iIOUN)
    kTempUnitOpen = -1

! open second file
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName2,FORM='formatted',STATUS='OLD', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine linemix, error reading'
        write (kStdErr,*) 'file (number 2) that has linemix parameters'
        WRITE(kStdErr,1070) iErr, caFName2
        CALL DoSTOP
    END IF
    kTempUnitOpen = 1

    read(iIOUN,*) iGasID,iNum,iISOTOPE,jLow,jHigh
    IF ((iGasID == 2) .AND. (iISOTOPE == iISO) .AND. &
    (jHigh == nint(dJU))) THEN
        dDiff0 = 1.0d0
    ELSE
        write(kStdErr,*) 'linemix reader : other bands'
        Call DoStop
    END IF

    IF (iNum > kHITRAN) THEN
        write(kStdErr,*) 'File has ',iNum,' line parameters for gas ',iGasID
        write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
        write(kStdErr,*) 'Please check kHITRAN in kcartaparam.f90 and fix'
        CALL DoStop
    END IF

    DO iI = 1,iNum
    !        read(iIOUN,*) jLow,daF2(iI),daY2(iI)
        read(iIOUN,*) jLow,daF2(iI),daY2(iI),caPorR_2(iI)
    END DO

    close (iIOUN)
    kTempUnitOpen = -1

! check that the line centers from the two files are in the same order
! at the same time, just fill in daYmix for kicks
    DO iI = 1,iNum
        dBlah = (daY2(iI)-daY1(iI))/(dT2 - dT1)    !slope between two  curves
        daYmix(iI) = daY1(iI) + dBlah*(dLTE-dT1)   !linear interpolation
    !        print *,dLTE,daF1(iI),daYmix(iI)*dP
        IF (dabs(daF1(iI) - daF2(iI)) > 1.0d-10) THEN
            write(kStdErr,*) 'Lines in files ',caFName1,caFName2,' are unequal!!'
            CALL DoStop
        END IF
    END DO

! now find the closest line to dLineShift, which is the input line
    iLine = 1
    dDiff0 = abs(dLineShift - daF1(iLine))
    DO iI = 1,iNum
        dDiff1 = abs(dLineShift - daF1(iI))
        IF (dDiff1 <= dDiff0) THEN
            iLine = iI
            dDiff0 = dDiff1
        END IF
    END DO
    IF (abs(dDiff0) >= 5.0d-3) THEN
    !!!  old code caused the program to stop
    ! aYmix(iLine) = 0.0d0
    ! TooFar = +1
    ! rint *,'large (df) : line center sent in = ',dLineShift,' smallest diff with lines in file ',abs(dDiff0)
    ! rint *,'dJU,iISO = ',dJU,iISO,' T1,T2 (in K) = ',dT1,dT2, 'first linemix file ',caFName1
    ! ALL DoStop
    !!! new code tries to see how bad rotation vib quantum number is, if > 100, this is a weak line so dont worry
        daYmix(iLine) = 0.0d0
        iTooFar = iJLowQR
    ! print *,'dLineCenter,diff = ',dLineShift,abs(dDiff0),' dJL,dJU,iISO = ',iJL,iJU,iISO,' JRotQuant ',iJLowQR
    END IF

!      write(kStdWarn,*) 'f0,df,iLine,ymix= ',sngl(daF1(iLine)),dDiff0,iLine,sngl(daYmix(iLine)*dP)
         
    1070 FORMAT('ERROR! number ',I5,' opening LINEMIX parameter file:',/,A80)
    1080 FORMAT('Opened LineMix parameter file ',A80)

    RETURN
    end SUBROUTINE linemix

!************************************************************************
! this subroutine reads in the LINEMIX parameters, from ONE big file
    SUBROUTINE linemixALL(dLineShift,daYmix,iTooFar,iNum,iLine,dJL,dJU,dJLowerQuantumRot,iISO,dLTE,dP)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: dLineShift   !line center
    DOUBLE PRECISION :: dLTE         !local kinetic temperature
    DOUBLE PRECISION :: dP           !layer pressure
    INTEGER :: iISO                  !isotope type
    DOUBLE PRECISION :: dJL,dJU            !!lower and upper quantum vibration numbers
    DOUBLE PRECISION :: dJLowerQuantumRot  !!lower           quantum rotation number
! output parameters
    DOUBLE PRECISION :: daYmix(kHITRAN)   !line mix coeffs
    INTEGER :: iLine                      !which one corresponds to dLineShift
    INTEGER :: iNum                       !number of Ymix lines for this band
    INTEGER :: iTooFar                    !is dLineShift too far from lines in files?
          
! local variables
    INTEGER :: iI,iJ,iErr,iIOUN,jLow,jHigh,iGasID,iISOTOPE,iJL,iJU,iJLowQR
    INTEGER :: iFound1,iFound2,iBand,iTemp,iSkip,iBandALL,iT1,iT2
    DOUBLE PRECISION :: daF1(kHITRAN),daF2(kHITRAN)
    DOUBLE PRECISION :: daY1(kHITRAN),daY2(kHITRAN)
    CHARACTER(120) :: caFnameALL
    CHARACTER(30) ::  caX
    DOUBLE PRECISION :: dDiff0,dDiff1,dT1,dT2,dBlah
    CHARACTER(1) :: caPorR_1(kHITRAN),caPorR_2(kHITRAN)
    CHARACTER(80) :: caSkip80
          
    iTooFar = -1         !!assume this line DOES exist in linemix files
    iJL = int(dJL)
    iJU = int(dJU)
    iJLowQR = int(dJLowerQuantumRot)

!! find bounding temperatures for dLTE
    IF (dLTE <= 160.0) THEN
        iFound1 = +1
        iT1 = 150
        iT2 = 160
    ELSEIF (dLTE >= 390.0) THEN
        iFound1 = +1
        iT1 = 390
        iT2 = 400
    ELSE
        iFound1 = -1
        iI = 1
        11 CONTINUE
        dT1 = 150.0d0 + (iI-1)*10.0d0
        dT2 = 150.0d0 + (iI)*10.0d0
        IF ((dLTE >= dT1) .AND. (dLTE < dT2) .AND. (iFound1 < 0)) THEN
            iFound1 = +1
        ELSE
            iI = iI + 1
            GOTO 11
        END IF
        iT1 = iI
        IT2 = IT1 + 1
    END IF
    IF (iFound1 < 0) THEN
        write(kStdErr,*) 'linemixALL : could not find bounding temperatures for input LTE temp ',dLTE
        CALL DoStop
    END IF
          
    IF ((iISO == 1) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFNameALL = 'band2350_linemix.dat'
        iBandALL   = 2350
    ELSEIF ((iISO == 2) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFNameALL = 'band2351_linemix.dat'
        iBandALL   = 2351
    ELSEIF ((iISO == 3) .AND. (nint(sngl(dJU)) == 9)) THEN
        caFNameALL = 'band2352_linemix.dat'
        iBandALL   = 2352
    ELSEIF ((iISO == 1) .AND. (nint(sngl(dJU)) == 16)) THEN
        caFNameALL = 'band2320_linemix.dat'
        iBandALL   = 2320
    ELSEIF ((iISO == 1) .AND. (nint(sngl(dJU)) == 24)) THEN
        caFNameALL = 'band2310_linemix.dat'
        iBandALL   = 2310
    ELSE
        write(kStdErr,*) 'Unidentified CO2 isotope for linemixing'
        write(kStdErr,*) 'iISO,dJU = ',iISO,dJU
        CALL DoStop
    END IF
             
    DO iI = 1,120
        caFNameALL(iI:iI) = ' '
    END DO
    caFNameALL = caLineMixDir
    caX        = '/all_linemixPR.dat'
! put everything together
    iI = 1
    10 CONTINUE
    IF ((caFnameALL(iI:iI) /= ' ') .AND. (iI < 120)) THEN
        iI = iI + 1
        GOTO 10
    END IF
    iJ = 1
    15 CONTINUE
    IF ((caX(iJ:iJ) /= ' ') .AND. (iJ < 30)) THEN
        iJ = iJ + 1
        GOTO 15
    END IF
    caFNameALL(iI:iI+iJ) = caX(1:iJ-1)
          
! open file and search for band, first temperature
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFNameALL,FORM='formatted',STATUS='OLD', &
    IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine linemix, error reading'
        write (kStdErr,*) 'file that has linemix parameters'
        WRITE(kStdErr,1070) iErr, caFNameALL
        CALL DoSTOP
    END IF
    kTempUnitOpen = 1

    iFound1 = -1
    DO WHILE (iFound1 < 0)
        read(iIOUN,*) iGasID,iNum,iISOTOPE,jLow,jHigh,iTemp,iBand
    !        print *,'A',iGasID,iNum,iISOTOPE,jLow,jHigh,iTemp,iBand,'X',iISO,dJU,iBandALL,iT1
        IF ((iGasID == 2) .AND. (iISOTOPE == iISO) .AND. &
        (jHigh == nint(dJU)) .AND. (iTemp == 150+(iT1-1)*10) .AND. (iBand == iBandALL)) THEN
            IF (iNum > kHITRAN) THEN
                write(kStdErr,*) 'File has ',iNum,' line parameters for gas ',iGasID
                write(kStdErr,*) 'Code can only handle ',kHITRAN,' line parameters '
                write(kStdErr,*) 'Please check kHITRAN in kcartaparam.f90 and fix'
                CALL DoStop
            END IF
            iFound1 = +1
            iFound1 = iTemp
            DO iI = 1,iNum
                read(iIOUN,*) jLow,daF1(iI),daY1(iI),caPorR_1(iI)
            END DO
        !! since we are at correct place, go read next chun, as it is at next Temperaturek!!!
            read(iIOUN,*) iGasID,iNum,iISOTOPE,jLow,jHigh,iTemp,iBand
        !          print *,'B',iGasID,iNum,iISOTOPE,jLow,jHigh,iTemp,iBand,'X',iISO,dJU,iBandALL,iT2
            IF ((iGasID == 2) .AND. (iISOTOPE == iISO) .AND. &
            (jHigh == nint(dJU)) .AND. (iTemp == 150+(iT2-1)*10) .AND. (iBand == iBandALL)) THEN
                iFound2 = +1
                iFound2 = iTemp
                DO iI = 1,iNum
                    read(iIOUN,*) jLow,daF2(iI),daY2(iI),caPorR_1(iI)
                END DO
            ELSE
                write(kStdErr,*) 'Very odd : iT2 should follow iT1 !!!!',iFound1,iFound2
                CALL DoStop
            END IF
        ELSE
            DO iI = 1,iNum
                read(iIOUN,*) caSkip80
            END DO
        END IF
    END DO
    close (iIOUN)
    kTempUnitOpen = -1

! check that the line centers from the two files are in the same order
! at the same time, just fill in daYmix for kicks
    DO iI = 1,iNum
        IF (dabs(daF1(iI) - daF2(iI)) > 1.0d-10) THEN
            write(kStdErr,*) 'Line centers for ',iT1,iT2,' are unequal!!'
            CALL DoStop
        END IF
        dBlah = (daY2(iI)-daY1(iI))/(dT2 - dT1)    !slope between two  curves
        daYmix(iI) = daY1(iI) + dBlah*(dLTE-dT1)   !linear interpolation
    !        print *,dLTE,daF1(iI),daYmix(iI)*dP
    END DO

! now find the closest line to dLineShift, which is the input line
    iLine = 1
    dDiff0 = abs(dLineShift - daF1(iLine))
    DO iI = 1,iNum
        dDiff1 = abs(dLineShift - daF1(iI))
        IF (dDiff1 <= dDiff0) THEN
            iLine = iI
            dDiff0 = dDiff1
        END IF
    END DO
    IF (abs(dDiff0) >= 5.0d-3) THEN
    !!!  old code caused the program to stop
    ! aYmix(iLine) = 0.0d0
    ! TooFar = +1
    ! rint *,'large (df) : line center sent in = ',dLineShift,' smallest diff with lines in file ',abs(dDiff0)
    ! rint *,'dJU,iISO = ',dJU,iISO,' T1,T2 (in K) = ',dT1,dT2, 'first linemix file ',caFName1
    ! ALL DoStop
    !!! new code tries to see how bad rotation vib quantum number is, if > 100, this is a weak line so dont worry
        daYmix(iLine) = 0.0d0
        iTooFar = iJLowQR
    ! print *,'dLineCenter,diff = ',dLineShift,abs(dDiff0),' dJL,dJU,iISO = ',iJL,iJU,iISO,' JRotQuant ',iJLowQR
    END IF

!      write(kStdWarn,*) 'f0,df,iLine,ymix= ',sngl(daF1(iLine)),dDiff0,iLine,sngl(daYmix(iLine)*dP)
         
    1070 FORMAT('ERROR! number ',I5,' opening LINEMIX parameter file:',/,A80)
    1080 FORMAT('Opened LineMix parameter file ',A80)

    RETURN
    end SUBROUTINE linemixALL

!************************************************************************
! this just figures out "tau" for birnbaum
    DOUBLE PRECISION FUNCTION tau2_birn(dP,dPP)

    include '../INCLUDE/kcartaparam.f90'

    DOUBLE PRECISION :: dP,dPP

! local var
    DOUBLE PRECISION :: tau2,dur(4),pt(4),ps(4),ratio(4),duration_self,pf
    DOUBLE PRECISION :: junk1(2),junk2(2)
    INTEGER :: iFr

    tau2 = 5.0d-3
          
    tau2 = 5.0d-3
    tau2 = 4.769324235354260D-003    !for CO2 ppmv = 370
    tau2 = 4.63-003    !for some reason, the run7code uses this ppmv

! copied from SPECTRA/CO2_COMMON/co2_param.m
    dur(1) = 10.46423654878965d-3
    dur(2) = 6.573033770334661d-3
    dur(3) = 5.614703762280506d-3
    dur(4) = 5.251501801911941d-3

    pt(1) = 2.6800d+01
    pt(2) = 1.6140d+02
    pt(3) = 5.6030d+02
    pt(4) = 9.6170d+02

    ps(1) = 2.6800d+01
    ps(2) = 2.6800d+01
    ps(3) = 2.6800d+01
    ps(4) = 2.6800d+01

    DO iFr = 1,4
        pt(iFr) = pt(iFr)/kAtm2mb
        ps(iFr) = ps(iFr)/kAtm2mb
        pf = pt(iFr) - ps(iFr)
        ratio(iFr) = pf/pt(iFr)
    END DO

    duration_self = dur(1)

!      CALL dspl(ratio,dur,4,(dP-dPP)/dP,tau2,1)
!      print *,'in tau2_birn : ',dPP,dP,(dP-dPP)/dP,tau2,(ratio(iFr),iFr=1,4)
!      tau2_birn = tau2

! this is for ifort compiler
    junk1(1) = (dP-dPP)/dP
    junk1(2) = (dP-dPP)/dP*1.0001
    junk2(1) = tau2
    junk2(2) = tau2
    CALL dspl(ratio,dur,4,junk1,junk2,2)
!      print *,'in tau2_birn : ',dPP,dP,(dP-dPP)/dP,tau2,(ratio(iFr),iFr=1,4)
    tau2_birn = junk2(1)
          
    RETURN
    END FUNCTION tau2_birn
          
!************************************************************************
! this subroutine finds Birnbaum
! copied directly from SPECTRA/FORTRANLINUX/birnbaum.f
! except that instead of assigning values to chi111,chi122 etc, it
! directly assigns the indices to chif,chip
    subroutine birnbaum(z,v,v0,w_tot,T_in,tau2,n_in)

! subroutine birnbaum(z,v,v0,w_tot,T_in,tau2,n_in)
! this is the birnbaum lineshape using lookup tables
! z    = results array
! v    = frequency array
! v0   = center freq
! T_in    = temperature
! tau2 = duration of collision
! n_in = no of input points

! see birn_lookup.m
! since w_tot is immaterial
! this does a scan of tau2 for 7 temps
! for ii=1:10
!   pow=1;
!   doc=1e-3*(1+(ii-1)*3);
!   y100(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),100,doc);
!   y150(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),150,doc);
!   y200(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),200,doc);
!   y250(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),250,doc);
!   y300(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),300,doc);
!   y350(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),350,doc);
!   y400(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),400,doc);
!   end

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: v(kMaxPtsBox),v0,T_in,tau2,w_tot
    integer :: n_in
! output params
    DOUBLE PRECISION :: z(kMaxPtsBox)

! local variables
    DOUBLE PRECISION :: dd
    integer :: npts,ndoc,ntemp
    parameter(npts=2001,ntemp=7,ndoc=10)

    integer ::           i,j,k,iTau,iTemp,nl,nh,nstp,IOUN,iErr
    DOUBLE PRECISION ::  tautau,fdif,fchi,df,chil,chiu,chif,chip
    DOUBLE PRECISION ::  vbc,vtc,dvc,temp(ntemp),doc(ndoc)

    DOUBLE PRECISION :: chi(npts,ndoc,ntemp) !!for T=100,150,200,250,300,350,400
    CHARACTER(80) :: caFName,caFName1

! we have 10 d of c parameters tau2=1e-3*(1,3,5,...,19) (ndoc=10)
! we have 7 temperatures : 100,150,200,250,300,350,400  (ntemp=7)
! the frequency spacing of the y's is 0.25
! the dnu's range from -250 to +250 in steps of 0.25    (npts=2001)
    data temp/100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0/
    data doc/1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0/
    data vbc,vtc,dvc/-250.0,+250.0,0.25/

    IF (n_in > kMaxPtsBox) THEN
        write(kStdErr,*) 'in birnbaum.f n_in > kMaxPtsBox'
        CALL DoStop
    END IF
    IF ((tau2 < 1e-3) .OR. (tau2 > 28e-3)) THEN
        write(kStdErr,*) 'ooops!! tau2 = ',tau2,'instead of 1e-3 <= x <= 28e-3'
        CALL DoStop
    END IF
    IF (T_in < 100.0) THEN
        write(kStdWarn,*) 'ooops!!! t = ',T_in,' instead of 100 <= T <= 450'
        write(kStdWarn,*) 'setting t = 100K'
        T_in = 100.0
    END IF
    IF (T_in > 400.0) THEN
        write(kStdWarn,*) 'ooops!!! t = ',T_in,' instead of 100 <= T <= 450'
        write(kStdWarn,*) 'setting t = 400K'
        T_in = 400.0
    END IF

! lookup table produced using birn_lookupNEW.m, which is a script file
! previous version read in the text file      include 'birn_lookup.dat'

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! now we load in bbiirrnn.dat
    caFName = caLineMixDir
    i = 1
    1 CONTINUE
    IF ((caFName(i:i) /= ' ') .AND. (i < 80)) THEN
        i = i + 1
        GOTO 1
    END IF
    caFName1 = 'bbiirrnn.dat'
    j = 1
    2 CONTINUE
    IF ((caFName1(j:j) /= ' ') .AND. (j < 80)) THEN
        j = j + 1
        GOTO 2
    END IF
    caFName(i:i+j-1) = caFName1(1:j)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    IOUN = kTempUnit
! this uses the asymmetric functions that Dave Tobin used in his thesis
! see birn_lookupNEW.m
    open(unit=IOUN,FILE=caFName,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine birnbaum_coarse, error reading file'
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    END IF
    1070 FORMAT('ERROR! number ',I5,' opening birnbaum parameter file:',/,A80)
          
    kTempUnitOpen = 1

! this uses the symmetric functions that Scott Hannon thinks to use
! see birn_lookupNEW2.m
!      open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn2.dat',STATUS='OLD',
!     $                                   FORM='UNFORMATTED')

    DO k = 1,ntemp
        read(IOUN) dd
        IF (abs(dd-temp(k)) >= 1e-2) THEN
            write(kStdErr,*) 'hmm : current temp should be ',temp(k),', not ',dd
            CALL DoStop
        END IF
        DO j=1,ndoc
            read(IOUN) (chi(i,j,k),i=1,npts)
        END DO
    END DO

    close(IOUN)
    kTempUnitOpen = -1

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    nl=1
    nh=n_in
    nstp=1

!c      T_in = 300.1d0
!c      tau2 = 4.001d-3
!c      v0 = 0.0d0
!c      nh = kMaxPts
!c      DO i = 1,kMaxPts
!c        v(i) = -200.0d0 + (400.0d0/kMaxPts)*(i-1)
!c      END DO

! find which temps (100,150,200,...,400) bracket the needed temp
    j=-1
    tautau=T_in
    iTemp=1
    i=1
    10 continue
    if ((j < 0) .AND. (tautau <= temp(i))) then
        j=1
        iTemp=i-1
    elseif (i < ntemp) then
        i=i+1
        goto 10
    else
        write(kStdErr,*) 'sorry : temp is outside 100:150:...:400'
        call DoStop
    end if
    if (iTemp < 1) iTemp = 1
          
! find which doc's (1,3,5,7,9,11,...,19) bracket the needed doc
    j=-1
    tautau=tau2/1e-3
    iTau=1
    i=1
    20 continue
    if ((j < 0) .AND. (tautau <= doc(i))) then
        j=1
        iTau=i-1
    elseif (i < ndoc) then
        i=i+1
        goto 20
    else
    !        print *,'sorry : doc is outside 1:3:5:7:9:11:..:19 e-3'
        write(kStdErr,*) 'sorry : doc is outside 1:4:7:..:28 e-3'
        call DoStop
    end if
    if (iTau < 1) iTau = 1

    do i=nl,nh,nstp
        fdif = (v(i) - v0)
        j = min((npts-2),int(fdif/dvc)) + 1
        j = j + 1000   !because we have chi100(1) ==> dnu = -250
        fchi = vbc + (j - 1)*dvc
        df=(fdif-fchi)/dvc

        k = iTemp
        IF (k == ntemp) k = ntemp - 1

    !          chi111=chi350(j,iTau)
    !          chi211=chi350(j+1,iTau)
    !          chi121=chi350(j,iTau+1)
    !          chi221=chi350(j+1,iTau+1)
    !          chi112=chi400(j,iTau)
    !          chi212=chi400(j+1,iTau)
    !          chi122=chi400(j,iTau+1)
    !          chi222=chi400(j+1,iTau+1)

    !       do the interpolation in Tau, at temp iTemp
    !        chil=chi111 + (chi211-chi111)*df
    !        chiu=chi121 + (chi221-chi121)*df
        chil=chi(j,iTau,k)   + (chi(j+1,iTau,k)-chi(j,iTau,k))*df
        chiu=chi(j,iTau+1,k) + (chi(j+1,iTau+1,k)-chi(j,iTau+1,k))*df
        chif = chil + (chiu-chil)*(tautau-doc(iTau))/(doc(iTau+1)-doc(iTau))

    !       do the interpolation in Tau, at temp iTemp+1
    !        chil=chi112 + (chi212-chi112)*df
    !        chiu=chi122 + (chi222-chi122)*df
        chil=chi(j,iTau,k+1)   + (chi(j+1,iTau,k+1)-chi(j,iTau,k+1))*df
        chiu=chi(j,iTau+1,k+1) + (chi(j+1,iTau+1,k+1)-chi(j,iTau+1,k+1))*df
        chip = chil + (chiu-chil)*(tautau-doc(iTau))/(doc(iTau+1)-doc(iTau))
         
    !        do the interpolation in temp
        z(i)=chif+(chip-chif)*(T_in-temp(iTemp))/(temp(iTemp+1)-temp(iTemp))
    enddo

    RETURN
    end subroutine birnbaum

!************************************************************************
! copied directly from SPECTRA/FORTRANLINUX/birnbaum.f
! uses Numerical Recipes ideas for 2d interpolation
! since ALL lines, for a given layer, should need same birn(T,tau,v-v0), the
! only difference is that v-v0 changes as line center v0 changes.
! IDEA : just output coarse (xx,zz) for ALL lines, and then finely interpolate
!        as necessary
    subroutine birnbaum_coarse(chiBirn,xBirn,iNptsBirn,T_in,tau2)

! subroutine birnbaum_coarse(zz,xBirn,T_in,tau2)
! this is the birnbaum lineshape using lookup tables
! chiBirn     = output results array
! xBirn       = output frequency array
! T_in     = temperature
! tau2     = duration of collision

! see birn_lookup.m
! since w_tot is immaterial
! this does a scan of tau2 for 7 temps
! for ii=1:10
!   pow=1;
!   doc=1e-3*(1+(ii-1)*3);
!   y100(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),100,doc);
!   y150(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),150,doc);
!   y200(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),200,doc);
!   y250(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),250,doc);
!   y300(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),300,doc);
!   y350(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),350,doc);
!   y400(ii,:)=birnbaum2(000:0.25:500,250,10^(-pow),400,doc);
!   end

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params
    DOUBLE PRECISION :: T_in,tau2
! output params
    DOUBLE PRECISION :: chiBirn(kMaxPts)
    DOUBLE PRECISION :: xBirn(kMaxPts)
    INTEGER :: iNptsBirn

! local variables
    DOUBLE PRECISION :: dd
    integer :: npts,ndoc,ntemp
    parameter(npts=2001,ntemp=7,ndoc=10)

    integer ::           i,j,k,iTau,iTemp,IOUN,iERR
    DOUBLE PRECISION ::  tautau,fdif,fchi,df
    DOUBLE PRECISION ::  vbc,vtc,dvc,temp(ntemp),doc(ndoc)

    DOUBLE PRECISION :: chi(npts,ndoc,ntemp) !!for T=100,150,200,250,300,350,400
    CHARACTER(80) :: caFName,caFName1
    DOUBLE PRECISION :: t,u,y1,y2,y3,y4,v(kMaxPts)

! we have 10 d of c parameters tau2=1e-3*(1,3,5,...,19) (ndoc=10)
! we have 7 temperatures : 100,150,200,250,300,350,400  (ntemp=7)
! the frequency spacing of the y's is 0.25
! the dnu's range from -250 to +250 in steps of 0.25    (npts=2001)
    data temp/100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0/
    data doc/1.0, 4.0, 7.0, 10.0, 13.0, 16.0, 19.0, 22.0, 25.0, 28.0/
    data vbc,vtc,dvc/-250.0,+250.0,0.25/

    IF ((tau2 < 1e-3) .OR. (tau2 > 28e-3)) THEN
        write(kStdErr,*) 'ooops!! tau2 = ',tau2,'instead of 1e-3 <= x <= 28e-3'
        write(kStdErr,*) 'temperature = ',T_in
        call DoStop
    END IF
    IF (T_in < 100.0) THEN
        write(kStdWarn,*) 'ooops!!! t = ',T_in,' instead of 100 <= T <= 450'
        write(kStdWarn,*) 'setting t = 100K'
        T_in = 100.0
    END IF
    IF (T_in > 400.0) THEN
        write(kStdWarn,*) 'ooops!!! t = ',T_in,' instead of 100 <= T <= 450'
        write(kStdWarn,*) 'setting t = 400K'
        T_in = 400.0
    END IF

! lookup table produced using birn_lookupNEW.m, which is a script file
! previous version read in the text file      include 'birn_lookup.dat'

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! load in bbiirrnn.dat
    caFName = caLineMixDir
    i = 1
    1 CONTINUE
    IF ((caFName(i:i) /= ' ') .AND. (i < 80)) THEN
        i = i + 1
        GOTO 1
    END IF
    caFName1 = 'bbiirrnn.dat'
    j = 1
    2 CONTINUE
    IF ((caFName1(j:j) /= ' ') .AND. (j < 80)) THEN
        j = j + 1
        GOTO 2
    END IF
    caFName(i:i+j-1) = caFName1(1:j)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    IOUN = kTempUnit
! this uses the asymmetric functions that Dave Tobin used in his thesis
! see birn_lookupNEW.m
    OPEN(unit=IOUN,FILE=caFName,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=iErr)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine birnbaum_coarse, error reading file'
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    END IF
    1070 FORMAT('ERROR! number ',I5,' opening birnbaum coarse parameter file:',/,A80)
     
    kTempUnitOpen = 1

! this uses the symmetric functions that Scott Hannon thinks to use
! see birn_lookupNEW2.m
!      open(unit=IOUN,FILE='../FORTRANLINUX/bbiirrnn2.dat',STATUS='OLD',
!     $                                   FORM='UNFORMATTED')

    DO k = 1,ntemp
        read(IOUN) dd
        IF (abs(dd-temp(k)) >= 1e-2) THEN
            write(kStdErr,*) 'hmm : current temp should be ',temp(k),', not ',dd
            CALL DoStop
        END IF
        DO j=1,ndoc
            read(IOUN) (chi(i,j,k),i=1,npts)
        END DO
    END DO

    close(IOUN)
    kTempUnitOpen = -1

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! find which temps (100,150,200,...,400) bracket the needed temp

!c      T_in = 300.1d0
!c     tau2 = 4.001d-3
!c      T_in = 200.1d0
!c      tau2 = 4.73001d-3
!c      v0 = 0.0d0
!c      nh = kMaxPts
!c      DO i = 1,kMaxPts
!c        v(i) = -200.0d0 + (400.0d0/kMaxPts)*(i-1)
!c      END DO

    j=-1
    tautau=T_in
    iTemp=1
    i=1
    10 continue
    if ((j < 0) .AND. (tautau <= temp(i))) then
        j=1
        iTemp=i-1
    elseif (i < ntemp) then
        i=i+1
        goto 10
    else
        write(kStdErr,*) 'sorry : temp is outside 100:150:...:400'
        call DoStop
    end if
    if (iTemp < 1) iTemp = 1
          
! find which doc's (1,3,5,7,9,11,...,19) bracket the needed doc
    j=-1
    tautau=tau2/1e-3
    iTau=1
    i=1
    20 continue
    if ((j < 0) .AND. (tautau <= doc(i))) then
        j=1
        iTau=i-1
    elseif (i < ndoc) then
        i=i+1
        goto 20
    else
    !        print *,'sorry : doc is outside 1:3:5:7:9:11:..:19 e-3'
        write(kStdErr,*) 'sorry : doc is outside 1:4:7:..:28 e-3'
        call DoStop
    end if
    if (iTau < 1) iTau = 1

    DO i = 1,npts
    ! chi(npts,ndoc,ntemp) !!for T=100,150,200,250,300,350,400
        y1 = chi(i,iTau,iTemp)
        y2 = chi(i,iTau+1,iTemp)
        y3 = chi(i,iTau+1,iTemp+1)
        y4 = chi(i,iTau,iTemp+1)
        t  = (tautau - doc(iTau))/(doc(iTau+1)-doc(iTau))
        u  = (T_in- temp(iTemp))/(temp(iTemp+1)-temp(iTemp))
        xBirn(i)   = vbc + dvc*(i-1)
        chiBirn(i) = (1.0d0 - t)*(1.0d0 - u)*y1 + t*(1.0d0 - u)*y2 + &
        t*u*y3 + (1.0d0 - t)*u*y4
    !c        print *,i,sngl(T_in),sngl(tautau),sngl(vbc + dvc*(i-1)),
    !c     $          sngl(chi(i,iTau,iTemp)),sngl(chiBirn(i))
    END DO

    iNptsBirn = npts

!c       stop

    RETURN
    end subroutine birnbaum_coarse

!************************************************************************
! uses data from birnbaum_coarse to quickly interp stuff to v
    subroutine birnbaum_interp(z,v,v0,n_in,chiBirn,xBirn,iNptsBirn)

! subroutine birnbaum(z,v,v0,n_in,chiBirn,xBirn,iNptsBirn)
! this is the birnbaum lineshape using lookup tables
! z       = results array
! v       = frequency array
! v0      = center freq
! n_in    = no of input points
! chiBirn,xBirn,iNptsBirn comes from birnbaum_coarse

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input params : input wavenumber array, linecenter
    INTEGER :: n_in
    DOUBLE PRECISION :: v(kMaxPtsBox),v0
! input parameters : the coarse birnbaum function
    DOUBLE PRECISION :: chiBirn(kMaxPts),xBirn(kMaxPts)
    integer :: iNptsBirn
! output params
    DOUBLE PRECISION :: z(kMaxPtsBox)

! local variables
    DOUBLE PRECISION :: dvc,vbc,vtc,df
          
    integer ::           i,j,nl,nh,nstp,npts
    DOUBLE PRECISION ::  fchi,fdif

! we have 10 d of c parameters tau2=1e-3*(1,3,5,...,19) (ndoc=10)
! we have 7 temperatures : 100,150,200,250,300,350,400  (ntemp=7)
! the frequency spacing of the y's is 0.25
! the dnu's range from -250 to +250 in steps of 0.25    (npts=2001)
    data vbc,vtc,dvc/-250.0,+250.0,0.25/

    if (n_in > kMaxPtsBox) THEN
        write(kStdErr,*) 'in birnbaum.f n_in > kMaxPtsBox'
        CALL DoStop
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    nl=1
    nh=n_in
    nstp=1

    npts = iNptsBirn
    DO i=nl,nh,nstp
        fdif = (v(i) - v0)
        j = min((npts-2),int(fdif/dvc)) + 1
        j = j + 1000     !!need offset because vbc = -250; 250/0.25 = 1000
        fchi = vbc + (j - 1)*dvc
        df=(fdif-fchi)/dvc
    ! o the interpolation in freq
        z(i) = chiBirn(j) + (chiBirn(j+1) - chiBirn(j))*df
    !        print *,i,sngl(v(i)),sngl(z(i))
    END DO

    RETURN
    end subroutine birnbaum_interp

!************************************************************************
