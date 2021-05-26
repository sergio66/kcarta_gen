c cp ../../NONLTE2/sergio/VT_48PROFILES/convertQV.f .
c this is for 48 LTE --> NLTE profiles
c this function reads in 
c a) the kinetic profiles from KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/pt
c b) the VT      profiles from KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/
c PERTURBS the 2349 cm-1 band by +1 K above 30 km
c merges them as necessary and puts the output into
c                  KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/sergio_merge
c************************************************************************

c this file takes in Manuel Lopez-Puertas files, and adds info to them so that
c they are in the same format as files required for Dave Edwards GENLN2
c The files can then be read in by subroutine  ReadGENLN2_NLTE_Profile

c same as convert.f except we try to compute the QV partition function using
c equation 17 of the JQSRT paper by Edwards, Lopez-Puertas, Gamache
c (v 59, no 3, pg 426 , 1998)

c 2/4/21  Trying to get things going for Manuel's latest profiles + disk
C 9/10/05 IMPORTANT UPDATE : use vt_md_new.dat instead of vt_md.dat
c      this "asks" for  51 profiles instead of 47 profiles

c compile with f77 -W -N114 -f -C convertQV_perturb.f; 
c or /usr/bin/g77 -fcase-lower -ffixed-line-length-120 convertQV_perturb.f
c or ifort -u -extend-source 132 convertQV_perturb.f
c   /bin/mv a.out convertQV_perturb.x

      IMPLICIT NONE 

      INTEGER iProf,iSol    !!! profile number, solar angle

      print *,'Enter profile number, solar zenith angle : '
      read *,iProf,iSol
      CALL process(iProf,iSol)

      END
c************************************************************************
c this function reads in 
c a) the kinetic profiles from KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/pt
c b) the VT      profiles from KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/
c merges them as necessary and puts the output into
c                  KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/sergio_merge

c a) the kinetic profiles from ../../NONLTE2/sergio/VT_48PROFILES/merge/pt -- is this similar to PUERTAS2021?
c                              or just use KCARTA/SRC/NONLTE2/sergio/VT_48PROFILES/merge/pt
c b) the VT      profiles from ../../NONLTE2/sergio/VT_48PROFILES/merge -- is this similar to ../PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/vt_airs_v3_2
c merges them as necessary and puts the output into
c                <<<<<<<<<<<              ../PUERTAS_DATA/merge/   >>>>>>>>>>     

      SUBROUTINE process(iProf,iSol)

      IMPLICIT NONE 

c input params
      INTEGER iProf,iSol    !!! profile number, solar angle

c local vars
      CHARACTER*132 caFname1,caFname2      !!! files to read
      CHARACTER*132 caFname                !!! file to output
      CHARACTER*132 caTemp

      INTEGER iN                          !!!from kinetic profile
      REAL raH(200),raP(200),raT(200)     !!!from kinetic profile
      INTEGER iaDaveEdwards(100)          !!!vib centers in vt_md_new.dat
      INTEGER iCount
      INTEGER iNV                         !!! from LP profile
      INTEGER iI,iJ,iK

c for QV
      REAL raaaVT(10,70,200)                      !!!ISO,band,T
      REAL raaBand(10,70)                        !!!ISO,band
      REAL raaQV(10,200)                         !!!the QVibPartFcns
      INTEGER iaQISO(10)
      
      DO iI = 1,10
        iaQISO(iI) = 0
        DO iJ = 1,70
          raaBand(iI,iJ) = -10.0
          DO iK = 1,200
            raaaVT(iI,iJ,iK) = -10.0
            END DO
          END DO
        END DO
      DO iI = 1,10
        DO iK = 1,200
          raaQV(iI,iK) = 1.0           !!!start out with the gnd state E = 0
                                       !!! and so (exp(-E/T) = 1)
          END DO
        END DO

      DO iN = 1,100
        iaDaveEdwards(iN) = -1
        END DO

c      caFname1 = '/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES/'
c      caFname2 = '/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES/'
c      caFname  = '/home/sergio/KCARTA/NONLTE2/sergio/VT_48PROFILES/'
c      caTemp = 'merge/'              !this is where Manuel's NLTE profs are
c      CALL addpath(caFname1,caTemp)  
c      caTemp = 'merge/pt/'           !the 48 regression profiles are here
c      CALL addpath(caFname2,caTemp)
c      caTemp = 'sergio_merge/'       !this is where the output file will be
c      CALL addpath(caFname ,caTemp)

c /home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/forKCARTA/ is basically 
c /home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/vt_airs_v3_2/
      caFname1 = '/home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/forKCARTA/'
      caFname2 = '/home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/forKCARTA/'
      caFname  = '/home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/forKCARTA_PerturbV1/'
      caTemp = 'merge/'              !this is where Manuel's NLTE profs are
      CALL addpath(caFname1,caTemp)  
      caTemp = 'merge/pt/'           !the 48 regression profiles are here
      CALL addpath(caFname2,caTemp)
      caTemp = 'sergio_merge/'       !this is where the output file will be
      CALL addpath(caFname ,caTemp)

      CALL addinfo(caFname1,iProf,iSol)
      CALL addinfo(caFname2,iProf,-1)
      CALL addinfo(caFname,iProf,iSol)

      print *,'Name of Manuel Lopez-Puertas input file = ',caFName1
      print *,'Name of file for kinetic temperature = ',caFName2
      print *,'Name of output file = ',caFName
      
      CALL ReadKineticProf(caFname2,raH,raP,raT,iN)
      CALL VibTemp(raH,raP,raT,iN,iProf,iSol,caFname1,caFName2,caFname,
     $                 iaDaveEdwards,-1,iNV,raaaVT,raaBand,iaQISO,raaQV,iCount)
      CALL ComputeQV(iNV,raaaVT,raaBand,iaQISO,raaQV)
      CALL VibTemp(raH,raP,raT,iN,iProf,iSol,caFname1,caFName2,caFname,
     $                 iaDaveEdwards,+1,iNV,raaaVT,raaBand,iaQISO,raaQV,iCount)
      RETURN
      END

************************************************************************
c this subroutine computes the QV function according to Edwards, Puertas, 
c Gamache eqn
      SUBROUTINE ComputeQV(iNV,raaaVT,raaBand,iaQISO,raaQV)

      IMPLICIT NONE

c input vars
      INTEGER iaQISO(10),iNV
      REAL raaaVT(10,70,200)                      !!!ISO,band,T
      REAL raaBand(10,70)                        !!!ISO,band
c output vars
      REAL raaQV(10,200)                         !!!the QVibPartFcns

c local vars
      INTEGER iI,iJ,iK
      REAL c2,rX,iD
      INTEGER iaDegenBand(10),iaaDegen_HIT_ID(10,51),iaaDegen_number(10,51)
      REAL raaDegen_Bandcenter(10,51)
      INTEGER find_degeneracy

      c2 = 1.4387863          !!!second radiative constant; see c1,c2 in kCARTA

      CALL load_co2_file(
     $  iaDegenBand,iaaDegen_HIT_ID,raaDegen_BandCenter,iaaDegen_number)

      DO iI = 1,10
        IF (iaQISO(iI) .GT. 0) THEN
          print *,' isotope, num of bands = ',iI,iaQISO(iI)
          DO iJ = 1,iaQISO(iI)
            iD = 1
            iD = find_degeneracy(iI,-1,raaBand(iI,iJ),
     $  iaDegenBand,iaaDegen_HIT_ID,raaDegen_BandCenter,iaaDegen_number)
            DO iK = 1,iNV
              rX = c2*raaBand(iI,iJ)/raaaVT(iI,iJ,iK)
              raaQV(iI,iK) = raaQV(iI,iK) + iD*exp(-rX)
              END DO 
            END DO
          END IF
        END DO

      RETURN
      END

c************************************************************************
c this function computes the degeneracy
      INTEGER FUNCTION find_degeneracy(iISO,iHIT_ID,rX,
     $  iaDegenBand,iaaDegen_HIT_ID,raaDegen_BandCenter,iaaDegen_number)

      IMPLICIT NONE

c input vars
      INTEGER iaDegenBand(10),iaaDegen_HIT_ID(10,51),iaaDegen_number(10,51)
      REAL raaDegen_Bandcenter(10,51),rX
      INTEGER iISO,iHIT_ID

c local vars 
      INTEGER iI,iJ,iAns,iFound

      iAns = 1    !default degen

      IF  ((iHIT_ID .GT. 0) .AND. (rX  .LT. 0)) THEN
        !!need to search on HITRAN iDs
        iFound = -1
        iJ = 1
 40     CONTINUE
        IF (iaaDegen_HIT_ID(iISO,iJ) .EQ. iHIT_ID) THEN
          iAns = iaaDegen_number(iISO,iJ)
          iFound = +1
          print *,'set     degen',iISO,iaaDegen_HIT_ID(iISO,iJ),
     $            raaDegen_BandCenter(iISO,iJ),iAns
          GOTO 50
        ELSEIF (iJ .LT. iaDegenBand(iISO)) THEN
          iJ = iJ + 1
          GOTO 40
          END IF
        END IF

      IF  ((iHIT_ID .LT. 0) .AND. (rX  .GT. 0)) THEN
        !!need to search on band centers
        iFound = -1
        iJ = 1
 45     CONTINUE
        IF (abs(raaDegen_Bandcenter(iISO,iJ)-rX) .LE. 0.1) THEN
          iAns = iaaDegen_number(iISO,iJ)
          iFound = +1
          print *,'set     degen',iISO,iaaDegen_HIT_ID(iISO,iJ),
     $            raaDegen_BandCenter(iISO,iJ),iAns
          GOTO 50
        ELSEIF (iJ .LT. iaDegenBand(iISO)) THEN
          iJ = iJ + 1
          GOTO 45
          END IF
        END IF

 50   CONTINUE

      IF (iFound .LT. 0) THEN
        print *,'default degen',iISO,iHIT_ID,rX,1
        END IF
      find_degeneracy = iAns

      RETURN
      END

c************************************************************************
c this function loads in the degeneracy file
      SUBROUTINE load_co2_file(
     $  iaDegenBand,iaaDegen_HIT_ID,raaDegen_BandCenter,iaaDegen_number)

      IMPLICIT NONE

c output vars
      INTEGER iaDegenBand(10),iaaDegen_HIT_ID(10,51),iaaDegen_number(10,51)
      REAL raaDegen_Bandcenter(10,51)

c local vars
      INTEGER iI,iJ,iK,iL,iIOUN,iCount,iISO,iHITRANID,iDegeneracy,iErr
      REAL rBandCenter
      CHARACTER*132 cafname,caStr

      DO iI = 1,10
        iaDegenBand(iI) = 0
        END DO
      DO iI = 1,10
        DO iJ = 1,51
          iaaDegen_HIT_ID(iI,iJ) = 0
          iaaDegen_number(iI,iJ) = 0
          raaDegen_BandCenter(iI,iJ) = 1.0e10
          END DO
        END DO

      caFname = 'co2_degeneracies.txt'

      iIOun = 11
      OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr)  
      IF (iErr .NE. 0) THEN  
        print *, 'error reading co2-degeneracies ',iErr,caFName
        stop
        END IF  
      
 40   CONTINUE
      READ (iIOUN,10,end=50) caStr
      IF (caStr(1:1) .EQ. '!') THEN
        GOTO 40
      ELSE
        READ(caStr,*)  iCount,iISO,iHITRANID,rBandCenter,iDegeneracy
!        print *,iCount,iISO,iHITRANID,rBandCenter,iDegeneracy
        iaDegenBand(iISO) = iaDegenBand(iISO) + 1
        iaaDegen_HIT_ID(iISO,iaDegenBand(iISO)) = iHITRANID
        raaDegen_BandCenter(iISO,iaDegenBand(iISO)) = rBandCenter
        iaaDegen_number(iISO,iaDegenBand(iISO)) = iDegeneracy
        GOTO 40
        END IF

 50   CONTINUE
      CLOSE(iIOUN)

 10   FORMAT(A80)

      print *,'read degeneracy table '
      DO iI = 1,10
        IF (iaDegenBand(ii) .GT. 0) THEN 
          print *,iI,iaDegenBand(ii)
          END IF
        END DO
      print *,' '

      RETURN
      END

c************************************************************************
c this subroutine reads in the raH,raP,raT from the kinetic profile
      SUBROUTINE ReadKineticProf(caFname2,raH,raP,raT,iN)

      IMPLICIT NONE

      CHARACTER*132 caFname2               !input kinetic profile file name
      REAL raH(200),raP(200),raT(200)     !outputs from kinetic profile
      INTEGER iN

c local vars 
      INTEGER iIOUN,iErr,iI
      CHARACTER*132 caStr

      iIOun = 11
      OPEN(UNIT=iIOun,FILE=caFName2,FORM='formatted',STATUS='OLD', 
     $     IOSTAT=iErr)  
      IF (iErr .NE. 0) THEN  
        print *, 'error reading kinetic profile ',iErr,caFName2
        stop
        END IF  

      read (iIOUN,10) caStr             !Mixer generated PT profiles
      read (iIOUN,10) caStr             !number of levels
      read (iIOUN,10) caStr             !dollar
      read (iIOUN,*)  iN                !number of levels

      read (iIOUN,10) caStr             !empty
      read (iIOUN,10) caStr             !altitude(km)
      read (iIOUN,10) caStr             !dollar
      read (iIOUN,*) (raH(iI),iI=1,iN)

      read (iIOUN,10) caStr             !empty
      read (iIOUN,10) caStr             !pressure (mb)
      read (iIOUN,10) caStr             !filename
      read (iIOUN,*) (raP(iI),iI=1,iN)

      read (iIOUN,10) caStr             !empty
      read (iIOUN,10) caStr             !kinetic temp (K)
      read (iIOUN,10) caStr             !filename
      read (iIOUN,*) (raT(iI),iI=1,iN)

      CLOSE(iIOUN)
 10   FORMAT(A80)

      print *,'read kinetic profile '
      RETURN
      END
c************************************************************************
c this subroutine reads in the vib temps, and outputs the kinetic+vibtemps
c in the format that GENLN2 and kCARTA like
      SUBROUTINE VibTemp(raH,raP,raT,iN,iProf,iSol,caFname1,caFName2,caFname,
     $                   iaDaveEdwards,iPrint,
     $                   iNV,raaaVT,raaBand,iaQISO,raaQV,iCount)

      IMPLICIT NONE

c input vars
      CHARACTER*132 caFname1               !input VT profile file name
      CHARACTER*132 caFname2               !input kinetic profile file name
      REAL raH(200),raP(200),raT(200)     !outputs from kinetic profile
      INTEGER iN,iProf,iSol
      INTEGER iPrint,iaDaveEdwards(100)   !print stuff to caFName?(N/Y = -1/+1)
                                          !how much vt_md stuff to put

c input/output vars for the vib part fcns
      REAL raaaVT(10,70,200)                      !!!ISO,band,T
      REAL raaBand(10,70)                        !!!ISO,band
      REAL raaQV(10,200)                         !!!the QVibPartFcns
      INTEGER iaQISO(10)
      INTEGER iNV                                !!!how many levels LP has

c output vars
      CHARACTER*132 caFname                !output profile file name
      INTEGER iCount                      !how many vib temp profiles to put in

c local vars 
      INTEGER iIOUN1,iIOUN,iErr,iI,iJ,iZero,iLocalIso
c used to read in LP's files
      INTEGER iNumPuertas,iLevel,iGasID_ISO,iAFGL
      CHARACTER*1 cDollar
      REAL rVibCenter,raVibTemp(200)
      CHARACTER*132 caStr
c interpolate the kinetic profile info, to the vibrational profile levels
      REAL raH2(200)            !heights used in the VT model
      REAL raP2(200),raT2(200)  !interp1(raH,raP,raH2,raP2)
                                !interp1(raH,raT,raH2,raT2)
      REAL raVibPart(200)       !dummy set to 1.0 for now
c used to read HITRAN/AFGL vib level info from vt_md.dat
      INTEGER iaIDAFGL(100),iaISO(100),iaIDISO(100),iaLevel(100),iaGasID(100)
      REAL raVib(100)

      DO iI = 1,200
        raVibPart(iI) = 1.0
        END DO

      IF (iPrint .GT. 0) THEN
        iIOun = 12
        OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='NEW', 
     $     IOSTAT=iErr)  
        IF (iErr .NE. 0) THEN  
          print *, 'error opening output profile ',iErr,caFName
          stop
          END IF  

        caStr = '! vibrational temperature profiles reformatted'
        write(iIOUN,13) caStr
        caStr = '! code : /home/sergio/KCARTA/NONLTE3_Feb2021/MAKE_vtPROF_sSOL_prf/convertQV.f'
        write(iIOUN,13) caStr
        caStr = '! Sergio Machado/Scott Hannon give levels (P,T) info'
        write(iIOUN,13) caStr
        caStr = '! for regression profile number ----->>> '
        write(iIOUN,12) caStr,iProf
        caStr = '! upto typical height 80 km '
        write(iIOUN,13) caStr
        caStr = '! Manuel Lopez-Puertas tacks on info to gives kinetic '
        write(iIOUN,13) caStr
        caSTr = '! profiles upto 200 km. Additionally, he computes '
        write(iIOUN,13) caStr
        caStr = '! NLTE Vib Temps for SZA ----->>>'
        write(iIOUN,12) caStr,iSol
        caStr = '! Files merged are (a) kinetic temperature = '
        write(iIOUN,13) caStr
        write(iIOUN,11) caFName2
        caStr = '!                  (b) vibrational temperature = '
        write(iIOUN,13) caStr
        write(iIOUN,11) caFName1
        END IF

      iIOun1 = 11
      OPEN(UNIT=iIOun1,FILE=caFName1,FORM='formatted',STATUS='OLD', 
     $     IOSTAT=iErr)  
      IF (iErr .NE. 0) THEN  
        print *, 'error reading vib profile ',iErr,caFName1
        stop
        END IF  

      read (iIOUN1,10) caStr             !dollar
      read (iIOUN1,*) iNV                !number of vib temp levels
      read (iIOUN1,10) caStr             !empty
      read (iIOUN1,10) caStr             !dollar
      read (iIOUN1,*) (raH2(iI),iI=1,iNV)!read the heights (km)
      CALL  vibstate(iaGasID,iaISO,iaIDISO,iaLevel,iaIDAFGL,raVib)

      IF (iPrint .GT. 0) THEN
        CALL rspl(raH,raT,iN,raH2,raT2,iNV)
        CALL rspl(raH,raP,iN,raH2,raP2,iNV)

        caStr = '! The vib temps are given for '
        write(iIOUN,13) caStr
        caStr = '! the following number of CO2 bands : iCount =  '
        write(iIOUN,12) caStr,iCount
        caStr = '! vib. states at following number of pressure levels : iNV = '
        write(iIOUN,12) caStr,iNV

        iJ = 0
        caStr = 
     $   '!   IL    MOL   IDMOL  ISO   IDISO  LEVEL   IDAFGL  ENERGY(cm-1)'
        write(iIOUN,13) caStr
        DO iI = 1,51
          IF (iaDaveEdwards(iI) .GT. 0) THEN
            iJ = iJ + 1
            write(iIOUN,20) iJ,iaGasID(iI),iaISO(iI),iaIDISO(iI),iaLevel(iI),
     $                  iaIDAFGL(iI),raVib(iI)
            END IF
          END DO

        iZero = 0
        caStr = '! '
        write(iIOUN,13) caStr
        caStr = '***'
        !!!dave edwards calls this FLAG,MODEL,NSET,NGAS,NLEV 
        !!!where FLAG  = *** 
        !!!      MODEL = atmosphere model number for data file 
        !!!      NSET  = NO OF VIBRATIONAL STATE DATA SETS FOR THIS MODEL 
        !!!      NGAS  = NO OF DIFFERENT GASES FOR WHICH THERE ARE T_VIB      
        !!!      NLEV  = NO OF PRESSURE LEVELS FOR THIS MODEL  
        !!! so read in flag, model,nset,ngas,nlev where flag = '***' 
        write(iIOUN,99)    caStr,1, iCount,1,iNV   
                                          ! 1,iNV gives # of NLTE gases and
                                          ! the number of levels
        write(iIOUN,*) 2,6         ! how many vib part fcns for the NLTE gas
        CALL printarray(iIOUN,iNV,raP2)  !!!print out press levels
        write(iIOUN,*) iZero       !dummy
        CALL printarray(iIOUN,iNV,raT2)  !!!print out kinetic temps

        caStr = '!QV   Vibrational Partition functions '
        write(iIOUN,13) caStr
        DO iI = 1,6
          DO iJ = 1,200
            raVibPart(iJ) =  raaQV(iI,iJ)
            END DO
          write(iIOUN,*) 2,iI
          CALL printarray(iIOUN,iNV,raVibPart)  !!!print out part fcn
          END DO
        END IF

      read (iIOUN1,10) caStr             !empty
      read (iIOUN1,10) caStr             !dollar
      read (iIOUN1,10) caStr             !F
      read (iIOUN1,10) caStr             !empty
      read (iIOUN1,10) caStr             !dollar
      read (iIOUN1,*)  iNumPuertas       !number of calcs done by Lopez-Puertas

      IF (iPrint .GT. 0) THEN
        caStr = '!TV   Vibrational Temperatures'
        write(iIOUN,13) caStr
        END IF
      iCount = 0

      DO iI = 1,iNumPuertas
        read (iIOUN1,10) caStr             !empty
        read(iIOUN1,*) cDollar,iLevel,rVibCenter
        read(iIOUN1,*) iGasID_ISO
        read(iIOUN1,*) iAFGL
        read(iIOUN1,*) (raVibTemp(iJ),iJ=1,iNV)

c see eg /home/sergio/KCARTA/NONLTE3_Feb2021/PUERTAS_DATA/AIRS_CO2_VTs_v3_ig2_Dusk_Mar2020/forKCARTA/merge/vt1_s0.prf
c heights are at 1 km intervals
c$    00011   2349.1433
c 21   
c 9    
        if ((iLevel .EQ. 00011) .and. (abs(rVibCenter-2349.1433) .LT. 1.0)) then
          DO iJ = 30,39
            raVibTemp(iJ) = raVibTemp(iJ) + (iJ*1.0-30.0)/10.0
          END DO
          DO iJ = 40,iNV
            raVibTemp(iJ) = raVibTemp(iJ) + 1.0
          END DO
        end if

        !!! store the info in matrices for QV computation
        iLocalIso = iGasID_ISO - 20   !!eg 21 == isotope 1 etc
        iaQISO(iLocalIso) = iaQISO(iLocalIso) + 1        
        raaBand(iLocalIso,iaQISO(iLocalIso)) = rVibCenter
        DO iJ = 1,iNV
          raaaVT(iLocalIso,iaQISO(iLocalIso),iJ) = raVibTemp(iJ)
          END DO
        !!! see if this band needs to be output
        CALL IdentifyPrint(iaGasID,iaISO,iaIDISO,iaLevel,iaIDAFGL,raVib,
     $                     iLevel,rVibCenter,iGasID_ISO,iAFGL,raVibTemp,
     $                     iCount,iI,iIOUN,iNV,iPrint,iaDaveEdwards)
        END DO
      print *,'matched up ',iCount,' out of ',iNumPuertas,'profiles '

      CLOSE(iIOUN1)
      IF (iPrint .GT. 0) THEN
        CLOSE(iIOUN)
        END IF

 10   FORMAT(A80)
 11   FORMAT('! ',A70)
 12   FORMAT(A70,I5)
 13   FORMAT(A70)
 20   FORMAT('! ',I4,' CO2   ',I5,'  ',I5,'  ',I5,'  ',I7,' ',I6,'   ',F11.5)
c 20   FORMAT('! ',I4,'   ',I5,'  ',I7,' ',I6,'   ',F11.5)
 99   FORMAT(A3,'   ',4(I3,' '))
      RETURN
      END
c************************************************************************
       SUBROUTINE rspl(XA,YA,N,XOUT,YOUT,NOUT) 
 
      IMPLICIT NONE

c this subroutine directly calls rsply2 and then rspline

C real version 
C      ----------------------------------------------------------------- 
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y 
C      XA  : I  : DOUB arr : x array(N) in increasing order  IN 
C      YA  : I  : DOUB arr : y array(N)                      IN 
C      N   : I  : INT      : number of points in arrays 
C      XOUT   : I  : DOUB ARR     : x points at which to evaluate spline 
C      YOUT   : O  : DOUB ARR     : y points from spline interpolation 
C      NOUT   : I  : INT          : number of points at which to spline 
C      ----------------------------------------------------------------- 
C 
C      Parameters 
       REAL XA(*),YA(*),XOUT(*),YOUT(*) 
       INTEGER N,NOUT 
 
       REAL Y2A(200),work(200),Yp1,yPn
       INTEGER I 

       yp1=1.0e16
       ypn=1.0e16

       CALL rSPLY2(XA,YA,N,Yp1,Ypn,Y2A,work) 
       DO I=1,NOUT 
         CALL rSPLIN(XA,YA,Y2A,N,XOUT(I),YOUT(I))  
         END DO 
 
       RETURN 
       END 
 
c************************************************************************ 
       SUBROUTINE rSPLIN(XA,YA,Y2A,N,X,Y)

      IMPLICIT NONE


C REAL version
C      -----------------------------------------------------------------
C      Uses Y2A from SPLY2 to do spline interpolation at X to get Y
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      Y2A : I  : DOUB arr : 2nd derivative of points
C      N   : I  : INT      : number of points in arrays
C      X   : I  : DOUB     : x point at which to evaluate spline
C      Y   : O  : DOUB     : y point from spline interpolation
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),Y2A(*),X,Y
       INTEGER N
C
C      Local Variables
       INTEGER K,KLO,KHI
       REAL A,B,H
C
C      -----------------------------------------------------------------
C
C      Determine between which pair of pints X falls (bisect loop)
       KLO=1
       KHI=N
 20    IF ( (KHI - KLO) .GT. 1) THEN
          K=(KHI + KLO)/2
          IF (XA(K) .GT. X) THEN
             KHI=K
          ELSE
             KLO=K
          ENDIF
          GOTO 20
       ENDIF
C
       H=XA(KHI) - XA(KLO)
       IF (H .LE. 0.0) THEN
          write(*,1010) KLO,KHI,XA(KLO),XA(KHI)
 1010     FORMAT('ERROR! rSPLINT: bad XA array.',/,
     $       'KLO=',I5,', KHI=',I5,', XA(KLO)=',1PE12.5,
     $       ', XA(KHI)=',E12.5,'. Quitting.')
       ENDIF
C
       A=(XA(KHI) - X)/H
       B=(X - XA(KLO))/H
C
       Y=A*YA(KLO) + B*YA(KHI) + ( Y2A(KLO)*(A**3 - A) +
     $    Y2A(KHI)*(B**3 - B) )*(H**2)/6.0
C
       RETURN
       END

c************************************************************************
       SUBROUTINE rSPLY2(XA,YA,N,YP1,YPN,Y2A,WORK)

      IMPLICIT NONE

C
C REAL version
C      -----------------------------------------------------------------
C      Calc 2nd derivative as preperation for SPLINT spline routine
C      XA  : I  : DOUB arr : x array(N) in increasing order
C      YA  : I  : DOUB arr : y array(N)
C      N   : I  : INT      : number of points in arrays
C      YP1 : I  : DOUB     : derivative of 1st point
C      YPN : I  : DOUB     : derivative of last point
C      Y2A : O  : DOUB arr : 2nd derivative array(N)
C      WORK: O  : DOUB arr : workspace array(N)
C      -----------------------------------------------------------------
C
C      Parameters
       REAL XA(*),YA(*),Y2A(*),YP1,YPN,WORK(*)
       INTEGER N
C
C      Local Variables
       INTEGER I,K
       REAL P,QN,SIG,UN
C
C      -----------------------------------------------------------------
C
C      Lower boundary
       IF (YP1 .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          Y2A(1)=0.0
          WORK(1)=0.0
       ELSE
C         Set to a specific first derivative
          Y2A(1)=-0.5
          WORK(1)=( 3.0/(XA(2) - XA(1)) )*( (YA(2) - YA(1))/
     $       (XA(2) - XA(1)) - YP1)
       ENDIF
C
C      Decomposition loop of the tridiagonal algorithm
       DO I=2,N-1
c this is from the progas code
c          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I))
          SIG=(XA(I) - XA(I-1))/(XA(I+1) - XA(I-1))
          P=SIG*Y2A(I-1) + 2.0
          Y2A(I)=(SIG - 1.0)/P
          WORK(I)=(YA(I+1) - YA(I))/(XA(I+1) - XA(I)) -
     $       (YA(I) - YA(I-1))/(XA(I) - XA(I-1))
          WORK(I)=( 6.0*WORK(I)/(XA(I+1) - XA(I-1)) -
     $       SIG*WORK(I-1) )/P
       ENDDO
C
C      Upper boundary
       IF (YPN .GT. 1.0E+15) THEN
C         "Natural" boundary condition
          QN=0.0
          UN=0.0
       ELSE
C         Set to a specific first derivative
          QN=0.5
          UN=( 3.0/(XA(N) - XA(N-1)) )*( YPN -
     $       (YA(N) - YA(N-1))/(XA(N) - XA(N-1)) )
       ENDIF
       Y2A(N)=(UN - QN*WORK(N-1))/(QN*Y2A(N-1) + 1.0)
C
C      Assign the other 2nd derivatives using the back-substitution
C      loop of the tridiagonal algorithm
       DO K=N-1,1,-1
          Y2A(K)=Y2A(K)*Y2A(K+1) + WORK(K)
       ENDDO
C
       RETURN
       END

c************************************************************************
c this function takes string s1 and appends string s2 to it
      SUBROUTINE addpath(caS1,caS2)

      IMPLICIT NONE

      CHARACTER*132 caS1,caS2

      INTEGER i1,i2

      i1 = 132
 10   CONTINUE
      IF ((caS1(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) THEN
        i1 = i1 - 1
        GOTO 10
        END IF

      i2 = 132
 20   CONTINUE
      IF ((caS2(i2:i2) .EQ. ' ') .AND. (i2 .GT. 1)) THEN
        i2 = i2 - 1
        GOTO 20
        END IF

      caS1(i1+1:i1+i2) = caS2(1:i2)

      RETURN
      END

c************************************************************************
c this function takes string s1 and adds info to it, based on iProf,iSol
      SUBROUTINE addinfo(caS1,iProf,iSol)

      IMPLICIT NONE

      CHARACTER*132 caS1
      INTEGER iProf,iSol

      INTEGER i1
      CHARACTER*3 caString3
      CHARACTER*2 caString2

      i1 = 132
 10   CONTINUE
      IF ((caS1(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) THEN
        i1 = i1 - 1
        GOTO 10
        END IF

      IF (iSol .GE. 0) THEN
        !!we are dealing with vib temp files, so we need "vt" names
        caS1(i1+1:i1+2) = 'vt'
        i1 = i1 + 2
      ELSE
        !!we are dealing with profile files, so we need "pt" names
        caS1(i1+1:i1+2) = 'pt'
        i1 = i1 + 2
        END IF

c now process iProf so that we end up with a right padded string 
c eg 2 ---> '2 ', 12 ---> '12' etc 
      WRITE(caString2,15) iProf
      IF (iProf .GE. 10) THEN
        caS1(i1+1:i1+2) = caString2(1:2)
        i1 = i1 + 2
      ELSE
        caS1(i1+1:i1+1) = caString2(2:2)
        i1 = i1 + 1
        END IF

      !!we are dealing with vib temp files, so we need to put in sol angle
      IF (iSol .GE. 0) THEN
        caS1(i1+1:i1+2) = '_s'
        i1 = i1 + 2
        WRITE(caString3,25) iSol
        IF (iSol .LT. 10) THEN 
          caS1(i1+1:i1+1) = caString3(3:3)
          i1 = i1 + 1
        ELSEIF (iSol .LT. 100) THEN
          caS1(i1+1:i1+2) = caString3(2:3)
          i1 = i1 + 2
        ELSEIF (iSol .LT. 1000) THEN
          caS1(i1+1:i1+3) = caString3(1:3)
          i1 = i1 + 3
          END IF
        END IF

      caS1(i1+1:i1+4) = '.prf'

 15   FORMAT(I2) 
 25   FORMAT(I3) 

      RETURN
      END
c************************************************************************
c this subroutine inputs/outputs vib state info
      SUBROUTINE vibstate(iaGasID,iaISO,iaIDISO,iaLevel,iaIDAFGL,raVib)

      INTEGER iaIDAFGL(100),iaISO(100),iaIDISO(100),iaLevel(100),iaGasID(100)
      REAL raVib(100)

c local vars 
      INTEGER iIOUN,iI,iL,iErr
      CHARACTER*132 caFName,caStr
      CHARACTER    cDollar
      CHARACTER*3  caCO2

!      caFName = 
!     $  '/home/sergio/KCARTA/SRCv1.16/NONLTE/sergio/VT_ScriptProfs/vt_md_new.dat'
      caFName = 
     $  '/home/sergio/KCARTA/NONLTE_PRODUCTION/VT_ScriptProfs/vt_md_new.dat'

      iIOun = 14
      OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr)  
      IF (iErr .NE. 0) THEN  
        print *, 'error reading vib profile ',iErr,caFName
        stop
        END IF  

      DO iI = 1,11
        !get and ignore first 11 lines of comments
        READ(iIOUN,10) caStr
        END DO
      
      DO iI = 1,51
        READ(iIOUN,*) cDollar,iL,caCO2,iaGasID(iI),iaISO(iI),iaIDISO(iI),
     $                iaLevel(iI),iaIDAFGL(iI),raVib(iI)
        END DO

 10   FORMAT(A80)

      CLOSE(iIOUN)
      RETURN
      END
c************************************************************************
c this subroutine prints out the arrays in a nice format
      SUBROUTINE printarray(iIOUN,iNV,raX)

      INTEGER iIOUN,iNV,iimod,iDiv
      REAL raX(200)

c local vars 
      INTEGER iM,iD,i1,i2,ia(5),iJ,iFive
      REAL ra(5)

      iFive = 5
      iM = iimod(iNv,iFive)
      iD = iDiv(iNv,iFive)

      ia(1) = -4
      ia(2) = -3
      ia(3) = -2
      ia(4) = -1
      ia(5) =  0
      DO i1 = 1,iD
        DO i2 = 1,5
          ia(i2) = ia(i2) + 5
          ra(i2) = raX(ia(i2))
          END DO
        write(iIOUN,15) (ra(iJ),iJ=1,5)
        END DO

      IF ((iM .GT. 0) .AND. (iM .LT. 5)) THEN
        IF (iM .EQ. 1) THEN
          ra(1) = raX(iNV)
          write(iIOUN,11) (ra(iJ),iJ=1,iM)
        ELSE IF (iM .EQ. 2) THEN
          ra(2) = raX(iNV)
          ra(1) = raX(iNV-1)
          write(iIOUN,12) (ra(iJ),iJ=1,iM)
        ELSE IF (iM .EQ. 3) THEN
          ra(3) = raX(iNV)
          ra(2) = raX(iNV-1)
          ra(1) = raX(iNV-2)
          write(iIOUN,13) (ra(iJ),iJ=1,iM)
        ELSE IF (iM .EQ. 4) THEN
          ra(4) = raX(iNV)
          ra(3) = raX(iNV-1)
          ra(2) = raX(iNV-2)
          ra(1) = raX(iNV-3)
          write(iIOUN,14) (ra(iJ),iJ=1,iM)
          END IF
        END IF

 11   FORMAT(5(E11.5,' '))
 12   FORMAT(5(E11.5,' '))
 13   FORMAT(5(E11.5,' '))
 14   FORMAT(5(E11.5,' '))
 15   FORMAT(5(E11.5,' '))

      RETURN
      END
c************************************************************************
c this function does the integer MOD eg 100 div 3 = 1, 10 div 2 = 0
c assuming i1,i2 > 0 
      INTEGER FUNCTION iimod(i1,i2) 
 
      IMPLICIT NONE 
 
      INTEGER i1,i2,iInt
 
c recall INT truncates the real number 
      iInt = INT((i1*1.0)/(i2*1.0)) 
      iimod = i1 - iInt*i2

      RETURN 
      END 
 
c************************************************************************
c this function does the integer DIVISION eg 100 div 3 = 33, 10 div 2 = 5 
c assuming i1,i2 > 0 
      INTEGER FUNCTION idiv(i1,i2) 
 
      IMPLICIT NONE 
 
      INTEGER i1,i2,iInt 
 
c recall INT truncates the real number 
      iInt = INT((i1*1.0)/(i2*1.0)) 
      idiv = iInt  
      RETURN 
      END 
 
c************************************************************************
c this subroutine goes thru the HITRAN/AFGL id's read from vt_md_new.dat,
c compares the current iAFGL,iLevel info from Lopez-Puertas; if there is
c a match, print it!
      SUBROUTINE IdentifyPrint(iaGasID,iaISO,iaIDISO,iaLevel,iaIDAFGL,raVib,
     $              iLevel,rVibCenter,iGasID_ISO,iAFGL,raVibTemp,
     $              iCount,iIn,iIOUN,iNV,iPrint,iaDaveEdwards)

c input params from Dave Edwards md_vt.dat
      INTEGER iaIDAFGL(100),iaISO(100),iaIDISO(100),iaLevel(100),iaGasID(100)
      REAL raVib(100)
      INTEGER iPrint,iaDaveEdwards(100),iCOunt
c input params from Lopez-Puertas files
      INTEGER iLevel,iGasID_ISO,iAFGL,iNV,iIN
      REAL rVibCenter,raVibTemp(200)
      INTEGER iIOUN

c local vars
      INTEGER iFound,i1,iGasID,iISO

      iGasID = 2
      iISO = iGasID_ISO - 20

c        read(iIOUN1,*) cDollar,iLevel,rVibCenter
c        read(iIOUN1,*) iGasID_ISO
c        read(iIOUN1,*) iAFGL

      iFound = -1
      i1 = 1
 10   CONTINUE
      IF (iaIDISO(i1) .EQ. iISO) THEN        !!!isotopes match
        IF (iaLevel(i1) .EQ. iLevel) THEN    !!!HITRAN upper levels match
          IF (iaIDAFGL(i1) .EQ. iAFGL) THEN  !!!AFGL numbers match
c            IF (abs(raVib(i1) - rVibCenter) .LE. 0.1) THEN !!!vib cntr match
            IF (abs(raVib(i1) - rVibCenter) .LE. 25.0) THEN !!!vib cntr match
              iFound = +1
              iaDaveEdwards(i1) = +1
              END IF
            END IF
          END IF
        END IF
      IF ((iFound .LT. 0) .AND. (i1 .LT. 51)) THEN
        i1 = i1 + 1   !!! try again
        GOTO 10
        END IF

      IF (iFound .GT. 0) THEN
        print *,'Match for    ',iLevel,rVibCenter,iAFGL,iGasID,iISO,' >>> ',i1
        iCOunt = iCount + 1
        IF (iPrint .GT. 0) THEN
          write(iIOUN,20) iCount,iGasID,iISO,iAFGL,rVibCenter
          CALL printarray(iIOUN,iNV,raVibTemp)
          END IF
      ELSE
        print *,'NO match for ',iLevel,rVibCenter,iAFGL,iGasID,iISO,' <<< ',i1
        END IF

 20   FORMAT(4(I3,'  '),F15.6)
      RETURN 
      END 
 
c************************************************************************
