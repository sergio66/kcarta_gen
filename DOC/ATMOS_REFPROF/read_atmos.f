c this program reads in atmos data files (from web) and produces new files
c for Matlab to read in
c compile with 
c f77 -o read_atmos.x -N109 -W -O3 read_atmos.f

      program read_atmos

      implicit none

      integer      iProf
      character*80 caFnameIN,caFout,caIDout
      real raaData(57,41)     !!there are 57 gases and 41 layers per gas
      character*20 caa20(57)  !!get the gasIDs and chemical names

      print *,'Enter profile number (1..12) '
      read *,iProf

      if (iProf .LT. 1 .OR. iProf .GT. 12) then
        print *,'invalid prof : should be between 1 -- 12'
        stop
        end if

      CALL getfilenames(iProf,caFnameIN,caFout,caIDout)
      print *,caFnameIN,caFout,caIDout

      CALL readdata(caFnameIN,raaData,caa20)
      CALL writedata(caFout,caIDout,raaData,caa20)
      end

c************************************************************************
      SUBROUTINE writedata(caFout,caIDout,raaData,caa20)

      implicit none

      character*80 caFout,caIDout
      real raaData(57,41)     !!there are 57 gases and 41 layers per gas
      character*20 caa20(57)  !!get the gasIDs and chemical names

c local vars
      INTEGER iG,iL,iIOUN,iErr

      iIOUN = 10
      OPEN(UNIT=iIOUN,FILE=caFout,STATUS='UNKNOWN',FORM='FORMATTED', 
     $    IOSTAT=iErr) 
      IF (iErr .NE. 0) THEN 
        WRITE(*,1080) iErr, caFout
        stop
        END IF

      DO iL = 1,41
        write(iIOUN,*) (raaData(iG,iL),iG=1,57)
        END DO
      CLOSE(iIOUN)

      iIOUN = 10
      OPEN(UNIT=iIOUN,FILE=caIDout,STATUS='UNKNOWN',FORM='FORMATTED', 
     $    IOSTAT=iErr) 
      IF (iErr .NE. 0) THEN 
        WRITE(*,1080) iErr, caFout
        stop
        END IF

      DO iG = 1,57
        write(iIOUN,20) caa20(iG)
        END DO
      CLOSE(iIOUN)

 20   FORMAT(A20)
 1080 FORMAT('ERROR! number ',I5,' opening profile file:',/,A80)       
      RETURN
      END

c************************************************************************
      SUBROUTINE parser(caStr,iG,caa20)

      implicit none

c input
      character*80 caStr
      integer iG
c output 
      character*20 caa20(57)  !!get the gasIDs and chemical names 
        
c local vars
      integer iI,iJ,iK
      character*20 caTemp

      iI = 1
      iJ = 0
      !! loop till we come to an integer
 10   CONTINUE
      if (caStr(iI:iI) .EQ. ' ') THEN
        iI = iI + 1
        GOTO 10
        END IF
      iJ = iI-1

      !!! now that we have come to an integer, loop till we leave the integer
 20   CONTINUE
      if (caStr(iI:iI) .NE. ' ') THEN
        iI = iI + 1
        GOTO 20
        END IF

      !! loop till we come to an char
 30   CONTINUE
      if (caStr(iI:iI) .EQ. ' ') THEN
        iI = iI + 1
        GOTO 30
        END IF

      !! loop till we leave the string
 40   CONTINUE
      if (caStr(iI:iI) .NE. ' ') THEN
        iI = iI + 1
        GOTO 40
        END IF
      
      DO iK = 1,20
        caTemp(iK:iK) = ' '
        END DO
      DO iK = iJ,iI
        caTemp(iK-iJ+1:iK-iJ+1) = caStr(iK:iK)
        END DO

      caa20(iG) = caTemp
      print *,caTemp

      RETURN
      END

c************************************************************************
      SUBROUTINE readdata(caFnameIN,raaData,caa20)

      implicit none

c input
      character*80 caFnameIN
c output 
      real raaData(57,41)     !!there are 57 gases and 41 layers per gas 
      character*20 caa20(57)  !!get the gasIDs and chemical names 

c local vars
      integer iG,iLay,iIOUN,iErr
      real raProf(41)
      character*80 caStr

      iIOUN = 10
      OPEN(UNIT=iIOUN,FILE=caFnameIN,STATUS='OLD',FORM='FORMATTED', 
     $    IOSTAT=iErr) 
      IF (iErr .NE. 0) THEN 
        WRITE(*,1080) iErr, caFnameIN 
        stop
        END IF

      DO iG = 1,57
        read(iIOUN,5020) caStr
        !!parse this to get the gas number and the chemical spieces
        CALL parser(caStr,iG,caa20)
        read(iIOUN,*) (raProf(iLay),iLay=1,41)        
        !turn from ppbv to ppmv
        DO iLay = 1,41
          raaData(iG,41-iLay+1) = raProf(iLay)*1e6
          END DO
        END DO
      CLOSE(iIOUN)

 5020 FORMAT(A80) 
 1080 FORMAT('ERROR! number ',I5,' opening profile file:',/,A80) 

      RETURN
      END

c************************************************************************
      SUBROUTINE getfilenames(iProf,caFnameIN,caFout,caIDout)

      implicit none

c input 
      integer iProf
c output
      character*80 caFnameIN,caFout,caIDout

c local vars
      character*2  ca2
      integer ii,i0,i1

      do ii = 1,80
        caFnameIN(ii:ii) = ' '
        caFout(ii:ii)    = ' '
        caIDout(ii:ii)   = ' '
        end do

      if (iProf .LT. 10) THEN
        caFnameIN = 'apriori_'
        caFout    = 'outfile_'
        caIDout   = 'idfile_'
        i0        = 0
        i1        = iProf       
      else
        caFnameIN = 'apriori_'
        caFout    = 'outfile_'
        caIDout   = 'idfile_'
        i0        = 1
        i1        = iProf-10       
        endif

      ca2 = CHAR(i0+48)//CHAR(i1+48)

      i1 = 80
 10   CONTINUE
      if ((caFnameIN(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) then
        i1 = i1 - 1
        GOTO 10
        END IF
      i1 = i1 + 1
      caFnameIN(i1:i1+1) = ca2
      i1 = i1 + 2
      caFnameIN(i1:i1+4) = '.dat'

 20   CONTINUE
      if ((caFout(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) then
        i1 = i1 - 1
        GOTO 20
        END IF
      i1 = i1 + 1
      caFout(i1:i1+1) = ca2
      i1 = i1 + 2
      caFout(i1:i1+4) = '.txt'

 30   CONTINUE
      if ((caIDout(i1:i1) .EQ. ' ') .AND. (i1 .GT. 1)) then
        i1 = i1 - 1
        GOTO 30
        END IF
      i1 = i1 + 1
      caIDout(i1:i1+1) = ca2
      i1 = i1 + 2
      caIDout(i1:i1+4) = '.txt'

      RETURN
      END
c************************************************************************
