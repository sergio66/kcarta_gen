c this file takes in a GENLN2 template file, and depending on the values 
c of the 2 variables read in, produces a copy that can be used by atmos.x
c it also produces an "inputfile" for atmos.x
c NOTE : IT IS ASSUMED THIS PROGRAM RUN IS SAME DIRECTORY AS ATMOS.X
c        FOR TESTING PURPOSES!

      IMPLICIT NONE

      INTEGER InfileIndex,iFreqIndex
      REAL r1,r2,rSurfaceTemp

      CHARACTER*80 caSubDir
      REAL raTLow(100)
      CHARACTER*130 caStr130
      INTEGER iIOIN,ii

c first get the name of the profile to be read in, to find the temperature
      print *,' Enter profile index of input file (1-48) : '
      read (5,*) InfileIndex
c then get the freq range (1-6)
      print *,'Enter frequency index 1-6 : '
      read (5,*) iFreqIndex
c get the name of the subdirectory that data files will go to
      print *,'Enter subdirectory to save genln2/radiance files '  
      print *,'    (w/o the begin/end slashes !! eg TEST, not /TEST/)'
      read (5,111) caSubDir
 111  FORMAT(A80)
      CALL rightpad(caSubDir)

      CALL GetSurfaceTemp(InfileIndex,rSurfaceTemp)

      CALL GetFreqRange(r1,r2,iFreqIndex)
      CALL Genln2(r1,r2,InfileIndex,iFreqIndex,rSurfaceTemp,caSubDir)
      
      STOP
      END
c************************************************************************
c this subroutine gets the freq range
      SUBROUTINE GetFreqRange(r1,r2,iFreqIndex)

      INTEGER iFreqIndex
      REAL r1,r2  

c now get the freq range
      IF (iFreqIndex .EQ. 1) THEN
        r1=755.0
        r2=780.0
c        r1=655.0
c        r2=680.0
      ELSE IF (iFreqIndex .EQ. 2) THEN
        r1=1005.0
        r2=1055.0
      ELSE IF (iFreqIndex .EQ. 3) THEN
        r1=1230.0
        r2=1255.0
      ELSE IF (iFreqIndex .EQ. 4) THEN
        r1=1430.0
        r2=1505.0
      ELSE IF (iFreqIndex .EQ. 5) THEN
        r1=1530.0
        r2=1605.0
      ELSE IF (iFreqIndex .EQ. 6) THEN
        r1=2355.0
        r2=2405.0
      ELSE IF (iFreqIndex .EQ. 7) THEN
        r1=2605.0
        r2=2630.0
      ELSE IF (iFreqIndex .EQ. 8) THEN   !!!these are special one run cases
        r1=1205.0
        r2=1805.0
        END IF  
      RETURN
      END

c************************************************************************
c this subroutine gets the surface temp
      SUBROUTINE GetSurfaceTemp(InfileIndex,rSurfaceTemp)

      CHARACTER*80 Inputfile
      CHARACTER*130 caStr130
      CHARACTER*2 ca2,ca2t

      INTEGER InFileIndex,iCnt,iIOIN,i1,iFind
      REAL r1,r2,r3,r4,r5,r6,rSurfaceTemp

      CALL DoClear(InputFile)
      Inputfile='../DATA/TestProf/prof'
      write(ca2t,4000)  InfileIndex
 4000 FORMAT(I2)
      ca2='  '
      IF (InfileIndex .LT. 10) THEN
        ca2(1:1)=ca2t(2:2)
      ELSE
        ca2(1:2)=ca2t(1:2)
        END IF
      CALL FindBlank(InputFile,iFind)
      InputFile(iFind:iFind+1)=ca2(1:2)

c get the surface temp info from the profile
      iCnt=0
      iIOIN=10
      OPEN(UNIT=iIOIN,file=Inputfile,status='old',form='formatted')
 40   CONTINUE 
      IF (iCnt .GT. 0) THEN
        GO TO 50
        END IF 
      READ(iIOIN,5013,END=50) caStr130
      IF ((iCnt .EQ. 0) .AND. (caStr130(1:1) .NE. '!') .AND. 
     $    (caStr130(1:1) .NE. 'r')) THEN
        iCnt=iCnt+1
        READ(caStr130,*) i1,r1,r2,r3,rSurfaceTemp,r4,r5,r6
        END IF
      GO TO 40
 50   CONTINUE
      CLOSE(iIOIN)
ccc      print *,'the surface temp is ',rSurfaceTemp

 5013 FORMAT(A130)
 5000 FORMAT(A80)
 
      RETURN
      END

c************************************************************************ 
c this subroutine clears a string
      SUBROUTINE doclear(caStr)

      CHARACTER*80 caStr
      INTEGER iInt

      do iINt=1,80
        caStr(iInt:iInt)=' '
        end do
  
      RETURN
      END
c************************************************************************ 
c this subroutine finds the first nonblank character in a string
      SUBROUTINE FindBlank(caStr,iFind)

      CHARACTER*80 caStr
      INTEGER iFind

      iFind=1
 10   CONTINUE
      IF ((caStr(iFind:iFind) .NE. ' ') .AND. (iFind .LT. 80)) THEN
        iFind=iFind+1
        GO TO 10
        END IF

      RETURN
      END
c************************************************************************
c this subroutine finds the posn of . in blah.dat
      SUBROUTINE FindPeriod(caStr,iFind)

      CHARACTER*80 caStr
      INTEGER iFind,iTruth

      iFind=0
      iTruth=-1
 10   CONTINUE
      IF ((iTruth .LT. 0) .AND. (iFind .LT. 76)) THEN
        iFind=iFind+1
        IF (caStr(iFind:iFind+3) .EQ. '.dat') THEN
          iTruth=1
          END IF
        GO TO 10
        END IF

      RETURN
      END
c************************************************************************
c this subroutine (1)reads in the template input file and writes out
c                    the correct version
c                 (2)creates an inputfile for atmos.x

      SUBROUTINE Genln2(r1,r2,InfileIndex,iFreqIndex,rSurfT,caSubDir)

      REAL r1,r2,rSurFT,pr1,pr2,angle,tspace,rr,height
      INTEGER InfileIndex,iFreqIndex,iAtm

      INTEGER iIOIN,iIOOUT,iFind,iFind2
      CHARACTER*80 Inputfile,OutputFile
      CHARACTER*80 ca80,caSubDir
      CHARACTER*1 c1
      CHARACTER*2 c2,c2t
c caSub1 === frequency start/stop points
c caSub2 === binary output file name for atmos.x
c caSub3 === input profile file name for the genln2 file read into kcarta.x
c caSub4 === radiance stuff
      CHARACTER*80 caSub1,caSub2,caSub3,caSub4,caPath
 
c get the name of the template genln2 file
      CALL DoClear(InputFile)
      InputFile='../DATA/Template/kcartainput_template'

      !depending ion iFreIndex = 1,2,3  we get a,b,c,d ... 
      c1= char(ichar('a')+(iFreqIndex-1)) 

c get the name of the output KCARTA file
      c2='  '
      WRITE(c2t,5008) InfileIndex
 5008 FORMAT(I2)
      IF (InfileIndex .LT. 10) THEN
        c2(1:1)=c2t(2:2)
      ELSE
        c2(1:2)=c2t(1:2)
        END IF
      CALL DoClear(OutputFile)
      OutPutFile='../'
      CALL FindBlank(OutPutFile,iFind)
      CALL FindBlank(caSubDir,iFind2)
      OutPutFile(iFind:iFind+iFind2-1)=caSubDir(1:iFind2-1)
      CALL FindBlank(OutPutFile,iFind)
      OutPutFile(iFind:iFind+4)='/test'

      CALL FindBlank(OutputFile,iFind)
      OutputFile(iFind:iFind)=c1(1:1)
      OutputFile(iFind+1:iFind+2)=c2(1:2)

      CALL FindBlank(OutputFile,iFind)
      OutputFile(iFind:iFind+2)='.ip'

c get the string caSub1 (*DEFGRD frequency wavenumber substitution)
      write(caSub1,5500) r1,r2
 5500 FORMAT(2(F10.5,'   '))

c get the string caSub2 (*OUTPUT binary file from atmos.x)
      CALL doclear(caSub2)
      caSub2='''../'
      CALL FindBlank(caSub2,iFind)
      CALL FindBlank(caSubDir,iFind2)
      caSub2(iFind:iFind+iFind2-1)=caSubDir(1:iFind2-1)
      CALL FindBlank(caSub2,iFind)
      caSub2(iFind:iFind+4)='/test'

      CALL FindBlank(caSub2,iFind)
      caSub2(iFind:iFind)=c1(1:1)
      caSub2(iFind+1:iFind+2)=c2(1:2)
      IF (InfileIndex .LT. 10) THEN  
        caSub2(iFind+2:iFind+6)='.dat'''
      ELSE 
        caSub2(iFind+3:iFind+7)='.dat'''
        END IF

c get the string caSub3 (filename for *PTHFIL)
      CALL doclear(caSub3)
      caSub3='''../'
      CALL FindBlank(caSub3,iFind)
      CALL FindBlank(caSubDir,iFind2)
      caSub3(iFind:iFind+iFind2-1)=caSubDir(1:iFind2-1)
      CALL FindBlank(caSub3,iFind)
      caSub3(iFind:iFind+8)='/testprof'

      CALL FindBlank(caSub3,iFind)
      caSub3(iFind:iFind+1)=c2(1:2)
      IF (InfileIndex .LT. 10) THEN  
        caSub3(iFind+1:iFind+1)=''''
      ELSE 
        caSub3(iFind+2:iFind+2)=''''
        END IF

      iIOIN=10
      iIOOUT=11
      OPEN(UNIT=iIOIN,file=Inputfile,status='old',form='formatted')
      OPEN(UNIT=iIOOUT,file=Outputfile,status='unknown',
     $     form='formatted')

c this looks for the special flags FLAG1..42, so it knows the next
c lines have to be substituted with the caSub1,caSub2 terms
c else it just copies the line over to the new file
 20   READ(iIOIN,5000,END=30) ca80
       IF (ca80(1:6) .EQ. '!FLAG1') THEN
        WRITE(iIOOUT,5000)ca80
        READ(iIOIN,5000,END=30) ca80
        WRITE(iIOOUT,5000)caSub1
      ELSE IF (ca80(1:6) .EQ. '!FLAG2') THEN
        WRITE(iIOOUT,5000)ca80
        READ(iIOIN,5000,END=30) ca80
        WRITE(iIOOUT,5000)caSub2
      ELSE IF (ca80(1:6) .EQ. '!FLAG3') THEN
        WRITE(iIOOUT,5000)ca80
        READ(iIOIN,5000,END=30) ca80
        WRITE(iIOOUT,5000)caSub3
      ELSE IF (ca80(1:6) .EQ. '!FLAG4') THEN
        WRITE(iIOOUT,5000)ca80
        READ(iIOIN,5000,END=30) ca80
        READ(ca80,*)iAtm,pr1,pr2,tspace,rr,angle,height
        rr=rSurfT
        write(caSub4,5555)iAtm,pr1,pr2,tspace,rr,angle,height
        WRITE(iIOOUT,5000)caSub4
      ELSE
        WRITE(iIOOUT,5000)ca80
        END IF
      GOTO 20
 5000 FORMAT(A80)
 5555 FORMAT(I3,' ',6(F10.5,' '))

30    CONTINUE 
      CLOSE(iIOIN)
      CLOSE(iIOOUT)

c no longer needed
ccc now create the input file for kcarta.x 
ccc      Inputfile='../SRC/inputfile'
ccc      CALL DoClear(InputFile)
ccc      InputFile='inputfile'
ccc      OPEN(UNIT=iIOIN,file=Inputfile,status='unknown',form='formatted')
ccc      CALL DoClear(caSub1)
ccc      CALL FindPeriod(caSub2,iFind)
cccc thus if caSub2 = '../testing.dat' then caSub1 --> ../testing (w/o quotes!)
ccc      caSub1(1:iFind-2)=caSub2(2:iFind-1)
cccc now caSub1 --> ../testing + .ip   (w/o quotes!)
ccc      CALL FindBlank(caSub1,iFind)
ccc      caSub1(iFind:iFInd+2)='.ip'
ccc      WRITE(iIOIN,5000) caSub1 
ccc      CLOSE(iIOIN)

      RETURN
      END
c************************************************************************
c this subroutine takes in the user supplied name, and strips off all the
c leftmost blanks, and outputs the string so that it is right padded with
c blanks. this makes appending substrings to the end a bit harder
c unfortunately, this is what the SGI's need
      SUBROUTINE rightpad(caFRootName)

      CHARACTER*80 caFRootName,caTempName
      INTEGER iR,iL,iInt

c find the "right" length of the input root name
      iR=len(caFRootName)
 11   continue
      IF (caFRootName(iR:iR) .eq. ' ') THEN
        iR=iR-1
        GO TO 11      
        END IF

c find the "left" length of the input root name
      iL=1
 12   continue
      IF (caFRootName(iL:iL) .eq. ' ') THEN
        iL=iL+1
        GO TO 12      
        END IF
c thus the entered word exists between iL:iR

c now rearrange the string, so that it is right padded with blanks
c this is basically equivalent to  ADJUSTR(caFRootName)
      DO iInt=1,80
        caTempName(iInt:iInt)=' '
        END DO
      caTempName(1:iR-iL+1)=caFRootName(iL:iR)
      caFRootName(1:80)=caTempName(1:80)

      RETURN
      END

c************************************************************************
