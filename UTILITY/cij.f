c this just figures out the cij lists
c compile with f77 -o cij.x -N109 -W cij.f
      implicit none

      character*80 caIn,caOut
      character*80 caLine
      character*5 caArray
      integer iIOUN1,iIOUN2,iErr,iI,iJ,iK,iChannel
      real    rF,rNoise,rFWHM,rCij

      caIn  = '/asl/data/airs/srf/clist.jpl/L2.chan_prop.2003.01.10.v6.6.1.anc'      
      caIn  = '/asl/data/airs/srf/clist.jpl/L2.I.channel_prop.v5.1.2.txt'
      caOut = '/home/sergio/MATLAB/CLOUD/cij.txt'

      iIOUN1 = 10
      iIOUN2 = 11
      OPEN(UNIT=iIOUN1,FILE=caIn,STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR) 
      OPEN(UNIT=iIOUN2,FILE=caOut,STATUS='UNKNOWN',FORM='FORMATTED',IOSTAT=IERR) 

      iI = 0
      iJ = 0
 10   CONTINUE
      READ(iIOUN1,50,END=1030) caLine
      print *,caLine
      IF (caLine(1:1) .EQ. '!') THEN
        iJ = iJ + 1
        print *,iJ
        GOTO 10
      ELSE
        iI = iI + 1
        READ (caLine,*) iK,rF,caArray,iChannel,rNoise,rFWHM,rCij
        write (iIOUN2,*) iK,rF,iChannel,rNoise,rFWHM,rCij
        GOTO 10
        END IF

 1030 CONTINUE
      CLOSE(iIOUN1) 
      CLOSE(iIOUN2) 

 50   format(a80)

      print *,' read in ',iI,' lines'
      end

