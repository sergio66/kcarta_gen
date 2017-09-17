      !ifort -extend-source 132 f77_nchar_to_120char.f90 -o f77_nchar_to_120char  
                                                                                
      PROGRAM f77_nchar_to_120char                                              
                                                                                
      IMPLICIT none                                                             
                                                                                
      INTEGER iIOUN1,iIOUN2,IERR,iCnt,iLong,iLen,xlen,iLenCut                   
      CHARACTER*80 caF1                                                         
      CHARACTER*81 caF2                                                         
      CHARACTER*132 caIn,caX                                                    
                                                                                
      INTEGER iI                                                                
                                                                                
      iLenCut = 100
      iLenCut = 72      
      iLenCut = 80
                                                                                
      print *,' input caF1 : '                                                  
      read *,caF1                                                               
      caF2 = 'x' // caF1                                                        
                                                                                
      iIOUN1 = 11                                                               
      iIOUN2 = 12                                                               
      OPEN(UNIT=iIOUN1,FILE=caF1,STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR)                                                          
      OPEN(UNIT=iIOUN2,FILE=caF2,STATUS='UNKNOWN',FORM='FORMATTED',IOSTAT=IERR)                                                          
                                                                                
      iCnt = 0                                                                  
      iLong = 0                                                                 
 101  CONTINUE                                                                  
      read(iIOUN1,'(A)',END=199) caIn                                           
      iCnt = iCnt + 1
      CALL dollarremove(caIn)
      CALL bangremove(caIn)
      caIn = trim(caIn)                                                         
      ilen = xlen(caIn)
!      print *,iCnt,iLen
      IF (iLen .GT. iLenCut) THEN                                               
!        write(*,'(I6,I6,A)') iCnt,iLen,caIn                                    
      iLong = iLong + 1                                                         
      CALL CutStr(iLenCut,iCnt,caIn,caX)
        ilen = xlen(caIn)      
        write(iIOUN2,'(A)') caIn(1:iLen)
        ilen = xlen(caX)      	
        write(iIOUN2,'(A)') caX(1:iLen)                                      
      ELSE                                                                      
        !!! no problemo                                                         
        write(iIOUN2,'(A)') caIn(1:iLen)                                     
      END IF                                                                    
                                                                                
      GOTO 101                                                                  
                                                                                
 199  CONTINUE                                                                  
      CLOSE(iIOUN1)                                                             
      CLOSE(iIOUN2)                                                             
                                                                                
      write(*,'(A,A,A,I4,A,I4,A,I3,A)') 'file ',trim(caF1),' had ',iCnt,' lines where',iLong,' were longer than ',iLenCut,' char'
      END                                                                       
                                                                                
!************************************************************************       
                                                                                
      SUBROUTINE dollarremove(caIn)

      IMPLICIT NONE

      CHARACTER*132 caIn
      
      INTEGER iLen,xLen,iI
      
      iLen = xlen(caIn)
      DO iI = 1,iLen
	IF (caIn(iI:iI) .EQ. '$') caIn(iI:iI) = 'c'
      END DO
      
      RETURN
      END

!************************************************************************       
                                                                                
      SUBROUTINE bangremove(caIn)

      IMPLICIT NONE

      CHARACTER*132 caIn,caX,ca0

      INTEGER iLen,xLen,iStop,iI,iRepeatBangFound,iCnt,iJ

      ca0 = caIn
      iStop = -1
      iCnt = 0
      DO WHILE (iStop .LT. 0)
        iCnt = iCnt + 1
        iLen = xlen(caIn)
	DO iI = 1,132
	  caX(iI:iI) = ' '
	END DO
	caX = caIn
	iRepeatBangFound = -1
	DO iI = 1,iLen
	  IF ((caX(iI:iI) .EQ. '!') .AND. (caX(iI+1:iI+1) .EQ. '!')) THEN
	    iRepeatBangFound = +1
	    iCnt = iCnt + 1
	    DO iJ = iI+1,132
	      caIn(iJ:iJ) = ' '
	    END DO
	    caIn(iI+1:131) = caX(iI+2:132)
	  END IF	  
	END DO
	IF (iRepeatBangFound .LT. 0) iStop = 1
      END DO
      
      RETURN
      END

!************************************************************************       
                                                                                
      INTEGER FUNCTION xlen(caStrx)                                             
                                                                                
      IMPLICIT NONE                                                             
                                                                                
      CHARACTER*132 caStrx                                                      
      INTEGER iI,iFound                                                         
                                                                                
      iFound = -1                                                               
      iI = len(caStrx)                                                          
      DO WHILE ((iI .GT. 1) .AND. (iFound .LT. 0))                              
        IF (caStrx(iI:iI) .NE. ' ') THEN                                        
          iFound = iI                                                           
        ELSE                                                                    
          iI = iI - 1                                                           
        END IF                                                                  
      END DO                                                                    
                                                                                
      xlen = iI                                                                 
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
!************************************************************************       
                                                                                
      SUBROUTINE CutStr(iLenCut,iCnt,ca1,ca2)                                   
                                                                                
      IMPLICIT NONE                                                             
                                                                                
      INTEGER iLenCut,iCnt                                                      
      CHARACTER*132 ca1,ca2                                                     
                                                                                
      CHARACTER*132 ca0                                                         
                                                                                
      INTEGER iI,iLen,xlen,iFound                                               
                                                                                
      ca0 = ca1                                                                 
      iLen = xlen(ca0)                                                          
                                                                                
      DO iI = 1,132                                                             
        ca2(iI:iI) = ' '                                                        
      END DO                                                                    
                                                                                
      iFound = -1                                                               
      iI = iLenCut-1                                                            
      DO WHILE ((iI .GT. 1) .AND. (iFound .LT. 0))                              
        IF ((ca1(iI:iI) .EQ. ' ') .OR. (ca1(iI:iI) .EQ. ',')) THEN              
          iFound = iI                                                           
        ELSE                                                                    
          iI = iI - 1                                                           
        END IF                                                                  
      END DO                                                                    
                                                                                
      ca2(1:7) = '     c '                                                      
      ca2(8:8+iLen-iFound-1) = ca1(iFound+1:iLen)                                   
                                                                                
      DO iI = iFound+1,132                                                      
        ca1(iI:iI) = ' '                                                        
      END DO                                                                    
                                                                                
      print *,' '                                                               
      write(*,'(A,I4,A,A,A)') 'Line ',iCnt,' Split the input string <',ca0(1:iLen),'> into the following'
      write(*,'(A,A,A)') ca1,' and ', ca2                                       
      print *,' '                                                               
                                                                                
      RETURN                                                                    
      END                                                                       
                                                                                
!************************************************************************       
