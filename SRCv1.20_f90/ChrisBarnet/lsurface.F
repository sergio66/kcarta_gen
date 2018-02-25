c***********************************************************************
      function lsurface ( numlev, pres, psurf, plow, phigh )
      implicit none
      integer*4  lsurface
c***********************************************************************
c
c        Routine to find the pointer to the bottom of atmosphere.
c
c***********************************************************************
      integer*4 numlev
      real*4    pres(*), Psurf, Plow, Phigh

C     local variables
c     ---------------

      integer *4 L

c***********************************************************************
c***********************************************************************
c
c     The bottom layer is at least 5 mb thick
c
c     Surface pressure between plow, phigh
c     ------------------------------------

      if ( psurf .gt. plow .and. psurf .le. phigh ) then

         do L = numlev, 1, -1
            if ( psurf .ge. pres(L-1)+5.0 ) then
               lsurface = L
               goto 990
            end if
         end do
      end if

c     surface pressure exceeded limits, see if it is within Pobs space
c     ----------------------------------------------------------------

      if(psurf.le.plow) then
        lsurface = 1
        do L = 1, numlev
          if(psurf.ge.pres(L) + 5.0) lsurface = L
        enddo
      else
        lsurface = numlev
        do L = numlev,1,-1
          if(psurf.le.pres(L) - 5.0) lsurface = L
        enddo
      endif

      print 100, plow, psurf, phigh, lsurface, pres(lsurface)

  990 return
  100 format('lsurface: ',f7.2,' <= (psurf=',f7.2,') <= ',f7.2,
     1    ' pres(',i3,')=',f7.2)
      end

