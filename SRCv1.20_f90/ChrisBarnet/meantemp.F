

      subroutine meantemp(numlev, Pobs, Psurf, Tair, Tavg)
      implicit none

c     input variables
c     ---------------
      integer*4 numlev  ! # of levels for this profile
      real*4  Psurf     ! surface pressure
      real*4  Pobs(*)   ! pressure of levels
      real*4  Tair(*)   ! temperature at Pobs(L)

c     output variables
c     ----------------

      real*4  Tavg(*)   ! average layer temperature

c     local variables
c     ---------------

      integer*4 L
      real*4    half
      parameter (half=0.5)
      real*4    Tsurfair   ! temperature at surface


      Tsurfair = Tair(numlev-1) + (Tair(numlev)-Tair(numlev-1))
     $    * (Psurf-Pobs(numlev-1))/(Pobs(numlev)-Pobs(numlev-1))
      Tavg(1) = Tair(1)
      do L = 2, numlev-1
        Tavg(L) = half*(Tair(L-1)+Tair(L))
      enddo
      Tavg(numlev) = half*(Tair(numlev-1)+Tsurfair)

      return
      end
