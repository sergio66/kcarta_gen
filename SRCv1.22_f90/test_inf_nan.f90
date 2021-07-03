      program test

use ieee_arithmetic
      real :: x(6)
      logical :: l(6)
      integer :: i

      x(1) = nan
      x(2) = inf
      x(1) = NaN
      x(2) = Infinity

      x(3) = 33.3
      x(4) = -1.0/0.0
      x(5) = 0.0/0.0

     print *,' nan inf 33.3  1/0 0/0 not-initialized'
     write(*,'(F12.4)') x

      l = ieee_is_nan(x);    print *,'nan test',l
      l = ieee_is_finite(x); print *,'finite test',l
      l = ieee_is_normal(x); print *,'normal test',l

!************************************************************************
      print *,' '
      print *,'TEST 1 '
      do i = 1,6
        if (ieee_is_nan(x(i))) then
          print *,'A1 nan',i,x(i)
        end if
!        if (.false. /= ieee_is_nan(x(i))) then
!          print *,'B finite, not nan',i,x(i)
!        end if

        if (ieee_is_finite(x(i))) then
          print *,'C1 finite',i,x(i)
        end if
!        if (.false. /= ieee_is_finite(x(i))) then
!          print *,'D finite, not finite',i,x(i)
!        end if

        if (ieee_is_normal(x(i))) then
          print *,'E1 normal',i,x(i)
        end if
!        if (.false. /= ieee_is_normal(x(i))) then
!          print *,'F normal, not normal',i,x(i)
!        end if

        print *,' '
      end do
      
!************************************************************************
      print *,' '
      print *,'TEST 2 '
      do i = 1,6
        if (.false. .eq. ieee_is_nan(x(i))) then
          print *,'B2 not nan but ok',i,x(i)
        end if

        if (.false. .eq. ieee_is_finite(x(i))) then
          print *,'D2 not finite but ok',i,x(i)
        end if

        if (.false. .eq. ieee_is_normal(x(i))) then
          print *,'F2 not normal but ok',i,x(i)
        end if

        print *,' '
      end do
      
!************************************************************************
      print *,' '
      print *,'TEST 3 '
      do i = 1,6
        if (.true. .EQ.  ieee_is_nan(x(i))) then
          print *,'B3 nan',i,x(i)
        end if

        if (.true. .EQ.  ieee_is_finite(x(i))) then
          print *,'D3 finite',i,x(i)
        end if

        if (.true. .EQ.  ieee_is_normal(x(i))) then
          print *,'F3 normal',i,x(i)
        end if

        print *,' '
      end do
      
!************************************************************************
      print *,' '
      print *,'TEST 4 '
      print *,' nan inf 33.3  1/0 0/0 not-initialized'
      write(*,'(F12.4)') x
      do i = 1,6
        if (.true. .EQ. ieee_is_nan(x(i))) then
          print *,'B4 nan',i,x(i)
        end if

        if (.not. ieee_is_finite(x(i))) then
          print *,'D4 infinite',i,x(i)
        end if

        if (.not. ieee_is_normal(x(i))) then
          print *,'F4 abnormal',i,x(i)
        end if

        print *,' '
      end do
      
!************************************************************************

      print *,' '
      print *,'TEST 5 '
      print *,' nan inf 33.3  1/0 0/0 not-initialized'
      write(*,'(F12.4)') x
      do i = 1,6
        if (.true. .EQ. is_badnum(x(i))) then
          print *,'ABCDEF5 is badnum',i,x(i)
        end if
        print *,' '
      end do

      end      
!************************************************************************

     logical function is_badnum(x)

     real :: x
     logical :: test

     test = .false.

     if (x > huge(x)) then
       print *, 'x is huge ',x
       test = .true.
     end if
     if (abs(x) > huge(x)) then
       print *, '+/-x is huge ',x
       test = .true.
     end if
       
     is_badnum = test

    RETURN
    end FUNCTION is_badnum

!************************************************************************
