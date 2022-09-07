!
! from NUM library
! subroutines of temperature may be used for further work
! new added useful functions : linspace, diag
!

module globals
   use input
   implicit none
   integer, parameter :: dp = kind(0.d0) ! double precision
   !
   ! Useful mathematical constants:
   !
   real(dp), parameter :: onethird = 1.d0/3.d0
   real(dp), parameter :: twothirds = 2.d0/3.d0
   real(dp), parameter :: threequarters = 3.d0/4.d0
   real(dp), parameter :: pi = 4*ATAN(1.d0)

   ! Small number to avoid divisions by zero
   real(dp), parameter :: eps = 1d-200

   ! Temperature Q10 corrections (for Q10=2 and Q10=1.5)
   real(dp) :: fTemp2, fTemp15
   real(dp), parameter:: Tref = 10. ! Reference temperature

   !
   ! Specification of what to do with HTL losses:
   !
   real(dp) :: fracHTL_to_N ! Half becomes urine that is routed back to N
   real(dp) :: fracHTL_to_POM ! Another half is fecal pellets that are routed back to the largest POM size class
   real(dp) :: rhoCN
   !real(dp), parameter :: fracHTL_to_N = 0.5 ! Half becomes urine that is routed back to N
   !real(dp), parameter :: fracHTL_to_POM = 0.5 ! Another half is fecal pellets that are routed back to the largest POM size class
   !real(dp), parameter :: rhoCN = 5.68

contains

   ! -----------------------------------------------
   ! Read in general parameters
   ! -----------------------------------------------
   subroutine read_namelist_general()
      integer :: file_unit, io_err

      namelist /input_general/ rhoCN, fracHTL_to_N, fracHTL_to_POM

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_general, iostat=io_err)
      call close_inputfile(file_unit, io_err)

   end subroutine read_namelist_general

   ! -----------------------------------------------
   ! Temperature Q10 function
   ! -----------------------------------------------
   function fTemp(Q10, T) result(f)
      real(dp), intent(in):: Q10, T
      real(dp):: f

      f = Q10**((T - Tref)/10.)
   end function fTemp

   ! -----------------------------------------------
   ! Update the temperature corrections only if T has changed
   ! -----------------------------------------------
   subroutine updateTemperature(T)
      real(dp), intent(in) :: T
      real(dp), save :: Told = -1000.

      if (T .ne. Told) then
         Told = T
         fTemp2 = fTemp(2.d0, T)
         fTemp15 = fTemp(1.5d0, T)
      end if
   end subroutine updateTemperature


  ! linspace function:
      function linspace(value_start, value_end, length) result(res)
         real(dp), intent(in)::value_start, value_end
         integer, intent(in)::length
         real(dp) ::dx, res(length)
         integer::i

         dx = (value_end - value_start)/(length - 1)
         res(1:length) = [(value_start + (i - 1)*dx, i=1, length)]
      end function linspace

      ! a simple function to create diagonal matrix, only for 1-D array
    function diag(A) result(res)
    real(dp), intent(in) :: A(:)
    real(dp) :: res(size(A),size(A))
    integer :: i
    res=0.d0
    do i = 1, size(A)
        res(i,i) = A(i)
    end do
    end function

end module globals
