!
! newly added useful functions : linspace, diag
!

module globals
   use input
   implicit none
   integer, parameter :: dp = kind(0.d0) ! double precision

   real(dp), parameter :: pi = 4*ATAN(1.d0)

   ! Small number to avoid divisions by zero
   real(dp), parameter :: eps = 1d-200

   ! Temperature Q10 corrections
   real(dp), parameter :: Q10=1.88d0, Q10m=1.88d0, Q10mPetrik=2.35d0
   real(dp) :: fTemp, fTempm, fTempdem, fTempmdem, fTempdem_shallow, fTempmdem_shallow
   real(dp),save :: fTempold=1.d0, fTempmold=1.d0, fTempdemold=1.d0, fTempmdemold=1.d0
   real(dp),save :: fTempdem_shallowold=1.d0, fTempmdem_shallowold=1.d0
   real(dp), parameter:: Tref = 10.d0 ! Reference temperature

   real(dp), allocatable :: fTempV(:), fTempmV(:)

contains

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
