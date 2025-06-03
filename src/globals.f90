!
! newly added useful functions : linspace, diag
!

module globals
   
   use fabm_types
   
   implicit none
   !integer, parameter :: rk = selected_real_kind(13) ! dp =  kind(0.d0) ! double precision

   real(rk), parameter :: pi = 4*ATAN(1._rk)

   ! Small number to avoid divisions by zero
   real(rk), parameter :: eps = 1.0e-200_rk

   ! Temperature Q10 corrections
   real(rk), parameter :: Q10=1.88_rk, Q10m=1.88_rk, Q10mPetrik=2.35_rk
   real(rk) :: fTemp, fTempm, fTempdem, fTempmdem, fTempdem_shallow, fTempmdem_shallow
   real(rk),save :: fTempold=1._rk, fTempmold=1._rk, fTempdemold=1._rk, fTempmdemold=1._rk
   real(rk),save :: fTempdem_shallowold=1._rk, fTempmdem_shallowold=1._rk
   real(rk), parameter:: Tref = 10._rk ! Reference temperature

   real(rk), allocatable :: fTempV(:), fTempmV(:)

contains

  ! linspace function:
      function linspace(value_start, value_end, length) result(res)
         real(rk), intent(in)::value_start, value_end
         integer, intent(in)::length
         real(rk) ::dx, res(length)
         integer::i

         dx = (value_end - value_start)/(length - 1)
         res(1:length) = [(value_start + (i - 1)*dx, i=1, length)]
      end function linspace

      ! a simple function to create diagonal matrix, only for 1-D array
    function diag(A) result(res)
    real(rk), intent(in) :: A(:)
    real(rk) :: res(size(A),size(A))
    integer :: i
    res=0._rk
    do i = 1, size(A)
        res(i,i) = A(i)
    end do
    end function

end module globals
