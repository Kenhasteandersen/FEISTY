!
! newly added useful functions : linspace, diag
!

module globals
   use inputnml
   implicit none
   integer, parameter :: dp = selected_real_kind(13) ! dp =  kind(0.d0) ! double precision

   real(dp), parameter :: pi = 4*ATAN(1._dp)

   ! Small number to avoid divisions by zero
   real(dp), parameter :: eps = 1.0e-200_dp

   ! Temperature Q10 corrections
   real(dp), parameter :: Q10=1.88_dp, Q10m=1.88_dp, Q10mPetrik=2.35_dp
   real(dp) :: fTemp, fTempm, fTempdem, fTempmdem, fTempdem_shallow, fTempmdem_shallow
   real(dp),save :: fTempold=1._dp, fTempmold=1._dp, fTempdemold=1._dp, fTempmdemold=1._dp
   real(dp),save :: fTempdem_shallowold=1._dp, fTempmdem_shallowold=1._dp
   real(dp), parameter:: Tref = 10._dp ! Reference temperature

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
    res=0._dp
    do i = 1, size(A)
        res(i,i) = A(i)
    end do
    end function

end module globals
