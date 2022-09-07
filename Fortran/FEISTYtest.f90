! for debug
program FEISTYtest
   use FEISTY
   use globals
   implicit none

   real(dp), allocatable:: u0(:), dudt(:)
   !real(dp),intent(out):: thetaF(nGrid,nGrid)
   !real(dp), allocatable:: flvl_r(:), mortpred_r(:), g_r(:)

   !call setupbasic(100.d0,5.d0)
  !call setupbasic2(100.d0,5.d0,9)
   !call setupVertical(80.d0) ! ?
   call setupsquid( 50.d0, 100.d0 , 6)

   allocate (u0(nGrid))
   allocate (dudt(nGrid))
   !allocate (flvl_r(nGrid))
   !allocate (mortpred_r(nGrid))
   !allocate (g_r(nGrid-nResources))

   u0(1) = 5.d0
   u0(2) = 5.d0
   u0(3) = 0.d0
   u0(4) = 0.d0
   u0(idxF:nGrid) = 0.0001d0
   dudt = 0.d0

!   u0(1) = 100.d0
!   u0(2) = 100.d0
!   u0(3) = 5.d0
!   u0(4) = 0.d0
!   u0(idxF:nGrid) = 1.d0
!   dudt = 0.d0

!! van Denderen et al., 2020
!   u0(1) = 0.5d0
!   u0(2) = 0.5d0
!   u0(3) = 0.5d0
!   u0(4) = 0.d0
!   u0(idxF:nGrid) = 0.0001d0
!   dudt = 0.d0

 !call calcderivatives(u0,dudt)
 call calcderivativesSquid(u0,dudt)

   !call getrates(u0, dudt,flvl_r,mortpred_r,g_r)

end program FEISTYtest
