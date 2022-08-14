! debug tool only works without input things?
program FEISTYtest
   use FEISTY
   use globals
   implicit none

   real(dp), allocatable:: u0(:), dudt(:)
   !real(dp), allocatable:: flvl_r(:), mortpred_r(:), g_r(:)

   call setupFEISTY(100.d0)

   allocate (u0(nGrid))
   allocate (dudt(nGrid))
   !allocate (flvl_r(nGrid))
   !allocate (mortpred_r(nGrid))
   !allocate (g_r(nGrid-nResources))

   u0(1) = 100.d0
   u0(2) = 100.d0
   u0(3) = 5.d0
   u0(4) = 0.d0
   u0(idxF:nGrid) = 1.d0
   dudt = 0.d0

   !call getrates(u0, dudt,flvl_r,mortpred_r,g_r)
   call calcderivatives(u0, dudt)

   !print *, dudt
   print *, 'hello world'

end program FEISTYtest
