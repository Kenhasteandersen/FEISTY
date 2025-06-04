! for debug
program FEISTYtest
   use setup
   use globals
   implicit none

   real(dp), allocatable:: u0(:), dudt(:)
   !real(dp),intent(out):: thetaF(nGrid,nGrid)
   !real(dp), allocatable:: flvl_r(:), mortpred_r(:), g_r(:)

   !call setupbasic(100._dp,100._dp,-1._dp,5._dp,100._dp,10._dp,8._dp)
  !call setupbasic2(100._dp,100._dp,-5._dp,5._dp,6,1500._dp,10._dp,10._dp,0.002_dp,0._dp,0.05_dp)
   !call setupVertical(80._dp,80._dp,-150._dp,-100._dp,100._dp,4,1500._dp,150._dp) !
     call setupVertical2(80._dp,80._dp,-150._dp,-100._dp,100._dp,9,10._dp,10._dp,10._dp,1500._dp,150._dp,0.002_dp,&
                       & 250._dp,1.5_dp,0._dp,0.05_dp)
   !call setupsquid( 50._dp, 100._dp , 6)
!   call setupVerticalGlobal(0.1_dp,0.1_dp,10._dp,2580._dp,219.4018_dp,&
!                         (/0._dp, 50._dp,130._dp,240._dp,380._dp, 50._dp,750._dp,980._dp&
!                ,1240._dp, 1530._dp, 1850._dp, 2200._dp, 2580._dp, 2990._dp, 3430._dp, 3900._dp, 4400._dp, 4930._dp, 5490._dp, 6080._dp/)&
!    ,(/-1.8593931_dp, -1.8593931_dp, -1.8547935_dp, -1.7973094_dp, -1.3357751_dp,  0.5777518_dp,  0.8374313_dp&
!,1.0768622_dp, 1.1814129_dp,  1.0766119_dp,  0.7584100_dp,  0.4179280_dp,  0.3816292_dp,  0._dp, 0._dp,  0._dp, 0._dp, 0._dp, 0._dp, 0._dp/)&
!,6,0.25_dp)

   allocate (u0(nGrid))
   allocate (dudt(nGrid))
   !allocate (flvl_r(nGrid))
   !allocate (mortpred_r(nGrid))
   !allocate (g_r(nGrid-nResources))

!   u0(1) = 5._dp
!   u0(2) = 5._dp
!   u0(3) = 0._dp
!   u0(4) = 0._dp
!   u0(idxF:nGrid) = 0.0001_dp
!   dudt = 0._dp

!   u0(1) = 100._dp
!   u0(2) = 100._dp
!   u0(3) = 5._dp
!   u0(4) = 0._dp
!   u0(idxF:nGrid) = 1._dp
!   dudt = 0._dp

! van Denderen et al., 2020
   u0(1) = 0.5_dp
   u0(2) = 0.5_dp
   u0(3) = 0.5_dp
   u0(4) = 0._dp
   u0(idxF:nGrid) = 0.0001_dp
   dudt = 0._dp


!call simulateEuler(u0,dudt, 100._dp, 0.01_dp)   !3stages dt=0.1   6stages 0.001 9stages 0.00001




 call calcderivatives(u0,dudt)
    call setupbasic(100._dp,100._dp,5._dp,-1._dp,100._dp,10._dp,9._dp)
     call calcderivatives(u0,dudt)
 !call calcderivativesSquid(u0,dudt)

   !call getrates(u0, dudt,flvl_r,mortpred_r,g_r)

end program FEISTYtest
