! for debug
program FEISTYtest
   use setup
   use globals
   implicit none

   real(rk), allocatable:: u0(:), dudt(:)
   !real(rk),intent(out):: thetaF(nGrid,nGrid)
   !real(rk), allocatable:: flvl_r(:), mortpred_r(:), g_r(:)

   !call setupbasic(100._rk,100._rk,-1._rk,5._rk,100._rk,10._rk,8._rk)
  !call setupbasic2(100._rk,100._rk,-5._rk,5._rk,6,1500._rk,10._rk,10._rk,0.002_rk,0._rk,0.05_rk)
   !call setupVertical(80._rk,80._rk,-150._rk,-100._rk,100._rk,4,1500._rk,150._rk) !
     call setupVertical2(80._rk,80._rk,-150._rk,-100._rk,100._rk,9,10._rk,10._rk,10._rk,1500._rk,150._rk,0.002_rk,&
                       & 250._rk,1.5_rk,0._rk,0.05_rk)
   !call setupsquid( 50._rk, 100._rk , 6)
!   call setupVerticalGlobal(0.1_rk,0.1_rk,10._rk,2580._rk,219.4018_rk,&
!                         (/0._rk, 50._rk,130._rk,240._rk,380._rk, 50._rk,750._rk,980._rk&
!                ,1240._rk, 1530._rk, 1850._rk, 2200._rk, 2580._rk, 2990._rk, 3430._rk, 3900._rk, 4400._rk, 4930._rk, 5490._rk, 6080._rk/)&
!    ,(/-1.8593931_rk, -1.8593931_rk, -1.8547935_rk, -1.7973094_rk, -1.3357751_rk,  0.5777518_rk,  0.8374313_rk&
!,1.0768622_rk, 1.1814129_rk,  1.0766119_rk,  0.7584100_rk,  0.4179280_rk,  0.3816292_rk,  0._rk, 0._rk,  0._rk, 0._rk, 0._rk, 0._rk, 0._rk/)&
!,6,0.25_rk)

   allocate (u0(nGrid))
   allocate (dudt(nGrid))
   !allocate (flvl_r(nGrid))
   !allocate (mortpred_r(nGrid))
   !allocate (g_r(nGrid-nResources))

!   u0(1) = 5._rk
!   u0(2) = 5._rk
!   u0(3) = 0._rk
!   u0(4) = 0._rk
!   u0(idxF:nGrid) = 0.0001_rk
!   dudt = 0._rk

!   u0(1) = 100._rk
!   u0(2) = 100._rk
!   u0(3) = 5._rk
!   u0(4) = 0._rk
!   u0(idxF:nGrid) = 1._rk
!   dudt = 0._rk

! van Denderen et al., 2020
   u0(1) = 0.5_rk
   u0(2) = 0.5_rk
   u0(3) = 0.5_rk
   u0(4) = 0._rk
   u0(idxF:nGrid) = 0.0001_rk
   dudt = 0._rk


!call simulateEuler(u0,dudt, 100._rk, 0.01_rk)   !3stages dt=0.1   6stages 0.001 9stages 0.00001




 call calcderivatives(u0,dudt)
    call setupbasic(100._rk,100._rk,5._rk,-1._rk,100._rk,10._rk,9._rk)
     call calcderivatives(u0,dudt)
 !call calcderivativesSquid(u0,dudt)

   !call getrates(u0, dudt,flvl_r,mortpred_r,g_r)

end program FEISTYtest
