! for debug
program FEISTYtest
   use FEISTY
   use globals
   implicit none

   real(dp), allocatable:: u0(:), dudt(:)
   !real(dp),intent(out):: thetaF(nGrid,nGrid)
   !real(dp), allocatable:: flvl_r(:), mortpred_r(:), g_r(:)

   call setupbasic(100.d0,5.d0,20.d0,9.d0)
  !call setupbasic2(100.d0,5.d0,9,10.d0)
   !call setupVertical(80.d0,3,4) ! ?
   !call setupsquid( 50.d0, 100.d0 , 6)
   !call setupVerticalGlobal(0.1d0,0.1d0,10.d0,2580.d0,219.4018d0,&
    !                     (/0.d0, 50.d0,130.d0,240.d0,380.d0, 50.d0,750.d0,980.d0&
   !             ,1240.d0, 1530.d0, 1850.d0, 2200.d0, 2580.d0, 2990.d0, 3430.d0, 3900.d0, 4400.d0, 4930.d0, 5490.d0, 6080.d0/)&
   ! ,(/-1.8593931d0, -1.8593931d0, -1.8547935d0, -1.7973094d0, -1.3357751d0,  0.5777518d0,  0.8374313d0&
!,1.0768622d0, 1.1814129d0,  1.0766119d0,  0.7584100d0,  0.4179280d0,  0.3816292d0,  0.d0, 0.d0,  0.d0, 0.d0, 0.d0, 0.d0, 0.d0/),6)

   allocate (u0(nGrid))
   allocate (dudt(nGrid))
   !allocate (flvl_r(nGrid))
   !allocate (mortpred_r(nGrid))
   !allocate (g_r(nGrid-nResources))

!   u0(1) = 5.d0
!   u0(2) = 5.d0
!   u0(3) = 0.d0
!   u0(4) = 0.d0
!   u0(idxF:nGrid) = 0.0001d0
!   dudt = 0.d0

   u0(1) = 100.d0
   u0(2) = 100.d0
   u0(3) = 5.d0
   u0(4) = 0.d0
   u0(idxF:nGrid) = 1.d0
   dudt = 0.d0

!! van Denderen et al., 2020
!   u0(1) = 0.5d0
!   u0(2) = 0.5d0
!   u0(3) = 0.5d0
!   u0(4) = 0.d0
!   u0(idxF:nGrid) = 0.0001d0
!   dudt = 0.d0


!call simulateEuler(u0,dudt, 100.d0, 0.01d0)   !3stages dt=0.1   6stages 0.001 9stages 0.00001




 call calcderivatives(u0,dudt)
 !call calcderivativesSquid(u0,dudt)

   !call getrates(u0, dudt,flvl_r,mortpred_r,g_r)

end program FEISTYtest
