! Handle fish groups

module fish
   use globals
   use spectrum
   use input

   implicit none

   private
!physiology:
   real(dp) :: h
   real(dp) :: nn
   real(dp) :: q
   real(dp) :: gamma
   real(dp) :: kk
   real(dp) :: p
   real(dp) :: epsAssim
   real(dp) :: epsRepro
!   real(dp), parameter ::  h = 20.d0            ! Max. consumption coefficient
!   real(dp), parameter ::  nn = -0.25d0         ! Max. consumption exponent
!   real(dp), parameter ::  q = -0.2d0           ! Clearance rate exponent
!   real(dp), parameter ::  gamma = 70.d0        ! Coef. for clearance rate
!   real(dp), parameter ::  kk = 0.011d0*365.d0  ! Metabolism coefficient   4.0150 in input file
!   real(dp), parameter ::  p = -0.175d0         ! Metabolism exponent
!   real(dp), parameter ::  epsAssim = 0.7d0     ! Assimilation efficiency
!   real(dp), parameter ::  epsRepro = 0.01d0    ! reproduction * recruitment efficiency

   real(dp) :: beta
   real(dp) :: sigma
   real(dp) :: mMin
!   real(dp), parameter ::  beta = 400.d0
!   real(dp), parameter ::  sigma = 1.3d0
!   real(dp), parameter ::  mMin = 0.001d0       !min fish mass (boundary of the grid) used for discretization

   real(dp) :: mMedium
   real(dp) :: mLarge
!   real(dp), parameter ::  mMedium = 10.d0      !meidum fish central mass used for feeding preference calc
!   real(dp), parameter ::  mLarge = 5000.d0     !large fish central mass used for feeding preference calc

!resources:
   real(dp) :: lbenk
   real(dp) :: szoog
   real(dp) :: lzoog
   real(dp) :: sbeng
   real(dp) :: lbeng
!   real(dp), parameter ::     lbenk = 0.d0              ! large benthos carry capacity
!   real(dp), parameter ::     szoog = 1.d0              ! small zooplankton growth rate
!   real(dp), parameter ::     lzoog = 1.d0               ! large zooplankton growth rate
!   real(dp), parameter ::     sbeng = 1.d0              ! small benthos growth rate
!   real(dp), parameter ::     lbeng = 0.d0               ! large benthos growth rate

   type, extends(typeSpectrum) :: spectrumFish

   contains
      procedure, pass :: initFish
      procedure :: calcfluxfish
   end type spectrumFish

   public h, nn, q, gamma, kk, p, epsAssim, epsRepro, beta, sigma, mMin, mMedium, mLarge,  &
          spectrumFish, &
          initFish, calcfluxfish, &
          lbenk, szoog, lzoog, sbeng, lbeng

contains

   subroutine initFish(this, n, mMax, mMature)
      class(spectrumFish), intent(inout) :: this
      integer, intent(in):: n
      real(dp), intent(in):: mMax, mMature

      call this%initSpectrum(n, mMin, mMax)

      allocate (this%nu(n))
      allocate (this%nupositive(n))
      allocate (this%Repro(n))
      allocate (this%g(n))

      allocate (this%psiMature(n))
      this%psiMature = (1.d0 + (this%m/mMature)**(-5.d0))**(-1.d0)  ! Maturity level

      this%beta = beta
      this%sigma = sigma
      this%epsAssim = epsAssim           ! Assimilation efficiency
      this%Cmax = h*(this%m**nn)         ! Maximum consumption rate
      this%V = gamma*this%m**q           ! Clearance rate
      this%metabolism = 0.2d0*this%Cmax  ! Standard metabolism

      this%mort0 = 0.1d0                 ! Intrinsic mortality
      this%mortF = 0.d0                  ! Fishing mortality

   end subroutine initFish

! calc Flux out and Flux in for each fish group
! u : the vector of state variables for each fish group
   subroutine calcfluxfish(this, u)
      class(spectrumFish), intent(inout):: this
      real(dp), intent(in):: u(this%n)
      real(dp)::kappa(this%n), ggamma(this%n)
      integer::i, ixPrev(this%n)

      this%nu = this%Eavail

!Flux out
      do i = 1, this%n
         this%nupositive(i) = max(0.d0, this%nu(i))

         kappa(i) = 1.d0 - this%psiMature(i)
         ggamma(i) = (kappa(i)*this%nupositive(i) - this%mort(i))/ &
                     (1.d0 - (1.d0/this%z(i))**(1.d0 - this%mort(i)/(kappa(i)*this%nupositive(i))))
         this%g(i) = kappa(i)*this%nupositive(i) !growth rate
         if (kappa(i) .eq. 0) then
            ggamma(i) = 0.d0 ! No growth of fully mature classes
         end if

         this%Jout(i) = ggamma(i)*u(i)                              ! Energy flow out of the current stage
         this%Repro(i) = (1 - kappa(i))*this%nupositive(i)*u(i)     ! Energy can be used for Reproduction
      end do

!Flux in
      ixPrev = [this%n, (i, i=1, (this%n - 1))]   ! index: the flux into the current grid is the flux out from the previous grid
      do i = 1, this%n
         this%Jin(i) = this%Jout(ixPrev(i))  ! the first grid will be overwrite in next step
      end do
      ! Reproduction:
      this%Jin(1) = epsRepro*(this%Jin(1) + sum(this%Repro))    ! overwrite first grid by reproduction data
   end subroutine calcfluxfish
!

end module fish
