! Handle fish groups

module fish
   use globals
   use spectrum

   implicit none

   private
! fish physiology:
   real(rk) :: h
   real(rk) :: hCepha   ! Note in input.nml
   real(rk) :: nn
   real(rk) :: q
   real(rk) :: gamma
   real(rk) :: kk
   real(rk) :: p
   real(rk) :: epsAssim
   real(rk) :: epsRepro
   real(rk) :: epst
!   real(rk), parameter ::  h = 20._rk            ! Max. consumption coefficient
!   real(rk), parameter ::  nn = -0.25_rk         ! Max. consumption exponent
!   real(rk), parameter ::  q = -0.2_rk           ! Clearance rate exponent
!   real(rk), parameter ::  gamma = 70._rk        ! Coef. for clearance rate
!   real(rk), parameter ::  kk = 0.011_rk*365._rk  ! Metabolism coefficient   4.0150 in input file
!   real(rk), parameter ::  p = -0.175_rk         ! Metabolism exponent
!   real(rk), parameter ::  epsAssim = 0.7_rk     ! Assimilation efficiency
!   real(rk), parameter ::  epsRepro = 0.01_rk    ! reproduction * recruitment efficiency
!   real(rk), parameter ::  epst = 0.1_rk         ! efficiency of benthos community
!   real(rk), parameter ::  hCepha = 28._rk / (epsAssim * (0.6_rk - 0.4_rk)) ! Max. consumption coefficient of squid

! for size spectrum
   real(rk) :: beta
   real(rk) :: betaCepha
   real(rk) :: sigma
   real(rk) :: mMin
!   real(rk), parameter ::  beta = 400._rk
!   real(rk), parameter ::  betaCepha = 50._rk ! squid
!   real(rk), parameter ::  sigma = 1.3_rk
!   real(rk) ::  mMin = 0.001_rk       ! min fish mass (boundary of the grid) used for discretization

! for feeding preference calc
   real(rk) :: mMedium
   real(rk) :: mLarge
!   real(rk), parameter ::  mMedium = 10._rk      !meidum fish central mass used for feeding preference calc
!   real(rk), parameter ::  mLarge = 5000._rk     !large fish central mass used for feeding preference calc

!resources:
   real(rk) :: lbenk
   real(rk) :: szoog
   real(rk) :: lzoog
   real(rk) :: sbeng
   real(rk) :: lbeng
!   real(rk), parameter ::     lbenk = 0._rk              ! large benthos carry capacity
!   real(rk), parameter ::     szoog = 1._rk              ! small zooplankton growth rate
!   real(rk), parameter ::     lzoog = 1._rk              ! large zooplankton growth rate
!   real(rk), parameter ::     sbeng = 1._rk              ! small benthos growth rate
!   real(rk), parameter ::     lbeng = 0._rk              ! large benthos growth rate

   type, extends(typeSpectrum) :: spectrumFish

   contains
      procedure, pass :: initFish
      !procedure :: calcfluxfish
   end type spectrumFish

   public h, hCepha , nn, q, gamma, kk, p, epsAssim, epsRepro, epst, beta, betaCepha, sigma, mMin, mMedium, mLarge,  &
          spectrumFish, &
          initFish,&! calcfluxfish, &
          lbenk, szoog, lzoog, sbeng, lbeng

contains

! initialize fish groups
   subroutine initFish(this, n, mMax, mMature)
      class(spectrumFish), intent(inout) :: this
      integer, intent(in):: n
      real(rk), intent(in):: mMax, mMature

      call this%initSpectrum(n, mMin, mMax) ! create mass spectrum

      allocate (this%nu(n))
      allocate (this%nupositive(n))
      allocate (this%Repro(n))
      allocate (this%g(n))

      allocate (this%psiMature(n))
      this%psiMature = (1._rk + (this%m/mMature)**(-5._rk))**(-1._rk)  ! Maturity level

      this%beta = beta
      this%sigma = sigma
      this%epsAssim = epsAssim           ! Assimilation efficiency
      this%Cmax = h*(this%m**nn)         ! Maximum consumption rate
      this%V = gamma*this%m**q           ! Clearance rate
      this%metabolism = 0.2_rk*this%Cmax  ! Standard metabolism

      this%mort0 = 0.1_rk                 ! Intrinsic mortality
      this%mortF = 0._rk                  ! Fishing mortality

      this%Cmaxsave = this%Cmax
      this%Vsave = this%V
      this%metabolismsave = this%metabolism

   end subroutine initFish

end module fish
