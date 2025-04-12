! Handle spectrum building

module spectrum
   use globals

   implicit none

   type, abstract :: typeSpectrum
      integer:: n
! Grid:
      real(rk), allocatable:: m(:)  ! Geometric center mass of size-bin
      real(rk), allocatable:: mLower(:)  ! Smallest size in the bin
      real(rk), allocatable:: mUpper(:)   ! Width of the bin
      real(rk), allocatable:: z(:) ! Ratio btw upper and lower size of bin

      real(rk), allocatable::psiMature(:)          ! Maturity level indicator
!feeding
      real(rk):: beta, sigma                       ! Pred:prey mass ratio and width
      real(rk):: epsAssim                          ! Assimilation efficiency
      real(rk), allocatable:: Enc(:)               ! Encounter
      real(rk), allocatable:: flvl(:)              ! Feeding level
      real(rk), allocatable :: Cmax(:)             ! Maximum consumption rate
      real(rk), allocatable ::  V(:)               ! Clearance rate
      real(rk), allocatable  :: metabolism(:)      ! Standard metabolism
      real(rk), allocatable  :: Eavail(:)          ! Available energy

!mortality
      real(rk), dimension(:), allocatable:: mort, mortpred, mort0, mortF
! total mortality=predation mortality+intrinsic mortality+fishing mortality :  mort=mortpred+mort0+mortF
!flux
      real(rk), dimension(:), allocatable:: Jin, Jout                 !
      real(rk), dimension(:), allocatable:: nu, nupositive, Repro, g  !nu = Eavail

      ! save for temperature (Feb 2024 added)
      real(rk), allocatable :: Cmaxsave(:)             !
      real(rk), allocatable :: Vsave(:)                !
      real(rk), allocatable :: metabolismsave(:)      !

   contains

      procedure, pass :: initSpectrum
      !procedure :: calcFeeding
   end type typeSpectrum

   type spectrumContainer
      class(typeSpectrum), allocatable :: spec
   end type spectrumContainer

!---------------------------------------------------------------------------
contains
! create mass spectrum
   subroutine initSpectrum(this, n, mMin, mMax)
      class(typeSpectrum) :: this
      integer, intent(in):: n
      real(rk), intent(in):: mMin, mMax

      this%n = n
      allocate (this%m(n))
      allocate (this%mLower(n))
      allocate (this%mUpper(n))
      allocate (this%z(n))

      call calcGrid()

      allocate (this%Enc(n))
      allocate (this%flvl(n))

      allocate (this%Cmax(n))
      allocate (this%V(n))
      allocate (this%metabolism(n))
      allocate (this%Eavail(n))

      allocate (this%mort(n))
      allocate (this%mortpred(n))
      allocate (this%mort0(n))
      allocate (this%mortF(n))

      allocate (this%Jin(n))
      allocate (this%Jout(n))

      allocate (this%Cmaxsave(n))
      allocate (this%Vsave(n))
      allocate (this%metabolismsave(n))

   contains
!
! Set up grids in terms of given minimum and maximum boundary masses
!
      subroutine calcGrid()
         real(rk):: mb(n + 1)

         mb = exp(linspace(log(mMin), log(mMax), n + 1))        ! boundary mass   grid of lower sizes, with the last being the upper size of the last cell
         this%mLower = mb(1:n)                                  ! lower boundary
         this%mUpper = mb(2:(n + 1))                            ! upper boundary
         this%z = this%mUpper/this%mLower                       ! The ratio between upper and lower sizes
         this%m = exp(log(this%mLower) + 0.5_rk*(log(this%z)))   ! Geometric mean center mass

      end subroutine calcGrid

   end subroutine initSpectrum


end module spectrum
