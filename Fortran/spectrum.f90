! Handle spectrum building

module spectrum
   use globals

   implicit none

   type, abstract :: typeSpectrum
      integer:: n
! Grid:
      real(dp), allocatable:: m(:)  ! Geometric center mass of size-bin
      real(dp), allocatable:: mLower(:)  ! Smallest size in the bin
      real(dp), allocatable:: mUpper(:)   ! Width of the bin
      real(dp), allocatable:: z(:) ! Ratio btw upper and lower size of bin

      real(dp), allocatable::psiMature(:)          ! Maturity level indicator
!feeding
      real(dp):: beta, sigma                       ! Pred:prey mass ratio and width
      real(dp):: epsAssim                          ! Assimilation efficiency
      real(dp), allocatable:: Enc(:)               ! Encounter
      real(dp), allocatable:: flvl(:)              ! Feeding level
      real(dp), allocatable :: Cmax(:)             ! Maximum consumption rate
      real(dp), allocatable ::  V(:)               ! Clearance rate
      real(dp), allocatable  :: metabolism(:)      ! Standard metabolism
      real(dp), allocatable  :: Eavail(:)          ! Available energy

!mortality
      real(dp), dimension(:), allocatable:: mort, mortpred, mort0, mortF
! total mortality=predation mortality+intrinsic mortality+fishing mortality :  mort=mortpred+mort0+mortF
!flux
      real(dp), dimension(:), allocatable:: Jin, Jout                 !
      real(dp), dimension(:), allocatable:: nu, nupositive, Repro, g  !nu = Eavail

   contains

      procedure, pass :: initSpectrum
      procedure :: calcFeeding
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
      real(dp), intent(in):: mMin, mMax

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

   contains
!
! Set up grids in terms of given minimum and maximum boundary masses
!
      subroutine calcGrid()
         real(dp):: mb(n + 1)

         mb = exp(linspace(log(mMin), log(mMax), n + 1))        ! boundary mass   grid of lower sizes, with the last being the upper size of the last cell
         this%mLower = mb(1:n)                                  ! lower boundary
         this%mUpper = mb(2:(n + 1))                            ! upper boundary
         this%z = this%mUpper/this%mLower                       ! The ratio between upper and lower sizes
         this%m = exp(log(this%mLower) + 0.5d0*(log(this%z)))   ! Geometric mean center mass

      end subroutine calcGrid

   end subroutine initSpectrum

! calc for each fish group
! in : available food of each group
! out:  Enc flvl Eavail for each fish group
   subroutine calcFeeding(this, F)
      class(typeSpectrum), intent(inout):: this
      real(dp), intent(in):: F(this%n)

      this%Enc = this%V*F                                ! Encounter rate
      this%flvl = this%Enc/(this%Cmax + this%Enc)        ! feeding level
      this%Eavail = this%epsAssim*this%flvl*this%Cmax - this%metabolism ! available energy

   end subroutine calcFeeding

end module spectrum
