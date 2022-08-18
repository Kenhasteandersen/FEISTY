!FEISTY
Module FEISTY
   use globals
   use input
   use spectrum
   use fish

   implicit none

   !integer, parameter :: idxR = 1 !index of resource 1:4

   ! Types parameter of fish species used in parametersAddGroup:
   !integer, parameter :: typeFish = 10
   integer, parameter :: typeResource = 100000000
   integer, parameter :: fishSmall = 1
   integer, parameter :: fishLarge = 10
   integer, parameter :: fishDemersal = 20

!spectrum
   integer:: nGroups       ! Number of fish spectrum groups
   integer:: iCurrentGroup ! The current group to be added
   integer:: nResources    ! Number of resource state variables
   integer:: idxF          ! First index into fish groups (=nResources+1)
   integer:: nGrid         ! Total number of grids including all kinds of Resources
   type(spectrumContainer), allocatable :: group(:)     ! Structure pointing to each group
   integer, dimension(:), allocatable :: ixStart, ixEnd ! Indices into u for each fish spectrum

   real(dp), allocatable:: theta(:, :)             ! Feeding preference matrix
   real(dp), dimension(:), allocatable:: upositive ! State variable must be positive
   real(dp), dimension(:), allocatable:: F         ! Available food

   real(dp) :: prod !primary production input from R

! assembled vectors including resources and fish
   real(dp), allocatable:: V(:)          ! assembled vector of clearance rate
   real(dp), allocatable:: Enc(:)        ! assembled vector of encounter rate
   real(dp), allocatable:: flvl(:)       ! assembled vector of feeding level
   real(dp), allocatable:: Cmax(:)       ! assembled vector of maximum consumption rate
   real(dp), allocatable:: mortpred(:)   ! assembled vector of predation mortality
   real(dp), allocatable:: mc(:)         ! assembled vector of central mass

! resource parameters    input from input file
   real(dp), allocatable:: K(:), rr(:)   ! Carrying capacity of resources and growth rate of resources
   real(dp) :: sbenk
   real(dp) :: lbenk
   real(dp) :: szoog
   real(dp) :: lzoog
   real(dp) :: sbeng
   real(dp) :: lbeng
!   real(dp) ::     sbenk =  5.d0               ! small benthos carry capacity
!   real(dp) ::     lbenk =  0.d0              ! large benthos carry capacity
!   real(dp) ::     szoog =  1.d0              ! small zooplankton growth rate
!   real(dp) ::     lzoog =  1.d0               ! large zooplankton growth rate
!   real(dp) ::     sbeng =  1.d0              ! small benthos growth rate
!   real(dp) ::     lbeng =  0.d0               ! large benthos growth rate

contains

   subroutine setupFEISTY(pprod)
      !integer::iFish
      real(dp), intent(in)::pprod

      call parametersInit(3, 2 + 3 + 3, 4, pprod)! (fish groups, total fish stages 2stages+3stages+3stages, 4 resources, pprod)
!      call parametersAddGroup(typeResource,) ! add resources

!      do
      call parametersAddGroup(fishSmall, 2, 2.50d2, 0.5d0) ! add fishes  !! Do you need these "fishSmall" etc variables
      call parametersAddGroup(fishLarge, 3, 1.25d5, 2.5d2)
      call parametersAddGroup(fishDemersal, 3, 1.25d5, 2.5d2)
!      end do

      call parametersFinalize()

   end subroutine setupFEISTY

! nnGroups: Fish group numbers
! nnGrid: all fish grids
   subroutine parametersInit(nnGroups, nnGrid, nnResources, pprod)
      integer, intent(in):: nnGroups, nnGrid, nnResources
      real(dp), intent(in)::pprod

      allocate (K(nResources))
      allocate (rr(nResources))

      nGroups = nnGroups                   ! fish size spectrum group numbers (species)
      iCurrentGroup = 0                    !
      nResources = nnResources             ! resource numbers
      nGrid = nnGrid + nnResources         ! total grid numbers   resources + total fish stages
      idxF = nResources + 1                ! fish grid begins at...
      prod = pprod                         ! primary production
      !
      ! Allocate variables:
      !
      if (allocated(upositive)) then
         deallocate (group)
         deallocate (ixStart)
         deallocate (ixEnd)
         deallocate (upositive)
         deallocate (F)
         deallocate (theta)
         !maybe more variables
         !
      end if

      allocate (group(nGroups))
      allocate (ixStart(nGroups))
      allocate (ixEnd(nGroups))
      allocate (upositive(nGrid))
      allocate (F(nGrid))
      allocate (theta(nGrid, nGrid))
      !

      call read_namelist_resources()    !load resources
      K = [prod, prod, sbenk, lbenk]    ! Carrying capacity of resources
      rr = [szoog, lzoog, sbeng, lbeng] ! growth rate of resources
   end subroutine parametersInit

!n:stages of a fish species
   subroutine parametersAddGroup(typeGroup, n, mMax, mMature)
      integer, intent(in) :: typeGroup, n      !number of stages
      real(dp), intent(in) :: mMax, mMature

      type(spectrumFish) :: specFish
      !
      ! Find the group number and grid location:
      !
      iCurrentGroup = iCurrentGroup + 1

      if (iCurrentGroup .eq. 1) then
         ixStart(iCurrentGroup) = idxF
      else
         ixStart(iCurrentGroup) = ixEnd(iCurrentGroup - 1) + 1
      end if

      ixEnd(iCurrentGroup) = ixStart(iCurrentGroup) + n - 1

      select case (typeGroup)
!case (typeResource)
!call initResource()
      case (fishSmall)
         call initFish(specFish, n, mMax, mMature, fishSmall)
         allocate (group(iCurrentGroup)%spec, source=specFish)

      case (fishLarge)
         call initFish(specFish, n, mMax, mMature, fishLarge)
         allocate (group(iCurrentGroup)%spec, source=specFish)

      case (fishDemersal)
         call initFish(specFish, n, mMax, mMature, fishDemersal)
         allocate (group(iCurrentGroup)%spec, source=specFish)
      end select

   end subroutine parametersAddGroup

   subroutine parametersFinalize()
      integer::iGroup, jGroup, i, j
        real(dp)::thetaS = 0.25d0 ! Medium fish pref for small zooplankton
        real(dp)::thetaA = 0.5d0  ! Large fish pref for medium forage fish
        real(dp)::thetaD = 0.75d0 ! Pref of large demersal on pelagic prey
!vectors:
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))

      theta = 0.d0 ! overwritten latter
      V = 0.d0
      Cmax = 0.d0

!Feeding preference matrix:    different from NUM theta & calcPhi
!assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formvector(spec, iGroup, V, Cmax, mc)
         end select
      end do
      mc(1:nResources) = [2d-06*sqrt(500d0), 1d-3*sqrt(500d0), 0.5d-03*sqrt(250000d0), 0.25d0*sqrt(500d0)] ! overwrite by resource mass
!basic feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            theta(i, j) = exp(-(log(mc(i)/(beta*mc(j))))**2/(2*sigma)**2)
            if (mc(j) .gt. mc(i)) then                                     !small can't eat large
               theta(i, j) = 0.d0
            end if
         end do
      end do

! further clac theta : feeding selection in terms of fish kinds and resources
      do iGroup = 1, nGroups

         select case (group(iGroup)%spec%typeGroup)

         case (fishSmall)

            theta(ixStart(iGroup):ixEnd(iGroup), 3:4) = 0.d0     ! pelagic has not feeding on benthos
            do jGroup = 1, nGroups
               if (group(jGroup)%spec%typeGroup .eq. fishDemersal) then
                  theta(ixStart(iGroup):ixEnd(iGroup), ixStart(jGroup):ixEnd(jGroup)) = 0.d0  !small pelagic has not feeding on demersals
               end if
            end do

         case (fishLarge)

            theta(ixStart(iGroup):ixEnd(iGroup), 3:4) = 0.d0            ! pelagic has not feeding on benthos

            do jGroup = 1, nGroups
               if (group(jGroup)%spec%typeGroup .eq. fishDemersal) then
                  do i = 1, group(jGroup)%spec%n
                  if (group(jGroup)%spec%m(i) > mMedium .AND. &
                      group(jGroup)%spec%m(i) < mLarge) then
                     theta(ixStart(iGroup):ixEnd(iGroup), ixStart(jGroup) + i - 1) = 0.d0 !large pelagic has not feeding on medium demersals
                  end if
                  end do
               end if
            end do

         case (fishDemersal)

            do i = 1, group(iGroup)%spec%n
               if (group(iGroup)%spec%m(i) < mMedium) then
                  theta(ixStart(iGroup) + i - 1, 3:4) = 0.d0        !Small demersals has not feeding on benthos
               end if

               if (group(iGroup)%spec%m(i) > mMedium .AND. &
                   group(iGroup)%spec%m(i) < mLarge) then
                  theta(ixStart(iGroup) + i - 1, 1:2) = 0.d0             ! Medium demersals has not feeding on zooplankton
                  theta(ixStart(iGroup) + i - 1, idxF:nGrid) = 0.d0      ! Medium demersals has not feeding on all fish (only eat benthos)
               end if
               ! Large demersals eat everything
            end do

         end select

      end do

!--------------------update---------------
theta=0.d0 ! reset
  ! Small pelagics:
  theta(5,1) = 1 ! Small ones eat only small zooplankton
  theta(6,1:10) = [thetaS, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0,0.d0, 1.d0]
  ! Large pelagics:
  theta(7,1) = 1
  theta(8,1:10) = [thetaS, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0,0.d0, 1.d0]
  theta(9,6) = thetaA ! medium forage fish
  theta(9,8) = 1 ! medium large pelagics
  ! Demersals:
  theta(10,1) = 1
  theta(11,3) = 1 ! Medium demersals eat benthos
  theta(12,6) = thetaA*thetaD  ! medium forage fish
  theta(12,8) = thetaD ! medium large pelagics
  theta(12,11) = 1 ! medium demersals






   end subroutine parametersFinalize

   subroutine formvector(this, iGroup, V, Cmax, mc) ! return assembled vectors containing values for all fish grid (no resources)
      integer, intent(in)::iGroup
      class(spectrumfish)::this
      real(dp), intent(out)::V(nGrid), Cmax(nGrid), mc(nGrid)

      V(ixStart(iGroup):ixEnd(iGroup)) = this%V
      Cmax(ixStart(iGroup):ixEnd(iGroup)) = this%Cmax
      mc(ixStart(iGroup):ixEnd(iGroup)) = this%m
   end subroutine formvector

!----------------------------------------------------------------------
!  Calculate the derivatives for all groups:
!  In:
!  u: the vector of state variables (all resources and all fish grids)
!  dudt: vector to hold the derivative (input and output)
!----------------------------------------------------------------------
   subroutine calcderivatives(u, dudt)
      real(dp), intent(in) :: u(nGrid)
      real(dp), intent(inout) :: dudt(nGrid)
      integer :: i, j, iGroup!, jGroup, ixi, ixj

      dudt = 0.d0 ! overwritten latter

      do i = 1, nGrid
         upositive(i) = max(0.d0, u(i)) ! keep state variables positive
      end do

      Enc = V*matmul(theta, upositive)  !a vector contain all Enc for all grid (resources+fish)
      flvl = Enc/(Cmax + Enc)           !a vector contain all flvl for all grid (resources+fish)

      do i = 1, nGrid
         if (isnan(flvl(i))) then
            flvl(i) = 0.d0
         end if
      end do
! Calc available food:
      do i = 1, nGrid                               !all components  R and Fish
         F(i) = 0.d0                             ! overwritten latter
         do j = 1, nGrid
            F(i) = F(i) + theta(i, j)*upositive(j) !matrix multiplication
         end do
      end do
!F=matmul(theta,upositive) !totally same as above

!Calc Encounter/feeding:
!return Enc flvl  Eavail for each group
      do iGroup = 1, nGroups
         call calcFeeding(group(iGroup)%spec, F(ixStart(iGroup):ixEnd(iGroup)))
      end do

! Mortality:
! Predation mortality, including all resources and fish grids    (different from NUM subroutine calcDerivativesUnicellulars)
      mortpred = matmul(transpose(theta), (flvl*Cmax/epsAssim*upositive/mc))

! Total mortality, only for each fish group
      do iGroup = 1, nGroups
         do i = 1, (ixEnd(iGroup) - ixStart(iGroup) + 1)
            group(iGroup)%spec%mort(i) = mortpred(ixStart(iGroup) + i - 1) + &  !predation mortality
                                         group(iGroup)%spec%mort0(i) + &        !intrinsic mortality
                                         group(iGroup)%spec%mortF(i)            !fishing mortality
         end do
      end do

! Flux and derivative assemble
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call calcfluxfish(spec, upositive(ixStart(iGroup):ixEnd(iGroup)))                          ! Flux out and Flux in
            call calcderiv(spec, upositive(1:nResources), upositive(ixStart(iGroup):ixEnd(iGroup)), &  ! Assemble derivatives
                           dudt(1:nResources), dudt(ixStart(iGroup):ixEnd(iGroup)))
         end select
      end do

   contains
!-------------------------------------------------------------------------
!  Assemble derivatives of resources and fish:
!  In:
!  R: the vector of state variables (all resources)
!  B: the vector of state variables (all fish grids)
!  dRdt: vector to hold the derivative of resources (input and output)
!  dBdt: vector to hold the derivative of fish (input and output)
!-------------------------------------------------------------------------
      subroutine calcderiv(this, R, B, dRdt, dBdt)
         class(spectrumfish)::this
         real(dp)::R(nResources), B(this%n)
         real(dp), intent(inout)::dRdt(nResources), dBdt(this%n)
         integer::i, j

!Fish:
         do i = 1, this%n
            dBdt(i) = this%Jin(i) - this%Jout(i) + (this%nu(i) - this%mort(i))*B(i) - this%Repro(i)
         end do
!Resources:
         do j = 1, nResources
            dRdt(j) = rr(j)*(K(j) - R(j)) - mortpred(j)*R(j)
         end do
      end subroutine calcderiv

   end subroutine calcderivatives
!----------------------------------------------------------------------------------------------------------
! Get rates from last step of calcderivatives
! in:
!      u: all state variables from the last time step
! out:
!      flvl_r vector of feeding level (resources + fish) resources are 0
!      mortpred_r vector of predation mortality (resources + fish)
!      g_r vector of growth rate (only fish grids)
! inout:    dudt: a vector holding all derivatives , probably not used
   subroutine getrates(u, dudt, flvl_r, mortpred_r, g_r)
      real(dp), intent(in) :: u(nGrid)
      real(dp), intent(inout):: dudt(nGrid) ! because dudt is set as intent(inout) in calcderivatives, overwritten in calcderivatives
      real(dp), intent(out) :: flvl_r(nGrid)
      real(dp), intent(out) :: mortpred_r(nGrid)
      real(dp), intent(out) :: g_r(nGrid - nResources)
      integer::iGroup

      call calcderivatives(u, dudt)
      flvl_r = flvl
      mortpred_r = mortpred

! output growth rate vector (only fish)
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call calcfluxfish(spec, u(ixStart(iGroup):ixEnd(iGroup)))   !get g of each fish group
            call assembleg(spec, iGroup, g_r)                           !put g in a vector group by group
         end select
      end do
   end subroutine getrates

! used for getrates: assemble growth rate vector for all fish
   subroutine assembleg(this, iGroup, g_r)
      class(spectrumFish) :: this
      integer, intent(in) :: iGroup
      real(dp), intent(out) :: g_r(nGrid - nResources)

      g_r(ixStart(iGroup) - nResources:ixEnd(iGroup) - nResources) = this%g
   end subroutine assembleg

   subroutine read_namelist_resources()
      integer :: file_unit, io_err

      namelist /input_resources/ sbenk, lbenk,&
                                  &szoog, lzoog, sbeng, lbeng

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_resources, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_resources

end module FEISTY
