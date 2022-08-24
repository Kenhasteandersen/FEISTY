!FEISTY
Module FEISTY
   use globals
   use input
   use spectrum
   use fish

   implicit none

   !integer, parameter :: idxR = 1 !index of resource 1:4

!spectrum
   integer:: nGroups       ! Number of fish spectrum groups
   integer:: iCurrentGroup ! The current group to be added
   integer:: nResources    ! Number of resource state variables
   integer:: idxF          ! First index into fish groups (=nResources+1)
   integer:: nGrid         ! Total number of grids including all kinds of Resources and fish
   type(spectrumContainer), allocatable :: group(:)     ! Structure pointing to each group
   integer, dimension(:), allocatable :: ixStart, ixEnd ! Indices into u for each fish spectrum

   real(dp), allocatable:: sizeprefer(:, :)        ! Size preference matrix
   real(dp), allocatable:: vertover(:, :)          ! Vertical overlap matrix
   real(dp), allocatable:: theta(:, :)             ! Feeding preference matrix
   real(dp), dimension(:), allocatable :: upositive ! State variable must be positive
   real(dp), dimension(:), allocatable :: F         ! Available food

! assembled vectors including resources and fish
   real(dp), allocatable:: V(:)          ! assembled vector of clearance rate
   real(dp), allocatable:: Enc(:)        ! assembled vector of encounter rate
   real(dp), allocatable:: flvl(:)       ! assembled vector of feeding level
   real(dp), allocatable:: Cmax(:)       ! assembled vector of maximum consumption rate
   real(dp), allocatable:: mortpred(:)   ! assembled vector of predation mortality
   real(dp), allocatable:: mc(:)         ! assembled vector of central mass
   real(dp), allocatable:: mL(:)         ! assembled vector of lower mass
   real(dp), allocatable:: mU(:)         ! assembled vector of upper mass

! resource parameters    input from input file
   real(dp), allocatable:: K(:), rr(:)   ! Carrying capacity of resources and growth rate of resources
!   real(dp) :: sbenk
!   real(dp) :: lbenk
!   real(dp) :: szoog
!   real(dp) :: lzoog
!   real(dp) :: sbeng
!   real(dp) :: lbeng
!   real(dp), parameter ::     sbenk = 5.d0               ! small benthos carry capacity
!   real(dp), parameter ::     lbenk = 0.d0              ! large benthos carry capacity
!   real(dp), parameter ::     szoog = 1.d0              ! small zooplankton growth rate
!   real(dp), parameter ::     lzoog = 1.d0               ! large zooplankton growth rate
!   real(dp), parameter ::     sbeng = 1.d0              ! small benthos growth rate
!   real(dp), parameter ::     lbeng = 0.d0               ! large benthos growth rate

! van Denderen et al., 2020

   real(dp), parameter ::     lbenk = 0.d0
   real(dp), parameter ::     szoog = 1.d0
   real(dp), parameter ::     lzoog = 1.d0
   real(dp), parameter ::     sbeng = 1.d0
   real(dp), parameter ::     lbeng = 0.d0

contains

   subroutine setupbasic(pprod, bprod) ! Petrik et al., 2019
      real(dp), intent(in)::pprod, bprod
      integer :: iGroup ! , i, j
      real(dp),parameter ::   thetaS = 0.25d0 ! Medium fish pref for small zooplankton
      real(dp),parameter ::   thetaA = 0.5d0  ! Large fish pref for medium forage fish
      real(dp),parameter ::   thetaD = 0.75d0 ! Pref of large demersal on pelagic prey

      ! call read_namelist_resources()    !load resources
      call parametersInit(3, 2 + 3 + 3, 4, pprod, bprod)! (fish groups, total fish stages 2stages+3stages+3stages, 4 resources, pprod)
      call parametersAddGroup(2, 2.50d2, 0.5d0) ! fishSmall,
      call parametersAddGroup(3, 1.25d5, 2.5d2) ! fishLarge,
      call parametersAddGroup(3, 1.25d5, 2.5d2)   ! fishDemersal,

      ! vectors:
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))
      allocate (mL(nGrid))
      allocate (mU(nGrid))

      theta = 0.d0 ! overwritten latter
      V = 0.d0
      Cmax = 0.d0

! Overwrite
      do iGroup = 1, nGroups
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p) ! overwrite matabolism

         group(iGroup)%spec%psiMature = 0.d0 ! reset
         group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5d0 ! only adults reproduce

         if (group(iGroup)%spec%n .eq. 2) then
            group(iGroup)%spec%mortF(group(iGroup)%spec%n) = 0.3d0 ! only adults have fishing mortality
         else if ((group(iGroup)%spec%n .eq. 3)) then
            group(iGroup)%spec%mortF(group(iGroup)%spec%n - 1) = 0.03d0 ! juvenile
            group(iGroup)%spec%mortF(group(iGroup)%spec%n) = 0.3d0 ! adult
         end if
      end do

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formvector(spec, iGroup, V, Cmax, mc, mL, mU)
         end select
      end do
      mc(1:nResources) = [2.d-06*sqrt(500.d0), 1.d-3*sqrt(500.d0), 0.5d-03*sqrt(250000.d0), 0.25d0*sqrt(500.d0)] ! overwrite by resource mass
      !mU(1:nResources) = [2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)] ! weight central size
      !mL(1:nResources) = [2e-06,0.001, 0.5e-03, 0.25] ! weight lower limit)

      ! Small pelagics:
      theta(5, 1) = 1.d0 ! Small ones eat only small zooplankton
      theta(6, 1:10) = [thetaS, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0]

      ! Large pelagics:
      theta(7, 1) = 1.d0
      theta(8, 1:10) = [thetaS, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0]
      theta(9, 6) = thetaA          ! medium forage fish
      theta(9, 8) = 1.d0            ! medium large pelagics

      ! Demersals:
      theta(10, 1) = 1.d0           ! larval demersals ear small zooplankton
      theta(11, 3) = 1.d0           ! medium demersals eat small benthos
      theta(12, 6) = thetaA*thetaD  ! medium forage fish
      theta(12, 8) = thetaD         ! medium large pelagics
      theta(12, 11) = 1.d0          ! medium demersals

   end subroutine setupbasic
!------------------------------------------------------------------------------------------------------------------------------

   subroutine setupbasic2(pprod, bprod, nStages) !
      real(dp), intent(in) :: pprod, bprod
      integer, intent(in) :: nStages
      integer :: iGroup, i, j
      real(dp), parameter :: thetaA = 0.5d0  ! Large fish pref for medium forage fish
      real(dp), parameter :: thetaD = 0.75d0 ! Pref of large demersal on pelagic prey

      ! call read_namelist_resources()    !load resources
      call parametersInit(3, nint(0.66d0*nStages) + nStages + nStages, 4, pprod, bprod)! (fish groups, total fish stages 2stages+3stages+3stages, 4 resources, pprod)
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, 0.5d0) ! fishSmall    nint: Returns the nearest integer to the argument.
      call parametersAddGroup(nStages, 1.25d5, 2.5d2) ! fishLarge (stages, max mass, mature mass)
      call parametersAddGroup(nStages, 1.25d5, 2.5d2) ! fishDemersal

      ! vectors:
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))
      allocate (mL(nGrid))
      allocate (mU(nGrid))

      theta = 0.d0 ! overwritten latter
      V = 0.d0
      Cmax = 0.d0

! Overwrite
      do iGroup = 1, nGroups
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p)!overwrite matabolism
      end do

! Feeding preference matrix:    different from NUM theta & calcPhi
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formvector(spec, iGroup, V, Cmax, mc, mL, mU)
         end select
      end do
      mc(1:nResources) = [2.d-06*sqrt(500.d0), 1.d-3*sqrt(500.d0), 0.5d-03*sqrt(250000.d0), 0.25d0*sqrt(500.d0)] ! overwrite by resource mass
      ! mU = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) ! weight central size
      ! mL = c(2e-06,0.001, 0.5e-03, 0.25) ! weight lower limit)

!basic feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            theta(i, j) = exp(-(log(mc(i)/(beta*mc(j))))**2/(2*sigma)**2)
            if (mc(j) .gt. mc(i)) then                   !small can't eat large
               theta(i, j) = 0.d0
            end if
         end do
      end do

! further clac theta : feeding selection in terms of fish kinds and resources
      ! Small pelagic
      theta(ixStart(1):ixEnd(1), 3:4) = 0.d0                  ! pelagic has not feeding on benthos
      do i = 1, group(3)%spec%n
      if (group(3)%spec%m(i) > mMedium .AND. &
          group(3)%spec%m(i) < mLarge) then
         theta(ixStart(1):ixEnd(1), ixStart(3) + i - 1) = 0.d0 !small pelagic has not feeding on medium demersals
      end if
      end do
      ! Large pelagic
      theta(ixStart(2):ixEnd(2), 3:4) = 0.d0                  ! pelagic has not feeding on benthos
      do i = 1, group(3)%spec%n
      if (group(3)%spec%m(i) > mMedium .AND. &
          group(3)%spec%m(i) < mLarge) then
         theta(ixStart(2):ixEnd(2), ixStart(3) + i - 1) = 0.d0 !large pelagic has not feeding on medium demersals
      end if
      end do
      ! Large pelagics have reduced feeding efficiency on small pelagics:
      theta(ixStart(2):ixEnd(2), ixStart(1):ixEnd(1)) = thetaA*theta(ixStart(2):ixEnd(2), ixStart(1):ixEnd(1))
      ! Demersals
      do i = 1, group(3)%spec%n
!         if (group(3)%spec%m(i) < mMedium) then
!            theta(ixStart(3) + i - 1, 3:4) = 0.d0             !Small demersals has not feeding on benthos
!         end if
!         if (group(3)%spec%m(i) > mMedium .AND. &
!             group(3)%spec%m(i) < mLarge) then
!            theta(ixStart(3) + i - 1, 1:2) = 0.d0             ! Medium demersals has not feeding on zooplankton
!            theta(ixStart(3) + i - 1, idxF:nGrid) = 0.d0      ! Medium demersals has not feeding on all fish (only eat benthos)
!         end if
         if (group(3)%spec%m(i) > mMedium .AND. &
             group(3)%spec%m(i) < mLarge) then
            theta(ixStart(3) + i - 1, ixStart(1):ixEnd(1)) = 0.d0    ! perhaps not correct ?
            theta(ixStart(3) + i - 1, ixStart(2):ixEnd(2)) = 0.d0    !
         end if

      end do
      ! Large demersals feed have reduced feeding effiiency on pelagic species:
      theta(ixStart(3):ixEnd(3), ixStart(1):ixEnd(1)) = thetaA*thetaD*theta(ixStart(3):ixEnd(3), ixStart(1):ixEnd(1))
      theta(ixStart(3):ixEnd(3), ixStart(2):ixEnd(2)) = thetaD*theta(ixStart(3):ixEnd(3), ixStart(2):ixEnd(2))

   end subroutine setupbasic2
!------------------------------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------------------------------
   subroutine setupVertical(pprod) ! Vertical overlap from MATLAB  van Denderen et al., 2020
      real(dp), intent(in) :: pprod !

! for theta calc
      real(dp) :: ssigma = 10.d0
      real(dp) :: tau = 10.d0
      real(dp), allocatable :: sigmap(:) ! width for each size class
      real(dp), parameter :: bottom = 1500.d0 ! total depth meter
      real(dp), parameter :: photic = 150.d0  ! photic zone depth
      real(dp), parameter :: mesop = 250.d0   ! depth ?
      real(dp), parameter :: bent = 150.d0  ! ?
      real(dp) :: bprod
      real(dp), parameter :: visual = 1.5d0 ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
      real(dp), dimension(int(bottom) + 1) :: xrange
      real(dp) :: dvm  ! vertical migration depth 650
      real(dp) :: xloc ! vertical location    will be overwritten again and again
      real(dp), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(dp), allocatable :: zp_n(:, :), zp_d(:, :), & ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :)        ! demersal night/ day
      real(dp) :: demmig ! ?
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(dp), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
      integer :: iGroup, i, j, ixjuv, ixadult

      bprod = 0.1d0*(bent*(bottom/photic)**(-0.86d0)) ! from matlab
      if (bprod .ge. bent*0.1d0) then
         bprod = bent*0.1d0
      end if

      call parametersInit(5, 2 + 2 + 3 + 3 + 3, 4, pprod, bprod)!
      call parametersAddGroup(2, 2.50d2, 0.5d0) ! fishSmall,
      call parametersAddGroup(2, 2.50d2, 0.5d0) ! fishMeso,
      call parametersAddGroup(3, 1.25d5, 2.5d2) ! fishLarge,
      call parametersAddGroup(3, 1.25d5, 2.5d2) ! fishBathy,
      call parametersAddGroup(3, 1.25d5, 2.5d2) ! fishDemersal,

! vectors and matrix:
      allocate (sigmap(nGrid))
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))
      allocate (mL(nGrid))
      allocate (mU(nGrid))
      allocate (vertover(nGrid, nGrid)) !

      vertover = 0.d0
      sizeprefer = 0.d0
      theta = 0.d0 ! overwritten latter
      V = 0.d0
      Cmax = 0.d0

! Overwrite
      do iGroup = 1, nGroups

         group(iGroup)%spec%metabolism = (0.2d0*h*group(iGroup)%spec%m**p)

         group(iGroup)%spec%psiMature = 0.d0 ! reset
         group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5d0! only adults reproduce

      end do
      group(nGroups)%spec%mortF(group(nGroups)%spec%n) = 0.5d0 ! only demersal adults have fishing mortality

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formvector(spec, iGroup, V, Cmax, mc, mL, mU)
         end select
      end do
      !from baseparameters.m
      mc(1:nResources) = [2.d-06*sqrt(500.d0), 1.d-3*sqrt(500.d0), 0.5d-03*sqrt(250000.d0), 0.25d0*sqrt(500.d0)] ! resource central mass
      mU(1:nResources) = [0.001d0, 0.5d0, 125.d0, 125.d0]  ! resource mass upper limit
      mL(1:nResources) = [2.d-6, 0.001d0, 0.5d-3, 0.25d0] ! resource mass lower limit
! basic feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            sizeprefer(i, j) = sqrt(pi/2.d0)*sigma*( &
                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)) &
                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)))
            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
         end do
      end do
!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0.d0, bottom, int(bottom) + 1)
      dvm = photic + 500.d0 ! 650.d0
      !  from matlab
      if (bottom .lt. (photic + 500.d0)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      else if (bottom .le. mesop) then
         dvm = 0.d0                   ! no migration in shallow habitats
      else
      end if

! first stages as juvenile/adult for predators
      ixjuv = 2     !minloc(abs(sizes-smat)); from matlab
      ixadult = 3   !minloc(abs(sizes-lmat));

! zooplankton night
      allocate (zp_n(size(xrange), 2))
      xloc = 0.d0 ! zoo on surface at night
      do i = 1, 2 !small zoo & large zoo
         zp_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(i)**2.d0)))* &
                      exp(-(((xrange - xloc)**2.d0)/(2.d0*sigmap(i)**2.d0)))
      end do
      zp_n = matmul(zp_n, diag(1.d0/sum(zp_n, 1)))

! zooplankton day (half at surface, half at dvm depth
      allocate (zp_d(size(xrange), 2))
      xloc = dvm !
      do i = 1, 2 !index: small zoo & large zoo
         zp_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(i)**2.d0)))* &
                      exp(-(((xrange - xloc)**2.d0)/(2.d0*sigmap(i)**2.d0)))
      end do
      zp_d = matmul(zp_d, diag(1.d0/sum(zp_d, 1)))
      zp_d = (zp_n + zp_d)/2.d0

! benthos small and large (at bottom with width sigma)
      allocate (bent_dn(size(xrange), 2))
      xloc = bottom
      do i = 1, 2 ! small benthos & large benthos
         bent_dn(:, i) = (1.d0/(sqrt(2.d0*pi*ssigma**2.d0)))*exp(-((xrange - xloc)**2.d0/(2.d0*ssigma**2.d0)))
         bent_dn(:, i) = bent_dn(:, i)/sum(bent_dn(:, i))
      end do

! small pelagic fish (day + night) always at surface
      allocate (ix(ixEnd(1) - ixStart(1) + 1))
      allocate (spel_dn(size(xrange), ixEnd(1) - ixStart(1) + 1))
      xloc = 0.d0
      ix = [(i, i=ixStart(1), ixEnd(1))]
      do i = 1, size(ix)
         spel_dn(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                         exp(-((xrange - xloc)**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      spel_dn = matmul(spel_dn, diag(1.d0/sum(spel_dn, 1)))

! meso pelagic night   at surface
      allocate (mpel_n(size(xrange), ixEnd(2) - ixStart(2) + 1))
      deallocate (ix)
      allocate (ix(ixEnd(2) - ixStart(2) + 1))
      mpel_n = spel_dn

! meso pelagic day (all at dvm)
      allocate (mpel_d(size(xrange), ixEnd(2) - ixStart(2) + 1))
      xloc = dvm
      ix = [(i, i=ixStart(2), ixEnd(2))]
      do i = 1, size(ix)
         mpel_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xloc)**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      mpel_d = matmul(mpel_d, diag(1.d0/sum(mpel_d, 1)))

! large pelagic fish night (all at surface)
      allocate (lpel_n(size(xrange), ixEnd(3) - ixStart(3) + 1))
      deallocate (ix)
      allocate (ix(ixEnd(3) - ixStart(3) + 1))
      xloc = 0.d0
      ix = [(i, i=ixStart(3), ixEnd(3))]
      do i = 1, size(ix)
         lpel_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xloc)**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      lpel_n = matmul(lpel_n, diag(1.d0/sum(lpel_n, 1)))

! large pelagic fish day (non-adult at surface   adult at dvm)
      allocate (lpel_d(size(xrange), ixEnd(3) - ixStart(3) + 1))
      allocate (xlocvec(ixEnd(3) - ixStart(3) + 1))
      xlocvec = 0.d0 ! initialization
      xlocvec(ixadult:size(xlocvec)) = dvm  !  non-adult at surface   adult at dvm
      do i = 1, size(ix)
         lpel_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      lpel_d = matmul(lpel_d, diag(1.d0/sum(lpel_d, 1)))
      lpel_d = (lpel_d + lpel_n)/2.d0

! bathypelagic night (adults in midwater, others at surface)
      allocate (bpel_n(size(xrange), ixEnd(4) - ixStart(4) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(4) - ixStart(4) + 1))
      allocate (ix(ixEnd(4) - ixStart(4) + 1))
      xlocvec = 0.d0 ! initialization
      xlocvec(ixadult:size(xlocvec)) = dvm  !  non-adult at surface   adult at dvm
      ix = [(i, i=ixStart(4), ixEnd(4))]
      do i = 1, size(ix)
         bpel_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      bpel_n = matmul(bpel_n, diag(1.d0/sum(bpel_n, 1)))

! bathypelagic day (all at dvm)
      allocate (bpel_d(size(xrange), ixEnd(4) - ixStart(4) + 1))
      xlocvec = dvm ! overwrite all elements by dvm
      do i = 1, size(ix)
         bpel_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      bpel_d = matmul(bpel_d, diag(1.d0/sum(bpel_d, 1)))

! demersal fish night
      allocate (dem_n(size(xrange), ixEnd(5) - ixStart(5) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(5) - ixStart(5) + 1))
      allocate (ix(ixEnd(5) - ixStart(5) + 1))
      xlocvec = 0.d0 ! initialization
      xlocvec(ixjuv:size(xlocvec)) = bottom  !  larvae at surface   juvenile and adult at bottom
      ix = [(i, i=ixStart(5), ixEnd(5))]
      do i = 1, size(ix)
         dem_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_n = matmul(dem_n, diag(1.d0/sum(dem_n, 1)))

! demersal fish day
      demmig = dvm ! ??? from matlab
      if ((bottom - dvm) .ge. 1200.d0) then
         demmig = dvm + (bottom - dvm - 1200.d0)
      else if ((bottom - dvm) .ge. 1500.d0) then
         demmig = bottom
      else
      end if
      allocate (dem_d(size(xrange), ixEnd(5) - ixStart(5) + 1))
      xlocvec(ixadult:size(xlocvec)) = dvm ! larvae at surface/ juvenile at bottom/ adult and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_d = matmul(dem_d, diag(1.d0/sum(dem_d, 1)))
! ?? from matlab
      ! if shallower than euphotic depth, adult demersals feed across-habitats
      if (bottom .le. photic) then
         dem_d = (dem_d + dem_n)/2.d0
         dem_n = dem_d
      end if

! calculate overlap during day
      allocate (depthDay(size(xrange), nGrid))
      allocate (test(size(xrange), nGrid))
      allocate (dayout(nGrid, nGrid))
      depthDay(:, 1:2) = zp_d ! resources
      depthDay(:, 3:4) = bent_dn ! resources
!do i=1,nGroups
      depthDay(:, ixStart(1):ixEnd(1)) = spel_dn
      depthDay(:, ixStart(2):ixEnd(2)) = mpel_d
      depthDay(:, ixStart(3):ixEnd(3)) = lpel_d
      depthDay(:, ixStart(4):ixEnd(4)) = bpel_d
      depthDay(:, ixStart(nGroups):ixEnd(nGroups)) = dem_d
!end do
      dayout = 0.d0
      test = 0.d0
      do i = 1, nGrid
         do j = 1, nGrid
            test(:, j) = min(depthDay(:, i), depthDay(:, j))
         end do
         dayout(:, i) = sum(test, 1)
      end do

! calculate overlap during night
      allocate (depthNight(size(xrange), nGrid))
      !test has already allocated
      allocate (nightout(nGrid, nGrid))
      depthNight(:, 1:2) = zp_n ! resources
      depthNight(:, 3:4) = bent_dn ! resources
      depthNight(:, ixStart(1):ixEnd(1)) = spel_dn
      depthNight(:, ixStart(2):ixEnd(2)) = mpel_n
      depthNight(:, ixStart(3):ixEnd(3)) = lpel_n
      depthNight(:, ixStart(4):ixEnd(4)) = bpel_n
      depthNight(:, ixStart(nGroups):ixEnd(nGroups)) = dem_n

      nightout = 0.d0
      test = 0.d0 !reset
      do i = 1, nGrid
         do j = 1, nGrid
            test(:, j) = min(depthNight(:, i), depthNight(:, j))
         end do
         nightout(:, i) = sum(test, 1)
      end do

!!! visual ability
! visual predation is good at light, bad in the dark
      visualpred = [(i, i=ixStart(1), ixEnd(1)), (i, i=ixStart(3), ixEnd(3))] ! small palegic 5 6 always at surface   large pelagic 9 10 11
      dayout(visualpred, :) = dayout(visualpred, :)*visual   ! predation enhance during day
      nightout(visualpred, :) = nightout(visualpred, :)*(2.d0 - visual) ! predation decrease at night

! pelagic predators have limited vision in twilight zone during day
      pelpred = [(i, i=ixStart(3), ixEnd(3))]   ! large pelagic   9 10 11
      pelpred = pelpred(ixadult:size(pelpred)) ! adult large pelagic  11  at dvm during day
      preytwi = [(i, i=ixStart(2), ixEnd(2)), (i, i=ixStart(4), ixEnd(4))] ! mesopelagic 7 8   bathypelagic 12 13 14
      dayout(pelpred, preytwi) = dayout(pelpred, preytwi)/visual*(2.d0 - visual)    ! /visual to restore and then *0.5

!!! average overlap during the whole day
      vertover = (dayout + nightout)*0.5d0
!!! calculate combined feeding preference matrix
      theta = sizeprefer*vertover

!! specific revision of feeding preference
!
      idx_be = [(i, i=idxF, ixStart(5) + (ixjuv - 2))]     ! all pelagic and larval demersals
      theta(idx_be, 3:4) = 0.d0      ! all pelagic and larval demersals do not eat benthos,
      ! only juvenile & adult demersals eat benthos
! small demersals are less preyed on
      idx_smd = [(i, i=ixStart(5) + (ixjuv - 1), ixStart(5) + (ixadult - 2))] !?????????????????????????
      theta(idx_be, idx_smd) = theta(idx_be, idx_smd)*0.25d0
! juvenile & adult demersals do not eat zooplankton
      theta(ixStart(5) + (ixjuv - 1):ixEnd(5), 1:2) = 0.d0
! provide benefit to forage and mesopelagic fish (predator avoidance)
      pred1 = [(i, i=ixStart(3) + (ixadult - 1), ixEnd(3))]
      pred2 = [(i, i=ixStart(4) + (ixadult - 1), ixEnd(4))]
      pred3 = [(i, i=ixStart(5) + (ixadult - 1), ixEnd(5))]
      prey1 = [(i, i=ixStart(1) + (ixjuv - 1), ixEnd(1))]
      prey2 = [(i, i=ixStart(2) + (ixjuv - 1), ixEnd(2))]
      idx_predat = [pred1, pred2, pred3]
      idx_prey = [prey1, prey2]
      theta(idx_predat, idx_prey) = theta(idx_predat, idx_prey)*0.5d0
   end subroutine setupVertical

! nnGroups: Fish group numbers
! nnGrid: all fish grids
   subroutine parametersInit(nnGroups, nnGrid, nnResources, pprod, bprod)
      integer, intent(in):: nnGroups, nnGrid, nnResources
      real(dp), intent(in):: pprod, bprod

      nGroups = nnGroups                   ! fish size spectrum group numbers (species)
      iCurrentGroup = 0                    !
      nResources = nnResources             ! resource numbers
      nGrid = nnGrid + nnResources         ! total grid numbers   resources + total fish stages
      idxF = nResources + 1                ! fish grid begins at...

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
         deallocate (sizeprefer)
         deallocate (vertover)
         deallocate (V)
         deallocate (Enc)
         deallocate (flvl)
         deallocate (Cmax)
         deallocate (mortpred)
         deallocate (mc)
         deallocate (mL)
         deallocate (mU)
         deallocate (K)
         deallocate (rr)
         !maybe more variables
         !
      end if

      allocate (group(nGroups))
      allocate (ixStart(nGroups))
      allocate (ixEnd(nGroups))
      allocate (upositive(nGrid))
      allocate (F(nGrid))
      allocate (theta(nGrid, nGrid))
      allocate (sizeprefer(nGrid, nGrid))
      !allocate (vertover(nGrid, nGrid)) move to theta part
      !

      ! define resources:
      K = [pprod, pprod, bprod, lbenk]    ! Carrying capacity of resources
      rr = [szoog, lzoog, sbeng, lbeng]   ! growth rate of resources
   end subroutine parametersInit

!n:stages of a fish species
   subroutine parametersAddGroup(n, mMax, mMature) ! parametersAddGroup(typeGroup, n, mMax, mMature)
      integer, intent(in) :: n      !number of stages
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

      call initFish(specFish, n, mMax, mMature)
      allocate (group(iCurrentGroup)%spec, source=specFish)

   end subroutine parametersAddGroup
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
            F(i) = F(i) + theta(i, j)*upositive(j) ! matrix multiplication    use upositive
         end do
      end do
!F=matmul(theta,upositive) ! totally same as above

!Calc Encounter/feeding:
!return Enc flvl  Eavail for each group
      do iGroup = 1, nGroups
         call calcFeeding(group(iGroup)%spec, F(ixStart(iGroup):ixEnd(iGroup)))
      end do

! Mortality:
! Predation mortality, including all resources and fish grids    (different from NUM subroutine calcDerivativesUnicellulars)
      ! mortpred = matmul(transpose(theta), (flvl*Cmax/epsAssim*upositive/mc)) ! Petrik et al., 2019

      !if (allocated(vertover))then ! overwritten  van Denderen et al., 2020
      mortpred = (Cmax*V/(Enc + Cmax)*upositive) ! temporarily store
      do i = 1, nGrid
         if (isnan(mortpred(i))) then            ! because some are nan?
            mortpred(i) = 0.d0
         end if
      end do
      mortpred = matmul(transpose(theta), mortpred)
      !end if

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
            call calcfluxfish(spec, upositive(ixStart(iGroup):ixEnd(iGroup)))     ! Flux out and Flux in  Petrik et al., 2019
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
         real(dp), intent(in)::R(nResources), B(this%n)
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

   subroutine formvector(this, iGroup, V, Cmax, mc, mL, mU) ! return assembled vectors containing values for all fish grid (no resources)
      integer, intent(in)::iGroup
      class(spectrumfish)::this
      real(dp), intent(out)::V(nGrid), Cmax(nGrid), mc(nGrid), mL(nGrid), mU(nGrid)

      V(ixStart(iGroup):ixEnd(iGroup)) = this%V
      Cmax(ixStart(iGroup):ixEnd(iGroup)) = this%Cmax
      mc(ixStart(iGroup):ixEnd(iGroup)) = this%m
      mL(ixStart(iGroup):ixEnd(iGroup)) = this%mLower
      mU(ixStart(iGroup):ixEnd(iGroup)) = this%mUpper
   end subroutine formvector

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

      call calcderivatives(u, dudt) ! get flvl mortpred
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

!   subroutine read_namelist_resources()
!      integer :: file_unit, io_err
!
!      namelist /input_resources/ sbenk, lbenk,&
!                                  &szoog, lzoog, sbeng, lbeng
!
!      call open_inputfile(file_unit, io_err)
!      read (file_unit, nml=input_resources, iostat=io_err)
!      call close_inputfile(file_unit, io_err)
!   end subroutine read_namelist_resources

end module FEISTY
