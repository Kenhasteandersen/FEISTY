!
! FEISTY model
! References: Petrik et al., 2019; van Denderen et al., 2020.
! The library follows MATLAB/R codes from Ken Haste Andersen; Pieter Daniël van Denderen; Rémy Denéchère; Daniel Ottmann Riera ...
! Programmed by Yixin Zhao, August 2022.
!
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

   real(dp) :: martin !  martin curve depth
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

!=======================================================================================================
! New from Karline Soetaert package   (Oct 2023 added)
         integer:: nFGrid        ! Total number of fish size grids

         integer :: Rtype   ! resource dynamics, 1=chemostat; 2=logistic

         ! vectors defined in all stages
         real(dp), allocatable:: epsAssim_vec(:)   ! assimilation efficiency
         real(dp), allocatable:: metabolism(:) ! basal respiration rate
         real(dp), allocatable:: mort0(:)      ! basal mortality rate
         real(dp), allocatable:: mortF(:)      ! fishing mortality rate
         ! vectors defined only in fish stages
         real(dp), allocatable:: z(:)          ! ratio of mass in fish size class
         real(dp), allocatable:: psiMature(:)  ! fraction mature in fish size class
         ! vectors defined in fish groups
         real(dp), allocatable:: epsRepro_vec(:)   ! efficiency of reproduction

         real(dp), allocatable:: grazing(:)    ! grazing rate
         real(dp), allocatable:: loss(:)       ! summed loss rates
         real(dp), allocatable:: mort(:)       ! summed mortaliy rates
         real(dp), allocatable:: eAvail(:)     ! available energy

         ! Defined in fish grids only
         real(dp), allocatable:: B(:), dBdt(:) ! state variables and derivatives
         real(dp), allocatable:: eplus(:), eFish(:) ! available energy
         real(dp), allocatable:: grow(:)       ! energy for growth [/yr]
         real(dp), allocatable:: gamma_vec(:)      ! growth to next stage
         real(dp), allocatable:: Repro(:)      ! reproduction rates
         real(dp), allocatable:: mortFish(:)   ! mortality
         real(dp), allocatable:: Fout(:), Fin(:) ! flux in and out of stages

         ! Defined in fish groups only
         real(dp), allocatable:: totMort(:)    ! Total mortality
         real(dp), allocatable:: totGrazing(:) ! Total grazing
         real(dp), allocatable:: totLoss(:)    ! Total losses
         real(dp), allocatable:: totRepro(:)   ! Total reproductive losses
         real(dp), allocatable:: totRecruit(:) ! Total recruitment
         real(dp), allocatable:: totBiomass(:) ! Total biomass

         ! Defined for resources
         real(dp), allocatable:: R(:), dRdt(:) ! state variables and derivatives
         real(dp), allocatable:: mortRes(:)    ! mortality rate pf resources

         ! Other
         real(dp):: tmp(12)
         logical :: feistyinitialised  = .FALSE.

contains
! ======================================
!  Different Setups
! ======================================

! --------------------------------------
! Setup according to Petrik et al., 2019
! --------------------------------------
   subroutine setupbasic(szprod,lzprod, bprod,Ts,Tb)
      real(dp), intent(in)::szprod,lzprod, bprod, Ts, Tb
      integer :: iGroup
! predation preference coefficient
      real(dp) :: thetaS
      real(dp) :: thetaA
      real(dp) :: thetaD
!      real(dp),parameter ::   thetaS = 0.25d0 ! Medium fish pref for small zooplankton
!      real(dp),parameter ::   thetaA = 0.5d0  ! Large fish pref for medium forage fish
!      real(dp),parameter ::   thetaD = 0.75d0 ! Pref of large demersal on pelagic prey

      call read_namelist_setupbasic()
      call parametersInit(3, 2 + 3 + 3, 4, szprod,lzprod, bprod) ! (fish groups, total fish stages 2stages+3stages+3stages, 4 resources, szprod,lzprod, bprod,none)
      call parametersAddGroup(2, 2.50d2, 0.5d0) ! fishSmall
      call parametersAddGroup(3, 1.25d5, 2.5d2) ! fishLarge
      call parametersAddGroup(3, 1.25d5, 2.5d2)   ! fishDemersal

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
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p) ! overwrite metabolism

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

! update temperature
call updateTemp(Ts,Tb)

  !all fish group
  !pelagic
do iGroup = 1, nGroups-1
    group(iGroup)%spec%V=group(iGroup)%spec%V*fTemp
    group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTemp
    group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempm
end do
  !demersal
    group(3)%spec%V=group(3)%spec%V*fTempdem
    group(3)%spec%Cmax=group(3)%spec%Cmax*fTempdem
    group(3)%spec%metabolism=group(3)%spec%metabolism*fTempmdem
  !vector
  !pelagic
V(ixStart(1):ixEnd(2))=V(ixStart(1):ixEnd(2))*fTemp
Cmax(ixStart(1):ixEnd(2))=Cmax(ixStart(1):ixEnd(2))*fTemp
  !demersal
V(ixStart(3):ixEnd(3))=V(ixStart(3):ixEnd(3))*fTempdem
Cmax(ixStart(3):ixEnd(3))=Cmax(ixStart(3):ixEnd(3))*fTempdem

call set2vec
Rtype=1

contains
   subroutine read_namelist_setupbasic()
      integer :: file_unit, io_err

      namelist /input_setupbasic/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
                                  & beta, sigma, mMin, &
                                  & mMedium, mLarge, &
                                  & lbenk, szoog, lzoog, sbeng, lbeng, &
                                  & thetaS, thetaA, thetaD

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_setupbasic, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_setupbasic

   end subroutine setupbasic

! --------------------------------------
! Setup by Ken.
! --------------------------------------
   subroutine setupbasic2(szprod,lzprod, bprod, nStages, Ts, Tb, etaMature) !
      real(dp), intent(in) :: szprod,lzprod, bprod, Ts, Tb
      real(dp), intent(in) :: etaMature ! Mature mass relative to asymptotic size default 0.25, original in van Denderen et al., 2021 was 0.002
      integer, intent(in) :: nStages
      integer :: iGroup, i, j
! predation preference coefficient
      real(dp) :: thetaA
      real(dp) :: thetaD
!      real(dp),parameter ::   thetaA = 0.5d0  ! Large fish pref for medium forage fish
!      real(dp),parameter ::   thetaD = 0.75d0 ! Pref of large demersal on pelagic prey

      call read_namelist_setupbasic2()    !
      call parametersInit(3, nint(0.66d0*nStages) + nStages + nStages, 4, szprod,lzprod, bprod) ! (fish groups, total fish stages, 4 resources,szprod,lzprod, bprod,none)
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, etaMature*2.50d2) ! fishSmall  original mature mass is 0.002 *2.50d2=0.5d0
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishLarge (stages, max mass, mature mass) 0.002 *1.25d5=2.5d2
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishDemersal

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
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p)!overwrite metabolism
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
      !mU = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) ! weight central size
      !mL = c(2e-06,0.001, 0.5e-03, 0.25) ! weight lower limit)

!!basic feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            theta(i, j) = exp(-(log(mc(i)/(beta*mc(j))))**2/(2*sigma)**2)
!            if (mc(j) .gt. mc(i)) then                   !small can't eat large
!               theta(i, j) = 0.d0
!            end if
!         end do
!      end do

! FROM NUM Andersen, K. H., & Visser, A. W. (2023). Appendix https://doi.org/10.1016/j.pocean.2023.102995
! feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            theta(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
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

         if (group(3)%spec%m(i) > mMedium .AND. &
             group(3)%spec%m(i) < mLarge) then
            theta(ixStart(3) + i - 1, 1:2) = 0.d0          ! Medium demersals has not feeding on zooplankton
            theta(ixStart(3) + i - 1, idxF:nGrid) = 0.d0   ! Medium demersals has not feeding on all fish (only eat benthos)
         end if

      end do
      ! Large demersals feed have reduced feeding effiiency on pelagic species:
      theta(ixStart(3):ixEnd(3), ixStart(1):ixEnd(1)) = thetaA*thetaD*theta(ixStart(3):ixEnd(3), ixStart(1):ixEnd(1))
      theta(ixStart(3):ixEnd(3), ixStart(2):ixEnd(2)) = thetaD*theta(ixStart(3):ixEnd(3), ixStart(2):ixEnd(2))

! update temperature
call updateTemp(Ts, Tb)
  !all fish group
  !pelagic
do iGroup = 1, nGroups-1
    group(iGroup)%spec%V=group(iGroup)%spec%V*fTemp
    group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTemp
    group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempm
end do
  !demersal
    group(3)%spec%V=group(3)%spec%V*fTempdem
    group(3)%spec%Cmax=group(3)%spec%Cmax*fTempdem
    group(3)%spec%metabolism=group(3)%spec%metabolism*fTempmdem

  !vector
  !pelagic
V(ixStart(1):ixEnd(2))=V(ixStart(1):ixEnd(2))*fTemp
Cmax(ixStart(1):ixEnd(2))=Cmax(ixStart(1):ixEnd(2))*fTemp
  !demersal
V(ixStart(3):ixEnd(3))=V(ixStart(3):ixEnd(3))*fTempdem
Cmax(ixStart(3):ixEnd(3))=Cmax(ixStart(3):ixEnd(3))*fTempdem

call set2vec
Rtype=1

contains
   subroutine read_namelist_setupbasic2()
      integer :: file_unit, io_err

      namelist /input_setupbasic2/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
                                  & beta, sigma, mMin, &
                                  & mMedium, mLarge, &
                                  & lbenk, szoog, lzoog, sbeng, lbeng, &
                                  & thetaA, thetaD

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_setupbasic2, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_setupbasic2

   end subroutine setupbasic2
! --------------------------------------
! Setup of vertical overlap from MATLAB (van Denderen et al., 2020) the simple run folder
! --------------------------------------
   subroutine setupVertical(szprod,lzprod, bent, nStages, region, bottom, photic, etaMature)
      real(dp), intent(in) :: szprod,lzprod, bottom, photic, bent, etaMature !  default bottom:1500m euphotic depth 150m
      integer, intent(in) :: nStages,region                                  ! Mature mass relative to asymptotic size default 0.25, original in van Denderen et al., 2021 was 0.002

! for theta calc
       real(dp) :: ssigma
       real(dp) :: tau
       !real(dp) :: bottom
       !real(dp) :: photic
       real(dp) :: mesop
       real(dp) :: visual
!      real(dp) :: ssigma = 10.d0
!      real(dp) :: tau = 10.d0
!      real(dp), parameter :: bottom = 1500.d0 ! total depth meter
!      real(dp), parameter :: photic = 150.d0  ! photic zone depth
!      real(dp), parameter :: mesop = 250.d0   ! depth ?
!      real(dp), parameter :: visual = 1.5d0 ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
      real(dp), allocatable :: sigmap(:) ! width for each size class
      !real(dp) :: bent ! for bprod calc
      real(dp) :: bprod
      real(dp), dimension(:), allocatable :: xrange
      real(dp) :: dvm  ! vertical migration depth photic + 500.d0
      real(dp) :: xloc ! vertical location    will be overwritten again and again
      real(dp), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(dp), allocatable :: zp_n(:, :), zp_d(:, :), & ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :)        ! demersal night/ day
      real(dp) :: demmig ! demersal migration
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(dp), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
      real(dp),allocatable :: sizes(:)
      integer :: iGroup, i, j, ixjuv, ixadult, nsize, matstageS, matstageL


      call read_namelist_setupvertical()
      allocate(xrange(int(bottom) + 1))

! calc bprod before initialization
      bprod = 0.1d0*(bent*(bottom/photic)**(-0.86d0)) ! from matlab
      if (bprod .ge. bent*0.1d0) then
         bprod = bent*0.1d0
      end if

      call parametersInit(5, nint(0.66d0*nStages) + nint(0.66d0*nStages) + nStages + nStages + nStages, 4, szprod,lzprod, bprod)!
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, etaMature*2.50d2) ! fishSmall  original mature mass is 0.002 *2.50d2=0.5d0
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, etaMature*2.50d2) ! fishMeso,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishLarge,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishBathy,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishDemersal,

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

         !group(iGroup)%spec%psiMature = 0.d0 ! reset
         !group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5d0! only adults reproduce

      end do

      !overwrite psiMature    from matlab simple run
!      nsize=nStages+1
!      allocate (sizes(nsize))
!      sizes = 10**(linspace(log10(mMin), log10(1.25d5), nsize)) !      mMin=0.001     mMax=1.25d5 predatory fish
!      matstageS = minloc(abs(sizes-0.5d0),dim=1)
!      matstageL = minloc(abs(sizes-2.5d2),dim=1)
!      group(1)%spec%psiMature(matstageS:group(1)%spec%n) = 0.5d0 ! fishSmall
!      group(2)%spec%psiMature(matstageS:group(2)%spec%n) = 0.5d0 ! fishMeso
!      group(3)%spec%psiMature(matstageL:group(3)%spec%n) = 0.5d0 ! fishLarge
!      group(4)%spec%psiMature(matstageL:group(4)%spec%n) = 0.5d0 ! fishBathy
!      group(5)%spec%psiMature(matstageL:group(5)%spec%n) = 0.5d0 ! fishDemersal

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
!! basic feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            sizeprefer(i, j) = sqrt(pi/2.d0)*sigma*( &
!                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)) &
!                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)))
!            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
!         end do
!      end do

! FROM NUM Andersen, K. H., & Visser, A. W. (2023). Appendix https://doi.org/10.1016/j.pocean.2023.102995
! feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            sizeprefer(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
            if (mc(j) .gt. mc(i)) then                   !small can't eat large
               sizeprefer(i, j) = 0.d0
            end if
         end do
      end do

!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0.d0, bottom, int(bottom) + 1)
      dvm = photic + 500.d0 ! 650.d0
      !  from matlab
      if (bottom .lt. (photic + 500.d0)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      end if
      if (bottom .le. mesop) then
         dvm = 0.d0                   ! no migration in shallow habitats
      end if

! first stages as juvenile/adult for predators
      !ixjuv = minloc(abs(sizes-0.5d0),dim=1) !from matlab
      !ixadult = minloc(abs(sizes-2.5d2),dim=1)

      ixjuv = minloc(abs(mL(ixStart(5):ixEnd(5))-0.5d0),dim=1) !predatory fish
      ixadult = minloc(abs(mL(ixStart(5):ixEnd(5))-etaMature*1.25d5),dim=1)


!  deallocate (sizes)!see above overwrite psimature

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
      end if
      if ((bottom - dvm) .ge. 1500.d0) then
         demmig = bottom
      end if
      allocate (dem_d(size(xrange), ixEnd(5) - ixStart(5) + 1))
      xlocvec(ixadult:size(xlocvec)) = dvm ! larvae at surface/ juvenile at bottom/ adult and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_d = matmul(dem_d, diag(1.d0/sum(dem_d, 1)))
! from matlab
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
      depthDay(:, ixStart(1):ixEnd(1)) = spel_dn
      depthDay(:, ixStart(2):ixEnd(2)) = mpel_d
      depthDay(:, ixStart(3):ixEnd(3)) = lpel_d
      depthDay(:, ixStart(4):ixEnd(4)) = bpel_d
      depthDay(:, ixStart(nGroups):ixEnd(nGroups)) = dem_d

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
      idx_smd = [(i, i=ixStart(5) + (ixjuv - 1), ixStart(5) + (ixadult - 2))] ! juvenile demersal is at bottom
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

! update temperature
call updateTempV(depthDay, depthNight, bottom, region)
! all fish group
do iGroup = 1, nGroups
    group(iGroup)%spec%V=group(iGroup)%spec%V*fTempV(ixStart(iGroup):ixEnd(iGroup))
    group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTempV(ixStart(iGroup):ixEnd(iGroup))
    group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempmV(ixStart(iGroup):ixEnd(iGroup))
end do
  !vector
V=V*fTempV
Cmax=Cmax*fTempV

call set2vec
Rtype=1

contains
   subroutine read_namelist_setupvertical()
      integer :: file_unit, io_err

      namelist /input_setupvertical/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
                                  & beta, sigma, mMin, &
                                  & lbenk, szoog, lzoog, sbeng, lbeng,&!bent,
                                  & mMedium, mLarge, &
                                  & ssigma, tau, mesop, visual

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_setupvertical, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_setupvertical

   end subroutine setupVertical
!
!setup for global offline coupling
!input: szprod: small zooplankton production
!       lzprod: large zooplankton production
!       bprod: benthos production
!       bottom: seafloor depth
!       photic: photic zone depth
!       Dgrid: depth grid (boundary), max 200 layers
!       Tprof: temperature profile (boundary), max 200 layers
!       nStages: size number
subroutine setupVerticalGlobal(szprod, lzprod, bprod, bottom, photic, Dgrid, Tprof, nStages, etaMature)
      real(dp), intent(in) :: szprod, lzprod, bprod, bottom, photic, Dgrid(:), Tprof(:), etaMature !
      integer, intent(in) :: nStages
       !real(dp) :: pprod=100.d0 !overwrite later

! for theta calc
       real(dp) :: ssigma
       real(dp) :: tau
       !real(dp) :: bottom
       !real(dp) :: photic
       real(dp) :: mesop
       real(dp) :: visual
!      real(dp) :: ssigma = 10.d0
!      real(dp) :: tau = 10.d0
!      real(dp), parameter :: bottom = 1500.d0 ! total depth meter
!      real(dp), parameter :: photic = 150.d0  ! photic zone depth
!      real(dp), parameter :: mesop = 250.d0   ! depth ?
!      real(dp), parameter :: visual = 1.5d0 ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
      real(dp), allocatable :: sigmap(:) ! width for each size class
      real(dp) :: bent ! for bprod calc
      !real(dp) :: bprod
      real(dp), dimension(:), allocatable :: xrange
      real(dp) :: dvm  ! vertical migration depth photic + 500.d0
      real(dp) :: xloc ! vertical location    will be overwritten again and again
      real(dp), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(dp), allocatable :: zp_n(:, :), zp_d(:, :), & ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :)        ! demersal night/ day
      real(dp) :: demmig ! demersal migration
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(dp), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
      real(dp),allocatable :: sizes(:)
      integer :: iGroup, i, j, ixjuv, ixadult, nsize, matstageS, matstageL


      call read_namelist_setupvertical()
      allocate(xrange(int(bottom) + 1))

!! calc bprod before initialization
!      bprod = 0.1d0*(bent*(bottom/photic)**(-0.86d0)) ! from matlab
!      if (bprod .ge. bent*0.1d0) then
!         bprod = bent*0.1d0
!      end if

      call parametersInit(5, nint(0.66d0*nStages) + nint(0.66d0*nStages) + nStages + nStages + nStages, 4, szprod,lzprod, bprod)!
      !overwrite K
      !K = [szprod, lzprod, bprod, lbenk]
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, etaMature*2.50d2) ! fishSmall  original mature mass is 0.002 *2.50d2=0.5d0
      call parametersAddGroup(nint(0.66d0*nStages), 2.50d2, etaMature*2.50d2) ! fishMeso,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishLarge,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishBathy,
      call parametersAddGroup(nStages, 1.25d5, etaMature*1.25d5) ! fishDemersal,
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
         !group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5d0! only adults reproduce

      end do

      !overwrite psiMature    from matlab simple run
!      nsize=nStages+1
!      allocate (sizes(nsize))
!      sizes = 10**(linspace(log10(mMin), log10(1.25d5), nsize)) !      mMin=0.001     mMax=1.25d5 predatory fish
!      matstageS = minloc(abs(sizes-0.5d0),dim=1)
!      matstageL = minloc(abs(sizes-2.5d2),dim=1)
!      group(1)%spec%psiMature(matstageS:group(1)%spec%n) = 0.5d0 ! fishSmall
!      group(2)%spec%psiMature(matstageS:group(2)%spec%n) = 0.5d0 ! fishMeso
!      group(3)%spec%psiMature(matstageL:group(3)%spec%n) = 0.5d0 ! fishLarge
!      group(4)%spec%psiMature(matstageL:group(4)%spec%n) = 0.5d0 ! fishBathy
!      group(5)%spec%psiMature(matstageL:group(5)%spec%n) = 0.5d0 ! fishDemersal

      !group(nGroups)%spec%mortF(group(nGroups)%spec%n) = 0.5d0 ! only demersal adults have fishing mortality

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
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            sizeprefer(i, j) = sqrt(pi/2.d0)*sigma*( &
!                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)) &
!                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)))
!            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
!         end do
!      end do

      do i = idxF, nGrid
         do j = 1, nGrid
            sizeprefer(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
            if (mc(j) .gt. mc(i)) then                   !small can't eat large
               sizeprefer(i, j) = 0.d0
            end if
         end do
      end do

!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0.d0, bottom, int(bottom) + 1)
      dvm = photic + 500.d0 ! 650.d0
      !  from matlab
      if (bottom .lt. (photic + 500.d0)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      end if
      if (bottom .le. mesop) then
         dvm = 0.d0                   ! no migration in shallow habitats
      end if

! first stages as juvenile/adult for predators
!      ixjuv = minloc(abs(sizes-0.5d0),dim=1) !from matlab
!      ixadult = minloc(abs(sizes-2.5d2),dim=1)


      ixjuv = minloc(abs(mL(ixStart(5):ixEnd(5))-0.5d0),dim=1) !predatory fish
      ixadult = minloc(abs(mL(ixStart(5):ixEnd(5))-etaMature*1.25d5),dim=1)

!  deallocate (sizes)!see above overwrite psimature

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
      end if
      if ((bottom - dvm) .ge. 1500.d0) then
         demmig = bottom
      end if
      allocate (dem_d(size(xrange), ixEnd(5) - ixStart(5) + 1))
      xlocvec(ixadult:size(xlocvec)) = dvm ! larvae at surface/ juvenile at bottom/ adult and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_d = matmul(dem_d, diag(1.d0/sum(dem_d, 1)))
! from matlab
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
      depthDay(:, ixStart(1):ixEnd(1)) = spel_dn
      depthDay(:, ixStart(2):ixEnd(2)) = mpel_d
      depthDay(:, ixStart(3):ixEnd(3)) = lpel_d
      depthDay(:, ixStart(4):ixEnd(4)) = bpel_d
      depthDay(:, ixStart(nGroups):ixEnd(nGroups)) = dem_d

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
      idx_smd = [(i, i=ixStart(5) + (ixjuv - 1), ixStart(5) + (ixadult - 2))] ! juvenile demersal is at bottom
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

! update temperature
call updateTempVG(Dgrid, Tprof, depthDay, depthNight, bottom)
! all fish group
do iGroup = 1, nGroups
    group(iGroup)%spec%V=group(iGroup)%spec%V*fTempV(ixStart(iGroup):ixEnd(iGroup))
    group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTempV(ixStart(iGroup):ixEnd(iGroup))
    group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempmV(ixStart(iGroup):ixEnd(iGroup))
end do
  !vector
V=V*fTempV
Cmax=Cmax*fTempV

contains
   subroutine read_namelist_setupvertical()
      integer :: file_unit, io_err

      namelist /input_setupvertical/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
                                  & beta, sigma, mMin, &
                                  & mMedium, mLarge, &
                                  & bent, lbenk, szoog, lzoog, sbeng, lbeng,&
                                  & ssigma, tau, mesop, visual

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_setupvertical, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_setupvertical

   end subroutine setupVerticalGlobal


! --------------------------------------
! Setup of vertical overlap (van Denderen et al., 2020) with squid Remy & Daniel Ottmann
! --------------------------------------
   subroutine setupsquid(szprod,lzprod, bottom, nStages)
      real(dp), intent(in) :: szprod,lzprod !
      real(dp), intent(in) :: bottom ! water depth default 1000.d0     revise input.nml
      integer, intent(in) :: nStages ! stage numbers

! for theta calc
       real(dp) :: ssigma
       real(dp) :: tau
       real(dp) :: photic
       real(dp) :: mesop
       real(dp) :: visual
       real(dp) :: S2P
!      real(dp) :: ssigma = 10.d0
!      real(dp) :: tau = 10.d0
!      real(dp), parameter :: photic = 150.d0  ! photic zone depth
!      real(dp), parameter :: mesop = 250.d0   ! depth ?
!      real(dp), parameter :: visual = 1.5d0   ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
!      real(dp), parameter :: S2P = 0.5d0        ! predation from Squid to pelagics
      real(dp), allocatable :: mortFi(:) ! default fishing intensity
      real(dp), allocatable :: maxfishi(:) ! vector of max mass of fish
      real(dp), allocatable :: sigmap(:) ! width for each size class
      real(dp) :: bprod
      real(dp) :: mat_const, smaxfish, lmaxfish, smat, lmat ! for maturity mass calc
      real(dp), dimension(:), allocatable ::  sizes, xrange
      real(dp) :: dvm  ! vertical migration depth photic + 500.d0
      real(dp) :: xloc ! vertical location    will be overwritten again and again
      real(dp), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(dp), allocatable :: zp_n(:, :), zp_d(:, :), &       ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :) , &    ! demersal night/ day
                               cph_n(:,:), cph_d(:,:)          ! suqid night/day
      real(dp) :: demmig ! ?
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(dp), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
      integer :: iGroup, i, j, ixjuv, ixadult

      call read_namelist_setupsquid()
      allocate(xrange(int(bottom) + 1))
      allocate(sizes(nStages+1))

  martin = min((bottom / photic)**(-0.86d0), 1.d0)

  mat_const = 0.28d0 ! Andersen 2019 pp45
  smaxfish = 2.50d2 ! max size of small pelagics   boundary
  lmaxfish = 1.25d5 ! max size of predator fish   boundary
  smat = smaxfish * mat_const !  weight at maturity forage/meso
  lmat = lmaxfish * mat_const ! weight at maturity predators
  sizes = exp(linspace(log(mMin), log(lmaxfish), nStages + 1))
  !sizes = 10.d0**(linspace(log10(mMin), log10(lmaxfish), nStages + 1))

!! calc bprod before initialization
!      bprod = 0.1d0*(bent*(bottom/photic)**(-0.86d0)) ! from matlab
!      if (bprod .ge. bent*0.1d0) then
!         bprod = bent*0.1d0
!      end if
      bprod=0.d0 ! no benthos production

      call parametersInit(6, 2*nint(0.66d0*nStages) + 4*nStages, 4, szprod,lzprod, bprod)!
      call parametersAddGroup(nint(0.66d0*nStages), smaxfish, smat) ! fishSmall
      call parametersAddGroup(nint(0.66d0*nStages), smaxfish, smat) ! fishMeso  (stages, max mass, mature mass)
      call parametersAddGroup(nStages, lmaxfish, lmat)              ! fishLarge
      call parametersAddGroup(nStages, lmaxfish, lmat)              ! fishDemersal
      mMin = 0.01d0  ! for squid   smallest size all cephalopod
      call parametersAddGroup(nStages, 3.5d3, 1.d10)                ! Squid  maturity mass is not used, psiMature will be overwritten later
      mMin = 0.001d0 ! restore for fish
      call parametersAddGroup(nStages, lmaxfish, lmat)              ! fishBathy

! vectors and matrix:
      allocate (mortFi(nGroups))
      allocate (maxfishi(nGroups))
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
      mortFi = [0.3d0, 0.d0, 0.3d0, 0.3d0, 0.3d0, 0.d0] ! default fishing intensity.
      maxfishi = [2.50d2, 2.50d2, 1.25d5, 1.25d5, 3.5d3, 1.25d5] ! vector of max mass of fish

! Overwrite--------------
      ! all
      do iGroup = 1, nGroups
         ! metabolism
         group(iGroup)%spec%metabolism = (0.2d0*h*group(iGroup)%spec%m**p)
         ! mortF
         group(iGroup)%spec%mortF = mortFi(iGroup) * &
                                    (1.d0 + (group(iGroup)%spec%m / &
                                    (group(iGroup)%spec%mUpper(group(iGroup)%spec%n) * 0.05d0))**(-3.d0))**(-1.d0)
         ! psiMature
         group(iGroup)%spec%psiMature = 0.d0 ! reset
         if (iGroup .le. 2) then ! small
         group(iGroup)%spec%psiMature = (1.d0 + (group(iGroup)%spec%m / smat)**(-5.d0))**(-1.d0) * &
                         (group(iGroup)%spec%m / maxfishi(iGroup))**(1.d0-(nn+1.d0))! from Daniel Ottmann Riera?
         else ! large
         group(iGroup)%spec%psiMature = (1.d0 + (group(iGroup)%spec%m / lmat)**(-5.d0))**(-1.d0) * &
                         (group(iGroup)%spec%m / maxfishi(iGroup))**(1.d0-(nn+1.d0))! from Daniel Ottmann Riera?
         end if

      end do
      ! squid
      group(5)%spec%psiMature = 0.d0
      group(5)%spec%Cmax = hCepha*(group(5)%spec%m**nn)         ! Maximum consumption rate of squid
      group(5)%spec%metabolism = (0.2d0*hCepha*group(5)%spec%m**p)

!  -------------------------
! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formvector(spec, iGroup, V, Cmax, mc, mL, mU)
         end select
      end do

      mc(1:nResources) = [2.d-06*sqrt(500.d0), 1.d-3*sqrt(500.d0), 0.5d-03*sqrt(250000.d0), 0.25d0*sqrt(500.d0)] ! resource central mass
      mU(1:nResources) = [0.001d0, 0.5d0, 125.d0, 125.d0]  ! resource mass upper limit
      mL(1:nResources) = [2.d-6, 0.001d0, 0.5d-3, 0.25d0] ! resource mass lower limit
! basic feeding preference matrix theta
      do i = idxF, nGrid ! all
         do j = 1, nGrid
            sizeprefer(i, j) = sqrt(pi/2.d0)*sigma*( &
                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)) &
                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2.d0)*sigma)))
            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
         end do
      end do
  ! overwrite squid preference
      do i = ixStart(5), ixEnd(5) ! squid
         do j = 1, nGrid
            sizeprefer(i, j) = sqrt(pi/2.d0)*sigma*( &
                               erf((log(mU(j)) - log(mc(i)/betaCepha))/(sqrt(2.d0)*sigma)) &
                               - erf((log(mL(j)) - log(mc(i)/betaCepha))/(sqrt(2.d0)*sigma)))
            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
         end do
      end do

!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0.d0, bottom, int(bottom) + 1)
      dvm = photic + 500.d0 ! 650.d0

      if (bottom .lt. (photic + 500.d0)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      end if
      if (bottom .le. mesop) then
         dvm = 0.d0                   ! no migration in shallow habitats
      end if

! first stages as juvenile/adult for predators
      ixjuv = minloc(abs(sizes-smat/15.d0),dim=1)
      ixadult = minloc(abs(sizes-lmat),dim=1)
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

! demersal fish night
      allocate (dem_n(size(xrange), ixEnd(4) - ixStart(4) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(4) - ixStart(4) + 1))
      allocate (ix(ixEnd(4) - ixStart(4) + 1))
      xlocvec = 0.d0 ! initialization
      xlocvec(ixjuv:size(xlocvec)) = bottom  !  larvae at surface   juvenile and adult at bottom
      ix = [(i, i=ixStart(4), ixEnd(4))]
      do i = 1, size(ix)
         dem_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_n = matmul(dem_n, diag(1.d0/sum(dem_n, 1)))

! demersal fish day
      demmig = dvm ! ??? from matlab
      if ((bottom - dvm) .ge. 1200.d0) then
         demmig = dvm + (bottom - dvm - 1200.d0)
         end if
      if ((bottom - dvm) .ge. 1500.d0) then
         demmig = bottom
      end if
      allocate (dem_d(size(xrange), ixEnd(4) - ixStart(4) + 1))
      xlocvec(ixadult:size(xlocvec)) = dvm ! larvae at surface/ juvenile at bottom/ adult and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      dem_d = matmul(dem_d, diag(1.d0/sum(dem_d, 1)))
!  from matlab
      ! if shallower than euphotic depth, adult demersals feed across-habitats
      if (bottom .le. photic) then
         dem_d = (dem_d + dem_n)/2.d0
         dem_n = dem_d
      end if

! squid night
      allocate (cph_n(size(xrange), ixEnd(5) - ixStart(5) + 1)) ! cephalopod night
      deallocate (xlocvec)
      deallocate (ix)
      cph_n=lpel_n

! squid day !................ need check
      allocate (cph_d(size(xrange), ixEnd(5) - ixStart(5) + 1)) ! cephalopod day
      allocate (xlocvec(ixEnd(5) - ixStart(5) + 1))
      allocate (ix(ixEnd(5) - ixStart(5) + 1))
      xlocvec = 0.d0 ! initialization
      if (bottom .gt. mesop) then
      xlocvec = dvm
      else                                  !!................ need check
      xlocvec = bottom
      end if
      ix = [(i, i=ixStart(5), ixEnd(5))]
      do i = 1, size(ix)
         cph_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                       exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      cph_d = matmul(cph_d, diag(1.d0/sum(cph_d, 1)))

! bathypelagic night (adults in midwater, others at surface)
      allocate (bpel_n(size(xrange), ixEnd(6) - ixStart(6) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(6) - ixStart(6) + 1))
      allocate (ix(ixEnd(6) - ixStart(6) + 1))
      xlocvec = 0.d0 ! initialization
      xlocvec(ixadult:size(xlocvec)) = dvm  !  non-adult at surface   adult at dvm
      ix = [(i, i=ixStart(6), ixEnd(6))]
      do i = 1, size(ix)
         bpel_n(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      bpel_n = matmul(bpel_n, diag(1.d0/sum(bpel_n, 1)))

! bathypelagic day (all at dvm)
      allocate (bpel_d(size(xrange), ixEnd(6) - ixStart(6) + 1))
      xlocvec = dvm ! overwrite all elements by dvm
      do i = 1, size(ix)
         bpel_d(:, i) = (1.d0/(sqrt(2.d0*pi*sigmap(ix(i))**2.d0)))* &
                        exp(-((xrange - xlocvec(i))**2.d0/(2.d0*sigmap(ix(i))**2.d0)))
      end do
      bpel_d = matmul(bpel_d, diag(1.d0/sum(bpel_d, 1)))

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
      depthDay(:, ixStart(4):ixEnd(4)) = dem_d
      depthDay(:, ixStart(5):ixEnd(5)) = cph_d
      depthDay(:, ixStart(nGroups):ixEnd(nGroups)) = bpel_d
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
      depthNight(:, ixStart(4):ixEnd(4)) = dem_n
      depthNight(:, ixStart(5):ixEnd(5)) = cph_n
      depthNight(:, ixStart(nGroups):ixEnd(nGroups)) = bpel_n

      nightout = 0.d0
      test = 0.d0 !reset
      do i = 1, nGrid
         do j = 1, nGrid
            test(:, j) = min(depthNight(:, i), depthNight(:, j))
         end do
         nightout(:, i) = sum(test, 1)
      end do

!!! visual ability                        !!................ need check...............
! visual predation is good at light, bad in the dark
      visualpred = [(i, i=ixStart(1), ixEnd(1)), (i, i=ixStart(3), ixEnd(3))] ! small palegic    large pelagic
      !dayout(visualpred, :) = dayout(visualpred, :)*visual   ! predation enhance during day
      !nightout(visualpred, :) = nightout(visualpred, :)*(2.d0 - visual) ! predation decrease at night

  deallocate(ix)
  allocate(ix(ixStart(5)))
  ix = [(i,i=1,ixStart(5))] ! ..................only 1st stage squid?
    dayout(visualpred, ix) = dayout(visualpred, ix) * visual   ! predation enhance during day
  deallocate(ix)
  allocate(ix(ixEnd(5)-ixStart(1)+1))
  ix = [(i,i=ixStart(1),ixEnd(5))]
    nightout(visualpred, ix) = nightout(visualpred, ix)*(2.d0 - visual) ! all prey are less prayed at night

! pelagic predators have limited vision in twilight zone during day
      pelpred = [(i, i=ixStart(3), ixEnd(3))]   ! large pelagic
      pelpred = pelpred(ixadult:size(pelpred)) ! adult large pelagic    at dvm during day
      preytwi = [(i, i=ixStart(2), ixEnd(2)), (i, i=ixStart(6), ixEnd(6))] ! mesopelagic    bathypelagic
      dayout(pelpred, preytwi) = dayout(pelpred, preytwi)/visual*(2.d0 - visual)    ! /visual to restore and then *0.5

! Squid    from R
  deallocate(ix)
  allocate(ix(ixEnd(5)-ixStart(1)+1))
  ix = [(i,i=ixStart(5),ixEnd(5))]

  if (bottom .gt. mesop) then ! In deep regions, cephalopods are hiding in the mesopelgaic regions during the day (less predation)
                              ! and are going up at night (less predation from visual predators) -> we multiply x 0.75
    dayout(visualpred, ix) = dayout(visualpred, ix) * 0.75d0
    nightout(visualpred, ix) = nightout(visualpred, ix) * 0.75d0
  end if

!!! average overlap during the whole day
      vertover = (dayout + nightout)*0.5d0
!!! calculate combined feeding preference matrix
      theta = sizeprefer*vertover

!! specific revision of feeding preference
!
      idx_be = [(i, i=idxF, ixStart(4) + (ixjuv - 2)), (i, i=ixStart(5), ixEnd(6)) ]  ! all pelagic and larval demersals & squid
      theta(idx_be, 3:4) = 0.d0      ! all pelagic and larval demersals do not eat benthos,
      ! only juvenile & adult demersals eat benthos
! small demersals are less prayed on
      idx_smd = [(i, i=ixStart(4) + (ixjuv - 1), ixStart(4) + (ixadult - 2))] !??
      theta(idx_be, idx_smd) = theta(idx_be, idx_smd)*0.25d0
! juvenile & adult demersals do not eat zooplankton
      theta(ixStart(4) + (ixjuv - 1):ixEnd(4), 1:2) = 0.d0

! provide benefit to forage and mesopelagic fish (predator avoidance)
      pred1 = [(i, i=ixStart(3) + (ixadult - 1), ixEnd(3))] ! large pelagic
      pred2 = [(i, i=ixStart(4) + (ixadult - 1), ixEnd(4))] ! demersal
      pred3 = [(i, i=ixStart(6) + (ixadult - 1), ixEnd(6))] ! bathypelagics
      prey1 = [(i, i=ixStart(1) + (ixjuv - 1), ixEnd(1))]   ! small pelagics
      prey2 = [(i, i=ixStart(2) + (ixjuv - 1), ixEnd(2))]   ! mesopelagic
      idx_predat = [pred1, pred2, pred3]
      idx_prey = [prey1, prey2]
      theta(idx_predat, idx_prey) = theta(idx_predat, idx_prey)*0.5d0

! adjustment of Squid
  deallocate(ix)
  allocate(ix(ixEnd(5)-ixStart(5)+1))
  ix = [(i,i=ixStart(5),ixEnd(5))]   !idx of squid
  theta(ix,ixStart(3):ixEnd(3))= theta(ix,ixStart(3):ixEnd(3)) * S2P  ! S2P=0.5d0
  theta(ix,ixStart(4):(ixStart(4)+ixjuv - 2))= theta(ix,ixStart(4):(ixStart(4)+ixjuv - 2)) * S2P ! larval demersals?
  theta(ix,prey1)= theta(ix,prey1) * S2P

contains
   subroutine read_namelist_setupsquid()
      integer :: file_unit, io_err

      namelist /input_setupsquid/ h, hCepha, nn, q, gamma, kk, p, epsAssim, epsRepro, epst, &
                                  & beta, betaCepha, sigma, mMin, &
                                  & mMedium, mLarge, &
                                  & lbenk, szoog, lzoog, sbeng, lbeng,&
                                  & ssigma, tau, photic, mesop, visual, S2P

      call open_inputfile(file_unit, io_err)
      read (file_unit, nml=input_setupsquid, iostat=io_err)
      call close_inputfile(file_unit, io_err)
   end subroutine read_namelist_setupsquid

   end subroutine setupsquid


! =============================
! Initialization
! =============================

!--------------------------------
! Initialize
! input:
! nnGroups: Fish group numbers
! nnGrid: all fish grids (no resources)
! nnResources: resource numbers
! szprod: small zooplankton carrying capacity
! lzprod: large zooplankton carrying capacity
! bprod: small benthos carrying capacity
! -------------------------------
   subroutine parametersInit(nnGroups, nnGrid, nnResources, szprod,lzprod, bprod)
      integer, intent(in):: nnGroups, nnGrid, nnResources
      real(dp), intent(in):: szprod,lzprod, bprod

      nGroups = nnGroups                   ! fish size spectrum group numbers (species)
      iCurrentGroup = 0
      nResources = nnResources             ! resource numbers
      nGrid = nnGrid + nnResources         ! total grid numbers   resources + total fish stages
      idxF = nResources + 1                ! fish grid begins at...
      nFGrid=nGrid-nResources
!
! Allocate/deallocate variables:
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

!-----------------Oct 2023 add-----------
         deallocate(epsAssim_vec)
         deallocate(metabolism)
         deallocate(mort0)
         deallocate(mortF)
         deallocate(z)
         deallocate(psiMature)
         deallocate(epsRepro_vec)
         deallocate(grazing)
         deallocate(loss)
         deallocate(mort)
         deallocate(EAvail)
         deallocate(B)
         deallocate(dBdt)
         deallocate(eplus)
         deallocate(eFish)
         deallocate(grow)
         deallocate(gamma_vec)
         deallocate(Repro)
         deallocate(mortFish)
         deallocate(Fout)
         deallocate(Fin)
         deallocate(totMort)
         deallocate(totGrazing)
         deallocate(totLoss)
         deallocate(totRepro)
         deallocate(totRecruit)
         deallocate(totBiomass)
         deallocate(R)
         deallocate(dRdt)
         deallocate(mortRes)
!-----------------------------------------

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

!-----------------Oct 2023 add-----------
      allocate(epsRepro_vec(nGroups))
      allocate(psiMature(nFGrid))
      allocate(z(nFGrid))
      allocate(epsAssim_vec(nGrid))
      allocate(metabolism(nGrid))
      allocate(mort(nGrid))
      allocate(mort0(nGrid))
      allocate(mortF(nGrid))
      allocate(R(nResources))
      allocate(dRdt(nResources))
      allocate(mortRes(nResources))
      allocate(grow(nFGrid))
      allocate(eplus(nFGrid))
      allocate(mortFish(nFGrid))
      allocate(eFish(nFGrid))
      allocate(B(nFGrid))
      allocate(dBdt(nFGrid))
      allocate(gamma_vec(nFGrid))
      allocate(Fout(nFGrid))
      allocate(Fin(nFGrid))
      allocate(Repro(nFGrid))
      allocate(Eavail(nGrid))
      allocate(grazing(nGrid))
      allocate(loss(nGrid))
      allocate(totMort(nGroups))
      allocate(totGrazing(nGroups))
      allocate(totLoss(nGroups))
      allocate(totRepro(nGroups))
      allocate(totRecruit(nGroups))
      allocate(totBiomass(nGroups))
!---------------------------------------

      ! define resources:
      K = [szprod, lzprod, bprod, lbenk]    ! Carrying capacity of resources [g m-2]]
      rr = [szoog, lzoog, sbeng, lbeng]   ! growth rate of resources       [yr-1]
   end subroutine parametersInit

! --------------------------
! Add a size spectrum group of fish
! input:
! n: stages of a fish species
! mMax: max fish size the of the species (boundary of the grid)
! mMature: mature size (middle point of the grid)
! --------------------------
   subroutine parametersAddGroup(n, mMax, mMature)
      integer, intent(in) :: n      !number of stages
      real(dp), intent(in) :: mMax, mMature

      type(spectrumFish) :: specFish
!
! define idx for different fish size group:
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

!
! =====================================
! derivative calculation
! =====================================
!
! ----------------------------------------------------------------------
!  Calculate the derivatives for all groups:
!  In:
!  u: the vector of state variables (all resources and all fish grids)
!  dudt: vector to hold the derivative (input and output)
! ----------------------------------------------------------------------

!Replaced by vectorized 'calcderivatives' in Karline Soetaert package below.




!-------------------------------------------------------------
! return assembled vectors containing values for all fish grid (no resources)
   subroutine formvector(this, iGroup, V, Cmax, mc, mL, mU)
      integer, intent(in)::iGroup
      class(spectrumfish)::this
      real(dp), intent(out)::V(nGrid), Cmax(nGrid), mc(nGrid), mL(nGrid), mU(nGrid)

      V(ixStart(iGroup):ixEnd(iGroup)) = this%V
      Cmax(ixStart(iGroup):ixEnd(iGroup)) = this%Cmax
      mc(ixStart(iGroup):ixEnd(iGroup)) = this%m
      mL(ixStart(iGroup):ixEnd(iGroup)) = this%mLower
      mU(ixStart(iGroup):ixEnd(iGroup)) = this%mUpper
   end subroutine formvector

   ! FROM NUM
    ! Calculate the interaction coefficient between two size groups.
    ! In:
    !   z : The predator:prey body mass ratio between the two groups
    !   beta: preferred predator:prey body mass ratio
    !   sigma: width of selection
    !   Delta: ratio between upper and lower body mass in size groups
    !
    function calcPhi(z, beta,sigma, Delta) result(res)
      real(dp), intent(in):: z,beta,sigma,Delta
      real(dp):: res, s

      if (beta .eq. 0.d0) then
         res = 0.d0 ! beta = 0 is interpreted as if the group is not feeding
      else
         s = 2*sigma*sigma
         res = max(0.d0, &
         (Sqrt(Delta)*(((exp(-Log((beta*Delta)/z)**2/s) - 2/exp(Log(z/beta)**2/s) + &
         exp(-Log((Delta*z)/beta)**2/s))*s)/2. - &
         (Sqrt(Pi)*Sqrt(s)*(Erf((-Log(beta*Delta) + Log(z))/Sqrt(s))*Log((beta*Delta)/z) + &
         2*Erf(Log(z/beta)/Sqrt(s))*Log(z/beta) + &
         Erf((Log(beta) - Log(Delta*z))/Sqrt(s))*Log((Delta*z)/beta)))/2.))/ &
         ((-1 + Delta)*Log(Delta)) )
      end if
    end function calcPhi

! =========================
! rates
! =========================

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
      integer :: iGroup

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


! Oct 2023
! assign parameter set from Yixin Zhao library to vectors in Karline Soetaert package
subroutine set2vec
 integer :: iGroup

epsAssim_vec=epsAssim
epsRepro_vec=epsRepro


metabolism=0.d0
mort0=0.d0
mortF=0.d0
psiMature=0.d0
z=0.d0

do iGroup=1,nGroups

metabolism( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%metabolism
mort0( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%mort0
mortF( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%mortF
psiMature( ixStart(iGroup)-nResources : ixEnd(iGroup)-nResources )    =   group(iGroup)%spec%psiMature
z( ixStart(iGroup)-nResources : ixEnd(iGroup)-nResources )   = group(iGroup)%spec%z

end do
end subroutine set2vec

!============================================================
! New stuff based on Karline Soetaert package (Oct 2023 added)
!============================================================

! =====================================
! Initialize:
! =====================================

!--------------------------------------
! read values
! from long vector passed from R
! update index to vector (ir)
!--------------------------------------

   subroutine getVec(vec, rpar, ir, n)
      implicit none
      integer,  intent(in)   :: n
      integer,  intent(inout):: ir
      real(dp), intent(out)  :: vec(n)
      real(dp), intent(in)   :: rpar(*)
      integer :: i

      do i = 1, n
        vec(i) = rpar(ir)
        ir=ir+1
      end do
   end subroutine getVec

!--------------------------------------
! allocate/deallocate variables
! pass values from R
!
! input:
! ipar: integer inputs from R
! rpar: double precision inputs from R
!
! output:
! allocated vectors/matrices
! indices to groups
! input values
!--------------------------------------

   subroutine allocfeisty(ipar, rpar)

      implicit none
      integer, intent(in):: ipar(*)
      real(dp), intent(in):: rpar(*)

      integer:: i, j, ii, ir


      ii = 4    !  first usable element in ipar

     !--------------------------
     ! number of groups, stages
     !--------------------------
      nGroups = ipar(ii)
      ii = ii+1

      nResources = ipar(ii)
      ii = ii+1

      idxF       = nResources + 1             ! fish grid begins at...

     !--------------------------
     ! indices to stages, indices
     !--------------------------
      if (allocated (ixStart))  deallocate (ixStart)
      allocate (ixStart(nGroups))

      if (allocated (ixEnd))    deallocate (ixEnd)
      allocate (ixEnd(nGroups))

      ! First and last indices within fish grid
      ixStart(1) = 1   +nResources
      do i = 2, nGroups
        ixstart(i) = ixstart(i-1) + ipar(ii)
        ixEnd(i-1) = ixstart(i)-1
        ii = ii+1
      end do
      ixEnd(nGroups) = ixstart(nGroups)+ipar(ii)-1
      ii = ii+1
      Rtype = ipar(ii)
      ii = ii+1

      nFGrid = ixEnd(nGroups) -nResources
      nGrid  = nFgrid + nResources
!CALL intpr ("nFgrid ", 7, nFgrid, 1)  ! to check if correctly passed

     !--------------------------
     ! resource parameters
     !--------------------------
      ir = ipar(1)+1   !  first usable element in rpar

      if (allocated (K))   deallocate (K)
      allocate (K(nResources))

      if (allocated (rr))  deallocate (rr)
      allocate (rr(nResources))

      call getVec(K,  rpar, ir, nResources)
      call getVec(rr, rpar, ir, nResources)

!     CALL intpr ("ir   ", 4, ir, 1)
!     tmp(:) = 0.d0
!      do i=1,4
!        tmp(i) = K(i)
!      end do
!      CALL dblepr ("k R ", 4, tmp, 12)
!      call Rexit("stop")

     !--------------------------
     ! fish group parameters
     !--------------------------

      if (allocated (epsRepro_vec)) deallocate (epsRepro_vec)
      allocate (epsRepro_vec(nGroups))

      call getVec(epsRepro_vec,  rpar, ir, nGroups)

     !--------------------------
     ! fish stage parameters
     !--------------------------

      if (allocated (psiMature))  deallocate (psiMature)
      allocate (psiMature(nFGrid))

      if (allocated (z))          deallocate (z)
      allocate (z(nFGrid))

      call getVec(psiMature,  rpar, ir, nFGrid)
      call getVec(z,          rpar, ir, nFGrid)

     !--------------------------
     ! parameters for resource+fish:
     !--------------------------

      if (allocated (theta))  deallocate (theta)
      allocate (theta(nGrid, nGrid))

      do i = 1, nGrid
        do j = 1, nGrid
          theta(i,j) = rpar(ir)
          ir = ir+1
        end do
      end do

      if (allocated (epsAssim_vec))    deallocate (epsAssim_vec)
      allocate (epsAssim_vec(nGrid))

      if (allocated (V))           deallocate (V)
      allocate (V(nGrid))

      if (allocated (Cmax))        deallocate (Cmax)
      allocate (Cmax(nGrid))

      if (allocated (metabolism))  deallocate (metabolism)
      allocate (metabolism(nGrid))

      if (allocated (mort))        deallocate (mort)
      allocate (mort(nGrid))

      if (allocated (mort0))       deallocate (mort0)
      allocate (mort0(nGrid))

      if (allocated (mortF))       deallocate (mortF)
      allocate (mortF(nGrid))

      call getVec(epsAssim_vec,    rpar, ir, nGrid)
      call getVec(V,           rpar, ir, nGrid)
      call getVec(Cmax,        rpar, ir, nGrid)
      call getVec(metabolism,  rpar, ir, nGrid)
      call getVec(mort0,       rpar, ir, nGrid)
      call getVec(mortF,       rpar, ir, nGrid)


     ! resource vectors to be calculated in R-code

      if (allocated (R))  deallocate (R)
      allocate (R(nResources))

      if (allocated (dRdt))  deallocate (dRdt)
      allocate (dRdt(nResources))

      if (allocated (mortRes))  deallocate (mortRes)
      allocate (mortRes(nResources))

     ! fish grid vectors to be calculated in R-code

      if (allocated (grow))  deallocate (grow)
      allocate (grow(nFGrid))

      if (allocated (eplus))  deallocate (eplus)
      allocate (eplus(nFGrid))

      if (allocated (mortFish))  deallocate (mortFish)
      allocate (mortFish(nFGrid))

      if (allocated (eFish))  deallocate (eFish)
      allocate (eFish(nFGrid))

      if (allocated (B))  deallocate (B)
      allocate (B(nFGrid))

      if (allocated (dBdt))  deallocate (dBdt)
      allocate (dBdt(nFGrid))

      if (allocated (gamma_vec))  deallocate (gamma_vec)
      allocate (gamma_vec(nFGrid))

      if (allocated (Fout))  deallocate (Fout)
      allocate (Fout(nFGrid))

      if (allocated (Fin))  deallocate (Fin)
      allocate (Fin(nFGrid))

      if (allocated (Repro))  deallocate (Repro)
      allocate (Repro(nFGrid))

     ! resource+fish vectors to be calculated in R-code

      if (allocated (flvl))  deallocate (flvl)
      allocate (flvl(nGrid))

      if (allocated (Enc))  deallocate (Enc)
      allocate (Enc(nGrid))

      if (allocated (Eavail))  deallocate (Eavail)
      allocate (Eavail(nGrid))

      if (allocated (mortpred))  deallocate (mortpred)
      allocate (mortpred(nGrid))

      if (allocated (grazing))  deallocate (grazing)
      allocate (grazing(nGrid))

      if (allocated (loss))  deallocate (loss)
      allocate (loss(nGrid))

     ! Fish group vectors to be calculated in R-code

      if (allocated (totMort)) deallocate (totMort)
      allocate (totMort(nGroups))

      if (allocated (totGrazing)) deallocate (totGrazing)
      allocate (totGrazing(nGroups))

      if (allocated (totLoss)) deallocate (totLoss)
      allocate (totLoss(nGroups))

      if (allocated (totRepro)) deallocate (totRepro)
      allocate (totRepro(nGroups))

      if (allocated (totRecruit)) deallocate (totRecruit)
      allocate (totRecruit(nGroups))

      if (allocated (totBiomass)) deallocate (totBiomass)
      allocate (totBiomass(nGroups))

   end subroutine allocfeisty

!--------------------------------------

   subroutine checknan(vec, n)
      integer, intent(in)::n
      real(dp), intent(inout) :: vec(n)
      integer :: i

      do i = 1, n
         if (isnan(vec(i))) vec(i) = 0.d0
      end do
   end subroutine checknan


! =====================================
! derivative calculation
! =====================================

! ----------------------------------------------------------------------
!  Calculate the derivatives for all groups:
!  In:
!  u: vector of state variables (all resources and fish grids, input)
!  dudt: vector to hold the derivative (input and output)
! ----------------------------------------------------------------------

  subroutine calcderivatives(u, dudt)

      real(dp), intent(in)    :: u(nGrid)
      real(dp), intent(inout) :: dudt(nGrid)

      integer :: i, j, ii, istart, istop, iGroup

! ----------------------------------------------------------------------

! ----------------------------------------------
! Feeding $ losses for resources and fish grids:
! ----------------------------------------------
      Enc  = V*matmul(theta, u)                         ! Encounter rates     [/yr]
      flvl = Enc/(Cmax + Enc)                           ! food limitation     [-]
      call checknan(flvl, nGrid)                        ! remove Nans

      Eavail  = epsAssim_vec*flvl*Cmax - metabolism         ! available energy    [/yr]

      grazing = Cmax * flvl*u                           ! grazing             [gWW/m2/yr]

      loss    = (1.d0-epsAssim_vec)*grazing + metabolism*u  ! total loss          [gWW/m2/yr]

! ----------------------------------------------
! Mortality for resources and fish grids:
! ----------------------------------------------

      mortpred = Cmax*V/(Enc + Cmax)*u
      call checknan(mortpred, nGrid)

      mortpred = matmul(transpose(theta), mortpred)      ! Predation mortality [/yr]

!             add basal and fishing mortality)
      mort = mortpred + mort0 + mortF                    ! Total mortality     [/yr]

! ----------------------------------------------
!  Flux out of the fish size group:
! ----------------------------------------------

      ! fish only data (fish grids)

      ii = 1
      do i = nResources+1, nGrid
        B(ii)        = u(i)                              ! fish stages          [g/m2]
        eFish(ii)    = Eavail(i)                         ! availabel energy     [/yr]
        eplus(ii)    = max(0d0, Eavail(i))               ! net growth rate      [/yr]
        mortFish(ii) = mort(i)                           ! total mortality      [/yr]
        ii = ii+1
      end do

      grow = (1.d0 - psiMature)*eplus                    ! energy  for growth   [/yr]

      gamma_vec = (grow - mortFish) /   &                    ! growth to next stage [/yr]
            (1d0 - (1/z)**(1d0-mortFish/grow) )

      call checknan(gamma_vec, nFGrid)    ! No growth of fully mature classes (grow=0)

      Fout = gamma_vec*B                                     ! flux out of stage    [g/m2/yr]

      Repro = psiMature*eplus*B                          ! reproduction         [g/m2/yr]

! ----------------------------------------------
! Flux into the size group
! ----------------------------------------------

      do i = 1, nGroups

      ! stages of this group in fish grid
        istart = ixStart(i) -nResources
        istop  = ixEnd(i) -nResources            ! last stage

        totRepro(i)    = repro(istart)
        totBiomass(i)  = B(istart)

        do j = istart+1, istop
          Fin(j)         = Fout(j-1)
          totRepro(i)    = totRepro(i) + Repro(j)
          totBiomass(i)  = totBiomass(i) + B(j)
        end do

        totRepro(i) = totRepro(i) + Fout(istop) ! growth out ouf final stage is reproduction
        Fin(istart) = epsRepro_vec(i)*totRepro(i)   ! reproduction

      ! stages of this group in total grid
        istart = ixStart(i) ! + nResources
        istop  = ixEnd(i)  ! + nResources

        totGrazing(i)  = 0d0
        totLoss(i)     = 0d0
        totMort(i)     = 0d0

        do j = istart, istop
          totGrazing(i)  = totGrazing(i) + grazing(j)
          totLoss(i)     = totLoss(i)    + loss(j)
          totMort(i)     = totMort(i)    + mort(j)*u(j)
        end do

      end do
      totRecruit   = totRepro*epsRepro_vec

! ----------------------------------------------
! Derivatives of fish:
! ----------------------------------------------
      dBdt = Fin - Fout + (eFish - mortFish)*B - Repro

! ----------------------------------------------
! Derivative of resources
! ----------------------------------------------
      ! resource only data
      do i = 1, nResources
        mortRes(i) = mort(i)              ! mortality rate [/year]
        R(i)       = u(i)                 ! resource [gWW/m2]
      enddo
      if (Rtype == 1) then
        dRdt = rr*(K-R) - mortRes*R       ! chemostat formulation
      else
        dRdt = rr*R*(1-R/K) - mortRes*R   ! logistic formulation
      end if

      do i = 1, nResources
        dudt(i) = dRdt(i)
      enddo

      do i = 1, nFGrid
        dudt(i+nResources) = dBdt(i)
      enddo

  end subroutine calcderivatives
!============================================================



end module FEISTY




!============================================================
! New stuff based on Karline Soetaert package (Oct 2023 added)
!============================================================

!==========================================================================
!==========================================================================
! interface subroutines between R and fortran for thefeisty model
!==========================================================================
!==========================================================================

  subroutine initfeisty (steadyparms)
    use feisty
    implicit none
    external steadyparms  ! not used
    real(dp) :: pars

       feistyinitialised = .FALSE.

   end subroutine initfeisty

   subroutine initfeistysetupbasic (steadyparms)
    use feisty
    implicit none
    external steadyparms  ! not used
    real(dp) :: parmsbasic(5)

       feistyinitialised = .FALSE.

       call steadyparms(5,parmsbasic)
       call setupbasic(parmsbasic(1),parmsbasic(2),parmsbasic(3),parmsbasic(4),parmsbasic(5))

       feistyinitialised = .TRUE.

   end subroutine initfeistysetupbasic

   subroutine initfeistysetupbasic2 (steadyparms)
    use feisty
    implicit none
    external steadyparms  ! not used
    real(dp) :: parmsbasic(7)

       feistyinitialised = .FALSE.

       call steadyparms(7,parmsbasic)
       call setupbasic2(parmsbasic(1),parmsbasic(2),parmsbasic(3),INT(parmsbasic(4)),parmsbasic(5),&
                       parmsbasic(6),parmsbasic(7))

       feistyinitialised = .TRUE.

   end subroutine initfeistysetupbasic2

   subroutine initfeistysetupVertical (steadyparms)
    use feisty
    implicit none
    external steadyparms  ! not used
    real(dp) :: parmsbasic(8)

       feistyinitialised = .FALSE.

       call steadyparms(8,parmsbasic)
       call setupVertical(parmsbasic(1),parmsbasic(2),parmsbasic(3),INT(parmsbasic(4)),INT(parmsbasic(5)),&
                       parmsbasic(6),parmsbasic(7),parmsbasic(8))

       feistyinitialised = .TRUE.

   end subroutine initfeistysetupVertical

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the feisty model
!==========================================================================
!==========================================================================

   subroutine runfeisty (neq, t, Conc, dConc, yout, ip)
    use feisty
    implicit none

    integer,  intent(in):: neq, ip(*)
    real(dp), intent(in):: t, conc(neq)
    real(dp), intent(inout):: yout(*)
    real(dp), intent(out):: dconc(neq)
    integer              :: i
!..........................................................................

      if (.NOT. feistyinitialised) then
         call allocfeisty(ip, yout)
         feistyinitialised = .TRUE.
      end if

      call calcderivatives(conc, dconc)
      CALL outfeisty(yout)

   end subroutine runfeisty

!==========================================================================

  subroutine outfeisty(yout)
    use feisty
    implicit none

    real(dp), intent(out):: yout(*)
    integer:: ir, i

!..........................................................................
    ir = 1
    do i = 1, nGrid
     yout(ir) = flvl(i)
     ir = ir + 1
    end do
    do i = 1, nGrid
     yout(ir) = mortpred(i)
     ir = ir + 1
    end do
    do i = 1, nFGrid
     yout(ir) = grow(i)
     ir = ir + 1
    end do
    do i = 1, nFGrid
     yout(ir) = Repro(i)
     ir = ir + 1
    end do
    do i = 1, nFGrid
     yout(ir) = Fin(i)
     ir = ir + 1
    end do
    do i = 1, nFGrid
     yout(ir) = Fout(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totMort(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totGrazing(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totLoss(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totRepro(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totRecruit(i)
     ir = ir + 1
    end do
    do i = 1, nGroups
     yout(ir) = totBiomass(i)
     ir = ir + 1
    end do

   end subroutine outfeisty
!============================================================
