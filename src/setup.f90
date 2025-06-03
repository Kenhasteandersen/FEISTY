!
! FEISTY model
! References: Petrik et al., 2019; van Denderen et al., 2020.
! The library follows MATLAB/R codes from Ken H. Andersen; P. Daniël van Denderen; Rémy Denéchère; Daniel Ottmann Riera ...
! Programmed by Yixin Zhao, August 2022.
!
Module setup
   use globals
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

   real(rk) :: martin !  martin curve depth
   real(rk), allocatable:: sizeprefer(:, :)        ! Size preference matrix
   real(rk), allocatable:: vertover(:, :)          ! Vertical overlap matrix
   real(rk), allocatable:: theta(:, :)             ! Feeding preference matrix
   real(rk), dimension(:), allocatable :: F         ! Available food

! assembled vectors including resources and fish
   real(rk), allocatable:: V(:)          ! assembled vector of clearance rate
   real(rk), allocatable:: Enc(:)        ! assembled vector of encounter rate
   real(rk), allocatable:: flvl(:)       ! assembled vector of feeding level
   real(rk), allocatable:: Cmax(:)       ! assembled vector of maximum consumption rate
   real(rk), allocatable:: mortpred(:)   ! assembled vector of predation mortality
   real(rk), allocatable:: mc(:)         ! assembled vector of central mass
   real(rk), allocatable:: mL(:)         ! assembled vector of lower mass
   real(rk), allocatable:: mU(:)         ! assembled vector of upper mass

! resource parameters    input from input file
   real(rk), allocatable:: K(:), rr(:)   ! Carrying capacity of resources and growth rate of resources
   
! predation preference coefficient
   real(rk) :: thetaS
   real(rk) :: thetaA
   real(rk) :: thetaD
!  real(rk),parameter ::   thetaS = 0.25_rk ! Medium fish pref for small zooplankton
!  real(rk),parameter ::   thetaA = 0.5_rk  ! Large fish pref for medium forage fish
!  real(rk),parameter ::   thetaD = 0.75_rk ! Pref of large demersal on pelagic prey

!=======================================================================================================
! New from Karline Soetaert package   (Oct 2023 added)
         integer:: nFGrid        ! Total number of fish size grids

         integer :: Rtype   ! resource dynamics, 1=chemostat; 2=logistic

         ! vectors defined in all stages
         real(rk), allocatable:: epsAssim_vec(:)   ! assimilation efficiency
         real(rk), allocatable:: metabolism(:) ! basal respiration rate
         real(rk), allocatable:: mort0(:)      ! basal mortality rate
         real(rk), allocatable:: mortF(:)      ! fishing mortality rate
         ! vectors defined only in fish stages
         real(rk), allocatable:: z(:)          ! ratio of mass in fish size class
         real(rk), allocatable:: psiMature(:)  ! fraction mature in fish size class
         ! vectors defined in fish groups
         real(rk), allocatable:: epsRepro_vec(:)   ! efficiency of reproduction

         real(rk), allocatable:: grazing(:)    ! grazing rate
         real(rk), allocatable:: loss(:)       ! summed loss rates
         real(rk), allocatable:: mort(:)       ! summed mortaliy rates
         real(rk), allocatable:: eAvail(:)     ! available energy

         ! Defined in fish grids only
         real(rk), allocatable:: B(:), dBdt(:) ! state variables and derivatives
         real(rk), allocatable:: eplus(:), eFish(:) ! available energy
         real(rk), allocatable:: grow(:)       ! energy for growth [/yr]
         real(rk), allocatable:: gamma_vec(:)      ! growth to next stage
         real(rk), allocatable:: Repro(:)      ! reproduction rates
         real(rk), allocatable:: mortFish(:)   ! mortality
         real(rk), allocatable:: Fout(:), Fin(:) ! flux in and out of stages

         ! Defined in fish groups only
         real(rk), allocatable:: totMort(:)    ! Total mortality
         real(rk), allocatable:: totGrazing(:) ! Total grazing
         real(rk), allocatable:: totLoss(:)    ! Total losses
         real(rk), allocatable:: totRepro(:)   ! Total reproductive losses
         real(rk), allocatable:: totRecruit(:) ! Total recruitment
         real(rk), allocatable:: totBiomass(:) ! Total biomass

         ! Defined for resources
         real(rk), allocatable:: R(:), dRdt(:) ! state variables and derivatives
         real(rk), allocatable:: mortRes(:)    ! mortality rate pf resources

         ! Other
         real(rk):: tmp(12)
         logical :: feistyinitialised  = .FALSE.

!===================================
!for temperature effects (Feb 2024)
         integer, allocatable :: pelgroupidx(:)
         integer, allocatable :: demgroupidx(:)
         real(rk) :: pelagicT
         real(rk) :: benthicT
         integer, allocatable :: smdemidx(:)
         integer, allocatable :: lgdemidx(:)
         integer, allocatable :: pelRidx(:)
         integer(rk), allocatable :: pelgrididx(:), allgrididx(:)
         real(rk), allocatable:: metabolismsave(:)
         real(rk), allocatable:: Vsave(:)
         real(rk), allocatable:: Cmaxsave(:)
         logical :: bET
         real(rk) :: Q10ET
         real(rk) :: Q10mET
         real(rk) :: depthET

!for time-series input (Dec 2024)
         real(rk), allocatable :: dr_fac_theta(:, :)             ! down-regulation factor matrix
         real(rk) :: szprod, lzprod
         real(rk) :: dr_fac_sz, dr_fac_lz
         real(rk), allocatable :: forcs(:)
         integer :: nforcs
         logical :: bTS
         real(rk) :: smzcsp, lgzcsp, smzcsp_dr, lgzcsp_dr


contains
! ======================================
!  Different Setups
! ======================================

! --------------------------------------
! Setup according to Petrik et al. (2019)
! --------------------------------------
   subroutine setupbasic(szprod, lzprod, bprodin, dfbot, depth, Ts, Tb)
      real(rk), intent(in)::szprod,lzprod, bprodin, dfbot, depth, Ts, Tb ! bprodin: benthic productivity, dfbot: detrital flux reaching the sea floor
      real(rk) :: bprod                                                  ! only one of them works, keep the unused arguments negative e.g., bprodin = 100._rk, dfbot = -1._rk
      integer :: iGroup

      !benthic productivity
      if (bprodin < 0._rk .and. dfbot<0._rk) then
        bprod = 5._rk
      else if (bprodin .gt. 0._rk .and. dfbot.gt.0._rk) then
        stop
      else
       if (bprodin .ge. 0._rk) then
        bprod=bprodin
       end if
       if (dfbot .ge. 0._rk) then
        bprod=dfbot*0.1_rk
       end if
      end if

      !call read_namelist_setupbasic()
      call parametersInit(3, 2 + 3 + 3, 4, szprod,lzprod, bprod) ! (fish groups, total fish stages 2stages+3stages+3stages, 4 resources, szprod,lzprod, bprod,none)
      call parametersAddGroup(2, 2.50e2_rk, 0.5_rk) ! fishSmall
      call parametersAddGroup(3, 1.25e5_rk, 2.5e2_rk) ! fishLarge
      call parametersAddGroup(3, 1.25e5_rk, 2.5e2_rk)   ! fishDemersal

      ! vectors:
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))
      allocate (mL(nGrid))
      allocate (mU(nGrid))

      theta = 0._rk ! overwritten latter

! Overwrite
      do iGroup = 1, nGroups
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p) ! overwrite metabolism
         group(iGroup)%spec%metabolismsave =   group(iGroup)%spec%metabolism

         group(iGroup)%spec%psiMature = 0._rk ! reset
         group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5_rk ! only adults reproduce

!         if (group(iGroup)%spec%n .eq. 2) then
!            group(iGroup)%spec%mortF(group(iGroup)%spec%n) = 0.3_rk ! only adults have fishing mortality
!         else if ((group(iGroup)%spec%n .eq. 3)) then
!            group(iGroup)%spec%mortF(group(iGroup)%spec%n - 1) = 0.03_rk ! juvenile
!            group(iGroup)%spec%mortF(group(iGroup)%spec%n) = 0.3_rk ! adult
!         end if
      end do

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formmassvector(spec, iGroup, mc, mL, mU)
         end select
      end do
      mc(1:nResources) = [2.e-06_rk*sqrt(500._rk), 1.e-3_rk*sqrt(500._rk), 0.5e-03_rk*sqrt(250000._rk), 0.25_rk*sqrt(500._rk)] ! overwrite by resource mass
      !mU(1:nResources) = [2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)] ! weight central size
      !mL(1:nResources) = [2e-06,0.001, 0.5e-03, 0.25] ! weight lower limit)

      ! Small pelagics:
      theta(5, 1) = 1._rk ! Small ones eat only small zooplankton
      theta(6, 1:10) = [thetaS, 1._rk, 0._rk, 0._rk, 1._rk, 0._rk, 1._rk, 0._rk, 0._rk, 1._rk]

      ! Large pelagics:
      theta(7, 1) = 1._rk
      theta(8, 1:10) = [thetaS, 1._rk, 0._rk, 0._rk, 1._rk, 0._rk, 1._rk, 0._rk, 0._rk, 1._rk]
      theta(9, 6) = thetaA          ! medium forage fish
      theta(9, 8) = 1._rk            ! medium large pelagics

      ! Demersals:
      theta(10, 1) = 1._rk           ! larval demersals ear small zooplankton
      theta(11, 3) = 1._rk           ! medium demersals eat small benthos
      if (depth .lt. 200._rk) then
      theta(12, 6) = thetaA*thetaD  ! medium forage fish
      theta(12, 8) = thetaD         ! medium large pelagics
      end if
      theta(12, 3) = 1._rk           ! large demersals eat small benthos
      theta(12, 11) = 1._rk          ! medium demersals

    ! update temperature
        if(allocated(pelRidx)) deallocate (pelRidx)
        pelRidx=[1,2]! hard coded, used in effective Temperature
        call updateTemp(Ts,Tb, depth, [1,2],2,[3],1)
        bET = .TRUE.

    !update vector V Cmax for T effects

    call set2vec
    Rtype=1
    
!contains
!   subroutine read_namelist_setupbasic()
!      integer :: file_unit, io_err
!
!      namelist /input_setupbasic/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
!                                  & beta, sigma, mMin, &
!                                  & mMedium, mLarge, &
!                                  & lbenk, szoog, lzoog, sbeng, lbeng, &
!                                  & thetaS, thetaA, thetaD
!
!      call open_inputfile(file_unit, io_err)
!      read (file_unit, nml=input_setupbasic, iostat=io_err)
!      call close_inputfile(file_unit, io_err)
!   end subroutine read_namelist_setupbasic
   
   end subroutine setupbasic

! --------------------------------------
! Setup by Ken H. Andersen based on Petrik et al. (2019).
! --------------------------------------
   subroutine setupbasic2(szprod,lzprod, bprodin, dfbot, nStages, depth, Ts, Tb, etaMature,Fmax,etaF,bETin) !
      real(rk), intent(in) :: szprod,lzprod, bprodin, dfbot, Ts, Tb, depth ! bprodin: benthic productivity, dfbot: detrital flux reaching the sea floor
                                                                           ! only one of them works, keep the unused arguments negative e.g., bprodin = 100._rk, dfbot = -1._rk

      real(rk), intent(in) :: etaMature ! Mature mass relative to asymptotic size default 0.25, original in van Denderen et al., 2021 was 0.002
      real(rk), intent(in) :: Fmax, etaF ! fishing mortality, etaF * asymptotic size =fish size with 50% fishing mortality
      integer, intent(in) :: nStages,bETin
      real(rk) :: bprod
      integer :: iGroup, i, j

      !benthic productivity
      if (bprodin < 0._rk .and. dfbot<0._rk) then
        bprod = 5._rk
      else if (bprodin .gt. 0._rk .and. dfbot.gt.0._rk) then
        stop
      else
       if (bprodin .ge. 0._rk) then
        bprod=bprodin
       end if
       if (dfbot .ge. 0._rk) then
        bprod=dfbot*0.1_rk
       end if
      end if

      !call read_namelist_setupbasic2()    !
      call parametersInit(3, nint(0.66_rk*nStages) + nStages + nStages, 4, szprod,lzprod, bprod) ! (fish groups, total fish stages, 4 resources,szprod,lzprod, bprod,none)
      call parametersAddGroup(nint(0.66_rk*nStages), 2.50e2_rk, etaMature*2.50e2_rk) ! fishSmall  original mature mass is 0.002 *2.50e2_rk=0.5_rk
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishLarge (stages, max mass, mature mass) 0.002 *1.25e5_rk=2.5e2_rk
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishDemersal

      ! vectors:
      allocate (V(nGrid))
      allocate (Enc(nGrid))
      allocate (flvl(nGrid))
      allocate (Cmax(nGrid))
      allocate (mortpred(nGrid))
      allocate (mc(nGrid))
      allocate (mL(nGrid))
      allocate (mU(nGrid))

      theta = 0._rk ! overwritten latter

! Overwrite
      do iGroup = 1, nGroups
         group(iGroup)%spec%metabolism = (kk*group(iGroup)%spec%m**p)!overwrite metabolism
         group(iGroup)%spec%metabolismsave = group(iGroup)%spec%metabolism
      end do

! fishing mortality
      call setFishing(Fmax,etaF)

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formmassvector(spec, iGroup, mc, mL, mU)
         end select
      end do
      mc(1:nResources) = [2.e-06_rk*sqrt(500._rk), 1.e-3_rk*sqrt(500._rk), 1.e-4_rk*sqrt(250000._rk), 0.25_rk*sqrt(500._rk)] ! overwrite by resource mass
      !mU = c(2e-06*sqrt(500), 0.001*sqrt(500), 0.5e-03*sqrt(250000), 0.25*sqrt(500)) ! weight central size
      !mL = c(2e-06,0.001, 0.5e-03, 0.25) ! weight lower limit)

!!basic feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            theta(i, j) = exp(-(log(mc(i)/(beta*mc(j))))**2/(2*sigma)**2)
            if (mc(j) .gt. mc(i)) then                   !small can't eat large
               theta(i, j) = 0._rk
            end if
         end do
      end do

! FROM NUM Andersen, K. H., & Visser, A. W. (2023). Appendix https://doi.org/10.1016/j.pocean.2023.102995
! feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            theta(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
!            if (mc(j) .gt. mc(i)) then                   !small can't eat large
!               theta(i, j) = 0._rk
!            end if
!         end do
!      end do

! further clac theta : feeding selection in terms of fish kinds and resources
      ! Small pelagic
      theta(ixStart(1):ixEnd(1), 3:4) = 0._rk                  ! pelagic has not feeding on benthos
      do i = 1, group(3)%spec%n
      if (group(3)%spec%m(i) > mMedium .AND. &
          group(3)%spec%m(i) < mLarge) then
         theta(ixStart(1):ixEnd(1), ixStart(3) + i - 1) = 0._rk !small pelagic has not feeding on medium demersals
      end if
      end do
      ! Large pelagic
      theta(ixStart(2):ixEnd(2), 3:4) = 0._rk                  ! pelagic has not feeding on benthos
      do i = 1, group(3)%spec%n
      if (group(3)%spec%m(i) > mMedium .AND. &
          group(3)%spec%m(i) < mLarge) then
         theta(ixStart(2):ixEnd(2), ixStart(3) + i - 1) = 0._rk !large pelagic has not feeding on medium demersals
      end if
      end do
      ! Large pelagics have reduced feeding efficiency on small pelagics:
      theta(ixStart(2):ixEnd(2), ixStart(1):ixEnd(1)) = thetaA*theta(ixStart(2):ixEnd(2), ixStart(1):ixEnd(1))
      ! Demersals
      do i = 1, group(3)%spec%n

        ! small
         if (group(3)%spec%m(i) .le. mMedium) then
            theta(ixStart(3) + i - 1, 3:4) = 0._rk          ! Small demersals has no feeding on benthos
         end if

        ! medium
         if (group(3)%spec%m(i) > mMedium .AND. &
             group(3)%spec%m(i) < mLarge) then
            theta(ixStart(3) + i - 1, 1:2) = 0._rk          ! Medium demersals has no feeding on zooplankton

            theta(ixStart(3) + i - 1, ixStart(1):ixEnd(2)) = 0._rk   ! Medium demersals only eat benthos and do cannibalism
            do j = 1, group(3)%spec%n
                if (group(3)%spec%m(j) .le. mMedium .OR. &
                      group(3)%spec%m(j) .ge. mLarge) then
                    theta(ixStart(3) + i - 1, ixStart(3) + j - 1) = 0._rk
                end if
            end do
         end if

        ! Large demersals feed have reduced feeding efficiency on pelagic species:
         if (group(3)%spec%m(i) .ge. mLarge) then
           if(depth .lt. 200) then
      theta(ixStart(3) + i - 1, ixStart(1):ixEnd(1)) = thetaA*thetaD*theta(ixStart(3) + i - 1, ixStart(1):ixEnd(1))
      theta(ixStart(3) + i - 1, ixStart(2):ixEnd(2)) = thetaD*theta(ixStart(3) + i - 1, ixStart(2):ixEnd(2))
!      theta(ixStart(3) + i - 1, 1:2) =
!              do j = 1, group(3)%spec%n
!              if (group(3)%spec%m(j) .le. mMedium) then
!      theta(ixStart(3) + i - 1, ixStart(3) + j - 1 ) =
!              endif
!              enddo

           else

      theta(ixStart(3) + i - 1, ixStart(1):ixEnd(1)) = 0
      theta(ixStart(3) + i - 1, ixStart(2):ixEnd(2)) = 0
      theta(ixStart(3) + i - 1, 1:2) = 0
              do j = 1, group(3)%spec%n
              if (group(3)%spec%m(j) .le. mMedium) then
      theta(ixStart(3) + i - 1, ixStart(3) + j - 1 ) = 0
              endif
              enddo
           endif
         end if

      end do

    !   theta(:,4) = 0._rk

    ! update temperature
        if(allocated(pelRidx)) deallocate (pelRidx)
        pelRidx=[1,2]! hard coded, used in effective Temperature
        call updateTemp(Ts, Tb, depth, [1,2],2,[3],1)
        if(bETin .eq. 1) bET = .TRUE.
        if(bETin .eq. 0) bET = .FALSE.

    call set2vec
    Rtype=1

!contains
   !subroutine read_namelist_setupbasic2()
   !   integer :: file_unit, io_err
   !
   !   namelist /input_setupbasic2/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
   !                               & beta, sigma, mMin, &
   !                               & mMedium, mLarge, &
   !                               & lbenk, szoog, lzoog, sbeng, lbeng, &
   !                               & thetaA, thetaD
   !
   !   call open_inputfile(file_unit, io_err)
   !   read (file_unit, nml=input_setupbasic2, iostat=io_err)
   !   call close_inputfile(file_unit, io_err)
   !end subroutine read_namelist_setupbasic2

   end subroutine setupbasic2
! --------------------------------------
! Setup of vertical overlap (van Denderen et al., 2020)
! --------------------------------------
   subroutine setupVertical(szprod,lzprod, bprodin, dfbot, dfpho, region, bottom, photic)
     !  default bottom:1500m euphotic depth 150m
      real(rk), intent(in) :: szprod,lzprod, bottom, bprodin, dfbot, dfpho, photic ! bprodin: benthic productivity, dfbot: detrital flux reaching the sea floor, dfpho: detrital flux out of the photic zone
                                                                                   ! only one of them works, keep the unused arguments negative e.g., bprodin = -1._rk, dfbot = -1._rk, dfpho = 100._rk
      integer, intent(in) :: region!, nStages                    ! Mature mass relative to asymptotic size default 0.25, original in van Denderen et al., 2021 was 0.002

! for theta calc
       real(rk) :: ssigma
       real(rk) :: tau
       !real(rk) :: bottom
       !real(rk) :: photic
       real(rk) :: shelfdepth
       real(rk) :: visual
!      real(rk) :: ssigma = 10._rk
!      real(rk) :: tau = 10._rk
!      real(rk), parameter :: bottom = 1500._rk ! total depth meter
!      real(rk), parameter :: photic = 150._rk  ! photic zone depth
!      real(rk), parameter :: shelfdepth = 250._rk   ! depth ?
!      real(rk), parameter :: visual = 1.5_rk ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
      real(rk), allocatable :: sigmap(:) ! width for each size class
      !real(rk) :: bent ! for bprod calc
      real(rk) :: bprod
      real(rk), dimension(:), allocatable :: xrange
      real(rk) :: dvm  ! vertical migration depth photic + 500._rk
      real(rk) :: xloc ! vertical location    will be overwritten again and again
      real(rk), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(rk), allocatable :: zp_n(:, :), zp_d(:, :), & ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :)        ! demersal night/ day
      real(rk) :: demmig ! demersal migration
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(rk), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
      real(rk),allocatable :: sizes(:)
      integer :: iGroup, i, j, ixmedium, ixlarge, nsize, matstageS, matstageL
      real(rk), parameter :: etaMature = 0.002_rk
      integer, parameter :: nStages = 6

      !call read_namelist_setupvertical()
      allocate(xrange(int(bottom) + 1))

! calc bprod before initialization
      if (bprodin < 0._rk .and. dfbot < 0._rk .and. dfpho < 0._rk) then
        !dfpho = 150._rk
        bprod = 0.1_rk*(150._rk*(bottom/photic)**(-0.86_rk)) ! from matlab
        if (bprod .ge. 150._rk*0.1_rk) then
           bprod = 150._rk*0.1_rk
        end if
      else if (bprodin .gt. 0._rk .and. dfbot .gt. 0._rk .and. dfpho .gt. 0._rk) then
        stop
      else
       if (bprodin >0._rk) then
        bprod=bprodin
       end if
       if (dfbot >0._rk) then
        bprod=dfbot*0.1_rk
       end if
       if (dfpho > 0._rk) then
        bprod = 0.1_rk*(dfpho*(bottom/photic)**(-0.86_rk)) ! from matlab
        if (bprod .ge. dfpho*0.1_rk) then
           bprod = dfpho*0.1_rk
        end if
       end if
      end if

      call parametersInit(5, nint(0.66_rk*nStages) + nint(0.66_rk*nStages) + nStages + nStages + nStages, 4, szprod,lzprod, bprod)!
      call parametersAddGroup(nint(0.66_rk*nStages), 2.50e2_rk, etaMature*2.50e2_rk) ! fishSmall  original mature mass is 0.002 *2.50e2_rk=0.5_rk
      call parametersAddGroup(nint(0.66_rk*nStages), 2.50e2_rk, etaMature*2.50e2_rk) ! fishMeso,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishLarge,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishBathy,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishDemersal,

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
      !allocate (vertover(nGrid, nGrid)) !

      vertover = 0._rk
      sizeprefer = 0._rk
      theta = 0._rk ! overwritten latter

! Overwrite
      do iGroup = 1, nGroups

         group(iGroup)%spec%metabolism = (0.2_rk*h*group(iGroup)%spec%m**p)
         group(iGroup)%spec%metabolismsave = group(iGroup)%spec%metabolism

         group(iGroup)%spec%psiMature = 0._rk ! reset
         !group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5_rk! only adults reproduce

      end do

      !overwrite psiMature    from matlab simple run
      nsize=nStages+1
      allocate (sizes(nsize))
      sizes = 10**(linspace(log10(mMin), log10(1.25e5_rk), nsize)) !      mMin=0.001     mMax=1.25e5_rk predatory fish
      matstageS = minloc(abs(sizes-0.5_rk),dim=1)! 0.002_rk*250._rk
      matstageL = minloc(abs(sizes-2.5e2_rk),dim=1)! 0.002_rk*125000._rk
      group(1)%spec%psiMature(matstageS:group(1)%spec%n) = 0.5_rk ! fishSmall
      group(2)%spec%psiMature(matstageS:group(2)%spec%n) = 0.5_rk ! fishMeso
      group(3)%spec%psiMature(matstageL:group(3)%spec%n) = 0.5_rk ! fishLarge
      group(4)%spec%psiMature(matstageL:group(4)%spec%n) = 0.5_rk ! fishBathy
      group(5)%spec%psiMature(matstageL:group(5)%spec%n) = 0.5_rk ! fishDemersal

      !group(nGroups)%spec%mortF(group(nGroups)%spec%n) = 0.5_rk ! only demersal adults have fishing mortality

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formmassvector(spec, iGroup,  mc, mL, mU)
         end select
      end do
      !from baseparameters.m
      mc(1:nResources) = [2.e-06_rk*sqrt(500._rk), 1.e-3_rk*sqrt(500._rk), 0.5e-03_rk*sqrt(250000._rk), 0.25_rk*sqrt(500._rk)] ! resource central mass
      mU(1:nResources) = [0.001_rk, 0.5_rk, 125._rk, 125._rk]  ! resource mass upper limit
      mL(1:nResources) = [2.e-6_rk, 0.001_rk, 0.5e-3_rk, 0.25_rk] ! resource mass lower limit
!! basic feeding preference matrix theta
      do i = idxF, nGrid
         do j = 1, nGrid
            sizeprefer(i, j) = sqrt(pi/2._rk)*sigma*( &
                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2._rk)*sigma)) &
                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2._rk)*sigma)))
            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
         end do
      end do

! FROM NUM Andersen, K. H., & Visser, A. W. (2023). Appendix https://doi.org/10.1016/j.pocean.2023.102995
! feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            sizeprefer(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
!            if (mc(j) .gt. mc(i)) then                   !small can't eat large
!               sizeprefer(i, j) = 0._rk
!            end if
!         end do
!      end do

!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0._rk, bottom, int(bottom) + 1)
      dvm = photic + 500._rk ! 650._rk
      !  from matlab
      if (bottom .lt. (photic + 500._rk)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      end if
      if (bottom .le. shelfdepth) then
         dvm = 0._rk                   ! no migration in shallow habitats
      end if

! first stages as medium/large for predators
      ixmedium = minloc(abs(sizes-0.5_rk),dim=1) ! 0.002_rk*250._rk
      ixlarge = minloc(abs(sizes-2.5e2_rk),dim=1) ! 0.002_rk*125000._rk

!      ixmedium = minloc(abs(mL(ixStart(5):ixEnd(5))-0.5_rk),dim=1) !predatory fish
!      ixlarge = minloc(abs(mL(ixStart(5):ixEnd(5))-etaMature*1.25e5_rk),dim=1)


  deallocate (sizes)!see above overwrite psimature

! zooplankton night
      allocate (zp_n(size(xrange), 2))
      xloc = 0._rk ! zoo on surface at night
      do i = 1, 2 !small zoo & large zoo
         zp_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(i)**2._rk)))* &
                      exp(-(((xrange - xloc)**2._rk)/(2._rk*sigmap(i)**2._rk)))
      end do
      zp_n = matmul(zp_n, diag(1._rk/sum(zp_n, 1)))

! zooplankton day (half at surface, half at dvm depth
      allocate (zp_d(size(xrange), 2))
      xloc = dvm !
      do i = 1, 2 !index: small zoo & large zoo
         zp_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(i)**2._rk)))* &
                      exp(-(((xrange - xloc)**2._rk)/(2._rk*sigmap(i)**2._rk)))
      end do
      zp_d = matmul(zp_d, diag(1._rk/sum(zp_d, 1)))
      zp_d = (zp_n + zp_d)/2._rk

! benthos small and large (at bottom with width sigma)
      allocate (bent_dn(size(xrange), 2))
      xloc = bottom
      do i = 1, 2 ! small benthos & large benthos
         bent_dn(:, i) = (1._rk/(sqrt(2._rk*pi*ssigma**2._rk)))*exp(-((xrange - xloc)**2._rk/(2._rk*ssigma**2._rk)))
         bent_dn(:, i) = bent_dn(:, i)/sum(bent_dn(:, i))
      end do

! small pelagic fish (day + night) always at surface
      allocate (ix(ixEnd(1) - ixStart(1) + 1))
      allocate (spel_dn(size(xrange), ixEnd(1) - ixStart(1) + 1))
      xloc = 0._rk
      ix = [(i, i=ixStart(1), ixEnd(1))]
      do i = 1, size(ix)
         spel_dn(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                         exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      spel_dn = matmul(spel_dn, diag(1._rk/sum(spel_dn, 1)))

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
         mpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      mpel_d = matmul(mpel_d, diag(1._rk/sum(mpel_d, 1)))

! large pelagic fish night (all at surface)
      allocate (lpel_n(size(xrange), ixEnd(3) - ixStart(3) + 1))
      deallocate (ix)
      allocate (ix(ixEnd(3) - ixStart(3) + 1))
      xloc = 0._rk
      ix = [(i, i=ixStart(3), ixEnd(3))]
      do i = 1, size(ix)
         lpel_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      lpel_n = matmul(lpel_n, diag(1._rk/sum(lpel_n, 1)))

! large pelagic fish day (non-large at surface   large at dvm)
      allocate (lpel_d(size(xrange), ixEnd(3) - ixStart(3) + 1))
      allocate (xlocvec(ixEnd(3) - ixStart(3) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixlarge:size(xlocvec)) = dvm  !  non-large at surface   large at dvm
      do i = 1, size(ix)
         lpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      lpel_d = matmul(lpel_d, diag(1._rk/sum(lpel_d, 1)))
      lpel_d = (lpel_d + lpel_n)/2._rk

! bathypelagic night (large in midwater, others at surface)
      allocate (bpel_n(size(xrange), ixEnd(4) - ixStart(4) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(4) - ixStart(4) + 1))
      allocate (ix(ixEnd(4) - ixStart(4) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixlarge:size(xlocvec)) = dvm  !  non-large at surface   large at dvm
      ix = [(i, i=ixStart(4), ixEnd(4))]
      do i = 1, size(ix)
         bpel_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      bpel_n = matmul(bpel_n, diag(1._rk/sum(bpel_n, 1)))

! bathypelagic day (all at dvm)
      allocate (bpel_d(size(xrange), ixEnd(4) - ixStart(4) + 1))
      xlocvec = dvm ! overwrite all elements by dvm
      do i = 1, size(ix)
         bpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      bpel_d = matmul(bpel_d, diag(1._rk/sum(bpel_d, 1)))

! demersal fish night
      allocate (dem_n(size(xrange), ixEnd(5) - ixStart(5) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(5) - ixStart(5) + 1))
      allocate (ix(ixEnd(5) - ixStart(5) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixmedium:size(xlocvec)) = bottom  !  small at surface   medium and large at bottom
      ix = [(i, i=ixStart(5), ixEnd(5))]
      do i = 1, size(ix)
         dem_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                       exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      dem_n = matmul(dem_n, diag(1._rk/sum(dem_n, 1)))

! demersal fish day
      demmig = dvm ! ??? from matlab
      if ((bottom - dvm) .ge. 1200._rk) then
         demmig = dvm + (bottom - dvm - 1200._rk)
      end if
      if ((bottom - dvm) .ge. 1500._rk) then
         demmig = bottom
      end if
      allocate (dem_d(size(xrange), ixEnd(5) - ixStart(5) + 1))
      xlocvec(ixlarge:size(xlocvec)) = demmig !=dvm? ! small at surface/ medium at bottom/ large and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                       exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      dem_d = matmul(dem_d, diag(1._rk/sum(dem_d, 1)))
! from matlab
      ! if shallower than euphotic depth, large demersals feed across-habitats
      if (bottom .le. photic) then
         dem_d = (dem_d + dem_n)/2._rk
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

      dayout = 0._rk
      test = 0._rk
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

      nightout = 0._rk
      test = 0._rk !reset
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
      nightout(visualpred, :) = nightout(visualpred, :)*(2._rk - visual) ! predation decrease at night

! pelagic predators have limited vision in twilight zone during day
      pelpred = [(i, i=ixStart(3), ixEnd(3))]   ! large pelagic   9 10 11
      pelpred = pelpred(ixlarge:size(pelpred)) ! large large pelagic  11  at dvm during day
      preytwi = [(i, i=ixStart(2), ixEnd(2)), (i, i=ixStart(4), ixEnd(4))] ! mesopelagic 7 8   bathypelagic 12 13 14
      dayout(pelpred, preytwi) = dayout(pelpred, preytwi)/visual*(2._rk - visual)    ! /visual to restore and then *0.5

!!! average overlap during the whole day
      vertover = (dayout + nightout)*0.5_rk
!!! calculate combined feeding preference matrix
      theta = sizeprefer*vertover

!! specific revision of feeding preference
!
      idx_be = [(i, i=idxF, ixStart(5) + (ixmedium - 2))]     ! all pelagic and larval demersals
      theta(idx_be, 3:4) = 0._rk      ! all pelagic and larval demersals do not eat benthos,
      ! only medium & large demersals eat benthos
! small demersals are less preyed on
      idx_smd = [(i, i=ixStart(5) + (ixmedium - 1), ixStart(5) + (ixlarge - 2))] ! medium demersal is at bottom
      theta(idx_be, idx_smd) = theta(idx_be, idx_smd)*0.25_rk
! medium & large demersals do not eat zooplankton
      theta(ixStart(5) + (ixmedium - 1):ixEnd(5), 1:2) = 0._rk
! provide benefit to forage and mesopelagic fish (predator avoidance)
      pred1 = [(i, i=ixStart(3) + (ixlarge - 1), ixEnd(3))]
      pred2 = [(i, i=ixStart(4) + (ixlarge - 1), ixEnd(4))]
      pred3 = [(i, i=ixStart(5) + (ixlarge - 1), ixEnd(5))]
      prey1 = [(i, i=ixStart(1) + (ixmedium - 1), ixEnd(1))]
      prey2 = [(i, i=ixStart(2) + (ixmedium - 1), ixEnd(2))]
      idx_predat = [pred1, pred2, pred3]
      idx_prey = [prey1, prey2]
      theta(idx_predat, idx_prey) = theta(idx_predat, idx_prey)*0.5_rk

    ! update temperature
    call updateTempV(depthDay, depthNight, bottom, region)
    ! all fish group
    do iGroup = 1, nGroups
        group(iGroup)%spec%V=group(iGroup)%spec%V*fTempV(ixStart(iGroup):ixEnd(iGroup))
        group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTempV(ixStart(iGroup):ixEnd(iGroup))
        group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempmV(ixStart(iGroup):ixEnd(iGroup))
    end do

      !vector
    call set2vec
    Rtype=1

!contains
!   subroutine read_namelist_setupvertical()
!      integer :: file_unit, io_err
!
!      namelist /input_setupvertical/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
!                                  & beta, sigma, mMin, &
!                                  & lbenk, szoog, lzoog, sbeng, lbeng,&!bent,
!                                  & mMedium, mLarge, &
!                                  & ssigma, tau, shelfdepth, visual
!
!      call open_inputfile(file_unit, io_err)
!      read (file_unit, nml=input_setupvertical, iostat=io_err)
!      call close_inputfile(file_unit, io_err)
!   end subroutine read_namelist_setupvertical

   end subroutine setupVertical

! --------------------------------------
! Revised setup of vertical overlap based on van Denderen et al. (2020)
! --------------------------------------
   subroutine setupVertical2(szprod,lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, bottom, photic, etaMature,&
                             shelfdepth, visual, Fmax, etaF)
     !  default bottom:1500m euphotic depth 150m
      real(rk), intent(in) :: szprod,lzprod, bottom, photic, bprodin, dfbot, dfpho, Tp, Tm, Tb, etaMature,Fmax,etaF ! bprodin: benthic productivity, dfbot: detrital flux reaching the sea floor, dfpho: detrital flux out of the photic zone
                                                                                                           ! only one of them works, keep the unused arguments negative e.g., bprodin = -1._rk, dfbot = -1._rk, dfpho = 100._rk
      integer, intent(in) :: nStages                                ! Mature mass relative to asymptotic size default 0.25, original in van Denderen et al., 2021 was 0.002
      real(rk), intent(in) :: shelfdepth  ! shelf region depth, meso- and bathy-elagic fish exist when water is deeper than this. default 250m
      real(rk), intent(in) :: visual      ! >1 visual predation primarily during the day, = 1 equal day and night. default 1.5

! for theta calc
       real(rk) :: ssigma
       real(rk) :: tau
       !real(rk) :: bottom
       !real(rk) :: photic
       !real(rk) :: mesop !change to shelfdepth
       !real(rk) :: visual
!      real(rk) :: ssigma = 10._rk
!      real(rk) :: tau = 10._rk
!      real(rk), parameter :: bottom = 1500._rk ! total depth meter
!      real(rk), parameter :: photic = 150._rk  ! photic zone depth
!      real(rk), parameter :: mesop = 250._rk   ! depth ?
!      real(rk), parameter :: visual = 1.5_rk ! scalar; >1 visual predation primarily during the day, = 1 equal day and night
      real(rk), allocatable :: sigmap(:) ! width for each size class
      !real(rk) :: bent ! for bprod calc
      real(rk) :: bprod
      real(rk), dimension(:), allocatable :: xrange
      real(rk) :: dvm  ! vertical migration depth photic + 500._rk
      real(rk) :: xloc ! vertical location    will be overwritten again and again
      real(rk), allocatable :: xlocvec(:) ! vertical location vector used for some species
      real(rk), allocatable :: zp_n(:, :), zp_d(:, :), & ! zooplankton day / night
                               bent_dn(:, :), spel_dn(:, :), & ! benthos day & night      small pelagic day & night
                               mpel_n(:, :), mpel_d(:, :), &   ! mesopelagic night / day
                               lpel_n(:, :), lpel_d(:, :), &   ! large pelagic night / day
                               bpel_n(:, :), bpel_d(:, :), &   ! bathypelagic night/ day
                               dem_n(:, :), dem_d(:, :)        ! demersal night/ day
      real(rk) :: demmig ! demersal migration
      integer, allocatable :: ix(:), idx_be(:), idx_smd(:), pred1(:), pred2(:), pred3(:), prey1(:), prey2(:), &
                              idx_predat(:), idx_prey(:)
      real(rk), allocatable :: depthDay(:, :), dayout(:, :), depthNight(:, :), nightout(:, :), test(:, :)
      integer, allocatable :: visualpred(:), pelpred(:), preytwi(:)
!      real(rk),allocatable :: sizes(:)
      integer :: iGroup, i, j, ixmedium, ixlarge!, nsize, matstageS, matstageL
      real(rk), allocatable :: Teff (:)
      real(rk) :: Tday, Tnight, Tdaylarge, Tdaynonlarge, &
                  Tsmall, Tmedium, Tnightlarge, Tnightnonlarge


      !call read_namelist_setupvertical()
      allocate(xrange(int(bottom) + 1))

! calc bprod before initialization
      if (bprodin < 0._rk .and. dfbot < 0._rk .and. dfpho < 0._rk) then
        !dfpho = 150._rk
        bprod = 0.1_rk*(150._rk*(bottom/photic)**(-0.86_rk)) ! from matlab
        if (bprod .ge. 150._rk*0.1_rk) then
           bprod = 150._rk*0.1_rk
        end if
      else if (bprodin .gt. 0._rk .and. dfbot .gt. 0._rk .and. dfpho .gt. 0._rk) then
        stop
      else
       if (bprodin >0._rk) then
        bprod=bprodin
       end if
       if (dfbot >0._rk) then
        bprod=dfbot*0.1_rk
       end if
       if (dfpho > 0._rk) then
        bprod = 0.1_rk*(dfpho*(bottom/photic)**(-0.86_rk)) ! from matlab
        if (bprod .ge. dfpho*0.1_rk) then
           bprod = dfpho*0.1_rk
        end if
       end if
      end if

      call parametersInit(5, nint(0.66_rk*nStages) + nint(0.66_rk*nStages) + nStages + nStages + nStages, 4, szprod,lzprod, bprod)!
      call parametersAddGroup(nint(0.66_rk*nStages), 2.50e2_rk, etaMature*2.50e2_rk) ! fishSmall  original mature mass is 0.002 *2.50e2_rk=0.5_rk
      call parametersAddGroup(nint(0.66_rk*nStages), 2.50e2_rk, etaMature*2.50e2_rk) ! fishMeso,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishLarge,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishBathy,
      call parametersAddGroup(nStages, 1.25e5_rk, etaMature*1.25e5_rk) ! fishDemersal,

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
      !allocate (vertover(nGrid, nGrid)) !

      vertover = 0._rk
      sizeprefer = 0._rk
      theta = 0._rk ! overwritten latter

! Overwrite
      do iGroup = 1, nGroups

         group(iGroup)%spec%metabolism = (0.2_rk*h*group(iGroup)%spec%m**p)
         group(iGroup)%spec%metabolismsave = group(iGroup)%spec%metabolism

         !group(iGroup)%spec%psiMature = 0._rk ! reset
         !group(iGroup)%spec%psiMature(group(iGroup)%spec%n) = 0.5_rk! only adults reproduce

      end do

      !overwrite psiMature    from matlab simple run
!      nsize=nStages+1
!      allocate (sizes(nsize))
!      sizes = 10**(linspace(log10(mMin), log10(1.25e5_rk), nsize)) !      mMin=0.001     mMax=1.25e5_rk predatory fish
!      matstageS = minloc(abs(sizes-0.5_rk),dim=1)
!      matstageL = minloc(abs(sizes-2.5e2_rk),dim=1)
!      group(1)%spec%psiMature(matstageS:group(1)%spec%n) = 0.5_rk ! fishSmall
!      group(2)%spec%psiMature(matstageS:group(2)%spec%n) = 0.5_rk ! fishMeso
!      group(3)%spec%psiMature(matstageL:group(3)%spec%n) = 0.5_rk ! fishLarge
!      group(4)%spec%psiMature(matstageL:group(4)%spec%n) = 0.5_rk ! fishBathy
!      group(5)%spec%psiMature(matstageL:group(5)%spec%n) = 0.5_rk ! fishDemersal


! fishing mortality
      !group(nGroups)%spec%mortF(group(nGroups)%spec%n) = 0.5_rk ! only demersal adults have fishing mortality

      call setFishing(Fmax,etaF)

! Feeding preference matrix:
! assemble vectors
      do iGroup = 1, nGroups
         select type (spec => group(iGroup)%spec)
         type is (spectrumfish)
            call formmassvector(spec, iGroup, mc, mL, mU)
         end select
      end do
      !from baseparameters.m
      mc(1:nResources) = [2.e-06_rk*sqrt(500._rk), 1.e-3_rk*sqrt(500._rk), 1.e-4_rk*sqrt(250000._rk), 0.25_rk*sqrt(500._rk)] ! resource central mass
      mU(1:nResources) = [0.001_rk, 0.5_rk, 25._rk, 125._rk]  ! resource mass upper limit
      mL(1:nResources) = [2.e-6_rk, 0.001_rk, 1.e-4_rk, 0.25_rk] ! resource mass lower limit
!! basic feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            sizeprefer(i, j) = sqrt(pi/2._rk)*sigma*( &
!                               erf((log(mU(j)) - log(mc(i)/beta))/(sqrt(2._rk)*sigma)) &
!                               - erf((log(mL(j)) - log(mc(i)/beta))/(sqrt(2._rk)*sigma)))
!            sizeprefer(i, j) = sizeprefer(i, j)/(log(mU(j)) - log(mL(j)))
!         end do
!      end do

      do i = idxF, nGrid
         do j = 1, nGrid
            sizeprefer(i, j) = exp(-(log(mc(i)/(beta*mc(j))))**2/(2*sigma)**2)
            if (mc(j) .gt. mc(i)) then                   !small can't eat large
               sizeprefer(i, j) = 0._rk
            end if
         end do
      end do

! FROM NUM Andersen, K. H., & Visser, A. W. (2023). Appendix https://doi.org/10.1016/j.pocean.2023.102995
! feeding preference matrix theta
!      do i = idxF, nGrid
!         do j = 1, nGrid
!            sizeprefer(i, j) = calcPhi(mc(i)/mc(j), beta, sigma,mU(i)/mL(i))
!            if (mc(j) .gt. mc(i)) then                   !small can't eat large
!               sizeprefer(i, j) = 0._rk
!            end if
!         end do
!      end do

!!!vertical overlap
      sigmap = ssigma + tau*log10(mc/mc(1)) ! width for each size class
      xrange = linspace(0._rk, bottom, int(bottom) + 1)
      dvm = photic + 500._rk ! 650._rk
      !  from matlab
      if (bottom .lt. (photic + 500._rk)) then
         dvm = bottom     ! migration to bottom in intermediate habitats
      end if
      if (bottom .le. shelfdepth) then
         dvm = 0._rk                   ! no migration in shallow habitats
      end if

! first stages as medium/large for predators
      !ixmedium = minloc(abs(sizes-0.5_rk),dim=1) !from matlab
      !ixlarge = minloc(abs(sizes-2.5e2_rk),dim=1)

!      ixmedium = minloc(abs(mL(ixStart(5):ixEnd(5))-etaMature*250._rk),dim=1) !predatory fish
!      ixlarge = minloc(abs(mL(ixStart(5):ixEnd(5))-etaMature*1.25e5_rk),dim=1)

      ixmedium = minloc(abs(mL(ixStart(5):ixEnd(5))-0.5_rk),dim=1) !predatory fish
      ixlarge = minloc(abs(mL(ixStart(5):ixEnd(5))-250._rk),dim=1)


!  deallocate (sizes)!see above overwrite psimature

! zooplankton night
      allocate (zp_n(size(xrange), 2))
      xloc = 0._rk ! zoo on surface at night
      do i = 1, 2 !small zoo & large zoo
         zp_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(i)**2._rk)))* &
                      exp(-(((xrange - xloc)**2._rk)/(2._rk*sigmap(i)**2._rk)))
      end do
      zp_n = matmul(zp_n, diag(1._rk/sum(zp_n, 1)))

! zooplankton day (half at surface, half at dvm depth
      allocate (zp_d(size(xrange), 2))
      xloc = dvm !
      do i = 1, 2 !index: small zoo & large zoo
         zp_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(i)**2._rk)))* &
                      exp(-(((xrange - xloc)**2._rk)/(2._rk*sigmap(i)**2._rk)))
      end do
      zp_d = matmul(zp_d, diag(1._rk/sum(zp_d, 1)))
      zp_d = (zp_n + zp_d)/2._rk

! benthos small and large (at bottom with width sigma)
      allocate (bent_dn(size(xrange), 2))
      xloc = bottom
      do i = 1, 2 ! small benthos & large benthos
         bent_dn(:, i) = (1._rk/(sqrt(2._rk*pi*ssigma**2._rk)))*exp(-((xrange - xloc)**2._rk/(2._rk*ssigma**2._rk)))
         bent_dn(:, i) = bent_dn(:, i)/sum(bent_dn(:, i))
      end do

! small pelagic fish (day + night) always at surface
      allocate (ix(ixEnd(1) - ixStart(1) + 1))
      allocate (spel_dn(size(xrange), ixEnd(1) - ixStart(1) + 1))
      xloc = 0._rk
      ix = [(i, i=ixStart(1), ixEnd(1))]
      do i = 1, size(ix)
         spel_dn(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                         exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      spel_dn = matmul(spel_dn, diag(1._rk/sum(spel_dn, 1)))

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
         mpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      mpel_d = matmul(mpel_d, diag(1._rk/sum(mpel_d, 1)))

! large pelagic fish night (all at surface)
      allocate (lpel_n(size(xrange), ixEnd(3) - ixStart(3) + 1))
      deallocate (ix)
      allocate (ix(ixEnd(3) - ixStart(3) + 1))
      xloc = 0._rk
      ix = [(i, i=ixStart(3), ixEnd(3))]
      do i = 1, size(ix)
         lpel_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xloc)**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      lpel_n = matmul(lpel_n, diag(1._rk/sum(lpel_n, 1)))

! large pelagic fish day (non-large at surface   large at dvm)
      allocate (lpel_d(size(xrange), ixEnd(3) - ixStart(3) + 1))
      allocate (xlocvec(ixEnd(3) - ixStart(3) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixlarge:size(xlocvec)) = dvm  !  non-large at surface   large at dvm
      do i = 1, size(ix)
         lpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      lpel_d = matmul(lpel_d, diag(1._rk/sum(lpel_d, 1)))
      lpel_d = (lpel_d + lpel_n)/2._rk

! bathypelagic night (larges in midwater, others at surface)
      allocate (bpel_n(size(xrange), ixEnd(4) - ixStart(4) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(4) - ixStart(4) + 1))
      allocate (ix(ixEnd(4) - ixStart(4) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixlarge:size(xlocvec)) = dvm  !  non-large at surface   large at dvm
      ix = [(i, i=ixStart(4), ixEnd(4))]
      do i = 1, size(ix)
         bpel_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      bpel_n = matmul(bpel_n, diag(1._rk/sum(bpel_n, 1)))

! bathypelagic day (all at dvm)
      allocate (bpel_d(size(xrange), ixEnd(4) - ixStart(4) + 1))
      xlocvec = dvm ! overwrite all elements by dvm
      do i = 1, size(ix)
         bpel_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                        exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      bpel_d = matmul(bpel_d, diag(1._rk/sum(bpel_d, 1)))

! demersal fish night
      allocate (dem_n(size(xrange), ixEnd(5) - ixStart(5) + 1))
      deallocate (xlocvec)
      deallocate (ix)
      allocate (xlocvec(ixEnd(5) - ixStart(5) + 1))
      allocate (ix(ixEnd(5) - ixStart(5) + 1))
      xlocvec = 0._rk ! initialization
      xlocvec(ixmedium:size(xlocvec)) = bottom  !  small at surface   medium and large at bottom
      ix = [(i, i=ixStart(5), ixEnd(5))]
      do i = 1, size(ix)
         dem_n(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                       exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      dem_n = matmul(dem_n, diag(1._rk/sum(dem_n, 1)))

! demersal fish day
      demmig = dvm ! ??? from matlab
      if ((bottom - dvm) .ge. 1200._rk) then
         demmig = dvm + (bottom - dvm - 1200._rk)
      end if
      if ((bottom - dvm) .ge. 1500._rk) then
         demmig = bottom
      end if
      allocate (dem_d(size(xrange), ixEnd(5) - ixStart(5) + 1))
      xlocvec(ixlarge:size(xlocvec)) = demmig !=dvm? ! small at surface/ medium at bottom/ large and middle
      do i = 1, size(ix)
         dem_d(:, i) = (1._rk/(sqrt(2._rk*pi*sigmap(ix(i))**2._rk)))* &
                       exp(-((xrange - xlocvec(i))**2._rk/(2._rk*sigmap(ix(i))**2._rk)))
      end do
      dem_d = matmul(dem_d, diag(1._rk/sum(dem_d, 1)))
! from matlab
      ! if shallower than euphotic depth, large demersals feed across-habitats
      if (bottom .le. photic) then
         dem_d = (dem_d + dem_n)/2._rk
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

      dayout = 0._rk
      test = 0._rk
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

      nightout = 0._rk
      test = 0._rk !reset
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
      nightout(visualpred, :) = nightout(visualpred, :)*(2._rk - visual) ! predation decrease at night

! pelagic predators have limited vision in twilight zone during day
      pelpred = [(i, i=ixStart(3), ixEnd(3))]   ! large pelagic   9 10 11
      pelpred = pelpred(ixlarge:size(pelpred)) ! large large pelagic  11  at dvm during day
      preytwi = [(i, i=ixStart(2), ixEnd(2)), (i, i=ixStart(4), ixEnd(4))] ! mesopelagic 7 8   bathypelagic 12 13 14
      dayout(pelpred, preytwi) = dayout(pelpred, preytwi)/visual*(2._rk - visual)    ! /visual to restore and then *0.5

!!! average overlap during the whole day
      vertover = (dayout + nightout)*0.5_rk
!!! calculate combined feeding preference matrix
      theta = sizeprefer*vertover

!! specific revision of feeding preference
!
      idx_be = [(i, i=idxF, ixStart(5) + (ixmedium - 2))]     ! all pelagic and larval demersals
      theta(idx_be, 3:4) = 0._rk      ! all pelagic and larval demersals do not eat benthos,
      ! only medium & large demersals eat benthos
! small demersals are less preyed on
      idx_smd = [(i, i=ixStart(5) + (ixmedium - 1), ixStart(5) + (ixlarge - 2))] ! medium demersal is at bottom
      theta(idx_be, idx_smd) = theta(idx_be, idx_smd)*0.25_rk
! medium & large demersals do not eat zooplankton
      theta(ixStart(5) + (ixmedium - 1):ixEnd(5), 1:2) = 0._rk
! provide benefit to forage and mesopelagic fish (predator avoidance)
      pred1 = [(i, i=ixStart(3) + (ixlarge - 1), ixEnd(3))]
      pred2 = [(i, i=ixStart(4) + (ixlarge - 1), ixEnd(4))]
      pred3 = [(i, i=ixStart(5) + (ixlarge - 1), ixEnd(5))]
      prey1 = [(i, i=ixStart(1) + (ixmedium - 1), ixEnd(1))]
      prey2 = [(i, i=ixStart(2) + (ixmedium - 1), ixEnd(2))]
      idx_predat = [pred1, pred2, pred3]
      idx_prey = [prey1, prey2]
      theta(idx_predat, idx_prey) = theta(idx_predat, idx_prey)*0.5_rk

! ====================
! update temperature
! ====================
      if (allocated (fTempV)) then
        deallocate (fTempV)
        deallocate (fTempmV)
        !deallocate (Teff)
      end if

      allocate (fTempV(nGrid))
      allocate (fTempmV(nGrid))
      allocate (Teff(nGrid))

      fTempV  = 0._rk
      fTempmV = 0._rk
      Teff   = 0._rk

! zooplankton (no use)
      Tday = (Tp + Tm) / 2  ! half surface half dvm  dvm = photic + 500
      if (dvm == bottom) Tday = (Tp + Tb) / 2 ! when bottom < (photic + 500)
      if (dvm == 0) Tday = Tp  ! when bottom <= shelfdepth
      Tnight = Tp  ! all surface
      Teff(1:2) = (Tday + Tnight) / 2
! benthos (no use)
      Teff(3:4) = Tb
! small pelagics
      deallocate (ix)
      allocate (ix(ixEnd(1) - ixStart(1) + 1))
      ix = [(i, i=ixStart(1), ixEnd(1))]
      Teff(ix) = Tp ! always surface
! mesopelagics
      deallocate (ix)
      allocate (ix(ixEnd(2) - ixStart(2) + 1))
      ix = [(i, i=ixStart(2), ixEnd(2))]
      Tday = Tm ! dvm
      if (dvm == bottom) Tday = Tb
      if (dvm == 0) Tday = Tp
      Tnight = Tp ! surface
      Teff(ix) = (Tday + Tnight) / 2
! large pelagics
      deallocate (ix)
      allocate (ix(ixEnd(3) - ixStart(3) + 1))
      ix = [(i, i=ixStart(3), ixEnd(3))]
      ! daytime large half at surface half at dvm
      Tdaylarge = (Tp + Tm) / 2
      if (dvm == bottom) Tdaylarge = (Tp + Tb) / 2
      if (dvm == 0) Tdaylarge = Tp
      Tdaynonlarge = Tp  ! non-large at surface at daytime
      Tnight = Tp        ! all at surface at night
      Teff(ix(ixlarge:size(ix))) = (Tdaylarge + Tnight) / 2 ! large
      Teff(ix(1:(ixlarge-1))) = (Tdaynonlarge + Tnight) / 2 ! non-large
! bathypelagics
      deallocate (ix)
      allocate (ix(ixEnd(4) - ixStart(4) + 1))
      ix = [(i, i=ixStart(4), ixEnd(4))]
      Tday = Tm ! all at dvm at daytime
      if (dvm == bottom) Tday = Tb
      if (dvm == 0)      Tday = Tp
      Tnightlarge = Tm ! large at dvm
      if (dvm == bottom) Tnightlarge = Tb
      if (dvm == 0)      Tnightlarge = Tp
      Tnightnonlarge = Tp ! non-large at surface at night
      Teff(ix(ixlarge:size(ix))) = (Tday + Tnightlarge) / 2 ! large
      Teff(ix(1:ixlarge-1)) = (Tday + Tnightnonlarge) / 2 ! non-large
! demersals
      deallocate (ix)
      allocate (ix(ixEnd(5) - ixStart(5) + 1))
      ix = [(i, i=ixStart(5), ixEnd(5))]
      Tsmall = Tp ! small always at surface
      Tmedium = Tb ! medium always at bottom
      ! large
      ! daytime
      Tdaylarge = Tm ! large at middle
      ! if the water is very deep large demersals always stay at the bottom
      if ((bottom - dvm) >= 1500) Tdaylarge = Tb
      ! if the water is very shallow large demersals migrate over the whole water column both day and night
      if (bottom <= photic) then
        Tdaylarge = (Tp + Tb) / 2
      end if
      ! nighttime
      Tnightlarge = Tb ! large at bottom if water is deep enough
      ! if the water is very shallow large demersals migrate over the whole water column both day and night
      if (bottom <= photic) then
        Tnightlarge = (Tp + Tb) / 2
      end if

      Teff(ix(1:ixmedium-1)) = Tsmall ! small
      Teff(ix(ixmedium:ixlarge-1)) = Tmedium ! medium
      Teff(ix(ixlarge:size(ix))) = (Tdaylarge + Tnightlarge) / 2 ! large

      fTempV = Q10**((Teff - 10._rk) / 10._rk)
      fTempmV = Q10m**((Teff - 10._rk) / 10._rk)

    ! all fish group
    do iGroup = 1, nGroups
        group(iGroup)%spec%V=group(iGroup)%spec%V*fTempV(ixStart(iGroup):ixEnd(iGroup))
        group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmax*fTempV(ixStart(iGroup):ixEnd(iGroup))
        group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolism*fTempmV(ixStart(iGroup):ixEnd(iGroup))
    end do

      !vector
    call set2vec
    Rtype=1

!contains
!   subroutine read_namelist_setupvertical()
!      integer :: file_unit, io_err
!
!      namelist /input_setupvertical2/ h, nn, q, gamma, kk, p, epsAssim, epsRepro, &
!                                  & beta, sigma, mMin, &
!                                  & lbenk, szoog, lzoog, sbeng, lbeng,&!bent,
!                                  & mMedium, mLarge, &
!                                  & ssigma, tau!, mesop, visual
!
!      call open_inputfile(file_unit, io_err)
!      read (file_unit, nml=input_setupvertical2, iostat=io_err)
!      call close_inputfile(file_unit, io_err)
!   end subroutine read_namelist_setupvertical

   end subroutine setupVertical2

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
      real(rk), intent(in):: szprod,lzprod, bprod

      nGroups = nnGroups                   ! fish size spectrum group numbers (species)
      iCurrentGroup = 0
      nResources = nnResources             ! resource numbers
      nGrid = nnGrid + nnResources         ! total grid numbers   resources + total fish stages
      idxF = nResources + 1                ! fish grid begins at...
      nFGrid=nGrid-nResources
!
! Allocate/deallocate variables:
!
        if (allocated(group)) deallocate(group)
        allocate(group(nGroups))

        if (allocated(ixStart)) deallocate(ixStart)
        allocate(ixStart(nGroups))

        if (allocated(ixEnd)) deallocate(ixEnd)
        allocate(ixEnd(nGroups))

        if (allocated(F)) deallocate(F)
        allocate(F(nGrid))

        if (allocated(theta)) deallocate(theta)
        allocate(theta(nGrid, nGrid))

        if (allocated(sizeprefer)) deallocate(sizeprefer)
        allocate(sizeprefer(nGrid, nGrid))

        if (allocated(vertover)) deallocate(vertover)
        allocate(vertover(nGrid, nGrid))

        if (allocated(V)) deallocate(V)
        if (allocated(Enc)) deallocate(Enc)
        if (allocated(flvl)) deallocate(flvl)
        if (allocated(Cmax)) deallocate(Cmax)
        if (allocated(mortpred)) deallocate(mortpred)
        if (allocated(mc)) deallocate(mc)
        if (allocated(mL)) deallocate(mL)
        if (allocated(mU)) deallocate(mU)

        if (allocated(K)) deallocate(K)
        allocate (K(nResources))

        if (allocated(rr)) deallocate(rr)
        allocate (rr(nResources))

      ! define resources:
      K = [szprod, lzprod, bprod, lbenk]    ! Carrying capacity of resources [g m-2]]
      rr = [szoog, lzoog, sbeng, lbeng]   ! growth rate of resources       [yr-1]

         !maybe more variables

!-----------------Oct 2023 add-----------
        if (allocated(epsAssim_vec)) deallocate(epsAssim_vec)
        allocate(epsAssim_vec(nGrid))

        if (allocated(metabolism)) deallocate(metabolism)
        allocate(metabolism(nGrid))

        if (allocated(mort0)) deallocate(mort0)
        allocate(mort0(nGrid))

        if (allocated(mortF)) deallocate(mortF)
        allocate(mortF(nGrid))

        if (allocated(z)) deallocate(z)
        allocate(z(nFGrid))

        if (allocated(psiMature)) deallocate(psiMature)
        allocate(psiMature(nFGrid))

        if (allocated(epsRepro_vec)) deallocate(epsRepro_vec)
        allocate(epsRepro_vec(nGroups))

        if (allocated(grazing)) deallocate(grazing)
        allocate(grazing(nGrid))

        if (allocated(loss)) deallocate(loss)
        allocate(loss(nGrid))

        if (allocated(mort)) deallocate(mort)
        allocate(mort(nGrid))

        if (allocated(Eavail)) deallocate(Eavail)
        allocate(Eavail(nGrid))

        if (allocated(B)) deallocate(B)
        allocate(B(nFGrid))

        if (allocated(dBdt)) deallocate(dBdt)
        allocate(dBdt(nFGrid))

        if (allocated(R)) deallocate(R)
        allocate(R(nResources))

        if (allocated(dRdt)) deallocate(dRdt)
        allocate(dRdt(nResources))

        if (allocated(mortRes)) deallocate(mortRes)
        allocate(mortRes(nResources))

        if (allocated(grow)) deallocate(grow)
        allocate(grow(nFGrid))

        if (allocated(eplus)) deallocate(eplus)
        allocate(eplus(nFGrid))

        if (allocated(mortFish)) deallocate(mortFish)
        allocate(mortFish(nFGrid))

        if (allocated(eFish)) deallocate(eFish)
        allocate(eFish(nFGrid))

        if (allocated(gamma_vec)) deallocate(gamma_vec)
        allocate(gamma_vec(nFGrid))

        if (allocated(Fout)) deallocate(Fout)
        allocate(Fout(nFGrid))

        if (allocated(Fin)) deallocate(Fin)
        allocate(Fin(nFGrid))

        if (allocated(Repro)) deallocate(Repro)
        allocate(Repro(nFGrid))

        if (allocated(totMort)) deallocate(totMort)
        allocate(totMort(nGroups))

        if (allocated(totGrazing)) deallocate(totGrazing)
        allocate(totGrazing(nGroups))

        if (allocated(totLoss)) deallocate(totLoss)
        allocate(totLoss(nGroups))

        if (allocated(totRepro)) deallocate(totRepro)
        allocate(totRepro(nGroups))

        if (allocated(totRecruit)) deallocate(totRecruit)
        allocate(totRecruit(nGroups))

        if (allocated(totBiomass)) deallocate(totBiomass)
        allocate(totBiomass(nGroups))

!-----------------Feb 2024 add-----------
        if (allocated(metabolismsave)) deallocate(metabolismsave)
        allocate(metabolismsave(nGrid))

        if (allocated(Vsave)) deallocate(Vsave)
        allocate(Vsave(nGrid))

        if (allocated(Cmaxsave)) deallocate(Cmaxsave)
        allocate(Cmaxsave(nGrid))

        bET=.FALSE.

!-----------------Jan 2025 add-----------
        bTS=.FALSE.

!---------------------------------------

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
      real(rk), intent(in) :: mMax, mMature

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

   subroutine setFishing(Fmax,etaF)
     real(rk), intent(in) :: Fmax, etaF
     integer :: iGroup
     real(rk),allocatable :: psiF(:)
     real(rk) :: mFishing

!     if(F.eq.0._rk) then
!       RETURN
!     end if

     do iGroup = 1,nGroups
        allocate(psiF(group(iGroup)%spec%n))

      mFishing = etaF*maxval(group(iGroup)%spec%mUpper)
      psiF =( 1 + (group(iGroup)%spec%m/mFishing)**(-3) )**(-1)
      group(iGroup)%spec%mortF=psiF*Fmax
      deallocate(psiF)
     end do
   end subroutine


!-------------------------------------------------------------
! return assembled vectors containing values for all fish grid (no resources)
! resources mass must be assigned manually
   subroutine formmassvector(this, iGroup, mc, mL, mU)
      integer, intent(in) :: iGroup
      class(spectrumfish) :: this
      real(rk), intent(out) :: mc(nGrid), mL(nGrid), mU(nGrid)

      !V(ixStart(iGroup):ixEnd(iGroup)) = this%V
      !Cmax(ixStart(iGroup):ixEnd(iGroup)) = this%Cmax
      mc(ixStart(iGroup):ixEnd(iGroup)) = this%m
      mL(ixStart(iGroup):ixEnd(iGroup)) = this%mLower
      mU(ixStart(iGroup):ixEnd(iGroup)) = this%mUpper
   end subroutine formmassvector

   ! FROM NUM
    ! Calculate the interaction coefficient between two size groups.
    ! In:
    !   z : The predator:prey body mass ratio between the two groups
    !   beta: preferred predator:prey body mass ratio
    !   sigma: width of selection
    !   Delta: ratio between upper and lower body mass in size groups
    !
    function calcPhi(z, beta,sigma, Delta) result(res)
      real(rk), intent(in):: z,beta,sigma,Delta
      real(rk):: res, s

      if (beta .eq. 0._rk) then
         res = 0._rk ! beta = 0 is interpreted as if the group is not feeding
      else
         s = 2*sigma*sigma
         res = max(0._rk, &
         (Sqrt(Delta)*(((exp(-Log((beta*Delta)/z)**2/s) - 2/exp(Log(z/beta)**2/s) + &
         exp(-Log((Delta*z)/beta)**2/s))*s)/2. - &
         (Sqrt(Pi)*Sqrt(s)*(Erf((-Log(beta*Delta) + Log(z))/Sqrt(s))*Log((beta*Delta)/z) + &
         2*Erf(Log(z/beta)/Sqrt(s))*Log(z/beta) + &
         Erf((Log(beta) - Log(Delta*z))/Sqrt(s))*Log((Delta*z)/beta)))/2.))/ &
         ((-1 + Delta)*Log(Delta)) )
      end if
    end function calcPhi


! Oct 2023
! assign parameter set from Yixin Zhao library to vectors in Karline Soetaert package
subroutine set2vec
 integer :: iGroup

    epsAssim_vec=epsAssim
    epsRepro_vec=epsRepro

    V=0._rk
    Cmax=0._rk
    metabolism=0._rk
    mort0=0._rk
    mortF=0._rk
    psiMature=0._rk
    z=0._rk

    metabolismsave=0._rk
    Vsave=0._rk
    Cmaxsave=0._rk

    do iGroup=1,nGroups

    V( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%V
    Cmax( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%Cmax
    metabolism( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%metabolism
    mort0( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%mort0
    mortF( ixStart(iGroup) :ixEnd(iGroup) )    =   group(iGroup)%spec%mortF
    psiMature( ixStart(iGroup)-nResources : ixEnd(iGroup)-nResources )    =   group(iGroup)%spec%psiMature
    z( ixStart(iGroup)-nResources : ixEnd(iGroup)-nResources )   = group(iGroup)%spec%z

    metabolismsave( ixStart(iGroup) :ixEnd(iGroup) )=group(iGroup)%spec%metabolismsave
    Vsave( ixStart(iGroup) :ixEnd(iGroup) )=group(iGroup)%spec%Vsave
    Cmaxsave( ixStart(iGroup) :ixEnd(iGroup) )=group(iGroup)%spec%Cmaxsave

    end do

end subroutine set2vec


!--------------------------------------------------------------------
!     Temperature
!--------------------------------------------------------------------


   ! -----------------------------------------------
   ! Temperature Q10 function
   ! -----------------------------------------------
   function calfTemp(Q10, T) result(f) !calculate temperature factor
      real(rk), intent(in):: Q10, T
      real(rk):: f

      f = Q10**((T - Tref)/10._rk)
   end function calfTemp

! update temperature effects on physiological rates, setupbasic (Petrik et al., 2019) and setupbasic2
subroutine updateTemp(Tp, Tb, depth, pelgroup, npelgroup, demgroup, ndemgroup)
      real(rk), intent(in) :: Tp, Tb, depth
      integer, intent(in) :: npelgroup,ndemgroup
      integer, intent(in) ::  pelgroup(npelgroup), demgroup(ndemgroup)
      real(rk) :: eT, lambda
      real(rk), save :: Toldp = -1000._rk
      real(rk), save :: Toldb = -1000._rk
      integer, allocatable:: smdemidx_storage(:), lgdemidx_storage(:)  !temporary storage
      integer:: i,ii,j,iGroup

      if (allocated (smdemidx_storage))  deallocate (smdemidx_storage)
      if (allocated (lgdemidx_storage))  deallocate (lgdemidx_storage)
      if (allocated (smdemidx))  deallocate (smdemidx)
      if (allocated (lgdemidx))  deallocate (lgdemidx)

      if (allocated (pelgroupidx))    deallocate (pelgroupidx)
      allocate (pelgroupidx(npelgroup))
      pelgroupidx=pelgroup
      if (allocated (demgroupidx))    deallocate (demgroupidx)
      allocate (demgroupidx(ndemgroup))
      demgroupidx=demgroup
      pelagicT=Tp
      benthicT=Tb

      if (Tp .ne. Toldp .OR. Tb .ne. Toldb) then
         Toldp = Tp
         Toldb = Tb
         fTemp = calfTemp(Q10, Tp)  !Q10=1.88 clearance rate,  maximum consumption rate
         fTempm = calfTemp(Q10mPetrik, Tp) !Q10m=2.35 for metabolism    Petrik
         fTempdem = calfTemp(Q10, Tb)  !for demersal
         fTempmdem = calfTemp(Q10mPetrik, Tb) !for demersal

          !lambda should depend on B but temporarily set as 0.5
          lambda = 0.5_rk
          eT     = Tp * lambda + Tb * (1-lambda)
          fTempdem_shallow  = calfTemp(Q10, eT)
          fTempmdem_shallow = calfTemp(Q10mPetrik, eT)
      end if

    !all fish group
      !pelagic
    if (npelgroup.gt.0) then
     do i = 1, size(pelgroup)

        iGroup=pelgroup(i)

        group(iGroup)%spec%V=group(iGroup)%spec%Vsave *fTemp
        group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmaxsave *fTemp
        group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolismsave *fTempm
     end do
    end if

      !demersal
   if (ndemgroup.gt.0) then
    do ii= 1, size(demgroup)
        iGroup=demgroup(ii)

      ! form small & large demersal index array
        smdemidx_storage=[(i, i = ixStart(iGroup), ixEnd(iGroup))]!temporary
        smdemidx_storage=smdemidx_storage(PACK([(i, i=1, SIZE(mc(ixStart(iGroup):ixEnd(iGroup))))], &
                                           & mc(ixStart(iGroup):ixEnd(iGroup)) .le. mMedium))
        smdemidx        =smdemidx_storage
        lgdemidx_storage=[(i, i = ixStart(iGroup), ixEnd(iGroup))]!temporary
        lgdemidx_storage=lgdemidx_storage(PACK([(i, i=1, SIZE(mc(ixStart(iGroup):ixEnd(iGroup))))], &
                                           & mc(ixStart(iGroup):ixEnd(iGroup)) .ge. mLarge))
        lgdemidx        =lgdemidx_storage
         if (ii.ne.1) then
          smdemidx = [smdemidx,smdemidx_storage]
          lgdemidx = [lgdemidx,lgdemidx_storage]
         end if

     do i = 1, group(iGroup)%spec%n

       if (group(iGroup)%spec%m(i) .le. mMedium) then
      !small
         group(iGroup)%spec%V(i)=group(iGroup)%spec%Vsave(i) *fTemp
         group(iGroup)%spec%Cmax(i)=group(iGroup)%spec%Cmaxsave(i) *fTemp
         group(iGroup)%spec%metabolism(i)=group(iGroup)%spec%metabolismsave(i) *fTempm

       elseif (group(iGroup)%spec%m(i) .gt. mMedium .and. group(iGroup)%spec%m(i) .lt. mLarge)then
      !medium
         group(iGroup)%spec%V(i)=group(iGroup)%spec%Vsave(i) *fTempdem
         group(iGroup)%spec%Cmax(i)=group(iGroup)%spec%Cmaxsave(i) *fTempdem
         group(iGroup)%spec%metabolism(i)=group(iGroup)%spec%metabolismsave(i) *fTempmdem
       elseif (group(iGroup)%spec%m(i) .ge. mLarge)then
      !large
          if (depth .lt. 200._rk) then
          group(iGroup)%spec%V(i)=group(iGroup)%spec%Vsave(i) *fTempdem_shallow
          group(iGroup)%spec%Cmax(i)=group(iGroup)%spec%Cmaxsave(i) *fTempdem_shallow
          group(iGroup)%spec%metabolism(i)=group(iGroup)%spec%metabolismsave(i) *fTempmdem_shallow
          else
          group(iGroup)%spec%V(i)=group(iGroup)%spec%Vsave(i) *fTempdem
          group(iGroup)%spec%Cmax(i)=group(iGroup)%spec%Cmaxsave(i) *fTempdem
          group(iGroup)%spec%metabolism(i)=group(iGroup)%spec%metabolismsave(i) *fTempmdem
          end if
       end if
     end do
    end do
   end if

!   fTempold=fTemp
!   fTempmold=fTempm
!   fTempdemold=fTempdem
!   fTempmdemold=fTempmdem
!   fTempdem_shallowold=fTempdem_shallow
!   fTempmdem_shallowold=fTempmdem_shallow

   !prepare index for effective temperature
   if (allocated (pelgrididx)) deallocate (pelgrididx)
   if (allocated (allgrididx)) deallocate (allgrididx)

! assemble pelgrididx vector
   pelgrididx = pelRidx ! all pelagic resources (zoop)
   do ii= 1, size(pelgroupidx) ! add pelagic fish
       iGroup=pelgroupidx(ii)
        pelgrididx=[int(pelgrididx,4),(i, i = ixStart(iGroup), ixEnd(iGroup))]
   end do
   pelgrididx=[int(pelgrididx,4),smdemidx] ! add small demersal fish (pelagic)
   ! Note pelgrididx is integer(8)

   allgrididx=[(i, i = 1, nResources), (j, j = ixStart(1), ixEnd(nGroups))]

   Q10ET=Q10
   Q10mET=Q10mPetrik
   depthET=depth

end subroutine updateTemp


subroutine updateET(u)
    real(rk), intent(in) :: u(nGrid)
    real(rk) :: eT, lambda
    integer :: i,ii
    integer(8), allocatable :: pelpreyidx(:), allpreyidx(:)

   do ii= 1,size(lgdemidx)
    i=lgdemidx(ii)

    if (allocated (pelpreyidx)) deallocate (pelpreyidx)
    if (allocated (allpreyidx)) deallocate (allpreyidx)

    pelpreyidx = pack(pelgrididx, theta(i, pelgrididx) /= 0._rk)
    allpreyidx = pack(allgrididx, theta(i, allgrididx) /= 0._rk)

    lambda = sum(u(pelpreyidx)) / (sum(u(allpreyidx)) + eps) ! Eq. 15     eps = 1e-200_rk in case NA values generated 0/0
    eT = pelagicT * lambda + benthicT * (1 - lambda)
    fTempdem_shallow  = calfTemp(Q10ET, eT)
    fTempmdem_shallow = calfTemp(Q10mET, eT)

    !update vectors
    V(i)          = Vsave(i) * fTempdem_shallow
    Cmax(i)       = Cmaxsave(i) * fTempdem_shallow
    metabolism(i) = metabolismsave(i) * fTempmdem_shallow

   end do

end subroutine

! update Temperature for vertical version van Denderen et al., 2021
subroutine updateTempV(depthDay, depthNight, bottom, region)
 real(rk), intent(in) :: depthDay(:, :), depthNight(:, :), bottom
 integer :: i,region
 real(rk), allocatable :: dist(:,:), TQ10(:), TQ10m(:), fTemp_stepV(:,:), fTempm_stepV(:,:)
 real(rk), save :: tempdata(5501,5) ! contain tempdata from van Denderen et al., 2021 + default(10 celcius)

    if (allocated (fTempV)) then
        deallocate (fTempV)
        deallocate (fTempmV)
    end if

allocate (dist(size(depthDay,1), size(depthDay,2)))
allocate (TQ10(int(bottom+1)))
allocate (TQ10m(int(bottom+1)))
allocate (fTemp_stepV(size(depthDay,1), size(depthDay,2)))
allocate (fTempm_stepV(size(depthDay,1), size(depthDay,2)))
allocate (fTempV(size(depthDay,2)))
allocate (fTempmV(size(depthDay,2)))

    dist = 0._rk
    TQ10 = 0._rk
    TQ10m = 0._rk
    fTemp_stepV = 0._rk
    fTempm_stepV = 0._rk
    fTempV = 0._rk
    fTempmV  = 0._rk

!open(unit=1,action='read', file=file_path_V,status="old")!C:/Users/Admin/Desktop/FEISTY-main/FEISTY-main
do i = 1,5501
    read(1,*) tempdata(i,1),tempdata(i,2),tempdata(i,3),tempdata(i,4) !depth 0-5500(no use), tropical, temperate, boreal, default(10 celcius)
end do
close(1)
tempdata(:,5)=10._rk !default temp, so no temp-effects

dist = (depthDay + depthNight)/2._rk
! region+1: 1+1 tropical, 2+1 temperate, 3+1 boreal, 4+1 default(10 celcius)
TQ10 =  Q10**((tempdata(1:int(bottom)+1 , (region+1))-10._rk)/10._rk)
TQ10m =  Q10m**((tempdata(1:int(bottom)+1 , (region+1))-10._rk)/10._rk)


do i=1,size(dist,2)
fTemp_stepV(:,i) = dist(:,i) * TQ10
fTempm_stepV(:,i) = dist(:,i) * TQ10m
end do

fTempV = sum(fTemp_stepV,1)
fTempmV = sum(fTempm_stepV,1)

end subroutine updateTempV



end module setup

!============================================================
