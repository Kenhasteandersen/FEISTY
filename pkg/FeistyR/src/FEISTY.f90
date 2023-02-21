! ==============================================================================
! the FEISTY model in fortran, to be linked to R
! References: Petrik et al., 2019; van Denderen et al., 2020.
! The library follows MATLAB/R codes from Ken Haste Andersen; 
!        Pieter Daniel van Denderen; R?my Den?ch?re; Daniel Ottmann Riera ...
! Programmed originally by Yixin Zhao, August 2022.
! 
! Modified for linking with R by Karline Soetaert, Januari 2023
! ==============================================================================

Module modulefeisty

   implicit none
   integer, parameter :: dp = kind(0.d0) ! double precision

! .....................................
! values and vectors passed from R:
! .....................................

! dimensions and indices
   integer:: nGroups       ! Number of fish spectrum groups
   integer:: nResources    ! Number of resource state variables
   integer:: idxF          ! First index into fish groups (=nResources+1)
   integer:: nFGrid        ! Total number of fish size grids
   integer:: nGrid         ! Total number of size grids in Resources and fish

   integer, allocatable :: ixStart(:) ! First index into u for each fish spectrum
   integer, allocatable :: ixEnd(:)   ! Last index into u for 
   
   integer :: Rtype   ! resource dynamics, 1=chemostat; 2=logistic

! interaction matrix
   real(dp), allocatable:: theta(:, :)   ! Feeding preference matrix

! vectors defined in all stages
   real(dp), allocatable:: epsAssim(:)   ! assimilation efficiency 
   real(dp), allocatable:: Cmax(:)       ! maximum consumption rate 
   real(dp), allocatable:: metabolism(:) ! basal respiration rate 
   real(dp), allocatable:: mort0(:)      ! basal mortality rate 
   real(dp), allocatable:: mortF(:)      ! fishing mortality rate

! vectors defined only in fish stages
   real(dp), allocatable:: z(:)          ! ratio of mass in fish size class
   real(dp), allocatable:: psiMature(:)  ! fraction mature in fish size class

! vectors defined in fish groups
   real(dp), allocatable:: epsRepro(:)   ! efficiency of reproduction 
   
! resource parameter vectors 
   real(dp), allocatable:: K(:)          ! Carrying capacity of resources
   real(dp), allocatable:: rr(:)         ! nudging term of resources

! .....................................
! local values:
! .....................................

! defined in all stages
   real(dp), allocatable:: V(:)          ! clearance rate 
   real(dp), allocatable:: Enc(:)        ! encounter rate 
   real(dp), allocatable:: flvl(:)       ! feeding level 
   real(dp), allocatable:: grazing(:)    ! grazing rate
   real(dp), allocatable:: loss(:)       ! summed loss rates
   real(dp), allocatable:: mort(:)       ! summed mortaliy rates
   real(dp), allocatable:: mortpred(:)   ! predation mortality
   real(dp), allocatable:: eAvail(:)     ! available energy 

! Defined in fish grids only   
   real(dp), allocatable:: B(:), dBdt(:) ! state variables and derivatives
   real(dp), allocatable:: eplus(:), eFish(:) ! available energy
   real(dp), allocatable:: grow(:)       ! energy for growth [/yr]
   real(dp), allocatable:: gamma(:)      ! growth to next stage 
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
      ixStart(1) = 1
      do i = 2, nGroups
        ixstart(i) = ixstart(i-1) + ipar(ii)
        ixEnd(i-1) = ixstart(i)-1
        ii = ii+1
      end do
      ixEnd(nGroups) = ixstart(nGroups)+ipar(ii)-1
      ii = ii+1
      Rtype = ipar(ii)
      ii = ii+1
      
      nFGrid = ixEnd(nGroups)
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

      if (allocated (epsRepro)) deallocate (epsRepro)
      allocate (epsRepro(nGroups))
      
      call getVec(epsRepro,  rpar, ir, nGroups)
      
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
      
      if (allocated (epsAssim))    deallocate (epsAssim)
      allocate (epsAssim(nGrid))
      
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

      call getVec(epsAssim,    rpar, ir, nGrid)
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
      
      if (allocated (gamma))  deallocate (gamma)
      allocate (gamma(nFGrid))
      
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

      Eavail  = epsAssim*flvl*Cmax - metabolism         ! available energy    [/yr]
      
      grazing = Cmax * flvl*u                           ! grazing             [gWW/m2/yr]
      
      loss    = (1.d0-epsAssim)*grazing + metabolism*u  ! total loss          [gWW/m2/yr]   

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
    
      gamma = (grow - mortFish) /   &                    ! growth to next stage [/yr]
            (1d0 - (1/z)**(1d0-mortFish/grow) )
      
      call checknan(gamma, nFGrid)    ! No growth of fully mature classes (grow=0)

      Fout = gamma*B                                     ! flux out of stage    [g/m2/yr]

      Repro = psiMature*eplus*B                          ! reproduction         [g/m2/yr]
 
! ----------------------------------------------
! Flux into the size group
! ----------------------------------------------

      do i = 1, nGroups
      
      ! stages of this group in fish grid
        istart = ixStart(i)          
        istop  = ixEnd(i)            ! last stage

        totRepro(i)    = repro(istart)
        totBiomass(i)  = B(istart)
        
        do j = istart+1, istop
          Fin(j)         = Fout(j-1)
          totRepro(i)    = totRepro(i) + Repro(j)
          totBiomass(i)  = totBiomass(i) + B(j)
        end do
        
        totRepro(i) = totRepro(i) + Fout(istop) ! growth out ouf final stage is reproduction
        Fin(istart) = epsRepro(i)*totRepro(i)   ! reproduction

      ! stages of this group in total grid
        istart = ixStart(i) + nResources         
        istop  = ixEnd(i)  + nResources          

        totGrazing(i)  = 0d0
        totLoss(i)     = 0d0
        totMort(i)     = 0d0

        do j = istart, istop
          totGrazing(i)  = totGrazing(i) + grazing(j)
          totLoss(i)     = totLoss(i)    + loss(j)
          totMort(i)     = totMort(i)    + mort(j)*u(j)
        end do
        
      end do
      totRecruit   = totRepro*epsrepro
      
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

end module modulefeisty


!==========================================================================
!==========================================================================
! interface subroutines between R and fortran for thefeisty model
!==========================================================================
!==========================================================================

  subroutine initfeisty (steadyparms)
    use modulefeisty
    implicit none
    external steadyparms  ! not used

       feistyinitialised = .FALSE.
           
   end subroutine initfeisty

!==========================================================================
!==========================================================================
! subroutine calculating the rate of change of
! the feisty model
!==========================================================================
!==========================================================================
      
   subroutine runfeisty (neq, t, Conc, dConc, yout, ip)
    use modulefeisty
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
    use modulefeisty
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
