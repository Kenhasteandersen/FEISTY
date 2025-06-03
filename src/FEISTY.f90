!
! FABM-FEISTY
! Yixin Zhao 2025
!
!--------------------------------------

   subroutine checknan(vec, n)
        use  setup
      integer, intent(in)::n
      real(rk), intent(inout) :: vec(n)
      integer :: i

      do i = 1, n
         if (isnan(vec(i))) vec(i) = 0._rk
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

  subroutine calcderivatives(uin, dudt)
    use  setup
      real(rk), intent(in)    :: uin(nGrid)
      real(rk), intent(inout) :: dudt(nGrid)
      real(rk):: u(nGrid)

      integer :: i, j, ii, istart, istop!, iGroup

! ----------------------------------------------------------------------
dudt=0._rk

! ----------------------------------------------
! Get proper state variable
! ----------------------------------------------
u=uin

if(bTS .eqv. .TRUE.) call allocfeisty_ts(u)

do i = 1, nGrid
  u(i) = max(0._rk , u(i))
end do

if(bET .eqv. .TRUE. .and. depthET .lt. 200) call updateET(u)

! Encounter rates
      Enc  = V*matmul(theta, u)                         ! Encounter rates     [/yr]

! ----------------------------------------------
! Mortality for resources and fish grids:
! ----------------------------------------------

      mortpred = Cmax*V/(Enc + Cmax)*u
      call checknan(mortpred, nGrid)

      mortpred = matmul(transpose(theta), mortpred)      ! Predation mortality [/yr]

! down-regulation in time-series input
if(bTS .eqv. .TRUE.)then
    dr_fac_theta = 1._rk
! original zooplankton consumption by fish
    smzcsp = mortpred(1)*u(1)
    lgzcsp = mortpred(2)*u(2)
! small zooplankton consumption cannot beyond the production
        if (mortpred(1)*u(1) > szprod) then
          dr_fac_sz = szprod / (mortpred(1)*u(1))
          do i = idxF, nGrid
            dr_fac_theta(i, 1) = dr_fac_sz
          end do
          mortpred(1) = dr_fac_sz * mortpred(1)
        end if
       !print*,(mortpred(1)*u(1))
! large zooplankton consumption cannot beyond the production
        if (mortpred(2)*u(2) > lzprod) then
         dr_fac_lz = lzprod / (mortpred(2)*u(2))
          do i = idxF, nGrid
            dr_fac_theta(i, 2) = dr_fac_lz
          end do
         mortpred(2) = dr_fac_lz * mortpred(2)
        end if
! down-regulated zooplankton consumption by fish
    smzcsp_dr = mortpred(1)*u(1)
    lgzcsp_dr = mortpred(2)*u(2)

! new Enc / (Cmax + original Enc)
    flvl = (V * (matmul((theta*dr_fac_theta), u))) /(Cmax + Enc) ! food limitation     [-]
else
! non-time-series input
    flvl = Enc/(Cmax + Enc)                                 ! food limitation     [-]
end if

! ----------------------------------------------
! Feeding $ losses for resources and fish grids:
! ----------------------------------------------
      !flvl = (V * (matmul((theta*dr_fac_theta), u))) /(Cmax + Enc)                           ! food limitation     [-]
      call checknan(flvl, nGrid)                        ! remove Nans

      Eavail  = epsAssim_vec*flvl*Cmax - metabolism         ! available energy    [/yr]

      grazing = Cmax * flvl*u                           ! grazing             [gWW/m2/yr]

      loss    = (1._rk-epsAssim_vec)*grazing + metabolism*u  ! Energy loss to environments  [gWW/m2/yr] Updated below.

! add basal and fishing mortality)
      mort = mortpred + mort0 + mortF                    ! Total mortality     [/yr]

! ----------------------------------------------
!  Flux out of the fish size group:
! ----------------------------------------------

      ! fish only data (fish grids)

      ii = 1
      do i = nResources+1, nGrid
        B(ii)        = u(i)                              ! fish stages          [g/m2]
        eFish(ii)    = Eavail(i)                         ! availabel energy     [/yr]
        eplus(ii)    = max(0._rk, Eavail(i))               ! net growth rate      [/yr]
        mortFish(ii) = mort(i)                           ! total mortality      [/yr]
        ii = ii+1
      end do

      grow = (1._rk - psiMature)*eplus                    ! energy  for growth   [/yr]

      gamma_vec = (grow - mortFish) /   &                    ! growth to next stage [/yr]
            (1._rk - (1/z)**(1._rk-mortFish/grow) )

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

      ! Add the waste energy in reproduction of each stages.
      ! Note it does not include the waste energy from last stage energy flux out, added in `totLoss` below.
        loss(ixStart(i):ixEnd(i)) = loss(ixStart(i):ixEnd(i)) + (1._rk-epsRepro_vec(i)) * Repro(istart:istop)

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

        totGrazing(i)  = 0._rk
        totLoss(i)     = 0._rk
        totMort(i)     = 0._rk

        do j = istart, istop
          totGrazing(i)  = totGrazing(i) + grazing(j)
          totLoss(i)     = totLoss(i)    + loss(j) ! updated below
          totMort(i)     = totMort(i)    + mort(j)*u(j)
        end do
      ! Add the waste energy in reproduction from flux out of the last stage of each functional group.
        totLoss(i)=totLoss(i) + (1._rk - epsRepro_vec(i)) * Fout(istop-nResources)

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

      if(bTS .eqv. .TRUE.)then
        dRdt = 0._rk
        dRdt(3) = rr(3)*(1-R(3)/K(3)) - mortRes(3)*R(3)   ! logistic formulation
      else
        if (Rtype == 1) then
          dRdt = rr*(K-R) - mortRes*R       ! chemostat formulation
        else
          dRdt = rr*R*(1-R/K) - mortRes*R   ! logistic formulation
        end if
      end if

      do i = 1, nResources
        dudt(i) = dRdt(i)
      enddo

      do i = 1, nFGrid
        dudt(i+nResources) = dBdt(i)
      enddo

  end subroutine calcderivatives

!============================================================
