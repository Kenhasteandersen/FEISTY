#include "fabm_driver.h"

module FEISTY_FABM

   use fabm_types                         ! FABM
   use fabm_particle                      ! FABM
   use fabm_builtin_depth_mapping         ! FABM
   use fabm_expressions                   ! FABM
   use feisty_fabm_vertical_distribution  ! FEISTY-controlled vertical distributions

   use globals                            ! FEISTY
   use spectrum                           ! FEISTY
   use fish                               ! FEISTY
   use setup                              ! FEISTY
   
   implicit none

   private

   type, extends(type_depth_integrated_particle), public :: type_feisty_fabm
      ! State variables for fish (small pelagics, mesopelagics, large pelagics, mid-water predators, and demersal fish):
      type (type_bottom_state_variable_id),         allocatable :: id_fish(:)
      ! The vertical distributios of the fish:
      type (type_vertical_distribution_id),         allocatable :: id_fish_w(:)!, id_smpel_w(:), id_mesopel_w(:), id_lgpel_w(:), id_midp_w(:), id_dem_w(:)
      type (type_vertical_distribution_id)                      :: dummy_w
      ! The benthos state variable:
      type (type_bottom_state_variable_id)                      :: id_benthos
      type (type_bottom_state_variable_id),         allocatable :: id_ben_c(:), id_ben_n(:), id_ben_p(:)
      !type (type_bottom_state_variable_id)                      :: id_large_benthos
      type (type_state_variable_id)                             :: id_det_c, id_det_n, id_det_p
      type (type_state_variable_id)                             :: id_nut_c, id_nut_n, id_nut_p

      ! Dependency IDs for external FABM variables (e.g., small and large zooplankton, nutrients, detritus):
      type (type_bottom_dependency_id),    allocatable          :: id_smzoo_c(:), id_smzoo_n(:), id_smzoo_p(:) !NIAZ: these should be groups of dependencies in order to accept models with more zooplankton groups
      type (type_bottom_dependency_id),    allocatable          :: id_lgzoo_c(:), id_lgzoo_n(:), id_lgzoo_p(:) !NIAZ: these should be groups of dependencies in order to accept models with more zooplankton groups
      type (type_bottom_dependency_id)                          :: id_zoo_c, id_zoo_n, id_zoo_p
      ! Excretion of nitrogen (n) and phosphorus (p), allocatable to the number of fish (nGrid-nResources):
      type (type_bottom_state_variable_id),    allocatable      :: id_excre_fish_n(:), id_excre_fish_p(:)
      ! Respiration (CO2 / DIC), allocatable to the number of fish (nGrid-nResources):
      type (type_bottom_state_variable_id),    allocatable      :: id_respiration_fish_c(:),id_respiration_fish_n(:),id_respiration_fish_p(:)
      ! Fecal pellets (for detritus / POM), allocatable to the number of fish (nGrid-nResources):
      type (type_bottom_state_variable_id),    allocatable      :: id_feces_fish_c(:), id_feces_fish_n(:), id_feces_fish_p(:)
      ! Carcasses, allocatable to the number of fish (nGrid-nResources):
      type (type_bottom_state_variable_id),    allocatable      :: id_carcasses_fish_c(:), id_carcasses_fish_n(:), id_carcasses_fish_p(:)
      ! Coupling; pointer to which model contains the zooplankton state variable that we will integrater over the water column:
      type (type_model_id),    allocatable                      :: id_smzoo_int(:), id_benprey(:) ! NIAZ: these should also be allocatable to the number of zooplankton models in each instance
      type (type_model_id),    allocatable                      :: id_lgzoo_int(:) ! NIAZ: these should also be allocatable to the number of zooplankton models in each instance
      type (type_model_id)                                      :: id_zoo_int
      
      type (type_dependency_id)                                 :: id_temp, id_depth
      type (type_horizontal_dependency_id)                      :: id_temp_vertmean_100m,id_temp_vertmean_500to1500m,id_bottom_depth

      type (type_bottom_diagnostic_variable_id)                 :: id_temp_vertmean_100m_diag
      type (type_bottom_diagnostic_variable_id), allocatable    :: id_fish_total_biomass(:)
      
      real (rk)                              :: b_ini
      real (rk) ,    allocatable             :: ben_availability(:), depth_array(:) ! Array of depth values for pre-computed profiles
      real (rk) ,    allocatable             :: theta_matrix(:,:,:) ! Pre-computed theta matrix (predator,prey,n_depth)
      real (rk) ,    allocatable             :: vertical_distribution_matrix_centers(:,:,:) ! Pre-computed vertical distribution matrix (depth_center,nFGrid,n_depth) - contains center values
      
      real(rk)                               :: shelfdepth
      real(rk)                               :: max_depth
      real(rk)                               :: photic_depth
      integer                                :: n_depth, n_zoo, n_s, n_l, n_ben
      character,    allocatable              :: zoo_size(:)
      character                              :: respiration_elements, feces_elements, carcass_elements
      real(rk)                               :: bgmort    ! Fish background mortality (per year)
      real(rk)                               :: gww_gc    ! gram wet weight to gram carbon conversion factor
      real(rk)                               :: qnc, qpc  ! Nitrogen:Carbon and Phosphorus:Carbon ratios
      !real(rk)                               :: mol_targetunit  ! mol to target unit conversion factor NIAZ: always the same when we couple through standard_variables, so we use a hardcoded conversion
      real(rk)                               :: gww_mmolc, gww_FISHn, gww_FISHp   ! gww to target C/N/P unit conversion factors

      logical                                :: internal_benthos, nutrient_in_c, nutrient_in_n, nutrient_in_p           ! BGC has C/N/P elements for nutrient coupling
      logical                                :: detritus_in_c, detritus_in_n, detritus_in_p           ! BGC has C/N/P elements for detritus coupling
      logical                                :: small_zooplankton_in_c, small_zooplankton_in_n, small_zooplankton_in_p  ! Small zooplankton has C/N/P elements
      logical                                :: large_zooplankton_in_c, large_zooplankton_in_n, large_zooplankton_in_p  ! Large zooplankton has C/N/P elements
      
   contains
      procedure :: initialize
      procedure :: do_bottom
      ! Reference model procedures here.
   end type type_feisty_fabm
   
   real(rk), parameter :: seconds_per_year = 365._rk*86400._rk
   real(rk), parameter :: d_per_s      = 1.0_rk/86400.0_rk
   real(rk), parameter :: trcmin       = 5e-12_rk            ! min tracer concentration NIAZ: Remove after new benthos implementation

contains

   subroutine initialize(self, configunit)
      class (type_feisty_fabm), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
      
      real(rk)           :: smz_ini, lgz_ini, smbent_ini, lgbent_ini
      real(rk)           :: szprod, lzprod, bprodin, dfbot, depth, Tp, Tb
      real(rk)           :: dfpho, Tm, etamature, visual, Fmax, etaF!, ssigma, tau, shelfdepth
      integer            :: nStages, bET_val
      integer            :: i,j,n_l,n_s
      character(len=100) :: i_str, j_str, size_number_str, fft_long_name,fft_short_name, ns_str, nl_str
      
      class (type_vertical_depth_range), pointer :: depth_distribution
      class (type_feisty_vertical_distribution), pointer :: feisty_vertical_distribution

      !call self%type_depth_integrated_particle%initialize(configunit)
      
      ! Register model parameters  
      ! Fish physiological parameters. See input.nml in the original FEISTY R package and Fortran library
      call self%get_parameter(h, 'h', 'g^nn yr-1', 'Max. consumption coefficient', default=20._rk)  
      call self%get_parameter(nn, 'nn', '-', 'Metabolic exponent', default=-0.25_rk)  
      call self%get_parameter(gamma, 'gamma', 'm2 g^q yr-1', 'Coef. for clearance rate', default=70._rk)  
      call self%get_parameter(q, 'q', '-', 'Clearance rate exponent', default=-0.2_rk)  
      call self%get_parameter(kk, 'kk', 'g^p yr-1', 'Metabolism coefficient', default=4._rk)  
      call self%get_parameter(p, 'p', '-', 'Metabolism exponent', default=-0.175_rk)  
      call self%get_parameter(epsAssim, 'epsAssim', '-', 'Assimilation efficiency', default=0.7_rk)
      call self%get_parameter(epsRepro, 'epsRepro', '-', 'Reproduction & recruitment efficiency', default=0.01_rk)  
      call self%get_parameter(self%bgmort, 'bgmort', 'yr-1', 'Fish background mortality', default=0.1_rk) 
      
      call self%get_parameter(beta, 'beta', '-', 'Beta parameter for size-based predation preference', default=400._rk)  
      call self%get_parameter(sigma, 'sigma', '-', 'Sigma parameter for size-based predation preference', default=1.3_rk)  
      call self%get_parameter(mMin, 'mMin', 'g', 'Minimum fish mass (boundary of the grid)', default=0.001_rk)  
      call self%get_parameter(mMedium, 'mMedium', 'g', 'Medium fish central mass for feeding preference calc', default=0.5_rk)   ! for vertical distribution 10
      call self%get_parameter(mLarge, 'mLarge', 'g', 'Large fish central mass for feeding preference calc', default=250._rk)     ! for vertical distribution 5000
      
          
      !call self%get_parameter(lbenk, 'lbenk', 'g m-2', 'Large benthos carry capacity', default=0._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(szoog, 'szoog', 'yr-1', 'Small zooplankton growth rate', default=1._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(lzoog, 'lzoog', 'yr-1', 'Large zooplankton growth rate', default=1._rk)  !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(sbeng, 'sbeng', 'yr-1', 'Small benthos growth rate', default=1._rk)      !NIAZ: deleted (no zoop parameters should be in the fish code)
      !call self%get_parameter(lbeng, 'lbeng', 'yr-1', 'Large benthos growth rate', default=0._rk)      !NIAZ: deleted (no zoop parameters should be in the fish code)
      ! predation preference coefficient
      call self%get_parameter(thetaS, 'thetaS', '-', 'Medium fish preference for small zooplankton', default=0.25_rk)  
      call self%get_parameter(thetaA, 'thetaA', '-', 'Large fish preference for medium forage fish', default=0.5_rk)  
      call self%get_parameter(thetaD, 'thetaD', '-', 'Preference of large demersal on pelagic prey', default=0.75_rk)  
      ! get initial values
      call self%get_parameter(smz_ini,'smz_ini', 'g m-2',  'initial small mesozooplankton biomass', default=1.e2_rk)  
      call self%get_parameter(lgz_ini, 'lgz_ini', 'g m-2', 'initial large mesozooplankton biomass', default=1.e2_rk)  
      call self%get_parameter(smbent_ini, 'smbent_ini', 'g m-2', 'initial small benthos biomass', default=5._rk)  
      call self%get_parameter(lgbent_ini, 'lgbent_ini', 'g m-2', 'initial large benthos biomass', default=0._rk)        
      call self%get_parameter(self%b_ini, 'b_ini', 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  
      ! get input values
      call self%get_parameter(szprod, 'szprod', 'g m-2 yr-1', 'small mesozooplankton production', default=100._rk)  
      call self%get_parameter(lzprod, 'lzprod', 'g m-2 yr-1', 'large mesozooplankton production', default=100._rk)  
      call self%get_parameter(bprodin, 'bprodin', 'g m-2 yr-1', 'benthos production', default=5._rk)  
      call self%get_parameter(dfbot, 'dfbot', 'g m-2 yr-1', 'detrital flux reaching the bottom', default=-1._rk)  
      call self%get_parameter(depth, 'depth', 'm', 'water column depth', default=100._rk)  
      call self%get_parameter(Tp, 'Tp', 'Celsius', 'pelagic layer averaged temperature', default=10._rk)  
      call self%get_parameter(Tb, 'Tb', 'Celsius', 'bottom layer depth temperature (500m - up to 1500m)', default=8._rk) 
      
      call self%get_parameter(dfpho, 'dfpho', 'g m-2 yr-1', 'detrital flux out of the photic zone', default=-1._rk)
      call self%get_parameter(nStages, 'nStages', '-', 'size number of large fish functional types', default=9)
      call self%get_parameter(Tm, 'Tm', 'Celsius', 'mid-water temperature', default=Tb)
      call self%get_parameter(self%photic_depth, 'photic_depth', 'm', 'photic zone depth', default=100._rk)
      call self%get_parameter(etamature, 'etamature', '-', 'the coefficient determines the fish size with a 50% maturity level', default=0.25_rk)
      call self%get_parameter(self%shelfdepth, 'shelfdepth', 'm', 'continental shelf depth', default=250._rk)
      call self%get_parameter(visual, 'visual', '-', 'the coefficient determines the visual ability of fish', default=1.5_rk) ! 1.5:visual predator or 1:non-visual predator. Be careful to use other values.
      call self%get_parameter(Fmax, 'Fmax', 'yr-1', 'Maximum fishing mortality', default=0._rk)
      call self%get_parameter(etaF, 'etaF', '-', 'the coefficient determines the fish size with 50% fishing selectivity', default=0.05_rk)
      call self%get_parameter(ssigma, 'ssigma', '-', 'width of initial vertical distribution', default=10._rk)
      call self%get_parameter(tau, 'tau', '-', 'increase in width', default=10._rk)
      
      call self%get_parameter(bET_val, 'bET_val', '-', 'whether turn on the effective temperature (integer 1 or 0)', default=1)
      
      call self%get_parameter(self%max_depth, 'max_depth', 'm', 'maximum depth of the model domain', default=5000._rk)
      call self%get_parameter(self%n_depth, 'n_depth', '-', 'number of depth levels for vertical distribution matrices', default=50)
      call self%get_parameter(self%gww_gc, 'gww_gc', '-', 'gram wet weight to gram carbon ratio', default=9._rk)
      call self%get_parameter(self%qnc, 'qnc', 'N/C', 'Nitrogen to Carbon ratio', default=16.0_rk/106.0_rk)
      call self%get_parameter(self%qpc, 'qpc', 'P/C', 'Phosphorus to Carbon ratio', default=1.0_rk/106.0_rk)
      !call self%get_parameter(self%mol_targetunit, 'mol_targetunit', '-', 'mol to target unit (in LTL model) ratio', default=1._rk/1000._rk)! default is mol m-3:mmol m-3 (used in LTL model)
      
      call self%get_parameter(self%feces_elements, 'feces_elements', '', 'elemental composition of feces e.g. "cnp"')
      call self%get_parameter(self%carcass_elements, 'carcass_elements', '', 'elemental composition of carcass e.g. "cnp"')
      call self%get_parameter(self%respiration_elements, 'respiration_elements', '', 'elemental composition of respiration e.g. "cnp"')
      
      ! Element switches indicating which BGC elements are available for coupling (C, N, P) !NIAZ: These should probably go.
      call self%get_parameter(self%nutrient_in_c, 'nutrient_in_c', '-', 'BGC has carbon (C) element for coupling', default=.false.)
      call self%get_parameter(self%nutrient_in_n, 'nutrient_in_n', '-', 'BGC has nitrogen (N) element for coupling', default=.true.)
      call self%get_parameter(self%nutrient_in_p, 'nutrient_in_p', '-', 'BGC has phosphorus (P) element for coupling', default=.true.)
      
      ! Element switches indicating which BGC elements are available for detritus coupling (C, N, P) !NIAZ: These should probably go.
      call self%get_parameter(self%detritus_in_c, 'detritus_in_c', '-', 'BGC has carbon (C) element for detritus coupling', default=.false.)
      call self%get_parameter(self%detritus_in_n, 'detritus_in_n', '-', 'BGC has nitrogen (N) element for detritus coupling', default=.true.)
      call self%get_parameter(self%detritus_in_p, 'detritus_in_p', '-', 'BGC has phosphorus (P) element for detritus coupling', default=.true.)
      
      ! Element switches indicating which elements are available in small zooplankton (C, N, P) !NIAZ: These should probably go.
      call self%get_parameter(self%small_zooplankton_in_c, 'small_zooplankton_in_c', '-', 'Small zooplankton has carbon (C) element', default=.false.)
      call self%get_parameter(self%small_zooplankton_in_n, 'small_zooplankton_in_n', '-', 'Small zooplankton has nitrogen (N) element', default=.true.)
      call self%get_parameter(self%small_zooplankton_in_p, 'small_zooplankton_in_p', '-', 'Small zooplankton has phosphorus (P) element', default=.true.)
      
      ! Element switches indicating which elements are available in large zooplankton (C, N, P) !NIAZ: These should probably go.
      call self%get_parameter(self%large_zooplankton_in_c, 'large_zooplankton_in_c', '-', 'Large zooplankton has carbon (C) element', default=.false.)
      call self%get_parameter(self%large_zooplankton_in_n, 'large_zooplankton_in_n', '-', 'Large zooplankton has nitrogen (N) element', default=.true.)
      call self%get_parameter(self%large_zooplankton_in_p, 'large_zooplankton_in_p', '-', 'Large zooplankton has phosphorus (P) element', default=.true.)
      
      ! Calculate gww_targetc: gww to target C unit conversion factor NIAZ: all this has to go
      self%gww_mmolc = self%gww_gc * 12.01_rk * 1._rk/1000._rk ! Redfield unit conversion from mmol C m-3 to gWW m-3
      
      ! Calculate gww_targetn and gww_targetp: gww to target N and P unit conversion factors
      self%gww_FISHn = self%gww_mmolc / self%qnc                   ! gww to target N unit
      self%gww_FISHp = self%gww_mmolc / self%qpc                   ! gww to target P unit
      
      ! Register model parameters and variables here.
      
      ! convert to per second
      h     = h/seconds_per_year !NIAZ change to scale factor at parameter registration - also, not used?
      gamma = gamma/seconds_per_year !NIAZ change to scale factor at parameter registration - also, not used?
      kk    = kk/seconds_per_year !NIAZ change to scale factor at parameter registration - also, not used?
      
      !select case (setupidx)
      !   case ()
      !call setupbasic(szprod, lzprod, bprodin, dfbot, depth, Tp, Tb)
            !call setupbasic2(szprod, lzprod, bprodin, dfbot,nStages, depth, Tp, Tb,etaMature,Fmax,etaF,bET_val)
      !   case()
       !call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,3,Tp,Tm,Tb,depth,self%photic_depth,etamature,shelfdepth,visual,Fmax,etaF)
       call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,nStages,Tp,Tm,Tb,depth,self%photic_depth,etamature,self%shelfdepth,visual,Fmax,etaF)
      !      
      
      !end select     
      
      ! Depth-integrated dependencies and coupling setup
      call self%register_vertical_distribution(self%dummy_w,'dummy')
      allocate(depth_distribution)
      call self%add_child(depth_distribution, 'dummy')
      call self%request_coupling(self%dummy_w, 'dummy'//'/'//'w')
       
      ! Get the number of zooplankton types. NIAZ: vectorized for more zooplankton groups
      call self%get_parameter(self%n_zoo, 'n_zoo', '', 'number of zooplankton prey')
      allocate(self%zoo_size(self%n_zoo))
      do i = 1, self%n_zoo ! find the number of large and small zooplankton models, register dependencies to both types, 
          write (i_str,'(i0)') i
          self%n_l = 0
          self%n_s = 0
          if (self%n_zoo > 1) then !if more than one zoopankton
              call self%get_parameter(self%zoo_size(i), 'zoo_size'//trim(i_str), '', 'small or large zooplankton? "s" for small, "l" for large') ! register whether current model is large or small
              if (self%zoo_size(i) == 's') then ! if small
                  self%n_s = self%n_s + 1 ! iterate counter
                  write (ns_str,'(i0)') self%n_s 
                  call self%register_model_dependency('small_zooplankton_target'//TRIM(ns_str)) ! register model dependency
              else if (self%zoo_size(i) == 'l') then ! if large
                  self%n_l = self%n_l + 1  ! iterate counter
                  write (nl_str,'(i0)') self%n_l 
                  call self%register_model_dependency('large_zooplankton_target'//TRIM(i_str)) ! register model dependency
              end if
          else if (self%n_zoo == 1) then ! if only one zooplankton
              call self%register_model_dependency('zooplankton_target')
              self%n_l = 1
              self%n_s = 1
          end if   
      end do
      
      ! iterate over small zooplankton models and set up couplings
      if (self%n_zoo > 1) then
          allocate(self%id_smzoo_c(self%n_s))
          allocate(self%id_smzoo_n(self%n_s))
          allocate(self%id_smzoo_p(self%n_s))
          allocate(self%id_smzoo_int(self%n_s))
          do i = 1, self%n_s ! for each small zooplankton type
              write (i_str,'(i0)') i
              call self%register_dependency(self%id_smzoo_c(i), 'small_zoo_c'//trim(i_str), 'mmol C m-2', 'depth-integrated small zooplankton carbon'//trim(i_str))
              call self%register_dependency(self%id_smzoo_n(i), 'small_zoo_n'//trim(i_str), 'mmol N m-2', 'depth-integrated small zooplankton nitrogen'//trim(i_str))
              call self%register_dependency(self%id_smzoo_p(i), 'small_zoo_p'//trim(i_str), 'mmol P m-2', 'depth-integrated small zooplankton phosphorus'//trim(i_str))
              call self%request_mapped_coupling_to_model(self%id_smzoo_c(i), 'small_zooplankton_c'//trim(i_str), standard_variables%total_carbon, id_w=self%dummy_w)
              call self%request_mapped_coupling_to_model(self%id_smzoo_n(i), 'small_zooplankton_n'//trim(i_str), standard_variables%total_nitrogen, id_w=self%dummy_w)
              call self%request_mapped_coupling_to_model(self%id_smzoo_p(i), 'small_zooplankton_p'//trim(i_str), standard_variables%total_phosphorus, id_w=self%dummy_w)
              call self%couplings%set_string('small_zooplankton_c'//trim(i_str), "small_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              call self%couplings%set_string('small_zooplankton_n'//trim(i_str), "small_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              call self%couplings%set_string('small_zooplankton_p'//trim(i_str), "small_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              
              call self%register_mapped_model_dependency(self%id_smzoo_int(i), 'small_zooplankton_target'//trim(i_str), proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w) ! NIAZ: is this acutally neccesary?
          enddo      
          
          ! iterate over large zooplankton models and set up couplings
          allocate(self%id_lgzoo_c(self%n_l))
          allocate(self%id_lgzoo_n(self%n_l))
          allocate(self%id_lgzoo_p(self%n_l))
          allocate(self%id_lgzoo_int(self%n_l))
          do i = 1, self%n_l ! for each small zooplankton type
              write (i_str,'(i0)') i
              call self%register_dependency(self%id_lgzoo_c(i), 'large_zoo_c'//trim(i_str), 'mmol C m-2', 'depth-integrated large zooplankton carbon'//trim(i_str))
              call self%register_dependency(self%id_lgzoo_n(i), 'large_zoo_n'//trim(i_str), 'mmol N m-2', 'depth-integrated large zooplankton nitrogen'//trim(i_str))
              call self%register_dependency(self%id_lgzoo_p(i), 'large_zoo_p'//trim(i_str), 'mmol P m-2', 'depth-integrated large zooplankton phosphorus'//trim(i_str))
              call self%request_mapped_coupling_to_model(self%id_lgzoo_c(i), 'large_zooplankton_c'//trim(i_str), standard_variables%total_carbon, id_w=self%dummy_w)
              call self%request_mapped_coupling_to_model(self%id_lgzoo_n(i), 'large_zooplankton_n'//trim(i_str), standard_variables%total_nitrogen, id_w=self%dummy_w)
              call self%request_mapped_coupling_to_model(self%id_lgzoo_p(i), 'large_zooplankton_p'//trim(i_str), standard_variables%total_phosphorus, id_w=self%dummy_w)
              call self%couplings%set_string('large_zooplankton_c'//trim(i_str), "large_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              call self%couplings%set_string('large_zooplankton_n'//trim(i_str), "large_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              call self%couplings%set_string('large_zooplankton_p'//trim(i_str), "large_zooplankton_target"//trim(i_str)) ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
              
              call self%register_mapped_model_dependency(self%id_lgzoo_int(i), 'large_zooplankton_target'//trim(i_str), proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w)
          enddo
      
          ! if only one zooplankton, set up couplings to that one.
      else
          call self%register_dependency(self%id_zoo_c, 'zoo_c', 'mmol C m-2', 'depth-integrated zooplankton carbon')
          call self%register_dependency(self%id_zoo_n, 'zoo_n', 'mmol N m-2', 'depth-integrated zooplankton nitrogen')
          call self%register_dependency(self%id_zoo_p, 'zoo_p', 'mmol P m-2', 'depth-integrated zooplankton phosphorus')
          call self%request_mapped_coupling_to_model(self%id_zoo_c, 'zooplankton_c', standard_variables%total_carbon, id_w=self%dummy_w)
          call self%request_mapped_coupling_to_model(self%id_zoo_n, 'zooplankton_n', standard_variables%total_nitrogen, id_w=self%dummy_w)
          call self%request_mapped_coupling_to_model(self%id_zoo_p, 'zooplankton_p', standard_variables%total_phosphorus, id_w=self%dummy_w)
          call self%couplings%set_string('zooplankton_c', "zooplankton_target") ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
          call self%couplings%set_string('zooplankton_n', "zooplankton_target") ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
          call self%couplings%set_string('zooplankton_p', "zooplankton_target") ! links every constituent dependency to the same prey model pointer, so you only need to point to it once in the yaml
          
          call self%register_mapped_model_dependency(self%id_zoo_int, 'zooplankton_target', proportional_change=.true., domain=domain_bottom,id_w=self%dummy_w)
      endif
      
      ! Determine whether benthos will be modelled internally or externally
      call self%get_parameter(self%internal_benthos, 'internal_benthos', '', 'benthos modelled internally (True) or from external zoobenthos model (False)')
      if (self%internal_benthos) then ! if no coupling to external zoobenthos prey, setup internal zoobenthos variables, and couple to detrital fluxes instead.
          !Register fake deposition pools, for POM/detritus to dump their sedimentation in, and register a dependency on its ROC, to extract the flux
          call self%register_dependency(self%id_det_flux_c,  'det_flux_c',   'mmol m-3 s-1',  'bottom flux '//trim(index))
          
          
          call self%register_diagnostic_variable(self%id_pool(i), 'pool'//trim(index),'kmol m-2', 'target bottom pool for pelagic tracer '//trim(index),source=source_constant, output=output_none, act_as_state_variable=.true.)
          
          call self%register_diagnostic_variable(self%id_flux_out(i), 'pool'//trim(index)//'flux','kmol m-2', 'flux to target bottom pool for pelagic tracer '//trim(index),source=source_do, output=output_instantaneous)
          
          call self%request_coupling(self%id_flux(i),'./pool'//trim(index)//'_sms_tot')
          
               
               
          call self%register_state_variable(self%id_benthos, 'benthos', 'g m-2', 'benthos biomass', initial_value=smbent_ini, minimum=0.0_rk)
          if (self%detritus_in_c) call self%register_state_dependency(self%id_det_c, 'detritus_c', 'mmol C m-3', 'detritus carbon reaching the bottom for driving benthos')
          if (self%detritus_in_n) call self%register_state_dependency(self%id_det_n, 'detritus_n', 'mmol N m-3', 'detritus nitrogen reaching the bottom for driving benthos')
          if (self%detritus_in_p) call self%register_state_dependency(self%id_det_p, 'detritus_p', 'mmol P m-3', 'detritus phosphorus reaching the bottom for driving benthos')
          if (self%nutrient_in_c) call self%register_state_dependency(self%id_nut_c, 'nutrient_c', 'mmol C m-3', 'nutrient carbon in the bottom layer')
          if (self%nutrient_in_n) call self%register_state_dependency(self%id_nut_n, 'nutrient_n', 'mmol N m-3', 'nutrient nitrogen in the bottom layer')
          if (self%nutrient_in_p) call self%register_state_dependency(self%id_nut_p, 'nutrient_p', 'mmol P m-3', 'nutrient phosphorus in the bottom layer')
          
          ! Add benthos to elemental aggregates, conditional on which BGC elements are used
          call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_benthos, scale_factor = 1._rk/self%gww_mmolc)     !if (self%nutrient_in_c) 
          call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_benthos, scale_factor = 1._rk/self%gww_FISHn)   !if (self%nutrient_in_n) 
          call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_benthos, scale_factor = 1._rk/self%gww_FISHp) !if (self%nutrient_in_p)
          
      else ! if external zoobenthos models are available, set up couplings to zoobenthos prey.
          call self%get_parameter(self%n_ben, 'n_ben', '', 'number of zoobenthos prey')
          
          allocate(self%ben_availability(self%n_ben))
          allocate(self%id_benprey(self%n_ben))
          allocate(self%id_ben_c(self%n_ben))
          allocate(self%id_ben_n(self%n_ben))
          allocate(self%id_ben_p(self%n_ben))
          do i = 1, self%n_ben ! find the number of zoobenthos models, register their individual availability, and set up automatic couplings to their c, n, and p constituents 
              write (i_str,'(i0)') i
              call self%get_parameter(self%ben_availability(i), 'ben_'//TRIM(i_str)//'_availability', '', 'availability of zoobenthos prey'//TRIM(i_str)//' (between 0 and 1)')
              call self%register_model_dependency(self%id_benprey(i), 'zoobenthos_target'//TRIM(i_str)) ! register model dependency
        
              ! register dependencies on zoobenthos C,N,P
              call self%register_state_dependency(self%id_ben_c(i), 'zoobenthos_c_'//trim(i_str), 'mmol m-3', 'zoobenthos '//trim(i_str)//' c biomass')
              call self%register_state_dependency(self%id_ben_n(i), 'zoobenthos_n_'//trim(i_str), 'mmol m-3', 'zoobenthos '//trim(i_str)//' n biomass')
              call self%register_state_dependency(self%id_ben_p(i), 'zoobenthos_p_'//trim(i_str), 'mmol m-3', 'zoobenthos '//trim(i_str)//' p biomass')
              
              call self%request_coupling_to_model(self%id_ben_c(i), self%id_benprey(i), standard_variables%total_carbon)
              call self%request_coupling_to_model(self%id_ben_n(i), self%id_benprey(i), standard_variables%total_nitrogen)
              call self%request_coupling_to_model(self%id_ben_p(i), self%id_benprey(i), standard_variables%total_phosphorus)
          enddo
      endif
      
      ! register model dependencies for recipient nutrient and POM/detritus pools
      if (index(self%respiration_elements, 'c') /= 0) call self%register_model_dependency("respiration_c_target")
      if (index(self%respiration_elements, 'n') /= 0) call self%register_model_dependency("respiration_n_target")
      if (index(self%respiration_elements, 'p') /= 0) call self%register_model_dependency("respiration_p_target")
      if (index(self%feces_elements, 'c') /= 0) call self%register_model_dependency("feces_c_target")
      if (index(self%feces_elements, 'n') /= 0) call self%register_model_dependency("feces_n_target")
      if (index(self%feces_elements, 'p') /= 0) call self%register_model_dependency("feces_p_target")
      if (index(self%carcass_elements, 'c') /= 0) call self%register_model_dependency("carcasses_c_target")
      if (index(self%carcass_elements, 'n') /= 0) call self%register_model_dependency("carcasses_n_target")
      if (index(self%carcass_elements, 'p') /= 0) call self%register_model_dependency("carcasses_p_target")

      allocate(self%id_fish(nGrid-nResources))
      allocate(self%id_fish_w(nGrid-nResources))! allocate size of vertical distribution of small pelagics     
      allocate(self%id_excre_fish_n(nGrid-nResources), self%id_excre_fish_p(nGrid-nResources))
      allocate(self%id_respiration_fish_c(nGrid-nResources),self%id_respiration_fish_n(nGrid-nResources),self%id_respiration_fish_p(nGrid-nResources))
      allocate(self%id_feces_fish_c(nGrid-nResources), self%id_feces_fish_n(nGrid-nResources), self%id_feces_fish_p(nGrid-nResources))
      allocate(self%id_carcasses_fish_c(nGrid-nResources), self%id_carcasses_fish_n(nGrid-nResources), self%id_carcasses_fish_p(nGrid-nResources))

      allocate(self%id_fish_total_biomass(nGroups))

      ! Prepare pre-computed matrices for depth-dependent fish vertical distribution and predation preference (MUST be done before creating child models)
      call initialize_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)  ! Read from file or compute theta/vertical distribution matrices

      !assign fish background mortality to FEISTY
      mort0(idxF:nGrid) = self%bgmort/seconds_per_year !NIAZ change to scale factor at parameter registration
      mortF(idxF:nGrid) = mortF/seconds_per_year !NIAZ change to scale factor at parameter registration

      ! registration of fish state variables and diagnostic variables
      do i = 1, nGroups
         write (i_str,'(i0)') i
         ! read functional type longnames from yaml
         call self%get_parameter(fft_long_name,  'fft_long_name'//trim(i_str),  units='', long_name='fish functional type long name',  default="fish_functional_type_"//trim(i_str)) !NIAZ: not used
         ! read functional type shortnames from yaml
         call self%get_parameter(fft_short_name, 'fft_short_name_'//trim(i_str), units='', long_name='fish functional type short name', default="fft_"//trim(i_str)) !NIAZ: not used
         ! read the initial value of each functional type from yaml
         call self%get_parameter(self%b_ini, 'b_ini_'//trim(i_str), 'g m-2', 'initial fish biomass of each size class', default=1.e-5_rk)  

         do j = ixStart(i)-nResources, ixEnd(i)-nResources
               write (j_str,'(i0)') j
               write (size_number_str,'(i0)') j+nResources-ixStart(i)+1
               ! assign the initial biomass to each fish size class of a functional type.
               !!change names!!
                call self%register_state_variable(self%id_fish(j), trim(fft_short_name)//'_size_'//trim(size_number_str), 'g m-2', trim(fft_short_name)//'_size_'//trim(size_number_str)//'_biomass', initial_value=self%b_ini, minimum=0.0_rk)
         end do
            
         call self%register_diagnostic_variable(self%id_fish_total_biomass(i), trim(fft_short_name)//'_totB', 'g m-2', trim(fft_short_name)//'_total_biomass')   
      end do 
      
      !sets up biogeochemical coupling and vertical distribution for each fish size class
      do i = 1, nGrid-nResources !NIAZ change this to a nested loop over functional types and size classes like in the previous loop, and use the functional type short name for the variable and coupling names
         write (i_str,'(i0)') i
         call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_fish(i), scale_factor = 1._rk/self%gww_mmolc)     !if (self%nutrient_in_c) NIAZ: remove the target scaling and use redfield.
         call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_fish(i), scale_factor = 1._rk/self%gww_FISHn)   !if (self%nutrient_in_n) NIAZ: remove the target scaling and use redfield.
         call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_fish(i), scale_factor = 1._rk/self%gww_FISHp) !if (self%nutrient_in_p) NIAZ: remove the target scaling and use redfield.
         
         ! Register vertical distribution for this fish size class
         call self%register_vertical_distribution(self%id_fish_w(i),'fish_'//trim(i_str))
         
         ! Create FEISTY-controlled vertical distribution model
         allocate(feisty_vertical_distribution)
         call self%add_child(feisty_vertical_distribution, 'feisty_vertical_distribution_'//trim(i_str), configunit=-1)
         call self%request_coupling(self%id_fish_w(i), 'feisty_vertical_distribution_'//trim(i_str)//'/'//'w')
         
         ! Set pointers to parent FEISTY model's data
         feisty_vertical_distribution%vertical_distribution_matrix_centers => self%vertical_distribution_matrix_centers
         feisty_vertical_distribution%depth_array => self%depth_array
         feisty_vertical_distribution%fish_index = i
                  
         ! Respiration - register state dependency and set up coupling
         if (index(self%respiration_elements, 'c') /= 0) then
             call self%register_state_dependency(self%id_respiration_fish_c(i),'respiration_c', 'mmol C m-2', 'respiration carbon')
             call self%request_mapped_coupling_to_model(self%id_respiration_fish_c(i),'respiration_fish_'//trim(i_str)//'_c',standard_variables%total_carbon, id_w=self%id_fish_w(i))
             call self%couplings%set_string('respiration_fish_'//trim(i_str)//'_c', "respiration_c_target") ! links the respiration of every fish to the same pool, so you only need to point to it once in the yaml
         endif
         
         if (index(self%respiration_elements, 'n') /= 0) then
             call self%register_state_dependency(self%id_respiration_fish_n(i),'respiration_n', 'mmol N m-2', 'respiration nitrogen')
             call self%request_mapped_coupling_to_model(self%id_respiration_fish_n(i),'respiration_fish_'//trim(i_str)//'_n',standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
             call self%couplings%set_string('respiration_fish_'//trim(i_str)//'_n', "respiration_n_target") ! links the respiration of every fish to the same pool, so you only need to point to it once in the yaml
         endif
         
         if (index(self%respiration_elements, 'p') /= 0) then
             call self%register_state_dependency(self%id_respiration_fish_p(i),'respiration_p', 'mmol P m-2', 'respiration phosphorus')
             call self%request_mapped_coupling_to_model(self%id_respiration_fish_p(i),'respiration_fish_'//trim(i_str)//'_p',standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
             call self%couplings%set_string('respiration_fish_'//trim(i_str)//'_p', "respiration_p_target") ! links the respiration of every fish to the same pool, so you only need to point to it once in the yaml
         end if

         ! Feces - register state dependency and set up coupling
         if (index(self%feces_elements, 'c') /= 0) then
            call self%register_state_dependency(self%id_feces_fish_c(i),'feces_c', 'mmol C m-2', 'feces carbon')
            call self%request_mapped_coupling_to_model(self%id_feces_fish_c(i), 'feces_fish_'//trim(i_str)//'_c',standard_variables%total_carbon, id_w=self%id_fish_w(i))
            call self%couplings%set_string('feces_fish_'//trim(i_str)//'_c', "feces_c_target")
         end if
         if (index(self%feces_elements, 'n') /= 0) then
            call self%register_state_dependency(self%id_feces_fish_n(i),'feces_n', 'mmol N m-2', 'feces nitrogen')
            call self%request_mapped_coupling_to_model(self%id_feces_fish_n(i), 'feces_fish_'//trim(i_str)//'_n',standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
            call self%couplings%set_string('feces_fish_'//trim(i_str)//'_n', "feces_n_target")
         end if
         if (index(self%feces_elements, 'p') /= 0) then
            call self%register_state_dependency(self%id_feces_fish_p(i),'feces_p', 'mmol P m-2', 'feces phosphorus')
            call self%request_mapped_coupling_to_model(self%id_feces_fish_p(i), 'feces_fish_'//trim(i_str)//'_p',standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
            call self%couplings%set_string('feces_fish_'//trim(i_str)//'_p', "feces_p_target")
         end if

         ! Carcasses - register state dependency and set up coupling
         if (index(self%carcass_elements, 'c') /= 0) then
            call self%register_state_dependency(self%id_carcasses_fish_c(i),'carcasses_c', 'mmol C m-2', 'carcasses carbon')
            call self%request_mapped_coupling_to_model(self%id_carcasses_fish_c(i), 'carcasses_fish_'//trim(i_str)//'_c',standard_variables%total_carbon, id_w=self%id_fish_w(i))
            call self%couplings%set_string('carcasses_fish_'//trim(i_str)//'_c', "carcasses_c_target")
         end if
         if (index(self%carcass_elements, 'n') /= 0) then
            call self%register_state_dependency(self%id_carcasses_fish_n(i),'carcasses_n', 'mmol N m-2', 'carcasses nitrogen')
            call self%request_mapped_coupling_to_model(self%id_carcasses_fish_n(i), 'carcasses_fish_'//trim(i_str)//'_n',standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
            call self%couplings%set_string('carcasses_fish_'//trim(i_str)//'_n', "carcasses_n_target")
         end if
         if (index(self%carcass_elements, 'p') /= 0) then
            call self%register_state_dependency(self%id_carcasses_fish_p(i),'carcasses_p', 'mmol P m-2', 'carcasses phosphorus')
            call self%request_mapped_coupling_to_model(self%id_carcasses_fish_p(i), 'carcasses_fish_'//trim(i_str)//'_p',standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
            call self%couplings%set_string('carcasses_fish_'//trim(i_str)//'_p', "carcasses_p_target")
         end if

        !call self%request_mapped_coupling_to_model(self%id_excre_fish_n(i), 'excretion_fish_'//trim(i_str),standard_variables%total_nitrogen, id_w=self%id_fish_w(i))
         !call self%request_mapped_coupling_to_model(self%id_excre_fish_p(i), 'excretion_fish_'//trim(i_str),standard_variables%total_phosphorus, id_w=self%id_fish_w(i))
        ! call self%couplings%set('excretion_fish_'//trim(i_str), "excretion") 
         
      end do
       
       !temperature and depth
       call self%register_dependency(self%id_temp, standard_variables%temperature)
       call self%register_dependency(self%id_temp_vertmean_100m, vertical_mean(self%id_temp, maximum_depth=100._rk))
       call self%register_dependency(self%id_temp_vertmean_500to1500m, vertical_mean(self%id_temp, minimum_depth=500._rk, maximum_depth=1500._rk))
       call self%register_diagnostic_variable(self%id_temp_vertmean_100m_diag, 'temp_vertmean_100m', 'degree_C', 'vertical mean temperature above 100 m')
       call self%register_dependency(self%id_bottom_depth, standard_variables%bottom_depth)
       call self%register_dependency(self%id_depth,     standard_variables%depth)
       
!   -----declared in setup.f90-----
      allocate (excretion(nGrid))
      allocate (respiration(nGrid))
      allocate (carcasses(nGrid))
      allocate (feces(nGrid))
      
   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_feisty_fabm), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: c, temp, prey_c, prey_n, prey_p, prey_s, w_int
      real(rk) :: ingestion_c, ingestion_n, ingestion_p, prey_loss_rate, p1,p2, net_growth,dflag,wd
      integer  :: iGroup, i, istate, iprey, closest_bottom_depth_idx
      real(rk),dimension(ixEnd(1)-ixStart(1)+1)  :: smpel
      real(rk),dimension(ixEnd(2)-ixStart(2)+1)  :: lgpel
      real(rk),dimension(ixEnd(3)-ixStart(3)+1)  :: dem
      real(rk),dimension(nGrid-nResources)       :: fish, rel_fish_cons, pred_morts
      !real(rk)                                   :: smzoo_c, smzoo_n, smzoo_p, lgzoo_c, lgzoo_n, lgzoo_p
      real(rk)                                   :: excess_c,excess_n,excess_p
      real(rk)                                   :: excess_c_smzoo,excess_n_smzoo, excess_p_smzoo
      real(rk)                                   :: excess_c_lgzoo,excess_n_lgzoo, excess_p_lgzoo
      real(rk)                                   :: zoo1, zoo2, benthos
      real(rk)                                   :: det_bot, det_bot_c, det_bot_n, det_bot_p
      real(rk)                                   :: det_bot_flux_gww, benthos_g_gww, benthos_loss_gww
      real(rk),dimension(nGrid)                  :: uin, dudt, mortpred_contri_zoo1, mortpred_contri_zoo2
      real(rk)                                   :: temp_100m,temp_bottom, temp_500to1500m, bottom_depth, depth
      real(rk)                                   :: photic_depth_local  ! local photic depth adjusted per grid cell
      real(rk),dimension(self%n_s)               :: smzoo_c, smzoo_n, smzoo_p, smzoo_cn, smzoo_cp, smzoo_c_rd, smzoo_n_rd, smzoo_p_rd, smzoo_ww_rf
      real(rk),dimension(self%n_l)               :: lgzoo_c, lgzoo_n, lgzoo_p, lgzoo_cn, lgzoo_cp, lgzoo_c_rd, lgzoo_n_rd, lgzoo_p_rd, lgzoo_ww_rf
      real(rk),dimension(1)                      :: zoo_c, zoo_n, zoo_p 
      real(rk),dimension(self%n_ben)             :: benthos_c, benthos_n, benthos_p, ben_ww_rf
      real(rk)                                   :: smzoo_c_tot, smzoo_n_tot, smzoo_p_tot, lgzoo_c_tot, lgzoo_n_tot, lgzoo_p_tot, smzoo_ww_tot, lgzoo_ww_tot, smzoo_mort_ww, lgzoo_mort_ww
      real(rk)                                   :: excess_c_zoo, excess_n_zoo, excess_p_zoo, ben_ww_tot, excess_c_ben, excess_n_ben, excess_p_ben, ben_mort_ww
      
      _BOTTOM_LOOP_BEGIN_

         ! Get depth-integrated predator biomass       
         do i = 1, nGrid-nResources
         _GET_BOTTOM_(self%id_fish(i), fish(i))
         end do      
         !print*,(sum(fish))
         
         ! get depth-integrated zooplankton
         if (self%n_zoo > 1) then ! if theres more than one zooplankton, get all biomasses and then sum up for each size group.
             do i = 1, self%n_s
                _GET_BOTTOM_(self%id_smzoo_c(i), smzoo_c(i)) ! NIAZ: These are in mmol m-3
                _GET_BOTTOM_(self%id_smzoo_n(i), smzoo_n(i)) ! NIAZ: These are in mmol m-3
                _GET_BOTTOM_(self%id_smzoo_p(i), smzoo_p(i)) ! NIAZ: These are in mmol m-3
             enddo
             ! calculate total small zoop. biomasss
             smzoo_c_tot = sum(smzoo_c)
             smzoo_n_tot = sum(smzoo_n)
             smzoo_p_tot = sum(smzoo_p)
                          
             ! calculate WW of limiting nutrient for each small zooplankton
             smzoo_ww_rf = min(merge(smzoo_c*self%gww_mmolc, huge(0._rk), smzoo_c /= 0._rk), &
                               merge(smzoo_n*self%gww_FISHn, huge(0._rk), smzoo_n /= 0._rk), &
                               merge(smzoo_p*self%gww_FISHp, huge(0._rk), smzoo_p /= 0._rk))
                          
             do i = 1, self%n_l
                _GET_BOTTOM_(self%id_lgzoo_c(i), lgzoo_c(i)) ! NIAZ: These are in mmol m-3
                _GET_BOTTOM_(self%id_lgzoo_n(i), lgzoo_n(i)) ! NIAZ: These are in mmol m-3
                _GET_BOTTOM_(self%id_lgzoo_p(i), lgzoo_p(i)) ! NIAZ: These are in mmol m-3
             enddo
                          
             ! calculate WW of limiting nutrient for each large zooplankton
             lgzoo_ww_rf = min(merge(lgzoo_c*self%gww_mmolc, huge(0._rk), lgzoo_c /= 0._rk), &
                               merge(lgzoo_n*self%gww_FISHn, huge(0._rk), lgzoo_n /= 0._rk), &
                               merge(lgzoo_p*self%gww_FISHp, huge(0._rk), lgzoo_p /= 0._rk))
             
         else ! if theres only one registered zooplankton prey, get its biomass, and then split evenly between small and large. NIAZ: Ask Ken about this!!!!
             _GET_BOTTOM_(self%id_zoo_c, zoo_c(1)) ! NIAZ: These are in mmol m-3
             _GET_BOTTOM_(self%id_zoo_n, zoo_n(1)) ! NIAZ: These are in mmol m-3
             _GET_BOTTOM_(self%id_zoo_p, zoo_p(1)) ! NIAZ: These are in mmol m-3
             
             !calculate WW of limiting nutrient
             smzoo_ww_rf = min(merge(zoo_c*self%gww_mmolc, huge(0._rk), zoo_c /= 0._rk), &
                               merge(zoo_n*self%gww_FISHn, huge(0._rk), zoo_n /= 0._rk), &
                               merge(zoo_p*self%gww_FISHp, huge(0._rk), zoo_p /= 0._rk)) / 2._rk
             lgzoo_ww_rf = min(merge(zoo_c*self%gww_mmolc, huge(0._rk), zoo_c /= 0._rk), &
                               merge(zoo_n*self%gww_FISHn, huge(0._rk), zoo_n /= 0._rk), &
                               merge(zoo_p*self%gww_FISHp, huge(0._rk), zoo_p /= 0._rk)) / 2._rk
         endif
         smzoo_ww_tot = sum(smzoo_ww_rf)
         lgzoo_ww_tot = sum(lgzoo_ww_rf)

         
         !get benthos (gww m-2)
         if (self%internal_benthos) then
             _GET_BOTTOM_(self%id_benthos,ben_ww_tot)
         else
             do i = 1,self%n_ben
                 _GET_BOTTOM_(self%id_ben_c(i),benthos_c(i))
                 _GET_BOTTOM_(self%id_ben_n(i),benthos_n(i))
                 _GET_BOTTOM_(self%id_ben_p(i),benthos_p(i))
             enddo
             ben_ww_rf = min(merge(benthos_c*self%ben_availability*self%gww_mmolc, huge(0._rk), benthos_c /= 0._rk), &
                             merge(benthos_n*self%ben_availability*self%gww_FISHn, huge(0._rk), benthos_n /= 0._rk), &
                             merge(benthos_p*self%ben_availability*self%gww_FISHp, huge(0._rk), benthos_p /= 0._rk))
             ben_ww_tot = sum(ben_ww_rf)
             !NIAZ: remember to check that we multiply by availability where appropriate
         endif
         
         
         ! Depth-averaged environmental dependencies
         _GET_BOTTOM_(self%id_bottom_depth, bottom_depth) 
         _GET_BOTTOM_(self%id_temp_vertmean_100m, temp_100m)
         _GET_(self%id_temp, temp_bottom)
         _GET_BOTTOM_(self%id_temp_vertmean_500to1500m,temp_500to1500m)
         if(bottom_depth .lt. 500._rk) temp_500to1500m = temp_bottom
         !print*,temp_100m
         _SET_BOTTOM_DIAGNOSTIC_(self%id_temp_vertmean_100m_diag, temp_100m)
         
         ! Adjust photic_depth if bottom is shallower than or equal to photic_depth
         ! Use local variable to avoid modifying shared member variable
         if (bottom_depth .le. self%photic_depth) then
            photic_depth_local = bottom_depth
         else
            photic_depth_local = self%photic_depth
         end if

         !update predator-prey preference matrix theta and vertical distribution
         ! Find the closest pre-computed depth profile
         closest_bottom_depth_idx = minloc(abs(self%depth_array - bottom_depth), dim=1)
         theta= self%theta_matrix(:,:,closest_bottom_depth_idx)
         
         ! From FEISTY
         dvm = photic_depth_local + 500._dp ! 650._dp 
         if (bottom_depth .lt. (photic_depth_local + 500._dp)) dvm = bottom_depth 
         if (bottom_depth .le. self%shelfdepth) dvm = 0._dp

         !update temperature effect
         call updateTempV2(temp_100m, temp_500to1500m, temp_bottom, dvm, bottom_depth, photic_depth_local, ixmedium, ixlarge)
         do iGroup = 1, nGroups
             group(iGroup)%spec%V=group(iGroup)%spec%Vsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%Cmax=group(iGroup)%spec%Cmaxsave*fTempV(ixStart(iGroup):ixEnd(iGroup))
             group(iGroup)%spec%metabolism=group(iGroup)%spec%metabolismsave*fTempmV(ixStart(iGroup):ixEnd(iGroup))
         end do
         call set2vec
         call scale_mortality_rates(self)

         ! FEISTY derivatives
         uin= [smzoo_ww_tot,lgzoo_ww_tot ,ben_ww_tot,0._rk ,fish] ! NIAZ: at this point, zooplankton and benthos biomasses (in WW) must be calculated and aggregated.
         call calcderivatives(uin, dudt)
         !print*,dudt
         !NIAZ: revive the following to inject excess_nutrients? 
         !do i = 1, nGrid
         ! uin(i) = max(0._rk , uin(i))
         !end do
         !! recalculate mortality contribution from each predator
         !mortpred_contri_zoo1= theta(:,1) * Cmax*V/(Enc + Cmax)*uin* uin(1)
         !call checknan(mortpred_contri_zoo1, nGrid)
         !!print*,mortpred_contri_zoo1
         !mortpred_contri_zoo2= theta(:,2) * Cmax*V/(Enc + Cmax)*uin* uin(2)
         !call checknan(mortpred_contri_zoo2, nGrid)
         !!print*,SUM(mortpred_contri_zoo1)/0.01201_rk/9_rk *16._rk/106._rk  

         !do i = 1, nGrid
         uin = max(0._rk , uin)
         !end do
         !! recalculate mortality contribution from each predator, in order to distribute excess nutrient excresion. NIAZ: nGrid-nResources
         pred_morts = Cmax(nResources+1:nGrid) * V(nResources+1:nGrid) / (Enc(nResources+1:nGrid) + Cmax(nResources+1:nGrid)) * uin(nResources+1:nGrid)
         rel_fish_cons = (theta(nResources+1:nGrid,1) + theta(nResources+1:nGrid,2) + theta(nResources+1:nGrid,3)) * pred_morts
         rel_fish_cons = rel_fish_cons/sum(rel_fish_cons)
         
         ! NIAZ: Zooplankton removal should be scaled to the relative WW of the most and least limiting prey nutrients - the difference is then excreted.
         smzoo_mort_ww = dudt(1)/smzoo_ww_tot
         lgzoo_mort_ww = dudt(2)/lgzoo_ww_tot
         
         !calculate total zooplankton mortality as dudt(1)/smzoo_ww_tot, then apply mortality to state variables of prey
         if (self%n_zoo > 1) then ! if multiple zooplankton prey types
             ! update small zooplankton state var ROCs
             do iprey = 1, self%n_s
                do istate = 1, size(self%id_smzoo_int(iprey)%bottom_state)
                   _GET_BOTTOM_(self%id_smzoo_int(iprey)%bottom_state(istate), p1)
                   _ADD_BOTTOM_SOURCE_(self%id_smzoo_int(iprey)%bottom_state(istate), smzoo_mort_ww * p1) ! NIAZ: this assumes dudt(1) is already negative.
                end do
             end do
             ! update large zooplankton state var ROCs
             do iprey = 1, self%n_l
                do istate = 1, size(self%id_lgzoo_int(iprey)%bottom_state)
                   _GET_BOTTOM_(self%id_lgzoo_int(iprey)%bottom_state(istate), p1)
                   _ADD_BOTTOM_SOURCE_(self%id_lgzoo_int(iprey)%bottom_state(istate), lgzoo_mort_ww * p1) ! NIAZ: this assumes dudt(1) is already negative.
                end do
             end do
             
             ! calculate excretion/respiration due to stochiometry
             excess_c_zoo = sum(smzoo_c*self%gww_mmolc * smzoo_mort_ww + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) ! gww biomass * specific loss rate - loss
             excess_n_zoo = sum(smzoo_n*self%gww_FISHn * smzoo_mort_ww + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) ! gww biomass * specific loss rate - loss
             excess_p_zoo = sum(smzoo_p*self%gww_FISHp * smzoo_mort_ww + dudt(1)*smzoo_ww_rf/smzoo_ww_tot) ! gww biomass * specific loss rate - loss
             
             excess_c_zoo = sum(excess_c_zoo + lgzoo_c*self%gww_mmolc * lgzoo_mort_ww + dudt(1)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_mmolc ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
             excess_n_zoo = sum(excess_n_zoo + lgzoo_n*self%gww_FISHn * lgzoo_mort_ww + dudt(1)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_FISHn ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
             excess_p_zoo = sum(excess_p_zoo + lgzoo_p*self%gww_FISHp * lgzoo_mort_ww + dudt(1)*lgzoo_ww_rf/lgzoo_ww_tot) / self%gww_FISHp ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
         else ! if only one zooplankton prey type
             ! Apply mortality to prey state variables
             do istate = 1, size(self%id_zoo_int%bottom_state)
                _GET_BOTTOM_(self%id_zoo_int%bottom_state(istate), p1)
                _ADD_BOTTOM_SOURCE_(self%id_zoo_int%bottom_state(istate), (lgzoo_mort_ww + smzoo_mort_ww)/2._rk * p1) ! NIAZ: this assumes dudt(1) is already negative.
             end do
             
             ! calculate excretion/respiration due to stochiometry
             excess_c_zoo = sum(zoo_c*self%gww_mmolc * (lgzoo_mort_ww + smzoo_mort_ww)/2._rk + (dudt(1)+dudt(2))) / self%gww_mmolc ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
             excess_n_zoo = sum(zoo_n*self%gww_FISHn * (lgzoo_mort_ww + smzoo_mort_ww)/2._rk + (dudt(1)+dudt(2))) / self%gww_FISHn ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
             excess_p_zoo = sum(zoo_p*self%gww_FISHp * (lgzoo_mort_ww + smzoo_mort_ww)/2._rk + (dudt(1)+dudt(2))) / self%gww_FISHp ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
         endif
      
         !calculate total zoobenthos mortality as dudt(1)/smzoo_ww_tot, then apply mortality to state variables of prey
         ben_mort_ww = dudt(3)/ben_ww_tot
         if (self%internal_benthos) then ! if coupling to detritus fluxes, running an internal benthos model.
!NIAZ
             ! benthos dynamics NIAZ: This need to be generalized and rewritten
             det_bot_c = 0._rk
             det_bot_n = 0._rk
             det_bot_p = 0._rk
             if (self%detritus_in_c) _GET_(self%id_det_c, det_bot_c)! mmol C m-3
             if (self%detritus_in_n) _GET_(self%id_det_n, det_bot_n)! mmol N m-3
             if (self%detritus_in_p) _GET_(self%id_det_p, det_bot_p)! mmol P m-3
             ! Use the most limiting element (smallest in gww m-3) for benthos growth calculation
             det_bot = minval([det_bot_c*self%gww_mmolc, det_bot_n*self%gww_FISHn, det_bot_p*self%gww_FISHp], &
                              mask=[self%detritus_in_c .and. det_bot_c*self%gww_mmolc .ne. 0._rk, &
                                    self%detritus_in_n .and. det_bot_n*self%gww_FISHn .ne. 0._rk, &
                                    self%detritus_in_p .and. det_bot_p*self%gww_FISHp .ne. 0._rk])
             if (det_bot .eq. 0._rk) det_bot = minval([det_bot_c*self%gww_mmolc, det_bot_n*self%gww_FISHn, det_bot_p*self%gww_FISHp], &
                                                       mask=[self%detritus_in_c, self%detritus_in_n, self%detritus_in_p])
             _GET_(self%id_depth,depth)
             dflag = 0.5_rk + sign(0.5_rk,det_bot - trcmin*self%gww_FISHn)
             wd = (6.0_rk+6.0e-2_rk*bottom_depth)*d_per_s*dflag
             det_bot_flux_gww=det_bot*wd*dflag !det_bot*w_d   gww m-3 to gww m-3 * m s-1 to gww/m2 s-1
             i = 3!,3!, 4
             benthos_g_gww = 0.1_rk*det_bot_flux_gww*(1-benthos/80)!gww/m2/s
             _ADD_BOTTOM_SOURCE_(self%id_benthos, benthos_g_gww + dudt(i) )!gww/m2  dRdt(3) = rr(3)*(1-R(3)/K(3)) - mortRes(3)*R(3) K is 80
             if (self%detritus_in_c) _ADD_BOTTOM_FLUX_(self%id_det_c, -det_bot_flux_gww/self%gww_mmolc) ! gww/m2 s-1 to mmol C/m2 s-1
             if (self%detritus_in_n) _ADD_BOTTOM_FLUX_(self%id_det_n, -det_bot_flux_gww/self%gww_FISHn) ! gww/m2 s-1 to mmol N/m2 s-1
             if (self%detritus_in_p) _ADD_BOTTOM_FLUX_(self%id_det_p, -det_bot_flux_gww/self%gww_FISHp) ! gww/m2 s-1 to mmol P/m2 s-1
             
             benthos_loss_gww = 0.1_rk*det_bot_flux_gww*(1-(1-benthos/80)) !gww, convert below
             if (self%nutrient_in_c) _ADD_BOTTOM_FLUX_(self%id_nut_c, (0.9_rk *det_bot_flux_gww + benthos_loss_gww) / self%gww_mmolc) ! gww/m2 s-1 to mmol C/m2 s-1
             if (self%nutrient_in_n) _ADD_BOTTOM_FLUX_(self%id_nut_n, (0.9_rk *det_bot_flux_gww + benthos_loss_gww) / self%gww_FISHn) ! gww/m2 s-1 to mmol N/m2 s-1
             if (self%nutrient_in_p) _ADD_BOTTOM_FLUX_(self%id_nut_p, (0.9_rk *det_bot_flux_gww + benthos_loss_gww) / self%gww_FISHp) ! gww/m2 s-1 to mmol P/m2 s-1
!NIAZ
         else ! if using external benthos model(s)
             do iprey = 1, self%n_ben
                 do istate = 1, size(self%id_benprey(iprey)%bottom_state)
                     _GET_BOTTOM_(self%id_benprey(iprey)%bottom_state(istate), p1)
                     _ADD_BOTTOM_SOURCE_(self%id_benprey(iprey)%bottom_state(istate), ben_mort_ww * p1) ! NIAZ: this assumes dudt(3) is already negative.
                 end do
                 
                 ! calculate excretion/respiration due to stochiometry
                 excess_c_ben = sum(benthos_c*self%ben_availability*self%gww_mmolc * ben_mort_ww + dudt(3)) / self%gww_mmolc ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
                 excess_n_ben = sum(benthos_n*self%ben_availability*self%gww_FISHn * ben_mort_ww + dudt(3)) / self%gww_FISHn ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
                 excess_p_ben = sum(benthos_p*self%ben_availability*self%gww_FISHp * ben_mort_ww + dudt(3)) / self%gww_FISHp ! gww biomass * specific loss rate - loss. NIAZ: in mmol m-3
             end do
         endif
         excess_c = excess_c_ben + excess_c_zoo
         excess_n = excess_n_ben + excess_n_zoo
         excess_p = excess_p_ben + excess_p_zoo
         !NIAZ: remember to add excess_ nutrients to their target pools:
         !  - calculate relative excess nutrient contributions from each fish, and add them to respiration fluxes.
         !      excess_x*rel_fish_cons
         
         do i = 1, nGrid-nResources
             ! fish dynamics
             _ADD_BOTTOM_SOURCE_(self%id_fish(i), dudt(i+nResources) )!gww/m2
             
             ! excretion/respiration
             if (index(self%respiration_elements, 'c') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_c(i), respiration(i+nResources) / self%gww_mmolc + excess_c*rel_fish_cons(i))!* gwwC /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to target C/m2
             endif
             if (index(self%respiration_elements, 'n') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_n(i), respiration(i+nResources) / self%gww_FISHn + excess_n*rel_fish_cons(i))
             endif
             if (index(self%respiration_elements, 'p') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_respiration_fish_p(i), respiration(i+nResources) / self%gww_FISHp + excess_p*rel_fish_cons(i))
             endif
             
             ! feces
             if (index(self%feces_elements, 'c') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_feces_fish_c(i), feces(i+nResources) / self%gww_mmolc)!gww/m2 to target C/m2
             endif
             if (index(self%feces_elements, 'n') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_feces_fish_n(i), feces(i+nResources) / self%gww_FISHn)!* self%qnc /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
             endif
             if (index(self%feces_elements, 'p') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_feces_fish_p(i), feces(i+nResources) / self%gww_FISHp)!gww/m2 to mmol P/m2
             endif
             
             ! carcasses
             if (index(self%carcass_elements, 'c') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_c(i), carcasses(i+nResources) / self%gww_mmolc)!gww/m2 to target C/m2
             endif
             if (index(self%carcass_elements, 'n') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_n(i), carcasses(i+nResources) / self%gww_FISHn)!* self%qnc /0.01201_rk/9_rk *16._rk/106._rk)!gww/m2 to mmol N/m2
             endif
             if (index(self%carcass_elements, 'p') /= 0) then
                 _ADD_BOTTOM_SOURCE_(self%id_carcasses_fish_p(i), carcasses(i+nResources) / self%gww_FISHp)!gww/m2 to mmol P/m2
             endif
         end do
         
         do i=1,nGroups            
            _SET_BOTTOM_DIAGNOSTIC_(self%id_fish_total_biomass(i), totBiomass(i))
         end do
      _BOTTOM_LOOP_END_
   end subroutine   

   subroutine initialize_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)
      class (type_feisty_fabm), intent(inout) :: self
      real(rk), intent(in) :: szprod, lzprod, bprodin, dfbot, dfpho
      integer, intent(in)  :: nStages
      real(rk), intent(in) :: Tp, Tm, Tb, etamature, visual, Fmax, etaF
      logical :: file_exists
      integer :: dims(3)

      call allocate_vertical_arrays(self) 
      self%theta_matrix = 0._rk
      self%vertical_distribution_matrix_centers = -999._rk

      inquire(file='matrix.dat', exist=file_exists)
      if (file_exists) then
         print*, "Reading matrices from existing file..."
         open(unit=10, file='matrix.dat', form='unformatted', access='stream', status='old')
         read(10) dims
         if (dims(1)==size(self%theta_matrix,1) .and. dims(2)==size(self%theta_matrix,2) .and. dims(3)==size(self%theta_matrix,3)) then
            print *, "Dimensions match.", dims, ". Reading matrices..."
            read(10) self%theta_matrix
            read(10) self%vertical_distribution_matrix_centers
            close(10)
            print *, "Matrices read successfully."
         else
            print *, "Dimension mismatch! Expected:", size(self%theta_matrix,1),size(self%theta_matrix,2),size(self%theta_matrix,3), " Got:", dims, "."
            close(10)
            file_exists = .false.
         end if
      end if

      if (.not. file_exists) then

         call build_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)
      end if

      print *, "Matrices (theta, vertical_distribution) ready for computation."
   end subroutine initialize_theta_matrices

   ! Allocate depth arrays for pre-computed matrices
   subroutine allocate_vertical_arrays(self)
      ! Creates a logarithmically-spaced depth grid for efficient lookup of depth-dependent fish behavior.
      class (type_feisty_fabm), intent(inout) :: self
      logical :: changed
      integer :: i_depth, j_depth
      logical :: is_unique

      ! Allocate depth arrays for pre-computed matrices
      allocate (self%depth_array(self%n_depth))  ! Depth sampling points (e.g., 50 depths from 1m to 5000m)
      allocate (self%theta_matrix(nGrid,nGrid,self%n_depth))  ! Predation preference matrix (predator×prey×depth)
      allocate (self%vertical_distribution_matrix_centers(int(self%max_depth),nFGrid,self%n_depth))  ! Vertical distributions (depth×fish_classes×depth_profile)

      ! Create logarithmically-spaced depth array: more points near surface, fewer in deep ocean
      self%depth_array = [(10._rk**(real(i_depth-1,rk)*(log10(self%max_depth))/real(self%n_depth-1,rk)), i_depth=1, self%n_depth)]
      self%depth_array = floor(self%depth_array)  ! Round down to integer depths
      self%depth_array(size(self%depth_array)) = self%max_depth  ! Ensure last value is exactly max_depth

      ! Remove duplicate depths that may arise from flooring
      do
         changed = .false.
         do i_depth = 1, self%n_depth - 1
            if (self%depth_array(i_depth) >= self%depth_array(i_depth+1)) then
               self%depth_array(i_depth+1) = min(self%depth_array(i_depth) + 1._rk, self%max_depth)  ! Increment by 1m
               changed = .true.
            end if
         end do
         if (.not. changed) exit  ! Stop when no more duplicates found
      end do

      ! Verify all depths are unique (quality check)
      is_unique = .true.
      do i_depth = 1, self%n_depth - 1
         do j_depth = i_depth + 1, self%n_depth
            if (self%depth_array(i_depth) == self%depth_array(j_depth)) then
               print*, "WARNING: depth_array still contains duplicates after removal loop!"
               is_unique = .false.
               exit
            end if
         end do
         if (.not. is_unique) exit
      end do
   end subroutine allocate_vertical_arrays

   subroutine build_theta_matrices(self, szprod, lzprod, bprodin, dfbot, dfpho, nStages, Tp, Tm, Tb, etamature, visual, Fmax, etaF)
      class (type_feisty_fabm), intent(inout) :: self
      real(rk), intent(in) :: szprod, lzprod, bprodin, dfbot, dfpho
      integer, intent(in)  :: nStages
      real(rk), intent(in) :: Tp, Tm, Tb, etamature, visual, Fmax, etaF
      real(rk) :: depth_local, photic_depth_local
      integer  :: i_depth

      print*, "Preparing the theta matrix data file..."
      do i_depth = 1, size(self%theta_matrix,3)
         depth_local = self%depth_array(i_depth)
         photic_depth_local = min(depth_local, self%photic_depth)
         if (depth_local .le. self%photic_depth) then
            call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,nStages,Tp,Tm,Tb,depth_local,depth_local,etamature,self%shelfdepth,visual,Fmax,etaF)
         else
            call setupVertical2(szprod,lzprod,bprodin,dfbot,dfpho,nStages,Tp,Tm,Tb,depth_local,photic_depth_local,etamature,self%shelfdepth,visual,Fmax,etaF)
         end if

         self%theta_matrix(:,:,i_depth)=theta
         call update_vertical_distribution_centers(self, i_depth)

         if(depth_local.le.self%shelfdepth) self%theta_matrix(ixStart(2):ixEnd(2),:,i_depth)=0._rk
         if(depth_local.le.self%shelfdepth) self%theta_matrix(ixStart(4):ixEnd(4),:,i_depth)=0._rk
      end do

      call save_theta_matrices(self)
      print *, "Matrices saved successfully."
   end subroutine build_theta_matrices

   subroutine update_vertical_distribution_centers(self, depth_index)
      class (type_feisty_fabm), intent(inout) :: self
      integer, intent(in) :: depth_index
      real(rk), allocatable :: vertical_dist_at_boundary(:,:)
      real(rk) :: column_sum(nFGrid)
      integer :: j, n_depth_bounds

      n_depth_bounds = size(depthDay,1)
      if (n_depth_bounds <= 1) return

      allocate(vertical_dist_at_boundary(n_depth_bounds, nFGrid))
      vertical_dist_at_boundary = (depthDay(:,idxF:nGrid) + depthNight(:,idxF:nGrid)) / 2.0_rk
      self%vertical_distribution_matrix_centers(1:n_depth_bounds-1,:,depth_index) = &
         (vertical_dist_at_boundary(1:n_depth_bounds-1,:) + vertical_dist_at_boundary(2:n_depth_bounds,:)) / 2.0_rk

      do j = 1, nFGrid
         column_sum(j) = sum(self%vertical_distribution_matrix_centers(1:n_depth_bounds-1,j,depth_index))
         if (column_sum(j) > 0.0_rk) then
            self%vertical_distribution_matrix_centers(1:n_depth_bounds-1,j,depth_index) = &
               self%vertical_distribution_matrix_centers(1:n_depth_bounds-1,j,depth_index) / column_sum(j)
         end if
      end do
      deallocate(vertical_dist_at_boundary)
   end subroutine update_vertical_distribution_centers

   subroutine save_theta_matrices(self)
      class (type_feisty_fabm), intent(in) :: self
      integer :: dims(3)

      print *, "Saving matrices to file..."
      dims = [size(self%theta_matrix,1),size(self%theta_matrix,2),size(self%theta_matrix,3)]
      open(unit=10, file='matrix.dat', form='unformatted', access='stream', status='replace')
      write(10) dims
      write(10) self%theta_matrix
      write(10) self%vertical_distribution_matrix_centers
      close(10)
   end subroutine save_theta_matrices

   subroutine scale_mortality_rates(self)
      class (type_feisty_fabm), intent(in) :: self
      if (self%bgmort .ne. mort0(nGrid)) then
         mort0(idxF:nGrid) = self%bgmort/seconds_per_year !NIAZ change to scale factor at parameter registration
      else
         mort0(idxF:nGrid) = mort0(idxF:nGrid)/seconds_per_year !NIAZ change to scale factor at parameter registration
      end if
      mortF(idxF:nGrid) = mortF(idxF:nGrid)/seconds_per_year !NIAZ change to scale factor at parameter registration
   end subroutine scale_mortality_rates

   ! Add model subroutines here.

end module
